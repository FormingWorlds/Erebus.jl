using Test
using Erebus
using Erebus.Config
using Erebus.Geometry
using Erebus.Physics
using Erebus.Particles
using Erebus.Numerics
using ExtendableSparse
using SparseArrays
using StaticArrays

@testset "Cold Lid Hydrofracture Venting & Ice Seals" begin
    @testset "compute_ice_sealed_permeability Invariants & Asymptotics" begin
        k0 = 1.0e-11
        T_freeze = 273.15
        dT_seal = 10.0
        k_min_ratio = 1.0e-6

        # Temperature at or above freezing: unsealed reference permeability
        @test isapprox(
            compute_ice_sealed_permeability(
                k0,
                T_freeze;
                T_freeze=T_freeze,
                delta_T_seal=dT_seal,
                k_min_ratio=k_min_ratio,
            ),
            k0;
            rtol=1e-12,
        )
        @test isapprox(
            compute_ice_sealed_permeability(
                k0, 300.0; T_freeze=T_freeze, delta_T_seal=dT_seal, k_min_ratio=k_min_ratio
            ),
            k0;
            rtol=1e-12,
        )

        # Temperature below freezing: exponential sealing
        T_cold1 = T_freeze - dT_seal # 1 e-folding drop
        expected_ratio1 = (1.0 - k_min_ratio) * exp(-1.0) + k_min_ratio
        k_sealed1 = compute_ice_sealed_permeability(
            k0, T_cold1; T_freeze=T_freeze, delta_T_seal=dT_seal, k_min_ratio=k_min_ratio
        )
        @test isapprox(k_sealed1, k0 * expected_ratio1; rtol=1e-10)
        @test k_sealed1 < k0

        # Monotonicity test: lower temperatures yield lower permeability
        T_cold2 = T_freeze - 2.0 * dT_seal
        k_sealed2 = compute_ice_sealed_permeability(
            k0, T_cold2; T_freeze=T_freeze, delta_T_seal=dT_seal, k_min_ratio=k_min_ratio
        )
        @test k_sealed2 < k_sealed1

        # Asymptotic cryogenic limit: T << T_freeze approaches k0 * k_min_ratio
        T_deep_freeze = 50.0 # 50 K surface
        k_sealed_deep = compute_ice_sealed_permeability(
            k0,
            T_deep_freeze;
            T_freeze=T_freeze,
            delta_T_seal=dT_seal,
            k_min_ratio=k_min_ratio,
        )
        expected_deep =
            k0 *
            ((1.0 - k_min_ratio) * exp(-(T_freeze - T_deep_freeze) / dT_seal) + k_min_ratio)
        @test isapprox(k_sealed_deep, expected_deep; rtol=1e-12)
        @test isapprox(k_sealed_deep, k0 * k_min_ratio; rtol=1e-3)
        @test k_sealed_deep >= k0 * k_min_ratio

        # DomainError guards on non-physical inputs
        @test_throws DomainError compute_ice_sealed_permeability(0.0, 200.0)
        @test_throws DomainError compute_ice_sealed_permeability(-1.0e-11, 200.0)
        @test_throws DomainError compute_ice_sealed_permeability(k0, 0.0)
        @test_throws DomainError compute_ice_sealed_permeability(k0, -50.0)
        @test_throws DomainError compute_ice_sealed_permeability(k0, 200.0; T_freeze=0.0)
        @test_throws DomainError compute_ice_sealed_permeability(
            k0, 200.0; delta_T_seal=0.0
        )
        @test_throws DomainError compute_ice_sealed_permeability(k0, 200.0; k_min_ratio=0.0)
        @test_throws DomainError compute_ice_sealed_permeability(k0, 200.0; k_min_ratio=1.5)
        @test_throws DomainError compute_ice_sealed_permeability(NaN, 200.0)
        @test_throws DomainError compute_ice_sealed_permeability(k0, NaN)
    end

    @testset "is_hydrofracture_breached Invariants" begin
        sigma_t = 1.0e7 # 10 MPa tensile strength

        # Compressive regime (Pt > Pf, Peff > 0): intact lid
        Peff_comp = 5.0e6
        @test !is_hydrofracture_breached(Peff_comp, sigma_t)
        @test !is_hydrofracture_breached(10.0e6, 5.0e6, sigma_t)

        # Mild tensile regime below rupture (-sigma_t < Peff <= 0): intact lid
        Peff_mild = -5.0e6
        @test !is_hydrofracture_breached(Peff_mild, sigma_t)
        @test !is_hydrofracture_breached(5.0e6, 10.0e6, sigma_t)

        # Rupture threshold (Peff == -sigma_t): breached
        Peff_thresh = -1.0e7
        @test is_hydrofracture_breached(Peff_thresh, sigma_t)
        @test is_hydrofracture_breached(5.0e6, 15.0e6, sigma_t)

        # Severe overpressure (Peff < -sigma_t): breached
        Peff_severe = -2.5e7
        @test is_hydrofracture_breached(Peff_severe, sigma_t)
        @test is_hydrofracture_breached(5.0e6, 30.0e6, sigma_t)

        # Non-finite and invalid inputs: safely return false
        @test !is_hydrofracture_breached(NaN, sigma_t)
        @test !is_hydrofracture_breached(Peff_severe, NaN)
        @test !is_hydrofracture_breached(Peff_severe, 0.0)
        @test !is_hydrofracture_breached(Peff_severe, -1.0e7)
    end

    @testset "Surface Boundary Hydrofracture Gating & Cryogenic Sealing Assembly" begin
        Nx, Ny = 9, 9
        Nx1, Ny1 = Nx + 1, Ny + 1
        coords = GridCoordinates(Nx, Ny; xsize=140_000.0, ysize=140_000.0)
        rplanet = 50_000.0
        xcenter = 70_000.0
        ycenter = 70_000.0
        P_amb = 10.0
        k_vent = 1.0e-11
        conductance_factor = 1.0
        Kcont = 1.0e20
        n_dof = 6 * Ny1 * Nx1

        tk_cold = fill(150.0, Ny1, Nx1) # Sub-freezing surface
        sigma_t_val = 1.0e7
        TEN = fill(sigma_t_val, Ny1, Nx1)
        pr_litho = fill(2.0e6, Ny1, Nx1)

        # Case 1: Hydrofracture-gated mode with sub-tensile pore pressure (Peff > -sigma_t)
        # Lid is intact: Robin conductance must NOT be added to L or R, S_vent must be zero
        pf_subtensile = fill(5.0e6, Ny1, Nx1) # Peff = 2 MPa - 5 MPa = -3 MPa > -10 MPa
        L_gated_closed = ExtendableSparseMatrix{Float64,Int64}(n_dof, n_dof)
        R_gated_closed = zeros(Float64, n_dof)
        S_vent_gated_closed = zeros(Float64, Ny1, Nx1)

        apply_venting_surface_boundary!(
            L_gated_closed,
            R_gated_closed,
            tk_cold,
            coords,
            rplanet,
            xcenter,
            ycenter,
            P_amb;
            k_vent=k_vent,
            conductance_factor=conductance_factor,
            mode=:hydrofracture_gated,
            pr=pr_litho,
            pf=pf_subtensile,
            TEN=TEN,
            Kcont=Kcont,
            S_vent_out=S_vent_gated_closed,
        )
        flush!(L_gated_closed)
        @test iszero(SparseArrays.nnz(L_gated_closed.cscmatrix))
        @test all(iszero, R_gated_closed)
        @test all(iszero, S_vent_gated_closed)

        # Case 2: Hydrofracture-gated mode with tensile failure (Peff <= -sigma_t)
        # Lid is breached: Robin conductance must be assembled, S_vent must be positive
        pf_breached = fill(20.0e6, Ny1, Nx1) # Peff = 2 MPa - 20 MPa = -18 MPa <= -10 MPa
        L_gated_open = ExtendableSparseMatrix{Float64,Int64}(n_dof, n_dof)
        R_gated_open = zeros(Float64, n_dof)
        S_vent_gated_open = zeros(Float64, Ny1, Nx1)

        apply_venting_surface_boundary!(
            L_gated_open,
            R_gated_open,
            tk_cold,
            coords,
            rplanet,
            xcenter,
            ycenter,
            P_amb;
            k_vent=k_vent,
            conductance_factor=conductance_factor,
            mode=:hydrofracture_gated,
            pr=pr_litho,
            pf=pf_breached,
            TEN=TEN,
            Kcont=Kcont,
            S_vent_out=S_vent_gated_open,
        )
        flush!(L_gated_open)
        @test SparseArrays.nnz(L_gated_open.cscmatrix) > 0
        @test any(R_gated_open .> 0.0)
        @test any(S_vent_gated_open .> 0.0)

        # Case 3: Darcy sink mode with cryogenic ice sealing enabled vs disabled
        # When unbreached, ice sealing must suppress surface venting rate by ~1e-6
        S_vent_unsealed = zeros(Float64, Ny1, Nx1)
        apply_venting_surface_boundary!(
            nothing,
            nothing,
            tk_cold,
            coords,
            rplanet,
            xcenter,
            ycenter,
            P_amb;
            k_vent=k_vent,
            conductance_factor=conductance_factor,
            mode=:darcy_sink,
            ice_sealing=false,
            pr=pr_litho,
            pf=pf_subtensile,
            TEN=TEN,
            S_vent_out=S_vent_unsealed,
        )

        S_vent_sealed = zeros(Float64, Ny1, Nx1)
        apply_venting_surface_boundary!(
            nothing,
            nothing,
            tk_cold,
            coords,
            rplanet,
            xcenter,
            ycenter,
            P_amb;
            k_vent=k_vent,
            conductance_factor=conductance_factor,
            mode=:darcy_sink,
            ice_sealing=true,
            t_freeze=273.15,
            dt_seal=10.0,
            k_seal_min_ratio=1.0e-6,
            pr=pr_litho,
            pf=pf_subtensile,
            TEN=TEN,
            S_vent_out=S_vent_sealed,
        )

        # Sealed venting must be strictly smaller and scaled by ice sealing permeability ratio
        max_unsealed = maximum(S_vent_unsealed)
        max_sealed = maximum(S_vent_sealed)
        expected_seal_ratio = compute_ice_sealed_permeability(
            1.0, 150.0; T_freeze=273.15, delta_T_seal=10.0, k_min_ratio=1.0e-6
        )
        @test max_unsealed > 0.0
        @test max_sealed > 0.0
        @test max_sealed < max_unsealed
        @test isapprox(max_sealed / max_unsealed, expected_seal_ratio; rtol=1e-5)

        # Case 4: Darcy sink mode with breached pressure - hydrofracture flag gating
        # When hydrofracture=false (default), Darcy sink retains reference k_vent even if Peff <= -sigma_t
        # When hydrofracture=true, permeability enhances to hydrofracture value
        S_vent_darcy_plain = zeros(Float64, Ny1, Nx1)
        apply_venting_surface_boundary!(
            nothing,
            nothing,
            tk_cold,
            coords,
            rplanet,
            xcenter,
            ycenter,
            P_amb;
            k_vent=k_vent,
            conductance_factor=conductance_factor,
            mode=:darcy_sink,
            hydrofracture=false,
            ice_sealing=false,
            pr=pr_litho,
            pf=pf_breached,
            TEN=TEN,
            S_vent_out=S_vent_darcy_plain,
        )

        S_vent_darcy_enhanced = zeros(Float64, Ny1, Nx1)
        apply_venting_surface_boundary!(
            nothing,
            nothing,
            tk_cold,
            coords,
            rplanet,
            xcenter,
            ycenter,
            P_amb;
            k_vent=k_vent,
            conductance_factor=conductance_factor,
            mode=:darcy_sink,
            hydrofracture=true,
            ice_sealing=false,
            pr=pr_litho,
            pf=pf_breached,
            TEN=TEN,
            S_vent_out=S_vent_darcy_enhanced,
        )

        max_plain = maximum(S_vent_darcy_plain)
        max_enhanced = maximum(S_vent_darcy_enhanced)
        @test max_plain > 0.0
        @test max_enhanced > max_plain
        @test isapprox(max_enhanced / max_plain, 100.0; rtol=1e-4)
    end

    @testset "compute_face_venting_permeability Invariants" begin
        k0 = 1.0e-11
        # 1. Unbreached, unsealed -> k0
        @test isapprox(
            compute_face_venting_permeability(k0, false, false, 200.0, 0.0, 1.0e7),
            k0;
            rtol=1e-12,
        )
        # 2. Unbreached, ice-sealed -> sealed k
        k_sealed = compute_face_venting_permeability(
            k0,
            false,
            true,
            200.0,
            0.0,
            1.0e7;
            t_freeze=273.15,
            dt_seal=10.0,
            k_seal_min_ratio=1.0e-6,
        )
        @test k_sealed < k0
        @test isapprox(
            k_sealed,
            compute_ice_sealed_permeability(
                k0, 200.0; T_freeze=273.15, delta_T_seal=10.0, k_min_ratio=1.0e-6
            );
            rtol=1e-12,
        )
        # 3. Breached -> enhanced hydrofracture permeability regardless of ice_sealing
        k_breached = compute_face_venting_permeability(
            k0,
            true,
            true,
            150.0,
            -18.0e6,
            10.0e6;
            kappa_frac=1.0e3,
            gamma_frac=1.0,
            k_frac_max=1.0e-9,
        )
        @test isapprox(k_breached, 1.0e-9; rtol=1e-12)
    end

    @testset "VentingConfig Validation on Ice Sealing Parameters" begin
        quick_toml = joinpath(@__DIR__, "..", "configs", "test_quick.toml")
        cfg_base = load_config(quick_toml)

        # Invalid t_freeze <= 0 or non-finite
        cfg_bad_tf = SimulationConfig(
            grid=cfg_base.grid,
            geometry=cfg_base.geometry,
            time=cfg_base.time,
            solver=cfg_base.solver,
            poroelasticity=cfg_base.poroelasticity,
            thermodynamics=cfg_base.thermodynamics,
            reaction=cfg_base.reaction,
            materials=cfg_base.materials,
            output=cfg_base.output,
            disk=cfg_base.disk,
            melting=cfg_base.melting,
            venting=VentingConfig(t_freeze=0.0),
        )
        @test_throws ArgumentError validate_config(cfg_bad_tf)

        # Invalid dt_seal <= 0
        cfg_bad_dt = SimulationConfig(
            grid=cfg_base.grid,
            geometry=cfg_base.geometry,
            time=cfg_base.time,
            solver=cfg_base.solver,
            poroelasticity=cfg_base.poroelasticity,
            thermodynamics=cfg_base.thermodynamics,
            reaction=cfg_base.reaction,
            materials=cfg_base.materials,
            output=cfg_base.output,
            disk=cfg_base.disk,
            melting=cfg_base.melting,
            venting=VentingConfig(dt_seal=-5.0),
        )
        @test_throws ArgumentError validate_config(cfg_bad_dt)

        # Invalid k_seal_min_ratio <= 0 or > 1.0
        cfg_bad_ratio0 = SimulationConfig(
            grid=cfg_base.grid,
            geometry=cfg_base.geometry,
            time=cfg_base.time,
            solver=cfg_base.solver,
            poroelasticity=cfg_base.poroelasticity,
            thermodynamics=cfg_base.thermodynamics,
            reaction=cfg_base.reaction,
            materials=cfg_base.materials,
            output=cfg_base.output,
            disk=cfg_base.disk,
            melting=cfg_base.melting,
            venting=VentingConfig(k_seal_min_ratio=0.0),
        )
        @test_throws ArgumentError validate_config(cfg_bad_ratio0)

        cfg_bad_ratio2 = SimulationConfig(
            grid=cfg_base.grid,
            geometry=cfg_base.geometry,
            time=cfg_base.time,
            solver=cfg_base.solver,
            poroelasticity=cfg_base.poroelasticity,
            thermodynamics=cfg_base.thermodynamics,
            reaction=cfg_base.reaction,
            materials=cfg_base.materials,
            output=cfg_base.output,
            disk=cfg_base.disk,
            melting=cfg_base.melting,
            venting=VentingConfig(k_seal_min_ratio=1.5),
        )
        @test_throws ArgumentError validate_config(cfg_bad_ratio2)
    end

    @testset "Episodic Hydrofracture Breaching & Lid Resealing Dynamics" begin
        # Simulate an episodic venting cycle at a surface cell:
        # Step 1: Pressure buildup triggers hydrofracturing (Peff <= -sigma_t) -> lid opens
        # Step 2: Venting releases pore fluid, lowering Pf (Peff > -sigma_t) -> lid reseals
        Nx, Ny = 7, 7
        coords = GridCoordinates(Nx, Ny; xsize=100_000.0, ysize=100_000.0)
        rplanet = 35_000.0
        xcenter = 50_000.0
        ycenter = 50_000.0
        P_amb = 1.0e-4 # space vacuum
        k_vent = 1.0e-11
        tk = fill(180.0, Ny + 1, Nx + 1)
        TEN = fill(1.0e7, Ny + 1, Nx + 1)
        pr = fill(2.0e6, Ny + 1, Nx + 1)

        # Pulse 1: Overpressure breach
        pf_pulse1 = fill(15.0e6, Ny + 1, Nx + 1) # Peff = 2 - 15 = -13 MPa <= -10 MPa
        S_vent_pulse1 = zeros(Float64, Ny + 1, Nx + 1)
        apply_venting_surface_boundary!(
            nothing,
            nothing,
            tk,
            coords,
            rplanet,
            xcenter,
            ycenter,
            P_amb;
            k_vent=k_vent,
            mode=:hydrofracture_gated,
            ice_sealing=true,
            pr=pr,
            pf=pf_pulse1,
            TEN=TEN,
            S_vent_out=S_vent_pulse1,
        )
        @test maximum(S_vent_pulse1) > 0.0

        # Pulse 2: Resealing after pressure dissipation
        pf_pulse2 = fill(8.0e6, Ny + 1, Nx + 1) # Peff = 2 - 8 = -6 MPa > -10 MPa
        S_vent_pulse2 = zeros(Float64, Ny + 1, Nx + 1)
        apply_venting_surface_boundary!(
            nothing,
            nothing,
            tk,
            coords,
            rplanet,
            xcenter,
            ycenter,
            P_amb;
            k_vent=k_vent,
            mode=:hydrofracture_gated,
            ice_sealing=true,
            pr=pr,
            pf=pf_pulse2,
            TEN=TEN,
            S_vent_out=S_vent_pulse2,
        )
        @test all(iszero, S_vent_pulse2)
    end

    @testset "Simulation Loop with Hydrofracture-Gated & Ice Sealed Venting" begin
        output_dir = mktempdir()
        try
            quick_toml = joinpath(@__DIR__, "..", "configs", "test_quick.toml")
            cfg = load_config(quick_toml)
            cfg_hf = SimulationConfig(
                grid=cfg.grid,
                geometry=cfg.geometry,
                time=TimeConfig(
                    dt_initial=cfg.time.dt_initial,
                    dt_longest=cfg.time.dt_longest,
                    dtcoefdn=cfg.time.dtcoefdn,
                    dtcoefup=cfg.time.dtcoefup,
                    dtstep=cfg.time.dtstep,
                    dxymax=cfg.time.dxymax,
                    vpratio=cfg.time.vpratio,
                    DTmax=cfg.time.DTmax,
                    start_time=cfg.time.start_time,
                    endtime=cfg.time.endtime,
                    start_step=1,
                    n_steps=2,
                ),
                solver=cfg.solver,
                poroelasticity=cfg.poroelasticity,
                thermodynamics=cfg.thermodynamics,
                reaction=cfg.reaction,
                materials=cfg.materials,
                output=OutputConfig(output_dir=output_dir, savematstep=2),
                disk=DiskConfig(
                    t_dispersal_myr=0.01,
                    dt_dispersal_myr=0.002,
                    p_amb_disk=10.0,
                    p_amb_space=1.0e-4,
                    dispersal_active=true,
                ),
                melting=cfg.melting,
                venting=VentingConfig(
                    active=true,
                    mode=:hydrofracture_gated,
                    k_vent=1.0e-11,
                    conductance_factor=1.0,
                    L_sublimation=2.83e6,
                    latent_cooling=true,
                    ice_sealing=true,
                    t_freeze=273.15,
                    dt_seal=10.0,
                    k_seal_min_ratio=1.0e-6,
                ),
            )
            Erebus.simulation_loop(cfg_hf; output_path=output_dir)

            files = readdir(output_dir)
            @test "output_00000.jld2" in files
            @test "output_00002.jld2" in files

            data2 = load_state(joinpath(output_dir, "output_00002.jld2"))
            @test data2["timestep"] == 2
            @test haskey(data2, "M_vent_total")
            @test isfinite(data2["M_vent_total"])
            @test data2["M_vent_total"] >= 0.0
            @test haskey(data2, "S_vent")
            @test all(isfinite, data2["S_vent"])
            @test all(data2["S_vent"] .>= 0.0)
            @test !any(isnan, data2["tk2"])
            @test !any(isinf, data2["tk2"])
            @test all(data2["tk2"] .> 0.0)
            @test !any(isnan, data2["pf"])
            @test !any(isinf, data2["pf"])
            @test !any(isnan, data2["pr"])
            @test !any(isinf, data2["pr"])
        finally
            rm(output_dir; recursive=true, force=true)
        end
    end

    @testset "NaN and non-positive temperature handling in boundary assembly" begin
        Nx, Ny = 5, 5
        coords = GridCoordinates(Nx, Ny; xsize=100_000.0, ysize=100_000.0)
        tk_nan = fill(150.0, Ny + 1, Nx + 1)
        tk_nan[3, 3] = NaN
        tk_nan[2, 2] = -50.0
        S_out = zeros(Float64, Ny + 1, Nx + 1)
        # Must execute without throwing DomainError
        apply_venting_surface_boundary!(
            nothing,
            nothing,
            tk_nan,
            coords,
            35_000.0,
            50_000.0,
            50_000.0,
            10.0;
            ice_sealing=true,
            S_vent_out=S_out,
        )
        @test all(isfinite, S_out)
        @test all(S_out .>= 0.0)
    end
end

using Test
using Erebus
using StaticArrays
using LinearAlgebra

include("test_helpers.jl")

@testset "Magma Ocean Multi-Component Degassing Physics" begin
    # Reference planet parameters (lunar-sized / Vesta-sized planetesimal)
    R_p = 50_000.0 # 50 km radius
    M_p = (4.0 / 3.0) * π * (R_p^3) * 3300.0 # ~1.73e18 kg
    g_surf = Erebus.G_GRAV * M_p / (R_p^2)
    T_mo = 1800.0 # Liquidus / super-liquidus temperature [K]
    M_melt = 0.8 * M_p # 80% magma ocean by mass

    # Test volatile inventories (representative carbonaceous chondrite volatile content)
    # H: 5000 ppm (as H2O equiv ~4.5 wt%)
    # C: 2000 ppm
    # N: 200 ppm
    # S: 10000 ppm (1 wt%)
    M_tot_H = M_p * 5.0e-3
    M_tot_C = M_p * 2.0e-3
    M_tot_N = M_p * 2.0e-4
    M_tot_S = M_p * 1.0e-2

    # -------------------------------------------------------------------------
    # 1. Equilibrium Magma Ocean Partitioning & Mass Conservation
    # -------------------------------------------------------------------------
    @testset "Equilibrium Volatile Partitioning & Mass Conservation" begin
        # Solve partitioning at neutral IW (delta_IW = 0.0)
        sol = solve_magma_ocean_volatile_partitioning(
            M_melt,
            M_tot_H,
            M_tot_C,
            M_tot_N,
            M_tot_S,
            R_p,
            g_surf,
            T_mo,
            0.0;
            carbon_active=true,
            sulfur_active=true,
        )

        @test sol isa NamedTuple
        @test sol.P_surf > 0.0
        @test isfinite(sol.P_surf)

        # Dalton's law: sum of partial pressures equals total surface pressure
        P_sum = sum(values(sol.p_i))
        @test isapprox(P_sum, sol.P_surf; rtol=1e-10)

        # Strict elemental mass conservation across melt and atmosphere:
        # 1. Hydrogen conservation (elemental H mass in melt + atmosphere)
        # H2O is 2.01588 / 18.01528 H by mass; H2 is 1.0 H by mass; CH4 is 4.032 / 16.042 H; NH3 is 3.024 / 17.031 H; H2S is 2.016 / 34.08 H
        M_atm_H = (
            get(sol.M_atm_i, :H2, 0.0) +
            get(sol.M_atm_i, :H2O, 0.0) * (2.01588 / 18.01528) +
            get(sol.M_atm_i, :CH4, 0.0) * (4.03176 / 16.04246) +
            get(sol.M_atm_i, :NH3, 0.0) * (3.02382 / 17.03052) +
            get(sol.M_atm_i, :H2S, 0.0) * (2.01588 / 34.08088)
        )
        M_melt_H = sol.M_melt_H
        @test isapprox(M_melt_H + M_atm_H, M_tot_H; rtol=1e-6)

        # 2. Carbon conservation
        M_atm_C = (
            get(sol.M_atm_i, :CO, 0.0) * (12.011 / 28.0101) +
            get(sol.M_atm_i, :CO2, 0.0) * (12.011 / 44.0095) +
            get(sol.M_atm_i, :CH4, 0.0) * (12.011 / 16.04246)
        )
        M_melt_C = sol.M_melt_C
        @test isapprox(M_melt_C + M_atm_C, M_tot_C; rtol=1e-6)

        # 3. Nitrogen conservation
        M_atm_N = (
            get(sol.M_atm_i, :N2, 0.0) + get(sol.M_atm_i, :NH3, 0.0) * (14.007 / 17.03052)
        )
        M_melt_N = sol.M_melt_N
        @test isapprox(M_melt_N + M_atm_N, M_tot_N; rtol=1e-6)

        # 4. Sulfur conservation
        M_atm_S = (
            get(sol.M_atm_i, :H2S, 0.0) * (32.06 / 34.08088) +
            get(sol.M_atm_i, :SO2, 0.0) * (32.06 / 64.066) +
            get(sol.M_atm_i, :S2, 0.0)
        )
        M_melt_S = sol.M_melt_S
        @test isapprox(M_melt_S + M_atm_S, M_tot_S; rtol=1e-6)

        # Monotonicity with total volatile inventory
        sol_2x = solve_magma_ocean_volatile_partitioning(
            M_melt,
            2.0 * M_tot_H,
            2.0 * M_tot_C,
            2.0 * M_tot_N,
            2.0 * M_tot_S,
            R_p,
            g_surf,
            T_mo,
            0.0;
            carbon_active=true,
            sulfur_active=true,
        )
        @test sol_2x.P_surf > sol.P_surf
        @test sol_2x.M_atm_tot > sol.M_atm_tot
        @test sol_2x.M_melt_tot > sol.M_melt_tot

        # Monotonicity with melt mass (more melt volume dissolves more volatiles, lowering surface pressure)
        sol_more_melt = solve_magma_ocean_volatile_partitioning(
            2.0 * M_melt,
            M_tot_H,
            M_tot_C,
            M_tot_N,
            M_tot_S,
            R_p,
            g_surf,
            T_mo,
            0.0;
            carbon_active=true,
            sulfur_active=true,
        )
        @test sol_more_melt.P_surf < sol.P_surf
        @test sol_more_melt.M_melt_tot > sol.M_melt_tot
    end

    # -------------------------------------------------------------------------
    # 2. Redox State Control on Degassed Volatile Speciation
    # -------------------------------------------------------------------------
    @testset "Redox Sensitivity on Outgassed Speciation" begin
        # Reduced case: delta_IW = -3.0 (enstatite/ordinary chondrite conditions)
        sol_red = solve_magma_ocean_volatile_partitioning(
            M_melt,
            M_tot_H,
            M_tot_C,
            M_tot_N,
            M_tot_S,
            R_p,
            g_surf,
            T_mo,
            -3.0;
            carbon_active=true,
            sulfur_active=true,
        )

        # Oxidized case: delta_IW = +2.0 (oxidized carbonaceous chondrite / terrestrial)
        sol_ox = solve_magma_ocean_volatile_partitioning(
            M_melt,
            M_tot_H,
            M_tot_C,
            M_tot_N,
            M_tot_S,
            R_p,
            g_surf,
            T_mo,
            2.0;
            carbon_active=true,
            sulfur_active=true,
        )

        # Reduced atmospheres must have higher CO/CO2 and H2/H2O ratios
        co_co2_red = sol_red.p_i[:CO] / max(sol_red.p_i[:CO2], 1e-30)
        co_co2_ox = sol_ox.p_i[:CO] / max(sol_ox.p_i[:CO2], 1e-30)
        @test co_co2_red > co_co2_ox

        h2_h2o_red = sol_red.p_i[:H2] / max(sol_red.p_i[:H2O], 1e-30)
        h2_h2o_ox = sol_ox.p_i[:H2] / max(sol_ox.p_i[:H2O], 1e-30)
        @test h2_h2o_red > h2_h2o_ox

        # Under reducing conditions, nitrogen dissolves chemically as nitride in melt
        # Hence, dissolved N in melt must be higher under reducing conditions
        @test sol_red.M_melt_N > sol_ox.M_melt_N
    end

    # -------------------------------------------------------------------------
    # 3. Asymptotic Limits & Boundary Cases
    # -------------------------------------------------------------------------
    @testset "Asymptotic Limits & Boundary Conditions" begin
        # Zero melt mass -> All volatiles in atmosphere
        sol_no_melt = solve_magma_ocean_volatile_partitioning(
            0.0,
            M_tot_H,
            M_tot_C,
            M_tot_N,
            M_tot_S,
            R_p,
            g_surf,
            T_mo,
            0.0;
            carbon_active=true,
            sulfur_active=true,
        )
        @test sol_no_melt.M_melt_tot ≈ 0.0 atol=1e-12
        @test isapprox(
            sol_no_melt.M_atm_tot, M_tot_H + M_tot_C + M_tot_N + M_tot_S; rtol=1e-4
        )

        # Zero volatile inventory -> Zero surface pressure
        sol_dry = solve_magma_ocean_volatile_partitioning(
            M_melt, 0.0, 0.0, 0.0, 0.0, R_p, g_surf, T_mo, 0.0
        )
        @test iszero(sol_dry.P_surf)
        @test iszero(sol_dry.M_atm_tot)

        # Error handling on invalid inputs
        @test_throws DomainError solve_magma_ocean_volatile_partitioning(
            -1.0, M_tot_H, M_tot_C, M_tot_N, M_tot_S, R_p, g_surf, T_mo, 0.0
        )
        @test_throws DomainError solve_magma_ocean_volatile_partitioning(
            M_melt, -1.0, M_tot_C, M_tot_N, M_tot_S, R_p, g_surf, T_mo, 0.0
        )
        @test_throws DomainError solve_magma_ocean_volatile_partitioning(
            M_melt, M_tot_H, M_tot_C, M_tot_N, M_tot_S, 0.0, g_surf, T_mo, 0.0
        )
        @test_throws DomainError solve_magma_ocean_volatile_partitioning(
            M_melt, M_tot_H, M_tot_C, M_tot_N, M_tot_S, R_p, 0.0, T_mo, 0.0
        )
        @test_throws DomainError solve_magma_ocean_volatile_partitioning(
            M_melt, M_tot_H, M_tot_C, M_tot_N, M_tot_S, R_p, g_surf, 0.0, 0.0
        )
    end

    # -------------------------------------------------------------------------
    # 4. Marker-Scale Magma Ocean Degassing & Crystallization Exsolution
    # -------------------------------------------------------------------------
    @testset "Lagrangian Marker Magma Ocean Degassing" begin
        # Setup mock marker arrays
        N_markers = 100
        xm = fill(0.0, N_markers)
        ym = fill(0.95 * R_p, N_markers) # Near-surface markers
        tm = fill(2, N_markers) # Rock markers
        tkm = fill(1800.0, N_markers) # Molten
        Fm = fill(0.8, N_markers) # 80% melt fraction
        Fm_old = fill(0.8, N_markers)

        # High volatile concentration in standard marker units:
        # XH2Om: wt% (2.0 wt%)
        # XCm: ppmw (1000.0 ppm)
        # XNm: ppmw (100.0 ppm)
        # XSm: ppmw (2000.0 ppm)
        XH2Om = fill(2.0, N_markers)
        XCm = fill(1000.0, N_markers)
        XNm = fill(100.0, N_markers)
        XSm = fill(2000.0, N_markers)
        rhosolidm_val = 3000.0
        marker_vol = ((2.0 * R_p)^2) / 10000.0 # Marker area in 2D

        cfg_degas = MagmaOceanDegassingConfig(;
            active=true,
            mode=:dynamic_flux,
            F_melt_threshold=0.40,
            degas_depth_fraction=0.90,
            efficiency=1.0,
        )

        dt = 1000.0 # 1000 s
        P_surf = 1.0e5 # 1 bar

        # Decompression degassing pass
        res_degas = degas_magma_ocean_markers!(
            xm,
            ym,
            tm,
            tkm,
            Fm,
            Fm_old,
            XH2Om,
            XCm,
            XNm,
            XSm,
            N_markers,
            dt,
            P_surf,
            R_p,
            cfg_degas;
            rho_solid=rhosolidm_val,
            marker_volume=marker_vol,
            delta_IW=0.0,
        )
        rates = res_degas.rates

        @test rates isa Dict{Symbol,Float64}
        @test haskey(rates, :H2O)
        @test rates[:H2O] > 0.0
        @test rates[:CO] > 0.0 || rates[:CO2] > 0.0
        @test rates[:N2] > 0.0
        @test rates[:H2S] > 0.0 || rates[:SO2] > 0.0 || rates[:S2] > 0.0
        @test res_degas.dM_2D isa Dict{Symbol,Float64}
        @test res_degas.dM_3D isa Dict{Symbol,Float64}
        @test res_degas.records isa Vector{TransferRecord}

        # Marker volatiles must have decreased due to degassing
        @test all(XH2Om .< 2.0)
        @test all(XCm .< 1000.0)
        @test all(XNm .< 100.0)
        @test all(XSm .< 2000.0)

        # Crystallization exsolution pass: melt fraction decreases from 0.8 to 0.2
        Fm_crystallizing = fill(0.2, N_markers)
        Fm_prev = fill(0.8, N_markers)
        XH2O_pre_cryst = copy(XH2Om)

        res_cryst = degas_magma_ocean_markers!(
            xm,
            ym,
            tm,
            tkm,
            Fm_crystallizing,
            Fm_prev,
            XH2Om,
            XCm,
            XNm,
            XSm,
            N_markers,
            dt,
            P_surf,
            R_p,
            cfg_degas;
            rho_solid=rhosolidm_val,
            marker_volume=marker_vol,
            delta_IW=0.0,
        )
        rates_cryst = res_cryst.rates

        @test rates_cryst[:H2O] > 0.0
        @test all(XH2Om .<= XH2O_pre_cryst)
    end

    # -------------------------------------------------------------------------
    # 5. Extreme Elemental Ratios Mass Conservation
    # -------------------------------------------------------------------------
    @testset "Extreme Elemental Ratios Mass Conservation" begin
        sol_ext = solve_magma_ocean_volatile_partitioning(
            M_melt, 1.0e10, 1.0e15, 1.0e8, 1.0e11, R_p, g_surf, T_mo, 3.0;
        )
        m_atm_N =
            get(sol_ext.M_atm_i, :N2, 0.0) +
            get(sol_ext.M_atm_i, :NH3, 0.0) * (14.007 / 17.03052)
        @test isapprox(m_atm_N + sol_ext.M_melt_N, 1.0e8; rtol=1e-5)
        m_atm_H =
            get(sol_ext.M_atm_i, :H2, 0.0) * 1.0 +
            get(sol_ext.M_atm_i, :H2O, 0.0) * (2.01588 / 18.01528) +
            get(sol_ext.M_atm_i, :CH4, 0.0) * (4.03176 / 16.04246) +
            get(sol_ext.M_atm_i, :NH3, 0.0) * (3.02382 / 17.03052) +
            get(sol_ext.M_atm_i, :H2S, 0.0) * (2.01588 / 34.08088)
        @test isapprox(m_atm_H + sol_ext.M_melt_H, 1.0e10; rtol=1e-5)
    end

    # -------------------------------------------------------------------------
    # 6. Magma Ocean Degassing Integration in simulation_loop
    # -------------------------------------------------------------------------
    @testset "Magma Ocean Degassing Runtime Integration in simulation_loop" begin
        quick_toml = joinpath(@__DIR__, "..", "configs", "test_quick.toml")
        cfg_base = load_config(quick_toml)

        for mode in (:dynamic_flux, :equilibrium)
            mktempdir() do output_dir
                cfg_run = SimulationConfig(;
                    grid=cfg_base.grid,
                    geometry=cfg_base.geometry,
                    time=TimeConfig(
                        dt_initial=cfg_base.time.dt_initial,
                        dt_longest=cfg_base.time.dt_longest,
                        dtcoefdn=cfg_base.time.dtcoefdn,
                        dtcoefup=cfg_base.time.dtcoefup,
                        dtstep=cfg_base.time.dtstep,
                        dxymax=cfg_base.time.dxymax,
                        vpratio=cfg_base.time.vpratio,
                        DTmax=cfg_base.time.DTmax,
                        start_time=cfg_base.time.start_time,
                        endtime=cfg_base.time.endtime,
                        start_step=1,
                        n_steps=2,
                    ),
                    solver=cfg_base.solver,
                    poroelasticity=cfg_base.poroelasticity,
                    thermodynamics=cfg_base.thermodynamics,
                    reaction=cfg_base.reaction,
                    materials=cfg_base.materials,
                    output=OutputConfig(output_dir=output_dir, savematstep=2),
                    disk=cfg_base.disk,
                    melting=cfg_base.melting,
                    venting=VentingConfig(active=false),
                    magma_degassing=MagmaOceanDegassingConfig(active=true, mode=mode),
                    atmosphere=AtmosphereConfig(active=true, mode=:guillot),
                )

                Erebus.simulation_loop(cfg_run; output_path=output_dir)

                files = readdir(output_dir)
                @test "output_00000.jld2" in files
                @test "output_00002.jld2" in files

                data2 = load_state(joinpath(output_dir, "output_00002.jld2"))
                @test data2["timestep"] == 2
                @test haskey(data2, "M_planet_val")
                @test data2["M_planet_val"] > 0.0
                @test haskey(data2, "atm_P_surf")
            end
        end
    end

    @testset "V1: Center-relative degassing coordinate invariance" begin
        cfg_degas = MagmaOceanDegassingConfig(;
            active=true,
            mode=:dynamic_flux,
            degas_depth_fraction=0.90,
            F_melt_threshold=0.40,
        )
        R_val = 50000.0
        xc_val = 70000.0
        yc_val = 70000.0
        n_m = 400
        th = range(0, 2π; length=n_m + 1)[1:n_m]

        # Ring at 0.95R around (70km, 70km)
        xm_c = xc_val .+ 0.95 * R_val .* cos.(th)
        ym_c = yc_val .+ 0.95 * R_val .* sin.(th)
        # Same ring around origin (0, 0)
        xm_o = 0.95 * R_val .* cos.(th)
        ym_o = 0.95 * R_val .* sin.(th)
        # Ring at 0.85R around (70km, 70km) (below degas_depth_fraction 0.90)
        xm_deep = xc_val .+ 0.85 * R_val .* cos.(th)
        ym_deep = yc_val .+ 0.85 * R_val .* sin.(th)

        tm_ring = fill(1, n_m)
        tkm_ring = fill(1900.0, n_m)
        Fm_ring = fill(1.0, n_m)
        Fm_old_ring = fill(1.0, n_m)

        mk_vols() = (fill(5.0, n_m), fill(2000.0, n_m), fill(500.0, n_m), fill(3000.0, n_m))

        w3d_c = [
            marker_out_of_plane_length(xm_c[m], ym_c[m], xc_val, yc_val) for m in 1:n_m
        ]
        w3d_o = [marker_out_of_plane_length(xm_o[m], ym_o[m], 0.0, 0.0) for m in 1:n_m]
        w3d_deep = [
            marker_out_of_plane_length(xm_deep[m], ym_deep[m], xc_val, yc_val) for
            m in 1:n_m
        ]

        XH_c, XC_c, XN_c, XS_c = mk_vols()
        res_c = degas_magma_ocean_markers!(
            xm_c,
            ym_c,
            tm_ring,
            tkm_ring,
            Fm_ring,
            Fm_old_ring,
            XH_c,
            XC_c,
            XN_c,
            XS_c,
            n_m,
            1.0e6,
            1.0e5,
            R_val,
            cfg_degas;
            xcenter=xc_val,
            ycenter=yc_val,
            w3d_m=w3d_c,
            rho_solid=3000.0,
            marker_volume=1.0e7,
        )

        XH_o, XC_o, XN_o, XS_o = mk_vols()
        res_o = degas_magma_ocean_markers!(
            xm_o,
            ym_o,
            tm_ring,
            tkm_ring,
            Fm_ring,
            Fm_old_ring,
            XH_o,
            XC_o,
            XN_o,
            XS_o,
            n_m,
            1.0e6,
            1.0e5,
            R_val,
            cfg_degas;
            xcenter=0.0,
            ycenter=0.0,
            w3d_m=w3d_o,
            rho_solid=3000.0,
            marker_volume=1.0e7,
        )

        XH_deep, XC_deep, XN_deep, XS_deep = mk_vols()
        res_deep = degas_magma_ocean_markers!(
            xm_deep,
            ym_deep,
            tm_ring,
            tkm_ring,
            Fm_ring,
            Fm_old_ring,
            XH_deep,
            XC_deep,
            XN_deep,
            XS_deep,
            n_m,
            1.0e6,
            1.0e5,
            R_val,
            cfg_degas;
            xcenter=xc_val,
            ycenter=yc_val,
            w3d_m=w3d_deep,
            rho_solid=3000.0,
            marker_volume=1.0e7,
        )

        # 1. Ring around (70km, 70km) degases at rate > 0
        @test any(!iszero, values(res_c.rates))
        # 2. Ring around (70km, 70km) and ring around origin give identical rates to 1e-12
        for sp in keys(res_c.rates)
            @test isapprox(res_c.rates[sp], res_o.rates[sp]; atol=1.0e-12, rtol=1.0e-12)
        end
        # 3. Ring at 0.85R (below degas_depth_fraction) returns zero rates
        @test all(iszero, values(res_deep.rates))
        @test isempty(res_deep.records)
    end

    @testset "PR 1b: V2 Crystallising marker melt-frame saturation" begin
        R_p = 50_000.0
        P_surf = 1.0e5
        T_melt = 1500.0
        d_IW = 0.0
        cfg = MagmaOceanDegassingConfig(;
            active=true,
            mode=:dynamic_flux,
            degas_depth_fraction=0.90,
            F_melt_threshold=0.40,
            water_As=0.40,
            efficiency=1.0,
        )

        spec = solve_chnos_speciation(P_surf, T_melt, d_IW)
        p_H2O_MPa = spec.p_H2O_Pa * 1.0e-6
        S_H2O = cfg.water_As * sqrt(p_H2O_MPa)
        S_N = compute_nitrogen_solubility_melt(spec.p_N2_Pa, d_IW).total_ppm
        S_C = compute_carbon_solubility_melt(P_surf, T_melt, d_IW).total_ppm
        S_S = compute_sulfur_solubility_melt(spec.p_S2_Pa, T_melt, d_IW)

        # 1. Crystallising marker: Fm_old = 0.8 -> Fm = 0.5
        xm = [0.0]
        ym = [0.95 * R_p]
        tm = [2]
        tkm = [T_melt]
        Fm = [0.5]
        Fm_old = [0.8]
        XH = [2.0 * S_H2O]
        XC = [2.0 * S_C]
        XN = [2.0 * S_N]
        XS = [2.0 * S_S]

        degas_magma_ocean_markers!(
            xm,
            ym,
            tm,
            tkm,
            Fm,
            Fm_old,
            XH,
            XC,
            XN,
            XS,
            1,
            1000.0,
            P_surf,
            R_p,
            cfg,
            T_melt;
            marker_volume=1.0,
        )

        @test isapprox(XH[1] / (Fm[1] * S_H2O), 1.0; rtol=1e-6)
        @test isapprox(XC[1] / (Fm[1] * S_C), 1.0; rtol=1e-6)
        @test isapprox(XN[1] / (Fm[1] * S_N), 1.0; rtol=1e-6)
        @test isapprox(XS[1] / (Fm[1] * S_S), 1.0; rtol=1e-6)

        # 2. Crystallising marker: Fm_old = 0.8 -> Fm = 0.25
        Fm2 = [0.25]
        Fm_old2 = [0.8]
        XH2 = [2.0 * S_H2O]
        XC2 = [2.0 * S_C]
        XN2 = [2.0 * S_N]
        XS2 = [2.0 * S_S]

        degas_magma_ocean_markers!(
            xm,
            ym,
            tm,
            tkm,
            Fm2,
            Fm_old2,
            XH2,
            XC2,
            XN2,
            XS2,
            1,
            1000.0,
            P_surf,
            R_p,
            cfg,
            T_melt;
            marker_volume=1.0,
        )

        @test isapprox(XH2[1] / (Fm2[1] * S_H2O), 1.0; rtol=1e-6)
        @test isapprox(XC2[1] / (Fm2[1] * S_C), 1.0; rtol=1e-6)
        @test isapprox(XN2[1] / (Fm2[1] * S_N), 1.0; rtol=1e-6)
        @test isapprox(XS2[1] / (Fm2[1] * S_S), 1.0; rtol=1e-6)
    end

    @testset "PR 1b: Steady supersaturated marker" begin
        R_p = 50_000.0
        P_surf = 1.0e5
        T_melt = 1500.0
        d_IW = 0.0
        cfg = MagmaOceanDegassingConfig(;
            active=true,
            mode=:dynamic_flux,
            degas_depth_fraction=0.90,
            F_melt_threshold=0.40,
            water_As=0.40,
            efficiency=1.0,
        )

        spec = solve_chnos_speciation(P_surf, T_melt, d_IW)
        S_H2O = cfg.water_As * sqrt(spec.p_H2O_Pa * 1.0e-6)

        # 1. Steady fully molten marker: Fm = Fm_old = 1.0 extracts down to S
        xm = [0.0]
        ym = [0.95 * R_p]
        tm = [2]
        tkm = [T_melt]
        Fm_1 = [1.0]
        Fm_old_1 = [1.0]
        XH_1 = [2.5 * S_H2O]
        XC_1 = [0.0]
        XN_1 = [0.0]
        XS_1 = [0.0]

        degas_magma_ocean_markers!(
            xm,
            ym,
            tm,
            tkm,
            Fm_1,
            Fm_old_1,
            XH_1,
            XC_1,
            XN_1,
            XS_1,
            1,
            1000.0,
            P_surf,
            R_p,
            cfg,
            T_melt;
            marker_volume=1.0,
        )
        @test isapprox(XH_1[1], S_H2O; rtol=1e-6)

        # 2. Steady partially molten marker: Fm = Fm_old = 0.5 with bulk above F*S extracts down to F*S
        Fm_half = [0.5]
        Fm_old_half = [0.5]
        XH_half = [0.75 * S_H2O]
        XC_half = [0.0]
        XN_half = [0.0]
        XS_half = [0.0]

        degas_magma_ocean_markers!(
            xm,
            ym,
            tm,
            tkm,
            Fm_half,
            Fm_old_half,
            XH_half,
            XC_half,
            XN_half,
            XS_half,
            1,
            1000.0,
            P_surf,
            R_p,
            cfg,
            T_melt;
            marker_volume=1.0,
        )
        @test isapprox(XH_half[1], 0.5 * S_H2O; rtol=1e-6)
    end

    @testset "PR 1b: Physical bounds and skip invariants" begin
        R_p = 50_000.0
        P_surf = 1.0e5
        T_melt = 1500.0
        cfg = MagmaOceanDegassingConfig(;
            active=true,
            mode=:dynamic_flux,
            degas_depth_fraction=0.90,
            F_melt_threshold=0.40,
            water_As=0.40,
            efficiency=1.0,
        )
        spec = solve_chnos_speciation(P_surf, T_melt, 0.0)
        S_H2O = cfg.water_As * sqrt(spec.p_H2O_Pa * 1.0e-6)

        # Extraction never exceeds bulk inventory for F in {0.011, 0.5, 1.0}
        for F_val in (0.011, 0.5, 1.0)
            xm = [0.0]
            ym = [0.95 * R_p]
            tm = [2]
            tkm = [T_melt]
            Fm = [F_val]
            Fm_old = [F_val]
            X_init = 5.0
            XH = [X_init]
            XC = [5000.0]
            XN = [500.0]
            XS = [5000.0]

            degas_magma_ocean_markers!(
                xm,
                ym,
                tm,
                tkm,
                Fm,
                Fm_old,
                XH,
                XC,
                XN,
                XS,
                1,
                1000.0,
                P_surf,
                R_p,
                cfg,
                T_melt;
                marker_volume=1.0,
            )
            @test 0.0 <= XH[1] <= X_init
            @test isapprox(XH[1], F_val * S_H2O; rtol=1e-6)
        end

        # F = 0 markers are skipped
        xm_zero = [0.0]
        ym_zero = [0.95 * R_p]
        tm_zero = [2]
        tkm_zero = [T_melt]
        Fm_zero = [0.0]
        Fm_old_zero = [0.0]
        XH_zero = [2.0]
        XC_zero = [1000.0]
        XN_zero = [100.0]
        XS_zero = [1000.0]
        res_zero = degas_magma_ocean_markers!(
            xm_zero,
            ym_zero,
            tm_zero,
            tkm_zero,
            Fm_zero,
            Fm_old_zero,
            XH_zero,
            XC_zero,
            XN_zero,
            XS_zero,
            1,
            1000.0,
            P_surf,
            R_p,
            cfg,
            T_melt;
            marker_volume=1.0,
        )
        @test isapprox(XH_zero[1], 2.0; rtol=1e-12)
        @test all(iszero, values(res_zero.rates))

        # Marker below saturation extracts nothing
        Fm_sub = [0.5]
        Fm_old_sub = [0.5]
        XH_sub = [0.25 * S_H2O]
        XC_sub = [0.0]
        XN_sub = [0.0]
        XS_sub = [0.0]
        res_sub = degas_magma_ocean_markers!(
            xm_zero,
            ym_zero,
            tm_zero,
            tkm_zero,
            Fm_sub,
            Fm_old_sub,
            XH_sub,
            XC_sub,
            XN_sub,
            XS_sub,
            1,
            1000.0,
            P_surf,
            R_p,
            cfg,
            T_melt;
            marker_volume=1.0,
        )
        @test isapprox(XH_sub[1], 0.25 * S_H2O; rtol=1e-12)
        @test all(iszero, values(res_sub.rates))
    end

    @testset "PR 1b: V7 partial_pressures_to_masses and mu_bar guard" begin
        amu = Erebus.SPECIES_AMU
        col_c = 1.0e10

        # All-zero partial pressures returns finite zeros (no NaN)
        p_zeros = Dict{Symbol,Float64}(:H2O => 0.0, :CO2 => 0.0, :N2 => 0.0)
        m_zero_p0 = Erebus.partial_pressures_to_masses(p_zeros, 0.0, col_c, amu)
        @test all(iszero, values(m_zero_p0))
        @test all(isfinite, values(m_zero_p0))

        m_zero_p_pos = Erebus.partial_pressures_to_masses(p_zeros, 1.0e5, col_c, amu)
        @test all(iszero, values(m_zero_p_pos))
        @test all(isfinite, values(m_zero_p_pos))

        # Negative P_total throws DomainError
        p_normal = Dict{Symbol,Float64}(:H2O => 1.0e4, :CO2 => 9.0e4)
        @test_throws DomainError Erebus.partial_pressures_to_masses(
            p_normal, -1.0, col_c, amu
        )

        # Negative species partial pressure throws DomainError
        p_neg = Dict{Symbol,Float64}(:H2O => -1.0, :CO2 => 1.0e5)
        @test_throws DomainError Erebus.partial_pressures_to_masses(
            p_neg, 1.0e5, col_c, amu
        )
    end

    @testset "PR 1b: Retention floor enforcement with linear melt blend" begin
        R_p = 50_000.0
        P_surf = 1.0e5
        T_melt = 1500.0
        cfg = MagmaOceanDegassingConfig(;
            active=true,
            mode=:dynamic_flux,
            degas_depth_fraction=0.90,
            F_melt_threshold=0.40,
            water_As=0.40,
            efficiency=1.0,
        )
        ret_cfg = RetentionConfig(;
            active=true,
            retention_law=:linear_melt_blend,
            h2o_retention_ppm=1000.0,
            carbon_retention_ppm=500.0,
            nitrogen_retention_ppm=50.0,
            sulfur_retention_ppm=500.0,
        )

        # Marker with F = 0.2: retention floor is 1000 * (1 - 0.2) = 800 ppmw = 0.08 wt%
        xm = [0.0]
        ym = [0.95 * R_p]
        tm = [2]
        tkm = [T_melt]
        Fm = [0.2]
        Fm_old = [0.2]
        XH = [0.20] # 2000 ppm bulk water
        XC = [1000.0]
        XN = [100.0]
        XS = [1000.0]

        degas_magma_ocean_markers!(
            xm,
            ym,
            tm,
            tkm,
            Fm,
            Fm_old,
            XH,
            XC,
            XN,
            XS,
            1,
            1000.0,
            P_surf,
            R_p,
            cfg,
            T_melt;
            marker_volume=1.0,
            retention_cfg=ret_cfg,
        )

        expected_ret_floor_wtpct = 1000.0 * (1.0 - 0.2) * 1.0e-4 # 0.08 wt%
        expected_ret_c_ppm = 500.0 * (1.0 - 0.2) # 400 ppm
        expected_ret_n_ppm = 50.0 * (1.0 - 0.2)  # 40 ppm
        expected_ret_s_ppm = 500.0 * (1.0 - 0.2) # 400 ppm

        @test isapprox(XH[1], expected_ret_floor_wtpct; rtol=1e-6)
        @test isapprox(XC[1], expected_ret_c_ppm; rtol=1e-6)
        @test isapprox(XN[1], expected_ret_n_ppm; rtol=1e-6)
        @test isapprox(XS[1], expected_ret_s_ppm; rtol=1e-6)
    end

    @testset "PR 1b: Reference temperature discrimination" begin
        R_p = 50_000.0
        P_surf = 1.0e5
        d_IW = 0.0
        cfg = MagmaOceanDegassingConfig(;
            active=true,
            mode=:dynamic_flux,
            degas_depth_fraction=0.90,
            F_melt_threshold=0.40,
            water_As=0.40,
            efficiency=1.0,
        )

        spec_1500 = solve_chnos_speciation(P_surf, 1500.0, d_IW)
        S_1500 = cfg.water_As * sqrt(spec_1500.p_H2O_Pa * 1.0e-6)

        spec_1200 = solve_chnos_speciation(P_surf, 1200.0, d_IW)
        S_1200 = cfg.water_As * sqrt(spec_1200.p_H2O_Pa * 1.0e-6)

        @test !isapprox(S_1500, S_1200; rtol=1e-3)

        # Test at 1500 K
        xm = [0.0]
        ym = [0.95 * R_p]
        tm = [2]
        tkm = [1800.0]
        Fm1 = [0.5]
        Fm_old1 = [0.5]
        XH1 = [2.0 * S_1500]
        XC1 = [0.0]
        XN1 = [0.0]
        XS1 = [0.0]
        degas_magma_ocean_markers!(
            xm,
            ym,
            tm,
            tkm,
            Fm1,
            Fm_old1,
            XH1,
            XC1,
            XN1,
            XS1,
            1,
            1000.0,
            P_surf,
            R_p,
            cfg,
            1500.0;
            marker_volume=1.0,
        )
        @test isapprox(XH1[1], Fm1[1] * S_1500; rtol=1e-6)

        # Test at 1200 K
        Fm2 = [0.5]
        Fm_old2 = [0.5]
        XH2 = [2.0 * S_1200]
        XC2 = [0.0]
        XN2 = [0.0]
        XS2 = [0.0]
        degas_magma_ocean_markers!(
            xm,
            ym,
            tm,
            tkm,
            Fm2,
            Fm_old2,
            XH2,
            XC2,
            XN2,
            XS2,
            1,
            1000.0,
            P_surf,
            R_p,
            cfg,
            1200.0;
            marker_volume=1.0,
        )
        @test isapprox(XH2[1], Fm2[1] * S_1200; rtol=1e-6)
    end
end

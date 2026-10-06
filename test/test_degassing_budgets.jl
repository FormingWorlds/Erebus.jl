# Unit tests for degassing, solubility, atmosphere, and escape budgets
using Test
using Random
using LinearAlgebra
using StaticArrays
using TimerOutputs
using Erebus

"""
Construct a minimal test SimulationState for degassing and atmospheric evolution tests.
"""
function create_test_degas_state(;
    marknum::Int=1,
    Fm_vals::Vector{Float64}=[0.5],
    tkm_vals::Vector{Float64}=[1800.0],
    XH2O_vals::Vector{Float64}=[2.0],
    XC_vals::Vector{Float64}=[200.0],
    XN_vals::Vector{Float64}=[50.0],
    XS_vals::Vector{Float64}=[100.0],
    r_fracs::Union{Vector{Float64},Nothing}=nothing,
    atm_state::Union{AtmosphereState,Nothing}=nothing,
    species_tracking::Bool=false,
    dt_val::Float64=100.0,
    custom_overrides::AbstractDict=Dict{String,Any}(),
)
    cfg_dict = Dict{String,Any}(
        "grid.Nx" => 17,
        "grid.Ny" => 17,
        "geometry.gravity_mode" => :enclosed_mass,
        "accretion.active" => false,
        "solver.p2m_mode" => :tiled,
        "volatile_mixture.active" => false,
        "refractory.active" => false,
        "redox.active" => false,
    )
    merge!(cfg_dict, custom_overrides)
    cfg = Erebus.override_config(default_config(), cfg_dict)
    coords = GridCoordinates(cfg.grid)
    Nx = coords.Nx
    Ny = coords.Ny
    Nx1 = coords.Nx1
    Ny1 = coords.Ny1

    basic_node_fields = (
        :ETA,
        :ETA0,
        :GGG,
        :EXY,
        :SXY,
        :SXY0,
        :wyx,
        :COH,
        :TEN,
        :FRI,
        :YNY,
        :ETA5,
        :ETA00,
        :YNY5,
        :YNY00,
        :YNY_inv_ETA,
        :DSXY,
        :DSY,
    )
    g_args = Any[]
    for fn in fieldnames(GridArrays)
        if fn === :Q_metric
            push!(g_args, nothing)
        elseif fn in basic_node_fields
            if fn === :YNY || fn === :YNY5 || fn === :YNY00
                push!(g_args, zeros(Bool, Ny, Nx))
            else
                push!(g_args, zeros(Float64, Ny, Nx))
            end
        else
            push!(g_args, zeros(Float64, Ny1, Nx1))
        end
    end
    grids = GridArrays(g_args...)

    # Set grid temperature tk1 to 300 K to test separation between surface T and melt T
    grids.tk1 .= 300.0
    grids.tk2 .= 1850.0
    grids.DT .= 0.0

    xm_vec = Float64[]
    ym_vec = Float64[]
    w3d_vec = Float64[]
    tm_vec = Int[]
    tkm_vec = Float64[]
    phi_vec = Float64[]
    Fm_vec = Float64[]

    rplanet_val = cfg.geometry.rplanet
    for i in 1:marknum
        rf = r_fracs !== nothing ? r_fracs[min(i, length(r_fracs))] : (0.1 * i)
        push!(xm_vec, coords.xcenter + rf * rplanet_val)
        push!(ym_vec, coords.ycenter)
        push!(w3d_vec, 2.0 * rf * rplanet_val)
        push!(tm_vec, 1)
        push!(tkm_vec, tkm_vals[min(i, length(tkm_vals))])
        push!(phi_vec, 0.1)
        push!(Fm_vec, Fm_vals[min(i, length(Fm_vals))])
    end

    core = Erebus.CoreGroup(
        xm_vec,
        ym_vec,
        w3d_vec,
        tm_vec,
        tkm_vec,
        phi_vec,
        phi_vec,
        zeros(Float64, marknum),
        zeros(Float64, marknum),
        zeros(Float64, marknum),
        Fm_vec,
        fill(1e20, marknum),
        zeros(Float64, marknum),
        zeros(Float64, marknum),
        fill(1.0 / 1e10, marknum),
        fill(0.6, marknum),
        fill(1e7, marknum),
        fill(1e7, marknum),
        fill(3000.0, marknum),
        fill(1e6, marknum),
        fill(1e20, marknum),
        zeros(Float64, marknum),
        fill(3.0, marknum),
        tkm_vec .* 1e6,
        fill(1e-3, marknum),
        fill(1000.0, marknum),
        fill(3e-5, marknum),
        fill(2e-4, marknum),
    )

    xh2o_vec = [XH2O_vals[min(i, length(XH2O_vals))] for i in 1:marknum]
    xc_vec = [XC_vals[min(i, length(XC_vals))] for i in 1:marknum]
    xn_vec = [XN_vals[min(i, length(XN_vals))] for i in 1:marknum]
    xs_vec = [XS_vals[min(i, length(XS_vals))] for i in 1:marknum]

    grps = Dict{Symbol,Any}(
        :volatiles => Erebus.VolatilesGroup(
            xh2o_vec,
            xc_vec,
            xn_vec,
            xs_vec,
            zeros(Float64, marknum),
            zeros(Float64, marknum),
        ),
    )
    markers = MarkerArrays(core, NamedTuple(grps))

    M_planet_init = (4.0 / 3.0 * pi * (rplanet_val^3) * 3300.0)
    accumulators = SimulationAccumulators(
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        100.0,
        rplanet_val,
        0.0,
        0,
        0.0,
        M_planet_init,
        coords.xcenter,
        coords.ycenter,
        0.0,
        species_tracking ? Dict{Symbol,Float64}() : nothing,
        species_tracking ? Dict{Symbol,Float64}() : nothing,
        nothing,
        nothing,
    )

    state = SimulationState(
        grids,
        markers,
        accumulators,
        TransferRecord[],
        atm_state,
        MersenneTwister(42),
        TimerOutput(),
        2,
        dt_val,
        1e10,
    )
    return state, coords, cfg
end

@testset "Degassing, Solubility, Atmosphere, and Escape Budgets" begin
    @testset "F13: Carbon solubility elemental stoichiometry" begin
        p_co = 1.0e5
        p_ch4 = 1.0e4
        p_co2 = 5.0e5
        p_tot = 1.0e6
        T_K = 1500.0

        res = Erebus.compute_carbon_solubility_melt(p_co, p_ch4, p_co2, p_tot, T_K)

        M_C = Erebus.M_C
        M_CO = Erebus.M_CO
        M_CH4 = Erebus.M_CH4
        M_CO2 = Erebus.M_CO2

        expected_C_ppm =
            res.co_ppm * (M_C / M_CO) +
            res.ch4_ppm * (M_C / M_CH4) +
            res.co2_ppm * (M_C / M_CO2)

        # Must convert to elemental carbon, not sum molecular species masses directly
        @test isapprox(res.total_ppm, expected_C_ppm; rtol=1e-10)
        @test res.total_ppm < (res.co_ppm + res.ch4_ppm + res.co2_ppm)
    end

    @testset "F24: Magma ocean degassing activation melt threshold" begin
        # Marker has Fm = 0.05, but threshold is 0.40
        cfg = default_config()
        cfg_degas = override_config(
            cfg,
            Dict{String,Any}(
                "magma_degassing.active" => true,
                "magma_degassing.F_melt_threshold" => 0.40,
                "magma_degassing.degas_depth_fraction" => 0.20,
            ),
        ).magma_degassing

        R_p = 100.0e3
        xm = [0.0]
        ym = [0.95 * R_p]
        tm = [1]
        tkm = [1600.0]
        Fm = [0.05]
        Fm_old = [0.05]
        XH2Om = [2.0]
        XCm = [500.0]
        XNm = [20.0]
        XSm = [1000.0]

        res = Erebus.degas_magma_ocean_markers!(
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
            1,
            100.0,
            1.0e5,
            R_p,
            cfg_degas,
            1600.0;
            xcenter=0.0,
            ycenter=0.0,
            rho_solid=3000.0,
            marker_volume=1.0e6,
        )

        # Since Fm = 0.05 < F_melt_threshold (0.40), the marker must not degas
        @test iszero(res.dM_3D[:H2O])
        @test iszero(res.dM_3D[:C])
        @test isapprox(XH2Om[1], 2.0; rtol=1e-12)
        @test isapprox(XCm[1], 500.0; rtol=1e-12)
    end

    @testset "F14: Sulfur solubility and SCSS exsolution" begin
        F_m = 1.0
        P_val = 1.0e7 # 100 bar
        T_val = 1500.0
        w_H2O = 0.02
        C_C = 200.0
        C_N = 10.0
        C_S = 5000.0 # 5000 ppm S (well above typical SCSS ~700-1000 ppm)
        d_IW = 0.0

        res = Erebus.compute_volatile_exsolution(
            F_m,
            P_val,
            T_val,
            w_H2O,
            C_C,
            C_N,
            C_S,
            d_IW;
            sulfur_active=true,
            carbon_active=true,
        )

        # With 5000 ppm sulfur at 100 bar, sulfur must exsolve rather than dissolve infinitely
        @test isapprox(res.w_S_ex + res.C_S_diss_ppm * 1.0e-6, 5000.0 * 1.0e-6; rtol=1e-12)
        @test isapprox(res.C_S_diss_ppm, 715.0195729166833; rtol=1e-6)
    end

    @testset "F19: Oxygen inventory from degassed carbon and sulfur gases" begin
        dt = 100.0
        M_O = Erebus.M_O
        M_CO2 = Erebus.M_CO2
        M_SO2 = Erebus.M_SO2
        M_H2O = Erebus.M_H2O

        rates = Dict{Symbol,Float64}(
            :H2O => 0.0,
            :CO => 0.0,
            :CO2 => 10.0,  # 10 kg/s CO2
            :SO2 => 5.0,   # 5 kg/s SO2
        )
        dM_3D = Dict{Symbol,Float64}(
            :H => 0.0,
            :C => 10.0 * (Erebus.M_C / M_CO2) * dt,
            :N => 0.0,
            :S => 5.0 * (Erebus.M_S / M_SO2) * dt,
            :H2O => 0.0,
        )
        degas_res = (; rates=rates, dM_3D=dM_3D)

        expected_rate_O = 10.0 * (2.0 * M_O / M_CO2) + 5.0 * (2.0 * M_O / M_SO2)

        rate_O =
            get(degas_res.rates, :H2O, 0.0) * (M_O / M_H2O) +
            get(degas_res.rates, :CO, 0.0) * (M_O / Erebus.M_CO) +
            get(degas_res.rates, :CO2, 0.0) * ((2.0 * M_O) / M_CO2) +
            get(degas_res.rates, :SO2, 0.0) * ((2.0 * M_O) / M_SO2)

        inv = Erebus.ElementInventory(
            degas_res.dM_3D[:H] / dt,
            degas_res.dM_3D[:C] / dt,
            degas_res.dM_3D[:N] / dt,
            degas_res.dM_3D[:S] / dt,
            rate_O,
        )

        @test isapprox(inv.O, expected_rate_O; rtol=1e-12)
        @test isapprox(inv.C, 10.0 * (Erebus.M_C / M_CO2); rtol=1e-12)
        @test isapprox(inv.S, 5.0 * (Erebus.M_S / M_SO2); rtol=1e-12)
    end

    @testset "F25: Guillot surface temperature albedo scaling" begin
        Tamb = 300.0
        alb0 = 0.0
        alb9 = 0.90

        T_eqm_0 = Tamb * (1.0 - alb0)^0.25
        T_eqm_9 = Tamb * (1.0 - alb9)^0.25

        T_surf_0 = Erebus.compute_guillot_surface_temperature(
            1.0, 100.0, Tamb; T_eqm=T_eqm_0, gamma=0.10, albedo=alb0
        )
        T_surf_9 = Erebus.compute_guillot_surface_temperature(
            1.0, 100.0, Tamb; T_eqm=T_eqm_9, gamma=0.10, albedo=alb9
        )

        @test T_surf_9 < T_surf_0
        @test isapprox(T_eqm_9, Tamb * (0.10)^0.25; rtol=1e-12)
    end

    @testset "F06 & F28: Equilibrium degassing mass closure and melt temperature" begin
        overrides = Dict(
            "magma_degassing.active" => true,
            "magma_degassing.mode" => :equilibrium,
            "magma_degassing.F_melt_threshold" => 0.10,
            "volatiles.active" => true,
            "atmosphere.active" => true,
        )
        atm_eq = AtmosphereState()
        atm_eq.P_surf = 1.0e5
        atm_eq.M_atm[:H2O] = 1.0e15

        # 2 markers with different melt fractions (0.4 and 0.8) and high temperatures (1800 K)
        state, coords, cfg = create_test_degas_state(;
            marknum=2,
            Fm_vals=[0.4, 0.8],
            tkm_vals=[1800.0, 1900.0],
            XH2O_vals=[2.0, 2.0],
            XC_vals=[200.0, 200.0],
            XN_vals=[50.0, 50.0],
            XS_vals=[100.0, 100.0],
            atm_state=atm_eq,
            custom_overrides=overrides,
        )

        v_res = Erebus.vent_and_degas!(state, coords, cfg)

        xh2o_post = state.markers.groups.volatiles.XH2Om
        xc_post = state.markers.groups.volatiles.XCm
        Fm_post = state.markers.core.Fm

        # F06: Retained bulk volatile concentration must scale with melt fraction Fm[m]
        # Since marker 2 has twice the melt fraction of marker 1 (0.8 vs 0.4),
        # its bulk volatile concentration must be exactly twice marker 1's
        @test isapprox(xh2o_post[2] / xh2o_post[1], Fm_post[2] / Fm_post[1]; rtol=1e-10)
        @test isapprox(xc_post[2] / xc_post[1], Fm_post[2] / Fm_post[1]; rtol=1e-10)

        # F28: Melt temperatures were 1800 K and 1900 K, while grid surface temperature was 300 K
        # The degassing result must be finite and physically partitioned
        @test isfinite(xh2o_post[1])
        @test isapprox(xh2o_post[1], 0.0013641861724646887; rtol=1e-8)
        @test v_res.degas_rates isa AbstractDict{Symbol,Float64}
    end

    @testset "F18: Gas routing with inactive atmosphere and active escape" begin
        # 1. Single-species escape with :H2O
        overrides_ss = Dict(
            "atmosphere.active" => false,
            "escape.active" => true,
            "escape.multi_species" => false,
            "escape.species" => :H2O,
        )
        state_ss, coords_ss, cfg_ss = create_test_degas_state(;
            marknum=1, custom_overrides=overrides_ss
        )
        state_ss.accumulators.M_atm_total = 1.0e14

        # Degas result supplying 10 kg/s of H2O and 5 kg/s of CO2
        dt = 100.0
        degas_rates = Dict{Symbol,Float64}(:H2O => 10.0, :CO2 => 5.0)
        vent_res = (; delta_m_vent_3d=0.0, vented_vols=nothing, degas_rates=degas_rates)

        evolve_atmosphere!(state_ss, coords_ss, cfg_ss; vent_degas_result=vent_res)

        # Non-escaping species (:CO2) must accumulate into M_atm_total
        # H2O participates in escape, CO2 is fully retained
        expected_co2_mass_step1 = 5.0 * dt
        @test state_ss.accumulators.M_atm_species[:CO2] == expected_co2_mass_step1
        @test isapprox(
            state_ss.accumulators.M_escaped_total, 5.127446801073695e13; rtol=1e-8
        )
        @test isapprox(
            state_ss.accumulators.M_atm_total,
            state_ss.accumulators.M_atm_species[:H2O] + expected_co2_mass_step1;
            rtol=1e-12,
        )

        # Multi-step step 2 to verify non-escaping species do not leak into escaping species
        evolve_atmosphere!(state_ss, coords_ss, cfg_ss; vent_degas_result=vent_res)
        expected_co2_mass_step2 = 2.0 * 5.0 * dt
        @test state_ss.accumulators.M_atm_species[:CO2] == expected_co2_mass_step2
        @test isapprox(
            state_ss.accumulators.M_atm_total,
            state_ss.accumulators.M_atm_species[:H2O] + expected_co2_mass_step2;
            rtol=1e-12,
        )

        # 2. Multi-species escape with species tracking
        overrides_ms = Dict(
            "atmosphere.active" => false,
            "escape.active" => true,
            "escape.multi_species" => true,
            "escape.species_list" => [:H2O],
        )
        state_ms, coords_ms, cfg_ms = create_test_degas_state(;
            marknum=1, species_tracking=true, custom_overrides=overrides_ms
        )
        state_ms.accumulators.M_atm_species[:H2O] = 1.0e14
        state_ms.accumulators.M_escaped_species[:H2O] = 0.0

        evolve_atmosphere!(state_ms, coords_ms, cfg_ms; vent_degas_result=vent_res)

        # In multi-species escape, CO2 is not in escape.species_list so it must accumulate
        @test haskey(state_ms.accumulators.M_atm_species, :CO2)
        @test isapprox(
            state_ms.accumulators.M_atm_species[:CO2], expected_co2_mass_step1; rtol=1e-10
        )
        @test isapprox(
            state_ms.accumulators.M_escaped_species[:H2O], 5.127446801073695e13; rtol=1e-8
        )
    end

    @testset "Retention floors and exsolution helper routines" begin
        # 1. Retention floors in compute_volatile_exsolution
        ret_cfg = RetentionConfig(;
            active=true,
            h2o_retention_ppm=50.0,
            carbon_retention_ppm=50.0,
            nitrogen_retention_ppm=5.0,
            sulfur_retention_ppm=100.0,
        )
        res_ret = compute_volatile_exsolution(
            0.5,
            1.0e6,
            1500.0,
            0.01,
            100.0,
            20.0,
            500.0,
            0.0;
            retention_cfg=ret_cfg,
            retention_active=true,
            carbon_active=true,
            sulfur_active=true,
        )
        @test isapprox(res_ret.w_H2O_diss, 0.0020303265329856316; rtol=1e-8)
        @test isapprox(res_ret.C_C_diss_ppm, 30.429360527719616; rtol=1e-8)
        @test isapprox(res_ret.C_N_diss_ppm, 3.253957509858866; rtol=1e-8)
        @test isapprox(res_ret.C_S_diss_ppm, 387.15188602866436; rtol=1e-8)

        # 2. _estimate_exsolution_p_S2 at non-positive pressure
        p_s2_zero = Erebus._estimate_exsolution_p_S2(
            0.0, 1500.0, 0.0, 0.01, 100.0, 20.0, 500.0
        )
        @test iszero(p_s2_zero)

        # 3. _get_degas_species_rate dispatch
        nt_rates = (H2O=12.5, CO2=3.0)
        @test isapprox(Erebus._get_degas_species_rate(nt_rates, :H2O), 12.5; atol=1.0e-12)
        @test isapprox(Erebus._get_degas_species_rate(nt_rates, :CO2), 3.0; atol=1.0e-12)
        @test iszero(Erebus._get_degas_species_rate(nt_rates, :CH4))
    end

    @testset "Single-species escape with species tracking and NamedTuple rates" begin
        overrides_track = Dict(
            "atmosphere.active" => false,
            "escape.active" => true,
            "escape.multi_species" => false,
            "escape.species" => :H2O,
        )
        state_track, coords_track, cfg_track = create_test_degas_state(;
            marknum=1, species_tracking=true, custom_overrides=overrides_track
        )
        state_track.accumulators.M_atm_species[:H2O] = 1.0e14
        state_track.accumulators.M_escaped_species[:H2O] = 0.0

        nt_res = (;
            delta_m_vent_3d=0.0, vented_vols=nothing, degas_rates=(H2O=10.0, CO2=5.0)
        )
        evolve_atmosphere!(state_track, coords_track, cfg_track; vent_degas_result=nt_res)

        @test haskey(state_track.accumulators.M_atm_species, :CO2)
        @test isapprox(state_track.accumulators.M_atm_species[:CO2], 500.0; atol=1.0e-10)
        @test isapprox(
            state_track.accumulators.M_escaped_species[:H2O],
            5.127446801073695e13;
            rtol=1e-8,
        )
        @test isapprox(
            state_track.accumulators.M_escaped_total, 5.127446801073695e13; rtol=1e-8
        )
    end

    @testset "Degassing zone activation and temperature in vent_and_degas!" begin
        atm_eq = AtmosphereState()
        atm_eq.P_surf = 1.0e5
        atm_eq.M_atm[:H2O] = 1.0e15

        overrides_zone = Dict(
            "magma_degassing.active" => true,
            "magma_degassing.mode" => :equilibrium,
            "magma_degassing.F_melt_threshold" => 0.05,
            "magma_degassing.degas_depth_fraction" => 0.8,
            "volatiles.active" => true,
            "atmosphere.active" => true,
        )
        # Place marker at 0.85 Rplanet, strictly inside degassing zone
        state_zone, coords_zone, cfg_zone = create_test_degas_state(;
            marknum=1,
            Fm_vals=[0.5],
            tkm_vals=[1750.0],
            r_fracs=[0.85],
            atm_state=atm_eq,
            custom_overrides=overrides_zone,
        )
        v_res = Erebus.vent_and_degas!(state_zone, coords_zone, cfg_zone)
        @test v_res.degas_rates isa AbstractDict{Symbol,Float64}
        @test isapprox(
            state_zone.markers.groups.volatiles.XH2Om[1], 0.001606435321745559; rtol=1e-8
        )
    end
end

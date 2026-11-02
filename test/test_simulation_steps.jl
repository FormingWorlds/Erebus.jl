using Test
using Random
using LinearAlgebra
using StaticArrays
using TimerOutputs
using JLD2
using Erebus

"""
Construct a minimal test SimulationState with 0 or 1 marker.
"""
function create_mock_simulation_state(;
    marknum::Int=0,
    r_marker::Float64=0.0,
    tm_val::Int=1,
    phi_val::Float64=0.1,
    tkm_val::Float64=300.0,
    gravity_mode::Symbol=:enclosed_mass,
    accretion_active::Bool=false,
    p2m_mode::Symbol=:tiled,
    volatile_mixture_active::Bool=false,
    refractory_active::Bool=false,
    redox_active::Bool=false,
    has_metal::Bool=false,
    Xfe_bulk_val::Float64=0.2,
    has_redox::Bool=false,
    has_volatiles::Bool=false,
    XH2O_val::Float64=1.0,
    XC_val::Float64=100.0,
    XN_val::Float64=50.0,
    XS_val::Float64=100.0,
    Fm_val::Float64=0.0,
    atm_state::Union{AtmosphereState,Nothing}=nothing,
    species_tracking::Bool=false,
    dt_val::Float64=1e9,
    custom_overrides::AbstractDict=Dict{String,Any}(),
)
    cfg_dict = Dict{String,Any}(
        "grid.Nx" => 17,
        "grid.Ny" => 17,
        "geometry.gravity_mode" => gravity_mode,
        "accretion.active" => accretion_active,
        "solver.p2m_mode" => p2m_mode,
        "volatile_mixture.active" => volatile_mixture_active,
        "refractory.active" => refractory_active,
        "redox.active" => redox_active,
    )
    merge!(cfg_dict, custom_overrides)
    cfg = Erebus.override_config(default_config(), cfg_dict)
    coords = GridCoordinates(cfg.grid)
    Nx = coords.Nx
    Ny = coords.Ny
    Nx1 = coords.Nx1
    Ny1 = coords.Ny1

    # Allocate GridArrays with appropriate staggered grid dimensions
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

    # Allocate CoreGroup
    if marknum == 0
        core = Erebus.CoreGroup(
            Float64[],
            Float64[],
            Float64[],
            Int[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
        )
    else
        xm_val = coords.xcenter + r_marker
        ym_val = coords.ycenter
        core = Erebus.CoreGroup(
            [xm_val],
            [ym_val],
            [1.0],
            [tm_val],
            [tkm_val],
            [phi_val],
            [phi_val],
            [0.0],
            [0.0],
            [0.0],
            [Fm_val],
            [1e20],
            [0.0],
            [0.0],
            [1.0 / 1e10],
            [0.6],
            [1e7],
            [1e7],
            [3000.0],
            [1e6],
            [1e20],
            [0.0],
            [3.0],
            [tkm_val * 1e6],
            [1e-3],
            [1000.0],
            [3e-5],
            [2e-4],
        )
    end
    grps = Dict{Symbol,Any}()
    if marknum > 0
        if has_metal
            grps[:metal] = Erebus.MetalGroup(
                [Xfe_bulk_val], [Xfe_bulk_val], [Xfe_bulk_val], [0.0], [0.0], [0.0], [0.0]
            )
        end
        if has_redox
            grps[:redox] = Erebus.RedoxGroup(
                [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0], [0.0]
            )
        end
        if has_volatiles
            grps[:volatiles] = Erebus.VolatilesGroup(
                [XH2O_val], [XC_val], [XN_val], [XS_val], [0.0], [0.0]
            )
        end
    end
    markers = MarkerArrays(core, NamedTuple(grps))

    rplanet_init = cfg.accretion.active ? cfg.accretion.R_initial : cfg.geometry.rplanet
    M_planet_init = if cfg.accretion.active
        cfg.accretion.M_initial
    else
        (4.0 / 3.0 * pi * (rplanet_init^3) * cfg.accretion.rho_bulk)
    end
    accumulators = SimulationAccumulators(
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        100.0,
        rplanet_init,
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
        Any[],
        atm_state,
        MersenneTwister(42),
        TimerOutput(),
        1,
        dt_val,
        1e10,
    )
    return state, coords, cfg
end

@testset "Step Functions Contract and Invariance" begin
    @testset "Step 1: accrete! Contract and Invariance" begin
        # 1. Empty marker set contract: returns nothing, state bitwise unchanged
        state_empty, coords, cfg = create_mock_simulation_state(;
            marknum=0, accretion_active=true
        )
        state_empty_copy = copy(state_empty)
        res = accrete!(state_empty, coords, cfg)
        @test res === nothing
        @test isequal(
            state_empty.accumulators.rplanet, state_empty_copy.accumulators.rplanet
        )
        @test isequal(
            state_empty.accumulators.M_planet_val,
            state_empty_copy.accumulators.M_planet_val,
        )
        @test isequal(
            state_empty.accumulators.M_accreted_total,
            state_empty_copy.accumulators.M_accreted_total,
        )
        @test isequal(state_empty.markers.core.xm, state_empty_copy.markers.core.xm)
        @test isequal(state_empty.grids.ETA, state_empty_copy.grids.ETA)

        # 2. One marker state at r = R (in air, tm=3): accretion active advances boundary
        R_p = cfg.accretion.R_initial
        state_one, coords, cfg_acc = create_mock_simulation_state(;
            marknum=1, r_marker=R_p * 1.0001, tm_val=3, accretion_active=true
        )
        state_one_copy = copy(state_one)
        accrete!(state_one, coords, cfg_acc)

        # Verify mutation under # Mutates: block
        dM_dt_acc = Erebus.compute_accretion_rate(
            state_one_copy.timesum,
            state_one_copy.accumulators.M_planet_val,
            state_one_copy.accumulators.rplanet,
            cfg_acc.accretion,
            cfg_acc.disk,
        )
        dM_expected = min(
            dM_dt_acc * state_one_copy.dt,
            cfg_acc.accretion.M_target - state_one_copy.accumulators.M_planet_val,
        )
        dr_expected = Erebus.compute_radius_increment(
            state_one_copy.accumulators.rplanet, dM_expected, cfg_acc.accretion.rho_bulk
        )
        @test isapprox(
            state_one.accumulators.rplanet,
            state_one_copy.accumulators.rplanet + dr_expected;
            rtol=1e-10,
        )
        @test isapprox(
            state_one.accumulators.M_planet_val,
            state_one_copy.accumulators.M_planet_val + dM_expected;
            rtol=1e-10,
        )
        @test isapprox(state_one.accumulators.M_accreted_total, dM_expected; rtol=1e-10)
        # Verify marker phase transition from air (3) to accreted mantle (2)
        @test state_one.markers.core.tm[1] == 2
        @test isapprox(
            state_one.markers.core.phim[1], cfg_acc.accretion.phi_accreted; atol=1e-6
        )
        @test isfinite(state_one.markers.core.tkm[1])
        # Non-mutated fields remain bitwise identical
        @test isequal(state_one.grids.ETA, state_one_copy.grids.ETA)
        @test isequal(state_one.grids.RHO, state_one_copy.grids.RHO)
        @test isequal(state_one.transfers, state_one_copy.transfers)
        @test length(state_one.markers.core.xm) == 1

        # 3. Physical edge case: interior marker at r = 0 remains rock tm=1
        state_center, coords, cfg_center = create_mock_simulation_state(;
            marknum=1, r_marker=0.0, tm_val=1, accretion_active=true
        )
        accrete!(state_center, coords, cfg_center)
        @test state_center.markers.core.tm[1] == 1

        # 4. Error contract: invalid accretion parameters
        @test_throws DomainError Erebus.compute_radius_increment(-1.0, 1.0, 3000.0)
        @test_throws DomainError Erebus.compute_radius_increment(1000.0, -1.0, 3000.0)

        # 5. Branch coverage: accretion inactive returns immediately
        state_inactive, coords, cfg_inact = create_mock_simulation_state(;
            marknum=1, r_marker=R_p * 1.0001, tm_val=3, accretion_active=false
        )
        state_inact_copy = copy(state_inactive)
        accrete!(state_inactive, coords, cfg_inact)
        @test isequal(
            state_inactive.accumulators.rplanet, state_inact_copy.accumulators.rplanet
        )
        @test isequal(
            state_inactive.accumulators.M_planet_val,
            state_inact_copy.accumulators.M_planet_val,
        )

        # 6. Branch coverage: volatile condensation in accretion shell
        state_vols, coords, cfg_vols = create_mock_simulation_state(;
            marknum=1,
            r_marker=R_p * 1.0001,
            tm_val=3,
            accretion_active=true,
            volatile_mixture_active=true,
            refractory_active=true,
        )
        accrete!(state_vols, coords, cfg_vols)
        @test state_vols.accumulators.rplanet > R_p
        @test state_vols.accumulators.M_planet_val > cfg_vols.accretion.M_initial
    end

    @testset "Step 2: radiogenic_heating! Contract and Closed-Form Assertion" begin
        # 1. Empty marker set contract: returns nothing, state bitwise unchanged
        state_empty, coords, cfg = create_mock_simulation_state(; marknum=0)
        state_empty_copy = copy(state_empty)
        res = radiogenic_heating!(state_empty, coords, cfg)
        @test res === nothing
        @test isequal(
            state_empty.markers.core.hrtotalm, state_empty_copy.markers.core.hrtotalm
        )
        @test isequal(state_empty.grids.HR, state_empty_copy.grids.HR)

        # 2. One marker state: rock tm=1, phi=0.0 at r = 0
        state_one, coords, cfg = create_mock_simulation_state(;
            marknum=1, r_marker=0.0, tm_val=1, phi_val=0.0
        )
        state_one_copy = copy(state_one)
        radiogenic_heating!(state_one, coords, cfg)

        # Mutates only hrtotalm on markers
        @test state_one.markers.core.hrtotalm != state_one_copy.markers.core.hrtotalm
        @test isequal(state_one.markers.core.tkm, state_one_copy.markers.core.tkm)
        @test isequal(state_one.grids.HR, state_one_copy.grids.HR)
        @test isequal(state_one.accumulators.rplanet, state_one_copy.accumulators.rplanet)
        @test length(state_one.markers.core.xm) == 1

        # 3. Closed-form verification: specific radiogenic power at time t
        t_decay = state_one.timesum
        f_al = cfg.thermodynamics.f_al
        ratio_al = cfg.thermodynamics.ratio_al
        E_al = cfg.thermodynamics.E_al
        tau_al = cfg.thermodynamics.t_half_al / log(2.0)
        rho_rock = cfg.materials.rhosolidm[1]
        rho_metal = cfg.coreformation.rho_metal
        phi_fe_ref =
            (X_FE_REF_CHONDRITE / rho_metal) /
            (X_FE_REF_CHONDRITE / rho_metal + (1.0 - X_FE_REF_CHONDRITE) / rho_rock)

        fac_al_step = -expm1(-state_one.dt / tau_al) / (state_one.dt / tau_al)
        Q_al_analytic =
            f_al * ratio_al * E_al * exp(-t_decay / tau_al) / tau_al * fac_al_step
        Q_al_silicate = Q_al_analytic / (1.0 - X_FE_REF_CHONDRITE)
        hr_volumetric_expected = (1.0 - phi_fe_ref) * Q_al_silicate * rho_rock

        @test isapprox(
            state_one.markers.core.hrtotalm[1], hr_volumetric_expected; rtol=1e-12
        )
        # 3-class discrimination guards
        # Exponent guard: wrong lifetime differs by orders of magnitude
        wrong_decay =
            f_al * ratio_al * E_al * exp(-t_decay / (2.0 * tau_al)) / (2.0 * tau_al) /
            (1.0 - X_FE_REF_CHONDRITE) *
            (1.0 - phi_fe_ref) *
            rho_rock
        @test abs(state_one.markers.core.hrtotalm[1] - wrong_decay) >
            1e-10 * hr_volumetric_expected
        # Scale and sign guard: order of magnitude between 1e-10 and 1e2 W/m^3
        @test 1e-10 < state_one.markers.core.hrtotalm[1] < 1e2

        # 4. Physical edge case: sticky air tm=3 produces zero radiogenic heating
        state_air, coords, cfg = create_mock_simulation_state(;
            marknum=1, r_marker=0.0, tm_val=3, phi_val=0.0
        )
        radiogenic_heating!(state_air, coords, cfg)
        @test iszero(state_air.markers.core.hrtotalm[1])

        # 5. Error contract: negative decay time throws DomainError
        @test_throws DomainError Erebus.Q_radiogenic(f_al, ratio_al, E_al, tau_al, -1.0)
        @test_throws DomainError Erebus.Q_radiogenic(f_al, ratio_al, E_al, -1.0, 100.0)
        state_neg = SimulationState(
            state_one.grids,
            state_one.markers,
            state_one.accumulators,
            state_one.transfers,
            state_one.atm,
            state_one.rng,
            state_one.timer,
            state_one.timestep,
            state_one.dt,
            -1.0,
        )
        @test_throws DomainError radiogenic_heating!(state_neg, coords, cfg)

        # 5b. Half-life decay discrimination: exp(-t_half / tau) == 0.5
        state_hl = SimulationState(
            state_one.grids,
            state_one.markers,
            state_one.accumulators,
            state_one.transfers,
            state_one.atm,
            state_one.rng,
            state_one.timer,
            state_one.timestep,
            state_one.dt,
            cfg.thermodynamics.t_half_al,
        )
        radiogenic_heating!(state_hl, coords, cfg)
        Q_hl_analytic = f_al * ratio_al * E_al * 0.5 / tau_al * fac_al_step
        Q_hl_silicate = Q_hl_analytic / (1.0 - X_FE_REF_CHONDRITE)
        @test isapprox(
            state_hl.markers.core.hrtotalm[1],
            (1.0 - phi_fe_ref) * Q_hl_silicate * rho_rock;
            rtol=1e-12,
        )

        # 6. Branch coverage: metal group contribution (phi_fe > 0 and phi_fe == 0)
        state_metal, coords, cfg_metal = create_mock_simulation_state(;
            marknum=1, r_marker=0.0, tm_val=1, phi_val=0.0, has_metal=true, Xfe_bulk_val=0.2
        )
        radiogenic_heating!(state_metal, coords, cfg_metal)
        # With has_metal=true and Xfe_bulk=0.2, silicate power is Q_al / (1 - 0.2)
        Q_al_silicate_metal = Q_al_analytic / (1.0 - cfg_metal.coreformation.Xfe_bulk)
        hr_expected_metal = (1.0 - 0.2) * Q_al_silicate_metal * rho_rock
        @test isapprox(state_metal.markers.core.hrtotalm[1], hr_expected_metal; rtol=1e-12)

        state_nometal, coords, cfg_nometal = create_mock_simulation_state(;
            marknum=1, r_marker=0.0, tm_val=1, phi_val=0.0, has_metal=true, Xfe_bulk_val=0.0
        )
        radiogenic_heating!(state_nometal, coords, cfg_nometal)
        hr_expected_nometal = (1.0 - 0.0) * Q_al_silicate_metal * rho_rock
        @test isapprox(
            state_nometal.markers.core.hrtotalm[1], hr_expected_nometal; rtol=1e-12
        )

        # 7. Radionuclide budget of a body with metal: CI bulk up to the CI iron content,
        # fixed concentrations per kg silicate and per kg metal above it.
        @testset "Radionuclide budget with a metal group" begin
            function marker_heating(;
                phi_marker, phi_body, al, fe, phim=0.0, has_metal=true
            )
                st, cr, c = create_mock_simulation_state(;
                    marknum=1,
                    r_marker=0.0,
                    tm_val=1,
                    phi_val=phim,
                    has_metal=has_metal,
                    Xfe_bulk_val=phi_marker,
                    custom_overrides=Dict{String,Any}(
                        "coreformation.Xfe_bulk" => phi_body,
                        "thermodynamics.hr_al" => al,
                        "thermodynamics.hr_fe" => fe,
                    ),
                )
                radiogenic_heating!(st, cr, c)
                th = c.thermodynamics
                q_al = Erebus.Q_radiogenic(
                    th.f_al, th.ratio_al, th.E_al, th.t_half_al / log(2.0), st.timesum;
                    dt=st.dt,
                )
                q_fe = Erebus.Q_radiogenic(
                    th.f_fe, th.ratio_fe, th.E_fe, th.t_half_fe / log(2.0), st.timesum;
                    dt=st.dt,
                )
                return st.markers.core.hrtotalm[1], q_al, q_fe, c
            end
            al_only(; kw...) = marker_heating(; al=true, fe=false, kw...)
            fe_only(; kw...) = marker_heating(; al=false, fe=true, kw...)
            w_ci = 0.182
            _, _, _, c0 = al_only(; phi_marker=0.1, phi_body=0.1)
            rho_s = c0.materials.rhosolidm[1]
            rho_m = c0.coreformation.rho_metal
            rho_b(phi) = phi * rho_m + (1.0 - phi) * rho_s
            w_of(phi) = phi * rho_m / rho_b(phi)

            # At the CI metal fraction the group changes nothing.
            phi_ci = (w_ci / rho_m) / (w_ci / rho_m + (1.0 - w_ci) / rho_s)
            for (al, fe) in ((true, false), (false, true), (true, true))
                hr_metal, _, _, _ = marker_heating(;
                    phi_marker=phi_ci, phi_body=phi_ci, al=al, fe=fe
                )
                hr_plain, _, _, _ = marker_heating(;
                    phi_marker=0.0, phi_body=phi_ci, al=al, fe=fe, has_metal=false
                )
                @test isapprox(hr_metal, hr_plain; rtol=1e-12)
            end

            # Up to the CI iron content the body is CI: both nuclides follow the bulk mass.
            for phi in (0.0, 0.015, 0.10)
                @test w_of(phi) < w_ci
                hr_al, q_al, _, _ = al_only(; phi_marker=phi, phi_body=phi)
                hr_fe, _, q_fe, _ = fe_only(; phi_marker=phi, phi_body=phi)
                @test isapprox(hr_al, rho_b(phi) * q_al; rtol=1e-12)
                @test isapprox(hr_fe, rho_b(phi) * q_fe; rtol=1e-12)
            end

            # Above the CI iron content the concentrations per phase stay at the CI value.
            phi_rich = 0.20
            @test w_of(phi_rich) > w_ci
            hr_al, q_al, _, _ = al_only(; phi_marker=phi_rich, phi_body=phi_rich)
            hr_fe, _, q_fe, _ = fe_only(; phi_marker=phi_rich, phi_body=phi_rich)
            @test isapprox(hr_al, (1.0 - phi_rich) * rho_s * q_al / 0.818; rtol=1e-12)
            @test isapprox(hr_fe, phi_rich * rho_m * q_fe / 0.182; rtol=1e-12)

            # Segregated markers keep the concentrations set by the body.
            phi_body = 0.05
            w_body = w_of(phi_body)
            hr_al_metal, _, _, _ = al_only(; phi_marker=1.0, phi_body=phi_body)
            hr_fe_metal, _, q_fe, _ = fe_only(; phi_marker=1.0, phi_body=phi_body)
            @test hr_al_metal == 0.0
            @test isapprox(hr_fe_metal, rho_m * q_fe / 0.182; rtol=1e-12)
            hr_al_rock, q_al, _, _ = al_only(; phi_marker=0.0, phi_body=phi_body)
            hr_fe_rock, _, _, _ = fe_only(; phi_marker=0.0, phi_body=phi_body)
            @test isapprox(hr_al_rock, rho_s * q_al / (1.0 - w_body); rtol=1e-12)
            @test isapprox(
                hr_fe_rock, rho_s * q_fe * (0.182 - w_body) / (0.182 * (1.0 - w_body));
                rtol=1e-12,
            )

            # Moving metal between 2 equal volumes keeps the total power of the pair.
            for (al, fe) in ((true, false), (false, true))
                hr_uniform, _, _, _ = marker_heating(;
                    phi_marker=phi_body, phi_body=phi_body, al=al, fe=fe
                )
                hr_low, _, _, _ = marker_heating(;
                    phi_marker=0.0, phi_body=phi_body, al=al, fe=fe
                )
                hr_high, _, _, _ = marker_heating(;
                    phi_marker=2.0 * phi_body, phi_body=phi_body, al=al, fe=fe
                )
                @test isapprox(hr_low + hr_high, 2.0 * hr_uniform; rtol=1e-12)
            end

            # Porosity scales the whole solid.
            hr_dry, _, _, _ = fe_only(; phi_marker=0.05, phi_body=0.05)
            hr_wet, _, _, _ = fe_only(; phi_marker=0.05, phi_body=0.05, phim=0.25)
            @test isapprox(hr_wet, 0.75 * hr_dry; rtol=1e-12)
        end
    end

    @testset "Step 3: interpolate_markers_to_grid! Contract and Invariance" begin
        # 1. Empty marker set contract: returns nothing, state bitwise unchanged
        state_empty, coords, cfg = create_mock_simulation_state(; marknum=0)
        state_empty_copy = copy(state_empty)
        res = interpolate_markers_to_grid!(state_empty, coords, cfg)
        @test res === nothing
        @test isequal(state_empty.grids.RHO, state_empty_copy.grids.RHO)
        @test isequal(state_empty.grids.ETA, state_empty_copy.grids.ETA)
        @test isequal(
            state_empty.accumulators.rplanet, state_empty_copy.accumulators.rplanet
        )

        # 2. One marker state at r = 0, rock tm=1, phi=0.1
        state_one, coords, cfg = create_mock_simulation_state(;
            marknum=1, r_marker=0.0, tm_val=1, phi_val=0.1, tkm_val=500.0
        )
        radiogenic_heating!(state_one, coords, cfg)
        state_one_copy = copy(state_one)
        interpolate_markers_to_grid!(state_one, coords, cfg)

        # Mutates grid properties
        @test state_one.grids.RHO != state_one_copy.grids.RHO
        @test state_one.grids.tk1 != state_one_copy.grids.tk1
        @test state_one.grids.tk2 != state_one_copy.grids.tk2
        @test isapprox(maximum(state_one.grids.RHO), 3000.0; rtol=0.2)
        # Non-mutated state fields remain bitwise identical
        @test isequal(state_one.accumulators.rplanet, state_one_copy.accumulators.rplanet)
        @test isequal(
            state_one.accumulators.M_planet_val, state_one_copy.accumulators.M_planet_val
        )
        @test isequal(state_one.transfers, state_one_copy.transfers)
        @test length(state_one.markers.core.xm) == 1

        # 3. Physical edge cases: low porosity (phi=0.1) vs high porosity (phi=0.8)
        state_fluid, coords, cfg = create_mock_simulation_state(;
            marknum=1, r_marker=0.0, tm_val=1, phi_val=0.8, tkm_val=500.0
        )
        interpolate_markers_to_grid!(state_fluid, coords, cfg)
        @test isapprox(maximum(state_fluid.grids.PHI), 0.8; atol=1e-6)
        @test maximum(state_fluid.grids.RHO) < maximum(state_one.grids.RHO)

        # 4. Branch coverage: Pre-allocated P2MTiledWorkspace
        ws_pre = Erebus.P2MTiledWorkspace(coords, 1, cfg.solver.tile_size)
        state_ws, coords_ws, cfg_ws = create_mock_simulation_state(;
            marknum=1, r_marker=0.0, tm_val=1, phi_val=0.1, tkm_val=500.0
        )
        interpolate_markers_to_grid!(state_ws, coords_ws, cfg_ws; p2m_workspace=ws_pre)
        @test isapprox(maximum(state_ws.grids.RHO), 3000.0; rtol=0.2)

        # 5. Branch coverage: Buffered fallback without thread buffers (serial execution)
        state_ser, coords_ser, cfg_ser = create_mock_simulation_state(;
            marknum=1,
            r_marker=0.0,
            tm_val=1,
            phi_val=0.1,
            tkm_val=500.0,
            p2m_mode=:buffered,
        )
        interpolate_markers_to_grid!(state_ser, coords_ser, cfg_ser; thread_buffers=nothing)
        @test isapprox(maximum(state_ser.grids.RHO), 3000.0; rtol=0.2)

        # 6. Branch coverage: Buffered p2m_mode with thread buffers
        tb = Erebus.allocate_thread_interpolation_buffers(
            max(2, Threads.nthreads()), coords
        )
        state_buf, coords_buf, cfg_buf = create_mock_simulation_state(;
            marknum=1,
            r_marker=0.0,
            tm_val=1,
            phi_val=0.1,
            tkm_val=500.0,
            p2m_mode=:buffered,
        )
        interpolate_markers_to_grid!(state_buf, coords_buf, cfg_buf; thread_buffers=tb)
        @test isapprox(maximum(state_buf.grids.RHO), 3000.0; rtol=0.2)

        # 7. Branch coverage: Redox active with redox marker group
        state_rdx, coords_rdx, cfg_rdx = create_mock_simulation_state(;
            marknum=1,
            r_marker=0.0,
            tm_val=1,
            phi_val=0.1,
            tkm_val=500.0,
            redox_active=true,
            has_redox=true,
        )
        interpolate_markers_to_grid!(state_rdx, coords_rdx, cfg_rdx)
        @test isapprox(maximum(state_rdx.grids.RHO), 3000.0; rtol=0.2)

        # 8. Accumulator propagation: External interp_arrays mutated in-place with non-zero WTPSUM
        interp_arrs = Erebus.setup_interpolated_properties(coords)
        interp_arrs[33][1, 1] = 999.0
        interpolate_markers_to_grid!(state_one, coords, cfg; interp_arrays=interp_arrs)
        WTPSUM_res = interp_arrs[33]
        @test isapprox(maximum(WTPSUM_res), 0.25; atol=1e-6)
        @test interp_arrs[33][1, 1] != 999.0

        # 9. Branch coverage: Metal group with partition inactive ensures Xfe_S_m is nothing
        state_mnp, coords_mnp, cfg_mnp = create_mock_simulation_state(;
            marknum=1,
            r_marker=0.0,
            tm_val=1,
            phi_val=0.1,
            tkm_val=1800.0,
            has_metal=true,
            Xfe_bulk_val=0.2,
        )
        @test !cfg_mnp.metal_partition.active
        interpolate_markers_to_grid!(state_mnp, coords_mnp, cfg_mnp)
        @test isapprox(maximum(state_mnp.grids.RHO), 3000.0; rtol=0.2)
    end

    @testset "Step 4: solve_gravity! Contract and Invariance" begin
        # 1. Empty marker set contract: returns nothing, state bitwise unchanged
        state_empty, coords, cfg = create_mock_simulation_state(; marknum=0)
        state_empty_copy = copy(state_empty)
        res = solve_gravity!(state_empty, coords, cfg)
        @test res === nothing
        @test isequal(state_empty.grids.gx, state_empty_copy.grids.gx)
        @test isequal(state_empty.grids.gy, state_empty_copy.grids.gy)
        @test isequal(state_empty.grids.FI, state_empty_copy.grids.FI)

        # 2. One marker state with enclosed_mass mode
        state_one, coords, cfg_enc = create_mock_simulation_state(;
            marknum=1, r_marker=10000.0, tm_val=1, phi_val=0.1, gravity_mode=:enclosed_mass
        )
        interpolate_markers_to_grid!(state_one, coords, cfg_enc)
        state_one_copy = copy(state_one)
        solve_gravity!(state_one, coords, cfg_enc)

        # Mutates FI, gx, gy
        @test state_one.grids.gx != state_one_copy.grids.gx ||
            state_one.grids.gy != state_one_copy.grids.gy ||
            state_one.grids.FI != state_one_copy.grids.FI
        # Non-mutated fields remain bitwise identical
        @test isequal(state_one.grids.RHO, state_one_copy.grids.RHO)
        @test isequal(state_one.markers.core.xm, state_one_copy.markers.core.xm)
        @test isequal(state_one.accumulators.rplanet, state_one_copy.accumulators.rplanet)

        # 3. Symmetry and attractive direction of gravity field
        j_c = clamp(Int(floor(coords.xcenter / coords.dx)) + 1, 1, coords.Nx1)
        i_c = clamp(Int(floor(coords.ycenter / coords.dy)) + 1, 1, coords.Ny)
        @test isfinite(state_one.grids.gx[i_c, j_c])
        j_right = clamp(j_c + 2, 1, coords.Nx1)
        j_left = clamp(j_c - 2, 1, coords.Nx1)
        @test signbit(state_one.grids.gx[i_c, j_right])
        @test !signbit(state_one.grids.gx[i_c, j_left])
        @test signbit(state_one.grids.FI[i_c, j_c])

        # 4. Mode support: poisson2d mode
        state_poi, coords_poi, cfg_poi = create_mock_simulation_state(;
            marknum=1, r_marker=10000.0, tm_val=1, phi_val=0.1, gravity_mode=:poisson2d
        )
        interpolate_markers_to_grid!(state_poi, coords_poi, cfg_poi)
        solve_gravity!(state_poi, coords_poi, cfg_poi)
        @test any(!iszero, state_poi.grids.FI)
        j_poi_c = clamp(
            Int(floor(coords_poi.xcenter / coords_poi.dx)) + 1, 1, coords_poi.Nx1
        )
        i_poi_c = clamp(
            Int(floor(coords_poi.ycenter / coords_poi.dy)) + 1, 1, coords_poi.Ny
        )
        j_poi_right = clamp(j_poi_c + 2, 1, coords_poi.Nx1)
        j_poi_left = clamp(j_poi_c - 2, 1, coords_poi.Nx1)
        @test signbit(state_poi.grids.gx[i_poi_c, j_poi_right])
        @test !signbit(state_poi.grids.gx[i_poi_c, j_poi_left])
        @test signbit(state_poi.grids.FI[i_poi_c, j_poi_c])

        # 5. Error contract: invalid gravity mode rejected by config validator
        @test_throws ArgumentError Erebus.override_config(
            cfg_enc, Dict("geometry.gravity_mode" => :invalid_mode)
        )

        # 6. Branch coverage: Pre-factorized F_grav pass-through
        RP_test, _ = Erebus.setup_gravitational_lse(coords_poi)
        LP_test = Erebus.assemble_gravitational_lse!(
            zeros(coords_poi.Ny1, coords_poi.Nx1), RP_test; coords=coords_poi
        )
        F_grav_fact = lu(LP_test.cscmatrix)
        state_fgrav, coords_fgrav, cfg_fgrav = create_mock_simulation_state(;
            marknum=1, r_marker=10000.0, tm_val=1, phi_val=0.1, gravity_mode=:poisson2d
        )
        interpolate_markers_to_grid!(state_fgrav, coords_fgrav, cfg_fgrav)
        solve_gravity!(state_fgrav, coords_fgrav, cfg_fgrav; F_grav=F_grav_fact)
        @test any(!iszero, state_fgrav.grids.FI)
    end

    @testset "Step 7: vent_and_degas! Contract and Invariance" begin
        # 1. Empty marker contract: returns zero vented mass, state bitwise unchanged
        state_empty, coords, cfg = create_mock_simulation_state(; marknum=0)
        state_empty_copy = copy(state_empty)
        res = vent_and_degas!(state_empty, coords, cfg)
        @test iszero(res.delta_m_vent_3d)
        @test res.vented_vols === nothing
        @test res.degas_rates === nothing
        @test isequal(state_empty.markers.core.xm, state_empty_copy.markers.core.xm)
        @test isequal(state_empty.grids.ETA, state_empty_copy.grids.ETA)
        @test isequal(
            state_empty.accumulators.M_vent_total,
            state_empty_copy.accumulators.M_vent_total,
        )

        # 2. Single marker state with inactive venting and inactive magma degassing
        state_one, coords, cfg_one = create_mock_simulation_state(;
            marknum=1, r_marker=1000.0, tm_val=1, phi_val=0.1, tkm_val=500.0
        )
        state_one.markers.core.phinewm[1] = 0.15
        state_one.markers.core.XWsolidm[1] = 0.05
        state_one_copy = copy(state_one)
        res_one = vent_and_degas!(state_one, coords, cfg_one)
        @test isapprox(state_one.markers.core.phim[1], 0.15; rtol=1e-12)
        @test isapprox(state_one.markers.core.XWsolidm0[1], 0.05; rtol=1e-12)
        @test isequal(state_one.grids.RHO, state_one_copy.grids.RHO)
        @test isequal(state_one.accumulators.rplanet, state_one_copy.accumulators.rplanet)
        @test res_one.degas_rates === nothing

        # 3. Single marker state with active venting and volatile tracking
        overrides_v = Dict("venting.active" => true, "volatiles.active" => true)
        state_vent, coords_v, cfg_v = create_mock_simulation_state(;
            marknum=1,
            r_marker=1000.0,
            tm_val=1,
            phi_val=0.25,
            tkm_val=500.0,
            has_volatiles=true,
            custom_overrides=overrides_v,
        )
        res_v = vent_and_degas!(state_vent, coords_v, cfg_v)
        @test isfinite(res_v.delta_m_vent_3d)
        @test isfinite(state_vent.accumulators.M_vent_total)

        # 4. Coreformation active with metal group advances Xfem0
        overrides_cf = Dict("coreformation.percolation_active" => true)
        state_cf, coords_cf, cfg_cf = create_mock_simulation_state(;
            marknum=1,
            r_marker=1000.0,
            tm_val=1,
            phi_val=0.1,
            tkm_val=1400.0,
            has_metal=true,
            Xfe_bulk_val=0.3,
            custom_overrides=overrides_cf,
        )
        state_cf.markers.groups.metal.Xfem[1] = 0.25
        state_cf.markers.groups.metal.Xfem0[1] = 0.10
        vent_and_degas!(state_cf, coords_cf, cfg_cf)
        @test isapprox(state_cf.markers.groups.metal.Xfem0[1], 0.25; rtol=1e-12)

        # 4b. Metal partition active without percolation also advances Xfem0
        overrides_mp = Dict("metal_partition.active" => true, "volatiles.active" => true)
        state_mp, coords_mp, cfg_mp = create_mock_simulation_state(;
            marknum=1,
            r_marker=1000.0,
            tm_val=1,
            phi_val=0.1,
            tkm_val=1400.0,
            has_metal=true,
            has_volatiles=true,
            Xfe_bulk_val=0.3,
            custom_overrides=overrides_mp,
        )
        state_mp.markers.groups.metal.Xfem[1] = 0.28
        state_mp.markers.groups.metal.Xfem0[1] = 0.12
        vent_and_degas!(state_mp, coords_mp, cfg_mp)
        @test isapprox(state_mp.markers.groups.metal.Xfem0[1], 0.28; rtol=1e-12)

        # 5. Magma degassing with dynamic flux mode
        overrides_df = Dict(
            "magma_degassing.active" => true,
            "magma_degassing.mode" => :dynamic_flux,
            "volatiles.active" => true,
            "atmosphere.active" => true,
        )
        atm_df = AtmosphereState()
        state_df, coords_df, cfg_df = create_mock_simulation_state(;
            marknum=1,
            r_marker=1000.0,
            tm_val=1,
            phi_val=0.1,
            tkm_val=1800.0,
            Fm_val=0.5,
            has_volatiles=true,
            XH2O_val=2.0,
            XC_val=200.0,
            XN_val=50.0,
            XS_val=100.0,
            atm_state=atm_df,
            custom_overrides=overrides_df,
        )
        res_df = vent_and_degas!(state_df, coords_df, cfg_df)
        @test res_df.degas_rates isa Erebus.ElementInventory
        @test isfinite(res_df.degas_rates.H)
        @test isfinite(res_df.degas_rates.C)

        # 6. Magma degassing with equilibrium mode
        overrides_eq = Dict(
            "magma_degassing.active" => true,
            "magma_degassing.mode" => :equilibrium,
            "volatiles.active" => true,
            "atmosphere.active" => true,
        )
        atm_eq = AtmosphereState()
        atm_eq.P_surf = 1e5
        atm_eq.M_atm[:H2O] = 1e16
        state_eq, coords_eq, cfg_eq = create_mock_simulation_state(;
            marknum=1,
            r_marker=1000.0,
            tm_val=1,
            phi_val=0.1,
            tkm_val=1800.0,
            Fm_val=0.5,
            has_volatiles=true,
            XH2O_val=2.0,
            XC_val=200.0,
            XN_val=50.0,
            XS_val=100.0,
            atm_state=atm_eq,
            custom_overrides=overrides_eq,
        )
        state_eq.grids.tk2 .= 1800.0
        vent_and_degas!(state_eq, coords_eq, cfg_eq)
        @test isfinite(state_eq.markers.groups.volatiles.XH2Om[1])
        @test any(r -> r.channel === :degassing, state_eq.transfers)
    end

    @testset "Step 8: evolve_atmosphere! Contract and Invariance" begin
        # 1. Empty marker contract: returns nothing, state bitwise unchanged
        state_empty, coords, cfg = create_mock_simulation_state(; marknum=0)
        state_empty_copy = copy(state_empty)
        res = evolve_atmosphere!(state_empty, coords, cfg)
        @test res === nothing
        @test isequal(
            state_empty.accumulators.M_atm_total, state_empty_copy.accumulators.M_atm_total
        )
        @test isequal(
            state_empty.accumulators.M_escaped_total,
            state_empty_copy.accumulators.M_escaped_total,
        )
        @test isequal(state_empty.grids.ETA, state_empty_copy.grids.ETA)

        # 2. Inactive atmosphere and escape: returns nothing and state unchanged
        state_one, coords, cfg_one = create_mock_simulation_state(;
            marknum=1, r_marker=1000.0, tm_val=1
        )
        state_one_copy = copy(state_one)
        res_one = evolve_atmosphere!(state_one, coords, cfg_one)
        @test res_one === nothing
        @test isequal(
            state_one.accumulators.M_atm_total, state_one_copy.accumulators.M_atm_total
        )
        @test isequal(
            state_one.accumulators.M_escaped_total,
            state_one_copy.accumulators.M_escaped_total,
        )

        # 3. Active atmosphere with coupled escape
        atm = AtmosphereState()
        atm.M_atm[:H2] = 1.0e16
        atm.M_atm[:CO2] = 2.0e16
        overrides_atm = Dict("atmosphere.active" => true, "escape.active" => true)
        state_atm, coords_atm, cfg_atm = create_mock_simulation_state(;
            marknum=1,
            r_marker=1000.0,
            tm_val=1,
            phi_val=0.1,
            tkm_val=300.0,
            atm_state=atm,
            species_tracking=true,
            custom_overrides=overrides_atm,
        )
        evolve_atmosphere!(state_atm, coords_atm, cfg_atm)
        @test isfinite(state_atm.accumulators.M_atm_total)
        @test isapprox(
            state_atm.accumulators.M_atm_total, sum(values(state_atm.atm.M_atm)); rtol=1e-12
        )
        @test isapprox(
            state_atm.accumulators.M_escaped_total,
            sum(values(state_atm.atm.M_escaped));
            rtol=1e-12,
        )
        @test haskey(state_atm.accumulators.M_atm_species, :H2)
        @test haskey(state_atm.accumulators.M_escaped_species, :H2)

        # 3b. Active atmosphere advance with non-empty vent_degas_result coupling
        atm_coupled = AtmosphereState()
        atm_coupled.M_atm[:H2O] = 1.0e15
        state_coupled, coords_c, cfg_c = create_mock_simulation_state(;
            marknum=1,
            r_marker=1000.0,
            tm_val=1,
            phi_val=0.1,
            tkm_val=300.0,
            atm_state=atm_coupled,
            species_tracking=true,
            custom_overrides=overrides_atm,
        )
        vent_res = (;
            delta_m_vent_3d=1.0e14,
            vented_vols=Dict(:H2O => 1.0e14),
            degas_rates=Erebus.ElementInventory(1.0e10, 0.0, 0.0, 0.0, 0.0),
        )
        evolve_atmosphere!(state_coupled, coords_c, cfg_c; vent_degas_result=vent_res)
        @test isapprox(
            state_coupled.accumulators.M_atm_total,
            sum(values(state_coupled.atm.M_atm));
            rtol=1e-12,
        )
        @test isapprox(
            state_coupled.accumulators.M_escaped_total,
            sum(values(state_coupled.atm.M_escaped));
            rtol=1e-12,
        )
        @test haskey(state_coupled.accumulators.M_atm_species, :H2O)

        # 4. Active escape only with multi-species branch
        overrides_esc_ms = Dict(
            "atmosphere.active" => false,
            "escape.active" => true,
            "escape.multi_species" => true,
        )
        state_esc_ms, coords_esc_ms, cfg_esc_ms = create_mock_simulation_state(;
            marknum=1,
            r_marker=1000.0,
            tm_val=1,
            phi_val=0.1,
            tkm_val=300.0,
            species_tracking=true,
            custom_overrides=overrides_esc_ms,
        )
        state_esc_ms.accumulators.M_atm_species[:H2O] = 1.0e15
        state_esc_ms.accumulators.M_escaped_species[:H2O] = 0.0
        evolve_atmosphere!(state_esc_ms, coords_esc_ms, cfg_esc_ms)
        @test isfinite(state_esc_ms.accumulators.M_atm_total)
        @test isfinite(state_esc_ms.accumulators.M_escaped_total)
        @test haskey(state_esc_ms.accumulators.M_atm_species, :H2O)

        # 5. Active escape only with single-species branch
        overrides_esc_ss = Dict(
            "atmosphere.active" => false,
            "escape.active" => true,
            "escape.multi_species" => false,
            "escape.species" => :H2O,
        )
        state_esc_ss, coords_esc_ss, cfg_esc_ss = create_mock_simulation_state(;
            marknum=1,
            r_marker=1000.0,
            tm_val=1,
            phi_val=0.1,
            tkm_val=300.0,
            custom_overrides=overrides_esc_ss,
        )
        state_esc_ss.accumulators.M_atm_total = 1.0e15
        state_esc_ss.accumulators.M_escaped_total = 0.0
        evolve_atmosphere!(state_esc_ss, coords_esc_ss, cfg_esc_ss)
        @test isfinite(state_esc_ss.accumulators.M_atm_total)
        @test isfinite(state_esc_ss.accumulators.M_escaped_total)
    end

    @testset "Step 9: advect_markers! Contract and Closed-Form RK4 Verification" begin
        # 1. Empty marker contract: returns nothing, state bitwise unchanged
        state_empty, coords, cfg = create_mock_simulation_state(; marknum=0)
        state_empty_copy = copy(state_empty)
        res = advect_markers!(state_empty, coords, cfg)
        @test res === nothing
        @test isequal(state_empty.markers.core.xm, state_empty_copy.markers.core.xm)
        @test isequal(state_empty.grids.vx, state_empty_copy.grids.vx)

        # 2. Closed-form RK4 verification on uniform velocity field: Δx = u0*dt, Δy = v0*dt
        dt_test = 1000.0
        state_rk4, coords_rk4, cfg_rk4 = create_mock_simulation_state(;
            marknum=1, r_marker=0.0, tm_val=1, phi_val=0.1, tkm_val=300.0, dt_val=dt_test
        )
        u0 = 0.1
        v0 = -0.2
        state_rk4.grids.vx .= u0
        state_rk4.grids.vy .= v0
        state_rk4.grids.vxf .= u0
        state_rk4.grids.vyf .= v0

        x0 = state_rk4.markers.core.xm[1]
        y0 = state_rk4.markers.core.ym[1]
        advect_markers!(state_rk4, coords_rk4, cfg_rk4)

        dx_actual = state_rk4.markers.core.xm[1] - x0
        dy_actual = state_rk4.markers.core.ym[1] - y0
        dx_expected = u0 * dt_test
        dy_expected = v0 * dt_test

        @test isapprox(dx_actual, dx_expected; rtol=1e-12)
        @test isapprox(dy_actual, dy_expected; rtol=1e-12)
        j_c = clamp(Int(floor(coords_rk4.xcenter / coords_rk4.dx)) + 1, 2, coords_rk4.Nx)
        i_c = clamp(Int(floor(coords_rk4.ycenter / coords_rk4.dy)) + 1, 2, coords_rk4.Ny)
        @test isapprox(state_rk4.grids.vxp[i_c, j_c], u0; rtol=1e-12)
        @test isapprox(state_rk4.grids.vyp[i_c, j_c], v0; rtol=1e-12)
        @test isapprox(
            state_rk4.markers.core.phinewm, state_rk4.markers.core.phim; rtol=1e-12
        )

        # 3. Contract: non-mutated fields remain bitwise identical
        @test isequal(state_rk4.accumulators.rplanet, cfg_rk4.geometry.rplanet)
        @test isequal(state_rk4.grids.RHO, zeros(coords_rk4.Ny1, coords_rk4.Nx1))
    end

    @testset "Step 10: replenish! Contract and Invariance" begin
        # 1. Empty marker contract: returns 0 and leaves markers empty
        state_empty, coords, cfg = create_mock_simulation_state(; marknum=0)
        state_empty_copy = copy(state_empty)
        n_added = replenish!(state_empty, coords, cfg)
        @test n_added == 0
        @test length(state_empty.markers) == 0
        @test isequal(state_empty.markers.core.xm, state_empty_copy.markers.core.xm)

        # 2. Under-populated cell replenishment trigger
        state_rep, coords_rep, cfg_rep = create_mock_simulation_state(;
            marknum=1, r_marker=0.0, tm_val=1, phi_val=0.1, tkm_val=300.0
        )
        count_before = length(state_rep.markers)
        step_buf1 = zeros(Float64, count_before)
        step_buf2 = zeros(Float64, count_before)
        buffers = (step_buf1, step_buf2)
        count_after = replenish!(state_rep, coords_rep, cfg_rep; step_start_buffers=buffers)

        @test count_after > count_before
        @test length(state_rep.markers) == count_after
        @test length(state_rep.markers.core.w3d_m) == count_after
        @test length(step_buf1) == count_after
        @test length(step_buf2) == count_after
        @test all(isfinite, state_rep.markers.core.xm)
        @test all(isfinite, state_rep.markers.core.ym)
        @test all(isfinite, state_rep.markers.core.w3d_m)
        @test all(w -> w >= 0.0, state_rep.markers.core.w3d_m)
        @test any(w -> w > 0.0, state_rep.markers.core.w3d_m)

        # 3. Geometric weight limit at center vs edge
        rplanet = cfg_rep.geometry.rplanet
        w_center = Erebus.marker_out_of_plane_length(
            coords_rep.xcenter, coords_rep.ycenter, coords_rep.xcenter, coords_rep.ycenter
        )
        w_edge = Erebus.marker_out_of_plane_length(
            coords_rep.xcenter + rplanet,
            coords_rep.ycenter,
            coords_rep.xcenter,
            coords_rep.ycenter,
        )
        @test isapprox(w_center, 0.0; atol=1e-12)
        @test isapprox(w_edge, 2.0 * rplanet; rtol=1e-12)

        # 4. Successive pass populates remaining sparse cells
        @test replenish!(state_rep, coords_rep, cfg_rep) >= count_after
    end

    @testset "radiogenic_heating! Step-Average and Mid-Step dt Refresh" begin
        yr = Erebus.SECONDS_PER_YEAR
        t0 = 2.0e6 * yr
        dt_step = 5.0e5 * yr

        cfg = SimulationConfig(;
            grid=GridConfig(; Nx=9, Ny=9, xsize=200_000.0, ysize=200_000.0),
            geometry=GeometryConfig(;
                rplanet=90_000.0, rcrust=90_000.0, xcenter=100_000.0, ycenter=100_000.0
            ),
            time=TimeConfig(;
                n_steps=1, dt_initial=500_000.0, dt_longest=1.0e7, DTmax=20.0, yearlength=yr
            ),
            thermodynamics=ThermalConfig(; hr_al=true, hr_fe=false, phim0=0.2),
            reaction=ReactionConfig(; active=false),
            solver=SolverConfig(; p2m_mode=:tiled),
            output=OutputConfig(;
                mode=:telemetry, savematstep=100, visstep=100, save_final=false
            ),
        )
        Erebus.Config.validate_config(cfg)
        state, coords, ws = mktempdir() do tmpdir
            st, co, w, tel_io, _ = Erebus.init_simulation(cfg; output_path=tmpdir)
            tel_io !== nothing && close(tel_io)
            return st, co, w
        end
        state.timesum = t0
        state.dt = dt_step

        # Find first interior planet marker (tm < 3)
        m_idx = findfirst(t -> t < 3, state.markers.core.tm)
        @test m_idx isa Int
        t_m = state.markers.core.tm[m_idx]
        phi_m = state.markers.core.phim[m_idx]
        rho_sol = cfg.materials.rhosolidm[t_m]
        rho_met = cfg.coreformation.rho_metal
        X_fe = Erebus.X_FE_REF_CHONDRITE
        v_fe = X_fe / rho_met
        phi_fe = clamp(v_fe / (v_fe + (1.0 - X_fe) / rho_sol), 0.0, 1.0)
        tau_al = cfg.thermodynamics.t_half_al / log(2.0)
        Q_al_0 =
            cfg.thermodynamics.f_al *
            cfg.thermodynamics.ratio_al *
            cfg.thermodynamics.E_al *
            exp(-t0 / tau_al) / tau_al
        fac_al_oracle = -expm1(-dt_step / tau_al) / (dt_step / tau_al)
        hr_start_oracle = (1.0 - phi_m) * (1.0 - phi_fe) * (Q_al_0 / (1.0 - X_fe)) * rho_sol
        hr_step_avg_oracle = hr_start_oracle * fac_al_oracle

        Erebus.radiogenic_heating!(state, coords, cfg)
        @test isapprox(state.markers.core.hrtotalm[m_idx], hr_step_avg_oracle; rtol=1e-13)
        @test !isapprox(state.markers.core.hrtotalm[m_idx], hr_start_oracle; rtol=0.1)

        # Zero-step instantaneous identity (dt = 0.0)
        state.dt = 0.0
        Erebus.radiogenic_heating!(state, coords, cfg)
        @test isapprox(state.markers.core.hrtotalm[m_idx], hr_start_oracle; rtol=1e-13)

        # DomainError on negative or NaN state.dt
        state.dt = -1.0
        @test_throws DomainError Erebus.radiogenic_heating!(state, coords, cfg)
        state.dt = NaN
        @test_throws DomainError Erebus.radiogenic_heating!(state, coords, cfg)

        # 1. Mid-step dt cut refresh invariance (:tiled) with pyrolysis kinetics, pyrolysis redox, and hydration backload
        cfg_cut = SimulationConfig(;
            grid=GridConfig(; Nx=9, Ny=9, xsize=200_000.0, ysize=200_000.0),
            geometry=GeometryConfig(;
                rplanet=90_000.0, rcrust=90_000.0, xcenter=100_000.0, ycenter=100_000.0
            ),
            time=TimeConfig(;
                n_steps=1, dt_initial=500_000.0, dt_longest=1.0e7, DTmax=20.0, yearlength=yr
            ),
            materials=MaterialConfig(;
                tkm0=Erebus.SVector{3,Float64}([1100.0, 1100.0, 200.0]),
                XWsolidm_init=Erebus.SVector{3,Float64}([0.01, 0.01, NaN]),
            ),
            thermodynamics=ThermalConfig(; hr_al=true, hr_fe=true, phim0=0.2),
            reaction=ReactionConfig(;
                active=true,
                hydration_active=true,
                hydration_mode=9,
                dtreaction_hydration=1.0e11,
            ),
            refractory=RefractoryConfig(; active=true, kinetics_active=true),
            redox=RedoxConfig(; active=true, pyrolysis_redox=true),
            coreformation=CoreFormationConfig(; percolation_active=true),
            solver=SolverConfig(; p2m_mode=:tiled, max_plastic_iterations=2),
            output=OutputConfig(;
                mode=:telemetry, savematstep=100, visstep=100, save_final=false
            ),
        )
        Erebus.Config.validate_config(cfg_cut)
        coords_cut = GridCoordinates(cfg_cut.grid)

        function _init_cut_state(dt_val::Float64)
            st, w = mktempdir() do tmpdir
                s, _, wk, tel_io, _ = Erebus.init_simulation(cfg_cut; output_path=tmpdir)
                tel_io !== nothing && close(tel_io)
                return s, wk
            end
            st.timestep = 1
            st.timesum = 0.0
            st.dt = dt_val
            fill!(st.markers.core.tkm, 1100.0)
            fill!(st.grids.tk1, 1100.0)
            fill!(st.grids.tk2, 1100.0)
            for m in 1:length(st.markers)
                if st.markers.core.tm[m] < 3
                    st.markers.core.phim[m] = 0.15 + 0.05 * sin(0.1 * m)
                    st.markers.core.phinewm[m] = st.markers.core.phim[m]
                    st.markers.core.XWsolidm0[m] = 0.01
                    st.markers.core.XWsolidm[m] = 0.01
                    st.markers.groups.metal.Xfe_bulk[m] = 0.08 + 0.02 * cos(0.1 * m)
                    st.markers.groups.hcnspo.X_refr_C_m[m] = 0.05
                end
            end
            return st, w
        end

        state_a, ws_a = _init_cut_state(dt_step)
        phim_ref = copy(state_a.markers.core.phim)
        xfe_ref = copy(state_a.markers.groups.metal.Xfe_bulk)

        Erebus.prepare_step_ambient!(state_a, coords_cut, cfg_cut, ws_a)
        Erebus.radiogenic_heating!(state_a, coords_cut, cfg_cut)
        DHP_a = Erebus.update_pyrolysis!(state_a, coords_cut, cfg_cut)
        # Pyrolysis kinetics and carbothermal redox mutate both phim and Xfe_bulk before line 305
        @test !isapprox(state_a.markers.core.phim, phim_ref; rtol=0.05)
        @test !isapprox(state_a.markers.groups.metal.Xfe_bulk, xfe_ref; rtol=0.05)

        Erebus.interpolate_markers_to_grid!(
            state_a,
            coords_cut,
            cfg_cut;
            p2m_workspace=ws_a.p2m,
            thread_buffers=ws_a.thread_buffers,
            interp_arrays=ws_a.interp_arrays,
        )
        Erebus.snapshot_step_start_inventories!(ws_a, state_a, cfg_cut)
        Erebus.update_step_start_pressures!(state_a)
        Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
            Erebus.solve_thermomechanical_iterations!(
                state_a, coords_cut, cfg_cut, ws_a; DHP_pyro=DHP_a, dt_step_initial=dt_step
            )
        end

        dt_cut = state_a.dt
        @test dt_cut < dt_step

        # Fresh step-start call directly at dt_cut
        state_b, ws_b = _init_cut_state(dt_cut)
        Erebus.prepare_step_ambient!(state_b, coords_cut, cfg_cut, ws_b)
        Erebus.radiogenic_heating!(state_b, coords_cut, cfg_cut)
        Erebus.update_pyrolysis!(state_b, coords_cut, cfg_cut)
        Erebus.interpolate_markers_to_grid!(
            state_b,
            coords_cut,
            cfg_cut;
            p2m_workspace=ws_b.p2m,
            thread_buffers=ws_b.thread_buffers,
            interp_arrays=ws_b.interp_arrays,
        )

        # Closed-form oracle at dt_cut using pre-pyrolysis step-start phim and Xfe_bulk on marker m_c
        m_c = findfirst(t -> t < 3, state_a.markers.core.tm)
        t_mc = state_a.markers.core.tm[m_c]
        tau_fe = cfg_cut.thermodynamics.t_half_fe / log(2.0)
        Q_al_t0 =
            cfg_cut.thermodynamics.f_al *
            cfg_cut.thermodynamics.ratio_al *
            cfg_cut.thermodynamics.E_al / tau_al
        Q_fe_t0 =
            cfg_cut.thermodynamics.f_fe *
            cfg_cut.thermodynamics.ratio_fe *
            cfg_cut.thermodynamics.E_fe / tau_fe
        fac_al_cut = -expm1(-dt_cut / tau_al) / (dt_cut / tau_al)
        fac_fe_cut = -expm1(-dt_cut / tau_fe) / (dt_cut / tau_fe)
        X_fe_cfg = cfg_cut.coreformation.Xfe_bulk
        hr_sol_cut =
            (Q_al_t0 * fac_al_cut / (1.0 - X_fe_cfg)) * cfg_cut.materials.rhosolidm[t_mc]
        hr_met_cut = (Q_fe_t0 * fac_fe_cut / X_fe_cfg) * cfg_cut.coreformation.rho_metal
        pf_mc = clamp(xfe_ref[m_c], 0.0, 1.0)
        hr_cut_oracle =
            (1.0 - phim_ref[m_c]) * ((1.0 - pf_mc) * hr_sol_cut + pf_mc * hr_met_cut)

        @test isapprox(state_a.markers.core.hrtotalm[m_c], hr_cut_oracle; rtol=1e-13)
        @test size(state_a.grids.HR) == (coords_cut.Ny1, coords_cut.Nx1)
        @test isapprox(state_a.grids.HR, state_b.grids.HR; rtol=1e-14, atol=1e-22)
        @test !isapprox(state_a.markers.core.phim, phim_ref; rtol=0.05)
        @test !isapprox(state_a.markers.groups.metal.Xfe_bulk, xfe_ref; rtol=0.05)

        # 2. Buffered P2M branch (p2m_mode = :buffered) and Xfe_bulk === nothing branch
        cfg_buf = SimulationConfig(;
            grid=GridConfig(; Nx=9, Ny=9, xsize=200_000.0, ysize=200_000.0),
            geometry=GeometryConfig(;
                rplanet=90_000.0, rcrust=90_000.0, xcenter=100_000.0, ycenter=100_000.0
            ),
            time=TimeConfig(;
                n_steps=1, dt_initial=500_000.0, dt_longest=1.0e7, DTmax=20.0, yearlength=yr
            ),
            thermodynamics=ThermalConfig(; hr_al=true, hr_fe=false, phim0=0.2),
            reaction=ReactionConfig(; active=false),
            coreformation=CoreFormationConfig(;
                percolation_active=false, settling_active=false
            ),
            redox=RedoxConfig(; active=false),
            solver=SolverConfig(; p2m_mode=:buffered, max_plastic_iterations=2),
            output=OutputConfig(;
                mode=:telemetry, savematstep=100, visstep=100, save_final=false
            ),
        )
        Erebus.Config.validate_config(cfg_buf)
        coords_buf = GridCoordinates(cfg_buf.grid)

        function _init_buf_state(dt_val::Float64)
            st, w = mktempdir() do tmpdir
                s, _, wk, tel_io, _ = Erebus.init_simulation(cfg_buf; output_path=tmpdir)
                tel_io !== nothing && close(tel_io)
                return s, wk
            end
            st.timestep = 1
            st.timesum = 0.0
            st.dt = dt_val
            for m in 1:length(st.markers)
                if st.markers.core.tm[m] < 3
                    st.markers.core.phim[m] = 0.15 + 0.05 * sin(0.1 * m)
                    st.markers.core.phinewm[m] = st.markers.core.phim[m]
                end
            end
            return st, w
        end

        state_buf_a, ws_buf_a = _init_buf_state(dt_step)
        @test !haskey(state_buf_a.markers.groups, :metal)
        Erebus.prepare_step_ambient!(state_buf_a, coords_buf, cfg_buf, ws_buf_a)
        Erebus.radiogenic_heating!(state_buf_a, coords_buf, cfg_buf)
        Erebus.interpolate_markers_to_grid!(
            state_buf_a,
            coords_buf,
            cfg_buf;
            p2m_workspace=ws_buf_a.p2m,
            thread_buffers=ws_buf_a.thread_buffers,
            interp_arrays=ws_buf_a.interp_arrays,
        )
        Erebus.snapshot_step_start_inventories!(ws_buf_a, state_buf_a, cfg_buf)
        Erebus.update_step_start_pressures!(state_buf_a)
        Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
            Erebus.solve_thermomechanical_iterations!(
                state_buf_a, coords_buf, cfg_buf, ws_buf_a; dt_step_initial=dt_step
            )
        end
        dt_buf_cut = state_buf_a.dt
        @test dt_buf_cut < dt_step

        state_buf_b, ws_buf_b = _init_buf_state(dt_buf_cut)
        Erebus.prepare_step_ambient!(state_buf_b, coords_buf, cfg_buf, ws_buf_b)
        Erebus.radiogenic_heating!(state_buf_b, coords_buf, cfg_buf)
        Erebus.interpolate_markers_to_grid!(
            state_buf_b,
            coords_buf,
            cfg_buf;
            p2m_workspace=ws_buf_b.p2m,
            thread_buffers=ws_buf_b.thread_buffers,
            interp_arrays=ws_buf_b.interp_arrays,
        )
        @test isapprox(state_buf_a.grids.HR, state_buf_b.grids.HR; rtol=1e-14, atol=1e-22)

        # 3. Exception safety (try/finally restoration) and guard contracts
        @test_throws DimensionMismatch Erebus._swap_step_start_vec!([1.0, 2.0], [1.0])
        phim_before_err = copy(state_a.markers.core.phim)
        xfe_before_err = copy(state_a.markers.groups.metal.Xfe_bulk)
        state_a.dt = -1.0
        @test_throws DomainError Erebus._update_radiogenic_heating_for_dt!(
            state_a, coords_cut, cfg_cut, ws_a
        )
        state_a.dt = dt_cut
        @test isapprox(state_a.markers.core.phim, phim_before_err; rtol=1e-15)
        @test isapprox(state_a.markers.groups.metal.Xfe_bulk, xfe_before_err; rtol=1e-15)

        # Missing step-start snapshot for the current timestep throws ArgumentError
        state_uninit, ws_uninit = _init_cut_state(dt_step)
        state_uninit.timestep = 7
        Erebus.radiogenic_heating!(state_uninit, coords_cut, cfg_cut)
        Erebus.interpolate_markers_to_grid!(
            state_uninit,
            coords_cut,
            cfg_cut;
            p2m_workspace=ws_uninit.p2m,
            thread_buffers=ws_uninit.thread_buffers,
            interp_arrays=ws_uninit.interp_arrays,
        )
        Erebus.snapshot_step_start_inventories!(ws_uninit, state_uninit, cfg_cut)
        Erebus.update_step_start_pressures!(state_uninit)
        err_msg_snap = "Radiogenic step-start state was not snapshotted for timestep 7 (hr_step_start_timestep=-1)"
        @test_throws ArgumentError Erebus._update_radiogenic_heating_for_dt!(
            state_uninit, coords_cut, cfg_cut, ws_uninit
        )
        @test_throws err_msg_snap Erebus._update_radiogenic_heating_for_dt!(
            state_uninit, coords_cut, cfg_cut, ws_uninit
        )
        @test_throws ArgumentError Base.CoreLogging.with_logger(
            Base.CoreLogging.NullLogger()
        ) do
            Erebus.solve_thermomechanical_iterations!(
                state_uninit, coords_cut, cfg_cut, ws_uninit; dt_step_initial=dt_step
            )
        end
        @test_throws err_msg_snap Base.CoreLogging.with_logger(
            Base.CoreLogging.NullLogger()
        ) do
            Erebus.solve_thermomechanical_iterations!(
                state_uninit, coords_cut, cfg_cut, ws_uninit; dt_step_initial=dt_step
            )
        end

        # Unbinned tiled P2M workspace throws DimensionMismatch before @inbounds loop
        saved_end_offset = ws_a.p2m.tile_offsets[end]
        try
            ws_a.p2m.tile_offsets[end] = 0
            @test_throws DimensionMismatch Erebus._update_radiogenic_heating_for_dt!(
                state_a, coords_cut, cfg_cut, ws_a
            )
        finally
            ws_a.p2m.tile_offsets[end] = saved_end_offset
        end

        # 4. Accretion-on with has_metal=true: direct snapshot contract and two-attempt DTmax cut
        dt_acc_init_yr = 40_000.0
        T0_acc = 250.0
        sv3_acc(a, b, c) = Erebus.SVector{3,Float64}([a, b, c])
        cfg_acc = SimulationConfig(;
            grid=GridConfig(; Nx=9, Ny=9, xsize=200_000.0, ysize=200_000.0),
            geometry=GeometryConfig(;
                rplanet=50_000.0, rcrust=50_000.0, xcenter=100_000.0, ycenter=100_000.0
            ),
            time=TimeConfig(;
                n_steps=1,
                start_time=0.0,
                dt_initial=dt_acc_init_yr,
                dt_longest=1.0e7,
                DTmax=0.5,
                dxymax=1.0,
                yearlength=yr,
            ),
            materials=MaterialConfig(;
                alphasolidm=sv3_acc(0.0, 0.0, 0.0),
                alphafluidm=sv3_acc(0.0, 0.0, 0.0),
                kphim0=sv3_acc(1.0e-30, 1.0e-30, 1.0e-30),
                tkm0=sv3_acc(T0_acc, T0_acc, T0_acc),
            ),
            thermodynamics=ThermalConfig(;
                hr_al=true,
                hr_fe=true,
                ratio_al=5.0e-3,
                t_half_al=5000.0 * yr,
                phim0=0.20,
                thermal_buoyancy=false,
                surface_radiation=false,
            ),
            disk=DiskConfig(; enabled=false, model=:fixed, t_ambient=T0_acc),
            accretion=AccretionConfig(;
                active=true,
                mode=:linear_radius,
                R_initial=50_000.0,
                R_target=85_000.0,
                M_target=1.0e25,
                t_start_myr=0.0,
                t_duration_myr=10.0,
                dR_dt_constant=1.0,
                h_impact=0.0,
                phi_accreted=0.35,
                Xfe_bulk_accreted=0.25,
            ),
            reaction=ReactionConfig(; active=false),
            melting=MeltingConfig(; active=false),
            coreformation=CoreFormationConfig(;
                percolation_active=false, settling_active=false, Xfe_bulk=0.15
            ),
            redox=RedoxConfig(; active=true),
            venting=VentingConfig(; active=false),
            solver=SolverConfig(; p2m_mode=:tiled, max_plastic_iterations=2),
            output=OutputConfig(;
                mode=:snapshots, savematstep=1, visstep=100, save_final=false
            ),
        )
        Erebus.Config.validate_config(cfg_acc)

        # Direct step-sequence check: prepare_step_ambient! -> accrete! -> radiogenic_heating! -> snapshot_radiogenic_step_start!
        st_dir, coords_dir, ws_dir = mktempdir() do tmp_d
            s, c, w, tel_io, _ = Erebus.init_simulation(cfg_acc; output_path=tmp_d)
            tel_io !== nothing && close(tel_io)
            return s, c, w
        end
        @test ws_dir.Xfe_bulk_hr_step_start isa Vector{Float64}
        st_dir.timestep = 1
        st_dir.timesum = 0.0
        st_dir.dt = dt_acc_init_yr * yr
        tm_pre_dir = copy(st_dir.markers.core.tm)
        Erebus.prepare_step_ambient!(st_dir, coords_dir, cfg_acc, ws_dir)
        Erebus.accrete!(st_dir, coords_dir, cfg_acc)
        acc_dir_idx = findall(
            m -> tm_pre_dir[m] == 3 && st_dir.markers.core.tm[m] == 2, eachindex(tm_pre_dir)
        )
        @test !isempty(acc_dir_idx)
        @test all(
            m -> ws_dir.phim_hr_step_start[m] == cfg_acc.poroelasticity.phimin, acc_dir_idx
        )
        @test all(m -> iszero(ws_dir.Xfe_bulk_hr_step_start[m]), acc_dir_idx)
        Erebus.radiogenic_heating!(st_dir, coords_dir, cfg_acc)
        Erebus.snapshot_radiogenic_step_start!(ws_dir, st_dir, cfg_acc)
        @test ws_dir.phim_hr_step_start == st_dir.markers.core.phim
        @test ws_dir.Xfe_bulk_hr_step_start == st_dir.markers.groups.metal.Xfe_bulk
        @test all(
            m -> ws_dir.phim_hr_step_start[m] == cfg_acc.accretion.phi_accreted, acc_dir_idx
        )
        @test all(
            m -> ws_dir.Xfe_bulk_hr_step_start[m] == cfg_acc.accretion.Xfe_bulk_accreted,
            acc_dir_idx,
        )
        # Measure Attempt 1's DTmax-reduced dt to verify Attempt 2 cuts dt again
        Erebus.interpolate_markers_to_grid!(
            st_dir,
            coords_dir,
            cfg_acc;
            p2m_workspace=ws_dir.p2m,
            thread_buffers=ws_dir.thread_buffers,
            interp_arrays=ws_dir.interp_arrays,
        )
        Erebus.solve_gravity!(
            st_dir, coords_dir, cfg_acc; F_grav=ws_dir.F_grav, RP=ws_dir.RP, SP=ws_dir.SP
        )
        Erebus.apply_surface_radiation!(st_dir, coords_dir, cfg_acc)
        Erebus.snapshot_step_start_inventories!(ws_dir, st_dir, cfg_acc)
        Erebus.update_step_start_pressures!(st_dir)
        Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
            Erebus.solve_thermomechanical_iterations!(
                st_dir, coords_dir, cfg_acc, ws_dir; dt_step_initial=dt_acc_init_yr * yr
            )
        end
        dt_attempt1 = st_dir.dt
        @test dt_attempt1 < dt_acc_init_yr * yr

        mktempdir() do tmpdir
            out_dir = tmpdir * "/"
            Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
                Erebus.simulation_loop(cfg_acc; output_path=out_dir)
            end
            d0 = JLD2.load(joinpath(tmpdir, "output_00000.jld2"))
            d1 = JLD2.load(joinpath(tmpdir, "output_00001.jld2"))
            dt_acc_final = d1["dt"]
            # Attempt 2 must cut dt below Attempt 1's dt so _update_radiogenic_heating_for_dt! runs on Attempt 2
            @test dt_acc_final < dt_attempt1
            accreted_idx = findall(
                m -> d0["tm"][m] == 3 && d1["tm"][m] == 2, eachindex(d0["tm"])
            )
            @test !isempty(accreted_idx)

            function _eval_fresh_post_acc(dt_val::Float64)
                st_tmp, coords_tmp, ws_tmp = mktempdir() do tmp_f
                    s, c, w, tel_io, _ = Erebus.init_simulation(cfg_acc; output_path=tmp_f)
                    tel_io !== nothing && close(tel_io)
                    return s, c, w
                end
                st_tmp.timestep = 1
                st_tmp.timesum = 0.0
                st_tmp.dt = dt_val
                for m in accreted_idx
                    st_tmp.markers.core.tm[m] = Int8(2)
                    st_tmp.markers.core.phim[m] = cfg_acc.accretion.phi_accreted
                    st_tmp.markers.core.phinewm[m] = cfg_acc.accretion.phi_accreted
                    st_tmp.markers.groups.metal.Xfe_bulk[m] =
                        cfg_acc.accretion.Xfe_bulk_accreted
                end
                Erebus.radiogenic_heating!(st_tmp, coords_tmp, cfg_acc)
                Erebus.interpolate_markers_to_grid!(
                    st_tmp,
                    coords_tmp,
                    cfg_acc;
                    p2m_workspace=ws_tmp.p2m,
                    thread_buffers=ws_tmp.thread_buffers,
                    interp_arrays=ws_tmp.interp_arrays,
                )
                return st_tmp
            end

            st_f = _eval_fresh_post_acc(dt_acc_final)
            st_pre_cut = _eval_fresh_post_acc(dt_attempt1)

            tau_al_acc = cfg_acc.thermodynamics.t_half_al / log(2.0)
            tau_fe_acc = cfg_acc.thermodynamics.t_half_fe / log(2.0)
            Q_al_0 =
                cfg_acc.thermodynamics.f_al *
                cfg_acc.thermodynamics.ratio_al *
                cfg_acc.thermodynamics.E_al / tau_al_acc
            Q_fe_0 =
                cfg_acc.thermodynamics.f_fe *
                cfg_acc.thermodynamics.ratio_fe *
                cfg_acc.thermodynamics.E_fe / tau_fe_acc
            X_fe_ref_acc = cfg_acc.coreformation.Xfe_bulk
            fac_al_final = -expm1(-dt_acc_final / tau_al_acc) / (dt_acc_final / tau_al_acc)
            fac_fe_final = -expm1(-dt_acc_final / tau_fe_acc) / (dt_acc_final / tau_fe_acc)
            hrs_crust =
                (Q_al_0 * fac_al_final / (1.0 - X_fe_ref_acc)) *
                cfg_acc.materials.rhosolidm[2]
            hrm_crust =
                (Q_fe_0 * fac_fe_final / X_fe_ref_acc) * cfg_acc.coreformation.rho_metal
            phi_acc = cfg_acc.accretion.phi_accreted
            xfe_acc = cfg_acc.accretion.Xfe_bulk_accreted
            hr_acc_oracle =
                (1.0 - phi_acc) * ((1.0 - xfe_acc) * hrs_crust + xfe_acc * hrm_crust)
            phi_pre_acc_air = cfg_acc.poroelasticity.phimin
            hr_acc_pre_accretion_air = (1.0 - phi_pre_acc_air) * hrs_crust
            hr_acc_stale_xfe_zero = (1.0 - phi_acc) * hrs_crust

            fac_al_pre_cut = -expm1(-dt_attempt1 / tau_al_acc) / (dt_attempt1 / tau_al_acc)
            fac_fe_pre_cut = -expm1(-dt_attempt1 / tau_fe_acc) / (dt_attempt1 / tau_fe_acc)
            hrs_crust_pre_cut =
                (Q_al_0 * fac_al_pre_cut / (1.0 - X_fe_ref_acc)) *
                cfg_acc.materials.rhosolidm[2]
            hrm_crust_pre_cut =
                (Q_fe_0 * fac_fe_pre_cut / X_fe_ref_acc) * cfg_acc.coreformation.rho_metal
            hr_acc_attempt2_pre_cut =
                (1.0 - phi_acc) *
                ((1.0 - xfe_acc) * hrs_crust_pre_cut + xfe_acc * hrm_crust_pre_cut)

            @test isapprox(d1["hrtotalm"][accreted_idx[1]], hr_acc_oracle; rtol=1e-14)
            @test !isapprox(
                d1["hrtotalm"][accreted_idx[1]], hr_acc_attempt2_pre_cut; rtol=1e-7
            )
            @test !isapprox(
                d1["hrtotalm"][accreted_idx[1]], hr_acc_pre_accretion_air; rtol=0.3
            )
            @test !isapprox(
                d1["hrtotalm"][accreted_idx[1]], hr_acc_stale_xfe_zero; rtol=0.1
            )
            @test isapprox(d1["hrtotalm"], st_f.markers.core.hrtotalm; rtol=1e-14)
            @test isapprox(d1["HR"], st_f.grids.HR; rtol=1e-14, atol=1e-22)
            @test !isapprox(d1["HR"], st_pre_cut.grids.HR; rtol=1e-7)
        end
    end
end

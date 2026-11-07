using Test
using Random
using LinearAlgebra
using Erebus

@testset "Simulation Modular Steps and Snapshots" begin
    @testset "StepSnapshot Round-Trip and In-Place Restoration" begin
        cfg = default_config()
        coords = GridCoordinates(cfg.grid)
        state, coords_actual, ws, telemetry_io, start_step = init_simulation(cfg)

        state.grids.pr .= 1.234e6
        state.grids.tk1 .= 450.0
        state.accumulators.rplanet = 1.05e5
        state.accumulators.M_planet_val = 5.67e19
        ws.YERRNOD[1] = 9.87e3

        snapshot = snapshot_step_state(state, coords_actual, ws)

        state.grids.pr .= 9.999e6
        state.grids.tk1 .= 800.0
        state.accumulators.rplanet = 2.0e5
        state.accumulators.M_planet_val = 1.0e20
        ws.YERRNOD[1] = 0.0

        coords_ref = Ref(coords_actual)
        restore_step_state!(state, coords_ref, ws, snapshot)

        @test isapprox(state.grids.pr[1, 1], 1.234e6; rtol=1e-12)
        @test isapprox(state.grids.tk1[1, 1], 450.0; rtol=1e-12)
        @test isapprox(state.accumulators.rplanet, 1.05e5; rtol=1e-12)
        @test isapprox(state.accumulators.M_planet_val, 5.67e19; rtol=1e-12)
        @test isapprox(ws.YERRNOD[1], 9.87e3; rtol=1e-12)

        state.grids.pr .= 3.456e6
        state.accumulators.rplanet = 1.15e5
        snapshot_step_state!(snapshot, state, coords_actual, ws)
        @test isapprox(snapshot.grids.pr[1, 1], 3.456e6; rtol=1e-12)
        @test isapprox(snapshot.scalars.rplanet_val, 1.15e5; rtol=1e-12)

        alloc_snap, alloc_rest =
            let snap_loc = snapshot,
                st_loc = state,
                co_loc = coords_actual,
                cref_loc = coords_ref,
                ws_loc = ws

                snapshot_step_state!(snap_loc, st_loc, co_loc, ws_loc)
                restore_step_state!(st_loc, cref_loc, ws_loc, snap_loc)
                a_s = @allocated snapshot_step_state!(snap_loc, st_loc, co_loc, ws_loc)
                a_r = @allocated restore_step_state!(st_loc, cref_loc, ws_loc, snap_loc)
                (a_s, a_r)
            end
        @test alloc_snap == 0
        @test alloc_rest == 0

        # Grid resize fallback branch in snapshot_step_state! (e.g. after telescoping)
        cfg_big = SimulationConfig(grid=GridConfig(Nx=21, Ny=21))
        state_big, coords_big, ws_big, _, _ = init_simulation(cfg_big)
        state_big.grids.pr .= 7.89e6
        ws_big.YERRNOD[end] = 42.0
        snapshot_step_state!(snapshot, state_big, coords_big, ws_big)
        @test size(snapshot.grids.pr) == size(state_big.grids.pr)
        @test isapprox(snapshot.grids.pr[1, 1], 7.89e6; rtol=1e-12)
        @test length(snapshot.YERRNOD) == length(ws_big.YERRNOD)
        @test isapprox(snapshot.YERRNOD[end], 42.0; rtol=1e-12)

        # Single-marker restore with optional groups
        cfg_grp = load_config("configs/core_formation_benchmark.toml")
        state_grp, _, _, _, _ = init_simulation(cfg_grp)
        snap_grp = copy(state_grp.markers)
        state_grp.markers.core.tkm[1] = 9999.0
        state_grp.markers.groups.metal.Xfe_bulk[1] = 0.99
        Erebus._restore_single_marker!(state_grp.markers, snap_grp, 1)
        @test isapprox(state_grp.markers.core.tkm[1], snap_grp.core.tkm[1]; rtol=1e-12)
        @test isapprox(
            state_grp.markers.groups.metal.Xfe_bulk[1],
            snap_grp.groups.metal.Xfe_bulk[1];
            rtol=1e-12,
        )
    end

    @testset "In-Place copy_grid_arrays! Accuracy" begin
        cfg = default_config()
        coords = GridCoordinates(cfg.grid)
        state1, coords1, ws1, _, _ = init_simulation(cfg)
        state2, coords2, ws2, _, _ = init_simulation(cfg)

        state1.grids.ETA .= 7.77e20
        state1.grids.tk2 .= 555.5
        state1.grids.vx .= 0.123

        copy_grid_arrays!(state2.grids, state1.grids)

        @test isapprox(state2.grids.ETA[1, 1], 7.77e20; rtol=1e-12)
        @test isapprox(state2.grids.tk2[1, 1], 555.5; rtol=1e-12)
        @test isapprox(state2.grids.vx[1, 1], 0.123; rtol=1e-12)
    end

    @testset "Simulation Step Ambient and Surface Radiation" begin
        cfg = default_config()
        coords = GridCoordinates(cfg.grid)
        state, coords_actual, ws, _, _ = init_simulation(cfg)

        _, P_amb_exp, _ = compute_ambient_conditions(state.timesum, cfg.disk)
        prepare_step_ambient!(state, coords_actual, cfg, ws)
        @test isapprox(state.accumulators.P_amb, P_amb_exp; rtol=1e-12)
        @test isapprox(maximum(ws.interp_arrays[1]), 0.0; atol=1e-12)

        apply_surface_radiation!(state, coords_actual, cfg)
        @test all(isfinite, state.grids.KX)
        @test all(isfinite, state.grids.KY)
    end

    @testset "Simulation Step Start Inventories and Pressure Snapshots" begin
        cfg = default_config()
        coords = GridCoordinates(cfg.grid)
        state, coords_actual, ws, _, _ = init_simulation(cfg)

        state.grids.pr .= 2.5e6
        state.grids.pf .= 1.8e6
        state.grids.ps .= 7.0e5
        update_step_start_pressures!(state)

        @test isapprox(state.grids.pr0[1, 1], 2.5e6; rtol=1e-12)
        @test isapprox(state.grids.pf0[1, 1], 1.8e6; rtol=1e-12)
        @test isapprox(state.grids.ps0[1, 1], 7.0e5; rtol=1e-12)

        snapshot_step_start_inventories!(ws, state, cfg)
        @test isa(ws.hydromech, HydromechanicalLSEWorkspace)
        @test isa(ws.thermal, ThermalLSEWorkspace)
    end

    @testset "Workspace Reset on Domain Expansion" begin
        cfg = default_config()
        coords = GridCoordinates(cfg.grid)
        state, coords_actual, ws, _, _ = init_simulation(cfg)

        ws.hydromech_cache = "dummy_cache"
        ws.thermal_cache = "dummy_cache"

        coords_new = compute_telescoped_coordinates(coords_actual)
        reset_workspaces_for_grid!(ws, coords_new, cfg, length(state.markers))

        @test ws.hydromech_cache === nothing
        @test ws.thermal_cache === nothing
        @test size(ws.fractured_cells) == (coords_new.Ny, coords_new.Nx)
        dof_per_node = cfg.solver.darcy_elimination ? 4 : 6
        @test length(ws.R) == coords_new.Ny1 * coords_new.Nx1 * dof_per_node
    end

    @testset "Diffuse and Update Markers Post-Convergence" begin
        cfg = default_config()
        coords = GridCoordinates(cfg.grid)
        state, coords_actual, ws, _, _ = init_simulation(cfg)

        state.grids.DT .= 5.0
        diffuse_and_update_markers!(state, coords_actual, cfg, ws)

        @test isapprox(state.grids.DT0[1, 1], 5.0; rtol=1e-12)
        @test all(isfinite, state.markers.core.etavpm)
        @test all(>(0.0), state.markers.core.etavpm)
    end

    @testset "Marker Temperature Range in the Post-Convergence Update" begin
        # Equal marker weights, so the heat budget below is a plain sum.
        cfg = override_config(
            default_config(), Dict{String,Any}("geometry.spherical_metric" => false)
        )
        state, coords_actual, ws, _, _ = init_simulation(cfg)
        core = state.markers.core
        rock = findfirst(<(3), core.tm)
        core.tkm[rock] = 20.0
        state.grids.tk2 .= 300.0
        state.grids.DT .= -5.0

        # The first step sets markers from the grid, so it takes no correction.
        state.timestep = 1
        tkm_first = copy(core.tkm)
        diffuse_and_update_markers!(state, coords_actual, cfg, ws)
        @test core.tkm ≈ tkm_first rtol = 1.0e-14

        # The projection of a real step sets the heat capacity of each rock marker.
        rocks = findall(<(3), core.tm)
        core.rhocptotalm[rocks] .= 3.0e6
        core.tkm[rocks] .= 500.0
        core.tkm[rock] = 20.0
        tkm_before = copy(core.tkm)
        state.timestep = 2
        state.grids.DT .= -5.0
        diffuse_and_update_markers!(state, coords_actual, cfg, ws)
        @test core.tkm[rock] ≈ 25.0 atol = 1.0e-9
        # The other rock markers take the cooling that this marker rejects.
        @test sum(core.tkm[rocks] .- tkm_before[rocks]) ≈ 0.0 atol = 1.0e-9
        @test sum(min.(core.tkm[rocks] .- tkm_before[rocks], 0.0)) ≈ -5.0 atol = 1.0e-9
    end

    @testset "Marker Temperature Range with Subgrid Relaxation" begin
        cfg = override_config(default_config(), Dict{String,Any}("solver.dsubgridt" => 1.0))
        # With no conduction in the increment the limit must not act: the first step skips
        # the limit, and the later step must give the same marker temperatures.
        function relaxed_tkm(timestep)
            state, coords_actual, ws, _, _ = init_simulation(cfg)
            core = state.markers.core
            rocks = findall(<(3), core.tm)
            core.rhocptotalm[rocks] .= 3.0e6
            core.tkm[rocks] .= 100.0
            state.grids.tk1 .= 300.0
            state.grids.tk2 .= 300.0
            state.grids.DT .= 0.0
            state.timestep = timestep
            diffuse_and_update_markers!(state, coords_actual, cfg, ws)
            return copy(core.tkm)
        end
        tkm_free = relaxed_tkm(1)
        @test relaxed_tkm(2) ≈ tkm_free rtol = 1.0e-12
        # The relaxation moved the markers, so the comparison holds a real change.
        @test !isapprox(tkm_free[tkm_free .< 200.0], fill(100.0, count(tkm_free .< 200.0)))
    end

    @testset "Marker Temperature Range under the Spherical Metric" begin
        cfg = override_config(
            default_config(), Dict{String,Any}("geometry.spherical_metric" => true)
        )
        state, coords_actual, ws, _, _ = init_simulation(cfg)
        core = state.markers.core
        rocks = findall(<(3), core.tm)
        core.rhocptotalm[rocks] .= 3.0e6
        core.tkm[rocks] .= 500.0
        rock = rocks[argmin(abs.(core.w3d_m[rocks] .- 0.5 * maximum(core.w3d_m[rocks])))]
        core.tkm[rock] = 20.0
        tkm_before = copy(core.tkm)
        state.grids.tk2 .= 300.0
        state.grids.DT .= -5.0
        state.timestep = 2
        diffuse_and_update_markers!(state, coords_actual, cfg, ws)
        change = core.tkm[rocks] .- tkm_before[rocks]
        @test core.tkm[rock] ≈ 25.0 atol = 1.0e-9
        # The heat closes only when each marker counts with its out-of-plane length.
        @test sum(core.w3d_m[rocks] .* change) ≈ 0.0 atol = 1.0e-9 * 5.0 * core.w3d_m[rock]
    end

    @testset "Enclosed Mass Gravity around an Off-Centre Planet" begin
        cfg = override_config(
            default_config(), Dict{String,Any}("geometry.xcenter" => 60_000.0)
        )
        state, coords_actual, ws, _, _ = init_simulation(cfg)
        @test cfg.geometry.gravity_mode === :enclosed_mass
        core = state.markers.core
        # The projection of a real step sets the density of each rock marker.
        core.rhototalm[core.tm .< 3] .= 3000.0
        solve_gravity!(state, coords_actual, cfg)
        # The result must match a grid whose own centre is the planet centre.
        coords_planet = GridCoordinates(cfg.grid; xcenter=60_000.0, ycenter=70_000.0)
        gx = zeros(size(state.grids.gx))
        gy = zeros(size(state.grids.gy))
        Erebus.compute_gravity_enclosed_mass!(
            gx,
            gy;
            xm=core.xm,
            ym=core.ym,
            rhototalm=core.rhototalm,
            tm=core.tm,
            coords=coords_planet,
            gravity_nr_factor=cfg.geometry.gravity_nr_factor,
            rplanet=state.accumulators.rplanet,
        )
        @test state.grids.gx ≈ gx rtol = 1.0e-12
        @test state.grids.gy ≈ gy rtol = 1.0e-12
        # The grid centre alone gives another field, so the comparison can fail.
        gx_grid = zeros(size(gx))
        Erebus.compute_gravity_enclosed_mass!(
            gx_grid,
            zeros(size(gy));
            xm=core.xm,
            ym=core.ym,
            rhototalm=core.rhototalm,
            tm=core.tm,
            coords=coords_actual,
            gravity_nr_factor=cfg.geometry.gravity_nr_factor,
            rplanet=state.accumulators.rplanet,
        )
        @test !isapprox(gx_grid, gx; rtol=1.0e-3)
    end

    @testset "Telescoping Trigger Evaluation" begin
        cfg = default_config()
        coords = GridCoordinates(cfg.grid)
        state, coords_actual, ws, _, _ = init_simulation(cfg)
        coords_ref = Ref(coords_actual)

        res = telescope_domain!(state, coords_ref, cfg, ws)
        @test res == false
        @test coords_ref[].Nx == coords_actual.Nx
    end

    @testset "SimulationState Delegation, Indexing, and Fallbacks" begin
        cfg = default_config()
        coords = GridCoordinates(cfg.grid)
        state, coords_actual, ws, _, _ = init_simulation(cfg)

        @test isapprox(state.rplanet, state.accumulators.rplanet; rtol=1e-12)
        @test isapprox(state[:rplanet], state.accumulators.rplanet; rtol=1e-12)
        @test isapprox(state["rplanet"], state.accumulators.rplanet; rtol=1e-12)
        @test state[:ETA] === state.grids.ETA
        @test state["ETA"] === state.grids.ETA
        @test state[:transfer_log] === state.transfers
        @test state[:S_vent] === state.grids.S_vent_grid

        @test haskey(state, :rplanet)
        @test haskey(state, "rplanet")
        @test haskey(state, :ETA)
        @test haskey(state, :transfer_log)
        @test haskey(state, :S_vent)
        @test !haskey(state, :nonexistent_field)

        state.rplanet = 1.234e5
        @test isapprox(state.accumulators.rplanet, 1.234e5; rtol=1e-12)

        @test_throws KeyError state[:nonexistent_field]
        @test_throws ErrorException state.nonexistent_field
        @test_throws ErrorException (state.nonexistent_field = 1.0)

        @test :rplanet in propertynames(state)
        @test :ETA in propertynames(state.grids)

        state_copied = copy(state)
        @test isapprox(state_copied.dt, state.dt; rtol=1e-12)
        @test any(p -> p.first === :rplanet, pairs(state))
        @test any(p -> p.first === :ETA, pairs(state.grids))

        nt_state = NamedTuple(state)
        @test haskey(nt_state, :grids)
        @test haskey(nt_state, :markers)

        nt_snap = snapshot_step_state((; a=[1.0, 2.0], b=3.0))
        @test isapprox(nt_snap.a[1], 1.0; rtol=1e-12)
        @test isapprox(nt_snap.b, 3.0; rtol=1e-12)

        dict_target = Dict(:a => [0.0, 0.0])
        restore_step_state!(dict_target, Dict(:a => [4.0, 5.0]))
        @test isapprox(dict_target[:a][1], 4.0; rtol=1e-12)
        @test isapprox(dict_target[:a][2], 5.0; rtol=1e-12)

        err_buf = IOBuffer()
        showerror(err_buf, Erebus.CheckpointError("test checkpoint error"))
        @test occursin("CheckpointError", String(take!(err_buf)))
    end

    @testset "Simulation Step Ambient and Surface Radiation with Atmosphere" begin
        cfg = default_config()
        cfg_atm = SimulationConfig(;
            (
                f => (f === :atmosphere ? AtmosphereConfig(active=true) : getfield(cfg, f))
                for f in fieldnames(SimulationConfig)
            )...,
        )
        coords = GridCoordinates(cfg_atm.grid)
        state_atm, coords_actual, ws, _, _ = init_simulation(cfg_atm)

        state_atm.atm.P_surf = 1.0e5
        state_atm.atm.T_surf_eq = 300.0
        state_atm.atm.tau_LW = 0.5

        _, P_amb_disk, _ = compute_ambient_conditions(state_atm.timesum, cfg_atm.disk)
        prepare_step_ambient!(state_atm, coords_actual, cfg_atm, ws)
        @test isapprox(state_atm.accumulators.P_amb, P_amb_disk + 1.0e5; rtol=1e-12)

        apply_surface_radiation!(state_atm, coords_actual, cfg_atm)
        @test all(isfinite, state_atm.grids.KX)
        @test all(isfinite, state_atm.grids.KY)
    end

    @testset "Refractory Organic Pyrolysis Step" begin
        cfg = default_config()
        coords = GridCoordinates(cfg.grid)
        state, coords_actual, ws, _, _ = init_simulation(cfg)

        res_inactive = Erebus.update_pyrolysis!(state, coords_actual, cfg)
        @test res_inactive === nothing

        marknum = length(state.markers)
        cfg_pyro = SimulationConfig(;
            (
                f => (
                    if f === :refractory
                        RefractoryConfig(active=true, kinetics_active=true)
                    else
                        getfield(cfg, f)
                    end
                ) for f in fieldnames(SimulationConfig)
            )...,
        )
        hcnspo_props = setup_marker_hcnspo_properties(
            marknum, cfg_pyro.volatile_mixture, cfg_pyro.refractory
        )
        hcnspo_props.X_refr_C_m .= 0.01
        hcnspo_props.X_refr_N_m .= 0.001
        hcnspo_props.X_refr_H_m .= 0.002
        vols_grp = VolatilesGroup(
            zeros(Float64, marknum),
            zeros(Float64, marknum),
            zeros(Float64, marknum),
            zeros(Float64, marknum),
            zeros(Float64, marknum),
            zeros(Float64, marknum),
        )
        state_pyro = SimulationState(
            state.grids,
            MarkerArrays(state.markers.core, (; hcnspo=hcnspo_props, volatiles=vols_grp)),
            state.accumulators,
            state.transfers,
            state.atm,
            state.rng,
            state.timer,
            state.timestep,
            state.dt,
            state.timesum,
        )
        state_pyro.markers.core.tkm .= 800.0
        res_active = Erebus.update_pyrolysis!(state_pyro, coords_actual, cfg_pyro)
        @test res_active isa Matrix{Float64}
        @test size(res_active) == (coords_actual.Ny1, coords_actual.Nx1)
        @test all(isfinite, res_active)
        @test any(>(0.0), vols_grp.XCm)
        @test any(>(0.0), vols_grp.XNm)
        @test any(>(0.0), vols_grp.X_graphite_m)
    end

    @testset "Post-processing time step reads the fluid displacement limit" begin
        cfg0 = default_config()
        state0, coords0, ws0, _, _ = init_simulation(cfg0)
        # Liquid pore water (above 273 K), so the fluid moves faster than the solid
        state0.markers.core.tkm .= 300.0
        interpolate_markers_to_grid!(
            state0,
            coords0,
            cfg0;
            p2m_workspace=ws0.p2m,
            thread_buffers=ws0.thread_buffers,
            interp_arrays=ws0.interp_arrays,
        )
        solve_gravity!(state0, coords0, cfg0; F_grav=ws0.F_grav, RP=ws0.RP, SP=ws0.SP)
        apply_surface_radiation!(state0, coords0, cfg0)
        snapshot_step_start_inventories!(ws0, state0, cfg0)
        Erebus.snapshot_radiogenic_step_start!(ws0, state0, cfg0)
        update_step_start_pressures!(state0)
        Erebus.assemble_and_solve_hydromechanical!(
            state0,
            coords0,
            cfg0,
            ws0;
            titer=1,
            iplast=1,
            cur_betasolid=0.0,
            cur_betafluid=0.0,
        )

        # Same state and solution, only the 2 displacement limits differ
        function postprocess_with_limits(dxymax, dxymax_fluid; rhocpfluid=nothing)
            tkw = (; (fn => getfield(cfg0.time, fn) for fn in fieldnames(TimeConfig))...)
            tcfg = TimeConfig(; merge(tkw, (; dxymax=dxymax, dxymax_fluid=dxymax_fluid))...)
            mat = cfg0.materials
            if rhocpfluid !== nothing
                MC = typeof(mat)
                mat = MC(
                    (
                        fn === :rhocpfluidm ? rhocpfluid : getfield(mat, fn) for
                        fn in fieldnames(MC)
                    )...,
                )
            end
            cfg = SimulationConfig(
                (
                    if fn === :time
                        tcfg
                    elseif fn === :materials
                        mat
                    else
                        getfield(cfg0, fn)
                    end for fn in fieldnames(SimulationConfig)
                )...
            )
            st = deepcopy(state0)
            w = deepcopy(ws0)
            Erebus.postprocess_hydromechanical_solution!(
                st,
                coords0,
                cfg,
                w;
                titer=1,
                iplast=1,
                dt_step_initial=st.dt,
                cur_betasolid=0.0,
                cur_betafluid=0.0,
            )
            return st
        end
        # Loose limits do not cut dt, so these are the velocities of the solution
        loose = postprocess_with_limits(1.0, 1.0)
        g = loose.grids
        @test maximum(abs, g.vxf) > 1.2 * maximum(abs, g.vx)
        dt_fluid = min(coords0.dx / maximum(abs, g.vxf), coords0.dy / maximum(abs, g.vyf))
        strict = postprocess_with_limits(0.05, 1.0e-4)
        @test strict.dt ≈ 1.0e-4 * dt_fluid rtol=1e-10
        @test postprocess_with_limits(0.05, 2.0e-4).dt ≈ 2.0e-4 * dt_fluid rtol=1e-10

        # A large fluid heat capacity in the step (the grid RHOCP stays): the heat limit binds
        cf = 1.0e13
        rate = Erebus.Numerics.fluid_advection_inflow_rate(
            g.tk1, g.qxD, g.qyD, g.RHOCP, _ -> cf; coords=coords0
        )
        @test 1.0 / rate < 0.1 * loose.dt
        heat = postprocess_with_limits(
            1.0, 1.0; rhocpfluid=typeof(cfg0.materials.rhocpfluidm)((cf, cf, cf))
        )
        @test heat.dt ≈ 1.0 / rate rtol=1e-10
    end

    @testset "Hydromechanical Assembly, Postprocessing, and Outer Iterations" begin
        cfg = default_config()
        coords = GridCoordinates(cfg.grid)
        state, coords_actual, ws, _, _ = init_simulation(cfg)

        interpolate_markers_to_grid!(
            state,
            coords_actual,
            cfg;
            p2m_workspace=ws.p2m,
            thread_buffers=ws.thread_buffers,
            interp_arrays=ws.interp_arrays,
        )
        solve_gravity!(state, coords_actual, cfg; F_grav=ws.F_grav, RP=ws.RP, SP=ws.SP)
        apply_surface_radiation!(state, coords_actual, cfg)
        snapshot_step_start_inventories!(ws, state, cfg)
        Erebus.snapshot_radiogenic_step_start!(ws, state, cfg)
        update_step_start_pressures!(state)

        Erebus.assemble_and_solve_hydromechanical!(
            state,
            coords_actual,
            cfg,
            ws;
            titer=1,
            iplast=1,
            cur_betasolid=0.0,
            cur_betafluid=0.0,
        )
        @test all(isfinite, ws.S)

        res_post = Erebus.postprocess_hydromechanical_solution!(
            state,
            coords_actual,
            cfg,
            ws;
            titer=1,
            iplast=1,
            dt_step_initial=state.dt,
            cur_betasolid=0.0,
            cur_betafluid=0.0,
        )
        @test haskey(res_post, :adjustment_ok)
        @test haskey(res_post, :aphimax)
        @test haskey(res_post, :n_flips_iter)

        # A state without the metric array cannot run with the option on.
        grids_no_metric = Erebus.GridArrays(
            (
                fn === :Q_metric ? nothing : getfield(state.grids, fn) for
                fn in fieldnames(Erebus.GridArrays)
            )...,
        )
        state_no_metric = SimulationState(
            (
                fn === :grids ? grids_no_metric : getfield(state, fn) for
                fn in fieldnames(typeof(state))
            )...,
        )
        @test_throws ArgumentError Erebus.solve_thermal_energy!(
            state_no_metric, coords_actual, cfg, ws; titer=1, DHP_pyro=nothing
        )
        res_therm = Erebus.solve_thermal_energy!(
            state, coords_actual, cfg, ws; titer=1, DHP_pyro=nothing
        )
        @test haskey(res_therm, :thermochemical_converged)
        @test haskey(res_therm, :maxDTcurrent)
        @test haskey(res_therm, :dt_next)

        res_outer = Erebus.solve_thermomechanical_iterations!(
            state, coords_actual, cfg, ws; DHP_pyro=nothing, dt_step_initial=state.dt
        )
        @test haskey(res_outer, :plastic_converged)
        @test haskey(res_outer, :thermochemical_converged)
        @test res_outer.plastic_converged == true

        cfg_noreact = override_config(cfg, Dict{String,Any}("reaction.active" => false))
        fill!(state.grids.DMP, 1.0e-9)
        fill!(state.grids.DQPF, 2.0e-9)
        fill!(state.grids.DHP, 3.0)
        res_noreact = Erebus.solve_thermomechanical_iterations!(
            state,
            coords_actual,
            cfg_noreact,
            ws;
            DHP_pyro=nothing,
            dt_step_initial=state.dt,
        )
        @test res_noreact.thermochemical_converged == true
        @test all(iszero, state.grids.DMP)
        @test all(iszero, state.grids.DQPF)
        @test all(iszero, state.grids.DHP)

        cfg_rdx_melt = override_config(
            load_config("configs/test_quick.toml"),
            Dict{String,Any}(
                "redox.active" => true,
                "melting.active" => true,
                "magma_transport.active" => true,
                "magma_transport.track_depletion" => true,
                "solver.max_plastic_iterations" => 2,
            ),
        )
        st_rm, coords_rm, ws_rm, tel_rm, _ = init_simulation(cfg_rdx_melt)
        tel_rm !== nothing && close(tel_rm)
        st_rm.markers.core.tkm .= 1700.0
        interpolate_markers_to_grid!(
            st_rm,
            coords_rm,
            cfg_rdx_melt;
            p2m_workspace=ws_rm.p2m,
            thread_buffers=ws_rm.thread_buffers,
            interp_arrays=ws_rm.interp_arrays,
        )
        snapshot_step_start_inventories!(ws_rm, st_rm, cfg_rdx_melt)
        Erebus.snapshot_radiogenic_step_start!(ws_rm, st_rm, cfg_rdx_melt)
        @test ws_rm.hr_step_start_timestep == st_rm.timestep
        @test ws_rm.nFe0_m_step_start isa Vector{Float64}
        @test isapprox(
            ws_rm.nFe2_m_step_start, st_rm.markers.groups.redox.nFe2_m; rtol=1.0e-12
        )
        # Mutate redox array before outer iteration to verify titer > 1 restoration
        nfe2_saved = copy(st_rm.markers.groups.redox.nFe2_m)
        st_rm.markers.groups.redox.nFe2_m .*= 0.5
        Erebus.solve_thermomechanical_iterations!(
            st_rm,
            coords_rm,
            cfg_rdx_melt,
            ws_rm;
            DHP_pyro=nothing,
            dt_step_initial=st_rm.dt,
        )
        @test isapprox(st_rm.markers.groups.redox.nFe2_m, nfe2_saved; rtol=1.0e-12)
    end

    @testset "Coupled Physics Simulation Step Integration" begin
        cfg_base = load_config("configs/test_quick.toml")
        cfg_coupled = override_config(
            cfg_base,
            Dict{String,Any}(
                "time.n_steps" => 1,
                "volatiles.active" => true,
                "melting.active" => true,
                "coreformation.percolation_active" => true,
                "coreformation.settling_active" => true,
                "coreformation.segregation_heating" => true,
                "metal_partition.active" => true,
                "magma_transport.active" => true,
                "venting.active" => true,
                "venting.latent_cooling" => true,
            ),
        )
        state_coupled = simulation_loop(cfg_coupled)
        @test state_coupled.timestep == 1
        @test all(isfinite, state_coupled.grids.Q_seg_grid)
        @test all(isfinite, state_coupled.grids.Q_lat_grid)
    end

    @testset "Stokes-Darcy Plastic Non-Convergence Branch" begin
        cfg_base = default_config()
        coords = GridCoordinates(cfg_base.grid)
        state, coords_actual, ws, _, _ = init_simulation(cfg_base)

        cfg_noconv = SimulationConfig(;
            (
                f => (
                    if f === :solver
                        SolverConfig(max_plastic_iterations=0)
                    else
                        getfield(cfg_base, f)
                    end
                ) for f in fieldnames(SimulationConfig)
            )...,
        )

        res_noconv = Erebus.solve_stokes_darcy!(
            state, coords_actual, cfg_noconv, ws; titer=1, dt_step_initial=state.dt
        )
        @test res_noconv.plastic_converged == false
        @test iszero(res_noconv.last_plastic_residual)
        @test res_noconv.n_flips_step_total == 0
    end

    @testset "Setup Simulation Checkpoint Error Handling" begin
        cfg = default_config()
        @test_throws Erebus.CheckpointError init_simulation(
            cfg; restart_from="nonexistent_checkpoint.jld2"
        )
        @test_throws Erebus.CheckpointError init_simulation(
            cfg; restart_from="invalid_path.jld2"
        )
    end

    @testset "Unconditional Marker Pressure Interpolation and Hydrostatic Init" begin
        cfg_dry = SimulationConfig(
            grid=GridConfig(Nx=25, Ny=25, xsize=140000.0, ysize=140000.0),
            geometry=GeometryConfig(
                rplanet=50000.0,
                rcrust=50000.0,
                xcenter=70000.0,
                ycenter=70000.0,
                gravity_mode=:enclosed_mass,
            ),
            reaction=ReactionConfig(active=false),
        )
        state, coords_actual, ws, _, _ = init_simulation(cfg_dry)
        core = state.markers.core
        marknum = length(state.markers)

        # 1. Step-0 hydrostatic pressure initialization for interior and sticky-air markers
        r_p = cfg_dry.geometry.rplanet
        p_amb = max(cfg_dry.geometry.psurface, cfg_dry.disk.p_amb_disk)
        for m in 1:marknum
            tm_m = core.tm[m]
            if tm_m < 3
                rho_s = cfg_dry.materials.rhosolidm[tm_m] * (1.0 - core.phim[m])
                r_sq =
                    (core.xm[m] - cfg_dry.geometry.xcenter)^2 +
                    (core.ym[m] - cfg_dry.geometry.ycenter)^2
                p_expected =
                    p_amb +
                    (2.0 / 3.0) * pi * Erebus.G_GRAV * (rho_s^2) * max(0.0, r_p^2 - r_sq)
                @test isapprox(core.pfm0[m], p_expected; rtol=1e-10)
            else
                @test iszero(core.pfm0[m])
            end
        end

        # Verify :poisson2d branch of _init_hydrostatic_marker_pressure!
        cfg_p2d = SimulationConfig(
            grid=GridConfig(Nx=15, Ny=15, xsize=140000.0, ysize=140000.0),
            geometry=GeometryConfig(
                rplanet=50000.0,
                rcrust=50000.0,
                xcenter=70000.0,
                ycenter=70000.0,
                gravity_mode=:poisson2d,
            ),
        )
        state_p2d, coords_p2d, ws_p2d, _, _ = init_simulation(cfg_p2d)
        core_p2d = state_p2d.markers.core
        for m in 1:length(state_p2d.markers)
            tm_m = core_p2d.tm[m]
            if tm_m < 3
                rho_s = cfg_p2d.materials.rhosolidm[tm_m] * (1.0 - core_p2d.phim[m])
                r_sq =
                    (core_p2d.xm[m] - cfg_p2d.geometry.xcenter)^2 +
                    (core_p2d.ym[m] - cfg_p2d.geometry.ycenter)^2
                p_expected = p_amb + pi * Erebus.G_GRAV * (rho_s^2) * max(0.0, r_p^2 - r_sq)
                @test isapprox(core_p2d.pfm0[m], p_expected; rtol=1e-10)
            else
                @test iszero(core_p2d.pfm0[m])
            end
        end
        Erebus.interpolate_markers_to_grid!(
            state_p2d,
            coords_p2d,
            cfg_p2d;
            p2m_workspace=ws_p2d.p2m,
            thread_buffers=ws_p2d.thread_buffers,
            interp_arrays=ws_p2d.interp_arrays,
        )
        alloc_proc_grav = let st_g = state_p2d, co_g = coords_p2d, ws_g = ws_p2d
            Erebus.process_gravitational_solution!(
                ws_g.SP, st_g.grids.FI, st_g.grids.gx, st_g.grids.gy; coords=co_g
            )
            @allocated Erebus.process_gravitational_solution!(
                ws_g.SP, st_g.grids.FI, st_g.grids.gx, st_g.grids.gy; coords=co_g
            )
        end
        if Base.JLOptions().code_coverage == 0
            @test alloc_proc_grav == 0
        end
        alloc_grav_p2d =
            let st_g = state_p2d, co_g = coords_p2d, cfg_g = cfg_p2d, ws_g = ws_p2d
                Erebus.solve_gravity!(
                    st_g, co_g, cfg_g; F_grav=ws_g.F_grav, RP=ws_g.RP, SP=ws_g.SP
                )
                @allocated Erebus.solve_gravity!(
                    st_g, co_g, cfg_g; F_grav=ws_g.F_grav, RP=ws_g.RP, SP=ws_g.SP
                )
            end
        @test alloc_grav_p2d < 16384

        # Initial factorization, in-place lu! reuse, and non-finite RHS guard in _solve_direct_umfpack!
        L_tiny = Erebus.SparseMatrixCSC{Float64,Int64}(Matrix{Float64}(I, 4, 4))
        S_tiny = zeros(Float64, 4)
        R_ones = ones(Float64, 4)
        cache_tiny = Erebus._solve_direct_umfpack!(S_tiny, L_tiny, R_ones, nothing)
        @test isapprox(S_tiny, R_ones; rtol=1e-12)
        L_tiny2 = 2.0 .* L_tiny
        cache_tiny2 = Erebus._solve_direct_umfpack!(S_tiny, L_tiny2, R_ones, cache_tiny)
        @test cache_tiny2 === cache_tiny
        @test Erebus._same_csc_pattern(cache_tiny, L_tiny2)
        @test isapprox(S_tiny, fill(0.5, 4); rtol=1e-12)
        L_tiny_offdiag = copy(L_tiny2)
        L_tiny_offdiag[1, 2] = 0.5
        @test !Erebus._same_csc_pattern(cache_tiny, L_tiny_offdiag)
        cache_tiny3 = Erebus._solve_direct_umfpack!(
            S_tiny, L_tiny_offdiag, R_ones, cache_tiny2
        )
        @test isapprox(L_tiny_offdiag * S_tiny, R_ones; rtol=1e-12)
        @test Erebus._same_csc_pattern(cache_tiny3, L_tiny_offdiag)
        R_nan = fill(NaN, 4)
        @test_throws ErrorException Erebus._solve_direct_umfpack!(
            S_tiny, L_tiny, R_nan, cache_tiny
        )

        # 2. Unconditional grid-to-marker pressure update when reaction.active == false
        fill!(core.pfm0, 0.0)
        state.timestep = 1
        Erebus.prepare_step_ambient!(state, coords_actual, cfg_dry, ws)
        Erebus.interpolate_markers_to_grid!(
            state,
            coords_actual,
            cfg_dry;
            p2m_workspace=ws.p2m,
            thread_buffers=ws.thread_buffers,
            interp_arrays=ws.interp_arrays,
        )
        Erebus.solve_gravity!(
            state, coords_actual, cfg_dry; F_grav=ws.F_grav, RP=ws.RP, SP=ws.SP
        )
        Erebus.update_step_start_pressures!(state)
        sd_res = Erebus.solve_stokes_darcy!(
            state, coords_actual, cfg_dry, ws; titer=1, dt_step_initial=state.dt
        )
        @test sd_res.plastic_converged == true
        @test isapprox(maximum(core.pfm0), maximum(state.grids.pf); rtol=0.08)

        # 3. Direct branch tests for update_marker_pressures! (pr fallback, zero early return, pfcoeff, cavitation)
        pf_zero = zeros(Float64, coords_actual.Ny1, coords_actual.Nx1)
        pr_zero = zeros(Float64, coords_actual.Ny1, coords_actual.Nx1)
        pfm_test = fill(123.0, marknum)
        Erebus.update_marker_pressures!(
            pfm_test,
            core.xm,
            core.ym,
            core.tm,
            pf_zero,
            pr_zero,
            marknum;
            coords=coords_actual,
        )
        @test all(x -> isapprox(x, 123.0; rtol=1e-12), pfm_test)

        pr_const = fill(2.0e6, coords_actual.Ny1, coords_actual.Nx1)
        fill!(pfm_test, 1.0e6)
        Erebus.update_marker_pressures!(
            pfm_test,
            core.xm,
            core.ym,
            core.tm,
            pf_zero,
            pr_const,
            marknum;
            coords=coords_actual,
            pfcoeff=0.4,
        )
        for m in 1:marknum
            if core.tm[m] < 3
                @test isapprox(pfm_test[m], 0.6 * 2.0e6 + 0.4 * 1.0e6; rtol=1e-12)
            else
                @test iszero(pfm_test[m])
            end
        end

        # 4. Verify perform_thermochemical_reaction! fallback to pfm0 when grid pf is all zeros
        fill!(state.grids.pf, 0.0)
        fill!(state.grids.tk2, 500.0)
        fill!(core.pfm0, 3.0e6)
        Erebus.perform_thermochemical_reaction!(
            state.grids.DMP,
            state.grids.DHP,
            zeros(Float64, coords_actual.Ny1, coords_actual.Nx1),
            zeros(Float64, coords_actual.Ny1, coords_actual.Nx1),
            zeros(Float64, coords_actual.Ny1, coords_actual.Nx1),
            state.grids.pf,
            state.grids.tk2,
            core.tm,
            core.xm,
            core.ym,
            core.XWsolidm0,
            core.XWsolidm,
            core.phim,
            core.phinewm,
            core.pfm0,
            marknum,
            state.dt,
            1,
            1;
            coords=coords_actual,
            DQPF=state.grids.DQPF,
            DQPFSUM=state.grids.DQPFSUM,
            cfg=ReactionConfig(active=true),
        )
        for m in 1:marknum
            if core.tm[m] < 3
                @test isapprox(core.pfm0[m], 3.0e6; rtol=1e-12)
            end
        end
    end

    @testset "Stokes-Darcy Adaptive Timestep Reduction Re-Solve Consistency" begin
        cfg_sd = SimulationConfig(
            grid=GridConfig(Nx=17, Ny=17, xsize=140000.0, ysize=140000.0),
            poroelasticity=PoroelasticConfig(betasolid=2.5e-11, betafluid=4.0e-10),
            reaction=ReactionConfig(active=false),
        )
        state, coords_actual, ws, _, _ = init_simulation(cfg_sd)
        state.timestep = 1
        Erebus.prepare_step_ambient!(state, coords_actual, cfg_sd, ws)
        Erebus.interpolate_markers_to_grid!(
            state,
            coords_actual,
            cfg_sd;
            p2m_workspace=ws.p2m,
            thread_buffers=ws.thread_buffers,
            interp_arrays=ws.interp_arrays,
        )
        Erebus.solve_gravity!(
            state, coords_actual, cfg_sd; F_grav=ws.F_grav, RP=ws.RP, SP=ws.SP
        )
        Erebus.update_step_start_pressures!(state)
        Erebus.solve_stokes_darcy!(
            state, coords_actual, cfg_sd, ws; titer=1, dt_step_initial=state.dt
        )

        # Advance to step 2 where poroelasticity is active (cur_betasolid > 0)
        state.timestep = 2
        Erebus.update_step_start_pressures!(state)
        # Perturb baseline pressure pr0 so poroelastic storage (pr - pr0)/dt is non-zero
        state.grids.pr0 .+= 5.0e5
        dt_huge = 1.0e15
        state.dt = dt_huge
        sd_res2 = Erebus.solve_stokes_darcy!(
            state, coords_actual, cfg_sd, ws; titer=1, dt_step_initial=dt_huge
        )
        @test sd_res2.plastic_converged == true
        dt_longest_sec = cfg_sd.time.dt_longest * cfg_sd.time.yearlength
        @test isapprox(min(state.dt, dt_longest_sec), state.dt; rtol=1e-12)

        # Verify discrete solid continuity equation on interior cells (3:Ny-1, 3:Nx-1)
        g = state.grids
        dx_val = coords_actual.dx
        dy_val = coords_actual.dy
        for j in 3:(coords_actual.Nx - 1), i in 3:(coords_actual.Ny - 1)
            div_vs =
                (g.vx[i, j] - g.vx[i, j - 1]) / dx_val +
                (g.vy[i, j] - g.vy[i - 1, j]) / dy_val
            rhs_cont = -g.PHI[i, j] * g.APHI[i, j] + g.DMP[i, j]
            @test isapprox(div_vs, rhs_cont; atol=1.0e-22, rtol=1.0e-10)
        end
    end
end

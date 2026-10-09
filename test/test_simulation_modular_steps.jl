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
        state_pyro = SimulationState(
            state.grids,
            MarkerArrays(state.markers.core, (; hcnspo=hcnspo_props)),
            state.accumulators,
            state.transfers,
            state.atm,
            state.rng,
            state.timer,
            state.timestep,
            state.dt,
            state.timesum,
        )
        res_active = Erebus.update_pyrolysis!(state_pyro, coords_actual, cfg_pyro)
        @test res_active isa Matrix{Float64}
        @test size(res_active) == (coords_actual.Ny1, coords_actual.Nx1)
        @test all(isfinite, res_active)
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
end

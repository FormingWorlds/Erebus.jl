
"""
Compute surface-mean oxygen fugacity (ΔIW) from near-surface silicate markers.

$(SIGNATURES)
"""
function compute_surface_mean_delta_iw(
    redox_props,
    tm,
    xm,
    ym,
    marknum::Int,
    xcenter::Float64,
    ycenter::Float64,
    rplanet_val::Float64,
    fallback::Float64,
)
    if redox_props === nothing || redox_props.deltaIW_m === nothing
        return fallback
    end
    surf_count = 0
    surf_diw = 0.0
    r_cut_sq = (0.8 * rplanet_val)^2
    r_max_sq = rplanet_val^2
    for m in 1:marknum
        if tm[m] < 3
            r_sq = (xm[m] - xcenter)^2 + (ym[m] - ycenter)^2
            if r_sq >= r_cut_sq && r_sq <= r_max_sq
                surf_diw += redox_props.deltaIW_m[m]
                surf_count += 1
            end
        end
    end
    if surf_count > 0
        return surf_diw / surf_count
    else
        delta_iw_arr = redox_props.deltaIW_m
        return isempty(delta_iw_arr) ? fallback : sum(delta_iw_arr) / length(delta_iw_arr)
    end
end

"""
Compute multi-species surface volatile venting rates [kg/s].

$(SIGNATURES)

# Arguments
- `cfg`: Simulation configuration.
- `delta_m_vent_3d`: 3D-equivalent pore fluid water mass vented from surface Darcy sink [kg].
- `vented_vols`: Drained mobile mineral volatiles NamedTuple or nothing.
- `dt`: Current timestep [s].
- `P_surf`: Surface boundary pressure [Pa].
- `T_surf`: Surface boundary temperature [K].
- `redox_props`: Redox properties struct or nothing.
- `marknum`: Number of active markers.
- `tm`: Marker material type array.
- `xm`: Marker x coordinates [m].
- `ym`: Marker y coordinates [m].
- `rplanet_val`: Planetesimal radius [m].
- `xcenter_val`: Planetesimal center x coordinate [m] (default: 0.0).
- `ycenter_val`: Planetesimal center y coordinate [m] (default: 0.0).

# Returns
- `Dict{Symbol,Float64}` mapping volatile species to surface venting rates [kg/s].
"""
function compute_surface_venting_rates(
    cfg::SimulationConfig,
    delta_m_vent_3d::Float64,
    vented_vols,
    dt::Float64,
    P_surf::Float64,
    T_surf::Float64,
    redox_props,
    marknum::Int,
    tm,
    xm,
    ym,
    rplanet_val::Float64,
    xcenter_val::Float64=0.0,
    ycenter_val::Float64=0.0,
)
    vent_rates = Dict{Symbol,Float64}(sp => 0.0 for sp in cfg.escape.species_list)
    (
        dt <= 0.0 || (
            !cfg.venting.active && (
                vented_vols === nothing ||
                !cfg.retention.active ||
                !cfg.retention.venting_drainage_active
            )
        )
    ) && return vent_rates

    drain_on = (
        vented_vols !== nothing &&
        cfg.retention.active &&
        cfg.retention.venting_drainage_active
    )
    m_pore_H2O = cfg.venting.active ? delta_m_vent_3d : 0.0
    m_mineral_H2O = if drain_on
        (
            if hasproperty(vented_vols, :M_vent_H2O_3d)
                vented_vols.M_vent_H2O_3d
            else
                throw(ArgumentError("vented_vols must contain :M_vent_H2O_3d"))
            end
        )
    else
        0.0
    end
    m_H2O_step = m_pore_H2O + m_mineral_H2O

    m_C_step = if drain_on
        (
            if hasproperty(vented_vols, :M_vent_C_3d)
                vented_vols.M_vent_C_3d
            else
                throw(ArgumentError("vented_vols must contain :M_vent_C_3d"))
            end
        )
    else
        0.0
    end
    m_N_step = if drain_on
        (
            if hasproperty(vented_vols, :M_vent_N_3d)
                vented_vols.M_vent_N_3d
            else
                throw(ArgumentError("vented_vols must contain :M_vent_N_3d"))
            end
        )
    else
        0.0
    end
    m_S_step = if drain_on
        (
            if hasproperty(vented_vols, :M_vent_S_3d)
                vented_vols.M_vent_S_3d
            else
                throw(ArgumentError("vented_vols must contain :M_vent_S_3d"))
            end
        )
    else
        0.0
    end

    (m_H2O_step + m_C_step + m_N_step + m_S_step) <= 0.0 && return vent_rates

    if cfg.volatiles.speciation_active
        fO2_delta_IW_vent = compute_surface_mean_delta_iw(
            redox_props,
            tm,
            xm,
            ym,
            marknum,
            xcenter_val,
            ycenter_val,
            rplanet_val,
            cfg.volatiles.fO2_delta_IW,
        )
        fO2_delta_IW_vent = clamp(fO2_delta_IW_vent, -50.0, 50.0)
        P_surf_eval = max(P_surf, 1.0)
        T_surf_eval = max(T_surf, 273.15)

        spec_dict = speciate_vented_volatiles(
            m_H2O_step,
            m_C_step,
            m_N_step,
            m_S_step,
            P_surf_eval,
            T_surf_eval,
            fO2_delta_IW_vent;
            graphite_saturation=cfg.volatiles.graphite_saturation,
        )
        for (sp, m_sp) in spec_dict
            vent_rates[sp] = get(vent_rates, sp, 0.0) + m_sp / dt
        end
    else
        vent_sp = cfg.venting.species
        if vent_sp === :H2
            vent_rates[:H2] =
                get(vent_rates, :H2, 0.0) + (m_pore_H2O * ((2.0 * M_H) / M_H2O)) / dt
        else
            vent_rates[vent_sp] = get(vent_rates, vent_sp, 0.0) + m_pore_H2O / dt
        end
        vent_rates[:H2O] = get(vent_rates, :H2O, 0.0) + m_mineral_H2O / dt
        # Stoichiometric conversion: elemental C to CO2 (M_CO2 / M_C)
        vent_rates[:CO2] = get(vent_rates, :CO2, 0.0) + (m_C_step * (M_CO2 / M_C)) / dt
        vent_rates[:N2] = get(vent_rates, :N2, 0.0) + m_N_step / dt
        # Stoichiometric conversion: elemental S to H2S (M_H2S / M_S)
        vent_rates[:H2S] = get(vent_rates, :H2S, 0.0) + (m_S_step * (M_H2S / M_S)) / dt
    end

    return vent_rates
end

"""
Main simulation loop: run calculations with timestepping.

$(SIGNATURES)

# Operator-Splitting Sequence
1. Accretion & mass addition (`accrete!`)
2. Radionuclide decay heating (`radiogenic_heating!`)
3. Marker-to-mesh interpolation (`interpolate_markers_to_grid!`)
4. Gravitational field solve (`solve_gravity!`)
5. Hydro-mechanical Stokes-Darcy solve (`assemble_hydromechanical_lse!`)
6. Thermal energy & convection solve
7. Fluid-silicate reactions, venting & degassing (`vent_and_degas!`)
8. Coupled atmosphere & escape evolution (`evolve_atmosphere!`)
9. Marker advection & pressure backtracking (`advect_markers!`)
10. Marker replenishment & weight updates (`replenish!`)
11. Telescoping domain doubling (`telescope_domain!`)

# Details
- `output_path`: Absolute path where to save simulation output files
- `restart_from`: Optional path to checkpoint JLD2 file to resume from
- `force_restart_config`: Allow configuration override on restart

# Returns
- `SimulationState`: Complete simulation state container
"""
function simulation_loop(
    cfg::SimulationConfig=default_config();
    output_path::String=cfg.output.output_dir,
    restart_from::AbstractString=cfg.output.restart_from,
    force_restart_config::Bool=false,
)
    if cfg.mpi.enable
        error("Distributed execution of simulation_loop is not available in this release.")
    end

    (state, coords, ws, telemetry_io, start_step_val) = init_simulation(
        cfg;
        output_path=output_path,
        restart_from=restart_from,
        force_restart_config=force_restart_config,
    )
    coords_ref = Ref(coords)
    dt_step_target = state.dt
    dt_reduced_by_maxDT = false
    dt_next = state.dt
    dt_longest = cfg.time.dt_longest * cfg.time.yearlength
    last_timestep = start_step_val

    p = Progress(
        cfg.time.n_steps;
        showspeed=true,
        dt=0.5,
        barglyphs=BarGlyphs('|', '█', ['▁', '▂', '▃', '▄', '▅', '▆', '▇'], ' ', '|'),
        barlen=10,
    )

    try
        for timestep in start_step_val:cfg.time.n_steps
            last_timestep = timestep
            state.timestep = timestep
            @timeit state.timer "step" begin
                timestep_begin = Dates.now()
                num_dt_reductions = 0
                snapshot = snapshot_step_state(state, coords_ref[], ws)
                dt_aphimax_step_max = 0.0
                n_flips_last = 0
                n_flips_step_total = 0

                while true
                    # Determine target timestep duration for this attempt
                    dt_step_target = dt_reduced_by_maxDT ? min(state.dt, dt_longest) : min(state.dt * cfg.time.dtcoefup, dt_longest)
                    state.dt = dt_step_target
                    dt_reduced_by_maxDT = false
                    dt_step_initial = state.dt

                    # Step preparation: ambient conditions and sticky air
                    prepare_step_ambient!(state, coords_ref[], cfg, ws)

                    # Step 1: Accretion
                    accrete!(state, coords_ref[], cfg)

                    # Step 2: Radioactive heating and pyrolysis
                    radiogenic_heating!(state, coords_ref[], cfg)
                    DHP_pyro = update_pyrolysis!(state, coords_ref[], cfg)

                    # Step 3: P2M interpolation
                    interpolate_markers_to_grid!(
                        state, coords_ref[], cfg;
                        p2m_workspace=ws.p2m,
                        thread_buffers=ws.thread_buffers,
                        interp_arrays=ws.interp_arrays,
                    )

                    # Step 4: Gravity solve
                    solve_gravity!(
                        state, coords_ref[], cfg;
                        F_grav=ws.F_grav,
                        RP=ws.RP,
                        SP=ws.SP,
                    )

                    # Step 4b: Surface radiation boundary condition
                    apply_surface_radiation!(state, coords_ref[], cfg)

                    # Step 4c: Step-start inventory snapshots and pressure baselines for outer iterations
                    snapshot_step_start_inventories!(ws, state, cfg)
                    update_step_start_pressures!(state)

                    # Steps 5 & 6: Thermomechanical outer iteration loop
                    res = solve_thermomechanical_iterations!(
                        state, coords_ref[], cfg, ws;
                        DHP_pyro=DHP_pyro,
                        dt_step_initial=dt_step_initial,
                    )
                    dt_next = res.dt_next
                    dt_reduced_by_maxDT = res.dt_reduced_by_maxDT
                    dt_aphimax_step_max = res.dt_aphimax_step_max
                    n_flips_last = res.n_flips_last
                    n_flips_step_total = res.n_flips_step_total

                    # Check plastic convergence for retry
                    if !res.plastic_converged
                        num_dt_reductions += 1
                        if num_dt_reductions > cfg.solver.max_dt_reductions
                            throw(PlasticConvergenceError(
                                timestep, res.last_plastic_residual, state.dt,
                                "Plastic iterations failed to converge after $(cfg.solver.max_dt_reductions) dt reductions",
                            ))
                        end
                        restore_step_state!(state, coords_ref, ws, snapshot)
                        dt_step_target = min(dt_step_target / 2.0, dt_next)
                        state.dt = dt_step_target
                        continue
                    end

                    break # Plastic iterations converged
                end

                # Post-convergence: marker viscoelasticity and subgrid diffusion
                diffuse_and_update_markers!(state, coords_ref[], cfg, ws)

                # Step 7: Venting and degassing
                vent_res = vent_and_degas!(
                    state, coords_ref[], cfg;
                    Fm_step_start=ws.Fm_step_start,
                )

                # Step 8: Coupled surface atmosphere and escape
                evolve_atmosphere!(state, coords_ref[], cfg; vent_degas_result=vent_res)

                # Step 9: Marker advection
                advect_markers!(state, coords_ref[], cfg)

                # Step 10: Replenishment
                replenish!(
                    state, coords_ref[], cfg;
                    mdis=ws.mdis,
                    mnum=ws.mnum,
                    randomized=random_markers,
                    step_start_buffers=(ws.Xfe_bulk_step_start, ws.Xfem_step_start, ws.F_extract_m_step_start, ws.Fm_step_start, ws.Xfe_H_m_step_start, ws.Xfe_C_m_step_start, ws.Xfe_N_m_step_start, ws.Xfe_S_m_step_start),
                )

                # Step 11: Telescoping domain doubling (post-convergence)
                telescope_domain!(state, coords_ref, cfg, ws)

                # Time progression
                state.timesum += state.dt

                # Step diagnostics, telemetry, checkpoints, progress
                advance_step_diagnostics!(
                    state, coords_ref[], cfg, ws;
                    output_path=output_path,
                    telemetry_io=telemetry_io,
                    dt_aphimax_step_max=dt_aphimax_step_max,
                    n_flips_last=n_flips_last,
                    n_flips_step_total=n_flips_step_total,
                    timestep_begin=timestep_begin,
                    progress_bar=p,
                )

                if dt_reduced_by_maxDT
                    state.dt = dt_next
                end
            end

            if state.timesum > cfg.time.endtime * cfg.time.yearlength
                break
            end
        end
    finally
        if telemetry_io !== nothing
            close(telemetry_io)
        end
    end

    return state
end

"""
Simulation loop overload for path or configuration file input.

$(SIGNATURES)

# Arguments
- `path_or_dir`: Path to `.toml` configuration file, or output directory path.
- `output_path`: Optional output path override.
"""
function simulation_loop(path_or_dir::String; output_path::String="")
    if endswith(path_or_dir, ".toml")
        cfg = load_config(path_or_dir)
        actual_output = isempty(output_path) ? cfg.output.output_dir : output_path
        return simulation_loop(cfg; output_path=actual_output)
    else
        actual_output = isempty(output_path) ? path_or_dir : output_path
        return simulation_loop(default_config(); output_path=actual_output)
    end
end

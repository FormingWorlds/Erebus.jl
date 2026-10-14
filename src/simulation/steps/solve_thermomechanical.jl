# Thermomechanical iteration driver coordinating Stokes-Darcy and thermal energy solves

"""
Restore marker redox speciation and bulk metal fraction from step-start snapshots on outer iterations.

$(SIGNATURES)

# Parameters
- `state`: Simulation state container.
- `cfg`: Simulation configuration struct.
- `ws`: Preallocated simulation workspaces.

# Returns
- `Nothing`: Restores marker redox vectors and `Xfe_bulk` in place.
"""
function _restore_step_start_redox!(
    state::SimulationState, cfg::SimulationConfig, ws::SimulationWorkspaces
)::Nothing
    if cfg.redox.active && haskey(state.markers.groups, :redox)
        redox = state.markers.groups.redox
        _sync_step_start_vec!(redox.nFe0_m, ws.nFe0_m_step_start)
        _sync_step_start_vec!(redox.nFe2_m, ws.nFe2_m_step_start)
        _sync_step_start_vec!(redox.nFe3_m, ws.nFe3_m_step_start)
        _sync_step_start_vec!(redox.deltaIW_m, ws.deltaIW_m_step_start)
        if haskey(state.markers.groups, :metal)
            _sync_step_start_vec!(
                state.markers.groups.metal.Xfe_bulk, ws.Xfe_bulk_step_start
            )
        end
    end
    return nothing
end

"""
Execute one outer iteration of fluid-silicate thermochemical hydration/dehydration or zero source grids.

$(SIGNATURES)

# Parameters
- `state`: Simulation state container.
- `coords`: Grid coordinates struct.
- `cfg`: Simulation configuration struct.
- `ws`: Preallocated simulation workspaces.
- `titer`: Current outer thermomechanical iteration index.

# Returns
- `Nothing`: Mutates `state.grids.DMP`, `state.grids.DHP`, and `state.grids.DQPF` in place.
"""
function _step_thermochemical_reaction!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig,
    ws::SimulationWorkspaces,
    titer::Int,
)::Nothing
    if cfg.reaction.active
        core = state.markers.core
        copyto!(core.phinewm, core.phim)
        DMPSUM = ws.interp_arrays[30]
        DHPSUM = ws.interp_arrays[31]
        WTPSUM = ws.interp_arrays[33]
        perform_thermochemical_reaction!(
            state.grids.DMP,
            state.grids.DHP,
            DMPSUM,
            DHPSUM,
            WTPSUM,
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
            length(state.markers),
            state.dt,
            state.timestep,
            titer;
            coords=coords,
            DQPF=state.grids.DQPF,
            DQPFSUM=state.grids.DQPFSUM,
            cfg=cfg.reaction,
        )
    else
        fill!(state.grids.DMP, 0.0)
        fill!(state.grids.DHP, 0.0)
        fill!(state.grids.DQPF, 0.0)
    end
    return nothing
end

"""
Coordinate thermochemical iterations, Stokes-Darcy iterations, and thermal energy updates.

$(SIGNATURES)
"""
function solve_thermomechanical_iterations!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig,
    ws::SimulationWorkspaces;
    DHP_pyro=nothing,
    dt_step_initial::Float64,
)
    dt_aphimax_step_max = 0.0
    plastic_converged = true
    thermochemical_converged = true
    last_plastic_residual = 0.0
    dt_reduced_by_maxDT = false
    dt_next = state.dt
    n_flips_last = 0
    n_flips_step_total = 0

    max_iters = cfg.solver.max_plastic_iterations
    for titer in 1:max_iters
        titer > 1 && _restore_step_start_redox!(state, cfg, ws)
        _step_thermochemical_reaction!(state, coords, cfg, ws, titer)
        if DHP_pyro !== nothing && size(DHP_pyro) == size(state.grids.DHP)
            state.grids.DHP .+= DHP_pyro
        end

        sd_res = solve_stokes_darcy!(
            state, coords, cfg, ws; titer=titer, dt_step_initial=dt_step_initial
        )
        dt_aphimax_step_max = max(dt_aphimax_step_max, sd_res.dt_aphimax_step_max)
        n_flips_last = sd_res.n_flips_last
        n_flips_step_total += sd_res.n_flips_step_total

        if !sd_res.plastic_converged
            plastic_converged = false
            last_plastic_residual = sd_res.last_plastic_residual
            break
        end

        th_res = solve_thermal_energy!(
            state, coords, cfg, ws; titer=titer, DHP_pyro=DHP_pyro
        )
        dt_next = th_res.dt_next
        dt_reduced_by_maxDT = dt_reduced_by_maxDT || th_res.dt_reduced_by_maxDT

        if th_res.thermochemical_converged
            thermochemical_converged = true
            break
        elseif titer == max_iters
            thermochemical_converged = false
            break
        else
            state.dt = dt_next
        end
    end

    if plastic_converged
        update_marker_pressures!(
            state.markers.core.pfm0,
            state.markers.core.xm,
            state.markers.core.ym,
            state.markers.core.tm,
            state.grids.pf,
            state.grids.pr,
            length(state.markers);
            coords=coords,
            p_cavitation=cfg.reaction.p_cavitation,
            pfcoeff=0.0,
        )
    end

    return (;
        plastic_converged=plastic_converged,
        last_plastic_residual=last_plastic_residual,
        thermochemical_converged=thermochemical_converged,
        dt_next=dt_next,
        dt_reduced_by_maxDT=dt_reduced_by_maxDT,
        dt_aphimax_step_max=dt_aphimax_step_max,
        n_flips_last=n_flips_last,
        n_flips_step_total=n_flips_step_total,
    )
end

# Thermomechanical iteration driver coordinating Stokes-Darcy and thermal energy solves

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
        if cfg.reaction.active
            core = state.markers.core
            (
                ETA0SUM, ETASUM, GGGSUM, SXYSUM, COHSUM, TENSUM, FRISUM, WTSUM,
                RHOXSUM, RHOFXSUM, KXSUM, PHIXSUM, RXSUM, WTXSUM,
                RHOYSUM, RHOFYSUM, KYSUM, PHIYSUM, RYSUM, WTYSUM,
                RHOSUM, RHOCPSUM, ALPHASUM, ALPHAFSUM, HRSUM, GGGPSUM,
                SXXSUM, TKSUM, PHISUM, DMPSUM, DHPSUM, XWSSUM, WTPSUM,
            ) = ws.interp_arrays

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
            fill!(state.grids.DHP, 0.0)
        end

        if DHP_pyro !== nothing
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

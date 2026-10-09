# Diagnostics, streaming telemetry, checkpoint serialization, and progress reporting

"""
Evaluate diagnostics, write telemetry rows, and serialize simulation checkpoints.

$(SIGNATURES)
"""
function advance_step_diagnostics!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig,
    ws::SimulationWorkspaces;
    output_path::String,
    telemetry_io=nothing,
    dt_aphimax_step_max::Float64=0.0,
    n_flips_last::Int=0,
    n_flips_step_total::Int=0,
    timestep_begin=nothing,
    progress_bar=nothing,
)
    timestep = state.timestep
    n_steps_val = cfg.time.n_steps
    savematstep_val = cfg.output.savematstep

    should_save_snapshot = if cfg.output.mode in (:snapshots, :both)
        savematstep_val > 0 && (timestep % savematstep_val == 0)
    elseif cfg.output.mode == :telemetry
        cfg.output.save_final && (timestep == n_steps_val)
    else
        false
    end

    need_telemetry = (
        telemetry_io !== nothing &&
        (
            (cfg.output.telemetrystep > 0 && timestep % cfg.output.telemetrystep == 0) ||
            timestep == n_steps_val
        )
    )

    if (cfg.metal_partition.active && state.markers.Xfe_bulk !== nothing) &&
        (need_telemetry || should_save_snapshot || cfg.hydrothermal.active)
        core_budgets = compute_core_volatile_budgets(
            state.markers.xm,
            state.markers.ym,
            state.markers.tm,
            state.markers.Xfe_bulk,
            state.markers.Xfe_H_m,
            state.markers.Xfe_C_m,
            state.markers.Xfe_N_m,
            state.markers.Xfe_S_m,
            length(state.markers);
            coords=coords,
            w3d_m=state.markers.w3d_m,
            xcenter=state.accumulators.xcenter,
            ycenter=state.accumulators.ycenter,
            rplanet=state.accumulators.rplanet,
            rho_metal=cfg.coreformation.rho_metal,
            core_radius_fraction=cfg.metal_partition.core_radius_fraction,
            phi_core_threshold=cfg.metal_partition.phi_core_threshold,
        )
        state.accumulators.core_budgets = core_budgets
        if (core_budgets !== nothing && core_budgets.M_core_metal > 0.0)
            state.accumulators.rcore =
                (
                    3.0 * core_budgets.M_core_metal /
                    (4.0 * π * cfg.coreformation.rho_metal)
                )^(1.0 / 3.0)
        else
            state.accumulators.rcore = 0.0
        end
    end

    if need_telemetry
        core_radius_current = state.accumulators.rcore
        tk2 = state.grids.tk2
        PHI = state.grids.PHI
        Fm = state.markers.Fm
        M_vent_all = (
            state.accumulators.M_vent_H2O_total +
            state.accumulators.M_vent_C_total +
            state.accumulators.M_vent_N_total +
            state.accumulators.M_vent_S_total
        )
        stream_telemetry_row!(
            telemetry_io,
            timestep,
            s_to_Ma(state.timesum; yearlength=cfg.time.yearlength),
            state.dt / cfg.time.yearlength,
            state.accumulators.rplanet,
            core_radius_current,
            maximum(tk2),
            sum(tk2) / length(tk2),
            maximum(PHI),
            sum(PHI) / length(PHI),
            M_vent_all,
            state.accumulators.M_vent_H2O_total,
            state.accumulators.M_atm_total,
            state.accumulators.M_escaped_total,
            Fm !== nothing ? maximum(Fm) : 0.0,
            Fm !== nothing ? (sum(Fm) / length(Fm)) : 0.0,
            dt_aphimax_step_max,
            n_flips_last,
            n_flips_step_total,
        )
    end

    if should_save_snapshot
        if cfg.phase_tracking.active &&
            cfg.phase_tracking.track_regional_modes &&
            state.markers.Xfe_bulk !== nothing
            state.accumulators.regional_mineral_modes = compute_regional_mineral_modes(
                state.markers.xm,
                state.markers.ym,
                state.markers.tm,
                state.markers.tkm,
                state.markers.Xfe_bulk,
                state.markers.Xfe_S_m,
                state.markers.Xfe_C_m,
                state.markers.Xfe_N_m,
                length(state.markers);
                coords=coords,
                w3d_m=state.markers.w3d_m,
                cfg=cfg.phase_tracking,
                xcenter=state.accumulators.xcenter,
                ycenter=state.accumulators.ycenter,
                rplanet=state.accumulators.rplanet,
                rho_metal=cfg.coreformation.rho_metal,
            )
        end
        save_state(
            output_path,
            state,
            coords,
            cfg,
        )
    end

    maxT = maximum(state.grids.tk2)
    marknum = length(state.markers)
    if timestep_begin !== nothing
        duration_str = Dates.canonicalize(Dates.CompoundPeriod(Dates.now() - timestep_begin))
        @info "timestep $timestep computed in $duration_str"
    else
        @info "timestep $timestep computed"
    end
    @info "total time = $(s_to_Ma(state.timesum; yearlength=cfg.time.yearlength)) Ma"
    @info "markers in use = $marknum"
    @info "max T = $maxT K"

    if progress_bar !== nothing
        endtime_val = cfg.time.endtime * cfg.time.yearlength
        showvals = () -> [
            (:timestep, timestep),
            (:markers, marknum),
            (:maxT_K, maxT),
            (:dt_s, state.dt),
            (:timesum_Ma, s_to_Ma(state.timesum; yearlength=cfg.time.yearlength)),
            (:to_go_Ma, s_to_Ma(endtime_val - state.timesum; yearlength=cfg.time.yearlength)),
        ]
        next!(progress_bar; showvalues=showvals)
    end

    return nothing
end

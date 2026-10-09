# Step preparation, ambient conditions, and baseline pressure updates

"""
Prepare ambient conditions, update sticky air markers, and reset interpolation buffers.

$(SIGNATURES)
"""
function prepare_step_ambient!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig,
    ws::Union{Nothing,SimulationWorkspaces}=nothing,
)
    if ws !== nothing
        (
            ETA0SUM, ETASUM, GGGSUM, SXYSUM, COHSUM, TENSUM, FRISUM, WTSUM,
            RHOXSUM, RHOFXSUM, KXSUM, PHIXSUM, RXSUM, WTXSUM,
            RHOYSUM, RHOFYSUM, KYSUM, PHIYSUM, RYSUM, WTYSUM,
            RHOSUM, RHOCPSUM, ALPHASUM, ALPHAFSUM, HRSUM, GGGPSUM,
            SXXSUM, TKSUM, PHISUM, DMPSUM, DHPSUM, XWSSUM, WTPSUM,
        ) = ws.interp_arrays

        reset_interpolated_properties!(
            ETA0SUM, ETASUM, GGGSUM, SXYSUM, COHSUM, TENSUM, FRISUM, WTSUM,
            RHOXSUM, RHOFXSUM, KXSUM, PHIXSUM, RXSUM, WTXSUM,
            RHOYSUM, RHOFYSUM, KYSUM, PHIYSUM, RYSUM, WTYSUM,
            RHOSUM, RHOCPSUM, ALPHASUM, ALPHAFSUM, HRSUM, GGGPSUM,
            SXXSUM, TKSUM, PHISUM, WTPSUM,
        )
    end

    T_amb, P_amb, _ = compute_ambient_conditions(state.timesum, cfg.disk)
    isfinite(T_amb) || throw(DomainError(T_amb, "Ambient disk temperature must be finite, got $T_amb"))

    P_atm = if cfg.atmosphere.active && state.atm !== nothing
        state.atm.P_surf
    elseif cfg.escape.active
        compute_surface_atmospheric_pressure(
            state.accumulators.M_atm_total,
            state.accumulators.M_planet_val,
            state.accumulators.rplanet,
        )
    else
        0.0
    end
    P_amb_eff = P_amb + P_atm
    state.accumulators.P_amb = P_amb_eff

    if cfg.disk.enabled || cfg.thermodynamics.surface_radiation
        core = state.markers.core
        marknum = length(state.markers)
        @threads :dynamic for m in 1:marknum
            if core.tm[m] >= 3
                core.tkm[m] = T_amb
            end
        end
    end

    return nothing
end

"""
Apply radiative thermal surface boundary condition to grid conductivity and temperature.

$(SIGNATURES)
"""
function apply_surface_radiation!(
    state::SimulationState, coords::GridCoordinates, cfg::SimulationConfig
)
    if cfg.thermodynamics.surface_radiation
        T_amb, _, _ = compute_ambient_conditions(state.timesum, cfg.disk)
        tau_LW_val = (cfg.atmosphere.active && state.atm !== nothing) ? state.atm.tau_LW : 0.0
        T_amb_rad = if cfg.atmosphere.active && state.atm !== nothing && state.atm.T_surf_eq > 0.0
            state.atm.T_surf_eq
        else
            T_amb
        end
        apply_radiative_surface_boundary!(
            state.grids.KX,
            state.grids.KY,
            state.grids.tk1,
            coords,
            state.accumulators.rplanet,
            state.accumulators.xcenter,
            state.accumulators.ycenter,
            T_amb_rad;
            emissivity=cfg.thermodynamics.emissivity,
            sigma_sb=cfg.thermodynamics.sigma_sb,
            marker_property_mode=marker_property_mode,
            phi=cfg.thermodynamics.phim0,
            kfluid=cfg.materials.kfluidm[2],
            tau_LW=tau_LW_val,
        )
    end
    return nothing
end

"""
Snapshot metal and magma inventories at attempt start for rollback safety.

$(SIGNATURES)
"""
function snapshot_step_start_inventories!(
    ws::SimulationWorkspaces, state::SimulationState, cfg::SimulationConfig
)
    coreformation_active_val =
        cfg.coreformation.percolation_active ||
        cfg.coreformation.settling_active ||
        cfg.metal_partition.active

    if coreformation_active_val && haskey(state.markers.groups, :metal)
        metal = state.markers.groups.metal
        if ws.Xfe_bulk_step_start !== nothing
            if length(ws.Xfe_bulk_step_start) != length(metal.Xfe_bulk)
                resize!(ws.Xfe_bulk_step_start, length(metal.Xfe_bulk))
            end
            copyto!(ws.Xfe_bulk_step_start, metal.Xfe_bulk)
        end
        if ws.Xfem_step_start !== nothing
            if length(ws.Xfem_step_start) != length(metal.Xfem)
                resize!(ws.Xfem_step_start, length(metal.Xfem))
            end
            copyto!(ws.Xfem_step_start, metal.Xfem)
        end
    end

    if cfg.metal_partition.active && haskey(state.markers.groups, :metal)
        metal = state.markers.groups.metal
        for (src, dst) in (
            (metal.Xfe_H_m, ws.Xfe_H_m_step_start),
            (metal.Xfe_C_m, ws.Xfe_C_m_step_start),
            (metal.Xfe_N_m, ws.Xfe_N_m_step_start),
            (metal.Xfe_S_m, ws.Xfe_S_m_step_start),
        )
            if src !== nothing && dst !== nothing
                if length(dst) != length(src)
                    resize!(dst, length(src))
                end
                copyto!(dst, src)
            end
        end
    end

    if cfg.magma_transport.active
        core = state.markers.core
        if ws.Fm_step_start !== nothing
            if length(ws.Fm_step_start) != length(core.Fm)
                resize!(ws.Fm_step_start, length(core.Fm))
            end
            copyto!(ws.Fm_step_start, core.Fm)
        end
        if haskey(state.markers.groups, :volatiles) && ws.F_extract_m_step_start !== nothing
            vols = state.markers.groups.volatiles
            if length(ws.F_extract_m_step_start) != length(vols.F_extract_m)
                resize!(ws.F_extract_m_step_start, length(vols.F_extract_m))
            end
            copyto!(ws.F_extract_m_step_start, vols.F_extract_m)
        end
    end

    return nothing
end

"""
Update step-start pressures on staggered grid once per timestep attempt.

$(SIGNATURES)
"""
function update_step_start_pressures!(state::SimulationState)
    state.grids.pr0 .= state.grids.pr
    state.grids.pf0 .= state.grids.pf
    state.grids.ps0 .= state.grids.ps
    return nothing
end

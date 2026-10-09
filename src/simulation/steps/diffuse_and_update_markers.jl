# Marker viscoplastic viscosity updates and subgrid stress/thermal diffusion step

"""
Interpolate viscoplastic viscosity to markers and diffuse subgrid stress and thermal increments.

$(SIGNATURES)
"""
function diffuse_and_update_markers!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig,
    ws::Union{Nothing,SimulationWorkspaces}=nothing,
)
    g = state.grids
    core = state.markers.core
    marknum = length(state.markers)

    g.DT0 .= g.DT

    @threads :dynamic for m in 1:marknum
        update_marker_viscosity!(
            m,
            core.xm,
            core.ym,
            core.tm,
            core.tkm,
            core.etatotalm,
            core.etavpm,
            g.YNY,
            g.YNY_inv_ETA;
            coords=coords,
            Fm=core.Fm,
            etamin=cfg.solver.etamin,
            etamax=cfg.solver.etamax,
            melting_active=cfg.melting.active,
            alpha_eta_val=cfg.melting.alpha_eta,
            phi_crit_val=cfg.melting.phi_crit,
            eta_melt_val=cfg.melting.eta_melt,
            tmsolidphase=cfg.thermodynamics.tmsolidphase,
            tmfluidphase=cfg.thermodynamics.tmfluidphase,
            etasolidm=cfg.materials.etasolidm,
            etasolidmm=cfg.materials.etasolidmm,
            etafluidm=cfg.materials.etafluidm,
            etafluidmm=cfg.materials.etafluidmm,
        )
    end

    interp_arrays =
        ws !== nothing ? ws.interp_arrays : setup_interpolated_properties(coords)
    (
        ETA0SUM,
        ETASUM,
        GGGSUM,
        SXYSUM,
        COHSUM,
        TENSUM,
        FRISUM,
        WTSUM,
        RHOXSUM,
        RHOFXSUM,
        KXSUM,
        PHIXSUM,
        RXSUM,
        WTXSUM,
        RHOYSUM,
        RHOFYSUM,
        KYSUM,
        PHIYSUM,
        RYSUM,
        WTYSUM,
        RHOSUM,
        RHOCPSUM,
        ALPHASUM,
        ALPHAFSUM,
        HRSUM,
        GGGPSUM,
        SXXSUM,
        TKSUM,
        PHISUM,
        DMPSUM,
        DHPSUM,
        XWSSUM,
        WTPSUM,
    ) = interp_arrays

    apply_subgrid_stress_diffusion!(
        core.xm,
        core.ym,
        core.tm,
        core.inv_gggtotalm,
        core.sxxm,
        core.sxym,
        g.SXX0,
        g.SXY0,
        g.DSXX,
        g.DSXY,
        SXXSUM,
        SXYSUM,
        WTPSUM,
        WTSUM,
        state.dt,
        marknum;
        coords=coords,
        dsubgrids=cfg.solver.dsubgrids,
    )

    update_marker_stress!(
        core.xm, core.ym, core.sxxm, core.sxym, g.DSXX, g.DSXY, marknum; coords=coords
    )

    apply_subgrid_temperature_diffusion!(
        core.xm,
        core.ym,
        core.tm,
        core.tkm,
        core.phim,
        g.tk1,
        g.DT,
        TKSUM,
        RHOCPSUM,
        state.dt,
        marknum,
        marker_property_mode;
        coords=coords,
        dsubgridt=cfg.solver.dsubgridt,
    )

    return nothing
end

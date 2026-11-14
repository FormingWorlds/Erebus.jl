# Marker viscoplastic viscosity updates and subgrid stress/thermal diffusion step

"""
Collect the out-of-plane stress arrays that the subgrid stress routines take as keywords.

$(SIGNATURES)

# Arguments
- `core`: `CoreGroup` holding the marker τ33 vector `s33m` [Pa]
- `g`: `GridArrays` holding `S330` and `DS33` at P nodes [Pa]
- `interp_arrays`: 34-tuple of interpolation accumulators from `setup_interpolated_properties`
- `on`: `cfg.geometry.spherical_flow`; with `false` every entry is `nothing`

# Returns
- `NamedTuple` with `s33m`, `S330`, `DS33` and `S33SUM`

# Throws
- `BoundsError` if `interp_arrays` has fewer than 34 elements
"""
function s33_stress_kwargs(core::CoreGroup, g::GridArrays, interp_arrays::Tuple, on::Bool)
    on || return (; s33m=nothing, S330=nothing, DS33=nothing, S33SUM=nothing)
    return (; s33m=core.s33m, S330=g.S330, DS33=g.DS33, S33SUM=interp_arrays[34])
end

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

    if Threads.nthreads() > 1
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
    else
        for m in 1:marknum
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

    s33_kw = s33_stress_kwargs(core, g, interp_arrays, cfg.geometry.spherical_flow)
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
        etavpm=core.etavpm,
        etasolidm=cfg.materials.etasolidm,
        s33_kw...,
    )

    update_marker_stress!(
        core.xm,
        core.ym,
        core.sxxm,
        core.sxym,
        g.DSXX,
        g.DSXY,
        marknum;
        coords=coords,
        s33m=s33_kw.s33m,
        DS33=s33_kw.DS33,
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

    # The first step sets markers from the grid, without an increment.
    if state.timestep > 1
        limit_marker_temperature_increment!(
            core.xm,
            core.ym,
            core.tm,
            core.tkm,
            core.rhocptotalm,
            g.DT,
            g.tk2,
            TKSUM,
            RHOCPSUM,
            marknum;
            coords=coords,
            DT_conduction=cfg.solver.dsubgridt > 0.0 ? g.DT0 : g.DT,
            w3d=cfg.geometry.spherical_metric ? core.w3d_m : nothing,
        )
    end

    return nothing
end

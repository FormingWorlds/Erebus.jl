# Marker property computation and particle-to-mesh (P2M) grid interpolation

"""
Pre-evaluate marker melt fractions and update marker redox speciation before grid interpolation.

$(SIGNATURES)

# Parameters
- `redox_props`: Marker redox property container.
- `core`: Core marker array container.
- `cfg`: Simulation configuration struct.
- `marknum`: Number of active markers.
- `Xfe_bulk`: Optional marker bulk metal volume fraction vector.
- `Xfem`: Optional marker molten metal volume fraction vector.
- `F_extract_m`: Optional marker cumulative extracted melt fraction vector.

# Returns
- `Nothing`: Mutates `core.Fm`, `redox_props`, and `Xfe_bulk` in place.
"""
function _preupdate_marker_redox!(
    redox_props,
    core,
    cfg::SimulationConfig,
    marknum::Int;
    Xfe_bulk=nothing,
    Xfem=nothing,
    F_extract_m=nothing,
)::Nothing
    if cfg.melting.active
        tm = core.tm
        tkm = core.tkm
        pfm0 = core.pfm0
        Fm = core.Fm
        for m in 1:marknum
            tm_m = tm[m]
            if tm_m < 3
                P_val = max(0.0, pfm0[m])
                F_m = compute_melt_fraction(
                    tkm[m],
                    P_val,
                    tm_m;
                    T_sol=cfg.melting.T_solidus[tm_m],
                    T_liq=cfg.melting.T_liquidus[tm_m],
                    dpdt=cfg.melting.dpdt_clapeyron,
                )
                if cfg.magma_transport.active &&
                    cfg.magma_transport.track_depletion &&
                    F_extract_m !== nothing
                    F_m = clamp(F_m - F_extract_m[m], 0.0, 1.0)
                end
                Fm[m] = F_m
            else
                Fm[m] = 0.0
            end
        end
    end
    update_marker_redox!(
        redox_props,
        core.tkm,
        core.pfm0,
        cfg.redox;
        Xfe_bulk=Xfe_bulk,
        Xfem=Xfem,
        XWsolidm=core.XWsolidm,
        Fm=core.Fm,
        tm=core.tm,
        rhosolidm=cfg.materials.rhosolidm,
        rhosolid=cfg.materials.rhosolidm[1],
        rho_metal=cfg.coreformation.rho_metal,
    )
    return nothing
end

"""
    interpolate_markers_to_grid!(state::SimulationState, coords::GridCoordinates, cfg::SimulationConfig; p2m_workspace=nothing, thread_buffers=nothing, interp_arrays=nothing, rock_sums=nothing)::Nothing

Interpolate marker properties to staggered Eulerian mesh and normalize nodal properties.

`rock_sums` (from `allocate_rock_thermal_sums`) holds the rock-marker sums of
`apply_fixed_surface_temperature!` under `thermodynamics.surface_temperature = :fixed`
and the face conductivities before the cut-cell scaling. `solve_thermal_energy!` reads
these from `ws.rock_sums`, so a time step passes `ws.rock_sums` here; with `nothing` the
pass allocates its own arrays.

# Mutates:
- `state.grids.ETA`
- `state.grids.ETA0`
- `state.grids.GGG`
- `state.grids.SXY0`
- `state.grids.COH`
- `state.grids.TEN`
- `state.grids.FRI`
- `state.grids.YNY`
- `state.grids.RHOX`
- `state.grids.RHOFX`
- `state.grids.KX`
- `state.grids.PHIX`
- `state.grids.RX`
- `state.grids.RHOY`
- `state.grids.RHOFY`
- `state.grids.KY`
- `state.grids.PHIY`
- `state.grids.RY`
- `state.grids.RHO`
- `state.grids.RHOCP`
- `state.grids.ALPHA`
- `state.grids.ALPHAF`
- `state.grids.HR`
- `state.grids.GGGP`
- `state.grids.SXX0`
- `state.grids.S330`
- `state.grids.tk1`
- `state.grids.tk2`
- `state.grids.PHI`
- `state.grids.BETAPHI`
- `state.markers.core` (constitutive properties: rhototalm, etatotalm, hrtotalm, ktotalm, etc.)
"""
function interpolate_markers_to_grid!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig;
    p2m_workspace=nothing,
    thread_buffers=nothing,
    interp_arrays=nothing,
    rock_sums=nothing,
)::Nothing
    # Empty marker set contract: return immediately with no mutations
    length(state.markers) == 0 && return nothing

    markers = state.markers
    core = markers.core
    grids = state.grids
    acc = state.accumulators
    marknum = length(markers)

    tm = core.tm
    xm = core.xm
    ym = core.ym
    tkm = core.tkm
    phim = core.phim
    pfm0 = core.pfm0
    XWsolidm0 = core.XWsolidm0
    rhototalm = core.rhototalm
    rhocptotalm = core.rhocptotalm
    etatotalm = core.etatotalm
    hrtotalm = core.hrtotalm
    ktotalm = core.ktotalm
    tkm_rhocptotalm = core.tkm_rhocptotalm
    etafluidcur_inv_kphim = core.etafluidcur_inv_kphim
    rhofluidcur = core.rhofluidcur
    alphasolidcur = core.alphasolidcur
    alphafluidcur = core.alphafluidcur
    etavpm = core.etavpm
    inv_gggtotalm = core.inv_gggtotalm
    sxym = core.sxym
    sxxm = core.sxxm
    s33m = cfg.geometry.spherical_flow ? core.s33m : nothing
    cohestotalm = core.cohestotalm
    tenstotalm = core.tenstotalm
    fricttotalm = core.fricttotalm
    Fm = core.Fm

    # Setup or reset interpolation accumulators
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
        S33SUM,
    ) = if interp_arrays !== nothing
        interp_arrays
    else
        setup_interpolated_properties(coords)
    end

    if interp_arrays !== nothing
        reset_interpolated_properties!(
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
            WTPSUM;
            S33SUM=S33SUM,
        )
    end

    # Optional marker groups
    has_metal = haskey(markers.groups, :metal)
    Xfe_bulk = has_metal ? markers.groups.metal.Xfe_bulk : nothing
    Xfem = has_metal ? markers.groups.metal.Xfem : nothing
    Xfe_H_m = has_metal ? markers.groups.metal.Xfe_H_m : nothing
    Xfe_C_m = has_metal ? markers.groups.metal.Xfe_C_m : nothing
    Xfe_N_m = has_metal ? markers.groups.metal.Xfe_N_m : nothing
    Xfe_S_m = has_metal ? markers.groups.metal.Xfe_S_m : nothing

    has_vols = haskey(markers.groups, :volatiles)
    XH2Om = has_vols ? markers.groups.volatiles.XH2Om : nothing
    XCm = has_vols ? markers.groups.volatiles.XCm : nothing
    XNm = has_vols ? markers.groups.volatiles.XNm : nothing
    XSm = has_vols ? markers.groups.volatiles.XSm : nothing
    F_extract_m = has_vols ? markers.groups.volatiles.F_extract_m : nothing

    has_phase = haskey(markers.groups, :phase)
    Xmin_troilite_m = has_phase ? markers.groups.phase.Xmin_troilite_m : nothing
    Xmin_schreibersite_m = has_phase ? markers.groups.phase.Xmin_schreibersite_m : nothing
    Xmin_cohenite_m = has_phase ? markers.groups.phase.Xmin_cohenite_m : nothing
    Xmin_graphite_m = has_phase ? markers.groups.phase.Xmin_graphite_m : nothing
    Xmin_nitride_m = has_phase ? markers.groups.phase.Xmin_nitride_m : nothing
    Xmin_metal_matrix_m = has_phase ? markers.groups.phase.Xmin_metal_matrix_m : nothing

    has_redox = haskey(markers.groups, :redox)
    redox_props = has_redox ? markers.groups.redox : nothing
    deltaIW_m = has_redox ? markers.groups.redox.deltaIW_m : nothing

    if cfg.redox.active && redox_props !== nothing
        _preupdate_marker_redox!(
            redox_props,
            core,
            cfg,
            marknum;
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            F_extract_m=F_extract_m,
        )
    end

    # compute_marker_properties! runs with compute_hr=false; hrtotalm is precomputed in radiogenic_heating!
    hrsolidm = hrfluidm = hrmetalm = @SVector zeros(3)

    thermal_buoyancy_val = cfg.thermodynamics.thermal_buoyancy
    alphafluid_val = cfg.materials.alphafluidm
    tmfluidphase_val = cfg.thermodynamics.tmfluidphase
    fluid_viscosity_mode_val = cfg.thermodynamics.fluid_viscosity_mode
    fluid_viscosity_Ea_val = cfg.thermodynamics.fluid_viscosity_Ea
    fluid_viscosity_T0_val = cfg.thermodynamics.fluid_viscosity_T0
    fluid_viscosity_eta0_val = cfg.thermodynamics.fluid_viscosity_eta0
    coreformation_active_val =
        cfg.coreformation.percolation_active || cfg.coreformation.settling_active

    phi_dry = dry_porosity(markers)
    use_tiled_p2m = cfg.solver.p2m_mode === :tiled
    use_threading =
        cfg.solver.p2m_mode === :buffered &&
        Threads.nthreads() > 1 &&
        thread_buffers !== nothing

    if use_tiled_p2m
        ws_p2m = if p2m_workspace !== nothing
            ensure_workspace_compatible(p2m_workspace, coords, marknum, cfg.solver.tile_size)
        else
            P2MTiledWorkspace(coords, marknum, cfg.solver.tile_size)
        end
        bin_markers_into_tiles!(ws_p2m, xm, ym, coords, marknum)
        hr_rates = (hrsolidm, hrfluidm, hrmetalm)
        interp_tup = (
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
            S33SUM,
        )
        for color in 1:4
            tiles = ws_p2m.tiles_by_color[color]
            if Threads.nthreads() > 1
                Threads.@threads :dynamic for t in tiles
                    _scatter_tiled_p2m_tile!(
                        t, ws_p2m, state, coords, cfg, interp_tup, hr_rates
                    )
                end
            else
                for t in tiles
                    _scatter_tiled_p2m_tile!(
                        t, ws_p2m, state, coords, cfg, interp_tup, hr_rates
                    )
                end
            end
        end
    elseif use_threading
        reset_thread_buffers!(thread_buffers)
        nchunks = length(thread_buffers)
        Threads.@threads :dynamic for c in 1:nchunks
            buf = thread_buffers[c]
            lo = (c - 1) * div(marknum, nchunks) + 1
            hi = c == nchunks ? marknum : c * div(marknum, nchunks)
            for m in lo:hi
                compute_marker_properties!(
                    m,
                    tm,
                    tkm,
                    rhototalm,
                    rhocptotalm,
                    etatotalm,
                    hrtotalm,
                    ktotalm,
                    tkm_rhocptotalm,
                    etafluidcur_inv_kphim,
                    hrsolidm,
                    hrfluidm,
                    phim,
                    XWsolidm0,
                    marker_property_mode,
                    rhofluidcur;
                    compute_hr=false,
                    thermal_buoyancy=thermal_buoyancy_val,
                    alphafluid=alphafluid_val,
                    tmfluidphase_val=tmfluidphase_val,
                    fluid_viscosity_mode=fluid_viscosity_mode_val,
                    fluid_viscosity_Ea=fluid_viscosity_Ea_val,
                    fluid_viscosity_T0=fluid_viscosity_T0_val,
                    fluid_viscosity_eta0=fluid_viscosity_eta0_val,
                    pm=pfm0,
                    Fm=Fm,
                    melting_active=cfg.melting.active,
                    magma_transport_active=cfg.magma_transport.active,
                    track_depletion=cfg.magma_transport.track_depletion,
                    F_extract_m=F_extract_m,
                    T_solidus_val=cfg.melting.T_solidus,
                    T_liquidus_val=cfg.melting.T_liquidus,
                    L_melt_val=cfg.melting.L_melt,
                    rho_melt_val=cfg.melting.rho_melt,
                    alpha_eta_val=cfg.melting.alpha_eta,
                    phi_crit_val=cfg.melting.phi_crit,
                    eta_melt_val=cfg.melting.eta_melt,
                    dpdt_clapeyron_val=cfg.melting.dpdt_clapeyron,
                    soft_turbulence=cfg.melting.soft_turbulence,
                    eta_fluid_silicate_val=cfg.melting.eta_fluid_silicate,
                    F_turb_start_val=cfg.melting.F_turb_start,
                    F_turb_end_val=cfg.melting.F_turb_end,
                    turb_exponent_val=cfg.melting.turb_exponent,
                    dT_turb_min_val=cfg.melting.dT_turb_min,
                    T_surface_ref_val=cfg.melting.T_surface_ref,
                    k_turb_cutoff_val=cfg.melting.k_turb_cutoff,
                    k_turb_floor_val=cfg.melting.k_turb_floor,
                    tmsolidphase_val=cfg.thermodynamics.tmsolidphase,
                    Xfe_bulk=Xfe_bulk,
                    Xfem=Xfem,
                    coreformation_active=coreformation_active_val,
                    hrmetalm=hrmetalm,
                    sulfur_fraction_val=cfg.coreformation.sulfur_fraction,
                    metal_density_mode_val=cfg.coreformation.metal_density_mode,
                    T_eutectic_val=cfg.coreformation.T_eutectic,
                    dT_metal_val=cfg.coreformation.dT_metal,
                    rho_metal_val=cfg.coreformation.rho_metal,
                    rho_metal_solid_val=cfg.coreformation.rho_metal_solid,
                    L_metal_val=cfg.coreformation.L_metal,
                    k_metal_val=cfg.coreformation.k_metal,
                    rhocp_metal_val=cfg.coreformation.rhocp_metal,
                    volatiles_active=cfg.volatiles.active,
                    magma_degassing_active=cfg.magma_degassing.active,
                    etamin=cfg.solver.etamin,
                    etamax=cfg.solver.etamax,
                    volatiles_cfg=cfg.volatiles,
                    retention_cfg=cfg.retention,
                    XH2Om=XH2Om,
                    XCm=XCm,
                    XNm=XNm,
                    XSm=XSm,
                    metal_partition_cfg=cfg.metal_partition,
                    Xfe_H_m=Xfe_H_m,
                    Xfe_C_m=Xfe_C_m,
                    Xfe_N_m=Xfe_N_m,
                    Xfe_S_m=Xfe_S_m,
                    phase_tracking_cfg=cfg.phase_tracking,
                    Xmin_troilite_m=Xmin_troilite_m,
                    Xmin_schreibersite_m=Xmin_schreibersite_m,
                    Xmin_cohenite_m=Xmin_cohenite_m,
                    Xmin_graphite_m=Xmin_graphite_m,
                    Xmin_nitride_m=Xmin_nitride_m,
                    Xmin_metal_matrix_m=Xmin_metal_matrix_m,
                    hydrothermal_active=cfg.hydrothermal.active,
                    hydrothermal_cfg=cfg.hydrothermal,
                    deltaIW_m=deltaIW_m,
                    xm=xm,
                    ym=ym,
                    coords=coords,
                    rplanet_val=acc.rplanet,
                    M_planet_val=acc.M_planet_val,
                    rcore_val=acc.rcore,
                    qxD_val=grids.qxD,
                    qyD_val=grids.qyD,
                    xcenter_val=acc.xcenter,
                    ycenter_val=acc.ycenter,
                    kphim0_val=cfg.materials.kphim0,
                    phim0_val=cfg.thermodynamics.phim0,
                    etafluidm_val=cfg.materials.etafluidm,
                    rhofluidm_val=cfg.materials.rhofluidm,
                    rhosolidm_val=cfg.materials.rhosolidm,
                    etasolidm_val=cfg.materials.etasolidm,
                    etasolidmm_val=cfg.materials.etasolidmm,
                    rhocpsolidm_val=cfg.materials.rhocpsolidm,
                    rhocpfluidm_val=cfg.materials.rhocpfluidm,
                    alphafluidm_val=cfg.materials.alphafluidm,
                    ksolidm_val=cfg.materials.ksolidm,
                    kfluidm_val=cfg.materials.kfluidm,
                    phimax_val=cfg.poroelasticity.phimax,
                    phi_dry=phi_dry,
                )
                @inbounds marker_to_basic_nodes!(
                    m,
                    xm[m],
                    ym[m],
                    etatotalm,
                    etavpm,
                    inv_gggtotalm,
                    sxym,
                    cohestotalm,
                    tenstotalm,
                    fricttotalm,
                    buf.ETA0SUM,
                    buf.ETASUM,
                    buf.GGGSUM,
                    buf.SXYSUM,
                    buf.COHSUM,
                    buf.TENSUM,
                    buf.FRISUM,
                    buf.WTSUM;
                    coords=coords,
                )
                @inbounds marker_to_vx_nodes!(
                    m,
                    xm[m],
                    ym[m],
                    rhototalm,
                    rhofluidcur,
                    ktotalm,
                    phim,
                    etafluidcur_inv_kphim,
                    buf.RHOXSUM,
                    buf.RHOFXSUM,
                    buf.KXSUM,
                    buf.PHIXSUM,
                    buf.RXSUM,
                    buf.WTXSUM;
                    coords=coords,
                )
                @inbounds marker_to_vy_nodes!(
                    m,
                    xm[m],
                    ym[m],
                    rhototalm,
                    rhofluidcur,
                    ktotalm,
                    phim,
                    etafluidcur_inv_kphim,
                    buf.RHOYSUM,
                    buf.RHOFYSUM,
                    buf.KYSUM,
                    buf.PHIYSUM,
                    buf.RYSUM,
                    buf.WTYSUM;
                    coords=coords,
                )
                @inbounds marker_to_p_nodes!(
                    m,
                    xm[m],
                    ym[m],
                    inv_gggtotalm,
                    sxxm,
                    rhototalm,
                    rhocptotalm,
                    alphasolidcur,
                    alphafluidcur,
                    hrtotalm,
                    phim,
                    tkm_rhocptotalm,
                    buf.GGGPSUM,
                    buf.SXXSUM,
                    buf.RHOSUM,
                    buf.RHOCPSUM,
                    buf.ALPHASUM,
                    buf.ALPHAFSUM,
                    buf.HRSUM,
                    buf.PHISUM,
                    buf.TKSUM,
                    buf.WTPSUM;
                    coords=coords,
                    s33m=s33m,
                    S33SUM=buf.S33SUM,
                )
            end
        end
        reduce_thread_buffers!(
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
            WTPSUM,
            thread_buffers;
            S33SUM=S33SUM,
        )
    else
        # Single-threaded serial scatter
        for m in 1:marknum
            compute_marker_properties!(
                m,
                tm,
                tkm,
                rhototalm,
                rhocptotalm,
                etatotalm,
                hrtotalm,
                ktotalm,
                tkm_rhocptotalm,
                etafluidcur_inv_kphim,
                hrsolidm,
                hrfluidm,
                phim,
                XWsolidm0,
                marker_property_mode,
                rhofluidcur;
                compute_hr=false,
                thermal_buoyancy=thermal_buoyancy_val,
                alphafluid=alphafluid_val,
                tmfluidphase_val=tmfluidphase_val,
                fluid_viscosity_mode=fluid_viscosity_mode_val,
                fluid_viscosity_Ea=fluid_viscosity_Ea_val,
                fluid_viscosity_T0=fluid_viscosity_T0_val,
                fluid_viscosity_eta0=fluid_viscosity_eta0_val,
                pm=pfm0,
                Fm=Fm,
                melting_active=cfg.melting.active,
                magma_transport_active=cfg.magma_transport.active,
                track_depletion=cfg.magma_transport.track_depletion,
                F_extract_m=F_extract_m,
                T_solidus_val=cfg.melting.T_solidus,
                T_liquidus_val=cfg.melting.T_liquidus,
                L_melt_val=cfg.melting.L_melt,
                rho_melt_val=cfg.melting.rho_melt,
                alpha_eta_val=cfg.melting.alpha_eta,
                phi_crit_val=cfg.melting.phi_crit,
                eta_melt_val=cfg.melting.eta_melt,
                dpdt_clapeyron_val=cfg.melting.dpdt_clapeyron,
                soft_turbulence=cfg.melting.soft_turbulence,
                eta_fluid_silicate_val=cfg.melting.eta_fluid_silicate,
                F_turb_start_val=cfg.melting.F_turb_start,
                F_turb_end_val=cfg.melting.F_turb_end,
                turb_exponent_val=cfg.melting.turb_exponent,
                dT_turb_min_val=cfg.melting.dT_turb_min,
                T_surface_ref_val=cfg.melting.T_surface_ref,
                k_turb_cutoff_val=cfg.melting.k_turb_cutoff,
                k_turb_floor_val=cfg.melting.k_turb_floor,
                tmsolidphase_val=cfg.thermodynamics.tmsolidphase,
                Xfe_bulk=Xfe_bulk,
                Xfem=Xfem,
                coreformation_active=coreformation_active_val,
                hrmetalm=hrmetalm,
                sulfur_fraction_val=cfg.coreformation.sulfur_fraction,
                metal_density_mode_val=cfg.coreformation.metal_density_mode,
                T_eutectic_val=cfg.coreformation.T_eutectic,
                dT_metal_val=cfg.coreformation.dT_metal,
                rho_metal_val=cfg.coreformation.rho_metal,
                rho_metal_solid_val=cfg.coreformation.rho_metal_solid,
                L_metal_val=cfg.coreformation.L_metal,
                k_metal_val=cfg.coreformation.k_metal,
                rhocp_metal_val=cfg.coreformation.rhocp_metal,
                volatiles_active=cfg.volatiles.active,
                magma_degassing_active=cfg.magma_degassing.active,
                etamin=cfg.solver.etamin,
                etamax=cfg.solver.etamax,
                volatiles_cfg=cfg.volatiles,
                retention_cfg=cfg.retention,
                XH2Om=XH2Om,
                XCm=XCm,
                XNm=XNm,
                XSm=XSm,
                metal_partition_cfg=cfg.metal_partition,
                Xfe_H_m=Xfe_H_m,
                Xfe_C_m=Xfe_C_m,
                Xfe_N_m=Xfe_N_m,
                Xfe_S_m=Xfe_S_m,
                phase_tracking_cfg=cfg.phase_tracking,
                Xmin_troilite_m=Xmin_troilite_m,
                Xmin_schreibersite_m=Xmin_schreibersite_m,
                Xmin_cohenite_m=Xmin_cohenite_m,
                Xmin_graphite_m=Xmin_graphite_m,
                Xmin_nitride_m=Xmin_nitride_m,
                Xmin_metal_matrix_m=Xmin_metal_matrix_m,
                hydrothermal_active=cfg.hydrothermal.active,
                hydrothermal_cfg=cfg.hydrothermal,
                deltaIW_m=deltaIW_m,
                xm=xm,
                ym=ym,
                coords=coords,
                rplanet_val=acc.rplanet,
                M_planet_val=acc.M_planet_val,
                rcore_val=acc.rcore,
                qxD_val=grids.qxD,
                qyD_val=grids.qyD,
                xcenter_val=acc.xcenter,
                ycenter_val=acc.ycenter,
                kphim0_val=cfg.materials.kphim0,
                phim0_val=cfg.thermodynamics.phim0,
                etafluidm_val=cfg.materials.etafluidm,
                rhofluidm_val=cfg.materials.rhofluidm,
                rhosolidm_val=cfg.materials.rhosolidm,
                etasolidm_val=cfg.materials.etasolidm,
                etasolidmm_val=cfg.materials.etasolidmm,
                rhocpsolidm_val=cfg.materials.rhocpsolidm,
                rhocpfluidm_val=cfg.materials.rhocpfluidm,
                alphafluidm_val=cfg.materials.alphafluidm,
                ksolidm_val=cfg.materials.ksolidm,
                kfluidm_val=cfg.materials.kfluidm,
                phimax_val=cfg.poroelasticity.phimax,
                phi_dry=phi_dry,
            )
            scatter_marker_to_master_grids!(
                m,
                xm[m],
                ym[m],
                coords,
                etatotalm,
                etavpm,
                inv_gggtotalm,
                sxym,
                cohestotalm,
                tenstotalm,
                fricttotalm,
                ETA0SUM,
                ETASUM,
                GGGSUM,
                SXYSUM,
                COHSUM,
                TENSUM,
                FRISUM,
                WTSUM,
                rhototalm,
                rhofluidcur,
                ktotalm,
                phim,
                etafluidcur_inv_kphim,
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
                sxxm,
                rhocptotalm,
                alphasolidcur,
                alphafluidcur,
                hrtotalm,
                tkm_rhocptotalm,
                GGGPSUM,
                SXXSUM,
                RHOSUM,
                RHOCPSUM,
                ALPHASUM,
                ALPHAFSUM,
                HRSUM,
                PHISUM,
                TKSUM,
                WTPSUM;
                s33m=s33m,
                S33SUM=S33SUM,
            )
        end
    end

    if cfg.volatiles.active && !cfg.magma_degassing.active
        copyto!(core.phinewm, core.phim)
    end

    # Normalize basic node properties
    compute_basic_node_properties!(
        ETA0SUM,
        ETASUM,
        GGGSUM,
        SXYSUM,
        COHSUM,
        TENSUM,
        FRISUM,
        WTSUM,
        grids.ETA0,
        grids.ETA,
        grids.GGG,
        grids.SXY0,
        grids.COH,
        grids.TEN,
        grids.FRI,
        grids.YNY,
    )

    # Normalize Vx properties
    compute_vx_node_properties!(
        RHOXSUM,
        RHOFXSUM,
        KXSUM,
        PHIXSUM,
        RXSUM,
        WTXSUM,
        grids.RHOX,
        grids.RHOFX,
        grids.KX,
        grids.PHIX,
        grids.RX,
    )

    # Normalize Vy properties
    compute_vy_node_properties!(
        RHOYSUM,
        RHOFYSUM,
        KYSUM,
        PHIYSUM,
        RYSUM,
        WTYSUM,
        grids.RHOY,
        grids.RHOFY,
        grids.KY,
        grids.PHIY,
        grids.RY,
    )

    # Normalize P-node properties
    compute_p_node_properties!(
        RHOSUM,
        RHOCPSUM,
        ALPHASUM,
        ALPHAFSUM,
        HRSUM,
        GGGPSUM,
        SXXSUM,
        TKSUM,
        PHISUM,
        WTPSUM,
        grids.RHO,
        grids.RHOCP,
        grids.ALPHA,
        grids.ALPHAF,
        grids.HR,
        grids.GGGP,
        grids.SXX0,
        grids.tk1,
        grids.PHI,
        grids.BETAPHI;
        S33SUM=cfg.geometry.spherical_flow ? S33SUM : nothing,
        S330=cfg.geometry.spherical_flow ? grids.S330 : nothing,
    )

    # Boundary conditions: the box edges copy tk1 after the fixed-surface override
    apply_fixed_surface_temperature!(state, coords, cfg; rock_sums=rock_sums)
    apply_insulating_boundary_conditions!(grids.tk1)
    grids.tk2 .= grids.tk1

    return nothing
end

"""
Allocate the 10 P grid arrays of `scatter_rock_thermal_sums!` and the 2 arrays of the
face conductivities before the cut-cell scaling.

$(SIGNATURES)

# Arguments
- `coords`: grid coordinates

# Returns
- `NamedTuple` of zero `Ny1` x `Nx1` arrays `WP`, `WPr`, `C`, `E`, `HRr`, `KPr`, `KXr`,
  `WXr`, `KYr`, `WYr`, `KXu`, `KYu`; `apply_fixed_surface_temperature!` sets `KXu` and
  `KYu` [W/(m K)] for the spherical metric term of the thermal solve

# Throws
- nothing
"""
function allocate_rock_thermal_sums(coords::GridCoordinates)
    sz = (coords.Ny1, coords.Nx1)
    return (;
        WP=zeros(sz),
        WPr=zeros(sz),
        C=zeros(sz),
        E=zeros(sz),
        HRr=zeros(sz),
        KPr=zeros(sz),
        KXr=zeros(sz),
        WXr=zeros(sz),
        KYr=zeros(sz),
        WYr=zeros(sz),
        KXu=zeros(sz),
        KYu=zeros(sz),
    )
end

"""
Sum the thermal properties of the rock markers (`tm < 3`) and the weights of all markers
at the P, Vx and Vy nodes.

$(SIGNATURES)

The weights are the bilinear weights of `marker_to_p_nodes!`, `marker_to_vx_nodes!` and
`marker_to_vy_nodes!`, so the sums do not depend on the scatter path.

# Arguments
- `sums`: result of `allocate_rock_thermal_sums`, reset and filled in place: `WP` (all
  weights), `WPr` (rock weights), `C` (rock heat capacity sum [J/(m^3 K)]), `E` (rock
  `rhocp T` sum [J/m^3]), `HRr` (rock heating sum [W/m^3]), `KPr` (rock conductivity sum
  at P nodes), `KXr`, `WXr`, `KYr`, `WYr` (rock conductivity and weight sums at Vx and Vy
  nodes); the face conductivities `KXu` and `KYu` are reset to 0
- `core`: marker core group (`xm`, `ym`, `tm`, `rhocptotalm`, `tkm_rhocptotalm`,
  `hrtotalm`, `ktotalm` of the step)
- `marknum`: number of markers in use
- `coords`: grid coordinates

# Returns
- `sums`

# Throws
- `ArgumentError`: if an array of `sums` does not have the size of the P grid
- `DomainError`: if a marker position is not finite
"""
function scatter_rock_thermal_sums!(sums, core, marknum::Int, coords::GridCoordinates)
    all(a -> size(a) == (coords.Ny1, coords.Nx1), values(sums)) ||
        throw(ArgumentError("the rock sums must have the size of the P grid"))
    foreach(a -> fill!(a, 0.0), values(sums))
    (; WP, WPr, C, E, HRr, KPr, KXr, WXr, KYr, WYr) = sums
    @unpack_coords coords xp yp xvx yvx xvy yvy dx dy
    @unpack_coords coords jmin_p jmax_p imin_p imax_p jmin_vx jmax_vx imin_vx imax_vx
    @unpack_coords coords jmin_vy jmax_vy imin_vy imax_vy
    @inbounds for m in 1:marknum
        x = core.xm[m]
        y = core.ym[m]
        i, j, w = fix_weights(
            x,
            y,
            xp_val,
            yp_val,
            dx_val,
            dy_val,
            jmin_p_val,
            jmax_p_val,
            imin_p_val,
            imax_p_val,
        )
        interpolate_add_to_grid!(i, j, w, 1.0, WP)
        core.tm[m] < 3 || continue
        k_m = core.ktotalm[m]
        interpolate_add_to_grid!(i, j, w, 1.0, WPr)
        interpolate_add_to_grid!(i, j, w, core.rhocptotalm[m], C)
        interpolate_add_to_grid!(i, j, w, core.tkm_rhocptotalm[m], E)
        interpolate_add_to_grid!(i, j, w, core.hrtotalm[m], HRr)
        interpolate_add_to_grid!(i, j, w, k_m, KPr)
        i, j, w = fix_weights(
            x,
            y,
            xvx_val,
            yvx_val,
            dx_val,
            dy_val,
            jmin_vx_val,
            jmax_vx_val,
            imin_vx_val,
            imax_vx_val,
        )
        interpolate_add_to_grid!(i, j, w, k_m, KXr)
        interpolate_add_to_grid!(i, j, w, 1.0, WXr)
        i, j, w = fix_weights(
            x,
            y,
            xvy_val,
            yvy_val,
            dx_val,
            dy_val,
            jmin_vy_val,
            jmax_vy_val,
            imin_vy_val,
            imax_vy_val,
        )
        interpolate_add_to_grid!(i, j, w, k_m, KYr)
        interpolate_add_to_grid!(i, j, w, 1.0, WYr)
    end
    return sums
end

"""
Give the free P nodes and the faces between them the thermal properties of their rock
markers.

$(SIGNATURES)

At a free node with rock, `RHOCP = C / WP` (rock heat capacity over the weight of all
markers, the heat that the rock markers take), `tk1 = E / C` and `HR = HRr / WP`. A free
node without rock becomes fixed when it lies within 1.5 cells of `rplanet`. A face between
2 free nodes, or between a free node and a fixed one, takes the rock conductivity of the
face, or the rock conductivity of the free node when the face has no rock. Faces between
2 fixed nodes do not change.

# Arguments
- `RHOCP`, `tk1`, `HR`: P node heat capacity [J/(m^3 K)], temperature [K] and heating
  [W/m^3], changed in place
- `KX`, `KY`: conductivity at Vx and Vy nodes [W/(m K)], changed in place
- `fixed`: Bool mask of the P grid, changed in place for nodes without rock
- `sums`: rock sums filled by `scatter_rock_thermal_sums!`
- `coords`: grid coordinates
- `xcenter`, `ycenter`, `rplanet`: planet centre and radius [m]

# Returns
- `Int`: the number of free nodes without rock that became fixed

# Throws
- `ArgumentError`: if a free node deeper than 1.5 cells below `rplanet` has no rock marker
"""
function apply_rock_only_thermal_properties!(
    RHOCP, tk1, HR, KX, KY, fixed, sums, coords::GridCoordinates; xcenter, ycenter, rplanet
)
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    r_shallow = rplanet - 1.5 * max(coords.dx, coords.dy)
    n_norock = 0
    @inbounds for j in 2:(Nx1 - 1), i in 2:(Ny1 - 1)
        fixed[i, j] && continue
        if sums.WPr[i, j] > 0.0 && sums.C[i, j] > 0.0
            RHOCP[i, j] = sums.C[i, j] / sums.WP[i, j]
            tk1[i, j] = sums.E[i, j] / sums.C[i, j]
            HR[i, j] = sums.HRr[i, j] / sums.WP[i, j]
        else
            r = hypot(coords.xp[j] - xcenter, coords.yp[i] - ycenter)
            r > r_shallow || throw(
                ArgumentError(
                    "P node ($i, $j) at r = $r m inside the planet has no rock marker"
                ),
            )
            fixed[i, j] = true
            n_norock += 1
        end
    end
    k_node(i, j) = sums.KPr[i, j] / sums.WPr[i, j]
    @inbounds for j in 2:(Nx1 - 2), i in 2:(Ny1 - 1)
        fixed[i, j] && fixed[i, j + 1] && continue
        if sums.WXr[i, j] > 0.0
            KX[i, j] = sums.KXr[i, j] / sums.WXr[i, j]
        else
            KX[i, j] = fixed[i, j] ? k_node(i, j + 1) : k_node(i, j)
        end
    end
    @inbounds for j in 2:(Nx1 - 1), i in 2:(Ny1 - 2)
        fixed[i, j] && fixed[i + 1, j] && continue
        if sums.WYr[i, j] > 0.0
            KY[i, j] = sums.KYr[i, j] / sums.WYr[i, j]
        else
            KY[i, j] = fixed[i, j] ? k_node(i + 1, j) : k_node(i, j)
        end
    end
    return n_norock
end

"""
Hold the P nodes outside the planet at the disk temperature
(`thermodynamics.surface_temperature = :fixed`).

$(SIGNATURES)

Builds `state.grids.fixed_tk` from the planet centre and radius of `state.accumulators`,
gives the free nodes and faces their rock-only properties
(`apply_rock_only_thermal_properties!`), keeps a copy of these face conductivities in
`KXu` and `KYu` of the sums, places the fixed value at the planet radius
(`apply_cut_cell_conductance!`) and sets `tk1` of the fixed nodes to the disk temperature
of the step. The counts go to `state.accumulators.fixed_surface_counts`. The first pass of
a run that makes free nodes without rock fixed logs a warning. Does nothing with
`:sticky_air`.

# Arguments
- `state`: simulation state after the scatter
- `coords`: grid coordinates
- `cfg`: simulation configuration
- `rock_sums`: arrays from `allocate_rock_thermal_sums`, reused in place, or `nothing` to
  allocate them

# Returns
- `nothing` with `:sticky_air`; else a `NamedTuple` with `n_fixed` (fixed P nodes),
  `n_norock` (free nodes without rock made fixed), `n_rock_out` (rock markers outside
  `rplanet`) and `n_norock_steps` (passes of the run with `n_norock > 0`)

# Throws
- `ArgumentError`: if the state has no `fixed_tk` array, or from
  `scatter_rock_thermal_sums!` and `apply_rock_only_thermal_properties!`
"""
function apply_fixed_surface_temperature!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig;
    rock_sums=nothing,
)
    cfg.thermodynamics.surface_temperature === :fixed || return nothing
    g = state.grids
    fixed = g.fixed_tk
    fixed === nothing && throw(
        ArgumentError("surface_temperature = :fixed but this state has no fixed_tk array"),
    )
    acc = state.accumulators
    geo = (; xcenter=acc.xcenter, ycenter=acc.ycenter, rplanet=acc.rplanet)
    T_fixed, _, _ = compute_ambient_conditions(state.timesum, cfg.disk)
    compute_fixed_surface_mask!(fixed, coords; geo...)
    marknum = length(state.markers)
    core = state.markers.core
    sums = rock_sums === nothing ? allocate_rock_thermal_sums(coords) : rock_sums
    scatter_rock_thermal_sums!(sums, core, marknum, coords)
    n_norock = apply_rock_only_thermal_properties!(
        g.RHOCP, g.tk1, g.HR, g.KX, g.KY, fixed, sums, coords; geo...
    )
    # The explicit metric term of the solve takes the faces before the cut-cell scaling
    copyto!(sums.KXu, g.KX)
    copyto!(sums.KYu, g.KY)
    apply_cut_cell_conductance!(g.KX, g.KY, fixed, coords; geo...)
    n_fixed = 0
    @inbounds for k in eachindex(fixed)
        fixed[k] || continue
        g.tk1[k] = T_fixed
        n_fixed += 1
    end
    n_rock_out = count(
        m ->
            core.tm[m] < 3 &&
            hypot(core.xm[m] - geo.xcenter, core.ym[m] - geo.ycenter) > geo.rplanet,
        1:marknum,
    )
    n_steps_before = acc.fixed_surface_counts.n_norock_steps
    if n_norock > 0 && n_steps_before == 0
        @warn "fixed surface temperature: $n_norock free P nodes within 1.5 cells of " *
            "rplanet have no rock marker and are held at the disk temperature"
    end
    counts = (;
        n_fixed, n_norock, n_rock_out, n_norock_steps=n_steps_before + (n_norock > 0)
    )
    acc.fixed_surface_counts = counts
    @debug "fixed surface temperature" T_fixed counts
    return counts
end

"""
Evaluate constitutive properties and scatter marker contributions for a single P2M tile `t`.

$(SIGNATURES)

# Parameters
- `t`: 1-based tile index in `ws.tile_offsets`
- `ws`: `P2MTiledWorkspace` holding CSR marker permutation and offsets
- `state`: `SimulationState` holding markers and master grids
- `coords`: `GridCoordinates`
- `cfg`: `SimulationConfig`
- `interp_arrays`: 34-tuple of master grid accumulation arrays
- `hr_rates`: tuple of precomputed radiogenic heating rates [W/kg]

# Returns
- `Nothing`
"""
function _scatter_tiled_p2m_tile!(
    t::Int,
    ws::P2MTiledWorkspace,
    state::SimulationState{G},
    coords::GridCoordinates,
    cfg::SimulationConfig,
    interp_arrays::NTuple{34,Matrix{Float64}},
    hr_rates::HR,
)::Nothing where {G,HR<:Tuple}
    lo = @inbounds ws.tile_offsets[t]
    hi = @inbounds ws.tile_offsets[t + 1] - 1
    lo > hi && return nothing

    markers = state.markers
    core = markers.core
    grids = state.grids
    acc = state.accumulators

    tm = core.tm
    xm = core.xm
    ym = core.ym
    tkm = core.tkm
    phim = core.phim
    pfm0 = core.pfm0
    XWsolidm0 = core.XWsolidm0
    rhototalm = core.rhototalm
    rhocptotalm = core.rhocptotalm
    etatotalm = core.etatotalm
    hrtotalm = core.hrtotalm
    ktotalm = core.ktotalm
    tkm_rhocptotalm = core.tkm_rhocptotalm
    etafluidcur_inv_kphim = core.etafluidcur_inv_kphim
    rhofluidcur = core.rhofluidcur
    alphasolidcur = core.alphasolidcur
    alphafluidcur = core.alphafluidcur
    etavpm = core.etavpm
    inv_gggtotalm = core.inv_gggtotalm
    sxym = core.sxym
    sxxm = core.sxxm
    s33m = cfg.geometry.spherical_flow ? core.s33m : nothing
    cohestotalm = core.cohestotalm
    tenstotalm = core.tenstotalm
    fricttotalm = core.fricttotalm
    Fm = core.Fm

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
        S33SUM,
    ) = interp_arrays

    (hrsolidm, hrfluidm, hrmetalm) = hr_rates

    has_metal = haskey(markers.groups, :metal)
    Xfe_bulk = has_metal ? markers.groups.metal.Xfe_bulk : nothing
    Xfem = has_metal ? markers.groups.metal.Xfem : nothing
    Xfe_H_m = has_metal ? markers.groups.metal.Xfe_H_m : nothing
    Xfe_C_m = has_metal ? markers.groups.metal.Xfe_C_m : nothing
    Xfe_N_m = has_metal ? markers.groups.metal.Xfe_N_m : nothing
    Xfe_S_m = has_metal ? markers.groups.metal.Xfe_S_m : nothing

    has_vols = haskey(markers.groups, :volatiles)
    XH2Om = has_vols ? markers.groups.volatiles.XH2Om : nothing
    XCm = has_vols ? markers.groups.volatiles.XCm : nothing
    XNm = has_vols ? markers.groups.volatiles.XNm : nothing
    XSm = has_vols ? markers.groups.volatiles.XSm : nothing
    F_extract_m = has_vols ? markers.groups.volatiles.F_extract_m : nothing

    has_phase = haskey(markers.groups, :phase)
    Xmin_troilite_m = has_phase ? markers.groups.phase.Xmin_troilite_m : nothing
    Xmin_schreibersite_m = has_phase ? markers.groups.phase.Xmin_schreibersite_m : nothing
    Xmin_cohenite_m = has_phase ? markers.groups.phase.Xmin_cohenite_m : nothing
    Xmin_graphite_m = has_phase ? markers.groups.phase.Xmin_graphite_m : nothing
    Xmin_nitride_m = has_phase ? markers.groups.phase.Xmin_nitride_m : nothing
    Xmin_metal_matrix_m = has_phase ? markers.groups.phase.Xmin_metal_matrix_m : nothing

    has_redox = haskey(markers.groups, :redox)
    deltaIW_m = has_redox ? markers.groups.redox.deltaIW_m : nothing
    phi_dry = dry_porosity(markers)

    thermal_buoyancy_val = cfg.thermodynamics.thermal_buoyancy
    alphafluid_val = cfg.materials.alphafluidm
    tmfluidphase_val = cfg.thermodynamics.tmfluidphase
    fluid_viscosity_mode_val = cfg.thermodynamics.fluid_viscosity_mode
    fluid_viscosity_Ea_val = cfg.thermodynamics.fluid_viscosity_Ea
    fluid_viscosity_T0_val = cfg.thermodynamics.fluid_viscosity_T0
    fluid_viscosity_eta0_val = cfg.thermodynamics.fluid_viscosity_eta0
    coreformation_active_val =
        cfg.coreformation.percolation_active || cfg.coreformation.settling_active

    for idx in lo:hi
        m = @inbounds ws.tile_markers[idx]
        compute_marker_properties!(
            m,
            tm,
            tkm,
            rhototalm,
            rhocptotalm,
            etatotalm,
            hrtotalm,
            ktotalm,
            tkm_rhocptotalm,
            etafluidcur_inv_kphim,
            hrsolidm,
            hrfluidm,
            phim,
            XWsolidm0,
            marker_property_mode,
            rhofluidcur;
            compute_hr=false,
            thermal_buoyancy=thermal_buoyancy_val,
            alphafluid=alphafluid_val,
            tmfluidphase_val=tmfluidphase_val,
            fluid_viscosity_mode=fluid_viscosity_mode_val,
            fluid_viscosity_Ea=fluid_viscosity_Ea_val,
            fluid_viscosity_T0=fluid_viscosity_T0_val,
            fluid_viscosity_eta0=fluid_viscosity_eta0_val,
            pm=pfm0,
            Fm=Fm,
            melting_active=cfg.melting.active,
            magma_transport_active=cfg.magma_transport.active,
            track_depletion=cfg.magma_transport.track_depletion,
            F_extract_m=F_extract_m,
            T_solidus_val=cfg.melting.T_solidus,
            T_liquidus_val=cfg.melting.T_liquidus,
            L_melt_val=cfg.melting.L_melt,
            rho_melt_val=cfg.melting.rho_melt,
            alpha_eta_val=cfg.melting.alpha_eta,
            phi_crit_val=cfg.melting.phi_crit,
            eta_melt_val=cfg.melting.eta_melt,
            dpdt_clapeyron_val=cfg.melting.dpdt_clapeyron,
            soft_turbulence=cfg.melting.soft_turbulence,
            eta_fluid_silicate_val=cfg.melting.eta_fluid_silicate,
            F_turb_start_val=cfg.melting.F_turb_start,
            F_turb_end_val=cfg.melting.F_turb_end,
            turb_exponent_val=cfg.melting.turb_exponent,
            dT_turb_min_val=cfg.melting.dT_turb_min,
            T_surface_ref_val=cfg.melting.T_surface_ref,
            k_turb_cutoff_val=cfg.melting.k_turb_cutoff,
            k_turb_floor_val=cfg.melting.k_turb_floor,
            tmsolidphase_val=cfg.thermodynamics.tmsolidphase,
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            coreformation_active=coreformation_active_val,
            hrmetalm=hrmetalm,
            sulfur_fraction_val=cfg.coreformation.sulfur_fraction,
            metal_density_mode_val=cfg.coreformation.metal_density_mode,
            T_eutectic_val=cfg.coreformation.T_eutectic,
            dT_metal_val=cfg.coreformation.dT_metal,
            rho_metal_val=cfg.coreformation.rho_metal,
            rho_metal_solid_val=cfg.coreformation.rho_metal_solid,
            L_metal_val=cfg.coreformation.L_metal,
            k_metal_val=cfg.coreformation.k_metal,
            rhocp_metal_val=cfg.coreformation.rhocp_metal,
            volatiles_active=cfg.volatiles.active,
            magma_degassing_active=cfg.magma_degassing.active,
            etamin=cfg.solver.etamin,
            etamax=cfg.solver.etamax,
            volatiles_cfg=cfg.volatiles,
            retention_cfg=cfg.retention,
            XH2Om=XH2Om,
            XCm=XCm,
            XNm=XNm,
            XSm=XSm,
            metal_partition_cfg=cfg.metal_partition,
            Xfe_H_m=Xfe_H_m,
            Xfe_C_m=Xfe_C_m,
            Xfe_N_m=Xfe_N_m,
            Xfe_S_m=Xfe_S_m,
            phase_tracking_cfg=cfg.phase_tracking,
            Xmin_troilite_m=Xmin_troilite_m,
            Xmin_schreibersite_m=Xmin_schreibersite_m,
            Xmin_cohenite_m=Xmin_cohenite_m,
            Xmin_graphite_m=Xmin_graphite_m,
            Xmin_nitride_m=Xmin_nitride_m,
            Xmin_metal_matrix_m=Xmin_metal_matrix_m,
            hydrothermal_active=cfg.hydrothermal.active,
            hydrothermal_cfg=cfg.hydrothermal,
            deltaIW_m=deltaIW_m,
            xm=xm,
            ym=ym,
            coords=coords,
            rplanet_val=acc.rplanet,
            M_planet_val=acc.M_planet_val,
            rcore_val=acc.rcore,
            qxD_val=grids.qxD,
            qyD_val=grids.qyD,
            xcenter_val=acc.xcenter,
            ycenter_val=acc.ycenter,
            kphim0_val=cfg.materials.kphim0,
            phim0_val=cfg.thermodynamics.phim0,
            etafluidm_val=cfg.materials.etafluidm,
            rhofluidm_val=cfg.materials.rhofluidm,
            rhosolidm_val=cfg.materials.rhosolidm,
            etasolidm_val=cfg.materials.etasolidm,
            etasolidmm_val=cfg.materials.etasolidmm,
            rhocpsolidm_val=cfg.materials.rhocpsolidm,
            rhocpfluidm_val=cfg.materials.rhocpfluidm,
            alphafluidm_val=cfg.materials.alphafluidm,
            ksolidm_val=cfg.materials.ksolidm,
            kfluidm_val=cfg.materials.kfluidm,
            phimax_val=cfg.poroelasticity.phimax,
            phi_dry=phi_dry,
        )
        scatter_marker_to_master_grids!(
            m,
            xm[m],
            ym[m],
            coords,
            etatotalm,
            etavpm,
            inv_gggtotalm,
            sxym,
            cohestotalm,
            tenstotalm,
            fricttotalm,
            ETA0SUM,
            ETASUM,
            GGGSUM,
            SXYSUM,
            COHSUM,
            TENSUM,
            FRISUM,
            WTSUM,
            rhototalm,
            rhofluidcur,
            ktotalm,
            phim,
            etafluidcur_inv_kphim,
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
            sxxm,
            rhocptotalm,
            alphasolidcur,
            alphafluidcur,
            hrtotalm,
            tkm_rhocptotalm,
            GGGPSUM,
            SXXSUM,
            RHOSUM,
            RHOCPSUM,
            ALPHASUM,
            ALPHAFSUM,
            HRSUM,
            PHISUM,
            TKSUM,
            WTPSUM;
            s33m=s33m,
            S33SUM=S33SUM,
        )
    end
    return nothing
end

# Marker property computation and particle-to-mesh (P2M) grid interpolation

"""
    interpolate_markers_to_grid!(state::SimulationState, coords::GridCoordinates, cfg::SimulationConfig; p2m_workspace=nothing, thread_buffers=nothing, interp_arrays=nothing)::Nothing

Interpolate marker properties to staggered Eulerian mesh and normalize nodal properties.

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
            WTPSUM,
        )
    end

    # Optional marker groups
    has_metal = haskey(markers.groups, :metal)
    Xfe_bulk = has_metal ? markers.groups.metal.Xfe_bulk : nothing
    Xfem = has_metal ? markers.groups.metal.Xfem : nothing
    has_metal_part = has_metal && cfg.metal_partition.active
    Xfe_H_m = has_metal_part ? markers.groups.metal.Xfe_H_m : nothing
    Xfe_C_m = has_metal_part ? markers.groups.metal.Xfe_C_m : nothing
    Xfe_N_m = has_metal_part ? markers.groups.metal.Xfe_N_m : nothing
    Xfe_S_m = has_metal_part ? markers.groups.metal.Xfe_S_m : nothing

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
        update_marker_redox!(
            redox_props,
            tkm,
            pfm0,
            cfg.redox;
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            XWsolidm=core.XWsolidm,
            rhosolid=cfg.materials.rhosolidm[1],
            rho_metal=cfg.coreformation.rho_metal,
        )
    end

    tau_al = cfg.thermodynamics.t_half_al / log(2.0)
    tau_fe = cfg.thermodynamics.t_half_fe / log(2.0)

    hrsolidm, hrfluidm, hrmetalm = calculate_radioactive_heating(
        cfg.thermodynamics.hr_al,
        cfg.thermodynamics.hr_fe,
        state.timesum;
        ratio_al=cfg.thermodynamics.ratio_al,
        E_al=cfg.thermodynamics.E_al,
        f_al=cfg.thermodynamics.f_al,
        tau_al=tau_al,
        ratio_fe=cfg.thermodynamics.ratio_fe,
        E_fe=cfg.thermodynamics.E_fe,
        f_fe=cfg.thermodynamics.f_fe,
        tau_fe=tau_fe,
        rho_metal=cfg.coreformation.rho_metal,
        rhosolidm=cfg.materials.rhosolidm,
        rhofluidm=cfg.materials.rhofluidm,
    )

    thermal_buoyancy_val = cfg.thermodynamics.thermal_buoyancy
    alphafluid_val = cfg.materials.alphafluidm
    tmfluidphase_val = cfg.thermodynamics.tmfluidphase
    fluid_viscosity_mode_val = cfg.thermodynamics.fluid_viscosity_mode
    fluid_viscosity_Ea_val = cfg.thermodynamics.fluid_viscosity_Ea
    fluid_viscosity_T0_val = cfg.thermodynamics.fluid_viscosity_T0
    fluid_viscosity_eta0_val = cfg.thermodynamics.fluid_viscosity_eta0
    coreformation_active_val =
        cfg.coreformation.percolation_active || cfg.coreformation.settling_active

    use_tiled_p2m = cfg.solver.p2m_mode === :tiled
    use_threading =
        cfg.solver.p2m_mode === :buffered &&
        Threads.nthreads() > 1 &&
        thread_buffers !== nothing

    if use_tiled_p2m
        ws = if p2m_workspace !== nothing
            ensure_workspace_compatible(p2m_workspace, coords, marknum, cfg.solver.tile_size)
        else
            P2MTiledWorkspace(coords, marknum, cfg.solver.tile_size)
        end
        bin_markers_into_tiles!(ws, xm, ym, coords, marknum)
        for color in 1:4
            tiles = ws.tiles_by_color[color]
            Threads.@threads :dynamic for t in tiles
                lo = ws.tile_offsets[t]
                hi = ws.tile_offsets[t + 1] - 1
                lo > hi && continue
                for idx in lo:hi
                    m = ws.tile_markers[idx]
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
                        WTPSUM,
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
            thread_buffers,
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
                WTPSUM,
            )
        end
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
        grids.BETAPHI,
    )

    # Boundary conditions
    apply_insulating_boundary_conditions!(grids.tk1)
    grids.tk2 .= grids.tk1

    return nothing
end

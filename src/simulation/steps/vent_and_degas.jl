# Hydrofracture venting, volatile drainage, and magma ocean degassing step

"""
    vent_and_degas!(state::SimulationState, coords::GridCoordinates, cfg::SimulationConfig; Fm_step_start=nothing)

Advance marker thermal state, porosity compaction, hydrofracture venting, and magma ocean volatile degassing.

# Mutates:
- `state.markers.core.XWsolidm0`
- `state.markers.core.phim`
- `state.markers.core.tkm`
- `state.markers.core.tm`
- `state.markers.core.w3d_m`
- `state.markers.volatiles` (if active)
- `state.markers.metal` (if active)
- `state.accumulators.M_vent_total`
- `state.accumulators.M_vent_H2O_total`
- `state.accumulators.M_vent_C_total`
- `state.accumulators.M_vent_N_total`
- `state.accumulators.M_vent_S_total`
- `state.transfers`
- `state.grids.S_vent_grid`

# Returns:
- `NamedTuple`: `(; delta_m_vent_3d, vented_vols, degas_rates)`
"""
function vent_and_degas!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig;
    Fm_step_start=nothing,
)
    marknum = length(state.markers)
    if marknum == 0
        return (; delta_m_vent_3d=0.0, vented_vols=nothing, degas_rates=nothing)
    end

    markers = state.markers
    grids = state.grids
    acc = state.accumulators
    transfer_log = state.transfers
    atm_state = state.atm

    xm = markers.core.xm
    ym = markers.core.ym
    tm = markers.core.tm
    tkm = markers.core.tkm
    phim = markers.core.phim
    phinewm = markers.core.phinewm
    XWsolidm = markers.core.XWsolidm
    XWsolidm0 = markers.core.XWsolidm0
    w3d_m = markers.core.w3d_m

    has_metal = haskey(markers.groups, :metal)
    Xfem = has_metal ? markers.groups.metal.Xfem : nothing
    Xfem0 = has_metal ? markers.groups.metal.Xfem0 : nothing
    Xfe_bulk = has_metal ? markers.groups.metal.Xfe_bulk : nothing
    Fm = markers.core.Fm
    has_redox = haskey(markers.groups, :redox)
    redox_props = has_redox ? markers.groups.redox : nothing

    has_vols = haskey(markers.groups, :volatiles)
    XH2Om = has_vols ? markers.groups.volatiles.XH2Om : nothing
    XCm = has_vols ? markers.groups.volatiles.XCm : nothing
    XNm = has_vols ? markers.groups.volatiles.XNm : nothing
    XSm = has_vols ? markers.groups.volatiles.XSm : nothing

    dt = state.dt
    timestep = state.timestep
    timesum = state.timesum
    rplanet_val = acc.rplanet
    xcenter_val = acc.xcenter
    ycenter_val = acc.ycenter
    M_planet_val = acc.M_planet_val
    P_amb_eff = acc.P_amb

    T_amb, _, _ = compute_ambient_conditions(timesum, cfg.disk)

    phimin_val = cfg.poroelasticity.phimin
    phimax_val = cfg.poroelasticity.phimax
    rhofluidm = cfg.materials.rhofluidm

    DT = grids.DT
    tk1 = grids.tk1
    tk2 = grids.tk2
    APHI = grids.APHI
    S_vent_grid = grids.S_vent_grid

    # Advance thermal state and porosity compaction
    XWsolidm0 .= XWsolidm
    phim .= phinewm
    if (
            cfg.coreformation.percolation_active ||
            cfg.coreformation.settling_active ||
            cfg.metal_partition.active
        ) &&
        Xfem !== nothing &&
        Xfem0 !== nothing
        Xfem0 .= Xfem
    end

    if length(w3d_m) != marknum
        resize!(markers.core.w3d_m, marknum)
        for m in 1:marknum
            markers.core.w3d_m[m] = marker_out_of_plane_length(
                xm[m], ym[m], xcenter_val, ycenter_val
            )
        end
        w3d_m = markers.core.w3d_m
    end

    marker_results = advance_marker_thermo_porosity_venting!(
        xm,
        ym,
        tm,
        tkm,
        phim,
        DT,
        tk2,
        APHI,
        dt,
        timestep,
        marknum;
        coords=coords,
        phimin=phimin_val,
        phimax=phimax_val,
        venting=cfg.venting.active,
        S_vent_grid=S_vent_grid,
        rhofluidcur=rhofluidm[2],
        ret_cfg=cfg.retention,
        XH2Om=cfg.volatiles.active ? XH2Om : nothing,
        XCm=cfg.volatiles.active ? XCm : nothing,
        XNm=cfg.volatiles.active ? XNm : nothing,
        XSm=cfg.volatiles.active ? XSm : nothing,
        rhosolid=cfg.materials.rhosolidm,
        Fm=Fm,
        redox_props=redox_props,
        w3d_m=w3d_m,
        Xfe_bulk=Xfe_bulk,
    )
    append!(transfer_log, marker_results.records)

    delta_m_vent_3d = marker_results.delta_m_vent_3d
    vented_vols = marker_results.vented_vols

    if cfg.venting.active
        acc.M_vent_total += delta_m_vent_3d

        if vented_vols !== nothing
            acc.M_vent_H2O_total += vented_vols.M_vent_H2O_3d
            acc.M_vent_C_total += vented_vols.M_vent_C_3d
            acc.M_vent_N_total += vented_vols.M_vent_N_3d
            acc.M_vent_S_total += vented_vols.M_vent_S_3d
        end
    end

    # Surface degassing from magma ocean
    degas_rates = if cfg.magma_degassing.active && XH2Om !== nothing && Fm !== nothing
        p_surf_mo =
            if (cfg.atmosphere.active && atm_state !== nothing && atm_state.P_surf > 0.0)
                atm_state.P_surf
            else
                P_amb_eff
            end
        fO2_diw_mo = compute_surface_mean_delta_iw(
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
        v_m = marker_area(coords)
        Fm_prev = Fm_step_start !== nothing ? Fm_step_start : Fm

        # Mass-weighted mean temperature over degassing-zone markers
        sum_t_mass = 0.0
        sum_degas_mass = 0.0
        r_degas_sq = (cfg.magma_degassing.degas_depth_fraction * rplanet_val)^2
        r_planet_sq = rplanet_val^2
        f_thresh = cfg.magma_degassing.F_melt_threshold
        rho_rock = cfg.materials.rhosolidm[1]

        @inbounds for m in 1:marknum
            if tm[m] >= 3
                continue
            end
            dx = xm[m] - xcenter_val
            dy = ym[m] - ycenter_val
            r_sq = dx * dx + dy * dy
            if r_sq > r_planet_sq
                continue
            end
            f_m = Fm[m]
            if (r_sq >= r_degas_sq) && (f_m >= f_thresh || f_m > 0.01) && (f_m > 0.0)
                w3d = w3d_m !== nothing ? w3d_m[m] : (2.0 * sqrt(r_sq))
                m_wt = rho_rock * v_m * w3d
                sum_t_mass += tkm[m] * m_wt
                sum_degas_mass += m_wt
            end
        end
        T_melt_ref = sum_degas_mass > 0.0 ? (sum_t_mass / sum_degas_mass) : 1500.0

        if cfg.magma_degassing.mode === :dynamic_flux
            degas_res = degas_magma_ocean_markers!(
                xm,
                ym,
                tm,
                tkm,
                Fm,
                Fm_prev,
                XH2Om,
                XCm,
                XNm,
                XSm,
                marknum,
                dt,
                p_surf_mo,
                rplanet_val,
                cfg.magma_degassing,
                T_melt_ref;
                xcenter=xcenter_val,
                ycenter=ycenter_val,
                w3d_m=w3d_m,
                rho_solid=cfg.materials.rhosolidm[1],
                marker_volume=v_m,
                delta_IW=fO2_diw_mo,
                retention_cfg=cfg.retention,
                step=timestep,
                Xfe_bulk=Xfe_bulk,
            )
            append!(transfer_log, degas_res.records)
            if cfg.atmosphere.active
                ElementInventory(
                    degas_res.dM_3D[:H] / dt,
                    degas_res.dM_3D[:C] / dt,
                    degas_res.dM_3D[:N] / dt,
                    degas_res.dM_3D[:S] / dt,
                    degas_res.dM_3D[:H2O] * (15.9994 / 18.01528) / dt,
                )
            else
                degas_res.rates
            end
        elseif cfg.atmosphere.active && atm_state !== nothing
            _degas_magma_ocean_equilibrium!(
                state, coords, cfg, atm_state, fO2_diw_mo, T_amb, v_m, marknum
            )
        else
            nothing
        end
    else
        nothing
    end

    return (;
        delta_m_vent_3d=delta_m_vent_3d, vented_vols=vented_vols, degas_rates=degas_rates
    )
end

"""
    _degas_magma_ocean_equilibrium!(
        state::SimulationState,
        coords::GridCoordinates,
        cfg::SimulationConfig,
        atm_state::AtmosphereState,
        fO2_diw_mo::Float64,
        T_amb::Float64,
        v_m::Float64,
        marknum::Int,
    )

Calculate equilibrium magma ocean volatile partitioning and deplete molten markers.
"""
function _degas_magma_ocean_equilibrium!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig,
    atm_state::AtmosphereState,
    fO2_diw_mo::Float64,
    T_amb::Float64,
    v_m::Float64,
    marknum::Int,
)
    markers = state.markers
    grids = state.grids
    acc = state.accumulators
    transfer_log = state.transfers

    xm = markers.core.xm
    ym = markers.core.ym
    tm = markers.core.tm
    w3d_m = markers.core.w3d_m
    Fm = markers.core.Fm

    has_vols = haskey(markers.groups, :volatiles)
    XH2Om = has_vols ? markers.groups.volatiles.XH2Om : nothing
    XCm = has_vols ? markers.groups.volatiles.XCm : nothing
    XNm = has_vols ? markers.groups.volatiles.XNm : nothing
    XSm = has_vols ? markers.groups.volatiles.XSm : nothing

    (XH2Om === nothing || Fm === nothing) && return nothing

    dt = state.dt
    timestep = state.timestep
    rplanet_val = acc.rplanet
    xcenter_val = acc.xcenter
    ycenter_val = acc.ycenter
    M_planet_val = acc.M_planet_val
    tk1 = grids.tk1

    T_int_val = compute_mean_surface_temperature(
        tk1, coords, rplanet_val, xcenter_val, ycenter_val; T_default=T_amb
    )
    m_melt_tot = 0.0
    m_H_melt = 0.0
    m_C_melt = 0.0
    m_N_melt = 0.0
    m_S_melt = 0.0
    for m in 1:marknum
        if tm[m] < 3 &&
            (((xm[m] - xcenter_val)^2 + (ym[m] - ycenter_val)^2) <= rplanet_val^2) &&
            Fm[m] >= cfg.magma_degassing.F_melt_threshold
            w3d = if w3d_m !== nothing
                w3d_m[m]
            else
                2.0 * hypot(xm[m] - xcenter_val, ym[m] - ycenter_val)
            end
            m_marker_3d = cfg.materials.rhosolidm[1] * v_m * w3d
            m_melt_tot += Fm[m] * m_marker_3d
            m_H_melt += (XH2Om[m] * 0.01) * (2.01588 / 18.01528) * m_marker_3d
            m_C_melt += (XCm[m] * 1.0e-6) * m_marker_3d
            m_N_melt += (XNm[m] * 1.0e-6) * m_marker_3d
            m_S_melt += (XSm[m] * 1.0e-6) * m_marker_3d
        end
    end

    # Atmospheric elemental inventories
    m_H_atm =
        get(atm_state.M_atm, :H2, 0.0) * 1.0 +
        get(atm_state.M_atm, :H2O, 0.0) * (2.01588 / 18.01528) +
        get(atm_state.M_atm, :CH4, 0.0) * (4.03176 / 16.04246) +
        get(atm_state.M_atm, :NH3, 0.0) * (3.02382 / 17.03052) +
        get(atm_state.M_atm, :H2S, 0.0) * (2.01588 / 34.08088)

    m_C_atm =
        get(atm_state.M_atm, :CO, 0.0) * (12.011 / 28.0101) +
        get(atm_state.M_atm, :CO2, 0.0) * (12.011 / 44.0095) +
        get(atm_state.M_atm, :CH4, 0.0) * (12.011 / 16.04246)

    m_N_atm =
        get(atm_state.M_atm, :N2, 0.0) * 1.0 +
        get(atm_state.M_atm, :NH3, 0.0) * (14.007 / 17.03052)

    m_S_atm =
        get(atm_state.M_atm, :H2S, 0.0) * (32.060 / 34.08088) +
        get(atm_state.M_atm, :SO2, 0.0) * (32.060 / 64.066) +
        get(atm_state.M_atm, :S2, 0.0) * 1.0

    m_H_tot = m_H_melt + m_H_atm
    m_C_tot = m_C_melt + m_C_atm
    m_N_tot = m_N_melt + m_N_atm
    m_S_tot = m_S_melt + m_S_atm

    if m_melt_tot > 0.0 && (m_H_tot + m_C_tot + m_N_tot + m_S_tot) > 0.0
        g_surf = GRAVITATIONAL_CONSTANT * M_planet_val / (rplanet_val^2)
        sol_eq = solve_magma_ocean_volatile_partitioning(
            m_melt_tot,
            m_H_tot,
            m_C_tot,
            m_N_tot,
            m_S_tot,
            rplanet_val,
            g_surf,
            T_int_val,
            fO2_diw_mo,
        )

        # Deplete molten markers according to residual melt volatile concentration
        new_XH2O_wtpct = (sol_eq.M_melt_H * (18.01528 / 2.01588) / m_melt_tot) * 100.0
        new_XC_ppm = (sol_eq.M_melt_C / m_melt_tot) * 1.0e6
        new_XN_ppm = (sol_eq.M_melt_N / m_melt_tot) * 1.0e6
        new_XS_ppm = (sol_eq.M_melt_S / m_melt_tot) * 1.0e6

        for m in 1:marknum
            if tm[m] < 3 &&
                (((xm[m] - xcenter_val)^2 + (ym[m] - ycenter_val)^2) <= rplanet_val^2) &&
                Fm[m] >= cfg.magma_degassing.F_melt_threshold
                old_XH2O = XH2Om[m]
                old_XC = XCm[m]
                old_XN = XNm[m]
                old_XS = XSm[m]
                XH2Om[m] = new_XH2O_wtpct
                XCm[m] = new_XC_ppm
                XNm[m] = new_XN_ppm
                XSm[m] = new_XS_ppm
                m_rock_2d = cfg.materials.rhosolidm[1] * v_m
                w3d = if w3d_m !== nothing
                    w3d_m[m]
                else
                    2.0 * hypot(xm[m] - xcenter_val, ym[m] - ycenter_val)
                end
                if old_XH2O != new_XH2O_wtpct
                    dm_h2o_2d = (old_XH2O - new_XH2O_wtpct) * 0.01 * m_rock_2d
                    dm_h_2d = dm_h2o_2d * (2.01588 / 18.01528)
                    push!(
                        transfer_log,
                        TransferRecord(
                            timestep,
                            :degassing,
                            :H,
                            m,
                            xm[m],
                            ym[m],
                            dm_h_2d,
                            dm_h_2d * w3d,
                        ),
                    )
                end
                if old_XC != new_XC_ppm
                    dm_c_2d = (old_XC - new_XC_ppm) * 1.0e-6 * m_rock_2d
                    push!(
                        transfer_log,
                        TransferRecord(
                            timestep,
                            :degassing,
                            :C,
                            m,
                            xm[m],
                            ym[m],
                            dm_c_2d,
                            dm_c_2d * w3d,
                        ),
                    )
                end
                if old_XN != new_XN_ppm
                    dm_n_2d = (old_XN - new_XN_ppm) * 1.0e-6 * m_rock_2d
                    push!(
                        transfer_log,
                        TransferRecord(
                            timestep,
                            :degassing,
                            :N,
                            m,
                            xm[m],
                            ym[m],
                            dm_n_2d,
                            dm_n_2d * w3d,
                        ),
                    )
                end
                if old_XS != new_XS_ppm
                    dm_s_2d = (old_XS - new_XS_ppm) * 1.0e-6 * m_rock_2d
                    push!(
                        transfer_log,
                        TransferRecord(
                            timestep,
                            :degassing,
                            :S,
                            m,
                            xm[m],
                            ym[m],
                            dm_s_2d,
                            dm_s_2d * w3d,
                        ),
                    )
                end
            end
        end

        rates = Dict{Symbol,Float64}()
        for (sp, m_atm_eq) in sol_eq.M_atm_i
            m_atm_cur = get(atm_state.M_atm, sp, 0.0)
            rates[sp] = max(0.0, m_atm_eq - m_atm_cur) / dt
        end
        return rates
    else
        return nothing
    end
end

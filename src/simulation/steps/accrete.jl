# Planetary accretion step: mass addition, impact heating, and boundary advance

"""
    accrete!(state::SimulationState, coords::GridCoordinates, cfg::SimulationConfig)::Nothing

Advance planetesimal accretion boundary, update marker inventory and mass accumulators.

# Mutates:
- `state.accumulators.rplanet`
- `state.accumulators.M_planet_val`
- `state.accumulators.M_accreted_total`
- `state.markers.core.tm`
- `state.markers.core.tkm`
- `state.markers.core.phim`
- `state.markers.core.XWsolidm`
"""
function accrete!(
    state::SimulationState, coords::GridCoordinates, cfg::SimulationConfig
)::Nothing
    # Empty marker set contract: return immediately with no mutations
    length(state.markers) == 0 && return nothing

    if !cfg.accretion.active
        return nothing
    end

    acc = state.accumulators
    markers = state.markers
    timesum = state.timesum
    dt = state.dt
    rplanet_val = acc.rplanet
    M_planet_val = acc.M_planet_val

    dM_dt_acc = compute_accretion_rate(
        timesum, M_planet_val, rplanet_val, cfg.accretion, cfg.disk
    )
    # Clamp mass increment so M_planet_val does not overshoot M_target
    dM_remain = max(0.0, cfg.accretion.M_target - M_planet_val)
    dM_acc = min(dM_dt_acc * dt, dM_remain)

    if dM_acc > 0.0 && rplanet_val < cfg.accretion.R_target
        dR_acc = compute_radius_increment(rplanet_val, dM_acc, cfg.accretion.rho_bulk)
        # Clamp radius increment so rplanet_val does not overshoot R_target
        dR_remain = max(0.0, cfg.accretion.R_target - rplanet_val)
        dR_acc = min(dR_acc, dR_remain)

        T_amb, P_amb, _ = compute_ambient_conditions(timesum, cfg.disk)
        isfinite(T_amb) ||
            throw(DomainError(T_amb, "Ambient disk temperature must be finite, got $T_amb"))

        T_acc = T_amb
        if cfg.accretion.h_impact > 0.0
            _, delta_T_imp = compute_impact_heating(
                M_planet_val,
                rplanet_val;
                h_impact=cfg.accretion.h_impact,
                c_p=cfg.accretion.cp_rock,
                v_inf=cfg.accretion.v_inf,
            )
            T_acc += delta_T_imp
        end

        XW_acc = cfg.accretion.XWsolid_dry
        H2O_acc = cfg.accretion.XH2O_dry_wtpct
        if cfg.accretion.snowline_coupling
            XW_acc, H2O_acc = evaluate_snowline_water_content(
                T_amb;
                T_snowline_cond=cfg.accretion.T_snowline_cond,
                XW_wet=cfg.accretion.XWsolid_wet,
                XW_dry=cfg.accretion.XWsolid_dry,
                H2O_wet_wtpct=cfg.accretion.XH2O_wet_wtpct,
                H2O_dry_wtpct=cfg.accretion.XH2O_dry_wtpct,
            )
        end

        cond_state_acc = if (cfg.volatile_mixture.active || cfg.refractory.active)
            evaluate_disk_volatile_condensation(
                T_amb,
                P_amb,
                cfg.volatile_mixture,
                cfg.refractory;
                P_ref=cfg.volatile_mixture.P_ref,
                alpha_P=cfg.volatile_mixture.alpha_P,
            )
        else
            nothing
        end

        if cfg.volatile_mixture.active && cond_state_acc !== nothing
            if cond_state_acc.condensed_H2O
                XW_acc = cfg.volatile_mixture.X_ice_H2O
                H2O_acc = cfg.volatile_mixture.X_ice_H2O * 100.0
            else
                XW_acc = cfg.accretion.XWsolid_dry
                H2O_acc = cfg.accretion.XH2O_dry_wtpct
            end
        end

        advance_accretion_boundary!(
            rplanet_val,
            dR_acc,
            markers;
            xcenter=acc.xcenter,
            ycenter=acc.ycenter,
            T_accreted=T_acc,
            phi_accreted=cfg.accretion.phi_accreted,
            XWsolid_accreted=XW_acc,
            Xfe_accreted=cfg.accretion.Xfe_bulk_accreted,
            current_time=timesum,
            XH2O_accreted=H2O_acc,
            XC_accreted=cfg.accretion.XC_accreted_ppm,
            XN_accreted=cfg.accretion.XN_accreted_ppm,
            XS_accreted=cfg.accretion.XS_accreted_ppm,
            cfg=cfg,
            disk_state=cond_state_acc,
        )

        acc.rplanet += dR_acc
        acc.M_planet_val += dM_acc
        acc.M_accreted_total += dM_acc
    end

    return nothing
end

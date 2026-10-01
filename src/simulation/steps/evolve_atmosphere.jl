# Coupled surface atmosphere and hydrodynamic escape evolution step

"""
    evolve_atmosphere!(state::SimulationState, coords::GridCoordinates, cfg::SimulationConfig; vent_degas_result=nothing)::Nothing

Evolve coupled surface atmosphere inventory, surface pressure/temperature, and hydrodynamic escape.

# Mutates:
- `state.atm`
- `state.accumulators.M_atm_total`
- `state.accumulators.M_escaped_total`
- `state.accumulators.M_atm_species`
- `state.accumulators.M_escaped_species`
"""
function evolve_atmosphere!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig;
    vent_degas_result=nothing,
)::Nothing
    marknum = length(state.markers)
    if marknum == 0 || (!cfg.atmosphere.active && !cfg.escape.active)
        return nothing
    end

    delta_m_vent_3d =
        vent_degas_result !== nothing ? vent_degas_result.delta_m_vent_3d : 0.0
    vented_vols = vent_degas_result !== nothing ? vent_degas_result.vented_vols : nothing
    degas_rates = vent_degas_result !== nothing ? vent_degas_result.degas_rates : nothing

    atm_state = state.atm
    acc = state.accumulators
    grids = state.grids
    markers = state.markers

    xm = markers.core.xm
    ym = markers.core.ym
    tm = markers.core.tm
    has_redox = haskey(markers.groups, :redox)
    redox_props = has_redox ? markers.groups.redox : nothing

    dt = state.dt
    timesum = state.timesum
    rplanet_val = acc.rplanet
    xcenter_val = acc.xcenter
    ycenter_val = acc.ycenter
    M_planet_val = acc.M_planet_val
    P_amb_eff = acc.P_amb

    tk1 = grids.tk1

    T_amb, _, w_disp = compute_ambient_conditions(timesum, cfg.disk)

    M_atm_species = acc.M_atm_species
    M_escaped_species = acc.M_escaped_species

    if cfg.atmosphere.active && atm_state !== nothing
        p_surf_val = max(atm_state.P_surf, P_amb_eff)
        T_surf_val = atm_state.T_surf_eq > 0.0 ? atm_state.T_surf_eq : T_amb
        vent_rates = compute_surface_venting_rates(
            cfg,
            delta_m_vent_3d,
            vented_vols,
            dt,
            p_surf_val,
            T_surf_val,
            redox_props,
            marknum,
            tm,
            xm,
            ym,
            rplanet_val,
            xcenter_val,
            ycenter_val,
        )

        c_s_disk = compute_sound_speed(T_amb)
        rho_disk_val = if (cfg.disk.enabled && c_s_disk > 0.0)
            max(0.0, (1.0 - w_disp) * cfg.disk.p_amb_disk) / (c_s_disk^2)
        else
            0.0
        end
        a_orb_val = cfg.disk.orbital_distance_au * AU_METERS
        M_star_val = cfg.disk.stellar_mass_msun * M_SUN_KG
        R_exo_val = max(cfg.escape.R_exobase, rplanet_val)
        T_int_val = compute_mean_surface_temperature(
            tk1, coords, rplanet_val, xcenter_val, ycenter_val; T_default=T_amb
        )

        evolve_coupled_atmosphere_step!(
            atm_state,
            vent_rates,
            dt,
            M_planet_val,
            rplanet_val,
            T_amb,
            cfg.atmosphere;
            rho_disk=rho_disk_val,
            c_s=c_s_disk,
            M_star=M_star_val,
            a_orb=a_orb_val,
            T_int=T_int_val,
            T_exobase=cfg.escape.T_exobase,
            R_exobase=R_exo_val,
            hydrodynamic=cfg.escape.hydrodynamic,
            gamma=cfg.escape.gamma,
            escape_active=cfg.escape.active,
            escape_cfg=cfg.escape,
            sim_time_s=timesum,
            degas_rates=degas_rates,
        )

        acc.M_atm_total = sum(values(atm_state.M_atm))
        acc.M_escaped_total = sum(values(atm_state.M_escaped))
        if M_atm_species !== nothing
            for (sp, val) in atm_state.M_atm
                M_atm_species[sp] = val
            end
        end
        if M_escaped_species !== nothing
            for (sp, val) in atm_state.M_escaped
                M_escaped_species[sp] = val
            end
        end
    elseif cfg.escape.active
        R_exo_val = max(cfg.escape.R_exobase, rplanet_val)
        T_surf_esc = compute_mean_surface_temperature(
            tk1, coords, rplanet_val, xcenter_val, ycenter_val; T_default=T_amb
        )

        vent_rates_esc = compute_surface_venting_rates(
            cfg,
            delta_m_vent_3d,
            vented_vols,
            dt,
            P_amb_eff,
            T_surf_esc,
            redox_props,
            marknum,
            tm,
            xm,
            ym,
            rplanet_val,
            xcenter_val,
            ycenter_val,
        )

        if cfg.escape.multi_species &&
            M_atm_species !== nothing &&
            M_escaped_species !== nothing
            for sp in cfg.escape.species_list
                m_sp = get_species_molecular_mass(sp)
                v_rate_sp = get(vent_rates_esc, sp, 0.0)
                prev_sp = get(M_atm_species, sp, 0.0)
                esc_sp = evolve_atmospheric_species_inventory(
                    prev_sp,
                    v_rate_sp,
                    dt,
                    M_planet_val,
                    rplanet_val,
                    cfg.escape.T_exobase,
                    m_sp;
                    R_exobase=R_exo_val,
                    gamma=cfg.escape.gamma,
                    hydrodynamic=cfg.escape.hydrodynamic,
                )
                M_atm_species[sp] = esc_sp.M_atm
                M_escaped_species[sp] =
                    get(M_escaped_species, sp, 0.0) + esc_sp.M_escaped_step
            end
            acc.M_atm_total = sum(values(M_atm_species))
            acc.M_escaped_total = sum(values(M_escaped_species))
        else
            M_vent_rate_eff = get(vent_rates_esc, cfg.escape.species, 0.0)
            esc_res = evolve_atmospheric_species_inventory(
                acc.M_atm_total,
                M_vent_rate_eff,
                dt,
                M_planet_val,
                rplanet_val,
                cfg.escape.T_exobase,
                get_species_molecular_mass(cfg.escape.species);
                R_exobase=R_exo_val,
                gamma=cfg.escape.gamma,
                hydrodynamic=cfg.escape.hydrodynamic,
            )
            acc.M_atm_total = esc_res.M_atm
            acc.M_escaped_total += esc_res.M_escaped_step
        end
    end

    return nothing
end

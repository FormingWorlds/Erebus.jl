# Radiogenic heating step: radioactive isotope decay power computation and marker heating

"""
    radiogenic_heating!(state::SimulationState, coords::GridCoordinates, cfg::SimulationConfig)::Nothing

Compute short-lived radionuclide (26Al and 60Fe) decay power and update marker radiogenic heating.

# Mutates:
- `state.markers.core.hrtotalm`
"""
function radiogenic_heating!(
    state::SimulationState, coords::GridCoordinates, cfg::SimulationConfig
)::Nothing
    # Empty marker set contract: return immediately with no mutations
    length(state.markers) == 0 && return nothing

    timesum = state.timesum
    (!isfinite(timesum) || timesum < 0.0) &&
        throw(DomainError(timesum, "Decay time must be non-negative, got $timesum"))
    dt = state.dt
    (!isfinite(dt) || dt < 0.0) &&
        throw(DomainError(dt, "Timestep dt must be finite and non-negative, got $dt"))

    hr_al = cfg.thermodynamics.hr_al
    hr_fe = cfg.thermodynamics.hr_fe

    tau_al = cfg.thermodynamics.t_half_al / log(2.0)
    tau_fe = cfg.thermodynamics.t_half_fe / log(2.0)

    has_metal = haskey(state.markers.groups, :metal)
    w_metal_ref = if has_metal
        metal_volume_to_mass_fraction(
            cfg.coreformation.Xfe_bulk,
            cfg.coreformation.rho_metal,
            cfg.materials.rhosolidm[1],
        )
    else
        X_FE_REF_CHONDRITE
    end

    hrsolidm, hrfluidm, hrmetalm = calculate_radioactive_heating(
        hr_al,
        hr_fe,
        timesum;
        dt=dt,
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
        w_metal_ref=w_metal_ref,
    )

    core = state.markers.core
    tm = core.tm
    phim = core.phim
    hrtotalm = core.hrtotalm
    Xfe_bulk = has_metal ? state.markers.groups.metal.Xfe_bulk : nothing

    @inbounds for m in 1:length(tm)
        t_m = tm[m]
        if t_m < 3
            if Xfe_bulk !== nothing
                phi_fe = clamp(Xfe_bulk[m], 0.0, 1.0)
            else
                rho_sol = cfg.materials.rhosolidm[t_m]
                rho_met = cfg.coreformation.rho_metal
                phi_fe = clamp(
                    (X_FE_REF_CHONDRITE / rho_met) /
                    (X_FE_REF_CHONDRITE / rho_met + (1.0 - X_FE_REF_CHONDRITE) / rho_sol),
                    0.0,
                    1.0,
                )
            end
            hr_solid_rock = (1.0 - phi_fe) * hrsolidm[t_m] + phi_fe * hrmetalm[t_m]
            hrtotalm[m] = (1.0 - phim[m]) * hr_solid_rock + phim[m] * hrfluidm[t_m]
        else
            # Sticky air / space
            hrtotalm[m] = 0.0
        end
    end

    return nothing
end

"""
Snapshot pre-pyrolysis marker porosity and bulk iron fraction for mid-step radiogenic heating refresh.

$(SIGNATURES)

Parameters
----------
ws : SimulationWorkspaces
    Preallocated solver workspaces and step-start buffers.
state : SimulationState
    Mutable simulation state container.
cfg : SimulationConfig
    Master simulation configuration.

Returns
-------
Nothing
    Copies `state.markers.core.phim` into `ws.phim_hr_step_start` and `Xfe_bulk` into `ws.Xfe_bulk_hr_step_start`.
"""
function snapshot_radiogenic_step_start!(
    ws::SimulationWorkspaces, state::SimulationState, cfg::SimulationConfig
)::Nothing
    (cfg.thermodynamics.hr_al || cfg.thermodynamics.hr_fe) || return nothing
    _sync_step_start_vec!(ws.phim_hr_step_start, state.markers.core.phim)
    if haskey(state.markers.groups, :metal)
        _sync_step_start_vec!(
            ws.Xfe_bulk_hr_step_start, state.markers.groups.metal.Xfe_bulk
        )
    end
    ws.hr_step_start_timestep = state.timestep
    return nothing
end

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
    timesum < 0.0 &&
        throw(DomainError(timesum, "Decay time must be non-negative, got $timesum"))

    hr_al = cfg.thermodynamics.hr_al
    hr_fe = cfg.thermodynamics.hr_fe

    tau_al = cfg.thermodynamics.t_half_al / log(2.0)
    tau_fe = cfg.thermodynamics.t_half_fe / log(2.0)

    has_metal = haskey(state.markers.groups, :metal)
    X_fe_ref_val = has_metal ? cfg.coreformation.Xfe_bulk : X_FE_REF_CHONDRITE

    hrsolidm, hrfluidm, hrmetalm = calculate_radioactive_heating(
        hr_al,
        hr_fe,
        timesum;
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
        X_fe_ref=X_fe_ref_val,
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

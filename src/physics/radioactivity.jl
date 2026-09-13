"""
Compute radiogenic heat production of isotope mixture.

$(SIGNATURES)

# Details

    - f: fraction of radioactive matter [atoms/kg]
    - ratio: initial ratio of radioactive to non-radioactive isotopes
    - E: heat energy [J]
    - tau: exp decay mean lifetime ``\\tau=\\frac{t_{1/2}}{\\log{2}}`` [s]
    - time: time elapsed since start of radioactive decay [s]

# Returns

    - Q: radiogenic heat production [W/kg]
"""
function Q_radiogenic(f, ratio, E, tau, time)
    if time < 0.0 || tau <= 0.0 || f < 0.0 || ratio < 0.0
        throw(
            DomainError(
                (f, ratio, E, tau, time),
                "Decay time, abundance, and ratio must be non-negative, and lifetime must be positive",
            ),
        )
    end
    return f * ratio * E * exp(-time * inv(tau)) * inv(tau)
end

"""
Compute radiogenic heat production of 26Al and 60Fe isotopes.

$(SIGNATURES)

# Details

    - al: true if radioactive isotope 26Al is present
    - fe: true if radioactive isotope 60Fe is present
    - timesum: time elapsed since initial conditions at start of simulation

# Keyword Arguments

    - ratio_al: initial isotopic ratio of 26Al/27Al
    - E_al: decay energy of 26Al [J]
    - f_al: mass fraction of 27Al in silicate matrix
    - tau_al: mean lifetime of 26Al [s]
    - ratio_fe: initial isotopic ratio of 60Fe/56Fe
    - E_fe: decay energy of 60Fe [J]
    - f_fe: mass fraction of 56Fe in metallic phase
    - tau_fe: mean lifetime of 60Fe [s]
    - rhosolidm: density vector of solid silicate phases [kg/m^3]
    - rhofluidm: density vector of pore fluid phases [kg/m^3]
    - rho_metal: density of metallic phase [kg/m^3] (default 5450.0)

# Returns

    - hrsolidm: radiogenic heat production of 26Al in silicate rock [W/m^3]
    - hrfluidm: radiogenic heat production in pore fluid [W/m^3] (zero: pore water carries no 26Al/60Fe)
    - hrmetalm: radiogenic heat production of 60Fe in metallic iron phase [W/m^3]
"""
function calculate_radioactive_heating(
    al,
    fe,
    timesum;
    ratio_al=ratio_al,
    E_al=E_al,
    f_al=f_al,
    tau_al=tau_al,
    ratio_fe=ratio_fe,
    E_fe=E_fe,
    f_fe=f_fe,
    tau_fe=tau_fe,
    rhosolidm=rhosolidm,
    rhofluidm=rhofluidm,
    rho_metal=5450.0,
)
    # 26Al: planet ✓, crust ✓, space × (lithophile, in silicate rock phase)
    if al
        # 26Al radiogenic heat production [W/kg]
        Q_al = Q_radiogenic(f_al, ratio_al, E_al, tau_al, timesum)
        # Solid phase 26Al radiogenic heat production [W/m^3]
        @inbounds hrsolidm = @SVector [Q_al * rhosolidm[1], Q_al * rhosolidm[2], 0.0]
    else
        hrsolidm = @SVector zeros(3)
    end
    # Fluid phase: pore fluid (water/ice) carries no radiogenic isotopes
    hrfluidm = @SVector zeros(3)
    # 60Fe: planet ✓, crust ✓, space × (siderophile, resides in metallic iron phase)
    if fe
        # 60Fe radiogenic heat production [W/kg]
        Q_fe = Q_radiogenic(f_fe, ratio_fe, E_fe, tau_fe, timesum)
        # Metallic phase 60Fe radiogenic heat production [W/m^3]
        # Silicate mantle (1) and crust (2) metal carries 60Fe; sticky air/space (3) carries zero.
        @inbounds hrmetalm = @SVector [Q_fe * rho_metal, Q_fe * rho_metal, 0.0]
    else
        hrmetalm = @SVector zeros(3)
    end
    return hrsolidm, hrfluidm, hrmetalm
end

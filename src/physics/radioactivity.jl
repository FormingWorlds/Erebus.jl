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
function Q_radiogenic(f, ratio, E, tau, time; dt=0.0)
    if !isfinite(time) ||
        !isfinite(tau) ||
        !isfinite(dt) ||
        time < 0.0 ||
        tau <= 0.0 ||
        f < 0.0 ||
        ratio < 0.0 ||
        dt < 0.0
        throw(
            DomainError(
                (f, ratio, E, tau, time, dt),
                "Decay time, timestep, abundance, and ratio must be finite and non-negative, and lifetime must be positive",
            ),
        )
    end
    Q_inst = f * ratio * E * exp(-time * inv(tau)) * inv(tau)
    return Q_inst * radiogenic_step_average_factor(dt, tau)
end

"""
Compute dimensionless step-average decay factor over `[t, t + dt]`.

$(SIGNATURES)

Parameters
----------
dt : Real
    Timestep duration [s]. Must be finite and non-negative.
tau : Real
    Mean decay lifetime `t_half / log(2)` [s]. Must be finite and positive.

Returns
-------
factor : AbstractFloat
    Step-average multiplier `(1 - exp(-dt / tau)) / (dt / tau)`, or `1` when `dt == 0`.

Raises
------
DomainError
    If `dt < 0`, `tau <= 0`, or either input is non-finite.
"""
@inline function radiogenic_step_average_factor(dt::Real, tau::Real)
    if !isfinite(dt) || !isfinite(tau) || dt < zero(dt) || tau <= zero(tau)
        throw(
            DomainError(
                (dt, tau),
                "Timestep dt must be finite and non-negative and lifetime tau must be finite and positive",
            ),
        )
    end
    T = promote_type(typeof(float(dt)), typeof(float(tau)))
    dt_f = T(dt)
    tau_f = T(tau)
    iszero(dt_f) && return one(T)
    x = dt_f * inv(tau_f)
    if x < T(1.0e-8)
        return evalpoly(x, (one(T), -inv(T(2)), inv(T(6)), -inv(T(24))))
    end
    return -expm1(-x) * inv(x)
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
    - f_al: abundance of 27Al [atoms/kg bulk rock]
    - tau_al: mean lifetime of 26Al [s]
    - ratio_fe: initial isotopic ratio of 60Fe/56Fe
    - E_fe: decay energy of 60Fe [J]
    - f_fe: abundance of 56Fe [atoms/kg bulk rock]
    - tau_fe: mean lifetime of 60Fe [s]
    - rhosolidm: density vector of solid silicate phases [kg/m^3]
    - rhofluidm: density vector of pore fluid phases [kg/m^3]
    - rho_metal: density of metallic phase [kg/m^3] (default 5450.0)
    - w_metal_ref: metal mass fraction of the body at its initial composition (default
      `X_FE_REF_CHONDRITE`, which places all CI iron in the metal)

# Returns

    - hrsolidm: radiogenic heat production of 26Al, and of 60Fe in oxidised iron, in silicate
      rock [W/m^3]
    - hrfluidm: radiogenic heat production in pore fluid [W/m^3] (zero: pore water carries no 26Al/60Fe)
    - hrmetalm: radiogenic heat production of 60Fe in metallic iron phase [W/m^3]
"""
function calculate_radioactive_heating(
    al,
    fe,
    timesum;
    dt=0.0,
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
    w_metal_ref=X_FE_REF_CHONDRITE,
)
    (isfinite(w_metal_ref) && 0.0 <= w_metal_ref <= 1.0) ||
        throw(DomainError(w_metal_ref, "Reference metal mass fraction must be in [0, 1]"))
    # f_al and f_fe are CI abundances per kg bulk rock. Metal up to the CI iron content comes
    # out of a CI bulk; metal above it adds iron to a body with CI silicate.
    w_ci = X_FE_REF_CHONDRITE
    w_metal_ci = min(w_metal_ref, w_ci)
    w_silicate_ci = 1.0 - w_metal_ci
    # 26Al: planet ✓, crust ✓, space × (lithophile, in silicate rock phase)
    # Concentrating lithophile 26Al into silicate phase [W/kg silicate]
    Q_silicate = 0.0
    if al
        Q_silicate +=
            Q_radiogenic(f_al, ratio_al, E_al, tau_al, timesum; dt=dt) / w_silicate_ci
    end
    # Fluid phase: pore fluid (water/ice) carries no radiogenic isotopes
    hrfluidm = @SVector zeros(3)
    # 60Fe: planet ✓, crust ✓, space × (follows iron atoms: metal, and oxidised iron in silicate)
    if fe
        # 60Fe radiogenic heat production [W/kg bulk rock]
        Q_fe = Q_radiogenic(f_fe, ratio_fe, E_fe, tau_fe, timesum; dt=dt)
        # f_fe counts 56Fe per kg bulk CI, so one kg of iron holds f_fe / w_ci.
        Q_fe_metal = Q_fe / w_ci
        Q_silicate += Q_fe_metal * (w_ci - w_metal_ci) / w_silicate_ci
        @inbounds hrmetalm = @SVector [Q_fe_metal * rho_metal, Q_fe_metal * rho_metal, 0.0]
    else
        hrmetalm = @SVector zeros(3)
    end
    # Solid silicate heat production [W/m^3]
    @inbounds hrsolidm = @SVector [
        Q_silicate * rhosolidm[1], Q_silicate * rhosolidm[2], 0.0
    ]
    return hrsolidm, hrfluidm, hrmetalm
end

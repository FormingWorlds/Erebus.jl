"""
    compute_drained_compressibility(betaphi, phi, betasolid; phimin=phimin, phimax=phimax)

Compute drained bulk compressibility of a porous medium.

# Details

    β_d = (β_ϕ + β_s) / (1 - ϕ)

where β_ϕ is pore compressibility [1/Pa], β_s is solid matrix compressibility [1/Pa],
and ϕ is porosity [-]. References: Biot (1941), Detournay & Cheng (1993), Gerya (2019).

# Arguments

    - betaphi: pore compressibility β_ϕ [1/Pa]
    - phi: porosity ϕ [-]
    - betasolid: solid matrix compressibility β_s [1/Pa]
    - phimin: minimum porosity limit
    - phimax: maximum porosity limit

# Returns

    - betadrained: drained bulk compressibility β_d [1/Pa]
"""
function compute_drained_compressibility(
    betaphi::T1, phi::T2, betasolid::T3; phimin::Real=phimin, phimax::Real=phimax
) where {T1<:Real,T2<:Real,T3<:Real}
    T = promote_type(T1, T2, T3)
    bphi = max(betaphi, zero(T))
    bsolid = max(betasolid, zero(T))
    phi_eff = clamp(phi, T(phimin), T(phimax))
    return (bphi + bsolid) / (one(T) - phi_eff)
end

"""
    compute_biot_willis_coefficient(betadrained, betasolid)

Compute Biot-Willis coefficient for poroelastic coupling.

# Details

    K_BW = 1 - β_s / β_d

For an incompressible solid matrix (β_s = 0), K_BW = 1.
For intact zero-porosity rock (β_d → β_s), K_BW → 0.
Physical bounds: K_BW ∈ [0, 1]. References: Biot (1941), Wang (2000).

# Arguments

    - betadrained: drained bulk compressibility β_d [1/Pa]
    - betasolid: solid matrix compressibility β_s [1/Pa]

# Returns

    - kbw: Biot-Willis coefficient K_BW [-]
"""
function compute_biot_willis_coefficient(
    betadrained::T1, betasolid::T2
) where {T1<:Real,T2<:Real}
    T = promote_type(T1, T2)
    if betasolid <= zero(T)
        return one(T)
    end
    if betadrained <= betasolid
        return zero(T)
    end
    return clamp(one(T) - betasolid / betadrained, zero(T), one(T))
end

"""
    compute_skempton_coefficient(betadrained, phi, betasolid, betafluid; phimin=phimin, phimax=phimax)

Compute Skempton coefficient B for pore pressure response to mean stress.

# Details

    B = (β_d - β_s) / (β_d - β_s + ϕ * (β_f - β_s))

For incompressible constituents (β_s = 0, β_f = 0), B = 1.
Physical bounds: B ∈ [0, 1]. References: Skempton (1954), Rice & Cleary (1976).

# Arguments

    - betadrained: drained bulk compressibility β_d [1/Pa]
    - phi: porosity ϕ [-]
    - betasolid: solid matrix compressibility β_s [1/Pa]
    - betafluid: pore fluid compressibility β_f [1/Pa]
    - phimin: minimum porosity limit
    - phimax: maximum porosity limit

# Returns

    - ksk: Skempton coefficient B [-]
"""
function compute_skempton_coefficient(
    betadrained::T1,
    phi::T2,
    betasolid::T3,
    betafluid::T4;
    phimin::Real=phimin,
    phimax::Real=phimax,
) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real}
    T = promote_type(T1, T2, T3, T4)
    if betasolid <= zero(T) && betafluid <= zero(T)
        return one(T)
    end
    bsolid = max(betasolid, zero(T))
    bfluid = max(betafluid, zero(T))
    phi_eff = clamp(phi, T(phimin), T(phimax))
    num = betadrained - bsolid
    denom = num + phi_eff * (bfluid - bsolid)
    if denom <= zero(T) || num <= zero(T)
        return one(T)
    end
    return clamp(num / denom, zero(T), one(T))
end

"""
    compute_rhofluid(T::Real, rho0::Real, alpha::Real, T0::Real; thermal_buoyancy::Bool = true)

Compute temperature-dependent pore fluid density with volumetric thermal expansion:
    ρ_f(T) = ρ_{f0} * max(0.1, 1.0 - α_f * (T - T_0))   for T > T_0

When `thermal_buoyancy = false`, returns reference density `rho0` unmodified.
Density is clamped to a lower bound of `0.1 * rho0` to prevent unphysical negative values
at extreme temperatures.

# Arguments

    - T: temperature [K]
    - rho0: reference fluid density at T0 [kg/m³]
    - alpha: fluid volumetric thermal expansion coefficient α_f [1/K]
    - T0: reference temperature [K]
    - thermal_buoyancy: enable or disable thermal expansion (default: true)

# Returns

    - rhof: temperature-dependent fluid density [kg/m³]
"""
function compute_rhofluid(
    T::Real, rho0::Real, alpha::Real, T0::Real; thermal_buoyancy::Bool=true
)
    if !isfinite(T) || !thermal_buoyancy || alpha <= 0.0 || T <= T0
        return Float64(rho0)
    end
    factor = max(0.1, 1.0 - alpha * (T - T0))
    return isnan(factor) ? Float64(rho0) : Float64(rho0 * factor)
end

"""
    compute_fluid_viscosity(T::Real, tm::Integer;
                            mode::Symbol = :arrhenius,
                            eta0::Real = 1.0e-3,
                            eta_ice::Real = 1.0e12,
                            eta_air::Real = 1.0e-3,
                            Ea::Real = 15.0e3,
                            T0::Real = 293.15,
                            tmfluidphase::Real = 273.0,
                            etamin::Real = 1.0e-5,
                            etamax::Real = 1.0e12)

Compute temperature-dependent dynamic fluid viscosity η_f(T) [Pa s].

For sticky air markers (`tm >= 3`), returns `eta_air`.
For sub-freezing rock markers (`T <= tmfluidphase`), returns `eta_ice`.
For non-finite or corrupt temperatures (`!isfinite(T)`), returns `eta_ice` to prevent runaway mobility.
For liquid fluid markers (`T > tmfluidphase`):
  - `:constant` mode: returns `eta0`.
  - `:arrhenius` mode:
      η_f(T) = η_0 * exp((E_a / R) * (1/T - 1/T_0))
    clamped to [etamin, etamax].

# Arguments

    - T: temperature [K]
    - tm: material type index (1: core, 2: crust, 3: sticky air)
    - mode: `:arrhenius` or `:constant` (default: `:arrhenius`)
    - eta0: reference liquid fluid viscosity at T0 [Pa s]
    - eta_ice: sub-freezing ice viscosity [Pa s]
    - eta_air: sticky air fluid viscosity [Pa s]
    - Ea: activation energy for fluid viscous flow [J/mol]
    - T0: reference temperature [K]
    - tmfluidphase: melting temperature [K]
    - etamin: minimum viscosity floor [Pa s]
    - etamax: maximum viscosity ceiling [Pa s]

# Returns

    - etafluid: dynamic fluid viscosity [Pa s]
"""
function compute_fluid_viscosity(
    T::Real,
    tm::Integer;
    mode::Symbol=:arrhenius,
    eta0::Real=1.0e-3,
    eta_ice::Real=1.0e12,
    eta_air::Real=1.0e-3,
    Ea::Real=15.0e3,
    T0::Real=293.15,
    tmfluidphase::Real=273.0,
    etamin::Real=1.0e-5,
    etamax::Real=1.0e12,
)
    if tm >= 3
        return Float64(eta_air)
    end
    if !isfinite(T)
        return Float64(eta_ice)
    end
    if T <= tmfluidphase
        return Float64(eta_ice)
    end
    if mode === :constant || Ea <= 0.0
        return Float64(eta0)
    elseif mode === :arrhenius
        log_ratio = (Ea / R_GAS) * (inv(T) - inv(T0))
        val = eta0 * exp(log_ratio)
        return clamp(val, Float64(etamin), Float64(etamax))
    else
        throw(
            ArgumentError(
                "Unknown fluid viscosity mode: $mode (expected :arrhenius or :constant)"
            ),
        )
    end
end

"""
    hydrofracture_overpressure_ramp(x::Real, delta::Real=0.0) -> Float64

Evaluate C¹ continuous overpressure regularisation ramp for dynamic hydrofracturing.

# Arguments
- `x::Real`: Normalised tensile fluid overpressure `(-Peff - sigma_t) / sigma_t` [-].
- `delta::Real`: Regularisation ramp half-width delta in [0, 1] [-] (default: 0.0).

# Returns
- `Float64`: Regularised effective overpressure value `s(x)`.
"""
@inline function hydrofracture_overpressure_ramp(x::Real, delta::Real=0.0)
    d = Float64(delta)
    if !isfinite(d) || d < 0.0
        throw(DomainError(delta, "ramp_width delta must be finite and >= 0"))
    end
    xv = Float64(x)
    if isnan(xv)
        return NaN
    elseif xv <= 0.0
        return 0.0
    elseif xv == Inf
        return Inf
    end
    if d == 0.0 || xv >= d
        return xv - 0.5 * d
    else
        return (xv^2) / (2.0 * d)
    end
end

"""
    compute_hydrofracture_factor(Peff::Real, sigma_t::Real;
                                 active::Bool = true,
                                 kappa_frac::Real = 1.0e3,
                                 gamma::Real = 1.0,
                                 max_factor::Real = Inf,
                                 ramp_width::Real = 0.0)

Compute dimensionless permeability enhancement factor from dynamic hydrofracturing.

When pore fluid pressure exceeds total confining pressure plus tensile strength
(Terzaghi effective pressure Peff = Pt - Pf <= -sigma_t), hydraulic tensile
fractures open and increase effective permeability:

    factor = 1.0 + kappa_frac * s(x)^gamma

clamped to [1.0, max_factor].
"""
@inline function compute_hydrofracture_factor(
    Peff::Real,
    sigma_t::Real;
    active::Bool=true,
    kappa_frac::Real=1.0e3,
    gamma::Real=1.0,
    max_factor::Real=Inf,
    ramp_width::Real=0.0,
)
    if !active || !isfinite(Peff) || !isfinite(sigma_t) || sigma_t <= 0.0
        return 1.0
    end
    overpressure = -Peff - sigma_t
    norm_overpressure = overpressure / sigma_t
    s = hydrofracture_overpressure_ramp(norm_overpressure, ramp_width)
    if s <= 0.0
        return 1.0
    end
    sg = gamma == 1.0 ? s : (gamma == 2.0 ? abs2(s) : (s ^ gamma))
    factor = 1.0 + kappa_frac * sg
    return clamp(Float64(factor), 1.0, Float64(max_factor))
end

"""
    compute_hydrofracture_permeability(kphi::Real, Peff::Real, sigma_t::Real;
                                       active::Bool = true,
                                       kappa_frac::Real = 1.0e3,
                                       gamma::Real = 1.0,
                                       kmax::Real = 1.0e-9,
                                       ramp_width::Real = 0.0)

Compute effective permeability k_eff [m²] with dynamic hydrofracturing enhancement.

When pore fluid pressure exceeds total confining pressure plus rock tensile strength:

    Peff = Pt - Pf <= -sigma_t

tensile microcracks open and enhance permeability according to:

    k_eff = min(kphi * compute_hydrofracture_factor(Peff, sigma_t; active, kappa_frac, gamma, ramp_width), kmax)

# Arguments
- `kphi`: baseline matrix permeability [m²]
- `Peff`: Terzaghi effective pressure Pt - Pf [Pa]
- `sigma_t`: rock tensile strength [Pa]
- `active`: enable hydrofracture enhancement (default: true)
- `kappa_frac`: dimensionless enhancement multiplier (default: 1.0e3)
- `gamma`: power-law exponent (default: 1.0)
- `kmax`: maximum permeability ceiling [m²] (default: 1.0e-9)
- `ramp_width`: C¹ regularisation ramp half-width [0, 1] [-] (default: 0.0)

# Returns
- `k_eff`: effective permeability [m²]
"""
function compute_hydrofracture_permeability(
    kphi::Real,
    Peff::Real,
    sigma_t::Real;
    active::Bool=true,
    kappa_frac::Real=1.0e3,
    gamma::Real=1.0,
    kmax::Real=1.0e-9,
    ramp_width::Real=0.0,
)
    if kphi < 0.0
        throw(
            DomainError(
                kphi, "Negative matrix permeability kphi is not physically allowed."
            ),
        )
    end
    if !active || !isfinite(Peff) || !isfinite(sigma_t) || sigma_t <= 0.0 || kphi == 0.0
        return Float64(kphi)
    end
    factor = compute_hydrofracture_factor(
        Peff,
        sigma_t;
        active=active,
        kappa_frac=kappa_frac,
        gamma=gamma,
        ramp_width=ramp_width,
    )
    k_enhanced = Float64(kphi * factor)
    kphi_f = Float64(kphi)
    kmax_f = Float64(kmax)
    return min(max(k_enhanced, kphi_f), max(kphi_f, kmax_f))
end

"""
    evaluate_hydrofracture_resistance(
        r_base::Real,
        Peff::Real,
        sigma_t::Real,
        phi_face::Real;
        kphim0::Union{Nothing,Real,AbstractVector}=nothing,
        phim0_val::Union{Nothing,Real}=nothing,
        phimin_val::Real=1.0e-4,
        kappa_frac::Real=1.0e3,
        gamma_frac::Real=1.0,
        k_frac_max::Real=1.0e-9,
        ramp_width::Real=0.0,
        rx_floor_prefactor::Real=1.0e-5,
    )

Evaluate hydrofracture-enhanced Darcy resistance R [Pa s / m²] using matrix permeability
derived from face porosity and reference permeability.
"""
function evaluate_hydrofracture_resistance(
    r_base::Real,
    Peff::Real,
    sigma_t::Real,
    phi_face::Real;
    kphim0::Union{Nothing,Real,AbstractVector}=nothing,
    phim0_val::Union{Nothing,Real}=nothing,
    phimin_val::Real=1.0e-4,
    kappa_frac::Real=1.0e3,
    gamma_frac::Real=1.0,
    k_frac_max::Real=1.0e-9,
    ramp_width::Real=0.0,
    rx_floor_prefactor::Real=1.0e-5,
)
    phi_eff = clamp(Float64(phi_face), Float64(phimin_val), 1.0 - 1.0e-6)
    kphim0_ref = if kphim0 !== nothing
        kphim0 isa AbstractVector ? Float64(kphim0[1]) : Float64(kphim0)
    else
        1.0e-13
    end
    phim0_eff = phim0_val !== nothing ? Float64(phim0_val) : 0.01
    kphi_val = max(0.0, kphi(kphim0_ref, phi_eff; phim0=phim0_eff))
    r_floor = rx_floor_prefactor / k_frac_max
    if kphi_val > 0.0
        keff = compute_hydrofracture_permeability(
            kphi_val,
            Peff,
            sigma_t;
            active=true,
            kappa_frac=kappa_frac,
            gamma=gamma_frac,
            kmax=k_frac_max,
            ramp_width=ramp_width,
        )
        return max(Float64(r_base) * (kphi_val / keff), r_floor)
    else
        ffrac = compute_hydrofracture_factor(
            Peff,
            sigma_t;
            active=true,
            kappa_frac=kappa_frac,
            gamma=gamma_frac,
            ramp_width=ramp_width,
        )
        return max(Float64(r_base) / ffrac, r_floor)
    end
end

"""
    phi_cold(P_bar)

Porosity of a granular silicate powder after cold isostatic pressing.

# Details

    φ_H(P) = 0.42 + 0.46 / ((P / 0.13 bar)^1.72 + 1)

Eq. 15 of Henke et al. (2012), the laboratory fit of Güttler et al. (2009, in Henke et
al. 2012), with the relative density D = 1 - φ, as used by Lichtenberg et al. (2016a)
Eq. 10. The marker step applies it to the highest pressure the marker has felt before it
starts to sinter.

# Arguments

    - P_bar: pressure P [bar], >= 0

# Returns

    - φ_H: porosity [-], in (0.42, 0.88]

# Throws

    - DomainError: `P_bar` is negative or not finite
"""
function phi_cold(P_bar::Real)
    require_nonneg_finite(P_bar, "pressure P_bar")
    return 0.42 + 0.46 / ((P_bar / 0.13)^1.72 + 1.0)
end

"""
    porosity_phi0(P_bar, phi_init; phi_cap=0.85)

Porosity of a rock marker before hot sintering: the smallest of the initial porosity,
the cold pressing porosity `phi_cold(P_bar)` and the upper cap `phi_cap`.

# Arguments

    - P_bar: pressure P [bar], >= 0
    - phi_init: initial porosity φ_init [-], in [0, 1]
    - phi_cap: upper porosity cap [-], in (0, 1) as `porosity.phi_cap` (i2elvis value 0.85)

# Returns

    - φ_0 = min(φ_init, φ_H(P), φ_cap) [-]

# Throws

    - DomainError: φ_init outside [0, 1], φ_cap outside (0, 1) or an invalid pressure
"""
function porosity_phi0(P_bar::Real, phi_init::Real; phi_cap::Real=0.85)
    require_unit_interval(phi_init, "initial porosity phi_init")
    (isfinite(phi_cap) && 0.0 < phi_cap < 1.0) ||
        throw(DomainError(phi_cap, "porosity cap phi_cap must be in (0, 1)"))
    return min(Float64(phi_init), phi_cold(P_bar), Float64(phi_cap))
end

"""
    grain_geometry(phi, phi_0, R0)

Radius of a truncated grain and radius of one of its contact faces after hot sintering
from porosity `phi_0` to porosity `phi`, for 8 contacts per grain.

# Details

Volume conservation of a sphere of radius R0 with 8 cut caps (Henke et al. 2012
Eqs. 34 to 37, after Kakar and Chaklader 1967 in Henke et al. 2012) gives, with
X = (1 - φ_0) / (1 - φ),

    R^3 - 3 R^2 R0 X^(1/3) + R0^3 X + R0^3 = 0,    a^2 = R^2 - R0^2 X^(2/3).

The physical root is the one in [R0, 2 R0 X^(1/3)], where the cubic decreases. It is
found by Newton iteration from R0 inside a bisection bracket. For X <= 1/3 the contact
faces overlap and the geometry has no root; X is then held at 1/3 + 1e-12.

# Arguments

    - phi: porosity φ [-], in [0, φ_0]
    - phi_0: porosity before sintering φ_0 [-], in [0, 1)
    - R0: initial grain radius [m], > 0

# Returns

    - (R, a): grain radius R [m] and contact face radius a [m]; a = 0 at φ = φ_0

# Throws

    - DomainError: φ or φ_0 outside [0, 1), φ > φ_0, or R0 not positive and finite
"""
function grain_geometry(phi::Real, phi_0::Real, R0::Real)
    require_positive_finite(R0, "grain radius R0")
    (isfinite(phi) && 0.0 <= phi < 1.0) ||
        throw(DomainError(phi, "porosity phi must be in [0, 1)"))
    (isfinite(phi_0) && 0.0 <= phi_0 < 1.0) ||
        throw(DomainError(phi_0, "porosity phi_0 must be in [0, 1)"))
    phi <= phi_0 || throw(DomainError(phi, "porosity phi must not exceed phi_0"))
    n_contacts = 8
    X = max((1.0 - phi_0) / (1.0 - phi), 1.0 / 3.0 + 1.0e-12)
    c = cbrt(X)
    # Eq. 36 with Eqs. 34, 35, 37: R0^3 = R^3 - Z/8 (2 R^3 - 3 R^2 R0 c + R0^3 X)
    f(R) = n_contacts / 8 * (2R^3 - 3R^2 * R0 * c + R0^3 * X) - R^3 + R0^3
    df(R) = n_contacts / 8 * (6R^2 - 6R * R0 * c) - 3R^2
    lo, hi = Float64(R0), 2.0 * R0 * c
    R = lo
    for _ in 1:200
        fR = f(R)
        iszero(fR) && break
        fR > 0.0 ? (lo = R) : (hi = R)
        dR = df(R)
        R_new = R - fR / dR
        if !(lo < R_new < hi) || dR >= 0.0
            R_new = 0.5 * (lo + hi)
        end
        converged = abs(R_new - R) <= 4.0 * eps(R)
        R = R_new
        (converged || hi - lo <= 4.0 * eps(hi)) && break
    end
    return R, sqrt(max(R^2 - R0^2 * c^2, 0.0))
end

"""
    contact_stress(P_bar, R, a)

Effective stress on the contact faces of the grains (Henke et al. 2012 Eqs. 40, 41).

# Details

    C_av = 2 sqrt(3) (R^2 - a^2),    σ = C_av P / (π a^2)

`C_av` is the mean cross-section of the cell of one grain in the hexagonal prismatic
packing. Early in sintering a << R and σ >> P.

# Arguments

    - P_bar: applied pressure P [bar], >= 0
    - R: grain radius [m], > 0
    - a: contact face radius [m], in (0, R]

# Returns

    - σ: contact stress [bar]

# Throws

    - DomainError: P < 0, R <= 0, a <= 0, a > R, or a value that is not finite
"""
function contact_stress(P_bar::Real, R::Real, a::Real)
    require_nonneg_finite(P_bar, "pressure P_bar")
    require_positive_finite(R, "grain radius R")
    require_positive_finite(a, "contact radius a")
    a <= R || throw(DomainError(a, "contact radius a must not exceed the grain radius R"))
    C_av = 2.0 * sqrt(3.0) * (R^2 - a^2)
    return C_av * P_bar / (pi * a^2)
end

"""
    sintering_rate(phi, sigma_bar, R, T; A=4.0e-5, E_act=85000.0, n=1.5)

Rate of porosity decrease by hot sintering (power-law creep of olivine of Schwenn and
Goetze 1978 in Henke et al. 2012, Eqs. 39, 43, 44, and Lichtenberg et al. 2016a Eq. 13).

# Details

    r = F(φ) A (1 - φ) σ^n / (100 R)^3 exp(-E_act / (R_gas T)),    dφ/dt = -r

with R in cm (100 R for R in m), σ in bar, R_gas = 1.9872 cal/(K mol) and
F(φ) = 10 φ for φ <= 0.1, else 1 (Eq. 43, pore closure stage). The rate is >= 0.

# Arguments

    - phi: porosity φ [-], in [0, 1]
    - sigma_bar: contact stress σ [bar], >= 0
    - R: grain radius [m], > 0
    - T: temperature [K], > 0
    - A: creep prefactor [cm^3 bar^-n s^-1]
    - E_act: activation energy [cal/mol]
    - n: stress exponent [-], > 0

# Returns

    - r: porosity decrease rate [1/s], >= 0

# Throws

    - DomainError: an argument outside its physical domain
"""
function sintering_rate(
    phi::Real,
    sigma_bar::Real,
    R::Real,
    T::Real;
    A::Real=4.0e-5,
    E_act::Real=85000.0,
    n::Real=1.5,
)
    require_unit_interval(phi, "porosity phi")
    require_nonneg_finite(sigma_bar, "contact stress sigma_bar")
    require_positive_finite(R, "grain radius R")
    require_positive_finite(T, "temperature T")
    require_nonneg_finite(A, "creep prefactor A")
    require_nonneg_finite(E_act, "activation energy E_act")
    require_positive_finite(n, "stress exponent n")
    R_gas = 1.9872
    F = phi <= 0.1 ? 10.0 * phi : 1.0
    return F * A * (1.0 - phi) * sigma_bar^n / (100.0 * R)^3 * exp(-E_act / (R_gas * T))
end

"""
    porosity_decrease_rate(phi, phi_0, P_bar, T; R0=1.0e-6, A=4.0e-5, E_act=85000.0, n=1.5)

Hot sintering rate at porosity `phi` for a rock that started sintering at `phi_0`:
`grain_geometry`, then `contact_stress`, then `sintering_rate`.

# Arguments

    - phi: porosity φ [-], in [0, φ_0]
    - phi_0: porosity before sintering φ_0 [-], in [0, 1)
    - P_bar: pressure P [bar], >= 0
    - T: temperature [K], > 0
    - R0: initial grain radius [m]
    - A, E_act, n: creep law constants (see `sintering_rate`)

# Returns

    - r: porosity decrease rate [1/s], >= 0; 0 for P = 0 or φ = 0; Inf at φ = φ_0 > 0
      with P > 0 (unsintered grains touch at points)

# Throws

    - DomainError: an argument outside its physical domain
"""
function porosity_decrease_rate(
    phi::Real,
    phi_0::Real,
    P_bar::Real,
    T::Real;
    R0::Real=1.0e-6,
    A::Real=4.0e-5,
    E_act::Real=85000.0,
    n::Real=1.5,
)
    require_nonneg_finite(P_bar, "pressure P_bar")
    require_positive_finite(T, "temperature T")
    R, a = grain_geometry(phi, phi_0, R0)
    (iszero(P_bar) || iszero(phi)) && return 0.0
    iszero(a) && return Inf
    sigma = contact_stress(P_bar, R, a)
    return sintering_rate(phi, sigma, R, T; A=A, E_act=E_act, n=n)
end

"""
Porosity decrease [-] below which a decrease is numerical, not sintering.

Near the onset of sintering (X = 1) the rounding of 1 - φ moves φ by a few 1e-16, a
decrease that the sintering law does not produce. A rock marker counts as sintered only
once its cumulative sintering decrease `dphi_sinter` exceeds this value, and
`sinter_porosity_step` needs a largest substep decrease `dphi_max` of at least this value.
"""
const DPHI_SINTERED = 1.0e-12

"""
Backward Euler porosity over one substep `h` [s]: the root of
g(φ) = φ - φ_old + h r(φ) in (0, φ_old] by bisection; g(0) < 0 since r(0) = 0.
The bisection keeps g(lo) < 0 <= g(hi) until lo and hi are adjacent floats (at most
1100 halvings reach the smallest subnormal from 1) and returns `hi`, the end whose
decrease φ_old - φ does not exceed h r(φ). A decrease below 1 ULP of φ_old therefore
gives φ_old, and so does h <= 0.
"""
function _backward_euler_porosity(phi_old::Float64, h::Float64, rate)
    (h <= 0.0 || h * rate(phi_old) <= 0.0) && return phi_old
    lo, hi = 0.0, phi_old
    for _ in 1:1100
        mid = 0.5 * (lo + hi)
        (mid <= lo || mid >= hi) && break
        if mid - phi_old + h * rate(mid) < 0.0
            lo = mid
        else
            hi = mid
        end
    end
    return hi
end

"""
    sinter_porosity_step(phi_n, phi_0, P_bar, T, dt; R0=1.0e-6, A=4.0e-5,
                         E_act=85000.0, n=1.5, dphi_max=0.01, max_substeps=1000)

Advance the porosity of one rock by hot sintering over a time step with sub-stepped
backward Euler.

# Details

Each substep solves φ^{k+1} = φ^k - h r(φ^{k+1}) for φ^{k+1} in (0, φ^k] by bisection.
The substep h adapts so that each substep lowers φ by at most `dphi_max`. After
`max_substeps - 1` substeps the last substep takes the remaining time and the step is
reported as capped. P = 0, φ^n = 0 or dt = 0 return φ^n. φ never increases.

# Arguments

    - phi_n: porosity at the start of the step [-], in [0, φ_0]
    - phi_0: porosity before sintering [-], in [0, 1)
    - P_bar: pressure [bar], >= 0
    - T: temperature [K], > 0
    - dt: time step [s], >= 0
    - R0, A, E_act, n: grain radius [m] and creep law constants (see `sintering_rate`)
    - dphi_max: largest porosity decrease per substep [-], >= `DPHI_SINTERED` (1e-12),
      so that a substep can step over the rounding of 1 - φ near φ_0
    - max_substeps: largest number of substeps, >= 1

# Returns

    - (phi, nsub, capped): porosity at the end of the step [-], number of substeps,
      and `true` when the substep cap was reached before the end of the step

# Throws

    - DomainError: an argument outside its physical domain, `dphi_max` below
      `DPHI_SINTERED`, or a rate so large that a rejected substep cannot shrink further
      (the value is a NamedTuple of the inputs)
    - ArgumentError: `max_substeps` < 1
"""
function sinter_porosity_step(
    phi_n::Real,
    phi_0::Real,
    P_bar::Real,
    T::Real,
    dt::Real;
    R0::Real=1.0e-6,
    A::Real=4.0e-5,
    E_act::Real=85000.0,
    n::Real=1.5,
    dphi_max::Real=0.01,
    max_substeps::Integer=1000,
)
    require_nonneg_finite(dt, "time step dt")
    (isfinite(dphi_max) && dphi_max >= DPHI_SINTERED) || throw(
        DomainError(dphi_max, "substep porosity change dphi_max must be >= DPHI_SINTERED"),
    )
    max_substeps >= 1 || throw(ArgumentError("max_substeps must be >= 1"))
    rate(p) = porosity_decrease_rate(p, phi_0, P_bar, T; R0=R0, A=A, E_act=E_act, n=n)
    rate(phi_n)
    phi = Float64(phi_n)
    (iszero(P_bar) || iszero(phi) || iszero(dt)) && return (phi, 0, false)
    t_left = Float64(dt)
    h = t_left
    nsub = 0
    capped = false
    while t_left > 0.0
        h = min(h, t_left)
        last = nsub + 1 == max_substeps
        capped = last && h < t_left
        h = last ? t_left : h
        phi_try = _backward_euler_porosity(phi, h, rate)
        dphi = phi - phi_try
        if !last && dphi > dphi_max
            # A rejected substep must shrink; a subnormal h times a factor above 0.5
            # rounds back to h, and an overflowing rate rejects every h
            h_new = h * max(0.1, 0.9 * dphi_max / dphi)
            (0.0 < h_new < h) || throw(
                DomainError(
                    (; phi_n, phi_0, P_bar, T, dt, R0, A, E_act, n),
                    "sintering rate overflows: a rejected substep does not shrink",
                ),
            )
            h = h_new
            continue
        end
        phi = phi_try
        t_left = h >= t_left ? 0.0 : t_left - h
        nsub += 1
        h *= min(2.0, 0.9 * dphi_max / max(dphi, floatmin(Float64)))
    end
    return (phi, nsub, capped)
end

"""
    k_porous(k, phi)

Thermal conductivity of porous silicate rock (Lichtenberg et al. 2016a Eq. 12, the
3-branch form of the i2elvis code, fitted by Henke et al. 2012 Eqs. 28 and 29).

# Details

    k1 = k exp(-φ / 0.08),    k2 = k exp(-1.2 - φ / 0.167)
    k_por = k1 (φ < 0.2),  (k1^4 + k2^4)^(1/4) (0.2 <= φ <= 0.4),  k2 (φ > 0.4)

with a floor of 1e-3 W/(m K). The form is not continuous at φ = 0.2 and 0.4.

# Arguments

    - k: conductivity of the dense solid [W/(m K)], >= 0
    - phi: porosity φ [-], in [0, 1]

# Returns

    - k_por: conductivity of the porous rock [W/(m K)], >= 1e-3

# Throws

    - DomainError: k negative or not finite, or φ outside [0, 1]
"""
function k_porous(k::Real, phi::Real)
    require_nonneg_finite(k, "solid conductivity k")
    require_unit_interval(phi, "porosity phi")
    k1 = k * exp(-phi / 0.08)
    k2 = k * exp(-1.2 - phi / 0.167)
    k_eff = if phi < 0.2
        k1
    elseif phi <= 0.4
        (k1^4 + k2^4)^0.25
    else
        k2
    end
    return max(k_eff, 1.0e-3)
end

"""
    dry_porosity_update(phi_d, P_bar, T, dt, pcfg::PorosityConfig; phi0_prev=1.0,
                        sintered=false)

Advance the dry porosity of one rock marker over one time step: cold pressing until the
marker starts to sinter, then hot sintering for `T >= pcfg.T_sinter`.

# Details

Before the onset of sintering (`sintered = false`), the porosity before sintering is
φ_0 = min(φ_0,prev, φ_init, φ_H(P), φ_cap). With `phi0_prev` the φ_0 of the previous
step, φ_0 belongs to the highest pressure the marker has felt, because φ_H decreases with
P; a pressure drop leaves it unchanged. Cold pressing sets φ <- min(φ, φ_0), so φ <= φ_0
before sintering starts.

Once the marker has sintered (`sintered = true`), φ_0 stays at `phi0_prev` and cold
pressing stops: φ_0 is the porosity of the initial packing of Henke et al. (2012, Eq. 37),
at which the contact faces have zero radius. The caller sets `sintered` once the
cumulative sintering decrease exceeds `DPHI_SINTERED` (1e-12): a smaller decrease is the
rounding of 1 - φ near the onset, not sintering. For a marker that starts to sinter at
φ = φ_0, as every marker of the setup does, the cumulative sintering decrease then equals
φ_0 - φ to within `DPHI_SINTERED`. For a marker that is below its φ_0 at the onset of
sintering (for example from `push_marker!`), φ_0 - φ exceeds the sintering decrease by
the gap φ_0 - φ at the onset of sintering, after any cold pressing. The i2elvis code
instead recomputes φ_0 from the current pressure every step.

Sintering lowers φ with `sinter_porosity_step`. Porosity never increases: a marker that
cools keeps its sintered porosity.

# Arguments

    - phi_d: dry porosity at the start of the step [-], in [0, 1]
    - P_bar: pressure [bar], >= 0
    - T: temperature [K], > 0
    - dt: time step [s], >= 0
    - pcfg: `PorosityConfig` with φ_init, φ_cap, R0, T_sinter, A, E_act, n
    - phi0_prev: φ_0 of the previous step [-], in [0, 1]; 1 means no pressure history,
      so a sintered marker needs phi0_prev < 1 after conversion to Float64
    - sintered: `true` once earlier steps have lowered φ by sintering by more than
      `DPHI_SINTERED`

# Returns

    - NamedTuple `(phi, phi0, dphi_sinter, capped)`: dry porosity at the end of the step,
      porosity before sintering φ_0, the decrease by sintering (not by cold pressing),
      and whether the substep cap of `sinter_porosity_step` was reached

# Throws

    - DomainError: an argument outside its physical domain (`dt` is checked at every
      temperature), or `sintered` with `phi0_prev` that is 1 in Float64 or
      `phi_d > phi0_prev`
"""
function dry_porosity_update(
    phi_d::Real,
    P_bar::Real,
    T::Real,
    dt::Real,
    pcfg::PorosityConfig;
    phi0_prev::Real=1.0,
    sintered::Bool=false,
)
    require_unit_interval(phi_d, "dry porosity phi_d")
    require_unit_interval(phi0_prev, "porosity before sintering phi0_prev")
    require_positive_finite(T, "temperature T")
    require_nonneg_finite(dt, "time step dt")
    if sintered
        Float64(phi0_prev) < 1.0 ||
            throw(DomainError(phi0_prev, "a sintered marker needs phi0_prev < 1"))
        phi_d <= phi0_prev ||
            throw(DomainError(phi_d, "a sintered marker needs phi_d <= phi0_prev"))
        require_nonneg_finite(P_bar, "pressure P_bar")
        phi_0 = Float64(phi0_prev)
    else
        phi_0 = min(
            Float64(phi0_prev), porosity_phi0(P_bar, pcfg.phi_init; phi_cap=pcfg.phi_cap)
        )
    end
    phi_p = min(Float64(phi_d), phi_0)
    T >= pcfg.T_sinter || return (phi=phi_p, phi0=phi_0, dphi_sinter=0.0, capped=false)
    phi_s, _, capped = sinter_porosity_step(
        phi_p,
        phi_0,
        P_bar,
        T,
        dt;
        R0=pcfg.grain_radius,
        A=pcfg.A,
        E_act=pcfg.E_act,
        n=pcfg.n,
    )
    phi_s = clamp(phi_s, 0.0, 1.0)
    return (phi=phi_s, phi0=phi_0, dphi_sinter=phi_p - phi_s, capped=capped)
end

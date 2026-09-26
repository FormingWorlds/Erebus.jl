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
        # Universal gas constant R [J/(mol K)]
        R_gas = 8.31446261815324
        log_ratio = (Ea / R_gas) * (inv(T) - inv(T0))
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
    compute_hydrofracture_factor(Peff::Real, sigma_t::Real;
                                 active::Bool = true,
                                 kappa_frac::Real = 1.0e3,
                                 gamma::Real = 1.0,
                                 max_factor::Real = Inf)

Compute dimensionless permeability enhancement factor from dynamic hydrofracturing.

When pore fluid pressure exceeds total confining pressure plus tensile strength
(Terzaghi effective pressure Peff = Pt - Pf <= -sigma_t), hydraulic tensile
fractures open and increase effective permeability:

    factor = 1.0 + kappa_frac * ((-Peff - sigma_t) / sigma_t)^gamma

clamped to [1.0, max_factor].
"""
function compute_hydrofracture_factor(
    Peff::Real,
    sigma_t::Real;
    active::Bool=true,
    kappa_frac::Real=1.0e3,
    gamma::Real=1.0,
    max_factor::Real=Inf,
)
    if !active || !isfinite(Peff) || !isfinite(sigma_t) || sigma_t <= 0.0
        return 1.0
    end
    overpressure = -Peff - sigma_t
    if overpressure <= 0.0
        return 1.0
    end
    norm_overpressure = overpressure / sigma_t
    factor = 1.0 + kappa_frac * (norm_overpressure ^ gamma)
    return clamp(Float64(factor), 1.0, Float64(max_factor))
end

"""
    compute_hydrofracture_permeability(kphi::Real, Peff::Real, sigma_t::Real;
                                       active::Bool = true,
                                       kappa_frac::Real = 1.0e3,
                                       gamma::Real = 1.0,
                                       kmax::Real = 1.0e-9)

Compute effective permeability k_eff [m²] with dynamic hydrofracturing enhancement.

When pore fluid pressure exceeds total confining pressure plus rock tensile strength:

    Peff = Pt - Pf <= -sigma_t

tensile microcracks open and enhance permeability according to:

    k_eff = min(kphi * compute_hydrofracture_factor(Peff, sigma_t; active, kappa_frac, gamma), kmax)

# Arguments
- `kphi`: baseline matrix permeability [m²]
- `Peff`: Terzaghi effective pressure Pt - Pf [Pa]
- `sigma_t`: rock tensile strength [Pa]
- `active`: enable hydrofracture enhancement (default: true)
- `kappa_frac`: dimensionless enhancement multiplier (default: 1.0e3)
- `gamma`: power-law exponent (default: 1.0)
- `kmax`: maximum permeability ceiling [m²] (default: 1.0e-9)

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
        Peff, sigma_t; active=active, kappa_frac=kappa_frac, gamma=gamma
    )
    k_enhanced = Float64(kphi * factor)
    kphi_f = Float64(kphi)
    kmax_f = Float64(kmax)
    return min(max(k_enhanced, kphi_f), max(kphi_f, kmax_f))
end

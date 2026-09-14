"""
Compute silicate partial melt fraction based on temperature, pressure, and rock type.

$(SIGNATURES)

# Details
- `T`: marker temperature [K]
- `P`: marker pressure [Pa] (fluid or lithostatic)
- `tm`: marker phase material type (1: core/mantle, 2: crust, 3: air)
- `T_sol`: solidus temperature [K]
- `T_liq`: liquidus temperature [K]
- `dpdt`: Clapeyron slope dT_sol/dP [K/Pa]

# Returns
- `F_m`: silicate melt mass/volume fraction in [0, 1]
"""
function compute_melt_fraction(
    T::Real, P::Real, tm::Integer; T_sol::Real=1400.0, T_liq::Real=1800.0, dpdt::Real=0.0
)
    if !isfinite(T) || T <= 0.0
        throw(DomainError(T, "Absolute temperature must be positive and finite"))
    end
    if !isfinite(P)
        throw(DomainError(P, "Pressure must be finite"))
    end
    if !isfinite(T_sol) || !isfinite(T_liq) || T_sol >= T_liq
        throw(
            DomainError(
                (T_sol, T_liq),
                "Solidus temperature must be finite and strictly less than liquidus temperature",
            ),
        )
    end
    if tm >= 3
        return 0.0
    end
    T_s = T_sol + dpdt * max(0.0, P)
    T_l = T_liq + dpdt * max(0.0, P)
    if T_s >= T_l
        throw(
            DomainError(
                (T_s, T_l),
                "Pressure-shifted solidus must be strictly less than shifted liquidus",
            ),
        )
    end
    if T <= T_s
        return 0.0
    elseif T >= T_l
        return 1.0
    else
        return (T - T_s) / (T_l - T_s)
    end
end

"""
Compute apparent volumetric heat capacity of silicate rock including latent heat of melting.

$(SIGNATURES)

# Details
- `T`: marker temperature [K]
- `P`: marker pressure [Pa]
- `rhocp_solid`: sensible solid volumetric heat capacity [J/(m^3 K)]
- `rho_solid`: solid density [kg/m^3]
- `tm`: marker phase material type (1: core/mantle, 2: crust, 3: air)
- `T_sol`: solidus temperature [K]
- `T_liq`: liquidus temperature [K]
- `L_melt`: latent heat of silicate melting [J/kg]
- `active`: boolean flag to enable latent heat addition
- `dpdt`: Clapeyron slope dT_sol/dP [K/Pa]

# Returns
- `rhocp_eff`: effective volumetric heat capacity [J/(m^3 K)]
"""
function rhocp_apparent_silicate(
    T::Real,
    P::Real,
    rhocp_solid::Real,
    rho_solid::Real,
    tm::Integer;
    T_sol::Real=1400.0,
    T_liq::Real=1800.0,
    L_melt::Real=4.0e5,
    active::Bool=true,
    dpdt::Real=0.0,
)
    if !isfinite(T) || T <= 0.0
        throw(DomainError(T, "Absolute temperature must be positive and finite"))
    end
    if !isfinite(P)
        throw(DomainError(P, "Pressure must be finite"))
    end
    if !isfinite(T_sol) || !isfinite(T_liq) || T_sol >= T_liq
        throw(
            DomainError(
                (T_sol, T_liq),
                "Solidus temperature must be finite and strictly less than liquidus temperature",
            ),
        )
    end
    if !active || tm >= 3
        return rhocp_solid
    end
    T_s = T_sol + dpdt * max(0.0, P)
    T_l = T_liq + dpdt * max(0.0, P)
    if T_s >= T_l
        throw(
            DomainError(
                (T_s, T_l),
                "Pressure-shifted solidus must be strictly less than shifted liquidus",
            ),
        )
    end
    if T_s < T < T_l
        dFdT = inv(T_l - T_s)
        return rhocp_solid + rho_solid * L_melt * dFdT
    else
        return rhocp_solid
    end
end

"""
Compute melt-weakened matrix viscosity and suspension transition.

$(SIGNATURES)

# Details
- `eta_solid`: solid rock matrix viscosity [Pa s]
- `F_m`: silicate melt fraction in [0, 1]
- `tm`: marker phase material type (1: core/mantle, 2: crust, 3: air)
- `alpha_eta`: melt weakening exponent (Costa et al. 2009; Gerya 2019 Section 16.6.2)
- `phi_crit`: rheologically critical melt fraction for disaggregation
- `eta_melt`: pure liquid melt viscosity limit [Pa s]
- `etamin`: lower viscosity clamp [Pa s]
- `etamax`: upper viscosity clamp [Pa s]

# Returns
- `eta_eff`: effective shear viscosity [Pa s]
"""
function compute_melt_weakened_viscosity(
    eta_solid::Real,
    F_m::Real,
    tm::Integer;
    alpha_eta::Real=28.0,
    phi_crit::Real=0.4,
    eta_melt::Real=10.0,
    etamin::Real=1.0e12,
    etamax::Real=1.0e23,
)
    if !isfinite(F_m)
        throw(DomainError(F_m, "Melt fraction must be finite"))
    end
    if !isfinite(eta_solid) || eta_solid <= 0.0
        throw(DomainError(eta_solid, "Solid viscosity must be positive and finite"))
    end
    if tm >= 3
        return clamp(eta_solid, etamin, etamax)
    end
    F_clamped = clamp(F_m, 0.0, 1.0)
    if iszero(F_clamped)
        return clamp(eta_solid, etamin, etamax)
    elseif F_clamped < phi_crit
        eta_weak = eta_solid * exp(-alpha_eta * F_clamped)
        return clamp(eta_weak, etamin, etamax)
    else
        # Critical disaggregation into magma ocean suspension
        eta_at_crit = eta_solid * exp(-alpha_eta * phi_crit)
        frac = (F_clamped - phi_crit) / (1.0 - phi_crit)
        log_eta = (1.0 - frac) * log(eta_at_crit) + frac * log(eta_melt)
        eta_susp = exp(log_eta)
        return clamp(eta_susp, etamin, etamax)
    end
end

"""
Compute silicate melt permeability through a compacting solid matrix.

$(SIGNATURES)

Calculates the permeable channel network conductivity according to the
McKenzie (1984) power-law formulation with residual melt retention threshold.

# Arguments
- `F_m::Real`: Silicate melt volume fraction in [0, 1]

# Keyword Arguments
- `k0::Real=1.0e-11`: Reference permeability [m^2]
- `phi0::Real=0.10`: Reference melt fraction [-]
- `n::Real=3.0`: Permeability power-law exponent [-]
- `phi_residual::Real=0.01`: Residual melt retention threshold [-]
- `phi_crit::Real=0.40`: Rheological critical melt fraction for matrix disaggregation [-]

# Returns
- `k_m::Float64`: Effective silicate melt permeability [m^2]
"""
function silicate_melt_permeability(
    F_m::Real;
    k0::Real=1.0e-11,
    phi0::Real=0.10,
    n::Real=3.0,
    phi_residual::Real=0.01,
    phi_crit::Real=0.40,
)
    if !isfinite(F_m) || F_m < 0.0 || F_m > 1.0
        throw(DomainError(F_m, "Silicate melt fraction must be finite and within [0, 1]"))
    end
    if !isfinite(k0) || k0 <= 0.0
        throw(DomainError(k0, "Reference permeability k0 must be positive and finite"))
    end
    if !isfinite(phi0) || phi0 <= 0.0 || phi0 > 1.0
        throw(DomainError(phi0, "Reference melt fraction phi0 must be within (0, 1]"))
    end
    if !isfinite(n) || n <= 0.0
        throw(DomainError(n, "Permeability exponent n must be positive and finite"))
    end
    if !isfinite(phi_residual) || phi_residual < 0.0 || phi_residual >= phi_crit
        throw(
            DomainError(
                phi_residual, "Residual threshold must satisfy 0 <= phi_residual < phi_crit"
            ),
        )
    end
    if !isfinite(phi_crit) || phi_crit <= 0.0 || phi_crit >= 1.0
        throw(
            DomainError(phi_crit, "Critical melt fraction phi_crit must be within (0, 1)")
        )
    end

    if F_m <= phi_residual
        return 0.0
    elseif F_m <= phi_crit
        phi_eff = (F_m - phi_residual) / phi0
        return Float64(k0 * (phi_eff^n))
    else
        phi_crit_eff = (phi_crit - phi_residual) / phi0
        k_crit = k0 * (phi_crit_eff^n)
        d_phi = (F_m - phi_crit) / phi0
        return Float64(k_crit + k0 * d_phi)
    end
end

"""
Compute buoyant silicate melt segregation velocity through solid matrix or crystal mush.

$(SIGNATURES)

Calculates the upward/outward relative migration velocity between buoyant liquid
silicate melt and solid silicate rock. Transitions smoothly from Darcy porous flow
in the percolation regime to Stokes settling in the crystal suspension regime.

# Arguments
- `F_m::Real`: Silicate melt volume fraction in [0, 1]
- `drho::Real`: Density contrast (rho_solid - rho_melt) [kg/m^3]
- `g_acc::Real`: Gravitational acceleration magnitude [m/s^2]
- `eta_melt::Real`: Liquid silicate melt dynamic viscosity [Pa s]

# Keyword Arguments
- `k_melt_ref::Real=1.0e-11`: Reference permeability [m^2]
- `phi0::Real=0.10`: Reference melt fraction [-]
- `perm_exponent::Real=3.0`: Permeability power-law exponent [-]
- `phi_residual::Real=0.01`: Residual melt retention threshold [-]
- `phi_crit::Real=0.40`: Rheological critical melt fraction [-]
- `r_grain::Real=1.0e-3`: Crystal grain radius for Stokes settling [m]
- `hindered_exponent::Real=2.0`: Richardson-Zaki hindrance exponent [-]
- `F_perc_end::Real=0.35`: Upper melt fraction of pure percolation regime [-]
- `F_settle_start::Real=0.45`: Lower melt fraction of pure suspension regime [-]

# Returns
- `v_seg::Float64`: Buoyant segregation velocity magnitude [m/s]
"""
function silicate_melt_segregation_velocity(
    F_m::Real,
    drho::Real,
    g_acc::Real,
    eta_melt::Real;
    k_melt_ref::Real=1.0e-11,
    phi0::Real=0.10,
    perm_exponent::Real=3.0,
    phi_residual::Real=0.01,
    phi_crit::Real=0.40,
    r_grain::Real=1.0e-3,
    hindered_exponent::Real=2.0,
    F_perc_end::Real=0.35,
    F_settle_start::Real=0.45,
)
    if !isfinite(F_m) || F_m < 0.0 || F_m > 1.0
        throw(DomainError(F_m, "Silicate melt fraction must be finite and within [0, 1]"))
    end
    if !isfinite(drho)
        throw(DomainError(drho, "Density contrast must be finite"))
    end
    if !isfinite(g_acc) || g_acc < 0.0
        throw(DomainError(g_acc, "Gravity acceleration must be non-negative and finite"))
    end
    if !isfinite(eta_melt) || eta_melt <= 0.0
        throw(DomainError(eta_melt, "Melt viscosity must be positive and finite"))
    end
    if !isfinite(r_grain) || r_grain <= 0.0
        throw(DomainError(r_grain, "Grain radius must be positive and finite"))
    end
    if !isfinite(hindered_exponent) || hindered_exponent < 0.0
        throw(DomainError(hindered_exponent, "Hindrance exponent must be non-negative"))
    end
    if !(0.0 < F_perc_end <= F_settle_start <= 1.0)
        throw(
            DomainError(
                (F_perc_end, F_settle_start),
                "Regime bounds must satisfy 0 < F_perc_end <= F_settle_start <= 1",
            ),
        )
    end

    if F_m <= phi_residual || drho <= 0.0 || g_acc <= 0.0
        return 0.0
    end

    # Darcy percolation velocity: v_perc = (k_m / (eta_m * F_m)) * drho * g
    k_m = silicate_melt_permeability(
        F_m;
        k0=k_melt_ref,
        phi0=phi0,
        n=perm_exponent,
        phi_residual=phi_residual,
        phi_crit=phi_crit,
    )
    v_perc = (k_m / (eta_melt * F_m)) * drho * g_acc

    # Stokes crystal suspension velocity: v_susp = (2 r^2 drho g / (9 eta_m)) * (F_m)^m
    v_stokes = (2.0 * r_grain^2 * drho * g_acc) / (9.0 * eta_melt)
    v_susp = v_stokes * (F_m^hindered_exponent)

    if F_m <= F_perc_end
        return Float64(v_perc)
    elseif F_m >= F_settle_start
        return Float64(v_susp)
    else
        # Smooth Hermite cubic interpolation between regimes
        xi = (F_m - F_perc_end) / (F_settle_start - F_perc_end)
        w = 3.0 * xi^2 - 2.0 * xi^3
        return Float64((1.0 - w) * v_perc + w * v_susp)
    end
end

"""
Compute gravitational potential energy dissipation heating during silicate melt ascent.

$(SIGNATURES)

Calculates the volumetric heat generation rate from buoyancy dissipation
as molten rock migrates relative to the solid mantle matrix.

# Arguments
- `F_m::Real`: Silicate melt volume fraction in [0, 1]
- `drho::Real`: Density contrast (rho_solid - rho_melt) [kg/m^3]
- `g_acc::Real`: Gravitational acceleration magnitude [m/s^2]
- `v_seg::Real`: Silicate melt segregation velocity magnitude [m/s]

# Returns
- `Q_diss::Float64`: Volumetric dissipation heat source [W/m^3]
"""
function silicate_melt_dissipation_heating(F_m::Real, drho::Real, g_acc::Real, v_seg::Real)
    if !isfinite(F_m) || F_m < 0.0 || F_m > 1.0
        throw(DomainError(F_m, "Silicate melt fraction must be finite and within [0, 1]"))
    end
    if !isfinite(drho) || !isfinite(g_acc) || !isfinite(v_seg)
        throw(DomainError((drho, g_acc, v_seg), "Inputs must be finite"))
    end
    if F_m <= 0.0 || drho <= 0.0 || g_acc <= 0.0 || v_seg <= 0.0
        return 0.0
    end
    return Float64(drho * g_acc * F_m * v_seg)
end

"""
Compute regularized effective thermal conductivity from soft turbulence.

$(SIGNATURES)

Blends conductive thermal conductivity and turbulent convective conductivity
smoothly in logarithmic space across a melt fraction transition window.

# Arguments
- `k_cond`: Conductive thermal conductivity [W/(m K)]
- `eta_num`: Numerical shear viscosity used in momentum solver [Pa s]
- `eta_fluid`: Physical fluid/magma viscosity [Pa s]
- `F_m`: Silicate melt fraction [0, 1]
- `T_marker`: Local marker temperature [K]
- `T_surface`: Reference surface/ambient temperature [K]

# Keyword Arguments
- `turb_exponent`: Power exponent for viscosity ratio (default 1/3 from Solomatov 2007)
- `F_start`: Lower boundary of transition window [0, 1] (default 0.30)
- `F_end`: Upper boundary of transition window [0, 1] (default 0.50)
- `dT_min`: Temperature contrast scale [K] (default 10.0)
- `k_floor`: Minimum thermal conductivity [W/(m K)] (default 1.0e-3)
- `k_cutoff`: Maximum thermal conductivity [W/(m K)] (default 1.0e6)

# Returns
- `k_eff`: Effective thermal conductivity [W/(m K)]
"""
function regularized_soft_turbulence_conductivity(
    k_cond::Real,
    eta_num::Real,
    eta_fluid::Real,
    F_m::Real,
    T_marker::Real,
    T_surface::Real;
    turb_exponent::Real=1.0 / 3.0,
    F_start::Real=0.30,
    F_end::Real=0.50,
    dT_min::Real=10.0,
    k_floor::Real=1.0e-3,
    k_cutoff::Real=1.0e6,
)
    if !isfinite(k_cond) ||
        k_cond <= 0.0 ||
        !isfinite(eta_num) ||
        eta_num <= 0.0 ||
        !isfinite(eta_fluid) ||
        eta_fluid <= 0.0
        throw(
            DomainError(
                (k_cond, eta_num, eta_fluid),
                "k_cond and viscosities must be positive and finite",
            ),
        )
    end
    if !isfinite(F_m) || !isfinite(T_marker) || !isfinite(T_surface)
        throw(
            DomainError(
                (F_m, T_marker, T_surface), "Melt fraction and temperatures must be finite"
            ),
        )
    end
    if F_end <= F_start
        throw(DomainError((F_start, F_end), "F_end must be strictly greater than F_start"))
    end
    if !isfinite(turb_exponent) || turb_exponent <= 0.0
        throw(DomainError(turb_exponent, "Turbulence exponent must be positive and finite"))
    end

    # Melt fraction weight factor using cubic smoothstep
    w_F = smoothstep(F_start, F_end, F_m)

    # Temperature contrast weight factor
    dT = max(0.0, T_marker - T_surface)
    w_T = clamp(dT / dT_min, 0.0, 1.0)^2

    w_total = w_F * w_T
    if w_total <= 0.0
        return max(k_cond, clamp(k_cond, k_floor, k_cutoff))
    end

    # Target turbulent conductivity from Solomatov (2007) scaling reduction
    k_turb_raw = k_cond * (eta_num / eta_fluid)^turb_exponent
    k_turb = max(k_cond, clamp(k_turb_raw, k_floor, k_cutoff))

    # Geometric blend in logarithmic space
    log_k = lerp(log10(k_cond), log10(k_turb), w_total)
    return max(k_cond, clamp(10.0^log_k, k_floor, k_cutoff))
end

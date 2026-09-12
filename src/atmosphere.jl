"""
Coupled 1D atmosphere, disk gas envelope, Guillot semi-grey radiation, and crossover escape physics.

This module models:
1. Gravitational gas envelope capture and Ormel et al. (2015) recycling limits.
2. Hydrodynamic boil-off mass loss during protoplanetary disk dispersal.
3. Multi-species atmospheric column optical depth and greenhouse blanketing.
4. Guillot (2010) semi-grey analytical radiative equilibrium profiles.
5. Greenhouse-attenuated surface heat transfer coefficients.
6. Zahnle & Kasting (1986) hydrodynamic crossover escape for multi-species outgassing.
"""

"""
Coupled 1D atmosphere and envelope dynamic state.

$(FIELDS)
"""
mutable struct AtmosphereState
    M_atm::Dict{Symbol,Float64}
    M_escaped::Dict{Symbol,Float64}
    P_surf::Float64
    T_surf_eq::Float64
    tau_LW::Float64
    M_env_bound::Float64
    F_net_rad::Float64
    h_rad_eff::Float64
end

function AtmosphereState(;
    M_atm::Dict{Symbol,Float64}=Dict{Symbol,Float64}(),
    M_escaped::Dict{Symbol,Float64}=Dict{Symbol,Float64}(),
    P_surf::Real=0.0,
    T_surf_eq::Real=0.0,
    tau_LW::Real=0.0,
    M_env_bound::Real=0.0,
    F_net_rad::Real=0.0,
    h_rad_eff::Real=0.0,
)
    return AtmosphereState(
        M_atm,
        M_escaped,
        Float64(P_surf),
        Float64(T_surf_eq),
        Float64(tau_LW),
        Float64(M_env_bound),
        Float64(F_net_rad),
        Float64(h_rad_eff),
    )
end

"""
    compute_gravitational_capture_radius(M::Real, M_star::Real, a::Real, c_s::Real)::Float64

Compute the gravitational capture radius of an embedded planetesimal, defined as the minimum of the Bondi radius and the Hill radius:
    R_cap = min(R_Bondi, R_Hill) = min(G * M / c_s^2, a * (M / (3 * M_star))^(1/3))

# Parameters
- `M`: Planetesimal mass [kg].
- `M_star`: Central stellar mass [kg].
- `a`: Semi-major axis / orbital separation [m].
- `c_s`: Disk gas sound speed [m/s].

# Returns
- `R_cap`: Gravitational capture radius [m].

# Raises
- `DomainError`: If any parameter is <= 0 or non-finite.
"""
function compute_gravitational_capture_radius(
    M::Real, M_star::Real, a::Real, c_s::Real
)::Float64
    M_val = Float64(M)
    M_star_val = Float64(M_star)
    a_val = Float64(a)
    c_s_val = Float64(c_s)

    if M_val <= 0.0 || !isfinite(M_val)
        throw(DomainError(M_val, "Planetesimal mass must be > 0 and finite"))
    end
    if M_star_val <= 0.0 || !isfinite(M_star_val)
        throw(DomainError(M_star_val, "Stellar mass must be > 0 and finite"))
    end
    if a_val <= 0.0 || !isfinite(a_val)
        throw(DomainError(a_val, "Semi-major axis must be > 0 and finite"))
    end
    if c_s_val <= 0.0 || !isfinite(c_s_val)
        throw(DomainError(c_s_val, "Sound speed must be > 0 and finite"))
    end

    R_Bondi = GRAVITATIONAL_CONSTANT * M_val / (c_s_val^2)
    R_Hill = a_val * cbrt(M_val / (3.0 * M_star_val))
    return min(R_Bondi, R_Hill)
end

"""
    compute_disk_envelope_mass(
        M::Real, R_planet::Real, R_cap::Real, rho_disk::Real, c_s::Real;
        f_rec::Real=0.10,
    )::Float64

Compute the bound gas envelope mass within capture radius R_cap under isothermal hydrostatic equilibrium, capped by the Ormel et al. (2015) recycling limit:
    M_env = min(M_iso, f_rec * (4π/3) * R_cap^3 * rho_disk)
where M_iso = 4π ∫_{R_planet}^{R_cap} r^2 rho_disk exp[(G*M/c_s^2)*(1/r - 1/R_cap)] dr.

# Parameters
- `M`: Planetesimal mass [kg].
- `R_planet`: Planetesimal radius [m].
- `R_cap`: Gravitational capture radius [m].
- `rho_disk`: Protoplanetary disk gas density [kg/m^3].
- `c_s`: Disk sound speed [m/s].

# Keywords
- `f_rec`: Recycling fraction cap (default: 0.10, Ormel et al. 2015).

# Returns
- `M_env`: Bound envelope mass [kg].

# Raises
- `DomainError`: If inputs are negative, zero, or non-finite.
"""
function compute_disk_envelope_mass(
    M::Real, R_planet::Real, R_cap::Real, rho_disk::Real, c_s::Real; f_rec::Real=0.10
)::Float64
    M_val = Float64(M)
    R_p = Float64(R_planet)
    R_c = Float64(R_cap)
    rho = Float64(rho_disk)
    cs = Float64(c_s)
    f_rec_val = Float64(f_rec)

    if M_val < 0.0 || !isfinite(M_val)
        throw(DomainError(M_val, "Planetesimal mass must be >= 0 and finite"))
    end
    if R_p < 0.0 || !isfinite(R_p)
        throw(DomainError(R_p, "Planetesimal radius must be >= 0 and finite"))
    end
    if R_c < 0.0 || !isfinite(R_c)
        throw(DomainError(R_c, "Capture radius must be >= 0 and finite"))
    end
    if rho < 0.0 || !isfinite(rho)
        throw(DomainError(rho, "Disk gas density must be >= 0 and finite"))
    end
    if cs <= 0.0 || !isfinite(cs)
        throw(DomainError(cs, "Sound speed must be > 0 and finite"))
    end
    if f_rec_val < 0.0 || !isfinite(f_rec_val)
        throw(DomainError(f_rec_val, "Recycling factor must be >= 0 and finite"))
    end

    if rho == 0.0 || R_p >= R_c || M_val == 0.0
        return 0.0
    end

    # 64-panel Simpson integration
    N_panels = 64
    dr = (R_c - R_p) / N_panels
    GM_cs2 = GRAVITATIONAL_CONSTANT * M_val / (cs^2)
    inv_Rc = 1.0 / R_c

    integral = 0.0
    for i in 0:N_panels
        r = R_p + i * dr
        w = if i == 0 || i == N_panels
            1.0 / 3.0
        elseif isodd(i)
            4.0 / 3.0
        else
            2.0 / 3.0
        end
        arg = clamp(GM_cs2 * (1.0 / r - inv_Rc), 0.0, 50.0)
        rho_r = rho * exp(arg)
        integral += w * (r^2) * rho_r * dr
    end

    M_iso = 4.0 * π * integral
    M_rec_cap = f_rec_val * (4.0 * π / 3.0) * (R_c^3) * rho
    return min(M_iso, M_rec_cap)
end

"""
    compute_atmospheric_optical_depth(
        M_atm::AbstractDict{Symbol,<:Real},
        R_planet::Real,
        opacities::AbstractDict{Symbol,<:Real};
        kappa_default::Real=1.0e-2,
    )::Float64

Compute the total longwave optical depth of an outgassed atmosphere over surface area 4π R_planet^2:
    τ_LW = (1 / 4π R_planet^2) * ∑_i κ_i M_{atm, i}

# Parameters
- `M_atm`: Dictionary of species atmospheric masses [kg].
- `R_planet`: Planetary surface radius [m].
- `opacities`: Dictionary of species mass absorption coefficients / opacities [m^2/kg].

# Keywords
- `kappa_default`: Fallback opacity [m^2/kg] for unlisted species (default: 1.0e-2).

# Returns
- `tau_LW`: Total infrared optical depth (dimensionless).

# Raises
- `DomainError`: If R_planet <= 0, any mass < 0, or any opacity < 0.
"""
function compute_atmospheric_optical_depth(
    M_atm::AbstractDict{Symbol,<:Real},
    R_planet::Real,
    opacities::AbstractDict{Symbol,<:Real};
    kappa_default::Real=1.0e-2,
)::Float64
    R_p = Float64(R_planet)
    if R_p <= 0.0 || !isfinite(R_p)
        throw(DomainError(R_p, "Planetary radius must be > 0 and finite"))
    end
    k_def = Float64(kappa_default)
    if k_def < 0.0 || !isfinite(k_def)
        throw(DomainError(k_def, "Default opacity must be >= 0 and finite"))
    end

    area = 4.0 * π * (R_p^2)
    sum_kappa_mass = 0.0
    for (sp, m) in M_atm
        m_val = Float64(m)
        if m_val < 0.0 || !isfinite(m_val)
            throw(
                DomainError(
                    m_val, "Atmospheric species mass for $sp must be >= 0 and finite"
                ),
            )
        end
        kap = Float64(get(opacities, sp, k_def))
        if kap < 0.0 || !isfinite(kap)
            throw(DomainError(kap, "Opacity for species $sp must be >= 0 and finite"))
        end
        sum_kappa_mass += kap * m_val
    end

    return sum_kappa_mass / area
end

"""
    compute_guillot_surface_temperature(
        tau_LW::Real, T_int::Real, T_irr::Real;
        T_eqm::Union{Real,Nothing}=nothing,
        gamma::Real=0.10, albedo::Real=0.20,
    )::Float64

Evaluate analytical surface temperature under semi-grey radiative equilibrium (Guillot 2010, Eq. 49):
    T_eqm^4 = (1 - albedo) * T_irr^4 / 4
    T^4(τ) = (3 * T_int^4 / 4) * (2/3 + τ) +
             (3 * T_eqm^4 / 4) * { 2/3 + 1 / (γ * √3) + (γ / √3 - 1 / (γ * √3)) * exp(-γ * τ * √3) }

# Parameters
- `tau_LW`: Infrared optical depth τ at the surface.
- `T_int`: Planetary internal effective temperature [K] (from interior heat flux F_int = σ T_int^4).
- `T_irr`: Irradiation / ambient stellar temperature [K].

# Keywords
- `T_eqm`: Planetary equilibrium temperature [K]. If supplied, T_eqm^4 is scaled by (1 - albedo).
- `gamma`: Ratio of visible/shortwave opacity to thermal/longwave opacity κ_vis / κ_th (default: 0.10).
- `albedo`: Bond albedo (default: 0.20).

# Returns
- `T_surf`: Radiative equilibrium surface temperature [K].

# Raises
- `DomainError`: If tau_LW < 0, T_int < 0, T_irr < 0, gamma <= 0, or albedo not in [0, 1).
"""
function compute_guillot_surface_temperature(
    tau_LW::Real,
    T_int::Real,
    T_irr::Real;
    T_eqm::Union{Real,Nothing}=nothing,
    gamma::Real=0.10,
    albedo::Real=0.20,
)::Float64
    tau = Float64(tau_LW)
    Tint = Float64(T_int)
    Tirr = Float64(T_irr)
    gam = Float64(gamma)
    alb = Float64(albedo)

    if tau < 0.0 || !isfinite(tau)
        throw(DomainError(tau, "Optical depth must be >= 0 and finite"))
    end
    if Tint < 0.0 || !isfinite(Tint)
        throw(DomainError(Tint, "Internal temperature must be >= 0 and finite"))
    end
    if Tirr < 0.0 || !isfinite(Tirr)
        throw(DomainError(Tirr, "Irradiation temperature must be >= 0 and finite"))
    end
    if gam <= 0.0 || !isfinite(gam)
        throw(DomainError(gam, "Opacity ratio gamma must be > 0 and finite"))
    end
    if alb < 0.0 || alb >= 1.0 || !isfinite(alb)
        throw(DomainError(alb, "Albedo must be in [0, 1) and finite"))
    end

    Teqm4 = if T_eqm !== nothing
        T_e = Float64(T_eqm)
        if T_e < 0.0 || !isfinite(T_e)
            throw(DomainError(T_e, "Equilibrium temperature must be >= 0 and finite"))
        end
        T_e^4
    else
        0.25 * (1.0 - alb) * (Tirr^4)
    end
    term1 = 0.75 * (Tint^4) * (2.0 / 3.0 + tau)

    sqrt3 = sqrt(3.0)
    inv_gam_sqrt3 = 1.0 / (gam * sqrt3)
    gam_over_sqrt3 = gam / sqrt3
    exp_term = exp(-gam * tau * sqrt3)
    bracket = 2.0 / 3.0 + inv_gam_sqrt3 + (gam_over_sqrt3 - inv_gam_sqrt3) * exp_term

    term2 = 0.75 * Teqm4 * bracket
    T4 = max(1.0, term1 + term2)
    return (T4)^0.25
end

"""
    compute_effective_radiation_htc(
        T_surf::Real, T_amb::Real, tau_LW::Real;
        emissivity::Real=0.9, sigma_sb::Real=5.670374419e-8,
    )::Float64

Compute the greenhouse-attenuated radiative heat transfer coefficient across the planetary surface:
    h_rad,eff = compute_radiation_htc(T_surf, T_amb; emissivity, sigma_sb) / (1 + 0.75 * tau_LW)

# Parameters
- `T_surf`: Surface temperature [K].
- `T_amb`: Ambient / skin temperature [K].
- `tau_LW`: Infrared optical depth (dimensionless).

# Keywords
- `emissivity`: Surface emissivity (default: 0.9).
- `sigma_sb`: Stefan-Boltzmann constant (default: 5.670374419e-8 W/(m^2 K^4)).

# Returns
- `h_rad_eff`: Effective heat transfer coefficient [W/(m^2 K)].

# Raises
- `DomainError`: If tau_LW < 0, temperatures are non-positive, or emissivity not in (0, 1].
"""
function compute_effective_radiation_htc(
    T_surf::Real,
    T_amb::Real,
    tau_LW::Real;
    emissivity::Real=0.9,
    sigma_sb::Real=5.670374419e-8,
)::Float64
    T_s = Float64(T_surf)
    T_a = Float64(T_amb)
    tau = Float64(tau_LW)

    if !isfinite(T_s) || T_s <= 0.0
        throw(DomainError(T_s, "Surface temperature must be > 0 and finite"))
    end
    if !isfinite(T_a) || T_a <= 0.0
        throw(DomainError(T_a, "Ambient temperature must be > 0 and finite"))
    end
    if tau < 0.0 || !isfinite(tau)
        throw(DomainError(tau, "Optical depth must be >= 0 and finite"))
    end
    if !(0.0 <= emissivity <= 1.0)
        throw(DomainError(emissivity, "Emissivity must be in [0.0, 1.0]"))
    end
    h_bare = compute_radiation_htc(T_s, T_a; emissivity=emissivity, sigma_sb=sigma_sb)
    return h_bare / (1.0 + 0.75 * tau)
end

"""
    compute_boiloff_rate(M_env::Real, M_env_target::Real, tau_boil::Real)::Float64

Compute the hydrodynamic boil-off mass loss rate of an unbound gas envelope:
    dM_boil / dt = max(0, M_env - M_env_target) / tau_boil

# Parameters
- `M_env`: Current envelope mass [kg].
- `M_env_target`: Equilibrium target envelope mass [kg].
- `tau_boil`: Boil-off relaxation timescale [s].

# Returns
- `dM_dt`: Hydrodynamic boil-off mass loss rate [kg/s].

# Raises
- `DomainError`: If masses are negative, tau_boil <= 0, or inputs are non-finite.
"""
function compute_boiloff_rate(M_env::Real, M_env_target::Real, tau_boil::Real)::Float64
    M_e = Float64(M_env)
    M_t = Float64(M_env_target)
    tau = Float64(tau_boil)

    if M_e < 0.0 || !isfinite(M_e)
        throw(DomainError(M_e, "Envelope mass must be >= 0 and finite"))
    end
    if M_t < 0.0 || !isfinite(M_t)
        throw(DomainError(M_t, "Target envelope mass must be >= 0 and finite"))
    end
    if tau <= 0.0 || !isfinite(tau)
        throw(DomainError(tau, "Boil-off timescale must be > 0 and finite"))
    end

    excess = max(0.0, M_e - M_t)
    return excess / tau
end

"""
    compute_crossover_mass(
        m_carrier::Real, T_exo::Real, Phi_carrier::Real, g::Real, X_carrier::Real;
        b_diff::Real=1.0e21,
    )::Float64

Compute the Zahnle & Kasting (1986) hydrodynamic crossover mass m_c for species dragged by an escaping light carrier:
    m_c = m_carrier + (k_B * T_exo * Phi_carrier) / (b_diff * g * X_carrier)

# Parameters
- `m_carrier`: Molecular mass of escaping carrier gas (e.g. H2) [kg].
- `T_exo`: Exobase / upper atmosphere temperature [K].
- `Phi_carrier`: Escape flux of carrier gas [molecules / (m^2 s)].
- `g`: Local gravitational acceleration [m/s^2].
- `X_carrier`: Mole fraction of carrier gas in escaping flow.

# Keywords
- `b_diff`: Binary diffusion parameter [m^-1 s^-1] (default: 1.0e21).

# Returns
- `m_c`: Crossover mass [kg]. Species with molecular mass m_j >= m_c cannot escape hydrodynamically.

# Raises
- `DomainError`: If inputs are non-positive or non-finite.
"""
function compute_crossover_mass(
    m_carrier::Real,
    T_exo::Real,
    Phi_carrier::Real,
    g::Real,
    X_carrier::Real;
    b_diff::Real=1.0e21,
)::Float64
    m_car = Float64(m_carrier)
    T = Float64(T_exo)
    Phi = Float64(Phi_carrier)
    grav = Float64(g)
    X = Float64(X_carrier)
    b = Float64(b_diff)

    if m_car <= 0.0 || !isfinite(m_car)
        throw(DomainError(m_car, "Carrier mass must be > 0 and finite"))
    end
    if T <= 0.0 || !isfinite(T)
        throw(DomainError(T, "Exobase temperature must be > 0 and finite"))
    end
    if Phi < 0.0 || !isfinite(Phi)
        throw(DomainError(Phi, "Carrier flux must be >= 0 and finite"))
    end
    if grav <= 0.0 || !isfinite(grav)
        throw(DomainError(grav, "Surface gravity must be > 0 and finite"))
    end
    if X <= 0.0 || !isfinite(X)
        throw(DomainError(X, "Carrier mole fraction must be > 0 and finite"))
    end
    if b <= 0.0 || !isfinite(b)
        throw(DomainError(b, "Binary diffusion coefficient must be > 0 and finite"))
    end

    if Phi == 0.0
        return m_car
    end

    drag_term = (K_BOLTZMANN * T * Phi) / (b * grav * X)
    return m_car + drag_term
end

"""
    compute_crossover_drag_fraction(m_species::Real, m_c::Real, m_carrier::Real)::Float64

Compute the hydrodynamic drag efficiency factor x_j for a heavier species dragged by escaping hydrogen:
    x_j = 1.0 - (m_species - m_carrier) / (m_c - m_carrier)   if m_carrier <= m_species < m_c
    x_j = 1.0                                                  if m_species <= m_carrier
    x_j = 0.0                                                  if m_species >= m_c

# Parameters
- `m_species`: Molecular mass of dragged species [kg].
- `m_c`: Crossover mass [kg].
- `m_carrier`: Molecular mass of escaping carrier species [kg].

# Returns
- `x_j`: Drag efficiency factor in [0, 1].

# Raises
- `DomainError`: If molecular masses are non-positive or m_c < m_carrier.
"""
function compute_crossover_drag_fraction(
    m_species::Real, m_c::Real, m_carrier::Real
)::Float64
    m_sp = Float64(m_species)
    mc = Float64(m_c)
    m_car = Float64(m_carrier)

    if m_sp <= 0.0 || !isfinite(m_sp)
        throw(DomainError(m_sp, "Species mass must be > 0 and finite"))
    end
    if m_car <= 0.0 || !isfinite(m_car)
        throw(DomainError(m_car, "Carrier mass must be > 0 and finite"))
    end
    if mc < m_car || !isfinite(mc)
        throw(DomainError(mc, "Crossover mass must be >= carrier mass and finite"))
    end

    if m_sp <= m_car
        return 1.0
    elseif m_sp >= mc || mc == m_car
        return 0.0
    else
        return 1.0 - (m_sp - m_car) / (mc - m_car)
    end
end

# -----------------------------------------------------------------------------
# Multispecies Escape Closure and Binary Diffusion Parameters
# -----------------------------------------------------------------------------

"""
Standard atomic and molecular weights for atmospheric escape species [amu].
References `SPECIES_AMU` from the physics module as the single source of truth.
"""
const SPECIES_AMU_ESCAPE = SPECIES_AMU

"""
Kinetic diameters for atmospheric escape species [pm].
"""
const SPECIES_DIAMETER_PM = Dict{Symbol,Float64}(
    :H => 260.0,
    :D => 265.0,
    :He => 260.0,
    :Ne => 275.0,
    :Ar => 340.0,
    :Kr => 360.0,
    :Xe => 396.0,
    :O => 275.0,
    :C => 307.5658,
    :N => 280.4276,
    :S => 325.6579,
    :Na => 410.6908,
    :Mg => 312.9934,
    :Si => 379.9342,
    :Fe => 441.4474,
    :H2 => 289.0,
    :H2O => 265.0,
    :CO2 => 330.0,
    :O2 => 346.0,
    :CH4 => 380.0,
    :N2 => 364.0,
    :CO => 376.0,
    :NH3 => 260.0,
    :H2S => 360.0,
    :SO2 => 360.0,
    :S2 => 400.0,
)

"""
Atomic binary diffusion parameters at T = 1000 K [m^-1 s^-1].
Tabulated from Attia & Lichtenberg (2026, Data S1).
"""
const ATOMIC_BINARY_DIFFUSION_1000K_SI = Dict{Tuple{Symbol,Symbol},Float64}(
    (:H, :He) => 1.6e22,
    (:H, :N) => 1.059e22,
    (:H, :C) => 9.655e21,
    (:H, :Ne) => 9.3e21,
    (:H, :Mg) => 9.2856e21,
    (:H, :O) => 9.0e21,
    (:H, :S) => 8.8454e21,
    (:H, :Kr) => 7.818e21,
    (:H, :Si) => 7.4245e21,
    (:H, :Xe) => 6.9685e21,
    (:H, :Na) => 6.7851e21,
    (:H, :Ar) => 6.5e21,
    (:He, :Ne) => 6.3014e21,
    (:H, :Fe) => 6.126e21,
    (:He, :O) => 6.0e21,
    (:He, :N) => 5.9149e21,
    (:He, :C) => 5.4609e21,
    (:He, :Mg) => 5.0078e21,
    (:He, :S) => 4.7109e21,
    (:He, :Ar) => 4.4906e21,
    (:He, :Si) => 3.9765e21,
    (:He, :Kr) => 3.8388e21,
    (:He, :Na) => 3.6699e21,
    (:He, :Fe) => 3.2054e21,
    (:He, :Xe) => 3.1869e21,
    (:N, :O) => 3.1296e21,
    (:O, :Ne) => 3.0e21,
    (:C, :O) => 2.9681e21,
    (:C, :N) => 2.5189e21,
    (:N, :Ne) => 2.4832e21,
    (:O, :Mg) => 2.4569e21,
    (:C, :Ne) => 2.3649e21,
    (:O, :S) => 2.2385e21,
    (:N, :Mg) => 2.1096e21,
    (:C, :Mg) => 1.9757e21,
    (:O, :Si) => 1.9267e21,
    (:Ne, :Mg) => 1.92e21,
    (:N, :S) => 1.895e21,
    (:O, :Na) => 1.8271e21,
    (:C, :S) => 1.82e21,
    (:O, :Ar) => 1.8e21,
    (:N, :Ar) => 1.7907e21,
    (:O, :Kr) => 1.7853e21,
    (:C, :Ar) => 1.7393e21,
    (:Ne, :S) => 1.7365e21,
    (:Ne, :Ar) => 1.6161e21,
    (:O, :Xe) => 1.5518e21,
    (:N, :Si) => 1.509e21,
    (:Ne, :Si) => 1.4998e21,
    (:N, :Kr) => 1.4703e21,
    (:C, :Si) => 1.4674e21,
    (:O, :Fe) => 1.4576e21,
    (:Ne, :Na) => 1.43e21,
    (:N, :Na) => 1.4276e21,
    (:C, :Kr) => 1.3927e21,
    (:C, :Na) => 1.3883e21,
    (:Mg, :S) => 1.3776e21,
    (:Ne, :Kr) => 1.3571e21,
    (:Mg, :Ar) => 1.3477e21,
    (:N, :Xe) => 1.2358e21,
    (:C, :Xe) => 1.2251e21,
    (:S, :Ar) => 1.1997e21,
    (:Ne, :Xe) => 1.1723e21,
    (:Mg, :Si) => 1.1607e21,
    (:N, :Fe) => 1.1536e21,
    (:C, :Fe) => 1.1405e21,
    (:Na, :Mg) => 1.1175e21,
    (:Ne, :Fe) => 1.1163e21,
    (:Si, :Ar) => 1.0633e21,
    (:Si, :S) => 1.0443e21,
    (:Na, :Ar) => 1.037e21,
    (:Mg, :Kr) => 1.0233e21,
    (:Na, :S) => 1.014e21,
    (:Na, :Si) => 9.0514e20,
    (:Ar, :Kr) => 8.9036e20,
    (:S, :Kr) => 8.886e20,
    (:Mg, :Xe) => 8.837e20,
    (:Mg, :Fe) => 8.5885e20,
    (:Si, :Kr) => 8.0112e20,
    (:Na, :Kr) => 7.9739e20,
    (:Ar, :Fe) => 7.6625e20,
    (:S, :Xe) => 7.6095e20,
    (:Ar, :Xe) => 7.6038e20,
    (:S, :Fe) => 7.575e20,
    (:Na, :Xe) => 6.989e20,
    (:Si, :Xe) => 6.9465e20,
    (:Si, :Fe) => 6.8976e20,
    (:Na, :Fe) => 6.8649e20,
    (:Fe, :Kr) => 5.4101e20,
    (:Kr, :Xe) => 4.9214e20,
    (:Fe, :Xe) => 4.5826e20,
)

"""
Molecular binary diffusion parameters at T = 1000 K [m^-1 s^-1].
Tabulated from Zahnle & Kasting (2023) and Marrero & Mason (1972).
"""
const MOLECULAR_BINARY_DIFFUSION_1000K_SI = Dict{Tuple{Symbol,Symbol},Float64}(
    (:H2O, :CO2) => 1.56e21,
    (:H2O, :O2) => 1.59e21,
    (:CO2, :O2) => 1.00e21,
    (:CO2, :He) => 3.56e21,
    (:CO2, :Ne) => 1.62e21,
    (:CO2, :Ar) => 1.00e21,
    (:CO2, :N2) => 1.04e21,
    (:H2, :CO2) => 4.09e21,
    (:H2, :H2O) => 4.80e21,
    (:H2, :N2) => 4.71e21,
    (:H2, :He) => 1.30e23,
    (:H2, :Ne) => 9.40e22,
    (:H2, :Ar) => 7.10e22,
    (:H2, :Kr) => 6.10e22,
    (:H2, :Xe) => 5.00e22,
)

"""
    _canonical_escape_symbol(s::Symbol)::Symbol

Normalize volatile species symbol `s` to canonical casing used in diffusion dictionaries.
"""
function _canonical_escape_symbol(s::Symbol)::Symbol
    u = uppercase(String(s))
    if u == "H"
        return :H
    elseif u == "D"
        return :D
    elseif u == "HE"
        return :He
    elseif u == "C"
        return :C
    elseif u == "N"
        return :N
    elseif u == "O"
        return :O
    elseif u == "NE"
        return :Ne
    elseif u == "NA"
        return :Na
    elseif u == "MG"
        return :Mg
    elseif u == "SI"
        return :Si
    elseif u == "S"
        return :S
    elseif u == "AR"
        return :Ar
    elseif u == "FE"
        return :Fe
    elseif u == "KR"
        return :Kr
    elseif u == "XE"
        return :Xe
    elseif u == "H2"
        return :H2
    elseif u == "H2O"
        return :H2O
    elseif u == "CO"
        return :CO
    elseif u == "CO2"
        return :CO2
    elseif u == "CH4"
        return :CH4
    elseif u == "N2"
        return :N2
    elseif u == "NH3"
        return :NH3
    elseif u == "O2"
        return :O2
    elseif u == "H2S"
        return :H2S
    elseif u == "SO2"
        return :SO2
    elseif u == "S2"
        return :S2
    else
        return s
    end
end

"""
    get_binary_diffusion_parameter(
        species_1::Symbol, species_2::Symbol, T_K::Real;
        b_anchor::Real=4.09e21,
    )::Float64

Evaluate the binary diffusion parameter b_ij(T) = n D_ij in m^-1 s^-1 for a gas pair.
Directly accesses tabulated values or applies Zahnle & Kasting (2023) Eq. (10) scaling.

# Parameters
- `species_1`: Symbol of first gas species.
- `species_2`: Symbol of second gas species.
- `T_K`: Temperature [K].

# Keywords
- `b_anchor`: Reference anchor parameter [m^-1 s^-1] (default: 4.09e21 for H2-CO2).

# Returns
- `b_ij`: Binary diffusion parameter [m^-1 s^-1].

# Raises
- `DomainError`: If temperature is non-positive or non-finite.
"""
function get_binary_diffusion_parameter(
    species_1::Symbol, species_2::Symbol, T_K::Real; b_anchor::Real=4.09e21
)::Float64
    T = Float64(T_K)
    if T <= 0.0 || !isfinite(T)
        throw(DomainError(T, "Temperature must be > 0 and finite"))
    end
    s1 = _canonical_escape_symbol(species_1)
    s2 = _canonical_escape_symbol(species_2)
    if s1 === s2
        return Inf
    end

    pair_key = (s1, s2)
    rev_key = (s2, s1)

    b_1000 = if haskey(ATOMIC_BINARY_DIFFUSION_1000K_SI, pair_key)
        ATOMIC_BINARY_DIFFUSION_1000K_SI[pair_key]
    elseif haskey(ATOMIC_BINARY_DIFFUSION_1000K_SI, rev_key)
        ATOMIC_BINARY_DIFFUSION_1000K_SI[rev_key]
    elseif haskey(MOLECULAR_BINARY_DIFFUSION_1000K_SI, pair_key)
        MOLECULAR_BINARY_DIFFUSION_1000K_SI[pair_key]
    elseif haskey(MOLECULAR_BINARY_DIFFUSION_1000K_SI, rev_key)
        MOLECULAR_BINARY_DIFFUSION_1000K_SI[rev_key]
    else
        m1 = get(SPECIES_AMU_ESCAPE, s1, nothing)
        m2 = get(SPECIES_AMU_ESCAPE, s2, nothing)
        d1 = get(SPECIES_DIAMETER_PM, s1, nothing)
        d2 = get(SPECIES_DIAMETER_PM, s2, nothing)
        if m1 === nothing || m2 === nothing || d1 === nothing || d2 === nothing
            return Float64(b_anchor) * (T / 1000.0)^0.75
        end
        m_a1 = SPECIES_AMU_ESCAPE[:H2]
        m_a2 = SPECIES_AMU_ESCAPE[:CO2]
        d_a1 = SPECIES_DIAMETER_PM[:H2]
        d_a2 = SPECIES_DIAMETER_PM[:CO2]
        rm_target = 1.0 / m1 + 1.0 / m2
        rm_anchor = 1.0 / m_a1 + 1.0 / m_a2
        d_sum_target = d1 + d2
        d_sum_anchor = d_a1 + d_a2
        Float64(b_anchor) * sqrt(rm_target / rm_anchor) * ((d_sum_anchor / d_sum_target)^2)
    end

    return b_1000 * (T / 1000.0)^0.75
end

"""
    assemble_binary_diffusion_matrix(
        species_list::AbstractVector{Symbol}, T_K::Real
    )::Matrix{Float64}

Assemble the symmetric N x N binary diffusion matrix b_ij in m^-1 s^-1 with Inf on diagonal.

# Parameters
- `species_list`: Vector of species symbols.
- `T_K`: Temperature [K].

# Returns
- `b_mat`: Symmetric N x N matrix [m^-1 s^-1].
"""
function assemble_binary_diffusion_matrix(
    species_list::AbstractVector{Symbol}, T_K::Real
)::Matrix{Float64}
    N = length(species_list)
    b_mat = fill(Inf, N, N)
    for j in 1:N
        for i in (j + 1):N
            bij = get_binary_diffusion_parameter(species_list[i], species_list[j], T_K)
            b_mat[i, j] = bij
            b_mat[j, i] = bij
        end
    end
    return b_mat
end

"""
    solve_fixed_active(
        phi::Real, X::AbstractVector{<:Real}, m::AbstractVector{<:Real},
        T::Real, g0::Real, b::AbstractMatrix{<:Real},
        active::Union{AbstractVector{Int},Set{Int}},
    )::Tuple{Vector{Float64},Float64}

Solve the multispecies linear closure equations on a fixed candidate escaping active set.
Uses two-sided diagonal matrix equilibration to ensure numerical stability across large dynamic ranges.

# Parameters
- `phi`: Total mass escape flux [kg / (m^2 s)].
- `X`: Base mole fractions summing to 1.
- `m`: Particle masses [kg].
- `T`: Exobase temperature [K].
- `g0`: Gravitational acceleration [m/s^2].
- `b`: Symmetric binary diffusion parameter matrix [m^-1 s^-1].
- `active`: Collection of 1-based indices of escaping species.

# Returns
- `(w, C)`: Drift variables w_j [m^-2 s^-1] and inverse scale height C [m^-1].
"""
function solve_fixed_active(
    phi::Real,
    X::AbstractVector{<:Real},
    m::AbstractVector{<:Real},
    T::Real,
    g0::Real,
    b::AbstractMatrix{<:Real},
    active::Union{AbstractVector{Int},Set{Int}},
)::Tuple{Vector{Float64},Float64}
    kT = K_BOLTZMANN * Float64(T)
    grav = Float64(g0)
    phi_val = Float64(phi)
    A = sort(collect(active))
    N = length(X)
    R = [k for k in 1:N if k ∉ active]
    nA = length(A)

    if nA == 1
        j = A[1]
        w = zeros(Float64, N)
        w[j] = phi_val / (Float64(m[j]) * Float64(X[j]))
        R_sum = isempty(R) ? 0.0 : sum(Float64(X[k]) / Float64(b[j, k]) for k in R)
        C = Float64(m[j]) * grav / kT + w[j] * R_sum
        return w, C
    end

    M = zeros(Float64, nA + 1, nA + 1)
    rhs = zeros(Float64, nA + 1)
    for (row, j) in enumerate(A)
        diag = 0.0
        for (col, i) in enumerate(A)
            if i != j
                term = Float64(X[i]) / Float64(b[i, j])
                M[row, col] += term
                diag += term
            end
        end
        for k in R
            diag += Float64(X[k]) / Float64(b[j, k])
        end
        M[row, row] -= diag
        M[row, nA + 1] = 1.0
        rhs[row] = Float64(m[j]) * grav / kT
    end
    for (col, j) in enumerate(A)
        M[nA + 1, col] = Float64(m[j]) * Float64(X[j])
    end
    rhs[nA + 1] = phi_val

    sol = M \ rhs
    w = zeros(Float64, N)
    for (col, j) in enumerate(A)
        w[j] = sol[col]
    end
    return w, sol[nA + 1]
end

"""
    solve_multispecies_escape_closure(
        phi::Real, X::AbstractVector{<:Real}, m::AbstractVector{<:Real},
        T::Real, g0::Real, b::AbstractMatrix{<:Real};
        return_diag::Bool=false,
    )

Solve the general convex multispecies hydrodynamic escape closure (Attia & Lichtenberg 2026).
Partitions total mass flux phi into non-negative individual species number fluxes Phi_j.

# Parameters
- `phi`: Total mass escape flux [kg / (m^2 s)].
- `X`: Species base mole fractions.
- `m`: Species molecular masses [kg].
- `T`: Exobase temperature [K].
- `g0`: Exobase gravitational acceleration [m/s^2].
- `b`: Symmetric binary diffusion parameter matrix [m^-1 s^-1].

# Keywords
- `return_diag`: If true, returns `(Phi, C, active_set)`.

# Returns
- `Phi`: Vector of species number escape fluxes [molecules / (m^2 s)].

# Raises
- `DomainError`: If phi < 0, T <= 0, or g0 <= 0.
- `DimensionMismatch`: If array lengths do not conform.
"""
function solve_multispecies_escape_closure(
    phi::Real,
    X::AbstractVector{<:Real},
    m::AbstractVector{<:Real},
    T::Real,
    g0::Real,
    b::AbstractMatrix{<:Real};
    return_diag::Bool=false,
)
    phi_val = Float64(phi)
    T_val = Float64(T)
    g0_val = Float64(g0)
    N = length(X)

    if phi_val < 0.0 || !isfinite(phi_val)
        throw(DomainError(phi_val, "Mass flux phi must be non-negative and finite"))
    end
    if T_val <= 0.0 || !isfinite(T_val)
        throw(DomainError(T_val, "Exobase temperature must be positive and finite"))
    end
    if g0_val <= 0.0 || !isfinite(g0_val)
        throw(DomainError(g0_val, "Gravitational acceleration must be positive and finite"))
    end
    if length(m) != N || size(b, 1) != N || size(b, 2) != N
        throw(DimensionMismatch("Input dimensions of X, m, b do not match"))
    end

    kT = K_BOLTZMANN * T_val
    if phi_val == 0.0
        Phi = zeros(Float64, N)
        if return_diag
            return Phi, minimum(m) * g0_val / kT, Set{Int}()
        end
        return Phi
    end

    active = Set(1:N)
    wscale = phi_val / minimum(m)
    visited = Set{Set{Int}}()
    best_feasible_Phi = zeros(Float64, N)
    best_feasible_C = minimum(m) * g0_val / kT
    best_feasible_active = Set{Int}()
    min_violation = Inf

    for _ in 1:(4 * N + 8)
        if active in visited
            if !isinf(min_violation)
                if return_diag
                    return best_feasible_Phi, best_feasible_C, best_feasible_active
                end
                return best_feasible_Phi
            end
        end
        push!(visited, copy(active))

        w, C = solve_fixed_active(phi_val, X, m, T_val, g0_val, b, active)
        neg = Set([j for j in active if w[j] < -1e-12 * wscale])
        if !isempty(neg)
            setdiff!(active, neg)
            if isempty(active)
                error("Active-set iteration produced empty active set")
            end
            continue
        end

        # Retention check for inactive species
        viol = 0
        worst = 0.0
        for k in 1:N
            if k in active
                continue
            end
            Rk =
                sum(Float64(X[i]) * w[i] / Float64(b[i, k]) for i in active) -
                (Float64(m[k]) * g0_val / kT - C)
            if Rk > 1e-12 * abs(Float64(m[k]) * g0_val / kT) && Rk > worst
                viol = k
                worst = Rk
            end
        end

        if worst < min_violation
            min_violation = worst
            w_clamp = max.(w, 0.0)
            best_feasible_Phi = [Float64(X[k]) * w_clamp[k] for k in 1:N]
            for k in 1:N
                if k ∉ active
                    best_feasible_Phi[k] = 0.0
                end
            end
            best_feasible_C = C
            best_feasible_active = copy(active)
        end

        if viol == 0
            w = max.(w, 0.0)
            Phi = [Float64(X[k]) * w[k] for k in 1:N]
            for k in 1:N
                if k ∉ active
                    Phi[k] = 0.0
                end
            end
            if return_diag
                return Phi, C, active
            end
            return Phi
        end
        push!(active, viol)
    end

    if !isinf(min_violation)
        if return_diag
            return best_feasible_Phi, best_feasible_C, best_feasible_active
        end
        return best_feasible_Phi
    end
    return error("Active-set iteration did not converge")
end

"""
    compute_escape_activation_threshold(
        X::AbstractVector{<:Real}, m::AbstractVector{<:Real},
        T::Real, g0::Real, b::AbstractMatrix{<:Real},
    )::Float64

Compute the threshold mass flux phi* at which the first heavier gas entrains with the escaping wind.

# Parameters
- `X`: Species mole fractions.
- `m`: Species masses [kg].
- `T`: Exobase temperature [K].
- `g0`: Gravitational acceleration [m/s^2].
- `b`: Symmetric binary diffusion parameter matrix [m^-1 s^-1].

# Returns
- `phi_star`: Activation threshold mass flux [kg / (m^2 s)].
"""
function compute_escape_activation_threshold(
    X::AbstractVector{<:Real},
    m::AbstractVector{<:Real},
    T::Real,
    g0::Real,
    b::AbstractMatrix{<:Real},
)::Float64
    T_val = Float64(T)
    g0_val = Float64(g0)
    kT = K_BOLTZMANN * T_val
    l = argmin(m)
    best = Inf
    N = length(X)
    for k in 1:N
        if k == l
            continue
        end
        denom =
            Float64(X[l]) / Float64(b[l, k]) +
            sum(Float64(X[kk]) / Float64(b[l, kk]) for kk in 1:N if kk != l; init=0.0)
        w_star = (Float64(m[k]) - Float64(m[l])) * g0_val / kT / denom
        best = min(best, Float64(m[l]) * Float64(X[l]) * w_star)
    end
    return best
end

"""
    evolve_coupled_atmosphere_step!(
        atm_state::AtmosphereState,
        vent_rates::AbstractDict{Symbol,<:Real},
        dt_s::Real,
        M_planet::Real,
        R_planet::Real,
        T_amb::Real,
        cfg::AtmosphereConfig;
        rho_disk::Real=0.0,
        c_s::Real=300.0,
        M_star::Real=1.98847e30,
        a_orb::Real=1.495978707e11,
        T_int::Real=T_amb,
        T_exobase::Real=T_amb,
        R_exobase::Real=R_planet,
        hydrodynamic::Bool=true,
        gamma::Real=1.4,
        escape_active::Bool=true,
    )::AtmosphereState

Advance the atmospheric species inventory, gas envelope capture/boil-off, radiative equilibrium, and multispecies hydrodynamic escape over time step dt_s with mass conservation.

# Parameters
- `atm_state`: Mutable `AtmosphereState` containing atmospheric species masses and cumulative escape.
- `vent_rates`: Dictionary of species venting rates [kg/s].
- `dt_s`: Simulation time step [s].
- `M_planet`: Planetesimal mass [kg].
- `R_planet`: Planetesimal radius [m].
- `T_amb`: Ambient temperature [K].
- `cfg`: Atmosphere configuration parameters (`AtmosphereConfig`).

# Keywords
- `rho_disk`: Protoplanetary disk gas density [kg/m^3] (default: 0.0).
- `c_s`: Disk sound speed [m/s] (default: 300.0).
- `M_star`: Central star mass [kg] (default: 1.98847e30).
- `a_orb`: Planetary semi-major axis [m] (default: 1.495978707e11).
- `T_int`: Interior temperature [K] (default: T_amb).
- `T_exobase`: Exobase temperature [K] (default: T_amb).
- `R_exobase`: Exobase radius [m] (default: R_planet).
- `hydrodynamic`: Whether to include hydrodynamic blow-off (default: true).
- `gamma`: Heat capacity ratio (default: 1.4).
- `escape_active`: Whether escape is active (default: true).

# Returns
- `atm_state`: Updated `AtmosphereState`.

# Raises
- `DomainError`: If planet properties or ambient temperatures are non-positive or non-finite.
"""
function evolve_coupled_atmosphere_step!(
    atm_state::AtmosphereState,
    vent_rates::AbstractDict{Symbol,<:Real},
    dt_s::Real,
    M_planet::Real,
    R_planet::Real,
    T_amb::Real,
    cfg::AtmosphereConfig;
    rho_disk::Real=0.0,
    c_s::Real=300.0,
    M_star::Real=1.98847e30,
    a_orb::Real=1.495978707e11,
    T_int::Real=T_amb,
    T_exobase::Real=T_amb,
    R_exobase::Real=R_planet,
    hydrodynamic::Bool=true,
    gamma::Real=1.4,
    escape_active::Bool=true,
)
    dt = Float64(dt_s)
    if !isfinite(dt)
        throw(DomainError(dt, "Time step dt_s must be finite"))
    end
    if dt <= 0.0
        return atm_state
    end

    M_p = Float64(M_planet)
    R_p = Float64(R_planet)
    Tamb = Float64(T_amb)
    T_int_actual = Float64(T_int)

    if M_p <= 0.0 || !isfinite(M_p)
        throw(DomainError(M_p, "Planet mass must be > 0 and finite"))
    end
    if R_p <= 0.0 || !isfinite(R_p)
        throw(DomainError(R_p, "Planet radius must be > 0 and finite"))
    end
    if Tamb <= 0.0 || !isfinite(Tamb)
        throw(DomainError(Tamb, "Ambient temperature must be > 0 and finite"))
    end
    if T_int_actual < 0.0 || !isfinite(T_int_actual)
        throw(DomainError(T_int_actual, "Internal temperature must be >= 0 and finite"))
    end

    T_exo = max(Tamb, Float64(T_exobase))
    R_exo = max(R_p, Float64(R_exobase))
    g_surf = GRAVITATIONAL_CONSTANT * M_p / (R_p^2)
    g_exo = GRAVITATIONAL_CONSTANT * M_p / (R_exo^2)
    area = 4.0 * π * (R_p^2)
    area_exo = 4.0 * π * (R_exo^2)

    # 1. Influx from interior venting
    for (sp, rate) in vent_rates
        M_influx = Float64(rate) * dt
        atm_state.M_atm[sp] = get(atm_state.M_atm, sp, 0.0) + M_influx
    end

    # 2. Disk envelope capture and boil-off if embedded in disk or clearing
    if (rho_disk > 0.0 || atm_state.M_env_bound > 0.0) && M_p > 0.0
        M_env_target = if rho_disk > 0.0 && c_s > 0.0
            R_cap = compute_gravitational_capture_radius(M_p, M_star, a_orb, c_s)
            compute_disk_envelope_mass(M_p, R_p, R_cap, rho_disk, c_s; f_rec=cfg.f_rec)
        else
            0.0
        end
        if M_env_target > atm_state.M_env_bound
            # Envelope growth from disk gas capture (assumed H2)
            dM_cap = M_env_target - atm_state.M_env_bound
            atm_state.M_env_bound = M_env_target
            atm_state.M_atm[:H2] = get(atm_state.M_atm, :H2, 0.0) + dM_cap
        elseif M_env_target < atm_state.M_env_bound
            # Envelope boil-off during disk dispersal
            dM_boil_rate = compute_boiloff_rate(
                atm_state.M_env_bound, M_env_target, cfg.tau_boil
            )
            dM_boil = min(atm_state.M_env_bound - M_env_target, dM_boil_rate * dt)
            dM_loss_actual = min(get(atm_state.M_atm, :H2, 0.0), dM_boil)
            atm_state.M_atm[:H2] = get(atm_state.M_atm, :H2, 0.0) - dM_loss_actual
            atm_state.M_escaped[:H2] = get(atm_state.M_escaped, :H2, 0.0) + dM_loss_actual
            atm_state.M_env_bound -= dM_loss_actual
        end
    end

    # 3. Hydrodynamic escape and multispecies closure (active once disk disperses)
    if escape_active && rho_disk <= 0.0
        present_species = Symbol[sp for (sp, mass) in atm_state.M_atm if mass > 0.0]
        if !isempty(present_species)
            m_species = [get_species_molecular_mass(sp) for sp in present_species]
            min_idx = argmin(m_species)
            carrier_sp = present_species[min_idx]
            m_carrier = m_species[min_idx]
            M_carrier = atm_state.M_atm[carrier_sp]

            # Unconstrained escape of lightest species via hydrodynamic blow-off
            esc_carrier = evolve_atmospheric_species_inventory(
                M_carrier,
                0.0,
                dt,
                M_p,
                R_p,
                T_exo,
                m_carrier;
                R_exobase=R_exo,
                hydrodynamic=hydrodynamic,
                gamma=gamma,
            )
            dM_esc_base = esc_carrier.M_escaped_step

            if !cfg.crossover_active || length(present_species) == 1
                atm_state.M_atm[carrier_sp] = esc_carrier.M_atm
                atm_state.M_escaped[carrier_sp] =
                    get(atm_state.M_escaped, carrier_sp, 0.0) + dM_esc_base
                if carrier_sp === :H2
                    atm_state.M_env_bound = min(atm_state.M_env_bound, atm_state.M_atm[:H2])
                end
            elseif dM_esc_base > 0.0 && dt > 0.0 && area_exo > 0.0
                total_moles = sum(
                    atm_state.M_atm[sp] / m_species[idx] for
                    (idx, sp) in enumerate(present_species)
                )
                X_vec = [
                    (atm_state.M_atm[sp] / m_species[idx]) / total_moles for
                    (idx, sp) in enumerate(present_species)
                ]
                phi_base = (dM_esc_base / dt) / area_exo
                b_mat = assemble_binary_diffusion_matrix(present_species, T_exo)

                Phi_vec = solve_multispecies_escape_closure(
                    phi_base, X_vec, m_species, T_exo, g_exo, b_mat
                )

                for (idx, sp) in enumerate(present_species)
                    Phi_j = Phi_vec[idx]
                    dM_esc_j = min(
                        atm_state.M_atm[sp], Phi_j * m_species[idx] * area_exo * dt
                    )
                    atm_state.M_atm[sp] -= dM_esc_j
                    atm_state.M_escaped[sp] = get(atm_state.M_escaped, sp, 0.0) + dM_esc_j
                end
                if haskey(atm_state.M_atm, :H2)
                    atm_state.M_env_bound = min(atm_state.M_env_bound, atm_state.M_atm[:H2])
                end
            end
        end
    end

    # 4. Update surface diagnostics: P_surf, tau_LW, T_surf_eq, h_rad_eff, F_net_rad
    M_tot = sum(values(atm_state.M_atm))
    atm_state.P_surf = compute_surface_atmospheric_pressure(M_tot, M_p, R_p)
    atm_state.tau_LW = compute_atmospheric_optical_depth(
        atm_state.M_atm, R_p, cfg.opacities; kappa_default=cfg.kappa_ir_default
    )

    T_calc = if cfg.mode === :guillot
        compute_guillot_surface_temperature(
            atm_state.tau_LW,
            T_int_actual,
            Tamb;
            T_eqm=Tamb,
            gamma=cfg.gamma_guillot,
            albedo=cfg.albedo,
        )
    elseif cfg.mode === :isothermal
        Tamb
    else
        # Simple Eddington grey approximation
        Tamb * (1.0 + 0.75 * atm_state.tau_LW)^0.25
    end

    atm_state.T_surf_eq = max(cfg.T_skin_floor, T_calc)
    atm_state.h_rad_eff = compute_effective_radiation_htc(
        atm_state.T_surf_eq, Tamb, atm_state.tau_LW
    )
    atm_state.F_net_rad = atm_state.h_rad_eff * (T_int_actual - Tamb)

    return atm_state
end

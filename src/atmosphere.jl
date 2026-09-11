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
        (1.0 - alb) * (T_e^4)
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

"""
    evolve_coupled_atmosphere_step!(
        atm_state::AtmosphereState,
        vent_rates::AbstractDict{Symbol,<:Real},
        dt_s::Real,
        M_planet::Real,
        R_planet::Real,
        T_amb::Real,
        cfg::AtmosphereConfig;
        P_disk::Real=0.0,
        rho_disk::Real=0.0,
        c_s::Real=300.0,
        M_star::Real=1.98847e30,
        a_orb::Real=1.495978707e11,
        T_int::Real=T_amb,
        T_exobase::Real=T_amb,
        R_exobase::Real=R_planet,
    )

Advance the atmospheric species inventory, gas envelope capture/boil-off, radiative equilibrium, and hydrodynamic crossover escape over time step dt_s with machine-precision mass conservation.
"""
function evolve_coupled_atmosphere_step!(
    atm_state::AtmosphereState,
    vent_rates::AbstractDict{Symbol,<:Real},
    dt_s::Real,
    M_planet::Real,
    R_planet::Real,
    T_amb::Real,
    cfg::AtmosphereConfig;
    P_disk::Real=0.0,
    rho_disk::Real=0.0,
    c_s::Real=300.0,
    M_star::Real=1.98847e30,
    a_orb::Real=1.495978707e11,
    T_int::Real=T_amb,
    T_exobase::Real=T_amb,
    R_exobase::Real=R_planet,
)
    dt = Float64(dt_s)
    if dt <= 0.0
        return atm_state
    end

    M_p = Float64(M_planet)
    R_p = Float64(R_planet)
    Tamb = Float64(T_amb)

    if M_p <= 0.0 || !isfinite(M_p)
        throw(DomainError(M_p, "Planet mass must be > 0 and finite"))
    end
    if R_p <= 0.0 || !isfinite(R_p)
        throw(DomainError(R_p, "Planet radius must be > 0 and finite"))
    end
    if Tamb <= 0.0 || !isfinite(Tamb)
        throw(DomainError(Tamb, "Ambient temperature must be > 0 and finite"))
    end

    T_exo = max(Tamb, Float64(T_exobase))
    R_exo = max(R_p, Float64(R_exobase))
    g_surf = GRAVITATIONAL_CONSTANT * M_p / (R_p^2)
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

    # 3. Hydrodynamic escape and crossover drag
    has_h2 = haskey(atm_state.M_atm, :H2) && atm_state.M_atm[:H2] > 0.0
    if has_h2
        carrier_sp = :H2
        m_carrier = get_species_molecular_mass(:H2)
        M_carrier = atm_state.M_atm[:H2]

        # Carrier mole fraction before escape step
        total_moles_pre = sum(
            m_curr / get_species_molecular_mass(sp) for
            (sp, m_curr) in atm_state.M_atm if m_curr > 0.0
        )
        X_carrier = if total_moles_pre > 0.0
            clamp((M_carrier / m_carrier) / total_moles_pre, 0.0, 1.0)
        else
            1.0
        end

        # Carrier escape flux via hydrodynamic blow-off
        esc_carrier = evolve_atmospheric_species_inventory(
            M_carrier,
            0.0,
            dt,
            M_p,
            R_p,
            T_exo,
            m_carrier;
            R_exobase=R_exo,
            hydrodynamic=true,
        )
        dM_esc_carrier = esc_carrier.M_escaped_step
        atm_state.M_atm[:H2] = esc_carrier.M_atm
        atm_state.M_escaped[:H2] = get(atm_state.M_escaped, :H2, 0.0) + dM_esc_carrier
        atm_state.M_env_bound = min(atm_state.M_env_bound, atm_state.M_atm[:H2])

        Phi_carrier = if dt > 0.0 && area_exo > 0.0
            (dM_esc_carrier / dt) / (m_carrier * area_exo)
        else
            0.0
        end

        if cfg.crossover_active && Phi_carrier > 0.0
            m_c = compute_crossover_mass(
                m_carrier, T_exo, Phi_carrier, g_surf, X_carrier; b_diff=cfg.b_diff_ref
            )
            for (sp, m_curr) in atm_state.M_atm
                if sp != :H2 && m_curr > 0.0
                    m_sp = get_species_molecular_mass(sp)
                    x_drag = compute_crossover_drag_fraction(m_sp, m_c, m_carrier)
                    # Zahnle & Kasting (1986): momentum coupling drags species proportional to carrier loss
                    dM_drag = if (M_carrier > 0.0 && x_drag > 0.0)
                        min(m_curr, dM_esc_carrier * (m_curr / M_carrier) * x_drag)
                    else
                        0.0
                    end
                    atm_state.M_atm[sp] = m_curr - dM_drag
                    atm_state.M_escaped[sp] = get(atm_state.M_escaped, sp, 0.0) + dM_drag
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

    T_int_actual = Float64(T_int)
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

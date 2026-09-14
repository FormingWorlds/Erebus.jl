# =============================================================================
# Atmospheric Jeans Kinetic Escape & Volatile Inventory Dynamics
# =============================================================================

# Fundamental physical constants
# Note: GRAVITATIONAL_CONSTANT is defined in constants.jl and test_constants.jl (CODATA 2018)

"""
Boltzmann constant k_B in J / K (SI exact definition).
"""
const BOLTZMANN_CONSTANT = 1.380649e-23

"""
Avogadro constant N_A in 1/mol (SI exact definition).
"""
const AVOGADRO_CONSTANT = 6.02214076e23

"""
Unified atomic mass unit (amu, Dalton) in kilograms [kg] (CODATA 2018).
"""
const ATOMIC_MASS_UNIT = 1.66053906660e-27

"""
Standard atomic and molecular weights for planetary volatile and escape species [amu].
Values follow IUPAC standard atomic weights (Meija et al. 2016).
"""
const SPECIES_AMU = Dict{Symbol,Float64}(
    :H => 1.008,
    :D => 2.0141,
    :He => 4.0026,
    :C => 12.011,
    :N => 14.007,
    :O => 15.999,
    :Ne => 20.1797,
    :Na => 22.990,
    :Mg => 24.305,
    :Si => 28.085,
    :S => 32.060,
    :Ar => 39.948,
    :Fe => 55.845,
    :Kr => 83.798,
    :Xe => 131.293,
    :H2 => 2.01588,
    :H2O => 18.01528,
    :CO => 28.0101,
    :CO2 => 44.0095,
    :CH4 => 16.04246,
    :N2 => 28.0134,
    :NH3 => 17.03052,
    :O2 => 31.998,
    :H2S => 34.08088,
    :SO2 => 64.066,
    :S2 => 64.130,
)

# Volatile molecular masses [kg] (Standard atomic weights divided by Avogadro constant)
"""
Molecular mass of water (H2O) in kilograms [kg], with molar mass 18.01528 g/mol.
"""
const MASS_H2O_KG = 2.991507e-26

"""
Molecular mass of molecular hydrogen (H2) in kilograms [kg], with molar mass 2.01588 g/mol.
"""
const MASS_H2_KG = 3.347447e-27

"""
Molecular mass of molecular nitrogen (N2) in kilograms [kg], with molar mass 28.01340 g/mol.
"""
const MASS_N2_KG = 4.651735e-26

"""
Molecular mass of ammonia (NH3) in kilograms [kg], with molar mass 17.03052 g/mol.
"""
const MASS_NH3_KG = 2.827986e-26

"""
Molecular mass of carbon monoxide (CO) in kilograms [kg], with molar mass 28.01010 g/mol.
"""
const MASS_CO_KG = 4.651187e-26

"""
Molecular mass of carbon dioxide (CO2) in kilograms [kg], with molar mass 44.00950 g/mol.
"""
const MASS_CO2_KG = 7.307950e-26

"""
Molecular mass of methane (CH4) in kilograms [kg], with molar mass 16.04246 g/mol.
"""
const MASS_CH4_KG = 2.663920e-26

"""
Molecular mass of hydrogen sulfide (H2S) in kilograms [kg], with molar mass 34.08088 g/mol.
"""
const MASS_H2S_KG = 5.659267e-26

"""
Molecular mass of disulfur (S2) in kilograms [kg], with molar mass 64.130 g/mol.
"""
const MASS_S2_KG = 1.064904e-25

"""
Molecular mass of sulfur dioxide (SO2) in kilograms [kg], with molar mass 64.066 g/mol.
"""
const MASS_SO2_KG = 1.063841e-25

"""
Upper threshold on the Jeans parameter λ above which kinetic effusion is numerically negligible.
"""
const JEANS_LAMBDA_CUTOFF = 100.0

"""
Lower threshold on the Jeans parameter λ below which escape is pure hydrodynamic sound-speed blow-off.
"""
const HYDRODYNAMIC_ESCAPE_LAMBDA_LOW = 1.0

"""
Upper threshold on the Jeans parameter λ above which escape transitions fully to kinetic effusion.
"""
const HYDRODYNAMIC_ESCAPE_LAMBDA_CUTOFF = 2.0

"""
    get_species_molecular_mass(species::Symbol)::Float64

Retrieve molecular mass in kilograms for standard planetary volatile species.

# Parameters
- `species`: Volatile species identifier (`:H2O`, `:H2`, `:N2`, `:NH3`, `:CO`, `:CO2`, `:CH4`, `:H2S`, `:S2`, `:SO2`, etc.).

# Returns
- `mass`: Molecular mass [kg].
"""
function get_species_molecular_mass(species::Symbol)::Float64
    s = Symbol(uppercase(String(species)))
    if s === :H2O
        return MASS_H2O_KG
    elseif s === :H2
        return MASS_H2_KG
    elseif s === :N2
        return MASS_N2_KG
    elseif s === :NH3
        return MASS_NH3_KG
    elseif s === :CO
        return MASS_CO_KG
    elseif s === :CO2
        return MASS_CO2_KG
    elseif s === :CH4
        return MASS_CH4_KG
    elseif s === :H2S
        return MASS_H2S_KG
    elseif s === :S2
        return MASS_S2_KG
    elseif s === :SO2
        return MASS_SO2_KG
    elseif s === :H
        return SPECIES_AMU[:H] * ATOMIC_MASS_UNIT
    elseif s === :D
        return SPECIES_AMU[:D] * ATOMIC_MASS_UNIT
    elseif s === :HE
        return SPECIES_AMU[:He] * ATOMIC_MASS_UNIT
    elseif s === :C
        return SPECIES_AMU[:C] * ATOMIC_MASS_UNIT
    elseif s === :N
        return SPECIES_AMU[:N] * ATOMIC_MASS_UNIT
    elseif s === :O
        return SPECIES_AMU[:O] * ATOMIC_MASS_UNIT
    elseif s === :NE
        return SPECIES_AMU[:Ne] * ATOMIC_MASS_UNIT
    elseif s === :NA
        return SPECIES_AMU[:Na] * ATOMIC_MASS_UNIT
    elseif s === :MG
        return SPECIES_AMU[:Mg] * ATOMIC_MASS_UNIT
    elseif s === :SI
        return SPECIES_AMU[:Si] * ATOMIC_MASS_UNIT
    elseif s === :S
        return SPECIES_AMU[:S] * ATOMIC_MASS_UNIT
    elseif s === :AR
        return SPECIES_AMU[:Ar] * ATOMIC_MASS_UNIT
    elseif s === :FE
        return SPECIES_AMU[:Fe] * ATOMIC_MASS_UNIT
    elseif s === :KR
        return SPECIES_AMU[:Kr] * ATOMIC_MASS_UNIT
    elseif s === :XE
        return SPECIES_AMU[:Xe] * ATOMIC_MASS_UNIT
    elseif s === :O2
        return SPECIES_AMU[:O2] * ATOMIC_MASS_UNIT
    else
        throw(
            ArgumentError(
                "Unknown species: $species. Supported species: :H2O, :H2, :N2, :NH3, :CO, :CO2, :CH4, :H2S, :S2, :SO2, :H, :D, :He, :C, :N, :O, :Ne, :Na, :Mg, :Si, :S, :Ar, :Fe, :Kr, :Xe, :O2",
            ),
        )
    end
end

"""
Compute planetary escape velocity at a given radial distance.

$(SIGNATURES)

    v_esc = sqrt(2 * G * M / r)

# Arguments
- `M_planet::Real`: Planetary mass [kg]
- `r::Real`: Radial distance from planetary center [m]

# Returns
- `v_esc`: Escape velocity [m/s]
"""
function compute_escape_velocity(M_planet::Real, r::Real)::Float64
    M_val = require_positive_finite(Float64(M_planet), "Planetary mass")
    r_val = require_positive_finite(Float64(r), "Radius")
    return sqrt(2.0 * GRAVITATIONAL_CONSTANT * M_val / r_val)
end

"""
Compute most probable thermal speed of a gas particle.

$(SIGNATURES)

    v_th = sqrt(2 * k_B * T / m)

# Arguments
- `T_K::Real`: Gas temperature [K]
- `m_species_kg::Real`: Particle molecular mass [kg]

# Returns
- `v_th`: Most probable thermal speed [m/s]
"""
function compute_thermal_velocity(T_K::Real, m_species_kg::Real)::Float64
    T_val = require_positive_finite(Float64(T_K), "Temperature")
    m_val = require_positive_finite(Float64(m_species_kg), "Molecular mass")
    return sqrt(2.0 * BOLTZMANN_CONSTANT * T_val / m_val)
end

"""
Compute dimensionless Jeans escape parameter λ at the exobase.

$(SIGNATURES)

    λ = (G * M * m) / (k_B * T_exo * R_exo) = (v_esc / v_th)^2

# Arguments
- `M_planet::Real`: Planetary mass [kg]
- `r_exo_m::Real`: Exobase radius [m]
- `T_exo_K::Real`: Exobase temperature [K]
- `m_species_kg::Real`: Particle molecular mass [kg]

# Returns
- `lambda`: Dimensionless Jeans parameter
"""
function compute_jeans_parameter(
    M_planet::Real, r_exo_m::Real, T_exo_K::Real, m_species_kg::Real
)::Float64
    M_val = require_positive_finite(Float64(M_planet), "Planetary mass")
    r_val = require_positive_finite(Float64(r_exo_m), "Exobase radius")
    T_val = require_positive_finite(Float64(T_exo_K), "Exobase temperature")
    m_val = require_positive_finite(Float64(m_species_kg), "Molecular mass")
    return (GRAVITATIONAL_CONSTANT * M_val * m_val) / (BOLTZMANN_CONSTANT * T_val * r_val)
end

"""
Compute Jeans kinetic escape particle number flux across the exobase.

$(SIGNATURES)

    Φ = n_exo * (v_th / (2 * sqrt(π))) * (1 + λ) * exp(-λ)

# Arguments
- `n_exo::Real`: Particle number density at exobase [m^-3]
- `T_exo_K::Real`: Exobase temperature [K]
- `m_species_kg::Real`: Particle molecular mass [kg]
- `lambda::Real`: Dimensionless Jeans parameter [-]

# Keyword Arguments
- `gamma::Real`: Adiabatic index for hydrodynamic blow-off regime (default: 1.4)
- `hydrodynamic::Bool`: Enable hydrodynamic sound speed escape when λ < 2.0 (default: true)

# Returns
- `flux`: Particle number escape flux [m^-2 s^-1]
"""
function compute_jeans_escape_flux(
    n_exo::Real,
    T_exo_K::Real,
    m_species_kg::Real,
    lambda::Real;
    gamma::Real=1.4,
    hydrodynamic::Bool=true,
)::Float64
    n_val = Float64(n_exo)
    if !isfinite(n_val)
        throw(DomainError(n_val, "Number density must be finite"))
    end
    if n_val <= 0.0
        return 0.0
    end
    lam_val = require_nonneg_finite(Float64(lambda), "Jeans parameter")
    # Strong gravitational retention underflow guard
    if lam_val > JEANS_LAMBDA_CUTOFF
        return 0.0
    end
    T_val = Float64(T_exo_K)
    m_val = Float64(m_species_kg)
    if hydrodynamic && lam_val < HYDRODYNAMIC_ESCAPE_LAMBDA_CUTOFF
        c_s = sqrt(Float64(gamma) * BOLTZMANN_CONSTANT * T_val / m_val)
        flux_hydro = n_val * c_s
        v_th = compute_thermal_velocity(T_val, m_val)
        effusion_factor = (1.0 + lam_val) * exp(-lam_val)
        flux_eff = (n_val * v_th / (2.0 * sqrt(π))) * effusion_factor
        w = smoothstep(
            HYDRODYNAMIC_ESCAPE_LAMBDA_LOW, HYDRODYNAMIC_ESCAPE_LAMBDA_CUTOFF, lam_val
        )
        return lerp(flux_hydro, flux_eff, w)
    end
    v_th = compute_thermal_velocity(T_val, m_val)
    effusion_factor = (1.0 + lam_val) * exp(-lam_val)
    return (n_val * v_th / (2.0 * sqrt(π))) * effusion_factor
end

"""
Compute global planetary Jeans mass loss rate for a given volatile species.

$(SIGNATURES)

    dM/dt = 4 * π * R_exo^2 * ρ_exo * (v_th / (2 * sqrt(π))) * (1 + λ) * exp(-λ)

# Arguments
- `M_planet::Real`: Planetary mass [kg]
- `R_exo_m::Real`: Exobase radius [m]
- `T_exo_K::Real`: Exobase temperature [K]
- `m_species_kg::Real`: Particle molecular mass [kg]
- `rho_exo_kg_m3::Real`: Gas mass density at exobase [kg/m^3]

# Keyword Arguments
- `gamma::Real`: Adiabatic index for hydrodynamic blow-off regime (default: 1.4)
- `hydrodynamic::Bool`: Enable hydrodynamic sound speed escape when λ < 2.0 (default: true)

# Returns
- `loss_rate`: Mass escape rate [kg/s]
"""
function compute_jeans_mass_loss_rate(
    M_planet::Real,
    R_exo_m::Real,
    T_exo_K::Real,
    m_species_kg::Real,
    rho_exo_kg_m3::Real;
    gamma::Real=1.4,
    hydrodynamic::Bool=true,
)::Float64
    rho_val = Float64(rho_exo_kg_m3)
    if !isfinite(rho_val)
        throw(DomainError(rho_val, "Exobase density must be finite"))
    end
    if rho_val <= 0.0
        return 0.0
    end
    R_val = require_positive_finite(Float64(R_exo_m), "Exobase radius")
    lam = compute_jeans_parameter(M_planet, R_val, T_exo_K, m_species_kg)
    if lam > JEANS_LAMBDA_CUTOFF
        return 0.0
    end
    area = 4.0 * π * R_val^2
    T_val = Float64(T_exo_K)
    m_val = Float64(m_species_kg)
    if hydrodynamic && lam < HYDRODYNAMIC_ESCAPE_LAMBDA_CUTOFF
        c_s = sqrt(Float64(gamma) * BOLTZMANN_CONSTANT * T_val / m_val)
        flux_hydro = rho_val * c_s
        v_th = compute_thermal_velocity(T_val, m_val)
        effusion_factor = (1.0 + lam) * exp(-lam)
        flux_eff = rho_val * (v_th / (2.0 * sqrt(π))) * effusion_factor
        w = smoothstep(
            HYDRODYNAMIC_ESCAPE_LAMBDA_LOW, HYDRODYNAMIC_ESCAPE_LAMBDA_CUTOFF, lam
        )
        mass_flux = lerp(flux_hydro, flux_eff, w)
        return area * mass_flux
    end
    v_th = compute_thermal_velocity(T_val, m_val)
    effusion_factor = (1.0 + lam) * exp(-lam)
    mass_flux = rho_val * (v_th / (2.0 * sqrt(π))) * effusion_factor
    return area * mass_flux
end

"""
Compute atmospheric barometric scale height.

$(SIGNATURES)

    H = (k_B * T) / (m * g)

# Arguments
- `M_planet::Real`: Planetary mass [kg]
- `R_planet::Real`: Planetary radius [m]
- `T_exo_K::Real`: Exobase temperature [K]
- `m_species_kg::Real`: Particle molecular mass [kg]

# Returns
- `H`: Scale height [m]
"""
function compute_atmospheric_scale_height(
    M_planet::Real, R_planet::Real, T_exo_K::Real, m_species_kg::Real
)::Float64
    M_val = require_positive_finite(Float64(M_planet), "Planetary mass")
    R_val = require_positive_finite(Float64(R_planet), "Planetary radius")
    T_val = require_positive_finite(Float64(T_exo_K), "Temperature")
    m_val = require_positive_finite(Float64(m_species_kg), "Molecular mass")
    g = GRAVITATIONAL_CONSTANT * M_val / (R_val^2)
    return (BOLTZMANN_CONSTANT * T_val) / (m_val * g)
end

"""
Compute surface atmospheric pressure from total atmospheric mass.

$(SIGNATURES)

    P_surf = (M_atm * g) / (4 * π * R_planet^2)

# Arguments
- `M_atm_total::Real`: Total atmospheric mass [kg]
- `M_planet::Real`: Planetary mass [kg]
- `R_planet::Real`: Planetary radius [m]

# Returns
- `P_surf`: Surface pressure [Pa]
"""
function compute_surface_atmospheric_pressure(
    M_atm_total::Real, M_planet::Real, R_planet::Real
)::Float64
    M_atm_val = require_nonneg_finite(Float64(M_atm_total), "Atmospheric mass")
    M_val = require_positive_finite(Float64(M_planet), "Planetary mass")
    R_val = require_positive_finite(Float64(R_planet), "Planetary radius")
    if M_atm_val == 0.0
        return 0.0
    end
    g = GRAVITATIONAL_CONSTANT * M_val / (R_val^2)
    area = 4.0 * π * R_val^2
    return M_atm_val * g / area
end

"""
Evolve atmospheric species mass inventory over a time step dt under venting and Jeans escape.

$(SIGNATURES)

Integrates the first-order ODE:
    dM_atm / dt = M_vent_rate - k_escape * M_atm

where k_escape is determined from Jeans kinetic flux and effective atmospheric column scale.
Mass is conserved to machine precision:
    M_atm(t + dt) + M_escaped_step == M_atm(t) + M_vent_rate * dt

# Arguments
- `M_atm_prev::Real`: Initial atmospheric mass of species [kg]
- `M_vent_rate::Real`: Venting mass influx rate from interior [kg/s]
- `dt_s::Real`: Time step duration [s]
- `M_planet::Real`: Planetary mass [kg]
- `R_planet::Real`: Planetary radius [m]
- `T_exo::Real`: Exobase temperature [K]
- `m_species::Real`: Molecular mass of volatile species [kg]
- `R_exobase::Real=R_planet`: Exobase radius for escape evaluation [m]

# Keyword Arguments
- `gamma::Real`: Adiabatic index for hydrodynamic blow-off regime (default: 1.4)
- `hydrodynamic::Bool`: Enable hydrodynamic blow-off when λ < 2.0 (default: true)

# Returns
- `NamedTuple`: `(; M_atm, M_escaped_step, escape_rate)`
"""
function evolve_atmospheric_species_inventory(
    M_atm_prev::Real,
    M_vent_rate::Real,
    dt_s::Real,
    M_planet::Real,
    R_planet::Real,
    T_exo::Real,
    m_species::Real;
    R_exobase::Real=R_planet,
    gamma::Real=1.4,
    hydrodynamic::Bool=true,
)::@NamedTuple{M_atm::Float64, M_escaped_step::Float64, escape_rate::Float64}
    M_prev = Float64(M_atm_prev)
    if !isfinite(M_prev) || M_prev < 0.0
        throw(
            DomainError(M_prev, "Initial atmospheric mass must be non-negative and finite")
        )
    end
    M_dot_vent = Float64(M_vent_rate)
    if !isfinite(M_dot_vent) || M_dot_vent < 0.0
        throw(DomainError(M_dot_vent, "Venting rate must be non-negative and finite"))
    end
    dt = Float64(dt_s)
    if !isfinite(dt) || dt < 0.0
        throw(DomainError(dt, "Time step dt must be non-negative and finite"))
    end
    if dt == 0.0
        return (M_atm=M_prev, M_escaped_step=0.0, escape_rate=0.0)
    end

    R_val = Float64(R_planet)
    if R_val <= 0.0 || !isfinite(R_val)
        throw(DomainError(R_val, "Planetary radius must be > 0 and finite"))
    end
    R_exo_val = Float64(R_exobase)
    if R_exo_val < R_val || !isfinite(R_exo_val)
        throw(
            DomainError(R_exo_val, "Exobase radius must be >= planetary radius and finite")
        )
    end

    lam = compute_jeans_parameter(M_planet, R_exo_val, T_exo, m_species)
    v_th = compute_thermal_velocity(T_exo, m_species)

    # Compute loss rate coefficient k_escape [s^-1]
    k_escape = if lam > JEANS_LAMBDA_CUTOFF
        0.0
    elseif hydrodynamic && lam < HYDRODYNAMIC_ESCAPE_LAMBDA_CUTOFF
        c_s = sqrt(Float64(gamma) * BOLTZMANN_CONSTANT * Float64(T_exo) / Float64(m_species))
        k_hydro = c_s / R_exo_val
        H = compute_atmospheric_scale_height(M_planet, R_exo_val, T_exo, m_species)
        effusion_factor = (1.0 + lam) * exp(-lam)
        k_eff = (v_th / (2.0 * sqrt(π) * H)) * effusion_factor
        w = smoothstep(
            HYDRODYNAMIC_ESCAPE_LAMBDA_LOW, HYDRODYNAMIC_ESCAPE_LAMBDA_CUTOFF, lam
        )
        lerp(k_hydro, k_eff, w)
    else
        H = compute_atmospheric_scale_height(M_planet, R_exo_val, T_exo, m_species)
        effusion_factor = (1.0 + lam) * exp(-lam)
        (v_th / (2.0 * sqrt(π) * H)) * effusion_factor
    end

    x = k_escape * dt
    # Analytical solution of dM/dt = S - k*M
    M_next = if x < 1.0e-6
        # Taylor expansion to prevent cancellation or divide-by-zero
        int_factor = dt * (1.0 - 0.5 * x + (x^2) / 6.0)
        M_prev * exp(-x) + M_dot_vent * int_factor
    elseif x > 40.0
        # Fully decoupled rapid escape regime
        M_dot_vent / k_escape
    else
        M_prev * exp(-x) + (M_dot_vent / k_escape) * (1.0 - exp(-x))
    end

    # Guarantee non-negative mass
    M_next = max(0.0, M_next)

    # Mass balance: escaped mass is exactly influx plus initial minus final
    total_available = M_prev + M_dot_vent * dt
    M_escaped = max(0.0, total_available - M_next)
    current_escape_rate = k_escape * M_next

    return (M_atm=M_next, M_escaped_step=M_escaped, escape_rate=current_escape_rate)
end

"""
Compute stellar extreme ultraviolet (XUV) flux at orbital separation.

$(SIGNATURES)

Follows the empirical power-law evolution of stellar high-energy radiation (Ribas et al. 2005; Tu et al. 2015):
- Saturated regime (t <= t_sat): F_XUV = F_xuv_1au_sat / d_au^2
- Decaying regime (t > t_sat): F_XUV = (F_xuv_1au_sat / d_au^2) * (t / t_sat)^(-beta)

# Arguments
- `t_yr`: Stellar/system age in years.
- `d_au`: Orbital semi-major axis in AU.

# Keyword Arguments
- `F_xuv_1au_sat`: Saturated XUV flux at 1 AU [W/m^2] (default: 1.361 W/m^2).
- `t_sat_yr`: Duration of saturated phase in years (default: 1.0e8 yr = 100 Myr).
- `beta`: Power-law decay exponent for t > t_sat (default: 1.23, Ribas et al. 2005).

# Returns
- `F_XUV`: Incident stellar XUV flux [W/m^2].
"""
function compute_stellar_xuv_flux(
    t_yr::Real, d_au::Real; F_xuv_1au_sat::Real=1.361, t_sat_yr::Real=1.0e8, beta::Real=1.23
)::Float64
    t_val = Float64(t_yr)
    d_val = Float64(d_au)
    F_sat = Float64(F_xuv_1au_sat)
    t_sat = Float64(t_sat_yr)
    b_val = Float64(beta)

    if t_val < 0.0 || !isfinite(t_val)
        throw(DomainError(t_val, "System age t_yr must be non-negative and finite"))
    end
    if d_val <= 0.0 || !isfinite(d_val)
        throw(DomainError(d_val, "Orbital distance d_au must be > 0 and finite"))
    end
    if F_sat < 0.0 || !isfinite(F_sat)
        throw(DomainError(F_sat, "F_xuv_1au_sat must be >= 0 and finite"))
    end
    if t_sat <= 0.0 || !isfinite(t_sat)
        throw(DomainError(t_sat, "t_sat_yr must be > 0 and finite"))
    end
    if b_val < 0.0 || !isfinite(b_val)
        throw(DomainError(b_val, "Decay exponent beta must be >= 0 and finite"))
    end

    decay_factor = if t_val <= t_sat
        1.0
    else
        (t_val / t_sat)^(-b_val)
    end
    return (F_sat / (d_val^2)) * decay_factor
end

"""
Compute the Roche lobe / tidal reduction factor for hydrodynamic escape.

$(SIGNATURES)

Accounts for the reduction of the gravitational potential barrier due to stellar tidal forces
(Erkaev et al. 2007; Lammer et al. 2009):
    K_tide = 1 - 3 / (2 * xi) + 1 / (2 * xi^3)
where xi = R_Hill / R_planet and R_Hill = a * (M_planet / (3 * M_star))^(1/3).

# Arguments
- `M_planet`: Planetary mass [kg].
- `M_star`: Central stellar mass [kg].
- `a_orb`: Orbital semi-major axis [m].
- `R_planet`: Planetary surface radius [m].

# Returns
- `K_tide`: Dimensionless tidal correction factor in (0, 1].
"""
function compute_roche_lobe_correction(
    M_planet::Real, M_star::Real, a_orb::Real, R_planet::Real
)::Float64
    Mp = Float64(M_planet)
    Ms = Float64(M_star)
    a = Float64(a_orb)
    Rp = Float64(R_planet)

    if Mp <= 0.0 || !isfinite(Mp)
        throw(DomainError(Mp, "Planetary mass must be > 0 and finite"))
    end
    if Ms <= 0.0 || !isfinite(Ms)
        throw(DomainError(Ms, "Stellar mass must be > 0 and finite"))
    end
    if a <= 0.0 || !isfinite(a)
        throw(DomainError(a, "Orbital distance must be > 0 and finite"))
    end
    if Rp <= 0.0 || !isfinite(Rp)
        throw(DomainError(Rp, "Planetary radius must be > 0 and finite"))
    end

    R_Hill = a * cbrt(Mp / (3.0 * Ms))
    if Rp >= R_Hill
        throw(DomainError(Rp, "Planetary radius $Rp must be less than Hill radius $R_Hill"))
    end

    xi = R_Hill / Rp
    K_tide = 1.0 - 1.5 / xi + 0.5 / (xi^3)
    return clamp(K_tide, 0.05, 1.0)
end

"""
Compute energy-limited hydrodynamic mass loss rate and exobase base mass flux.

$(SIGNATURES)

Evaluates the photoevaporation escape rate under stellar XUV irradiation (Watson et al. 1981; Erkaev et al. 2007):
    dM_XUV / dt = (epsilon * π * R_XUV^2 * F_XUV * R_planet) / (G * M_planet * K_tide)
    phi_XUV = (dM_XUV / dt) / (4 * π * R_XUV^2) = (epsilon * F_XUV * R_planet) / (4 * G * M_planet * K_tide)

# Arguments
- `M_planet`: Planetary mass [kg].
- `R_planet`: Planetary surface radius [m].
- `F_xuv`: Incident stellar XUV flux [W/m^2].

# Keyword Arguments
- `epsilon`: Efficiency factor for converting incident XUV energy into hydrodynamic expansion (default: 0.15).
- `R_xuv`: Effective absorption radius for XUV photons [m] (default: R_planet).
- `K_tide`: Tidal potential correction factor (default: 1.0).

# Returns
- `@NamedTuple{M_dot_xuv::Float64, phi_xuv::Float64}`:
  - `M_dot_xuv`: Total energy-limited mass loss rate [kg/s].
  - `phi_xuv`: Base mass flux per unit exobase area [kg / (m^2 s)].
"""
function compute_energy_limited_escape_flux(
    M_planet::Real,
    R_planet::Real,
    F_xuv::Real;
    epsilon::Real=0.15,
    R_xuv::Real=R_planet,
    K_tide::Real=1.0,
)::@NamedTuple{M_dot_xuv::Float64, phi_xuv::Float64}
    Mp = Float64(M_planet)
    Rp = Float64(R_planet)
    Fxuv = Float64(F_xuv)
    eps = Float64(epsilon)
    Rxuv = Float64(R_xuv)
    Ktide = Float64(K_tide)

    if Mp <= 0.0 || !isfinite(Mp)
        throw(DomainError(Mp, "Planetary mass must be > 0 and finite"))
    end
    if Rp <= 0.0 || !isfinite(Rp)
        throw(DomainError(Rp, "Planetary radius must be > 0 and finite"))
    end
    if Fxuv < 0.0 || !isfinite(Fxuv)
        throw(DomainError(Fxuv, "XUV flux must be >= 0 and finite"))
    end
    if eps <= 0.0 || !isfinite(eps)
        throw(DomainError(eps, "Escape efficiency epsilon must be > 0 and finite"))
    end
    if Rxuv <= 0.0 || !isfinite(Rxuv)
        throw(DomainError(Rxuv, "R_xuv must be > 0 and finite"))
    end
    if Ktide <= 0.0 || !isfinite(Ktide)
        throw(DomainError(Ktide, "K_tide must be > 0 and finite"))
    end

    if Fxuv == 0.0
        return (M_dot_xuv=0.0, phi_xuv=0.0)
    end

    pot_well = (GRAVITATIONAL_CONSTANT * Mp * Ktide) / Rp
    absorbed_power = eps * π * (Rxuv^2) * Fxuv
    M_dot = absorbed_power / pot_well
    area_xuv = 4.0 * π * (Rxuv^2)
    phi = M_dot / area_xuv

    return (M_dot_xuv=M_dot, phi_xuv=phi)
end

"""
Compute species and isotopic fractionation factors from hydrodynamic escape fluxes.

$(SIGNATURES)

Evaluates the relative fractionation factor between species pairs:
    alpha_{ij} = (Phi_i / X_i) / (Phi_j / X_j) = w_i / w_j

# Arguments
- `Phi`: Number fluxes [molecules / (m^2 s)].
- `X`: Mole fractions.
- `species`: List of species symbols.

# Returns
- `Dict{Tuple{Symbol,Symbol},Float64}`: Fractionation factors alpha_{ij}.
"""
function compute_escape_fractionation_factors(
    Phi::AbstractVector{<:Real}, X::AbstractVector{<:Real}, species::AbstractVector{Symbol}
)::Dict{Tuple{Symbol,Symbol},Float64}
    N = length(species)
    factors = Dict{Tuple{Symbol,Symbol},Float64}()
    w = zeros(Float64, N)
    for i in 1:N
        w[i] = X[i] > 0.0 ? Float64(Phi[i]) / Float64(X[i]) : 0.0
    end
    for j in 1:N
        for i in 1:N
            if i != j
                pair = (species[i], species[j])
                alpha = if w[j] > 0.0
                    w[i] / w[j]
                elseif w[i] > 0.0
                    Inf
                else
                    0.0
                end
                factors[pair] = alpha
            end
        end
    end
    return factors
end

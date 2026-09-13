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
    M_val = Float64(M_planet)
    if M_val <= 0.0 || !isfinite(M_val)
        throw(DomainError(M_val, "Planetary mass must be > 0 and finite"))
    end
    r_val = Float64(r)
    if r_val <= 0.0 || !isfinite(r_val)
        throw(DomainError(r_val, "Radius must be > 0 and finite"))
    end
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
    T_val = Float64(T_K)
    if T_val <= 0.0 || !isfinite(T_val)
        throw(DomainError(T_val, "Temperature must be > 0 and finite"))
    end
    m_val = Float64(m_species_kg)
    if m_val <= 0.0 || !isfinite(m_val)
        throw(DomainError(m_val, "Molecular mass must be > 0 and finite"))
    end
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
    M_val = Float64(M_planet)
    if M_val <= 0.0 || !isfinite(M_val)
        throw(DomainError(M_val, "Planetary mass must be > 0 and finite"))
    end
    r_val = Float64(r_exo_m)
    if r_val <= 0.0 || !isfinite(r_val)
        throw(DomainError(r_val, "Exobase radius must be > 0 and finite"))
    end
    T_val = Float64(T_exo_K)
    if T_val <= 0.0 || !isfinite(T_val)
        throw(DomainError(T_val, "Exobase temperature must be > 0 and finite"))
    end
    m_val = Float64(m_species_kg)
    if m_val <= 0.0 || !isfinite(m_val)
        throw(DomainError(m_val, "Molecular mass must be > 0 and finite"))
    end
    return (GRAVITATIONAL_CONSTANT * M_val * m_val) / (BOLTZMANN_CONSTANT * T_val * r_val)
end

"""
Compute Jeans kinetic escape particle number flux across the exobase.

$(SIGNATURES)

    Φ_Jeans = (n_exo * v_th) / (2 * sqrt(π)) * (1 + λ) * exp(-λ)

# Arguments
- `n_exo::Real`: Particle number density at exobase [m^-3]
- `T_exo_K::Real`: Exobase temperature [K]
- `m_species_kg::Real`: Particle molecular mass [kg]
- `lambda::Real`: Dimensionless Jeans parameter

# Returns
- `Phi_Jeans`: Kinetic escape number flux [m^-2 s^-1]
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
    lam_val = Float64(lambda)
    if !isfinite(lam_val) || lam_val < 0.0
        throw(DomainError(lam_val, "Jeans parameter must be non-negative and finite"))
    end
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
        s = clamp(
            (lam_val - HYDRODYNAMIC_ESCAPE_LAMBDA_LOW) /
            (HYDRODYNAMIC_ESCAPE_LAMBDA_CUTOFF - HYDRODYNAMIC_ESCAPE_LAMBDA_LOW),
            0.0,
            1.0,
        )
        w = s * s * (3.0 - 2.0 * s)
        return (1.0 - w) * flux_hydro + w * flux_eff
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
    R_val = Float64(R_exo_m)
    if R_val <= 0.0 || !isfinite(R_val)
        throw(DomainError(R_val, "Exobase radius must be > 0 and finite"))
    end
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
        s = clamp(
            (lam - HYDRODYNAMIC_ESCAPE_LAMBDA_LOW) /
            (HYDRODYNAMIC_ESCAPE_LAMBDA_CUTOFF - HYDRODYNAMIC_ESCAPE_LAMBDA_LOW),
            0.0,
            1.0,
        )
        w = s * s * (3.0 - 2.0 * s)
        mass_flux = (1.0 - w) * flux_hydro + w * flux_eff
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
- `T_exo_K::Real`: Atmospheric temperature [K]
- `m_species_kg::Real`: Particle molecular mass [kg]

# Returns
- `H`: Scale height [m]
"""
function compute_atmospheric_scale_height(
    M_planet::Real, R_planet::Real, T_exo_K::Real, m_species_kg::Real
)::Float64
    M_val = Float64(M_planet)
    if M_val <= 0.0 || !isfinite(M_val)
        throw(DomainError(M_val, "Planetary mass must be > 0 and finite"))
    end
    R_val = Float64(R_planet)
    if R_val <= 0.0 || !isfinite(R_val)
        throw(DomainError(R_val, "Planetary radius must be > 0 and finite"))
    end
    T_val = Float64(T_exo_K)
    if T_val <= 0.0 || !isfinite(T_val)
        throw(DomainError(T_val, "Temperature must be > 0 and finite"))
    end
    m_val = Float64(m_species_kg)
    if m_val <= 0.0 || !isfinite(m_val)
        throw(DomainError(m_val, "Molecular mass must be > 0 and finite"))
    end
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
    M_atm_val = Float64(M_atm_total)
    if !isfinite(M_atm_val) || M_atm_val < 0.0
        throw(DomainError(M_atm_val, "Atmospheric mass must be non-negative and finite"))
    end
    M_val = Float64(M_planet)
    if M_val <= 0.0 || !isfinite(M_val)
        throw(DomainError(M_val, "Planetary mass must be > 0 and finite"))
    end
    R_val = Float64(R_planet)
    if R_val <= 0.0 || !isfinite(R_val)
        throw(DomainError(R_val, "Planetary radius must be > 0 and finite"))
    end
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
        s = clamp(
            (lam - HYDRODYNAMIC_ESCAPE_LAMBDA_LOW) /
            (HYDRODYNAMIC_ESCAPE_LAMBDA_CUTOFF - HYDRODYNAMIC_ESCAPE_LAMBDA_LOW),
            0.0,
            1.0,
        )
        w = s * s * (3.0 - 2.0 * s)
        (1.0 - w) * k_hydro + w * k_eff
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

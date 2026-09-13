"""
    compute_disk_temperature(t_seconds::Real, cfg::DiskConfig)::Float64

Compute ambient protoplanetary disk midplane temperature T_disk [K] at time t [s].

Supports three operating modes in `cfg.model`:
- `:fixed`: Constant ambient temperature `cfg.t_ambient`.
- `:monotonic`: Power-law viscous clearing (Lynden-Bell & Pringle 1974;
  Johansen et al. 2015).
- `:class1_to_class2`: Cold-to-hot-to-cold evolution from disk buildup
  through peak accretion to flared irradiation floor
  (Drążkowska & Dullemond 2018; Lichtenberg et al. 2021; Williams et al. 2026;
  `:class0_to_class2` supported as alias).

# Arguments
- `t_seconds`: Simulation time [s]
- `cfg`: Protoplanetary disk configuration struct (`DiskConfig`)
- `orbital_distance_au`: Optional orbital distance override [AU]
  (defaults to `cfg.orbital_distance_au`)
- `stellar_mass_msun`: Optional host star mass override [M_sun]
  (defaults to `cfg.stellar_mass_msun`)

# Returns
- `T_disk`: Ambient disk temperature [K]
"""
function compute_disk_temperature(
    t_seconds::Real,
    cfg::DiskConfig;
    orbital_distance_au::Real=cfg.orbital_distance_au,
    stellar_mass_msun::Real=cfg.stellar_mass_msun,
)::Float64
    if !cfg.enabled || cfg.model === :fixed
        return Float64(cfg.t_ambient)
    end

    t_sec_nonneg = max(0.0, Float64(t_seconds))
    t_Myr = t_sec_nonneg / (1.0e6 * (365.25 * 86400.0))
    r_au = Float64(orbital_distance_au)
    m_star = Float64(stellar_mass_msun)

    T_irr = cfg.t_irr_1au * (m_star ^ cfg.p_m_irr) * (r_au ^ (-cfg.q_irr))
    T_peak = cfg.t_peak_1au * (m_star ^ cfg.p_m_visc) * (r_au ^ (-cfg.q_visc))

    T_visc_excess4 = max(0.0, T_peak^4 - T_irr^4)

    if cfg.model === :monotonic
        t_visc = cfg.t_visc_0_myr * (m_star ^ cfg.p_m_visc_decay)
        decay = (1.0 + t_Myr / t_visc) ^ (-cfg.gamma)
        T4 = T_irr^4 + T_visc_excess4 * decay
        return max(cfg.t_cloud, T4 ^ 0.25)
    elseif cfg.model === :class1_to_class2 || cfg.model === :class0_to_class2
        t_peak = cfg.t_peak_time_1au_myr * (m_star ^ cfg.p_m_t) * (r_au ^ cfg.p_r_t)
        tau_star = 0.8 * t_peak

        if t_Myr <= 0.0
            f_acc = 0.0
            g_star = 0.0
        else
            x = t_Myr / t_peak
            f_acc =
                (1.0 + cfg.alpha / cfg.gamma) * (x ^ cfg.alpha) /
                (1.0 + (cfg.alpha / cfg.gamma) * (x ^ (cfg.alpha + cfg.gamma)))
            g_star = 1.0 - exp(-t_Myr / tau_star)
        end

        T_eff_irr4 = cfg.t_cloud^4 + (T_irr^4 - cfg.t_cloud^4) * g_star
        T4 = T_eff_irr4 + T_visc_excess4 * f_acc
        return max(cfg.t_cloud, T4 ^ 0.25)
    else
        throw(ArgumentError("Unknown disk model: $(cfg.model)"))
    end
end

"""
    compute_snowline_radius(
        t_seconds::Real, cfg::DiskConfig;
        T_sub::Real=170.0, r_min::Real=0.05, r_max::Real=100.0,
        tol::Real=1e-4, max_iter::Int=50,
        stellar_mass_msun::Real=cfg.stellar_mass_msun
    )::Float64

Compute heliocentric water snowline radius [AU] at time `t_seconds` where ambient
disk temperature equals volatile sublimation temperature `T_sub` (default: 170.0 K).
Returns `r_min` if the entire disk is below `T_sub`, or `r_max` if the disk remains
above `T_sub`. Employs a coarse-to-fine radial scan from `r_max` inward to guarantee
locating the outermost snowline under non-monotonic profiles.

# Arguments
- `t_seconds`: Simulation time [s]
- `cfg`: Protoplanetary disk configuration struct (`DiskConfig`)
- `T_sub`: Volatile sublimation threshold temperature [K] (default: 170.0 K)
- `r_min`: Minimum search radius [AU] (default: 0.05 AU)
- `r_max`: Maximum search radius [AU] (default: 100.0 AU)
- `tol`: Absolute convergence tolerance in orbital distance [AU] (default: 1e-4 AU)
- `max_iter`: Maximum bisection iterations (default: 50)
- `stellar_mass_msun`: Host star mass [M_sun] (default: `cfg.stellar_mass_msun`)

# Returns
- `r_snow`: Water snowline orbital distance [AU]
"""
function compute_snowline_radius(
    t_seconds::Real,
    cfg::DiskConfig;
    T_sub::Real=170.0,
    r_min::Real=0.05,
    r_max::Real=100.0,
    tol::Real=1e-4,
    max_iter::Int=50,
    stellar_mass_msun::Real=cfg.stellar_mass_msun,
)::Float64
    T_outer = compute_disk_temperature(
        t_seconds, cfg; orbital_distance_au=r_max, stellar_mass_msun=stellar_mass_msun
    )
    if T_outer >= T_sub
        return Float64(r_max)
    end

    # Scan radially inward from r_max to r_min across log-spaced intervals
    # to locate the outermost crossing bracket [r_lo, r_hi]
    n_coarse = 200
    log_rmax = log(Float64(r_max))
    log_rmin = log(Float64(r_min))
    dlog_r = (log_rmax - log_rmin) / n_coarse

    r_hi = Float64(r_max)
    r_lo = Float64(r_min)
    found_bracket = false

    for step in 1:n_coarse
        r_step = exp(log_rmax - step * dlog_r)
        T_step = compute_disk_temperature(
            t_seconds, cfg; orbital_distance_au=r_step, stellar_mass_msun=stellar_mass_msun
        )
        if T_step >= T_sub
            r_lo = r_step
            r_hi = exp(log_rmax - (step - 1) * dlog_r)
            found_bracket = true
            break
        end
    end

    if !found_bracket
        return Float64(r_min)
    end

    for _ in 1:max_iter
        r_mid = 0.5 * (r_lo + r_hi)
        T_mid = compute_disk_temperature(
            t_seconds, cfg; orbital_distance_au=r_mid, stellar_mass_msun=stellar_mass_msun
        )
        if abs(T_mid - T_sub) < tol || (r_hi - r_lo) < tol
            return r_mid
        end
        if T_mid >= T_sub
            r_lo = r_mid
        else
            r_hi = r_mid
        end
    end
    return 0.5 * (r_lo + r_hi)
end

"""
    compute_radiation_htc(
        T_surf::Real, T_amb::Real;
        emissivity::Real=0.9, sigma_sb::Real=5.670374419e-8
    )::Float64

Compute linearized Stefan-Boltzmann radiative heat transfer coefficient h_rad [W/(m² K)]:

    h_rad = ε * σ_SB * (T_surf² + T_amb²) * (T_surf + T_amb)

such that h_rad * (T_surf - T_amb) = ε * σ_SB * (T_surf⁴ - T_amb⁴).

# Arguments
- `T_surf`: Surface temperature [K]
- `T_amb`: Ambient disk temperature [K]
- `emissivity`: Surface thermal emissivity in [0, 1]
- `sigma_sb`: Stefan-Boltzmann constant [W/(m² K⁴)]

# Returns
- `h_rad`: Linearized radiative heat transfer coefficient [W/(m² K)]
"""
function compute_radiation_htc(
    T_surf::Real, T_amb::Real; emissivity::Real=0.9, sigma_sb::Real=5.670374419e-8
)::Float64
    if !(0.0 <= emissivity <= 1.0)
        throw(DomainError(emissivity, "Emissivity must be in [0.0, 1.0]"))
    end
    if !isfinite(T_surf) || !isfinite(T_amb) || T_surf <= 0.0 || T_amb <= 0.0
        return 0.0
    end
    T_s = Float64(T_surf)
    T_a = Float64(T_amb)
    return Float64(emissivity) * Float64(sigma_sb) * (T_s^2 + T_a^2) * (T_s + T_a)
end

"""
    compute_disk_dispersal_weight(time_seconds::Real; t_dispersal_myr::Real=3.0, dt_dispersal_myr::Real=0.1)::Float64

Compute smooth sigmoid transition weight w_disp in [0, 1] representing the fraction
of circumstellar gas disk cleared at time `time_seconds`:

    w_disp = 1 / (1 + exp(-(t_Myr - t_dispersal_myr) / dt_dispersal_myr))

# Arguments
- `time_seconds`: Simulation time [s]
- `t_dispersal_myr`: Epoch of disk gas dispersal [Myr] (default: 3.0 Myr)
- `dt_dispersal_myr`: Characteristic duration of dispersal transition [Myr] (default: 0.1 Myr)

# Returns
- `w_disp`: Sigmoid dispersal weight in [0, 1]
"""
function compute_disk_dispersal_weight(
    time_seconds::Real; t_dispersal_myr::Real=3.0, dt_dispersal_myr::Real=0.1
)::Float64
    if dt_dispersal_myr <= 0.0 || !isfinite(dt_dispersal_myr)
        throw(DomainError(dt_dispersal_myr, "dt_dispersal_myr must be > 0 and finite"))
    end
    t_sec_nonneg = max(0.0, Float64(time_seconds))
    t_Myr = t_sec_nonneg / (1.0e6 * (365.25 * 86400.0))
    t_disp = Float64(t_dispersal_myr)
    dt_disp = Float64(dt_dispersal_myr)
    arg = clamp((t_Myr - t_disp) / dt_disp, -100.0, 100.0)
    w = 1.0 / (1.0 + exp(-arg))
    return clamp(w, 0.0, 1.0)
end

"""
    compute_solar_equilibrium_temperature(orbital_distance_au::Real; albedo::Real=0.06, stellar_luminosity_lsun::Real=1.0)::Float64

Compute vacuum solar radiation equilibrium temperature T_eq [K] at heliocentric distance
`orbital_distance_au` [AU] assuming fast planetary rotation or uniform spherical emission:

    T_eq = ((1 - A) * L_star / (16 * π * σ_SB * d²))^(1/4)

# Arguments
- `orbital_distance_au`: Heliocentric orbital distance [AU]
- `albedo`: Bond albedo in [0, 1) (default: 0.06 for dark carbonaceous planetesimals)
- `stellar_luminosity_lsun`: Host star luminosity in solar units [L_sun] (default: 1.0)

# Returns
- `T_eq`: Solar radiation equilibrium temperature [K]
"""
function compute_solar_equilibrium_temperature(
    orbital_distance_au::Real; albedo::Real=0.06, stellar_luminosity_lsun::Real=1.0
)::Float64
    r_au = Float64(orbital_distance_au)
    if r_au <= 0.0 || !isfinite(r_au)
        throw(DomainError(r_au, "orbital_distance_au must be > 0 and finite"))
    end
    A = Float64(albedo)
    if !(0.0 <= A < 1.0) || !isfinite(A)
        throw(DomainError(A, "albedo must be in [0.0, 1.0)"))
    end
    L_sun = 3.828e26 * Float64(stellar_luminosity_lsun)
    sigma_sb = 5.670374419e-8
    d_m = r_au * 1.495978707e11
    F_sun = L_sun / (4.0 * π * d_m^2)
    T_eq4 = (1.0 - A) * F_sun / (4.0 * sigma_sb)
    return T_eq4^0.25
end

"""
    compute_ambient_conditions(time_seconds::Real, cfg::DiskConfig)::Tuple{Float64,Float64,Float64}

Compute evolving ambient temperature T_amb [K], ambient pressure P_amb [Pa], and disk
dispersal weight w_disp at time `time_seconds`. Transitions smoothly from nebular disk
conditions to solar radiative equilibrium and space vacuum upon disk gas clearing.

# Arguments
- `time_seconds`: Simulation time [s]
- `cfg`: Protoplanetary disk configuration struct (`DiskConfig`)

# Returns
- `(T_amb, P_amb, w_disp)`: Ambient temperature [K], ambient pressure [Pa], and dispersal weight in [0, 1]
"""
function compute_ambient_conditions(
    time_seconds::Real, cfg::DiskConfig
)::Tuple{Float64,Float64,Float64}
    w_disp = if cfg.dispersal_active
        compute_disk_dispersal_weight(
            time_seconds;
            t_dispersal_myr=cfg.t_dispersal_myr,
            dt_dispersal_myr=cfg.dt_dispersal_myr,
        )
    else
        0.0
    end
    T_disk = compute_disk_temperature(time_seconds, cfg)
    T_eq = if isfinite(cfg.t_eq_custom) && cfg.t_eq_custom > 0.0
        Float64(cfg.t_eq_custom)
    else
        compute_solar_equilibrium_temperature(cfg.orbital_distance_au; albedo=cfg.albedo)
    end
    T_amb = (1.0 - w_disp) * T_disk + w_disp * T_eq
    P_amb = (1.0 - w_disp) * cfg.p_amb_disk + w_disp * cfg.p_amb_space
    return (T_amb, P_amb, w_disp)
end

"""
    compute_ice_vapor_pressure(T::Real; P0::Real=611.66, T0::Real=273.16, L_sub::Real=2.83e6, Rv::Real=461.5)::Float64

Compute water ice sublimation equilibrium vapor pressure P_sat,ice [Pa] at temperature `T` [K]
using the integrated Clausius-Clapeyron relation:

    P_sat,ice = P0 * exp(-(L_sub / Rv) * (1/T - 1/T0))

anchored at the water triple point (T0 = 273.16 K, P0 = 611.66 Pa).

# Arguments
- `T`: Temperature [K]
- `P0`: Triple-point water vapor pressure [Pa] (default: 611.66 Pa)
- `T0`: Triple-point temperature [K] (default: 273.16 K)
- `L_sub`: Latent heat of ice sublimation [J/kg] (default: 2.83e6 J/kg)
- `Rv`: Specific gas constant for water vapor [J/(kg K)] (default: 461.5 J/(kg K))

# Returns
- `P_sat`: Equilibrium ice sublimation vapor pressure [Pa]

# Notes
- For temperatures below the triple point (`T < T0`), vapor pressure follows ice sublimation
  via the integrated Clausius-Clapeyron relation with latent heat `L_sub = 2.83e6 J/kg`.
- Above the triple point (`T >= T0`), saturation vapor pressure over liquid water follows
  the Arden Buck (1981) formulation up to the critical point (`T_crit = 647.096 K`), clamped
  at water critical pressure `P_crit = 22.064 MPa`.
"""
function compute_water_vapor_pressure(
    T::Real; P0::Real=611.66, T0::Real=273.16, L_sub::Real=2.83e6, Rv::Real=461.5
)::Float64
    T_val = Float64(T)
    if T_val <= 0.0 || !isfinite(T_val)
        throw(DomainError(T_val, "Temperature must be > 0 and finite"))
    end
    P0_val = Float64(P0)
    if P0_val <= 0.0 || !isfinite(P0_val)
        throw(DomainError(P0_val, "Triple-point pressure P0 must be > 0 and finite"))
    end
    T0_val = Float64(T0)
    if T0_val <= 0.0 || !isfinite(T0_val)
        throw(DomainError(T0_val, "Triple-point temperature T0 must be > 0 and finite"))
    end
    Rv_val = Float64(Rv)
    if Rv_val <= 0.0 || !isfinite(Rv_val)
        throw(DomainError(Rv_val, "Gas constant Rv must be > 0 and finite"))
    end
    L_sub_val = Float64(L_sub)

    if T_val <= T0_val
        return P0_val * exp(-(L_sub_val / Rv_val) * (1.0 / T_val - 1.0 / T0_val))
    end

    # Liquid water saturation vapor pressure via Arden Buck (1981)
    T_crit = 647.096
    P_crit = 22.064e6
    if T_val >= T_crit
        return P_crit
    end
    Tc = T_val - 273.15
    p_buck = 611.21 * exp((18.678 - Tc / 234.5) * (Tc / (Tc + 257.14)))
    return min(P_crit, max(P0_val, p_buck))
end

const compute_ice_vapor_pressure = compute_water_vapor_pressure

"""
Compute equilibrium saturation vapor pressure [Pa] for volatile species at temperature T_K.

Supported species:
- `:H2O`: Water ice / liquid sublimation & vapor pressure (Clausius-Clapeyron)
- `:CO2`: Carbon dioxide sublimation vapor pressure
- `:CH4`: Methane vapor pressure
- `:CO`: Carbon monoxide vapor pressure
- `:N2`: Molecular nitrogen vapor pressure
- `:H2`: Molecular hydrogen (hyper-volatile / supercritical at T >= 33 K)
- `:NH3`: Ammonia vapor pressure
- `:H2S`: Hydrogen sulfide vapor pressure
- `:SO2`: Sulfur dioxide vapor pressure
- `:S2`: Diatomic sulfur vapor pressure
"""
function compute_species_vapor_pressure(species::Symbol, T_K::Real)::Float64
    T_val = Float64(T_K)
    if T_val <= 0.0 || !isfinite(T_val)
        throw(DomainError(T_val, "Temperature must be > 0 and finite"))
    end
    if species === :H2O
        return compute_water_vapor_pressure(T_val)
    elseif species === :CO2
        # Triple point: T0 = 216.58 K, P0 = 5.18e5 Pa, L_sub = 5.71e5 J/kg, Rv = 188.92 J/(kg K)
        arg = clamp(-(5.71e5 / 188.92) * (1.0 / T_val - 1.0 / 216.58), -100.0, 100.0)
        return 5.18e5 * exp(arg)
    elseif species === :CH4
        # Non-condensible / supercritical above Tc = 190.6 K; triple point: T0 = 90.69 K, P0 = 1.17e4 Pa
        if T_val >= 190.6
            return 0.0
        end
        arg = clamp(-(5.10e5 / 518.3) * (1.0 / T_val - 1.0 / 90.69), -100.0, 100.0)
        return 1.17e4 * exp(arg)
    elseif species === :CO
        # Non-condensible / supercritical above Tc = 132.9 K; triple point: T0 = 68.15 K, P0 = 1.54e4 Pa
        if T_val >= 132.9
            return 0.0
        end
        arg = clamp(-(2.97e5 / 296.8) * (1.0 / T_val - 1.0 / 68.15), -100.0, 100.0)
        return 1.54e4 * exp(arg)
    elseif species === :N2
        # Non-condensible / supercritical above Tc = 126.2 K; triple point: T0 = 63.15 K, P0 = 1.25e4 Pa
        if T_val >= 126.2
            return 0.0
        end
        arg = clamp(-(2.54e5 / 296.8) * (1.0 / T_val - 1.0 / 63.15), -100.0, 100.0)
        return 1.25e4 * exp(arg)
    elseif species === :H2
        # Non-condensible / supercritical above Tc = 33.1 K; triple point: T0 = 13.8 K, P0 = 7.04e3 Pa
        if T_val >= 33.0
            return 0.0
        end
        arg = clamp(-(4.54e5 / 4124.0) * (1.0 / T_val - 1.0 / 13.8), -100.0, 100.0)
        return 7.04e3 * exp(arg)
    elseif species === :NH3
        # Triple point: T0 = 195.4 K, P0 = 6060.0 Pa, L_sub = 1.70e6 J/kg, Rv = 488.2 J/(kg K)
        arg = clamp(-(1.70e6 / 488.2) * (1.0 / T_val - 1.0 / 195.4), -100.0, 100.0)
        return 6060.0 * exp(arg)
    elseif species === :H2S
        # Triple point: T0 = 187.6 K, P0 = 2.32e4 Pa, L_sub = 6.98e5 J/kg, Rv = 243.9 J/(kg K)
        arg = clamp(-(6.98e5 / 243.9) * (1.0 / T_val - 1.0 / 187.6), -100.0, 100.0)
        return 2.32e4 * exp(arg)
    elseif species === :SO2
        # Triple point: T0 = 197.7 K, P0 = 1670.0 Pa, L_sub = 5.25e5 J/kg, Rv = 129.8 J/(kg K)
        arg = clamp(-(5.25e5 / 129.8) * (1.0 / T_val - 1.0 / 197.7), -100.0, 100.0)
        return 1670.0 * exp(arg)
    elseif species === :S2
        # Boiling / reference: T0 = 718.0 K, P0 = 1.0e5 Pa, L_sub = 1.45e6 J/kg, Rv = 129.6 J/(kg K)
        arg = clamp(-(1.45e6 / 129.6) * (1.0 / T_val - 1.0 / 718.0), -100.0, 100.0)
        return 1.0e5 * exp(arg)
    else
        throw(ArgumentError("Unknown volatile species: $species"))
    end
end

"""
    compute_venting_pressure(T_surf::Real, P_amb::Real; species::Symbol=:H2O, P0::Real=611.66, T0::Real=273.16, L_sub::Real=2.83e6, Rv::Real=461.5)::Float64

Compute effective boundary venting fluid pressure P_vent [Pa] at a planetesimal surface:

    P_vent = max(P_amb, P_sat(T_surf))

Enforces the physical cold-trap and boiling constraints: if ambient nebular gas pressure exceeds
saturation vapor pressure at surface temperatures, ambient gas confines pore fluid;
if ambient pressure drops below saturation vapor pressure, boiling or flash sublimation
sets the effective boundary venting pressure.

# Arguments
- `T_surf`: Planetesimal surface temperature [K]
- `P_amb`: Ambient surrounding gas pressure [Pa]

# Keyword Arguments
- `species`: Volatile species identifier (default: `:H2O`)

# Returns
- `P_vent`: Effective venting boundary pressure [Pa]
"""
function compute_venting_pressure(
    T_surf::Real,
    P_amb::Real;
    species::Symbol=:H2O,
    P0::Real=611.66,
    T0::Real=273.16,
    L_sub::Real=2.83e6,
    Rv::Real=461.5,
)::Float64
    P_sat = if species === :H2O
        compute_water_vapor_pressure(T_surf; P0=P0, T0=T0, L_sub=L_sub, Rv=Rv)
    else
        compute_species_vapor_pressure(species, T_surf)
    end
    return max(Float64(P_amb), P_sat)
end

"""
    compute_ice_sealed_permeability(
        k0::Real, T::Real;
        T_freeze::Real=273.15, delta_T_seal::Real=10.0, k_min_ratio::Real=1.0e-6
    )::Float64

Compute effective rock permeability k_sealed [m²] reduced by cryogenic pore ice freezing:

    k_sealed = k0 * [ (1 - r_min) * exp(-(T_freeze - T) / ΔT_seal) + r_min ]

For temperatures at or above freezing (`T >= T_freeze`), pore ice melts and permeability
equals `k0`. For sub-freezing temperatures (`T < T_freeze`), pore ice blocks pore throats
and reduces permeability exponentially toward the floor ratio `r_min = k_min_ratio`.

# Arguments
- `k0`: Reference unsealed permeability [m²]
- `T`: Local rock temperature [K]
- `T_freeze`: Water freezing temperature [K] (default: 273.15 K)
- `delta_T_seal`: Temperature scale for ice sealing [K] (default: 10.0 K)
- `k_min_ratio`: Minimum residual permeability floor ratio in (0, 1] (default: 1.0e-6)

# Returns
- `k_sealed`: Effective sealed permeability [m²]
"""
function compute_ice_sealed_permeability(
    k0::Real,
    T::Real;
    T_freeze::Real=273.15,
    delta_T_seal::Real=10.0,
    k_min_ratio::Real=1.0e-6,
)::Float64
    k0_val = Float64(k0)
    if k0_val <= 0.0 || !isfinite(k0_val)
        throw(DomainError(k0_val, "Reference permeability k0 must be > 0 and finite"))
    end
    T_val = Float64(T)
    if T_val <= 0.0 || !isfinite(T_val)
        throw(DomainError(T_val, "Temperature T must be > 0 and finite"))
    end
    T_frz = Float64(T_freeze)
    if T_frz <= 0.0 || !isfinite(T_frz)
        throw(DomainError(T_frz, "T_freeze must be > 0 and finite"))
    end
    dT_seal = Float64(delta_T_seal)
    if dT_seal <= 0.0 || !isfinite(dT_seal)
        throw(DomainError(dT_seal, "delta_T_seal must be > 0 and finite"))
    end
    r_min = Float64(k_min_ratio)
    if !(0.0 < r_min <= 1.0) || !isfinite(r_min)
        throw(DomainError(r_min, "k_min_ratio must be in (0, 1] and finite"))
    end

    if T_val >= T_frz
        return k0_val
    end

    arg = (T_frz - T_val) / dT_seal
    factor = (1.0 - r_min) * exp(-arg) + r_min
    return clamp(k0_val * factor, k0_val * r_min, k0_val)
end

"""
    is_hydrofracture_breached(Peff::Real, sigma_t::Real)::Bool
    is_hydrofracture_breached(Pt::Real, Pf::Real, sigma_t::Real)::Bool

Assess whether hydraulic tensile failure criterion (Peff <= -sigma_t) is satisfied,
breaching the rock matrix or cryogenic ice lid.

# Arguments
- `Peff`: Terzaghi effective pressure Pt - Pf [Pa]
- `Pt`: Total confining pressure [Pa]
- `Pf`: Pore fluid pressure [Pa]
- `sigma_t`: Rock tensile strength [Pa]

# Returns
- `breached::Bool`: `true` if hydraulic tensile fractures open, `false` otherwise.
  Non-finite inputs or non-positive `sigma_t <= 0.0` return `false` as an invalid-input safety guard.
"""
function is_hydrofracture_breached(Peff::Real, sigma_t::Real)::Bool
    if !isfinite(Peff) || !isfinite(sigma_t) || sigma_t <= 0.0
        return false
    end
    return Float64(Peff) <= -Float64(sigma_t)
end

function is_hydrofracture_breached(Pt::Real, Pf::Real, sigma_t::Real)::Bool
    return is_hydrofracture_breached(Pt - Pf, sigma_t)
end

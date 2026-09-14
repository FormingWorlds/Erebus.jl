# Erebus Synthetic Meteorite Suite
# References:
# - Van Schmus, W. R. & Wood, J. A. (1967), Geochim. Cosmochim. Acta, 31, 747-765
# - Yang, J. & Goldstein, J. I. (2006), Geochim. Cosmochim. Acta, 70, 3197-3215
# - Trieloff, M. et al. (2003), Nature, 422, 502-506
# - Huss, G. R. et al. (2006), Meteorites and the Early Solar System II, 567-586

const SEC_PER_MYR_MET = 3.15576e13
const SEC_PER_YEAR_MET = 3.15576e7

"""
Synthetic meteorite sample record for planetesimal thermal and petrologic history.

$(FIELDS)
"""
struct SyntheticMeteoriteSample
    depth_m::Float64
    radius_m::Float64
    T_peak_K::Float64
    F_melt_peak::Float64
    cooling_rate_773K_K_per_Myr::Float64
    petrologic_type::Symbol
    X_water_final::Float64
    X_refr_C_final::Float64
end

"""
Compute cooling rate between two temperature-time points.

$(SIGNATURES)

# Arguments
- `t1`: Initial time [s]
- `T1`: Initial temperature [K]
- `t2`: Final time [s]
- `T2`: Final temperature [K]

# Keyword Arguments
- `units`: Output units (`:K_per_Myr` or `:K_per_yr`, default: `:K_per_Myr`)

# Returns
- `Float64`: Cooling rate -dT/dt in requested units
"""
function compute_cooling_rate(
    t1::Real, T1::Real, t2::Real, T2::Real; units::Symbol=:K_per_Myr
)::Float64
    dt = Float64(t2 - t1)
    dt > 0.0 || throw(DomainError(dt, "t2 must be strictly greater than t1"))
    dT = Float64(T1 - T2)
    cr_per_sec = dT / dt
    if units === :K_per_Myr
        return cr_per_sec * SEC_PER_MYR_MET
    elseif units === :K_per_yr
        return cr_per_sec * SEC_PER_YEAR_MET
    else
        throw(ArgumentError("units must be :K_per_Myr or :K_per_yr, got $units"))
    end
end

"""
Compute cooling rate at a target closure temperature from temperature-time history.

$(SIGNATURES)

# Arguments
- `times_s`: Vector of simulation times [s] (monotonically increasing)
- `temps_K`: Vector of temperatures [K]
- `T_closure_K`: Target closure temperature [K] (default: 773.15 K / 500 °C)

# Keyword Arguments
- `units`: Output unit (`:K_per_Myr` or `:K_per_yr`, default: `:K_per_Myr`)

# Returns
- `Float64`: Cooling rate -dT/dt (> 0 for cooling) in requested units, or `NaN` if closure not crossed during cooling.
"""
function compute_cooling_rate(
    times_s::AbstractVector{<:Real},
    temps_K::AbstractVector{<:Real},
    T_closure_K::Real=773.15;
    units::Symbol=:K_per_Myr,
)::Float64
    n = length(times_s)
    n == length(temps_K) ||
        throw(DimensionMismatch("times_s and temps_K must have equal length"))
    n >= 2 || throw(ArgumentError("Need at least 2 points to compute cooling rate"))

    T_clos = Float64(T_closure_K)
    T_clos > 0.0 || throw(DomainError(T_clos, "T_closure_K must be positive"))

    idx_max = 1
    T_max = Float64(temps_K[1])
    for i in 2:n
        T_i = Float64(temps_K[i])
        if T_i > T_max
            T_max = T_i
            idx_max = i
        end
    end

    if T_max < T_clos
        return NaN
    end

    for i in idx_max:(n - 1)
        T_a = Float64(temps_K[i])
        T_b = Float64(temps_K[i + 1])
        t_a = Float64(times_s[i])
        t_b = Float64(times_s[i + 1])

        if T_a >= T_clos && T_b <= T_clos && t_b > t_a
            return compute_cooling_rate(t_a, T_a, t_b, T_b; units=units)
        end
    end

    return NaN
end

"""
Compute analytical conductive cooling rate proxy for a spherical body.

$(SIGNATURES)

Approximates central or depth-dependent cooling rate at temperature T_closure
for a conductive sphere of radius R with thermal diffusivity kappa.

# Arguments
- `R_m`: Body radius [m]
- `depth_m`: Depth below surface [m]

# Keyword Arguments
- `kappa`: Thermal diffusivity [m^2/s] (default: 1.0e-6)
- `delta_T`: Initial peak excess temperature [K] (default: 500.0)
- `units`: Output units (`:K_per_Myr` or `:K_per_yr`, default: `:K_per_Myr`)

# Returns
- `Float64`: Analytical cooling rate proxy [K/Myr]
"""
function compute_conductive_cooling_rate(
    R_m::Real,
    depth_m::Real;
    kappa::Real=1.0e-6,
    delta_T::Real=500.0,
    units::Symbol=:K_per_Myr,
)::Float64
    R = Float64(R_m)
    R > 0.0 || throw(DomainError(R, "R_m must be positive"))
    d = Float64(depth_m)
    (0.0 <= d <= R) || throw(DomainError(d, "depth_m must be in [0, R_m]"))
    kap = Float64(kappa)
    kap > 0.0 || throw(DomainError(kap, "kappa must be positive"))
    dT = Float64(delta_T)
    dT > 0.0 || throw(DomainError(dT, "delta_T must be positive"))

    cr_core_per_sec = (pi^2 * kap * dT) / (R^2)
    rel_d = clamp(d / R, 0.05, 1.0)
    geom_factor = sqrt(1.0 / rel_d)
    cr_per_sec = cr_core_per_sec * geom_factor

    if units === :K_per_Myr
        return cr_per_sec * SEC_PER_MYR_MET
    elseif units === :K_per_yr
        return cr_per_sec * SEC_PER_YEAR_MET
    else
        throw(ArgumentError("units must be :K_per_Myr or :K_per_yr, got $units"))
    end
end

"""
Compute metallographic cooling rate and apparent taenite Ni proxy (Yang & Goldstein 2006).

$(SIGNATURES)

Calculates the metallographic cooling rate proxy based on Ni-in-taenite diffusion
closure at ~773 K (500 °C).

# Arguments
- `cooling_rate_773K_K_per_Myr`: Thermal cooling rate at 773 K [K/Myr]

# Returns
- `NamedTuple`: `(; cooling_rate_metallographic_K_per_Myr, taenite_central_Ni_wtpct)`
"""
function compute_metallographic_cooling_rate_proxy(cooling_rate_773K_K_per_Myr::Real)
    cr = Float64(cooling_rate_773K_K_per_Myr)
    if isnan(cr) || cr <= 0.0
        return (; cooling_rate_metallographic_K_per_Myr=NaN, taenite_central_Ni_wtpct=NaN)
    end
    log_cr = log10(max(cr, 1.0e-4))
    ni_wtpct = clamp(25.0 - 5.0 * log_cr, 10.0, 45.0)
    return (; cooling_rate_metallographic_K_per_Myr=cr, taenite_central_Ni_wtpct=ni_wtpct)
end

"""
Classify synthetic meteorite sample into petrologic type (Van Schmus & Wood 1967; Huss et al. 2006).

$(SIGNATURES)

# Arguments
- `T_peak_K`: Peak metamorphic/magmatic temperature reached [K]
- `F_melt_peak`: Peak silicate melt fraction [-]

# Keyword Arguments
- `X_water`: Retained water mass fraction [-] (default: 0.0)

# Returns
- `Symbol`: Petrologic classification
  - `:achondrite`: Extensively melted/differentiated (F_melt > 0.50)
  - `:primitive_achondrite`: Partially melted (0.10 < F_melt <= 0.50)
  - `Symbol("Type 1")`: Aqueously altered, low temperature (T_peak < 420 K, X_water > 0.05)
  - `Symbol("Type 2")`: Moderately aqueously altered (T_peak < 550 K, X_water > 0.01)
  - `Symbol("Type 3")`: Pristine / unmetamorphosed (T_peak < 873.15 K)
  - `Symbol("Type 4")`: Low thermal metamorphism (873.15 K <= T_peak < 973.15 K)
  - `Symbol("Type 5")`: Intermediate thermal metamorphism (973.15 K <= T_peak < 1073.15 K)
  - `Symbol("Type 6")`: High thermal metamorphism (1073.15 K <= T_peak < 1223.15 K)
  - `Symbol("Type 7")`: Incipient melting / transitional (T_peak >= 1223.15 K, F_melt <= 0.10)
"""
function classify_petrologic_type(
    T_peak_K::Real, F_melt_peak::Real; X_water::Real=0.0
)::Symbol
    T = Float64(T_peak_K)
    T >= 0.0 || throw(DomainError(T, "T_peak_K must be non-negative"))
    F = Float64(F_melt_peak)
    (0.0 <= F <= 1.0) || throw(DomainError(F, "F_melt_peak must be in [0, 1]"))
    w = Float64(X_water)
    (0.0 <= w <= 1.0) || throw(DomainError(w, "X_water must be in [0, 1]"))

    if F > 0.50
        return :achondrite
    elseif F > 0.10
        return :primitive_achondrite
    end

    if T < 420.0 && w > 0.05
        return Symbol("Type 1")
    elseif T < 550.0 && w > 0.01
        return Symbol("Type 2")
    elseif T < 873.15
        return Symbol("Type 3")
    elseif T < 973.15
        return Symbol("Type 4")
    elseif T < 1073.15
        return Symbol("Type 5")
    elseif T < 1223.15
        return Symbol("Type 6")
    else
        return Symbol("Type 7")
    end
end

"""
Generate synthetic meteorite sample suite from planetesimal radial/marker data.

$(SIGNATURES)

# Arguments
- `r_m`: Marker or shell radius vector [m]
- `T_peak_m`: Peak temperature vector [K]
- `F_melt_peak_m`: Peak melt fraction vector
- `cooling_rates_K_per_Myr`: Cooling rate at 773 K vector [K/Myr]
- `R_planet_m`: Planetesimal outer radius [m]

# Keyword Arguments
- `X_water_m`: Retained water fraction vector (optional)
- `X_refr_C_m`: Retained refractory carbon fraction vector (optional)

# Returns
- `Vector{SyntheticMeteoriteSample}`
"""
function generate_synthetic_meteorite_suite(
    r_m::AbstractVector{<:Real},
    T_peak_m::AbstractVector{<:Real},
    F_melt_peak_m::AbstractVector{<:Real},
    cooling_rates_K_per_Myr::AbstractVector{<:Real},
    R_planet_m::Real;
    X_water_m::Union{Nothing,AbstractVector{<:Real}}=nothing,
    X_refr_C_m::Union{Nothing,AbstractVector{<:Real}}=nothing,
)::Vector{SyntheticMeteoriteSample}
    n = length(r_m)
    (
        length(T_peak_m) == n &&
        length(F_melt_peak_m) == n &&
        length(cooling_rates_K_per_Myr) == n
    ) || throw(DimensionMismatch("Input vectors must have identical length"))
    R_p = Float64(R_planet_m)
    R_p > 0.0 || throw(DomainError(R_p, "R_planet_m must be positive"))

    samples = Vector{SyntheticMeteoriteSample}(undef, n)

    for i in 1:n
        r_val = Float64(r_m[i])
        depth_val = max(0.0, R_p - r_val)
        T_pk = Float64(T_peak_m[i])
        F_pk = Float64(F_melt_peak_m[i])
        cr_val = Float64(cooling_rates_K_per_Myr[i])
        w_h2o = X_water_m !== nothing ? Float64(X_water_m[i]) : 0.0
        w_c = X_refr_C_m !== nothing ? Float64(X_refr_C_m[i]) : 0.0

        p_type = classify_petrologic_type(T_pk, F_pk; X_water=w_h2o)

        samples[i] = SyntheticMeteoriteSample(
            depth_val, r_val, T_pk, F_pk, cr_val, p_type, w_h2o, w_c
        )
    end

    return samples
end


"""
Calculate Euclidean distance between two point coordinates.

$(SIGNATURES)

# Details

    - x1: x-coordinate of point 1 [m]
    - y1: y-coordinate of point 1 [m]
    - x2: x-coordinate of point 2 [m]
    - y2: y-coordinate of point 2 [m]

# Returns

    - Euclidean distance between point 1 and point 2 [m]
"""
function distance(x1, y1, x2, y2)
    return sqrt(abs2(x1-x2) + abs2(y1-y2))
end

"""
Compute convex combination of fluid and solid properties to get total property.

$(SIGNATURES)

# Details

    - fluid: fluid properties
    - solid: solid properties
    - ϕ: porosity (fraction of fluid)

# Returns

    - total: computed total property
"""
function total(solid, fluid, ϕ)
    if !(0.0 <= ϕ <= 1.0)
        throw(DomainError(ϕ, "Porosity must be in [0, 1]"))
    end
    return solid * (1.0 - ϕ) + fluid * ϕ
end

"""
Compute total thermal conductivity of two-phase material.

$(SIGNATURES)

# Details

    - ksolid: solid thermal conductivity [W/m/K]
    - kfluid: fluid thermal conductivity [W/m/K]
    - phi: porosity (fraction of fluid)

# Returns

    - ktotal: total thermal conductivity of mixed phase [W/m/K]
"""
function ktotal(ksolid, kfluid, phi)
    if !(0.0 <= phi <= 1.0)
        throw(DomainError(phi, "Porosity must be in [0, 1]"))
    end
    if ksolid < 0.0 || kfluid < 0.0
        throw(DomainError((ksolid, kfluid), "Thermal conductivities must be non-negative"))
    end
    return (
        sqrt(
            ksolid * kfluid / 2 +
            ((ksolid * (3.0 * phi - 2.0) + kfluid * (1.0 - 3.0 * phi))^2) * inv(16.0),
        ) - 0.25 * (ksolid * (3.0 * phi - 2.0) + kfluid * (1.0 - 3.0 * phi))
    )
end

"""
Compute porosity-dependent permeability (eqn 16.64 in Gerya (2019)).

$(SIGNATURES)

# Details

    - kphim0m: standard (reference) permeability (of marker type) [m^2]
    - phimm: actual (marker) porosity

# Returns

    - kphim: empirical porosity-dependent permeability [m^2]
"""
function kphi(kphim0m, phimm)
    if !(0.0 <= phimm < 1.0)
        throw(DomainError(phimm, "Porosity must be in [0, 1)"))
    end
    if kphim0m < 0.0
        throw(DomainError(kphim0m, "Reference permeability must be non-negative"))
    end
    # phim0 is a global constant defined independent of material type
    return kphim0m * (phimm * inv(phim0))^3.0 * ((1.0 - phimm) * inv(1.0 - phim0))^-2.0
end

"""
Compute inverse of porosity-dependent permeability (eqn 16.64 in Gerya (2019))
times current fluid viscosity.

$(SIGNATURES)

# Details

    - kϕᵣ: reference permeability [m^2]
    - ϕ: current porosity
    - ηᶠcur: current fluid viscosity [Pa s]

# Returns

    - etafluidcur_inv_kphi: inverse empirical porosity-dependent permeability
                            times current fluid viscosity
"""
function ηᶠcur_inv_kᵠ(kϕᵣ, ϕ, ηᶠcur)
    if !(0.0 < ϕ < 1.0)
        throw(DomainError(ϕ, "Porosity must be strictly between 0 and 1"))
    end
    if kϕᵣ <= 0.0 || ηᶠcur < 0.0
        throw(
            DomainError(
                (kϕᵣ, ηᶠcur), "Permeability must be positive and viscosity non-negative"
            ),
        )
    end
    return ηᶠcur * inv(kϕᵣ) * (phim0 * inv(ϕ))^3.0 * ((1.0 - ϕ) * inv(1.0 - phim0))^2.0
end

"""
Compute total rocky marker viscosity based on temperature and material type.

$(SIGNATURES)

# Details

    - tkmm: marker temperature [K]
    - tmm: marker type [1, 2]

# Returns

    - etatotal: rocky marker temperature-dependent total viscosity
"""
function etatotal_rocks(tkmm, tmm; etamin::Real=1.0e12)
    if tkmm <= 0.0
        throw(DomainError(tkmm, "Absolute temperature must be positive"))
    end
    @inbounds etasolidcur = ifelse(tkmm > tmsolidphase, etasolidmm[tmm], etasolidm[tmm])
    @inbounds etafluidcur = ifelse(tkmm > tmfluidphase, etafluidmm[tmm], etafluidm[tmm])
    return max(etamin, etasolidcur, etafluidcur)
end

"""
Compute volumetric isobaric heat capacity of H₂O (fluid phase)
based on temperature.

$(SIGNATURES)

# Details

        - T: temperature [K]
        - mode: marker property computation mode
            - 1: dynamic, based on (Touloukian, 1970; Hobbs, 1974;
                 Travis and Schubert, 2005)
            - 9: constant parameter rhocpfluidm

# Returns

        - ρᶠCₚᶠ: volumetric isobaric heat capacity of fluid
"""
function compute_rhocpfluidm(T, mode)
    if mode == 1
        if T <= 0.0
            throw(DomainError(T, "Absolute temperature must be positive"))
        end
        if T < tmfluidphase - 5.0
            ρᶠCₚᶠ = ρH₂Oᶠⁱ * 7.67T
        elseif T < tmfluidphase
            ρᶠCₚᶠ = ρH₂Oᶠⁱ * (7.67T + 0.1Lᶠ)
        elseif T < tmfluidphase + 5.0
            ρᶠCₚᶠ = ρH₂Oᶠ * (4200.0 + 0.1Lᶠ)
        elseif T < 410.0
            ρᶠCₚᶠ = ρH₂Oᶠ * 4200.0
        else
            ρᶠCₚᶠ = ρH₂Oᶠ * (-4.67e4 + 333T - 0.731T^2 + 5.4e-4T^3)
        end
    elseif mode == 9
        @inbounds ρᶠCₚᶠ = rhocpfluidm[1]
    else
        throw(ArgumentError("unknown mode $mode"))
    end
    return ρᶠCₚᶠ
end # function compute_rhocpfluidm

"""
Compute thermal conductivity of silicate (solid phase) based on temperature.

$(SIGNATURES)

# Details

        - T: temperature [K]
        - mode:
            - 1: dynamic, based on (Gerya, 2019)
            - 9: constant parameter ksolidm

# Returns

        - kᶠ: thermal conductivity of solid
"""
function compute_ksolidm(T, mode)
    if mode == 1
        if T <= 0.0
            throw(DomainError(T, "Absolute temperature must be positive"))
        end
        kˢ = 0.73 + 1293.0 / (T + 77.0)
    elseif mode == 9
        @inbounds kˢ = ksolidm[1]
    else
        throw(ArgumentError("unknown mode $mode"))
    end
    return kˢ
end # function compute_ksolidm

"""
Compute thermal conductivity of H₂O (fluid phase) based on temperature.

$(SIGNATURES)

# Details

        - T: temperature [K]
        - mode:
            - 1: dynamic, based on (Touloukian, 1970; Hobbs, 1974;
                 Grimm & Mcsween, 1989; Bland & Travis, 2017)
            - 9: constant parameter kfluidm

# Returns

        - kᶠ: thermal conductivity of fluid
"""
function compute_kfluidm(T, mode)
    if mode == 1
        if T <= 0.0
            throw(DomainError(T, "Absolute temperature must be positive"))
        end
        if T < tmfluidphase
            kᶠ = 0.465 + 488.0 / T
        elseif T < 410.0
            kᶠ = -0.581 + 6.34e-3T - 7.93e-6T^2
        else
            kᶠ = -0.142 + 4.12e-3T - 5.01e-6T^2
        end
    elseif mode == 9
        @inbounds kᶠ = kfluidm[1]
    else
        throw(ArgumentError("unknown mode $mode"))
    end
    return kᶠ
end # function compute_kfluidm

"""
    compute_spherical_metric_heat_source!(Q_metric, tk, KX, KY, coords;
                                          xcenter, ycenter, rplanet, reg_cells=0.5)

Compute geometric metric volumetric heat source Q_metric [W/m³] on P-nodes
to account for 3D spherical divergence on a 2D Cartesian grid:

    Q_metric = (k / r_eff²) * [ (x - xc) * (∂T/∂x) + (y - yc) * (∂T/∂y) ]

inside r <= rplanet, and 0 in sticky air.

# Arguments
- `Q_metric`: Output matrix [W/m³] of size (Ny1, Nx1)
- `tk`: Temperature field [K] of size (Ny1, Nx1)
- `KX`: Thermal conductivity at Vx nodes [W/(m K)]
- `KY`: Thermal conductivity at Vy nodes [W/(m K)]
- `coords`: Grid coordinate descriptors
- `xcenter`: Horizontal center coordinate [m]
- `ycenter`: Vertical center coordinate [m]
- `rplanet`: Planetesimal radius [m]
- `reg_cells`: Radial regularization parameter in grid cell units
"""
function compute_spherical_metric_heat_source!(
    Q_metric::AbstractMatrix{Float64},
    tk::AbstractMatrix{Float64},
    KX::AbstractMatrix{Float64},
    KY::AbstractMatrix{Float64},
    coords::GridCoordinates;
    xcenter::Real,
    ycenter::Real,
    rplanet::Real,
    reg_cells::Real=0.5,
)
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    dx = coords.dx
    dy = coords.dy
    inv_2dx = inv(2.0 * dx)
    inv_2dy = inv(2.0 * dy)
    eps_r = reg_cells * min(dx, dy)
    eps_r2 = eps_r^2
    rplanet2 = rplanet^2

    fill!(Q_metric, 0.0)

    @inbounds for j in 2:(Nx1 - 1)
        xj = coords.xp[j]
        dx_c = xj - xcenter
        for i in 2:(Ny1 - 1)
            yi = coords.yp[i]
            dy_c = yi - ycenter
            r2 = dx_c^2 + dy_c^2
            if r2 <= rplanet2
                k_node = 0.25 * (KX[i, j - 1] + KX[i, j] + KY[i - 1, j] + KY[i, j])
                dT_dx = (tk[i, j + 1] - tk[i, j - 1]) * inv_2dx
                dT_dy = (tk[i + 1, j] - tk[i - 1, j]) * inv_2dy
                r_eff2 = r2 + eps_r2
                Q_metric[i, j] = (k_node / r_eff2) * (dx_c * dT_dx + dy_c * dT_dy)
            end
        end
    end
    return Q_metric
end

"""
    lerp(a, b, t)

Linear interpolation between `a` and `b` with blend factor `t`.
"""
@inline lerp(a, b, t) = a * (1.0 - t) + b * t

"""
    smoothstep(x0, x1, x)

Cubic smoothstep Hermite interpolation between 0 and 1 for `x` clamped to `[x0, x1]`.
"""
@inline function smoothstep(x0::Real, x1::Real, x::Real)
    x0 == x1 && return x >= x1 ? 1.0 : 0.0
    xi = clamp((x - x0) / (x1 - x0), 0.0, 1.0)
    return xi * xi * (3.0 - 2.0 * xi)
end

"""
    require_positive_finite(v, name)

Validate that `v` is positive and finite, otherwise throw `DomainError`.
"""
@inline function require_positive_finite(v::Real, name::AbstractString)
    (isfinite(v) && v > 0.0) ||
        throw(DomainError(v, string(name, " must be > 0 and finite")))
    return v
end

"""
    require_nonneg_finite(v, name)

Validate that `v` is non-negative and finite, otherwise throw `DomainError`.
"""
@inline function require_nonneg_finite(v::Real, name::AbstractString)
    (isfinite(v) && v >= 0.0) ||
        throw(DomainError(v, string(name, " must be non-negative and finite")))
    return v
end

"""
    require_unit_interval(v, name)

Validate that `v` is in unit interval [0, 1] and finite, otherwise throw `DomainError`.
"""
@inline function require_unit_interval(v::Real, name::AbstractString)
    (isfinite(v) && 0.0 <= v <= 1.0) ||
        throw(DomainError(v, string(name, " must be in [0, 1] and finite")))
    return v
end

"""
    compute_l3d_metric(rplanet::Real)

Compute the scaling factor `L_3D` used to convert a 2D out-of-plane extruded mass inventory [kg/m]
into a 3D spherical inventory [kg].

Derived from the ratio of spherical volume to 2D disk cross-sectional area:
`L_3D = (4/3 π R^3) / (π R^2) = (4/3) R`
"""
function compute_l3d_metric(rplanet::Real)
    rplanet <= 0.0 && throw(DomainError(rplanet, "Planet radius must be positive"))
    return (4.0 / 3.0) * rplanet
end

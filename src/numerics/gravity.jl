"""
Assemble the LHS sparse coefficient matrix and fill RHS coefficient vector
of the Poisson equation to be solved for the gravitational potential Φ.

$(SIGNATURES)

# Details

    - RHO: density at P nodes
    - RP: right hand side coefficient vector
    - coords: grid coordinates
    - LP: optional ExtendableSparseMatrix buffer to reuse

# Returns

    - LP: LHS sparse coefficient matrix

"""
function assemble_gravitational_lse!(
    RHO, RP; coords::GridCoordinates=default_grid_coordinates(), LP=nothing
)
    Ny1, Nx1 = size(RHO)
    @unpack_coords coords dx dy xp yp
    xc_val = coords.xcenter
    yc_val = coords.ycenter
    r_limit = min(xc_val, yc_val)

    L = if LP === nothing
        ExtendableSparseMatrix(Nx1 * Ny1, Nx1 * Ny1)
    else
        if !isempty(LP.cscmatrix.nzval)
            nonzeros(LP.cscmatrix) .= zero(0.0)
        end
        LP
    end
    # reset RHS coefficient vector
    RP .= zero(0.0)
    G_coeff = 4.0 * 2.0 * inv(3.0) * π * G
    # iterate over P nodes
    for j in 1:1:Nx1, i in 1:1:Ny1
        # define global index in algebraic space
        gk = (j - 1) * Ny1 + i
        # decide if external / boundary points
        @inbounds if is_gravitational_boundary(
            i, j, Ny1, Nx1, xp_val, yp_val, xc_val, yc_val, r_limit
        )
            # boundary condition: ϕ = 0
            updateindex!(L, +, 1.0, gk, gk)
        else
            # internal points: 2D Poisson equation: gravitational potential Φ
            updateindex!(L, +, inv(dx_val^2), gk, gk - Ny1) # Φ₁
            updateindex!(L, +, inv(dy_val^2), gk, gk - 1) # Φ₂
            updateindex!(L, +, -2.0 * (inv(dx_val^2) + inv(dy_val^2)), gk, gk) # Φ₃
            updateindex!(L, +, inv(dy_val^2), gk, gk + 1) # Φ₄
            updateindex!(L, +, inv(dx_val^2), gk, gk + Ny1) # Φ₅
            @inbounds RP[gk] = G_coeff * RHO[i, j]
        end
    end
    flush!(L)
    return L
end

"""
Assemble right-hand side coefficient vector for gravitational Poisson equation.

$(SIGNATURES)

# Details

    - RHO: density at P nodes
    - RP: right hand side coefficient vector to populate
    - coords: grid coordinates

# Returns

    - RP
"""
function assemble_gravitational_rhs!(
    RHO, RP; coords::GridCoordinates=default_grid_coordinates()
)
    Ny1, Nx1 = size(RHO)
    @unpack_coords coords xp yp
    xc_val = coords.xcenter
    yc_val = coords.ycenter
    r_limit = min(xc_val, yc_val)

    RP .= zero(0.0)
    G_coeff = 4.0 * 2.0 * inv(3.0) * π * G
    for j in 1:1:Nx1, i in 1:1:Ny1
        @inbounds if !is_gravitational_boundary(
            i, j, Ny1, Nx1, xp_val, yp_val, xc_val, yc_val, r_limit
        )
            gk = (j - 1) * Ny1 + i
            RP[gk] = G_coeff * RHO[i, j]
        end
    end
    return RP
end

"""
Process gravitational potential solution vector to output physical observables.

$(SIGNATURES)

# Details

    - FI: gravitational potential
    - gx: x-component of gravitational acceleration
    - gy: y-component of gravitational acceleration

# Returns

    - nothing
"""
function process_gravitational_solution!(
    SP, FI, gx, gy; coords::GridCoordinates=default_grid_coordinates()
)
    Ny1, Nx1 = size(FI)
    Nx_val = Nx1 - 1
    Ny_val = Ny1 - 1
    inv_dx = inv(coords.dx)
    inv_dy = inv(coords.dy)
    @inbounds for j in 1:Nx1, i in 1:Ny1
        FI[i, j] = SP[(j - 1) * Ny1 + i]
    end
    @inbounds for j in 1:Nx_val, i in 1:Ny1
        gx[i, j] = -(FI[i, j + 1] - FI[i, j]) * inv_dx
    end
    @inbounds for j in 1:Nx1, i in 1:Ny_val
        gy[i, j] = -(FI[i + 1, j] - FI[i, j]) * inv_dy
    end
    @inbounds for i in 1:Ny1
        gx[i, Nx1] = gx[i, Nx_val]
    end
    @inbounds for j in 1:Nx1
        gy[Ny1, j] = gy[Ny_val, j]
    end
    return nothing
end

"""
Compute gravity solution in P nodes to obtain
gravitational accelerations gx for Vx nodes, gy for Vy nodes.

$(SIGNATURES)

# Details

    - SP: solution vector
    - RP: right hand side vector
    - RHO: density at P nodes
    - FI: gravity potential at P nodes
    - gx: x gravitational acceleration at Vx nodes
    - gy: y gravitational acceleration at Vy nodes

# Returns

- nothing
"""
function compute_gravity_solution!(
    SP, RP, RHO, FI, gx, gy; coords::GridCoordinates=default_grid_coordinates()
)
    LP = assemble_gravitational_lse!(RHO, RP; coords=coords)
    SP .= LP \ RP
    process_gravitational_solution!(SP, FI, gx, gy; coords=coords)
    return nothing
end # function compute_gravity_solution!

"""
Evaluate piecewise-linear radial gravitational acceleration magnitude `g(r)` with linear core regularization below `r1`.

$(SIGNATURES)

# Parameters
- `r`: radial distance from body center [m]
- `r1`: outer radius of innermost radial shell [m]
- `r_max`: maximum radial shell boundary [m]
- `dr`: radial shell thickness [m]
- `Nr`: number of radial shells
- `g1`: gravitational acceleration at `r1` [m/s^2]
- `M_tot`: total enclosed planetary mass [kg]
- `g_bins`: radial shell gravitational acceleration vector [m/s^2]

# Returns
- `Float64`: radial gravitational acceleration magnitude [m/s^2]
"""
@inline function _eval_radial_gravity(
    r::Float64,
    r1::Float64,
    r_max::Float64,
    dr::Float64,
    Nr::Int,
    g1::Float64,
    M_tot::Float64,
    g_bins::Vector{Float64},
)::Float64
    if r <= 0.0
        return 0.0
    elseif r < r1
        return g1 * (r / r1)
    elseif r >= r_max
        return G * M_tot / (r^2)
    else
        k = clamp(floor(Int, r / dr), 1, Nr - 1)
        rk_prev = k * dr
        frac = (r - rk_prev) / dr
        return @inbounds g_bins[k] + frac * (g_bins[k + 1] - g_bins[k])
    end
end

"""
Evaluate piecewise-linear radial gravitational potential `phi(r)` with parabolic core regularization below `r1`.

$(SIGNATURES)

# Parameters
- `r`: radial distance from body center [m]
- `r1`: outer radius of innermost radial shell [m]
- `r_max`: maximum radial shell boundary [m]
- `dr`: radial shell thickness [m]
- `Nr`: number of radial shells
- `g1`: gravitational acceleration at `r1` [m/s^2]
- `phi0`: gravitational potential at body center `r = 0` [J/kg]
- `M_tot`: total enclosed planetary mass [kg]
- `phi_bins`: radial shell gravitational potential vector [J/kg]

# Returns
- `Float64`: gravitational potential [J/kg]
"""
@inline function _eval_radial_potential(
    r::Float64,
    r1::Float64,
    r_max::Float64,
    dr::Float64,
    Nr::Int,
    g1::Float64,
    phi0::Float64,
    M_tot::Float64,
    phi_bins::Vector{Float64},
)::Float64
    if r >= r_max
        return -G * M_tot / r
    elseif r <= 0.0
        return phi0
    elseif r < r1
        return @inbounds phi_bins[1] - 0.5 * g1 * (r1 - r^2 / r1)
    else
        k = clamp(floor(Int, r / dr), 1, Nr - 1)
        rk_prev = k * dr
        frac = (r - rk_prev) / dr
        return @inbounds phi_bins[k] + frac * (phi_bins[k + 1] - phi_bins[k])
    end
end

"""
Compute gravitational acceleration using the 3D spherical enclosed-mass formulation.

$(SIGNATURES)

# Arguments
- `gx::AbstractMatrix{Float64}`: Horizontal gravitational acceleration at Vx nodes [m/s²].
- `gy::AbstractMatrix{Float64}`: Vertical gravitational acceleration at Vy nodes [m/s²].

# Keywords
- `xm::AbstractVector{Float64}`: Horizontal marker coordinates [m].
- `ym::AbstractVector{Float64}`: Vertical marker coordinates [m].
- `rhototalm::AbstractVector{Float64}`: Marker total density [kg/m³].
- `tm::AbstractVector{<:Integer}`: Marker material phase tag (`tm < 3` for rock/metal, `tm >= 3` for sticky air).
- `coords::GridCoordinates`: Grid coordinate arrays and geometry.
- `gravity_nr_factor::Int = 4`: Radial bin resolution multiplier (`Nr = gravity_nr_factor * Nx`).
- `rplanet::Union{Nothing,Real} = nothing`: Optional planet radius [m].
- `FI::Union{Nothing,AbstractMatrix{Float64}} = nothing`: Optional gravitational potential at P nodes [J/kg].

# Returns
- `Tuple{Vector{Float64}, Vector{Float64}}`: `(r_bins, g_bins)` radial coordinates and gravity profile.
"""
function compute_gravity_enclosed_mass!(
    gx::AbstractMatrix{Float64},
    gy::AbstractMatrix{Float64};
    xm::AbstractVector{Float64},
    ym::AbstractVector{Float64},
    rhototalm::AbstractVector{Float64},
    tm::AbstractVector{<:Integer},
    coords::GridCoordinates=default_grid_coordinates(),
    gravity_nr_factor::Int=4,
    rplanet::Union{Nothing,Real}=nothing,
    FI::Union{Nothing,AbstractMatrix{Float64}}=nothing,
)
    xc = coords.xcenter
    yc = coords.ycenter
    Nx_val = coords.Nx
    dxm_val = coords.dxm
    dym_val = coords.dym
    Am = dxm_val * dym_val

    # Determine maximum radius for enclosed radial bins
    r_max = if rplanet !== nothing
        Float64(rplanet)
    else
        coords.xsize / 2.0
    end
    for m in eachindex(xm, ym, tm)
        if @inbounds tm[m] < 3
            r_m = hypot(xm[m] - xc, ym[m] - yc)
            if r_m > r_max
                r_max = r_m
            end
        end
    end

    Nr = max(1, gravity_nr_factor * Nx_val)
    dr = r_max / Nr
    r1 = dr

    # Accumulate 3D shell mass increments: M_3D(<r) = sum rho_m * Am * 2 * r_m
    bin_mass = zeros(Float64, Nr)
    for m in eachindex(xm, ym, rhototalm, tm)
        if @inbounds tm[m] < 3
            r_m = hypot(xm[m] - xc, ym[m] - yc)
            k = clamp(ceil(Int, r_m / dr), 1, Nr)
            @inbounds bin_mass[k] += rhototalm[m] * Am * 2.0 * r_m
        end
    end

    M_enc = cumsum(bin_mass)
    M_tot = M_enc[end]

    # Fill empty center bins from discrete marker sampling
    k0 = findfirst(>(0.0), M_enc)
    if k0 !== nothing && k0 > 1
        r_k0 = k0 * dr
        rho_bar = M_enc[k0] / ((4.0 / 3.0) * π * (r_k0^3))
        for k in 1:(k0 - 1)
            rk = k * dr
            M_enc[k] = (4.0 / 3.0) * π * rho_bar * (rk^3)
        end
    end

    r_bins = [k * dr for k in 1:Nr]
    g_bins = zeros(Float64, Nr)
    for k in 1:Nr
        g_bins[k] = G * M_enc[k] / (r_bins[k]^2)
    end
    g1 = g_bins[1]

    # Apply radial gravity vector to Vx nodes (gx) and Vy nodes (gy)
    Ny1, Nx1 = size(gx)
    @unpack_coords coords xvx yvx xvy yvy xp yp

    @inbounds for j in 1:Nx1, i in 1:Ny1
        dx_val = xvx_val[j] - xc
        dy_val = yvx_val[i] - yc
        r_node = hypot(dx_val, dy_val)
        if r_node == 0.0
            gx[i, j] = 0.0
        else
            g_mag = _eval_radial_gravity(r_node, r1, r_max, dr, Nr, g1, M_tot, g_bins)
            gx[i, j] = -g_mag * (dx_val / r_node)
        end
    end

    @inbounds for j in 1:Nx1, i in 1:Ny1
        dx_val = xvy_val[j] - xc
        dy_val = yvy_val[i] - yc
        r_node = hypot(dx_val, dy_val)
        if r_node == 0.0
            gy[i, j] = 0.0
        else
            g_mag = _eval_radial_gravity(r_node, r1, r_max, dr, Nr, g1, M_tot, g_bins)
            gy[i, j] = -g_mag * (dy_val / r_node)
        end
    end

    # Optional: populate gravitational potential FI at P nodes
    if FI !== nothing
        phi_rmax = -G * M_tot / r_max
        phi_bins = zeros(Float64, Nr)
        phi_bins[Nr] = phi_rmax
        for k in (Nr - 1):-1:1
            g_avg = 0.5 * (g_bins[k] + g_bins[k + 1])
            phi_bins[k] = phi_bins[k + 1] - g_avg * dr
        end
        phi0 = phi_bins[1] - 0.5 * g1 * r1

        @inbounds for j in 1:Nx1, i in 1:Ny1
            r_node = hypot(xp_val[j] - xc, yp_val[i] - yc)
            FI[i, j] = _eval_radial_potential(
                r_node, r1, r_max, dr, Nr, g1, phi0, M_tot, phi_bins
            )
        end
    end

    return r_bins, g_bins
end

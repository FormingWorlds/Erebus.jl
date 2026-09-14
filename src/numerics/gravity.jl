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
function assemble_gravitational_lse!(RHO, RP; coords=nothing, LP=nothing)
    Ny1, Nx1 = size(RHO)
    @unpack_coords coords dx dy xp yp
    xc_val = coords === nothing ? xcenter : coords.xcenter
    yc_val = coords === nothing ? ycenter : coords.ycenter
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
function assemble_gravitational_rhs!(RHO, RP; coords=nothing)
    Ny1, Nx1 = size(RHO)
    @unpack_coords coords xp yp
    xc_val = coords === nothing ? xcenter : coords.xcenter
    yc_val = coords === nothing ? ycenter : coords.ycenter
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
function process_gravitational_solution!(SP, FI, gx, gy; coords=nothing)
    Ny1, Nx1 = size(FI)
    Nx_val = Nx1 - 1
    Ny_val = Ny1 - 1
    @unpack_coords coords dx dy
    FI .= reshape(SP, Ny1, Nx1)
    @inbounds gx[:, 1:Nx_val] .= -diff(FI, dims=2) ./ dx_val
    @inbounds gy[1:Ny_val, :] .= -diff(FI, dims=1) ./ dy_val
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
function compute_gravity_solution!(SP, RP, RHO, FI, gx, gy; coords=nothing)
    LP = assemble_gravitational_lse!(RHO, RP; coords=coords)
    SP .= LP \ RP
    process_gravitational_solution!(SP, FI, gx, gy; coords=coords)
    return nothing
end # function compute_gravity_solution!

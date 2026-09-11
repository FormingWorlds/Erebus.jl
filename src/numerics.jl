
"""
Set up gravitational linear system of equations structures.

$(SIGNATURES)

# Details

    - nothing

# Returns 

    - RP: gravitational linear system of equations: RHS vector
    - SP: gravitational linear system of equations: solution vector
"""
function setup_gravitational_lse(Nx1::Int=Nx1, Ny1::Int=Ny1)
    RP = Vector{Float64}(undef, Ny1*Nx1)
    SP = Vector{Float64}(undef, Ny1*Nx1)
    return RP, SP
end
function setup_gravitational_lse(coords::GridCoordinates)
    return setup_gravitational_lse(coords.Nx1, coords.Ny1)
end

"""
Set up hydromechanical linear system of equations structures.

$(SIGNATURES)

# Details

    - nothing

# Returns 

    - R: hydromechanical linear system of equations: RHS vector
    - S: hydromechanical linear system of equations: solution vector
"""
function setup_hydromechanical_lse(Nx1::Int=Nx1, Ny1::Int=Ny1)
    R = Vector{Float64}(undef, Ny1*Nx1*6)
    S = Vector{Float64}(undef, Ny1*Nx1*6)
    return R, S
end
function setup_hydromechanical_lse(coords::GridCoordinates)
    return setup_hydromechanical_lse(coords.Nx1, coords.Ny1)
end

"""
Set up thermal linear system of equations structures.

$(SIGNATURES)

# Details

    - nothing

# Returns 

    - RT: thermal linear system of equations: RHS vector
    - ST: thermal linear system of equations: solution vector
"""
function setup_thermal_lse(Nx1::Int=Nx1, Ny1::Int=Ny1)
    RT = Vector{Float64}(undef, Ny1*Nx1)
    ST = Vector{Float64}(undef, Ny1*Nx1)
    return RT, ST
end
setup_thermal_lse(coords::GridCoordinates) = setup_thermal_lse(coords.Nx1, coords.Ny1)

"""
Initialize `iparm` parameters of Pardiso MKL solver.

$(SIGNATURES)

# Details

    - ps: Instance of pardiso solver
    - iparms_dict: dictionary of iparm parameters

# Returns

    - nothing
"""
function initialize_pardiso!(pardiso_solver, iparms_dict)
    set_msglvl!(pardiso_solver, Pardiso.MESSAGE_LEVEL_OFF)
    set_matrixtype!(pardiso_solver, Pardiso.REAL_NONSYM)
    set_nprocs!(pardiso_solver, cache_kwargs.nprocs)
    for (i, v) in iparms_dict
        set_iparm!(pardiso_solver, i+1, v)
    end
    return set_phase!(pardiso_solver, Pardiso.ANALYSIS)
end

"""
Compute viscosities, stresses, and density gradients
for hydromechanical solver.

$(SIGNATURES)

# Details

## In

    - ETA: viscoplastic viscosity at basic nodes
    - ETAP: viscosity at P nodes
    - GGG: shear modulus at basic nodes
    - GGGP: shear modulus at P nodes
    - SXY0: σ₀xy XY stress at basic nodes
    - SXX0:σ₀xy XY stress at basic nodes
    - RHOX: density at Vx nodes
    - RHOY: density at Vy nodes
    - dt: time step

## Out 

    - ETAcomp: computational viscosity at basic nodes
    - ETAPcomp: computational viscosity at P nodes
    - SXYcomp: previous XY stresses at basic nodes
    - SXXcomp: previous XX stresses at P nodes
    - SYYcomp: previous YY stresses at P nodes
    - dRHOXdx: density gradient at Vx nodes in x direction
    - dRHOXdy: density gradient at Vx nodes in y direction
    - dRHOYdx: density gradient at Vy nodes in x direction
    - dRHOYdy: density gradient at Vy nodes in y direction

# Returns

    - nothing
"""
function get_viscosities_stresses_density_gradients!(
    ETA,
    ETAP,
    GGG,
    GGGP,
    SXY0,
    SXX0,
    RHOX,
    RHOY,
    dt,
    ETAcomp,
    ETAPcomp,
    SXYcomp,
    SXXcomp,
    SYYcomp,
    dRHOXdx,
    dRHOXdy,
    dRHOYdx,
    dRHOYdy;
    coords=nothing,
)
    Ny1, Nx1 = size(RHOX)
    Nx_val = Nx1 - 1
    Ny_val = Ny1 - 1
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    # computational viscosity
    @views @. ETAcomp = ETA*GGG*dt / (GGG*dt + ETA)
    @views @. ETAPcomp = ETAP*GGGP*dt / (GGGP*dt + ETAP)
    # previous stresses
    @views @. SXYcomp = SXY0*ETA / (GGG*dt+ETA)
    @views @. SXXcomp = SXX0*ETAP / (GGGP*dt+ETAP)
    @views @. SYYcomp = -SXX0*ETAP / (GGGP*dt+ETAP)
    # for erroneously undersized (Ny, Nx) SSX0, SSX
    # @views @. SXXcomp = (
    # SXX0*ETAP[1:Ny, 1:Nx] / (GGGP[1:Ny, 1:Nx]*dt + ETAP[1:Ny, 1:Nx])
    # )
    # @views @. SYYcomp = (
    # -SXX0*ETAP[1:Ny, 1:Nx] / (GGGP[1:Ny, 1:Nx]*dt+ETAP[1:Ny, 1:Nx])
    # )
    # density gradients
    @inbounds begin
        @views @. dRHOXdx[:, 2:Nx_val] =
            0.5 * (RHOX[:, 3:Nx1]-RHOX[:, 1:(Nx1 - 2)]) * inv(dx_val)
        @views @. dRHOXdy[2:Ny_val, :] =
            0.5 * (RHOX[3:Ny1, :]-RHOX[1:(Ny1 - 2), :]) * inv(dy_val)
        @views @. dRHOYdx[:, 2:Nx_val] =
            0.5 * (RHOY[:, 3:Nx1]-RHOY[:, 1:(Nx1 - 2)]) * inv(dx_val)
        @views @. dRHOYdy[2:Ny_val, :] =
            0.5 * (RHOY[3:Ny1, :]-RHOY[1:(Ny1 - 2), :]) * inv(dy_val)
    end # @inbounds
    return nothing
end # function get_viscosities_stresses_density_gradients!

"""
Check if point (i, j) lies on or outside domain boundary for gravitational Poisson solver.

$(SIGNATURES)

# Details

    - i: y-index
    - j: x-index
    - Ny1: total number of P nodes in y
    - Nx1: total number of P nodes in x
    - xp_val: x coordinates of P nodes
    - yp_val: y coordinates of P nodes
    - xc_val: center coordinate in x
    - yc_val: center coordinate in y
    - r_limit: radius limit for circular planetary domain

# Returns

    - true if node is boundary or outside domain, false if interior
"""
@inline function is_gravitational_boundary(
    i::Integer,
    j::Integer,
    Ny1::Integer,
    Nx1::Integer,
    xp_val,
    yp_val,
    xc_val,
    yc_val,
    r_limit,
)
    return (
        i == 1 ||
        i == Ny1 ||
        j == 1 ||
        j == Nx1 ||
        distance(xp_val[j], yp_val[i], xc_val, yc_val) > r_limit
    )
end

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
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    xp_val = coords === nothing ? xp : coords.xp
    yp_val = coords === nothing ? yp : coords.yp
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
    xp_val = coords === nothing ? xp : coords.xp
    yp_val = coords === nothing ? yp : coords.yp
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
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
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

"""
Assemble hydromechanical system of equations.

$(SIGNATURES)

# Details

    - ETA: viscosity at basic nodes
    - ETAP: viscosity at P nodes
    - GGG: shear modulus at basic nodes
    - GGGP: shear modulus at P nodes
    - SXY0: previous XY stress at basic nodes
    - SXX0: previous XX stress at P nodes
    - RHOX: total density at Vx nodes
    - RHOY: total density at Vy nodes
    - RHOFX: fluid density at Vx nodes
    - RHOFY: fluid density at Vy nodes
    - RX: ηfluid/Kϕ at Vx nodes
    - RY: ηfluid/Kϕ at Vy nodes
    - ETAPHI: bulk viscosity at P nodes
    - BETAPHI: bulk compressibility at P nodes
    - PHI: porosity at P nodes
    - gx: x gravitational acceleration at Vx nodes
    - gy: y gravitational acceleration at Vy nodes
    - pr0: previous total pressure at P nodes
    - pf0: previous fluid pressure at P nodes
    - DMP: mass transfer term at P nodes
    - dt: time step
    - R: vector to store RHS coefficients
    - coords: grid coordinates
    - L: optional ExtendableSparseMatrix buffer to reuse

# Returns

    - L: LHS coefficient matrix (SparseMatrixCSC)
"""
function assemble_hydromechanical_lse!(
    ETA,
    ETAP,
    GGG,
    GGGP,
    SXY0,
    SXX0,
    RHOX,
    RHOY,
    RHOFX,
    RHOFY,
    RX,
    RY,
    ETAPHI,
    BETAPHI,
    PHI,
    gx,
    gy,
    pr0,
    pf0,
    DMP,
    dt,
    R;
    betasolid=betasolid,
    betafluid=betafluid,
    phimin=phimin,
    phimax=phimax,
    hydrofracture::Bool=false,
    pr=nothing,
    pf=nothing,
    TEN=nothing,
    KX=nothing,
    KY=nothing,
    kappa_frac::Real=1.0e3,
    gamma_frac::Real=1.0,
    k_frac_max::Real=1.0e-9,
    coords=nothing,
    L=nothing,
    venting::Bool=false,
    venting_mode::Symbol=:darcy_sink,
    k_vent::Real=1.0e-11,
    conductance_factor::Real=1.0,
    ice_sealing::Bool=false,
    t_freeze::Real=273.15,
    dt_seal::Real=10.0,
    k_seal_min_ratio::Real=1.0e-6,
    rplanet::Real=50000.0,
    xcenter::Real=70000.0,
    ycenter::Real=70000.0,
    P_amb::Real=10.0,
    venting_species::Symbol=:H2O,
    tk=nothing,
    eta_fluid_surf::Real=1.0e-3,
    L_sub::Real=2.83e6,
    S_vent_out=nothing,
)
    Ny1, Nx1 = size(ETAP)
    Nx_val = Nx1 - 1
    Ny_val = Ny1 - 1
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy

    # initialize or reuse LHS sparse coefficient matrix
    L = if L === nothing
        ExtendableSparseMatrix(Nx1 * Ny1 * 6, Nx1 * Ny1 * 6)
    else
        if !isempty(L.cscmatrix.nzval)
            nonzeros(L.cscmatrix) .= zero(0.0)
        end
        L
    end
    # reset RHS coefficient vector
    R .= 0.0
    @inbounds begin
        for j in 1:1:Nx1, i in 1:1:Ny1
            # define global indices in algebraic space
            kvx = ((j-1)*Ny1 + i-1) * 6 + 1 # Vx solid
            kvy = kvx + 1 # Vy solid
            kpm = kvx + 2 # P total
            kqx = kvx + 3 # qx Darcy
            kqy = kvx + 4 # qy Darcy
            kpf = kvx + 5 # P fluid
            # Vx equation
            if i==1 || i==Ny1 || j==1 || j==Nx_val || j==Nx1
                # Vx equation external points: boundary conditions
                # all locations: ghost unknowns Vx₃=0 -> 1.0⋅Vx[i,j]=0.0
                updateindex!(L, +, 1.0, kvx, kvx)
                # R[kvx] = 0.0 # already done with initialization
                # left boundary
                if j == 1
                    R[kvx] = vxleft
                end
                # right boundary
                if j == Nx_val
                    R[kvx] = vxright
                end
                # top boundary
                if i==1 && 1<j<Nx_val
                    updateindex!(L, +, bctop, kvx, kvx+6)
                end
                # bottom boundary
                if i==Ny1 && 1<j<Nx_val
                    updateindex!(L, +, bcbottom, kvx, kvx-6)
                end
            else
                # Vx momentum internal stencil (see docs/src/explanations/discretization_numerics.md)
                # computational viscosity
                ETA₁ =
                    ETA[i - 1, j] * GGG[i - 1, j] * dt / (GGG[i - 1, j]*dt + ETA[i - 1, j])
                ETA₂ = ETA[i, j] * GGG[i, j] * dt / (GGG[i, j]*dt + ETA[i, j])
                ETAP₁ = ETAP[i, j] * GGGP[i, j] * dt / (GGGP[i, j]*dt + ETAP[i, j])
                ETAP₂ =
                    ETAP[i, j + 1] * GGGP[i, j + 1] * dt /
                    (GGGP[i, j + 1]*dt + ETAP[i, j + 1])
                # previous stresses
                SXY₁ = SXY0[i - 1, j] * ETA[i - 1, j] / (GGG[i - 1, j]*dt + ETA[i - 1, j])
                SXY₂ = SXY0[i, j] * ETA[i, j] / (GGG[i, j]*dt + ETA[i, j])
                SXX₁ = SXX0[i, j] * ETAP[i, j] / (GGGP[i, j]*dt + ETAP[i, j])
                SXX₂ =
                    SXX0[i, j + 1] * ETAP[i, j + 1] / (GGGP[i, j + 1]*dt + ETAP[i, j + 1])
                # density gradients
                ∂RHO∂x = 0.5 * (RHOX[i, j + 1] - RHOX[i, j - 1]) * inv(dx_val)
                ∂RHO∂y = 0.5 * (RHOX[i + 1, j] - RHOX[i - 1, j]) * inv(dy_val)
                # LHS coefficient matrix
                updateindex!(L, +, ETAP₁/dx_val^2, kvx, kvx-6*Ny1) # Vx₁
                updateindex!(L, +, ETA₁/dy_val^2, kvx, kvx-6) # Vx₂
                updateindex!(
                    L,
                    +,
                    (
                        -(ETAP₁+ETAP₂) * inv(dx_val^2) - (ETA₁+ETA₂) * inv(dy_val^2) -
                        ∂RHO∂x * gx[i, j] * dt
                    ),
                    kvx,
                    kvx,
                ) # Vx₃
                updateindex!(L, +, ETA₂/dy_val^2, kvx, kvx+6) # Vx₄
                updateindex!(L, +, ETAP₂/dx_val^2, kvx, kvx+6*Ny1) # Vx₅
                updateindex!(
                    L,
                    +,
                    (
                        ETAP₁ * inv(dx_val) * inv(dy_val) -
                        ETA₂ * inv(dx_val) * inv(dy_val) - ∂RHO∂y * gx[i, j] * dt * 0.25
                    ),
                    kvx,
                    kvy,
                ) # Vy₂
                updateindex!(
                    L,
                    +,
                    (
                        -ETAP₂ * inv(dx_val) * inv(dy_val) +
                        ETA₂ * inv(dx_val) * inv(dy_val) - ∂RHO∂y * gx[i, j] * dt * 0.25
                    ),
                    kvx,
                    kvy+6*Ny1,
                ) # Vy₄
                updateindex!(
                    L,
                    +,
                    (
                        -ETAP₁ * inv(dx_val) * inv(dy_val) +
                        ETA₁ * inv(dx_val) * inv(dy_val) - ∂RHO∂y * gx[i, j] * dt * 0.25
                    ),
                    kvx,
                    kvy-6,
                ) # Vy₁
                updateindex!(
                    L,
                    +,
                    (
                        ETAP₂ * inv(dx_val) * inv(dy_val) -
                        ETA₁ * inv(dx_val) * inv(dy_val) - ∂RHO∂y * gx[i, j] * dt * 0.25
                    ),
                    kvx,
                    kvy+6*Ny1-6,
                ) # Vy₃
                updateindex!(L, +, Kcont*inv(dx_val), kvx, kpm) # P₁
                updateindex!(L, +, -Kcont*inv(dx_val), kvx, kpm+6*Ny1) # P₂
                # RHS coefficient vector
                R[kvx] = (
                    -RHOX[i, j] * gx[i, j] - (SXY₂-SXY₁) * inv(dy_val) -
                    (SXX₂-SXX₁) * inv(dx_val)
                )
            end # Vx equation
            # Vy equation
            if i==1 || i==Ny_val || i==Ny1 || j==1 || j==Nx1
                # Vy equation external points: boundary conditions
                # all locations: ghost unknowns Vy₃=0 -> 1.0⋅Vy[i,j]=0.0
                updateindex!(L, +, 1.0, kvy, kvy)
                # R[kvy] = 0.0 # already done with initialization
                # top boundary
                if i == 1
                    R[kvy] = vytop
                end
                # bottom boundary
                if i == Ny_val
                    R[kvy] = vybottom
                end
                # left boundary
                if j==1 && 1<i<Ny_val
                    updateindex!(L, +, bcleft, kvy, kvy+6*Ny1)
                end
                # right boundary
                if j==Nx1 && 1<i<Ny_val
                    updateindex!(L, +, bcright, kvy, kvy-6*Ny1)
                end
            else
                # Vy momentum internal stencil (see docs/src/explanations/discretization_numerics.md)
                # computational viscosity
                ETA₁ =
                    ETA[i, j - 1] * GGG[i, j - 1] * dt / (GGG[i, j - 1]*dt + ETA[i, j - 1])
                ETA₂ = ETA[i, j] * GGG[i, j] * dt / (GGG[i, j]*dt + ETA[i, j])
                ETAP₁ = ETAP[i, j] * GGGP[i, j] * dt / (GGGP[i, j]*dt + ETAP[i, j])
                ETAP₂ =
                    ETAP[i + 1, j] * GGGP[i + 1, j] * dt /
                    (GGGP[i + 1, j]*dt + ETAP[i + 1, j])
                # previous stresses
                SXY₁ = SXY0[i, j - 1] * ETA[i, j - 1] / (GGG[i, j - 1]*dt + ETA[i, j - 1])
                SXY₂ = SXY0[i, j] * ETA[i, j] / (GGG[i, j]*dt + ETA[i, j])
                SYY₁ = -SXX0[i, j] * ETAP[i, j] / (GGGP[i, j]*dt + ETAP[i, j])
                SYY₂ =
                    -SXX0[i + 1, j] * ETAP[i + 1, j] / (GGGP[i + 1, j]*dt + ETAP[i + 1, j])
                # density gradients
                ∂RHO∂x = 0.5 * (RHOY[i, j + 1]-RHOY[i, j - 1]) / dx_val
                ∂RHO∂y = 0.5 * (RHOY[i + 1, j]-RHOY[i - 1, j]) / dy_val
                # LHS coefficient matrix
                updateindex!(L, +, ETA₁/dx_val^2, kvy, kvy-6*Ny1) # Vy₁
                updateindex!(L, +, ETAP₁/dy_val^2, kvy, kvy-6) # Vy₂
                updateindex!(
                    L,
                    +,
                    (
                        -(ETAP₁+ETAP₂) * inv(dy_val^2) - (ETA₁+ETA₂) * inv(dx_val^2) -
                        ∂RHO∂y * gy[i, j] * dt
                    ),
                    kvy,
                    kvy,
                ) # Vy₃
                updateindex!(L, +, ETAP₂ * inv(dy_val^2), kvy, kvy+6) # Vy₄
                updateindex!(L, +, ETA₂ * inv(dx_val^2), kvy, kvy+6*Ny1) # Vy₅
                updateindex!(
                    L,
                    +,
                    (
                        ETAP₁ * inv(dx_val) * inv(dy_val) -
                        ETA₂ * inv(dx_val) * inv(dy_val) - ∂RHO∂x * gy[i, j] * dt * 0.25
                    ),
                    kvy,
                    kvx,
                ) # Vx₃
                updateindex!(
                    L,
                    +,
                    (
                        -ETAP₂ * inv(dx_val) * inv(dy_val) +
                        ETA₂ * inv(dx_val) * inv(dy_val) - ∂RHO∂x * gy[i, j] * dt * 0.25
                    ),
                    kvy,
                    kvx+6,
                ) # Vx₄
                updateindex!(
                    L,
                    +,
                    (
                        -ETAP₁ * inv(dx_val) * inv(dy_val) +
                        ETA₁ * inv(dx_val) * inv(dy_val) - ∂RHO∂x * gy[i, j] * dt * 0.25
                    ),
                    kvy,
                    kvx-6*Ny1,
                ) # Vx₁
                updateindex!(
                    L,
                    +,
                    (
                        ETAP₂ * inv(dx_val) * inv(dy_val) -
                        ETA₁ * inv(dx_val) * inv(dy_val) - ∂RHO∂x * gy[i, j] * dt * 0.25
                    ),
                    kvy,
                    kvx+6-6*Ny1,
                ) # Vx₂
                updateindex!(L, +, Kcont*inv(dy_val), kvy, kpm) # P₁
                updateindex!(L, +, -Kcont*inv(dy_val), kvy, kpm+6) # P₂
                R[kvy] = (
                    -RHOY[i, j] * gy[i, j] - (SXY₂-SXY₁) * inv(dx_val) -
                    (SYY₂-SYY₁) * inv(dy_val)
                ) # RHS
            end # Vy equation
            # P equation
            if i==1 || i==Ny1 || j==1 || j==Nx1
                # P equation external points: boundary conditions
                # all locations: ghost unknowns P=0 -> 1.0⋅P[i,j]=0.0
                updateindex!(L, +, 1.0, kpm, kpm)
                # R[kpm] = 0.0 # already done with initialization
                # elseif i==j==2
            elseif (
                (i==2 && 2<=j<=Nx_val) ||
                (j==2 && 2<i<Ny_val) ||
                (i==Ny_val && 2<=j<=Nx_val) ||
                (j==Nx_val && 2<i<Ny_val)
            )
                # Ptotal/Pfluid real pressure boundary condition 'anchor'
                updateindex!(L, +, Kcont, kpm, kpm)
                R[kpm] = psurface
            else
                # Solid continuity internal stencil: ∂Vx/∂x + ∂Vy/∂y = 0
                updateindex!(L, +, -1.0/dx_val, kpm, kvx-6*Ny1) # Vx₁
                updateindex!(L, +, 1.0/dx_val, kpm, kvx) # Vx₂
                updateindex!(L, +, -1.0/dy_val, kpm, kvy-6) # Vy₁
                updateindex!(L, +, 1.0/dy_val, kpm, kvy) # Vy₂
                # Poroelastic continuity stencils based on simple3anpfl.m (Taras Gerya, pers. comm.)
                betadrained = compute_drained_compressibility(
                    BETAPHI[i, j], PHI[i, j], betasolid; phimin=phimin, phimax=phimax
                )
                kbw = compute_biot_willis_coefficient(betadrained, betasolid)
                ksk = compute_skempton_coefficient(
                    betadrained,
                    PHI[i, j],
                    betasolid,
                    betafluid;
                    phimin=phimin,
                    phimax=phimax,
                )

                # LHS coefficient matrix
                updateindex!(
                    L,
                    +,
                    Kcont * (inv(ETAPHI[i, j]) / (1.0 - PHI[i, j]) + betadrained / dt),
                    kpm,
                    kpm,
                ) # P: Ptotal
                updateindex!(
                    L,
                    +,
                    -Kcont *
                    (inv(ETAPHI[i, j]) / (1.0 - PHI[i, j]) + betadrained * kbw / dt),
                    kpm,
                    kpf,
                ) # P: Pfluid
                # RHS coefficient vector
                R[kpm] = (betadrained * (pr0[i, j] - kbw * pf0[i, j]) / dt + DMP[i, j])
            end # P equation
            # qxDarcy equation
            if i==1 || i==Ny1 || j==1 || j==Nx_val || j==Nx1
                # qxDarcy equation external points: boundary conditions
                # all locations: ghost unknowns qyD = 0 -> 1.0⋅qxD[i, j] = 0.0
                updateindex!(L, +, 1.0, kqx, kqx)
                # R[kqx] = 0.0 # already done with initialization
                # top boundary
                if i==1 && 1<j<Nx_val
                    updateindex!(L, +, bcftop, kqx, kqx+6)
                end
                # bottom boundary
                if i==Ny1 && 1<j<Nx_val
                    updateindex!(L, +, bcfbottom, kqx, kqx-6)
                end
            else
                # x-Darcy flux internal stencil: (η_f/k_ϕx) * qxD + ∂P/∂x = ρ_f * gx
                # See Stencil Topologies in docs/src/explanations/discretization_numerics.md
                # LHS coefficient matrix
                rx_val = RX[i, j]
                if hydrofracture && pr !== nothing && pf !== nothing && TEN !== nothing
                    Peff_x = 0.5 * (pr[i, j] + pr[i, j + 1] - pf[i, j] - pf[i, j + 1])
                    sigma_t_x = 0.5 * (TEN[i, j] + TEN[i - 1, j])
                    kphi_x = (KX !== nothing) ? KX[i, j] : 0.0
                    if kphi_x > 0.0
                        keff_x = compute_hydrofracture_permeability(
                            kphi_x,
                            Peff_x,
                            sigma_t_x;
                            active=true,
                            kappa_frac=kappa_frac,
                            gamma=gamma_frac,
                            kmax=k_frac_max,
                        )
                        rx_val = RX[i, j] * (kphi_x / keff_x)
                    else
                        ffrac_x = compute_hydrofracture_factor(
                            Peff_x,
                            sigma_t_x;
                            active=true,
                            kappa_frac=kappa_frac,
                            gamma=gamma_frac,
                        )
                        rx_floor = 1.0e-5 / k_frac_max
                        rx_val = max(RX[i, j] / ffrac_x, rx_floor)
                    end
                end
                updateindex!(L, +, rx_val, kqx, kqx) # qxD
                updateindex!(L, +, -Kcont*inv(dx_val), kqx, kpf) # P₁
                updateindex!(L, +, Kcont*inv(dx_val), kqx, kpf+6*Ny1) # P₂
                # RHS coefficient vector
                R[kqx] = RHOFX[i, j] * gx[i, j]
            end # qxDarcy equation
            # qyDarcy equation
            if i==1 || i==Ny_val || i==Ny1 || j==1 || j==Nx1
                # qyDarcy equation external points: boundary conditions
                # all locations: ghost unknowns qyD = 0 -> 1.0⋅qyD[i, j] = 0.0
                updateindex!(L, +, 1.0, kqy, kqy)
                # R[kqy] = 0.0 # already done with initialization
                # left boundary
                if j==1 && 1<i<Ny_val
                    updateindex!(L, +, bcfleft, kqy, kqy+6*Ny1)
                end
                # right boundary
                if j==Nx1 && 1<i<Ny_val
                    updateindex!(L, +, bcfright, kqy, kqy-6*Ny1)
                end
            else
                # y-Darcy flux internal stencil: (η_f/k_ϕy) * qyD + ∂P/∂y = ρ_f * gy
                # LHS coefficient matrix
                ry_val = RY[i, j]
                if hydrofracture && pr !== nothing && pf !== nothing && TEN !== nothing
                    Peff_y = 0.5 * (pr[i, j] + pr[i + 1, j] - pf[i, j] - pf[i + 1, j])
                    sigma_t_y = 0.5 * (TEN[i, j] + TEN[i, j - 1])
                    kphi_y = (KY !== nothing) ? KY[i, j] : 0.0
                    if kphi_y > 0.0
                        keff_y = compute_hydrofracture_permeability(
                            kphi_y,
                            Peff_y,
                            sigma_t_y;
                            active=true,
                            kappa_frac=kappa_frac,
                            gamma=gamma_frac,
                            kmax=k_frac_max,
                        )
                        ry_val = RY[i, j] * (kphi_y / keff_y)
                    else
                        ffrac_y = compute_hydrofracture_factor(
                            Peff_y,
                            sigma_t_y;
                            active=true,
                            kappa_frac=kappa_frac,
                            gamma=gamma_frac,
                        )
                        ry_floor = 1.0e-5 / k_frac_max
                        ry_val = max(RY[i, j] / ffrac_y, ry_floor)
                    end
                end
                updateindex!(L, +, ry_val, kqy, kqy) # qyD
                updateindex!(L, +, -Kcont*inv(dy_val), kqy, kpf) # P₁
                updateindex!(L, +, Kcont*inv(dy_val), kqy, kpf+6) # P₂
                # RHS coefficient vector
                R[kqy] = RHOFY[i, j] * gy[i, j]
            end # qyDarcy equation
            # Ptotal/Pfluid equation 
            if i==1 || i==Ny1 || j==1 || j==Nx1
                # Ptotal/Pfluid equation external points: boundary conditions
                # all locations: ghost unknowns P = 0 -> 1.0⋅P[i, j] = 0.0
                updateindex!(L, +, 1.0, kpf, kpf)
                # R[kpf] = 0.0 # already done with initialization
                # elseif i==j==2
            elseif (
                (i==2 && 2<=j<=Nx_val) ||
                (j==2 && 2<i<Ny_val) ||
                (i==Ny_val && 2<=j<=Nx_val) ||
                (j==Nx_val && 2<i<Ny_val)
            )
                # Ptotal/Pfluid real pressure boundary condition 'anchor'
                updateindex!(L, +, Kcont, kpf, kpf)
                R[kpf] = psurface
            else
                # Fluid continuity internal stencil: ∂qxD/∂x + ∂qyD/∂y - (Pt - Pf)/η_ϕ = 0
                # LHS coefficient matrix
                updateindex!(L, +, -inv(dx_val), kpf, kqx-6*Ny1) # qxD₁
                updateindex!(L, +, inv(dx_val), kpf, kqx) # qxD₂
                updateindex!(L, +, -inv(dy_val), kpf, kqy-6) # qyD₁
                updateindex!(L, +, inv(dy_val), kpf, kqy) # qyD₂

                # LHS coefficient matrix
                updateindex!(
                    L,
                    +,
                    -Kcont *
                    (inv(ETAPHI[i, j]) / (1.0 - PHI[i, j]) + betadrained * kbw / dt),
                    kpf,
                    kpm,
                ) # Ptotal
                updateindex!(
                    L,
                    +,
                    Kcont *
                    (inv(ETAPHI[i, j]) / (1.0 - PHI[i, j]) + betadrained * kbw / ksk / dt),
                    kpf,
                    kpf,
                ) # Pfluid
                # RHS coefficient vector
                R[kpf] = -betadrained * kbw * (pr0[i, j] - (1.0 / ksk) * pf0[i, j]) / dt
            end # Ptotal/Pfluid equation
        end # for j=1:1:Nx1, i=1:1:Ny1
    end # @inbounds 

    if venting && tk !== nothing && coords !== nothing
        apply_venting_surface_boundary!(
            L,
            R,
            tk,
            coords,
            rplanet,
            xcenter,
            ycenter,
            P_amb;
            species=venting_species,
            k_vent=k_vent,
            conductance_factor=conductance_factor,
            mode=venting_mode,
            hydrofracture=hydrofracture,
            ice_sealing=ice_sealing,
            t_freeze=t_freeze,
            dt_seal=dt_seal,
            k_seal_min_ratio=k_seal_min_ratio,
            kappa_frac=kappa_frac,
            gamma_frac=gamma_frac,
            k_frac_max=k_frac_max,
            pr=pr,
            pf=pf,
            TEN=TEN,
            PHI=PHI,
            phimin=phimin,
            dt=dt,
            eta_fluid_surf=eta_fluid_surf,
            L_sub=L_sub,
            Kcont=Kcont,
            S_vent_out=S_vent_out,
        )
    end

    flush!(L) # finalize CSC matrix
    # return L
    return L.cscmatrix
end # function assemble_hydromechanical_lse!

"""
Process hydromechanical solution vector to output physical observables.

$(SIGNATURES)

# Details

    - S: hydromechanical solution vector
    - vx: solid velocity at Vx nodes
    - vy: solid velocity at Vy nodes
    - pr: total pressure at P nodes
    - qxD: qx-Darcy flux at Vx nodes
    - qyD: qy-Darcy flux at Vy nodes
    - pf: fluid pressure at P nodes

# Returns

    - nothing
"""
function process_hydromechanical_solution!(S, vx, vy, pr, qxD, qyD, pf; coords=nothing)
    Ny1, Nx1 = size(vx)
    S_mat = reshape(S, (:, Ny1, Nx1))
    @inbounds begin
        @views @. vx = S_mat[1, :, :]
        @views @. vy = S_mat[2, :, :]
        @views @. pr = S_mat[3, :, :] .* Kcont
        @views @. qxD = S_mat[4, :, :]
        @views @. qyD = S_mat[5, :, :]
        @views @. pf = S_mat[6, :, :] .* Kcont
    end # @inbounds
    # Δp = 0.25 * (pf[2, 2]+pf[2, Nx]+pf[Ny, 2]+pf[Ny, Nx]) - psurface
    # pr .-= Δp
    # pf .-= Δp
    return nothing
end # function process_hydromechanical_solution!

"""
Recompute bulk viscosity at P nodes.

# Details

    - ETA: viscoplastic viscosity at basic nodes
    - ETAP: viscosity at P nodes
    - ETAPHI: bulk viscosity at P nodes
    - PHI: porosity at P nodes
    - etaphikoef: coefficient: shear viscosity -> compaction viscosity

# Returns

    - nothing
"""
function recompute_bulk_viscosity!(ETA, ETAP, ETAPHI, PHI, etaphikoef)
    @inbounds begin
        @views @. ETAP[2:(end - 1), 2:(end - 1)] =
            4.0 / (
                inv(ETA[1:(end - 1), 1:(end - 1)]) +
                inv(ETA[2:end, 1:(end - 1)]) +
                inv(ETA[1:(end - 1), 2:end]) +
                inv(ETA[2:end, 2:end])
            )
        @views @. ETAPHI = etaphikoef * ETAP * inv(PHI)
    end # @inbounds
    return nothing
end

"""
Compute porosity coefficient Aϕ = Dln[(1-ϕ)/ϕ]/Dt

$(SIGNATURES)

# Details

## In

    - ETAPHI: bulk viscosity at P Nodes
    - BETAPHI: bulk compressibility at P nodes
    - PHI: porosity at P Nodes
    - pr: total pressure at P nodes
    - pf: fluid pressure at P nodes
    - pr0: previous step total pressure at P nodes
    - pf0: previous step fluid pressure at P nodes
    - dt: time step

## Out

    - APHI: porosity coefficient at P nodes

# Returns

    - aphimax: maximum absolute porosity coefficient
"""
function compute_Aϕ!(
    APHI,
    ETAPHI,
    BETAPHI,
    PHI,
    pr,
    pf,
    pr0,
    pf0,
    dt;
    coords=nothing,
    betasolid=betasolid,
    phimin=phimin,
    phimax=phimax,
    S_vent::Union{AbstractMatrix{Float64},Nothing}=nothing,
)
    # APHI .= 0.0
    Ny1, Nx1 = size(APHI)
    Nx = Nx1 - 1
    Ny = Ny1 - 1
    @inbounds begin
        for j in 2:Nx, i in 2:Ny
            betadrained = compute_drained_compressibility(
                BETAPHI[i, j], PHI[i, j], betasolid; phimin=phimin, phimax=phimax
            )
            kbw = compute_biot_willis_coefficient(betadrained, betasolid)
            compaction = (
                (pr[i, j] - pf[i, j]) / (ETAPHI[i, j] * (1.0 - PHI[i, j])) +
                betadrained * ((pr[i, j] - pr0[i, j]) - kbw * (pf[i, j] - pf0[i, j])) / dt
            )
            if S_vent !== nothing
                s_v = S_vent[i, j]
                if s_v > 0.0
                    # Subtract venting sink to avoid double drainage with explicit marker sink.
                    compaction = max(0.0, compaction - s_v)
                end
            end
            APHI[i, j] = compaction / PHI[i, j]
        end
        return maximum(abs, @view APHI[2:Ny, 2:Nx]) # includes [2, 2] anchor abberation
    end # @inbounds
    # return maximum(abs, APHI[3:Ny-1, 3:Nx-1]) # no abberation
end # function compute_Aϕ!

"""
Compute current fluid velocities.

$(SIGNATURES)

# Details

## In

    - PHIX: porosity at Vx nodes
    - PHIY: porosity at Vy nodes
    - qxD: qx-Darcy flux at Vx nodes
    - qyD: qy-Darcy flux at Vy nodes
    - vx: solid velocity at Vx nodes
    - vy: solid velocity at Vy nodes

## Out 

    - vxf: fluid vx velocity at Vx nodes
    - vyf: fluid vy velocity at Vy nodes

# Returns

    - nothing
"""
function compute_fluid_velocities!(PHIX, PHIY, qxD, qyD, vx, vy, vxf, vyf; coords=nothing)
    Ny1, Nx1 = size(vxf)
    Nx = Nx1 - 1
    Ny = Ny1 - 1
    @inbounds begin
        # vx velocity
        @views @. vxf[2:Ny, 1:Nx] = qxD[2:Ny, 1:Nx] / PHIX[2:Ny, 1:Nx]
        # top boundary
        @views @. vxf[1, :] = -bcftop * vxf[2, :]
        # bottom boundary
        @views @. vxf[Ny1, :] = -bcfbottom * vxf[Ny, :]
        # vy velocity
        @views @. vyf[1:Ny, 2:Nx] = qyD[1:Ny, 2:Nx] / PHIY[1:Ny, 2:Nx]
        # left boundary
        @views @. vyf[:, 1] = -bcfleft * vyf[:, 2]
        # right boundary
        @views @. vyf[:, Nx1] = -bcfright * vyf[:, Nx]
        # adding solid velocity
        @views @. vxf += vx
        @views @. vyf += vy
    end # @inbounds

    # for j=1:1:Nx, i=2:1:Ny
    #     vxf[i, j] = qxD[i, j]*inv(PHIX[i,j]) + vx[i, j]
    # end
    # @views @. vxf[1, :] = -bcftop*vxf[2, :]    
    # @views @. vxf[Ny1, :] = -bcfbottom*vxf[Ny, :]
    # for j=2:1:Nx, i=1:1:Ny
    #     vyf[i,j] = qyD[i,j]*inv(PHIY[i,j]) + vy[i,j]
    # end
    # @views @. vyf[:, 1] = -bcfleft*vyf[:, 2]    
    # @views @. vyf[:, Nx1] = -bcfright*vyf[:, Nx]     
    return nothing
end # function compute_fluid_velocities!

"""
Compute velocity-/displacement-limited time step.

$(SIGNATURES)

# Details

    - vx: solid vx velocity at Vx nodes
    - vy: solid vy velocity at Vy nodes
    - vxf: fluid vx velocity at Vx nodes
    - vyf: fluid vy velocity at Vy nodes
    - dt: current time step
    - aphimax: maximum observed porosity coefficient
   
# Returns

    - dt: displacement time step
"""
function compute_displacement_timestep(
    vx,
    vy,
    vxf,
    vyf,
    dt,
    aphimax;
    coords=nothing,
    dx_val=coords === nothing ? dx : coords.dx,
    dy_val=coords === nothing ? dy : coords.dy,
    dxymax_val=dxymax,
    dphimax_val=dphimax,
)
    maxvx = maximum(abs, vx)
    maxvy = maximum(abs, vy)
    maxvxf = maximum(abs, vxf)
    maxvyf = maximum(abs, vyf)
    @info "dt before velocity limitations = $dt s"
    dt = ifelse(dt*maxvx > dxymax_val*dx_val, dxymax_val*dx_val*inv(maxvx), dt)
    @info "dt after vx limitation = $dt s"
    dt = ifelse(dt*maxvy > dxymax_val*dy_val, dxymax_val*dy_val*inv(maxvy), dt)
    @info "dt after vy limitation = $dt s"
    dt = ifelse(dt*maxvxf > dxymax_val*dx_val, dxymax_val*dx_val*inv(maxvxf), dt)
    @info "dt after vxf limitation = $dt s"
    dt = ifelse(dt*maxvyf > dxymax_val*dy_val, dxymax_val*dy_val*inv(maxvyf), dt)
    @info "dt after vyf limitation = $dt s"
    dt = ifelse(dt*aphimax > dphimax_val, dphimax_val*inv(aphimax), dt)
    @info "dt after aphimax limitation = $dt s"
    return dt
end # function compute_displacement_timestep

"""
    compute_adaptive_timestep(
        vx, vy, vxf, vyf, dt, aphimax;
        coords=nothing,
        dx_val=coords === nothing ? dx : coords.dx,
        dy_val=coords === nothing ? dy : coords.dy,
        dxymax_val=dxymax,
        dphimax_val=dphimax,
        maxDTcurrent=0.0,
        DTmax_val=DTmax,
        dt_longest_val=dt_longest,
        dt_min=1.0,
    )

Compute multi-criterion adaptive timestep constrained by velocity CFL, porosity compaction,
thermal variation, and stability bounds.
"""
function compute_adaptive_timestep(
    vx,
    vy,
    vxf,
    vyf,
    dt,
    aphimax;
    coords=nothing,
    dx_val=coords === nothing ? dx : coords.dx,
    dy_val=coords === nothing ? dy : coords.dy,
    dxymax_val=dxymax,
    dphimax_val=dphimax,
    dt_ref=nothing,
    maxDTcurrent=0.0,
    DTmax_val=DTmax,
    dt_longest_val=dt_longest,
    dt_min=1.0,
    max_v_seg::Real=0.0,
    max_subcycles::Integer=2000,
    cfl_settling::Real=0.5,
)
    dt_cand = compute_displacement_timestep(
        vx,
        vy,
        vxf,
        vyf,
        dt,
        aphimax;
        coords=coords,
        dx_val=dx_val,
        dy_val=dy_val,
        dxymax_val=dxymax_val,
        dphimax_val=dphimax_val,
    )
    ref_dt = dt_ref === nothing ? dt : dt_ref
    if maxDTcurrent > DTmax_val && maxDTcurrent > 0.0
        dt_cand = min(dt_cand, ref_dt * (DTmax_val * inv(maxDTcurrent)))
    end
    if max_v_seg > 0.0
        min_dx = min(dx_val, dy_val)
        dt_cfl_seg = cfl_settling * min_dx / max_v_seg
        dt_seg_bound = max_subcycles * dt_cfl_seg
        dt_cand = min(dt_cand, dt_seg_bound)
    end
    dt_cand = clamp(dt_cand, dt_min, dt_longest_val)
    return dt_cand
end

"""
Compute stress, stress change, and strain rate components.

$(SIGNATURES)

# Details

## In

    - vx: solid vx velocity at Vx nodes
    - vy: solid vy velocity at Vy nodes
    - ETA: viscosity at basic nodes
    - GGG: shear modulus at basic nodes
    - ETAP: viscosity at P nodes
    - GGGP: shear modulus at P nodes
    - SXX0: previous time step σ₀′xx at P nodes
    - SXY0: previous time step σ₀xy at basic nodes
    - dt: computational time step

## Out

    - EXX: ϵxx at P nodes
    - EXY: ϵxy at basic nodes
    - SXX: σ′xx P nodes
    - SXY: σxy at basic nodes
    - DSXX: stress change Δσ′xx at P nodes
    - DSXY: stress change Δσxy at basic nodes
    - EII: second strain rate invariant ϵᴵᴵ at P nodes
    - SII: second stress invariant σᴵᴵ at P nodes

# Returns
    
        - nothing
"""
function compute_stress_strainrate!(
    vx,
    vy,
    ETA,
    GGG,
    ETAP,
    GGGP,
    SXX0,
    SXY0,
    EXX,
    EXY,
    SXX,
    SXY,
    DSXX,
    DSXY,
    EII,
    SII,
    dt;
    coords=nothing,
)
    Ny, Nx = size(EXY)
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    @inbounds begin
        # ϵxy, σxy, Δσxy at basic nodes
        for j in 1:1:Nx, i in 1:1:Ny
            EXY[i, j] =
                0.5 * ((vx[i + 1, j]-vx[i, j]) / dy_val + (vy[i, j + 1]-vy[i, j])/dx_val)
            SXY[i, j] = (
                2*ETA[i, j]*EXY[i, j]*GGG[i, j]*dt / (GGG[i, j]*dt+ETA[i, j]) +
                SXY0[i, j]*ETA[i, j] / (GGG[i, j]*dt+ETA[i, j])
            )
            DSXY[i, j] = SXY[i, j] - SXY0[i, j]
        end
        # ϵxx, σ′xx, Δσ'xx and Eᴵᴵ, Sᴵᴵ at P nodes
        for j in 2:1:Nx, i in 2:1:Ny
            EXX[i, j] =
                0.5 * ((vx[i, j]-vx[i, j - 1]) / dx_val - (vy[i, j]-vy[i - 1, j]) / dy_val)
            SXX[i, j] = (
                2*ETAP[i, j]*EXX[i, j]*GGGP[i, j]*dt / (GGGP[i, j]*dt+ETAP[i, j]) +
                SXX0[i, j]*ETAP[i, j] / (GGGP[i, j]*dt+ETAP[i, j])
            )
            DSXX[i, j] = SXX[i, j] - SXX0[i, j]
            EII[i, j] = sqrt(EXX[i, j]^2 + grid_average(i-1, j-1, EXY)^2)
            SII[i, j] = sqrt(SXX[i, j]^2 + grid_average(i-1, j-1, SXY)^2)
        end
    end # @inbounds        
    # # ϵxy, σxy, Δσxy at basic nodes
    # EXY .= 0.5.*(diff(vx, dims=1)[:, 1:Nx]./dy .+ diff(vy, dims=2)[1:Ny, :]./dx)
    # @. SXY = 2*ETA*EXY*GGG*dt/(GGG*dt+ETA) + SXY0*ETA/(GGG*dt+ETA)
    # @. DSXY = SXY - SXY0
    # # ϵxx, σ′xx at P nodes
    # # @. DIVV[2:end, 2:end] = 
    # #     diff(vx, dims=2)[2:end, :]/dx + diff(vy, dims=1)[:, 2:end]/dy
    # EXX[2:Ny1, 2:Nx1] .= ( 
    #     0.5 .* (
    #         diff(vx, dims=2)[2:Ny1, :]./dx .- diff(vy, dims=1)[:, 2:Nx1]./dy
    #     )
    # )
    # @. SXX = 2.0*ETAP*EXX*GGGP*dt/(GGGP*dt+ETAP) + SXX0*ETAP/(GGGP*dt+ETAP)
    # @. DSXX = SXX - SXX0
    # @. EII[2:Ny, 2:Nx] = sqrt(
    #     EXX[2:Ny, 2:Nx]^2 + (
    #         (
    #             EXY[2:Ny, 2:Nx]
    #             +EXY[1:Ny-1,2:Nx]
    #             +EXY[2:Ny,1:Nx-1]
    #             +EXY[1:Ny-1,1:Nx-1]
    #         )/4.0
    #     )^2
    # )
    # @. SII[2:Ny, 2:Nx] = sqrt(
    #     SXX[2:Ny, 2:Nx]^2 + (
    #         (
    #             SXY[2:Ny, 2:Nx]
    #             +SXY[1:Ny-1,2:Nx]
    #             +SXY[2:Ny,1:Nx-1]
    #             +SXY[1:Ny-1,1:Nx-1]
    #         )/4.0
    #     )^2
    # )
    return nothing
end # function compute_stress_strainrate!

"""
Apply symmetry to P node observables.

$(SIGNATURES)

# Details

    - SXX: σ′xx at P nodes
    - APHI: Aϕ = Dln[(1-ϕ)/ϕ]/Dt at P nodes
    - PHI: porosity at P nodes
    - pr: total pressure at P nodes
    - pf: fluid pressure at P nodes
    - ps: solid pressure at P nodes

# Returns

    - nothing
"""
function symmetrize_p_node_observables!(SXX, APHI, PHI, pr, pf, ps)
    Ny1, Nx1 = size(SXX)
    Nx = Nx1 - 1
    Ny = Ny1 - 1
    # top boundary
    @inbounds @views @. begin
        SXX[1, 2:Nx] = SXX[2, 2:Nx]
        APHI[1, 2:Nx] = APHI[2, 2:Nx]
        PHI[1, 2:Nx] = PHI[2, 2:Nx]
        pr[1, 2:Nx] = pr[2, 2:Nx]
        pf[1, 2:Nx] = pf[2, 2:Nx]
        # bottom boundary
        SXX[Ny1, 2:Nx] = SXX[Ny, 2:Nx]
        APHI[Ny1, 2:Nx] = APHI[Ny, 2:Nx]
        PHI[Ny1, 2:Nx] = PHI[Ny, 2:Nx]
        pr[Ny1, 2:Nx] = pr[Ny, 2:Nx]
        pf[Ny1, 2:Nx] = pf[Ny, 2:Nx]
        # left boundary
        SXX[:, 1] = SXX[:, 2]
        APHI[:, 1] = APHI[:, 2]
        PHI[:, 1] = PHI[:, 2]
        pr[:, 1] = pr[:, 2]
        pf[:, 1] = pf[:, 2]
        # right boundary
        SXX[:, Nx1] = SXX[:, Nx]
        APHI[:, Nx1] = APHI[:, Nx]
        PHI[:, Nx1] = PHI[:, Nx]
        pr[:, Nx1] = pr[:, Nx]
        pf[:, Nx1] = pf[:, Nx]
        # solid pressure
        ps = (pr-pf*PHI) * inv(1-PHI)
    end
    return nothing
end # function symmetrize_p_node_observables!

"""
Compute nodal adjustment and return plastic iterations completeness status.

$(SIGNATURES)

# Details

    - ETA: viscoplastic viscosity at basic nodes
    - ETA0: previous time step viscoplastic viscosity at basic nodes
    - ETA5: plastic iterations viscoplastic viscosity at basic nodes
    - GGG: shear modulus at basic nodes
    - SXX: σ′xx at P nodes
    - SXY: σxy at basic nodes
    - pr: total pressure at P nodes
    - pf: fluid pressure at P nodes
    - COH: compressive strength at basic nodes 
    - TEN: tensile strength at basic nodes 
    - FRI: friction at basic nodes
    - YNY: plastic yielding status at basic nodes 
    - YNY5: plastic iterations plastic yielding status at basic nodes
    - YERRNOD: vector of summed yielding errors of nodes over plastic iterations
    - DSY: (SIIB-syield) at basic nodes
    - dt: time set
    - iplast: plastic iteration step 

# Returns

    - plastic_iterations_complete: true if plastic iterations complete
"""
function compute_nodal_adjustment!(
    ETA,
    ETA0,
    ETA5,
    GGG,
    SXX,
    SXY,
    pr,
    pf,
    COH,
    TEN,
    FRI,
    YNY,
    YNY5,
    YERRNOD,
    DSY,
    dt,
    iplast,
)
    # reset / setup
    Ny, Nx = size(ETA)
    ETA5 .= ETA0
    YNY5 .= 0
    DSY .= 0.0
    ynpl = 0
    ddd = 0.0
    @inbounds begin
        for j in 1:1:Nx, i in 1:1:Ny
            # second stress invariant at basic nodes
            SIIB = sqrt(SXY[i, j]^2 + grid_average(i, j, SXX)^2)
            # second invariant for purely elastic stress buildup at basic nodes
            siiel = SIIB * (GGG[i, j]*dt+ETA[i, j]) / ETA[i, j]
            # interpolate total and fluid pressure at basic nodes
            prB = grid_average(i, j, pr)
            pfB = grid_average(i, j, pf)
            # Yielding stress: confined and tensile fracture.
            # Note: uses Terzaghi effective stress (prB - pfB) for frictional shear and tensile failure,
            # following standard rock mechanics (Terzaghi 1943, Handin et al. 1963).
            syieldc = COH[i, j] + FRI[i, j] * (prB-pfB)
            syieldt = TEN[i, j] + (prB-pfB)
            # non-negative yielding stress requirement
            syield = max(min(syieldc, syieldt), 0.0)
            # update error for previous yielding nodes
            ynn = false
            if YNY[i, j] > 0
                ynn = true
                DSY[i, j] = SIIB - syield
                ddd += DSY[i, j]^2
                ynpl += 1
            end
            # correcting viscosity for yielding
            if syield < siiel
                # update viscosity for basic node
                etapl = dt * GGG[i, j] * syield/(siiel-syield)
                if etapl < ETA0[i, j]
                    # recompute nodal viscosity, apply min/max viscosity cutoffs
                    ETA5[i, j] = etapl^(1.0-etawt) * ETA[i, j]^etawt
                    if ETA5[i, j] > etamax
                        ETA5[i, j] = etamax
                    elseif ETA5[i, j] < etamin
                        ETA5[i, j] = etamin
                    end
                    # mark yielding nodes
                    YNY5[i, j] = 1
                    # update error for new yielding nodes
                    if ynn == false
                        DSY[i, j] = SIIB - syield
                        ddd += DSY[i, j]^2
                        ynpl += 1
                    end
                end
            end
        end
        if ynpl > 0
            YERRNOD[iplast] = sqrt(ddd/ynpl)
        end
        # return plastic iteration completeness
        @info "end plastic iter $iplast: ynpl=$ynpl, YERRNOD=$(YERRNOD[iplast])"
        return ynpl==0 || YERRNOD[iplast]<yerrmax || iplast==nplast
    end # @inbounds
end # function compute_nodal_adjustment!

"""
Compare two arrays of identical sizes element-wise and fill a third array with the larger value if positive and zero otherwise.

$(SIGNATURES)

# Details

    - A: first array
    - B: second array
    - C: result array

# Returns

    - nothing
"""
function positive_max!(A, B, C)
    @inbounds for i in eachindex(A)
        C[i] = max(0, ifelse(A[i] > B[i], A[i], B[i]))
    end
    return nothing
end # function positive_max

"""
Decide next pass plastic iteration time step, viscoplastic viscosity,
and basic node yielding status.

$(SIGNATURES)

# Details:

    - ETA: viscoplastic viscosity at basic nodes
    - ETA5: plastic iterations viscoplastic viscosity at basic nodes
    - ETA00: previous time step viscoplastic viscosity at basic nodes
    - YNY: plastic yielding status at basic nodes 
    - YNY5: plastic iterations plastic yielding status at basic nodes
    - YNY00: previous time step plastic yielding status at basic nodes
    - YNY_inv_ETA: inverse of plastic viscosity at yielding basic nodes
    - dt: current time step
    - iplast: current plastic iteration counter

# Returns

    - dt: adjusted next time step
"""
function finalize_plastic_iteration_pass!(
    ETA, ETA5, ETA00, YNY, YNY5, YNY00, YNY_inv_ETA, dt, iplast
)
    if iplast % dtstep == 0
        # dtstep plastic iterations performed without reaching targets:
        # decrease time step and reset to previous viscoplastic viscosity
        dt *= dtcoefdn
        @info "reducing dt due to plastic iteration limit: dt=$dt s"
        ETA .= ETA00
        YNY .= YNY00
    else
        # perform next plastic iteration pass with new viscoplastic viscosity
        ETA .= ETA5
        YNY .= YNY5
    end
    @views @. YNY_inv_ETA = YNY / ETA
    return dt
end # function finalize_plastic_iteration_pass

"""
Assemble the LHS sparse coefficient matrix and fill RHS coefficient vector
of the energy conservation (heat) equation.

$(SIGNATURES)

# Details

	- tk1: current temperature at P nodes
	- RHOCP: volumetric heat capacity at P nodes  
	- KX: thermal conductivity at Vx nodes
	- KY: thermal conductivity at Vy nodes 
	- HR: radioactive heating at P nodes
	- HA: adiabatic heating at P nodes 
	- HS: shear heating at P nodes
    - DHP: latent heating (HL) at P nodes
    - RT: thermal RHS coefficient vector
    - dt: current time step length
    - coords: grid coordinates
    - LT: optional ExtendableSparseMatrix buffer to reuse

# Returns

    - LT: LHS sparse coefficient matrix
"""
function assemble_thermal_lse!(
    tk1,
    RHOCP,
    KX,
    KY,
    HR,
    HA,
    HS,
    DHP,
    RT,
    dt;
    coords=nothing,
    LT=nothing,
    Q_metric=nothing,
    Q_lat=nothing,
    Q_seg=nothing,
)
    Ny1, Nx1 = size(tk1)
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    # fresh or reusable LHS coefficient matrix
    LT = if LT === nothing
        ExtendableSparseMatrix(Ny1 * Nx1, Ny1 * Nx1)
    else
        if !isempty(LT.cscmatrix.nzval)
            nonzeros(LT.cscmatrix) .= zero(0.0)
        end
        LT
    end
    # reset RHS coefficient vector
    RT .= zero(0.0)
    # compose global thermal matrix LT and coefficient vector RT
    @inbounds begin
        for j in 1:1:Nx1, i in 1:1:Ny1
            # define global index in algebraic space
            gk = (j-1)*Ny1 + i
            # External points
            if i==1 || i==Ny1 || j==1 || j==Nx1
                # thermal equation external points: boundary conditions
                # all locations: ghost unknowns T₃=0 -> 1.0⋅T[i,j]=0.0
                updateindex!(LT, +, 1.0, gk, gk)
                # R[gk] = 0.0 # already done with initialization
                # left boundary: ∂T/∂x=0
                if j == 1
                    updateindex!(LT, +, -1.0, gk, gk+Ny1)
                end
                # right boundary: ∂T/∂x=0
                if j == Nx1
                    updateindex!(LT, +, -1.0, gk, gk-Ny1)
                end
                # top inner boundary: ∂T/∂y=0
                if i==1 && 1<j<Nx1
                    updateindex!(LT, +, -1.0, gk, gk+1)
                end
                # bottom inner boundary: ∂T/∂y=0
                if i==Ny1 && 1<j<Nx1
                    updateindex!(LT, +, -1.0, gk, gk-1)
                end
            else
                # internal points: 2D thermal equation (conservative formulation)
                # extract thermal conductivities
                Kx₁ = KX[i, j - 1]
                Kx₂ = KX[i, j]
                Ky₁ = KY[i - 1, j]
                Ky₂ = KY[i, j]
                # fill system of equations: LHS
                updateindex!(LT, +, -Kx₁*inv(dx_val^2), gk, gk-Ny1) # T₁
                updateindex!(LT, +, -Ky₁*inv(dy_val^2), gk, gk-1) # T₂
                updateindex!(
                    LT,
                    +,
                    (RHOCP[i, j]/dt + (Kx₁+Kx₂)*inv(dx_val^2) + (Ky₁+Ky₂)*inv(dy_val^2)),
                    gk,
                    gk,
                ) # T₃
                updateindex!(LT, +, -Ky₂*inv(dy_val^2), gk, gk+1) # T₄
                updateindex!(LT, +, -Kx₂*inv(dx_val^2), gk, gk+Ny1) # T₅
                # fill system of equations: RHS
                RT[gk] = (
                    RHOCP[i, j]/dt*tk1[i, j] + HR[i, j] + HA[i, j] + HS[i, j] + DHP[i, j]
                )
                if Q_metric !== nothing
                    RT[gk] += Q_metric[i, j]
                end
                if Q_lat !== nothing
                    RT[gk] += Q_lat[i, j]
                end
                if Q_seg !== nothing
                    RT[gk] += Q_seg[i, j]
                end
            end
        end
    end # @inbounds

    flush!(LT) # finalize CSC matrix
    return LT
end # function assemble_thermal_lse!

"""
    apply_radiative_surface_boundary!(
        KX, KY, tk, coords, rplanet, xcenter, ycenter, T_amb;
        emissivity=0.9, sigma_sb=5.670374419e-8,
        k_rock=nothing, marker_property_mode=1
    )

Apply linearized Stefan-Boltzmann radiative cooling at planetesimal-sticky air
interface faces.

Modifies interface conductivities in KX and KY using harmonic series resistance:

    1 / U_eff = (Δ / (2 * k_bulk)) + (1 / h_rad)
    k_face = U_eff * Δ = 2 * k_bulk * (h_rad * Δ) / (2 * k_bulk + h_rad * Δ)

where `h_rad = compute_radiation_htc(T_surf, T_amb; emissivity, sigma_sb)`, and `k_bulk`
is local rock thermal conductivity: evaluated dynamically via
`compute_ksolidm(T_surf, marker_property_mode)` when `k_rock === nothing`, or given by
constant scalar `k_rock` when specified.

# Arguments
- `KX`: Horizontal conductivity array at Vx nodes [W/(m K)]
- `KY`: Vertical conductivity array at Vy nodes [W/(m K)]
- `tk`: Temperature field [K] at P nodes
- `coords`: Grid coordinate descriptors
- `rplanet`: Planetesimal radius [m]
- `xcenter`: Center horizontal coordinate [m]
- `ycenter`: Center vertical coordinate [m]
- `T_amb`: Ambient disk temperature [K]
- `emissivity`: Radiative emissivity
- `sigma_sb`: Stefan-Boltzmann constant
- `k_rock`: Rock thermal conductivity [W/(m K)] override
  (default: `nothing` for temperature-dependent calculation)
- `marker_property_mode`: Marker property regime index (default: 1)
- `tau_LW`: Infrared optical depth of atmosphere for greenhouse blanketing (default: 0.0)
"""
function apply_radiative_surface_boundary!(
    KX::AbstractMatrix{Float64},
    KY::AbstractMatrix{Float64},
    tk::AbstractMatrix{Float64},
    coords::GridCoordinates,
    rplanet::Real,
    xcenter::Real,
    ycenter::Real,
    T_amb::Real;
    emissivity::Real=0.9,
    sigma_sb::Real=5.670374419e-8,
    k_rock::Union{Real,Nothing}=nothing,
    marker_property_mode::Int=1,
    phi::Real=0.0,
    kfluid::Real=50.0,
    tau_LW::Real=0.0,
)
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    dx = coords.dx
    dy = coords.dy
    rplanet2 = rplanet^2

    # Horizontal faces KX[i, j] between P(i, j) and P(i, j+1)
    @inbounds for j in 1:(Nx1 - 1)
        xj1 = coords.xp[j] - xcenter
        xj2 = coords.xp[j + 1] - xcenter
        for i in 1:Ny1
            yi = coords.yp[i] - ycenter
            r1_sq = xj1^2 + yi^2
            r2_sq = xj2^2 + yi^2
            is_rock1 = r1_sq <= rplanet2
            is_rock2 = r2_sq <= rplanet2
            if is_rock1 != is_rock2
                T_surf = is_rock1 ? tk[i, j] : tk[i, j + 1]
                h_rad = compute_effective_radiation_htc(
                    T_surf, T_amb, tau_LW; emissivity=emissivity, sigma_sb=sigma_sb
                )
                k_rad = h_rad * dx
                k_bulk = if k_rock !== nothing
                    Float64(k_rock)
                elseif phi > 0.0
                    ktotal(
                        compute_ksolidm(T_surf, marker_property_mode),
                        Float64(kfluid),
                        Float64(phi),
                    )
                else
                    compute_ksolidm(T_surf, marker_property_mode)
                end
                KX[i, j] = (2.0 * k_bulk * k_rad) / (2.0 * k_bulk + k_rad)
            end
        end
    end

    # Vertical faces KY[i, j] between P(i, j) and P(i+1, j)
    @inbounds for j in 1:Nx1
        xj = coords.xp[j] - xcenter
        for i in 1:(Ny1 - 1)
            yi1 = coords.yp[i] - ycenter
            yi2 = coords.yp[i + 1] - ycenter
            r1_sq = xj^2 + yi1^2
            r2_sq = xj^2 + yi2^2
            is_rock1 = r1_sq <= rplanet2
            is_rock2 = r2_sq <= rplanet2
            if is_rock1 != is_rock2
                T_surf = is_rock1 ? tk[i, j] : tk[i + 1, j]
                h_rad = compute_effective_radiation_htc(
                    T_surf, T_amb, tau_LW; emissivity=emissivity, sigma_sb=sigma_sb
                )
                k_rad = h_rad * dy
                k_bulk = if k_rock !== nothing
                    Float64(k_rock)
                elseif phi > 0.0
                    ktotal(
                        compute_ksolidm(T_surf, marker_property_mode),
                        Float64(kfluid),
                        Float64(phi),
                    )
                else
                    compute_ksolidm(T_surf, marker_property_mode)
                end
                KY[i, j] = (2.0 * k_bulk * k_rad) / (2.0 * k_bulk + k_rad)
            end
        end
    end
    return nothing
end

"""
    compute_face_venting_permeability(
        k_v, breached, ice_sealing, T_surf, peff, sigma_t;
        t_freeze=273.15, dt_seal=10.0, k_seal_min_ratio=1.0e-6,
        kappa_frac=1.0e3, gamma_frac=1.0, k_frac_max=1.0e-9,
    )

Compute effective rock face permeability at the planetesimal surface, accounting for
tensile hydrofracture breaching or cryogenic pore ice sealing.

# Arguments
- `k_v::Real`: Reference matrix permeability [m^2].
- `breached::Bool`: Whether tensile failure has ruptured the rock lid.
- `ice_sealing::Bool`: Whether cryogenic pore ice sealing is active below freezing.
- `T_surf::Real`: Local surface rock temperature [K].
- `peff::Real`: Terzaghi effective stress `P_t - P_f` [Pa].
- `sigma_t::Real`: Rock tensile strength [Pa].

# Keywords
- `t_freeze::Real`: Water freezing temperature [K] (default: 273.15).
- `dt_seal::Real`: Temperature sealing interval [K] (default: 10.0).
- `k_seal_min_ratio::Real`: Minimum residual cryogenic permeability ratio (default: 1.0e-6).
- `kappa_frac::Real`: Hydrofracture multiplier (default: 1.0e3).
- `gamma_frac::Real`: Hydrofracture power-law exponent (default: 1.0).
- `k_frac_max::Real`: Maximum fractured permeability ceiling [m^2] (default: 1.0e-9).

# Returns
- Effective face permeability [m^2].
"""
function compute_face_venting_permeability(
    k_v::Real,
    breached::Bool,
    ice_sealing::Bool,
    T_surf::Real,
    peff::Real,
    sigma_t::Real;
    species::Symbol=:H2O,
    t_freeze::Real=273.15,
    dt_seal::Real=10.0,
    k_seal_min_ratio::Real=1.0e-6,
    kappa_frac::Real=1.0e3,
    gamma_frac::Real=1.0,
    k_frac_max::Real=1.0e-9,
)
    if breached
        return compute_hydrofracture_permeability(
            k_v,
            peff,
            sigma_t;
            active=true,
            kappa_frac=kappa_frac,
            gamma=gamma_frac,
            kmax=k_frac_max,
        )
    elseif ice_sealing && (species === :H2O || species === :water)
        return compute_ice_sealed_permeability(
            k_v,
            T_surf;
            T_freeze=t_freeze,
            delta_T_seal=dt_seal,
            k_min_ratio=k_seal_min_ratio,
        )
    else
        return Float64(k_v)
    end
end

"""
    apply_venting_surface_boundary!(
        L, R, tk, coords, rplanet, xcenter, ycenter, P_amb;
        k_vent=1.0e-11, conductance_factor=1.0, mode=:darcy_sink,
        hydrofracture=false,
        ice_sealing=false, t_freeze=273.15, dt_seal=10.0, k_seal_min_ratio=1.0e-6,
        kappa_frac=1.0e3, gamma_frac=1.0, k_frac_max=1.0e-9,
        pr=nothing, pf=nothing, TEN=nothing, PHI=nothing, phimin=1.0e-4, dt=1.0e10,
        eta_fluid_surf=1.0e-3, L_sub=2.83e6, Kcont=1.0e20, S_vent_out=nothing
    )

Apply permeable venting sink boundary condition at rock-air interface faces (`r = rplanet`).
Computes local venting pressure `P_vent = max(P_amb, P_sat,ice(T_surf))` and assembles Robin
conductance into the fluid continuity row (scaled by `Kcont`).

If `mode === :hydrofracture_gated`, venting requires `pr`, `pf`, and `TEN` to evaluate tensile
failure `Peff <= -sigma_t`. If any pressure array is missing or the lid is unbreached, the face
remains closed (`is_open = false`).
If `ice_sealing === true`, sub-freezing rock faces (`T_surf < t_freeze`) experience exponential
pore ice permeability sealing during unbreached porous flow (`:darcy_sink` mode).
When overpressure breaches the lid (`Peff <= -sigma_t` and `hydrofracture === true` or
`mode === :hydrofracture_gated`), enhanced hydrofracture permeability opens.
Only outward venting is permitted (`pf > P_vent`), and venting is fluid-limited (`phi > phimin`).
"""
function apply_venting_surface_boundary!(
    L,
    R::Union{AbstractVector{Float64},Nothing},
    tk::AbstractMatrix{Float64},
    coords::GridCoordinates,
    rplanet::Real,
    xcenter::Real,
    ycenter::Real,
    P_amb::Real;
    species::Symbol=:H2O,
    k_vent::Real=1.0e-11,
    conductance_factor::Real=1.0,
    mode::Symbol=:darcy_sink,
    hydrofracture::Bool=false,
    ice_sealing::Bool=false,
    t_freeze::Real=273.15,
    dt_seal::Real=10.0,
    k_seal_min_ratio::Real=1.0e-6,
    kappa_frac::Real=1.0e3,
    gamma_frac::Real=1.0,
    k_frac_max::Real=1.0e-9,
    pr::Union{AbstractMatrix{Float64},Nothing}=nothing,
    pf::Union{AbstractMatrix{Float64},Nothing}=nothing,
    TEN::Union{AbstractMatrix{Float64},Nothing}=nothing,
    PHI::Union{AbstractMatrix{Float64},Nothing}=nothing,
    phimin::Real=1.0e-4,
    dt::Real=1.0e10,
    eta_fluid_surf::Real=1.0e-3,
    L_sub::Real=2.83e6,
    Kcont::Real=1.0e20,
    S_vent_out::Union{AbstractMatrix{Float64},Nothing}=nothing,
)
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    dx = coords.dx
    dy = coords.dy
    if dx <= 0.0 || dy <= 0.0 || !isfinite(dx) || !isfinite(dy)
        throw(DomainError((dx, dy), "Grid spacing must be > 0 and finite"))
    end
    eta_f = Float64(eta_fluid_surf)
    if eta_f <= 0.0 || !isfinite(eta_f)
        throw(DomainError(eta_f, "eta_fluid_surf must be > 0 and finite"))
    end
    rplanet2 = Float64(rplanet)^2
    p_amb_val = Float64(P_amb)
    k_v = Float64(k_vent)
    c_factor = Float64(conductance_factor)
    kcont_val = Float64(Kcont)
    dt_val = max(Float64(dt), 1.0e-12)
    phimin_val = Float64(phimin)

    if S_vent_out !== nothing
        S_vent_out .= 0.0
    end

    # Horizontal faces between P(i, j) and P(i, j+1)
    @inbounds for j in 1:(Nx1 - 1)
        xj1 = coords.xp[j] - xcenter
        xj2 = coords.xp[j + 1] - xcenter
        for i in 1:Ny1
            yi = coords.yp[i] - ycenter
            r1_sq = xj1^2 + yi^2
            r2_sq = xj2^2 + yi^2
            is_rock1 = r1_sq <= rplanet2
            is_rock2 = r2_sq <= rplanet2
            if is_rock1 != is_rock2
                i_rock = i
                j_rock = is_rock1 ? j : (j + 1)

                # Skip domain boundary ghost and anchor nodes
                if i_rock < 2 || i_rock > Ny1 - 1 || j_rock < 2 || j_rock > Nx1 - 1
                    continue
                end

                breached = false
                peff = 0.0
                sigma_t = 0.0
                if (mode === :hydrofracture_gated || hydrofracture) &&
                    pr !== nothing &&
                    pf !== nothing &&
                    TEN !== nothing
                    peff = pr[i_rock, j_rock] - pf[i_rock, j_rock]
                    sigma_t = TEN[i_rock, j_rock]
                    breached = is_hydrofracture_breached(peff, sigma_t)
                end

                is_open = true
                if mode === :hydrofracture_gated && !breached
                    is_open = false
                end

                T_raw = tk[i_rock, j_rock]
                T_surf = isfinite(T_raw) ? max(T_raw, 1.0e-3) : 1.0e-3
                P_vent = compute_venting_pressure(
                    T_surf, p_amb_val; species=species, L_sub=L_sub
                )

                # One-sided venting condition: pore fluid must exceed venting pressure
                if pf !== nothing
                    pf_cur = pf[i_rock, j_rock]
                    if pf_cur <= P_vent
                        is_open = false
                    end
                end

                # Fluid-availability limit: no venting from dry rock
                phi_avail = 1.0
                if PHI !== nothing
                    phi_rock = PHI[i_rock, j_rock]
                    phi_avail = max(0.0, phi_rock - phimin_val)
                    if phi_avail <= 0.0
                        is_open = false
                    end
                end

                if is_open
                    k_face = compute_face_venting_permeability(
                        k_v,
                        breached,
                        ice_sealing,
                        T_surf,
                        peff,
                        sigma_t;
                        species=species,
                        t_freeze=t_freeze,
                        dt_seal=dt_seal,
                        k_seal_min_ratio=k_seal_min_ratio,
                        kappa_frac=kappa_frac,
                        gamma_frac=gamma_frac,
                        k_frac_max=k_frac_max,
                    )

                    C_face = (k_face / (eta_f * dx^2)) * c_factor
                    kpf = ((j_rock - 1) * Ny1 + i_rock - 1) * 6 + 6

                    # If pf and PHI are known, check and cap Darcy rate by available fluid
                    C_face_eff = C_face
                    S_vent_actual = 0.0
                    if pf !== nothing
                        pf_cur = pf[i_rock, j_rock]
                        S_darcy = C_face * (pf_cur - P_vent)
                        if PHI !== nothing
                            S_max = phi_avail / dt_val
                            if S_darcy > S_max && S_darcy > 0.0
                                C_face_eff = C_face * (S_max / S_darcy)
                                S_vent_actual = S_max
                            else
                                S_vent_actual = max(0.0, S_darcy)
                            end
                        else
                            S_vent_actual = max(0.0, S_darcy)
                        end
                    end

                    if L !== nothing && R !== nothing
                        updateindex!(L, +, kcont_val * C_face_eff, kpf, kpf)
                        R[kpf] += C_face_eff * P_vent
                    end

                    if S_vent_out !== nothing
                        S_vent_out[i_rock, j_rock] += S_vent_actual
                    end
                end
            end
        end
    end

    # Vertical faces between P(i, j) and P(i+1, j)
    @inbounds for j in 1:Nx1
        xj = coords.xp[j] - xcenter
        for i in 1:(Ny1 - 1)
            yi1 = coords.yp[i] - ycenter
            yi2 = coords.yp[i + 1] - ycenter
            r1_sq = xj^2 + yi1^2
            r2_sq = xj^2 + yi2^2
            is_rock1 = r1_sq <= rplanet2
            is_rock2 = r2_sq <= rplanet2
            if is_rock1 != is_rock2
                i_rock = is_rock1 ? i : (i + 1)
                j_rock = j

                # Skip domain boundary ghost and anchor nodes
                if i_rock < 2 || i_rock > Ny1 - 1 || j_rock < 2 || j_rock > Nx1 - 1
                    continue
                end

                breached = false
                peff = 0.0
                sigma_t = 0.0
                if (mode === :hydrofracture_gated || hydrofracture) &&
                    pr !== nothing &&
                    pf !== nothing &&
                    TEN !== nothing
                    peff = pr[i_rock, j_rock] - pf[i_rock, j_rock]
                    sigma_t = TEN[i_rock, j_rock]
                    breached = is_hydrofracture_breached(peff, sigma_t)
                end

                is_open = true
                if mode === :hydrofracture_gated && !breached
                    is_open = false
                end

                T_raw = tk[i_rock, j_rock]
                T_surf = isfinite(T_raw) ? max(T_raw, 1.0e-3) : 1.0e-3
                P_vent = compute_venting_pressure(
                    T_surf, p_amb_val; species=species, L_sub=L_sub
                )

                # One-sided venting condition: pore fluid must exceed venting pressure
                if pf !== nothing
                    pf_cur = pf[i_rock, j_rock]
                    if pf_cur <= P_vent
                        is_open = false
                    end
                end

                # Fluid-availability limit: no venting from dry rock
                phi_avail = 1.0
                if PHI !== nothing
                    phi_rock = PHI[i_rock, j_rock]
                    phi_avail = max(0.0, phi_rock - phimin_val)
                    if phi_avail <= 0.0
                        is_open = false
                    end
                end

                if is_open
                    k_face = compute_face_venting_permeability(
                        k_v,
                        breached,
                        ice_sealing,
                        T_surf,
                        peff,
                        sigma_t;
                        species=species,
                        t_freeze=t_freeze,
                        dt_seal=dt_seal,
                        k_seal_min_ratio=k_seal_min_ratio,
                        kappa_frac=kappa_frac,
                        gamma_frac=gamma_frac,
                        k_frac_max=k_frac_max,
                    )

                    C_face = (k_face / (eta_f * dy^2)) * c_factor
                    kpf = ((j_rock - 1) * Ny1 + i_rock - 1) * 6 + 6

                    # If pf and PHI are known, check and cap Darcy rate by available fluid
                    C_face_eff = C_face
                    S_vent_actual = 0.0
                    if pf !== nothing
                        pf_cur = pf[i_rock, j_rock]
                        S_darcy = C_face * (pf_cur - P_vent)
                        if PHI !== nothing
                            S_max = phi_avail / dt_val
                            if S_darcy > S_max && S_darcy > 0.0
                                C_face_eff = C_face * (S_max / S_darcy)
                                S_vent_actual = S_max
                            else
                                S_vent_actual = max(0.0, S_darcy)
                            end
                        else
                            S_vent_actual = max(0.0, S_darcy)
                        end
                    end

                    if L !== nothing && R !== nothing
                        updateindex!(L, +, kcont_val * C_face_eff, kpf, kpf)
                        R[kpf] += C_face_eff * P_vent
                    end

                    if S_vent_out !== nothing
                        S_vent_out[i_rock, j_rock] += S_vent_actual
                    end
                end
            end
        end
    end
    return nothing
end

"""
Perform thermal iterations to time step thermal field at P nodes.

$(SIGNATURES)

# Details

    - tk0: previous temperature at P nodes 
	- tk1: current temperature at P nodes
	- tk2: next temperature at P nodes 
	- DT: calculated temperature difference at P nodes 
	- DT0: previous calculated temperature difference at P nodes
	- RHOCP: volumetric heat capacity at P nodes  
	- KX: thermal conductivity at Vx nodes
	- KY: thermal conductivity at Vy nodes 
	- HR: radioactive heating at P nodes
	- HA: adiabatic heating at P nodes 
	- HS: shear heating at P nodes
    - DHP: latent heating (HL) at P nodes
    - RT: thermal RHS coefficient vector
    - ST: thermal solution vector
	- dt: computational time step

# Returns

    - nothing
"""
function perform_thermal_iterations!(
    tk0,
    tk1,
    tk2,
    DT,
    DT0,
    RHOCP,
    KX,
    KY,
    HR,
    HA,
    HS,
    DHP,
    RT,
    ST,
    dt;
    coords=nothing,
    Q_metric=nothing,
    Q_lat=nothing,
    DTmax::Real=DTmax,
)
    # set up thermal iterations
    Ny1, Nx1 = size(tk1)
    tk0 .= tk1
    dtt = dt
    dttsum = 0.0
    titer = 1
    # perform thermal iterations until reaching time limit
    while dttsum < dt
        # fresh LHS coefficient matrix
        LT = assemble_thermal_lse!(
            tk1,
            RHOCP,
            KX,
            KY,
            HR,
            HA,
            HS,
            DHP,
            RT,
            dtt;
            coords=coords,
            Q_metric=Q_metric,
            Q_lat=Q_lat,
        )

        # solve system of equations
        ST .= LT \ RT # implicit: flush!(LT)
        # reshape solution vector to 2D array
        tk2 .= reshape(ST, Ny1, Nx1)
        # compute ΔT
        DT .= tk2 .- tk1
        if titer == 1
            # during first thermal iteration pass:
            # apply thermal timestepping stability condition
            maxDTcurrent = maximum(abs, DT)
            if maxDTcurrent > DTmax
                dtt *= DTmax * inv(maxDTcurrent)
            else
                dttsum += dtt
                tk1 .= tk2
            end
        else
            # second+ thermal iteration passes:
            # update dttsum and adjust timestep
            dttsum += dtt
            tk1 .= tk2
            dtt = min(dtt, dt - dttsum)
        end
        # increase thermal iteration counter
        titer += 1
    end
    # finalize overall temperature change and advance temperature field
    DT .= tk2 .- tk0
    DT0 .= DT
    return nothing
end # function perform_thermal_iterations!

"""
Decide next pass thermochemical iteration time step.

$(SIGNATURES)

# Details:

    - maxDTcurrent: maximum temperature difference between current and
                    previous time step
    - dt: current time step duration
    - titer: current thermochemical iteration counter

# Returns

    - dt: adjusted next time step
"""
function finalize_thermochemical_iteration_pass(maxDTcurrent, dt, titer)
    if titer == 1
        if maxDTcurrent > DTmax
            dt *= (DTmax * inv(maxDTcurrent))
            @info "titer 1: reducing dt due to maxDT: dt=$dt s"
        end
    end
    return dt
end # function finalize_thermochemical_iteration_pass

"""
Assess outcome of thermochemical iteration and return thermochemical iterations
completeness status.

$(SIGNATURES)

# Details:

    - DMP: mass transfer term at P nodes
    - pf: fluid pressure at P nodes
    - pf0: previous time step fluid pressure at P nodes
    - titer: current thermochemical iteration counter

# Returns

    - dt: adjusted next time step
"""
function compute_thermochemical_iteration_outcome(DMP, pf, pf0, titer; pferrmax=1.0e5)
    pferrcur = maximum(abs, pf - pf0)
    DMPmax = maximum(abs, DMP)
    @info "end thermochemical iter $titer" pferrcur DMPmax
    return pferrcur < pferrmax && (titer > 2 || DMPmax <= 0.0)
end # function compute_thermochemical_iteration_outcome

"""
Apply iron metal segregation via sub-cycled conservative drift-flux transport.

$(SIGNATURES)

Solves the conservative drift-flux transport equation for molten iron metal
percolation through a solid silicate matrix and Stokes settling through a magma ocean.
Subcycles the explicit finite-volume transport step using a local CFL criterion.
Guarantees mass conservation of total metal to machine precision when input marker
metal fractions satisfy 0 <= Xfe_bulk <= phi_pack.

References:
- Stevenson (1990), Fluid dynamics of core formation.
- Deguen et al. (2014), Earth Planet. Sci. Lett., 391, 274-287.
- Lichtenberg et al. (2019, 2021), Science / JGR Planets.

# Arguments
- `xm::AbstractVector{Float64}`: Marker x-coordinates [m]
- `ym::AbstractVector{Float64}`: Marker y-coordinates [m]
- `tm::AbstractVector{<:Integer}`: Marker material phase type
- `tkm::AbstractVector{Float64}`: Marker temperature [K]
- `phim::AbstractVector{Float64}`: Marker silicate melt fraction / porosity [-]
- `Xfe_bulk::AbstractVector{Float64}`: Marker bulk metal volume fraction [-]
- `Xfem::AbstractVector{Float64}`: Marker molten metal fraction [-]
- `marknum::Integer`: Number of markers
- `dt::Real`: Timestep duration [s]
- `cfg_core::CoreFormationConfig`: Core formation configuration parameters

# Keyword Arguments
- `coords=nothing`: `GridCoordinates` domain geometry struct
- `xcenter::Real=coords !== nothing ? coords.xcenter : 70000.0`: Planet center x [m]
- `ycenter::Real=coords !== nothing ? coords.ycenter : 70000.0`: Planet center y [m]
- `rplanet::Real=50000.0`: Planet radius [m]
- `g_surf::Real=0.1`: Reference surface gravity magnitude [m/s^2]
- `gx::Union{Nothing,AbstractMatrix{Float64}}=nothing`: Optional x-gravity on grid [m/s^2]
- `gy::Union{Nothing,AbstractMatrix{Float64}}=nothing`: Optional y-gravity on grid [m/s^2]
- `Q_seg_grid::Union{Nothing,AbstractMatrix{Float64}}=nothing`: Optional grid to accumulate dissipation heating [W/m^3]
- `rho_silicate::Real=3000.0`: Reference silicate rock density [kg/m^3]
- `eta_silicate::Real=1.0e18`: Reference silicate rock dynamic viscosity [Pa s]
- `ETA::Union{Nothing,AbstractMatrix{Float64}}=nothing`: Matrix viscosity on grid [Pa s]
- `Fm::Union{Nothing,AbstractVector{Float64}}=nothing`: Marker silicate melt fraction [-]
- `T_solidus_silicate::Real=1400.0`: Reference silicate solidus temperature [K]
- `T_liquidus_silicate::Real=1800.0`: Reference silicate liquidus temperature [K]

# Returns
- NamedTuple `(; max_v_seg, n_subcycles, dt_sub, total_dissipation_energy)` where `total_dissipation_energy` is in [J/m].
"""
function apply_metal_segregation!(
    xm::AbstractVector{Float64},
    ym::AbstractVector{Float64},
    tm::AbstractVector{<:Integer},
    tkm::AbstractVector{Float64},
    phim::AbstractVector{Float64},
    Xfe_bulk::AbstractVector{Float64},
    Xfem::AbstractVector{Float64},
    marknum::Integer,
    dt::Real,
    cfg_core::CoreFormationConfig;
    coords=nothing,
    xcenter::Real=coords !== nothing ? coords.xcenter : 70000.0,
    ycenter::Real=coords !== nothing ? coords.ycenter : 70000.0,
    rplanet::Real=50000.0,
    g_surf::Real=0.1,
    gx::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    gy::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    Q_seg_grid::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    rho_silicate::Real=3000.0,
    eta_silicate::Real=1.0e18,
    ETA::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    Fm::Union{Nothing,AbstractVector{Float64}}=nothing,
    T_solidus_silicate::Real=1400.0,
    T_liquidus_silicate::Real=1800.0,
    Xfe_H_m::Union{Nothing,AbstractVector{Float64}}=nothing,
    Xfe_C_m::Union{Nothing,AbstractVector{Float64}}=nothing,
    Xfe_N_m::Union{Nothing,AbstractVector{Float64}}=nothing,
    Xfe_S_m::Union{Nothing,AbstractVector{Float64}}=nothing,
    cfg_partition::Union{Nothing,MetalPartitionConfig}=nothing,
)
    if (!cfg_core.percolation_active && !cfg_core.settling_active) ||
        dt <= 0.0 ||
        marknum <= 0
        return (; max_v_seg=0.0, n_subcycles=0, dt_sub=0.0, total_dissipation_energy=0.0)
    end

    track_volatiles = (
        cfg_partition !== nothing &&
        cfg_partition.active &&
        Xfe_H_m !== nothing &&
        Xfe_C_m !== nothing &&
        Xfe_N_m !== nothing &&
        Xfe_S_m !== nothing
    )

    # Validate input marker bounds
    @inbounds for m in 1:marknum
        if tm[m] < 3
            rmark = distance(xm[m], ym[m], xcenter, ycenter)
            if rmark <= rplanet
                xfe = Xfe_bulk[m]
                if !isfinite(xfe) || xfe < 0.0 || xfe > cfg_core.phi_pack + 1.0e-7
                    throw(
                        DomainError(
                            xfe,
                            "Marker bulk metal fraction must be finite, non-negative, and <= phi_pack",
                        ),
                    )
                end
            end
        end
    end

    Nx_val = coords !== nothing ? coords.Nx : 32
    Ny_val = coords !== nothing ? coords.Ny : 32
    dx_val = if coords !== nothing
        coords.dx
    else
        (coords !== nothing ? coords.xsize / Nx_val : 4375.0)
    end
    dy_val = if coords !== nothing
        coords.dy
    else
        (coords !== nothing ? coords.ysize / Ny_val : 4375.0)
    end

    # Allocate cell accumulations
    M_fe_cell = zeros(Float64, Ny_val, Nx_val)
    M_rock_markers = zeros(Int, Ny_val, Nx_val)
    v_seg_cell = zeros(Float64, Ny_val, Nx_val)
    phi_m_cell = zeros(Float64, Ny_val, Nx_val)
    F_m_cell = zeros(Float64, Ny_val, Nx_val)
    Xfem_cell = zeros(Float64, Ny_val, Nx_val)
    g_acc_cell = zeros(Float64, Ny_val, Nx_val)
    cap_cell = zeros(Float64, Ny_val, Nx_val)
    phi_fe_cell = zeros(Float64, Ny_val, Nx_val)
    T_cell = zeros(Float64, Ny_val, Nx_val)
    drho_cell = zeros(Float64, Ny_val, Nx_val)

    M_fe_H_cell = track_volatiles ? zeros(Float64, Ny_val, Nx_val) : zeros(Float64, 0, 0)
    M_fe_C_cell = track_volatiles ? zeros(Float64, Ny_val, Nx_val) : zeros(Float64, 0, 0)
    M_fe_N_cell = track_volatiles ? zeros(Float64, Ny_val, Nx_val) : zeros(Float64, 0, 0)
    M_fe_S_cell = track_volatiles ? zeros(Float64, Ny_val, Nx_val) : zeros(Float64, 0, 0)

    # Bin markers into grid cells
    @inbounds for m in 1:marknum
        if tm[m] < 3
            rmark = distance(xm[m], ym[m], xcenter, ycenter)
            if rmark <= rplanet
                j_c = clamp(Int(floor(xm[m] / dx_val)) + 1, 1, Nx_val)
                i_c = clamp(Int(floor(ym[m] / dy_val)) + 1, 1, Ny_val)
                fe_m = Xfe_bulk[m]
                M_fe_cell[i_c, j_c] += fe_m
                M_rock_markers[i_c, j_c] += 1
                phi_m_cell[i_c, j_c] += Xfem[m]
                T_cell[i_c, j_c] += tkm[m]
                F_m_val = if Fm !== nothing
                    Fm[m]
                elseif tkm[m] >= T_solidus_silicate &&
                    T_liquidus_silicate > T_solidus_silicate
                    clamp(
                        (tkm[m] - T_solidus_silicate) /
                        (T_liquidus_silicate - T_solidus_silicate),
                        0.0,
                        1.0,
                    )
                else
                    0.0
                end
                F_m_cell[i_c, j_c] += F_m_val
                cap_cell[i_c, j_c] += max(cfg_core.phi_pack - fe_m, 0.0)

                if track_volatiles
                    M_fe_H_cell[i_c, j_c] += fe_m * Xfe_H_m[m]
                    M_fe_C_cell[i_c, j_c] += fe_m * Xfe_C_m[m]
                    M_fe_N_cell[i_c, j_c] += fe_m * Xfe_N_m[m]
                    M_fe_S_cell[i_c, j_c] += fe_m * Xfe_S_m[m]
                end
            end
        end
    end

    # Normalize cell averages
    @inbounds for j in 1:Nx_val, i in 1:Ny_val
        n_m = M_rock_markers[i, j]
        if n_m > 0
            phi_fe_cell[i, j] = M_fe_cell[i, j] / n_m
            phi_m_cell[i, j] /= n_m
            F_m_cell[i, j] /= n_m
            T_cell[i, j] /= n_m
            m_bulk = phi_fe_cell[i, j]
            Xfem_cell[i, j] =
                m_bulk > 0.0 ? clamp(phi_m_cell[i, j] / m_bulk, 0.0, 1.0) : 0.0
        end
    end

    # Compute cell segregation velocities
    @inbounds for j in 1:Nx_val, i in 1:Ny_val
        n_m = M_rock_markers[i, j]
        if n_m == 0
            continue
        end
        xc = (j - 0.5) * dx_val
        yc = (i - 0.5) * dy_val
        rc = distance(xc, yc, xcenter, ycenter)
        if rc > rplanet
            continue
        end

        g_acc = if gx !== nothing && gy !== nothing && i <= size(gx, 1) && j <= size(gx, 2)
            g_mag = sqrt(gx[i, j]^2 + gy[i, j]^2)
            g_mag > 0.0 ? g_mag : g_surf * min(rc / rplanet, 1.0)
        else
            g_surf * min(rc / rplanet, 1.0)
        end
        g_acc_cell[i, j] = g_acc

        phi_m = phi_m_cell[i, j]
        F_m = F_m_cell[i, j]

        rho_metal_eff = if cfg_core.metal_density_mode !== :constant
            w_S_val =
                if track_volatiles &&
                    cfg_partition.dynamic_sulfur_density &&
                    M_fe_cell[i, j] > 0.0
                    clamp((M_fe_S_cell[i, j] / M_fe_cell[i, j]) * 1.0e-6, 0.0, 0.40)
                else
                    cfg_core.sulfur_fraction
                end
            compute_liquid_metal_density(
                w_S_val; T=max(T_cell[i, j], 100.0), law=cfg_core.metal_density_mode
            )
        else
            cfg_core.rho_metal
        end
        drho = max(rho_metal_eff - rho_silicate, 1.0)
        drho_cell[i, j] = drho

        if phi_m > 0.0 && g_acc > 0.0 && drho > 0.0
            eta_matrix = if ETA !== nothing && i <= size(ETA, 1) && j <= size(ETA, 2)
                ETA[i, j]
            else
                eta_silicate
            end
            eta_susp = compute_melt_weakened_viscosity(
                eta_matrix, F_m, 1; phi_crit=0.4, eta_melt=10.0, etamin=0.1, etamax=1.0e20
            )
            r_drop = if cfg_core.droplet_size_mode === :fixed
                cfg_core.droplet_diameter_fixed / 2.0
            elseif cfg_core.droplet_size_mode === :capillary_mean ||
                cfg_core.droplet_size_mode === :bond_mean ||
                cfg_core.droplet_size_mode === :weber_mean
                # Gravity-capillary (Bond) balance: d = sqrt(We_crit * sigma / (drho * g))
                d_cap = sqrt(
                    cfg_core.We_crit * cfg_core.sigma_metal_silicate /
                    max(drho * g_acc, 1.0e-8),
                )
                clamp(d_cap / 2.0, 1.0e-4, 5.0e-2)
            else # :weber_turbulent
                v_est = stokes_settling_velocity(
                    cfg_core.droplet_diameter_fixed / 2.0,
                    drho,
                    max(g_acc, 1.0e-5),
                    eta_susp,
                )
                v_rel = max(v_est, 1.0e-6)
                d_weber = weber_equilibrium_diameter(
                    rho_silicate,
                    v_rel,
                    cfg_core.sigma_metal_silicate;
                    We_crit=cfg_core.We_crit,
                )
                clamp(d_weber / 2.0, 1.0e-4, 5.0e-2)
            end

            v_seg_cell[i, j] = metal_segregation_velocity(
                phi_m,
                F_m,
                drho,
                g_acc,
                eta_susp;
                percolation_active=cfg_core.percolation_active,
                settling_active=cfg_core.settling_active,
                k_metal_ref=cfg_core.k_metal_ref,
                eta_metal=cfg_core.eta_metal,
                phi_crit_perc=cfg_core.phi_crit_perc,
                phi_residual=cfg_core.phi_residual,
                phi0=cfg_core.phi0,
                perm_exponent=cfg_core.perm_exponent,
                r_drop=r_drop,
                hindered_exponent=cfg_core.hindered_exponent,
                phi_pack=cfg_core.phi_pack,
                hadamard_rybczynski=cfg_core.hadamard_rybczynski,
                F_settle_start=cfg_core.F_settle_start,
                F_perc_end=cfg_core.F_perc_end,
            )
        end
    end

    max_v = maximum(v_seg_cell)
    if max_v <= 0.0
        return (; max_v_seg=0.0, n_subcycles=0, dt_sub=0.0, total_dissipation_energy=0.0)
    end

    # CFL calculation and subcycling
    dt_cfl = cfg_core.cfl_settling * min(dx_val, dy_val) / max_v
    n_sub_raw = Int(ceil(dt / dt_cfl))
    if n_sub_raw > cfg_core.max_subcycles
        @warn "CFL subcycling requires $n_sub_raw steps, capped at max_subcycles $(cfg_core.max_subcycles); segregation transport throttled" maxlog=10
    end
    n_sub = clamp(n_sub_raw, 1, cfg_core.max_subcycles)
    dt_sub = dt / n_sub

    # Working copy of cell metal mass for subcycling
    m_fe = copy(M_fe_cell)
    total_diss_energy = 0.0

    m_fe_H = track_volatiles ? copy(M_fe_H_cell) : zeros(Float64, 0, 0)
    m_fe_C = track_volatiles ? copy(M_fe_C_cell) : zeros(Float64, 0, 0)
    m_fe_N = track_volatiles ? copy(M_fe_N_cell) : zeros(Float64, 0, 0)
    m_fe_S = track_volatiles ? copy(M_fe_S_cell) : zeros(Float64, 0, 0)

    flux_H_x = track_volatiles ? zeros(Float64, Ny_val, Nx_val - 1) : zeros(Float64, 0, 0)
    flux_H_y = track_volatiles ? zeros(Float64, Ny_val - 1, Nx_val) : zeros(Float64, 0, 0)
    flux_C_x = track_volatiles ? zeros(Float64, Ny_val, Nx_val - 1) : zeros(Float64, 0, 0)
    flux_C_y = track_volatiles ? zeros(Float64, Ny_val - 1, Nx_val) : zeros(Float64, 0, 0)
    flux_N_x = track_volatiles ? zeros(Float64, Ny_val, Nx_val - 1) : zeros(Float64, 0, 0)
    flux_N_y = track_volatiles ? zeros(Float64, Ny_val - 1, Nx_val) : zeros(Float64, 0, 0)
    flux_S_x = track_volatiles ? zeros(Float64, Ny_val, Nx_val - 1) : zeros(Float64, 0, 0)
    flux_S_y = track_volatiles ? zeros(Float64, Ny_val - 1, Nx_val) : zeros(Float64, 0, 0)

    # Pre-allocated arrays for subcycling fluxes and limiters
    req_flux_x = zeros(Float64, Ny_val, Nx_val - 1)
    req_flux_y = zeros(Float64, Ny_val - 1, Nx_val)
    flux_x = zeros(Float64, Ny_val, Nx_val - 1)
    flux_y = zeros(Float64, Ny_val - 1, Nx_val)
    outflow_tot = zeros(Float64, Ny_val, Nx_val)
    inflow_tot = zeros(Float64, Ny_val, Nx_val)
    alpha_out = ones(Float64, Ny_val, Nx_val)
    alpha_in = ones(Float64, Ny_val, Nx_val)

    # Subcycling loop
    for _ in 1:n_sub
        fill!(outflow_tot, 0.0)
        fill!(inflow_tot, 0.0)
        fill!(req_flux_x, 0.0)
        fill!(req_flux_y, 0.0)

        # 1. Compute unscaled requested fluxes across East-West faces
        @inbounds for j in 1:(Nx_val - 1)
            xf = j * dx_val
            for i in 1:Ny_val
                yf = (i - 0.5) * dy_val
                dxf = xf - xcenter
                dyf = yf - ycenter
                rf = sqrt(dxf^2 + dyf^2)
                if rf > rplanet || rf < 1.0e-3
                    continue
                end

                # Determine transport direction
                nx =
                    if gx !== nothing &&
                        gy !== nothing &&
                        i <= size(gx, 1) &&
                        j <= size(gx, 2)
                        gx_f = gx[i, j]
                        gy_f =
                            0.5 *
                            (gy[i, j] + (j + 1 <= size(gy, 2) ? gy[i, j + 1] : gy[i, j]))
                        g_f = sqrt(gx_f^2 + gy_f^2)
                        g_f > 1.0e-10 ? gx_f / g_f : -dxf / rf
                    else
                        -dxf / rf
                    end

                vf = 0.5 * (v_seg_cell[i, j] + v_seg_cell[i, j + 1])
                uf = vf * nx

                donor_j = uf > 0.0 ? j : j + 1
                rec_j = uf > 0.0 ? j + 1 : j

                n_donor = M_rock_markers[i, donor_j]
                n_rec = M_rock_markers[i, rec_j]
                if n_donor == 0 || n_rec == 0
                    continue
                end

                X_donor = m_fe[i, donor_j] / n_donor
                X_mob = max(X_donor - cfg_core.phi_residual, 0.0) * Xfem_cell[i, donor_j]
                m_avail = X_mob * n_donor

                fx = abs(uf) * (dt_sub / dx_val) * m_avail
                fx_req = uf > 0.0 ? fx : -fx
                req_flux_x[i, j] = fx_req

                if fx_req > 0.0
                    outflow_tot[i, j] += fx_req
                    inflow_tot[i, j + 1] += fx_req
                else
                    outflow_tot[i, j + 1] += -fx_req
                    inflow_tot[i, j] += -fx_req
                end
            end
        end

        # 2. Compute unscaled requested fluxes across North-South faces
        @inbounds for i in 1:(Ny_val - 1)
            yf = i * dy_val
            for j in 1:Nx_val
                xf = (j - 0.5) * dx_val
                dxf = xf - xcenter
                dyf = yf - ycenter
                rf = sqrt(dxf^2 + dyf^2)
                if rf > rplanet || rf < 1.0e-3
                    continue
                end

                ny =
                    if gx !== nothing &&
                        gy !== nothing &&
                        i <= size(gy, 1) &&
                        j <= size(gy, 2)
                        gy_f = gy[i, j]
                        gx_f =
                            0.5 *
                            (gx[i, j] + (i + 1 <= size(gx, 1) ? gx[i + 1, j] : gx[i, j]))
                        g_f = sqrt(gx_f^2 + gy_f^2)
                        g_f > 1.0e-10 ? gy_f / g_f : -dyf / rf
                    else
                        -dyf / rf
                    end

                vf = 0.5 * (v_seg_cell[i, j] + v_seg_cell[i + 1, j])
                wf = vf * ny

                donor_i = wf > 0.0 ? i : i + 1
                rec_i = wf > 0.0 ? i + 1 : i

                n_donor = M_rock_markers[donor_i, j]
                n_rec = M_rock_markers[rec_i, j]
                if n_donor == 0 || n_rec == 0
                    continue
                end

                X_donor = m_fe[donor_i, j] / n_donor
                X_mob = max(X_donor - cfg_core.phi_residual, 0.0) * Xfem_cell[donor_i, j]
                m_avail = X_mob * n_donor

                fy = abs(wf) * (dt_sub / dy_val) * m_avail
                fy_req = wf > 0.0 ? fy : -fy
                req_flux_y[i, j] = fy_req

                if fy_req > 0.0
                    outflow_tot[i, j] += fy_req
                    inflow_tot[i + 1, j] += fy_req
                else
                    outflow_tot[i + 1, j] += -fy_req
                    inflow_tot[i, j] += -fy_req
                end
            end
        end

        # 3. Multi-dimensional flux limiters per cell
        @inbounds for j in 1:Nx_val, i in 1:Ny_val
            n_m = M_rock_markers[i, j]
            if n_m > 0
                X_c = m_fe[i, j] / n_m
                m_avail = max(X_c - cfg_core.phi_residual, 0.0) * Xfem_cell[i, j] * n_m
                m_cap = max(cfg_core.phi_pack - X_c, 0.0) * n_m
                alpha_out[i, j] = if outflow_tot[i, j] > m_avail && m_avail > 0.0
                    m_avail / outflow_tot[i, j]
                else
                    (outflow_tot[i, j] > m_avail ? 0.0 : 1.0)
                end
                alpha_in[i, j] = if inflow_tot[i, j] > m_cap && m_cap > 0.0
                    m_cap / inflow_tot[i, j]
                else
                    (inflow_tot[i, j] > m_cap ? 0.0 : 1.0)
                end
            else
                alpha_out[i, j] = 0.0
                alpha_in[i, j] = 0.0
            end
        end

        # 4. Scale fluxes by joint donor-receiver limiters
        @inbounds for j in 1:(Nx_val - 1), i in 1:Ny_val
            fx_req = req_flux_x[i, j]
            if iszero(fx_req)
                flux_x[i, j] = 0.0
            else
                donor_j = fx_req > 0.0 ? j : j + 1
                rec_j = fx_req > 0.0 ? j + 1 : j
                lim = min(alpha_out[i, donor_j], alpha_in[i, rec_j])
                flux_x[i, j] = fx_req * lim
            end
        end

        @inbounds for j in 1:Nx_val, i in 1:(Ny_val - 1)
            fy_req = req_flux_y[i, j]
            if iszero(fy_req)
                flux_y[i, j] = 0.0
            else
                donor_i = fy_req > 0.0 ? i : i + 1
                rec_i = fy_req > 0.0 ? i + 1 : i
                lim = min(alpha_out[donor_i, j], alpha_in[rec_i, j])
                flux_y[i, j] = fy_req * lim
            end
        end

        if track_volatiles
            @inbounds for j in 1:(Nx_val - 1), i in 1:Ny_val
                fx = flux_x[i, j]
                if iszero(fx)
                    flux_H_x[i, j] = 0.0
                    flux_C_x[i, j] = 0.0
                    flux_N_x[i, j] = 0.0
                    flux_S_x[i, j] = 0.0
                else
                    donor_j = fx > 0.0 ? j : j + 1
                    m_d = m_fe[i, donor_j]
                    if m_d > 0.0
                        flux_H_x[i, j] = fx * (m_fe_H[i, donor_j] / m_d)
                        flux_C_x[i, j] = fx * (m_fe_C[i, donor_j] / m_d)
                        flux_N_x[i, j] = fx * (m_fe_N[i, donor_j] / m_d)
                        flux_S_x[i, j] = fx * (m_fe_S[i, donor_j] / m_d)
                    else
                        flux_H_x[i, j] = 0.0
                        flux_C_x[i, j] = 0.0
                        flux_N_x[i, j] = 0.0
                        flux_S_x[i, j] = 0.0
                    end
                end
            end

            @inbounds for j in 1:Nx_val, i in 1:(Ny_val - 1)
                fy = flux_y[i, j]
                if iszero(fy)
                    flux_H_y[i, j] = 0.0
                    flux_C_y[i, j] = 0.0
                    flux_N_y[i, j] = 0.0
                    flux_S_y[i, j] = 0.0
                else
                    donor_i = fy > 0.0 ? i : i + 1
                    m_d = m_fe[donor_i, j]
                    if m_d > 0.0
                        flux_H_y[i, j] = fy * (m_fe_H[donor_i, j] / m_d)
                        flux_C_y[i, j] = fy * (m_fe_C[donor_i, j] / m_d)
                        flux_N_y[i, j] = fy * (m_fe_N[donor_i, j] / m_d)
                        flux_S_y[i, j] = fy * (m_fe_S[donor_i, j] / m_d)
                    else
                        flux_H_y[i, j] = 0.0
                        flux_C_y[i, j] = 0.0
                        flux_N_y[i, j] = 0.0
                        flux_S_y[i, j] = 0.0
                    end
                end
            end
        end

        # 5. Conservative update of cell metal masses
        @inbounds for j in 1:Nx_val, i in 1:Ny_val
            F_w = (j > 1) ? flux_x[i, j - 1] : 0.0
            F_e = (j < Nx_val) ? flux_x[i, j] : 0.0
            F_n = (i > 1) ? flux_y[i - 1, j] : 0.0
            F_s = (i < Ny_val) ? flux_y[i, j] : 0.0
            m_fe[i, j] += (F_w - F_e + F_n - F_s)

            if track_volatiles
                F_H_w = (j > 1) ? flux_H_x[i, j - 1] : 0.0
                F_H_e = (j < Nx_val) ? flux_H_x[i, j] : 0.0
                F_H_n = (i > 1) ? flux_H_y[i - 1, j] : 0.0
                F_H_s = (i < Ny_val) ? flux_H_y[i, j] : 0.0
                m_fe_H[i, j] = max(0.0, m_fe_H[i, j] + (F_H_w - F_H_e + F_H_n - F_H_s))

                F_C_w = (j > 1) ? flux_C_x[i, j - 1] : 0.0
                F_C_e = (j < Nx_val) ? flux_C_x[i, j] : 0.0
                F_C_n = (i > 1) ? flux_C_y[i - 1, j] : 0.0
                F_C_s = (i < Ny_val) ? flux_C_y[i, j] : 0.0
                m_fe_C[i, j] = max(0.0, m_fe_C[i, j] + (F_C_w - F_C_e + F_C_n - F_C_s))

                F_N_w = (j > 1) ? flux_N_x[i, j - 1] : 0.0
                F_N_e = (j < Nx_val) ? flux_N_x[i, j] : 0.0
                F_N_n = (i > 1) ? flux_N_y[i - 1, j] : 0.0
                F_N_s = (i < Ny_val) ? flux_N_y[i, j] : 0.0
                m_fe_N[i, j] = max(0.0, m_fe_N[i, j] + (F_N_w - F_N_e + F_N_n - F_N_s))

                F_S_w = (j > 1) ? flux_S_x[i, j - 1] : 0.0
                F_S_e = (j < Nx_val) ? flux_S_x[i, j] : 0.0
                F_S_n = (i > 1) ? flux_S_y[i - 1, j] : 0.0
                F_S_s = (i < Ny_val) ? flux_S_y[i, j] : 0.0
                m_fe_S[i, j] = max(0.0, m_fe_S[i, j] + (F_S_w - F_S_e + F_S_n - F_S_s))
            end
        end

        # 6. Gravitational potential energy dissipation heating
        @inbounds for j in 1:Nx_val, i in 1:Ny_val
            n_m = M_rock_markers[i, j]
            v_s = v_seg_cell[i, j]
            if n_m > 0 && v_s > 0.0
                phi_m_curr = (m_fe[i, j] / n_m) * Xfem_cell[i, j]
                Q_diss = segregation_dissipation_heating(
                    min(phi_m_curr, 1.0), drho_cell[i, j], g_acc_cell[i, j], v_s
                )
                total_diss_energy += Q_diss * (dx_val * dy_val) * dt_sub
                if Q_seg_grid !== nothing
                    dQ = 0.25 * Q_diss * (dt_sub / dt)
                    if i <= size(Q_seg_grid, 1) && j <= size(Q_seg_grid, 2)
                        Q_seg_grid[i, j] += dQ
                    end
                    if i <= size(Q_seg_grid, 1) && (j + 1) <= size(Q_seg_grid, 2)
                        Q_seg_grid[i, j + 1] += dQ
                    end
                    if (i + 1) <= size(Q_seg_grid, 1) && j <= size(Q_seg_grid, 2)
                        Q_seg_grid[i + 1, j] += dQ
                    end
                    if (i + 1) <= size(Q_seg_grid, 1) && (j + 1) <= size(Q_seg_grid, 2)
                        Q_seg_grid[i + 1, j + 1] += dQ
                    end
                end
            end
        end
    end

    # Distribute net cell mass changes to markers in each cell
    initial_sum = 0.0
    @inbounds for m in 1:marknum
        initial_sum += Xfe_bulk[m]
    end

    @inbounds for m in 1:marknum
        if tm[m] < 3
            rmark = distance(xm[m], ym[m], xcenter, ycenter)
            if rmark <= rplanet
                j_c = clamp(Int(floor(xm[m] / dx_val)) + 1, 1, Nx_val)
                i_c = clamp(Int(floor(ym[m] / dy_val)) + 1, 1, Ny_val)
                n_m = M_rock_markers[i_c, j_c]
                if n_m > 0
                    m_target = m_fe[i_c, j_c]
                    m_init = M_fe_cell[i_c, j_c]
                    dm_cell = m_target - m_init
                    dX = 0.0
                    if dm_cell > 0.0
                        c_tot = cap_cell[i_c, j_c]
                        if c_tot > 0.0
                            frac_gain = min(dm_cell / c_tot, 1.0)
                            dX = frac_gain * max(cfg_core.phi_pack - Xfe_bulk[m], 0.0)
                            Xfe_bulk[m] = clamp(Xfe_bulk[m] + dX, 0.0, cfg_core.phi_pack)
                        else
                            Xfe_bulk[m] = clamp(Xfe_bulk[m], 0.0, cfg_core.phi_pack)
                        end
                    elseif dm_cell < 0.0
                        if m_init > 0.0
                            scale_loss = max(m_target / m_init, 0.0)
                            Xfe_bulk[m] = clamp(
                                Xfe_bulk[m] * scale_loss, 0.0, cfg_core.phi_pack
                            )
                        else
                            Xfe_bulk[m] = 0.0
                        end
                    else
                        Xfe_bulk[m] = clamp(Xfe_bulk[m], 0.0, cfg_core.phi_pack)
                    end

                    if track_volatiles
                        X_new = Xfe_bulk[m]
                        if X_new > 0.0 && m_target > 0.0
                            Xfe_H_m[m] = m_fe_H[i_c, j_c] / m_target
                            Xfe_C_m[m] = m_fe_C[i_c, j_c] / m_target
                            Xfe_N_m[m] = m_fe_N[i_c, j_c] / m_target
                            Xfe_S_m[m] = m_fe_S[i_c, j_c] / m_target
                        else
                            Xfe_H_m[m] = 0.0
                            Xfe_C_m[m] = 0.0
                            Xfe_N_m[m] = 0.0
                            Xfe_S_m[m] = 0.0
                        end
                    end
                end
            end
        end
    end

    # Enforce floating point conservation without creating out-of-bounds markers
    final_sum = 0.0
    @inbounds for m in 1:marknum
        final_sum += Xfe_bulk[m]
    end

    diff_sum = initial_sum - final_sum
    if abs(diff_sum) > 1.0e-12 * initial_sum
        eligible_count = 0
        @inbounds for m in 1:marknum
            if tm[m] < 3 && distance(xm[m], ym[m], xcenter, ycenter) <= rplanet
                if diff_sum > 0.0 && Xfe_bulk[m] < cfg_core.phi_pack
                    eligible_count += 1
                elseif diff_sum < 0.0 && Xfe_bulk[m] > 0.0
                    eligible_count += 1
                end
            end
        end
        if eligible_count > 0
            corr = diff_sum / eligible_count
            @inbounds for m in 1:marknum
                if tm[m] < 3 && distance(xm[m], ym[m], xcenter, ycenter) <= rplanet
                    if diff_sum > 0.0 && Xfe_bulk[m] < cfg_core.phi_pack
                        Xfe_bulk[m] = min(Xfe_bulk[m] + corr, cfg_core.phi_pack)
                    elseif diff_sum < 0.0 && Xfe_bulk[m] > 0.0
                        Xfe_bulk[m] = max(Xfe_bulk[m] + corr, 0.0)
                    end
                end
            end
        end
    end

    # Keep molten metal fraction Xfem consistent with newly segregated bulk metal
    @inbounds for m in 1:marknum
        if tm[m] < 3
            F_fe = compute_metal_melt_fraction(
                tkm[m]; T_eutectic=cfg_core.T_eutectic, dT_metal=cfg_core.dT_metal
            )
            Xfem[m] = Xfe_bulk[m] * F_fe
        else
            Xfem[m] = 0.0
        end
    end

    return (;
        max_v_seg=max_v,
        n_subcycles=n_sub,
        dt_sub=dt_sub,
        total_dissipation_energy=total_diss_energy,
    )
end

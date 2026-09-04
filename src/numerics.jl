
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
    # @timeit to "setup_gravitational_lse()" begin
    RP = Vector{Float64}(undef, Ny1*Nx1)
    SP = Vector{Float64}(undef, Ny1*Nx1)
    # end # @timeit to "setup_gravitational_lse()"
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
    # @timeit to "setup_hydromechanical_lse()" begin
    R = Vector{Float64}(undef, Ny1*Nx1*6)
    S = Vector{Float64}(undef, Ny1*Nx1*6)
    # end # @timeit to "setup_hydromechanical_lse()"
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
    # @timeit to "setup_thermal_lse()" begin
    RT = Vector{Float64}(undef, Ny1*Nx1)
    ST = Vector{Float64}(undef, Ny1*Nx1)
    # end # @timeit to "setup_thermal_lse()"
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

    # @timeit to "get_viscosities_stresses_density_gradients!()" begin
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
    # end # @timeit to "get_viscosities_stresses_density_gradients!()"
end # function get_viscosities_stresses_density_gradients!

"""
Assemble the LHS sparse coefficient matrix and fill RHS coefficient vector
of the Poisson equation to be solved for the gravitational potential Φ.

$(SIGNATURES)

# Details

    - RHO: density at P nodes
    - RP: right hand side coefficient vector

# Returns

    - LP: LHS sparse coefficient matrix

"""
function assemble_gravitational_lse!(RHO, RP; coords=nothing)
    Ny1, Nx1 = size(RHO)
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    xp_val = coords === nothing ? xp : coords.xp
    yp_val = coords === nothing ? yp : coords.yp
    xc_val = coords === nothing ? xcenter : coords.xcenter
    yc_val = coords === nothing ? ycenter : coords.ycenter

    # fresh LHS sparse coefficient matrix
    LP = ExtendableSparseMatrix(Nx1*Ny1, Nx1*Ny1)
    # reset RHS coefficient vector
    RP .= zero(0.0)
    # iterate over P nodes
    for j in 1:1:Nx1, i in 1:1:Ny1
        # define global index in algebraic space
        gk = (j-1) * Ny1 + i
        # decide if external / boundary points
        @inbounds if (
            i==1 ||
            i==Ny1 ||
            j==1 ||
            j==Nx1 ||
            distance(xp_val[j], yp_val[i], xc_val, yc_val) > min(xc_val, yc_val)
        )
            # boundary condition: ϕ = 0
            updateindex!(LP, +, 1.0, gk, gk)
        else
            # internal points: 2D Poisson equation: gravitational potential Φ
            updateindex!(LP, +, inv(dx_val^2), gk, gk-Ny1) # Φ₁
            updateindex!(LP, +, inv(dy_val^2), gk, gk-1) # Φ₂
            updateindex!(LP, +, -2.0*(inv(dx_val^2)+inv(dy_val^2)), gk, gk) # Φ₃
            updateindex!(LP, +, inv(dy_val^2), gk, gk+1) # Φ₄
            updateindex!(LP, +, inv(dx_val^2), gk, gk+Ny1) # Φ₅
            @inbounds RP[gk] = 4.0 * 2.0 * inv(3.0) * π * G * RHO[i, j]
        end
    end
    return LP
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

# Returns

    - L: LHS coefficient matrix
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
)
    Ny1, Nx1 = size(ETAP)
    Nx_val = Nx1 - 1
    Ny_val = Ny1 - 1
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy

    # initialize LHS sparse coefficient matrix
    L = ExtendableSparseMatrix(Nx1*Ny1*6, Nx1*Ny1*6)
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
                # Vx equation internal points: x-Stokes
                #
                #                           kvx-6
                #                            Vx₂
                #                             |
                #               kvy-6     ETA(i-1,j)   kvy+6⋅Ny1-6
                #                Vy₁      GGG(i-1,j)     Vy₃
                #                 *       SXY0(i-1,j)     *
                #                           basic₁
                #                            ETA₁                       
                #                            SXY₁
                #               ETAP(i,j)     |      ETAP(i,j+1)
                #               GGGP(i,j)     |      GGGP(i,j+1) 
                #   kvx-6⋅Ny1   SXX0(i,j)    kvx     SXX0(i,j+1)  kvx+6⋅Ny1
                #     Vx₁---------P₁---------Vx₃---------P₂---------Vx₅
                #                kpm          |        kpm+6⋅Ny1
                #               ETAP₁         |        ETAP₂
                #               SXX₁          |        SXX₂
                #                          ETA(i,j) 
                #                kvy       GGG(i,j)     kvy+6⋅Ny1
                #                Vy₂       SXY0(i,j)      Vy₄
                #                 *        basic₂          * 
                #                            ETA₂ 
                #                            SXY₂
                #                             |
                #                           kvx+6
                #                            Vx₄
                #                             *
                #
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
                # Vy equation internal points: y-Stokes
                #
                #                           kvy-6
                #                            Vy₂
                #                             |
                #                          ETAP(i,j)
                #                          GGGP(i,j)
                #             kvx-6⋅Ny1    SXX0(i,j)     kvx
                #                Vx₁          P₁         Vx₃
                #                 *         ETAP₁         *
                #                            SYY₁
                #               ETA(i,j-1)   kpm       ETA(i,j)
                #               GGG(i,j-1)    |        GGG(i,j)
                #   kvy-6⋅Ny1   SXY0(i,j-1)  kvy       SXY0(i,j)  kvy+6⋅Ny1
                #     Vy₁-------basic₁-------Vv₃-------basic₂-------Vy₅
                #               ETA₁          |        ETA₂     
                #               SXY₁          |        SXY₂
                #                            kpm+6
                #                         ETAP(i+1,j)
                #                         GGGP(i+1,j)
                #          kvx-6⋅Ny1+6    SXX0(i+1,j)   kvx+6
                #                Vx₂          P₂        Vx₄
                #                 *         ETAP₂        *  
                #                            SYY₂
                #                             |
                #                           kvy+6
                #                            Vy₄
                #
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
                # P equation internal points: continuity equation: ∂Vx/∂x+∂Vy/∂y=0
                #
                #                 kvy-6
                #                  Vy₁
                #                   |
                #                   |
                #      kvx-6⋅Ny1   kpm       kvx
                #        Vx₁--------P--------Vx₂
                #                   |
                #                   |
                #                  kvy
                #                  Vy₂
                #
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
                # qxDarcy equation internal points: x-Darcy equation:
                # ηfluid/kᵩx⋅qxDarcy + ∂P/∂x = ρfluid⋅gx
                #
                #        P₁--------qxD--------P₂
                #       kpf        kqx     kpf+6⋅Ny1
                #
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
                # qyDarcy equation internal points: y-Darcy equation:
                # ηfluid/kᵩy⋅qyDarcy + ∂P/∂y = ρfluid⋅gy
                #
                #                  P₁
                #                 kpf
                #                  |
                #                 qyD
                #                 kqy
                #                  |
                #                  P₂
                #                kpf+6
                #
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
                # Ptotal/Pfluid equation internal points: continuity equation:
                # ∂qxD/∂x + ∂qyD/∂y - (Ptotal-Pfluid)/ηϕ = 0.0
                #
                #                 qyD₁
                #                kqy-6
                #                  |
                #       qxD₁-------P-------qxD₂
                #    kqx-6⋅Ny1    kpf      kqx
                #                  |
                #                 qyD₂
                #                 kqy
                #
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
    flush!(L) # finalize CSC matrix
    # end # @timeit to "assemble_hydromechanical_lse()"
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
    # @timeit to "process_hydromechanical_solution!()" begin
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
    # end # @timeit to "process_hydromechanical_solution!()"
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
    # @timeit to "recompute_bulk_viscosity!" begin
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
    # end # @timeit to "recompute_bulk_viscosity!"
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
)
    # @timeit to "compute_Aϕ!()" begin
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
            APHI[i, j] = compaction / PHI[i, j]
        end
        return maximum(abs, @view APHI[2:Ny, 2:Nx]) # includes [2, 2] anchor abberation
    end # @inbounds
    # return maximum(abs, APHI[3:Ny-1, 3:Nx-1]) # no abberation
    # end # @timeit to "compute_Aϕ!()"
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
    # @timeit to "compute_fluid_velocities!()" begin
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
    # end # @timeit to "compute_fluid_velocities!()"
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
    # @timeit to "compute_stress_strainrate!()" begin
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
    # end # @timeit to "compute_stress_strainrate!()"
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
    # @timeit to "symmetrize_p_node_observables!()" begin
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
    # end # @timeit to "symmetrize_p_node_observables!()"
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
    # @timeit to "compute_nodal_adjustment!()" begin
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
    # end # @timeit to "compute_nodal_adjustment!()
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
    # @timeit to "positive_max!()" begin
    @inbounds for i in eachindex(A)
        C[i] = max(0, ifelse(A[i] > B[i], A[i], B[i]))
    end
    # end # @timeit to "positive_max!()"
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
    # @timeit to "finalize_plastic_iteration_pass!()" begin
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
    # end # @timeit to "finalize_plastic_iteration_pass!()"
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

# Returns

    - LT: LHS sparse coefficient matrix
"""
function assemble_thermal_lse!(tk1, RHOCP, KX, KY, HR, HA, HS, DHP, RT, dt; coords=nothing)
    Ny1, Nx1 = size(tk1)
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    # fresh LHS coefficient matrix
    LT = ExtendableSparseMatrix(Ny1*Nx1, Ny1*Nx1)
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
            end
        end
    end # @inbounds
    flush!(LT) # finalize CSC matrix
    return LT
end # function assemble_thermal_lse!

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
    tk0, tk1, tk2, DT, DT0, RHOCP, KX, KY, HR, HA, HS, DHP, RT, ST, dt; coords=nothing
)
    # @timeit to "perform_thermal_iterations!" begin
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
            tk1, RHOCP, KX, KY, HR, HA, HS, DHP, RT, dtt; coords=coords
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
            end
        else
            # second+ thermal iteration passes:
            # update dttsum and adjust timestep
            dttsum += dtt
            dtt = min(dtt, dt-dttsum)
        end
        # increase thermal iteration counter
        titer += 1
    end
    # finalize overall temperature change and advance temperature field
    DT .= tk2 .- tk0
    DT0 .= DT
    # end # @timeit to "perform_thermal_iterations!"
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
    # @timeit to "finalize_thermochemical_iteration_pass" begin
    if titer == 1
        if maxDTcurrent > DTmax
            dt *= (DTmax * inv(maxDTcurrent))
            @info "titer 1: reducing dt due to maxDT: dt=$dt s"
        end
    end
    # end # @timeit to "finalize_thermochemical_iteration_pass"
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
function compute_thermochemical_iteration_outcome(DMP, pf, pf0, titer)
    # @timeit to "compute_thermochemical_iteration_outcome" begin
    pferrcur = maximum(abs, pf-pf0)
    DMPmax = maximum(abs, DMP)
    @info "end thermochemical iter $titer" pferrcur DMPmax
    return pferrcur < pferrmax && (titer>2||DMPmax<=0.0)
    # end # @timeit to "compute_thermochemical_iteration_outcome"
end # function compute_thermochemical_iteration_outcome

# Iterative Krylov solvers, preconditioners, and matrix-free operators for Stokes-Darcy systems.

"""
    MatrixFreeStokesDarcyOperator{T<:AbstractFloat}

Matrix-free linear operator for the condensed 4-variable Stokes-Darcy system.
Evaluates matrix-vector products directly from grid property arrays without
allocating or assembling a global sparse matrix.

# Parameters
- `Ny1::Int`: Number of vertical grid nodes.
- `Nx1::Int`: Number of horizontal grid nodes.
- `dx::T`: Horizontal grid spacing [m].
- `dy::T`: Vertical grid spacing [m].
- `Nx_val::Int`: Number of horizontal grid cells.
- `Ny_val::Int`: Number of vertical grid cells.
- `ETA::Matrix{T}`: Shear viscosity at basic nodes [Pa s].
- `ETAP::Matrix{T}`: Normal viscosity at pressure nodes [Pa s].
- `GGG::Matrix{T}`: Shear modulus at basic nodes [Pa].
- `GGGP::Matrix{T}`: Shear modulus at pressure nodes [Pa].
- `RHOX::Matrix{T}`: Density at horizontal velocity nodes [kg/m^3].
- `RHOY::Matrix{T}`: Density at vertical velocity nodes [kg/m^3].
- `RHOFX::Matrix{T}`: Fluid density at horizontal velocity nodes [kg/m^3].
- `RHOFY::Matrix{T}`: Fluid density at vertical velocity nodes [kg/m^3].
- `RX::Matrix{T}`: Darcy hydraulic drag in x-direction [Pa s/m^2].
- `RY::Matrix{T}`: Darcy hydraulic drag in y-direction [Pa s/m^2].
- `ETAPHI::Matrix{T}`: Bulk viscosity at pressure nodes [Pa s].
- `BETAPHI::Matrix{T}`: Porosity-dependent compressibility [1/Pa].
- `PHI::Matrix{T}`: Porosity [-].
- `gx::Matrix{T}`: Horizontal gravity acceleration [m/s^2].
- `gy::Matrix{T}`: Vertical gravity acceleration [m/s^2].
- `dt::T`: Current timestep [s].
- `betasolid::T`: Solid rock compressibility [1/Pa].
- `betafluid::T`: Fluid compressibility [1/Pa].
- `phimin::T`: Minimum porosity cutoff [-].
- `phimax::T`: Maximum porosity cutoff [-].
- `Kcont::T`: Continuity equation scaling parameter [-].
- `bctop::T`: Top boundary condition velocity factor [-].
- `bcbottom::T`: Bottom boundary condition velocity factor [-].
- `bcleft::T`: Left boundary condition velocity factor [-].
- `bcright::T`: Right boundary condition velocity factor [-].
"""
struct MatrixFreeStokesDarcyOperator{T<:AbstractFloat,M<:AbstractMatrix{T}}
    Ny1::Int
    Nx1::Int
    dx::T
    dy::T
    Nx_val::Int
    Ny_val::Int
    ETA::M
    ETAP::M
    GGG::M
    GGGP::M
    RHOX::M
    RHOY::M
    RHOFX::M
    RHOFY::M
    RX::M
    RY::M
    ETAPHI::M
    BETAPHI::M
    PHI::M
    gx::M
    gy::M
    dt::T
    betasolid::T
    betafluid::T
    phimin::T
    phimax::T
    Kcont::T
    bctop::T
    bcbottom::T
    bcleft::T
    bcright::T
    bc_north::Bool
    bc_south::Bool
    bc_west::Bool
    bc_east::Bool
end

function MatrixFreeStokesDarcyOperator{T,M}(
    Ny1::Int,
    Nx1::Int,
    dx::T,
    dy::T,
    Nx_val::Int,
    Ny_val::Int,
    ETA::M,
    ETAP::M,
    GGG::M,
    GGGP::M,
    RHOX::M,
    RHOY::M,
    RHOFX::M,
    RHOFY::M,
    RX::M,
    RY::M,
    ETAPHI::M,
    BETAPHI::M,
    PHI::M,
    gx::M,
    gy::M,
    dt::T,
    betasolid::T,
    betafluid::T,
    phimin::T,
    phimax::T,
    Kcont::T,
    bctop::T,
    bcbottom::T,
    bcleft::T,
    bcright::T;
    bc_north::Bool=true,
    bc_south::Bool=true,
    bc_west::Bool=true,
    bc_east::Bool=true,
) where {T<:AbstractFloat,M<:AbstractMatrix{T}}
    return MatrixFreeStokesDarcyOperator{T,M}(
        Ny1,
        Nx1,
        dx,
        dy,
        Nx_val,
        Ny_val,
        ETA,
        ETAP,
        GGG,
        GGGP,
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
        dt,
        betasolid,
        betafluid,
        phimin,
        phimax,
        Kcont,
        bctop,
        bcbottom,
        bcleft,
        bcright,
        bc_north,
        bc_south,
        bc_west,
        bc_east,
    )
end

function MatrixFreeStokesDarcyOperator{T}(
    Ny1::Int,
    Nx1::Int,
    dx::T,
    dy::T,
    Nx_val::Int,
    Ny_val::Int,
    ETA::M,
    ETAP::M,
    GGG::M,
    GGGP::M,
    RHOX::M,
    RHOY::M,
    RHOFX::M,
    RHOFY::M,
    RX::M,
    RY::M,
    ETAPHI::M,
    BETAPHI::M,
    PHI::M,
    gx::M,
    gy::M,
    dt::T,
    betasolid::T,
    betafluid::T,
    phimin::T,
    phimax::T,
    Kcont::T,
    bctop::T,
    bcbottom::T,
    bcleft::T,
    bcright::T;
    bc_north::Bool=true,
    bc_south::Bool=true,
    bc_west::Bool=true,
    bc_east::Bool=true,
) where {T<:AbstractFloat,M<:AbstractMatrix{T}}
    return MatrixFreeStokesDarcyOperator{T,M}(
        Ny1,
        Nx1,
        dx,
        dy,
        Nx_val,
        Ny_val,
        ETA,
        ETAP,
        GGG,
        GGGP,
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
        dt,
        betasolid,
        betafluid,
        phimin,
        phimax,
        Kcont,
        bctop,
        bcbottom,
        bcleft,
        bcright,
        bc_north,
        bc_south,
        bc_west,
        bc_east,
    )
end

"""
    MatrixFreeStokesDarcyOperator(ETA, ETAP, GGG, GGGP, RHOX, RHOY, RHOFX, RHOFY, RX, RY, ETAPHI, BETAPHI, PHI, gx, gy, dt; coords, kwargs...)

Construct a `MatrixFreeStokesDarcyOperator` from grid arrays and geometry.
"""
function MatrixFreeStokesDarcyOperator(
    ETA::M,
    ETAP::M,
    GGG::M,
    GGGP::M,
    RHOX::M,
    RHOY::M,
    RHOFX::M,
    RHOFY::M,
    RX::M,
    RY::M,
    ETAPHI::M,
    BETAPHI::M,
    PHI::M,
    gx::M,
    gy::M,
    dt::Real;
    coords::GridCoordinates,
    betasolid::Real=Erebus.betasolid,
    betafluid::Real=Erebus.betafluid,
    phimin::Real=Erebus.phimin,
    phimax::Real=Erebus.phimax,
    Kcont::Real=Erebus.Kcont,
    bctop::Real=Erebus.bctop,
    bcbottom::Real=Erebus.bcbottom,
    bcleft::Real=Erebus.bcleft,
    bcright::Real=Erebus.bcright,
    bc_north::Bool=true,
    bc_south::Bool=true,
    bc_west::Bool=true,
    bc_east::Bool=true,
) where {T<:AbstractFloat,M<:AbstractMatrix{T}}
    return MatrixFreeStokesDarcyOperator{T,M}(
        coords.Ny1,
        coords.Nx1,
        T(coords.dx),
        T(coords.dy),
        coords.Nx,
        coords.Ny,
        ETA,
        ETAP,
        GGG,
        GGGP,
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
        T(dt),
        T(betasolid),
        T(betafluid),
        T(phimin),
        T(phimax),
        T(Kcont),
        T(bctop),
        T(bcbottom),
        T(bcleft),
        T(bcright),
        bc_north,
        bc_south,
        bc_west,
        bc_east,
    )
end

Base.size(op::MatrixFreeStokesDarcyOperator) = (op.Ny1 * op.Nx1 * 4, op.Ny1 * op.Nx1 * 4)
Base.size(op::MatrixFreeStokesDarcyOperator, d::Int) = d in (1, 2) ? op.Ny1 * op.Nx1 * 4 : 1
Base.eltype(::MatrixFreeStokesDarcyOperator{T}) where {T} = T

@inline function is_boundary_vx(
    i::Int,
    j::Int,
    Ny_val::Int,
    Nx_val::Int,
    Ny1::Int,
    Nx1::Int;
    bc_north::Bool=true,
    bc_south::Bool=true,
    bc_west::Bool=true,
    bc_east::Bool=true,
)
    return i == 1 || i == Ny1 || j == 1 || j == Nx1 || (bc_east && j == Nx_val)
end

@inline function is_boundary_vy(
    i::Int,
    j::Int,
    Ny_val::Int,
    Nx_val::Int,
    Ny1::Int,
    Nx1::Int;
    bc_north::Bool=true,
    bc_south::Bool=true,
    bc_west::Bool=true,
    bc_east::Bool=true,
)
    return i == 1 || i == Ny1 || j == 1 || j == Nx1 || (bc_south && i == Ny_val)
end

@inline function is_boundary_p(
    i::Int,
    j::Int,
    Ny_val::Int,
    Nx_val::Int,
    Ny1::Int,
    Nx1::Int;
    bc_north::Bool=true,
    bc_south::Bool=true,
    bc_west::Bool=true,
    bc_east::Bool=true,
)
    return i == 1 ||
           i == Ny1 ||
           j == 1 ||
           j == Nx1 ||
           (bc_north && i == 2 && 2 <= j <= Nx_val) ||
           (bc_west && j == 2 && 2 < i < Ny_val) ||
           (bc_south && i == Ny_val && 2 <= j <= Nx_val) ||
           (bc_east && j == Nx_val && 2 < i < Ny_val)
end

"""
    maxwell_effective_viscosity(eta::T, g::T, dt::T) where {T<:Real}

Compute Maxwell viscoelastic effective viscosity avoiding intermediate overflow.
"""
@inline function maxwell_effective_viscosity(eta::T, g::T, dt::T) where {T<:Real}
    g_dt = g * dt
    denom = g_dt + eta
    return iszero(denom) ? zero(T) : eta * (g_dt / denom)
end

"""
    LinearAlgebra.mul!(y, op::MatrixFreeStokesDarcyOperator, x)

Evaluate matrix-free operator application y = A * x for the condensed 4-variable system.
"""
function LinearAlgebra.mul!(
    y::AbstractVector, op::MatrixFreeStokesDarcyOperator, x::AbstractVector
)
    backend = KernelAbstractions.get_backend(y)
    backend_x = KernelAbstractions.get_backend(x)
    backend == backend_x || throw(
        ArgumentError(
            "Backend mismatch in mul!: y is on $(typeof(backend)) while x is on $(typeof(backend_x))",
        ),
    )
    backend_op = KernelAbstractions.get_backend(op.ETA)
    backend == backend_op || throw(
        ArgumentError(
            "Backend mismatch in mul!: operator arrays are on $(typeof(backend_op)) while vectors are on $(typeof(backend))",
        ),
    )
    if !(backend isa KernelAbstractions.CPU)
        return mul_device!(y, op, x; backend=backend)
    end
    Ny1 = op.Ny1
    Nx1 = op.Nx1
    dx_val = op.dx
    dy_val = op.dy
    Nx_val = op.Nx_val
    Ny_val = op.Ny_val
    dt = op.dt

    x_mat = reshape(x, (4, Ny1, Nx1))
    y_mat = reshape(y, (4, Ny1, Nx1))

    Kcont = op.Kcont
    bctop = op.bctop
    bcbottom = op.bcbottom
    bcleft = op.bcleft
    bcright = op.bcright
    betasolid = op.betasolid
    betafluid = op.betafluid
    phimin = op.phimin
    phimax = op.phimax
    ETA = op.ETA
    ETAP = op.ETAP
    GGG = op.GGG
    GGGP = op.GGGP
    RHOX = op.RHOX
    RHOY = op.RHOY
    RHOFX = op.RHOFX
    RHOFY = op.RHOFY
    RX = op.RX
    RY = op.RY
    ETAPHI = op.ETAPHI
    BETAPHI = op.BETAPHI
    PHI = op.PHI
    gx = op.gx
    gy = op.gy

    # Node-disjoint updates allow thread-parallel execution over columns.
    Threads.@threads for j in 1:Nx1
        @inbounds for i in 1:Ny1
            pt = evaluate_stokes_darcy_point(
                i,
                j,
                x_mat,
                Ny1,
                Nx1,
                dx_val,
                dy_val,
                Nx_val,
                Ny_val,
                dt,
                Kcont,
                bctop,
                bcbottom,
                bcleft,
                bcright,
                betasolid,
                betafluid,
                phimin,
                phimax,
                ETA,
                ETAP,
                GGG,
                GGGP,
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
                op.bc_north,
                op.bc_south,
                op.bc_west,
                op.bc_east,
            )
            y_mat[1, i, j] = pt[1]
            y_mat[2, i, j] = pt[2]
            y_mat[3, i, j] = pt[3]
            y_mat[4, i, j] = pt[4]
        end
    end
    return y
end

"""
    LinearAlgebra.mul!(y, op::MatrixFreeStokesDarcyOperator, x, alpha, beta)

Evaluate 5-argument matrix-vector product y = alpha * A * x + beta * y.
"""
function LinearAlgebra.mul!(
    y::AbstractVector,
    op::MatrixFreeStokesDarcyOperator,
    x::AbstractVector,
    alpha::Number,
    beta::Number,
)
    if iszero(beta)
        LinearAlgebra.mul!(y, op, x)
        if !isone(alpha)
            y .*= alpha
        end
    else
        tmp = similar(y)
        LinearAlgebra.mul!(tmp, op, x)
        y .= alpha .* tmp .+ beta .* y
    end
    return y
end

"""
    compute_operator_diagonal(op::MatrixFreeStokesDarcyOperator{T}) where {T}

Compute diagonal entries of the matrix-free Stokes-Darcy operator in O(N) operations.

# Returns
- `d::Vector{T}`: Vector of diagonal entries of length `Ny1 * Nx1 * 4`.
"""
function compute_operator_diagonal(op::MatrixFreeStokesDarcyOperator{T}) where {T}
    backend_op = KernelAbstractions.get_backend(op.ETA)
    if !(backend_op isa KernelAbstractions.CPU)
        return compute_operator_diagonal_device(op; backend=backend_op)
    end
    Ny1 = op.Ny1
    Nx1 = op.Nx1
    dx_val = op.dx
    dy_val = op.dy
    Nx_val = op.Nx_val
    Ny_val = op.Ny_val
    dt = op.dt
    Kcont = op.Kcont
    bcn = op.bc_north
    bcs = op.bc_south
    bcw = op.bc_west
    bce = op.bc_east

    d = zeros(T, Ny1 * Nx1 * 4)
    d_mat = reshape(d, (4, Ny1, Nx1))

    for j in 1:Nx1, i in 1:Ny1
        # Vx diagonal
        if is_boundary_vx(
            i,
            j,
            Ny_val,
            Nx_val,
            Ny1,
            Nx1;
            bc_north=bcn,
            bc_south=bcs,
            bc_west=bcw,
            bc_east=bce,
        )
            d_mat[1, i, j] = one(T)
        else
            ETA1 = maxwell_effective_viscosity(op.ETA[i - 1, j], op.GGG[i - 1, j], dt)
            ETA2 = maxwell_effective_viscosity(op.ETA[i, j], op.GGG[i, j], dt)
            ETAP1 = maxwell_effective_viscosity(op.ETAP[i, j], op.GGGP[i, j], dt)
            ETAP2 = maxwell_effective_viscosity(op.ETAP[i, j + 1], op.GGGP[i, j + 1], dt)
            dRHOdx = 0.5 * (op.RHOX[i, j + 1] - op.RHOX[i, j - 1]) * inv(dx_val)
            d_mat[1, i, j] =
                -(ETAP1 + ETAP2) * inv(dx_val^2) - (ETA1 + ETA2) * inv(dy_val^2) -
                dRHOdx * op.gx[i, j] * dt
        end

        # Vy diagonal
        if is_boundary_vy(
            i,
            j,
            Ny_val,
            Nx_val,
            Ny1,
            Nx1;
            bc_north=bcn,
            bc_south=bcs,
            bc_west=bcw,
            bc_east=bce,
        )
            d_mat[2, i, j] = one(T)
        else
            ETA1 = maxwell_effective_viscosity(op.ETA[i, j - 1], op.GGG[i, j - 1], dt)
            ETA2 = maxwell_effective_viscosity(op.ETA[i, j], op.GGG[i, j], dt)
            ETAP1 = maxwell_effective_viscosity(op.ETAP[i, j], op.GGGP[i, j], dt)
            ETAP2 = maxwell_effective_viscosity(op.ETAP[i + 1, j], op.GGGP[i + 1, j], dt)
            dRHOdy = 0.5 * (op.RHOY[i + 1, j] - op.RHOY[i - 1, j]) * inv(dy_val)
            d_mat[2, i, j] =
                -(ETA1 + ETA2) * inv(dx_val^2) - (ETAP1 + ETAP2) * inv(dy_val^2) -
                dRHOdy * op.gy[i, j] * dt
        end

        # Pt diagonal
        if i == 1 || i == Ny1 || j == 1 || j == Nx1
            d_mat[3, i, j] = one(T)
        elseif (
            (bcn && i == 2 && 2 <= j <= Nx_val) ||
            (bcw && j == 2 && 2 < i < Ny_val) ||
            (bcs && i == Ny_val && 2 <= j <= Nx_val) ||
            (bce && j == Nx_val && 2 < i < Ny_val)
        )
            d_mat[3, i, j] = Kcont
        else
            betadrained = compute_drained_compressibility(
                op.BETAPHI[i, j],
                op.PHI[i, j],
                op.betasolid;
                phimin=op.phimin,
                phimax=op.phimax,
            )
            d_mat[3, i, j] =
                Kcont * (inv(op.ETAPHI[i, j]) / (1.0 - op.PHI[i, j]) + betadrained / dt)
        end

        # Pf diagonal
        if i == 1 || i == Ny1 || j == 1 || j == Nx1
            d_mat[4, i, j] = one(T)
        elseif (
            (bcn && i == 2 && 2 <= j <= Nx_val) ||
            (bcw && j == 2 && 2 < i < Ny_val) ||
            (bcs && i == Ny_val && 2 <= j <= Nx_val) ||
            (bce && j == Nx_val && 2 < i < Ny_val)
        )
            d_mat[4, i, j] = Kcont
        else
            betadrained = compute_drained_compressibility(
                op.BETAPHI[i, j],
                op.PHI[i, j],
                op.betasolid;
                phimin=op.phimin,
                phimax=op.phimax,
            )
            kbw = compute_biot_willis_coefficient(betadrained, op.betasolid)
            ksk = compute_skempton_coefficient(
                betadrained,
                op.PHI[i, j],
                op.betasolid,
                op.betafluid;
                phimin=op.phimin,
                phimax=op.phimax,
            )

            diag_pf =
                Kcont *
                (inv(op.ETAPHI[i, j]) / (1.0 - op.PHI[i, j]) + betadrained * kbw / ksk / dt)
            if 1 < i < Ny1 && 1 < j < Nx_val
                diag_pf += Kcont / (dx_val^2 * op.RX[i, j])
            end
            if 1 < i < Ny1 && 2 < j <= Nx_val
                diag_pf += Kcont / (dx_val^2 * op.RX[i, j - 1])
            end
            if 1 < j < Nx1 && 1 < i < Ny_val
                diag_pf += Kcont / (dy_val^2 * op.RY[i, j])
            end
            if 1 < j < Nx1 && 2 < i <= Ny_val
                diag_pf += Kcont / (dy_val^2 * op.RY[i - 1, j])
            end
            d_mat[4, i, j] = diag_pf
        end
    end

    return d
end

"""
    AbstractStokesDarcyPreconditioner

Supertype for preconditioners in Stokes-Darcy iterative solves.
"""
abstract type AbstractStokesDarcyPreconditioner end

"""
    DiagonalPreconditioner{T<:AbstractFloat} <: AbstractStokesDarcyPreconditioner

Point-Jacobi diagonal preconditioner. Stores inverse diagonal entries.
"""
struct DiagonalPreconditioner{T<:AbstractFloat} <: AbstractStokesDarcyPreconditioner
    inv_diag::Vector{T}
end

"""
    BlockSchurPreconditioner{T<:AbstractFloat} <: AbstractStokesDarcyPreconditioner

Decoupled velocity-pressure block preconditioner. Stores inverse diagonal
factors partitioned by physical field.
"""
struct BlockSchurPreconditioner{T<:AbstractFloat} <: AbstractStokesDarcyPreconditioner
    inv_diag::Vector{T}
    dof_stride::Int
end

# LinearAlgebra interface for DiagonalPreconditioner
function LinearAlgebra.ldiv!(
    y::AbstractVector, P::DiagonalPreconditioner, x::AbstractVector
)
    @inbounds @simd for i in eachindex(y, x, P.inv_diag)
        y[i] = P.inv_diag[i] * x[i]
    end
    return y
end

function LinearAlgebra.ldiv!(P::DiagonalPreconditioner, x::AbstractVector)
    @inbounds @simd for i in eachindex(x, P.inv_diag)
        x[i] *= P.inv_diag[i]
    end
    return x
end

function LinearAlgebra.mul!(y::AbstractVector, P::DiagonalPreconditioner, x::AbstractVector)
    return LinearAlgebra.ldiv!(y, P, x)
end

function LinearAlgebra.mul!(
    y::AbstractVector,
    P::DiagonalPreconditioner,
    x::AbstractVector,
    alpha::Number,
    beta::Number,
)
    if iszero(beta)
        LinearAlgebra.mul!(y, P, x)
        if !isone(alpha)
            y .*= alpha
        end
    else
        @inbounds for i in eachindex(y, x, P.inv_diag)
            y[i] = alpha * P.inv_diag[i] * x[i] + beta * y[i]
        end
    end
    return y
end

# LinearAlgebra interface for BlockSchurPreconditioner
function LinearAlgebra.ldiv!(
    y::AbstractVector, P::BlockSchurPreconditioner, x::AbstractVector
)
    @inbounds @simd for i in eachindex(y, x, P.inv_diag)
        y[i] = P.inv_diag[i] * x[i]
    end
    return y
end

function LinearAlgebra.ldiv!(P::BlockSchurPreconditioner, x::AbstractVector)
    @inbounds @simd for i in eachindex(x, P.inv_diag)
        x[i] *= P.inv_diag[i]
    end
    return x
end

function LinearAlgebra.mul!(
    y::AbstractVector, P::BlockSchurPreconditioner, x::AbstractVector
)
    return LinearAlgebra.ldiv!(y, P, x)
end

function LinearAlgebra.mul!(
    y::AbstractVector,
    P::BlockSchurPreconditioner,
    x::AbstractVector,
    alpha::Number,
    beta::Number,
)
    if iszero(beta)
        LinearAlgebra.mul!(y, P, x)
        if !isone(alpha)
            y .*= alpha
        end
    else
        @inbounds for i in eachindex(y, x, P.inv_diag)
            y[i] = alpha * P.inv_diag[i] * x[i] + beta * y[i]
        end
    end
    return y
end

"""
    build_diagonal_preconditioner(A; tol=1.0e-30)

Build a `DiagonalPreconditioner` from a sparse matrix or matrix-free operator.
"""
function build_diagonal_preconditioner(A::AbstractMatrix{T}; tol::Real=1.0e-30) where {T}
    d = diag(A)
    inv_d = [abs(v) > tol ? inv(v) : one(T) for v in d]
    return DiagonalPreconditioner{T}(inv_d)
end

function build_diagonal_preconditioner(
    op::MatrixFreeStokesDarcyOperator{T}; tol::Real=1.0e-30
) where {T}
    d = compute_operator_diagonal(op)
    inv_d = [abs(v) > tol ? inv(v) : one(T) for v in d]
    return DiagonalPreconditioner{T}(inv_d)
end

"""
    build_block_schur_preconditioner(A; dof_stride=4, tol=1.0e-30)

Build a `BlockSchurPreconditioner` from a sparse matrix or matrix-free operator.
"""
function build_block_schur_preconditioner(
    A::AbstractMatrix{T};
    coords::Union{GridCoordinates,Nothing}=nothing,
    dof_stride::Int=4,
    tol::Real=1.0e-30,
) where {T}
    d = diag(A)
    inv_d = [abs(v) > tol ? inv(v) : one(T) for v in d]
    if (dof_stride == 4 || dof_stride == 6) && coords !== nothing
        Ny1 = coords.Ny1
        Nx1 = coords.Nx1
        Nx_val = coords.Nx
        Ny_val = coords.Ny
        if length(d) == Ny1 * Nx1 * dof_stride
            inv_d_mat = reshape(inv_d, (dof_stride, Ny1, Nx1))
            d_mat = reshape(d, (dof_stride, Ny1, Nx1))
            dx2 = T(coords.dx)^2
            dy2 = T(coords.dy)^2
            Kcont = T(Erebus.Kcont)
            for j in 2:Nx_val, i in 2:Ny_val
                dvx1 =
                    abs(d_mat[1, i, j - 1]) > tol ? inv(abs(d_mat[1, i, j - 1])) : zero(T)
                dvx2 = abs(d_mat[1, i, j]) > tol ? inv(abs(d_mat[1, i, j])) : zero(T)
                dvy1 =
                    abs(d_mat[2, i - 1, j]) > tol ? inv(abs(d_mat[2, i - 1, j])) : zero(T)
                dvy2 = abs(d_mat[2, i, j]) > tol ? inv(abs(d_mat[2, i, j])) : zero(T)
                schur_pt =
                    d_mat[3, i, j] + Kcont * ((dvx1 + dvx2) / dx2 + (dvy1 + dvy2) / dy2)
                inv_d_mat[3, i, j] = abs(schur_pt) > tol ? inv(schur_pt) : one(T)
            end
        end
    end
    return BlockSchurPreconditioner{T}(inv_d, dof_stride)
end

function build_block_schur_preconditioner(
    op::MatrixFreeStokesDarcyOperator{T};
    coords::Union{GridCoordinates,Nothing}=nothing,
    dof_stride::Int=4,
    tol::Real=1.0e-30,
) where {T}
    d = compute_operator_diagonal(op)
    inv_d = [abs(v) > tol ? inv(v) : one(T) for v in d]
    if dof_stride == 4
        inv_d_mat = reshape(inv_d, (4, op.Ny1, op.Nx1))
        d_mat = reshape(d, (4, op.Ny1, op.Nx1))
        dx2 = op.dx^2
        dy2 = op.dy^2
        Kcont = op.Kcont
        for j in 2:op.Nx_val, i in 2:op.Ny_val
            dvx1 = abs(d_mat[1, i, j - 1]) > tol ? inv(abs(d_mat[1, i, j - 1])) : zero(T)
            dvx2 = abs(d_mat[1, i, j]) > tol ? inv(abs(d_mat[1, i, j])) : zero(T)
            dvy1 = abs(d_mat[2, i - 1, j]) > tol ? inv(abs(d_mat[2, i - 1, j])) : zero(T)
            dvy2 = abs(d_mat[2, i, j]) > tol ? inv(abs(d_mat[2, i, j])) : zero(T)
            schur_pt = d_mat[3, i, j] + Kcont * ((dvx1 + dvx2) / dx2 + (dvy1 + dvy2) / dy2)
            inv_d_mat[3, i, j] = abs(schur_pt) > tol ? inv(schur_pt) : one(T)
        end
    end
    return BlockSchurPreconditioner{T}(inv_d, dof_stride)
end

"""
    solve_hydromechanical_iterative!(
        A, b, S;
        method = :gmres,
        rtol = 1.0e-6,
        atol = 1.0e-10,
        maxiter = 200,
        restart = 50,
        preconditioner = :block_schur,
        P = nothing,
    )

Solve the linear system `A * S = b` using a preconditioned Krylov subspace iterative method.

# Parameters
- `A`: Linear operator (`SparseMatrixCSC` or `MatrixFreeStokesDarcyOperator`).
- `b`: Right-hand side vector.
- `S`: Output solution vector (updated in place).
- `method::Symbol`: Krylov method (`:gmres`, `:fgmres`, or `:bicgstab`).
- `rtol::Real`: Relative convergence tolerance.
- `atol::Real`: Absolute convergence tolerance.
- `maxiter::Int`: Maximum iteration count.
- `restart::Int`: Restart subspace size for GMRES and FGMRES.
- `preconditioner::Symbol`: Preconditioner type (`:block_schur`, `:diagonal`, or `:none`).
- `P`: Optional user-supplied preconditioner.

# Returns
- `(S, stats)`: Solution vector and Krylov convergence statistics struct.
"""
function solve_hydromechanical_iterative!(
    A,
    b::AbstractVector{T},
    S::AbstractVector{T};
    coords::Union{GridCoordinates,Nothing}=nothing,
    method::Symbol=:gmres,
    rtol::Real=1.0e-6,
    atol::Real=1.0e-10,
    maxiter::Int=200,
    restart::Int=50,
    preconditioner::Symbol=:block_schur,
    mg_levels::Int=4,
    mg_pre_smooth::Int=2,
    mg_post_smooth::Int=2,
    mg_omega::Real=0.67,
    mg_smoother::Symbol=:damped_jacobi,
    P=nothing,
) where {T<:AbstractFloat}
    # Validate arguments early
    if method ∉ (:gmres, :fgmres, :bicgstab)
        throw(
            ArgumentError(
                "Unknown Krylov method: $(method). Choose :gmres, :fgmres, or :bicgstab."
            ),
        )
    end
    if P === nothing && preconditioner ∉ (:block_schur, :diagonal, :none, :multigrid)
        throw(ArgumentError("Unknown preconditioner: $(preconditioner)"))
    end

    # Construct preconditioner if not supplied
    prec = if P !== nothing
        P
    elseif preconditioner == :multigrid
        if A isa MatrixFreeStokesDarcyOperator
            build_multigrid_preconditioner(
                A;
                levels=mg_levels,
                pre_smooth=mg_pre_smooth,
                post_smooth=mg_post_smooth,
                omega=mg_omega,
                smoother=mg_smoother,
            )
        else
            @warn "preconditioner=:multigrid is only supported for MatrixFreeStokesDarcyOperator; falling back to :block_schur"
            stride = if coords !== nothing && coords.Ny1 * coords.Nx1 > 0
                Int(size(A, 1) ÷ (coords.Ny1 * coords.Nx1))
            else
                4
            end
            build_block_schur_preconditioner(A; coords=coords, dof_stride=stride)
        end
    elseif preconditioner == :block_schur
        stride = if A isa MatrixFreeStokesDarcyOperator
            4
        elseif coords !== nothing && coords.Ny1 * coords.Nx1 > 0
            Int(size(A, 1) ÷ (coords.Ny1 * coords.Nx1))
        else
            4
        end
        build_block_schur_preconditioner(A; coords=coords, dof_stride=stride)
    elseif preconditioner == :diagonal
        build_diagonal_preconditioner(A)
    elseif preconditioner == :none
        nothing
    end

    K = LinearSolve.Krylov

    sol, stats = if method == :gmres
        if prec !== nothing
            K.gmres(
                A,
                b;
                N=prec,
                rtol=T(rtol),
                atol=T(atol),
                itmax=maxiter,
                restart=true,
                memory=restart,
            )
        else
            K.gmres(
                A,
                b;
                rtol=T(rtol),
                atol=T(atol),
                itmax=maxiter,
                restart=true,
                memory=restart,
            )
        end
    elseif method == :fgmres
        if prec !== nothing
            K.fgmres(
                A,
                b;
                N=prec,
                rtol=T(rtol),
                atol=T(atol),
                itmax=maxiter,
                restart=true,
                memory=restart,
            )
        else
            K.fgmres(
                A,
                b;
                rtol=T(rtol),
                atol=T(atol),
                itmax=maxiter,
                restart=true,
                memory=restart,
            )
        end
    elseif method == :bicgstab
        if prec !== nothing
            K.bicgstab(A, b; N=prec, rtol=T(rtol), atol=T(atol), itmax=maxiter)
        else
            K.bicgstab(A, b; rtol=T(rtol), atol=T(atol), itmax=maxiter)
        end
    else
        throw(
            ArgumentError(
                "Unknown Krylov method: $(method). Choose :gmres, :fgmres, or :bicgstab.",
            ),
        )
    end

    if !stats.solved
        @warn "Krylov solver $(method) did not converge within $(stats.niter) iterations: $(stats.status)"
    end

    S .= sol
    return S, stats
end

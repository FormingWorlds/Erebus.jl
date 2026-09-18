# Hardware acceleration and device-agnostic kernels via KernelAbstractions.jl.

using KernelAbstractions
using LinearAlgebra

"""
    to_device(backend::KernelAbstractions.Backend, a::AbstractArray{T}) where {T}

Transfer array `a` to the device managed by `backend`.
Returns `a` directly when `backend` is a CPU backend.
"""
function to_device(backend::KernelAbstractions.CPU, a::AbstractArray)
    return a isa Array ? a : collect(a)
end

function to_device(backend::KernelAbstractions.Backend, a::AbstractArray{T}) where {T}
    d_a = KernelAbstractions.allocate(backend, T, size(a))
    KernelAbstractions.copyto!(backend, d_a, a)
    KernelAbstractions.synchronize(backend)
    return d_a
end

"""
    to_host(a::AbstractArray{T}) where {T}

Transfer device array `a` to standard host memory.
Returns `a` directly when already resident on host.
"""
function to_host(a::AbstractArray{T}) where {T}
    backend = KernelAbstractions.get_backend(a)
    if backend isa KernelAbstractions.CPU
        return a isa Array ? a : collect(a)
    else
        h_a = Array{T}(undef, size(a))
        KernelAbstractions.copyto!(backend, h_a, a)
        KernelAbstractions.synchronize(backend)
        return h_a
    end
end

"""
    to_device(backend::KernelAbstractions.Backend, op::MatrixFreeStokesDarcyOperator{T}) where {T}

Transfer all property arrays of `op` to the device managed by `backend`.
"""
function to_device(
    backend::KernelAbstractions.Backend, op::MatrixFreeStokesDarcyOperator{T}
) where {T}
    if backend isa KernelAbstractions.CPU && op.ETA isa Matrix{T}
        return op
    end
    ETA_dev = to_device(backend, op.ETA)
    ETAP_dev = to_device(backend, op.ETAP)
    GGG_dev = to_device(backend, op.GGG)
    GGGP_dev = to_device(backend, op.GGGP)
    RHOX_dev = to_device(backend, op.RHOX)
    RHOY_dev = to_device(backend, op.RHOY)
    RHOFX_dev = to_device(backend, op.RHOFX)
    RHOFY_dev = to_device(backend, op.RHOFY)
    RX_dev = to_device(backend, op.RX)
    RY_dev = to_device(backend, op.RY)
    ETAPHI_dev = to_device(backend, op.ETAPHI)
    BETAPHI_dev = to_device(backend, op.BETAPHI)
    PHI_dev = to_device(backend, op.PHI)
    gx_dev = to_device(backend, op.gx)
    gy_dev = to_device(backend, op.gy)

    return MatrixFreeStokesDarcyOperator{T,typeof(ETA_dev)}(
        op.Ny1,
        op.Nx1,
        op.dx,
        op.dy,
        op.Nx_val,
        op.Ny_val,
        ETA_dev,
        ETAP_dev,
        GGG_dev,
        GGGP_dev,
        RHOX_dev,
        RHOY_dev,
        RHOFX_dev,
        RHOFY_dev,
        RX_dev,
        RY_dev,
        ETAPHI_dev,
        BETAPHI_dev,
        PHI_dev,
        gx_dev,
        gy_dev,
        op.dt,
        op.betasolid,
        op.betafluid,
        op.phimin,
        op.phimax,
        op.Kcont,
        op.bctop,
        op.bcbottom,
        op.bcleft,
        op.bcright,
        op.bc_north,
        op.bc_south,
        op.bc_west,
        op.bc_east,
    )
end

"""
    to_device(backend::KernelAbstractions.Backend, lvl::MultigridLevel{T}) where {T}

Transfer level operators and pre-allocated buffers to device memory.
"""
function to_device(backend::KernelAbstractions.Backend, lvl::MultigridLevel{T}) where {T}
    if backend isa KernelAbstractions.CPU && lvl.x_buf isa Vector{T}
        return lvl
    end
    op_dev = to_device(backend, lvl.op)
    inv_diag_dev = to_device(backend, lvl.inv_diag)
    x_buf_dev = to_device(backend, lvl.x_buf)
    r_buf_dev = to_device(backend, lvl.r_buf)
    res_buf_dev = to_device(backend, lvl.res_buf)

    return MultigridLevel{T,typeof(x_buf_dev),typeof(op_dev)}(
        lvl.level,
        lvl.Nx,
        lvl.Ny,
        lvl.Nx1,
        lvl.Ny1,
        lvl.dx,
        lvl.dy,
        op_dev,
        inv_diag_dev,
        x_buf_dev,
        r_buf_dev,
        res_buf_dev,
    )
end

"""
    to_device(backend::KernelAbstractions.Backend, hierarchy::StaggeredGridHierarchy{T}) where {T}

Transfer all grid levels in `hierarchy` to device memory.
"""
function to_device(
    backend::KernelAbstractions.Backend, hierarchy::StaggeredGridHierarchy{T}
) where {T}
    levels_dev = [to_device(backend, lvl) for lvl in hierarchy.levels]
    return StaggeredGridHierarchy{T,eltype(levels_dev)}(levels_dev)
end

"""
    to_host(op::MatrixFreeStokesDarcyOperator{T}) where {T}

Transfer all property arrays of `op` back to host memory.
"""
function to_host(op::MatrixFreeStokesDarcyOperator{T}) where {T}
    backend = KernelAbstractions.get_backend(op.ETA)
    if backend isa KernelAbstractions.CPU && op.ETA isa Matrix{T}
        return op
    end
    ETA_h = to_host(op.ETA)
    return MatrixFreeStokesDarcyOperator{T,typeof(ETA_h)}(
        op.Ny1,
        op.Nx1,
        op.dx,
        op.dy,
        op.Nx_val,
        op.Ny_val,
        ETA_h,
        to_host(op.ETAP),
        to_host(op.GGG),
        to_host(op.GGGP),
        to_host(op.RHOX),
        to_host(op.RHOY),
        to_host(op.RHOFX),
        to_host(op.RHOFY),
        to_host(op.RX),
        to_host(op.RY),
        to_host(op.ETAPHI),
        to_host(op.BETAPHI),
        to_host(op.PHI),
        to_host(op.gx),
        to_host(op.gy),
        op.dt,
        op.betasolid,
        op.betafluid,
        op.phimin,
        op.phimax,
        op.Kcont,
        op.bctop,
        op.bcbottom,
        op.bcleft,
        op.bcright,
        op.bc_north,
        op.bc_south,
        op.bc_west,
        op.bc_east,
    )
end

"""
    to_host(lvl::MultigridLevel{T}) where {T}

Transfer level operators and working buffers back to host memory.
"""
function to_host(lvl::MultigridLevel{T}) where {T}
    backend = KernelAbstractions.get_backend(lvl.x_buf)
    if backend isa KernelAbstractions.CPU &&
        lvl.x_buf isa Vector{T} &&
        lvl.op.ETA isa Matrix{T}
        return lvl
    end
    op_host = to_host(lvl.op)
    inv_diag_host = to_host(lvl.inv_diag)
    x_buf_host = to_host(lvl.x_buf)
    r_buf_host = to_host(lvl.r_buf)
    res_buf_host = to_host(lvl.res_buf)

    return MultigridLevel{T,typeof(x_buf_host),typeof(op_host)}(
        lvl.level,
        lvl.Nx,
        lvl.Ny,
        lvl.Nx1,
        lvl.Ny1,
        lvl.dx,
        lvl.dy,
        op_host,
        inv_diag_host,
        x_buf_host,
        r_buf_host,
        res_buf_host,
    )
end

"""
    to_host(hierarchy::StaggeredGridHierarchy{T}) where {T}

Transfer all grid levels in `hierarchy` back to host memory.
"""
function to_host(hierarchy::StaggeredGridHierarchy{T}) where {T}
    levels_host = [to_host(lvl) for lvl in hierarchy.levels]
    return StaggeredGridHierarchy{T,eltype(levels_host)}(levels_host)
end

# -----------------------------------------------------------------------------
# Point Stencil Evaluation
# -----------------------------------------------------------------------------

@inline function evaluate_stokes_darcy_point(
    i::Int,
    j::Int,
    x_mat,
    Ny1::Int,
    Nx1::Int,
    dx_val::T,
    dy_val::T,
    Nx_val::Int,
    Ny_val::Int,
    dt::T,
    Kcont::T,
    bctop::T,
    bcbottom::T,
    bcleft::T,
    bcright::T,
    betasolid::T,
    betafluid::T,
    phimin::T,
    phimax::T,
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
    bc_north::Bool=true,
    bc_south::Bool=true,
    bc_west::Bool=true,
    bc_east::Bool=true,
) where {T<:AbstractFloat}
    # Equation 1: Solid x-velocity (Vx)
    val_vx = zero(T)
    if is_boundary_vx(
        i,
        j,
        Ny_val,
        Nx_val,
        Ny1,
        Nx1;
        bc_north=bc_north,
        bc_south=bc_south,
        bc_west=bc_west,
        bc_east=bc_east,
    )
        val_vx = x_mat[1, i, j]
        if bc_north && i == 1 && 1 < j < Nx_val
            val_vx += bctop * x_mat[1, i + 1, j]
        end
        if bc_south && i == Ny1 && 1 < j < Nx_val
            val_vx += bcbottom * x_mat[1, i - 1, j]
        end
    else
        ETA1 = maxwell_effective_viscosity(ETA[i - 1, j], GGG[i - 1, j], dt)
        ETA2 = maxwell_effective_viscosity(ETA[i, j], GGG[i, j], dt)
        ETAP1 = maxwell_effective_viscosity(ETAP[i, j], GGGP[i, j], dt)
        ETAP2 = maxwell_effective_viscosity(ETAP[i, j + 1], GGGP[i, j + 1], dt)
        dRHOdx = T(0.5) * (RHOX[i, j + 1] - RHOX[i, j - 1]) * inv(dx_val)
        dRHOdy = T(0.5) * (RHOX[i + 1, j] - RHOX[i - 1, j]) * inv(dy_val)

        val_vx =
            (ETAP1 / dx_val^2) * x_mat[1, i, j - 1] +
            (ETA1 / dy_val^2) * x_mat[1, i - 1, j] +
            (
                -(ETAP1 + ETAP2) * inv(dx_val^2) - (ETA1 + ETA2) * inv(dy_val^2) -
                dRHOdx * gx[i, j] * dt
            ) * x_mat[1, i, j] +
            (ETA2 / dy_val^2) * x_mat[1, i + 1, j] +
            (ETAP2 / dx_val^2) * x_mat[1, i, j + 1] +
            (
                ETAP1 * inv(dx_val) * inv(dy_val) - ETA2 * inv(dx_val) * inv(dy_val) -
                dRHOdy * gx[i, j] * dt * T(0.25)
            ) * x_mat[2, i, j] +
            (
                -ETAP2 * inv(dx_val) * inv(dy_val) + ETA2 * inv(dx_val) * inv(dy_val) -
                dRHOdy * gx[i, j] * dt * T(0.25)
            ) * x_mat[2, i, j + 1] +
            (
                -ETAP1 * inv(dx_val) * inv(dy_val) + ETA1 * inv(dx_val) * inv(dy_val) -
                dRHOdy * gx[i, j] * dt * T(0.25)
            ) * x_mat[2, i - 1, j] +
            (
                ETAP2 * inv(dx_val) * inv(dy_val) - ETA1 * inv(dx_val) * inv(dy_val) -
                dRHOdy * gx[i, j] * dt * T(0.25)
            ) * x_mat[2, i - 1, j + 1] +
            (Kcont * inv(dx_val)) * x_mat[3, i, j] -
            (Kcont * inv(dx_val)) * x_mat[3, i, j + 1]
    end

    # Equation 2: Solid y-velocity (Vy)
    val_vy = zero(T)
    if is_boundary_vy(
        i,
        j,
        Ny_val,
        Nx_val,
        Ny1,
        Nx1;
        bc_north=bc_north,
        bc_south=bc_south,
        bc_west=bc_west,
        bc_east=bc_east,
    )
        val_vy = x_mat[2, i, j]
        if bc_west && j == 1 && 1 < i < Ny_val
            val_vy += bcleft * x_mat[2, i, j + 1]
        end
        if bc_east && j == Nx1 && 1 < i < Ny_val
            val_vy += bcright * x_mat[2, i, j - 1]
        end
    else
        ETA1 = maxwell_effective_viscosity(ETA[i, j - 1], GGG[i, j - 1], dt)
        ETA2 = maxwell_effective_viscosity(ETA[i, j], GGG[i, j], dt)
        ETAP1 = maxwell_effective_viscosity(ETAP[i, j], GGGP[i, j], dt)
        ETAP2 = maxwell_effective_viscosity(ETAP[i + 1, j], GGGP[i + 1, j], dt)
        dRHOdx = T(0.5) * (RHOY[i, j + 1] - RHOY[i, j - 1]) * inv(dx_val)
        dRHOdy = T(0.5) * (RHOY[i + 1, j] - RHOY[i - 1, j]) * inv(dy_val)

        val_vy =
            (ETA1 / dx_val^2) * x_mat[2, i, j - 1] +
            (ETAP1 / dy_val^2) * x_mat[2, i - 1, j] +
            (
                -(ETA1 + ETA2) * inv(dx_val^2) - (ETAP1 + ETAP2) * inv(dy_val^2) -
                dRHOdy * gy[i, j] * dt
            ) * x_mat[2, i, j] +
            (ETAP2 / dy_val^2) * x_mat[2, i + 1, j] +
            (ETA2 / dx_val^2) * x_mat[2, i, j + 1] +
            (
                ETAP1 * inv(dx_val) * inv(dy_val) - ETA2 * inv(dx_val) * inv(dy_val) -
                dRHOdx * gy[i, j] * dt * T(0.25)
            ) * x_mat[1, i, j] +
            (
                -ETAP2 * inv(dx_val) * inv(dy_val) + ETA2 * inv(dx_val) * inv(dy_val) -
                dRHOdx * gy[i, j] * dt * T(0.25)
            ) * x_mat[1, i + 1, j] +
            (
                -ETAP1 * inv(dx_val) * inv(dy_val) + ETA1 * inv(dx_val) * inv(dy_val) -
                dRHOdx * gy[i, j] * dt * T(0.25)
            ) * x_mat[1, i, j - 1] +
            (
                ETAP2 * inv(dx_val) * inv(dy_val) - ETA1 * inv(dx_val) * inv(dy_val) -
                dRHOdx * gy[i, j] * dt * T(0.25)
            ) * x_mat[1, i + 1, j - 1] +
            (Kcont * inv(dy_val)) * x_mat[3, i, j] -
            (Kcont * inv(dy_val)) * x_mat[3, i + 1, j]
    end

    # Equation 3: Total pressure (Pt)
    val_pt = zero(T)
    if i == 1 || i == Ny1 || j == 1 || j == Nx1
        val_pt = x_mat[3, i, j]
    elseif (
        (bc_north && i == 2 && 2 <= j <= Nx_val) ||
        (bc_west && j == 2 && 2 < i < Ny_val) ||
        (bc_south && i == Ny_val && 2 <= j <= Nx_val) ||
        (bc_east && j == Nx_val && 2 < i < Ny_val)
    )
        val_pt = Kcont * x_mat[3, i, j]
    else
        betadrained = compute_drained_compressibility(
            BETAPHI[i, j], PHI[i, j], betasolid; phimin=phimin, phimax=phimax
        )
        kbw = compute_biot_willis_coefficient(betadrained, betasolid)

        val_pt =
            (-one(T) / dx_val) * x_mat[1, i, j - 1] +
            (one(T) / dx_val) * x_mat[1, i, j] +
            (-one(T) / dy_val) * x_mat[2, i - 1, j] +
            (one(T) / dy_val) * x_mat[2, i, j] +
            (Kcont * (inv(ETAPHI[i, j]) / (one(T) - PHI[i, j]) + betadrained / dt)) *
            x_mat[3, i, j] +
            (-Kcont * (inv(ETAPHI[i, j]) / (one(T) - PHI[i, j]) + betadrained * kbw / dt)) *
            x_mat[4, i, j]
    end

    # Equation 4: Fluid pressure (Pf) with condensed Darcy divergence
    val_pf = zero(T)
    if i == 1 || i == Ny1 || j == 1 || j == Nx1
        val_pf = x_mat[4, i, j]
    elseif (
        (bc_north && i == 2 && 2 <= j <= Nx_val) ||
        (bc_west && j == 2 && 2 < i < Ny_val) ||
        (bc_south && i == Ny_val && 2 <= j <= Nx_val) ||
        (bc_east && j == Nx_val && 2 < i < Ny_val)
    )
        val_pf = Kcont * x_mat[4, i, j]
    else
        betadrained = compute_drained_compressibility(
            BETAPHI[i, j], PHI[i, j], betasolid; phimin=phimin, phimax=phimax
        )
        kbw = compute_biot_willis_coefficient(betadrained, betasolid)
        ksk = compute_skempton_coefficient(
            betadrained, PHI[i, j], betasolid, betafluid; phimin=phimin, phimax=phimax
        )

        val_pf =
            (-Kcont * (inv(ETAPHI[i, j]) / (one(T) - PHI[i, j]) + betadrained * kbw / dt)) *
            x_mat[3, i, j]
        diag_pf =
            Kcont *
            (inv(ETAPHI[i, j]) / (one(T) - PHI[i, j]) + betadrained * kbw / ksk / dt)

        if 1 < i < Ny1 && 1 < j < Nx_val
            rx2 = RX[i, j]
            coeff_x2 = Kcont / (dx_val^2 * rx2)
            diag_pf += coeff_x2
            val_pf -= coeff_x2 * x_mat[4, i, j + 1]
        end
        if 1 < i < Ny1 && 2 < j <= Nx_val
            rx1 = RX[i, j - 1]
            coeff_x1 = Kcont / (dx_val^2 * rx1)
            diag_pf += coeff_x1
            val_pf -= coeff_x1 * x_mat[4, i, j - 1]
        end
        if 1 < j < Nx1 && 1 < i < Ny_val
            ry2 = RY[i, j]
            coeff_y2 = Kcont / (dy_val^2 * ry2)
            diag_pf += coeff_y2
            val_pf -= coeff_y2 * x_mat[4, i + 1, j]
        end
        if 1 < j < Nx1 && 2 < i <= Ny_val
            ry1 = RY[i - 1, j]
            coeff_y1 = Kcont / (dy_val^2 * ry1)
            diag_pf += coeff_y1
            val_pf -= coeff_y1 * x_mat[4, i - 1, j]
        end

        val_pf += diag_pf * x_mat[4, i, j]
    end

    return val_vx, val_vy, val_pt, val_pf
end

# -----------------------------------------------------------------------------
# KernelAbstractions Kernels
# -----------------------------------------------------------------------------

@kernel function stokes_darcy_matvec_kernel!(
    y_mat,
    @Const(x_mat),
    Ny1::Int,
    Nx1::Int,
    dx_val,
    dy_val,
    Nx_val::Int,
    Ny_val::Int,
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
    @Const(ETA),
    @Const(ETAP),
    @Const(GGG),
    @Const(GGGP),
    @Const(RHOX),
    @Const(RHOY),
    @Const(RHOFX),
    @Const(RHOFY),
    @Const(RX),
    @Const(RY),
    @Const(ETAPHI),
    @Const(BETAPHI),
    @Const(PHI),
    @Const(gx),
    @Const(gy),
)
    i, j = @index(Global, NTuple)
    if i <= Ny1 && j <= Nx1
        vx, vy, pt, pf = evaluate_stokes_darcy_point(
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
        )
        @inbounds begin
            y_mat[1, i, j] = vx
            y_mat[2, i, j] = vy
            y_mat[3, i, j] = pt
            y_mat[4, i, j] = pf
        end
    end
end

@kernel function stokes_darcy_diag_kernel!(
    d_mat,
    Ny1::Int,
    Nx1::Int,
    dx_val,
    dy_val,
    Nx_val::Int,
    Ny_val::Int,
    dt,
    Kcont,
    betasolid,
    betafluid,
    phimin,
    phimax,
    @Const(ETA),
    @Const(ETAP),
    @Const(GGG),
    @Const(GGGP),
    @Const(RHOX),
    @Const(RHOY),
    @Const(RX),
    @Const(RY),
    @Const(ETAPHI),
    @Const(BETAPHI),
    @Const(PHI),
    @Const(gx),
    @Const(gy),
)
    i, j = @index(Global, NTuple)
    if i <= Ny1 && j <= Nx1
        T = eltype(d_mat)

        # Vx diagonal
        d_vx = one(T)
        if !is_boundary_vx(i, j, Ny_val, Nx_val, Ny1, Nx1)
            ETA1 = maxwell_effective_viscosity(ETA[i - 1, j], GGG[i - 1, j], dt)
            ETA2 = maxwell_effective_viscosity(ETA[i, j], GGG[i, j], dt)
            ETAP1 = maxwell_effective_viscosity(ETAP[i, j], GGGP[i, j], dt)
            ETAP2 = maxwell_effective_viscosity(ETAP[i, j + 1], GGGP[i, j + 1], dt)
            dRHOdx = T(0.5) * (RHOX[i, j + 1] - RHOX[i, j - 1]) * inv(dx_val)
            d_vx =
                -(ETAP1 + ETAP2) * inv(dx_val^2) - (ETA1 + ETA2) * inv(dy_val^2) -
                dRHOdx * gx[i, j] * dt
        end

        # Vy diagonal
        d_vy = one(T)
        if !is_boundary_vy(i, j, Ny_val, Nx_val, Ny1, Nx1)
            ETA1 = maxwell_effective_viscosity(ETA[i, j - 1], GGG[i, j - 1], dt)
            ETA2 = maxwell_effective_viscosity(ETA[i, j], GGG[i, j], dt)
            ETAP1 = maxwell_effective_viscosity(ETAP[i, j], GGGP[i, j], dt)
            ETAP2 = maxwell_effective_viscosity(ETAP[i + 1, j], GGGP[i + 1, j], dt)
            dRHOdy = T(0.5) * (RHOY[i + 1, j] - RHOY[i - 1, j]) * inv(dy_val)
            d_vy =
                -(ETA1 + ETA2) * inv(dx_val^2) - (ETAP1 + ETAP2) * inv(dy_val^2) -
                dRHOdy * gy[i, j] * dt
        end

        # Pt diagonal
        d_pt = one(T)
        if i == 1 || i == Ny1 || j == 1 || j == Nx1
            d_pt = one(T)
        elseif (
            (i == 2 && 2 <= j <= Nx_val) ||
            (j == 2 && 2 < i < Ny_val) ||
            (i == Ny_val && 2 <= j <= Nx_val) ||
            (j == Nx_val && 2 < i < Ny_val)
        )
            d_pt = Kcont
        else
            betadrained = compute_drained_compressibility(
                BETAPHI[i, j], PHI[i, j], betasolid; phimin=phimin, phimax=phimax
            )
            d_pt = Kcont * (inv(ETAPHI[i, j]) / (one(T) - PHI[i, j]) + betadrained / dt)
        end

        # Pf diagonal
        d_pf = one(T)
        if i == 1 || i == Ny1 || j == 1 || j == Nx1
            d_pf = one(T)
        elseif (
            (i == 2 && 2 <= j <= Nx_val) ||
            (j == 2 && 2 < i < Ny_val) ||
            (i == Ny_val && 2 <= j <= Nx_val) ||
            (j == Nx_val && 2 < i < Ny_val)
        )
            d_pf = Kcont
        else
            betadrained = compute_drained_compressibility(
                BETAPHI[i, j], PHI[i, j], betasolid; phimin=phimin, phimax=phimax
            )
            kbw = compute_biot_willis_coefficient(betadrained, betasolid)
            ksk = compute_skempton_coefficient(
                betadrained, PHI[i, j], betasolid, betafluid; phimin=phimin, phimax=phimax
            )
            diag_pf =
                Kcont *
                (inv(ETAPHI[i, j]) / (one(T) - PHI[i, j]) + betadrained * kbw / ksk / dt)
            if 1 < i < Ny1 && 1 < j < Nx_val
                diag_pf += Kcont / (dx_val^2 * RX[i, j])
            end
            if 1 < i < Ny1 && 2 < j <= Nx_val
                diag_pf += Kcont / (dx_val^2 * RX[i, j - 1])
            end
            if 1 < j < Nx1 && 1 < i < Ny_val
                diag_pf += Kcont / (dy_val^2 * RY[i, j])
            end
            if 1 < j < Nx1 && 2 < i <= Ny_val
                diag_pf += Kcont / (dy_val^2 * RY[i - 1, j])
            end
            d_pf = diag_pf
        end

        @inbounds begin
            d_mat[1, i, j] = d_vx
            d_mat[2, i, j] = d_vy
            d_mat[3, i, j] = d_pt
            d_mat[4, i, j] = d_pf
        end
    end
end

@kernel function restrict_4var_kernel!(
    rc_mat,
    @Const(rf_mat),
    Ny_c::Int,
    Nx_c::Int,
    Ny_f::Int,
    Nx_f::Int,
    Ny1_c::Int,
    Nx1_c::Int,
    Ny1_f::Int,
    Nx1_f::Int,
)
    ci, cj = @index(Global, NTuple)
    if ci <= Ny1_c && cj <= Nx1_c
        T = eltype(rc_mat)

        # 1. Gather Vx at coarse node (ci, cj)
        val_vx = zero(T)
        for fj in (2 * cj - 2):(2 * cj)
            if 2 <= fj <= Nx_f
                wj = if isodd(fj)
                    ((fj + 1) ÷ 2 == cj ? one(T) : zero(T))
                else
                    (fj ÷ 2 == cj || (fj ÷ 2 + 1) == cj ? T(0.5) : zero(T))
                end
                if wj > zero(T)
                    for fi in (2 * ci - 2):(2 * ci + 1)
                        if 2 <= fi <= Ny_f
                            if !is_boundary_vx(fi, fj, Ny_f, Nx_f, Ny1_f, Nx1_f)
                                i_mid = isodd(fi) ? (fi + 1) ÷ 2 : fi ÷ 2
                                i_adj = isodd(fi) ? i_mid - 1 : i_mid + 1
                                wi = if i_mid == ci
                                    T(0.75)
                                else
                                    (i_adj == ci ? T(0.25) : zero(T))
                                end
                                if wi > zero(T)
                                    val_vx += T(0.25) * wi * wj * rf_mat[1, fi, fj]
                                end
                            end
                        end
                    end
                end
            end
        end

        # 2. Gather Vy at coarse node (ci, cj)
        val_vy = zero(T)
        for fi in (2 * ci - 2):(2 * ci)
            if 2 <= fi <= Ny_f
                wi = if isodd(fi)
                    ((fi + 1) ÷ 2 == ci ? one(T) : zero(T))
                else
                    (fi ÷ 2 == ci || (fi ÷ 2 + 1) == ci ? T(0.5) : zero(T))
                end
                if wi > zero(T)
                    for fj in (2 * cj - 2):(2 * cj + 1)
                        if 2 <= fj <= Nx_f
                            if !is_boundary_vy(fi, fj, Ny_f, Nx_f, Ny1_f, Nx1_f)
                                j_mid = isodd(fj) ? (fj + 1) ÷ 2 : fj ÷ 2
                                j_adj = isodd(fj) ? j_mid - 1 : j_mid + 1
                                wj = if j_mid == cj
                                    T(0.75)
                                else
                                    (j_adj == cj ? T(0.25) : zero(T))
                                end
                                if wj > zero(T)
                                    val_vy += T(0.25) * wi * wj * rf_mat[2, fi, fj]
                                end
                            end
                        end
                    end
                end
            end
        end

        # 3. Gather Pt and Pf at coarse cell (ci, cj)
        val_pt = zero(T)
        val_pf = zero(T)
        if !is_boundary_p(ci, cj, Ny_c, Nx_c, Ny1_c, Nx1_c)
            for dfj in 0:1, dfi in 0:1
                fi = 2 * ci - 1 + dfi
                fj = 2 * cj - 1 + dfj
                if 2 <= fi <= Ny_f &&
                    2 <= fj <= Nx_f &&
                    !is_boundary_p(fi, fj, Ny_f, Nx_f, Ny1_f, Nx1_f)
                    val_pt += T(0.25) * rf_mat[3, fi, fj]
                    val_pf += T(0.25) * rf_mat[4, fi, fj]
                end
            end
        end

        @inbounds begin
            rc_mat[1, ci, cj] = val_vx
            rc_mat[2, ci, cj] = val_vy
            rc_mat[3, ci, cj] = val_pt
            rc_mat[4, ci, cj] = val_pf
        end
    end
end

@kernel function prolongate_4var_kernel!(
    ef_mat,
    @Const(ec_mat),
    Ny_f::Int,
    Nx_f::Int,
    Ny_c::Int,
    Nx_c::Int,
    Ny1_f::Int,
    Nx1_f::Int,
    Ny1_c::Int,
    Nx1_c::Int,
)
    fi, fj = @index(Global, NTuple)
    if fi <= Ny1_f && fj <= Nx1_f
        T = eltype(ef_mat)

        # 1. Vx Prolongation
        if 2 <= fi <= Ny_f &&
            2 <= fj <= Nx_f &&
            !is_boundary_vx(fi, fj, Ny_f, Nx_f, Ny1_f, Nx1_f)
            i_mid = isodd(fi) ? (fi + 1) ÷ 2 : fi ÷ 2
            i_adj = isodd(fi) ? i_mid - 1 : i_mid + 1
            if isodd(fj)
                j_mid = (fj + 1) ÷ 2
                ef_mat[1, fi, fj] += T(0.75) * ec_mat[1, i_mid, j_mid]
                if 1 <= i_adj <= Ny1_c
                    ef_mat[1, fi, fj] += T(0.25) * ec_mat[1, i_adj, j_mid]
                end
            else
                j_lt = fj ÷ 2
                j_rt = j_lt + 1
                ef_mat[1, fi, fj] += T(0.375) * ec_mat[1, i_mid, j_lt]
                ef_mat[1, fi, fj] += T(0.375) * ec_mat[1, i_mid, j_rt]
                if 1 <= i_adj <= Ny1_c
                    ef_mat[1, fi, fj] += T(0.125) * ec_mat[1, i_adj, j_lt]
                    ef_mat[1, fi, fj] += T(0.125) * ec_mat[1, i_adj, j_rt]
                end
            end
        end

        # 2. Vy Prolongation
        if 2 <= fi <= Ny_f &&
            2 <= fj <= Nx_f &&
            !is_boundary_vy(fi, fj, Ny_f, Nx_f, Ny1_f, Nx1_f)
            j_mid = isodd(fj) ? (fj + 1) ÷ 2 : fj ÷ 2
            j_adj = isodd(fj) ? j_mid - 1 : j_mid + 1
            if isodd(fi)
                i_mid = (fi + 1) ÷ 2
                ef_mat[2, fi, fj] += T(0.75) * ec_mat[2, i_mid, j_mid]
                if 1 <= j_adj <= Nx1_c
                    ef_mat[2, fi, fj] += T(0.25) * ec_mat[2, i_mid, j_adj]
                end
            else
                i_dn = fi ÷ 2
                i_up = i_dn + 1
                ef_mat[2, fi, fj] += T(0.375) * ec_mat[2, i_dn, j_mid]
                ef_mat[2, fi, fj] += T(0.375) * ec_mat[2, i_up, j_mid]
                if 1 <= j_adj <= Nx1_c
                    ef_mat[2, fi, fj] += T(0.125) * ec_mat[2, i_dn, j_adj]
                    ef_mat[2, fi, fj] += T(0.125) * ec_mat[2, i_up, j_adj]
                end
            end
        end

        # 3. Pt and Pf Prolongation
        if 2 <= fi <= Ny_f &&
            2 <= fj <= Nx_f &&
            !is_boundary_p(fi, fj, Ny_f, Nx_f, Ny1_f, Nx1_f)
            ci = (fi + 1) ÷ 2
            cj = (fj + 1) ÷ 2
            if !is_boundary_p(ci, cj, Ny_c, Nx_c, Ny1_c, Nx1_c)
                ef_mat[3, fi, fj] += ec_mat[3, ci, cj]
                ef_mat[4, fi, fj] += ec_mat[4, ci, cj]
            end
        end
    end
end

@kernel function diag_scale_4var_kernel!(
    x_mat, @Const(b_mat), @Const(inv_diag_mat), Ny1::Int, Nx1::Int
)
    i, j = @index(Global, NTuple)
    if i <= Ny1 && j <= Nx1
        @inbounds for comp in 1:4
            x_mat[comp, i, j] = inv_diag_mat[comp, i, j] * b_mat[comp, i, j]
        end
    end
end

@kernel function smooth_jacobi_update_kernel!(
    x_mat, @Const(b_mat), @Const(res_mat), @Const(inv_diag_mat), omega, Ny1::Int, Nx1::Int
)
    i, j = @index(Global, NTuple)
    if i <= Ny1 && j <= Nx1
        @inbounds for v in 1:4
            x_mat[v, i, j] +=
                omega * inv_diag_mat[v, i, j] * (b_mat[v, i, j] - res_mat[v, i, j])
        end
    end
end

@kernel function smooth_jacobi_velocity_kernel!(
    x_mat,
    @Const(b_mat),
    @Const(res_mat),
    @Const(inv_diag_mat),
    omega,
    Ny_val::Int,
    Nx_val::Int,
    Ny1::Int,
    Nx1::Int,
)
    i, j = @index(Global, NTuple)
    if i <= Ny1 && j <= Nx1
        @inbounds begin
            if !is_boundary_vx(i, j, Ny_val, Nx_val, Ny1, Nx1)
                rvx = b_mat[1, i, j] - res_mat[1, i, j]
                x_mat[1, i, j] += omega * inv_diag_mat[1, i, j] * rvx
            end
            if !is_boundary_vy(i, j, Ny_val, Nx_val, Ny1, Nx1)
                rvy = b_mat[2, i, j] - res_mat[2, i, j]
                x_mat[2, i, j] += omega * inv_diag_mat[2, i, j] * rvy
            end
        end
    end
end

@kernel function compute_velocity_res_kernel!(
    r_mat, @Const(b_mat), @Const(res_mat), Ny_val::Int, Nx_val::Int, Ny1::Int, Nx1::Int
)
    i, j = @index(Global, NTuple)
    if i <= Ny1 && j <= Nx1
        @inbounds begin
            T = eltype(r_mat)
            r_mat[1, i, j] = if !is_boundary_vx(i, j, Ny_val, Nx_val, Ny1, Nx1)
                (b_mat[1, i, j] - res_mat[1, i, j])
            else
                zero(T)
            end
            r_mat[2, i, j] = if !is_boundary_vy(i, j, Ny_val, Nx_val, Ny1, Nx1)
                (b_mat[2, i, j] - res_mat[2, i, j])
            else
                zero(T)
            end
            r_mat[3, i, j] = zero(T)
            r_mat[4, i, j] = zero(T)
        end
    end
end

@kernel function smooth_jacobi_darcy_kernel!(
    x_mat,
    @Const(b_mat),
    @Const(res_mat),
    @Const(inv_diag_mat),
    omega,
    Ny_val::Int,
    Nx_val::Int,
    Ny1::Int,
    Nx1::Int,
)
    i, j = @index(Global, NTuple)
    if i <= Ny1 && j <= Nx1
        @inbounds begin
            if !is_boundary_p(i, j, Ny_val, Nx_val, Ny1, Nx1)
                b_val = size(b_mat, 1) == 1 ? b_mat[1, i, j] : b_mat[4, i, j]
                rpf = b_val - res_mat[4, i, j]
                x_mat[4, i, j] += omega * inv_diag_mat[4, i, j] * rpf
            end
        end
    end
end

@kernel function compute_darcy_res_kernel!(
    r_mat, @Const(b_mat), @Const(res_mat), Ny_val::Int, Nx_val::Int, Ny1::Int, Nx1::Int
)
    i, j = @index(Global, NTuple)
    if i <= Ny1 && j <= Nx1
        @inbounds begin
            T = eltype(r_mat)
            r_mat[1, i, j] = zero(T)
            r_mat[2, i, j] = zero(T)
            r_mat[3, i, j] = zero(T)
            r_mat[4, i, j] = if !is_boundary_p(i, j, Ny_val, Nx_val, Ny1, Nx1)
                b_val = size(b_mat, 1) == 1 ? b_mat[1, i, j] : b_mat[4, i, j]
                (b_val - res_mat[4, i, j])
            else
                zero(T)
            end
        end
    end
end

@kernel function smooth_redblack_velocity_kernel!(
    x_mat,
    @Const(b_mat),
    @Const(res_mat),
    @Const(inv_diag_mat),
    omega,
    Ny_val::Int,
    Nx_val::Int,
    Ny1::Int,
    Nx1::Int,
    parity::Int,
)
    i, j = @index(Global, NTuple)
    if i <= Ny1 && j <= Nx1
        if (i + j) % 2 == parity
            @inbounds begin
                if !is_boundary_vx(i, j, Ny_val, Nx_val, Ny1, Nx1)
                    rvx = b_mat[1, i, j] - res_mat[1, i, j]
                    x_mat[1, i, j] += omega * inv_diag_mat[1, i, j] * rvx
                end
                if !is_boundary_vy(i, j, Ny_val, Nx_val, Ny1, Nx1)
                    rvy = b_mat[2, i, j] - res_mat[2, i, j]
                    x_mat[2, i, j] += omega * inv_diag_mat[2, i, j] * rvy
                end
            end
        end
    end
end

@kernel function smooth_redblack_darcy_kernel!(
    x_mat,
    @Const(b_mat),
    @Const(res_mat),
    @Const(inv_diag_mat),
    omega,
    Ny_val::Int,
    Nx_val::Int,
    Ny1::Int,
    Nx1::Int,
    parity::Int,
)
    i, j = @index(Global, NTuple)
    if i <= Ny1 && j <= Nx1
        if (i + j) % 2 == parity
            @inbounds begin
                if !is_boundary_p(i, j, Ny_val, Nx_val, Ny1, Nx1)
                    b_val = size(b_mat, 1) == 1 ? b_mat[1, i, j] : b_mat[4, i, j]
                    rpf = b_val - res_mat[4, i, j]
                    x_mat[4, i, j] += omega * inv_diag_mat[4, i, j] * rpf
                end
            end
        end
    end
end

@kernel function copy_2var_to_4var_kernel!(x4_mat, @Const(x2_mat), Ny1::Int, Nx1::Int)
    i, j = @index(Global, NTuple)
    if i <= Ny1 && j <= Nx1
        @inbounds begin
            x4_mat[1, i, j] = x2_mat[1, i, j]
            x4_mat[2, i, j] = x2_mat[2, i, j]
            x4_mat[3, i, j] = zero(eltype(x4_mat))
            x4_mat[4, i, j] = zero(eltype(x4_mat))
        end
    end
end

@kernel function copy_4var_to_2var_kernel!(x2_mat, @Const(x4_mat), Ny1::Int, Nx1::Int)
    i, j = @index(Global, NTuple)
    if i <= Ny1 && j <= Nx1
        @inbounds begin
            x2_mat[1, i, j] = x4_mat[1, i, j]
            x2_mat[2, i, j] = x4_mat[2, i, j]
        end
    end
end

@kernel function copy_1var_to_4var_kernel!(x4_mat, @Const(x1_mat), Ny1::Int, Nx1::Int)
    i, j = @index(Global, NTuple)
    if i <= Ny1 && j <= Nx1
        @inbounds begin
            x4_mat[1, i, j] = zero(eltype(x4_mat))
            x4_mat[2, i, j] = zero(eltype(x4_mat))
            x4_mat[3, i, j] = zero(eltype(x4_mat))
            x4_mat[4, i, j] = x1_mat[i, j]
        end
    end
end

@kernel function copy_4var_to_1var_kernel!(x1_mat, @Const(x4_mat), Ny1::Int, Nx1::Int)
    i, j = @index(Global, NTuple)
    if i <= Ny1 && j <= Nx1
        @inbounds begin
            x1_mat[i, j] = x4_mat[4, i, j]
        end
    end
end

# -----------------------------------------------------------------------------
# Device API Wrappers
# -----------------------------------------------------------------------------

"""
    mul_device!(y, op, x; backend=KernelAbstractions.get_backend(x), workgroupsize=(16, 16))

Evaluate matrix-free Stokes-Darcy operator application using `KernelAbstractions.jl`.
"""
function mul_device!(
    y::AbstractVector{T},
    op::MatrixFreeStokesDarcyOperator{T},
    x::AbstractVector{T};
    backend::KernelAbstractions.Backend=KernelAbstractions.get_backend(x),
    workgroupsize::Tuple{Int,Int}=(16, 16),
) where {T<:AbstractFloat}
    Ny1 = op.Ny1
    Nx1 = op.Nx1
    x_mat = reshape(x, (4, Ny1, Nx1))
    y_mat = reshape(y, (4, Ny1, Nx1))

    kernel! = stokes_darcy_matvec_kernel!(backend, workgroupsize)
    kernel!(
        y_mat,
        x_mat,
        Ny1,
        Nx1,
        op.dx,
        op.dy,
        op.Nx_val,
        op.Ny_val,
        op.dt,
        op.Kcont,
        op.bctop,
        op.bcbottom,
        op.bcleft,
        op.bcright,
        op.betasolid,
        op.betafluid,
        op.phimin,
        op.phimax,
        op.ETA,
        op.ETAP,
        op.GGG,
        op.GGGP,
        op.RHOX,
        op.RHOY,
        op.RHOFX,
        op.RHOFY,
        op.RX,
        op.RY,
        op.ETAPHI,
        op.BETAPHI,
        op.PHI,
        op.gx,
        op.gy;
        ndrange=(Ny1, Nx1),
    )
    KernelAbstractions.synchronize(backend)
    return y
end

"""
    compute_operator_diagonal_device(op; backend=KernelAbstractions.get_backend(op.ETA), workgroupsize=(16, 16))

Compute diagonal entries of `op` on device without sparse matrix allocation.
"""
function compute_operator_diagonal_device(
    op::MatrixFreeStokesDarcyOperator{T};
    backend::KernelAbstractions.Backend=KernelAbstractions.get_backend(op.ETA),
    workgroupsize::Tuple{Int,Int}=(16, 16),
) where {T<:AbstractFloat}
    Ny1 = op.Ny1
    Nx1 = op.Nx1
    d_mat = KernelAbstractions.allocate(backend, T, (4, Ny1, Nx1))

    kernel! = stokes_darcy_diag_kernel!(backend, workgroupsize)
    kernel!(
        d_mat,
        Ny1,
        Nx1,
        op.dx,
        op.dy,
        op.Nx_val,
        op.Ny_val,
        op.dt,
        op.Kcont,
        op.betasolid,
        op.betafluid,
        op.phimin,
        op.phimax,
        op.ETA,
        op.ETAP,
        op.GGG,
        op.GGGP,
        op.RHOX,
        op.RHOY,
        op.RX,
        op.RY,
        op.ETAPHI,
        op.BETAPHI,
        op.PHI,
        op.gx,
        op.gy;
        ndrange=(Ny1, Nx1),
    )
    KernelAbstractions.synchronize(backend)
    return vec(d_mat)
end

"""
    restrict_4var_device!(r_c, r_f, Ny_c, Nx_c, Ny_f, Nx_f; backend=..., workgroupsize=(16, 16))

Conservative volume-weighted restriction on device.
"""
function restrict_4var_device!(
    r_c::AbstractVector{T},
    r_f::AbstractVector{T},
    Ny_c::Int,
    Nx_c::Int,
    Ny_f::Int,
    Nx_f::Int;
    backend::KernelAbstractions.Backend=KernelAbstractions.get_backend(r_c),
    workgroupsize::Tuple{Int,Int}=(16, 16),
) where {T<:AbstractFloat}
    Ny1_c = Ny_c + 1
    Nx1_c = Nx_c + 1
    Ny1_f = Ny_f + 1
    Nx1_f = Nx_f + 1
    rc_mat = reshape(r_c, (4, Ny1_c, Nx1_c))
    rf_mat = reshape(r_f, (4, Ny1_f, Nx1_f))

    kernel! = restrict_4var_kernel!(backend, workgroupsize)
    kernel!(
        rc_mat,
        rf_mat,
        Ny_c,
        Nx_c,
        Ny_f,
        Nx_f,
        Ny1_c,
        Nx1_c,
        Ny1_f,
        Nx1_f;
        ndrange=(Ny1_c, Nx1_c),
    )
    KernelAbstractions.synchronize(backend)
    return r_c
end

"""
    prolongate_4var_device!(e_f, e_c, Ny_f, Nx_f, Ny_c, Nx_c; backend=..., workgroupsize=(16, 16))

Continuous staggered prolongation on device.
"""
function prolongate_4var_device!(
    e_f::AbstractVector{T},
    e_c::AbstractVector{T},
    Ny_f::Int,
    Nx_f::Int,
    Ny_c::Int,
    Nx_c::Int;
    backend::KernelAbstractions.Backend=KernelAbstractions.get_backend(e_f),
    workgroupsize::Tuple{Int,Int}=(16, 16),
) where {T<:AbstractFloat}
    Ny1_c = Ny_c + 1
    Nx1_c = Nx_c + 1
    Ny1_f = Ny_f + 1
    Nx1_f = Nx_f + 1
    ef_mat = reshape(e_f, (4, Ny1_f, Nx1_f))
    ec_mat = reshape(e_c, (4, Ny1_c, Nx1_c))

    kernel! = prolongate_4var_kernel!(backend, workgroupsize)
    kernel!(
        ef_mat,
        ec_mat,
        Ny_f,
        Nx_f,
        Ny_c,
        Nx_c,
        Ny1_f,
        Nx1_f,
        Ny1_c,
        Nx1_c;
        ndrange=(Ny1_f, Nx1_f),
    )
    KernelAbstractions.synchronize(backend)
    return e_f
end

"""
    smooth_velocity_device!(x, b, op, inv_diag, res_buf; backend=..., omega=0.67, iterations=2, smoother=:damped_jacobi, workgroupsize=(16, 16))

Relax velocity components (Vx, Vy) on device with boundary condition enforcement.
"""
function smooth_velocity_device!(
    x::AbstractVector{T},
    b::AbstractVector{T},
    op::MatrixFreeStokesDarcyOperator{T},
    inv_diag::AbstractVector{T},
    res_buf::AbstractVector{T};
    backend::KernelAbstractions.Backend=KernelAbstractions.get_backend(x),
    omega::Real=0.67,
    iterations::Int=2,
    smoother::Symbol=:damped_jacobi,
    workgroupsize::Tuple{Int,Int}=(16, 16),
) where {T<:AbstractFloat}
    smoother in (:damped_jacobi, :redblack_gauss_seidel) || throw(
        ArgumentError(
            "Unsupported smoother $smoother. Supported: :damped_jacobi, :redblack_gauss_seidel",
        ),
    )
    Ny1 = op.Ny1
    Nx1 = op.Nx1
    x_mat = reshape(x, (4, Ny1, Nx1))
    b_stride = length(b) == 2 * Ny1 * Nx1 ? 2 : 4
    b_mat = reshape(b, (b_stride, Ny1, Nx1))
    res_mat = reshape(res_buf, (4, Ny1, Nx1))
    inv_mat = reshape(inv_diag, (4, Ny1, Nx1))
    omega_T = T(omega)

    if smoother == :damped_jacobi
        kernel! = smooth_jacobi_velocity_kernel!(backend, workgroupsize)
        for _ in 1:iterations
            mul_device!(res_buf, op, x; backend=backend, workgroupsize=workgroupsize)
            kernel!(
                x_mat,
                b_mat,
                res_mat,
                inv_mat,
                omega_T,
                op.Ny_val,
                op.Nx_val,
                Ny1,
                Nx1;
                ndrange=(Ny1, Nx1),
            )
            KernelAbstractions.synchronize(backend)
        end
    elseif smoother == :redblack_gauss_seidel
        kernel! = smooth_redblack_velocity_kernel!(backend, workgroupsize)
        for _ in 1:iterations
            # Red sweep (parity 0)
            mul_device!(res_buf, op, x; backend=backend, workgroupsize=workgroupsize)
            kernel!(
                x_mat,
                b_mat,
                res_mat,
                inv_mat,
                omega_T,
                op.Ny_val,
                op.Nx_val,
                Ny1,
                Nx1,
                0;
                ndrange=(Ny1, Nx1),
            )
            KernelAbstractions.synchronize(backend)

            # Black sweep (parity 1)
            mul_device!(res_buf, op, x; backend=backend, workgroupsize=workgroupsize)
            kernel!(
                x_mat,
                b_mat,
                res_mat,
                inv_mat,
                omega_T,
                op.Ny_val,
                op.Nx_val,
                Ny1,
                Nx1,
                1;
                ndrange=(Ny1, Nx1),
            )
            KernelAbstractions.synchronize(backend)
        end
    end
    return x
end

"""
    smooth_darcy_device!(x, b, op, inv_diag, res_buf; backend=..., omega=0.67, iterations=2, smoother=:damped_jacobi, workgroupsize=(16, 16))

Relax Darcy fluid pressure (Pf) on device with boundary condition enforcement.
"""
function smooth_darcy_device!(
    x::AbstractVector{T},
    b::AbstractVector{T},
    op::MatrixFreeStokesDarcyOperator{T},
    inv_diag::AbstractVector{T},
    res_buf::AbstractVector{T};
    backend::KernelAbstractions.Backend=KernelAbstractions.get_backend(x),
    omega::Real=0.67,
    iterations::Int=2,
    smoother::Symbol=:damped_jacobi,
    workgroupsize::Tuple{Int,Int}=(16, 16),
) where {T<:AbstractFloat}
    smoother in (:damped_jacobi, :redblack_gauss_seidel) || throw(
        ArgumentError(
            "Unsupported smoother $smoother. Supported: :damped_jacobi, :redblack_gauss_seidel",
        ),
    )
    Ny1 = op.Ny1
    Nx1 = op.Nx1
    x_mat = reshape(x, (4, Ny1, Nx1))
    b_stride = length(b) == Ny1 * Nx1 ? 1 : 4
    b_mat = reshape(b, (b_stride, Ny1, Nx1))
    res_mat = reshape(res_buf, (4, Ny1, Nx1))
    inv_mat = reshape(inv_diag, (4, Ny1, Nx1))
    omega_T = T(omega)

    if smoother == :damped_jacobi
        kernel! = smooth_jacobi_darcy_kernel!(backend, workgroupsize)
        for _ in 1:iterations
            mul_device!(res_buf, op, x; backend=backend, workgroupsize=workgroupsize)
            kernel!(
                x_mat,
                b_mat,
                res_mat,
                inv_mat,
                omega_T,
                op.Ny_val,
                op.Nx_val,
                Ny1,
                Nx1;
                ndrange=(Ny1, Nx1),
            )
            KernelAbstractions.synchronize(backend)
        end
    elseif smoother == :redblack_gauss_seidel
        kernel! = smooth_redblack_darcy_kernel!(backend, workgroupsize)
        for _ in 1:iterations
            # Red sweep (parity 0)
            mul_device!(res_buf, op, x; backend=backend, workgroupsize=workgroupsize)
            kernel!(
                x_mat,
                b_mat,
                res_mat,
                inv_mat,
                omega_T,
                op.Ny_val,
                op.Nx_val,
                Ny1,
                Nx1,
                0;
                ndrange=(Ny1, Nx1),
            )
            KernelAbstractions.synchronize(backend)

            # Black sweep (parity 1)
            mul_device!(res_buf, op, x; backend=backend, workgroupsize=workgroupsize)
            kernel!(
                x_mat,
                b_mat,
                res_mat,
                inv_mat,
                omega_T,
                op.Ny_val,
                op.Nx_val,
                Ny1,
                Nx1,
                1;
                ndrange=(Ny1, Nx1),
            )
            KernelAbstractions.synchronize(backend)
        end
    end
    return x
end

"""
    smooth_damped_jacobi_device!(x, b, op, inv_diag, res_buf; backend=..., omega=0.67, iterations=2, workgroupsize=(16, 16))

Device-agnostic relaxation sweeps on velocity unknowns.
"""
function smooth_damped_jacobi_device!(
    x::AbstractVector{T},
    b::AbstractVector{T},
    op::MatrixFreeStokesDarcyOperator{T},
    inv_diag::AbstractVector{T},
    res_buf::AbstractVector{T};
    backend::KernelAbstractions.Backend=KernelAbstractions.get_backend(x),
    omega::Real=0.67,
    iterations::Int=2,
    workgroupsize::Tuple{Int,Int}=(16, 16),
) where {T<:AbstractFloat}
    Ny1 = op.Ny1
    Nx1 = op.Nx1
    x_mat = reshape(x, (4, Ny1, Nx1))
    b_mat = reshape(b, (4, Ny1, Nx1))
    res_mat = reshape(res_buf, (4, Ny1, Nx1))
    inv_diag_mat = reshape(inv_diag, (4, Ny1, Nx1))

    update_k! = smooth_jacobi_update_kernel!(backend, workgroupsize)
    om = T(omega)

    for _ in 1:iterations
        mul_device!(res_buf, op, x; backend=backend, workgroupsize=workgroupsize)
        update_k!(x_mat, b_mat, res_mat, inv_diag_mat, om, Ny1, Nx1; ndrange=(Ny1, Nx1))
        KernelAbstractions.synchronize(backend)
    end
    return x
end

"""
    v_cycle_velocity_device!(xv, bv, hierarchy, level_idx=1; backend=..., smoother=:damped_jacobi, pre_smooth=2, post_smooth=2, omega=0.67, workgroupsize=(16, 16))

Execute recursive geometric multigrid V-cycle for Stokes velocity on device.
"""
function v_cycle_velocity_device!(
    xv::AbstractVector{T},
    bv::AbstractVector{T},
    hierarchy::StaggeredGridHierarchy{T},
    level_idx::Int=1;
    backend::KernelAbstractions.Backend=KernelAbstractions.get_backend(xv),
    smoother::Symbol=:damped_jacobi,
    pre_smooth::Int=2,
    post_smooth::Int=2,
    omega::Real=0.67,
    workgroupsize::Tuple{Int,Int}=(16, 16),
) where {T<:AbstractFloat}
    levels = hierarchy.levels
    curr = levels[level_idx]
    num_levels = length(levels)
    Ny1 = curr.Ny1
    Nx1 = curr.Nx1

    is_2var = (length(xv) == 2 * Ny1 * Nx1)
    if is_2var
        x4 = curr.x_buf
        r4 = curr.r_buf
        copy_k! = copy_2var_to_4var_kernel!(backend, workgroupsize)
        xv_mat = reshape(xv, (2, Ny1, Nx1))
        bv_mat = reshape(bv, (2, Ny1, Nx1))
        x4_mat = reshape(x4, (4, Ny1, Nx1))
        r4_mat = reshape(r4, (4, Ny1, Nx1))
        copy_k!(x4_mat, xv_mat, Ny1, Nx1; ndrange=(Ny1, Nx1))
        copy_k!(r4_mat, bv_mat, Ny1, Nx1; ndrange=(Ny1, Nx1))
        KernelAbstractions.synchronize(backend)
        x_use = x4
        b_use = r4
    else
        x_use = xv
        b_use = bv
    end

    smooth_velocity_device!(
        x_use,
        b_use,
        curr.op,
        curr.inv_diag,
        curr.res_buf;
        backend=backend,
        smoother=smoother,
        omega=omega,
        iterations=pre_smooth,
        workgroupsize=workgroupsize,
    )

    if level_idx == num_levels
        coarse_sweeps = max(pre_smooth + post_smooth, 16)
        smooth_velocity_device!(
            x_use,
            b_use,
            curr.op,
            curr.inv_diag,
            curr.res_buf;
            backend=backend,
            smoother=smoother,
            omega=omega,
            iterations=coarse_sweeps,
            workgroupsize=workgroupsize,
        )
        if is_2var
            copy_back_k! = copy_4var_to_2var_kernel!(backend, workgroupsize)
            xv_mat = reshape(xv, (2, Ny1, Nx1))
            x4_mat = reshape(x_use, (4, Ny1, Nx1))
            copy_back_k!(xv_mat, x4_mat, Ny1, Nx1; ndrange=(Ny1, Nx1))
            KernelAbstractions.synchronize(backend)
        end
        return xv
    end

    mul_device!(curr.res_buf, curr.op, x_use; backend=backend, workgroupsize=workgroupsize)
    r_k! = compute_velocity_res_kernel!(backend, workgroupsize)
    r_mat = reshape(curr.r_buf, (4, Ny1, Nx1))
    b_mat = reshape(b_use, (4, Ny1, Nx1))
    res_mat = reshape(curr.res_buf, (4, Ny1, Nx1))
    r_k!(r_mat, b_mat, res_mat, curr.Ny, curr.Nx, Ny1, Nx1; ndrange=(Ny1, Nx1))
    KernelAbstractions.synchronize(backend)

    next_lvl = levels[level_idx + 1]
    restrict_4var_device!(
        next_lvl.r_buf,
        curr.r_buf,
        next_lvl.Ny,
        next_lvl.Nx,
        curr.Ny,
        curr.Nx;
        backend=backend,
        workgroupsize=workgroupsize,
    )
    fill!(next_lvl.x_buf, zero(T))

    v_cycle_velocity_device!(
        next_lvl.x_buf,
        next_lvl.r_buf,
        hierarchy,
        level_idx + 1;
        backend=backend,
        smoother=smoother,
        pre_smooth=pre_smooth,
        post_smooth=post_smooth,
        omega=omega,
        workgroupsize=workgroupsize,
    )

    prolongate_4var_device!(
        x_use,
        next_lvl.x_buf,
        curr.Ny,
        curr.Nx,
        next_lvl.Ny,
        next_lvl.Nx;
        backend=backend,
        workgroupsize=workgroupsize,
    )

    smooth_velocity_device!(
        x_use,
        b_use,
        curr.op,
        curr.inv_diag,
        curr.res_buf;
        backend=backend,
        smoother=smoother,
        omega=omega,
        iterations=post_smooth,
        workgroupsize=workgroupsize,
    )

    if is_2var
        copy_back_k! = copy_4var_to_2var_kernel!(backend, workgroupsize)
        xv_mat = reshape(xv, (2, Ny1, Nx1))
        x4_mat = reshape(x_use, (4, Ny1, Nx1))
        copy_back_k!(xv_mat, x4_mat, Ny1, Nx1; ndrange=(Ny1, Nx1))
        KernelAbstractions.synchronize(backend)
    end
    return xv
end

"""
    v_cycle_darcy_device!(xpf, bpf, hierarchy, level_idx=1; backend=..., smoother=:damped_jacobi, pre_smooth=2, post_smooth=2, omega=0.67, workgroupsize=(16, 16))

Execute recursive geometric multigrid V-cycle for Darcy fluid pressure on device.
"""
function v_cycle_darcy_device!(
    xpf::AbstractVector{T},
    bpf::AbstractVector{T},
    hierarchy::StaggeredGridHierarchy{T},
    level_idx::Int=1;
    backend::KernelAbstractions.Backend=KernelAbstractions.get_backend(xpf),
    smoother::Symbol=:damped_jacobi,
    pre_smooth::Int=2,
    post_smooth::Int=2,
    omega::Real=0.67,
    workgroupsize::Tuple{Int,Int}=(16, 16),
) where {T<:AbstractFloat}
    levels = hierarchy.levels
    curr = levels[level_idx]
    num_levels = length(levels)
    Ny1 = curr.Ny1
    Nx1 = curr.Nx1

    is_1var = (length(xpf) == Ny1 * Nx1)
    if is_1var
        x4 = curr.x_buf
        r4 = curr.r_buf
        copy_k! = copy_1var_to_4var_kernel!(backend, workgroupsize)
        xpf_mat = reshape(xpf, (Ny1, Nx1))
        bpf_mat = reshape(bpf, (Ny1, Nx1))
        x4_mat = reshape(x4, (4, Ny1, Nx1))
        r4_mat = reshape(r4, (4, Ny1, Nx1))
        copy_k!(x4_mat, xpf_mat, Ny1, Nx1; ndrange=(Ny1, Nx1))
        copy_k!(r4_mat, bpf_mat, Ny1, Nx1; ndrange=(Ny1, Nx1))
        KernelAbstractions.synchronize(backend)
        x_use = x4
        b_use = r4
    else
        x_use = xpf
        b_use = bpf
    end

    smooth_darcy_device!(
        x_use,
        b_use,
        curr.op,
        curr.inv_diag,
        curr.res_buf;
        backend=backend,
        smoother=smoother,
        omega=omega,
        iterations=pre_smooth,
        workgroupsize=workgroupsize,
    )

    if level_idx == num_levels
        coarse_sweeps = max(pre_smooth + post_smooth, 16)
        smooth_darcy_device!(
            x_use,
            b_use,
            curr.op,
            curr.inv_diag,
            curr.res_buf;
            backend=backend,
            smoother=smoother,
            omega=omega,
            iterations=coarse_sweeps,
            workgroupsize=workgroupsize,
        )
        if is_1var
            copy_back_k! = copy_4var_to_1var_kernel!(backend, workgroupsize)
            xpf_mat = reshape(xpf, (Ny1, Nx1))
            x4_mat = reshape(x_use, (4, Ny1, Nx1))
            copy_back_k!(xpf_mat, x4_mat, Ny1, Nx1; ndrange=(Ny1, Nx1))
            KernelAbstractions.synchronize(backend)
        end
        return xpf
    end

    mul_device!(curr.res_buf, curr.op, x_use; backend=backend, workgroupsize=workgroupsize)
    r_k! = compute_darcy_res_kernel!(backend, workgroupsize)
    r_mat = reshape(curr.r_buf, (4, Ny1, Nx1))
    b_mat = reshape(b_use, (4, Ny1, Nx1))
    res_mat = reshape(curr.res_buf, (4, Ny1, Nx1))
    r_k!(r_mat, b_mat, res_mat, curr.Ny, curr.Nx, Ny1, Nx1; ndrange=(Ny1, Nx1))
    KernelAbstractions.synchronize(backend)

    next_lvl = levels[level_idx + 1]
    restrict_4var_device!(
        next_lvl.r_buf,
        curr.r_buf,
        next_lvl.Ny,
        next_lvl.Nx,
        curr.Ny,
        curr.Nx;
        backend=backend,
        workgroupsize=workgroupsize,
    )
    fill!(next_lvl.x_buf, zero(T))

    v_cycle_darcy_device!(
        next_lvl.x_buf,
        next_lvl.r_buf,
        hierarchy,
        level_idx + 1;
        backend=backend,
        smoother=smoother,
        pre_smooth=pre_smooth,
        post_smooth=post_smooth,
        omega=omega,
        workgroupsize=workgroupsize,
    )

    prolongate_4var_device!(
        x_use,
        next_lvl.x_buf,
        curr.Ny,
        curr.Nx,
        next_lvl.Ny,
        next_lvl.Nx;
        backend=backend,
        workgroupsize=workgroupsize,
    )

    smooth_darcy_device!(
        x_use,
        b_use,
        curr.op,
        curr.inv_diag,
        curr.res_buf;
        backend=backend,
        smoother=smoother,
        omega=omega,
        iterations=post_smooth,
        workgroupsize=workgroupsize,
    )

    if is_1var
        copy_back_k! = copy_4var_to_1var_kernel!(backend, workgroupsize)
        xpf_mat = reshape(xpf, (Ny1, Nx1))
        x4_mat = reshape(x_use, (4, Ny1, Nx1))
        copy_back_k!(xpf_mat, x4_mat, Ny1, Nx1; ndrange=(Ny1, Nx1))
        KernelAbstractions.synchronize(backend)
    end
    return xpf
end

"""
    apply_multigrid_vcycle_device!(e, r, hierarchy; backend=..., smoother=:damped_jacobi, pre_smooth=2, post_smooth=2, omega=0.67, workgroupsize=(16, 16))

Execute decoupled geometric multigrid V-cycle on device memory.
"""
function apply_multigrid_vcycle_device!(
    e::AbstractVector{T},
    r::AbstractVector{T},
    hierarchy::StaggeredGridHierarchy{T};
    backend::KernelAbstractions.Backend=KernelAbstractions.get_backend(e),
    smoother::Symbol=:damped_jacobi,
    pre_smooth::Int=2,
    post_smooth::Int=2,
    omega::Real=0.67,
    workgroupsize::Tuple{Int,Int}=(16, 16),
) where {T<:AbstractFloat}
    smoother in (:damped_jacobi, :redblack_gauss_seidel) || throw(
        ArgumentError(
            "Unsupported smoother $smoother. Supported: :damped_jacobi, :redblack_gauss_seidel",
        ),
    )
    curr = hierarchy.levels[1]
    Ny1 = curr.Ny1
    Nx1 = curr.Nx1

    scale_k! = diag_scale_4var_kernel!(backend, workgroupsize)
    e_mat = reshape(e, (4, Ny1, Nx1))
    r_mat = reshape(r, (4, Ny1, Nx1))
    inv_mat = reshape(curr.inv_diag, (4, Ny1, Nx1))
    scale_k!(e_mat, r_mat, inv_mat, Ny1, Nx1; ndrange=(Ny1, Nx1))
    KernelAbstractions.synchronize(backend)

    v_cycle_velocity_device!(
        e,
        r,
        hierarchy,
        1;
        backend=backend,
        smoother=smoother,
        pre_smooth=pre_smooth,
        post_smooth=post_smooth,
        omega=omega,
        workgroupsize=workgroupsize,
    )

    v_cycle_darcy_device!(
        e,
        r,
        hierarchy,
        1;
        backend=backend,
        smoother=smoother,
        pre_smooth=pre_smooth,
        post_smooth=post_smooth,
        omega=omega,
        workgroupsize=workgroupsize,
    )

    return e
end

"""
    apply_multigrid_vcycle!(e, r, hierarchy; kwargs...)

Execute a geometric multigrid V-cycle on `hierarchy`.
"""
function apply_multigrid_vcycle!(
    e::AbstractVector{T},
    r::AbstractVector{T},
    hierarchy::StaggeredGridHierarchy{T};
    kwargs...,
) where {T<:AbstractFloat}
    return apply_multigrid_vcycle_device!(e, r, hierarchy; kwargs...)
end

"""
    to_device(backend::KernelAbstractions.Backend, P::MultigridPreconditioner{T}) where {T}

Transfer a `MultigridPreconditioner` to device memory.
"""
function to_device(
    backend::KernelAbstractions.Backend, P::MultigridPreconditioner{T}
) where {T}
    h_dev = to_device(backend, P.hierarchy)
    return MultigridPreconditioner(h_dev, P.pre_smooth, P.post_smooth, P.omega, P.smoother)
end

"""
    to_host(P::MultigridPreconditioner{T}) where {T}

Transfer a `MultigridPreconditioner` back to host memory.
"""
function to_host(P::MultigridPreconditioner{T}) where {T}
    return MultigridPreconditioner(
        to_host(P.hierarchy), P.pre_smooth, P.post_smooth, P.omega, P.smoother
    )
end

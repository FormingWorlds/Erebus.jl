# Geometric Multigrid (GMG) solvers, transfer operators, and preconditioners for Stokes-Darcy systems.

using LinearAlgebra

"""
    MultigridLevel{T<:AbstractFloat}

Grid level representation within a geometric multigrid hierarchy.

# Parameters
- `level::Int`: Index of this grid level (1 = finest).
- `Nx::Int`: Number of horizontal cells.
- `Ny::Int`: Number of vertical cells.
- `Nx1::Int`: Number of horizontal nodes (`Nx + 1`).
- `Ny1::Int`: Number of vertical nodes (`Ny + 1`).
- `dx::T`: Horizontal grid spacing [m].
- `dy::T`: Vertical grid spacing [m].
- `op::MatrixFreeStokesDarcyOperator{T}`: Matrix-free linear operator at this level.
- `inv_diag::Vector{T}`: Precomputed inverse diagonal elements for smoothers.
- `x_buf::Vector{T}`: Pre-allocated solution and correction vector.
- `r_buf::Vector{T}`: Pre-allocated residual and right-hand side vector.
- `res_buf::Vector{T}`: Pre-allocated operator evaluation buffer.
"""
struct MultigridLevel{
    T<:AbstractFloat,V<:AbstractVector{T},O<:MatrixFreeStokesDarcyOperator{T}
}
    level::Int
    Nx::Int
    Ny::Int
    Nx1::Int
    Ny1::Int
    dx::T
    dy::T
    op::O
    inv_diag::V
    x_buf::V
    r_buf::V
    res_buf::V
end

function MultigridLevel(
    level::Int,
    Nx::Int,
    Ny::Int,
    Nx1::Int,
    Ny1::Int,
    dx::T,
    dy::T,
    op::O,
    inv_diag::V,
    x_buf::V,
    r_buf::V,
) where {T<:AbstractFloat,V<:AbstractVector{T},O<:MatrixFreeStokesDarcyOperator{T}}
    res_buf = similar(x_buf)
    fill!(res_buf, zero(T))
    return MultigridLevel{T,V,O}(
        level, Nx, Ny, Nx1, Ny1, dx, dy, op, inv_diag, x_buf, r_buf, res_buf
    )
end

function MultigridLevel{T}(
    level::Int,
    Nx::Int,
    Ny::Int,
    Nx1::Int,
    Ny1::Int,
    dx::T,
    dy::T,
    op::O,
    inv_diag::V,
    x_buf::V,
    r_buf::V,
    res_buf::V,
) where {T<:AbstractFloat,V<:AbstractVector{T},O<:MatrixFreeStokesDarcyOperator{T}}
    return MultigridLevel{T,V,O}(
        level, Nx, Ny, Nx1, Ny1, dx, dy, op, inv_diag, x_buf, r_buf, res_buf
    )
end

"""
    StaggeredGridHierarchy{T<:AbstractFloat, L<:MultigridLevel{T}}

Multi-level geometric grid hierarchy for staggered-grid Stokes-Darcy systems.

# Parameters
- `levels::Vector{L}`: Ordered list of grid levels from finest to coarsest.
"""
struct StaggeredGridHierarchy{T<:AbstractFloat,L<:MultigridLevel{T}}
    levels::Vector{L}
end

function StaggeredGridHierarchy(levels::Vector{L}) where {T,L<:MultigridLevel{T}}
    return StaggeredGridHierarchy{T,L}(levels)
end

function StaggeredGridHierarchy{T}(levels::Vector{L}) where {T,L<:MultigridLevel{T}}
    return StaggeredGridHierarchy{T,L}(levels)
end

"""
    MultigridPreconditioner{T<:AbstractFloat, H<:StaggeredGridHierarchy{T}} <: AbstractStokesDarcyPreconditioner

Geometric multigrid preconditioner applying decoupled V-cycles on velocity and fluid pressure.

# Parameters
- `hierarchy::H`: Multigrid level hierarchy.
- `pre_smooth::Int`: Number of pre-smoothing relaxation sweeps.
- `post_smooth::Int`: Number of post-smoothing relaxation sweeps.
- `omega::T`: Relaxation damping parameter (e.g. 2/3 for damped Jacobi).
- `smoother::Symbol`: Relaxation method (`:damped_jacobi` or `:redblack_gauss_seidel`).
"""
struct MultigridPreconditioner{T<:AbstractFloat,H<:StaggeredGridHierarchy{T}} <:
       AbstractStokesDarcyPreconditioner
    hierarchy::H
    pre_smooth::Int
    post_smooth::Int
    omega::T
    smoother::Symbol
end

# Backward-compatible constructor accepting optional dof_stride keyword/positional
function MultigridPreconditioner(
    hierarchy::H,
    pre_smooth::Int,
    post_smooth::Int,
    omega::Real,
    smoother::Symbol,
    dof_stride::Int=4,
) where {T<:AbstractFloat,H<:StaggeredGridHierarchy{T}}
    dof_stride == 4 || throw(
        ArgumentError(
            "MultigridPreconditioner requires dof_stride == 4, got $(dof_stride)"
        ),
    )
    return MultigridPreconditioner{T,H}(
        hierarchy, pre_smooth, post_smooth, T(omega), smoother
    )
end

"""
    restrict_4var!(r_c, r_f, Ny_c, Nx_c, Ny_f, Nx_f)

Restrict 4-variable Stokes-Darcy residual from fine grid to coarse grid with boundary guards.
Exact adjoint of prolongation scaled by 1/4 (R = 1/4 P^T).
"""
function restrict_4var!(
    r_c::AbstractVector{T},
    r_f::AbstractVector{T},
    Ny_c::Int,
    Nx_c::Int,
    Ny_f::Int,
    Nx_f::Int,
) where {T<:AbstractFloat}
    Ny1_c = Ny_c + 1
    Nx1_c = Nx_c + 1
    Ny1_f = Ny_f + 1
    Nx1_f = Nx_f + 1
    rc_mat = reshape(r_c, (4, Ny1_c, Nx1_c))
    rf_mat = reshape(r_f, (4, Ny1_f, Nx1_f))

    fill!(r_c, zero(T))

    @inbounds for fj in 2:Nx_f, fi in 2:Ny_f
        # 1. Restrict Vx (Adjoint of Vx prolongation scaled by 1/4)
        if !is_boundary_vx(fi, fj, Ny_f, Nx_f, Ny1_f, Nx1_f)
            rf_val = 0.25 * rf_mat[1, fi, fj]
            i_mid = isodd(fi) ? (fi + 1) ÷ 2 : fi ÷ 2
            i_adj = isodd(fi) ? i_mid - 1 : i_mid + 1
            if isodd(fj)
                j_mid = (fj + 1) ÷ 2
                rc_mat[1, i_mid, j_mid] += 0.75 * rf_val
                if 1 <= i_adj <= Ny1_c
                    rc_mat[1, i_adj, j_mid] += 0.25 * rf_val
                end
            else
                j_lt = fj ÷ 2
                j_rt = j_lt + 1
                rc_mat[1, i_mid, j_lt] += 0.5 * 0.75 * rf_val
                rc_mat[1, i_mid, j_rt] += 0.5 * 0.75 * rf_val
                if 1 <= i_adj <= Ny1_c
                    rc_mat[1, i_adj, j_lt] += 0.5 * 0.25 * rf_val
                    rc_mat[1, i_adj, j_rt] += 0.5 * 0.25 * rf_val
                end
            end
        end

        # 2. Restrict Vy (Adjoint of Vy prolongation scaled by 1/4)
        if !is_boundary_vy(fi, fj, Ny_f, Nx_f, Ny1_f, Nx1_f)
            rf_val = 0.25 * rf_mat[2, fi, fj]
            j_mid = isodd(fj) ? (fj + 1) ÷ 2 : fj ÷ 2
            j_adj = isodd(fj) ? j_mid - 1 : j_mid + 1
            if isodd(fi)
                i_mid = (fi + 1) ÷ 2
                rc_mat[2, i_mid, j_mid] += 0.75 * rf_val
                if 1 <= j_adj <= Nx1_c
                    rc_mat[2, i_mid, j_adj] += 0.25 * rf_val
                end
            else
                i_dn = fi ÷ 2
                i_up = i_dn + 1
                rc_mat[2, i_dn, j_mid] += 0.5 * 0.75 * rf_val
                rc_mat[2, i_up, j_mid] += 0.5 * 0.75 * rf_val
                if 1 <= j_adj <= Nx1_c
                    rc_mat[2, i_dn, j_adj] += 0.5 * 0.25 * rf_val
                    rc_mat[2, i_up, j_adj] += 0.5 * 0.25 * rf_val
                end
            end
        end

        # 3. Restrict Pressures Pt and Pf (Volume-weighted 4-point average)
        if !is_boundary_p(fi, fj, Ny_f, Nx_f, Ny1_f, Nx1_f)
            ci = (fi + 1) ÷ 2
            cj = (fj + 1) ÷ 2
            if !is_boundary_p(ci, cj, Ny_c, Nx_c, Ny1_c, Nx1_c)
                rc_mat[3, ci, cj] += 0.25 * rf_mat[3, fi, fj]
                rc_mat[4, ci, cj] += 0.25 * rf_mat[4, fi, fj]
            end
        end
    end

    return r_c
end

"""
    prolongate_4var!(e_f, e_c, Ny_f, Nx_f, Ny_c, Nx_c)

Prolongate coarse-grid correction `e_c` to fine grid `e_f` in place (`e_f .+= P * e_c`).
"""
function prolongate_4var!(
    e_f::AbstractVector{T},
    e_c::AbstractVector{T},
    Ny_f::Int,
    Nx_f::Int,
    Ny_c::Int,
    Nx_c::Int,
) where {T<:AbstractFloat}
    Ny1_c = Ny_c + 1
    Nx1_c = Nx_c + 1
    Ny1_f = Ny_f + 1
    Nx1_f = Nx_f + 1
    ef_mat = reshape(e_f, (4, Ny1_f, Nx1_f))
    ec_mat = reshape(e_c, (4, Ny1_c, Nx1_c))

    @inbounds for fj in 2:Nx_f, fi in 2:Ny_f
        # 1. Prolongate Vx (normal to vertical faces)
        if !is_boundary_vx(fi, fj, Ny_f, Nx_f, Ny1_f, Nx1_f)
            i_mid = isodd(fi) ? (fi + 1) ÷ 2 : fi ÷ 2
            i_adj = isodd(fi) ? i_mid - 1 : i_mid + 1
            if isodd(fj)
                j_mid = (fj + 1) ÷ 2
                ef_mat[1, fi, fj] += 0.75 * ec_mat[1, i_mid, j_mid]
                if 1 <= i_adj <= Ny1_c
                    ef_mat[1, fi, fj] += 0.25 * ec_mat[1, i_adj, j_mid]
                end
            else
                j_lt = fj ÷ 2
                j_rt = j_lt + 1
                ef_mat[1, fi, fj] += 0.5 * 0.75 * ec_mat[1, i_mid, j_lt]
                ef_mat[1, fi, fj] += 0.5 * 0.75 * ec_mat[1, i_mid, j_rt]
                if 1 <= i_adj <= Ny1_c
                    ef_mat[1, fi, fj] += 0.5 * 0.25 * ec_mat[1, i_adj, j_lt]
                    ef_mat[1, fi, fj] += 0.5 * 0.25 * ec_mat[1, i_adj, j_rt]
                end
            end
        end

        # 2. Prolongate Vy (normal to horizontal faces)
        if !is_boundary_vy(fi, fj, Ny_f, Nx_f, Ny1_f, Nx1_f)
            j_mid = isodd(fj) ? (fj + 1) ÷ 2 : fj ÷ 2
            j_adj = isodd(fj) ? j_mid - 1 : j_mid + 1
            if isodd(fi)
                i_mid = (fi + 1) ÷ 2
                ef_mat[2, fi, fj] += 0.75 * ec_mat[2, i_mid, j_mid]
                if 1 <= j_adj <= Nx1_c
                    ef_mat[2, fi, fj] += 0.25 * ec_mat[2, i_mid, j_adj]
                end
            else
                i_dn = fi ÷ 2
                i_up = i_dn + 1
                ef_mat[2, fi, fj] += 0.5 * 0.75 * ec_mat[2, i_dn, j_mid]
                ef_mat[2, fi, fj] += 0.5 * 0.75 * ec_mat[2, i_up, j_mid]
                if 1 <= j_adj <= Nx1_c
                    ef_mat[2, fi, fj] += 0.5 * 0.25 * ec_mat[2, i_dn, j_adj]
                    ef_mat[2, fi, fj] += 0.5 * 0.25 * ec_mat[2, i_up, j_adj]
                end
            end
        end

        # 3. Prolongate Pressures Pt and Pf (cell-centered piecewise constant injection)
        if !is_boundary_p(fi, fj, Ny_f, Nx_f, Ny1_f, Nx1_f)
            ci = (fi + 1) ÷ 2
            cj = (fj + 1) ÷ 2
            if !is_boundary_p(ci, cj, Ny_c, Nx_c, Ny1_c, Nx1_c)
                ef_mat[3, fi, fj] += ec_mat[3, ci, cj]
                ef_mat[4, fi, fj] += ec_mat[4, ci, cj]
            end
        end
    end

    return e_f
end

"""
    smooth_damped_jacobi!(x, b, op, inv_diag, res_buf; omega=0.67, iterations=2)

Perform damped Jacobi relaxation sweeps: `x <- x + omega * D^{-1} * (b - A * x)`.
"""
function smooth_damped_jacobi!(
    x::AbstractVector{T},
    b::AbstractVector{T},
    op::MatrixFreeStokesDarcyOperator{T},
    inv_diag::AbstractVector{T},
    res_buf::AbstractVector{T};
    omega::Real=0.67,
    iterations::Int=2,
) where {T<:AbstractFloat}
    omega_T = T(omega)
    for _ in 1:iterations
        LinearAlgebra.mul!(res_buf, op, x)
        @inbounds @simd for i in eachindex(res_buf, b, inv_diag, x)
            r = b[i] - res_buf[i]
            x[i] += omega_T * inv_diag[i] * r
        end
    end
    return x
end

"""
    smooth_velocity!(x4, b, op, inv_d, res_buf; omega=0.67, iterations=2, smoother=:damped_jacobi)

Perform relaxation sweeps on the velocity unknowns with boundary condition enforcement.
"""
function smooth_velocity!(
    x4::AbstractVector{T},
    b::AbstractVector{T},
    op::MatrixFreeStokesDarcyOperator{T},
    inv_d::AbstractVector{T},
    res_buf::AbstractVector{T};
    omega::Real=0.67,
    iterations::Int=2,
    smoother::Symbol=:damped_jacobi,
) where {T<:AbstractFloat}
    Ny1 = op.Ny1
    Nx1 = op.Nx1
    x4_mat = reshape(x4, (4, Ny1, Nx1))
    inv_mat = reshape(inv_d, (4, Ny1, Nx1))
    omega_T = T(omega)
    b_stride = length(b) == 2 * Ny1 * Nx1 ? 2 : 4
    b_mat = reshape(b, (b_stride, Ny1, Nx1))

    if smoother == :redblack_gauss_seidel
        for _ in 1:iterations
            LinearAlgebra.mul!(res_buf, op, x4)
            res_mat = reshape(res_buf, (4, Ny1, Nx1))
            @inbounds for j in 1:Nx1, i in 1:Ny1
                if iseven(i + j)
                    if !is_boundary_vx(i, j, op.Ny_val, op.Nx_val, Ny1, Nx1)
                        rvx = b_mat[1, i, j] - res_mat[1, i, j]
                        x4_mat[1, i, j] += omega_T * inv_mat[1, i, j] * rvx
                    end
                    if !is_boundary_vy(i, j, op.Ny_val, op.Nx_val, Ny1, Nx1)
                        rvy = b_mat[2, i, j] - res_mat[2, i, j]
                        x4_mat[2, i, j] += omega_T * inv_mat[2, i, j] * rvy
                    end
                end
            end
            LinearAlgebra.mul!(res_buf, op, x4)
            res_mat = reshape(res_buf, (4, Ny1, Nx1))
            @inbounds for j in 1:Nx1, i in 1:Ny1
                if isodd(i + j)
                    if !is_boundary_vx(i, j, op.Ny_val, op.Nx_val, Ny1, Nx1)
                        rvx = b_mat[1, i, j] - res_mat[1, i, j]
                        x4_mat[1, i, j] += omega_T * inv_mat[1, i, j] * rvx
                    end
                    if !is_boundary_vy(i, j, op.Ny_val, op.Nx_val, Ny1, Nx1)
                        rvy = b_mat[2, i, j] - res_mat[2, i, j]
                        x4_mat[2, i, j] += omega_T * inv_mat[2, i, j] * rvy
                    end
                end
            end
        end
    else
        for _ in 1:iterations
            LinearAlgebra.mul!(res_buf, op, x4)
            res_mat = reshape(res_buf, (4, Ny1, Nx1))
            @inbounds for j in 1:Nx1, i in 1:Ny1
                if !is_boundary_vx(i, j, op.Ny_val, op.Nx_val, Ny1, Nx1)
                    rvx = b_mat[1, i, j] - res_mat[1, i, j]
                    x4_mat[1, i, j] += omega_T * inv_mat[1, i, j] * rvx
                end
                if !is_boundary_vy(i, j, op.Ny_val, op.Nx_val, Ny1, Nx1)
                    rvy = b_mat[2, i, j] - res_mat[2, i, j]
                    x4_mat[2, i, j] += omega_T * inv_mat[2, i, j] * rvy
                end
            end
        end
    end
    return x4
end

"""
    smooth_darcy!(x4, bpf, op, inv_d, res_buf; omega=0.67, iterations=2, smoother=:damped_jacobi)

Perform relaxation sweeps on the Darcy fluid pressure unknowns with boundary condition enforcement.
"""
function smooth_darcy!(
    x4::AbstractVector{T},
    bpf::AbstractVector{T},
    op::MatrixFreeStokesDarcyOperator{T},
    inv_d::AbstractVector{T},
    res_buf::AbstractVector{T};
    omega::Real=0.67,
    iterations::Int=2,
    smoother::Symbol=:damped_jacobi,
) where {T<:AbstractFloat}
    Ny1 = op.Ny1
    Nx1 = op.Nx1
    x4_mat = reshape(x4, (4, Ny1, Nx1))
    inv_mat = reshape(inv_d, (4, Ny1, Nx1))
    omega_T = T(omega)
    b_stride = length(bpf) == Ny1 * Nx1 ? 1 : 4
    b_mat = reshape(bpf, (b_stride, Ny1, Nx1))
    b_idx = b_stride == 1 ? 1 : 4

    if smoother == :redblack_gauss_seidel
        for _ in 1:iterations
            LinearAlgebra.mul!(res_buf, op, x4)
            res_mat = reshape(res_buf, (4, Ny1, Nx1))
            @inbounds for j in 1:Nx1, i in 1:Ny1
                if iseven(i + j) && !is_boundary_p(i, j, op.Ny_val, op.Nx_val, Ny1, Nx1)
                    rpf = b_mat[b_idx, i, j] - res_mat[4, i, j]
                    x4_mat[4, i, j] += omega_T * inv_mat[4, i, j] * rpf
                end
            end
            LinearAlgebra.mul!(res_buf, op, x4)
            res_mat = reshape(res_buf, (4, Ny1, Nx1))
            @inbounds for j in 1:Nx1, i in 1:Ny1
                if isodd(i + j) && !is_boundary_p(i, j, op.Ny_val, op.Nx_val, Ny1, Nx1)
                    rpf = b_mat[b_idx, i, j] - res_mat[4, i, j]
                    x4_mat[4, i, j] += omega_T * inv_mat[4, i, j] * rpf
                end
            end
        end
    else
        for _ in 1:iterations
            LinearAlgebra.mul!(res_buf, op, x4)
            res_mat = reshape(res_buf, (4, Ny1, Nx1))
            @inbounds for j in 1:Nx1, i in 1:Ny1
                if !is_boundary_p(i, j, op.Ny_val, op.Nx_val, Ny1, Nx1)
                    rpf = b_mat[b_idx, i, j] - res_mat[4, i, j]
                    x4_mat[4, i, j] += omega_T * inv_mat[4, i, j] * rpf
                end
            end
        end
    end
    return x4
end

"""
    v_cycle_velocity!(xv, bv, hierarchy, level_idx=1; pre_smooth=2, post_smooth=2, omega=0.67, smoother=:damped_jacobi)

Execute one recursive geometric multigrid V-cycle on the decoupled Stokes velocity operator.
"""
function v_cycle_velocity!(
    xv::AbstractVector{T},
    bv::AbstractVector{T},
    hierarchy::StaggeredGridHierarchy{T,L},
    level_idx::Int=1;
    pre_smooth::Int=2,
    post_smooth::Int=2,
    omega::Real=0.67,
    smoother::Symbol=:damped_jacobi,
) where {T<:AbstractFloat,L<:MultigridLevel{T}}
    levels = hierarchy.levels
    curr = levels[level_idx]
    num_levels = length(levels)
    Ny1 = curr.Ny1
    Nx1 = curr.Nx1

    is_2var = (length(xv) == 2 * Ny1 * Nx1)
    if is_2var
        x4 = curr.x_buf
        r4 = curr.r_buf
        fill!(x4, zero(T))
        fill!(r4, zero(T))
        x4_mat = reshape(x4, (4, Ny1, Nx1))
        r4_mat = reshape(r4, (4, Ny1, Nx1))
        xv_mat = reshape(xv, (2, Ny1, Nx1))
        bv_mat = reshape(bv, (2, Ny1, Nx1))
        x4_mat[1:2, :, :] .= xv_mat
        r4_mat[1:2, :, :] .= bv_mat
        b_use = r4
        x_use = x4
    else
        x_use = xv
        b_use = bv
    end

    Ny_val = curr.Ny
    Nx_val = curr.Nx

    # 1. Pre-smoothing
    smooth_velocity!(
        x_use,
        b_use,
        curr.op,
        curr.inv_diag,
        curr.res_buf;
        omega=omega,
        iterations=pre_smooth,
        smoother=smoother,
    )

    if level_idx == num_levels
        coarse_sweeps = max(pre_smooth + post_smooth, 16)
        smooth_velocity!(
            x_use,
            b_use,
            curr.op,
            curr.inv_diag,
            curr.res_buf;
            omega=omega,
            iterations=coarse_sweeps,
            smoother=smoother,
        )
        if is_2var
            xv_mat = reshape(xv, (2, Ny1, Nx1))
            x4_mat = reshape(x_use, (4, Ny1, Nx1))
            xv_mat .= x4_mat[1:2, :, :]
        end
        return xv
    end

    # 2. Residual computation: r = b - A * x directly into curr.r_buf
    LinearAlgebra.mul!(curr.res_buf, curr.op, x_use)
    res_mat = reshape(curr.res_buf, (4, Ny1, Nx1))
    r_mat = reshape(curr.r_buf, (4, Ny1, Nx1))
    b_mat = reshape(b_use, (4, Ny1, Nx1))

    @inbounds for j in 1:Nx1, i in 1:Ny1
        if !is_boundary_vx(i, j, Ny_val, Nx_val, Ny1, Nx1)
            r_mat[1, i, j] = b_mat[1, i, j] - res_mat[1, i, j]
        else
            r_mat[1, i, j] = zero(T)
        end
        if !is_boundary_vy(i, j, Ny_val, Nx_val, Ny1, Nx1)
            r_mat[2, i, j] = b_mat[2, i, j] - res_mat[2, i, j]
        else
            r_mat[2, i, j] = zero(T)
        end
        r_mat[3, i, j] = zero(T)
        r_mat[4, i, j] = zero(T)
    end

    # 3. Restrict residual to next level buffer
    next_level = levels[level_idx + 1]
    restrict_4var!(
        next_level.r_buf, curr.r_buf, next_level.Ny, next_level.Nx, curr.Ny, curr.Nx
    )

    # 4. Recursive V-cycle
    fill!(next_level.x_buf, zero(T))
    v_cycle_velocity!(
        next_level.x_buf,
        next_level.r_buf,
        hierarchy,
        level_idx + 1;
        pre_smooth=pre_smooth,
        post_smooth=post_smooth,
        omega=omega,
        smoother=smoother,
    )

    # 5. Prolongate correction directly into x_use
    prolongate_4var!(
        x_use, next_level.x_buf, curr.Ny, curr.Nx, next_level.Ny, next_level.Nx
    )

    # 6. Post-smoothing
    smooth_velocity!(
        x_use,
        b_use,
        curr.op,
        curr.inv_diag,
        curr.res_buf;
        omega=omega,
        iterations=post_smooth,
        smoother=smoother,
    )

    if is_2var
        xv_mat = reshape(xv, (2, Ny1, Nx1))
        x4_mat = reshape(x_use, (4, Ny1, Nx1))
        xv_mat .= x4_mat[1:2, :, :]
    end

    return xv
end

"""
    v_cycle_darcy!(xpf, bpf, hierarchy, level_idx=1; pre_smooth=2, post_smooth=2, omega=0.67, smoother=:damped_jacobi)

Execute one recursive geometric multigrid V-cycle on the decoupled Darcy fluid pressure operator.
"""
function v_cycle_darcy!(
    xpf::AbstractVector{T},
    bpf::AbstractVector{T},
    hierarchy::StaggeredGridHierarchy{T,L},
    level_idx::Int=1;
    pre_smooth::Int=2,
    post_smooth::Int=2,
    omega::Real=0.67,
    smoother::Symbol=:damped_jacobi,
) where {T<:AbstractFloat,L<:MultigridLevel{T}}
    levels = hierarchy.levels
    curr = levels[level_idx]
    num_levels = length(levels)
    Ny1 = curr.Ny1
    Nx1 = curr.Nx1

    is_1var = (length(xpf) == Ny1 * Nx1)
    if is_1var
        x4 = curr.x_buf
        r4 = curr.r_buf
        fill!(x4, zero(T))
        fill!(r4, zero(T))
        x4_mat = reshape(x4, (4, Ny1, Nx1))
        r4_mat = reshape(r4, (4, Ny1, Nx1))
        xpf_mat = reshape(xpf, (Ny1, Nx1))
        bpf_mat = reshape(bpf, (Ny1, Nx1))
        x4_mat[4, :, :] .= xpf_mat
        r4_mat[4, :, :] .= bpf_mat
        b_use = r4
        x_use = x4
    else
        x_use = xpf
        b_use = bpf
    end

    Ny_val = curr.Ny
    Nx_val = curr.Nx

    # 1. Pre-smoothing
    smooth_darcy!(
        x_use,
        b_use,
        curr.op,
        curr.inv_diag,
        curr.res_buf;
        omega=omega,
        iterations=pre_smooth,
        smoother=smoother,
    )

    if level_idx == num_levels
        coarse_sweeps = max(pre_smooth + post_smooth, 16)
        smooth_darcy!(
            x_use,
            b_use,
            curr.op,
            curr.inv_diag,
            curr.res_buf;
            omega=omega,
            iterations=coarse_sweeps,
            smoother=smoother,
        )
        if is_1var
            xpf_mat = reshape(xpf, (Ny1, Nx1))
            x4_mat = reshape(x_use, (4, Ny1, Nx1))
            xpf_mat .= x4_mat[4, :, :]
        end
        return xpf
    end

    # 2. Residual computation: r = b - A * x directly into curr.r_buf
    LinearAlgebra.mul!(curr.res_buf, curr.op, x_use)
    res_mat = reshape(curr.res_buf, (4, Ny1, Nx1))
    r_mat = reshape(curr.r_buf, (4, Ny1, Nx1))
    b_mat = reshape(b_use, (4, Ny1, Nx1))

    @inbounds for j in 1:Nx1, i in 1:Ny1
        r_mat[1, i, j] = zero(T)
        r_mat[2, i, j] = zero(T)
        r_mat[3, i, j] = zero(T)
        if !is_boundary_p(i, j, Ny_val, Nx_val, Ny1, Nx1)
            r_mat[4, i, j] = b_mat[4, i, j] - res_mat[4, i, j]
        else
            r_mat[4, i, j] = zero(T)
        end
    end

    # 3. Restrict residual to next level buffer
    next_level = levels[level_idx + 1]
    restrict_4var!(
        next_level.r_buf, curr.r_buf, next_level.Ny, next_level.Nx, curr.Ny, curr.Nx
    )

    # 4. Recursive V-cycle
    fill!(next_level.x_buf, zero(T))
    v_cycle_darcy!(
        next_level.x_buf,
        next_level.r_buf,
        hierarchy,
        level_idx + 1;
        pre_smooth=pre_smooth,
        post_smooth=post_smooth,
        omega=omega,
        smoother=smoother,
    )

    # 5. Prolongate correction directly into x_use
    prolongate_4var!(
        x_use, next_level.x_buf, curr.Ny, curr.Nx, next_level.Ny, next_level.Nx
    )

    # 6. Post-smoothing
    smooth_darcy!(
        x_use,
        b_use,
        curr.op,
        curr.inv_diag,
        curr.res_buf;
        omega=omega,
        iterations=post_smooth,
        smoother=smoother,
    )

    if is_1var
        xpf_mat = reshape(xpf, (Ny1, Nx1))
        x4_mat = reshape(x_use, (4, Ny1, Nx1))
        xpf_mat .= x4_mat[4, :, :]
    end

    return xpf
end

"""
    v_cycle!(x, b, hierarchy, level_idx=1; pre_smooth=2, post_smooth=2, omega=0.67, smoother=:damped_jacobi)

Execute one recursive geometric multigrid V-cycle on `hierarchy`.
"""
function v_cycle!(
    x::AbstractVector{T},
    b::AbstractVector{T},
    hierarchy::StaggeredGridHierarchy{T,L},
    level_idx::Int=1;
    pre_smooth::Int=2,
    post_smooth::Int=2,
    omega::Real=0.67,
    smoother::Symbol=:damped_jacobi,
) where {T<:AbstractFloat,L<:MultigridLevel{T}}
    levels = hierarchy.levels
    curr = levels[level_idx]

    # Base diagonal scaling across all components (Schur Pt scaling)
    @inbounds @simd for i in eachindex(x, b, curr.inv_diag)
        x[i] = curr.inv_diag[i] * b[i]
    end

    # Apply decoupled V-cycles: velocity on [vx, vy], darcy on Pf
    v_cycle_velocity!(
        x,
        b,
        hierarchy,
        level_idx;
        pre_smooth=pre_smooth,
        post_smooth=post_smooth,
        omega=omega,
        smoother=smoother,
    )

    # Darcy on Pf
    v_cycle_darcy!(
        x,
        b,
        hierarchy,
        level_idx;
        pre_smooth=pre_smooth,
        post_smooth=post_smooth,
        omega=omega,
        smoother=smoother,
    )

    return x
end

# LinearAlgebra interface for MultigridPreconditioner
function LinearAlgebra.ldiv!(
    y::AbstractVector, P::MultigridPreconditioner, x::AbstractVector
)
    backend_y = KernelAbstractions.get_backend(y)
    backend_x = KernelAbstractions.get_backend(x)
    backend_y == backend_x || throw(
        ArgumentError(
            "Backend mismatch in ldiv!: y is on $(typeof(backend_y)) while x is on $(typeof(backend_x))",
        ),
    )
    backend_p = KernelAbstractions.get_backend(P.hierarchy.levels[1].op.ETA)
    backend_y == backend_p || throw(
        ArgumentError(
            "Backend mismatch in ldiv!: preconditioner hierarchy arrays are on $(typeof(backend_p)) while vectors are on $(typeof(backend_y))",
        ),
    )
    if !(backend_y isa KernelAbstractions.CPU)
        fill!(y, zero(eltype(y)))
        return apply_multigrid_vcycle_device!(
            y,
            x,
            P.hierarchy;
            backend=backend_y,
            smoother=P.smoother,
            pre_smooth=P.pre_smooth,
            post_smooth=P.post_smooth,
            omega=P.omega,
        )
    end
    fill!(y, zero(eltype(y)))
    v_cycle!(
        y,
        x,
        P.hierarchy,
        1;
        pre_smooth=P.pre_smooth,
        post_smooth=P.post_smooth,
        omega=P.omega,
        smoother=P.smoother,
    )
    return y
end

function LinearAlgebra.ldiv!(P::MultigridPreconditioner, x::AbstractVector)
    tmp = similar(x)
    LinearAlgebra.ldiv!(tmp, P, x)
    x .= tmp
    return x
end

function LinearAlgebra.mul!(
    y::AbstractVector, P::MultigridPreconditioner, x::AbstractVector
)
    return LinearAlgebra.ldiv!(y, P, x)
end

function LinearAlgebra.mul!(
    y::AbstractVector,
    P::MultigridPreconditioner,
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
        tmp = similar(y)
        LinearAlgebra.mul!(tmp, P, x)
        y .= alpha .* tmp .+ beta .* y
    end
    return y
end

"""
    coarsen_property_harmonic(A_f::Matrix{T}, Ny1_c::Int, Nx1_c::Int; tol=1.0e-30) where {T}

Coarsen a 2D property array using harmonic cell averaging.
"""
function coarsen_property_harmonic(
    A_f::Matrix{T}, Ny1_c::Int, Nx1_c::Int; tol::Real=1.0e-30
) where {T<:AbstractFloat}
    A_c = zeros(T, Ny1_c, Nx1_c)
    Ny1_f, Nx1_f = size(A_f)

    @inbounds for j in 1:Nx1_c, i in 1:Ny1_c
        if i == 1 || i == Ny1_c || j == 1 || j == Nx1_c
            fi = clamp(2 * i - 1, 1, Ny1_f)
            fj = clamp(2 * j - 1, 1, Nx1_f)
            A_c[i, j] = A_f[fi, fj]
        else
            fi = 2 * i - 1
            fj = 2 * j - 1
            v1 = A_f[fi - 1, fj - 1]
            v2 = A_f[fi - 1, fj]
            v3 = A_f[fi, fj - 1]
            v4 = A_f[fi, fj]
            inv_sum =
                (abs(v1) > tol ? inv(v1) : one(T)) +
                (abs(v2) > tol ? inv(v2) : one(T)) +
                (abs(v3) > tol ? inv(v3) : one(T)) +
                (abs(v4) > tol ? inv(v4) : one(T))
            A_c[i, j] = inv(0.25 * inv_sum)
        end
    end
    return A_c
end

"""
    coarsen_property_arithmetic(A_f::Matrix{T}, Ny1_c::Int, Nx1_c::Int) where {T}

Coarsen a 2D property array using arithmetic cell averaging.
"""
function coarsen_property_arithmetic(
    A_f::Matrix{T}, Ny1_c::Int, Nx1_c::Int
) where {T<:AbstractFloat}
    A_c = zeros(T, Ny1_c, Nx1_c)
    Ny1_f, Nx1_f = size(A_f)

    @inbounds for j in 1:Nx1_c, i in 1:Ny1_c
        if i == 1 || i == Ny1_c || j == 1 || j == Nx1_c
            fi = clamp(2 * i - 1, 1, Ny1_f)
            fj = clamp(2 * j - 1, 1, Nx1_f)
            A_c[i, j] = A_f[fi, fj]
        else
            fi = 2 * i - 1
            fj = 2 * j - 1
            A_c[i, j] =
                0.25 *
                (A_f[fi - 1, fj - 1] + A_f[fi - 1, fj] + A_f[fi, fj - 1] + A_f[fi, fj])
        end
    end
    return A_c
end

"""
    compute_multigrid_inv_diag(op; tol=1.0e-30)

Compute inverse diagonal scaling for multigrid levels, including Schur complement scaling for Pt.
"""
function compute_multigrid_inv_diag(
    op::MatrixFreeStokesDarcyOperator{T}; tol::Real=1.0e-30
) where {T<:AbstractFloat}
    backend_op = KernelAbstractions.get_backend(op.ETA)
    d_raw = compute_operator_diagonal(op)
    d_host = to_host(d_raw)
    inv_d_host = [abs(v) > tol ? inv(v) : one(T) for v in d_host]
    inv_d_mat = reshape(inv_d_host, (4, op.Ny1, op.Nx1))
    d_mat = reshape(d_host, (4, op.Ny1, op.Nx1))
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
    return to_device(backend_op, inv_d_host)
end

"""
    build_staggered_multigrid_hierarchy(op_fine; max_levels=4)

Build a `StaggeredGridHierarchy` from a fine-grid `MatrixFreeStokesDarcyOperator`.
"""
function build_staggered_multigrid_hierarchy(
    op_fine::MatrixFreeStokesDarcyOperator{T,M}; max_levels::Int=4
) where {T<:AbstractFloat,M<:AbstractMatrix{T}}
    backend_op = KernelAbstractions.get_backend(op_fine.ETA)
    if !(backend_op isa KernelAbstractions.CPU)
        op_host = to_host(op_fine)
        h_host = build_staggered_multigrid_hierarchy(op_host; max_levels=max_levels)
        return to_device(backend_op, h_host)
    end
    LType = MultigridLevel{T,Vector{T},MatrixFreeStokesDarcyOperator{T,M}}
    levels = LType[]

    # Level 1: fine grid
    inv_d_fine = compute_multigrid_inv_diag(op_fine)
    d_fine_len = op_fine.Ny1 * op_fine.Nx1 * 4
    push!(
        levels,
        MultigridLevel(
            1,
            op_fine.Nx_val,
            op_fine.Ny_val,
            op_fine.Nx1,
            op_fine.Ny1,
            op_fine.dx,
            op_fine.dy,
            op_fine,
            inv_d_fine,
            zeros(T, d_fine_len),
            zeros(T, d_fine_len),
            zeros(T, d_fine_len),
        ),
    )

    curr_op = op_fine
    curr_level = 1

    # Recursively build coarser levels while dimensions are even and >= 8
    while curr_level < max_levels &&
          curr_op.Nx_val >= 8 &&
          curr_op.Ny_val >= 8 &&
          iseven(curr_op.Nx_val) &&
          iseven(curr_op.Ny_val)
        Nx_c = curr_op.Nx_val ÷ 2
        Ny_c = curr_op.Ny_val ÷ 2
        Nx1_c = Nx_c + 1
        Ny1_c = Ny_c + 1
        dx_c = curr_op.dx * 2
        dy_c = curr_op.dy * 2

        ETA_c = coarsen_property_harmonic(curr_op.ETA, Ny1_c, Nx1_c)
        ETAP_c = coarsen_property_harmonic(curr_op.ETAP, Ny1_c, Nx1_c)
        GGG_c = coarsen_property_harmonic(curr_op.GGG, Ny1_c, Nx1_c)
        GGGP_c = coarsen_property_harmonic(curr_op.GGGP, Ny1_c, Nx1_c)
        RHOX_c = coarsen_property_arithmetic(curr_op.RHOX, Ny1_c, Nx1_c)
        RHOY_c = coarsen_property_arithmetic(curr_op.RHOY, Ny1_c, Nx1_c)
        RHOFX_c = coarsen_property_arithmetic(curr_op.RHOFX, Ny1_c, Nx1_c)
        RHOFY_c = coarsen_property_arithmetic(curr_op.RHOFY, Ny1_c, Nx1_c)
        RX_c = coarsen_property_harmonic(curr_op.RX, Ny1_c, Nx1_c)
        RY_c = coarsen_property_harmonic(curr_op.RY, Ny1_c, Nx1_c)
        ETAPHI_c = coarsen_property_harmonic(curr_op.ETAPHI, Ny1_c, Nx1_c)
        BETAPHI_c = coarsen_property_arithmetic(curr_op.BETAPHI, Ny1_c, Nx1_c)
        PHI_c = coarsen_property_arithmetic(curr_op.PHI, Ny1_c, Nx1_c)
        gx_c = coarsen_property_arithmetic(curr_op.gx, Ny1_c, Nx1_c)
        gy_c = coarsen_property_arithmetic(curr_op.gy, Ny1_c, Nx1_c)

        op_c = MatrixFreeStokesDarcyOperator{T}(
            Ny1_c,
            Nx1_c,
            dx_c,
            dy_c,
            Nx_c,
            Ny_c,
            ETA_c,
            ETAP_c,
            GGG_c,
            GGGP_c,
            RHOX_c,
            RHOY_c,
            RHOFX_c,
            RHOFY_c,
            RX_c,
            RY_c,
            ETAPHI_c,
            BETAPHI_c,
            PHI_c,
            gx_c,
            gy_c,
            curr_op.dt,
            curr_op.betasolid,
            curr_op.betafluid,
            curr_op.phimin,
            curr_op.phimax,
            curr_op.Kcont,
            curr_op.bctop,
            curr_op.bcbottom,
            curr_op.bcleft,
            curr_op.bcright,
        )

        inv_d_c = compute_multigrid_inv_diag(op_c)
        d_c_len = Ny1_c * Nx1_c * 4

        curr_level += 1
        push!(
            levels,
            MultigridLevel(
                curr_level,
                Nx_c,
                Ny_c,
                Nx1_c,
                Ny1_c,
                dx_c,
                dy_c,
                op_c,
                inv_d_c,
                zeros(T, d_c_len),
                zeros(T, d_c_len),
                zeros(T, d_c_len),
            ),
        )
        curr_op = op_c
    end

    if length(levels) < max_levels
        @warn "Staggered grid hierarchy stopped at level $(length(levels)) with coarsest resolution $(levels[end].Nx)x$(levels[end].Ny) cells"
    end

    return StaggeredGridHierarchy{T,LType}(levels)
end

"""
    build_multigrid_preconditioner(op_fine; levels=4, pre_smooth=2, post_smooth=2, omega=0.67, smoother=:damped_jacobi, dof_stride=4)

Build a `MultigridPreconditioner` wrapping a geometric multigrid hierarchy.
"""
function build_multigrid_preconditioner(
    op_fine::MatrixFreeStokesDarcyOperator{T,M};
    levels::Int=4,
    pre_smooth::Int=2,
    post_smooth::Int=2,
    omega::Real=0.67,
    smoother::Symbol=:damped_jacobi,
    dof_stride::Int=4,
) where {T<:AbstractFloat,M<:AbstractMatrix{T}}
    dof_stride == 4 || throw(
        ArgumentError(
            "MultigridPreconditioner requires dof_stride == 4, got $(dof_stride)"
        ),
    )
    hierarchy = build_staggered_multigrid_hierarchy(op_fine; max_levels=levels)
    return MultigridPreconditioner(hierarchy, pre_smooth, post_smooth, T(omega), smoother)
end

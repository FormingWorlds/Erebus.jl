# Unit and integration tests for geometric multigrid hierarchy, transfer operators, and solvers.

using Test
using LinearAlgebra
using Erebus

@testset "Geometric Multigrid (GMG) Solvers and Preconditioners" begin
    # 1. Setup synthetic physical properties on 64x64 grid
    Nx = 64
    Ny = 64
    Nx1 = Nx + 1
    Ny1 = Ny + 1
    xsize = 100_000.0
    ysize = 100_000.0
    coords = GridCoordinates(Nx, Ny; xsize=xsize, ysize=ysize)
    dx = coords.dx
    dy = coords.dy

    ETA = fill(1.0e19, Ny1, Nx1)
    ETAP = fill(1.0e19, Ny1, Nx1)
    GGG = fill(1.0e10, Ny1, Nx1)
    GGGP = fill(1.0e10, Ny1, Nx1)
    RHOX = fill(3300.0, Ny1, Nx1)
    RHOY = fill(3300.0, Ny1, Nx1)
    RHOFX = fill(1000.0, Ny1, Nx1)
    RHOFY = fill(1000.0, Ny1, Nx1)
    RX = fill(1.0e10, Ny1, Nx1)
    RY = fill(1.0e10, Ny1, Nx1)
    ETAPHI = fill(1.0e19, Ny1, Nx1)
    BETAPHI = fill(1.0e-11, Ny1, Nx1)
    PHI = fill(0.01, Ny1, Nx1)
    gx = zeros(Ny1, Nx1)
    gy = fill(-9.81, Ny1, Nx1)
    dt = 1000.0 * 365.25 * 86400.0

    op_fine = MatrixFreeStokesDarcyOperator(
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
        dt;
        coords=coords,
    )

    @testset "Hierarchy Construction and Geometric Scaling" begin
        hierarchy = build_staggered_multigrid_hierarchy(op_fine; max_levels=4)
        @test length(hierarchy.levels) == 4

        # Verify geometric dimension coarsening
        @test hierarchy.levels[1].Nx == 64
        @test hierarchy.levels[1].Ny == 64
        @test hierarchy.levels[2].Nx == 32
        @test hierarchy.levels[2].Ny == 32
        @test hierarchy.levels[3].Nx == 16
        @test hierarchy.levels[3].Ny == 16
        @test hierarchy.levels[4].Nx == 8
        @test hierarchy.levels[4].Ny == 8

        # Verify grid spacing scaling
        @test isapprox(hierarchy.levels[1].dx, dx)
        @test isapprox(hierarchy.levels[2].dx, 2.0 * dx)
        @test isapprox(hierarchy.levels[3].dx, 4.0 * dx)
        @test isapprox(hierarchy.levels[4].dx, 8.0 * dx)

        # Verify work buffers
        for lvl in hierarchy.levels
            expected_len = lvl.Ny1 * lvl.Nx1 * 4
            @test length(lvl.x_buf) == expected_len
            @test length(lvl.r_buf) == expected_len
            @test length(lvl.res_buf) == expected_len
            @test length(lvl.inv_diag) == expected_len
            @test all(isfinite, lvl.inv_diag)
            @test all(!iszero, lvl.inv_diag)
        end
    end

    @testset "Restriction and Prolongation Invariants" begin
        hierarchy = build_staggered_multigrid_hierarchy(op_fine; max_levels=3)
        lvl1 = hierarchy.levels[1]
        lvl2 = hierarchy.levels[2]

        # 1. Constant pressure field restriction preservation in interior
        rf = zeros(lvl1.Ny1 * lvl1.Nx1 * 4)
        rf_mat = reshape(rf, (4, lvl1.Ny1, lvl1.Nx1))
        rf_mat[3, :, :] .= 1.5e6
        rf_mat[4, :, :] .= 0.8e6

        rc = zeros(lvl2.Ny1 * lvl2.Nx1 * 4)
        restrict_4var!(rc, rf, lvl2.Ny, lvl2.Nx, lvl1.Ny, lvl1.Nx)
        rc_mat = reshape(rc, (4, lvl2.Ny1, lvl2.Nx1))

        # Interior coarse pressures must match constant value
        @test isapprox(rc_mat[3, 4, 4], 1.5e6; atol=1e-10)
        @test isapprox(rc_mat[4, 4, 4], 0.8e6; atol=1e-10)

        # 2. Prolongation coverage of boundary-adjacent interior nodes (row 2, column 2)
        ec_ones = ones(lvl2.Ny1 * lvl2.Nx1 * 4)
        ef_cov = zeros(lvl1.Ny1 * lvl1.Nx1 * 4)
        prolongate_4var!(ef_cov, ec_ones, lvl1.Ny, lvl1.Nx, lvl2.Ny, lvl2.Nx)
        ef_cov_mat = reshape(ef_cov, (4, lvl1.Ny1, lvl1.Nx1))

        # Verify fine row 2 and column 2 are covered for velocity interior nodes
        @test ef_cov_mat[1, 2, 2] > 0.0
        @test ef_cov_mat[2, 2, 2] > 0.0
        # Verify interior pressure nodes are covered
        @test ef_cov_mat[3, 5, 5] > 0.0
        @test ef_cov_mat[4, 5, 5] > 0.0

        # 3. Exact adjointness: <R u, v>_c == 1/4 <u, P v>_f for arbitrary vectors across full domain
        u_f = randn(lvl1.Ny1 * lvl1.Nx1 * 4)
        v_c = randn(lvl2.Ny1 * lvl2.Nx1 * 4)

        Ru = zeros(lvl2.Ny1 * lvl2.Nx1 * 4)
        restrict_4var!(Ru, u_f, lvl2.Ny, lvl2.Nx, lvl1.Ny, lvl1.Nx)

        Pv = zeros(lvl1.Ny1 * lvl1.Nx1 * 4)
        prolongate_4var!(Pv, v_c, lvl1.Ny, lvl1.Nx, lvl2.Ny, lvl2.Nx)

        dot_c = dot(Ru, v_c)
        dot_f = 0.25 * dot(u_f, Pv)
        @test isapprox(dot_c, dot_f; rtol=1e-12, atol=1e-14)
    end

    @testset "Damped Jacobi and Red-Black Smoothing" begin
        hierarchy = build_staggered_multigrid_hierarchy(op_fine; max_levels=2)
        lvl1 = hierarchy.levels[1]

        # 1. Damped Jacobi high-frequency attenuation
        b = zeros(lvl1.Ny1 * lvl1.Nx1 * 4)
        x_jacobi = zeros(lvl1.Ny1 * lvl1.Nx1 * 4)
        x_mat = reshape(x_jacobi, (4, lvl1.Ny1, lvl1.Nx1))
        for j in 3:(lvl1.Nx1 - 2), i in 3:(lvl1.Ny1 - 2)
            x_mat[1, i, j] = sin(16 * pi * i / lvl1.Ny1) * sin(16 * pi * j / lvl1.Nx1)
        end

        res_buf = similar(b)
        LinearAlgebra.mul!(res_buf, lvl1.op, x_jacobi)
        initial_norm = norm(res_buf)
        @test initial_norm > 0.0

        smooth_damped_jacobi!(
            x_jacobi, b, lvl1.op, lvl1.inv_diag, res_buf; omega=0.67, iterations=5
        )
        LinearAlgebra.mul!(res_buf, lvl1.op, x_jacobi)
        final_norm_jacobi = norm(res_buf)
        @test final_norm_jacobi < initial_norm

        # 2. Red-Black Gauss-Seidel velocity relaxation
        x_rb = copy(x_jacobi)
        x_mat_rb = reshape(x_rb, (4, lvl1.Ny1, lvl1.Nx1))
        for j in 3:(lvl1.Nx1 - 2), i in 3:(lvl1.Ny1 - 2)
            x_mat_rb[1, i, j] = sin(16 * pi * i / lvl1.Ny1) * sin(16 * pi * j / lvl1.Nx1)
        end
        LinearAlgebra.mul!(res_buf, lvl1.op, x_rb)
        initial_norm_rb = norm(res_buf)

        Erebus.smooth_velocity!(
            x_rb,
            b,
            lvl1.op,
            lvl1.inv_diag,
            res_buf;
            omega=0.67,
            iterations=5,
            smoother=:redblack_gauss_seidel,
        )
        LinearAlgebra.mul!(res_buf, lvl1.op, x_rb)
        final_norm_rb = norm(res_buf)
        @test final_norm_rb < initial_norm_rb
    end

    @testset "Decoupled Velocity and Darcy V-Cycles" begin
        hierarchy = build_staggered_multigrid_hierarchy(op_fine; max_levels=2)
        lvl1 = hierarchy.levels[1]

        # 1. Velocity V-Cycle attenuation
        bv = randn(2 * lvl1.Ny1 * lvl1.Nx1)
        bv_mat = reshape(bv, (2, lvl1.Ny1, lvl1.Nx1))
        for j in 1:lvl1.Nx1, i in 1:lvl1.Ny1
            if Erebus.is_boundary_vx(i, j, lvl1.Ny, lvl1.Nx, lvl1.Ny1, lvl1.Nx1)
                bv_mat[1, i, j] = 0.0
            end
            if Erebus.is_boundary_vy(i, j, lvl1.Ny, lvl1.Nx, lvl1.Ny1, lvl1.Nx1)
                bv_mat[2, i, j] = 0.0
            end
        end

        xv = zeros(2 * lvl1.Ny1 * lvl1.Nx1)
        x4 = zeros(4 * lvl1.Ny1 * lvl1.Nx1)
        x4_mat = reshape(x4, (4, lvl1.Ny1, lvl1.Nx1))
        y4 = zeros(4 * lvl1.Ny1 * lvl1.Nx1)
        y4_mat = reshape(y4, (4, lvl1.Ny1, lvl1.Nx1))

        x4_mat[1:2, :, :] .= bv_mat
        LinearAlgebra.mul!(y4, lvl1.op, x4)
        r0_v = norm(bv)

        v_cycle_velocity!(
            xv,
            bv,
            hierarchy,
            1;
            pre_smooth=3,
            post_smooth=3,
            omega=0.67,
            smoother=:damped_jacobi,
        )
        x4_mat[1:2, :, :] .= reshape(xv, (2, lvl1.Ny1, lvl1.Nx1))
        LinearAlgebra.mul!(y4, lvl1.op, x4)
        rf_v = norm(bv_mat .- y4_mat[1:2, :, :])
        @test rf_v < r0_v

        # Red-Black GS variant
        xv_rb = zeros(2 * lvl1.Ny1 * lvl1.Nx1)
        v_cycle_velocity!(
            xv_rb,
            bv,
            hierarchy,
            1;
            pre_smooth=3,
            post_smooth=3,
            omega=0.67,
            smoother=:redblack_gauss_seidel,
        )
        x4_mat[1:2, :, :] .= reshape(xv_rb, (2, lvl1.Ny1, lvl1.Nx1))
        LinearAlgebra.mul!(y4, lvl1.op, x4)
        rf_v_rb = norm(bv_mat .- y4_mat[1:2, :, :])
        @test rf_v_rb < r0_v

        # 2. Darcy Fluid Pressure V-Cycle attenuation
        bpf = randn(lvl1.Ny1 * lvl1.Nx1)
        bpf_mat = reshape(bpf, (lvl1.Ny1, lvl1.Nx1))
        for j in 1:lvl1.Nx1, i in 1:lvl1.Ny1
            if Erebus.is_boundary_p(i, j, lvl1.Ny, lvl1.Nx, lvl1.Ny1, lvl1.Nx1)
                bpf_mat[i, j] = 0.0
            end
        end

        xpf = zeros(lvl1.Ny1 * lvl1.Nx1)
        r0_pf = norm(bpf)
        v_cycle_darcy!(
            xpf,
            bpf,
            hierarchy,
            1;
            pre_smooth=3,
            post_smooth=3,
            omega=0.67,
            smoother=:damped_jacobi,
        )
        x4_mat .= 0.0
        x4_mat[4, :, :] .= reshape(xpf, (lvl1.Ny1, lvl1.Nx1))
        LinearAlgebra.mul!(y4, lvl1.op, x4)
        rf_pf = norm(bpf_mat .- y4_mat[4, :, :])
        @test rf_pf < r0_pf
    end

    @testset "MultigridPreconditioner Interface" begin
        hierarchy = build_staggered_multigrid_hierarchy(op_fine; max_levels=2)

        for sm in (:damped_jacobi, :redblack_gauss_seidel)
            P_mg = MultigridPreconditioner(hierarchy, 2, 2, 0.67, sm)

            v_test = randn(op_fine.Ny1 * op_fine.Nx1 * 4)
            y_test = zeros(length(v_test))

            LinearAlgebra.ldiv!(y_test, P_mg, v_test)
            @test norm(y_test) > 0.0
            @test all(isfinite, y_test)

            y_mul = zeros(length(v_test))
            LinearAlgebra.mul!(y_mul, P_mg, v_test)
            @test norm(y_mul - y_test) < 1.0e-14

            y_5arg = copy(v_test)
            LinearAlgebra.mul!(y_5arg, P_mg, v_test, 2.0, 3.0)
            @test isapprox(y_5arg, 2.0 .* y_test .+ 3.0 .* v_test; atol=1e-12)
        end
    end

    @testset "End-to-End Multigrid FGMRES Solve and Mesh Independence" begin
        # 32x32 Grid Solve
        coords_32 = GridCoordinates(32, 32; xsize=xsize, ysize=ysize)
        op_32 = MatrixFreeStokesDarcyOperator(
            fill(1.0e19, coords_32.Ny1, coords_32.Nx1),
            fill(1.0e19, coords_32.Ny1, coords_32.Nx1),
            fill(1.0e10, coords_32.Ny1, coords_32.Nx1),
            fill(1.0e10, coords_32.Ny1, coords_32.Nx1),
            fill(3300.0, coords_32.Ny1, coords_32.Nx1),
            fill(3300.0, coords_32.Ny1, coords_32.Nx1),
            fill(1000.0, coords_32.Ny1, coords_32.Nx1),
            fill(1000.0, coords_32.Ny1, coords_32.Nx1),
            fill(1.0e10, coords_32.Ny1, coords_32.Nx1),
            fill(1.0e10, coords_32.Ny1, coords_32.Nx1),
            fill(1.0e19, coords_32.Ny1, coords_32.Nx1),
            fill(1.0e-11, coords_32.Ny1, coords_32.Nx1),
            fill(0.01, coords_32.Ny1, coords_32.Nx1),
            zeros(coords_32.Ny1, coords_32.Nx1),
            fill(-9.81, coords_32.Ny1, coords_32.Nx1),
            dt;
            coords=coords_32,
        )

        x_true_32 = zeros(op_32.Ny1 * op_32.Nx1 * 4)
        x_mat_32 = reshape(x_true_32, (4, op_32.Ny1, op_32.Nx1))
        for j in 1:op_32.Nx1, i in 1:op_32.Ny1
            if !Erebus.is_boundary_vx(i, j, 32, 32, op_32.Ny1, op_32.Nx1)
                x_mat_32[1, i, j] =
                    sin(2 * pi * i / op_32.Ny1) * cos(2 * pi * j / op_32.Nx1)
            end
            if !Erebus.is_boundary_vy(i, j, 32, 32, op_32.Ny1, op_32.Nx1)
                x_mat_32[2, i, j] =
                    cos(2 * pi * i / op_32.Ny1) * sin(2 * pi * j / op_32.Nx1)
            end
            if !Erebus.is_boundary_p(i, j, 32, 32, op_32.Ny1, op_32.Nx1)
                x_mat_32[3, i, j] = 1.0e6 + 1.0e5 * sin(pi * i / op_32.Ny1)
                x_mat_32[4, i, j] = 1.0e5 * cos(pi * j / op_32.Nx1)
            end
        end
        b_32 = zeros(length(x_true_32))
        LinearAlgebra.mul!(b_32, op_32, x_true_32)
        S_32 = zeros(length(x_true_32))

        _, stats_32 = solve_hydromechanical_iterative!(
            op_32,
            b_32,
            S_32;
            method=:fgmres,
            preconditioner=:multigrid,
            mg_levels=3,
            rtol=1.0e-5,
            atol=1.0e-10,
            maxiter=100,
            restart=50,
        )
        @test stats_32.solved
        @test stats_32.niter <= 40

        # 64x64 Grid Solve
        x_true_64 = zeros(op_fine.Ny1 * op_fine.Nx1 * 4)
        x_mat_64 = reshape(x_true_64, (4, op_fine.Ny1, op_fine.Nx1))
        for j in 1:op_fine.Nx1, i in 1:op_fine.Ny1
            if !Erebus.is_boundary_vx(
                i, j, op_fine.Ny_val, op_fine.Nx_val, op_fine.Ny1, op_fine.Nx1
            )
                x_mat_64[1, i, j] =
                    sin(2 * pi * i / op_fine.Ny1) * cos(2 * pi * j / op_fine.Nx1)
            end
            if !Erebus.is_boundary_vy(
                i, j, op_fine.Ny_val, op_fine.Nx_val, op_fine.Ny1, op_fine.Nx1
            )
                x_mat_64[2, i, j] =
                    cos(2 * pi * i / op_fine.Ny1) * sin(2 * pi * j / op_fine.Nx1)
            end
            if !Erebus.is_boundary_p(
                i, j, op_fine.Ny_val, op_fine.Nx_val, op_fine.Ny1, op_fine.Nx1
            )
                x_mat_64[3, i, j] = 1.0e6 + 1.0e5 * sin(pi * i / op_fine.Ny1)
                x_mat_64[4, i, j] = 1.0e5 * cos(pi * j / op_fine.Nx1)
            end
        end

        b_64 = zeros(length(x_true_64))
        LinearAlgebra.mul!(b_64, op_fine, x_true_64)

        S_64 = zeros(length(x_true_64))
        _, stats_64 = solve_hydromechanical_iterative!(
            op_fine,
            b_64,
            S_64;
            method=:fgmres,
            preconditioner=:multigrid,
            mg_levels=3,
            rtol=1.0e-5,
            atol=1.0e-10,
            maxiter=100,
            restart=50,
        )
        @test stats_64.solved
        @test stats_64.niter <= 50

        # Verify mesh-independent iteration scaling ratio
        iter_ratio = stats_64.niter / stats_32.niter
        @test iter_ratio < 2.0

        # Red-Black Gauss-Seidel solve verification
        S_rb = zeros(length(x_true_32))
        _, stats_rb = solve_hydromechanical_iterative!(
            op_32,
            b_32,
            S_rb;
            method=:fgmres,
            preconditioner=:multigrid,
            mg_levels=3,
            mg_smoother=:redblack_gauss_seidel,
            rtol=1.0e-5,
            atol=1.0e-10,
            maxiter=100,
            restart=50,
        )
        @test stats_rb.solved
        @test stats_rb.niter <= 40

        # Viscosity contrast test (10^3 ratio)
        ETA_jump = fill(1.0e19, coords_32.Ny1, coords_32.Nx1)
        ETA_jump[1:(coords_32.Ny1 ÷ 2), :] .= 1.0e22
        op_jump = MatrixFreeStokesDarcyOperator(
            ETA_jump,
            ETA_jump,
            fill(1.0e10, coords_32.Ny1, coords_32.Nx1),
            fill(1.0e10, coords_32.Ny1, coords_32.Nx1),
            fill(3300.0, coords_32.Ny1, coords_32.Nx1),
            fill(3300.0, coords_32.Ny1, coords_32.Nx1),
            fill(1000.0, coords_32.Ny1, coords_32.Nx1),
            fill(1000.0, coords_32.Ny1, coords_32.Nx1),
            fill(1.0e10, coords_32.Ny1, coords_32.Nx1),
            fill(1.0e10, coords_32.Ny1, coords_32.Nx1),
            fill(1.0e19, coords_32.Ny1, coords_32.Nx1),
            fill(1.0e-11, coords_32.Ny1, coords_32.Nx1),
            fill(0.01, coords_32.Ny1, coords_32.Nx1),
            zeros(coords_32.Ny1, coords_32.Nx1),
            fill(-9.81, coords_32.Ny1, coords_32.Nx1),
            dt;
            coords=coords_32,
        )
        b_jump = zeros(length(x_true_32))
        LinearAlgebra.mul!(b_jump, op_jump, x_true_32)
        S_jump = zeros(length(x_true_32))
        _, stats_jump = solve_hydromechanical_iterative!(
            op_jump,
            b_jump,
            S_jump;
            method=:fgmres,
            preconditioner=:multigrid,
            mg_levels=3,
            rtol=1.0e-5,
            atol=1.0e-10,
            maxiter=150,
            restart=50,
        )
        @test stats_jump.solved
        res_jump = zeros(length(x_true_32))
        LinearAlgebra.mul!(res_jump, op_jump, S_jump)
        @test norm(b_jump - res_jump) / norm(b_jump) < 1.0e-5
    end
end

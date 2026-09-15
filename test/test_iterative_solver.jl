using Erebus
using ExtendableSparse
using LinearAlgebra
using SparseArrays
using Test

@testset "Iterative Solvers and Matrix-Free Operators" begin
    cfg = load_config("configs/test_quick.toml")
    coords = GridCoordinates(
        cfg.grid.Nx, cfg.grid.Ny; xsize=cfg.grid.xsize, ysize=cfg.grid.ysize
    )
    Nx1 = coords.Nx1
    Ny1 = coords.Ny1
    dx_val = coords.dx
    dy_val = coords.dy

    # Material properties for benchmark
    ETA = fill(1.0e19, Ny1, Nx1)
    ETAP = fill(1.0e19, Ny1, Nx1)
    GGG = fill(1.0e10, Ny1, Nx1)
    GGGP = fill(1.0e10, Ny1, Nx1)
    SXY0 = zeros(Ny1, Nx1)
    SXX0 = zeros(Ny1, Nx1)
    RHOX = fill(3000.0, Ny1, Nx1)
    RHOY = fill(3000.0, Ny1, Nx1)
    RHOFX = fill(1000.0, Ny1, Nx1)
    RHOFY = fill(1000.0, Ny1, Nx1)
    RX = fill(1.0e8, Ny1, Nx1)
    RY = fill(1.0e8, Ny1, Nx1)
    ETAPHI = fill(1.0e19, Ny1, Nx1)
    BETAPHI = fill(1.0e-11, Ny1, Nx1)
    PHI = fill(0.1, Ny1, Nx1)
    gx = zeros(Ny1, Nx1)
    gy = fill(1.0, Ny1, Nx1)
    pr0 = fill(1.0e6, Ny1, Nx1)
    pf0 = fill(1.0e6, Ny1, Nx1)
    DMP = zeros(Ny1, Nx1)
    dt = 1.0e10

    R4 = zeros(Ny1 * Nx1 * 4)
    L4 = assemble_hydromechanical_4var_lse!(
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
        R4;
        coords=coords,
    )

    op = MatrixFreeStokesDarcyOperator(
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

    @testset "Matrix-Free Operator Properties and Stencil" begin
        dof_total = Ny1 * Nx1 * 4
        @test size(op) == (dof_total, dof_total)
        @test size(op, 1) == dof_total
        @test size(op, 2) == dof_total
        @test size(op, 3) == 1
        @test eltype(op) == Float64

        # 3-arg mul! vs assembled CSC matrix
        x_rand = rand(dof_total)
        y_csc = L4 * x_rand
        y_op = zeros(dof_total)
        mul!(y_op, op, x_rand)
        rel_mul_err = norm(y_op - y_csc) / norm(y_csc)
        @test rel_mul_err < 1.0e-15

        # 5-arg mul! (y = alpha * A * x + beta * y)
        alpha_val = 1.75
        beta_val = 0.5
        y_5arg_csc = alpha_val .* (L4 * x_rand) .+ beta_val .* y_csc
        y_5arg_op = copy(y_csc)
        mul!(y_5arg_op, op, x_rand, alpha_val, beta_val)
        rel_5arg_err = norm(y_5arg_op - y_5arg_csc) / norm(y_5arg_csc)
        @test rel_5arg_err < 1.0e-15

        # Diagonal extraction exactness
        d_csc = diag(L4)
        d_op = compute_operator_diagonal(op)
        rel_d_err = norm(d_op - d_csc) / norm(d_csc)
        @test rel_d_err < 1.0e-15
    end

    @testset "Preconditioner Construction and Application" begin
        # Diagonal Preconditioner
        P_diag = build_diagonal_preconditioner(op)
        @test isa(P_diag, DiagonalPreconditioner)
        @test length(P_diag.inv_diag) == Ny1 * Nx1 * 4

        v_test = rand(Ny1 * Nx1 * 4)
        y_test = zeros(Ny1 * Nx1 * 4)
        mul!(y_test, P_diag, v_test)
        @test norm(y_test - P_diag.inv_diag .* v_test) < 1.0e-14

        y_5arg = copy(v_test)
        mul!(y_5arg, P_diag, v_test, 2.0, 3.0)
        @test norm(y_5arg - (2.0 .* P_diag.inv_diag .* v_test .+ 3.0 .* v_test)) < 1.0e-14

        # In-place ldiv!
        v_copy = copy(v_test)
        ldiv!(P_diag, v_copy)
        @test norm(v_copy - P_diag.inv_diag .* v_test) < 1.0e-14

        # Block-Schur Preconditioner
        P_schur_mat = build_block_schur_preconditioner(L4; coords=coords)
        P_schur_op = build_block_schur_preconditioner(op)
        @test isa(P_schur_mat, BlockSchurPreconditioner)
        @test isa(P_schur_op, BlockSchurPreconditioner)
        @test norm(P_schur_mat.inv_diag - P_schur_op.inv_diag) /
              norm(P_schur_mat.inv_diag) < 1.0e-14
    end

    @testset "Krylov Solvers Convergence and Solution Accuracy" begin
        S_direct = L4 \ R4

        # GMRES with matrix-free operator and Block-Schur preconditioner
        S_gmres = zeros(Ny1 * Nx1 * 4)
        _, stats_gmres = solve_hydromechanical_iterative!(
            op,
            R4,
            S_gmres;
            coords=coords,
            method=:gmres,
            preconditioner=:block_schur,
            rtol=1.0e-6,
            atol=1.0e-10,
            maxiter=500,
            restart=50,
        )
        @test stats_gmres.solved
        res_gmres = norm(L4 * S_gmres - R4) / norm(R4)
        @test res_gmres < 1.0e-5
        rel_sol_gmres = norm(S_gmres - S_direct) / norm(S_direct)
        @test rel_sol_gmres < 1.0e-2

        # FGMRES with matrix-free operator
        S_fgmres = zeros(Ny1 * Nx1 * 4)
        _, stats_fgmres = solve_hydromechanical_iterative!(
            op,
            R4,
            S_fgmres;
            coords=coords,
            method=:fgmres,
            preconditioner=:block_schur,
            rtol=1.0e-6,
            atol=1.0e-10,
            maxiter=500,
            restart=50,
        )
        @test stats_fgmres.solved
        res_fgmres = norm(L4 * S_fgmres - R4) / norm(R4)
        @test res_fgmres < 1.0e-5

        # BiCGStab with matrix-free operator
        S_bicg = zeros(Ny1 * Nx1 * 4)
        _, stats_bicg = solve_hydromechanical_iterative!(
            op,
            R4,
            S_bicg;
            coords=coords,
            method=:bicgstab,
            preconditioner=:block_schur,
            rtol=1.0e-6,
            atol=1.0e-10,
            maxiter=500,
        )
        @test stats_bicg.solved
        res_bicg = norm(L4 * S_bicg - R4) / norm(R4)
        @test res_bicg < 1.0e-5

        # Iterative on assembled CSC matrix
        S_csc_iter = zeros(Ny1 * Nx1 * 4)
        _, stats_csc = solve_hydromechanical_iterative!(
            L4,
            R4,
            S_csc_iter;
            coords=coords,
            method=:gmres,
            preconditioner=:diagonal,
            rtol=1.0e-6,
            atol=1.0e-10,
            maxiter=500,
            restart=50,
        )
        @test stats_csc.solved
        res_csc = norm(L4 * S_csc_iter - R4) / norm(R4)
        @test res_csc < 1.0e-5
    end

    @testset "6-variable Block Schur Preconditioner" begin
        # Test that dof_stride == 6 modifies Schur complement for pressure without throwing
        R6 = zeros(Ny1 * Nx1 * 6)
        L6 = assemble_hydromechanical_lse!(
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
            R6;
            coords=coords,
        )
        prec6 = build_block_schur_preconditioner(L6; coords=coords, dof_stride=6)
        @test prec6.dof_stride == 6
        @test length(prec6.inv_diag) == Ny1 * Nx1 * 6
        # Apply preconditioner
        y6 = zeros(Ny1 * Nx1 * 6)
        x6 = ones(Ny1 * Nx1 * 6)
        LinearAlgebra.ldiv!(y6, prec6, x6)
        @test all(isfinite, y6)
        @test !all(iszero, y6)
    end

    @testset "Iterative Solver Error Handling and Non-Convergence" begin
        S_dummy = zeros(Ny1 * Nx1 * 4)
        @test_throws ArgumentError solve_hydromechanical_iterative!(
            op, R4, S_dummy; method=:unknown_method
        )
        @test_throws ArgumentError solve_hydromechanical_iterative!(
            op, R4, S_dummy; preconditioner=:invalid_precond
        )

        # Test that tight tolerance with 1 iteration reports unsolved
        S_fail = zeros(Ny1 * Nx1 * 4)
        _, stats_fail = solve_hydromechanical_iterative!(
            op,
            R4,
            S_fail;
            coords=coords,
            method=:gmres,
            rtol=1.0e-20,
            atol=1.0e-20,
            maxiter=1,
        )
        @test stats_fail.solved == false
    end
end

using Erebus
using Erebus.Numerics
using ExtendableSparse
using LinearAlgebra
using LinearSolve
using Random
using SparseArrays
using Test

@testset "Linear Solver Optimization & Factorization Caching" begin
    coords = GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
    Ny1, Nx1 = coords.Ny1, coords.Nx1

    @testset "Gravitational Poisson pre-factorization equivalence" begin
        # Test unified boundary predicate
        r_limit = min(coords.xcenter, coords.ycenter)
        @test is_gravitational_boundary(
            1, 1, Ny1, Nx1, coords.xp, coords.yp, coords.xcenter, coords.ycenter, r_limit
        )
        @test is_gravitational_boundary(
            Ny1, 1, Ny1, Nx1, coords.xp, coords.yp, coords.xcenter, coords.ycenter, r_limit
        )
        @test is_gravitational_boundary(
            1, Nx1, Ny1, Nx1, coords.xp, coords.yp, coords.xcenter, coords.ycenter, r_limit
        )
        @test is_gravitational_boundary(
            Ny1,
            Nx1,
            Ny1,
            Nx1,
            coords.xp,
            coords.yp,
            coords.xcenter,
            coords.ycenter,
            r_limit,
        )
        # Center node must be internal
        mid_y, mid_x = Ny1 ÷ 2, Nx1 ÷ 2
        @test !is_gravitational_boundary(
            mid_y,
            mid_x,
            Ny1,
            Nx1,
            coords.xp,
            coords.yp,
            coords.xcenter,
            coords.ycenter,
            r_limit,
        )

        # Verify LP matrix invariance to density field
        RHO_zero = zeros(Ny1, Nx1)
        RP_dummy = zeros(Ny1 * Nx1)
        LP_zero = assemble_gravitational_lse!(RHO_zero, RP_dummy; coords=coords)

        RHO = fill(3000.0, Ny1, Nx1)
        # Add spatial perturbation
        RHO[5:12, 5:12] .+= 500.0

        RP_orig = zeros(Ny1 * Nx1)
        LP_orig = assemble_gravitational_lse!(RHO, RP_orig; coords=coords)
        @test LP_zero.cscmatrix == LP_orig.cscmatrix
        SP_orig = LP_orig \ RP_orig

        # Test assemble_gravitational_rhs!
        RP_rhs = zeros(Ny1 * Nx1)
        assemble_gravitational_rhs!(RHO, RP_rhs; coords=coords)
        @test RP_rhs == RP_orig

        # Test pre-factorized solve
        F_grav = lu(LP_orig.cscmatrix)
        SP_fact = F_grav \ RP_rhs
        @test SP_fact == SP_orig

        # Test physical observables gx, gy
        FI_orig = zeros(Ny1, Nx1)
        gx_orig = zeros(Ny1, Nx1)
        gy_orig = zeros(Ny1, Nx1)
        process_gravitational_solution!(SP_orig, FI_orig, gx_orig, gy_orig; coords=coords)

        FI_fact = zeros(Ny1, Nx1)
        gx_fact = zeros(Ny1, Nx1)
        gy_fact = zeros(Ny1, Nx1)
        process_gravitational_solution!(SP_fact, FI_fact, gx_fact, gy_fact; coords=coords)

        @test FI_fact == FI_orig
        @test gx_fact == gx_orig
        @test gy_fact == gy_orig

        # Verify LP buffer reuse and ensure no entry doubling on repeated calls
        LP_buf = ExtendableSparseMatrix(Nx1 * Ny1, Nx1 * Ny1)
        LP_reused1 = assemble_gravitational_lse!(RHO, RP_rhs; coords=coords, LP=LP_buf)
        @test LP_reused1.cscmatrix == LP_orig.cscmatrix
        LP_reused2 = assemble_gravitational_lse!(RHO, RP_rhs; coords=coords, LP=LP_buf)
        @test LP_reused2.cscmatrix == LP_orig.cscmatrix
    end

    @testset "Hydromechanical LSE in-place assembly equivalence" begin
        ETA = fill(1e20, Ny1, Nx1)
        ETAP = fill(1e20, Ny1, Nx1)
        GGG = fill(1e10, Ny1, Nx1)
        GGGP = fill(1e10, Ny1, Nx1)
        SXY0 = zeros(Ny1, Nx1)
        SXX0 = zeros(Ny1, Nx1)
        RHOX = fill(3000.0, Ny1, Nx1)
        RHOY = fill(3000.0, Ny1, Nx1)
        RHOFX = fill(1000.0, Ny1, Nx1)
        RHOFY = fill(1000.0, Ny1, Nx1)
        RX = fill(1e-12, Ny1, Nx1)
        RY = fill(1e-12, Ny1, Nx1)
        ETAPHI = fill(1e20, Ny1, Nx1)
        BETAPHI = zeros(Ny1, Nx1)
        PHI = fill(0.1, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = fill(-1.0, Ny1, Nx1)
        pr0 = zeros(Ny1, Nx1)
        pf0 = zeros(Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)
        dt = 1e11

        R_fresh = zeros(Nx1 * Ny1 * 6)
        L_fresh = assemble_hydromechanical_lse!(
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
            R_fresh;
            coords=coords,
        )

        # In-place assembly with pre-allocated buffer
        L_buf = ExtendableSparseMatrix(Nx1 * Ny1 * 6, Nx1 * Ny1 * 6)
        R_buf = zeros(Nx1 * Ny1 * 6)
        L_reused = assemble_hydromechanical_lse!(
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
            R_buf;
            coords=coords,
            L=L_buf,
        )
        @test L_reused == L_fresh
        @test R_buf == R_fresh

        # Second in-place pass with mutated properties
        ETA .+= 5e19
        R_fresh2 = zeros(Nx1 * Ny1 * 6)
        L_fresh2 = assemble_hydromechanical_lse!(
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
            R_fresh2;
            coords=coords,
        )
        R_buf2 = zeros(Nx1 * Ny1 * 6)
        L_reused2 = assemble_hydromechanical_lse!(
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
            R_buf2;
            coords=coords,
            L=L_buf,
        )
        @test L_reused2 == L_fresh2
        @test R_buf2 == R_fresh2
    end

    @testset "Thermal LSE in-place assembly equivalence" begin
        tk1 = fill(300.0, Ny1, Nx1)
        RHOCP = fill(2e6, Ny1, Nx1)
        KX = fill(2.0, Ny1, Nx1)
        KY = fill(2.0, Ny1, Nx1)
        HR = fill(1e-7, Ny1, Nx1)
        HA = zeros(Ny1, Nx1)
        HS = zeros(Ny1, Nx1)
        DHP = zeros(Ny1, Nx1)
        RT_fresh = zeros(Ny1 * Nx1)
        dt = 1e11

        LT_fresh = assemble_thermal_lse!(
            tk1, RHOCP, KX, KY, HR, HA, HS, DHP, RT_fresh, dt; coords=coords
        )

        LT_buf = ExtendableSparseMatrix(Ny1 * Nx1, Ny1 * Nx1)
        RT_buf = zeros(Ny1 * Nx1)
        LT_reused = assemble_thermal_lse!(
            tk1, RHOCP, KX, KY, HR, HA, HS, DHP, RT_buf, dt; coords=coords, LT=LT_buf
        )
        @test LT_reused.cscmatrix == LT_fresh.cscmatrix
        @test RT_buf == RT_fresh

        # Second in-place pass with mutated properties
        tk1 .+= 20.0
        KX .+= 0.5
        RT_fresh2 = zeros(Ny1 * Nx1)
        LT_fresh2 = assemble_thermal_lse!(
            tk1, RHOCP, KX, KY, HR, HA, HS, DHP, RT_fresh2, dt; coords=coords
        )
        RT_buf2 = zeros(Ny1 * Nx1)
        LT_reused2 = assemble_thermal_lse!(
            tk1, RHOCP, KX, KY, HR, HA, HS, DHP, RT_buf2, dt; coords=coords, LT=LT_buf
        )
        @test LT_reused2.cscmatrix == LT_fresh2.cscmatrix
        @test RT_buf2 == RT_fresh2
    end

    @testset "LinearSolve symbolic factorization caching across Picard iterations" begin
        ETA = fill(1e20, Ny1, Nx1)
        ETAP = fill(1e20, Ny1, Nx1)
        GGG = fill(1e10, Ny1, Nx1)
        GGGP = fill(1e10, Ny1, Nx1)
        SXY0 = zeros(Ny1, Nx1)
        SXX0 = zeros(Ny1, Nx1)
        RHOX = fill(3000.0, Ny1, Nx1)
        RHOY = fill(3000.0, Ny1, Nx1)
        RHOFX = fill(1000.0, Ny1, Nx1)
        RHOFY = fill(1000.0, Ny1, Nx1)
        RX = fill(1e-12, Ny1, Nx1)
        RY = fill(1e-12, Ny1, Nx1)
        ETAPHI = fill(1e20, Ny1, Nx1)
        BETAPHI = zeros(Ny1, Nx1)
        PHI = fill(0.1, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = fill(-1.0, Ny1, Nx1)
        pr0 = zeros(Ny1, Nx1)
        pf0 = zeros(Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)
        dt = 1e11

        L_buf = ExtendableSparseMatrix(Nx1 * Ny1 * 6, Nx1 * Ny1 * 6)
        R = zeros(Nx1 * Ny1 * 6)
        L = assemble_hydromechanical_lse!(
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
            coords=coords,
            L=L_buf,
        )

        prob = LinearProblem(L, R)
        cache = init(prob, UMFPACKFactorization(; reuse_symbolic=true))
        sol1 = solve!(cache)
        @test LinearSolve.SciMLBase.successful_retcode(sol1)
        @test sol1.retcode == LinearSolve.SciMLBase.ReturnCode.Success

        sol1_ref = LinearSolve.solve(LinearProblem(L, R), UMFPACKFactorization()).u
        @test maximum(abs.(sol1.u - sol1_ref)) == 0.0

        # Simulate 3 successive Picard iterations with property updates
        for iter in 2:4
            ETA .+= 1e19 * iter
            assemble_hydromechanical_lse!(
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
                coords=coords,
                L=L_buf,
            )
            cache.A = L
            cache.b = R
            sol_cur = solve!(cache)
            @test LinearSolve.SciMLBase.successful_retcode(sol_cur)
            sol_ref = LinearSolve.solve(LinearProblem(L, R), UMFPACKFactorization()).u
            @test maximum(abs.(sol_cur.u - sol_ref)) < 1e-12
        end
    end
end

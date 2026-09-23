# Unit tests for GPU acceleration via KernelAbstractions.jl

using Test
using LinearAlgebra
using Random
using KernelAbstractions
using Erebus

@testset "GPU Acceleration via KernelAbstractions" begin
    # 1. Setup synthetic physical properties on 32x32 grid
    Nx = 32
    Ny = 32
    Nx1 = Nx + 1
    Ny1 = Ny + 1
    xsize = 100_000.0
    ysize = 100_000.0
    coords = GridCoordinates(Nx, Ny; xsize=xsize, ysize=ysize)

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

    backend = KernelAbstractions.CPU()

    @testset "Device Data Migration (to_device / to_host)" begin
        v = randn(Ny1 * Nx1 * 4)
        v_dev = Erebus.to_device(backend, v)
        @test size(v_dev) == size(v)
        @test isapprox(Erebus.to_host(v_dev), v; rtol=1e-14, atol=1e-14)

        op_dev = Erebus.to_device(backend, op)
        @test op_dev.Ny1 == op.Ny1
        @test op_dev.Nx1 == op.Nx1
        @test isapprox(Erebus.to_host(op_dev.ETA), op.ETA)
    end

    @testset "Matrix-Free Stencil Kernel Equivalence" begin
        x = randn(Ny1 * Nx1 * 4)
        y_cpu = zeros(Ny1 * Nx1 * 4)
        y_dev = zeros(Ny1 * Nx1 * 4)

        mul!(y_cpu, op, x)
        Erebus.mul_device!(y_dev, op, x; backend=backend)

        @test isapprox(y_cpu, y_dev; rtol=1e-12, atol=1e-14)
        @test norm(y_dev) > 0.0
        @test all(isfinite, y_dev)
    end

    @testset "Operator Diagonal Kernel Equivalence" begin
        d_cpu = compute_operator_diagonal(op)
        d_dev = Erebus.compute_operator_diagonal_device(op; backend=backend)

        @test isapprox(d_cpu, d_dev; rtol=1e-12, atol=1e-14)
        @test all(isfinite, d_dev)
        @test length(d_dev) == length(d_cpu)
    end

    @testset "Restriction and Prolongation Kernels" begin
        Ny_c = Ny ÷ 2
        Nx_c = Nx ÷ 2
        Ny1_c = Ny_c + 1
        Nx1_c = Nx_c + 1

        u_f = randn(Ny1 * Nx1 * 4)
        v_c = randn(Ny1_c * Nx1_c * 4)

        Ru_cpu = zeros(Ny1_c * Nx1_c * 4)
        Ru_dev = zeros(Ny1_c * Nx1_c * 4)
        restrict_4var!(Ru_cpu, u_f, Ny_c, Nx_c, Ny, Nx)
        Erebus.restrict_4var_device!(Ru_dev, u_f, Ny_c, Nx_c, Ny, Nx; backend=backend)
        @test isapprox(Ru_cpu, Ru_dev; rtol=1e-12, atol=1e-14)

        Pv_cpu = zeros(Ny1 * Nx1 * 4)
        Pv_dev = zeros(Ny1 * Nx1 * 4)
        prolongate_4var!(Pv_cpu, v_c, Ny, Nx, Ny_c, Nx_c)
        Erebus.prolongate_4var_device!(Pv_dev, v_c, Ny, Nx, Ny_c, Nx_c; backend=backend)
        @test isapprox(Pv_cpu, Pv_dev; rtol=1e-12, atol=1e-14)

        # Adjointness of device kernels
        dot_c = dot(Ru_dev, v_c)
        dot_f = 0.25 * dot(u_f, Pv_dev)
        @test isapprox(dot_c, dot_f; rtol=1e-12, atol=1e-14)
    end

    @testset "Relaxation Smoother Kernels" begin
        x_jacobi = randn(Ny1 * Nx1 * 4)
        b = randn(Ny1 * Nx1 * 4)
        inv_diag = inv.(compute_operator_diagonal(op))
        res_buf = zeros(Ny1 * Nx1 * 4)

        x_ref = copy(x_jacobi)
        smooth_damped_jacobi!(x_ref, b, op, inv_diag, res_buf; iterations=2)

        x_dev = copy(x_jacobi)
        Erebus.smooth_damped_jacobi_device!(
            x_dev, b, op, inv_diag, res_buf; backend=backend, iterations=2
        )
        @test isapprox(x_ref, x_dev; rtol=1e-12, atol=1e-14)
        @test all(isfinite, x_dev)
        @test norm(x_dev - x_jacobi) > 0.0

        # Decoupled velocity and Darcy smoothers
        x_vel_ref = copy(x_jacobi)
        Erebus.smooth_velocity!(x_vel_ref, b, op, inv_diag, res_buf; iterations=2)
        x_vel_dev = copy(x_jacobi)
        Erebus.smooth_velocity_device!(
            x_vel_dev, b, op, inv_diag, res_buf; backend=backend, iterations=2
        )
        @test isapprox(x_vel_ref, x_vel_dev; rtol=1e-12, atol=1e-14)

        x_darcy_ref = copy(x_jacobi)
        Erebus.smooth_darcy!(x_darcy_ref, b, op, inv_diag, res_buf; iterations=2)
        x_darcy_dev = copy(x_jacobi)
        Erebus.smooth_darcy_device!(
            x_darcy_dev, b, op, inv_diag, res_buf; backend=backend, iterations=2
        )
        @test isapprox(x_darcy_ref, x_darcy_dev; rtol=1e-12, atol=1e-14)
    end

    @testset "Multigrid Hierarchy Device Migration and V-Cycle" begin
        hierarchy = build_staggered_multigrid_hierarchy(op; max_levels=3)
        hierarchy_dev = Erebus.to_device(backend, hierarchy)

        @test length(hierarchy_dev.levels) == 3
        @test hierarchy_dev.levels[1].Ny1 == hierarchy.levels[1].Ny1

        b = randn(Ny1 * Nx1 * 4)
        x_cpu = zeros(Ny1 * Nx1 * 4)
        x_dev = zeros(Ny1 * Nx1 * 4)

        v_cycle!(x_cpu, b, hierarchy; smoother=:damped_jacobi)
        Erebus.apply_multigrid_vcycle_device!(
            x_dev, b, hierarchy_dev; backend=backend, smoother=:damped_jacobi
        )
        @test isapprox(x_cpu, x_dev; rtol=1e-12, atol=1e-14)
        @test all(isfinite, x_dev)
    end

    @testset "Heterogeneous Property Field Equivalence" begin
        ETA_het = 1.0e19 .* exp.(randn(Ny1, Nx1) .* 0.5)
        ETAP_het = copy(ETA_het)
        GGG_het = fill(1e10, Ny1, Nx1)
        GGGP_het = copy(GGG_het)
        RHOX_het = 3300.0 .+ 50.0 .* randn(Ny1, Nx1)
        RHOY_het = copy(RHOX_het)
        RHOFX_het = fill(1000.0, Ny1, Nx1)
        RHOFY_het = copy(RHOFX_het)
        RX_het = 1.0e10 .* exp.(randn(Ny1, Nx1) .* 0.5)
        RY_het = copy(RX_het)
        ETAPHI_het = copy(ETA_het)
        BETAPHI_het = fill(1e-11, Ny1, Nx1)
        PHI_het = fill(0.01, Ny1, Nx1)
        gx_het = zeros(Ny1, Nx1)
        gy_het = fill(-9.81, Ny1, Nx1)

        op_het = MatrixFreeStokesDarcyOperator(
            ETA_het,
            ETAP_het,
            GGG_het,
            GGGP_het,
            RHOX_het,
            RHOY_het,
            RHOFX_het,
            RHOFY_het,
            RX_het,
            RY_het,
            ETAPHI_het,
            BETAPHI_het,
            PHI_het,
            gx_het,
            gy_het,
            dt;
            coords=coords,
        )

        x_het = randn(Ny1 * Nx1 * 4)
        y_cpu_het = zeros(Ny1 * Nx1 * 4)
        y_dev_het = zeros(Ny1 * Nx1 * 4)

        mul!(y_cpu_het, op_het, x_het)
        Erebus.mul_device!(y_dev_het, op_het, x_het; backend=backend)
        @test isapprox(y_cpu_het, y_dev_het; rtol=1e-12, atol=1e-14)

        h_het = build_staggered_multigrid_hierarchy(op_het; max_levels=3)
        h_dev_het = Erebus.to_device(backend, h_het)
        b_het = randn(Ny1 * Nx1 * 4)
        xc_het = zeros(Ny1 * Nx1 * 4)
        xd_het = zeros(Ny1 * Nx1 * 4)

        v_cycle!(xc_het, b_het, h_het; smoother=:damped_jacobi)
        Erebus.apply_multigrid_vcycle_device!(
            xd_het, b_het, h_dev_het; backend=backend, smoother=:damped_jacobi
        )
        @test isapprox(xc_het, xd_het; rtol=1e-12, atol=1e-14)
        @test all(isfinite, xd_het)
    end

    @testset "Float32 Precision and Type Stability" begin
        ETA_f32 = fill(1.0f19, Ny1, Nx1)
        ETAP_f32 = fill(1.0f19, Ny1, Nx1)
        GGG_f32 = fill(1.0f10, Ny1, Nx1)
        GGGP_f32 = fill(1.0f10, Ny1, Nx1)
        RHOX_f32 = fill(3300.0f0, Ny1, Nx1)
        RHOY_f32 = fill(3300.0f0, Ny1, Nx1)
        RHOFX_f32 = fill(1000.0f0, Ny1, Nx1)
        RHOFY_f32 = fill(1000.0f0, Ny1, Nx1)
        RX_f32 = fill(1.0f10, Ny1, Nx1)
        RY_f32 = fill(1.0f10, Ny1, Nx1)
        ETAPHI_f32 = fill(1.0f19, Ny1, Nx1)
        BETAPHI_f32 = fill(1.0f-11, Ny1, Nx1)
        PHI_f32 = fill(0.01f0, Ny1, Nx1)
        gx_f32 = zeros(Float32, Ny1, Nx1)
        gy_f32 = fill(-9.81f0, Ny1, Nx1)
        dt_f32 = 1000.0f0 * 365.25f0 * 86400.0f0

        op_f32 = MatrixFreeStokesDarcyOperator(
            ETA_f32,
            ETAP_f32,
            GGG_f32,
            GGGP_f32,
            RHOX_f32,
            RHOY_f32,
            RHOFX_f32,
            RHOFY_f32,
            RX_f32,
            RY_f32,
            ETAPHI_f32,
            BETAPHI_f32,
            PHI_f32,
            gx_f32,
            gy_f32,
            dt_f32;
            coords=coords,
        )

        x_f32 = randn(Float32, Ny1 * Nx1 * 4)
        y_f32 = zeros(Float32, Ny1 * Nx1 * 4)

        Erebus.mul_device!(y_f32, op_f32, x_f32; backend=backend)
        @test y_f32 isa Vector{Float32}
        @test all(isfinite, y_f32)
        @test norm(y_f32) > 0.0f0

        d_f32 = Erebus.compute_operator_diagonal_device(op_f32; backend=backend)
        @test d_f32 isa Vector{Float32}
        @test all(isfinite, d_f32)

        y_cpu_f32 = zeros(Float32, Ny1 * Nx1 * 4)
        mul!(y_cpu_f32, op_f32, x_f32)
        @test isapprox(y_cpu_f32, y_f32; rtol=1.0f-5, atol=1.0f-6)

        d_cpu_f32 = compute_operator_diagonal(op_f32)
        @test isapprox(d_cpu_f32, d_f32; rtol=1.0f-5, atol=1.0f-6)

        cd_f32 = Erebus.compute_drained_compressibility(1.0f-11, 0.2f0, 1.0f-12)
        cd_f64 = Erebus.compute_drained_compressibility(1.0e-11, 0.2, 1.0e-12)
        @test cd_f32 isa Float32
        @test isapprox(Float64(cd_f32), cd_f64; rtol=1.0e-6)

        alpha_f32 = Erebus.compute_biot_willis_coefficient(1.0f-11, 2.0f-11)
        alpha_f64 = Erebus.compute_biot_willis_coefficient(1.0e-11, 2.0e-11)
        @test alpha_f32 isa Float32
        @test isapprox(Float64(alpha_f32), alpha_f64; rtol=1.0e-6)

        pt_eval = Erebus.evaluate_stokes_darcy_point(
            10,
            10,
            reshape(x_f32, (4, Ny1, Nx1)),
            op_f32.Ny1,
            op_f32.Nx1,
            op_f32.dx,
            op_f32.dy,
            op_f32.Nx_val,
            op_f32.Ny_val,
            op_f32.dt,
            op_f32.Kcont,
            op_f32.bctop,
            op_f32.bcbottom,
            op_f32.bcleft,
            op_f32.bcright,
            op_f32.betasolid,
            op_f32.betafluid,
            op_f32.phimin,
            op_f32.phimax,
            op_f32.ETA,
            op_f32.ETAP,
            op_f32.GGG,
            op_f32.GGGP,
            op_f32.RHOX,
            op_f32.RHOY,
            op_f32.RHOFX,
            op_f32.RHOFY,
            op_f32.RX,
            op_f32.RY,
            op_f32.ETAPHI,
            op_f32.BETAPHI,
            op_f32.PHI,
            op_f32.gx,
            op_f32.gy,
        )
        @test pt_eval isa NTuple{4,Float32}
        @test all(isfinite, pt_eval)
    end

    @testset "Red-Black Gauss-Seidel Smoother Equivalence" begin
        h_rb = build_staggered_multigrid_hierarchy(op; max_levels=3)
        h_rb_dev = Erebus.to_device(backend, h_rb)
        b_rb = randn(Ny1 * Nx1 * 4)

        x_rb_cpu = zeros(Ny1 * Nx1 * 4)
        x_rb_dev = zeros(Ny1 * Nx1 * 4)
        v_cycle!(x_rb_cpu, b_rb, h_rb; smoother=:redblack_gauss_seidel)
        Erebus.apply_multigrid_vcycle_device!(
            x_rb_dev, b_rb, h_rb_dev; backend=backend, smoother=:redblack_gauss_seidel
        )
        @test isapprox(x_rb_cpu, x_rb_dev; rtol=1.0e-12, atol=1.0e-14)
        @test all(isfinite, x_rb_dev)

        # Unsupported smoother must throw ArgumentError
        @test_throws ArgumentError Erebus.apply_multigrid_vcycle_device!(
            x_rb_dev, b_rb, h_rb_dev; backend=backend, smoother=:invalid_smoother
        )
    end

    @testset "Reduced-DOF Sub-System V-Cycles" begin
        h_sub = build_staggered_multigrid_hierarchy(op; max_levels=3)
        h_sub_dev = Erebus.to_device(backend, h_sub)

        # 2-variable velocity V-cycle
        xv_cpu = randn(2 * Ny1 * Nx1)
        bv_2 = randn(2 * Ny1 * Nx1)
        xv_dev = copy(xv_cpu)
        v_cycle_velocity!(xv_cpu, bv_2, h_sub, 1; pre_smooth=2, post_smooth=2)
        Erebus.v_cycle_velocity_device!(
            xv_dev, bv_2, h_sub_dev, 1; backend=backend, pre_smooth=2, post_smooth=2
        )
        @test isapprox(xv_cpu, xv_dev; rtol=1.0e-12, atol=1.0e-14)

        # 1-variable Darcy pressure V-cycle
        xpf_cpu = randn(Ny1 * Nx1)
        bpf_1 = randn(Ny1 * Nx1)
        xpf_dev = copy(xpf_cpu)
        v_cycle_darcy!(xpf_cpu, bpf_1, h_sub, 1; pre_smooth=2, post_smooth=2)
        Erebus.v_cycle_darcy_device!(
            xpf_dev, bpf_1, h_sub_dev, 1; backend=backend, pre_smooth=2, post_smooth=2
        )
        @test isapprox(xpf_cpu, xpf_dev; rtol=1.0e-12, atol=1.0e-14)
    end

    @testset "Multigrid Preconditioner ldiv! and Memory Allocations" begin
        precond = build_multigrid_preconditioner(op; levels=3)
        @test isconcretetype(typeof(precond))
        @test isconcretetype(typeof(precond.hierarchy))

        b_vec = randn(Ny1 * Nx1 * 4)
        y_vec = zeros(Ny1 * Nx1 * 4)
        ldiv!(y_vec, precond, b_vec)
        @test all(isfinite, y_vec)
        @test norm(y_vec) > 0.0

        # Allocation verification: concrete CPU preconditioner avoids dynamic dispatch regression
        alloc_ldiv = @allocated ldiv!(y_vec, precond, b_vec)
        @test alloc_ldiv < 1024

        # Device preconditioner transfer
        p_dev = Erebus.to_device(backend, precond)
        y_dev = zeros(Ny1 * Nx1 * 4)
        ldiv!(y_dev, p_dev, b_vec)
        @test isapprox(y_vec, y_dev; rtol=1.0e-12, atol=1.0e-14)

        # Bidirectional to_host transfer
        p_host = Erebus.to_host(p_dev)
        y_host = zeros(Ny1 * Nx1 * 4)
        ldiv!(y_host, p_host, b_vec)
        @test isapprox(y_vec, y_host; rtol=1.0e-12, atol=1.0e-14)
    end

    @testset "SubArray View Execution and Backend Consistency" begin
        x_pad = randn(Ny1 * Nx1 * 4 + 20)
        y_pad = zeros(Ny1 * Nx1 * 4 + 20)
        x_view = view(x_pad, 1:(Ny1 * Nx1 * 4))
        y_view = view(y_pad, 1:(Ny1 * Nx1 * 4))

        y_ref = zeros(Ny1 * Nx1 * 4)
        mul!(y_ref, op, copy(x_view))
        mul!(y_view, op, x_view)
        @test isapprox(y_ref, y_view; rtol=1.0e-12, atol=1.0e-14)

        # Preconditioner with SubArray views
        precond = build_multigrid_preconditioner(op; levels=3)
        y_prec_ref = zeros(Ny1 * Nx1 * 4)
        ldiv!(y_prec_ref, precond, copy(x_view))
        fill!(y_view, 0.0)
        ldiv!(y_view, precond, x_view)
        @test isapprox(y_prec_ref, y_view; rtol=1.0e-12, atol=1.0e-14)

        # Automatic redirection on device operator
        op_dev = Erebus.to_device(backend, op)
        d_dev = compute_operator_diagonal(op_dev)
        d_cpu = compute_operator_diagonal(op)
        @test isapprox(d_cpu, d_dev; rtol=1.0e-12, atol=1.0e-14)

        # Build hierarchy directly from device operator
        h_from_dev = build_staggered_multigrid_hierarchy(op_dev; max_levels=3)
        @test length(h_from_dev.levels) == 3

        # Non-identity to_host reconstruction with SubArray fields
        sub_view = view(fill(1.0, 40, 40), 1:Ny1, 1:Nx1)
        S = typeof(sub_view)
        op_sub = Erebus.MatrixFreeStokesDarcyOperator{Float64,S}(
            op.Ny1,
            op.Nx1,
            op.dx,
            op.dy,
            op.Nx_val,
            op.Ny_val,
            sub_view,
            sub_view,
            sub_view,
            sub_view,
            sub_view,
            sub_view,
            sub_view,
            sub_view,
            sub_view,
            sub_view,
            sub_view,
            sub_view,
            sub_view,
            sub_view,
            sub_view,
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
        )
        op_h = Erebus.to_host(op_sub)
        @test op_h.ETA isa Matrix{Float64}
        @test size(op_h.ETA) == (Ny1, Nx1)
        @test op_h.Ny1 == op.Ny1
        @test op_h.Nx1 == op.Nx1
    end

    @testset "Package Extension Syntax and Guards" begin
        ext_dir = joinpath(@__DIR__, "..", "ext")
        for ext_file in ["ErebusCUDAExt.jl", "ErebusAMDGPUExt.jl", "ErebusMetalExt.jl"]
            path = joinpath(ext_dir, ext_file)
            @test isfile(path)
            content = read(path, String)
            expr = Meta.parse("begin\n" * content * "\nend")
            @test expr isa Expr
        end
    end
end

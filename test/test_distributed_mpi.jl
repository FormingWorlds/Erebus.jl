# Unit and integration tests for distributed memory domain decomposition via MPI.jl

using Test
using LinearAlgebra
using Random
using MPI
using Erebus

if !MPI.Initialized()
    MPI.Init()
end

@testset "Distributed Memory & MPI Domain Decomposition" begin
    @testset "Serial Topology Fallbacks & Partitioning" begin
        topo = SerialTopology()
        @test !is_distributed(topo)
        @test !is_distributed(42)

        # 1D index partitioning tests across procs
        r1 = partition_indices(10, 3, 1)
        r2 = partition_indices(10, 3, 2)
        r3 = partition_indices(10, 3, 3)
        @test r1 == 1:4
        @test r2 == 5:7
        @test r3 == 8:10
        @test length(r1) + length(r2) + length(r3) == 10

        # Boundary checks for partition_indices
        @test_throws ArgumentError partition_indices(10, 0, 1)
        @test_throws ArgumentError partition_indices(10, 3, 4)

        # Fallback operations on SerialTopology
        A = fill(1.0, 4, 4)
        exchange_halos!(topo, A)
        @test A[1, 1] ≈ 1.0
        @test A[4, 4] ≈ 1.0

        v1 = [1.0, 2.0, 3.0]
        v2 = [4.0, 5.0, 6.0]
        @test distributed_dot(topo, v1, v2) ≈ 32.0
        @test distributed_norm(topo, v1) ≈ sqrt(14.0)

        mig = migrate_markers!(topo, [1.0], [2.0])
        @test mig == (0, 0)
        @test is_distributed(topo) == false
    end

    @testset "MPI Configuration & Validation" begin
        cfg_def = default_config()
        @test !cfg_def.mpi.enable
        @test cfg_def.mpi.halo_width == 1
        @test cfg_def.mpi.px == 0
        @test cfg_def.mpi.py == 0

        # Valid custom MPI config
        cfg_custom = SimulationConfig(;
            mpi=MPIConfig(; enable=true, px=2, py=2, halo_width=1)
        )
        @test cfg_custom.mpi.enable
        @test cfg_custom.mpi.px == 2
        @test cfg_custom.mpi.py == 2
        validate_config(cfg_custom)

        # Serialization to dict
        d = config_to_dict(cfg_custom)
        @test haskey(d, "mpi")
        @test d["mpi"]["enable"] == true
        @test d["mpi"]["px"] == 2
        @test d["mpi"]["py"] == 2

        # Invalid MPI configs
        @test_throws ArgumentError validate_config(
            SimulationConfig(; mpi=MPIConfig(; enable=true, px=-1))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; mpi=MPIConfig(; enable=true, py=-1))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; mpi=MPIConfig(; enable=true, halo_width=0))
        )
    end

    @testset "DistributedGridTopology2D Single Process" begin
        comm = MPI.COMM_SELF
        topo = DistributedGridTopology2D(
            comm,
            32,
            32;
            px=1,
            py=1,
            halo_width=1,
            xmin=0.0,
            xmax=100.0,
            ymin=0.0,
            ymax=100.0,
        )
        @test is_distributed(topo)
        @test topo.rank == 0
        @test topo.size == 1
        @test topo.px == 1
        @test topo.py == 1
        @test topo.coords == (0, 0)
        @test topo.neighbors.north == MPI.PROC_NULL
        @test topo.neighbors.south == MPI.PROC_NULL
        @test topo.neighbors.east == MPI.PROC_NULL
        @test topo.neighbors.west == MPI.PROC_NULL
        @test topo.dx ≈ 100.0 / 32
        @test topo.dy ≈ 100.0 / 32
        @test topo.xmin_loc ≈ 0.0
        @test topo.xmax_loc ≈ 100.0
        @test topo.ymin_loc ≈ 0.0
        @test topo.ymax_loc ≈ 100.0

        # Parameter validation checks
        @test_throws ArgumentError DistributedGridTopology2D(comm, 32, 32; px=64, py=1)
        @test_throws ArgumentError DistributedGridTopology2D(comm, 32, 32; px=1, py=64)
        @test_throws ArgumentError DistributedGridTopology2D(comm, 32, 32; px=2, py=1)

        # Automatic factorization when one coordinate is unspecified
        topo_auto_y = DistributedGridTopology2D(comm, 32, 32; px=1, py=0)
        @test topo_auto_y.px == 1
        @test topo_auto_y.py == 1

        buf = HaloBuffer{Float64}(4, topo.ny1_loc, topo.nx1_loc; halo_width=1)
        @test size(buf.send_east) == (4, topo.ny1_loc, 1)
        @test size(buf.send_north) == (4, 1, topo.nx1_loc)
        @test size(buf.work_x) == (4, topo.ny1_loc, topo.nx1_loc)
    end

    @testset "Distributed Dot and Norm on Comm Self" begin
        comm = MPI.COMM_SELF
        topo = DistributedGridTopology2D(comm, 16, 16; px=1, py=1)
        x = [1.0, 2.0, 3.0, 4.0]
        y = [2.0, 0.5, 1.0, -1.0]

        d_val = distributed_dot(topo, x, y)
        n_val = distributed_norm(topo, x)
        @test d_val ≈ dot(x, y)
        @test n_val ≈ norm(x)

        dist_x = DistributedVector(x, topo, 4)
        dist_y = DistributedVector(y, topo, 4)
        @test is_distributed(dist_x)
        @test length(dist_x) == 4
        @test dist_x[2] ≈ 2.0
        dist_x[2] = 5.0
        @test dist_x[2] ≈ 5.0
        dist_x[2] = 2.0
        @test dot(dist_x, dist_y) ≈ dot(x, y)
        @test norm(dist_x) ≈ norm(x)

        dist_sim = similar(dist_x)
        @test length(dist_sim) == 4
        dist_copy = copy(dist_x)
        @test dist_copy[1] ≈ dist_x[1]
        @test dist_copy[4] ≈ dist_x[4]

        # 2D Matrix reduction with halo exclusion
        M1 = fill(2.0, 6, 6)
        M2 = fill(3.0, 6, 6)
        @test distributed_dot(topo, M1, M2) ≈ 6.0 * 16
        @test distributed_norm(topo, M1) ≈ sqrt(4.0 * 16)
    end

    @testset "Distributed Operator Application Single Process" begin
        comm = MPI.COMM_SELF
        Nx = 16
        Ny = 16
        Nx1 = Nx + 1
        Ny1 = Ny + 1
        coords = GridCoordinates(Nx, Ny; xsize=100_000.0, ysize=100_000.0)

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

        local_op = MatrixFreeStokesDarcyOperator(
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
        topo = DistributedGridTopology2D(comm, Ny, Nx; px=1, py=1)
        dist_op = DistributedStokesDarcyOperator(local_op, topo)

        @test is_distributed(dist_op)
        @test size(dist_op) == size(local_op)
        @test eltype(dist_op) == Float64

        rng = MersenneTwister(42)
        x = randn(rng, Ny1 * Nx1 * 4)
        x_orig = copy(x)
        y_ref = zeros(Ny1 * Nx1 * 4)
        y_dist = zeros(Ny1 * Nx1 * 4)

        mul!(y_ref, local_op, x)
        mul!(y_dist, dist_op, x)
        @test y_dist ≈ y_ref
        @test norm(y_dist) ≈ norm(y_ref)

        # Assert mul! does NOT mutate input vector x
        @test x == x_orig

        # 5-argument mul!
        y_dist5 = zeros(Ny1 * Nx1 * 4)
        mul!(y_dist5, dist_op, x, 2.0, 0.0)
        @test y_dist5 ≈ 2.0 .* y_ref
        @test norm(y_dist5) ≈ 2.0 * norm(y_ref)
        @test x == x_orig
    end

    @testset "Marker Particle Migration Single Process" begin
        comm = MPI.COMM_SELF
        topo = DistributedGridTopology2D(
            comm, 16, 16; px=1, py=1, xmin=0.0, xmax=100.0, ymin=0.0, ymax=100.0
        )
        xm = [10.0, 25.0, 50.0, 90.0]
        ym = [15.0, 30.0, 60.0, 85.0]
        tm = [1, 2, 1, 3]
        tkm = [280.0, 290.0, 300.0, 310.0]

        # In single-process mode, all markers remain in local domain
        out_cnt, in_cnt = migrate_markers!(topo, xm, ym, tm, tkm)
        @test out_cnt == 0
        @test in_cnt == 0
        @test length(xm) == 4
        @test length(tm) == 4
        @test xm[1] ≈ 10.0
        @test tm[4] == 3

        # NamedTuple interface
        markers_nt = (xm=xm, ym=ym, tm=tm, tkm=tkm)
        out_nt, in_nt = migrate_markers!(topo, markers_nt)
        @test out_nt == 0
        @test in_nt == 0
        @test length(markers_nt.xm) == 4
        @test markers_nt.tkm[2] ≈ 290.0
    end

    function run_mpi_test_with_timeout(cmd::Cmd; timeout_secs::Real=90.0)
        proc = run(cmd, wait=false)
        timed_out = false
        timer = Timer(timeout_secs) do _
            timed_out = true
            kill(proc)
        end
        wait(proc)
        close(timer)
        !timed_out || error("MPI command timed out after $(timeout_secs)s: $cmd")
        return proc
    end

    @testset "Multi-Process MPI Integration (2 Ranks)" begin
        script = joinpath(@__DIR__, "mpi_worker_tests.jl")
        cmd = `$(MPI.mpiexec()) -n 2 $(Base.julia_cmd()) --project=$(normpath(joinpath(@__DIR__, ".."))) $script 2`
        p = run_mpi_test_with_timeout(cmd; timeout_secs=90.0)
        @test success(p)
        @test p.exitcode == 0
    end

    @testset "Multi-Process MPI Integration (4 Ranks 2x2 Grid)" begin
        script = joinpath(@__DIR__, "mpi_worker_tests.jl")
        cmd = `$(MPI.mpiexec()) -n 4 $(Base.julia_cmd()) --project=$(normpath(joinpath(@__DIR__, ".."))) $script 4`
        p = run_mpi_test_with_timeout(cmd; timeout_secs=90.0)
        @test success(p)
        @test p.exitcode == 0
    end

    @testset "Simulation Loop MPI Guard" begin
        cfg_mpi = SimulationConfig(; mpi=MPIConfig(; enable=true))
        @test cfg_mpi.mpi.enable
        @test_throws ErrorException simulation_loop(cfg_mpi)
    end
end

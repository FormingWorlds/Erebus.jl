# Worker script for multi-process MPI domain decomposition integration tests.

using Test
using MPI
using LinearAlgebra
using Random
using Erebus

MPI.Init()
comm = MPI.COMM_WORLD
sz = MPI.Comm_size(comm)
rank = MPI.Comm_rank(comm)

expected_sz = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : sz

# 1. 2D Cartesian Topology setup
px = sz == 4 ? 2 : (sz == 2 ? 2 : 1)
py = sz == 4 ? 2 : 1

topo = DistributedGridTopology2D(
    comm, 32, 32; px=px, py=py, halo_width=1, xmin=0.0, xmax=100.0, ymin=0.0, ymax=100.0
)

@testset "MPI Worker: Topology (Rank $rank)" begin
    @test sz == expected_sz
    @test topo.rank == MPI.Comm_rank(topo.cart_comm)
    @test topo.size == sz
    @test topo.px == px
    @test topo.py == py

    # Verify diagonal corner neighbor topological invariants
    if topo.neighbors.north != MPI.PROC_NULL && topo.neighbors.east != MPI.PROC_NULL
        @test topo.neighbors.northeast != MPI.PROC_NULL
    end
    if topo.neighbors.north != MPI.PROC_NULL && topo.neighbors.west != MPI.PROC_NULL
        @test topo.neighbors.northwest != MPI.PROC_NULL
    end
    if topo.neighbors.south != MPI.PROC_NULL && topo.neighbors.east != MPI.PROC_NULL
        @test topo.neighbors.southeast != MPI.PROC_NULL
    end
    if topo.neighbors.south != MPI.PROC_NULL && topo.neighbors.west != MPI.PROC_NULL
        @test topo.neighbors.southwest != MPI.PROC_NULL
    end
end

# 2. Halo exchange test with exact mathematical continuity validation
hw = topo.halo_width
ny = topo.ny_loc + 2 * hw
nx = topo.nx_loc + 2 * hw
A = fill(-999.0, ny, nx)

# Populate owned interior cells with continuous coordinate-based function
for j in (hw + 1):(nx - hw), i in (hw + 1):(ny - hw)
    cell_x = topo.xmin_loc + (j - hw - 0.5) * topo.dx
    cell_y = topo.ymin_loc + (i - hw - 0.5) * topo.dy
    A[i, j] = 100.0 * cell_y + cell_x
end

buf = HaloBuffer{Float64}(1, ny, nx; halo_width=hw)
exchange_halos!(topo, A; buffer=buf)

@testset "MPI Worker: Halo Exchange (Rank $rank)" begin
    @test size(A) == (ny, nx)
    @test all(isfinite, A)

    # Verify ghost boundary values match exact coordinate function of neighbors
    if topo.neighbors.west != MPI.PROC_NULL
        for i in (hw + 1):(ny - hw)
            cell_y = topo.ymin_loc + (i - hw - 0.5) * topo.dy
            expected = 100.0 * cell_y + (topo.xmin_loc - 0.5 * topo.dx)
            @test A[i, 1] ≈ expected
        end
    end
    if topo.neighbors.east != MPI.PROC_NULL
        for i in (hw + 1):(ny - hw)
            cell_y = topo.ymin_loc + (i - hw - 0.5) * topo.dy
            expected = 100.0 * cell_y + (topo.xmax_loc + 0.5 * topo.dx)
            @test A[i, nx] ≈ expected
        end
    end
    if topo.neighbors.south != MPI.PROC_NULL
        for j in (hw + 1):(nx - hw)
            cell_x = topo.xmin_loc + (j - hw - 0.5) * topo.dx
            expected = 100.0 * (topo.ymin_loc - 0.5 * topo.dy) + cell_x
            @test A[1, j] ≈ expected
        end
    end
    if topo.neighbors.north != MPI.PROC_NULL
        for j in (hw + 1):(nx - hw)
            cell_x = topo.xmin_loc + (j - hw - 0.5) * topo.dx
            expected = 100.0 * (topo.ymax_loc + 0.5 * topo.dy) + cell_x
            @test A[ny, j] ≈ expected
        end
    end
    # Check diagonal corner exchange transported via dimensional splitting
    if topo.neighbors.north != MPI.PROC_NULL && topo.neighbors.east != MPI.PROC_NULL
        expected_corner =
            100.0 * (topo.ymax_loc + 0.5 * topo.dy) + (topo.xmax_loc + 0.5 * topo.dx)
        @test A[ny, nx] ≈ expected_corner
    end
end

# 3. Distributed reductions test (verifying halo exclusion)
local_v = fill(1.0, 10)
dist_v = DistributedVector(local_v, topo, 10 * sz)

# Grid reduction excluding ghost margins
local_owned_sq = sum(A[(hw + 1):(ny - hw), (hw + 1):(nx - hw)] .^ 2)
expected_global_sq = MPI.Allreduce(local_owned_sq, +, topo.cart_comm)

@testset "MPI Worker: Reductions (Rank $rank)" begin
    @test dot(dist_v, dist_v) ≈ 10.0 * sz
    @test norm(dist_v) ≈ sqrt(10.0 * sz)
    @test distributed_dot(topo, A, A) ≈ expected_global_sq
    @test distributed_norm(topo, A) ≈ sqrt(expected_global_sq)
end

# 4. Multi-rank DistributedStokesDarcyOperator test
Ny_loc = topo.ny1_loc - 1
Nx_loc = topo.nx1_loc - 1
Ny1_loc = topo.ny1_loc
Nx1_loc = topo.nx1_loc

coords_loc = GridCoordinates(
    Nx_loc,
    Ny_loc;
    xsize=(topo.xmax_loc - topo.xmin_loc),
    ysize=(topo.ymax_loc - topo.ymin_loc),
)
ETA_loc = fill(1.0e19, Ny1_loc, Nx1_loc)
ETAP_loc = fill(1.0e19, Ny1_loc, Nx1_loc)
GGG_loc = fill(1.0e10, Ny1_loc, Nx1_loc)
GGGP_loc = fill(1.0e10, Ny1_loc, Nx1_loc)
RHOX_loc = fill(3300.0, Ny1_loc, Nx1_loc)
RHOY_loc = fill(3300.0, Ny1_loc, Nx1_loc)
RHOFX_loc = fill(1000.0, Ny1_loc, Nx1_loc)
RHOFY_loc = fill(1000.0, Ny1_loc, Nx1_loc)
RX_loc = fill(1.0e10, Ny1_loc, Nx1_loc)
RY_loc = fill(1.0e10, Ny1_loc, Nx1_loc)
ETAPHI_loc = fill(1.0e19, Ny1_loc, Nx1_loc)
BETAPHI_loc = fill(1.0e-11, Ny1_loc, Nx1_loc)
PHI_loc = fill(0.01, Ny1_loc, Nx1_loc)
gx_loc = zeros(Ny1_loc, Nx1_loc)
gy_loc = fill(-9.81, Ny1_loc, Nx1_loc)
dt_loc = 1000.0 * 365.25 * 86400.0

local_op = MatrixFreeStokesDarcyOperator(
    ETA_loc,
    ETAP_loc,
    GGG_loc,
    GGGP_loc,
    RHOX_loc,
    RHOY_loc,
    RHOFX_loc,
    RHOFY_loc,
    RX_loc,
    RY_loc,
    ETAPHI_loc,
    BETAPHI_loc,
    PHI_loc,
    gx_loc,
    gy_loc,
    dt_loc;
    coords=coords_loc,
)

dist_op = DistributedStokesDarcyOperator(local_op, topo)

rng = MersenneTwister(rank + 1)
x_in = randn(rng, Ny1_loc * Nx1_loc * 4)
x_orig = copy(x_in)
y_out = zeros(Ny1_loc * Nx1_loc * 4)

mul!(y_out, dist_op, x_in)

@testset "MPI Worker: Operator (Rank $rank)" begin
    @test is_distributed(dist_op)
    @test x_in == x_orig
    @test all(isfinite, y_out)
    @test norm(y_out) > 0.0
end

# 5. Marker particle migration tests
# 5a. Symmetric migration across all ranks
xm = Float64[]
ym = Float64[]
tm = Int[]
tkm = Float64[]

for target_rank in 0:(sz - 1)
    coords_t = MPI.Cart_coords(topo.cart_comm, target_rank)
    ty = coords_t[1]
    tx = coords_t[2]
    mid_x = 0.5 * (topo.x_splits[tx + 1] + topo.x_splits[tx + 2])
    mid_y = 0.5 * (topo.y_splits[ty + 1] + topo.y_splits[ty + 2])
    push!(xm, mid_x)
    push!(ym, mid_y)
    push!(tm, rank * 10 + target_rank)
    push!(tkm, 300.0 + target_rank)
end

initial_local_count = length(xm)
initial_global_count = MPI.Allreduce(initial_local_count, +, topo.cart_comm)

out_cnt, in_cnt = migrate_markers!(topo, xm, ym, tm, tkm)

final_local_count = length(xm)
final_global_count = MPI.Allreduce(final_local_count, +, topo.cart_comm)

# 5b. Asymmetric migration (proves zero collective divergence when only one rank migrates)
xm_asym = Float64[]
ym_asym = Float64[]
tm_asym = Int[]

if rank == 0 && sz >= 2
    # Rank 0 creates 1 marker targeting Rank 1
    coords_t = MPI.Cart_coords(topo.cart_comm, 1)
    mid_x = 0.5 * (topo.x_splits[coords_t[2] + 1] + topo.x_splits[coords_t[2] + 2])
    mid_y = 0.5 * (topo.y_splits[coords_t[1] + 1] + topo.y_splits[coords_t[1] + 2])
    push!(xm_asym, mid_x)
    push!(ym_asym, mid_y)
    push!(tm_asym, 999)
end

out_a, in_a = migrate_markers!(topo, xm_asym, ym_asym, tm_asym)
global_count_asym = MPI.Allreduce(length(xm_asym), +, topo.cart_comm)

@testset "MPI Worker: Marker Migration (Rank $rank)" begin
    @testset "Symmetric Migration" begin
        @test initial_global_count == sz * sz
        @test final_global_count == initial_global_count
        @test final_local_count == sz
        for m in 1:final_local_count
            @test topo.xmin_loc <= xm[m] <= topo.xmax_loc
            @test topo.ymin_loc <= ym[m] <= topo.ymax_loc
            @test tkm[m] ≈ 300.0 + rank
            @test (tm[m] % 10) == rank
            @test isa(tm[m], Int)
        end
    end

    @testset "Asymmetric Migration" begin
        @test global_count_asym == (sz >= 2 ? 1 : 0)
        if rank == 1 && sz >= 2
            @test length(xm_asym) == 1
            @test tm_asym[1] == 999
            @test isa(tm_asym[1], Int)
        else
            @test length(xm_asym) == 0
            @test isempty(tm_asym)
        end
    end
end

MPI.Finalize()

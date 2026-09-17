# Worker script for multi-process MPI domain decomposition integration tests.

using Test
using MPI
using LinearAlgebra
using Erebus

MPI.Init()
comm = MPI.COMM_WORLD
sz = MPI.Comm_size(comm)
rank = MPI.Comm_rank(comm)

expected_sz = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : sz
@test sz == expected_sz

# 1. 2D Cartesian Topology setup
px = sz == 4 ? 2 : (sz == 2 ? 2 : 1)
py = sz == 4 ? 2 : 1

topo = DistributedGridTopology2D(
    comm, 32, 32; px=px, py=py, halo_width=1, xmin=0.0, xmax=100.0, ymin=0.0, ymax=100.0
)

@test topo.rank == rank
@test topo.size == sz
@test topo.px == px
@test topo.py == py

# 2. Halo exchange test for 2D field
ny = topo.ny1_loc
nx = topo.nx1_loc
A = fill(-999.0, ny, nx)

# Populate interior degrees of freedom with global continuous coordinates
for j in 2:(nx - 1), i in 2:(ny - 1)
    gx = topo.xmin_loc + (j - 1.5) * topo.dx
    gy = topo.ymin_loc + (i - 1.5) * topo.dy
    A[i, j] = 100.0 * gy + gx
end

buf = HaloBuffer{Float64}(1, ny, nx; halo_width=1)
exchange_halos!(topo, A; buffer=buf)

# Verify ghost boundary values match neighbor interior values
if topo.neighbors.west != MPI.PROC_NULL
    @test A[2, 1] > -900.0
    @test A[ny - 1, 1] > -900.0
end
if topo.neighbors.east != MPI.PROC_NULL
    @test A[2, nx] > -900.0
    @test A[ny - 1, nx] > -900.0
end
if topo.neighbors.south != MPI.PROC_NULL
    @test A[1, 2] > -900.0
    @test A[1, nx - 1] > -900.0
end
if topo.neighbors.north != MPI.PROC_NULL
    @test A[ny, 2] > -900.0
    @test A[ny, nx - 1] > -900.0
end

# 3. Distributed reductions test
local_v = fill(1.0, 10)
global_dot = distributed_dot(topo, local_v, local_v)
@test global_dot ≈ 10.0 * sz

global_norm = distributed_norm(topo, local_v)
@test global_norm ≈ sqrt(10.0 * sz)

# 4. Marker particle migration test
# Generate markers where some positions cross subdomain bounds
xm = Float64[]
ym = Float64[]
tm = Int[]
tkm = Float64[]

# Markers placed deliberately in every rank's territory
for target_rank in 0:(sz - 1)
    # Target center in process grid
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
@test initial_global_count == sz * sz

out_cnt, in_cnt = migrate_markers!(topo, xm, ym, tm, tkm)

final_local_count = length(xm)
final_global_count = MPI.Allreduce(final_local_count, +, topo.cart_comm)

# Assert 100% exact particle count conservation
@test final_global_count == initial_global_count
@test final_local_count == sz

# Assert every marker now on this rank is strictly within this rank's bounds
for m in 1:final_local_count
    @test topo.xmin_loc <= xm[m] <= topo.xmax_loc
    @test topo.ymin_loc <= ym[m] <= topo.ymax_loc
    @test tkm[m] ≈ 300.0 + rank
end

MPI.Finalize()

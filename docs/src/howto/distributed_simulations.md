# Distributed Computations via MPI

This guide explains how to configure and use distributed-memory primitives in `Erebus.jl` through the `MPI.jl` extension.

---

## Overview

Planetary thermal models at extreme spatial resolutions exceed the memory limit of a single compute node. To enable multi-node scaling, `Erebus.jl` provides 2D Cartesian domain decomposition and distributed memory primitives.

Key distributed components include:
- `DistributedGridTopology2D`: Maps subdomains to an MPI process grid via `MPI.Cart_create`.
- `HaloBuffer`: Stores pre-allocated buffers for non-blocking ghost cell exchange via `MPI.Isend` and `MPI.Irecv!`.
- `exchange_halos!`: Transports boundary data between adjacent subdomains in two dimensional steps.
- `DistributedStokesDarcyOperator`: Wraps local matrix-free operators and evaluates stencils at partition seams.
- `DistributedVector`: Stores partitioned degree-of-freedom vectors with global reductions (`distributed_dot`, `distributed_norm`).
- `migrate_markers!`: Transports marker particles that cross subdomain boundaries using `MPI.Alltoallv!`.

---

## Configuration Settings

The `[mpi]` configuration section controls parallel execution settings:

```toml
[grid]
xsize = 140000.0
ysize = 140000.0
Nx = 1024
Ny = 1024

[solver]
hydromech_solver = "matrix_free"
darcy_elimination = true
preconditioner = "block_schur"
krylov_method = "fgmres"
krylov_rtol = 1.0e-6

[mpi]
enable = true
px = 2
py = 2
halo_width = 1
```

Parameters:
- `enable`: Activates distributed memory settings (default: `false`).
- `px`: Number of subdomain partitions in the horizontal direction. Set to `0` for automatic factorization.
- `py`: Number of subdomain partitions in the vertical direction. Set to `0` for automatic factorization.
- `halo_width`: Width of the ghost cell layer (default: `1`).

---

## Using Distributed Primitives

Distributed functionality is provided via the `ErebusMPIExt` package extension, loaded automatically when `using MPI` is present in your environment.

### Setting Up a Distributed Topology

```julia
using MPI
using Erebus

MPI.Init()
comm = MPI.COMM_WORLD

# Create a 2D Cartesian topology over a 1024x1024 grid
topo = DistributedGridTopology2D(comm, 1024, 1024; px=2, py=2, halo_width=1)

println("Rank $(topo.rank): Local grid size $(topo.ny1_loc) x $(topo.nx1_loc)")
```

### Performing Halo Exchanges

Allocate a `HaloBuffer` and update ghost margins without memory allocations:

```julia
# Allocate field with halo margins
A = zeros(topo.ny1_loc + 2 * topo.halo_width, topo.nx1_loc + 2 * topo.halo_width)

# Pre-allocate communication buffers
buf = HaloBuffer{Float64}(1, size(A, 1), size(A, 2); halo_width=topo.halo_width)

# Exchange ghost cell values with cardinal and diagonal neighbors
exchange_halos!(topo, A; buffer=buf)
```

### Distributed Operator Evaluation

Wrap a local `MatrixFreeStokesDarcyOperator` for parallel Krylov iterations:

```julia
# Wrap local operator with topology and communication buffer
dist_op = DistributedStokesDarcyOperator(local_op, topo)

# Allocate distributed solution and source vectors
x = DistributedVector(local_x, topo, global_length)
y = similar(x)

# Evaluate matrix-vector product without mutating x
mul!(y, dist_op, x)

# Compute global inner product and norm
dot_val = distributed_dot(topo, x, y)
norm_val = distributed_norm(topo, x)
```

### Migrating Marker Particles

When marker particles move across subdomain bounds, route them to their owning MPI rank:

```julia
# Markers coordinates and associated properties
xm = [10.0, 50000.0, 120000.0]
ym = [15.0, 60000.0, 130000.0]
material_type = [1, 2, 1]

# Migrate particles across ranks
migrated_out, migrated_in = migrate_markers!(topo, xm, ym, material_type)
```

---

## Current Status and Integration Roadmap

The distributed memory primitives in `ErebusMPIExt` implement domain decomposition, halo communication, distributed operator evaluation, and marker migration. Full multi-node orchestration of the top-level time integration loop (`simulation_loop`) for cluster nodes is under active development.

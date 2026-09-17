# Distributed Simulations via MPI

This guide explains how to set up, configure, and run distributed-memory simulations on multiple compute nodes with `MPI.jl`.

---

## Overview

Planetary thermal models at extreme spatial resolutions, such as $2048 \times 2048$ nodes or larger, exceed the memory limit of a single server. To resolve fine structures, `Erebus.jl` partitions the domain into subdomains and distributes them among MPI ranks.

Key distributed components include:
- A 2D Cartesian process grid, organized through `DistributedGridTopology2D`, maps subdomains to MPI ranks.
- Pre-allocated buffers, managed in `HaloBuffer`, handle non-blocking ghost cell updates between adjacent ranks via `exchange_halos!`.
- Matrix-free operator evaluation, wrapped in `DistributedStokesDarcyOperator`, evaluates local stencils in parallel with message passing.
- Global dot products and norms, provided by `DistributedVector`, enable Krylov solvers to iterate on distributed vectors.
- Particle migration, handled by `migrate_markers!`, routes moving particles to destination ranks, maintaining exact marker counts.

---

## Configuration Settings

To run on multiple ranks, enable the `[mpi]` section in your `.toml` configuration file:

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

Here, `px` and `py` define the number of process subdivisions in the horizontal and vertical directions, respectively. If set to `0`, the solver factors the total MPI process count automatically. The parameter `halo_width` specifies the ghost cell margin, which defaults to 1.

---

## Running Multi-Process Jobs

### Local Multi-Process Testing

To test distributed simulations on a local workstation, invoke `mpiexec`:

```bash
mpiexec -n 4 julia --project=. -e '
using MPI
using Erebus

MPI.Init()
cfg = load_config("configs/distributed_example.toml")
run_simulation(cfg)
MPI.Finalize()
'
```

### HPC Cluster Deployment (Slurm)

On clusters managed by Slurm, such as Habrok, submit a batch script that requests multiple nodes and allocates compute cores:

```bash
#!/bin/bash
#SBATCH --job-name=erebus-distributed
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --partition=regular
#SBATCH --mem=64G

module load OpenMPI/4.1.4-GCC-11.3.0
module load Julia/1.10.0

srun julia --project=. -e '
using MPI
using Erebus

MPI.Init()
cfg = load_config("configs/cluster_2048.toml")
run_simulation(cfg)
MPI.Finalize()
'
```

---

## Practical Guidelines

When setting up a distributed calculation, choose `px` and `py` such that each subdomain remains close to a square, which minimizes boundary surface area relative to subdomain volume. Furthermore, maintain at least 16 markers per cell in every subdomain, which ensures smooth interpolation when markers move between neighbor ranks.

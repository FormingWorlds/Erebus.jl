# Installation Guide

This guide provides instructions for installing and setting up `Erebus.jl` on your local workstation or cluster.

---

## Prerequisites

- **Julia**: Version 1.10 or higher. You can download and install Julia via [juliaup](https://github.com/JuliaLang/juliaup):
  ```bash
  curl -fsSL https://install.julialang.org | sh
  ```
- **Git**: For version control and cloning the repository.
- **Linear Algebra Acceleration**:
  - macOS: Uses Apple Accelerate BLAS automatically.
  - Linux: Uses OpenBLAS or Intel MKL.

---

## Installing Erebus.jl

Clone the repository from GitHub:

```bash
git clone https://github.com/FormingWorlds/Erebus.jl.git
cd Erebus.jl
```

Instantiate the package dependencies using Julia's package manager:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

This installs all required direct and transitive dependencies:
- `StaticArrays.jl`: Compile-time sized static arrays for computational stencils.
- `LinearSolve.jl`: High-performance sparse linear solvers (UMFPACK, Pardiso).
- `ExtendableSparse.jl`: Fast dynamic sparse matrix assembly.
- `TOML.jl`: Standard library TOML configuration parser.
- `JLD2.jl`: HDF5-compatible binary checkpoint storage.
- `ArgParse.jl`: Command-line argument parsing.

---

## Running the Test Suite

To verify that your installation functions properly and all numerical routines pass baseline assertions:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

All unit tests (geometry, physics, particles, numerics, simulation, config, and integration) should pass cleanly.

---

## Environment and Multi-Threading

Erebus uses multi-threading for particle operations and staggered grid interpolation. For optimal performance, start Julia with multiple execution threads:

```bash
export JULIA_NUM_THREADS=auto
julia --project=.
```

Alternatively, pass the `-t` flag directly:

```bash
julia --project=. -t 8
```

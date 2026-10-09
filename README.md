# Erebus.jl

[![Build Status](https://github.com/FormingWorlds/Erebus.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/FormingWorlds/Erebus.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage Status](https://github.com/FormingWorlds/Erebus.jl/actions/workflows/Coverage.yml/badge.svg?branch=main)](https://github.com/FormingWorlds/Erebus.jl/actions/workflows/Coverage.yml?query=branch%3Amain)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://formingworlds.github.io/Erebus.jl/dev)
[![codecov](https://codecov.io/gh/FormingWorlds/Erebus.jl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/FormingWorlds/Erebus.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

`Erebus.jl` models the thermal and physical history of porous planetesimals in the early Solar System. It solves rock and pore fluid flow on a staggered grid with Lagrangian markers.

---

## Table of Contents

- [Capabilities](#capabilities)
- [Physical and Numerical Limits](#physical-and-numerical-limits)
- [Experimental Modules](#experimental-modules)
- [Setup](#setup)
- [Quickstart](#quickstart)
  - [Run from Command Line](#run-from-command-line)
  - [Run from Julia](#run-from-julia)
  - [Resume from Checkpoint](#resume-from-checkpoint)
  - [Run Parameter Sweeps](#run-parameter-sweeps)
  - [Benchmark Suite](#benchmark-suite)
- [Execution Pipeline](#execution-pipeline)
- [Tests and Quality Checks](#tests-and-quality-checks)
- [Contributing and Code Style](#contributing-and-code-style)
- [Documentation](#documentation)
- [License](#license)

---

## Capabilities

- **Marker-in-Cell Grid**: Solves 2D Cartesian flow for planetesimals in sticky air. The out-of-plane length $L(r) = 2r$ acts as a spherical proxy for mass, heat, and volatile budgets.
- **Two-Phase Stokes-Darcy Flow**: Solves viscous and plastic rock matrix flow coupled to Darcy fluid flow with Biot poroelastic compressibility and Drucker-Prager failure.
- **Implicit Heat Transport**: Solves heat conduction, fluid flow, radiogenic heat ($^{26}\text{Al}$, $^{60}\text{Fe}$), reaction and melt latent heats, shear heating, and turbulent heat flow.
- **Serpentine Reactions and Venting**: Tracks rock hydration and dehydration, surface vents, and dynamic hydrofracture permeability when pore fluid pressure exceeds rock strength.
- **Silicate Melt and Core Formation**: Models rock melting, melt ascent, and liquid metal settling to a central core via Darcy flow and Stokes rain.
- **Volatile Chemistry**: Conserves elemental budgets (H, C, N, S, O) with ice mixtures, organics, melt solubility, equilibrium gas release, metal-silicate partition, and Iron-Wüstite redox buffers.
- **Atmosphere and Escape**: Couples a 1D Guillot radiative profile with gas chemistry, photo-evaporative loss, and Jeans escape.
- **Accretion and Grid Growth**: Models body growth from pebble and embryo impacts, with domain doubling at constant cell size.
- **Self-Gravity**: Evaluates gravity fields from a 2D Poisson solve or a 3D enclosed-mass radial profile.
- **Configuration Files and Sweeps**: Reads TOML files, writes restart files, streams run metrics, and executes Latin Hypercube parameter sweeps.

---

## Physical and Numerical Limits

- **2D Geometry Proxy**: Solves momentum and mass conservation in 2D Cartesian coordinates. The out-of-plane length $L(r) = 2r$ enters volume integrals, not spatial derivatives. Flow patterns show 2D slab flow rather than 3D spherical flow.
- **Gravity Field**: The default 2D Poisson solver yields the field of an infinite cylinder. For differentiated bodies, the enclosed-mass profile gives a closer match to spherical gravity.
- **Sticky-Air Boundary**: Exterior sticky air has high viscosity ($10^{16}\text{ Pa s}$) compared to silicate melt ($10^{12}\text{ Pa s}$ floor). This forms a rigid lid rather than a free surface.
- **Linear Solvers**: The direct sparse solver (UMFPACK) is the production path for grids up to $256^2$. Iterative solvers, multigrid, GPU kernels, and MPI paths remain experimental.
- **Melt Transport**: Silicate melt moves by a kinematic drift step rather than as a separate Darcy fluid phase in the matrix momentum equations.
- **Pore Fluid Regime**: The code tracks one aqueous fluid phase. Compaction without venting loses pore fluid mass.
- **Benchmark Provenance**: Documentation figures state their provenance class. Many figures show Python reference models or isolated module outputs rather than full 2D simulation runs. Automated unit and reference tests in `test/` check individual solver components.

---

## Experimental Modules

The MPI extension (`ErebusMPIExt`) and GPU kernels (`ErebusCUDAExt`, `ErebusMetalExt`, `ErebusAMDGPUExt`) are experimental. They run in dedicated tests and do not enter production runs.

---

## Setup

Install Julia 1.12 or later. Add `Erebus.jl` through the package tool:

```julia
using Pkg
Pkg.add(url="https://github.com/FormingWorlds/Erebus.jl.git")
```

For local work:

```bash
git clone https://github.com/FormingWorlds/Erebus.jl.git
cd Erebus.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

---

## Quickstart

### Run from Command Line

Run runs with `launch.jl` and a TOML file:

```bash
# Run the 2D hydrothermal benchmark (32x32 cells, 100 km diameter body)
julia --project launch.jl configs/hydrothermal_benchmark.toml -o output_hydrothermal/

# Run a quick test run (5 timesteps)
julia --project launch.jl configs/test_quick.toml -o output_test/
```

### Run from Julia

Load a config, check parameters, and launch the run:

```julia
using Erebus

# Load and validate configuration
cfg = load_config("configs/test_quick.toml")

# Execute simulation loop
state = simulation_loop(cfg)
```

### Resume from Checkpoint

Resume a run from a saved JLD2 file with the `--restart` (`-r`) flag:

```bash
julia --project launch.jl configs/hydrothermal_benchmark.toml -o output_hydrothermal/ -r output_hydrothermal/output_00008.jld2
```

### Run Parameter Sweeps

Execute parallel parameter sweeps with Latin Hypercube sampling:

```bash
julia --project tools/run_ensemble.jl configs/test_ensemble_sweep.toml
```

### Benchmark Suite

Build verification figures and benchmark plots:

```bash
./benchmarks/run_all.sh
```

---

## Execution Pipeline

`Erebus.jl` couples staggered finite-difference grids with Lagrangian markers in a multi-stage loop:

![Erebus.jl Architecture and Execution Flow](docs/src/assets/erebus_architecture_flowchart.svg)

---

## Tests and Quality Checks

Run the test suite locally:

```bash
julia --project=. -t 4 test/runtests.jl
```

The test suite covers:
- Unit tests: grid coordinates, interpolation weights, constitutive laws, and TOML validation.
- Physical benchmarks: Terzaghi 1D consolidation, 2D hydrothermal convection, Stefan moving front, and radioactive decay.
- Mutation tests: discriminating physical checks on ten key solver functions.
- Reference baselines: regression checks for mass, heat, and volatile totals on standard setups.
- Restart tests: exact bitwise checks for resumed runs on all markers, grid fields, atmosphere, and random seed states.

---

## Contributing and Code Style

Contributions are welcome. `Erebus.jl` enforces the [BlueStyle](https://github.com/invenia/BlueStyle) code formatting convention:

1. Format code before committing:
   ```julia
   using JuliaFormatter
   format(".", BlueStyle())
   ```
2. Verify that all tests pass:
   ```bash
   julia --project=. test/runtests.jl
   ```
3. Open a pull request against `main`. Continuous Integration verifies the test suite (group `all`) on Julia 1.12 and 1.13 (`ubuntu-latest`), quick simulation on `macos-latest` (`macos-test`), coverage, documentation (`docs`), code style formatting (`format`), and dedicated checks for architecture ratchet (`check_architecture`), performance budget (`check_budget`), bitwise determinism (`check_determinism`), and test quality standards (`test-quality`).

---

## Documentation

Read complete tutorials, how-to guides, equations, and benchmarks online:

[https://formingworlds.github.io/Erebus.jl/dev](https://formingworlds.github.io/Erebus.jl/dev)

The documentation uses the Diataxis structure:
- **Tutorials**: Hands-on exercises for core formation and planetary growth.
- **How-To Guides**: Practical guides for config files, checkpoints, and parameter sweeps.
- **Explanations**: Physics of fluid flow, porous mechanics, heat flow, and gas chemistry.
- **Reference**: Config schema, benchmark validation records, and bibliography.

---

## License

`Erebus.jl` is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

# Erebus.jl

[![Build Status](https://github.com/FormingWorlds/Erebus.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/FormingWorlds/Erebus.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage Status](https://github.com/FormingWorlds/Erebus.jl/actions/workflows/Coverage.yml/badge.svg?branch=main)](https://github.com/FormingWorlds/Erebus.jl/actions/workflows/Coverage.yml?query=branch%3Amain)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://formingworlds.github.io/Erebus.jl/dev)
[![codecov](https://codecov.io/gh/FormingWorlds/Erebus.jl/branch/main/graph/badge.svg)](https://app.codecov.io/gh/FormingWorlds/Erebus.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

`Erebus.jl` is a two-dimensional hydro-thermo-mechanical-chemical (HTMC) geodynamic modeling code for planetesimal evolution. It simulates coupled two-phase fluid-solid flow, visco-elasto-plastic rock deformation, radiogenic heating, fluid thermal buoyancy, temperature-dependent Darcy percolation, and dynamic hydrofracturing in early solar system planetesimals.

---

## Table of Contents

- [Key Capabilities](#key-capabilities)
- [Installation](#installation)
- [Quickstart](#quickstart)
  - [Run from Configuration](#run-from-configuration)
  - [Resume from Checkpoint](#resume-from-checkpoint)
  - [Post-Processing and Plotting](#post-processing-and-plotting)
- [Code Architecture](#code-architecture)
- [Verification and Testing](#verification-and-testing)
- [Contributing and Code Style](#contributing-and-code-style)
- [References and Citations](#references-and-citations)
- [License](#license)

---

## Key Capabilities

- Solves coupled two-phase Stokes-Darcy fluid-solid equations on a staggered Eulerian grid.
- Incorporates poroelastic compressibility from Biot theory for solid and fluid phases.
- Models visco-elasto-plastic matrix rheology with Drucker-Prager yielding.
- Couples fluid thermal expansion and temperature-dependent Arrhenius viscosity to Darcy flux.
- Enhances Darcy permeability via dynamic hydrofracturing under fluid overpressure.
- Computes radiogenic decay heat from 26Al and 60Fe and silicate reaction kinetics.
- Tracks physical properties with Marker-in-Cell advection (Runge-Kutta 4th-order).
- Configures geometry, solvers, and physics through validated TOML files.

---

## Installation

Ensure you have Julia 1.10 or later installed. Install `Erebus.jl` via the Julia package manager:

```julia
using Pkg
Pkg.add(url="https://github.com/FormingWorlds/Erebus.jl.git")
```

For local development:

```bash
git clone https://github.com/FormingWorlds/Erebus.jl.git
cd Erebus.jl
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

---

## Quickstart

### Run from Configuration

Launch simulations directly from the command line using structured TOML input files:

```bash
# Run the 2D hydrothermal circulation benchmark (32x32 cells, 100 km diameter planetesimal)
julia --project launch.jl configs/hydrothermal_benchmark.toml -o output_hydrothermal/

# Run a quick test simulation
julia --project launch.jl configs/test_quick.toml -o output_test/
```

Available pre-configured setups in `configs/`:

- `configs/hydrothermal_benchmark.toml`: 2D hydrothermal circulation benchmark (32x32 cells, 100 km diameter planetesimal, 26Al decay heating, Darcy thermal buoyancy, Arrhenius fluid viscosity, dynamic hydrofracturing).
- `configs/default.toml`: Default reference configuration with baseline physical properties.
- `configs/test_quick.toml`: Short-duration configuration (5 timesteps) for rapid execution and sanity checks.


### Resume from Checkpoint

`Erebus.jl` saves simulation states to portable JLD2 archives. Resume an interrupted or extended run using the `--restart` (`-r`) flag:

```bash
julia --project launch.jl configs/hydrothermal_benchmark.toml -o output_hydrothermal/ -r output_hydrothermal/output_00008.jld2
```

### Post-Processing and Plotting

Helper scripts in `scripts/` produce publication-ready figures and animations from output files:

```bash
# Export benchmark simulation data to JSON
julia --project scripts/export_hydrothermal_data.jl output_hydrothermal/

# Generate the 4-panel hydrothermal circulation benchmark figure
python3 scripts/generate_hydrothermal_benchmark.py output_hydrothermal/benchmark_plot_data.json

# Generate comprehensive multi-panel diagnostic plots (requires PyPlot and StatsBase)
julia --project scripts/generate_plots.jl output_hydrothermal/

# Generate time-evolution MP4 animations (requires Plots)
julia --project scripts/generate_animations.jl output_hydrothermal/
```

---

## Code Architecture

`Erebus.jl` divides numerical responsibilities across focused submodules:

| Submodule | Responsibilities |
| :--- | :--- |
| `Erebus.Geometry` | Staggered grid coordinate vectors, basic/velocity/pressure node indexing, and domain boundary masks. |
| `Erebus.Particles` | Marker initialization, material property mapping, Marker-in-Cell interpolation, and RK4 advection. |
| `Erebus.Physics` | Constitutive equations, Arrhenius fluid viscosity, Darcy buoyancy, hydrofracture criteria, 26Al/60Fe heat production, and reaction kinetics. |
| `Erebus.Numerics` | Discretization and linear assembly of coupled Stokes-Darcy hydromechanics, Poisson gravity, and thermal energy equations. |
| `Erebus.Simulation` | Global timestep progression, non-linear Picard and plastic iterations, checkpoint serialization, and state resumption. |
| `Erebus.Config` | Strongly typed configuration structs (`SimulationConfig`), default generation, bounds validation, and TOML serialization. |

---

## Verification and Testing

Execute the automated test suite locally:

```bash
julia --project=. test/runtests.jl
```

The test suite covers:

- Unit verification: Grid metrics, particle interpolation kernels, constitutive relations, and TOML schema validators.
- Terzaghi 1D consolidation benchmark: Verifies poroelastic compressibility against the analytical Fourier series solution.
- 2D hydrothermal benchmark: Verifies coupled Darcy thermal buoyancy, fluid viscosity transitions, and hydrofracturing.
- Checkpoint restart continuity: Asserts exact numerical state matching across checkpoint save and resume boundaries.

Consult the [online documentation](https://formingworlds.github.io/Erebus.jl/dev) for tutorials and mathematical derivations.

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
3. Open a pull request against `main`. Continuous Integration verifies test execution on Julia 1.10 and 1.11, and enforces formatting and documentation builds on Julia 1.10.

---

## References and Citations

If you use `Erebus.jl` in your research, please cite:

- Hubmann, B. (2022). *Hydrology of Planetesimals*. Master's thesis, ETH Zurich. [DOI: 10.5281/zenodo.7058229](https://doi.org/10.5281/zenodo.7058229).
- Gerya, T. (2019). *Introduction to Numerical Geodynamic Modelling* (2nd ed.). Cambridge University Press. [DOI: 10.1017/9781316534243](https://doi.org/10.1017/9781316534243).

### Acknowledgements

`Erebus.jl` was originally created by Beat Hubmann (ETH Zurich). Ongoing development is supported by:

- European Research Council (ERC) Starting Grant *MagmaWorlds* (grant agreement no. 101219807).
- Dutch Research Council (NWO) NWA *PRELIFE* (grant no. NWA.1630.23.013).


---

## License

`Erebus.jl` is licensed under the Apache License 2.0. See [LICENSE](LICENSE) for details.

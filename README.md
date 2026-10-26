<h1 align="center">Erebus.jl</h1>

<p align="center">
  <a href="https://github.com/FormingWorlds/Erebus.jl/actions/workflows/CI.yml">
    <img src="https://github.com/FormingWorlds/Erebus.jl/actions/workflows/CI.yml/badge.svg?branch=main" alt="CI">
  </a>
  <a href="https://github.com/FormingWorlds/Erebus.jl/actions/workflows/nightly.yml">
    <img src="https://github.com/FormingWorlds/Erebus.jl/actions/workflows/nightly.yml/badge.svg?branch=main" alt="Nightly">
  </a>
  <a href="https://app.codecov.io/gh/FormingWorlds/Erebus.jl">
    <img src="https://codecov.io/gh/FormingWorlds/Erebus.jl/branch/main/graph/badge.svg" alt="Coverage">
  </a>
  <a href="https://formingworlds.github.io/Erebus.jl/dev/">
    <img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Documentation">
  </a>
  <a href="https://www.apache.org/licenses/LICENSE-2.0">
    <img src="https://img.shields.io/badge/License-Apache_2.0-blue.svg" alt="License">
  </a>
  <a href="https://github.com/invenia/BlueStyle">
    <img src="https://img.shields.io/badge/code%20style-blue-4495d1.svg" alt="Code Style: Blue">
  </a>
</p>

`Erebus.jl` is a Julia code for the thermal and chemical history of porous planetesimals. It solves 3 coupled processes: the deformation of solid rock (Stokes flow), the flow of pore fluid through the rock (Darcy flow), and heat transport. Radioactive decay, accretion, and segregation generate heat. It uses a staggered finite-difference grid and Lagrangian markers (particles that carry material properties). The code simulates ice melting, hydrothermal circulation, serpentinization reactions, metallic core segregation, and volatile degassing.

## Get started

Install Julia 1.12 or later. Clone the repository and instantiate dependencies:

```bash
git clone https://github.com/FormingWorlds/Erebus.jl.git
cd Erebus.jl
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

Run a 2-step test simulation from the command line:

```bash
julia --project launch.jl configs/test_quick.toml -o output_test/
```

The run writes JLD2 simulation checkpoints to `output_test`. To call `Erebus.jl` from your own Julia scripts without cloning the repository, install the package directly:

```julia
using Pkg
Pkg.add(url="https://github.com/FormingWorlds/Erebus.jl.git")
```

The [quickstart tutorial](https://formingworlds.github.io/Erebus.jl/dev/tutorials/quickstart/) shows how to read the output.

## Documentation

Full documentation: **[formingworlds.github.io/Erebus.jl/dev](https://formingworlds.github.io/Erebus.jl/dev/)**

- [Model Overview](https://formingworlds.github.io/Erebus.jl/dev/explanations/model_overview/): physical scope, governing equations, and execution stages
- [Model Limits](https://formingworlds.github.io/Erebus.jl/dev/explanations/discretization_numerics/#Physical-and-Numerical-Limits): physical and numerical limits of the 2-D solver
- [Running Tests](https://formingworlds.github.io/Erebus.jl/dev/howto/installation/#running-the-test-suite): test execution commands and validation assertions
- [Quickstart Tutorial](https://formingworlds.github.io/Erebus.jl/dev/tutorials/quickstart/): initial simulation run and output inspection
- [Installation Guide](https://formingworlds.github.io/Erebus.jl/dev/howto/installation/): package setup, environment flags, and multi-threading
- [Configuration Schema](https://formingworlds.github.io/Erebus.jl/dev/reference/config_schema/): complete parameter reference for TOML configuration files
- [Benchmark Suite](https://formingworlds.github.io/Erebus.jl/dev/reference/benchmarks/): physical verifications against analytical and reference solutions
- [Contributing Guide](https://formingworlds.github.io/Erebus.jl/dev/community/contributing/): code style, tests, and architecture checks

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

## Community and citation

Ask questions on the [discussions page](https://github.com/orgs/FormingWorlds/discussions) or report issues on [GitHub](https://github.com/FormingWorlds/Erebus.jl/issues). If you use `Erebus.jl` in your research, cite the literature listed in the [bibliography](https://formingworlds.github.io/Erebus.jl/dev/reference/bibliography/).

## License

`Erebus.jl` is licensed under the Apache License 2.0. See the [LICENSE](LICENSE) file for details.

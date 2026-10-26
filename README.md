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

**Two-dimensional thermo-hydro-mechanical evolution of early planetesimals.**

`Erebus.jl` is a Julia code for the thermal, chemical, and physical history of porous planetesimals. It couples solid matrix Stokes deformation, Darcy pore fluid flow, and radiogenic heat transport on a staggered finite-difference grid with Lagrangian markers. The code simulates ice melting, hydrothermal circulation, serpentinization reactions, metallic core segregation, and volatile degassing.

## Get started

Install `Erebus.jl` through the Julia package manager:

```julia
using Pkg
Pkg.add(url="https://github.com/FormingWorlds/Erebus.jl.git")
```

Run a test simulation with the command line interface:

```bash
julia --project launch.jl configs/test_quick.toml -o output_test/
```

See the [quickstart tutorial](https://formingworlds.github.io/Erebus.jl/dev/tutorials/quickstart/) for output inspection and next steps.

## Documentation

Full documentation: **[formingworlds.github.io/Erebus.jl/dev](https://formingworlds.github.io/Erebus.jl/dev/)**

- [Model Overview](https://formingworlds.github.io/Erebus.jl/dev/explanations/model_overview/): physical scope, governing equations, and execution stages
- [Quickstart Tutorial](https://formingworlds.github.io/Erebus.jl/dev/tutorials/quickstart/): initial simulation run and output inspection
- [Installation Guide](https://formingworlds.github.io/Erebus.jl/dev/howto/installation/): package setup, environment flags, and multi-threading
- [Configuration Schema](https://formingworlds.github.io/Erebus.jl/dev/reference/config_schema/): complete parameter reference for TOML configuration files
- [Benchmark Suite](https://formingworlds.github.io/Erebus.jl/dev/reference/benchmarks/): physical verifications against analytical and reference solutions
- [Contributing Guide](https://formingworlds.github.io/Erebus.jl/dev/community/contributing/): code style, test gates, and architecture ratchets

## Community and citation

Ask questions on the [discussions page](https://github.com/orgs/FormingWorlds/discussions) or report issues on [GitHub](https://github.com/FormingWorlds/Erebus.jl/issues). If you use `Erebus.jl` in your research, cite the literature listed in the [bibliography](https://formingworlds.github.io/Erebus.jl/dev/reference/bibliography/).

## License

`Erebus.jl` is licensed under the Apache License 2.0. See the [LICENSE](LICENSE) file for details.

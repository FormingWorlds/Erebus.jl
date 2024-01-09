# Erebus.jl

[![Build Status](https://github.com/formingworlds/Erebus.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/formingworlds/Erebus.jl/actions/workflows/CI.yml?query=branch%3Amain)
<!-- [![Stable Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://formingworlds.github.io/Erebus.jl/stable) -->
[![Dev Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://formingworlds.github.io/Erebus.jl/dev)
[![Test Coverage](https://codecov.io/gh/formingworlds/Erebus.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/formingworlds/Erebus.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

## Table of Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Introduction
`Erebus.jl` is a Julia package for simulating the geophysical evolution of planetesimals, focusing on two-phase fluid flow. It takes into account factors such as composition, mineralogical and hydrodynamic properties, temperature curves, and reaction rates[^1][^2].

## Installation
Ensure you have a current version of Julia installed on your system. You can add `Erebus.jl` using Julia's native package manager as follows:

```julia
using Pkg
Pkg.add(url="https://github.com/FormingWorlds/Erebus.jl.git")
```

## Usage

### Simulation parameters
Modify the simulation parameters in [```src/constants.jl```](src/constants.jl) as needed. For automated unit/CI testing, refer to [```src/test_constants.jl```](src/test_constants.jl). For production simulation runs, start with the parameters in [```src/prod_constants.jl```](src/prod_constants.jl).

### Complete simulation runs
Initiate a simulation run with ```julia -O3 --tauto launch.jl /PATH/TO/OUTPUT/```. The `-O3` and `--tauto` flags are used for optimization and automatic threading respectively. All intermediate simulation data is written into HDF-5 compatible ```.jld2``` files in the output directory for later analysis and plotting. The log file ```Erebus_run.log``` in the output directory contains relevant simulation progress and parameters.

### Manual operations and experiments
After installation, you can import and manually use the package in your Julia scripts as follows:

```julia
using Erebus

# Your code here
```

Refer to the linked documentation for detailed information on the internal functions and features provided by `Erebus.jl`.

### Plotting results

Use the helper script ```generate_plots.jl``` to generate plots by typing ```julia generate_plots.jl /PATH/TO/OUTPUT/```. The resulting ```.pdf``` plot files are stored in the same directory as the simulation output. Generate animations using the helper script ```generate_animations.jl``` with the command ```julia generate_animations.jl /PATH/TO/OUTPUT/```. The resulting ```.mp4``` animation files are stored in the same directory as the simulation output.

Note: These helper scripts are independent of the ```Erebus.jl``` package and are provided as-is. They have their own dependencies which need to be fulfilled manually based on your Julia installation.

## Contributing
Contributions to `Erebus.jl` are welcome. If you find a bug or have an idea for an improvement or new feature, please open an issue to discuss it. If you'd like to contribute code, please fork the repository, make your changes, and submit a pull request.

## License
`Erebus.jl` is licensed under the Apache 2.0 license. Please see the `LICENSE` file for more details.

Citations:
[^1]: For more information on the simulation code and proof-of-concept experiments, refer to the report https://doi.org/10.5281/zenodo.7058229
[^2]: For more information on the overall ansatz, refer to the book https://doi.org/10.1017/9781316534243

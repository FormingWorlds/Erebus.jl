# HydrologyPlanetesimals

[![CI](https://github.com/BeatHubmann/HydrologyPlanetesimals.jl/workflows/CI/badge.svg?branch=main)](https://github.com/BeatHubmann/HydrologyPlanetesimals.jl/actions/workflows/CI.yml)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://BeatHubmann.github.io/HydrologyPlanetesimals.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://BeatHubmann.github.io/HydrologyPlanetesimals.jl/dev)
[![Coverage](https://codecov.io/gh/BeatHubmann/HydrologyPlanetesimals.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/BeatHubmann/HydrologyPlanetesimals.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)


Julia package to simulate the hydrology in planetesimals[^1].

A MSc project at ETH Zurich under the supervision of Taras Gerya (ETHZ) and [Tim Lichtenberg](https://github.com/timlichtenberg) (Oxford).

[^1]: Based on: Gerya, T. (2019). Introduction to Numerical Geodynamic Modelling (2nd ed.). Cambridge: Cambridge University Press, [doi:10.1017/9781316534243](https://doi.org/10.1017/9781316534243).

## Usage

### Installation

At the Julia prompt, add with ```using Pkg; pkg"add https://github.com/BeatHubmann/HydrologyPlanetesimals.jl"```.

### Parameters

Adjust the simulation parameters in ```src/constants.jl``` if necessary.

### Running the simulation

Launch with ```julia -O3 --math-mode=fast launch.jl /PATH/TO/OUTPUT```.

### Plotting results

Launch with ```julia plots.jl /PATH/TO/OUTPUT```; ```.mp4``` files are stored in same directory.

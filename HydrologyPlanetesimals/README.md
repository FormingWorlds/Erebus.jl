# HydrologyPlanetesimals

[![CI](https://github.com/BeatHubmann/HydrologyPlanetesimals.jl/workflows/CI/badge.svg?branch=main)](https://github.com/BeatHubmann/HydrologyPlanetesimals.jl/actions/workflows/CI.yml)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://BeatHubmann.github.io/HydrologyPlanetesimals.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://BeatHubmann.github.io/HydrologyPlanetesimals.jl/dev)
[![Coverage](https://codecov.io/gh/BeatHubmann/HydrologyPlanetesimals.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/BeatHubmann/HydrologyPlanetesimals.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7058229.svg)](https://doi.org/10.5281/zenodo.7058229)

Julia package to simulate the hydrology in planetesimals[^1].

A MSc project at ETH Zurich under the supervision of Taras Gerya (ETH Zurich) and [Tim Lichtenberg](https://github.com/timlichtenberg) (University of Oxford).

[^1]: Based on: Gerya, T. (2019). Introduction to Numerical Geodynamic Modelling (2nd ed.). Cambridge: Cambridge University Press, [doi:10.1017/9781316534243](https://doi.org/10.1017/9781316534243).

## Usage

### Installation

At the Julia prompt, add with ```using Pkg; pkg"add https://github.com/BeatHubmann/HydrologyPlanetesimals.jl"```.

### Parameters

Adjust the simulation parameters in ```src/constants.jl``` if necessary. A set of parameters suitable for unit testing is given at ```src/test_constants.jl```.

### Running the simulation

Launch with ```julia -O3 --tauto launch.jl /PATH/TO/OUTPUT/```. All intermediate data is written into HDF-5 compatible ```.jld2``` files in the output directory for further analysis and plotting.

Relevant simulation progress and used parameters are written to the log file ```HydrologyPlanetesimals_run.log``` in the output directory.

### Plotting results

Generate plots with ```julia generate_plots.jl /PATH/TO/OUTPUT/```; ```.pdf``` plot files are stored in same directory as the simulation output.
Generate animations with ```julia generate_animations.jl /PATH/TO/OUTPUT/```; ```.mp4``` animation files are stored in same directory as the simulation output.

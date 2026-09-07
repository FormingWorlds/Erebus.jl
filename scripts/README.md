# Benchmark and Diagnostic Plotting Scripts

This directory contains standalone plotting scripts that generate benchmark, verification, and diagnostic figures for `Erebus.jl` documentation in `docs/src/assets/`.

## Master Regeneration Script

To regenerate all documentation figures in a single command, run:

```bash
./scripts/replot_all_benchmarks.sh
```

## Figure Catalog

| Documentation Asset | Script | Description |
|:---|:---|:---|
| `hcns_solubility_benchmark.png` | `generate_hcns_solubility_benchmark.py` | 4-panel multi-volatile solubility suite for H₂O, H₂, CO, CO₂, CH₄, N₂, S, graphite saturation, and SCSS ceilings. |
| `volatile_solubility_benchmark.png` | `generate_volatile_solubility_benchmark.py` | 4-panel diagnostic figure for water solubility, IW buffer oxygen fugacity, nitrogen dissolution, and organic devolatilization. |
| `hydrofracture_venting_benchmark.png` | `generate_hydrofracture_venting_benchmark.py` | 4-panel benchmark for cryogenic ice sealing, phase space regime boundaries, saw-tooth overpressure cycles, and episodic mass fluxes. |
| `cold_surface_venting_benchmark.png` | `generate_cold_venting_benchmark.py` | Multi-panel verification for cold lid hydrofracture and episodic venting dynamics. |
| `jeans_escape_benchmark.png` | `generate_jeans_escape_benchmark.py` | 4-panel diagnostic verification for thermal Jeans parameter, effusion suppression, loss timescale, and dynamic inventory partitioning across body sizes. |
| `terzaghi_benchmark.png` / `.svg` | `generate_terzaghi_benchmark.py` | Analytical 1D Fourier series consolidation benchmark compared with numerical Stokes-Darcy solution and pointwise relative discretization error. |
| `poroelastic_verification.png` / `.svg` | `generate_poroelastic_benchmark.py` | Poroelastic constitutive limits for Biot-Willis coupling ($K_{\mathrm{BW}}$) and Skempton coefficient ($B$) against grain and fluid compressibility limits. |
| `darcy_buoyancy_verification.png` / `.svg` | `generate_buoyancy_benchmark.py` | Verification of two-phase buoyancy coupling and Darcy fluid velocity scaling. |
| `fluid_viscosity_temperature.png` / `.svg` | `generate_viscosity_benchmark.py` | Temperature-dependent dynamic water viscosity formulation and physical limits. |
| `hydrofracture_verification.png` / `.svg` | `generate_hydrofracture_benchmark.py` | Non-linear tensile overpressure hydrofracture permeability transition and saturation ceiling. |
| `disk_temperature_multidistance_multimass.png` | `generate_disk_temperature_plots.py` | Protoplanetary disk surface temperature evolution across radial distances and stellar host masses. |
| `hydrothermal_circulation_benchmark.png` / `.svg` | `generate_hydrothermal_benchmark.py` | Verification of 2D coupled Stokes-Darcy hydrothermal convection cells. |
| `hydrothermal_grid_convergence.png` | `compare_grid_convergence.py` | Numerical grid convergence comparison across mesh resolutions ($16\times 16$, $32\times 32$, $64\times 64$, $128\times 128$). |
| `hydrothermal_porosity_sweep.png` | `compare_porosity_sweep.py` | Hydrothermal circulation regime transitions across initial matrix porosity values. |
| `hydrothermal_reaction_*.png` | `plot_reaction_benchmark.py` | Serpentinization reaction progress, fluid consumption, and heat generation fields. |

## Running Individual Scripts

Each Python script is self-contained and uses standard scientific libraries (`numpy`, `matplotlib`):

```bash
python3 scripts/generate_hcns_solubility_benchmark.py
python3 scripts/generate_jeans_escape_benchmark.py
python3 scripts/generate_hydrofracture_venting_benchmark.py
python3 scripts/generate_terzaghi_benchmark.py
python3 scripts/generate_poroelastic_benchmark.py
```

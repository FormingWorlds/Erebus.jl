# Erebus.jl Benchmark Suite

This directory contains standalone benchmark exporters and figure scripts for `Erebus.jl`.
The suite produces diagnostic plots and datasets shown in the documentation.

## Code Structure

The benchmark suite separates numerical data output from plotting:

1. **Julia Benchmark Exporters**:
   - `export_core_formation_benchmark.jl`: 1D core formation and metal segregation solver.
   - `export_degassing_benchmark.jl`: Magma ocean gas release against analytical Burnham solubility models.
   - `export_gravity_two_layer_benchmark.jl`: Gravity fields in differentiated core-mantle structure.
   - `export_hydrofracture_ramp_benchmark.jl`: Hydrofracture permeability response under fluid overpressure ramps.
   - `export_jeans_effusion_benchmark.jl`: Atmospheric Jeans escape and kinetic effusion fluxes.
   - `export_radiogenic_decay_benchmark.jl`: Short-lived radioactive decay (26Al, 60Fe) and half-life fits.
   - `export_thermal_slab_benchmark.jl`: 2D transient heat conduction against analytical Fourier series.
   - `generate_soft_turbulence_benchmarks.jl`: 1D planetesimal magma ocean thermal model with soft turbulence.

   Data from these scripts are written to `output_files/` in structured JSON format.

2. **Python Plotting and Diagram Scripts**:
   - Python scripts (`generate_*_benchmarks.py`, `render_*_movie.py`) load the JSON datasets or evaluate analytical reference models.
   - Procedural diagram scripts (`generate_architecture_diagram.py`, `generate_reservoir_diagram.py`) draw light, brand-free vector diagrams.
   - Plots, diagrams, and animations are saved directly to `docs/src/assets/`.

## Provenance Taxonomy

All figures from this suite follow a four-tier grouping:

- **Class A (2D Simulation Output)**: Figures from multi-dimensional `Erebus.jl` simulation runs.
- **Class B (Julia Library Exporter / 1D Benchmark Solver)**: Figures from Julia exporters that evaluate compiled `Erebus.jl` modules or 1D finite-difference benchmark solvers.
- **Class C (Analytical / Empirical Reference)**: Figures displaying closed-form mathematical formulas or published fits evaluated in Python for checks against discrete simulation behavior.
- **Class D (Schematic Diagram)**: Procedural vector graphics drawn by Python scripts (`generate_architecture_diagram.py` and `generate_reservoir_diagram.py`) showing system layouts and reservoir networks.

Captions and documentation pages state the provenance tier for every figure.

## Reproduction Steps

### Environment Setup

Julia benchmark exporters require the packages defined in `benchmarks/Project.toml`:

```bash
julia --project=benchmarks -e 'using Pkg; Pkg.instantiate()'
```

Python scripts require Python 3 with `numpy`, `matplotlib`, and `Pillow`.

### Execution

To run all Julia exporters and Python figure scripts in order:

```bash
./benchmarks/run_all.sh
```

To run an individual exporter and its paired plot script:

```bash
julia --project=benchmarks benchmarks/export_thermal_slab_benchmark.jl
python3 benchmarks/generate_thermal_slab_benchmark.py
```

## Quantitative Tolerances

| Benchmark | Reference Standard | Tolerance |
|:---|:---|:---|
| Degassing | Burnham (1979); Dixon et al. (1995) | Max relative error $< 10^{-6}$ |
| Two-Layer Gravity | Analytical shell theorem | Surface gravity relative error $< 10^{-3}$ |
| Hydrofracture Ramp | Rubin (1995); Hubmann (2022) | Pointwise formula match $< 10^{-6}$ |
| Jeans Effusion | Jeans (1925); Chamberlain (1963) | Relative error $< 10^{-6}$ |
| Radiogenic Decay | Tang & Dauphas (2012) | Half-life error $< 10^{-10}$ |
| Thermal Slab | Carslaw & Jaeger (1959) | $L_2$ relative error $< 10^{-3}$, energy drift $< 10^{-12}$ |

## Visual Style

All benchmark figures follow a clean, light neutral scientific style:
- Figure facecolor: `#FFFFFF`
- Axes facecolor: `#FFFFFF`
- Neutral palette: `#1F4E79` (Slate Blue), `#C0392B` (Crimson), `#27AE60` (Green), `#D97706` (Amber), `#2C3E50` (Charcoal)
- Tick marks and axis spines: `#222222` with linewidth 0.8 pt
- Grid lines: light dotted `#E0E0E0` with linewidth 0.6 pt

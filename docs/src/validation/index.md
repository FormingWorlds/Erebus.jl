# Validation

## What this section checks

Validation means that we run the code on a case where we know the right answer from another source, and we check that the code gives that answer. The other source can be an exact mathematical solution, a second independent solver, or a published measurement.

The pages in [Physics Function Reference](../reference/benchmarks/index.md) mostly check one function at a time. The pages here check how the parts work together, in most cases with the complete `simulation_loop`, where a fault in any step of the time loop changes the result.

Each case on a page has the same layout:

- why the check exists, which names the fault it catches;
- what it tests, with the setup and the expected result;
- how we know the right answer, which names the source of the reference;
- why it matters, which names the model results that depend on it;
- a figure, a note on how to read it, and a box with the test details.

## Pages

| Page | Part of the model | Reference |
|---|---|---|
| [Spherical Heat Conduction](thermal_conduction.md) | Heat conduction, the spherical geometry term, radiogenic heating in the time loop | Exact solution, backward-Euler series, 1-D finite-volume solver |
| [Spherical Darcy Flow](spherical_darcy_flow.md) | Darcy pressure and flux, the spherical face weights of the fluid continuity rows, a vented rate equal to the produced rate | Exact steady pressure and flux of a uniform source in a sphere and a disk, with uniform permeability and with a low-permeability core |
| [Spherical Compaction](spherical_compaction.md) | The compaction term of the solid and the fluid continuity rows, their spherical face weights, a vented rate equal to the rate that compaction squeezes out | Exact compaction-length profile of the effective pressure and the solid velocity in a sphere and a disk |
| [Spherical Pressure Diffusion](spherical_pressure_diffusion.md) | The poroelastic storage terms of the solid and the fluid continuity rows over 20 time steps, the spherical face weights of the fluid row | Backward-Euler series of pore pressure diffusion in a rigid sphere and disk with a vented surface |
| [Spherical Momentum](spherical_momentum.md) | How the flow of the solid changes in a slice of a sphere: the weights on the stress, the out-of-plane stress $S_{33}$ and its hoop term, the volumetric trace of the stress, and the stored `S33`, `EII`, `SII` and heating | Exact shear flow around the centre with and without `spherical_flow`, a manufactured flow towards and away from the centre for the momentum equation of the spherical annulus, and a radial compaction flow of a compressible rock for the 3-D equations of a sphere |
| [Porous Shell](porous_shell.md) | The dry porosity model in the time loop: cold pressing, sintering, the porous conductivity and heat, and the depth of the porous shell | Independent 1-D finite-volume model of a sphere with the same laws, at the same time step |
| [Melting Shell](melting_shell.md) | Silicate melting with latent heat in the porous shell case: the apparent heat capacity of a marker, the centre peak temperature and the radius of the molten centre | The same 1-D model with the latent heat as an apparent heat capacity, at the same time step |

## How the results are stored

Each validation page has a record file in `benchmarks/records/`. The record holds the measured values, the tolerances, a hash of the configuration, and the last commit that changed `src/`. The table on the page comes from the record, and `tools/check_doc_numbers.jl` fails if the two differ.

`Pkg.test` runs each case again and compares it with its record.
The Darcy flow, the compaction, the pressure diffusion and the momentum tests skip the finest grid and use its recorded values.
`tools/check_validation.jl` runs every case of its table on every grid and compares it with its record, and the CI job `check_validation` runs it on every pull request.
To write a new record after an intended physics change, run the export tool for that record.
The tool also writes the table on the page.
It refuses to run with uncommitted changes in `src/` or `benchmarks/` (the Darcy flow, the compaction, the pressure diffusion, the momentum, the porous shell and the melting shell tools also check `tools/`).
It also refuses to write a value that leaves its regression band (the allowed range around the recorded value), unless you pass `--accept-physics-change`, and a value that fails its reference tolerance in all cases.

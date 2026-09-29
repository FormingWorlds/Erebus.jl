# Output Inspection and Checkpoint Resumption

`Erebus.jl` saves simulation checkpoints in the HDF5-compatible binary format [JLD2.jl](https://github.com/JuliaIO/JLD2.jl). This guide shows how to inspect output fields, post-process data, and resume simulations from saved checkpoints.

---

## Output File Structure

Output files are stored in the directory configured under `[output]` (`output_dir = "output"`). At each checkpoint interval (`savematstep`), a file named `output_<timestep:05d>.jld2` is written (for example, `output_00010.jld2`).

Each checkpoint file stores:
- Staggered grid fields:
  - `pr`: Total mixture pressure on P nodes [Pa]
  - `pf`: Fluid pore pressure on P nodes [Pa]
  - `vx`: Horizontal solid velocity on staggered Vx nodes [m/s]
  - `vy`: Vertical solid velocity on staggered Vy nodes [m/s]
  - `qxD`: Horizontal Darcy fluid flux on staggered Vx nodes [m/s]
  - `qyD`: Vertical Darcy fluid flux on staggered Vy nodes [m/s]
  - `tk1`: Temperature array on grid nodes [K]
  - `PHI`: Porosity field on P nodes [-]
  - `SXX`: Deviatoric normal stress on P nodes [Pa], where $\sigma_{yy}' = -\sigma_{xx}'$
  - `SXY`: Deviatoric shear stress on shear nodes [Pa]
- Progression state variables:
  - `timesum`: Total elapsed physical time [s]
  - `dt`: Current computational timestep [s]
  - `marknum`: Number of active Lagrangian markers in domain [-]
- Lagrangian marker property arrays:
  - Coordinate locations ($x_m, y_m$)
  - Marker temperature ($T_m$), pressure, composition, and phase fractions

---

## Loading Checkpoints in Julia

To read and inspect variables from a saved checkpoint:

```julia
using JLD2

# Open checkpoint file
file_path = "output/output_00010.jld2"
data = jldopen(file_path, "r")

# Read fields
pr = data["pr"]
pf = data["pf"]
temperature = data["tk1"]
porosity = data["PHI"]
timesum = data["timesum"]

close(data)

println("Elapsed time: ", timesum / (365.25 * 24 * 3600 * 1e6), " Ma")
println("Max temperature: ", maximum(temperature), " K")
println("Min/Max pore pressure: ", extrema(pf), " Pa")
```

---

## Calculating Effective Stress

Terzaghi effective pressure $P_{\text{eff}} = P_t - P_f$ determines matrix compaction and shear strength:

```julia
using JLD2

data = jldopen("output/output_00010.jld2", "r")
pr = data["pr"]
pf = data["pf"]
close(data)

# Compute Terzaghi effective pressure
peff = pr .- pf

println("Minimum effective pressure: ", minimum(peff), " Pa")
println("Pore fluid overpressured cells (peff < 0): ", count(peff .< 0))
```

---

## Resuming Simulations from Checkpoints

`Erebus.jl` supports exact checkpoint resumption. You can continue execution from any saved checkpoint without recomputing earlier timesteps.

### 1. Via Command Line

Pass the `--restart` flag followed by the checkpoint path:

```bash
julia --project=. launch.jl configs/hydrothermal_benchmark.toml --restart output_hydrothermal/output_00010.jld2
```

### 2. Via TOML Configuration

Specify `restart_from` in the `[output]` section of the configuration file:

```toml
[output]
output_dir   = "output_hydrothermal"
restart_from = "output_hydrothermal/output_00010.jld2"
savematstep  = 1

[time]
n_steps = 50  # Advance through step 50
```

### 3. In Julia Scripts

Pass `restart_from` as a keyword argument to `run_simulation`:

```julia
using Erebus

cfg = load_config("configs/hydrothermal_benchmark.toml")
run_simulation(cfg; restart_from="output_hydrothermal/output_00010.jld2")
```

When resuming, the solver reloads all grid fields and marker positions directly into memory, bypasses marker seeding, and resumes stepping forward from `timestep = loaded_step + 1`.

---

## Reproducibility and Restart

Erebus provides a determinism contract for simulation reproducibility and restart operations.

### Determinism Contract

Single-threaded runs and multi-threaded runs guarantee bitwise reproducibility under fixed seed values:
- Pseudo-random generation: setting `cfg.solver.seed` initializes the generator. Identical seeds generate identical initial marker positions and property assignments. Different seeds generate different realizations.
- Thread count invariance: particle-to-mesh interpolation (`p2m_mode = :buffered`) uses 16 fixed logical chunks. The accumulation order is invariant to `Threads.nthreads()`. Simulations run with 2 threads or 4 threads produce bitwise identical results.
- Process isolation: when running parameter sweeps with `run_ensemble`, each ensemble member executes in an isolated worker process with an independent seed (`spec.seed + i`).

### Checkpoint Contents and Schema

Binary checkpoint files (`.jld2`) store complete simulation state:
- Grid fields include total pressure `pr`, fluid pressure `pf`, velocities `vx` and `vy`, Darcy fluxes `qxD` and `qyD`, temperature `tk1`, porosity `PHI`, and deviatoric stresses `SXX` and `SXY`.
- Marker arrays store coordinates `xm` and `ym`, temperature `tkm`, phase assignments, and volatile inventories.
- Progress metrics record physical time `timesum`, timestep duration `dt`, active marker count `marknum`, and current step number.
- State records save progress metrics and grid and marker field data.


### Restart Semantics and Config Differences

When you resume from a checkpoint:
- The solver loads grid arrays, marker arrays, and progress metrics directly into memory.
- Marker generation and seeding routines do not run; marker replenishment draws use a reseeded generator (`MersenneTwister(cfg.solver.seed)`).
- Time integration resumes from `timestep = loaded_step + 1`.
- Grid resolution (`Nx`, `Ny`) and domain geometry (`xsize`, `ysize`) must match the checkpoint, or the solver raises a `DimensionMismatch`.
- Execution parameters can differ between the checkpoint and restart configuration:
  - Output intervals (`savematstep`, `visstep`)
  - End conditions (`time.n_steps`, `time.endtime`)
- Pass the `--force-restart-config` command-line option to apply configuration overrides on restart.

### Recomputed Fields Upon Resume

On the first timestep after restart, Erebus recomputes derived quantities and solver workspace arrays:
- Linear system matrices and preconditioners assemble fresh from the loaded state.
- Coordinate arrays and boundary condition stencils rebuild from geometry parameters.
- Hydrofracture and reaction rate limiters initialize from current field values.


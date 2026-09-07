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
  - `SXX`: Deviatoric normal stress on P nodes [Pa] (with $\sigma_{yy}' = -\sigma_{xx}'$)
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

When resuming, the solver reloads all grid fields and marker positions directly into memory, bypasses marker seeding, and resumes stepping forward from `timestep = loaded_step + 1`. The resumed run matches an uninterrupted continuous execution within relative tolerance `rtol ≈ 1e-8`.

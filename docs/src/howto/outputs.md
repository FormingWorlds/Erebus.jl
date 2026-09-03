# Output Inspection and Post-Processing

`Erebus.jl` saves simulation checkpoints in the HDF5-compatible binary format [JLD2.jl](https://github.com/JuliaIO/JLD2.jl).

---

## Output File Structure

Output files are stored in the directory configured under `[output]` (`output_dir = "output"`). At each checkpoint interval (`savematstep`), a file named `output_<timestep:05d>.jld2` is written (for example, `output_00010.jld2`).

Each file contains the primary staggered grid and marker fields:
- `pr`: Total mixture pressure array on P nodes [Pa]
- `pf`: Fluid pore pressure array on P nodes [Pa]
- `vx`: Horizontal solid velocity on staggered Vx nodes [m/s]
- `vy`: Vertical solid velocity on staggered Vy nodes [m/s]
- `qxD`: Horizontal Darcy fluid flux on staggered Vx nodes [m/s]
- `qyD`: Vertical Darcy fluid flux on staggered Vy nodes [m/s]
- `tk1`: Temperature array on grid nodes [K]
- `PHI`: Porosity field on P nodes [-]
- `SXX`: Deviatoric normal stress on basic nodes [Pa] (with $\sigma_{yy}' = -\sigma_{xx}'$)
- `SXY`: Deviatoric shear stress on shear nodes [Pa]
- `timesum`: Total elapsed physical time [s]
- `dt`: Current computational timestep [s]
- `marknum`: Number of active Lagrangian markers in domain [-]

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

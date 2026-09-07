# How to Configure Simulations

This guide shows how to configure simulation runs in `Erebus.jl` using `.toml` configuration files.

For the complete list of parameters, data types, physical units, and default values, see the [Configuration Schema Reference](../reference/config_schema.md).

---

## 1. Configure a Hydrothermal Benchmark Simulation

To set up a 2D hydrothermal circulation benchmark with radiogenic heating and water phase changes, create a `.toml` configuration file with the following core sections:

```toml
[grid]
xsize = 140000.0  # Horizontal domain width [m]
ysize = 140000.0  # Vertical domain height [m]
Nx    = 33        # Grid points in x (>= 3)
Ny    = 33        # Grid points in y (>= 3)

[geometry]
rplanet  = 50000.0  # Planetesimal radius [m]
rcrust   = 50000.0  # Crustal boundary radius [m]
xcenter  = 70000.0  # Center x-coordinate [m]
ycenter  = 70000.0  # Center y-coordinate [m]
psurface = 1000.0   # Surface pressure anchor [Pa]

[time]
dt_initial = 3168.80878  # Initial timestep [Julian yr] (1e11 s)
dt_longest = 3168.80878  # Maximum timestep [Julian yr] (1e11 s)
start_time = 2.25e6      # Start time [Julian yr] (2.25 Ma after CAIs)
endtime    = 15.0e6      # End time [Julian yr] (15.0 Ma)
start_step = 1           # Starting step index
n_steps    = 10          # Total steps to run

[thermodynamics]
hr_al            = true   # Enable 26Al radiogenic decay heating
hr_fe            = false  # Disable 60Fe decay heating
ratio_al         = 5.0e-5 # Initial 26Al/27Al isotope ratio
thermal_buoyancy = true   # Enable Darcy thermal buoyancy
```

Load and execute the configuration in Julia:

```julia
using Erebus

cfg = load_config("configs/hydrothermal_benchmark.toml")
run_simulation(cfg)
```

---

## 2. Configure Cold Surface Venting and Atmospheric Loss

To model volatile drainage across a cold planetesimal lid and couple it to kinetic escape into space, add `[venting]` and `[escape]` sections:

```toml
[venting]
active      = true                  # Enable surface boundary venting
mode        = "hydrofracture_gated" # Only vent when pore pressure breaches rock strength
k_vent      = 1.0e-11               # Surface boundary permeability [m^2]
ice_sealing = true                  # Cryogenic permeability reduction below freezing
t_freeze    = 273.15                # Freezing temperature [K]
dt_seal     = 10.0                  # Freezing transition width [K]

[escape]
active    = true   # Enable kinetic atmospheric escape
species   = "H2O"  # Volatile species ("H2O", "CO2", "N2", "CH4", "CO", "H2")
T_exobase = 200.0  # Exobase temperature [K]
```

When `ice_sealing = true`, cryogenic pore ice reduces matrix permeability below $273.15\text{ K}$. Venting activates only when pore fluid pressure breaches the cold lid.

---

## 3. Scale Grid Resolution Dynamically

`Erebus.jl` dynamically allocates coordinates and field arrays from the input configuration. You can change domain sizes and grid resolutions without recompiling:

```toml
[grid]
xsize = 200000.0  # Enlarged domain width [m]
ysize = 200000.0  # Enlarged domain height [m]
Nx    = 65        # Higher horizontal resolution
Ny    = 65        # Higher vertical resolution
```

The code instantiates `GridCoordinates(cfg.grid)` at runtime. Grid resolution must satisfy $N_x \ge 3$ and $N_y \ge 3$. Domain sizes must be strictly positive ($xsize > 0$, $ysize > 0$).

---

## 4. Manage Checkpoints and Restart Simulations

To save simulation state at regular intervals, configure the `[output]` section:

```toml
[output]
output_dir  = "output_run1"  # Directory for JLD2 output files
savematstep = 10             # Save checkpoint every 10 timesteps
visstep     = 1              # Save visualization outputs every step
```

To resume a simulation from a saved checkpoint, specify `restart_from`:

```toml
[output]
output_dir   = "output_run1"
restart_from = "output_run1/output_00010.jld2"
savematstep  = 10
visstep      = 1

[time]
n_steps = 50  # Advance through step 50
```

When `restart_from` is non-empty, `run_simulation` reloads marker distributions, field arrays, and time counters directly from the JLD2 checkpoint. The simulation automatically resumes from the saved checkpoint step (`timestep + 1`).

---

## 5. Model Protoplanetary Disk Thermal Evolution

To simulate a planetesimal embedded in an evolving protoplanetary disk, configure the `[disk]` section:

```toml
[disk]
enabled             = true                 # Enable disk temperature evolution
model               = "class1_to_class2"   # Two-stage accretion-to-clearing model
orbital_distance_au = 2.5                  # Planetesimal semi-major axis [AU]
stellar_mass_msun   = 1.0                  # Central star mass [Solar masses]
t_dispersal_myr     = 3.0                  # Disk dispersal epoch [Myr]
dt_dispersal_myr    = 0.5                  # Dispersal transition duration [Myr]
p_amb_disk          = 1.0                  # Nebular gas pressure [Pa]
p_amb_space         = 1.0e-4               # Post-dispersal vacuum pressure [Pa]
albedo              = 0.06                 # Post-dispersal Bond albedo
```

During the simulation, the external boundary conditions smoothly transition from gas-rich disk conditions to solar radiative equilibrium in vacuum space.

---

## Further Reading

- [Configuration Schema Reference](../reference/config_schema.md): Exhaustive parameter definitions, data types, and units.
- [Running Simulations](running.md): Instructions for command-line execution and batch processing.
- [Outputs and Checkpoints](outputs.md): How to analyze JLD2 files and inspect marker data.

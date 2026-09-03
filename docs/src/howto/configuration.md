# Configuration Guide

`Erebus.jl` uses structured `.toml` configuration files to define all physical, geometric, solver, and storage parameters. This allows running parameter space explorations and reproducing simulation runs without modifying source code.

---

## Configuration File Sections

A simulation configuration is organized into eight logical groups:

### 1. `[grid]`
Defines the numerical domain dimensions and basic grid resolution.

```toml
[grid]
xsize = 14000.0   # Horizontal domain size [m]
ysize = 14000.0   # Vertical domain size [m]
Nx = 15           # Basic grid resolution in x
Ny = 15           # Basic grid resolution in y
```

*Note: In the current version, grid resolution `Nx` and `Ny` are coupled to compile-time static stencils. Changes to grid dimensions require updating constants and recompiling.*

### 2. `[geometry]`
Defines the spherical planetesimal geometry and location within the Cartesian computational domain.

```toml
[geometry]
rplanet = 5000.0   # Total planetesimal radius [m]
rcrust = 4800.0    # Crust radius [m] (transition from core/mantle to crust)
xcenter = 7000.0   # Horizontal center coordinate [m]
ycenter = 7000.0   # Vertical center coordinate [m]
psurface = 1.0e+3  # Surface ambient pressure anchor [Pa]
```

### 3. `[time]`
Controls the timestepping and temporal integration.

```toml
[time]
dt_initial = 1.0e+11       # Initial timestep [s]
dt_longest = 1.0e+11       # Maximum allowable timestep [s]
dtcoefdn = 0.5             # Timestep reduction factor upon non-convergence
dtcoefup = 1.2             # Timestep growth factor upon rapid convergence
dtstep = 200               # Iteration interval before increasing timestep
dxymax = 0.05              # Maximum marker displacement per timestep [grid units]
vpratio = 0.333333333333   # Marker velocity weighting (staggered vs node)
DTmax = 20.0               # Maximum temperature change per step [K]
start_time = 7.10046e+13   # Initial time [s] (e.g. 2.25 Ma after CAIs)
endtime = 4.73364e+14       # Final simulation time [s] (e.g. 15 Ma)
start_step = 1             # Starting step counter
n_steps = 10               # Total timesteps to compute
```

### 4. `[solver]`
Controls the Picard iteration loop, yielding tolerances, and sparse matrix solver.

```toml
[solver]
titermax = 10000     # Maximum global thermo-mechanical iterations
nplast = 100000      # Maximum plastic yielding iterations
yerrmax = 100.0      # Yielding tolerance level
etawt = 0.0          # Viscosity relaxation weight
dphimax = 100.01     # Maximum allowable porosity change ratio per step
seed = 42            # Random seed for marker spatial distribution
use_pardiso = false  # Set true to use Pardiso solver, false for SuiteSparse UMFPACK
etaphikoef = 1.0     # Bulk viscosity scaling factor
etamin = 1.0e+12     # Lower viscosity cutoff [Pa s]
etamax = 1.0e+23     # Upper viscosity cutoff [Pa s]
```

### 5. `[poroelasticity]`
Sets solid matrix and fluid pore compressibility and porosity bounds.

```toml
[poroelasticity]
betasolid = 2.5e-11   # Solid silicate matrix compressibility [1/Pa]
betafluid = 4.0e-10   # Pore fluid (water) compressibility [1/Pa]
phimin = 1.0e-4       # Minimum porosity threshold floor [-]
phimax = 0.9999       # Maximum porosity threshold ceiling [-]
```

### 6. `[thermodynamics]`
Controls radiogenic isotope abundances, half-lives, decay energies, and phase change parameters.

```toml
[thermodynamics]
hr_al = true           # Enable 26Al decay heating
hr_fe = false          # Enable 60Fe decay heating
ratio_al = 5.0e-5      # Initial 26Al/27Al isotope ratio at CAI formation
ratio_fe = 1.0e-6      # Initial 60Fe/56Fe isotope ratio at CAI formation
tmsolidphase = 1416.0  # Silicate solidus melting temperature [K]
tmfluidphase = 273.0   # H2O ice melting temperature [K]
Lᶠ = 333.55e+3         # Latent heat of melting for water ice [J/kg]
phim0 = 0.2            # Standard reference porosity [-]
thermal_buoyancy = true # Enable Darcy thermal buoyancy
fluid_viscosity_mode = "arrhenius" # Viscosity mode: "arrhenius" or "constant"
fluid_viscosity_Ea = 15.0e+3       # Activation energy [J/mol]
fluid_viscosity_T0 = 293.15        # Reference temperature [K]
fluid_viscosity_eta0 = 1.0e-3      # Reference viscosity at T0 [Pa s]
```

### 7. `[materials]`
Specifies material phase properties on the 3-phase staggered grid:
- Index 1: Planetesimal core / mantle (`rmark <= rcrust`)
- Index 2: Porous silicate crust / rock (`rcrust < rmark < rplanet`)
- Index 3: Sticky air / space (`rmark >= rplanet`)

*Note: In the current release, eight material arrays (`rhosolidm`, `rhofluidm`, `etasolidm`, `etasolidmm`, `etafluidm`, `etafluidmm`, `ksolidm`, `kfluidm`) and radii (`rplanet`, `rcrust`) are locked to compiled constants in `src/constants.jl`. Modifying them requires recompilation. The other 10 material arrays can be configured freely in TOML.*

```toml
[materials]
rhosolidm = [3300.0, 3300.0, 1.0]
rhofluidm = [1000.0, 1000.0, 1.0]
etasolidm = [1.0e+19, 1.0e+19, 1.0e+16]
etasolidmm = [1.0e+19, 1.0e+19, 1.0e+16]
etafluidm = [1.0e+12, 1.0e+12, 1.0e-03]
etafluidmm = [1.0e-03, 1.0e-03, 1.0e-03]
rhocpsolidm = [3.3e+06, 3.3e+06, 3.0e+06]
rhocpfluidm = [1.0e+06, 1.0e+06, 3.0e+06]
alphasolidm = [3.0e-05, 3.0e-05, 0.0]
alphafluidm = [5.0e-05, 5.0e-05, 0.0]
ksolidm = [3.0, 3.0, 3000.0]
kfluidm = [50.0, 50.0, 3000.0]
gggsolidm = [1.0e+10, 1.0e+10, 1.0e+10]
frictsolidm = [0.6, 0.6, 0.0]
cohessolidm = [1.0e+08, 1.0e+08, 1.0e+08]
tenssolidm = [6.0e+07, 6.0e+07, 6.0e+07]
kphim0 = [1.0e-13, 1.0e-13, 1.0e-17]
tkm0 = [170.0, 170.0, 170.0]
```

### 8. `[output]`
Configures file output paths and checkpoint cadences.

```toml
[output]
output_dir = "output"   # Directory for simulation JLD2 checkpoints
savematstep = 10        # Save full checkpoint every N timesteps
visstep = 1             # Visualization output frequency
```

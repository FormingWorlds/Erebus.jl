# Tutorial: Planetesimal Differentiation

This tutorial walks through modeling the complete differentiation of an early planetesimal in `Erebus.jl`. The simulation traces the evolution of a cold, homogeneous planetesimal through radiogenic heating, ice melting, hydrothermal circulation, surface hydrofracture venting, silicate rock melting, volatile degassing and atmospheric escape, and Fe-FeS metallic core segregation.

---

## Physical Overview

We model a planetesimal with radius $R_{\text{planet}} = 50\,000\text{ m}$ ($100\text{ km}$ diameter) in a $140\,000\text{ m} \times 140\,000\text{ m}$ Cartesian domain. The initial planetesimal represents a primordial carbonaceous mixture:
- **Silicate rock matrix**: $\sim 50\text{ vol}\%$
- **Water ice in pore network**: $\sim 30\text{ vol}\%$ ($\phi_0 = 0.30$)
- **Dispersed Fe-FeS metal**: $\sim 20\text{ vol}\%$ ($X_{\text{fe,bulk}} = 0.20$)
- **Initial temperature**: $T_0 = 150\text{ K}$ (or preheated to $1350\text{ K}$ in core segregation benchmarks)

```text
    Initial State (t = 0 Ma)              Differentiated State (t = 3.5 Ma)
   +------------------------+           +------------------------+
   |                        |           |       Sticky Air       |
   |      Cold Mixture      |   ====>   |  +------------------+  |
   |   (Ice + Rock + Metal) |           |  | Primordial Crust |  |
   |        T = 150 K       |           |  |  Silicate Mantle |  |
   |                        |           |  |   (Fe-FeS Core)  |  |
   +------------------------+           +--+------------------+--+
```

---

## The Evolutionary Sequence

### 1. Early Radiogenic Heating ($0.0 \text{ to } 0.3\text{ Ma}$)

Heating is driven by the short-lived radionuclides $^{26}\text{Al}$ ($t_{1/2} = 0.717\text{ Ma}$) and $^{60}\text{Fe}$ ($t_{1/2} = 2.6\text{ Ma}$):
- Rock matrix heating: $^{26}\text{Al}$ decay power deposits directly into the solid rock matrix.
- Metallic phase heating: $^{60}\text{Fe}$ decay power deposits into the metallic iron phase ($X_{\text{fe,bulk}}$).

Because the interior is thermally insulating ($k_{\text{rock}} \approx 2.5\text{ W/(m}\cdot\text{K)}$), internal temperatures rise rapidly while the surface radiates to space.

### 2. Pore Ice Melting and Hydrothermal Circulation ($0.3 \text{ to } 0.8\text{ Ma}$)

When the deep interior exceeds $T = 273.15\text{ K}$, pore ice melts to liquid water:
- Dynamic viscosity drops by 15 orders of magnitude, from ice ($\eta_{\text{ice}} \approx 10^{12}\text{ Pa}\cdot\text{s}$) to liquid water ($\eta_{\text{water}} \approx 10^{-3}\text{ Pa}\cdot\text{s}$).
- Thermal expansion of water ($\alpha_f \approx 5 \times 10^{-5}\text{ K}^{-1}$) produces Darcy thermal buoyancy, driving convective circulation in the permeable interior.

### 3. Hydrofracture Venting and Cryovolcanism ($0.5 \text{ to } 1.5\text{ Ma}$)

Thermal expansion and prograde dehydration reactions generate high pore fluid pressures ($P_f$). When the effective pressure exceeds the tensile strength of the overlying cold rock lid ($P_t - P_f \le -\sigma_t$), hydraulic fractures breach the cold lid:
- Overpressured fluids discharge through surface fractures.
- Sublimation latent cooling ($L_{\text{sub}} \approx 2.83\times 10^6\text{ J/kg}$) buffers surface temperatures against runaway heating.

### 4. Silicate Melting and Magma Rheology ($1.0 \text{ to } 2.5\text{ Ma}$)

As temperatures surpass the silicate solidus ($T > 1400\text{ K}$), partial melting of the rock matrix initiates:
- Solid matrix viscosity decreases by orders of magnitude via the Costa et al. (2009) rheology formulation.
- A central magma ocean develops beneath the conductive crust.

### 5. Multi-Species Volatile Degassing and Atmospheric Escape

Dissolved volatiles partition between silicate melt and vapor phases:
- Water exsolves according to the Burnham-Dixon square-root solubility law.
- Under reducing conditions ($\Delta\text{IW} \le 0$), nitrogen dissolves chemically as nitride ($\text{N}^{3-}$).
- Homogeneous speciation computes equilibrium abundances for ten gas species ($\text{H}_2$, $\text{H}_2\text{O}$, $\text{CO}$, $\text{CO}_2$, $\text{CH}_4$, $\text{N}_2$, $\text{NH}_3$, $\text{H}_2\text{S}$, $\text{S}_2$, $\text{SO}_2$).
- Low-gravity escape strips light volatiles through a cubic Hermite transition between sound-speed hydrodynamic blow-off and kinetic Jeans effusion.

### 6. Metallic Core Formation ($1.5 \text{ to } 3.5\text{ Ma}$)

When temperatures reach the Fe-FeS eutectic melting point ($T_{\text{eutectic}} = 1213\text{ K}$):
1. **Porous Percolation ($F_{\text{silicate}} < F_{\text{settle\_start}} = 0.40$)**: Liquid Fe-FeS drains downward through permeable pore channels in the solid silicate matrix.
2. **Transition Regime ($0.40 \le F_{\text{silicate}} \le 0.50$)**: A $C^1$-continuous cubic Hermite polynomial smoothly blends percolation and settling velocities as melt pockets coalesce.
3. **Stokes Droplet Settling ($F_{\text{silicate}} > F_{\text{perc\_end}} = 0.50$)**: In regions where silicates melt extensively into a magma ocean, metal droplets separate and settle rapidly at terminal Stokes velocities modified by Richardson-Zaki hindrance.
4. **Weber Droplet Equilibrium**: Droplet diameters equilibrate between surface tension and hydrodynamic shear forces.
5. **Core Compaction**: Settling metal accumulates at the planetary center, packing to $\phi_{\text{pack}} = 0.65$ and releasing gravitational potential energy as dissipation heating.

---

## Simulation Configuration

The full 2D evolutionary history of a planetesimal from a cold primordial state ($T_0 = 150\text{ K}$) to complete differentiation spans several million years ($0 \text{ to } 3.5\text{ Ma}$) and requires thousands of computational timesteps on high-performance computing clusters.

To demonstrate the numerical formulation without long execution times, this tutorial uses a focused core formation benchmark (`configs/core_formation_benchmark.toml`). This benchmark initializes the interior in a preheated state ($T = 1350\text{ K}$ above the Fe-FeS eutectic at $t = 2.25\text{ Ma}$) and integrates 5 representative timesteps ($\sim 16\text{ kyr}$) to isolate and verify the Darcy percolation, Stokes settling, and core compaction dynamics.

Below is the complete benchmark configuration:

```toml
[grid]
xsize = 140000.0        # horizontal domain size [m]
ysize = 140000.0        # vertical domain size [m]
Nx = 33                  # basic grid resolution in x
Ny = 33                  # basic grid resolution in y

[geometry]
rplanet = 50000.0       # planetesimal radius [m]
rcrust = 50000.0        # crust radius [m]
xcenter = 70000.0       # horizontal center [m]
ycenter = 70000.0       # vertical center [m]
psurface = 1.0e+3        # surface pressure anchor [Pa]

[time]
dt_initial = 3168.80878      # initial computational timestep [yr] (~1e11 s)
dt_longest = 3168.80878      # maximum allowed computational timestep [yr] (~1e11 s)
dtcoefdn = 0.5               # coefficient to decrease computational timestep
dtcoefup = 1.2               # coefficient to increase computational timestep
dtstep = 200                 # iterations before changing timestep
dxymax = 0.05                # max marker displacement per timestep [grid units]
vpratio = 0.3333333333333333 # weight of averaged velocity for moving markers
DTmax = 20.0                 # max temperature change per timestep [K]
start_time = 2.25e6          # initial simulation time [yr] (2.25 Ma)
endtime = 15.0e6             # maximum simulation time [yr] (15.0 Ma)
start_step = 1               # starting timestep counter
n_steps = 5                  # total number of timesteps to run

[materials]
tkm0 = [1350.0, 1350.0, 170.0]  # planet interior preheated above Fe-FeS eutectic (1213 K)

[thermodynamics]
hr_al = true             # 26Al decay heating active in solid phase
hr_fe = true             # 60Fe decay heating active in solid phase

[melting]
active = true            # silicate melt tracking active

[coreformation]
percolation_active = true
settling_active = true
sulfur_fraction = 0.31
metal_density_mode = "sanloup2000"
rho_metal = 5450.0
rho_metal_solid = 5700.0
eta_metal = 1.0e-2
k_metal = 40.0
rhocp_metal = 4.0e6
Xfe_bulk = 0.20
phi_pack = 0.65
T_eutectic = 1213.0
dT_metal = 50.0
k_metal_ref = 1.0e-9
perm_exponent = 3.0
phi_crit_perc = 0.05
phi_residual = 0.02
phi0 = 0.1
droplet_size_mode = "capillary_mean"
droplet_diameter_fixed = 5.0e-3
sigma_metal_silicate = 1.0
We_crit = 10.0
hindered_exponent = 4.5
hadamard_rybczynski = false
F_settle_start = 0.40
F_perc_end = 0.50
segregation_heating = true
cfl_settling = 0.5
max_subcycles = 2000

[output]
output_dir = "output_core_benchmark"
savematstep = 1
```

---

## Running the Simulation

Execute the benchmark simulation in Julia:

```julia
using Erebus

# Load configuration
cfg = load_config("configs/core_formation_benchmark.toml")

# Run coupled simulation
run_simulation(cfg)
```

---

## Analyzing Results

Examine the benchmark checkpoint using `load_state`:

```julia
using Erebus

# Load checkpoint
data = load_state("output_core_benchmark/output_00005.jld2")

# Temperature field on staggered grid [K]
T_grid = data["tk2"]
max_T = maximum(T_grid)
println("Peak interior temperature: ", round(max_T, digits=1), " K")

# Marker metal volume fraction
Xfe = data["Xfe_bulk"]
println("Max core metal volume fraction: ", round(maximum(Xfe), digits=3))

# Marker porosity
phim = data["phim"]
println("Mean marker porosity: ", round(sum(phim) / length(phim), digits=3))
```

In full-scale production simulations integrated through the entire 3.5 Ma differentiation history, the planetesimal evolves into three distinct physical zones:
1. **Central Metallic Core**: Dense Fe-FeS metal pool packed to $\phi_{\text{fe}} \approx 0.65$.
2. **Depleted Silicate Mantle**: High-temperature partially molten silicate shell stripped of metallic iron.
3. **Conductive Primordial Crust**: Cold, non-melted exterior retaining primordial ice, phyllosilicates, and fine-grained metal.

---

### Full Multi-Million Year Production Runs

To simulate the complete 0 to 3.5 Ma sequence from a primordial cold start:
- Set `materials.tkm0 = [150.0, 150.0, 150.0]` for an initial cold mixture.
- Set `time.start_time = 0.0` and `time.endtime = 3.5e6` (years).
- Set `time.n_steps = 10000` with adaptive timestepping enabled.
- Enable `[venting]` (`active = true`, `mode = :hydrofracture_gated`) and `[escape]` (`active = true`, `multi_species = true`).

---

## Related Documentation

- [Iron Core Formation & Metal Segregation (Explanations)](../explanations/core_formation.md)
- [Iron Core Formation Verification (Validation)](../validation/core_formation.md)
- [Silicate Melting & Soft Turbulence (Explanations)](../explanations/rock_melting.md)
- [Degassing & Cold Venting (Explanations)](../explanations/degassing_and_venting.md)
- [Configuration Schema Reference](../reference/config_schema.md)


# Model Overview

`Erebus.jl` models the internal physical and geochemical evolution of icy and rocky planetesimals during the first tens of millions of years of Solar System history.

---

## The Physical Scenario and Evolution Phases

### How does a planetesimal evolve from a cold porous body into a differentiated world?

Planetesimals evolve across 3 physical phases driven by radiogenic decay heating: Phase I (cold accretion), Phase II (hydrothermal alteration), and Phase III (magma ocean and core formation).

Figure 2 shows these physical phases and links them to the active processes in Erebus.

![Erebus.jl Planetesimal Evolution in Three Phases](../assets/erebus_planetesimal_stages.svg)

*Figure 2: Planetesimal evolution across three physical phases. Sectors show successive phases of 1 body; sector widths do not scale with time. Left and right cards show physical processes that Erebus solves. Active modules depend on configuration switches; all modules are shown active. Flux arrows are schematic. (IW: Iron-Wüstite buffer; QFM: Quartz-Fayalite-Magnetite buffer; XUV: extreme-ultraviolet irradiation). Schematic (hand-built diagram, not simulation output).*

In Figure 2:
- **Phase I (Cold Accretion)**: Low temperatures preserve porous ice and pristine rock. Accretion flux $\dot{M}_\mathrm{acc}$ delivers pebbles and volatile ices to the surface.
- **Phase II (Hydrothermal Alteration)**: Radiogenic decay warms the interior. Ice melts into pore water. Darcy percolation $\mathbf{q}_\mathrm{D}$ drives hydrothermal circulation, serpentinization, and dehydration reactions through permeable crust.
- **Phase III (Magma Ocean & Core Segregation)**: Deep heating drives silicate melting above the silicate solidus ($T_\mathrm{s}$). Above the critical rheological transition ($F_\mathrm{m} > F_\mathrm{crit} \approx 0.40$), the mixture behaves as a liquid magma ocean. Liquid metal droplets segregate to form a metallic core, releasing segregation heating. Magma ocean surface heat flux $F_\mathrm{rad}$ radiates into space.

### Physical Scenario Sequence

Planetesimals form in the protoplanetary disk through gravitational collapse of pebble swarms or streaming instabilities.
Bodies that accrete during the first 2 to 3 million years after Calcium-Aluminum-rich Inclusions (CAIs) incorporate short-lived radioactive isotopes: $^{26}\text{Al}$ (half-life $t_{1/2} \approx 0.717\text{ Ma}$) and $^{60}\text{Fe}$ ($t_{1/2} \approx 2.62\text{ Ma}$).
As radionuclides decay, volumetric heating drives 7 physical processes:

1. **Ice Melting (Phase I to II)**:
   Initial temperatures of $\approx 170\text{ K}$ rise.
   At $273\text{ K}$, primordial water ice melts, absorbing latent heat ($L^f \approx 333\text{ kJ/kg}$) and generating mobile pore water in the porous silicate matrix.

2. **Hydrothermal Circulation (Phase II)**:
   Liquid water percolates through the pore network under Darcy flow.
   Flow is driven by thermal buoyancy and compaction-induced fluid pressure gradients.

3. **Clay Dehydration and Gas Generation (Phase II)**:
   At elevated temperatures ($T \approx 500\text{ to } 900\text{ K}$), hydrous phyllosilicates (such as serpentine) undergo thermal dehydration reactions, releasing structurally bound water into the pore fluid.

4. **Matrix Compaction and Pore Overpressure (Phase II)**:
   Viscous, elastic, and plastic deformation of the silicate matrix reduces porosity $\phi$.
   When matrix compaction exceeds fluid escape via Darcy percolation, pore fluid pressure $P_\mathrm{f}$ rises toward total pressure $P_\mathrm{t}$.

5. **Failure and Hydrofracture (Phase II)**:
   When Terzaghi effective stress becomes tensile ($\sigma_\text{eff} \le -\sigma_\text{tensile}$), the silicate matrix fractures.
   This creates high-permeability pathways that vent pore fluids to the planetesimal surface.

6. **Silicate Rock Melting and Magma Ocean (Phase III)**:
   In bodies with early accretion, decay heating raises temperatures past the silicate solidus ($T_\mathrm{s} \approx 1400\text{ K}$).
   Silicate rock melts into a crystal-liquid mixture.
   Above critical melt fraction $F_\mathrm{crit} \approx 0.40$ (`melting.phi_crit`), the mixture behaves as a liquid magma ocean.
   Soft turbulence convection transports heat outward to the surface, where the magma ocean radiates heat flux $F_\mathrm{rad}$ into space.

7. **Iron Core Formation and Metal Segregation (Phase III)**:
   Metal means liquid Fe-FeS alloy (iron plus iron sulfide).
   At temperatures above the Fe-FeS eutectic ($T_\text{eut} \approx 1213\text{ K}$), dense metallic liquid melts before silicate rock melts completely.
   Metal moves by one of two methods, depending on the silicate melt fraction $F_\mathrm{m}$.
   Below 0.40, metal trickles through solid rock via porous percolation.
   Above 0.50, metal droplets sink through magma via Stokes settling.
   Between 0.40 and 0.50, the code smoothly blends the two methods (`coreformation.F_settle_start` and `coreformation.F_perc_end`).
   Liquid metal droplets segregate to form a metallic core, releasing gravitational dissipation heating that accelerates core growth.

### Summary of Physical Processes by Evolution Phase

Table 2 lists the symbols in Figure 2.

| Phase | Symbol | Meaning | Unit | Where in the code |
|:---|:---|:---|:---|:---|
| **Phase I** | $\dot{M}_\mathrm{acc}$ | Pebble and embryo accretion mass rate | $\text{kg/s}$ | `accrete!` (`src/simulation/steps/accrete.jl`) |
| | $T_\mathrm{amb}$ | Protoplanetary disk ambient gas temperature | $\text{K}$ | `prepare_step_ambient!` (`src/simulation/steps/prepare_step.jl`) |
| **Phase II** | $\mathbf{q}_\mathrm{D}$ | Darcy pore fluid velocity vector | $\text{m/s}$ | `reconstruct_darcy_fluxes!` (`src/numerics/stokes_darcy.jl`) |
| | $^{26}\text{Al}, ^{60}\text{Fe}$ | Short-lived radionuclides driving decay heating | dimensionless | `thermodynamics.f_al`, `thermodynamics.f_fe` (`src/config.jl`) |
| **Phase III** | $\mathbf{J}_\mathrm{Fe}$ | Metal segregation volume flux vector | $\text{m/s}$ | `apply_metal_segregation!` (`src/numerics/segregation.jl`) |
| | $D_i$ | Metal-silicate trace element partition coefficients | dimensionless | `apply_metal_segregation!` (`src/numerics/segregation.jl`) |
| | $F_\mathrm{m}$ | Silicate melt fraction | dimensionless | `compute_melt_fraction` (`src/physics/melting_convection.jl`) |
| | $F_\mathrm{crit}$ | Critical rheological melt fraction threshold | dimensionless | `melting.phi_crit` (`src/config.jl`) |
| | $k_\mathrm{eff}$ | Turbulent convective thermal conductivity | $\text{W/}(\text{m}\cdot\text{K})$ | `regularized_soft_turbulence_conductivity` (`src/physics/melting_convection.jl`) |
| | $F_\mathrm{rad}$ | Surface radiative heat flux into space | $\text{W/m}^2$ | `apply_surface_radiation!` (`src/simulation/steps/prepare_step.jl`) |
| | $T(\tau_\mathrm{IR})$ | Atmospheric temperature at infrared optical depth $\tau_\mathrm{IR}$ | $\text{K}$ | `compute_guillot_surface_temperature` (`src/atmosphere.jl`) |
| | $\tau_\mathrm{IR}$ | Thermal infrared optical depth | dimensionless | `compute_guillot_surface_temperature` (`src/atmosphere.jl`) |

---

## Global Architecture and Execution Pipeline

### What does the code do in one time step?

The simulation loop advances physical time by executing 10 operator-splitting stages per time step.
Each stage updates temperatures, mechanical flows, or geochemical states.
Figure 1 shows the 10 stages and their governing balance equations.

### Execution Pipeline and Governing Equations

![Erebus.jl Execution Pipeline and Governing Equations](../assets/erebus_execution_pipeline.svg)

*Figure 1: Execution pipeline and governing equations in Erebus. The center column shows the full run, where the boxes inside the frame form 1 time step. Side cards show governing equations and physical modules for each stage. The configuration file switches each module on or off; all modules are shown active. The outer iteration repeats stages 4, 5, and 6 until fields converge. If the plastic Stokes-Darcy iteration fails to converge, the step restarts from stage 1 with a reduced time step. Schematic (hand-built diagram, not simulation output).*

Each marker carries the out-of-plane spherical weight $w_\mathrm{3D} = 2r$; see [Geometry Convention](#geometry-convention) for derivation and limits.

### Numerical Execution Pipeline Stages

The simulation loop advances physical time in 10 sequential stages:

1. **Ambient & Accretion (`prepare_step_ambient!`, `accrete!`)**:
   The code computes disk ambient temperature and pressure.
   It deposits accreted pebbles or embryos onto the body.
   It advances planetary radius $R(t)$ and updates sticky air boundaries.

2. **Radiogenic Heating (`radiogenic_heating!`, `update_pyrolysis!`)**:
   The code computes analytical decay rates of short-lived radionuclides ($^{26}\text{Al}$, $^{60}\text{Fe}$).
   It adds volumetric heat sources to markers.
   It tracks volatile release from organic carbon pyrolysis.

3. **Marker Mapping & Gravity (`interpolate_markers_to_grid!`, `solve_gravity!`, `apply_surface_radiation!`)**:
   The code interpolates marker properties to staggered grid nodes.
   It computes silicate melt fractions $F_\mathrm{m}$.
   It solves the Poisson equation for gravitational potential $\Phi$ and gravity $\mathbf{g}$.
   It applies radiative surface boundary conditions.

4. **Hydration Reactions (`perform_thermochemical_reaction!`, outer iteration)**:
   The code solves reversible serpentinization and dehydration reactions through time-relaxed equilibrium with timescale $\Delta t_\mathrm{r}$.
   The code defines free energy change $\Delta G$ for the reaction $\mathrm{W_{sil}} \to \mathrm{D_{sil}} + \mathrm{H_2O}$ (dehydration).
   It updates mineral hydration states, fluid production, and reaction latent heat.

5. **Hydro-Mechanical Stokes-Darcy (`assemble_hydromechanical_lse!`, outer iteration)**:
   The solver assembles and solves the linear system for solid velocities, total pressure, fluid Darcy flux, and fluid pressure.
   It updates matrix compaction and plastic yielding.

6. **Thermal & Segregation (`solve_thermal_energy!`, `apply_metal_segregation!`)**:
   The code solves the thermal energy equation with conduction and latent heat.
   In Stage 9, marker advection transports matrix and Darcy pore fluid heat.
   Turbulent convection and optional hydrothermal circulation enter through effective thermal conductivity $k_\mathrm{eff}$.
   The code computes metal droplet settling and silicate melt ascent using conservative drift-flux transport.
   It adds heat $H_\mathrm{seg}$ released by metal segregation.

7. **Marker Stress & Subgrid Diffusion (`diffuse_and_update_markers!`)**:
   The code interpolates grid stresses to markers and applies Maxwell viscoelastic relaxation.
   It updates marker dynamic viscosity.
   When $d_s > 0$ (`solver.dsubgrids`), it applies subgrid stress diffusion; $d_s = 0$ (default) switches stress diffusion off.
   The code also applies subgrid temperature diffusion (`solver.dsubgridt`).

8. **Fluid Venting & Atmosphere (`vent_and_degas!`, `evolve_atmosphere!`)**:
   The code computes hydrofracture fluid venting and volatile exsolution from magma.
   Venting drainage rate $S_\mathrm{vent}$ is computed from surface Darcy flux in Stage 5 (`apply_venting_surface_boundary!`) and applied to drain marker pores in Stage 8.
   It calls `advance_marker_thermo_porosity_venting!` to update marker temperatures and compaction porosity.
   It couples vented volatiles into surface atmospheric pressure and equilibrium temperature.
   It computes photo-evaporative and Jeans escape losses to space.

9. **Marker Advection, Replenish & Domain Doubling (`advect_markers!`, `replenish!`, `telescope_domain!`)**:
   The code advects markers with a 4th-order Runge-Kutta scheme (RK4).
   It transports matrix heat and Darcy pore fluid heat along marker paths.
   It updates nodal total and fluid pressures.
   It replenishes markers in empty cells.
   The code doubles domain size when the body radius exceeds 70% of domain half-width ($0.70 \cdot x_\mathrm{size}/2$; see [Telescoping Domain](telescoping_domain.md)).

10. **Diagnostics & Output (`advance_step_diagnostics!`)**:
    After telescoping, the loop advances physical time $t \leftarrow t + \Delta t$.
    The code logs conservation metrics, saves JLD2 checkpoints, and streams telemetry records.

### Outer Iteration and Timestep Control

Stages 4, 5, and 6 form an outer iteration loop that repeats up to $n_\mathrm{iter}$ (`solver.max_plastic_iterations`):
- **Convergence criteria**: The loop stops when all 3 conditions are true:
  1. The fluid pressure error is below `reaction.pferrmax`.
  2. The loop ran more than 2 iterations, or the reaction rate is zero (`DMP == 0`).
  3. The maximum temperature change is at most $\Delta T_\mathrm{max}$ (`time.DTmax`).
  If the loop does not converge within `solver.max_plastic_iterations`, the step continues with `thermochemical_converged = false`.
- **Plastic convergence and retries**: The code saves a snapshot (a copy of the state at the start of the step) before the outer loop. If plastic iterations fail to converge, the step restores the snapshot, sets $\Delta t = \min(\Delta t/2, \Delta t_\mathrm{next})$, and restarts from Stage 1. The simulation stops after `solver.max_dt_reductions` reductions.
- **Accretion timestep resync**: If `cfg.accretion.active` is true and $\Delta t$ changed during the step, the code restores the snapshot and restarts the step once at the reduced $\Delta t$ to preserve accretion rate synchronization.

---

## Model Capabilities and Scope

### Which physical and numerical capabilities does Erebus include?

Erebus simulates the coupled thermo-hydro-mechanical evolution of porous planetesimals. Key capabilities include:

- Solves 2-D Cartesian flow for planetesimals in sticky air.
- Solves viscous and plastic rock matrix flow coupled to Darcy fluid flow with Biot poroelastic compressibility and Drucker-Prager failure.
- Solves heat conduction, fluid flow, radiogenic heat ($^{26}\text{Al}$, $^{60}\text{Fe}$), reaction and melt latent heats, shear heating, and turbulent heat flow.
- Tracks rock hydration and dehydration, surface vents, and dynamic hydrofracture permeability when pore fluid pressure exceeds rock strength.
- Models rock melting, melt ascent, and liquid metal settling to a central core via Darcy flow and Stokes settling.
- Conserves elemental budgets (H, C, N, S, O) with ice mixtures, organics, melt solubility, equilibrium gas release, metal-silicate partition, and Iron-Wüstite redox buffers.
- Couples a 1D Guillot radiative profile with gas chemistry, photo-evaporative loss, and Jeans escape.
- Models body growth from pebble and embryo impacts, with domain doubling at constant cell size.
- Computes gravity fields from a 2-D Poisson solve or an enclosed-mass radial profile.
- Reads TOML files, writes restart files, streams run metrics, and executes Latin Hypercube parameter sweeps.
- Supports exact bitwise checks for resumed runs on all markers, grid fields, atmosphere, and random seed states.

---

## Symbols and Governing Equations by Stage

### What do the symbols and equations in Figure 1 represent?

The tables below list the symbols in Figure 1, grouped by execution stage.

#### Stage 1: Ambient & Accretion
Module: `src/simulation/steps/prepare_step.jl`, `src/simulation/steps/accrete.jl`

| Symbol | Meaning | Unit | Where in the code |
|:---|:---|:---|:---|
| $w_\mathrm{3D}$ | Out-of-plane cylindrical-to-spherical marker weighting factor ($2r$) | $\text{m}$ | `marker_out_of_plane_length` (`src/geometry.jl`) |
| $r$ | Radial distance from planetesimal center | $\text{m}$ | `marker_out_of_plane_length` (`src/geometry.jl`) |
| $A$ | Cross-sectional computational area | $\text{m}^2$ | `marker_area` (`src/geometry.jl`) |
| $R$ | Outer planetesimal radius | $\text{m}$ | `geometry.rplanet` (`src/config.jl`) |

#### Stage 2: Radiogenic Heating
Module: `src/physics/radioactivity.jl`, `src/simulation/steps/radiogenic_heating.jl`

| Symbol | Meaning | Unit | Where in the code |
|:---|:---|:---|:---|
| $Q_i$ | Specific radiogenic decay heating for isotope $i$ | $\text{W/kg}$ | `Q_radiogenic` (`src/physics/radioactivity.jl`) |
| $f_i$ | Number of host-element atoms per kg of bulk rock | $\text{kg}^{-1}$ | `thermodynamics.f_al`, `thermodynamics.f_fe` (`src/config.jl`) |
| $r_i$ | Initial isotopic abundance ratio ($^{26}\text{Al}/^{27}\text{Al}$, $^{60}\text{Fe}/^{56}\text{Fe}$) | dimensionless | `thermodynamics.ratio_al`, `thermodynamics.ratio_fe` (`src/config.jl`) |
| $E_i$ | Decay energy per decay event | $\text{J}$ | `thermodynamics.E_al`, `thermodynamics.E_fe` (`src/config.jl`) |
| $t_{1/2, i}$ | Isotope half-life | $\text{s}$ | `thermodynamics.t_half_al`, `thermodynamics.t_half_fe` (`src/config.jl`) |
| $t$ | Elapsed physical time since CAI formation | $\text{s}$ | `radiogenic_heating!` (`src/simulation/steps/radiogenic_heating.jl`) |
| $i$ | Radionuclide system index ($^{26}\text{Al}$, $^{60}\text{Fe}$) | dimensionless | `radiogenic_heating!` (`src/simulation/steps/radiogenic_heating.jl`) |

#### Stage 3: Marker Mapping & Gravity
Module: `src/simulation/steps/interpolate_markers_to_grid.jl`, `src/numerics/gravity.jl`, `src/physics/melting_convection.jl`, `src/config.jl`

| Symbol | Meaning | Unit | Where in the code |
|:---|:---|:---|:---|
| $\Phi$ | Gravitational potential | $\text{m}^2/\text{s}^2$ | `solve_gravity!` (`src/numerics/gravity.jl`) |
| $G$ | Newtonian gravitational constant | $\text{m}^3/(\text{kg}\cdot\text{s}^2)$ | `constants.jl` |
| $\rho$ | Bulk mixture density | $\text{kg/m}^3$ | `interpolate_markers_to_grid!` (`src/simulation/steps/interpolate_markers_to_grid.jl`) |
| $F_\mathrm{m}$ | Silicate melt fraction | dimensionless | `compute_melt_fraction` (`src/physics/melting_convection.jl`) |
| $T$ | Temperature | $\text{K}$ | `interpolate_markers_to_grid!` (`src/simulation/steps/interpolate_markers_to_grid.jl`) |
| $T_\mathrm{s}$ | Silicate solidus temperature | $\text{K}$ | `melting.T_solidus` (`src/config.jl`) |
| $T_\mathrm{l}$ | Silicate liquidus temperature | $\text{K}$ | `melting.T_liquidus` (`src/config.jl`) |
| $\sigma_\mathrm{SB}$ | Stefan-Boltzmann constant | $\text{W/}(\text{m}^2\cdot\text{K}^4)$ | `thermodynamics.sigma_sb` (`src/config.jl`) |

#### Stage 4: Hydration Reactions
Module: `src/physics/serpentinization.jl`, `src/constants.jl`

| Symbol | Meaning | Unit | Where in the code |
|:---|:---|:---|:---|
| $\mathrm{W_\mathrm{sil}}$ | Hydrous silicate mineral phase (serpentine) | dimensionless | `perform_thermochemical_reaction!` (`src/physics/serpentinization.jl`) |
| $\mathrm{D_\mathrm{sil}}$ | Anhydrous silicate mineral phase (olivine) | dimensionless | `perform_thermochemical_reaction!` (`src/physics/serpentinization.jl`) |
| $\mathrm{H_2O}$ | Pore water phase | dimensionless | `perform_thermochemical_reaction!` (`src/physics/serpentinization.jl`) |
| $\Delta t_\mathrm{r}$ | Thermochemical reaction relaxation timescale | $\text{s}$ | `compute_Δtreaction` (`src/physics/serpentinization.jl`) |
| $\Delta G$ | Gibbs free energy change of reaction | $\text{J/mol}$ | `compute_gibbs_free_energy` (`src/physics/serpentinization.jl`) |
| $\Delta H$ | Standard molar enthalpy change | $\text{J/mol}$ | `compute_gibbs_free_energy` (`src/physics/serpentinization.jl`) |
| $\Delta S$ | Standard molar entropy change | $\text{J/}(\text{mol}\cdot\text{K})$ | `compute_gibbs_free_energy` (`src/physics/serpentinization.jl`) |
| $\Delta V$ | Solid-fluid reaction molar volume change | $\text{m}^3/\text{mol}$ | `compute_gibbs_free_energy` (`src/physics/serpentinization.jl`) |
| $P_\mathrm{f}$ | Pore fluid pressure | $\text{Pa}$ | `perform_thermochemical_reaction!` (`src/physics/serpentinization.jl`) |
| $R_\mathrm{gas}$ | Universal gas constant | $\text{J/}(\text{mol}\cdot\text{K})$ | `constants.jl` |
| $X_\mathrm{D}, X_\mathrm{W}$ | Molar fractions of dry and wet silicate phases | dimensionless | `perform_thermochemical_reaction!` (`src/physics/serpentinization.jl`) |

#### Outer Iteration Loop (Stages 4–6)
Module: `src/simulation/loop.jl`, `src/config.jl`

| Symbol | Meaning | Unit | Where in the code |
|:---|:---|:---|:---|
| $n_\mathrm{iter}$ | Maximum plastic iterations per step | dimensionless | `solver.max_plastic_iterations` (`src/config.jl`) |
| $\text{iter}$ | Outer iteration loop counter | dimensionless | `simulation_loop` (`src/simulation/loop.jl`) |

#### Stage 5: Hydro-Mechanical Stokes-Darcy
Module: `src/numerics/stokes_darcy.jl`, `src/simulation/steps/solve_stokes_darcy.jl`

| Symbol | Meaning | Unit | Where in the code |
|:---|:---|:---|:---|
| $\mathbf{v}_\mathrm{s}$ | Solid rock matrix velocity vector | $\text{m/s}$ | `assemble_hydromechanical_lse!` (`src/numerics/stokes_darcy.jl`) |
| $P_\mathrm{t}$ | Total pressure (solid plus fluid) | $\text{Pa}$ | `assemble_hydromechanical_lse!` (`src/numerics/stokes_darcy.jl`) |
| $P_\mathrm{f}$ | Pore fluid pressure | $\text{Pa}$ | `assemble_hydromechanical_lse!` (`src/numerics/stokes_darcy.jl`) |
| $\eta_\phi$ | Bulk compaction viscosity | $\text{Pa}\cdot\text{s}$ | `assemble_hydromechanical_lse!` (`src/numerics/stokes_darcy.jl`) |
| $\phi$ | Porosity (pore volume fraction) | dimensionless | `assemble_hydromechanical_lse!` (`src/numerics/stokes_darcy.jl`) |
| $\beta_\mathrm{d}$ | Drained matrix compressibility | $\text{Pa}^{-1}$ | `assemble_hydromechanical_lse!` (`src/numerics/stokes_darcy.jl`) |
| $K_\mathrm{BW}$ | Biot-Willis poroelastic coefficient | dimensionless | `assemble_hydromechanical_lse!` (`src/numerics/stokes_darcy.jl`) |
| $\Delta t$ | Computational time step | $\text{s}$ | `assemble_hydromechanical_lse!` (`src/numerics/stokes_darcy.jl`) |
| $\Delta P_\mathrm{t}$ | Total pressure increment | $\text{Pa}$ | `assemble_hydromechanical_lse!` (`src/numerics/stokes_darcy.jl`) |
| $\Delta P_\mathrm{f}$ | Pore fluid pressure increment | $\text{Pa}$ | `assemble_hydromechanical_lse!` (`src/numerics/stokes_darcy.jl`) |
| $\Gamma_\mathrm{r}$ | Reaction volume source rate ($DMP$) | $\text{s}^{-1}$ | `assemble_hydromechanical_lse!` (`src/numerics/stokes_darcy.jl`) |
| $\mathbf{q}_\mathrm{D}$ | Darcy fluid discharge velocity vector | $\text{m/s}$ | `reconstruct_darcy_fluxes!` (`src/numerics/stokes_darcy.jl`) |
| $k$ | Dynamic matrix permeability | $\text{m}^2$ | `assemble_hydromechanical_lse!` (`src/numerics/stokes_darcy.jl`) |
| $\eta_\mathrm{f}$ | Pore fluid dynamic viscosity | $\text{Pa}\cdot\text{s}$ | `assemble_hydromechanical_lse!` (`src/numerics/stokes_darcy.jl`) |
| $\rho_\mathrm{f}$ | Pore fluid density | $\text{kg/m}^3$ | `assemble_hydromechanical_lse!` (`src/numerics/stokes_darcy.jl`) |
| $\mathbf{g}$ | Gravitational acceleration vector | $\text{m/s}^2$ | `assemble_hydromechanical_lse!` (`src/numerics/stokes_darcy.jl`) |

#### Stage 6: Thermal & Segregation
Module: `src/numerics/thermal.jl`, `src/simulation/steps/solve_thermal_energy.jl`, `src/numerics/segregation.jl`

| Symbol | Meaning | Unit | Where in the code |
|:---|:---|:---|:---|
| $\rho c_{p,\mathrm{app}}$ | Apparent volumetric heat capacity including melting | $\text{J/}(\text{m}^3\cdot\text{K})$ | `rhocp_apparent_silicate` (`src/physics/melting_convection.jl`) |
| $\mathrm{D}/\mathrm{D}t$ | Material derivative along marker trajectory | $\text{s}^{-1}$ | `solve_thermal_energy!` (`src/simulation/steps/solve_thermal_energy.jl`) |
| $T$ | Temperature | $\text{K}$ | `solve_thermal_energy!` (`src/simulation/steps/solve_thermal_energy.jl`) |
| $k_\mathrm{eff}$ | Effective thermal conductivity with convective scaling | $\text{W/}(\text{m}\cdot\text{K})$ | `regularized_soft_turbulence_conductivity` (`src/physics/melting_convection.jl`) |
| $H_\mathrm{r}$ | Volumetric radiogenic heating rate | $\text{W/m}^3$ | `assemble_thermal_lse!` (`src/numerics/thermal.jl`) |
| $H_\mathrm{s}$ | Volumetric shear heating dissipation | $\text{W/m}^3$ | `assemble_thermal_lse!` (`src/numerics/thermal.jl`) |
| $H_\mathrm{a}$ | Volumetric adiabatic work heating | $\text{W/m}^3$ | `assemble_thermal_lse!` (`src/numerics/thermal.jl`) |
| $H_\mathrm{seg}$ | Volumetric segregation dissipation heating | $\text{W/m}^3$ | `assemble_thermal_lse!` (`src/numerics/thermal.jl`) |
| $H_\mathrm{L}$ | Latent heat of reactions, crystallization, and venting | $\text{W/m}^3$ | `assemble_thermal_lse!` (`src/numerics/thermal.jl`) |
| $\mathbf{J}_\mathrm{Fe}$ | Metal segregation volume flux vector | $\text{m/s}$ | `apply_metal_segregation!` (`src/numerics/segregation.jl`) |
| $X_\mathrm{Fe}$ | Local metal volume fraction on markers | dimensionless | `apply_metal_segregation!` (`src/numerics/segregation.jl`) |
| $\phi_\mathrm{res}$ | Residual trapped metal volume fraction floor | dimensionless | `coreformation.phi_residual` (`src/config.jl`) |
| $\mathbf{v}_\mathrm{seg}$ | Metal segregation velocity | $\text{m/s}$ | `metal_segregation_velocity` (`src/numerics/segregation.jl`) |
| $v_\mathrm{perc}$ | Porous percolation metal segregation velocity | $\text{m/s}$ | `metal_segregation_velocity` (`src/numerics/segregation.jl`) |
| $k_\mathrm{Fe}$ | Permeability of molten metal through silicate matrix | $\text{m}^2$ | `metal_segregation_velocity` (`src/numerics/segregation.jl`) |
| $\eta_\mathrm{Fe}$ | Liquid metal dynamic viscosity | $\text{Pa}\cdot\text{s}$ | `metal_segregation_velocity` (`src/numerics/segregation.jl`) |
| $\Delta\rho$ | Metal-silicate density contrast ($\rho_\mathrm{Fe} - \rho_\mathrm{sil}$) | $\text{kg/m}^3$ | `metal_segregation_velocity` (`src/numerics/segregation.jl`) |
| $g$ | Magnitude of gravitational acceleration | $\text{m/s}^2$ | `metal_segregation_velocity` (`src/numerics/segregation.jl`) |

#### Stage 7: Marker Stress & Subgrid Diffusion
Module: `src/simulation/steps/diffuse_and_update_markers.jl`, `src/config.jl`

| Symbol | Meaning | Unit | Where in the code |
|:---|:---|:---|:---|
| $\Delta \sigma_m$ | Subgrid stress diffusion increment on marker $m$ | $\text{Pa}$ | `diffuse_and_update_markers!` (`src/simulation/steps/diffuse_and_update_markers.jl`) |
| $\sigma_\mathrm{grid}$ | Deviatoric stress interpolated from grid nodes | $\text{Pa}$ | `diffuse_and_update_markers!` (`src/simulation/steps/diffuse_and_update_markers.jl`) |
| $\sigma_m$ | Marker deviatoric stress state | $\text{Pa}$ | `diffuse_and_update_markers!` (`src/simulation/steps/diffuse_and_update_markers.jl`) |
| $d_s$ | Subgrid stress diffusion damping parameter | dimensionless | `solver.dsubgrids` (`src/config.jl`) |
| $\Delta t$ | Computational time step | $\text{s}$ | `diffuse_and_update_markers!` (`src/simulation/steps/diffuse_and_update_markers.jl`) |
| $\tau_\mathrm{M}$ | Maxwell viscoelastic relaxation timescale | $\text{s}$ | `diffuse_and_update_markers!` (`src/simulation/steps/diffuse_and_update_markers.jl`) |

#### Stage 8: Fluid Venting & Atmosphere
Module: `src/numerics/thermal.jl`, `src/atmosphere.jl`, `src/simulation/steps/vent_and_degas.jl`

| Symbol | Meaning | Unit | Where in the code |
|:---|:---|:---|:---|
| $T(\tau_\mathrm{IR})$ | Atmospheric temperature at infrared optical depth $\tau_\mathrm{IR}$ | $\text{K}$ | `compute_guillot_surface_temperature` (`src/atmosphere.jl`) |
| $\tau_\mathrm{IR}$ | Thermal infrared optical depth coordinate | dimensionless | `compute_guillot_surface_temperature` (`src/atmosphere.jl`) |
| $T_\mathrm{int}$ | Interior effective radiating temperature | $\text{K}$ | `compute_guillot_surface_temperature` (`src/atmosphere.jl`) |
| $T_\mathrm{eq}$ | Planetary equilibrium irradiation temperature | $\text{K}$ | `compute_guillot_surface_temperature` (`src/atmosphere.jl`) |
| $f_\gamma(\tau_\mathrm{IR})$ | Semi-grey atmospheric temperature shape factor | dimensionless | `compute_guillot_surface_temperature` (`src/atmosphere.jl`) |
| $S_\mathrm{vent}$ | Surface venting drainage rate (computed in Stage 5, applied in Stage 8) | $\text{s}^{-1}$ | `apply_venting_surface_boundary!` (`src/numerics/thermal.jl`) |
| $q_{D,\mathrm{face}}$ | Outward Darcy flux at venting boundary face | $\text{m/s}$ | `apply_venting_surface_boundary!` (`src/numerics/thermal.jl`) |
| $\Delta x$ | Grid cell spatial spacing | $\text{m}$ | `apply_venting_surface_boundary!` (`src/numerics/thermal.jl`) |
| $\phi$ | Local rock porosity | dimensionless | `apply_venting_surface_boundary!` (`src/numerics/thermal.jl`) |
| $\phi_\mathrm{min}$ | Minimum residual porosity floor | dimensionless | `poroelasticity.phimin` (`src/config.jl`) |
| $\Delta t$ | Computational time step | $\text{s}$ | `apply_venting_surface_boundary!` (`src/numerics/thermal.jl`) |

#### Stage 9: Marker Advection, Replenish & Domain Doubling
Module: `src/simulation/steps/advect_markers.jl`, `src/simulation/steps/replenish.jl`, `src/simulation/steps/telescope_domain.jl`

| Symbol | Meaning | Unit | Where in the code |
|:---|:---|:---|:---|
| $\mathbf{x}_m$ | Marker coordinate position vector in 2-D | $\text{m}$ | `advect_markers!` (`src/simulation/steps/advect_markers.jl`) |
| $t$ | Elapsed simulation physical time | $\text{s}$ | `advect_markers!` (`src/simulation/steps/advect_markers.jl`) |
| $\mathbf{v}(\mathbf{x}_m)$ | Velocity field interpolated to marker position | $\text{m/s}$ | `advect_markers!` (`src/simulation/steps/advect_markers.jl`) |
| $\text{RK4}$ | Fourth-order Runge-Kutta integration scheme | dimensionless | `advect_markers!` (`src/simulation/steps/advect_markers.jl`) |
| $r_\mathrm{max}$ | Maximum radial distance of active markers | $\text{m}$ | `telescope_domain!` (`src/simulation/steps/telescope_domain.jl`) |
| $R_\mathrm{tel}$ | Domain doubling threshold ($0.70 \cdot x_\mathrm{size}/2$) | $\text{m}$ | `telescoping.r_threshold_fraction` (`src/config.jl`) |

#### Stage 10: Diagnostics & Output
Module: `src/simulation/steps/diagnostics_and_output.jl`, `src/config.jl`

| Symbol | Meaning | Unit | Where in the code |
|:---|:---|:---|:---|
| $t$ | Elapsed simulation physical time | $\text{s}$ | `simulation_loop` (`src/simulation/loop.jl`) |
| $t_\mathrm{end}$ | Simulation termination physical time | $\text{s}$ | `time.endtime` (`src/config.jl`) |
| $\text{step}$ | Current simulation time step index | dimensionless | `simulation_loop` (`src/simulation/loop.jl`) |
| $n_\mathrm{steps}$ | Maximum allowed simulation time step count | dimensionless | `time.n_steps` (`src/config.jl`) |

### Standard Abbreviations

| Abbreviation | Meaning |
|:---|:---|
| **IW** | Iron-Wüstite oxygen fugacity buffer (equilibrium between metallic iron and wüstite FeO). |
| **QFM** | Quartz-Fayalite-Magnetite oxygen fugacity buffer (silicate-oxide redox reference). |
| **XUV** | Extreme-ultraviolet stellar irradiation driving photo-evaporative atmospheric escape. |
| **CAI** | Calcium-Aluminum-rich Inclusion (earliest dated Solar System refractory solids). |
| **RK4** | Fourth-order Runge-Kutta numerical integration scheme. |
| **UMFPACK** | Unsymmetric MultiFrontal sparse direct solver library used for production 2-D solves. |

---

## Thermo-Hydro-Mechanical Coupling

### How do mechanics, fluids, and thermodynamics exchange properties across stages?

The 10 operator stages couple through staggered grid fields and Lagrangian markers:

1. **Hydration to Mechanics**:
   Serpentinization reactions produce or consume pore fluid, injecting volume source rate $\Gamma_\mathrm{r} = DMP$ into the Stokes-Darcy fluid mass balance.
   Dehydration generates pore overpressure ($P_\mathrm{f} \to P_\mathrm{t}$), which triggers plastic yielding and dynamic hydrofracture.

2. **Mechanics to Heat**:
   Plastic shear deformation dissipates mechanical energy into shear heating $H_\mathrm{s} = \boldsymbol{\sigma}' : \dot{\boldsymbol{\epsilon}}_\mathrm{pl}$.
   Volumetric compaction generates adiabatic heating $H_\mathrm{a} = -P_\mathrm{t} (\nabla \cdot \mathbf{v}_\mathrm{s})$.

3. **Thermodynamics to Mechanics**:
   Temperature fields control solid rock viscosities via Arrhenius activation laws and mineral reaction rates via Gibbs free energy $\Delta G(T, P_\mathrm{f})$.
   Silicate melting reduces matrix viscosity by up to 10 orders of magnitude across the rheological transition $F_\mathrm{crit} \approx 0.40$.

4. **Fluids to Atmosphere**:
   Darcy fluid flux reaching the planetesimal surface vents pore water and volatiles into the atmospheric reservoir, supplying gases that govern greenhouse warming and surface pressure.

---

## Geometry Convention and Cylindrical-to-Spherical Mapping

### How does a 2-D cross-section represent a 3-D spherical planetesimal?

`Erebus.jl` solves conservation laws on a 2-D Cartesian grid $(x, y) \in [0, W_\mathrm{box}] \times [0, H_\mathrm{box}]$ with coordinates centered at $(x_c, y_c)$.
Radial distance from the center is:

$$r = \sqrt{(x - x_c)^2 + (y - y_c)^2}$$

For spherically symmetric distributions, mass increments and volatile transfers map between 2-D planar sums and 3-D spherical inventories via:

$$\Delta M = \sum_m \rho_m A_m w_\mathrm{3D}(r_m) \Delta C_m, \quad w_\mathrm{3D}(r_m) = 2 r_m, \quad A_m = \Delta x_m \cdot \Delta y_m$$

where $\rho_m$ is marker density, $A_m$ is marker differential area, and $\Delta C_m$ is the dimensionless mass fraction or phase change increment.

### Weighting Limits

The out-of-plane weighting function $w_\mathrm{3D}(r) = 2r$ exhibits 3 characteristic limits:
1. **Planetary surface limit ($r \to R$)**: Near the outer planetary radius, $w_\mathrm{3D}(R) = 2R$. Surface flux and atmospheric exchange scale with diameter $2R$.
2. **Finite shell mean ($r \in [r_\text{inner}, R]$)**: For an outer spherical shell, the ratio of integrated 3-D shell mass to 2-D planar shell mass equals $\frac{4}{3} \frac{R^3 - r_\text{inner}^3}{R^2 - r_\text{inner}^2}$. For a near-surface shell spanning $0.9R \le r \le R$, this ratio equals approximately $1.902 R$.
3. **Full spherical volume mean**: Integrating over the entire circular section yields full spherical volume $\frac{4}{3} \pi R^3$, corresponding to an area-weighted mean out-of-plane thickness of $\frac{4}{3} R$.

### Dynamic Flow and Conservation

Hydromechanical Stokes flow, Darcy filtration, and marker transport operate in the 2-D Cartesian cross-section.
Because 2-D Cartesian divergence-free velocity fields ($\nabla \cdot \mathbf{v} = 0$) do not preserve the axisymmetric volume element $r \,\mathrm{d}r$, planar advection does not automatically conserve the 3-D weighted inventory $\sum_m \rho_m A_m (2 r_m)$.
The helper functions `marker_out_of_plane_length` and `marker_area` provide geometric primitives to evaluate 3-D spherical projections alongside 2-D cross-sectional quantities.

Thermal conduction incorporates an optional radial metric term ($H_\mathrm{metric}$) to account for spherical divergence in conductive heat flow when `cfg.geometry.spherical_metric` is active.
It also applies spherical metric weighting ($w_\mathrm{3D}$) to marker segregation.

### Gravitational Potential

The 2-D Poisson equation uses the source term $\frac{8}{3}\pi G\rho$.
For a uniform sphere, this term yields exact gravity at every radius inside the body ($r \le R$).
Outside the body ($r > R$), 2-D cylindrical gravity drops as $1/r$ (potential $\propto \ln r$), not as $1/r^2$.
For differentiated bodies with dense metallic cores, select `geometry.gravity_mode = :enclosed_mass` (see [Telescoping Domain](telescoping_domain.md)).

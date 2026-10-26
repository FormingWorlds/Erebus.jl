# Scientific Context: Early Planetesimal Evolution

`Erebus.jl` models the internal physical and geochemical evolution of icy and rocky planetesimals during the first tens of millions of years of Solar System history.

---

## The Physical Scenario and Evolution Phases

### How does a planetesimal evolve from a cold porous body into a differentiated world?

Planetesimals evolve across 3 physical phases driven by radiogenic decay heating: Phase I (cold accretion), Phase II (hydrothermal alteration), and Phase III (magma ocean and core formation).

Figure 2 shows these physical phases and links them to the active processes in Erebus.

![Erebus.jl Planetesimal Evolution in Three Stages](../assets/erebus_planetesimal_stages.svg)

*Figure 2: Planetesimal evolution across three physical phases. Sectors show successive phases of 1 body; sector widths do not scale with time. Left and right cards show physical processes that Erebus solves. Active modules depend on configuration switches; all modules are shown active. Flux arrows are schematic. Schematic (hand-built diagram, not simulation output).*

In Figure 2:
- **Phase I (Cold Accretion)**: Low temperatures preserve porous ice and pristine rock. Accretion flux $\dot{M}_\mathrm{acc}$ delivers pebbles and volatile ices to the surface.
- **Phase II (Hydrothermal Alteration)**: Radiogenic decay warms the interior. Ice melts into pore water. Darcy percolation $\mathbf{q}_\mathrm{D}$ drives hydrothermal circulation, serpentinization, and dehydration reactions through permeable crust.
- **Phase III (Magma Ocean & Core Segregation)**: Deep heating drives silicate melting above the critical rheological transition ($F_\mathrm{m} > \phi_\mathrm{crit} \approx 0.40$). Liquid metal droplets segregate to form a metallic core, releasing segregation heating $H_\mathrm{seg}$. Magma ocean surface heat flux $F_\mathrm{MO}$ radiates into space and drives atmospheric escape.

### Physical Scenario Sequence

Planetesimals form in the protoplanetary disk through gravitational collapse of pebble swarms or streaming instabilities.
Bodies that accrete during the first 2 to 3 million years after Calcium-Aluminum-rich Inclusions (CAIs) incorporate short-lived radioactive isotopes: $^{26}\text{Al}$ (half-life $t_{1/2} \approx 0.717\text{ Ma}$) and $^{60}\text{Fe}$ ($t_{1/2} \approx 2.62\text{ Ma}$).
As radionuclides decay, volumetric heating drives 7 physical stages:

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
   When Terzaghi effective stress becomes tensile ($\\sigma_\text{eff} \le -\sigma_\text{tensile}$), the silicate matrix fractures.
   This creates high-permeability pathways that vent pore fluids to the planetesimal surface.

6. **Silicate Rock Melting and Magma Ocean (Phase III)**:
   In bodies with early accretion, decay heating raises temperatures past the silicate solidus ($T_\mathrm{s} \approx 1400\text{ K}$).
   Silicate rock melts into a crystal-liquid mixture.
   Above critical melt fraction $F_\mathrm{crit} \approx 0.40$ (`cfg.melting.phi_crit`), the mixture behaves as a liquid magma ocean.
   Soft turbulence convection transports heat outward to the surface, where the magma ocean radiates heat flux $F_\mathrm{MO}$ into space.

7. **Iron Core Formation and Metal Segregation (Phase III)**:
   Metal means liquid Fe-FeS alloy (iron plus iron sulfide).
   At temperatures above the Fe-FeS eutectic ($T_\text{eut} \approx 1213\text{ K}$), dense metallic liquid melts before silicate rock melts completely.
   In partially molten rock ($F_\mathrm{m} < F_\mathrm{perc,end} = 0.50$), metal trickles downward through solid silicates via porous percolation.
   In magma oceans ($F_\mathrm{m} \ge F_\mathrm{settle,start} = 0.40$), metal droplets rain downward via Stokes settling, releasing gravitational dissipation heating that accelerates core growth.

### Summary of Physical Processes by Evolution Phase

Table 2 lists key symbols in Figure 2 and links each process to code files.

| Phase | Symbol | Description | Unit | Code Location |
|:---|:---|:---|:---|:---|
| **Phase I** | $\dot{M}_\mathrm{acc}$ | Pebble and embryo accretion mass flux | $\text{kg/s}$ | `src/simulation/steps/accrete.jl` (`accrete!`) |
| | $T_\mathrm{amb}$ | Protoplanetary disk ambient temperature | $\text{K}$ | `src/simulation/steps/prepare_step.jl` (`prepare_step_ambient!`) |
| | $P_\mathrm{amb}$ | Protoplanetary disk ambient pressure | $\text{Pa}$ | `src/simulation/steps/prepare_step.jl` (`prepare_step_ambient!`) |
| **Phase II** | $\mathbf{q}_\mathrm{D}$ | Darcy pore fluid filtration velocity vector | $\text{m/s}$ | `src/numerics/stokes_darcy.jl` (`assemble_hydromechanical_lse!`) |
| | $\Delta G$ | Gibbs free energy change of serpentine dehydration | $\text{J/mol}$ | `src/physics/serpentinization.jl` (`compute_gibbs_free_energy`) |
| | $S_\mathrm{vent}$ | Surface venting volumetric drainage rate | $\text{s}^{-1}$ | `src/numerics/thermal.jl` (`apply_venting_surface_boundary!`) |
| | $\tau_\mathrm{IR}$ | Thermal infrared optical depth coordinate | dimensionless | `src/atmosphere.jl` (`compute_guillot_surface_temperature`) |
| **Phase III** | $F_\mathrm{m}$ | Degree of silicate partial melting | dimensionless | `src/physics/melting_convection.jl` (`compute_melt_fraction`) |
| | $\mathbf{J}_\mathrm{Fe}$ | Metal segregation volume flux vector | $\text{m/s}$ | `src/numerics/segregation.jl` (`apply_metal_segregation!`) |
| | $k_\mathrm{eff}$ | Effective thermal conductivity with soft turbulence | $\text{W/}(\text{m}\cdot\text{K})$ | `src/numerics/thermal.jl` (`assemble_thermal_lse!`) |
| | $F_\mathrm{MO}$ | Magma ocean radiative surface heat flux | $\text{W/m}^2$ | `src/simulation/steps/prepare_step.jl` (`prepare_step_ambient!`) |

---

## Global Architecture and Execution Pipeline

### What does the code do in one time step?

The simulation loop advances physical time by executing 10 operator-splitting stages per time step.
Each stage updates temperatures, mechanical flows, or geochemical states.
Figure 1 shows the 10 stages and their governing balance equations.

### Execution Pipeline and Governing Equations

![Erebus.jl Execution Pipeline and Governing Equations](../assets/erebus_execution_pipeline.svg)

*Figure 1: Execution pipeline and governing equations in Erebus. The center column shows the full run, where the boxes inside the frame form 1 time step. Side cards show governing equations and physical modules for each stage. The configuration file switches each module on or off; all modules are shown active. The outer iteration repeats stages 4, 5, and 6 until fields converge. If the plastic Stokes-Darcy iteration does not converge, the step restarts from stage 1 with dt = min(dt/2, dt_next); the run stops after solver.max_dt_reductions reductions. Schematic (hand-built diagram, not simulation output).*

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
   Free energy change $\Delta G$ is written for the forward dehydration direction.
   It updates mineral hydration states, fluid production, and reaction latent heat.

5. **Hydro-Mechanical Stokes-Darcy (`assemble_hydromechanical_lse!`, outer iteration)**:
   The solver assembles and solves the linear system for solid velocities, total pressure, fluid Darcy flux, and fluid pressure.
   It updates matrix compaction and plastic yielding.

6. **Thermal & Segregation (`solve_thermal_energy!`, `apply_metal_segregation!`, `apply_silicate_melt_segregation!`, outer iteration)**:
   The code solves the thermal energy equation with conduction and latent heat.
   Advective heat transport occurs through marker advection in Stage 9.
   Hydrothermal convection and soft-turbulence transport enter through effective thermal conductivity $k_\mathrm{eff}$.
   The code computes metal droplet settling and silicate melt ascent using conservative drift-flux transport.
   It adds heat $H_\mathrm{seg}$ released by metal and melt segregation.

7. **Marker Stress & Diffusion (`diffuse_and_update_markers!`)**:
   The code updates viscoelastic stress relaxation on markers after the outer iteration loop converges.
   It conducts subgrid stress diffusion between grid nodes and markers.

8. **Fluid Venting & Atmosphere (`vent_and_degas!`, `evolve_atmosphere!`)**:
   The code computes hydrofracture fluid venting and volatile exsolution from magma.
   It calls `advance_marker_thermo_porosity_venting!` to update marker temperatures and compaction porosity.
   It couples vented volatiles into surface atmospheric pressure and equilibrium temperature.
   It computes hydrodynamic and kinetic escape losses to space.

9. **Marker Advection & Doubling (`advect_markers!`, `replenish!`, `telescope_domain!`)**:
   The code advects markers with a 4th-order Runge-Kutta scheme (RK4).
   It updates nodal total and fluid pressures along marker paths.
   It replenishes markers in empty cells.
   The code doubles domain size when the body radius exceeds 70% of domain half-width (see [Telescoping Domain](telescoping_domain.md)).

10. **Diagnostics & Output (`advance_step_diagnostics!`)**:
    After telescoping, the loop advances physical time $t \leftarrow t + \Delta t$.
    The code logs conservation metrics, saves JLD2 checkpoints, and streams telemetry records.

### Outer Iteration and Timestep Control

Stages 4, 5, and 6 form an outer iteration loop that repeats up to `cfg.solver.max_plastic_iterations`:
- **Plastic convergence**: If plastic iterations do not converge, the step restores the snapshot, sets $\Delta t = \min(\Delta t/2, \Delta t_\mathrm{next})$, and restarts from Stage 1. The simulation stops after `cfg.solver.max_dt_reductions` reductions.
- **Thermochemical convergence**: The loop checks fluid pressure error against `cfg.solver.pferrmax`, iteration count $> 2$ (or zero reaction rate $DMP = 0$), and temperature increment $\max(\Delta T) \le \Delta T_\mathrm{max}$. If thermochemical iterations do not converge within `max_plastic_iterations`, the step continues with `thermochemical_converged = false`.
- **Accretion timestep resync**: If `cfg.accretion.active` is true and $\Delta t$ was reduced during the step, the code restores the snapshot and restarts the step once at the reduced $\Delta t$ to preserve accretion rate synchronization.

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

The tables below list every mathematical symbol in Figure 1, grouped by execution stage.

#### Stage 1: Ambient & Accretion
Module: `src/simulation/steps/prepare_step.jl`, `src/simulation/steps/accrete.jl`

| Symbol | Description | Unit | Code Variable |
|:---|:---|:---|:---|
| $\dot{M}_\mathrm{acc}$ | Pebble and embryo accretion mass rate | $\text{kg/s}$ | `accrete!` |
| $T_\mathrm{amb}$ | Protoplanetary disk ambient gas temperature | $\text{K}$ | `prepare_step_ambient!` |
| $P_\mathrm{amb}$ | Protoplanetary disk ambient gas pressure | $\text{Pa}$ | `prepare_step_ambient!` |
| $R(t)$ | Planetesimal physical outer radius | $\text{m}$ | `state.accumulators.rplanet` |

#### Stage 2: Radiogenic Heating
Module: `src/physics/radioactivity.jl`, `src/simulation/steps/radiogenic_heating.jl`

| Symbol | Description | Unit | Code Variable |
|:---|:---|:---|:---|
| $H_\mathrm{r}$ | Volumetric radiogenic heating rate | $\text{W/m}^3$ | `radiogenic_heating!` |
| $Q_i$ | Specific decay heating for isotope $i$ ($^{26}\text{Al}$, $^{60}\text{Fe}$) | $\text{W/kg}$ | `Q_radiogenic` |
| $f_i$ | Mass fraction of element host in rock | dimensionless | `cfg.radioactivity.f_al`, `cfg.radioactivity.f_fe` |
| $r_i$ | Initial isotopic abundance ratio ($^{26}\text{Al}/^{27}\text{Al}$, $^{60}\text{Fe}/^{56}\text{Fe}$) | dimensionless | `cfg.radioactivity.r_al`, `cfg.radioactivity.r_fe` |
| $E_i$ | Decay energy per decay event | $\text{J}$ | `cfg.radioactivity.E_al`, `cfg.radioactivity.E_fe` |
| $\tau_i$ | Mean isotope lifetime ($\tau_i = t_{1/2} / \ln 2$) | $\text{s}$ | `radioactivity.jl` |
| $t$ | Elapsed time since CAI formation | $\text{s}$ | `state.timesum` |
| $X_\mathrm{ref}$ | Reference chondritic silicate or metal fraction | dimensionless | `cfg.coreformation.Xfe_bulk` |

#### Stage 3: Marker Mapping & Gravity
Module: `src/simulation/steps/interpolate_markers_to_grid.jl`, `src/numerics/gravity.jl`, `src/physics/melting_convection.jl`, `src/config.jl`

| Symbol | Description | Unit | Code Variable |
|:---|:---|:---|:---|
| $\Phi$ | Gravitational potential | $\text{m}^2/\text{s}^2$ | `ws.F_grav` |
| $G$ | Newtonian gravitational constant | $\text{m}^3/(\text{kg}\cdot\text{s}^2)$ | `constants.jl` (`G_GRAV`) |
| $\rho$ | Bulk mixture density | $\text{kg/m}^3$ | `state.RHO` |
| $\mathbf{g}$ | Gravitational acceleration vector ($\mathbf{g} = -\nabla\Phi$) | $\text{m/s}^2$ | `state.gx`, `state.gy` |
| $F_\mathrm{grav}$ | Radial gravity field diagnostic label | $\text{m/s}^2$ | `ws.F_grav` |
| $F_\mathrm{m}$ | Silicate melt fraction | dimensionless | `compute_melt_fraction` |
| $T_\mathrm{s}$ | Silicate solidus temperature | $\text{K}$ | `cfg.melting.T_solidus` |
| $T_\mathrm{l}$ | Silicate liquidus temperature | $\text{K}$ | `cfg.melting.T_liquidus` |
| $\sigma_\mathrm{SB}$ | Stefan-Boltzmann constant | $\text{W/}(\text{m}^2\cdot\text{K}^4)$ | `cfg.thermodynamics.sigma_sb` |

#### Stage 4: Hydration Reactions
Module: `src/physics/serpentinization.jl`, `src/constants.jl`

| Symbol | Description | Unit | Code Variable |
|:---|:---|:---|:---|
| $\mathrm{W_\mathrm{sil}}$ | Hydrous silicate mineral phase (serpentine) | dimensionless | `serpentinization.jl` |
| $\mathrm{D_\mathrm{sil}}$ | Anhydrous silicate mineral phase (olivine) | dimensionless | `serpentinization.jl` |
| $\mathrm{H_2O}$ | Pore water phase | dimensionless | `serpentinization.jl` |
| $\Delta t_\mathrm{r}$ | Thermochemical reaction relaxation timescale | $\text{s}$ | `constants.jl` (`Δtreaction`) |
| $\Delta G$ | Gibbs free energy change of reaction | $\text{J/mol}$ | `compute_gibbs_free_energy` |
| $\Delta H$ | Reaction standard enthalpy change | $\text{J/mol}$ | `serpentinization.jl` |
| $\Delta S$ | Reaction standard entropy change | $\text{J/}(\text{mol}\cdot\text{K})$ | `serpentinization.jl` |
| $\Delta V$ | Solid-fluid reaction molar volume change | $\text{m}^3/\text{mol}$ | `serpentinization.jl` |
| $P_\mathrm{f}$ | Pore fluid pressure | $\text{Pa}$ | `state.pf` |
| $R_\mathrm{gas}$ | Universal gas constant | $\text{J/}(\text{mol}\cdot\text{K})$ | `constants.jl` (`R_GAS`) |
| $X_\mathrm{D}, X_\mathrm{W}$ | Molar fractions of dry and wet silicate phases | dimensionless | `serpentinization.jl` |

#### Stage 5: Hydro-Mechanical Stokes-Darcy
Module: `src/numerics/stokes_darcy.jl`, `src/simulation/steps/solve_stokes_darcy.jl`

| Symbol | Description | Unit | Code Variable |
|:---|:---|:---|:---|
| $\mathbf{v}_\mathrm{s}$ | Solid rock matrix velocity vector | $\text{m/s}$ | `state.vx`, `state.vy` |
| $P_\mathrm{t}$ | Total pressure (solid plus fluid) | $\text{Pa}$ | `state.pr` |
| $P_\mathrm{f}$ | Pore fluid pressure | $\text{Pa}$ | `state.pf` |
| $\eta_\phi$ | Porosity-dependent bulk compaction viscosity | $\text{Pa}\cdot\text{s}$ | `stokes_darcy.jl` |
| $\phi$ | Porosity (pore volume fraction) | dimensionless | `state.PHIX`, `state.PHIY` |
| $\beta_\mathrm{d}$ | Dynamic drained matrix compressibility | $\text{Pa}^{-1}$ | `stokes_darcy.jl` |
| $K_\mathrm{BW}$ | Biot-Willis poroelastic coefficient | dimensionless | `stokes_darcy.jl` |
| $\Gamma_\mathrm{r}$ | Reaction volume source rate ($DMP$) | $\text{s}^{-1}$ | `stokes_darcy.jl` (`DMP`) |
| $\Delta t$ | Computational time step | $\text{s}$ | `state.dt` |
| $\mathbf{q}_\mathrm{D}$ | Darcy fluid discharge velocity vector | $\text{m/s}$ | `state.qxD`, `state.qyD` |
| $k$ | Dynamic matrix permeability | $\text{m}^2$ | `state.KX`, `state.KY` |
| $\eta_\mathrm{f}$ | Pore fluid dynamic viscosity | $\text{Pa}\cdot\text{s}$ | `state.eta_f` |
| $\rho_\mathrm{f}$ | Pore fluid density | $\text{kg/m}^3$ | `state.rho_f` |
| $\mathbf{g}$ | Gravitational acceleration vector | $\text{m/s}^2$ | `state.gx`, `state.gy` |
| $\boldsymbol{\sigma}'$ | Deviatoric stress tensor | $\text{Pa}$ | `state.SXX`, `state.SXY` |

#### Stage 6: Thermal & Segregation
Module: `src/numerics/thermal.jl`, `src/simulation/steps/solve_thermal_energy.jl`, `src/numerics/segregation.jl`

| Symbol | Description | Unit | Code Variable |
|:---|:---|:---|:---|
| $\rho c_{p,\mathrm{app}}$ | Apparent volumetric heat capacity including melting | $\text{J/}(\text{m}^3\cdot\text{K})$ | `rhocp_apparent_silicate` |
| $\mathrm{D}/\mathrm{D}t$ | Material derivative along marker trajectory | $\text{s}^{-1}$ | `assemble_thermal_lse!` |
| $k_\mathrm{eff}$ | Effective thermal conductivity with convective scaling | $\text{W/}(\text{m}\cdot\text{K})$ | `assemble_thermal_lse!` |
| $H_\mathrm{s}$ | Volumetric shear heating dissipation | $\text{W/m}^3$ | `assemble_thermal_lse!` |
| $H_\mathrm{a}$ | Volumetric adiabatic work heating | $\text{W/m}^3$ | `assemble_thermal_lse!` |
| $H_\mathrm{seg}$ | Volumetric segregation dissipation heating | $\text{W/m}^3$ | `assemble_thermal_lse!` |
| $H_\mathrm{L}$ | Reaction enthalpy, melting latent heat, and venting cooling | $\text{W/m}^3$ | `assemble_thermal_lse!` |
| $H_\mathrm{metric}$ | Optional spherical metric conduction term | $\text{W/m}^3$ | `assemble_thermal_lse!` |
| $\mathbf{J}_\mathrm{Fe}$ | Metal segregation volume flux vector | $\text{m/s}$ | `apply_metal_segregation!` |
| $k_\mathrm{Fe}$ | Permeability of molten metal through silicate matrix | $\text{m}^2$ | `apply_metal_segregation!` |
| $\eta_\mathrm{Fe}$ | Liquid metal dynamic viscosity | $\text{Pa}\cdot\text{s}$ | `apply_metal_segregation!` |
| $\mathbf{v}_\mathrm{seg}$ | Metal segregation velocity (percolation or settling) | $\text{m/s}$ | `metal_segregation_velocity` |
| $\Delta\rho$ | Metal-silicate density contrast ($\rho_\mathrm{Fe} - \rho_\mathrm{sil}$) | $\text{kg/m}^3$ | `segregation.jl` |

#### Stage 7: Marker Stress & Diffusion
Module: `src/simulation/steps/diffuse_and_update_markers.jl`, `src/config.jl`

| Symbol | Description | Unit | Code Variable |
|:---|:---|:---|:---|
| $\Delta \sigma_m$ | Subgrid stress diffusion increment on marker $m$ | $\text{Pa}$ | `diffuse_and_update_markers!` |
| $\sigma_\mathrm{grid}$ | Interpolated deviatoric stress from grid nodes | $\text{Pa}$ | `diffuse_and_update_markers!` |
| $\sigma_m$ | Marker deviatoric stress state | $\text{Pa}$ | `diffuse_and_update_markers!` |
| $d_s$ | Subgrid stress diffusion damping parameter | dimensionless | `cfg.solver.dsubgrids` |
| $\tau_\mathrm{M}$ | Maxwell viscoelastic relaxation timescale | $\text{s}$ | `diffuse_and_update_markers!` |

#### Stage 8: Fluid Venting & Atmosphere
Module: `src/numerics/thermal.jl`, `src/atmosphere.jl`, `src/simulation/steps/vent_and_degas.jl`

| Symbol | Description | Unit | Code Variable |
|:---|:---|:---|:---|
| $T(\tau_\mathrm{IR})$ | Atmospheric temperature at infrared optical depth $\tau_\mathrm{IR}$ | $\text{K}$ | `compute_guillot_surface_temperature` |
| $\tau_\mathrm{IR}$ | Thermal infrared optical depth coordinate | dimensionless | `compute_guillot_surface_temperature` |
| $T_\mathrm{int}$ | Interior effective radiating temperature | $\text{K}$ | `compute_guillot_surface_temperature` |
| $T_\mathrm{eq}$ | Planetary equilibrium irradiation temperature | $\text{K}$ | `compute_guillot_surface_temperature` |
| $f_\gamma(\tau_\mathrm{IR})$ | Semi-grey atmospheric temperature profile shape factor | dimensionless | `compute_guillot_surface_temperature` |
| $\gamma$ | Visible-to-thermal opacity ratio | dimensionless | `compute_guillot_surface_temperature` |
| $S_\mathrm{vent}$ | Surface venting volumetric drainage rate | $\text{s}^{-1}$ | `apply_venting_surface_boundary!` |
| $\phi_\mathrm{min}$ | Minimum residual porosity floor | dimensionless | `cfg.poroelasticity.phimin` |

#### Stage 9: Marker Advection & Doubling
Module: `src/simulation/steps/advect_markers.jl`, `src/simulation/steps/replenish.jl`, `src/simulation/steps/telescope_domain.jl`

| Symbol | Description | Unit | Code Variable |
|:---|:---|:---|:---|
| $\mathbf{x}_m$ | Marker position vector in 2-D coordinates | $\text{m}$ | `state.markers.xm`, `state.markers.ym` |
| $\mathbf{v}(\mathbf{x}_m)$ | Solid velocity interpolated to marker position | $\text{m/s}$ | `advect_markers!` |
| $\text{RK4}$ | Fourth-order Runge-Kutta time integration scheme | dimensionless | `advect_markers!` |
| $R_\mathrm{tel}$ | Domain doubling telescoping threshold ($0.70 W_\mathrm{box}$) | $\text{m}$ | `telescope_domain!` |

#### Stage 10: Diagnostics & Output
Module: `src/simulation/steps/diagnostics_and_output.jl`, `src/config.jl`

| Symbol | Description | Unit | Code Variable |
|:---|:---|:---|:---|
| $n_\mathrm{iter}$ | Maximum plastic iterations per step | dimensionless | `cfg.solver.max_plastic_iterations` |
| $t_\mathrm{end}$ | Simulation termination physical time | $\text{yr (config) / s (internal)}$ | `cfg.time.endtime` |
| $n_\mathrm{steps}$ | Maximum allowed simulation time steps | dimensionless | `cfg.time.n_steps` |

### Standard Abbreviations

| Abbreviation | Full Name and Physical Meaning |
|:---|:---|
| **IW** | Iron-Wüstite oxygen fugacity buffer (equilibrium between metallic iron and wüstite FeO). |
| **QFM** | Quartz-Fayalite-Magnetite oxygen fugacity buffer (silicate-oxide redox reference). |
| **XUV** | Extreme-ultraviolet stellar irradiation driving photo-evaporative atmospheric escape. |
| **CAI** | Calcium-Aluminum-rich Inclusion (earliest dated Solar System refractory solids). |
| **THM** | Thermo-Hydro-Mechanical (coupled thermal, porous fluid, and matrix deformation equations). |
| **LSE** | Linear System of Equations (assembled sparse matrix problem). |
| **RK4** | Fourth-order Runge-Kutta numerical integration scheme. |
| **UMFPACK** | Unsymmetric MultiFrontal sparse direct solver library used for production 2-D solves. |

---

## Thermo-Hydro-Mechanical Coupling

### Which solvers share state in each time step?

The simulation couples 5 numerical modules in a unified sequence:
- The thermal solver calculates conduction, radiogenic decay, latent heat buffering, and segregation heating. Solid advective heat transport occurs via marker movement. Hydrothermal heat transport enters through effective thermal conductivity $k_\mathrm{eff}$.
- The hydromechanical solver simultaneously solves Stokes solid matrix deformation, Darcy fluid flux, and poroelastic volume changes in a monolithic linear system.
- Marker-in-Cell advects material phases, temperature, composition, and porosity without numerical diffusion across moving boundaries.
- Melting and soft turbulence routines evaluate pressure-dependent silicate melting, apparent heat capacity latent heat buffering, suspension rheology, and convective conductivity.
- Magma and metal drift-flux solvers track conservative migration of silicate melt and molten iron droplets, computing gravitational dissipation heating and mantle depletion.

---

## Geometry Convention

### How does a 2-D cross-section represent a 3-D sphere?

`Erebus.jl` models the planetesimal as a 2-D Cartesian cross-section passing through planetary center $(x_\text{center}, y_\text{center})$.

To evaluate 3-D integral quantities (mass, component inventories, and volatile degassing) from planar marker positions, each marker $m$ at distance $r_m = \sqrt{(x_m - x_\text{center})^2 + (y_m - y_\text{center})^2}$ carries an effective out-of-plane spherical integration length $w_\mathrm{3D}(r_m)$:

$$w_\mathrm{3D}(r_m) = 2 r_m$$

Under spherical symmetry, each concentric ring at radius $r_m$ represents a spherical shell of radius $r_m$.
Integrating over a circular disk of radius $R$ with differential marker area $\mathrm{d}A = \mathrm{d}x \,\mathrm{d}y$ recovers the exact volume of a sphere:

$$V_\text{3D} = \int_\text{disk} w_\mathrm{3D}(r) \,\mathrm{d}A = \int_0^{2\pi} \mathrm{d}\theta \int_0^R (2r) \, r \,\mathrm{d}r = 4\pi \int_0^R r^2 \,\mathrm{d}r = \frac{4}{3} \pi R^3$$

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

### Gravitational Potential

The 2-D Poisson equation uses the source term $\frac{8}{3}\pi G\rho$.
For a uniform sphere, this term yields exact gravity at every radius inside the body ($r \le R$).
Outside the body ($r > R$), 2-D cylindrical potential drops as $1/r$ rather than $1/r^2$.
For differentiated bodies with dense metallic cores, select `geometry.gravity_mode = :enclosed_mass` (see [Telescoping Domain](telescoping_domain.md)).

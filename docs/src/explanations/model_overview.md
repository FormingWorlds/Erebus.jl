# Scientific Context: Early Planetesimal Evolution

`Erebus.jl` models the internal physical and geochemical evolution of icy and rocky planetesimals during the first tens of millions of years of Solar System history.

---

## The Physical Scenario

Planetesimals formed in the solar protoplanetary disk through the gravitational collapse of pebble swarms or streaming instabilities.
Bodies that accreted during the first 2 to 3 million years after the formation of Calcium-Aluminum-rich Inclusions (CAIs) incorporated substantial quantities of short-lived radioactive isotopes, predominantly $^{26}\text{Al}$ (half-life $t_{1/2} \approx 0.717\text{ Ma}$) and $^{60}\text{Fe}$ ($t_{1/2} \approx 2.62\text{ Ma}$).

As these radionuclides decayed, they generated intense internal volumetric heating:

### 1. Ice Melting
Initial temperatures of $\approx 170\text{ K}$ rose rapidly.
Once internal temperatures reached $273\text{ K}$, primordial water ice melted, absorbing latent heat ($L^f \approx 333\text{ kJ/kg}$) and generating mobile pore water in the porous silicate matrix.

### 2. Hydrothermal Circulation
Liquid water percolated through the interconnected pore network under Darcy flow, driven by thermal buoyancy and compaction-induced fluid pressure gradients.

### 3. Clay Dehydration and Gas Generation
At elevated temperatures ($T \approx 500\text{ to } 900\text{ K}$), hydrous phyllosilicates, such as serpentine, underwent thermal dehydration reactions, releasing bound hydroxyl groups as free supercritical fluid.

### 4. Matrix Compaction and Pore Overpressure
Viscous, elastic, and plastic deformation of the silicate matrix reduced porosity ($\phi$).
When the compaction rate exceeded the rate of fluid escape via Darcy percolation, pore fluid pressure ($P_\mathrm{f}$) rose to match or exceed the total solid pressure ($P_\mathrm{t}$).

### 5. Failure and Hydrofracture
When the Terzaghi effective stress became tensile ($\sigma_{\text{eff}} \le -\sigma_{\text{tensile}}$), the silicate matrix experienced hydrofracturing, creating high-permeability pathways that vented fluids to the planetesimal surface.

### 6. Silicate Rock Melting, Magma Ascent, and Magma Ocean Formation
In planetesimals that accreted early, decay heating drives temperatures past the silicate solidus ($T_\mathrm{s} \approx 1400\text{ K}$).
Above the solidus, silicate rock partially melts into a crystal-liquid mush.
Above the critical melt fraction ($F_\mathrm{crit} \approx 0.40$, configured via `melting.phi_crit`), the mush behaves as a liquid magma ocean.
Liquid Fe-FeS metal trickles downward by porous percolation before silicate rock melts (`coreformation.percolation_active`).
Buoyant silicate melt segregates outward via porous Darcy percolation and hindered Stokes crystal settling.
This ascent releases gravitational shear dissipation and latent heat of crystallization.
Sub-grid soft turbulence convection transports heat outward to the surface, buffering interior temperatures and governing core segregation.

### 7. Iron Core Formation and Metal Segregation
At temperatures exceeding the Fe-FeS eutectic ($T_{\text{eut}} \approx 1213\text{ K}$), dense metallic liquid melts and segregates inward toward the center.
Metal means liquid Fe-FeS alloy (iron plus iron sulfide).
In partially molten rock, metal trickles downward by porous Darcy percolation.
In magma oceans, metal droplets rain downward via Stokes settling, releasing gravitational potential energy that accelerates core formation.

---

## Global Architecture & Execution Flow

### What does the code do in one time step?
The simulation loop runs 10 stages in each time step.
Each stage updates temperatures, mechanical flows, or geochemical states.
Stages 4, 5, and 6 repeat inside an outer iteration until mechanical and thermal fields converge.
Figure 1 shows the 10 stages and their governing balance equations.

### Execution Pipeline and Governing Equations

![Erebus.jl Execution Pipeline and Governing Equations](../assets/erebus_execution_pipeline.svg)

*Figure 1: Execution pipeline and governing equations in Erebus. The center column shows the full run, where the boxes inside the frame form 1 time step. Side cards show governing equations and physical modules for each stage. The config file switches each module on or off; all modules are shown active. The outer iteration repeats stages 4, 5, and 6 until they agree. If the plastic iteration does not converge, the code restores the step snapshot and halves $\Delta t$. Thermal non-convergence reduces $\Delta t$ inside the iteration and does not restart the step.*

Each marker carries the out-of-plane spherical weight $w_\mathrm{3D} = 2r$ (`src/geometry.jl`, `marker_out_of_plane_length`).
Integrating over a disk of radius $R$ gives the sphere volume $\frac{4}{3}\pi R^3$.
The 2-D Poisson equation uses the source term $\frac{8}{3}\pi G\rho$.
For a uniform sphere, this term yields exact gravity at every radius.
Differentiated bodies can select `geometry.gravity_mode = :enclosed_mass` (see [Telescoping Domain](telescoping_domain.md)).

Volumetric radiogenic heating $H_\mathrm{r}$ scales radionuclide decay powers by phase volume fractions and phase densities (`src/simulation/steps/radiogenic_heating.jl`, `radiogenic_heating!`):
$$H_\mathrm{r} = (1 - \phi) \left[ (1 - \phi_\mathrm{Fe})\,\rho_\mathrm{sil}\,\frac{Q_{^{26}\mathrm{Al}}}{1 - X_\mathrm{ref}} + \phi_\mathrm{Fe}\,\rho_\mathrm{Fe}\,\frac{Q_{^{60}\mathrm{Fe}}}{X_\mathrm{ref}} \right]$$
Without metal markers, the code sets $X_\mathrm{ref} = 0.182$ (`src/constants.jl`, `X_FE_REF_CHONDRITE`) and computes $\phi_\mathrm{Fe}$ from silicate and metal densities.
With metal markers active, $X_\mathrm{ref}$ uses `cfg.coreformation.Xfe_bulk` and $\phi_\mathrm{Fe}$ is the marker metal volume fraction `Xfe_bulk[m]`.

Table 1 defines the mathematical symbols in Figure 1.

*Table 1: Symbols in Figure 1 equation cards.*

| Card / Stage | Symbol | Meaning | Unit | Code location |
| :--- | :--- | :--- | :--- | :--- |
| **1. Ambient & Accretion** | $w_\mathrm{3D}$ | Out-of-plane spherical weight per marker | $\text{m}$ | `src/geometry.jl` (`marker_out_of_plane_length`) |
| | $r$ | Radial distance from planetesimal center | $\text{m}$ | `src/geometry.jl` (`marker_out_of_plane_length`) |
| | $A$ | Cross-sectional area in 2-D planar slice | $\text{m}^2$ | `src/simulation/steps/accrete.jl` (`accrete!`) |
| | $R$ | Planetesimal outer radius | $\text{m}$ | `src/simulation/steps/accrete.jl` (`accrete!`) |
| **2. Radiogenic Heating** | $Q_i$ | Specific radiogenic heat production rate per kg of bulk rock | $\text{W/kg}$ | `src/physics/radioactivity.jl` (`Q_radiogenic`) |
| | $f_i$ | Atoms of stable reference isotope ($^{27}\text{Al}$, $^{56}\text{Fe}$) per kg of bulk rock | $\text{atoms/kg}$ | `src/physics/radioactivity.jl` (`Q_radiogenic`) |
| | $r_i$ | Initial isotopic abundance ratio at CAI formation | 1 | `src/physics/radioactivity.jl` (`Q_radiogenic`) |
| | $E_i$ | Nuclear decay energy | $\text{J}$ | `src/physics/radioactivity.jl` (`Q_radiogenic`) |
| | $\tau_{1/2, i}$ | Radionuclide half-life | $\text{s}$ | `src/physics/radioactivity.jl` (`Q_radiogenic`) |
| | $t$ | Time since CAI formation (tracked by `timesum`, starts at `time.start_time`) | $\text{s}$ | `src/physics/radioactivity.jl` (`Q_radiogenic`) |
| | $H_\mathrm{r}$ | Volumetric radiogenic heating rate | $\text{W/m}^3$ | `src/simulation/steps/radiogenic_heating.jl` (`radiogenic_heating!`) |
| | $\phi$ | Matrix porosity | 1 | `src/simulation/steps/radiogenic_heating.jl` (`radiogenic_heating!`) |
| | $\phi_\mathrm{Fe}$ | Metal volume fraction in solid matrix | 1 | `src/simulation/steps/radiogenic_heating.jl` (`radiogenic_heating!`) |
| | $X_\mathrm{ref}$ | Reference Fe mass fraction (`X_FE_REF_CHONDRITE` or `Xfe_bulk`) | 1 | `src/constants.jl`, `src/config.jl` |
| | $\rho_\mathrm{sil}$ | Silicate grain density | $\text{kg/m}^3$ | `src/physics/radioactivity.jl` (`calculate_radioactive_heating`) |
| | $\rho_\mathrm{Fe}$ | Iron metal density | $\text{kg/m}^3$ | `src/physics/radioactivity.jl` (`calculate_radioactive_heating`) |
| **3. Marker Mapping & Gravity** | $\Phi$ | Gravitational potential | $\text{m}^2/\text{s}^2$ | `src/simulation/steps/solve_gravity.jl` (`solve_gravity!`) |
| | $G$ | Gravitational constant | $\text{m}^3/(\text{kg}\cdot\text{s}^2)$ | `src/constants.jl` (`G_GRAV`) |
| | $\rho$ | Total bulk mass density | $\text{kg/m}^3$ | `src/simulation/steps/solve_gravity.jl` (`solve_gravity!`) |
| | $F_\mathrm{m}$ | Equilibrium silicate melt mass fraction | 1 | `src/physics/melting_convection.jl` (`compute_melt_fraction`) |
| | $T_\mathrm{s}, T_\mathrm{l}$ | Solidus and liquidus temperatures | $\text{K}$ | `src/physics/melting_convection.jl` (`compute_melt_fraction`) |
| | $\sigma_\mathrm{SB}$ | Stefan-Boltzmann constant | $\text{W/}(\text{m}^2\cdot\text{K}^4)$ | `src/constants.jl` (`SIGMA_SB`) |
| **4. Hydration Reactions** | $\mathrm{W_{sil}}$ | Hydrous silicate phase (serpentine) | 1 | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| | $\mathrm{D_{sil}}$ | Dry silicate phase (olivine) | 1 | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| | $\mathrm{H_2O}$ | Pore fluid water | 1 | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| | $\Delta G$ | Molar Gibbs free energy change of dehydration (forward direction) | $\text{J/mol}$ | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| | $\Delta H$ | Standard molar reaction enthalpy | $\text{J/mol}$ | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| | $T$ | Temperature | $\text{K}$ | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| | $\Delta S$ | Standard molar reaction entropy | $\text{J/(mol}\cdot\text{K)}$ | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| | $P_\mathrm{f}$ | Pore fluid pressure | $\text{Pa}$ | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| | $\Delta V$ | Standard molar volume change | $\text{m}^3/\text{mol}$ | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| | $R_\mathrm{gas}$ | Universal gas constant | $\text{J/(mol}\cdot\text{K)}$ | `src/constants.jl` (`R_GAS`) |
| | $X_\mathrm{D}, X_\mathrm{W}$ | Molar fractions of dry and hydrous silicate | 1 | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| | $\Delta t_\mathrm{r}$ | Reaction kinetic relaxation timescale | $\text{s}$ | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| **5. Hydro-Mechanical Stokes-Darcy** | $\mathbf{v}_\mathrm{s}$ | Solid matrix velocity vector | $\text{m/s}$ | `src/numerics/stokes_darcy.jl` (`assemble_hydromechanical_lse!`) |
| | $P_\mathrm{t}$ | Total solid pressure | $\text{Pa}$ | `src/numerics/stokes_darcy.jl` (`assemble_hydromechanical_lse!`) |
| | $\eta_\phi$ | Bulk compaction viscosity | $\text{Pa}\cdot\text{s}$ | `src/numerics/stokes_darcy.jl` (`assemble_hydromechanical_lse!`) |
| | $\beta_\mathrm{d}$ | Drained matrix compressibility | $\text{Pa}^{-1}$ | `src/numerics/stokes_darcy.jl` (`assemble_hydromechanical_lse!`) |
| | $K_\mathrm{BW}$ | Biot-Willis poroelastic coefficient | 1 | `src/numerics/stokes_darcy.jl` (`assemble_hydromechanical_lse!`) |
| | $B$ | Skempton pore pressure coefficient | 1 | `src/numerics/stokes_darcy.jl` (`assemble_hydromechanical_lse!`) |
| | $\Delta t$ | Computational time step | $\text{s}$ | `src/numerics/stokes_darcy.jl` (`assemble_hydromechanical_lse!`) |
| | $\mathbf{q}_\mathrm{D}$ | Darcy fluid discharge velocity | $\text{m/s}$ | `src/numerics/stokes_darcy.jl` (`assemble_hydromechanical_lse!`) |
| | $k$ | Dynamic matrix permeability | $\text{m}^2$ | `src/numerics/stokes_darcy.jl` (`assemble_hydromechanical_lse!`) |
| | $\eta_\mathrm{f}$ | Pore fluid dynamic viscosity | $\text{Pa}\cdot\text{s}$ | `src/numerics/stokes_darcy.jl` (`assemble_hydromechanical_lse!`) |
| | $\rho_\mathrm{f}$ | Pore fluid density | $\text{kg/m}^3$ | `src/numerics/stokes_darcy.jl` (`assemble_hydromechanical_lse!`) |
| | $\mathbf{g}$ | Gravitational acceleration vector | $\text{m/s}^2$ | `src/numerics/stokes_darcy.jl` (`assemble_hydromechanical_lse!`) |
| | $\boldsymbol{\sigma}'$ | Deviatoric stress tensor | $\text{Pa}$ | `src/numerics/stokes_darcy.jl` (`assemble_hydromechanical_lse!`) |
| **6. Thermal & Segregation** | $\rho c_{p,\mathrm{app}}$ | Apparent volumetric heat capacity including melting latent heat | $\text{J/}(\text{m}^3\cdot\text{K})$ | `src/physics/melting_convection.jl` (`rhocp_apparent_silicate`) |
| | $\mathrm{D}/\mathrm{D}t$ | Material derivative along marker path | $\text{s}^{-1}$ | `src/numerics/thermal.jl` (`assemble_thermal_lse!`) |
| | $k_\mathrm{eff}$ | Effective thermal conductivity with convective scaling | $\text{W/}(\text{m}\cdot\text{K})$ | `src/numerics/thermal.jl` (`assemble_thermal_lse!`) |
| | $H_\mathrm{s}$ | Volumetric shear heating dissipation | $\text{W/m}^3$ | `src/numerics/thermal.jl` (`assemble_thermal_lse!`) |
| | $H_\mathrm{a}$ | Volumetric adiabatic work heating | $\text{W/m}^3$ | `src/numerics/thermal.jl` (`assemble_thermal_lse!`) |
| | $H_\mathrm{seg}$ | Volumetric segregation dissipation heating | $\text{W/m}^3$ | `src/numerics/thermal.jl`, `src/numerics/segregation.jl` |
| | $H_\mathrm{L}$ | Reaction/pyrolysis enthalpy ($DHP$) plus crystallization latent heat ($Q_\mathrm{lat}$) | $\text{W/m}^3$ | `src/numerics/thermal.jl` (`assemble_thermal_lse!`) |
| | $H_\mathrm{metric}$ | Optional spherical metric conduction term | $\text{W/m}^3$ | `src/numerics/thermal.jl` (`assemble_thermal_lse!`) |
| | $\mathbf{J}_\mathrm{Fe}$ | Molten iron settling drift-flux vector | $\text{kg/}(\text{m}^2\cdot\text{s})$ | `src/numerics/segregation.jl` (`apply_metal_segregation!`) |
| | $k_\mathrm{Fe}$ | Metal permeability in silicate matrix | $\text{m}^2$ | `src/numerics/segregation.jl` (`apply_metal_segregation!`) |
| | $\eta_\mathrm{Fe}$ | Liquid metal dynamic viscosity | $\text{Pa}\cdot\text{s}$ | `src/numerics/segregation.jl` (`apply_metal_segregation!`) |
| **7 & 9. Markers & Domain** | $\mathbf{x}_m$ | Marker position vector | $\text{m}$ | `src/simulation/steps/advect_markers.jl` (`advect_markers!`) |
| | $\mathbf{v}$ | Interpolated velocity vector | $\text{m/s}$ | `src/simulation/steps/advect_markers.jl` (`advect_markers!`) |
| | $\Delta \sigma_m$ | Subgrid stress diffusion increment | $\text{Pa}$ | `src/simulation/steps/diffuse_and_update_markers.jl` (`diffuse_and_update_markers!`) |
| | $d_s$ | Subgrid stress diffusion damping parameter | 1 | `src/config.jl` (`dsubgrids`) |
| | $\tau_\mathrm{M}$ | Maxwell viscoelastic relaxation timescale | $\text{s}$ | `src/simulation/steps/diffuse_and_update_markers.jl` |
| **8. Venting & Atmosphere** | $T(\tau_\mathrm{IR})$ | Atmospheric temperature at infrared optical depth $\tau_\mathrm{IR}$ | $\text{K}$ | `src/atmosphere.jl` (`compute_guillot_surface_temperature`) |
| | $\tau_\mathrm{IR}$ | Thermal infrared optical depth coordinate | 1 | `src/atmosphere.jl` (`compute_guillot_surface_temperature`) |
| | $T_\mathrm{int}$ | Interior effective temperature | $\text{K}$ | `src/atmosphere.jl` (`compute_guillot_surface_temperature`) |
| | $T_\mathrm{eq}$ | Planetary equilibrium temperature | $\text{K}$ | `src/atmosphere.jl` (`compute_guillot_surface_temperature`) |
| | $f_\gamma(\tau_\mathrm{IR})$ | Semi-grey atmospheric temperature profile | 1 | `src/atmosphere.jl` (`compute_guillot_surface_temperature`) |
| | $\gamma$ | Visible-to-thermal opacity ratio | 1 | `src/atmosphere.jl` (`compute_guillot_surface_temperature`) |
| | $S_\mathrm{vent}$ | Surface venting volumetric drainage rate | $\text{s}^{-1}$ | `src/simulation/steps/vent_and_degas.jl` (`vent_and_degas!`) |
| | $\phi_\mathrm{min}$ | Minimum residual porosity floor | 1 | `src/config.jl` |
| **Controls & Bounds** | $n_\mathrm{iter}$ | Maximum plastic iterations per step | 1 | `src/config.jl` (`solver.max_plastic_iterations`) |
| | $t_\mathrm{end}$ | Simulation termination physical time | $\text{s}$ | `src/config.jl` (`time.endtime`) |
| | $n_\mathrm{steps}$ | Maximum allowed simulation time steps | 1 | `src/config.jl` (`time.n_steps`) |

### Numerical Execution Pipeline

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
   The code solves reversible serpentinization and dehydration reactions through time-relaxed equilibrium with reaction timescale $\Delta t_\mathrm{r}$.
   The free energy change $\Delta G$ is written for the forward dehydration direction.
   It updates mineral hydration states, fluid production, and reaction latent heat.

5. **Hydro-Mechanical Stokes-Darcy (`assemble_hydromechanical_lse!`, outer iteration)**:
   The solver assembles and solves the linear system for solid velocities, total pressure, fluid Darcy flux, and fluid pressure.
   It updates matrix compaction and plastic yielding.

6. **Thermal Energy & Metal Segregation (`solve_thermal_energy!`, `apply_metal_segregation!`, `apply_silicate_melt_segregation!`, outer iteration)**:
   The code solves the thermal energy equation with conduction and latent heat.
   Advective heat transport occurs through marker advection in Stage 9.
   Hydrothermal convection and soft-turbulence transport enter through the effective thermal conductivity $k_\mathrm{eff}$.
   The code computes metal droplet settling and silicate melt ascent using conservative drift-flux transport.
   It adds the heat $H_\mathrm{seg}$ that metal and melt segregation release (`src/numerics/thermal.jl`, `src/numerics/segregation.jl`).

7. **Marker Viscoelasticity & Subgrid Diffusion (`diffuse_and_update_markers!`)**:
   The code updates viscoelastic stress relaxation on markers after the outer iteration loop converges.
   It conducts subgrid thermal diffusion between grid nodes and markers.

8. **Fluid Venting & Coupled Atmosphere (`vent_and_degas!`, `evolve_atmosphere!`)**:
   The code computes hydrofracture fluid venting and volatile exsolution from magma.
   It couples vented volatiles into surface atmospheric pressure and equilibrium temperature.
   It computes hydrodynamic and kinetic escape losses to space.

9. **Marker Advection, Replenish & Domain Doubling (`advect_markers!`, `replenish!`, `telescope_domain!`)**:
   The code advects markers with a 4th-order Runge-Kutta scheme.
   It updates nodal total and fluid pressures along marker paths.
   It replenishes markers in empty cells.
   The code doubles the domain size when the body radius exceeds 70% of the domain half-width (telescoping, see [Telescoping Domain](telescoping_domain.md)).

10. **Diagnostics & Output (`advance_step_diagnostics!`)**:
    After telescoping, the loop advances physical time $t \leftarrow t + \Delta t$.
    The code logs conservation metrics, saves JLD2 checkpoints, and streams telemetry records.

**Where this is checked.** The execution pipeline and solver convergence loops are tested in [test_simulation_modular_steps.jl](https://github.com/FormingWorlds/Erebus.jl/blob/main/test/test_simulation_modular_steps.jl), [test_integration.jl](https://github.com/FormingWorlds/Erebus.jl/blob/main/test/test_integration.jl), [test_simulation.jl](https://github.com/FormingWorlds/Erebus.jl/blob/main/test/test_simulation.jl), and [test_conservation_closure.jl](https://github.com/FormingWorlds/Erebus.jl/blob/main/test/test_conservation_closure.jl).

---

## Planetesimal Evolution Across Physical Stages

### How does a planetesimal evolve from a cold porous body into a differentiated world?
Planetesimals undergo 3 successive physical stages as radiogenic decay heats their interiors.
Figure 2 summarizes these stages and links them to the physical processes in Erebus.

![Erebus.jl Planetesimal Evolution in Three Stages](../assets/erebus_planetesimal_stages.svg)

*Figure 2: Planetesimal evolution across 3 physical stages. Sectors show successive stages of 1 body; sector widths do not scale with time. Left and right cards show physical processes that Erebus solves, marked with their stage tag. Active modules depend on configuration switches; all modules are shown active. Arrow sizes carry no scale.*

In Figure 2:
- **Stage I (Cold Accretion)**: Low temperatures preserve porous ice and pristine rock. Accretion flux $\dot{M}_\mathrm{acc}$ delivers pebbles and volatile ices to the surface.
- **Stage II (Hydrothermal Alteration)**: Radiogenic decay warms the interior. Ice melts into pore water. Darcy flow $\mathbf{q}_\mathrm{D}$ carries pore water through the permeable crust. Melted ice reacts with rock (serpentinization). Higher temperatures later drive dehydration.
- **Stage III (Magma Ocean & Core Formation)**: Silicate rock melts above the solidus. Above the critical threshold ($F_\mathrm{m} > F_\mathrm{crit} \approx 0.40$, configured via `melting.phi_crit`), the mush becomes a magma ocean. Liquid Fe-FeS metal percolates downward before rock melting and rains downward through magma oceans as drift-flux $\mathbf{J}_\mathrm{Fe}$ to form a metallic core. This descent releases gravitational segregation dissipation heat $H_\mathrm{seg}$ (`src/numerics/thermal.jl`, `src/numerics/segregation.jl`). The planetesimal surface emits radiative heat flux $F_\mathrm{rad}$ into space.

Table 2 defines the physical symbols in Figure 2.

*Table 2: Symbols in Figure 2 evolution diagram.*

| Symbol | Meaning | Unit | Code location |
| :--- | :--- | :--- | :--- |
| $\dot{M}_\mathrm{acc}$ | Pebble and embryo accretion mass flux | $\text{kg/s}$ | `src/simulation/steps/accrete.jl` (`accrete!`) |
| $T_\mathrm{amb}$ | Protoplanetary disk ambient surface temperature | $\text{K}$ | `src/simulation/steps/prepare_step.jl` (`prepare_step_ambient!`) |
| $\mathbf{q}_\mathrm{D}$ | Darcy hydrothermal pore water filtration vector | $\text{m/s}$ | `src/numerics/stokes_darcy.jl` (`assemble_hydromechanical_lse!`) |
| $F_\mathrm{crit}$ | Rheological melt fraction threshold for disaggregation | 1 | `src/config.jl` (`melting.phi_crit`) |
| $\mathbf{J}_\mathrm{Fe}$ | Liquid Fe-FeS metal settling drift-flux vector | $\text{kg/}(\text{m}^2\cdot\text{s})$ | `src/numerics/segregation.jl` (`apply_metal_segregation!`) |
| $F_\mathrm{rad}$ | Surface radiative energy emission flux into space | $\text{W/m}^2$ | `src/simulation/steps/prepare_step.jl` (`apply_surface_radiation!`) |
| $D_i$ | Metal-silicate geochemical partition coefficient for trace element $i$ | 1 | `src/physics/metal_partitioning.jl` |
| $\tau_\mathrm{IR}$ | Infrared atmospheric optical depth coordinate | 1 | `src/atmosphere.jl` (`compute_guillot_surface_temperature`) |
| $R_\mathrm{core}$ | Metallic core boundary radius | $\text{m}$ | `src/numerics/segregation.jl` |
| $R_\mathrm{mo}$ | Magma ocean outer boundary radius | $\text{m}$ | `src/physics/melting_convection.jl` |
| $R_\mathrm{cr}$ | Planetesimal solid crust/lid outer surface radius | $\text{m}$ | `src/geometry.jl` |
| $R_\mathrm{halo}$ | Computational sticky-air buffer outer boundary radius | $\text{m}$ | `src/geometry.jl` |

**Where this is checked.** The [Benchmark Reference Suite](../reference/benchmarks/index.md) documents analytical benchmarks for these physical regimes. Tests in [test_hydrothermal_convection.jl](https://github.com/FormingWorlds/Erebus.jl/blob/main/test/test_hydrothermal_convection.jl), [test_core_formation.jl](https://github.com/FormingWorlds/Erebus.jl/blob/main/test/test_core_formation.jl), and [test_magma_transport.jl](https://github.com/FormingWorlds/Erebus.jl/blob/main/test/test_magma_transport.jl) verify these solutions.

---

## Thermo-Hydro-Mechanical Coupling

`Erebus.jl` resolves these interacting regimes through a fully coupled numerical framework:
- The thermal solver calculates conduction, radioactive decay, latent heat buffering, fluid advective heat transport, and segregation dissipation heating.
- The hydromechanical solver simultaneously solves coupled Stokes solid deformation, Darcy fluid flux, and poroelastic volume changes in a monolithic linear system.
- Marker-in-Cell advects material phases, temperature, composition, and porosity without numerical diffusion across moving boundaries.
- Melting and soft turbulence routines evaluate pressure-dependent silicate melting, apparent heat capacity latent heat buffering, suspension rheology, and regularized sub-grid convective conductivity.
- The magma transport drift-flux solver tracks conservative outward buoyant migration of silicate melt, couples Darcy percolation and hindered Stokes crystal settling, models subsolidus crystallization and dissipation heating, and tracks mantle depletion.
- The metal segregation drift-flux solver tracks conservative downward migration of molten iron, couples percolation and Stokes droplet settling through the rheological transition, and computes gravitational dissipation heating.

---

## Geometry Convention

`Erebus.jl` models the planetesimal as a 2-D Cartesian cross-section passing through the planetary center $(x_{\text{center}}, y_{\text{center}})$.

### Two-Dimensional Form and Three-Dimensional Mapping

The computational domain represents a planar slice through a spherical body.
To evaluate three-dimensional integral quantities (such as total mass, component inventories, and integrated volatile degassing) from planar marker positions, each Lagrangian marker $m$ at distance $r_m = \sqrt{(x_m - x_{\text{center}})^2 + (y_m - y_{\text{center}})^2}$ carries an effective out-of-plane spherical integration length $w_\mathrm{3D}(r_m)$:

$$w_\mathrm{3D}(r_m) = 2 r_m$$

Under the assumption of spherical symmetry, this geometric mapping reinterprets each concentric ring at radius $r_m$ as a spherical shell of radius $r_m$.
Integrating over a circular disk of radius $R$ in the 2-D Cartesian cross-section with differential marker area $A_m = \Delta x_m \Delta y_m$ recovers the exact volume of a sphere:

$$V_{\text{3D}} = \int_{\text{disk}} w_\mathrm{3D}(r) \, \mathrm{d}A = \int_0^{2\pi} \mathrm{d}\theta \int_0^R (2r) \, r \, \mathrm{d}r = 4\pi \int_0^R r^2 \, \mathrm{d}r = \frac{4}{3} \pi R^3$$

For spherically symmetric or radially averaged distributions, mass increments and volatile transfers map between 2-D planar sums and 3-D spherical inventories via:

$$\Delta M = \sum_m \rho_m A_m w_\mathrm{3D}(r_m) \Delta C_m, \quad w_\mathrm{3D}(r_m) = 2 r_m, \quad A_m = \Delta x_m \cdot \Delta y_m$$

where $\rho_m$ is the marker density, $A_m$ is the marker differential area, and $\Delta C_m$ is the dimensionless mass fraction or phase change increment.

### Weighting Limits

The out-of-plane weighting function $w_\mathrm{3D}(r) = 2r$ exhibits three characteristic limits in planetary structures:
1. **Planetary surface limit ($r \to R$)**: Near the outer planetary radius, $w_\mathrm{3D}(R) = 2R$. Surface flux and atmospheric exchange calculations scale with out-of-plane diameter $2R$.
2. **Finite shell mean ($r \in [r_{\text{inner}}, R]$)**: For an outer spherical shell, the ratio of integrated 3-D shell mass to 2-D planar shell mass equals $\frac{4}{3} \frac{R^3 - r_{\text{inner}}^3}{R^2 - r_{\text{inner}}^2}$. For a near-surface shell spanning $0.9R \le r \le R$, this ratio equals approximately $1.902 R$.
3. **Full spherical volume mean**: Integrating over the entire circular section yields the full spherical volume $\frac{4}{3} \pi R^3$, corresponding to an area-weighted mean out-of-plane thickness of $\frac{4}{3} R$.

### Dynamic Flow and Conservation

Hydromechanical Stokes flow, Darcy filtration, and marker transport operate in the 2-D Cartesian cross-section.
Because 2-D Cartesian divergence-free velocity fields ($\nabla \cdot \mathbf{v} = 0$) do not preserve the axisymmetric volume element $r \, \mathrm{d}r$, planar advection does not automatically conserve the 3-D weighted inventory $\sum_m \rho_m A_m (2 r_m)$.
The helper functions `marker_out_of_plane_length` and `marker_area` provide the geometric primitives to evaluate 3-D spherical projections alongside 2-D cross-sectional quantities.

Thermal conduction incorporates an optional radial geometric metric term to account for spherical divergence in conductive heat flow in the cross-section.

### Gravitational Source Term

The 2-D Poisson equation uses the source term $\frac{8}{3}\pi G\rho$.
For a uniform sphere, this term yields exact gravity at every radius.
For differentiated bodies with core formation, select `geometry.gravity_mode = :enclosed_mass` (see [Telescoping Domain](telescoping_domain.md)).

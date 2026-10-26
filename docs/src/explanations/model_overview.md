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
At elevated temperatures ($T \approx 500\text{--}900\text{ K}$), hydrous phyllosilicates, such as serpentine, underwent thermal dehydration reactions, releasing bound hydroxyl groups as free supercritical fluid.

### 4. Matrix Compaction and Pore Overpressure
Viscous, elastic, and plastic deformation of the silicate matrix reduced porosity ($\phi$).
When the compaction rate exceeded the rate of fluid escape via Darcy percolation, pore fluid pressure ($P_f$) rose to match or exceed the lithostatic solid pressure ($P_t$).

### 5. Failure and Hydrofracture
When the Terzaghi effective stress became tensile ($\sigma_{\text{eff}} \le -\sigma_{\text{tensile}}$), the silicate matrix experienced hydrofracturing, creating high-permeability pathways that vented fluids to the planetesimal surface.

### 6. Silicate Rock Melting, Magma Ascent, and Magma Ocean Formation
In planetesimals that accreted early, decay heating drove temperatures past the silicate solidus ($T_{\text{sol}} \approx 1400\text{ K}$).
Rock melted into a crystal-liquid mush and transitioned into a magma ocean above the rheological threshold ($F_{\text{crit}} \approx 0.40$).
Buoyant silicate melt segregates outward via porous Darcy percolation and hindered Stokes crystal settling, releasing gravitational shear dissipation and latent heat of crystallization.
Sub-grid soft turbulence convection transported heat outward to the surface, buffering interior temperatures and governing core segregation.

### 7. Iron Core Formation and Metal Segregation
At temperatures exceeding the Fe-FeS eutectic ($T_{\text{eut}} \approx 1213\text{ K}$), dense metallic liquid melts and segregates inward toward the center.
In partially molten rock, metal trickles downward by porous Darcy percolation.
In magma oceans, metal droplets rain downward via Stokes settling, releasing gravitational potential energy that accelerates core formation.

---

## Global Architecture & Execution Flow

How does the simulation loop advance physical time, and which governing equations does it solve?
The simulation loop couples Eulerian staggered grids with Lagrangian markers.
In each time step, the solver executes 10 sequential stages to model thermal, mechanical, chemical, and atmospheric processes.

### Execution Pipeline and Governing Equations

Figure 1 shows the 10-stage execution pipeline and the governing equations of the simulation loop.

![Erebus.jl Execution Pipeline and Governing Equations](../assets/erebus_execution_pipeline.svg)

*Figure 1: Class D (Schematic Diagram): Erebus execution pipeline and governing equations. The center column shows 1 time step in simulation call order. Side cards show governing equations and physical modules for each step. Module activation depends on configuration switches; all modules are shown active. Outer iterations solve coupled hydration, Stokes-Darcy fluid-matrix mechanics, and thermal energy equations. The non-convergence path restores state and halves the time step.*

Table 1 defines the symbols in the Figure 1 equation cards.
Each 2-D cell is weighted by $w_\mathrm{3D} = 2r$, so the slice volume equals the sphere volume $\frac{4}{3}\pi R^3$ (`src/geometry.jl`, `marker_out_of_plane_length`).
Volumetric radiogenic heating $H_\mathrm{r}$ weights isotope decay powers by phase mass fractions and phase densities (`src/simulation/steps/radiogenic_heating.jl`, `radiogenic_heating!`).
$$H_\mathrm{r} = (1 - \phi) [(1 - X_\mathrm{Fe}) \rho_\mathrm{sil} Q_{^{26}\mathrm{Al}} + X_\mathrm{Fe} \rho_\mathrm{Fe} Q_{^{60}\mathrm{Fe}}]$$
The 2-D Poisson equation uses factor $\frac{8}{3}\pi G\rho$ to match the interior gravity of a 3-D sphere (see [Telescoping Domain](telescoping_domain.md)).

*Table 1: Symbols in Figure 1 equation cards.*

| Symbol | Meaning | Unit | Code location |
| :--- | :--- | :--- | :--- |
| $w_\mathrm{3D}$ | Out-of-plane cylindrical integration width | $\text{m}$ | `src/geometry.jl` (`marker_out_of_plane_length`) |
| $r$ | Radial distance from planetesimal center | $\text{m}$ | `src/geometry.jl` (`marker_out_of_plane_length`) |
| $A$ | Cross-sectional area in 2-D planar slice | $\text{m}^2$ | `src/simulation/steps/accrete.jl` (`accrete!`) |
| $R$ | Planetesimal outer radius | $\text{m}$ | `src/simulation/steps/accrete.jl` (`accrete!`) |
| $Q_i$ | Specific radiogenic heat production rate of isotope $i$ | $\text{W/kg}$ | `src/physics/radioactivity.jl` (`Q_radiogenic`) |
| $f_i$ | Reference isotope mass abundance | $\text{atoms/kg}$ | `src/physics/radioactivity.jl` (`Q_radiogenic`) |
| $r_i$ | Initial isotopic abundance ratio | - | `src/physics/radioactivity.jl` (`Q_radiogenic`) |
| $E_i$ | Nuclear decay energy | $\text{J}$ | `src/physics/radioactivity.jl` (`Q_radiogenic`) |
| $\tau_i$ | Mean lifetime of radionuclide $i$ | $\text{s}$ | `src/physics/radioactivity.jl` (`Q_radiogenic`) |
| $t$ | Elapsed physical simulation time | $\text{s}$ | `src/physics/radioactivity.jl` (`Q_radiogenic`) |
| $H_\mathrm{r}$ | Volumetric radiogenic heating rate | $\text{W/m}^3$ | `src/simulation/steps/radiogenic_heating.jl` (`radiogenic_heating!`) |
| $\phi$ | Matrix porosity | - | `src/simulation/steps/radiogenic_heating.jl` (`radiogenic_heating!`) |
| $X_\mathrm{Fe}$ | Iron metal mass fraction | - | `src/simulation/steps/radiogenic_heating.jl` (`radiogenic_heating!`) |
| $\rho_\mathrm{sil}$ | Silicate grain density | $\text{kg/m}^3$ | `src/physics/radioactivity.jl` (`calculate_radioactive_heating`) |
| $\rho_\mathrm{Fe}$ | Iron metal density | $\text{kg/m}^3$ | `src/physics/radioactivity.jl` (`calculate_radioactive_heating`) |
| $\mathrm{W_{sil}}$ | Hydrous silicate phase (serpentine) | - | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| $\mathrm{D_{sil}}$ | Dry silicate phase (olivine) | - | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| $\mathrm{H_2O}$ | Pore fluid water | - | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| $\Delta G$ | Molar Gibbs free energy change of reaction | $\text{J/mol}$ | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| $\Delta H$ | Standard molar reaction enthalpy | $\text{J/mol}$ | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| $T$ | Temperature | $\text{K}$ | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| $\Delta S$ | Standard molar reaction entropy | $\text{J/(mol}\cdot\text{K)}$ | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| $P_\mathrm{f}$ | Pore fluid pressure | $\text{Pa}$ | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| $\Delta V$ | Standard molar volume change | $\text{m}^3/\text{mol}$ | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| $R_\mathrm{gas}$ | Universal gas constant | $\text{J/(mol}\cdot\text{K)}$ | `src/constants.jl` (`R_GAS`) |
| $X_\mathrm{D}, X_\mathrm{W}$ | Molar fractions of dry and hydrous silicate | - | `src/physics/serpentinization.jl` (`perform_thermochemical_reaction!`) |
| $T(\tau)$ | Atmospheric temperature at optical depth $\tau$ | $\text{K}$ | `src/atmosphere.jl` (`compute_guillot_surface_temperature`) |
| $\tau$ | Thermal optical depth coordinate | - | `src/atmosphere.jl` (`compute_guillot_surface_temperature`) |
| $T_\mathrm{int}$ | Interior effective temperature | $\text{K}$ | `src/atmosphere.jl` (`compute_guillot_surface_temperature`) |
| $T_\mathrm{eq}$ | Planetary equilibrium temperature | $\text{K}$ | `src/atmosphere.jl` (`compute_guillot_surface_temperature`) |
| $f_\gamma(\tau)$ | Semi-grey atmospheric temperature profile | - | `src/atmosphere.jl` (`compute_guillot_surface_temperature`) |
| $\gamma$ | Visible-to-thermal opacity ratio | - | `src/atmosphere.jl` (`compute_guillot_surface_temperature`) |
| $\Phi$ | Gravitational potential proxy | $\text{m}^2/\text{s}^2$ | `src/simulation/steps/solve_gravity.jl` (`solve_gravity!`) |
| $G$ | Gravitational constant | $\text{m}^3/(\text{kg}\cdot\text{s}^2)$ | `src/constants.jl` (`G_GRAV`) |
| $\rho$ | Total bulk mass density | $\text{kg/m}^3$ | `src/simulation/steps/solve_gravity.jl` (`solve_gravity!`) |
| $F_\mathrm{m}$ | Equilibrium silicate melt mass fraction | - | `src/physics/melting_convection.jl` (`compute_melt_fraction`) |
| $T_\mathrm{s}, T_\mathrm{l}$ | Solidus and liquidus temperatures | $\text{K}$ | `src/physics/melting_convection.jl` (`compute_melt_fraction`) |
| $\boldsymbol{\sigma}'$ | Deviatoric stress tensor | $\text{Pa}$ | `src/numerics/stokes_darcy.jl` (`assemble_hydromechanical_lse!`) |
| $P$ | Total mixture pressure | $\text{Pa}$ | `src/numerics/stokes_darcy.jl` (`assemble_hydromechanical_lse!`) |
| $\mathbf{g}$ | Gravitational acceleration vector | $\text{m/s}^2$ | `src/numerics/stokes_darcy.jl` (`assemble_hydromechanical_lse!`) |
| $\mathbf{v}_\mathrm{s}$ | Solid matrix velocity vector | $\text{m/s}$ | `src/numerics/stokes_darcy.jl` (`assemble_hydromechanical_lse!`) |
| $\rho c_p$ | Volumetric heat capacity | $\text{J/(m}^3\cdot\text{K)}$ | `src/numerics/thermal.jl` (`assemble_thermal_lse!`) |
| $k_\mathrm{eff}$ | Effective thermal conductivity | $\text{W/(m}\cdot\text{K)}$ | `src/numerics/thermal.jl` (`assemble_thermal_lse!`) |
| $H_\mathrm{s}$ | Volumetric shear heating dissipation | $\text{W/m}^3$ | `src/numerics/thermal.jl` (`assemble_thermal_lse!`) |
| $H_\mathrm{a}$ | Volumetric adiabatic work heating | $\text{W/m}^3$ | `src/numerics/thermal.jl` (`assemble_thermal_lse!`) |
| $Q_\mathrm{seg}$ | Segregation dissipation heating rate | $\text{W/m}^3$ | `src/numerics/thermal.jl` (`assemble_thermal_lse!`) |
| $H_\mathrm{L}$ | Latent heat of reactions and melting | $\text{W/m}^3$ | `src/numerics/thermal.jl` (`assemble_thermal_lse!`) |
| $\mathbf{x}_m$ | Marker position vector | $\text{m}$ | `src/simulation/steps/advect_markers.jl` (`advect_markers!`) |
| $\mathbf{v}$ | Interpolated velocity vector | $\text{m/s}$ | `src/simulation/steps/advect_markers.jl` (`advect_markers!`) |


### Numerical Execution Pipeline

The simulation loop advances physical time in 10 sequential stages:

1. **Accretion & Ambient (`prepare_step_ambient!`, `accrete!`)**:
   The code computes disk ambient temperature and pressure.
   It deposits accreted pebbles or embryos onto the body.
   It advances planetary radius $R(t)$ and updates sticky air boundaries.

2. **Heating & Pyrolysis (`radiogenic_heating!`, `update_pyrolysis!`)**:
   The code computes analytical decay rates of short-lived radionuclides ($^{26}\text{Al}$, $^{60}\text{Fe}$).
   It adds volumetric heat sources to markers.
   It tracks volatile release from organic carbon pyrolysis.

3. **Marker-in-Cell Mapping & Gravity (`interpolate_markers_to_grid!`, `solve_gravity!`, `apply_surface_radiation!`)**:
   The code interpolates marker properties to staggered grid nodes.
   It solves the Poisson equation for gravitational potential $\Phi$ and gravity $\mathbf{g}$.
   It applies radiative surface boundary conditions.

4. **Hydration Reaction (`perform_thermochemical_reaction!`)**:
   The code solves reversible serpentinization and dehydration reactions based on Gibbs free energy minimization.
   It updates mineral hydration states, fluid production, and reaction latent heat.

5. **Hydro-Mechanical Stokes-Darcy Solve (`assemble_hydromechanical_lse!`)**:
   The solver assembles and solves the linear system for solid velocities, total pressure, fluid Darcy flux, and fluid pressure.
   It updates matrix compaction and plastic yielding.

6. **Thermal Energy & Metal Segregation (`solve_thermal_energy!`, `apply_metal_segregation!`, `apply_silicate_melt_segregation!`)**:
   The solver solves the thermal energy equation with conduction, advection, and latent heat.
   It computes metal droplet settling and silicate melt ascent using conservative drift-flux transport.
   It adds gravitational segregation dissipation heating $Q_\mathrm{seg}$ (`src/numerics/thermal.jl`, `src/numerics/segregation.jl`).

7. **Marker Viscoelasticity & Subgrid Diffusion (`diffuse_and_update_markers!`)**:
   The code updates viscoelastic stress relaxation on markers after the outer iteration loop converges.
   It conducts subgrid thermal diffusion between grid nodes and markers.

8. **Fluid Venting & Coupled Atmosphere (`vent_and_degas!`, `evolve_atmosphere!`)**:
   The code computes hydrofracture fluid venting and volatile exsolution from magma.
   It couples vented volatiles into surface atmospheric pressure and equilibrium temperature.
   It computes hydrodynamic and kinetic escape losses to space.

9. **Marker Advection & Replenishment (`advect_markers!`, `replenish!`, `telescope_domain!`)**:
   The code advects markers with a 4th-order Runge-Kutta scheme.
   It updates nodal total and fluid pressures along marker paths.
   It replenishes markers in empty cells and doubles the domain when the body radius exceeds 70% of domain half-width.

10. **Diagnostics & Telemetry (`advance_step_diagnostics!`)**:
    The code advances physical time $t \leftarrow t + \Delta t$.
    It logs conservation metrics, saves JLD2 checkpoints, and streams telemetry records.

**Where this is checked.** The execution pipeline and solver convergence loops are tested in `test/test_integration.jl`, `test/test_simulation.jl`, and `test/test_conservation_closure.jl`.

---

## Planetesimal Evolution Across Physical Stages

How does a planetesimal evolve from a cold dust aggregate into a differentiated body?
Planetesimals undergo 3 successive physical stages as radiogenic decay heats their interiors.
Figure 2 summarizes these stages and links them to the physical processes in Erebus.

![Erebus.jl Planetesimal Evolution in Three Stages](../assets/erebus_planetesimal_stages.svg)

*Figure 2: Class D (Schematic Diagram): Planetesimal evolution across 3 physical stages. Sectors show successive stages of 1 body; sector widths do not scale with time. Left and right cards show physical processes that Erebus solves. Active modules depend on configuration switches; all modules are shown active. Flux arrows are schematic.*

In Figure 2:
- **Stage I (Cold Accretion)**: Low temperatures preserve porous ice and pristine rock. Accretion flux $\dot{M}_\mathrm{acc}$ delivers pebbles and volatile ices to the surface.
- **Stage II (Hydrothermal Alteration)**: Radiogenic decay warms the interior. Ice melts into pore water. Darcy percolation $\mathbf{q}_\mathrm{D}$ drives hydrothermal circulation, serpentinization, and dehydration reactions through permeable crust.
- **Stage III (Magma Ocean & Core Formation)**: Heating drives silicate melting above the critical rheological threshold ($F_\mathrm{m} > F_\mathrm{crit} \approx 0.40$, configured via `phi_crit`). Molten iron drift-flux $\mathbf{J}_\mathrm{Fe}$ descends to form a metallic core. This descent releases gravitational segregation dissipation heating $Q_\mathrm{seg}$ (`src/numerics/thermal.jl`, `src/numerics/segregation.jl`). Surface radiative heat flux $F_\mathrm{rad}$ radiates into space and drives atmospheric escape.

**Where this is checked.** The [Benchmark Reference Suite](../reference/benchmarks/index.md) documents analytical benchmarks for these physical regimes. Tests in `test/test_hydrothermal_convection.jl`, `test/test_core_formation.jl`, and `test/test_magma_transport.jl` verify these solutions.

## Thermo-Hydro-Mechanical Coupling

`Erebus.jl` resolves these interacting regimes through a fully coupled numerical framework:
- Thermal solver calculates conduction, radioactive decay, latent heat buffering, fluid advective heat transport, and segregation dissipation heating.
- Hydromechanical solver simultaneously solves coupled Stokes solid deformation, Darcy fluid flux, and poroelastic volume changes in a monolithic linear system.
- Marker-in-Cell advects material phases, temperature, composition, and porosity without numerical diffusion across moving boundaries.
- Melting and soft turbulence routines evaluate pressure-dependent silicate melting, apparent heat capacity latent heat buffering, suspension rheology, and regularized sub-grid convective conductivity.
- Magma transport drift-flux solver tracks conservative outward buoyant migration of silicate melt, couples Darcy percolation and hindered Stokes crystal settling, models subsolidus crystallization and dissipation heating, and tracks mantle depletion.
- Metal segregation drift-flux solver tracks conservative downward migration of molten iron, couples percolation and Stokes droplet settling through the rheological transition, and computes gravitational dissipation heating.

---

## Geometry Convention

Erebus.jl models the planetesimal as a 2D Cartesian cross-section passing through the planetary center $(x_{\text{center}}, y_{\text{center}})$.

### Two-Dimensional Form and Three-Dimensional Mapping

The computational domain represents a planar slice through a spherical body. To evaluate three-dimensional integral quantities (such as total mass, component inventories, and integrated volatile degassing) from planar marker positions, each Lagrangian marker $m$ at distance $r_m = \sqrt{(x_m - x_{\text{center}})^2 + (y_m - y_{\text{center}})^2}$ carries an effective out-of-plane spherical integration length $L(r_m)$:

\[
L(r_m) = 2 r_m
\]

Under the assumption of spherical symmetry, this geometric mapping reinterprets each concentric ring at radius $r_m$ as a spherical shell of radius $r_m$. Integrating over a circular disk of radius $R$ in the 2D Cartesian cross-section with differential marker area $A_m = \Delta x_m \Delta y_m$ recovers the exact volume of a sphere:

\[
V_{\text{3D}} = \int_{\text{disk}} L(r) \, dA = \int_0^{2\pi} d\theta \int_0^R (2r) \, r \, dr = 4\pi \int_0^R r^2 \, dr = \frac{4}{3} \pi R^3
\]

For spherically symmetric or radially averaged distributions, mass increments and volatile transfers map between 2D planar sums and 3D spherical inventories via:

\[
\Delta M = \sum_m \rho_m A_m L(r_m) \Delta C_m, \quad L(r_m) = 2 r_m, \quad A_m = \Delta x_m \cdot \Delta y_m
\]

where $\rho_m$ is the marker density, $A_m$ is the marker differential area, and $\Delta C_m$ is the dimensionless mass fraction or phase change increment.

### Weighting Limits

The out-of-plane weighting function $L(r) = 2r$ exhibits three characteristic limits in planetary structures:
1. **Planetary surface limit ($r \to R$)**: Near the outer planetary radius, $L(R) = 2R$. Surface flux and atmospheric exchange calculations scale with out-of-plane diameter $2R$.
2. **Finite shell mean ($r \in [r_{\text{inner}}, R]$)**: For an outer spherical shell, the ratio of integrated 3D shell mass to 2D planar shell mass equals $\frac{4}{3} \frac{R^3 - r_{\text{inner}}^3}{R^2 - r_{\text{inner}}^2}$. For a near-surface shell spanning $0.9R \le r \le R$, this ratio evaluates to approximately $1.902 R$.
3. **Full spherical volume mean**: Integrating over the entire circular section yields the full spherical volume $\frac{4}{3} \pi R^3$, corresponding to an area-weighted mean out-of-plane thickness of $\frac{4}{3} R$.

### Dynamic Flow and Conservation

Hydromechanical Stokes flow, Darcy filtration, and marker transport operate in the 2D Cartesian cross-section. Because 2D Cartesian divergence-free velocity fields ($\nabla \cdot \mathbf{v} = 0$) do not preserve the axisymmetric volume element $r \, dr$, planar advection does not automatically conserve the 3D-weighted inventory $\sum_m \rho_m A_m (2 r_m)$. The helper functions `marker_out_of_plane_length` and `marker_area` provide the geometric primitives to evaluate 3D spherical projections alongside 2D cross-sectional quantities.

Thermal conduction incorporates an optional radial geometric metric term to account for spherical divergence in conductive heat flow in the cross-section.

### Gravitational Source Term

The Poisson solver evaluates gravitational potential from the 2D density field using an effective source term $(8/3) \pi G \rho$. This source coefficient reproduces the exact radial gravitational acceleration at the surface of a uniform-density sphere of radius $R$.



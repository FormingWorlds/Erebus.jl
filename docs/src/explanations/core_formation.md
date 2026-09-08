# Iron Core Formation and Metal Segregation

This section documents the physical theory, regime transitions, governing continuum equations, and numerical implementation of iron core formation in `Erebus.jl`.

---

## Physical Background and Motivation

Planetesimal differentiation is the primary planetary process separating metallic iron-nickel-sulfur liquids from silicate rock to form a central metallic core. In the early solar system, internal heating from short-lived radionuclides ($^{26}\text{Al}$ and $^{60}\text{Fe}$) heats planetesimal interiors.

Because Fe-FeS mixtures melt at a eutectic temperature ($T_{\text{eutectic}} \approx 1213\text{ K}$) well below the silicate solidus ($T_{\text{sol,sil}} \approx 1400\text{ K}$), metal liquefies while the silicate matrix remains largely or entirely solid. As heating continues, progressive melting of silicates transforms the mechanical environment from a rigid or deforming porous solid rock into a disaggregated, convective magma ocean.

The transport of dense molten metal through this changing silicate environment spans two fundamental regimes:
1. **Porous Darcy percolation**: Dense molten metal trickles downward through the pore channels of a solid or partially molten silicate crystalline matrix.
2. **Stokes droplet settling**: Molten metal breaks into discrete liquid droplets that rain downward through a liquid silicate magma ocean.

`Erebus.jl` couples these two transport mechanisms into a unified drift-flux formulation across the entire melting history of the planetesimal.

---

## Metal Melting Thermodynamics

The metallic phase consists of an iron-nickel alloy with light elements, primarily sulfur. Liquid metal formation is modeled as a partial melting process parameterized by the eutectic temperature $T_{\text{eutectic}}$ and a melting interval $\Delta T_{\text{metal}}$:

$$F_{\text{fe}}(T) = \begin{cases}
0, & T \le T_{\text{eutectic}} \\
\frac{T - T_{\text{eutectic}}}{\Delta T_{\text{metal}}}, & T_{\text{eutectic}} < T < T_{\text{eutectic}} + \Delta T_{\text{metal}} \\
1, & T \ge T_{\text{eutectic}} + \Delta T_{\text{metal}}
\end{cases}$$

For a marker with bulk metal volume fraction $X_{\text{fe,bulk}}$, the local molten metal volume fraction $\phi_m$ is:

$$\phi_m = X_{\text{fe,bulk}} \cdot F_{\text{fe}}(T)$$

When $\phi_m > 0$, the marker density, thermal conductivity, and volumetric heat capacity are updated via volume-weighted mixture rules:

$$\rho_{\text{eff}} = (1 - \phi_m) \rho_{\text{sil}} + \phi_m \rho_{\text{metal}}$$

$$k_{\text{eff}} = (1 - \phi_m) k_{\text{sil}} + \phi_m k_{\text{metal}}$$

$$(\rho C_p)_{\text{eff}} = (1 - \phi_m) (\rho C_p)_{\text{sil}} + \phi_m (\rho C_p)_{\text{metal}}$$

---

## Segregation Regimes

```
Solid Silicate Matrix                      Partially Molten Silicate                  Magma Ocean
(F_m < F_settle_start)                     (F_settle_start <= F_m <= F_perc_end)     (F_m > F_perc_end)
+------------------------+                 +------------------------+                 +------------------------+
| Darcy Percolation      |                 | Regime Transition      |                 | Stokes Droplet         |
| v_perc ~ k_fe / eta_fe |  ============>  | Cubic Hermite          |  ============>  | Settling               |
| Permeability power law |                 | Blending of v_seg      |                 | v_settle ~ d^2 Delta_rho|
+------------------------+                 +------------------------+                 +------------------------+
```

The segregation velocity magnitude $v_{\text{seg}}$ depends on the degree of silicate melting $F_m \in [0, 1]$.

### Regime 1: Porous Darcy Percolation ($F_m < F_{\text{settle\_start}}$)

In solid or low-melt silicate rock, molten iron moves by gravity-driven percolation through an interconnected network of pores and grain boundaries.

Because the interfacial energy between molten iron and solid silicates is high (dihedral wetting angle $\theta > 60^\circ$ at low pressures), metal melt forms isolated pockets until a critical percolation threshold $\phi_{\text{crit,perc}} \approx 0.05$ is exceeded. Above this threshold, interconnected channels permit flow.

The interstitial pore segregation velocity $v_{\text{perc}}$ for mobile metal transport is:

$$v_{\text{perc}} = \frac{k_{\text{metal}}(\phi_m)}{\phi_m \, \eta_{\text{metal}}} \Delta\rho \, g \left(\frac{\phi_m - \phi_{\text{residual}}}{\phi_m}\right)$$

where:
- $\Delta\rho = \rho_{\text{metal}} - \rho_{\text{silicate}}$ is the positive density contrast between metal and silicate.
- $\eta_{\text{metal}}$ is the dynamic viscosity of molten iron-sulfur alloy ($\sim 10^{-2}\text{ Pa}\cdot\text{s}$).
- $g$ is the local gravitational acceleration directed radially inward.
- $\phi_{\text{residual}}$ is the residual trapped threshold below which metal is immobilized in matrix pores.

The effective permeability for molten metal $k_{\text{metal}}(\phi_m)$ follows a modified Kozeny-Carman formulation:

$$k_{\text{metal}}(\phi_m) = \begin{cases}
0, & \phi_m \le \phi_{\text{crit,perc}} \\
k_{\text{ref}} \left(\frac{\phi_m - \phi_{\text{crit,perc}}}{\phi_0}\right)^n \left(\frac{1 - (\phi_m - \phi_{\text{crit,perc}})}{1 - \phi_0}\right)^{-2}, & \phi_{\text{crit,perc}} < \phi_m < \phi_{\text{pack}}
\end{cases}$$

where $\phi_{\text{crit,perc}}$ is the percolation threshold ($0.05$), $k_{\text{ref}}$ is the reference permeability ($10^{-9}\text{ m}^2$), $\phi_0$ is the reference porosity ($0.10$), and $n$ is the permeability exponent ($n = 3$).

### Regime 2: Stokes Droplet Settling ($F_m > F_{\text{perc\_end}}$)

When silicate melting exceeds the rheological transition ($F_m \ge 0.40\text{ to }0.50$), the solid silicate framework breaks down. The environment becomes a low-viscosity suspension or liquid magma ocean. In this regime, molten metal can no longer percolate through pore channels; instead, it emulsifies into liquid droplets that settle under gravity.

The terminal settling velocity of an isolated spherical metal droplet of radius $r_{\text{drop}}$ in a fluid of dynamic viscosity $\eta_{\text{susp}}$ is given by the Stokes-Hadamard-Rybczynski formulation:

$$v_{\text{stokes}} = \frac{2}{9} \frac{\Delta\rho \, g \, r_{\text{drop}}^2}{\eta_{\text{susp}}} \cdot C_{\text{HR}}$$

The Hadamard-Rybczynski mobility factor $C_{\text{HR}}$ accounts for internal circulation within the liquid metal droplet:

$$C_{\text{HR}} = \frac{3 \eta_{\text{metal}} + 3 \eta_{\text{susp}}}{2 \eta_{\text{susp}} + 3 \eta_{\text{metal}}}$$

- For a rigid sphere or when surface active contaminants immobilize the interface ($\eta_{\text{metal}} \gg \eta_{\text{susp}}$), $C_{\text{HR}} \to 1$.
- For a clean liquid droplet settling through a viscous liquid ($\eta_{\text{susp}} \gg \eta_{\text{metal}}$), $C_{\text{HR}} \to 1.5$.

#### Droplet Size Determination

Droplet size is controlled by the balance between disruptive hydrodynamic shear forces and restorative surface tension forces, characterized by the droplet Weber number $\text{We} = \rho_{\text{sil}} v^2 d / \sigma$. Three modes are available:

1. `:fixed`: Constant prescribed droplet diameter $d_{\text{fixed}} = 2 r_{\text{drop}}$ (default $1.0\text{ cm}$).
2. `:weber_mean`: Gravity-capillary equilibrium balance where maximum droplet size is limited by Rayleigh-Taylor/Kelvin-Helmholtz instabilities:
   $$d_{\text{weber}} = \sqrt{\frac{\text{We}_{\text{crit}} \, \sigma}{\Delta\rho \, g}}$$
3. `:weber_turbulent`: Dynamic breakup based on the estimated settling velocity:
   $$d_{\text{weber}} = \frac{\text{We}_{\text{crit}} \, \sigma}{\rho_{\text{sil}} \, v_{\text{est}}^2}$$

#### Hindered Settling and Maximum Packing

At high metal volume fractions, droplet-droplet hydrodynamic interactions reduce the settling velocity. This hindrance is modeled via the Richardson-Zaki power law:

$$h(\phi_m) = \left(1 - \frac{\phi_m}{\phi_{\text{pack}}}\right)^m$$

where $m$ is the hindrance exponent (default $4.5$), and $\phi_{\text{pack}}$ is the maximum packing fraction (default $0.65$). Furthermore, as metal ponds toward maximum packing ($\phi_m \to \phi_{\text{pack}}$), settling terminates as droplets touch and coalesce.

### Transition Regime ($F_{\text{settle\_start}} \le F_m \le F_{\text{perc\_end}}$)

Between the onset of droplet settling ($F_{\text{settle\_start}} = 0.40$) and complete matrix breakdown ($F_{\text{perc\_end}} = 0.50$), both processes occur concurrently as melt pockets coalesce into liquid pools. `Erebus.jl` blends the two velocity limits using a $C^1$-continuous cubic Hermite polynomial weight $\xi(F_m)$:

$$\theta = \frac{F_m - F_{\text{settle\_start}}}{F_{\text{perc\_end}} - F_{\text{settle\_start}}}$$

$$\xi(F_m) = 3 \theta^2 - 2 \theta^3$$

$$v_{\text{seg}}(F_m) = (1 - \xi) \, v_{\text{perc}} + \xi \, v_{\text{settle}}$$

---

## Continuum Transport: Drift-Flux Formulation

Metal segregation is implemented as an explicit drift-flux transport step operating on top of the bulk hydrodynamic advection of the rock-fluid mixture.

The conservation of total metal mass is:

$$\frac{\partial \phi_m}{\partial t} + \nabla \cdot (\phi_m \mathbf{v}_{\text{matrix}}) + \nabla \cdot (\phi_m \mathbf{v}_{\text{seg}}) = 0$$

In the operator-split architecture of `Erebus.jl`:
1. Bulk advection by matrix flow $\nabla \cdot (\phi_m \mathbf{v}_{\text{matrix}})$ is handled by the marker Runge-Kutta advection scheme (`advect_markers!`).
2. Relative drift $\nabla \cdot (\phi_m \mathbf{v}_{\text{seg}})$ is solved on the Eulerian staggered grid by `apply_metal_segregation!`.

### Finite-Volume Discretization and Subcycling

Because local Stokes settling velocities in low-viscosity magma oceans can reach $10^{-2}\text{ m/s}$ (traversing a grid cell in minutes to hours), while the outer geodynamic timestep $dt$ spans years, the drift-flux update is subcycled with an internal CFL-constrained sub-timestep:

$$\Delta t_{\text{cfl}} = \text{CFL} \cdot \frac{\min(\Delta x, \Delta y)}{\max |v_{\text{seg}}|}$$

$$N_{\text{sub}} = \min\left(\left\lceil \frac{\Delta t}{\Delta t_{\text{cfl}}} \right\rceil, N_{\text{max}}\right), \quad \Delta t_{\text{sub}} = \frac{\Delta t}{N_{\text{sub}}}$$

During each subcycle:
1. Unscaled requested fluxes across cell faces are computed using donor-cell upwinding based on the gravity direction.
2. Multi-dimensional flux limiters prevent donor cells from evacuating below residual fraction $\phi_{\text{res}}$ and receiver cells from overfilling beyond maximum packing $\phi_{\text{pack}}$.
3. The cell metal mass is updated conservatively.
4. Net cell mass increments are distributed back to Lagrangian markers within each cell.
5. Floating-point roundoff conservation correction maintains total metal mass conservation to machine precision ($< 10^{-12}$).

---

## Gravitational Dissipation Heating

Core formation releases gravitational potential energy as heavy iron drops toward the planetesimal center. The volumetric energy dissipation rate $Q_{\text{seg}}$ [$\text{W/m}^3$] is:

$$Q_{\text{seg}} = \Delta\rho \, g \, v_{\text{seg}} \, \phi_m$$

This volumetric dissipation acts as a thermal source in the planetary energy balance:

$$\rho C_p \frac{DT}{dt} = \nabla \cdot (K \nabla T) + H_{\text{radio}} + Q_{\text{seg}} + Q_{\text{shear}}$$

In large planetesimals and protoplanets ($R > 100\text{ km}$), segregation heating can raise internal temperatures by tens to hundreds of Kelvins, accelerating core formation and establishing a positive feedback loop toward complete planetary melting.

---

## Source Code Architecture

| Physical Component | Source File | Key Functions |
|:---|:---|:---|
| Parameter definition | `src/config.jl` | `CoreFormationConfig` |
| Melt fraction & velocity | `src/physics.jl` | `compute_metal_melt_fraction`, `metal_segregation_velocity`, `segregation_dissipation_heating` |
| Marker tracking & properties | `src/particles.jl` | `compute_marker_properties!`, `setup_marker_metal_properties`, `replenish_markers!` |
| Conservative transport | `src/numerics.jl` | `apply_metal_segregation!`, `assemble_thermal_lse!` |
| Main simulation integration | `src/simulation.jl` | Timestep loop sequence |

---

## References

- Deguen, R., Olson, P., & Cardin, P. (2011). Driving magma ocean convection by core formation. *Earth and Planetary Science Letters*, 310(3-4), 303-313.
- Deguen, R., Landeau, M., & Olson, P. (2014). Tunnelling of liquid metal through a magma ocean. *Earth and Planetary Science Letters*, 391, 274-287.
- Lichtenberg, T., Golabek, G. J., Burn, R., Meyer, M. R., Alibert, Y., Gerya, T. V., & Mordasini, C. (2019). A water budget divide in developing exoplanetary systems. *Nature Astronomy*, 3(4), 307-313.
- Lichtenberg, T., Bower, D. J., Hammond, M., Boukrouche, R., Sanan, P., Tsai, S. M., & Pierrehumbert, R. T. (2021). Vertically resolved magma ocean-protoatmosphere evolution. *Journal of Geophysical Research: Planets*, 126(2), e2020JE006711.
- Monteux, J., Ricard, Y., Coltice, N., Dubuffet, F., & Aguilar, M. (2009a). A model of metal-silicate separation on growing planets. *Geophysical Journal International*, 179(1), 515-526.
- Monteux, J., Jellinek, A. M., & Buffett, B. A. (2009b). Heating of the early Earth by core formation: Physical mechanisms and thermal impact. *Journal of Geophysical Research*, 114(B6), B06404.
- Rubie, D. C., Melosh, H. J., Reid, J. E., Liebske, C., & Righter, K. (2003). Mechanisms of metal-silicate equilibration in the terrestrial magma ocean. *Earth and Planetary Science Letters*, 205(3-4), 239-255.
- Stevenson, D. J. (1990). Fluid dynamics of core formation. In *Origin of the Earth* (pp. 231-249). Oxford University Press.
- Yoshino, T., Walter, M. J., & Katsura, T. (2003). Core formation in planetesimals triggered by permeable flow. *Nature*, 422(6928), 154-157.

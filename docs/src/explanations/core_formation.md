# Iron Core Formation and Metal Segregation

This section documents the physical theory, regime transitions, governing continuum equations, and numerical implementation of iron core formation in `Erebus.jl`.

---

## Physical Background and Motivation

Planetesimal differentiation is the primary planetary process separating metallic iron-nickel-sulfur liquids from silicate rock to form a central metallic core. In the early solar system, internal heating from short-lived radionuclides ($^{26}\text{Al}$ and $^{60}\text{Fe}$) heats planetesimal interiors.

Because Fe-FeS mixtures melt at a eutectic temperature ($T_{\text{eutectic}} \approx 1213\text{ K}$) well below the silicate solidus ($T_{\text{sol,sil}} \approx 1400\text{ K}$), metal liquefies while the silicate matrix remains largely or entirely solid. As heating continues, progressive melting of silicates transforms the mechanical environment from a rigid or deforming porous solid rock into a disaggregated, convective magma ocean.

The transport of dense molten metal through this changing silicate environment spans two fundamental regimes:
1. **Porous Darcy percolation**: Dense molten metal trickles downward through the pore channels of a solid or partially molten silicate crystalline matrix.
2. **Stokes droplet settling**: Molten metal breaks into discrete liquid droplets that rain downward through a liquid silicate magma ocean.

`Erebus.jl` couples these two transport mechanisms into a unified drift-flux formulation throughout the melting history of the planetesimal.

---

## Metal Melting Thermodynamics

The metallic phase consists of an iron-nickel alloy with light elements, primarily sulfur. Liquid metal formation is modeled as a partial melting process parameterized by the eutectic temperature $T_{\text{eutectic}}$ and a melting interval $\Delta T_{\text{metal}}$:

$$F_{\text{fe}}(T) = \begin{cases}
0, & T \le T_{\text{eutectic}} \\
\frac{T - T_{\text{eutectic}}}{\Delta T_{\text{metal}}}, & T_{\text{eutectic}} < T < T_{\text{eutectic}} + \Delta T_{\text{metal}} \\
1, & T \ge T_{\text{eutectic}} + \Delta T_{\text{metal}}
\end{cases}$$

For a marker with bulk metal volume fraction $X_{\text{fe,bulk}}$, the local molten metal volume fraction $\phi_m$ available for segregation transport is:

$$\phi_m = X_{\text{fe,bulk}} \cdot F_{\text{fe}}(T)$$

The local metallic phase density transitions from solid ($\rho_{\text{metal,solid}}$) to liquid ($\rho_{\text{metal,liquid}}$) across the melting interval:

$$\rho_{\text{metal}}(T) = (1 - F_{\text{fe}}) \rho_{\text{metal,solid}} + F_{\text{fe}} \rho_{\text{metal,liquid}}$$

### Liquid Metal Density Equations of State (`metal_density_mode`)

The molten metal density $\rho_{\text{metal,liquid}}$ is governed by the light-element alloy chemistry, configured via `metal_density_mode`:

1. `:sanloup2000` (default): Fe-FeS liquid alloy formulation from Sanloup et al. (2000) parameterized by the sulfur mass fraction $w_S \in [0, 0.40]$ (`sulfur_fraction`, default $0.31$ for eutectic composition):
   $$\rho_0(w_S) = 7020.0 - 5050.0 \cdot w_S \quad [\text{kg/m}^3]$$
   Coupled with thermal expansion ($\alpha_m = 10^{-4}\text{ K}^{-1}$) and isothermal compressibility ($K_T = 65\text{ GPa}$):
   $$\rho_{\text{liquid}}(T, P) = \rho_0(w_S) \left[1 - \alpha_m (T - T_0) + \frac{P}{K_T}\right]$$
2. `:morard2014`: Fe-S liquid alloy formulation from Morard et al. (2014):
   $$\rho_0(w_S) = 7020.0 \cdot (1.0 - 0.72 \cdot w_S) \quad [\text{kg/m}^3]$$
3. `:constant`: Fixed prescribed liquid metal density $\rho_{\text{metal}}$ (default $5450\text{ kg/m}^3$).

The marker composite density and thermal conductivity blend the host rock or ice matrix with metal across all temperatures using bulk metal volume fraction $X_{\text{fe,bulk}}$:

$$\rho_{\text{eff}} = (1 - X_{\text{fe,bulk}}) \rho_{\text{sil}} + X_{\text{fe,bulk}} \rho_{\text{metal}}(T)$$

$$k_{\text{eff}} = (1 - X_{\text{fe,bulk}}) k_{\text{sil}} + X_{\text{fe,bulk}} k_{\text{metal}}$$

During metal melting ($T_{\text{eutectic}} \le T \le T_{\text{eutectic}} + \Delta T_{\text{metal}}$), the apparent volumetric heat capacity of the metal incorporates the latent heat of melting $L_{\text{metal}}$ using the solid metal reference density:

$$(\rho C_p)_{\text{eff,metal}} = (\rho C_p)_{\text{metal}} + \rho_{\text{metal,solid}} \frac{L_{\text{metal}}}{\Delta T_{\text{metal}}}$$

$$(\rho C_p)_{\text{eff}} = (1 - X_{\text{fe,bulk}}) (\rho C_p)_{\text{sil}} + X_{\text{fe,bulk}} (\rho C_p)_{\text{eff,metal}}$$

Solid and liquid metal share baseline conductivity $k_{\text{metal}}$ and unbuffered heat capacity $(\rho C_p)_{\text{metal}}$.

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

| Symbol | Description | Units |
|:---|:---|:---|
| $\Delta\rho$ | Density contrast: $\rho_{\text{metal}} - \rho_{\text{silicate}}$ | $\text{kg/m}^3$ |
| $\eta_{\text{metal}}$ | Dynamic viscosity of molten iron-sulfur alloy ($\sim 10^{-2}\text{ Pa}\cdot\text{s}$) | $\text{Pa}\cdot\text{s}$ |
| $g$ | Local gravitational acceleration directed radially inward | $\text{m/s}^2$ |
| $\phi_{\text{residual}}$ | Residual trapped threshold below which metal is immobilized | - |

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

Droplet size is controlled by the balance between disruptive hydrodynamic shear forces and restorative surface tension forces. Five modes are supported:

1. `:capillary_mean` (default): Gravity-capillary equilibrium balance where maximum droplet size is limited by Rayleigh-Taylor instabilities:
   $$d_{\text{cap}} = \sqrt{\frac{\text{We}_{\text{crit}} \, \sigma}{\Delta\rho \, g}}$$
2. `:bond_mean`: Equivalent alias to `:capillary_mean`.
3. `:weber_mean`: Equivalent alias to `:capillary_mean`.
4. `:weber_turbulent`: Dynamic breakup based on estimated relative settling velocity:
   $$d_{\text{turb}} = \frac{\text{We}_{\text{crit}} \, \sigma}{\rho_{\text{sil}} \, v_{\text{est}}^2}$$
5. `:fixed`: Constant prescribed droplet diameter $d_{\text{fixed}} = 2 r_{\text{drop}}$ (default $0.5\text{ cm} = 5.0\times 10^{-3}\text{ m}$).

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

## Metal-Silicate Volatile Partitioning

During planetesimal differentiation and core formation, volatile elements (H, C, N, and S) partition between molten metallic iron alloys and silicate melts. The thermodynamic distribution is quantified by the metal-silicate partition coefficient $D_i^{\text{met/sil}}$:

$$D_i^{\text{met/sil}} = \frac{C_i^{\text{metal}}}{C_i^{\text{silicate}}}$$

where $C_i^{\text{metal}}$ is the concentration of element $i$ in the molten metal phase [ppmw] and $C_i^{\text{silicate}}$ is the concentration dissolved in the coexisting silicate melt [ppmw].

### Empirical Parameterizations

`Erebus.jl` parameterizes $D_i^{\text{met/sil}}(T, P, \Delta\text{IW}, w_S)$ as functions of temperature $T$ [K], pressure $P$ [Pa], oxygen fugacity $\Delta\text{IW}$ relative to the iron-wüstite buffer, and alloy sulfur mass fraction $w_S$:

1. **Carbon Partitioning (`model_carbon`)**:
   - `:grewal2019` (default): Grewal et al. (2019b) parameterization accounting for the strong suppression of carbon siderophile affinity by dissolved sulfur:
     $$\log_{10} D_C = 1.80 + \frac{2200}{T} - 1.5 \times 10^{-8} \frac{P}{T} - 0.25 \, \Delta\text{IW} + 4.2 \ln(1 - X_S)$$
     where $X_S$ is the mole fraction of sulfur in the Fe-S liquid alloy. In sulfur-free metallic iron, carbon is strongly siderophile ($D_C \approx 2000 - 3000$). At the Fe-FeS eutectic ($w_S \approx 0.31$, $X_S \approx 0.44$), $D_C$ drops sharply to $\approx 15 - 30$.
   - `:fischer2020`: High-pressure parameterization from Fischer et al. (2020):
     $$\log_{10} D_C = 1.50 + \frac{2500}{T} - 1.2 \times 10^{-8} \frac{P}{T} - 0.20 \, \Delta\text{IW}$$
   - `:constant`: Fixed prescribed value `D_C_const`.

2. **Nitrogen Partitioning (`model_nitrogen`)**:
   - `:grewal2019` (default): Grewal et al. (2019a, 2019b) parameterization:
     $$\log_{10} D_N = 0.85 + \frac{1200}{T} - 0.25 \, \Delta\text{IW} + 0.60 \ln(1 - X_S)$$
     Nitrogen is moderately siderophile ($D_N \approx 20 - 50$). Because the sulfur interaction term ($0.60 \ln(1 - X_S)$) is much smaller than that for carbon ($4.2 \ln(1 - X_S)$), nitrogen partitioning is comparatively insensitive to sulfur content. Consequently, core formation in sulfur-rich planetesimals lowers the metallic C/N ratio, generating superchondritic C/N in the residual silicate mantle.
   - `:constant`: Fixed prescribed value `D_N_const`.

3. **Hydrogen Partitioning (`model_hydrogen`)**:
   - `:clesi2018` (default): Clesi et al. (2018) low-pressure parameterization:
     $$\log_{10} D_H = -0.80 + \frac{300}{T} + 5.0 \times 10^{-8} \frac{P}{T} + 0.05 \, \Delta\text{IW}$$
     In the low-pressure planetesimal regime ($P < 1\text{ GPa}$), hydrogen is moderately siderophile to lithophile ($D_H \approx 0.1 - 1.0$). Stoichiometric conversion connects silicate water concentration to elemental hydrogen via $f_H = (2 \times 1.00794 / 18.01528) \times 10^4 \approx 1118.98\text{ ppmw H}$ per $1\text{ wt}\% \text{ H}_2\text{O}$.
   - `:constant`: Fixed prescribed value `D_H_const`.

4. **Sulfur Partitioning (`model_sulfur`)**:
   - `:boujibar2014` (default): Boujibar et al. (2014) parameterization:
     $$\log_{10} D_S = 2.80 - \frac{800}{T} + 1.0 \times 10^{-10} P - 0.20 \, \Delta\text{IW}$$
     Sulfur partitions strongly into metallic liquid ($D_S \approx 100 - 500$), concentrating primordial sulfur into the segregated metallic core.
   - `:constant`: Fixed prescribed value `D_S_const`.

### Elemental Mass Conservation

When molten metal ($F_{\text{fe}} > 0$) coexists with silicate melt ($F_{\text{melt}} > 0$), elemental volatile mass is conserved across both phases:

$$M_{i,\text{total}} = m_{\text{sil}} C_{i,\text{sil}} + m_{\text{met}} C_{i,\text{met}}$$

where $m_{\text{sil}} = \phi_{\text{sil}} \rho_{\text{sil}}$ and $m_{\text{met}} = \phi_{\text{fe}} F_{\text{fe}} \rho_{\text{met}}$. The thermodynamic equilibrium concentration in the silicate melt is:

$$C_{i,\text{sil}}^{\text{eq}} = \frac{M_{i,\text{total}}}{m_{\text{sil}} + D_i m_{\text{met}}}$$

Kinetic exchange advances toward equilibrium with rate fraction $\alpha_{\text{eq}} \in [0, 1]$ (`equilibration_rate`):

$$\Delta C_{i,\text{sil}} = \alpha_{\text{eq}} \left(C_{i,\text{sil}}^{\text{eq}} - C_{i,\text{sil}}\right)$$

$$\Delta C_{i,\text{met}} = -\Delta C_{i,\text{sil}} \left(\frac{m_{\text{sil}}}{m_{\text{met}}}\right)$$

This formulation conserves total volatile mass to floating-point precision on every marker.

---

## Advective Core Segregation Transport of Volatiles

As molten metallic droplets segregate downward under gravity, metal-hosted volatiles are advected alongside the metallic mass flux. Within `apply_metal_segregation!`, the volatile mass flux across cell face $(i, j)$ is computed using upwind donor-cell concentrations:

$$F_{i, k}^x = F_{\text{fe}, x} \cdot \left(\frac{M_{\text{fe}, k}}{M_{\text{fe}}}\right)_{\text{donor}}$$

$$F_{i, k}^y = F_{\text{fe}, y} \cdot \left(\frac{M_{\text{fe}, k}}{M_{\text{fe}}}\right)_{\text{donor}}$$

where $k \in \{H, C, N, S\}$ and the donor cell is selected by the sign of the metallic mass flux. Cell volatile inventories update conservatively:

$$M_{\text{fe}, k}^{n+1} = M_{\text{fe}, k}^n + \Delta t \left(F_{w, k} - F_{e, k} + F_{n, k} - F_{s, k}\right)$$

### Dynamic Sulfur Density Feedback

When `dynamic_sulfur_density = true`, the local sulfur content of the metallic alloy $w_S = X_{\text{fe}, S} \cdot 10^{-6}$ is evaluated dynamically on each marker and grid cell. The local alloy density feeds back into droplet buoyancy $\Delta\rho = \rho_{\text{metal}}(w_S) - \rho_{\text{silicate}}$:

$$\rho_0(w_S) = 7020.0 - 5050.0 \cdot w_S \quad [\text{kg/m}^3]$$

Sulfur-poor metal differentiates faster due to its higher density contrast, whereas sulfur-rich metal exhibits lower settling velocities and longer segregation timescales.

### Integrated Core Volatile Budgets

Integrated core mass and volatile budgets are evaluated using `compute_core_volatile_budgets`:

$$M_{\text{core}, k} = \sum_{m \in \text{core}} \phi_{\text{fe}, m} \, \rho_{\text{metal}} \, \left(X_{\text{fe}, k, m} \cdot 10^{-6}\right)$$

The resulting mean core volatile concentrations $w_{\text{core}, k} = M_{\text{core}, k} / M_{\text{core}, \text{metal}}$ are compared with empirical concentrations measured in magmatic iron meteorites (groups IIAB, IIIAB, IVA, IVB).

---

## Normative Accessory Mineral Tracking and Meteorite Diagnostics

While the binary Fe-FeS eutectic is ~1261 K at low pressure, minor nickel, phosphorus, and carbon depress initial melting to ~1213 K (940 °C; Goldstein et al., 2009; Chabot and Drake, 1999). Below this temperature ($T \le T_{\text{eutectic}}$), minor and volatile elements in solid metallic iron exsolve into accessory phases:
1. Troilite ($\text{FeS}$): Formed from sulfur via molar conversion $w_{\text{troilite}} = w_S \cdot (M_{\text{FeS}} / M_S) \approx 2.742 \cdot w_S$.
2. Schreibersite ($(\text{Fe,Ni})_3\text{P}$): Formed from phosphorus with nickel molar fraction $x_{\text{Ni}} = 0.25$, yielding $w_{\text{schreibersite}} \approx 6.478 \cdot w_P$.
3. Cohenite ($(\text{Fe,Ni})_3\text{C}$) and graphite ($\text{C}$): Formed from carbon up to the carbide saturation limit ($w_{C,\text{max}} \approx 0.0667$). Excess carbon precipitates as elemental graphite.
4. Nitrides ($\text{Fe}_4\text{N}$, $\text{CrN}$, or $\text{TiN}$): Formed from nitrogen according to the selected mode.
5. Metal matrix: Residual solid metallic alloy.

### Thermal Dissolution

During heating through the eutectic interval ($T_{\text{eutectic}} \le T \le T_{\text{eutectic}} + \Delta T_{\text{transition}}$), solid accessory phases dissolve into metallic liquid:

$$X_{m, i} = F_{\text{solid}}(T) \cdot w_i^{\text{stoich}}$$

where $F_{\text{solid}}(T) = 1.0 - (T - T_{\text{eutectic}}) / \Delta T_{\text{transition}}$.

### Regional Classification

Marker phase distributions are integrated across radial domains (core, mantle, crust). Simulation outputs are classified into meteorite affinities:
- Magmatic differentiated bodies ($f_{\text{molten,core}} \ge 0.80$, $f_{\text{metal,core}} \ge 0.40$, and $f_{\text{solid,crust}} \le 0.005$): Corresponds to groups IIIAB, IVA, and IVB.
- Primitive incomplete bodies ($f_{\text{molten,core}} \le 0.60$ and $f_{\text{solid,crust}} \ge 0.01$): Corresponds to IAB complexes and winonaites.
- Transitional bodies: Intermediate segregation states.

---

## Source Code Architecture

| Physical Component | Source File | Key Functions |
|:---|:---|:---|
| Parameter definition | `src/config.jl` | `CoreFormationConfig`, `MetalPartitionConfig`, `PhaseTrackingConfig` |
| Partition coefficients | `src/physics.jl` | `compute_metal_silicate_partition_coefficient`, `compute_metal_silicate_partition_coefficients` |
| Volatile equilibration | `src/physics.jl` | `equilibrate_metal_silicate_volatiles!` |
| Accessory minerals | `src/physics.jl` | `compute_troilite_stoichiometry`, `compute_schreibersite_stoichiometry`, `compute_cohenite_graphite_stoichiometry`, `compute_nitride_stoichiometry`, `compute_normative_mineral_assemblage` |
| Regional modes | `src/physics.jl` | `compute_regional_mineral_modes` |
| Core budget integration | `src/physics.jl` | `compute_core_volatile_budgets` |
| Melt fraction and velocity | `src/physics.jl` | `compute_metal_melt_fraction`, `metal_segregation_velocity`, `segregation_dissipation_heating` |
| Marker tracking and properties | `src/particles.jl` | `compute_marker_properties!`, `setup_marker_metal_properties`, `setup_marker_metal_volatile_properties`, `setup_marker_phase_tracking_properties`, `replenish_markers!` |
| Conservative transport | `src/numerics.jl` | `apply_metal_segregation!`, `assemble_thermal_lse!` |
| Main simulation integration | `src/simulation.jl` | Timestep loop sequence and checkpoint persistence |

---

## References

- Benedix, G. K., McCoy, T. J., Keil, K., & Bogard, D. D. (2000). A petrologic and geochemical study of winonaites: Implications for trace element behavior during primitive achondrite differentiation. *Geochimica et Cosmochimica Acta*, 64(14), 2535-2553.
- Boujibar, A., Andrault, D., Bolfan-Casanova, N., Bouhifd, M. A., & Kawamoto, T. (2014). Metal-silicate partitioning of sulphur, new experimental constraints by EMPA and SIMS. *Earth and Planetary Science Letters*, 391, 42-54.
- Chabot, N. L., & Drake, M. J. (1999). Crystallization of magmatic iron meteorites: The role of phosphorus and sulfur. *Meteoritics & Planetary Science*, 34(2), 235-246.
- Clesi, V., Bouhifd, M. A., Bolfan-Casanova, N., Manthilake, G., Schiavi, F., Kawamoto, T., & Andrault, D. (2018). Low hydrogen contents in Earth's core. *Science Advances*, 4(3), e1701876.
- Deguen, R., Olson, P., & Cardin, P. (2011). Experiments on turbulent metal-silicate mixing in a magma ocean. *Earth and Planetary Science Letters*, 310(3-4), 303-313.
- Deguen, R., Landeau, M., & Olson, P. (2014). Turbulent metal-silicate mixing, fragmentation, and equilibration in magma oceans. *Earth and Planetary Science Letters*, 391, 274-287.
- Fischer, R. A., Cottrell, E., Hauri, E., Lee, K. K. M., & Le Voyer, M. (2020). The partitioning of carbon and oxygen between core and mantle in the early Earth. *Proceedings of the National Academy of Sciences*, 117(16), 8743-8749.
- Goldstein, J. I., Scott, E. R. D., & Chabot, N. L. (2009). Iron meteorites: Crystallization, thermal history, parent bodies, and origin. *Chemie der Erde - Geochemistry*, 69(4), 293-325.
- Grewal, D. S., Dasgupta, R., Sun, C., Tsuno, K., & Costin, G. (2019a). Delivery of carbon, nitrogen, and sulfur to the silicate Earth by a planetary merger. *Science Advances*, 5(1), eaau3669.
- Grewal, D. S., Dasgupta, R., & Farnell, A. (2019b). The speciation of carbon, nitrogen, and water in magma oceans and its effect on volatile partitioning between metal and silicate. *Geochimica et Cosmochimica Acta*, 251, 87-115.
- Lichtenberg, T., Golabek, G. J., Burn, R., Meyer, M. R., Alibert, Y., Gerya, T. V., & Mordasini, C. (2019). A water budget dichotomy of rocky protoplanets from 26Al-heating. *Nature Astronomy*, 3(4), 307-313.
- Lichtenberg, T., Bower, D. J., Hammond, M., Boukrouche, R., Sanan, P., Tsai, S. M., & Pierrehumbert, R. T. (2021). Vertically resolved magma ocean-protoatmosphere evolution. *Journal of Geophysical Research: Planets*, 126(2), e2020JE006711.
- Monteux, J., Ricard, Y., Coltice, N., Dubuffet, F., & Aguilar, M. (2009a). A model of metal-silicate separation on growing planets. *Geophysical Journal International*, 179(1), 515-526.
- Monteux, J., Jellinek, A. M., & Buffett, B. A. (2009b). Heating of the early Earth by core formation: Physical mechanisms and thermal impact. *Journal of Geophysical Research*, 114(B6), B06404.
- Rubie, D. C., Melosh, H. J., Reid, J. E., Liebske, C., & Righter, K. (2003). Mechanisms of metal-silicate equilibration in the terrestrial magma ocean. *Earth and Planetary Science Letters*, 205(3-4), 239-255.
- Stevenson, D. J. (1990). Fluid dynamics of core formation. In *Origin of the Earth* (pp. 231-249). Oxford University Press.
- Yoshino, T., Walter, M. J., & Katsura, T. (2003). Core formation in planetesimals triggered by permeable flow. *Nature*, 422(6928), 154-157.


# Iron Core Formation and Metal Segregation

How does molten iron leave silicate rock to collect in a central core, and how does core formation alter the thermal evolution and volatile budget of a planetesimal?
Iron-nickel-sulfur alloy starts to melt near $1213\text{ K}$, below the silicate solidus near $1400\text{ K}$.
Dense metallic liquid first percolates downward through pores in solid rock.
Once silicate melting disrupts the crystalline matrix, metal settles as droplets.
This segregation releases gravitational energy as heat, and sinking metal and rising silicate exchange sensible heat.
It also carries hydrogen (H), carbon (C), nitrogen (N), and sulfur (S) toward the core and sets the core radius reported in output files.

---

## Physical Background and Motivation

Planetesimal differentiation separates metallic Fe-Ni-S liquids from silicate rock to form a central metallic core.
In the early solar system, short-lived radionuclides ($^{26}\text{Al}$ and $^{60}\text{Fe}$) heat planetesimal interiors ([Lichtenberg et al., 2019](../reference/bibliography.md#Primary-References)).

Fe-Ni-S alloy starts to melt near $T_{\text{eutectic}} \approx 1213\text{ K}$ (below the pure Fe-FeS eutectic near $1261\text{ K}$ due to Ni, P, and C), below the silicate solidus ($T_{\text{sol,sil}} \approx 1400\text{ K}$).
Progressive silicate melting then turns the rigid rock matrix into a convective magma ocean.

Dense molten metal moves through this environment in two regimes:
1. **Porous Darcy percolation**: Liquid metal trickles downward along silicate grain boundaries and pores.
2. **Stokes droplet settling**: Molten metal forms liquid droplets that sink through magma.

`Erebus.jl` couples both regimes in one drift-flux formulation during planetesimal melting.

---

## Metal Melting Thermodynamics

The metallic phase consists of an iron-nickel alloy with light elements, primarily sulfur.
Liquid metal formation is modeled as a partial melting process parameterized by the eutectic temperature $T_{\text{eutectic}}$ and a melting interval $\Delta T_{\text{metal}}$:

$$F_{\text{fe}}(T) = \begin{cases}
0, & T \le T_{\text{eutectic}} \\
\frac{T - T_{\text{eutectic}}}{\Delta T_{\text{metal}}}, & T_{\text{eutectic}} < T < T_{\text{eutectic}} + \Delta T_{\text{metal}} \\
1, & T \ge T_{\text{eutectic}} + \Delta T_{\text{metal}}
\end{cases}$$

For a marker with bulk metal volume fraction $X_{\text{fe,bulk}}$, the local molten metal volume fraction $\phi_m$ available for segregation transport is:

$$\phi_m = X_{\text{fe,bulk}} \cdot F_{\text{fe}}(T)$$

**Known limit.** After each segregation call, marker molten fraction resets to $X_{\text{fe,bulk}} \cdot F_{\text{fe}}(T)$ from local temperature without tracking refractory Fe enrichment in solid residue.
Repeated steps at constant sub-liquidus temperature therefore slowly drain residual metal toward $\phi_{\text{trap}} / F_{\text{fe}}(T)$, where $\phi_{\text{trap}}$ is the capillary trap floor (defined in Regime 1).

The local metallic phase density transitions from solid ($\rho_{\text{metal,solid}}$) to liquid ($\rho_{\text{metal,liquid}}$) over the melting interval:

$$\rho_{\text{metal}}(T) = (1 - F_{\text{fe}}) \rho_{\text{metal,solid}} + F_{\text{fe}} \rho_{\text{metal,liquid}}$$

### Liquid Metal Density Equations of State (`metal_density_mode`)

`metal_density_mode` selects the liquid metal density $\rho_{\text{metal,liquid}}$.
Both alloy laws define a reference density $\rho_0(w_S)$ at $T_0 = 1500\text{ K}$ from sulfur mass fraction $w_S \in [0, 0.40]$ (`sulfur_fraction`, default $0.31$):

1. `:sanloup2000` (default; [Sanloup et al., 2000](../reference/bibliography.md#Core-Formation-and-Metal-Silicate-Segregation)):
   $$\rho_0(w_S) = 7020.0 - 5050.0 \cdot w_S \quad [\text{kg/m}^3]$$
2. `:morard2014` ([Morard et al., 2014](../reference/bibliography.md#Core-Formation-and-Metal-Silicate-Segregation)):
   $$\rho_0(w_S) = 7020.0 \cdot (1.0 - 0.72 \cdot w_S) \quad [\text{kg/m}^3]$$
3. `:constant`: Fixed liquid density $\rho_{\text{metal}}$ (default $5450\text{ kg/m}^3$).

Both alloy laws apply the thermal-expansion and compression factor with $\alpha_m = 10^{-4}\text{ K}^{-1}$ and isothermal bulk modulus $K_T = 65\text{ GPa}$:

$$\rho_{\text{liquid}}(T, P) = \rho_0(w_S) \max\left(0.5,\ 1 - \alpha_m (T - T_0) + \frac{P}{K_T}\right)$$

Marker properties pass marker pressure $P$, whereas segregation cells evaluate $P = 0$.

Marker density and thermal conductivity blend rock-plus-pore silicate matrix ($\rho_{\text{sil}}$, $k_{\text{sil}}$) and metal with $X_{\text{fe,bulk}}$:

$$\rho_{\text{eff}} = (1 - X_{\text{fe,bulk}}) \rho_{\text{sil}} + X_{\text{fe,bulk}} \rho_{\text{metal}}(T)$$

$$k_{\text{eff}} = (1 - X_{\text{fe,bulk}}) k_{\text{sil}} + X_{\text{fe,bulk}} k_{\text{metal}}$$

Within $T_{\text{eutectic}} \le T \le T_{\text{eutectic}} + \Delta T_{\text{metal}}$, metal heat capacity includes latent heat $L_{\text{metal}}$ weighted by $\rho_{\text{metal,solid}}$:

$$(\rho C_p)_{\text{eff,metal}} = (\rho C_p)_{\text{metal}} + \rho_{\text{metal,solid}} \frac{L_{\text{metal}}}{\Delta T_{\text{metal}}}$$

$$(\rho C_p)_{\text{eff}} = (1 - X_{\text{fe,bulk}}) (\rho C_p)_{\text{sil}} + X_{\text{fe,bulk}} (\rho C_p)_{\text{eff,metal}}$$

Solid and liquid metal share baseline $k_{\text{metal}}$ and $(\rho C_p)_{\text{metal}}$.

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

Because the interfacial energy between molten iron and solid silicates is high (dihedral wetting angle $\theta > 60^\circ$ at low pressures), metal melt forms isolated pockets until a critical percolation threshold $\phi_{\text{crit,perc}} \approx 0.05$ is exceeded.
Above this threshold, interconnected channels permit flow ([Yoshino et al., 2003](../reference/bibliography.md#Core-Formation-and-Metal-Silicate-Segregation)).

The interstitial pore segregation velocity $v_{\text{perc}}$ for mobile metal transport follows Darcy's law:

$$v_{\text{perc}} = \frac{k_{\text{metal}}(\phi_m)}{\phi_m \, \eta_{\text{metal}}} \Delta\rho \, g$$

for $\phi_m > \max(\phi_{\text{crit,perc}}, \phi_{\text{trap}})$, where $\Delta\rho = \rho_{\text{metal}} - \rho_{\text{sil,eff}}(F_m)$ uses the melt-blended silicate density $\rho_{\text{sil,eff}} = (1 - F_m)\rho_{\text{sil,solid}} + F_m\rho_{\text{sil,melt}}$.

| Symbol | Description | Units |
|:---|:---|:---|
| $\Delta\rho$ | Density contrast: $\rho_{\text{metal}} - \rho_{\text{sil,eff}}$ | $\text{kg/m}^3$ |
| $\eta_{\text{metal}}$ | Dynamic viscosity of molten iron-sulfur alloy ($\sim 10^{-2}\text{ Pa}\cdot\text{s}$) | $\text{Pa}\cdot\text{s}$ |
| $g$ | Local gravitational acceleration directed radially inward | $\text{m/s}^2$ |
| $\phi_{\text{trap}}$ | Effective trapped metal threshold reduced smoothly from $\phi_{\text{residual}}$ to zero over $[F_{\text{settle\_start}}, F_{\text{perc\_end}}]$ | - |

The effective permeability for molten metal $k_{\text{metal}}(\phi_m)$ follows a modified Kozeny-Carman formulation:

$$k_{\text{metal}}(\phi_m) = \begin{cases}
0, & \phi_m \le \phi_{\text{crit,perc}} \\
k_{\text{ref}} \left(\frac{\phi_m - \phi_{\text{crit,perc}}}{\phi_0}\right)^n \left(\frac{1 - (\phi_m - \phi_{\text{crit,perc}})}{1 - \phi_0}\right)^{-2}, & \phi_{\text{crit,perc}} < \phi_m < \phi_{\text{pack}}
\end{cases}$$

where $\phi_{\text{crit,perc}}$ is the percolation threshold ($0.05$), $k_{\text{ref}}$ is the reference permeability ($10^{-9}\text{ m}^2$), $\phi_0$ is the reference porosity ($0.10$), and $n$ is the permeability exponent ($n = 3$).
`effective_metal_trap_fraction` lowers $\phi_{\text{trap}}$ smoothly from $\phi_{\text{residual}}$ to $0$ over $[F_{\text{settle\_start}}, F_{\text{perc\_end}}]$ via a cubic Hermite weight.
This reduction applies only when droplet settling is active.

### Regime 2: Stokes Droplet Settling ($F_m > F_{\text{perc\_end}}$)

Above $F_{\text{perc\_end}}$ (`[coreformation] F_perc_end`, default $0.50$), metal moves by droplet settling alone.
Silicate rheology uses a separate threshold, $\phi_{\text{crit}}$ (`[melting] phi_crit`, default $0.40$); it shares a default value with $F_{\text{settle\_start}}$ but configures matrix breakdown independently.

`compute_silicate_suspension_viscosity` evaluates exponential melt weakening $\eta_{\text{solid}} \exp(-\alpha_\eta F_m)$ below $\phi_{\text{crit}}$.
Over melt-fraction width `F_blend_width` ($0.10$) above $\phi_{\text{crit}}$, it blends in log space onto the crystal-suspension relation ([Roscoe, 1952](../reference/bibliography.md#Silicate-Melting-and-Magma-Rheology); [Solomatov, 2015](../reference/bibliography.md#Silicate-Melting-and-Magma-Rheology)):

$$\eta_{\text{susp}}(F_m) = \eta_{\text{melt}} \left(1 - \frac{1 - F_m}{1 - \phi_{\text{crit}}}\right)^{-2.5}$$

Here $\eta_{\text{solid}}$ is melt-free silicate viscosity [$\text{Pa}\cdot\text{s}$], $\alpha_\eta$ is the weakening exponent (`[melting] alpha_eta`), and $\eta_{\text{melt}}$ is silicate melt viscosity (`[melting] eta_melt`) [$\text{Pa}\cdot\text{s}$].
$\eta_{\text{susp}}$ clamps to $[\eta_{\min}, \eta_{\max}]$, with $\eta_{\min} = \min(0.01\text{ Pa}\cdot\text{s}, \eta_{\text{melt}})$ and `[solver] etamax`.
Settling uses this suspension viscosity rather than the matrix law in [Silicate melting](rock_melting.md) and `[solver] etamin`.

Terminal speed of an isolated droplet of radius $r_{\text{drop}}$ follows Stokes-Hadamard-Rybczynski flow:

$$v_{\text{stokes}} = \frac{2}{9} \frac{\Delta\rho \, g \, r_{\text{drop}}^2}{\eta_{\text{susp}}} \cdot C_{\text{HR}}$$

The factor $C_{\text{HR}} = (3 \eta_{\text{metal}} + 3 \eta_{\text{susp}}) / (2 \eta_{\text{susp}} + 3 \eta_{\text{metal}})$ activates only when `hadamard_rybczynski = true` (default `false` sets $C_{\text{HR}} = 1$).
The regime blend uses $v_{\text{settle}} = v_{\text{stokes}} \, h(\phi_m)$, where $h(\phi_m)$ is the hindered-settling factor defined below.

#### Droplet Size Determination

Droplet size depends on interfacial surface tension $\sigma$ (`sigma_metal_silicate` [N/m]) and critical Weber number $\text{We}_{\text{crit}}$ (`We_crit`) ([Stevenson, 1990](../reference/bibliography.md#Core-Formation-and-Metal-Silicate-Segregation); [Rubie et al., 2003](../reference/bibliography.md#Core-Formation-and-Metal-Silicate-Segregation); [Deguen et al., 2011](../reference/bibliography.md#Core-Formation-and-Metal-Silicate-Segregation); [Deguen et al., 2014](../reference/bibliography.md#Core-Formation-and-Metal-Silicate-Segregation)):

1. `:capillary_mean` (default): Gravity-capillary balance against Rayleigh-Taylor breakup:
   $$d_{\text{cap}} = \sqrt{\frac{\text{We}_{\text{crit}} \, \sigma}{\Delta\rho \, g}}$$
2. `:bond_mean`: Equivalent alias to `:capillary_mean`.
3. `:weber_mean`: Equivalent alias to `:capillary_mean`.
4. `:weber_turbulent`: Dynamic breakup from Stokes terminal speed estimate $v_{\text{est}}$ evaluated at `droplet_diameter_fixed / 2`:
   $$d_{\text{turb}} = \frac{\text{We}_{\text{crit}} \, \sigma}{\rho_{\text{sil,eff}} \, v_{\text{est}}^2}$$
5. `:fixed`: Constant prescribed droplet diameter $d_{\text{fixed}} = 2 r_{\text{drop}}$ (default $0.5\text{ cm} = 5.0\times 10^{-3}\text{ m}$).

#### Hindered Settling and Maximum Packing

At high metal volume fractions, droplet-droplet hydrodynamic interactions reduce the settling velocity.
This hindrance follows the [Richardson and Zaki (1954)](../reference/bibliography.md#Silicate-Melting-and-Magma-Rheology) power law:

$$h(\phi_m) = \left(1 - \frac{\phi_m}{\phi_{\text{pack}}}\right)^m$$

where $m$ is the hindrance exponent (default $4.5$), and $\phi_{\text{pack}}$ is maximum packing fraction (default $0.65$).
Settling flux $q_{\text{settle}}(\phi) = v_{\text{stokes},0} \, \phi \, (1 - \phi/\phi_{\text{pack}})^m$ [m/s] uses unhindered Stokes speed $v_{\text{stokes},0}$ and peaks at $\phi^* = \phi_{\text{pack}}/(m + 1)$.
Each active face evaluates the scalar Godunov flux on $[ \phi_R, \phi_D ]$ using donor mobility $v_{0,D}$ (unhindered donor Stokes speed) to avoid oscillations above $\phi^*$.

### Transition Regime ($F_{\text{settle\_start}} \le F_m \le F_{\text{perc\_end}}$)

Between the onset of droplet settling ($F_{\text{settle\_start}} = 0.40$) and complete matrix breakdown ($F_{\text{perc\_end}} = 0.50$), both processes occur concurrently as melt pockets coalesce into liquid pools.
`Erebus.jl` blends the two velocity limits using a $C^1$-continuous cubic Hermite polynomial weight $\xi(F_m)$ with normalized coordinate $s$:

$$s = \frac{F_m - F_{\text{settle\_start}}}{F_{\text{perc\_end}} - F_{\text{settle\_start}}}$$

$$\xi(F_m) = 3 s^2 - 2 s^3$$

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

Metal can cross several grid cells in one outer step, so the solver subcycles under the Courant-Friedrichs-Lewy factor CFL (`cfl_settling`, default $0.5$) using characteristic wave speed $v_{\text{char}} = \max(v_{\text{seg}}, |\partial q/\partial \phi_m|)$, where $q$ is metal flux [m/s]:

$$\Delta t_{\text{cfl}} = \text{CFL} \cdot \frac{\min(\Delta x, \Delta y)}{\max_{\text{active faces}} v_{\text{char}}}$$

Each substep takes $\Delta t_{\text{sub}} = \min(t_{\text{rem}}, \Delta t_{\text{cfl}})$, where $t_{\text{rem}}$ is remaining outer step time, so no substep exceeds $\Delta t_{\text{cfl}}$.

**Known limit.** If subcycling reaches `max_subcycles` (default `2000`) with $t_{\text{rem}} > 0$, it stops, logs a warning (`maxlog=10`), and freezes transport over $t_{\text{rem}}$.

!!! details "Subcycle operations in `apply_metal_segregation!`"
    1. Weight interior markers by $w_{3\text{D},m} = 2 r_m$ (`spherical_metric = true`) or $1$ (2D planar) to accumulate cell weight $W_{i,j}$, bulk metal $M_{\text{fe},i,j}$, and immutable solid metal $m_{\text{solid},i,j} = M_{\text{fe},i,j} - M_{\text{molten},i,j}$.
    2. Evaluate active-face rates along gravity from upwind percolation plus the Godunov settling flux.
    3. Limit donor outflows above $m_{\text{solid}} + \phi_{\text{trap}} W_{\text{donor}}$ (keeping liquid capillary floor $\phi_{\text{trap}}$ alongside solid metal) and cap receiver inflows at $\phi_{\text{pack}} W_{\text{recv}}$.
    4. Update cell metal and volatile inventories conservatively.
    5. Map negative cell changes onto molten markers proportionally to molten metal fraction $X_{\text{fem},m}$, and positive changes proportionally to $\phi_{\text{pack}} - X_{\text{fe,bulk},m}$.
    6. Correct roundoff to conserve $\sum_m X_{\text{fe,bulk},m} w_{3\text{D},m}$ within $10^{-12}$.

**Where this is checked.** [Core formation benchmark](../reference/benchmarks/core_formation.md).

---

## Energetics of Core Segregation

### Gravitational Dissipation

How much gravitational energy does core formation release?
Sinking metal releases gravitational energy as heat ([Monteux et al., 2009](../reference/bibliography.md#Core-Formation-and-Metal-Silicate-Segregation)).
Volumetric dissipation $Q_{\text{diss}}$ [$\text{W/m}^3$] uses realized face fluxes $J_{\text{realized}}$ [m/s], so immobilized metal in packed cores releases zero heat:

$$Q_{\text{diss}} = \Delta\rho \, g \, J_{\text{realized}}$$

The mean temperature rise from this heat scales as $R^2$.
It is about $0.01$ to $0.1\text{ K}$ at $R = 50\text{ km}$, and a rough estimate gives of order $10\text{ K}$ at $R = 500\text{ km}$, with larger local heating along descent paths.

### Counterflow Sensible Heat

How does counterflow move enthalpy between sinking metal and rising silicate?
At each active face $(D \to R)$, descending metal $\Delta V$ at donor substep temperature $T^*_D$ [K] replaces equal silicate rising at receiver substep temperature $T^*_R$ [K], shifting sensible enthalpy by:

$$\Delta E_{D \to R} = \left[(\rho C_p)_m T^*_D - (\rho C_p)_s T^*_R\right] \Delta V$$

Each substep updates enthalpy $E_i = C_i T^*_i$ and recovers $T^*_i$ from updated sensible capacity $C_i$ [$\text{J/(m}^3\text{ K)}$].
After subcycling, $w_{c,i}/w_{\text{full},i}$ normalizes step heat using initial capacity $C_i^{\text{init}}$ and temperature $T_i^{\text{init}}$, where $w_{c,i}$ is active rock weight and $w_{\text{full},i}$ is expected full-cell weight:

$$Q_{\text{sens}, i} = \frac{w_{c,i}}{w_{\text{full},i}} C_i^{\text{init}} \frac{T^*_i - T_i^{\text{init}}}{\Delta t}$$

The thermal solver applies the combined source $Q_{\text{seg}} = Q_{\text{diss}} + Q_{\text{sens}}$:

$$\rho C_p \frac{DT}{dt} = \nabla \cdot (K \nabla T) + H_{\text{radio}} + Q_{\text{seg}} + Q_{\text{shear}}$$

**Known limit.** Counterflow conserves bare sensible enthalpy using $(\rho C_p)_m$ and $(\rho C_p)_s$, excluding metal and silicate latent heat across faces with differing melt fractions.

---

## Metal-Silicate Volatile Partitioning

During planetesimal differentiation and core formation, volatile elements (H, C, N, and S) partition between molten metallic iron alloys and silicate melts.
The thermodynamic distribution is quantified by the metal-silicate partition coefficient $D_i^{\text{met/sil}}$:

$$D_i^{\text{met/sil}} = \frac{C_i^{\text{metal}}}{C_i^{\text{silicate}}}$$

where $C_i^{\text{metal}}$ is the concentration of element $i$ in the molten metal phase [ppmw] and $C_i^{\text{silicate}}$ is the concentration dissolved in the coexisting silicate melt [ppmw].

### Empirical Parameterizations

`Erebus.jl` uses empirical fits after the cited studies for $D_i^{\text{met/sil}}(T, P, \Delta\text{IW}, w_S)$ with temperature $T$ [K], pressure $P$ [Pa], oxygen fugacity $\Delta\text{IW}$, and sulfur mass fraction $w_S$:

1. **Carbon (`model_carbon`)**:
   - `:grewal2019` (default; [Grewal et al., 2019a](../reference/bibliography.md#Metal-Silicate-Volatile-Partitioning-and-Geochemistry)) captures carbon suppression by sulfur mole fraction $X_S$:
     $$\log_{10} D_C = 1.80 + \frac{2200}{T} - 1.5 \times 10^{-8} \frac{P}{T} - 0.25 \, \Delta\text{IW} + 4.2 \ln(1 - X_S)$$
   - `:fischer2020` ([Fischer et al., 2020](../reference/bibliography.md#Metal-Silicate-Volatile-Partitioning-and-Geochemistry)):
     $$\log_{10} D_C = 1.50 + \frac{2500}{T} - 1.2 \times 10^{-8} \frac{P}{T} - 0.20 \, \Delta\text{IW}$$
   - `:constant`: Prescribed value `D_C_const`.

2. **Nitrogen (`model_nitrogen`)**:
   - `:grewal2019` (default; [Grewal et al., 2019a](../reference/bibliography.md#Metal-Silicate-Volatile-Partitioning-and-Geochemistry); [Grewal et al., 2019b](../reference/bibliography.md#Metal-Silicate-Volatile-Partitioning-and-Geochemistry)):
     $$\log_{10} D_N = 0.85 + \frac{1200}{T} - 0.25 \, \Delta\text{IW} + 0.60 \ln(1 - X_S)$$
     Because nitrogen depends less strongly on $X_S$ than carbon does, sulfur-rich core segregation lowers metallic C/N and raises silicate mantle C/N.
   - `:constant`: Prescribed value `D_N_const`.

3. **Hydrogen (`model_hydrogen`)**:
   - `:clesi2018` (default; [Clesi et al., 2018](../reference/bibliography.md#Metal-Silicate-Volatile-Partitioning-and-Geochemistry)):
     $$\log_{10} D_H = -0.80 + \frac{300}{T} + 5.0 \times 10^{-8} \frac{P}{T} + 0.05 \, \Delta\text{IW}$$
     Silicate water mass percent converts to elemental hydrogen via $1\text{ wt}\% \text{ H}_2\text{O} \approx 1118.98\text{ ppmw H}$.
   - `:constant`: Prescribed value `D_H_const`.

4. **Sulfur (`model_sulfur`)**:
   - `:boujibar2014` (default; [Boujibar et al., 2014](../reference/bibliography.md#Metal-Silicate-Volatile-Partitioning-and-Geochemistry)):
     $$\log_{10} D_S = 2.80 - \frac{800}{T} + 1.0 \times 10^{-10} P - 0.20 \, \Delta\text{IW}$$
   - `:constant`: Prescribed value `D_S_const`.

The table evaluates default partition coefficients at $T = 1500\text{ K}$ and $P = 0.1\text{ GPa}$.
At the Fe-FeS eutectic ($w_S = 0.31$, $X_S \approx 0.44$), sulfur lowers $D_C$ by a factor of $\approx 268$.

| Element | $w_S$ | $D$ at $\Delta\text{IW} = 0$ | $D$ at $\Delta\text{IW} = -2$ |
|:---|:---|:---|:---|
| C | $0$ | $\approx 1844$ | $\approx 5830$ |
| C | $0.31$ | $\approx 6.9$ | $\approx 21.8$ |
| N | $0$ | $\approx 45$ | $\approx 141$ |
| N | $0.31$ | $\approx 20$ | $\approx 64$ |
| H | any | $\approx 0.25$ | $\approx 0.20$ |
| S | any | $\approx 189$ | $\approx 475$ |

**Where this is checked.** [Core geochemistry benchmark](../reference/benchmarks/core_geochemistry.md).

### Elemental Mass Conservation

When molten metal ($F_{\text{fe}} > 0$) coexists with silicate melt ($F_m > 0$), elemental volatile mass is conserved between both phases:

$$M_{i,\text{total}} = m_{\text{sil}} C_{i,\text{sil}} + m_{\text{met}} C_{i,\text{met}}$$

where $m_{\text{sil}} = (1 - X_{\text{fe,bulk}}) \rho_{\text{sil}}$ and interacting metal mass $m_{\text{met}} = X_{\text{fe,bulk}} \rho_{\text{met}}$ use bulk metal volume fraction $X_{\text{fe,bulk}}$.
The equilibrium concentration in silicate melt is:

$$C_{i,\text{sil\_melt}}^{\text{eq}} = \frac{M_{i,\text{total}}}{m_{\text{sil}} F_{\text{melt,eff}} + D_i m_{\text{met}}}$$

where $F_{\text{melt,eff}} = \max(F_m, 10^{-4})$ prevents division by zero at incipient melting.
Equilibrium bulk concentrations follow $C_{i,\text{sil}}^{\text{eq}} = F_m C_{i,\text{sil\_melt}}^{\text{eq}}$ and $C_{i,\text{met}}^{\text{eq}} = D_i C_{i,\text{sil\_melt}}^{\text{eq}}$.

Kinetic exchange advances toward equilibrium with rate fraction $\alpha_{\text{eq}} \in [0, 1]$ (`equilibration_rate`):

$$\Delta C_{i,\text{sil}} = \alpha_{\text{eq}} \left(C_{i,\text{sil}}^{\text{eq}} - C_{i,\text{sil}}\right)$$

$$\Delta C_{i,\text{met}} = -\Delta C_{i,\text{sil}} \left(\frac{m_{\text{sil}}}{m_{\text{met}}}\right)$$

This formulation conserves total volatile mass to floating-point precision on every marker.

---

## Advective Core Segregation Transport of Volatiles

As molten metallic droplets sink under gravity, dissolved volatiles move with the metallic flux.
Within `apply_metal_segregation!`, volatile flux $F_k$ at each active cell face uses donor-cell metal volatile concentration $M_{\text{fe},k}/M_{\text{fe}}$ [ppmw]:

$$F_k = F_{\text{fe}} \cdot \left(\frac{M_{\text{fe}, k}}{M_{\text{fe}}}\right)_{\text{donor}}$$

where $F_{\text{fe}}$ is the realized metal flux, $k \in \{\text{H}, \text{C}, \text{N}, \text{S}\}$, and flux sign selects the donor cell.
Cell inventories update over west, east, north, and south faces ($F_{w,k}, F_{e,k}, F_{n,k}, F_{s,k}$) during each substep $\Delta t_{\text{sub}}$:

$$M_{\text{fe}, k}^{n+1} = M_{\text{fe}, k}^n + \Delta t_{\text{sub}} \left(F_{w, k} - F_{e, k} + F_{n, k} - F_{s, k}\right)$$

### Dynamic Sulfur Density Feedback

With `dynamic_sulfur_density = true`, each marker and cell computes alloy sulfur content $w_S = X_{\text{fe},S} \cdot 10^{-6}$ from metal sulfur concentration $X_{\text{fe},S}$ [ppmw], clamped to $[0, 0.40]$.
The local $w_S$ enters the selected `metal_density_mode` alloy law to set $\Delta\rho = \rho_{\text{metal}}(w_S) - \rho_{\text{sil,eff}}(F_m)$; under `:constant`, alloy density stays fixed.
Sulfur-poor metal is denser and segregates faster, whereas sulfur-rich metal sinks more slowly.

### Integrated Core Volatile Budgets and Core Radius

`compute_core_volatile_budgets` integrates over core markers whose bulk metal fraction reaches $\phi_{\text{cut}}$ (`[metal_partition] phi_core_threshold`, default $0.40$).
Each marker $m$ carries 3D volume $V_m = A_m (2 r_m)$, where $A_m$ is its 2D cell share [$\text{m}^2$] and $r_m$ is radial distance from center [m]:

$$M_{\text{core}, k} = \sum_{X_{\text{fe,bulk}, m} \ge \phi_{\text{cut}}} X_{\text{fe,bulk}, m} \, \rho_{\text{metal}} \, V_m \, \left(X_{\text{fe}, k, m} \cdot 10^{-6}\right)$$

$$M_{\text{core}, \text{metal}} = \sum_{X_{\text{fe,bulk}, m} \ge \phi_{\text{cut}}} X_{\text{fe,bulk}, m} \, \rho_{\text{metal}} \, V_m$$

Mean core volatile mass fractions $w_{\text{core}, k} = M_{\text{core}, k} / M_{\text{core}, \text{metal}}$ allow comparison with magmatic iron meteorites (groups IIAB, IIIAB, IVA, IVB).

When core segregation or partitioning is active, `state.accumulators.rcore` stores $(3 V_{\text{core}} / 4\pi)^{1/3}$ from core marker volume $V_{\text{core}}$.
This update runs on telemetry, snapshot, and active hydrothermal steps.
Hydrothermal convection uses `rcore` as the inner core boundary to suppress porous water flow inside segregated metal cores ($r_m \le r_{\text{core}}$).
Separately, `compute_core_radius` converts a compact metallic mass $M_{\text{core,metal}}$ at density $\rho_{\text{metal}}$ and $\phi_{\text{pack}}$ to an equivalent radius; the time loop does not call this helper.

---

## Normative Accessory Mineral Tracking and Meteorite Diagnostics

While the binary Fe-FeS eutectic is ~1261 K at low pressure, minor nickel, phosphorus, and carbon depress initial melting to ~1213 K (940 °C; [Goldstein et al., 2009](../reference/bibliography.md#Normative-Accessory-Minerals-and-Meteorite-Petrology); [Chabot and Drake, 1999](../reference/bibliography.md#Normative-Accessory-Minerals-and-Meteorite-Petrology)).
Below this temperature ($T \le T_{\text{eutectic}}$), minor and volatile elements in solid metallic iron exsolve into accessory phases:
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

Marker phase distributions are integrated over core, mantle, and crust domains.
`Erebus.jl` assigns meteorite affinity classes following [Benedix et al. (1998)](../reference/bibliography.md#Normative-Accessory-Minerals-and-Meteorite-Petrology) and [Benedix et al. (2000)](../reference/bibliography.md#Normative-Accessory-Minerals-and-Meteorite-Petrology) using internal numerical model thresholds:
- Magmatic differentiated bodies ($f_{\text{molten,core}} \ge 0.80$, $f_{\text{metal,core}} \ge 0.40$, and $f_{\text{solid,crust}} \le 0.005$): Corresponds to groups IIIAB, IVA, and IVB.
- Primitive incomplete bodies ($f_{\text{molten,core}} \le 0.60$ and $f_{\text{solid,crust}} \ge 0.01$): Corresponds to IAB complexes and winonaites.
- Transitional bodies: Intermediate segregation states.

---

## Source Code Architecture

| Physical Component | Source File | Key Functions |
|:---|:---|:---|
| Parameter definition | `src/config.jl` | `CoreFormationConfig`, `MetalPartitionConfig`, `PhaseTrackingConfig` |
| Partition coefficients | `src/physics/metal_partitioning.jl` | `compute_metal_silicate_partition_coefficient`, `compute_metal_silicate_partition_coefficients` |
| Volatile equilibration | `src/physics/metal_partitioning.jl` | `equilibrate_metal_silicate_volatiles!` |
| Accessory minerals | `src/physics/normative_minerals.jl` | `compute_troilite_stoichiometry`, `compute_schreibersite_stoichiometry`, `compute_cohenite_graphite_stoichiometry`, `compute_nitride_stoichiometry`, `compute_normative_mineral_assemblage` |
| Regional modes | `src/physics/normative_minerals.jl` | `compute_regional_mineral_modes` |
| Core budget integration | `src/physics/metal_partitioning.jl` | `compute_core_volatile_budgets` |
| Melt fraction, rheology, and velocity | `src/physics/core_dynamics.jl` | `compute_metal_melt_fraction`, `compute_silicate_suspension_viscosity`, `effective_metal_trap_fraction`, `metal_segregation_velocity`, `segregation_dissipation_heating`, `compute_core_radius` |
| Marker tracking and properties | `src/particles.jl` | `compute_marker_properties!`, `setup_marker_metal_properties`, `setup_marker_metal_volatile_properties`, `setup_marker_phase_tracking_properties`, `replenish_markers!` |
| Conservative transport | `src/numerics/segregation.jl`, `src/numerics/thermal.jl` | `apply_metal_segregation!`, `assemble_thermal_lse!` |
| Main simulation integration | `src/simulation/loop.jl` | Timestep loop sequence and checkpoint persistence |

For parameter defaults and bounds, see the [Configuration Schema Reference](../reference/config_schema.md#coreformation).
For full citations and DOIs, see the [Bibliography and Citations](../reference/bibliography.md).

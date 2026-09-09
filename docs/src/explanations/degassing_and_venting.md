# Volatile Degassing, Cold Surface Venting, and Atmospheric Escape

This section documents the physical principles, governing equations, and numerical implementations for volatile release from planetesimal interiors, cold surface venting, disk dispersal transitions, and atmospheric loss in `Erebus.jl`.

---

## Physical Overview

Early planetesimals accrete volatile-rich components, including water ice and hydrated phyllosilicates (e.g., serpentine, saponite), prior to the establishment of global surface magma oceans. Internal heating from short-lived radionuclides ($^{26}\text{Al}$ and $^{60}\text{Fe}$) drives internal temperatures upward, triggering prograde metamorphic dehydration reactions at $400\text{ to }900\text{ K}$.

When internal dehydration occurs beneath a cold, solid crust, fluid pressure accumulates in the pore network. If pore pressure exceeds confining pressure and rock tensile strength, hydraulic fractures open, routing fluids toward the surface. Unlike magma ocean degassing (where volatiles partition between liquid silicate melt and an overlying vapor envelope according to pressure-dependent solubility laws), early volatile release across a cold solid lid is governed by porous Darcy drainage, cold-trap flash sublimation, and surface boundary conductance.

---

## Cold Surface Venting Mechanics

### 1. Leaky Robin Surface Drainage

At the permeable planetesimal-nebula boundary ($r = R_{\text{planet}}$), pore fluid escapes into the ambient environment. The surface venting flux $\mathbf{q}_{\text{vent}}$ [$\text{m/s}$] is modeled as a leaky Robin boundary condition:

$$q_{\text{vent}} = \frac{k_{\text{vent}}}{\eta_f} \frac{P_f - P_{\text{vent}}}{\Delta}$$

| Symbol | Description | Units | Default |
|:---|:---|:---|:---|
| $q_{\text{vent}}$ | Outward venting discharge flux | $\text{m/s}$ | - |
| $k_{\text{vent}}$ | Surface venting boundary permeability | $\text{m}^2$ | $1.0\times 10^{-11}$ |
| $\eta_f$ | Fluid dynamic viscosity at surface temperature | $\text{Pa}\cdot\text{s}$ | $1.0\times 10^{-3}$ |
| $P_f$ | Local pore fluid pressure in boundary cell | $\text{Pa}$ | - |
| $P_{\text{vent}}$ | Effective surface venting pressure | $\text{Pa}$ | - |
| $\Delta$ | Local grid spacing ($\Delta = \sqrt{\Delta x \Delta y}$) | $\text{m}$ | - |

In the staggered finite-difference hydromechanical system, this flux manifests as an effective fluid volumetric loss rate $S_{\text{vent}}$ [$\text{s}^{-1}$] on boundary rock cells:

$$S_{\text{vent}} = \frac{C_{\text{face}}}{\Delta} \max(0, P_f - P_{\text{vent}})$$

where $C_{\text{face}} = \frac{k_{\text{vent}}}{\eta_f \Delta} \times f_{\text{conductance}}$ is the face conductance.

### 2. Venting Activation Modes and Cold Lid Rupture

`Erebus.jl` supports two venting activation modes:

1. `:darcy_sink`: Porous drainage whenever pore fluid pressure exceeds venting pressure ($P_f > P_{\text{vent}}$).
2. `:hydrofracture_gated`: Venting activates only at tensile failure zones where Terzaghi effective pressure satisfies:
   $$P_{\text{eff}} = P_t - P_f \le -\sigma_t \iff P_f \ge P_t + \sigma_t$$
   where $P_t$ is total mechanical mixture pressure and $\sigma_t$ is rock tensile strength [$\text{Pa}$]. In this mode, an intact, non-fractured crust prevents venting until internal fluid overpressure breaches the lid.

### 3. Cryogenic Pore Ice Permeability Sealing

In cold nebular or post-dispersal environments ($T_{\text{surf}} < 273.15\text{ K}$), pore fluid freezes into water ice, clogging pores and reducing matrix permeability:

$$k_{\text{eff}}(T) = k_v \cdot \left[(1 - r_{\text{min}}) \exp\left(-\frac{T_{\text{freeze}} - T}{\Delta T_{\text{seal}}}\right) + r_{\text{min}}\right] \quad (T < T_{\text{freeze}})$$

| Parameter | Description | Value | Units |
|:---|:---|:---|:---|
| $T_{\text{freeze}}$ | Liquid-solid transition temperature | $273.15$ | $\text{K}$ |
| $\Delta T_{\text{seal}}$ | Sealing temperature interval | $10.0$ | $\text{K}$ |
| $r_{\text{min}}$ | Residual cryogenic permeability ratio | $1.0\times 10^{-6}$ | - |

When overpressures breach the lid ($P_{\text{eff}} \le -\sigma_t$), macroscopic hydrofractures cut through rock and ice, providing high-permeability pathways:

$$k_{\text{frac}}(P_{\text{eff}}) = \min\left(k_{\text{max}}, k_v \left[1 + \kappa_{\text{frac}} \left(\frac{-P_{\text{eff}} - \sigma_t}{\sigma_t}\right)^\gamma\right]\right)$$

This mechanism produces episodic cryovolcanic venting cycles: internal dehydration inflates pore pressure, breaches the cold lid, rapidly discharges volatiles, drops fluid pressure, and reseals the crust.

---

## Surface Ice Cold-Trap Thermodynamics

In cold nebular or vacuum environments, ambient pressure $P_{\text{amb}}$ is extremely low ($10^{-4}\text{ to }10\text{ Pa}$). When venting fluid arrives at a cold surface ($T_{\text{surf}} < 273.16\text{ K}$), liquid water flash-freezes into ice or sublimates directly into vapor.

The equilibrium vapor pressure of water ice $P_{\text{sat,ice}}$ [$\text{Pa}$] is computed via the integrated Clausius-Clapeyron equation anchored at the water triple point:

$$P_{\text{sat,ice}}(T) = P_0 \exp\left[-\frac{L_{\text{sub}}}{R_v}\left(\frac{1}{T} - \frac{1}{T_0}\right)\right]$$

| Parameter | Description | Value | Units |
|:---|:---|:---|:---|
| $T_0$ | Water triple-point temperature | $273.16$ | $\text{K}$ |
| $P_0$ | Water triple-point vapor pressure | $611.66$ | $\text{Pa}$ |
| $L_{\text{sub}}$ | Latent heat of ice sublimation | $2.83\times 10^6$ | $\text{J/kg}$ |
| $R_v$ | Specific gas constant for water vapor | $461.5$ | $\text{J/(kg}\cdot\text{K)}$ |

The effective boundary venting pressure $P_{\text{vent}}$ represents the higher of the ambient nebular/space pressure and the local ice sublimation pressure:

$$P_{\text{vent}} = \max\left(P_{\text{amb}}, P_{\text{sat,ice}}(T_{\text{surf}})\right)$$

- At temperatures $T_{\text{surf}} < 150\text{ K}$, $P_{\text{sat,ice}} \ll 1\text{ Pa}$; the cold trap retains ice, and $P_{\text{vent}}$ is bounded below by $P_{\text{amb}}$.
- At temperatures $T_{\text{surf}} > 232\text{ K}$, $P_{\text{sat,ice}} > 10\text{ Pa}$; sublimating vapor establishes the local boundary pressure, driving free venting into the ambient medium.
- At temperatures $T_{\text{surf}} \ge 273.16\text{ K}$, saturation vapor pressure over liquid water follows the Arden Buck formulation up to water critical pressure ($22.064\text{ MPa}$ at $647.1\text{ K}$).

---

## Protoplanetary Disk Dispersal & Solar Equilibrium

Planetesimals form embedded within gas-rich protoplanetary disks, where ambient temperature $T_{\text{disk}}(t)$ and gas pressure $P_{\text{amb,disk}}$ are maintained by disk accretion heating and stellar irradiation. Over $1\text{ to }5\text{ Myr}$, photoevaporation and disk clearing disperse nebular gas.

### 1. Dispersal Weight Evolution

The transition from nebular immersion to vacuum space is parameterized by a smooth sigmoid function centered at dispersal time $t_{\text{disp}}$:

$$w_{\text{disp}}(t) = \frac{1}{1 + \exp\left[-\frac{t - t_{\text{disp}}}{\Delta t_{\text{disp}}}\right]}$$

where $t_{\text{disp}}$ is the dispersal epoch [$\text{s}$] and $\Delta t_{\text{disp}}$ is the transition timescale [$\text{s}$].

### 2. Solar Radiative Equilibrium

Following disk clearing, the ambient thermal boundary condition transitions to solar radiative equilibrium $T_{\text{eq}}$ [$\text{K}$]:

$$T_{\text{eq}} = \left[\frac{(1 - A) L_\odot}{16 \pi \sigma_{\text{SB}} d^2}\right]^{1/4}$$

| Symbol | Description | Value / Source | Units |
|:---|:---|:---|:---|
| $A$ | Surface Bond albedo | $0.06$ (default, configurable) | - |
| $L_\odot$ | Solar luminosity | $3.828\times 10^{26}$ | $\text{W}$ |
| $\sigma_{\text{SB}}$ | Stefan-Boltzmann constant | $5.670374419\times 10^{-8}$ | $\text{W/(m}^2\cdot\text{K}^4)$ |
| $d$ | Orbital distance from Sun | $d_{\text{au}} \times 1.495978707\times 10^{11}$ | $\text{m}$ |

At $2.7\text{ AU}$ with default Bond albedo $A = 0.06$, the solar equilibrium temperature is approximately $166.8\text{ K}$ ($165.0\text{ K}$ for $A = 0.1$).

### 3. Ambient State Evolution

The time-dependent ambient temperature $T_{\text{amb}}(t)$ and pressure $P_{\text{amb}}(t)$ applied to surface radiation and venting boundaries are:

$$T_{\text{amb}}(t) = (1 - w_{\text{disp}}(t)) T_{\text{disk}}(t) + w_{\text{disp}}(t) T_{\text{eq}}$$

$$P_{\text{amb}}(t) = (1 - w_{\text{disp}}(t)) P_{\text{amb,disk}} + w_{\text{disp}}(t) P_{\text{amb,space}}$$

where $P_{\text{amb,space}} = 1.0\times 10^{-4}\text{ Pa}$ represents the interplanetary space pressure floor.

---

## Coupled Porosity Drainage and Sublimation Latent Cooling

### 1. Marker Porosity Depletion

Vented fluid drains directly from Lagrangian rock markers residing within surface boundary cells:

$$\frac{D\phi_m}{Dt} = -S_{\text{vent}}(\mathbf{x}_m)$$

$$\phi_m(t + \Delta t) = \max\left(\phi_{\text{min}}, \phi_m(t) - S_{\text{vent}}(\mathbf{x}_m) \Delta t\right)$$

The cumulative vented fluid mass $M_{\text{vent}}$ [$\text{kg}$] is tracked by summing drained marker fluid:

$$\Delta M_{\text{vent}} = \sum_{m \in \text{vent}} \rho_f \cdot \left(\phi_{m,\text{old}} - \phi_{m,\text{new}}\right) \cdot V_{\text{marker}}$$

where $V_{\text{marker}} = (x_{\text{size}} y_{\text{size}}) / N_{\text{markers}}$.

### 2. Sublimation Latent Heat Sink

Venting of volatile vapor carries latent heat of phase change $L_{\text{sub}}$, extracting thermal energy from the surface rock:

$$Q_{\text{lat}} = - L_{\text{sub}} \cdot \rho_f \cdot S_{\text{vent}} \quad [\text{W/m}^3]$$

This volumetric sink is incorporated into the right-hand side of the thermal energy balance equation:

$$\rho C_p \frac{\partial T}{\partial t} = \nabla \cdot (k \nabla T) + H_R + H_A + H_S + Q_{\text{lat}}$$

where $H_R$ is radiogenic heating, $H_A$ is adiabatic work, and $H_S$ is shear heating. Sublimation cooling acts as a strong thermostatic buffer, stabilizing surface crust against runaway heating while venting continues.

---

## Multi-Species Volatile Solubility and Organic Devolatilization

During planetesimal differentiation and internal melting, volatile elements partition between crystalline minerals, liquid silicate melt, and free hydrothermal pore fluids.

### 1. Silicate Melt Water Solubility Law

At crustal pressures ($P_f \le 100\text{ MPa}$), water dissolves primarily as hydroxyl ($\text{OH}^-$) ions in liquid silicate melt. The equilibrium saturation concentration follows the square-root law of Burnham (1979) and Dixon et al. (1995):

$$w_{\text{sat}}^{\text{H}_2\text{O}} = A_s \sqrt{\max(0, P_f \times 10^{-6})} \quad [\text{wt}\%]$$

where $A_s \approx 0.40\text{ wt}\%/\text{MPa}^{0.5}$ is a representative baseline coefficient for basaltic melts and $P_f$ is pore fluid pressure [$\text{Pa}$].

### 2. Multi-Species Nitrogen Solubility Under Reducing Conditions

Nitrogen exhibits a dual dissolution mechanism in silicate melts (Libourel et al. 2003; Boulliung et al. 2020). The current parameterization is isothermal at reference magmatic conditions ($T \approx 1673\text{ K}$):

1. **Physical Molecular Dissolution ($\text{N}_2$)**: At oxidizing conditions, nitrogen dissolves as molecular dinitrogen, governed by Henry's law:
   $$w_{\text{phys}}^{\text{N}} = K_h \cdot f_{\text{N}_2} \quad [\text{ppm}]$$
   where $K_h = 0.40\text{ ppm/bar}$ is an illustrative baseline constant and $f_{\text{N}_2} = P_f \times 10^{-5}\text{ bar}$.

2. **Chemical Nitride Dissolution ($\text{N}^{3-}$)**: Under the reducing conditions typical of unoxidized planetesimal interiors ($\Delta\text{IW} \le 0$), nitrogen dissolves chemically as nitride ions:
   $$\frac{1}{2}\text{N}_2\text{(g)} + \frac{3}{2}\text{O}^{2-}\text{(melt)} \rightleftharpoons \text{N}^{3-}\text{(melt)} + \frac{3}{4}\text{O}_2\text{(g)}$$
   The chemical nitride saturation concentration scales inversely with oxygen fugacity:
   $$w_{\text{chem}}^{\text{N}} = (C_{\text{nitride}} \times 10^4) \sqrt{f_{\text{N}_2}} \cdot \left(\frac{f_{\text{O}_2}}{f_{\text{O}_2}^{\text{IW}}}\right)^{-3/4} \quad [\text{ppm}]$$
   where $C_{\text{nitride}} = 1.0\times 10^{-3}\text{ wt}\%/\text{bar}^{0.5}$ is the baseline chemical nitride capacity, and $f_{\text{O}_2} / f_{\text{O}_2}^{\text{IW}} = 10^{\Delta\text{IW}}$.

3. **Total Nitrogen Capacity**:
   $$w_{\text{total}}^{\text{N}} = w_{\text{phys}}^{\text{N}} + w_{\text{chem}}^{\text{N}} \quad [\text{ppm}]$$

When $\Delta\text{IW}$ decreases by $2$ (i.e., $f_{\text{O}_2}$ drops by two orders of magnitude, for example from $\text{IW}$ to $\text{IW}-2$), chemical nitride solubility increases by a factor of $10^{2.0 \times 0.75} = 10^{1.5} \approx 31.62$, rendering the melt a major reservoir for nitrogen storage during interior magma ocean episodes.

### 3. Iron-Wüstite Oxygen Fugacity Buffer

The 1-bar oxygen fugacity along the iron-wüstite buffer is evaluated following empirical calibrations (e.g. O'Neill 1988; Campbell et al. 2009):

$$\log_{10}(f_{\text{O}_2} [\text{bar}]) = 6.541 - \frac{28164}{T} + \Delta\text{IW}$$

where $T$ is local rock temperature in Kelvin. The default $\Delta\text{IW} = 0$ corresponds to neutral IW, while planetesimal interiors use reduced offsets such as $\Delta\text{IW} = -1.0$.

### 4. Primordial Organic Nitrogen Devolatilization

Carbonaceous chondrite parent bodies contain up to several hundred parts per million of primordial macromolecular organic nitrogen. During prograde radiogenic metamorphism ($400\text{ to }700\text{ K}$), thermal breakdown releases ammonia ($\text{NH}_3$) and molecular nitrogen ($\text{N}_2$) into the hydrothermal pore network.

This devolatilization yield $y(T) \in [0, 1]$ is parameterized as a continuous logistic transition:

$$y(T) = \frac{1}{1 + \exp\left[-\frac{T - T_{\text{devol}}}{\Delta T}\right]}$$

where $T_{\text{devol}} = 550.0\text{ K}$ is the characteristic devolatilization midpoint temperature and $\Delta T = 50.0\text{ K}$ is the thermal transition scale.

### 5. Multi-Component H-C-N-S Volatile Solubility and Saturation Limits

Beyond water and nitrogen, the complete volatile inventory includes carbon and sulfur species:

- Carbon dissolves as neutral CO and CH₄ in reduced melts (Armstrong et al. 2015; Ardia et al. 2013; Yoshioka et al. 2019) and as carbonate in oxidized melts (Dixon et al. 1995). Carbon gas fugacities are physically limited by graphite saturation ($a_{\text{C}} = 1$) via the heterogeneous equilibria of French (1966) and Holloway et al. (1992).
- Sulfur dissolves as sulfide (S²⁻) under reducing conditions and sulfate (SO₄²⁻) under oxidizing conditions (Boulliung & Wood 2022, 2023; Gaillard et al. 2022). Total dissolved sulfur is bounded by the Sulfur Content at Sulfide Saturation (O'Neill & Mavrogenes 2002; Smythe et al. 2017). When the sulfide saturation ceiling is enabled, dissolved sulfur is capped at this saturation limit.
- Homogeneous chemical equilibrium for the ten-species system (H₂, H₂O, CO, CO₂, CH₄, N₂, NH₃, H₂S, S₂, SO₂) includes graphite saturation clipping and exact elemental quadratic partitioning for nitrogen and sulfur.

The multi-species H-C-N-S solubility, speciation, and organic devolatilization routines are implemented in the physics module of `Erebus.jl` (`Erebus.Physics`). The simulation pipeline couples marker volatile exsolution to multi-species gas speciation, species-resolved atmospheric accumulation, and kinetic escape.

For benchmark comparisons across all four elemental systems, see the validation chapter [Multi-Species H-C-N-S Volatile Solubility, Speciation, and Saturation Ceilings](../validation/hcns_solubility.md).

---

## Thermodynamic Volatile Retention Floors and Low-Temperature Vent Drainage Coupling

### 1. Solid Matrix Volatile Retention Floors in Nominally Anhydrous Minerals

In classical melt solubility formulations (e.g., Burnham, 1979; Dixon et al., 1995), volatile solubility in silicate melt scales with pressure as $S_{\text{eq}} \propto \sqrt{P}$ or $S_{\text{eq}} \propto P$. Under decompression toward near-vacuum surface environments ($P \to 0$), these classical parameterizations predict that all dissolved volatiles exsolve into the fluid or vapor phase. However, laboratory analyses of natural mantle xenoliths, meteorites, and high-pressure experiments demonstrate that crystalline silicates retain trace volatile concentrations in crystal lattice defects within nominally anhydrous minerals (NAMs: olivine, pyroxenes) and refractory carbonaceous grains (Hirschmann et al., 2006; Peslier et al., 2017; Shcheka et al., 2006; Hirschmann, 2018; Li et al., 2013).

To capture this physical behavior, `Erebus.jl` incorporates thermodynamic volatile retention floors. For each volatile species (H2O, C, N, S), a temperature-dependent retention floor $C_{\text{ret}}(T)$ defines the minimum volatile concentration retained in the solid matrix. Below the reference solidus temperature $T_{\text{solidus}}$, volatiles are locked in the crystalline lattice at nominal concentration $C_{\text{floor}}$:

- **Nominally Anhydrous Minerals (NAMs) Exponential Decay** (`:nams_exponential`, default):
  $$C_{\text{ret}}(T) = \begin{cases} C_{\text{floor}}, & T \le T_{\text{solidus}} \\ C_{\text{floor}} \exp\left(-\frac{T - T_{\text{solidus}}}{\Delta T_{\text{ret}}}\right), & T > T_{\text{solidus}} \end{cases}$$
  where $\Delta T_{\text{ret}}$ is the characteristic temperature scale for melt extraction of lattice-bound volatiles.

- **Linear Melt Blend** (`:linear_melt_blend`):
  $$C_{\text{ret}}(T) = C_{\text{floor}} \left[1 - \text{clamp}\left(\frac{T - T_{\text{solidus}}}{\Delta T_{\text{ret}}}, 0, 1\right)\right]$$

- **Constant Floor** (`:constant_floor`):
  $$C_{\text{ret}}(T) = C_{\text{floor}}$$

### 2. Decompression Exsolution Clamping and Vacuum Preservation Invariant

When volatile exsolution is evaluated in marker particles, the total bulk volatile inventory $C_{\text{bulk}}$ partitions between the immovable retention floor $C_{\text{ret}}(T)$ and the mobile volatile concentration $C_{\text{mob}}$:

$$C_{\text{mob}} = \max\left(0, C_{\text{bulk}} - C_{\text{ret}}(T)\right)$$

The equilibrium melt solubility $S_{\text{eq}}(P, T)$ acts exclusively on this mobile fraction. The resulting dissolved concentration $C_{\text{sol}}$ and exsolved volatile concentration $C_{\text{exs}}$ are given by:

$$C_{\text{sol}} = C_{\text{ret}}(T) + \min\left(C_{\text{mob}}, S_{\text{eq}}\right)$$
$$C_{\text{exs}} = \max\left(0, C_{\text{mob}} - S_{\text{eq}}\right)$$

This formulation enforces two fundamental physical invariants:
1. **Vacuum Retention Floor Invariant**: As $P \to 0$ (implying $S_{\text{eq}} \to 0$), the exsolved concentration satisfies $C_{\text{exs}} \to C_{\text{mob}}$, while the dissolved concentration asymptotes to $C_{\text{sol}} \to C_{\text{ret}}(T) > 0$. Unphysical total degassing of cold crustal planetesimal rocks is strictly prevented.
2. **Sub-Floor Invariance**: When $C_{\text{bulk}} \le C_{\text{ret}}(T)$, mobile volatile concentration $C_{\text{mob}} = 0$, yielding $C_{\text{sol}} = C_{\text{bulk}}$ and $C_{\text{exs}} = 0$ identically, regardless of the local ambient pressure.

### 3. Low-Temperature Hydrothermal Venting and Mobile Volatile Drainage

During early planetesimal evolution, fluid release through porous Darcy flow and hydrofracture networks vents fluids into space or the surrounding nebula. In `Erebus.jl`, when surface venting is active ($S_{\text{vent}} > 0$) on marker particles, mobile dissolved volatiles drain at a rate proportional to the local venting volumetric sink rate:

$$\frac{d C_{\text{mob}}}{dt} = -S_{\text{vent}} \cdot C_{\text{mob}} \cdot \chi_{\text{vent}}$$

where $\chi_{\text{vent}}$ is the volatile venting extraction efficiency factor. Over a discrete computational timestep $\Delta t$, this linear kinetic equation integrates analytically:

$$C_{\text{mob}}(t + \Delta t) = C_{\text{mob}}(t) \exp\left(-S_{\text{vent}} \chi_{\text{vent}} \Delta t\right)$$

The updated bulk volatile concentration is reconstructed by adding the retained floor:

$$C_{\text{bulk}}(t + \Delta t) = C_{\text{ret}}(T) + C_{\text{mob}}(t + \Delta t)$$

The incremental mass of volatiles extracted from marker $m$ during the timestep is:

$$\Delta M_{\text{vent}, m} = \rho_m V_m \left[C_{\text{bulk}, m}(t) - C_{\text{bulk}, m}(t + \Delta t)\right]$$

These increments are scaled from 2D Cartesian cross-sectional geometry to 3D spherical geometry using the planetary metric factor $L_{\text{3D}} = 2 R_{\text{planet}}$, and accumulated into cumulative species inventories ($M_{\text{vent}}^{\text{H2O}}$, $M_{\text{vent}}^{\text{C}}$, $M_{\text{vent}}^{\text{N}}$, $M_{\text{vent}}^{\text{S}}$). This mechanism couples low-temperature hydrothermal fluid discharge directly to marker chemistry while preserving solid retention floors against exhaustive depletion.

---

## Atmospheric Accumulation and Jeans Kinetic Escape

Volatiles released through cold surface venting or magma degassing collect above the solid surface, forming a transient or steady-state atmosphere. For low-mass planetesimals, thermal effusion (Jeans escape) strips this vapor envelope to space.

### 1. Escape Velocity and Thermal Velocity

The gravitational escape velocity $v_{\text{esc}}$ [$\text{m/s}$] at planetary radius $R$ [$\text{m}$] for a body of mass $M$ [$\text{kg}$] is:

$$v_{\text{esc}} = \sqrt{\frac{2 G M}{R}}$$

The most probable Maxwellian thermal velocity $v_{\text{th}}$ [$\text{m/s}$] of a gas species with molecular mass $m$ [$\text{kg}$] at exobase temperature $T_{\text{exo}}$ [$\text{K}$] is:

$$v_{\text{th}} = \sqrt{\frac{2 k_B T_{\text{exo}}}{m}}$$

where $G = 6.67430\times 10^{-11}\text{ m}^3/(\text{kg}\cdot\text{s}^2)$ and $k_B = 1.380649\times 10^{-23}\text{ J/K}$.

### 2. Jeans Parameter and Kinetic Escape Flux

The dimensionless Jeans parameter $\lambda$ represents the ratio of gravitational binding energy to thermal kinetic energy at the exobase:

$$\lambda = \left(\frac{v_{\text{esc}}}{v_{\text{th}}}\right)^2 = \frac{G M m}{k_B T_{\text{exo}} R_{\text{exo}}}$$

The classic kinetic Jeans escape flux $\Phi_{\text{Jeans}}$ [$\text{molecules}/(\text{m}^2\cdot\text{s})$] across the exobase radius $R_{\text{exo}}$ is given by Jeans (1925):

$$\Phi_{\text{Jeans}} = \frac{n_{\text{exo}} v_{\text{th}}}{2 \sqrt{\pi}} (1 + \lambda) \exp(-\lambda)$$

where $n_{\text{exo}}$ is the number density of the species at the exobase [$\text{m}^{-3}$].

When $\lambda \ll 1$ (typical for small asteroids with $R < 100\text{ km}$ and $v_{\text{esc}} < 100\text{ m/s}$), thermal velocities exceed the escape velocity. Escape operates in the rapid effusion regime, and vented gases depart into the interplanetary medium within hours to weeks. Conversely, when $\lambda \gg 10$ (massive planetary embryos or giant planets), the exponential factor $\exp(-\lambda)$ suppresses kinetic escape, and vented volatiles accumulate into an enduring atmosphere.

### 3. Integrated Atmospheric Mass Loss Rate

Relating the exobase density $n_{\text{exo}}$ to the total atmospheric inventory $M_{\text{atm}}$ [$\text{kg}$] through the atmospheric scale height $H = k_B T / (m g)$, the total planetary mass loss rate $\dot{M}_{\text{escape}}$ [$\text{kg/s}$] can be written in linear relaxation form:

$$\dot{M}_{\text{escape}} = 4 \pi R_{\text{exo}}^2 m \Phi_{\text{Jeans}} = k_{\text{escape}} M_{\text{atm}}$$

where $k_{\text{escape}} = \frac{v_{\text{th}}}{2 \sqrt{\pi} H} (1 + \lambda) \exp(-\lambda)$ [$\text{s}^{-1}$] is the effective escape rate coefficient. The exobase density closure $M_{\text{atm}} \approx 4 \pi R_{\text{exo}}^2 \rho_{\text{exo}} H$ represents an upper bound on loss for bound atmospheres ($\lambda > 1$), because true exobase density falls below the column-averaged density.

### 4. Atmospheric Mass Conservation and Evolution

The time evolution of atmospheric mass subject to surface venting flux $\dot{M}_{\text{vent}}$ [$\text{kg/s}$] and kinetic escape is governed by:

$$\frac{d M_{\text{atm}}}{dt} = \dot{M}_{\text{vent}} - k_{\text{escape}} M_{\text{atm}}$$

In the 2D Cartesian cross-sectional domain, marker fluid drainage is evaluated per unit out-of-plane length ($[\text{kg/m}]$). Before coupling with the 3D spherical atmosphere, this 2D surface boundary flux is scaled to 3D by the ratio of spherical surface area to 2D circular boundary perimeter: $L_{\text{3D}} = A_{\text{sphere}} / P_{\text{2D}} = \frac{4 \pi R_{\text{planet}}^2}{2 \pi R_{\text{planet}}} = 2 R_{\text{planet}}$ [m]; this ensures that atmospheric inventory $M_{\text{atm}}$ has units of kilograms and surface pressure evaluates in true Pascals.

For a constant computational timestep $\Delta t$, the analytical solution yields:

$$M_{\text{atm}}(t + \Delta t) = M_{\text{atm}}(t) \exp(-x) + \frac{\dot{M}_{\text{vent}}}{k_{\text{escape}}} \left[1 - \exp(-x)\right]$$

where $x = k_{\text{escape}} \Delta t$. For small loss rates ($x < 10^{-6}$), numerical precision is preserved via a Taylor expansion:

$$\frac{1 - \exp(-x)}{k_{\text{escape}}} = \Delta t \left(1 - \frac{x}{2} + \frac{x^2}{6}\right)$$

The cumulative mass lost to space during the step is determined from exact conservation:

$$\Delta M_{\text{escaped}} = M_{\text{atm}}(t) + \dot{M}_{\text{vent}} \Delta t - M_{\text{atm}}(t + \Delta t)$$

### 5. Surface Pressure Boundary Feedback

The accumulated atmospheric mass exerts a downward hydrostatic column pressure at the planetesimal surface:

$$P_{\text{atm}} = \frac{M_{\text{atm}} g}{4 \pi R_{\text{planet}}^2}$$

where $g = G M_{\text{planet}} / R_{\text{planet}}^2$. In the coupled hydromechanical system, this atmospheric pressure contributes to the effective ambient boundary pressure:

$$P_{\text{amb,eff}} = P_{\text{amb}} + P_{\text{atm}}$$

When substantial atmospheres accumulate, $P_{\text{amb,eff}}$ opposes ongoing boiling and venting, naturally throttling further surface volatile discharge. In the explicit time-advancement scheme of `Erebus.jl`, surface atmospheric pressure $P_{\text{atm}}$ is evaluated from the previous step atmospheric inventory, lagging the hydromechanical solve by one timestep in an operator-split fashion. Within the 2D Stokes-Darcy simulation loop, marker fluid drainage and volatile exsolution feed species-resolved atmospheric inventories ($M_{\text{atm}}^{\text{species}}$) for the 10-species gas mixture. Multi-species kinetic escape is evaluated for each volatile species using its molecular mass, with a cubic Hermite transition between sound-speed blow-off and kinetic Jeans effusion, tracking cumulative escaped masses ($M_{\text{escaped}}^{\text{species}}$) across timesteps.


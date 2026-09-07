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
- At temperatures $T_{\text{surf}} > 200\text{ K}$, $P_{\text{sat,ice}} > 10\text{ Pa}$; sublimating vapor establishes the local boundary pressure, driving free venting into the ambient medium.

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


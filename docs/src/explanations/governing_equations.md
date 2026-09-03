# Governing Equations

This section states the complete mathematical system of partial differential equations and constitutive relations solved in `Erebus.jl`.

---

## Hydromechanical Conservation Equations

The mechanical response of a two-phase deformable porous medium saturated with fluid is governed by the coupled Stokes-Darcy equations (Gerya, 2019).

### 1. Solid Momentum Conservation
In the creeping-flow limit (negligible inertia, $\text{Re} \ll 1$), momentum conservation for the combined solid-fluid mixture is:

$$\frac{\partial \sigma_{ij}'}{\partial x_j} - \frac{\partial P_t}{\partial x_i} + \rho_{\text{total}} g_i = 0$$

| Symbol | Description | Units |
|:---|:---|:---|
| $\sigma_{ij}'$ | Deviatoric stress tensor of the solid matrix | $\text{Pa}$ |
| $P_t$ | Total mechanical mixture pressure | $\text{Pa}$ |
| $\rho_{\text{total}}$ | Bulk mixture density: $(1 - \phi)\rho_s + \phi\rho_f$ | $\text{kg/m}^3$ |
| $g_i$ | Gravitational acceleration vector | $\text{m/s}^2$ |

### 2. Solid Mass Continuity
Conservation of solid mass accounting for volumetric compaction is:

$$\nabla \cdot \mathbf{v}_s = - \frac{D \ln(1 - \phi)}{Dt}$$

where $\mathbf{v}_s$ is solid matrix velocity [$\text{m/s}$], $\phi$ is porosity [-], and $D/Dt$ is the material derivative.

### 3. Darcy Pore Fluid Flux and Thermal Buoyancy
Percolation of pore fluid relative to the solid matrix follows Darcy's law:

$$\mathbf{q}_D = - \frac{k_\phi}{\eta_f} \left(\nabla P_f - \rho_f(T) \mathbf{g}\right)$$

| Symbol | Description | Units |
|:---|:---|:---|
| $\mathbf{q}_D$ | Darcy fluid discharge flux: $\phi (\mathbf{v}_f - \mathbf{v}_s)$ | $\text{m/s}$ |
| $k_\phi$ | Porosity-dependent permeability | $\text{m}^2$ |
| $\eta_f$ | Dynamic fluid viscosity | $\text{Pa}\cdot\text{s}$ |
| $P_f$ | Pore fluid pressure | $\text{Pa}$ |
| $\rho_f(T)$ | Temperature-dependent pore fluid density | $\text{kg/m}^3$ |

Fluid density incorporates volumetric thermal expansion $\alpha_f$ [$\text{K}^{-1}$] above the reference melting temperature $T_{\text{melt}}$:

$$\rho_f(T) = \rho_{f0} \max\left(0.1, 1.0 - \alpha_f (T - T_{\text{melt}})\right) \quad (T > T_{\text{melt}})$$

The temperature gradient creates buoyancy forces $(\rho_f(T) - \rho_{f0})\mathbf{g}$ that drive hydrothermal convection through permeable silicate crust. Fluid density is referenced to the melting point $T_{\text{melt}} = 273.0\text{ K}$ ($\rho_{f0} = 1000\text{ kg/m}^3$), while thermal expansion is applied to the pore fluid density $\rho_f(T)$ to drive relative Darcy percolation; solid grain density $\rho_s$ remains compositionally dependent on hydration state $X_W^s$.

Dynamic fluid viscosity $\eta_f(T)$ decreases with temperature following an Arrhenius relation for liquid water:

$$\eta_f(T) = \begin{cases} \eta_{\text{ice}} & T \le T_{\text{melt}} \\ \max\left(\eta_{\text{min}}, \min\left(\eta_{\text{max}}, \eta_{f0} \exp\left[\frac{E_a}{R} \left(\frac{1}{T} - \frac{1}{T_0}\right)\right]\right)\right) & T > T_{\text{melt}} \end{cases}$$

where viscosity is referenced to $T_0 = 293.15\text{ K}$ ($20^\circ\text{C}$, $\eta_{f0} = 1.0\times 10^{-3}\text{ Pa}\cdot\text{s}$) with activation energy $E_a = 15.0\text{ kJ/mol}$ and gas constant $R$. This relation reproduces liquid water viscosity within $12\%$ across $273\text{--}373\text{ K}$ and provides an effective approximation across the liquid hydrothermal range ($273\text{--}600\text{ K}$), where fluid mobility $k_\phi / \eta_f(T)$ increases by a factor of $5\times\text{--}20\times$ (clamped to a maximum factor of $\eta_{f0}/\eta_{\text{min}} = 100\times$).

### 4. Coupled Fluid Mass Continuity
Conservation of pore fluid mass with fluid compressibility is:

$$\nabla \cdot \mathbf{q}_D = \frac{D \ln \phi}{Dt} + \phi \beta_f \frac{D P_f}{Dt}$$

where $\beta_f$ is the fluid compressibility [$\text{Pa}^{-1}$]. In the discrete numerical system, the total mechanical continuity equation also incorporates the volumetric expansion rate $\text{DMP}$ [$1/\text{s}$] from mineral dehydration reactions.

---

## Poroelastic Constitutive Equations

`Erebus.jl` couples solid matrix and fluid compressibilities using linear poroelasticity theory:

### 1. Drained and Bulk Compressibility
The drained compressibility $\beta_d$ represents skeletal compliance under drained conditions ($P_f = 0$):

$$\beta_d = \frac{\beta_\phi + \beta_s}{1 - \phi}$$

where $\beta_\phi = \phi / G$ is the elastic pore compliance determined by porosity $\phi$ and the solid shear elastic modulus $G$ [$\text{Pa}$], and $\beta_s$ is the intrinsic solid grain compressibility [$\text{Pa}^{-1}$].

### 2. Biot-Willis Coefficient
The Biot-Willis coefficient $K_{\text{BW}}$ scales the coupling between total mechanical pressure and pore fluid pressure:

$$K_{\text{BW}} = 1 - \frac{\beta_s}{\beta_d}$$

In the incompressible solid grain limit ($\beta_s \to 0$), $K_{\text{BW}} \to 1$.

### 3. Skempton Coefficient
The Skempton coefficient $B$ quantifies the undrained pore pressure response to changes in total confining pressure:

$$B = \frac{\beta_d - \beta_s}{\beta_d - \beta_s + \phi (\beta_f - \beta_s)}$$

---

## Terzaghi Effective Stress and Failure Criteria

### 1. Effective Stress Decomposition
The Terzaghi effective stress tensor $\boldsymbol{\sigma}'_{\text{eff}}$ governs matrix deformation:

$$\sigma_{ij}'^{\text{eff}} = \sigma_{ij}' - P_{\text{eff}} \delta_{ij}$$

where $P_{\text{eff}} = P_t - P_f$ is the Terzaghi effective pressure, related to the solid skeleton pressure $P_s$ and total mixture pressure $P_t = (1 - \phi) P_s + \phi P_f$ by $P_{\text{eff}} = (1 - \phi)(P_s - P_f)$.

### 2. Mohr-Coulomb Plasticity and Tensile Failure
The yield stress $\sigma_{\text{yield}}$ is evaluated using a linear Drucker-Prager / Mohr-Coulomb criterion:

$$\sigma_{\text{yield}} = C + \mu P_{\text{eff}}$$

where $C$ is cohesion [$\text{Pa}$] and $\mu$ is the internal friction coefficient. If the effective pressure drops below the tensile strength $\sigma_t$:

$$P_{\text{eff}} \le -\sigma_t$$

tensile hydrofracturing occurs, limiting the sustainable pore overpressure.

---

## Thermal Energy Conservation

Heat transport across the planetesimal is governed by the energy conservation equation:

$$\rho_{\text{total}} c_p \frac{\partial T}{\partial t} = \nabla \cdot (k \nabla T) + Q_{\text{radiogenic}} + Q_{\text{adiabatic}} + Q_{\text{shear}} + \text{DHP}$$

| Symbol | Description | Units |
|:---|:---|:---|
| $c_p$ | Specific heat capacity (with melting latent heat buffer) | $\text{J}/(\text{kg}\cdot\text{K})$ |
| $k$ | Bulk thermal conductivity | $\text{W}/(\text{m}\cdot\text{K})$ |
| $Q_{\text{adiabatic}}$ | Adiabatic compression/decompression heating | $\text{W/m}^3$ |
| $Q_{\text{shear}}$ | Viscous shear and compaction dissipation | $\text{W/m}^3$ |
| $\text{DHP}$ | Mineral dehydration reaction enthalpy transfer rate | $\text{W/m}^3$ |
| $Q_{\text{radiogenic}}$ | Total radiogenic heat production rate: $Q_{\text{al}} + Q_{\text{fe}}$ | $\text{W/m}^3$ |

### Radiogenic Decay Kinetics
Volumetric heat production from $^{26}\text{Al}$ and $^{60}\text{Fe}$ decay is computed with the decay rate constant $\lambda = 1/\tau = \ln(2)/t_{1/2}$:

$$Q_{\text{al}}(t) = f_{\text{al}} \left(\frac{^{26}\text{Al}}{^{27}\text{Al}}\right)_0 E_{\text{al}} \frac{1}{\tau_{\text{al}}} \exp\left(-\frac{t}{\tau_{\text{al}}}\right) \rho_s$$

$$Q_{\text{fe}}(t) = f_{\text{fe}} \left(\frac{^{60}\text{Fe}}{^{56}\text{Fe}}\right)_0 E_{\text{fe}} \frac{1}{\tau_{\text{fe}}} \exp\left(-\frac{t}{\tau_{\text{fe}}}\right) \rho_f$$

where $\tau$ is the mean lifetime of each radioactive isotope, and $\rho_s$, $\rho_f$ are solid and fluid densities.

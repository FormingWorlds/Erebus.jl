# Governing Equations

This section states the complete mathematical system of partial differential equations and constitutive relations solved in `Erebus.jl`.

---

## Hydromechanical Conservation Equations

The mechanical response of a two-phase deformable porous medium saturated with fluid is governed by the coupled Stokes-Darcy equations (Gerya, 2019).

### 1. Solid Momentum Conservation
In the creeping-flow limit (negligible inertia, $\text{Re} \ll 1$), momentum conservation for the combined solid-fluid mixture is:

$$\frac{\partial \sigma_{ij}'}{\partial x_j} - \frac{\partial P_t}{\partial x_i} + \rho_{\text{total}} g_i = 0$$

where:
- $\sigma_{ij}'$ is the deviatoric stress tensor of the solid matrix [$\text{Pa}$]
- $P_t$ is the total mechanical pressure [$\text{Pa}$]
- $\rho_{\text{total}} = (1 - \phi) \rho_s + \phi \rho_f$ is the mixture bulk density [$\text{kg/m}^3$]
- $g_i$ is the gravitational acceleration vector [$\text{m/s}^2$]

### 2. Solid Mass Continuity
Conservation of solid mass accounting for volumetric compaction is:

$$\nabla \cdot \mathbf{v}_s = - \frac{D \ln(1 - \phi)}{Dt}$$

where $\mathbf{v}_s$ is the solid matrix velocity [$\text{m/s}$], $\phi$ is the porosity [-], and $D/Dt = \partial/\partial t + \mathbf{v}_s \cdot \nabla$ represents the material derivative following the solid skeleton.

### 3. Darcy Pore Fluid Flux
Percolation of the pore fluid relative to the solid matrix follows Darcy's law:

$$\mathbf{q}_D = - \frac{k_\phi}{\eta_f} \left(\nabla P_f - \rho_f \mathbf{g}\right)$$

where:
- $\mathbf{q}_D = \phi (\mathbf{v}_f - \mathbf{v}_s)$ is the Darcy fluid discharge flux [$\text{m/s}$]
- $k_\phi$ is the porosity-dependent permeability [$\text{m}^2$]
- $\eta_f$ is the dynamic fluid viscosity [$\text{Pa}\cdot\text{s}$]
- $P_f$ is the pore fluid pressure [$\text{Pa}$]
- $\rho_f$ is the pore fluid density [$\text{kg/m}^3$]

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

where:
- $c_p$ is the specific heat capacity [$\text{J}/(\text{kg}\cdot\text{K})$], with water ice melting latent heat ($L^f$) modeled as an effective heat capacity peak near $273\text{ K}$
- $k$ is the thermal conductivity [$\text{W}/(\text{m}\cdot\text{K})$]
- $Q_{\text{adiabatic}}$ is adiabatic compression/decompression heating [$\text{W/m}^3$]
- $Q_{\text{shear}}$ is viscous shear and compaction dissipation [$\text{W/m}^3$]
- $\text{DHP}$ is the mineral dehydration reaction enthalpy transfer rate [$\text{W/m}^3$]
- $Q_{\text{radiogenic}}$ is the total radiogenic heat production rate [$\text{W/m}^3$].

### Radiogenic Decay Kinetics
Volumetric heat production from $^{26}\text{Al}$ and $^{60}\text{Fe}$ decay is computed with the decay rate constant $\lambda = 1/\tau = \ln(2)/t_{1/2}$:

$$Q_{\text{al}}(t) = f_{\text{al}} \left(\frac{^{26}\text{Al}}{^{27}\text{Al}}\right)_0 E_{\text{al}} \frac{1}{\tau_{\text{al}}} \exp\left(-\frac{t}{\tau_{\text{al}}}\right) \rho_s$$

$$Q_{\text{fe}}(t) = f_{\text{fe}} \left(\frac{^{60}\text{Fe}}{^{56}\text{Fe}}\right)_0 E_{\text{fe}} \frac{1}{\tau_{\text{fe}}} \exp\left(-\frac{t}{\tau_{\text{fe}}}\right) \rho_f$$

where $\tau$ is the mean lifetime of each radioactive isotope, and $\rho_s$, $\rho_f$ are solid and fluid densities.

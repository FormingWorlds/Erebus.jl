# Crustal Sill Cooling and Magma-Hydrothermal Coupling Benchmarks

This page documents the verification, and analytical validation, of crustal sill cooling, sensible heat enthalpy advection, latent heat crystallization, and hydrothermal convection coupling in `Erebus.jl`.

---

## 1. Jaeger (1957) 1D Analytical Sill Cooling Benchmark

### Physical Setup and Governing Equation

We model the conductive dissipation of heat from an intrusive igneous sheet, or sill, of thickness $2b$, emplaced instantaneously at initial temperature $T_0$, into cold country rock at initial ambient temperature $T_c$.

Heat transfer in the solid host, and in the cooling intrusion, follows the 1D transient diffusion equation:

$$\frac{\partial T}{\partial t} = \kappa \frac{\partial^2 T}{\partial y^2}$$

where $\kappa = k / (\rho c_p)$ is the constant thermal diffusivity [$\text{m}^2/\text{s}$].

The exact analytical solution (Jaeger 1957), as a function of distance $y$ from the sill center line ($y = 0$) at time $t$, is:

$$T(y, t) = T_c + \frac{T_0 - T_c}{2} \left[ \text{erf}\left(\frac{b - y}{2\sqrt{\kappa t}}\right) + \text{erf}\left(\frac{b + y}{2\sqrt{\kappa t}}\right) \right]$$

At the center of the sill, where $y = 0$:

$$T(0, t) = T_c + (T_0 - T_c) \text{erf}\left(\frac{b}{2\sqrt{\kappa t}}\right)$$

At the intrusion contact, where $y = \pm b$ as $t \to 0$:

$$T(\pm b, 0^+) = \frac{T_0 + T_c}{2}$$

### Benchmark Parameters

| Parameter | Symbol | Value | Units |
| :--- | :--- | :--- | :--- |
| Sill Half-Thickness | $b$ | 50.0 | m |
| Initial Magma Temperature | $T_0$ | 1400.0 | K |
| Country Rock Temperature | $T_c$ | 400.0 | K |
| Rock Density | $\rho$ | 2800.0 | $\text{kg/m}^3$ |
| Specific Heat Capacity | $c_p$ | 1000.0 | $\text{J/(kg K)}$ |
| Thermal Conductivity | $k$ | 2.8 | $\text{W/(m K)}$ |
| Thermal Diffusivity | $\kappa$ | $1.0 \times 10^{-6}$ | $\text{m}^2/\text{s}$ |
| Simulation Domain Width | $L$ | 400.0 | m |
| Target Cooling Time | $t_{\text{total}}$ | $5.0 \times 10^8$ | s ($\approx 15.8$ yr) |

### Numerical Verification

The benchmark drives the production `assemble_thermal_lse!` solver over a 1D column with backward Euler time integration on a mesh of $N = 201$ nodes ($\Delta y = 2.0$ m).
At time $t = 5.0 \times 10^8$ s, the characteristic thermal diffusion length, $2\sqrt{\kappa t} \approx 44.7$ m, spans 22.4 grid cells.
The test asserts:
- Relative $L_2$ error norm, $\|T_{\text{num}} - T_{\text{ana}}\|_2 / \|T_{\text{ana}}\|_2$, is less than $0.1\%$ ($< 1.0 \times 10^{-3}$).
- Maximum point-wise absolute error, $\max |T_{\text{num}} - T_{\text{ana}}|$, is less than $1.0$ K throughout the domain.
- Centerline temperature agreement at $y = 0$ is within 1.0 K (numerical difference $< 0.75$ K).

---

## 2. Stefan Latent Heat Buffering Benchmark

### Solidification Phase Change Formulation

When ponded magma cools between the liquidus, $T_{\text{liq}}$, and solidus, $T_{\text{sol}}$, crystallization releases latent heat of fusion, $L_m$ [$\text{J/kg}$].
The phase change is characterized by the dimensionless Stefan number:

$$\text{Ste} = \frac{c_p (T_{\text{liq}} - T_{\text{sol}})}{L_m}$$

In the mushy crystallization interval, latent heat release increases the apparent heat capacity:

$$c_{p,\text{apparent}} = c_p + \frac{L_m}{T_{\text{liq}} - T_{\text{sol}}} = c_p \left( 1 + \frac{1}{\text{Ste}} \right)$$

### Verification Results

The test evaluates `rhocp_apparent_silicate` directly across subsolidus, mushy, and superliquidus intervals with silicate parameters ($L_m = 4.0 \times 10^5$ J/kg, $c_p = 1000$ J/(kg K), $\Delta T = 200$ K):
- Dimensionless Stefan number $\text{Ste} = 0.5$.
- Inside the mushy interval, the apparent heat capacity buffering factor is $1 + \text{Ste}^{-1} = 3.0$.
- Outside the mushy interval, the apparent heat capacity equals baseline sensible heat capacity ($1.0 \times \rho c_p$).

---

## 3. Sensible Heat Enthalpy Transport and Energy Conservation

### Pairwise Flux Discretization

Melt segregation transports sensible heat advectively through cell faces using relative temperature differences between donor and receiver cells:

$$\mathbf{H}_{\text{flux}} = \mathbf{q}_m \rho_m c_{p,m} (T_{\text{donor}} - T_{\text{rec}})$$

This advective formulation computes the sensible transport term $-\rho_m c_{p,m} (\mathbf{q}_m \cdot \nabla T)$ without spurious divergence artifacts on fixed-capacity Eulerian grids.
The heating increment is deposited directly onto interior staggered grid nodes $(i+1, j+1)$ to prevent boundary ghost node leakage.

### Conservation Verification

Tests in `test/test_sill_cooling.jl` verify:
1. **Isothermal Zero Sensible Heating**: In an isothermal domain ($\nabla T = 0$), sensible heating is identically zero ($Q_{\text{seg}} \equiv 0$ to machine precision, $< 10^{-12} \text{ W/m}^3$), preventing unphysical heating from non-solenoidal melt ponding.
2. **Ghost-Free Nodal Deposition**: All sensible heating is deposited onto interior nodes ($2 \le i \le Ny$, $2 \le j \le Nx$), leaving boundary ghost rims strictly zero.
3. **Thermal Stratification Extraction**: Upward segregation extracts heat from the deep, hot interior ($Q_{\text{sens}} < 0$), and deposits it in the shallow, cooler lithosphere ($Q_{\text{sens}} > 0$).

---

## 4. Hydrothermal Convection and Contact Metamorphism Coupling

### Convection Enhancement Above Emplaced Sills

Cooling crustal sills heat overlying porous rock.
When rock temperature exceeds $T_{\text{surface}} + \Delta T_{\text{min}}$ in permeable aquifers ($\phi > \phi_{\text{start}}$), porous Rayleigh-Darcy convection triggers ($Ra_m > Ra_{m,\text{crit}} = 4\pi^2$):

$$k_{\text{eff}} = \text{Nu} \cdot k_{\text{cond}} > k_{\text{cond}}$$

Enhanced effective thermal conductivity accelerates heat removal from the underlying sill, relative to conductive diffusion.
The effective convective layer thickness, $H_{\text{eff}}$, scales the local Rayleigh number, providing consistent convective closure for shallow intrusive sheets.

### Contact Metamorphism and Dehydration Aureoles

In hydrated host rock, such as serpentine or chlorite, conductive heat from the sill drives host rock temperatures above the thermodynamic dehydration equilibrium ($T > T_{\text{eq}} = \Delta H / \Delta S$).
Contact metamorphic dehydration produces water:

$$\text{DQPF} = \frac{\Gamma_{\text{water}}}{\rho_f} > 0$$

In low-permeability rock, this dehydration fluid source creates dynamic pore fluid overpressures ($\Delta P = P_f - P_r > 50\text{ kPa}$), driving hydrofracturing, and fluid venting, into overlying porous layers.

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

The simulation solves the 1D diffusion equation, using a second-order Crank-Nicolson implicit scheme, on a mesh of $N = 201$ nodes ($\Delta y = 2.0$ m).
At time $t = 5.0 \times 10^8$ s, the characteristic thermal diffusion length, $2\sqrt{\kappa t} \approx 44.7$ m, spans more than 22 grid points.
The test asserts:
- Relative $L_2$ error norm, $\|T_{\text{num}} - T_{\text{ana}}\|_2 / \|T_{\text{ana}}\|_2$, is less than $0.1\%$ ($< 1.0 \times 10^{-3}$).
- Maximum point-wise absolute error, $\max |T_{\text{num}} - T_{\text{ana}}|$, is less than $1.0$ K throughout the domain.
- Agreement at the sill center, $T(0, t_{\text{total}})$, is better than 0.5 K.

---

## 2. Stefan Latent Heat Buffering Benchmark

### Solidification Phase Change Formulation

When ponded magma cools between the liquidus, $T_{\text{liq}}$, and solidus, $T_{\text{sol}}$, crystallization releases latent heat of fusion, $L_m$ [$\text{J/kg}$].
The phase change is characterized by the dimensionless Stefan number:

$$\text{Ste} = \frac{c_p (T_{\text{liq}} - T_{\text{sol}})}{L_m}$$

In the mushy crystallization interval, latent heat release increases the apparent heat capacity:

$$c_{p,\text{apparent}} = c_p + \frac{L_m}{T_{\text{liq}} - T_{\text{sol}}} = c_p \left( 1 + \frac{1}{\text{Ste}} \right)$$

### Verification Results

With typical silicate melt parameters ($L_m = 4.0 \times 10^5$ J/kg, $c_p = 1000$ J/(kg K), $\Delta T = 200$ K):
- $\text{Ste} = 0.5$.
- The effective thermal buffering factor is $1 + \text{Ste}^{-1} = 3.0$.
- The duration required to cool through the crystallization interval increases by a factor of 3.0, relative to sensible cooling alone.

---

## 3. Sensible Heat Enthalpy Transport and Energy Conservation

### Pairwise Flux Discretization

Melt segregation transports sensible enthalpy advectively through cell faces:

$$\mathbf{H}_{\text{flux}} = \mathbf{q}_m \rho_m c_{p,m} T_{\text{donor}}$$

The volumetric net heating rate, deposited into the energy solver, is:

$$Q_{\text{sens}} = -\nabla \cdot \mathbf{H}_{\text{flux}}$$

### Conservation Verification

Tests in `test/test_sill_cooling.jl` verify:
1. **Machine-Precision Telescoping Cancellation**: Every face enthalpy flux is subtracted from the donor cell, and added to the receiver cell. The domain-integrated net sensible energy satisfies $\sum \Delta H = 0.0$ to machine precision.
2. **Nodal Source Balance**: Interpolation to staggered grid nodes preserves zero net energy injection: $\int Q_{\text{seg}} \, dV = 0.0$, remaining below $10^{-10}$ W.
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

In hydrated host rock, such as serpentine or chlorite, conductive heat from the sill drives host rock temperatures above the devolatilization limit ($T > 650$ K).
Contact metamorphic dehydration produces water:

$$\text{DQPF} = \frac{\Gamma_{\text{water}}}{\rho_f} > 0$$

In low-permeability rock, this dehydration fluid source creates dynamic pore fluid overpressures ($\Delta P = P_f - P_r > 50\text{ kPa}$), driving hydrofracturing, and fluid venting, into overlying porous layers.

# Fluid Viscosity and Phase Transition Validation

This module validates temperature-dependent fluid viscosity and phase state transitions for pore water.

---

## Governing Formulation

Above the melting point ($T \ge T_{\text{melt}} = 273.0\text{ K}$), liquid water viscosity follows the Arrhenius activation law normalized at reference temperature $T_0 = 293.15\text{ K}$:

$$\eta_f(T) = \eta_0 \exp\left[ \frac{E_a}{R} \left( \frac{1}{T} - \frac{1}{T_0} \right) \right]$$

Below the melting point ($T < T_{\text{melt}}$), pore water freezes into solid ice, switching to high-viscosity solid rheology:

$$\eta_f(T) = \eta_{\text{ice}} \approx 10^{12}\text{ Pa}\cdot\text{s}$$

---

## Literature Anchors

- **Hubmann, B. (2022)**. *Hydrology of Planetesimals*. Master's thesis, ETH Zurich.  
  [https://doi.org/10.5281/zenodo.7058229](https://doi.org/10.5281/zenodo.7058229) (Equation 2.10)
- **Gerya, T. (2019)**. *Introduction to Numerical Geodynamic Modelling* (2nd ed.). Cambridge University Press.  
  [https://doi.org/10.1017/9781316534243](https://doi.org/10.1017/9781316534243)

---

## Invariants and Physical Limits

1. **Liquid Monotonicity**:
   In the liquid phase ($T \ge T_{\text{melt}}$), viscosity decreases strictly monotonically with temperature: $\partial \eta_f / \partial T < 0$.
2. **Phase Boundary Discontinuity**:
   Crossing $T = T_{\text{melt}}$ produces a sharp, physically bounded transition from liquid Darcy flow ($\eta_f \sim 10^{-3}\text{ Pa}\cdot\text{s}$) to immobile frozen pore ice ($\eta_f \approx 10^{12}\text{ Pa}\cdot\text{s}$).
3. **Strict Positivity**:
   $\eta_f(T) > 0$ for all physical temperatures $T > 0\text{ K}$.
4. **Physical Bounds**:
   Liquid viscosity is clamped within physical bounds: $\eta_{\text{min}} \le \eta_f \le \eta_{\text{max}}$.
5. **Phase Limiting Behavior**:
   Non-positive absolute temperatures ($T \le 0\text{ K}$) or non-finite values evaluate to the frozen ice regime ($\eta_{\text{ice}} = 10^{12}\text{ Pa}\cdot\text{s}$).

---

## Parameterization Benchmarking

The temperature-dependent fluid viscosity formulation approximates liquid water behavior under planetesimal hydrothermal conditions:

$$\eta_f(T) = \eta_{f0} \exp\left[\frac{E_a}{R} \left(\frac{1}{T} - \frac{1}{T_0}\right)\right]$$

with reference viscosity $\eta_{f0} = 1.0\times 10^{-3}\text{ Pa}\cdot\text{s}$ at $T_0 = 293.15\text{ K}$ ($20^\circ\text{C}$) and activation energy $E_a = 15.0\text{ kJ/mol}$.

This single-activation energy relation closely tracks experimental liquid water data (IAPWS / NIST standard tables) within about $13\%$ in the sub-boiling regime $T \in [273, 373]\text{ K}$. At higher temperatures ($T \approx 470\text{ to }570\text{ K}$), liquid water curvature follows a Vogel-Fulcher-Tammann profile, where the single Arrhenius fit underestimates viscosity by $\approx 29\text{ to }43\%$ before entering the supercritical regime ($T > 647\text{ K}$).

![Temperature-Dependent Fluid Viscosity Benchmarking](../assets/fluid_viscosity_temperature.png)

*Figure 1: Benchmarking of temperature-dependent pore fluid viscosity $\eta_f(T)$ in Erebus. (a) Dynamic fluid viscosity over the range $T \in [270, 650]\text{ K}$ on a logarithmic scale, comparing the default Arrhenius model ($E_a = 15.0\text{ kJ/mol}$, blue curve) against experimental liquid water measurements from IAPWS/NIST standards (red circles). (b) Hydrothermal Darcy mobility enhancement factor $\eta_{f0} / \eta_f(T)$ illustrating the $5\times\text{ to }24\times$ increase in fluid percolation speed as interior temperatures rise in liquid hydrothermal conditions ($273\text{ to }600\text{ K}$).*

---

## Verification Test Suite

- `test/test_physics.jl`:
  - `@testset "etatotal_rocks(): phase transition and viscosity bounds"`
  - `@testset "temperature-dependent fluid viscosity and compute_fluid_viscosity()"`

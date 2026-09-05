# Fluid Viscosity and Phase Transition Validation

This module validates temperature-dependent fluid viscosity and phase state transitions for pore water.

---

## Governing Formulation

Above the melting point ($T \ge T_{\text{melt}} = 273.15\text{ K}$), liquid water viscosity follows the Arrhenius activation law:

$$\eta_f(T) = \eta_0 \exp\left( \frac{E_a}{R T} \right)$$

Below the melting point ($T < T_{\text{melt}}$), pore water freezes into solid ice, switching to high-viscosity solid rheology:

$$\eta_f(T) = \eta_{\text{ice}} \approx 10^{14}\text{ Pa}\cdot\text{s}$$

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
   Crossing $T = T_{\text{melt}}$ produces a sharp, physically bounded transition from liquid Darcy flow ($\eta_f \sim 10^{-3}\text{ Pa}\cdot\text{s}$) to immobile frozen pore ice ($\eta_f \sim 10^{14}\text{ Pa}\cdot\text{s}$).
3. **Strict Positivity**:
   $\eta_f(T) > 0$ for all physical temperatures $T > 0\text{ K}$.
4. **Physical Bounds**:
   Liquid viscosity is clamped within physical bounds: $\eta_{\text{min}} \le \eta_f \le \eta_{\text{max}}$.
5. **Error Contract**:
   Passing non-positive absolute temperatures ($T \le 0\text{ K}$) throws `DomainError`.

---

## Verification Test Suite

- `test/test_physics.jl`: `@testset "etatotal_rocks()"`

# Thermal Conduction and Metric Geometry Validation

This module validates the effective thermal conductivity of porous rock-fluid mixtures and the geometric metric divergence term for spherical heat conduction in 2D Cartesian coordinates.

---

## 1. Bulk Thermal Conductivity of Porous Mixtures

### Governing Formulation
The effective bulk thermal conductivity $k_{\text{total}}$ of a porous rock matrix saturated with fluid (water or ice) is computed using the quadratic mixture formulation:

$$k_{\text{total}} = \sqrt{\frac{k_s k_f}{2} + \frac{[k_s(3\phi - 2) + k_f(1 - 3\phi)]^2}{16}} - \frac{k_s(3\phi - 2) + k_f(1 - 3\phi)}{4}$$

where $k_s$ is matrix solid conductivity, $k_f$ is pore fluid conductivity, and $\phi$ is porosity.

### Literature Anchors
- **Gerya, T. (2019)**. *Introduction to Numerical Geodynamic Modelling* (2nd ed.). Cambridge University Press.  
  [https://doi.org/10.1017/9781316534243](https://doi.org/10.1017/9781316534243) (Section 10.3)
- **Hubmann, B. (2022)**. *Hydrology of Planetesimals*. Master's thesis, ETH Zurich.  
  [https://doi.org/10.5281/zenodo.7058229](https://doi.org/10.5281/zenodo.7058229) (Equation 2.14)

### Invariants and Limits
1. **Solid Matrix Limit**: As $\phi \to 0$, $k_{\text{total}} \to k_s$.
2. **Pure Fluid Limit**: As $\phi \to 1$, $k_{\text{total}} \to k_f$.
3. **Monotonicity**: For $k_s > k_f$, $\partial k_{\text{total}} / \partial \phi < 0$.
4. **Physical Bounding**: $\min(k_s, k_f) \le k_{\text{total}} \le \max(k_s, k_f)$ for all $\phi \in [0, 1]$.
5. **Error Contract**: Passing $\phi < 0$, $\phi > 1$, $k_s < 0$, or $k_f < 0$ throws `DomainError`.

### Verification Test Suite
- `test/test_physics.jl`: `@testset "ktotal(): physical limits and non-linear mixture invariants"`

---

## 2. Spherical Geometric Metric Divergence (2D Cartesian)

### Governing Formulation
Radial heat conduction in a 3D spherically symmetric body simulated on a 2D Cartesian grid requires adding a geometric curvature heat source term:

$$Q_{\text{metric}} = -\frac{k}{r} \frac{\partial T}{\partial r} = \frac{k}{r_{\text{eff}}^2} \left[ (x - x_c)\frac{\partial T}{\partial x} + (y - y_c)\frac{\partial T}{\partial y} \right]$$

with core regularization radius $r_{\text{eff}} = \sqrt{(x - x_c)^2 + (y - y_c)^2 + \epsilon_r^2}$.

### Analytical Benchmark
For steady-state conduction in a sphere of radius $R$ with uniform volumetric heating $H$ and surface temperature $T_s$:
- 3D Spherical analytical solution: $\Delta T_{\text{sph}} = T_{\text{center}} - T_s = \frac{H R^2}{6 k}$
- 2D Cylindrical / planar solution: $\Delta T_{\text{cyl}} = \frac{H R^2}{4 k}$
- Core temperature ratio:

$$\frac{\Delta T_{\text{sph}}}{\Delta T_{\text{cyl}}} = \frac{4}{6} = \frac{2}{3} \approx 0.6667$$

### Verification Test Suite
- `test/test_geometry_radiation.jl`: `@testset "Analytical 1D Steady-State Spherical Conduction Benchmark"`  
  Verifies that the numerical 2D Cartesian solve with $Q_{\text{metric}}$ recovers the analytical $2/3$ temperature gradient ratio within 1% accuracy.
- `test/test_geometry_radiation.jl`: `@testset "Spherical Geometric Metric Divergence"`  
  Verifies core regularization, directional sign consistency ($Q_{\text{metric}} < 0$ for outwardly decreasing temperature), and sticky-air vanishing.

---

## 3. Darcy Thermal Buoyancy and Fluid Equation of State

### Governing Formulation
The fluid thermal expansion equation of state and the Darcy thermal buoyancy driving force are verified against theoretical scaling:

$$\mathbf{q}_D = - \frac{k_\phi}{\eta_f} \left(\nabla P_f - \rho_f(T) \mathbf{g}\right)$$

In a hydrostatic vertical column with temperature contrast $\Delta T = T - T_{\text{melt}}$, the upward buoyant Darcy discharge velocity scaling is:

$$|q_{yD}| = \frac{k_\phi}{\eta_f} \rho_{f0} \alpha_f \Delta T g$$

### Parameterization Behavior

![Darcy Thermal Buoyancy Verification](../assets/darcy_buoyancy_verification.png)

*Figure 1: Thermal buoyancy and fluid equation of state verification in Erebus. (a) Temperature-dependent fluid density $\rho_f(T)$ over the range $T \in [240, 700]\text{ K}$ displaying sub-freezing ice density ($\rho_{\text{ice}} = 917\text{ kg/m}^3$), liquid water density at $T_{\text{melt}} = 273.0\text{ K}$ ($\rho_{\text{water}} = 1000\text{ kg/m}^3$), and linear thermal expansion above melting. Curves compare the code baseline ($\alpha_f = 5\times 10^{-5}\text{ K}^{-1}$) against ambient water and hydrothermal regimes. (b) Upward buoyant Darcy discharge velocity $|q_{yD}|$ as a function of thermal contrast $\Delta T$ for representative crustal permeabilities ($k_\phi \in [10^{-14}, 10^{-12}]\text{ m}^2$) at the code baseline $\alpha_f = 5\times 10^{-5}\text{ K}^{-1}$.*

### Verification Test Suite
- `test/test_physics.jl`: Thermal expansion and buoyancy driving forces


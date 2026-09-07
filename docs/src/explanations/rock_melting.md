# Silicate Rock Melting and Magma Ocean Dynamics

This page explains the physical formulations, constitutive equations, and numerical regularizations for silicate rock melting and sub-grid magma ocean convection in `Erebus.jl`.

---

## Physical Context: Planetesimal Differentiation

Planetesimals that accreted during the first 1.5 million years of Solar System history incorporated substantial abundances of short-lived radioactive isotopes, primarily $^{26}\text{Al}$ ($t_{1/2} \approx 0.717\text{ Ma}$) and $^{60}\text{Fe}$ ($t_{1/2} \approx 2.62\text{ Ma}$).
Volumetric decay heating in bodies larger than approximately 20 to 30 km in radius drove internal temperatures past the silicate solidus ($T_{\text{sol}} \approx 1400\text{ K}$).

As silicate minerals melt, the mechanical and thermal state of the planetesimal undergoes a dramatic transition:
1. **Partially Molten Silicate Mush:**
   Between the solidus and the rheologically critical melt fraction, ($\phi_{\text{crit}} \approx 0.40$), solid silicate grains remain in physical contact.
   Permeability increases, and buoyant silicate melt or immiscible metallic liquid can segregate through porous percolation.
2. **Rheological Disaggregation and Magma Ocean Transition:**
   When the melt fraction exceeds $\phi_{\text{crit}}$, solid grains lose physical contact and become suspended in liquid silicate melt.
   Viscosity drops by ten to fifteen orders of magnitude, from solid mantle rock ($\approx 10^{18}\text{ to }10^{21}\text{ Pa}\cdot\text{s}$) down to liquid silicate magma ($\approx 10^{-2}\text{ to }10^2\text{ Pa}\cdot\text{s}$).
3. **Vigorous Turbulent Convection:**
   In the liquid suspension regime, thermal convection operates at extreme Rayleigh numbers ($\text{Ra} \sim 10^{15}\text{ to }10^{20}$).
   Vigorous convective motion rapidly transports core heat outward to the surface, buffering interior temperatures and governing core formation.

Resolving both viscous porous compaction and turbulent magma ocean circulation on a single computational grid presents a fundamental challenge.
`Erebus.jl` addresses this challenge by combining continuum porous-media mechanics with a regularized sub-grid soft turbulence model.

---

## Thermodynamics of Silicate Melting

### 1. Pressure-Dependent Melting Interval

Planetesimal interior pressure $P$ alters the thermodynamic equilibrium between solid and liquid silicates.
Following linear Clapeyron slopes, the solidus $T_s(P)$ and liquidus $T_l(P)$ shift upward with lithostatic pressure:

$$T_s(P) = T_{s,0} + \frac{dT}{dP} P$$

$$T_l(P) = T_{l,0} + \frac{dT}{dP} P$$

| Parameter | Description | Default Value | Units |
|:---|:---|:---|:---|
| $T_{s,0}$ | Silicate solidus temperature at zero pressure | $1400.0$ | $\text{K}$ |
| $T_{l,0}$ | Silicate liquidus temperature at zero pressure | $1800.0$ | $\text{K}$ |
| $\frac{dT}{dP}$ | Clapeyron melting slope with pressure | $0.0$ ($1.2\times 10^{-7}$ when enabled) | $\text{K/Pa}$ ($120\text{ K/GPa}$) |

Shifting both the solidus and the liquidus by the identical slope $dT/dP$ preserves the melting temperature interval:

$$\Delta T_{\text{melt}} = T_l(P) - T_s(P) = T_{l,0} - T_{s,0} = 400.0\text{ K}$$

Preserving $\Delta T_{\text{melt}}$ at all depths prevents artificial distortion of the phase-change interval under lithostatic overburden.

### 2. Silicate Melt Fraction

The equilibrium silicate melt fraction $F_m(T, P)$ is evaluated by linear interpolation between the local solidus and liquidus:

$$F_m(T, P) = \text{clamp}\left( \frac{T - T_s(P)}{T_l(P) - T_s(P)}, 0.0, 1.0 \right)$$

For sticky air and non-rock boundary phases (marker phase type $tm \ge 3$), $F_m$ is fixed to zero.
Input validation asserts finite, strictly positive temperatures and valid pressure states before evaluation.

### 3. Latent Heat Buffering via Apparent Heat Capacity

Phase transitions absorb or release latent heat of fusion $L_m \approx 4.0\times 10^5\text{ J/kg}$.
`Erebus.jl` incorporates latent heat buffering through the apparent heat capacity method (Gerya, 2019, Section 16.6):

$$\rho c_{p,\text{eff}} = \rho c_p + \rho_s L_m \frac{\partial F_m}{\partial T}$$

Because $\Delta T_{\text{melt}} = T_l(P) - T_s(P)$ remains constant under the parallel Clapeyron shift, the derivative in the melting interval is constant:

$$\frac{\partial F_m}{\partial T} = \begin{cases} \frac{1}{T_l - T_s} & \text{if } T_s(P) \le T \le T_l(P) \\ 0.0 & \text{otherwise} \end{cases}$$

Integrating the effective heat capacity over the solidus-liquidus interval yields the exact total latent heat enthalpy:

$$\int_{T_s}^{T_l} \left(\rho c_{p,\text{eff}} - \rho c_p\right) \, dT = \rho_s L_m \int_{T_s}^{T_l} \frac{1}{T_l - T_s} \, dT = \rho_s L_m$$

This formulation satisfies total energy conservation without requiring non-linear enthalpy iterations.

---

## Magma Suspension Rheology

The effective viscosity of partially molten rock depends strongly on the silicate melt fraction $F_m$.
`Erebus.jl` implements a two-regime rheological model based on Costa et al. (2009) and Gerya (2019, Section 16.6.2).

### 1. Partially Molten Regime ($F_m < \phi_{\text{crit}}$)

When the melt fraction remains below the rheologically critical threshold $\phi_{\text{crit}} = 0.40$, solid grains form a continuous framework.
Melt accumulation weakens the solid framework through exponential softening:

$$\eta(F_m) = \eta_{\text{solid}} \exp\left(-\alpha_\eta F_m\right)$$

where $\eta_{\text{solid}}$ is the unmolten rock matrix viscosity, and $\alpha_\eta = 28.0$ is the melt weakening exponent.
At $F_m = 0.40$, this exponential reduction decreases viscosity by a factor of $\exp(-11.2) \approx 1.37\times 10^{-5}$.

### 2. Disaggregated Suspension Regime ($F_m \ge \phi_{\text{crit}}$)

When $F_m \ge \phi_{\text{crit}}$, solid grains detach from one another, transforming the aggregate into a liquid suspension.
Viscosity transitions from the disaggregation value $\eta(\phi_{\text{crit}})$ down to the viscosity of pure liquid silicate melt $\eta_{\text{melt}} \approx 10\text{ Pa}\cdot\text{s}$.
`Erebus.jl` employs continuous log-linear interpolation over the suspension interval:

$$\log_{10} \eta(F_m) = \log_{10} \eta(\phi_{\text{crit}}) - \frac{F_m - \phi_{\text{crit}}}{1.0 - \phi_{\text{crit}}} \left[ \log_{10} \eta(\phi_{\text{crit}}) - \log_{10} \eta_{\text{melt}} \right]$$

This formulation guarantees $C^0$ continuity at the disaggregation boundary $F_m = \phi_{\text{crit}}$.

### 3. Solver Viscosity Clamping

To maintain numerical stability in the implicit Stokes-Darcy velocity-pressure solver, the mechanical viscosity is bounded by admissible limits:

$$\eta_{\text{eff}} = \max\left(\eta_{\min}, \min\left(\eta_{\max}, \eta(F_m)\right)\right)$$

The default numerical lower limit is $\eta_{\min} = 10^{12}\text{ Pa}\cdot\text{s}$.
This floor prevents ill-conditioning of the discrete Stokes matrix operator while allowing realistic flow patterns in solid and partially molten regions.

---

## Sub-Grid Soft Turbulence Model

### 1. Physical Motivation

Liquid silicate magma possesses a very low dynamic viscosity ($\eta_{\text{fluid}} \sim 10^{-1}\text{ to }10^2\text{ Pa}\cdot\text{s}$) and density $\rho \sim 2800\text{ kg/m}^3$.
In a molten planetesimal interior ($R \approx 50\text{ km}$), thermal convection operates at Rayleigh numbers exceeding $10^{16}$.
The associated thermal boundary layer thickness scales as $\delta \sim R \, \text{Ra}^{-1/3} \sim 0.2\text{ m}$.

Direct numerical simulation of such thin boundary layers on a planetary grid with kilometer-scale cell resolution ($\Delta x \sim 1\text{ km}$) is impossible.
Clamping mechanical viscosity to $\eta_{\min} = 10^{12}\text{ Pa}\cdot\text{s}$ suppresses resolved convective velocities by many orders of magnitude.
Without a sub-grid model, resolved advection cannot carry the heat flux, causing unphysical super-liquidus overheating.

### 2. Convective Heat Flux Parameterization

Following Solomatov (2007) and boundary-layer scaling theory for turbulent convection (Kraichnan 1962; Solomatov & Stevenson 1993), the convective heat flux $F_{\text{conv}}$ scales inversely with fluid viscosity:

$$F_{\text{conv}} \propto \left(\frac{1}{\eta}\right)^\beta$$

Because the numerical Stokes solver enforces an artificial viscosity floor, $\eta_{\text{num}} \ge \eta_{\min} \gg \eta_{\text{fluid}}$, the missing convective heat transport is represented as an effective turbulent thermal conductivity $k_{\text{turb}}$:

$$k_{\text{turb}} = k_{\text{cond}} \left(\frac{\eta_{\text{num}}}{\eta_{\text{fluid}}}\right)^\beta$$

| Parameter | Symbol | Default Value | Physical Meaning |
|:---|:---|:---|:---|
| Background conductivity | $k_{\text{cond}}$ | $3.0\text{ W}/(\text{m}\cdot\text{K})$ | Solid silicate lattice conductivity |
| Numerical viscosity | $\eta_{\text{num}}$ | $\ge 10^{12}\text{ Pa}\cdot\text{s}$ | Viscosity evaluated on the numerical grid |
| Magma fluid viscosity | $\eta_{\text{fluid}}$ | $100.0\text{ Pa}\cdot\text{s}$ | True physical liquid silicate viscosity |
| Turbulent scaling exponent | $\beta$ | $1/3 \approx 0.3333$ | Solomatov (2007) soft turbulence exponent |

The scaling exponent $\beta = 1/3$ derives from asymptotic boundary layer theory in the soft turbulence regime, where convective heat transport is governed by boundary layer instability.
An exponent of $\beta = 1/2$ corresponds to classical laminar boundary layer parameterizations, such as that employed in `i2elvis` (Gerya 2019).
`Erebus.jl` supports arbitrary user-configured exponents via the `turb_exponent` setting in TOML configurations.

---

## Regularization and Singularity Elimination

Discontinuous step thresholds at marker state transitions produce numerical artifacts. When conductivity jumps discontinuously by three orders of magnitude at $F_m = 0.40$, the spatial flux derivative develops a Dirac $\delta$-function spike.
This step change introduces three severe numerical pathologies:
1. **Dirac Flux Spikes:**
   A step change in conductivity between adjacent grid cells generates a discontinuous spatial derivative $\nabla \cdot (k \nabla T)$, producing non-physical temperature oscillations.
2. **Conductivity Suppression at Melting Onset:**
   If $\eta_{\text{num}}$ weakens before $\eta_{\text{fluid}}$ is adjusted, the ratio $\eta_{\text{num}} / \eta_{\text{fluid}}$ can fall below unity, causing $k_{\text{turb}} < k_{\text{cond}}$.
   This artificially suppresses heat conduction at the exact onset of melting.
3. **Cold Boundary Singularities:**
   Markers near the cold planetesimal surface that retain high melt fractions can calculate large convective heat conductivities, draining surface energy unnaturally.

`Erebus.jl` eliminates these pathologies with a threefold regularization scheme.

### 1. Cubic Smoothstep Blending

Conductivity is blended continuously over a finite melt transition interval $[F_{\text{start}}, F_{\text{end}}] = [0.30, 0.50]$ using a Hermite cubic smoothstep function:

$$\xi = \text{clamp}\left( \frac{F_m - F_{\text{start}}}{F_{\text{end}} - F_{\text{start}}}, 0.0, 1.0 \right)$$

$$w_F(\xi) = 3 \xi^2 - 2 \xi^3$$

The function $w_F(\xi)$ satisfies $w_F(0) = 0$, $w_F(1) = 1$, and possesses vanishing first derivatives at both boundaries:

$$\left.\frac{dw_F}{d\xi}\right|_{\xi=0} = 0, \qquad \left.\frac{dw_F}{d\xi}\right|_{\xi=1} = 0$$

This property yields $C^1$ smoothness of the effective thermal conductivity at both solid and liquid transition thresholds.

### 2. Logarithmic Geometric Interpolation

Because thermal conductivity varies over several orders of magnitude between conduction ($3\text{ W}/(\text{m}\cdot\text{K})$) and turbulent convection ($10^3\text{ to }10^5\text{ W}/(\text{m}\cdot\text{K})$), linear interpolation produces steep gradients.
`Erebus.jl` blends conductivity logarithmically:

$$\log_{10} k_{\text{eff}} = (1 - w_{\text{total}}) \log_{10} k_{\text{cond}} + w_{\text{total}} \log_{10} k_{\text{turb}}$$

$$k_{\text{eff}} = 10^{\log_{10} k_{\text{eff}}}$$

Strict monotonicity is enforced by clamping $k_{\text{eff}} = \max(k_{\text{cond}}, k_{\text{eff}})$; this condition guarantees that convective enhancement never reduces conductivity below background solid conduction.

### 3. Viscosity Matching in the Mush Interval

To prevent the unphysical drop where $\eta_{\text{num}} < \eta_{\text{fluid}}$ during the initial stages of melting, the reference fluid viscosity $\eta_{\text{fluid}}$ is blended from the current matrix viscosity down to the magma viscosity:

$$\log_{10} \eta_{\text{fluid}}(F_m) = (1 - \xi) \log_{10} \eta_{\text{matrix}} + \xi \log_{10} \eta_{\text{melt}}$$

This formulation maintains $\eta_{\text{num}} / \eta_{\text{fluid}} \ge 1.0$ throughout the entire transition window.

### 4. Quadratic Surface Thermal Contrast Weighting

To prevent convective heat enhancement from penetrating cold planetesimal crust or radiating surfaces, a thermal contrast weight $w_T$ modulates the transition:

$$\Delta T = \max(0.0, T - T_{\text{surface}})$$

$$w_T = \left[ \text{clamp}\left( \frac{\Delta T}{\Delta T_{\min}}, 0.0, 1.0 \right) \right]^2$$

The quadratic exponent guarantees that $\partial w_T / \partial T \to 0$ continuously as $T \to T_{\text{surface}}$, which eliminates gradient kinks at the surface boundary.
The complete transition weight is:

$$w_{\text{total}} = w_F(\xi) \cdot w_T$$

When $\Delta T \to 0$ or $F_m \le F_{\text{start}}$, $w_{\text{total}} \to 0$, recovering exact solid lattice conductivity $k_{\text{cond}}$.

---

## Numerical Integration Architecture

The melting and soft turbulence system operates on both Lagrangian markers and Eulerian grid cells:
1. **Marker Update:**
   Marker temperatures and lithostatic pressures determine $F_m$ for each marker.
   Marker viscosities are computed using the suspension rheology model.
2. **Marker-to-Grid Interpolation:**
   Harmonic and arithmetic interpolation schemes transfer marker viscosities and densities to staggered velocity and pressure nodes.
3. **Cell Conductivity Evaluation:**
   Effective thermal conductivities $k_{\text{eff}}$ are evaluated at staggered cell faces using the regularized soft turbulence formulation.
4. **Implicit Energy Solver:**
   The parabolic thermal energy equation is assembled into a linear system with variable apparent heat capacity and variable thermal conductivity, and is solved implicitly to maintain numerical stability.

---

## Configuration Controls

Melting and turbulence parameters are set in the `[melting]` table of simulation configuration files. Key controls include:

- `active`: Enables silicate melting thermodynamics and apparent heat capacity buffering.
- `soft_turbulence`: Enables regularized sub-grid soft turbulence convection.
- `turb_exponent`: Convective heat flux scaling exponent ($\beta = 1/3$ for Solomatov 2007 soft turbulence; $\beta = 1/2$ for laminar boundary layers).
- `eta_fluid_silicate`: Dynamic viscosity of liquid silicate magma ($100.0\text{ Pa}\cdot\text{s}$).
- `F_turb_start` and `F_turb_end`: Melt fraction bounds for smoothstep conductivity blending ($[0.30, 0.50]$).
- `T_solidus` and `T_liquidus`: Phase-dependent solidus and liquidus temperatures.
- `dpdt_clapeyron`: Pressure-dependent melting slope $dT/dP$ (default $0.0\text{ K/Pa}$; typical $1.2\times 10^{-7}\text{ K/Pa}$).

For complete schema details, default values, and data types, see the [Configuration Schema Reference](../reference/config_schema.md#melting).

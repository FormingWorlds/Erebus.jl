# Silicate Rock Melting and Magma Rheology

This page documents and validates the silicate rock melting formulation, latent heat buffering, magma rheology, and sub-grid soft turbulence model in `Erebus.jl`.

## Theoretical Formulation

The thermodynamics of silicate melting (linear melt fraction $F_m$, apparent heat capacity $\rho c_{p,\text{eff}}$, latent heat buffering), Costa et al. (2009) rheological weakening, and Solomatov (2007) sub-grid soft turbulence scaling are derived in detail in [Silicate Melting & Soft Turbulence](../explanations/rock_melting.md).

Key constitutive formulations validated here include:

- **Linear Melt Fraction ($F_m$):** Interpolated across the Clapeyron-shifted solidus $T_s(P)$ and liquidus $T_l(P)$:
  $$F_m(T, P) = \text{clamp}\left(\frac{T - T_s(P)}{T_l(P) - T_s(P)}, 0.0, 1.0\right)$$
- **Apparent Heat Capacity ($c_{p,\text{eff}}$):** Conserving latent heat $L_m = 4.0\times 10^5\text{ J/kg}$ over the melting interval:
  $$\rho c_{p,\text{eff}} = \rho c_p + \rho_s L_m \frac{\partial F_m}{\partial T}$$
- **Melt-Weakened Rheology ($\eta$):** Exponential softening below $\phi_{\text{crit}} = 0.40$ and suspension breakdown above $\phi_{\text{crit}}$ (Costa et al., 2009; Gerya, 2019, Section 16.6.2):
  $$\eta(F_m) = \begin{cases} \eta_{\text{solid}} \exp(-\alpha_\eta F_m), & F_m \le \phi_{\text{crit}} \\ \eta(\phi_{\text{crit}}) \left(\frac{\eta_{\text{melt}}}{\eta(\phi_{\text{crit}})}\right)^{\frac{F_m - \phi_{\text{crit}}}{1 - \phi_{\text{crit}}}}, & F_m > \phi_{\text{crit}} \end{cases}$$
- **Sub-Grid Soft Turbulence ($k_{\text{turb}}$):** Scaling effective conductivity by $\text{Nu} \sim \text{Ra}^{1/3}$ in liquid magma (Solomatov, 2007):
  $$k_{\text{turb}} = k_{\text{cond}} \left(\frac{\eta_{\text{num}}}{\eta_{\text{fluid}}}\right)^\beta, \quad \beta = 1/3$$

---

## Planetesimal Magma Ocean Solidification Benchmark

The benchmark simulates the cooling and solidification of a 50 km radius planetesimal using a 1D implicit spherical finite-volume solver, coupled with the regularized soft turbulence parameterization and apparent heat capacity formulation. The radial field is revolved into a circular cross-section for 2D visualization and video animation.

The interior possesses an initial 30 km radius molten core at 1850 K ($F_m = 1.0$, super-liquidus), enclosed by a 20 km conductive solid crust that tapers to a surface temperature of 300 K. Sticky air outside the planetesimal displays in pure white.

The benchmark tracks thermal relaxation, sub-grid convective heat transport, phase change enthalpy, and solidification front retreat over 50 kyr.

The benchmark verifies four physical mechanisms:
1. **Convective Heat Extraction**: Solomatov soft turbulence enhances thermal conductivity by more than three orders of magnitude in molten silicate, cooling the interior efficiently.
2. **Phase Boundary Transition**: The cubic smoothstep transitions conductivity smoothly without numerical spikes through the mush shell ($F_m \in [0.30, 0.50]$).
3. **Enthalpy Buffering**: Latent heat release slows solidification at the solidus-liquidus interface ($T \in [1400, 1800]\text{ K}$).
4. **Front Propagation & Energy Balance**: Heat discharged through the surface balances internal energy loss over the cooling trajectory.

### Benchmark Results

The benchmark tracks thermal state, phase fraction, and convective heat transport:

![Planetesimal Magma Ocean Benchmark Summary](../assets/magma_ocean_cooling_benchmark.png)

- **(a) Revolved Thermal Field Snapshot ($t = 15\text{ kyr}$):** Planetesimal center sits at origin $(0, 0)\text{ km}$. The interior core cools to 1720 K while maintaining a sharp boundary layer beneath the solid crust. The solidus contour (1400 K) and liquidus contour (1800 K) mark the crystallization zone.
- **(b) Core Thermal Quenching:** In the baseline conduction model (Soft Turb. OFF), central core temperature stays at 1850 K for 50 kyr because conductive diffusion through 50 km requires $\sim 25\text{ Myr}$. With regularized soft turbulence active (Soft Turb. ON), core temperature drops below liquidus (1800 K) in 2 kyr and reaches 1660 K at 50 kyr.
- **(c) Magma Ocean Solidification Front:** Tracks the radial retreat of the rheological breakdown front ($F_m = 0.40$). Turbulent mixing delivers heat to the front, controlling the freezing velocity.
- **(d) Radial Temperature Profiles:** Profiles at $t \in [0, 5, 15, 50]\text{ kyr}$ display convective flattening in the core ($r < 35\text{ km}$) and steep conductive gradients in the outer crust ($r \in [35, 50]\text{ km}$).
- **(e) Convective Conductivity Profiles:** Effective thermal conductivity reaches $k_{\text{eff}} \approx 7 \times 10^3\text{ W/(m K)}$ in the liquid core, declining smoothly to $k_{\text{cond}} = 3.0\text{ W/(m K)}$ within the mush layer without artificial jumps.
- **(f) Planetary Heat Loss ($q_{\text{surf}}$):** Surface heat flux starts at $0.23\text{ W/m}^2$, sustaining heat discharge through the conductive lid.

---

### Planetesimal Magma Ocean Solidification Video (Revolved Spherical Model)

The animation below displays the benchmark ($N_r = 128$ radial cells, revolved onto a 128x128 visualization mesh, 50 km planetesimal) over 50 kyr. The panels show temperature $T$ (left), silicate melt fraction $F_m$ (center), and effective thermal conductivity $k_{\text{eff}}$ on a logarithmic scale (right). Color limits stay fixed and normalized in all frames.

![Planetesimal Magma Ocean Solidification Animation](../assets/magma_ocean_cooling_128.gif)

---

### Grid Convergence (32, 64, 128, and 256 cells)

To test spatial convergence, simulations compare four radial grid resolutions: $N_r = 32$ ($\Delta r = 1.56\text{ km}$), $N_r = 64$ ($\Delta r = 0.78\text{ km}$), $N_r = 128$ ($\Delta r = 0.39\text{ km}$), and $N_r = 256$ ($\Delta r = 0.20\text{ km}$).

![Grid Convergence Comparison](../assets/magma_ocean_grid_convergence.png)

Metrics demonstrate spatial convergence:
- **Thermal Match:** Core temperature at 15 kyr reaches 1765.2 K at $N_r = 32$, 1740.1 K at $N_r = 64$, 1721.4 K at $N_r = 128$, and 1704.8 K at $N_r = 256$. The relative difference between $N_r = 128$ and $N_r = 256$ is 0.97%.
- **Solidification Front Match:** The crystallization front radius at 15 kyr converges to $r_{\text{melt}} = 37.1\text{ km}$, differing by less than one grid cell width between all resolutions.
- **Monotonic Progression:** Central temperature decreases monotonically toward the continuum limit as grid resolution is refined.

---

### Conductivity Regularization and Singularity Elimination
 
Discontinuous step thresholds at marker state transitions produce numerical artifacts. When conductivity jumps discontinuously by three orders of magnitude at $F_m = 0.40$, the spatial flux derivative develops a Dirac $\delta$-function spike.
 
`Erebus.jl` resolves this issue with regularized geometric blending:
 
![Conductivity Regularization](../assets/magma_ocean_regularization.png)
 
The regularization provides three improvements:
1. **$C^1$ Smoothness in Mush Interval:** The cubic smoothstep provides continuous first derivatives $d(\log_{10} k_{\text{eff}})/dF_m$ throughout the melting interval $[F_{\text{start}}, F_{\text{end}}]$. This eliminates singular flux spikes.
2. **Viscosity Matching:** Blending $\eta_{\text{fluid}}$ from matrix viscosity down to liquid silicate viscosity prevents the conductivity dip at melting onset.
3. **Surface Boundary Weighting:** The quadratic thermal contrast weight $w_T = [\text{clamp}(\Delta T / \Delta T_{\text{min}}, 0, 1)]^2$ forces $k_{\text{eff}} \to k_{\text{cond}}$ smoothly as $\Delta T \to 0$, preventing isothermal boundary artifacts.
 
---
 
## Analytical Verification
 
The implementation is verified against the following benchmarks:
1. **Enthalpy Conservation**: Numerical integration of apparent heat capacity over the solidus-liquidus interval recovers the theoretical latent heat energy within $10^{-6}$ relative error in unit tests (midpoint Riemann sum).
2. **Rheological Invariants**: Viscosity remains monotonic with $F_m$ for fixed solid viscosity, strictly positive, and continuous at the transition $F_m = \phi_{\text{crit}}$.
3. **Soft Turbulence Monotonicity & Asymptotics**: Effective thermal conductivity matches $k_{\text{cond}}$ exactly for $F_m \le F_{\text{start}}$ or isothermal conditions, recovers $k_{\text{turb}}$ for $F_m \ge F_{\text{end}}$, and satisfies $k_{\text{eff}} \ge k_{\text{cond}}$ monotonically in $F_m$ over the transition window for constant viscosities.
4. **Planetesimal Magma Ocean Solidification Benchmark**: Effective convective conductivity transports core heat to the surface, buffering interior temperatures and preventing runaway super-liquidus overheating.

---

## Configurations

Model setup files live in `configs/`:
- `magma_ocean_cooling_turb_on_128.toml` (High-resolution benchmark, soft turbulence enabled, 128x128)
- `magma_ocean_cooling_turb_off_128.toml` (Baseline benchmark, conduction only, 128x128)
- `magma_ocean_cooling_turb_on_64.toml` (Medium-resolution benchmark, 64x64)
- `magma_ocean_cooling_turb_on_32.toml` (Fast benchmark, 32x32)

---

## Verification Test Suite

- `test/test_melting.jl`:
  - `@testset "Silicate Rock Melting & Magma Rheology"`
  - `@testset "compute_melt_fraction: Analytical Limits & Monotonicity"`
  - `@testset "rhocp_apparent_silicate: Latent Heat & Conservation Integral"`
  - `@testset "compute_melt_weakened_viscosity: Rheological Transition"`
- `test/test_soft_turbulence.jl`:
  - `@testset "Sub-Grid Soft Turbulence & Regularized Conductivity"`
  - `@testset "Physics: regularized_soft_turbulence_conductivity"`
  - `@testset "Grid Interpolation: KX & KY receive enhanced conductivity"`
  - `@testset "Mini-Simulation Execution with Soft Turbulence"`
- `test/test_magma_transport.jl`:
  - `@testset "Magma Transport Configuration & Validation"`
  - `@testset "Silicate Melt Permeability Formulation"`
  - `@testset "Silicate Melt Segregation Velocity & Regimes"`
  - `@testset "Silicate Melt Gravitational Dissipation Heating"`
  - `@testset "Marker Magma Allocation & Depletion Properties"`
  - `@testset "Two-Phase Melt Segregation Operator: Mass Conservation & Ascent"`
  - `@testset "Two-Phase Melt Segregation: Thermal Dissipation & Crystallization Latent Heat"`
  - `@testset "Two-Phase Melt Segregation: Neutral & Negative Buoyancy Cutoff"`

---

## Two-Phase Silicate Melt Segregation Verification

The two-phase silicate melt segregation, porous percolation, hindered settling, and conservative drift-flux transport implementation is verified by dedicated unit and integration tests:

1. **Analytical Regime Limits:**
   - Darcy percolation velocity matches $v_{\text{perc}} = (k_\phi / \eta_{\text{melt}} F_m) \Delta\rho g$ exactly when $F_m \le F_{\text{perc\_end}}$.
   - Hindered settling velocity matches $v_{\text{susp}} = v_{\text{Stokes}} F_m^n$ exactly when $F_m \ge F_{\text{settle\_start}}$, and recovers unhindered Stokes settling when $F_m \to 1.0$.
   - The Hermite cubic blend smoothly interpolates between Darcy percolation and Stokes crystal settling without discontinuities or negative derivatives.
   - Permeability and segregation velocities evaluate to exactly zero when melt fraction is at or below the residual threshold $\phi_{\text{residual}}$.
   - Segregation velocity vanishes under neutral buoyancy ($\Delta\rho = 0$) or negative buoyancy ($\Delta\rho < 0$).

2. **Machine-Precision Mass Conservation:**
   - Conservative finite-volume donor-receiver flux limiting conserves total silicate melt mass to machine precision across all subcycles ($< 10^{-12}$ relative drift).
   - Clamping prevents donor cell melt fraction from becoming negative ($F_m \ge 0$) and receiver cells from exceeding the packing ceiling ($F_m \le \phi_{\text{pack}}$).

3. **Thermal and Depletion Coupling:**
   - Gravitational shear dissipation heating rates satisfy $\Psi = \Delta\rho g F_m v_{\text{seg}}$ and accumulate non-negative thermal energy into the energy solver.
   - Subsolidus crystallization latent heat releases $Q_{\text{lat}} = \dot{M}_{\text{cryst}} L_m$, buffering temperature drops when buoyant melt rises into subsolidus crustal regions.
   - Mantle depletion tracking accurately records cumulative extracted melt on Lagrangian markers and prevents remelting of depleted residues.


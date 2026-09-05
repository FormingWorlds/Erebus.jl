# Verification and Validation

This page summarizes the benchmarks and tests used to ensure numerical accuracy and physical fidelity in `Erebus.jl`.

---

## 1. 1D Terzaghi Analytical Consolidation Benchmark

For the coupled Stokes-Darcy formulation, the primary verification test is the one-dimensional Terzaghi consolidation problem.

- **Setup**: Saturated porous column of height $H$, with draining boundary anchors ($P_f = p_{\text{surface}} = 1000\text{ Pa}$), subjected to excess pore pressure dissipation.
- **Analytical Solution**: Fourier series representation of two-way pore pressure dissipation over dimensionless time factor $T_v = c_v t / H^2$.
- **Result**: The numerical pore pressure field, computed by the monolithic solver, matches the analytical Fourier series solution at all column depths, with a relative error $< 3.5\%$. See the [Terzaghi Consolidation Tutorial](../tutorials/terzaghi_consolidation.md) for the full derivation and code.

![1D Terzaghi Consolidation Benchmark Verification](../assets/terzaghi_benchmark.png)

*Figure 1: Benchmark verification of the coupled Stokes-Darcy solver against the analytical 1D Terzaghi consolidation solution. (a) Pore fluid pressure dissipation along column height $y \in [0, H]$ at three successive timesteps, comparing the analytical Fourier series (dashed curves) with Erebus numerical results (dots). (b) Pointwise relative error $|p_f^{\text{num}} - p_f^{\text{ana}}| / p_0$ showing that numerical discretization error remains strictly below $3.5\%$ throughout the column.*

---

## 2. Poroelastic Constitutive Limits

In `src/physics.jl`, the constitutive poroelastic functions are verified against physical asymptotic limits:

1. **Incompressible Solid Skeleton ($\beta_s \to 0$)**:
   $$\lim_{\beta_s \to 0} K_{\text{BW}} = 1, \quad \lim_{\beta_s \to 0} B = \frac{\beta_\phi}{\beta_\phi + \phi(1 - \phi)\beta_f}$$
   Verified over porosity values $\phi \in [\phi_{\text{min}}, \phi_{\text{max}}]$.

2. **Incompressible Pore Fluid ($\beta_f \to 0$)**:
   As fluid compressibility approaches zero, Skempton coefficient, $B$, approaches its undrained upper bound:
   $$\lim_{\beta_f \to 0} B = \min\left(1, \frac{\beta_d - \beta_s}{\beta_d - (1 + \phi)\beta_s}\right) = 1$$
   where the code clamps the theoretical ratio to $[0, 1]$, to enforce the physical upper bound.

3. **Porosity Bounding Guarantees**:
   Constitutive routines clamp porosity to $[\phi_{\text{min}}, \phi_{\text{max}}]$, to prevent singular division when $\phi \to 0$, or $\phi \to 1$.

![Poroelastic Constitutive Limits Verification](../assets/poroelastic_verification.png)

*Figure 2: Theoretical behavior of derived poroelastic coefficients in Erebus. (a) Biot-Willis coefficient $K_{\text{BW}}$ as a function of porosity $\phi$ for varied solid grain compressibility $\beta_s$ to confirm asymptotic convergence toward unity ($K_{\text{BW}} \equiv 1$) in the incompressible solid grain limit. (b) Skempton pore pressure coefficient $B$ as a function of fluid compressibility $\beta_f$ for representative porosity values to display undrained response transitions.*

---

## 3. Matrix Block Consistency Tests

In `assemble_hydromechanical_lse!`, the individual coupling sub-blocks are tested independently against theoretical values:
- **Solid continuity diagonal ($L[P_t, P_t]$)**: Matches $+K_{\text{cont}} \left(\frac{1}{(1-\phi)\eta_\phi} + \frac{\beta_d}{\Delta t}\right)$.
- **Poroelastic cross-coupling ($L[P_t, P_f] = L[P_f, P_t]$)**: Matches $-K_{\text{cont}} \left(\frac{1}{(1-\phi)\eta_\phi} + \frac{\beta_d K_{\text{BW}}}{\Delta t}\right)$.
- **Fluid continuity diagonal ($L[P_f, P_f]$)**: Matches $+K_{\text{cont}} \left(\frac{1}{(1-\phi)\eta_\phi} + \frac{\beta_d K_{\text{BW}}}{B \Delta t}\right)$.

---

## 4. Radiogenic Energy and Mass Conservation

- **Radioactive Decay Kinetics**: Exponential decay of $^{26}\text{Al}$ and $^{60}\text{Fe}$ is verified against analytical integration, $\int_0^t Q(t') dt'$, to confirm exact conservation of energy release.
- **Dehydration Mass Balance**: Fluid released from decomposing hydrated silicates matches matrix mass loss ($DMP$).

---

## 5. Darcy Thermal Buoyancy Verification

The fluid thermal expansion equation of state, and the Darcy thermal buoyancy driving force, are verified against theoretical scaling:

$$\mathbf{q}_D = - \frac{k_\phi}{\eta_f} \left(\nabla P_f - \rho_f(T) \mathbf{g}\right)$$

In a hydrostatic vertical column, with temperature contrast $\Delta T = T - T_{\text{melt}}$, the upward buoyant Darcy discharge velocity scaling is:

$$|q_{yD}| = \frac{k_\phi}{\eta_f} \rho_{f0} \alpha_f \Delta T g$$

![Darcy Thermal Buoyancy Verification](../assets/darcy_buoyancy_verification.png)

*Figure 3: Thermal buoyancy and fluid equation of state verification in Erebus. (a) Temperature-dependent fluid density $\rho_f(T)$ over the range $T \in [240, 700]\text{ K}$ displaying sub-freezing ice density ($\rho_{\text{ice}} = 917\text{ kg/m}^3$), liquid water density at $T_{\text{melt}} = 273.0\text{ K}$ ($\rho_{\text{water}} = 1000\text{ kg/m}^3$), and linear thermal expansion above melting. Curves compare the code baseline ($\alpha_f = 5\times 10^{-5}\text{ K}^{-1}$) against ambient water and hydrothermal regimes. (b) Upward buoyant Darcy discharge velocity $|q_{yD}|$ as a function of thermal contrast $\Delta T$ for representative crustal permeabilities ($k_\phi \in [10^{-14}, 10^{-12}]\text{ m}^2$) at the code baseline $\alpha_f = 5\times 10^{-5}\text{ K}^{-1}$.*

---

## 6. Temperature-Dependent Fluid Viscosity Benchmarking

The temperature-dependent fluid viscosity formulation approximates liquid water behavior under planetesimal hydrothermal conditions:

$$\eta_f(T) = \eta_{f0} \exp\left[\frac{E_a}{R} \left(\frac{1}{T} - \frac{1}{T_0}\right)\right]$$

with reference viscosity $\eta_{f0} = 1.0\times 10^{-3}\text{ Pa}\cdot\text{s}$ at $T_0 = 293.15\text{ K}$ ($20^\circ\text{C}$) and activation energy $E_a = 15.0\text{ kJ/mol}$.

This single-activation energy relation closely tracks experimental liquid water data (IAPWS / NIST standard tables) within about $13\%$ in the sub-boiling regime $T \in [273, 373]\text{ K}$. At higher temperatures ($T \approx 470\text{ to }570\text{ K}$), liquid water curvature follows a Vogel-Fulcher-Tammann profile, where the single Arrhenius fit underestimates viscosity by $\approx 29\text{ to }43\%$ before entering the supercritical regime ($T > 647\text{ K}$).

![Temperature-Dependent Fluid Viscosity Benchmarking](../assets/fluid_viscosity_temperature.png)

*Figure 4: Benchmarking of temperature-dependent pore fluid viscosity $\eta_f(T)$ in Erebus. (a) Dynamic fluid viscosity over the range $T \in [270, 650]\text{ K}$ on a logarithmic scale, comparing the default Arrhenius model ($E_a = 15.0\text{ kJ/mol}$, blue curve) against experimental liquid water measurements from IAPWS/NIST standards (red circles). (b) Hydrothermal Darcy mobility enhancement factor $\eta_{f0} / \eta_f(T)$ illustrating the $5\times\text{ to }24\times$ increase in fluid percolation speed as interior temperatures rise in liquid hydrothermal conditions ($273\text{ to }600\text{ K}$).*

---

## 7. Dynamic Hydrofracture Permeability Verification

When pore fluid pressure exceeds total confining pressure plus the rock tensile strength, the dynamic hydrofracture formulation enhances Darcy permeability:

$$P_{\text{eff}} = P_t - P_f \le -\sigma_t$$

The effective permeability scaling:

$$k_\phi^{\text{eff}} = \min\left(k_\phi \cdot \left[1 + \kappa_{\text{frac}} \left(\frac{\max(0, -P_{\text{eff}} - \sigma_t)}{\sigma_t}\right)^\gamma\right], k_{\text{max}}\right)$$

is verified in intact and fractured regimes:

![Dynamic Hydrofracturing Verification](../assets/hydrofracture_verification.png)

*Figure 5: Verification of dynamic hydrofracturing permeability enhancement in Erebus. (a) Effective permeability $k_\phi^{\text{eff}}$ as a function of Terzaghi effective stress $P_{\text{eff}} = P_t - P_f$ for compressive ($P_{\text{eff}} > 0$), intact tensile ($-\sigma_t < P_{\text{eff}} \le 0$), and hydrofractured ($P_{\text{eff}} \le -\sigma_t$) regimes for representative matrix permeabilities ($k_0 \in [10^{-16}, 10^{-14}]\text{ m}^2$) at tensile strength $\sigma_t = 10\text{ MPa}$. (b) Permeability enhancement factor $k_{\text{eff}} / k_0$ as a function of normalized overpressure for scaling exponents $\gamma \in \{0.5, 1.0, 2.0\}$ at $\kappa_{\text{frac}} = 10^3$.*

---

## 8. 1D Stefan Moving-Boundary Dehydration Benchmark

In the equilibrium limit, the code verifies moving dehydration fronts, comparing marker reaction progress against the analytical Stefan similarity solution. During reaction progress, silicate dehydration absorbs latent heat, $L_{\text{eff}}$, at reaction temperature $T_{\text{eq}}$:

$$D^s + \text{H}_2\text{O} \rightleftharpoons W^s$$

where $W^s$ is hydrated silicate (serpentine), and $D^s$ is anhydrous silicate. In the equilibrium limit, when $\Delta t \ge \Delta t_{\text{reaction}}$, minerals dehydrate at phase-equilibrium temperature $T_{\text{eq}} = (\Delta H_{WD} + P_f \Delta V_{WD}) / \Delta S_{WD}$.

The 1D planar Stefan problem sets wall temperature, $T_{\text{wall}} > T_{\text{eq}}$, at boundary $y = 0$. The interface moves at position $s(t) = 2 \lambda \sqrt{\kappa t}$, where thermal diffusivity is $\kappa = k / (\rho c_p)$, and $\lambda$ solves the transcendental eigenvalue equation:

$$\lambda \exp(\lambda^2) \text{erf}(\lambda) = \frac{\text{Ste}}{\sqrt{\pi}}$$

with Stefan number, $\text{Ste} = c_p (T_{\text{wall}} - T_{\text{eq}}) / L_{\text{eff}}$, and effective latent heat, $L_{\text{eff}} = \Delta X_W^s \Delta H_{WD} / (M_D + M_{H_2O})$.

Energy balance at the moving interface couples thermal conduction to latent heat sink terms:

$$-k \left. \frac{\partial T}{\partial y} \right|_{s^-} = \rho L_{\text{eff}} \frac{ds}{dt}$$

In the numerical benchmark, `perform_thermochemical_reaction!` operates on markers across a 1D column. Markers behind the moving front ($y_m < s_{\text{ana}}$), where $T > T_{\text{eq}}$, dehydrate, releasing fluid, and absorbing latent heat ($\text{DHP} < 0$). Markers ahead of the front ($y_m > s_{\text{ana}}$), where $T \le T_{\text{eq}}$, remain hydrated. The dehydration front position, extracted from marker compositions, matches the analytical Stefan front position, $s(t)$, within marker spacing.

---

## 9. Two-Way Hydration-Dehydration Reaction Pathways and Poroelastic Coupling

The thermochemical reaction code supports reversible, two-way phase changes between dry silicate, pore fluid, and hydrous minerals:

1. Forward exothermic hydration (serpentinization): consumes pore water, reduces rock pore space ($\Delta\phi < 0$), and releases latent heat ($\text{DHP} > 0$).
2. Reverse endothermic dehydration: expels lattice water into pore space, creates pore space ($\Delta\phi > 0$), and absorbs heat ($\text{DHP} < 0$).

Reaction direction follows the sign of chemical affinity:

$$\Delta G_{WD} = \Delta H_{WD} - T \Delta S_{WD} + P_f \Delta V_{WD} + R_G T \ln\left(\frac{X_D^s}{X_W^s}\right)$$

Two kinetic regimes set the characteristic reaction timescale $\Delta t_{\text{reaction}}$:

- Hydration regime, when $X_{W,\text{eq}}^s > X_W^s$: Gaussian temperature dependence, with peak rate at $T = 543\text{ K}$ ($270^\circ\text{C}$).
- Dehydration regime, when $X_{W,\text{eq}}^s < X_W^s$: Arrhenius kinetics, active at high temperatures ($T > 600\text{ K}$).

Poroelastic volume change enters the two-phase mass balance through the compaction or dilation rate:

$$\text{DMP} = \frac{1 - R_V}{\Delta t}$$

where $R_V$ is the solid-plus-fluid volume ratio. The local fluid mass exchange diagnostic:

$$\text{DQPF} = \frac{\Gamma_{\text{water}}}{\rho_f}$$

tracks fluid production or consumption on the grid. In a closed box, the combined water mass of solid mineral lattice and pore fluid remains conserved:

$$\Delta m_{\text{mineral,water}} + \Delta m_{\text{pore,water}} = 0$$

# Silicate Rock Melting and Magma Ocean Dynamics

How does `Erebus.jl` treat rock that melts, and how does it carry the heat of a convecting magma ocean that its grid cannot resolve?
Rock melts linearly between a solidus and a liquidus, and the latent heat enters the heat equation as a larger apparent heat capacity inside that interval.
Melt lowers the rock viscosity, while the Stokes-Darcy flow solver keeps viscosity above a numerical floor, so the grid cannot resolve thin convective boundary layers in liquid magma.
A regularized effective thermal conductivity on molten markers stands for this sub-grid convective heat transport.
The sections below give the constitutive equations, the smoothing regularizations, and the configuration parameters.

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
`Erebus.jl` incorporates latent heat buffering through the apparent heat capacity method ([Gerya, 2019](../reference/bibliography.md#Primary-References), Section 16.6):

$$\rho c_{p,\text{eff}} = \rho c_p + \rho_s L_m \frac{\partial F_m}{\partial T}$$

Because $\Delta T_{\text{melt}} = T_l(P) - T_s(P)$ remains constant under the parallel Clapeyron shift, the derivative in the melting interval is constant:

$$\frac{\partial F_m}{\partial T} = \begin{cases} \frac{1}{T_l - T_s} & \text{if } T_s(P) < T < T_l(P) \\ 0.0 & \text{otherwise} \end{cases}$$

Integrating the effective heat capacity over the solidus-liquidus interval yields the exact total latent heat enthalpy:

$$\int_{T_s}^{T_l} \left(\rho c_{p,\text{eff}} - \rho c_p\right) \, dT = \rho_s L_m \int_{T_s}^{T_l} \frac{1}{T_l - T_s} \, dT = \rho_s L_m$$

This integral is exact, so in the limit of small steps a cell that crosses the whole interval takes up the full latent heat.
Erebus takes the heat capacity at the temperature of the start of each step and does no non-linear enthalpy iterations.
A cell that enters or leaves the interval within a step therefore takes the latent heat one step late, and its total latent heat can differ from $\rho_s L_m$ by up to the latent heat of one step.
This lag is an enthalpy error of first order in the time step, as the section on [thermal subcycles](discretization_numerics.md#3.-Thermal-Subcycles-and-Temperature-Change-Limit) explains.

**Where this is checked.**
Latent heat buffering in the complete time loop is checked on the [Melting Shell](../validation/melting_shell.md) validation page.

---

## Magma Suspension Rheology

The effective viscosity of partially molten rock depends strongly on the silicate melt fraction $F_m$.
`Erebus.jl` implements a two-regime rheological model based on [Costa et al. (2009)](../reference/bibliography.md#Silicate-Melting-and-Magma-Rheology) and [Gerya (2019, Section 16.6.2)](../reference/bibliography.md#Primary-References).

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

Following [Solomatov (2007)](../reference/bibliography.md#Silicate-Melting-and-Magma-Rheology), [Solomatov (2015)](../reference/bibliography.md#Silicate-Melting-and-Magma-Rheology), and boundary-layer scaling theory for turbulent convection ([Kraichnan, 1962](../reference/bibliography.md#Hydrothermal-Convection-and-Fluid-Transport); [Solomatov and Stevenson, 1993](../reference/bibliography.md#Silicate-Melting-and-Magma-Rheology); [Siggia, 1994](../reference/bibliography.md#Silicate-Melting-and-Magma-Rheology)), the convective heat flux $F_{\text{conv}}$ scales inversely with fluid viscosity:

$$F_{\text{conv}} \propto \left(\frac{1}{\eta}\right)^\beta$$

Because the numerical Stokes solver enforces an artificial viscosity floor, $\eta_{\text{num}} \ge \eta_{\min} \gg \eta_{\text{fluid}}$, the missing convective heat transport is represented on each rock marker as an effective turbulent thermal conductivity $k_{\text{turb}}$, clamped to admissible bounds $[k_{\text{floor}}, k_{\text{cutoff}}]$ and bounded below by $k_{\text{cond}}$:

$$k_{\text{turb}} = \max\left(k_{\text{cond}}, \, \text{clamp}\left(k_{\text{cond}} \left(\frac{\eta_{\text{num}}}{\eta_{\text{fluid}}}\right)^\beta, \, k_{\text{floor}}, \, k_{\text{cutoff}}\right)\right)$$

| Parameter | Symbol | Default Value | Physical Meaning |
|:---|:---|:---|:---|
| Background conductivity | $k_{\text{cond}}$ | $3.0\text{ W}/(\text{m}\cdot\text{K})$ at $\phi_d = 0$ | Marker porous rock conductivity (`k_porous(k_solid, phi_dry)` on `henke2012`, or `ktotal` with pore fluid) |
| Numerical marker viscosity | $\eta_{\text{num}}$ | $\ge 10^{12}\text{ Pa}\cdot\text{s}$ | Solver-clamped melt-weakened rock marker viscosity (`etatotalm`, using `etasolidmm` above `tmsolidphase`) |
| Reference fluid viscosity | $\eta_{\text{fluid}}$ | $\eta_{\text{num}} \to 100.0\text{ Pa}\cdot\text{s}$ | Mush-blended silicate liquid viscosity (`melting.eta_fluid_silicate` above $F_{\text{end}}$) |
| Turbulent scaling exponent | $\beta$ | $1/3 \approx 0.3333$ | Soft-turbulence Rayleigh scaling exponent ([Solomatov, 2007](../reference/bibliography.md#Silicate-Melting-and-Magma-Rheology); [Solomatov, 2015](../reference/bibliography.md#Silicate-Melting-and-Magma-Rheology)) |
| Conductivity bounds | $[k_{\text{floor}}, k_{\text{cutoff}}]$ | $[10^{-3}, 10^6]\text{ W}/(\text{m}\cdot\text{K})$ | Numerical floor (`melting.k_turb_floor`) and cap (`melting.k_turb_cutoff`) |

The scaling exponent $\beta = 1/3$ derives from asymptotic boundary layer theory in the soft turbulence regime, where convective heat transport is governed by boundary layer instability.
An exponent of $\beta = 1/2$ corresponds to square-root viscosity parameterizations, such as that employed in `i2elvis` ([Gerya, 2019](../reference/bibliography.md#Primary-References)).
`Erebus.jl` supports user-configured exponents via `melting.turb_exponent`.

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

$$\xi_F = \text{clamp}\left( \frac{F_m - F_{\text{start}}}{F_{\text{end}} - F_{\text{start}}}, 0.0, 1.0 \right)$$

$$w_F(\xi_F) = 3 \xi_F^2 - 2 \xi_F^3$$

The weight $w_F(\xi_F)$ satisfies $w_F(0) = 0$, $w_F(1) = 1$, and has vanishing first derivatives at both endpoints:

$$\left.\frac{dw_F}{d\xi_F}\right|_{\xi_F=0} = 0, \qquad \left.\frac{dw_F}{d\xi_F}\right|_{\xi_F=1} = 0$$

This gives $C^1$ smoothness of $k_{\text{eff}}(F_m)$ at the onset $F_m = F_{\text{start}} = 0.30$.
At three higher melt fractions the viscosities change slope: $F_m = \phi_{\text{crit}} = 0.40$, the $\eta_{\min}$ clamp near $F_m \approx 0.498$ at default parameters, and $F_m = F_{\text{end}} = 0.50$.
At these points $k_{\text{eff}}(F_m)$ stays continuous ($C^0$) and non-decreasing, with bounded piecewise-continuous derivatives.
The validation pages use "mush" in a narrower sense than the planetesimal differentiation overview above: there, the mush is rock whose peak melt fraction exceeds $F_{\text{start}}$, where the soft-turbulence conductivity begins to rise, and the mush radius is its outer boundary.

### 2. Logarithmic Geometric Interpolation

Because thermal conductivity varies over several orders of magnitude between conduction ($3\text{ W}/(\text{m}\cdot\text{K})$) and turbulent convection ($10^3\text{ to }10^5\text{ W}/(\text{m}\cdot\text{K})$), linear interpolation produces steep gradients.
`Erebus.jl` blends conductivity logarithmically and enforces both the floor/cutoff bounds and the lower bound $k_{\text{cond}}$:

$$\log_{10} \tilde{k} = (1 - w_{\text{total}}) \log_{10} k_{\text{cond}} + w_{\text{total}} \log_{10} k_{\text{turb}}$$

$$k_{\text{eff}} = \max\left(k_{\text{cond}}, \, \text{clamp}\left(10^{\log_{10} \tilde{k}}, \, k_{\text{floor}}, \, k_{\text{cutoff}}\right)\right)$$

This clamp ensures that convective enhancement never reduces marker conductivity below $k_{\text{cond}}$.

### 3. Viscosity Matching in the Mush Interval

To prevent $\eta_{\text{num}} < \eta_{\text{fluid}}$ at the onset of melting, the reference fluid viscosity $\eta_{\text{fluid}}$ equals the solver-clamped melt-weakened viscosity $\eta_{\text{num}}$ (`etasolidcur`) for $F_m \le F_{\text{start}}$, and blends logarithmically to `melting.eta_fluid_silicate` over $(F_{\text{start}}, F_{\text{end}}]$ using $\xi_F$ from Section 1:

$$\ln \eta_{\text{fluid}}(F_m) = (1 - \xi_F) \ln \eta_{\text{num}}(F_m) + \xi_F \ln \eta_{\text{fluid,silicate}}$$

This formulation maintains $\eta_{\text{num}} / \eta_{\text{fluid}} \ge 1.0$ throughout the entire transition window.

### 4. Quadratic Surface Thermal Contrast Weighting

To prevent convective heat enhancement from acting on markers near the surface reference temperature (for example when segregated melt reaches cold surface layers), a thermal contrast weight $w_T$ modulates the transition relative to $T_{\text{surface,ref}}$ (`melting.T_surface_ref`, default $300.0\text{ K}$) over `melting.dT_turb_min` ($\Delta T_{\min} = 10.0\text{ K}$):

$$\Delta T = \max(0.0, T - T_{\text{surface,ref}})$$

$$w_T = \left[ \text{clamp}\left( \frac{\Delta T}{\Delta T_{\min}}, 0.0, 1.0 \right) \right]^2$$

The quadratic exponent guarantees that $\partial w_T / \partial T \to 0$ continuously as $T \to T_{\text{surface,ref}}$, which eliminates gradient kinks at the surface boundary.
The complete transition weight is:

$$w_{\text{total}} = w_F(\xi_F) \cdot w_T$$

When $w_{\text{total}} = 0$ ($\Delta T = 0$ or $F_m \le F_{\text{start}}$), the routine returns $\max(k_{\text{cond}}, \text{clamp}(k_{\text{cond}}, k_{\text{floor}}, k_{\text{cutoff}}))$, which equals $k_{\text{cond}}$ whenever $k_{\text{cond}} \ge k_{\text{floor}}$.

**Where this is checked.**
Regularized soft-turbulence conductivity in the complete time loop is checked on the [Magma Ocean Conductivity](../validation/magma_ocean_conductivity.md) validation page.

---

## Two-Phase Silicate Melt Segregation and Magma Ascent

When internal temperatures exceed the silicate solidus, silicate rock partially melts.
Liquid silicate melt is less dense than the crystalline silicate matrix ($\Delta\rho = \rho_{\text{solid}} - \rho_{\text{melt}} > 0$).
This density contrast generates a buoyant body force that drives melt upward and outward toward the surface, while the solid crystalline matrix compacts or settles inward.

### 1. Governing Migration Regimes

Silicate melt migration spans two physical regimes separated by the rheologically critical melt fraction $\phi_{\text{crit}} \approx 0.40$:

1. **Porous Darcy Percolation ($F_m \le F_{\text{perc\_end}}$):**
   Below the disaggregation threshold, solid silicate grains form a continuous, compacting framework.
   Buoyant melt percolates through interconnected channels following Darcy's law ([McKenzie, 1984](../reference/bibliography.md#Silicate-Melting-and-Magma-Rheology)):

   $$\mathbf{v}_{\text{perc}} = -\frac{k_\phi}{\eta_{\text{melt}} F_m} \Delta\rho \, \mathbf{g}$$

   where the matrix permeability $k_\phi$ follows a McKenzie power-law relation modified with a residual melt retention threshold $\phi_{\text{residual}}$ and a post-disaggregation linear branch:

   $$k_\phi = \begin{cases}
   0, & F_m \le \phi_{\text{residual}} \\
   k_{\text{ref}} \left[ \frac{F_m - \phi_{\text{residual}}}{\phi_0 - \phi_{\text{residual}}} \right]^n, & \phi_{\text{residual}} < F_m \le \phi_{\text{crit}} \\
   k_{\text{crit}} + k_{\text{ref}} \left( \frac{F_m - \phi_{\text{crit}}}{\phi_0 - \phi_{\text{residual}}} \right), & F_m > \phi_{\text{crit}}
   \end{cases}$$

   When $F_m \le \phi_{\text{residual}}$, permeability vanishes ($k_\phi = 0$), retaining a trapped melt fraction within the crystalline matrix.

2. **Hindered Stokes Crystal Settling ($F_m \ge F_{\text{settle\_start}}$):**
   When the melt fraction exceeds the disaggregation threshold, the solid framework breaks down into a crystal suspension.
   Isolated silicate crystals settle downward through liquid magma under gravity, while buoyant melt ascends.
   The upward melt segregation velocity equals the downward crystal settling flux scaled by the [Richardson and Zaki (1954)](../reference/bibliography.md#Silicate-Melting-and-Magma-Rheology) hindered settling relation:

   $$\mathbf{v}_{\text{susp}} = \mathbf{v}_{\text{Stokes}} F_m^n$$

   $$\mathbf{v}_{\text{Stokes}} = -\frac{2 r_{\text{grain}}^2 \Delta\rho \, \mathbf{g}}{9 \eta_{\text{melt}}}$$

   As melt fraction approaches unity ($F_m \to 1.0$), hindered settling approaches the unhindered Stokes velocity limit.

3. **Smoothstep Regime Transition:**
   Across the intermediate mush window $[F_{\text{perc\_end}}, F_{\text{settle\_start}}]$, the segregation velocity transitions smoothly using a $C^1$ cubic Hermite smoothstep blend:

   $$\xi = \text{clamp}\left( \frac{F_m - F_{\text{perc\_end}}}{F_{\text{settle\_start}} - F_{\text{perc\_end}}}, 0.0, 1.0 \right)$$

   $$w(\xi) = 3 \xi^2 - 2 \xi^3$$

   $$\mathbf{v}_{\text{seg}} = (1 - w(\xi)) \mathbf{v}_{\text{perc}} + w(\xi) \mathbf{v}_{\text{susp}}$$

### 2. Operator-Split Conservative Drift-Flux Transport

Resolving simultaneous Stokes matrix flow, aqueous fluid Darcy flow, and silicate melt migration in a single monolithic matrix leads to ill-conditioning because of the 20-order-of-magnitude viscosity contrast between liquid magma and solid mantle rock.
`Erebus.jl` solves silicate melt segregation via an operator-split conservative drift-flux formulation on Eulerian cell volumes:

$$\frac{\partial (\rho_m F_m)}{\partial t} + \nabla \cdot \left(\rho_m F_m \mathbf{v}_{\text{seg}}\right) = -\dot{M}_{\text{cryst}}$$

To preserve 3D spherical volume-weighted melt mass $\sum_m F_m w_{3\text{D},m}$ to machine precision ($< 10^{-12}$ relative drift) and maintain numerical monotonicity:
- Each marker has a spherical shell weight $w_{3\text{D},m} = 2 r_m$ when `spherical_metric = true`, with cell weight $W_{i,j} = \sum_{m \in (i, j)} w_{3\text{D},m}$ and cell melt inventory $M_{\text{melt},i,j} = \sum_{m \in (i, j)} F_m w_{3\text{D},m}$.
- Explicit subcycling subdivides each global timestep $\Delta t$ into $N_{\text{sub}}$ subcycles constrained by the Courant-Friedrichs-Lewy (CFL) limit:

  $$\Delta t_{\text{sub}} \le C_{\text{cfl}} \min\left( \frac{\Delta x}{|v_{x,\text{seg}}|}, \frac{\Delta y}{|v_{y,\text{seg}}|} \right)$$

- Multidimensional donor-receiver flux limiters clamp outgoing fluxes so no donor cell drains below $\phi_{\text{residual}} W_{\text{donor}}$ or exceeds the physical packing ceiling $\phi_{\text{pack}} W_{\text{recv}}$.
- Subcycle flux increments are accumulated across staggered velocity faces and redistributed proportionally to $w_{3\text{D},m}$ among Lagrangian markers within each cell.

### 3. Energetics and Phase Coupling

Silicate melt migration couples to the thermal energy equation through three mechanisms accumulated onto staggered pressure nodes $(i+1, j+1)$:

1. **Gravitational Energy Dissipation:**
   Segregation releases gravitational potential energy as heat through viscous shear dissipation:

   $$\Psi_{\text{seg}} = \Delta\rho \, g \, F_m \, \|\mathbf{v}_{\text{seg}}\|$$

   This volumetric heating rate is accumulated into the Eulerian thermal source grid $Q_{\text{seg}}$ and passed to the implicit energy solver.

2. **Advective Sensible Heat Transport:**
   When `sensible_heat_transport = true`, ascending melt transfers sensible enthalpy $\Delta E_{\text{sens}} = \Delta m_{\text{melt}} c_{p,\text{melt}} (T_{\text{donor}} - T_{\text{recv}})$ across cell faces, crediting $+\Delta Q_{\text{sens}}$ to the receiver cell and $-\Delta Q_{\text{sens}}$ to the donor cell so net sensible energy is conserved.

3. **Subsolidus Crystallization Latent Heat:**
   When buoyant melt migrates into cooler regions where temperature drops below the local silicate solidus or thermodynamic equilibrium melt fraction ($T < T_{\text{sol}}$), excess liquid melt crystallizes up to the local thermal capacity ceiling $\Delta T_{\max} = T_{\text{liq}} - T_{\text{cell}}$.
   Crystallization releases latent heat of fusion $L_m$:

   $$\dot{Q}_{\text{lat}} = \rho_m \dot{F}_{\text{freeze}} L_m$$

   During subcycling, the local cell temperature $T_{\text{cell}}$ is updated dynamically by the latent temperature rise $\Delta T_{\text{lat}} = \dot{Q}_{\text{lat}} \Delta t / (\rho_s c_{p,s})$ so cumulative crystallization across subcycles never exceeds the sensible heat deficit to reach the liquidus.

### 4. Mantle Depletion and Enrichment Tracking

When melt segregates outward, the residual solid mantle becomes depleted in melt while shallower receiver regions become enriched.
`Erebus.jl` tracks cumulative melt extraction and arrival on each Lagrangian marker $m$ through the signed state variable $F_{\text{extract}}[m] \in [-1, 1]$ (positive for melt-depleted donor markers, negative for melt-enriched receiver markers).
The net local melt fraction is bounded by $[0.0, 1.0]$:

$$F_m[m] = \text{clamp}\left(F_{\text{eq}}(T, P) - F_{\text{extract}}[m], \, 0.0, \, 1.0\right)$$

This bidirectional tracking preserves segregated melt on receiver markers across subsequent thermal and mechanical steps while preventing remelting of depleted donor residues at constant temperature.

---

## Numerical Integration Architecture

The melting and soft turbulence system operates on both Lagrangian markers and Eulerian grid cells:
1. **Marker Update:**
   Marker temperatures and lithostatic pressures determine $F_m$, the solver-clamped melt-weakened viscosity `etatotalm`, the apparent heat capacity `rhocptotalm`, and the regularized conductivity `ktotalm` on each rock marker.
2. **Marker-to-Grid Interpolation:**
   Harmonic and arithmetic interpolation schemes transfer marker viscosities and densities to staggered velocity and pressure nodes.
   Marker conductivities `ktotalm` are averaged arithmetically onto staggered `KX` and `KY` faces.
   When `surface_temperature = :fixed`, the uncut rock-only face conductivities `KXu` and `KYu` are also retained for the spherical metric source.
3. **Implicit Energy Solver and Subgrid Relaxation:**
   The parabolic thermal energy equation is assembled into a linear system with variable apparent heat capacity and face conductivities and solved implicitly, followed by subgrid marker-node thermal relaxation using each marker's `ktotalm` and `rhocptotalm`.

---

## Configuration Controls

### Silicate Melting Parameters (`[melting]`)

Melting and turbulence parameters are set in the `[melting]` table:
- `active`: Enables silicate melting thermodynamics and apparent heat capacity buffering.
- `soft_turbulence`: Enables regularized sub-grid soft turbulence convection.
- `turb_exponent`: Convective heat flux scaling exponent ($\beta = 1/3$ for [Solomatov (2007)](../reference/bibliography.md#Silicate-Melting-and-Magma-Rheology) and [Solomatov (2015)](../reference/bibliography.md#Silicate-Melting-and-Magma-Rheology) soft turbulence; $\beta = 1/2$ for the square-root viscosity law of `i2elvis`).
- `eta_fluid_silicate`: Dynamic viscosity of liquid silicate magma ($100.0\text{ Pa}\cdot\text{s}$).
- `F_turb_start` and `F_turb_end`: Melt fraction bounds for smoothstep conductivity blending ($[0.30, 0.50]$).
- `T_surface_ref` and `dT_turb_min`: Reference surface temperature ($300.0\text{ K}$) and ramp width ($10.0\text{ K}$) of the thermal contrast weight $w_T$.
- `k_turb_floor` and `k_turb_cutoff`: Lower ($10^{-3}\text{ W}/(\text{m}\cdot\text{K})$) and upper ($10^6\text{ W}/(\text{m}\cdot\text{K})$) bounds on the effective conductivity.
- `T_solidus` and `T_liquidus`: Phase-dependent solidus and liquidus temperatures.
- `dpdt_clapeyron`: Pressure-dependent melting slope $dT/dP$ (default $0.0\text{ K/Pa}$; typical $1.2\times 10^{-7}\text{ K/Pa}$).

### Two-Phase Magma Transport Parameters (`[magma_transport]`)

Two-phase buoyant melt segregation is configured in the `[magma_transport]` table:
- `active`: Activates buoyant silicate melt segregation advection.
- `F_perc_end` and `F_settle_start`: Melt fraction bounds separating Darcy percolation and hindered Stokes crystal settling ($[0.35, 0.45]$).
- `k_melt_ref`: Reference matrix permeability ($1.0\times 10^{-11}\text{ m}^2$).
- `phi_residual`: Residual melt fraction retained in pores ($0.01$).
- `cfl_melt`: Sub-cycling CFL Courant number ($0.5$).

For complete schema details, default values, and data types, see the [Configuration Schema Reference](../reference/config_schema.md#melting).

# Telescoping Computational Domain for Growth to Lunar Mass

This page documents the mathematical formulation, grid transformation invariants, marker buffer replenishment mechanics, and physical verification benchmarks for the telescoping domain engine in `Erebus.jl`. The engine doubles spatial domain dimensions during planetary growth, maintaining constant grid cell resolution $dx$ from early planetesimals ($R \sim 20\text{ km}$) to lunar-mass bodies ($R \sim 1,737\text{ km}$).

---

## 1. Physical and Numerical Motivation

Simulating planetesimal evolution from seed bodies to protoplanets presents a severe multiscale spatial challenge:

1. **Spatial Scale Range:** Planetesimal seeds initiate accretion at radii $R \sim 20 - 50\text{ km}$, while oligarchic and pebble growth can grow embryos to lunar mass ($R \approx 1,737\text{ km}$, $M \approx 7.35 \times 10^{22}\text{ kg}$). This represents a factor of 35 to 87 increase in radius and 4 to 6 orders of magnitude in mass.
2. **Resolution Trade-Offs on Static Grids:** 
   - A static grid sized to contain the final lunar body ($x_{\text{size}} \approx 5,000\text{ km}$) with modest node counts ($N_x = 101$) yields cell resolution $dx \approx 50\text{ km}$. On such a grid, an initial seed of $R = 25\text{ km}$ spans less than a single cell.
   - Maintaining $dx = 1\text{ km}$ on a static box of $5,000\text{ km}$ requires $N_x = 5001$ nodes. The resulting 2D linear system has $N_{\text{dof}} \approx 2.5 \times 10^7$ degrees of freedom per timestep, which is computationally prohibitive for million-year evolutionary runs.
3. **Free Surface Boundary Proximity:** In the marker-in-cell method, the planetary surface is represented as an internal free surface bounded by low-viscosity sticky air (Gerya, 2019). If the planetary boundary approaches the outer computational boundary, artificial boundary stresses and spurious traction forces contaminate internal convective circulation and compaction flow.
4. **Telescoping Domain Solution:** The telescoping domain doubles domain dimensions ($x_{\text{size}}^{\text{new}} = 2 x_{\text{size}}$) and basic node counts ($N_x^{\text{new}} = 2(N_x - 1) + 1$) whenever the body radius exceeds a predefined fraction of the domain half-width ($R > 0.70 \cdot x_{\text{size}}/2$). This preserves exact cell spacing $dx = \text{const}$, preserves all interior marker positions and thermochemical invariants, and replenishes the newly created outer volume with sticky air markers.

---

## 2. Theoretical Formulation and Invariants

The mathematical principles of coordinate doubling, radial distance invariance $r_m^{\text{new}} = r_m^{\text{old}}$, odd-grid parity alignment ($N_x = 2k + 1$), and sticky-air buffer replenishment are derived in detail in [Telescoping Computational Domains](../explanations/telescoping_domain.md).

Key discrete transformations validated on this page include:

- **Invariant Grid Resolution ($dx = \text{const}$):**
  $$x_{\text{size}}^{\text{new}} = 2 \, x_{\text{size}}, \quad N_x^{\text{new}} = 2(N_x - 1) + 1 \implies dx^{\text{new}} = dx$$
- **Trigger Criterion:**
  $$R(t) > f_{\text{threshold}} \cdot \frac{x_{\text{size}}}{2}, \quad f_{\text{threshold}} = 0.70$$
- **Marker Coordinate Translation:**
  $$x_m^{\text{new}} = x_m^{\text{old}} + \frac{x_{\text{size}}}{2}, \quad y_m^{\text{new}} = y_m^{\text{old}} + \frac{y_{\text{size}}}{2} \implies r_m^{\text{new}} = r_m^{\text{old}}$$
- **Odd-Parity Grid Remapping:**
  $$j_{\text{off}} = \frac{N_x^{\text{new}} - N_x^{\text{old}}}{2} = k, \quad \Delta x_{\text{shift}} = j_{\text{off}} \, dx$$
  $$A^{\text{new}}[i_{\text{off}} + i, \, j_{\text{off}} + j] = A^{\text{old}}[i, j]$$
- **Sticky-Air Buffer Replenishment:** Injection of $n_{\text{sub}}$ neutral markers per cell in outer cells of area $3 A_{\text{old}}$ with phase type $tm = 3$ and ambient thermal properties.

---

## 3. Physical Conservation Invariants

Telescoping domain doubling satisfies strict physical conservation laws:

| Invariant | Discrete Conservation Law | Physical Guarantee |
|:---|:---|:---|
| **Solid Mass** | $M_{\text{solid}} = \sum_{m=1}^{N} m_m \, \delta_{tm_m \in \{1, 2\}} = \text{const}$ | No rock mass is created or destroyed during grid doubling. Preserved to machine precision ($< 10^{-15}$). |
| **Thermal Energy** | $E_{\text{thermal}} = \sum_{m=1}^{N} (\rho c_p V)_m \, T_m \, \delta_{tm_m \in \{1, 2\}} = \text{const}$ | Planetary internal heat content is strictly conserved. Added buffer markers carry independent ambient energy. |
| **Volatile Inventory** | $M_{\text{vol}} = \sum_{m=1}^{N} m_m \, X_{m,\text{species}} = \text{const}$ | Total inventories of $\text{H}_2\text{O}$, $\text{C}$, $\text{N}$, and $\text{S}$ are invariant. |
| **Core Metal Volume** | $V_{\text{metal}} = \sum_{m=1}^{N} V_m \, X_{\text{fe},m} = \text{const}$ | Metallic iron volume and segregation budgets are preserved. |
| **Radial Distance** | $r_m^{\text{new}} = r_m^{\text{old}} \quad \forall \; m \le N_{\text{old}}$ | Radial profiles of temperature, melt fraction, and density remain identical. |

---

## 4. Gravitational Poisson Solver Re-factorization

Self-gravitational potential $\Phi$ satisfies the 2D Poisson equation on the discrete staggered grid (Gerya, 2019):

$$\nabla^2 \Phi = \frac{8}{3} \pi G \rho_{\text{total}}$$

where the coefficient $8/3 \pi G = (2/3) \times 4 \pi G$ accounts for the 2D geometric approximation of spherical self-gravity. When grid node counts change from $(N_y, N_x)$ to $(N_y^{\text{new}}, N_x^{\text{new}})$:

1. **Operator Re-assembly:** The sparse 5-point discrete 2D Laplacian operator $L \in \mathbb{R}^{N_{\text{nodes}} \times N_{\text{nodes}}}$ is reassembled for the new dimensions:
   $$L_{(i,j),(i,j)} = -\left(\frac{2}{dx^2} + \frac{2}{dy^2}\right)$$
   $$L_{(i,j),(i\pm 1, j)} = \frac{1}{dy^2}, \qquad L_{(i,j),(i, j\pm 1)} = \frac{1}{dx^2}$$
2. **Dirichlet Boundary Re-indexing:** Boundary nodes along the computational box boundary and nodes outside the inscribed circle enforce homogeneous Dirichlet conditions:
   $$\Phi_{\partial\Omega} = 0$$
3. **LU Re-factorization:** Sparse LU factorization ($P L Q = L_U U_U$) is recomputed via UMFPACK or Pardiso. Re-factorization occurs once per domain doubling event, incurring negligible overhead during an evolutionary run.
4. **Acceleration Update:** Gravitational acceleration components $(g_x, g_y) = -\nabla \Phi$ are re-evaluated on the new grid to provide continuous gravitational body forces for Darcy percolation and Stokes diapirism.

---

## 5. Architectural Scaling and Profile Invariance Benchmarks

The 4-panel benchmark suite illustrates the spatial hierarchy, radial invariance, gravitational potential, and accretion trajectory of the telescoping domain engine:

![Telescoping domain benchmark](../assets/telescoping_domain_benchmark.png)

*Figure 1: Telescoping domain benchmark suite for growth from planetesimal seed ($R = 40\text{ km}$) to lunar radius ($R = 1,737\text{ km}$). Panel (a): Nested computational domain hierarchy across successive doubling levels ($x_{\text{size}} = 140, 280, 560\text{ km}$) with sticky-air buffer regions and the 70% domain threshold circle. Panel (b): Planetesimal radial temperature profile invariance $T(r)$, showing identical radial coordinates for continuum profiles and Lagrangian markers before and after domain translation. Panel (c): Gravitational potential $\Phi(r)$ and gravitational acceleration $g(r)$ across the enlarged computational domain, enforcing homogeneous Dirichlet boundary conditions ($\Phi = 0$ at the domain boundary). Panel (d): Planetesimal accretion growth trajectory to lunar mass ($R_{\text{lunar}} = 1,737\text{ km}$), displaying discrete domain doubling events triggered whenever $R(t) > 0.70 \cdot (x_{\text{size}} / 2)$.*

---

## 6. Configuration Example

The TOML configuration snippet below enables the telescoping domain engine:

```toml
[telescoping]
active = true
r_threshold_fraction = 0.70      # Trigger doubling when R > 0.70 * (xsize / 2)
max_telescope_levels = 10        # Maximum allowable doubling levels
target_radius = 1737000.0        # Target final radius (Moon: 1,737 km)
buffer_markers_per_cell = 4      # Sticky-air markers injected per outer cell

[accretion]
active = true
mode = "pebble_auto"
M_initial = 1.0e17               # Initial seed mass [kg]
R_initial = 20000.0              # Initial seed radius [m] (20 km)
M_target = 7.35e22               # Lunar mass [kg]
R_target = 1737000.0             # Lunar radius [m]
rho_bulk = 3340.0                # Bulk density [kg/m^3]
t_start_myr = 0.1                # Accretion onset time [Myr]
t_duration_myr = 3.0             # Accretion duration [Myr]
```

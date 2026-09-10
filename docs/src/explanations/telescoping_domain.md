# Physical and Numerical Principles of Telescoping Computational Domains

This section explains the physical principles, multiscale numerical challenges, staggered-grid transformation mechanics, and conservation invariants underlying telescoping computational domains in `Erebus.jl`.

---

## 1. The Multiscale Spatial Challenge in Planetary Accretion

Simulating planetesimal accretion from initial planetesimal seeds to lunar-mass planetary embryos encompasses a growth of several orders of magnitude in both mass and radius:

- **Initial Planetesimal Seeds:** Planetesimals forming via streaming instability or gravitational collapse initiate growth with characteristic radii $R \approx 20 - 50\text{ km}$ and masses $M \approx 10^{17} - 10^{18}\text{ kg}$.
- **Planetary Embryos and Moon-Sized Protoplanets:** Oligarchic planetesimal collisions and pebble accretion drive rapid accretion toward lunar dimensions ($R \approx 1,737\text{ km}$, $M \approx 7.35 \times 10^{22}\text{ kg}$).

### Limitations of Static Spatial Domains

In traditional numerical geodynamic simulations, the computational domain remains fixed throughout the calculation. For planetesimal growth over multiple spatial scales, static domains introduce an intractable dilemma:

1. **Resolution Degradation (Large Fixed Box):** If the box is sized to encompass the final lunar embryo ($x_{\text{size}} \ge 4,000\text{ km}$) with a computationally tractable grid ($N_x = 101$), the cell spacing is $dx \approx 40\text{ km}$. On such a grid, an initial seed of $R = 25\text{ km}$ is represented by a single cell. Consequently, interior convection, Darcy fluid percolation, core metal segregation, and cold lid formation cannot be resolved.
2. **Prohibitive Computational Cost (High-Resolution Static Box):** If high spatial resolution ($dx = 1\text{ km}$) is maintained over a static $4,000\text{ km}$ domain, the grid requires $N_x = N_y = 4,001$ nodes ($1.6 \times 10^7$ grid cells). In 2D hydromechanical solvers solving coupled Stokes, Darcy, and Poisson equations, the resulting linear systems require tens of millions of degrees of freedom. Solving these systems at each timestep over million-year evolutionary spans exceeds accessible compute budgets.
3. **Boundary Artifacts (Small Fixed Box):** If the domain is sized to resolve the initial planetesimal ($x_{\text{size}} = 140\text{ km}$), expanding planetary crust rapidly encroaches upon the outer domain boundaries. In the marker-in-cell technique, free planetary surfaces require a surrounding buffer of low-viscosity sticky air (Gerya, 2019). When the planetary surface approaches within a few cells of the computational boundary, artificial boundary traction, mirror stresses, and boundary reflection distort internal circulation and mantle compaction.

---

## 2. Mechanics of Dynamic Domain Doubling

The telescoping domain method resolves this multiscale dilemma by expanding the computational domain dynamically as the planetesimal accretes:

```
+---------------------------+        +-------------------------------------------------------+
|                           |        |                                                       |
|        Sticky Air         |        |                      Sticky Air                       |
|          (tm = 3)         |        |                       (tm = 3)                        |
|                           |        |                                                       |
|       +-----------+       |        |             +---------------------------+             |
|       | Planetary |       |  --->  |             |                           |             |
|       |   Body    |       |        |             |      Planetary Body       |             |
|       | (tm = 1,2)|       |        |             |         (tm = 1,2)        |             |
|       +-----------+       |        |             +---------------------------+             |
|                           |        |                                                       |
|       Level 0: Box 1X     |        |                    Level 1: Box 2X                    |
+---------------------------+        +-------------------------------------------------------+
```

### Invariant Grid Cell Spacing

When the planetesimal radius $R(t)$ exceeds 70% of the domain half-width ($R > 0.70 \cdot x_{\text{size}} / 2$), the domain executes a doubling transformation. The physical box dimensions and basic grid node counts double according to:

$$x_{\text{size}}^{\text{new}} = 2 \, x_{\text{size}}, \qquad y_{\text{size}}^{\text{new}} = 2 \, y_{\text{size}}$$

$$N_x^{\text{new}} = 2(N_x - 1) + 1, \qquad N_y^{\text{new}} = 2(N_y - 1) + 1$$

Because the number of grid intervals $N_x - 1$ doubles concurrently with physical box length, the spatial resolution is strictly invariant:

$$dx^{\text{new}} = \frac{x_{\text{size}}^{\text{new}}}{N_x^{\text{new}} - 1} = \frac{2 \, x_{\text{size}}}{2(N_x - 1)} = dx$$

This ensures that numerical diffusion, numerical dispersion, and boundary layer resolution remain constant through telescoping events.

---

## 3. Staggered Grid and Lagrangian Marker Invariance

### Centering Transformation and Grid Parity

To preserve planetary symmetry, the physical center of the planetesimal shifts from $(x_c^{\text{old}}, y_c^{\text{old}})$ to the geometric center of the doubled domain:

$$x_c^{\text{new}} = \frac{x_{\text{size}}^{\text{new}}}{2}, \qquad y_c^{\text{new}} = \frac{y_{\text{size}}^{\text{new}}}{2}$$

The translation vector is equal to the half-width of the original domain:

$$\Delta x_{\text{shift}} = x_c^{\text{new}} - x_c^{\text{old}} = \frac{x_{\text{size}}}{2}, \qquad \Delta y_{\text{shift}} = y_c^{\text{new}} - y_c^{\text{old}} = \frac{y_{\text{size}}}{2}$$

Exact alignment between continuous marker translations and discrete grid node blocks requires odd basic grid dimensions $N_x = 2k + 1$ and $N_y = 2m + 1$. Under this condition, the discrete node offset is $j_{\text{off}} = (N_x^{\text{new}} - N_x) / 2 = k$. Because $\Delta x_{\text{shift}} = (N_x - 1) dx / 2 = k \, dx = j_{\text{off}} \, dx$, marker positions and grid cell centers align with zero offset. If $N_x$ were even, a half-cell offset ($dx / 2$) would arise between marker centers and grid node blocks. For this reason, `Erebus.jl` requires and validates odd node counts for telescoping domains.

### Lagrangian Marker Distance Preservation

All existing Lagrangian markers undergo pure translation:

$$\mathbf{x}_m^{\text{new}} = \mathbf{x}_m^{\text{old}} + \Delta \mathbf{x}_{\text{shift}}$$

Because both marker coordinates and planetary center translate by identical displacement vectors, the radial distance of every marker from the center remains unchanged:

$$r_m^{\text{new}} = \|\mathbf{x}_m^{\text{new}} - \mathbf{x}_c^{\text{new}}\| = \|(\mathbf{x}_m^{\text{old}} + \Delta \mathbf{x}_{\text{shift}}) - (\mathbf{x}_c^{\text{old}} + \Delta \mathbf{x}_{\text{shift}})\| = \|\mathbf{x}_m^{\text{old}} - \mathbf{x}_c^{\text{old}}\| = r_m^{\text{old}}$$

This mathematical identity guarantees that radial lithostatic pressure, internal temperature profiles, melt fraction contours, and chemical layering are unaffected by domain doubling.

---

## 4. Sticky-Air Free Surface Mechanics

### Traction-Free Surface Condition

The marker-in-cell method simulates internal free surfaces using the sticky-air approach (Gerya, 2019; Crameri et al., 2012):

1. **Viscosity Contrast:** The sticky-air layer is assigned an effective viscosity $\eta_{\text{air}} \ll \eta_{\text{rock}}$, typically $\eta_{\text{air}} \sim 10^{16}\text{ Pa s}$ compared to silicate lithospheric viscosities $\eta_{\text{rock}} \sim 10^{20} - 10^{24}\text{ Pa s}$.
2. **Density Contrast:** Sticky air carries a nominal density $\rho_{\text{air}} = 1.0\text{ kg/m}^3 \ll \rho_{\text{rock}} \sim 3000\text{ kg/m}^3$.
3. **Traction Cancellation:** At the irregular interface between rock markers ($tm = 1, 2$) and sticky-air markers ($tm = 3$), shear stresses and normal stresses relax to near-zero values ($\sigma_{nt} \approx 0, \sigma_{nn} \approx 0$), accurately mimicking unconfined atmospheric or vacuum planetary surfaces.

### Buffer Marker Replenishment

When the computational box doubles in each dimension, the total area quadruples ($A_{\text{new}} = 4 A_{\text{old}}$). The outer annular zone of area $3 A_{\text{old}}$ contains no preexisting markers. To maintain uniform marker resolution, new sticky-air markers are synthesized:

- Each outer grid cell receives $n_{\text{buffer}} = \text{buffer\_markers\_per\_cell}$ markers placed on a regular sub-grid.
- Newly created buffer markers receive neutral ambient properties ($tm = 3$, $T = T_{\text{ambient}}$, $hr = 0$, zero volatiles, zero metal).
- Preexisting solid rock and core markers are unmodified; this preserves exact physical conservation.

---

## 5. Discrete Operator Re-factorization

Grid doubling alters the sparsity structure and dimensions of discrete differential operators:

### Gravitational Poisson Equation

Self-gravitational potential $\Phi$ satisfies the 2D Poisson equation with cylindrical/spherical geometric correction factor $2/3$ (Gerya, 2019):

$$\nabla^2 \Phi = \frac{8}{3} \pi G \rho_{\text{total}}$$

where $8/3 \pi G = (2/3) \times 4 \pi G$. On the expanded grid:
1. The 5-point discrete Laplacian matrix $L$ is reassembled with dimensions $(N_x^{\text{new}} N_y^{\text{new}}) \times (N_x^{\text{new}} N_y^{\text{new}})$.
2. Boundary potentials are updated to enforce homogeneous Dirichlet conditions along the computational box boundary and outside the inscribed circle:
   $$\Phi_{\partial\Omega} = 0$$
3. The sparse matrix $L$ is refactored via sparse LU decomposition (UMFPACK or Pardiso). Because domain doubling occurs only a few times throughout a multimillion-year simulation, the one-time factorization cost is negligible compared to regular timestepping.

---

## 6. Conservation Invariant Proofs

The telescoping domain algorithm guarantees exact conservation of fundamental physical quantities:

### 1. Solid Mass Invariance

The total solid mass is the sum over all rock and core markers:

$$M_{\text{solid}} = \sum_{m=1}^{N_{\text{old}}} m_m \, \delta_{tm_m \in \{1, 2\}}$$

Because added buffer markers carry phase type $tm = 3$ (sticky air) and original markers are preserved without modification:

$$M_{\text{solid}}^{\text{new}} = \sum_{m=1}^{N_{\text{new}}} m_m \, \delta_{tm_m \in \{1, 2\}} = \sum_{m=1}^{N_{\text{old}}} m_m \, \delta_{tm_m \in \{1, 2\}} + \sum_{m=N_{\text{old}}+1}^{N_{\text{new}}} m_m \, \underbrace{\delta_{tm_m \in \{1, 2\}}}_{= 0} = M_{\text{solid}}$$

Solid mass is conserved to machine precision.

### 2. Internal Thermal Energy Invariance

The thermal energy of the planetary body is defined as:

$$E_{\text{thermal}} = \sum_{m=1}^{N_{\text{old}}} (\rho c_p V)_m \, T_m \, \delta_{tm_m \in \{1, 2\}}$$

By identical reasoning, because newly injected markers carry $tm = 3$:

$$E_{\text{thermal}}^{\text{new}} = E_{\text{thermal}}$$

### 3. Volatile and Metallic Elemental Invariants

Total species inventories (water, carbon, nitrogen, sulfur, metallic iron) satisfy:

$$\mathcal{M}_{\text{species}} = \sum_{m=1}^{N} m_m \, X_{m,\text{species}}$$

Because buffer markers are initialized with zero volatile fractions and zero metallic iron ($X_{\text{species}} = 0$ for all newly injected markers), elemental inventories remain invariant during domain doubling.

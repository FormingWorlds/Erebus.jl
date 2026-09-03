# Discretization and Numerical Methods

This page explains the numerical algorithms and spatial discretization techniques employed in `Erebus.jl`.

---

## The Staggered Finite-Difference Grid

`Erebus.jl` uses a 2D Cartesian staggered grid arrangement to avoid spurious pressure-velocity checkerboard oscillations:

```text
    (i,j)---------Vy(i,j+1/2)---------(i,j+1)
      |                                  |
      |                                  |
   Vx(i+1/2,j)       P, Pf(i+1/2,j+1/2)  Vx(i+1/2,j+1)
      |              T, phi              |
      |                                  |
    (i+1,j)-------Vy(i+1,j+1/2)-------(i+1,j+1)
```

- **P Nodes (Cell Centers)**: Total solid pressure $P_t$, fluid pore pressure $P_f$, temperature $T$, porosity $\phi$, normal deviatoric stress $\sigma_{xx}, \sigma_{yy}$, and material phase fractions.
- **Vx Nodes (Vertical Cell Faces)**: Horizontal solid velocity $v_x$ and horizontal Darcy fluid discharge $q_{xD}$.
- **Vy Nodes (Horizontal Cell Faces)**: Vertical solid velocity $v_y$ and vertical Darcy fluid discharge $q_{yD}$.
- **Basic Nodes (Cell Vertices)**: Shear deviatoric stress $\sigma_{xy}$, shear strain rate $\dot{\varepsilon}_{xy}$, and effective shear viscosity $\eta$.

---

## The Marker-in-Cell (MIC) Method

Material properties, temperature, hydration state, and porosity are carried on millions of Lagrangian markers that advect through the Eulerian grid:

1. **Marker to Grid Interpolation**:
   Marker properties are projected onto adjacent staggered grid nodes using 4-point bilinear weighting:
   $$w_{ij} = \left(1 - \frac{|\Delta x|}{\Delta x_{\text{cell}}}\right) \left(1 - \frac{|\Delta y|}{\Delta y_{\text{cell}}}\right)$$

2. **Subgrid Diffusion**:
   To prevent unphysical short-wavelength oscillations, subgrid diffusion damping is applied to temperature and stress before interpolation.

3. **Marker Advection**:
   Markers move with the solid velocity field using a 4th-order Runge-Kutta scheme or 2nd-order predictor-corrector scheme. Marker velocity is weighted between interpolated staggered velocities and local nodal velocities (`vpratio = 1/3`).

---

## Monolithic Hydro-Mechanical System Assembly

The coupled Stokes-Darcy equations are assembled into a single monolithic sparse linear system:

$$\mathbf{L} \mathbf{u} = \mathbf{R}$$

where the unknown solution vector $\mathbf{u}$ contains five degrees of freedom per grid cell:
$$\mathbf{u} = \begin{bmatrix} v_x \\ v_y \\ P_t \\ q_{xD} \\ q_{yD} \\ P_f \end{bmatrix}$$

- **$v_x, v_y$ Rows**: Discretized momentum conservation in the $x$ and $y$ directions.
- **$P_t$ Row**: Discretized solid continuity with poroelastic compaction coupling.
- **$q_{xD}, q_{yD}$ Rows**: Discretized Darcy flux relations.
- **$P_f$ Row**: Discretized fluid continuity with matrix and fluid compressibility terms.

The sparse system is solved using direct sparse LU factorization (UMFPACK via `LinearSolve.jl` or Pardiso via `Pardiso.jl`).

---

## Non-Linear Iterations (Picard Loop)

Because effective viscosities $\eta(\dot{\varepsilon}_{\text{II}}, P_{\text{eff}})$, permeabilities $k_\phi(\phi)$, and bulk viscosities $\eta_\phi(\phi)$ depend non-linearly on the state variables, each timestep executes nested Picard iteration loops:

1. **Global Thermochemical Loop (`titer`)**:
   Updates temperature, phase changes, radiogenic decay, and mass/enthalpy source terms ($DMP$, $DHP$).
2. **Plastic Yielding Loop (`iplast`)**:
   Checks Mohr-Coulomb stress invariants $\sigma_{\text{II}}$ against yield surfaces, adjusting effective viscosity and bulk compaction viscosity until stress changes drop below `yerrmax`.

---

## Boundary Conditions

- **Mechanical**: Free-slip or no-slip solid boundary conditions on domain walls.
- **Hydraulic**: Draining ($P_f = 0$) or impermeable ($\partial P_f / \partial n = 0$) boundaries.
- **Sticky Air Layer**: The domain includes a low-density, low-viscosity buffer layer representing open space above the planetesimal surface, allowing the free surface of the planetesimal to deform naturally.

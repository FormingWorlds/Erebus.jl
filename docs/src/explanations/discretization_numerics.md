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

- **P Nodes (Cell Centers)**: Total solid pressure $P_t$, fluid pore pressure $P_f$, temperature $T$, porosity $\phi$, normal deviatoric stress $\sigma_{xx}$ (with $\sigma_{yy} = -\sigma_{xx}$ from the 2D traceless deviatoric relation), and material phase fractions.
- **Vx Nodes (Vertical Cell Faces)**: Horizontal solid velocity $v_x$ and horizontal Darcy fluid discharge $q_{xD}$.
- **Vy Nodes (Horizontal Cell Faces)**: Vertical solid velocity $v_y$ and vertical Darcy fluid discharge $q_{yD}$.
- **Basic Nodes (Cell Vertices)**: Shear deviatoric stress $\sigma_{xy}$, shear strain rate $\dot{\varepsilon}_{xy}$, and effective shear viscosity $\eta$.

---

## The Marker-in-Cell (MIC) Method

Material properties, temperature, hydration state, and porosity are carried on thousands of Lagrangian markers that advect through the Eulerian grid:

1. **Marker to Grid Interpolation**:
   Marker properties are projected onto adjacent staggered grid nodes using 4-point bilinear weighting:
   $$w_{ij} = \left(1 - \frac{|\Delta x|}{\Delta x_{\text{cell}}}\right) \left(1 - \frac{|\Delta y|}{\Delta y_{\text{cell}}}\right)$$

2. **Subgrid Diffusion**:
   Subgrid diffusion damping for temperature and stress is available as a configurable option (`dsubgridt`, `dsubgrids`), set to zero (disabled) by default.

3. **Marker Advection**:
   Markers move with the solid velocity field using a 4th-order Runge-Kutta scheme (`move_markers_rk4!`) with 2nd-order continuity-based parabolic spatial velocity corrections.

---

## Monolithic Hydro-Mechanical System Assembly

The coupled Stokes-Darcy equations are assembled into a single monolithic sparse linear system:

$$\mathbf{L} \mathbf{u} = \mathbf{R}$$

where the unknown solution vector $\mathbf{u}$ contains six scalar degrees of freedom per grid cell (two velocities, total pressure, two Darcy fluxes, and fluid pressure):
$$\mathbf{u} = \begin{bmatrix} v_x \\ v_y \\ P_t \\ q_{xD} \\ q_{yD} \\ P_f \end{bmatrix}$$

- **$v_x, v_y$ Rows**: Discretized momentum conservation in the $x$ and $y$ directions.
- **$P_t$ Row**: Discretized total continuity with poroelastic compaction coupling.
- **$q_{xD}, q_{yD}$ Rows**: Discretized Darcy flux relations.
- **$P_f$ Row**: Discretized fluid continuity with matrix and fluid compressibility terms.

The sparse system is solved using direct sparse LU factorization (UMFPACK via `LinearSolve.jl` or Pardiso via `Pardiso.jl`).

### Stencil Topologies

#### 1. Horizontal Momentum ($v_x$ Stencil)

```text
                           kvx-6
                            Vx₂
                             |
               kvy-6     ETA(i-1,j)   kvy+6⋅Ny1-6
                Vy₁      GGG(i-1,j)     Vy₃
                 *       SXY0(i-1,j)     *
                           basic₁
                            ETA₁                       
                            SXY₁
               ETAP(i,j)     |      ETAP(i,j+1)
               GGGP(i,j)     |      GGGP(i,j+1) 
   kvx-6⋅Ny1   SXX0(i,j)    kvx     SXX0(i,j+1)  kvx+6⋅Ny1
     Vx₁---------P₁---------Vx₃---------P₂---------Vx₅
                kpm          |        kpm+6⋅Ny1
               ETAP₁         |        ETAP₂
               SXX₁          |        SXX₂
                          ETA(i,j) 
                kvy       GGG(i,j)     kvy+6⋅Ny1
                Vy₂       SXY0(i,j)      Vy₄
                 *        basic₂          * 
                            ETA₂ 
                            SXY₂
                             |
                           kvx+6
                            Vx₄
                             *
```

#### 2. Vertical Momentum ($v_y$ Stencil)

```text
                           kvy-6
                            Vy₂
                             |
                          ETAP(i,j)
                          GGGP(i,j)
             kvx-6⋅Ny1    SXX0(i,j)     kvx
                Vx₁          P₁         Vx₃
                 *         ETAP₁         *
                            SYY₁
               ETA(i,j-1)   kpm       ETA(i,j)
               GGG(i,j-1)    |        GGG(i,j)
   kvy-6⋅Ny1   SXY0(i,j-1)  kvy       SXY0(i,j)  kvy+6⋅Ny1
     Vy₁-------basic₁-------Vy₃-------basic₂-------Vy₅
               ETA₁          |        ETA₂     
               SXY₁          |        SXY₂
                            kpm+6
                         ETAP(i+1,j)
                         GGGP(i+1,j)
          kvx-6⋅Ny1+6    SXX0(i+1,j)   kvx+6
                Vx₂          P₂        Vx₄
                 *         ETAP₂        *  
                            SYY₂
                             |
                           kvy+6
                            Vy₄
```

#### 3. Solid Continuity ($P_t$ Stencil)

```text
                 kvy-6
                  Vy₁
                   |
                   |
      kvx-6⋅Ny1   kpm       kvx
        Vx₁--------P--------Vx₂
                   |
                   |
                  kvy
                  Vy₂
```

#### 4. Darcy Fluid Continuity ($P_f$ Stencil)

```text
                 qyD₁
                kqy-6
                  |
       qxD₁-------P-------qxD₂
    kqx-6⋅Ny1    kpf      kqx
                  |
                 qyD₂
                 kqy
```

#### 5. Horizontal Darcy Flux ($q_x^D$ Stencil)

```text
       P₁--------qxD--------P₂
      kpf        kqx     kpf+6⋅Ny1
```

#### 6. Vertical Darcy Flux ($q_y^D$ Stencil)

```text
       P₁  kpf
       |
      qyD  kqy
       |
       P₂  kpf+6
```

---

## Analytical Darcy Elimination (4-Variable System)

To reduce computational overhead and memory footprint, the linear momentum equations for Darcy filtration:

$$R_x q_{xD} + \frac{\partial P_f}{\partial x} = \rho_f g_x$$

$$R_y q_{yD} + \frac{\partial P_f}{\partial y} = \rho_f g_y$$

can be substituted directly into the fluid mass conservation equation. Because the Darcy flux degrees of freedom appear with strictly diagonal drag coefficients $R_x = \eta_f / k_{\phi, x}$ and $R_y = \eta_f / k_{\phi, y}$, the fluxes $q_{xD}$ and $q_{yD}$ are eliminated analytically at $\mathcal{O}(N)$ operations prior to matrix assembly or operator evaluation:

$$q_{xD, ij} = \frac{1}{R_{x, ij}} \left(\rho_{f, ij} g_{x, ij} - \frac{P_{f, i, j+1} - P_{f, ij}}{\Delta x}\right)$$

$$q_{yD, ij} = \frac{1}{R_{y, ij}} \left(\rho_{f, ij} g_{y, ij} - \frac{P_{f, i+1, j} - P_{f, ij}}{\Delta y}\right)$$

Substituting these expressions into $\nabla \cdot \mathbf{q}_D$ yields a discrete 5-point Laplacian operator for $P_f$:

$$\nabla \cdot \mathbf{q}_D = \frac{q_{xD, ij} - q_{xD, i, j-1}}{\Delta x} + \frac{q_{yD, ij} - q_{yD, i-1, j}}{\Delta y}$$

This condenses the Stokes-Darcy system from six unknowns per cell ($v_x, v_y, P_t, q_{xD}, q_{yD}, P_f$) to four unknowns per cell ($v_x, v_y, P_t, P_f$):

$$\mathbf{u}_4 = \begin{bmatrix} v_x \\ v_y \\ P_t \\ P_f \end{bmatrix}$$

Following the solution of the 4-variable linear system, the Darcy velocity fields $q_{xD}$ and $q_{yD}$ are reconstructed in $\mathcal{O}(N)$ time with machine-precision fidelity via `reconstruct_darcy_fluxes!`.

---

## Matrix-Free Operators and Preconditioned Krylov Solvers

At grid resolutions exceeding $1024 \times 1024$ ($> 4 \times 10^6$ nodes, $> 1.6 \times 10^7$ degrees of freedom), storing a global sparse matrix in Compressed Sparse Column (CSC) format requires gigabytes of memory, and sparse direct factorization requires tens to hundreds of gigabytes of RAM.

`Erebus.jl` provides matrix-free operator evaluation (`MatrixFreeStokesDarcyOperator`) and iterative Krylov solvers (`solve_hydromechanical_iterative!`):

1. **Matrix-Free Operator Evaluation**:
   The operator action $\mathbf{y} = \mathbf{A} \mathbf{x}$ is evaluated directly on 2D grid property arrays without allocating sparse matrix indices or non-zero value arrays. Thread-parallel column execution yields high cache locality and eliminates memory allocation during solver iterations.

2. **Block-Schur Preconditioning**:
   Decouples the velocity and pressure blocks using a Schur complement approximation. The inverse diagonal of the velocity momentum block scales the velocity fields, while a discrete Laplacian approximation scales the pressure Schur complement block to maintain fast convergence over diverse permeability and viscosity regimes.

3. **Geometric Multigrid (GMG) Preconditioning**:
   At high grid resolutions, geometric multigrid provides fast error damping across spatial scales:
   - **Grid Hierarchy**: The grid spacing doubles at each coarser level ($\Delta x_{l+1} = 2 \Delta x_l, \Delta y_{l+1} = 2 \Delta y_l$). Rheological viscosities and hydraulic drag coarsen with harmonic averaging. Densities and body forces coarsen with arithmetic averaging.
   - **Restriction and Prolongation**: Volume-weighted conservative restriction ($R = \frac{1}{4} P^T$) transfers fine residuals to coarse grids. Continuous staggered prolongation prolongates coarse corrections to fine grids with boundary condition guards.
   - **Decoupled V-Cycles**: Damped Jacobi or Red-Black Gauss-Seidel relaxation smoothers operate separately on the elliptic velocity momentum block ($v_x, v_y$) and the elliptic fluid Darcy pressure block ($P_f$). Total solid pressure ($P_t$) scales with the discrete Schur complement diagonal.

4. **Krylov Solvers**:
   Supports Flexible GMRES (`fgmres`), restarted GMRES (`gmres`), and stabilized Bi-conjugate Gradient (`bicgstab`) through `LinearSolve.jl`.

5. **Hardware Acceleration (GPU via `KernelAbstractions.jl`)**:
   Device-agnostic parallel kernels support multi-threaded CPUs and hardware accelerators (NVIDIA GPUs via `CUDA.jl`, Apple Silicon via `Metal.jl`, and AMD GPUs via `AMDGPU.jl`):
   - **Operator Evaluation (`mul_device!`)**: Computes point stencils for momentum, continuity, and Darcy flow directly on device memory with configurable workgroup tiling (default 16x16).
   - **Operator Diagonal (`compute_operator_diagonal_device`)**: Evaluates diagonal elements in parallel on device for Jacobi relaxation and Schur complement scaling.
   - **Adjoint Restriction (`restrict_4var_device!`)**: Formulated as a coarse-grid gather operation where each thread computes one coarse cell from sixteen fine staggered nodes. This eliminates race conditions and atomic write instructions on GPUs.
   - **Continuous Prolongation (`prolongate_4var_device!`)**: Bilinear interpolation for velocity nodes ($v_x, v_y$) and piecewise-constant injection for pressure nodes ($P_t, P_f$) matching staggered boundary conditions.
   - **Relaxation Smoothers (`smooth_velocity_device!`, `smooth_darcy_device!`)**: Executes in-place relaxation sweeps (damped Jacobi or Red-Black Gauss-Seidel) using pre-allocated working buffers on each level to eliminate buffer reallocations during V-cycles. Note that on Julia's host CPU backend `KernelAbstractions.jl` allocates task partition contexts during launch, whereas hardware GPU backends enqueue directly to device streams without host heap allocations.
   - **Precision Considerations on Apple Silicon**: Apple Metal GPUs do not provide hardware double-precision (`Float64`) arithmetic. Simulations deployed to Apple Metal hardware must configure single-precision floating point (`Float32`). Transferring `Float64` arrays or operators to Metal devices raises an informative `ArgumentError`.

6. **2D Domain Decomposition and Distributed Memory (MPI)**:
   For extreme-resolution simulations on High-Performance Computing (HPC) clusters, `Erebus.jl` integrates domain decomposition via `MPI.jl` (v0.20):
   - **2D Cartesian Process Grid (`DistributedGridTopology2D`)**: Organizes compute nodes into a 2D Cartesian mesh ($P_y \times P_x$) using `MPI.Cart_create` with automatic surface-to-volume ratio optimization ($P_y \approx P_x \approx \sqrt{P}$). Identifies cardinal (North, South, East, West) and diagonal corner neighbors.
   - **Asynchronous Halo Exchange (`exchange_halos!`)**: Non-blocking peer-to-peer ghost cell communication using pre-allocated buffers (`HaloBuffer`) with non-blocking `MPI.Isend` and `MPI.Irecv!` followed by `MPI.Waitall`. Dimensional 2-step exchange ensures corners are communicated without diagonal message overhead.
   - **Distributed Matrix-Free Operator (`DistributedStokesDarcyOperator`)**: Evaluates $\mathbf{y} = \mathbf{A} \mathbf{x}$ for distributed subdomains with overlapped interior computation and boundary halo exchange.
   - **Distributed Krylov Solvers (`DistributedVector`)**: Overloaded inner products (`distributed_dot`) and norms (`distributed_norm`) using `MPI.Allreduce` enable standard Krylov solvers (`fgmres`, `gmres`, `bicgstab`) to scale over distributed compute nodes.
   - **Lagrangian Marker Migration (`migrate_markers!`)**: Markers crossing subdomain boundaries are routed to destination ranks via `MPI.Alltoall` count exchange and `MPI.Alltoallv!` payload transfer, enforcing 100% exact particle count and mass conservation among ranks.

---

## Non-Linear Iterations (Picard Loop)

Because effective viscosities $\eta(\dot{\varepsilon}_{\text{II}}, P_{\text{eff}})$, permeabilities $k_\phi(\phi)$, and bulk viscosities $\eta_\phi(\phi)$ depend non-linearly on the state variables, each timestep executes nested Picard iteration loops:

1. **Global Thermochemical Loop (`titer`)**:
   Updates temperature, phase changes, radiogenic decay, and mass/enthalpy source terms ($DMP$, $DHP$).
2. **Plastic Yielding Loop (`iplast`)**:
   Checks Mohr-Coulomb stress invariants $\sigma_{\text{II}}$ against yield surfaces, adjusting effective viscosity and bulk compaction viscosity until stress changes drop below `yerrmax`.

---

## Iteration Control and Numerical Limiters

To maintain numerical stability, `Erebus.jl` limits timesteps, plastic iterations, and field rates:

### 1. Plastic Yielding and Timestep Cut Retry

The plastic loop iterates until the maximum relative yield error drops below `yerrmax`:

$$\text{YERRNOD} = \max_i \frac{|\sigma_{\text{yield}} - \sigma_{\text{II}}|}{\sigma_{\text{yield}}} < \text{yerrmax}$$

If yielding nodes persist after `max_plastic_iterations` (default: `10000`):
- The solver rejects the candidate timestep.
- The system restores state from a start-of-step snapshot.
- The computational timestep halves ($dt \leftarrow 0.5 dt$).
- The step retries from the restored state.

If the run reaches `max_dt_reductions` (default: `5`) step cuts without plastic balance, the code throws a `PlasticConvergenceError`.

### 2. Porosity Change Rate Limiter

Rapid pore collapse or opening can destabilize Darcy fluid flow. The maximum displacement timestep limits relative porosity change per step:

$$\Delta t_{\phi} \le \frac{d\phi_{\text{max}}}{\max_{\text{interior}} |a\phi|}$$

where `dphimax` defaults to `0.1` (10% relative porosity change per step), and $a\phi$ is the fluid continuity divergence source term. The check omits the fixed pressure anchor cell $(i=2, j=2)$ to avoid false timestep limits at boundary cells.

### 3. Thermal Subcycles and Temperature Change Limit

Within the thermal solver, thermal subcycles advance with local timestep $dt_t \le dt$. On every substep, the solver tracks the maximum temperature change:

$$\Delta T_{\text{max}} = \max |T - T_{\text{old}}|$$

If $\Delta T_{\text{max}} > \text{DTmax}$ (default: `20.0 K`), the subcycle timestep scales down:

$$dt_t \leftarrow dt_t \frac{\text{DTmax}}{\Delta T_{\text{max}}}$$

This prevents thermal spikes from concentrated radiogenic decay or rapid phase change.

### 4. Hydrofracture Permeability Bounds

Dynamic hydrofracturing increases matrix permeability when pore fluid pressure exceeds the minimum compressive stress plus tensile rock strength:

$$k_{\text{eff}} = \min\left( \max(k_{\text{enhanced}}, k_\phi), \max(k_\phi, k_{\text{frac\_max}}) \right)$$

This two-sided clamp ensures that enhanced permeability cannot fall below matrix permeability $k_\phi$ or exceed ceiling $k_{\text{frac\_max}}$.

---

## Boundary Conditions

- **Mechanical**: Free-slip solid boundary conditions on the outer domain boundaries.
- **Hydraulic**: Draining ($P_f = p_{\text{surface}} = 1000\text{ Pa}$) pore pressure boundary anchors on outer walls.
- **Sticky Air Layer**: The domain includes a low-density, low-viscosity buffer layer representing open space above the planetesimal surface, allowing the free surface of the planetesimal to deform naturally.

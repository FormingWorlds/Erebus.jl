# Spherical Momentum

## What this page checks

This page checks the momentum rows of the hydromechanical solve of `Erebus.jl` with the option `geometry.spherical_flow = true`.
With the option on, the 2-D grid stands for a slice through a sphere, and the momentum rows get metric terms that make the stress balance that of the sphere.
The function `_add_spherical_momentum!` in `src/numerics/stokes_darcy.jl` adds these terms, and both assembly paths of the solver call it (lines 701 and 1515).
Validation means that we run the code on a case where we know the right answer from another source, and we check that the code gives that answer.
If these terms were wrong, the flow of a deforming planetesimal interior and the stresses that drive its fracturing would be wrong.

The page checks 3 groups of terms, each with an effect of order 1 on the flow:

- the weighted stress fluxes of the stress divergence, which every shear flow reaches;
- the out-of-plane stress $S_{33}$ and its hoop term, which only a flow with a radial velocity reaches;
- the stress update that writes `grids.S33` after the solve.

## The momentum equation of the slice

With the option on, the momentum rows solve

$$\frac{1}{w} \partial_j (w \tau_{ij}) - \tau_{33} \frac{\partial_i w}{w} - \partial_i P_t + f_i = 0, \qquad i, j \in \{x, y\}, \qquad w = 2r,$$

with the body force $f_i = \rho g_i$ and the distance $r$ from the centre of the sphere.
The stress is the 3-D deviatoric stress $\tau_{ij} = 2 \eta (\dot\varepsilon_{ij} - \delta_{ij} \nabla_3 \cdot \mathbf{u} / 3)$ with the out-of-plane strain rate $\dot\varepsilon_{33} = \mathbf{u} \cdot \nabla w / w$.
This is the momentum equation of the spherical annulus of Hernlund and Tackley (2008), their Eqs. (2) and (3) with $d = 2$, written in Cartesian components (see [Discretization and Numerics](../explanations/discretization_numerics.md)).
With the option off the rows solve the 2-D equation $\partial_j \tau_{ij} - \partial_i P_t + f_i = 0$.

The 2-D rows hold the 2-D terms, and `_add_spherical_momentum!` adds the difference:

- line 754 computes $\dot\varepsilon_{33}$ at a P node as the weighted divergence of the continuity row minus the 2-D divergence, with the face multipliers $m_F = w_F / w_P$;
- line 756 computes $S_{33} = \tau_{33} = \frac{2\eta}{3} (2 \dot\varepsilon_{33} - \dot\varepsilon_{xx} - \dot\varepsilon_{yy})$;
- lines 789 to 792 (x row) and 817 to 820 (y row) add $(m_F - 1)$ times the stress flux through each face, the difference between $\frac{1}{w} \partial_j (w \tau_{ij})$ and $\partial_j \tau_{ij}$;
- lines 793 and 794 (x row) and 821 and 822 (y row) add the $S_{33}$ part of the weighted divergence, with $\tau_{xx} = S_{xx} - S_{33} / 2$, and the hoop term $\tau_{33} \partial_i w / w$ with $\partial_x w / w = (m_E - m_W) / \Delta x$ from lines 784 and 811.

After the solve, `_spherical_stress_update!` in `src/numerics/plasticity.jl` (line 1229) writes `grids.S33` from the solved velocity.

## The case

The case uses a sphere of radius $R = 50$ km in a square box of 140 km with the centre in the middle of the box.
`Nx` is the number of grid nodes on each side of the box; it is odd, so a basic node lies on the centre.
The whole box has the shear viscosity $\eta = 10^{20}$ Pa s, with no air layer.
The matrix is nearly incompressible: the ratio of the shear to the compaction viscosity is $\eta / \zeta = 10^{-6}$.
The shear modulus is $10^{35}$ Pa, so the elastic term leaves the viscosity unchanged, and the stress history is 0.
The fluid density is 0, and the case sets the body force $f = \rho g$ directly in the gravity fields `gx` and `gy`.
The total pressure is held at $P_\mathrm{amb} = 100$ Pa on the box edge.
Every exact flow below is 0 for $r \ge R$, so the free-slip walls and the pressure ring of the box edge see the exact solution.

The case builds the state with `init_simulation`, sets the grid fields directly and calls `solve_stokes_darcy!` once on timestep 1.
With the force scale $F_0 = 100$ N/m³ the velocity unit of case A is $F_0 R^2 / \eta = 2.5 \times 10^{-9}$ m/s, and its largest velocity is a few $10^{-12}$ m/s.

## Case A: azimuthal shear, exact with and without the option

**Why this check exists.**
The weighted stress fluxes $(m_F - 1)$ could be missing, have the wrong sign or use the weight of the wrong face.
For an azimuthal shear flow they change the balance of the shear stress at order 1, so each of these faults changes the flow.

**What it tests.**
The force $f_\theta = F_0 \rho (1 - \rho^2)^2 (1 - 32 \rho^2 / 5 + 39 \rho^4 / 5)$ with $\rho = r / R$ drives the azimuthal flow $\mathbf{u} = r \Omega(r) \mathbf{e}_\theta$.
We compare the velocity at all interior nodes of the box with the exact flow for `spherical_flow = true` and `false`, at `Nx = 33`, 65 and 129.
The L2 error is the relative discrete L2 error over all interior velocity nodes with equal weights.
We also compare the rotation rate at the centre and the L2 norms of the 2 profiles, in the sections below.

**How we know the right answer.**
For an azimuthal flow the divergence is 0, $\tau_{33} = 0$, the pressure is uniform and the only stress is $\tau_{r\theta} = \eta r \Omega'$.
Hernlund and Tackley (2008) give the angular momentum equation of the annulus as their Eq. (3).
With $\partial_\phi = 0$ and $v_r = 0$ it reduces with $d = 2$ (the sphere, option on) to

$$\frac{1}{r^3} \frac{d}{dr} \left( r^3 \tau_{r\theta} \right) + f_\theta = 0,$$

and with $d = 1$ (the cylinder, which for this flow is the 2-D equation of the option off) to $\frac{1}{r^2} \frac{d}{dr} (r^2 \tau_{r\theta}) + f_\theta = 0$.
The force has $\int_0^R s^2 f_\theta \, ds = \int_0^R s^3 f_\theta \, ds = 0$, so $\Omega' = 0$ at $r = R$ in both forms and the flow is 0 outside.
The exact solutions are

$$\Omega_\mathrm{on} = \frac{F_0 R}{\eta} \frac{(1 - \rho^2)^4 (1 - 6 \rho^2)}{100}, \qquad \Omega_\mathrm{off} = \frac{F_0 R}{\eta} \frac{(1 - \rho^2)^4 (3 - 13 \rho^2)}{200}.$$

The test checks both closed forms against their equation and the 2 moments of the force with exact rational arithmetic, and the Float64 forms against finite differences.

The option-on equation is also the equation of the full 3-D sphere for this flow.
The 3-D flow $\mathbf{u} = \Omega(r) \, r \sin\vartheta \, \mathbf{e}_\varphi$ with the polar angle $\vartheta$ and the force $f_\varphi = \sin\vartheta \, f_\theta(r)$ has the only stress component $\tau_{r\varphi} = \eta \, r \sin\vartheta \, \Omega'(r)$, because $\tau_{\vartheta\varphi} = \eta \sin\vartheta \, \partial_\vartheta \Omega = 0$.
The $\varphi$ component of the stress divergence in spherical coordinates is then $\frac{1}{r^3} \partial_r (r^3 \tau_{r\varphi})$, and the factor $\sin\vartheta$ cancels against that of the force.
So the equator of the 3-D toroidal flow, and every other latitude, satisfies the option-on equation above, and case A tests the slice against a flow of the full sphere.

**Why it matters.**
Every run with `spherical_flow = true` uses the weighted stress fluxes in each momentum row.

## Case A0: no rotation at the centre

Case A0 has $\Omega = -\frac{F_0 R}{\eta} \frac{\rho^2 (1 - \rho^2)^4}{10}$ with the force $f_\theta = F_0 \rho (1 - \rho^2)^2 (1 - 46 \rho^2 / 5 + 13 \rho^4)$ of the option-on equation, so $\Omega(0) = 0$.
The test checks its closed form in the same way as case A.
Case A0 runs in the test only, at `Nx = 33` and 65, and has no record row; it shows the cause of the slow centre convergence below.

## Case B: radial flow and the out-of-plane stress

**Why this check exists.**
The out-of-plane stress $S_{33}$ and its hoop term are 0 in case A, so case A cannot see them.
The prefactor $2/3$ of $S_{33}$, the factor 2 of $\dot\varepsilon_{33}$, the factor $1/2$ of the hoop term and the factor $1/2$ of the $S_{33}$ flux could be wrong, and the stress update could write a wrong `grids.S33`.

**What it tests.**
Case B is a manufactured flow with a radial velocity.
With $X = (x - x_c) / R$, $Y = (y - y_c) / R$ and the potential $\phi = X Y (1 - \rho^2)^6$, the velocity is $r \mathbf{u} = \nabla \times (r^3 \phi \, \mathbf{e}_z)$ in units of $U_0 = 10^{-9}$ m/s:

$$u_x = 3 Y \phi + \rho^2 \partial_Y \phi, \qquad u_y = -(3 X \phi + \rho^2 \partial_X \phi).$$

It has $\partial_j (r u_j) = 0$, so the 3-D divergence of the slice is 0, and $\tau_{33} = 2 \eta u_r / r$ is of order $\eta U_0 / R$.
The flow and its first 4 derivatives are 0 at $r = R$.
The body force is $f = -D(\mathbf{u})$, with $D$ the stress divergence of the option-on equation above, computed from the exact polynomials in rational arithmetic and not from the code.
We compare the velocity with the exact flow over the box and over the nodes with $r < R$, and `grids.S33` with the exact $\tau_{33}$ over the P nodes with $r < R$, at `Nx = 33`, 65 and 129.
Case B has no exact solution with the option off.

**How we know the right answer.**
The force is a polynomial in $X$ and $Y$ plus a remainder of degree 3 divided by $\rho^2$, which is at most a constant times $\rho$ and makes the force 0 at the centre.
The test evaluates the force a second time by finite differences of the exact velocity through the equation in Cartesian components, and a third time from Eqs. (2) to (7) of Hernlund and Tackley (2008) with $d = 2$ in polar coordinates; both agree with the force to a relative $10^{-6}$, and the form with $d = 1$ does not.
Case B therefore checks the documented momentum equation of the spherical annulus; the case does not check it against the full 3-D equations on a sphere.

**Why it matters.**
Every flow with a radial component, such as the compaction of a planetesimal interior or a plume, reaches the $S_{33}$ terms.

```@raw html
<div class="validation-fig-light">
```
![Azimuthal velocity of case A against radius on the axes through the centre at Nx = 65 with the exact profiles with and without the option, and the errors of the record against Nx with a second-order guide.](../assets/validation_momentum_profile_light.svg)
```@raw html
</div>
<div class="validation-fig-dark">
```
![Azimuthal velocity of case A against radius on the axes through the centre at Nx = 65 with the exact profiles with and without the option, and the errors of the record against Nx with a second-order guide.](../assets/validation_momentum_profile_dark.svg)
```@raw html
</div>
```

**How to read the plot.**
Panel (a) shows the azimuthal velocity of case A against $r / R$ at `Nx = 65` on the 2 axes through the centre, where the flow is purely azimuthal.
The solid line is the exact solution with the option on and the dashed line the exact solution with the option off.
The points are grid nodes.
A pass means that the points with the option on follow the solid line and the points with the option off follow the dashed line.
Panel (b) shows the errors of the record against `Nx` on logarithmic axes, with a dotted line of second order.
A pass means that the errors stay below their limits and fall along the dotted line, except the error of the centre node, which the next section explains.

## The rotation rate at the centre

The rotation rate at the centre is 1/100 with the option on and 3/200 without, in units of $F_0 R / \eta$, so the exact ratio is 2/3.
We measure it in 2 ways.
The fit takes $\Omega = v_\theta / r$ at the velocity nodes on the 2 axes through the centre with $0.18 \le \rho \le 0.85$ and fits a polynomial of degree 5 in $\rho^2$ by least squares.
The exact $\Omega$ has degree 5 in $\rho^2$, so the fit has no bias.
The fitted centre value converges at order 2.06 from `Nx = 33` to 65 in both settings, and the ratio of the fitted values is 2/3 within 0.08 % at `Nx = 65`.
The centre node value is half the vorticity at the basic node on the centre.

**Known limit.**
With the option on, the error of the centre node falls only at order 0.45 from `Nx = 33` to 65 and at order 0.74 from 65 to 129 (from the record rows `centre_node_a_on`), while without the option it falls at second order.
The record therefore holds the centre node error as a bound only, with no order and no ratio row.
The cause is the centre of the slice.
The weight $w = 2 \sqrt{x^2 + y^2}$ has a cone at $r = 0$, and within a few cells of the centre the truncation error of the weighted divergence does not fall with the grid spacing, as the docstring of `flow_metric_face_weights` states.
A flow with $\Omega(0) \ne 0$ is a rigid rotation near the centre, which has no strain in the continuum.
On the grid the stress update gives this rotation the out-of-plane strain rate $\dot\varepsilon_{33} = 0.011 \, \Omega$ at the P nodes at $r = 1.58 \Delta x$, the same on every grid, while away from the centre it falls at second order.
Case A has a total pressure error of about 1100 Pa next to the centre with the option on that does not fall with the grid, against a stress scale $F_0 R = 5 \times 10^6$ Pa.
Case A0, with $\Omega(0) = 0$, converges at second order in the L2 error, in the largest velocity error and in the largest pressure error.
The error stays near the centre: the L2 error of case A, which sums over the whole box, and the fitted centre value, which reads the nodes with $\rho \ge 0.18$, converge at second order.

## The L2 norms of the 2 profiles

The ratio of the L2 norms of the case A profiles with and without the option, $\sqrt{\int_0^1 \rho \, u_\mathrm{on}^2 \, d\rho / \int_0^1 \rho \, u_\mathrm{off}^2 \, d\rho}$, is exactly $\sqrt{8 / 11} = 0.853$, from the closed forms in rational arithmetic.
The discrete ratio over the velocity nodes with $r < R$ is compared with it.

## Conditions of the case

The exact solutions hold only if the case meets 6 conditions.
Each condition is a test assertion, so a change of the set-up that breaks one fails the test and not only the comparison.

- C1, incompressible matrix: after the solve, on every interior P node, $\eta / \zeta \le 10^{-6}$.
- C2, viscous limit: on every basic node and every P node the shear modulus times $\Delta t$ is at least $10^6$ times the viscosity.
- C3, body force: the solve leaves `gx` and `gy` unchanged.
- C4, force moments: the moments $\int s^2 f_\theta$ and $\int s^3 f_\theta$ of the case A force are 0 in exact rational arithmetic.
- C5, compact flow: the exact velocity of every case is 0 on $r = R$, checked in rational arithmetic at rational points of the circle and to round-off in Float64.
- C6, one solve: the time step after the solve equals $\Delta t$ to a relative $10^{-12}$, with one plastic pass that converges.

!!! details "Test details"
    - Reference: the closed forms of case A with and without the option; the manufactured case B with its force from the polynomials of the option-on equation; $\tau_{33} = 2 \eta u_r / r$; the exact centre rotation ratio 2/3; the exact profile norm ratio $\sqrt{8 / 11}$; equal fields of the 6-variable and the 4-variable solver paths.
    - Tolerance: each error at most about 1.4 times its value on the committed solver (the table below); each order between `Nx = 65` and `Nx = 129` equal to 2 within 0.35; the fitted centre ratio equal to 2/3 within a relative $8 \times 10^{-4}$ at `Nx = 65` and $2.5 \times 10^{-4}$ at `Nx = 129`; the profile norm ratio equal to $\sqrt{8 / 11}$ within a relative $7.5 \times 10^{-4}$ and $2.1 \times 10^{-4}$; the difference of the 2 solver paths at most $10^{-7}$, with the pressure difference scaled by $F_0 R$ (case A) or $\eta U_0 / R$ (case B).
    - Regression band: the allowed range around the recorded value is `PHYSICS_REGRESSION_BANDS` in `benchmarks/validation_spherical_momentum.jl`.
      It is $2 \times 10^{-3}$, $5 \times 10^{-4}$ and $1.5 \times 10^{-4}$ on each error at `Nx = 33`, 65 and 129, 0.05 on each order, $2 \times 10^{-4}$ and $10^{-4}$ on each ratio at `Nx = 65` and 129, and $10^{-7}$ on the solver path difference.
      Measured against the exact flow with the option on, the run with the option off has the L2 error 0.63 at `Nx = 33`; the test checks that this run leaves the bands of the option-on rows by more than twice the band.
    - Test id: `l2_a_{on,off}_nx{33,65,129}`, `l2_b_nx{33,65,129}`, `l2_inside_b_nx{33,65,129}`, `s33_b_nx{33,65,129}`, `centre_fit_a_{on,off}_nx{33,65,129}`, `centre_node_a_on_nx{33,65,129}`, `l2_order_a_{on,off}_65_129`, `l2_order_b_65_129`, `l2_inside_order_b_65_129`, `s33_order_b_65_129`, `centre_fit_order_a_{on,off}_65_129`, `centre_fit_ratio_nx{65,129}`, `profile_ratio_nx{65,129}`, `elim_{a,b}_nx33`.
      `Pkg.test` runs `test/test_validation_spherical_momentum.jl` at `Nx = 33` and 65 and pins the L2 error of case A with the option at `Nx = 65`, 8.6724558541e-3, to a relative $10^{-6}$.
    - Command: `julia --project=. tools/check_validation.jl --check`.

## What this case does not test

- The volumetric trace $\nabla_3 \cdot \mathbf{u} / 3$ in $S_{33}$: with $\eta / \zeta = 10^{-6}$ the divergence is close to 0, so a fault in the trace changes the flow far less than the grid error does.
  The operator test `Smooth non-radial flow with variable viscosity` in `test/test_spherical_momentum.jl` uses a flow with a nonzero divergence and checks the rows against the continuum stress divergence that includes the trace.
- The visco-elastic factors and the stress history `SXX0`, `SXY0` and `S330`: the shear modulus is $10^{35}$ Pa and the history is 0.
  The tests `Out-of-plane stress history in the momentum right-hand side`, `Stress history of all components in both momentum rows` and `Pressure stays unweighted and the visco-elastic factor scales the stress` in `test/test_spherical_momentum.jl` cover them.
- A variable viscosity `ETA` and `ETAP`, a grid with $\Delta x \ne \Delta y$ and a centre off the grid nodes: the case has a uniform viscosity, a square grid and the centre on a basic node.
  The operator tests of `test/test_spherical_momentum.jl` use a variable viscosity, $\Delta x \ne \Delta y$ and a centre between the nodes.
- A compressible matrix: there is no case with a large $\eta / \zeta$.
- Gravity of the planet, buoyancy, Darcy flow, plasticity and the time loop: the case sets the force directly, has no fluid density and solves once.
- How the code computes the viscosity and the force from the markers: the case sets the grid fields directly.

## All recorded values

The table comes from the record and must match it.
To make the figure again, run these 2 commands from the repository root:

```
julia --project=. scripts/export_validation_momentum_data.jl <data_dir>
python3 scripts/generate_validation_momentum_figures.py <data_dir>
```

The second command needs NumPy and Matplotlib.

```@raw html
<!-- VALIDATION_RECORD_START: spherical_momentum -->
```

**Validation record**: `benchmarks/records/spherical_momentum.toml` (Source commit `b80632a7`, config SHA-256 `f31389014dd1`)

| Case | Metric | Grid `Nx` | Measured | Reference | Tolerance | Units | Status |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| `A` | Velocity L2 relative error over the box, case A, option on (Nx=33) | 33 | `3.7452e-02` | `0.0000e+00` | `<= 5.30e-02` | dimensionless | **PASS** |
| `A` | Velocity L2 relative error over the box, case A, option on (Nx=65) | 65 | `8.6725e-03` | `0.0000e+00` | `<= 1.22e-02` | dimensionless | **PASS** |
| `A` | Velocity L2 relative error over the box, case A, option on (Nx=129) | 129 | `2.1555e-03` | `0.0000e+00` | `<= 3.10e-03` | dimensionless | **PASS** |
| `A` | Velocity L2 relative error over the box, case A, option off (Nx=33) | 33 | `3.7361e-02` | `0.0000e+00` | `<= 5.30e-02` | dimensionless | **PASS** |
| `A` | Velocity L2 relative error over the box, case A, option off (Nx=65) | 65 | `8.5455e-03` | `0.0000e+00` | `<= 1.20e-02` | dimensionless | **PASS** |
| `A` | Velocity L2 relative error over the box, case A, option off (Nx=129) | 129 | `2.1173e-03` | `0.0000e+00` | `<= 3.00e-03` | dimensionless | **PASS** |
| `B` | Velocity L2 relative error over the box, case B, option on (Nx=33) | 33 | `4.0409e-02` | `0.0000e+00` | `<= 5.70e-02` | dimensionless | **PASS** |
| `B` | Velocity L2 relative error over the box, case B, option on (Nx=65) | 65 | `9.8283e-03` | `0.0000e+00` | `<= 1.38e-02` | dimensionless | **PASS** |
| `B` | Velocity L2 relative error over the box, case B, option on (Nx=129) | 129 | `2.4411e-03` | `0.0000e+00` | `<= 3.50e-03` | dimensionless | **PASS** |
| `B` | Velocity L2 relative error from the nodes with r < R, case B, option on (Nx=33) | 33 | `4.0340e-02` | `0.0000e+00` | `<= 5.70e-02` | dimensionless | **PASS** |
| `B` | Velocity L2 relative error from the nodes with r < R, case B, option on (Nx=65) | 65 | `9.8063e-03` | `0.0000e+00` | `<= 1.38e-02` | dimensionless | **PASS** |
| `B` | Velocity L2 relative error from the nodes with r < R, case B, option on (Nx=129) | 129 | `2.4355e-03` | `0.0000e+00` | `<= 3.50e-03` | dimensionless | **PASS** |
| `B` | S33 L2 relative error against tau_33 = 2 eta u_r / r for r < R, case B, option on (Nx=33) | 33 | `4.1855e-02` | `0.0000e+00` | `<= 5.90e-02` | dimensionless | **PASS** |
| `B` | S33 L2 relative error against tau_33 = 2 eta u_r / r for r < R, case B, option on (Nx=65) | 65 | `1.0839e-02` | `0.0000e+00` | `<= 1.52e-02` | dimensionless | **PASS** |
| `B` | S33 L2 relative error against tau_33 = 2 eta u_r / r for r < R, case B, option on (Nx=129) | 129 | `2.7375e-03` | `0.0000e+00` | `<= 3.90e-03` | dimensionless | **PASS** |
| `A` | abs(Omega_fit(0) / Omega(0) - 1), fit of degree 5 in rho^2 for 0.18 <= rho <= 0.85, case A, option on (Nx=33) | 33 | `3.5554e-02` | `0.0000e+00` | `<= 5.00e-02` | dimensionless | **PASS** |
| `A` | abs(Omega_fit(0) / Omega(0) - 1), fit of degree 5 in rho^2 for 0.18 <= rho <= 0.85, case A, option on (Nx=65) | 65 | `8.5180e-03` | `0.0000e+00` | `<= 1.20e-02` | dimensionless | **PASS** |
| `A` | abs(Omega_fit(0) / Omega(0) - 1), fit of degree 5 in rho^2 for 0.18 <= rho <= 0.85, case A, option on (Nx=129) | 129 | `2.1618e-03` | `0.0000e+00` | `<= 3.10e-03` | dimensionless | **PASS** |
| `A` | abs(Omega_fit(0) / Omega(0) - 1), fit of degree 5 in rho^2 for 0.18 <= rho <= 0.85, case A, option off (Nx=33) | 33 | `3.3180e-02` | `0.0000e+00` | `<= 4.70e-02` | dimensionless | **PASS** |
| `A` | abs(Omega_fit(0) / Omega(0) - 1), fit of degree 5 in rho^2 for 0.18 <= rho <= 0.85, case A, option off (Nx=65) | 65 | `7.9723e-03` | `0.0000e+00` | `<= 1.12e-02` | dimensionless | **PASS** |
| `A` | abs(Omega_fit(0) / Omega(0) - 1), fit of degree 5 in rho^2 for 0.18 <= rho <= 0.85, case A, option off (Nx=129) | 129 | `1.9917e-03` | `0.0000e+00` | `<= 2.80e-03` | dimensionless | **PASS** |
| `A` | abs(Omega_node / Omega(0) - 1) at the centre basic node, case A, option on (Nx=33) | 33 | `9.7009e-03` | `0.0000e+00` | `<= 1.36e-02` | dimensionless | **PASS** |
| `A` | abs(Omega_node / Omega(0) - 1) at the centre basic node, case A, option on (Nx=65) | 65 | `7.0838e-03` | `0.0000e+00` | `<= 1.00e-02` | dimensionless | **PASS** |
| `A` | abs(Omega_node / Omega(0) - 1) at the centre basic node, case A, option on (Nx=129) | 129 | `4.2314e-03` | `0.0000e+00` | `<= 6.00e-03` | dimensionless | **PASS** |
| `A` | Convergence order log2(e(Nx=65) / e(Nx=129)) of the velocity L2 error, case A, option on | 129 | `2.0084e+00` | `2.0000e+00` | `atol=3.50e-01, rtol=0.00e+00` | dimensionless | **PASS** |
| `A` | Convergence order log2(e(Nx=65) / e(Nx=129)) of the velocity L2 error, case A, option off | 129 | `2.0129e+00` | `2.0000e+00` | `atol=3.50e-01, rtol=0.00e+00` | dimensionless | **PASS** |
| `B` | Convergence order log2(e(Nx=65) / e(Nx=129)) of the velocity L2 error, case B, option on | 129 | `2.0094e+00` | `2.0000e+00` | `atol=3.50e-01, rtol=0.00e+00` | dimensionless | **PASS** |
| `B` | Convergence order log2(e(Nx=65) / e(Nx=129)) of the velocity L2 error for r < R, case B, option on | 129 | `2.0095e+00` | `2.0000e+00` | `atol=3.50e-01, rtol=0.00e+00` | dimensionless | **PASS** |
| `B` | Convergence order log2(e(Nx=65) / e(Nx=129)) of the S33 L2 error, case B, option on | 129 | `1.9853e+00` | `2.0000e+00` | `atol=3.50e-01, rtol=0.00e+00` | dimensionless | **PASS** |
| `A` | Convergence order log2(e(Nx=65) / e(Nx=129)) of the fitted centre rotation error, case A, option on | 129 | `1.9783e+00` | `2.0000e+00` | `atol=3.50e-01, rtol=0.00e+00` | dimensionless | **PASS** |
| `A` | Convergence order log2(e(Nx=65) / e(Nx=129)) of the fitted centre rotation error, case A, option off | 129 | `2.0010e+00` | `2.0000e+00` | `atol=3.50e-01, rtol=0.00e+00` | dimensionless | **PASS** |
| `A` | Fitted centre rotation ratio, option on to off, case A (Nx=65) | 65 | `6.6703e-01` | `6.6667e-01` | `atol=0.00e+00, rtol=8.00e-04` | dimensionless | **PASS** |
| `A` | Fitted centre rotation ratio, option on to off, case A (Nx=129) | 129 | `6.6678e-01` | `6.6667e-01` | `atol=0.00e+00, rtol=2.50e-04` | dimensionless | **PASS** |
| `A` | Velocity L2 norm ratio over r < R, option on to off, case A (Nx=65) | 65 | `8.5235e-01` | `8.5280e-01` | `atol=0.00e+00, rtol=7.50e-04` | dimensionless | **PASS** |
| `A` | Velocity L2 norm ratio over r < R, option on to off, case A (Nx=129) | 129 | `8.5268e-01` | `8.5280e-01` | `atol=0.00e+00, rtol=2.10e-04` | dimensionless | **PASS** |
| `solver` | Largest relative difference of vx, vy, pr between the 6- and 4-variable paths, case A, option on (Nx=33) | 33 | `1.3014e-10` | `0.0000e+00` | `<= 1.00e-07` | dimensionless | **PASS** |
| `solver` | Largest relative difference of vx, vy, pr between the 6- and 4-variable paths, case B, option on (Nx=33) | 33 | `2.1994e-10` | `0.0000e+00` | `<= 1.00e-07` | dimensionless | **PASS** |

```@raw html
<!-- VALIDATION_RECORD_END: spherical_momentum -->
```

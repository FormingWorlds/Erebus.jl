# Verification and Validation

This page summarizes the benchmarks and tests used to ensure numerical accuracy and physical fidelity in `Erebus.jl`.

---

## 1. 1D Terzaghi Analytical Consolidation Benchmark

The primary verification test for the coupled Stokes-Darcy formulation is the one-dimensional Terzaghi consolidation problem.

- **Setup**: A saturated porous column subjected to instantaneous compressive surface load with a free-draining top surface ($P_f = 0$) and impermeable base.
- **Analytical Solution**: Fourier series representation of pore pressure dissipation over dimensionless time factor $T_v = c_v t / h^2$.
- **Result**: The numerical pore pressure field computed by the monolithic solver matches the analytical Fourier series solution across all column depths with a relative error $< 3.3\%$. See the [Terzaghi Consolidation Tutorial](../tutorials/terzaghi_consolidation.md) for the full derivation and code.

---

## 2. Poroelastic Constitutive Limits

The constitutive poroelastic functions in `src/physics.jl` are verified against physical asymptotic limits:

1. **Incompressible Solid Skeleton ($\beta_s \to 0$)**:
   $$\lim_{\beta_s \to 0} K_{\text{BW}} = 1, \quad \lim_{\beta_s \to 0} B = \frac{\beta_\phi}{\beta_\phi + \phi \beta_f}$$
   Verified across porosity values $\phi \in [\phi_{\text{min}}, \phi_{\text{max}}]$.

2. **Incompressible Pore Fluid ($\beta_f \to 0$)**:
   As fluid compressibility approaches zero, Skempton coefficient $B$ approaches its undrained upper bound:
   $$\lim_{\beta_f \to 0} B = \frac{\beta_d - \beta_s}{\beta_d - (1 - \phi)\beta_s} \approx 1$$

3. **Porosity Bounding Guarantees**:
   Constitutive routines clamp porosity to $[\phi_{\text{min}}, \phi_{\text{max}}]$ to prevent singular division when $\phi \to 0$ or $\phi \to 1$.

---

## 3. Matrix Block Consistency Tests

The individual coupling sub-blocks in `assemble_hydromechanical_lse!` are tested independently:
- **Solid continuity diagonal ($L[P_t, P_t]$)**: Matches $-(1 - K_{\text{BW}}) / (\Delta t / \beta_s)$.
- **Fluid continuity cross-term ($L[P_f, P_t]$)**: Matches $K_{\text{BW}} / \Delta t$.
- **Fluid continuity diagonal ($L[P_f, P_f]$)**: Matches $-\beta_{\text{fluid,eff}} / \Delta t$.

---

## 4. Radiogenic Energy and Mass Conservation

- **Radioactive Decay Kinetics**: Exponential decay of $^{26}\text{Al}$ and $^{60}\text{Fe}$ is verified against analytical integration $\int_0^t Q(t') dt'$, demonstrating exact conservation of energy release.
- **Dehydration Mass Balance**: Fluid released from decomposing hydrated silicates matches matrix mass loss ($DMP$).

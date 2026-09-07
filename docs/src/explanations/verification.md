# Verification and Benchmarking Strategy

This page explains the verification hierarchy, numerical consistency checks, and conservation audits in `Erebus.jl`.

For complete physical derivations, parameter tables, and benchmark figures, see the dedicated chapters in the [Validation section](../validation/index.md).

---

## The Verification Hierarchy

`Erebus.jl` couples thermo-mechanical solid deformation, Darcy fluid percolation, mineral phase transitions, and volatile loss. To establish simulation fidelity across this multi-physics chain, verification proceeds through four sequential tiers:

```
[ Tier 1: Analytical Closed-Form Solutions ]
                  │
[ Tier 2: Asymptotic Constitutive Limits   ]
                  │
[ Tier 3: Discrete Matrix Consistency      ]
                  │
[ Tier 4: Global Conservation Audits       ]
```

---

## 1. Analytical Closed-Form Solutions

Analytical benchmarks compare numerical outputs against exact mathematical solutions:

- **1D Terzaghi Consolidation**: Verifies coupled Stokes-Darcy dissipation against the classical Fourier series solution in a consolidating porous column. Pointwise relative errors remain below $3.5\%$. See the [Terzaghi Consolidation Tutorial](../tutorials/terzaghi_consolidation.md).
- **1D Stefan Moving-Boundary Front**: Verifies endothermic dehydration front propagation against the transcendental Neumann-Stefan similarity solution. The front location matches the analytical interface within grid resolution. See [Hydrothermal Reactions](../validation/hydrothermal_reactions.md).
- **Radionuclide Decay Kinetics**: Confirms analytic integration of $^{26}\text{Al}$ and $^{60}\text{Fe}$ heat release over multi-million-year timescales. See [Radionuclide Decay](../validation/radionuclides.md).
- **Jeans Kinetic Atmospheric Loss**: Verifies analytic time-integrated atmospheric mass loss and surface pressure evolution across light and heavy volatile species. See [Jeans Atmospheric Escape](../validation/jeans_escape.md).

---

## 2. Asymptotic Constitutive Limits

Constitutive parameterizations are tested against theoretical asymptotic bounds:

- **Incompressible Solid Skeleton ($\beta_s \to 0$)**: The Biot-Willis coefficient satisfies $\lim_{\beta_s \to 0} K_{\text{BW}} = 1$. The Skempton coefficient converges to $B = \frac{\beta_\phi}{\beta_\phi + \phi(1 - \phi)\beta_f}$.
- **Incompressible Pore Fluid ($\beta_f \to 0$)**: The Skempton coefficient converges to its undrained upper bound: $\lim_{\beta_f \to 0} B = 1$.
- **Porosity Safeguards**: Constitutive routines clamp porosity to $[\phi_{\text{min}}, \phi_{\text{max}}]$ to prevent numerical singularities during compaction and fluid dilation.
- **Solubility Thresholds**: Dissolved volatile concentrations vanish continuously at zero pressure ($w_{\text{sat}} \to 0$ as $P_f \to 0$). Dissolved species respect graphite and sulfide saturation ceilings.

For full derivations and asymptotic plots, see [Permeability & Hydrofracture](../validation/permeability.md) and [H-C-N-S Volatile Solubility](../validation/hcns_solubility.md).

---

## 3. Discrete Operator and Matrix Consistency

Discrete finite-difference operators are tested independently before assembling global systems:

- **Staggered-Grid Operator Symmetry**: Discrete divergence and gradient operators satisfy discrete adjoint properties on staggered cell faces and centers.
- **Stokes-Darcy Schur Coupling**: Coupling sub-blocks in `assemble_hydromechanical_lse!` satisfy cross-coupling consistency ($L[P_t, P_f] = L[P_f, P_t]$).
- **Nonlinear Picard Convergence**: Visco-elasto-plastic yielding iterations monitor Euclidean norm residuals until mechanical yielding errors fall below user-specified tolerances.

For discretization details and finite-difference stencils, see [Discretization and Numerics](discretization_numerics.md).

---

## 4. Global Conservation Audits

Simulations must satisfy exact physical conservation across all time steps:

- **Energy Balance**: Integrated radiogenic decay, latent heat absorption/release, and surface radiative loss balance total planetary internal energy changes.
- **Water Mass Conservation**: In closed systems, mineral lattice water plus pore fluid water remains strictly constant:
  $$\Delta m_{\text{mineral,water}} + \Delta m_{\text{pore,water}} = 0$$
- **Volatile Inventory Accounting**: Planetary volatile mass balances satisfy:
  $$M_{\text{initial}} = M_{\text{interior}}(t) + M_{\text{atm}}(t) + M_{\text{escaped}}(t)$$

---

## Validation Benchmark Directory

Detailed benchmarks, literature anchors, governing equations, and reproduction figures are organized into dedicated validation chapters:

| Benchmark Area | Focus | Reference Chapter |
|:---|:---|:---|
| **Thermal Conduction** | Porous thermal conductivity and spherical metric divergence | [Thermal Conduction](../validation/thermal_conduction.md) |
| **Pore Permeability** | Kozeny-Carman flow and overpressure hydrofracturing | [Permeability](../validation/permeability.md) |
| **Radionuclide Heating** | $^{26}\text{Al}$ and $^{60}\text{Fe}$ volumetric decay heating | [Radionuclides](../validation/radionuclides.md) |
| **Hydrothermal Fluid** | Arrhenius water viscosity and hydrothermal phase limits | [Fluid Viscosity](../validation/fluid_viscosity.md) |
| **Surface Radiation** | Stefan-Boltzmann cooling and protoplanetary disk dispersal | [Surface Radiation](../validation/disk_radiation.md) |
| **Metamorphic Reactions** | Serpentine hydration-dehydration kinetics and latent heat | [Hydrothermal Reactions](../validation/hydrothermal_reactions.md) |
| **Magma Convection** | Silicate melting, suspension rheology, and soft turbulence | [Silicate Melting](../validation/rock_melting.md) |
| **Surface Drainage** | Leaky Robin venting and Clausius-Clapeyron cold trap | [Cold Surface Venting](../validation/cold_surface_venting.md) |
| **Cryogenic Sealing** | Pore ice clogging and episodic lid hydrofracture | [Hydrofracture Venting](../validation/hydrofracture_venting.md) |
| **Volatile Speciation** | Multi-species H-C-N-S solubility, graphite, and SCSS limits | [H-C-N-S Solubility](../validation/hcns_solubility.md) |
| **Atmospheric Escape** | Maxwellian effusion, scale height, and pressure feedback | [Jeans Escape](../validation/jeans_escape.md) |

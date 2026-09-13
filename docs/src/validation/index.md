# Validation and Physical Benchmarks

The physical formulations in `Erebus.jl` are anchored against peer-reviewed literature, analytical limits, and laboratory measurements. This section documents the verification hierarchy, numerical consistency checks, conservation audits, and physical benchmarks governing simulation fidelity.

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

- **1D Terzaghi Consolidation**: Verifies coupled Stokes-Darcy dissipation against the classical Fourier series solution in a consolidating porous column. Pointwise relative errors remain below $3.5\%$. See [1D Terzaghi Consolidation Benchmark](terzaghi_consolidation.md).
- **1D Stefan Moving-Boundary Front**: Verifies endothermic dehydration front propagation against the transcendental Neumann-Stefan similarity solution. The front location matches the analytical interface within grid resolution. See [Hydrothermal Reactions](hydrothermal_reactions.md).
- **Radionuclide Decay Kinetics**: Confirms analytic integration of $^{26}\text{Al}$ and $^{60}\text{Fe}$ heat release over multi-million-year timescales. See [Radionuclide Decay](radionuclides.md).
- **Jeans Kinetic Atmospheric Loss**: Verifies analytic time-integrated atmospheric mass loss and surface pressure evolution across light and heavy volatile species. See [Jeans Atmospheric Escape](jeans_escape.md).

---

## 2. Asymptotic Constitutive Limits

Constitutive parameterizations are tested against theoretical asymptotic bounds:

- **Incompressible Solid Skeleton ($\beta_s \to 0$)**: The Biot-Willis coefficient satisfies $\lim_{\beta_s \to 0} K_{\text{BW}} = 1$. The Skempton coefficient converges to $B = \frac{\beta_\phi}{\beta_\phi + \phi(1 - \phi)\beta_f}$.
- **Incompressible Pore Fluid ($\beta_f \to 0$)**: The Skempton coefficient converges to its undrained upper bound: $\lim_{\beta_f \to 0} B = 1$.
- **Porosity Safeguards**: Constitutive routines clamp porosity to $[\phi_{\text{min}}, \phi_{\text{max}}]$ to prevent numerical singularities during compaction and fluid dilation.
- **Solubility Thresholds**: Dissolved volatile concentrations vanish continuously at zero pressure ($w_{\text{sat}} \to 0$ as $P_f \to 0$). Dissolved species respect graphite and sulfide saturation ceilings.

For full derivations and asymptotic plots, see [Permeability & Hydrofracture](permeability.md) and [H-C-N-S Volatile Solubility](hcns_solubility.md).

---

## 3. Discrete Operator and Matrix Consistency

Discrete finite-difference operators are tested independently before assembling global systems:

- **Staggered-Grid Operator Symmetry**: Discrete divergence and gradient operators satisfy discrete adjoint properties on staggered cell faces and centers.
- **Stokes-Darcy Schur Coupling**: Coupling sub-blocks in `assemble_hydromechanical_lse!` satisfy cross-coupling consistency ($L[P_t, P_f] = L[P_f, P_t]$).
- **Nonlinear Picard Convergence**: Visco-elasto-plastic yielding iterations monitor Euclidean norm residuals until mechanical yielding errors fall below user-specified tolerances.

For discretization details and finite-difference stencils, see [Discretization & Numerics](../explanations/discretization_numerics.md).

---

## 4. Global Conservation Audits

Simulations must satisfy exact physical conservation across all time steps:

- **Energy Balance**: Integrated radiogenic decay, latent heat absorption/release, and surface radiative loss balance total planetary internal energy changes.
- **Water Mass Conservation**: In closed systems, mineral lattice water plus pore fluid water remains strictly constant:
  $$\Delta m_{\text{mineral,water}} + \Delta m_{\text{pore,water}} = 0$$
- **Volatile Inventory Accounting**: Planetary volatile mass balances satisfy:
  $$M_{\text{initial}} = M_{\text{interior}}(t) + M_{\text{atm}}(t) + M_{\text{escaped}}(t)$$

---

## Validation Matrix

| Physical Process | Governing Theory | Primary Reference | Verification Test Suite |
|:---|:---|:---|:---|
| **1D Terzaghi Consolidation** | Poroelastic excess pore pressure consolidation in porous column | Terzaghi (1925); Wang (2000) | `test/test_numerics.jl` |
| **Thermal Conduction & Geometry** | Two-phase porous conductivity and 3D spherical metric divergence | Gerya (2019); Hubmann (2022) | `test/test_geometry_radiation.jl`, `test/test_physics.jl` |
| **Permeability & Hydrofracturing** | Kozeny-Carman flow and Terzaghi effective overpressure failure | Carman (1937); Terzaghi (1925); Wang (2000) | `test/test_physics.jl`, `test/test_numerics.jl` |
| **Radionuclide Decay** | Short-lived radioactive heating ($^{26}\text{Al}$, $^{60}\text{Fe}$) | Russell et al. (1996); Tachibana & Huss (2003); Lichtenberg et al. (2019) | `test/test_physics.jl` |
| **Fluid Viscosity & Phase Transitions** | Arrhenius water viscosity and hydrothermal phase limits | Hubmann (2022); Gerya (2019) | `test/test_physics.jl` |
| **Surface Radiation & Disk Evolution** | Stefan-Boltzmann boundary and protoplanetary disk clearing | Chiang & Goldreich (1997); Drążkowska & Dullemond (2018); Williams et al. (2026) | `test/test_geometry_radiation.jl` |
| **Hydrothermal Reactions** | Hydration and dehydration kinetics (Arrhenius) with latent heat and mass coupling | Hubmann (2022); Gerya (2019) | `test/test_reaction_pathways.jl` |
| **Silicate Rock Melting & Magma Convection** | Linear melt fraction, apparent heat capacity, melt-weakened rheology, and sub-grid soft turbulence | Gerya (2019); Costa et al. (2009); Solomatov (2007) | `test/test_melting.jl`, `test/test_soft_turbulence.jl` |
| **Cold Surface Venting & Ice Sealing** | Darcy Robin leaky drainage, Clausius-Clapeyron cold trap, and cryogenic pore ice sealing | Hubmann (2022); Gerya (2019) | `test/test_venting_thermodynamics.jl`, `test/test_venting_darcy_sink.jl`, `test/test_venting_integration.jl`, `test/test_hydrofracture_venting.jl` |
| **H-C-N-S Volatile Solubility & Speciation** | Multi-species (H, C, N, S) solubility laws, graphite and sulfide saturation limits, organic devolatilization | Burnham (1979); Dixon et al. (1995); Dasgupta et al. (2022); Boulliung & Wood (2022) | `test/test_volatile_solubility.jl`, `test/test_volatile_solubility_hcns.jl` |
| **Volatile Retention Floors & Vent Drainage** | Nominally anhydrous mineral retention floors, vacuum exsolution clamping, and low-temperature venting drainage coupling | Hirschmann et al. (2006); Peslier et al. (2017); Shcheka et al. (2006); Hirschmann (2018); Li et al. (2013) | `test/test_volatile_retention.jl` |
| **Redox Buffers & Electron Accounting** | Solid-oxide linear buffers (Frost 1991), graphite CCO inversion, and extensive electron budget conservation | Frost (1991); Campbell et al. (2009); Evans (2006, 2012) | `test/test_redox.jl` |
| **Jeans Kinetic Atmospheric Escape** | Maxwellian effusion, scale height relaxation, mass conservation, and surface pressure feedback | Jeans (1925); Catling & Kasting (2017) | `test/test_jeans_escape.jl` |
| **Iron Core Formation** | Fe-FeS percolation Darcy flow, Stokes droplet settling, regime handover, dissipation heating, and mass conservation | Yoshino et al. (2003); Rubie et al. (2003, 2015); Monteux et al. (2009); Lichtenberg et al. (2019) | `test/test_core_formation.jl` |
| **Core Geochemistry & Volatiles** | Siderophile partitioning (H, C, N, S), sulfur suppression of carbon, and dynamic alloy density | Grewal et al. (2019a, 2019b); Clesi et al. (2018); Boujibar et al. (2014) | `test/test_core_volatile_partitioning.jl` |
| **Normative Accessory Minerals** | Stoichiometric sub-eutectic mineral exsolution, thermal eutectic dissolution, and meteorite classification | Benedix et al. (2000); Goldstein et al. (2009); Chabot & Drake (1999) | `test/test_normative_accessory_minerals.jl` |
| **Hydrothermal Subgrid Convection** | Porous Rayleigh-Darcy scaling, boundary-layer free-fluid Rayleigh scaling, smoothstep porosity transition, and cell-Péclet damping | Horton & Rogers (1945); Lapwood (1948); Kraichnan (1962); Elder (1967) | `test/test_hydrothermal_convection.jl` |
| **Planetesimal Accretion & Impact Heating** | Bondi and Hill pebble accretion, Safronov gravitational focusing, exact 3D spherical mapping, impact heating, and radiogenic clock inheritance | Safronov (1972); Ormel & Klahr (2010); Lambrechts & Johansen (2012); Lichtenberg et al. (2019) | `test/test_accretion.jl`, `test/test_config.jl` |
| **Telescoping Domain** | Constant cell spacing coordinate doubling, marker distance preservation, sticky-air buffer generation, and physical conservation | Gerya (2019); Crameri et al. (2012) | `test/test_telescoping.jl` |
| **Volatiles & Refractory Mixtures** | Multi-snowline condensation, ammonia eutectic freezing depression, and refractory organic pyrolysis | Bergin et al. (2026); Kama et al. (2019) | `test/test_volatile_mixtures.jl` |
| **Multi-Stage Accretion Sequence** | Planetesimal collisions before settling pebble accretion, pebble isolation mass, and smoothstep transitions | Visser & Ormel (2016); Liu et al. (2019); Lambrechts et al. (2014) | `test/test_multistage_accretion.jl` |

---

## Standards of Verification

Each validation module fulfills three requirements:
1. **Literature Grounding**: Every parameterization maps to a verified scientific publication with a valid DOI.
2. **Analytical Limits**: Numerical implementations are tested against closed-form mathematical limits (e.g. pure solid limit, infinite-time decay, steady-state conduction).
3. **Physical Invariants**: Unit tests assert physical conservation (energy, mass), boundedness ($T > 0\text{ K}$, $\phi \in [0, 1]$), and monotonicity.

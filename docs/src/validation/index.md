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

## Figure Provenance Taxonomy

All validation and benchmark figures in the `Erebus.jl` documentation are classified into four provenance tiers:

- **Class A (2D Simulation Output)**: Figures generated directly from multi-dimensional `Erebus.jl` simulation loop solves (`simulation_loop` or output checkpoint files).
- **Class B (Julia Library Exporter / 1D Benchmark Solver)**: Figures generated from Julia benchmark scripts that execute compiled `Erebus.jl` modules or 1D finite-difference benchmark solvers exporting structured data.
- **Class C (Analytical / Empirical Reference Formulation)**: Figures illustrating closed-form analytical solutions, theoretical limits, or literature empirical parameterizations evaluated in Python for comparison with discrete simulation behavior.
- **Class D (Schematic Diagram)**: Procedural vector graphics illustrating computational pipelines, reservoir architectures, or coordinate systems.

Each validation page provides a standardized summary block detailing diagnostic targets, comparison standards, figure provenance, and automated test coverage.

---

## 1. Analytical Closed-Form Solutions

Analytical benchmarks compare numerical outputs against exact mathematical solutions:

- **1D Terzaghi Consolidation**: Verifies coupled Stokes-Darcy dissipation against the classical Fourier series solution in a consolidating porous column. Pointwise relative errors remain below $3.5\%$. See [1D Terzaghi Consolidation Benchmark](terzaghi_consolidation.md).
- **1D Stefan Moving-Boundary Front**: Verifies endothermic dehydration front propagation against the transcendental Neumann-Stefan similarity solution. The front location matches the analytical interface to grid resolution. See [Hydrothermal Reactions](hydrothermal_reactions.md).
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

| Physical Process | Governing Theory | Primary Reference | Provenance | Verification Test Suite | Quantitative Tolerance |
|:---|:---|:---|:---:|:---|:---|
| **1D Terzaghi Consolidation** | Poroelastic excess pore pressure consolidation in porous column | Terzaghi (1925); Wang (2000) | Class B | `test/test_numerics.jl` | Pointwise relative error $< 3.5\%$ throughout column |
| **Thermal Conduction & Geometry** | Two-phase porous conductivity and 3D spherical metric divergence | Gerya (2019); Hubmann (2022); Carslaw & Jaeger (1959) | Class C / B | `test/test_geometry_radiation.jl`, `test/test_physics.jl` | 2D Fourier temperature diffusion $L_2 < 1.0\times 10^{-3}$; energy drift $< 10^{-12}$ |
| **Permeability & Hydrofracturing** | Kozeny-Carman flow and Terzaghi effective overpressure failure | Carman (1937); Terzaghi (1925); Wang (2000) | Class C / B | `test/test_physics.jl`, `test/test_hydrofracture_stability.jl` | Biot-Willis $K_{\text{BW}}$ and Skempton $B$ match $< 10^{-12}$; ramp match $< 10^{-6}$ |
| **Radionuclide Decay** | Short-lived radioactive heating ($^{26}\text{Al}$, $^{60}\text{Fe}$) | Russell et al. (1996); Tachibana & Huss (2003); Lichtenberg et al. (2019) | Class B | `test/test_physics.jl` | Radioactive power decay matches exponential law $< 10^{-10}$ |
| **Fluid Viscosity & Phase Transitions** | Arrhenius water viscosity and hydrothermal phase limits | Hubmann (2022); Gerya (2019) | Class C | `test/test_physics.jl` | Arrhenius viscosity matches formula $< 10^{-12}$ |
| **Surface Radiation & Disk Evolution** | Stefan-Boltzmann boundary and protoplanetary disk clearing | Chiang & Goldreich (1997); Drążkowska & Dullemond (2018) | Class C | `test/test_geometry_radiation.jl` | $T_{\text{disk}}(r, t)$ match analytical formula $< 10^{-12}$ |
| **Hydrothermal Reactions** | Hydration and dehydration kinetics (Arrhenius) with latent heat and mass coupling | Hubmann (2022); Gerya (2019); Neumann-Stefan solution | Class A | `test/test_reaction_pathways.jl`, `test/test_stefan_benchmark.jl` | Reaction front position resolves to grid cell size $\Delta y$; mass conservation $\Delta m_{\text{solid}} + \Delta m_{\text{pore}} = 0$ closed to $10^{-12}$ |
| **Silicate Rock Melting & Magma Convection** | Linear melt fraction, apparent heat capacity, melt-weakened rheology, and sub-grid soft turbulence | Gerya (2019); Costa et al. (2009); Solomatov (2007) | Class B (1D FD) | `test/test_melting.jl`, `test/test_soft_turbulence.jl` | Apparent heat capacity integral $\int c_p^{\text{eff}} dT = c_p \Delta T + L_m$ closed to $10^{-12}$; $F_m \in [0, 1]$ exact |
| **Cold Surface Venting & Ice Sealing** | Darcy Robin leaky drainage, Clausius-Clapeyron cold trap, and cryogenic pore ice sealing | Hubmann (2022); Washburn (1924) | Class C | `test/test_venting_thermodynamics.jl`, `test/test_venting_darcy_sink.jl`, `test/test_venting_integration.jl` | Exponential decay match $< 10^{-5}$; $P_{\text{surf}} \ge P_{\text{amb}}$ |
| **H-C-N-S Volatile Solubility & Speciation** | Multi-species (H, C, N, S) solubility laws, graphite and sulfide saturation limits | Burnham (1979); Dixon et al. (1995); Dasgupta et al. (2022) | Class C | `test/test_volatile_solubility.jl`, `test/test_volatile_solubility_hcns.jl` | Exsolved volatile mass matches analytical saturation limits $< 10^{-6}$ |
| **Volatile Retention Floors & Vent Drainage** | Nominally anhydrous mineral retention floors, vacuum exsolution clamping, and venting drainage | Hirschmann et al. (2006); Peslier et al. (2017); Shcheka et al. (2006) | Class C | `test/test_volatile_retention.jl` | Retained concentrations strictly respect floor $w_{\text{ret}} \ge w_{\text{floor}}$ |
| **Redox Buffers & Electron Accounting** | Solid-oxide linear buffers (Frost 1991), graphite CCO inversion, and extensive electron budget conservation | Frost (1991); Campbell et al. (2009); Evans (2012) | Class B | `test/test_redox.jl`, `test/test_redox_metal.jl` | Electron budget conservation $\Delta n_{e^-} = 0$ closed to $10^{-12}$; buffer $\log_{10} f_{\text{O2}}$ match $< 10^{-4}$ |
| **Jeans Kinetic Atmospheric Escape** | Maxwellian effusion, scale height relaxation, mass conservation, and surface pressure feedback | Jeans (1925); Catling & Kasting (2017) | Class C / B | `test/test_jeans_escape.jl` | Effusion flux matches kinetic rate $< 10^{-6}$; exobase scaling exact $< 10^{-12}$ |
| **Iron Core Formation** | Fe-FeS percolation Darcy flow, Stokes droplet settling, regime handover, dissipation heating | Rubie et al. (2015); Monteux et al. (2009); Lichtenberg et al. (2019) | Class B (1D FD) | `test/test_core_formation.jl`, `test/test_reference_runs.jl` | Analytical 1D gravity $L_2 < 10^{-4}$; core metal conservation $< 10^{-12}$ |
| **Core Geochemistry & Volatiles** | Siderophile partitioning (H, C, N, S), sulfur suppression of carbon, and dynamic alloy density | Grewal et al. (2019a, 2019b); Clesi et al. (2018); Boujibar et al. (2014) | Class C | `test/test_core_volatile_partitioning.jl`, `test/test_redox_metal.jl` | Partition mass balance closed to $10^{-12}$; $D_i$ match equations $< 10^{-6}$ |
| **Normative Accessory Minerals** | Stoichiometric sub-eutectic mineral exsolution, thermal eutectic dissolution, and meteorite classification | Benedix et al. (2000); Goldstein et al. (2009); Chabot & Drake (1999) | Class C | `test/test_normative_accessory_minerals.jl` | Stoichiometric mass balance $\sum m_{\text{minerals}} = m_{\text{parcel}}$ closed to $10^{-12}$ |
| **Hydrothermal Subgrid Convection** | Porous Rayleigh-Darcy scaling, boundary-layer free-fluid Rayleigh scaling, smoothstep porosity transition | Horton & Rogers (1945); Lapwood (1948); Elder (1967) | Class C | `test/test_hydrothermal_convection.jl` | $Ra_{\text{crit}} = 4\pi^2$ to $0.5\%$; $Nu \sim Ra/Ra_{\text{crit}}$ to $1.0\%$ |
| **Planetesimal Accretion & Impact Heating** | Bondi and Hill pebble accretion, Safronov gravitational focusing, exact 3D spherical mapping | Safronov (1972); Ormel & Klahr (2010); Lambrechts & Johansen (2012) | Class C | `test/test_accretion.jl` | Pebble capture cross-section matches 2D/3D Hill and Bondi limits $< 10^{-6}$ |
| **Telescoping Domain** | Constant cell spacing coordinate doubling, marker distance preservation, sticky-air buffer generation | Gerya (2019); Crameri et al. (2012) | Class C | `test/test_telescoping.jl` | Marker coordinates relative to planet center invariant under domain doubling $< 10^{-14}\text{ m}$ |
| **Volatiles & Refractory Mixtures** | Multi-snowline condensation, ammonia eutectic freezing depression, and refractory organic pyrolysis | Bergin et al. (2026); Kama et al. (2019) | Class C | `test/test_volatile_mixtures.jl` | Condensation temperatures and eutectic freezing curves match analytical models $< 10^{-6}$ |
| **Multi-Stage Accretion Sequence** | Planetesimal collisions before settling pebble accretion, pebble isolation mass, and smoothstep transitions | Visser & Ormel (2016); Liu et al. (2019); Lambrechts et al. (2014) | Class C | `test/test_multistage_accretion.jl`, `test/test_accretion.jl` | Mass growth rate $\dot{M}$ matches analytical focusing equations $< 10^{-6}$ |
| **Coupled 1D Atmosphere & Degassing** | Semi-grey radiative equilibrium, greenhouse blanketing, and multi-species hydrodynamic escape | Guillot (2010); Dixon et al. (1995); Gaillard & Scaillet (2014); Attia & Lichtenberg (2026) | Class D / C / B | `test/test_atmosphere.jl`, `test/test_volatile_solubility.jl` | Radiative profile $L_\infty < 10^{-4}$; CHNOS speciation equilibrium $< 10^{-5}$ |
| **Dehydration-Darcy Coupling** | Poroelastic overpressure coupling, 4-variable condensed assembly, and unified surface venting | McKenzie (1984); Young et al. (1999); Hubmann (2022) | Class B | `test/test_dehydration_darcy_coupling.jl`, `test/test_atmosphere.jl` | Fluid pressure equivalence $\|P_f^{(4)} - P_f^{(6)}\|_\infty / \|P_f^{(6)}\|_\infty < 10^{-8}$ |
| **Magma Compaction & Sills** | Matrix compaction, dynamic compaction pressure, and 1D sill conductive cooling | McKenzie (1984); Scott & Stevenson (1984); Jaeger (1957) | Class B | `test/test_magma_transport.jl`, `test/test_sill_cooling.jl` | Dynamic compaction pressure error $< 10^{-7}$; sill cooling $L_2 < 1.0\times 10^{-3}$ |

---

## Standards of Verification

Each validation module fulfills three requirements:
1. **Literature Grounding**: Every parameterization maps to a verified scientific publication with a valid citation.
2. **Analytical Limits**: Numerical implementations are tested against closed-form mathematical limits (e.g. pure solid limit, infinite-time decay, steady-state conduction).
3. **Physical Invariants**: Unit tests assert physical conservation (energy, mass), boundedness ($T > 0\text{ K}$, $\phi \in [0, 1]$), and monotonicity.

---

## Validation and Provenance Summary

| Attribute | Specification |
|:---|:---|
| **Target Physics / Diagnostic** | Multi-physics validation suite spanning analytical solutions, constitutive limits, discrete operators, and global conservation audits |
| **Reference Standard** | Terzaghi (1925); Carslaw & Jaeger (1959); Gerya (2019); Hubmann (2022); Lichtenberg et al. (2019, 2021) |
| **Figure Provenance** | Suite-wide index covering Class A, Class B, Class C, and Class D figures |
| **Generating Script** | `benchmarks/run_all.sh` |
| **Automated Verification Test** | Full test suite in `test/` |
| **Quantitative Tolerance** | Pointwise analytical matches $< 3.5\%$; conservation closures closed to $10^{-12}$; zero drift in invariants |

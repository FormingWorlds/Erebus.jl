# Dehydration-Darcy Fluid Overpressure Coupling and Unified Surface Venting Validation

This module validates the two-way coupling between prograde metamorphic dehydration fluid generation, the Stokes-Darcy fluid continuity equation, unified non-double-counted surface water venting, and dynamic equilibrium C-H-N-O-S gas speciation under variable oxygen fugacity in `Erebus.jl`.

---

## Theoretical Formulation

The physical theory of hydrothermal reactions, two-phase Stokes-Darcy porous flow, surface venting, and volatile speciation is detailed in [Volatile Degassing, Cold Surface Venting, and Atmospheric Escape](../explanations/degassing_and_venting.md).

Key constitutive relations validated on this page include:

- **Stokes-Darcy Fluid Continuity with Dehydration Source Term ($\Delta Q^f$):**
  $$\nabla \cdot \mathbf{v}_D = \Delta Q^f - \frac{1 - \phi}{\rho_s} \frac{D\rho_s}{Dt} - \frac{\phi}{\rho_f} \frac{D\rho_f}{Dt}$$
  where $\mathbf{v}_D = -\frac{k}{\eta_f} \nabla P_f$ is the Darcy filtration velocity and $\Delta Q^f = \text{DQPF}$ [$\text{s}^{-1}$] is the volumetric fluid production rate from mineral dehydration.

- **Condensed Hydromechanical Assembly (4-Variable and 6-Variable):**
  In both the classical 6-variable formulation ($v_x, v_y, P_t, v_{x}^D, v_{y}^D, P_f$) and the condensed 4-variable formulation ($v_x, v_y, P_t, P_f$), dehydration fluid production injects into the fluid continuity residual at internal pressure cell centers:
  $$R_{\text{fluid}}(i, j) = \text{DQPF}(i, j) - \beta_d K_{\text{bw}} \frac{P_t^n(i, j) - \frac{1}{K_{\text{sk}}} P_f^n(i, j)}{\Delta t}$$

- **Unified Surface Venting Water Mass Budget:**
  Surface water venting couples both Darcy pore fluid drainage ($m_{\text{pore}}$) and mobile mineral volatile drainage ($m_{\text{mineral}}$) additively, preserving distinct physical reservoirs on markers:
  $$m_{\text{H}_2\text{O}} = m_{\text{pore}} + m_{\text{mineral}}$$
  where $m_{\text{pore}} = \Delta M_{\text{vent}} \cdot L_{\text{3D}}$ when Darcy surface venting is active, and $m_{\text{mineral}} = M_{\text{vent}}^{\text{H}_2\text{O}} \cdot L_{\text{3D}}$ when retention drainage is active.

- **Dynamic Gas Speciation at the Surface:**
  Vented volatile elemental masses ($m_{\text{H}_2\text{O}}, m_{\text{C}}, m_{\text{N}}, m_{\text{S}}$) are partitioned into chemical equilibrium gas species ($\text{H}_2, \text{H}_2\text{O}, \text{CO}, \text{CO}_2, \text{CH}_4, \text{N}_2, \text{NH}_3, \text{H}_2\text{S}, \text{S}_2, \text{SO}_2$) evaluated at local surface conditions:
  $$\mathbf{m}_{\text{spec}} = \text{speciate\_vented\_volatiles}(m_{\text{H}_2\text{O}}, m_{\text{C}}, m_{\text{N}}, m_{\text{S}}, P_{\text{surf}}, T_{\text{surf}}, \Delta\text{IW})$$

---

## Invariants and Physical Limits

1. **Fluid Overpressure Positivity:** In a closed or low-permeability domain, positive dehydration volume production ($\text{DQPF} > 0$) strictly generates positive fluid overpressure ($\Delta P_f > 0$) relative to the background boundary pressure.
2. **Formulation Equivalence:** The condensed 4-variable formulation and the 6-variable formulation produce identical fluid pressure fields under dehydration forcing:
   $$\frac{\|P_f^{(4)} - P_f^{(6)}\|_\infty}{\|P_f^{(6)}\|_\infty} < 10^{-6}$$
3. **Outward Darcy Flux Divergence:** The divergence of outward Darcy filtration balances the fluid production rate and transient storage:
   $$\int_{\Omega} \nabla \cdot \mathbf{q}_D \, d\Omega = \int_{\Omega} \text{DQPF} \, d\Omega - \int_{\Omega} S_{\text{storage}} \, d\Omega$$
   approaching exact equivalence in the high-permeability quasi-steady limit.
4. **Water and Elemental Mass Conservation:** When dynamic speciation is inactive (`speciation_active = false`), the cumulative atmospheric water inventory matches the sum of cumulative vented pore fluid mass and drained mineral water mass:
   $$M_{\text{atm}}(\text{H}_2\text{O}) = M_{\text{vent}}^{\text{total}} + M_{\text{vent, mineral}}^{\text{total}}$$
   When dynamic speciation is active without graphite precipitation (`speciation_active = true`, `graphite_saturation = false`), elemental moles of hydrogen ($n_{\text{H}}$), carbon ($n_{\text{C}}$), nitrogen ($n_{\text{N}}$), and sulfur ($n_{\text{S}}$) are strictly conserved across equilibrium gas species:
   $$\sum_{\text{species}} \nu_{i, \text{sp}} n_{\text{sp}}^{\text{gas}} = n_i^{\text{vent}}$$
   When graphite saturation occurs, carbon excess precipitates as solid graphite while mobile H, N, and S partition into gas species.
5. **Redox-Dependent Gas Partitioning:** Under reducing conditions ($\Delta\text{IW} < 0$), dynamic speciation generates finite molecular hydrogen ($\text{H}_2 > 0$) and carbon monoxide ($\text{CO} > 0$), providing the light driving species necessary for hydrodynamic and crossover escape.

---

## Literature Anchors

- **McKenzie, D. (1984)**. The generation and compaction of partially molten rock. *Journal of Petrology*, 25(3), 713-754.
- **Young, E. D., Ash, R. D., England, P., & Rumble, D. (1999)**. Fluid flow in carbonaceous chondrite parent bodies and the origin of magnetites. *Science*, 286(5443), 1331-1335.
- **Fu, R. R., & Elkins-Tanton, L. T. (2014)**. The early thermal evolution of planetesimals: Implications for differentiated asteroids and carbonaceous chondrite parent bodies. *Earth and Planetary Science Letters*, 390, 128-137.
- **French, B. M. (1966)**. Some geological implications of equilibrium between graphite and a C-H-O gas at high temperatures and pressures. *Reviews of Geophysics*, 4(2), 223-253.
- **Holloway, J. R., Pan, V., & Gudmundsson, G. (1992)**. High-pressure fluid-absent melting in mantle carbonatite systems. *European Journal of Mineralogy*, 4(1), 105-114.
- **Gaillard, F., Bouhifd, M. A., Furi, E., Malavergne, V., & Marrocchi, Y. (2022)**. The diverse paths to planetary atmospheres: A volatile perspective. *Space Science Reviews*, 218(8), 65.

---

## Verification Test Suite
 
- `test/test_dehydration_darcy_coupling.jl`:
  - `@testset "ReactionConfig Schema & Serialization"`
  - `@testset "Analytical Dehydration Overpressure Pulse (Zero-Permeability Limit)"`
  - `@testset "4-Variable Condensed vs 6-Variable Solution Equivalence Under DQPF"`
  - `@testset "Outward Darcy Filtration Driven by Dehydration Overpressure"`
  - `@testset "Unified Surface Venting Water Mass Balance & Speciation"`
  - `@testset "Dynamic Equilibrium Gas Speciation Under Reducing Surface Venting"`

- `test/test_atmosphere.jl`:
  - `@testset "Simulation Loop Integration with Coupled Atmosphere"` verifies additive atmospheric H2O budget.

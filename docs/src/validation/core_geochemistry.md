# Metal-Silicate Volatile Partitioning and Core Geochemistry

This page documents the physical formulation, thermodynamic parameterizations, and numerical validation of metal-silicate volatile partitioning during iron core formation in `Erebus.jl`. The model couples empirical partition coefficients $D_i^{\text{met/sil}}$ for H, C, N, and S with local phase equilibration, dynamic sulfur density feedback, and conservative advective transport.

---

## 1. Physical Motivation

Volatile elements (hydrogen, carbon, nitrogen, and sulfur) play a decisive role in determining the planetary atmospheric budget, internal oxidation state, and long-term habitability of differentiated rocky worlds. During early planetesimal differentiation, metallic iron-nickel-sulfur liquids separate from coexisting silicate melts. Siderophile ("metal-loving") and chalcophile ("sulfur-loving") volatiles preferentially dissolve into the descending metallic liquid, sequestering volatile inventories into the planetesimal core.

Magmatic iron meteorites (such as groups IIAB, IIIAB, IVA, and IVB) represent the fragmented metallic cores of differentiated planetesimals accreted in the first few million years of Solar System history. Geochemical analyses of these meteorites demonstrate distinct volatile signatures:

1. **Group IIAB**: Formed in sulfur-rich parent bodies ($w_S \approx 10 - 17\text{ wt}\%$). They exhibit relatively low carbon concentrations ($100 - 350\text{ ppmw}$) despite chondritic precursors, producing sub-chondritic $(C/N)_{\text{metal}}$ ratios.
2. **Group IIIAB**: Formed at moderate sulfur contents ($w_S \approx 5 - 10\text{ wt}\%$), displaying intermediate carbon ($250 - 550\text{ ppmw}$) and nitrogen ($10 - 25\text{ ppmw}$) contents.
3. **Group IVA**: Highly depleted in sulfur ($w_S \approx 1 - 2.5\text{ wt}\%$) and moderately depleted in other volatiles ($C \approx 30 - 150\text{ ppmw}$, $N \approx 2 - 8\text{ ppmw}$).
4. **Group IVB**: Extremely refractory and volatile-poor ($w_S < 0.5\text{ wt}\%$, $C < 60\text{ ppmw}$, $N < 3\text{ ppmw}$); these concentrations indicate high-temperature condensation or early catastrophic degassing.

A realistic geochemical model of core formation must reproduce these geochemical patterns, specifically capturing how dissolved sulfur in metallic liquid alters the chemical activity and partitioning behavior of light elements.

---

## 2. Theoretical Formulation

The thermodynamics of metal-silicate volatile partitioning ($D_i^{\text{met/sil}}$ parameterizations for H, C, N, and S), phase equilibration, normative mineral crystallization, and conservative volatile drift-flux transport are derived in detail in [Iron Core Formation and Metal Segregation](../explanations/core_formation.md).

Key constitutive relations validated on this page include:

- **Nernst Partition Coefficient ($D_i^{\text{met/sil}}$):**
  $$D_i^{\text{met/sil}} = \frac{C_i^{\text{metal}}}{C_i^{\text{silicate}}}$$
- **Sulfur-Dependent Carbon and Nitrogen Partitioning (Grewal et al., 2019a, 2019b):**
  $$\log_{10} D_C = 1.80 + \frac{2200}{T} - 1.5 \times 10^{-8} \frac{P}{T} - 0.25 \, \Delta\text{IW} + 4.2 \ln(1 - X_S)$$
  $$\log_{10} D_N = 0.85 + \frac{1200}{T} - 0.25 \, \Delta\text{IW} + 0.60 \ln(1 - X_S)$$
- **Hydrogen (Clesi et al., 2018) and Sulfur (Boujibar et al., 2014) Partitioning:**
  $$\log_{10} D_H = -0.80 + \frac{300}{T} + 5.0 \times 10^{-8} \frac{P}{T} + 0.05 \, \Delta\text{IW}$$
  $$\log_{10} D_S = 2.80 - \frac{800}{T} + 1.0 \times 10^{-10} P - 0.20 \, \Delta\text{IW}$$
- **Phase Equilibration in the Silicate Melt Frame:**
  Metal-silicate volatile partitioning occurs between liquid metal and molten silicate:
  $$C_{i,\text{sil\_melt}} = \frac{C_{i,\text{sil\_bulk}}}{F_{\text{melt}}}$$
  For an unconstrained system with total volatile mass $M_{i,\text{tot}} = m_{\text{sil}} C_{i,\text{sil\_bulk}} + m_{\text{met}} C_{i,\text{met}}$, the thermodynamic equilibrium concentrations satisfy:
  $$C_{i,\text{sil\_melt}}^{\text{eq}} = \frac{M_{i,\text{tot}}}{m_{\text{sil}} F_{\text{melt}} + m_{\text{met}} D_i}, \quad C_{i,\text{met}}^{\text{eq}} = D_i \cdot C_{i,\text{sil\_melt}}^{\text{eq}}, \quad C_{i,\text{sil\_bulk}}^{\text{eq}} = F_{\text{melt}} \cdot C_{i,\text{sil\_melt}}^{\text{eq}}$$
- **Physical Saturation Ceilings and Four-Case Resolution:**
  Both reservoirs possess physical saturation limits:
  - Silicate melt ceilings ($C_{i,\text{sil\_melt\_max}}$): $1.0 \times 10^6\text{ ppmw}$ for C, N, and S; $100.0\text{ wt}\%$ for $\text{H}_2\text{O}$ ($= 100.0 \times f_H\text{ ppmw H}$).
  - Metal alloy ceilings ($C_{i,\text{met\_max}}$): $7.0 \times 10^4\text{ ppmw}$ for C, $4.0 \times 10^4\text{ ppmw}$ for N, $3.65 \times 10^5\text{ ppmw}$ for S (Fe-FeS eutectic), and $1.0 \times 10^4\text{ ppmw}$ for H.
  The solver resolves saturation across four mutually exclusive regimes:
  1. *Unconstrained*: Neither ceiling binds; concentrations follow the Nernst law.
  2. *Metal saturation only*: Metal alloy saturates at $C_{i,\text{met\_max}}$; the remainder resides in the silicate melt.
  3. *Silicate saturation only*: Silicate melt saturates at $C_{i,\text{sil\_melt\_max}}$; excess volatile mass partitions into the metallic liquid up to $C_{i,\text{met\_max}}$.
  4. *Dual saturation ($M_{i,\text{tot}} > M_{i,\text{sil\_max}} + M_{i,\text{met\_max}}$)*: Liquid metal saturates at $C_{i,\text{met\_max}}$, while excess volatile mass is retained within the silicate array. An atomic warning counter (`METAL_SILICATE_CAP_WARNING_COUNTER`) records the event in telemetry.
- **Strict Mass Conservation via Mirrored Writes:**
  Kinetic relaxation with rate $\alpha_{\text{eq}} \in [0, 1]$ updates the metal concentration by $\Delta C_{i,\text{met}} = \alpha_{\text{eq}} (C_{i,\text{met}}^{\text{eq}} - C_{i,\text{met}})$. Silicate concentration updates mirror the applied metal change:
  $$\Delta C_{i,\text{sil\_bulk}} = -\Delta C_{i,\text{met}} \left(\frac{m_{\text{met}}}{m_{\text{sil}}}\right)$$
  This guarantees whole-marker elemental mass invariance to machine precision ($< 10^{-12}$) across all saturation regimes.
- **Dynamic Sulfur Density Feedback:**
  $$\rho_{\text{metal}}(w_S) = 7020.0 - 5050.0 \cdot w_S \quad [\text{kg/m}^3]$$

---

## 3. Benchmark Validation

![Metal-Silicate Volatile Partitioning and Core Geochemistry Benchmark](../assets/core_geochemistry_benchmark.png)

*Figure: Metal-silicate volatile partitioning and core geochemistry benchmark in Erebus.jl. (a) Oxygen fugacity sensitivity of partition coefficients $D_i^{\text{met/sil}}$ over the range $\Delta\text{IW} \in [-4, 0]$ at $T = 1600\text{ K}$, $P = 0.1\text{ GPa}$, and $w_S = 0.05$. (b) Suppression of carbon partition coefficient $D_C$ by dissolved sulfur in metallic liquid for $w_S \in [0, 0.31]$ compared to nitrogen $D_N$, showing a steep drop in $D_C / D_N$ from $>100$ in sulfur-free metal down to $\sim 1$ at the Fe-FeS eutectic. (c) Temperature dependence of partition coefficients from $1300\text{ K}$ to $2200\text{ K}$ at $\Delta\text{IW} = -2.0$ and $w_S = 0.10$. (d) Core carbon versus nitrogen concentrations predicted by Erebus.jl over varied oxygen fugacities compared against empirical fields for magmatic iron meteorite groups (IIAB, IIIAB, IVA, and IVB). (e) Integrated core volatile delivery timeline during runaway core formation in a $R = 50\text{ km}$ planetesimal. (f) Total planetary elemental mass allocation among segregated core, retained silicate mantle, and degassed/vented losses.*

### Analysis of Benchmark Results

1. **Redox Sensitivity (Panel a)**: Siderophile volatiles exhibit distinct sensitivities to ambient oxygen fugacity. As $\Delta\text{IW}$ increases from $-4$ (highly reducing) to $0$ (oxidizing), both $D_C$ and $D_N$ decrease systematically due to lower siderophile affinity under oxidizing conditions. Sulfur partitioning remains strongly favorable to metal throughout the range ($D_S > 100$), while hydrogen remains largely lithophile ($D_H \approx 0.2 - 0.5$).
2. **Sulfur Suppression of Carbon (Panel b)**: Dissolved sulfur strongly suppresses carbon partitioning. In sulfur-free liquid iron, $D_C \approx 2800$, whereas at the Fe-FeS eutectic ($w_S \approx 0.31$), $D_C$ falls to $\approx 18$. In contrast, nitrogen affinity declines modestly from $D_N \approx 36$ to $\approx 14$. This differential behavior drives the $D_C / D_N$ ratio from $78$ down to $1.3$. This shift explains why sulfur-rich parent bodies (group IIAB) retain carbon in their silicate mantles while sulfur-poor bodies sequester carbon efficiently into the metallic core.
3. **Thermal Sensitivity (Panel c)**: Carbon and nitrogen partition coefficients decrease with increasing temperature, whereas sulfur partition coefficients increase following Boujibar et al. (2014). Hydrogen exhibits weak temperature dependence.
4. **Iron Meteorite Geochemical Matching (Panel d)**: Model equilibrium tracks over varied oxygen fugacity ($\Delta\text{IW} \in [-3.5, -1.0]$) and sulfur contents overlap directly with the measured compositions of magmatic iron meteorites. Sulfur-rich simulations ($w_S = 0.15$) naturally reproduce the low-C, moderate-N signature of IIAB irons, whereas intermediate sulfur models ($w_S = 0.07$) match the IIIAB group. Refractory groups IVA and IVB correspond to volatile-depleted starting inventories.
5. **Core Segregation Dynamics (Panel e)**: Metal-hosted volatiles migrate into the central core concurrently with metallic iron during runaway differentiation ($t \sim 1.05 - 1.6\text{ Ma}$). Over $90\%$ of total core volatile delivery occurs within the interval of peak magma ocean settling. Core volatile inventories are integrated using canonical per-marker 3D volume weighting $V_m = A_m L(r_m) = A_m (2 r_m)$.
6. **Planetary Budget Allocation (Panel f)**: In closed-system differentiation with core mass fraction $x_{\text{met}} \approx 0.22$, sulfur, carbon, and nitrogen are heavily sequestered into the metallic core ($>95\%$). Hydrogen is lithophile ($D_H \approx 0.19$), retaining the majority of its budget in the silicate mantle ($74\%$) or degassing into the atmosphere ($21\%$), with only $5\%$ entering the metallic core.

---

## 4. Source Code Architecture

| Component | Source File | Functions & Structs |
|:---|:---|:---|
| Configuration Schema | `src/config/metal_partition.jl` | `MetalPartitionConfig`, `validate_config` |
| Partition Thermodynamics | `src/physics/metal_partitioning.jl` | `compute_metal_silicate_partition_coefficient`, `compute_metal_silicate_partition_coefficients` |
| Marker Phase Equilibration | `src/physics/metal_partitioning.jl` | `equilibrate_metal_silicate_volatiles!`, `get_metal_silicate_cap_warning_count`, `reset_metal_silicate_cap_warning_count!` |
| Core Inventory Integration | `src/physics/metal_partitioning.jl` | `compute_core_volatile_budgets` |
| Marker Arrays & Properties | `src/particles.jl` | `setup_marker_metal_volatile_properties`, `compute_marker_properties!`, `replenish_markers!` |
| Advective Transport | `src/numerics/darcy.jl` | `apply_metal_segregation!` |
| Simulation Integration | `src/simulation/step.jl` | Caching, equilibration calls, and checkpoint persistence |

---

## 5. References

- Boujibar, A., Andrault, D., Bolfan-Casanova, N., Bouhifd, M. A., & Kawamoto, T. (2014). Metal-silicate partitioning of sulphur, new experimental constraints by EMPA and SIMS. *Earth and Planetary Science Letters*, 391, 42-54. [https://doi.org/10.1016/j.epsl.2014.01.021](https://doi.org/10.1016/j.epsl.2014.01.021)
- Clesi, V., Bouhifd, M. A., Bolfan-Casanova, N., Manthilake, G., Schiavi, F., Kawamoto, T., & Andrault, D. (2018). Low hydrogen contents in Earth's core. *Science Advances*, 4(3), e1701876. [https://doi.org/10.1126/sciadv.1701876](https://doi.org/10.1126/sciadv.1701876)
- Fischer, R. A., Cottrell, E., Hauri, E., Lee, K. K. M., & Le Voyer, M. (2020). The partitioning of carbon and oxygen between core and mantle in the early Earth. *Proceedings of the National Academy of Sciences*, 117(16), 8743-8749. [https://doi.org/10.1073/pnas.1919930117](https://doi.org/10.1073/pnas.1919930117)
- Grewal, D. S., Dasgupta, R., Sun, C., Tsuno, K., & Costin, G. (2019a). Delivery of carbon, nitrogen, and sulfur to the silicate Earth by a planetary merger. *Science Advances*, 5(1), eaau3669. [https://doi.org/10.1126/sciadv.aau3669](https://doi.org/10.1126/sciadv.aau3669)
- Grewal, D. S., Dasgupta, R., & Farnell, A. (2019b). The speciation of carbon, nitrogen, and water in magma oceans and its effect on volatile partitioning between metal and silicate. *Geochimica et Cosmochimica Acta*, 251, 87-115. [https://doi.org/10.1016/j.gca.2019.02.009](https://doi.org/10.1016/j.gca.2019.02.009)

# Normative Accessory Mineral Tracking and Meteorite Diagnostics

This page documents the physical formulation, stoichiometric parameterizations, and numerical validation of normative accessory mineral tracking and meteorite parent body diagnostics in `Erebus.jl`. The model tracks sub-eutectic stoichiometric mineral allocation (troilite $\text{FeS}$, schreibersite $(\text{Fe,Ni})_3\text{P}$, cohenite $(\text{Fe,Ni})_3\text{C}$, graphite $\text{C}$, and nitrides $\text{Fe}_4\text{N}/\text{CrN}/\text{TiN}$) in metallic phases, models eutectic phase dissolution near $T_{\text{eutectic}} \approx 1213\text{ K}$, and computes regional modal distributions to classify meteorite parent body affinities.

---

## 1. Physical Motivation

The physical and chemical signatures of iron and stony-iron meteorites preserve records of planetary differentiation, core segregation, and thermal metamorphism in early planetesimals:

1. Magmatic iron meteorites (groups IIAB, IIIAB, IVA, IVB) represent fractional crystallization products of fully segregated metallic cores. In these bodies, extensive internal heating melted both metallic alloy and silicates. Molten metal collected into a central core, leaving silicate mantles and crusts depleted in metallic iron and accessory phases (Goldstein et al., 2009; Chabot and Drake, 1999).
2. Primitive iron and achondrite complexes (the IAB complex and winonaites) contain abundant angular silicate inclusions, high carbon abundances (cohenite and graphite), schreibersite, and troilite veins. Their textures and geochemistry indicate incomplete differentiation, partial melting, and limited core segregation in partially melted parent bodies (Benedix et al., 2000).

Tracking the spatial and thermal distribution of accessory phases provides a direct geochemical diagnostic link between numerical geodynamic simulations and laboratory meteorite petrology.

---

## 2. Mathematical Formulation

### Stoichiometric Accessory Phase Allocation

In solid metallic iron-nickel alloys below the eutectic temperature ($T \le T_{\text{eutectic}}$), minor and volatile elements (S, P, C, and N) exsolve into stoichiometric accessory minerals.

#### Troilite ($\text{FeS}$)

All available sulfur in the metallic phase forms stoichiometric troilite:

$$w_{\text{troilite}} = w_S \cdot \left(\frac{M_{\text{FeS}}}{M_S}\right) \approx 2.74162 \cdot w_S$$

$$w_{\text{Fe,troilite}} = w_S \cdot \left(\frac{M_{\text{Fe}}}{M_S}\right) \approx 1.74162 \cdot w_S$$

where $M_S = 32.065\text{ g/mol}$, $M_{\text{Fe}} = 55.845\text{ g/mol}$, and $M_{\text{FeS}} = 87.910\text{ g/mol}$.

#### Schreibersite ($(\text{Fe,Ni})_3\text{P}$)

Phosphorus combines with metallic iron and nickel to form schreibersite. With nickel fraction in metal $x_{\text{Ni}} \in [0, 1]$ (default $0.25$):

$$M_{\text{metal,avg}} = (1 - x_{\text{Ni}}) M_{\text{Fe}} + x_{\text{Ni}} M_{\text{Ni}}$$

$$M_{\text{schreibersite}} = 3 \cdot M_{\text{metal,avg}} + M_P$$

$$w_{\text{schreibersite}} = w_P \cdot \left(\frac{M_{\text{schreibersite}}}{M_P}\right) \approx 6.478 \cdot w_P$$

where $M_{\text{Ni}} = 58.6934\text{ g/mol}$ and $M_P = 30.97376\text{ g/mol}$.

#### Cohenite ($(\text{Fe,Ni})_3\text{C}$) and Graphite ($\text{C}$)

Carbon combines with iron to form cohenite up to the carbide saturation threshold $w_{C,\text{max}} \approx 0.0667$ ($6.67\text{ wt}\%$ C, corresponding to stoichiometric $\text{Fe}_3\text{C}$):

$$w_{\text{cohenite}} = \begin{cases}
w_C \cdot \left(\frac{3 M_{\text{Fe}} + M_C}{M_C}\right) \approx 14.948 \cdot w_C, & w_C \le w_{C,\text{max}} \\
w_{C,\text{max}} \cdot \left(\frac{3 M_{\text{Fe}} + M_C}{M_C}\right), & w_C > w_{C,\text{max}}
\end{cases}$$

$$w_{\text{graphite}} = \begin{cases}
0.0, & w_C \le w_{C,\text{max}} \\
w_C - w_{C,\text{max}}, & w_C > w_{C,\text{max}}
\end{cases}$$

Excess carbon above carbide saturation precipitates as elemental graphite, consistent with petrographic observations in primitive IAB meteorites.

#### Nitrides ($\text{Fe}_4\text{N}$, $\text{CrN}$, $\text{TiN}$)

Nitrogen allocates to accessory nitrides based on the selected `nitride_mode`:

1. `:roaldite` ($\text{Fe}_4\text{N}$, default): $w_{\text{nitride}} = w_N \cdot (4 M_{\text{Fe}} + M_N) / M_N \approx 16.948 \cdot w_N$.
2. `:carlsbergite` ($\text{CrN}$): $w_{\text{nitride}} = w_N \cdot (M_{\text{Cr}} + M_N) / M_N \approx 4.712 \cdot w_N$.
3. `:osbornite` ($\text{TiN}$): $w_{\text{nitride}} = w_N \cdot (M_{\text{Ti}} + M_N) / M_N \approx 4.417 \cdot w_N$.

#### Residual Metal Matrix and Mass Normalization

The remaining fraction of the solid metallic phase forms the metallic iron-nickel matrix:

$$w_{\text{matrix}} = \max\left(0.0, 1.0 - (w_{\text{troilite}} + w_{\text{schreibersite}} + w_{\text{cohenite}} + w_{\text{graphite}} + w_{\text{nitride}})\right)$$

When the sum of accessory mineral mass fractions exceeds $1.0$ (for extreme volatile enrichments), mineral fractions scale by $1.0 / \sum w_i$ so their total equals $1.0$, and $w_{\text{matrix}} = 0.0$.

---

### Thermal Eutectic Phase Dissolution

At elevated temperatures, solid accessory phases dissolve into eutectic metallic liquid. While the binary Fe-FeS eutectic temperature at low planetary pressure is approximately $1261\text{ K}$ ($988\ ^\circ\text{C}$), minor additions of nickel, phosphorus, and carbon depress the initial melting temperature to $\approx 1213\text{ K}$ ($940\ ^\circ\text{C}$; Chabot and Drake, 1999; Goldstein et al., 2009). The solid metal fraction $F_{\text{solid}}(T)$ is parameterized around eutectic temperature $T_{\text{eutectic}}$ (default $1213.0\text{ K}$) and transition interval $\Delta T_{\text{transition}}$ (default $50.0\text{ K}$):

$$F_{\text{solid}}(T) = \begin{cases}
1.0, & T \le T_{\text{eutectic}} \\
1.0 - \frac{T - T_{\text{eutectic}}}{\Delta T_{\text{transition}}}, & T_{\text{eutectic}} < T < T_{\text{eutectic}} + \Delta T_{\text{transition}} \\
0.0, & T \ge T_{\text{eutectic}} + \Delta T_{\text{transition}}
\end{cases}$$

The abundance of each solid accessory mineral on Lagrangian markers scales with the solid metal fraction:

$$X_{m, i} = F_{\text{solid}}(T) \cdot w_i^{\text{stoich}}$$

where $i \in \{\text{troilite}, \text{schreibersite}, \text{cohenite}, \text{graphite}, \text{nitride}, \text{metal matrix}\}$.

---

### Regional Modal Distributions and Meteorite Classification

Planetesimal markers are grouped into three concentric zones based on normalized radial distance $r / R_{\text{planet}}$:
- Core region: $r / R_{\text{planet}} \le r_{\text{core,norm}}$ (default $0.50$).
- Mantle region: $r_{\text{core,norm}} < r / R_{\text{planet}} < r_{\text{mantle,norm}}$ (default $0.85$).
- Crust region: $r / R_{\text{planet}} \ge r_{\text{mantle,norm}}$.

For each zone, volume-weighted mean abundances of accessory phases and metallic melt fraction are computed. The resulting planetary state is classified into diagnostic meteorite categories:

- Magmatic differentiated signature ($f_{\text{molten,core}} \ge 0.80$ and $f_{\text{metal,core}} \ge 0.40$): Represents bodies that underwent complete segregation into a liquid metallic core, characteristic of group IIIAB, IVA, and IVB iron meteorites.
- Primitive incomplete signature ($f_{\text{molten,core}} \le 0.60$ and $f_{\text{solid,crust}} \ge 0.01$): Represents bodies with partial melting, retained crustal accessory minerals, and incomplete metal segregation, characteristic of IAB complexes and winonaites.
- Transitional signature: Intermediate differentiation states with partial core segregation.

---

## 3. Benchmark Validation

![Normative Accessory Mineral Tracking and Meteorite Diagnostics Benchmark](../assets/mineral_assemblage_benchmark.png)

*Figure: Benchmark results for normative accessory mineral tracking and meteorite diagnostics in Erebus.jl. Panel (a) shows stoichiometric accessory mineral conversion factors as a function of precursor element content in metallic alloy. Panel (b) shows thermal eutectic phase dissolution showing progressive melting of troilite, schreibersite, and cohenite through the eutectic transition ($T_{\text{eutectic}} = 1213\text{ K}$, $\Delta T = 50\text{ K}$). Panel (c) shows the radial mineral assemblage in a differentiated 50 km radius planetesimal, displaying the molten metallic core, depleted mantle, and accessory-rich primitive crust. Panel (d) shows meteorite parent body diagnostic regimes comparing core melt fraction against crustal accessory retention, with petrologic fields for magmatic iron groups (IIIAB, IVA, IVB) and primitive complexes (IAB, winonaites).*

### Analysis of Benchmark Results

1. Panel (a) shows that linear stoichiometric factors convert precursor element concentrations into accessory phase abundances. Troilite conversion matches the exact molar ratio $M_{\text{FeS}} / M_S \approx 2.742$. Cohenite tracks carbon content up to $6.67\text{ wt}\%$ C; beyond this saturation point, excess carbon precipitates as graphite.
2. Panel (b) illustrates that as temperature rises above $1213\text{ K}$, solid accessory phases dissolve proportionally into molten Fe-FeS liquid alloy. Above $1263\text{ K}$, solid accessory minerals are completely dissolved into metallic liquid.
3. Panel (c) demonstrates that in a differentiated planetesimal with a conductive cooling lid, the hot interior ($T > 1263\text{ K}$) segregates a molten metallic core ($r \le 25\text{ km}$). The conductive crust ($r \ge 42.5\text{ km}$) remains below $1213\text{ K}$, preserving pristine troilite, schreibersite, and cohenite.
4. Panel (d) shows that model trajectories separate magmatic irons (high core melt fraction, zero crustal accessory retention) from primitive IAB and winonaite complexes (low core melt fraction, high crustal accessory retention), establishing a quantitative framework for matching simulation outputs to meteorite collections.

### Modeling Simplifications and Limitations

- Shared dissolution interval: All accessory minerals (troilite, schreibersite, cohenite, graphite, and nitrides) dissolve across the shared temperature interval $[T_{\text{eutectic}}, T_{\text{eutectic}} + \Delta T_{\text{transition}}]$. In multicomponent metallic systems, refractory graphite and carbides exhibit higher thermal stability than sulfides and dissolve according to composition-dependent liquidus curves.
- Linear melt fraction parameterization: Solid metal fraction scales linearly across the melting interval rather than following non-linear thermodynamic lever-rule trajectories.

---

## 4. Source Code Architecture

| Component | Source File | Functions and Structs |
|:---|:---|:---|
| Configuration Schema | `src/config.jl` | `PhaseTrackingConfig`, `validate_config` |
| Stoichiometry and Melting | `src/physics.jl` | `compute_troilite_stoichiometry`, `compute_schreibersite_stoichiometry`, `compute_cohenite_graphite_stoichiometry`, `compute_nitride_stoichiometry`, `compute_normative_mineral_assemblage` |
| Regional Diagnostics | `src/physics.jl` | `compute_regional_mineral_modes` |
| Marker Phase Properties | `src/particles.jl` | `setup_marker_phase_tracking_properties`, `compute_marker_properties!`, `replenish_markers!` |
| Simulation Integration | `src/simulation.jl` | Marker allocation, periodic diagnostic logging, and checkpoint persistence |

---

## 5. References

- Benedix, G. K., McCoy, T. J., Keil, K., & Bogard, D. D. (2000). A petrologic and geochemical study of winonaites: Implications for trace element behavior during primitive achondrite differentiation. *Geochimica et Cosmochimica Acta*, 64(14), 2535-2553. [https://doi.org/10.1016/S0016-7037(00)00383-5](https://doi.org/10.1016/S0016-7037(00)00383-5)
- Chabot, N. L., & Drake, M. J. (1999). Crystallization of magmatic iron meteorites: The role of phosphorus and sulfur. *Meteoritics & Planetary Science*, 34(2), 235-246. [https://doi.org/10.1111/j.1945-5100.1999.tb01749.x](https://doi.org/10.1111/j.1945-5100.1999.tb01749.x)
- Goldstein, J. I., Scott, E. R. D., & Chabot, N. L. (2009). Iron meteorites: Crystallization, thermal history, parent bodies, and origin. *Chemie der Erde - Geochemistry*, 69(4), 293-325. [https://doi.org/10.1016/j.chemer.2009.01.002](https://doi.org/10.1016/j.chemer.2009.01.002)

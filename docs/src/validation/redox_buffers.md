# Redox Buffers and Evans (2012) Electron Budget Accounting

This page documents the thermodynamic formulations, buffer conversions, and electron conservation framework implemented in `src/redox.jl` of `Erebus.jl`.

---

## 1. Physical Motivation

Oxygen fugacity ($f_{\text{O2}}$) controls volatile speciation, mineral stability, metal-silicate partitioning, and outgassing compositions in planetesimals and protoplanets. In early Solar System materials, redox states span a wide dynamic range, from highly reduced enstatite chondrites ($\Delta\text{IW} \approx -7$ to $-4$) to moderately reduced ordinary and carbonaceous chondrites ($\Delta\text{IW} \approx -3$ to $0$), up to oxidized terrestrial magmas ($\Delta\text{IW} \approx +3$ to $+5$, or near the Quartz-Fayalite-Magnetite buffer).

In previous versions of `Erebus.jl`, redox state was parameterized as a single static scalar ($\Delta\text{IW}$). However, geochemical processes require:
1. Seamless bidirectional translation between conventional petrologic buffers (Iron-Wüstite, Quartz-Fayalite-Magnetite, Nickel-Nickel Oxide, Magnetite-Hematite, Wüstite-Magnetite, Quartz-Iron-Fayalite, and Graphite-CO-CO2).
2. Dynamic local buffer determination from phase assemblages across individual grid cells.
3. Rigorous conservation of oxidation-reduction potential during open- and closed-system thermochemical evolution (serpentinization, core formation, and gas venting), following the extensive electron budget framework of Evans (2012).

---

## 2. Mathematical Formulation

### 2.1 Solid-Oxide Redox Buffers

For linear solid-oxide buffers, absolute oxygen fugacity is calculated following Frost (1991):

$$\log_{10} f_{\text{O2}} = \frac{A}{T} + B + C \frac{P_{\text{bar}} - 1}{T}$$

where $T$ is temperature in Kelvin, $P_{\text{bar}}$ is pressure in bar ($P_{\text{bar}} = P_{\text{Pa}} \times 10^{-5}$), and $A, B, C$ are calibrated thermodynamic constants.

| Buffer | Reaction | $A$ | $B$ | $C$ | Reference |
|:---|:---|:---|:---|:---|:---|
| **IW** (Iron-Wüstite) | $2\text{Fe} + \text{O}_2 \rightleftharpoons 2\text{FeO}$ | $-28164.0$ | $6.541$ | $0.0$ | Campbell et al. (2009); O'Neill (1988) |
| **QFM** (Quartz-Fayalite-Magnetite) | $3\text{Fe}_2\text{SiO}_4 + \text{O}_2 \rightleftharpoons 2\text{Fe}_3\text{O}_4 + 3\text{SiO}_2$ | $-25096.3$ | $8.735$ | $0.110$ | Frost (1991); O'Neill (1987) |
| **NNO** (Nickel-Bunsenite) | $2\text{Ni} + \text{O}_2 \rightleftharpoons 2\text{NiO}$ | $-24930.0$ | $9.360$ | $0.046$ | Frost (1991); O'Neill & Pownceby (1993) |
| **MH** (Magnetite-Hematite) | $4\text{Fe}_3\text{O}_4 + \text{O}_2 \rightleftharpoons 6\text{Fe}_2\text{O}_3$ | $-25497.5$ | $14.330$ | $0.019$ | Frost (1991); Chou (1978) |
| **WM** (Wüstite-Magnetite) | $6\text{Fe}_{1-x}\text{O} + \text{O}_2 \rightleftharpoons 2\text{Fe}_3\text{O}_4$ | $-32807.0$ | $13.012$ | $0.083$ | Frost (1991) |
| **QIF** (Quartz-Iron-Fayalite) | $2\text{Fe} + \text{SiO}_2 + \text{O}_2 \rightleftharpoons \text{Fe}_2\text{SiO}_4$ | $-29435.7$ | $7.391$ | $0.044$ | Frost (1991) |

### 2.2 Graphite-CO-CO2 Buffer (CCO)

The CCO equilibrium depends on total gas pressure $P_{\text{bar}} = f_{\text{CO}} + f_{\text{CO2}}$. The equilibrium constants for:

$$\text{C} + \frac{1}{2}\text{O}_2 \rightleftharpoons \text{CO}, \quad \log_{10} K_{\text{CO}} = \frac{5785.0}{T} + 4.545$$

$$\text{C} + \text{O}_2 \rightleftharpoons \text{CO}_2, \quad \log_{10} K_{\text{CO2}} = \frac{20590.0}{T} - 0.043$$

give the quadratic equation in $x = \sqrt{f_{\text{O2}}}$:

$$K_{\text{CO2}} x^2 + K_{\text{CO}} x - P_{\text{bar}} = 0$$

with unique positive physical solution:

$$x = \frac{-K_{\text{CO}} + \sqrt{K_{\text{CO}}^2 + 4 K_{\text{CO2}} P_{\text{bar}}}}{2 K_{\text{CO2}}}$$

$$\log_{10} f_{\text{O2}} = 2 \log_{10}(x)$$

### 2.3 Evans (2012) Electron Budget Framework

Oxygen fugacity is an intensive variable that depends on temperature, pressure, and local mineral assemblages. In contrast, the redox budget ($RB$) is an extensive, conserved quantity representing the moles of electrons required to bring an elemental assemblage to a specified reference state:

$$RB = \sum_i n_i \nu_i \quad [\text{mol } e^-]$$

where $n_i$ is the molar inventory of species $i$, and $\nu_i$ is the number of electrons required to bring species $i$ to the reference oxidation state:

$$\nu_i = z_i - z_{\text{ref}, i}$$

In `Erebus.jl`, two standard reference states from Evans (2012) are supported:
1. **Mantle Reference State ($\text{M}$)**: $\text{Fe}^{2+}, \text{C}^0, \text{S}^{2-}, \text{H}^+, \text{O}^{2-}, \text{P}^{5+}$. In this state, the background mantle minerals have $RB = 0$.
2. **Crust Reference State ($\text{C}$)**: $\text{Fe}^{3+}, \text{C}^{4+}, \text{S}^{6+}, \text{H}^+, \text{O}^{2-}, \text{P}^{5+}$.

The specific redox budget is normalized by total system mass $M$:

$$RB_M = \frac{RB}{M} \quad [\text{mol } e^- / \text{kg}]$$

---

## 3. Conservation Theorems and Verification

### 3.1 Serpentinization Reaction

Hydrothermal alteration oxidizes ferrous iron in olivine/pyroxene to ferric iron in magnetite while reducing aqueous protons to molecular hydrogen:

$$3\text{FeO} + \text{H}_2\text{O} \to \text{Fe}_3\text{O}_4 + \text{H}_2$$

Under the mantle reference state:
- $3\text{FeO}$ ($\text{Fe}^{2+}$) has $\nu = 0$.
- $\text{H}_2\text{O}$ has $\nu = 0$.
- $\text{Fe}_3\text{O}_4$ contains $2\text{Fe}^{3+}$ ($\nu = +1$ each, total $+2$).
- $\text{H}_2$ contains $2\text{H}^0$ ($\nu = -1$ each, total $-2$).

Total electron balance is identically zero: $\Delta RB = (+2) + (-2) = 0$.

### 3.2 Core Segregation

Segregation of metallic liquid into an isolated core concentrates reduced metallic phases ($\text{Fe}^0, \text{Fe}_3\text{C}, \text{P}_{\text{phosphide}}$):

$$RB_{\text{bulk}} = RB_{\text{mantle}} + RB_{\text{core}}$$

Because the core sequesters strongly negative electron budget ($RB_{\text{core}} < 0$), the residual silicate mantle becomes oxidized ($RB_{\text{mantle}} > RB_{\text{bulk}}$).

### 3.3 Degassing and Venting

Loss of reduced volatile gases ($\text{H}_2, \text{CO}, \text{CH}_4$) through porous vents removes electron equivalents from the residual planetesimal:

$$RB_{\text{rock, initial}} = RB_{\text{rock, final}} + RB_{\text{gas, vented}}$$

---

## 4. Verification Test Suite

The redox engine is validated in `test/test_redox.jl`:

| Testset | Target Invariant | Tolerance / Assertion |
|:---|:---|:---|
| `Redox Buffer Coefficients and Equilibrium Physics` | Exact legacy IW reproduction; Petrologic order (MH > NNO > QFM > WM > IW > QIF); 3-class guards | `isapprox(atol=1e-12)`; sign and scale bounds |
| `Graphite CCO Buffer Inversion Physics` | Inversion recovers total gas pressure $f_{\text{CO}} + f_{\text{CO2}} = P_{\text{bar}}$ | `isapprox(rtol=1e-6)` |
| `Bidirectional Buffer Translation and Invariants` | Exact round-trip identity; Triangle closure IW $\to$ QFM $\to$ NNO $\to$ IW | `isapprox(atol=1e-12)` |
| `Local Controlling Buffer Regime Selection` | Phase-dependent buffer selection | Symbolic equality (`===`) |
| `Evans 2012 Redox Budget Electron Accounting` | Electron conservation across serpentinization, core formation, and gas venting | `isapprox(atol=1e-12)` |

---

## 5. References

- Campbell, A. J., Danielson, L., Righter, K., Seagle, C. T., Wang, Y., & Prakapenka, V. B. (2009). High pressure effects on the iron-wüstite and nickel-nickel oxide oxygen fugacity buffers. *Earth and Planetary Science Letters*, 286(3-4), 556-564. [https://doi.org/10.1016/j.epsl.2009.07.022](https://doi.org/10.1016/j.epsl.2009.07.022)
- Evans, K. A. (2006). Redox decoupling and redox budgets: Conceptual tools for the study of earth systems. *Geology*, 34(6), 489-492. [https://doi.org/10.1130/G22472.1](https://doi.org/10.1130/G22472.1)
- Evans, K. A. (2012). The redox budget of subduction zones. *Earth-Science Reviews*, 113(1-2), 11-32. [https://doi.org/10.1016/j.earscirev.2012.03.003](https://doi.org/10.1016/j.earscirev.2012.03.003)
- Frost, B. R. (1991). Introduction to oxygen fugacity and its petrologic importance. In D. H. Lindsley (Ed.), *Oxide Minerals: Petrologic and Magnetic Significance* (Reviews in Mineralogy, Vol. 25, pp. 1-9). Mineralogical Society of America. [https://doi.org/10.1515/9781501508684-004](https://doi.org/10.1515/9781501508684-004)

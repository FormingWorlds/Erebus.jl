# Permeability and Hydrofracture Validation

This module validates the porosity-dependent matrix permeability and the effective stress-dependent hydrofracturing enhancement.

---

## 1. Kozeny-Carman Porosity-Permeability Relation

### Governing Formulation
Pore fluid percolation through the solid planetesimal matrix follows Darcy's law with permeability $k(\phi)$ governed by the Kozeny-Carman relationship:

$$k(\phi) = k_{\phi 0} \left(\frac{\phi}{\phi_0}\right)^3 \left(\frac{1 - \phi_0}{1 - \phi}\right)^2$$

where $k_{\phi 0}$ is reference permeability at reference porosity $\phi_0$.

### Literature Anchors
- **Carman, P. C. (1937)**. Fluid flow through granular beds. *Transactions of the Institution of Chemical Engineers*, 15, 150-166.  
  [https://doi.org/10.1016/S0263-8762(97)80003-2](https://doi.org/10.1016/S0263-8762(97)80003-2)
- **Hubmann, B. (2022)**. *Hydrology of Planetesimals*. Master's thesis, ETH Zurich.  
  [https://doi.org/10.5281/zenodo.7058229](https://doi.org/10.5281/zenodo.7058229) (Equation 2.11)
- **Gerya, T. (2019)**. *Introduction to Numerical Geodynamic Modelling* (2nd ed.). Cambridge University Press.  
  [https://doi.org/10.1017/9781316534243](https://doi.org/10.1017/9781316534243)

### Invariants and Limits
1. **Zero Porosity Limit**: As $\phi \to 0$, $k(\phi) \to 0$ (impermeable solid).
2. **Reference Consistency**: At $\phi = \phi_0$, $k(\phi_0) = k_{\phi 0}$.
3. **Monotonicity**: Permeability is strictly monotonically increasing with porosity $\phi$ for all $\phi \in (0, 1)$.
4. **Error Contract**: Passing $\phi < 0$ or $\phi \ge 1$ throws `DomainError`.

### Verification Test Suite
- `test/test_physics.jl`: `@testset "kphi()"`

---

## 2. Terzaghi Effective Overpressure and Hydrofracturing

### Governing Formulation
When pore fluid pressure $P_f$ exceeds solid confining pressure $P_s$ plus tensile strength $\sigma_T$, rock fractures and increases effective permeability:

$$\Delta P_{\text{eff}} = P_f - P_s - \sigma_T$$

$$k_{\text{eff}} = \min\left( k(\phi) \left[ 1 + \left(\frac{\max(0, \Delta P_{\text{eff}})}{\sigma_{\text{scale}}}\right)^\alpha \right], k_{\text{max}} \right)$$

### Literature Anchors
- **Terzaghi, K. (1925)**. *Erdbaumechanik auf bodenphysikalischer Grundlage*. Franz Deuticke, Leipzig.
- **Wang, H. F. (2000)**. *Theory of Linear Poroelasticity with Applications to Geomechanics and Hydrogeology*. Princeton University Press.

### Invariants and Limits
1. **No Overpressure**: When $P_f \le P_s + \sigma_T$, $\Delta P_{\text{eff}} \le 0$ and $k_{\text{eff}} = k(\phi)$.
2. **Strict Upper Bound**: $k_{\text{eff}} \le k_{\text{max}}$ under arbitrarily high fluid overpressure.
3. **Monotonic Enhancement**: $k_{\text{eff}}$ increases monotonically with pore fluid overpressure $\Delta P_{\text{eff}} > 0$.

### Verification Test Suite
- `test/test_physics.jl`: Hydrofracturing permeability bounds
- `test/test_numerics.jl`: Stokes-Darcy coupled fluid-matrix pressure solve

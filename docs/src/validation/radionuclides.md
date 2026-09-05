# Radionuclide Heating Validation

This module validates the radiogenic heating rate calculations from short-lived radionuclides ($^{26}\text{Al}$ and $^{60}\text{Fe}$).

---

## Governing Formulation

Radioactive decay generates volumetric heating $Q(t)$ in the solid rock and fluid phases:

$$Q(t) = Q_0 \exp\left(-\frac{t}{\tau}\right) = Q_0 \left(\frac{1}{2}\right)^{t / t_{1/2}}$$

with mean lifetime $\tau = t_{1/2} / \ln 2$. The initial volumetric power density at time of CAI formation ($t = 0$) is:

$$Q_0 = \rho f_m \left(\frac{^{26}\text{Al}}{^{27}\text{Al}}\right)_0 E_{\text{al}} \frac{1}{\tau_{\text{al}}}$$

where $\rho$ is phase density, $f_m$ is elemental mass abundance, and $E_{\text{al}}$ is decay energy per atom.

---

## Literature Anchors

- **Russell, S. S., Srinivasan, G., Huss, G. R., & Wasserburg, G. J. (1996)**. Evidence for widespread $^{26}\text{Al}$ in the solar nebula and constraints for nebula time scales. *Science*, 273(5276), 757-762.  
  [https://doi.org/10.1126/science.273.5276.757](https://doi.org/10.1126/science.273.5276.757)
- **Tachibana, S., & Huss, G. R. (2003)**. The initial abundance of $^{60}\text{Fe}$ in the Solar System. *The Astrophysical Journal*, 588(1), L41-L44.  
  [https://doi.org/10.1086/374597](https://doi.org/10.1086/374597)
- **Lichtenberg, T., Golabek, G. J., Burn, R., Meyer, M. R., Alibert, Y., Gerya, T. V., & Mordasini, C. (2019)**. A water budget dichotomy of rocky protoplanets from 26Al-heating. *Nature Astronomy*, 3(4), 307-313.  
  [https://doi.org/10.1038/s41550-018-0688-5](https://doi.org/10.1038/s41550-018-0688-5)

---

## Invariants and Analytical Limits

1. **Total Energy Conservation Closure**:
   The time-integrated radiogenic energy released per unit volume over all time equals the total isotopic energy budget:

   $$\int_0^\infty Q(t) dt = Q_0 \tau = \rho f_m \left(\frac{^{26}\text{Al}}{^{27}\text{Al}}\right)_0 E_{\text{al}}$$

2. **Half-Life Consistency**:
   At $t = t_{1/2}$, $Q(t_{1/2}) = 0.5\, Q_0$ to floating point precision.
   At $t = 2 t_{1/2}$, $Q(2 t_{1/2}) = 0.25\, Q_0$.
3. **Strict Monotonicity**:
   $dQ/dt < 0$ for all $t \ge 0$.
4. **Positivity**:
   $Q(t) > 0$ for all finite $t \ge 0$.
5. **Error Contract**:
   Passing negative times ($t < 0$) or negative isotope ratios throws `DomainError`.

---

## Verification Test Suite

- `test/test_physics.jl`: `@testset "Q_radiogenic(): half-life and conservation closure"`

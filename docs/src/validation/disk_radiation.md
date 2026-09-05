# Surface Radiation and Disk Evolution Validation

This module validates the Stefan-Boltzmann surface radiation boundary condition and ambient protoplanetary disk thermal evolution models.

---

## 1. Stefan-Boltzmann Surface Radiative Boundary

### Governing Formulation
The boundary heat flux between the planetesimal surface and the surrounding disk gas is:

$$F_{\text{rad}} = \epsilon \sigma_{\text{SB}} \left( T_{\text{surf}}^4 - T_{\text{disk}}^4 \right)$$

Linearized as $F_{\text{rad}} = h_{\text{rad}} (T_{\text{surf}} - T_{\text{disk}})$ with:

$$h_{\text{rad}} = \epsilon \sigma_{\text{SB}} \left( T_{\text{surf}}^2 + T_{\text{disk}}^2 \right) \left( T_{\text{surf}} + T_{\text{disk}} \right)$$

At boundary faces, effective harmonic series conductivity couples internal conduction to radiative transfer:

$$k_{\text{interface}} = \frac{2 k_{\text{rock}} h_{\text{rad}} \Delta}{2 k_{\text{rock}} + h_{\text{rad}} \Delta}$$

### Literature Anchors
- **Gerya, T. (2019)**. *Introduction to Numerical Geodynamic Modelling* (2nd ed.). Cambridge University Press.  
  [https://doi.org/10.1017/9781316534243](https://doi.org/10.1017/9781316534243)

### Invariants and Limits
1. **Thermal Equilibrium**: When $T_{\text{surf}} = T_{\text{disk}}$, $F_{\text{rad}} = 0$ exactly.
2. **Directional Heat Transfer**: $T_{\text{surf}} > T_{\text{disk}} \implies F_{\text{rad}} > 0$ (cooling); $T_{\text{surf}} < T_{\text{disk}} \implies F_{\text{rad}} < 0$ (heating).
3. **Harmonic Bounding**: $k_{\text{interface}} \le 2 k_{\text{rock}}$ under arbitrarily intense radiation ($h_{\text{rad}} \to \infty$).
4. **Error Contract**: Passing non-positive absolute temperatures ($T \le 0\text{ K}$) or unphysical emissivity ($\epsilon \notin [0, 1]$) throws `DomainError`.

---

## 2. Protoplanetary Disk Temperature Evolution

### Governing Formulations
1. **Monotonic Viscous Clearing (`:monotonic`)**:

   $$T_{\text{disk}}(t, r, M_\star) = \left[ T_{\text{irr}}^4 + \max(0, T_{\text{peak}}^4 - T_{\text{irr}}^4) \left( 1 + \frac{t}{t_{\text{visc}}} \right)^{-\gamma} \right]^{1/4}$$

2. **Two-Stage Accretion Heating (`:class1_to_class2` or `:class0_to_class2`)**:

   $$T_{\text{disk}}(t, r, M_\star) = \left[ T_{\text{eff, irr}}(t)^4 + \max(0, T_{\text{peak}}^4 - T_{\text{irr}}^4) f_{\text{acc}}(t) \right]^{1/4}$$

### Literature Anchors
- **Chiang, E. I., & Goldreich, P. (1997)**. Spectral energy distributions of T Tauri stars with passive circumstellar disks. *The Astrophysical Journal*, 490(1), 368-376.  
  [https://doi.org/10.1086/304869](https://doi.org/10.1086/304869)
- **Drążkowska, J., & Dullemond, C. P. (2018)**. Planetesimal formation during protoplanetary disk buildup. *Astronomy & Astrophysics*, 614, A62.  
  [https://doi.org/10.1051/0004-6361/201732221](https://doi.org/10.1051/0004-6361/201732221)
- **Williams, J., Krijt, S., Drążkowska, J., & Lichtenberg, T. (2026)**. Planetesimal formation across the stellar mass spectrum and its influence on exoplanet-inherited volatile budgets. *Monthly Notices of the Royal Astronomical Society*, 551(3), stag1510.  
  [https://doi.org/10.1093/mnras/stag1510](https://doi.org/10.1093/mnras/stag1510)
- **Lynden-Bell, D., & Pringle, J. E. (1974)**. The evolution of viscous discs and the origin of the nebular variables. *Monthly Notices of the Royal Astronomical Society*, 168(3), 603-637.  
  [https://doi.org/10.1093/mnras/168.3.603](https://doi.org/10.1093/mnras/168.3.603)

### Invariants and Limits
1. **Asymptotic Convergence**: For both models, $\lim_{t \to \infty} T_{\text{disk}}(t, r) = T_{\text{irr}}(r)$.
2. **Molecular Cloud Floor**: $T_{\text{disk}} \ge T_{\text{cloud}} = 30.0\text{ K}$ across all radial distances and times.
3. **Viscous Peak Heating**: For Model 2, $T_{\text{disk}}(t_{\text{peak}}) = T_{\text{peak}}$ when $T_{\text{peak}} \ge T_{\text{irr}}$.
4. **Snowline Migration**: The water snowline ($T = 170\text{ K}$) expands outward during accretion peak and retreats inward during viscous clearing.

---

## Verification Test Suite

- `test/test_geometry_radiation.jl`:
  - `@testset "Stefan-Boltzmann Surface Radiation Physics"`
  - `@testset "Radiative Boundary Application on Staggered Grid"`
  - `@testset "Disk Temperature Model 1: Monotonic Viscous Decay"`
  - `@testset "Disk Temperature Model 2: Cold -> Hot -> Cold (Option 2A)"`
  - `@testset "Water Snowline Dynamics Across Evolutionary Stages"`

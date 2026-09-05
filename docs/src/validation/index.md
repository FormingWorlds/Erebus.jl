# Validation and Physical Benchmarks

The physical formulations in `Erebus.jl` are anchored against peer-reviewed literature, analytical limits, and laboratory measurements. This section provides an inventory of the physical benchmarks and verification tests governing simulation fidelity.

---

## Validation Matrix

| Physical Process | Governing Theory | Primary Reference | Verification Test Suite |
|:---|:---|:---|:---|
| **Thermal Conduction & Geometry** | Two-phase porous conductivity and 3D spherical metric divergence | Gerya (2019); Hubmann (2022) | `test/test_geometry_radiation.jl`, `test/test_physics.jl` |
| **Permeability & Hydrofracturing** | Kozeny-Carman flow and Terzaghi effective overpressure failure | Carman (1937); Terzaghi (1925); Wang (2000) | `test/test_physics.jl`, `test/test_numerics.jl` |
| **Radionuclide Decay** | Short-lived radioactive heating ($^{26}\text{Al}$, $^{60}\text{Fe}$) | Russell et al. (1996); Tachibana & Huss (2003); Lichtenberg et al. (2019) | `test/test_physics.jl` |
| **Fluid Viscosity & Phase Transitions** | Arrhenius water viscosity and hydrothermal phase limits | Hubmann (2022); Gerya (2019) | `test/test_physics.jl` |
| **Surface Radiation & Disk Evolution** | Stefan-Boltzmann boundary and protoplanetary disk clearing | Chiang & Goldreich (1997); Drążkowska & Dullemond (2018); Williams et al. (2026) | `test/test_geometry_radiation.jl` |

---

## Standards of Verification

Each validation module fulfills three requirements:
1. **Literature Grounding**: Every parameterization maps to a verified scientific publication with a valid DOI.
2. **Analytical Limits**: Numerical implementations are tested against closed-form mathematical limits (e.g. pure solid limit, infinite-time decay, steady-state conduction).
3. **Physical Invariants**: Unit tests assert physical conservation (energy, mass), boundedness ($T > 0\text{ K}$, $\phi \in [0, 1]$), and monotonicity.

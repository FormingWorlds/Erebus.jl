# Erebus.jl Project State & Plan

## 1. What Has Been Done So Far

We have addressed major physics flaws related to mass conservation, thermodynamic equilibrium, and volatile reference frames. These issues were initially identified by an adversarial review (Opus) and have now been fully resolved.

### Completed Fixes & Implementations:
- **Venting Mass Creation (Atmosphere Boundary)**
  - *Location*: `src/simulation/loop.jl`
  - *Fix*: Removed artificial stoichiometric multipliers (`44.0095 / 12.011` and `34.08 / 32.06`) for `:CO2` and `:H2S` when speciation is inactive. This prevents elemental mass from being artificially "synthesized" into molecular mass as it vents to the atmosphere.
- **Thermodynamic Hard Capping (Magma Ocean Degassing)**
  - *Location*: `src/physics/magma_ocean_degassing.jl`
  - *Fix*: Replaced independent, per-species elemental scaling caps (`scale_H`, `scale_C`, etc.) with a uniform `scale_global`. This ensures that when the massive gas inventory must be scaled to fit solver bounds, the H/C/N/S equilibrium ratios and Dalton's Law are strictly preserved.
- **Bulk vs. Melt Extraction (Magma Ocean Degassing)**
  - *Location*: `src/physics/magma_ocean_degassing.jl`
  - *Fix*: Scaled bulk marker concentrations by `1.0 / F_curr` to evaluate physical supersaturation strictly in the melt volume, then scaled the extracted mass back by `F_curr`.
- **Metal Partitioning Reference Frame**
  - *Location*: `src/physics/metal_partitioning.jl`
  - *Fix*: Corrected the computation of total available volatile mass (`M_tot`) by evaluating the true melt concentration `C_sil_melt = C_sil_bulk / F_melt_val` across the C, N, S, and H equilibration loops.
- **Test Suite Modernization**
  - *Locations*: `test/test_core_volatile_partitioning.jl`, `test/test_atmosphere.jl`, `test/test_dehydration_darcy_coupling.jl`
  - *Fix*: Updated test assertions and stripped out hardcoded mass inflation multipliers so that they align with the newly rigorous mass conservation logic.
- **Verification**
  - The entire test suite (`Pkg.test("Erebus")`) was run successfully, with all 403,315 tests passing.

## 2. Current State

- **Code Tip**: All fixes from the previous adversarial review are currently on the local `fix-architecture` branch (uncommitted).
- **Ongoing Check**: A new adversarial review (using four independent Claude Opus 4.8 agents) is currently executing in the background against this latest code tip. The lenses are targeting correctness, completeness/test-quality, design, and physics/chemistry conservation.

## 3. What Is Planned Further

### Immediate Next Steps:
1. **Opus 4.8 Adversarial Review Results**:
   - Wait for the ongoing Opus 4.8 background task to finish.
   - Review and consolidate the findings.
   - Address and implement fixes for any new issues uncovered by these Opus 4.8 agents.
2. **Commit and Pull Request**:
   - Commit the fully verified fixes.
   - Request explicit sign-off from the user to open a Pull Request (per the standing PR rules).
   - Once approved, draft the PR title and description using the `git-conventions` and `ai-check` skills.

### Backlog & Future Tasks:
3. **Address Solver Scalability**: 
   - Adjust the massive `Kcont` penalty term or adapt the matrix-free operator to better support stable hydrofracture numerics.
4. **Documentation Consolidation and Verification**:
   - Consolidate redundant explanations in the documentation and align it with the newly corrected volumetric scaling ratios.
5. **Repair Submodule Facades**:
   - Assess and repair the valid Julia facade pattern (pending user confirmation).

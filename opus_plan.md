I now have everything needed. Here is my analysis report.

---

# Analysis: three mass-conservation bugs in Erebus.jl

Read-only analysis. All line numbers are from the current `fix-architecture` working tree. Every claim below comes from reading the source; I ran no Julia and executed no test, so the failure magnitudes are analytic estimates, not measured. I flag confidence per item.

## Summary table

| # | Bug | File | Key lines | Physical error |
|---|-----|------|-----------|----------------|
| 1 | Dalton's-law / partial-pressure→mass conversion omits mean molar mass | `src/physics/magma_ocean_degassing.jl` | 223, 276, 430, 499 (+ elemental sums 279-296, 434-453) | Per-species atmospheric mass set ∝ partial pressure with one universal coefficient; distorts the elemental (H/C/N/S) split. |
| 2 | Silicate interacting mass uses bulk (not molten) silicate | `src/physics/metal_partitioning.jl` | 306 | Silicate volatile capacity inflated by 1/F_melt when partially molten; asymmetric with the metal side (line 305). |
| 3 | 2D→3D volume extrusion uses cylinder length 2Rp | `src/physics/magma_ocean_degassing.jl` (618, 723-726); `src/simulation/loop.jl` (3299, 3312, 3529) | Extruding the πR² disk by 2Rp gives 2πR³ vs. true sphere (4/3)πR³: a 1.5× overcount, inconsistent with the core budgets that use (4/3)πR³. |

---

## Bug 1: partial pressure → species mass conversion (Dalton's-law misapplication)

**Confidence: high on the physics; the exact preferred fix form is a design choice.**

### Location
`src/physics/magma_ocean_degassing.jl`, function `solve_magma_ocean_volatile_partitioning`.

- `col_coeff = area / grav` (`4πR²/g`) at line 146 is correct for the *total* column: `P_surf = M_atm_tot / col_coeff`.
- The error is applying that same coefficient per species:
  - line 223: `m_atm_dict = Dict(k => v * col_coeff …)` (M_m==0 branch)
  - line 276: `m_atm_dict = Dict(k => v * col_coeff …)` (inside `eval_partitioning_at_pressure`)
  - lines 430, 499: `final_p_i = Dict(k => v / col_coeff …)` (the inverse, same wrong coefficient)
- The distorted per-species masses then feed the elemental sums at lines 279-296 and 434-453, which drive the bisection target `m_tot_calc` (line 343) and the converged `P_surf`.

### Why it is wrong
Hydrostatics fixes the **total**: `P_surf = g·M_atm_tot/A`, independent of composition. For an individual species, mass is proportional to **moles**, and moles (not mass) are proportional to partial pressure:

- moles of species i in the column = `p_i · A / (g · μ̄)`, where `μ̄ = Σ p_i M_i / Σ p_i` is the mean molar mass.
- therefore `M_i = col_coeff · p_i · (M_i_molar / μ̄)`.

The code uses `M_i = col_coeff · p_i`, i.e. it drops the `M_i_molar/μ̄` factor. Equivalently, when it later divides by the species molar mass to get elements, it divides by each species' own `m_i` instead of by `μ̄`. The ratio (code/correct) for each species' contribution is `μ̄/m_i`.

Consequence: the **total** atmospheric mass is still `col_coeff·P_surf` (so a total-mass conservation test still passes), but the **distribution across H, C, N, S is biased**. In an H2/H2O-rich atmosphere (`μ̄` a few g/mol), heavy carriers such as CO2 (44 g/mol) have their elemental contribution scaled by roughly `μ̄/44`, an order-of-magnitude distortion. Because the bisection converges `m_tot_calc` to the true total using these biased elemental masses, the converged `P_surf` and melt/atmosphere split are physically wrong, not merely rescaled.

### Reference implementation already in the repo
`speciate_vented_volatiles` in `src/physics/solubility_speciation.jl` (lines 1559-1627) does it correctly: it forms mole fractions `y_i = p_i/p_sum`, gets total gas moles from elemental conservation, then `mass_i = N_gas · y_i · M_i_molar`. The degassing solver should mirror that logic.

### Proposed fix
In every place that maps partial pressures to species masses (and its inverse), introduce `μ̄` and weight by molar mass. Concrete form for `eval_partitioning_at_pressure` (replace line 276):

```julia
molar = Dict{Symbol,Float64}(
    :H2 => 2.01588e-3, :H2O => 18.01528e-3, :CO => 28.0101e-3,
    :CO2 => 44.0095e-3, :CH4 => 16.04246e-3, :N2 => 28.0134e-3,
    :NH3 => 17.03052e-3, :H2S => 34.08088e-3, :S2 => 64.12e-3, :SO2 => 64.066e-3,
)
p_tot = sum(values(p_dict))               # already rescaled to P_trial above
mu_bar = p_tot > 0.0 ? sum(p_dict[k] * molar[k] for k in keys(p_dict)) / p_tot : 0.0
m_atm_dict = Dict{Symbol,Float64}(
    k => (mu_bar > 0.0 ? v * col_coeff * (molar[k] / mu_bar) : 0.0) for (k, v) in p_dict
)
```

This preserves `Σ M_i = col_coeff · P_surf` (verify: `Σ col_coeff p_i M_i/μ̄ = col_coeff/μ̄ · Σ p_i M_i = col_coeff/μ̄ · P·μ̄ = col_coeff·P`) while giving correct per-species masses.

Apply the same molar-mass weighting to:
- the M_m==0 branch (line 223): `m_atm_dict = Dict(k => v * col_coeff * molar[k]/μ̄ …)`.
- the final reconstruction (lines 430 and 499): the inverse map must invert the *same* weighted relation, i.e. `p_i = M_i · μ̄ / (col_coeff · M_i_molar)`. Simplest: keep `p_i` from the converged `best_res.p_dict` directly and derive `M_i` from it with the weighted formula, rather than round-tripping through `/col_coeff`.

Define the `molar` dict once at module or function scope to avoid duplication (BlueStyle / DRY).

### Verification the caller should run
1. RED test first (per the repo TDD rule): construct a two-species case where `μ̄ ≠ m_i` (e.g. pure H2O vs. pure CO2 at fixed total moles) and assert the returned `M_atm_i` matches the mole-based masses from `speciate_vented_volatiles`. The current code fails this; the fix passes.
2. Independent check: confirm `sum(values(M_atm_i)) ≈ col_coeff * P_surf` still holds after the fix (total-mass invariant).
3. Confirm the existing conservation tests in `test/test_magma_degassing.jl` (e.g. the no-melt total at line 195) still pass; the per-species assertions there recompute elements with the same ratios, so they test total, not the split.

---

## Bug 2: inflated silicate volatile capacity (`m_sil`)

**Confidence: high.**

### Location
`src/physics/metal_partitioning.jl`, function `equilibrate_metal_silicate_volatiles!`, lines 305-306:

```julia
m_met = Xfem[m] * max(Float64(rho_metal), 100.0)      # molten metal volume fraction  → OK
m_sil = phi_sil * max(Float64(rho_silicate), 100.0)    # BULK silicate volume fraction → BUG
```

with `phi_sil = max(1 - phi_fe, 0)` (line 292) and `F_melt_val` already available (line 294).

### Why it is wrong
`Xfem[m]` is documented (line 243) as the **molten** metal volume fraction, so `m_met` correctly counts only molten metal. The silicate side uses `phi_sil`, the **total** silicate volume fraction, including solid silicate that does not chemically equilibrate with molten metal. The routine's own docstring (line 229-231) and guard (line 295, which requires `F_melt_val > 0`) state the process is equilibration between molten metal and **silicate melt**. The two phases are treated asymmetrically.

`m_sil` and `m_met` are the lever-rule weights in `C_sil_eq = M_tot/(m_sil + D·m_met)` (lines 333-334, 370-371, 407-408, 445-446). Overstating `m_sil` by `1/F_melt` (when `F_melt < 1`) biases the equilibrium toward the silicate: it inflates the silicate's apparent capacity to hold H, C, N, S and suppresses metal uptake, exactly the reported symptom.

### Proposed fix
Restrict the silicate mass to the molten fraction, symmetric with the metal side (replace line 306):

```julia
m_sil = phi_sil * F_melt_val * max(Float64(rho_silicate), 100.0)
```

The `phi_fe/phi_sil/F_fe/F_melt > 0` guard at line 295 already protects against a zero `m_sil` and the `m_sil <= 0.0` return at line 307 remains a valid backstop.

### Caveat to check
The routine redistributes per-marker *concentrations*; changing `m_sil` changes the lever weight and hence the equilibrium concentrations. Confirm that whatever consumes `XH2Om/XCm/XNm/XSm` downstream (the bulk-silicate transport and the atmosphere/venting coupling) interprets those arrays consistently as molten-silicate concentrations. If any downstream step maps them back over the *bulk* silicate mass, the coupling, not just this routine, must use the molten mass for a closed elemental budget. I did not trace the full downstream path; the caller should.

### Verification
1. RED test: partially molten marker (`F_melt = 0.5`, `phi_sil = phi_fe = 0.5`, `D` large). Assert the equilibrium metal concentration is higher (more siderophile drawdown) with the fix than with the current bulk-silicate weight, and that total elemental mass on the marker is unchanged (the routine is a two-reservoir redistribution, so `m_sil·C_sil + m_met·C_met` must be invariant under the α-step for each element).
2. Check `test/test_core_formation.jl` and `test/test_core_volatile_partitioning.jl` for assertions that assume the current `F_melt`-independent weight; those may go RED and need updating to the correct physics.

---

## Bug 3: cylinder (2Rp) 2D→3D extrusion vs. spherical geometry

**Confidence: high that 2Rp is a surface-area factor misapplied to a volume integral; medium on the single best replacement, because a uniform extrusion cannot reproduce a sphere exactly and the repo already carries a different (shell-equal-volume) convention.**

### Locations
- `src/physics/magma_ocean_degassing.jl:618`: `L_3D = 2.0 * Rp` and lines 723-726 multiply the volume-integrated extracted masses (`tot_ex_*`, summed over interior markers at lines 706-718) by `L_3D`.
- `src/simulation/loop.jl:3299, 3312, 3529`: `L_3D_equiv = 2.0 * rplanet_val` for venting/atmosphere coupling.
- Documented as `L_3D = 2 R_planet` in `docs/src/validation/volatile_retention.md:84`.

### Why it is wrong (and internally inconsistent)
The comment at `loop.jl:3298` states the intent: `4πR² / 2πR = 2R`. That is a **surface-flux** conversion (3D sphere surface / 2D circle perimeter) and is defensible for a quantity defined per unit length of the surface perimeter. But:

- `degas_magma_ocean_markers!` sums `tot_ex_*` over **interior/volume** markers (a magma-ocean shell between `r_degas` and `Rp`, lines 628-720), not a surface line. For a volume integral the correct 2D→3D map relates cross-sectional area (πR²) to volume ((4/3)πR³), giving an effective thickness `(4/3)R`, not `2R`.
- Extruding the πR² disk by `2R` yields a cylinder of volume `2πR³`, versus the true sphere `(4/3)πR³`: an overcount of `2πR³ / ((4/3)πR³) = 1.5`, i.e. **+50%**.
- This is inconsistent with the same module's core/mineral budgets, which use the true sphere: `compute_core_volatile_budgets` and `compute_regional_mineral_modes` set `V_tot = (4/3)π R³` and `V_m = V_tot/N_planet` (`metal_partitioning.jl:561, 1054`), and with the planetary mass `M_planet_val = (4/3)π R³ ρ_bulk` (`loop.jl:233`). The degassed/vented volatile inventory is therefore inflated by ~1.5× relative to the body it is supposed to be conserved against.

### Two fix options (the caller must pick the convention)

**Option A — minimal, keep the extrusion form, make it volume-consistent.** For the volume integrals (degassing, and the venting sink if that too is a volume-integrated bulk drainage), replace the length with `(4/3)·Rp`:

```julia
L_3D = (4.0 / 3.0) * Rp   # extrude πR² cross-section to (4/3)πR³ sphere-equivalent volume
```

This makes the degassed mass match the spherical planet volume in the aggregate. It is still only exact in the aggregate: a uniform extrusion does not reproduce the radial distribution of a revolved axisymmetric body (each marker at radius r should be weighted by its revolution, Pappus 2πr). Note that the `loop.jl:3298` surface-flux use of `2R` may be *correct as-is* if that specific term truly is a surface flux (`compute_surface_venting_rates` at line 3315), while the bulk drainage at `loop.jl:3300-3307` is a volume term that needs `(4/3)R`. These two uses must be separated; do not blanket-replace every `2.0*rplanet_val`.

**Option B — adopt the equal-volume marker convention used elsewhere (preferred for a closed budget).** Replace the per-area × length extrusion with the same `V_m = (4/3)πR³/N_planet` scheme that `compute_core_volatile_budgets` uses, so degassing, venting, and core budgets share one geometry. This removes the standing 1.5× inconsistency at the source and gives a single auditable convention. It is a larger change (the degassing sum currently works in per-out-of-plane-length units via `m_marker = rho_solid * marker_volume`).

Whichever option, `docs/src/validation/volatile_retention.md:84` and the `loop.jl:3298` comment must be corrected in the same change (docs assert geometry that would then be wrong).

### Verification
1. Analytic check: build a uniform planetesimal, degas a known bulk concentration, and confirm the summed 3D degassed mass equals `concentration × (4/3)πR³ρ` to the marker-discretisation tolerance. The current `2R` version overshoots by ~1.5×.
2. Cross-module conservation: sum (atmosphere + melt-retained + core + mantle + crust) elemental masses and confirm closure against the initial planetary inventory `(4/3)πR³ρ_bulk × w_bulk`. This is the decisive check and will expose whether venting's `2R` term is a genuine surface flux or a mislabelled volume term.
3. Separate the `loop.jl` uses: confirm by reading `compute_surface_venting_rates` whether its input is a per-perimeter surface flux (keep `2R`) or a bulk mass (needs `(4/3)R`). I did not fully trace that function.

---

## Cross-cutting notes for the caller

- **Independent verification is essential.** I established mechanism by reading, not by running. Before treating any of these as confirmed, run the RED tests above and the cross-module elemental-closure check; the repo's own rules require breaking the behaviour and observing the failure, not accepting a code-reading argument.
- **Interaction between the bugs.** They compound in the same elemental budget: bug 1 mis-splits atmospheric elements, bug 2 mis-splits metal/silicate, bug 3 mis-scales the absolute inventory by ~1.5×. A closure test will only pass once all three are consistent, so fix and test them together, then re-run the full closure.
- **Docs and tests to update:** `docs/src/validation/volatile_retention.md:84`, the `loop.jl:3298` comment, `test/test_magma_degassing.jl` (per-species split), and the core-formation tests. `test/test_dehydration_darcy_coupling.jl` passes `L_3D_equiv` as a free parameter (1.2), so it will not pin the geometry.
- **Not independently verified by me:** the downstream consumers of the silicate concentration arrays (bug 2 caveat), and whether `compute_surface_venting_rates` is a surface flux or a volume term (bug 3, `loop.jl` uses). Trace both before applying.

This is a lead for the caller to verify and a proposal to apply, not a completed fix.

---
description: Erebus test quality deep-dive. Anti-happy-path patterns, physical invariant families, 3-class discrimination guards, TDD RED-GREEN-REFACTOR cycle, and validation certification.
---

# Erebus Test Quality and Scientific Rigor Rules

Erebus simulates hydrothermal circulation, planetesimal interior thermodynamics, and early volatile evolution. The test suite is held to physics-grade rigor. Tests exist to catch real bugs. A test that asserts the wrong thing, or that passes for the wrong reason, is worse than no test because it creates false confidence.

These rules define what scientific test rigor means in Erebus.jl.

---

## 1. The Iron Law of Test-Driven Development (TDD)

Every new physical feature, boundary condition, solver modification, and bug fix in Erebus MUST follow the strict RED-GREEN-REFACTOR development cycle:

```
NO PRODUCTION CODE WITHOUT A FAILING TEST FIRST
```

### Phase 1: RED (Write the Failing Test)
- Write the minimal test expressing the physical invariant, analytical limit, boundary condition, or error contract before writing any implementation code.
- Execute the test and watch it fail.
- Verify that the test fails for the expected physical or mathematical reason, not because of a syntax error.

### Phase 2: GREEN (Make the Test Pass)
- Write the minimal implementation code to satisfy the test contract.
- Re-run the test suite and verify that the test passes.

### Phase 3: REFACTOR (Harden and Clean)
- Enforce BlueStyle formatting with `JuliaFormatter`.
- Ensure clean physical variable names and dimensional consistency.
- Eliminate duplication while preserving all test assertions.
- Verify that no performance regressions or extra memory allocations were introduced.

---

## 2. The Three Slogans (Invariants of Intent)

### Slogan 1: Physical Invariant Mandate
Every unit test on a physics routine must assert at least one physical invariant from the four families:
1. **Conservation**: Mass closure, energy conservation ($LHS = RHS$), flux continuity.
2. **Positivity and Boundedness**: Temperature $T > 0\text{ K}$, fluid pressure $P > 0\text{ Pa}$, porosity $\phi \in [0, 1]$, melt fraction in $[0, 1]$, finite floating-point values (no `NaN` or `Inf`).
3. **Monotonicity and Symmetry**: Pressure increases with depth ($dP/dr \le 0$), fluid viscosity decreases with temperature ($d\eta/dT < 0$), geometric flux operators preserve coordinate symmetries.
4. **Analytical Limits with 3-Class Discrimination Guards**: Pinned values against closed-form solutions with explicit guards for exponent, sign, and scale errors.

### Slogan 2: Zero Tautological Tests and Anti-Happy-Path Discipline
- Never test a function by evaluating a copy of its source equation with random numbers.
- Every test must exercise edge cases and boundary limits.
- Every function with a physical domain must enforce an error contract (`DomainError` or `ArgumentError`) when supplied unphysical inputs, and tests must verify this behavior with `@test_throws`.

### Slogan 3: Literature Grounding and Anti-Reference-Hallucination
- Every physical parameterization or equation must map to an authoritative reference with a verified DOI in `docs/src/validation/`.
- Test docstrings must cite specific equation, table, or figure numbers from the original publications.
- Never cite papers from memory. All citations must be verified through Crossref or OpenAlex.

---

## 3. Anti-Happy-Path Rules for Every Test

Every new test set or test function MUST include:

1. **At least one physical edge case**:
   - Boundary values: $\phi = 0.0$ (pure rock limit), $\phi = 1.0$ (pure fluid limit), $T = T_{\text{solidus}}$, $r = 0$ (center), $r = R_{\text{planet}}$ (surface), $r \gg R_{\text{planet}}$ (ambient disk).
2. **At least one path exercising the error contract**:
   - If a function has a physical domain ($T > 0$, $0 \le \phi \le 1$, $k > 0$, $\rho > 0$), passing unphysical inputs must throw `DomainError` or `ArgumentError`.
   - Tests must assert the exception: `@test_throws DomainError Erebus.ktotal(-1.0, 1.0, 0.5)`.
   - For closed-form algebra without exceptions, test the degenerate limit behavior (e.g. zero radioactive abundance yields zero radiogenic heating).
3. **Assertion values not trivially derivable from the implementation**:
   - Assert physical properties or analytical limits.
   - Do not write tests of the form `@test f(x, y) ≈ x * y + 2` where `x * y + 2` is copied from `src/physics.jl`.

### Forbidden Patterns (Flagged by Test Linter)

1. **Single-assert testsets**:
   - Every testset must contain at least 2 assertions. The second assertion typically pins the invariant that the first assertion tests broadly.
2. **Standalone weak assertions**:
   - `@test result !== nothing`
   - `@test result > 0` (when not part of a 3-class discrimination guard)
   - `@test length(result) > 0`
   - `@test typeof(result) == Float64`
3. **Float comparison with exact equality (`==`)**:
   - Never compare floating-point values using `==`. Always use `isapprox(val, expected; rtol=..., atol=...)` or `val ≈ expected rtol=...`.
4. **Untagged tests without physical descriptions**:
   - Testsets must document what physical scenario or equation is being tested.

---

## 4. Discriminating Numeric Values (3-Class Guards)

When a test pins a specific numerical value against an analytical solution or benchmark, it must include guards for the three most common bug classes:

1. **Exponent or Prefactor Error**:
   Assert that plausible wrong powers (e.g. $T^3$ instead of $T^4$, missing factor of 2 or $\pi$) fail the test:
   `@test abs(val - wrong_val) > threshold`
2. **Sign Error**:
   Assert the physical sign explicitly:
   `@test val > 0.0` or `@test val < 0.0`
3. **Scale or Unit Error**:
   Assert order-of-magnitude bounds to catch unit conversion bugs (e.g. bar vs Pa, years vs seconds, km vs m):
   `@test 1.0e-4 < val < 1.0e-2`

### Canonical Example

```julia
@testset "Stefan-Boltzmann Surface Radiation Physics" begin
    # Test at T_surf = 350 K, T_amb = 170 K with emissivity = 0.9
    T_s = 350.0
    T_a = 170.0
    eps = 0.9
    sigma = 5.670374419e-8
    h_rad = Erebus.compute_radiation_htc(T_s, T_a; emissivity=eps, sigma_sb=sigma)

    # Primary pin: h_rad * (T_s - T_a) == eps * sigma * (T_s^4 - T_a^4)
    flux_linear = h_rad * (T_s - T_a)
    flux_exact = eps * sigma * (T_s^4 - T_a^4)
    @test isapprox(flux_linear, flux_exact; rtol=1e-12)

    # 1. Limiting derivative behavior: as T_amb -> T_surf, h_rad -> 4 * eps * sigma * T^3
    T_close = T_s + 1e-6
    h_close = Erebus.compute_radiation_htc(T_s, T_close; emissivity=eps, sigma_sb=sigma)
    h_tangent = 4.0 * eps * sigma * (T_s^3)
    @test isapprox(h_close, h_tangent; rtol=1e-5)

    # 2. Sign and positivity guard
    @test h_rad > 0.0

    # 3. Scale guard: order of magnitude is ~3 W/(m² K)
    @test 1.0 < h_rad < 10.0

    # 4. Truncation and non-positive temperature handling: returns 0.0 smooth cutoff
    @test iszero(Erebus.compute_radiation_htc(0.0, 170.0))
    @test iszero(Erebus.compute_radiation_htc(300.0, -10.0))
    @test iszero(Erebus.compute_radiation_htc(NaN, 170.0))

    # 5. Error contract: unphysical emissivity throws DomainError
    @test_throws DomainError Erebus.compute_radiation_htc(T_s, T_a; emissivity=-0.1)
    @test_throws DomainError Erebus.compute_radiation_htc(T_s, T_a; emissivity=1.5)
end
```

---

## 5. Validation Documentation and Citation Grounding

Every physics module in `src/` must have a corresponding validation anchor page in `docs/src/validation/`:

- `docs/src/validation/thermal_conduction.md`: Porous rock mixture models and geometric metric divergence.
- `docs/src/validation/permeability.md`: Porosity-permeability relations and hydrofracture bounds.
- `docs/src/validation/radionuclides.md`: Short-lived radionuclide heating and decay constants.
- `docs/src/validation/fluid_viscosity.md`: Temperature-dependent viscosity and phase change boundaries.
- `docs/src/validation/disk_radiation.md`: Protoplanetary disk evolution and surface radiative boundaries.

Each validation page must list:
1. Full literature reference with verified DOI.
2. The specific equation number or figure from the paper.
3. The exact testset in `test/` verifying the formulation.

---

## 6. Test Quality Linter (`tools/check_test_quality.jl`)

A native Julia AST linter enforces these rules in CI. The linter parses all test files using `Meta.parse`:

- Verifies that testsets contain >= 2 assertions (`min_asserts`).
- Blocks standalone weak assertions (`weak_assert`).
- Flags float comparisons using `==` or `.==` (`float_equality`).
- Compares violation counts against `tools/test_quality_baseline.json` as a one-way ratchet: violation counts cannot increase.

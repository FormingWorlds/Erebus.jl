# Contributing Guide

Contributions to `Erebus.jl` are welcome. This guide outlines development practices, test verification, and code style standards.

---

## Development Workflow

1. **Fork or Branch**:
   Create a dedicated feature branch for your changes:
   ```bash
   git checkout -b feature/my-enhancement
   ```

2. **Environment Setup**:
   Ensure dependencies are instantiated:
   ```bash
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```

3. **Code Style**:
   - Follow the [BlueStyle](https://github.com/invenia/BlueStyle) code formatting convention.
   - Format code using `JuliaFormatter.jl` before submitting:
     ```julia
     using JuliaFormatter
     format(".", BlueStyle())
     ```
   - Continuous Integration automatically validates BlueStyle formatting on every pull request.
   - Use `DocStringExtensions` for docstrings with `$(SIGNATURES)` and `$(FIELDS)`.
   - Prefer type annotations and `StaticArrays` for performance-critical inner loops.


4. **Testing**:
   Run the test suite before submitting:
   ```bash
   julia --project=. -e 'using Pkg; Pkg.test()'
   ```
   All tests across `Geometry`, `Physics`, `Particles`, `Numerics`, `Config`, `Simulation`, and `Integration` must pass.

5. **Pull Requests**:
   Submit your pull request against the `main` branch with a concise description of what changed, why, and how it was verified.

---

## Testing Standards and Quality Gates

Erebus enforces automated quality ratchets across code structure, test quality, configuration schema, and numerical regressions.

### Test Execution Groups

The test suite divides into two tiers using the `EREBUS_TEST_GROUP` environment variable:
- `unit`: Fast unit tests, analytical verifications, and physics checks (default).
  ```bash
  julia --project=. -e 'ENV["EREBUS_TEST_GROUP"] = "unit"; include("test/runtests.jl")'
  ```
- `integration`: Full simulation loops, checkpoint restarts, multi-step coupled runs, and numerical reference runs.
  ```bash
  julia --project=. -e 'ENV["EREBUS_TEST_GROUP"] = "integration"; include("test/runtests.jl")'
  ```
- `all`: Executes both unit and integration suites.
  ```bash
  julia --project=. -e 'ENV["EREBUS_TEST_GROUP"] = "all"; include("test/runtests.jl")'
  ```

### Static Architecture Ratchet

The repository enforces architectural invariants with `tools/check_architecture.jl`:
- Function line span and positional argument count cannot exceed tracked baselines in `tools/architecture_baseline.json`.
- Global mutable collections (such as `Dict` or `Set`) at top level are forbidden.
- Random number generation requires explicit RNG passing; bare `rand()` calls outside authorized sampling routines are rejected.
Verify architecture compliance locally:
```bash
julia --project=. tools/check_architecture.jl --check
```

### Test Quality Ratchet

Test assertions must verify concrete quantitative invariants rather than trivial structural checks:
- No float comparisons with equality (`==` or `.==`). Use `isapprox` or `≈` with explicit physical tolerances.
- No weak assertions such as testing `!== nothing`, `length(x) > 0`, or bare non-negativity against zero when an analytical expectation exists.
- Leaf `@testset` blocks must contain at least 2 assertions.
Verify test quality compliance locally:
```bash
julia --project=. tools/check_test_quality.jl --check
```

### Configuration Schema Ratchet

Every field in `src/config.jl` must exist in `docs/src/reference/config_schema.md` with its type, default value, and physical description:
```bash
julia --project=. tools/check_config_schema.jl --check
```

### Mutation Testing Suite

Key physics and numerical routines in `test/test_mutation.jl` are verified against mutations:
- Stokes-Darcy continuity coupling
- Thermal diffusion direction
- Volatile solubility scaling
- Jeans kinetic escape velocity
- Metal-silicate volatile mass conservation and partition coefficients
- Surface venting drainage budgets
- Radiogenic heating lifetimes
- Marker-to-grid bilinear weight partition of unity
- Silicate melt segregation buoyancy
- Plastic yielding thresholds
Every mutant must fail its invariant check to confirm test discrimination.

### Numerical Reference Runs

Integration tests in `test/test_reference_runs.jl` verify physical solutions against baseline outputs stored in `test/data/reference_runs.json`:
- Reference simulations track total mass, total thermal energy, peak temperature, core radius, and exsolved volatile mass across multiple timesteps.
- Solutions must agree with reference baselines within a relative tolerance of 1.0e-3.

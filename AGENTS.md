# Erebus Development Guidelines

Erebus is a Julia package for 2D thermochemical evolution, hydrothermal circulation, and volatile transport in porous planetesimals.

## Rules Discovery
Testing standards and deep-dive test rules live in [`.github/.claude/rules/erebus-tests.md`](.claude/rules/erebus-tests.md). When opening, creating, or editing any test or physics routine, follow those rules.

## The Iron Law of Test-Driven Development (TDD)
All new features, bug fixes, physical parameterizations, and solver adjustments in Erebus MUST follow the RED-GREEN-REFACTOR cycle:
1. RED. Write the minimal failing test first (asserting physical invariants, analytical limits, or error contracts) and watch it fail.
2. GREEN. Write minimal physical code to make the test pass.
3. REFACTOR. Apply BlueStyle formatting, eliminate duplication, and verify zero performance regression.

## Development Commands

```bash
# Run full test suite
julia --project -e 'using Pkg; Pkg.test()'

# Run a single test file
julia --project -e 'using Test; include("test/test_geometry_radiation.jl")'

# Check code formatting (BlueStyle)
julia -e 'using JuliaFormatter; if !format(".", BlueStyle(); overwrite=false, verbose=true) error("Code style check failed!") end'

# Format code automatically
julia -e 'using JuliaFormatter; format(".", BlueStyle())'

# Build documentation
julia --project=docs docs/make.jl

# Check test quality linter
julia --project tools/check_test_quality.jl --check
```

## Testing Standards
- Tests must verify physical laws rather than duplicating equations from `src/` with random numbers.
- Assert Conservation (mass, energy), Positivity ($T > 0$, $P > 0$, $\phi \in [0, 1]$), Monotonicity, or Analytical Limits with 3-class discrimination guards.
- Enforce input domain contracts by validating inputs and throwing `DomainError` or `ArgumentError` on unphysical parameters.
- Never compare floating-point numbers with `==`. Use `isapprox` or `≈` with explicit tolerances.

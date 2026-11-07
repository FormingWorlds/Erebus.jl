# Erebus Development Guidelines

Erebus is a Julia package for 2D thermochemical evolution, hydrothermal circulation, and volatile transport in porous planetesimals.

## Rules Discovery
Testing standards and deep-dive test rules live in [`.github/.claude/rules/erebus-tests.md`](.github/.claude/rules/erebus-tests.md). When opening, creating, or editing any test or physics routine, follow those rules.

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
julia --project=tools tools/check_test_quality.jl --check
```

## Testing Standards
- Tests must verify physical laws rather than duplicating equations from `src/` with random numbers.
- Assert Conservation (mass, energy), Positivity ($T > 0$, $P > 0$, $\phi \in [0, 1]$), Monotonicity, or Analytical Limits with 3-class discrimination guards.
- Enforce input domain contracts by validating inputs and throwing `DomainError` or `ArgumentError` on unphysical parameters.
- Never compare floating-point numbers with `==`. Use `isapprox` or `≈` with explicit tolerances.

## Code Style and Linting
- Format all Julia code with JuliaFormatter 2.14.5 in BlueStyle (`.JuliaFormatter.toml`, line length 92). The CI job "Code Style - BlueStyle" fails on any file that the formatter changes.
- Keep `julia --project=tools tools/check_architecture.jl --check` passing: no mutable `const` globals outside its allowlist, argument budgets, RNG and coordinate guards, function line spans. Never raise a baseline to make a check pass.
- Keep `julia --project=tools tools/check_test_quality.jl --check` passing: no weak assertions, no float `==`, at least 2 assertions per leaf testset.
- A function that changes one of its arguments has a name that ends in `!`. Compare floats with `isapprox` or `≈` and an explicit tolerance.
- Every new function has a docstring with its parameters (and units), return value and errors, in the style of the file (`# Parameters`, `# Returns`, `# Raises` in `src/`).
- Keep inline comments to 1 to 3 lines; a longer explanation belongs in the docstring.
- Comments, docstrings and docs describe the current code. No change history, dates, issue, PR or commit references.

## Documentation
- The docs are Documenter.jl pages in `docs/src/` (tutorials, howto, explanations, reference, validation). Build them with `julia --project=docs docs/make.jl` and check a change on the built pages.
- Every new physics option or default gets an explanation of its equations, assumptions and limits.
- Formatting: one sentence per line in new or edited prose; inline math in `$...$`, display math in a `math` block; no list item that starts with math (`tools/check_docs_formatting.jl` checks the source and the built HTML); no em or en dashes.
- References: cite in the text as (Author, Year), (Author and Author, Year) or (Author et al., Year). Every cited work has an entry in `docs/src/reference/bibliography.md` with authors, year, title, journal and a `https://doi.org/` link where one exists; `tools/verify_references.jl` checks that every DOI resolves. Cite a work only for a claim that it makes, and only after you have checked the source.
- Numbers: a number in the docs that a test asserts is registered in `docs/src/reference/benchmarks/number_map.toml` and checked by `tools/check_doc_numbers.jl`. Validation tables come from the records in `benchmarks/records/` through their export tools, not from hand edits.

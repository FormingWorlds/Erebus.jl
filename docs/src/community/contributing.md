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

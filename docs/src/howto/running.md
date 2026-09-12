# Running Simulations

This guide covers running `Erebus.jl` simulations from the command-line interface (CLI) and from within the Julia REPL.

---

## Command-Line Interface (CLI)

The package provides a built-in CLI entry point via `run_simulation()`.

### Running with a TOML Configuration File

To execute a simulation from a configuration file:

```bash
julia --project=. -e 'using Erebus; run_simulation("configs/default.toml")'
```

Alternatively, invoke the CLI script directly with command-line arguments:

```bash
julia --project=. launch.jl configs/default.toml
```

### Overriding the Output Directory

To run with a configuration file while redirecting the output to a custom directory:

```bash
julia --project=. -e 'using Erebus; simulation_loop("configs/default.toml", output_path="custom_run_01")'
```

Via CLI flags:

```bash
julia --project=. launch.jl configs/default.toml --output_path custom_run_01
```

### Displaying Solver Timers

To print detailed performance diagnostics at the end of the simulation:

```bash
julia --project=. launch.jl configs/default.toml --show_timer true
```

---

## Interactive Julia REPL Workflow

You can load and execute simulations interactively within Julia:

```julia
using Erebus

# Load and validate configuration
cfg = load_config("configs/default.toml")

# Inspect specific parameters
println("Timesteps: ", cfg.time.n_steps)
println("Solid compressibility: ", cfg.poroelasticity.betasolid)

# Run the simulation loop
simulation_loop(cfg)
```

---

## Output and Log Monitoring

During simulation execution, Erebus emits structured log messages:
- Initial grid dimensions, domain sizes, and physical parameters.
- Thermochemical iteration progress (`titer`) and hydromechanical solver iterations (`iplast`).
- Timestep duration (`dt`), CFL displacement bounds, and convergence metrics (`pferrcur`, `DMPmax`).
- Runtime duration per timestep and total elapsed simulation time (in Ma).

Binary checkpoints (`.jld2`) are saved to the designated `output_dir` at intervals determined by `savematstep`.

---

## Performance Profiling and Benchmarking

To benchmark execution runtime and memory allocation across repeated runs:

```bash
julia --project tools/profile_simulation.jl configs/test_quick.toml
```

The profiling harness:
1. Executes a warm-up run to compile all methods.
2. Executes repeated benchmark runs (default: 3 iterations).
3. Reports minimum and mean wall-clock execution time, plus mean memory allocation and mean garbage collection duration.
4. Saves structured JSON metrics to `output_files/profiling_baseline.json`.

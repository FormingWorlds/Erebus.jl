# Parameter Space Exploration

`Erebus.jl` allows systematic parameter space explorations by programmatically modifying parameters, writing configuration overlays, or constructing custom `SimulationConfig` instances.

---

## Partial Configuration Overlays

`load_config` supports partial TOML strings. Omitted sections and keys automatically inherit standard default values.

```julia
using Erebus

# Define only the parameters you wish to explore
custom_toml = """
[poroelasticity]
betasolid = 5.0e-11
betafluid = 8.0e-10

[thermodynamics]
hr_fe = true
ratio_al = 1.0e-4

[time]
n_steps = 20

[output]
output_dir = "sweep_run_high_heating"
"""

cfg = load_config(custom_toml)
simulation_loop(cfg)
```

---

## Parameter Sweeps

To execute a series of parameter variations across a grid of compressibility or heating values:

```julia
using Erebus

# Define parameter grid
compressibilities = [1.0e-11, 2.5e-11, 5.0e-11, 1.0e-10]

for (idx, beta_s) in enumerate(compressibilities)
    run_dir = "run_beta_$(idx)"
    
    # Construct modified configuration
    cfg = SimulationConfig(
        poroelasticity = PoroelasticConfig(
            betasolid = beta_s,
            betafluid = 4.0e-10
        ),
        time = TimeConfig(n_steps = 5),
        output = OutputConfig(output_dir = run_dir)
    )
    
    # Validate before running
    validate_config(cfg)
    
    println("Launching run $idx: betasolid = $beta_s -> $run_dir")
    simulation_loop(cfg)
end
```

---

## Saving and Loading Generated Configurations

Configurations can be serialized to disk for archival or cluster deployment:

```julia
using Erebus

cfg = default_config()
# Modify parameters
custom_cfg = SimulationConfig(
    poroelasticity = PoroelasticConfig(betasolid = 3.0e-11, betafluid = 6.0e-10),
    time = TimeConfig(n_steps = 50),
    output = OutputConfig(output_dir = "cluster_job_01")
)

# Save to TOML file
save_config("job_config.toml", custom_cfg)

# Reload and run
reloaded_cfg = load_config("job_config.toml")
simulation_loop(reloaded_cfg)
```

---

## Ensemble Sweeps via Command Line

Run multiple simulations with varied parameter combinations (an ensemble). The script samples parameters using Latin Hypercube sampling, where each parameter interval is divided into equal segments and sampled once.

To run parameter sweeps in parallel across worker processes:

```bash
julia -p 4 --project tools/run_ensemble.jl configs/test_ensemble_sweep.toml --workers 4
```

The script parses parameter ranges from `configs/test_ensemble_sweep.toml`, generates sample configurations, and executes ensemble simulations in parallel. Output files are saved to the directories specified in the sweep configuration file.

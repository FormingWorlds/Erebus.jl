# Quickstart Tutorial

In this tutorial, you will execute a fast 2-step simulation of a porous planetesimal using the pre-configured `configs/test_quick.toml` input file.

---

## Step 1: Examine the Test Configuration

View the contents of `configs/test_quick.toml`:

```toml
[time]
n_steps = 2
dt_initial = 1.0e+11
dt_longest = 1.0e+11

[output]
output_dir = "output_test"
savematstep = 2
visstep = 1
```

This configuration runs 2 timesteps of $1.0 \times 10^{11}\text{ s}$ each ($\approx 3,170\text{ years}$ per step), using the standard $15 \times 15$ staggered grid with $3,136$ markers.

---

## Step 2: Run the Simulation

Launch the simulation from the Julia REPL:

```julia
using Erebus

# Load quick test configuration
cfg = load_config("configs/test_quick.toml")

# Run the simulation loop
simulation_loop(cfg)
```

The solver outputs information about each iteration:

```text
[ Info: Simulation layout Nx = 15 Ny = 15 xsize = 14000.0 ysize = 14000.0 marknum = 3136
[ Info: ********** begin timestep 1 - dt = 1.0e11 s **********
[ Info: thermochemical iter 1 - hydromechanical iter 1
[ Info: starting hydro-mechanical solver 1-1
[ Info: finished hydro-mechanical solver 1-1
[ Info: timestep 1 computed in 6 milliseconds
[ Info: total time = 2.253 Ma
[ Info: ********** begin timestep 2 - dt = 1.0e11 s **********
...
[ Info: timestep 2 computed in 15 milliseconds
[ Info: total time = 2.256 Ma
```

---

## Step 3: Inspect the Output

Check that the checkpoint was generated in `output_test`:

```julia
using JLD2

# Open the checkpoint
data = jldopen("output_test/step_2.jld2", "r")
println("Saved keys: ", keys(data))
println("Max temperature: ", maximum(data["tk"]), " K")
println("Solid velocity vx range: ", extrema(data["vx"]), " m/s")
close(data)
```

You have successfully completed a full simulation run with `Erebus.jl`.

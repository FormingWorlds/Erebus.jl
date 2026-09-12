#!/usr/bin/env julia
# tools/profile_simulation.jl
# In-tree profiling harness for Erebus.jl
# Measures runtime and memory allocation baselines for simulation execution.

using Dates
using Erebus
using JSON
using Printf

function run_profile(config_path::String; n_runs::Int=3)
    isfile(config_path) || error("Config file not found: $config_path")
    cfg = load_config(config_path)

    # Warm-up run to compile all Julia methods
    println("Performing warm-up run...")
    mktempdir() do tmpdir
        Erebus.run_simulation(cfg; output_path=tmpdir)
    end

    # Profile runs
    println("Executing $(n_runs) profiling runs on $(config_path)...")
    times = Float64[]
    allocs = Int64[]
    gctimes = Float64[]

    for i in 1:n_runs
        mktempdir() do tmpdir
            res = @timed Erebus.run_simulation(cfg; output_path=tmpdir)
            push!(times, res.time)
            push!(allocs, res.bytes)
            push!(gctimes, res.gctime)
            @printf(
                "  Run %d: %.3f s, %.2f MB allocated (gc: %.3f s)\n",
                i,
                res.time,
                res.bytes / (1024 * 1024),
                res.gctime
            )
        end
    end

    mean_time = sum(times) / n_runs
    min_time = minimum(times)
    mean_alloc = sum(allocs) / n_runs
    mean_gc = sum(gctimes) / n_runs

    profile_data = Dict(
        "timestamp" => string(now()),
        "config" => config_path,
        "grid_size" => [cfg.grid.Nx, cfg.grid.Ny],
        "n_steps" => cfg.time.n_steps,
        "n_runs" => n_runs,
        "min_wall_time_s" => min_time,
        "mean_wall_time_s" => mean_time,
        "mean_alloc_bytes" => mean_alloc,
        "mean_alloc_mb" => mean_alloc / (1024 * 1024),
        "mean_gc_time_s" => mean_gc,
        "all_wall_times_s" => times,
        "all_alloc_bytes" => allocs,
    )

    return profile_data
end

function main()
    config_file = length(ARGS) >= 1 ? ARGS[1] : "configs/test_quick.toml"
    output_dir = "output_files"
    mkpath(output_dir)
    output_file = joinpath(output_dir, "profiling_baseline.json")

    println("==========================================================")
    println(" Erebus.jl Simulation Profiling Harness")
    println("==========================================================")
    println("Config: $config_file")

    data = run_profile(config_file; n_runs=3)

    open(output_file, "w") do f
        JSON.print(f, data, 2)
    end

    println("----------------------------------------------------------")
    @printf("Min wall time:   %.3f s\n", data["min_wall_time_s"])
    @printf("Mean wall time:  %.3f s\n", data["mean_wall_time_s"])
    @printf("Mean allocated:  %.2f MB\n", data["mean_alloc_mb"])
    @printf("Mean GC time:    %.3f s\n", data["mean_gc_time_s"])
    println("Baseline metrics saved to: $output_file")
    println("==========================================================")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

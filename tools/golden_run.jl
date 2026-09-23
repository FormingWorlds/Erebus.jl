#!/usr/bin/env julia
# tools/golden_run.jl
# Golden-run execution harness for Erebus.jl
# Runs 5 steps deterministically and writes complete state to JLD2.

using Pkg
const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(ROOT_DIR; io=devnull)

@eval using Dates
@eval using Erebus
@eval using JLD2
@eval using LinearAlgebra
@eval using Random
@eval using TOML

function run_golden(config_path::String, out_jld2_path::String)
    isfile(config_path) || error("Config file not found: $config_path")

    # Enforce deterministic single-thread execution
    LinearAlgebra.BLAS.set_num_threads(1)

    cfg_base = load_config(config_path)
    overrides = Dict{String,Any}(
        "time.n_steps" => 5,
        "solver.p2m_mode" => :tiled,
        "solver.seed" => 42,
    )
    cfg = Erebus.override_config(cfg_base, overrides)

    # Seed RNGs deterministically for golden run reproducibility
    Random.seed!(Erebus.rgen, cfg.solver.seed)
    Random.seed!(cfg.solver.seed)

    mktempdir() do tmpdir
        res = Erebus.simulation_loop(cfg; output_path=tmpdir)
        jldopen(out_jld2_path, "w") do f
            # 1. Marker arrays
            for (k, v) in pairs(res.markers)
                if v !== nothing
                    f[string(k)] = v
                end
            end

            # 2. Grid fields
            for (k, v) in pairs(res.grids)
                if v !== nothing
                    f[string(k)] = v
                end
            end

            # 3. Atmosphere state
            if res.atm !== nothing
                f["atm"] = Dict{String,Any}(
                    "P_surf" => res.atm.P_surf,
                    "T_surf_eq" => res.atm.T_surf_eq,
                    "M_atm" => res.atm.M_atm,
                    "M_escaped" => res.atm.M_escaped,
                    "tau_LW" => res.atm.tau_LW,
                    "M_env_bound" => res.atm.M_env_bound,
                    "F_net_rad" => res.atm.F_net_rad,
                    "h_rad_eff" => res.atm.h_rad_eff,
                )
            else
                f["atm"] = Dict{String,Any}()
            end

            # 4. Time accumulators
            f["timesum"] = res.timesum
            f["dt"] = res.dt
            f["timestep"] = res.timestep
        end
    end
    println("Golden run written successfully to: $out_jld2_path")
    return out_jld2_path
end

function main()
    if length(ARGS) < 2
        println(stderr, "Usage: julia tools/golden_run.jl <config.toml> <out.jld2>")
        exit(1)
    end
    config_path = ARGS[1]
    out_jld2_path = ARGS[2]
    run_golden(config_path, out_jld2_path)
    exit(0)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

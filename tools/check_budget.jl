#!/usr/bin/env julia
# tools/check_budget.jl
# Performance budget verification harness for Erebus.jl.
# Measures @allocated and wall-clock over steps 3-5 after 2 warm-up steps.

using Pkg
const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(ROOT_DIR; io=devnull)

@eval using Dates
@eval using Erebus
@eval using JSON
@eval using Printf
@eval using TOML

const BASELINE_PATH = joinpath(@__DIR__, "budget_baseline.json")

const SHIPPED_CONFIGS = [
    "configs/test_quick.toml",
    "configs/hydrothermal_benchmark.toml",
    "configs/core_formation_benchmark.toml",
    "configs/lunar_growth_tutorial.toml",
    "configs/magma_ocean_cooling_turb_on_32.toml",
]

function measure_config_budget(config_path::String)
    full_path = isabspath(config_path) ? config_path : joinpath(ROOT_DIR, config_path)
    isfile(full_path) || error("Config file not found: $full_path")

    cfg_base = load_config(full_path)

    mktempdir() do tmpdir
        # 1. Warm-up steps 1 to 2
        cfg_warmup = Erebus.override_config(
            cfg_base,
            Dict{String,Any}(
                "time.n_steps" => 2,
                "output.output_dir" => tmpdir,
                "output.savematstep" => 2,
                "solver.p2m_mode" => :tiled,
                "solver.seed" => 42,
            ),
        )
        Erebus.simulation_loop(cfg_warmup; output_path=tmpdir)

        ckpt2 = joinpath(tmpdir, "output_00002.jld2")
        isfile(ckpt2) || error("Warm-up checkpoint not found: $ckpt2")

        # 2. Measured steps 3 to 5
        cfg_eval = Erebus.override_config(
            cfg_base,
            Dict{String,Any}(
                "time.n_steps" => 5,
                "time.start_step" => 3,
                "output.output_dir" => tmpdir,
                "output.restart_from" => ckpt2,
                "output.savematstep" => 1000,
                "solver.p2m_mode" => :tiled,
                "solver.seed" => 42,
            ),
        )

        # Force GC before measurement
        GC.gc()

        wall_seconds = 0.0
        allocated_bytes = @allocated begin
            t0 = time_ns()
            Erebus.simulation_loop(cfg_eval; output_path=tmpdir)
            wall_seconds = (time_ns() - t0) / 1.0e9
        end

        return (; wall_seconds, allocated_bytes)
    end
end

function write_baseline(baseline_data::Dict{String,Any})
    open(BASELINE_PATH, "w") do io
        JSON.print(io, baseline_data, 4)
        return println(io)
    end
    return println("Saved budget baseline to $BASELINE_PATH")
end

function read_baseline()
    if !isfile(BASELINE_PATH)
        return Dict{String,Any}()
    end
    return JSON.parsefile(BASELINE_PATH)
end

function main()
    mode = length(ARGS) >= 1 ? ARGS[1] : "--check"

    if mode == "--baseline"
        data = Dict{String,Any}()
        for cfg_rel in SHIPPED_CONFIGS
            println("Profiling budget baseline for $cfg_rel...")
            res = measure_config_budget(cfg_rel)
            @printf(
                "  %s: %.3f s, %d bytes\n", cfg_rel, res.wall_seconds, res.allocated_bytes
            )
            data[cfg_rel] = Dict(
                "wall_seconds" => res.wall_seconds, "allocated_bytes" => res.allocated_bytes
            )
        end
        write_baseline(data)
        exit(0)
    elseif mode == "--check"
        baseline = read_baseline()
        if isempty(baseline)
            println(
                stderr,
                "No baseline file found at $BASELINE_PATH. Run with --baseline first.",
            )
            exit(2)
        end

        has_regression = false
        println("=== Performance Budget Verification ===")
        for cfg_rel in SHIPPED_CONFIGS
            haskey(baseline, cfg_rel) || continue
            base = baseline[cfg_rel]
            res = measure_config_budget(cfg_rel)

            time_ratio = res.wall_seconds / max(1e-6, base["wall_seconds"])
            alloc_ratio = res.allocated_bytes / max(1, base["allocated_bytes"])

            # If jitter caused threshold breach, retry once and take best sample
            if time_ratio > 1.10 || alloc_ratio > 1.10
                res2 = measure_config_budget(cfg_rel)
                res = (;
                    wall_seconds=min(res.wall_seconds, res2.wall_seconds),
                    allocated_bytes=min(res.allocated_bytes, res2.allocated_bytes),
                )
                time_ratio = res.wall_seconds / max(1e-6, base["wall_seconds"])
                alloc_ratio = res.allocated_bytes / max(1, base["allocated_bytes"])
            end

            status_time = time_ratio <= 1.10 ? "PASS" : "FAIL (REGRESSION)"
            status_alloc = alloc_ratio <= 1.10 ? "PASS" : "FAIL (REGRESSION)"

            @printf(
                "%s:\n  Time: %.3fs vs base %.3fs (ratio %.2f) [%s]\n  Alloc: %d B vs base %d B (ratio %.2f) [%s]\n",
                cfg_rel,
                res.wall_seconds,
                base["wall_seconds"],
                time_ratio,
                status_time,
                res.allocated_bytes,
                base["allocated_bytes"],
                alloc_ratio,
                status_alloc,
            )

            if time_ratio > 1.10 || alloc_ratio > 1.10
                has_regression = true
            end
        end

        if has_regression
            println(
                stderr, "Performance budget regression detected (>10% threshold exceeded)."
            )
            exit(1)
        else
            println("All configurations satisfied performance budget.")
            exit(0)
        end
    else
        println(stderr, "Usage: julia tools/check_budget.jl [--baseline | --check]")
        exit(1)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

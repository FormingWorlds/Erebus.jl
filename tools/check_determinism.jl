#!/usr/bin/env julia
# tools/check_determinism.jl
# CI determinism and idempotence verification harness for Erebus.jl.
# Runs test_quick.toml for 5 steps twice and asserts bitwise equality.

using Pkg
const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(ROOT_DIR; io=devnull)

@eval using Erebus
@eval using JLD2
@eval using LinearAlgebra
@eval include(joinpath(normpath(joinpath(@__DIR__, "..")), "test", "golden_helpers.jl"))
@eval using .GoldenHelpers
@eval include(joinpath(normpath(joinpath(@__DIR__, "..")), "tools", "golden_run.jl"))

function check_determinism()
    config_path = joinpath(ROOT_DIR, "configs", "test_quick.toml")
    isfile(config_path) || error("Config file not found: $config_path")

    mktempdir() do tmpdir
        out1 = joinpath(tmpdir, "run1.jld2")
        out2 = joinpath(tmpdir, "run2.jld2")

        println("Executing Determinism Run 1 (5 steps)...")
        run_golden(config_path, out1)

        println("Executing Determinism Run 2 (5 steps)...")
        run_golden(config_path, out2)

        println("Comparing states bitwise...")
        state1 = JLD2.load(out1)
        state2 = JLD2.load(out2)

        try
            compare_golden(state1, state2)
            println(
                "Bitwise determinism check PASSED: both runs produced identical fields across all markers and grids.",
            )
            return true
        catch err
            println(stderr, "Bitwise determinism check FAILED:")
            showerror(stderr, err)
            println(stderr)
            return false
        end
    end
end

function main()
    success = check_determinism()
    return exit(success ? 0 : 1)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

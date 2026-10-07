#!/usr/bin/env julia
# tools/run_ensemble.jl
# Command-line driver for running Erebus ensemble parameter sweeps.

using ArgParse
using Erebus
using TOML

function parse_cli()
    s = ArgParseSettings(; description="Run Erebus ensemble simulation sweeps.")
    @add_arg_table! s begin
        "spec_file"
        help = "Path to sweep TOML specification file"
        required = true
        "--workers", "-w"
        help = "Number of worker threads / instances"
        arg_type = Int
        default = 1
    end
    return parse_args(s)
end

function main()
    args = parse_cli()
    spec_path = args["spec_file"]
    isfile(spec_path) || error("Sweep spec file not found: $spec_path")

    spec = load_ensemble_config(spec_path)
    catalog = run_ensemble(spec; max_workers=args["workers"])
    println("Completed ensemble sweep with $(length(catalog)) members.")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

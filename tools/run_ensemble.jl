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

    raw = TOML.parsefile(spec_path)
    haskey(raw, "base") || error("Missing [base] section in $spec_path")
    haskey(raw["base"], "config") || error("Missing 'config' in [base] section")

    base_cfg_path = raw["base"]["config"]
    isfile(base_cfg_path) || error("Base configuration file not found: $base_cfg_path")
    base_cfg = load_config(base_cfg_path)

    ens_sec = get(raw, "ensemble", Dict{String,Any}())
    out_dir = get(ens_sec, "output_dir", "ensemble_output")
    method = Symbol(get(ens_sec, "method", "lhs"))
    num_samples = get(ens_sec, "num_samples", 10)
    seed = get(ens_sec, "seed", 42)
    params = get(raw, "parameters", Dict{String,Any}())

    spec = EnsembleSweepSpec(
        base_cfg;
        output_dir=out_dir,
        sampling_method=method,
        num_samples=num_samples,
        parameters=params,
        seed=seed,
    )

    catalog = run_ensemble(spec; max_workers=args["workers"])
    println("Completed ensemble sweep with $(length(catalog)) members.")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

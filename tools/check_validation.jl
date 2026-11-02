#!/usr/bin/env julia
# Re-run the spherical conduction validation benchmark and verify reproducibility against the committed TOML record.

using Erebus
using TOML

const DEFAULT_SPHERICAL_RECORD = normpath(
    joinpath(@__DIR__, "..", "benchmarks", "records", "spherical_conduction.toml")
)

"""
    run_validation_check(record_path::String=DEFAULT_SPHERICAL_RECORD; bench=nothing) -> Bool

Execute the spherical conduction validation benchmark and verify schema and metric reproducibility
against `record_path`. Requires `SphericalConductionValidation` (`benchmarks/validation_spherical_conduction.jl`)
to be loaded in `Main` before calling (loaded automatically when `tools/check_validation.jl` is run as a script).

Parameters
----------
record_path : String, optional
    Path to the TOML validation record file.
bench : Any, optional
    Pre-evaluated benchmark tuple, or `nothing` to run `SphericalConductionValidation.evaluate_spherical_conduction_benchmark()`.

Returns
-------
Bool
    `true` when all validation and reproducibility checks pass.

Raises
------
DomainError
    If `isempty(record_path)`.
"""
function run_validation_check(
    record_path::String=DEFAULT_SPHERICAL_RECORD; bench=nothing
)::Bool
    isempty(record_path) && throw(DomainError(record_path, "Record path must be non-empty"))
    if !isfile(record_path)
        println(stderr, "Validation record file not found: $record_path")
        return false
    end
    raw = TOML.parsefile(record_path)
    ok_schema, schema_errs = SphericalConductionValidation.validate_record_dict(
        raw; require_doc_pages=false, check_doc_sync=false
    )
    if !ok_schema
        for e in schema_errs
            println(stderr, "[schema] $e")
        end
        return false
    end

    bench_res = if bench === nothing
        SphericalConductionValidation.evaluate_spherical_conduction_benchmark()
    else
        bench
    end
    ok_repro, repro_errs = SphericalConductionValidation.check_validation_reproducibility(
        bench_res, record_path
    )
    if !ok_repro
        for e in repro_errs
            println(stderr, "[reproducibility] $e")
        end
        return false
    end
    println("check_validation: all metrics reproduced and verified against $record_path")
    return true
end

"""
    run_validation_cli(args::Vector{String}=ARGS; bench=nothing) -> Int

Run the command-line validation check for the spherical conduction record. Requires
`SphericalConductionValidation` (`benchmarks/validation_spherical_conduction.jl`) to be loaded
in `Main` before calling.

Parameters
----------
args : Vector{String}, optional
    Command-line arguments (`--check` and optional `record_path`).
bench : Any, optional
    Pre-evaluated benchmark tuple, or `nothing` to evaluate on demand.

Returns
-------
Int
    Exit code (`0` on pass, `1` on validation failure, `2` on invalid CLI flag).
"""
function run_validation_cli(args::Vector{String}=ARGS; bench=nothing)::Int
    record_path = DEFAULT_SPHERICAL_RECORD
    for arg in args
        if arg == "--check"
            continue
        elseif startswith(arg, "--")
            println(
                stderr,
                "Usage: julia --project=. tools/check_validation.jl [--check] [record_path]",
            )
            return 2
        else
            record_path = arg
        end
    end
    ok = run_validation_check(record_path; bench=bench)
    return ok ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_conduction.jl"))
    exit(Base.invokelatest(run_validation_cli, ARGS))
end

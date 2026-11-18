#!/usr/bin/env julia
# Re-run the validation benchmarks of VALIDATION_CASES and verify reproducibility against the committed TOML records.

using Erebus
using TOML

const VALIDATION_RECORDS_DIR = normpath(joinpath(@__DIR__, "..", "benchmarks", "records"))
const VALIDATION_BENCHMARKS_DIR = normpath(joinpath(@__DIR__, "..", "benchmarks"))
const DEFAULT_SPHERICAL_RECORD = joinpath(
    VALIDATION_RECORDS_DIR, "spherical_conduction.toml"
)
const DEFAULT_DARCY_RECORD = joinpath(VALIDATION_RECORDS_DIR, "spherical_darcy_flow.toml")
const DEFAULT_COMPACTION_RECORD = joinpath(
    VALIDATION_RECORDS_DIR, "spherical_compaction.toml"
)

# Ordered table of the validation records: the no-argument run checks them in this order.
const VALIDATION_CASES = (
    (;
        id="spherical_conduction",
        module_name=:SphericalConductionValidation,
        file="validation_spherical_conduction.jl",
        record=DEFAULT_SPHERICAL_RECORD,
        evaluate=:evaluate_spherical_conduction_benchmark,
        bench_keyword=:bench,
    ),
    (;
        id="spherical_darcy_flow",
        module_name=:SphericalDarcyFlowValidation,
        file="validation_spherical_darcy_flow.jl",
        record=DEFAULT_DARCY_RECORD,
        evaluate=:evaluate_spherical_darcy_flow_benchmark,
        bench_keyword=:darcy_bench,
    ),
    (;
        id="spherical_compaction",
        module_name=:SphericalCompactionValidation,
        file="validation_spherical_compaction.jl",
        record=DEFAULT_COMPACTION_RECORD,
        evaluate=:evaluate_spherical_compaction_benchmark,
        bench_keyword=:compaction_bench,
    ),
)

"""
    validation_case(id::AbstractString) -> Union{NamedTuple,Nothing}

Return the entry of `VALIDATION_CASES` with the record id `id`.

Parameters
----------
id : AbstractString
    Record id, the `metadata.id` of a validation record.

Returns
-------
Union{NamedTuple,Nothing}
    The table entry, or `nothing` when no entry has this id.
"""
function validation_case(id::AbstractString)
    idx = findfirst(c -> c.id == id, VALIDATION_CASES)
    return idx === nothing ? nothing : VALIDATION_CASES[idx]
end

"""
    validation_module(case::NamedTuple; host::Module=@__MODULE__) -> Module

Return the benchmark module of a `VALIDATION_CASES` entry from `host`. A module not yet
defined in `host` is loaded there from `benchmarks/` first, so the check runs from a script
and from a test that loaded only some of the modules.

Parameters
----------
case : NamedTuple
    Entry of `VALIDATION_CASES`.
host : Module, optional
    Module that holds the benchmark modules (default: the module that includes this file).

Returns
-------
Module
    The benchmark module.
"""
function validation_module(case::NamedTuple; host::Module=@__MODULE__)::Module
    if !Base.invokelatest(isdefined, host, case.module_name)
        Base.include(host, joinpath(VALIDATION_BENCHMARKS_DIR, case.file))
    end
    return Base.invokelatest(getfield, host, case.module_name)
end

"""
    _error_report(err, backtrace) -> String

Return the message of `err` followed by its stack trace.

Parameters
----------
err : Exception
    Caught exception.
backtrace : Vector
    Backtrace of the catch, from `catch_backtrace()`.

Returns
-------
String
    Message and stack trace as `showerror` prints them.
"""
_error_report(err, backtrace) = sprint(showerror, err, backtrace)

"""
    run_record_check(case::NamedTuple, record_path::String; bench=nothing, raw=nothing) -> Bool

Run the benchmark of a `VALIDATION_CASES` entry (or use `bench`) and verify the schema and
the reproducibility of the record `record_path`. A record that is not valid TOML and an
exception in the benchmark or in the reproducibility check are reported, the latter with its
stack trace, as a failed check; an `InterruptException` stops the check.

Parameters
----------
case : NamedTuple
    Entry of `VALIDATION_CASES`.
record_path : String
    Path to the TOML validation record file.
bench : Any, optional
    Pre-evaluated benchmark tuple, or `nothing` to run the `evaluate` function of the case.
raw : Union{Dict{String,Any},Nothing}, optional
    The record already parsed from `record_path`, or `nothing` to parse it here.

Returns
-------
Bool
    `true` when all validation and reproducibility checks pass.

Raises
------
DomainError
    If `isempty(record_path)`.
InterruptException
    If the run is interrupted.
"""
function run_record_check(
    case::NamedTuple,
    record_path::String;
    bench=nothing,
    raw::Union{Dict{String,Any},Nothing}=nothing,
)::Bool
    isempty(record_path) && throw(DomainError(record_path, "Record path must be non-empty"))
    if !isfile(record_path)
        println(stderr, "Validation record file not found: $record_path")
        return false
    end
    if raw === nothing
        raw = try
            TOML.parsefile(record_path)
        catch err
            err isa InterruptException && rethrow()
            println(stderr, "[schema] Cannot parse $record_path: $(sprint(showerror, err))")
            return false
        end
    end
    mod = validation_module(case)
    # A module loaded by `validation_module` is newer than this method: read its bindings
    # in the latest world.
    binding(name) = Base.invokelatest(getfield, mod, name)
    ok_schema, schema_errs = Base.invokelatest(
        binding(:validate_record_dict), raw; require_doc_pages=false, check_doc_sync=false
    )
    if !ok_schema
        for e in schema_errs
            println(stderr, "[schema] $e")
        end
        return false
    end
    ok_repro, repro_errs = try
        bench_res = bench === nothing ? Base.invokelatest(binding(case.evaluate)) : bench
        Base.invokelatest(
            binding(:check_validation_reproducibility), bench_res, record_path
        )
    catch err
        err isa InterruptException && rethrow()
        (false, ["$(case.id) check stopped: $(_error_report(err, catch_backtrace()))"])
    end
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
    run_validation_check(record_path::String=DEFAULT_SPHERICAL_RECORD; bench=nothing) -> Bool

Execute the spherical conduction validation benchmark and verify schema and metric
reproducibility against `record_path` (see `run_record_check`).

Parameters
----------
record_path : String, optional
    Path to the TOML validation record file.
bench : Any, optional
    Pre-evaluated benchmark tuple, or `nothing` to run
    `SphericalConductionValidation.evaluate_spherical_conduction_benchmark()`.

Returns
-------
Bool
    `true` when all validation and reproducibility checks pass.

Raises
------
DomainError
    If `isempty(record_path)`.
InterruptException
    If the run is interrupted.
"""
function run_validation_check(
    record_path::String=DEFAULT_SPHERICAL_RECORD; bench=nothing
)::Bool
    return run_record_check(
        validation_case("spherical_conduction"), record_path; bench=bench
    )
end

"""
    run_darcy_validation_check(record_path::String=DEFAULT_DARCY_RECORD; bench=nothing) -> Bool

Execute the spherical Darcy flow validation benchmark and verify schema and metric
reproducibility against `record_path` (see `run_record_check`).

Parameters
----------
record_path : String, optional
    Path to the TOML validation record file.
bench : Any, optional
    Pre-evaluated benchmark tuple, or `nothing` to run
    `SphericalDarcyFlowValidation.evaluate_spherical_darcy_flow_benchmark()`.

Returns
-------
Bool
    `true` when all validation and reproducibility checks pass.

Raises
------
DomainError
    If `isempty(record_path)`.
InterruptException
    If the run is interrupted.
"""
function run_darcy_validation_check(
    record_path::String=DEFAULT_DARCY_RECORD; bench=nothing
)::Bool
    return run_record_check(
        validation_case("spherical_darcy_flow"), record_path; bench=bench
    )
end

"""
    run_compaction_validation_check(record_path::String=DEFAULT_COMPACTION_RECORD; bench=nothing) -> Bool

Execute the spherical compaction validation benchmark and verify schema and metric
reproducibility against `record_path` (see `run_record_check`).

Parameters
----------
record_path : String, optional
    Path to the TOML validation record file.
bench : Any, optional
    Pre-evaluated benchmark tuple, or `nothing` to run
    `SphericalCompactionValidation.evaluate_spherical_compaction_benchmark()`.

Returns
-------
Bool
    `true` when all validation and reproducibility checks pass.

Raises
------
DomainError
    If `isempty(record_path)`.
InterruptException
    If the run is interrupted.
"""
function run_compaction_validation_check(
    record_path::String=DEFAULT_COMPACTION_RECORD; bench=nothing
)::Bool
    return run_record_check(
        validation_case("spherical_compaction"), record_path; bench=bench
    )
end

"""
    run_validation_cli(args::Vector{String}=ARGS; bench=nothing, darcy_bench=nothing, compaction_bench=nothing) -> Int

Run the command-line validation check. Without a record path it checks every record of
`VALIDATION_CASES` in order and collects garbage after each case. With a record path it
checks that record only, with the entry of `VALIDATION_CASES` whose id is the
`metadata.id` of the record; a missing or unreadable file and an id outside the table are
reported as errors. An argument other than `--check` that starts with `-`, and a second
record path, are usage errors.

Parameters
----------
args : Vector{String}, optional
    Command-line arguments (`--check` and at most one `record_path`).
bench : Any, optional
    Pre-evaluated conduction benchmark tuple, or `nothing` to evaluate on demand.
darcy_bench : Any, optional
    Pre-evaluated Darcy flow benchmark tuple, or `nothing` to evaluate on demand.
compaction_bench : Any, optional
    Pre-evaluated compaction benchmark tuple, or `nothing` to evaluate on demand.

Returns
-------
Int
    Exit code (`0` on pass, `1` on validation failure, `2` on invalid command line).

Raises
------
InterruptException
    If the run is interrupted.
"""
function run_validation_cli(
    args::Vector{String}=ARGS; bench=nothing, darcy_bench=nothing, compaction_bench=nothing
)::Int
    benches = (; bench, darcy_bench, compaction_bench)
    record_path = nothing
    for arg in args
        arg == "--check" && continue
        if startswith(arg, "-") || record_path !== nothing
            println(
                stderr,
                "Usage: julia --project=. tools/check_validation.jl [--check] [record_path]",
            )
            return 2
        end
        record_path = arg
    end
    if record_path === nothing
        all_ok = true
        for case in VALIDATION_CASES
            ok = run_record_check(
                case, case.record; bench=getfield(benches, case.bench_keyword)
            )
            all_ok &= ok
            GC.gc()
        end
        return all_ok ? 0 : 1
    end
    if !isfile(record_path)
        println(stderr, "Validation record file not found: $record_path")
        return 1
    end
    raw = try
        TOML.parsefile(record_path)
    catch err
        err isa InterruptException && rethrow()
        println(stderr, "[schema] Cannot parse $record_path: $(sprint(showerror, err))")
        return 1
    end
    meta = get(raw, "metadata", Dict{String,Any}())
    rec_id = meta isa AbstractDict ? string(get(meta, "id", "")) : ""
    case = validation_case(rec_id)
    if case === nothing
        println(stderr, "[schema] Unknown record id '$rec_id' in $record_path")
        return 1
    end
    ok = run_record_check(
        case, record_path; bench=getfield(benches, case.bench_keyword), raw=raw
    )
    return ok ? 0 : 1
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(Base.invokelatest(run_validation_cli, ARGS))
end

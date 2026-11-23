#!/usr/bin/env julia
# Export the porous shell validation record TOML from live benchmark measurements.

module PorousShellExport

using Dates
using Printf
using TOML

const DEFAULT_EXPORT_RECORD_PATH = normpath(
    joinpath(@__DIR__, "..", "benchmarks", "records", "porous_shell.toml")
)

"""
    last_src_commit(repo_root::String=normpath(joinpath(@__DIR__, ".."))) -> String

Return the 40-character lowercase hex SHA-1 of the most recent commit that modified `src/`.

Parameters
----------
repo_root : String, optional
    Path to the repository root directory.

Returns
-------
String
    40-character lowercase hex commit SHA-1 from `git log -n 1 --format=%H -- src`.

Raises
------
ErrorException
    If git does not return a 40-character lowercase hex SHA-1.
"""
function last_src_commit(repo_root::String=normpath(joinpath(@__DIR__, "..")))::String
    sha = readchomp(`git -C $repo_root log -n 1 --format=%H -- src`)
    !occursin(r"^[0-9a-f]{40}$", sha) &&
        error("Invalid src/ commit SHA returned by git: '$sha'")
    return sha
end

"""
    check_clean_source_tree(repo_root::String=normpath(joinpath(@__DIR__, "..")); record="benchmarks/records/porous_shell.toml") -> Tuple{Bool,String}

Verify that `src/`, `benchmarks/` and `tools/` have no uncommitted modifications before exporting a validation record.
The record file itself is excluded.

Parameters
----------
repo_root : String, optional
    Path to the repository root directory.
record : String, optional
    Path of the record file relative to `repo_root`, left out of the check (default
    `benchmarks/records/porous_shell.toml`).

Returns
-------
Tuple{Bool,String}
    `true` and empty string when `src/`, `benchmarks/` and `tools/` are clean; `false` and `git status --porcelain` output otherwise.

Raises
------
ProcessFailedException
    If `git status` fails, for example because `repo_root` is not a git repository.
"""
function check_clean_source_tree(
    repo_root::String=normpath(joinpath(@__DIR__, ".."));
    record::String="benchmarks/records/porous_shell.toml",
)::Tuple{Bool,String}
    dirty = strip(
        readchomp(
            `git -C $repo_root status --porcelain -- src benchmarks tools ":(exclude)$record"`,
        ),
    )
    return (isempty(dirty), String(dirty))
end

"""
    metric_status(value, reference, atol, rtol, comparison) -> String

Return `"PASS"` when `value` meets its tolerance against `reference`, `"FAIL"` otherwise.

Parameters
----------
value : Float64
    Measured value.
reference : Float64
    Reference value; `0.0` for `comparison = "le_atol"`.
atol : Float64
    Absolute tolerance (upper bound for `le_atol`).
rtol : Float64
    Relative tolerance (used by `isapprox`).
comparison : String
    `"isapprox"` or `"le_atol"`; any other string gives `"FAIL"`.

Returns
-------
String
    `"PASS"` or `"FAIL"`.
"""
function metric_status(
    value::Float64, reference::Float64, atol::Float64, rtol::Float64, comparison::String
)::String
    passed = if comparison == "isapprox"
        isapprox(value, reference; atol=atol, rtol=rtol)
    elseif comparison == "le_atol"
        iszero(reference) && isfinite(value) && 0.0 <= value <= atol
    else
        false
    end
    return passed ? "PASS" : "FAIL"
end

"""
    metric_specs(metrics, refs) -> Vector{NamedTuple}

Return the record entry of every metric, in the order of `EXPECTED_METRIC_IDS`, with the
limits of `PorousShellValidation.expected_metric_limits`: the shell depth within one cell
or 15 % of the 1-D depth, whichever is larger, the maximum temperature error at most 0.10
and the energy residual at most 1e-9; all are fixed before any run.

Parameters
----------
metrics : Dict{String,Float64}
    Measured metric values by id.
refs : Dict{String,Float64}
    Reference values by id, from `expected_metric_references`.

Returns
-------
Vector{NamedTuple}
    One entry per metric with `id`, `case`, `description`, `nx`, `value`, `reference`,
    `atol`, `rtol`, `comparison` and `units`.

Raises
------
KeyError
    If `metrics` or `refs` lacks one of the metric ids.
ErrorException
    If `PorousShellValidation` is not loaded.
"""
function metric_specs(metrics::Dict{String,Float64}, refs::Dict{String,Float64})
    PSV = porous_shell_validation_module()
    nx = PSV.NX_REFERENCE
    limits = PSV.expected_metric_limits(; d_shell_ref=refs["shell_depth_nx$(nx)"])
    texts = Dict(
        "shell_depth_nx$(nx)" => (
            "Median shell depth of the rock markers against the 1-D shell depth (Nx=$(nx))",
            "m",
        ),
        "tmax_rel_l2_nx$(nx)" => (
            "Relative L2 error of the binned peak temperature rise against the 1-D reference (Nx=$(nx))",
            "dimensionless",
        ),
        "energy_lse_residual_nx$(nx)" =>
            ("Scaled energy residual of the last thermal step (Nx=$(nx))", "dimensionless"),
    )
    return [
        (;
            id,
            case=PSV.CASE_ID,
            description=texts[id][1],
            nx,
            value=metrics[id],
            reference=refs[id],
            atol=limits[id].atol,
            rtol=limits[id].rtol,
            comparison=limits[id].comparison,
            units=texts[id][2],
        ) for id in PSV.EXPECTED_METRIC_IDS
    ]
end

"""
    porous_shell_validation_module(hosts=(parentmodule(@__MODULE__), Main)) -> Module

Return the `PorousShellValidation` module (`benchmarks/validation_porous_shell.jl`) of the
first module in `hosts` that defines it: by default the module that includes this file, or
`Main` when that module has none.

Parameters
----------
hosts : Tuple, optional
    Modules searched in order (default: the module that includes this file, then `Main`).

Returns
-------
Module
    The validation module.

Raises
------
ErrorException
    If no module of `hosts` defines `PorousShellValidation`.
"""
function porous_shell_validation_module(
    hosts::Tuple=(parentmodule(@__MODULE__), Main)
)::Module
    for host in hosts
        if Base.invokelatest(isdefined, host, :PorousShellValidation)
            return Base.invokelatest(getfield, host, :PorousShellValidation)
        end
    end
    return error(
        "PorousShellValidation is not loaded; include benchmarks/validation_porous_shell.jl first",
    )
end

"""
    write_porous_shell_record(record_path::String=DEFAULT_EXPORT_RECORD_PATH; src_commit=nothing, bench=nothing, accept_physics_change=false) -> NamedTuple

Run the porous shell benchmark (or use a pre-evaluated `bench`), validate the candidate
record with `validate_record_dict`, verify that no metric leaves its physics regression band around
the existing `record_path` (unless `accept_physics_change=true`), and write `record_path`.
Requires `PorousShellValidation` (`benchmarks/validation_porous_shell.jl`), see
`porous_shell_validation_module`.

Parameters
----------
record_path : String, optional
    Output path for the TOML validation record.
src_commit : Union{String,Nothing}, optional
    40-character lowercase hex SHA-1 of the last commit modifying `src/`. When `nothing`, queries `last_src_commit()`.
bench : Any, optional
    Result of `PorousShellValidation.evaluate_porous_shell_benchmark()`.
accept_physics_change : Bool, optional
    Allow overwriting an existing record when a metric moves outside `PHYSICS_REGRESSION_BANDS` around the existing value (default `false`).

Returns
-------
NamedTuple
    Evaluated benchmark summary.

Raises
------
DomainError
    If `isempty(record_path)` or `src_commit` is not a 40-character lowercase hex SHA-1.
ErrorException
    If a metric leaves its regression band around the existing record without `accept_physics_change=true`, if `validate_record_dict` fails, or if `PorousShellValidation` is not loaded.
KeyError
    If `bench.metrics` lacks one of the metric ids.
"""
function write_porous_shell_record(
    record_path::String=DEFAULT_EXPORT_RECORD_PATH;
    src_commit::Union{String,Nothing}=nothing,
    bench=nothing,
    accept_physics_change::Bool=false,
)
    PSV = porous_shell_validation_module()
    isempty(record_path) && throw(DomainError(record_path, "Record path must be non-empty"))
    commit_str = src_commit === nothing ? last_src_commit() : src_commit
    !occursin(r"^[0-9a-f]{40}$", commit_str) && throw(
        DomainError(
            commit_str, "src_commit must be a 40-character lowercase hexadecimal SHA-1"
        ),
    )

    bench_res = bench === nothing ? PSV.evaluate_porous_shell_benchmark() : bench
    metrics = bench_res.metrics

    if isfile(record_path) && !accept_physics_change
        ok_bands, band_errs = PSV.check_physics_regression_bands(metrics, record_path)
        if !ok_bands
            error(
                "Refusing to overwrite $record_path without --accept-physics-change: " *
                join(band_errs, "; "),
            )
        end
    end

    cfg = PSV.make_case_config()
    cfg_sha = PSV.canonical_config_sha256(cfg)
    refs = PSV.expected_metric_references(;
        d_shell_ref=PSV.bench_shell_reference(bench_res)
    )
    specs = metric_specs(metrics, refs)
    ts_utc = Dates.format(now(Dates.UTC), "yyyy-mm-ddTHH:MM:SSZ")

    io = IOBuffer()
    println(
        io,
        "# Validation record for the porous shell of a planetesimal with the dry porosity laws of Henke et al. (2012).",
    )
    println(io, "# Generated by tools/export_validation_porous_shell.jl.")
    println(
        io,
        "# metadata.src_commit records the latest commit that modified src/ (git log -n 1 --format=%H -- src).",
    )
    println(io)
    println(io, "[metadata]")
    println(io, "id = \"porous_shell\"")
    println(io, "title = \"Porous Shell Validation\"")
    println(io, "page = \"validation/porous_shell.md\"")
    println(io, "src_commit = \"$(commit_str)\"")
    println(io, "timestamp_utc = \"$(ts_utc)\"")
    println(io, "julia_version = \"$(VERSION)\"")
    println(io, "config_sha256 = \"$(cfg_sha)\"")
    println(io, "reference_type = \"1d_fv\"")
    println(
        io,
        "reference_citation = \"One-dimensional spherical finite-volume solution with the same time step: backward Euler conduction with the porous conductivity of Lichtenberg et al. (2016a, Eq. 12), 26Al heating, cold pressing and hot sintering of Henke et al. (2012)\"",
    )
    println(io)
    cfg_sum = PSV.expected_config_summary(cfg)
    println(io, "[config]")
    println(io, "Nx = $(cfg_sum["Nx"])")
    for key in ("Nr_reference", "n_steps", "seed")
        println(io, "$(key) = $(Int(cfg_sum[key]))")
    end
    for key in (
        "dt_yr",
        "t_form_yr",
        "phi_init",
        "R_planet_m",
        "L_box_m",
        "ratio_al",
        "k_rock_W_mK",
        "k_air_W_mK",
        "rhocp_rock_J_m3K",
        "rhocp_air_J_m3K",
        "T_space_K",
        "shell_dphi",
        "shell_depth_rel_bound",
    )
        println(io, "$(key) = $(repr(Float64(cfg_sum[key])))")
    end
    println(io, "porosity_model = \"$(cfg_sum["porosity_model"])\"")
    println(io, "spherical_metric = $(cfg_sum["spherical_metric"])")

    for m in specs
        status_str = metric_status(
            Float64(m.value),
            Float64(m.reference),
            Float64(m.atol),
            Float64(m.rtol),
            String(m.comparison),
        )
        println(io)
        println(io, "[[metrics]]")
        println(io, "id = \"$(m.id)\"")
        println(io, "case = \"$(m.case)\"")
        println(io, "description = \"$(m.description)\"")
        println(io, "nx = $(m.nx)")
        @printf(io, "value = %.16e\n", m.value)
        @printf(io, "reference = %.16e\n", m.reference)
        @printf(io, "atol = %.16e\n", m.atol)
        @printf(io, "rtol = %.16e\n", m.rtol)
        println(io, "comparison = \"$(m.comparison)\"")
        println(io, "units = \"$(m.units)\"")
        println(io, "status = \"$(status_str)\"")
    end

    toml_str = String(take!(io))
    candidate_dict = TOML.parse(toml_str)
    ok_valid, valid_errs = PSV.validate_record_dict(
        candidate_dict; require_doc_pages=false, check_doc_sync=false
    )
    if !ok_valid
        error(
            "Refusing to write invalid validation record to $record_path: " *
            join(valid_errs, "; "),
        )
    end

    mkpath(dirname(record_path))
    write(record_path, toml_str)
    if abspath(record_path) == abspath(DEFAULT_EXPORT_RECORD_PATH)
        PSV.update_record_doc_block!(candidate_dict) ||
            @warn "Documentation page has no record table to update"
    end
    return bench_res
end

"""
    run_export_cli(args::Vector{String}=ARGS; src_commit=nothing, bench=nothing, check_clean_tree=true, dirty_status=nothing, repo_root=normpath(joinpath(@__DIR__, ".."))) -> Int

Run the command-line export tool for the porous shell validation record. Requires
`PorousShellValidation` (`benchmarks/validation_porous_shell.jl`), see
`porous_shell_validation_module`. An argument other than `--accept-physics-change` that
starts with `-`, and a second record path, are usage errors. An error of the export other
than a `DomainError`, including a failed `git status` of the clean-tree check, is printed
with its stack trace and gives the exit code `1`.

Parameters
----------
args : Vector{String}, optional
    Command-line arguments (`[--accept-physics-change]` and at most one `[record_path]`).
src_commit : Union{String,Nothing}, optional
    Explicit `src/` commit SHA-1, or `nothing` to query `last_src_commit(repo_root)`.
bench : Any, optional
    Pre-evaluated benchmark tuple, or `nothing` to evaluate on demand.
check_clean_tree : Bool, optional
    Refuse to export when `src/`, `benchmarks/` or `tools/` has uncommitted changes (default `true`).
dirty_status : Union{String,Nothing}, optional
    Replacement for the `git status --porcelain` output of the clean-tree check.
repo_root : String, optional
    Repository root directory checked for `src/`, `benchmarks/` and `tools/` cleanliness.

Returns
-------
Int
    Exit code (`0` on success, `1` on export refusal, `2` on invalid command line).

Raises
------
DomainError
    If the record path is empty or `src_commit` is not a 40-character lowercase hex SHA-1.
InterruptException
    If the run is interrupted.
"""
function run_export_cli(
    args::Vector{String}=ARGS;
    src_commit::Union{String,Nothing}=nothing,
    bench=nothing,
    check_clean_tree::Bool=true,
    dirty_status::Union{String,Nothing}=nothing,
    repo_root::String=normpath(joinpath(@__DIR__, "..")),
)::Int
    accept_physics_change = false
    record_path = nothing
    for arg in args
        if arg == "--accept-physics-change"
            accept_physics_change = true
            continue
        end
        if startswith(arg, "-") || record_path !== nothing
            println(
                stderr,
                "Usage: julia --project=. tools/export_validation_porous_shell.jl [--accept-physics-change] [record_path]",
            )
            return 2
        end
        record_path = arg
    end
    record_path = something(record_path, DEFAULT_EXPORT_RECORD_PATH)

    try
        if check_clean_tree
            ok_clean, dirty_out = if dirty_status === nothing
                check_clean_source_tree(repo_root)
            else
                (isempty(strip(dirty_status)), String(strip(dirty_status)))
            end
            if !ok_clean
                println(
                    stderr,
                    "Refusing to export validation record with uncommitted changes in src/, benchmarks/ or tools/:\n$dirty_out",
                )
                return 1
            end
        end
        commit_val = src_commit === nothing ? last_src_commit(repo_root) : src_commit
        write_porous_shell_record(
            record_path;
            src_commit=commit_val,
            bench=bench,
            accept_physics_change=accept_physics_change,
        )
        println("Wrote $record_path")
        return 0
    catch err
        (err isa DomainError || err isa InterruptException) && rethrow()
        println(stderr, sprint(showerror, err, catch_backtrace()))
        return 1
    end
end

"""
    export_main(args::Vector{String}=ARGS; kwargs...) -> Int

Run `run_export_cli` as the script does: a `DomainError` (an empty record path or an invalid
`src_commit`) is printed and gives the exit code `1`.

Parameters
----------
args : Vector{String}, optional
    Command-line arguments, as for `run_export_cli`.
kwargs
    Keyword arguments passed on to `run_export_cli`.

Returns
-------
Int
    Exit code of `run_export_cli`, or `1` after a `DomainError`.

Raises
------
InterruptException
    If the run is interrupted.
"""
function export_main(args::Vector{String}=ARGS; kwargs...)::Int
    return try
        run_export_cli(args; kwargs...)
    catch err
        err isa DomainError || rethrow()
        println(stderr, sprint(showerror, err))
        1
    end
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "..", "benchmarks", "validation_porous_shell.jl"))
    exit(Base.invokelatest(PorousShellExport.export_main, ARGS))
end

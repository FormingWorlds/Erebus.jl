#!/usr/bin/env julia
# Export the melting shell validation records: the CI case from a live run, and the suite
# record of the case at Nx 129 from the result files of scripts/run_validation_melting_shell.jl.
# Usage: julia --project=. tools/export_validation_melting_shell.jl [--accept-physics-change] [record_path]
#        julia --project=. tools/export_validation_melting_shell.jl --suite <run_dir> [record_path]

isdefined(@__MODULE__, :PorousShellExport) ||
    include(joinpath(@__DIR__, "export_validation_porous_shell.jl"))

module MeltingShellExport

using Dates
using Printf
using TOML
using ..PorousShellExport: check_clean_source_tree, last_src_commit, metric_status

const DEFAULT_RECORD_PATH = normpath(
    joinpath(@__DIR__, "..", "benchmarks", "records", "melting_shell.toml")
)
const DEFAULT_SUITE_RECORD_PATH = normpath(
    joinpath(@__DIR__, "..", "benchmarks", "records", "melting_shell_suite.toml")
)
const RECORD_NOTE = "The melt radius row is a regression band: the runs with L = 0 and with 1.5 L stay inside it, so it does not separate them; the centre temperature row does."

"""
    melting_shell_validation_module(hosts=(parentmodule(@__MODULE__), Main)) -> Module

Return the `MeltingShellValidation` module of the first module in `hosts` that defines it.

# Throws

    - `ErrorException`: if no module of `hosts` defines it
"""
function melting_shell_validation_module(
    hosts::Tuple=(parentmodule(@__MODULE__), Main)
)::Module
    for host in hosts
        if Base.invokelatest(isdefined, host, :MeltingShellValidation)
            return Base.invokelatest(getfield, host, :MeltingShellValidation)
        end
    end
    return error(
        "MeltingShellValidation is not loaded; include benchmarks/validation_melting_shell.jl first",
    )
end

"""
    metric_description(id) -> String

The text of the record table for a metric id of the melting shell records.

# Throws

    - `ArgumentError`: if `id` has no known prefix
"""
function metric_description(id::AbstractString)::String
    MSV = melting_shell_validation_module()
    startswith(id, "enthalpy_residual_") &&
        return "Relative enthalpy residual of the 1-D reference over the run, the lag of the apparent heat capacity"
    nx = MSV.metric_nx(id)
    startswith(id, "centre_tmax_") &&
        return "Median peak temperature of the rock markers in the centre bin against the 1-D peak temperature (Nx=$(nx))"
    startswith(id, "melt_radius_") &&
        return "Melt radius against the 1-D melt radius, limit 1.5 cells, a regression band (Nx=$(nx))"
    startswith(id, "tmax_rel_l2_hot_") &&
        return "Relative L2 error of the binned peak temperature rise where the 1-D peak temperature is at least 1100 K (Nx=$(nx))"
    startswith(id, "energy_lse_residual_") &&
        return "Scaled energy residual of the last thermal step (Nx=$(nx))"
    return throw(ArgumentError("unknown metric id '$id'"))
end

"""
    record_text(; id, title, citation, comment, src_commit, config_sha256, summary, keys,
                metrics) -> String

The TOML text of a melting shell record: header comments, `[metadata]`, the `[config]`
entries of `summary` in the order of `keys`, and one `[[metrics]]` table per entry of
`metrics` (`(; id, case, nx, value, reference)`), with the description of
`metric_description`, the limit of `metric_limit`, the units of `metric_units` and the
status of `metric_status`.
"""
function record_text(;
    id, title, citation, comment, src_commit, config_sha256, summary, keys, metrics
)::String
    MSV = melting_shell_validation_module()
    io = IOBuffer()
    println(io, "# ", comment)
    println(io, "# ", RECORD_NOTE)
    println(
        io,
        "# metadata.src_commit records the latest commit that modified src/ (git log -n 1 --format=%H -- src).",
    )
    println(io)
    println(io, "[metadata]")
    println(io, "id = \"$id\"")
    println(io, "title = \"$title\"")
    println(io, "page = \"validation/melting_shell.md\"")
    println(io, "src_commit = \"$src_commit\"")
    println(
        io, "timestamp_utc = \"$(Dates.format(now(Dates.UTC), "yyyy-mm-ddTHH:MM:SSZ"))\""
    )
    println(io, "julia_version = \"$(VERSION)\"")
    println(io, "config_sha256 = \"$config_sha256\"")
    println(io, "reference_type = \"1d_fv\"")
    println(io, "reference_citation = \"$citation\"")
    println(io)
    println(io, "[config]")
    for key in keys
        v = summary[key]
        text = if v isa AbstractString || v isa Bool || v isa Integer
            sprint(TOML.print, Dict(key => v))
        else
            "$key = $(repr(Float64(v)))\n"
        end
        print(io, text)
    end
    for m in metrics
        lim = MSV.metric_limit(m.id)
        status = metric_status(
            Float64(m.value), Float64(m.reference), lim.atol, lim.rtol, lim.comparison
        )
        println(io)
        println(io, "[[metrics]]")
        println(io, "id = \"$(m.id)\"")
        println(io, "case = \"$(m.case)\"")
        println(io, "description = \"$(metric_description(m.id))\"")
        println(io, "nx = $(m.nx)")
        @printf(io, "value = %.16e\n", m.value)
        @printf(io, "reference = %.16e\n", m.reference)
        @printf(io, "atol = %.16e\n", lim.atol)
        @printf(io, "rtol = %.16e\n", lim.rtol)
        println(io, "comparison = \"$(lim.comparison)\"")
        println(io, "units = \"$(MSV.metric_units(m.id))\"")
        println(io, "status = \"$status\"")
    end
    return String(take!(io))
end

# Validate a candidate record text and write it, with its table on the page for a default path.
function _write_validated(text::String, record_path::String, default_path::String)
    MSV = melting_shell_validation_module()
    data = TOML.parse(text)
    ok, errs = MSV.PSV.validate_record_dict(
        data; require_doc_pages=false, check_doc_sync=false
    )
    ok || error(
        "Refusing to write invalid validation record to $record_path: " * join(errs, "; "),
    )
    mkpath(dirname(record_path))
    write(record_path, text)
    if abspath(record_path) == abspath(default_path)
        MSV.PSV.update_record_doc_block!(data) ||
            @warn "Documentation page has no record table to update"
    end
    return data
end

# A 40-character lowercase hex SHA-1, or a DomainError.
function _checked_commit(commit::String)::String
    occursin(r"^[0-9a-f]{40}$", commit) || throw(
        DomainError(
            commit, "src_commit must be a 40-character lowercase hexadecimal SHA-1"
        ),
    )
    return commit
end

"""
    write_melting_shell_record(record_path=DEFAULT_RECORD_PATH; src_commit=nothing,
                               bench=nothing, accept_physics_change=false) -> NamedTuple

Run the CI case (or use `bench`, an `evaluate_melting_shell_benchmark` result), check the
regression bands against an existing record unless `accept_physics_change`, validate the
candidate record and write it; for the default path also its table on the page.

# Returns

The benchmark result.

# Throws

    - `DomainError`: if `record_path` is empty or `src_commit` is not a SHA-1
    - `ErrorException`: if a metric leaves its band, or the record fails
      `validate_record_dict`
"""
function write_melting_shell_record(
    record_path::String=DEFAULT_RECORD_PATH;
    src_commit::Union{String,Nothing}=nothing,
    bench=nothing,
    accept_physics_change::Bool=false,
)
    MSV = melting_shell_validation_module()
    isempty(record_path) && throw(DomainError(record_path, "Record path must be non-empty"))
    commit = _checked_commit(something(src_commit, last_src_commit()))
    b = bench === nothing ? MSV.evaluate_melting_shell_benchmark() : bench
    if isfile(record_path) && !accept_physics_change
        ok, errs = MSV.check_physics_regression_bands(b.metrics, record_path)
        ok || error(
            "Refusing to overwrite $record_path without --accept-physics-change: " *
            join(errs, "; "),
        )
    end
    cfg = MSV.make_case_config()
    refs = MSV.expected_metric_references(; ref=b.reference)
    metrics = [
        (;
            id,
            case=MSV.CASE_ID,
            nx=MSV.metric_nx(id),
            value=b.metrics[id],
            reference=refs[id],
        ) for id in MSV.EXPECTED_METRIC_IDS
    ]
    text = record_text(;
        id=MSV.RECORD_ID,
        title="Melting Shell Validation",
        citation="One-dimensional spherical finite-volume solution with the same time step: the porous shell reference with silicate latent heat as an apparent heat capacity at the temperature of the step start (solve_porous_sphere in benchmarks/validation_porous_shell.jl)",
        comment="Validation record of the melting shell: silicate melting with latent heat in the porous shell case. Generated by tools/export_validation_melting_shell.jl.",
        src_commit=commit,
        config_sha256=MSV.canonical_config_sha256(cfg),
        summary=MSV.expected_config_summary(cfg),
        keys=MSV.PSV.REQUIRED_CONFIG_KEYS[MSV.RECORD_ID],
        metrics,
    )
    _write_validated(text, record_path, DEFAULT_RECORD_PATH)
    return b
end

"""
    read_suite_result(run_dir) -> NamedTuple

Read the result row of the run at `NX_SUITE` with the nominal latent heat from `run_dir`:
`header.csv` and `result_nx<NX_SUITE>.csv` of `scripts/run_validation_melting_shell.jl`.
The row must hold the grid, the latent heat and the cell size of the case, the pinned 1-D
references of `SUITE_REFERENCE_VALUES` and an energy residual within its limit.

# Returns

`(; T_centre, r_melt, tmax_rel_l2_hot)` [K, m, -].

# Throws

    - `ArgumentError`: if a file is missing, the header differs from
      `SUITE_RESULT_HEADER`, or the row does not match the case
"""
function read_suite_result(run_dir::String)
    MSV = melting_shell_validation_module()
    rows_of(path) = [split(strip(l), ",") for l in eachline(path) if !isempty(strip(l))]
    header_path = joinpath(run_dir, "header.csv")
    isfile(header_path) || throw(ArgumentError("missing $header_path"))
    header = String.(only(rows_of(header_path)))
    header == collect(MSV.SUITE_RESULT_HEADER) ||
        throw(ArgumentError("$header_path has columns $header"))
    nx = MSV.NX_SUITE
    path = joinpath(run_dir, "result_nx$(nx).csv")
    isfile(path) || throw(ArgumentError("missing $path"))
    rows = rows_of(path)
    (length(rows) == 1 && length(rows[1]) == length(header)) ||
        throw(ArgumentError("$path must hold 1 row of $(length(header)) fields"))
    row = rows[1]
    num(name) = parse(Float64, row[findfirst(==(name), header)])
    e_res = num("energy_res")
    matches =
        row[1] == "nx$(nx)" &&
        parse(Int, row[2]) == nx &&
        isapprox(num("L_melt"), MSV.L_MELT; rtol=1.0e-12) &&
        isapprox(num("dx"), MSV.PSV.L_BOX / (nx - 1); rtol=1.0e-12) &&
        isapprox(
            num("T_centre_ref"),
            MSV.suite_reference_value("centre_tmax_nx$nx");
            rtol=1.0e-12,
        ) &&
        isapprox(
            num("r_melt_ref"), MSV.suite_reference_value("melt_radius_nx$nx"); rtol=1.0e-12
        ) &&
        isfinite(e_res) &&
        0.0 <= e_res <= MSV.ENERGY_RESIDUAL_LIMIT
    matches || throw(ArgumentError("$path does not match the suite case at Nx = $nx"))
    return (;
        T_centre=num("T_centre"),
        r_melt=num("r_melt"),
        tmax_rel_l2_hot=num("tmax_rel_l2_hot"),
    )
end

"""
    write_suite_record(run_dir; record_path=DEFAULT_SUITE_RECORD_PATH, src_commit=nothing)
        -> Dict{String,Any}

Read the result of the run at `NX_SUITE` from `run_dir`, write the suite record and, for
the default path, its table on the page. A candidate record that fails
`validate_record_dict` is not written.

# Throws

    - `DomainError`: if `record_path` is empty or `src_commit` is not a SHA-1
    - `ArgumentError`: from `read_suite_result`
    - `ErrorException`: if the record fails `validate_record_dict`
"""
function write_suite_record(
    run_dir::String;
    record_path::String=DEFAULT_SUITE_RECORD_PATH,
    src_commit::Union{String,Nothing}=nothing,
)
    MSV = melting_shell_validation_module()
    isempty(record_path) && throw(DomainError(record_path, "Record path must be non-empty"))
    commit = _checked_commit(something(src_commit, last_src_commit()))
    res = read_suite_result(run_dir)
    values = (res.T_centre, res.r_melt, res.tmax_rel_l2_hot)
    metrics = [
        (;
            id,
            case=MSV.CASE_ID,
            nx=MSV.NX_SUITE,
            value=v,
            reference=MSV.suite_reference_value(id),
        ) for (id, v) in zip(MSV.suite_metric_ids(), values)
    ]
    text = record_text(;
        id=MSV.SUITE_RECORD_ID,
        title="Melting Shell Suite Validation",
        citation="One-dimensional spherical finite-volume solution of the run with the same time step and number of steps (solve_porous_sphere in benchmarks/validation_porous_shell.jl, with latent heat)",
        comment="Validation record of the melting shell at Nx = 129. Generated by tools/export_validation_melting_shell.jl from the result files of scripts/run_validation_melting_shell.jl.",
        src_commit=commit,
        config_sha256=MSV.suite_config_sha256(),
        summary=MSV.expected_suite_config_summary(),
        keys=MSV.PSV.REQUIRED_CONFIG_KEYS[MSV.SUITE_RECORD_ID],
        metrics,
    )
    return _write_validated(text, record_path, DEFAULT_SUITE_RECORD_PATH)
end

"""
    run_export_cli(args=ARGS; src_commit=nothing, bench=nothing, check_clean_tree=true,
                   repo_root=normpath(joinpath(@__DIR__, ".."))) -> Int

Run the export from the command line: without `--suite` the CI record, with
`--suite <run_dir>` the suite record; then at most one record path. An unknown option or
a second path is a usage error. An error of the export other than a `DomainError` is
printed with its stack trace.

# Returns

Exit code: 0 on success, 1 on a refused or failed export, 2 on a usage error.

# Throws

    - `DomainError`: if `src_commit` is not a SHA-1
    - `InterruptException`: if the run is interrupted
"""
function run_export_cli(
    args::Vector{String}=ARGS;
    src_commit::Union{String,Nothing}=nothing,
    bench=nothing,
    check_clean_tree::Bool=true,
    repo_root::String=normpath(joinpath(@__DIR__, "..")),
)::Int
    usage = "Usage: julia --project=. tools/export_validation_melting_shell.jl [--accept-physics-change] [record_path] | --suite <run_dir> [record_path]"
    accept, suite_dir, record_path = false, nothing, nothing
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--accept-physics-change"
            accept = true
        elseif a == "--suite" && i < length(args) && suite_dir === nothing
            suite_dir = args[i + 1]
            i += 1
        elseif startswith(a, "-") || record_path !== nothing
            println(stderr, usage)
            return 2
        else
            record_path = a
        end
        i += 1
    end
    default = suite_dir === nothing ? DEFAULT_RECORD_PATH : DEFAULT_SUITE_RECORD_PATH
    record_path = something(record_path, default)
    try
        if check_clean_tree
            rel = relpath(record_path, repo_root)
            ok, dirty = check_clean_source_tree(repo_root; record=rel)
            if !ok
                println(
                    stderr,
                    "Refusing to export validation record with uncommitted changes in src/, benchmarks/ or tools/:\n$dirty",
                )
                return 1
            end
        end
        commit = something(src_commit, last_src_commit(repo_root))
        if suite_dir === nothing
            write_melting_shell_record(
                record_path; src_commit=commit, bench, accept_physics_change=accept
            )
        else
            write_suite_record(suite_dir; record_path, src_commit=commit)
        end
        println("Wrote $record_path")
        return 0
    catch err
        (err isa DomainError || err isa InterruptException) && rethrow()
        println(stderr, sprint(showerror, err, catch_backtrace()))
        return 1
    end
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "..", "benchmarks", "validation_melting_shell.jl"))
    exit(Base.invokelatest(MeltingShellExport.run_export_cli, ARGS))
end

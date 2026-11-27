#!/usr/bin/env julia
# Export the magma-ocean conductivity CI record and Nx = 129 suite record.
# Usage: julia --project=. tools/export_validation_magma_ocean_conductivity.jl [--accept-physics-change | --suite <run_dir>] [record_path]

isdefined(@__MODULE__, :PorousShellExport) ||
    include(joinpath(@__DIR__, "export_validation_porous_shell.jl"))

module MagmaOceanConductivityExport

using Dates
using Printf
using TOML
using ..PorousShellExport: check_clean_source_tree, last_src_commit, metric_status

const DEFAULT_RECORD_PATH = normpath(
    joinpath(@__DIR__, "..", "benchmarks", "records", "magma_ocean_conductivity.toml")
)
const DEFAULT_SUITE_RECORD_PATH = normpath(
    joinpath(@__DIR__, "..", "benchmarks", "records", "magma_ocean_conductivity_suite.toml")
)
const RECORD_NOTE = "The mush and melt radius rows are regression bands (1.5 cells); the centre temperature and plateau drop rows separate convectively buffered runs from unbuffered or mis-scaled exponents."

"""
    magma_ocean_validation_module(hosts=(parentmodule(@__MODULE__), Main)) -> Module

Return the `MagmaOceanConductivityValidation` module of the first module in `hosts`.
"""
function magma_ocean_validation_module(
    hosts::Tuple=(parentmodule(@__MODULE__), Main)
)::Module
    for host in hosts
        if Base.invokelatest(isdefined, host, :MagmaOceanConductivityValidation)
            return Base.invokelatest(getfield, host, :MagmaOceanConductivityValidation)
        end
    end
    return error(
        "MagmaOceanConductivityValidation is not loaded; include benchmarks/validation_magma_ocean_conductivity.jl first",
    )
end

"""
    metric_description(id) -> String

Human-readable description of a metric id for the magma-ocean conductivity records.
"""
function metric_description(id::AbstractString)::String
    MOSV = magma_ocean_validation_module()
    startswith(id, "enthalpy_residual_") &&
        return "Relative enthalpy residual of the 1-D reference over the run"
    nx = MOSV.metric_nx(id)
    startswith(id, "centre_tmax_") &&
        return "Median peak temperature of rock markers in the centre bin against the 1-D peak temperature (Nx=$(nx))"
    startswith(id, "plateau_drop_") &&
        return "Inner plateau temperature drop from the centre bin to 0.5 R_p (Nx=$(nx))"
    startswith(id, "mush_radius_") &&
        return "Mush radius (F_max > 0.30) against the 1-D mush radius, limit 1.5 cells, a regression band (Nx=$(nx))"
    startswith(id, "melt_radius_") &&
        return "Melt radius (F_max > 0) against the 1-D melt radius, limit 1.5 cells, a regression band (Nx=$(nx))"
    startswith(id, "tmax_rel_l2_hot_") &&
        return "Relative L2 error of the binned peak temperature rise where the 1-D peak temperature is at least 1100 K (Nx=$(nx))"
    startswith(id, "energy_lse_residual_") &&
        return "Scaled energy residual of the last thermal step (Nx=$(nx))"
    return throw(ArgumentError("unknown metric id '$id'"))
end

"""
    record_text(; id, title, citation, comment, src_commit, config_sha256, summary, keys,
                metrics) -> String

Format a magma-ocean conductivity validation record as TOML text.
"""
function record_text(;
    id, title, citation, comment, src_commit, config_sha256, summary, keys, metrics
)::String
    MOSV = magma_ocean_validation_module()
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
    println(io, "page = \"validation/magma_ocean_conductivity.md\"")
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
        lim = MOSV.metric_limit(m.id)
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
        println(io, "units = \"$(MOSV.metric_units(m.id))\"")
        println(io, "status = \"$status\"")
    end
    return String(take!(io))
end

"""
    _write_validated(text, record_path, default_path) -> Dict{String,Any}

Validate TOML `text`, write it to `record_path`, and update the documentation block when writing `default_path`.
"""
function _write_validated(text::String, record_path::String, default_path::String)
    MOSV = magma_ocean_validation_module()
    data = TOML.parse(text)
    ok, errs = MOSV.PSV.validate_record_dict(
        data; require_doc_pages=false, check_doc_sync=false
    )
    ok || error(
        "Refusing to write invalid validation record to $record_path: " * join(errs, "; "),
    )
    mkpath(dirname(record_path))
    write(record_path, text)
    if abspath(record_path) == abspath(default_path)
        MOSV.PSV.update_record_doc_block!(data) ||
            @warn "Documentation page has no record table to update"
    end
    return data
end

"""
    _checked_commit(commit) -> String

Verify that `commit` is a 40-character lowercase hexadecimal SHA-1 digest.
"""
function _checked_commit(commit::String)::String
    occursin(r"^[0-9a-f]{40}$", commit) || throw(
        DomainError(
            commit, "src_commit must be a 40-character lowercase hexadecimal SHA-1"
        ),
    )
    return commit
end

"""
    write_magma_ocean_record(record_path=DEFAULT_RECORD_PATH; src_commit=nothing,
                             bench=nothing, accept_physics_change=false) -> NamedTuple

Run (or accept `bench`) and write the `Nx = 65` CI validation record.
"""
function write_magma_ocean_record(
    record_path::String=DEFAULT_RECORD_PATH;
    src_commit::Union{String,Nothing}=nothing,
    bench=nothing,
    accept_physics_change::Bool=false,
)
    MOSV = magma_ocean_validation_module()
    isempty(record_path) && throw(DomainError(record_path, "Record path must be non-empty"))
    commit = _checked_commit(something(src_commit, last_src_commit()))
    b = bench === nothing ? MOSV.evaluate_magma_ocean_benchmark() : bench
    if isfile(record_path) && !accept_physics_change
        ok, errs = MOSV.check_physics_regression_bands(b.metrics, record_path)
        ok || error(
            "Refusing to overwrite $record_path without --accept-physics-change: " *
            join(errs, "; "),
        )
    end
    cfg = MOSV.make_case_config()
    refs = MOSV.expected_metric_references(; ref=b.reference)
    metrics = [
        (;
            id,
            case=MOSV.CASE_ID,
            nx=MOSV.metric_nx(id),
            value=b.metrics[id],
            reference=refs[id],
        ) for id in MOSV.EXPECTED_METRIC_IDS
    ]
    text = record_text(;
        id=MOSV.RECORD_ID,
        title="Magma Ocean Effective Conductivity Validation",
        citation="One-dimensional spherical finite-volume solution with the same time step: porous shell reference with silicate latent heat and regularized soft-turbulence conductivity (solve_porous_sphere in benchmarks/validation_porous_shell.jl)",
        comment="Validation record of convecting melt effective conductivity in a porous planetesimal. Generated by tools/export_validation_magma_ocean_conductivity.jl.",
        src_commit=commit,
        config_sha256=MOSV.canonical_config_sha256(cfg),
        summary=MOSV.expected_config_summary(cfg),
        keys=MOSV.PSV.REQUIRED_CONFIG_KEYS[MOSV.RECORD_ID],
        metrics,
    )
    _write_validated(text, record_path, DEFAULT_RECORD_PATH)
    return b
end

"""
    read_suite_result(run_dir) -> NamedTuple

Read the `Nx = 129` result row from `run_dir` (`header.csv` and `result_nx129.csv`).
"""
function read_suite_result(run_dir::String)
    MOSV = magma_ocean_validation_module()
    rows_of(path) = [split(strip(l), ",") for l in eachline(path) if !isempty(strip(l))]
    header_path = joinpath(run_dir, "header.csv")
    isfile(header_path) || throw(ArgumentError("missing $header_path"))
    header = String.(only(rows_of(header_path)))
    header == collect(MOSV.SUITE_RESULT_HEADER) ||
        throw(ArgumentError("$header_path has columns $header"))
    nx = MOSV.NX_SUITE
    path = joinpath(run_dir, "result_nx$(nx).csv")
    isfile(path) || throw(ArgumentError("missing $path"))
    rows = rows_of(path)
    (length(rows) == 1 && length(rows[1]) == length(header)) ||
        throw(ArgumentError("$path must hold 1 row of $(length(header)) fields"))
    row = rows[1]
    col(name) = row[findfirst(==(name), header)]
    num(name) = parse(Float64, col(name))
    e_res = num("energy_res")
    matches =
        row[1] == "nx$(nx)" &&
        parse(Int, row[2]) == nx &&
        col("soft_turbulence") == "true" &&
        isapprox(num("turb_exponent"), MOSV.TURB_EXP_REF; rtol=1.0e-12) &&
        isapprox(num("dx"), MOSV.PSV.L_BOX / (nx - 1); rtol=1.0e-12) &&
        isapprox(
            num("T_centre_ref"),
            MOSV.suite_reference_value("centre_tmax_nx$nx");
            rtol=1.0e-12,
        ) &&
        isapprox(
            num("r_mush_ref"), MOSV.suite_reference_value("mush_radius_nx$nx"); rtol=1.0e-12
        ) &&
        isapprox(
            num("r_melt_ref"), MOSV.suite_reference_value("melt_radius_nx$nx"); rtol=1.0e-12
        ) &&
        isfinite(e_res) &&
        0.0 <= e_res <= MOSV.ENERGY_RESIDUAL_LIMIT
    matches || throw(ArgumentError("$path does not match the suite case at Nx = $nx"))
    return (;
        T_centre=num("T_centre"),
        plateau_drop=num("plateau_drop"),
        r_mush=num("r_mush"),
        r_melt=num("r_melt"),
        tmax_rel_l2_hot=num("tmax_rel_l2_hot"),
    )
end

"""
    write_suite_record(run_dir; record_path=DEFAULT_SUITE_RECORD_PATH, src_commit=nothing)
        -> Dict{String,Any}

Write the `Nx = 129` suite validation record from `run_dir`.
"""
function write_suite_record(
    run_dir::String;
    record_path::String=DEFAULT_SUITE_RECORD_PATH,
    src_commit::Union{String,Nothing}=nothing,
)
    MOSV = magma_ocean_validation_module()
    isempty(record_path) && throw(DomainError(record_path, "Record path must be non-empty"))
    commit = _checked_commit(something(src_commit, last_src_commit()))
    res = read_suite_result(run_dir)
    values = (res.T_centre, res.plateau_drop, res.r_mush, res.r_melt, res.tmax_rel_l2_hot)
    metrics = [
        (;
            id,
            case=MOSV.CASE_ID,
            nx=MOSV.NX_SUITE,
            value=v,
            reference=MOSV.suite_reference_value(id),
        ) for (id, v) in zip(MOSV.suite_metric_ids(), values)
    ]
    text = record_text(;
        id=MOSV.SUITE_RECORD_ID,
        title="Magma Ocean Effective Conductivity Suite Validation",
        citation="One-dimensional spherical finite-volume solution with the same time step and number of steps (solve_porous_sphere in benchmarks/validation_porous_shell.jl, with latent heat and soft turbulence)",
        comment="Validation record of convecting melt effective conductivity at Nx = 129. Generated by tools/export_validation_magma_ocean_conductivity.jl from scripts/run_validation_magma_ocean_conductivity.jl.",
        src_commit=commit,
        config_sha256=MOSV.suite_config_sha256(),
        summary=MOSV.expected_suite_config_summary(),
        keys=MOSV.PSV.REQUIRED_CONFIG_KEYS[MOSV.SUITE_RECORD_ID],
        metrics,
    )
    return _write_validated(text, record_path, DEFAULT_SUITE_RECORD_PATH)
end

"""
    run_export_cli(args=ARGS; src_commit=nothing, bench=nothing, check_clean_tree=true,
                   repo_root=normpath(joinpath(@__DIR__, ".."))) -> Int

Command-line entry point for exporting CI and suite magma-ocean conductivity records.
"""
function run_export_cli(
    args::Vector{String}=ARGS;
    src_commit::Union{String,Nothing}=nothing,
    bench=nothing,
    check_clean_tree::Bool=true,
    repo_root::String=normpath(joinpath(@__DIR__, "..")),
)::Int
    usage = "Usage: julia --project=. tools/export_validation_magma_ocean_conductivity.jl [--accept-physics-change] [record_path] | --suite <run_dir> [record_path]"
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
            write_magma_ocean_record(
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
    include(
        joinpath(@__DIR__, "..", "benchmarks", "validation_magma_ocean_conductivity.jl")
    )
    exit(Base.invokelatest(MagmaOceanConductivityExport.run_export_cli, ARGS))
end

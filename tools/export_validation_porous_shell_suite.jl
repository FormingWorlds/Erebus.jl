#!/usr/bin/env julia
# Export the porous shell suite record TOML from the result files of the full runs.
# Usage: julia --project=. tools/export_validation_porous_shell_suite.jl <suite_dir> [record_path]

isdefined(@__MODULE__, :PorousShellExport) ||
    include(joinpath(@__DIR__, "export_validation_porous_shell.jl"))

module PorousShellSuiteExport

using Dates
using Printf
using TOML
using ..PorousShellExport:
    check_clean_source_tree, last_src_commit, metric_status, porous_shell_validation_module

const DEFAULT_SUITE_RECORD_PATH = normpath(
    joinpath(@__DIR__, "..", "benchmarks", "records", "porous_shell_suite.toml")
)
const SUITE_RECORD_RELPATH = "benchmarks/records/porous_shell_suite.toml"

"""
    read_suite_results(suite_dir) -> Vector{NamedTuple}

Read the result rows of the full runs from `suite_dir`, as written by
`scripts/run_validation_porous_shell_suite.jl`: `header.csv` and one `result_<tag>.csv`
per case of `SUITE_CASES`. Each row must describe its case: the tag, `Nx`, the formation
time, the initial porosity, the number of steps, the air conductivity (`-1` for the
matched air) and the cell size `L_BOX / (Nx - 1)`. The energy residual `energy_res` of
each row must be finite and in `[0, ENERGY_RESIDUAL_LIMIT]`, its 1-D depth `d_ref` must
equal the pinned depth of its case (`suite_reference_depth`), and the folder must hold no
`result_<tag>.csv` of a tag outside `SUITE_CASES`.

Parameters
----------
suite_dir : String
    Folder with the result files.

Returns
-------
Vector{NamedTuple}
    One `(; tag, nx, t_form_myr, phi_init, n_steps, k_air, d_shell, d_ref, tmax_rel_l2)`
    per case, in the order of `SUITE_CASES`; depths in m, `t_form_myr` in Myr.

Raises
------
ArgumentError
    If the header differs from `SUITE_RESULT_HEADER`, a result file is missing or does not
    hold one row of the header length, a row does not match its case or has an energy
    residual above `ENERGY_RESIDUAL_LIMIT`, or a result file has a tag outside
    `SUITE_CASES`.
ErrorException
    If `PorousShellValidation` is not loaded.
"""
function read_suite_results(suite_dir::String)
    PSV = porous_shell_validation_module()
    rows_of(path) = [split(strip(l), ",") for l in eachline(path) if !isempty(strip(l))]
    header_path = joinpath(suite_dir, "header.csv")
    isfile(header_path) || throw(ArgumentError("missing $header_path"))
    header = String.(only(rows_of(header_path)))
    header == collect(PSV.SUITE_RESULT_HEADER) ||
        throw(ArgumentError("$header_path has columns $header"))
    col(name) = findfirst(==(name), header)
    tags = Set(String(c[1]) for c in PSV.SUITE_CASES)
    for f in readdir(suite_dir)
        mt = match(r"^result_(.*)\.csv$", f)
        (mt === nothing || mt.captures[1] in tags) || throw(
            ArgumentError(
                "$(joinpath(suite_dir, f)) has the tag '$(mt.captures[1])', not in SUITE_CASES",
            ),
        )
    end
    out = NamedTuple[]
    for (tag, nx, tf, phi, n, k_air) in PSV.SUITE_CASES
        path = joinpath(suite_dir, "result_$(tag).csv")
        isfile(path) || throw(ArgumentError("missing $path"))
        rows = rows_of(path)
        (length(rows) == 1 && length(rows[1]) == length(header)) ||
            throw(ArgumentError("$path must hold 1 row of $(length(header)) fields"))
        row = rows[1]
        num(name) = parse(Float64, row[col(name)])
        k_rec = num("k_air")
        matches =
            row[col("tag")] == tag &&
            parse(Int, row[col("Nx")]) == nx &&
            parse(Int, row[col("n_steps")]) == n &&
            isapprox(num("t_form_myr"), tf; rtol=1.0e-12) &&
            isapprox(num("phi_init"), phi; rtol=1.0e-12) &&
            isapprox(num("dx"), PSV.L_BOX / (nx - 1); rtol=1.0e-12) &&
            (k_air === nothing ? k_rec < 0.0 : isapprox(k_rec, k_air; rtol=1.0e-12))
        matches || throw(ArgumentError("$path does not match the case $tag of SUITE_CASES"))
        e_res, e_lim = num("energy_res"), PSV.ENERGY_RESIDUAL_LIMIT
        (isfinite(e_res) && 0.0 <= e_res <= e_lim) ||
            throw(ArgumentError("$path has energy_res $e_res outside [0, $e_lim]"))
        d_ref, d_pin = num("d_ref"), PSV.suite_reference_depth(tag)
        isapprox(d_ref, d_pin; rtol=1.0e-12) || throw(
            ArgumentError("$path has the 1-D depth $d_ref m, not the pinned $d_pin m")
        )
        push!(
            out,
            (;
                tag,
                nx,
                t_form_myr=tf,
                phi_init=phi,
                n_steps=n,
                k_air,
                d_shell=num("d_shell"),
                d_ref=num("d_ref"),
                tmax_rel_l2=num("tmax_rel_l2"),
            ),
        )
    end
    return out
end

"""
    suite_metric_specs(results) -> Vector{NamedTuple}

Return the record entry of every suite metric, in the order of `suite_metric_ids`: per run
the shell depth against the 1-D depth with the tolerance `max(dx, 0.15 d_ref)`, and for
`Nx >= SUITE_TMAX_NX_MIN` the relative L2 error of the peak temperature, at most 0.10;
the units come from `metric_units`.

Parameters
----------
results : Vector{NamedTuple}
    Result rows from `read_suite_results`.

Returns
-------
Vector{NamedTuple}
    One entry per metric with `id`, `case`, `description`, `nx`, `value`, `reference`,
    `atol`, `rtol`, `comparison` and `units`.

Raises
------
ErrorException
    If `PorousShellValidation` is not loaded.
"""
function suite_metric_specs(results)
    PSV = porous_shell_validation_module()
    specs = NamedTuple[]
    for r in results
        air = r.k_air === nothing ? "matched air" : @sprintf("air k = %g W/(m K)", r.k_air)
        run = "t_form $(r.t_form_myr) Myr, phi_init $(r.phi_init), $(air), $(r.n_steps) steps of 10 kyr"
        entries = [(
            "shell_depth_$(r.tag)",
            "Median shell depth against the 1-D shell depth, $run",
            r.d_shell,
            r.d_ref,
        )]
        r.nx >= PSV.SUITE_TMAX_NX_MIN && push!(
            entries,
            (
                "tmax_rel_l2_$(r.tag)",
                "Relative L2 error of the binned peak temperature rise, $run",
                r.tmax_rel_l2,
                0.0,
            ),
        )
        for (id, description, value, reference) in entries
            lim = PSV.suite_metric_limit(id, r.nx, reference)
            push!(
                specs,
                (;
                    id,
                    case=r.tag,
                    description,
                    nx=r.nx,
                    value,
                    reference,
                    atol=lim.atol,
                    rtol=lim.rtol,
                    comparison=lim.comparison,
                    units=PSV.metric_units(id),
                ),
            )
        end
    end
    return specs
end

"""
    suite_record_text(results; src_commit) -> String

Return the TOML text of the suite record for the result rows `results`.

Parameters
----------
results : Vector{NamedTuple}
    Result rows from `read_suite_results`.
src_commit : String
    40-character lowercase hex SHA-1 of the last commit that modified `src/`.

Returns
-------
String
    The record text.

Raises
------
ErrorException
    If `PorousShellValidation` is not loaded.
"""
function suite_record_text(results; src_commit::String)::String
    PSV = porous_shell_validation_module()
    io = IOBuffer()
    println(
        io,
        "# Validation record of the full porous shell runs: Erebus against the 1-D reference to the shell freeze.",
    )
    println(
        io,
        "# Generated by tools/export_validation_porous_shell_suite.jl from the result files of scripts/run_validation_porous_shell_suite.jl.",
    )
    println(
        io,
        "# metadata.src_commit records the latest commit that modified src/ (git log -n 1 --format=%H -- src).",
    )
    println(io)
    println(io, "[metadata]")
    println(io, "id = \"$(PSV.SUITE_RECORD_ID)\"")
    println(io, "title = \"Porous Shell Suite Validation\"")
    println(io, "page = \"validation/porous_shell.md\"")
    println(io, "src_commit = \"$(src_commit)\"")
    println(
        io, "timestamp_utc = \"$(Dates.format(now(Dates.UTC), "yyyy-mm-ddTHH:MM:SSZ"))\""
    )
    println(io, "julia_version = \"$(VERSION)\"")
    println(io, "config_sha256 = \"$(PSV.suite_config_sha256())\"")
    println(io, "reference_type = \"1d_fv\"")
    println(
        io,
        "reference_citation = \"One-dimensional spherical finite-volume solution of each run with the same time step and number of steps (solve_porous_sphere in benchmarks/validation_porous_shell.jl)\"",
    )
    println(io)
    println(io, "[config]")
    summary = PSV.expected_suite_config_summary()
    for key in PSV.REQUIRED_CONFIG_KEYS[PSV.SUITE_RECORD_ID]
        v = summary[key]
        text = if v isa AbstractVector{<:AbstractString} || v isa AbstractString
            sprint(TOML.print, Dict(key => v))
        elseif v isa AbstractVector{<:AbstractFloat}
            "$key = [" * join(repr.(Float64.(v)), ", ") * "]\n"
        elseif v isa AbstractFloat
            "$key = $(repr(Float64(v)))\n"
        else
            "$key = $(v isa AbstractVector ? "[" * join(v, ", ") * "]" : v)\n"
        end
        print(io, text)
    end
    for m in suite_metric_specs(results)
        status = metric_status(
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
        println(io, "status = \"$(status)\"")
    end
    return String(take!(io))
end

"""
    write_suite_record(suite_dir; record_path=DEFAULT_SUITE_RECORD_PATH, src_commit=nothing) -> Dict{String,Any}

Read the result files of `suite_dir`, write the suite record to `record_path` and, for the
default path, its table on the documentation page. A candidate record that fails
`validate_record_dict` (a value outside its tolerance among others) is not written.

Parameters
----------
suite_dir : String
    Folder with `header.csv` and the `result_<tag>.csv` files.
record_path : String, optional
    Output path of the record (default `benchmarks/records/porous_shell_suite.toml`).
src_commit : Union{String,Nothing}, optional
    40-character lowercase hex SHA-1, or `nothing` to query `last_src_commit()`.

Returns
-------
Dict{String,Any}
    The parsed record.

Raises
------
DomainError
    If `record_path` is empty or `src_commit` is not a 40-character lowercase hex SHA-1.
ArgumentError
    From `read_suite_results`.
ErrorException
    If the candidate record fails `validate_record_dict`, or `PorousShellValidation` is not
    loaded.
"""
function write_suite_record(
    suite_dir::String;
    record_path::String=DEFAULT_SUITE_RECORD_PATH,
    src_commit::Union{String,Nothing}=nothing,
)
    PSV = porous_shell_validation_module()
    isempty(record_path) && throw(DomainError(record_path, "Record path must be non-empty"))
    commit = something(src_commit, last_src_commit())
    occursin(r"^[0-9a-f]{40}$", commit) || throw(
        DomainError(
            commit, "src_commit must be a 40-character lowercase hexadecimal SHA-1"
        ),
    )
    text = suite_record_text(read_suite_results(suite_dir); src_commit=commit)
    data = TOML.parse(text)
    ok, errs = PSV.validate_record_dict(data; require_doc_pages=false, check_doc_sync=false)
    ok ||
        error("Refusing to write invalid suite record to $record_path: " * join(errs, "; "))
    mkpath(dirname(record_path))
    write(record_path, text)
    if abspath(record_path) == abspath(DEFAULT_SUITE_RECORD_PATH)
        PSV.update_record_doc_block!(data) ||
            @warn "Documentation page has no suite record table to update"
    end
    return data
end

"""
    run_suite_export_cli(args::Vector{String}=ARGS; src_commit=nothing, check_clean_tree=true, dirty_status=nothing, repo_root=normpath(joinpath(@__DIR__, ".."))) -> Int

Run the command-line export of the suite record: `<suite_dir> [record_path]`. Uncommitted
changes in `src/`, `benchmarks/` or `tools/` (other than the suite record) refuse the
export, as in `run_export_cli` of `tools/export_validation_porous_shell.jl`. An error of
the export is printed and gives the exit code `1`.

Parameters
----------
args : Vector{String}, optional
    Command-line arguments: the folder of the result files and at most one record path.
src_commit : Union{String,Nothing}, optional
    Explicit `src/` commit SHA-1, or `nothing` to query `last_src_commit(repo_root)`.
check_clean_tree : Bool, optional
    Refuse to export when `src/`, `benchmarks/` or `tools/` has uncommitted changes
    (default `true`).
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
InterruptException
    If the run is interrupted.
"""
function run_suite_export_cli(
    args::Vector{String}=ARGS;
    src_commit::Union{String,Nothing}=nothing,
    check_clean_tree::Bool=true,
    dirty_status::Union{String,Nothing}=nothing,
    repo_root::String=normpath(joinpath(@__DIR__, "..")),
)::Int
    if !(1 <= length(args) <= 2) || any(startswith("-"), args)
        println(
            stderr,
            "Usage: julia --project=. tools/export_validation_porous_shell_suite.jl <suite_dir> [record_path]",
        )
        return 2
    end
    record_path = length(args) == 2 ? args[2] : DEFAULT_SUITE_RECORD_PATH
    try
        if check_clean_tree
            ok_clean, dirty_out = if dirty_status === nothing
                check_clean_source_tree(repo_root; record=SUITE_RECORD_RELPATH)
            else
                (isempty(strip(dirty_status)), String(strip(dirty_status)))
            end
            if !ok_clean
                println(
                    stderr,
                    "Refusing to export the suite record with uncommitted changes in src/, benchmarks/ or tools/:\n$dirty_out",
                )
                return 1
            end
        end
        commit = something(src_commit, last_src_commit(repo_root))
        write_suite_record(args[1]; record_path, src_commit=commit)
        println("Wrote $record_path")
        return 0
    catch err
        err isa InterruptException && rethrow()
        println(stderr, sprint(showerror, err))
        return 1
    end
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "..", "benchmarks", "validation_porous_shell.jl"))
    exit(Base.invokelatest(PorousShellSuiteExport.run_suite_export_cli, ARGS))
end

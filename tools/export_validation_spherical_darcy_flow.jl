#!/usr/bin/env julia
# Export the spherical Darcy flow validation record TOML from live benchmark measurements.

module SphericalDarcyFlowExport

using Dates
using Printf
using TOML
using Erebus

const DEFAULT_EXPORT_RECORD_PATH = normpath(
    joinpath(@__DIR__, "..", "benchmarks", "records", "spherical_darcy_flow.toml")
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
    check_clean_source_tree(repo_root::String=normpath(joinpath(@__DIR__, ".."))) -> Tuple{Bool,String}

Verify that `src/` and `benchmarks/` have no uncommitted modifications before exporting a validation record.
The record file itself is excluded.

Parameters
----------
repo_root : String, optional
    Path to the repository root directory.

Returns
-------
Tuple{Bool,String}
    `true` and empty string when `src/` and `benchmarks/` are clean; `false` and `git status --porcelain` output otherwise.
"""
function check_clean_source_tree(
    repo_root::String=normpath(joinpath(@__DIR__, ".."))
)::Tuple{Bool,String}
    dirty = strip(
        readchomp(
            `git -C $repo_root status --porcelain -- src benchmarks ":(exclude)benchmarks/records/spherical_darcy_flow.toml"`,
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
        reference == 0.0 && isfinite(value) && 0.0 <= value <= atol
    else
        false
    end
    return passed ? "PASS" : "FAIL"
end

"""
    metric_specs(metrics, refs) -> Vector{NamedTuple}

Return the record entry of every metric, in the order of `EXPECTED_METRIC_IDS`.

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
"""
function metric_specs(metrics::Dict{String,Float64}, refs::Dict{String,Float64})
    specs = NamedTuple[]
    function add!(id, case, description, nx, atol, comparison, units)
        return push!(
            specs,
            (;
                id,
                case,
                description,
                nx,
                value=metrics[id],
                reference=refs[id],
                atol,
                rtol=0.0,
                comparison,
                units,
            ),
        )
    end
    l2_atol = Dict(33 => 0.20, 65 => 0.08, 129 => 0.05)
    for (tag, label) in (("off", "2-D disk"), ("on", "3-D sphere"))
        for nx in (33, 65, 129)
            add!(
                "profile_l2_inner_$(tag)_nx$(nx)",
                "profile",
                "Excess pressure L2 relative error for r <= 0.9 R, $(label) (Nx=$(nx))",
                nx,
                l2_atol[nx],
                "le_atol",
                "dimensionless",
            )
        end
    end
    for (tag, label) in (("off", "2-D disk"), ("on", "3-D sphere"))
        add!(
            "profile_order_$(tag)_65_129",
            "profile",
            "Convergence order log2(L2(Nx=65) / L2(Nx=129)), $(label)",
            129,
            0.3,
            "isapprox",
            "dimensionless",
        )
    end
    add!(
        "centre_ratio_on_off_nx65",
        "profile",
        "Centre excess pressure ratio, sphere to disk (Nx=65)",
        65,
        0.01,
        "isapprox",
        "dimensionless",
    )
    for (tag, label) in (("off", "2-D disk"), ("on", "3-D sphere"))
        add!(
            "vent_closure_$(tag)_nx65",
            "closure",
            "Vented to produced fluid rate, abs(ratio - 1), $(label) (Nx=65)",
            65,
            1.0e-8,
            "le_atol",
            "dimensionless",
        )
    end
    add!(
        "solid_speed_max_nx65",
        "closure",
        "Largest solid velocity component, both settings (Nx=65)",
        65,
        1.0e-12,
        "le_atol",
        "m/s",
    )
    return specs
end

"""
    write_spherical_darcy_flow_record(record_path::String=DEFAULT_EXPORT_RECORD_PATH; src_commit=nothing, bench=nothing, accept_physics_change=false) -> NamedTuple

Run the spherical Darcy flow benchmark (or use a pre-evaluated `bench`), validate the candidate
record with `validate_record_dict`, verify that no metric leaves its physics regression band around
the existing `record_path` (unless `accept_physics_change=true`), and write `record_path`.
Requires `SphericalDarcyFlowValidation` (`benchmarks/validation_spherical_darcy_flow.jl`) in `Main`.

Parameters
----------
record_path : String, optional
    Output path for the TOML validation record.
src_commit : Union{String,Nothing}, optional
    40-character lowercase hex SHA-1 of the last commit modifying `src/`. When `nothing`, queries `last_src_commit()`.
bench : Any, optional
    Result of `SphericalDarcyFlowValidation.evaluate_spherical_darcy_flow_benchmark()`.
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
    If a metric leaves its regression band around the existing record without `accept_physics_change=true`, or if `validate_record_dict` fails.
"""
function write_spherical_darcy_flow_record(
    record_path::String=DEFAULT_EXPORT_RECORD_PATH;
    src_commit::Union{String,Nothing}=nothing,
    bench=nothing,
    accept_physics_change::Bool=false,
)
    SDF = Main.SphericalDarcyFlowValidation
    isempty(record_path) && throw(DomainError(record_path, "Record path must be non-empty"))
    commit_str = src_commit === nothing ? last_src_commit() : src_commit
    !occursin(r"^[0-9a-f]{40}$", commit_str) && throw(
        DomainError(
            commit_str, "src_commit must be a 40-character lowercase hexadecimal SHA-1"
        ),
    )

    bench_res = if bench === nothing
        SDF.evaluate_spherical_darcy_flow_benchmark()
    else
        bench
    end

    if isfile(record_path) && !accept_physics_change
        ok_bands, band_errs = SDF.check_physics_regression_bands(
            bench_res.metrics, record_path
        )
        if !ok_bands
            error(
                "Refusing to overwrite $record_path without --accept-physics-change: " *
                join(band_errs, "; "),
            )
        end
    end

    cfg = SDF.make_config(SDF.NX_REFERENCE; spherical_flow=true)
    cfg_sha = SDF.canonical_config_sha256(cfg)
    specs = metric_specs(bench_res.metrics, SDF.expected_metric_references())
    ts_utc = Dates.format(now(Dates.UTC), "yyyy-mm-ddTHH:MM:SSZ")

    io = IOBuffer()
    println(
        io, "# Validation record for steady Darcy flow from a uniform source in a sphere."
    )
    println(io, "# Generated by tools/export_validation_spherical_darcy_flow.jl.")
    println(
        io,
        "# metadata.src_commit records the latest commit that modified src/ (git log -n 1 --format=%H -- src).",
    )
    println(io)
    println(io, "[metadata]")
    println(io, "id = \"spherical_darcy_flow\"")
    println(io, "title = \"Spherical Darcy Flow Validation\"")
    println(io, "page = \"validation/spherical_darcy_flow.md\"")
    println(io, "src_commit = \"$(commit_str)\"")
    println(io, "timestamp_utc = \"$(ts_utc)\"")
    println(io, "julia_version = \"$(VERSION)\"")
    println(io, "config_sha256 = \"$(cfg_sha)\"")
    println(io, "reference_type = \"analytical\"")
    println(
        io,
        "reference_citation = \"Exact steady Darcy pressure of a uniform fluid source in a sphere (d = 3) and a disk (d = 2) with fixed surface pressure\"",
    )
    println(io)
    cfg_sum = SDF.expected_config_summary(cfg)
    println(io, "[config]")
    println(io, "Nx = $(cfg_sum["Nx"])")
    for key in (
        "R_planet_m",
        "L_box_m",
        "mu_f_Pa_s",
        "k_m2",
        "S_source_per_s",
        "phi",
        "eta_solid_Pa_s",
        "k_vent_m2",
        "T_surf_K",
        "P_amb_Pa",
        "inner_fraction",
    )
        println(io, "$(key) = $(repr(Float64(cfg_sum[key])))")
    end
    println(io, "spherical_metric = $(cfg_sum["spherical_metric"])")
    println(io, "darcy_elimination = $(cfg_sum["darcy_elimination"])")

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
    ok_valid, valid_errs = SDF.validate_record_dict(
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
        SDF.update_record_doc_block!(candidate_dict) ||
            @warn "Documentation page has no record table to update"
    end
    return bench_res
end

"""
    run_export_cli(args::Vector{String}=ARGS; src_commit=nothing, bench=nothing, check_clean_tree=true, dirty_status=nothing, repo_root=normpath(joinpath(@__DIR__, ".."))) -> Int

Run the command-line export tool for the spherical Darcy flow validation record. Requires
`SphericalDarcyFlowValidation` (`benchmarks/validation_spherical_darcy_flow.jl`) to be loaded
in `Main` before calling.

Parameters
----------
args : Vector{String}, optional
    Command-line arguments (`[--accept-physics-change]` and optional `[record_path]`).
src_commit : Union{String,Nothing}, optional
    Explicit `src/` commit SHA-1, or `nothing` to query `last_src_commit(repo_root)`.
bench : Any, optional
    Pre-evaluated benchmark tuple, or `nothing` to evaluate on demand.
check_clean_tree : Bool, optional
    Refuse to export when `src/` or `benchmarks/` has uncommitted changes (default `true`).
dirty_status : Union{String,Nothing}, optional
    Replacement for the `git status --porcelain` output of the clean-tree check.
repo_root : String, optional
    Repository root directory checked for `src/` and `benchmarks/` cleanliness.

Returns
-------
Int
    Exit code (`0` on success, `1` on export refusal, `2` on invalid CLI flag).

Raises
------
DomainError
    If the record path is empty or `src_commit` is not a 40-character lowercase hex SHA-1.
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
    record_path = DEFAULT_EXPORT_RECORD_PATH
    for arg in args
        if arg == "--accept-physics-change"
            accept_physics_change = true
        elseif startswith(arg, "--")
            println(
                stderr,
                "Usage: julia --project=. tools/export_validation_spherical_darcy_flow.jl [--accept-physics-change] [record_path]",
            )
            return 2
        else
            record_path = arg
        end
    end

    if check_clean_tree
        ok_clean, dirty_out = if dirty_status === nothing
            check_clean_source_tree(repo_root)
        else
            (isempty(strip(dirty_status)), String(strip(dirty_status)))
        end
        if !ok_clean
            println(
                stderr,
                "Refusing to export validation record with uncommitted changes in src/ or benchmarks/:\n$dirty_out",
            )
            return 1
        end
    end

    try
        commit_val = src_commit === nothing ? last_src_commit(repo_root) : src_commit
        write_spherical_darcy_flow_record(
            record_path;
            src_commit=commit_val,
            bench=bench,
            accept_physics_change=accept_physics_change,
        )
        println("Wrote $record_path")
        return 0
    catch err
        err isa DomainError && rethrow(err)
        println(stderr, sprint(showerror, err))
        return 1
    end
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_darcy_flow.jl"))
    exit(Base.invokelatest(SphericalDarcyFlowExport.run_export_cli, ARGS))
end

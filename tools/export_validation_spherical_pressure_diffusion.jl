#!/usr/bin/env julia
# Export the spherical pressure diffusion validation record TOML from live benchmark measurements.

module SphericalPressureDiffusionExport

using Dates
using Printf
using TOML

const DEFAULT_EXPORT_RECORD_PATH = normpath(
    joinpath(@__DIR__, "..", "benchmarks", "records", "spherical_pressure_diffusion.toml")
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

Verify that `src/`, `benchmarks/` and `tools/` have no uncommitted modifications before exporting a validation record.
The record file itself is excluded.

Parameters
----------
repo_root : String, optional
    Path to the repository root directory.

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
    repo_root::String=normpath(joinpath(@__DIR__, ".."))
)::Tuple{Bool,String}
    dirty = strip(
        readchomp(
            `git -C $repo_root status --porcelain -- src benchmarks tools ":(exclude)benchmarks/records/spherical_pressure_diffusion.toml"`,
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

Return the record entry of every metric, in the order of `EXPECTED_METRIC_IDS`. The bounds
of the centre and profile errors and of the vent pressure are about 1.4 times the values
measured on the committed solver. The solid velocity metrics are in m/s; every other
metric is dimensionless.

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
"""
function metric_specs(metrics::Dict{String,Float64}, refs::Dict{String,Float64})
    specs = NamedTuple[]
    function add!(
        id, case, description, nx, atol, comparison; rtol=0.0, units="dimensionless"
    )
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
                rtol,
                comparison,
                units,
            ),
        )
    end
    labels = (("off", "2-D disk"), ("on", "3-D sphere"))
    # Bounds per row prefix, setting and grid
    bounds = Dict(
        "centre_err_m5" => Dict(
            "off" => Dict(33 => 0.0027, 65 => 7.3e-4, 129 => 3.2e-4),
            "on" => Dict(33 => 0.0064, 65 => 0.0018, 129 => 7.9e-4),
        ),
        "centre_err_m20" => Dict(
            "off" => Dict(33 => 0.069, 65 => 0.027, 129 => 0.014),
            "on" => Dict(33 => 0.12, 65 => 0.048, 129 => 0.025),
        ),
        "l2_m20" => Dict(
            "off" => Dict(33 => 0.18, 65 => 0.075, 129 => 0.04),
            "on" => Dict(33 => 0.25, 65 => 0.108, 129 => 0.058),
        ),
    )
    error_rows = (
        ("centre_err_m5", "centre", "abs(p / p_BE - 1) at the centre cell, m = 5"),
        ("centre_err_m20", "centre", "abs(p / p_BE - 1) at the centre cell, m = 20"),
        (
            "l2_m20",
            "profile",
            "Fluid overpressure L2 relative error for r <= 0.9 R, m = 20",
        ),
    )
    for (prefix, case, text) in error_rows, (tag, label) in labels, nx in (33, 65, 129)
        add!(
            "$(prefix)_$(tag)_nx$(nx)",
            case,
            "$(text), $(label) (Nx=$(nx))",
            nx,
            bounds[prefix][tag][nx],
            "le_atol",
        )
    end
    for (short, case, text) in
        (("centre", "centre", "centre error"), ("l2", "profile", "L2 profile error")),
        (tag, label) in labels

        add!(
            "$(short)_order_m20_$(tag)_65_129",
            case,
            "Convergence order log2(e(Nx=65) / e(Nx=129)) of the $(text) at m = 20, $(label)",
            129,
            0.35,
            "isapprox",
        )
    end
    for (nx, rtol) in ((65, 0.022), (129, 0.012))
        add!(
            "centre_ratio_m20_nx$(nx)",
            "centre",
            "Centre fluid overpressure ratio, sphere to disk, m = 20 (Nx=$(nx))",
            nx,
            0.0,
            "isapprox";
            rtol=rtol,
        )
    end
    for (tag, label) in labels, nx in (33, 65, 129)
        add!(
            "pt_kbw_$(tag)_nx$(nx)",
            "storage",
            "Largest abs((Pt - Pt0) - kbw (Pf - Pf0)) / P0 over the rock, m = 20, $(label) (Nx=$(nx))",
            nx,
            1.0e-9,
            "le_atol",
        )
    end
    for (tag, label) in labels
        add!(
            "s_invariance_$(tag)_nx33",
            "storage",
            "Largest difference of (Pf - P_amb) / P0 between bs = 1e-10 and 1e-11 Pa^-1, m = 20, $(label) (Nx=33)",
            33,
            1.0e-9,
            "le_atol",
        )
    end
    for (tag, label) in labels
        add!(
            "elim_$(tag)_nx33",
            "solver",
            "Largest relative difference of pf, pr, qxD, qyD between the 6- and 4-variable paths, m = 20, $(label) (Nx=33)",
            33,
            1.0e-10,
            "le_atol",
        )
    end
    for (tag, label, atol) in (("off", "2-D disk", 1.4e-3), ("on", "3-D sphere", 9.0e-4))
        add!(
            "vent_dev_$(tag)_nx65",
            "vent",
            "Largest abs(Pf - P_amb) / P0 over the vent cells, m = 20, $(label) (Nx=65)",
            65,
            atol,
            "le_atol",
        )
    end
    for (tag, label) in labels
        add!(
            "vent_closure_$(tag)_nx65",
            "vent",
            "Vented to stored fluid volume rate, abs(ratio - 1), m = 20, $(label) (Nx=65)",
            65,
            1.0e-8,
            "le_atol",
        )
    end
    for (tag, label) in labels
        add!(
            "v_max_$(tag)_nx65",
            "rigid",
            "Largest solid velocity component [m/s] over the 20 steps, $(label) (Nx=65)",
            65,
            1.0e-15,
            "le_atol";
            units="m/s",
        )
    end
    return specs
end

"""
    diffusion_validation_module(hosts=(parentmodule(@__MODULE__), Main)) -> Module

Return the `SphericalPressureDiffusionValidation` module (`benchmarks/validation_spherical_pressure_diffusion.jl`)
of the first module in `hosts` that defines it: by default the module that includes this
file, or `Main` when that module has none.

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
    If no module of `hosts` defines `SphericalPressureDiffusionValidation`.
"""
function diffusion_validation_module(hosts::Tuple=(parentmodule(@__MODULE__), Main))::Module
    for host in hosts
        if Base.invokelatest(isdefined, host, :SphericalPressureDiffusionValidation)
            return Base.invokelatest(getfield, host, :SphericalPressureDiffusionValidation)
        end
    end
    return error(
        "SphericalPressureDiffusionValidation is not loaded; include benchmarks/validation_spherical_pressure_diffusion.jl first",
    )
end

"""
    write_spherical_pressure_diffusion_record(record_path::String=DEFAULT_EXPORT_RECORD_PATH; src_commit=nothing, bench=nothing, accept_physics_change=false) -> NamedTuple

Run the spherical pressure diffusion benchmark (or use a pre-evaluated `bench`), validate the candidate
record with `validate_record_dict`, verify that no metric leaves its physics regression band around
the existing `record_path` (unless `accept_physics_change=true`), and write `record_path`.
Requires `SphericalPressureDiffusionValidation` (`benchmarks/validation_spherical_pressure_diffusion.jl`), see
`diffusion_validation_module`.

Parameters
----------
record_path : String, optional
    Output path for the TOML validation record.
src_commit : Union{String,Nothing}, optional
    40-character lowercase hex SHA-1 of the last commit modifying `src/`. When `nothing`, queries `last_src_commit()`.
bench : Any, optional
    Result of `SphericalPressureDiffusionValidation.evaluate_spherical_pressure_diffusion_benchmark()`.
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
    If a metric leaves its regression band around the existing record without `accept_physics_change=true`, if `validate_record_dict` fails, or if `SphericalPressureDiffusionValidation` is not loaded.
KeyError
    If `bench.metrics` lacks one of the metric ids.
"""
function write_spherical_pressure_diffusion_record(
    record_path::String=DEFAULT_EXPORT_RECORD_PATH;
    src_commit::Union{String,Nothing}=nothing,
    bench=nothing,
    accept_physics_change::Bool=false,
)
    SPDV = diffusion_validation_module()
    isempty(record_path) && throw(DomainError(record_path, "Record path must be non-empty"))
    commit_str = src_commit === nothing ? last_src_commit() : src_commit
    !occursin(r"^[0-9a-f]{40}$", commit_str) && throw(
        DomainError(
            commit_str, "src_commit must be a 40-character lowercase hexadecimal SHA-1"
        ),
    )

    bench_res = if bench === nothing
        SPDV.evaluate_spherical_pressure_diffusion_benchmark()
    else
        bench
    end

    if isfile(record_path) && !accept_physics_change
        ok_bands, band_errs = SPDV.check_physics_regression_bands(
            bench_res.metrics, record_path
        )
        if !ok_bands
            error(
                "Refusing to overwrite $record_path without --accept-physics-change: " *
                join(band_errs, "; "),
            )
        end
    end

    cfg = SPDV.make_config(SPDV.NX_REFERENCE; spherical_flow=true)
    cfg_sha = SPDV.canonical_config_sha256(cfg)
    specs = metric_specs(bench_res.metrics, SPDV.expected_metric_references())
    ts_utc = Dates.format(now(Dates.UTC), "yyyy-mm-ddTHH:MM:SSZ")

    io = IOBuffer()
    println(
        io,
        "# Validation record for transient pore pressure diffusion in a rigid porous sphere with a vented surface.",
    )
    println(io, "# Generated by tools/export_validation_spherical_pressure_diffusion.jl.")
    println(
        io,
        "# metadata.src_commit records the latest commit that modified src/ (git log -n 1 --format=%H -- src).",
    )
    println(io)
    println(io, "[metadata]")
    println(io, "id = \"spherical_pressure_diffusion\"")
    println(io, "title = \"Spherical Pressure Diffusion Validation\"")
    println(io, "page = \"validation/spherical_pressure_diffusion.md\"")
    println(io, "src_commit = \"$(commit_str)\"")
    println(io, "timestamp_utc = \"$(ts_utc)\"")
    println(io, "julia_version = \"$(VERSION)\"")
    println(io, "config_sha256 = \"$(cfg_sha)\"")
    println(io, "reference_type = \"analytical\"")
    println(
        io,
        "reference_citation = \"Backward-Euler Fourier series (sphere, d = 3) and Fourier-Bessel series (disk, d = 2) of pore pressure diffusion with a fixed surface pressure, rigid-matrix storage coefficient S = phi bf + (1 - phi) bphi bs / (bphi + bs)\"",
    )
    println(io)
    cfg_sum = SPDV.expected_config_summary(cfg)
    println(io, "[config]")
    println(io, "Nx = $(cfg_sum["Nx"])")
    for key in (
        "R_planet_m",
        "L_box_m",
        "mu_f_Pa_s",
        "k_m2",
        "phi",
        "beta_phi_per_Pa",
        "beta_s_per_Pa",
        "beta_s_alt_per_Pa",
        "beta_f_per_Pa",
        "storage_per_Pa",
        "kbw",
        "P0_Pa",
        "P_amb_Pa",
        "k_vent_m2",
        "T_surf_K",
        "eta_rock_Pa_s",
        "shear_modulus_Pa",
        "tau_end",
        "dt_s",
        "inner_fraction",
    )
        println(io, "$(key) = $(repr(Float64(cfg_sum[key])))")
    end
    println(io, "n_steps = $(cfg_sum["n_steps"])")
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
    ok_valid, valid_errs = SPDV.validate_record_dict(
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
        SPDV.update_record_doc_block!(candidate_dict) ||
            @warn "Documentation page has no record table to update"
    end
    return bench_res
end

"""
    run_export_cli(args::Vector{String}=ARGS; src_commit=nothing, bench=nothing, check_clean_tree=true, dirty_status=nothing, repo_root=normpath(joinpath(@__DIR__, ".."))) -> Int

Run the command-line export tool for the spherical pressure diffusion validation record. Requires
`SphericalPressureDiffusionValidation` (`benchmarks/validation_spherical_pressure_diffusion.jl`), see
`diffusion_validation_module`. An argument other than `--accept-physics-change` that
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
                "Usage: julia --project=. tools/export_validation_spherical_pressure_diffusion.jl [--accept-physics-change] [record_path]",
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
        write_spherical_pressure_diffusion_record(
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
    include(
        joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_pressure_diffusion.jl")
    )
    exit(Base.invokelatest(SphericalPressureDiffusionExport.export_main, ARGS))
end

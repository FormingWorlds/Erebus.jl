#!/usr/bin/env julia
# Export the spherical conduction validation record TOML from live benchmark measurements.

using Dates
using Printf
using TOML
using Erebus

const DEFAULT_EXPORT_RECORD_PATH = normpath(
    joinpath(@__DIR__, "..", "benchmarks", "records", "spherical_conduction.toml")
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
            `git -C $repo_root status --porcelain -- src benchmarks ":(exclude)benchmarks/records/spherical_conduction.toml"`,
        ),
    )
    return (isempty(dirty), String(dirty))
end

function _compute_metric_status(
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
    write_spherical_conduction_record(record_path::String=DEFAULT_EXPORT_RECORD_PATH; src_commit=nothing, bench=nothing, accept_physics_change=false) -> NamedTuple

Run the spherical conduction benchmark suite (or use a pre-evaluated `bench`), validate the candidate
record with `validate_record_dict`, verify that no metric leaves its physics regression band around
the existing `record_path` (unless `accept_physics_change=true`), and write `record_path`.

Parameters
----------
record_path : String, optional
    Output path for the TOML validation record.
src_commit : Union{String,Nothing}, optional
    40-character lowercase hex SHA-1 of the last commit modifying `src/`. When `nothing`, queries `last_src_commit()`.
bench : Any, optional
    Pre-evaluated benchmark result from `SphericalConductionValidation.evaluate_spherical_conduction_benchmark()`.
accept_physics_change : Bool, optional
    Allow overwriting an existing record when a metric moves outside `PHYSICS_REGRESSION_BANDS` around the existing value (default `false`).

Returns
-------
NamedTuple
    Evaluated benchmark summary tuple.

Raises
------
DomainError
    If `isempty(record_path)` or `src_commit` is not a 40-character lowercase hex SHA-1.
ErrorException
    If a metric leaves its regression band around the existing record without `accept_physics_change=true`, or if `validate_record_dict` fails.
"""
function write_spherical_conduction_record(
    record_path::String=DEFAULT_EXPORT_RECORD_PATH;
    src_commit::Union{String,Nothing}=nothing,
    bench=nothing,
    accept_physics_change::Bool=false,
)
    isempty(record_path) && throw(DomainError(record_path, "Record path must be non-empty"))
    commit_str = src_commit === nothing ? last_src_commit() : src_commit
    !occursin(r"^[0-9a-f]{40}$", commit_str) && throw(
        DomainError(
            commit_str, "src_commit must be a 40-character lowercase hexadecimal SHA-1"
        ),
    )

    bench_res = if bench === nothing
        SphericalConductionValidation.evaluate_spherical_conduction_benchmark()
    else
        bench
    end

    if isfile(record_path) && !accept_physics_change
        ok_bands, band_errs = SphericalConductionValidation.check_physics_regression_bands(
            bench_res.metrics, record_path
        )
        if !ok_bands
            error(
                "Refusing to overwrite $record_path without --accept-physics-change: " *
                join(band_errs, "; "),
            )
        end
    end

    cfg65 = SphericalConductionValidation.make_spherical_conduction_config(;
        Nx=SphericalConductionValidation.NX_FINE
    )
    cfg_sha = SphericalConductionValidation.canonical_config_sha256(cfg65)
    exp_refs = SphericalConductionValidation.expected_metric_references(bench_res)
    ts_utc = Dates.format(now(Dates.UTC), "yyyy-mm-ddTHH:MM:SSZ")

    metric_specs = [
        (
            id="steady_center_dT_ratio_3d_2d_nx33",
            case="steady",
            description="Steady 3D-to-2D centre warming ratio (Nx=33)",
            nx=33,
            value=bench_res.metrics["steady_center_dT_ratio_3d_2d_nx33"],
            reference=exp_refs["steady_center_dT_ratio_3d_2d_nx33"],
            atol=0.01,
            rtol=0.015,
            comparison="isapprox",
            units="dimensionless",
        ),
        (
            id="steady_center_dT_ratio_3d_2d_nx65",
            case="steady",
            description="Steady 3D-to-2D centre warming ratio (Nx=65)",
            nx=65,
            value=bench_res.metrics["steady_center_dT_ratio_3d_2d_nx65"],
            reference=exp_refs["steady_center_dT_ratio_3d_2d_nx65"],
            atol=0.01,
            rtol=0.015,
            comparison="isapprox",
            units="dimensionless",
        ),
        (
            id="steady_profile_l2_rel_error_nx33",
            case="steady",
            description="Steady radial profile L2 relative error (r <= R, Nx=33)",
            nx=33,
            value=bench_res.metrics["steady_profile_l2_rel_error_nx33"],
            reference=exp_refs["steady_profile_l2_rel_error_nx33"],
            atol=0.15,
            rtol=0.0,
            comparison="le_atol",
            units="dimensionless",
        ),
        (
            id="steady_profile_l2_rel_error_nx65",
            case="steady",
            description="Steady radial profile L2 relative error (r <= R, Nx=65)",
            nx=65,
            value=bench_res.metrics["steady_profile_l2_rel_error_nx65"],
            reference=exp_refs["steady_profile_l2_rel_error_nx65"],
            atol=0.08,
            rtol=0.0,
            comparison="le_atol",
            units="dimensionless",
        ),
        (
            id="steady_energy_conservation_residual_nx65",
            case="steady",
            description="Steady discrete energy conservation residual (Nx=65)",
            nx=65,
            value=bench_res.metrics["steady_energy_conservation_residual_nx65"],
            reference=exp_refs["steady_energy_conservation_residual_nx65"],
            atol=1.0e-11,
            rtol=0.0,
            comparison="le_atol",
            units="dimensionless",
        ),
        (
            id="transient_be_vs_1dfv_rel_diff",
            case="transient",
            description="Step-averaged analytical BE series vs 1D spherical FV relative difference",
            nx=0,
            value=bench_res.metrics["transient_be_vs_1dfv_rel_diff"],
            reference=exp_refs["transient_be_vs_1dfv_rel_diff"],
            atol=1.0e-5,
            rtol=0.0,
            comparison="le_atol",
            units="dimensionless",
        ),
        (
            id="transient_center_dT_nx65",
            case="transient",
            description="Transient centre temperature rise dT(r_min) vs step-averaged BE reference (Nx=65)",
            nx=65,
            value=bench_res.metrics["transient_center_dT_nx65"],
            reference=exp_refs["transient_center_dT_nx65"],
            atol=25.0,
            rtol=0.05,
            comparison="isapprox",
            units="K",
        ),
        (
            id="transient_profile_l2_rel_error_nx65",
            case="transient",
            description="Transient radial profile L2 relative error vs step-averaged BE reference (r <= R, Nx=65)",
            nx=65,
            value=bench_res.metrics["transient_profile_l2_rel_error_nx65"],
            reference=exp_refs["transient_profile_l2_rel_error_nx65"],
            atol=0.10,
            rtol=0.0,
            comparison="le_atol",
            units="dimensionless",
        ),
        (
            id="transient_profile_l2_convergence_ratio_65_33",
            case="transient",
            description="Transient radial profile L2 convergence ratio L2(Nx=65) / L2(Nx=33)",
            nx=65,
            value=bench_res.metrics["transient_profile_l2_convergence_ratio_65_33"],
            reference=exp_refs["transient_profile_l2_convergence_ratio_65_33"],
            atol=0.60,
            rtol=0.0,
            comparison="le_atol",
            units="dimensionless",
        ),
        (
            id="transient_max_air_shell_temperature_nx65",
            case="transient",
            description="Maximum sticky-air temperature in r >= R + dx (Nx=65)",
            nx=65,
            value=bench_res.metrics["transient_max_air_shell_temperature_nx65"],
            reference=exp_refs["transient_max_air_shell_temperature_nx65"],
            atol=0.01,
            rtol=0.0,
            comparison="isapprox",
            units="K",
        ),
    ]

    io = IOBuffer()
    println(
        io, "# Validation record for spherical conduction and step-averaged 26Al heating."
    )
    println(io, "# Generated by tools/export_validation_spherical_conduction.jl.")
    println(
        io,
        "# metadata.src_commit records the latest commit that modified src/ (git log -n 1 --format=%H -- src).",
    )
    println(io)
    println(io, "[metadata]")
    println(io, "id = \"spherical_conduction\"")
    println(io, "title = \"Spherical Conduction and Radiogenic Heating Validation\"")
    println(io, "page = \"validation/thermal_conduction.md\"")
    println(io, "src_commit = \"$(commit_str)\"")
    println(io, "timestamp_utc = \"$(ts_utc)\"")
    println(io, "julia_version = \"$(VERSION)\"")
    println(io, "config_sha256 = \"$(cfg_sha)\"")
    println(io, "reference_type = \"analytical_and_1d_fv\"")
    println(
        io,
        "reference_citation = \"Exact spherical Poisson solution and step-averaged backward-Euler Dirichlet eigenmode series cross-verified against a 1D spherical finite-volume solver\"",
    )
    println(io)
    cfg_sum = SphericalConductionValidation.expected_config_summary(cfg65)
    println(io, "[config]")
    println(io, "Nx_coarse = $(cfg_sum["Nx_coarse"])")
    println(io, "Nx_fine = $(cfg_sum["Nx_fine"])")
    println(io, "n_steps = $(cfg_sum["n_steps"])")
    println(io, "dt_yr = $(repr(Float64(cfg_sum["dt_yr"])))")
    println(io, "R_planet_m = $(repr(Float64(cfg_sum["R_planet_m"])))")
    println(io, "L_box_m = $(repr(Float64(cfg_sum["L_box_m"])))")
    println(io, "k_W_mK = $(repr(Float64(cfg_sum["k_W_mK"])))")
    println(io, "rho_kg_m3 = $(repr(Float64(cfg_sum["rho_kg_m3"])))")
    println(io, "cp_J_kgK = $(repr(Float64(cfg_sum["cp_J_kgK"])))")
    println(io, "T_surf_K = $(repr(Float64(cfg_sum["T_surf_K"])))")
    println(io, "H_steady_W_m3 = $(repr(Float64(cfg_sum["H_steady_W_m3"])))")
    println(io, "spherical_metric = $(cfg_sum["spherical_metric"])")
    println(io, "reg_cells = $(repr(Float64(cfg_sum["reg_cells"])))")
    println(io, "disk_enabled = $(cfg_sum["disk_enabled"])")

    for m in metric_specs
        status_str = _compute_metric_status(
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
    ok_valid, valid_errs = SphericalConductionValidation.validate_record_dict(
        candidate_dict; require_doc_pages=false
    )
    if !ok_valid
        error(
            "Refusing to write invalid validation record to $record_path: " *
            join(valid_errs, "; "),
        )
    end

    mkpath(dirname(record_path))
    write(record_path, toml_str)
    return bench_res
end

"""
    run_export_cli(args::Vector{String}=ARGS; src_commit=nothing, bench=nothing, check_clean_tree=true, repo_root=normpath(joinpath(@__DIR__, ".."))) -> Int

Run the command-line export tool for the spherical conduction validation record. Requires
`SphericalConductionValidation` (`benchmarks/validation_spherical_conduction.jl`) to be loaded
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
repo_root : String, optional
    Repository root directory checked for `src/` and `benchmarks/` cleanliness.

Returns
-------
Int
    Exit code (`0` on success, `1` on export refusal, `2` on invalid CLI flag).
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
                "Usage: julia --project=. tools/export_validation_spherical_conduction.jl [--accept-physics-change] [record_path]",
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
        write_spherical_conduction_record(
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

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_conduction.jl"))
    exit(Base.invokelatest(run_export_cli, ARGS))
end

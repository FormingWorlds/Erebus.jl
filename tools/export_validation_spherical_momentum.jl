#!/usr/bin/env julia
# Export the spherical momentum validation record TOML from live benchmark measurements.

module SphericalMomentumExport

using Dates
using Printf
using TOML

const DEFAULT_EXPORT_RECORD_PATH = normpath(
    joinpath(@__DIR__, "..", "benchmarks", "records", "spherical_momentum.toml")
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
ProcessFailedException
    If `git log` fails, for example because `repo_root` is not a git repository.
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
            `git -C $repo_root status --porcelain -- src benchmarks tools ":(exclude)benchmarks/records/spherical_momentum.toml"`,
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
of the errors are about 1.4 times the values measured on the committed solver (case `:C`:
1.4 times rounded up to 2 significant digits), the order rows accept `2 +- 0.35`, the ratio rows a relative deviation of about 1.4 times the
measured one, and the elimination rows the fixed bound `1e-7`. Every metric is
dimensionless.

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
    function add!(id, case, description, nx, atol, comparison; rtol=0.0)
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
                units="dimensionless",
            ),
        )
    end
    # Bounds per row prefix and grid
    bounds = Dict(
        "l2_a_on" => Dict(33 => 0.053, 65 => 0.0122, 129 => 0.0031),
        "l2_a_off" => Dict(33 => 0.053, 65 => 0.012, 129 => 0.003),
        "l2_b" => Dict(33 => 0.057, 65 => 0.0138, 129 => 0.0035),
        "l2_inside_b" => Dict(33 => 0.057, 65 => 0.0138, 129 => 0.0035),
        "s33_b" => Dict(33 => 0.059, 65 => 0.0152, 129 => 0.0039),
        "centre_fit_a_on" => Dict(33 => 0.05, 65 => 0.012, 129 => 0.0031),
        "centre_fit_a_off" => Dict(33 => 0.047, 65 => 0.0112, 129 => 0.0028),
        "centre_node_a_on" => Dict(33 => 0.0136, 65 => 0.01, 129 => 0.006),
        "l2_c" => Dict(33 => 0.025, 65 => 0.0062, 129 => 0.0016),
        "l2_inside_c" => Dict(33 => 0.025, 65 => 0.0062, 129 => 0.0016),
        "s33_c" => Dict(33 => 0.025, 65 => 0.006, 129 => 0.0015),
        "pt_c" => Dict(33 => 0.011, 65 => 0.0026, 129 => 0.00064),
        "pf_c" => Dict(33 => 0.024, 65 => 0.0056, 129 => 0.0015),
        "eii_c" => Dict(33 => 0.018, 65 => 0.0045, 129 => 0.0012),
        "sii_c" => Dict(33 => 0.018, 65 => 0.0045, 129 => 0.0012),
        "hs_c" => Dict(33 => 0.024, 65 => 0.0058, 129 => 0.0015),
    )
    error_rows = (
        ("l2_a_on", "A", "Velocity L2 relative error over the box, case A, option on"),
        ("l2_a_off", "A", "Velocity L2 relative error over the box, case A, option off"),
        ("l2_b", "B", "Velocity L2 relative error over the box, case B, option on"),
        (
            "l2_inside_b",
            "B",
            "Velocity L2 relative error from the nodes with r < R, case B, option on",
        ),
        (
            "s33_b",
            "B",
            "S33 L2 relative error against tau_33 = 2 eta u_r / r for r < R, case B, option on",
        ),
        (
            "centre_fit_a_on",
            "A",
            "abs(Omega_fit(0) / Omega(0) - 1), fit of degree 5 in rho^2 for 0.18 <= rho <= 0.85, case A, option on",
        ),
        (
            "centre_fit_a_off",
            "A",
            "abs(Omega_fit(0) / Omega(0) - 1), fit of degree 5 in rho^2 for 0.18 <= rho <= 0.85, case A, option off",
        ),
        (
            "centre_node_a_on",
            "A",
            "abs(Omega_node / Omega(0) - 1) at the centre basic node, case A, option on",
        ),
        ("l2_c", "C", "Velocity L2 relative error over the box, case C, option on"),
        (
            "l2_inside_c",
            "C",
            "Velocity L2 relative error from the nodes with r < R, case C, option on",
        ),
        (
            "s33_c",
            "C",
            "S33 L2 relative error against tau_33 = 2 eta (u_r / r - div3 / 3) for r < R, case C, option on",
        ),
        (
            "pt_c",
            "C",
            "Total pressure L2 error over the L2 norm of P_t - p_s for r < R, case C, option on",
        ),
        (
            "pf_c",
            "C",
            "Fluid pressure L2 error over the L2 norm of P_f - p_s for r < R, case C, option on",
        ),
        (
            "eii_c",
            "C",
            "EII L2 relative error against the exact second strain rate invariant for r < R, case C, option on",
        ),
        (
            "sii_c",
            "C",
            "SII L2 relative error against 2 eta E_II for r < R, case C, option on",
        ),
        (
            "hs_c",
            "C",
            "HS L2 relative error against the exact deviatoric, compaction and Darcy heating for r < R, case C, option on",
        ),
    )
    for (prefix, case, text) in error_rows, nx in (33, 65, 129)
        add!(
            "$(prefix)_nx$(nx)",
            case,
            "$(text) (Nx=$(nx))",
            nx,
            bounds[prefix][nx],
            "le_atol",
        )
    end
    order_rows = (
        ("l2_order_a_on", "A", "velocity L2 error, case A, option on"),
        ("l2_order_a_off", "A", "velocity L2 error, case A, option off"),
        ("l2_order_b", "B", "velocity L2 error, case B, option on"),
        ("l2_inside_order_b", "B", "velocity L2 error for r < R, case B, option on"),
        ("s33_order_b", "B", "S33 L2 error, case B, option on"),
        ("centre_fit_order_a_on", "A", "fitted centre rotation error, case A, option on"),
        ("centre_fit_order_a_off", "A", "fitted centre rotation error, case A, option off"),
        ("l2_order_c", "C", "velocity L2 error, case C, option on"),
        ("l2_inside_order_c", "C", "velocity L2 error for r < R, case C, option on"),
        ("s33_order_c", "C", "S33 L2 error, case C, option on"),
        ("pt_order_c", "C", "total pressure L2 error, case C, option on"),
        ("pf_order_c", "C", "fluid pressure L2 error, case C, option on"),
        ("eii_order_c", "C", "EII L2 error, case C, option on"),
        ("sii_order_c", "C", "SII L2 error, case C, option on"),
        ("hs_order_c", "C", "HS L2 error, case C, option on"),
    )
    for (prefix, case, text) in order_rows
        add!(
            "$(prefix)_65_129",
            case,
            "Convergence order log2(e(Nx=65) / e(Nx=129)) of the $(text)",
            129,
            0.35,
            "isapprox",
        )
    end
    for (nx, rtol) in ((65, 8.0e-4), (129, 2.5e-4))
        add!(
            "centre_fit_ratio_nx$(nx)",
            "A",
            "Fitted centre rotation ratio, option on to off, case A (Nx=$(nx))",
            nx,
            0.0,
            "isapprox";
            rtol=rtol,
        )
    end
    for (nx, rtol) in ((65, 7.5e-4), (129, 2.1e-4))
        add!(
            "profile_ratio_nx$(nx)",
            "A",
            "Velocity L2 norm ratio over r < R, option on to off, case A (Nx=$(nx))",
            nx,
            0.0,
            "isapprox";
            rtol=rtol,
        )
    end
    for (tag, case) in (("a", "A"), ("b", "B"), ("c", "C"))
        add!(
            "elim_$(tag)_nx33",
            "solver",
            "Largest relative difference of vx, vy, pr, pf between the 6- and 4-variable paths, case $(case), option on (Nx=33)",
            33,
            1.0e-7,
            "le_atol",
        )
    end
    return specs
end

"""
    momentum_validation_module(hosts=(parentmodule(@__MODULE__), Main)) -> Module

Return the `SphericalMomentumValidation` module (`benchmarks/validation_spherical_momentum.jl`)
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
    If no module of `hosts` defines `SphericalMomentumValidation`.
"""
function momentum_validation_module(hosts::Tuple=(parentmodule(@__MODULE__), Main))::Module
    for host in hosts
        if Base.invokelatest(isdefined, host, :SphericalMomentumValidation)
            return Base.invokelatest(getfield, host, :SphericalMomentumValidation)
        end
    end
    return error(
        "SphericalMomentumValidation is not loaded; include benchmarks/validation_spherical_momentum.jl first",
    )
end

"""
    write_spherical_momentum_record(record_path::String=DEFAULT_EXPORT_RECORD_PATH; src_commit=nothing, bench=nothing, accept_physics_change=false) -> NamedTuple

Run the spherical momentum benchmark (or use a pre-evaluated `bench`), validate the candidate
record with `validate_record_dict`, verify that no metric leaves its physics regression band around
the existing `record_path` (unless `accept_physics_change=true`), and write `record_path`.
Requires `SphericalMomentumValidation` (`benchmarks/validation_spherical_momentum.jl`), see
`momentum_validation_module`.

Parameters
----------
record_path : String, optional
    Output path for the TOML validation record.
src_commit : Union{String,Nothing}, optional
    40-character lowercase hex SHA-1 of the last commit modifying `src/`. When `nothing`, queries `last_src_commit()`.
bench : Any, optional
    Result of `SphericalMomentumValidation.evaluate_spherical_momentum_benchmark()`.
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
    If a metric leaves its regression band around the existing record without `accept_physics_change=true`, if `validate_record_dict` fails, or if `SphericalMomentumValidation` is not loaded.
KeyError
    If `bench.metrics` lacks one of the metric ids.
ProcessFailedException
    From `last_src_commit` when `src_commit` is `nothing` and `git log` fails.
"""
function write_spherical_momentum_record(
    record_path::String=DEFAULT_EXPORT_RECORD_PATH;
    src_commit::Union{String,Nothing}=nothing,
    bench=nothing,
    accept_physics_change::Bool=false,
)
    SMV = momentum_validation_module()
    isempty(record_path) && throw(DomainError(record_path, "Record path must be non-empty"))
    commit_str = src_commit === nothing ? last_src_commit() : src_commit
    !occursin(r"^[0-9a-f]{40}$", commit_str) && throw(
        DomainError(
            commit_str, "src_commit must be a 40-character lowercase hexadecimal SHA-1"
        ),
    )

    bench_res = if bench === nothing
        SMV.evaluate_spherical_momentum_benchmark()
    else
        bench
    end

    if isfile(record_path) && !accept_physics_change
        ok_bands, band_errs = SMV.check_physics_regression_bands(
            bench_res.metrics, record_path
        )
        if !ok_bands
            error(
                "Refusing to overwrite $record_path without --accept-physics-change: " *
                join(band_errs, "; "),
            )
        end
    end

    cfg = SMV.make_config(SMV.NX_REFERENCE; spherical_flow=true)
    cfg_sha = SMV.canonical_config_sha256(cfg)
    specs = metric_specs(bench_res.metrics, SMV.expected_metric_references())
    ts_utc = Dates.format(now(Dates.UTC), "yyyy-mm-ddTHH:MM:SSZ")

    io = IOBuffer()
    println(
        io,
        "# Validation record for steady Stokes flow driven by a body force in a sphere, against exact solutions of the spherical momentum equation.",
    )
    println(io, "# Generated by tools/export_validation_spherical_momentum.jl.")
    println(
        io,
        "# metadata.src_commit records the latest commit that modified src/ (git log -n 1 --format=%H -- src).",
    )
    println(io)
    println(io, "[metadata]")
    println(io, "id = \"spherical_momentum\"")
    println(io, "title = \"Spherical Momentum Validation\"")
    println(io, "page = \"validation/spherical_momentum.md\"")
    println(io, "src_commit = \"$(commit_str)\"")
    println(io, "timestamp_utc = \"$(ts_utc)\"")
    println(io, "julia_version = \"$(VERSION)\"")
    println(io, "config_sha256 = \"$(cfg_sha)\"")
    println(io, "reference_type = \"analytical\"")
    println(
        io,
        "reference_citation = \"Closed forms of azimuthal shear u = r Omega(r) e_theta for the momentum equation with and without the spherical metric terms, a manufactured flow r u = curl(r^3 X Y (1 - rho^2)^6 e_z) for the documented spherical-slice momentum equation, and a radial compaction flow u = U0 R grad (1 - rho^2)^5 with q = -u for the 3-D sphere equations of a compressible matrix\"",
    )
    println(io)
    cfg_sum = SMV.expected_config_summary(cfg)
    println(io, "[config]")
    println(io, "Nx = $(cfg_sum["Nx"])")
    for key in (
        "R_planet_m",
        "L_box_m",
        "eta_Pa_s",
        "eta_zeta_ratio",
        "phi",
        "k_m2",
        "mu_f_Pa_s",
        "rho_kg_m3",
        "rho_f_kg_m3",
        "shear_modulus_Pa",
        "F0_N_m3",
        "U0_m_s",
        "dt_yr",
        "P_amb_Pa",
        "fit_rho_min",
        "fit_rho_max",
        "eta_zeta_ratio_c",
        "P_amb_c_Pa",
    )
        println(io, "$(key) = $(repr(Float64(cfg_sum[key])))")
    end
    println(io, "b_power = $(cfg_sum["b_power"])")
    println(io, "c_power = $(cfg_sum["c_power"])")
    println(io, "fit_degree = $(cfg_sum["fit_degree"])")
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
    ok_valid, valid_errs = SMV.validate_record_dict(
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
        SMV.update_record_doc_block!(candidate_dict) ||
            @warn "Documentation page has no record table to update"
    end
    return bench_res
end

"""
    run_export_cli(args::Vector{String}=ARGS; src_commit=nothing, bench=nothing, check_clean_tree=true, dirty_status=nothing, repo_root=normpath(joinpath(@__DIR__, ".."))) -> Int

Run the command-line export tool for the spherical momentum validation record. Requires
`SphericalMomentumValidation` (`benchmarks/validation_spherical_momentum.jl`), see
`momentum_validation_module`. An argument other than `--accept-physics-change` that
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
                "Usage: julia --project=. tools/export_validation_spherical_momentum.jl [--accept-physics-change] [record_path]",
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
        write_spherical_momentum_record(
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
    include(joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_momentum.jl"))
    exit(Base.invokelatest(SphericalMomentumExport.export_main, ARGS))
end

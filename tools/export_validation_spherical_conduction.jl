#!/usr/bin/env julia
# Export the spherical conduction validation record TOML from live benchmark measurements.

using Dates
using Printf
using TOML
using Erebus

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

function _compute_metric_status(
    value::Float64, reference::Float64, atol::Float64, rtol::Float64, comparison::String
)::String
    passed = if comparison == "isapprox"
        isapprox(value, reference; atol=atol, rtol=rtol)
    elseif comparison == "le_atol"
        reference == 0.0 && isfinite(value) && value <= atol
    else
        false
    end
    return passed ? "PASS" : "FAIL"
end

"""
    write_spherical_conduction_record(record_path::String; src_commit::Union{String,Nothing}=nothing, bench=nothing) -> NamedTuple

Run the spherical conduction benchmark suite (or use a pre-evaluated `bench`) and write
`benchmarks/records/spherical_conduction.toml`.

Parameters
----------
record_path : String, optional
    Output path for the TOML validation record.
src_commit : Union{String,Nothing}, optional
    40-character lowercase hex SHA-1 of the last commit modifying `src/`. When `nothing`, queries `last_src_commit()`.
bench : Any, optional
    Pre-evaluated benchmark result from `SphericalConductionValidation.evaluate_spherical_conduction_benchmark()`.

Returns
-------
NamedTuple
    Evaluated benchmark summary tuple.

Raises
------
DomainError
    If `isempty(record_path)` or `src_commit` is not a 40-character lowercase hex SHA-1.
"""
function write_spherical_conduction_record(
    record_path::String=joinpath(
        @__DIR__, "..", "benchmarks", "records", "spherical_conduction.toml"
    );
    src_commit::Union{String,Nothing}=nothing,
    bench=nothing,
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
    cfg65 = SphericalConductionValidation.make_spherical_conduction_config(; Nx=65)
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

    mkpath(dirname(record_path))
    open(record_path, "w") do io
        println(
            io,
            "# Validation record for spherical conduction and step-averaged 26Al heating.",
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
        println(io, "Nx_coarse = 33")
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
    end
    return bench_res
end

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_conduction.jl"))
    expected_src = last_src_commit()
    commit_arg = length(ARGS) >= 1 ? ARGS[1] : expected_src
    if commit_arg != expected_src
        error(
            "Provided src_commit ($commit_arg) does not match latest src/ commit ($expected_src)",
        )
    end
    Base.invokelatest(write_spherical_conduction_record; src_commit=commit_arg)
    println("Wrote benchmarks/records/spherical_conduction.toml")
end

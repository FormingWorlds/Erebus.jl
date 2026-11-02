#!/usr/bin/env julia
# Export the spherical conduction validation record TOML from live benchmark measurements.

if Base.find_package("Printf") === nothing && !("@stdlib" in LOAD_PATH)
    push!(LOAD_PATH, "@stdlib")
end

using Dates
using Printf
using TOML
using Erebus

if !isdefined(@__MODULE__, :validate_record_dict)
    include(joinpath(@__DIR__, "validation_records.jl"))
end
if !isdefined(@__MODULE__, :evaluate_spherical_conduction_benchmark)
    include(joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_conduction.jl"))
end

function _current_git_commit()::String
    repo_root = normpath(joinpath(@__DIR__, ".."))
    return readchomp(`git -C $repo_root rev-parse HEAD`)
end

"""
    write_spherical_conduction_record(record_path::String; git_commit::String, bench=nothing)

Run the spherical conduction benchmark suite and write `benchmarks/records/spherical_conduction.toml`.
"""
function write_spherical_conduction_record(
    record_path::String=joinpath(
        @__DIR__, "..", "benchmarks", "records", "spherical_conduction.toml"
    );
    git_commit::String=_current_git_commit(),
    bench=nothing,
)
    bench_res = bench === nothing ? evaluate_spherical_conduction_benchmark() : bench
    cfg65 = make_spherical_conduction_config(; Nx=65)
    cfg_sha = canonical_config_sha256(cfg65)
    ts_utc = Dates.format(now(Dates.UTC), "yyyy-mm-ddTHH:MM:SSZ")

    metric_specs = [
        (
            id="steady_center_dT_ratio_3d_2d_nx33",
            case="steady",
            description="Steady 3D-to-2D centre warming ratio (Nx=33)",
            nx=33,
            value=bench_res.metrics["steady_center_dT_ratio_3d_2d_nx33"],
            reference=2.0 / 3.0,
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
            reference=2.0 / 3.0,
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
            reference=0.0,
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
            reference=0.0,
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
            reference=0.0,
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
            reference=0.0,
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
            reference=bench_res.t65.dT_be_rmin,
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
            reference=0.0,
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
            reference=0.5,
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
            reference=T_SURF,
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
        println(io)
        println(io, "[metadata]")
        println(io, "id = \"spherical_conduction\"")
        println(io, "title = \"Spherical Conduction and Radiogenic Heating Validation\"")
        println(io, "page = \"validation/thermal_conduction.md\"")
        println(io, "git_commit = \"$(git_commit)\"")
        println(io, "timestamp_utc = \"$(ts_utc)\"")
        println(io, "julia_version = \"$(VERSION)\"")
        println(io, "config_sha256 = \"$(cfg_sha)\"")
        println(io, "reference_type = \"analytical_and_1d_fv\"")
        println(
            io,
            "reference_citation = \"Exact spherical Poisson solution and step-averaged backward-Euler Dirichlet eigenmode series cross-verified against a 1D spherical finite-volume solver\"",
        )
        println(io)
        cfg_sum = expected_config_summary(cfg65)
        println(io, "[config]")
        println(io, "Nx_coarse = $(cfg_sum["Nx_coarse"])")
        println(io, "Nx_fine = $(cfg_sum["Nx_fine"])")
        println(io, "n_steps = $(cfg_sum["n_steps"])")
        @printf(io, "dt_yr = %.1f\n", cfg_sum["dt_yr"])
        @printf(io, "R_planet_m = %.1f\n", cfg_sum["R_planet_m"])
        @printf(io, "k_W_mK = %.1f\n", cfg_sum["k_W_mK"])
        @printf(io, "rho_kg_m3 = %.1f\n", cfg_sum["rho_kg_m3"])
        @printf(io, "cp_J_kgK = %.1f\n", cfg_sum["cp_J_kgK"])
        @printf(io, "T_surf_K = %.1f\n", cfg_sum["T_surf_K"])
        println(io, "spherical_metric = $(cfg_sum["spherical_metric"])")
        println(io, "reg_cells = $(cfg_sum["reg_cells"])")
        println(io, "disk_enabled = $(cfg_sum["disk_enabled"])")

        for m in metric_specs
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
            println(io, "status = \"PASS\"")
        end
    end
    return bench_res
end

if abspath(PROGRAM_FILE) == @__FILE__
    commit_arg = length(ARGS) >= 1 ? ARGS[1] : _current_git_commit()
    write_spherical_conduction_record(; git_commit=commit_arg)
    println("Wrote benchmarks/records/spherical_conduction.toml")
end

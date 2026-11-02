#!/usr/bin/env julia
# Write the data of the spherical conduction validation figures as CSV files into a folder.
# Usage: julia --project=. scripts/export_validation_conduction_data.jl <output_dir>

using Erebus

include(joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_conduction.jl"))
const SV = SphericalConductionValidation

"""
    write_columns(path, header, columns)

Write equal-length numeric `columns` to the CSV file `path` with the given `header` names.
"""
function write_columns(
    path::String, header::Vector{String}, columns::Vector{<:AbstractVector}
)
    open(path, "w") do io
        println(io, join(header, ","))
        for i in eachindex(columns[1])
            println(io, join((repr(Float64(c[i])) for c in columns), ","))
        end
    end
    return path
end

"""
    export_all(outdir)

Run the steady and transient cases, the marker seed runs and the heating change, and write the data.
"""
function export_all(outdir::String)
    mkpath(outdir)
    node_header = ["r_m", "dT_sim_K", "dT_ref_K"]
    for nx in (33, 65)
        s = SV.run_steady_spherical_conduction(nx)
        write_columns(
            joinpath(outdir, "steady_$(nx).csv"),
            node_header,
            [s.r_nodes, s.dT_sim_nodes, s.dT_ref_nodes],
        )
        t = SV.run_transient_spherical_conduction(nx)
        write_columns(
            joinpath(outdir, "transient_$(nx).csv"),
            node_header,
            [t.r_nodes, t.dT_sim_nodes, t.dT_ref_nodes],
        )
    end
    radii = collect(range(0.0, SV.R_PLANET; length=201))
    write_columns(
        joinpath(outdir, "reference_curves.csv"),
        ["r_m", "dT_steady_K", "dT_be_K", "dT_cont_K"],
        [
            radii,
            [SV.H_STEADY * (SV.R_PLANET^2 - r^2) / (6.0 * SV.K_COND) for r in radii],
            [SV.analytical_T_BE(r) - SV.T_SURF for r in radii],
            [SV.analytical_T_continuous(r) - SV.T_SURF for r in radii],
        ],
    )
    r_fv, t_fv = SV.solve_1d_spherical_fv()
    write_columns(joinpath(outdir, "fv_1d.csv"), ["r_m", "dT_K"], [r_fv, t_fv .- SV.T_SURF])
    seeds = [42, 1, 2, 3, 4, 101, 202, 303, 404]
    runs = [SV.run_transient_spherical_conduction(65; seed=s) for s in seeds]
    write_columns(
        joinpath(outdir, "seeds.csv"),
        ["seed", "dT_center_K", "l2_rel"],
        [seeds, [c.dT_center for c in runs], [c.l2_rel for c in runs]],
    )
    h = SV.run_transient_spherical_conduction(65; f_al=1.02 * SV.F_AL)
    write_columns(
        joinpath(outdir, "heating_plus_2pct.csv"),
        ["dT_center_K", "l2_rel"],
        [[h.dT_center], [h.l2_rel]],
    )
    return outdir
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("Usage: export_validation_conduction_data.jl <output_dir>")
    export_all(ARGS[1])
end

#!/usr/bin/env julia
# Write the data of the spherical compaction validation figure as CSV files into a folder.
# Usage: julia --project=. scripts/export_validation_compaction_data.jl <output_dir>

using Erebus
using TOML

include(joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_compaction.jl"))
const SCV = SphericalCompactionValidation
const DEFAULT_RECORD = joinpath(
    @__DIR__, "..", "benchmarks", "records", "spherical_compaction.toml"
)

"""
    write_columns(path, header, columns)

Write equal-length numeric `columns` to the CSV file `path` with the given `header` names.

Parameters
----------
path : String
    Output file path.
header : Vector{String}
    One name per column.
columns : Vector{<:AbstractVector}
    Numeric columns of equal length, one per name in `header`.

Returns
-------
String
    `path`.

Raises
------
ArgumentError
    If `columns` is empty, `header` and `columns` differ in length, or the columns differ in
    length. The file is then left as it was.
"""
function write_columns(
    path::String, header::Vector{String}, columns::Vector{<:AbstractVector}
)
    isempty(columns) && throw(ArgumentError("at least one column is needed"))
    length(header) == length(columns) || throw(
        ArgumentError("$(length(header)) names for $(length(columns)) columns in $path")
    )
    n = length(columns[1])
    all(c -> length(c) == n, columns) ||
        throw(ArgumentError("the columns of $path differ in length"))
    open(path, "w") do io
        println(io, join(header, ","))
        for i in eachindex(columns[1])
            println(io, join((repr(Float64(c[i])) for c in columns), ","))
        end
    end
    return path
end

"""
    node_profile(N; spherical_flow) -> NamedTuple

Solve the case with `SCV.run_case` and return the radius [m] and the computed and exact
effective pressure `Pe = Pt - Pf` [Pa] of the P nodes with `r <= R_PLANET` that the error
measures read (`SCV.run_case` compares the nodes with `r <= 0.9 R_PLANET`).

Parameters
----------
N : Int
    Basic grid resolution in x and y, at least `SCV.NX_MIN`.
spherical_flow : Bool
    Value of `geometry.spherical_flow`.

Returns
-------
NamedTuple
    `r`, `pe` and `pe_exact`: vectors over the P nodes inside the planet.

Raises
------
ArgumentError
    If `N < SCV.NX_MIN`.
"""
function node_profile(N::Int; spherical_flow::Bool)
    run = SCV.run_case(N; spherical_flow=spherical_flow)
    d = spherical_flow ? 3 : 2
    r, pe, pe_exact = Float64[], Float64[], Float64[]
    for j in 3:(length(run.xp) - 2), i in 3:(length(run.yp) - 2)
        ri = hypot(run.xp[j] - SCV.X_CENTRE, run.yp[i] - SCV.Y_CENTRE)
        ri <= SCV.R_PLANET || continue
        push!(r, ri)
        push!(pe, run.pe[i, j])
        push!(pe_exact, SCV.exact_pe(ri, d))
    end
    return (; r, pe, pe_exact)
end

"""
    export_all(outdir; record_path=DEFAULT_RECORD)

Run the node profiles at `Nx = SCV.NX_REFERENCE`, write the exact curves, and copy the
profile and centre errors of the record at every grid in `SCV.NX_LIST`.

Parameters
----------
outdir : String
    Output folder, created if missing.
record_path : String, optional
    Validation record TOML file with the ids `profile_l2_inner_{off,on}_nx<N>` and
    `centre_err_{off,on}_nx<N>`.

Returns
-------
String
    `outdir`, holding `profile_off_65.csv`, `profile_on_65.csv`, `reference_curves.csv` and
    `convergence.csv`.

Raises
------
KeyError
    If the record lacks an error id for a grid in `SCV.NX_LIST`.
"""
function export_all(outdir::String; record_path::String=DEFAULT_RECORD)
    mkpath(outdir)
    for (sf, tag) in ((false, "off"), (true, "on"))
        p = node_profile(SCV.NX_REFERENCE; spherical_flow=sf)
        write_columns(
            joinpath(outdir, "profile_$(tag)_$(SCV.NX_REFERENCE).csv"),
            ["r_m", "pe_Pa", "pe_exact_Pa"],
            [p.r, p.pe, p.pe_exact],
        )
    end
    radii = collect(range(0.0, SCV.R_PLANET; length=201))
    write_columns(
        joinpath(outdir, "reference_curves.csv"),
        ["r_m", "pe_exact_off_Pa", "pe_exact_on_Pa"],
        [radii, [SCV.exact_pe(r, 2) for r in radii], [SCV.exact_pe(r, 3) for r in radii]],
    )
    values = Dict(
        m["id"] => Float64(m["value"]) for m in TOML.parsefile(record_path)["metrics"]
    )
    nx = collect(SCV.NX_LIST)
    write_columns(
        joinpath(outdir, "convergence.csv"),
        ["nx", "l2_off", "l2_on", "centre_off", "centre_on"],
        [
            nx,
            [values["profile_l2_inner_off_nx$(n)"] for n in nx],
            [values["profile_l2_inner_on_nx$(n)"] for n in nx],
            [values["centre_err_off_nx$(n)"] for n in nx],
            [values["centre_err_on_nx$(n)"] for n in nx],
        ],
    )
    return outdir
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("Usage: export_validation_compaction_data.jl <output_dir>")
    export_all(ARGS[1])
end

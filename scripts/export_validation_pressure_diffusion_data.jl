#!/usr/bin/env julia
# Write the data of the spherical pressure diffusion validation figure as CSV files into a folder.
# Usage: julia --project=. scripts/export_validation_pressure_diffusion_data.jl <output_dir>

using Erebus
using TOML

include(
    joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_pressure_diffusion.jl")
)
const SPDV = SphericalPressureDiffusionValidation
const DEFAULT_RECORD = joinpath(
    @__DIR__, "..", "benchmarks", "records", "spherical_pressure_diffusion.toml"
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
    cell_profile(run) -> NamedTuple

Return the radius `x = r / R_PLANET` and the nondimensional fluid overpressure
`(Pf - P_AMB) / P0` after 5 and 20 steps of the interior rock P nodes of a `SPDV.run_case`
result, the nodes that the error measures read.

Parameters
----------
run : NamedTuple
    Result of `SPDV.run_case` with `xp`, `yp`, `pf_m5` and `pf`.

Returns
-------
NamedTuple
    `x`, `p_m5` and `p_m20`: vectors over the rock P nodes.
"""
function cell_profile(run)
    x, p5, p20 = Float64[], Float64[], Float64[]
    for j in 2:(length(run.xp) - 1), i in 2:(length(run.yp) - 1)
        r = hypot(run.xp[j] - SPDV.X_CENTRE, run.yp[i] - SPDV.Y_CENTRE)
        r <= SPDV.R_PLANET || continue
        push!(x, r / SPDV.R_PLANET)
        push!(p5, (run.pf_m5[i, j] - SPDV.P_AMB) / SPDV.P0)
        push!(p20, (run.pf[i, j] - SPDV.P_AMB) / SPDV.P0)
    end
    return (; x, p_m5=p5, p_m20=p20)
end

"""
    export_all(outdir; record_path=DEFAULT_RECORD, runs=nothing, nx=SPDV.NX_REFERENCE)

Write the cell profiles at `Nx = nx` after 5 and 20 steps, the backward-Euler series
curves, and the absolute centre and profile errors at 20 steps of the record at every grid
in `SPDV.NX_LIST`.

Parameters
----------
outdir : String
    Output folder, created if missing.
record_path : String, optional
    Validation record TOML file with the ids `l2_m20_{off,on}_nx<N>` and
    `centre_err_m20_{off,on}_nx<N>`.
runs : Union{AbstractDict,Nothing}, optional
    `run_case` results keyed by `(N, spherical_flow, :base)` with the runs at `nx`, or
    `nothing` to run them here (default).
nx : Int, optional
    Grid of the cell profiles (default `SPDV.NX_REFERENCE`, the grid of the figure).

Returns
-------
String
    `outdir`, holding `profile_off_<nx>.csv`, `profile_on_<nx>.csv`,
    `reference_curves.csv` and `convergence.csv`.

Raises
------
KeyError
    If the record lacks an error id for a grid in `SPDV.NX_LIST`, or `runs` lacks a run
    at `nx`.
ArgumentError
    If `runs` is `nothing` and `SPDV.run_case` rejects `nx`.
"""
function export_all(
    outdir::String;
    record_path::String=DEFAULT_RECORD,
    runs=nothing,
    nx::Int=SPDV.NX_REFERENCE,
)
    mkpath(outdir)
    for (sf, tag) in ((false, "off"), (true, "on"))
        run = if runs === nothing
            SPDV.run_case(nx; spherical_flow=sf)
        else
            runs[(nx, sf, :base)]
        end
        p = cell_profile(run)
        write_columns(
            joinpath(outdir, "profile_$(tag)_$(nx).csv"),
            ["x", "p_m5", "p_m20"],
            [p.x, p.p_m5, p.p_m20],
        )
    end
    xs = collect(range(0.0, 1.0; length=201))
    write_columns(
        joinpath(outdir, "reference_curves.csv"),
        ["x", "p_be_m5_off", "p_be_m5_on", "p_be_m20_off", "p_be_m20_on"],
        [
            xs,
            SPDV.exact_p_cells(xs, 2, 5),
            SPDV.exact_p_cells(xs, 3, 5),
            SPDV.exact_p_cells(xs, 2, 20),
            SPDV.exact_p_cells(xs, 3, 20),
        ],
    )
    values = Dict(
        m["id"] => Float64(m["value"]) for m in TOML.parsefile(record_path)["metrics"]
    )
    grids = collect(SPDV.NX_LIST)
    write_columns(
        joinpath(outdir, "convergence.csv"),
        ["nx", "l2_off", "l2_on", "centre_off", "centre_on"],
        [
            grids,
            [values["l2_m20_off_nx$(n)"] for n in grids],
            [values["l2_m20_on_nx$(n)"] for n in grids],
            [values["centre_err_m20_off_nx$(n)"] for n in grids],
            [values["centre_err_m20_on_nx$(n)"] for n in grids],
        ],
    )
    return outdir
end

"""
    export_data_main(args=ARGS; kwargs...)

Run the script: write the figure data into the one output folder that `args` names.

Parameters
----------
args : Vector{String}, optional
    Command-line arguments, exactly one output folder that is not empty or blank and does
    not start with `-` (default `ARGS`).
kwargs
    Keyword arguments passed on to `export_all`.

Returns
-------
String
    The output folder.

Raises
------
ArgumentError
    If `args` does not hold exactly one entry, or the entry is empty, blank or starts with
    `-`, before any file is written.
"""
function export_data_main(args::Vector{String}=ARGS; kwargs...)
    usage = "Usage: export_validation_pressure_diffusion_data.jl <output_dir>"
    length(args) == 1 || throw(ArgumentError(usage))
    isempty(strip(args[1])) &&
        throw(ArgumentError("the output folder must be non-empty. $usage"))
    startswith(args[1], "-") && throw(ArgumentError("unknown option '$(args[1])'. $usage"))
    return export_all(args[1]; kwargs...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    export_data_main(ARGS)
end

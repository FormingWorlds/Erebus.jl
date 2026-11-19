#!/usr/bin/env julia
# Write the data of the spherical momentum validation figure as CSV files into a folder.
# Usage: julia --project=. scripts/export_validation_momentum_data.jl <output_dir>

using Erebus
using TOML

include(joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_momentum.jl"))
const SMV = SphericalMomentumValidation
const DEFAULT_RECORD = joinpath(
    @__DIR__, "..", "benchmarks", "records", "spherical_momentum.toml"
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
    axis_profile(run) -> NamedTuple

Return the azimuthal velocity of a case `:A` run on the 2 axes through the centre with
`0 < r < R`: `u_theta = sign(x - xc) vy` on the vy nodes of the row `y = yc` and
`-sign(y - yc) vx` on the vx nodes of the column `x = xc`, where the flow is purely
azimuthal, in units of `SMV.velocity_scale(:A)`.

Parameters
----------
run : NamedTuple
    Result of `SMV.run_case` with `vx`, `vy` and `grid`; the centre lies on a basic node.

Returns
-------
NamedTuple
    `rho`: `r / R`; `u`: `u_theta`, dimensionless; sorted by `rho`.
"""
function axis_profile(run)
    g = run.grid
    v = SMV.velocity_scale(:A)
    bx = argmin(abs.(g.xvx .- SMV.X_CENTRE))
    ay = argmin(abs.(g.yvy .- SMV.Y_CENTRE))
    rho, u = Float64[], Float64[]
    for j in 2:(g.Nx1 - 1)
        x = (g.xvy[j] - SMV.X_CENTRE) / SMV.R_PLANET
        0 < abs(x) < 1 || continue
        push!(rho, abs(x))
        push!(u, sign(x) * run.vy[ay, j] / v)
    end
    for i in 2:(g.Ny1 - 1)
        y = (g.yvx[i] - SMV.Y_CENTRE) / SMV.R_PLANET
        0 < abs(y) < 1 || continue
        push!(rho, abs(y))
        push!(u, -sign(y) * run.vx[i, bx] / v)
    end
    order = sortperm(rho)
    return (; rho=rho[order], u=u[order])
end

"""
    export_all(outdir; record_path=DEFAULT_RECORD, runs=nothing, nx=SMV.NX_REFERENCE)

Write the case `:A` axis profiles at `Nx = nx` with the option on and off, the exact
profiles `u_theta = rho Omega(rho)`, and the errors of the record at every grid in
`SMV.NX_LIST`.

Parameters
----------
outdir : String
    Output folder, created if missing.
record_path : String, optional
    Validation record TOML file with the ids `l2_a_on_nx<N>`, `l2_a_off_nx<N>`,
    `l2_b_nx<N>`, `s33_b_nx<N>`, `centre_fit_a_on_nx<N>` and `centre_node_a_on_nx<N>`.
runs : Union{AbstractDict,Nothing}, optional
    `run_case` results keyed by `(N, :A, spherical_flow, :base)` with the runs at `nx`, or
    `nothing` to run them here (default).
nx : Int, optional
    Grid of the profiles (default `SMV.NX_REFERENCE`, the grid of the figure).

Returns
-------
String
    `outdir`, holding `profile_on_<nx>.csv`, `profile_off_<nx>.csv`,
    `reference_curves.csv` and `convergence.csv`.

Raises
------
KeyError
    If the record lacks an error id for a grid in `SMV.NX_LIST`, or `runs` lacks a run at
    `nx`.
ArgumentError
    If `runs` is `nothing` and `SMV.run_case` rejects `nx`.
"""
function export_all(
    outdir::String;
    record_path::String=DEFAULT_RECORD,
    runs=nothing,
    nx::Int=SMV.NX_REFERENCE,
)
    mkpath(outdir)
    for (sf, tag) in ((true, "on"), (false, "off"))
        run = if runs === nothing
            SMV.run_case(nx; case=:A, spherical_flow=sf)
        else
            runs[(nx, :A, sf, :base)]
        end
        p = axis_profile(run)
        write_columns(
            joinpath(outdir, "profile_$(tag)_$(nx).csv"), ["rho", "u"], [p.rho, p.u]
        )
    end
    rhos = collect(range(0.0, 1.0; length=201))
    write_columns(
        joinpath(outdir, "reference_curves.csv"),
        ["rho", "u_on", "u_off"],
        [
            rhos,
            [r * SMV.omega_hat(:A, r; spherical_flow=true) for r in rhos],
            [r * SMV.omega_hat(:A, r; spherical_flow=false) for r in rhos],
        ],
    )
    values = Dict(
        m["id"] => Float64(m["value"]) for m in TOML.parsefile(record_path)["metrics"]
    )
    grids = collect(SMV.NX_LIST)
    names = ("l2_a_on", "l2_a_off", "l2_b", "s33_b", "centre_fit_a_on", "centre_node_a_on")
    write_columns(
        joinpath(outdir, "convergence.csv"),
        ["nx", names...],
        [grids, ([values["$(name)_nx$(n)"] for n in grids] for name in names)...],
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
    usage = "Usage: export_validation_momentum_data.jl <output_dir>"
    length(args) == 1 || throw(ArgumentError(usage))
    isempty(strip(args[1])) &&
        throw(ArgumentError("the output folder must be non-empty. $usage"))
    startswith(args[1], "-") && throw(ArgumentError("unknown option '$(args[1])'. $usage"))
    return export_all(args[1]; kwargs...)
end

if abspath(PROGRAM_FILE) == @__FILE__
    export_data_main(ARGS)
end

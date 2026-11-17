#!/usr/bin/env julia
# Write the data of the spherical Darcy flow validation figure as CSV files into a folder.
# Usage: julia --project=. scripts/export_validation_darcy_flow_data.jl <output_dir>

using Erebus
using Logging

include(joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_darcy_flow.jl"))
const SDF = SphericalDarcyFlowValidation

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
BoundsError
    If a column is shorter than the first one.
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
    exact_excess_pressure(r, spherical_flow) -> Float64

Return the exact excess fluid pressure `p - P_AMB` [Pa] at radius `r` [m]: `d = 3` for
`spherical_flow = true` and `d = 2` otherwise.

Parameters
----------
r : Real
    Radius [m].
spherical_flow : Bool
    Value of `geometry.spherical_flow`.

Returns
-------
Float64
    Exact excess pressure [Pa].
"""
function exact_excess_pressure(r::Real, spherical_flow::Bool)::Float64
    d = spherical_flow ? 3.0 : 2.0
    return SDF.MU_F * SDF.S_SOURCE * (SDF.R_PLANET^2 - r^2) / (2d * SDF.K_PERM)
end

"""
    node_profile(N; spherical_flow) -> NamedTuple

Solve the case as `SDF.run_case` does and return the radius [m] and the computed and exact
excess pressure `p - P_AMB` [Pa] of the P nodes inside the planet that `SDF.run_case` compares.

Parameters
----------
N : Int
    Basic grid resolution in x and y, at least `SDF.NX_MIN`.
spherical_flow : Bool
    Value of `geometry.spherical_flow`.

Returns
-------
NamedTuple
    `r`, `dp` and `dp_exact`: vectors over the P nodes inside the planet.

Raises
------
ArgumentError
    If `N < SDF.NX_MIN`.
"""
function node_profile(N::Int; spherical_flow::Bool)
    cfg = SDF.make_config(N; spherical_flow=spherical_flow)
    state, coords, ws, tio, _ = with_logger(NullLogger()) do
        return Erebus.init_simulation(cfg; output_path=mktempdir())
    end
    tio === nothing || close(tio)
    state.timestep = 1
    state.accumulators.P_amb = SDF.P_AMB
    SDF.set_fields!(state, coords)
    with_logger(NullLogger()) do
        return Erebus.solve_stokes_darcy!(
            state, coords, cfg, ws; titer=1, dt_step_initial=state.dt
        )
    end
    r, dp, dp_exact = Float64[], Float64[], Float64[]
    for j in 3:(coords.Nx1 - 2), i in 3:(coords.Ny1 - 2)
        ri = hypot(coords.xp[j] - SDF.X_CENTRE, coords.yp[i] - SDF.Y_CENTRE)
        ri <= SDF.R_PLANET || continue
        push!(r, ri)
        push!(dp, state.grids.pf[i, j] - SDF.P_AMB)
        push!(dp_exact, exact_excess_pressure(ri, spherical_flow))
    end
    return (; r, dp, dp_exact)
end

"""
    export_all(outdir)

Run the node profiles at `Nx = 65`, the profile error at every grid in `SDF.NX_LIST`, and the
exact curves, and write the data.

Parameters
----------
outdir : String
    Output folder, created if missing.

Returns
-------
String
    `outdir`, holding `profile_off_65.csv`, `profile_on_65.csv`, `reference_curves.csv` and
    `convergence.csv`.
"""
function export_all(outdir::String)
    mkpath(outdir)
    for (sf, tag) in ((false, "off"), (true, "on"))
        p = node_profile(SDF.NX_REFERENCE; spherical_flow=sf)
        write_columns(
            joinpath(outdir, "profile_$(tag)_$(SDF.NX_REFERENCE).csv"),
            ["r_m", "dp_Pa", "dp_exact_Pa"],
            [p.r, p.dp, p.dp_exact],
        )
    end
    radii = collect(range(0.0, SDF.R_PLANET; length=201))
    write_columns(
        joinpath(outdir, "reference_curves.csv"),
        ["r_m", "dp_exact_off_Pa", "dp_exact_on_Pa"],
        [
            radii,
            [exact_excess_pressure(r, false) for r in radii],
            [exact_excess_pressure(r, true) for r in radii],
        ],
    )
    nx = collect(SDF.NX_LIST)
    err = Dict(
        sf => [SDF.run_case(n; spherical_flow=sf).l2_inner for n in nx] for
        sf in (false, true)
    )
    write_columns(
        joinpath(outdir, "convergence.csv"),
        ["nx", "l2_off", "l2_on"],
        [nx, err[false], err[true]],
    )
    return outdir
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("Usage: export_validation_darcy_flow_data.jl <output_dir>")
    export_all(ARGS[1])
end

#!/usr/bin/env julia
# Write the data of the spherical Darcy flow validation figure as CSV files into a folder.
# Usage: julia --project=. scripts/export_validation_darcy_flow_data.jl <output_dir>

using Erebus
using Logging
using TOML

include(joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_darcy_flow.jl"))
const SDF = SphericalDarcyFlowValidation
const DEFAULT_RECORD = joinpath(
    @__DIR__, "..", "benchmarks", "records", "spherical_darcy_flow.toml"
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
    node_profile(N; spherical_flow, jump=false) -> NamedTuple

Solve the case as `SDF.run_case` does (with the core permeability `SDF.K_CORE` for `r <
SDF.A_JUMP` when `jump` is `true`) and return the radius [m] and the computed and exact
excess pressure `p - P_AMB` [Pa] of all P nodes with `r <= R_PLANET` (`SDF.run_case` compares
the nodes with `r <= 0.9 R_PLANET`).

Parameters
----------
N : Int
    Basic grid resolution in x and y, at least `SDF.NX_MIN`.
spherical_flow : Bool
    Value of `geometry.spherical_flow`.
jump : Bool, optional
    Whether the core has the permeability `SDF.K_CORE` (default `false`).

Returns
-------
NamedTuple
    `r`, `dp` and `dp_exact`: vectors over the P nodes inside the planet.

Raises
------
ArgumentError
    If `N < SDF.NX_MIN`.
"""
function node_profile(N::Int; spherical_flow::Bool, jump::Bool=false)
    cfg = SDF.make_config(N; spherical_flow=spherical_flow)
    state, coords = mktempdir() do out
        state, coords, ws, tio, _ = with_logger(NullLogger()) do
            return Erebus.init_simulation(cfg; output_path=out)
        end
        tio === nothing || close(tio)
        state.timestep = 1
        state.accumulators.P_amb = SDF.P_AMB
        SDF.set_fields!(state, coords; jump=jump)
        with_logger(NullLogger()) do
            return Erebus.solve_stokes_darcy!(
                state, coords, cfg, ws; titer=1, dt_step_initial=state.dt
            )
        end
        return state, coords
    end
    r, dp, dp_exact = Float64[], Float64[], Float64[]
    for j in 3:(coords.Nx1 - 2), i in 3:(coords.Ny1 - 2)
        ri = hypot(coords.xp[j] - SDF.X_CENTRE, coords.yp[i] - SDF.Y_CENTRE)
        ri <= SDF.R_PLANET || continue
        push!(r, ri)
        push!(dp, state.grids.pf[i, j] - SDF.P_AMB)
        push!(dp_exact, SDF.exact_excess_pressure(ri, spherical_flow ? 3 : 2; jump=jump))
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
    `outdir`, holding `profile_off_65.csv`, `profile_on_65.csv`, `reference_curves.csv`,
    `convergence.csv` and the files of `export_jump`.
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
            [SDF.exact_excess_pressure(r, 2) for r in radii],
            [SDF.exact_excess_pressure(r, 3) for r in radii],
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
    export_jump(outdir)
    return outdir
end

"""
    export_jump(outdir; record_path=DEFAULT_RECORD)

Write the permeability-jump data: the node profiles at `Nx = SDF.NX_REFERENCE`, the exact
curves with the jump, and the profile errors of the record at every grid in `SDF.NX_LIST`.

Parameters
----------
outdir : String
    Output folder, created if missing.
record_path : String, optional
    Validation record TOML file with the ids `profile_l2_inner_jump_{off,on}_nx<N>`.

Returns
-------
String
    `outdir`, holding `jump_profile_off_65.csv`, `jump_profile_on_65.csv`,
    `jump_reference_curves.csv` and `jump_convergence.csv`.

Raises
------
KeyError
    If the record lacks a jump error id for a grid in `SDF.NX_LIST`.
"""
function export_jump(outdir::String; record_path::String=DEFAULT_RECORD)
    mkpath(outdir)
    for (sf, tag) in ((false, "off"), (true, "on"))
        p = node_profile(SDF.NX_REFERENCE; spherical_flow=sf, jump=true)
        write_columns(
            joinpath(outdir, "jump_profile_$(tag)_$(SDF.NX_REFERENCE).csv"),
            ["r_m", "dp_Pa", "dp_exact_Pa"],
            [p.r, p.dp, p.dp_exact],
        )
    end
    radii = collect(range(0.0, SDF.R_PLANET; length=201))
    write_columns(
        joinpath(outdir, "jump_reference_curves.csv"),
        ["r_m", "dp_exact_off_Pa", "dp_exact_on_Pa"],
        [
            radii,
            [SDF.exact_excess_pressure(r, 2; jump=true) for r in radii],
            [SDF.exact_excess_pressure(r, 3; jump=true) for r in radii],
        ],
    )
    values = Dict(
        m["id"] => Float64(m["value"]) for m in TOML.parsefile(record_path)["metrics"]
    )
    nx = collect(SDF.NX_LIST)
    write_columns(
        joinpath(outdir, "jump_convergence.csv"),
        ["nx", "l2_off", "l2_on"],
        [
            nx,
            [values["profile_l2_inner_jump_off_nx$(n)"] for n in nx],
            [values["profile_l2_inner_jump_on_nx$(n)"] for n in nx],
        ],
    )
    return outdir
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("Usage: export_validation_darcy_flow_data.jl <output_dir>")
    export_all(ARGS[1])
end

using Test
using Erebus
using StaticArrays
using Random

"""
Assert that `validate_config` rejects an invalid configuration with `ArgumentError`.
"""
macro reject_config(args...)
    if length(args) == 1
        arg = args[1]
        if Meta.isexpr(arg, :(=)) || Meta.isexpr(arg, :kw)
            cfg_expr = :(SimulationConfig(; $(args...)))
        else
            cfg_expr = arg
        end
    else
        cfg_expr = :(SimulationConfig(; $(args...)))
    end
    return quote
        @test_throws ArgumentError validate_config($(esc(cfg_expr)))
    end
end

if !isdefined(@__MODULE__, :G)
    include("../src/constants.jl")
end

"""
    test_grid_coordinates(; Nx=33, Ny=33, xsize=140_000.0, ysize=140_000.0)

Construct fresh grid coordinates for test fixtures.
"""
function test_grid_coordinates(; Nx=33, Ny=33, xsize=140_000.0, ysize=140_000.0)
    return GridCoordinates(Nx, Ny; xsize=xsize, ysize=ysize)
end

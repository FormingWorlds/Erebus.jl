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

if !isdefined(@__MODULE__, :rgen)
    const rgen = MersenneTwister(42)
end

if !isdefined(@__MODULE__, :G)
    include("../src/constants.jl")
end

# Standard test fixture grid coordinates
const coords_test_default = default_grid_coordinates()
const Nx = coords_test_default.Nx
const Ny = coords_test_default.Ny
const dx = coords_test_default.dx
const dy = coords_test_default.dy
const xsize = coords_test_default.xsize
const ysize = coords_test_default.ysize
const xcenter = coords_test_default.xcenter
const ycenter = coords_test_default.ycenter
const x = coords_test_default.x
const y = coords_test_default.y
const xp = coords_test_default.xp
const yp = coords_test_default.yp
const xvx = coords_test_default.xvx
const yvx = coords_test_default.yvx
const xvy = coords_test_default.xvy
const yvy = coords_test_default.yvy
const Nx1 = coords_test_default.Nx1
const Ny1 = coords_test_default.Ny1
const start_marknum = coords_test_default.start_marknum
const jmin_basic = coords_test_default.jmin_basic
const imin_basic = coords_test_default.imin_basic
const jmax_basic = coords_test_default.jmax_basic
const imax_basic = coords_test_default.imax_basic
const jmin_vx = coords_test_default.jmin_vx
const imin_vx = coords_test_default.imin_vx
const jmax_vx = coords_test_default.jmax_vx
const imax_vx = coords_test_default.imax_vx
const jmin_vy = coords_test_default.jmin_vy
const imin_vy = coords_test_default.imin_vy
const jmax_vy = coords_test_default.jmax_vy
const imax_vy = coords_test_default.imax_vy
const jmin_p = coords_test_default.jmin_p
const imin_p = coords_test_default.imin_p
const jmax_p = coords_test_default.jmax_p
const imax_p = coords_test_default.imax_p
const Nxmc = coords_test_default.Nxmc
const Nymc = coords_test_default.Nymc
const Nxm = coords_test_default.Nxm
const Nym = coords_test_default.Nym
const dxm = coords_test_default.dxm
const dym = coords_test_default.dym
const xxm = coords_test_default.xxm
const yym = coords_test_default.yym
const jmin_m = coords_test_default.jmin_m
const imin_m = coords_test_default.imin_m
const jmax_m = coords_test_default.jmax_m
const imax_m = coords_test_default.imax_m
const vxleft = 0.0
const vxright = 0.0
const vytop = 0.0
const vybottom = 0.0

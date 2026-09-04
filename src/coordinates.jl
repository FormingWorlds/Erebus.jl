"""
Coordinate and spatial domain geometry for dynamic grid resolutions in Erebus.jl.
"""

using DocStringExtensions

"""
Coordinate arrays, grid dimensions, and spacing parameters for staggered grid discretization.

$(FIELDS)
"""
Base.@kwdef struct GridCoordinates
    Nx::Int
    Ny::Int
    Nx1::Int
    Ny1::Int
    xsize::Float64
    ysize::Float64
    dx::Float64
    dy::Float64
    xcenter::Float64
    ycenter::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    xvx::Vector{Float64}
    yvx::Vector{Float64}
    xvy::Vector{Float64}
    yvy::Vector{Float64}
    xp::Vector{Float64}
    yp::Vector{Float64}
    jmin_basic::Int
    imin_basic::Int
    jmax_basic::Int
    imax_basic::Int
    jmin_vx::Int
    imin_vx::Int
    jmax_vx::Int
    imax_vx::Int
    jmin_vy::Int
    imin_vy::Int
    jmax_vy::Int
    imax_vy::Int
    jmin_p::Int
    imin_p::Int
    jmax_p::Int
    imax_p::Int
    Nxmc::Int
    Nymc::Int
    Nxm::Int
    Nym::Int
    dxm::Float64
    dym::Float64
    start_marknum::Int
    xxm::Vector{Float64}
    yym::Vector{Float64}
    jmin_m::Int
    imin_m::Int
    jmax_m::Int
    imax_m::Int
end

"""
    GridCoordinates(Nx::Int, Ny::Int; xsize=140_000.0, ysize=140_000.0, Nxmc=4, Nymc=4)

Construct a `GridCoordinates` object for arbitrary grid resolutions and domain sizes.
"""
function GridCoordinates(
    Nx::Int,
    Ny::Int;
    xsize::Float64=140_000.0,
    ysize::Float64=140_000.0,
    Nxmc::Int=4,
    Nymc::Int=4,
)
    Nx >= 3 || throw(ArgumentError("Grid Nx must be >= 3, got $Nx"))
    Ny >= 3 || throw(ArgumentError("Grid Ny must be >= 3, got $Ny"))
    xsize > 0.0 || throw(ArgumentError("Domain xsize must be > 0, got $xsize"))
    ysize > 0.0 || throw(ArgumentError("Domain ysize must be > 0, got $ysize"))
    Nxmc >= 1 || throw(ArgumentError("Nxmc must be >= 1, got $Nxmc"))
    Nymc >= 1 || throw(ArgumentError("Nymc must be >= 1, got $Nymc"))

    Nx1 = Nx + 1
    Ny1 = Ny + 1
    dx = xsize / (Nx - 1)
    dy = ysize / (Ny - 1)
    xcenter = xsize / 2.0
    ycenter = ysize / 2.0

    x = [j for j in 0.0:dx:xsize]
    y = [i for i in 0.0:dy:ysize]
    xvx = [j for j in 0.0:dx:(xsize + dx)]
    yvx = [i for i in (-dy / 2.0):dy:(ysize + dy / 2.0)]
    xvy = [j for j in (-dx / 2.0):dx:(xsize + dx / 2.0)]
    yvy = [i for i in 0.0:dy:(ysize + dy)]
    xp = [j for j in (-dx / 2.0):dx:(xsize + dx / 2.0)]
    yp = [i for i in (-dy / 2.0):dy:(ysize + dy / 2.0)]

    Nxm = (Nx - 1) * Nxmc
    Nym = (Ny - 1) * Nymc
    dxm = xsize / Nxm
    dym = ysize / Nym
    start_marknum = Nxm * Nym
    xxm = [j for j in (dxm / 2.0):dxm:(xsize - dxm / 2.0)]
    yym = [i for i in (dym / 2.0):dym:(ysize - dym / 2.0)]

    return GridCoordinates(;
        Nx=Nx,
        Ny=Ny,
        Nx1=Nx1,
        Ny1=Ny1,
        xsize=xsize,
        ysize=ysize,
        dx=dx,
        dy=dy,
        xcenter=xcenter,
        ycenter=ycenter,
        x=x,
        y=y,
        xvx=xvx,
        yvx=yvx,
        xvy=xvy,
        yvy=yvy,
        xp=xp,
        yp=yp,
        jmin_basic=1,
        imin_basic=1,
        jmax_basic=Nx - 1,
        imax_basic=Ny - 1,
        jmin_vx=1,
        imin_vx=1,
        jmax_vx=Nx - 1,
        imax_vx=Ny,
        jmin_vy=1,
        imin_vy=1,
        jmax_vy=Nx,
        imax_vy=Ny - 1,
        jmin_p=1,
        imin_p=1,
        jmax_p=Nx,
        imax_p=Ny,
        Nxmc=Nxmc,
        Nymc=Nymc,
        Nxm=Nxm,
        Nym=Nym,
        dxm=dxm,
        dym=dym,
        start_marknum=start_marknum,
        xxm=xxm,
        yym=yym,
        jmin_m=1,
        imin_m=1,
        jmax_m=Nxm - 1,
        imax_m=Nym - 1,
    )
end

"""
    GridCoordinates(cfg::GridConfig; Nxmc=4, Nymc=4)

Construct a `GridCoordinates` object from a `GridConfig`.
"""
function GridCoordinates(cfg::GridConfig; Nxmc::Int=4, Nymc::Int=4)
    return GridCoordinates(
        cfg.Nx, cfg.Ny; xsize=cfg.xsize, ysize=cfg.ysize, Nxmc=Nxmc, Nymc=Nymc
    )
end

"""
    default_grid_coordinates()

Construct the baseline `GridCoordinates` corresponding to compiled default constants.
"""
function default_grid_coordinates()
    return GridCoordinates(Nx, Ny; xsize=xsize, ysize=ysize, Nxmc=Nxmc, Nymc=Nymc)
end

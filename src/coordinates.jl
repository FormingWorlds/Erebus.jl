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
    xcenter::Union{Nothing,Float64}=nothing,
    ycenter::Union{Nothing,Float64}=nothing,
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
    xc = xcenter !== nothing ? xcenter : xsize / 2.0
    yc = ycenter !== nothing ? ycenter : ysize / 2.0

    x = collect(range(0.0, xsize; length=Nx))
    y = collect(range(0.0, ysize; length=Ny))
    xvx = collect(range(0.0, xsize + dx; length=Nx1))
    yvx = collect(range(-dy / 2.0, ysize + dy / 2.0; length=Ny1))
    xvy = collect(range(-dx / 2.0, xsize + dx / 2.0; length=Nx1))
    yvy = collect(range(0.0, ysize + dy; length=Ny1))
    xp = collect(range(-dx / 2.0, xsize + dx / 2.0; length=Nx1))
    yp = collect(range(-dy / 2.0, ysize + dy / 2.0; length=Ny1))

    Nxm = (Nx - 1) * Nxmc
    Nym = (Ny - 1) * Nymc
    dxm = xsize / Nxm
    dym = ysize / Nym
    start_marknum = Nxm * Nym
    xxm = collect(range(dxm / 2.0, xsize - dxm / 2.0; length=Nxm))
    yym = collect(range(dym / 2.0, ysize - dym / 2.0; length=Nym))

    return GridCoordinates(;
        Nx=Nx,
        Ny=Ny,
        Nx1=Nx1,
        Ny1=Ny1,
        xsize=xsize,
        ysize=ysize,
        dx=dx,
        dy=dy,
        xcenter=xc,
        ycenter=yc,
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
    GridCoordinates(cfg::GridConfig; xcenter=nothing, ycenter=nothing, Nxmc=4, Nymc=4)

Construct a `GridCoordinates` object from a `GridConfig`.
"""
function GridCoordinates(
    cfg::GridConfig;
    xcenter::Union{Nothing,Float64}=nothing,
    ycenter::Union{Nothing,Float64}=nothing,
    Nxmc::Int=4,
    Nymc::Int=4,
)
    return GridCoordinates(
        cfg.Nx,
        cfg.Ny;
        xsize=cfg.xsize,
        ysize=cfg.ysize,
        xcenter=xcenter,
        ycenter=ycenter,
        Nxmc=Nxmc,
        Nymc=Nymc,
    )
end

"""
    GridCoordinates(cfg::SimulationConfig; Nxmc=4, Nymc=4)

Construct a `GridCoordinates` object from a `SimulationConfig`.
"""
function GridCoordinates(cfg::SimulationConfig; Nxmc::Int=4, Nymc::Int=4)
    return GridCoordinates(
        cfg.grid;
        xcenter=cfg.geometry.xcenter,
        ycenter=cfg.geometry.ycenter,
        Nxmc=Nxmc,
        Nymc=Nymc,
    )
end

"""
    default_grid_coordinates()

Construct the baseline `GridCoordinates` with standard default grid dimensions.
"""
function default_grid_coordinates()
    return GridCoordinates(33, 33; xsize=140_000.0, ysize=140_000.0, Nxmc=4, Nymc=4)
end

"""
    @unpack_coords coords fields...

Unpack grid coordinate fields from `coords::GridCoordinates`.
Defines `<field>_val` for each symbol in `fields`.
"""
macro unpack_coords(coords, fields...)
    c_var = gensym("coords")
    assigns = map(fields) do f
        val = Symbol(f, "_val")
        return :($(esc(val)) = getfield($c_var, $(QuoteNode(f))))
    end
    return quote
        $c_var = $(esc(coords))
        $(assigns...)
    end
end

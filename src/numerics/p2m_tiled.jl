"""
Cell-tiled marker-to-mesh interpolation structures and routines.

Divides the domain into spatial tiles and uses a 4-color disjoint partitioning
to allow thread-parallel scatter directly into master grid arrays without
atomics, locks, or per-thread duplicate grids.
"""

"""
Workspace for cell-tiled marker-to-mesh interpolation.

Preallocates marker sorting and tile binning structures to achieve zero allocations
during simulation timestepping.

$(FIELDS)
"""
mutable struct P2MTiledWorkspace
    tile_size::Int
    Ncx::Int
    Ncy::Int
    Ntx::Int
    Nty::Int
    Ntiles::Int
    tile_counts::Vector{Int}
    tile_offsets::Vector{Int}
    tile_markers::Vector{Int}
    marker_tiles::Vector{Int}
    tiles_by_color::Vector{Vector{Int}}
end

"""
    P2MTiledWorkspace(coords, marknum, tile_size=4)

Allocate a new cell-tiled marker-to-mesh workspace for the given grid and marker count.
"""
function P2MTiledWorkspace(coords::GridCoordinates, marknum::Int, tile_size::Int=4)
    tile_size >= 2 || throw(ArgumentError("tile_size must be >= 2, got $tile_size"))
    Ncx = coords.Nx - 1
    Ncy = coords.Ny - 1
    Ntx = cld(Ncx, tile_size)
    Nty = cld(Ncy, tile_size)
    Ntiles = Ntx * Nty
    tile_counts = zeros(Int, Ntiles)
    tile_offsets = zeros(Int, Ntiles + 1)
    tile_markers = zeros(Int, marknum)
    marker_tiles = zeros(Int, marknum)
    tiles_by_color = [Int[] for _ in 1:4]
    for tj in 1:Ntx, ti in 1:Nty
        t = (tj - 1) * Nty + ti
        # 4-color parity pattern in (ti, tj)
        color = ((ti - 1) % 2) * 2 + ((tj - 1) % 2) + 1
        push!(tiles_by_color[color], t)
    end
    return P2MTiledWorkspace(
        tile_size,
        Ncx,
        Ncy,
        Ntx,
        Nty,
        Ntiles,
        tile_counts,
        tile_offsets,
        tile_markers,
        marker_tiles,
        tiles_by_color,
    )
end

"""
    bin_markers_into_tiles!(ws, xm, ym, coords, marknum)

Sort marker indices into tiles using a two-pass linear counting sort.
Reuses preallocated workspace buffers to avoid allocations.
"""
function bin_markers_into_tiles!(
    ws::P2MTiledWorkspace,
    xm::AbstractVector{<:Real},
    ym::AbstractVector{<:Real},
    coords::GridCoordinates,
    marknum::Int,
)
    if length(ws.tile_markers) < marknum
        resize!(ws.tile_markers, marknum)
        resize!(ws.marker_tiles, marknum)
    end

    fill!(ws.tile_counts, 0)
    x1 = coords.x[1]
    y1 = coords.y[1]
    inv_dx = 1.0 / coords.dx
    inv_dy = 1.0 / coords.dy
    Ncx = ws.Ncx
    Ncy = ws.Ncy
    Nty = ws.Nty
    ts = ws.tile_size

    # Pass 1: compute tile index for each marker and accumulate histogram
    @inbounds for m in 1:marknum
        cj = clamp(unsafe_trunc(Int, (xm[m] - x1) * inv_dx) + 1, 1, Ncx)
        ci = clamp(unsafe_trunc(Int, (ym[m] - y1) * inv_dy) + 1, 1, Ncy)
        tj = div(cj - 1, ts) + 1
        ti = div(ci - 1, ts) + 1
        t = (tj - 1) * Nty + ti
        ws.marker_tiles[m] = t
        ws.tile_counts[t] += 1
    end

    # Pass 2: prefix sum to compute offsets
    ws.tile_offsets[1] = 1
    @inbounds for t in 1:ws.Ntiles
        ws.tile_offsets[t + 1] = ws.tile_offsets[t] + ws.tile_counts[t]
    end

    # Pass 3: populate sorted marker indices
    copyto!(ws.tile_counts, 1, ws.tile_offsets, 1, ws.Ntiles)
    @inbounds for m in 1:marknum
        t = ws.marker_tiles[m]
        pos = ws.tile_counts[t]
        ws.tile_markers[pos] = m
        ws.tile_counts[t] = pos + 1
    end

    return nothing
end

"""
    scatter_marker_to_master_grids!(m, xm, ym, coords, ...)

Scatter properties of marker `m` to Basic, Vx, Vy, and P master nodes.
"""
@inline function scatter_marker_to_master_grids!(
    m::Int,
    xmm::Real,
    ymm::Real,
    coords::GridCoordinates,
    etatotalm,
    etavpm,
    inv_gggtotalm,
    sxym,
    cohestotalm,
    tenstotalm,
    fricttotalm,
    ETA0SUM,
    ETASUM,
    GGGSUM,
    SXYSUM,
    COHSUM,
    TENSUM,
    FRISUM,
    WTSUM,
    rhototalm,
    rhofluidcur,
    ktotalm,
    phim,
    etafluidcur_inv_kphim,
    RHOXSUM,
    RHOFXSUM,
    KXSUM,
    PHIXSUM,
    RXSUM,
    WTXSUM,
    RHOYSUM,
    RHOFYSUM,
    KYSUM,
    PHIYSUM,
    RYSUM,
    WTYSUM,
    sxxm,
    rhocptotalm,
    alphasolidcur,
    alphafluidcur,
    hrtotalm,
    tkm_rhocptotalm,
    GGGPSUM,
    SXXSUM,
    RHOSUM,
    RHOCPSUM,
    ALPHASUM,
    ALPHAFSUM,
    HRSUM,
    PHISUM,
    TKSUM,
    WTPSUM,
)
    @inbounds marker_to_basic_nodes!(
        m,
        xmm,
        ymm,
        etatotalm,
        etavpm,
        inv_gggtotalm,
        sxym,
        cohestotalm,
        tenstotalm,
        fricttotalm,
        ETA0SUM,
        ETASUM,
        GGGSUM,
        SXYSUM,
        COHSUM,
        TENSUM,
        FRISUM,
        WTSUM;
        coords=coords,
    )
    @inbounds marker_to_vx_nodes!(
        m,
        xmm,
        ymm,
        rhototalm,
        rhofluidcur,
        ktotalm,
        phim,
        etafluidcur_inv_kphim,
        RHOXSUM,
        RHOFXSUM,
        KXSUM,
        PHIXSUM,
        RXSUM,
        WTXSUM;
        coords=coords,
    )
    @inbounds marker_to_vy_nodes!(
        m,
        xmm,
        ymm,
        rhototalm,
        rhofluidcur,
        ktotalm,
        phim,
        etafluidcur_inv_kphim,
        RHOYSUM,
        RHOFYSUM,
        KYSUM,
        PHIYSUM,
        RYSUM,
        WTYSUM;
        coords=coords,
    )
    @inbounds marker_to_p_nodes!(
        m,
        xmm,
        ymm,
        inv_gggtotalm,
        sxxm,
        rhototalm,
        rhocptotalm,
        alphasolidcur,
        alphafluidcur,
        hrtotalm,
        phim,
        tkm_rhocptotalm,
        GGGPSUM,
        SXXSUM,
        RHOSUM,
        RHOCPSUM,
        ALPHASUM,
        ALPHAFSUM,
        HRSUM,
        PHISUM,
        TKSUM,
        WTPSUM;
        coords=coords,
    )
    return nothing
end

"""
    ensure_workspace_compatible(ws, coords, marknum, tile_size=4)

Ensure that `ws` is allocated and compatible with the current `coords` and `marknum`.
Re-allocates or resizes when coordinates or marker counts change.
"""
function ensure_workspace_compatible(
    ws::Union{Nothing,P2MTiledWorkspace},
    coords::GridCoordinates,
    marknum::Int,
    tile_size::Int=4,
)
    if ws === nothing || (coords.Nx - 1 != ws.Ncx) || (coords.Ny - 1 != ws.Ncy)
        return P2MTiledWorkspace(coords, marknum, tile_size)
    end
    if length(ws.tile_markers) < marknum
        resize!(ws.tile_markers, marknum)
        resize!(ws.marker_tiles, marknum)
    end
    return ws
end

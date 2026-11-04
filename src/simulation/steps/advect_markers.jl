# Marker advection and pressure backtracking step

"""
    advect_markers!(state::SimulationState, coords::GridCoordinates, cfg::SimulationConfig; XWSSUM=nothing, WTPSUM=nothing)::Nothing

Advect Lagrangian markers with 4th-order Runge-Kutta scheme and backtrace nodal pressures.

# Mutates:
- `state.markers.core.phinewm`
- `state.markers.core.xm`
- `state.markers.core.ym`
- `state.markers.core.sxym`
- `state.markers.core.sxxm`
- `state.grids.XWS`
- `state.grids.vxp`
- `state.grids.vyp`
- `state.grids.vxpf`
- `state.grids.vypf`
- `state.grids.wyx`
- `state.grids.pr`
- `state.grids.pr0`
- `state.grids.ps`
- `state.grids.ps0`
- `state.grids.pf`
- `state.grids.pf0`
"""
function advect_markers!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig;
    XWSSUM=nothing,
    WTPSUM=nothing,
)::Nothing
    marknum = length(state.markers)
    marknum == 0 && return nothing

    markers = state.markers
    grids = state.grids

    xm = markers.core.xm
    ym = markers.core.ym
    phinewm = markers.core.phinewm
    sxym = markers.core.sxym
    sxxm = markers.core.sxxm
    XWsolidm0 = markers.core.XWsolidm0

    vx = grids.vx
    vy = grids.vy
    phim = markers.core.phim
    vxf = grids.vxf
    vyf = grids.vyf
    vxp = grids.vxp
    vyp = grids.vyp
    vxpf = grids.vxpf
    vypf = grids.vypf
    wyx = grids.wyx
    pr = grids.pr
    pr0 = grids.pr0
    ps = grids.ps
    ps0 = grids.ps0
    pf = grids.pf
    pf0 = grids.pf0
    XWS = grids.XWS

    dt = state.dt

    # 1. Synchronize porosity state
    phinewm .= phim

    # 2. Interpolate melt composition from markers to P nodes
    xws_sum = XWSSUM !== nothing ? XWSSUM : zeros(Float64, coords.Ny1, coords.Nx1)
    wtp_sum = WTPSUM !== nothing ? WTPSUM : zeros(Float64, coords.Ny1, coords.Nx1)
    update_p_nodes_melt_composition!(
        xm, ym, XWsolidm0, XWS, xws_sum, wtp_sum, marknum; coords=coords
    )

    # 3. Compute velocities in P nodes
    compute_velocities!(vx, vy, vxf, vyf, vxp, vyp, vxpf, vypf; coords=coords)

    # 4. Compute rotation rate in basic nodes
    compute_rotation_rate!(vx, vy, wyx; coords=coords)

    # 5. Move markers with 4th-order Runge-Kutta scheme
    move_markers_rk4!(xm, ym, sxym, sxxm, vx, vy, wyx, marknum, dt; coords=coords)

    # 6. Backtrack nodal total and fluid pressures
    backtrace_pressures_rk4!(pr, pr0, ps, ps0, pf, pf0, vx, vy, vxf, vyf, dt; coords=coords)

    return nothing
end

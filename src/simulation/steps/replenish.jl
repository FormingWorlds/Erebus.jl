# Marker replenishment and geometric weight calculation step

"""
    replenish!(
        state::SimulationState,
        coords::GridCoordinates,
        cfg::SimulationConfig;
        mdis=nothing,
        mnum=nothing,
        randomized::Bool=random_markers,
        step_start_buffers=nothing,
    )::Int

Replenish depleted grid cells with new markers and recalculate out-of-plane geometric weights.

# Mutates:
- `state.markers`
- `state.markers.core.w3d_m`

# Returns:
- `Int`: Updated total marker count.
"""
function replenish!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig;
    mdis=nothing,
    mnum=nothing,
    randomized::Bool=random_markers,
    step_start_buffers=nothing,
)::Int
    marknum = length(state.markers)
    marknum == 0 && return 0

    markers = state.markers
    acc = state.accumulators
    rng = state.rng
    xcenter_val = acc.xcenter
    ycenter_val = acc.ycenter

    mdis_buf, mnum_buf = if mdis !== nothing && mnum !== nothing
        (mdis, mnum)
    else
        helpers = setup_marker_geometry_helpers(coords)
        (mdis !== nothing ? mdis : helpers[1], mnum !== nothing ? mnum : helpers[2])
    end

    new_marknum = replenish_markers!(
        markers, mdis_buf, mnum_buf; randomized=randomized, coords=coords, cfg=cfg, rng=rng
    )

    if step_start_buffers !== nothing
        for buf in step_start_buffers
            if buf !== nothing && length(buf) != new_marknum
                resize!(buf, new_marknum)
            end
        end
    end

    resize!(markers.core.w3d_m, new_marknum)
    xm = markers.core.xm
    ym = markers.core.ym
    for m in 1:new_marknum
        markers.core.w3d_m[m] = marker_out_of_plane_length(
            xm[m], ym[m], xcenter_val, ycenter_val
        )
    end

    assert_marker_arrays_invariants(markers, new_marknum)
    return new_marknum
end

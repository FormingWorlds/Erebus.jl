# Thermomechanical iteration driver coordinating Stokes-Darcy and thermal energy solves

"""
Restore marker redox speciation and bulk metal fraction from step-start snapshots on outer iterations.

$(SIGNATURES)

# Parameters
- `state`: Simulation state container.
- `cfg`: Simulation configuration struct.
- `ws`: Preallocated simulation workspaces.

# Returns
- `Nothing`: Restores marker redox vectors and `Xfe_bulk` in place.
"""
function _restore_step_start_redox!(
    state::SimulationState, cfg::SimulationConfig, ws::SimulationWorkspaces
)::Nothing
    if cfg.redox.active && haskey(state.markers.groups, :redox)
        redox = state.markers.groups.redox
        _sync_step_start_vec!(redox.nFe0_m, ws.nFe0_m_step_start)
        _sync_step_start_vec!(redox.nFe2_m, ws.nFe2_m_step_start)
        _sync_step_start_vec!(redox.nFe3_m, ws.nFe3_m_step_start)
        _sync_step_start_vec!(redox.deltaIW_m, ws.deltaIW_m_step_start)
        if haskey(state.markers.groups, :metal)
            _sync_step_start_vec!(
                state.markers.groups.metal.Xfe_bulk, ws.Xfe_bulk_step_start
            )
        end
    end
    return nothing
end

"""
Execute one outer iteration of fluid-silicate thermochemical hydration/dehydration or zero source grids.

$(SIGNATURES)

# Parameters
- `state`: Simulation state container.
- `coords`: Grid coordinates struct.
- `cfg`: Simulation configuration struct.
- `ws`: Preallocated simulation workspaces.
- `titer`: Current outer thermomechanical iteration index.
- `dt_prev`: Time step of the previous reaction call; the previous iterate is scaled to `state.dt`.

# Returns
- `Nothing`: Mutates `state.grids.DMP`, `state.grids.DHP`, and `state.grids.DQPF` in place.
"""
function _step_thermochemical_reaction!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig,
    ws::SimulationWorkspaces,
    titer::Int;
    dt_prev::Real=state.dt,
)::Nothing
    if cfg.reaction.active
        core = state.markers.core
        copyto!(core.phinewm, core.phim)
        DMPSUM = ws.interp_arrays[30]
        DHPSUM = ws.interp_arrays[31]
        WTPSUM = ws.interp_arrays[33]
        perform_thermochemical_reaction!(
            state.grids.DMP,
            state.grids.DHP,
            DMPSUM,
            DHPSUM,
            WTPSUM,
            state.grids.pf,
            state.grids.tk2,
            core.tm,
            core.xm,
            core.ym,
            core.XWsolidm0,
            core.XWsolidm,
            core.phim,
            core.phinewm,
            core.pfm0,
            length(state.markers),
            state.dt,
            state.timestep,
            titer;
            coords=coords,
            DQPF=state.grids.DQPF,
            DQPFSUM=state.grids.DQPFSUM,
            cfg=cfg.reaction,
            dt_prev=dt_prev,
        )
    else
        fill!(state.grids.DMP, 0.0)
        fill!(state.grids.DHP, 0.0)
        fill!(state.grids.DQPF, 0.0)
    end
    return nothing
end

"""
Swap contents of `vec` and `step_start_vec` in place when both vectors are non-nothing.

$(SIGNATURES)

Parameters
----------
vec : Union{Nothing,AbstractVector}
    Active marker property vector or `nothing`.
step_start_vec : Union{Nothing,AbstractVector}
    Step-start snapshot vector or `nothing`.

Returns
-------
Nothing
    Swaps elements of `vec` and `step_start_vec` in place.

Raises
------
DimensionMismatch
    If `vec` and `step_start_vec` have different lengths.
"""
@inline function _swap_step_start_vec!(vec, step_start_vec)::Nothing
    (vec === nothing || step_start_vec === nothing) && return nothing
    if length(vec) != length(step_start_vec)
        throw(
            DimensionMismatch(
                "Marker vector length $(length(vec)) does not match step-start snapshot length $(length(step_start_vec))",
            ),
        )
    end
    @inbounds for m in eachindex(vec, step_start_vec)
        vec[m], step_start_vec[m] = step_start_vec[m], vec[m]
    end
    return nothing
end

"""
Recompute step-averaged marker radiogenic heating and re-project `HR` onto P-nodes for a new `state.dt`.

$(SIGNATURES)

Parameters
----------
state : SimulationState
    Mutable simulation state container.
coords : GridCoordinates
    Staggered grid coordinate specification.
cfg : SimulationConfig
    Master simulation configuration.
ws : SimulationWorkspaces
    Preallocated solver workspaces and step-start buffers.

Returns
-------
Nothing
    Updates `state.markers.core.hrtotalm` and `state.grids.HR` in place using pre-pyrolysis step-start `phim` and `Xfe_bulk`.
"""
function _update_radiogenic_heating_for_dt!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig,
    ws::SimulationWorkspaces,
)::Nothing
    marknum = length(state.markers)
    marknum == 0 && return nothing
    (cfg.thermodynamics.hr_al || cfg.thermodynamics.hr_fe) || return nothing
    if ws.hr_step_start_timestep != state.timestep || ws.phim_hr_step_start === nothing
        throw(
            ArgumentError(
                "Radiogenic step-start state was not snapshotted for timestep $(state.timestep) (hr_step_start_timestep=$(ws.hr_step_start_timestep))",
            ),
        )
    end

    core = state.markers.core
    Xfe_bulk = if haskey(state.markers.groups, :metal)
        state.markers.groups.metal.Xfe_bulk
    else
        nothing
    end

    _swap_step_start_vec!(core.phim, ws.phim_hr_step_start)
    try
        _swap_step_start_vec!(Xfe_bulk, ws.Xfe_bulk_hr_step_start)
        try
            radiogenic_heating!(state, coords, cfg)
        finally
            _swap_step_start_vec!(Xfe_bulk, ws.Xfe_bulk_hr_step_start)
        end
    finally
        _swap_step_start_vec!(core.phim, ws.phim_hr_step_start)
    end

    HRSUM = ws.interp_arrays[25]
    WTPSUM = ws.interp_arrays[33]
    fill!(HRSUM, 0.0)
    fill!(WTPSUM, 0.0)

    @unpack_coords coords xp yp dx dy jmin_p jmax_p imin_p imax_p
    xm = core.xm
    ym = core.ym
    hrtotalm = core.hrtotalm

    if cfg.solver.p2m_mode == :tiled && ws.p2m !== nothing
        ws_p2m = ws.p2m
        if ws_p2m.Ncx != coords.Nx - 1 ||
            ws_p2m.Ncy != coords.Ny - 1 ||
            length(ws_p2m.tile_markers) < marknum ||
            length(ws_p2m.tile_offsets) != ws_p2m.Ntiles + 1 ||
            ws_p2m.tile_offsets[ws_p2m.Ntiles + 1] - 1 != marknum
            throw(
                DimensionMismatch(
                    "Tiled P2M workspace binning does not match grid ($(coords.Ny)x$(coords.Nx)) or marker count ($marknum)",
                ),
            )
        end
        @inbounds for color in 1:4
            for t in ws_p2m.tiles_by_color[color]
                lo = ws_p2m.tile_offsets[t]
                hi = ws_p2m.tile_offsets[t + 1] - 1
                for idx in lo:hi
                    m = ws_p2m.tile_markers[idx]
                    i, j, weights = fix_weights(
                        xm[m],
                        ym[m],
                        xp_val,
                        yp_val,
                        dx_val,
                        dy_val,
                        jmin_p_val,
                        jmax_p_val,
                        imin_p_val,
                        imax_p_val,
                    )
                    interpolate_add_to_grid!(i, j, weights, hrtotalm[m], HRSUM)
                    interpolate_add_to_grid!(i, j, weights, 1.0, WTPSUM)
                end
            end
        end
    else
        @inbounds for m in 1:marknum
            i, j, weights = fix_weights(
                xm[m],
                ym[m],
                xp_val,
                yp_val,
                dx_val,
                dy_val,
                jmin_p_val,
                jmax_p_val,
                imin_p_val,
                imax_p_val,
            )
            interpolate_add_to_grid!(i, j, weights, hrtotalm[m], HRSUM)
            interpolate_add_to_grid!(i, j, weights, 1.0, WTPSUM)
        end
    end

    HR = state.grids.HR
    Ny1, Nx1 = size(WTPSUM)
    @inbounds for j in 1:Nx1, i in 1:Ny1
        wt = WTPSUM[i, j]
        if wt > 0.0
            HR[i, j] = HRSUM[i, j] * inv(wt)
        end
    end
    fill_zero_weight_nodes!(HR, WTPSUM; background_val=0.0)
    return nothing
end

"""
Coordinate thermochemical iterations, Stokes-Darcy iterations, and thermal energy updates.

$(SIGNATURES)
"""
function solve_thermomechanical_iterations!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig,
    ws::SimulationWorkspaces;
    DHP_pyro=nothing,
    dt_step_initial::Float64,
)
    ws.hydromech_cache = nothing
    ws.thermal_cache = nothing
    dt_aphimax_step_max = 0.0
    plastic_converged = true
    thermochemical_converged = true
    last_plastic_residual = 0.0
    dt_reduced_by_maxDT = false
    dt_next = state.dt
    dt_hr = dt_step_initial
    n_flips_last = 0
    n_flips_step_total = 0

    max_iters = cfg.solver.max_plastic_iterations
    dt_rx = state.dt
    n_titers = 0
    for titer in 1:max_iters
        n_titers = titer
        titer > 1 && _restore_step_start_redox!(state, cfg, ws)
        _step_thermochemical_reaction!(state, coords, cfg, ws, titer; dt_prev=dt_rx)
        dt_rx = state.dt
        if DHP_pyro !== nothing && size(DHP_pyro) == size(state.grids.DHP)
            state.grids.DHP .+= DHP_pyro
        end

        sd_res = solve_stokes_darcy!(
            state, coords, cfg, ws; titer=titer, dt_step_initial=dt_step_initial
        )
        dt_aphimax_step_max = max(dt_aphimax_step_max, sd_res.dt_aphimax_step_max)
        n_flips_last = sd_res.n_flips_last
        n_flips_step_total += sd_res.n_flips_step_total

        if !sd_res.plastic_converged
            plastic_converged = false
            last_plastic_residual = sd_res.last_plastic_residual
            break
        end

        if state.dt != dt_hr
            _update_radiogenic_heating_for_dt!(state, coords, cfg, ws)
            dt_hr = state.dt
        end
        th_res = solve_thermal_energy!(
            state, coords, cfg, ws; titer=titer, DHP_pyro=DHP_pyro
        )
        dt_next = th_res.dt_next
        dt_reduced_by_maxDT = dt_reduced_by_maxDT || th_res.dt_reduced_by_maxDT

        # A reaction from before a large dt cut must be recomputed before the step is accepted
        reaction_stale = cfg.reaction.active && dt_rx > 2.0 * state.dt
        if th_res.thermochemical_converged && !reaction_stale
            thermochemical_converged = true
            break
        elseif titer == max_iters
            thermochemical_converged = false
            break
        else
            state.dt = dt_next
        end
    end

    if plastic_converged
        update_marker_pressures!(
            state.markers.core.pfm0,
            state.markers.core.xm,
            state.markers.core.ym,
            state.markers.core.tm,
            state.grids.pf,
            state.grids.pr,
            length(state.markers);
            coords=coords,
            p_cavitation=cfg.reaction.p_cavitation,
            pfcoeff=0.0,
        )
    end

    return (;
        plastic_converged=plastic_converged,
        last_plastic_residual=last_plastic_residual,
        thermochemical_converged=thermochemical_converged,
        dt_next=dt_next,
        dt_reduced_by_maxDT=dt_reduced_by_maxDT,
        dt_aphimax_step_max=dt_aphimax_step_max,
        n_flips_last=n_flips_last,
        n_flips_step_total=n_flips_step_total,
        n_titers=n_titers,
    )
end

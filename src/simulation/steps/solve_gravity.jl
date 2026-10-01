# Gravitational potential and acceleration solve

"""
    solve_gravity!(state::SimulationState, coords::GridCoordinates, cfg::SimulationConfig; F_grav=nothing, RP=nothing, SP=nothing)::Nothing

Solve gravitational potential and acceleration on staggered grid nodes.

# Mutates:
- `state.grids.FI`
- `state.grids.gx`
- `state.grids.gy`
"""
function solve_gravity!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig;
    F_grav=nothing,
    RP=nothing,
    SP=nothing,
)::Nothing
    # Empty marker set contract: return immediately with no mutations
    length(state.markers) == 0 && return nothing

    grids = state.grids
    gx = grids.gx
    gy = grids.gy
    FI = grids.FI
    mode = cfg.geometry.gravity_mode

    if mode === :enclosed_mass
        compute_gravity_enclosed_mass!(
            gx,
            gy;
            xm=state.markers.xm,
            ym=state.markers.ym,
            rhototalm=state.markers.rhototalm,
            tm=state.markers.tm,
            coords=coords,
            gravity_nr_factor=cfg.geometry.gravity_nr_factor,
            rplanet=cfg.geometry.rplanet,
            FI=FI,
        )
    elseif mode === :poisson2d
        RP_vec, SP_vec = if RP !== nothing && SP !== nothing
            RP, SP
        else
            setup_gravitational_lse(coords)
        end
        F_fact = if F_grav !== nothing
            F_grav
        else
            LP = assemble_gravitational_lse!(
                zeros(coords.Ny1, coords.Nx1), RP_vec; coords=coords
            )
            lu(LP.cscmatrix)
        end
        assemble_gravitational_rhs!(grids.RHO, RP_vec; coords=coords)
        SP_vec = F_fact \ RP_vec
        process_gravitational_solution!(SP_vec, FI, gx, gy; coords=coords)
    else
        throw(
            ArgumentError(
                "Unknown gravity_mode: $(mode). Must be :enclosed_mass or :poisson2d"
            ),
        )
    end

    return nothing
end

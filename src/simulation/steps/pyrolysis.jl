# Refractory organic devolatilization and pyrolysis step

"""
Evaluate refractory organic devolatilization and compute pyrolysis enthalpy sink.

$(SIGNATURES)
"""
function update_pyrolysis!(
    state::SimulationState, coords::GridCoordinates, cfg::SimulationConfig
)::Union{Nothing,Matrix{Float64}}
    hcnspo_props = get(state.markers.groups, :hcnspo, nothing)
    if !cfg.refractory.active || !cfg.refractory.kinetics_active || hcnspo_props === nothing
        return nothing
    end
    core = state.markers.core
    metal = get(state.markers.groups, :metal, nothing)
    redox_props = get(state.markers.groups, :redox, nothing)
    vols = get(state.markers.groups, :volatiles, nothing)
    dhp_p = zeros(Float64, coords.Ny1, coords.Nx1)
    phim_target = cfg.reaction.active ? core.phinewm : core.phim
    update_marker_pyrolysis!(
        core.tkm,
        state.dt,
        phim_target,
        hcnspo_props.X_refr_C_m,
        hcnspo_props.X_refr_N_m,
        hcnspo_props.X_refr_H_m,
        cfg.refractory;
        tm=core.tm,
        rhosolidm=cfg.materials.rhosolidm,
        xm=core.xm,
        ym=core.ym,
        coords=coords,
        DHP=dhp_p,
        redox_props=redox_props,
        redox_cfg=cfg.redox,
        Xfem=(metal !== nothing ? metal.Xfem : nothing),
        Xfe_bulk=(metal !== nothing ? metal.Xfe_bulk : nothing),
        rho_metal=cfg.coreformation.rho_metal,
        XH2Om=(vols !== nothing ? vols.XH2Om : nothing),
        XCm=(vols !== nothing ? vols.XCm : nothing),
        XNm=(vols !== nothing ? vols.XNm : nothing),
        XSm=(vols !== nothing ? vols.XSm : nothing),
        X_graphite_m=(vols !== nothing ? vols.X_graphite_m : nothing),
    )
    copyto!(cfg.reaction.active ? core.phim : core.phinewm, phim_target)
    state.grids.DHP .+= dhp_p
    return dhp_p
end

# Refractory organic devolatilization and pyrolysis step

"""
Evaluate refractory organic devolatilization and compute pyrolysis enthalpy sink.

$(SIGNATURES)
"""
function update_pyrolysis!(
    state::SimulationState, coords::GridCoordinates, cfg::SimulationConfig
)::Union{Nothing,Matrix{Float64}}
    hcnspo_props = haskey(state.markers.groups, :hcnspo) ? state.markers.groups.hcnspo : nothing
    if !cfg.refractory.active || !cfg.refractory.kinetics_active || hcnspo_props === nothing
        return nothing
    end

    core = state.markers.core
    metal = haskey(state.markers.groups, :metal) ? state.markers.groups.metal : nothing
    redox_props = haskey(state.markers.groups, :redox) ? state.markers.groups.redox : nothing

    Xfem = metal !== nothing ? metal.Xfem : nothing
    Xfe_bulk = metal !== nothing ? metal.Xfe_bulk : nothing
    dhp_p = zeros(Float64, coords.Ny1, coords.Nx1)

    update_marker_pyrolysis!(
        core.tkm,
        state.dt,
        core.phim,
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
        Xfem=Xfem,
        Xfe_bulk=Xfe_bulk,
        rho_metal=cfg.coreformation.rho_metal,
    )

    return dhp_p
end

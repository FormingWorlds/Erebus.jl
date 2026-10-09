# Simulation initialization and domain setup step

"""
Unpack core marker array references from CoreGroup.

$(SIGNATURES)
"""
function _unpack_core_arrays(core::CoreGroup)
    return (
        core.xm,
        core.ym,
        core.w3d_m,
        core.tm,
        core.tkm,
        core.phim,
        core.phinewm,
        core.pfm0,
        core.XWsolidm,
        core.XWsolidm0,
        core.Fm,
        core.etavpm,
        core.sxxm,
        core.sxym,
        core.inv_gggtotalm,
        core.fricttotalm,
        core.cohestotalm,
        core.tenstotalm,
        core.rhototalm,
        core.rhocptotalm,
        core.etatotalm,
        core.hrtotalm,
        core.ktotalm,
        core.tkm_rhocptotalm,
        core.etafluidcur_inv_kphim,
        core.rhofluidcur,
        core.alphasolidcur,
        core.alphafluidcur,
    )
end

"""
Extract optional marker group arrays from MarkerArrays container.

$(SIGNATURES)
"""
function _extract_optional_marker_arrays(
    markers::MarkerArrays, cfg::SimulationConfig, marknum::Int, magma_active_val::Bool
)
    metal = haskey(markers.groups, :metal) ? markers.groups.metal : nothing
    volatiles = haskey(markers.groups, :volatiles) ? markers.groups.volatiles : nothing
    phase = haskey(markers.groups, :phase) ? markers.groups.phase : nothing
    accretion = haskey(markers.groups, :accretion) ? markers.groups.accretion : nothing

    Xfem = metal !== nothing ? metal.Xfem : nothing
    Xfem0 = metal !== nothing ? metal.Xfem0 : nothing
    Xfe_bulk = metal !== nothing ? metal.Xfe_bulk : nothing
    Xfe_bulk_step_start = Xfe_bulk !== nothing ? zeros(Float64, marknum) : nothing
    Xfem_step_start = Xfem !== nothing ? zeros(Float64, marknum) : nothing

    F_extract_m = if volatiles !== nothing
        volatiles.F_extract_m
    else
        (magma_active_val ? setup_marker_magma_properties(marknum)[1] : nothing)
    end
    F_extract_m_step_start = F_extract_m !== nothing ? zeros(Float64, marknum) : nothing
    Fm_step_start = magma_active_val ? zeros(Float64, marknum) : nothing

    XH2Om = volatiles !== nothing ? volatiles.XH2Om : nothing
    XCm = volatiles !== nothing ? volatiles.XCm : nothing
    XNm = volatiles !== nothing ? volatiles.XNm : nothing
    XSm = volatiles !== nothing ? volatiles.XSm : nothing
    X_graphite_m = volatiles !== nothing ? volatiles.X_graphite_m : nothing

    has_metal_part = metal !== nothing && cfg.metal_partition.active
    Xfe_H_m = has_metal_part ? metal.Xfe_H_m : nothing
    Xfe_C_m = has_metal_part ? metal.Xfe_C_m : nothing
    Xfe_N_m = has_metal_part ? metal.Xfe_N_m : nothing
    Xfe_S_m = has_metal_part ? metal.Xfe_S_m : nothing
    Xfe_H_m_step_start = Xfe_H_m !== nothing ? zeros(Float64, marknum) : nothing
    Xfe_C_m_step_start = Xfe_C_m !== nothing ? zeros(Float64, marknum) : nothing
    Xfe_N_m_step_start = Xfe_N_m !== nothing ? zeros(Float64, marknum) : nothing
    Xfe_S_m_step_start = Xfe_S_m !== nothing ? zeros(Float64, marknum) : nothing

    Xmin_troilite_m = phase !== nothing ? phase.Xmin_troilite_m : nothing
    Xmin_schreibersite_m = phase !== nothing ? phase.Xmin_schreibersite_m : nothing
    Xmin_cohenite_m = phase !== nothing ? phase.Xmin_cohenite_m : nothing
    Xmin_graphite_m = phase !== nothing ? phase.Xmin_graphite_m : nothing
    Xmin_nitride_m = phase !== nothing ? phase.Xmin_nitride_m : nothing
    Xmin_metal_matrix_m = phase !== nothing ? phase.Xmin_metal_matrix_m : nothing

    t_accreted = accretion !== nothing ? accretion.t_accreted : nothing
    hcnspo_props = haskey(markers.groups, :hcnspo) ? markers.groups.hcnspo : nothing
    redox_props = haskey(markers.groups, :redox) ? markers.groups.redox : nothing

    return (;
        Xfem,
        Xfem0,
        Xfe_bulk,
        Xfe_bulk_step_start,
        Xfem_step_start,
        F_extract_m,
        F_extract_m_step_start,
        Fm_step_start,
        XH2Om,
        XCm,
        XNm,
        XSm,
        X_graphite_m,
        Xfe_H_m,
        Xfe_C_m,
        Xfe_N_m,
        Xfe_S_m,
        Xfe_H_m_step_start,
        Xfe_C_m_step_start,
        Xfe_N_m_step_start,
        Xfe_S_m_step_start,
        Xmin_troilite_m,
        Xmin_schreibersite_m,
        Xmin_cohenite_m,
        Xmin_graphite_m,
        Xmin_nitride_m,
        Xmin_metal_matrix_m,
        t_accreted,
        hcnspo_props,
        redox_props,
    )
end

"""
Set up simulation output path and grid coordinates.

$(SIGNATURES)
"""
function setup_simulation_grid_and_coords(cfg::SimulationConfig; output_path::String)
    formatted_path = endswith(output_path, "/") ? output_path : output_path * "/"
    isdir(formatted_path) || mkpath(formatted_path)
    coords = GridCoordinates(cfg.grid)
    return coords, formatted_path
end

"""
Set up or restore simulation state, markers, and grid arrays.

$(SIGNATURES)
"""
function setup_simulation_markers_and_state(
    cfg::SimulationConfig,
    coords::GridCoordinates;
    output_path::String,
    restart_from::AbstractString,
    force_restart_config::Bool,
)
    timestep, dt, timesum, marknum, hrsolidm, hrfluidm, YERRNOD = setup_dynamic_simulation_parameters(
        cfg; coords=coords
    )

    rplanet_val = cfg.accretion.active ? cfg.accretion.R_initial : cfg.geometry.rplanet
    rcrust_val = if cfg.accretion.active
        min(cfg.geometry.rcrust, cfg.accretion.R_initial)
    else
        cfg.geometry.rcrust
    end
    M_planet_val = if cfg.accretion.active
        cfg.accretion.M_initial
    elseif cfg.escape.active
        cfg.escape.M_planet
    else
        (4.0 / 3.0 * pi * (rplanet_val^3) * cfg.accretion.rho_bulk)
    end
    xcenter_val = cfg.geometry.xcenter
    ycenter_val = cfg.geometry.ycenter

    rng = MersenneTwister(cfg.solver.seed)
    timer = TimerOutput()
    transfer_log = TransferRecord[]

    atm_state = if cfg.atmosphere.active
        AtmosphereState(
            Dict{Symbol,Float64}(sp => 0.0 for sp in cfg.escape.species_list),
            Dict{Symbol,Float64}(sp => 0.0 for sp in cfg.escape.species_list),
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        )
    else
        nothing
    end

    is_restart = !isempty(restart_from)
    start_step_val = cfg.time.start_step
    coords_actual = coords

    if is_restart
        (ckpt_state, ckpt_coords, cfg_saved) = load_simulation_state(
            restart_from; force_restart_config=force_restart_config, current_cfg=cfg
        )
        coords_actual = ckpt_coords
        grids = copy(ckpt_state.grids)
        markers = copy(ckpt_state.markers)
        accumulators = copy(ckpt_state.accumulators)
        transfers = deepcopy(ckpt_state.transfers)
        atm_state = ckpt_state.atm !== nothing ? copy(ckpt_state.atm) : atm_state
        rng = copy(ckpt_state.rng)
        timer = copy(ckpt_state.timer)
        start_step_val = ckpt_state.timestep + 1
        dt = ckpt_state.dt
        timesum = ckpt_state.timesum

        state = SimulationState(
            grids,
            markers,
            accumulators,
            transfers,
            atm_state,
            rng,
            timer,
            ckpt_state.timestep,
            dt,
            timesum,
        )
        @info "Resumed simulation from checkpoint: $restart_from at timestep $(start_step_val - 1)"
    else
        (
            ETA, ETA0, GGG, EXY, SXY, SXY0, wyx, COH, TEN, FRI, YNY,
            RHOX, RHOFX, KX, PHIX, vx, vxf, RX, qxD, gx,
            RHOY, RHOFY, KY, PHIY, vy, vyf, RY, qyD, gy,
            RHO, RHOCP, ALPHA, ALPHAF, HR, HA, HS, ETAP, GGGP,
            EXX, SXX, SXX0, tk1, tk2, DT, DT0, vxp, vyp, vxpf, vypf,
            pr, pf, ps, pr0, pf0, ps0, ETAPHI, BETAPHI, PHI, APHI, FI, DMP, DHP, XWS
        ) = setup_staggered_grid_properties(coords; rng=rng)

        (ETA5, ETA00, YNY5, YNY00, YNY_inv_ETA, DSXY, DSY, EII, SII, DSXX, tk0) = setup_staggered_grid_properties_helpers(
            coords; rng=rng
        )

        Q_metric = cfg.geometry.spherical_metric ? zeros(Float64, coords.Ny1, coords.Nx1) : nothing
        DQPF = zeros(Float64, coords.Ny1, coords.Nx1)
        DQPFSUM = zeros(Float64, coords.Ny1, coords.Nx1)
        S_vent_grid = zeros(Float64, coords.Ny1, coords.Nx1)
        Q_lat_grid = zeros(Float64, coords.Ny1, coords.Nx1)
        Q_seg_grid = zeros(Float64, coords.Ny1, coords.Nx1)

        grids = GridArrays(
            ETA, ETA0, GGG, EXY, SXY, SXY0, wyx, COH, TEN, FRI, YNY,
            RHOX, RHOFX, KX, PHIX, vx, vxf, RX, qxD, gx,
            RHOY, RHOFY, KY, PHIY, vy, vyf, RY, qyD, gy,
            RHO, RHOCP, ALPHA, ALPHAF, HR, HA, HS, ETAP, GGGP,
            EXX, SXX, SXX0, tk1, tk2, DT, DT0, vxp, vyp, vxpf, vypf,
            pr, pf, ps, pr0, pf0, ps0, ETAPHI, BETAPHI, PHI, APHI, FI, DMP, DHP, XWS,
            ETA5, ETA00, YNY5, YNY00, YNY_inv_ETA, DSXY, DSY, EII, SII, DSXX, tk0,
            DQPF, DQPFSUM, S_vent_grid, Q_lat_grid, Q_seg_grid, Q_metric,
        )

        markers = init_marker_arrays(marknum, cfg, coords; rng=rng, initial_time=timesum)
        core = markers.core
        (;
            Xfem, Xfem0, Xfe_bulk, Xfe_bulk_step_start, Xfem_step_start,
            F_extract_m, F_extract_m_step_start, Fm_step_start,
            XH2Om, XCm, XNm, XSm, X_graphite_m,
            Xfe_H_m, Xfe_C_m, Xfe_N_m, Xfe_S_m,
            Xmin_troilite_m, Xmin_schreibersite_m, Xmin_cohenite_m,
            Xmin_graphite_m, Xmin_nitride_m, Xmin_metal_matrix_m,
            t_accreted, hcnspo_props, redox_props,
        ) = _extract_optional_marker_arrays(markers, cfg, marknum, cfg.magma_transport.active)

        define_markers!(
            markers;
            coords=coords,
            xcenter_val=xcenter_val,
            ycenter_val=ycenter_val,
            rplanet_val=rplanet_val,
            rcrust_val=rcrust_val,
            XWsolidm_init_val=cfg.materials.XWsolidm_init,
            phim0_val=cfg.thermodynamics.phim0,
            Xfe_bulk_val=cfg.coreformation.Xfe_bulk,
            T_eutectic_val=cfg.coreformation.T_eutectic,
            dT_metal_val=cfg.coreformation.dT_metal,
            tkm0_val=cfg.materials.tkm0,
            gggsolidm_val=cfg.materials.gggsolidm,
            frictsolidm_val=cfg.materials.frictsolidm,
            cohessolidm_val=cfg.materials.cohessolidm,
            tenssolidm_val=cfg.materials.tenssolidm,
            rhosolidm_val=cfg.materials.rhosolidm,
            rhofluidm_val=cfg.materials.rhofluidm,
            etasolidm_val=cfg.materials.etasolidm,
            rhocpsolidm_val=cfg.materials.rhocpsolidm,
            alphasolidm_val=cfg.materials.alphasolidm,
            alphafluidm_val=cfg.materials.alphafluidm,
            ksolidm_val=cfg.materials.ksolidm,
            start_hrsolidm_val=hrsolidm,
            phimin_val=cfg.poroelasticity.phimin,
            rng=rng,
        )

        core.XWsolidm .= core.XWsolidm0
        core.phinewm .= core.phim
        _init_fresh_hcnspo_markers!(hcnspo_props, core.tm, marknum, cfg)

        resize!(core.w3d_m, marknum)
        for m in 1:marknum
            core.w3d_m[m] = marker_out_of_plane_length(core.xm[m], core.ym[m], xcenter_val, ycenter_val)
        end

        M_atm_species = if cfg.escape.multi_species || cfg.atmosphere.active || cfg.escape.active
            Dict{Symbol,Float64}(sp => 0.0 for sp in cfg.escape.species_list)
        else
            nothing
        end
        M_escaped_species = if cfg.escape.multi_species || cfg.atmosphere.active || cfg.escape.active
            Dict{Symbol,Float64}(sp => 0.0 for sp in cfg.escape.species_list)
        else
            nothing
        end

        accumulators = SimulationAccumulators(
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            Float64(cfg.disk.p_amb_disk),
            Float64(rplanet_val),
            0.0,
            0,
            0.0,
            Float64(M_planet_val),
            Float64(xcenter_val),
            Float64(ycenter_val),
            0.0,
            M_atm_species,
            M_escaped_species,
            nothing,
            nothing,
        )

        state = SimulationState(
            grids,
            markers,
            accumulators,
            transfer_log,
            atm_state,
            rng,
            timer,
            0,
            dt,
            timesum,
        )

        if cfg.output.mode != :telemetry && cfg.output.savematstep > 0
            save_state(output_path, state, coords, cfg)
        end
    end

    return state, coords_actual, start_step_val, is_restart
end

"""
Allocate workspace memory containers for simulation solvers.

$(SIGNATURES)
"""
function setup_simulation_workspaces(
    cfg::SimulationConfig, coords::GridCoordinates, marknum::Int
)
    darcy_elim_val = cfg.solver.darcy_elimination
    dof_per_node_val = darcy_elim_val ? 4 : 6
    R, S = setup_hydromechanical_lse(coords; dof_per_node=dof_per_node_val)
    hydromech = HydromechanicalLSEWorkspace(coords; dof_per_node=dof_per_node_val)

    RT, ST = setup_thermal_lse(coords)
    thermal = ThermalLSEWorkspace(coords)

    RP, SP = setup_gravitational_lse(coords)
    F_grav = if cfg.geometry.gravity_mode === :poisson2d
        LP = assemble_gravitational_lse!(zeros(coords.Ny1, coords.Nx1), RP; coords=coords)
        lu(LP.cscmatrix)
    else
        nothing
    end

    coreformation_active_val =
        cfg.coreformation.percolation_active ||
        cfg.coreformation.settling_active ||
        cfg.metal_partition.active
    metal_segregation = if coreformation_active_val
        MetalSegregationWorkspace(
            coords.Ny, coords.Nx; track_volatiles=cfg.metal_partition.active
        )
    else
        nothing
    end

    magma_segregation = if cfg.magma_transport.active
        MagmaSegregationWorkspace(coords.Ny, coords.Nx)
    else
        nothing
    end

    use_tiled_p2m = cfg.solver.p2m_mode == :tiled
    p2m = use_tiled_p2m ? P2MTiledWorkspace(coords, marknum, cfg.solver.tile_size) : nothing

    thread_buffers = if (!use_tiled_p2m && Threads.nthreads() > 1)
        allocate_thread_interpolation_buffers(16, coords)
    else
        nothing
    end

    interp_arrays = setup_interpolated_properties(coords)
    mdis, mnum = setup_marker_geometry_helpers(coords)
    YERRNOD = zeros(Float64, cfg.solver.max_plastic_iterations)
    fractured_cells = zeros(Bool, coords.Ny, coords.Nx)
    fractured_cells_prev = zeros(Bool, coords.Ny, coords.Nx)

    has_metal = coreformation_active_val
    has_metal_part = cfg.metal_partition.active
    magma_on = cfg.magma_transport.active
    Xfe_bulk_step_start = has_metal ? zeros(Float64, marknum) : nothing
    Xfem_step_start = has_metal ? zeros(Float64, marknum) : nothing
    Fm_step_start = magma_on ? zeros(Float64, marknum) : nothing
    F_extract_m_step_start = magma_on ? zeros(Float64, marknum) : nothing
    Xfe_H_m_step_start = has_metal_part ? zeros(Float64, marknum) : nothing
    Xfe_C_m_step_start = has_metal_part ? zeros(Float64, marknum) : nothing
    Xfe_N_m_step_start = has_metal_part ? zeros(Float64, marknum) : nothing
    Xfe_S_m_step_start = has_metal_part ? zeros(Float64, marknum) : nothing

    return SimulationWorkspaces(
        hydromech,
        thermal,
        metal_segregation,
        magma_segregation,
        p2m,
        thread_buffers,
        interp_arrays,
        nothing,
        nothing,
        R,
        S,
        RT,
        ST,
        RP,
        SP,
        F_grav,
        mdis,
        mnum,
        YERRNOD,
        fractured_cells,
        fractured_cells_prev,
        Xfe_bulk_step_start,
        Xfem_step_start,
        Fm_step_start,
        F_extract_m_step_start,
        Xfe_H_m_step_start,
        Xfe_C_m_step_start,
        Xfe_N_m_step_start,
        Xfe_S_m_step_start,
    )
end

"""
Initialize simulation state, coordinates, workspaces, and telemetry stream.

$(SIGNATURES)
"""
function init_simulation(
    cfg::SimulationConfig=default_config();
    output_path::String=cfg.output.output_dir,
    restart_from::AbstractString=cfg.output.restart_from,
    force_restart_config::Bool=false,
)
    coords, formatted_output_path = setup_simulation_grid_and_coords(cfg; output_path=output_path)
    state, coords_actual, start_step_val, is_restart = setup_simulation_markers_and_state(
        cfg, coords;
        output_path=formatted_output_path,
        restart_from=restart_from,
        force_restart_config=force_restart_config,
    )
    ws = setup_simulation_workspaces(cfg, coords_actual, length(state.markers))
    telemetry_io = if (cfg.output.mode in (:telemetry, :both))
        init_telemetry(formatted_output_path, cfg.output.telemetry_file; append=is_restart)
    else
        nothing
    end
    return state, coords_actual, ws, telemetry_io, start_step_val
end

"""
Reconstruct marker arrays from a checkpoint tuple and core group.

# Parameters
- `core`: Core marker properties group.
- `optional_tuple`: NamedTuple containing optional marker field arrays.
- `marknum`: Number of markers.

# Returns
- Reconstructed `MarkerArrays` struct.
"""
function _reconstruct_checkpoint_marker_arrays(
    core::CoreGroup, optional_tuple::NamedTuple, marknum::Int
)
    group_pairs = Pair{Symbol,Any}[]
    o = optional_tuple
    if o.Xfem !== nothing
        push!(
            group_pairs,
            :metal => MetalGroup(
                o.Xfem,
                o.Xfem0,
                o.Xfe_bulk,
                o.Xfe_H_m !== nothing ? o.Xfe_H_m : zeros(Float64, marknum),
                o.Xfe_C_m !== nothing ? o.Xfe_C_m : zeros(Float64, marknum),
                o.Xfe_N_m !== nothing ? o.Xfe_N_m : zeros(Float64, marknum),
                o.Xfe_S_m !== nothing ? o.Xfe_S_m : zeros(Float64, marknum),
            ),
        )
    end
    if o.XH2Om !== nothing || o.F_extract_m !== nothing
        push!(
            group_pairs,
            :volatiles => VolatilesGroup(
                o.XH2Om !== nothing ? o.XH2Om : zeros(Float64, marknum),
                o.XCm !== nothing ? o.XCm : zeros(Float64, marknum),
                o.XNm !== nothing ? o.XNm : zeros(Float64, marknum),
                o.XSm !== nothing ? o.XSm : zeros(Float64, marknum),
                o.X_graphite_m !== nothing ? o.X_graphite_m : zeros(Float64, marknum),
                o.F_extract_m !== nothing ? o.F_extract_m : zeros(Float64, marknum),
            ),
        )
    end
    if o.redox_props !== nothing
        rp = o.redox_props
        push!(
            group_pairs,
            :redox => RedoxGroup(
                rp.nFe0_m,
                rp.nFe2_m,
                rp.nFe3_m,
                rp.deltaIW_m,
                rp.nC_graphite_m,
                rp.nCO_m,
                rp.nCO2_m,
                rp.nCH4_m,
            ),
        )
    end
    if o.hcnspo_props !== nothing
        hp = o.hcnspo_props
        push!(
            group_pairs,
            :hcnspo => HcnspoGroup(
                hp.X_ice_H2O_m,
                hp.X_ice_NH3_m,
                hp.X_ice_CO2_m,
                hp.X_ice_CO_m,
                hp.X_ice_CH4_m,
                hp.X_ice_N2_m,
                hp.X_ice_H2S_m,
                hp.X_ice_PH3_m,
                hp.X_refr_C_m,
                hp.X_refr_S_m,
                hp.X_refr_N_m,
                hp.X_refr_P_m,
                hp.X_refr_H_m,
            ),
        )
    end
    if o.Xmin_troilite_m !== nothing
        push!(
            group_pairs,
            :phase => PhaseGroup(
                o.Xmin_troilite_m,
                o.Xmin_schreibersite_m,
                o.Xmin_cohenite_m,
                o.Xmin_graphite_m,
                o.Xmin_nitride_m,
                o.Xmin_metal_matrix_m,
            ),
        )
    end
    if o.t_accreted !== nothing
        push!(group_pairs, :accretion => AccretionGroup(o.t_accreted))
    end
    return MarkerArrays(core, NamedTuple(group_pairs))
end


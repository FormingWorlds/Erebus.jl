"""
Compute surface-mean oxygen fugacity (ΔIW) from near-surface silicate markers.

$(SIGNATURES)
"""
function compute_surface_mean_delta_iw(
    redox_props,
    tm,
    xm,
    ym,
    marknum::Int,
    xcenter::Float64,
    ycenter::Float64,
    rplanet_val::Float64,
    fallback::Float64,
)
    if redox_props === nothing || redox_props.deltaIW_m === nothing
        return fallback
    end
    surf_count = 0
    surf_diw = 0.0
    r_cut_sq = (0.8 * rplanet_val)^2
    r_max_sq = rplanet_val^2
    for m in 1:marknum
        if tm[m] < 3
            r_sq = (xm[m] - xcenter)^2 + (ym[m] - ycenter)^2
            if r_sq >= r_cut_sq && r_sq <= r_max_sq
                surf_diw += redox_props.deltaIW_m[m]
                surf_count += 1
            end
        end
    end
    if surf_count > 0
        return surf_diw / surf_count
    else
        delta_iw_arr = redox_props.deltaIW_m
        return isempty(delta_iw_arr) ? fallback : sum(delta_iw_arr) / length(delta_iw_arr)
    end
end

"""
Compute multi-species surface volatile venting rates [kg/s].

$(SIGNATURES)

# Arguments
- `cfg`: Simulation configuration.
- `delta_m_vent_3d`: 3D-equivalent pore fluid water mass vented from surface Darcy sink [kg].
- `vented_vols`: Drained mobile mineral volatiles NamedTuple or nothing.
- `L_3D_equiv`: 3D geometric cross-section factor [m].
- `dt`: Current timestep [s].
- `P_surf`: Surface boundary pressure [Pa].
- `T_surf`: Surface boundary temperature [K].
- `redox_props`: Redox properties struct or nothing.
- `marknum`: Number of active markers.
- `tm`: Marker material type array.
- `xm`: Marker x coordinates [m].
- `ym`: Marker y coordinates [m].
- `rplanet_val`: Planetesimal radius [m].
- `xcenter_val`: Planetesimal center x coordinate [m] (default: 0.0).
- `ycenter_val`: Planetesimal center y coordinate [m] (default: 0.0).

# Returns
- `Dict{Symbol,Float64}` mapping volatile species to surface venting rates [kg/s].
"""
function compute_surface_venting_rates(
    cfg::SimulationConfig,
    delta_m_vent_3d::Float64,
    vented_vols,
    L_3D_equiv::Float64,
    dt::Float64,
    P_surf::Float64,
    T_surf::Float64,
    redox_props,
    marknum::Int,
    tm,
    xm,
    ym,
    rplanet_val::Float64,
    xcenter_val::Float64=0.0,
    ycenter_val::Float64=0.0,
)
    vent_rates = Dict{Symbol,Float64}(sp => 0.0 for sp in cfg.escape.species_list)
    (
        dt <= 0.0 || (
            !cfg.venting.active && (
                vented_vols === nothing ||
                !cfg.retention.active ||
                !cfg.retention.venting_drainage_active
            )
        )
    ) && return vent_rates

    drain_on = (
        vented_vols !== nothing &&
        cfg.retention.active &&
        cfg.retention.venting_drainage_active
    )
    m_pore_H2O = cfg.venting.active ? delta_m_vent_3d : 0.0
    m_mineral_H2O = drain_on ? vented_vols.M_vent_H2O * L_3D_equiv : 0.0
    m_H2O_step = m_pore_H2O + m_mineral_H2O

    m_C_step = drain_on ? vented_vols.M_vent_C * L_3D_equiv : 0.0
    m_N_step = drain_on ? vented_vols.M_vent_N * L_3D_equiv : 0.0
    m_S_step = drain_on ? vented_vols.M_vent_S * L_3D_equiv : 0.0

    (m_H2O_step + m_C_step + m_N_step + m_S_step) <= 0.0 && return vent_rates

    if cfg.volatiles.speciation_active
        fO2_delta_IW_vent = compute_surface_mean_delta_iw(
            redox_props,
            tm,
            xm,
            ym,
            marknum,
            xcenter_val,
            ycenter_val,
            rplanet_val,
            cfg.volatiles.fO2_delta_IW,
        )
        fO2_delta_IW_vent = clamp(fO2_delta_IW_vent, -50.0, 50.0)
        P_surf_eval = max(P_surf, 1.0)
        T_surf_eval = max(T_surf, 273.15)

        spec_dict = speciate_vented_volatiles(
            m_H2O_step,
            m_C_step,
            m_N_step,
            m_S_step,
            P_surf_eval,
            T_surf_eval,
            fO2_delta_IW_vent;
            graphite_saturation=cfg.volatiles.graphite_saturation,
        )
        for (sp, m_sp) in spec_dict
            vent_rates[sp] = get(vent_rates, sp, 0.0) + m_sp / dt
        end
    else
        vent_sp = cfg.venting.species
        vent_rates[vent_sp] = get(vent_rates, vent_sp, 0.0) + m_pore_H2O / dt
        vent_rates[:H2O] = get(vent_rates, :H2O, 0.0) + m_mineral_H2O / dt
        # Stoichiometric conversion: elemental C to CO2 (44.0095 / 12.011)
        vent_rates[:CO2] = get(vent_rates, :CO2, 0.0) + (m_C_step * (44.0095 / 12.011)) / dt
        vent_rates[:N2] = get(vent_rates, :N2, 0.0) + m_N_step / dt
        # Stoichiometric conversion: elemental S to H2S (34.08 / 32.06)
        vent_rates[:H2S] = get(vent_rates, :H2S, 0.0) + (m_S_step * (34.08 / 32.06)) / dt
    end

    return vent_rates
end

"""
Main simulation loop: run calculations with timestepping.

$(SIGNATURES)

# Details

    - output_path: Absolute path where to save simulation output files
    - restart_from: Optional path to checkpoint JLD2 file to resume from

# Returns
    
    - nothing
"""
function simulation_loop(
    cfg::SimulationConfig=default_config();
    output_path::String=cfg.output.output_dir,
    restart_from::AbstractString=cfg.output.restart_from,
)
    if cfg.mpi.enable
        error(
            "Distributed multi-node orchestration for simulation_loop is scheduled for Milestone 4/5.",
        )
    end
    output_path = endswith(output_path, "/") ? output_path : output_path * "/"
    isdir(output_path) || mkpath(output_path)
    coords = GridCoordinates(cfg.grid)

    # -------------------------------------------------------------------------
    # set up dynamic simulation parameters from given static parameters"
    # -------------------------------------------------------------------------
    timestep, dt, timesum, marknum, hrsolidm, hrfluidm, YERRNOD = setup_dynamic_simulation_parameters(
        cfg; coords=coords
    )

    # Extract dynamic simulation control parameters from cfg
    n_steps_val = cfg.time.n_steps
    start_step_val = cfg.time.start_step
    savematstep_val = cfg.output.savematstep
    titermax_val = cfg.solver.titermax
    use_pardiso_val = cfg.solver.use_pardiso
    etaphikoef_val = cfg.solver.etaphikoef
    betasolid_val = cfg.poroelasticity.betasolid
    betafluid_val = cfg.poroelasticity.betafluid
    phimin_val = cfg.poroelasticity.phimin
    phimax_val = cfg.poroelasticity.phimax
    hydrofracture_val = cfg.poroelasticity.hydrofracture
    kappa_frac_val = cfg.poroelasticity.kappa_frac
    gamma_frac_val = cfg.poroelasticity.gamma_frac
    k_frac_max_val = cfg.poroelasticity.k_frac_max
    dt_longest_val = cfg.time.dt_longest * cfg.time.yearlength
    endtime_val = cfg.time.endtime * cfg.time.yearlength
    dtcoefup_val = cfg.time.dtcoefup
    hr_al_val = cfg.thermodynamics.hr_al
    hr_fe_val = cfg.thermodynamics.hr_fe
    ratio_al_val = cfg.thermodynamics.ratio_al
    ratio_fe_val = cfg.thermodynamics.ratio_fe
    E_al_val = cfg.thermodynamics.E_al
    f_al_val = cfg.thermodynamics.f_al
    tau_al_val = cfg.thermodynamics.t_half_al / log(2.0)
    E_fe_val = cfg.thermodynamics.E_fe
    f_fe_val = cfg.thermodynamics.f_fe
    tau_fe_val = cfg.thermodynamics.t_half_fe / log(2.0)
    thermal_buoyancy_val = cfg.thermodynamics.thermal_buoyancy
    tmfluidphase_val = cfg.thermodynamics.tmfluidphase
    alphafluid_val = cfg.materials.alphafluidm
    fluid_viscosity_mode_val = cfg.thermodynamics.fluid_viscosity_mode
    fluid_viscosity_Ea_val = cfg.thermodynamics.fluid_viscosity_Ea
    fluid_viscosity_T0_val = cfg.thermodynamics.fluid_viscosity_T0
    fluid_viscosity_eta0_val = cfg.thermodynamics.fluid_viscosity_eta0
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
    M_accreted_total = 0.0
    xcenter_val = cfg.geometry.xcenter
    ycenter_val = cfg.geometry.ycenter
    psurface_val = cfg.geometry.psurface
    spherical_metric_val = cfg.geometry.spherical_metric
    metric_reg_val = cfg.geometry.metric_regularization_cells
    surface_radiation_val = cfg.thermodynamics.surface_radiation
    emissivity_val = cfg.thermodynamics.emissivity
    sigma_sb_val = cfg.thermodynamics.sigma_sb
    phim0_val = cfg.thermodynamics.phim0
    kfluidm_val = cfg.materials.kfluidm
    disk_enabled_val = cfg.disk.enabled
    reaction_active_val = cfg.reaction.active
    delta_H_val = cfg.reaction.delta_H
    delta_S_val = cfg.reaction.delta_S
    pfcoeff_val = cfg.reaction.pfcoeff
    pferrmax_val = cfg.reaction.pferrmax
    dtr_hyd_val = cfg.reaction.dtreaction_hydration
    dtr_deh_val = cfg.reaction.dtreaction_dehydration
    melting_active_val = cfg.melting.active
    T_solidus_val = cfg.melting.T_solidus
    T_liquidus_val = cfg.melting.T_liquidus
    L_melt_val = cfg.melting.L_melt
    rho_melt_val = cfg.melting.rho_melt
    alpha_eta_val = cfg.melting.alpha_eta
    phi_crit_val = cfg.melting.phi_crit
    eta_melt_val = cfg.melting.eta_melt
    dpdt_clapeyron_val = cfg.melting.dpdt_clapeyron
    soft_turbulence_val = cfg.melting.soft_turbulence
    eta_fluid_silicate_val = cfg.melting.eta_fluid_silicate
    F_turb_start_val = cfg.melting.F_turb_start
    F_turb_end_val = cfg.melting.F_turb_end
    turb_exponent_val = cfg.melting.turb_exponent
    dT_turb_min_val = cfg.melting.dT_turb_min
    T_surface_ref_val = cfg.melting.T_surface_ref
    k_turb_cutoff_val = cfg.melting.k_turb_cutoff
    k_turb_floor_val = cfg.melting.k_turb_floor
    coreformation_active_val =
        cfg.coreformation.percolation_active ||
        cfg.coreformation.settling_active ||
        cfg.metal_partition.active
    Xfe_bulk_val = cfg.coreformation.Xfe_bulk
    T_eutectic_val = cfg.coreformation.T_eutectic
    dT_metal_val = cfg.coreformation.dT_metal
    rho_metal_val = cfg.coreformation.rho_metal
    rho_metal_solid_val = cfg.coreformation.rho_metal_solid
    sulfur_fraction_val = cfg.coreformation.sulfur_fraction
    metal_density_mode_val = cfg.coreformation.metal_density_mode
    L_metal_val = cfg.coreformation.L_metal
    k_metal_val = cfg.coreformation.k_metal
    rhocp_metal_val = cfg.coreformation.rhocp_metal
    max_v_seg_prev = 0.0

    nthreads = Threads.nthreads()

    use_threading = nthreads > 1
    num_buffers = nthreads
    use_tiled_p2m = cfg.solver.p2m_mode == :tiled
    thread_buffers = if (!use_tiled_p2m && use_threading)
        allocate_thread_interpolation_buffers(num_buffers, coords)
    else
        nothing
    end
    p2m_workspace =
        use_tiled_p2m ? P2MTiledWorkspace(coords, marknum, cfg.solver.tile_size) : nothing

    @info "Simulation layout" coords.Nx coords.Ny coords.xsize coords.dx coords.dy coords.ysize rplanet_val rcrust_val marknum nthreads
    @info(
        "Parameters",
        random_markers,
        marker_property_mode,
        hr_al,
        hr_fe,
        reaction_active = reaction_active_val,
        reaction_rate_coeff_mode,
        log_completion_rate,
        t_half_al,
        ratio_al,
        E_al,
        f_al,
        t_half_fe,
        ratio_fe,
        E_fe,
        f_fe,
        rhosolidm,
        rhofluidm,
        etasolidm,
        etasolidmm,
        etafluidm,
        etafluidmm,
        rhocpsolidm,
        rhocpfluidm,
        alphasolidm,
        alphafluidm,
        ksolidm,
        kfluidm,
        gggsolidm,
        frictsolidm,
        cohessolidm,
        tenssolidm,
        kphim0,
        tkm0,
        etaphikoef,
        αη,
        tmsolidphase,
        tmfluidphase,
        phim0=phim0_val,
        phimin,
        phimax,
        ΔHWD = delta_H_val,
        ΔSWD = delta_S_val,
        ΔVWD,
        dtreaction_hydration = dtr_hyd_val,
        dtreaction_dehydration = dtr_deh_val,
        pfcoeff = pfcoeff_val,
        pferrmax = pferrmax_val,
        start_time,
        start_step,
        endtime,
        dsubgrids,
        dsubgridt,
        dt_longest,
        dphimax,
        dxymax,
        vpratio,
        seed
    )
    @info "Solver" use_pardiso BLAS.get_config() BLAS.get_num_threads()

    # -------------------------------------------------------------------------
    # set up staggered grid"
    # -------------------------------------------------------------------------
    (ETA, ETA0, GGG, EXY, SXY, SXY0, wyx, COH, TEN, FRI, YNY, RHOX, RHOFX, KX, PHIX, vx, vxf, RX, qxD, gx, RHOY, RHOFY, KY, PHIY, vy, vyf, RY, qyD, gy, RHO, RHOCP, ALPHA, ALPHAF, HR, HA, HS, ETAP, GGGP, EXX, SXX, SXX0, tk1, tk2, DT, DT0, vxp, vyp, vxpf, vypf, pr, pf, ps, pr0, pf0, ps0, ETAPHI, BETAPHI, PHI, APHI, FI, DMP, DHP, XWS) = setup_staggered_grid_properties(
        coords
    )
    (ETA5, ETA00, YNY5, YNY00, YNY_inv_ETA, DSXY, DSY, EII, SII, DSXX, tk0) = setup_staggered_grid_properties_helpers(
        coords
    )
    Q_metric = spherical_metric_val ? zeros(Float64, coords.Ny1, coords.Nx1) : nothing
    DQPF = zeros(Float64, coords.Ny1, coords.Nx1)
    DQPFSUM = zeros(Float64, coords.Ny1, coords.Nx1)

    # -------------------------------------------------------------------------
    # set up markers and state (from checkpoint or fresh definition)
    # -------------------------------------------------------------------------
    mdis, mnum = setup_marker_geometry_helpers(coords)
    M_vent_total = 0.0
    M_vent_H2O_total = 0.0
    M_vent_C_total = 0.0
    M_vent_N_total = 0.0
    M_vent_S_total = 0.0
    M_atm_total = 0.0
    M_escaped_total = 0.0
    S_vent_grid = zeros(Float64, coords.Ny1, coords.Nx1)
    Q_lat_grid = zeros(Float64, coords.Ny1, coords.Nx1)
    Q_seg_grid = zeros(Float64, coords.Ny1, coords.Nx1)
    Xfem = nothing
    Xfem0 = nothing
    Xfe_bulk = nothing
    Xfe_bulk_step_start = nothing
    Xfem_step_start = nothing
    XH2Om = nothing
    XCm = nothing
    XNm = nothing
    XSm = nothing
    Xfe_H_m = nothing
    Xfe_C_m = nothing
    Xfe_N_m = nothing
    Xfe_S_m = nothing
    Xfe_H_m_step_start = nothing
    Xfe_C_m_step_start = nothing
    Xfe_N_m_step_start = nothing
    Xfe_S_m_step_start = nothing
    core_budgets = nothing
    rcore_val = 0.0
    Xmin_troilite_m = nothing
    Xmin_schreibersite_m = nothing
    Xmin_cohenite_m = nothing
    Xmin_graphite_m = nothing
    Xmin_nitride_m = nothing
    Xmin_metal_matrix_m = nothing
    regional_mineral_modes = nothing
    magma_active_val = cfg.magma_transport.active
    F_extract_m = nothing
    F_extract_m_step_start = nothing
    Fm_step_start = nothing
    M_atm_species = if cfg.escape.multi_species || cfg.atmosphere.active
        Dict{Symbol,Float64}(sp => 0.0 for sp in cfg.escape.species_list)
    else
        nothing
    end
    M_escaped_species = if cfg.escape.multi_species || cfg.atmosphere.active
        Dict{Symbol,Float64}(sp => 0.0 for sp in cfg.escape.species_list)
    else
        nothing
    end
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
    telescope_level = 0
    is_restart = !isempty(restart_from)
    telemetry_io = nothing
    if is_restart
        ckpt = load_state(restart_from)
        if haskey(ckpt, "telescope_level")
            telescope_level = Int(ckpt["telescope_level"])
        end
        if haskey(ckpt, "M_vent_total")
            M_vent_total = Float64(ckpt["M_vent_total"])
        end
        if haskey(ckpt, "M_vent_H2O_total")
            M_vent_H2O_total = Float64(ckpt["M_vent_H2O_total"])
        end
        if haskey(ckpt, "M_vent_C_total")
            M_vent_C_total = Float64(ckpt["M_vent_C_total"])
        end
        if haskey(ckpt, "M_vent_N_total")
            M_vent_N_total = Float64(ckpt["M_vent_N_total"])
        end
        if haskey(ckpt, "M_vent_S_total")
            M_vent_S_total = Float64(ckpt["M_vent_S_total"])
        end
        if haskey(ckpt, "M_atm_total")
            M_atm_total = Float64(ckpt["M_atm_total"])
        end
        if haskey(ckpt, "M_escaped_total")
            M_escaped_total = Float64(ckpt["M_escaped_total"])
        end
        if cfg.escape.multi_species || cfg.atmosphere.active
            if haskey(ckpt, "M_atm_species")
                raw_atm = ckpt["M_atm_species"]
                M_atm_species = Dict{Symbol,Float64}(
                    Symbol(k) => Float64(v) for (k, v) in pairs(raw_atm)
                )
            end
            if haskey(ckpt, "M_escaped_species")
                raw_esc = ckpt["M_escaped_species"]
                M_escaped_species = Dict{Symbol,Float64}(
                    Symbol(k) => Float64(v) for (k, v) in pairs(raw_esc)
                )
            end
        end
        if cfg.atmosphere.active && atm_state !== nothing
            if haskey(ckpt, "atm_M_atm")
                for (k, v) in pairs(ckpt["atm_M_atm"])
                    atm_state.M_atm[Symbol(k)] = Float64(v)
                end
            end
            if haskey(ckpt, "atm_M_escaped")
                for (k, v) in pairs(ckpt["atm_M_escaped"])
                    atm_state.M_escaped[Symbol(k)] = Float64(v)
                end
            end
            if haskey(ckpt, "atm_P_surf")
                atm_state.P_surf = Float64(ckpt["atm_P_surf"])
            end
            if haskey(ckpt, "atm_T_surf_eq")
                atm_state.T_surf_eq = Float64(ckpt["atm_T_surf_eq"])
            end
            if haskey(ckpt, "atm_tau_LW")
                atm_state.tau_LW = Float64(ckpt["atm_tau_LW"])
            end
            if haskey(ckpt, "atm_M_env_bound")
                atm_state.M_env_bound = Float64(ckpt["atm_M_env_bound"])
            end
            if haskey(ckpt, "atm_F_net_rad")
                atm_state.F_net_rad = Float64(ckpt["atm_F_net_rad"])
            end
            if haskey(ckpt, "atm_h_rad_eff")
                atm_state.h_rad_eff = Float64(ckpt["atm_h_rad_eff"])
            end
        end
        if cfg.telescoping.active &&
            haskey(ckpt, "Nx") &&
            haskey(ckpt, "Ny") &&
            (ckpt["Nx"] != coords.Nx || ckpt["Ny"] != coords.Ny)
            coords = GridCoordinates(
                ckpt["Nx"],
                ckpt["Ny"];
                xsize=Float64(ckpt["xsize"]),
                ysize=Float64(ckpt["ysize"]),
                Nxmc=coords.Nxmc,
                Nymc=coords.Nymc,
            )
            xcenter_val = coords.xcenter
            ycenter_val = coords.ycenter
            (ETA, ETA0, GGG, EXY, SXY, SXY0, wyx, COH, TEN, FRI, YNY, RHOX, RHOFX, KX, PHIX, vx, vxf, RX, qxD, gx, RHOY, RHOFY, KY, PHIY, vy, vyf, RY, qyD, gy, RHO, RHOCP, ALPHA, ALPHAF, HR, HA, HS, ETAP, GGGP, EXX, SXX, SXX0, tk1, tk2, DT, DT0, vxp, vyp, vxpf, vypf, pr, pf, ps, pr0, pf0, ps0, ETAPHI, BETAPHI, PHI, APHI, FI, DMP, DHP, XWS) = setup_staggered_grid_properties(
                coords
            )
            (ETA5, ETA00, YNY5, YNY00, YNY_inv_ETA, DSXY, DSY, EII, SII, DSXX, tk0) = setup_staggered_grid_properties_helpers(
                coords
            )
            Q_metric =
                spherical_metric_val ? zeros(Float64, coords.Ny1, coords.Nx1) : nothing
            DQPF = zeros(Float64, coords.Ny1, coords.Nx1)
            DQPFSUM = zeros(Float64, coords.Ny1, coords.Nx1)
            mdis, mnum = setup_marker_geometry_helpers(coords)
            S_vent_grid = zeros(Float64, coords.Ny1, coords.Nx1)
            Q_lat_grid = zeros(Float64, coords.Ny1, coords.Nx1)
            Q_seg_grid = zeros(Float64, coords.Ny1, coords.Nx1)
            if !use_tiled_p2m && use_threading
                thread_buffers = allocate_thread_interpolation_buffers(num_buffers, coords)
            end
        else
            if haskey(ckpt, "Nx") && haskey(ckpt, "Ny")
                (ckpt["Nx"] == coords.Nx && ckpt["Ny"] == coords.Ny) || throw(
                    DimensionMismatch(
                        "Checkpoint grid size ($(ckpt["Nx"])x$(ckpt["Ny"])) does not match current simulation grid size ($(coords.Nx)x$(coords.Ny))",
                    ),
                )
            end
            if haskey(ckpt, "xsize") && haskey(ckpt, "ysize")
                (ckpt["xsize"] == coords.xsize && ckpt["ysize"] == coords.ysize) || throw(
                    DimensionMismatch(
                        "Checkpoint domain size ($(ckpt["xsize"])x$(ckpt["ysize"])) does not match current simulation domain size ($(coords.xsize)x$(coords.ysize))",
                    ),
                )
            end
        end
        start_step_val = ckpt["timestep"] + 1
        dt = ckpt["dt"]
        timesum = ckpt["timesum"]
        marknum = ckpt["marknum"]
        n_steps_val = cfg.time.n_steps
        if n_steps_val < start_step_val
            @warn "Restart checkpoint timestep ($(ckpt["timestep"])) >= target n_steps ($n_steps_val). No timesteps will be executed."
        end

        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0, Fm) = setup_marker_properties(
            marknum, coords
        )
        (rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm, tkm_rhocptotalm, etafluidcur_inv_kphim, inv_gggtotalm, fricttotalm, cohestotalm, tenstotalm, rhofluidcur, alphasolidcur, alphafluidcur) = setup_marker_properties_helpers(
            marknum
        )

        # Restore staggered grid arrays
        ETA .= ckpt["ETA"]
        ETA0 .= ckpt["ETA0"]
        GGG .= ckpt["GGG"]
        EXY .= ckpt["EXY"]
        SXY .= ckpt["SXY"]
        SXY0 .= ckpt["SXY0"]
        wyx .= ckpt["wyx"]
        COH .= ckpt["COH"]
        TEN .= ckpt["TEN"]
        FRI .= ckpt["FRI"]
        YNY .= ckpt["YNY"]
        RHOX .= ckpt["RHOX"]
        RHOFX .= ckpt["RHOFX"]
        KX .= ckpt["KX"]
        PHIX .= ckpt["PHIX"]
        vx .= ckpt["vx"]
        vxf .= ckpt["vxf"]
        RX .= ckpt["RX"]
        qxD .= ckpt["qxD"]
        gx .= ckpt["gx"]
        RHOY .= ckpt["RHOY"]
        RHOFY .= ckpt["RHOFY"]
        KY .= ckpt["KY"]
        PHIY .= ckpt["PHIY"]
        vy .= ckpt["vy"]
        vyf .= ckpt["vyf"]
        RY .= ckpt["RY"]
        qyD .= ckpt["qyD"]
        gy .= ckpt["gy"]
        RHO .= ckpt["RHO"]
        RHOCP .= ckpt["RHOCP"]
        ALPHA .= ckpt["ALPHA"]
        ALPHAF .= ckpt["ALPHAF"]
        HR .= ckpt["HR"]
        HA .= ckpt["HA"]
        HS .= ckpt["HS"]
        ETAP .= ckpt["ETAP"]
        GGGP .= ckpt["GGGP"]
        EXX .= ckpt["EXX"]
        SXX .= ckpt["SXX"]
        SXX0 .= ckpt["SXX0"]
        tk1 .= ckpt["tk1"]
        tk2 .= ckpt["tk2"]
        if haskey(ckpt, "DT0")
            DT0 .= ckpt["DT0"]
        end
        pr .= ckpt["pr"]
        pf .= ckpt["pf"]
        ps .= ckpt["ps"]
        pr0 .= ckpt["pr0"]
        pf0 .= ckpt["pf0"]
        ps0 .= ckpt["ps0"]
        ETAPHI .= ckpt["ETAPHI"]
        BETAPHI .= ckpt["BETAPHI"]
        PHI .= ckpt["PHI"]
        APHI .= ckpt["APHI"]
        FI .= ckpt["FI"]
        DMP .= ckpt["DMP"]
        DHP .= ckpt["DHP"]
        XWS .= ckpt["XWS"]

        # Restore marker properties
        xm .= ckpt["xm"]
        ym .= ckpt["ym"]
        tm .= ckpt["tm"]
        tkm .= ckpt["tkm"]
        sxxm .= ckpt["sxxm"]
        sxym .= ckpt["sxym"]
        etavpm .= ckpt["etavpm"]
        phim .= ckpt["phim"]
        phinewm .= phim
        rhototalm .= ckpt["rhototalm"]
        rhocptotalm .= ckpt["rhocptotalm"]
        etatotalm .= ckpt["etatotalm"]
        hrtotalm .= ckpt["hrtotalm"]
        ktotalm .= ckpt["ktotalm"]
        tkm_rhocptotalm .= ckpt["tkm_rhocptotalm"]
        etafluidcur_inv_kphim .= ckpt["etafluidcur_inv_kphim"]
        inv_gggtotalm .= ckpt["inv_gggtotalm"]
        fricttotalm .= ckpt["fricttotalm"]
        cohestotalm .= ckpt["cohestotalm"]
        tenstotalm .= ckpt["tenstotalm"]
        rhofluidcur .= ckpt["rhofluidcur"]
        alphasolidcur .= ckpt["alphasolidcur"]
        alphafluidcur .= ckpt["alphafluidcur"]
        XWsolidm0 .= ckpt["XWsolidm0"]
        XWsolidm .= XWsolidm0
        if coreformation_active_val || hr_fe_val
            Xfem, Xfem0, Xfe_bulk = setup_marker_metal_properties(marknum)
            if haskey(ckpt, "Xfe_bulk")
                Xfe_bulk .= ckpt["Xfe_bulk"]
            else
                for m in 1:marknum
                    rmark = distance(xm[m], ym[m], xcenter_val, ycenter_val)
                    if rmark < rplanet_val
                        Xfe_bulk[m] = Xfe_bulk_val
                    end
                end
            end
            if haskey(ckpt, "Xfem")
                Xfem .= ckpt["Xfem"]
            end
            if haskey(ckpt, "Xfem0")
                Xfem0 .= ckpt["Xfem0"]
            end
            Xfe_bulk_step_start = zeros(Float64, marknum)
            Xfem_step_start = zeros(Float64, marknum)
        end
        if cfg.volatiles.active
            if haskey(ckpt, "XH2Om")
                XH2Om = Vector{Float64}(ckpt["XH2Om"])
                XCm = Vector{Float64}(ckpt["XCm"])
                XNm = Vector{Float64}(ckpt["XNm"])
                XSm = Vector{Float64}(ckpt["XSm"])
            else
                (XH2Om, XCm, XNm, XSm) = setup_marker_volatile_properties(
                    marknum;
                    initial_water_wtpct=cfg.volatiles.initial_water_wtpct,
                    initial_carbon_ppm=cfg.volatiles.initial_carbon_ppm,
                    initial_nitrogen_ppm=cfg.volatiles.initial_nitrogen_ppm,
                    initial_sulfur_ppm=cfg.volatiles.initial_sulfur_ppm,
                )
            end
        end
        if cfg.metal_partition.active
            if haskey(ckpt, "Xfe_H_m")
                Xfe_H_m = Vector{Float64}(ckpt["Xfe_H_m"])
                Xfe_C_m = Vector{Float64}(ckpt["Xfe_C_m"])
                Xfe_N_m = Vector{Float64}(ckpt["Xfe_N_m"])
                Xfe_S_m = Vector{Float64}(ckpt["Xfe_S_m"])
            else
                (Xfe_H_m, Xfe_C_m, Xfe_N_m, Xfe_S_m) = setup_marker_metal_volatile_properties(
                    marknum;
                    initial_h_ppm=cfg.metal_partition.initial_metal_h_ppm,
                    initial_c_ppm=cfg.metal_partition.initial_metal_c_ppm,
                    initial_n_ppm=cfg.metal_partition.initial_metal_n_ppm,
                    initial_s_ppm=cfg.metal_partition.initial_metal_s_ppm,
                )
            end
            Xfe_H_m_step_start = zeros(Float64, marknum)
            Xfe_C_m_step_start = zeros(Float64, marknum)
            Xfe_N_m_step_start = zeros(Float64, marknum)
            Xfe_S_m_step_start = zeros(Float64, marknum)
        end
        if cfg.phase_tracking.active
            if haskey(ckpt, "Xmin_troilite_m")
                Xmin_troilite_m = Vector{Float64}(ckpt["Xmin_troilite_m"])
                Xmin_schreibersite_m = Vector{Float64}(ckpt["Xmin_schreibersite_m"])
                Xmin_cohenite_m = Vector{Float64}(ckpt["Xmin_cohenite_m"])
                Xmin_graphite_m = Vector{Float64}(ckpt["Xmin_graphite_m"])
                Xmin_nitride_m = Vector{Float64}(ckpt["Xmin_nitride_m"])
                Xmin_metal_matrix_m = Vector{Float64}(ckpt["Xmin_metal_matrix_m"])
            else
                phase_arrays = setup_marker_phase_tracking_properties(
                    marknum, cfg.phase_tracking
                )
                Xmin_troilite_m = phase_arrays.Xmin_troilite_m
                Xmin_schreibersite_m = phase_arrays.Xmin_schreibersite_m
                Xmin_cohenite_m = phase_arrays.Xmin_cohenite_m
                Xmin_graphite_m = phase_arrays.Xmin_graphite_m
                Xmin_nitride_m = phase_arrays.Xmin_nitride_m
                Xmin_metal_matrix_m = phase_arrays.Xmin_metal_matrix_m
            end
            if haskey(ckpt, "regional_mineral_modes")
                regional_mineral_modes = ckpt["regional_mineral_modes"]
            end
        end
        if cfg.accretion.active
            if haskey(ckpt, "t_accreted")
                t_accreted = Vector{Float64}(ckpt["t_accreted"])
            else
                t_accreted = setup_marker_accretion_properties(
                    marknum, cfg.accretion; initial_time=timesum
                )
            end
            if haskey(ckpt, "M_accreted_total")
                M_accreted_total = Float64(ckpt["M_accreted_total"])
            end
            if haskey(ckpt, "M_planet_val")
                M_planet_val = Float64(ckpt["M_planet_val"])
            end
            if haskey(ckpt, "rplanet")
                rplanet_val = Float64(ckpt["rplanet"])
            end
        else
            t_accreted = nothing
        end
        if haskey(ckpt, "rcore")
            rcore_val = Float64(ckpt["rcore"])
        end
        hcnspo_props = if cfg.volatile_mixture.active || cfg.refractory.active
            hp = setup_marker_hcnspo_properties(marknum, cfg.volatile_mixture, cfg.refractory)
            for k in keys(hp)
                k_str = string(k)
                if haskey(ckpt, k_str)
                    getfield(hp, k) .= ckpt[k_str]
                end
            end
            hp
        else
            nothing
        end
        redox_props = if cfg.redox.active
            rp = setup_marker_redox_properties(marknum, cfg.redox; initial_xfe_bulk=Xfe_bulk)
            for k in keys(rp)
                k_str = string(k)
                if haskey(ckpt, k_str)
                    getfield(rp, k) .= ckpt[k_str]
                end
            end
            rp
        else
            nothing
        end
        if magma_active_val
            F_extract_m = if haskey(ckpt, "F_extract_m")
                Vector{Float64}(ckpt["F_extract_m"])
            else
                zeros(Float64, marknum)
            end
            F_extract_m_step_start = zeros(Float64, marknum)
            Fm_step_start = zeros(Float64, marknum)
        end
        @info "Resumed simulation from checkpoint: $restart_from at timestep $(start_step_val-1) (running to $n_steps_val)"
    else
        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0, Fm) = setup_marker_properties(
            marknum, coords
        )
        (rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm, tkm_rhocptotalm, etafluidcur_inv_kphim, inv_gggtotalm, fricttotalm, cohestotalm, tenstotalm, rhofluidcur, alphasolidcur, alphafluidcur) = setup_marker_properties_helpers(
            marknum
        )
        if coreformation_active_val || hr_fe_val
            Xfem, Xfem0, Xfe_bulk = setup_marker_metal_properties(marknum)
            Xfe_bulk_step_start = zeros(Float64, marknum)
            Xfem_step_start = zeros(Float64, marknum)
        end
        if magma_active_val
            F_extract_m = setup_marker_magma_properties(marknum)[1]
            F_extract_m_step_start = zeros(Float64, marknum)
            Fm_step_start = zeros(Float64, marknum)
        end
        if cfg.volatiles.active
            (XH2Om, XCm, XNm, XSm) = setup_marker_volatile_properties(
                marknum;
                initial_water_wtpct=cfg.volatiles.initial_water_wtpct,
                initial_carbon_ppm=cfg.volatiles.initial_carbon_ppm,
                initial_nitrogen_ppm=cfg.volatiles.initial_nitrogen_ppm,
                initial_sulfur_ppm=cfg.volatiles.initial_sulfur_ppm,
            )
        end
        if cfg.metal_partition.active
            (Xfe_H_m, Xfe_C_m, Xfe_N_m, Xfe_S_m) = setup_marker_metal_volatile_properties(
                marknum;
                initial_h_ppm=cfg.metal_partition.initial_metal_h_ppm,
                initial_c_ppm=cfg.metal_partition.initial_metal_c_ppm,
                initial_n_ppm=cfg.metal_partition.initial_metal_n_ppm,
                initial_s_ppm=cfg.metal_partition.initial_metal_s_ppm,
            )
            Xfe_H_m_step_start = zeros(Float64, marknum)
            Xfe_C_m_step_start = zeros(Float64, marknum)
            Xfe_N_m_step_start = zeros(Float64, marknum)
            Xfe_S_m_step_start = zeros(Float64, marknum)
        end
        if cfg.phase_tracking.active
            phase_arrays = setup_marker_phase_tracking_properties(
                marknum, cfg.phase_tracking
            )
            Xmin_troilite_m = phase_arrays.Xmin_troilite_m
            Xmin_schreibersite_m = phase_arrays.Xmin_schreibersite_m
            Xmin_cohenite_m = phase_arrays.Xmin_cohenite_m
            Xmin_graphite_m = phase_arrays.Xmin_graphite_m
            Xmin_nitride_m = phase_arrays.Xmin_nitride_m
            Xmin_metal_matrix_m = phase_arrays.Xmin_metal_matrix_m
        end
        t_accreted = setup_marker_accretion_properties(
            marknum, cfg.accretion; initial_time=timesum
        )
        hcnspo_props = if cfg.volatile_mixture.active || cfg.refractory.active
            setup_marker_hcnspo_properties(marknum, cfg.volatile_mixture, cfg.refractory)
        else
            nothing
        end
        redox_props = if cfg.redox.active
            setup_marker_redox_properties(
                marknum, cfg.redox; initial_xfe_bulk=Xfe_bulk, tkm=tkm, pfm=pfm0
            )
        else
            nothing
        end
        define_markers!(
            xm,
            ym,
            tm,
            phim,
            etavpm,
            rhototalm,
            rhocptotalm,
            etatotalm,
            hrtotalm,
            ktotalm,
            tkm,
            inv_gggtotalm,
            fricttotalm,
            cohestotalm,
            tenstotalm,
            rhofluidcur,
            alphasolidcur,
            alphafluidcur,
            XWsolidm0;
            coords=coords,
            xcenter_val=xcenter_val,
            ycenter_val=ycenter_val,
            rplanet_val=rplanet_val,
            rcrust_val=rcrust_val,
            XWsolidm_init_val=cfg.materials.XWsolidm_init,
            phim0_val=phim0_val,
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            Xfem0=Xfem0,
            Xfe_bulk_val=Xfe_bulk_val,
            T_eutectic_val=T_eutectic_val,
            dT_metal_val=dT_metal_val,
            tkm0_val=cfg.materials.tkm0,
        )
        # copy thermodynamic marker properties to next generation for initial setup
        XWsolidm .= XWsolidm0
        phinewm .= phim

        if hcnspo_props !== nothing
            for m in 1:marknum
                if tm[m] == 2
                    if cfg.volatile_mixture.active
                        hcnspo_props.X_ice_H2O_m[m] = cfg.volatile_mixture.X_ice_H2O
                        hcnspo_props.X_ice_NH3_m[m] = cfg.volatile_mixture.X_ice_NH3
                        hcnspo_props.X_ice_CO2_m[m] = cfg.volatile_mixture.X_ice_CO2
                        hcnspo_props.X_ice_CO_m[m] = cfg.volatile_mixture.X_ice_CO
                        hcnspo_props.X_ice_CH4_m[m] = cfg.volatile_mixture.X_ice_CH4
                        hcnspo_props.X_ice_N2_m[m] = cfg.volatile_mixture.X_ice_N2
                        hcnspo_props.X_ice_H2S_m[m] = cfg.volatile_mixture.X_ice_H2S
                        hcnspo_props.X_ice_PH3_m[m] = cfg.volatile_mixture.X_ice_PH3
                    end
                    if cfg.refractory.active
                        hcnspo_props.X_refr_C_m[m] = cfg.refractory.f_refr_C
                        hcnspo_props.X_refr_S_m[m] = cfg.refractory.f_refr_S
                        hcnspo_props.X_refr_N_m[m] = cfg.refractory.f_refr_N
                        hcnspo_props.X_refr_P_m[m] = cfg.refractory.f_refr_P
                        hcnspo_props.X_refr_H_m[m] = cfg.refractory.f_refr_H
                    end
                end
            end
        end

        # save initial state
        if cfg.output.mode != :telemetry
            save_state(
                output_path,
                0,
                dt,
                timesum,
                marknum,
                ETA,
                ETA0,
                GGG,
                EXY,
                SXY,
                SXY0,
                wyx,
                COH,
                TEN,
                FRI,
                YNY,
                RHOX,
                RHOFX,
                KX,
                PHIX,
                vx,
                vxf,
                RX,
                qxD,
                gx,
                RHOY,
                RHOFY,
                KY,
                PHIY,
                vy,
                vyf,
                RY,
                qyD,
                gy,
                RHO,
                RHOCP,
                ALPHA,
                ALPHAF,
                HR,
                HA,
                HS,
                ETAP,
                GGGP,
                EXX,
                SXX,
                SXX0,
                tk1,
                tk2,
                vxp,
                vyp,
                vxpf,
                vypf,
                pr,
                pf,
                ps,
                pr0,
                pf0,
                ps0,
                ETAPHI,
                BETAPHI,
                PHI,
                APHI,
                FI,
                ETA5,
                ETA00,
                YNY5,
                YNY00,
                YNY_inv_ETA,
                DSXY,
                EII,
                SII,
                DSXX,
                DMP,
                DHP,
                DQPF,
                XWS,
                XWsolidm0,
                xm,
                ym,
                tm,
                tkm,
                sxxm,
                sxym,
                etavpm,
                phim,
                rhototalm,
                rhocptotalm,
                etatotalm,
                hrtotalm,
                ktotalm,
                tkm_rhocptotalm,
                etafluidcur_inv_kphim,
                inv_gggtotalm,
                fricttotalm,
                cohestotalm,
                tenstotalm,
                rhofluidcur,
                alphasolidcur,
                alphafluidcur;
                coords=coords,
                phim0_val=phim0_val,
                M_vent_total=M_vent_total,
                M_vent_H2O_total=M_vent_H2O_total,
                M_vent_C_total=M_vent_C_total,
                M_vent_N_total=M_vent_N_total,
                M_vent_S_total=M_vent_S_total,
                M_atm_total=M_atm_total,
                M_escaped_total=M_escaped_total,
                P_amb=cfg.disk.p_amb_disk,
                S_vent=S_vent_grid,
                Xfem=Xfem,
                Xfem0=Xfem0,
                Xfe_bulk=Xfe_bulk,
                M_atm_species=M_atm_species,
                M_escaped_species=M_escaped_species,
                XH2Om=XH2Om,
                XCm=XCm,
                XNm=XNm,
                XSm=XSm,
                Xfe_H_m=Xfe_H_m,
                Xfe_C_m=Xfe_C_m,
                Xfe_N_m=Xfe_N_m,
                Xfe_S_m=Xfe_S_m,
                core_budgets=core_budgets,
                Xmin_troilite_m=Xmin_troilite_m,
                Xmin_schreibersite_m=Xmin_schreibersite_m,
                Xmin_cohenite_m=Xmin_cohenite_m,
                Xmin_graphite_m=Xmin_graphite_m,
                Xmin_nitride_m=Xmin_nitride_m,
                Xmin_metal_matrix_m=Xmin_metal_matrix_m,
                regional_mineral_modes=regional_mineral_modes,
                DT0=DT0,
                rplanet=rplanet_val,
                rcore=rcore_val,
                t_accreted=t_accreted,
                M_accreted_total=cfg.accretion.active ? M_accreted_total : nothing,
                M_planet_val=M_planet_val,
                telescope_level=telescope_level,
                hcnspo_props=hcnspo_props,
                atm_state=atm_state,
            )
        end
    end

    # ---------------------------------------------------------------------
    # set up interpolation arrays"
    # ---------------------------------------------------------------------
    (ETA0SUM, ETASUM, GGGSUM, SXYSUM, COHSUM, TENSUM, FRISUM, WTSUM, RHOXSUM, RHOFXSUM, KXSUM, PHIXSUM, RXSUM, WTXSUM, RHOYSUM, RHOFYSUM, KYSUM, PHIYSUM, RYSUM, WTYSUM, RHOSUM, RHOCPSUM, ALPHASUM, ALPHAFSUM, HRSUM, GGGPSUM, SXXSUM, TKSUM, PHISUM, DMPSUM, DHPSUM, XWSSUM, WTPSUM) = setup_interpolated_properties(
        coords
    )

    # -------------------------------------------------------------------------
    # set up of matrices for global grav/thermal/hydromechanical solvers"
    # -------------------------------------------------------------------------
    # hydromechanical solver
    darcy_elim_val = cfg.solver.darcy_elimination
    dof_per_node_val = darcy_elim_val ? 4 : 6
    R, S = setup_hydromechanical_lse(coords; dof_per_node=dof_per_node_val)
    hydromech_sol = nothing
    hydromech_ws = HydromechanicalLSEWorkspace(coords; dof_per_node=dof_per_node_val)
    L_hydromech = hydromech_ws.L
    hydromech_cache = nothing
    # thermal solver
    RT, ST = setup_thermal_lse(coords)
    thermal_ws = ThermalLSEWorkspace(coords)
    LT_thermal = thermal_ws.LT
    thermal_cache = nothing
    # gravitational solver
    RP, SP = setup_gravitational_lse(coords)
    # precompute gravitational Poisson operator (invariant across timesteps)
    LP = assemble_gravitational_lse!(zeros(coords.Ny1, coords.Nx1), RP; coords=coords)
    F_grav = lu(LP.cscmatrix)
    # Pardiso MKL solver
    pardiso_solver = nothing
    pardiso_last_nnz = 0
    if use_pardiso_val
        pardiso_solver = Pardiso.MKLPardisoSolver()
        initialize_pardiso!(pardiso_solver, iparms_dict)
    end

    # Pre-allocated segregation workspaces
    metal_seg_ws = if coreformation_active_val
        MetalSegregationWorkspace(
            coords.Ny, coords.Nx; track_volatiles=cfg.metal_partition.active
        )
    else
        nothing
    end
    magma_seg_ws = if magma_active_val
        MagmaSegregationWorkspace(coords.Ny, coords.Nx)
    else
        nothing
    end

    # -------------------------------------------------------------------------
    # iterate timesteps"
    # -------------------------------------------------------------------------
    generate_showvalues(timestep, marknum, maxT, dt, timesum) =
        () -> [
            (:timestep, timestep),
            (:marknum, marknum),
            (:maxT_K, maxT),
            (:dt_s, dt),
            (:timesum_Ma, s_to_Ma(timesum)),
            (:to_go_Ma, s_to_Ma(endtime_val - timesum)),
        ]
    p = Progress(
        n_steps_val;
        showspeed=true,
        dt=0.5,
        barglyphs=BarGlyphs('|', '█', ['▁', '▂', '▃', '▄', '▅', '▆', '▇'], ' ', '|'),
        barlen=10,
    )
    if (cfg.output.mode in (:telemetry, :both))
        telemetry_io = init_telemetry(
            output_path, cfg.output.telemetry_file; append=is_restart
        )
    end
    try
        for timestep in start_step_val:1:n_steps_val
            timestep_begin = now()
            # ---------------------------------------------------------------------
            # reset interpolation arrays
            # ---------------------------------------------------------------------
            reset_interpolated_properties!(
                ETA0SUM,
                ETASUM,
                GGGSUM,
                SXYSUM,
                COHSUM,
                TENSUM,
                FRISUM,
                WTSUM,
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
                RHOSUM,
                RHOCPSUM,
                ALPHASUM,
                ALPHAFSUM,
                HRSUM,
                GGGPSUM,
                SXXSUM,
                TKSUM,
                PHISUM,
                WTPSUM,
            )

            # ---------------------------------------------------------------------
            # compute ambient conditions and update sticky air markers
            # ---------------------------------------------------------------------
            T_amb, P_amb, w_disp = compute_ambient_conditions(timesum, cfg.disk)
            isfinite(T_amb) || throw(
                DomainError(T_amb, "Ambient disk temperature must be finite, got $T_amb"),
            )
            P_atm = if cfg.atmosphere.active && atm_state !== nothing
                atm_state.P_surf
            elseif cfg.escape.active
                compute_surface_atmospheric_pressure(M_atm_total, M_planet_val, rplanet_val)
            else
                0.0
            end
            P_amb_eff = P_amb + P_atm
            if disk_enabled_val || surface_radiation_val
                @threads :static for m in 1:marknum
                    if tm[m] >= 3
                        tkm[m] = T_amb
                    end
                end
            end

            # ---------------------------------------------------------------------
            # planetesimal accretion engine: mass addition, heating, boundary advance
            # ---------------------------------------------------------------------
            if cfg.accretion.active
                dM_dt_acc = compute_accretion_rate(
                    timesum, M_planet_val, rplanet_val, cfg.accretion, cfg.disk
                )
                # Clamp mass increment so M_planet_val does not overshoot M_target
                dM_remain = max(0.0, cfg.accretion.M_target - M_planet_val)
                dM_acc = min(dM_dt_acc * dt, dM_remain)

                if dM_acc > 0.0 && rplanet_val < cfg.accretion.R_target
                    dR_acc = compute_radius_increment(
                        rplanet_val, dM_acc, cfg.accretion.rho_bulk
                    )
                    # Clamp radius increment so rplanet_val does not overshoot R_target
                    dR_remain = max(0.0, cfg.accretion.R_target - rplanet_val)
                    dR_acc = min(dR_acc, dR_remain)
                    T_acc = T_amb
                    if cfg.accretion.h_impact > 0.0
                        _, delta_T_imp = compute_impact_heating(
                            M_planet_val,
                            rplanet_val;
                            h_impact=cfg.accretion.h_impact,
                            c_p=cfg.accretion.cp_rock,
                            v_inf=cfg.accretion.v_inf,
                        )
                        T_acc += delta_T_imp
                    end
                    XW_acc = cfg.accretion.XWsolid_dry
                    H2O_acc = cfg.accretion.XH2O_dry_wtpct
                    if cfg.accretion.snowline_coupling
                        XW_acc, H2O_acc = evaluate_snowline_water_content(
                            T_amb;
                            T_snowline_cond=cfg.accretion.T_snowline_cond,
                            XW_wet=cfg.accretion.XWsolid_wet,
                            XW_dry=cfg.accretion.XWsolid_dry,
                            H2O_wet_wtpct=cfg.accretion.XH2O_wet_wtpct,
                            H2O_dry_wtpct=cfg.accretion.XH2O_dry_wtpct,
                        )
                    end

                    cond_state_acc =
                        if (cfg.volatile_mixture.active || cfg.refractory.active)
                            evaluate_disk_volatile_condensation(
                                T_amb,
                                P_amb,
                                cfg.volatile_mixture,
                                cfg.refractory;
                                P_ref=cfg.volatile_mixture.P_ref,
                                alpha_P=cfg.volatile_mixture.alpha_P,
                            )
                        else
                            nothing
                        end

                    if cfg.volatile_mixture.active && cond_state_acc !== nothing
                        if cond_state_acc.condensed_H2O
                            XW_acc = cfg.volatile_mixture.X_ice_H2O
                            H2O_acc = cfg.volatile_mixture.X_ice_H2O * 100.0
                        else
                            XW_acc = cfg.accretion.XWsolid_dry
                            H2O_acc = cfg.accretion.XH2O_dry_wtpct
                        end
                    end

                    advance_accretion_boundary!(
                        rplanet_val,
                        dR_acc,
                        xm,
                        ym,
                        tm,
                        tkm,
                        phim,
                        XWsolidm0,
                        Xfe_bulk,
                        Xfem;
                        xcenter=xcenter_val,
                        ycenter=ycenter_val,
                        T_accreted=T_acc,
                        phi_accreted=cfg.accretion.phi_accreted,
                        XWsolid_accreted=XW_acc,
                        Xfe_accreted=cfg.accretion.Xfe_bulk_accreted,
                        t_accreted=t_accreted,
                        current_time=timesum,
                        XWsolidm=XWsolidm,
                        phinewm=phinewm,
                        XH2Om=cfg.volatiles.active ? XH2Om : nothing,
                        XCm=cfg.volatiles.active ? XCm : nothing,
                        XNm=cfg.volatiles.active ? XNm : nothing,
                        XSm=cfg.volatiles.active ? XSm : nothing,
                        XH2O_accreted=H2O_acc,
                        XC_accreted=cfg.accretion.XC_accreted_ppm,
                        XN_accreted=cfg.accretion.XN_accreted_ppm,
                        XS_accreted=cfg.accretion.XS_accreted_ppm,
                        hcnspo_props=hcnspo_props,
                        disk_state=cond_state_acc,
                    )
                    rplanet_val += dR_acc
                    M_planet_val += dM_acc
                    M_accreted_total += dM_acc
                end
            end

            # ---------------------------------------------------------------------
            # telescoping domain expansion
            # ---------------------------------------------------------------------
            if cfg.telescoping.active && should_telescope_domain(
                rplanet_val, coords, cfg.telescoping; level=telescope_level
            )
                @info "Triggering domain telescoping expansion" level=telescope_level + 1 rplanet=rplanet_val old_xsize=coords.xsize new_xsize=2.0 *
                                                                                                                                               coords.xsize
                old_coords = coords
                coords = compute_telescoped_coordinates(old_coords)
                xcenter_val = coords.xcenter
                ycenter_val = coords.ycenter

                # Remap staggered grid variables
                ETA = remap_staggered_grid_array(
                    ETA, (coords.Ny, coords.Nx); background_val=cfg.materials.etasolidm[3]
                )
                ETA0 = remap_staggered_grid_array(
                    ETA0, (coords.Ny, coords.Nx); background_val=cfg.materials.etasolidm[3]
                )
                GGG = remap_staggered_grid_array(
                    GGG, (coords.Ny, coords.Nx); background_val=cfg.materials.gggsolidm[3]
                )
                EXY = remap_staggered_grid_array(EXY, (coords.Ny, coords.Nx))
                SXY = remap_staggered_grid_array(SXY, (coords.Ny, coords.Nx))
                SXY0 = remap_staggered_grid_array(SXY0, (coords.Ny, coords.Nx))
                wyx = remap_staggered_grid_array(wyx, (coords.Ny, coords.Nx))
                COH = remap_staggered_grid_array(
                    COH, (coords.Ny, coords.Nx); background_val=cfg.materials.cohessolidm[3]
                )
                TEN = remap_staggered_grid_array(
                    TEN, (coords.Ny, coords.Nx); background_val=cfg.materials.tenssolidm[3]
                )
                FRI = remap_staggered_grid_array(
                    FRI, (coords.Ny, coords.Nx); background_val=cfg.materials.frictsolidm[3]
                )
                YNY = remap_staggered_grid_array(YNY, (coords.Ny, coords.Nx))
                ETA5 = remap_staggered_grid_array(
                    ETA5, (coords.Ny, coords.Nx); background_val=cfg.materials.etasolidm[3]
                )
                ETA00 = remap_staggered_grid_array(
                    ETA00, (coords.Ny, coords.Nx); background_val=cfg.materials.etasolidm[3]
                )
                YNY5 = remap_staggered_grid_array(YNY5, (coords.Ny, coords.Nx))
                YNY00 = remap_staggered_grid_array(YNY00, (coords.Ny, coords.Nx))
                YNY_inv_ETA = remap_staggered_grid_array(
                    YNY_inv_ETA, (coords.Ny, coords.Nx)
                )
                DSXY = remap_staggered_grid_array(DSXY, (coords.Ny, coords.Nx))
                DSY = remap_staggered_grid_array(DSY, (coords.Ny, coords.Nx))

                RHOX = remap_staggered_grid_array(
                    RHOX,
                    (coords.Ny1, coords.Nx1);
                    background_val=cfg.materials.rhosolidm[3],
                )
                RHOFX = remap_staggered_grid_array(
                    RHOFX,
                    (coords.Ny1, coords.Nx1);
                    background_val=cfg.materials.rhofluidm[3],
                )
                KX = remap_staggered_grid_array(
                    KX, (coords.Ny1, coords.Nx1); background_val=cfg.materials.ksolidm[3]
                )
                PHIX = remap_staggered_grid_array(
                    PHIX, (coords.Ny1, coords.Nx1); background_val=cfg.poroelasticity.phimin
                )
                vx = remap_staggered_grid_array(vx, (coords.Ny1, coords.Nx1))
                vxf = remap_staggered_grid_array(vxf, (coords.Ny1, coords.Nx1))
                RX = remap_staggered_grid_array(RX, (coords.Ny1, coords.Nx1))
                qxD = remap_staggered_grid_array(qxD, (coords.Ny1, coords.Nx1))
                gx = remap_staggered_grid_array(gx, (coords.Ny1, coords.Nx1))

                RHOY = remap_staggered_grid_array(
                    RHOY,
                    (coords.Ny1, coords.Nx1);
                    background_val=cfg.materials.rhosolidm[3],
                )
                RHOFY = remap_staggered_grid_array(
                    RHOFY,
                    (coords.Ny1, coords.Nx1);
                    background_val=cfg.materials.rhofluidm[3],
                )
                KY = remap_staggered_grid_array(
                    KY, (coords.Ny1, coords.Nx1); background_val=cfg.materials.ksolidm[3]
                )
                PHIY = remap_staggered_grid_array(
                    PHIY, (coords.Ny1, coords.Nx1); background_val=cfg.poroelasticity.phimin
                )
                vy = remap_staggered_grid_array(vy, (coords.Ny1, coords.Nx1))
                vyf = remap_staggered_grid_array(vyf, (coords.Ny1, coords.Nx1))
                RY = remap_staggered_grid_array(RY, (coords.Ny1, coords.Nx1))
                qyD = remap_staggered_grid_array(qyD, (coords.Ny1, coords.Nx1))
                gy = remap_staggered_grid_array(gy, (coords.Ny1, coords.Nx1))

                RHO = remap_staggered_grid_array(
                    RHO, (coords.Ny1, coords.Nx1); background_val=cfg.materials.rhosolidm[3]
                )
                RHOCP = remap_staggered_grid_array(
                    RHOCP,
                    (coords.Ny1, coords.Nx1);
                    background_val=cfg.materials.rhocpsolidm[3],
                )
                ALPHA = remap_staggered_grid_array(
                    ALPHA,
                    (coords.Ny1, coords.Nx1);
                    background_val=cfg.materials.alphasolidm[3],
                )
                ALPHAF = remap_staggered_grid_array(
                    ALPHAF,
                    (coords.Ny1, coords.Nx1);
                    background_val=cfg.materials.alphafluidm[3],
                )
                HR = remap_staggered_grid_array(HR, (coords.Ny1, coords.Nx1))
                HA = remap_staggered_grid_array(HA, (coords.Ny1, coords.Nx1))
                HS = remap_staggered_grid_array(HS, (coords.Ny1, coords.Nx1))
                ETAP = remap_staggered_grid_array(
                    ETAP,
                    (coords.Ny1, coords.Nx1);
                    background_val=cfg.materials.etasolidm[3],
                )
                GGGP = remap_staggered_grid_array(
                    GGGP,
                    (coords.Ny1, coords.Nx1);
                    background_val=cfg.materials.gggsolidm[3],
                )
                EXX = remap_staggered_grid_array(EXX, (coords.Ny1, coords.Nx1))
                SXX = remap_staggered_grid_array(SXX, (coords.Ny1, coords.Nx1))
                SXX0 = remap_staggered_grid_array(SXX0, (coords.Ny1, coords.Nx1))
                tk1 = remap_staggered_grid_array(
                    tk1, (coords.Ny1, coords.Nx1); background_val=cfg.materials.tkm0[3]
                )
                tk2 = remap_staggered_grid_array(
                    tk2, (coords.Ny1, coords.Nx1); background_val=cfg.materials.tkm0[3]
                )
                DT = remap_staggered_grid_array(DT, (coords.Ny1, coords.Nx1))
                DT0 = remap_staggered_grid_array(DT0, (coords.Ny1, coords.Nx1))
                vxp = remap_staggered_grid_array(vxp, (coords.Ny1, coords.Nx1))
                vyp = remap_staggered_grid_array(vyp, (coords.Ny1, coords.Nx1))
                vxpf = remap_staggered_grid_array(vxpf, (coords.Ny1, coords.Nx1))
                vypf = remap_staggered_grid_array(vypf, (coords.Ny1, coords.Nx1))
                pr = remap_staggered_grid_array(pr, (coords.Ny1, coords.Nx1))
                pf = remap_staggered_grid_array(pf, (coords.Ny1, coords.Nx1))
                ps = remap_staggered_grid_array(ps, (coords.Ny1, coords.Nx1))
                pr0 = remap_staggered_grid_array(pr0, (coords.Ny1, coords.Nx1))
                pf0 = remap_staggered_grid_array(pf0, (coords.Ny1, coords.Nx1))
                ps0 = remap_staggered_grid_array(ps0, (coords.Ny1, coords.Nx1))
                ETAPHI = remap_staggered_grid_array(
                    ETAPHI,
                    (coords.Ny1, coords.Nx1);
                    background_val=cfg.materials.etasolidm[3],
                )
                BETAPHI = remap_staggered_grid_array(BETAPHI, (coords.Ny1, coords.Nx1))
                PHI = remap_staggered_grid_array(
                    PHI, (coords.Ny1, coords.Nx1); background_val=cfg.poroelasticity.phimin
                )
                APHI = remap_staggered_grid_array(APHI, (coords.Ny1, coords.Nx1))
                FI = remap_staggered_grid_array(FI, (coords.Ny1, coords.Nx1))
                DMP = remap_staggered_grid_array(DMP, (coords.Ny1, coords.Nx1))
                DHP = remap_staggered_grid_array(DHP, (coords.Ny1, coords.Nx1))
                XWS = remap_staggered_grid_array(XWS, (coords.Ny1, coords.Nx1))
                EII = remap_staggered_grid_array(EII, (coords.Ny1, coords.Nx1))
                SII = remap_staggered_grid_array(SII, (coords.Ny1, coords.Nx1))
                DSXX = remap_staggered_grid_array(DSXX, (coords.Ny1, coords.Nx1))
                tk0 = remap_staggered_grid_array(
                    tk0, (coords.Ny1, coords.Nx1); background_val=cfg.materials.tkm0[3]
                )
                if Q_metric !== nothing
                    Q_metric = remap_staggered_grid_array(
                        Q_metric, (coords.Ny1, coords.Nx1)
                    )
                end
                DQPF = remap_staggered_grid_array(DQPF, (coords.Ny1, coords.Nx1))
                DQPFSUM = remap_staggered_grid_array(DQPFSUM, (coords.Ny1, coords.Nx1))
                S_vent_grid = remap_staggered_grid_array(
                    S_vent_grid, (coords.Ny1, coords.Nx1)
                )
                Q_lat_grid = remap_staggered_grid_array(
                    Q_lat_grid, (coords.Ny1, coords.Nx1)
                )
                Q_seg_grid = remap_staggered_grid_array(
                    Q_seg_grid, (coords.Ny1, coords.Nx1)
                )

                # Helper geometries and buffers
                mdis, mnum = setup_marker_geometry_helpers(coords)
                (ETA0SUM, ETASUM, GGGSUM, SXYSUM, COHSUM, TENSUM, FRISUM, WTSUM, RHOXSUM, RHOFXSUM, KXSUM, PHIXSUM, RXSUM, WTXSUM, RHOYSUM, RHOFYSUM, KYSUM, PHIYSUM, RYSUM, WTYSUM, RHOSUM, RHOCPSUM, ALPHASUM, ALPHAFSUM, HRSUM, GGGPSUM, SXXSUM, TKSUM, PHISUM, DMPSUM, DHPSUM, XWSSUM, WTPSUM) = setup_interpolated_properties(
                    coords
                )
                if !use_tiled_p2m && use_threading
                    thread_buffers = allocate_thread_interpolation_buffers(
                        num_buffers, coords
                    )
                end

                # Telescope marker arrays
                marknum = telescope_marker_arrays!(
                    xm,
                    ym,
                    tm,
                    tkm,
                    sxxm,
                    sxym,
                    etavpm,
                    phim,
                    phinewm,
                    pfm0,
                    XWsolidm,
                    XWsolidm0,
                    Fm,
                    rhototalm,
                    rhocptotalm,
                    etatotalm,
                    hrtotalm,
                    ktotalm,
                    inv_gggtotalm,
                    fricttotalm,
                    cohestotalm,
                    tenstotalm,
                    rhofluidcur,
                    alphasolidcur,
                    alphafluidcur,
                    tkm_rhocptotalm,
                    etafluidcur_inv_kphim;
                    old_coords=old_coords,
                    new_coords=coords,
                    T_ambient=cfg.materials.tkm0[3],
                    phi_ambient=cfg.poroelasticity.phimin,
                    buffer_markers_per_cell=cfg.telescoping.buffer_markers_per_cell,
                    materials=cfg.materials,
                    Xfem=Xfem,
                    Xfem0=Xfem0,
                    Xfe_bulk=Xfe_bulk,
                    XH2Om=cfg.volatiles.active ? XH2Om : nothing,
                    XCm=cfg.volatiles.active ? XCm : nothing,
                    XNm=cfg.volatiles.active ? XNm : nothing,
                    XSm=cfg.volatiles.active ? XSm : nothing,
                    Xfe_H_m=Xfe_H_m,
                    Xfe_C_m=Xfe_C_m,
                    Xfe_N_m=Xfe_N_m,
                    Xfe_S_m=Xfe_S_m,
                    Xmin_troilite_m=Xmin_troilite_m,
                    Xmin_schreibersite_m=Xmin_schreibersite_m,
                    Xmin_cohenite_m=Xmin_cohenite_m,
                    Xmin_graphite_m=Xmin_graphite_m,
                    Xmin_nitride_m=Xmin_nitride_m,
                    Xmin_metal_matrix_m=Xmin_metal_matrix_m,
                    t_accreted=t_accreted,
                    hcnspo_props=hcnspo_props,
                    F_extract_m=F_extract_m,
                )

                # Re-initialize linear solvers and Poisson operator for new grid size
                R, S = setup_hydromechanical_lse(coords; dof_per_node=dof_per_node_val)
                hydromech_sol = nothing
                hydromech_ws = HydromechanicalLSEWorkspace(
                    coords; dof_per_node=dof_per_node_val
                )
                L_hydromech = hydromech_ws.L
                hydromech_cache = nothing

                RT, ST = setup_thermal_lse(coords)
                thermal_ws = ThermalLSEWorkspace(coords)
                LT_thermal = thermal_ws.LT
                thermal_cache = nothing

                RP, SP = setup_gravitational_lse(coords)
                LP = assemble_gravitational_lse!(
                    zeros(coords.Ny1, coords.Nx1), RP; coords=coords
                )
                F_grav = lu(LP.cscmatrix)

                if use_pardiso_val && pardiso_solver !== nothing
                    set_phase!(pardiso_solver, Pardiso.RELEASE_ALL)
                    pardiso(pardiso_solver)
                    pardiso_last_nnz = 0
                end

                telescope_level += 1
                if metal_seg_ws !== nothing
                    metal_seg_ws = MetalSegregationWorkspace(
                        coords.Ny, coords.Nx; track_volatiles=cfg.metal_partition.active
                    )
                end
                if magma_seg_ws !== nothing
                    magma_seg_ws = MagmaSegregationWorkspace(coords.Ny, coords.Nx)
                end
                @info "Telescoping complete" level=telescope_level Nx=coords.Nx Ny=coords.Ny marknum=marknum
            end

            # ---------------------------------------------------------------------
            # calculate radioactive heating
            # ---------------------------------------------------------------------
            hrsolidm, hrfluidm, hrmetalm = calculate_radioactive_heating(
                hr_al_val,
                hr_fe_val,
                timesum;
                ratio_al=ratio_al_val,
                E_al=E_al_val,
                f_al=f_al_val,
                tau_al=tau_al_val,
                ratio_fe=ratio_fe_val,
                E_fe=E_fe_val,
                f_fe=f_fe_val,
                tau_fe=tau_fe_val,
                rho_metal=rho_metal_val,
            )

            # ---------------------------------------------------------------------
            # compute marker properties and interpolate to staggered grid
            # ---------------------------------------------------------------------
            if cfg.redox.active && redox_props !== nothing
                update_marker_redox!(
                    redox_props, tkm, pfm0, cfg.redox; Xfem=Xfem, XWsolidm=XWsolidm
                )
            end

            if use_tiled_p2m
                p2m_workspace = ensure_workspace_compatible(
                    p2m_workspace, coords, marknum, cfg.solver.tile_size
                )
                bin_markers_into_tiles!(p2m_workspace, xm, ym, coords, marknum)
                for color in 1:4
                    tiles = p2m_workspace.tiles_by_color[color]
                    Threads.@threads :dynamic for t in tiles
                        lo = p2m_workspace.tile_offsets[t]
                        hi = p2m_workspace.tile_offsets[t + 1] - 1
                        lo > hi && continue
                        for idx in lo:hi
                            m = p2m_workspace.tile_markers[idx]
                            compute_marker_properties!(
                                m,
                                tm,
                                tkm,
                                rhototalm,
                                rhocptotalm,
                                etatotalm,
                                hrtotalm,
                                ktotalm,
                                tkm_rhocptotalm,
                                etafluidcur_inv_kphim,
                                hrsolidm,
                                hrfluidm,
                                phim,
                                XWsolidm0,
                                marker_property_mode,
                                rhofluidcur;
                                thermal_buoyancy=thermal_buoyancy_val,
                                alphafluid=alphafluid_val,
                                tmfluidphase_val=tmfluidphase_val,
                                fluid_viscosity_mode=fluid_viscosity_mode_val,
                                fluid_viscosity_Ea=fluid_viscosity_Ea_val,
                                fluid_viscosity_T0=fluid_viscosity_T0_val,
                                fluid_viscosity_eta0=fluid_viscosity_eta0_val,
                                pm=pfm0,
                                Fm=Fm,
                                melting_active=melting_active_val,
                                magma_transport_active=magma_active_val,
                                track_depletion=cfg.magma_transport.track_depletion,
                                F_extract_m=F_extract_m,
                                T_solidus_val=T_solidus_val,
                                T_liquidus_val=T_liquidus_val,
                                L_melt_val=L_melt_val,
                                rho_melt_val=rho_melt_val,
                                alpha_eta_val=alpha_eta_val,
                                phi_crit_val=phi_crit_val,
                                eta_melt_val=eta_melt_val,
                                dpdt_clapeyron_val=dpdt_clapeyron_val,
                                soft_turbulence=soft_turbulence_val,
                                eta_fluid_silicate_val=eta_fluid_silicate_val,
                                F_turb_start_val=F_turb_start_val,
                                F_turb_end_val=F_turb_end_val,
                                turb_exponent_val=turb_exponent_val,
                                dT_turb_min_val=dT_turb_min_val,
                                T_surface_ref_val=T_surface_ref_val,
                                k_turb_cutoff_val=k_turb_cutoff_val,
                                k_turb_floor_val=k_turb_floor_val,
                                Xfe_bulk=Xfe_bulk,
                                Xfem=Xfem,
                                coreformation_active=coreformation_active_val,
                                hrmetalm=hrmetalm,
                                sulfur_fraction_val=sulfur_fraction_val,
                                metal_density_mode_val=metal_density_mode_val,
                                T_eutectic_val=T_eutectic_val,
                                dT_metal_val=dT_metal_val,
                                rho_metal_val=rho_metal_val,
                                rho_metal_solid_val=rho_metal_solid_val,
                                L_metal_val=L_metal_val,
                                k_metal_val=k_metal_val,
                                rhocp_metal_val=rhocp_metal_val,
                                volatiles_active=cfg.volatiles.active,
                                volatiles_cfg=cfg.volatiles,
                                retention_cfg=cfg.retention,
                                XH2Om=XH2Om,
                                XCm=XCm,
                                XNm=XNm,
                                XSm=XSm,
                                metal_partition_cfg=cfg.metal_partition,
                                Xfe_H_m=Xfe_H_m,
                                Xfe_C_m=Xfe_C_m,
                                Xfe_N_m=Xfe_N_m,
                                Xfe_S_m=Xfe_S_m,
                                phase_tracking_cfg=cfg.phase_tracking,
                                Xmin_troilite_m=Xmin_troilite_m,
                                Xmin_schreibersite_m=Xmin_schreibersite_m,
                                Xmin_cohenite_m=Xmin_cohenite_m,
                                Xmin_graphite_m=Xmin_graphite_m,
                                Xmin_nitride_m=Xmin_nitride_m,
                                Xmin_metal_matrix_m=Xmin_metal_matrix_m,
                                hydrothermal_active=cfg.hydrothermal.active,
                                hydrothermal_cfg=cfg.hydrothermal,
                                deltaIW_m=if redox_props !== nothing
                                    redox_props.deltaIW_m
                                else
                                    nothing
                                end,
                                xm=xm,
                                ym=ym,
                                coords=coords,
                                rplanet_val=rplanet_val,
                                M_planet_val=M_planet_val,
                                rcore_val=rcore_val,
                                qxD_val=qxD,
                                qyD_val=qyD,
                                xcenter_val=xcenter_val,
                                ycenter_val=ycenter_val,
                            )
                            scatter_marker_to_master_grids!(
                                m,
                                xm[m],
                                ym[m],
                                coords,
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
                        end
                    end
                end
            elseif use_threading
                reset_thread_buffers!(thread_buffers)
                nchunks = length(thread_buffers)
                Threads.@threads :static for c in 1:nchunks
                    buf = thread_buffers[c]
                    lo = (c - 1) * div(marknum, nchunks) + 1
                    hi = c == nchunks ? marknum : c * div(marknum, nchunks)
                    for m in lo:hi
                        compute_marker_properties!(
                            m,
                            tm,
                            tkm,
                            rhototalm,
                            rhocptotalm,
                            etatotalm,
                            hrtotalm,
                            ktotalm,
                            tkm_rhocptotalm,
                            etafluidcur_inv_kphim,
                            hrsolidm,
                            hrfluidm,
                            phim,
                            XWsolidm0,
                            marker_property_mode,
                            rhofluidcur;
                            thermal_buoyancy=thermal_buoyancy_val,
                            alphafluid=alphafluid_val,
                            tmfluidphase_val=tmfluidphase_val,
                            fluid_viscosity_mode=fluid_viscosity_mode_val,
                            fluid_viscosity_Ea=fluid_viscosity_Ea_val,
                            fluid_viscosity_T0=fluid_viscosity_T0_val,
                            fluid_viscosity_eta0=fluid_viscosity_eta0_val,
                            pm=pfm0,
                            Fm=Fm,
                            melting_active=melting_active_val,
                            magma_transport_active=magma_active_val,
                            track_depletion=cfg.magma_transport.track_depletion,
                            F_extract_m=F_extract_m,
                            T_solidus_val=T_solidus_val,
                            T_liquidus_val=T_liquidus_val,
                            L_melt_val=L_melt_val,
                            rho_melt_val=rho_melt_val,
                            alpha_eta_val=alpha_eta_val,
                            phi_crit_val=phi_crit_val,
                            eta_melt_val=eta_melt_val,
                            dpdt_clapeyron_val=dpdt_clapeyron_val,
                            soft_turbulence=soft_turbulence_val,
                            eta_fluid_silicate_val=eta_fluid_silicate_val,
                            F_turb_start_val=F_turb_start_val,
                            F_turb_end_val=F_turb_end_val,
                            turb_exponent_val=turb_exponent_val,
                            dT_turb_min_val=dT_turb_min_val,
                            T_surface_ref_val=T_surface_ref_val,
                            k_turb_cutoff_val=k_turb_cutoff_val,
                            k_turb_floor_val=k_turb_floor_val,
                            Xfe_bulk=Xfe_bulk,
                            Xfem=Xfem,
                            coreformation_active=coreformation_active_val,
                            hrmetalm=hrmetalm,
                            sulfur_fraction_val=sulfur_fraction_val,
                            metal_density_mode_val=metal_density_mode_val,
                            T_eutectic_val=T_eutectic_val,
                            dT_metal_val=dT_metal_val,
                            rho_metal_val=rho_metal_val,
                            rho_metal_solid_val=rho_metal_solid_val,
                            L_metal_val=L_metal_val,
                            k_metal_val=k_metal_val,
                            rhocp_metal_val=rhocp_metal_val,
                            volatiles_active=cfg.volatiles.active,
                            volatiles_cfg=cfg.volatiles,
                            retention_cfg=cfg.retention,
                            XH2Om=XH2Om,
                            XCm=XCm,
                            XNm=XNm,
                            XSm=XSm,
                            metal_partition_cfg=cfg.metal_partition,
                            Xfe_H_m=Xfe_H_m,
                            Xfe_C_m=Xfe_C_m,
                            Xfe_N_m=Xfe_N_m,
                            Xfe_S_m=Xfe_S_m,
                            phase_tracking_cfg=cfg.phase_tracking,
                            Xmin_troilite_m=Xmin_troilite_m,
                            Xmin_schreibersite_m=Xmin_schreibersite_m,
                            Xmin_cohenite_m=Xmin_cohenite_m,
                            Xmin_graphite_m=Xmin_graphite_m,
                            Xmin_nitride_m=Xmin_nitride_m,
                            Xmin_metal_matrix_m=Xmin_metal_matrix_m,
                            hydrothermal_active=cfg.hydrothermal.active,
                            hydrothermal_cfg=cfg.hydrothermal,
                            deltaIW_m=if redox_props !== nothing
                                redox_props.deltaIW_m
                            else
                                nothing
                            end,
                            xm=xm,
                            ym=ym,
                            coords=coords,
                            rplanet_val=rplanet_val,
                            M_planet_val=M_planet_val,
                            rcore_val=rcore_val,
                            qxD_val=qxD,
                            qyD_val=qyD,
                            xcenter_val=xcenter_val,
                            ycenter_val=ycenter_val,
                        )
                        @inbounds marker_to_basic_nodes!(
                            m,
                            xm[m],
                            ym[m],
                            etatotalm,
                            etavpm,
                            inv_gggtotalm,
                            sxym,
                            cohestotalm,
                            tenstotalm,
                            fricttotalm,
                            buf.ETA0SUM,
                            buf.ETASUM,
                            buf.GGGSUM,
                            buf.SXYSUM,
                            buf.COHSUM,
                            buf.TENSUM,
                            buf.FRISUM,
                            buf.WTSUM;
                            coords=coords,
                        )
                        @inbounds marker_to_vx_nodes!(
                            m,
                            xm[m],
                            ym[m],
                            rhototalm,
                            rhofluidcur,
                            ktotalm,
                            phim,
                            etafluidcur_inv_kphim,
                            buf.RHOXSUM,
                            buf.RHOFXSUM,
                            buf.KXSUM,
                            buf.PHIXSUM,
                            buf.RXSUM,
                            buf.WTXSUM;
                            coords=coords,
                        )
                        @inbounds marker_to_vy_nodes!(
                            m,
                            xm[m],
                            ym[m],
                            rhototalm,
                            rhofluidcur,
                            ktotalm,
                            phim,
                            etafluidcur_inv_kphim,
                            buf.RHOYSUM,
                            buf.RHOFYSUM,
                            buf.KYSUM,
                            buf.PHIYSUM,
                            buf.RYSUM,
                            buf.WTYSUM;
                            coords=coords,
                        )
                        @inbounds marker_to_p_nodes!(
                            m,
                            xm[m],
                            ym[m],
                            inv_gggtotalm,
                            sxxm,
                            rhototalm,
                            rhocptotalm,
                            alphasolidcur,
                            alphafluidcur,
                            hrtotalm,
                            phim,
                            tkm_rhocptotalm,
                            buf.GGGPSUM,
                            buf.SXXSUM,
                            buf.RHOSUM,
                            buf.RHOCPSUM,
                            buf.ALPHASUM,
                            buf.ALPHAFSUM,
                            buf.HRSUM,
                            buf.PHISUM,
                            buf.TKSUM,
                            buf.WTPSUM;
                            coords=coords,
                        )
                    end
                end
                reduce_thread_buffers!(
                    ETA0SUM,
                    ETASUM,
                    GGGSUM,
                    SXYSUM,
                    COHSUM,
                    TENSUM,
                    FRISUM,
                    WTSUM,
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
                    RHOSUM,
                    RHOCPSUM,
                    ALPHASUM,
                    ALPHAFSUM,
                    HRSUM,
                    GGGPSUM,
                    SXXSUM,
                    TKSUM,
                    PHISUM,
                    WTPSUM,
                    thread_buffers,
                )
            else
                for m in 1:1:marknum
                    compute_marker_properties!(
                        m,
                        tm,
                        tkm,
                        rhototalm,
                        rhocptotalm,
                        etatotalm,
                        hrtotalm,
                        ktotalm,
                        tkm_rhocptotalm,
                        etafluidcur_inv_kphim,
                        hrsolidm,
                        hrfluidm,
                        phim,
                        XWsolidm0,
                        marker_property_mode,
                        rhofluidcur;
                        thermal_buoyancy=thermal_buoyancy_val,
                        alphafluid=alphafluid_val,
                        tmfluidphase_val=tmfluidphase_val,
                        fluid_viscosity_mode=fluid_viscosity_mode_val,
                        fluid_viscosity_Ea=fluid_viscosity_Ea_val,
                        fluid_viscosity_T0=fluid_viscosity_T0_val,
                        fluid_viscosity_eta0=fluid_viscosity_eta0_val,
                        pm=pfm0,
                        Fm=Fm,
                        melting_active=melting_active_val,
                        magma_transport_active=magma_active_val,
                        track_depletion=cfg.magma_transport.track_depletion,
                        F_extract_m=F_extract_m,
                        T_solidus_val=T_solidus_val,
                        T_liquidus_val=T_liquidus_val,
                        L_melt_val=L_melt_val,
                        rho_melt_val=rho_melt_val,
                        alpha_eta_val=alpha_eta_val,
                        phi_crit_val=phi_crit_val,
                        eta_melt_val=eta_melt_val,
                        dpdt_clapeyron_val=dpdt_clapeyron_val,
                        soft_turbulence=soft_turbulence_val,
                        eta_fluid_silicate_val=eta_fluid_silicate_val,
                        F_turb_start_val=F_turb_start_val,
                        F_turb_end_val=F_turb_end_val,
                        turb_exponent_val=turb_exponent_val,
                        dT_turb_min_val=dT_turb_min_val,
                        T_surface_ref_val=T_surface_ref_val,
                        k_turb_cutoff_val=k_turb_cutoff_val,
                        k_turb_floor_val=k_turb_floor_val,
                        Xfe_bulk=Xfe_bulk,
                        Xfem=Xfem,
                        coreformation_active=coreformation_active_val,
                        hrmetalm=hrmetalm,
                        sulfur_fraction_val=sulfur_fraction_val,
                        metal_density_mode_val=metal_density_mode_val,
                        T_eutectic_val=T_eutectic_val,
                        dT_metal_val=dT_metal_val,
                        rho_metal_val=rho_metal_val,
                        rho_metal_solid_val=rho_metal_solid_val,
                        L_metal_val=L_metal_val,
                        k_metal_val=k_metal_val,
                        rhocp_metal_val=rhocp_metal_val,
                        volatiles_active=cfg.volatiles.active,
                        volatiles_cfg=cfg.volatiles,
                        retention_cfg=cfg.retention,
                        XH2Om=XH2Om,
                        XCm=XCm,
                        XNm=XNm,
                        XSm=XSm,
                        metal_partition_cfg=cfg.metal_partition,
                        Xfe_H_m=Xfe_H_m,
                        Xfe_C_m=Xfe_C_m,
                        Xfe_N_m=Xfe_N_m,
                        Xfe_S_m=Xfe_S_m,
                        phase_tracking_cfg=cfg.phase_tracking,
                        Xmin_troilite_m=Xmin_troilite_m,
                        Xmin_schreibersite_m=Xmin_schreibersite_m,
                        Xmin_cohenite_m=Xmin_cohenite_m,
                        Xmin_graphite_m=Xmin_graphite_m,
                        Xmin_nitride_m=Xmin_nitride_m,
                        Xmin_metal_matrix_m=Xmin_metal_matrix_m,
                        hydrothermal_active=cfg.hydrothermal.active,
                        hydrothermal_cfg=cfg.hydrothermal,
                        deltaIW_m=redox_props !== nothing ? redox_props.deltaIW_m : nothing,
                        xm=xm,
                        ym=ym,
                        coords=coords,
                        rplanet_val=rplanet_val,
                        M_planet_val=M_planet_val,
                        rcore_val=rcore_val,
                        qxD_val=qxD,
                        qyD_val=qyD,
                        xcenter_val=xcenter_val,
                        ycenter_val=ycenter_val,
                    )
                    # interpolate marker properties to basic nodes
                    @inbounds marker_to_basic_nodes!(
                        m,
                        xm[m],
                        ym[m],
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
                    # interpolate marker properties to Vx nodes
                    @inbounds marker_to_vx_nodes!(
                        m,
                        xm[m],
                        ym[m],
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
                    # interpolate marker properties to Vy nodes
                    @inbounds marker_to_vy_nodes!(
                        m,
                        xm[m],
                        ym[m],
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
                    # interpolate marker properties to P nodes
                    @inbounds marker_to_p_nodes!(
                        m,
                        xm[m],
                        ym[m],
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
                end # for m=1:1:marknum
            end

            # ---------------------------------------------------------------------
            # compute physical properties of basic nodes
            # ---------------------------------------------------------------------
            compute_basic_node_properties!(
                ETA0SUM,
                ETASUM,
                GGGSUM,
                SXYSUM,
                COHSUM,
                TENSUM,
                FRISUM,
                WTSUM,
                ETA0,
                ETA,
                GGG,
                SXY0,
                COH,
                TEN,
                FRI,
                YNY,
            )

            # ---------------------------------------------------------------------
            # compute physical properties of Vx nodes
            # ---------------------------------------------------------------------
            compute_vx_node_properties!(
                RHOXSUM, RHOFXSUM, KXSUM, PHIXSUM, RXSUM, WTXSUM, RHOX, RHOFX, KX, PHIX, RX
            )

            # ---------------------------------------------------------------------
            # compute physical properties of Vy nodes
            # ---------------------------------------------------------------------
            compute_vy_node_properties!(
                RHOYSUM, RHOFYSUM, KYSUM, PHIYSUM, RYSUM, WTYSUM, RHOY, RHOFY, KY, PHIY, RY
            )

            # ---------------------------------------------------------------------
            # compute physical properties of P nodes
            # ---------------------------------------------------------------------
            compute_p_node_properties!(
                RHOSUM,
                RHOCPSUM,
                ALPHASUM,
                ALPHAFSUM,
                HRSUM,
                GGGPSUM,
                SXXSUM,
                TKSUM,
                PHISUM,
                WTPSUM,
                RHO,
                RHOCP,
                ALPHA,
                ALPHAF,
                HR,
                GGGP,
                SXX0,
                tk1,
                PHI,
                BETAPHI,
            )

            # ---------------------------------------------------------------------
            # apply thermal boundary conditions for interpolated temperature
            # ---------------------------------------------------------------------
            apply_insulating_boundary_conditions!(tk1)
            # Initialize tk2 from tk1 so thermochemical iterations before
            # the thermal Poisson solve receive valid physical temperatures (T > 0).
            tk2 .= tk1

            # ---------------------------------------------------------------------
            # compute gravity solution
            # compute gravitational acceleration
            # ---------------------------------------------------------------------
            assemble_gravitational_rhs!(RHO, RP; coords=coords)
            SP = F_grav \ RP
            process_gravitational_solution!(SP, FI, gx, gy; coords=coords)

            # ---------------------------------------------------------------------
            # probe increasing computational timestep
            # ---------------------------------------------------------------------
            dt = min(dt*dtcoefup_val, dt_longest_val)
            dt_step_initial = dt
            maxDTcurrent = maximum(abs, DT0)
            @info "\n\n ********** begin timestep $timestep - dt = $dt s **********"

            # ---------------------------------------------------------------------
            # apply surface radiation boundary condition
            # ---------------------------------------------------------------------
            if surface_radiation_val
                tau_LW_val = if cfg.atmosphere.active && atm_state !== nothing
                    atm_state.tau_LW
                else
                    0.0
                end
                T_amb_rad =
                    if cfg.atmosphere.active &&
                        atm_state !== nothing &&
                        atm_state.T_surf_eq > 0.0
                        atm_state.T_surf_eq
                    else
                        T_amb
                    end
                apply_radiative_surface_boundary!(
                    KX,
                    KY,
                    tk1,
                    coords,
                    rplanet_val,
                    xcenter_val,
                    ycenter_val,
                    T_amb_rad;
                    emissivity=emissivity_val,
                    sigma_sb=sigma_sb_val,
                    marker_property_mode=marker_property_mode,
                    phi=phim0_val,
                    kfluid=kfluidm_val[2],
                    tau_LW=tau_LW_val,
                )
            end

            # Snapshot metal inventory at timestep start for idempotent thermochemical iterations
            if coreformation_active_val
                if Xfe_bulk !== nothing && Xfe_bulk_step_start !== nothing
                    if length(Xfe_bulk_step_start) != length(Xfe_bulk)
                        resize!(Xfe_bulk_step_start, length(Xfe_bulk))
                    end
                    copyto!(Xfe_bulk_step_start, Xfe_bulk)
                end
                if Xfem !== nothing && Xfem_step_start !== nothing
                    if length(Xfem_step_start) != length(Xfem)
                        resize!(Xfem_step_start, length(Xfem))
                    end
                    copyto!(Xfem_step_start, Xfem)
                end
            end
            if cfg.metal_partition.active
                if Xfe_H_m !== nothing && Xfe_H_m_step_start !== nothing
                    if length(Xfe_H_m_step_start) != length(Xfe_H_m)
                        resize!(Xfe_H_m_step_start, length(Xfe_H_m))
                    end
                    copyto!(Xfe_H_m_step_start, Xfe_H_m)
                end
                if Xfe_C_m !== nothing && Xfe_C_m_step_start !== nothing
                    if length(Xfe_C_m_step_start) != length(Xfe_C_m)
                        resize!(Xfe_C_m_step_start, length(Xfe_C_m))
                    end
                    copyto!(Xfe_C_m_step_start, Xfe_C_m)
                end
                if Xfe_N_m !== nothing && Xfe_N_m_step_start !== nothing
                    if length(Xfe_N_m_step_start) != length(Xfe_N_m)
                        resize!(Xfe_N_m_step_start, length(Xfe_N_m))
                    end
                    copyto!(Xfe_N_m_step_start, Xfe_N_m)
                end
                if Xfe_S_m !== nothing && Xfe_S_m_step_start !== nothing
                    if length(Xfe_S_m_step_start) != length(Xfe_S_m)
                        resize!(Xfe_S_m_step_start, length(Xfe_S_m))
                    end
                    copyto!(Xfe_S_m_step_start, Xfe_S_m)
                end
            end

            # Snapshot magma inventory at timestep start for idempotent thermochemical iterations
            if magma_active_val
                if Fm !== nothing && Fm_step_start !== nothing
                    if length(Fm_step_start) != length(Fm)
                        resize!(Fm_step_start, length(Fm))
                    end
                    copyto!(Fm_step_start, Fm)
                end
                if F_extract_m !== nothing && F_extract_m_step_start !== nothing
                    if length(F_extract_m_step_start) != length(F_extract_m)
                        resize!(F_extract_m_step_start, length(F_extract_m))
                    end
                    copyto!(F_extract_m_step_start, F_extract_m)
                end
            end

            # ---------------------------------------------------------------------
            # perform thermochemical iterations (outer iteration loop)
            # ---------------------------------------------------------------------
            for titer in 1:1:titermax_val
                # perform thermochemical reaction
                if reaction_active_val
                    perform_thermochemical_reaction!(
                        DMP,
                        DHP,
                        DMPSUM,
                        DHPSUM,
                        WTPSUM,
                        pf,
                        tk2,
                        tm,
                        xm,
                        ym,
                        XWsolidm0,
                        XWsolidm,
                        phim,
                        phinewm,
                        pfm0,
                        marknum,
                        dt,
                        timestep,
                        titer;
                        coords=coords,
                        DQPF=DQPF,
                        DQPFSUM=DQPFSUM,
                        cfg=cfg.reaction,
                    )
                end

                if cfg.refractory.active &&
                    cfg.refractory.kinetics_active &&
                    hcnspo_props !== nothing
                    update_marker_pyrolysis!(
                        tkm,
                        dt,
                        phim,
                        hcnspo_props.X_refr_C_m,
                        hcnspo_props.X_refr_N_m,
                        hcnspo_props.X_refr_H_m,
                        cfg.refractory;
                        xm=xm,
                        ym=ym,
                        coords=coords,
                        DHP=DHP,
                        rhosolid=cfg.materials.rhosolidm,
                        redox_props=redox_props,
                        redox_cfg=cfg.redox,
                        Xfem=Xfem,
                    )
                end

                # -----------------------------------------------------------------
                # perform hydromechanical/plastic iterations (inner iteration loop)
                # -----------------------------------------------------------------

                # save initial viscosity, yielding nodes
                ETA00 .= ETA
                YNY00 .= YNY
                cur_betasolid = timestep == 1 ? 0.0 : betasolid_val
                cur_betafluid = timestep == 1 ? 0.0 : betafluid_val
                if timestep == 1
                    # no elastic compaction during first timestep
                    BETAPHI .= 0.0
                end
                # advance pressure generation inside thermochemical iteration
                pr0 .= pr
                pf0 .= pf

                # perform plastic iterations
                for iplast in 1:1:titermax_val
                    @info("thermochemical iter $titer - hydromechanical iter $iplast")
                    # recompute bulk viscosity at pressure nodes
                    recompute_bulk_viscosity!(ETA, ETAP, ETAPHI, PHI, etaphikoef_val)
                    fill!(S_vent_grid, 0.0)
                    # assemble hydromechanical system of equations
                    if darcy_elim_val
                        L = assemble_hydromechanical_4var_lse!(
                            ETA,
                            ETAP,
                            GGG,
                            GGGP,
                            SXY0,
                            SXX0,
                            RHOX,
                            RHOY,
                            RHOFX,
                            RHOFY,
                            RX,
                            RY,
                            ETAPHI,
                            BETAPHI,
                            PHI,
                            gx,
                            gy,
                            pr0,
                            pf0,
                            DMP,
                            dt,
                            R;
                            coords=coords,
                            betasolid=cur_betasolid,
                            betafluid=cur_betafluid,
                            phimin=phimin_val,
                            phimax=phimax_val,
                            hydrofracture=hydrofracture_val,
                            pr=pr,
                            pf=pf,
                            TEN=TEN,
                            KX=KX,
                            KY=KY,
                            kappa_frac=kappa_frac_val,
                            gamma_frac=gamma_frac_val,
                            k_frac_max=k_frac_max_val,
                            L=L_hydromech,
                            venting=cfg.venting.active,
                            venting_mode=cfg.venting.mode,
                            k_vent=cfg.venting.k_vent,
                            conductance_factor=cfg.venting.conductance_factor,
                            ice_sealing=cfg.venting.ice_sealing,
                            t_freeze=cfg.venting.t_freeze,
                            dt_seal=cfg.venting.dt_seal,
                            k_seal_min_ratio=cfg.venting.k_seal_min_ratio,
                            rplanet=rplanet_val,
                            xcenter=xcenter_val,
                            ycenter=ycenter_val,
                            P_amb=P_amb_eff,
                            venting_species=cfg.venting.species,
                            tk=tk1,
                            eta_fluid_surf=etafluidmm[2],
                            L_sub=cfg.venting.L_sublimation,
                            S_vent_out=S_vent_grid,
                            DQPF=DQPF,
                            fluid_overpressure_coupling=cfg.reaction.fluid_overpressure_coupling,
                            workspace=hydromech_ws,
                        )
                    else
                        L = assemble_hydromechanical_lse!(
                            ETA,
                            ETAP,
                            GGG,
                            GGGP,
                            SXY0,
                            SXX0,
                            RHOX,
                            RHOY,
                            RHOFX,
                            RHOFY,
                            RX,
                            RY,
                            ETAPHI,
                            BETAPHI,
                            PHI,
                            gx,
                            gy,
                            pr0,
                            pf0,
                            DMP,
                            dt,
                            R;
                            coords=coords,
                            betasolid=cur_betasolid,
                            betafluid=cur_betafluid,
                            phimin=phimin_val,
                            phimax=phimax_val,
                            hydrofracture=hydrofracture_val,
                            pr=pr,
                            pf=pf,
                            TEN=TEN,
                            KX=KX,
                            KY=KY,
                            kappa_frac=kappa_frac_val,
                            gamma_frac=gamma_frac_val,
                            k_frac_max=k_frac_max_val,
                            L=L_hydromech,
                            venting=cfg.venting.active,
                            venting_mode=cfg.venting.mode,
                            k_vent=cfg.venting.k_vent,
                            conductance_factor=cfg.venting.conductance_factor,
                            ice_sealing=cfg.venting.ice_sealing,
                            t_freeze=cfg.venting.t_freeze,
                            dt_seal=cfg.venting.dt_seal,
                            k_seal_min_ratio=cfg.venting.k_seal_min_ratio,
                            rplanet=rplanet_val,
                            xcenter=xcenter_val,
                            ycenter=ycenter_val,
                            P_amb=P_amb_eff,
                            venting_species=cfg.venting.species,
                            tk=tk1,
                            eta_fluid_surf=etafluidmm[2],
                            L_sub=cfg.venting.L_sublimation,
                            S_vent_out=S_vent_grid,
                            DQPF=DQPF,
                            fluid_overpressure_coupling=cfg.reaction.fluid_overpressure_coupling,
                            workspace=hydromech_ws,
                        )
                    end
                    # solve hydromechanical system of equations
                    @info "starting hydro-mechanical solver $titer-$iplast"
                    if darcy_elim_val
                        hydromech_ws.pr_presolve .= pr
                        hydromech_ws.pf_presolve .= pf
                    end
                    if cfg.solver.hydromech_solver == :matrix_free
                        op_mf = MatrixFreeStokesDarcyOperator(
                            ETA,
                            ETAP,
                            GGG,
                            GGGP,
                            RHOX,
                            RHOY,
                            RHOFX,
                            RHOFY,
                            RX,
                            RY,
                            ETAPHI,
                            BETAPHI,
                            PHI,
                            gx,
                            gy,
                            dt;
                            coords=coords,
                            betasolid=cur_betasolid,
                            betafluid=cur_betafluid,
                            phimin=phimin_val,
                            phimax=phimax_val,
                        )
                        _, stats = solve_hydromechanical_iterative!(
                            op_mf,
                            R,
                            S;
                            coords=coords,
                            method=cfg.solver.krylov_method,
                            rtol=cfg.solver.krylov_rtol,
                            atol=cfg.solver.krylov_atol,
                            maxiter=cfg.solver.krylov_maxiter,
                            restart=cfg.solver.krylov_restart,
                            preconditioner=cfg.solver.preconditioner,
                            mg_levels=cfg.solver.mg_levels,
                            mg_pre_smooth=cfg.solver.mg_pre_smooth,
                            mg_post_smooth=cfg.solver.mg_post_smooth,
                            mg_smoother=cfg.solver.mg_smoother,
                            mg_omega=cfg.solver.mg_omega,
                        )
                        if !stats.solved
                            error(
                                "Matrix-free Krylov solver $(cfg.solver.krylov_method) failed to converge within $(stats.niter) iterations (status: $(stats.status))",
                            )
                        end
                    elseif cfg.solver.hydromech_solver == :iterative
                        _, stats = solve_hydromechanical_iterative!(
                            L,
                            R,
                            S;
                            coords=coords,
                            method=cfg.solver.krylov_method,
                            rtol=cfg.solver.krylov_rtol,
                            atol=cfg.solver.krylov_atol,
                            maxiter=cfg.solver.krylov_maxiter,
                            restart=cfg.solver.krylov_restart,
                            preconditioner=cfg.solver.preconditioner,
                            mg_levels=cfg.solver.mg_levels,
                            mg_pre_smooth=cfg.solver.mg_pre_smooth,
                            mg_post_smooth=cfg.solver.mg_post_smooth,
                            mg_smoother=cfg.solver.mg_smoother,
                            mg_omega=cfg.solver.mg_omega,
                        )
                        if !stats.solved
                            error(
                                "Iterative Krylov solver $(cfg.solver.krylov_method) failed to converge within $(stats.niter) iterations (status: $(stats.status))",
                            )
                        end
                    elseif use_pardiso_val && pardiso_solver !== nothing
                        L_csc = get_matrix(pardiso_solver, L, :N)
                        current_nnz = nnz(L_csc)
                        if current_nnz != pardiso_last_nnz
                            set_phase!(
                                pardiso_solver, Pardiso.ANALYSIS_NUM_FACT_SOLVE_REFINE
                            )
                            pardiso(pardiso_solver, S, L_csc, R)
                            pardiso_last_nnz = current_nnz
                        else
                            set_phase!(pardiso_solver, Pardiso.NUM_FACT_SOLVE_REFINE)
                            pardiso(pardiso_solver, S, L_csc, R)
                        end
                    else
                        if hydromech_cache === nothing
                            hydromech_prob = LinearProblem(L, R)
                            hydromech_cache = init(
                                hydromech_prob, UMFPACKFactorization(; reuse_symbolic=true)
                            )
                        else
                            hydromech_cache.A = L
                            hydromech_cache.b = R
                        end
                        hydromech_sol = solve!(hydromech_cache)
                        if !LinearSolve.SciMLBase.successful_retcode(hydromech_sol) ||
                            !all(isfinite, hydromech_sol.u)
                            error(
                                "LinearSolve failed with retcode $(hydromech_sol.retcode) or produced non-finite values",
                            )
                        end
                        S .= hydromech_sol.u
                    end
                    @info "finished hydro-mechanical solver $titer-$iplast"
                    # process hydromechanical solution
                    if darcy_elim_val
                        process_hydromechanical_4var_solution!(
                            S, vx, vy, pr, pf; coords=coords
                        )
                        reconstruct_darcy_fluxes!(
                            qxD,
                            qyD,
                            pf,
                            RHOFX,
                            RHOFY,
                            RX,
                            RY,
                            gx,
                            gy,
                            coords;
                            hydrofracture=hydrofracture_val,
                            pr=hydromech_ws.pr_presolve,
                            pf_eff=hydromech_ws.pf_presolve,
                            TEN=TEN,
                            KX=KX,
                            KY=KY,
                            kappa_frac=kappa_frac_val,
                            gamma_frac=gamma_frac_val,
                            k_frac_max=k_frac_max_val,
                        )
                    else
                        process_hydromechanical_solution!(
                            S, vx, vy, pr, qxD, qyD, pf; coords=coords
                        )
                    end

                    # compute Aϕ = Dln[(1-PHI)/PHI]/Dt
                    aphimax = compute_Aϕ!(
                        APHI,
                        ETAPHI,
                        BETAPHI,
                        PHI,
                        pr,
                        pf,
                        pr0,
                        pf0,
                        dt;
                        coords=coords,
                        betasolid=cur_betasolid,
                        phimin=phimin_val,
                        phimax=phimax_val,
                        S_vent=cfg.venting.active ? S_vent_grid : nothing,
                    )

                    # compute fluid velocities
                    compute_fluid_velocities!(
                        PHIX, PHIY, qxD, qyD, vx, vy, vxf, vyf; coords=coords
                    )

                    # adapt timestep for displacement and multiple criteria
                    dt = compute_adaptive_timestep(
                        vx,
                        vy,
                        vxf,
                        vyf,
                        dt,
                        aphimax;
                        coords=coords,
                        dxymax_val=dxymax,
                        dphimax_val=dphimax,
                        dt_ref=dt_step_initial,
                        maxDTcurrent=maxDTcurrent,
                        DTmax_val=DTmax,
                        dt_longest_val=dt_longest_val,
                        max_v_seg=max_v_seg_prev,
                        max_subcycles=cfg.coreformation.max_subcycles,
                        cfl_settling=cfg.coreformation.cfl_settling,
                        DQPF=DQPF,
                        cfl_reaction=cfg.reaction.cfl_reaction,
                        dphi_reaction_max=cfg.reaction.dphi_reaction_max,
                    )

                    # compute stresses, stress changes and strain rate components
                    compute_stress_strainrate!(
                        vx,
                        vy,
                        ETA,
                        GGG,
                        ETAP,
                        GGGP,
                        SXX0,
                        SXY0,
                        EXX,
                        EXY,
                        SXX,
                        SXY,
                        DSXX,
                        DSXY,
                        EII,
                        SII,
                        dt;
                        coords=coords,
                    )

                    # recompute Dln[(1-PHI)/PHI]/Dt
                    _ = compute_Aϕ!(
                        APHI,
                        ETAPHI,
                        BETAPHI,
                        PHI,
                        pr,
                        pf,
                        pr0,
                        pf0,
                        dt;
                        coords=coords,
                        betasolid=cur_betasolid,
                        phimin=phimin_val,
                        phimax=phimax_val,
                        S_vent=cfg.venting.active ? S_vent_grid : nothing,
                    )
                    # symmetrize P node observables
                    symmetrize_p_node_observables!(SXX, APHI, PHI, pr, pf, ps)
                    # consider saving nodal stress changes - RMK: not required
                    # DSXX0 .= DSXX
                    # DSXY0 .= DSXY

                    # nodal adjustment
                    if compute_nodal_adjustment!(
                        ETA,
                        ETA0,
                        ETA5,
                        GGG,
                        SXX,
                        SXY,
                        pr,
                        pf,
                        COH,
                        TEN,
                        FRI,
                        YNY,
                        YNY5,
                        YERRNOD,
                        DSY,
                        dt,
                        iplast,
                    )
                        # exit plastic iterations loop    
                        break
                    else
                        # prepare next pass of plastic iteration 
                        dt = finalize_plastic_iteration_pass!(
                            ETA, ETA5, ETA00, YNY, YNY5, YNY00, YNY_inv_ETA, dt, iplast
                        )
                    end
                end # for iplast=1:1:nplast

                # Refresh venting drainage rate using converged fluid pressure
                if cfg.venting.active
                    apply_venting_surface_boundary!(
                        nothing,
                        nothing,
                        tk1,
                        coords,
                        rplanet_val,
                        xcenter_val,
                        ycenter_val,
                        P_amb_eff;
                        species=cfg.venting.species,
                        k_vent=cfg.venting.k_vent,
                        conductance_factor=cfg.venting.conductance_factor,
                        mode=cfg.venting.mode,
                        hydrofracture=hydrofracture_val,
                        ice_sealing=cfg.venting.ice_sealing,
                        t_freeze=cfg.venting.t_freeze,
                        dt_seal=cfg.venting.dt_seal,
                        k_seal_min_ratio=cfg.venting.k_seal_min_ratio,
                        kappa_frac=kappa_frac_val,
                        gamma_frac=gamma_frac_val,
                        k_frac_max=k_frac_max_val,
                        pr=pr,
                        pf=pf,
                        TEN=TEN,
                        PHI=PHI,
                        phimin=phimin_val,
                        dt=dt,
                        eta_fluid_surf=etafluidmm[2],
                        L_sub=cfg.venting.L_sublimation,
                        S_vent_out=S_vent_grid,
                    )
                end

                # ------------------------------------------------------------------
                # compute shear heating HS in P nodes
                # ------------------------------------------------------------------
                compute_shear_heating!(
                    HS,
                    ETA,
                    SXY,
                    ETAP,
                    SXX,
                    RX,
                    RY,
                    qxD,
                    qyD,
                    PHI,
                    ETAPHI,
                    pr,
                    pf;
                    hydrofracture=hydrofracture_val,
                    TEN=TEN,
                    KX=KX,
                    KY=KY,
                    kappa_frac=kappa_frac_val,
                    gamma_frac=gamma_frac_val,
                    k_frac_max=k_frac_max_val,
                    coords=coords,
                )

                # ------------------------------------------------------------------
                # no pressure changes for the first time step
                # ------------------------------------------------------------------
                if timestep == 1
                    pr0 .= pr
                    pf0 .= pf
                    ps0 .= ps
                end

                # ------------------------------------------------------------------
                # compute adiabatic heating HA in P nodes
                # ------------------------------------------------------------------
                compute_adiabatic_heating!(
                    HA, tk1, ALPHA, ALPHAF, PHI, vx, vy, vxf, vyf, ps, pf; coords=coords
                )

                # ------------------------------------------------------------------
                # compute spherical metric heat source Q_metric
                # ------------------------------------------------------------------
                if spherical_metric_val
                    compute_spherical_metric_heat_source!(
                        Q_metric,
                        tk1,
                        KX,
                        KY,
                        coords;
                        xcenter=xcenter_val,
                        ycenter=ycenter_val,
                        rplanet=rplanet_val,
                        reg_cells=metric_reg_val,
                    )
                end

                # ------------------------------------------------------------------
                # iron core formation segregation
                # ------------------------------------------------------------------
                fill!(Q_seg_grid, 0.0)
                if coreformation_active_val && Xfe_bulk !== nothing && Xfem !== nothing
                    if Xfe_bulk_step_start !== nothing
                        if length(Xfe_bulk) != length(Xfe_bulk_step_start)
                            resize!(Xfe_bulk, length(Xfe_bulk_step_start))
                        end
                        copyto!(Xfe_bulk, Xfe_bulk_step_start)
                    end
                    if Xfem_step_start !== nothing
                        if length(Xfem) != length(Xfem_step_start)
                            resize!(Xfem, length(Xfem_step_start))
                        end
                        copyto!(Xfem, Xfem_step_start)
                    end
                    if cfg.metal_partition.active
                        if Xfe_H_m_step_start !== nothing && Xfe_H_m !== nothing
                            if length(Xfe_H_m) != length(Xfe_H_m_step_start)
                                resize!(Xfe_H_m, length(Xfe_H_m_step_start))
                            end
                            copyto!(Xfe_H_m, Xfe_H_m_step_start)
                        end
                        if Xfe_C_m_step_start !== nothing && Xfe_C_m !== nothing
                            if length(Xfe_C_m) != length(Xfe_C_m_step_start)
                                resize!(Xfe_C_m, length(Xfe_C_m_step_start))
                            end
                            copyto!(Xfe_C_m, Xfe_C_m_step_start)
                        end
                        if Xfe_N_m_step_start !== nothing && Xfe_N_m !== nothing
                            if length(Xfe_N_m) != length(Xfe_N_m_step_start)
                                resize!(Xfe_N_m, length(Xfe_N_m_step_start))
                            end
                            copyto!(Xfe_N_m, Xfe_N_m_step_start)
                        end
                        if Xfe_S_m_step_start !== nothing && Xfe_S_m !== nothing
                            if length(Xfe_S_m) != length(Xfe_S_m_step_start)
                                resize!(Xfe_S_m, length(Xfe_S_m_step_start))
                            end
                            copyto!(Xfe_S_m, Xfe_S_m_step_start)
                        end
                    end
                    seg_res = apply_metal_segregation!(
                        xm,
                        ym,
                        tm,
                        tkm,
                        phim,
                        Xfe_bulk,
                        Xfem,
                        marknum,
                        dt,
                        cfg.coreformation;
                        coords=coords,
                        xcenter=xcenter_val,
                        ycenter=ycenter_val,
                        rplanet=rplanet_val,
                        gx=gx,
                        gy=gy,
                        Q_seg_grid=if cfg.coreformation.segregation_heating
                            Q_seg_grid
                        else
                            nothing
                        end,
                        rho_silicate=rhosolidm[1],
                        eta_silicate=etasolidm[1],
                        ETA=ETA,
                        Fm=Fm,
                        T_solidus_silicate=cfg.melting.T_solidus[1],
                        T_liquidus_silicate=cfg.melting.T_liquidus[1],
                        cfg_partition=cfg.metal_partition,
                        Xfe_H_m=Xfe_H_m,
                        Xfe_C_m=Xfe_C_m,
                        Xfe_N_m=Xfe_N_m,
                        Xfe_S_m=Xfe_S_m,
                        workspace=metal_seg_ws,
                    )
                    max_v_seg_prev = seg_res.max_v_seg
                end

                # ------------------------------------------------------------------
                # solve temperature equation
                # ------------------------------------------------------------------
                if cfg.venting.active && cfg.venting.latent_cooling
                    @. Q_lat_grid = -cfg.venting.L_sublimation * rhofluidm[2] * S_vent_grid
                else
                    fill!(Q_lat_grid, 0.0)
                end

                # ------------------------------------------------------------------
                # silicate melt magma segregation
                # ------------------------------------------------------------------
                if magma_active_val && Fm !== nothing
                    if Fm_step_start !== nothing
                        if length(Fm) != length(Fm_step_start)
                            resize!(Fm, length(Fm_step_start))
                        end
                        copyto!(Fm, Fm_step_start)
                    end
                    if F_extract_m_step_start !== nothing && F_extract_m !== nothing
                        if length(F_extract_m) != length(F_extract_m_step_start)
                            resize!(F_extract_m, length(F_extract_m_step_start))
                        end
                        copyto!(F_extract_m, F_extract_m_step_start)
                    end
                    magma_res = apply_silicate_melt_segregation!(
                        xm,
                        ym,
                        tm,
                        tkm,
                        Fm,
                        marknum,
                        dt,
                        cfg.magma_transport;
                        coords=coords,
                        xcenter=xcenter_val,
                        ycenter=ycenter_val,
                        rplanet=rplanet_val,
                        gx=gx,
                        gy=gy,
                        Q_seg_grid=if cfg.magma_transport.segregation_heating ||
                            cfg.magma_transport.sensible_heat_transport
                            Q_seg_grid
                        else
                            nothing
                        end,
                        Q_lat_grid=if cfg.magma_transport.latent_crystallization
                            Q_lat_grid
                        else
                            nothing
                        end,
                        rho_silicate=rhosolidm[1],
                        rho_melt=cfg.melting.rho_melt,
                        eta_silicate=etasolidm[1],
                        ETA=ETA,
                        T_solidus_silicate=cfg.melting.T_solidus[1],
                        T_liquidus_silicate=cfg.melting.T_liquidus[1],
                        L_melt=cfg.melting.L_melt,
                        F_extract_m=F_extract_m,
                        vx=vx,
                        vy=vy,
                        pr=pr,
                        XH2Om=cfg.volatiles.active ? XH2Om : nothing,
                        XCm=cfg.volatiles.active ? XCm : nothing,
                        XNm=cfg.volatiles.active ? XNm : nothing,
                        XSm=cfg.volatiles.active ? XSm : nothing,
                        phim=phim,
                        cfg_volatiles=cfg.volatiles,
                        workspace=magma_seg_ws,
                    )
                end

                Q_seg_val =
                    if (
                        coreformation_active_val && cfg.coreformation.segregation_heating
                    ) || (
                        magma_active_val && (
                            cfg.magma_transport.segregation_heating ||
                            cfg.magma_transport.sensible_heat_transport
                        )
                    )
                        Q_seg_grid
                    else
                        nothing
                    end
                # assemble thermal system of equations 
                LT = assemble_thermal_lse!(
                    tk1,
                    RHOCP,
                    KX,
                    KY,
                    HR,
                    HA,
                    HS,
                    DHP,
                    RT,
                    dt;
                    coords=coords,
                    LT=LT_thermal,
                    Q_metric=Q_metric,
                    Q_lat=Q_lat_grid,
                    Q_seg=Q_seg_val,
                    workspace=thermal_ws,
                )
                # solve thermal system of equations
                if thermal_cache === nothing
                    thermal_prob = LinearProblem(LT.cscmatrix, RT)
                    thermal_cache = init(
                        thermal_prob, UMFPACKFactorization(; reuse_symbolic=true)
                    )
                else
                    thermal_cache.A = LT.cscmatrix
                    thermal_cache.b = RT
                end
                thermal_sol = solve!(thermal_cache)
                if !LinearSolve.SciMLBase.successful_retcode(thermal_sol) ||
                    !all(isfinite, thermal_sol.u)
                    error(
                        "Thermal solver failed: retcode=$(thermal_sol.retcode), finite=$(all(isfinite, thermal_sol.u))",
                    )
                end
                ST = thermal_sol.u
                # reshape solution vector to 2D array
                tk2 .= reshape(ST, coords.Ny1, coords.Nx1)
                # compute ΔT
                @. DT = tk2 - tk1
                maxDTcurrent = maximum(abs, DT)
                @info "max DT = $maxDTcurrent K"
                # prepare next pass of thermochemical iteration
                dt = finalize_thermochemical_iteration_pass(maxDTcurrent, dt, titer)
                # evaluate iteration outcome
                if compute_thermochemical_iteration_outcome(
                    DMP, pf, pf0, titer; pferrmax=cfg.reaction.pferrmax
                )
                    # exit thermochemical iterations loop
                    break
                end
            end # for titer=1:1:ntiter

            # ---------------------------------------------------------------------
            # advance temperature generation
            # ---------------------------------------------------------------------
            DT0 .= DT

            # ---------------------------------------------------------------------
            # interpolate updated viscoplastic viscosity to markers
            # ---------------------------------------------------------------------
            @threads :static for m in 1:1:marknum
                update_marker_viscosity!(
                    m,
                    xm,
                    ym,
                    tm,
                    tkm,
                    etatotalm,
                    etavpm,
                    YNY,
                    YNY_inv_ETA;
                    coords=coords,
                    Fm=Fm,
                    melting_active=melting_active_val,
                    alpha_eta_val=alpha_eta_val,
                    phi_crit_val=phi_crit_val,
                    eta_melt_val=eta_melt_val,
                )
            end

            # ---------------------------------------------------------------------
            # apply subgrid stress diffusion to markers
            # ---------------------------------------------------------------------
            apply_subgrid_stress_diffusion!(
                xm,
                ym,
                tm,
                inv_gggtotalm,
                sxxm,
                sxym,
                SXX0,
                SXY0,
                DSXX,
                DSXY,
                SXXSUM,
                SXYSUM,
                WTPSUM,
                WTSUM,
                dt,
                marknum;
                coords=coords,
                dsubgrids=dsubgrids,
            )

            # ---------------------------------------------------------------------
            # interpolate DSXX, DSXY to markers
            # ---------------------------------------------------------------------
            update_marker_stress!(xm, ym, sxxm, sxym, DSXX, DSXY, marknum; coords=coords)

            # ---------------------------------------------------------------------
            # apply subgrid temperature diffusion on markers,
            # compute DTsubgrid
            # ---------------------------------------------------------------------
            apply_subgrid_temperature_diffusion!(
                xm,
                ym,
                tm,
                tkm,
                phim,
                tk1,
                DT,
                TKSUM,
                RHOCPSUM,
                dt,
                marknum,
                marker_property_mode;
                coords=coords,
                dsubgridt=dsubgridt,
            )

            # ---------------------------------------------------------------------
            # advance marker temperature, compaction, venting, and volatile drainage
            # ---------------------------------------------------------------------
            XWsolidm0 .= XWsolidm
            phim .= phinewm
            if coreformation_active_val && Xfem !== nothing && Xfem0 !== nothing
                Xfem0 .= Xfem
            end

            marker_results = advance_marker_thermo_porosity_venting!(
                xm,
                ym,
                tm,
                tkm,
                phim,
                DT,
                tk2,
                APHI,
                dt,
                timestep,
                marknum;
                coords=coords,
                phimin=phimin_val,
                phimax=phimax_val,
                venting=cfg.venting.active,
                S_vent_grid=S_vent_grid,
                rhofluidcur=rhofluidm[2],
                ret_cfg=cfg.retention,
                XH2Om=cfg.volatiles.active ? XH2Om : nothing,
                XCm=cfg.volatiles.active ? XCm : nothing,
                XNm=cfg.volatiles.active ? XNm : nothing,
                XSm=cfg.volatiles.active ? XSm : nothing,
                rhosolid=cfg.materials.rhosolidm,
                Fm=Fm,
                redox_props=redox_props,
            )

            delta_m_vent = marker_results.delta_m_vent
            delta_m_vent_3d = 0.0
            vented_vols = marker_results.vented_vols

            if cfg.venting.active
                # Area flux geometric scaling: 4πR² / 2πR = 2 * R_planet
                L_3D_equiv = 2.0 * rplanet_val
                delta_m_vent_3d = delta_m_vent * L_3D_equiv
                M_vent_total += delta_m_vent_3d

                if vented_vols !== nothing
                    M_vent_H2O_total += vented_vols.M_vent_H2O * L_3D_equiv
                    M_vent_C_total += vented_vols.M_vent_C * L_3D_equiv
                    M_vent_N_total += vented_vols.M_vent_N * L_3D_equiv
                    M_vent_S_total += vented_vols.M_vent_S * L_3D_equiv
                end
            end

            if cfg.atmosphere.active && atm_state !== nothing
                L_3D_equiv = 2.0 * rplanet_val
                p_surf_val = max(atm_state.P_surf, P_amb_eff)
                T_surf_val = atm_state.T_surf_eq > 0.0 ? atm_state.T_surf_eq : T_amb
                vent_rates = compute_surface_venting_rates(
                    cfg,
                    delta_m_vent_3d,
                    vented_vols,
                    L_3D_equiv,
                    dt,
                    p_surf_val,
                    T_surf_val,
                    redox_props,
                    marknum,
                    tm,
                    xm,
                    ym,
                    rplanet_val,
                    xcenter_val,
                    ycenter_val,
                )

                c_s_disk = compute_sound_speed(T_amb)
                rho_disk_val = if (cfg.disk.enabled && c_s_disk > 0.0)
                    max(0.0, (1.0 - w_disp) * cfg.disk.p_amb_disk) / (c_s_disk^2)
                else
                    0.0
                end
                a_orb_val = cfg.disk.orbital_distance_au * AU_METERS
                M_star_val = cfg.disk.stellar_mass_msun * M_SUN_KG
                R_exo_val = max(cfg.escape.R_exobase, rplanet_val)
                T_int_val = compute_mean_surface_temperature(
                    tk1, coords, rplanet_val, xcenter_val, ycenter_val; T_default=T_amb
                )

                degas_rates =
                    if cfg.magma_degassing.active && XH2Om !== nothing && Fm !== nothing
                        p_surf_mo = atm_state.P_surf > 0.0 ? atm_state.P_surf : P_amb_eff
                        fO2_diw_mo = compute_surface_mean_delta_iw(
                            redox_props,
                            tm,
                            xm,
                            ym,
                            marknum,
                            xcenter_val,
                            ycenter_val,
                            rplanet_val,
                            cfg.volatiles.fO2_delta_IW,
                        )
                        v_m = (coords.xsize * coords.ysize) / max(1, marknum)
                        Fm_prev = Fm_step_start !== nothing ? Fm_step_start : Fm

                        if cfg.magma_degassing.mode === :dynamic_flux
                            degas_magma_ocean_markers!(
                                xm,
                                ym,
                                tm,
                                tkm,
                                Fm,
                                Fm_prev,
                                XH2Om,
                                XCm,
                                XNm,
                                XSm,
                                marknum,
                                dt,
                                p_surf_mo,
                                rplanet_val,
                                cfg.magma_degassing;
                                rho_solid=cfg.materials.rhosolidm[1],
                                marker_volume=v_m,
                                delta_IW=fO2_diw_mo,
                                retention_cfg=cfg.retention,
                            )
                        else
                            # Equilibrium partitioning mode across molten magma ocean
                            m_melt_tot = 0.0
                            m_H_melt = 0.0
                            m_C_melt = 0.0
                            m_N_melt = 0.0
                            m_S_melt = 0.0
                            m_marker =
                                cfg.materials.rhosolidm[1] * v_m * (2.0 * rplanet_val)
                            for m in 1:marknum
                                if tm[m] < 3 &&
                                    (
                                        (
                                            (xm[m] - xcenter_val)^2 +
                                            (ym[m] - ycenter_val)^2
                                        ) <= rplanet_val^2
                                    ) &&
                                    Fm[m] >= cfg.magma_degassing.F_melt_threshold
                                    m_melt_tot += Fm[m] * m_marker
                                    m_H_melt +=
                                        (XH2Om[m] * 0.01) * (2.01588 / 18.01528) * m_marker
                                    m_C_melt += (XCm[m] * 1.0e-6) * m_marker
                                    m_N_melt += (XNm[m] * 1.0e-6) * m_marker
                                    m_S_melt += (XSm[m] * 1.0e-6) * m_marker
                                end
                            end

                            # Atmospheric elemental inventories
                            m_H_atm =
                                get(atm_state.M_atm, :H2, 0.0) * 1.0 +
                                get(atm_state.M_atm, :H2O, 0.0) * (2.01588 / 18.01528) +
                                get(atm_state.M_atm, :CH4, 0.0) * (4.03176 / 16.04246) +
                                get(atm_state.M_atm, :NH3, 0.0) * (3.02382 / 17.03052) +
                                get(atm_state.M_atm, :H2S, 0.0) * (2.01588 / 34.08088)

                            m_C_atm =
                                get(atm_state.M_atm, :CO, 0.0) * (12.011 / 28.0101) +
                                get(atm_state.M_atm, :CO2, 0.0) * (12.011 / 44.0095) +
                                get(atm_state.M_atm, :CH4, 0.0) * (12.011 / 16.04246)

                            m_N_atm =
                                get(atm_state.M_atm, :N2, 0.0) * 1.0 +
                                get(atm_state.M_atm, :NH3, 0.0) * (14.007 / 17.03052)

                            m_S_atm =
                                get(atm_state.M_atm, :H2S, 0.0) * (32.060 / 34.08088) +
                                get(atm_state.M_atm, :SO2, 0.0) * (32.060 / 64.066) +
                                get(atm_state.M_atm, :S2, 0.0) * 1.0

                            m_H_tot = m_H_melt + m_H_atm
                            m_C_tot = m_C_melt + m_C_atm
                            m_N_tot = m_N_melt + m_N_atm
                            m_S_tot = m_S_melt + m_S_atm

                            if m_melt_tot > 0.0 &&
                                (m_H_tot + m_C_tot + m_N_tot + m_S_tot) > 0.0
                                g_surf =
                                    GRAVITATIONAL_CONSTANT * M_planet_val / (rplanet_val^2)
                                sol_eq = solve_magma_ocean_volatile_partitioning(
                                    m_melt_tot,
                                    m_H_tot,
                                    m_C_tot,
                                    m_N_tot,
                                    m_S_tot,
                                    rplanet_val,
                                    g_surf,
                                    T_int_val,
                                    fO2_diw_mo,
                                )

                                # Deplete molten markers according to residual melt volatile concentration
                                new_XH2O_wtpct =
                                    (sol_eq.M_melt_H * (18.01528 / 2.01588) / m_melt_tot) *
                                    100.0
                                new_XC_ppm = (sol_eq.M_melt_C / m_melt_tot) * 1.0e6
                                new_XN_ppm = (sol_eq.M_melt_N / m_melt_tot) * 1.0e6
                                new_XS_ppm = (sol_eq.M_melt_S / m_melt_tot) * 1.0e6

                                for m in 1:marknum
                                    if tm[m] < 3 &&
                                        (
                                            (
                                                (xm[m] - xcenter_val)^2 +
                                                (ym[m] - ycenter_val)^2
                                            ) <= rplanet_val^2
                                        ) &&
                                        Fm[m] >= cfg.magma_degassing.F_melt_threshold
                                        XH2Om[m] = new_XH2O_wtpct
                                        XCm[m] = new_XC_ppm
                                        XNm[m] = new_XN_ppm
                                        XSm[m] = new_XS_ppm
                                    end
                                end

                                rates = Dict{Symbol,Float64}()
                                for (sp, m_atm_eq) in sol_eq.M_atm_i
                                    m_atm_cur = get(atm_state.M_atm, sp, 0.0)
                                    rates[sp] = max(0.0, m_atm_eq - m_atm_cur) / dt
                                end
                                rates
                            else
                                nothing
                            end
                        end
                    else
                        nothing
                    end

                evolve_coupled_atmosphere_step!(
                    atm_state,
                    vent_rates,
                    dt,
                    M_planet_val,
                    rplanet_val,
                    T_amb,
                    cfg.atmosphere;
                    rho_disk=rho_disk_val,
                    c_s=c_s_disk,
                    M_star=M_star_val,
                    a_orb=a_orb_val,
                    T_int=T_int_val,
                    T_exobase=cfg.escape.T_exobase,
                    R_exobase=R_exo_val,
                    hydrodynamic=cfg.escape.hydrodynamic,
                    gamma=cfg.escape.gamma,
                    escape_active=cfg.escape.active,
                    escape_cfg=cfg.escape,
                    sim_time_s=timesum,
                    degas_rates=degas_rates,
                )

                M_atm_total = sum(values(atm_state.M_atm))
                M_escaped_total = sum(values(atm_state.M_escaped))
                if M_atm_species !== nothing
                    for (sp, val) in atm_state.M_atm
                        M_atm_species[sp] = val
                    end
                end
                if M_escaped_species !== nothing
                    for (sp, val) in atm_state.M_escaped
                        M_escaped_species[sp] = val
                    end
                end
            elseif cfg.escape.active
                L_3D_equiv = 2.0 * rplanet_val
                R_exo_val = max(cfg.escape.R_exobase, rplanet_val)
                T_surf_esc = compute_mean_surface_temperature(
                    tk1, coords, rplanet_val, xcenter_val, ycenter_val; T_default=T_amb
                )

                vent_rates_esc = compute_surface_venting_rates(
                    cfg,
                    delta_m_vent_3d,
                    vented_vols,
                    L_3D_equiv,
                    dt,
                    P_amb_eff,
                    T_surf_esc,
                    redox_props,
                    marknum,
                    tm,
                    xm,
                    ym,
                    rplanet_val,
                    xcenter_val,
                    ycenter_val,
                )

                if cfg.escape.multi_species &&
                    M_atm_species !== nothing &&
                    M_escaped_species !== nothing
                    for sp in cfg.escape.species_list
                        m_sp = get_species_molecular_mass(sp)
                        v_rate_sp = get(vent_rates_esc, sp, 0.0)
                        prev_sp = get(M_atm_species, sp, 0.0)
                        esc_sp = evolve_atmospheric_species_inventory(
                            prev_sp,
                            v_rate_sp,
                            dt,
                            M_planet_val,
                            rplanet_val,
                            cfg.escape.T_exobase,
                            m_sp;
                            R_exobase=R_exo_val,
                            gamma=cfg.escape.gamma,
                            hydrodynamic=cfg.escape.hydrodynamic,
                        )
                        M_atm_species[sp] = esc_sp.M_atm
                        M_escaped_species[sp] =
                            get(M_escaped_species, sp, 0.0) + esc_sp.M_escaped_step
                    end
                    M_atm_total = sum(values(M_atm_species))
                    M_escaped_total = sum(values(M_escaped_species))
                else
                    M_vent_rate_eff = get(vent_rates_esc, cfg.escape.species, 0.0)
                    esc_res = evolve_atmospheric_species_inventory(
                        M_atm_total,
                        M_vent_rate_eff,
                        dt,
                        M_planet_val,
                        rplanet_val,
                        cfg.escape.T_exobase,
                        get_species_molecular_mass(cfg.escape.species);
                        R_exobase=R_exo_val,
                        gamma=cfg.escape.gamma,
                        hydrodynamic=cfg.escape.hydrodynamic,
                    )
                    M_atm_total = esc_res.M_atm
                    M_escaped_total += esc_res.M_escaped_step
                end
            end
            phinewm .= phim

            # ---------------------------------------------------------------------
            # interpolate melt composition from markers to P nodes
            # --------------------------------------------------------------------- 
            update_p_nodes_melt_composition!(
                xm, ym, XWsolidm0, XWS, XWSSUM, WTPSUM, marknum; coords=coords
            )

            # ---------------------------------------------------------------------
            # compute velocity in P nodes,
            # compute fluid velocity in P nodes including boundary conditions
            # ---------------------------------------------------------------------
            compute_velocities!(vx, vy, vxf, vyf, vxp, vyp, vxpf, vypf; coords=coords)

            # ---------------------------------------------------------------------
            # compute rotation rate in basic nodes
            # ---------------------------------------------------------------------
            compute_rotation_rate!(vx, vy, wyx; coords=coords)

            # ---------------------------------------------------------------------
            # move markers with RK4
            # ---------------------------------------------------------------------
            move_markers_rk4!(
                xm,
                ym,
                tm,
                tkm,
                phim,
                sxym,
                sxxm,
                vx,
                vy,
                vxf,
                vyf,
                wyx,
                tk2,
                marknum,
                dt,
                marker_property_mode;
                coords=coords,
            )

            # ---------------------------------------------------------------------
            # backtrack P nodes: Ptotal with RK4,
            # backtrack P nodes: Pfluid with RK4
            # ---------------------------------------------------------------------
            backtrace_pressures_rk4!(
                pr, pr0, ps, ps0, pf, pf0, vx, vy, vxf, vyf, dt; coords=coords
            )

            # ---------------------------------------------------------------------
            # replenish sparse areas with additional markers
            # ---------------------------------------------------------------------
            marknum = replenish_markers!(
                xm,
                ym,
                tm,
                tkm,
                phim,
                sxxm,
                sxym,
                etavpm,
                phinewm,
                pfm0,
                XWsolidm,
                XWsolidm0,
                rhototalm,
                rhocptotalm,
                etatotalm,
                hrtotalm,
                ktotalm,
                inv_gggtotalm,
                fricttotalm,
                cohestotalm,
                tenstotalm,
                rhofluidcur,
                alphasolidcur,
                alphafluidcur,
                tkm_rhocptotalm,
                etafluidcur_inv_kphim,
                mdis,
                mnum;
                Fm=Fm,
                randomized=random_markers,
                coords=coords,
                Xfem=Xfem,
                Xfem0=Xfem0,
                Xfe_bulk=Xfe_bulk,
                XH2Om=XH2Om,
                XCm=XCm,
                XNm=XNm,
                XSm=XSm,
                Xfe_H_m=Xfe_H_m,
                Xfe_C_m=Xfe_C_m,
                Xfe_N_m=Xfe_N_m,
                Xfe_S_m=Xfe_S_m,
                Xmin_troilite_m=Xmin_troilite_m,
                Xmin_schreibersite_m=Xmin_schreibersite_m,
                Xmin_cohenite_m=Xmin_cohenite_m,
                Xmin_graphite_m=Xmin_graphite_m,
                Xmin_nitride_m=Xmin_nitride_m,
                Xmin_metal_matrix_m=Xmin_metal_matrix_m,
                t_accreted=t_accreted,
                hcnspo_props=hcnspo_props,
                F_extract_m=F_extract_m,
            )
            if t_accreted !== nothing && length(t_accreted) != marknum
                resize!(t_accreted, marknum)
            end
            if hcnspo_props !== nothing
                for prop in values(hcnspo_props)
                    if length(prop) != marknum
                        resize!(prop, marknum)
                    end
                end
            end
            if coreformation_active_val
                if Xfe_bulk_step_start !== nothing && length(Xfe_bulk_step_start) != marknum
                    resize!(Xfe_bulk_step_start, marknum)
                end
                if Xfem_step_start !== nothing && length(Xfem_step_start) != marknum
                    resize!(Xfem_step_start, marknum)
                end
            end
            if magma_active_val
                if F_extract_m_step_start !== nothing &&
                    length(F_extract_m_step_start) != marknum
                    resize!(F_extract_m_step_start, marknum)
                end
                if Fm_step_start !== nothing && length(Fm_step_start) != marknum
                    resize!(Fm_step_start, marknum)
                end
            end
            if cfg.metal_partition.active
                if Xfe_H_m_step_start !== nothing && length(Xfe_H_m_step_start) != marknum
                    resize!(Xfe_H_m_step_start, marknum)
                end
                if Xfe_C_m_step_start !== nothing && length(Xfe_C_m_step_start) != marknum
                    resize!(Xfe_C_m_step_start, marknum)
                end
                if Xfe_N_m_step_start !== nothing && length(Xfe_N_m_step_start) != marknum
                    resize!(Xfe_N_m_step_start, marknum)
                end
                if Xfe_S_m_step_start !== nothing && length(Xfe_S_m_step_start) != marknum
                    resize!(Xfe_S_m_step_start, marknum)
                end
            end

            # ---------------------------------------------------------------------
            # update timesum
            # ---------------------------------------------------------------------
            timesum += dt
            timestep_end = now()

            # ---------------------------------------------------------------------
            # save data evaluation and core budget update
            # ---------------------------------------------------------------------
            should_save_snapshot = if cfg.output.mode in (:snapshots, :both)
                timestep % savematstep_val == 0
            elseif cfg.output.mode == :telemetry
                cfg.output.save_final && (timestep == n_steps_val)
            else
                false
            end

            need_telemetry = (
                telemetry_io !== nothing &&
                (timestep % cfg.output.telemetrystep == 0 || timestep == n_steps_val)
            )

            if (cfg.metal_partition.active && Xfe_bulk !== nothing) &&
                (need_telemetry || should_save_snapshot || cfg.hydrothermal.active)
                core_budgets = compute_core_volatile_budgets(
                    xm,
                    ym,
                    tm,
                    Xfe_bulk,
                    Xfe_H_m,
                    Xfe_C_m,
                    Xfe_N_m,
                    Xfe_S_m,
                    marknum;
                    xcenter=xcenter_val,
                    ycenter=ycenter_val,
                    rplanet=rplanet_val,
                    rho_metal=cfg.coreformation.rho_metal,
                    core_radius_fraction=cfg.metal_partition.core_radius_fraction,
                    phi_core_threshold=cfg.metal_partition.phi_core_threshold,
                )
                if (core_budgets !== nothing && core_budgets.M_core_metal > 0.0)
                    rcore_val =
                        (
                            3.0 * core_budgets.M_core_metal /
                            (4.0 * π * cfg.coreformation.rho_metal)
                        )^(1.0 / 3.0)
                else
                    rcore_val = 0.0
                end
            end

            # ---------------------------------------------------------------------
            # streaming telemetry record
            # ---------------------------------------------------------------------
            if need_telemetry
                core_radius_current = rcore_val
                stream_telemetry_row!(
                    telemetry_io,
                    timestep,
                    s_to_Ma(timesum),
                    dt / cfg.time.yearlength,
                    rplanet_val,
                    core_radius_current,
                    maximum(tk2),
                    sum(tk2) / length(tk2),
                    maximum(PHI),
                    sum(PHI) / length(PHI),
                    M_vent_H2O_total + M_vent_C_total + M_vent_N_total + M_vent_S_total,
                    M_vent_H2O_total,
                    M_atm_total,
                    M_escaped_total,
                    Fm !== nothing ? maximum(Fm) : 0.0,
                    Fm !== nothing ? (sum(Fm) / length(Fm)) : 0.0,
                )
            end

            # ---------------------------------------------------------------------
            # save data for analysis and visualization
            # ---------------------------------------------------------------------
            if should_save_snapshot
                if cfg.phase_tracking.active &&
                    cfg.phase_tracking.track_regional_modes &&
                    Xfe_bulk !== nothing
                    regional_mineral_modes = compute_regional_mineral_modes(
                        xm,
                        ym,
                        tm,
                        tkm,
                        Xfe_bulk,
                        Xfe_S_m,
                        Xfe_C_m,
                        Xfe_N_m,
                        marknum;
                        cfg=cfg.phase_tracking,
                        xcenter=xcenter_val,
                        ycenter=ycenter_val,
                        rplanet=rplanet_val,
                        rho_metal=cfg.coreformation.rho_metal,
                    )
                end
                save_state(
                    output_path,
                    timestep,
                    dt,
                    timesum,
                    marknum,
                    ETA,
                    ETA0,
                    GGG,
                    EXY,
                    SXY,
                    SXY0,
                    wyx,
                    COH,
                    TEN,
                    FRI,
                    YNY,
                    RHOX,
                    RHOFX,
                    KX,
                    PHIX,
                    vx,
                    vxf,
                    RX,
                    qxD,
                    gx,
                    RHOY,
                    RHOFY,
                    KY,
                    PHIY,
                    vy,
                    vyf,
                    RY,
                    qyD,
                    gy,
                    RHO,
                    RHOCP,
                    ALPHA,
                    ALPHAF,
                    HR,
                    HA,
                    HS,
                    ETAP,
                    GGGP,
                    EXX,
                    SXX,
                    SXX0,
                    tk1,
                    tk2,
                    vxp,
                    vyp,
                    vxpf,
                    vypf,
                    pr,
                    pf,
                    ps,
                    pr0,
                    pf0,
                    ps0,
                    ETAPHI,
                    BETAPHI,
                    PHI,
                    APHI,
                    FI,
                    ETA5,
                    ETA00,
                    YNY5,
                    YNY00,
                    YNY_inv_ETA,
                    DSXY,
                    EII,
                    SII,
                    DSXX,
                    DMP,
                    DHP,
                    DQPF,
                    XWS,
                    XWsolidm0,
                    xm,
                    ym,
                    tm,
                    tkm,
                    sxxm,
                    sxym,
                    etavpm,
                    phim,
                    rhototalm,
                    rhocptotalm,
                    etatotalm,
                    hrtotalm,
                    ktotalm,
                    tkm_rhocptotalm,
                    etafluidcur_inv_kphim,
                    inv_gggtotalm,
                    fricttotalm,
                    cohestotalm,
                    tenstotalm,
                    rhofluidcur,
                    alphasolidcur,
                    alphafluidcur;
                    coords=coords,
                    phim0_val=phim0_val,
                    M_vent_total=M_vent_total,
                    M_vent_H2O_total=M_vent_H2O_total,
                    M_vent_C_total=M_vent_C_total,
                    M_vent_N_total=M_vent_N_total,
                    M_vent_S_total=M_vent_S_total,
                    M_atm_total=M_atm_total,
                    M_escaped_total=M_escaped_total,
                    P_amb=P_amb_eff,
                    S_vent=S_vent_grid,
                    Xfem=Xfem,
                    Xfem0=Xfem0,
                    Xfe_bulk=Xfe_bulk,
                    M_atm_species=M_atm_species,
                    M_escaped_species=M_escaped_species,
                    XH2Om=XH2Om,
                    XCm=XCm,
                    XNm=XNm,
                    XSm=XSm,
                    Xfe_H_m=Xfe_H_m,
                    Xfe_C_m=Xfe_C_m,
                    Xfe_N_m=Xfe_N_m,
                    Xfe_S_m=Xfe_S_m,
                    core_budgets=core_budgets,
                    Xmin_troilite_m=Xmin_troilite_m,
                    Xmin_schreibersite_m=Xmin_schreibersite_m,
                    Xmin_cohenite_m=Xmin_cohenite_m,
                    Xmin_graphite_m=Xmin_graphite_m,
                    Xmin_nitride_m=Xmin_nitride_m,
                    Xmin_metal_matrix_m=Xmin_metal_matrix_m,
                    regional_mineral_modes=regional_mineral_modes,
                    DT0=DT0,
                    rplanet=rplanet_val,
                    rcore=rcore_val,
                    t_accreted=t_accreted,
                    M_accreted_total=cfg.accretion.active ? M_accreted_total : nothing,
                    M_planet_val=M_planet_val,
                    telescope_level=telescope_level,
                    hcnspo_props=hcnspo_props,
                    redox_props=redox_props,
                    atm_state=atm_state,
                )
            end
            # ---------------------------------------------------------------------
            #  save old stresses - RMK: not used anywhere in code
            # ---------------------------------------------------------------------
            #  sxxm00 = sxxm 
            #  sxym00 = sxym    

            # ---------------------------------------------------------------------
            # update progress indicators
            # ---------------------------------------------------------------------
            maxT = maximum(tk2)
            @info "timestep $timestep computed in $(
            Dates.canonicalize(
                Dates.CompoundPeriod(timestep_end-timestep_begin)
            )
        )"
            @info "total time = $(s_to_Ma(timesum)) Ma"
            @info "markers in use = $marknum"
            @info "max T = $maxT K"
            next!(p; showvalues=generate_showvalues(timestep, marknum, maxT, dt, timesum))

            # ---------------------------------------------------------------------
            # finish timestep
            # ---------------------------------------------------------------------
            if timesum > endtime_val
                break
            end
        end # for timestep = startstep:1:n_steps
    finally
        if telemetry_io !== nothing
            close(telemetry_io)
        end
        if use_pardiso_val && pardiso_solver !== nothing
            set_phase!(pardiso_solver, Pardiso.RELEASE_ALL)
            pardiso(pardiso_solver)
        end
    end
end # function simulation loop

"""
Simulation loop overload for path or configuration file input.

$(SIGNATURES)

# Arguments
- `path_or_dir`: Path to `.toml` configuration file, or output directory path.
- `output_path`: Optional output path override.
"""
function simulation_loop(path_or_dir::String; output_path::String="")
    if endswith(path_or_dir, ".toml")
        cfg = load_config(path_or_dir)
        actual_output = isempty(output_path) ? cfg.output.output_dir : output_path
        return simulation_loop(cfg; output_path=actual_output)
    else
        actual_output = isempty(output_path) ? path_or_dir : output_path
        return simulation_loop(default_config(); output_path=actual_output)
    end
end

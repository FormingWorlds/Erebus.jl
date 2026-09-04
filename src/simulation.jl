
"""
Convert seconds to Ma (millions of years).

$(SIGNATURES)

# Details

    - s: period in seconds

# Returns

    - Ma: period in millions of years
"""
function s_to_Ma(s)
    return s / (yearlength * 1e6)
end

"""
Set up and initialize dynamic simulation parameters.

$(SIGNATURES)

# Details

    - nothing

# Returns

    - timestep: simulation starting time step count
    - dt: simulation initial computational time step [s]
    - timesum: simulation starting time [s]
    - marknum: initial number of markers
    - hrsolidm: initial radiogenic heat production solid phase
    - hrfluidm: initial radiogenic heat production fluid phase
    - YERRNOD: vector of summed yielding errors of nodes over plastic iterations
"""
function setup_dynamic_simulation_parameters(
    cfg::SimulationConfig=default_config(); coords::Union{Nothing,GridCoordinates}=nothing
)
    # timestep counter (current), init to startstep
    timestep::Int64 = cfg.time.start_step
    # computational timestep (current), init to dt_longest [s]
    dt::Float64 = cfg.time.dt_longest
    # time sum (current), init to starttime [s]
    timesum::Float64 = cfg.time.start_time
    # current number of markers, init to startmarknum
    marknum::Int64 = coords === nothing ? start_marknum : coords.start_marknum
    # radiogenic heat production solid phase
    hrsolidm::SVector{3,Float64} = start_hrsolidm
    # radiogenic heat production fluid phase
    hrfluidm::SVector{3,Float64} = start_hrfluidm
    # nodes yielding error vector of plastic iterations
    YERRNOD::Vector{Float64} = zeros(Float64, cfg.solver.nplast)
    return timestep, dt, timesum, marknum, hrsolidm, hrfluidm, YERRNOD
end # function setup_dynamic_simulation_parameters()

"""
Save simulation state to JLD2 output file named after current timestep.

$(SIGNATURES)

# Details

    - output_path: absolute path to output directory
    - timestep: current time step number
    - dt: time step
    - timesum: total simulation time
    - marknum: number of markers
    - ETA... : simulation state variables

# Returns

    - nothing
"""
function save_state(
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
    coords::Union{Nothing,GridCoordinates}=nothing,
)
    # @timeit to "save_state" begin
    fid = output_path * "output_" * lpad(timestep, 5, "0") * ".jld2"
    Nx_val = coords === nothing ? Nx : coords.Nx
    Ny_val = coords === nothing ? Ny : coords.Ny
    Nx1_val = coords === nothing ? Nx1 : coords.Nx1
    Ny1_val = coords === nothing ? Ny1 : coords.Ny1
    Nxm_val = coords === nothing ? Nxm : coords.Nxm
    Nym_val = coords === nothing ? Nym : coords.Nym
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    dxm_val = coords === nothing ? dxm : coords.dxm
    dym_val = coords === nothing ? dym : coords.dym
    x_val = coords === nothing ? x : coords.x
    y_val = coords === nothing ? y : coords.y
    xvx_val = coords === nothing ? xvx : coords.xvx
    yvx_val = coords === nothing ? yvx : coords.yvx
    xvy_val = coords === nothing ? xvy : coords.xvy
    yvy_val = coords === nothing ? yvy : coords.yvy
    xp_val = coords === nothing ? xp : coords.xp
    yp_val = coords === nothing ? yp : coords.yp
    xxm_val = coords === nothing ? xxm : coords.xxm
    yym_val = coords === nothing ? yym : coords.yym
    xsize_val = coords === nothing ? xsize : coords.xsize
    ysize_val = coords === nothing ? ysize : coords.ysize
    xcenter_val = coords === nothing ? xcenter : coords.xcenter
    ycenter_val = coords === nothing ? ycenter : coords.ycenter
    jldsave(
        fid;
        timestep,
        dt,
        Δtreaction,
        reaction_rate_coeff_mode,
        marker_property_mode,
        timesum,
        marknum,
        phim0,
        XWsolidm_init,
        ratio_al,
        t_half_al,
        dsubgrids,
        dsubgridt,
        hr_al,
        hr_fe,
        rplanet,
        rcrust,
        psurface,
        xsize=xsize_val,
        ysize=ysize_val,
        xcenter=xcenter_val,
        ycenter=ycenter_val,
        Nx=Nx_val,
        Ny=Ny_val,
        Nx1=Nx1_val,
        Ny1=Ny1_val,
        Nxm=Nxm_val,
        Nym=Nym_val,
        dx=dx_val,
        dy=dy_val,
        dxm=dxm_val,
        dym=dym_val,
        x=x_val,
        y=y_val,
        xvx=xvx_val,
        yvx=yvx_val,
        xvy=xvy_val,
        yvy=yvy_val,
        xp=xp_val,
        yp=yp_val,
        xxm=xxm_val,
        yym=yym_val,
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
        alphafluidcur,
    )
    # end # @timeit to "save_state"
    return nothing
end

"""
Load simulation state from a JLD2 checkpoint archive.

$(SIGNATURES)

# Details

    - checkpoint_path: absolute or relative path to JLD2 checkpoint archive

# Returns

    - checkpoint_data: dictionary containing saved state arrays and progression parameters
"""
function load_state(checkpoint_path::AbstractString)
    isfile(checkpoint_path) ||
        throw(ArgumentError("Checkpoint file does not exist: $checkpoint_path"))
    return JLD2.load(checkpoint_path)
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
    dt_longest_val = cfg.time.dt_longest
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
    rplanet_val = cfg.geometry.rplanet
    rcrust_val = cfg.geometry.rcrust
    xcenter_val = cfg.geometry.xcenter
    ycenter_val = cfg.geometry.ycenter
    psurface_val = cfg.geometry.psurface

    nthreads = Threads.nthreads()
    use_threading = nthreads > 1
    num_buffers = nthreads
    thread_buffers =
        use_threading ? allocate_thread_interpolation_buffers(num_buffers, coords) : nothing

    @info "Simulation layout" coords.Nx coords.Ny coords.xsize coords.dx coords.dy coords.ysize rplanet_val rcrust_val marknum nthreads
    @info(
        "Parameters",
        random_markers,
        marker_property_mode,
        hr_al,
        hr_fe,
        reaction_active,
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
        XWsolidm_init,
        etaphikoef,
        αη,
        tmsolidphase,
        tmfluidphase,
        phim0,
        phimin,
        phimax,
        ΔHWD,
        ΔSWD,
        ΔVWD,
        Δtreaction,
        pfcoeff,
        pferrmax,
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

    # -------------------------------------------------------------------------
    # set up markers and state (from checkpoint or fresh definition)
    # -------------------------------------------------------------------------
    mdis, mnum = setup_marker_geometry_helpers(coords)
    is_restart = !isempty(restart_from)
    if is_restart
        ckpt = load_state(restart_from)
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
        start_step_val = ckpt["timestep"] + 1
        dt = ckpt["dt"]
        timesum = ckpt["timesum"]
        marknum = ckpt["marknum"]
        n_steps_val = cfg.time.n_steps
        if n_steps_val < start_step_val
            @warn "Restart checkpoint timestep ($(ckpt["timestep"])) >= target n_steps ($n_steps_val). No timesteps will be executed."
        end

        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0) = setup_marker_properties(
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
        @info "Resumed simulation from checkpoint: $restart_from at timestep $(start_step_val-1) (running to $n_steps_val)"
    else
        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0) = setup_marker_properties(
            marknum, coords
        )
        (rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm, tkm_rhocptotalm, etafluidcur_inv_kphim, inv_gggtotalm, fricttotalm, cohestotalm, tenstotalm, rhofluidcur, alphasolidcur, alphafluidcur) = setup_marker_properties_helpers(
            marknum
        )
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
        )
        # copy thermodynamic marker properties to next generation for initial setup
        XWsolidm .= XWsolidm0
        phinewm .= phim

        # save initial state
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
        )
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
    R, S = setup_hydromechanical_lse(coords)
    hydromech_sol = nothing
    # thermal solver
    RT, ST = setup_thermal_lse(coords)
    # gravitational solver
    RP, SP = setup_gravitational_lse(coords)
    # Pardiso MKL solver
    if use_pardiso
        pardiso_solver = Pardiso.MKLPardisoSolver()
        initialize_pardiso!(pardiso_solver, iparms_dict)
        # else
        #     F = lu(fdrand(Nx1*Ny1*6, 1, 1, matrixtype=ExtendableSparseMatrix))
    end

    # end # @timeit to "simulation_loop setup"

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
            (:to_go_Ma, s_to_Ma(endtime-timesum)),
        ]
    p = Progress(
        n_steps_val;
        showspeed=true,
        dt=0.5,
        barglyphs=BarGlyphs('|', '█', ['▁', '▂', '▃', '▄', '▅', '▆', '▇'], ' ', '|'),
        barlen=10,
    )
    for timestep in start_step_val:1:n_steps_val
        # @timeit to "set up interpolation arrays" begin
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
        # end # @timeit to "set up interpolation arrays" 

        # ---------------------------------------------------------------------
        # calculate radioactive heating
        # ---------------------------------------------------------------------
        hrsolidm, hrfluidm = calculate_radioactive_heating(
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
        )

        # ---------------------------------------------------------------------
        # compute marker properties and interpolate to staggered grid
        # ---------------------------------------------------------------------
        if use_threading
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

        # ---------------------------------------------------------------------
        # compute gravity solution
        # compute gravitational acceleration
        # ---------------------------------------------------------------------
        LP = assemble_gravitational_lse!(RHO, RP; coords=coords)
        #     @timeit to "solve gravitational LSE" begin
        SP = LP \ RP
        #     end # @timeit to "solve gravitational LSE"
        process_gravitational_solution!(SP, FI, gx, gy; coords=coords)

        # ---------------------------------------------------------------------
        # probe increasing computational timestep
        # ---------------------------------------------------------------------
        dt = min(dt*dtcoefup_val, dt_longest_val)
        dt_step_initial = dt
        maxDTcurrent = maximum(abs, DT0)
        @info "\n\n ********** begin timestep $timestep - dt = $dt s **********"

        # ---------------------------------------------------------------------
        # perform thermochemical iterations (outer iteration loop)
        # ---------------------------------------------------------------------
        for titer in 1:1:titermax_val
            #     @timeit to "thermochemical iteration (outer)" begin
            # perform thermochemical reaction
            if reaction_active
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
                )
            end

            # -----------------------------------------------------------------
            # perform hydromechanical/plastic iterations (inner iteration loop)
            # -----------------------------------------------------------------

            # save initial viscosity, yielding nodes
            #     @timeit to "save initial viscosity, yielding nodes" begin
            ETA00 .= ETA
            YNY00 .= YNY
            cur_betasolid = timestep == 1 ? 0.0 : betasolid_val
            cur_betafluid = timestep == 1 ? 0.0 : betafluid_val
            if timestep == 1
                # no elastic compaction during first timestep
                BETAPHI .= 0.0
            end
            #     end # @timeit to "save initial viscosity, yielding nodes"

            #     @timeit to "advance pressure generation" begin
            # advance pressure generation inside thermochemical iteration
            pr0 .= pr
            pf0 .= pf
            #     end # @timeit to "advance pressure generation"

            # perform plastic iterations
            for iplast in 1:1:titermax_val
                #     @timeit to "plastic iteration (inner)" begin
                @info("thermochemical iter $titer - hydromechanical iter $iplast")
                # recompute bulk viscosity at pressure nodes
                recompute_bulk_viscosity!(ETA, ETAP, ETAPHI, PHI, etaphikoef_val)
                # assemble hydromechanical system of equations
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
                )
                # solve hydromechanical system of equations
                @info "starting hydro-mechanical solver $titer-$iplast"
                #     @timeit to "solve hydromechanical system" begin
                if use_pardiso_val
                    set_phase!(pardiso_solver, Pardiso.ANALYSIS_NUM_FACT_SOLVE_REFINE)
                    pardiso(pardiso_solver, S, get_matrix(pardiso_solver, L, :N), R)
                    set_phase!(pardiso_solver, Pardiso.RELEASE_ALL)
                    pardiso(pardiso_solver, S, L, R)
                else
                    hydromech_prob = LinearProblem(L, R)
                    hydromech_sol = LinearSolve.solve(
                        hydromech_prob, UMFPACKFactorization()
                    )
                    S = hydromech_sol.u
                end
                #     end # @timeit to "solve hydromechanical system"
                @info "finished hydro-mechanical solver $titer-$iplast"
                # obtain hydromechanical observables from solution
                process_hydromechanical_solution!(
                    S, vx, vy, pr, qxD, qyD, pf; coords=coords
                )

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
                # end # @timeit to "plastic iteration"
            end # for iplast=1:1:nplast

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
            # solve temperature equation
            # ------------------------------------------------------------------
            # assemble thermal system of equations 
            LT = assemble_thermal_lse!(
                tk1, RHOCP, KX, KY, HR, HA, HS, DHP, RT, dt; coords=coords
            )
            # solve thermal system of equations
            ST = LT \ RT
            # reshape solution vector to 2D array
            tk2 .= reshape(ST, coords.Ny1, coords.Nx1)
            # compute ΔT
            @. DT = tk2 - tk1
            maxDTcurrent = maximum(abs, DT)
            @info "max DT = $maxDTcurrent K"
            # prepare next pass of thermochemical iteration
            dt = finalize_thermochemical_iteration_pass(maxDTcurrent, dt, titer)
            # evaluate iteration outcome
            if compute_thermochemical_iteration_outcome(DMP, pf, pf0, titer)
                # exit thermochemical iterations loop
                break
            end
            #     end # @timeit to "thermochemical iteration (outer)"
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
                m, xm, ym, tm, tkm, etatotalm, etavpm, YNY, YNY_inv_ETA; coords=coords
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
        )

        # ---------------------------------------------------------------------
        # interpolate DT to markers
        # ---------------------------------------------------------------------
        update_marker_temperature!(xm, ym, tkm, DT, tk2, timestep, marknum; coords=coords)

        # ---------------------------------------------------------------------
        # advance marker melt composition and porosity generation,
        # update porosity on markers for compaction,
        # update next marker porosity generation
        # ---------------------------------------------------------------------
        XWsolidm0 .= XWsolidm
        phim .= phinewm
        update_marker_porosity!(
            xm,
            ym,
            tm,
            phim,
            APHI,
            dt,
            marknum;
            phimin=phimin_val,
            phimax=phimax_val,
            coords=coords,
        )
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
            randomized=random_markers,
            coords=coords,
        )

        # ---------------------------------------------------------------------
        # update timesum
        # ---------------------------------------------------------------------
        timesum += dt
        timestep_end = now()

        # ---------------------------------------------------------------------
        #  save data for analysis and visualization
        # ---------------------------------------------------------------------
        if timestep % savematstep_val == 0
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
        if timesum > endtime
            break
        end
    end # for timestep = startstep:1:n_steps
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

"""
Parse command line arguments and feed them to the main function.

$(SIGNATURES)

# Details:
    
    - nothing

# Returns

    - parsed_args: parsed command line arguments
"""
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "config_or_output"
        help = "path to TOML configuration file (.toml) or output directory"
        default = "output"
        "--output_path", "-o"
        help = "output path for simulation data (overrides config output_dir if provided)"
        default = ""
        "--restart", "-r"
        help = "path to JLD2 checkpoint file to resume from"
        default = ""
        "--show_timer"
        help = "show timing results?"
        arg_type = Bool
        default = false
    end
    return parse_args(s)
end

"""
Runs the simulation with the given parameters.

$(SIGNATURES)

# Details

    - nothing

# Returns

    - nothing 
"""
function run_simulation(config_or_output::AbstractString="")
    if isempty(config_or_output)
        parsed_args = parse_commandline()
        target = parsed_args["config_or_output"]
        cli_output = parsed_args["output_path"]
        cli_restart = parsed_args["restart"]
        show_timer = parsed_args["show_timer"]
    else
        target = config_or_output
        cli_output = ""
        cli_restart = ""
        show_timer = false
    end

    cfg = if endswith(target, ".toml")
        load_config(target)
    elseif !isempty(target)
        def = default_config()
        SimulationConfig(;
            grid=def.grid,
            geometry=def.geometry,
            time=def.time,
            solver=def.solver,
            poroelasticity=def.poroelasticity,
            thermodynamics=def.thermodynamics,
            materials=def.materials,
            output=OutputConfig(; output_dir=target),
        )
    else
        default_config()
    end

    if !isempty(cli_restart)
        cfg = SimulationConfig(;
            grid=cfg.grid,
            geometry=cfg.geometry,
            time=cfg.time,
            solver=cfg.solver,
            poroelasticity=cfg.poroelasticity,
            thermodynamics=cfg.thermodynamics,
            materials=cfg.materials,
            output=OutputConfig(;
                output_dir=cfg.output.output_dir,
                savematstep=cfg.output.savematstep,
                visstep=cfg.output.visstep,
                restart_from=cli_restart,
            ),
        )
    end

    actual_output = isempty(cli_output) ? cfg.output.output_dir : cli_output
    actual_output = endswith(actual_output, "/") ? actual_output : actual_output * "/"
    mkpath(actual_output)

    io = open(actual_output * "Erebus_run.log", "w+")
    logger = SimpleLogger(io)
    global_logger(logger)
    if show_timer
        reset_timer!(to)
    end
    @info "=========== Erebus simulation run ==========="
    @info "system information: Apple=$(Sys.isapple()) Linux=$(Sys.islinux()) Win=$(Sys.iswindows())" Sys.cpu_info()
    @info "writing results to $actual_output"
    t1 = now()
    @info "start time = $t1"
    simulation_loop(cfg; output_path=actual_output, restart_from=cfg.output.restart_from)
    t2 = now()
    @info "end time = $t2"
    @info "total run time = $(Dates.canonicalize(
        Dates.CompoundPeriod(t2-t1)))"
    if show_timer
        show(to)
    end
    return close(io)
end

"""
    run_simulation(cfg::SimulationConfig; restart_from::AbstractString = "", output_path::AbstractString = "")

Execute a simulation using a pre-loaded `SimulationConfig` object.
"""
function run_simulation(
    cfg::SimulationConfig; restart_from::AbstractString="", output_path::AbstractString=""
)
    actual_restart = isempty(restart_from) ? cfg.output.restart_from : restart_from
    actual_output = isempty(output_path) ? cfg.output.output_dir : output_path
    return simulation_loop(cfg; output_path=actual_output, restart_from=actual_restart)
end

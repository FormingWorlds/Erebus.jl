
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
function setup_dynamic_simulation_parameters()
     # timestep counter (current), init to startstep
     timestep::Int64 = start_step
     # computational timestep (current), init to dt_longest [s]
     dt::Float64 = dt_longest
     # time sum (current), init to starttime [s]
     timesum::Float64 = start_time
     # current number of markers, init to startmarknum
     marknum::Int64 = start_marknum
     # radiogenic heat production solid phase
     hrsolidm::SVector{3, Float64} = start_hrsolidm
     # radiogenic heat production fluid phase
     hrfluidm::SVector{3, Float64} = start_hrfluidm
     # nodes yielding error vector of plastic iterations
     YERRNOD::Vector{Float64} = zeros(Float64, nplast) 
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
    alphafluidcur
    )
# @timeit to "save_state" begin
    fid = output_path * "output_" * lpad(timestep, 5, "0") * ".jld2"
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
        xsize,
        ysize,
        xcenter,
        ycenter,  
        Nx,
        Ny,
        Nx1,
        Ny1,
        Nxm,
        Nym,
        dx,
        dy,
        dxm,
        dym,
        x,
        y,
        xvx,
        yvx,
        xvy,
        yvy,
        xp,
        yp,
        xxm,
        yym,
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
        alphafluidcur
    )
# end # @timeit to "save_state"
    return nothing
end

"""
Main simulation loop: run calculations with timestepping.

$(SIGNATURES)

# Details

    - output_path: Absolute path where to save simulation output files

# Returns
    
    - nothing
"""
function simulation_loop(output_path)
    output_path = endswith(output_path, "/") ? output_path : output_path * "/"
# @timeit to "simulation_loop setup" begin
    # -------------------------------------------------------------------------
    # set up dynamic simulation parameters from given static parameters"
    # -------------------------------------------------------------------------
    timestep,
    dt,
    timesum,
    marknum,
    hrsolidm,
    hrfluidm,
    YERRNOD = setup_dynamic_simulation_parameters()

    @info "Simulation layout" Nx Ny xsize dx dy ysize rplanet rcrust marknum
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
    (
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
        DT,
        DT0,
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
        DMP,
        DHP,
        XWS
    ) = setup_staggered_grid_properties()
    (
        ETA5,
        ETA00,
        YNY5,
        YNY00,
        YNY_inv_ETA,
        DSXY,
        DSY,
        EII,
        SII,
        DSXX,
        tk0
    ) = setup_staggered_grid_properties_helpers()

    # -------------------------------------------------------------------------
    # set up markers"
    # -------------------------------------------------------------------------
    mdis, mnum = setup_marker_geometry_helpers()
    (
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
        XWsolidm0
    ) = setup_marker_properties(marknum)
    (
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
        alphafluidcur
    ) = setup_marker_properties_helpers(marknum)
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
        XWsolidm0
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
        alphafluidcur
    )

    # ---------------------------------------------------------------------
     # set up interpolation arrays"
    # ---------------------------------------------------------------------
    (
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
        DMPSUM,
        DHPSUM,
        XWSSUM,
        WTPSUM
    ) = setup_interpolated_properties()

    # -------------------------------------------------------------------------
    # set up of matrices for global grav/thermal/hydromechanical solvers"
    # -------------------------------------------------------------------------
    # hydromechanical solver
    R, S = setup_hydromechanical_lse()
    hydromech_sol = nothing
    # thermal solver
    RT, ST = setup_thermal_lse()
    # gravitational solver
    RP, SP = setup_gravitational_lse()
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
    generate_showvalues(timestep, marknum, maxT, dt, timesum) = () -> [
        (:timestep, timestep),
        (:marknum, marknum),
        (:maxT_K, maxT),
        (:dt_s, dt),
        (:timesum_Ma, s_to_Ma(timesum)),
        (:to_go_Ma, s_to_Ma(endtime-timesum))
    ]
    p = Progress(
        n_steps;
        showspeed=true,
        dt=0.5,
        barglyphs=BarGlyphs(
            '|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',), barlen=10)
    for timestep = start_step:1:n_steps
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
            WTPSUM
        )
# end # @timeit to "set up interpolation arrays" 

        # ---------------------------------------------------------------------
        # calculate radioactive heating
        # ---------------------------------------------------------------------
        hrsolidm, hrfluidm = calculate_radioactive_heating(
            hr_al, hr_fe, timesum)

        # ---------------------------------------------------------------------
        # compute marker properties and interpolate to staggered grid
        # ---------------------------------------------------------------------
        for m=1:1:marknum
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
                marker_property_mode
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
                WTSUM
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
                WTXSUM
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
                WTYSUM
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
                WTPSUM
            )
        end # for m=1:1:marknum

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
            YNY
        )

        # ---------------------------------------------------------------------
        # compute physical properties of Vx nodes
        # ---------------------------------------------------------------------
        compute_vx_node_properties!(
            RHOXSUM,
            RHOFXSUM,
            KXSUM,
            PHIXSUM,
            RXSUM,
            WTXSUM,
            RHOX,
            RHOFX,
            KX,
            PHIX,
            RX
        )

        # ---------------------------------------------------------------------
        # compute physical properties of Vy nodes
        # ---------------------------------------------------------------------
        compute_vy_node_properties!(
            RHOYSUM,
            RHOFYSUM,
            KYSUM,
            PHIYSUM,
            RYSUM,
            WTYSUM,
            RHOY,
            RHOFY,
            KY,
            PHIY,
            RY
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
            BETAPHI
        )

        # ---------------------------------------------------------------------
        # apply thermal boundary conditions for interpolated temperature
        # ---------------------------------------------------------------------
        apply_insulating_boundary_conditions!(tk1)

        # ---------------------------------------------------------------------
        # compute gravity solution
        # compute gravitational acceleration
        # ---------------------------------------------------------------------
        LP = assemble_gravitational_lse!(RHO, RP)
#     @timeit to "solve gravitational LSE" begin
        SP = LP \ RP
#     end # @timeit to "solve gravitational LSE"
        process_gravitational_solution!(SP, FI, gx, gy)

        # ---------------------------------------------------------------------
        # probe increasing computational timestep
        # ---------------------------------------------------------------------
        dt = min(dt*dtcoefup, dt_longest)
        @info "\n\n ********** begin timestep $timestep - dt = $dt s **********"

        # ---------------------------------------------------------------------
        # perform thermochemical iterations (outer iteration loop)
        # ---------------------------------------------------------------------
        for titer=1:1:titermax
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
                    titer
                )
            end

            # -----------------------------------------------------------------
            # perform hydromechanical/plastic iterations (inner iteration loop)
            # -----------------------------------------------------------------

            # save initial viscosity, yielding nodes
#     @timeit to "save initial viscosity, yielding nodes" begin
            ETA00 .= ETA
            YNY00 .= YNY
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
            for iplast=1:1:titermax
#     @timeit to "plastic iteration (inner)" begin
                @info(
                    "thermochemical iter $titer - hydromechanical iter $iplast")
                # recompute bulk viscosity at pressure nodes
                recompute_bulk_viscosity!(ETA, ETAP, ETAPHI, PHI, etaphikoef)
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
                    betasolid = betasolid,
                    betafluid = betafluid
                )
                # solve hydromechanical system of equations
                @info "starting hydro-mechanical solver $titer-$iplast"
#     @timeit to "solve hydromechanical system" begin
                if use_pardiso
                    set_phase!(
                        pardiso_solver, Pardiso.ANALYSIS_NUM_FACT_SOLVE_REFINE)
                    pardiso(
                        pardiso_solver,
                        S,
                        get_matrix(pardiso_solver, L, :N),
                        R
                    )
                    set_phase!(pardiso_solver, Pardiso.RELEASE_ALL)
                    pardiso(pardiso_solver, S, L, R)
                else
                    hydromech_prob = LinearProblem(L, R)
                    hydromech_sol = LinearSolve.solve(hydromech_prob, UMFPACKFactorization())
                    S = hydromech_sol.u
                end
#     end # @timeit to "solve hydromechanical system"
                @info "finished hydro-mechanical solver $titer-$iplast"
                # obtain hydromechanical observables from solution
                process_hydromechanical_solution!(
                    S,
                    vx,
                    vy,
                    pr,
                    qxD,
                    qyD,
                    pf
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
                    betasolid = betasolid
                )

                # compute fluid velocities
                compute_fluid_velocities!(
                    PHIX,
                    PHIY,
                    qxD,
                    qyD,
                    vx,
                    vy,
                    vxf,
                    vyf
                )

                # adapt timestep for displacement
                dt = compute_displacement_timestep(
                    vx,
                    vy,
                    vxf,
                    vyf,
                    dt,
                    aphimax
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
                    dt
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
                    betasolid = betasolid
                )
                # symmetrize P node observables
                symmetrize_p_node_observables!(
                    SXX,
                    APHI,
                    PHI,
                    pr,
                    pf,
                    ps
                )
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
                    iplast
                )
                    # exit plastic iterations loop    
                    break 
                else
                    # prepare next pass of plastic iteration 
                    dt = finalize_plastic_iteration_pass!(
                        ETA,
                        ETA5,
                        ETA00,
                        YNY,
                        YNY5,
                        YNY00,
                        YNY_inv_ETA,
                        dt,
                        iplast
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
                pf
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
                HA, tk1, ALPHA, ALPHAF, PHI, vx, vy, vxf, vyf, ps, pf)

            # ------------------------------------------------------------------
            # solve temperature equation
            # ------------------------------------------------------------------
            # assemble thermal system of equations 
            LT = assemble_thermal_lse!(
                tk1, RHOCP, KX, KY, HR, HA, HS, DHP, RT, dt)
            # solve thermal system of equations
            ST = LT \ RT
            # reshape solution vector to 2D array
            tk2 .= reshape(ST, Ny1, Nx1)
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
        @threads for m = 1:1:marknum
            update_marker_viscosity!(
                m, xm, ym, tm, tkm, etatotalm, etavpm, YNY, YNY_inv_ETA)
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
            marknum
        )

        # ---------------------------------------------------------------------
        # interpolate DSXX, DSXY to markers
        # ---------------------------------------------------------------------
        update_marker_stress!(xm, ym, sxxm, sxym, DSXX, DSXY, marknum)

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
            marker_property_mode
        )

        # ---------------------------------------------------------------------
        # interpolate DT to markers
        # ---------------------------------------------------------------------
        update_marker_temperature!(xm, ym, tkm, DT, tk2, timestep, marknum)

        # ---------------------------------------------------------------------
        # advance marker melt composition and porosity generation,
        # update porosity on markers for compaction,
        # update next marker porosity generation
        # ---------------------------------------------------------------------
        XWsolidm0 .= XWsolidm
        phim .= phinewm
        update_marker_porosity!(xm, ym, tm, phim, APHI, dt, marknum)
        phinewm .= phim

        # ---------------------------------------------------------------------
        # interpolate melt composition from markers to P nodes
        # --------------------------------------------------------------------- 
        update_p_nodes_melt_composition!(
            xm, ym, XWsolidm0, XWS, XWSSUM, WTPSUM, marknum)

        # ---------------------------------------------------------------------
        # compute velocity in P nodes,
        # compute fluid velocity in P nodes including boundary conditions
        # ---------------------------------------------------------------------
        compute_velocities!(vx, vy, vxf, vyf, vxp, vyp, vxpf, vypf)

        # ---------------------------------------------------------------------
        # compute rotation rate in basic nodes
        # ---------------------------------------------------------------------
        compute_rotation_rate!(vx, vy, wyx)
        
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
            marker_property_mode
        )

        # ---------------------------------------------------------------------
        # backtrack P nodes: Ptotal with RK4,
        # backtrack P nodes: Pfluid with RK4
        # ---------------------------------------------------------------------
        backtrace_pressures_rk4!(
            pr, pr0, ps, ps0, pf, pf0, vx, vy, vxf, vyf, dt)

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
            randomized=random_markers
        )

        # ---------------------------------------------------------------------
        # update timesum
        # ---------------------------------------------------------------------
        timesum += dt
        timestep_end = now() 

        # ---------------------------------------------------------------------
        #  save data for analysis and visualization
        # ---------------------------------------------------------------------
        if timestep % savematstep == 0
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
                alphafluidcur
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
        next!(p; showvalues = generate_showvalues(
            timestep, marknum, maxT, dt, timesum))

        # ---------------------------------------------------------------------
        # finish timestep
        # ---------------------------------------------------------------------
        if timesum > endtime
            break
        end
    end # for timestep = startstep:1:n_steps
end # function simulation loop

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
        "output_path"
            help = "output path for simulation data"
            required = true
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
function run_simulation()
    parsed_args = parse_commandline()
    output_path = parsed_args["output_path"]
    show_timer = parsed_args["show_timer"]
    mkpath(output_path)
    io = open(output_path * "Erebus_run.log", "w+")
    logger = SimpleLogger(io)
    global_logger(logger)
    if show_timer
        reset_timer!(to)
    end
    @info "=========== Erebus simulation run ==========="
    @info "system information: Apple=$(Sys.isapple()) Linux=$(Sys.islinux()) Win=$(Sys.iswindows())" Sys.cpu_info()
    @info "writing results to $output_path"
    t1 = now()
    @info "start time = $t1"
    simulation_loop(output_path)
    t2 = now()
    @info "end time = $t2"
    @info "total run time = $(Dates.canonicalize(
        Dates.CompoundPeriod(t2-t1)))"
    if show_timer
        show(to)
    end
    close(io)
end

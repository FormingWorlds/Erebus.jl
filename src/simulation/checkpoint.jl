
"""
Convert seconds to Ma (millions of years).

$(SIGNATURES)

# Details

    - s: period in seconds

# Returns

    - Ma: period in millions of years
"""
function s_to_Ma(s::Real; yearlength::Real=31557600.0)
    return s / (yearlength * 1e6)
end

"""
    init_telemetry(output_dir, filename="telemetry.csv")

Initialize telemetry CSV file with header row and return an open writable IO stream.
"""
function init_telemetry(
    output_dir::AbstractString, filename::AbstractString="telemetry.csv"; append::Bool=false
)
    mkpath(output_dir)
    filepath = joinpath(output_dir, filename)
    if append && isfile(filepath)
        return open(filepath, "a")
    end
    io = open(filepath, "w")
    header = join(
        [
            "step",
            "time_Ma",
            "dt_yr",
            "rplanet",
            "rcore",
            "T_peak",
            "T_mean",
            "phi_max",
            "phi_mean",
            "M_outgassed_total",
            "M_vent_H2O",
            "M_atm",
            "M_escaped",
            "F_melt_max",
            "F_melt_mean",
        ],
        ",",
    )
    println(io, header)
    flush(io)
    return io
end

"""
    stream_telemetry_row!(io, step, time_Ma, dt_yr, rplanet, rcore, T_peak, T_mean,
                          phi_max, phi_mean, M_outgassed, M_vent_H2O, M_atm, M_escaped,
                          F_melt_max, F_melt_mean)

Stream a single row of scalar diagnostic telemetry to the given IO stream.
"""
function stream_telemetry_row!(
    io::IO,
    step::Integer,
    time_Ma::Real,
    dt_yr::Real,
    rplanet::Real,
    rcore::Real,
    T_peak::Real,
    T_mean::Real,
    phi_max::Real,
    phi_mean::Real,
    M_outgassed::Real,
    M_vent_H2O::Real,
    M_atm::Real,
    M_escaped::Real,
    F_melt_max::Real,
    F_melt_mean::Real,
)
    println(
        io,
        string(
            step,
            ",",
            time_Ma,
            ",",
            dt_yr,
            ",",
            rplanet,
            ",",
            rcore,
            ",",
            T_peak,
            ",",
            T_mean,
            ",",
            phi_max,
            ",",
            phi_mean,
            ",",
            M_outgassed,
            ",",
            M_vent_H2O,
            ",",
            M_atm,
            ",",
            M_escaped,
            ",",
            F_melt_max,
            ",",
            F_melt_mean,
        ),
    )
    flush(io)
    return nothing
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
    # computational timestep (current), init to dt_initial [s]
    dt::Float64 = cfg.time.dt_initial * cfg.time.yearlength
    # time sum (current), init to start_time [s]
    timesum::Float64 = cfg.time.start_time * cfg.time.yearlength
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
    coords::Union{Nothing,GridCoordinates}=nothing,
    phim0_val=phim0,
    M_vent_total::Real=0.0,
    M_vent_H2O_total::Real=0.0,
    M_vent_C_total::Real=0.0,
    M_vent_N_total::Real=0.0,
    M_vent_S_total::Real=0.0,
    M_atm_total::Real=0.0,
    M_escaped_total::Real=0.0,
    P_amb::Real=10.0,
    S_vent::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    Xfem=nothing,
    Xfem0=nothing,
    Xfe_bulk=nothing,
    M_atm_species::Union{Nothing,Dict{Symbol,Float64}}=nothing,
    M_escaped_species::Union{Nothing,Dict{Symbol,Float64}}=nothing,
    XH2Om=nothing,
    XCm=nothing,
    XNm=nothing,
    XSm=nothing,
    Xfe_H_m=nothing,
    Xfe_C_m=nothing,
    Xfe_N_m=nothing,
    Xfe_S_m=nothing,
    core_budgets=nothing,
    Xmin_troilite_m=nothing,
    Xmin_schreibersite_m=nothing,
    Xmin_cohenite_m=nothing,
    Xmin_graphite_m=nothing,
    Xmin_nitride_m=nothing,
    Xmin_metal_matrix_m=nothing,
    regional_mineral_modes=nothing,
    DT0::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    t_accreted=nothing,
    M_accreted_total=nothing,
    M_planet_val=nothing,
    rplanet::Union{Nothing,Real}=nothing,
    rcore::Union{Nothing,Real}=nothing,
    telescope_level::Union{Nothing,Integer}=nothing,
    hcnspo_props=nothing,
    redox_props=nothing,
    atm_state::Union{Nothing,AtmosphereState}=nothing,
    F_extract_m=nothing,
    cfg::Union{Nothing,SimulationConfig}=nothing,
)
    fid = output_path * "output_" * lpad(timestep, 5, "0") * ".jld2"
    @unpack_coords coords Nx Ny Nx1 Ny1 Nxm Nym dx dy dxm dym
    @unpack_coords coords x y xvx yvx xvy yvy xp yp xxm yym xsize ysize xcenter ycenter
    jldsave(
        fid;
        timestep,
        dt,
        Δtreaction,
        reaction_rate_coeff_mode,
        marker_property_mode,
        timesum,
        marknum,
        phim0=phim0_val,
        M_vent_total,
        M_vent_H2O_total,
        M_vent_C_total,
        M_vent_N_total,
        M_vent_S_total,
        M_atm_total,
        M_escaped_total,
        P_amb,
        S_vent=S_vent === nothing ? zeros(Float64, Ny1_val, Nx1_val) : S_vent,
        ratio_al,
        t_half_al,
        dsubgrids=cfg === nothing ? 0.0 : cfg.solver.dsubgrids,
        dsubgridt=cfg === nothing ? 0.0 : cfg.solver.dsubgridt,
        hr_al=cfg === nothing ? true : cfg.thermodynamics.hr_al,
        hr_fe=cfg === nothing ? false : cfg.thermodynamics.hr_fe,
        rplanet=rplanet !== nothing ? Float64(rplanet) : 50000.0,
        rcore=rcore !== nothing ? Float64(rcore) : 0.0,
        telescope_level=telescope_level !== nothing ? Int(telescope_level) : 0,
        rcrust=cfg === nothing ? 50000.0 : cfg.geometry.rcrust,
        psurface=cfg === nothing ? 1000.0 : cfg.geometry.psurface,
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
        DT0=DT0 === nothing ? zeros(Float64, Ny1_val, Nx1_val) : DT0,
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
        alphafluidcur,
        (F_extract_m !== nothing ? (; F_extract_m) : (;))...,
        (Xfem !== nothing ? (; Xfem, Xfem0, Xfe_bulk) : (;))...,
        (XH2Om !== nothing ? (; XH2Om, XCm, XNm, XSm) : (;))...,
        (Xfe_H_m !== nothing ? (; Xfe_H_m, Xfe_C_m, Xfe_N_m, Xfe_S_m) : (;))...,
        (core_budgets !== nothing ? (; core_budgets) : (;))...,
        (M_atm_species !== nothing ? (; M_atm_species, M_escaped_species) : (;))...,
        (
            if Xmin_troilite_m !== nothing
                (;
                    Xmin_troilite_m,
                    Xmin_schreibersite_m,
                    Xmin_cohenite_m,
                    Xmin_graphite_m,
                    Xmin_nitride_m,
                    Xmin_metal_matrix_m,
                )
            else
                (;)
            end
        )...,
        (regional_mineral_modes !== nothing ? (; regional_mineral_modes) : (;))...,
        (t_accreted !== nothing ? (; t_accreted) : (;))...,
        (M_accreted_total !== nothing ? (; M_accreted_total) : (;))...,
        (M_planet_val !== nothing ? (; M_planet_val) : (;))...,
        (hcnspo_props !== nothing ? (; hcnspo_props...) : (;))...,
        (redox_props !== nothing ? (; redox_props...) : (;))...,
        (
            if atm_state !== nothing
                (;
                    atm_M_atm=atm_state.M_atm,
                    atm_M_escaped=atm_state.M_escaped,
                    atm_P_surf=atm_state.P_surf,
                    atm_T_surf_eq=atm_state.T_surf_eq,
                    atm_tau_LW=atm_state.tau_LW,
                    atm_M_env_bound=atm_state.M_env_bound,
                    atm_F_net_rad=atm_state.F_net_rad,
                    atm_h_rad_eff=atm_state.h_rad_eff,
                )
            else
                (;)
            end
        )...,
    )
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

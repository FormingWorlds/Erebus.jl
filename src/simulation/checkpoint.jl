
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
            "dt_aphimax_max",
            "n_flips_last",
            "n_flips_total",
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
                          F_melt_max, F_melt_mean, dt_aphimax_max)

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
    dt_aphimax_max::Real=0.0,
    n_flips_last::Integer=0,
    n_flips_total::Integer=0,
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
            ",",
            dt_aphimax_max,
            ",",
            n_flips_last,
            ",",
            n_flips_total,
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
    cfg::SimulationConfig=default_config();
    coords::GridCoordinates=default_grid_coordinates(),
)
    # timestep counter (current), init to startstep
    timestep::Int64 = cfg.time.start_step
    # computational timestep (current), init to dt_initial [s]
    dt::Float64 = cfg.time.dt_initial * cfg.time.yearlength
    # time sum (current), init to start_time [s]
    timesum::Float64 = cfg.time.start_time * cfg.time.yearlength
    # current number of markers, init to startmarknum
    marknum::Int64 = coords.start_marknum
    # radiogenic heat production solid phase
    hrsolidm::SVector{3,Float64} = start_hrsolidm
    # radiogenic heat production fluid phase
    hrfluidm::SVector{3,Float64} = start_hrfluidm
    # nodes yielding error vector of plastic iterations
    YERRNOD::Vector{Float64} = zeros(Float64, cfg.solver.max_plastic_iterations)
    return timestep, dt, timesum, marknum, hrsolidm, hrfluidm, YERRNOD
end # function setup_dynamic_simulation_parameters()

"""
    _collect_checkpoint_marker_dict(markers; kwargs...)

Collect marker arrays and optional simulation properties into a dictionary for checkpointing.
"""
function _collect_checkpoint_marker_dict(
    markers::Union{Nothing,MarkerArrays};
    xm=nothing,
    ym=nothing,
    tm=nothing,
    tkm=nothing,
    sxxm=nothing,
    sxym=nothing,
    etavpm=nothing,
    phim=nothing,
    rhototalm=nothing,
    rhocptotalm=nothing,
    etatotalm=nothing,
    hrtotalm=nothing,
    ktotalm=nothing,
    tkm_rhocptotalm=nothing,
    etafluidcur_inv_kphim=nothing,
    inv_gggtotalm=nothing,
    fricttotalm=nothing,
    cohestotalm=nothing,
    tenstotalm=nothing,
    rhofluidcur=nothing,
    alphasolidcur=nothing,
    alphafluidcur=nothing,
    XWsolidm0=nothing,
    F_extract_m=nothing,
    Xfem=nothing,
    Xfem0=nothing,
    Xfe_bulk=nothing,
    XH2Om=nothing,
    XCm=nothing,
    XNm=nothing,
    XSm=nothing,
    Xfe_H_m=nothing,
    Xfe_C_m=nothing,
    Xfe_N_m=nothing,
    Xfe_S_m=nothing,
    core_budgets=nothing,
    M_atm_species=nothing,
    M_escaped_species=nothing,
    Xmin_troilite_m=nothing,
    Xmin_schreibersite_m=nothing,
    Xmin_cohenite_m=nothing,
    Xmin_nitride_m=nothing,
    Xmin_metal_matrix_m=nothing,
    Xmin_graphite_m=nothing,
    X_graphite_m=nothing,
    regional_mineral_modes=nothing,
    t_accreted=nothing,
    M_accreted_total=nothing,
    M_planet_val=nothing,
    hcnspo_props=nothing,
    redox_props=nothing,
)
    dict = Dict{Symbol,Any}()
    if markers !== nothing
        for fn in fieldnames(CoreGroup)
            dict[fn] = getfield(markers.core, fn)
        end
        for grp in values(markers.groups)
            for fn in fieldnames(typeof(grp))
                dict[fn] = getfield(grp, fn)
            end
        end
    else
        dict[:xm] = xm
        dict[:ym] = ym
        dict[:tm] = tm
        dict[:tkm] = tkm
        dict[:sxxm] = sxxm
        dict[:sxym] = sxym
        dict[:etavpm] = etavpm
        dict[:phim] = phim
        dict[:rhototalm] = rhototalm
        dict[:rhocptotalm] = rhocptotalm
        dict[:etatotalm] = etatotalm
        dict[:hrtotalm] = hrtotalm
        dict[:ktotalm] = ktotalm
        dict[:tkm_rhocptotalm] = tkm_rhocptotalm
        dict[:etafluidcur_inv_kphim] = etafluidcur_inv_kphim
        dict[:inv_gggtotalm] = inv_gggtotalm
        dict[:fricttotalm] = fricttotalm
        dict[:cohestotalm] = cohestotalm
        dict[:tenstotalm] = tenstotalm
        dict[:rhofluidcur] = rhofluidcur
        dict[:alphasolidcur] = alphasolidcur
        dict[:alphafluidcur] = alphafluidcur
        dict[:XWsolidm0] = XWsolidm0
        F_extract_m !== nothing && (dict[:F_extract_m] = F_extract_m)
        Xfem !== nothing && (dict[:Xfem]=Xfem; dict[:Xfem0]=Xfem0; dict[:Xfe_bulk]=Xfe_bulk)
        XH2Om !== nothing &&
            (dict[:XH2Om]=XH2Om; dict[:XCm]=XCm; dict[:XNm]=XNm; dict[:XSm]=XSm)
        Xfe_H_m !== nothing && (
            dict[:Xfe_H_m]=Xfe_H_m;
            dict[:Xfe_C_m]=Xfe_C_m;
            dict[:Xfe_N_m]=Xfe_N_m;
            dict[:Xfe_S_m]=Xfe_S_m
        )
        Xmin_troilite_m !== nothing && (
            dict[:Xmin_troilite_m]=Xmin_troilite_m;
            dict[:Xmin_schreibersite_m]=Xmin_schreibersite_m;
            dict[:Xmin_cohenite_m]=Xmin_cohenite_m;
            dict[:Xmin_nitride_m]=Xmin_nitride_m;
            dict[:Xmin_metal_matrix_m]=Xmin_metal_matrix_m
        )
        Xmin_graphite_m !== nothing && (dict[:Xmin_graphite_m] = Xmin_graphite_m)
        X_graphite_m !== nothing && (dict[:X_graphite_m] = X_graphite_m)
        t_accreted !== nothing && (dict[:t_accreted] = t_accreted)
        if hcnspo_props !== nothing
            for fn in fieldnames(typeof(hcnspo_props))
                dict[fn] = getfield(hcnspo_props, fn)
            end
        end
        if redox_props !== nothing
            for fn in fieldnames(typeof(redox_props))
                dict[fn] = getfield(redox_props, fn)
            end
        end
    end
    regional_mineral_modes !== nothing &&
        (dict[:regional_mineral_modes] = regional_mineral_modes)
    core_budgets !== nothing && (dict[:core_budgets] = core_budgets)
    M_atm_species !== nothing &&
        (dict[:M_atm_species]=M_atm_species; dict[:M_escaped_species]=M_escaped_species)
    M_accreted_total !== nothing && (dict[:M_accreted_total] = M_accreted_total)
    M_planet_val !== nothing && (dict[:M_planet_val] = M_planet_val)
    return dict
end

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
    coords::GridCoordinates=default_grid_coordinates(),
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
    X_graphite_m=nothing,
    markers::Union{Nothing,MarkerArrays}=nothing,
    cfg::Union{Nothing,SimulationConfig}=nothing,
    transfer_log=nothing,
)
    marker_props = _collect_checkpoint_marker_dict(
        markers;
        xm=xm,
        ym=ym,
        tm=tm,
        tkm=tkm,
        sxxm=sxxm,
        sxym=sxym,
        etavpm=etavpm,
        phim=phim,
        rhototalm=rhototalm,
        rhocptotalm=rhocptotalm,
        etatotalm=etatotalm,
        hrtotalm=hrtotalm,
        ktotalm=ktotalm,
        tkm_rhocptotalm=tkm_rhocptotalm,
        etafluidcur_inv_kphim=etafluidcur_inv_kphim,
        inv_gggtotalm=inv_gggtotalm,
        fricttotalm=fricttotalm,
        cohestotalm=cohestotalm,
        tenstotalm=tenstotalm,
        rhofluidcur=rhofluidcur,
        alphasolidcur=alphasolidcur,
        alphafluidcur=alphafluidcur,
        XWsolidm0=XWsolidm0,
        F_extract_m=F_extract_m,
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
        core_budgets=core_budgets,
        M_atm_species=M_atm_species,
        M_escaped_species=M_escaped_species,
        Xmin_troilite_m=Xmin_troilite_m,
        Xmin_schreibersite_m=Xmin_schreibersite_m,
        Xmin_cohenite_m=Xmin_cohenite_m,
        Xmin_nitride_m=Xmin_nitride_m,
        Xmin_metal_matrix_m=Xmin_metal_matrix_m,
        Xmin_graphite_m=Xmin_graphite_m,
        X_graphite_m=X_graphite_m,
        regional_mineral_modes=regional_mineral_modes,
        t_accreted=t_accreted,
        M_accreted_total=M_accreted_total,
        M_planet_val=M_planet_val,
        hcnspo_props=hcnspo_props,
        redox_props=redox_props,
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
        marker_props...,
        (transfer_log !== nothing ? (; transfer_log) : (;))...,
        (
            if atm_state !== nothing
                (;
                    atm_elem=atm_state.elem,
                    atm_species=atm_state.species,
                    atm_escaped=atm_state.escaped,
                    atm_dO_buffer=atm_state.dO_buffer,
                    atm_log10_fO2=atm_state.log10_fO2,
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

"""
    compare_restart_configs(cfg_saved::SimulationConfig, cfg_current::SimulationConfig)::Vector{String}

Compare configuration of saved checkpoint against restart configuration.
Differences in `[output]`, `time.{n_steps, endtime, start_step}`, and `solver.seed` are allowed.
All other differences return as `\"section.field\"` strings.
"""
function compare_restart_configs(
    cfg_saved::SimulationConfig, cfg_current::SimulationConfig
)::Vector{String}
    d_saved = config_to_dict(cfg_saved)
    d_curr = config_to_dict(cfg_current)
    diffs = String[]

    all_secs = sort(collect(union(keys(d_saved), keys(d_curr))))
    for sec in all_secs
        sec == "output" && continue
        s_saved = get(d_saved, sec, Dict{String,Any}())
        s_curr = get(d_curr, sec, Dict{String,Any}())
        all_flds = sort(collect(union(keys(s_saved), keys(s_curr))))
        for fld in all_flds
            if (sec == "time" && fld in ("n_steps", "endtime", "start_step")) ||
                (sec == "solver" && fld == "seed")
                continue
            end
            val_saved = get(s_saved, fld, nothing)
            val_curr = get(s_curr, fld, nothing)
            if !isequal(val_saved, val_curr)
                push!(diffs, "$sec.$fld")
            end
        end
    end
    return diffs
end

"""
    save_state(output_dir::AbstractString, state::SimulationState, coords::GridCoordinates, cfg::SimulationConfig)

Serialize simulation state, grid coordinates, and configuration to a JLD2 checkpoint archive.
Emits checkpoint schema version 2.
"""
function save_state(
    output_dir::AbstractString,
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig,
)
    filepath = if endswith(output_dir, ".jld2")
        output_dir
    else
        joinpath(output_dir, "output_" * lpad(state.timestep, 5, "0") * ".jld2")
    end
    mkpath(dirname(filepath))

    grid_dict = Dict{String,Any}(
        string(fn) => getfield(state.grids, fn) for fn in fieldnames(GridArrays)
    )

    marker_dict = Dict{String,Any}()
    for fn in fieldnames(CoreGroup)
        marker_dict[string(fn)] = getfield(state.markers.core, fn)
    end
    for grp in values(state.markers.groups)
        for fn in fieldnames(typeof(grp))
            marker_dict[string(fn)] = getfield(grp, fn)
        end
    end

    coord_dict = Dict{String,Any}(
        string(fn) => getfield(coords, fn) for fn in fieldnames(GridCoordinates)
    )

    acc = state.accumulators

    jldopen(filepath, "w") do f
        f["schema_version"] = 2
        f["cfg"] = cfg
        f["coords"] = coords
        f["timestep"] = state.timestep
        f["dt"] = state.dt
        f["timesum"] = state.timesum
        f["marknum"] = length(state.markers)
        f["rng"] = copy(state.rng)
        f["transfers"] = deepcopy(state.transfers)
        f["transfer_log"] = deepcopy(state.transfers)
        f["S_vent"] = state.grids.S_vent_grid
        f["accumulators"] = copy(state.accumulators)
        clean_timer = copy(state.timer)
        empty!(clean_timer.timer_stack)
        f["timer"] = clean_timer

        for fn in (
            :M_vent_total,
            :M_vent_H2O_total,
            :M_vent_C_total,
            :M_vent_N_total,
            :M_vent_S_total,
            :M_atm_total,
            :M_escaped_total,
            :P_amb,
            :rplanet,
            :rcore,
            :telescope_level,
            :M_accreted_total,
            :M_planet_val,
            :max_v_seg_prev,
        )
            f[string(fn)] = getfield(acc, fn)
        end
        f["planet_xcenter"] = acc.xcenter
        f["planet_ycenter"] = acc.ycenter
        acc.M_atm_species !== nothing && (f["M_atm_species"] = acc.M_atm_species)
        acc.M_escaped_species !== nothing &&
            (f["M_escaped_species"] = acc.M_escaped_species)
        acc.core_budgets !== nothing && (f["core_budgets"] = acc.core_budgets)
        acc.regional_mineral_modes !== nothing &&
            (f["regional_mineral_modes"] = acc.regional_mineral_modes)

        for (k, v) in grid_dict
            f[k] = v
        end
        for (k, v) in marker_dict
            f[k] = v
        end
        for (k, v) in coord_dict
            f[k] = v
        end

        if state.atm !== nothing
            f["atm_state"] = state.atm
            for fn in (
                :elem,
                :species,
                :escaped,
                :dO_buffer,
                :log10_fO2,
                :M_atm,
                :M_escaped,
                :P_surf,
                :T_surf_eq,
                :tau_LW,
                :M_env_bound,
                :F_net_rad,
                :h_rad_eff,
            )
                f["atm_" * string(fn)] = getfield(state.atm, fn)
            end
        end
    end

    return filepath
end

"""
    load_simulation_state(path::AbstractString; force_restart_config::Bool=false, current_cfg::Union{Nothing,SimulationConfig}=nothing)

Load simulation state, coordinates, and saved configuration from a schema version 2 checkpoint archive.
Returns `(state::SimulationState, coords::GridCoordinates, cfg_saved::SimulationConfig)`.

# Throws
- `CheckpointError` on missing file, `schema_version < 2`, missing keys, or configuration mismatch.
"""
function load_simulation_state(
    path::AbstractString;
    force_restart_config::Bool=false,
    current_cfg::Union{Nothing,SimulationConfig}=nothing,
)
    isfile(path) || throw(CheckpointError("Checkpoint file does not exist: $path"))

    data = try
        JLD2.load(path)
    catch err
        throw(CheckpointError("Failed to open checkpoint archive '$path': $err"))
    end

    schema_ver = get(data, "schema_version", 0)
    if schema_ver < 2
        throw(
            CheckpointError(
                "Unsupported checkpoint schema version: got $schema_ver, expected >= 2"
            ),
        )
    end

    for req in ("cfg", "timestep", "dt", "timesum", "rng", "transfers")
        haskey(data, req) ||
            throw(CheckpointError("Missing required checkpoint key: '$req'"))
    end

    cfg_saved = if data["cfg"] isa SimulationConfig
        data["cfg"]
    else
        throw(CheckpointError("Invalid 'cfg' payload in checkpoint archive"))
    end

    if current_cfg !== nothing
        diffs = compare_restart_configs(cfg_saved, current_cfg)
        if !isempty(diffs)
            # Grid geometry and dimension mismatches cannot be overridden (unsupported regridding)
            grid_diffs = filter(
                d -> startswith(d, "grid.") || startswith(d, "geometry."), diffs
            )
            if !isempty(grid_diffs)
                throw(
                    CheckpointError(
                        "Unsupported grid geometry override: " *
                        join(grid_diffs, ", ") *
                        ". Grid dimensions and geometry cannot be changed on restart.",
                    ),
                )
            end
            if !force_restart_config
                throw(
                    CheckpointError(
                        "Configuration mismatch between checkpoint and restart config: " *
                        join(diffs, ", ") *
                        ". Pass --force-restart-config to override.",
                    ),
                )
            else
                for d in diffs
                    @warn "Restart configuration override: $d"
                end
            end
        end
    end

    effective_cfg = current_cfg !== nothing ? current_cfg : cfg_saved

    coords = if haskey(data, "coords") && data["coords"] isa GridCoordinates
        data["coords"]
    else
        GridCoordinates(effective_cfg.grid)
    end

    grid_vals = Any[]
    for fn in fieldnames(GridArrays)
        sfn = string(fn)
        if !haskey(data, sfn)
            if fn === :Q_metric
                push!(grid_vals, nothing)
            else
                throw(CheckpointError("Missing required grid array in checkpoint: '$sfn'"))
            end
        else
            push!(grid_vals, data[sfn])
        end
    end
    grids = GridArrays(grid_vals...)

    core_vals = Any[]
    for fn in fieldnames(CoreGroup)
        sfn = string(fn)
        haskey(data, sfn) || throw(
            CheckpointError("Missing required core marker array in checkpoint: '$sfn'")
        )
        push!(core_vals, data[sfn])
    end
    core = CoreGroup(core_vals...)

    group_pairs = Pair{Symbol,Any}[]
    marknum = length(core.xm)

    # Optional marker groups
    if effective_cfg.metal_partition.active || haskey(data, "Xfem")
        metal_vals = Any[]
        for k in fieldnames(MetalGroup)
            sk = string(k)
            v = get(data, sk, nothing)
            push!(metal_vals, v !== nothing ? v : zeros(Float64, marknum))
        end
        push!(group_pairs, :metal => MetalGroup(metal_vals...))
    end

    if effective_cfg.volatiles.active || haskey(data, "XH2Om")
        vol_vals = Any[]
        for k in fieldnames(VolatilesGroup)
            sk = string(k)
            v = get(data, sk, nothing)
            push!(vol_vals, v !== nothing ? v : zeros(Float64, marknum))
        end
        push!(group_pairs, :volatiles => VolatilesGroup(vol_vals...))
    end

    if effective_cfg.redox.active || haskey(data, "nFe0_m")
        redox_vals = Any[]
        for k in fieldnames(RedoxGroup)
            sk = string(k)
            v = get(data, sk, nothing)
            push!(redox_vals, v !== nothing ? v : zeros(Float64, marknum))
        end
        push!(group_pairs, :redox => RedoxGroup(redox_vals...))
    end

    if haskey(data, "X_ice_H2O_m") ||
        effective_cfg.volatile_mixture.active ||
        effective_cfg.refractory.active
        hcnspo_vals = Any[]
        for k in fieldnames(HcnspoGroup)
            sk = string(k)
            v = get(data, sk, nothing)
            push!(hcnspo_vals, v !== nothing ? v : zeros(Float64, marknum))
        end
        push!(group_pairs, :hcnspo => HcnspoGroup(hcnspo_vals...))
    end

    if effective_cfg.phase_tracking.active || haskey(data, "Xmin_troilite_m")
        phase_vals = Any[]
        for k in fieldnames(PhaseGroup)
            sk = string(k)
            v = get(data, sk, nothing)
            push!(phase_vals, v !== nothing ? v : zeros(Float64, marknum))
        end
        push!(group_pairs, :phase => PhaseGroup(phase_vals...))
    end

    if effective_cfg.accretion.active || haskey(data, "t_accreted")
        v = get(data, "t_accreted", zeros(Float64, marknum))
        push!(group_pairs, :accretion => AccretionGroup(v))
    end

    markers = MarkerArrays(core, NamedTuple(group_pairs))

    accumulators =
        if haskey(data, "accumulators") && data["accumulators"] isa SimulationAccumulators
            copy(data["accumulators"])
        else
            SimulationAccumulators(
                Float64(get(data, "M_vent_total", 0.0)),
                Float64(get(data, "M_vent_H2O_total", 0.0)),
                Float64(get(data, "M_vent_C_total", 0.0)),
                Float64(get(data, "M_vent_N_total", 0.0)),
                Float64(get(data, "M_vent_S_total", 0.0)),
                Float64(get(data, "M_atm_total", 0.0)),
                Float64(get(data, "M_escaped_total", 0.0)),
                Float64(get(data, "P_amb", 10.0)),
                Float64(get(data, "rplanet", 50000.0)),
                Float64(get(data, "rcore", 0.0)),
                Int(get(data, "telescope_level", 0)),
                Float64(get(data, "M_accreted_total", 0.0)),
                Float64(get(data, "M_planet_val", 0.0)),
                Float64(get(data, "planet_xcenter", get(data, "xcenter", coords.xcenter))),
                Float64(get(data, "planet_ycenter", get(data, "ycenter", coords.ycenter))),
                Float64(get(data, "max_v_seg_prev", 0.0)),
                get(data, "M_atm_species", nothing),
                get(data, "M_escaped_species", nothing),
                get(data, "core_budgets", nothing),
                get(data, "regional_mineral_modes", nothing),
            )
        end

    atm = if haskey(data, "atm_state") && data["atm_state"] isa AtmosphereState
        copy(data["atm_state"])
    else
        nothing
    end

    transfers = haskey(data, "transfers") ? deepcopy(data["transfers"]) : TransferRecord[]
    rng = data["rng"]
    timer = if haskey(data, "timer") && data["timer"] isa TimerOutput
        t = copy(data["timer"])
        empty!(t.timer_stack)
        t
    else
        TimerOutput()
    end
    timestep = Int(data["timestep"])
    dt = Float64(data["dt"])
    timesum = Float64(data["timesum"])

    state = SimulationState(
        grids, markers, accumulators, transfers, atm, rng, timer, timestep, dt, timesum
    )

    return (state, coords, cfg_saved)
end

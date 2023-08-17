module Erebus

using ArgParse
using Base.Threads
using Dates
using DocStringExtensions
using ExtendableSparse
using JLD2
using LinearAlgebra
using LinearSolve
using Logging
using ProgressMeter
using Random
using SparseArrays
using StaticArrays
using TimerOutputs

export run_simulation

# include("constants.jl")
include("test_constants.jl")

if use_pardiso
    using Pardiso
# else
    # using MKL
    # BLAS.set_num_threads(4)
end

const to = TimerOutput()
const rgen = MersenneTwister(seed)

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
Set up staggered grid properties for basic, Vx, Vy, and P nodes.

$(SIGNATURES)

# Details

    - randomized: fill in random values for grid properties instead of zeros

# Returns
    - ETA : viscoplastic viscosity at basic nodes [Pa*s]
    - ETA0 : viscous viscosity at basic nodes [Pa*s]
    - GGG : shear modulus at basic nodes [Pa]
    - EXY : ϵxy at basic nodes [1/s]
    - SXY : σxy at basic nodes [1/s]
    - SXY0 : σ₀xy at basic nodes [1/s]
    - wyx : rotation rate at basic nodes [1/s]
    - COH : compressive strength at basic nodes [Pa]
    - TEN : tensile strength at basic nodes [Pa]
    - FRI : friction at basic nodes
    - YNY : plastic yielding node property at basic nodes
    - RHOX : density at Vx nodes [kg/m^3]
    - RHOFX : fluid density at Vx nodes [kg/m^3]
    - KX : thermal conductivity at Vx nodes [W/m/K]
    - PHIX : porosity at Vx nodes
    - vx : solid vx-velocity at Vx nodes [m/s]
    - vxf : fluid vx-velocity at Vx nodes [m/s]
    - RX : etafluid/kphi ratio at Vx nodes [kg m⁻³s⁻¹]
    - qxD : qx-darcy flux at Vx nodes [m/s]
    - gx : gx-gravity at Vx nodes [m/s^2]
    - RHOY : density at Vy nodes [kg/m^3]
    - RHOFY : fluid density at Vy nodes [kg/m^3]
    - KY : thermal conductivity at Vy nodes [W/m/K]
    - PHIY : porosity at Vy nodes
    - vy : solid vy-velocity at Vy nodes [m/s]
    - vyf : fluid vy-velocity at Vy nodes [m/s]
    - RY : etafluid/kphi ratio at Vy nodes [kg m⁻³s⁻¹]
    - qyD : qy-Darcy flux at Vy nodes [m/s]
    - gy : gy-gravity at Vy nodes [m/s^2]
    - RHO : density at P nodes [kg/m^3]
    - RHOCP : volumetric heat capacity at P nodes [J/m^3/K]
    - ALPHA : thermal expansion at P nodes [K⁻¹]
    - ALPHAF : fluid thermal expansion at P nodes [K⁻¹]
    - HR : radioactive heating at P nodes [W/m^3]
    - HA : adiabatic heating at P nodes [W/m^3]
    - HS : shear heating at P nodes [W/m^3]
    - ETAP : viscosity at P nodes [Pa*s]
    - GGGP : shear modulus at P nodes [Pa]
    - EXX : ϵxx at P nodes [1/s]
    - SXX : σ′xx at P nodes [1/s]
    - SXX0 : σ₀′xx at P nodes [1/s]
    - tk1 : current temperature at P nodes [K]
    - tk2 : next temperature at P nodes [K]
    - DT: temperature difference at P nodes [K]
    - DT0: previous temperature difference at P nodes [K]
    - vxp : solid vx in pressure nodes at P nodes [m/s]
    - vyp : solid vy in pressure nodes at P nodes [m/s]
    - vxpf : fluid vx in pressure nodes at P nodes [m/s]
    - vypf : fluid vy in pressure nodes at P nodes [m/s]
    - pr : total pressure at P nodes [Pa]
    - pf : fluid pressure at P nodes [Pa]
    - ps : solid pressure at P nodes [Pa]
    - pr0 : previous total pressure at P nodes [Pa]
    - pf0 : previous fluid pressure at P nodes [Pa]
    - ps0 : previous solid pressure at P nodes [Pa]
    - ETAPHI : bulk viscosity at P nodes [Pa*s]
    - BETAPHI : bulk compresibility at P nodes [Pa*s]
    - PHI : porosity at P nodes
    - APHI : Dln at P nodes [(1-ϕ)/ϕ]/Dt
    - FI : gravity potential at P nodes [J/kg]
    - DMP: mass transfer term at P nodes
    - DHP: enthalpy transfer/latent heating term at P nodes
    - XWS: wet solid fraction at P nodes
"""
function setup_staggered_grid_properties(; randomized=false)
    # basic nodes
    # viscoplastic viscosity [Pa*s]
    ETA = randomized ? rand(rgen, Ny, Nx)*1e16 : zeros(Ny, Nx)
    # viscous viscosity [Pa*s]
    ETA0 = randomized ? rand(rgen, Ny, Nx)*1e16 : zeros(Ny, Nx)
    # shear modulus [Pa]
    GGG = randomized ? rand(rgen, Ny, Nx)*1e10 : zeros(Ny, Nx)
    # epsilonxy [1/s]
    EXY = randomized ? rand(rgen, Ny, Nx)*2e-13.-1e-13 : zeros(Ny, Nx)
    # σxy [1/s]
    SXY = randomized ? rand(rgen, Ny, Nx)*1e4 : zeros(Ny, Nx)
    # σ₀xy [1/s]
    SXY0 = randomized ? rand(rgen, Ny, Nx)*1e4 : zeros(Ny, Nx)
    # rotation rate [1/s]
    wyx = randomized ? rand(rgen, Ny, Nx)*2e-14.-1e-14 : zeros(Ny, Nx)
    # compressive strength [Pa]
    COH = randomized ? rand(rgen, Ny, Nx)*1e8 : zeros(Ny, Nx)
    # tensile strength [Pa]
    TEN = randomized ? rand(rgen, Ny, Nx)*1e8 : zeros(Ny, Nx)
    # friction
    FRI = randomized ? rand(rgen, Ny, Nx) : zeros(Ny, Nx)
    # plastic yielding node property
    YNY = randomized ? rand(rgen, Bool, Ny, Nx) : zeros(Bool, Ny, Nx)
    # Vx nodes
    # density [kg/m^3]
    RHOX = randomized ? rand(rgen, Ny1, Nx1)*1e4 : zeros(Ny1, Nx1)
    # fluid density [kg/m^3]
    RHOFX = randomized ? rand(rgen, Ny1, Nx1)*1e4 : zeros(Ny1, Nx1)
    # thermal conductivity [W/m/K]
    KX = randomized ? rand(rgen, Ny1, Nx1)*1e4 : zeros(Ny1, Nx1)
    # porosity
    PHIX = randomized ? rand(rgen, Ny1, Nx1) : zeros(Ny1, Nx1)
    # solid vx-velocity [m/s]
    vx = randomized ? rand(rgen, Ny1, Nx1)*2e-9.-1e-9 : zeros(Ny1, Nx1)
    # fluid vx-velocity [m/s]
    vxf = randomized ? rand(rgen, Ny1, Nx1)*2e-9.-1e-9 : zeros(Ny1, Nx1)
    # etafluid/kphi ratio [kg m⁻³s⁻¹]
    RX = randomized ? rand(rgen, Ny1, Nx1)*1e39 : zeros(Ny1, Nx1)
    # qx-darcy flux [m/s]
    qxD = randomized ? rand(rgen, Ny1, Nx1)*2e-10.-1e-10 : zeros(Ny1, Nx1)
    # gx-gravity [m/s^2]
    gx = randomized ? rand(rgen, Ny1, Nx1)*2e-1.-1e-1 : zeros(Ny1, Nx1)
    # Vy nodes
    # density [kg/m^3]
    RHOY = randomized ? rand(rgen, Ny1, Nx1)*1e4 : zeros(Ny1, Nx1)
    # fluid density [kg/m^3]
    RHOFY = randomized ? rand(rgen, Ny1, Nx1)*1e4 : zeros(Ny1, Nx1)
    # thermal conductivity [W/m/K]
    KY = randomized ? rand(rgen, Ny1, Nx1)*1e4 : zeros(Ny1, Nx1)
    # porosity
    PHIY = randomized ? rand(rgen, Ny1, Nx1) : zeros(Ny1, Nx1)
    # solid vy-velocity [m/s]
    vy = randomized ? rand(rgen, Ny1, Nx1)*2e-9.-1e-9 : zeros(Ny1, Nx1)
    # fluid vy-velocity [m/s]
    vyf = randomized ? rand(rgen, Ny1, Nx1)*2e-9.-1e-9 : zeros(Ny1, Nx1)
    # etafluid/kphi ratio [kg m⁻³s⁻¹]
    RY = randomized ? rand(rgen, Ny1, Nx1)*1e39 : zeros(Ny1, Nx1)
    # qy-Darcy flux [m/s]
    qyD = randomized ? rand(rgen, Ny1, Nx1)*2e-10.-1e-10 : zeros(Ny1, Nx1)
    # gy-gravity [m/s^2]
    gy = randomized ? rand(rgen, Ny1, Nx1)*2e-1.-1e-1 : zeros(Ny1, Nx1)
    # P nodes
    # density [kg/m^3]
    RHO = randomized ? rand(rgen, Ny1, Nx1)*1e4 : zeros(Ny1, Nx1)
    # volumetric heat capacity [J/m^3/K]
    RHOCP = randomized ? rand(rgen, Ny1, Nx1)*1e6 : zeros(Ny1, Nx1)
    # thermal expansion [J/m^3/K]
    ALPHA = randomized ? rand(rgen, Ny1, Nx1)*1e-4 : zeros(Ny1, Nx1)
    # fluid thermal expansion [J/m^3/K]
    ALPHAF = randomized ? rand(rgen, Ny1, Nx1)*1e-4 : zeros(Ny1, Nx1)
    # radioactive heating [W/m^3]
    HR = randomized ? rand(rgen, Ny1, Nx1)*1e-3 : zeros(Ny1, Nx1)
    # adiabatic heating [W/m^3]
    HA = randomized ? rand(rgen, Ny1, Nx1)*1e-10 : zeros(Ny1, Nx1)
    # shear heating [W/m^3]
    HS = randomized ? rand(rgen, Ny1, Nx1)*1e-9 : zeros(Ny1, Nx1)
    # viscosity [Pa*s]
    ETAP = randomized ? rand(rgen, Ny1, Nx1)*1e16 : zeros(Ny1, Nx1)
    # shear modulus [Pa]
    GGGP = randomized ? rand(rgen, Ny1, Nx1)*1e10 : zeros(Ny1, Nx1)
    # EPSILONxx [1/s]
    EXX = randomized ? rand(rgen, Ny1, Nx1)*2e-12.-1e-12 : zeros(Ny1, Nx1)
    # σ′xx [1/s]
    SXX = randomized ? rand(rgen, Ny1, Nx1)*2e3-1e3 : zeros(Ny1, Nx1)
    # σ₀′ (SIGMA0'xx) [1/s]
    SXX0 = randomized ? rand(rgen, Ny1, Nx1)*2e3-1e3 : zeros(Ny1, Nx1)
    # current temperature [K]
    tk1 = randomized ? rand(rgen, Ny1, Nx1)*1e3 : zeros(Ny1, Nx1)
    # next temperature [K]
    tk2 = randomized ? rand(rgen, Ny1, Nx1)*1e3 : zeros(Ny1, Nx1)
    # temperature difference at P nodes [K]
    DT = randomized ? rand(rgen, Ny1, Nx1)*2e2.-1e2 : zeros(Ny1, Nx1)
    # previous temperature difference at P nodes [K]
    DT0 = randomized ? rand(rgen, Ny1, Nx1)*2e2.-1e2 : zeros(Ny1, Nx1)
    # solid vx in pressure nodes [m/s]
    vxp = randomized ? rand(rgen, Ny1, Nx1)*2e-9.-1e-9 : zeros(Ny1, Nx1)
    # solid vy in pressure nodes [m/s]
    vyp = randomized ? rand(rgen, Ny1, Nx1)*2e-9.-1e-9 : zeros(Ny1, Nx1)
    # fluid vx in pressure nodes [m/s]
    vxpf = randomized ? rand(rgen, Ny1, Nx1)*2e-9.-1e-9 : zeros(Ny1, Nx1)
    # fluid vy in pressure nodes [m/s]
    vypf = randomized ? rand(rgen, Ny1, Nx1)*2e-9.-1e-9 : zeros(Ny1, Nx1)
    # total pressure [Pa]
    pr = randomized ? rand(rgen, Ny1, Nx1)*1e4 : zeros(Ny1, Nx1)
    # fluid pressure [Pa]
    pf = randomized ? rand(rgen, Ny1, Nx1)*1e4 : zeros(Ny1, Nx1)
    # solid pressure [Pa]
    ps = randomized ? rand(rgen, Ny1, Nx1)*1e4 : zeros(Ny1, Nx1)
    # previous total pressure [Pa]
    pr0 = randomized ? rand(rgen, Ny1, Nx1)*1e4 : zeros(Ny1, Nx1)
    # previous fluid pressure [Pa]
    pf0 = randomized ? rand(rgen, Ny1, Nx1)*1e4 : zeros(Ny1, Nx1)
    # previous solid pressure [Pa]
    ps0 = randomized ? rand(rgen, Ny1, Nx1)*1e4 : zeros(Ny1, Nx1)
    # bulk viscosity [Pa*s]
    ETAPHI = randomized ? rand(rgen, Ny1, Nx1)*1e14 : zeros(Ny1, Nx1)
    # bulk compressibility [Pa*s]
    BETAPHI = randomized ? rand(rgen, Ny1, Nx1)*1e-10 : zeros(Ny1, Nx1)
    # porosity
    PHI = randomized ? rand(rgen, Ny1, Nx1) : zeros(Ny1, Nx1)
    # Dln[(1-ϕ)/ϕ]/Dt
    APHI = randomized ? rand(rgen, Ny1, Nx1)*2e-12.-1e-12 : zeros(Ny1, Nx1)
    # gravity potential [J/kg]
    FI = randomized ? rand(rgen, Ny1, Nx1)*2e2.=1e2 : zeros(Ny1, Nx1)
    # mass transfer term
    DMP = randomized ? rand(rgen, Ny1, Nx1) : zeros(Ny1, Nx1)
    # enthalpy transfer/latent heating term
    DHP = randomized ? rand(rgen, Ny1, Nx1) : zeros(Ny1, Nx1)
    # wet solid fraction
    XWS = randomized ? rand(rgen, Ny1, Nx1) : zeros(Ny1, Nx1)
    return (
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
    )
end # function setup_staggered_grid_properties()

"""
Set up additional helper staggered grid properties to facilitate computations.

$(SIGNATURES)

# Details

    - randomized: fill in random values for grid properties instead of zeros

# Returns

    - ETA5: plastic iterations viscoplastic viscosity at basic nodes [Pa⋅s]
    - ETA00: previous viscous viscosity at basic nodes [Pa⋅s]
    - YNY5: plastic iterations plastic yielding node property at basic nodes
    - YNY00: previous plastic yielding node property at basic nodes
    - DSXY: stress change Δσxy at basic nodes [Pa]
    - DSY: (SIIB-syield) at basic nodes
    - EII :second strain rate invariant at P nodes [1/s]
    - SII :second stress invariant at P nodes [Pa]
    - DSXX :stress change Δσ′xx at P nodes [Pa]
    - tk0: previous temperature at P nodes [K]
"""
function setup_staggered_grid_properties_helpers(;randomized=false)
    # basic nodes
    # plastic iterations viscoplastic viscosity at basic nodes [Pa⋅s]
    ETA5 = randomized ? rand(rgen, Ny, Nx)*1e16 : zeros(Ny, Nx)
    # previous viscous viscosity at basic nodes [Pa⋅s]
    ETA00 = randomized ? rand(rgen, Ny, Nx)*1e16 : zeros(Ny, Nx)
    # plastic iterations plastic yielding node property at basic nodes
    YNY5 = randomized ? rand(rgen, Bool, Ny, Nx) : zeros(Bool, Ny, Nx)
    # previous plastic yielding node property at basic nodes
    YNY00 = randomized ? rand(rgen, Bool, Ny, Nx) : zeros(Bool, Ny, Nx)
    # inverse viscoplastic viscosity at yielding basic nodes [1/(Pa⋅s)]
    YNY_inv_ETA = randomized ? rand(rgen, Ny, Nx)*1e-16 : zeros(Ny, Nx)
    # stress change Δσxy at basic nodes [Pa]
    DSXY = randomized ? rand(rgen, Ny, Nx)*2e3.-1e3 : zeros(Ny, Nx)
    # (SIIB-syield) at basic nodes
    DSY = randomized ? rand(rgen, Ny, Nx)*2e3.-1e3 : zeros(Ny, Nx)
    # second strain rate invariant at P nodes [1/s]
    EII = randomized ? rand(rgen, Ny1, Nx1)*1e-12 : zeros(Ny1, Nx1)
    # second stress invariant at P nodes [Pa]
    SII = randomized ? rand(rgen, Ny1, Nx1)*1e3 : zeros(Ny1, Nx1)
    # stress change Δσ′xx at P nodes [Pa]
    DSXX = randomized ? rand(rgen, Ny1, Nx1)*2e3.-1e3 : zeros(Ny1, Nx1)
    # previous temperature at P nodes [K]
    tk0 = randomized ? rand(rgen, Ny1, Nx1)*1e3 : zeros(Ny1, Nx1)
    return (
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
    )
end # function setup_staggered_grid_properties_helpers()

"""
Set up geodesic and physical properties of the set of markers.

$(SIGNATURES)

# Details

    - randomized: fill in random values for grid properties instead of zeros

# Returns

    - xm : horizontal marker coordinate [m]
    - ym : vertical marker coordinate [m]
    - tm : marker material type
    - tkm : marker temperature [K]
    - sxxm : marker σ′xx [Pa]
    - sxym : marker σxy [Pa]
    - etavpm : marker viscoplastic viscosity [Pa]
    - phim : marker porosity
"""
function setup_marker_properties(marknum; randomized=false)
    # horizontal marker coordinate [m]
    xm = randomized ? rand(rgen, -dx:0.1:xsize+dx, marknum) : zeros(marknum)
    # vertical marker coordinate [m]
    ym = randomized ? rand(rgen, -dy:0.1:ysize+dy, marknum) : zeros(marknum)
    # marker material type
    tm = randomized ? rand(rgen, 1:3, marknum) : zeros(Int, marknum)
    # marker temperature [K]
    tkm = randomized ? rand(rgen, 273:300, marknum) : zeros(marknum)
    # marker σ′xx [Pa]
    sxxm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker σxy [Pa]
    sxym = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker viscoplastic viscosity [Pa]
    etavpm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker porosity
    phim = randomized ? rand(rgen, marknum) : zeros(marknum)
    # reacted marker porosity
    phinewm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # previous marker fluid pressure
    pfm0 = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker melt molar fraction
    XWsolidm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # previous marker melt molar fraction
    XWsolidm0 = randomized ? rand(rgen, marknum) : zeros(marknum)
    return (
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
    )
end # function setup_marker_properties()

"""
Set up additional helper marker properties to facility comptuations.

$(SIGNATURES)

# Details

    - randomized: fill in random values for grid properties instead of zeros

# Returns

    - rhototalm: total density of markers
    - rhocptotalm : total volumetric heat capacity of markers
    - etatotalm: total viscosity of markers
    - hrtotalm: total radiogenic heat production of markers
    - ktotalm: total thermal conductivity of markers
    - tkm_rhocptotalm: total thermal energy of markers
    - etafluidcur_inv_kphim: fluid viscosity over permeability of markers
    - inv_gggtotalm: inverse of total shear modulus of markers
    - fricttotalm: total friction coefficient of markers
    - cohestotalm: total compressive strength of markers
    - tenstotalm: total tensile strength of markers
    - rhofluidcur: fluid density of markers
    - alphasolidcur: solid thermal expansion coefficient of markers
    - alphafluidcur: fluid thermal expansion coefficient of markers
"""
function setup_marker_properties_helpers(marknum; randomized=false)
    # marker total density
    rhototalm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker total volumetric heat capacity
    rhocptotalm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker total viscosity
    etatotalm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker total radiogenic heat production
    hrtotalm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker total thermal conductivity
    ktotalm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker total thermal energy
    tkm_rhocptotalm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker fluid viscosity over permeability
    etafluidcur_inv_kphim = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker inverse of total shear modulus
    inv_gggtotalm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker total friction coefficient
    fricttotalm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker total compressive strength
    cohestotalm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker total tensile strength
    tenstotalm = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker fluid density
    rhofluidcur = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker solid thermal expansion coefficient
    alphasolidcur = randomized ? rand(rgen, marknum) : zeros(marknum)
    # marker fluid thermal expansion coefficient
    alphafluidcur = randomized ? rand(rgen, marknum) : zeros(marknum)
    return (
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
end # function setup_marker_properties_helpers()

"""
Set up additional marker geometry helpers to facilitate marker handling.

$(SIGNATURES)

# Details

    - nothing

# Returns

    - mdis: minimum distance of marker launch anchor points to nearest marker
    - mnum: number of marker nearest to marker launch anchor positions
"""
function setup_marker_geometry_helpers()
    mdis = fill(mdis_init, Nym, Nxm)
    mnum = zeros(Int, Nym, Nxm)
    return mdis, mnum
end

"""
Define initial set of markers according to model parameters

$(SIGNATURES)

# Details

    - xm: x coordinates of markers
    - ym: y coordinates of markers
    - tm: material type of markers
    - phim: porosity of markers
    - etavpm: matrix viscosity of markers
    - rhototalm: total density of markers
    - rhocptotalm: total volumetric heat capacity of markers
    - etatotalm: total viscosity of markers
    - hrtotalm: total radiogenic heat production of markers
    - ktotalm: total thermal conductivity of markers
    - tkm: temperature of markers 
    - inv_gggtotalm: inverse of total shear modulus of markers
    - fricttotalm: total friction coefficient of markers
    - cohestotalm: total compressive strength of markers
    - tenstotalm: total tensile strength of markers
    - rhofluidcur: fluid density of markers
    - alphasolidcur: solid thermal expansion coefficient of markers
    - alphafluidcur: fluid thermal expansion coefficient of markers
    - XWsolidm0: previous wet solid molar fraction of markers
    - randomized: uniformly random-distribute marker x/y positions within cells
                  and randomly set initial marker porosity 

# Returns

    - nothing
"""
function define_markers!(
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
    randomized=random_markers
)
    for jm=1:1:Nxm, im=1:1:Nym
        # calculate marker counter
        m = (jm-1) * Nym + im
        # define marker coordinates
        xm[m] = dxm/2 + (jm-1) * dxm 
        ym[m] = dym/2 + (im-1) * dym 
        # random marker position within cell
        if randomized
            xm[m] += (rand(rgen)-0.5) * dxm
            ym[m] += (rand(rgen)-0.5) * dym
        end
        # primary marker properties 
        rmark = distance(xm[m], ym[m], xcenter, ycenter)
        if rmark < rplanet
            # planet
            tm[m] = ifelse(rmark>rcrust, 2, 1)
            # porosity
            phim[m] = phim0
            if randomized
                phim[m] += phim0 * (rand(rgen)-0.5)
            end
            # matrix viscosity
            etavpm[m] = etasolidm[tm[m]] # * exp(-αη*phim[m]) # ∇! CHANGE!!!
            # wet solid molar fraction
            XWsolidm0[m] = XWsolidm_init[tm[m]]
            if randomized
                XWsolidm0[m] += XWsolidm_init[tm[m]] * (rand(rgen)-0.5)
            end
        else
            # sticky space ("air") [to have internal free surface]
            tm[m] = 3
            # porosity
            phim[m] = phimin
            # matrix viscosity
            etavpm[m] = etasolidm[tm[m]]
            # static properties for air markers
            rhototalm[m] = rhosolidm[tm[m]]
            rhocptotalm[m] = rhocpsolidm[tm[m]]
            etatotalm[m] = etasolidm[tm[m]]
            hrtotalm[m] = start_hrsolidm[tm[m]]
            ktotalm[m] = ksolidm[tm[m]]           
        end
        # common initialisations for all marker types
        tkm[m] = tkm0[tm[m]]
        inv_gggtotalm[m] = inv(gggsolidm[tm[m]])
        fricttotalm[m] = frictsolidm[tm[m]]
        cohestotalm[m] = cohessolidm[tm[m]]
        tenstotalm[m] = tenssolidm[tm[m]]
        rhofluidcur[m] = rhofluidm[tm[m]]
        alphasolidcur[m] = alphasolidm[tm[m]]
        alphafluidcur[m] = alphafluidm[tm[m]]
    end
    return nothing
end

"""
Compute properties of given marker and save them to corresponding arrays.

$(SIGNATURES)

# Details

    - m: marker number
    - tm: type of markers
    - tkm: temperature of markers
    - rhototalm: total density of markers
    - rhocptotalm: total volumetric heat capacity of markers
    - etatotalm: total viscosity of markers
    - hrtotalm: total radiogenic heat production of markers
    - ktotalm: total thermal conductivity of markers
    - tkm_rhocptotalm: total thermal energy of markers
    - etafluidcur_inv_kphim: (fluid viscosity)/permeability of markers
    - hrsolidm: vector of radiogenic heat production of solid materials
    - hrfluidm: vector of radiogenic heat production of fluid materials
    - phim: porosity of markers
    - XWˢm₀: previous wet solid molar fraction of markers
    - mode: marker property computation mode
        - 1: dynamic, based on (Touloukian, 1970; Hobbs, 1974;
             Travis and Schubert, 2005)
        - 2: constant parameter rhocpfluidm

# Returns

    - nothing
"""
function compute_marker_properties!(
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
    XWˢm₀,
    mode
)
@timeit to "compute_marker_properties!" begin
    if tm[m] < 3
        # rocks
        XDˢm₀ = 1.0 - XWˢm₀[m]
        rhosolidm0 = (MD + MH₂O*XWˢm₀[m]) / (VDˢ*XDˢm₀+ VWˢ*XWˢm₀[m]) # (16.161)
        rhofluidm0 = ifelse(
            tkm[m]>tmfluidphase, ρH₂Oᶠ, ρH₂Oᶠⁱ) # (16.162)
        rhototalm[m] = total(rhosolidm0, rhofluidm0, phim[m])
        rhocptotalm[m] = total(
            rhocpsolidm[tm[m]], compute_rhocpfluidm(tkm[m], mode), phim[m])
        etasolidcur = ifelse(
            tkm[m]>tmsolidphase, etasolidmm[tm[m]], etasolidm[tm[m]])
        etafluidcur = ifelse(
            tkm[m]>tmfluidphase, etafluidmm[tm[m]], etafluidm[tm[m]])
        etatotalm[m] = max(etamin, etasolidcur, etafluidcur)
        hrtotalm[m] = total(hrsolidm[tm[m]], hrfluidm[tm[m]], phim[m])
        ktotalm[m] = ktotal(
            compute_ksolidm(tkm[m], mode),
            compute_kfluidm(tkm[m], mode),
            phim[m]
        )
    else
        # sticky air
        etafluidcur = etafluidm[tm[m]]
    end
    # common for rocks and air
    tkm_rhocptotalm[m] = tkm[m] * rhocptotalm[m]
    # kphim[m] = kphi(kphim0[tm[m]], phim[m])
    # etafluidcur_inv_kphim[m] = etafluidcur[m] * inv(kphim[m])
    etafluidcur_inv_kphim[m] = ηᶠcur_inv_kᵠ(
        kphim0[tm[m]], phim[m], etafluidcur)
end # @timeit to "compute_marker_properties!"
    return nothing
end # function compute_marker_properties!

"""
Update marker viscoplastic viscosity based on basic node yielding status and
marker temperature- and material-based viscosity.

$(SIGNATURES)

# Details

    - m: marker number
    - xm: x-coordinates of markers
    - ym: y-coordinates of markers
    - tm: type of markers
    - tkm: temperature of markers
    - etatotalm: total viscosity of markers
    - etavpm: matrix viscosity of markers
    - YNY: plastic yielding node property at basic nodes
    - YNY_inv_ETA: inverse viscoplastic viscosity at yielding basic nodes

# Returns

    -nothing
"""
function update_marker_viscosity!(
    m, xm, ym, tm, tkm, etatotalm, etavpm, YNY, YNY_inv_ETA)
    @inbounds i, j, weights = fix_weights(
        xm[m],
        ym[m],
        x,
        y,
        dx,
        dy,
        jmin_basic,
        jmax_basic,
        imin_basic,
        imax_basic
    )
    @inbounds if tm[m] < 3
        # rocks: update etatotalm[m] based on current marker temperature
        @inbounds etatotalm[m] = etatotal_rocks(tkm[m], tm[m]) # * exp(-αη*phim[m]) # ∇! CHANGE!!!
    # else
        # air: constant etatotalm[m]=etasolidm[tm[m]] as initialized
        # pass
    end
    if any(grid_vector(i, j, YNY))
        interpolate_to_marker!(m, i, j, weights, etavpm, YNY_inv_ETA)
        @inbounds etavpm[m] = inv(etavpm[m])
        @inbounds etavpm[m] = ifelse(
            etavpm[m]>etatotalm[m], etatotalm[m], etavpm[m])
    else
        @inbounds etavpm[m] = etatotalm[m]
    end
    return nothing
end

"""
Set up properties to be interpolated from markers to staggered grid.

$(SIGNATURES)

# Details

    - nothing

# Returns

    - ETA0SUM: interpolation of ETA0 at basic nodes
    - ETASUM: interpolation of ETA at basic nodes
    - GGGSUM: interpolation of GGG at basic nodes
    - SXYSUM: interpolation of SXY at basic nodes
    - COHSUM: interpolation of COH at basic nodes
    - TENSUM: interpolation of TEN at basic nodes
    - FRISUM: interpolation of FRI at basic nodes
    - WTSUM: interpolation weights at basic nodes
    - RHOXSUM: interpolation of RHOX at Vx nodes
    - RHOFXSUM: interpolation of RHOFX at Vx nodes
    - KXSUM: interpolation of KX at Vx nodes
    - PHIXSUM: interpolation of PHIX at Vx nodes
    - RXSUM: interpolation of RX at Vx nodes
    - WTXSUM: interpolation weights at Vx nodes
    - RHOYSUM: interpolation of RHOY at Vy nodes
    - RHOFYSUM: interpolation of RHOFY at Vy nodes
    - KYSUM: interpolation of KY at Vy nodes
    - PHIYSUM: interpolation of PHIX at Vy nodes
    - RYSUM: interpolation of RY at Vy nodes
    - WTYSUM: interpolation weights at Vy nodes
    - RHOSUM: interpolation of RHO at P nodes
    - RHOCPSUM: interpolation of RHOCP at P nodes
    - ALPHASUM: interpolation of ALPHA at P nodes
    - ALPHAFSUM: interpolation of ALPHAF at P nodes
    - HRSUM: interpolation of HR at P nodes
    - GGGPSUM: interpolation of GGGP at P nodes
    - SXXSUM: interpolation of SXX at P nodes
    - TKSUM: interpolation of TK at P nodes
    - PHISUM: interpolation of PHI at P nodes
    - DMPSUM: interpolation of DMP at P nodes
    - DHPSUM: interpolation of DHP at P nodes
    - XWSSUM: interpolation of XWS at P nodes
    - WTPSUM: interpolation weights at P nodes
"""
function setup_interpolated_properties()
    # basic nodes
    ETA0SUM = zeros(Ny, Nx)
    ETASUM = zeros(Ny, Nx)
    GGGSUM = zeros(Ny, Nx)
    SXYSUM = zeros(Ny, Nx)
    COHSUM = zeros(Ny, Nx)
    TENSUM = zeros(Ny, Nx)
    FRISUM = zeros(Ny, Nx)
    WTSUM = zeros(Ny, Nx)
    # Vx nodes
    RHOXSUM = zeros(Ny1, Nx1)
    RHOFXSUM = zeros(Ny1, Nx1)
    KXSUM = zeros(Ny1, Nx1)
    PHIXSUM = zeros(Ny1, Nx1)
    RXSUM = zeros(Ny1, Nx1)
    WTXSUM = zeros(Ny1, Nx1)
    # Vy nodes
    RHOYSUM = zeros(Ny1, Nx1)
    RHOFYSUM = zeros(Ny1, Nx1)
    KYSUM = zeros(Ny1, Nx1)
    PHIYSUM = zeros(Ny1, Nx1)
    RYSUM = zeros(Ny1, Nx1)
    WTYSUM = zeros(Ny1, Nx1)
    # P Nodes
    RHOSUM = zeros(Ny1, Nx1)
    RHOCPSUM = zeros(Ny1, Nx1)
    ALPHASUM = zeros(Ny1, Nx1)
    ALPHAFSUM = zeros(Ny1, Nx1)
    HRSUM = zeros(Ny1, Nx1)
    GGGPSUM = zeros(Ny1, Nx1)
    SXXSUM = zeros(Ny1, Nx1)
    TKSUM = zeros(Ny1, Nx1)
    PHISUM = zeros(Ny1, Nx1)
    DMPSUM = zeros(Ny1, Nx1)
    DHPSUM = zeros(Ny1, Nx1)
    XWSSUM = zeros(Ny1, Nx1)
    WTPSUM = zeros(Ny1, Nx1)
    return (
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
    )
end

"""
Reset properties to be interpolated from markers to staggered grid.

$(SIGNATURES)

# Details   

    - ETA0SUM: interpolation of ETA0 at basic nodes
    - ETASUM: interpolation of ETA at basic nodes
    - GGGSUM: interpolation of GGG at basic nodes
    - SXYSUM: interpolation of SXY at basic nodes
    - COHSUM: interpolation of COH at basic nodes
    - TENSUM: interpolation of TEN at basic nodes
    - FRISUM: interpolation of FRI at basic nodes
    - WTSUM: interpolation weights at basic nodes
    - RHOXSUM: interpolation of RHOX at Vx nodes
    - RHOFXSUM: interpolation of RHOFX at Vx nodes
    - KXSUM: interpolation of KX at Vx nodes
    - PHIXSUM: interpolation of PHIX at Vx nodes
    - RXSUM: interpolation of RX at Vx nodes
    - WTXSUM: interpolation weights at Vx nodes
    - RHOYSUM: interpolation of RHOY at Vy nodes
    - RHOFYSUM: interpolation of RHOFY at Vy nodes
    - KYSUM: interpolation of KY at Vy nodes
    - PHIYSUM: interpolation of PHIX at Vy nodes
    - RYSUM: interpolation of RY at Vy nodes
    - WTYSUM: interpolation weights at Vy nodes
    - RHOSUM: interpolation of RHO at P nodes
    - RHOCPSUM: interpolation of RHOCP at P nodes
    - ALPHASUM: interpolation of ALPHA at P nodes
    - ALPHAFSUM: interpolation of ALPHAF at P nodes
    - HRSUM: interpolation of HR at P nodes
    - GGGPSUM: interpolation of GGGP at P nodes
    - SXXSUM: interpolation of SXX at P nodes
    - TKSUM: interpolation of TK at P nodes
    - PHISUM: interpolation of PHI at P nodes
    - WTPSUM: interpolation weights at P nodes

# Returns

    - nothing
"""
function reset_interpolated_properties!(
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
        # basic nodes
        ETA0SUM .= zero(0.0)
        ETASUM .= zero(0.0)
        GGGSUM .= zero(0.0)
        SXYSUM .= zero(0.0)
        COHSUM .= zero(0.0)
        TENSUM .= zero(0.0)
        FRISUM .= zero(0.0)
        WTSUM .= zero(0.0)
        # Vx nodes
        RHOXSUM .= zero(0.0)
        RHOFXSUM .= zero(0.0)
        KXSUM .= zero(0.0)
        PHIXSUM .= zero(0.0)
        RXSUM .= zero(0.0)
        WTXSUM .= zero(0.0)
        # Vy nodes
        RHOYSUM .= zero(0.0)
        RHOFYSUM .= zero(0.0)
        KYSUM .= zero(0.0)
        PHIYSUM .= zero(0.0)
        RYSUM .= zero(0.0)
        WTYSUM .= zero(0.0)
        # P Nodes
        RHOSUM .= zero(0.0)
        RHOCPSUM .= zero(0.0)
        ALPHASUM .= zero(0.0)
        ALPHAFSUM .= zero(0.0)
        HRSUM .= zero(0.0)
        GGGPSUM .= zero(0.0)
        SXXSUM .= zero(0.0)
        TKSUM .= zero(0.0)
        PHISUM .= zero(0.0)
        WTPSUM .= zero(0.0)
    return nothing
end

"""
Reset thermochemical properties to be interpolated to staggered grid.

$(SIGNATURES)

# Details   

    - DMPSUM: interpolation of DMP at P nodes
    - DHPSUM: interpolation of DHP at P nodes
    - WTPSUM: interpolation weights at P nodes
    
# Returns

    - nothing
"""
function reset_thermochemical_properties!(
    DMPSUM,
    DHPSUM,
    WTPSUM
)
        # P Nodes
        DMPSUM .= zero(0.0)
        DHPSUM .= zero(0.0)
        WTPSUM .= zero(0.0)
    return nothing
end

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
Calculate Euclidean distance between two point coordinates.

$(SIGNATURES)

# Details

    - x1: x-coordinate of point 1 [m]
    - y1: y-coordinate of point 1 [m]
    - x2: x-coordinate of point 2 [m]
    - y2: y-coordinate of point 2 [m]

# Returns

    - Euclidean distance between point 1 and point 2 [m]
"""
function distance(x1, y1, x2, y2)
    return sqrt(abs2(x1-x2) + abs2(y1-y2))
end

"""
Compute convex combination of fluid and solid properties to get total property.

$(SIGNATURES)

# Details

    - fluid: fluid properties
    - solid: solid properties
    - ϕ: porosity (fraction of fluid)

# Returns

    - total: computed total property
"""
function total(solid, fluid, ϕ)
    return solid*(1.0-ϕ) + fluid*ϕ
end

"""
Get a 4-vector of values from a grid 4-stencil anchored at (i, j) in 
column-major order.

$(SIGNATURES)

# Details

    - i: top left grid node column index
    - j: top left grid node row index
    - grid: data from which to build Vector

# Returns

    - grid_vector: 4-vector of values
        [grid[i, j], grid[i+1, j], grid[i, j+1], grid[i+1, j+1]]
"""
function grid_vector(i, j, grid)
    @inbounds return @SVector [
        grid[i, j], grid[i+1, j], grid[i, j+1], grid[i+1, j+1]
    ]
end

"""
Compute inner product of two 4-vectors of reals.

$(SIGNATURES)

# Details

    - v1: first 4-vector
    - v2: second 4-vector

# Returns

    - inner product of v1 and v2
"""
function dot4(v1, v2)
    @inbounds return v1[1]*v2[1] + v1[2]*v2[2] + v1[3]*v2[3] + v1[4]*v2[4]
end

"""
Get the arithmetic average of values from a grid 4-stencil anchored at (i, j).

$(SIGNATURES)

# Details

    - i: top left grid node column index
    - j: top left grid node row index
    - grid: data from which to get the average value

# Returns

    - grid_average: (grid[i, j]+grid[i+1, j]+grid[i, j+1]+grid[i+1, j+1]) / 4
"""
function grid_average(i, j, grid)
    # return sum(grid_vector(i, j, grid)) * inv(length(grid_vector(i, j, grid)))
    @inbounds return 0.25 * (
        grid[i, j]+grid[i+1, j]+grid[i, j+1]+grid[i+1, j+1])
end

"""
Add a RK4 stage velocity to the RK4 velocity vector.

$(SIGNATURES)

# Details:

    - vrk4: current RK4 velocity vector
    - v: RK4 velocity of stage `rk` to be added to velocity vector
    - rk: RK4 stage number

# Returns

    - vrk4: updated RK4 velocity vector
"""
function add_vrk4(vrk4, v, rk)
    if rk == 1
        return vrk4 + @SVector [v, 0.0, 0.0, 0.0]
    elseif rk == 2
        return vrk4 + @SVector [0.0, v, 0.0, 0.0]
    elseif rk == 3
        return vrk4 + @SVector [0.0, 0.0, v, 0.0]
    elseif rk == 4
        return vrk4 + @SVector [0.0, 0.0, 0.0, v]
    else
        return vrk4
    end
end

"""
Compute total thermal conductivity of two-phase material.

$(SIGNATURES)

# Details

    - ksolid: solid thermal conductivity [W/m/K]
    - kfluid: fluid thermal conductivity [W/m/K]
    - phi: fraction of solid

# Returns

    - ktotal: total thermal conductivity of mixed phase [W/m/K]
"""
function ktotal(ksolid, kfluid, phi)
    return (
        sqrt(
            ksolid * kfluid/2
            + ((ksolid*(3.0*phi-2.0) + kfluid*(1.0-3.0*phi))^2)*inv(16.0)
        )
        -0.25 * (ksolid*(3.0*phi-2.0) + kfluid*(1.0-3.0*phi))
    )
end


"""
Compute porosity-dependent permeability (eqn 16.64 in Gerya (2019)).

$(SIGNATURES)

# Details

    - kphim0m: standard (reference) permeability (of marker type) [m^2]
    - phimm: actual (marker) porosity

# Returns

    - kphim: empirical porosity-dependent permeability [m^2]
"""
function kphi(kphim0m, phimm)
    # phim0 is a global constant defined independent of material type
    return kphim0m * (phimm*inv(phim0))^3.0 * ((1.0-phimm)*inv(1.0-phim0))^-2.0
end

"""
Compute inverse of porosity-dependent permeability (eqn 16.64 in Gerya (2019)) 
times current fluid viscosity.

$(SIGNATURES)

# Details

    - kϕᵣ: reference permeability [m^2]
    - ϕ: current porosity
    - ηᶠcur: current fluid viscosity [Pa s]

# Returns

    - etafluidcur_inv_kphi: inverse empirical porosity-dependent permeability 
                            times current fluid viscosity
"""
function ηᶠcur_inv_kᵠ(kϕᵣ, ϕ, ηᶠcur)
    return ηᶠcur * inv(kϕᵣ) * (phim0*inv(ϕ))^3.0 * ((1.0-ϕ)*inv(1.0-phim0))^2.0
end

"""
Compute radiogenic heat production of isotope mixture.

$(SIGNATURES)

# Details

    - f: fraction of radioactive matter [atoms/kg]
    - ratio: initial ratio of radioactive to non-radioactive isotopes
    - E: heat energy [J]
    - tau: exp decay mean lifetime ``\\tau=\\frac{t_{1/2}}{\\log{2}}`` [s]
    - time: time elapsed since start of radioactive decay [s]

# Returns

    - Q: radiogenic heat production [W/kg]
"""
function Q_radiogenic(f, ratio, E, tau, time)
    return f * ratio * E * exp(-time*inv(tau)) * inv(tau)
end

"""
Compute total rocky marker viscosity based on temperature and material type.

$(SIGNATURES)

# Details

    - tkmm: marker temperature [K]
    - tmm: marker type [1, 2]

# Returns
    
    - etatotal: rocky marker temperature-dependent total viscosity 
"""
function etatotal_rocks(tkmm, tmm)
    @inbounds etasolidcur = ifelse(
        tkmm>tmsolidphase, etasolidmm[tmm], etasolidm[tmm])
    @inbounds etafluidcur = ifelse(
        tkmm>tmfluidphase, etafluidmm[tmm], etafluidm[tmm])
    return max(etamin, etasolidcur, etafluidcur)
end

"""
Compute radiogenic heat production of 26Al and 60Fe isotopes.

$(SIGNATURES)

# Details

    - al: true if radioactive isotope 26Al is present
    - fe: true if radioactive isotope 60Fe is present
    - timesum: time elapsed since initial conditions at start of simulation

# Returns

    - hrsolidm: radiogenic heat production of 26Al [W/m^3]
    - hrfluidm: radiogenic heat production of 60Fe [W/m^3]
"""
function calculate_radioactive_heating(al, fe, timesum)
    #26Al: planet ✓, crust ✓, space ×
    if al
        # 26Al radiogenic heat production [W/kg]
        Q_al = Q_radiogenic(f_al, ratio_al, E_al, tau_al, timesum)
        # Solid phase 26Al radiogenic heat production [W/m^3]
        @inbounds hrsolidm = @SVector [
            Q_al*rhosolidm[1], Q_al*rhosolidm[2], 0.0]
    else
        hrsolidm = @SVector zeros(3)
    end    
    #60Fe: planet ✓, crust ×, space ×
    if fe
        # 60Fe radiogenic heat production [W/kg]
        Q_fe = Q_radiogenic(f_fe, ratio_fe, E_fe, tau_fe, timesum)
        # Fluid phase 60Fe radiogenic heat production [W/m^3]
        @inbounds hrfluidm = @SVector [Q_fe*rhofluidm[1], 0.0, 0.0]
    else
        hrfluidm = @SVector zeros(3)
    end
    return hrsolidm, hrfluidm
end

"""
Compute top and left grid nodes indices and bilinear interpolation weigths to
nearest four grid nodes for given (x, y) position and grid axes.

$(SIGNATURES)

# Details

    - x: x-position [m]
    - y: y-position [m]
    - x_axis: x-grid reference axis array [m]
    - y_axis: y-grid reference axis array [m]
    - dx: x-grid axis mesh width [m]
    - dy: y-grid axis mesh width [m]
    - jmin: minimum assignable index on x-grid axis (basic/Vx/Vy/P)
    - jmax: maximum assignable index on x-grid axis (basic/Vx/Vy/P)
    - imin: minimum assignable index on y-grid axis (basic/Vx/Vy/P)
    - imax: maximum assignable index on y-grid axis (basic/Vx/Vy/P)

# Returns
    - i: top (with reference to y) node index on y-grid axis
    - j: left (with reference to x) node index on x-grid axis
    - bilinear_weights: vector of 4 bilinear interpolation weights to
      nearest four grid nodes:
        [wtmij  : i  , j   node,
        wtmi1j : i+1, j   node,
        wtmij1 : i  , j+1 node,
        wtmi1j1: i+1, j+1 node]
"""
function fix_weights(x, y, x_axis, y_axis, dx, dy, jmin, jmax, imin, imax)
    i, j, dxmj, dymi = fix_distances(
        x, y, x_axis, y_axis, dx, dy, jmin, jmax, imin, imax)
    return i, j, SVector(
            (1.0-dymi/dy) * (1.0-dxmj/dx),
            (dymi/dy) * (1.0-dxmj/dx),
            (1.0-dymi/dy) * (dxmj/dx),
            (dymi/dy) * (dxmj/dx)
        ) 
end # function fix_weights

"""
Compute top and left grid nodes indices (i, j) and x- and y-distances to that 
grid node (i, j) for given (x, y) position and grid axes.

$(SIGNATURES)

# Details

    - x: x-position [m]
    - y: y-position [m]
    - x_axis: x-grid reference axis array [m]
    - y_axis: y-grid reference axis array [m]
    - dx: x-grid axis mesh width [m]
    - dy: y-grid axis mesh width [m]
    - jmin: minimum assignable index on x-grid axis (basic/Vx/Vy/P)
    - jmax: maximum assignable index on x-grid axis (basic/Vx/Vy/P)
    - imin: minimum assignable index on y-grid axis (basic/Vx/Vy/P)
    - imax: maximum assignable index on y-grid axis (basic/Vx/Vy/P)

# Returns
    - i: top (with reference to y) node index on y-grid axis
    - j: left (with reference to x) node index on x-grid axis
    - dxmj: x-distance from (x, y) point to (i, j) node
    - dymi: y-distance from (x, y) point to (i, j) node
"""
function fix_distances(x, y, x_axis, y_axis, dx, dy, jmin, jmax, imin, imax)
    @inbounds begin
        i, j = fix(x, y, x_axis, y_axis, dx, dy, jmin, jmax, imin, imax)
        dxmj = x - x_axis[j]
        dymi = y - y_axis[i]
    end # @inbounds
    return i, j, dxmj, dymi
end # function fix_distances

"""
Compute top and left grid nodes indices (i, j) and x- and y-distances to that 
grid node (i, j) for given (x, y) position and grid axes.

$(SIGNATURES)

# Details

    - x: x-position [m]
    - y: y-position [m]
    - x_axis: x-grid reference axis array [m]
    - y_axis: y-grid reference axis array [m]
    - dx: x-grid axis mesh width [m]
    - dy: y-grid axis mesh width [m]
    - jmin: minimum assignable index on x-grid axis (basic/Vx/Vy/P)
    - jmax: maximum assignable index on x-grid axis (basic/Vx/Vy/P)
    - imin: minimum assignable index on y-grid axis (basic/Vx/Vy/P)
    - imax: maximum assignable index on y-grid axis (basic/Vx/Vy/P)

# Returns
    - i: top (with reference to y) node index on y-grid axis
    - j: left (with reference to x) node index on x-grid axis
"""
function fix(x, y, x_axis, y_axis, dx, dy, jmin, jmax, imin, imax)
    @inbounds j = unsafe_trunc(Int, (x-x_axis[1])*inv(dx)) + 1
    @inbounds i = unsafe_trunc(Int, (y-y_axis[1])*inv(dy)) + 1
    if j < jmin
        j = jmin
    elseif j > jmax
        j = jmax
    end
    if i < imin
        i = imin
    elseif i > imax
        i = imax
    end
    return i, j
end # function fix

"""
Reduce a 3D (i, j, k) array along its third (k) axis by addition and write the result
into (i, j, 1) without reallocating the array's memory.

$(SIGNATURES)

# Details

    - A: 3D array [i, j, k]

# Returns

    - nothing
"""
function reduce_add_3darray!(A)
    ii = axes(A, 1)
    ij = axes(A, 2)
    ik = axes(A, 3)
    for k in ik[(begin+1):end], j in ij, i in ii 
        @inbounds A[i, j, 1] += A[i, j, k]
    end
    return nothing
end

"""
Interpolate a property to neareast four nodes on a given grid location
using given bilinear interpolation weights.

# Details

    - i: top (with reference to y) node index on vertical y-grid axis
    - j: left (with reference to x) node index on horizontal x-grid axis
    - weights: vector of 4 bilinear interpolation weights to
      nearest four grid nodes:
        [wtmij  : i  , j   node,
        wtmi1j : i+1, j   node,
        wtmij1 : i  , j+1 node,
        wtmi1j1: i+1, j+1 node]
    - property: property to be interpolated to grid using weights
    - grid: threaded grid array on which to interpolate property

# Returns

    - nothing
"""
function interpolate_add_to_grid!(i, j, weights, property, grid)
@timeit to "interpolate_add_to_grid!" begin
    @inbounds grid[i, j] += property * weights[1]
    @inbounds grid[i+1, j] += property * weights[2]
    @inbounds grid[i, j+1] += property * weights[3]
    @inbounds grid[i+1, j+1] += property * weights[4]
end # @timeit to "interpolate_add_to_grid!"
    return nothing
end # function interpolate_add_to_grid!

"""
Interpolate a property from nearest the four nodes on a given grid to a marker.

# Details

    - m: number of marker to interpolate to
    - i: top (with reference to y) node index on vertical y-grid axis
    - j: left (with reference to x) node index on horizontal x-grid axis
    - weights: vector of 4 bilinear interpolation weights to
    nearest four grid nodes:
        [wtmij  : i  , j   node,
        wtmi1j : i+1, j   node,
        wtmij1 : i  , j+1 node,
        wtmi1j1: i+1, j+1 node]
    - marker_property: marker property array into which to interpolate
    - property_grid: grid whose property is to be interpolated to marker
    - m: marker 
"""
function interpolate_to_marker!(m, i, j, weights, marker_property, grid)
    @inbounds marker_property[m] = dot4(grid_vector(i, j, grid), weights)
    return nothing
end # function interpolate_to_marker


"""
Interpolate a property from nearest the four nodes on a given grid to a marker
and add it to the markers property.

# Details

    - m: number of marker to interpolate to
    - i: top (with reference to y) node index on vertical y-grid axis
    - j: left (with reference to x) node index on horizontal x-grid axis
    - weights: vector of 4 bilinear interpolation weights to
    nearest four grid nodes:
        [wtmij  : i  , j   node,
        wtmi1j : i+1, j   node,
        wtmij1 : i  , j+1 node,
        wtmi1j1: i+1, j+1 node]
    - marker_property: marker property array into which to interpolate and add
    - property_grid: grid whose property is to be interpolated to marker
    - m: marker 
"""
function interpolate_add_to_marker!(m, i, j, weights, marker_property, grid)
    @inbounds marker_property[m] += dot4(grid_vector(i, j, grid), weights)
    return nothing
end # function interpolate_add_to_marker

"""
Interpolate selected marker properties to basic nodes.

$(SIGNATURES)

# Details

    - m: marker number
    - xmm: marker's x-position [m]
    - ymm: marker's y-position [m]
    - etatotalm: total viscosity of markers
    - etavpm: matrix viscosity of markers
    - inv_gggtotalm: inverse of total shear modulus of markers
    - sxym: marker σxy [Pa]
    - cohestotalm: total compressive strength of markers
    - tenstotalm: total tensile strength of markers 
    - fricttotalm: total friction coefficient of markers
    - ETA0SUM: viscous viscosity interpolated to basic nodes
    - ETASUM: viscoplastic viscosity interpolated to basic nodes
    - GGGSUM: shear modulus interpolated to basic nodes
    - SXYSUM: σxy shear stress interpolated to basic nodes 
    - COHSUM: compressive strength interpolated to basic nodes
    - TENSUM: tensile strength interpolated to basic nodes
    - FRISUM: friction interpolated to basic nodes 
    - WTSUM: weight array for bilinear interpolation to basic nodes

# Returns

    - nothing
"""
function marker_to_basic_nodes!(
    m,
    xmm,
    ymm,
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
    i, j, weights = fix_weights(
        xmm,
        ymm,
        x,
        y,
        dx,
        dy,
        jmin_basic,
        jmax_basic,
        imin_basic,
        imax_basic
    )
    @inbounds begin
        interpolate_add_to_grid!(i, j, weights, etatotalm[m], ETA0SUM)
        interpolate_add_to_grid!(i, j, weights, etavpm[m], ETASUM)
        interpolate_add_to_grid!(i, j, weights, inv_gggtotalm[m], GGGSUM)
        interpolate_add_to_grid!(i, j, weights, sxym[m], SXYSUM)
        interpolate_add_to_grid!(i, j, weights, cohestotalm[m], COHSUM)
        interpolate_add_to_grid!(i, j, weights, tenstotalm[m], TENSUM)
        interpolate_add_to_grid!(i, j, weights, fricttotalm[m], FRISUM)
        interpolate_add_to_grid!(i, j, weights, one(1.0), WTSUM)
    end # @inbounds
    return nothing
end

"""
Interpolate selected marker properties to Vx nodes.

$(SIGNATURES)

# Details

    - m: marker number
    - xmm: marker's x-position [m]
    - ymm: marker's y-position [m]
    - rhototalm: total density of markers
    - rhofluidcur: fluid density of markers
    - ktotalm: total thermal conductivity of markers
    - phim: porosity of markers
    - etafluidcur_inv_kphim: fluid viscosity over permeability of markers
    - RHOXSUM: density interpolated to Vx nodes
    - RHOFXSUM: fluid density interpolated to Vx nodes
    - KXSUM: thermal conductivity interpolated to Vx nodes
    - PHIXSUM: porosity interpolated to Vx nodes
    - RXSUM: ηfluid/kϕ interpolated to Vx nodes
    - WTXSUM: weight for bilinear interpolation to Vx nodes

# Returns

    - nothing
"""
function marker_to_vx_nodes!(
    m,
    xmm,
    ymm,
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
    i, j, weights = fix_weights(
        xmm,
        ymm,
        xvx,
        yvx,
        dx,
        dy,
        jmin_vx,
        jmax_vx,
        imin_vx,
        imax_vx
    )
    @inbounds begin
        interpolate_add_to_grid!(i, j, weights, rhototalm[m], RHOXSUM)
        interpolate_add_to_grid!(i, j, weights, rhofluidcur[m], RHOFXSUM)
        interpolate_add_to_grid!(i, j, weights, ktotalm[m], KXSUM)
        interpolate_add_to_grid!(i, j, weights, phim[m], PHIXSUM)
        interpolate_add_to_grid!(i, j, weights, etafluidcur_inv_kphim[m], RXSUM)
        interpolate_add_to_grid!(i, j, weights, one(1.0), WTXSUM)
    end # @inbounds
    return nothing
end

"""
Interpolate selected marker properties to Vy nodes.

$(SIGNATURES)

# Details

    - m: marker number
    - xmm: marker's x-position [m]
    - ymm: marker's y-position [m]
    - rhototalm: total density of markers
    - rhofluidcur: fluid density of markers
    - ktotalm: total thermal conductivity of markers
    - phim: porosity of markers
    - etafluidcur_inv_kphim: fluid viscosity over permeability of markers
    - RHOYSUM: density interpolated to Vy nodes
    - RHOFYSUM: fluid density interpolated to Vy nodes
    - KYSUM: thermal conductivity interpolated to Vy nodes
    - PHIYSUM: porosity interpolated to Vy nodes
    - RYSUM: ηfluid/kϕ interpolated to Vy nodes
    - WTYSUM: weight for bilinear interpolation to Vy nodes

# Returns

    - nothing
"""
function marker_to_vy_nodes!(
    m,
    xmm,
    ymm,
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
    i, j, weights = fix_weights(
        xmm,
        ymm,
        xvy,
        yvy,
        dx,
        dy,
        jmin_vy,
        jmax_vy,
        imin_vy,
        imax_vy
    )
    @inbounds begin
        interpolate_add_to_grid!(i, j, weights, rhototalm[m], RHOYSUM)
        interpolate_add_to_grid!(i, j, weights, rhofluidcur[m], RHOFYSUM)
        interpolate_add_to_grid!(i, j, weights, ktotalm[m], KYSUM)
        interpolate_add_to_grid!(i, j, weights, phim[m], PHIYSUM)
        interpolate_add_to_grid!(i, j, weights, etafluidcur_inv_kphim[m], RYSUM)
        interpolate_add_to_grid!(i, j, weights, one(1.0), WTYSUM)
    end # @inbounds
    return nothing
end

"""
Interpolate selected marker properties to P nodes.

$(SIGNATURES)

# Details

    - m: marker number
    - xmm: marker x-position [m]
    - ymm: marker y-position [m]
	- inv_gggtotalm: inverse of total shear modulus of markers
	- sxxm: marker σ′xx [Pa]
	- rhototalm: total density of markers
	- rhocptotalm: total volumetric heat capacity of markers 
	- alphasolidcur: solid thermal expansion coefficient of markers 
	- alphafluidcur: fluid thermal expansion coefficient of markers
	- hrtotalm: total radiogenic heat production of markers
	- phim:  marker porosity
	- tkm_rhocptotalm: total thermal energy of markers
    - GGGPSUM: shear modulus interpolated to P nodes
    - SXXSUM: σ'xx interpolated to P nodes
    - RHOSUM: density interpolated to P nodes
    - RHOCPSUM: volumetric heat capacity interpolated to P nodes
    - ALPHASUM: thermal expansion interpolated to P nodes
    - ALPHAFSUM: fluid thermal expansion interpolated to P nodes
    - HRSUM: radioactive heating interpolated to P nodes
    - PHISUM: porosity interpolated to P nodes
    - TKSUM: heat capacity interpolated to P nodes
    - WTPSUM: weight for bilinear interpolation to P nodes

# Returns

    - nothing
"""
function marker_to_p_nodes!(
    m,
    xmm,
    ymm,
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
    i, j, weights = fix_weights(
        xmm,
        ymm,
        xp,
        yp,
        dx,
        dy,
        jmin_p,
        jmax_p,
        imin_p,
        imax_p
    )
    @inbounds begin
        interpolate_add_to_grid!(i, j, weights, inv_gggtotalm[m], GGGPSUM)
        interpolate_add_to_grid!(i, j, weights, sxxm[m], SXXSUM)
        interpolate_add_to_grid!(i, j, weights, rhototalm[m], RHOSUM)
        interpolate_add_to_grid!(i, j, weights, rhocptotalm[m], RHOCPSUM)
        interpolate_add_to_grid!(i, j, weights, alphasolidcur[m], ALPHASUM)
        interpolate_add_to_grid!(i, j, weights, alphafluidcur[m], ALPHAFSUM)
        interpolate_add_to_grid!(i, j, weights, hrtotalm[m], HRSUM)
        interpolate_add_to_grid!(i, j, weights, phim[m], PHISUM)
        interpolate_add_to_grid!(i, j, weights, tkm_rhocptotalm[m], TKSUM)
        interpolate_add_to_grid!(i, j, weights, one(1.0), WTPSUM)
    end # @inbounds
    return nothing
end

"""
Interpolate marker wet silicate (solid) fraction to P nodes.

$(SIGNATURES)

# Details

    - m: marker number
    - xmm: marker x-position [m]
    - ymm: marker y-position [m]
    - XWsolidm0: previous marker wet silicate (solid) molar fraction
    - XWSSUM: wet silicate (solid) molar fraction interpolated to P nodes
    - WTPSUM: weight for bilinear interpolation to P nodes

# Returns

    - nothing
"""
function molarfraction_marker_to_p_nodes!(
    m,
    xmm,
    ymm,
    XWsolidm0,
    XWSSUM,
    WTPSUM
)
    i, j, weights = fix_weights(
        xmm,
        ymm,
        xp,
        yp,
        dx,
        dy,
        jmin_p,
        jmax_p,
        imin_p,
        imax_p
    )
    @inbounds begin
        interpolate_add_to_grid!(i, j, weights, XWsolidm0[m], XWSSUM)
        interpolate_add_to_grid!(i, j, weights, one(1.0), WTPSUM)
    end # @inbounds
    return nothing
end

"""
Compute P nodes wet silicate (solid) fraction from markers.

$(SIGNATURES)

# Details

    - xm: x-coordinate of markers
    - ym: y-coordinate of markers
    - XWsolidm0: previous marker wet silicate (solid) molar fraction
    - XWS: wet silicate (solid) molar fraction at P nodes
    - XWSSUM: wet silicate (solid) molar fraction interpolated to P nodes
    - WTPSUM: weight for bilinear interpolation to P nodes
    - marknum: current total number of markers in use

# Returns

    - nothing
"""
function update_p_nodes_melt_composition!(
    xm, ym, XWsolidm0, XWS, XWSSUM, WTPSUM, marknum)
    XWSSUM .= zero(0.0)
    WTPSUM .= zero(0.0)
    for m=1:1:marknum
        molarfraction_marker_to_p_nodes!(
            m, xm[m], ym[m], XWsolidm0, XWSSUM, WTPSUM)
    end
    compute_molarfraction!(XWSSUM, WTPSUM, XWS)
    return nothing
end

"""
Compute properties of basic nodes based on interpolation arrays.

$(SIGNATURES)

# Details

    - ETA0SUM: ETA0 interpolation array
    - ETASUM: ETA interpolation array
    - GGGSUM: GGG interpolation array
    - SXYSUM: SXY interpolation array
    - COHSUM: COH interpolation array
    - TENSUM: TEN interpolation array
    - FRISUM: FRI interpolation array
    - WTSUM: WT interpolation array
    - ETA0: ETA0 basic node array
    - ETA: ETA basic node array
    - GGG: GGG basic node array
    - SXY0: SXY basic node array
    - COH: COH basic node array
    - TEN: TEN basic node array
    - FRI: FRI basic node array
    - YNY: YNY basic node array

# Returns

    - nothing

"""
function compute_basic_node_properties!(
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
@timeit to "compute_basic_node_properties!" begin
    @inbounds begin
        for j=1:1:Nx, i=1:1:Ny
            if WTSUM[i, j] > 0.0 
                ETA0[i, j] = ETA0SUM[i, j] * inv(WTSUM[i, j])
                ETA[i, j] = ETASUM[i, j] * inv(WTSUM[i, j])
                if ETA[i, j] < ETA0[i, j]
                    YNY[i, j] = true
                end
                GGG[i, j] = inv(GGGSUM[i, j]) * WTSUM[i, j]
                SXY0[i, j] = SXYSUM[i, j] * inv(WTSUM[i, j])
                COH[i, j] = COHSUM[i, j] * inv(WTSUM[i, j])
                TEN[i, j] = TENSUM[i, j] * inv(WTSUM[i, j])
                FRI[i, j] = FRISUM[i, j] * inv(WTSUM[i, j])
            end
        end 
    end # @inbounds
end # @timeit to "compute_basic_node_properties!"
    return nothing
end # function compute_basic_node_properties!

"""
Compute properties of Vx nodes based on interpolation arrays.

$(SIGNATURES)

# Details

    - RHOXSUM: RHOX interpolation array
    - RHOFXSUM: RHOFX interpolation array
    - KXSUM: KX interpolation array
    - PHIXSUM: PHIX interpolation array
    - RXSUM: RX interpolation array
    - WTXSUM: WTX interpolation array
    - RHOX: RHOX Vx node array
    - RHOFX: RHOFX Vx node array
    - KX: KX Vx node array
    - PHIX: PHIX Vx node array
    - RX: RX Vx node array

# Returns

    - nothing

"""
function compute_vx_node_properties!(
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
@timeit to "compute_vx_node_properties!" begin
    @inbounds begin
        for j=1:1:Nx1, i=1:1:Ny1
            if WTXSUM[i, j] > 0.0 
                RHOX[i, j] = RHOXSUM[i, j] * inv(WTXSUM[i, j])
                RHOFX[i, j] = RHOFXSUM[i, j] * inv(WTXSUM[i, j])
                KX[i, j] = KXSUM[i, j] * inv(WTXSUM[i, j])
                PHIX[i, j] = PHIXSUM[i, j] * inv(WTXSUM[i, j])
                RX[i, j] = RXSUM[i, j] * inv(WTXSUM[i, j])
            end
        end
    end # @inbounds
end # @timeit to "compute_vx_node_properties!"
    return nothing
end # function compute_vx_node_properties!

"""
Compute properties of Vy nodes based on interpolation arrays.

$(SIGNATURES)

# Details

    - RHOYSUM: RHOY interpolation array
    - RHOFYSUM: RHOFY interpolation array
    - KYSUM: KY interpolation array
    - PHIYSUM: PHIY interpolation array
    - RYSUM: RY interpolation array
    - WTYSUM: WTY interpolation array
    - RHOY: RHOY Vy node array
    - RHOFY: RHOFY Vy node array
    - KY: KY Vy node array
    - PHIY: PHIY Vy node array
    - RY: RY Vy node array

# Returns

    - nothing

"""
function compute_vy_node_properties!(
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
@timeit to "compute_vy_node_properties!" begin
    @inbounds begin
        for j=1:1:Nx1, i=1:1:Ny1
            if WTYSUM[i, j] > 0.0 
                RHOY[i, j] = RHOYSUM[i, j] * inv(WTYSUM[i, j])
                RHOFY[i, j] = RHOFYSUM[i, j] * inv(WTYSUM[i, j])
                KY[i, j] = KYSUM[i, j] * inv(WTYSUM[i, j])
                PHIY[i, j] = PHIYSUM[i, j] * inv(WTYSUM[i, j])
                RY[i, j] = RYSUM[i, j] * inv(WTYSUM[i, j])
            end
        end
    end # @inbounds
end # @timeit to "compute_vy_node_properties!"
    return nothing
end # function compute_vy_node_properties!

"""
Compute properties of P nodes based on interpolation arrays.

$(SIGNATURES)

# Details

    - GGGPSUM: GGGP interpolation array
    - SXX0SUM: SXX0 interpolation array
    - RHOSUM: RHO interpolation array
    - RHOCPSUM: RHOCP interpolation array
    - ALPHASUM: ALPHA interpolation array
    - ALPHAFSUM: ALPHAF interpolation array
    - HRSUM: HR interpolation array
    - PHISUM: PHI interpolation array
    - TKSUM: TK interpolation array
    - WTPSUM: WTP interpolation array
    - GGGP: GGGP P node array
    - SXX0: SXX0 P node array
    - RHO: RHO P node array
    - RHOCP: RHOCP P node array
    - ALPHA: ALPHA P node array
    - ALPHAF: ALPHAF P node array
    - HR: HR P node array
    - PHI: PHI P node array
    - BETAPHI: BETAPHI P node array
    - tk1: tk1 P node array

# Returns

    - nothing

"""
function compute_p_node_properties!(
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
@timeit to "compute_p_node_properties!" begin
    @inbounds begin
        for j=1:1:Nx1, i=1:1:Ny1
            if WTPSUM[i, j] > 0.0
                RHO[i, j] = RHOSUM[i, j] * inv(WTPSUM[i, j])
                RHOCP[i, j] = RHOCPSUM[i, j] * inv(WTPSUM[i, j])
                ALPHA[i, j] = ALPHASUM[i, j] * inv(WTPSUM[i, j])
                ALPHAF[i, j] = ALPHAFSUM[i, j] * inv(WTPSUM[i, j])
                HR[i, j] = HRSUM[i, j] * inv(WTPSUM[i, j])
                GGGP[i, j] = inv(GGGPSUM[i, j]) * WTPSUM[i, j]
                SXX0[i, j] = SXXSUM[i, j] * inv(WTPSUM[i, j])
                tk1[i, j] = TKSUM[i, j] * inv(RHOCPSUM[i, j])
                PHI[i, j] = PHISUM[i, j] * inv(WTPSUM[i, j])
                BETAPHI[i, j] = inv(GGGP[i, j]) * PHI[i, j]
            end
        end
    end # @inbounds
end # @timeit to "compute_p_node_properties!"
    return nothing
end # function compute_p_node_properties!

"""
Compute thermodynamic properties at P nodes based on interpolation arrays.

$(SIGNATURES)

# Details

    - DMPSUM: DMP interpolation array
    - DHPSUM: DHP interpolation array
    - WTPSUM: WTP interpolation array
    - DMP: mass transfer term at P nodes
    - DHP: enthalpy transfer/latent heating term at P nodes

# Returns

    - nothing
"""
function compute_thermodynamic_xfer!(DMPSUM, DHPSUM, WTPSUM, DMP, DHP)
@timeit to "compute_thermodynamic_xfer!" begin
    @inbounds begin
        for j=1:1:Nx1, i=1:1:Ny1
            if WTPSUM[i, j] > 0.0 
                DMP[i, j] = DMPSUM[i, j] * inv(WTPSUM[i, j])
                DHP[i, j] = DHPSUM[i, j] * inv(WTPSUM[i, j])
            else
                DMP[i, j] = DHP[i, j] = zero(0.0)
            end
        end 
    end # @inbounds
end # @timeit to "compute_thermodynamic_xfer!"
    return nothing
end # function compute_thermodynamic_xfer!

"""
Compute wet solid molar fraction at P nodes based on interpolation arrays.

$(SIGNATURES)

# Details

    - XWSSUM: XWX interpolation array
    - WTPSUM: WTP interpolation array
    - XWS: wet solid molar fraction at P nodes

# Returns

    - nothing
"""
function compute_molarfraction!(XWSSUM, WTPSUM, XWS)
@timeit to "compute_molarfraction!" begin
    @inbounds begin
        for j=1:1:Nx1, i=1:1:Ny1
            if WTPSUM[i, j] > 0.0 
                XWS[i, j] = XWSSUM[i, j] * inv(WTPSUM[i, j])
            else
                XWS[i, j] = zero(0.0)
            end
        end 
    end # @inbounds
end # @timeit to "compute_molarfraction!"
    return nothing
end # function compute_molarfraction!

"""
Compute volumetric isobaric heat capacity of H₂O (fluid phase)
based on temperature.

$(SIGNATURES)

# Details
    
        - T: temperature [K]
        - mode: marker property computation mode
            - 1: dynamic, based on (Touloukian, 1970; Hobbs, 1974;
                 Travis and Schubert, 2005)
            - 2: constant parameter rhocpfluidm

# Returns
    
        - ρᶠCₚᶠ: volumetric isobaric heat capacity of fluid
"""
function compute_rhocpfluidm(T, mode)
@timeit to "compute_rhocpfluidm" begin
    if mode == 1 
        if T < tmfluidphase-5.0
            ρᶠCₚᶠ = ρH₂Oᶠⁱ * 7.67T 
        elseif T < tmfluidphase
            ρᶠCₚᶠ = ρH₂Oᶠⁱ * (7.67T + 0.1Lᶠ)
        elseif T < tmfluidphase+5.0
            ρᶠCₚᶠ = ρH₂Oᶠ * (4200.0 + 0.1Lᶠ)
        elseif T < 410.0
            ρᶠCₚᶠ = ρH₂Oᶠ * 4200.0
        else
            ρᶠCₚᶠ = ρH₂Oᶠ * (-4.67e4 + 333T - 0.731T^2 + 5.4e-4T^3) 
        end
    elseif mode == 9
        @inbounds ρᶠCₚᶠ = rhocpfluidm[1]
    else
        throw("unknown mode $mode") 
    end
end # @timeit to "compute_rhocpfluidm"
    return ρᶠCₚᶠ
end # function compute_rhocpfluidm

"""
Compute thermal conductivity of silicate (solid phase) based on temperature.

$(SIGNATURES)

# Details
    
        - T: temperature [K]
        - mode:
            - 1: dynamic, based on (Gerya, 2019)
            - 9: constant parameter ksolidm

# Returns
    
        - kᶠ: thermal conductivity of solid
"""
function compute_ksolidm(T, mode)
@timeit to "compute_ksolidm" begin
    if mode == 1
        kˢ = 0.73 + 1293.0/(T+77.0)
    elseif mode == 9
        @inbounds kˢ = ksolidm[1]
    else
        throw("unknown mode $mode") 
    end
end # @timeit to "compute_ksolidm"
    return kˢ
end # function compute_ksolidm

"""
Compute thermal conductivity of H₂O (fluid phase) based on temperature.

$(SIGNATURES)

# Details
    
        - T: temperature [K]
        - mode:
            - 1: dynamic, based on (Touloukian, 1970; Hobbs, 1974;
                 Grimm & Mcsween, 1989; Bland & Travis, 2017)
            - 9: constant parameter ksolidm

# Returns
    
        - kᶠ: thermal conductivity of fluid
"""
function compute_kfluidm(T, mode)
@timeit to "compute_kfluidm" begin
    if mode == 1
        if T < tmfluidphase
            kᶠ = 0.465 + 488.0/T
        elseif T < 410
            kᶠ = -0.581 + 6.34e-3T - 7.93e-6T^2
        else
            kᶠ = -0.142 + 4.12e-3T - 5.01e-6T^2
        end
    elseif mode == 9
        @inbounds kᶠ = kfluidm[1]
    else
        throw("unknown mode $mode") 
    end
end # @timeit to "compute_kfluidm"
    return kᶠ
end # function compute_kfluidm

"""
Compute dehydration reaction time Δtreaction based on temperature and porosity 
according to selected method:

$(SIGNATURES)

# Details

    - T: temperature
    - ϕ: porosity
    - mode:
        - 1: Gaussian form reaction rate coefficient, based on (Martin & Fyfe,
             1970; Emmanuel & Berkowitz, 2006; Iyer et al., 2012). 
        - 2: pseudo-Arrhenius form reaction rate coefficient, based on
             (Bland & Travis, 2017).
        - 3: Arrhenius form reaction coefficent, based on (Travis et al., 2018).
        - 9: constant parameter Δtreaction
    
# Returns

    - Δtreaction: dehydration reaction time
"""
function compute_Δtreaction(T, ϕ, mode)
@timeit to "compute_Δtreaction" begin
    if mode == 1
        Δtr = -log_completion_rate / (A_I*ϕ) * exp(b_I*(T-c_I)^2)
    elseif mode == 2
        Δtr = -log_completion_rate / (Sxo_B*ϕ) * 2.0^((To_B-T)/Tscl_B)
    elseif mode == 3
        Δtr = -log_completion_rate / (Sxo_B*ϕ) * exp(Ea_T / RG * (1.0/T - 1.0/To_T))
    elseif mode == 9
        Δtr = Δtreaction
    else
        throw("unknown mode $mode")
    end
end # @timeit to "compute_Δtreaction"
    return Δtr
end # function compute_dtreaction

"""
Compute solid phase thermal conductivity for given temperature and method.

$(SIGNATURES)

# Details
- T: temperature
- method:
    - 1:
    - 2:

# Returns

- dtreaction
"""

"""
Apply insulating boundary conditions to given array:

[x x x x x x; x a b c d x; x e f g h x; x x x x x x]

becomes

[a a b c d d; a a b c d d; e e f g h h; e e f g h h]
     
# Details

    - t: array to apply insulating boundary conditions to

# Returns

    - nothing
"""
function apply_insulating_boundary_conditions!(t)
@timeit to "apply_insulating_boundary_conditions!" begin
    Nyy, Nxx = size(t)
    if Nyy>2 && Nxx>2
        @inbounds begin
            # upper boundary
            @views @. t[1, 2:Nxx-1] = t[2, 2:Nxx-1]
            # lower boundary
            @views @. t[Nyy, 2:Nxx-1] = t[Nyy-1, 2:Nxx-1]
            # left boundary
            @views @. t[:, 1] = t[:, 2]
            # right boundary
            @views @. t[:, Nxx] = t[:, Nxx-1]
        end # @inbounds
    end
end # @timeit to "apply_insulating_boundary_conditions!"
    return nothing
end

"""
Compute molar Gibbs free energy for single dehydration reaction
Wsilicate = Dsilicate + H₂O (16.144)

$(SIGNATURES)

# Details

    - T: temperature
    - pf: fluid pressure
    - XDˢ: molar fraction of dry solid
    - XWˢ: molar fraction of wet solid
    - Δt: timestep size
    - Δtr: total reaction time Δtreaction

# Returns

    - ΔGWD: molar Gibbs free energy for single dehydration reaction (16.165a/b).
"""
function compute_gibbs_free_energy(T, pf, XDˢ, XWˢ, Δt, Δtr)
@timeit to "compute_gibbs_free_energy" begin
    # compute incomplete reaction for short timestep Δt < Δtreaction
    if Δt < Δtr
        # compute ΔG for dehydration reaction (16.145), (16.165b)
        ΔGWD = (ΔHWD - T*ΔSWD + pf*ΔVWD + RG*T*log(XDˢ/XWˢ)) * (1.0 - Δt/Δtr)
    else
        # Δt ≥ Δtreaction (16.165a)
        ΔGWD = zero(0.0)    
    end
end # @timeit to "compute_gibbs_free_energy"
    return ΔGWD
end # function compute_gibbs_free_energy

"""
Compute relative enthalpy of system for single dehydration reaction
Wsilicate = Dsilicate + H₂O (16.144).

$(SIGNATURES)

# Details

    - T: temperature
    - pf: fluid pressure
    - XDsolid: molar fraction of dry solid
    - XWsolid: molar fraction of wet solid
    - Δt: timestep size

# Returns

    - Hᵗ: relative enthalpy of system for single dehydration reaction (16.163)
"""
function compute_relative_enthalpy(Xsolid, XWsolid)
    @timeit to "compute_relative_enthalpy" return -Xsolid * XWsolid * ΔHWD / (MD+MH₂O)
end # function compute_relative_enthalpy

"""
Compute dehydration reaction constant (16.151).

$(SIGNATURES)

# Details

    - T: temperature
    - pf: fluid pressure
    - ΔGWD: Gibbs free energy for dehydration reaction

# Returns

    - KWD: dehydration reaction constant (16.151)
"""
function compute_reaction_constant(T, pf, ΔGWD)
    # compute reaction constant (16.151)
    @timeit to "compute_reaction_constant" return exp(-(ΔHWD - T*ΔSWD + ΔVWD*pf - ΔGWD) / (RG*T))
end # function compute_reaction_constant


"""
Perform hydrothermomechanical iterations to time step thermal field at P nodes.

$(SIGNATURES)

# Details

    - DMP: mass transfer term at P nodes
    - DHP: enthalpy transfer/latent heating term at P nodes
    - DMPSUM: interpolation of DMP (mass transfer term) at P nodes
    - DHPSUM: interpolation of DHP (enthalpy transfer term) at P nodes
    - WTPSUM: interpolation weights at P nodes 
    - pf: fluid pressure at P nodes
    - tk2: next temperature at P nodes 
    - tm: type of markers
    - xm: x-coordinate of markers
    - ym: y-coordinate of markers
    - XWˢm₀: previous marker wet silicate (solid) fraction
    - XWˢm: current marker wet silicate (solid) fraction
    - phim: current marker porosity
    - phinewm: next generation marker porosity
    - pfm₀: previous marker fluid pressure
    - marknum: current total number of markers
    - Δt: current time step length
    - timestep: current time step
    - titer: current thermochemical iteration number

# Returns

    - nothing
"""
function perform_thermochemical_reaction!(
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
    XWˢm₀,
    XWˢm,
    phim,
    phinewm,
    pfm₀,
    marknum,
    Δt,
    timestep,
    titer
)
@timeit to "perform_thermochemical_reaction!" begin
    # reset interpolation arrays
    reset_thermochemical_properties!(DMPSUM, DHPSUM, WTPSUM)
    # iterate over markers
    @inbounds begin
        for m=1:1:marknum
            # for 
            if tm[m] < 3
                # for rocks only
                i, j, weights = fix_weights(
                    xm[m],
                    ym[m],
                    xp,
                    yp,
                    dx,
                    dy,
                    jmin_p,
                    jmax_p,
                    imin_p,
                    imax_p
                )
                # interpolate temperature from P nodes
                tknm = dot4(grid_vector(i, j, tk2), weights)
                # interpolate non-negative fluid pressure from P nodes
                pfnm = max(zero(0.0), dot4(grid_vector(i, j, pf), weights))
                # factor in previous iteration marker fluid pressure
                if titer > 2
                    pfnm = pfnm*(1.0-pfcoeff) + pfm₀[m]*pfcoeff
                end
                # store current marker fluid pressure for next iteration
                pfm₀[m] = pfnm
                # compute bulk composition of solid and fluid system:
                # compute previous dry solid molar fraction (16.146)
                XDˢm₀ = 1.0 - XWˢm₀[m]
                # get fluid molar volume
                VH₂O = ifelse(tknm>tmfluidphase, VH₂Oᶠ, VH₂Oᶠⁱ)
                # compute previous fluid molar fraction (16.164)
                Xᶠ₀ = phim[m]*(XWˢm₀[m]*VWˢ + XDˢm₀*VDˢ) / (
                    (1.0-phim[m])*VH₂O +
                    phim[m] * (XWˢm₀[m]*VWˢ + XDˢm₀*VDˢ)
                )
                # compute previous equilibrium solid molar fraction (16.150)
                Xˢ₀ = 1.0 - Xᶠ₀
                # compute previous water molar fraction (16.147)
                XH₂Oᵗ = (XWˢm₀[m]*Xˢ₀ + Xᶠ₀) / (1.0 + XWˢm₀[m]*Xˢ₀)
                # compute dry solid molar fraction (16.149)
                XDᵗ = 1.0 - XH₂Oᵗ
                # compute previous solid density (16.161)
                ρˢ₀ = (MD + MH₂O*XWˢm₀[m]) / (VDˢ*XDˢm₀ + VWˢ*XWˢm₀[m])
                # compute previous fluid density (16.162)
                ρᶠ₀ = ifelse(tknm>tmfluidphase, ρH₂Oᶠ, ρH₂Oᶠⁱ)
                # compute Δtreaction
                Δtr = compute_Δtreaction(
                    tknm, phim[m], reaction_rate_coeff_mode)
                # compute previous relative enthalpy of the system (16.163)
                Hᵗ₀ = compute_relative_enthalpy(Xˢ₀, XWˢm₀[m])
                # compute previous ΔG for dehydration reaction (16.165a/b)
                ΔGWD₀ = compute_gibbs_free_energy(
                    tknm, pfnm, XDˢm₀, XWˢm₀[m], Δt, Δtr) 
                # compute dehydration reaction constant (16.151)
                KWD = compute_reaction_constant(tknm, pfnm, ΔGWD₀)
                # compute reacted wet solid molar fraction (16.152)
                XWˢm₁ = inv(KWD + 1.0)
                # compute reacted dry solid molar fraction (16.153)
                XDˢm₁ = 1.0 - XWˢm₁
                # compute reacted total solid molar fraction (16.154)
                Xˢ₁ = XDᵗ / (1.0 - XDᵗ*XWˢm₁)
                # compute reacted fluid molar fraction (16.155)
                Xᶠ₁ = 1.0 - Xˢ₁
                # only process fluid-bearing rocks
                if 0.0 < Xᶠ₁ < 1.0
                    # compute reacted equilibrium porosity (16.156)
                    ϕ₁ = Xᶠ₁*VH₂O / (Xᶠ₁*VH₂O + Xˢ₁*(XWˢm₁*VWˢ+XDˢm₁*VDˢ))
                    # compute equilibrium solid density (16.161)
                    ρˢ₁ = (MD + MH₂O*XWˢm₁) / (VDˢ*XDˢm₁ + VWˢ*XWˢm₁)
                    # compute equilibrium fluid density (16.162)
                    ρᶠ₁ = ifelse(tknm>tmfluidphase, ρH₂Oᶠ, ρH₂Oᶠⁱ)
                    # compute equilibrium relative enthalpy of the system
                    # (16.163)
                    Hᵗ₁ = compute_relative_enthalpy(Xˢ₁, XWˢm₁)
                    # compute enthalpy change
                    ΔHᵗ = Hᵗ₁ - Hᵗ₀
                    # compute previous-to-reacted-equilibrium volume ratio
                    # (16.106)
                    RV = (ρˢ₁*(1.0-ϕ₁) + ρᶠ₁*ϕ₁) / (
                        ρˢ₀*(1.0-phim[m]) + ρᶠ₀*phim[m])
                    # compute mass transfer rate (16.103)
                    Γmass = (ρˢ₀*RV*(1.0-phim[m]) - ρˢ₁*(1.0-ϕ₁)) / Δt
                    # compute mass transfer term (16.112e)
                    ΔMm = (1.0-RV) / Δt
                    # compute enthalpy transfer/latent heating term (16.113)
                    ΔHm = Γmass * ΔHᵗ
                    # update wet solid (melt) molar fraction
                    XWˢm[m]=XWˢm₁
                    # update porosity
                    phinewm[m] = ϕ₁
                    # backload properties during first timestep
                    if timestep==1
                        XWˢm₀[m] = XWˢm[m]
                        phim[m] = phinewm[m]
                    end
                    # interpolate mass, enthalpy transfer terms to P nodes
                    interpolate_add_to_grid!(i, j, weights, ΔMm, DMPSUM)
                    interpolate_add_to_grid!(i, j, weights, ΔHm, DHPSUM)
                    interpolate_add_to_grid!(i, j, weights, one(1.0), WTPSUM)
                end
            end # if tm[m] < 3
        end # for m=1:1:marknum
        # compute thermodynamic properties at P nodes
        compute_thermodynamic_xfer!(DMPSUM, DHPSUM, WTPSUM, DMP, DHP)
    end # @inbounds
    @info "min/max mass transfer term" extrema(DMP)
    @info "min/max enthalpy transfer term" extrema(DHP)
end # @timeit to "perform_thermochemical_reaction!"
    return nothing
end # function perform_thermochemical_reaction!

"""
Compute gravity solution in P nodes to obtain
gravitational accelerations gx for Vx nodes, gy for Vy nodes.

$(SIGNATURES)

# Details

    - SP: solution vector
    - RP: right hand side vector
    - RHO: density at P nodes
    - FI: gravity potential at P nodes
    - gx: x gravitational acceleration at Vx nodes
    - gy: y gravitational acceleration at Vy nodes

# Returns

- nothing
"""
function compute_gravity_solution!(SP, RP, RHO, FI, gx, gy)
@timeit to "compute_gravity_solution!" begin
    # fresh LHS sparse coefficient matrix
    @timeit to "set up sparse matrix" LP = ExtendableSparseMatrix(Nx1*Ny1, Nx1*Ny1)
    # reset RHS coefficient vector
    RP .= 0.0
    # iterate over P nodes
    @timeit to "build system" begin
    for j=1:1:Nx1, i=1:1:Ny1
        # define global index in algebraic space
        gk = (j-1) * Ny1 + i
        # decide if external / boundary points
        @inbounds if (
            i==1 ||
            i==Ny1 ||
            j==1 ||
            j==Nx1 ||
            distance(xp[j], yp[i], xcenter, ycenter) > xcenter
        )
            # boundary condition: ϕ = 0
            updateindex!(LP, +, 1.0, gk, gk)
            # RP[gk] = 0.0 # already done at initialization
        else
            # internal points: 2D Poisson equation: gravitational potential Φ
            # ∂²Φ/∂x² + ∂²Φ/∂y² = 4KπGρ with K=2/3 for spherical 2D (11.10)
            #
            #           Φ₂
            #           |
            #           |
            #    Φ₁-----Φ₃-----Φ₅
            #           |
            #           |
            #           Φ₄
            #
            # density gradients
            # dRHOdx = (RHO[i, j+1]-RHO[i, j-1]) / 2 / dx
            # dRHOdy = (RHO[i+1, j]-RHO[i-1, j]) / 2 / dy
            # fill system of equations: LHS (11.11)
            updateindex!(LP, +, inv(dx^2), gk, gk-Ny1) # Φ₁
            updateindex!(LP, +, inv(dy^2), gk, gk-1) # Φ₂
            updateindex!(LP, +, -2.0*(inv(dx^2)+inv(dy^2)), gk, gk) # Φ₃
            updateindex!(LP, +, inv(dy^2), gk, gk+1) # Φ₄
            updateindex!(LP, +, inv(dx^2), gk, gk+Ny1) # Φ₅
            # fill system of equations: RHS (11.11)
            @inbounds RP[gk] = 4.0 * 2.0 * inv(3.0) * π * G * RHO[i, j]
        end
    end
    end # @timeit to "build system"
    @timeit to "solve system" begin
    # solve system of equations
    SP .= LP \ RP # implicit: flush!(LP)
    end # @timeit to "solve system"
    # reshape solution vector to 2D array
    @timeit to "reshape solution" begin
    FI .= reshape(SP, Ny1, Nx1)
    end # @timeit to "reshape solution"
    @timeit to "compute accelerations" begin
    # gx = -∂ϕ/∂x (11.12)
    @inbounds gx[:, 1:Nx] .= -diff(FI, dims=2) ./ dx
    # gy = -∂ϕ/∂y (11.13)   
    @inbounds gy[1:Ny, :] .= -diff(FI, dims=1) ./ dy
    end # @timeit to "compute accelerations"
end # @timeit to "compute_gravity_solution!"
    return nothing
end # function compute_gravity_solution!


"""
Assemble the LHS sparse coefficient matrix and fill RHS coefficient vector
of the Poisson equation to be solved for the gravitational potential Φ.

$(SIGNATURES)

# Details

    - RHO: density at P nodes
    - RP: right hand side coefficient vector

# Returns

    - LP: LHS sparse coefficient matrix

"""
function assemble_gravitational_lse!(RHO, RP)
    @timeit to "assemble_gravitational_lse!" begin
        # fresh LHS sparse coefficient matrix
        @timeit to "setup sparse LHS" LP = ExtendableSparseMatrix(Nx1*Ny1, Nx1*Ny1)
        # reset RHS coefficient vector
        @timeit to "reset RHS" RP .= zero(0.0)
        # iterate over P nodes
        @timeit to "build system" begin
        for j=1:1:Nx1, i=1:1:Ny1
            # define global index in algebraic space
            gk = (j-1) * Ny1 + i
            # decide if external / boundary points
            @inbounds if (
                i==1 ||
                i==Ny1 ||
                j==1 ||
                j==Nx1 ||
                distance(xp[j], yp[i], xcenter, ycenter) > xcenter
            )
                # boundary condition: ϕ = 0
                updateindex!(LP, +, 1.0, gk, gk)
                # RP[gk] = 0.0 # already done at initialization
            else
                # internal points: 2D Poisson equation: gravitational potential Φ
                # ∂²Φ/∂x² + ∂²Φ/∂y² = 4KπGρ with K=2/3 for spherical 2D (11.10)
                #
                #           Φ₂
                #           |
                #           |
                #    Φ₁-----Φ₃-----Φ₅
                #           |
                #           |
                #           Φ₄
                #
                # density gradients
                # dRHOdx = (RHO[i, j+1]-RHO[i, j-1]) / 2 / dx
                # dRHOdy = (RHO[i+1, j]-RHO[i-1, j]) / 2 / dy
                # fill system of equations: LHS (11.11)
                updateindex!(LP, +, inv(dx^2), gk, gk-Ny1) # Φ₁
                updateindex!(LP, +, inv(dy^2), gk, gk-1) # Φ₂
                updateindex!(LP, +, -2.0*(inv(dx^2)+inv(dy^2)), gk, gk) # Φ₃
                updateindex!(LP, +, inv(dy^2), gk, gk+1) # Φ₄
                updateindex!(LP, +, inv(dx^2), gk, gk+Ny1) # Φ₅
                # fill system of equations: RHS (11.11)
                @inbounds RP[gk] = 4.0 * 2.0 * inv(3.0) * π * G * RHO[i, j]
            end
        end
        end # @timeit to "build system"
    end # @timeit to "assemble_gravitational_lse!"
    return LP
end

"""
Process gravitational potential solution vector to output physical observables.

$(SIGNATURES)

# Details

    - FI: gravitational potential
    - gx: x-component of gravitational acceleration
    - gy: y-component of gravitational acceleration

# Returns

    - nothing
"""
function process_gravitational_solution!(SP, FI, gx, gy)
@timeit to "process gravitational solution" begin
    @timeit to "reshape solution" begin
    FI .= reshape(SP, Ny1, Nx1)
    end # @timeit to "reshape solution"
    @timeit to "compute accelerations" begin
    # gx = -∂ϕ/∂x (11.12)
    @inbounds gx[:, 1:Nx] .= -diff(FI, dims=2) ./ dx
    # gy = -∂ϕ/∂y (11.13)   
    @inbounds gy[1:Ny, :] .= -diff(FI, dims=1) ./ dy
    end # @timeit to "compute accelerations"
end # @timeit to "process gravitational solution"
    return nothing
end

"""
Compute viscosities, stresses, and density gradients
for hydromechanical solver.

$(SIGNATURES)

# Details

## In

    - ETA: viscoplastic viscosity at basic nodes
    - ETAP: viscosity at P nodes
    - GGG: shear modulus at basic nodes
    - GGGP: shear modulus at P nodes
    - SXY0: σ₀xy XY stress at basic nodes
    - SXX0:σ₀xy XY stress at basic nodes
    - RHOX: density at Vx nodes
    - RHOY: density at Vy nodes
    - dt: time step

## Out 

    - ETAcomp: computational viscosity at basic nodes
    - ETAPcomp: computational viscosity at P nodes
    - SXYcomp: previous XY stresses at basic nodes
    - SXXcomp: previous XX stresses at P nodes
    - SYYcomp: previous YY stresses at P nodes
    - dRHOXdx: density gradient at Vx nodes in x direction
    - dRHOXdy: density gradient at Vx nodes in y direction
    - dRHOYdx: density gradient at Vy nodes in x direction
    - dRHOYdy: density gradient at Vy nodes in y direction

# Returns

    - nothing
"""
function get_viscosities_stresses_density_gradients!(
    ETA,
    ETAP,
    GGG,
    GGGP,
    SXY0,
    SXX0,
    RHOX,
    RHOY,
    dt,
    ETAcomp,
    ETAPcomp,
    SXYcomp,
    SXXcomp,
    SYYcomp,
    dRHOXdx,
    dRHOXdy,
    dRHOYdx,
    dRHOYdy
)
@timeit to "get_viscosities_stresses_density_gradients!()" begin
    # computational viscosity
    @views @. ETAcomp = ETA*GGG*dt / (GGG*dt + ETA)
    @views @. ETAPcomp = ETAP*GGGP*dt / (GGGP*dt + ETAP)
    # previous stresses
    @views @. SXYcomp = SXY0*ETA / (GGG*dt+ETA)
    @views @. SXXcomp = SXX0*ETAP / (GGGP*dt+ETAP)
    @views @. SYYcomp = -SXX0*ETAP / (GGGP*dt+ETAP)
    # for erroneously undersized (Ny, Nx) SSX0, SSX
    # @views @. SXXcomp = (
        # SXX0*ETAP[1:Ny, 1:Nx] / (GGGP[1:Ny, 1:Nx]*dt + ETAP[1:Ny, 1:Nx])
    # )
    # @views @. SYYcomp = (
        # -SXX0*ETAP[1:Ny, 1:Nx] / (GGGP[1:Ny, 1:Nx]*dt+ETAP[1:Ny, 1:Nx])
    # )
    # density gradients
    @inbounds begin
    @views @. dRHOXdx[:, 2:Nx] = 0.5*(RHOX[:, 3:Nx1]-RHOX[:, 1:Nx1-2]) * inv(dx)
    @views @. dRHOXdy[2:Ny, :] = 0.5*(RHOX[3:Ny1, :]-RHOX[1:Ny1-2, :]) * inv(dy)
    @views @. dRHOYdx[:, 2:Nx] = 0.5*(RHOY[:, 3:Nx1]-RHOY[:, 1:Nx1-2]) * inv(dx)
    @views @. dRHOYdy[2:Ny, :] = 0.5*(RHOY[3:Ny1, :]-RHOY[1:Ny1-2, :]) * inv(dy)
    end # @inbounds
    return nothing
end # @timeit to "get_viscosities_stresses_density_gradients!()"
end # function get_viscosities_stresses_density_gradients!

"""
Set up hydromechanical linear system of equations structures.

$(SIGNATURES)

# Details

    - nothing

# Returns 

    - R: hydromechanical linear system of equations: RHS vector
    - S: hydromechanical linear system of equations: solution vector
"""
function setup_hydromechanical_lse()
@timeit to "setup_hydromechanical_lse()" begin
    # R = zeros(Ny1*Nx1*6)
    R = Vector{Float64}(undef, Ny1*Nx1*6)
    S = Vector{Float64}(undef, Ny1*Nx1*6)
end # @timeit to "setup_hydromechanical_lse()"
    return R, S
end

"""
Set up thermal linear system of equations structures.

$(SIGNATURES)

# Details

    - nothing

# Returns 

    - RT: thermal linear system of equations: RHS vector
    - ST: thermal linear system of equations: solution vector
"""
function setup_thermal_lse()
@timeit to "setup_thermal_lse()" begin
    # RT = zeros(Ny1*Nx1)
    RT = Vector{Float64}(undef, Ny1*Nx1)
    ST = Vector{Float64}(undef, Ny1*Nx1)
end # @timeit to "setup_thermal_lse()"
    return RT, ST
end

"""
Set up gravitational linear system of equations structures.

$(SIGNATURES)

# Details

    - nothing

# Returns 

    - RP: gravitational linear system of equations: RHS vector
    - SP: gravitational linear system of equations: solution vector
"""
function setup_gravitational_lse()
@timeit to "setup_gravitational_lse()" begin
    # RP = zeros(Ny1*Nx1)
    RP = Vector{Float64}(undef, Ny1*Nx1)
    SP = Vector{Float64}(undef, Ny1*Nx1)
end # @timeit to "setup_gravitational_lse()"
    return RP, SP
end

"""
Initialize `iparm` parameters of Pardiso MKL solver.

$(SIGNATURES)

# Details

    - ps: Instance of pardiso solver
    - iparms: dictionary of iparm parameters

# Returns

    - nothing
"""
function initialize_pardiso!(pardiso_solver, iparms)
    set_msglvl!(pardiso_solver, Pardiso.MESSAGE_LEVEL_OFF)
    set_matrixtype!(pardiso_solver, Pardiso.REAL_NONSYM)
    set_nprocs!(pardiso_solver, parse(Int32, ENV["OMP_NUM_THREADS"]))
    for (i, v) in iparms
        set_iparm!(pardiso_solver, i+1, v)
    end
end

"""
Assemble hydromechanical system of equations.

$(SIGNATURES)

# Details

    - ETA: viscosity at basic nodes
    - ETAP: viscosity at P nodes
    - GGG: shear modulus at basic nodes
    - GGGP: shear modulus at P nodes
    - SXY0: previous XY stress at basic nodes
    - SXX0: previous XX stress at P nodes
    - RHOX: total density at Vx nodes
    - RHOY: total density at Vy nodes
    - RHOFX: fluid density at Vx nodes
    - RHOFY: fluid density at Vy nodes
    - RX: ηfluid/Kϕ at Vx nodes
    - RY: ηfluid/Kϕ at Vy nodes
    - ETAPHI: bulk viscosity at P nodes
    - BETAPHI: bulk compressibility at P nodes
    - PHI: porosity at P nodes
    - gx: x gravitational acceleration at Vx nodes
    - gy: y gravitational acceleration at Vy nodes
    - pr0: previous total pressure at P nodes
    - pf0: previous fluid pressure at P nodes
    - DMP: mass transfer term at P nodes
    - dt: time step
    - R: vector to store RHS coefficients

# Returns

    - L: LHS coefficient matrix
"""
function assemble_hydromechanical_lse!(
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
    R
)
@timeit to "assemble_hydromechanical_lse()" begin
    # initialize LHS sparse coefficient matrix
    L = ExtendableSparseMatrix(Nx1*Ny1*6, Nx1*Ny1*6)
    # reset RHS coefficient vector
    R .= 0.0
    @inbounds begin
    for j=1:1:Nx1, i=1:1:Ny1
        # define global indices in algebraic space
        kvx = ((j-1)*Ny1 + i-1) * 6 + 1 # Vx solid
        kvy = kvx + 1 # Vy solid
        kpm = kvx + 2 # P total
        kqx = kvx + 3 # qx Darcy
        kqy = kvx + 4 # qy Darcy
        kpf = kvx + 5 # P fluid
        # Vx equation
        if i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1
            # Vx equation external points: boundary conditions
            # all locations: ghost unknowns Vx₃=0 -> 1.0⋅Vx[i,j]=0.0
            updateindex!(L, +, 1.0, kvx, kvx)
            # R[kvx] = 0.0 # already done with initialization
            # left boundary
            if j == 1 
                R[kvx] = vxleft
            end
            # right boundary
            if j == Nx 
                R[kvx] = vxright
            end
            # top boundary
            if i==1 && 1<j<Nx
                updateindex!(L, +, bctop, kvx, kvx+6)
            end
            # bottom boundary
            if i==Ny1 && 1<j<Nx
                updateindex!(L, +, bcbottom, kvx, kvx-6)
            end
        else
            # Vx equation internal points: x-Stokes
            #
            #                           kvx-6
            #                            Vx₂
            #                             |
            #               kvy-6     ETA(i-1,j)   kvy+6⋅Ny1-6
            #                Vy₁      GGG(i-1,j)     Vy₃
            #                 *       SXY0(i-1,j)     *
            #                           basic₁
            #                            ETA₁                       
            #                            SXY₁
            #               ETAP(i,j)     |      ETAP(i,j+1)
            #               GGGP(i,j)     |      GGGP(i,j+1) 
            #   kvx-6⋅Ny1   SXX0(i,j)    kvx     SXX0(i,j+1)  kvx+6⋅Ny1
            #     Vx₁---------P₁---------Vx₃---------P₂---------Vx₅
            #                kpm          |        kpm+6⋅Ny1
            #               ETAP₁         |        ETAP₂
            #               SXX₁          |        SXX₂
            #                          ETA(i,j) 
            #                kvy       GGG(i,j)     kvy+6⋅Ny1
            #                Vy₂       SXY0(i,j)      Vy₄
            #                 *        basic₂          * 
            #                            ETA₂ 
            #                            SXY₂
            #                             |
            #                           kvx+6
            #                            Vx₄
            #                             *
            #
            # computational viscosity
            ETA₁ = ETA[i-1, j] * GGG[i-1, j]*dt / (GGG[i-1, j]*dt + ETA[i-1, j])
            ETA₂ = ETA[i, j] * GGG[i, j]*dt / (GGG[i, j]*dt + ETA[i, j])
            ETAP₁ = ETAP[i, j] * GGGP[i, j]*dt / (GGGP[i, j]*dt + ETAP[i, j])
            ETAP₂ = ETAP[i, j+1] * GGGP[i, j+1]*dt / (
                GGGP[i, j+1]*dt + ETAP[i, j+1])
            # previous stresses
            SXY₁ = SXY0[i-1, j]* ETA[i-1, j] / (GGG[i-1, j]*dt + ETA[i-1, j])
            SXY₂ = SXY0[i, j] * ETA[i, j] / (GGG[i, j]*dt + ETA[i, j])
            SXX₁ = SXX0[i, j] * ETAP[i, j] / (GGGP[i, j]*dt + ETAP[i, j])
            SXX₂ = SXX0[i, j+1] * ETAP[i, j+1] / (
                GGGP[i, j+1]*dt + ETAP[i, j+1])
            # density gradients
            ∂RHO∂x = 0.5 * (RHOX[i, j+1] - RHOX[i, j-1]) * inv(dx)
            ∂RHO∂y = 0.5 * (RHOX[i+1, j] - RHOX[i-1, j]) * inv(dy)
            # LHS coefficient matrix
            updateindex!(L, +, ETAP₁/dx^2, kvx, kvx-6*Ny1) # Vx₁
            updateindex!(L, +, ETA₁/dy^2, kvx, kvx-6) # Vx₂
            updateindex!(
                L,
                +,
                (
                    -(ETAP₁+ETAP₂) * inv(dx^2)
                    -(ETA₁+ETA₂) * inv(dy^2)
                    -∂RHO∂x * gx[i, j] * dt
                ),
                kvx,
                kvx
            ) # Vx₃
            updateindex!(L, +, ETA₂/dy^2, kvx, kvx+6) # Vx₄
            updateindex!(L, +, ETAP₂/dx^2, kvx, kvx+6*Ny1) # Vx₅
            updateindex!(
                L,
                +,
                (
                    ETAP₁ * inv(dx) * inv(dy)
                    -ETA₂ * inv(dx) * inv(dy)
                    -∂RHO∂y * gx[i, j] * dt * 0.25
                ),
                kvx,
                kvy
            ) # Vy₂
            updateindex!(
                L,
                +,
                (
                    -ETAP₂ * inv(dx) * inv(dy)
                    +ETA₂ * inv(dx) * inv(dy)
                    -∂RHO∂y * gx[i, j] * dt * 0.25
                ),
                kvx,
                kvy+6*Ny1
            ) # Vy₄
            updateindex!(
                L,
                +,
                (
                    -ETAP₁ * inv(dx) * inv(dy)
                    +ETA₁ * inv(dx) * inv(dy)
                    -∂RHO∂y * gx[i, j] * dt * 0.25
                ),
                kvx,
                kvy-6
            ) # Vy₁
            updateindex!(
                L,
                +,
                (
                    ETAP₂ * inv(dx) * inv(dy)
                    -ETA₁ * inv(dx) * inv(dy)
                    -∂RHO∂y * gx[i, j] * dt * 0.25
                ),
                kvx,
                kvy+6*Ny1-6
            ) # Vy₃
            updateindex!(L, +, Kcont*inv(dx), kvx, kpm) # P₁
            updateindex!(L, +, -Kcont*inv(dx), kvx, kpm+6*Ny1) # P₂
            # RHS coefficient vector
            R[kvx] = (
                -RHOX[i, j] * gx[i, j]
                -(SXY₂-SXY₁) * inv(dy)
                -(SXX₂-SXX₁) * inv(dx)
            )
        end # Vx equation
        # Vy equation
        if i==1 || i==Ny || i==Ny1 || j==1 || j==Nx1
            # Vy equation external points: boundary conditions
            # all locations: ghost unknowns Vy₃=0 -> 1.0⋅Vy[i,j]=0.0
            updateindex!(L, +, 1.0, kvy, kvy)
            # R[kvy] = 0.0 # already done with initialization
            # top boundary
            if i == 1 
                R[kvy] = vytop
            end
            # bottom boundary
            if i == Ny 
                R[kvy] = vybottom
            end
            # left boundary
            if j==1 && 1<i<Ny
                updateindex!(L, +, bcleft, kvy, kvy+6*Ny1)
            end
            # right boundary
            if j==Nx1 && 1<i<Ny
                updateindex!(L, +, bcright, kvy, kvy-6*Ny1)
            end
        else
            # Vy equation internal points: y-Stokes
            #
            #                           kvy-6
            #                            Vy₂
            #                             |
            #                          ETAP(i,j)
            #                          GGGP(i,j)
            #             kvx-6⋅Ny1    SXX0(i,j)     kvx
            #                Vx₁          P₁         Vx₃
            #                 *         ETAP₁         *
            #                            SYY₁
            #               ETA(i,j-1)   kpm       ETA(i,j)
            #               GGG(i,j-1)    |        GGG(i,j)
            #   kvy-6⋅Ny1   SXY0(i,j-1)  kvy       SXY0(i,j)  kvy+6⋅Ny1
            #     Vy₁-------basic₁-------Vv₃-------basic₂-------Vy₅
            #               ETA₁          |        ETA₂     
            #               SXY₁          |        SXY₂
            #                            kpm+6
            #                         ETAP(i+1,j)
            #                         GGGP(i+1,j)
            #          kvx-6⋅Ny1+6    SXX0(i+1,j)   kvx+6
            #                Vx₂          P₂        Vx₄
            #                 *         ETAP₂        *  
            #                            SYY₂
            #                             |
            #                           kvy+6
            #                            Vy₄
            #
            # computational viscosity
            ETA₁ = ETA[i, j-1] * GGG[i, j-1]*dt / (GGG[i, j-1]*dt + ETA[i, j-1])
            ETA₂ = ETA[i, j] * GGG[i, j]*dt / (GGG[i, j]*dt + ETA[i, j])
            ETAP₁ = ETAP[i, j] * GGGP[i, j]*dt / (GGGP[i, j]*dt + ETAP[i, j])
            ETAP₂ = ETAP[i+1, j] * GGGP[i+1, j]*dt / (
                GGGP[i+1, j]*dt + ETAP[i+1, j])
            # previous stresses
            SXY₁ = SXY0[i, j-1] * ETA[i, j-1] / (GGG[i, j-1]*dt + ETA[i, j-1])
            SXY₂ = SXY0[i, j] * ETA[i, j] / (GGG[i, j]*dt + ETA[i, j])
            SYY₁ = -SXX0[i, j] * ETAP[i, j] / (GGGP[i, j]*dt + ETAP[i, j])
            SYY₂ = -SXX0[i+1, j] * ETAP[i+1, j] / (
                GGGP[i+1, j]*dt + ETAP[i+1, j])
            # density gradients
            ∂RHO∂x = 0.5 * (RHOY[i, j+1]-RHOY[i, j-1]) / dx
            ∂RHO∂y = 0.5 * (RHOY[i+1, j]-RHOY[i-1, j]) / dy
            # LHS coefficient matrix
            updateindex!(L, +, ETA₁/dx^2, kvy, kvy-6*Ny1) # Vy₁
            updateindex!(L, +, ETAP₁/dy^2, kvy, kvy-6) # Vy₂
            updateindex!(
                L,
                +,
                (
                    -(ETAP₁+ETAP₂) * inv(dy^2)
                    -(ETA₁+ETA₂) * inv(dx^2)
                    -∂RHO∂y * gy[i, j] * dt
                ),
                kvy,
                kvy
            ) # Vy₃
            updateindex!(L, +, ETAP₂ * inv(dy^2), kvy, kvy+6) # Vy₄
            updateindex!(L, +, ETA₂ * inv(dx^2), kvy, kvy+6*Ny1) # Vy₅
            updateindex!(
                L,
                +,
                (
                    ETAP₁ * inv(dx) * inv(dy)
                    -ETA₂ * inv(dx) * inv(dy)
                    -∂RHO∂x * gy[i, j] * dt * 0.25
                ),
                kvy,
                kvx
            ) # Vx₃
            updateindex!(
                L,
                +,
                (
                    -ETAP₂ * inv(dx) * inv(dy)
                    +ETA₂ * inv(dx) * inv(dy)
                    -∂RHO∂x * gy[i, j] * dt * 0.25
                ),
                kvy,
                kvx+6
            ) # Vx₄
            updateindex!(
                L,
                +,
                (
                    -ETAP₁ * inv(dx) * inv(dy)
                    +ETA₁ * inv(dx) * inv(dy)
                    -∂RHO∂x * gy[i, j] * dt * 0.25
                ),
                kvy,
                kvx-6*Ny1
            ) # Vx₁
            updateindex!(
                L,
                +,
                (
                    ETAP₂ * inv(dx) * inv(dy)
                    -ETA₁ * inv(dx) * inv(dy)
                    -∂RHO∂x * gy[i, j] * dt * 0.25
                ),
                kvy,
                kvx+6-6*Ny1
            ) # Vx₂
            updateindex!(L, +, Kcont*inv(dy), kvy, kpm) # P₁
            updateindex!(L, +, -Kcont*inv(dy), kvy, kpm+6) # P₂
            R[kvy] = (
                -RHOY[i, j] * gy[i, j]
                -(SXY₂-SXY₁) * inv(dx)
                -(SYY₂-SYY₁) * inv(dy)
            ) # RHS
        end # Vy equation
        # P equation
        if i==1 || i==Ny1 || j==1 || j==Nx1
            # P equation external points: boundary conditions
            # all locations: ghost unknowns P=0 -> 1.0⋅P[i,j]=0.0
            updateindex!(L, +, 1.0, kpm, kpm)
            # R[kpm] = 0.0 # already done with initialization
        # elseif i==j==2
        elseif (
            (i==2 && 2<=j<=Nx)
            || (j==2 && 2<i<Ny)
            || (i==Ny && 2<=j<=Nx)
            || (j==Nx && 2<i<Ny)
        )
            # Ptotal/Pfluid real pressure boundary condition 'anchor'
            updateindex!(L, +, Kcont, kpm, kpm)
            R[kpm] = psurface
        else
            # P equation internal points: continuity equation: ∂Vx/∂x+∂Vy/∂y=0
            #
            #                 kvy-6
            #                  Vy₁
            #                   |
            #                   |
            #      kvx-6⋅Ny1   kpm       kvx
            #        Vx₁--------P--------Vx₂
            #                   |
            #                   |
            #                  kvy
            #                  Vy₂
            #
            updateindex!(L, +, -1.0/dx, kpm, kvx-6*Ny1) # Vx₁
            updateindex!(L, +, 1.0/dx, kpm, kvx) # Vx₂
            updateindex!(L, +, -1.0/dy, kpm, kvy-6) # Vy₁
            updateindex!(L, +, 1.0/dy, kpm, kvy) # Vy₂
            # LHS coefficient matrix
            updateindex!(
                L,
                +,
                Kcont/(1-PHI[i, j]) * (inv(ETAPHI[i, j])+BETAPHI[i, j]/dt),
                kpm,
                kpm
            ) # P: Ptotal
            updateindex!(
                L,
                +,
                -Kcont/(1-PHI[i, j]) * (inv(ETAPHI[i, j])+BETAPHI[i, j]/dt),
                kpm,
                kpf
            ) # P: Pfluid
            # RHS coefficient vector
            R[kpm] = (
                (pr0[i,j]-pf0[i,j]) / (1-PHI[i,j])
                * BETAPHI[i,j]/dt
                + DMP[i, j]
            )
        end # P equation
        # qxDarcy equation
        if i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1
            # qxDarcy equation external points: boundary conditions
            # all locations: ghost unknowns qyD = 0 -> 1.0⋅qxD[i, j] = 0.0
            updateindex!(L, +, 1.0, kqx, kqx)
            # R[kqx] = 0.0 # already done with initialization
            # top boundary
            if i==1 && 1<j<Nx
                updateindex!(L, +, bcftop, kqx, kqx+6)
            end
            # bottom boundary
            if i==Ny1 && 1<j<Nx
                updateindex!(L, +, bcfbottom, kqx, kqx-6)
            end
        else
            # qxDarcy equation internal points: x-Darcy equation:
            # ηfluid/kᵩx⋅qxDarcy + ∂P/∂x = ρfluid⋅gx
            #
            #        P₁--------qxD--------P₂
            #       kpf        kqx     kpf+6⋅Ny1
            #
            # LHS coefficient matrix
            updateindex!(L, +, RX[i, j], kqx, kqx) # qxD
            updateindex!(L, +, -Kcont*inv(dx), kqx, kpf) # P₁
            updateindex!(L, +, Kcont*inv(dx), kqx, kpf+6*Ny1) # P₂
            # RHS coefficient vector
            R[kqx] = RHOFX[i, j] * gx[i, j]
        end # qxDarcy equation
        # qyDarcy equation
        if i==1 || i==Ny || i==Ny1 || j==1 || j==Nx1
            # qyDarcy equation external points: boundary conditions
            # all locations: ghost unknowns qyD = 0 -> 1.0⋅qyD[i, j] = 0.0
            updateindex!(L, +, 1.0, kqy, kqy)
            # R[kqy] = 0.0 # already done with initialization
            # left boundary
            if j==1 && 1<i<Ny
                updateindex!(L, +, bcfleft, kqy, kqy+6*Ny1)
            end
            # right boundary
            if j==Nx1 && 1<i<Ny
                updateindex!(L, +, bcfright, kqy, kqy-6*Ny1)
            end 
        else
            # qyDarcy equation internal points: y-Darcy equation:
            # ηfluid/kᵩy⋅qyDarcy + ∂P/∂y = ρfluid⋅gy
            #
            #                  P₁
            #                 kpf
            #                  |
            #                 qyD
            #                 kqy
            #                  |
            #                  P₂
            #                kpf+6
            #
            # LHS coefficient matrix
            updateindex!(L, +, RY[i, j], kqy, kqy) # qyD
            updateindex!(L, +, -Kcont*inv(dy), kqy, kpf) # P₁
            updateindex!(L, +, Kcont*inv(dy), kqy, kpf+6) # P₂
            # RHS coefficient vector
            R[kqy] = RHOFY[i, j] * gy[i, j]
        end # qyDarcy equation
        # Ptotal/Pfluid equation 
        if i==1 || i==Ny1 || j==1 || j==Nx1
            # Ptotal/Pfluid equation external points: boundary conditions
            # all locations: ghost unknowns P = 0 -> 1.0⋅P[i, j] = 0.0
            updateindex!(L, +, 1.0, kpf, kpf)
            # R[kpf] = 0.0 # already done with initialization
        # elseif i==j==2
        elseif (
            (i==2 && 2<=j<=Nx)
            || (j==2 && 2<i<Ny)
            || (i==Ny && 2<=j<=Nx)
            || (j==Nx && 2<i<Ny)
        )
            # Ptotal/Pfluid real pressure boundary condition 'anchor'
            updateindex!(L, +, Kcont, kpf, kpf)
            R[kpf] = psurface
        else
            # Ptotal/Pfluid equation internal points: continuity equation:
            # ∂qxD/∂x + ∂qyD/∂y - (Ptotal-Pfluid)/ηϕ = 0.0
            #
            #                 qyD₁
            #                kqy-6
            #                  |
            #       qxD₁-------P-------qxD₂
            #    kqx-6⋅Ny1    kpf      kqx
            #                  |
            #                 qyD₂
            #                 kqy
            #
            # LHS coefficient matrix
            updateindex!(L, +, -inv(dx), kpf, kqx-6*Ny1) # qxD₁
            updateindex!(L, +, inv(dx), kpf, kqx) # qxD₂
            updateindex!(L, +, -inv(dy), kpf, kqy-6) # qyD₁
            updateindex!(L, +, inv(dy), kpf, kqy) # qyD₂
            updateindex!(
                L,
                +,
                -Kcont/(1.0-PHI[i, j]) * (1.0/ETAPHI[i, j]+BETAPHI[i, j]/dt),
                 kpf,
                 kpm
            ) # Ptotal
            updateindex!(
                L,
                +,
                Kcont/(1.0-PHI[i, j]) * (1.0/ETAPHI[i, j]+BETAPHI[i, j]/dt),
                kpf,
                kpf
            ) # Pfluid
            # RHS coefficient vector
            R[kpf] = -(pr0[i, j]-pf0[i, j]) / (1-PHI[i, j]) * BETAPHI[i, j]/dt
        end # Ptotal/Pfluid equation
    end # for j=1:1:Nx1, i=1:1:Ny1
    end # @inbounds 
    flush!(L) # finalize CSC matrix
end # @timeit to "assemble_hydromechanical_lse()"
    return L
end # function assemble_hydromechanical_lse!

"""
Process hydromechanical solution vector to output physical observables.

$(SIGNATURES)

# Details

    - S: hydromechanical solution vector
    - vx: solid velocity at Vx nodes
    - vy: solid velocity at Vy nodes
    - pr: total pressure at P nodes
    - qxD: qx-Darcy flux at Vx nodes
    - qyD: qy-Darcy flux at Vy nodes
    - pf: fluid pressure at P nodes

# Returns

    - nothing
"""
function process_hydromechanical_solution!(
    S,
    vx,
    vy,
    pr,
    qxD,
    qyD,
    pf
)
@timeit to "process_hydromechanical_solution!()" begin
    S_mat = reshape(S, (:, Ny1, Nx1))
    @inbounds begin
        @views @. vx = S_mat[1, :, :]
        @views @. vy = S_mat[2, :, :]
        @views @. pr = S_mat[3, :, :] .* Kcont
        @views @. qxD = S_mat[4, :, :]
        @views @. qyD = S_mat[5, :, :]
        @views @. pf = S_mat[6, :, :] .* Kcont
    end # @inbounds
    # Δp = 0.25 * (pf[2, 2]+pf[2, Nx]+pf[Ny, 2]+pf[Ny, Nx]) - psurface
    # pr .-= Δp
    # pf .-= Δp
end # @timeit to "process_hydromechanical_solution!()"
    return nothing
end # function process_hydromechanical_solution!

"""
Recompute bulk viscosity at P nodes.

# Details

    - ETA: viscoplastic viscosity at basic nodes
    - ETAP: viscosity at P nodes
    - ETAPHI: bulk viscosity at P nodes
    - PHI: porosity at P nodes
    - etaphikoef: coefficient: shear viscosity -> compaction viscosity

# Returns

    - nothing
"""
function recompute_bulk_viscosity!(ETA, ETAP, ETAPHI, PHI, etaphikoef)
@timeit to "recompute_bulk_viscosity!" begin
    @inbounds begin  
        @views @. ETAP[2:end-1, 2:end-1] = 4.0 / (
            inv(ETA[1:end-1, 1:end-1]) +
            inv(ETA[2:end, 1:end-1]) +
            inv(ETA[1:end-1, 2:end]) +
            inv(ETA[2:end, 2:end])
        )
        @views @. ETAPHI = etaphikoef * ETAP * inv(PHI)
    end # @inbounds
end # @timeit to "recompute_bulk_viscosity!"
    return nothing
end

"""
Compute porosity coefficient Aϕ = Dln[(1-ϕ)/ϕ]/Dt

$(SIGNATURES)

# Details

## In

    - ETAPHI: bulk viscosity at P Nodes
    - BETAPHI: bulk compressibility at P nodes
    - PHI: porosity at P Nodes
    - pr: total pressure at P nodes
    - pf: fluid pressure at P nodes
    - pr0: previous step total pressure at P nodes
    - pf0: previous step fluid pressure at P nodes
    - dt: time step

## Out

    - APHI: porosity coefficient at P nodes

# Returns

    - aphimax: maximum absolute porosity coefficient
"""
function compute_Aϕ!(APHI, ETAPHI, BETAPHI, PHI, pr, pf, pr0, pf0, dt)
@timeit to "compute_Aϕ!()" begin
    # APHI .= 0.0
    @inbounds begin
    @views @. APHI[2:Ny, 2:Nx] = (
        ((pr[2:Ny, 2:Nx]-pf[2:Ny, 2:Nx])/ETAPHI[2:Ny, 2:Nx]
        + (
            (pr[2:Ny, 2:Nx]-pr0[2:Ny, 2:Nx])-(pf[2:Ny, 2:Nx]-pf0[2:Ny, 2:Nx])
        )/dt*BETAPHI[2:Ny, 2:Nx]) / (1-PHI[2:Ny, 2:Nx]) / PHI[2:Ny, 2:Nx]
    )
    return maximum(abs, APHI[2:Ny, 2:Nx]) # includes [2, 2] anchor abberation
    end # @inbounds
    # return maximum(abs, APHI[3:Ny-1, 3:Nx-1]) # no abberation
end # @timeit to "compute_Aϕ!()"
end # function compute_Aϕ!

"""
Compute current fluid velocities.

$(SIGNATURES)

# Details

## In

    - PHIX: porosity at Vx nodes
    - PHIY: porosity at Vy nodes
    - qxD: qx-Darcy flux at Vx nodes
    - qyD: qy-Darcy flux at Vy nodes
    - vx: solid velocity at Vx nodes
    - vy: solid velocity at Vy nodes

## Out 

    - vxf: fluid vx velocity at Vx nodes
    - vyf: fluid vy velocity at Vy nodes

# Returns

    - nothing
"""
function compute_fluid_velocities!(
    PHIX,
    PHIY,
    qxD,
    qyD,
    vx,
    vy,
    vxf,
    vyf
)
@timeit to "compute_fluid_velocities!()" begin
    @inbounds begin
        # vx velocity
        @views @. vxf[2:Ny, 1:Nx] = qxD[2:Ny, 1:Nx] / PHIX[2:Ny, 1:Nx]
        # top boundary
        @views @. vxf[1, :] = -bcftop  * vxf[2, :]
        # bottom boundary
        @views @. vxf[Ny1, :] = -bcfbottom * vxf[Ny, :]
        # vy velocity
        @views @. vyf[1:Ny, 2:Nx] = qyD[1:Ny, 2:Nx] / PHIY[1:Ny, 2:Nx]
        # left boundary
        @views @. vyf[:, 1] = -bcfleft * vyf[:, 2]
        # right boundary
        @views @. vyf[:, Nx1] = -bcfright * vyf[:, Nx]
        # adding solid velocity
        @views @. vxf += vx
        @views @. vyf += vy
    end # @inbounds

     # for j=1:1:Nx, i=2:1:Ny
    #     vxf[i, j] = qxD[i, j]*inv(PHIX[i,j]) + vx[i, j]
    # end
    # @views @. vxf[1, :] = -bcftop*vxf[2, :]    
    # @views @. vxf[Ny1, :] = -bcfbottom*vxf[Ny, :]
    # for j=2:1:Nx, i=1:1:Ny
    #     vyf[i,j] = qyD[i,j]*inv(PHIY[i,j]) + vy[i,j]
    # end
    # @views @. vyf[:, 1] = -bcfleft*vyf[:, 2]    
    # @views @. vyf[:, Nx1] = -bcfright*vyf[:, Nx]     
end # @timeit to "compute_fluid_velocities!()"
    return nothing
end # function compute_fluid_velocities!

"""
Compute velocity-/displacement-limited time step.

$(SIGNATURES)

# Details

    - vx: solid vx velocity at Vx nodes
    - vy: solid vy velocity at Vy nodes
    - vxf: fluid vx velocity at Vx nodes
    - vyf: fluid vy velocity at Vy nodes
    - dt: current time step
    - aphimax: maximum observed porosity coefficient
   
# Returns

    - dt: displacement time step
"""
function compute_displacement_timestep(
    vx,
    vy,
    vxf,
    vyf,
    dt,
    aphimax
)
@timeit to "compute_displacement_timestep()" begin
    maxvx = maximum(abs, vx)
    maxvy = maximum(abs, vy)
    maxvxf = maximum(abs, vxf)
    maxvyf = maximum(abs, vyf)    
    @info "dt before velocity limitations = $dt s"
    dt = ifelse(dt*maxvx > dxymax*dx, dxymax*dx*inv(maxvx), dt)
    @info "dt after vx limitation = $dt s"
    dt = ifelse(dt*maxvy > dxymax*dy, dxymax*dy*inv(maxvy), dt)
    @info "dt after vy limitation = $dt s"
    dt = ifelse(dt*maxvxf > dxymax*dx, dxymax*dx*inv(maxvxf), dt)
    @info "dt after vxf limitation = $dt s"
    dt = ifelse(dt*maxvyf > dxymax*dy, dxymax*dy*inv(maxvyf), dt)
    @info "dt after vyf limitation = $dt s"
    dt = ifelse(dt*aphimax > dphimax, dphimax*inv(aphimax), dt)
    @info "dt after aphimax limitation = $dt s"
end # @timeit to "compute_displacement_timestep()"
    return dt
end # function compute_displacement_timestep

"""
Compute stress, stress change, and strain rate components.

$(SIGNATURES)

# Details

## In

    - vx: solid vx velocity at Vx nodes
    - vy: solid vy velocity at Vy nodes
    - ETA: viscosity at basic nodes
    - GGG: shear modulus at basic nodes
    - ETAP: viscosity at P nodes
    - GGGP: shear modulus at P nodes
    - SXX0: previous time step σ₀′xx at P nodes
    - SXY0: previous time step σ₀xy at basic nodes
    - dt: computational time step

## Out

    - EXX: ϵxx at P nodes
    - EXY: ϵxy at basic nodes
    - SXX: σ′xx P nodes
    - SXY: σxy at basic nodes
    - DSXX: stress change Δσ′xx at P nodes
    - DSXY: stress change Δσxy at basic nodes
    - EII: second strain rate invariant ϵᴵᴵ at P nodes
    - SII: second stress invariant σᴵᴵ at P nodes

# Returns
    
        - nothing
"""
function compute_stress_strainrate!(
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
@timeit to "compute_stress_strainrate!()" begin
        @inbounds begin
            # ϵxy, σxy, Δσxy at basic nodes
            for j=1:1:Nx, i=1:1:Ny
                EXY[i, j] = 0.5 * (
                    (vx[i+1, j]-vx[i, j]) / dy
                    +(vy[i, j+1]-vy[i, j])/dx
                )
                SXY[i,j] = (
                    2*ETA[i, j]*EXY[i, j]*GGG[i,j]*dt / (
                        GGG[i, j]*dt+ETA[i, j]
                    )  + SXY0[i, j]*ETA[i, j] / (GGG[i, j]*dt+ETA[i, j])
                )
                DSXY[i, j] = SXY[i, j] - SXY0[i, j]
            end
            # ϵxx, σ′xx, Δσ'xx and Eᴵᴵ, Sᴵᴵ at P nodes
            for j=2:1:Nx, i=2:1:Ny
                EXX[i, j]= 0.5 *(
                    (vx[i, j]-vx[i, j-1]) / dx
                    -(vy[i, j]-vy[i-1, j]) / dy
                )
                SXX[i, j] = (
                    2*ETAP[i, j]*EXX[i, j]*GGGP[i, j]*dt / (
                        GGGP[i, j]*dt+ETAP[i, j]
                    ) + SXX0[i, j]*ETAP[i, j] / (GGGP[i, j]*dt+ETAP[i, j])
                )
                DSXX[i, j] = SXX[i, j] - SXX0[i, j]
                EII[i, j] = sqrt(EXX[i, j]^2 + grid_average(i-1, j-1, EXY)^2)
                SII[i, j] = sqrt(SXX[i, j]^2 + grid_average(i-1, j-1, SXY)^2)
            end
        end # @inbounds        
    # # ϵxy, σxy, Δσxy at basic nodes
    # EXY .= 0.5.*(diff(vx, dims=1)[:, 1:Nx]./dy .+ diff(vy, dims=2)[1:Ny, :]./dx)
    # @. SXY = 2*ETA*EXY*GGG*dt/(GGG*dt+ETA) + SXY0*ETA/(GGG*dt+ETA)
    # @. DSXY = SXY - SXY0
    # # ϵxx, σ′xx at P nodes
    # # @. DIVV[2:end, 2:end] = 
    # #     diff(vx, dims=2)[2:end, :]/dx + diff(vy, dims=1)[:, 2:end]/dy
    # EXX[2:Ny1, 2:Nx1] .= ( 
    #     0.5 .* (
    #         diff(vx, dims=2)[2:Ny1, :]./dx .- diff(vy, dims=1)[:, 2:Nx1]./dy
    #     )
    # )
    # @. SXX = 2.0*ETAP*EXX*GGGP*dt/(GGGP*dt+ETAP) + SXX0*ETAP/(GGGP*dt+ETAP)
    # @. DSXX = SXX - SXX0
    # @. EII[2:Ny, 2:Nx] = sqrt(
    #     EXX[2:Ny, 2:Nx]^2 + (
    #         (
    #             EXY[2:Ny, 2:Nx]
    #             +EXY[1:Ny-1,2:Nx]
    #             +EXY[2:Ny,1:Nx-1]
    #             +EXY[1:Ny-1,1:Nx-1]
    #         )/4.0
    #     )^2
    # )
    # @. SII[2:Ny, 2:Nx] = sqrt(
    #     SXX[2:Ny, 2:Nx]^2 + (
    #         (
    #             SXY[2:Ny, 2:Nx]
    #             +SXY[1:Ny-1,2:Nx]
    #             +SXY[2:Ny,1:Nx-1]
    #             +SXY[1:Ny-1,1:Nx-1]
    #         )/4.0
    #     )^2
    # )
end # @timeit to "compute_stress_strainrate!()"
    return nothing
end # function compute_stress_strainrate!

"""
Apply symmetry to P node observables.

$(SIGNATURES)

# Details

    - SXX: σ′xx at P nodes
    - APHI: Aϕ = Dln[(1-ϕ)/ϕ]/Dt at P nodes
    - PHI: porosity at P nodes
    - pr: total pressure at P nodes
    - pf: fluid pressure at P nodes
    - ps: solid pressure at P nodes

# Returns

    - nothing
"""
function symmetrize_p_node_observables!(
    SXX,
    APHI,
    PHI,
    pr,
    pf,
    ps
)
@timeit to "symmetrize_p_node_observables!()" begin
    # top boundary
    @inbounds @views @. begin
        SXX[1, 2:Nx] = SXX[2, 2:Nx]
        APHI[1, 2:Nx] = APHI[2, 2:Nx]    
        PHI[1, 2:Nx] = PHI[2, 2:Nx]    
        pr[1, 2:Nx] = pr[2, 2:Nx]    
        pf[1, 2:Nx] = pf[2, 2:Nx]    
        # bottom boundary
        SXX[Ny1, 2:Nx] = SXX[Ny, 2:Nx]
        APHI[Ny1, 2:Nx] = APHI[Ny, 2:Nx]    
        PHI[Ny1, 2:Nx] = PHI[Ny, 2:Nx]    
        pr[Ny1, 2:Nx] = pr[Ny, 2:Nx]    
        pf[Ny1, 2:Nx] = pf[Ny, 2:Nx]    
        # left boundary
        SXX[:, 1] = SXX[:, 2]
        APHI[:, 1] = APHI[:, 2]    
        PHI[:, 1] = PHI[:, 2]    
        pr[:, 1] = pr[:, 2]    
        pf[:, 1] = pf[:, 2]    
        # right boundary
        SXX[:, Nx1] = SXX[:, Nx]
        APHI[:, Nx1] = APHI[:, Nx]    
        PHI[:, Nx1] = PHI[:, Nx]    
        pr[:, Nx1] = pr[:, Nx]    
        pf[:, Nx1] = pf[:, Nx]
        # solid pressure
        ps = (pr-pf*PHI) * inv(1-PHI)
    end
end # @timeit to "symmetrize_p_node_observables!()"
    return nothing
    end # function symmetrize_p_node_observables!

"""
Compute nodal adjustment and return plastic iterations completeness status.

$(SIGNATURES)

# Details

    - ETA: viscoplastic viscosity at basic nodes
    - ETA0: previous time step viscoplastic viscosity at basic nodes
    - ETA5: plastic iterations viscoplastic viscosity at basic nodes
    - GGG: shear modulus at basic nodes
    - SXX: σ′xx at P nodes
    - SXY: σxy at basic nodes
    - pr: total pressure at P nodes
    - pf: fluid pressure at P nodes
    - COH: compressive strength at basic nodes 
    - TEN: tensile strength at basic nodes 
    - FRI: friction at basic nodes
    - YNY: plastic yielding status at basic nodes 
    - YNY5: plastic iterations plastic yielding status at basic nodes
    - YERRNOD: vector of summed yielding errors of nodes over plastic iterations
    - DSY: (SIIB-syield) at basic nodes
    - dt: time set
    - iplast: plastic iteration step 

# Returns

    - plastic_iterations_complete: true if plastic iterations complete
"""
function compute_nodal_adjustment!(
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
@timeit to "compute_nodal_adjustment!()" begin
    # reset / setup
    ETA5 .= ETA0
    YNY5 .= 0
    DSY .= 0.0
    ynpl = 0
    ddd = 0.0
    @inbounds begin
    for j=1:1:Nx, i=1:1:Ny
        # second stress invariant at basic nodes
        SIIB = sqrt(SXY[i, j]^2 + grid_average(i, j, SXX)^2) 
        # second invariant for purely elastic stress buildup at basic nodes
        siiel = SIIB * (GGG[i, j]*dt+ETA[i, j]) / ETA[i, j]
        # interpolate total and fluid pressure at basic nodes
        prB = grid_average(i, j, pr)
        pfB = grid_average(i, j, pf)
        # yielding stress: confined fracture
        syieldc = COH[i, j] + FRI[i, j] * (prB-pfB)
        # yielding stress: tensile fracture
        syieldt = TEN[i, j] + (prB-pfB)
        # non-negative yielding stress requirement
        syield = max(min(syieldc, syieldt), 0.0)
        # update error for previous yielding nodes
        ynn = false
        if YNY[i, j] > 0
            DSY[i, j] = SIIB - syield
            ddd += DSY[i, j]^2
            ynpl += 1
        end
        # correcting viscosity for yielding
        if syield < siiel
            # update viscosity for basic node
            etapl = dt * GGG[i,j]*syield/(siiel-syield)
            if etapl < ETA0[i, j]
                # recompute nodal viscosity, apply min/max viscosity cutoffs
                ETA5[i, j] = etapl^(1.0-etawt) * ETA[i, j]^etawt
                if ETA5[i, j] > etamax
                    ETA5[i, j] = etamax
                    elseif ETA5[i, j] < etamin
                        ETA5[i, j] = etamin
                end
                # mark yielding nodes
                YNY5[i, j] = 1
                # update error for new yielding nodes
                if ynn == false
                    DSY[i, j] = SIIB - syield
                    ddd += DSY[i, j]^2
                    ynpl += 1
                end  
            end
        end
    end
    if ynpl > 0
        YERRNOD[iplast] = sqrt(ddd/ynpl)
    end
    # return plastic iteration completeness
    @info "end plastic iter $iplast: ynpl=$ynpl, YERRNOD=$(YERRNOD[iplast])"
    return ynpl==0 || YERRNOD[iplast]<yerrmax || iplast==nplast
    end # @inbounds
end # @timeit to "compute_nodal_adjustment!()
end # function compute_nodal_adjustment!

"""
Compare two arrays of identical sizes element-wise and fill a third array with the larger value if positive and zero otherwise.

$(SIGNATURES)

# Details

    - A: first array
    - B: second array
    - C: result array

# Returns

    - nothing
"""
function positive_max!(A, B, C)
@timeit to "positive_max!()" begin
    @inbounds for i in eachindex(A)
        C[i] = max(0, ifelse(A[i] > B[i], A[i], B[i]))
    end
end # @timeit to "positive_max!()"
    return nothing
end # function positive_max

"""
Decide next pass plastic iteration time step, viscoplastic viscosity,
and basic node yielding status.

$(SIGNATURES)

# Details:

    - ETA: viscoplastic viscosity at basic nodes
    - ETA5: plastic iterations viscoplastic viscosity at basic nodes
    - ETA00: previous time step viscoplastic viscosity at basic nodes
    - YNY: plastic yielding status at basic nodes 
    - YNY5: plastic iterations plastic yielding status at basic nodes
    - YNY00: previous time step plastic yielding status at basic nodes
    - YNY_inv_ETA: inverse of plastic viscosity at yielding basic nodes
    - dt: current time step
    - iplast: current plastic iteration counter

# Returns

    - dt: adjusted next time step
"""
function finalize_plastic_iteration_pass!(
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
@timeit to "finalize_plastic_iteration_pass!()" begin
    if iplast % dtstep == 0
        # dtstep plastic iterations performed without reaching targets:
        # decrease time step and reset to previous viscoplastic viscosity
        dt *= dtcoefdn
        @info "reducing dt due to plastic iteration limit: dt=$dt s"
        ETA .= ETA00
        YNY .= YNY00
    else
        # perform next plastic iteration pass with new viscoplastic viscosity
        ETA .= ETA5
        YNY .= YNY5
    end
    @views @. YNY_inv_ETA = YNY / ETA
    return dt
end # @timeit to "finalize_plastic_iteration_pass!()"
end # function finalize_plastic_iteration_pass

"""
Decide next pass thermochemical iteration time step.

$(SIGNATURES)

# Details:

    - maxDTcurrent: maximum temperature difference between current and
                    previous time step
    - dt: current time step duration
    - titer: current thermochemical iteration counter

# Returns

    - dt: adjusted next time step
"""
function finalize_thermochemical_iteration_pass(maxDTcurrent, dt, titer)
@timeit to "finalize_thermochemical_iteration_pass" begin
    if titer == 1
        if maxDTcurrent > DTmax
            dt *= (DTmax * inv(maxDTcurrent))
            @info "titer 1: reducing dt due to maxDT: dt=$dt s"
        end
    end
end # @timeit to "finalize_thermochemical_iteration_pass"
    return dt
end # function finalize_thermochemical_iteration_pass

"""
Assess outcome of thermochemical iteration and return thermochemical iterations
completeness status.

$(SIGNATURES)

# Details:

    - DMP: mass transfer term at P nodes
    - pf: fluid pressure at P nodes
    - pf0: previous time step fluid pressure at P nodes
    - titer: current thermochemical iteration counter

# Returns

    - dt: adjusted next time step
"""
function compute_thermochemical_iteration_outcome(DMP, pf, pf0, titer)
@timeit to "compute_thermochemical_iteration_outcome" begin
    pferrcur = maximum(abs, pf-pf0)
    DMPmax = maximum(abs, DMP)
    @info "end thermochemical iter $titer" pferrcur DMPmax
    return pferrcur < pferrmax && (titer>2||DMPmax<=0.0)
end # @timeit to "compute_thermochemical_iteration_outcome"
end # function compute_thermochemical_iteration_outcome

"""
Apply xy (basic grid) and xx (P grid) subgrid stress diffusion to markers.

$(SIGNATURES)

# Details

    - xm: marker x-coordinates
    - ym: marker y-coordinates
    - tm: marker type
    - inv_gggtotalm: inverse of total shear modulus of markers
    - sxxm: marker σ′xx [Pa]
    - sxym: marker σxy [Pa]
    - SXX0: σ₀′xx at P nodes [1/s]
    - SXY0: σ₀xy at basic nodes [1/s]
    - DSXX: stress change Δσ′xx at P nodes
    - DSXY: stress change Δσxy at basic nodes
    - SXXSUM: interpolation of SXX at P nodes
    - SXYSUM: interpolation of SXY at basic nodes
    - WTPSUM: interpolation weights at P nodes
    - WTSUM: interpolation weights at basic nodes
    - dt: computational time step
    - marknum: total number of markers in use

# Returns
    
    - nothing
"""
function apply_subgrid_stress_diffusion!(
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
@timeit to "apply_subgrid_stress_diffusion!" begin
    # only perform subgrid stress diffusion if enabled by dsubgrids > 0
    if dsubgrids == 0.0
        return nothing
    end 
    # fix etam[tm[m]] RMK: It's a temporary fix, not yet implemented in source
    etam = @SVector ones(3)
    # reset interpolation arrays
    SXXSUM .= 0.0
    WTPSUM .= 0.0
    SXYSUM .= 0.0
    WTSUM .= 0.0
    # iterate over markers
    for m=1:1:marknum
        @inbounds i_p, j_p, weights_p = fix_weights(
            xm[m], ym[m], xp, yp, dx, dy, jmin_p, jmax_p, imin_p, imax_p)
        @inbounds i_basic, j_basic, weights_basic = fix_weights(
            xm[m],
            ym[m],
            x,
            y,
            dx,
            dy,
            jmin_basic,
            jmax_basic,
            imin_basic,
            imax_basic
        )
        # σ₀′xx at P nodes
        # compute marker-node σ′xx difference
        @inbounds δσxxm₀ = sxxm[m] - dot4(grid_vector(i_p, j_p, SXX0), weights_p)
        # time-relax σ′xx difference
        @inbounds δσxxm₀ *= (
            exp(-dsubgrids*dt/(etam[tm[m]]*inv_gggtotalm[m])) - 1.0) 
        # correct marker stress
        @inbounds sxxm[m] += δσxxm₀
        # update subgrid diffusion on P nodes
        interpolate_add_to_grid!(i_p, j_p, weights_p, δσxxm₀, SXXSUM)
        interpolate_add_to_grid!(i_p, j_p, weights_p, 1.0, WTPSUM)
        # σ₀xy at basic nodes
        # compute marker-node σxy difference
        @inbounds δσxy₀ = sxym[m] - dot4(
            grid_vector(i_basic, j_basic, SXY0), weights_basic)
        # time-relax σxy difference
        @inbounds δσxy₀ *= (
            exp(-dsubgrids*dt/(etam[tm[m]]*inv_gggtotalm[m])) - 1.0)
        # correct marker stress
        @inbounds sxym[m] += δσxy₀
        # update subgrid diffusion on basic nodes
        interpolate_add_to_grid!(i_basic, j_basic, weights_basic, δσxy₀, SXYSUM)
        interpolate_add_to_grid!(i_basic, j_basic, weights_basic, 1.0, WTSUM)
    end
    # compute DSXXsubgrid and update DSXX at inner P nodes
    @views @. DSXX[2:Ny, 2:Nx][WTPSUM[2:Ny, 2:Nx]>0.0] -= (
        SXXSUM[2:Ny, 2:Nx][WTPSUM[2:Ny, 2:Nx]>0.0] /
        WTPSUM[2:Ny, 2:Nx][WTPSUM[2:Ny, 2:Nx]>0.0]
    )
    # compute DSXYsubgrid and update DSXY at all basic nodes
    @views @. DSXY[WTSUM[:, :]>0.0] -= SXYSUM[:, :][WTSUM[:, :]>0.0] /
        WTSUM[:, :][WTSUM[:, :]>0.0]
end # @timeit to "apply_subgrid_stress_diffusion!"
    return nothing
end # function apply_subgrid_stress_diffusion!

"""
Update marker stress based on xy (basic grid) and xx (P grid) stress changes.

$(SIGNATURES)

# Details

    - xm: x-coordinates of markers
    - ym: y-coordinates of markers
    - sxxm: marker σ′xx [Pa]
    - sxym: marker σxy [Pa]
    - DSXX: stress change Δσ′xx at P nodes
    - DSXY: stress change Δσxy at basic nodes
    - marknum: total number of markers in use

# Returns

    - nothing
"""
function update_marker_stress!(xm, ym, sxxm, sxym, DSXX, DSXY, marknum)
@timeit to "update_marker_stress!" begin
    @threads for m=1:1:marknum    
        @inbounds i_p, j_p, weights_p = fix_weights(
            xm[m],
            ym[m],
            xp,
            yp,
            dx,
            dy,
            2,
            Nx-1,
            2,
            Ny-1
        )
        @inbounds i_basic, j_basic, weights_basic = fix_weights(
            xm[m],
            ym[m],
            x,
            y,
            dx,
            dy,
            1,
            Nx-1,
            1,
            Ny-1
        )
    # interpolate updated DSXX, DSXY back to markers
        interpolate_add_to_marker!(m, i_p, j_p, weights_p, sxxm, DSXX)
        interpolate_add_to_marker!(
            m, i_basic, j_basic, weights_basic, sxym, DSXY)
    end
end # @timeit to "update_marker_stress!"
    return nothing
end # function update_marker_stress!

"""
Compute shear heating based on basic (temperature) and P grids.

$(SIGNATURES)

# Details

    - HS: shear heating
    - ETA: viscoplastic viscosity at basic nodes
    - SXY: σ₀xy XY stress at basic nodes
    - ETAP: viscosity at P nodes
    - SXX: σ₀xy XY stress at basic nodes
    - RX: ηfluid/Kϕ at Vx nodes
    - RY: ηfluid/Kϕ at Vy nodes
    - qxD: qx-Darcy flux at Vx nodes
    - qyD: qy-Darcy flux at Vy nodes
    - PHI: porosity at P nodes
    - ETAPHI: bulk viscosity at P nodes
    - pr: total pressure at P nodes
    - pf: fluid pressure at P nodes

# Returns

    - nothing
"""
function compute_shear_heating!(
    HS, ETA, SXY, ETAP, SXX, RX, RY, qxD, qyD, PHI, ETAPHI, pr, pf)
@timeit to "compute_shear_heating!" begin
    for j=2:1:Nx, i=2:1:Ny
        # average SXY⋅EXY
        @timeit to "average SXY⋅EXY" SXYEXY = 0.25 * sum(
            grid_vector(i-1, j-1, SXY).^2 ./ grid_vector(i-1, j-1, ETA))
        # compute shear heating HS
        @timeit to "compute shear heating HS" @inbounds HS[i, j] = (
            SXX[i, j]^2 / ETAP[i,j]
            + SXYEXY
            + (pr[i, j]-pf[i, j])^2 / (1-PHI[i, j]) / ETAPHI[i, j]
            + 0.5 * (RX[i, j-1]*qxD[i, j-1]^2 + RX[i, j]*qxD[i, j]^2)
            + 0.5 * (RY[i-1, j]*qyD[i-1, j]^2 + RY[i, j]*qyD[i, j]^2)
        )
    end
end # @timeit to "compute_shear_heating!" 
    return nothing
end # function compute_shear_heating!

"""
Compute adiabatic heating based on basic (temperature) and P grids.

$(SIGNATURES)

# Details

    - HA: adiabatic heating at P nodes
    - tk1: previous temperature at P nodes
    - ALPHA: thermal expansion coefficient at P nodes
    - ALPHAF: fluid thermal expansion coefficient at P nodes
    - PHI: porosity at P nodes
    - vx: solid vx-velocity at Vx nodes
    - vy: solid vy-velocity at Vy nodes
    - vxf: fluid vx-velocity at Vx nodes
    - vyf: fluid vy-velocity at Vy nodes
    - ps: solid pressure at P nodes
    - pf: fluid pressure at P nodes

# Returns

    - nothing
"""
function compute_adiabatic_heating!(
    HA, tk1, ALPHA, ALPHAF, PHI, vx, vy, vxf, vyf, ps, pf)
@timeit to "compute_adiabatic_heating!" begin
    @inbounds begin
        for j=2:1:Nx, i=2:1:Ny
            # indirect calculation of DP/Dt ≈ (∂P/∂x)⋅vx + (∂P/∂y)⋅vy (eq. 9.23)
            # average vy, vx, vxf, vyf
            VXP = 0.5 * (vx[i, j]+vx[i, j-1])
            VYP = 0.5 * (vy[i, j]+vy[i-1, j])
            VXFP = 0.5 * (vxf[i, j]+vxf[i, j-1])
            VYFP = 0.5 * (vyf[i, j]+vyf[i-1, j])
            # evaluate DPsolid/Dt with upwind differences
            if VXP < 0.0
                dpsdx = (ps[i, j]-ps[i, j-1]) * inv(dx)
            else
                dpsdx = (ps[i, j+1]-ps[i, j]) * inv(dx)
            end
            if VYP < 0.0
                dpsdy = (ps[i, j]-ps[i-1, j]) * inv(dy)
            else
                dpsdy = (ps[i+1, j]-ps[i, j]) * inv(dy)
            end
            dpsdt = VXP*dpsdx + VYP*dpsdy
            # evaluate DPfluid/Dt with upwind differences
            if VXFP > 0.0
                dpfdx = (pf[i, j]-pf[i, j-1]) * inv(dx)
            else
                dpfdx = (pf[i, j+1]-pf[i, j]) * inv(dx)
            end
            if VYFP > 0.0
                dpfdy = (pf[i, j]-pf[i-1, j]) * inv(dy)
            else
                dpfdy = (pf[i+1, j]-pf[i, j]) * inv(dy)
            end
            dpfdt = VXFP*dpsdx + VYFP*dpsdy
            # Hₐ = (1-ϕ)Tαˢ⋅DPˢ/Dt + ϕTαᶠ⋅DPᶠ/Dt (eq. 9.23)
            HA[i, j] = (
                (1-PHI[i, j]) * tk1[i, j] * ALPHA[i, j] * dpsdt
                + PHI[i, j] * tk1[i, j] * ALPHAF[i, j] * dpfdt
            )
        end
    end # @inbounds
end # @timeit to "compute_adiabatic_heating!"
end # function compute_adiabatic_heating!

"""
Assemble the LHS sparse coefficient matrix and fill RHS coefficient vector
of the energy conservation (heat) equation.

$(SIGNATURES)

# Details

	- tk1: current temperature at P nodes
	- RHOCP: volumetric heat capacity at P nodes  
	- KX: thermal conductivity at Vx nodes
	- KY: thermal conductivity at Vy nodes 
	- HR: radioactive heating at P nodes
	- HA: adiabatic heating at P nodes 
	- HS: shear heating at P nodes
    - DHP: latent heating (HL) at P nodes
    - RT: thermal RHS coefficient vector
	- dt: current time step length

# Returns

    - LT: LHS sparse coefficient matrix
"""
function assemble_thermal_lse!(tk1, RHOCP, KX, KY, HR, HA, HS, DHP, RT, dt)
@timeit to "assemble_thermal_lse!" begin
    # fresh LHS coefficient matrix
    LT = ExtendableSparseMatrix(Ny1*Nx1, Ny1*Nx1)
    # reset RHS coefficient vector
    RT .= zero(0.0)
    # compose global thermal matrix LT and coefficient vector RT
    @inbounds begin
        for j=1:1:Nx1, i=1:1:Ny1
            # define global index in algebraic space
            gk = (j-1)*Ny1 + i
            # External points
            if i==1 || i==Ny1 || j==1 || j==Nx1
                # thermal equation external points: boundary conditions
                # all locations: ghost unknowns T₃=0 -> 1.0⋅T[i,j]=0.0
                updateindex!(LT, +, 1.0, gk, gk)
                # R[gk] = 0.0 # already done with initialization
                # left boundary: ∂T/∂x=0
                if j == 1
                    updateindex!(LT, +, -1.0, gk, gk+Ny1)
                end
                # right boundary: ∂T/∂x=0
                if j == Nx1
                    updateindex!(LT, +, -1.0, gk, gk-Ny1)
                end
                # top inner boundary: ∂T/∂y=0
                if i==1 && 1<j<Nx1 
                    updateindex!(LT, +, -1.0, gk, gk+1)
                end
                # bottom inner boundary: ∂T/∂y=0
                if i==Ny1 && 1<j<Nx1 
                    updateindex!(LT, +, -1.0, gk, gk-1)
                end
            else
                # internal points: 2D thermal equation (conservative formulation)
                # ρCₚ∂/∂t = k∇²T + Hᵣ + Hₛ + Hₐ + Hₗ
                #         = -∂qᵢ/∂xᵢ + Hᵣ + Hₛ + Hₐ + Hₗ (16.54, 16.97) 
                # in case of purely advective heat transport: ∂T/∂t+⃗v⋅∇T = 0
                #
                #                      gk-1
                #                        T₂
                #                        |
                #                       i-1,j
                #                       Ky₁
                #                       qy₁
                #                        |
                #         gk-Ny1 i,j-1  gk     i,j  gk+Ny1
                #           T₁----Kx₁----T₃----Kx₂----T₅
                #                 qx₁    |     qx₂
                #                       i,j
                #                       Ky₂
                #                       qy₂
                #                        |
                #                      gk+1
                #                        T₄
                #
                # extract thermal conductivities
                Kx₁ = KX[i, j-1] 
                Kx₂ = KX[i, j]
                Ky₁ = KY[i-1, j]
                Ky₂ = KY[i, j]
                # fill system of equations: LHS
                updateindex!(LT, +, -Kx₁*inv(dx^2), gk, gk-Ny1) # T₁
                updateindex!(LT, +, -Ky₁*inv(dy^2), gk, gk-1) # T₂
                updateindex!(LT, +, (
                    RHOCP[i, j]/dt+(Kx₁+Kx₂)*inv(dx^2)+(Ky₁+Ky₂)*inv(dy^2)),
                    gk,
                    gk
                ) # T₃
                updateindex!(LT, +, -Ky₂*inv(dy^2), gk, gk+1) # T₄
                updateindex!(LT, +, -Kx₂*inv(dx^2), gk, gk+Ny1) # T₅
                # fill system of equations: RHS
                RT[gk] = (
                    RHOCP[i, j]/dt*tk1[i, j]
                    + HR[i, j]
                    + HA[i, j]
                    + HS[i, j]
                    + DHP[i, j]
                )
            end
        end
    end # @inbounds
    flush!(LT) # finalize CSC matrix
end # @timeit to "assemble_thermal_lse!"
    return LT
end # function assemble_thermal_lse!

"""
Perform thermal iterations to time step thermal field at P nodes.

$(SIGNATURES)

# Details

    - tk0: previous temperature at P nodes 
	- tk1: current temperature at P nodes
	- tk2: next temperature at P nodes 
	- DT: calculated temperature difference at P nodes 
	- DT0: previous calculated temperature difference at P nodes
	- RHOCP: volumetric heat capacity at P nodes  
	- KX: thermal conductivity at Vx nodes
	- KY: thermal conductivity at Vy nodes 
	- HR: radioactive heating at P nodes
	- HA: adiabatic heating at P nodes 
	- HS: shear heating at P nodes
    - DHP: latent heating (HL) at P nodes
    - RT: thermal RHS coefficient vector
    - ST: thermal solution vector
	- dt: computational time step

# Returns

    - nothing
"""
function perform_thermal_iterations!(
    tk0, tk1, tk2, DT, DT0, RHOCP, KX, KY, HR, HA, HS, DHP, RT, ST, dt)
@timeit to "perform_thermal_iterations!" begin
    # set up thermal iterations
    tk0 .= tk1
    dtt = dt
    dttsum = 0.0
    titer = 1
    # perform thermal iterations until reaching time limit
    while dttsum < dt
        # fresh LHS coefficient matrix
        LT = assemble_thermal_lse!(tk1, RHOCP, KX, KY, HR, HA, HS, DHP, RT, dtt)
        # solve system of equations
        ST .= LT \ RT # implicit: flush!(LT)
        # reshape solution vector to 2D array
        tk2 .= reshape(ST, Ny1, Nx1)
        # compute ΔT
        DT .= tk2 .- tk1
        if titer == 1
            # during first thermal iteration pass:
            # apply thermal timestepping stability condition
            maxDTcurrent = maximum(abs, DT)
            if maxDTcurrent > DTmax
                dtt *= DTmax * inv(maxDTcurrent)
            else
                dttsum += dtt
            end
        else
            # second+ thermal iteration passes:
            # update dttsum and adjust timestep
            dttsum += dtt
            dtt = min(dtt, dt-dttsum)
        end
        # increase thermal iteration counter
        titer += 1
    end
    # finalize overall temperature change and advance temperature field
    DT .= tk2 .- tk0
    DT0 .= DT
end # @timeit to "perform_thermal_iterations!"
    return nothing
end # function perform_thermal_iterations!

"""
Apply subgrid temperature diffusion to markers.

$(SIGNATURES)

# Details

    - xm: marker x-coordinates 
	- ym: marker y-coordinates
	- tm: marker type
	- tkm: marker temperature
    - phim: marker porosity
    - tk1: current temperature at P nodes
	- DT: temperature change at P nodes
	- TKSUM: interpolation of TK at P nodes 
	- RHOCPSUM: interpolation of RHOCP at P nodes
    - dt: computational time step
	- marknum: total number of markers in use
    - mode: marker property computation mode
        - 1: dynamic, based on (Touloukian, 1970; Hobbs, 1974;
             Travis and Schubert, 2005)
        - 2: constant parameter rhocpfluidm

# Returns

    - nothing
"""
function apply_subgrid_temperature_diffusion!(
    xm, ym, tm, tkm, phim, tk1, DT, TKSUM, RHOCPSUM, dt, marknum, mode)
@timeit to "apply_subgrid_temperature_diffusion!" begin
    # only perform subgrid temperature diffusion if enabled by dsubgridt > 0
    if dsubgridt == 0.0
        return nothing
    end
    # reset interpolation arrays
    TKSUM .= 0.0
    RHOCPSUM .= 0.0
    # iterate over markers
    for m=1:1:marknum
        @inbounds i, j, weights = fix_weights(
            xm[m], ym[m], xp, yp, dx, dy, jmin_p, jmax_p, imin_p, imax_p)
        # compute marker-node temperature difference
        @inbounds δtkm = tkm[m] - dot4(grid_vector(i, j, tk1), weights)
        # compute marker properties
        @inbounds if tm[m] < 3
            # rocks
            @inbounds rhocptotalm = total(
                rhocpsolidm[tm[m]], compute_rhocpfluidm(tkm[m], mode), phim[m])
            @inbounds ktotalm[m] = ktotal(
                compute_ksolidm(tkm[m], mode),
                compute_kfluidm(tkm[m], mode),
                phim[m]
            )
        else
            # sticky air
            @inbounds rhocptotalm = rhocpsolidm[tm[m]]
            @inbounds ktotalm = ksolidm[tm[m]]
        end
        # time-relax δtkm difference
        δtkm *= (
            exp(
                -dsubgridt*ktotalm*dt/rhocptotalm*( 2.0*(inv(dx^2)+inv(dy^2)))
            ) - 1.0
        )
        # correct marker temperature
        @inbounds tkm[m] += δtkm
        # update subgrid diffusion on P nodes
        interpolate_add_to_grid!(i, j, weights, δtkm*rhocptotalm, TKSUM)
        interpolate_add_to_grid!(i, j, weights, rhocptotalm, RHOCPSUM)
    end
    # compute DTsubgrid=TKSUM/RHOCPSUM and update temperature field at P nodes
    for j=1:1:Nx1, i=1:1:Ny1
        if RHOCPSUM[i, j] > 0.0
            DT[i, j] -= TKSUM[i, j] / RHOCPSUM[i, j]
        end
    end
end # @timeit to "apply_subgrid_temperature_diffusion!"
    return nothing
end # function apply_subgrid_temperature_diffusion! 

"""
Update marker temperature based on P grid temperature changes.

$(SIGNATURES)

# Details

    - xm: x-coordinates of markers
    - ym: y-coordinates of markers
    - tkm: marker temperature
    - DT: temperature change at P nodes
    - tk2: next temperature at P nodes
    - timestep: current time step
    - marknum: total number of markers in use

# Returns

    - nothing
# """
function update_marker_temperature!(xm, ym, tkm, DT, tk2, timestep, marknum)
@timeit to "update_marker_temperature!" begin
    if timestep == 1
        # interpolate tk2 to markers instead of DT for first time step        
        @threads for m=1:1:marknum
            @inbounds i, j, weights = fix_weights(
                xm[m], ym[m], xp, yp, dx, dy, jmin_p, jmax_p, imin_p, imax_p)
            interpolate_to_marker!(m, i, j, weights, tkm, tk2)
        end
    else
        # interpolate and apply DT to markers for subsequent time steps
        @threads for m=1:1:marknum
            @inbounds i, j, weights = fix_weights(
                xm[m], ym[m], xp, yp, dx, dy, jmin_p, jmax_p, imin_p, imax_p)
            interpolate_add_to_marker!(m, i, j, weights, tkm, DT)
        end
    end
end # @timeit to "update_marker_temperature!"
    return nothing
end # function update_marker_temperature!

"""
Update marker porosity for compaction based on Dln[(1-ϕ)/ϕ]/Dt at P grid.

$(SIGNATURES)

# Details

    - xm: x-coordinates of markers
    - ym: y-coordinates of markers
    - tm: marker type
    - phim: marker porosity
    - APHI: Dln[(1-ϕ)/ϕ]/Dt at P nodes
    - dt: computational time step
    - marknum: total number of markers in use

# Returns

    - nothing
# """
function update_marker_porosity!(xm, ym, tm, phim, APHI, dt, marknum)
@timeit to "update_marker_porosity!" begin
    # update porosity for compaction
    @inbounds begin
        @threads for m=1:1:marknum
            if tm[m] < 3
                # rocks
                i, j, weights = fix_weights(
                    xm[m],
                    ym[m],
                    xp,
                    yp,
                    dx,
                    dy,
                    jmin_p,
                    jmax_p,
                    imin_p,
                    imax_p
                )
                # compute Dln[(1-ϕ)/ϕ]/Dt at marker
                aphim = dot4(grid_vector(i, j, APHI), weights)
                # update marker porosity
                phim[m] = max(
                    phimin,
                    min(
                        phimax,
                        phim[m] / ((1.0-phim[m])*exp(aphim*dt) + phim[m])
                    )
                )
            end
        end
    end # @inbounds
end # @timeit to "update_marker_porosity!"
    return nothing
end # function update_marker_porosity!

""" 
Compute solid velocities, fluid velocities at P nodes.

$(SIGNATURES)

# Details

    - vx: solid vx-velocity at Vx nodes
    - vy: solid vy-velocity at Vy nodes
    - vxf: fluid vx-velocity at Vx nodes
    - vyf: fluid vy-velocity at Vy nodes
    - vxp: solid vx-velocity at P nodes
    - vyp: solid vy-velocity at P nodes
    - vxpf: fluid vx-velocity at P nodes
    - vypf: fluid vy-velocity at P nodes

# Returns

    - nothing
"""
function compute_velocities!(vx, vy, vxf, vyf, vxp, vyp, vxpf, vypf)
@timeit to "compute_velocities!" begin
    @inbounds begin
        # compute solid velocities at P nodes
        for j=2:1:Nx, i=2:1:Ny
            vxp[i, j] = 0.5 * (vx[i, j] + vx[i, j-1]) 
            vyp[i, j] = 0.5 * (vy[i, j] + vy[i-1, j])
            vxpf[i, j] = 0.5 * (vxf[i, j] + vxf[i, j-1]) 
            vypf[i, j] = 0.5 * (vyf[i, j] + vyf[i-1, j])
        end
        # apply boundary conditions
        # vxp
        # top: free slip
        @views @. vxp[1, 2:Nx-1] =- bctop * vxp[2, 2:Nx-1]    
        # bottom: free slip
        @views @. vxp[Ny1, 2:Nx-1] =- bcbottom * vxp[Ny, 2:Nx-1]    
        # left
        @views @. vxp[:, 1] = 2.0*vxleft - vxp[:, 2]
        # right
        @views @. vxp[:, Nx1] = 2.0*vxright - vxp[:, Nx]
        # vyp
        # left: free slip
        @views @. vyp[2:Ny-1, 1] =- bcleft * vyp[2:Ny-1, 2]    
        # right: free slip
        @views @. vyp[2:Ny-1, Nx1] =- bcright * vyp[2:Ny-1, Nx]
        # top
        @views @. vyp[1, :] = 2.0*vytop - vyp[2, :]
        # bottom
        @views @. vyp[Ny1, :] = 2.0*vybottom - vyp[Ny, :]
        # vxpf
        # top: free slip
        @views @. vxpf[1, 2:Nx-1] =- bcftop * vxpf[2, 2:Nx-1]    
        # bottom: free slip
        @views @. vxpf[Ny1, 2:Nx-1] =- bcfbottom * vxpf[Ny, 2:Nx-1]    
        # left
        @views @. vxpf[:,1] = 2.0*vxleft - vxpf[:, 2]
        # right
        @views @. vxpf[:, Nx1] = 2.0*vxright - vxpf[:, Nx]
        # vypf
        # left: free slip
        @views @. vypf[2:Ny-1, 1] =- bcfleft * vypf[2:Ny-1, 2]    
        # right: free slip
        @views @. vypf[2:Ny-1, Nx1] =- bcfright * vypf[2:Ny-1, Nx]
        # top
        @views @. vypf[1, :] = 2.0*vytop - vypf[2, :]
        # bottom
        @views @. vypf[Ny1,:] = 2.0*vybottom - vypf[Ny, :]
    end # @inbounds
end # @timeit to "compute_velocities!"
    return nothing
end # function compute_velocities!

"""
Compute rotatation rate in basic nodes based on velocity derivatives
at Vx and Vy nodes.

$(SIGNATURES)

# Details

    - vx: solid vx-velocity at Vx nodes
    - vy: solid vy-velocity at Vy nodes
    - wyx: rotation rate at basic nodes
  
# Returns

    - nothing
"""
function compute_rotation_rate!(vx, vy, wyx)
@timeit to "compute_rotation_rate!" begin
    # compute rotation rate ωyx=1/2[∂Vy/∂x-∂Vx/∂y] at basic nodes
    for j=1:1:Nx, i=1:1:Ny
        @inbounds wyx[i, j] = 0.5 * (
        (vy[i, j+1]-vy[i, j])*inv(dx) - (vx[i+1, j]-vx[i, j])*inv(dy)
        )
    end
end # @timeit to "compute_rotation_rate!"
    return nothing
end # function compute_rotation_rate!

"""
Move markers using classic Runge-Kutta integration (RK4) taking into account
two-phase flow and solid-fluid temperatures.

$(SIGNATURES)

# Details

    - xm: x-coordinate of markers
    - ym: y-coordinate of markers
    - tm: type of markers 
    - tkm: temperature of markers
    - phim: porosity of markers  
    - sxxm: marker σ′xx
    - sxym: marker σxy
    - vx: solid vx-velocity at Vx nodes
    - vy: solid vy-velocity at Vy nodes
    - vxf: fluid vx-velocity at Vx nodes
    - vyf: fluid vy-velocity at Vy nodes
    - wyx: rotation rate at basic nodes
    - tk2: next temperature at P nodes
    - marknum: number of markers in use
    - dt: computational time step
    - mode: marker property computation mode
        - 1: dynamic, based on (Touloukian, 1970; Hobbs, 1974;
             Travis and Schubert, 2005)
        - 2: constant parameter rhocpfluidm

# Returns

    - nothing
"""
function move_markers_rk4!(
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
    mode
)
@timeit to "move_markers_rk4!" begin
    @inbounds begin
        for m=1:1:marknum
            xmm = xmrk4 = xm[m]
            ymm = ymrk4 = ym[m]        
            i, j, weights = fix_weights(
                xmrk4,
                ymrk4,
                xp,
                yp,
                dx,
                dy,
                jmin_p,
                jmax_p,
                imin_p,
                imax_p
            )
            # interpolate local solid temperature from P grid
            tksm₀ = dot4(grid_vector(i, j, tk2), weights)
            i, j, weights = fix_weights(
                xmrk4,
                ymrk4,
                x,
                y,
                dx,
                dy,
                jmin_basic,
                jmax_basic,
                imin_basic,
                imax_basic
            )
            # interpolate local rotation rate from basic grid
            ωm = dot4(grid_vector(i, j, wyx), weights)
            # incremental rotation angle
            θ = dt * ωm
            # compute analytic stress rotation using σ′′xx = -σ′′yy
            sxxm₁ = sxxm[m]*cos(θ)^2 - sxxm[m]*sin(θ)^2-sxym[m]*sin(2.0*θ)
            sxym₁ = sxxm[m]*sin(2.0*θ) + sxym[m]*cos(2.0*θ)
            # update stresses
            sxxm[m] = sxxm₁
            sxym[m] = sxym₁
            # setup RK4 scheme 
            # RK4 coordinate positions A, B, C, D
            # RK4 velocities va, vb, vc, vd
            vxrk4 = @SVector zeros(4)
            vyrk4 = @SVector zeros(4)
            # advance marker using RK4 scheme on solid velocity
            for rk=1:1:4
                # interpolate vx
                i, j, dxmj, dymi = fix_distances(
                    xmrk4,
                    ymrk4,
                    xvx,
                    yvx,
                    dx,
                    dy,
                    jmin_vx,
                    jmax_vx,
                    imin_vx,
                    imax_vx
                )
                # compute vx velocity for left and right of current cell
                vxm₁₃ = vx[i, j]*(1.0-dxmj/dx) + vx[i, j+1]*dxmj/dx
                vxm₂₄ = vx[i+1, j]*(1.0-dxmj/dx) + vx[i+1, j+1]*dxmj/dx
                # compute second order vx velocity corrections
                if dxmj/dx >= 0.5 
                    # in right half of cell but not at right edge of grid
                    if j < Nx-1
                        vxm₁₃ += 0.5*((dxmj/dx-0.5)^2) * (
                            vx[i, j] - 2.0*vx[i, j+1] + vx[i, j+2])
                        vxm₂₄ += 0.5*((dxmj/dx-0.5)^2) * (
                            vx[i+1, j] - 2.0*vx[i+1, j+1] + vx[i+1, j+2])
                    end
                else
                    # in left half of cell but not at left edge of grid
                    if j > 1
                        vxm₁₃ += 0.5*((dxmj/dx-0.5)^2) * (
                            vx[i, j-1] - 2.0*vx[i, j] + vx[i, j+1])
                        vxm₂₄ += 0.5*((dxmj/dx-0.5)^2) * (
                            vx[i+1, j-1] - 2.0*vx[i+1, j] + vx[i+1, j+1])
                    end
                end
                # compute current RK step vx
                vxrk4 = add_vrk4(vxrk4, vxm₁₃*(1.0-dymi/dy) + vxm₂₄*dymi/dy, rk)
                # interpolate vy
                i, j, dxmj, dymi = fix_distances(
                    xmrk4,
                    ymrk4,
                    xvy,
                    yvy,
                    dx,
                    dy,
                    jmin_vy,
                    jmax_vy,
                    imin_vy,
                    imax_vy
                )
                # compute vy velocity for top and bottom of current cell
                vym₁₂ = vy[i, j]*(1.0-dymi/dy) + vy[i+1, j]*dymi/dy
                vym₃₄ = vy[i, j+1]*(1.0-dymi/dy) + vy[i+1, j+1]*dymi/dy
                # compute second order vy velocity corrections
                if dymi/dy >= 0.5
                    # in bottom half of cell but not at bottom edge of grid
                    if i < Ny-1
                        vym₁₂ += 0.5*((dymi/dy-0.5)^2) * (
                            vy[i, j] - 2.0*vy[i+1, j] + vy[i+2, j])
                        vym₃₄ += 0.5*((dymi/dy-0.5)^2) * (
                            vy[i, j+1] - 2.0*vy[i+1, j+1] + vy[i+2, j+1])
                    end      
                else
                    # in top half of cell but not at top edge of grid
                    if i > 1
                        vym₁₂ += 0.5*((dymi/dy-0.5)^2) * (
                            vy[i-1, j] - 2.0*vy[i, j] + vy[i+1, j])
                        vym₃₄ += 0.5*((dymi/dy-0.5)^2) * (
                            vy[i-1, j+1] - 2.0*vy[i, j+1] + vy[i+1, j+1])
                    end
                end
                # compute current RK step vy
                vyrk4 = add_vrk4(vyrk4, vym₁₂*(1.0-dxmj/dx) + vym₃₄*dxmj/dx, rk)
                # calculate next RK step x and y positions if not at final step
                if rk < 4
                    xmrk4 = xmm + dt*crk4[rk]*vxrk4[rk]
                    ymrk4 = ymm + dt*crk4[rk]*vyrk4[rk]
                end
            end # RK4 solid velocity loop
            # advance marker using RK4 solid velocity
            xm[m] += dt * dot4(brk4, vxrk4)
            ym[m] += dt * dot4(brk4, vyrk4)
            # reset RK4 scheme for fluid velocity backtracing
            xmm = xmrk4 = xm[m]
            ymm = ymrk4 = ym[m]      
            vxrk4 = @SVector zeros(4)
            vyrk4 = @SVector zeros(4)
            # backtrack marker using RK4 scheme on fluid velocity
            for rk=1:1:4
                # interpolate vxf
                i, j, dxmj, dymi = fix_distances(
                    xmrk4,
                    ymrk4,
                    xvx,
                    yvx,
                    dx,
                    dy,
                    jmin_vx,
                    jmax_vx,
                    imin_vx,
                    imax_vx
                )
                # compute vxf velocity for left and right of current cell
                vxfm₁₃ = vxf[i, j]*(1.0-dxmj/dx) + vxf[i, j+1]*dxmj/dx
                vxfm₂₄ = vxf[i+1, j]*(1.0-dxmj/dx) + vxf[i+1, j+1]*dxmj/dx
                # compute second order vxf velocity corrections
                if dxmj/dx >= 0.5 
                    # in right half of cell but not at right edge of grid
                    if j < Nx-1
                        vxfm₁₃ += 0.5*((dxmj/dx-0.5)^2) * (
                            vxf[i, j] - 2.0*vxf[i, j+1] + vxf[i, j+2])
                        vxfm₂₄ += 0.5*((dxmj/dx-0.5)^2) * (
                            vxf[i+1, j] - 2.0*vxf[i+1, j+1] + vxf[i+1, j+2])
                    end
                else
                    # in left half of cell but not at left edge of grid
                    if j > 1
                        vxfm₁₃ += 0.5*((dxmj/dx-0.5)^2) * (
                            vxf[i, j-1] - 2.0*vx[i, j] + vxf[i, j+1])
                        vxfm₂₄ += 0.5*((dxmj/dx-0.5)^2) * (
                            vxf[i+1, j-1] - 2.0*vxf[i+1, j] + vxf[i+1, j+1])
                    end
                end
                # compute current RK step vxf
                vxrk4 = add_vrk4(
                    vxrk4, vxfm₁₃*(1.0-dymi/dy) + vxfm₂₄*dymi/dy, rk)
                # interpolate vyf
                i, j, dxmj, dymi = fix_distances(
                    xmrk4,
                    ymrk4,
                    xvy,
                    yvy,
                    dx,
                    dy,
                    jmin_vy,
                    jmax_vy,
                    imin_vy,
                    imax_vy
                )
                # compute vyf velocity for top and bottom of current cell
                vyfm₁₂ = vyf[i, j]*(1.0-dymi/dy) + vyf[i+1, j]*dymi/dy
                vyfm₃₄ = vyf[i, j+1]*(1.0-dymi/dy) + vyf[i+1, j+1]*dymi/dy
                # compute second order vyf velocity corrections
                if dymi/dy >= 0.5
                    # in bottom half of cell but not at bottom edge of grid
                    if i < Ny-1
                        vyfm₁₂ += 0.5*((dymi/dy-0.5)^2) * (
                            vyf[i, j] - 2.0*vyf[i+1, j] + vyf[i+2, j])
                        vyfm₃₄ += 0.5*((dymi/dy-0.5)^2) * (
                            vyf[i, j+1] - 2.0*vyf[i+1, j+1] + vyf[i+2, j+1])
                    end
                else
                    # in top half of cell but not at top edge of grid
                    if i > 1
                        vyfm₁₂ += 0.5*((dymi/dy-0.5)^2) * (
                            vyf[i-1, j] - 2.0*vyf[i, j] + vyf[i+1, j])
                        vyfm₃₄ += 0.5*((dymi/dy-0.5)^2) * (
                            vyf[i-1, j+1] - 2.0*vyf[i, j+1] + vyf[i+1, j+1])
                    end
                end
                # compute current RK step vyf
                vyrk4 = add_vrk4(
                    vyrk4, vyfm₁₂*(1.0-dxmj/dx) + vyfm₃₄*dxmj/dx, rk)
                # calculate next RK step x and y positions if not at final step
                if rk < 4
                    xmrk4 = xmm - dt*crk4[rk]*vxrk4[rk]
                    ymrk4 = ymm - dt*crk4[rk]*vyrk4[rk]
                end
            end # RK4 fluid velocity loop
            # backtrace marker using RK4 fluid velocity
            xmrk4 = xmm - dt*dot4(brk4, vxrk4)
            ymrk4 = ymm - dt*dot4(brk4, vyrk4)
            # interpolate fluid temperature at backtraced marker position
            i, j, weights = fix_weights(
                xmrk4,
                ymrk4,
                xp,
                yp,
                dx,
                dy,
                jmin_p,
                jmax_p,
                imin_p,
                imax_p,
            )
            # interpolate backtraced local fluid temperature from P grid
            tkfm₀ = dot4(grid_vector(i, j, tk2), weights)
            # compute marker fluid-solid temperature difference
            δtkfsm = tkfm₀ - tksm₀
            # correct marker temperature
            
            tkm[m] = (
                (1.0-phim[m])*tkm[m]*rhocpsolidm[tm[m]]
                    + phim[m]*(tkm[m]+δtkfsm)*compute_rhocpfluidm(tkm[m], mode)
            ) / (
                (1.0-phim[m])*rhocpsolidm[tm[m]]
                + phim[m]*compute_rhocpfluidm(tkm[m], mode)
            )
        end # marker loop
    end # @inbounds   
end # timeit to "move_markers_rk4!"
    return nothing
end # function move_markers_rk4!

"""
Backtrack pressure nodes using classic Runge-Kutta integration (RK4) 
to update total, solid, and fluid pressure under consideration of 
two-phase flow velocities.

$(SIGNATURES)

# Details

    - pr: total pressure at P nodes
    - pr0: previous time step total pressure at P nodes
    - ps: solid pressure at P nodes
    - ps0: previous time step solid pressure at P nodes
    - pf: fluid pressure at P nodes
    - pf0: previous time step fluid pressure at P nodes
    - vx: solid vx-velocity at Vx nodes
    - vy: solid vy-velocity at Vy nodes 
    - vxf: fluid vx-velocity at Vx nodes 
    - vyf: fluid vy-velocity at Vy nodes
    - dt: computational time step

# Returns

    - nothing
"""
function backtrace_pressures_rk4!(
    pr, pr0, ps, ps0, pf, pf0, vx, vy, vxf, vyf, dt)
@timeit to "backtrace_pressures_rk4!" begin
    @inbounds begin
        # setup RK4 scheme
        vxm = zeros(4)
        vym = zeros(4)
        # advance pressure generation
        pr0 .= pr
        ps0 .= ps
        # backtrace P nodes: total and solid pressure
        for jj=2:1:Nx, ii=2:1:Ny
            xA = xcur = xp[jj]
            yA = ycur =yp[ii]
            for rk=1:1:4
                # interpolate vx
                i, j, dxmj, dymi = fix_distances(
                    xcur,
                    ycur,
                    xvx,
                    yvx,
                    dx,
                    dy,
                    jmin_vx,
                    jmax_vx,
                    imin_vx,
                    imax_vx
                )
                # compute vx velocity for left and right of current cell
                vxm₁₃ = vx[i, j]*(1.0-dxmj/dx) + vx[i, j+1]*dxmj/dx
                vxm₂₄ = vx[i+1, j]*(1.0-dxmj/dx) + vx[i+1, j+1]*dxmj/dx
                # compute second order vx velocity corrections
                if dxmj/dx >= 0.5
                    if j < Nx-1
                        vxm₁₃ += 0.5*((dxmj/dx-0.5)^2) * (
                            vx[i, j] - 2.0*vx[i, j+1] + vx[i, j+2])
                        vxm₂₄ += 0.5*((dxmj/dx-0.5)^2) * (
                            vx[i+1, j] - 2.0*vx[i+1, j+1] + vx[i+1, j+2])
                    end
                else
                    if j > 1
                        vxm₁₃ += 0.5*((dxmj/dx-0.5)^2) * (
                            vx[i, j-1] - 2.0*vx[i, j] + vx[i, j+1])
                        vxm₂₄ += 0.5*((dxmj/dx-0.5)^2) * (
                            vx[i+1, j-1] - 2.0*vx[i+1, j] + vx[i+1, j+1])
                    end
                end
                # compute current RK step vx
                vxm[rk] = (1.0-dymi/dy)*vxm₁₃ + (dymi/dy)*vxm₂₄
                # interpolate vy
                i, j, dxmj, dymi = fix_distances(
                    xcur,
                    ycur,
                    xvy,
                    yvy,
                    dx,
                    dy,
                    jmin_vy,
                    jmax_vy,
                    imin_vy,
                    imax_vy
                )
                # compute vy velocity for left and right of current cell
                vym₁₂ = vy[i, j]*(1.0-dymi/dy) + vy[i+1, j]*dymi/dy
                vym₃₄ = vy[i, j+1]*(1.0-dymi/dy) + vy[i+1, j+1]*dymi/dy
                # compute second order vy velocity corrections
                if dymi/dy >= 0.5
                    if i < Ny-1
                        vym₁₂ += 0.5*((dymi/dy-0.5)^2) * (
                            vy[i, j] - 2.0*vy[i+1, j] + vy[i+2, j])
                        vym₃₄ += 0.5*((dymi/dy-0.5)^2) * (
                            vy[i, j+1] - 2.0*vy[i+1, j+1] + vy[i+2, j+1])
                    end      
                else
                    if i > 1
                        vym₁₂=vym₁₂+1/2*((dymi/dy-0.5)^2)*(vy[i-1,j]-2*vy[i,j]+vy[i+1,j])
                        vym₃₄=vym₃₄+1/2*((dymi/dy-0.5)^2)*(vy[i-1,j+1]-2*vy[i,j+1]+vy[i+1,j+1])
                    end
                end
                # compute current RK step vy
                vym[rk] = (1.0-dxmj/dx)*vym₁₂ + (dxmj/dx)*vym₃₄
                # calculate next RK step x and y positions if not at final 
                if rk==1 || rk==2
                    xcur=xA - 0.5*dt*vxm[rk]
                    ycur=yA - 0.5*dt*vym[rk]
                elseif rk==3
                    xcur = xA - dt*vxm[rk]
                    ycur = yA - dt*vym[rk]
                end
            end
            # compute effective velocity using RK4
            xcur = xA - dt*1//6*(vxm[1]+2*vxm[2]+2*vxm[3]+vxm[4])
            ycur = yA - dt*1//6*(vym[1]+2*vym[2]+2*vym[3]+vym[4])
            i, j, weights = fix_weights(
                xcur,
                ycur,
                xp,
                yp,
                dx,
                dy,
                jmin_p,
                jmax_p,
                imin_p,
                imax_p
            )
            pr0[ii,jj] = dot4(grid_vector(i, j, pr), weights)
            ps0[ii,jj] = dot4(grid_vector(i, j, ps), weights)
        end # jj, ii total and solid pressure loop
        # backtrace P nodes: fluid pressure
        pf0 .= pf
        for jj=2:1:Nx, ii=2:1:Ny
            # Save initial nodal coordinates
            xA = xcur = xp[jj]
            yA = ycur = yp[ii]
            for rk=1:1:4
                # interpolate vx
                i, j, dxmj, dymi = fix_distances(
                    xcur,
                    ycur,
                    xvx,
                    yvx,
                    dx,
                    dy,
                    jmin_vx,
                    jmax_vx,
                    imin_vx,
                    imax_vx
                )
                # compute vx velocity for left and right of current cell
                vxm₁₃ = vxf[i, j]*(1.0-dxmj/dx) + vxf[i, j+1]*dxmj/dx
                vxm₂₄ = vxf[i+1, j]*(1.0-dxmj/dx) + vxf[i+1, j+1]*dxmj/dx
                # compute second order vx velocity corrections
                if dxmj/dx>=0.5
                    if j < Nx-1
                        vxm₁₃ += 0.5*((dxmj/dx-0.5)^2) * (
                            vxf[i, j] - 2.0*vxf[i, j+1] + vxf[i, j+2])
                        vxm₂₄ += 0.5*((dxmj/dx-0.5)^2) * (
                            vxf[i+1, j] - 2.0*vxf[i+1, j+1] + vxf[i+1, j+2])
                    end
                else
                    if j > 1
                        vxm₁₃ += 0.5*((dxmj/dx-0.5)^2) * (
                            vxf[i, j-1] - 2.0*vxf[i, j] + vxf[i, j+1])
                        vxm₂₄ += 0.5*((dxmj/dx-0.5)^2) * (
                            vxf[i+1, j-1] - 2.0*vxf[i+1, j] + vxf[i+1, j+1])
                    end
                end
                # compute current RK step vx
                vxm[rk] = (1.0-dymi/dy)*vxm₁₃ + (dymi/dy)*vxm₂₄
                # interpolate vy
                i, j, dxmj, dymi = fix_distances(
                    xcur,
                    ycur,
                    xvy,
                    yvy,
                    dx,
                    dy,
                    jmin_vy,
                    jmax_vy,
                    imin_vy,
                    imax_vy
                )
                # compute vy velocity for top and bottom of current cell
                vym₁₂ = vyf[i, j]*(1.0-dymi/dy) + vyf[i+1, j]*dymi/dy
                vym₃₄ = vyf[i, j+1]*(1.0-dymi/dy) + vyf[i+1, j+1]*dymi/dy
                # compute second order vy velocity corrections
                if dymi/dy >= 0.5
                    if i < Ny-1
                        vym₁₂ += 0.5*((dymi/dy-0.5)^2) * (
                            vyf[i, j] - 2.0*vyf[i+1, j] + vyf[i+2, j])
                        vym₃₄ += 0.5*((dymi/dy-0.5)^2) * (
                            vyf[i, j+1] - 2.0*vyf[i+1, j+1] + vyf[i+2, j+1])
                    end      
                else
                    if i > 1
                        vym₁₂ += 0.5*((dymi/dy-0.5)^2) * (
                            vyf[i-1, j] - 2.0*vyf[i, j] + vyf[i+1, j])
                        vym₃₄ += 0.5*((dymi/dy-0.5)^2) * (
                            vyf[i-1, j+1] - 2.0*vyf[i, j+1] + vyf[i+1, j+1])
                    end
                end
                # compute current RK step vy
                vym[rk] = (1.0-dxmj/dx)*vym₁₂ + (dxmj/dx)*vym₃₄
                # calculate next RK step x and y positions if not at final 
                if rk==1 || rk==2
                    xcur = xA - 0.5*dt*vxm[rk]
                    ycur = yA - 0.5*dt*vym[rk]
                elseif rk==3
                    xcur = xA - dt*vxm[rk]
                    ycur = yA - dt*vym[rk]
                end
            end
            # compute effective velocity using RK4
            xcur = xA - dt*1//6*(vxm[1]+2*vxm[2]+2*vxm[3]+vxm[4])
            ycur = yA - dt*1//6*(vym[1]+2*vym[2]+2*vym[3]+vym[4])
            i, j, weights = fix_weights(
                xcur,
                ycur,
                xp,
                yp,
                dx,
                dy,
                jmin_p,
                jmax_p,
                imin_p,
                imax_p
            )
            pf0[ii,jj] = dot4(grid_vector(i, j, pf), weights)
        end
    end # inbounds
end # @timeit to "backtrace_pressures_rk4!"
    return nothing
end # function backtrace_pressures_rk4!

"""
Update marker population geometry status given a marker number and 
nearest top/left marker grid point.

$(SIGNATURES)

# Details

    - m: marker number
    - i: vertical index of top/left marker grid point
    - j: horizontal index of top/left marker grid point
    - xm: horizontal x-position of markers
	- ym: vertical y-position of markers
    - mdis: minimum distance of marker launch anchor points to nearest marker
    - mnum: number of marker nearest to marker launch anchor positions

# Returns

    - nothing
"""
function update_marker_population_geometry!(m, i, j, xm, ym, mdis, mnum)
@timeit to "update_marker_population_geometry!" begin
    @inbounds begin
        dismij = distance(xm[m], ym[m], xxm[j], yym[i])
        dismi1j = distance(xm[m], ym[m], xxm[j], yym[i+1])
        dismij1 = distance(xm[m], ym[m], xxm[j+1], yym[i])
        dismi1j1 = distance(xm[m], ym[m], xxm[j+1], yym[i+1])
        if dismij < mdis[i, j]
            mdis[i, j] = dismij
            mnum[i, j] = m
        end
        if dismi1j < mdis[i+1, j]
            mdis[i+1, j] = dismi1j
            mnum[i+1, j] = m
        end
        if dismij1 < mdis[i, j+1]
            mdis[i, j+1] = dismij1
            mnum[i, j+1] = m
        end
        if dismi1j1 < mdis[i+1, j+1]
            mdis[i+1, j+1] = dismi1j1
            mnum[i+1, j+1] = m
        end
    end # @inbounds
end # @timeit to "update_marker_population_geometry!"
    return nothing
end

"""
Add markers to populate currently sparsely filled grid areas.

$(SIGNATURES)

# Details

    - xm: horizontal x-position of markers
	- ym: vertical y-position of markers
	- tm: type of markers
	- tkm: temperature of markers 
	- phim: porosity of markers 
	- sxxm: marker σ′xx of markers
	- sxym: σxy of markers
	- etavpm: viscoplastic viscosity of markers
    - phinewm: reacted porosity of markers
    - pfm0: previous fluid pressure of markers
    - XWsolidm: melt molar fraction of markers
    - XWsolidm0: previous melt molar fraction of markers
    - rhototalm: total density of markers
    - rhocptotalm : total volumetric heat capacity of markers
    - etatotalm: total viscosity of markers
    - hrtotalm: total radiogenic heat production of markers
    - ktotalm: total thermal conductivity of markers
    - inv_gggtotalm: inverse of total shear modulus of markers
    - fricttotalm: total friction coefficient of markers
    - cohestotalm: total compressive strength of markers
    - tenstotalm: total tensile strength of markers
    - rhofluidcur: fluid density of markers
    - alphasolidcur: solid thermal expansion coefficient of markers
    - alphafluidcur: fluid thermal expansion coefficient of markers
    - tkm_rhocptotalm: total thermal energy of markers
    - etafluidcur_inv_kphim: fluid viscosity over permeability of markers
	- mdis: minimum distance of marker launch anchor points to nearest marker 
	- mnum: number of marker nearest to marker launch anchor positions

# Returns

    - marknum: updated number of markers in use
"""
function replenish_markers!(
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
@timeit to "replenish_markers!" begin
    # reset marker population geometry tracker
    mdis .= mdis_init
    mnum .= 0
    # establish marker distribution
    @inbounds begin
        for m=1:1:length(xm)
            i, j = fix(
                xm[m],
                ym[m],
                xxm,
                yym,
                dxm,
                dym,
                jmin_m,
                jmax_m,
                imin_m,
                imax_m
            )
            update_marker_population_geometry!(m, i, j, xm, ym, mdis, mnum)
        end 
        dii = 5 * Nymc
        djj = 5 * Nxmc
        for j=1:1:Nxm, i=1:1:Nym
            if mnum[i, j] == 0
                for jj=max(j-djj, 1):1:min(j+djj, Nxm)
                    for ii=max(i-dii, 1):1:min(i+dii, Nym)
                        if mnum[ii, jj] > 0
                            m = mnum[ii, jj]
                            dismij = distance(xm[m], ym[m], xxm[j], yym[i])
                            if dismij < mdis[i, j]
                                mdis[i, j] = dismij
                                mnum[i, j] = -m
                            end
                        end
                    end
                end 
                # add new marker            
                if mnum[i, j] < 0
                    # add marker
                    if randomized 
                        # production runs
                        push!(xm, xxm[j] + (rand(rgen)-0.5)*dxm)
                        push!(ym, yym[i] + (rand(rgen)-0.5)*dym)
                    else
                        # for testing
                        push!(xm, xxm[j])
                        push!(ym, yym[i])
                    end
                    # copy marker properties
                    m = -mnum[i,j]
                    push!(tm, tm[m])
                    push!(tkm, tkm[m])
                    push!(phim, phim[m])
                    push!(sxxm, sxxm[m])
                    push!(sxym, sxym[m])
                    push!(etavpm, etavpm[m])
                    push!(phinewm, phinewm[m])
                    push!(pfm0, pfm0[m])
                    push!(XWsolidm, XWsolidm[m])
                    push!(XWsolidm0, XWsolidm0[m])
                    push!(rhototalm, rhototalm[m])
                    push!(rhocptotalm, rhocptotalm[m])
                    push!(etatotalm, etatotalm[m])
                    push!(hrtotalm, hrtotalm[m])
                    push!(ktotalm, ktotalm[m])
                    push!(inv_gggtotalm, inv_gggtotalm[m])
                    push!(fricttotalm, fricttotalm[m])
                    push!(cohestotalm, cohestotalm[m])
                    push!(tenstotalm, tenstotalm[m])
                    push!(rhofluidcur, rhofluidcur[m])
                    push!(alphasolidcur, alphasolidcur[m])
                    push!(alphafluidcur, alphafluidcur[m])
                    push!(tkm_rhocptotalm, tkm_rhocptotalm[m])
                    push!(etafluidcur_inv_kphim, etafluidcur_inv_kphim[m])
                end
            end
        end  
    end # @inbounds  
    return length(xm)
end # @timeit to "replenish_markers!"
end # function replenish_markers!

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
@timeit to "save_state" begin
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
end # @timeit to "save_state"
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
@timeit to "simulation_loop setup" begin
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
    # thermal solver
    RT, ST = setup_thermal_lse()
    # gravitational solver
    RP, SP= setup_gravitational_lse()
    # Pardiso MKL solver
    if use_pardiso
        pardiso_solver = MKLPardisoSolver()
        initialize_pardiso!(pardiso_solver, iparms)
    # else
    #     F = lu(fdrand(Nx1*Ny1*6, 1, 1, matrixtype=ExtendableSparseMatrix))
    end

end # @timeit to "simulation_loop setup"

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
@timeit to "set up interpolation arrays" begin
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
end # @timeit to "set up interpolation arrays" 

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
    @timeit to "solve gravitational LSE" begin
        SP = LP \ RP
    end # @timeit to "solve gravitational LSE"
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
    @timeit to "thermochemical iteration (outer)" begin
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
    @timeit to "save initial viscosity, yielding nodes" begin
            ETA00 .= ETA
            YNY00 .= YNY
            if timestep == 1
                # no elastic compaction during first timestep
                BETAPHI .= 0.0
            end
    end # @timeit to "save initial viscosity, yielding nodes"

    @timeit to "advance pressure generation" begin
            # advance pressure generation inside thermochemical iteration
            pr0 .= pr
            pf0 .= pf
    end # @timeit to "advance pressure generation"

            # perform plastic iterations
            for iplast=1:1:titermax
    @timeit to "plastic iteration (inner)" begin
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
                    R
                )
                # solve hydromechanical system of equations
                @info "starting hydro-mechanical solver $titer-$iplast"
    @timeit to "solve hydromechanical system" begin
                if use_pardiso
                    # S = solve(pardiso_solver, L.cscmatrix, R)
                    set_phase!(
                        pardiso_solver, Pardiso.ANALYSIS_NUM_FACT_SOLVE_REFINE)
                    pardiso(
                        pardiso_solver,
                        S,
                        get_matrix(pardiso_solver, L.cscmatrix, :N),
                        R
                    )
                    set_phase!(pardiso_solver, Pardiso.RELEASE_ALL)
                    pardiso(pardiso_solver, S, L.cscmatrix, R)
                else
                    # lu!(F, L)
                    # S = F \ R
                    S = L \ R
                    # prob = LinearProblem(L, R)
                    # S = solve(prob)
                end
    end # @timeit to "solve hydromechanical system"
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
                    dt
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
                    dt
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
    end # timeit to "plastic iteration"
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
            @timeit to "solve thermal system" ST = LT \ RT
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
    end # @timeit to "thermochemical iteration (outer)"
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

end # module Erebus

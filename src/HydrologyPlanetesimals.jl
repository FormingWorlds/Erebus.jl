module HydrologyPlanetesimals

using Base.Threads
using DocStringExtensions
using ExtendableSparse
using MAT
using Parameters
using ProgressMeter
using SparseArrays
using StaticArrays
using TimerOutputs

export run_simulation

const to = TimerOutput()

"""
Static parameters: Grids, markers, switches, constants, etc. which remain
constant throughout the simulation.

$(TYPEDFIELDS)
"""
@with_kw struct StaticParameters
    # radioactive switches
    "radioactive heating from 26Al active"
    hr_al::Bool = true
    "radioactive heating from 60Fe active"	
    hr_fe::Bool = true
    # model size, geometry, and resolution
    "horizontal model size [m]"
    xsize::Float64 = 140_000.0
    "vertical model size [m]"
    ysize::Float64 = 140_000.0
    "horizontal center of model"
    xcenter::Float64 = xsize / 2
    "vertical center of model"
    ycenter::Float64 = ysize / 2  
    "basic grid resolution in x direction (horizontal)"
    Nx::Int64 = 141
    "basic grid resolution in y direction (vertical)"	
    Ny::Int64 = 141
    "Vx, Vy, P grid resolution in x direction (horizontal)"
    Nx1::Int64 = Nx + 1
    "Vx/Vy/P grid resolution in y direction (vertical)"
    Ny1::Int64 = Ny + 1
    "horizontal grid step [m]"
    dx::Float64 = xsize / (Nx-1)
    "vertical grid step [m]"
    dy::Float64 = ysize / (Ny-1)
    # basic grid min/max assignables indices
    "minimum assignable basic grid index in x direction"
    jmin_basic::Int64 = 1
    "minimum assignable basic grid index in y direction"
    imin_basic::Int64 = 1
    "maximum assignable basic grid index in x direction"
    jmax_basic::Int64 = Nx - 1
    "maximum assignable basic grid index in y direction"
    imax_basic::Int64 = Ny - 1
    # Vx grid min/max assignables indices
    "minimum assignable Vx grid index in x direction"
    jmin_vx::Int64 = 1
    "minimum assignable Vx grid index in y direction"
    imin_vx::Int64 = 1
    "maximum assignable Vx grid index in x direction"
    jmax_vx::Int64 = Nx - 1
    "maximum assignable Vx grid index in y direction"
    imax_vx::Int64 = Ny
    # Vy grid min/max assignables indices
    "minimum assignable Vy grid index in x direction"
    jmin_vy::Int64 = 1
    "minimum assignable Vy grid index in y direction"
    imin_vy::Int64 = 1
    "maximum assignable Vy grid index in x direction"
    jmax_vy::Int64 = Nx
    "maximum assignable Vy grid index in y direction"
    imax_vy::Int64 = Ny - 1
    # P grid min/max assignables indices
    "minimum assignable P grid index in x direction"
    jmin_p::Int64 = 1
    "minimum assignable P grid index in y direction"
    imin_p::Int64 = 1
    "maximum assignable P grid index in x direction"
    jmax_p::Int64 = Nx
    "maximum assignable P grid index in y direction"
    imax_p::Int64 = Ny
    # planetary parameters
    "planetary radius [m]"
    rplanet::Float64 = 50_000.0
    "crust radius [m]"
    rcrust::Float64 = 48_000.0
    "surface pressure [Pa]"
    psurface::Float64 = 1e+3
    # marker count and initial spacing
    "number of markers per cell in horizontal direction"
    Nxmc::Int64 = 4
    "number of markers per cell in vertical direction"
    Nymc::Int64 = 4
    "marker grid resolution in horizontal direction"
    Nxm::Int64 = (Nx - 1) * Nxmc
    "marker grid resolution in vertical direction"
    Nym::Int64 = (Ny - 1) * Nymc
    "marker grid step in horizontal direction"
    dxm::Float64 = xsize / Nxm
    "marker grid step in vertical direction"
    dym::Float64 = ysize / Nym
    "number of markers at start"
    start_marknum::Int64 = Nxm * Nym
    # physical constants
    "gravitational constant [m^3*kg^-1*s^-2]"
    G::Float64 = 6.672e-11
    "scaled pressure"    
    pscale::Float64 = 1e+23 / dx
    # materials properties:              planet      crust       space
    "solid Density [kg/m^3]"
    rhosolidm::SVector{3, Float64}      = [ 3300.0    , 3300.0    ,    1.0    ]
    "fluid density [kg/m^3]"	
    rhofluidm::SVector{3, Float64}      = [ 7000.0    , 7000.0    , 1000.0    ]
    "solid viscosity [Pa*s]"
    etasolidm::SVector{3, Float64}      = [    1.0e+16,    1.0e+16,    1.0e+14]
    "molten solid viscosity [Pa*s]"
    etasolidmm::SVector{3, Float64}     = [    1.0e+14,    1.0e+14,    1.0e+14]
    "fluid viscosity [Pa*s]"
    etafluidm::SVector{3, Float64}      = [    1.0e-02,    1.0e-02,    1.0e+12]
    "molten fluid viscosity [Pa*s]"
    etafluidmm::SVector{3, Float64}     = [    1.0e-02,    1.0e-02,    1.0e+12]
    "solid volumetric heat capacity [kg/m^3]"
    rhocpsolidm::SVector{3, Float64}    = [    3.3e+06,    3.3e+06,    3.0e+06]
    "fluid volumetric heat capacity [kg/m^3]"
    rhocpfluidm::SVector{3, Float64}    = [    7.0e+06,    7.0e+06,    3.0e+06]
    "solid thermal expansion [1/K]"
    alphasolidm::SVector{3, Float64}    = [    3.0e-05,    3.0e-05,    0.0    ]
    "fluid thermal expansion [1/K]"
    alphafluidm::SVector{3, Float64}    = [    5.0e-05,    5.0e-05,    0.0    ]
    "solid thermal conductivity [W/m/K]"
    ksolidm::SVector{3, Float64}        = [    3.0    ,    3.0    , 3000.0    ]
    "fluid thermal conductivity [W/m/K]"
    kfluidm::SVector{3, Float64}        = [   50.0    ,   50.0    , 3000.0    ]
    "solid radiogenic heat production [W/m^3]"
    start_hrsolidm::SVector{3, Float64} = [    0.0    ,    0.0    ,    0.0    ]
    "fluid radiogenic heat production [W/m^3]"
    start_hrfluidm::SVector{3, Float64} = [    0.0    ,    0.0    ,    0.0    ]
    "solid shear modulus [Pa]"
    gggsolidm::SVector{3, Float64}      = [    1.0e+10,    1.0e+10,    1.0e+10]
    "solid friction coefficient"
    frictsolidm::SVector{3, Float64}    = [    0.6    ,    0.6    ,    0.0    ]
    "solid compressive strength [Pa]"
    cohessolidm::SVector{3, Float64}    = [    1.0e+08,    1.0e+08,    1.0e+08]
    "solid tensile strength [Pa]"
    tenssolidm ::SVector{3, Float64}    = [    6.0e+07,    6.0e+07,    6.0e+07]
    "standard permeability [m^2]"
    kphim0::SVector{3, Float64}         = [    1.0e-13,    1.0e-13,    1.0e-17]
    "initial temperature [K]"
    tkm0::SVector{3, Float64}           = [  300.0    ,  300.0    ,  273.0    ]
    "Coefficient to compute compaction viscosity from shear viscosity"
    etaphikoef::Float64 = 1e-4
    # 26Al decay
    "26Al half life [s]"
    t_half_al::Float64 = 717000 * 31540000
    "26Al decay constant"
    tau_al::Float64 = t_half_al / log(2)
    "initial ratio of 26Al and 27Al isotopes"
    ratio_al::Float64 = 5.0e-5
    "E 26Al [J]"
    E_al::Float64 = 5.0470e-13
    "26Al atoms/kg"
    f_al::Float64 = 1.9e23
    # 60Fe decay
    "60Fe half life [s]"	
    t_half_fe::Float64 = 2620000 * 31540000
    "60Fe decay constant"
    tau_fe::Float64 = t_half_fe / log(2)
    "initial ratio of 60Fe and 56Fe isotopes"	
    ratio_fe::Float64 = 1e-6
    "E 60Fe [J]"	
    E_fe::Float64 = 4.34e-13
    "60Fe atoms/kg"	
    f_fe::Float64 = 1.957e24
    # melting temperatures
    "silicate melting temperature [K]"
    tmsilicate::Float64 = 1e+6
    "iron melting temperature [K]"
    tmiron::Float64 = 1273 
    # porosities
    "standard Fe fraction [porosity]"
    phim0::Float64 = 0.2
    "min porosity"	
    phimin::Float64 = 1e-4
    "max porosity"
    phimax::Float64 = 1 - phimin            
    # mechanical boundary conditions: free slip=-1 / no slip=1
    "mechanical boundary condition left"
    bcleft::Float64 = -1
    "mechanical boundary condition right"
    bcright::Float64 = -1
    "mechanical boundary condition top"
    bctop::Float64 = -1
    "mechanical boundary condition bottom"
    bcbottom::Float64 = -1
    # hydraulic boundary conditions: free slip=-1 / no slip=1
    "hydraulic boundary condition left"
    bcfleft::Float64 = -1
    "hydraulic boundary condition right"
    bcfright::Float64 = -1
    "hydraulic boundary condition top"
    bcftop::Float64 = -1
    "hydraulic boundary condition bottom"
    bcfbottom::Float64 = -1
    # extension/shortening velocities
    "shortening strain rate"
    strainrate::Float64 = 0e-13
    "x extension/shortening velocity left"
    vxleft::Float64 = strainrate * xsize / 2
    "x extension/shortening velocity right"
    vxright::Float64= -strainrate * xsize / 2
    "y extension/shortening velocity top"
    vytop::Float64 = - strainrate * ysize / 2
    "y extension/shortening velocity bottom"
    vybottom::Float64 = strainrate * ysize / 2
    # timestepping parameters
    "mat filename"
    nname::String = "madcph_"
    ".mat storage periodicity"
    savematstep::Int64 = 50
    "Maximal computational timestep [s]"
    dtelastic::Float64 = 1e+11 
    "Coefficient to decrease computational timestep"
    dtkoef::Float64 = 2 
    "Coefficient to increase computational timestep"
    dtkoefup::Float64 = 1.1 
    "Number of iterations before changing computational timestep"
    dtstep::Int64 = 200 
    "Max marker movement per time step [grid steps]"
    dxymax::Float64 = 0.05 
    "Weight of averaged velocity for moving markers"
    vpratio::Float64 = 1 / 3 
    "Max temperature change per time step [K]"
    DTmax::Float64 = 20 
    "Subgrid temperature diffusion parameter"
    dsubgridt::Float64 = 0 
    "Subgrid stress diffusion parameter"
    dsubgrids::Float64 = 0
    "length of year [s]"
    yearlength::Float64 = 365.25 * 24 * 3600
    "Time sum (start) [s]"
    start_time::Float64 = 1e6 * yearlength 
    "Time sum (end) [s]"
    endtime::Float64 = 15 * 1000000 * yearlength
    "Lower viscosity cut-off [Pa s]"	
    etamin::Float64 = 1e+12 
    "Upper viscosity cut-off [Pa s]"
    etamax::Float64 = 1e+23 
    "Number of plastic iterations"
    nplast::Int64 = 100000
    "Periodicity of visualization"
    visstep::Int64 = 1 
    "Tolerance level for yielding error()"
    yerrmax::Float64 = 1e+2 
    "Weight for old viscosity"
    etawt::Float64 = 0 
    "max porosity ratio change per time step"
    dphimax::Float64 = 0.01
    "starting timestep"
    start_step::Int64 = 1
    "number of timesteps to run"
    nsteps::Int64 = 30000 
end


"""

$(SIGNATURES)

# Details

    - sp: static simulation parameters

# Returns

    - timestep: simulation starting time step count
    - dt: simulation initial computational time step 
    - timesum: simulation starting time
    - marknum: initial number of markers
    - hrsolidm: initial radiogenic heat production solid phase
    - hrfluidm: initial radiogenic heat production fluid phase
    - YERRNOD: vector of summed yielding errors of nodes over plastic iterations
"""
function setup_dynamic_simulation_parameters(sp)
    @unpack start_step,
        dtelastic,
        start_time,
        start_marknum,
        start_hrsolidm,
        start_hrfluidm,
        nplast = sp
     # timestep counter (current), init to startstep
     timestep = start_step
     # computational timestep (current), init to dtelastic [s]
     dt = dtelastic
     # time sum (current), init to starttime [s]
     timesum = start_time
     # current number of markers, init to startmarknum
     marknum = start_marknum
     # radiogenic heat production solid phase
     hrsolidm = start_hrsolidm
     # radiogenic heat production fluid phase
     hrfluidm = start_hrfluidm
     # nodes yielding error vector of plastic iterations
     YERRNOD = zeros(Float64, nplast) 
    return timestep, dt, timesum, marknum, hrsolidm, hrfluidm, YERRNOD
end # function setup_dynamic_simulation_parameters()


"""
Set up staggered grid geometry with basic, Vx, Vy, and P nodes.

$(SIGNATURES)

# Details

    - sp: static simulation parameters

# Returns

    - x:
    - y:
    - xvx:
    - yvx:
    - xvy:
    - yvy:
    - xp:
    - yp:
"""
function setup_staggered_grid_geometry(sp)
    @unpack Nx, Ny, Nx1, Ny1, dx, dy, xsize, ysize = sp
    # basic nodes
    # x: horizontal coordinates of basic grid points [m]
    x = SVector{Nx, Float64}([j for j = 0:dx:xsize])
    # y: vertical coordinates of basic grid points [m]
    y = SVector{Ny, Float64}([j for j = 0:dy:ysize])
    # Vx nodes
    # xvx: horizontal coordinates of vx grid points [m]
    xvx = SVector{Ny1, Float64}([j for j = 0:dx:xsize+dy])
    # yvx: vertical coordinates of vx grid points [m]
    yvx = SVector{Nx1, Float64}([i for i = -dy/2:dy:ysize+dy/2])
    # Vy nodes
    # xvy: horizontal coordinates of vy grid points [m]
    xvy = SVector{Nx1, Float64}([j for j = -dx/2:dx:xsize+dx/2])
    # yvy: vertical coordinates of vy grid points [m]
    yvy = SVector{Ny1, Float64}([i for i = 0:dy:ysize+dy])
    # P nodes
    # xp: horizontal coordinates of p grid points [m]
    xp = SVector{Nx1, Float64}([j for j = -dx/2:dx:xsize+dx/2])
    # yp: vertical coordinates of p grid points [m]
    yp = SVector{Ny1, Float64}([i for i = -dy/2:dy:ysize+dy/2])
    return x, y, xvx, yvx, xvy, yvy, xp, yp
end # function setup_staggered_grid()


"""
Set up staggered grid properties for basic, Vx, Vy, and P nodes.

$(SIGNATURES)

# Details

    - sp: static simulation parameters
    - randomized: fill in random values for grid properties instead of zeros

# Returns

    - 
"""
function setup_staggered_grid_properties(sp; randomized=false)
    @unpack Nx, Ny, Nx1, Ny1 = sp
    # basic nodes
    # viscoplastic viscosity [Pa*s]
    ETA = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
    # viscous viscosity [Pa*s]
    ETA0 = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
    # shear modulus [Pa]
    GGG = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
    # epsilonxy [1/s]
    EXY = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
    # σxy [1/s]
    SXY = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
    # σ₀xy [1/s]
    SXY0 = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
    # rotation rate [1/s]
    wyx = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
    # compressive strength [Pa]
    COH = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
    # tensile strength [Pa]
    TEN = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
    # friction
    FRI = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
    # plastic yielding node property
    YNY = randomized ? rand(Bool, Ny, Nx) : zeros(Bool, Ny, Nx)
    # Vx nodes
    # density [kg/m^3]
    RHOX = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # fluid density [kg/m^3]
    RHOFX = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # thermal conductivity [W/m/K]
    KX = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # porosity
    PHIX = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # solid vx-velocity [m/s]
    vx = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # fluid vx-velocity [m/s]
    vxf = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # etafluid/kphi ratio [m^2]
    RX = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # qx-darcy flux [m/s]
    qxD = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # gx-gravity [m/s^2]
    gx = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # Vy nodes
    # density [kg/m^3]
    RHOY = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # fluid density [kg/m^3]
    RHOFY = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # thermal conductivity [W/m/K]
    KY = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # porosity
    PHIY = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # solid vy-velocity [m/s]
    vy = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # fluid vy-velocity [m/s]
    vyf = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # etafluid/kphi ratio [m^2]
    RY = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # qy-Darcy flux [m/s]
    qyD = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # gy-gravity [m/s^2]
    gy = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # P nodes
    # density [kg/m^3]
    RHO = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # volumetric heat capacity [J/m^3/K]
    RHOCP = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # thermal expansion [J/m^3/K]
    ALPHA = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # fluid thermal expansion [J/m^3/K]
    ALPHAF = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # radioactive heating [W/m^3]
    HR = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # adiabatic heating [W/m^3]
    HA = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # shear heating [W/m^3]
    HS = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # viscosity [Pa*s]
    ETAP = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # shear modulus [Pa]
    GGGP = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # EPSILONxx [1/s]
    EXX = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # σ′xx [1/s]
    SXX = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # σ₀′ (SIGMA0'xx) [1/s]
    SXX0 = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # previous temperature [K]
    tk1 = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # next temperature [K]
    tk2 = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # solid vx in pressure nodes [m/s]
    vxp = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # solid vy in pressure nodes [m/s]
    vyp = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # fluid vx in pressure nodes [m/s]
    vxpf = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # fluid vy in pressure nodes [m/s]
    vypf = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # total pressure [Pa]
    pr = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # fluid pressure [Pa]
    pf = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # solid pressure [Pa]
    ps = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # previous total pressure [Pa]
    pr0 = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # previous fluid pressure [Pa]
    pf0 = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # previous solid pressure [Pa]
    ps0 = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # bulk viscosity [Pa*s]
    ETAPHI = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # bulk compresibility [Pa*s]
    BETTAPHI = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # porosity
    PHI = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # Dln[(1-ϕ)/ϕ]/Dt
    APHI = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # gravity potential [J/kg]
    FI = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
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
        BETTAPHI,
        PHI,
        APHI,
        FI
    )
end # function setup_staggered_grid_properties()


"""
Define initial set of markers according to model parameters

$(SIGNATURES)

# Details

    - xm: array of x coordinates of markers
    - ym: array of y coordinates of markers
    - tm: array of material type of markers
    - phim: array of porosity of markers
    - etavpm: array of matrix viscosity of markers
    - rhototalm: array of total density of markers
    - rhocptotalm: array of total volumetric heat capacity of markers
    - etatotalm: array of total viscosity of markers
    - hrtotalm: array of total radiogenic heat production of markers
    - ktotalm: array of total thermal conductivity of markers
    - etafluidcur: array of fluid viscosity of markers
    - tkm: array of temperature of markers 
    - inv_gggtotalm: array of inverse of total shear modulus of markers
    - fricttotalm: array of total friction coefficient of markers
    - cohestotalm: array of total compressive strength of markers
    - tenstotalm: array of total tensile strength of markers
    - rhofluidcur: array of fluid density of markers
    - alphasolidcur: array of solid thermal expansion coefficient of markers
    - alphafluidcur: array of fluid thermal expansion coefficient of markers
    - sp: static simulation parameters

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
    etafluidcur,
    tkm,
    inv_gggtotalm,
    fricttotalm,
    cohestotalm,
    tenstotalm,
    rhofluidcur,
    alphasolidcur,
    alphafluidcur,
    sp
)
    @unpack xsize,
    ysize,
    xcenter,
    ycenter,
    Nxm,
    Nym,
    dxm,
    dym,
    rplanet,
    rcrust,
    phim0,
    phimin,
    etasolidm,
    rhosolidm,
    rhocpsolidm,
    etasolidm,
    etafluidm,
    tkm0,
    gggsolidm,
    frictsolidm,
    cohessolidm,
    tenssolidm,
    rhofluidm,
    start_hrsolidm,
    alphasolidm,
    alphafluidm,
    ksolidm = sp

    for jm=1:1:Nxm, im=1:1:Nym
        # calculate marker counter
        m = (jm-1) * Nym + im
        # define marker coordinates
        xm[m] = dxm/2 + (jm-1) * dxm + (rand()-0.5) * dxm
        ym[m] = dym/2 + (im-1) * dym + (rand()-0.5) * dym
        # primary marker properties 
        rmark = distance(xm[m], ym[m], xcenter, ycenter)
        if rmark < rplanet
            # planet
            tm[m] = ifelse(rmark>rcrust, 2, 1)
            # porosity
            phim[m] = phim0 * (1.0 + (rand()-0.5))
            # matrix viscosity
            etavpm[m] = etasolidm[tm[m]] # *exp(-28*phim[m])
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
            etafluidcur[m] = etafluidm[tm[m]]
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
    - tm: array of type of markers
    - tkm: array of temperature of markers
    - rhototalm: array of total density of markers
    - rhocptotalm: array of total volumetric heat capacity of markers
    - etasolidcur: array of solid viscosity of markers
    - etafluidcur: array of fluid viscosity of markers
    - etatotalm: array of total viscosity of markers
    - hrtotalm: array of total radiogenic heat production of markers
    - ktotalm: array of total thermal conductivity of markers
    - tkm_rhocptotalm: array of total thermal energy of markers
    - etafluidcur_inv_kphim: array of (fluid viscosity)/permeability of markers
    - phim: array of porosity of markers
    - hrsolidm: vector of radiogenic heat production of solid materials
    - hrfluidm: vector of radiogenic heat production of fluid materials
    - sp: static simulation parameters

# Returns

    - nothing
"""
function compute_marker_properties!(
    m,
    tm,
    tkm,
    rhototalm,
    rhocptotalm,
    etasolidcur,
    etafluidcur,
    etatotalm,
    hrtotalm,
    ktotalm,
    tkm_rhocptotalm,
    etafluidcur_inv_kphim,
    hrsolidm,
    hrfluidm,
    phim,
    sp::StaticParameters
)
    @unpack rhosolidm,
        rhofluidm,
        rhocpsolidm,
        rhocpfluidm,
        tmsilicate,
        tmiron,
        etamin,
        etasolidm,
        etasolidmm,
        etafluidm,
        etafluidmm,
        ksolidm,
        kfluidm,
        kphim0,
        phim0 = sp
# @timeit to "compute_marker_properties!" begin
    if tm[m] < 3
        # rocks
        rhototalm[m] = total(rhosolidm[tm[m]], rhofluidm[tm[m]], phim[m])
        rhocptotalm[m] = total(
            rhocpsolidm[tm[m]], rhocpfluidm[tm[m]], phim[m])
        etasolidcur[m] = ifelse(
            tkm[m]>tmsilicate, etasolidmm[tm[m]], etasolidm[tm[m]])
        etafluidcur[m] = ifelse(
            tkm[m]>tmiron, etafluidmm[tm[m]], etafluidm[tm[m]])
        etatotalm[m] = max(etamin, etasolidcur[m], etafluidcur[m])
        hrtotalm[m] = total(hrsolidm[tm[m]], hrfluidm[tm[m]], phim[m])
        ktotalm[m] = ktotal(ksolidm[tm[m]], kfluidm[tm[m]], phim[m])
    # else
        # air
        # pass  
    end
    # # common for rocks and air
    tkm_rhocptotalm[m] = tkm[m] * rhocptotalm[m]
    # kphim[m] = kphi(kphim0[tm[m]], phim0, phim[m])
    etafluidcur_inv_kphim[m] = etafluidcur[m]/kphi(
        kphim0[tm[m]], phim0, phim[m])
# end # @timeit to "compute_marker_properties!"
    return nothing
end # function compute_marker_properties!


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
    - phi: fraction of fluid

# Returns

    - total: computed total property
"""
function total(solid, fluid, phi)
    return solid * (1.0-phi) + fluid * phi
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
    return ((ksolid * kfluid/2 + ((ksolid * (3*phi-2)
                                 + kfluid * (1.0-3.0*phi))^2)/16)^0.5
            - (ksolid*(3.0*phi-2.0) + kfluid*(1.0-3.0*phi))/4)
end


"""
Compute iron porosity-dependent permeability.

$(SIGNATURES)

# Details

    - kphim0: standard permeability [m^2]
    - phim: actual (marker) porosity
    - phim0: standard iron fraction (porosity)

# Returns

    - kphim: iron porosity-dependent permeability [m^2]
"""
function kphi(kphim0, phim, phim0)
    return kphim0 * (phim/phim0)^3 / ((1.0-phim)/(1.0-phim0))^2
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
    return f * ratio * E * exp(-time/tau) / tau
end


"""
Compute radiogenic heat production of 26Al and 60Fe isotopes.

$(SIGNATURES)

# Details
    - timesum: time elapsed since initial conditions at start of simulation
    - sp: static simulation parameters

# Returns

    - hrsolidm: radiogenic heat production of 26Al [W/m^3]
    - hrfluidm: radiogenic heat production of 60Fe [W/m^3]
"""
function calculate_radioactive_heating(timesum, sp::StaticParameters)
    @unpack hr_al,
        f_al,
        ratio_al,
        E_al,
        tau_al,
        hr_fe,
        f_fe,
        ratio_fe,
        E_fe,
        tau_fe,
        rhosolidm,
        rhofluidm = sp
    #26Al: planet ✓, crust ✓, space ×
    if hr_al
        # 26Al radiogenic heat production [W/kg]
        Q_al = Q_radiogenic(f_al, ratio_al, E_al, tau_al, timesum)
        # Solid phase 26Al radiogenic heat production [W/m^3]
        hrsolidm = @SVector [Q_al*rhosolidm[1], Q_al*rhosolidm[2], 0.0]
    else
        hrsolidm = @SVector zeros(3)#[0.0, 0.0, 0.0]
    end    
    #60Fe: planet ✓, crust ×, space ×
    if hr_fe
        # 60Fe radiogenic heat production [W/kg]
        Q_fe = Q_radiogenic(f_fe, ratio_fe, E_fe, tau_fe, timesum)
        # Fluid phase 60Fe radiogenic heat production [W/m^3]
        hrfluidm = @SVector [Q_fe*rhofluidm[1], 0.0, 0.0]
    else
        hrfluidm = @SVector zeros(3)#[0.0, 0.0, 0.0]
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
# @timeit to "fix_weights" begin
    @inbounds begin
    j = min(max(trunc(Int, (x-x_axis[1])/dx)+1, jmin), jmax)
    i = min(max(trunc(Int, (y-y_axis[1])/dy)+1, imin), imax)
    dxmj = x - x_axis[j]
    dymi = y - y_axis[i]
    end # @inbounds
    weights = SVector(
        (1.0-dymi/dy) * (1.0-dxmj/dx),
        (dymi/dy) * (1.0-dxmj/dx),
        (1.0-dymi/dy) * (dxmj/dx),
        (dymi/dy) * (dxmj/dx)
    )
# end # @timeit to "fix_weights"
    return i, j, weights
end # function fix_weights


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
function interpolate_to_grid!(i, j, weights, property, grid)
# @timeit to "interpolate_to_grid!" begin
    grid[i, j, threadid()] += property * weights[1]
    grid[i+1, j, threadid()] += property * weights[2]
    grid[i, j+1, threadid()] += property * weights[3]
    grid[i+1, j+1, threadid()] += property * weights[4]
# end # @timeit to "interpolate_to_grid!"
    return nothing
end # function interpolate_to_grid!


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
@timeit to "interpolate_to_marker!()" begin
    marker_property[m] = (
        grid[i, j] * weights[1]
        + grid[i+1, j] * weights[2]
        + grid[i, j+1] * weights[3]
        + grid[i+1, j+1] * weights[4]
    )
end # @timeit to "interpolate_to_marker!()"
    return nothing
end # function interpolate_to_marker


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
# @timeit to "compute_basic_node_properties!" begin
# @timeit to "reduce" begin
    ETA0SUM = reduce(+, ETA0SUM, dims=3)[:, :, 1]
    ETASUM = reduce(+, ETASUM, dims=3)[:, :, 1]
    GGGSUM = reduce(+, GGGSUM, dims=3)[:, :, 1]
    SXYSUM = reduce(+, SXYSUM, dims=3)[:, :, 1]
    COHSUM = reduce(+, COHSUM, dims=3)[:, :, 1]
    TENSUM = reduce(+, TENSUM, dims=3)[:, :, 1]
    FRISUM = reduce(+, FRISUM, dims=3)[:, :, 1]
    WTSUM = reduce(+, WTSUM, dims=3)[:, :, 1]
# end # @timeit to "reduce"
# @timeit to "compute" begin
    ETA0[WTSUM.>0.0] .= ETA0SUM[WTSUM.>0.0] ./ WTSUM[WTSUM.>0.0]
    ETA[WTSUM.>0.0] .= ETASUM[WTSUM.>0.0] ./ WTSUM[WTSUM.>0.0]
    YNY[ETA.<ETA0] .= true
    GGG[WTSUM.>0.0] .= GGGSUM[WTSUM.>0.0] .\ WTSUM[WTSUM.>0.0]
    SXY0[WTSUM.>0.0] .= SXYSUM[WTSUM.>0.0] ./ WTSUM[WTSUM.>0.0]
    COH[WTSUM.>0.0] .= COHSUM[WTSUM.>0.0] ./ WTSUM[WTSUM.>0.0]
    TEN[WTSUM.>0.0] .= TENSUM[WTSUM.>0.0] ./ WTSUM[WTSUM.>0.0]
    FRI[WTSUM.>0.0] .= FRISUM[WTSUM.>0.0] ./ WTSUM[WTSUM.>0.0]
# end # @timeit to "compute"
# end # @timeit to "compute_basic_node_properties!"
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
# @timeit to "compute_vx_node_properties!" begin
# @timeit to "reduce" begin
    RHOXSUM = reduce(+, RHOXSUM, dims=3)[:, :, 1]
    RHOFXSUM = reduce(+, RHOFXSUM, dims=3)[:, :, 1]
    KXSUM = reduce(+, KXSUM, dims=3)[:, :, 1]
    PHIXSUM = reduce(+, PHIXSUM, dims=3)[:, :, 1]
    RXSUM = reduce(+, RXSUM, dims=3)[:, :, 1]
    WTXSUM = reduce(+, WTXSUM, dims=3)[:, :, 1]    
# end # @timeit to "reduce"
# @timeit to "compute" begin
    WTXSUM[WTXSUM .<= 0.0] .= Inf
    RHOX .= RHOXSUM ./ WTXSUM
    RHOFX .= RHOFXSUM ./ WTXSUM
    KX .= KXSUM ./ WTXSUM
    PHIX .= PHIXSUM ./ WTXSUM
    RX .= RXSUM ./ WTXSUM
# end # @timeit to "compute"
# end # @timeit to "compute_vx_node_properties!"
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
# @timeit to "compute_vy_node_properties!" begin
# @timeit to "reduce" begin
    RHOYSUM = reduce(+, RHOYSUM, dims=3)[:, :, 1]
    RHOFYSUM = reduce(+, RHOFYSUM, dims=3)[:, :, 1]
    KYSUM = reduce(+, KYSUM, dims=3)[:, :, 1]
    PHIYSUM = reduce(+, PHIYSUM, dims=3)[:, :, 1]
    RYSUM = reduce(+, RYSUM, dims=3)[:, :, 1]
    WTYSUM = reduce(+, WTYSUM, dims=3)[:, :, 1]    
# end # @timeit to "reduce"
# @timeit to "compute" begin
    WTYSUM[WTYSUM .<= 0.0] .= Inf
    RHOY .= RHOYSUM ./ WTYSUM
    RHOFY .= RHOFYSUM ./ WTYSUM
    KY .= KYSUM ./ WTYSUM
    PHIY .= PHIYSUM ./ WTYSUM
    RY .= RYSUM ./ WTYSUM
# end # @timeit to "compute"
# end # @timeit to "compute_vy_node_properties!"
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
    - BETTAPHI: BETTAPHI P node array
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
    BETTAPHI
)
# @timeit to "compute_p_node_properties!" begin
# @timeit to "reduce" begin
    RHOSUM = reduce(+, RHOSUM, dims=3)[:, :, 1]
    RHOCPSUM = reduce(+, RHOCPSUM, dims=3)[:, :, 1]
    ALPHASUM = reduce(+, ALPHASUM, dims=3)[:, :, 1]
    ALPHAFSUM = reduce(+, ALPHAFSUM, dims=3)[:, :, 1]
    HRSUM = reduce(+, HRSUM, dims=3)[:, :, 1]
    GGGPSUM = reduce(+, GGGPSUM, dims=3)[:, :, 1]
    SXXSUM = reduce(+, SXXSUM, dims=3)[:, :, 1]
    TKSUM = reduce(+, TKSUM, dims=3)[:, :, 1]
    PHISUM = reduce(+, PHISUM, dims=3)[:, :, 1]
    WTPSUM = reduce(+, WTPSUM, dims=3)[:, :, 1]
# end # @timeit to "reduce"
# @timeit to "compute" begin
    RHO[WTPSUM.>0.0] .= RHOSUM[WTPSUM.>0.0] ./ WTPSUM[WTPSUM.>0.0]
    RHOCP[WTPSUM.>0.0] .= RHOCPSUM[WTPSUM.>0.0] ./ WTPSUM[WTPSUM.>0.0]
    ALPHA[WTPSUM.>0.0] .= ALPHASUM[WTPSUM.>0.0] ./ WTPSUM[WTPSUM.>0.0]
    ALPHAF[WTPSUM.>0.0] .= ALPHAFSUM[WTPSUM.>0.0] ./ WTPSUM[WTPSUM.>0.0]
    HR[WTPSUM.>0.0] .= HRSUM[WTPSUM.>0.0] ./ WTPSUM[WTPSUM.>0.0]
    GGGP[WTPSUM.>0.0] .= GGGPSUM[WTPSUM.>0.0] .\ WTPSUM[WTPSUM.>0.0]
    SXX0[WTPSUM.>0.0] .= SXXSUM[WTPSUM.>0.0] ./ WTPSUM[WTPSUM.>0.0]
    tk1[WTPSUM.>0.0] .= TKSUM[WTPSUM.>0.0] ./ RHOCPSUM[WTPSUM.>0.0]
    PHI[WTPSUM.>0.0] .= PHISUM[WTPSUM.>0.0] ./ WTPSUM[WTPSUM.>0.0]
    BETTAPHI[WTPSUM.>0.0] .= GGGP[WTPSUM.>0.0] .\ PHI[WTPSUM.>0.0]
# end # @timeit to "compute"
# end # @timeit to "compute_p_node_properties!"
    return nothing
end # function compute_p_node_properties!


"""
Apply insulating boundary conditions to given array.

[x x x x x x        [a a b c d d

 x a b c d x         a a b c d d

 x e f g h x   ->    e e f g h h

 x x x x x x]        e e f g h h]
 
# Details

    - t: array to apply insulating boundary conditions to

# Returns

    - nothing
"""
function apply_insulating_boundary_conditions!(t)
# @timeit to "apply_insulating_boundary_conditions!" begin
    Ny, Nx = size(t)
    if Ny>2 && Nx>2
        # upper boundary
        t[1, 2:Nx-1] .= t[2, 2:Nx-1]
        # lower boundary
        t[Ny, 2:Nx-1] .= t[Ny-1, 2:Nx-1]
        # left boundary
        t[:, 1] .= t[:, 2]
        # right boundary
        t[:, Nx] .= t[:, Nx-1]
    end
# end # @timeit to "apply_insulating_boundary_conditions!"
    return nothing
end


"""
Compute gravity solution in P nodes to obtain
gravitational accelerations gx for Vx nodes, gy for Vy nodes.

$(SIGNATURES)

# Details

    - SP: solution vector
    - RP: right hand side vector
    - RHO: density at P nodes
    - xp: horizontal position of P nodes 
    - yp: vertical position of P nodes
    - gx: x gravitational acceleration at Vx nodes
    - gy: y gravitational acceleration at Vy nodes
    - sp: simulation parameters

# Returns

- nothing
"""
# function compute_gravity_solution!(LP, RP, RHO, xp, yp, gx, gy, sp)
function compute_gravity_solution!(SP, RP, RHO, xp, yp, gx, gy, sp)
# @timeit to "compute_gravity_solution!" begin
    @unpack Nx,
        Ny,
        Nx1,
        Ny1,
        xcenter,
        ycenter,
        dx,
        dy,
        G = sp
    # fresh LHS sparse coefficient matrix
    LP = ExtendableSparseMatrix(Nx1*Ny1, Nx1*Ny1)
    # iterate over P nodes
    # @timeit to "build system" begin
    for j=1:1:Nx1, i=1:1:Ny1
        # define global index in algebraic space
        gk = (j-1) * Ny1 + i
        # decide if external / boundary points
        if (
            i==1 ||
            i==Ny1 ||
            j==1 ||
            j==Nx1 ||
            distance(xp[j], yp[i], xcenter, ycenter) > xcenter
        )
            # boundary condition: ϕ = 0
            updateindex!(LP, +, 1.0, gk, gk)
            RP[gk] = 0.0
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
            updateindex!(LP, +, 1.0/dx^2, gk, gk-Ny1) # Φ₁
            updateindex!(LP, +, 1.0/dy^2, gk, gk-1) # Φ₂
            updateindex!(LP, +, -2.0/dx^2 -2.0/dy^2, gk, gk) # Φ₃
            updateindex!(LP, +, 1.0/dy^2, gk, gk+1) # Φ₄
            updateindex!(LP, +, 1.0/dx^2, gk, gk+Ny1) # Φ₅
            # fill system of equations: RHS (11.11)
            RP[gk] = 4.0 * 2.0/3.0 * π * G * RHO[i, j]
        end
    end
    # end # @timeit to "build system"
    # @timeit to "solve system" begin
    # solve system of equations
    SP .= LP \ RP # implicit: flush!(LP)
    # end # @timeit to "solve system"
    # reshape solution vector to 2D array
    # @timeit to "reshape solution" begin
    ϕ = reshape(SP, Ny1, Nx1)
    # end # @timeit to "reshape solution"
    # @timeit to "compute accelerations" begin
    # gx = -∂ϕ/∂x (11.12)
    gx[:, 1:Nx] .= -diff(ϕ, dims=2) ./ dx
    # gy = -∂ϕ/∂y (11.13)   
    gy[1:Ny, :] .= -diff(ϕ, dims=1) ./ dy
    # end # @timeit to "compute accelerations"
# end # @timeit to "compute_gravity_solution!"
    return nothing
end # function compute_gravity_solution!


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
    - SXX0: σ₀′xx XX stress at P nodes
    - RHOX: density at Vx nodes
    - RHOY: density at Vy nodes
    - dx: horizontal grid spacing
    - dy: vertical grid spacing
    - dt: time step
    - Nx: number of horizontal basic grid points
    - Ny: number of vertical basic grid points
    - Nx1: number of horizontal Vx/Vy/P grid points
    - Ny1: number of vertical Vx/Vy/P grid points

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
    dx,
    dy,
    dt,
    Nx,
    Ny,
    Nx1,
    Ny1,
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
    @views @. SXYcomp = SXY0*ETA / (GGG*dt + ETA)
    @views @. SXXcomp = (
        SXX0*ETAP[1:Ny, 1:Nx] / (GGGP[1:Ny, 1:Nx]*dt + ETAP[1:Ny, 1:Nx])
    )
    @views @. SYYcomp = (
        -SXX0*ETAP[1:Ny, 1:Nx] / (GGGP[1:Ny, 1:Nx]*dt+ETAP[1:Ny, 1:Nx])
    )
    # density gradients
    @views @. dRHOXdx[:, 2:Nx] = (RHOX[:, 3:Nx1]-RHOX[:, 1:Nx1-2]) / 2 / dx
    @views @. dRHOXdy[2:Ny, :] = (RHOX[3:Ny1, :]-RHOX[1:Ny1-2, :]) / 2 / dy
    @views @. dRHOYdx[:, 2:Nx] = (RHOY[:, 3:Nx1]-RHOY[:, 1:Nx1-2]) / 2 / dx
    @views @. dRHOYdy[2:Ny, :] = (RHOY[3:Ny1, :]-RHOY[1:Ny1-2, :]) / 2 / dy
    return nothing
end # @timeit to "get_viscosities_stresses_density_gradients!()"
end # function get_viscosities_stresses_density_gradients!


"""
Assemble hydromechanical system of equations.

$(SIGNATURES)

# Details

    - ETAcomp: computational viscosity at basic nodes
    - ETAPcomp: computational viscosity at P nodes
    - SXYcomp: computational previous XY stress at basic nodes
    - SXXcomp: computational previous XX stress at P nodes
    - SYYcomp: computational previous YY stress at P nodes
    - dRHOXdx: total density gradient in x direction at Vx nodes
    - dRHOXdy: total density gradient in y direction at Vx nodes
    - dRHOYdx: total density gradient in x direction at Vy nodes
    - dRHOYdy: total density gradient in y direction at Vy nodes
    - RHOX: total density at Vx nodes
    - RHOY: total density at Vy nodes
    - RHOFX: fluid density at Vx nodes
    - RHOFY: fluid density at Vy nodes
    - RX: ηfluid/Kϕ at Vx nodes
    - RY: ηfluid/Kϕ at Vy nodes
    - ETAPHI: bulk viscosity at P nodes
    - BETTAPHI: bulk compressibility at P nodes
    - PHI: porosity at P nodes
    - gx: x gravitational acceleration at Vx nodes
    - gy: y gravitational acceleration at Vy nodes
    - pr0: previous total pressure at P nodes
    - pf0: previous fluid pressure at P nodes
    - dt: time step
    - L: ExtendableSparse matrix to store LHS coefficients
    - R: vector to store RHS coefficients
    - sp: simulation parameters

# Returns

    - nothing
"""
function assemble_hydromechanical_lse!(
    ETAcomp,
    ETAPcomp,
    SXYcomp,
    SXXcomp,
    SYYcomp,
    dRHOXdx,
    dRHOXdy,
    dRHOYdx,
    dRHOYdy,
    RHOX,
    RHOY,
    RHOFX,
    RHOFY,
    RX,
    RY,
    ETAPHI,
    BETTAPHI,
    PHI,
    gx,
    gy,
    pr0,
    pf0,
    dt,
    L,
    R,
    sp
)
@timeit to "assemble_hydromechanical_lse()" begin
    @unpack Nx,
        Ny,
        Nx1,
        Ny1,
        dx,
        dy,
        vxleft,
        vxright,
        vytop,
        vybottom,
        bctop,
        bcbottom,
        bcleft,
        bcright,
        bcftop,
        bcfbottom,
        bcfleft,
        bcfright,
        pscale,
        psurface,
        etaphikoef = sp
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
            #                       kvx-6
            #                        Vx₂
            #                         |
            #             kvy-6   ETA(i-1,j) kvy+6⋅Ny1-6
            #              Vy₁    GGG(i-1,j)   Vy₃
            #                    SXY0(i-1,j)
            #                        ETA₁ (ETAcomp)
            #                       basic₁                       
            #                         |
            #             GGGP(i,j)   |    GGGP(i,j+1)
            #             ETAP(i,j)   |    ETAP(i,j+1) 
            #   kvx-6⋅Ny1 SXX0(i,j)  kvx   SXX0(i,j+1) kvx+6⋅Ny1
            #     Vx₁-------P₁-------Vx₃-------P₂-------Vx₅
            #              kpm        |     kpm+6⋅Ny1
            #             ETAPcomp    |    ETAPcomp    
            #                         |
            #              kvy     ETA(i,j)  kvy+6⋅Ny1
            #              Vy₂     GGG(i,j)    Vy₄
            #                     SXY0(i,j)
            #                        ETA₂ (ETAcomp)
            #                       basic₂
            #                         |
            #                       kvx+6
            #                        Vx₄
            #
            updateindex!(L, +, ETAPcomp[i, j]/dx^2, kvx, kvx-6*Ny1) # Vx₁
            updateindex!(L, +, ETAcomp[i-1, j]/dy^2, kvx, kvx-6) # Vx₂
            updateindex!(
                L,
                +,
                (
                    -(ETAPcomp[i, j]+ETAPcomp[i, j+1])/dx^2
                    -(ETAcomp[i-1, j]+ETAcomp[i,j])/dy^2
                    -dRHOXdx[i, j]*gx[i, j]*dt
                ),
                kvx,
                kvx
            ) # Vx₃
            updateindex!(L, +, ETAcomp[i, j]/dy^2, kvx, kvx+6) # Vx₄
            updateindex!(L, +, ETAPcomp[i, j+1]/dx^2, kvx, kvx+6*Ny1) # Vx₅
            updateindex!(
                L,
                +,
                (
                    ETAPcomp[i, j]/dx/dy
                    -ETAcomp[i, j]/dx/dy 
                    -dRHOXdy[i, j]*gx[i, j]*dt/4
                ),
                kvx,
                kvy
            ) # Vy₂
            updateindex!(
                L,
                +,
                (
                    -ETAPcomp[i, j+1]/dx/dy
                    +ETAcomp[i, j]/dx/dy
                    -dRHOXdy[i, j]*gx[i, j]*dt/4
                ),
                kvx,
                kvy+6*Ny1
            ) # Vy₄
            updateindex!(
                L,
                +,
                (
                    -ETAPcomp[i, j]/dx/dy
                    +ETAcomp[i-1, j]/dx/dy
                    -dRHOXdx[i, j]*gy[i, j]*dt/4
                ),
                kvx,
                kvy-6
            ) # Vy₁
            updateindex!(
                L,
                +,
                (
                    ETAPcomp[i, j+1]/dx/dy
                    -ETAcomp[i-1, j]/dx/dy
                    -dRHOXdx[i, j]*gy[i, j]*dt/4
                ),
                kvx,
                kvy+6*Ny1-6
            ) # Vy₃
            updateindex!(L, +, pscale/dx, kvx, kpm) # P₁
            updateindex!(L, +, -pscale/dx, kvx, kpm+6*Ny1) # P₂
            R[kvx] = (
                -RHOX[i, j]*gx[i, j]
                -(SXYcomp[i, j]- SXYcomp[i-1, j])/dy
                -(SXXcomp[i, j+1]-SXXcomp[i, j])/dx
            ) # RHS
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
            #                       kvy-6
            #                        Vy₂
            #                         |
            #                      SYY0(i,j)
            #                      ETAP(i,j)
            #           kvx-6⋅Ny1  GGGP(i,j)   kvx
            #              Vx₁        P₁       Vx₃
            #                        kpm 
            #                         |       
            #   kvy-6⋅Ny1 ETA(i,j-1) kvy     ETA(i,j)  kvy+6⋅Ny1
            #     Vy₁------ETA₁------Vv₃ ------ETA₂-------Vy₅
            #             GGG(i,j-1)  |      GGG(i,j)     
            #            SXY0(i,j-1) kpm+6  SXY0(i,j)
            #                         P₂
            #        kvx-6⋅Ny1+6  ETAP(i+1,j)  kvx+6
            #              Vx₂    GGGP(i+1,j)  Vx₄
            #                     SYY0(i+1,j)
            #                         |
            #                       kvy+6
            #                        Vy₄
            #
            updateindex!(L, +, ETAcomp[i, j-1]/dx^2, kvy, kvy-6*Ny1) # Vy₁
            updateindex!(L, +, ETAPcomp[i, j]/dy^2, kvy, kvy-6) # Vy₂
            updateindex!(
                L,
                +,
                (
                    -(ETAPcomp[i, j]+ETAPcomp[i+1, j])/dy^2
                    -(ETAcomp[i, j-1]+ETAcomp[i,j])/dx^2
                    -dRHOYdy[i, j]*gy[i, j]*dt
                ),
                kvy,
                kvy
            ) # Vy₃
            updateindex!(L, +, ETAPcomp[i+1, j]/dy^2, kvy, kvy+6) # Vy₄
            updateindex!(L, +, ETAcomp[i, j]/dx^2, kvy, kvy+6*Ny1) # Vy₅
            updateindex!(
                L,
                +,
                (
                    ETAPcomp[i, j]/dx/dy
                    -ETAcomp[i, j]/dx/dy
                    -dRHOYdx[i, j]*gy[i, j]*dt/4
                ),
                kvy,
                kvx
            ) # Vx₃
            updateindex!(
                L,
                +,
                (
                    -ETAPcomp[i+1, j]/dx/dy
                    +ETAcomp[i, j]/dx/dy
                    -dRHOYdx[i, j]*gy[i, j]*dt/4
                ),
                kvy,
                kvx+6
            ) # Vx₄
            updateindex!(
                L,
                +,
                (
                    -ETAPcomp[i, j]/dx/dy
                    +ETAcomp[i, j-1]/dx/dy
                    -dRHOYdy[i, j]*gx[i, j]*dt/4
                ),
                kvy,
                kvx-6*Ny1
            ) # Vx₁
            updateindex!(
                L,
                +,
                (
                    ETAPcomp[i+1, j]/dx/dy
                    -ETAcomp[i, j-1]/dx/dy
                    -dRHOYdy[i, j]*gx[i, j]*dt/4
                ),
                kvy,
                kvx-6*Ny1+6
            ) # Vx₂
            updateindex!(L, +, pscale/dy, kvy, kpm) # P₁
            updateindex!(L, +, -pscale/dy, kvy, kpm+6) # P₂
            R[kvy] = (
                -RHOY[i, j]*gy[i, j]
                -(SYYcomp[i+1, j]-SYYcomp[i, j])/dy
                -(SXYcomp[i, j]-SXYcomp[i, j-1])/dx
            ) # RHS
        end # Vy equation
        # P equation
        if i==1 || i==Ny1 || j==1 || j==Nx1
            # P equation external points: boundary conditions
            # all locations: ghost unknowns P=0 -> 1.0⋅P[i,j]=0.0
            updateindex!(L, +, 1.0, kpm, kpm)
            # R[kpm] = 0.0 # already done with initialization
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
            updateindex!(
                L,
                +,
                pscale/(1-PHI[i, j]) * (inv(ETAPHI[i, j])+BETTAPHI[i, j]/dt),
                kpm,
                kpm
            ) # P: Ptotal
            updateindex!(
                L,
                +,
                -pscale/(1-PHI[i, j]) * (inv(ETAPHI[i, j])+BETTAPHI[i, j]/dt),
                kpm,
                kpf
            ) # P: Pfluid
            R[kpm] = (pr0[i,j]-pf0[i,j]) / (1-PHI[i,j]) * BETTAPHI[i,j]/dt
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
            updateindex!(L, +, RX[i, j], kqx, kqx) # qxD
            updateindex!(L, +, -pscale/dx, kqx, kpf) # P₁
            updateindex!(L, +, pscale/dx, kqx, kqx+6*Ny1) # P₂
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
            updateindex!(L, +, RY[i, j], kqy, kqy) # qyD
            updateindex!(L, +, -pscale/dy, kqy, kpf) # P₁
            updateindex!(L, +, pscale/dy, kqy, kpf+6) # P₂
            R[kqy] = RHOFY[i, j] * gy[i, j]
        end # qyDarcy equation
        # Ptotal/Pfluid equation 
        if i==1 || i==Ny1 || j==1 || j==Nx1
            # Ptotal/Pfluid equation external points: boundary conditions
            # all locations: ghost unknowns P = 0 -> 1.0⋅P[i, j] = 0.0
            updateindex!(L, +, 1.0, kpf, kpf)
            # R[kpf] = 0.0 # already done with initialization
        elseif i==j==2
            # Ptotal/Pfluid real pressure boundary condition 'anchor'
            updateindex!(L, +, pscale, kpf, kpf)
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
            updateindex!(L, +, -1.0/dx, kpf, kqx-6*Ny1) # qxD₁
            updateindex!(L, +, 1.0/dx, kpf, kqx) # qxD₂
            updateindex!(L, +, -1.0/dy, kpf, kqy-6) # qyD₁
            updateindex!(L, +, 1.0/dy, kpf, kqy) # qyD₂
            updateindex!(
                L,
                +,
                -pscale/(1-PHI[i, j]) * (inv(ETAPHI[i, j])+BETTAPHI[i, j]/dt),
                 kpf,
                 kpm
            ) # Ptotal
            updateindex!(
                L,
                +,
                pscale/(1-PHI[i, j]) * (inv(ETAPHI[i, j])-BETTAPHI[i, j]/dt),
                kpf,
                kpf
            ) # Pfluid
            R[kpf] = -(pr0[i, j]-pf0[i, j]) / (1-PHI[i, j]) * BETTAPHI[i, j]/dt
        end # Ptotal/Pfluid equation
    end # for j=1:1:Nx1, i=1:1:Ny1
    flush!(L) # finalize CSC matrix
end # @timeit to "assemble_hydromechanical_lse()"
    return nothing
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
    - pscale: scaled pressure
    - Nx1: number of Vx/Vy/P nodes in horizontal x-direction
    - Ny1: number of Vx/Vy/P nodes in vertical y-direction

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
    pf,
    pscale,
    Nx1,
    Ny1
)
@timeit to "process_hydromechanical_solution!()" begin
    S = reshape(S, (:, Ny1, Nx1))
    vx .= S[1, :, :]
    vy .= S[2, :, :]
    pr .= S[3, :, :] .* pscale
    qxD .= S[4, :, :]
    qyD .= S[5, :, :]
    pf .= S[6, :, :] .* pscale
end # @timeit to "process_hydromechanical_solution!()"
    return nothing
end # function process_hydromechanical_solution!


"""
Recompute bulk viscosity at P nodes.

# Details

    - ETA: viscoplastic viscosity at basic nodes
    - ETAP: viscosity at P Nodes
    - ETAPHI: bulk viscosity at P Nodes
    - PHI: porosity at P Nodes
    - etaphikoef: coefficient: shear viscosity -> compaction viscosity

# Returns

    - nothing
"""
function recompute_bulk_viscosity!(ETA, ETAP, ETAPHI, PHI, etaphikoef)
@timeit to "recompute_bulk_viscosity!" begin    
    @views @. ETAP[2:end-1, 2:end-1] = 4.0 / (
        inv(ETA[1:end-1, 1:end-1]) +
        inv(ETA[2:end, 1:end-1]) +
        inv(ETA[1:end-1, 2:end]) +
        inv(ETA[2:end, 2:end])
    )
    @. ETAPHI = etaphikoef * ETAP / PHI
end # @timeit to "recompute_bulk_viscosity!"
    return nothing
end


"""
Compute porosity coefficient Aϕ = Dln[(1-ϕ)/ϕ]/Dt

$(SIGNATURES)

# Details

## In

    - ETAPHI: bulk viscosity at P Nodes
    - BETTAPHI: bulk compressibility at P nodes
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
function compute_Aϕ!(APHI, ETAPHI, BETTAPHI, PHI, pr, pf, pr0, pf0, dt)
@timeit to "compute_Aϕ!()" begin
    # reset APHI
    APHI .= 0.0
    @views @. APHI = (
        ((pr-pf)/ETAPHI + ((pr-pr0)-(pf-pf0))/dt*BETTAPHI) / (1-PHI) / PHI
    )
    return maximum(abs, APHI)
end # @timeit to "compute_Aϕ!()"
end # function compute_Aϕ!


"""
Compute current fluid velocity.

$(SIGNATURES)

# Details

## In

    - PHIX: porosity at Vx nodes
    - PHIY: porosity at Vy nodes
    - qxD: qx-Darcy flux at Vx nodes
    - qyD: qy-Darcy flux at Vy nodes
    - vx: solid velocity at Vx nodes
    - vy: solid velocity at Vy nodes
    - sp: simulation parameters

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
    vyf,
    sp
)
@timeit to "compute_fluid_velocities!()" begin
    @unpack Nx, Ny, Nx1, Ny1, bcftop, bcfbottom, bcfleft, bcfright = sp
    # vx velocity
    @. vxf[2:Ny, 1:Nx] = qxD[2:Ny, 1:Nx] / PHIX[2:Ny, 1:Nx]
    # top boundary
    vxf[1, :] = -bcftop  * vxf[2, :]
    # bottom boundary
    vxf[Ny1, :] = -bcfbottom * vxf[Ny, :]
    # vy velocity
    @. vyf[1:Ny, 2:Nx] = qyD[1:Ny, 2:Nx] / PHIY[1:Ny, 2:Nx]
    # left boundary
    vyf[:, 1] = -bcfleft * vyf[:, 2]
    # right boundary
    vyf[:, Nx1] = -bcfright * vyf[:, Nx]
    # adding solid velocity
    vxf .+= vx
    vyf .+= vy
end # @timeit to "compute_fluid_velocities!()"
    return nothing
end # function compute_fluid_velocities!


"""
Compute displacement time step.

$(SIGNATURES)

# Details

    - vx: solid vx velocity at Vx nodes
    - vy: solid vy velocity at Vy nodes
    - vxf: fluid vx velocity at Vx nodes
    - vyf: fluid vy velocity at Vy nodes
    - dt: current time step
    - dx: x-grid horizontal spacing
    - dy: y-grid vertical spacing
    - dxymax: maximum allowed grid spacing movement of markers per time step
    - aphimax: maximum observed porosity coefficient
    - dphimax: maximum allowed porosity ratio change per time step
   
# Returns

    - dtm: displacement time step
"""
function compute_displacement_timestep(
    vx,
    vy,
    vxf,
    vyf,
    dt,
    dx,
    dy,
    dxymax,
    aphimax,
    dphimax
)
@timeit to "compute_displacement_timestep()" begin
    maxvx = maximum(abs, vx)
    maxvy = maximum(abs, vy)
    maxvxf = maximum(abs, vxf)
    maxvyf = maximum(abs, vyf)    
    dtm = ifelse(dt*maxvx > dxymax*dx, dxymax*dx/maxvx, dt)
    dtm = ifelse(dtm*maxvy > dxymax*dy, dxymax*dy/maxvy, dtm)
    dtm = ifelse(dtm*maxvxf > dxymax*dx, dphimax*dx/maxvxf, dtm)
    dtm = ifelse(dtm*maxvyf > dxymax*dy, dphimax*dy/maxvyf, dtm)
    dtm = ifelse(dtm*aphimax > dphimax, dphimax/aphimax, dtm)
end # @timeit to "compute_displacement_timestep()"
    return dtm
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
    - dtm: displacement time step
    - sp: simulation parameters 

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
    dtm,
    sp
)
@timeit to "compute_stress_strainrate!()" begin
    @unpack Nx, Ny, Nx1, Ny1, dx, dy = sp
    # ϵxy, σxy, Δσxy at basic nodes
    EXY .= 0.5.*(diff(vx, dims=1)[:, 1:Nx]./dy .+ diff(vy, dims=2)[1:Ny, :]./dx)
    @. SXY = 2*ETA*EXY*GGG*dtm/(GGG*dtm+ETA) + SXY0*ETA/(GGG*dtm+ETA)
    @. DSXY = SXY - SXY0
    # ϵxx, σ′xx at P nodes
    # @. DIVV[2:end, 2:end] = 
    #     diff(vx, dims=2)[2:end, :]/dx + diff(vy, dims=1)[:, 2:end]/dy
    EXX[2:Ny1, 2:Nx1] .= ( 
        0.5 .* (
            diff(vx, dims=2)[2:Ny1, :]./dx .- diff(vy, dims=1)[:, 2:Nx1]./dy
        )
    )
    @. SXX = 2.0*ETAP*EXX*GGGP*dtm/(GGGP*dtm+ETAP) + SXX0*ETAP/(GGGP*dtm+ETAP)
    @. DSXX = SXX - SXX0
    @. EII[2:Ny, 2:Nx] = sqrt(
        EXX[2:Ny, 2:Nx]^2 + (
            (
                EXY[2:Ny, 2:Nx]
                +EXY[1:Ny-1,2:Nx]
                +EXY[2:Ny,1:Nx-1]
                +EXY[1:Ny-1,1:Nx-1]
            )/4.0
        )^2
    )
    @. SII[2:Ny, 2:Nx] = sqrt(
        SXX[2:Ny, 2:Nx]^2 + (
            (
                SXY[2:Ny, 2:Nx]
                +SXY[1:Ny-1,2:Nx]
                +SXY[2:Ny,1:Nx-1]
                +SXY[1:Ny-1,1:Nx-1]
            )/4.0
        )^2
    )
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
    - Nx: number of basic nodes in horizontal x direction
    - Ny: number of basic nodes in vertical y direction
    - Nx1: number of Vx/Vy/P nodes in horizontal x direction
    - Ny1: number of Vx/Vy/P nodes in vertical y direction

# Returns

    - nothing
"""
function symmetrize_p_node_observables!(
    SXX,
    APHI,
    PHI,
    pr,
    pf,
    ps,
    Nx,
    Ny,
    Nx1,
    Ny1
)
@timeit to "symmetrize_p_node_observables!()" begin
    # top boundary
    @. begin
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
    ps = (pr-pf*PHI) / (1-PHI)
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
    - SIIB: second stress invariant at basic nodes
    - siiel: second invariant for purely elastic stress buildup at basic nodes 
    - prB: interpolated total pressure at basic nodes 
    - pfB: interpolated fluid pressure at basic nodes 
    - syieldc: confined fractures yielding stress at basic nodes 
    - syieldt: tensile fractures yielding stress at basic nodes 
    - syield: non-negative maximum yielding stress at basic nodes
    - etapl: stress-based viscoplastic viscosity at basic nodes 
    - YNY: plastic yielding status at basic nodes 
    - YNY5: plastic iterations plastic yielding status at basic nodes
    - YERRNOD: vector of summed yielding errors of nodes over plastic iterations
    - DSY: (SIIB-syield) at basic nodes
    - YNPL: plastic iterations plastic yielding status at basic nodes
    - DDD: plastic iterations (SIIB-syield)² at basic nodes
    - dt: time set
    - iplast: plastic iteration step 
    - sp: static simulation parameters

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
    SIIB,
    siiel,
    prB,
    pfB,
    syieldc,
    syieldt,
    syield,
    etapl,
    YNY,
    YNY5,
    YERRNOD,
    DSY,
    YNPL,
    DDD,
    dt,
    iplast,
    sp
)
@timeit to "compute_nodal_adjustment!()" begin
    @unpack Nx, Ny, Nx1, Ny1, etamin, etamax, etawt, yerrmax, nplast = sp
    # reset / setup
    @. begin
    ETA5 = copy(ETA0)
    YNY5 = 0
    DSY = 0.0
    YNPL = 0
    DDD = 0.0
    end
    # second stress invariant at basic nodes
    @views @. SIIB = sqrt(
        SXY^2 + (
            0.25 * (
                SXX[1:Ny, 1:Nx]
                +SXX[2:Ny1, 1:Nx]
                +SXX[1:Ny, 2:Nx1]
                +SXX[2:Ny1, 2:Nx1]
            )
        )^2
    )
    # second invariant for purely elastic stress buildup at basic nodes
    @views @. siiel = SIIB * (GGG*dt+ETA) / ETA
    # interpolate total and fluid pressure at basic nodes
    @views @. prB = 0.25 * (
        pr[1:Ny, 1:Nx] + pr[2:Ny1, 1:Nx] + pr[1:Ny, 2:Nx1] + pr[2:Ny1, 2:Nx1]
    )   
    @views @. pfB = 0.25 * (
        pf[1:Ny, 1:Nx] + pf[2:Ny1, 1:Nx] + pf[1:Ny, 2:Nx1] + pf[2:Ny1, 2:Nx1]
    )
    # yielding stress: confined fracture
    @views @. syieldc = COH + FRI * (prB-pfB)
    # yielding stress: tensile fracture
    @views @. syieldt = TEN + (prB-pfB)
    # non-negative maximum yielding stress
    positive_max!(syieldc, syieldt, syield)
    # update error for previous yielding nodes
    @views @. DSY[YNY>0] = SIIB[YNY>0] - syield[YNY>0]
    @views @. DDD[YNY>0] += DSY[YNY>0]^2
    @views @. YNPL[YNY>0] .= 1
    # new viscosity for basic nodes
    @views @. etapl = dt * GGG * syield / (siiel-syield)
    # correcting viscosity for yielding
    # recompute nodal viscosity, apply min/max viscosity cutoffs
    @views @. ETA5[syield<siiel && etapl<ETA0] = (
        etapl[syield<siiel && etapl<ETA0]^(1.0-etawt)
        * ETA[syield<siiel && etapl<ETA0]^etawt
    )
    @views @. ETA5[syield<siiel && etapl<ETA0 && ETA5>etamax] = etamax
    @views @. ETA5[syield<siiel && etapl<ETA0 && ETA5<etamin] = etamin
    # mark yielding nodes
    @views @. YNY5[syield<siiel && etapl<ETA0] = 1
    # update error for new yielding nodes
    @views @. DSY[syield<siiel && etapl<ETA0 && YNPL==0] = (
        SIIB[syield<siiel && etapl<ETA0 && YNPL==0]
         -syield[syield<siiel && etapl<ETA0 && YNPL==0]
    )
    @views @. DDD[syield<siiel && etapl<ETA0 && YNPL==0] +=  
        DSY[syield<siiel && etapl<ETA0 && YNPL==0]^2
    @views @. YNPL[syield<siiel && etapl<ETA0 && YNPL==0] = 1
    if sum(YNPL) > 0
        YERRNOD[iplast] = sqrt(sum(DDD)/sum(YNPL))
    end
    # return plastic iteration completeness
    return sum(YNPL)==0 || YERRNOD[iplast]<yerrmax || iplast==nplast
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
    @inbounds for i in eachindex(A)
        C[i] = max(0, ifelse(A[i] > B[i], A[i], B[i]))
    end
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
    - dt: current time step
    - dtkoef: coefficient by which to decrease computational time step
    - dtstep: minimum of plastic iterations before adjusting time step
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
    dt,
    dtkoef,
    dtstep,
    iplast
)
@timeit to "finalize_plastic_iteration_pass!()" begin
    if iplast % dtstep == 0
        # dtstep plastic iterations performed without reaching targets:
        # decrease time step and reset to previous viscoplastic viscosity
        dt /= dtkoef
        ETA .= copy(ETA00)
        YNY .= copy(YNY00)
    else
        # perform next plastic iteration pass with new viscoplastic viscosity
        ETA .= copy(ETA5)
        YNY .= copy(YNY5)
    end
    return dt
end # @timeit to "finalize_plastic_iteration_pass!()"
    end # function finalize_plastic_iteration_pass


"""
Main simulation loop: run calculations with timestepping.

$(SIGNATURES)

# Details

    - markers: arrays containing all marker properties
    - sp: static simulation parameters

# Returns
    
    - nothing
"""
function simulation_loop(sp::StaticParameters)
    # -------------------------------------------------------------------------
    # unpack static simulation parameters
    # -------------------------------------------------------------------------
    @unpack xsize, ysize,
        Nx, Ny,
        Nx1, Ny1,
        dx, dy,
        jmin_basic, jmax_basic,
        imin_basic, imax_basic,
        jmin_vx, jmax_vx,
        imin_vx, imax_vx,
        jmin_vy, jmax_vy,
        imin_vy, imax_vy,
        jmin_p, jmax_p,
        imin_p, imax_p,
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
        start_hrsolidm,
        start_hrfluidm,
        gggsolidm,
        frictsolidm,
        cohessolidm,
        tenssolidm,
        kphim0,
        etaphikoef,
        phim0,
        tmsilicate,
        tmiron,
        etamin,
        nplast,
        dtelastic,
        dtkoef,
        dtkoefup,
        start_step,
        nsteps,
        start_time, 
        endtime,
        start_marknum = sp

# @timeit to "simulation_loop setup" begin
    # -------------------------------------------------------------------------
    # set up dynamic simulation parameters from given static parameters
    # -------------------------------------------------------------------------
    timestep,
    dt,
    timesum,
    marknum,
    hrsolidm,
    hrfluidm,
    YERRNOD = setup_dynamic_simulation_parameters(sp)
    
    # -------------------------------------------------------------------------
    # set up staggered grid
    # -------------------------------------------------------------------------
    x, y, xvx, yvx, xvy, yvy, xp, yp = setup_staggered_grid_geometry(sp)
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
        BETTAPHI,
        PHI,
        APHI,
        FI
    ) = setup_staggered_grid_properties(sp)

    # -------------------------------------------------------------------------
    # set up markers
    # -------------------------------------------------------------------------
    # primary marker arrays: initialized at beginning
    # horizontal marker coordinate [m]
    xm = zeros(Float64, marknum)
    # vertical marker coordinate [m]
    ym = zeros(Float64, marknum)
    # marker material type
    tm = zeros(Int8, marknum)
    # marker temperature [K]
    tkm = zeros(Float64, marknum)
    # marker σ′xx [Pa]
    sxxm = zeros(Float64, marknum)
    # marker σxy [Pa]
    sxym = zeros(Float64, marknum)
    # marker porosity ϕ
    phim = zeros(Float64, marknum)
    # marker viscoplastic viscosity [Pa]
    etavpm = zeros(Float64, marknum)
    # marker total density
    rhototalm = zeros(Float64, marknum)
    # marker total volumetric heat capacity
    rhocptotalm = zeros(Float64, marknum)
    # marker solid viscosity
    etasolidcur = zeros(Float64, marknum)
    # marker fluid viscosity
    etafluidcur = zeros(Float64, marknum)
    # marker total viscosity
    etatotalm = zeros(Float64, marknum)
    # marker total radiogenic heat production
    hrtotalm = zeros(Float64, marknum)
    # marker total thermal conductivity
    ktotalm = zeros(Float64, marknum)
    # marker total thermal energy
    tkm_rhocptotalm = zeros(Float64, marknum)
    # marker fluid viscosity over permeability
    etafluidcur_inv_kphim = zeros(Float64, marknum)
    # marker temperature 
    tkm = zeros(Float64, marknum)
    # marker inverse of total shear modulus
    inv_gggtotalm = zeros(Float64, marknum)
    # marker total friction coefficient
    fricttotalm = zeros(Float64, marknum)
    # marker total compressive strength
    cohestotalm = zeros(Float64, marknum)
    # marker total tensile strength
    tenstotalm = zeros(Float64, marknum)
    # marker fluid density
    rhofluidcur = zeros(Float64, marknum)
    # marker solid thermal expansion coefficient
    alphasolidcur = zeros(Float64, marknum)
    # marker fluid thermal expansion coefficient
    alphafluidcur = zeros(Float64, marknum)
    
    # define initial markers: coordinates, material type, and properties    
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
        etafluidcur,
        tkm,
        inv_gggtotalm,
        fricttotalm,
        cohestotalm,
        tenstotalm,
        rhofluidcur,
        alphasolidcur,
        alphafluidcur,
        sp
    )

    # -------------------------------------------------------------------------
    # set up of matrices for global gravity/thermal/hydromechanical solutions
    # -------------------------------------------------------------------------
    # hydromechanical solution: LHS coefficient matrix
    # L = ExtendableSparseMatrix(Nx1*Ny1*6, Nx1*Ny1*6)
    # hydromechanical solution: RHS vector
    R = zeros(Float64, Nx1*Ny1*6)
    # thermal solution: LHS coefficient matrix
    LT = ExtendableSparseMatrix(Nx1*Ny1, Nx1*Ny1)
    # thermal solution: RHS vector
    RT = zeros(Float64, Nx1*Ny1)
    # gravity solution: LHS coefficient matrix
    # LP = ExtendableSparseMatrix(Nx1*Ny1, Nx1*Ny1)
    # gravity solution: RHS vector
    RP = zeros(Float64, Nx1*Ny1)
    # gravity solution: solution vector (->matrix)
    SP = zeros(Float64, Nx1*Ny1)
# end # @timeit to "simulation_loop setup"

    # -------------------------------------------------------------------------
    # iterate timesteps   
    # -------------------------------------------------------------------------
    nsteps = 10 # <======= remove for production
    p = Progress(
        nsteps,
        dt=0.5,
        barglyphs=BarGlyphs(
            '|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',), barlen=10)
    for timestep = start_step:1:nsteps
# @timeit to "set up interpolation arrays" begin
        # ---------------------------------------------------------------------
        # set up interpolation arrays
        # ---------------------------------------------------------------------
        # basic nodes
        ETA0SUM = zeros(Ny, Nx, nthreads())
        ETASUM = zeros(Ny, Nx, nthreads())
        GGGSUM = zeros(Ny, Nx, nthreads())
        SXYSUM = zeros(Ny, Nx, nthreads())
        COHSUM = zeros(Ny, Nx, nthreads())
        TENSUM = zeros(Ny, Nx, nthreads())
        FRISUM = zeros(Ny, Nx, nthreads())
        WTSUM = zeros(Ny, Nx, nthreads())
        # Vx nodes
        RHOXSUM = zeros(Ny1, Nx1, nthreads())
        RHOFXSUM = zeros(Ny1, Nx1, nthreads())
        KXSUM = zeros(Ny1, Nx1, nthreads())
        PHIXSUM = zeros(Ny1, Nx1, nthreads())
        RXSUM = zeros(Ny1, Nx1, nthreads())
        WTXSUM = zeros(Ny1, Nx1, nthreads())
        # Vy nodes
        RHOYSUM = zeros(Ny1, Nx1, nthreads())
        RHOFYSUM = zeros(Ny1, Nx1, nthreads())
        KYSUM = zeros(Ny1, Nx1, nthreads())
        PHIYSUM = zeros(Ny1, Nx1, nthreads())
        RYSUM = zeros(Ny1, Nx1, nthreads())
        WTYSUM = zeros(Ny1, Nx1, nthreads())
        # P Nodes
        RHOSUM = zeros(Ny1, Nx1, nthreads())
        RHOCPSUM = zeros(Ny1, Nx1, nthreads())
        ALPHASUM = zeros(Ny1, Nx1, nthreads())
        ALPHAFSUM = zeros(Ny1, Nx1, nthreads())
        HRSUM = zeros(Ny1, Nx1, nthreads())
        GGGPSUM = zeros(Ny1, Nx1, nthreads())
        SXXSUM = zeros(Ny1, Nx1, nthreads())
        TKSUM = zeros(Ny1, Nx1, nthreads())
        PHISUM = zeros(Ny1, Nx1, nthreads())
        WTPSUM = zeros(Ny1, Nx1, nthreads())
# end # @timeit to "set up interpolation arrays" 

        # ---------------------------------------------------------------------
        # calculate radioactive heating
        # ---------------------------------------------------------------------
        hrsolidm, hrfluidm = calculate_radioactive_heating(timesum, sp)

        # ---------------------------------------------------------------------
        # computer marker properties and interpolate to staggered grid nodes
        # ---------------------------------------------------------------------
        @threads for m=1:1:marknum
            compute_marker_properties!(
                m,
                tm,
                tkm,
                rhototalm,
                rhocptotalm,
                etasolidcur,
                etafluidcur,
                etatotalm,
                hrtotalm,
                ktotalm,
                tkm_rhocptotalm,
                etafluidcur_inv_kphim,
                hrsolidm,
                hrfluidm,
                phim,
                sp
            )
                    
            # interpolate marker properties to basic nodes
            i, j, weights = fix_weights(
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
            # ETA0SUM: viscous viscosity interpolated to basic nodes
            interpolate_to_grid!(i, j, weights, etatotalm[m], ETA0SUM)
            # ETASUM: viscoplastic viscosity interpolated to basic nodes
            interpolate_to_grid!(i, j, weights, etavpm[m], ETASUM)
            # GGGSUM: shear modulus interpolated to basic nodes
            interpolate_to_grid!(i, j, weights, inv_gggtotalm[m], GGGSUM)
            # SXYSUM: σxy shear stress interpolated to basic nodes
            interpolate_to_grid!(i, j, weights, sxym[m], SXYSUM)
            # COHSUM: compressive strength interpolated to basic nodes
            interpolate_to_grid!(i, j, weights, cohestotalm[m], COHSUM)
            # TENSUM: tensile strength interpolated to basic nodes
            interpolate_to_grid!(i, j, weights, tenstotalm[m], TENSUM)
            # FRISUM: friction  interpolated to basic nodes
            interpolate_to_grid!(i, j, weights, fricttotalm[m], FRISUM)
            # WTSUM: weight array for bilinear interpolation to basic nodes
            interpolate_to_grid!(i, j, weights, 1.0, WTSUM)

            # interpolate marker properties to Vx nodes
            i, j, weights = fix_weights(
                xm[m],
                ym[m],
                xvx,
                yvx,
                dx,
                dy,
                jmin_vx,
                jmax_vx,
                imin_vx,
                imax_vx
            )
            # RHOXSUM: density interpolated to Vx nodes
            interpolate_to_grid!(i, j, weights, rhototalm[m], RHOXSUM)
            # RHOFXSUM: fluid density interpolated to Vx nodes
            interpolate_to_grid!(i, j, weights, rhofluidcur[m], RHOFXSUM)
            # KXSUM: thermal conductivity interpolated to Vx nodes
            interpolate_to_grid!(i, j, weights, ktotalm[m], KXSUM)
            # PHIXSUM: porosity interpolated to Vx nodes
            interpolate_to_grid!(i, j, weights, phim[m], PHIXSUM)
            # RXSUM: ηfluid/kϕ interpolated to Vx nodes
            interpolate_to_grid!(i, j, weights, etafluidcur_inv_kphim[m], RXSUM)
            # WTXSUM: weight for bilinear interpolation to Vx nodes
            interpolate_to_grid!(i, j, weights, 1.0, WTXSUM)

            # interpolate marker properties to Vy nodes
            i, j, weights = fix_weights(
                xm[m],
                ym[m],
                xvy,
                yvy,
                dx,
                dy,
                jmin_vy,
                jmax_vy,
                imin_vy,
                imax_vy
            )
            # RHOYSUM: density interpolated to Vy nodes
            interpolate_to_grid!(i, j, weights, rhototalm[m], RHOYSUM)
            # RHOFYSUM: fluid density interpolated to Vy nodes
            interpolate_to_grid!(i, j, weights, rhofluidcur[m], RHOFYSUM)
            # KYSUM: thermal conductivity interpolated to Vy nodes
            interpolate_to_grid!(i, j, weights, ktotalm[m], KYSUM)
            # PHIYSUM: porosity interpolated to Vy nodes
            interpolate_to_grid!(i, j, weights, phim[m], PHIYSUM)
            # RYSUM: ηfluid/kϕ interpolated to Vy nodes
            interpolate_to_grid!(i, j, weights, etafluidcur_inv_kphim[m], RYSUM)
            # WTYSUM: weight for bilinear interpolation to Vy nodes
            interpolate_to_grid!(i, j, weights, 1.0, WTYSUM)
            
            # interpolate marker properties to P nodes
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
            # GGGPSUM: shear modulus interpolated to P nodes
            interpolate_to_grid!(i, j, weights, inv_gggtotalm[m], GGGPSUM)
            # SXXSUM: σ'xx interpolated to P nodes
            interpolate_to_grid!(i, j, weights, sxxm[m], SXXSUM)
            # RHOSUM: density interpolated to P nodes
            interpolate_to_grid!(i, j, weights, rhototalm[m], RHOSUM)
            # RHOCPSUM: volumetric heat capacity interpolated to P nodes
            interpolate_to_grid!(i, j, weights, rhocptotalm[m], RHOCPSUM)
            # ALPHASUM: thermal expansion interpolated to P nodes
            interpolate_to_grid!(i, j, weights, alphasolidcur[m], ALPHASUM)
            # ALPHAFSUM: fluid thermal expansion interpolated to P nodes
            interpolate_to_grid!(i, j, weights, alphafluidcur[m], ALPHAFSUM)
            # HRSUM: radioactive heating interpolated to P nodes
            interpolate_to_grid!(i, j, weights, hrtotalm[m], HRSUM)
            # PHISUM: porosity interpolated to P nodes
            interpolate_to_grid!(i, j, weights, phim[m], PHISUM)
            # TKSUM: heat capacity interpolated to P nodes
            interpolate_to_grid!(i, j, weights, tkm_rhocptotalm[m], TKSUM)
            # WTPSUM: weight for bilinear interpolation to P nodes
            interpolate_to_grid!(i, j, weights, 1.0, WTPSUM)
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
            BETTAPHI
        )

        # ---------------------------------------------------------------------
        # applying thermal boundary conditions for interpolated temperature
        # ---------------------------------------------------------------------
        apply_insulating_boundary_conditions!(tk1)

        # ---------------------------------------------------------------------
        # compute gravity solution
        # compute gravitational acceleration
        # ---------------------------------------------------------------------
        # compute_gravity_solution!(LP, RP, RHO, xp, yp, gx, gy, sp)
        compute_gravity_solution!(SP, RP, RHO, xp, yp, gx, gy, sp)

        # ---------------------------------------------------------------------
        # # probe increasing computational timestep
        # ---------------------------------------------------------------------
        dt = min(dt*dtkoefup, dtelastic)

        # ---------------------------------------------------------------------
        # # save initial viscosity, yielding nodes
        # ---------------------------------------------------------------------
        ETA00, ETA = ETA, ETA00
        YNY00, YNY = YNY, YNY00

# Mon 18/Tue 19/Wed 20       
        # ---------------------------------------------------------------------
        # # perform plastic iterations
        # ---------------------------------------------------------------------
        # 768-1316
        if timestep == 1
            # no elastic compaction during first timestep
            BETTAPHI .= 0.0
        end
        # perform plastic iterations
        for iplast=1:1:nplast
            # recompute bulk viscosity at pressure nodes
            recompute_bulk_viscosity!(ETA, ETAP, ETAPHI, PHI, etaphikoef)
            # compute computational viscosities, stresses, and density gradients
            get_viscosities_stresses_density_gradients!(
                ETA,
                ETAP,
                GGG,
                GGGP,
                SXY0,
                SXX0,
                RHOX,
                RHOY,
                dx,
                dy,
                dt,
                Nx,
                Ny,
                Nx1,
                Ny1,
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
            # assemble hydromechanical system of equations
            assemble_hydromechanical_lse!(
                ETAcomp,
                ETAPcomp,
                SXYcomp,
                SXXcomp,
                SYYcomp,
                dRHOXdx,
                dRHOXdy,
                dRHOYdx,
                dRHOYdy,
                RHOX,
                RHOY,
                RHOFX,
                RHOFY,
                RX,
                RY,
                ETAPHI,
                BETTAPHI,
                PHI,
                gx,
                gy,
                pr0,
                pf0,
                dt,
                L,
                R,
                sp
            )
            # solve hydromechanical system of equations
@timeit to "solve system" begin
            S = L \ R
end # @timeit to "solve system"
            # obtain hydromechanical observables from solution
            process_hydromechanical_solution!(
                S,
                vx,
                vy,
                pr,
                qxD,
                qyD,
                pf,
                pscale,
                Nx1,
                Ny1
            )
            # compute Aϕ = Dln[(1-PHI)/PHI]/Dt
            aphimax = compute_Aϕ!(
                APHI,
                ETAPHI,
                BETTAPHI,
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
                vyf,
                sp
            )
            # define displacement timestep dtm
            dtm = compute_displacement_timestep(
                vx,
                vy,
                vxf,
                vyf,
                dt,
                dx,
                dy,
                dxymax,
                aphimax,
                dphimax
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
                dtm,
                sp
            )
            # recompute Dln[(1-PHI)/PHI]/Dt
            aphimax = compute_Aϕ!(
                APHI,
                ETAPHI,
                BETTAPHI,
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
                ps,
                Nx,
                Ny,
                Nx1,
                Ny1
            )
            # save nodal stress changes - RMK: not required in code
            # DSXX0 = copy(DSXX)
            # DSXY0 = copy(DSXY)
            # nodal adjustment
            if compute_nodal_adjustment(
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
                SIIB,
                siiel,
                prB,
                pfB,
                syieldc,
                syieldt,
                syield,
                etapl,
                YNY,
                YNY5,
                YERRNOD,
                DSY,
                YNPL,
                DDD,
                dt,
                iplast,
                sp
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
                    dt,
                    dtkoef,
                    dtstep,
                    iplast
                )
            end
        end # for iplast=1:1:nplast

# Thu 21
        # ---------------------------------------------------------------------
        # # interpolate updated viscoplastic viscosity to markers
        # ---------------------------------------------------------------------
        # 1320-1369
        # for m = 1:1:marknum
        #     # ~50 lines MATLAB 
        # end


# Thu 21
        # ---------------------------------------------------------------------
        # # apply subgrid stress diffusion to markers
        # ---------------------------------------------------------------------
        # 1374-1467
        # for m = 1:1:marknum
        #     # ~100 lines MATLAB 
        # end


# Fri 22
        # ---------------------------------------------------------------------
        # # compute DSXXsubgrid, DSXYsubgrid
        # ---------------------------------------------------------------------
        # 1471-1492
        # compute_dsxx_dsxy_subgrids!(sp, dp, interp_arrays)


# Fri 22
        # ---------------------------------------------------------------------
        # # interpolate DSXX, DSXY to markers
        # ---------------------------------------------------------------------
        # 1495-1547
        # for m = 1:1:marknum
        #     # ~50 lines MATLAB 
        # end

# Fri 22
        # ---------------------------------------------------------------------
        # # compute shear heating HS in P nodes
        # ---------------------------------------------------------------------
        # 1551-1563
        # compute_HS_p_nodes!(sp, dp, interp_arrays)


# Fri 22
        # ---------------------------------------------------------------------
        # # compute adiabatic heating HA in P nodes
        # ---------------------------------------------------------------------
        # 1567-1614
        # compute_HA_p_nodes!(sp, dp, interp_arrays)


# Sat 23
        # ---------------------------------------------------------------------
        # # perform thermal iterations
        # ---------------------------------------------------------------------
        # 1618-1725
        # # ~100 lines MATLAB

# Sat 23 
        # ---------------------------------------------------------------------
        # # apply subgrid temperature diffusion on markers
        # ---------------------------------------------------------------------
        # 1729-1786
        # for m = 1:1:marknum
        #     # ~50 lines MATLAB
        # end


# Sat 23
        # ---------------------------------------------------------------------
        # # compute DTsubgrid
        # ---------------------------------------------------------------------
        # 1787-1799
        # compute_DT_subgrid!(sp, dp, interp_arrays)


# Mon 25
        # ---------------------------------------------------------------------
        # # interpolate DT to markers
        # ---------------------------------------------------------------------
        # 1803-1833
        # for m = 1:1:marknum
        #     # ~30 lines MATLAB
        # end


# Mon 25
        # ---------------------------------------------------------------------
        # # update porosity on markers
        # ---------------------------------------------------------------------
        # 1838-1873
        # for m = 1:1:marknum
        #     # ~30 lines MATLAB
        # end


# Mon 25
        # ---------------------------------------------------------------------
        # # compute fluid velocity in P nodes including boundary conditions
        # ---------------------------------------------------------------------
        # 1877-1905
        # compute_v_fluid_p_nodes(sp, dp, interp_arrays)


# Mon 25
        # ---------------------------------------------------------------------
        # # compute velocity in P nodes
        # ---------------------------------------------------------------------
        # 1909-1937
        # compute_v_p_nodes!(sp, dp, interp_arrays)


# Mon 25
        # ---------------------------------------------------------------------
        # # compute rotation rate in basic nodes
        # ---------------------------------------------------------------------
        # 1940-1944
        # compute_ω_basic_nodes!(sp, dp, interp_arrays)


# Tue 26
        # ---------------------------------------------------------------------
        # # move markers with RK4
        # ---------------------------------------------------------------------
        # 1947-2227
        # for m = 1:1:marknum
        #     # ~300 lines MATLAB
        # end


# Wed 26
        # ---------------------------------------------------------------------
        # # backtrack P nodes: Ptotal with RK4
        # ---------------------------------------------------------------------
        # 2231-2360

# Thu 27
        # ---------------------------------------------------------------------
        # # backtrack P nodes: Pfluid with RK1/2/3
        # ---------------------------------------------------------------------
        # 2364-2487

# Fri 28
        # ---------------------------------------------------------------------
        # # replenish sparse areas with additional markers
        # ---------------------------------------------------------------------
        # 2489-2604
        # for m = 1:1:marknum
        #     # ~100 lines MATLAB
        # end

# Fri 28
        # ---------------------------------------------------------------------
        # # update timesum
        # ---------------------------------------------------------------------

# Fri 28
        # ---------------------------------------------------------------------
        # # save data for analysis and visualization
        # ---------------------------------------------------------------------

# Fri 28
        # ---------------------------------------------------------------------
        # # save old stresses - RMK: not used anywhere in code ?
        # ---------------------------------------------------------------------
        # # sxxm00 = sxxm 
        # # sxym00 = sxym    

# Fri 28
        # ---------------------------------------------------------------------
        # finish timestep
        # ---------------------------------------------------------------------
        next!(p)
        # if timestep % 20 == 0
        #     println("timestep: ", timestep)
        # end

        if timesum > endtime
            break
        end

    end # for timestep = startstep:1:nsteps
end # function simulation loop


"""
Runs the simulation with the given parameters.

# Details

    - xsize: size of the domain in horizontal x direction [m]
    - ysize: size of the domain in vertical y direction [m]
    - rplanet: radius of the planet [m]
    - rcrust: radius of the crust [m]
    - Nx: number of basic grid nodes in horizontal x direction
    - Ny: number of basic grid nodes in vertical y direction
    - Nxmc: initial number of markers per cell in horizontal x direction
    - Nymc: initial number of markers per cell in vertical y direction

# Returns

    - exit code
"""
function run_simulation(
    xsize=140_000.0,
    ysize=140_000.0,
    rplanet=50_000.0,
    rcrust=48_000.0,
    Nx=141,
    Ny=141,
    Nxmc=4,
    Nymc=4
)
    reset_timer!(to)
    sp = StaticParameters(
        xsize=xsize,
        ysize=ysize,
        rplanet=rplanet,
        rcrust=rcrust,
        Nx=Nx,
        Ny=Ny,
        Nxmc=Nxmc,
        Nymc=Nymc
        )
    simulation_loop(sp)
    show(to)
end

end # module HydrologyPlanetesimals
module HydrologyPlanetesimals

using Base.Threads
using DocStringExtensions
using ExtendableSparse
using LinearAlgebra
using MAT
using Parameters
using ProgressMeter
using SparseArrays
using StaticArrays
using TimerOutputs

export run_simulation

include("constants.jl")
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
    # basic nodes
    "horizontal coordinates of basic grid points [m]"
    x = SVector{Nx, Float64}([j for j = 0:dx:xsize])
    "vertical coordinates of basic grid points [m]"
    y = SVector{Ny, Float64}([i for i = 0:dy:ysize])
    # Vx nodes
    "horizontal coordinates of vx grid points [m]"
    xvx = SVector{Ny1, Float64}([j for j = 0:dx:xsize+dy])
    "vertical coordinates of vx grid points [m]"
    yvx = SVector{Nx1, Float64}([i for i = -dy/2:dy:ysize+dy/2])
    # Vy nodes
    "horizontal coordinates of vy grid points [m]"
    xvy = SVector{Nx1, Float64}([j for j = -dx/2:dx:xsize+dx/2])
    "vertical coordinates of vy grid points [m]"
    yvy = SVector{Ny1, Float64}([i for i = 0:dy:ysize+dy])
    # P nodes
    "horizontal coordinates of p grid points [m]"
    xp = SVector{Nx1, Float64}([j for j = -dx/2:dx:xsize+dx/2])
    "vertical coordinates of p grid points [m]"
    yp = SVector{Ny1, Float64}([i for i = -dy/2:dy:ysize+dy/2])
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
    "horizontal coordinates of marker grid/launch anchor points [m]"
    xxm = SVector{Nxm, Float64}([j for j = dxm/2:dxm:xsize-dxm/2])
    "vertical coordinates of marker grid/launch anchor points [m]"
    yym = SVector{Nym, Float64}([i for i = dym/2:dym:ysize-dym/2])
    "initialization distance of nearest marker to launch anchor point [m]"
    mdis_init = 1.0e30
    "number of markers at start"
    start_marknum::Int64 = Nxm * Nym
    # marker grid min/max assignables indices
    "minimum assignable marker grid index in x direction"
    jmin_m::Int64 = 1
    "minimum assignable marker grid index in y direction"
    imin_m::Int64 = 1
    "maximum assignable marker grid index in x direction"
    jmax_m::Int64 = Nxm - 1
    "maximum assignable marker grid index in y direction"
    imax_m::Int64 = Nym - 1
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
    # Runge-Kutta integration parameters
    "bⱼ Butcher coefficients for RK4"
    brk4::SVector{4, Float64} = [1/6, 2/6, 2/6, 1/6]
    "cⱼ Butcher coefficients for RK4"
    crk4::SVector{3, Float64} = [0.5, 0.5, 1.0] 
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
Set up and initialize dynamic simulation parameters.

$(SIGNATURES)

# Details

    - sp: static simulation parameters

# Returns

    - timestep: simulation starting time step count
    - dt: simulation initial computational time step [s]
    - timesum: simulation starting time [s]
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
     timestep::Int64 = start_step
     # computational timestep (current), init to dtelastic [s]
     dt::Float64 = dtelastic
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
 
function setup_dynamic_simulation_parameters()
     # timestep counter (current), init to startstep
     timestep::Int64 = start_step
     # computational timestep (current), init to dtelastic [s]
     dt::Float64 = dtelastic
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

    - sp: static simulation parameters
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
    - RX : etafluid/kphi ratio at Vx nodes [m^2]
    - qxD : qx-darcy flux at Vx nodes [m/s]
    - gx : gx-gravity at Vx nodes [m/s^2]
    - RHOY : density at Vx nodes [kg/m^3]
    - RHOFY : fluid density at Vx nodes [kg/m^3]
    - KY : thermal conductivity at Vx nodes [W/m/K]
    - PHIY : porosity at Vx nodes
    - vy : solid vy-velocity at Vx nodes [m/s]
    - vyf : fluid vy-velocity at Vx nodes [m/s]
    - RY : etafluid/kphi ratio at Vx nodes [m^2]
    - qyD : qy-Darcy flux at Vx nodes [m/s]
    - gy : gy-gravity at Vx nodes [m/s^2]
    - RHO : density at P nodes [kg/m^3]
    - RHOCP : volumetric heat capacity at P nodes [J/m^3/K]
    - ALPHA : thermal expansion at P nodes [J/m^3/K]
    - ALPHAF : fluid thermal expansion at P nodes [J/m^3/K]
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
    - BETTAPHI : bulk compresibility at P nodes [Pa*s]
    - PHI : porosity at P nodes
    - APHI : Dlnat P nodes [(1-ϕ)/ϕ]/Dt
    - FI : gravity potential at P nodes [J/kg]
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
    # current temperature [K]
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

# function setup_staggered_grid_properties(; randomized=false)
#     # basic nodes
#     # viscoplastic viscosity [Pa*s]
#     ETA = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
#     # viscous viscosity [Pa*s]
#     ETA0 = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
#     # shear modulus [Pa]
#     GGG = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
#     # epsilonxy [1/s]
#     EXY = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
#     # σxy [1/s]
#     SXY = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
#     # σ₀xy [1/s]
#     SXY0 = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
#     # rotation rate [1/s]
#     wyx = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
#     # compressive strength [Pa]
#     COH = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
#     # tensile strength [Pa]
#     TEN = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
#     # friction
#     FRI = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
#     # plastic yielding node property
#     YNY = randomized ? rand(Bool, Ny, Nx) : zeros(Bool, Ny, Nx)
#     # Vx nodes
#     # density [kg/m^3]
#     RHOX = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # fluid density [kg/m^3]
#     RHOFX = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # thermal conductivity [W/m/K]
#     KX = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # porosity
#     PHIX = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # solid vx-velocity [m/s]
#     vx = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # fluid vx-velocity [m/s]
#     vxf = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # etafluid/kphi ratio [m^2]
#     RX = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # qx-darcy flux [m/s]
#     qxD = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # gx-gravity [m/s^2]
#     gx = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # Vy nodes
#     # density [kg/m^3]
#     RHOY = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # fluid density [kg/m^3]
#     RHOFY = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # thermal conductivity [W/m/K]
#     KY = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # porosity
#     PHIY = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # solid vy-velocity [m/s]
#     vy = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # fluid vy-velocity [m/s]
#     vyf = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # etafluid/kphi ratio [m^2]
#     RY = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # qy-Darcy flux [m/s]
#     qyD = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # gy-gravity [m/s^2]
#     gy = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # P nodes
#     # density [kg/m^3]
#     RHO = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # volumetric heat capacity [J/m^3/K]
#     RHOCP = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # thermal expansion [J/m^3/K]
#     ALPHA = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # fluid thermal expansion [J/m^3/K]
#     ALPHAF = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # radioactive heating [W/m^3]
#     HR = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # adiabatic heating [W/m^3]
#     HA = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # shear heating [W/m^3]
#     HS = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # viscosity [Pa*s]
#     ETAP = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # shear modulus [Pa]
#     GGGP = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # EPSILONxx [1/s]
#     EXX = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # σ′xx [1/s]
#     SXX = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # σ₀′ (SIGMA0'xx) [1/s]
#     SXX0 = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # current temperature [K]
#     tk1 = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # next temperature [K]
#     tk2 = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # solid vx in pressure nodes [m/s]
#     vxp = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # solid vy in pressure nodes [m/s]
#     vyp = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # fluid vx in pressure nodes [m/s]
#     vxpf = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # fluid vy in pressure nodes [m/s]
#     vypf = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # total pressure [Pa]
#     pr = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # fluid pressure [Pa]
#     pf = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # solid pressure [Pa]
#     ps = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # previous total pressure [Pa]
#     pr0 = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # previous fluid pressure [Pa]
#     pf0 = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # previous solid pressure [Pa]
#     ps0 = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # bulk viscosity [Pa*s]
#     ETAPHI = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # bulk compresibility [Pa*s]
#     BETTAPHI = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # porosity
#     PHI = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # Dln[(1-ϕ)/ϕ]/Dt
#     APHI = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # gravity potential [J/kg]
#     FI = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     return (
#         ETA,
#         ETA0,
#         GGG,
#         EXY,
#         SXY,
#         SXY0,
#         wyx,
#         COH,
#         TEN,
#         FRI,
#         YNY,
#         RHOX,
#         RHOFX,
#         KX,
#         PHIX,
#         vx,
#         vxf,
#         RX,
#         qxD,
#         gx,
#         RHOY,
#         RHOFY,
#         KY,
#         PHIY,
#         vy,
#         vyf,
#         RY,
#         qyD,
#         gy,
#         RHO,
#         RHOCP,
#         ALPHA,
#         ALPHAF,
#         HR,
#         HA,
#         HS,
#         ETAP,
#         GGGP,
#         EXX,
#         SXX,
#         SXX0,       
#         tk1,
#         tk2,
#         vxp,
#         vyp,
#         vxpf,
#         vypf,
#         pr,
#         pf,
#         ps,
#         pr0,
#         pf0,
#         ps0,
#         ETAPHI,
#         BETTAPHI,
#         PHI,
#         APHI,
#         FI
#     )
# end # function setup_staggered_grid_properties()

"""
Set up additional helper staggered grid properties to facilitate computations.

$(SIGNATURES)

# Details

    - sp: static simulation parameters
    - randomized: fill in random values for grid properties instead of zeros

# Returns

    - ETA5: plastic iterations viscoplastic viscosity at basic nodes [Pa⋅s]
    - ETA00: previous viscous viscosity at basic nodes [Pa⋅s]
    - YNY5: plastic iterations plastic yielding node property at basic nodes
    - YNY00: previous plastic yielding node property at basic nodes
    - DSXY: stress change Δσxy at basic nodes [Pa]
    - ETAcomp: computational viscosity at basic nodes
    - SXYcomp: computational previous XY stress at basic nodes
    - dRHOXdx: total density gradient in x direction at Vx nodes
    - dRHOXdy: total density gradient in y direction at Vx nodes
    - dRHOYdx: total density gradient in x direction at Vy nodes
    - dRHOYdy: total density gradient in y direction at Vy nodes
    - ETAPcomp: computational viscosity at P nodes
    - SXXcomp: computational previous XX stress at P nodes
    - SYYcomp: computational previous YY stress at P nodes
    - EII :second strain rate invariant at P nodes [1/s]
    - SII :second stress invariant at P nodes [Pa]
    - DSXX :stress change Δσ′xx at P nodes [Pa]
    - tk0: previous temperature at P nodes [K]
"""
function setup_staggered_grid_properties_helpers(sp; randomized=false)
    @unpack Nx, Ny, Nx1, Ny1 = sp
    # basic nodes
    # plastic iterations viscoplastic viscosity at basic nodes [Pa⋅s]
    ETA5 = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
    # previous viscous viscosity at basic nodes [Pa⋅s]
    ETA00 = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
    # plastic iterations plastic yielding node property at basic nodes
    YNY5 = randomized ? rand(Bool, Ny, Nx) : zeros(Bool, Ny, Nx)
    # previous plastic yielding node property at basic nodes
    YNY00 = randomized ? rand(Bool, Ny, Nx) : zeros(Bool, Ny, Nx)
    # inverse viscoplastic viscosity at yielding basic nodes [1/(Pa⋅s)]
    YNY_inv_ETA = randomized ? rand(Bool, Ny, Nx) : zeros(Bool, Ny, Nx)
    # stress change Δσxy at basic nodes [Pa]
    DSXY = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
    # computational viscosity at basic nodes
    ETAcomp = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
    # computational previous xy stress at basic nodes
    SXYcomp = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
    # Vx nodes
    # total density gradient in x direction at Vx nodes
    dRHOXdx = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # total density gradient in y direction at Vx nodes
    dRHOXdy = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # Vy nodes
    # total density gradient in x direction at Vy nodes
    dRHOYdx = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # total density gradient in y direction at Vy nodes
    dRHOYdy = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # P nodes
    # computational viscosity at P nodes
    ETAPcomp = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # computational previous xx stress at P nodes
    SXXcomp = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # computational previous yy stress at P nodes
    SYYcomp = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # second strain rate invariant at P nodes [1/s]
    EII = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # second stress invariant at P nodes [Pa]
    SII = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # stress change Δσ′xx at P nodes [Pa]
    DSXX = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    # previous temperature at P nodes [K]
    tk0 = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
    return (
        ETA5,
        ETA00,
        YNY5,
        YNY00,
        YNY_inv_ETA,
        DSXY,
        ETAcomp,
        SXYcomp,
        dRHOXdx,
        dRHOXdy,
        dRHOYdx,
        dRHOYdy,
        ETAPcomp,
        SXXcomp,
        SYYcomp,
        EII,
        SII,
        DSXX,
        tk0
    )
end # function setup_staggered_grid_properties_helpers()

# function setup_staggered_grid_properties_helpers(;randomized=false)
#     # basic nodes
#     # plastic iterations viscoplastic viscosity at basic nodes [Pa⋅s]
#     ETA5 = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
#     # previous viscous viscosity at basic nodes [Pa⋅s]
#     ETA00 = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
#     # plastic iterations plastic yielding node property at basic nodes
#     YNY5 = randomized ? rand(Bool, Ny, Nx) : zeros(Bool, Ny, Nx)
#     # previous plastic yielding node property at basic nodes
#     YNY00 = randomized ? rand(Bool, Ny, Nx) : zeros(Bool, Ny, Nx)
#     # inverse viscoplastic viscosity at yielding basic nodes [1/(Pa⋅s)]
#     YNY_inv_ETA = randomized ? rand(Bool, Ny, Nx) : zeros(Bool, Ny, Nx)
#     # stress change Δσxy at basic nodes [Pa]
#     DSXY = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
#     # computational viscosity at basic nodes
#     ETAcomp = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
#     # computational previous xy stress at basic nodes
#     SXYcomp = randomized ? rand(Ny, Nx) : zeros(Ny, Nx)
#     # Vx nodes
#     # total density gradient in x direction at Vx nodes
#     dRHOXdx = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # total density gradient in y direction at Vx nodes
#     dRHOXdy = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # Vy nodes
#     # total density gradient in x direction at Vy nodes
#     dRHOYdx = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # total density gradient in y direction at Vy nodes
#     dRHOYdy = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # P nodes
#     # computational viscosity at P nodes
#     ETAPcomp = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # computational previous xx stress at P nodes
#     SXXcomp = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # computational previous yy stress at P nodes
#     SYYcomp = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # second strain rate invariant at P nodes [1/s]
#     EII = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # second stress invariant at P nodes [Pa]
#     SII = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # stress change Δσ′xx at P nodes [Pa]
#     DSXX = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     # previous temperature at P nodes [K]
#     tk0 = randomized ? rand(Ny1, Nx1) : zeros(Ny1, Nx1)
#     return (
#         ETA5,
#         ETA00,
#         YNY5,
#         YNY00,
#         YNY_inv_ETA,
#         DSXY,
#         ETAcomp,
#         SXYcomp,
#         dRHOXdx,
#         dRHOXdy,
#         dRHOYdx,
#         dRHOYdy,
#         ETAPcomp,
#         SXXcomp,
#         SYYcomp,
#         EII,
#         SII,
#         DSXX,
#         tk0
#     )
# end # function setup_staggered_grid_properties_helpers()


"""
Set up geodesic and physical properties of the set of markers.

$(SIGNATURES)

# Details

    - sp: static simulation parameters
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
function setup_marker_properties(sp; randomized=false)
    @unpack start_marknum, dx, dy, xsize, ysize = sp
    marknum = start_marknum
    # horizontal marker coordinate [m]
    xm = randomized ? rand(-dx:0.1:xsize+dx, marknum) : zeros(marknum)
    # vertical marker coordinate [m]
    ym = randomized ? rand(-dy:0.1:ysize+dy, marknum) : zeros(marknum)
    # marker material type
    tm = randomized ? rand(1:3, marknum) : zeros(Int, marknum)
    # marker temperature [K]
    tkm = randomized ? rand(273:300, marknum) : zeros(marknum)
    # marker σ′xx [Pa]
    sxxm = randomized ? rand(marknum) : zeros(marknum)
    # marker σxy [Pa]
    sxym = randomized ? rand(marknum) : zeros(marknum)
    # marker viscoplastic viscosity [Pa]
    etavpm = randomized ? rand(marknum) : zeros(marknum)
    # marker porosity
    phim = randomized ? rand(marknum) : zeros(marknum)
    return xm, ym, tm, tkm, sxxm, sxym, etavpm, phim
end # function setup_marker_properties()

# function setup_marker_properties(marknum; randomized=false)
#     # horizontal marker coordinate [m]
#     xm = randomized ? rand(-dx:0.1:xsize+dx, marknum) : zeros(marknum)
#     # vertical marker coordinate [m]
#     ym = randomized ? rand(-dy:0.1:ysize+dy, marknum) : zeros(marknum)
#     # marker material type
#     tm = randomized ? rand(1:3, marknum) : zeros(Int, marknum)
#     # marker temperature [K]
#     tkm = randomized ? rand(273:300, marknum) : zeros(marknum)
#     # marker σ′xx [Pa]
#     sxxm = randomized ? rand(marknum) : zeros(marknum)
#     # marker σxy [Pa]
#     sxym = randomized ? rand(marknum) : zeros(marknum)
#     # marker viscoplastic viscosity [Pa]
#     etavpm = randomized ? rand(marknum) : zeros(marknum)
#     # marker porosity
#     phim = randomized ? rand(marknum) : zeros(marknum)
#     return xm, ym, tm, tkm, sxxm, sxym, etavpm, phim
# end # function setup_marker_properties()


"""
Set up additional helper marker properties to facility comptuations.

$(SIGNATURES)

# Details

    - sp: static simulation parameters
    - randomized: fill in random values for grid properties instead of zeros

# Returns

    - rhototalm: total density of markers
    - rhocptotalm : total volumetric heat capacity of markers
    - etasolidcur: solid viscosity of markers
    - etafluidcur: fluid viscosity of markers
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
function setup_marker_properties_helpers(sp; randomized=false)
    @unpack start_marknum = sp
    marknum = start_marknum
    # marker total density
    rhototalm = randomized ? rand(marknum) : zeros(marknum)
    # marker total volumetric heat capacity
    rhocptotalm = randomized ? rand(marknum) : zeros(marknum)
    # marker solid viscosity
    etasolidcur = randomized ? rand(marknum) : zeros(marknum)
    # marker fluid viscosity
    etafluidcur = randomized ? rand(marknum) : zeros(marknum)
    # marker total viscosity
    etatotalm = randomized ? rand(marknum) : zeros(marknum)
    # marker total radiogenic heat production
    hrtotalm = randomized ? rand(marknum) : zeros(marknum)
    # marker total thermal conductivity
    ktotalm = randomized ? rand(marknum) : zeros(marknum)
    # marker total thermal energy
    tkm_rhocptotalm = randomized ? rand(marknum) : zeros(marknum)
    # marker fluid viscosity over permeability
    etafluidcur_inv_kphim = randomized ? rand(marknum) : zeros(marknum)
    # marker inverse of total shear modulus
    inv_gggtotalm = randomized ? rand(marknum) : zeros(marknum)
    # marker total friction coefficient
    fricttotalm = randomized ? rand(marknum) : zeros(marknum)
    # marker total compressive strength
    cohestotalm = randomized ? rand(marknum) : zeros(marknum)
    # marker total tensile strength
    tenstotalm = randomized ? rand(marknum) : zeros(marknum)
    # marker fluid density
    rhofluidcur = randomized ? rand(marknum) : zeros(marknum)
    # marker solid thermal expansion coefficient
    alphasolidcur = randomized ? rand(marknum) : zeros(marknum)
    # marker fluid thermal expansion coefficient
    alphafluidcur = randomized ? rand(marknum) : zeros(marknum)
    return (
        rhototalm,
        rhocptotalm,
        etasolidcur,
        etafluidcur,
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

# function setup_marker_properties_helpers(marknum; randomized=false)
#     # marker total density
#     rhototalm = randomized ? rand(marknum) : zeros(marknum)
#     # marker total volumetric heat capacity
#     rhocptotalm = randomized ? rand(marknum) : zeros(marknum)
#     # marker solid viscosity
#     etasolidcur = randomized ? rand(marknum) : zeros(marknum)
#     # marker fluid viscosity
#     etafluidcur = randomized ? rand(marknum) : zeros(marknum)
#     # marker total viscosity
#     etatotalm = randomized ? rand(marknum) : zeros(marknum)
#     # marker total radiogenic heat production
#     hrtotalm = randomized ? rand(marknum) : zeros(marknum)
#     # marker total thermal conductivity
#     ktotalm = randomized ? rand(marknum) : zeros(marknum)
#     # marker total thermal energy
#     tkm_rhocptotalm = randomized ? rand(marknum) : zeros(marknum)
#     # marker fluid viscosity over permeability
#     etafluidcur_inv_kphim = randomized ? rand(marknum) : zeros(marknum)
#     # marker inverse of total shear modulus
#     inv_gggtotalm = randomized ? rand(marknum) : zeros(marknum)
#     # marker total friction coefficient
#     fricttotalm = randomized ? rand(marknum) : zeros(marknum)
#     # marker total compressive strength
#     cohestotalm = randomized ? rand(marknum) : zeros(marknum)
#     # marker total tensile strength
#     tenstotalm = randomized ? rand(marknum) : zeros(marknum)
#     # marker fluid density
#     rhofluidcur = randomized ? rand(marknum) : zeros(marknum)
#     # marker solid thermal expansion coefficient
#     alphasolidcur = randomized ? rand(marknum) : zeros(marknum)
#     # marker fluid thermal expansion coefficient
#     alphafluidcur = randomized ? rand(marknum) : zeros(marknum)
#     return (
#         rhototalm,
#         rhocptotalm,
#         etasolidcur,
#         etafluidcur,
#         etatotalm,
#         hrtotalm,
#         ktotalm,
#         tkm_rhocptotalm,
#         etafluidcur_inv_kphim,
#         inv_gggtotalm,
#         fricttotalm,
#         cohestotalm,
#         tenstotalm,
#         rhofluidcur,
#         alphasolidcur,
#         alphafluidcur
#     )
# end # function setup_marker_properties_helpers()

"""
Set up additional marker geometry helpers to facilitate marker handling.

$(SIGNATURES)

# Details

    - sp: static simulation parameters

# Returns

    - mdis: minimum distance of marker launch anchor points to nearest marker
    - mnum: number of marker nearest to marker launch anchor positions
    # - mtyp: type of marker nearest to marker launch anchor positions
    # - mpor: porosity of marker nearest to marker launch anchor positions
"""
function setup_marker_geometry_helpers(sp)
    @unpack Nxm, Nym, mdis_init = sp
    mdis = fill(mdis_init, Nym, Nxm)
    mnum = zeros(Int, Nym, Nxm)
    # mtyp = zeros(Int, Nym, Nxm)
    # mpor = zeros(Nym, Nxm)
    return mdis, mnum#, mtyp, mpor 
end

function setup_marker_geometry_helpers()
    mdis = fill(mdis_init, Nym, Nxm)
    mnum = zeros(Int, Nym, Nxm)
    # mtyp = zeros(Int, Nym, Nxm)
    # mpor = zeros(Nym, Nxm)
    return mdis, mnum#, mtyp, mpor 
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
    - etafluidcur: fluid viscosity of markers
    - tkm: temperature of markers 
    - inv_gggtotalm: inverse of total shear modulus of markers
    - fricttotalm: total friction coefficient of markers
    - cohestotalm: total compressive strength of markers
    - tenstotalm: total tensile strength of markers
    - rhofluidcur: fluid density of markers
    - alphasolidcur: solid thermal expansion coefficient of markers
    - alphafluidcur: fluid thermal expansion coefficient of markers
    - sp: static simulation parameters
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
    etafluidcur,
    tkm,
    inv_gggtotalm,
    fricttotalm,
    cohestotalm,
    tenstotalm,
    rhofluidcur,
    alphasolidcur,
    alphafluidcur,
    sp;
    randomized=true
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
        xm[m] = dxm/2 + (jm-1) * dxm 
        ym[m] = dym/2 + (im-1) * dym 
        # random marker position within cell
        if randomized
            xm[m] += (rand()-0.5) * dxm
            ym[m] += (rand()-0.5) * dym
        end
        # primary marker properties 
        rmark = distance(xm[m], ym[m], xcenter, ycenter)
        if rmark < rplanet
            # planet
            tm[m] = ifelse(rmark>rcrust, 2, 1)
            # porosity
            phim[m] = phim0
            if randomized
                phim[m] += phim0 * (rand()-0.5)
            end
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
    - tm: type of markers
    - tkm: temperature of markers
    - rhototalm: total density of markers
    - rhocptotalm: total volumetric heat capacity of markers
    - etasolidcur: solid viscosity of markers
    - etafluidcur: fluid viscosity of markers
    - etatotalm: total viscosity of markers
    - hrtotalm: total radiogenic heat production of markers
    - ktotalm: total thermal conductivity of markers
    - tkm_rhocptotalm: total thermal energy of markers
    - etafluidcur_inv_kphim: (fluid viscosity)/permeability of markers
    - phim: porosity of markers
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
    sp
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
    - sp: static simulation parameters

# Returns

    -nothing
"""
function update_marker_viscosity!(
    m, xm, ym, tm, tkm, etatotalm, etavpm, YNY, YNY_inv_ETA, sp)
@timeit to "update_marker_viscosity!" begin
    @unpack x, y, dx, dy, jmin_basic, jmax_basic, imin_basic, imax_basic = sp
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
    if tm[m] < 3
        # rocks: update etatotalm[m] based on current marker temperature
        etatotalm[m] = etatotal_rocks(tkm[m], tm[m], sp)
    # else
        # air: constant etatotalm[m] as initialized
        # pass
    end
    if YNY[i,j] || YNY[i+1,j] || YNY[i,j+1] || YNY[i+1,j+1]
        interpolate_to_marker!(m, i, j, weights, etavpm, YNY_inv_ETA)
        etavpm[m] = inv(etavpm[m])
        etavpm[m] = ifelse(etavpm[m]>etatotalm[m], etatotalm[m], etavpm[m])
    else
        etavpm[m] = etatotalm[m]
    end
end # @timeit to "update_marker_viscosity!"
    return nothing
end


"""
Set up properties to be interpolated from markers to staggered grid.

$(SIGNATURES)

# Details

    - sp: static simulation parameters

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
    - WTPSUM: interpolation weights at P nodes
"""
function setup_interpolated_properties(sp)
    @unpack Nx, Ny, Nx1, Ny1 = sp
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
        WTPSUM
    )
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
Get a 4-vector of values from grid in row-major order.

$(SIGNATURES)

# Details

    - i: top left grid node column index
    - j: top left grid node row index
    - grid: data from which to build Vector

# Returns

    -grid_vector: 4-vector of values
    [grid[i, j], grid[i+1, j], grid[i, j+1], grid[i+1, j+1]]
"""
function grid_vector(i, j, grid)
    return @inbounds @SVector [
        grid[i, j], grid[i+1, j], grid[i, j+1], grid[i+1, j+1]
    ]
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
            + ((ksolid*(3.0*phi-2.0) + kfluid*(1.0-3.0*phi))^2)/16.0
        )
        -0.25 * (ksolid*(3.0*phi-2.0) + kfluid*(1.0-3.0*phi))
    )
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
Compute total rocky marker viscosity based on temperature and material type.

$(SIGNATURES)

# Details

    - tkmm: marker temperature [K]
    - tmm: marker type [1, 2]
    - sp: simulation parameters

# Returns
    
    - etatotal: rocky marker temperature-dependent total viscosity 
"""
function etatotal_rocks(tkmm, tmm, sp)
    @unpack tmsilicate,
        tmiron,
        etasolidm,
        etasolidmm,
        etafluidm,
        etafluidmm,
        etamin = sp
    @inbounds etasolidcur = ifelse(
        tkmm>tmsilicate, etasolidmm[tmm], etasolidm[tmm])
    @inbounds etafluidcur = ifelse(
        tkmm>tmiron, etafluidmm[tmm], etafluidm[tmm])
    return max(etamin, etasolidcur, etafluidcur)
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
# @timeit to "fix_distances" begin
    @inbounds begin
        i, j = fix(x, y, x_axis, y_axis, dx, dy, jmin, jmax, imin, imax)
        dxmj = x - x_axis[j]
        dymi = y - y_axis[i]
    end # @inbounds
    return i, j, dxmj, dymi
# end # @timeit to "fix_distances"
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
# @timeit to "fix" begin
    j = unsafe_trunc(Int, (x-x_axis[1])/dx) + 1
    i = unsafe_trunc(Int, (y-y_axis[1])/dy) + 1
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
# end # @timeit to "fix"
end # function fix

function fix(x, y, x_axis, y_axis, jmin, jmax, imin, imax)
    j = unsafe_trunc(Int, (x-x_axis[1])/dx) + 1
    i = unsafe_trunc(Int, (y-y_axis[1])/dy) + 1
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
end

"""
Reduce a 3D (i, j, k) along its third (k) axis by addition and write the result
into (i, j, 1) without reallocating the array's memory.

$(SIGNATURES)

# Details

    - A: 3D array [i, j, k]

# Returns

    - nothing
"""
function reduce_add_3darray!(A)
    for k in 2:size(A, 3), j in 1:size(A, 2), i in 1:size(A, 1)
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
# @timeit to "interpolate_add_to_grid!" begin
    grid[i, j, threadid()] += property * weights[1]
    grid[i+1, j, threadid()] += property * weights[2]
    grid[i, j+1, threadid()] += property * weights[3]
    grid[i+1, j+1, threadid()] += property * weights[4]
# end # @timeit to "interpolate_add_to_grid!"
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
# @timeit to "interpolate_to_marker!()" begin
    @inbounds marker_property[m] = dot(grid_vector(i, j, grid), weights)
# end # @timeit to "interpolate_to_marker!()"
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
# @timeit to "interpolate_add_to_marker!()" begin
    @inbounds marker_property[m] += dot(grid_vector(i, j, grid), weights)
# end # @timeit to "interpolate_add_to_marker!()"
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
    - sp: static simulation parameters

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
    WTSUM,
    sp
)
    @unpack dx, dy, x, y, jmin_basic, jmax_basic, imin_basic, imax_basic = sp
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
    interpolate_add_to_grid!(i, j, weights, etatotalm[m], ETA0SUM)
    interpolate_add_to_grid!(i, j, weights, etavpm[m], ETASUM)
    interpolate_add_to_grid!(i, j, weights, inv_gggtotalm[m], GGGSUM)
    interpolate_add_to_grid!(i, j, weights, sxym[m], SXYSUM)
    interpolate_add_to_grid!(i, j, weights, cohestotalm[m], COHSUM)
    interpolate_add_to_grid!(i, j, weights, tenstotalm[m], TENSUM)
    interpolate_add_to_grid!(i, j, weights, fricttotalm[m], FRISUM)
    interpolate_add_to_grid!(i, j, weights, 1.0, WTSUM)
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
    - sp: static simulation parameters

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
    WTXSUM,
    sp
)
    @unpack dx, dy, xvx, yvx, jmin_vx, jmax_vx, imin_vx, imax_vx = sp
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
    interpolate_add_to_grid!(i, j, weights, rhototalm[m], RHOXSUM)
    interpolate_add_to_grid!(i, j, weights, rhofluidcur[m], RHOFXSUM)
    interpolate_add_to_grid!(i, j, weights, ktotalm[m], KXSUM)
    interpolate_add_to_grid!(i, j, weights, phim[m], PHIXSUM)
    interpolate_add_to_grid!(i, j, weights, etafluidcur_inv_kphim[m], RXSUM)
    interpolate_add_to_grid!(i, j, weights, 1.0, WTXSUM)
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
    - sp: static simulation parameters

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
    WTYSUM,
    sp
)
    @unpack dx, dy, xvy, yvy, jmin_vy, jmax_vy, imin_vy, imax_vy = sp
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
    interpolate_add_to_grid!(i, j, weights, rhototalm[m], RHOYSUM)
    interpolate_add_to_grid!(i, j, weights, rhofluidcur[m], RHOFYSUM)
    interpolate_add_to_grid!(i, j, weights, ktotalm[m], KYSUM)
    interpolate_add_to_grid!(i, j, weights, phim[m], PHIYSUM)
    interpolate_add_to_grid!(i, j, weights, etafluidcur_inv_kphim[m], RYSUM)
    interpolate_add_to_grid!(i, j, weights, 1.0, WTYSUM)
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
    - sp: static simulation parameters

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
    WTPSUM,
    sp
)
    @unpack dx, dy, xp, yp, jmin_p, jmax_p, imin_p, imax_p = sp
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
    interpolate_add_to_grid!(i, j, weights, inv_gggtotalm[m], GGGPSUM)
    interpolate_add_to_grid!(i, j, weights, sxxm[m], SXXSUM)
    interpolate_add_to_grid!(i, j, weights, rhototalm[m], RHOSUM)
    interpolate_add_to_grid!(i, j, weights, rhocptotalm[m], RHOCPSUM)
    interpolate_add_to_grid!(i, j, weights, alphasolidcur[m], ALPHASUM)
    interpolate_add_to_grid!(i, j, weights, alphafluidcur[m], ALPHAFSUM)
    interpolate_add_to_grid!(i, j, weights, hrtotalm[m], HRSUM)
    interpolate_add_to_grid!(i, j, weights, phim[m], PHISUM)
    interpolate_add_to_grid!(i, j, weights, tkm_rhocptotalm[m], TKSUM)
    interpolate_add_to_grid!(i, j, weights, 1.0, WTPSUM)
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
    - SXX0:σ₀xy XY stress at basic nodes
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
    @views @. dRHOXdx[:, 2:Nx] = (RHOX[:, 3:Nx1]-RHOX[:, 1:Nx1-2]) / 2 / dx
    @views @. dRHOXdy[2:Ny, :] = (RHOX[3:Ny1, :]-RHOX[1:Ny1-2, :]) / 2 / dy
    @views @. dRHOYdx[:, 2:Nx] = (RHOY[:, 3:Nx1]-RHOY[:, 1:Nx1-2]) / 2 / dx
    @views @. dRHOYdy[2:Ny, :] = (RHOY[3:Ny1, :]-RHOY[1:Ny1-2, :]) / 2 / dy
    return nothing
end # @timeit to "get_viscosities_stresses_density_gradients!()"
end # function get_viscosities_stresses_density_gradients!


"""
Set up hydromechanical linear system of equations structures.

$(SIGNATURES)

# Details

    - sp: static simulation parameters

# Returns 

    _ L: hydromechanical linear system of equations: LHS coefficient matrix    
    - R: hydromechanical linear system of equations: RHS vector
    - S: hydromechanical linear system of equations: solution vector
"""
function setup_hydromechanical_lse(sp)
    @unpack Nx1, Ny1 = sp
    L = ExtendableSparseMatrix(Ny1*Nx1*6, Ny1*Nx1*6)
    R = zeros(Ny1*Nx1*6)
    S = Vector{Float64}(undef, Ny1*Nx1*6)
    return L, R, S
end


"""
Set up thermal linear system of equations structures.

$(SIGNATURES)

# Details

    - sp: static simulation parameters

# Returns 

    _ LT: thermal linear system of equations: LHS coefficient matrix
    - RT: thermal linear system of equations: RHS vector
    - ST: thermal linear system of equations: solution vector
"""
function setup_thermal_lse(sp)
    @unpack Nx1, Ny1 = sp
    LT = ExtendableSparseMatrix(Ny1*Nx1, Ny1*Nx1)
    RT = zeros(Ny1*Nx1)
    ST = Vector{Float64}(undef, Ny1*Nx1)
    return LT, RT, ST
end

"""
Set up gravitational linear system of equations structures.

$(SIGNATURES)

# Details

    - sp: static simulation parameters

# Returns 

    _ LP: gravitational linear system of equations: LHS coefficient matrix
    - RP: gravitational linear system of equations: RHS vector
    - SP: gravitational linear system of equations: solution vector
"""
function setup_gravitational_lse(sp)
    @unpack Nx1, Ny1 = sp
    LP = ExtendableSparseMatrix(Ny1*Nx1, Ny1*Nx1)
    RP = zeros(Ny1*Nx1)
    SP = Vector{Float64}(undef, Ny1*Nx1)
    return LP, RP, SP
end


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
    - ETAP: viscosity at P nodes
    - ETAPHI: bulk viscosity at P nodes
    - PHI: porosity at P nodes
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
    - sp: static simulation parameters

## Out

    - APHI: porosity coefficient at P nodes

# Returns

    - aphimax: maximum absolute porosity coefficient
"""
function compute_Aϕ!(APHI, ETAPHI, BETTAPHI, PHI, pr, pf, pr0, pf0, dt, sp)
@timeit to "compute_Aϕ!()" begin
    @unpack Nx, Ny = sp
    # reset APHI
    APHI .= 0.0
    @views @. APHI[2:Ny, 2:Nx] = (
        ((pr[2:Ny, 2:Nx]-pf[2:Ny, 2:Nx])/ETAPHI[2:Ny, 2:Nx]
        + (
            (pr[2:Ny, 2:Nx]-pr0[2:Ny, 2:Nx])-(pf[2:Ny, 2:Nx]-pf0[2:Ny, 2:Nx])
        )/dt*BETTAPHI[2:Ny, 2:Nx]) / (1-PHI[2:Ny, 2:Nx]) / PHI[2:Ny, 2:Nx]
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
    - YNY_inv_ETA: inverse of plastic viscosity at yielding basic nodes
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
    YNY_inv_ETA,
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
    @views @. YNY_inv_ETA = YNY / ETA
    return dt
end # @timeit to "finalize_plastic_iteration_pass!()"
    end # function finalize_plastic_iteration_pass


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
    - dtm: displacement time step
    - marknum: total number of markers in use
    - sp: static simulation parameters  

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
    dtm,
    marknum,
    sp
)
@timeit to "apply_subgrid_stress_diffusion!" begin
    @unpack Nx,
        Ny,
        dx,
        dy,
        x,
        y,
        jmin_basic,
        jmax_basic,
        imin_basic,
        imax_basic,
        xp,
        yp,
        jmin_p,
        jmax_p,
        imin_p,
        imax_p,
        dsubgrids = sp
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
    @threads for m=1:1:marknum
        i_p, j_p, weights_p = fix_weights(
            xm[m], ym[m], xp, yp, dx, dy, jmin_p, jmax_p, imin_p, imax_p)
        i_basic, j_basic, weights_basic = fix_weights(
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
        δσxxm₀ = sxxm[m] - dot(grid_vector(i_p, j_p, SXX0), weights_p)
        # time-relax σ′xx difference
        δσxxm₀ *= (exp(-dsubgrids*dtm/(etam[tm[m]]*inv_gggtotalm[m])) - 1.0) 
        # correct marker stress
        sxxm[m] += δσxxm₀
        # update subgrid diffusion on P nodes
        interpolate_add_to_grid!(i_p, j_p, weights_p, δσxxm₀, SXXSUM)
        interpolate_add_to_grid!(i_p, j_p, weights_p, 1.0, WTPSUM)
        # σ₀xy at basic nodes
        # compute marker-node σxy difference
        δσxy₀ = sxym[m] - dot(
            grid_vector(i_basic, j_basic, SXY0), weights_basic)
        # time-relax σxy difference
        δσxy₀ *= (exp(-dsubgrids*dtm/(etam[tm[m]]*inv_gggtotalm[m])) - 1.0)
        # correct marker stress
        sxym[m] += δσxy₀
        # update subgrid diffusion on basic nodes
        interpolate_add_to_grid!(i_basic, j_basic, weights_basic, δσxy₀, SXYSUM)
        interpolate_add_to_grid!(i_basic, j_basic, weights_basic, 1.0, WTSUM)
    end
    # reduce interpolation arrays
    SXXSUM = reduce(+, SXXSUM, dims=3)
    WTPSUM = reduce(+, WTPSUM, dims=3)
    SXYSUM = reduce(+, SXYSUM, dims=3)
    WTSUM = reduce(+, WTSUM, dims=3)

    # compute DSXXsubgrid and update DSXX at inner P nodes
    @views @. DSXX[2:Ny, 2:Nx][WTPSUM[2:Ny, 2:Nx, 1]>0.0] -= (
        SXXSUM[2:Ny, 2:Nx, 1][WTPSUM[2:Ny, 2:Nx, 1]>0.0] /
        WTPSUM[2:Ny, 2:Nx, 1][WTPSUM[2:Ny, 2:Nx, 1]>0.0]
    )
    # compute DSXYsubgrid and update DSXY at all basic nodes
    @views @. DSXY[WTSUM[:, :, 1]>0.0] -= SXYSUM[:, :, 1][WTSUM[:, :, 1]>0.0] /
        WTSUM[:, :, 1][WTSUM[:, :, 1]>0.0]
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
    - sp: static simulation parameters

# Returns

    - nothing
"""
function update_marker_stress!(xm, ym, sxxm, sxym, DSXX, DSXY, marknum, sp)
@timeit to "update_marker_stress!" begin
    @unpack dx,
        dy,
        x,
        y,
        jmin_basic,
        jmax_basic,
        imin_basic,
        imax_basic,
        xp,
        yp,
        jmin_p,
        jmax_p,
        imin_p,
        imax_p = sp
    @threads for m=1:1:marknum    
        i_p, j_p, weights_p = fix_weights(
            xm[m],
            ym[m],
            xp,
            yp,
            dx,
            dy,
            jmin_p+1,
            jmax_p-1,
            imin_p+1,
            imax_p-1
        )
        i_basic, j_basic, weights_basic = fix_weights(
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
    - sp: static simulation parameters

# Returns

    - nothing
"""
function compute_shear_heating!(
    HS, SXYEXY, ETA, SXY, ETAP, SXX, RX, RY, qxD, qyD, PHI, ETAPHI, pr, pf, sp)
# @timeit to "compute_shear_heating!" begin
    @unpack Nx, Ny, Nx1, Ny1 = sp
    # average SXY⋅EXY
    @inbounds begin
    @views @. SXYEXY[2:Ny, 2:Nx] = 0.25 * (
        SXY[1:Ny-1, 1:Nx-1]^2/ETA[1:Ny-1, 1:Nx-1]
        + SXY[2:Ny, 1:Nx-1]^2/ETA[2:Ny, 1:Nx-1]
        + SXY[1:Ny-1, 2:Nx]^2/ETA[1:Ny-1, 2:Nx]
        + SXY[2:Ny, 2:Nx]^2/ETA[2:Ny, 2:Nx]
    )
    # compute shear heating HS
    @views @. HS[2:Ny, 2:Nx] = (
        SXX[2:Ny, 2:Nx]^2 / ETAP[2:Ny, 2:Nx]
        + SXYEXY[2:Ny, 2:Nx]
        + (
            (pr[2:Ny, 2:Nx]-pf[2:Ny, 2:Nx])^2
            / (1-PHI[2:Ny, 2:Nx])
            / ETAPHI[2:Ny, 2:Nx]
        )
        + 0.5 * (
            RX[2:Ny, 1:Nx-1] * qxD[2:Ny, 1:Nx-1]^2
            + RX[2:Ny, 2:Nx] * qxD[2:Ny, 2:Nx]^2
        )
        + 0.5 * (
            RY[1:Ny-1, 2:Nx] * qyD[1:Ny-1, 2:Nx]^2
            + RY[2:Ny, 2:Nx] * qyD[2:Ny, 2:Nx]^2
        )
    )
    end # @inbounds
# end # @timeit to "compute_shear_heating!" 
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
    - sp: static simulation parameters

# Returns

    - nothing
"""
function compute_adiabatic_heating!(
    HA, tk1, ALPHA, ALPHAF, PHI, vx, vy, vxf, vyf, ps, pf, sp)
# @timeit to "compute_adiabatic_heating!" begin
    @unpack Nx, Ny, Nx1, Ny1, dx, dy = sp
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
                dpsdx = (ps[i,j]-ps[i,j-1]) / dx
            else
                dpsdx = (ps[i,j+1]-ps[i,j]) / dx
            end
            if VYP < 0.0
                dpsdy = (ps[i,j]-ps[i-1,j]) / dy
            else
                dpsdy = (ps[i+1,j]-ps[i,j]) / dy
            end
            dpsdt = VXP*dpsdx + VYP*dpsdy
            # evaluate DPfluid/Dt with upwind differences
            if VXFP > 0.0
                dpfdx = (pf[i,j]-pf[i,j-1]) / dx
            else
                dpfdx = (pf[i,j+1]-pf[i,j]) / dx
            end
            if VYFP > 0.0
                dpfdy = (pf[i,j]-pf[i-1,j]) / dy
            else
                dpfdy = (pf[i+1,j]-pf[i,j]) / dy
            end
            dpfdt = VXFP*dpsdx + VYFP*dpsdy
            # Hₐ = (1-ϕ)Tαˢ⋅DPˢ/Dt + ϕTαᶠ⋅DPᶠ/Dt (eq. 9.23)
            HA[i, j] = (
                (1-PHI[i, j]) * tk1[i, j] * ALPHA[i, j] * dpsdt
                + PHI[i, j] * tk1[i, j] * ALPHAF[i, j] * dpfdt
            )
        end
    end # @inbounds
# end # @timeit to "compute_adiabatic_heating!"
end # function compute_adiabatic_heating!

# RMK: vectorized version below: untested, but slower anyway
# function compute_adiabatic_heating!(
#     HA, VXP, VYP, VXFP, VYFP, ∂ps∂xₗ, ∂ps∂xᵣ, ∂ps∂yₗ, ∂ps∂yᵣ, ∂pf∂xₗ, ∂pf∂xᵣ, ∂pf∂yₗ, ∂pf∂yᵣ, dpsdx,dpsdy, dpfdx, dpfdy, dpsdt, dpfdt, tk1, ALPHA, ALPHAF, PHI, vx, vy, vxf, vyf, ps, pf, sp)
# # @timeit to "compute_adiabatic_heating!" begin
#     @unpack Nx, Ny, Nx1, Ny1, dx, dy = sp
#     # indirect calculation of ∂p∂t
#     # average vx, vy, vxf, vyf
#     @views @. VXP[2:Ny, 2:Nx] = 0.5 * (vx[2:Ny, 2:Nx] + vx[2:Ny, 1:Nx-1])
#     @views @. VYP[2:Ny, 2:Nx] = 0.5 * (vy[2:Ny, 2:Nx] + vy[1:Ny-1, 2:Nx])
#     @views @. VXFP[2:Ny, 2:Nx] = 0.5 * (vxf[2:Ny, 2:Nx] + vxf[2:Ny, 1:Nx-1])
#     @views @. VYFP[2:Ny, 2:Nx] = 0.5 * (vyf[2:Ny, 2:Nx] + vyf[1:Ny-1, 2:Nx])
#     # left and right upwind ∂p∂x, ∂p∂y, ∂pf∂x, ∂pf∂y
#     @views @. ∂ps∂xₗ[2:Ny, 2:Nx] = (ps[2:Ny, 2:Nx]-ps[2:Ny, 1:Nx-1]) / dx
#     @views @. ∂ps∂xᵣ[2:Ny, 2:Nx] = (ps[2:Ny, 3:Nx1]-ps[2:Ny, 2:Nx]) / dx
#     @views @. ∂ps∂yₗ[2:Ny, 2:Nx] = (ps[2:Ny, 2:Nx]-ps[1:Ny-1, 2:Nx]) / dy
#     @views @. ∂ps∂yᵣ[2:Ny, 2:Nx] = (ps[3:Ny1, 2:Nx]-ps[2:Ny, 2:Nx]) / dy
#     @views @. ∂pf∂xₗ[2:Ny, 2:Nx] = (pf[2:Ny, 2:Nx]-pf[2:Ny, 1:Nx-1]) / dx
#     @views @. ∂pf∂xᵣ[2:Ny, 2:Nx] = (pf[2:Ny, 3:Nx1]-pf[2:Ny, 2:Nx]) / dx
#     @views @. ∂pf∂yₗ[2:Ny, 2:Nx] = (pf[2:Ny, 2:Nx]-pf[1:Ny-1, 2:Nx]) / dy
#     @views @. ∂pf∂yᵣ[2:Ny, 2:Nx] = (pf[3:Ny1, 2:Nx]-pf[2:Ny, 2:Nx]) / dy
#     # evaluate DPsolid/Dt, DPfluid/Dt
#     @views @. dpsdx = ∂ps∂xₗ
#     @views @. dpsdx[2:Ny, 2:Nx][VXP[2:Ny, 2:Nx]>=0] = ∂ps∂xᵣ[2:Ny, 2:Nx][VXP[2:Ny, 2:Nx]>=0]
#     @views @. dpsdy = ∂ps∂yₗ
#     @views @. dpsdy[2:Ny, 2:Nx][VYP[2:Ny, 2:Nx]>=0] = ∂ps∂yᵣ[2:Ny, 2:Nx][VYP[2:Ny, 2:Nx]>=0]
#     @views @. dpsdt = VXP*dpsdx + VYP*dpsdy
#     @views @. dpfdx = ∂pf∂xₗ
#     @views @. dpfdx[2:Ny, 2:Nx][VXFP[2:Ny, 2:Nx]<=0] = ∂pf∂xᵣ[2:Ny, 2:Nx][VXFP[2:Ny, 2:Nx]<=0]
#     @views @. dpfdy = ∂pf∂yₗ
#     @views @. dpfdy[2:Ny, 2:Nx][VYFP[2:Ny, 2:Nx]<=0] = ∂pf∂yᵣ[2:Ny, 2:Nx][VYFP[2:Ny, 2:Nx]<=0]
#     @views @. dpfdt = VXFP*dpfdx + VYFP*dpfdy
#     # compute adiabatic heating HA
#     @views @. HA[2:Ny, 2:Nx] = (
#         (1-PHI[2:Ny, 2:Nx])*tk1[2:Ny, 2:Nx]*ALPHA[2:Ny, 2:Nx]*dpsdt[2:Ny, 2:Nx]
#         + PHI[2:Ny, 2:Nx]*tk1[2:Ny, 2:Nx]*ALPHAF[2:Ny, 2:Nx] * dpfdt[2:Ny, 2:Nx]
#     )    
# # end # @timeit to "compute_adiabatic_heating!"
#     return nothing
# end # function compute_adiabatic_heating!


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
	- dtm: displacement time step
	- sp: static simulation parameters

# Returns

    - nothing
"""
function perform_thermal_iterations!(
    tk0, tk1, tk2, DT, DT0, RHOCP, KX, KY, HR, HA, HS, dtm, sp)
# @timeit to "assemble_hydromechanical_lse!" begin
    @unpack Nx,
        Ny,
        Nx1,
        Ny1,
        dx,
        dy,
        DTmax = sp
    # set up thermal iterations
    @. tk0 = tk1
    dtt = dtm
    dttsum = 0.0
    titer = 1
    RT = zeros(Ny1*Nx1)
    ST = zeros(Ny1*Nx1)
    # perform thermal iterations until reaching time limit
    while dttsum < dtm
        # reset LSE 
        LT = ExtendableSparseMatrix(Ny1*Nx1, Ny1*Nx1)
        @. RT = 0.0
        # compose global thermal matrix LT and coefficient vector RT
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
                # ρCₚ∂/∂t = k∇²T + Hᵣ + Hₛ + Hₐ
                #         = -∂qᵢ/∂xᵢ + Hᵣ + Hₛ + Hₐ (eq 10.9) 
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
                # fill system of equations: LHS (10.9)
                updateindex!(LT, +, -Kx₁/dx^2, gk, gk-Ny1) # T₁
                updateindex!(LT, +, -Ky₁/dy^2, gk, gk-1) # T₂
                updateindex!(LT, +, (
                    RHOCP[i, j]/dtt
                    + (Kx₁+Kx₂)/dx^2
                    + (Ky₁+Ky₂)/dy^2),
                    gk,
                    gk
                ) # T₃
                updateindex!(LT, +, -Ky₂/dy^2, gk, gk+1) # T₄
                updateindex!(LT, +, -Kx₂/dx^2, gk, gk+Ny1) # T₅
                # fill system of equations: RHS (10.9)
                RT[gk] = (
                    RHOCP[i, j]/dtt*tk1[i, j] + HR[i, j] + HA[i, j] + HS[i, j])

            end
        end
        # solve system of equations
        ST .= LT \ RT # implicit: flush!(LT)
        # reshape solution vector to 2D array
        tk2 .= reshape(ST, Ny1, Nx1)
        # compute ΔT
        @. DT = tk2 - tk1
        if titer == 1
            # during first thermal iteration pass:
            # apply thermal timestepping stability condition
            maxDTcurrent = maximum(abs, DT)
            if maxDTcurrent > DTmax
                dtt *= DTmax / maxDTcurrent
            else
                dttsum += dtt
            end
        else
            # second+ thermal iteration passes:
            # update dttsum and adjust timestep
            dttsum += dtt
            dtt = min(dtt, dtm-dttsum)
        end
        # increase thermal iteration counter
        titer += 1
    end
    # finalize overall temperature change and advance temperature field
    @. DT = tk2 - tk0
    @. DT0 = DT
# end # @timeit to "perform_thermal_iterations!"
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
    - dtm: displacement time step
	- marknum: total number of markers in use
	- sp: static simulation parameters

# Returns

    - nothing
"""
function apply_subgrid_temperature_diffusion!(
    xm, ym, tm, tkm, phim, tk1, DT, TKSUM, RHOCPSUM, dtm, marknum, sp)
# @timeit to "apply_subgrid_temperature_diffusion!" begin
    @unpack Nx1,
        Ny1,
        dx,
        dy,
        xp,
        yp,
        jmin_p,
        jmax_p,
        imin_p,
        imax_p,
        dsubgridt,
        rhocpsolidm,
        rhocpfluidm,
        ksolidm,
        kfluidm = sp
    # only perform subgrid temperature diffusion if enabled by dsubgridt > 0
    if dsubgridt == 0.0
        return nothing
    end
    # reset interpolation arrays
    TKSUM .= 0.0
    RHOCPSUM .= 0.0
    # iterate over markers
    @threads for m=1:1:marknum
        i, j, weights = fix_weights(
            xm[m], ym[m], xp, yp, dx, dy, jmin_p, jmax_p, imin_p, imax_p)
        # compute marker-node temperature difference
        δtkm = tkm[m] - dot(grid_vector(i, j, tk1), weights)
        # compute marker properties
        if tm[m] < 3
            # rocks
            rhocptotalm = total(rhocpsolidm[tm[m]], rhocpfluidm[tm[m]], phim[m])
            ktotalm = ktotal(ksolidm[tm[m]], kfluidm[tm[m]], phim[m])
        else
            # sticky air
            rhocptotalm = rhocpsolidm[tm[m]]
            ktotalm = ksolidm[tm[m]]
        end
        # time-relax δtkm difference
        δtkm *= (
            exp(-dsubgridt*ktotalm*dtm/rhocptotalm*(2.0/dx^2+2.0/dy^2)) - 1.0)
        # correct marker temperature
        tkm[m] += δtkm
        # update subgrid diffusion on P nodes
        interpolate_add_to_grid!(i, j, weights, δtkm*rhocptotalm, TKSUM)
        interpolate_add_to_grid!(i, j, weights, rhocptotalm, RHOCPSUM)
    end
    # reduce interpolation arrays
    reduce_add_3darray!(TKSUM)
    reduce_add_3darray!(RHOCPSUM)
    # compute DTsubgrid=TKSUM/RHOCPSUM and update temperature field at P nodes
    for j=1:1:Nx1, i=1:1:Ny1
        if RHOCPSUM[i, j, 1] > 0.0
            DT[i, j] -= TKSUM[i, j, 1] / RHOCPSUM[i, j, 1]
        end
    end
# end # @timeit to "apply_subgrid_temperature_diffusion!"
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
    - sp: static simulation parameters

# Returns

    - nothing
# """
function update_marker_temperature!(xm, ym, tkm, DT, tk2, timestep, marknum, sp)
# @timeit to "update_marker_temperature!" begin
    @unpack dx,
        dy,
        xp,
        yp,
        jmin_p,
        jmax_p,
        imin_p,
        imax_p = sp
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
# end # @timeit to "update_marker_temperature!"
    return nothing
end # function update_marker_temperature!


"""
Update marker porosity based on Dln[(1-ϕ)/ϕ]/Dt at P grid.

$(SIGNATURES)

# Details

    - xm: x-coordinates of markers
    - ym: y-coordinates of markers
    - tm: marker type
    - phim: marker porosity
    - APHI: Dln[(1-ϕ)/ϕ]/Dt at P nodes
    - dtm: displacement time step
    - marknum: total number of markers in use
    - sp: static simulation parameters

# Returns

    - nothing
# """
function update_marker_porosity!(xm, ym, tm, phim, APHI, dtm, marknum, sp)
# @timeit to "update_marker_porosity!" begin
    @unpack dx,
        dy,
        xp,
        yp,
        jmin_p,
        jmax_p,
        imin_p,
        imax_p,
        phimin,
        phimax = sp
    # interpolate and apply DT to markers for subsequent time steps
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
                aphim = dot(grid_vector(i, j, APHI), weights)
                # update marker porosity
                phim[m] = max(
                    phimin,
                    min(
                        phimax,
                        phim[m]/((1-phim[m])*exp(aphim*dtm)+phim[m])
                    )
                )
            end
        end
    end # @inbounds
# end # @timeit to "update_marker_porosity!"
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
function compute_velocities!(vx, vy, vxf, vyf, vxp, vyp, vxpf, vypf, sp)
@timeit to "compute_velocities!" begin
    @unpack Nx,
        Ny,
        Nx1,
        Ny1,
        dx,
        dy,
        bctop,
        bcbottom,
        bcleft,
        bcright,
        bcftop,
        bcfbottom,
        bcfleft,
        bcfright,
        vxleft,
        vxright,
        vytop,
        vybottom = sp
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
function compute_rotation_rate!(vx, vy, wyx, sp)
@timeit to "compute_rotation_rate!" begin
    @unpack Nx, Ny, dx, dy = sp
    # compute rotation rate ωyx=1/2[∂Vy/∂x-∂Vx/∂y] at basic nodes
    for j=1:1:Nx, i=1:1:Ny
        @inbounds wyx[i, j] = 0.5 * (
        (vy[i, j+1]-vy[i, j])/dx - (vx[i+1, j]-vx[i, j])/dy
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
    - tm: 
    - tkm: 
    - phim: 
    - sxym: 
    - sxxm: 
    - vx: 
    - vy: 
    - vxf: 
    - vyf: 
    - wyx: 
    - tk2: 
    - marknum: 
    - dtm: 
    - sp: static simulation parameters

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
	dtm,
	sp
)
# @timeit to "move_markers_rk4!" begin
    @unpack Nx,
        Ny,
        dx,
        dy,
        x,
        y,
        xvx,
        yvx,
        xvy,
        yvy,
        xp,
        yp,
        jmin_basic,
        jmax_basic,
        imin_basic,
        imax_basic,
        jmin_vx,
        jmax_vx,
        imin_vx,
        imax_vx,
        jmin_vy,
        jmax_vy,
        imin_vy,
        imax_vy,
        jmin_p,
        jmax_p,
        imin_p,
        imax_p,
        brk4,
        crk4,
        rhocpsolidm,
        rhocpfluidm = sp
    # @threads for m=1:1:marknum        
    for m=1:1:marknum        
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
        # interpolate local solid temperature from P grid
        tkms₀ = dot(grid_vector(i, j, tk2), weights)
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
        # interpolate local rotation rate from basic grid
        ωm = dot(grid_vector(i, j, wyx), weights)
        # incremental rotation angle
        θ = dtm * ωm
        # compute analytic stress rotation using σ′′xx = -σ′′yy
        @inbounds sxxm₁ = sxxm[m]*cos(θ)^2 - sxxm[m]*sin(θ)^2-sxym[m]*sin(2.0*θ)
        @inbounds sxym₁ = sxxm[m]*sin(2.0*θ) + sxym[m]*cos(2.0*θ)
        # update stresses
        @inbounds sxxm[m] = sxxm₁
        @inbounds sxym[m] = sxym₁
        # setup RK4 scheme 
        # RK4 coordinate positions A, B, C, D
        @inbounds xmrk4 = @MVector [xm[m], 0.0, 0.0, 0.0]
        @inbounds ymrk4 = @MVector [ym[m], 0.0, 0.0, 0.0]
        # RK4 velocities va, vb, vc, vd
        vxrk4 = @MVector zeros(4)
        vyrk4 = @MVector zeros(4)
        # advance marker using RK4 scheme on solid velocity
        for rk=1:1:4
            # interpolate vx
            @inbounds i, j, dxmj, dymi = fix_distances(
                xmrk4[rk],
                ymrk4[rk],
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
            @inbounds vxm₁₃ = vx[i, j]*(1.0-dxmj/dx) + vx[i, j+1]*dxmj/dx
            @inbounds vxm₂₄ = vx[i+1, j]*(1.0-dxmj/dx) + vx[i+1, j+1]*dxmj/dx
            # compute second order vx velocity corrections
            if dxmj/dx >= 0.5 
                # in right half of cell but not at right edge of grid
                if j < Nx-1
                    @inbounds vxm₁₃ += 0.5*((dxmj/dx-0.5)^2) * (
                        vx[i, j] - 2.0*vx[i, j+1] + vx[i, j+2])
                    @inbounds vxm₂₄ += 0.5*((dxmj/dx-0.5)^2) * (
                        vx[i+1, j] - 2.0*vx[i+1, j+1] + vx[i+1, j+2])
                end
            else
                # in left half of cell but not at left edge of grid
                if j > 1
                    @inbounds vxm₁₃ += 0.5*((dxmj/dx-0.5)^2) * (
                        vx[i, j-1] - 2.0*vx[i, j] + vx[i, j+1])
                    @inbounds vxm₂₄ += 0.5*((dxmj/dx-0.5)^2) * (
                        vx[i+1, j-1] - 2.0*vx[i+1, j] + vx[i+1, j+1])
                end
            end
            # compute current RK step vx
            @inbounds vxrk4[rk] = vxm₁₃*(1.0-dymi/dy) + vxm₂₄*dymi/dy
            # interpolate vy
            @inbounds i, j, dxmj, dymi = fix_distances(
                xmrk4[rk],
                ymrk4[rk],
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
            @inbounds vym₁₂ = vy[i, j]*(1.0-dymi/dy) + vy[i+1, j]*dymi/dy
            @inbounds vym₃₄ = vy[i, j+1]*(1.0-dymi/dy) + vy[i+1, j+1]*dymi/dy
            # compute second order vy velocity corrections
            if dymi/dy >= 0.5
                # in bottom half of cell but not at bottom edge of grid
                if i < Ny-1
                    @inbounds vym₁₂ += 0.5*((dymi/dy-0.5)^2) * (
                        vy[i, j] - 2.0*vy[i+1, j] + vy[i+2, j])
                    @inbounds vym₃₄ += 0.5*((dymi/dy-0.5)^2) * (
                        vy[i, j+1] - 2.0*vy[i+1, j+1] + vy[i+2, j+1])
                end      
            else
                # in top half of cell but not at top edge of grid
                if i > 1
                    @inbounds vym₁₂ += 0.5*((dymi/dy-0.5)^2) * (
                        vy[i-1, j] - 2.0*vy[i, j] + vy[i+1, j])
                    @inbounds vym₃₄ += 0.5*((dymi/dy-0.5)^2) * (
                        vy[i-1, j+1] - 2.0*vy[i, j+1] + vy[i+1, j+1])
                end
            end
            # compute current RK step vy
            @inbounds vyrk4[rk] = vym₁₂*(1.0-dxmj/dx) + vym₃₄*dxmj/dx
            # calculate next RK step x and y positions if not at final step
            if rk < 4
                @inbounds xmrk4[rk+1] = xmrk4[1] + dtm*crk4[rk]*vxrk4[rk]
                @inbounds ymrk4[rk+1] = ymrk4[1] + dtm*crk4[rk]*vyrk4[rk]
            end
        end # RK4 solid velocity loop
        # advance marker using RK4 solid velocity
        @inbounds xm[m] += dtm * dot(brk4, vxrk4)
        @inbounds ym[m] += dtm * dot(brk4, vyrk4)
        # reset RK4 scheme for fluid velocity backtracing
        @inbounds xmrk4[1] = xm[m]
        @inbounds ymrk4[1] = ym[m]
        # backtrack marker using RK4 scheme on fluid velocity
        for rk=1:1:4
            # interpolate vxf
            @inbounds i, j, dxmj, dymi = fix_distances(
                xmrk4[rk],
                ymrk4[rk],
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
            @inbounds vxfm₁₃ = vxf[i, j]*(1.0-dxmj/dx) + vxf[i, j+1]*dxmj/dx
            @inbounds vxfm₂₄ = vxf[i+1, j]*(1.0-dxmj/dx) + vxf[i+1, j+1]*dxmj/dx
            # compute second order vxf velocity corrections
            if dxmj/dx >= 0.5 
                # in right half of cell but not at right edge of grid
                if j < Nx-1
                    @inbounds vxfm₁₃ += 0.5*((dxmj/dx-0.5)^2) * (
                        vxf[i, j] - 2.0*vxf[i, j+1] + vxf[i, j+2])
                    @inbounds vxfm₂₄ += 0.5*((dxmj/dx-0.5)^2) * (
                        vxf[i+1, j] - 2.0*vxf[i+1, j+1] + vxf[i+1, j+2])
                end
            else
                # in left half of cell but not at left edge of grid
                if j > 1
                    @inbounds vxfm₁₃ += 0.5*((dxmj/dx-0.5)^2) * (
                        vxf[i, j-1] - 2.0*vx[i, j] + vxf[i, j+1])
                    @inbounds vxfm₂₄ += 0.5*((dxmj/dx-0.5)^2) * (
                        vxf[i+1, j-1] - 2.0*vxf[i+1, j] + vxf[i+1, j+1])
                end
            end
            # compute current RK step vxf
            @inbounds vxrk4[rk] = vxfm₁₃*(1.0-dymi/dy) + vxfm₂₄*dymi/dy
            # interpolate vyf
            @inbounds i, j, dxmj, dymi = fix_distances(
                xmrk4[rk],
                ymrk4[rk],
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
            @inbounds vyfm₁₂ = vyf[i, j]*(1.0-dymi/dy) + vyf[i+1, j]*dymi/dy
            @inbounds vyfm₃₄ = vyf[i, j+1]*(1.0-dymi/dy) + vyf[i+1, j+1]*dymi/dy
            # compute second order vyf velocity corrections
            if dymi/dy >= 0.5
                # in bottom half of cell but not at bottom edge of grid
                if i < Ny-1
                    @inbounds vyfm₁₂ += 0.5*((dymi/dy-0.5)^2) * (
                        vyf[i, j] - 2.0*vyf[i+1, j] + vyf[i+2, j])
                    @inbounds vyfm₃₄ += 0.5*((dymi/dy-0.5)^2) * (
                        vyf[i, j+1] - 2.0*vyf[i+1, j+1] + vyf[i+2, j+1])
                end
            else
                # in top half of cell but not at top edge of grid
                if i > 1
                    @inbounds vyfm₁₂ += 0.5*((dymi/dy-0.5)^2) * (
                        vyf[i-1, j] - 2.0*vyf[i, j] + vyf[i+1, j])
                    @inbounds vyfm₃₄ += 0.5*((dymi/dy-0.5)^2) * (
                        vyf[i-1, j+1] - 2.0*vyf[i, j+1] + vyf[i+1, j+1])
                end
            end
            # compute current RK step vyf
            @inbounds vyrk4[rk] = vyfm₁₂*(1.0-dxmj/dx) + vyfm₃₄*dxmj/dx
            # calculate next RK step x and y positions if not at final step
            if rk < 4
                @inbounds xmrk4[rk+1] = xmrk4[1] - dtm*crk4[rk]*vxrk4[rk]
                @inbounds ymrk4[rk+1] = ymrk4[1] - dtm*crk4[rk]*vyrk4[rk]
            end
        end # RK4 fluid velocity loop
        # backtrace marker using RK4 fluid velocity
        @inbounds xmrk4[1] -= dtm * dot(brk4, vxrk4)
        @inbounds ymrk4[1] -= dtm * dot(brk4, vyrk4)
        # interpolate fluid temperature at backtraced marker position
        i, j, weights = fix_weights(
            xmrk4[1],
            ymrk4[1],
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
        tkmf₀ = dot(grid_vector(i, j, tk2), weights)
        # compute marker fluid-solid temperature difference
        δtkmfs = tkmf₀ - tkms₀
        # correct marker temperature
        @inbounds tkm[m] = (
            (1.0-phim[m])*tkm[m]*rhocpsolidm[tm[m]]
                + phim[m]*(tkm[m]+δtkmfs)*rhocpfluidm[tm[m]]
        ) / ((1-phim[m])*rhocpsolidm[tm[m]] + phim[m]*rhocpfluidm[tm[m]])
    end # marker loop
# end # timeit to "move_markers_rk4!"
    return nothing
end # function move_markers_rk4!


"""
Backtrack pressure nodes using classic Runge-Kutta integration (RK4) 
to update total, solid, and fluid pressure under consideration of 
two-phase flow velocities.

$(SIGNATURES)

# Details

    - pr:
    - pr0:
    - ps:
    - ps0:
    - pf:
    - pf0:
    - vx: 
    - vy: 
    - vxf: 
    - vyf:
    - sp: static simulation parameters

# Returns

    - nothing
"""
function backtrace_pressures_rk4!(
   pr, pr0, ps, ps0, pf, pf0, vx, vy, vxf, vyf, dtm, sp)
# @timeit to "backtrace_pressures_rk4!" begin
    @unpack Nx,
        Ny,
        dx,
        dy,
        xsize,
        ysize,
        xvx,
        yvx,
        xvy,
        yvy,
        xp,
        yp,
        jmin_vx,
        jmax_vx,
        imin_vx,
        imax_vx,
        jmin_vy,
        jmax_vy,
        imin_vy,
        imax_vy,
        jmin_p,
        jmax_p,
        imin_p,
        imax_p,
        brk4,
        crk4 = sp
    # advance pressure generations
    pr0 .= pr
    ps0 .= ps
    pf0 .= pf
    # # setup RK4 scheme
    # xrk4 = @MVector zeros(4)
    # yrk4 = @MVector zeros(4)
    # # RK4 velocities va, vb, vc, vd
    # vxrk4 = @MVector zeros(4)
    # vyrk4 = @MVector zeros(4)
    # # backtrace P nodes: total and solid pressure
    # for jj=2:1:Nx, ii=2:1:Ny
    #     @inbounds xrk4[1] = xp[jj]
    #     @inbounds yrk4[1] = yp[ii]
        for jj=2:1:Nx
            for ii=2:1:Ny
        # setup RK4 scheme
        xrk4 = @MVector [xp[jj], 0.0, 0.0, 0.0]
        yrk4 = @MVector [yp[ii], 0.0, 0.0, 0.0]
        # RK4 velocities va, vb, vc, vd
        vxrk4 = @MVector zeros(4)
        vyrk4 = @MVector zeros(4)
        # backtrace P node using RK4 scheme on solid velocity
        for rk=1:1:4
            # interpolate vx
            i, j, dxmj, dymi = fix_distances(
                xrk4[rk],
                yrk4[rk],
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
            @inbounds vxm₁₃ = vx[i, j]*(1.0-dxmj/dx) + vx[i, j+1]*dxmj/dx
            @inbounds vxm₂₄ = vx[i+1, j]*(1.0-dxmj/dx) + vx[i+1, j+1]*dxmj/dx
            # compute second order vx velocity corrections
            if dxmj/dx >= 0.5 
                # in right half of cell but not at right edge of grid
                if j < Nx-1
                    @inbounds vxm₁₃ += 0.5*((dxmj/dx-0.5)^2) * (
                        vx[i, j] - 2.0*vx[i, j+1] + vx[i, j+2])
                    @inbounds vxm₂₄ += 0.5*((dxmj/dx-0.5)^2) * (
                        vx[i+1, j] - 2.0*vx[i+1, j+1] + vx[i+1, j+2])
                end
            else
                # in left half of cell but not at left edge of grid
                if j > 1
                    @inbounds vxm₁₃ += 0.5*((dxmj/dx-0.5)^2) * (
                        vx[i, j-1] - 2.0*vx[i, j] + vx[i, j+1])
                    @inbounds vxm₂₄ += 0.5*((dxmj/dx-0.5)^2) * (
                        vx[i+1, j-1] - 2.0*vx[i+1, j] + vx[i+1, j+1])
                end
            end
            # compute current RK step vx
            @inbounds vxrk4[rk] = vxm₁₃*(1.0-dymi/dy) + vxm₂₄*dymi/dy
            # interpolate vy
            @inbounds i, j, dxmj, dymi = fix_distances(
                xrk4[rk],
                yrk4[rk],
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
            @inbounds vym₁₂ = vy[i, j]*(1.0-dymi/dy) + vy[i+1, j]*dymi/dy
            @inbounds vym₃₄ = vy[i, j+1]*(1.0-dymi/dy) + vy[i+1, j+1]*dymi/dy
            # compute second order vy velocity corrections
            if dymi/dy >= 0.5
                # in bottom half of cell but not at bottom edge of grid
                if i < Ny-1
                    @inbounds vym₁₂ += 0.5*((dymi/dy-0.5)^2) * (
                        vy[i, j] - 2.0*vy[i+1, j] + vy[i+2, j])
                    @inbounds vym₃₄ += 0.5*((dymi/dy-0.5)^2) * (
                        vy[i, j+1] - 2.0*vy[i+1, j+1] + vy[i+2, j+1])
                end      
            else
                # in top half of cell but not at top edge of grid
                if i > 1
                    @inbounds vym₁₂ += 0.5*((dymi/dy-0.5)^2) * (
                        vy[i-1, j] - 2.0*vy[i, j] + vy[i+1, j])
                    @inbounds vym₃₄ += 0.5*((dymi/dy-0.5)^2) * (
                        vy[i-1, j+1] - 2.0*vy[i, j+1] + vy[i+1, j+1])
                end
            end
            # compute current RK step vy
            @inbounds vyrk4[rk] = vym₁₂*(1.0-dxmj/dx) + vym₃₄*dxmj/dx
            # calculate next RK step x and y positions if not at final step
            if rk < 4
                @inbounds xrk4[rk+1] = xrk4[1] - dtm*crk4[rk]*vxrk4[rk]
                @inbounds yrk4[rk+1] = yrk4[1] - dtm*crk4[rk]*vyrk4[rk]
            end
        end # RK4 solid velocity loop
        # backtrace P node using RK4 solid velocity
        @inbounds xrk4[1] -= dtm * dot(brk4, vxrk4)
        @inbounds yrk4[1] -= dtm * dot(brk4, vyrk4)
        # interpolate total and solid pressure at backtraced P nodes
        @inbounds i, j, weights = fix_weights(
            xrk4[1],
            yrk4[1],
            xp,
            yp,
            dx,
            dy,
            jmin_p,
            jmax_p,
            imin_p,
            imax_p
        )
        @inbounds pr0[ii, jj] = dot(grid_vector(i, j, pr), weights)
        @inbounds ps0[ii, jj] = dot(grid_vector(i, j, ps), weights)
    end # ii total and solid pressure loop
    end # jj total and solid pressure loop
    # backtrace P nodes: fluid pressure
    # for jj=2:1:Nx, ii=2:1:Ny
    #     xrk4[1] = xp[ii]
    #     yrk4[1] = yp[jj]
    for jj=2:1:Nx
        @threads for ii=2:1:Ny
        # setup RK4 scheme
        xrk4 = @MVector [xp[jj], 0.0, 0.0, 0.0]
        yrk4 = @MVector [yp[ii], 0.0, 0.0, 0.0]
        # RK4 velocities va, vb, vc, vd
        vxrk4 = @MVector zeros(4)
        vyrk4 = @MVector zeros(4)
        # backtrace P node using RK4 scheme on fluid velocity
        for rk=1:1:4
            # interpolate vxf
            @inbounds i, j, dxmj, dymi = fix_distances(
                xrk4[rk],
                yrk4[rk],
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
            @inbounds vxfm₁₃ = vxf[i, j]*(1.0-dxmj/dx) + vxf[i, j+1]*dxmj/dx
            @inbounds vxfm₂₄ = vxf[i+1, j]*(1.0-dxmj/dx) + vxf[i+1, j+1]*dxmj/dx
            # compute second order vxf velocity corrections
            if dxmj/dx >= 0.5 
                # in right half of cell but not at right edge of grid
                if j < Nx-1
                    @inbounds vxfm₁₃ += 0.5*((dxmj/dx-0.5)^2) * (
                        vxf[i, j] - 2.0*vxf[i, j+1] + vxf[i, j+2])
                    @inbounds vxfm₂₄ += 0.5*((dxmj/dx-0.5)^2) * (
                        vxf[i+1, j] - 2.0*vxf[i+1, j+1] + vxf[i+1, j+2])
                end
            else
                # in left half of cell but not at left edge of grid
                if j > 1
                    @inbounds vxfm₁₃ += 0.5*((dxmj/dx-0.5)^2) * (
                        vxf[i, j-1] - 2.0*vx[i, j] + vxf[i, j+1])
                    @inbounds vxfm₂₄ += 0.5*((dxmj/dx-0.5)^2) * (
                        vxf[i+1, j-1] - 2.0*vxf[i+1, j] + vxf[i+1, j+1])
                end
            end
            # compute current RK step vxf
            @inbounds vxrk4[rk] = vxfm₁₃*(1.0-dymi/dy) + vxfm₂₄*dymi/dy
            # interpolate vyf
            @inbounds i, j, dxmj, dymi = fix_distances(
                xrk4[rk],
                yrk4[rk],
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
            @inbounds vyfm₁₂ = vyf[i, j]*(1.0-dymi/dy) + vyf[i+1, j]*dymi/dy
            @inbounds vyfm₃₄ = vyf[i, j+1]*(1.0-dymi/dy) + vyf[i+1, j+1]*dymi/dy
            # compute second order vyf velocity corrections
            if dymi/dy >= 0.5
                # in bottom half of cell but not at bottom edge of grid
                if i < Ny-1
                    @inbounds vyfm₁₂ += 0.5*((dymi/dy-0.5)^2) * (
                        vyf[i, j] - 2.0*vyf[i+1, j] + vyf[i+2, j])
                    @inbounds vyfm₃₄ += 0.5*((dymi/dy-0.5)^2) * (
                        vyf[i, j+1] - 2.0*vyf[i+1, j+1] + vyf[i+2, j+1])
                end
            else
                # in top half of cell but not at top edge of grid
                if i > 1
                    @inbounds vyfm₁₂ += 0.5*((dymi/dy-0.5)^2) * (
                        vyf[i-1, j] - 2.0*vyf[i, j] + vyf[i+1, j])
                    @inbounds vyfm₃₄ += 0.5*((dymi/dy-0.5)^2) * (
                        vyf[i-1, j+1] - 2.0*vyf[i, j+1] + vyf[i+1, j+1])
                end
            end
            # compute current RK step vyf
            @inbounds vyrk4[rk] = vyfm₁₂*(1.0-dxmj/dx) + vyfm₃₄*dxmj/dx
            # calculate next RK step x and y positions if not at final step
            if rk < 4
                @inbounds xrk4[rk+1] = xrk4[1] - dtm*crk4[rk]*vxrk4[rk]
                @inbounds yrk4[rk+1] = yrk4[1] - dtm*crk4[rk]*vyrk4[rk]
            end
        end # RK4 fluid velocity loop
        # backtrace P node using RK4 fluid velocity
        @inbounds xrk4[1] -= dtm * dot(brk4, vxrk4)
        @inbounds yrk4[1] -= dtm * dot(brk4, vyrk4)
        # interpolate fluid pressure at backtraced P nodes
        @inbounds i, j, weights = fix_weights(
            xrk4[1],
            yrk4[1],
            xp,
            yp,
            dx,
            dy,
            jmin_p,
            jmax_p,
            imin_p,
            imax_p
        )
        @inbounds pf0[ii, jj] = dot(grid_vector(i, j, pf), weights)
    end # ii fluid pressure loop
    end # jj fluid pressure loop
# end # timeit to "backtrace_pressures_rk4!"
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
	- etavpm: viscoplastic viscosity  
	- mdis: minimum distance of marker launch anchor points to nearest marker 
	- mnum: number of marker nearest to marker launch anchor positions

# Returns

    - nothing
"""
function replenish_markers!(
    xm, ym, tm, tkm, phim, sxxm, sxym, etavpm, mdis, mnum; randomized=true)
# @timeit to "replenish_markers!" begin
    # reset marker population geometry tracker
    mdis .= mdis_init
    mnum .= 0
    # establish marker distribution
    for m=1:1:length(xm)
        i, j = fix(
            xm[m], ym[m], xxm, yym, dxm, dym, jmin_m, jmax_m, imin_m, imax_m)
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
                    push!(xm, xxm[j] + (rand()-0.5)*dxm)
                    push!(ym, yym[i] + (rand()-0.5)*dym)
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
            end
        end
    end    
    return nothing
    # return length(xm)
# end # timeit to "replenish_markers!"
end # function replenish_markers!


"""
Main simulation loop: run calculations with timestepping.

$(SIGNATURES)

# Details

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
        x, y,
        xvx, yvx,
        xvy, yvy,
        xp, yp,
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
        dsubgrids,
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
    # x, y, xvx, yvx, xvy, yvy, xp, yp = setup_staggered_grid_geometry(sp)
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
    (
        ETA5,
        ETA00,
        YNY5,
        YNY00,
        YNY_inv_ETA,
        DSXY,
        ETAcomp,
        SXYcomp,
        dRHOXdx,
        dRHOXdy,
        dRHOYdx,
        dRHOYdy,
        ETAPcomp,
        SXXcomp,
        SYYcomp,
        EII,
        SII,
        DSXX
    ) = setup_staggered_grid_properties_helpers(sp)


    # -------------------------------------------------------------------------
    # set up markers
    # -------------------------------------------------------------------------
    mdis, mnum, mtyp, mpor = setup_marker_geometry_helpers(sp)
    xm, ym, tm, tkm, sxxm, sxym, etavpm, phim = setup_marker_properties(sp)
    (
        rhototalm,
        rhocptotalm,
        etasolidcur,
        etafluidcur,
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
    ) = setup_marker_properties_helpers(sp)
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
    # set up of matrices for global gravity/thermal/hydromechanical solvers
    # -------------------------------------------------------------------------
    # hydromechanical solver
    L, R, S = setup_hydromechanical_lse(sp)
    # thermal solver
    LT, RT, ST = setup_thermal_lse(sp)
    # gravitational solver
    LP, RP, SP= setup_gravitational_lse(sp)
# end # @timeit to "simulation_loop setup"

    # -------------------------------------------------------------------------
    # iterate timesteps   
    # -------------------------------------------------------------------------
    nsteps = 10 # <======= RMK: remove for production
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
            WTPSUM
        ) = setup_interpolated_properties(sp)
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
            marker_to_basic_nodes!(
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
                WTSUM,
                sp
            )
            # interpolate marker properties to Vx nodes
            marker_to_vx_nodes!(
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
                WTXSUM,
                sp
            )
            # interpolate marker properties to Vy nodes
            marker_to_vy_nodes!(
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
                WTYSUM,
                sp
            )     
            # interpolate marker properties to P nodes
            marker_to_p_nodes!(
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
                WTPSUM,
                sp
            )
        end # @threads for m=1:1:marknum

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
                dt,
                sp
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
                dt,
                sp
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
                    YNY_inv_ETA,
                    dt,
                    dtkoef,
                    dtstep,
                    iplast
                )
            end
        end # for iplast=1:1:nplast

        # ---------------------------------------------------------------------
        # # interpolate updated viscoplastic viscosity to markers
        # ---------------------------------------------------------------------
        @threads for m = 1:1:marknum
            update_marker_viscosity!(
                m, xm, ym, tm, tkm, etatotalm, etavpm, YNY, YNY_inv_ETA, sp)
        end

        # ---------------------------------------------------------------------
        # # apply subgrid stress diffusion to markers
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
            dtm,
            marknum,
            sp
        )

        # ---------------------------------------------------------------------
        # interpolate DSXX, DSXY to markers
        # ---------------------------------------------------------------------
        update_marker_stress!(xm, ym, sxxm, sxym, DSXX, DSXY, marknum, sp)

        # ---------------------------------------------------------------------
        # compute shear heating HS in P nodes
        # ---------------------------------------------------------------------
        compute_shear_heating!(
            HS,
            SXYEXY,
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
            pf,
            sp
        )

        # ---------------------------------------------------------------------
        # compute adiabatic heating HA in P nodes
        # perform thermal iterations
        # ---------------------------------------------------------------------
        if timestep==1
            # no pressure changes for the first timestep
            pr0 .= pr
            pf0 .= pf
            ps0 .= pf
        end
        perform_thermal_iterations!(
            tk0, tk1, tk2, DT, DT0, RHOCP, KX, KY, HR, HA, HS, dtm, sp)

        # ---------------------------------------------------------------------
        # apply subgrid temperature diffusion on markers
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
            dtm,
            marknum,
            sp
        )

        # ---------------------------------------------------------------------
        # interpolate DT to markers
        # ---------------------------------------------------------------------
        update_marker_temperature!(xm, ym, tkm, DT, tk2, timestep, marknum, sp)

        # ---------------------------------------------------------------------
        # update porosity on markers
        # ---------------------------------------------------------------------
        update_marker_porosity!(xm, ym, tm, phim, APHI, dtm, marknum, sp)

        # ---------------------------------------------------------------------
        # compute velocity in P nodes
        # compute fluid velocity in P nodes including boundary conditions
        # ---------------------------------------------------------------------
        compute_velocities!(vx, vy, vxf, vyf, vxp, vyp, vxpf, vypf, sp)

        # ---------------------------------------------------------------------
        # compute rotation rate in basic nodes
        # ---------------------------------------------------------------------
        compute_rotation_rate!(vx, vy, wyx, sp)

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
            dtm,
            sp
        )

        # ---------------------------------------------------------------------
        # backtrack P nodes: Ptotal with RK4
        # backtrack P nodes: Pfluid with RK4
        # ---------------------------------------------------------------------
        backtrace_pressures_rk4!(
            pr, pr0, ps, ps0, pf, pf0, vx, vy, vxf, vyf, dtm, sp)

        # ---------------------------------------------------------------------
        # # replenish sparse areas with additional markers
        # ---------------------------------------------------------------------
        # 2489-2604
        # for m = 1:1:marknum
        #     # ~100 lines MATLAB
        # end

        # ---------------------------------------------------------------------
        # # update timesum
        # ---------------------------------------------------------------------

        # ---------------------------------------------------------------------
        # # save data for analysis and visualization
        # ---------------------------------------------------------------------

        # ---------------------------------------------------------------------
        # # save old stresses - RMK: not used anywhere in code ?
        # ---------------------------------------------------------------------
        # # sxxm00 = sxxm 
        # # sxym00 = sxym    

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

$(SIGNATURES)

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
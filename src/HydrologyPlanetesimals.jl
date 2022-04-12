module HydrologyPlanetesimals

using Base.Threads
using SparseArrays
using MAT
using DocStringExtensions
using Parameters
using StaticArrays
using ProgressMeter
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
    xsize::Float64
    "vertical model size [m]"
    ysize::Float64
    "horizontal center of model"
    xcenter::Float64 = xsize / 2
    "vertical center of model"
    ycenter::Float64 = ysize / 2  
    "basic grid resolution in x direction (horizontal)"
    Nx::Int 
    "basic grid resolution in y direction (vertical)"	
    Ny::Int
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
    rplanet::Int64
    "crust radius [m]"
    rcrust::Int64
    "surface pressure [Pa]"
    psurface::Float64 = 1e+3
    # marker count and initial spacing
    "number of markers per cell in horizontal direction"
    Nxmc::Int
    "number of markers per cell in vertical direction"
    Nymc::Int
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
    sp::StaticParameters
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
Computers properties of given marker and saves them to corresponding arrays.

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
@timeit to "compute_marker_properties!" begin
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
end
    return nothing
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
    return (ksolid * kfluid/2 + ((ksolid * (3*phi-2)
                                 + kfluid * (1.0-3.0*phi))^2)/16)^0.5
            - (ksolid*(3.0*phi-2.0) + kfluid*(1.0-3.0*phi))/4
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
    - sp: static simulation parameters
    - dp: dynamic simulation parameters

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
        hrsolidm = @SVector [0.0, 0.0, 0.0]
    end    
    #60Fe: planet ✓, crust ×, space ×
    if hr_fe
        # 60Fe radiogenic heat production [W/kg]
        Q_fe = Q_radiogenic(f_fe, ratio_fe, E_fe, tau_fe, timesum)
        # Fluid phase 60Fe radiogenic heat production [W/m^3]
        hrfluidm = @SVector [Q_fe*rhofluidm[1], 0.0, 0.0]
    else
        hrfluidm = @SVector [0.0, 0.0, 0.0]
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
@timeit to "fix_weights" begin
    @inbounds j = trunc(Int, (x - x_axis[1]) / dx) + 1
    @inbounds i = trunc(Int, (y - y_axis[1]) / dy) + 1
    @inbounds dxmj = x - x_axis[min(max(j, jmin), jmax)]
    @inbounds dymi = y - y_axis[min(max(i, imin), imax)]
    weights = SVector(
        (1.0-dymi/dy) * (1.0-dxmj/dx),
        (dymi/dy) * (1.0-dxmj/dx),
        (1.0-dymi/dy) * (dxmj/dx),
        (dymi/dy) * (dxmj/dx)
        )
end
    return i, j, weights
end


"""
Interpolate a property to neareast four nodes on a given grid location
using given bilinear interpolation weights.

# Details

    - i: top (with reference to y) node index on y-grid axis
    - j: left (with reference to x) node index on x-grid axis
    - weights: vector of 4 bilinear interpolation weights to
      nearest four grid nodes:
        [wtmij  : i  , j   node,
        wtmi1j : i+1, j   node,
        wtmij1 : i  , j+1 node,
        wtmi1j1: i+1, j+1 node]
    - property: property to be interpolated to grid using weights
    - grid: threaded grid array on which to interpolate property
"""
function interpolate!(i, j, weights, property, grid)
@timeit to "interpolate!" begin
    grid[i, j, threadid()] += property * weights[1]
    grid[i+1, j, threadid()] += property * weights[2]
    grid[i, j+1, threadid()] += property * weights[3]
    grid[i+1, j+1, threadid()] += property * weights[4]
end
end


"""
Compute properties of basic nodes based on interpolation arrays.

$(SIGNATURES)

# Detail

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
    - YNY: YNY basic node array
    - GGG: GGG basic node array
    - SXY0: SXY basic node array
    - COH: COH basic node array
    - TEN: TEN basic node array
    - FRI: FRI basic node array

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
    YNY,
    GGG,
    SXY0,
    COH,
    TEN,
    FRI
)
@timeit to "compute_basic_node_properties!" begin
@timeit to "reduce" begin
    ETA0SUM = reduce(+, ETA0SUM, dims=3)[:, :, 1]
    ETASUM = reduce(+, ETASUM, dims=3)[:, :, 1]
    GGGSUM = reduce(+, GGGSUM, dims=3)[:, :, 1]
    SXYSUM = reduce(+, SXYSUM, dims=3)[:, :, 1]
    COHSUM = reduce(+, COHSUM, dims=3)[:, :, 1]
    TENSUM = reduce(+, TENSUM, dims=3)[:, :, 1]
    FRISUM = reduce(+, FRISUM, dims=3)[:, :, 1]
    WTSUM = reduce(+, WTSUM, dims=3)[:, :, 1]
end # @timeit to "reduce"
@timeit to "compute" begin
    WTSUM[WTSUM .== 0.0] .= 1.0
    ETA0 .= ETA0SUM ./ WTSUM
    ETA .= ETASUM ./ WTSUM
    YNY[ETA .< ETA0] .= true
    GGG .= GGGSUM .\ WTSUM
    SXY0 .= SXYSUM ./ WTSUM
    COH .= COHSUM ./ WTSUM
    TEN .= TENSUM ./ WTSUM
    FRI .= FRISUM ./ WTSUM
end # @timeit to "compute"
end # @timeit to "compute_basic_node_properties!"
end # function compute_basic_node_properties!


"""
Compute properties of Vx nodes based on interpolation arrays.

$(SIGNATURES)

# Detail

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
@timeit to "reduce" begin
    RHOXSUM = reduce(+, RHOXSUM, dims=3)[:, :, 1]
    RHOFXSUM = reduce(+, RHOFXSUM, dims=3)[:, :, 1]
    KXSUM = reduce(+, KXSUM, dims=3)[:, :, 1]
    PHIXSUM = reduce(+, PHIXSUM, dims=3)[:, :, 1]
    RXSUM = reduce(+, RXSUM, dims=3)[:, :, 1]
    WTXSUM = reduce(+, WTXSUM, dims=3)[:, :, 1]    
end # @timeit to "reduce"
@timeit to "compute" begin
    WTXSUM[WTXSUM .== 0.0] .= 1.0
    RHOX .= RHOXSUM ./ WTXSUM
    RHOFX .= RHOFXSUM ./ WTXSUM
    KX .= KXSUM ./ WTXSUM
    PHIX .= PHIXSUM ./ WTXSUM
    RX .= RXSUM ./ WTXSUM
end # @timeit to "compute"
end # @timeit to "compute_vx_node_properties!"
end # function compute_vx_node_properties!


"""
Compute properties of Vy nodes based on interpolation arrays.

$(SIGNATURES)

# Detail

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
@timeit to "reduce" begin
    RHOYSUM = reduce(+, RHOYSUM, dims=3)[:, :, 1]
    RHOFYSUM = reduce(+, RHOFYSUM, dims=3)[:, :, 1]
    KYSUM = reduce(+, KYSUM, dims=3)[:, :, 1]
    PHIYSUM = reduce(+, PHIYSUM, dims=3)[:, :, 1]
    RYSUM = reduce(+, RYSUM, dims=3)[:, :, 1]
    WTYSUM = reduce(+, WTYSUM, dims=3)[:, :, 1]    
end # @timeit to "reduce"
@timeit to "compute" begin
    WTYSUM[WTYSUM .== 0.0] .= 1.0
    RHOY .= RHOYSUM ./ WTYSUM
    RHOFY .= RHOFYSUM ./ WTYSUM
    KY .= KYSUM ./ WTYSUM
    PHIY .= PHIYSUM ./ WTYSUM
    RY .= RYSUM ./ WTYSUM
end # @timeit to "compute"
end # @timeit to "compute_vy_node_properties!"
end # function compute_vy_node_properties!


"""
Compute properties of P nodes based on interpolation arrays.

$(SIGNATURES)

# Detail

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
    GGGP,
    SXX0,
    RHO,
    RHOCP,
    ALPHA,
    ALPHAF,
    HR,
    PHI,
    BETTAPHI,
    tk1  
)
@timeit to "compute_p_node_properties!" begin
@timeit to "reduce" begin
    GGGPSUM = reduce(+, GGGPSUM, dims=3)[:, :, 1]
    SXXSUM = reduce(+, SXXSUM, dims=3)[:, :, 1]
    RHOSUM = reduce(+, RHOSUM, dims=3)[:, :, 1]
    RHOCPSUM = reduce(+, RHOCPSUM, dims=3)[:, :, 1]
    ALPHASUM = reduce(+, ALPHASUM, dims=3)[:, :, 1]
    ALPHAFSUM = reduce(+, ALPHAFSUM, dims=3)[:, :, 1]
    HRSUM = reduce(+, HRSUM, dims=3)[:, :, 1]
    PHISUM = reduce(+, PHISUM, dims=3)[:, :, 1]
    TKSUM = reduce(+, TKSUM, dims=3)[:, :, 1]
    WTPSUM = reduce(+, WTPSUM, dims=3)[:, :, 1]
end # @timeit to "reduce"
@timeit to "compute" begin
    GGGP .= GGGPSUM .\ WTPSUM
    SXX0 .= SXXSUM ./ WTPSUM
    RHO .= RHOSUM ./ WTPSUM
    RHOCP .= RHOCPSUM ./ WTPSUM
    ALPHA .= ALPHASUM ./ WTPSUM
    ALPHAF .= ALPHAFSUM ./ WTPSUM
    HR .= HRSUM ./ WTPSUM
    PHI .= PHISUM ./ WTPSUM
    BETTAPHI .= GGGP .\ PHI
    tk1 .= TKSUM ./ RHOCPSUM
end # @timeit to "compute"
end # @timeit to "compute_p_node_properties!"
end # function compute_p_node_properties!


"""
Main simulation loop: run calculations with timestepping.

$(SIGNATURES)

# Detail

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
    start_step,
    nsteps,
    start_time, 
    endtime,
    start_marknum = sp

@timeit to "simulation_loop setup" begin
    # -------------------------------------------------------------------------
    # set up dynamic simulation parameters from given static parameters
    # -------------------------------------------------------------------------
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
    # Yielding error of nodes
    YERRNOD = zeros(Float64, nplast) 
   

    # -------------------------------------------------------------------------
    # set up staggered grid
    # -------------------------------------------------------------------------
    # basic nodes
    # grid geometry
    # x: horizontal coordinates of basic grid points [m]
    # x = @SVector [j for j = 0:dx:xsize] # should work but doesn't
    x = SVector{Nx, Float64}([j for j = 0:dx:xsize])
    # y: vertical coordinates of basic grid points [m]
    # y = @SVector [i for i = 0:dy:ysize]
    y = SVector{Ny, Float64}([j for j = 0:dy:ysize])
    # physical node properties
    # viscoplastic viscosity, Pa*s
    ETA = zeros(Float64, Ny, Nx)
    # viscous viscosity, Pa*s
    ETA0 = zeros(Float64, Ny, Nx)
    # shear modulus, Pa
    GGG = zeros(Float64, Ny, Nx)
    # epsilonxy, 1/s
    EXY = zeros(Float64, Ny, Nx)
    # sigma0xy, 1/s
    SXY0 = zeros(Float64, Ny, Nx)
    # rotation rate, 1/s
    WYX = zeros(Float64, Ny, Nx)
    # compressive strength, Pa
    COH = zeros(Float64, Ny, Nx)
    # tensile strength, Pa
    TEN = zeros(Float64, Ny, Nx)
    # friction
    FRI = zeros(Float64, Ny, Nx)
    # plastic yielding mark
    YNY = zeros(Bool, Ny, Nx)

    # Vx nodes
    # grid geometry
    # xvx: horizontal coordinates of vx grid points [m]
    xvx = SVector{Ny1, Float64}([j for j = 0:dx:xsize+dy])
    # yvx: vertical coordinates of vx grid points [m]
    yvx = SVector{Nx1, Float64}([i for i = -dy/2:dy:ysize+dy/2])
    # physical node properties
    # density [kg/m^3]
    RHOX = zeros(Float64, Ny1, Nx1)
    # fluid density [kg/m^3]
    RHOFX = zeros(Float64, Ny1, Nx1)
    # thermal conductivity [W/m/K]
    KX = zeros(Float64, Ny1, Nx1)
    # porosity
    PHIX = zeros(Float64, Ny1, Nx1)
    # solid vx-velocity [m/s]
    vx = zeros(Float64, Ny1, Nx1)
    # fluid vx-velocity [m/s]
    vxf = zeros(Float64, Ny1, Nx1)
    # etafluid/kphi ratio [m^2]
    RX = zeros(Float64, Ny1, Nx1)
    # qx-darcy flux [m/s]
    qxD = zeros(Float64, Ny1, Nx1)
    # gx-gravity [m/s^2]
    gx = zeros(Float64, Ny1, Nx1)

    # Vy nodes
    # grid geometry
    # xvy: horizontal coordinates of vy grid points [m]
    xvy = SVector{Nx1, Float64}([j for j = -dx/2:dx:xsize+dx/2])
    # yvy: vertical coordinates of vy grid points [m]
    yvy = SVector{Ny1, Float64}([i for i = 0:dy:ysize+dy])
    # physical node properties
    # "density [kg/m^3]"
    RHOY = zeros(Float64, Ny1, Nx1)
    # "fluid density [kg/m^3]"
    RHOFY = zeros(Float64, Ny1, Nx1)
    # "thermal conductivity [W/m/K]"
    KY = zeros(Float64, Ny1, Nx1)
    # "porosity"
    PHIY = zeros(Float64, Ny1, Nx1)
    # "solid vy-velocity [m/s]"
    vy = zeros(Float64, Ny1, Nx1)
    # "fluid vy-velocity [m/s]"
    vyf = zeros(Float64, Ny1, Nx1)
    # "etafluid/kphi ratio [m^2]"
    RY = zeros(Float64, Ny1, Nx1)
    # "qy-darcy flux [m/s]"
    qyD = zeros(Float64, Ny1, Nx1)
    # "gy-gravity [m/s^2]"
    gy = zeros(Float64, Ny1, Nx1)

    # P nodes
    # grid geometry
    # xp: horizontal coordinates of p grid points [m]
    xp = SVector{Nx1, Float64}([j for j = -dx/2:dx:xsize+dx/2])
    # yp: vertical coordinates of p grid points [m]
    yp = SVector{Ny1, Float64}([i for i = -dy/2:dy:ysize+dy/2])
    # physical node properties
    # density [kg/m^3]
    RHO = zeros(Float64, Ny1, Nx1)
    # volumetric heat capacity [J/m^3/K]
    RHOCP = zeros(Float64, Ny1, Nx1)
    # thermal expansion [J/m^3/K]
    ALPHA = zeros(Float64, Ny1, Nx1)
    # fluid thermal expansion [J/m^3/K]
    ALPHAF = zeros(Float64, Ny1, Nx1)
    # radioactive heating [W/m^3]
    HR = zeros(Float64, Ny1, Nx1)
    # adiabatic heating [W/m^3]
    HA = zeros(Float64, Ny1, Nx1)
    # shear heating [W/m^3]
    HS = zeros(Float64, Ny1, Nx1)
    # viscosity [Pa*s]
    ETAP = zeros(Float64, Ny1, Nx1)
    # shear modulus [Pa]
    GGGP = zeros(Float64, Ny1, Nx1)
    # EPSILONxx [1/s]
    EXX = zeros(Float64, Ny1, Nx1)
    # SIGMA'xx [1/s]
    SXX = zeros(Float64, Ny1, Nx1)
    # SIGMA0'xx [1/s]
    SXX0 = zeros(Float64, Ny1, Nx1)
    # old temperature [K]
    tk1 = zeros(Float64, Ny1, Nx1)
    # new temperature [K]
    tk2 = zeros(Float64, Ny1, Nx1)
    # solid Vx in pressure nodes [m/s]
    vxp = zeros(Float64, Ny1, Nx1)
    # solid Vy in pressure nodes [m/s]
    vyp = zeros(Float64, Ny1, Nx1)
    # fluid Vx in pressure nodes [m/s]
    vxpf = zeros(Float64, Ny1, Nx1)
    # fluid Vy in pressure nodes [m/s]
    vypf = zeros(Float64, Ny1, Nx1)
    # total pressure [Pa]
    pr = zeros(Float64, Ny1, Nx1)
    # fluid pressure [Pa]
    pf = zeros(Float64, Ny1, Nx1)
    # solid pressure [Pa]
    ps = zeros(Float64, Ny1, Nx1)
    # old total pressure [Pa]
    pr0 = zeros(Float64, Ny1, Nx1)
    # old fluid pressure [Pa]
    pf0 = zeros(Float64, Ny1, Nx1)
    # old solid pressure [Pa]
    ps0 = zeros(Float64, Ny1, Nx1)
    # bulk viscosity [Pa*s]
    ETAPHI = zeros(Float64, Ny1, Nx1)
    # bulk compresibility [Pa*s]
    BETTAPHI = zeros(Float64, Ny1, Nx1)
    # porosity
    PHI = zeros(Float64, Ny1, Nx1)
    # Dln[(1-PHI)/PHI]/Dt
    APHI = zeros(Float64, Ny1, Nx1)
    # gravity potential [J/kg]
    FI = zeros(Float64, Ny1, Nx1)


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
    L = spzeros(Nx1*Ny1*6, Nx1*Ny1*6)
    # hydromechanical solution: RHS Vector
    R = zeros(Float64, Nx1*Ny1*6)
    # thermal solution: LHS coefficient matrix
    LT = spzeros(Nx1*Ny1, Nx1*Ny1)
    # thermal solution: RHS Vector
    RT = zeros(Float64, Nx1*Ny1)
    # gravity solution: LHS coefficient matrix
    LP = spzeros(Nx1*Ny1, Nx1*Ny1)
    # gravity solution: RHS Vector
    RP = zeros(Float64, Nx1*Ny1)
end

    # -------------------------------------------------------------------------
    # iterate timesteps   
    # -------------------------------------------------------------------------
    nsteps = 100 # <======= remove for production
    p = Progress(
        nsteps,
        dt=0.5,
        barglyphs=BarGlyphs(
            '|','█', ['▁' ,'▂' ,'▃' ,'▄' ,'▅' ,'▆', '▇'],' ','|',), barlen=10)
    for timestep = start_step:1:nsteps
@timeit to "set up interpolation arrays" begin
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
        GGGPSUM = zeros(Ny1, Nx1, nthreads())
        SXXSUM = zeros(Ny1, Nx1, nthreads())
        RHOSUM = zeros(Ny1, Nx1, nthreads())
        RHOCPSUM = zeros(Ny1, Nx1, nthreads())
        ALPHASUM = zeros(Ny1, Nx1, nthreads())
        ALPHAFSUM = zeros(Ny1, Nx1, nthreads())
        HRSUM = zeros(Ny1, Nx1, nthreads())
        PHISUM = zeros(Ny1, Nx1, nthreads())
        TKSUM = zeros(Ny1, Nx1, nthreads())
        WTPSUM = zeros(Ny1, Nx1, nthreads())
end

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
            interpolate!(i, j, weights, etatotalm[m], ETA0SUM)
            # ETASUM: viscoplastic viscosity interpolated to basic nodes
            interpolate!(i, j, weights, etavpm[m], ETASUM)
            # GGGSUM: shear modulus interpolated to basic nodes
            interpolate!(i, j, weights, inv_gggtotalm[m], GGGSUM)
            # SXYSUM: σxy shear stress interpolated to basic nodes
            interpolate!(i, j, weights, sxym[m], SXYSUM)
            # COHSUM: compressive strength interpolated to basic nodes
            interpolate!(i, j, weights, cohestotalm[m], COHSUM)
            # TENSUM: tensile strength interpolated to basic nodes
            interpolate!(i, j, weights, tenstotalm[m], TENSUM)
            # FRISUM: friction  interpolated to basic nodes
            interpolate!(i, j, weights, fricttotalm[m], FRISUM)
            # WTSUM: weight array for bilinear interpolation to basic nodes
            interpolate!(i, j, weights, 1.0, WTSUM)

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
            interpolate!(i, j, weights, rhototalm[m], RHOXSUM)
            # RHOFXSUM: fluid density interpolated to Vx nodes
            interpolate!(i, j, weights, rhofluidcur[m], RHOFXSUM)
            # KXSUM: thermal conductivity interpolated to Vx nodes
            interpolate!(i, j, weights, ktotalm[m], KXSUM)
            # PHIXSUM: porosity interpolated to Vx nodes
            interpolate!(i, j, weights, phim[m], PHIXSUM)
            # RXSUM: ηfluid/kϕ interpolated to Vx nodes
            interpolate!(i, j, weights, etafluidcur_inv_kphim[m], RXSUM)
            # WTXSUM: weight for bilinear interpolation to Vx nodes
            interpolate!(i, j, weights, 1.0, WTXSUM)

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
            interpolate!(i, j, weights, rhototalm[m], RHOYSUM)
            # RHOFYSUM: fluid density interpolated to Vy nodes
            interpolate!(i, j, weights, rhofluidcur[m], RHOFYSUM)
            # KYSUM: thermal conductivity interpolated to Vy nodes
            interpolate!(i, j, weights, ktotalm[m], KYSUM)
            # PHIYSUM: porosity interpolated to Vy nodes
            interpolate!(i, j, weights, phim[m], PHIYSUM)
            # RYSUM: ηfluid/kϕ interpolated to Vy nodes
            interpolate!(i, j, weights, etafluidcur_inv_kphim[m], RYSUM)
            # WTYSUM: weight for bilinear interpolation to Vy nodes
            interpolate!(i, j, weights, 1.0, WTYSUM)
            
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
            interpolate!(i, j, weights, inv_gggtotalm[m], GGGPSUM)
            # SXXSUM: σ'xx interpolated to P nodes
            interpolate!(i, j, weights, sxxm[m], SXXSUM)
            # RHOSUM: density interpolated to P nodes
            interpolate!(i, j, weights, rhototalm[m], RHOSUM)
            # RHOCPSUM: volumetric heat capacity interpolated to P nodes
            interpolate!(i, j, weights, rhocptotalm[m], RHOCPSUM)
            # ALPHASUM: thermal expansion interpolated to P nodes
            interpolate!(i, j, weights, alphasolidcur[m], ALPHASUM)
            # ALPHAFSUM: fluid thermal expansion interpolated to P nodes
            interpolate!(i, j, weights, alphafluidcur[m], ALPHAFSUM)
            # HRSUM: radioactive heating interpolated to P nodes
            interpolate!(i, j, weights, hrtotalm[m], HRSUM)
            # PHISUM: porosity interpolated to P nodes
            interpolate!(i, j, weights, phim[m], PHISUM)
            # TKSUM: heat capacity interpolated to P nodes
            interpolate!(i, j, weights, tkm_rhocptotalm[m], TKSUM)
            # WTPSUM: weight for bilinear interpolation to P nodes
            interpolate!(i, j, weights, 1.0, WTPSUM)
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
            YNY,
            GGG,
            SXY0,
            COH,
            TEN,
            FRI 
        )


        # ---------------------------------------------------------------------
        # # compute physical properties of Vx nodes
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
        # # compute physical properties of Vy nodes
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
        # # compute physical properties of P nodes
        # ---------------------------------------------------------------------
        compute_p_node_properties!(
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
            GGGP,
            SXX0,
            RHO,
            RHOCP,
            ALPHA,
            ALPHAF,
            HR,
            PHI,
            BETTAPHI,
            tk1  
        )


# Tue 12        
        # ---------------------------------------------------------------------
        # # applying thermal boundary conditions for interpolated temperature
        # ---------------------------------------------------------------------
        # apply_thermal_bc(sp, dp, tk1)


# Wed 13        
        # ---------------------------------------------------------------------
        # # compute gravity solution
        # ---------------------------------------------------------------------
        # compute_gravity_solution!(sp, dp, interp_arrays)


# Wed 13        
        # ---------------------------------------------------------------------
        # # compute gravitational acceleration
        # ---------------------------------------------------------------------
        # compute_grav_accel!(sp, dp, interp_arrays)


# Wed 13        
        # ---------------------------------------------------------------------
        # # probe increasing computational timestep
        # ---------------------------------------------------------------------
        # dt = min(dt*dtkoefup, dtelastic)


# Wed 13        
        # ---------------------------------------------------------------------
        # # save initial viscosity, yielding nodes
        # ---------------------------------------------------------------------


# Thu 14        
        # ---------------------------------------------------------------------
        # # perform plastic iterations
        # ---------------------------------------------------------------------
        # for iplast = 1:1:nplast
        #     # ~600 lines MATLAB
        # end


# Thu 14
        # ---------------------------------------------------------------------
        # # interpolate updated viscoplastic viscosity to markers
        # ---------------------------------------------------------------------
        # for m = 1:1:marknum
        #     # ~50 lines MATLAB 
        # end


# Thu 14
        # ---------------------------------------------------------------------
        # # apply subgrid stress diffusion to markers
        # ---------------------------------------------------------------------
        # for m = 1:1:marknum
        #     # ~100 lines MATLAB 
        # end


# Thu 14
        # ---------------------------------------------------------------------
        # # compute DSXXsubgrid, DSXYsubgrid
        # ---------------------------------------------------------------------
        # compute_dsxx_dsxy_subgrids!(sp, dp, interp_arrays)


# Thu 14
        # ---------------------------------------------------------------------
        # # interpolate DSXX, DSXY to markers
        # ---------------------------------------------------------------------
        # for m = 1:1:marknum
        #     # ~50 lines MATLAB 
        # end

# Fri 15
        # ---------------------------------------------------------------------
        # # compute shear heating HS in P nodes
        # ---------------------------------------------------------------------
        # compute_HS_p_nodes!(sp, dp, interp_arrays)


# Fri 15
        # ---------------------------------------------------------------------
        # # compute adiabatic heating HA in P nodes
        # ---------------------------------------------------------------------
        # compute_HA_p_nodes!(sp, dp, interp_arrays)


# Fri 15
        # ---------------------------------------------------------------------
        # # perform thermal iterations
        # ---------------------------------------------------------------------
        # # ~100 lines MATLAB

# Sat 16
        # ---------------------------------------------------------------------
        # # apply subgrid temperature diffusion on markers
        # ---------------------------------------------------------------------
        # for m = 1:1:marknum
        #     # ~50 lines MATLAB
        # end


# Sat 16
        # ---------------------------------------------------------------------
        # ---------------------------------------------------------------------
        # # compute DTsubgrid
        # ---------------------------------------------------------------------
        # compute_DT_subgrid!(sp, dp, interp_arrays)


# Sat 16
        # ---------------------------------------------------------------------
        # ---------------------------------------------------------------------
        # # interpolate DT to markers
        # ---------------------------------------------------------------------
        # for m = 1:1:marknum
        #     # ~30 lines MATLAB
        # end


# Sat 16
        # ---------------------------------------------------------------------
        # ---------------------------------------------------------------------
        # # update porosity on markers
        # ---------------------------------------------------------------------
        # for m = 1:1:marknum
        #     # ~30 lines MATLAB
        # end


# Mon 18
        # ---------------------------------------------------------------------
        # # compute fluid velocity in P nodes including boundary conditions
        # ---------------------------------------------------------------------
        # compute_v_fluid_p_nodes(sp, dp, interp_arrays)


# Mon 18
        # ---------------------------------------------------------------------
        # # compute velocity in P nodes
        # ---------------------------------------------------------------------
        # compute_v_p_nodes!(sp, dp, interp_arrays)


# Mon 18
        # ---------------------------------------------------------------------
        # # compute rotation rate in basic nodes
        # ---------------------------------------------------------------------
        # compute_ω_basic_nodes!(sp, dp, interp_arrays)


# Tue 19
        # ---------------------------------------------------------------------
        # # move markers with RK4
        # ---------------------------------------------------------------------
        # for m = 1:1:marknum
        #     # ~300 lines MATLAB
        # end


# Wed 20
        # ---------------------------------------------------------------------
        # # backtrack P nodes: Ptotal with RK4
        # ---------------------------------------------------------------------


# Thu 21
        # ---------------------------------------------------------------------
        # # backtrack P nodes: Pfluid with RK1/2/3
        # ---------------------------------------------------------------------


# Fri 22
        # ---------------------------------------------------------------------
        # # replenish sparse areas with additional markers
        # ---------------------------------------------------------------------
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
    xsize=140000.0,
    ysize=140000.0,
    rplanet=50000.0,
    rcrust=48000.0,
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
module HydrologyPlanetesimals

using Base.Threads
using SparseArrays
using MAT
using DocStringExtensions
using Parameters
using StaticArrays

export StaticParameters, simulation_loop


""""
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
    xsize::Float64 = 140000
    "vertical model size [m]"
    ysize::Float64 = 140000
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
    rplanet::Int64 = 50000
    "crust radius [m]"
    rcrust::Int64 = 48000
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
    startmarknum::Int64 = Nxm * Nym
    # physical constants
    "gravitational constant [m^3*kg^-1*s^-2]"
    G::Float64 = 6.672e-11
    "scaled pressure"    
    pscale::Float64 = 1e+23 / dx
    # materials properties:              planet      crust       space
    "solid Density [kg/m^3]"
    rhosolidm::SVector{3, Float64}   = [ 3300.0    , 3300.0    ,    1.0    ]
    "fluid density [kg/m^3]"	
    rhofluidm::SVector{3, Float64}   = [ 7000.0    , 7000.0    , 1000.0    ]
    "solid viscosity [Pa*s]"
    etasolidm::SVector{3, Float64}   = [    1.0e+16,    1.0e+16,    1.0e+14]
    "molten solid viscosity [Pa*s]"
    etasolidmm::SVector{3, Float64}  = [    1.0e+14,    1.0e+14,    1.0e+14]
    "fluid viscosity [Pa*s]"
    etafluidm::SVector{3, Float64}   = [    1.0e-02,    1.0e-02,    1.0e+12]
    "molten fluid viscosity [Pa*s]"
    etafluidmm::SVector{3, Float64}  = [    1.0e-02,    1.0e-02,    1.0e+12]
    "solid volumetric heat capacity [kg/m^3]"
    rhocpsolidm::SVector{3, Float64} = [    3.3e+06,    3.3e+06,    3.0e+06]
    "fluid volumetric heat capacity [kg/m^3]"
    rhocpfluidm::SVector{3, Float64} = [    7.0e+06,    7.0e+06,    3.0e+06]
    "solid thermal expansion [1/K]"
    alphasolidm::SVector{3, Float64} = [    3.0e-05,    3.0e-05,    0.0    ]
    "fluid thermal expansion [1/K]"
    alphafluidm::SVector{3, Float64} = [    5.0e-05,    5.0e-05,    0.0    ]
    "solid thermal conductivity [W/m/K]"
    ksolidm::SVector{3, Float64}     = [    3.0    ,    3.0    , 3000.0    ]
    "fluid thermal conductivity [W/m/K]"
    kfluidm::SVector{3, Float64}     = [   50.0    ,   50.0    , 3000.0    ]
    # "solid radiogenic heat production [W/m^3]"
    # hrsolidm::Array{Float64}    = [    0.0    ,    0.0    ,    0.0    ]
    # "fluid radiogenic heat production [W/m^3]"
    # hrfluidm::Array{Float64}    = [    0.0    ,    0.0    ,    0.0    ]
    "solid shear modulus [Pa]"
    gggsolidm::SVector{3, Float64}   = [    1.0e+10,    1.0e+10,    1.0e+10]
    "solid friction coefficient"
    frictsolidm::SVector{3, Float64} = [    0.6    ,    0.6    ,    0.0    ]
    "solid compressive strength [Pa]"
    cohessolidm::SVector{3, Float64} = [    1.0e+08,    1.0e+08,    1.0e+08]
    "solid tensile strength [Pa]"
    tenssolidm ::SVector{3, Float64} = [    6.0e+07,    6.0e+07,    6.0e+07]
    "standard permeability [m^2]"
    kphim0::SVector{3, Float64}      = [    1.0e-13,    1.0e-13,    1.0e-17]
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
    starttime::Float64 = 1e6 * yearlength 
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
    "Yielding error of nodes"
    YERRNOD::Array{Float64} = zeros(1, nplast) 
    "Weight for old viscosity"
    etawt::Float64 = 0 
    "max porosity ratio change per time step"
    dphimax::Float64 = 0.01
    "starting timestep"
    startstep::Int64 = 1
    "number of timesteps to run"
    nsteps::Int64 = 30000 
end


"""
Marker properties: Fixed and calculated during timestepping

$(TYPEDFIELDS)
"""
@with_kw struct MarkerArrays
    # original marker properties 
    "horizontal coordinates [m]"
    xm::Vector{Float64}
    "vertical coordinates [m]"
    ym::Vector{Float64}
    "material type"
    tm::Vector{Int8}
    "marker temperature [K]"
    tkm::Vector{Float64}
    "SIGMA'xx [Pa]"
    sxxm::Vector{Float64}
    "SIGMAxy [Pa]"
    sxym::Vector{Float64}
    "Visco-plastic viscosity [Pa]"
    etavpm::Vector{Float64}
    "Marker porosity"
    phim::Vector{Float64}
    # fixed marker properties used during timestepping calculations
    # RMK: omitted, as only used once - reconsider?
    # marker properties calculated during timestepping
    "kphim"
    kphim::Vector{Float64}
    "rhototalm"
    rhototalm::Vector{Float64}
    "rhocptotalm"
    rhocptotalm::Vector{Float64}
    "etatotalm"
    etatotalm::Vector{Float64}
    "hrtotalm"
    hrtotalm::Vector{Float64}
    "ktotalm"
    ktotalm::Vector{Float64}
    "gggtotalm"
    gggtotalm::Vector{Float64}
    "fricttotalm"
    fricttotalm::Vector{Float64}
    "cohestotalm"
    cohestotalm::Vector{Float64}
    "tenstotalm"
    tenstotalm::Vector{Float64}
    "etafluidcur"
    etafluidcur::Vector{Float64}
    "rhofluidcur"
    rhofluidcur::Vector{Float64}
    "inner constructor"
    MarkerArrays(marknum) = new(
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum),
        zeros(marknum)
    )
end


"""
Initialize markers according to model parameters

$(SIGNATURES)

# Details

    - ma: arrays containing marker properties
    - sp: static simulation parameters
    - dp: dynamic simulation parameters

# Returns

    - nothing
"""
function define_markers!(
    xm,
    ym,
    tm,
    tkm,
    phim,
    etavpm,
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
    phimin,
    etasolidm,
    phim0 = sp

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
            if rmark > rcrust
                # crust
                tm[m] = 2
            else
                # mantle
                tm[m] = 1
            end
            # temperature
            tkm[m] = 300
            # porosity
            phim[m] = phim0 * (1.0 + (rand()-0.5))
            # matrix viscosity
            etavpm[m] = etasolidm[tm[m]] # *exp(-28*phim[m])
        else
            # sticky space ("air") [to have internal free surface]
            # space
            tm[m] = 3
            # temperature
            tkm[m] = 273
            # porosity
            phim[m] = phimin
            # matrix viscosity
            etavpm[m] = etasolidm[tm[m]]
        end
    end
    return nothing
end


"""
Compute static marker properties which stay constant during simulation.
Runs once during initialization.

$(SIGNATURES)

# Details 

    - m: marker counter of marker whose static properties are to be computed
    - sp: static simulation parameters

# Returns

    - nothing
"""	
function compute_static_marker_properties!(
    m,
    tm,
    rhototalm,
    rhocptotalm,
    etatotalm,
    hrtotalm,
    ktotalm,
    gggtotalm,
    fricttotalm,
    cohestotalm,
    tenstotalm,
    etafluidcur,
    rhofluidcur,
    hrsolidm,
    sp::StaticParameters
)
    @unpack rhosolidm,
    rhocpsolidm,
    etasolidm,
    ksolidm,
    gggsolidm,
    frictsolidm,
    cohessolidm,
    tenssolidm,
    etafluidm,
    rhofluidm = sp

    # static secondary marker properties
    if tm[m] < 3
        # rocks
        # pass
    else
        # air
        rhototalm[m] = rhosolidm[tm[m]]
        rhocptotalm[m] = rhocpsolidm[tm[m]]
        etatotalm[m] = etasolidm[tm[m]]
        hrtotalm[m] = hrsolidm[tm[m]]
        ktotalm[m] = ksolidm[tm[m]]
    end
    # common for rocks and air
    inv_gggtotalm[m] = inv(gggsolidm[tm[m]])
    fricttotalm[m] = frictsolidm[tm[m]]
    cohestotalm[m] = cohessolidm[tm[m]]
    tenstotalm[m] = tenssolidm[tm[m]]
    etafluidcur[m] = etafluidm[tm[m]]
    rhofluidcur[m] = rhofluidm[tm[m]]

    return nothing
end


"""
Compute dynamic marker properties.
Runs every simulation loop (i.e. time step).

$(SIGNATURES)

# Detail

    - m: marker index of marker whose parameters are to be computed
    - ma: marker arrays containing marker properties to be updated
    - sp: static simulation parameters
    - dp: dynamic simulation parameters

# Returns
    - nothing
"""
function compute_dynamic_marker_params!(
    m::Int64,
    ma::MarkerArrays,
    sp::StaticParameters,
    dp::DynamicParameters
)
    # @unpack tm,
        # tkm,
        # phim,
        # rhototalm,
        # rhocptotalm,
        # etatotalm,
        # hrtotalm,
        # ktotalm,
        # kphim = ma
    # @unpack hrsolidm, hrfluidm = dp
    @unpack rhosolidm,
        rhofluidm,
        rhocpsolidm,
        rhocpfluidm,
        tmiron,
        tmsilicate,
        etamin,
        etasolidmm,
        etasolidm,
        etafluidmm,
        etafluidm,
        kphim0,
        phim0,
        ksolidm,
        kfluidm  = sp
    
    if tm[m] < 3
        # rocks
        rhototalm[m] = total(rhosolidm[tm[m]], rhofluidm[tm[m]], phim[m])
        rhocptotalm[m] = total(rhocpsolidm[tm[m]], rhocpfluidm[tm[m]], phim[m])
        etatotalm[m] = etatotal_rock(
            tkm[m],
            tmsilicate,
            tmiron,
            etamin,
            etasolidm[tm[m]],
            etasolidmm[tm[m]],
            etafluidm[tm[m]],
            etafluidmm[tm[m]]
            )
        hrtotalm[m] = total(hrsolidm[tm[m]], hrfluidm[tm[m]], phim[m])
        ktotalm[m] = ktotal(ksolidm[tm[m]], kfluidm[tm[m]], phim[m])
        
    else
        # air
    
    end
    # common for rocks and air
    kphim[m] = kphi(kphim0[tm[m]], phim0, phim[m])

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
Compute temperature-dependent total viscosity for iron-containing silicate rock.

$(SIGNATURES)

# Details

    - tk: temperature [K]
    - tmsilicate: melting temperature of silicate in [K]
    - tmiron: melting temperature of iron in K [K]
    - etamin: minimum viscosity [Pa s]
    - etasolidm: solid viscosity [Pa s]
    - etasolidmm: molten solid viscosity [Pa s]
    - etafluidm: fluid viscosity [Pa s]
    - etafluidmm: molten fluid viscosity [Pa s]

# Returns

    - etatotal: temperature-dependent total viscosity [Pa s]
"""	
function etatotal_rock(
    tk,
    tmsilicate,
    tmiron,
    etamin,
    etasolidm,
    etasolidmm,
    etafluidm,
    etafluidmm
    )
    return max(
        etamin,
        tk > tmsilicate ? etasolidmm : etasolidm,
        tk > tmiron ? etafluidmm : etafluidm
        )
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
function calculate_radioactive_heating(
    sp::StaticParameters, dp::DynamicParameters
    )
    @unpack hr_al, f_al, ratio_al, E_al, tau_al, hr_fe, f_fe, ratio_fe, E_fe,
        tau_fe, rhosolidm, rhofluidm = sp
    @unpack timesum = dp
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
    @inbounds j = trunc(Int, (x - x_axis[1]) / dx) + 1
    @inbounds i = trunc(Int, (y - y_axis[1]) / dy) + 1
    j = min(max(j, jmin), jmax)
    i = min(max(i, imin), imax)
    # ij = [SVector(i, j), SVector(i+1, j), SVector(i, j+1), SVector(i+1, j+1)]
    @inbounds dxmj = x - x_axis[j]
    @inbounds dymi = y - y_axis[i]
    weights = SVector(
        (1.0-dymi/dy) * (1.0-dxmj/dx),
        (dymi/dy) * (1.0-dxmj/dx),
        (1.0-dymi/dy) * (dxmj/dx),
        (dymi/dy) * (dxmj/dx)
        )
    return i, j, weights
end


"""
Interpolate marker properties to basic nodes.

$(SIGNATURES)

# Details

    - m: index of marker whose properties are to be interpolated to nodes
    - mrk: arrays containing all marker properties
    - i: top node index of marker m on basic grid
    - j: left node index of marker m on basic grid
    - weights: bilinear interpolation weights to four neighbor nodes of marker m
    - ETA0SUM: viscous viscosity array interpolated to basic nodes
    - ETASUM: viscoplastic viscosity array interpolated to basic nodes
    - GGGSUM: shear modulus array interpolated to basic nodes
    - SXYSUM: σxy shear stress array interpolated to basic nodes
    - COHSUM: copmressive strength array interpolated to basic nodes
    - TENSUM: tensile strength array interpolated to basic nodes
    - FRISUM: friction array interpolated to basic nodes
    - WTSUM: weight array for bilinear interpolation to basic nodes

# Returns

    -nothing
"""
function interpolate_basic_nodes!(
    m,
    mrk,
    i,
    j,
    weights,
    ETA0SUM,
    ETASUM,
    GGGSUM,
    SXYSUM,
    COHSUM,
    TENSUM,
    FRISUM,
    WTSUM
)
    ETA0SUM[i, j, threadid()] += mrk.etatotalm[m] * weights[1]
    ETA0SUM[i+1, j, threadid()] += mrk.etatotalm[m] * weights[2]
    ETA0SUM[i, j+1, threadid()] += mrk.etatotalm[m] * weights[3]
    ETA0SUM[i+1, j+1, threadid()] += mrk.etatotalm[m] * weights[4]

    ETASUM[i, j, threadid()] += mrk.etavpm[m] * weights[1]
    ETASUM[i+1, j, threadid()] += mrk.etavpm[m] * weights[2]
    ETASUM[i, j+1, threadid()] += mrk.etavpm[m] * weights[3]
    ETASUM[i+1, j+1, threadid()] += mrk.etavpm[m] * weights[4]

    GGGSUM[i, j, threadid()] += inv(mrk.gggtotalm[m]) * weights[1]
    GGGSUM[i+1, j, threadid()] += inv(mrk.gggtotalm[m]) * weights[2]
    GGGSUM[i, j+1, threadid()] += inv(mrk.gggtotalm[m]) * weights[3]
    GGGSUM[i+1, j+1, threadid()] += inv(mrk.gggtotalm[m]) * weights[4]

    SXYSUM[i, j, threadid()] += mrk.sxym[m] * weights[1]
    SXYSUM[i+1, j, threadid()] += mrk.sxym[m] * weights[2]
    SXYSUM[i, j+1, threadid()] += mrk.sxym[m] * weights[3]
    SXYSUM[i+1, j+1, threadid()] += mrk.sxym[m] * weights[4]

    COHSUM[i, j, threadid()] += mrk.cohestotalm[m] * weights[1]
    COHSUM[i+1, j, threadid()] += mrk.cohestotalm[m] * weights[2]
    COHSUM[i, j+1, threadid()] += mrk.cohestotalm[m] * weights[3]
    COHSUM[i+1, j+1, threadid()] += mrk.cohestotalm[m] * weights[4]

    TENSUM[i, j, threadid()] += mrk.tenstotalm[m] * weights[1]
    TENSUM[i+1, j, threadid()] += mrk.tenstotalm[m] * weights[2]
    TENSUM[i, j+1, threadid()] += mrk.tenstotalm[m] * weights[3]
    TENSUM[i+1, j+1, threadid()] += mrk.tenstotalm[m] * weights[4]

    FRISUM[i, j, threadid()] += mrk.fricttotalm[m] * weights[1]
    FRISUM[i+1, j, threadid()] += mrk.fricttotalm[m] * weights[2]
    FRISUM[i, j+1, threadid()] += mrk.fricttotalm[m] * weights[3]
    FRISUM[i+1, j+1, threadid()] += mrk.fricttotalm[m] * weights[4]

    WTSUM[i, j, threadid()] += weights[1]
    WTSUM[i+1, j, threadid()] += weights[2]
    WTSUM[i, j+1, threadid()] += weights[3]
    WTSUM[i+1, j+1, threadid()] += weights[4]

    return nothing
end


"""
Interpolate marker properties to Vx nodes.

$(SIGNATURES)

# Details

    - m: index of marker whose properties are to be interpolated to nodes
    - mrk: arrays containing all marker properties
    - i: top node index of marker m on Vx grid
    - j: left node index of marker m on Vx grid
    - weights: bilinear interpolation weights to four neighbor nodes of marker m
    - RHOXSUM: density array interpolated to Vx nodes
    - RHOFXSUM: fluid density array interpolated to Vx nodes
    - KXSUM: thermal conductivity array interpolated to Vx nodes
    - PHIXSUM: porosity array interpolated to Vx nodes
    - RXSUM: ηfluid/kϕ array interpolated to Vx nodes
    - WTXSUM: weight array for bilinear interpolation to Vx nodes

# Returns

    -nothing
"""
function interpolate_vx_nodes!(
    m,
    mrk,
    i,
    j,
    weights,
    RHOXSUM,
    RHOFXSUM,
    KXSUM,
    PHIXSUM,
    RXSUM,
    WTXSUM
)
    RHOXSUM[i, j, threadid()] += mrk.rhototalm[m] * weights[1]
    RHOXSUM[i+1, j, threadid()] += mrk.rhototalm[m] * weights[2]
    RHOXSUM[i, j+1, threadid()] += mrk.rhototalm[m] * weights[3]
    RHOXSUM[i+1, j+1, threadid()] += mrk.rhototalm[m] * weights[4]

    RHOFXSUM[i, j, threadid()] += mrk.rhofluidcur[m] * weights[1]
    RHOFXSUM[i+1, j, threadid()] += mrk.rhofluidcur[m] * weights[2]
    RHOFXSUM[i, j+1, threadid()] += mrk.rhofluidcur[m] * weights[3]
    RHOFXSUM[i+1, j+1, threadid()] += mrk.rhofluidcur[m] * weights[4]

    KXSUM[i, j, threadid()] += mrk.ktotalm[m] * weights[1]
    KXSUM[i+1, j, threadid()] += mrk.ktotalm[m] * weights[2]
    KXSUM[i, j+1, threadid()] += mrk.ktotalm[m] * weights[3]
    KXSUM[i+1, j+1, threadid()] += mrk.ktotalm[m] * weights[4]

    PHIXSUM[i, j, threadid()] += mrk.phim[m] * weights[1]
    PHIXSUM[i+1, j, threadid()] += mrk.phim[m] * weights[2]
    PHIXSUM[i, j+1, threadid()] += mrk.phim[m] * weights[3]
    PHIXSUM[i+1, j+1, threadid()] += mrk.phim[m] * weights[4]

    RXSUM[i, j, threadid()] += mrk.etafluidcur[m] / mrk.kphim[m] * weights[1]
    RXSUM[i+1, j, threadid()] += mrk.etafluidcur[m] / mrk.kphim[m] * weights[2]
    RXSUM[i, j+1, threadid()] += mrk.etafluidcur[m] / mrk.kphim[m] * weights[3]
    RXSUM[i+1, j+1, threadid()] += mrk.etafluidcur[m] /mrk.kphim[m] * weights[4]

    WTXSUM[i, j, threadid()] += weights[1]
    WTXSUM[i+1, j, threadid()] += weights[2]
    WTXSUM[i, j+1, threadid()] += weights[3]
    WTXSUM[i+1, j+1, threadid()] += weights[4]

    return nothing 
end


"""
Interpolate marker properties to Vy nodes.

$(SIGNATURES)

# Details

    - m: index of marker whose properties are to be interpolated to nodes
    - mrk: arrays containing all marker properties
    - i: top node index of marker m on Vy grid
    - j: left node index of marker m on Vy grid
    - weights: bilinear interpolation weights to four neighbor nodes of marker m
    - RHOYSUM: density array interpolated to Vy nodes
    - RHOFYSUM: fluid density array interpolated to Vy nodes
    - KYSUM: thermal conductivity array interpolated to Vy nodes
    - PHIYSUM: porosity array interpolated to Vy nodes
    - RYSUM: ηfluid/kϕ array interpolated to Vy nodes
    - WTYSUM: weight array for bilinear interpolation to Vy nodes

# Returns

    -nothing
"""
function interpolate_vy_nodes!(
    m,
    mrk,
    i,
    j,
    weights,
    RHOYSUM,
    RHOFYSUM,
    KYSUM,
    PHIYSUM,
    RYSUM,
    WTYSUM
)
    RHOYSUM[i, j, threadid()] += mrk.rhototalm[m] * weights[1]
    RHOYSUM[i+1, j, threadid()] += mrk.rhototalm[m] * weights[2]
    RHOYSUM[i, j+1, threadid()] += mrk.rhototalm[m] * weights[3]
    RHOYSUM[i+1, j+1, threadid()] += mrk.rhototalm[m] * weights[4]

    RHOFYSUM[i, j, threadid()] += mrk.rhofluidcur[m] * weights[1]
    RHOFYSUM[i+1, j, threadid()] += mrk.rhofluidcur[m] * weights[2]
    RHOFYSUM[i, j+1, threadid()] += mrk.rhofluidcur[m] * weights[3]
    RHOFYSUM[i+1, j+1, threadid()] += mrk.rhofluidcur[m] * weights[4]

    KYSUM[i, j, threadid()] += mrk.ktotalm[m] * weights[1]
    KYSUM[i+1, j, threadid()] += mrk.ktotalm[m] * weights[2]
    KYSUM[i, j+1, threadid()] += mrk.ktotalm[m] * weights[3]
    KYSUM[i+1, j+1, threadid()] += mrk.ktotalm[m] * weights[4]

    PHIYSUM[i, j, threadid()] += mrk.phim[m] * weights[1]
    PHIYSUM[i+1, j, threadid()] += mrk.phim[m] * weights[2]
    PHIYSUM[i, j+1, threadid()] += mrk.phim[m] * weights[3]
    PHIYSUM[i+1, j+1, threadid()] += mrk.phim[m] * weights[4]

    RYSUM[i, j, threadid()] += mrk.etafluidcur[m] / mrk.kphim[m] * weights[1]
    RYSUM[i+1, j, threadid()] += mrk.etafluidcur[m] / mrk.kphim[m] * weights[2]
    RYSUM[i, j+1, threadid()] += mrk.etafluidcur[m] / mrk.kphim[m] * weights[3]
    RYSUM[i+1, j+1, threadid()] += mrk.etafluidcur[m] /mrk.kphim[m] * weights[4]

    WTYSUM[i, j, threadid()] += weights[1]
    WTYSUM[i+1, j, threadid()] += weights[2]
    WTYSUM[i, j+1, threadid()] += weights[3]
    WTYSUM[i+1, j+1, threadid()] += weights[4]

    return nothing
end


"""
Interpolate marker properties to P nodes.

$(SIGNATURES)

# Details

    - m: index of marker whose properties are to be interpolated to nodes
    - mrk: arrays containing all marker properties
    - i: top node index of marker m on P grid
    - j: left node index of marker m on P grid
    - weights: bilinear interpolation weights to four neighbor nodes of marker m
    - GGGPSUM: shear modulus array interpolated to P nodes
    - SXXSUM: σ'xx array interpolated to P nodes
    - RHOSUM: density array interpolated to P nodes
    - RHOCPSUM: volumetric heat capacity array interpolated to P nodes
    - ALPHASUM: thermal expansion array interpolated to P nodes
    - ALPHAFSUM: fluid thermal expansion array interpolated to P nodes
    - HRSUM: radioactive heating array interpolated to P nodes
    - TKSUM: temperature array interpolated to P nodes
    - PHISUM: porosity array interpolated to P nodes
    - WTPSUM: weight array for bilinear interpolation to P nodes

# Returns

    -nothing
"""
function interpolate_p_nodes!(
    m,
    mrk,
    i,
    j,
    weights,
    GGGPSUM,
    SXXSUM,
    RHOSUM,
    RHOCPSUM,
    ALPHASUM,
    ALPHAFSUM,
    HRSUM,
    TKSUM,
    PHISUM,
    WTPSUM
)
    GGGPSUM[i, j, threadid()] += 1.0 / mrk.gggtotalm[m] * weights[1]
    GGGPSUM[i+1, j, threadid()] += 1.0 / mrk.gggtotalm[m] * weights[2]
    GGGPSUM[i, j+1, threadid()] += 1.0 / mrk.gggtotalm[m] * weights[3]
    GGGPSUM[i+1, j+1, threadid()] += 1.0 / mrk.gggtotalm[m] * weights[4]

    SXXSUM[i, j, threadid()] += mrk.sxxm[m] * weights[1]
    SXXSUM[i+1, j, threadid()] += mrk.sxxm[m] * weights[2]
    SXXSUM[i, j+1, threadid()] += mrk.sxxm[m] * weights[3]
    SXXSUM[i+1, j+1, threadid()] += mrk.sxxm[m] * weights[4]

    RHOSUM[i, j, threadid()] += mrk.rhototalm[m] * weights[1]
    RHOSUM[i+1, j, threadid()] += mrk.rhototalm[m] * weights[2]
    RHOSUM[i, j+1, threadid()] += mrk.rhototalm[m] * weights[3]
    RHOSUM[i+1, j+1, threadid()] += mrk.rhototalm[m] * weights[4]

    RHOCPSUM[i, j, threadid()] += mrk.rhocptotalm[m] * weights[1]
    RHOCPSUM[i+1, j, threadid()] += mrk.rhocptotalm[m] * weights[2]
    RHOCPSUM[i, j+1, threadid()] += mrk.rhocptotalm[m] * weights[3]
    RHOCPSUM[i+1, j+1, threadid()] += mrk.rhocptotalm[m] * weights[4]

    ALPHASUM[i, j, threadid()] += mrk.alphasolidm[m] * weights[1]
    ALPHASUM[i+1, j, threadid()] += mrk.alphasolidm[m] * weights[2]
    ALPHASUM[i, j+1, threadid()] += mrk.alphasolidm[m] * weights[3]
    ALPHASUM[i+1, j+1, threadid()] += mrk.alphasolidm[m] * weights[4]

    ALPHAFSUM[i, j, threadid()] += mrk.alphafluidm[m] * weights[1]
    ALPHAFSUM[i+1, j, threadid()] += mrk.alphafluidm[m] * weights[2]
    ALPHAFSUM[i, j+1, threadid()] += mrk.alphafluidm[m] * weights[3]
    ALPHAFSUM[i+1, j+1, threadid()] += mrk.alphafluidm[m] * weights[4]

    HRSUM[i, j, threadid()] += mrk.hrtotalm[m] * weights[1]
    HRSUM[i+1, j, threadid()] += mrk.hrtotalm[m] * weights[2]
    HRSUM[i, j+1, threadid()] += mrk.hrtotalm[m] * weights[3]
    HRSUM[i+1, j+1, threadid()] += mrk.hrtotalm[m] * weights[4]

    TKSUM[i, j, threadid()] += mrk.tkm[m] * mrk.rhocptotalm[m] * weights[1]
    TKSUM[i+1, j, threadid()] += mrk.tkm[m] * mrk.rhocptotalm[m] * weights[2]
    TKSUM[i, j+1, threadid()] += mrk.tkm[m] * mrk.rhocptotalm[m] * weights[3]
    TKSUM[i+1, j+1, threadid()] += mrk.tkm[m] * mrk.rhocptotalm[m] * weights[4]

    PHISUM[i, j, threadid()] += mrk.phim[m] * weights[1]
    PHISUM[i+1, j, threadid()] += mrk.phim[m] * weights[2]
    PHISUM[i, j+1, threadid()] += mrk.phim[m] * weights[3]
    PHISUM[i+1, j+1, threadid()] += mrk.phim[m] * weights[4]

    WTPSUM[i, j, threadid()] += weights[1]
    WTPSUM[i+1, j, threadid()] += weights[2]
    WTPSUM[i, j+1, threadid()] += weights[3]
    WTPSUM[i+1, j+1, threadid()] += weights[4]

    return nothing
end


"""
Main simulation loop: run calculations with timestepping.

$(SIGNATURES)

# Detail

    - markers: arrays containing all marker properties
    - sp: static simulation parameters

# Returns
    
    - nothing
"""
function simulation_loop(markers::MarkerArrays, sp::StaticParameters)
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
    dtelastic,
    startstep,
    nsteps,
    starttime, 
    endtime,
    startmarknum = sp

    
    # -------------------------------------------------------------------------
    # set up dynamic simulation parameters from given static parameters
    # -------------------------------------------------------------------------
    # timestep counter (current), init to startstep
    timestep = startstep
    # computational timestep (current), init to dtelastic [s]
    dt = dtelastic
    # time sum (current), init to starttime [s]
    timesum = starttime
    # current number of markers, init to startmarknum
    marknum = startmarknum
    # radiogenic heat production solid phase, init to zero
    hrsolidm = SVector{3, Float64}(zeros(3))
    # radiogenic heat production fluid phase, init to zero
    hrfluidm = SVector{3, Float64}(zeros(3))
   

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
    # plastic yielding mark, 1=yes,0=no
    YNY = zeros(Int8, Ny, Nx)

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
    tm = zeros(Float64, marknum)
    # marker temperature [K]
    tkm = zeros(Float64, marknum)
    # marker σ′xx [Pa]
    sxxm = zeros(Float64, marknum)
    # marker σxy [Pa]
    sxym = zeros(Float64, marknum)
    # marker viscoplastic viscosity [Pa]
    etavpm = zeros(Float64, marknum)
    # marker porosity ϕ
    phim = zeros(Float64, marknum)
    # # secondary marker arrays: derived properties calculated during time steps
    # # kphim
    # kphim = zeros(Float64, marknum)
    # # rhototalm
    # rhototalm = zeros(Float64, marknum)
    # # rhocptotalm
    # rhocptotalm = zeros(Float64, marknum)
    # # etatotalm
    # etatotalm = zeros(Float64, marknum)
    # # hrtotalm
    # hrtotalm = zeros(Float64, marknum)
    # # ktotalm
    # ktotalm = zeros(Float64, marknum)
    # # 1 / gggtotalm
    # inv_gggtotalm = zeros(Float64, marknum)
    # # fricttotalm
    # fricttotalm = zeros(Float64, marknum)
    # # cohestotalm
    # cohestotalm = zeros(Float64, marknum)
    # # tenstotalm
    # tenstotalm = zeros(Float64, marknum)
    # # etafluidcur
    # etafluidcur = zeros(Float64, marknum)
    # # rhofluidcur
    # rhofluidcur = zeros(Float64, marknum)

    # define markers: coordinates, temperature, and material type    
    define_markers!(xm, ym, tm, tkm, phim, etavpm, sp)


    # -------------------------------------------------------------------------
    # set up of matrices for global gravity/thermal/hydromechanical solutions
    # -------------------------------------------------------------------------
    # hydromechanical solution: LHS coefficient matrix
    L = SparseMatrixCSC{Float64, Int64}(Nx1*Ny1*6, Nx1*Ny1*6)
    # hydromechanical solution: RHS Vector
    R = zeros(Float64, Nx1*Ny1*6)
    # thermal solution: LHS coefficient matrix
    LT = SparseMatrixCSC{Float64, Int64}(Nx1*Ny1, Nx1*Ny1)
    # thermal solution: RHS Vector
    RT = zeros(Float64, Nx1*Ny1)
    # gravity solution: LHS coefficient matrix
    LP = SparseMatrixCSC{Float64, Int64}(Nx1*Ny1, Nx1*Ny1)
    # gravity solution: RHS Vector
    RP = zeros(Float64, Nx1*Ny1)


    # -------------------------------------------------------------------------
    # iterate timesteps   
    # -------------------------------------------------------------------------
    for timestep = startstep:1:100
    # for timestep = startstep:1:nsteps

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
        TKSUM = zeros(Ny1, Nx1, nthreads())
        PHISUM = zeros(Ny1, Nx1, nthreads())
        WTPSUM = zeros(Ny1, Nx1, nthreads())


        # ---------------------------------------------------------------------
        # calculate radioactive heating
        # ---------------------------------------------------------------------
        hrsolidm, hrfluidm = calculate_radioactive_heating(sp, dp)

        
        # ---------------------------------------------------------------------
        # computer marker properties and interpolate to staggered grid nodes
        # ---------------------------------------------------------------------
        @threads for m = 1:1:marknum
            # compute marker properties 
            if tm[m] < 3
                # rocks
                rhototalm[m] = total(rhosolidm[tm[m]], rhofluidm[tm[m]], phim[m])
                rhocptotalm[m] = total(rhocpsolidm[tm[m]], rhocpfluidm[tm[m]], phim[m])
                etatotalm[m] = etatotal_rock(
                    tkm[m],
                    tmsilicate,
                    tmiron,
                    etamin,
                    etasolidm[tm[m]],
                    etasolidmm[tm[m]],
                    etafluidm[tm[m]],
                    etafluidmm[tm[m]]
                    )
                hrtotalm[m] = total(hrsolidm[tm[m]], hrfluidm[tm[m]], phim[m])
                ktotalm[m] = ktotal(ksolidm[tm[m]], kfluidm[tm[m]], phim[m])
                
            else
                # air
            
            end
            # common for rocks and air
            kphim[m] = kphi(kphim0[tm[m]], phim0, phim[m])
            

            # interpolate marker properties to basic nodes
            i, j, weights = fix_weights(
                markers.xm[m],
                markers.ym[m],
                x,
                y,
                dx,
                dy,
                jmin_basic,
                jmax_basic,
                imin_basic,
                imax_basic
            )
            interpolate_basic_nodes!(
                m,
                markers,
                i,
                j,
                weights,
                ETA0SUM,
                ETASUM,
                GGGSUM,
                SXYSUM,
                COHSUM,
                TENSUM,
                FRISUM,
                WTSUM
            )

            # # interpolate marker properties to Vx nodes
            i, j, weights = fix_weights(
                markers.xm[m],
                markers.ym[m],
                xvx,
                yvx,
                dx,
                dy,
                jmin_vx,
                jmax_vx,
                imin_vx,
                imax_vx
            )
            interpolate_vx_nodes!(
                m,
                markers,
                i,
                j,
                weights,
                RHOXSUM,
                RHOFXSUM,
                KXSUM,
                PHIXSUM,
                RXSUM,
                WTXSUM
            )

            # interpolate marker properties to Vy nodes
            i, j, weights = fix_weights(
                markers.xm[m],
                markers.ym[m],
                xvy,
                yvy,
                dx,
                dy,
                jmin_vy,
                jmax_vy,
                imin_vy,
                imax_vy
            )
            interpolate_vy_nodes!(
                m,
                markers,
                i,
                j,
                weights,
                RHOYSUM,
                RHOFYSUM,
                KYSUM,
                PHIYSUM,
                RYSUM,
                WTYSUM
            )

            # interpolate marker properties to P nodes
            i, j, weights = fix_weights(
                markers.xm[m],
                markers.ym[m],
                xp,
                yp,
                dx,
                dy,
                jmin_p,
                jmax_p,
                imin_p,
                imax_p
            )
            interpolate_p_nodes!(
                m,
                markers,
                i,
                j,
                weights,
                GGGPSUM,
                SXXSUM,
                RHOSUM,
                RHOCPSUM,
                ALPHASUM,
                ALPHAFSUM,
                HRSUM,
                TKSUM,
                PHISUM,
                WTPSUM
            )
        end


        # reduce interpolation arrays
        # ETA = reduce(+, WTPSUM, dims=3)



        # ---------------------------------------------------------------------
        # compute physical properties of basic nodes
        # ---------------------------------------------------------------------
        # compute_properties_basic_nodes!(sp, dp, interp_arrays)


        # ---------------------------------------------------------------------
        # # compute physical properties of Vx nodes
        # ---------------------------------------------------------------------
        # compute_properties_vx_nodes!(sp, dp, interp_arrays)


        # ---------------------------------------------------------------------
        # # compute physical properties of Vy nodes
        # ---------------------------------------------------------------------
        # compute_properties_vy_nodes!(sp, dp, interp_arrays)


        # ---------------------------------------------------------------------
        # # compute physical properties of P nodes
        # ---------------------------------------------------------------------
        # compute_properties_p_nodes!(sp, dp, interp_arrays)


        # ---------------------------------------------------------------------
        # # applying thermal boundary conditions for interpolated temperature
        # ---------------------------------------------------------------------
        # apply_thermal_bc(sp, dp, tk1)


        # ---------------------------------------------------------------------
        # # compute gravity solution
        # ---------------------------------------------------------------------
        # compute_gravity_solution!(sp, dp, interp_arrays)


        # ---------------------------------------------------------------------
        # # compute gravitational acceleration
        # ---------------------------------------------------------------------
        # compute_grav_accel!(sp, dp, interp_arrays)


        # ---------------------------------------------------------------------
        # # probe increasing computational timestep
        # ---------------------------------------------------------------------
        # dt = min(dt*dtkoefup, dtelastic)


        # ---------------------------------------------------------------------
        # # perform plastic iterations
        # ---------------------------------------------------------------------
        # for iplast = 1:1:nplast
        #     # ~600 lines MATLAB
        # end


        # ---------------------------------------------------------------------
        # # interpolate updated viscoplastic viscosity to markers
        # ---------------------------------------------------------------------
        # for m = 1:1:marknum
        #     # ~50 lines MATLAB 
        # end


        # ---------------------------------------------------------------------
        # # apply subgrid stress diffusion to markers
        # ---------------------------------------------------------------------
        # for m = 1:1:marknum
        #     # ~100 lines MATLAB 
        # end


        # ---------------------------------------------------------------------
        # # compute DSXXsubgrid, DSXYsubgrid
        # ---------------------------------------------------------------------
        # compute_dsxx_dsxy_subgrids!(sp, dp, interp_arrays)


        # ---------------------------------------------------------------------
        # # interpolate DSXX, DSXY to markers
        # ---------------------------------------------------------------------
        # for m = 1:1:marknum
        #     # ~50 lines MATLAB 
        # end


        # ---------------------------------------------------------------------
        # # compute shear heating HS in P nodes
        # ---------------------------------------------------------------------
        # compute_HS_p_nodes!(sp, dp, interp_arrays)


        # ---------------------------------------------------------------------
        # # compute adiabatic heating HA in P nodes
        # ---------------------------------------------------------------------
        # compute_HA_p_nodes!(sp, dp, interp_arrays)


        # ---------------------------------------------------------------------
        # # perform thermal iterations
        # ---------------------------------------------------------------------
        # # ~100 lines MATLAB


        # ---------------------------------------------------------------------
        # # apply subgrid temperature diffusion on markers
        # ---------------------------------------------------------------------
        # for m = 1:1:marknum
        #     # ~50 lines MATLAB
        # end


        # ---------------------------------------------------------------------
        # # compute DTsubgrid
        # ---------------------------------------------------------------------
        # compute_DT_subgrid!(sp, dp, interp_arrays)


        # ---------------------------------------------------------------------
        # # interpolate DT to markers
        # ---------------------------------------------------------------------
        # for m = 1:1:marknum
        #     # ~30 lines MATLAB
        # end


        # ---------------------------------------------------------------------
        # # update porosity on markers
        # ---------------------------------------------------------------------
        # for m = 1:1:marknum
        #     # ~30 lines MATLAB
        # end


        # ---------------------------------------------------------------------
        # # compute fluid velocity in P nodes including boundary conditions
        # ---------------------------------------------------------------------
        # compute_v_fluid_p_nodes(sp, dp, interp_arrays)


        # ---------------------------------------------------------------------
        # # compute velocity in P nodes
        # ---------------------------------------------------------------------
        # compute_v_p_nodes!(sp, dp, interp_arrays)


        # ---------------------------------------------------------------------
        # # compute rotation rate in basic nodes
        # ---------------------------------------------------------------------
        # compute_ω_basic_nodes!(sp, dp, interp_arrays)


        # ---------------------------------------------------------------------
        # # move markers with RK4
        # ---------------------------------------------------------------------
        # for m = 1:1:marknum
        #     # ~300 lines MATLAB
        # end


        # ---------------------------------------------------------------------
        # # backtrack P nodes: Ptotal with RK4
        # ---------------------------------------------------------------------


        # ---------------------------------------------------------------------
        # # backtrack P nodes: Pfluid with RK1/2/3
        # ---------------------------------------------------------------------


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
        if timestep % 20 == 0
            println("timestep: ", timestep)
        end

        if timesum > endtime
            break
        end

    end # for timestep = startstep:1:nsteps
end # function simulation loop

end # module HydrologyPlanetesimals

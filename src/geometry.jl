
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
    - BETAPHI : bulk compresibility at P nodes [1/Pa]
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
# @timeit to "apply_insulating_boundary_conditions!" begin
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
# end # @timeit to "apply_insulating_boundary_conditions!"
    return nothing
end

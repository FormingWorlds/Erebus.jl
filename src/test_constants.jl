 # radioactive switches
 # radioactive heating from 26Al active
 const hr_al = true
 # radioactive heating from 60Fe active	
 const hr_fe = true
 # model size, geometry, and resolution
 # horizontal model size [m]
 const xsize = 14_000.0
 # vertical model size [m]
 const ysize = 14_000.0
 # horizontal center of model
 const xcenter = xsize / 2
 # vertical center of model
 const ycenter = ysize / 2  
 # basic grid resolution in x direction (horizontal)
 const Nx = 15
 # basic grid resolution in y direction (vertical)	
 const Ny = 15
 # Vx, Vy, P grid resolution in x direction (horizontal)
 const Nx1 = Nx + 1
 # Vx/Vy/P grid resolution in y direction (vertical)
 const Ny1 = Ny + 1
 # horizontal grid step [m]
 const dx = xsize / (Nx-1)
 # vertical grid step [m]
 const dy = ysize / (Ny-1)
 # basic nodes
 # horizontal coordinates of basic grid points [m]
#  const x = SVector{Nx, Float64}([j for j = 0:dx:xsize])
 const x = [j for j = 0:dx:xsize]
 # vertical coordinates of basic grid points [m]
#  const y = SVector{Ny, Float64}([i for i = 0:dy:ysize])
 const y = [i for i = 0:dy:ysize]
 # Vx nodes
 # horizontal coordinates of vx grid points [m]
#  const xvx = SVector{Ny1, Float64}([j for j = 0:dx:xsize+dy])
 const xvx = [j for j = 0:dx:xsize+dy]
 # vertical coordinates of vx grid points [m]
#  const yvx = SVector{Nx1, Float64}([i for i = -dy/2:dy:ysize+dy/2])
 const yvx = [i for i = -dy/2:dy:ysize+dy/2]
 # Vy nodes
 # horizontal coordinates of vy grid points [m]
#  const xvy = SVector{Nx1, Float64}([j for j = -dx/2:dx:xsize+dx/2])
 const xvy = [j for j = -dx/2:dx:xsize+dx/2]
 # vertical coordinates of vy grid points [m]
#  const yvy = SVector{Ny1, Float64}([i for i = 0:dy:ysize+dy])
 const yvy = [i for i = 0:dy:ysize+dy]
 # P nodes
 # horizontal coordinates of p grid points [m]
#  const xp = SVector{Nx1, Float64}([j for j = -dx/2:dx:xsize+dx/2])
 const xp = [j for j = -dx/2:dx:xsize+dx/2]
 # vertical coordinates of p grid points [m]
#  const yp = SVector{Ny1, Float64}([i for i = -dy/2:dy:ysize+dy/2])
 const yp = [i for i = -dy/2:dy:ysize+dy/2]
 # basic grid min/max assignables indices
 # minimum assignable basic grid index in x direction
 const jmin_basic = 1
 # minimum assignable basic grid index in y direction
 const imin_basic = 1
 # maximum assignable basic grid index in x direction
 const jmax_basic = Nx - 1
 # maximum assignable basic grid index in y direction
 const imax_basic = Ny - 1
 # Vx grid min/max assignables indices
 # minimum assignable Vx grid index in x direction
 const jmin_vx = 1
 # minimum assignable Vx grid index in y direction
 const imin_vx = 1
 # maximum assignable Vx grid index in x direction
 const jmax_vx = Nx - 1
 # maximum assignable Vx grid index in y direction
 const imax_vx = Ny
 # Vy grid min/max assignables indices
 # minimum assignable Vy grid index in x direction
 const jmin_vy = 1
 # minimum assignable Vy grid index in y direction
 const imin_vy = 1
 # maximum assignable Vy grid index in x direction
 const jmax_vy = Nx
 # maximum assignable Vy grid index in y direction
 const imax_vy = Ny - 1
 # P grid min/max assignables indices
 # minimum assignable P grid index in x direction
 const jmin_p = 1
 # minimum assignable P grid index in y direction
 const imin_p = 1
 # maximum assignable P grid index in x direction
 const jmax_p = Nx
 # maximum assignable P grid index in y direction
 const imax_p = Ny
 # planetary parameters
 # planetary radius [m]
 const rplanet = 5_000.0
 # crust radius [m]
 const rcrust = 4_800.0
 # surface pressure [Pa]
 const psurface = 1e+3
 # marker count and initial spacing
 # number of markers per cell in horizontal direction
 const Nxmc = 4
 # number of markers per cell in vertical direction
 const Nymc = 4
 # marker grid resolution in horizontal direction
 const Nxm = (Nx - 1) * Nxmc
 # marker grid resolution in vertical direction
 const Nym = (Ny - 1) * Nymc
 # marker grid step in horizontal direction
 const dxm = xsize / Nxm
 # marker grid step in vertical direction
 const dym = ysize / Nym
 # horizontal coordinates of marker grid/launch anchor points [m]
#  const xxm = SVector{Nxm, Float64}([j for j = dxm/2:dxm:xsize-dxm/2])
 const xxm = [j for j = dxm/2:dxm:xsize-dxm/2]
 # vertical coordinates of marker grid/launch anchor points [m]
 const yym = [i for i = dym/2:dym:ysize-dym/2]
 # initialization distance of nearest marker to launch anchor point [m]
 const mdis_init = 1.0e30
 # number of markers at start
 const start_marknum = Nxm * Nym
 # marker grid min/max assignables indices
 # minimum assignable marker grid index in x direction
 const jmin_m = 1
 # minimum assignable marker grid index in y direction
 const imin_m = 1
 # maximum assignable marker grid index in x direction
 const jmax_m = Nxm - 1
 # maximum assignable marker grid index in y direction
 const imax_m = Nym - 1
 # physical constants
 # gravitational constant [m^3*kg^-1*s^-2]
 const G = 6.672e-11
 # scaled pressure    
 const pscale = 1.0e+23 / dx
 # materials properties:              planet      crust       space
 # solid Density [kg/m^3]
 const rhosolidm = SVector{3, Float64}(    [3300.0    , 3300.0    ,    1.0    ])
 # fluid density [kg/m^3]	
 const rhofluidm = SVector{3, Float64}(    [7000.0    , 7000.0    ,    1.0    ])
 # solid viscosity [Pa*s]
 const etasolidm = SVector{3, Float64}(    [   1.0e+16,    1.0e+16,    1.0e+14])
 # molten solid viscosity [Pa*s]
 const etasolidmm = SVector{3, Float64}(   [   1.0e+14,    1.0e+14,    1.0e+14])
 # fluid viscosity [Pa*s]
 const etafluidm = SVector{3, Float64}(    [   1.0e-02,    1.0e-02,    1.0e+12])
 # molten fluid viscosity [Pa*s]
 const etafluidmm = SVector{3, Float64}(   [   1.0e-02,    1.0e-02,    1.0e+12])
 # solid volumetric heat capacity [kg/m^3]
 const rhocpsolidm = SVector{3, Float64}(  [   3.3e+06,    3.3e+06,    3.0e+06])
 # fluid volumetric heat capacity [kg/m^3]
 const rhocpfluidm = SVector{3, Float64}(  [   7.0e+06,    7.0e+06,    3.0e+06])
 # solid thermal expansion [1/K]
 const alphasolidm = SVector{3, Float64}(  [   3.0e-05,    3.0e-05,    0.0    ])
 # fluid thermal expansion [1/K]
 const alphafluidm = SVector{3, Float64}(  [   5.0e-05,    5.0e-05,    0.0    ])
 # solid thermal conductivity [W/m/K]
 const ksolidm = SVector{3, Float64}(      [   3.0    ,    3.0    , 3000.0    ])
 # fluid thermal conductivity [W/m/K]
 const kfluidm = SVector{3, Float64}(      [  50.0    ,   50.0    , 3000.0    ])
 # solid radiogenic heat production [W/m^3]
 const start_hrsolidm = SVector{3,Float64}([   0.0    ,    0.0    ,    0.0    ])
 # fluid radiogenic heat production [W/m^3]
 const start_hrfluidm = SVector{3,Float64}([   0.0    ,    0.0    ,    0.0    ])
 # solid shear modulus [Pa]
 const gggsolidm = SVector{3, Float64}(    [   1.0e+10,    1.0e+10,    1.0e+10])
 # solid friction coefficient
 const frictsolidm = SVector{3, Float64}(  [   0.6    ,    0.6    ,    0.0    ])
 # solid compressive strength [Pa]
 const cohessolidm = SVector{3, Float64}(  [   1.0e+08,    1.0e+08,    1.0e+08])
 # solid tensile strength [Pa]
 const tenssolidm  = SVector{3, Float64}(  [   6.0e+07,    6.0e+07,    6.0e+07])
 # standard permeability [m^2]
 const kphim0 = SVector{3, Float64}(       [   1.0e-13,    1.0e-13,    1.0e-17])
 # initial temperature [K]
 const tkm0 = SVector{3, Float64}(         [ 300.0    ,  300.0    ,  273.0    ])
 # Coefficient to compute compaction viscosity from shear viscosity
 const etaphikoef = 1e-4
 # 26Al decay
 # 26Al half life [s]
 const t_half_al = 717000 * 31540000
 # 26Al decay constant
 const tau_al = t_half_al / log(2)
 # initial ratio of 26Al and 27Al isotopes
 const ratio_al = 5.0e-5
 # E 26Al [J]
 const E_al = 5.0470e-13
 # 26Al atoms/kg
 const f_al = 1.9e23
 # 60Fe decay
 # 60Fe half life [s]	
 const t_half_fe = 2620000 * 31540000
 # 60Fe decay constant
 const tau_fe = t_half_fe / log(2)
 # initial ratio of 60Fe and 56Fe isotopes	
 const ratio_fe = 1e-6
 # E 60Fe [J]	
 const E_fe = 4.34e-13
 # 60Fe atoms/kg	
 const f_fe = 1.957e24
 # melting temperatures
 # silicate melting temperature [K]
 const tmsilicate = 1e+6
 # iron melting temperature [K]
 const tmiron = 1273 
 # porosities
 # standard Fe fraction [porosity]
 const phim0 = 0.2
 # min porosity	
 const phimin = 1e-4
 # max porosity
 const phimax = 1 - phimin            
 # mechanical boundary conditions: free slip=-1 / no slip=1
 # mechanical boundary condition left
 const bcleft = -1
 # mechanical boundary condition right
 const bcright = -1
 # mechanical boundary condition top
 const bctop = -1
 # mechanical boundary condition bottom
 const bcbottom = -1
 # hydraulic boundary conditions: free slip=-1 / no slip=1
 # hydraulic boundary condition left
 const bcfleft = -1
 # hydraulic boundary condition right
 const bcfright = -1
 # hydraulic boundary condition top
 const bcftop = -1
 # hydraulic boundary condition bottom
 const bcfbottom = -1
 # extension/shortening velocities
 # shortening strain rate
 const strainrate = 0e-13
 # x extension/shortening velocity left
 const vxleft = strainrate * xsize / 2
 # x extension/shortening velocity right
 const vxright= -strainrate * xsize / 2
 # y extension/shortening velocity top
 const vytop = - strainrate * ysize / 2
 # y extension/shortening velocity bottom
 const vybottom = strainrate * ysize / 2
 # Runge-Kutta integration parameters
 # bⱼ Butcher coefficients for RK4
 const brk4 = SVector{4, Float64}([1/6, 2/6, 2/6, 1/6])
 # cⱼ Butcher coefficients for RK4
 const crk4 = SVector{3, Float64}([0.5, 0.5, 1.0])
 # timestepping parameters
 # output storage periodicity
 const savematstep = 10
 # Maximal computational timestep [s]
 const dtelastic = 1e+11 
 # Coefficient to decrease computational timestep
 const dtkoef = 2 
 # Coefficient to increase computational timestep
 const dtkoefup = 1.1 
 # Number of iterations before changing computational timestep
 const dtstep = 200 
 # Max marker movement per time step [grid steps]
 const dxymax = 0.05 
 # Weight of averaged velocity for moving markers
 const vpratio = 1 / 3 
 # Max temperature change per time step [K]
 const DTmax = 20 
 # Subgrid temperature diffusion parameter
 const dsubgridt = 0 
 # Subgrid stress diffusion parameter
 const dsubgrids = 0
 # length of year [s]
 const yearlength = 365.25 * 24 * 3600
 # Time sum (start) [s]
 const start_time = 1e6 * yearlength 
 # Time sum (end) [s]
 const endtime = 15 * 1000000 * yearlength
 # Lower viscosity cut-off [Pa s]	
 const etamin = 1e+12 
 # Upper viscosity cut-off [Pa s]
 const etamax = 1e+23 
 # Number of plastic iterations
 const nplast = 100000
 # Periodicity of visualization
 const visstep = 1 
 # Tolerance level for yielding error()
 const yerrmax = 1e+2 
 # Weight for old viscosity
 const etawt = 0 
 # max porosity ratio change per time step
 const dphimax = 0.01
 # starting timestep
 const start_step = 1
 # number of timesteps to run
 const nsteps = 20
#  const nsteps = 30000 
# radioactive switches
# radioactive heating from 26Al active
const hr_al = true
# radioactive heating from 60Fe active	
const hr_fe = true
# planetary parameters
# planetary radius [m]
const rplanet = 50_000.0
# crust radius [m]
const rcrust = 48_000.0
# surface pressure [Pa]
const psurface = 1.0e+3
# model size, geometry, and resolution
# horizontal model size [m]
const xsize = 140_000.0
# vertical model size [m]
const ysize = 140_000.0
# horizontal center of model
const xcenter = xsize / 2
# vertical center of model
const ycenter = ysize / 2  
# basic grid resolution in x direction (horizontal)
const Nx = 141
# basic grid resolution in y direction (vertical)	
const Ny = 141
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
const x = SVector{Nx, Float64}([j for j = 0:dx:xsize])
#  const x = [j for j = 0:dx:xsize]
# vertical coordinates of basic grid points [m]
const y = SVector{Ny, Float64}([i for i = 0:dy:ysize])
#  const y = [i for i = 0:dy:ysize]
# Vx nodes
# horizontal coordinates of vx grid points [m]
const xvx = SVector{Ny1, Float64}([j for j = 0:dx:xsize+dy])
#  const xvx = [j for j = 0:dx:xsize+dy]
# vertical coordinates of vx grid points [m]
const yvx = SVector{Nx1, Float64}([i for i = -dy/2:dy:ysize+dy/2])
#  const yvx = [i for i = -dy/2:dy:ysize+dy/2]
# Vy nodes
# horizontal coordinates of vy grid points [m]
const xvy = SVector{Nx1, Float64}([j for j = -dx/2:dx:xsize+dx/2])
#  const xvy = [j for j = -dx/2:dx:xsize+dx/2]
# vertical coordinates of vy grid points [m]
const yvy = SVector{Ny1, Float64}([i for i = 0:dy:ysize+dy])
#  const yvy = [i for i = 0:dy:ysize+dy]
# P nodes
# horizontal coordinates of p grid points [m]
const xp = SVector{Nx1, Float64}([j for j = -dx/2:dx:xsize+dx/2])
#  const xp = [j for j = -dx/2:dx:xsize+dx/2]
# vertical coordinates of p grid points [m]
const yp = SVector{Ny1, Float64}([i for i = -dy/2:dy:ysize+dy/2])
#  const yp = [i for i = -dy/2:dy:ysize+dy/2]
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
const xxm = SVector{Nxm, Float64}([j for j = dxm/2:dxm:xsize-dxm/2])
#  const xxm = [j for j = dxm/2:dxm:xsize-dxm/2]
# vertical coordinates of marker grid/launch anchor points [m]
const yym = SVector{Nxm, Float64}([i for i = dym/2:dym:ysize-dym/2])
#  const yym = [i for i = dym/2:dym:ysize-dym/2]
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
# marker randomized positions and porosity for testing
const random_markers = false
# physical constants
# gravitational constant [m^3*kg^-1*s^-2]
const G = 6.672e-11
# scaled pressure    
# const pscale = 1.0e+23 * inv(dx)
# pressure scaling coefficient (eqn 7.19-7.21 in Gerya(2019))
# const Kcont = 2.0 * 1.0e15 * inv(dx+dy)
const Kcont = 1.0e20
# ------------------------------------------------------------------------------
# materials properties:                       planet      crust       space
# all identical except fluid density -------------------------------------------
#  # solid Density [kg/m^3]
#  const rhosolidm = SVector{3, Float64}(    [3300.0    , 3300.0    ,    1.0    ])
#  # fluid density [kg/m^3]	
#  const rhofluidm = SVector{3, Float64}(    [7000.0    , 7000.0    ,    1.0    ])
#  # solid viscosity [Pa*s]
#  const etasolidm = SVector{3, Float64}(    [   1.0e+16,    1.0e+16,    1.0e+16])
#  # molten solid viscosity [Pa*s]
#  const etasolidmm = SVector{3, Float64}(   [   1.0e+14,    1.0e+14,    1.0e+14])
#  # fluid viscosity [Pa*s]
#  const etafluidm = SVector{3, Float64}(    [   1.0e-02,    1.0e-02,    1.0e-02])
#  # molten fluid viscosity [Pa*s]
#  const etafluidmm = SVector{3, Float64}(   [   1.0e-02,    1.0e-02,    1.0e-02])
#  # solid volumetric heat capacity [kg/m^3]
#  const rhocpsolidm = SVector{3, Float64}(  [   3.3e+06,    3.3e+06,    3.3e+06])
#  # fluid volumetric heat capacity [kg/m^3]
#  const rhocpfluidm = SVector{3, Float64}(  [   7.0e+06,    7.0e+06,    7.0e+06])
#  # solid thermal expansion [1/K]
#  const alphasolidm = SVector{3, Float64}(  [   3.0e-05,    3.0e-05,    3.0e-05])
#  # fluid thermal expansion [1/K]
#  const alphafluidm = SVector{3, Float64}(  [   5.0e-05,    5.0e-05,    5.0e-05])
#  # solid thermal conductivity [W/m/K]
#  const ksolidm = SVector{3, Float64}(      [   3.0    ,    3.0    ,    3.0    ])
#  # fluid thermal conductivity [W/m/K]
#  const kfluidm = SVector{3, Float64}(      [  50.0    ,   50.0    ,   50.0    ])
#  # solid radiogenic heat production [W/m^3]
#  const start_hrsolidm = SVector{3,Float64}([   0.0    ,    0.0    ,    0.0    ])
#  # fluid radiogenic heat production [W/m^3]
#  const start_hrfluidm = SVector{3,Float64}([   0.0    ,    0.0    ,    0.0    ])
#  # solid shear modulus [Pa]
#  const gggsolidm = SVector{3, Float64}(    [   1.0e+10,    1.0e+10,    1.0e+10])
#  # solid friction coefficient
#  const frictsolidm = SVector{3, Float64}(  [   0.6    ,    0.6    ,    0.0    ])
#  # solid compressive strength [Pa]
#  const cohessolidm = SVector{3, Float64}(  [   1.0e+08,    1.0e+08,    1.0e+08])
#  # solid tensile strength [Pa]
#  const tenssolidm  = SVector{3, Float64}(  [   6.0e+07,    6.0e+07,    6.0e+07])
#  # standard permeability [m^2]
#  const kphim0 = SVector{3, Float64}(       [   1.0e-13,    1.0e-13,    1.0e-13])
#  # initial temperature [K]
#  const tkm0 = SVector{3, Float64}(         [ 300.0    ,  300.0    ,  300.0    ])
# fluid Fe ---------------------------------------------------------------------
# # solid Density [kg/m^3]
# const rhosolidm = SVector{3, Float64}(    [3300.0    , 3300.0    ,    1.0    ])
# # fluid density [kg/m^3]	
# const rhofluidm = SVector{3, Float64}(    [7000.0    , 7000.0    ,    1.0    ])
# # solid viscosity [Pa*s]
# const etasolidm = SVector{3, Float64}(    [   1.0e+16,    1.0e+16,    1.0e+14])
# # molten solid viscosity [Pa*s]
# const etasolidmm = SVector{3, Float64}(   [   1.0e+14,    1.0e+14,    1.0e+14])
# # fluid viscosity [Pa*s]
# const etafluidm = SVector{3, Float64}(    [   1.0e-02,    1.0e-02,    1.0e+12])
# # molten fluid viscosity [Pa*s]
# const etafluidmm = SVector{3, Float64}(   [   1.0e-02,    1.0e-02,    1.0e+12])
# # solid volumetric heat capacity [kg/m^3]
# const rhocpsolidm = SVector{3, Float64}(  [   3.3e+06,    3.3e+06,    3.0e+06])
# # fluid volumetric heat capacity [kg/m^3]
# const rhocpfluidm = SVector{3, Float64}(  [   7.0e+06,    7.0e+06,    3.0e+06])
# # solid thermal expansion [1/K]
# const alphasolidm = SVector{3, Float64}(  [   3.0e-05,    3.0e-05,    0.0    ])
# # fluid thermal expansion [1/K]
# const alphafluidm = SVector{3, Float64}(  [   5.0e-05,    5.0e-05,    0.0    ])
# # solid thermal conductivity [W/m/K]
# const ksolidm = SVector{3, Float64}(      [   3.0    ,    3.0    , 3000.0    ])
# # fluid thermal conductivity [W/m/K]
# const kfluidm = SVector{3, Float64}(      [  50.0    ,   50.0    , 3000.0    ])
# # solid radiogenic heat production [W/m^3]
# const start_hrsolidm = SVector{3,Float64}([   0.0    ,    0.0    ,    0.0    ])
# # fluid radiogenic heat production [W/m^3]
# const start_hrfluidm = SVector{3,Float64}([   0.0    ,    0.0    ,    0.0    ])
# # solid shear modulus [Pa]
# const gggsolidm = SVector{3, Float64}(    [   1.0e+10,    1.0e+10,    1.0e+10])
# # solid friction coefficient
# const frictsolidm = SVector{3, Float64}(  [   0.6    ,    0.6    ,    0.0    ])
# # solid compressive strength [Pa]
# const cohessolidm = SVector{3, Float64}(  [   1.0e+08,    1.0e+08,    1.0e+08])
# # solid tensile strength [Pa]
# const tenssolidm  = SVector{3, Float64}(  [   6.0e+07,    6.0e+07,    6.0e+07])
# # standard permeability [m^2]
# const kphim0 = SVector{3, Float64}(       [   1.0e-13,    1.0e-13,    1.0e-17])
# # initial temperature [K]
# const tkm0 = SVector{3, Float64}(         [ 300.0    ,  300.0    ,  273.0    ])
# # Coefficient to compute compaction viscosity from shear viscosity
# const etaphikoef = 1e-4
# fluid H₂O -------------------------------------------------------------------
# solid Density [kg/m^3]
const rhosolidm = SVector{3, Float64}(    [3300.0    , 3300.0    ,    1.0    ])
# fluid density [kg/m^3]	
const rhofluidm = SVector{3, Float64}(    [1000.0    , 1000.0    ,    1.0    ])
# solid viscosity [Pa*s]
const etasolidm = SVector{3, Float64}(    [   1.0e+19,    1.0e+19,    1.0e+16])
# molten solid viscosity [Pa*s]
const etasolidmm = SVector{3, Float64}(   [   1.0e+19,    1.0e+19,    1.0e+16])
# fluid viscosity [Pa*s]
const etafluidm = SVector{3, Float64}(    [   1.0e-02,    1.0e-02,    1.0e+12])
# molten fluid viscosity [Pa*s]
const etafluidmm = SVector{3, Float64}(   [   1.0e-02,    1.0e-02,    1.0e+12])
# solid volumetric heat capacity [kg/m^3]
const rhocpsolidm = SVector{3, Float64}(  [   3.3e+06,    3.3e+06,    3.0e+06])
# fluid volumetric heat capacity [kg/m^3]
const rhocpfluidm = SVector{3, Float64}(  [   1.0e+06,    1.0e+06,    3.0e+06])
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
const etaphikoef = 1
# ------------------------------------------------------------------------------
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
const strainrate = 0.0e-13
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
const brk4 = SVector{4, Rational{Int64}}([1//6, 2//6, 2//6, 1//6])
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
const start_time = 1.0e6 * yearlength 
# Time sum (end) [s]
const endtime = 15.0e6 * yearlength
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
 const n_steps = 10#30_000 
# using MKL Pardiso solver
const use_pardiso = false
# MKL Pardiso solver IPARM control parameters -> ATTN: zero-indexed as in docs:
# https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/sparse-solver-routines/onemkl-pardiso-parallel-direct-sparse-solver-iface/pardiso-iparm-parameter.html
const iparms = Dict([
    (0, 1), # in: do not use default
    (1, 2), # in: nested dissection from METIS
    (2, 0), # in: reserved, set to zero
    (3, 0), # in: no CGS/CG iterations
    (4, 0), # in: no user permutation
    (5, 0), # in: write solution on x or RHS b
    (6, 0), # out: number of iterative refinement steps performed
    (7, 20), # in: maximum number of iterative refinement steps
    (8, 0), # in: tolerance level for relative residual, only with iparm[23]=1
    (9, 13), # in: pivoting perturbation
    (10, 1), # in: scaling vectors
    (11, 1), # in: solve AX=B (no transpose, conjugate transpose): CSC matrix
    (12, 1), # in: improved accuracy using (non-)symmetric weighted matching
    (13, 0), # out: number of perturbed pivots
    (14, 0), # out: peak memory on symbolic factorization
    (15, 0), # out: permanent memory on symbolic factorization
    (16, 0), # out: size of factors/peak memory on symbolic factorization
    (17, -1), # in/out: report number of non-zero elements in the factors
    (18, -1), # in/out: report the number of FLOPs to factor matrix A
    (19, 0), # out: report CG/CGS diagnostics, iterations
    (20, 0), # in: pivoting for symmetric indefinite matrices
    (21, 0), # out: inertia: number of positive eigenvalues
    (22, 0), # out: inertia: number of negative eigenvalues
    (23, 10), # in: parallel factorization control, REQ: iparm[10]==iparm[12]==0
    (24, 0), # in: parallel forward/backward solve control
    (25, 0), # reserved, set to zero
    (26, 1), # in: matrix checker: checks ia, ja sorting order
    (27, 0), # in: single or double precision
    (28, 0), # reserved, set to zero
    (29, 0), # out: number of zero or negative pivots
    (30, 0), # in: partial solve
    (31, 0), # reserved, set to zero
    (32, 0), # reserved, set to zero
    (33, 0), # in: optimal number of OpenMP threads for CNR mode
    (34, 0), # in: one- or zero-based indexing of columns and rows
    (35, 0), # in/out: Schur complement matrix computation control
    (36, 0), # in: format for matrix storage: CSR or BSR or VBSR
    (37, 0), # reserved, set to zero
    (38, 0), # in: enable low-rank update for for multiple similar matrices
    (39, 0), # reserved, set to zero
    (40, 0), # reserved, set to zero
    (41, 0), # reserved, set to zero
    (42, 0), # in: compute diagonal of inverse matrix
    (43, 0), # reserved, set to zero 
    (44, 0), # reserved, set to zero 
    (45, 0), # reserved, set to zero 
    (46, 0), # reserved, set to zero 
    (47, 0), # reserved, set to zero 
    (48, 0), # reserved, set to zero 
    (49, 0), # reserved, set to zero 
    (50, 0), # reserved, set to zero 
    (51, 0), # reserved, set to zero 
    (52, 0), # reserved, set to zero 
    (53, 0), # reserved, set to zero 
    (54, 0), # reserved, set to zero 
    (55, 0), # in: diagonal and pivoting control
    (56, 0), # reserved, set to zero 
    (57, 0), # reserved, set to zero 
    (58, 0), # reserved, set to zero 
    (59, 0), # in: in-core (IC) or out-of-core (OOC) PARDISO mode
    (60, 0), # reserved, set to zero 
    (61, 0), # reserved, set to zero 
    (62, 0), # out: size of the minimum OOC memory for factorization and sol 
    (63, 0), # reserved, set to zero 
])
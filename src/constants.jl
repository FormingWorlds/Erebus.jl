# radioactive switches
# solid phase radioactive heating from 26Al active
const hr_al = true
# fluid phase radioactive heating from 60Fe active	
const hr_fe = false
# planetary parameters
# planetary radius [m]
const rplanet = 50_000.0
# crust radius [m]
const rcrust = 50_000.0
# surface pressure [Pa]
const psurface = 1.0e+3
# initialization distance of nearest marker to launch anchor point [m]
const mdis_init = 1.0e30
# marker randomized positions and porosity for testing
const random_markers = true
# const random_markers = false
# physical constants
# gravitational constant [m³*kg⁻¹*s⁻²] (CODATA 2018)
const G = 6.67430e-11
const GRAVITATIONAL_CONSTANT = G
const G_GRAV = G

# Boltzmann constant [J/K] (CODATA 2018 exact definition)
const BOLTZMANN_CONSTANT = 1.380649e-23
const K_BOLTZMANN = BOLTZMANN_CONSTANT

# Avogadro constant [1/mol] (CODATA 2018 exact definition)
const AVOGADRO_CONSTANT = 6.02214076e23

# Unified atomic mass unit (amu, Dalton) [kg] (CODATA 2018)
const ATOMIC_MASS_UNIT = 1.66053906660e-27

# Proton mass [kg] (CODATA 2018)
const M_PROTON = 1.67262192e-27

# Universal gas constant [J/(mol*K)] (CODATA 2018 exact definition: R = N_A * k_B)
const R_GAS = 8.31446261815324
const RG = R_GAS

# Astronomical unit [m] (IAU 2012 exact definition)
const AU_METERS = 1.495978707e11

# Solar mass [kg] (IAU standard nominal solar mass)
const M_SUN_KG = 1.98847e30

# Julian year length [s] (IAU standard: 365.25 days * 86,400 s/day)
const SECONDS_PER_YEAR = 31_557_600.0
const SEC_PER_YEAR = SECONDS_PER_YEAR

# Molar masses [kg/mol] (IUPAC 2013 standard atomic weights / conventional values)
const M_H = 0.001008
const M_C = 0.012011
const M_N = 0.014007
const M_O = 0.015999
const M_S = 0.03206
const M_Fe = 0.055845

# Compound species molar masses [kg/mol] (exact stoichiometric sums of constituent elements)
const M_H2O = 2.0 * M_H + M_O
const MH₂O = M_H2O
const MH2O = M_H2O
const M_CO = M_C + M_O
const M_CO2 = M_C + 2.0 * M_O
const M_CH4 = M_C + 4.0 * M_H
const M_N2 = 2.0 * M_N
const M_NH3 = M_N + 3.0 * M_H
const M_H2S = 2.0 * M_H + M_S
const M_S2 = 2.0 * M_S
const M_SO2 = M_S + 2.0 * M_O
const M_FeS = M_Fe + M_S
# scaled pressure    
# pressure scaling coefficient (eqn 7.19-7.21 in Gerya(2019))
# const Kcont = 2.0 * 1.0e15 * inv(dx+dy)
const Kcont = 1.0e20
# planetesimals: fluid phase H₂O -----------------------------------------------
# solid density [kg/m³]
const rhosolidm = SVector{3,Float64}([3300.0, 3300.0, 1.0])
# fluid density [kg/m³]	
const rhofluidm = SVector{3,Float64}([1000.0, 1000.0, 1.0])
# solid viscosity [Pa*s]
const etasolidm = SVector{3,Float64}([1.0e+19, 1.0e+19, 1.0e+16])
# molten solid viscosity [Pa*s]
const etasolidmm = SVector{3,Float64}([1.0e+19, 1.0e+19, 1.0e+16])
# fluid viscosity [Pa*s]
const etafluidm = SVector{3,Float64}([1.0e+12, 1.0e+12, 1.0e-03])
# molten fluid viscosity [Pa*s]
const etafluidmm = SVector{3,Float64}([1.0e-03, 1.0e-03, 1.0e-03])
# solid volumetric heat capacity [kg/m³]
const rhocpsolidm = SVector{3,Float64}([3.3e+06, 3.3e+06, 3.0e+06])
# fluid volumetric heat capacity [kg/m³]
const rhocpfluidm = SVector{3,Float64}([1.0e+06, 1.0e+06, 3.0e+06])
# solid thermal expansion [1/K]
const alphasolidm = SVector{3,Float64}([3.0e-05, 3.0e-05, 0.0])
# fluid thermal expansion [1/K]
const alphafluidm = SVector{3,Float64}([5.0e-05, 5.0e-05, 0.0])
# solid thermal conductivity [W/m/K]
const ksolidm = SVector{3,Float64}([3.0, 3.0, 3000.0])
# fluid thermal conductivity [W/m/K]
const kfluidm = SVector{3,Float64}([50.0, 50.0, 3000.0])
# solid radiogenic heat production [W/m³]
const start_hrsolidm = SVector{3,Float64}([0.0, 0.0, 0.0])
# fluid radiogenic heat production [W/m³]
const start_hrfluidm = SVector{3,Float64}([0.0, 0.0, 0.0])
# solid shear modulus [Pa]
const gggsolidm = SVector{3,Float64}([1.0e+10, 1.0e+10, 1.0e+10])
# solid friction coefficient
const frictsolidm = SVector{3,Float64}([0.6, 0.6, 0.0])
# solid compressive strength [Pa]
const cohessolidm = SVector{3,Float64}([1.0e+08, 1.0e+08, 1.0e+08])
# solid tensile strength [Pa]
const tenssolidm = SVector{3,Float64}([6.0e+07, 6.0e+07, 6.0e+07])
# solid matrix compressibility [1/Pa]
const betasolid = 2.5e-11
# fluid compressibility [1/Pa]
const betafluid = 4.0e-10
# standard permeability [m^2]
const kphim0 = SVector{3,Float64}([1.0e-13, 1.0e-13, 1.0e-17])
# initial temperature [K]
const tkm0 = SVector{3,Float64}([170.0, 170.0, 170.0])
# initial wet solid molar fraction
# marker property mode (1: dynamic calculations, 9: static parameters)
const marker_property_mode = 9
# coefficient to compute compaction viscosity from shear viscosity
const etaphikoef = 1
# melt-weakening coefficient (16.67)
const αη = 28.0
# ------------------------------------------------------------------------------
# 26Al decay
# 26Al half life [s]
const t_half_al = 717_000 * 31_540_000
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
const t_half_fe = 2_620_000 * 31_540_000
# 60Fe decay constant
const tau_fe = t_half_fe / log(2)
# initial ratio of 60Fe and 56Fe isotopes	
const ratio_fe = 1.0e-6
# E 60Fe [J]	
const E_fe = 4.34e-13
# 60Fe atoms/kg	
const f_fe = 1.957e24
# melting
# solid phase (silicate) melting temperature [K]
const tmsolidphase = 1416.0#1.0e+6
# fluid phase (H₂O) melting temperature [K]
const tmfluidphase = 273.0
# fluid phase (H₂O) heat of fusion of ice [J/kg]
const Lᶠ = 333.55e3
# porosities
# standard H₂O fraction [porosity]
const phim0 = 0.2
# min porosity	
const phimin = 1.0e-4
# max porosity
const phimax = 1.0 - phimin
# thermodynamic parameters: silicate dehydration reaction Wˢ = Dˢ + H₂O
# molar mass of dry silicate [kg/mol]
const MD = 0.120
# density of dry silicate [kg/m³]
const ρDˢ = 3300.0
# density of wet silicate [kg/m³]
const ρWˢ = 2600.0
# density of fluid (liquid H₂O) [kg/m³]
const ρH₂Oᶠ = 1000.0
# density of fluid ice (frozen H₂O) [kg/m³]
const ρH₂Oᶠⁱ = 917.0
# molar volume of dry solid [m³/mol]
const VDˢ = MD / ρDˢ
# molar volume of wet solid [m³/mol]
const VWˢ = (MD+MH₂O) / ρWˢ
# molar volume of fluid (liquid H₂O) [m³/mol]
const VH₂Oᶠ = MH₂O / ρH₂Oᶠ
# molar volume of fluid ice (frozen H₂O) [m³/mol]
const VH₂Oᶠⁱ = MH₂O / ρH₂Oᶠⁱ
# enthalpy change for dehydration of the wet silicate [J/mol]
const ΔHWD = 40000.0
# entropy change for dehydration of the wet silicate [J/K/mol]
const ΔSWD = 60.0
# volume change for dehydration of the wet silicate [m³/mol]
const ΔVWD = VDˢ + VH₂Oᶠ - VWˢ
# coefficient of pressure from previous hydrothermomechanical iteration
const pfcoeff = 0.5
# error limit to exit thermochemical iterations
const pferrmax = 1.0e+5
# reaction activation switch
const reaction_active = true
# time to run dehydration reaction to completion [s]
const Δtreaction = 1.0e+10
# log reaction completion rate ln(ρend/ρstart)
const log_completion_rate = log(0.01)
# reaction constant mode (1: [Martin & Fyfe, 1970; Emmanuel & Berkowitz, 2006;
# Iyer et al., 2012], 2: [Bland & Travis, 2017], 9: constant Δtreaction)
const reaction_rate_coeff_mode = 1
# reaction constant parameters mode 1 [Iyer et al., 2012]
# (A: kinetic coefficient, b: kinetic coefficient, c: kinetic coefficient)
A_I = 1.0e-11;
b_I = 2.5e-4;
c_I = 543.0
# reaction constant parameters mode 2 [Bland & Travis, 2017]
# (Sxo_B: reaction rate at ref T, Tscl_B: empirical scaling factor, To_B: reaction ref T)
Sxo_B = 2.0e-11;
Tscl_B = 10.0;
To_B = 293.0
# reaction constant parameters mode 3 [Travis et al., 2018]
# (Sxo_T: reaction rate at ref T, To_T: reaction ref T, Ea_T: reaction activation energy)
Sxo_T = 2.0e-11;
To_T = 293.0;
Ea_T = 63.8e3
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
# shortening strain rate
const strainrate = 0.0e-13
# Runge-Kutta integration parameters
# bⱼ Butcher coefficients for RK4
const brk4 = SVector{4,Rational{Int64}}([1//6, 2//6, 2//6, 1//6])
# cⱼ Butcher coefficients for RK4
const crk4 = SVector{3,Float64}([0.5, 0.5, 1.0])
# timestepping parameters
# output storage periodicity
# const savematstep = 10
# longest allowed computational timestep [s]
# const dt_longest = 1.0e+11
# coefficient to decrease computational timestep
# const dtcoefdn = 0.5
# coefficient to increase computational timestep
# const dtcoefup = 1.2
# number of iterations before changing computational timestep
# const dtstep = 200
# max marker movement per time step [grid steps]
# const dxymax = 0.05
# weight of averaged velocity for moving markers
# const vpratio = 1 / 3
# max temperature change per time step [K]
# const DTmax = 20.0
# subgrid temperature diffusion parameter
# const dsubgridt = 0.0
# subgrid stress diffusion parameter
# const dsubgrids = 0.0
# length of year [s]
# const yearlength = 365.25 * 24 * 3600
# time sum (start) [s]
# const start_time = 2.25e6 * yearlength
# time sum (end) [s]
# const endtime = 15.0e6 * yearlength
# lower viscosity cut-off [Pa s]	
# const etamin = 1e+12
# upper viscosity cut-off [Pa s]
# const etamax = 1e+23
# maximum number of plastic iterations
# const nplast = 100_000
# maximum number of global iterations
# const titermax = 10_000
# periodicity of visualization
# const visstep = 1
# tolerance level for yielding error()
# const yerrmax = 1e+2
# weight for old viscosity
# const etawt = 0.0
# max porosity ratio change per time step
# const dphimax = 100.01
# const dphimax = 0.01
# starting timestep
# const start_step = 1
# maximum number of timesteps to run
# const n_steps = 10
# # const n_steps = 100 
# const n_steps = 30_000 
# random number generator seed
# const seed = 42
# using MKL Pardiso solver
# const use_pardiso = false
# MKL Pardiso solver IPARM control parameters -> ∇ATTN: zero-indexed as in docs:
# https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/sparse-solver-routines/onemkl-pardiso-parallel-direct-sparse-solver-iface/pardiso-iparm-parameter.html
const iparms_dict = Dict([
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
# MKL Pardiso solver IPARM control parameters formatted for LinearSolve.jl
const iparms = collect(
    (key + 1, iparms_dict[key]) for key in sort!(collect(keys(iparms_dict)))
)
# LinearSolve.jl solver keyword arguments
const cache_kwargs = (;
    nprocs=4, verbose=true, abstol=1e-8, reltol=1e-8, maxiter=30, iparm=iparms
)

# Volatile species set in thermodynamic equilibrium surface speciation
const SPECIATION_SPECIES = Set{Symbol}([
    :H2, :H2O, :CO, :CO2, :CH4, :N2, :NH3, :H2S, :S2, :SO2
])

"""
Standard atomic and molecular weights for planetary volatile and escape species [amu].
Values follow IUPAC standard atomic weights (Meija et al. 2016).
"""
const SPECIES_AMU = Dict{Symbol,Float64}(
    :H => M_H * 1000.0,
    :D => 2.0141,
    :He => 4.0026,
    :C => M_C * 1000.0,
    :N => M_N * 1000.0,
    :O => M_O * 1000.0,
    :Ne => 20.1797,
    :Na => 22.990,
    :Mg => 24.305,
    :Si => 28.085,
    :S => M_S * 1000.0,
    :Ar => 39.948,
    :Fe => M_Fe * 1000.0,
    :Kr => 83.798,
    :Xe => 131.293,
    :H2 => 2.0 * M_H * 1000.0,
    :H2O => M_H2O * 1000.0,
    :CO => M_CO * 1000.0,
    :CO2 => M_CO2 * 1000.0,
    :CH4 => M_CH4 * 1000.0,
    :N2 => M_N2 * 1000.0,
    :NH3 => M_NH3 * 1000.0,
    :O2 => 2.0 * M_O * 1000.0,
    :H2S => M_H2S * 1000.0,
    :SO2 => M_SO2 * 1000.0,
    :S2 => M_S2 * 1000.0,
)

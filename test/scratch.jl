using Base.Threads
using SparseArrays
using MAT
using DocStringExtensions
using Parameters
using StaticArrays
using BenchmarkTools
using TimerOutputs



using HydrologyPlanetesimals
run_simulation()


xsize=140000.0
ysize=140000.0
rplanet=50000.0
rcrust=48000.0
Nx=141
Ny=141
Nxmc=4
Nymc=4

sp = HydrologyPlanetesimals.StaticParameters(
    xsize=xsize,
    ysize=ysize,
    rplanet=rplanet,
    rcrust=rcrust,
    Nx=Nx,
    Ny=Ny,
    Nxmc=Nxmc,
    Nymc=Nymc
    )

 # basic nodes
# grid geometry
# x: horizontal coordinates of basic grid points [m]
x = SVector{sp.Nx, Float64}([j for j = 0:sp.dx:sp.xsize])
# y: vertical coordinates of basic grid points [m]
y = SVector{sp.Ny, Float64}([j for j = 0:sp.dy:sp.ysize])
# physical node properties
# viscoplastic viscosity, Pa*s
ETA = zeros(Float64, sp.Ny, sp.Nx)
# viscous viscosity, Pa*s
ETA0 = zeros(Float64, sp.Ny, sp.Nx)
# shear modulus, Pa
GGG = zeros(Float64, sp.Ny, sp.Nx)
# epsilonxy, 1/s
EXY = zeros(Float64, sp.Ny, sp.Nx)
# sigma0xy, 1/s
SXY0 = zeros(Float64, sp.Ny, sp.Nx)
# rotation rate, 1/s
WYX = zeros(Float64, sp.Ny, sp.Nx)
# compressive strength, Pa
COH = zeros(Float64, sp.Ny, sp.Nx)
# tensile strength, Pa
TEN = zeros(Float64, sp.Ny, sp.Nx)
# friction
FRI = zeros(Float64, sp.Ny, sp.Nx)
# plastic yielding mark, 1=yes,0=no
YNY = zeros(Int8, sp.Ny, sp.Nx)



const ETA0SUM = zeros(sp.Ny, sp.Nx, nthreads())
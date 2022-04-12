# using Base.Threads
# using SparseArrays
# using MAT
# using DocStringExtensions
# using Parameters
# using StaticArrays
# using BenchmarkTools
# using TimerOutputs



using HydrologyPlanetesimals
run_simulation()


# xsize=140000.0
# ysize=140000.0
# rplanet=50000.0
# rcrust=48000.0
# Nx=141
# Ny=141
# Nxmc=4
# Nymc=4

# sp = HydrologyPlanetesimals.StaticParameters(
#     xsize=xsize,
#     ysize=ysize,
#     rplanet=rplanet,
#     rcrust=rcrust,
#     Nx=Nx,
#     Ny=Ny,
#     Nxmc=Nxmc,
#     Nymc=Nymc
#     )

#  # basic nodes
# # grid geometry
# # x: horizontal coordinates of basic grid points [m]
# x = SVector{sp.Nx, Float64}([j for j = 0:sp.dx:sp.xsize])
# # y: vertical coordinates of basic grid points [m]
# y = SVector{sp.Ny, Float64}([j for j = 0:sp.dy:sp.ysize])
# # physical node properties
# # viscoplastic viscosity, Pa*s
# ETA = zeros(Float64, sp.Ny, sp.Nx)
# # viscous viscosity, Pa*s
# ETA0 = zeros(Float64, sp.Ny, sp.Nx)
# # shear modulus, Pa
# GGG = zeros(Float64, sp.Ny, sp.Nx)
# # epsilonxy, 1/s
# EXY = zeros(Float64, sp.Ny, sp.Nx)
# # sigma0xy, 1/s
# SXY0 = zeros(Float64, sp.Ny, sp.Nx)
# # rotation rate, 1/s
# WYX = zeros(Float64, sp.Ny, sp.Nx)
# # compressive strength, Pa
# COH = zeros(Float64, sp.Ny, sp.Nx)
# # tensile strength, Pa
# TEN = zeros(Float64, sp.Ny, sp.Nx)
# # friction
# FRI = zeros(Float64, sp.Ny, sp.Nx)
# # plastic yielding mark, 1=yes,0=no
# YNY = zeros(Int8, sp.Ny, sp.Nx)



# const ETA0SUM = zeros(sp.Ny, sp.Nx, nthreads())

# compiled marker properties used in the simulation

# tm # static
# phim # dynamic
# etavpm # static / optionally dynamic
# tkm # dynamic 
# 1/gggtotalm # static
# fricttotalm # static
# cohestotalm # static
# tenstotalm # static
# rhofluidcur # static
# alphasolidm # static
# alphafluidm # static

# rhototalm # dynamic (rock) / static (air)
# rhocptotalm # dynamic (rock) / static (air)
# tkm*rhocptotalm # dynamic
# etatotalm # dynamic / static (air)
# hrtotalm # dynamic (rock) / static (air)
# ktotalm # dynamic (rock) / static (air)
# etafluidcur/kphim # dynamic

# sxxm # dynamic
# sxym # dynamic




function apply_insulating_boundary_conditions!(t)
    # @timeit to "apply_insulating_boundary_conditions!" begin
        Ny, Nx = size(t)
        # upper boundary
        t[1, 2:Nx-1] .= t[2, 2:Nx-1]
        # lower boundary
        t[Ny, 2:Nx-1] .= t[Ny-1, 2:Nx-1]
        # left boundary
        t[:, 1] .= t[:, 2]
        # right boundary
        t[:, Nx] .= t[:, Nx-1]
    # end # @timeit to "apply_insulating_boundary_conditions!"
        return nothing
    end
    
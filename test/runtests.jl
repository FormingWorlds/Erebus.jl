using Erebus
using ExtendableSparse
using JLD2
using LinearSolve
using Random
using StaticArrays
using Test

include("../src/test_constants.jl")
# include("../src/constants.jl")
const rgen = MersenneTwister(seed)

@testset verbose=true "Erebus.jl" begin
    include("test_geometry.jl")
    include("test_physics.jl")
    include("test_particles.jl")
    include("test_numerics.jl")
    include("test_simulation.jl")
    include("test_integration.jl")
end

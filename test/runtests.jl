using Erebus
using ExtendableSparse
using JLD2
using LinearSolve
using Random
using StaticArrays
using Test

include("../src/constants.jl")
include("test_helpers.jl")
const rgen = MersenneTwister(42)

test_group = get(ENV, "EREBUS_TEST_GROUP", "all")

unit_tests = [
    "test_config.jl",
    "test_geometry.jl",
    "test_dynamic_grid.jl",
    "test_threading.jl",
    "test_solver_optimization.jl",
    "test_darcy_elimination.jl",
    "test_iterative_solver.jl",
    "test_multigrid.jl",
    "test_gpu_acceleration.jl",
    "test_distributed_mpi.jl",
    "test_physics.jl",
    "test_particles.jl",
    "test_numerics.jl",
    "test_simulation.jl",
    "test_workspace.jl",
    "test_numerical_accelerations.jl",
    "test_geometry_radiation.jl",
    "test_reaction_pathways.jl",
    "test_stefan_benchmark.jl",
    "test_thermal_slab.jl",
    "test_melting.jl",
    "test_soft_turbulence.jl",
    "test_magma_transport.jl",
    "test_sill_cooling.jl",
    "test_venting_thermodynamics.jl",
    "test_venting_darcy_sink.jl",
    "test_hydrofracture_venting.jl",
    "test_volatile_solubility.jl",
    "test_volatile_solubility_hcns.jl",
    "test_volatile_retention.jl",
    "test_redox.jl",
    "test_meteorite_suite.jl",
    "test_jeans_escape.jl",
    "test_core_formation.jl",
    "test_core_volatile_partitioning.jl",
    "test_normative_accessory_minerals.jl",
    "test_hydrothermal_convection.jl",
    "test_accretion.jl",
    "test_multistage_accretion.jl",
    "test_volatile_mixtures.jl",
    "test_atmosphere.jl",
    "test_dehydration_darcy_coupling.jl",
    "test_magma_degassing.jl",
    "test_xuv_escape.jl",
    "test_telescoping.jl",
    "test_p2m_tiled.jl",
    "test_telemetry.jl",
    "test_tooling.jl",
]

integration_tests = [
    "test_venting_integration.jl",
    "test_ensemble.jl",
    "test_tutorial_lunar_growth.jl",
    "test_integration.jl",
    "test_conservation_closure.jl",
    "test_switch_sensitivity.jl",
]

all_test_files = filter(f -> endswith(f, ".jl") && f != "runtests.jl", readdir(@__DIR__))
listed_tests = Set(vcat(unit_tests, integration_tests))
for f in all_test_files
    if f ∉ listed_tests &&
        f != "test_helpers.jl" &&
        f != "golden_helpers.jl" &&
        f != "test_constants.jl" &&
        f != "mpi_worker_tests.jl"
        error(
            "Test file $f is not listed in runtests.jl (neither unit nor integration group).",
        )
    end
end

files_to_run = String[]
if test_group == "all"
    append!(files_to_run, unit_tests)
    append!(files_to_run, integration_tests)
elseif test_group == "unit"
    append!(files_to_run, unit_tests)
elseif test_group == "integration"
    append!(files_to_run, integration_tests)
else
    error("Unknown EREBUS_TEST_GROUP: $test_group")
end

if isempty(files_to_run)
    error("No tests found to run!")
end

@testset verbose=true "Erebus.jl" begin
    for f in files_to_run
        include(f)
    end
end

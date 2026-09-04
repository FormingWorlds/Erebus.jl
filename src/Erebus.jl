module Erebus

using ArgParse
using Base.Threads
using Dates
using DocStringExtensions
using ExtendableSparse
using JLD2
using LinearAlgebra
using LinearSolve
using Logging
using ProgressMeter
using Random
using SparseArrays
using StaticArrays
using TimerOutputs
using TOML

export run_simulation
export Config, GridConfig, GeometryConfig, TimeConfig, SolverConfig,
       PoroelasticConfig, ThermalConfig, MaterialConfig, OutputConfig,
       SimulationConfig, default_config, load_config, validate_config, save_config
export Geometry, Physics, Particles, Numerics, Simulation

include("constants.jl")
# include("test_constants.jl")

if use_pardiso
    using Pardiso
else
    if Sys.isapple()
        using AppleAccelerate
    else
        # using MKL
    end
end

const to = TimerOutput()
const rgen = MersenneTwister(seed)

# Core modular components
include("config.jl")
include("geometry.jl")
include("physics.jl")
include("particles.jl")
include("numerics.jl")
include("simulation.jl")

# Submodule interfaces
module Geometry
    import ..Erebus: setup_staggered_grid_properties, setup_staggered_grid_properties_helpers,
                     grid_vector, dot4, grid_average, apply_insulating_boundary_conditions!
    export setup_staggered_grid_properties, setup_staggered_grid_properties_helpers,
           grid_vector, dot4, grid_average, apply_insulating_boundary_conditions!
end

module Physics
    import ..Erebus: distance, total, ktotal, kphi, ηᶠcur_inv_kᵠ, Q_radiogenic, etatotal_rocks,
                     calculate_radioactive_heating, compute_rhocpfluidm, compute_ksolidm, compute_kfluidm,
                     compute_Δtreaction, compute_gibbs_free_energy, compute_relative_enthalpy,
                     compute_reaction_constant, compute_thermodynamic_xfer!, perform_thermochemical_reaction!,
                     compute_shear_heating!, compute_adiabatic_heating!,
                     compute_drained_compressibility, compute_biot_willis_coefficient, compute_skempton_coefficient,
                     compute_rhofluid, compute_fluid_viscosity,
                     compute_hydrofracture_factor, compute_hydrofracture_permeability
    export distance, total, ktotal, kphi, ηᶠcur_inv_kᵠ, Q_radiogenic, etatotal_rocks,
           calculate_radioactive_heating, compute_rhocpfluidm, compute_ksolidm, compute_kfluidm,
           compute_Δtreaction, compute_gibbs_free_energy, compute_relative_enthalpy,
           compute_reaction_constant, compute_thermodynamic_xfer!, perform_thermochemical_reaction!,
           compute_shear_heating!, compute_adiabatic_heating!,
           compute_drained_compressibility, compute_biot_willis_coefficient, compute_skempton_coefficient,
           compute_rhofluid, compute_fluid_viscosity,
           compute_hydrofracture_factor, compute_hydrofracture_permeability
end

module Particles
    import ..Erebus: setup_marker_properties, setup_marker_properties_helpers, setup_marker_geometry_helpers,
                     define_markers!, compute_marker_properties!, update_marker_viscosity!,
                     setup_interpolated_properties, reset_interpolated_properties!, reset_thermochemical_properties!,
                     fix_weights, fix_distances, fix, reduce_add_3darray!, interpolate_add_to_grid!,
                     interpolate_to_marker!, interpolate_add_to_marker!, marker_to_basic_nodes!,
                     marker_to_vx_nodes!, marker_to_vy_nodes!, marker_to_p_nodes!,
                     molarfraction_marker_to_p_nodes!, update_p_nodes_melt_composition!,
                     compute_basic_node_properties!, compute_vx_node_properties!, compute_vy_node_properties!,
                     compute_p_node_properties!, compute_molarfraction!, add_vrk4,
                     compute_velocities!, compute_rotation_rate!, move_markers_rk4!, backtrace_pressures_rk4!,
                     update_marker_population_geometry!, replenish_markers!, apply_subgrid_stress_diffusion!,
                     update_marker_stress!, apply_subgrid_temperature_diffusion!, update_marker_temperature!,
                     update_marker_porosity!
    export setup_marker_properties, setup_marker_properties_helpers, setup_marker_geometry_helpers,
           define_markers!, compute_marker_properties!, update_marker_viscosity!,
           setup_interpolated_properties, reset_interpolated_properties!, reset_thermochemical_properties!,
           fix_weights, fix_distances, fix, reduce_add_3darray!, interpolate_add_to_grid!,
           interpolate_to_marker!, interpolate_add_to_marker!, marker_to_basic_nodes!,
           marker_to_vx_nodes!, marker_to_vy_nodes!, marker_to_p_nodes!,
           molarfraction_marker_to_p_nodes!, update_p_nodes_melt_composition!,
           compute_basic_node_properties!, compute_vx_node_properties!, compute_vy_node_properties!,
           compute_p_node_properties!, compute_molarfraction!, add_vrk4,
           compute_velocities!, compute_rotation_rate!, move_markers_rk4!, backtrace_pressures_rk4!,
           update_marker_population_geometry!, replenish_markers!, apply_subgrid_stress_diffusion!,
           update_marker_stress!, apply_subgrid_temperature_diffusion!, update_marker_temperature!,
           update_marker_porosity!
end

module Numerics
    import ..Erebus: setup_gravitational_lse, setup_hydromechanical_lse, setup_thermal_lse,
                     initialize_pardiso!, get_viscosities_stresses_density_gradients!,
                     assemble_gravitational_lse!, process_gravitational_solution!, compute_gravity_solution!,
                     assemble_hydromechanical_lse!, process_hydromechanical_solution!, recompute_bulk_viscosity!,
                     compute_Aϕ!, compute_fluid_velocities!, compute_displacement_timestep,
                     compute_stress_strainrate!, symmetrize_p_node_observables!, compute_nodal_adjustment!,
                     positive_max!, finalize_plastic_iteration_pass!, assemble_thermal_lse!,
                     perform_thermal_iterations!, finalize_thermochemical_iteration_pass,
                     compute_thermochemical_iteration_outcome
    export setup_gravitational_lse, setup_hydromechanical_lse, setup_thermal_lse,
           initialize_pardiso!, get_viscosities_stresses_density_gradients!,
           assemble_gravitational_lse!, process_gravitational_solution!, compute_gravity_solution!,
           assemble_hydromechanical_lse!, process_hydromechanical_solution!, recompute_bulk_viscosity!,
           compute_Aϕ!, compute_fluid_velocities!, compute_displacement_timestep,
           compute_stress_strainrate!, symmetrize_p_node_observables!, compute_nodal_adjustment!,
           positive_max!, finalize_plastic_iteration_pass!, assemble_thermal_lse!,
           perform_thermal_iterations!, finalize_thermochemical_iteration_pass,
           compute_thermochemical_iteration_outcome
end

module Simulation
    import ..Erebus: s_to_Ma, setup_dynamic_simulation_parameters, save_state, simulation_loop,
                     parse_commandline, run_simulation
    export s_to_Ma, setup_dynamic_simulation_parameters, save_state, simulation_loop,
           parse_commandline, run_simulation
end

module Config
    import ..Erebus: GridConfig, GeometryConfig, TimeConfig, SolverConfig,
                     PoroelasticConfig, ThermalConfig, MaterialConfig, OutputConfig,
                     SimulationConfig, default_config, load_config, validate_config, save_config
    export GridConfig, GeometryConfig, TimeConfig, SolverConfig,
           PoroelasticConfig, ThermalConfig, MaterialConfig, OutputConfig,
           SimulationConfig, default_config, load_config, validate_config, save_config
end

end # module Erebus

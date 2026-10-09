# Pre-allocated memory workspaces for segregation drift-flux solvers

"""
Reusable pre-allocated memory buffers for metal segregation calculations.

$(FIELDS)
"""
mutable struct MetalSegregationWorkspace
    Ny::Int
    Nx::Int
    M_fe_cell::Matrix{Float64}
    M_rock_markers::Matrix{Int}
    v_seg_cell::Matrix{Float64}
    phi_m_cell::Matrix{Float64}
    F_m_cell::Matrix{Float64}
    Xfem_cell::Matrix{Float64}
    g_acc_cell::Matrix{Float64}
    cap_cell::Matrix{Float64}
    phi_fe_cell::Matrix{Float64}
    T_cell::Matrix{Float64}
    drho_cell::Matrix{Float64}
    M_fe_H_cell::Matrix{Float64}
    M_fe_C_cell::Matrix{Float64}
    M_fe_N_cell::Matrix{Float64}
    M_fe_S_cell::Matrix{Float64}
    m_fe::Matrix{Float64}
    m_fe_H::Matrix{Float64}
    m_fe_C::Matrix{Float64}
    m_fe_N::Matrix{Float64}
    m_fe_S::Matrix{Float64}
    flux_H_x::Matrix{Float64}
    flux_H_y::Matrix{Float64}
    flux_C_x::Matrix{Float64}
    flux_C_y::Matrix{Float64}
    flux_N_x::Matrix{Float64}
    flux_N_y::Matrix{Float64}
    flux_S_x::Matrix{Float64}
    flux_S_y::Matrix{Float64}
    req_flux_x::Matrix{Float64}
    req_flux_y::Matrix{Float64}
    flux_x::Matrix{Float64}
    flux_y::Matrix{Float64}
    outflow_tot::Matrix{Float64}
    inflow_tot::Matrix{Float64}
    alpha_out::Matrix{Float64}
    alpha_in::Matrix{Float64}

    function MetalSegregationWorkspace(Ny::Integer, Nx::Integer; track_volatiles::Bool=true)
        Ny_val = Int(Ny)
        Nx_val = Int(Nx)
        return new(
            Ny_val,
            Nx_val,
            zeros(Float64, Ny_val, Nx_val),
            zeros(Int, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            track_volatiles ? zeros(Float64, Ny_val, Nx_val) : zeros(Float64, 0, 0),
            track_volatiles ? zeros(Float64, Ny_val, Nx_val) : zeros(Float64, 0, 0),
            track_volatiles ? zeros(Float64, Ny_val, Nx_val) : zeros(Float64, 0, 0),
            track_volatiles ? zeros(Float64, Ny_val, Nx_val) : zeros(Float64, 0, 0),
            zeros(Float64, Ny_val, Nx_val),
            track_volatiles ? zeros(Float64, Ny_val, Nx_val) : zeros(Float64, 0, 0),
            track_volatiles ? zeros(Float64, Ny_val, Nx_val) : zeros(Float64, 0, 0),
            track_volatiles ? zeros(Float64, Ny_val, Nx_val) : zeros(Float64, 0, 0),
            track_volatiles ? zeros(Float64, Ny_val, Nx_val) : zeros(Float64, 0, 0),
            track_volatiles ? zeros(Float64, Ny_val, Nx_val - 1) : zeros(Float64, 0, 0),
            track_volatiles ? zeros(Float64, Ny_val - 1, Nx_val) : zeros(Float64, 0, 0),
            track_volatiles ? zeros(Float64, Ny_val, Nx_val - 1) : zeros(Float64, 0, 0),
            track_volatiles ? zeros(Float64, Ny_val - 1, Nx_val) : zeros(Float64, 0, 0),
            track_volatiles ? zeros(Float64, Ny_val, Nx_val - 1) : zeros(Float64, 0, 0),
            track_volatiles ? zeros(Float64, Ny_val - 1, Nx_val) : zeros(Float64, 0, 0),
            track_volatiles ? zeros(Float64, Ny_val, Nx_val - 1) : zeros(Float64, 0, 0),
            track_volatiles ? zeros(Float64, Ny_val - 1, Nx_val) : zeros(Float64, 0, 0),
            zeros(Float64, Ny_val, Nx_val - 1),
            zeros(Float64, Ny_val - 1, Nx_val),
            zeros(Float64, Ny_val, Nx_val - 1),
            zeros(Float64, Ny_val - 1, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            ones(Float64, Ny_val, Nx_val),
            ones(Float64, Ny_val, Nx_val),
        )
    end
end

"""
Reusable pre-allocated memory buffers for silicate melt segregation calculations.

$(FIELDS)
"""
mutable struct MagmaSegregationWorkspace
    Ny::Int
    Nx::Int
    M_melt_cell::Matrix{Float64}
    M_rock_markers::Matrix{Int}
    v_seg_cell::Matrix{Float64}
    Fm_cell::Matrix{Float64}
    T_cell::Matrix{Float64}
    g_acc_cell::Matrix{Float64}
    cap_cell::Matrix{Float64}
    drho_cell::Matrix{Float64}
    m_melt::Matrix{Float64}
    req_flux_x::Matrix{Float64}
    req_flux_y::Matrix{Float64}
    flux_x::Matrix{Float64}
    flux_y::Matrix{Float64}
    outflow_tot::Matrix{Float64}
    inflow_tot::Matrix{Float64}
    alpha_out::Matrix{Float64}
    alpha_in::Matrix{Float64}
    P_comp::Matrix{Float64}
    div_v::Matrix{Float64}
    H_flux_x::Matrix{Float64}
    H_flux_y::Matrix{Float64}

    function MagmaSegregationWorkspace(Ny::Integer, Nx::Integer)
        Ny_val = Int(Ny)
        Nx_val = Int(Nx)
        return new(
            Ny_val,
            Nx_val,
            zeros(Float64, Ny_val, Nx_val),
            zeros(Int, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val - 1),
            zeros(Float64, Ny_val - 1, Nx_val),
            zeros(Float64, Ny_val, Nx_val - 1),
            zeros(Float64, Ny_val - 1, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            ones(Float64, Ny_val, Nx_val),
            ones(Float64, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val),
            zeros(Float64, Ny_val, Nx_val - 1),
            zeros(Float64, Ny_val - 1, Nx_val),
        )
    end
end

"""
Reusable pre-allocated memory workspace for hydromechanical linear system assembly and solution.

$(FIELDS)
"""
mutable struct HydromechanicalLSEWorkspace
    Ny1::Int
    Nx1::Int
    L::ExtendableSparseMatrix{Float64,Int64}
    is_initialized::Bool
    dof_per_node::Int
    pr_presolve::Matrix{Float64}
    pf_presolve::Matrix{Float64}
    rx_eff::Matrix{Float64}
    ry_eff::Matrix{Float64}
    rx_eff_prev::Matrix{Float64}
    ry_eff_prev::Matrix{Float64}
    pf_prev_iter::Matrix{Float64}

    function HydromechanicalLSEWorkspace(
        Ny1::Integer, Nx1::Integer; dof_per_node::Integer=6
    )
        Ny1_val = Int(Ny1)
        Nx1_val = Int(Nx1)
        dof_val = Int(dof_per_node)
        dim = Ny1_val * Nx1_val * dof_val
        L = ExtendableSparseMatrix(dim, dim)
        return new(
            Ny1_val,
            Nx1_val,
            L,
            false,
            dof_val,
            zeros(Float64, Ny1_val, Nx1_val),
            zeros(Float64, Ny1_val, Nx1_val),
            zeros(Float64, Ny1_val, Nx1_val),
            zeros(Float64, Ny1_val, Nx1_val),
            zeros(Float64, Ny1_val, Nx1_val),
            zeros(Float64, Ny1_val, Nx1_val),
            zeros(Float64, Ny1_val, Nx1_val),
        )
    end
end
function HydromechanicalLSEWorkspace(coords::GridCoordinates; dof_per_node::Integer=6)
    return HydromechanicalLSEWorkspace(coords.Ny1, coords.Nx1; dof_per_node=dof_per_node)
end

"""
Reusable pre-allocated memory workspace for thermal linear system assembly and solution.

$(FIELDS)
"""
mutable struct ThermalLSEWorkspace
    Ny1::Int
    Nx1::Int
    LT::ExtendableSparseMatrix{Float64,Int64}
    is_initialized::Bool

    function ThermalLSEWorkspace(Ny1::Integer, Nx1::Integer)
        Ny1_val = Int(Ny1)
        Nx1_val = Int(Nx1)
        dim = Ny1_val * Nx1_val
        LT = ExtendableSparseMatrix(dim, dim)
        return new(Ny1_val, Nx1_val, LT, false)
    end
end
ThermalLSEWorkspace(coords::GridCoordinates) = ThermalLSEWorkspace(coords.Ny1, coords.Nx1)

"""
Unified container for pre-allocated simulation linear system and marker workspaces.

$(FIELDS)
"""
mutable struct SimulationWorkspaces
    hydromech::HydromechanicalLSEWorkspace
    thermal::ThermalLSEWorkspace
    metal_segregation::Union{Nothing,MetalSegregationWorkspace}
    magma_segregation::Union{Nothing,MagmaSegregationWorkspace}
    p2m::Union{Nothing,P2MTiledWorkspace}
    thread_buffers::Union{Nothing,Vector{ThreadInterpolationBuffers}}
    interp_arrays::NTuple{33,Matrix{Float64}}
    hydromech_cache::Any
    thermal_cache::Any
    R::Vector{Float64}
    S::Vector{Float64}
    RT::Vector{Float64}
    ST::Vector{Float64}
    RP::Vector{Float64}
    SP::Vector{Float64}
    F_grav::Any
    mdis::Matrix{Float64}
    mnum::Matrix{Int}
    YERRNOD::Vector{Float64}
    fractured_cells::Matrix{Bool}
    fractured_cells_prev::Matrix{Bool}
    Xfe_bulk_step_start::Union{Nothing,Vector{Float64}}
    Xfem_step_start::Union{Nothing,Vector{Float64}}
    Fm_step_start::Union{Nothing,Vector{Float64}}
    F_extract_m_step_start::Union{Nothing,Vector{Float64}}
    Xfe_H_m_step_start::Union{Nothing,Vector{Float64}}
    Xfe_C_m_step_start::Union{Nothing,Vector{Float64}}
    Xfe_N_m_step_start::Union{Nothing,Vector{Float64}}
    Xfe_S_m_step_start::Union{Nothing,Vector{Float64}}
end

"""
Reset workspace arrays and solver caches when grid coordinates change.

$(SIGNATURES)
"""
function reset_workspaces_for_grid!(
    ws::SimulationWorkspaces,
    coords::GridCoordinates,
    cfg::SimulationConfig,
    marknum::Int,
)
    darcy_elim_val = cfg.solver.darcy_elimination
    dof_per_node_val = darcy_elim_val ? 4 : 6
    ws.R, ws.S = setup_hydromechanical_lse(coords; dof_per_node=dof_per_node_val)
    ws.hydromech = HydromechanicalLSEWorkspace(coords; dof_per_node=dof_per_node_val)
    ws.hydromech_cache = nothing

    ws.RT, ws.ST = setup_thermal_lse(coords)
    ws.thermal = ThermalLSEWorkspace(coords)
    ws.thermal_cache = nothing

    ws.RP, ws.SP = setup_gravitational_lse(coords)
    if cfg.geometry.gravity_mode === :poisson2d
        LP = assemble_gravitational_lse!(zeros(coords.Ny1, coords.Nx1), ws.RP; coords=coords)
        ws.F_grav = lu(LP.cscmatrix)
    else
        ws.F_grav = nothing
    end

    if ws.metal_segregation !== nothing
        ws.metal_segregation = MetalSegregationWorkspace(
            coords.Ny, coords.Nx; track_volatiles=cfg.metal_partition.active
        )
    end
    if ws.magma_segregation !== nothing
        ws.magma_segregation = MagmaSegregationWorkspace(coords.Ny, coords.Nx)
    end

    if cfg.solver.p2m_mode == :tiled
        ws.p2m = P2MTiledWorkspace(coords, marknum, cfg.solver.tile_size)
    elseif Threads.nthreads() > 1
        ws.thread_buffers = allocate_thread_interpolation_buffers(16, coords)
    end

    ws.mdis, ws.mnum = setup_marker_geometry_helpers(coords)
    ws.interp_arrays = setup_interpolated_properties(coords)
    ws.fractured_cells = zeros(Bool, coords.Ny, coords.Nx)
    ws.fractured_cells_prev = zeros(Bool, coords.Ny, coords.Nx)

    return ws
end

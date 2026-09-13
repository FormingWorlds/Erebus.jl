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

    function HydromechanicalLSEWorkspace(Ny1::Integer, Nx1::Integer)
        Ny1_val = Int(Ny1)
        Nx1_val = Int(Nx1)
        dim = Ny1_val * Nx1_val * 6
        L = ExtendableSparseMatrix(dim, dim)
        return new(Ny1_val, Nx1_val, L, false)
    end
end
function HydromechanicalLSEWorkspace(coords::GridCoordinates)
    return HydromechanicalLSEWorkspace(coords.Ny1, coords.Nx1)
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

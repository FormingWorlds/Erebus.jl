
"""
Set up gravitational linear system of equations structures.

$(SIGNATURES)

# Details

    - nothing

# Returns 

    - RP: gravitational linear system of equations: RHS vector
    - SP: gravitational linear system of equations: solution vector
"""
function setup_gravitational_lse(Nx1::Int=Nx1, Ny1::Int=Ny1)
    RP = Vector{Float64}(undef, Ny1*Nx1)
    SP = Vector{Float64}(undef, Ny1*Nx1)
    return RP, SP
end
function setup_gravitational_lse(coords::GridCoordinates)
    return setup_gravitational_lse(coords.Nx1, coords.Ny1)
end

"""
Set up hydromechanical linear system of equations structures.

$(SIGNATURES)

# Details

    - nothing

# Returns 

    - R: hydromechanical linear system of equations: RHS vector
    - S: hydromechanical linear system of equations: solution vector
"""
function setup_hydromechanical_lse(Nx1::Int=Nx1, Ny1::Int=Ny1)
    R = Vector{Float64}(undef, Ny1*Nx1*6)
    S = Vector{Float64}(undef, Ny1*Nx1*6)
    return R, S
end
function setup_hydromechanical_lse(coords::GridCoordinates)
    return setup_hydromechanical_lse(coords.Nx1, coords.Ny1)
end

"""
Set up thermal linear system of equations structures.

$(SIGNATURES)

# Details

    - nothing

# Returns 

    - RT: thermal linear system of equations: RHS vector
    - ST: thermal linear system of equations: solution vector
"""
function setup_thermal_lse(Nx1::Int=Nx1, Ny1::Int=Ny1)
    RT = Vector{Float64}(undef, Ny1*Nx1)
    ST = Vector{Float64}(undef, Ny1*Nx1)
    return RT, ST
end
setup_thermal_lse(coords::GridCoordinates) = setup_thermal_lse(coords.Nx1, coords.Ny1)

"""
Initialize `iparm` parameters of Pardiso MKL solver.

$(SIGNATURES)

# Details

    - ps: Instance of pardiso solver
    - iparms_dict: dictionary of iparm parameters

# Returns

    - nothing
"""
function initialize_pardiso!(pardiso_solver, iparms_dict)
    set_msglvl!(pardiso_solver, Pardiso.MESSAGE_LEVEL_OFF)
    set_matrixtype!(pardiso_solver, Pardiso.REAL_NONSYM)
    set_nprocs!(pardiso_solver, cache_kwargs.nprocs)
    for (i, v) in iparms_dict
        set_iparm!(pardiso_solver, i+1, v)
    end
    return set_phase!(pardiso_solver, Pardiso.ANALYSIS)
end


"""
Set up gravitational linear system of equations structures.

$(SIGNATURES)

# Details

    - nothing

# Returns 

    - RP: gravitational linear system of equations: RHS vector
    - SP: gravitational linear system of equations: solution vector
"""
function setup_gravitational_lse(Nx1::Int, Ny1::Int)
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
function setup_hydromechanical_lse(Nx1::Int, Ny1::Int; dof_per_node::Int=6)
    R = Vector{Float64}(undef, Ny1*Nx1*dof_per_node)
    S = Vector{Float64}(undef, Ny1*Nx1*dof_per_node)
    return R, S
end
function setup_hydromechanical_lse(coords::GridCoordinates; dof_per_node::Int=6)
    return setup_hydromechanical_lse(coords.Nx1, coords.Ny1; dof_per_node=dof_per_node)
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
function setup_thermal_lse(Nx1::Int, Ny1::Int)
    RT = Vector{Float64}(undef, Ny1*Nx1)
    ST = Vector{Float64}(undef, Ny1*Nx1)
    return RT, ST
end
setup_thermal_lse(coords::GridCoordinates) = setup_thermal_lse(coords.Nx1, coords.Ny1)

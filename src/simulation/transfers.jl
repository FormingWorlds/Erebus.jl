# src/simulation/transfers.jl
# TransferRecord struct for mass and volatile transfers between markers and reservoirs.

"""
    TransferRecord

Record of mass transfer between markers and planetary reservoirs.

# Fields
- `step::Int`: Simulation step index.
- `channel::Symbol`: Transfer channel (e.g. `:degassing`, `:venting`).
- `element::Symbol`: Elemental species (e.g. `:H`, `:C`, `:N`, `:S`).
- `m::Int`: Marker index.
- `x::Float64`: Marker horizontal coordinate [m].
- `y::Float64`: Marker vertical coordinate [m].
- `dM2::Float64`: 2D planar mass transfer increment [kg/m].
- `dM3::Float64`: 3D spherical mass transfer increment [kg].
"""
struct TransferRecord
    step::Int
    channel::Symbol
    element::Symbol
    m::Int
    x::Float64
    y::Float64
    dM2::Float64
    dM3::Float64
end

function TransferRecord(
    step::Integer,
    channel::Symbol,
    element::Symbol,
    m::Integer,
    x::Real,
    y::Real,
    dM2::Real,
    dM3::Real,
)
    return TransferRecord(
        Int(step),
        channel,
        element,
        Int(m),
        Float64(x),
        Float64(y),
        Float64(dM2),
        Float64(dM3),
    )
end

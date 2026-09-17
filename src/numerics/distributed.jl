# Distributed memory abstraction and domain decomposition interfaces.

using LinearAlgebra

"""
    AbstractGridTopology

Abstract supertype for 2D domain decomposition topologies.
"""
abstract type AbstractGridTopology end

"""
    SerialTopology <: AbstractGridTopology

Default topology representing single-node execution without domain decomposition.
"""
struct SerialTopology <: AbstractGridTopology end

"""
    is_distributed(topo::AbstractGridTopology) -> Bool
    is_distributed(x) -> Bool

Query whether `topo` or object `x` is configured for distributed execution.
"""
is_distributed(::SerialTopology) = false
is_distributed(::AbstractGridTopology) = false
is_distributed(::Any) = false

"""
    partition_indices(n::Integer, num_procs::Integer, proc_idx::Integer) -> UnitRange{Int}

Compute the 1-based local index range `istart:iend` for a 1D dimension of size `n`
partitioned among `num_procs` processes, for process `proc_idx` (1-based, 1 <= proc_idx <= num_procs).
"""
function partition_indices(
    n::Integer, num_procs::Integer, proc_idx::Integer
)::UnitRange{Int}
    num_procs >= 1 || throw(ArgumentError("num_procs must be >= 1, got $num_procs"))
    1 <= proc_idx <= num_procs ||
        throw(ArgumentError("proc_idx must be in 1:$num_procs, got $proc_idx"))
    base = div(n, num_procs)
    rem_val = mod(n, num_procs)
    start_idx = (proc_idx - 1) * base + min(proc_idx - 1, rem_val) + 1
    len = base + (proc_idx <= rem_val ? 1 : 0)
    return start_idx:(start_idx + len - 1)
end

# 2D Cartesian distributed grid topology representation.
struct DistributedGridTopology2D{C} <: AbstractGridTopology
    comm::C
    cart_comm::C
    rank::Int
    size::Int
    px::Int
    py::Int
    coords::Tuple{Int,Int}
    neighbors::NamedTuple{
        (:north, :south, :east, :west, :northeast, :northwest, :southeast, :southwest),
        NTuple{8,Int},
    }
    Ny_global::Int
    Nx_global::Int
    Ny1_global::Int
    Nx1_global::Int
    istart::Int
    iend::Int
    jstart::Int
    jend::Int
    ny_loc::Int
    nx_loc::Int
    ny1_loc::Int
    nx1_loc::Int
    halo_width::Int
    xmin_loc::Float64
    xmax_loc::Float64
    ymin_loc::Float64
    ymax_loc::Float64
    dx::Float64
    dy::Float64
    x_splits::Vector{Float64}
    y_splits::Vector{Float64}
end

is_distributed(::DistributedGridTopology2D) = true

# Preallocated memory buffers for zero-allocation asynchronous halo exchange.
mutable struct HaloBuffer{T,R}
    nvars::Int
    ny::Int
    nx::Int
    halo_width::Int
    send_east::Array{T,3}
    recv_east::Array{T,3}
    send_west::Array{T,3}
    recv_west::Array{T,3}
    send_north::Array{T,3}
    recv_north::Array{T,3}
    send_south::Array{T,3}
    recv_south::Array{T,3}
    reqs_x::Vector{R}
    reqs_y::Vector{R}
    work_x::Array{T,3}
end

# Distributed Stokes-Darcy operator with halo communication.
struct DistributedStokesDarcyOperator{T<:AbstractFloat,M<:AbstractMatrix{T},H<:HaloBuffer}
    local_op::MatrixFreeStokesDarcyOperator{T,M}
    topology::DistributedGridTopology2D
    halo_buffer::H
end

is_distributed(::DistributedStokesDarcyOperator) = true
Base.eltype(::DistributedStokesDarcyOperator{T}) where {T} = T
Base.size(op::DistributedStokesDarcyOperator) = Base.size(op.local_op)
Base.size(op::DistributedStokesDarcyOperator, d::Int) = Base.size(op.local_op, d)

# Distributed vector representation for Krylov solvers.
struct DistributedVector{T,V<:AbstractVector{T}} <: AbstractVector{T}
    local_vec::V
    topology::DistributedGridTopology2D
    global_length::Int
end

is_distributed(::DistributedVector) = true
Base.size(v::DistributedVector) = Base.size(v.local_vec)
Base.length(v::DistributedVector) = Base.length(v.local_vec)
Base.IndexStyle(::Type{<:DistributedVector}) = IndexLinear()
Base.getindex(v::DistributedVector, i::Int) = v.local_vec[i]
Base.setindex!(v::DistributedVector, val, i::Int) = (v.local_vec[i] = val)
function Base.similar(v::DistributedVector)
    return DistributedVector(similar(v.local_vec), v.topology, v.global_length)
end
function Base.similar(v::DistributedVector, ::Type{S}) where {S}
    return DistributedVector(similar(v.local_vec, S), v.topology, v.global_length)
end
function Base.copy(v::DistributedVector)
    return DistributedVector(copy(v.local_vec), v.topology, v.global_length)
end
Base.fill!(v::DistributedVector, val) = (fill!(v.local_vec, val); v)

"""
    exchange_halos!(topo::AbstractGridTopology, fields...; kwargs...)

Exchange ghost cell / halo margins between neighbor subdomains.
Serial topology performs a no-op.
"""
exchange_halos!(::SerialTopology, fields...; kwargs...) = nothing
exchange_halos!(::Any, fields...; kwargs...) = nothing

"""
    migrate_markers!(topo::AbstractGridTopology, args...; kwargs...)

Migrate marker particles that cross subdomain boundaries to the owning process ranks.
Serial topology performs a no-op and returns `(0, 0)` for `(migrated_out, migrated_in)`.
"""
migrate_markers!(::SerialTopology, args...; kwargs...) = (0, 0)
migrate_markers!(::Any, args...; kwargs...) = (0, 0)

"""
    distributed_dot(topo::AbstractGridTopology, x, y)

Compute global inner product for distributed subdomains.
Serial topology falls back to `LinearAlgebra.dot(x, y)`.
"""
distributed_dot(::SerialTopology, x, y) = dot(x, y)
distributed_dot(::Any, x, y) = dot(x, y)

"""
    distributed_norm(topo::AbstractGridTopology, x)

Compute global Euclidean norm for distributed subdomains.
Serial topology falls back to `LinearAlgebra.norm(x)`.
"""
distributed_norm(::SerialTopology, x) = norm(x)
distributed_norm(::Any, x) = norm(x)

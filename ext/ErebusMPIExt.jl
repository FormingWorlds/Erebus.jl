module ErebusMPIExt

using Erebus
using MPI
using LinearAlgebra

"""
    DistributedGridTopology2D(comm::MPI.Comm, Ny_global::Int, Nx_global::Int; kwargs...)

Organize MPI processes into a 2D Cartesian grid for domain decomposition.
"""
function Erebus.DistributedGridTopology2D(
    comm::MPI.Comm,
    Ny_global::Int,
    Nx_global::Int;
    px::Int=0,
    py::Int=0,
    halo_width::Int=1,
    xmin::Real=0.0,
    xmax::Real=100_000.0,
    ymin::Real=0.0,
    ymax::Real=100_000.0,
)
    halo_width >= 1 || throw(ArgumentError("halo_width must be >= 1, got $halo_width"))
    total_size = MPI.Comm_size(comm)

    # Automatic 2D factorization minimizing surface-to-volume ratio.
    if px <= 0 && py <= 0
        best_diff = typemax(Int)
        best_py = 1
        best_px = total_size
        for p_cand in 1:floor(Int, sqrt(total_size))
            if mod(total_size, p_cand) == 0
                q_cand = div(total_size, p_cand)
                if abs(q_cand - p_cand) < best_diff
                    best_diff = abs(q_cand - p_cand)
                    best_py = p_cand
                    best_px = q_cand
                end
            end
        end
        py = best_py
        px = best_px
    elseif px <= 0 && py > 0
        mod(total_size, py) == 0 ||
            throw(ArgumentError("comm size $total_size is not divisible by py ($py)"))
        px = div(total_size, py)
    elseif py <= 0 && px > 0
        mod(total_size, px) == 0 ||
            throw(ArgumentError("comm size $total_size is not divisible by px ($px)"))
        py = div(total_size, px)
    end
    px * py == total_size ||
        throw(ArgumentError("px * py ($px * $py) must match comm size $total_size"))
    px <= Nx_global || throw(ArgumentError("px ($px) cannot exceed Nx_global ($Nx_global)"))
    py <= Ny_global || throw(ArgumentError("py ($py) cannot exceed Ny_global ($Ny_global)"))

    # Dimension 0 is Y (rows), dimension 1 is X (cols)
    dims = Cint[py, px]
    periodic = Cint[0, 0]
    cart_comm = MPI.Cart_create(comm, dims; periodic=periodic, reorder=true)
    rank = MPI.Comm_rank(cart_comm)
    coords_c = MPI.Cart_coords(cart_comm, rank)
    cy = Int(coords_c[1])
    cx = Int(coords_c[2])

    (south_rank, north_rank) = MPI.Cart_shift(cart_comm, 0, 1)
    (west_rank, east_rank) = MPI.Cart_shift(cart_comm, 1, 1)

    # Compute diagonal corner neighbor ranks (cy+1 is North, cy-1 is South, cx+1 is East, cx-1 is West).
    nw = if (cy < py - 1 && cx > 0)
        MPI.Cart_rank(cart_comm, Cint[cy + 1, cx - 1])
    else
        MPI.PROC_NULL
    end
    ne = if (cy < py - 1 && cx < px - 1)
        MPI.Cart_rank(cart_comm, Cint[cy + 1, cx + 1])
    else
        MPI.PROC_NULL
    end
    sw = (cy > 0 && cx > 0) ? MPI.Cart_rank(cart_comm, Cint[cy - 1, cx - 1]) : MPI.PROC_NULL
    se = if (cy > 0 && cx < px - 1)
        MPI.Cart_rank(cart_comm, Cint[cy - 1, cx + 1])
    else
        MPI.PROC_NULL
    end

    neighbors = (
        north=Int(north_rank),
        south=Int(south_rank),
        east=Int(east_rank),
        west=Int(west_rank),
        northeast=Int(ne),
        northwest=Int(nw),
        southeast=Int(se),
        southwest=Int(sw),
    )

    Ny1_global = Ny_global + 1
    Nx1_global = Nx_global + 1
    dx = (Float64(xmax) - Float64(xmin)) / Float64(Nx_global)
    dy = (Float64(ymax) - Float64(ymin)) / Float64(Ny_global)

    # Partition domain cells.
    irange = partition_indices(Ny_global, py, cy + 1)
    jrange = partition_indices(Nx_global, px, cx + 1)

    istart = irange.start
    iend = irange.stop
    jstart = jrange.start
    jend = jrange.stop

    ny_loc = length(irange)
    nx_loc = length(jrange)
    ny1_loc = ny_loc + 1
    nx1_loc = nx_loc + 1

    # Compute spatial partitions across all ranks for O(1) marker routing.
    x_splits = zeros(Float64, px + 1)
    x_splits[1] = Float64(xmin)
    for p in 1:px
        jr = partition_indices(Nx_global, px, p)
        x_splits[p + 1] = Float64(xmin) + jr.stop * dx
    end

    y_splits = zeros(Float64, py + 1)
    y_splits[1] = Float64(ymin)
    for p in 1:py
        ir = partition_indices(Ny_global, py, p)
        y_splits[p + 1] = Float64(ymin) + ir.stop * dy
    end

    xmin_loc = x_splits[cx + 1]
    xmax_loc = x_splits[cx + 2]
    ymin_loc = y_splits[cy + 1]
    ymax_loc = y_splits[cy + 2]

    return DistributedGridTopology2D{typeof(comm)}(
        comm,
        cart_comm,
        rank,
        total_size,
        px,
        py,
        (cy, cx),
        neighbors,
        Ny_global,
        Nx_global,
        Ny1_global,
        Nx1_global,
        istart,
        iend,
        jstart,
        jend,
        ny_loc,
        nx_loc,
        ny1_loc,
        nx1_loc,
        halo_width,
        xmin_loc,
        xmax_loc,
        ymin_loc,
        ymax_loc,
        dx,
        dy,
        x_splits,
        y_splits,
    )
end

"""
    HaloBuffer{T}(nvars::Int, ny::Int, nx::Int; halo_width::Int=1) where {T}

Construct pre-allocated buffers for non-blocking halo exchanges.
"""
function Erebus.HaloBuffer{T}(nvars::Int, ny::Int, nx::Int; halo_width::Int=1) where {T}
    hw = halo_width
    return HaloBuffer{T,MPI.Request}(
        nvars,
        ny,
        nx,
        hw,
        zeros(T, nvars, ny, hw),
        zeros(T, nvars, ny, hw),
        zeros(T, nvars, ny, hw),
        zeros(T, nvars, ny, hw),
        zeros(T, nvars, hw, nx),
        zeros(T, nvars, hw, nx),
        zeros(T, nvars, hw, nx),
        zeros(T, nvars, hw, nx),
        MPI.Request[],
        MPI.Request[],
        zeros(T, nvars, ny, nx),
    )
end

"""
    exchange_halos!(topo::DistributedGridTopology2D, A; buffer=nothing)

Perform non-blocking halo exchange across 2D subdomain boundaries.
"""
function Erebus.exchange_halos!(
    topo::DistributedGridTopology2D,
    A::AbstractArray{T,3};
    buffer::Union{Nothing,HaloBuffer{T,MPI.Request}}=nothing,
) where {T}
    nvars, ny, nx = size(A)
    hw = topo.halo_width
    buf = if buffer !== nothing
        buffer
    else
        Erebus.HaloBuffer{T}(nvars, ny, nx; halo_width=hw)
    end

    # Step 1: Exchange along X (West and East)
    empty!(buf.reqs_x)
    if topo.neighbors.east != MPI.PROC_NULL
        @views buf.send_east .= A[:, :, (nx - 2 * hw + 1):(nx - hw)]
        push!(
            buf.reqs_x,
            MPI.Irecv!(buf.recv_east, topo.cart_comm; source=topo.neighbors.east, tag=101),
        )
        push!(
            buf.reqs_x,
            MPI.Isend(buf.send_east, topo.cart_comm; dest=topo.neighbors.east, tag=102),
        )
    end
    if topo.neighbors.west != MPI.PROC_NULL
        @views buf.send_west .= A[:, :, (hw + 1):(2 * hw)]
        push!(
            buf.reqs_x,
            MPI.Irecv!(buf.recv_west, topo.cart_comm; source=topo.neighbors.west, tag=102),
        )
        push!(
            buf.reqs_x,
            MPI.Isend(buf.send_west, topo.cart_comm; dest=topo.neighbors.west, tag=101),
        )
    end
    if !isempty(buf.reqs_x)
        MPI.Waitall(buf.reqs_x)
    end
    if topo.neighbors.east != MPI.PROC_NULL
        @views A[:, :, (nx - hw + 1):nx] .= buf.recv_east
    end
    if topo.neighbors.west != MPI.PROC_NULL
        @views A[:, :, 1:hw] .= buf.recv_west
    end

    # Step 2: Exchange along Y (South and North, includes corner halo data)
    empty!(buf.reqs_y)
    if topo.neighbors.north != MPI.PROC_NULL
        @views buf.send_north .= A[:, (ny - 2 * hw + 1):(ny - hw), :]
        push!(
            buf.reqs_y,
            MPI.Irecv!(
                buf.recv_north, topo.cart_comm; source=topo.neighbors.north, tag=201
            ),
        )
        push!(
            buf.reqs_y,
            MPI.Isend(buf.send_north, topo.cart_comm; dest=topo.neighbors.north, tag=202),
        )
    end
    if topo.neighbors.south != MPI.PROC_NULL
        @views buf.send_south .= A[:, (hw + 1):(2 * hw), :]
        push!(
            buf.reqs_y,
            MPI.Irecv!(
                buf.recv_south, topo.cart_comm; source=topo.neighbors.south, tag=202
            ),
        )
        push!(
            buf.reqs_y,
            MPI.Isend(buf.send_south, topo.cart_comm; dest=topo.neighbors.south, tag=201),
        )
    end
    if !isempty(buf.reqs_y)
        MPI.Waitall(buf.reqs_y)
    end
    if topo.neighbors.north != MPI.PROC_NULL
        @views A[:, (ny - hw + 1):ny, :] .= buf.recv_north
    end
    if topo.neighbors.south != MPI.PROC_NULL
        @views A[:, 1:hw, :] .= buf.recv_south
    end

    return A
end

function Erebus.exchange_halos!(
    topo::DistributedGridTopology2D,
    A::AbstractMatrix{T};
    buffer::Union{Nothing,HaloBuffer{T,MPI.Request}}=nothing,
) where {T}
    A_3d = reshape(A, (1, size(A, 1), size(A, 2)))
    Erebus.exchange_halos!(topo, A_3d; buffer=buffer)
    return A
end

"""
    DistributedStokesDarcyOperator(local_op, topology)

Wrap a local `MatrixFreeStokesDarcyOperator` with halo communication for distributed execution.
"""
function Erebus.DistributedStokesDarcyOperator(
    local_op::MatrixFreeStokesDarcyOperator{T,M}, topology::DistributedGridTopology2D
) where {T,M}
    bc_north = (topology.neighbors.north == MPI.PROC_NULL)
    bc_south = (topology.neighbors.south == MPI.PROC_NULL)
    bc_west = (topology.neighbors.west == MPI.PROC_NULL)
    bc_east = (topology.neighbors.east == MPI.PROC_NULL)

    configured_local_op = MatrixFreeStokesDarcyOperator{T,M}(
        local_op.Ny1,
        local_op.Nx1,
        local_op.dx,
        local_op.dy,
        local_op.Nx_val,
        local_op.Ny_val,
        local_op.ETA,
        local_op.ETAP,
        local_op.GGG,
        local_op.GGGP,
        local_op.RHOX,
        local_op.RHOY,
        local_op.RHOFX,
        local_op.RHOFY,
        local_op.RX,
        local_op.RY,
        local_op.ETAPHI,
        local_op.BETAPHI,
        local_op.PHI,
        local_op.gx,
        local_op.gy,
        local_op.dt,
        local_op.betasolid,
        local_op.betafluid,
        local_op.phimin,
        local_op.phimax,
        local_op.Kcont,
        local_op.bctop,
        local_op.bcbottom,
        local_op.bcleft,
        local_op.bcright,
        bc_north,
        bc_south,
        bc_west,
        bc_east,
    )
    buf = Erebus.HaloBuffer{T}(
        4, configured_local_op.Ny1, configured_local_op.Nx1; halo_width=topology.halo_width
    )
    return DistributedStokesDarcyOperator{T,M,typeof(buf)}(
        configured_local_op, topology, buf
    )
end

function LinearAlgebra.mul!(
    y::AbstractVector, op::DistributedStokesDarcyOperator, x::AbstractVector
)
    copyto!(op.halo_buffer.work_x, x)
    Erebus.exchange_halos!(op.topology, op.halo_buffer.work_x; buffer=op.halo_buffer)
    LinearAlgebra.mul!(y, op.local_op, vec(op.halo_buffer.work_x))
    return y
end

function LinearAlgebra.mul!(
    y::AbstractVector,
    op::DistributedStokesDarcyOperator,
    x::AbstractVector,
    alpha::Number,
    beta::Number,
)
    if iszero(beta)
        LinearAlgebra.mul!(y, op, x)
        if !isone(alpha)
            y .*= alpha
        end
    else
        tmp = similar(y)
        LinearAlgebra.mul!(tmp, op, x)
        y .= alpha .* tmp .+ beta .* y
    end
    return y
end

function Erebus.distributed_dot(
    topo::DistributedGridTopology2D, x::AbstractVector, y::AbstractVector
)
    local_sum = dot(x, y)
    return MPI.Allreduce(local_sum, +, topo.cart_comm)
end

function Erebus.distributed_norm(topo::DistributedGridTopology2D, x::AbstractVector)
    local_sq = dot(x, x)
    global_sq = MPI.Allreduce(local_sq, +, topo.cart_comm)
    return sqrt(max(zero(global_sq), global_sq))
end

function Erebus.distributed_dot(
    topo::DistributedGridTopology2D, A::AbstractMatrix{T}, B::AbstractMatrix{T}
) where {T}
    hw = topo.halo_width
    ny, nx = size(A)
    local_sum = zero(T)
    @inbounds for j in (hw + 1):(nx - hw), i in (hw + 1):(ny - hw)
        local_sum += A[i, j] * B[i, j]
    end
    return MPI.Allreduce(local_sum, +, topo.cart_comm)
end

function Erebus.distributed_norm(
    topo::DistributedGridTopology2D, A::AbstractMatrix{T}
) where {T}
    return sqrt(max(zero(T), Erebus.distributed_dot(topo, A, A)))
end

function Erebus.distributed_dot(
    topo::DistributedGridTopology2D, A::AbstractArray{T,3}, B::AbstractArray{T,3}
) where {T}
    hw = topo.halo_width
    nvars, ny, nx = size(A)
    local_sum = zero(T)
    @inbounds for j in (hw + 1):(nx - hw), i in (hw + 1):(ny - hw), v in 1:nvars
        local_sum += A[v, i, j] * B[v, i, j]
    end
    return MPI.Allreduce(local_sum, +, topo.cart_comm)
end

function Erebus.distributed_norm(
    topo::DistributedGridTopology2D, A::AbstractArray{T,3}
) where {T}
    return sqrt(max(zero(T), Erebus.distributed_dot(topo, A, A)))
end

function LinearAlgebra.dot(x::DistributedVector, y::DistributedVector)
    return Erebus.distributed_dot(x.topology, x.local_vec, y.local_vec)
end

function LinearAlgebra.norm(x::DistributedVector)
    return Erebus.distributed_norm(x.topology, x.local_vec)
end

# Route marker particle (x, y) to destination MPI rank.
function _locate_marker_rank(x::Real, y::Real, topo::DistributedGridTopology2D)::Int
    px = topo.px
    py = topo.py
    xs = topo.x_splits
    ys = topo.y_splits

    jx = 1
    if x <= xs[1]
        jx = 1
    elseif x >= xs[end]
        jx = px
    else
        jx = clamp(searchsortedlast(xs, Float64(x)), 1, px)
    end

    iy = 1
    if y <= ys[1]
        iy = 1
    elseif y >= ys[end]
        iy = py
    else
        iy = clamp(searchsortedlast(ys, Float64(y)), 1, py)
    end

    return Int(MPI.Cart_rank(topo.cart_comm, Cint[iy - 1, jx - 1]))
end

"""
    migrate_markers!(topo::DistributedGridTopology2D, xm, ym, property_arrays...)

Migrate markers that leave local spatial boundaries to destination ranks via Alltoall.
Returns `(migrated_out, migrated_in)`.
"""
function Erebus.migrate_markers!(
    topo::DistributedGridTopology2D, xm::Vector{T}, ym::Vector{T}, properties::Vector...;
) where {T<:Real}
    N = length(xm)
    sz = topo.size
    my_rank = topo.rank

    # 1. Classify markers by destination rank.
    target_ranks = Vector{Int}(undef, N)
    send_counts = zeros(Cint, sz)
    num_outgoing = 0

    for m in 1:N
        dest = _locate_marker_rank(xm[m], ym[m], topo)
        target_ranks[m] = dest
        if dest != my_rank
            send_counts[dest + 1] += 1
            num_outgoing += 1
        end
    end

    # 2. Exchange counts across all ranks.
    recv_counts = zeros(Cint, sz)
    MPI.Alltoall!(MPI.UBuffer(send_counts, 1), MPI.UBuffer(recv_counts, 1), topo.cart_comm)
    num_incoming = Int(sum(recv_counts))

    # Collective synchronization: check globally if any rank has outgoing markers
    global_outgoing = MPI.Allreduce(num_outgoing, +, topo.cart_comm)
    if global_outgoing == 0
        return (0, 0)
    end

    # 3. Group outgoing indices by destination rank.
    grouped_indices = [Int[] for _ in 1:sz]
    for m in 1:N
        dest = target_ranks[m]
        if dest != my_rank
            push!(grouped_indices[dest + 1], m)
        end
    end

    keep_mask = [target_ranks[m] == my_rank for m in 1:N]

    # 4. Migrate coordinate arrays xm and ym.
    send_xm = Vector{T}(undef, num_outgoing)
    send_ym = Vector{T}(undef, num_outgoing)
    pos = 1
    for r in 1:sz
        for idx in grouped_indices[r]
            send_xm[pos] = xm[idx]
            send_ym[pos] = ym[idx]
            pos += 1
        end
    end

    recv_xm = Vector{T}(undef, num_incoming)
    recv_ym = Vector{T}(undef, num_incoming)
    MPI.Alltoallv!(
        MPI.VBuffer(send_xm, send_counts), MPI.VBuffer(recv_xm, recv_counts), topo.cart_comm
    )
    MPI.Alltoallv!(
        MPI.VBuffer(send_ym, send_counts), MPI.VBuffer(recv_ym, recv_counts), topo.cart_comm
    )

    xm_kept = xm[keep_mask]
    ym_kept = ym[keep_mask]
    empty!(xm)
    empty!(ym)
    append!(xm, xm_kept)
    append!(xm, recv_xm)
    append!(ym, ym_kept)
    append!(ym, recv_ym)

    # 5. Migrate all associated property arrays.
    for prop in properties
        PType = eltype(prop)
        send_prop = Vector{PType}(undef, num_outgoing)
        p_pos = 1
        for r in 1:sz
            for idx in grouped_indices[r]
                send_prop[p_pos] = prop[idx]
                p_pos += 1
            end
        end
        recv_prop = Vector{PType}(undef, num_incoming)
        MPI.Alltoallv!(
            MPI.VBuffer(send_prop, send_counts),
            MPI.VBuffer(recv_prop, recv_counts),
            topo.cart_comm,
        )
        prop_kept = prop[keep_mask]
        empty!(prop)
        append!(prop, prop_kept)
        append!(prop, recv_prop)
    end

    return (num_outgoing, num_incoming)
end

function Erebus.migrate_markers!(topo::DistributedGridTopology2D, markers::NamedTuple;)
    hasproperty(markers, :xm) && hasproperty(markers, :ym) ||
        throw(ArgumentError("markers NamedTuple must contain :xm and :ym"))
    xm = markers.xm
    ym = markers.ym
    other_props = [getproperty(markers, k) for k in keys(markers) if k !== :xm && k !== :ym]
    return Erebus.migrate_markers!(topo, xm, ym, other_props...)
end

end # module ErebusMPIExt

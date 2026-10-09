# Concrete state and grid containers for Erebus planetesimal simulations.

"""
Exception thrown when checkpoint validation, serialization, or restoration fails.
"""
struct CheckpointError <: Exception
    msg::String
end

Base.showerror(io::IO, e::CheckpointError) = print(io, "CheckpointError: ", e.msg)

"""
Unified container for all 80 staggered Eulerian mesh arrays.
"""
struct GridArrays
    ETA::Matrix{Float64}
    ETA0::Matrix{Float64}
    GGG::Matrix{Float64}
    EXY::Matrix{Float64}
    SXY::Matrix{Float64}
    SXY0::Matrix{Float64}
    wyx::Matrix{Float64}
    COH::Matrix{Float64}
    TEN::Matrix{Float64}
    FRI::Matrix{Float64}
    YNY::Union{Matrix{Bool},Matrix{Float64}}
    RHOX::Matrix{Float64}
    RHOFX::Matrix{Float64}
    KX::Matrix{Float64}
    PHIX::Matrix{Float64}
    vx::Matrix{Float64}
    vxf::Matrix{Float64}
    RX::Matrix{Float64}
    qxD::Matrix{Float64}
    gx::Matrix{Float64}
    RHOY::Matrix{Float64}
    RHOFY::Matrix{Float64}
    KY::Matrix{Float64}
    PHIY::Matrix{Float64}
    vy::Matrix{Float64}
    vyf::Matrix{Float64}
    RY::Matrix{Float64}
    qyD::Matrix{Float64}
    gy::Matrix{Float64}
    RHO::Matrix{Float64}
    RHOCP::Matrix{Float64}
    ALPHA::Matrix{Float64}
    ALPHAF::Matrix{Float64}
    HR::Matrix{Float64}
    HA::Matrix{Float64}
    HS::Matrix{Float64}
    ETAP::Matrix{Float64}
    GGGP::Matrix{Float64}
    EXX::Matrix{Float64}
    SXX::Matrix{Float64}
    SXX0::Matrix{Float64}
    tk1::Matrix{Float64}
    tk2::Matrix{Float64}
    DT::Matrix{Float64}
    DT0::Matrix{Float64}
    vxp::Matrix{Float64}
    vyp::Matrix{Float64}
    vxpf::Matrix{Float64}
    vypf::Matrix{Float64}
    pr::Matrix{Float64}
    pf::Matrix{Float64}
    ps::Matrix{Float64}
    pr0::Matrix{Float64}
    pf0::Matrix{Float64}
    ps0::Matrix{Float64}
    ETAPHI::Matrix{Float64}
    BETAPHI::Matrix{Float64}
    PHI::Matrix{Float64}
    APHI::Matrix{Float64}
    FI::Matrix{Float64}
    DMP::Matrix{Float64}
    DHP::Matrix{Float64}
    XWS::Matrix{Float64}
    ETA5::Matrix{Float64}
    ETA00::Matrix{Float64}
    YNY5::Union{Matrix{Bool},Matrix{Float64}}
    YNY00::Union{Matrix{Bool},Matrix{Float64}}
    YNY_inv_ETA::Matrix{Float64}
    DSXY::Matrix{Float64}
    DSY::Matrix{Float64}
    EII::Matrix{Float64}
    SII::Matrix{Float64}
    DSXX::Matrix{Float64}
    tk0::Matrix{Float64}
    DQPF::Matrix{Float64}
    DQPFSUM::Matrix{Float64}
    S_vent_grid::Matrix{Float64}
    Q_lat_grid::Matrix{Float64}
    Q_seg_grid::Matrix{Float64}
    Q_metric::Union{Nothing,Matrix{Float64}}
end

Base.keys(::GridArrays) = fieldnames(GridArrays)
Base.pairs(g::GridArrays) = (fn => getfield(g, fn) for fn in fieldnames(GridArrays))
Base.haskey(::GridArrays, s::Symbol) = hasfield(GridArrays, s)
Base.haskey(::GridArrays, s::AbstractString) = hasfield(GridArrays, Symbol(s))
Base.getindex(g::GridArrays, s::Symbol) = getfield(g, s)
Base.getindex(g::GridArrays, s::AbstractString) = getfield(g, Symbol(s))
Base.propertynames(::GridArrays) = fieldnames(GridArrays)

function Base.copy(g::GridArrays)
    vals = Any[]
    for fn in fieldnames(GridArrays)
        v = getfield(g, fn)
        push!(vals, v === nothing ? nothing : copy(v))
    end
    return GridArrays(vals...)
end

"""
Copy all grid arrays in-place from source to destination.

$(SIGNATURES)
"""
@generated function copy_grid_arrays!(dst::GridArrays, src::GridArrays)
    exprs = Expr[]
    for fn in fieldnames(GridArrays)
        push!(exprs, quote
            let v_src = getfield(src, $(QuoteNode(fn))),
                v_dst = getfield(dst, $(QuoteNode(fn)))
                if v_src !== nothing && v_dst !== nothing
                    copyto!(v_dst, v_src)
                end
            end
        end)
    end
    return Expr(:block, exprs..., :(return dst))
end

"""
Diagnostic and physical mass accumulators across global simulation evolution.
"""
mutable struct SimulationAccumulators
    M_vent_total::Float64
    M_vent_H2O_total::Float64
    M_vent_C_total::Float64
    M_vent_N_total::Float64
    M_vent_S_total::Float64
    M_atm_total::Float64
    M_escaped_total::Float64
    P_amb::Float64
    rplanet::Float64
    rcore::Float64
    telescope_level::Int
    M_accreted_total::Float64
    M_planet_val::Float64
    xcenter::Float64
    ycenter::Float64
    max_v_seg_prev::Float64
    M_atm_species::Union{Nothing,Dict{Symbol,Float64}}
    M_escaped_species::Union{Nothing,Dict{Symbol,Float64}}
    core_budgets::Any
    regional_mineral_modes::Any
end

function Base.copy(a::SimulationAccumulators)
    return SimulationAccumulators(
        a.M_vent_total,
        a.M_vent_H2O_total,
        a.M_vent_C_total,
        a.M_vent_N_total,
        a.M_vent_S_total,
        a.M_atm_total,
        a.M_escaped_total,
        a.P_amb,
        a.rplanet,
        a.rcore,
        a.telescope_level,
        a.M_accreted_total,
        a.M_planet_val,
        a.xcenter,
        a.ycenter,
        a.max_v_seg_prev,
        a.M_atm_species === nothing ? nothing : copy(a.M_atm_species),
        a.M_escaped_species === nothing ? nothing : copy(a.M_escaped_species),
        a.core_budgets === nothing ? nothing : deepcopy(a.core_budgets),
        a.regional_mineral_modes === nothing ? nothing : deepcopy(a.regional_mineral_modes),
    )
end

"""
Complete physical, numerical, and diagnostic simulation state.
"""
mutable struct SimulationState{G<:NamedTuple,R<:Random.AbstractRNG,T<:AbstractVector}
    grids::GridArrays
    markers::MarkerArrays{G}
    accumulators::SimulationAccumulators
    transfers::T
    atm::Union{Nothing,AtmosphereState}
    rng::R
    timer::TimerOutput
    timestep::Int
    dt::Float64
    timesum::Float64
end

function Base.propertynames(::SimulationState)
    return (
        fieldnames(SimulationState)...,
        fieldnames(SimulationAccumulators)...,
        fieldnames(GridArrays)...,
    )
end

function Base.getproperty(s::SimulationState, sym::Symbol)
    if sym in fieldnames(SimulationState)
        return getfield(s, sym)
    elseif hasfield(SimulationAccumulators, sym)
        return getfield(getfield(s, :accumulators), sym)
    elseif hasfield(GridArrays, sym)
        return getfield(getfield(s, :grids), sym)
    else
        error("type SimulationState has no field or delegated property '$sym'")
    end
end

function Base.copy(s::SimulationState)
    return SimulationState(
        copy(s.grids),
        copy(s.markers),
        copy(s.accumulators),
        deepcopy(s.transfers),
        s.atm === nothing ? nothing : copy(s.atm),
        copy(s.rng),
        copy(s.timer),
        s.timestep,
        s.dt,
        s.timesum,
    )
end

function Base.keys(::SimulationState)
    return (
        fieldnames(SimulationState)...,
        fieldnames(SimulationAccumulators)...,
        fieldnames(GridArrays)...,
        :transfer_log,
        :S_vent,
    )
end

function Base.haskey(s::SimulationState, sym::Symbol)
    return (
        hasfield(SimulationState, sym) ||
        hasfield(SimulationAccumulators, sym) ||
        hasfield(GridArrays, sym) ||
        sym === :transfer_log ||
        sym === :S_vent
    )
end
Base.haskey(s::SimulationState, str::AbstractString) = haskey(s, Symbol(str))

function Base.getindex(s::SimulationState, sym::Symbol)
    if sym === :transfer_log
        return getfield(s, :transfers)
    elseif sym === :S_vent
        return getfield(getfield(s, :grids), :S_vent_grid)
    elseif hasfield(SimulationState, sym)
        return getfield(s, sym)
    elseif hasfield(SimulationAccumulators, sym)
        return getfield(getfield(s, :accumulators), sym)
    elseif hasfield(GridArrays, sym)
        return getfield(getfield(s, :grids), sym)
    else
        throw(KeyError(sym))
    end
end
Base.getindex(s::SimulationState, str::AbstractString) = s[Symbol(str)]

Base.pairs(s::SimulationState) = (k => s[k] for k in keys(s))

function Base.NamedTuple(s::SimulationState)
    return (;
        markers=s.markers,
        grids=s.grids,
        atm=s.atm,
        transfers=s.transfers,
        timesum=s.timesum,
        dt=s.dt,
        timestep=s.timestep,
    )
end

function Base.setproperty!(s::SimulationState, sym::Symbol, val)
    if sym in fieldnames(SimulationState)
        return setfield!(s, sym, val)
    elseif hasfield(SimulationAccumulators, sym)
        return setfield!(getfield(s, :accumulators), sym, val)
    else
        error("type SimulationState has no field or delegated mutable property '$sym'")
    end
end

"""
Snapshot of simulation state at timestep start for plastic retry rollback.

$(FIELDS)
"""
struct StepSnapshot{G<:NamedTuple}
    grids::GridArrays
    markers::MarkerArrays{G}
    scalars::NamedTuple
    YERRNOD::Vector{Float64}
    coords::GridCoordinates
end

"""
Capture simulation state snapshot at timestep start.

$(SIGNATURES)
"""
function snapshot_step_state(
    state::SimulationState, coords::GridCoordinates, ws::SimulationWorkspaces
)
    return StepSnapshot(
        copy(state.grids),
        copy(state.markers),
        (;
            marknum=length(state.markers),
            M_planet_val=state.accumulators.M_planet_val,
            M_accreted_total=state.accumulators.M_accreted_total,
            telescope_level=state.accumulators.telescope_level,
            rplanet_val=state.accumulators.rplanet,
            xcenter_val=state.accumulators.xcenter,
            ycenter_val=state.accumulators.ycenter,
        ),
        copy(ws.YERRNOD),
        coords,
    )
end

"""
Snapshot step state generic fallback for NamedTuple or dict state representations.

$(SIGNATURES)
"""
function snapshot_step_state(state)
    if hasproperty(state, :markers) && state.markers isa MarkerArrays
        return (;
            markers=copy(state.markers),
            arrays=if hasproperty(state, :arrays)
                deepcopy(state.arrays)
            else
                (hasproperty(state, :grids) ? deepcopy(state.grids) : (;))
            end,
            scalars=deepcopy(state.scalars),
        )
    else
        return deepcopy(state)
    end
end

"""
Restore simulation state in-place from a step-start snapshot.

$(SIGNATURES)
"""
function restore_step_state!(
    state::SimulationState,
    coords_ref::Ref{GridCoordinates},
    ws::SimulationWorkspaces,
    snapshot::StepSnapshot,
)
    copy_grid_arrays!(state.grids, snapshot.grids)
    restore_marker_arrays!(state.markers, snapshot.markers)
    acc = state.accumulators
    acc.M_planet_val = snapshot.scalars.M_planet_val
    acc.M_accreted_total = snapshot.scalars.M_accreted_total
    acc.telescope_level = snapshot.scalars.telescope_level
    acc.rplanet = snapshot.scalars.rplanet_val
    acc.xcenter = snapshot.scalars.xcenter_val
    acc.ycenter = snapshot.scalars.ycenter_val
    copyto!(ws.YERRNOD, snapshot.YERRNOD)
    coords_ref[] = snapshot.coords
    return state
end

"""
Restore simulation state generic fallback for NamedTuple representations.

$(SIGNATURES)
"""
function restore_step_state!(target, source)
    for (k, v) in pairs(source)
        if v isa AbstractArray && haskey(target, k)
            tgt = target[k]
            if tgt isa AbstractArray
                if tgt isa Vector && length(tgt) != length(v)
                    resize!(tgt, length(v))
                end
                copyto!(tgt, v)
            end
        end
    end
    return target
end

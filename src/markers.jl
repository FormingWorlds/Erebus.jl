# Concrete marker groups and unified MarkerArrays container.

"""
Foundational marker arrays present in all planetesimal simulations.
"""
struct CoreGroup
    xm::Vector{Float64}
    ym::Vector{Float64}
    w3d_m::Vector{Float64}
    tm::Vector{Int}
    tkm::Vector{Float64}
    phim::Vector{Float64}
    phinewm::Vector{Float64}
    pfm0::Vector{Float64}
    XWsolidm::Vector{Float64}
    XWsolidm0::Vector{Float64}
    Fm::Vector{Float64}
    etavpm::Vector{Float64}
    sxxm::Vector{Float64}
    sxym::Vector{Float64}
    inv_gggtotalm::Vector{Float64}
    fricttotalm::Vector{Float64}
    cohestotalm::Vector{Float64}
    tenstotalm::Vector{Float64}
    rhototalm::Vector{Float64}
    rhocptotalm::Vector{Float64}
    etatotalm::Vector{Float64}
    hrtotalm::Vector{Float64}
    ktotalm::Vector{Float64}
    tkm_rhocptotalm::Vector{Float64}
    etafluidcur_inv_kphim::Vector{Float64}
    rhofluidcur::Vector{Float64}
    alphasolidcur::Vector{Float64}
    alphafluidcur::Vector{Float64}
end

"""
Marker arrays for core segregation and metal-silicate partitioning.
"""
struct MetalGroup
    Xfem::Vector{Float64}
    Xfem0::Vector{Float64}
    Xfe_bulk::Vector{Float64}
    Xfe_H_m::Vector{Float64}
    Xfe_C_m::Vector{Float64}
    Xfe_N_m::Vector{Float64}
    Xfe_S_m::Vector{Float64}
end

"""
Marker arrays for volatile speciation, transport, and magma degassing.
"""
struct VolatilesGroup
    XH2Om::Vector{Float64}
    XCm::Vector{Float64}
    XNm::Vector{Float64}
    XSm::Vector{Float64}
    X_graphite_m::Vector{Float64}
    F_extract_m::Vector{Float64}
end

"""
Marker arrays for oxygen fugacity buffering and Evans (2012) electron accounting.
"""
struct RedoxGroup
    nFe0_m::Vector{Float64}
    nFe2_m::Vector{Float64}
    nFe3_m::Vector{Float64}
    deltaIW_m::Vector{Float64}
    nC_graphite_m::Vector{Float64}
    nCO_m::Vector{Float64}
    nCO2_m::Vector{Float64}
    nCH4_m::Vector{Float64}
end

"""
Marker arrays for multi-species ices and refractory phases.
"""
struct HcnspoGroup
    X_ice_H2O_m::Vector{Float64}
    X_ice_NH3_m::Vector{Float64}
    X_ice_CO2_m::Vector{Float64}
    X_ice_CO_m::Vector{Float64}
    X_ice_CH4_m::Vector{Float64}
    X_ice_N2_m::Vector{Float64}
    X_ice_H2S_m::Vector{Float64}
    X_ice_PH3_m::Vector{Float64}
    X_refr_C_m::Vector{Float64}
    X_refr_S_m::Vector{Float64}
    X_refr_N_m::Vector{Float64}
    X_refr_P_m::Vector{Float64}
    X_refr_H_m::Vector{Float64}
end

"""
Marker arrays for sub-eutectic accessory mineral exsolution modes.
"""
struct PhaseGroup
    Xmin_troilite_m::Vector{Float64}
    Xmin_schreibersite_m::Vector{Float64}
    Xmin_cohenite_m::Vector{Float64}
    Xmin_graphite_m::Vector{Float64}
    Xmin_nitride_m::Vector{Float64}
    Xmin_metal_matrix_m::Vector{Float64}
end

"""
Marker arrays for pebble accretion timing.
"""
struct AccretionGroup
    t_accreted::Vector{Float64}
end

function Base.NamedTuple(
    g::T
) where {
    T<:Union{
        CoreGroup,MetalGroup,VolatilesGroup,RedoxGroup,HcnspoGroup,PhaseGroup,AccretionGroup
    },
}
    return NamedTuple{fieldnames(T)}(Tuple(getfield(g, name) for name in fieldnames(T)))
end

function Base.merge(
    a::NamedTuple,
    g::Union{
        CoreGroup,MetalGroup,VolatilesGroup,RedoxGroup,HcnspoGroup,PhaseGroup,AccretionGroup
    },
)
    return merge(a, NamedTuple(g))
end

function Base.keys(
    g::Union{
        CoreGroup,MetalGroup,VolatilesGroup,RedoxGroup,HcnspoGroup,PhaseGroup,AccretionGroup
    },
)
    return fieldnames(typeof(g))
end

function Base.values(
    g::Union{
        CoreGroup,MetalGroup,VolatilesGroup,RedoxGroup,HcnspoGroup,PhaseGroup,AccretionGroup
    },
)
    return Tuple(getfield(g, fn) for fn in fieldnames(typeof(g)))
end

function Base.iterate(
    g::Union{
        CoreGroup,MetalGroup,VolatilesGroup,RedoxGroup,HcnspoGroup,PhaseGroup,AccretionGroup
    },
    state...,
)
    return iterate(values(g), state...)
end

"""
Unified parametric marker container holding core and active optional groups.
"""
struct MarkerArrays{G<:NamedTuple}
    core::CoreGroup
    groups::G
end

Base.length(m::MarkerArrays) = length(m.core.xm)

@inline function Base.getproperty(m::MarkerArrays, sym::Symbol)
    if sym === :core
        return getfield(m, :core)
    elseif sym === :groups
        return getfield(m, :groups)
    elseif hasfield(CoreGroup, sym)
        return getfield(getfield(m, :core), sym)
    end
    grps = getfield(m, :groups)
    for grp in values(grps)
        if hasfield(typeof(grp), sym)
            return getfield(grp, sym)
        end
    end
    error("type MarkerArrays has no field $(sym)")
end

function Base.propertynames(m::MarkerArrays, private::Bool=false)
    prop_list = Symbol[:core, :groups]
    append!(prop_list, fieldnames(CoreGroup))
    for grp in values(m.groups)
        append!(prop_list, fieldnames(typeof(grp)))
    end
    return tuple(prop_list...)
end

function Base.haskey(m::MarkerArrays, sym::Symbol)
    sym === :core && return true
    sym === :groups && return true
    hasfield(CoreGroup, sym) && return true
    for grp in values(m.groups)
        hasfield(typeof(grp), sym) && return true
    end
    return false
end

Base.haskey(m::MarkerArrays, key::AbstractString) = haskey(m, Symbol(key))
Base.getindex(m::MarkerArrays, sym::Symbol) = getproperty(m, sym)
Base.getindex(m::MarkerArrays, key::AbstractString) = getproperty(m, Symbol(key))

"""
    all_marker_array_names(m::MarkerArrays)

Return a vector of all array symbols present in `m.core` and active groups.
"""
function all_marker_array_names(m::MarkerArrays)
    names = Symbol[]
    for fn in fieldnames(CoreGroup)
        push!(names, fn)
    end
    for grp in values(m.groups)
        for fn in fieldnames(typeof(grp))
            push!(names, fn)
        end
    end
    return names
end

function Base.keys(m::MarkerArrays)
    return tuple(all_marker_array_names(m)...)
end

function Base.pairs(m::MarkerArrays)
    return (fn => getproperty(m, fn) for fn in all_marker_array_names(m))
end

function Base.iterate(m::MarkerArrays, state...)
    return iterate(pairs(m), state...)
end

"""
    Base.resize!(m::MarkerArrays, new_len::Integer)

Resize all vectors in core and active groups to `new_len`.
"""
function Base.resize!(m::MarkerArrays, new_len::Integer)
    for fn in fieldnames(CoreGroup)
        resize!(getfield(m.core, fn), new_len)
    end
    for grp in values(m.groups)
        for fn in fieldnames(typeof(grp))
            resize!(getfield(grp, fn), new_len)
        end
    end
    return m
end

"""
    Base.copy(m::MarkerArrays)

Return a deep copy of `MarkerArrays` with independent vector allocations.
"""
function Base.copy(m::MarkerArrays)
    core_copied = CoreGroup([copy(getfield(m.core, fn)) for fn in fieldnames(CoreGroup)]...)
    group_pairs = Pair{Symbol,Any}[]
    for (k, grp) in pairs(m.groups)
        T = typeof(grp)
        grp_copied = T([copy(getfield(grp, fn)) for fn in fieldnames(T)]...)
        push!(group_pairs, k => grp_copied)
    end
    groups_copied = NamedTuple(group_pairs)
    return MarkerArrays(core_copied, groups_copied)
end

"""
    restore_marker_arrays!(target::MarkerArrays, source::MarkerArrays)

Copy all vector contents from `source` to `target` in-place, resizing target arrays if needed.
"""
function restore_marker_arrays!(target::MarkerArrays, source::MarkerArrays)
    for fn in fieldnames(CoreGroup)
        tgt_vec = getfield(target.core, fn)
        src_vec = getfield(source.core, fn)
        if length(tgt_vec) != length(src_vec)
            resize!(tgt_vec, length(src_vec))
        end
        copyto!(tgt_vec, src_vec)
    end
    for (k, src_grp) in pairs(source.groups)
        if haskey(target.groups, k)
            tgt_grp = target.groups[k]
            T = typeof(src_grp)
            for fn in fieldnames(T)
                tgt_vec = getfield(tgt_grp, fn)
                src_vec = getfield(src_grp, fn)
                if length(tgt_vec) != length(src_vec)
                    resize!(tgt_vec, length(src_vec))
                end
                copyto!(tgt_vec, src_vec)
            end
        end
    end
    return target
end

"""
    push_marker!(m::MarkerArrays; kwargs...)

Append a single marker entry to core and present groups.
"""
function push_marker!(m::MarkerArrays; kwargs...)
    kw_dict = Dict{Symbol,Any}(kwargs)
    for fn in fieldnames(CoreGroup)
        vec = getfield(m.core, fn)
        val = get(kw_dict, fn, fn === :tm ? 1 : 0.0)
        push!(vec, val)
    end
    for grp in values(m.groups)
        for fn in fieldnames(typeof(grp))
            vec = getfield(grp, fn)
            val = get(kw_dict, fn, 0.0)
            push!(vec, val)
        end
    end
    return m
end

"""
    serialize_marker_arrays(m::MarkerArrays)

Export marker arrays to a dictionary with string keys for checkpointing.
"""
function serialize_marker_arrays(m::MarkerArrays)
    dict = Dict{String,Any}()
    for fn in fieldnames(CoreGroup)
        dict[string(fn)] = copy(getfield(m.core, fn))
    end
    for grp in values(m.groups)
        for fn in fieldnames(typeof(grp))
            dict[string(fn)] = copy(getfield(grp, fn))
        end
    end
    return dict
end

"""
    init_marker_arrays(marknum::Integer, cfg::SimulationConfig, coords::GridCoordinates; initial_time::Real=0.0, rng::AbstractRNG=Random.default_rng())

Construct a `MarkerArrays` container sizing core and configured active groups to `marknum`.
"""
function init_marker_arrays(
    marknum::Integer,
    cfg::SimulationConfig,
    coords::GridCoordinates;
    initial_time::Real=0.0,
    rng::AbstractRNG=Random.default_rng(),
)
    (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0, Fm) = setup_marker_properties(
        marknum, coords; rng=rng
    )
    (rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm, tkm_rhocptotalm, etafluidcur_inv_kphim, inv_gggtotalm, fricttotalm, cohestotalm, tenstotalm, rhofluidcur, alphasolidcur, alphafluidcur) = setup_marker_properties_helpers(
        marknum; rng=rng
    )

    w3d_m = Vector{Float64}(undef, marknum)
    for m in 1:marknum
        w3d_m[m] = marker_out_of_plane_length(xm[m], ym[m], coords.xcenter, coords.ycenter)
    end

    core = CoreGroup(
        xm,
        ym,
        w3d_m,
        tm,
        tkm,
        phim,
        phinewm,
        pfm0,
        XWsolidm,
        XWsolidm0,
        Fm,
        etavpm,
        sxxm,
        sxym,
        inv_gggtotalm,
        fricttotalm,
        cohestotalm,
        tenstotalm,
        rhototalm,
        rhocptotalm,
        etatotalm,
        hrtotalm,
        ktotalm,
        tkm_rhocptotalm,
        etafluidcur_inv_kphim,
        rhofluidcur,
        alphasolidcur,
        alphafluidcur,
    )

    group_pairs = Pair{Symbol,Any}[]

    # Metal group
    if cfg.coreformation.percolation_active ||
        cfg.coreformation.settling_active ||
        cfg.metal_partition.active ||
        cfg.thermodynamics.hr_fe
        Xfem, Xfem0, Xfe_bulk = setup_marker_metal_properties(marknum)
        Xfe_H_m, Xfe_C_m, Xfe_N_m, Xfe_S_m = if cfg.metal_partition.active
            setup_marker_metal_volatile_properties(
                marknum;
                initial_h_ppm=cfg.metal_partition.initial_metal_h_ppm,
                initial_c_ppm=cfg.metal_partition.initial_metal_c_ppm,
                initial_n_ppm=cfg.metal_partition.initial_metal_n_ppm,
                initial_s_ppm=cfg.metal_partition.initial_metal_s_ppm,
            )
        else
            (
                zeros(Float64, marknum),
                zeros(Float64, marknum),
                zeros(Float64, marknum),
                zeros(Float64, marknum),
            )
        end
        push!(
            group_pairs,
            :metal => MetalGroup(Xfem, Xfem0, Xfe_bulk, Xfe_H_m, Xfe_C_m, Xfe_N_m, Xfe_S_m),
        )
    end

    # Volatiles group
    if cfg.volatiles.active || cfg.magma_degassing.active || cfg.magma_transport.active
        XH2Om, XCm, XNm, XSm = if cfg.volatiles.active || cfg.magma_degassing.active
            setup_marker_volatile_properties(
                marknum;
                initial_water_wtpct=cfg.volatiles.initial_water_wtpct,
                initial_carbon_ppm=cfg.volatiles.initial_carbon_ppm,
                initial_nitrogen_ppm=cfg.volatiles.initial_nitrogen_ppm,
                initial_sulfur_ppm=cfg.volatiles.initial_sulfur_ppm,
            )
        else
            (
                zeros(Float64, marknum),
                zeros(Float64, marknum),
                zeros(Float64, marknum),
                zeros(Float64, marknum),
            )
        end
        X_graphite_m = zeros(Float64, marknum)
        F_extract_m = setup_marker_magma_properties(marknum)[1]
        push!(
            group_pairs,
            :volatiles => VolatilesGroup(XH2Om, XCm, XNm, XSm, X_graphite_m, F_extract_m),
        )
    end

    # Redox group
    if cfg.redox.active
        initial_xfe =
            haskey(Dict(group_pairs), :metal) ? Dict(group_pairs)[:metal].Xfe_bulk : nothing
        rp = setup_marker_redox_properties(
            marknum, cfg.redox; initial_xfe_bulk=initial_xfe, tkm=tkm, pfm=pfm0
        )
        push!(
            group_pairs,
            :redox => RedoxGroup(
                rp.nFe0_m,
                rp.nFe2_m,
                rp.nFe3_m,
                rp.deltaIW_m,
                rp.nC_graphite_m,
                rp.nCO_m,
                rp.nCO2_m,
                rp.nCH4_m,
            ),
        )
    end

    # HCN-S-P-O group
    if cfg.volatile_mixture.active || cfg.refractory.active
        hp = setup_marker_hcnspo_properties(marknum, cfg.volatile_mixture, cfg.refractory)
        push!(
            group_pairs,
            :hcnspo => HcnspoGroup(
                hp.X_ice_H2O_m,
                hp.X_ice_NH3_m,
                hp.X_ice_CO2_m,
                hp.X_ice_CO_m,
                hp.X_ice_CH4_m,
                hp.X_ice_N2_m,
                hp.X_ice_H2S_m,
                hp.X_ice_PH3_m,
                hp.X_refr_C_m,
                hp.X_refr_S_m,
                hp.X_refr_N_m,
                hp.X_refr_P_m,
                hp.X_refr_H_m,
            ),
        )
    end

    # Phase tracking group
    if cfg.phase_tracking.active
        phases = setup_marker_phase_tracking_properties(marknum, cfg.phase_tracking)
        push!(
            group_pairs,
            :phase => PhaseGroup(
                phases.Xmin_troilite_m,
                phases.Xmin_schreibersite_m,
                phases.Xmin_cohenite_m,
                phases.Xmin_graphite_m,
                phases.Xmin_nitride_m,
                phases.Xmin_metal_matrix_m,
            ),
        )
    end

    # Accretion group
    if cfg.accretion.active
        t_acc = setup_marker_accretion_properties(
            marknum, cfg.accretion; initial_time=initial_time
        )
        if t_acc !== nothing
            push!(group_pairs, :accretion => AccretionGroup(t_acc))
        end
    end

    groups = NamedTuple(group_pairs)
    return MarkerArrays(core, groups)
end

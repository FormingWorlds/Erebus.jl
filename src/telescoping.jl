"""
Telescoping domain algorithm for expanding spatial computational domains
during planetesimal accretion to lunar mass in Erebus.jl.

Maintains constant cell resolution dx across domain doubling events,
shifts existing markers to preserve physical radial distances from the
planetesimal center, seeds outer buffer cells with sticky air markers,
and triggers re-factorization of gravitational Poisson solvers.
"""

using DocStringExtensions
using LinearAlgebra
using SparseArrays

"""
    should_telescope_domain(rplanet::Real, coords::GridCoordinates, cfg::TelescopingConfig; level::Integer=0)::Bool

Evaluates whether the planetesimal radius has exceeded the domain threshold fraction
to trigger a telescoping domain doubling event.

# Arguments
- `rplanet`: Current planetesimal radius [m].
- `coords`: Current spatial grid coordinates.
- `cfg`: Telescoping domain configuration.
- `level`: Current telescoping expansion level (0-indexed).

# Returns
- `true` if accretion has grown the body beyond `r_threshold_fraction * (xsize / 2)` and
  maximum telescoping levels have not been reached.
"""
function should_telescope_domain(
    rplanet::Real, coords::GridCoordinates, cfg::TelescopingConfig; level::Integer=0
)::Bool
    if !isfinite(rplanet) || rplanet < 0.0
        throw(DomainError(rplanet, "rplanet must be non-negative and finite"))
    end
    if !cfg.active || level >= cfg.max_telescope_levels
        return false
    end
    half_domain = coords.xsize / 2.0
    r_threshold = cfg.r_threshold_fraction * half_domain
    return Float64(rplanet) > r_threshold
end

function should_telescope_domain(
    rplanet::Real, coords::GridCoordinates, cfg::SimulationConfig; level::Integer=0
)::Bool
    return should_telescope_domain(rplanet, coords, cfg.telescoping; level=level)
end

"""
    compute_telescoped_coordinates(coords::GridCoordinates)::GridCoordinates

Constructs a doubled spatial grid domain while preserving constant cell spacing `dx` and `dy`.
New node counts satisfy `Nx_new = 2 * (Nx - 1) + 1`.

$(SIGNATURES)
"""
function compute_telescoped_coordinates(coords::GridCoordinates)::GridCoordinates
    if iseven(coords.Nx) || iseven(coords.Ny)
        throw(
            ArgumentError(
                "Telescoping domain requires odd Nx and Ny for staggered grid centering (got Nx=$(coords.Nx), Ny=$(coords.Ny))",
            ),
        )
    end
    Nx_new = 2 * (coords.Nx - 1) + 1
    Ny_new = 2 * (coords.Ny - 1) + 1
    xsize_new = 2.0 * coords.xsize
    ysize_new = 2.0 * coords.ysize
    return GridCoordinates(
        Nx_new, Ny_new; xsize=xsize_new, ysize=ysize_new, Nxmc=coords.Nxmc, Nymc=coords.Nymc
    )
end

"""
    remap_staggered_grid_array(old_arr::AbstractMatrix{T}, new_dims::Tuple{Int,Int}; background_val::Real=zero(T)) where {T}

Remaps a 2D staggered grid array into a larger telescoped grid dimension by centering
the original array and setting outer buffer cells to `background_val`.

$(SIGNATURES)
"""
function remap_staggered_grid_array(
    old_arr::AbstractMatrix{T}, new_dims::Tuple{Int,Int}; background_val::Real=zero(T)
)::Matrix{T} where {T}
    Ny_old, Nx_old = size(old_arr)
    Ny_new, Nx_new = new_dims
    if Ny_new < Ny_old || Nx_new < Nx_old
        throw(
            ArgumentError(
                "Target grid dimensions $new_dims must be >= original dimensions $(size(old_arr))",
            ),
        )
    end
    if isodd(Ny_new - Ny_old) || isodd(Nx_new - Nx_old)
        throw(
            ArgumentError(
                "Dimension difference ($Ny_new - $Ny_old, $Nx_new - $Nx_old) must be even for symmetric centering",
            ),
        )
    end
    ioff = (Ny_new - Ny_old) ÷ 2
    joff = (Nx_new - Nx_old) ÷ 2

    new_arr = fill(T(background_val), Ny_new, Nx_new)
    @views new_arr[(ioff + 1):(ioff + Ny_old), (joff + 1):(joff + Nx_old)] .= old_arr
    return new_arr
end

"""
    push_redox_marker!(redox_props, deltaIW_ambient::Float64)

Append sticky-air marker defaults to all active arrays in `redox_props`.

# Parameters
- `redox_props`: Redox properties container.
- `deltaIW_ambient`: Ambient delta IW oxygen fugacity buffer value.
"""
function push_redox_marker!(redox_props, deltaIW_ambient::Float64)
    redox_props === nothing && return nothing
    for fn in fieldnames(RedoxGroup)
        if hasproperty(redox_props, fn)
            arr = getproperty(redox_props, fn)
            if arr !== nothing
                push!(arr, fn === :deltaIW_m ? deltaIW_ambient : 0.0)
            end
        end
    end
    return nothing
end

@inline push_if_not_nothing!(arr, val=0.0) = (arr !== nothing && push!(arr, val); nothing)

"""
    telescope_marker_arrays!(
        xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0,
        XWsolidm, XWsolidm0, Fm, rhototalm, rhocptotalm, etatotalm,
        hrtotalm, ktotalm, inv_gggtotalm, fricttotalm, cohestotalm,
        tenstotalm, rhofluidcur, alphasolidcur, alphafluidcur,
        tkm_rhocptotalm, etafluidcur_inv_kphim;
        old_coords::GridCoordinates,
        new_coords::GridCoordinates,
        T_ambient::Real=250.0,
        phi_ambient::Real=0.35,
        buffer_markers_per_cell::Integer=4,
        Xfem=nothing,
        Xfem0=nothing,
        Xfe_bulk=nothing,
        XH2Om=nothing,
        XCm=nothing,
        XNm=nothing,
        XSm=nothing,
        Xfe_H_m=nothing,
        Xfe_C_m=nothing,
        Xfe_N_m=nothing,
        Xfe_S_m=nothing,
        Xmin_troilite_m=nothing,
        Xmin_schreibersite_m=nothing,
        Xmin_cohenite_m=nothing,
        Xmin_graphite_m=nothing,
        Xmin_nitride_m=nothing,
        Xmin_metal_matrix_m=nothing,
        t_accreted=nothing,
    )::Int

Translates existing markers by the domain offset `(shift_x, shift_y)` to preserve
exact physical radial distance from the planetesimal center, then populates newly
created outer buffer cells with ambient sticky air markers (`tm = 3`).

# Returns
- `new_marknum`: Total number of active markers following buffer replenishment.
"""
function telescope_marker_arrays!(
    xm::AbstractVector{<:Real},
    ym::AbstractVector{<:Real},
    tm::AbstractVector{<:Integer},
    tkm::AbstractVector{<:Real},
    sxxm::AbstractVector{<:Real},
    sxym::AbstractVector{<:Real},
    etavpm::AbstractVector{<:Real},
    phim::AbstractVector{<:Real},
    phinewm::AbstractVector{<:Real},
    pfm0::AbstractVector{<:Real},
    XWsolidm::AbstractVector{<:Real},
    XWsolidm0::AbstractVector{<:Real},
    Fm::AbstractVector{<:Real},
    rhototalm::AbstractVector{<:Real},
    rhocptotalm::AbstractVector{<:Real},
    etatotalm::AbstractVector{<:Real},
    hrtotalm::AbstractVector{<:Real},
    ktotalm::AbstractVector{<:Real},
    inv_gggtotalm::AbstractVector{<:Real},
    fricttotalm::AbstractVector{<:Real},
    cohestotalm::AbstractVector{<:Real},
    tenstotalm::AbstractVector{<:Real},
    rhofluidcur::AbstractVector{<:Real},
    alphasolidcur::AbstractVector{<:Real},
    alphafluidcur::AbstractVector{<:Real},
    tkm_rhocptotalm::AbstractVector{<:Real},
    etafluidcur_inv_kphim::AbstractVector{<:Real};
    old_coords::GridCoordinates,
    new_coords::GridCoordinates,
    T_ambient::Real=250.0,
    phi_ambient::Real=0.35,
    buffer_markers_per_cell::Integer=4,
    Xfem::Union{Nothing,AbstractVector{<:Real}}=nothing,
    Xfem0::Union{Nothing,AbstractVector{<:Real}}=nothing,
    Xfe_bulk::Union{Nothing,AbstractVector{<:Real}}=nothing,
    XH2Om::Union{Nothing,AbstractVector{<:Real}}=nothing,
    XCm::Union{Nothing,AbstractVector{<:Real}}=nothing,
    XNm::Union{Nothing,AbstractVector{<:Real}}=nothing,
    XSm::Union{Nothing,AbstractVector{<:Real}}=nothing,
    Xfe_H_m::Union{Nothing,AbstractVector{<:Real}}=nothing,
    Xfe_C_m::Union{Nothing,AbstractVector{<:Real}}=nothing,
    Xfe_N_m::Union{Nothing,AbstractVector{<:Real}}=nothing,
    Xfe_S_m::Union{Nothing,AbstractVector{<:Real}}=nothing,
    Xmin_troilite_m::Union{Nothing,AbstractVector{<:Real}}=nothing,
    Xmin_schreibersite_m::Union{Nothing,AbstractVector{<:Real}}=nothing,
    Xmin_cohenite_m::Union{Nothing,AbstractVector{<:Real}}=nothing,
    Xmin_graphite_m::Union{Nothing,AbstractVector{<:Real}}=nothing,
    Xmin_nitride_m::Union{Nothing,AbstractVector{<:Real}}=nothing,
    Xmin_metal_matrix_m::Union{Nothing,AbstractVector{<:Real}}=nothing,
    t_accreted::Union{Nothing,AbstractVector{<:Real}}=nothing,
    materials::Union{Nothing,MaterialConfig}=nothing,
    hcnspo_props=nothing,
    F_extract_m::Union{Nothing,AbstractVector{<:Real}}=nothing,
    w3d_m::Union{Nothing,AbstractVector{<:Real}}=nothing,
    X_graphite_m::Union{Nothing,AbstractVector{<:Real}}=nothing,
    redox_props=nothing,
    deltaIW_ambient::Real=0.0,
    cfg::Union{Nothing,SimulationConfig}=nothing,
)::Int
    if iseven(old_coords.Nx) || iseven(old_coords.Ny)
        throw(
            ArgumentError(
                "Telescoping domain requires odd Nx and Ny for staggered grid centering (got Nx=$(old_coords.Nx), Ny=$(old_coords.Ny))",
            ),
        )
    end
    if iseven(new_coords.Nx) || iseven(new_coords.Ny)
        throw(
            ArgumentError(
                "Telescoping domain requires odd Nx and Ny for new staggered grid centering (got Nx=$(new_coords.Nx), Ny=$(new_coords.Ny))",
            ),
        )
    end
    if length(ym) != length(xm)
        throw(
            DimensionMismatch(
                "Coordinate vector lengths do not match: length(xm)=$(length(xm)), length(ym)=$(length(ym))",
            ),
        )
    end
    shift_x = new_coords.xcenter - old_coords.xcenter
    shift_y = new_coords.ycenter - old_coords.ycenter

    N_old = length(xm)
    for m in 1:N_old
        @inbounds xm[m] += shift_x
        @inbounds ym[m] += shift_y
    end

    Nx_cells_new = new_coords.Nx - 1
    Ny_cells_new = new_coords.Ny - 1
    Nx_cells_old = old_coords.Nx - 1
    Ny_cells_old = old_coords.Ny - 1
    joff_cells = (Nx_cells_new - Nx_cells_old) ÷ 2
    ioff_cells = (Ny_cells_new - Ny_cells_old) ÷ 2

    dx = new_coords.dx
    dy = new_coords.dy

    nx_sub = isqrt(buffer_markers_per_cell)
    ny_sub = nx_sub
    if nx_sub * ny_sub != buffer_markers_per_cell
        for d in nx_sub:-1:1
            if buffer_markers_per_cell % d == 0
                nx_sub = d
                ny_sub = buffer_markers_per_cell ÷ d
                break
            end
        end
    end
    @assert nx_sub * ny_sub == buffer_markers_per_cell
    dx_sub = dx / nx_sub
    dy_sub = dy / ny_sub

    T_amb = Float64(T_ambient)
    phi_amb = Float64(phi_ambient)

    # Use materials configuration for sticky-air phase (index 3) when provided
    rho_air = materials !== nothing ? materials.rhosolidm[3] : 1.0
    rhofluid_air = materials !== nothing ? materials.rhofluidm[3] : 1.0
    eta_air = materials !== nothing ? materials.etasolidm[3] : 1.0e16
    rhocp_air = materials !== nothing ? materials.rhocpsolidm[3] : 3.0e6
    k_air = materials !== nothing ? materials.ksolidm[3] : 3000.0
    inv_ggg_air = materials !== nothing ? inv(materials.gggsolidm[3]) : 1.0e-10
    frict_air = materials !== nothing ? materials.frictsolidm[3] : 0.0
    cohes_air = materials !== nothing ? materials.cohessolidm[3] : 1.0e8
    tens_air = materials !== nothing ? materials.tenssolidm[3] : 6.0e7
    alphasolid_air = materials !== nothing ? materials.alphasolidm[3] : 0.0
    alphafluid_air = materials !== nothing ? materials.alphafluidm[3] : 0.0
    tkm_rhocp_air = T_amb * rhocp_air

    for j in 1:Nx_cells_new
        for i in 1:Ny_cells_new
            is_inner =
                (joff_cells < j <= joff_cells + Nx_cells_old) &&
                (ioff_cells < i <= ioff_cells + Ny_cells_old)
            if !is_inner
                cell_x0 = (j - 1) * dx
                cell_y0 = (i - 1) * dy
                for iy in 1:ny_sub
                    for ix in 1:nx_sub
                        x_marker = cell_x0 + (ix - 0.5) * dx_sub
                        y_marker = cell_y0 + (iy - 0.5) * dy_sub

                        push!(xm, x_marker)
                        push!(ym, y_marker)
                        push!(tm, 3)
                        push!(tkm, T_amb)
                        push!(sxxm, 0.0)
                        push!(sxym, 0.0)
                        push!(etavpm, eta_air)
                        push!(phim, phi_amb)
                        push!(phinewm, phi_amb)
                        push!(pfm0, 0.0)
                        push!(XWsolidm, 0.0)
                        push!(XWsolidm0, 0.0)
                        push!(Fm, 0.0)
                        push!(rhototalm, rho_air)
                        push!(rhocptotalm, rhocp_air)
                        push!(etatotalm, eta_air)
                        push!(hrtotalm, 0.0)
                        push!(ktotalm, k_air)
                        push!(inv_gggtotalm, inv_ggg_air)
                        push!(fricttotalm, frict_air)
                        push!(cohestotalm, cohes_air)
                        push!(tenstotalm, tens_air)
                        push!(rhofluidcur, rhofluid_air)
                        push!(alphasolidcur, alphasolid_air)
                        push!(alphafluidcur, alphafluid_air)
                        push!(tkm_rhocptotalm, tkm_rhocp_air)
                        push!(etafluidcur_inv_kphim, 1.0e14)

                        push_if_not_nothing!(Xfem)
                        push_if_not_nothing!(Xfem0)
                        push_if_not_nothing!(Xfe_bulk)
                        push_if_not_nothing!(XH2Om)
                        push_if_not_nothing!(XCm)
                        push_if_not_nothing!(XNm)
                        push_if_not_nothing!(XSm)
                        push_if_not_nothing!(Xfe_H_m)
                        push_if_not_nothing!(Xfe_C_m)
                        push_if_not_nothing!(Xfe_N_m)
                        push_if_not_nothing!(Xfe_S_m)
                        push_if_not_nothing!(Xmin_troilite_m)
                        push_if_not_nothing!(Xmin_schreibersite_m)
                        push_if_not_nothing!(Xmin_cohenite_m)
                        push_if_not_nothing!(Xmin_graphite_m)
                        push_if_not_nothing!(Xmin_nitride_m)
                        push_if_not_nothing!(Xmin_metal_matrix_m)
                        push_if_not_nothing!(t_accreted)
                        push_if_not_nothing!(F_extract_m)
                        if w3d_m !== nothing
                            push!(
                                w3d_m,
                                marker_out_of_plane_length(
                                    x_marker,
                                    y_marker,
                                    new_coords.xcenter,
                                    new_coords.ycenter,
                                ),
                            )
                        end
                        push_if_not_nothing!(X_graphite_m)
                        push_redox_marker!(redox_props, Float64(deltaIW_ambient))
                        if hcnspo_props !== nothing
                            for prop in values(hcnspo_props)
                                push_if_not_nothing!(prop)
                            end
                        end
                    end
                end
            end
        end
    end

    return length(xm)
end

"""
    telescope_marker_arrays!(
        markers::MarkerArrays;
        old_coords::GridCoordinates,
        new_coords::GridCoordinates,
        cfg::Union{Nothing,SimulationConfig}=nothing,
        buffer_markers_per_cell::Integer=4,
    )::Int

Extend all active arrays in `markers.core` and `markers.groups` following domain telescoping.
Ensures array length invariants are maintained across all groups.
"""
function telescope_marker_arrays!(
    markers::MarkerArrays;
    old_coords::GridCoordinates,
    new_coords::GridCoordinates,
    cfg::Union{Nothing,SimulationConfig}=nothing,
    buffer_markers_per_cell::Integer=4,
)::Int
    c = markers.core
    grps = markers.groups
    mats = cfg !== nothing ? cfg.materials : nothing

    Xfem = haskey(grps, :metal) ? grps[:metal].Xfem : nothing
    Xfem0 = haskey(grps, :metal) ? grps[:metal].Xfem0 : nothing
    Xfe_bulk = haskey(grps, :metal) ? grps[:metal].Xfe_bulk : nothing
    Xfe_H_m = haskey(grps, :metal) ? grps[:metal].Xfe_H_m : nothing
    Xfe_C_m = haskey(grps, :metal) ? grps[:metal].Xfe_C_m : nothing
    Xfe_N_m = haskey(grps, :metal) ? grps[:metal].Xfe_N_m : nothing
    Xfe_S_m = haskey(grps, :metal) ? grps[:metal].Xfe_S_m : nothing

    XH2Om = haskey(grps, :volatiles) ? grps[:volatiles].XH2Om : nothing
    XCm = haskey(grps, :volatiles) ? grps[:volatiles].XCm : nothing
    XNm = haskey(grps, :volatiles) ? grps[:volatiles].XNm : nothing
    XSm = haskey(grps, :volatiles) ? grps[:volatiles].XSm : nothing
    X_graphite_m = haskey(grps, :volatiles) ? grps[:volatiles].X_graphite_m : nothing
    F_extract_m = haskey(grps, :volatiles) ? grps[:volatiles].F_extract_m : nothing

    redox_props = haskey(grps, :redox) ? grps[:redox] : nothing
    hcnspo_props = haskey(grps, :hcnspo) ? grps[:hcnspo] : nothing

    Xmin_troilite_m = haskey(grps, :phase) ? grps[:phase].Xmin_troilite_m : nothing
    Xmin_schreibersite_m =
        haskey(grps, :phase) ? grps[:phase].Xmin_schreibersite_m : nothing
    Xmin_cohenite_m = haskey(grps, :phase) ? grps[:phase].Xmin_cohenite_m : nothing
    Xmin_graphite_m = haskey(grps, :phase) ? grps[:phase].Xmin_graphite_m : nothing
    Xmin_nitride_m = haskey(grps, :phase) ? grps[:phase].Xmin_nitride_m : nothing
    Xmin_metal_matrix_m = haskey(grps, :phase) ? grps[:phase].Xmin_metal_matrix_m : nothing

    t_accreted = haskey(grps, :accretion) ? grps[:accretion].t_accreted : nothing

    d_IW_amb = if cfg !== nothing && hasproperty(cfg, :volatiles)
        cfg.volatiles.fO2_delta_IW
    else
        0.0
    end

    T_amb =
        if cfg !== nothing &&
            hasproperty(cfg, :materials) &&
            length(cfg.materials.tkm0) >= 3
            cfg.materials.tkm0[3]
        else
            250.0
        end

    phi_amb = if cfg !== nothing && hasproperty(cfg, :poroelasticity)
        cfg.poroelasticity.phimin
    else
        0.35
    end

    new_marknum = telescope_marker_arrays!(
        c.xm,
        c.ym,
        c.tm,
        c.tkm,
        c.sxxm,
        c.sxym,
        c.etavpm,
        c.phim,
        c.phinewm,
        c.pfm0,
        c.XWsolidm,
        c.XWsolidm0,
        c.Fm,
        c.rhototalm,
        c.rhocptotalm,
        c.etatotalm,
        c.hrtotalm,
        c.ktotalm,
        c.inv_gggtotalm,
        c.fricttotalm,
        c.cohestotalm,
        c.tenstotalm,
        c.rhofluidcur,
        c.alphasolidcur,
        c.alphafluidcur,
        c.tkm_rhocptotalm,
        c.etafluidcur_inv_kphim;
        old_coords=old_coords,
        new_coords=new_coords,
        T_ambient=T_amb,
        phi_ambient=phi_amb,
        buffer_markers_per_cell=buffer_markers_per_cell,
        materials=mats,
        w3d_m=c.w3d_m,
        Xfem=Xfem,
        Xfem0=Xfem0,
        Xfe_bulk=Xfe_bulk,
        XH2Om=XH2Om,
        XCm=XCm,
        XNm=XNm,
        XSm=XSm,
        Xfe_H_m=Xfe_H_m,
        Xfe_C_m=Xfe_C_m,
        Xfe_N_m=Xfe_N_m,
        Xfe_S_m=Xfe_S_m,
        Xmin_troilite_m=Xmin_troilite_m,
        Xmin_schreibersite_m=Xmin_schreibersite_m,
        Xmin_cohenite_m=Xmin_cohenite_m,
        Xmin_graphite_m=Xmin_graphite_m,
        Xmin_nitride_m=Xmin_nitride_m,
        Xmin_metal_matrix_m=Xmin_metal_matrix_m,
        t_accreted=t_accreted,
        hcnspo_props=hcnspo_props,
        F_extract_m=F_extract_m,
        X_graphite_m=X_graphite_m,
        redox_props=redox_props,
        deltaIW_ambient=d_IW_amb,
        cfg=cfg,
    )
    assert_marker_arrays_invariants(markers, new_marknum)
    return new_marknum
end

"""
    assert_marker_arrays_invariants(markers::MarkerArrays, marknum::Integer)

Assert that all marker vectors in `markers.core` and all active groups in `markers.groups`
have length equal to `marknum`. Throw `DimensionMismatch` if any vector length deviates.
"""
function assert_marker_arrays_invariants(markers::MarkerArrays, marknum::Integer)
    for fn in fieldnames(CoreGroup)
        arr = getfield(markers.core, fn)
        if length(arr) != marknum
            throw(
                DimensionMismatch(
                    "Marker core array :$fn has length $(length(arr)), expected marknum=$marknum",
                ),
            )
        end
    end
    for (gname, grp) in pairs(markers.groups)
        for fn in fieldnames(typeof(grp))
            arr = getfield(grp, fn)
            if length(arr) != marknum
                throw(
                    DimensionMismatch(
                        "Marker group :$gname array :$fn has length $(length(arr)), expected marknum=$marknum",
                    ),
                )
            end
        end
    end
    return true
end

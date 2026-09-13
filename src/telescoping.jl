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
)::Int
    if iseven(old_coords.Nx) || iseven(old_coords.Ny)
        throw(
            ArgumentError(
                "Telescoping domain requires odd Nx and Ny for staggered grid centering (got Nx=$(old_coords.Nx), Ny=$(old_coords.Ny))",
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

                        if Xfem !== nothing
                            push!(Xfem, 0.0)
                        end
                        if Xfem0 !== nothing
                            push!(Xfem0, 0.0)
                        end
                        if Xfe_bulk !== nothing
                            push!(Xfe_bulk, 0.0)
                        end
                        if XH2Om !== nothing
                            push!(XH2Om, 0.0)
                        end
                        if XCm !== nothing
                            push!(XCm, 0.0)
                        end
                        if XNm !== nothing
                            push!(XNm, 0.0)
                        end
                        if XSm !== nothing
                            push!(XSm, 0.0)
                        end
                        if Xfe_H_m !== nothing
                            push!(Xfe_H_m, 0.0)
                        end
                        if Xfe_C_m !== nothing
                            push!(Xfe_C_m, 0.0)
                        end
                        if Xfe_N_m !== nothing
                            push!(Xfe_N_m, 0.0)
                        end
                        if Xfe_S_m !== nothing
                            push!(Xfe_S_m, 0.0)
                        end
                        if Xmin_troilite_m !== nothing
                            push!(Xmin_troilite_m, 0.0)
                        end
                        if Xmin_schreibersite_m !== nothing
                            push!(Xmin_schreibersite_m, 0.0)
                        end
                        if Xmin_cohenite_m !== nothing
                            push!(Xmin_cohenite_m, 0.0)
                        end
                        if Xmin_graphite_m !== nothing
                            push!(Xmin_graphite_m, 0.0)
                        end
                        if Xmin_nitride_m !== nothing
                            push!(Xmin_nitride_m, 0.0)
                        end
                        if Xmin_metal_matrix_m !== nothing
                            push!(Xmin_metal_matrix_m, 0.0)
                        end
                        if t_accreted !== nothing
                            push!(t_accreted, 0.0)
                        end
                        if F_extract_m !== nothing
                            push!(F_extract_m, 0.0)
                        end
                        if hcnspo_props !== nothing
                            for prop in values(hcnspo_props)
                                push!(prop, 0.0)
                            end
                        end
                    end
                end
            end
        end
    end

    return length(xm)
end

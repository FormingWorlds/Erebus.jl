using ArgParse
using DocStringExtensions
using JLD2
using LinearAlgebra
using Printf
using ProgressMeter
using Random 
using StaticArrays
using Statistics
using StatsBase
import PyPlot; const plt = PyPlot


macro gprintf(fmt::String)
    :((io::IO, arg) -> Printf.@printf(io, $fmt, arg))
end

const fi = @gprintf "%d"
const f2p = @gprintf "%.2f"
const f3p = @gprintf "%.3f"
const f6p = @gprintf "%.6f"
# jit_printf(fmt) = @eval @gprintf($fmt)

const seed = 42
const rgen = MersenneTwister(seed)

"""
Get x- and y-ticks in scientific notation for given x- and y- ranges.

$(SIGNATURES)

# Details

    - x: x-range
    - y: y-range
    - n: number of ticks

# Returns

    - xticks: (ticks, ticklabels) for x-axis
    - yticks: (ticks, ticklabels) for y-axis
"""
function label_ticks(x, y; n=10)
    xticks = range(x[1], x[end], n)
    yticks = range(y[1], y[end], n)
    xlabels = map(t -> @sprintf("%.1E", t), xticks)
    ylabels = map(t -> @sprintf("%.1E", t), yticks)
    return (xticks, xlabels), (yticks, ylabels)
end    

"""
Get x- and y- vectors as an analogue to the 2D function `meshgrid()` as known
from another famous language.

$(SIGNATURES)

# Details

    - x: x-range vector
    - y: y-range vector

# Returns

    - (X, Y): a tuple of vectors describing the 2D meshgrid
"""
function meshgrid(x, y)
    return repeat(x, outer=length(y)), repeat(y, inner=length(x))
end

"""
Obtain extrema and labeled colorbar ticks for given scalar field.

$(SIGNATURES)

# Details

    - grid: 2D array containing the scalar field
    - n_ticks: number of ticks on the colorbar

# Returns

    - (limits, (cbar_ticks, cbar_lables)): tuple of extrema and colorbar ticks
"""
function colorbar_properties(grid, n_ticks=5)
    limits = get_extrema(grid)
    cbartickrange = range(limits..., n_ticks)
    cbarticklabels = map(t -> @sprintf("%.1E", t), cbartickrange)
    return limits, (cbartickrange, cbarticklabels)
end

"""
Obtain an array's extrema, but return (0, 1) in case the extrema range is
invalid for plotting.

$(SIGNATURES)

# Details

    - grid: 2D array containing the scalar field whose extrema are sought

# Returns

    - (min, max): tuple of extrema
"""
function get_extrema(grid)
    a, b = extrema(grid)
    if a==b ||  (a==-Inf && b==Inf)
        return zero(0.0), one(1.0)
    else
        return a, b
    end
end

"""
Compute longitudinal average of input field for given radius.

$(SIGNATURES)

# Details

    - f: input field
    - r: radius of the radial average [matrix cell units]
    - xcenter: x-coordinate of the center of the radial average [matrix
        cell units]
    - ycenter: y-coordinate of the center of the radial average [matrix
        cell units]

# Returns

    - ravg: radial average of the field at radius r
"""
function ravg(f, r)
    ysize, xsize = size(f)
    ycenter = trunc(Int, ysize / 2)
    xcenter = trunc(Int, xsize / 2)
    r = min(r, ycenter, xcenter, ysize-ycenter, xsize-xcenter)
    circle = zeros(Bool, ysize, xsize)
    for θ = 0:0.01:2π
        @inbounds circle[
            round(Int, ycenter + r*cos(θ)),
            round(Int, xcenter + r*sin(θ))
        ] = true
        @inbounds return mean(f[circle])
    end
end

"""
Compute radial velocity field from center of given input x- and y-velocity
fields.

$(SIGNATURES)

# Details

    - vx: x-velocity field
    - vy: y-velocity field

# Returns

    -vrad: radial velocity field
"""
function vrad(vx, vy)
    @assert size(vx) == size(vy)
    ysize, xsize = size(vx)
    vrad = zero(vx)
    ycenter = trunc(Int, ysize / 2)
    xcenter = trunc(Int, xsize / 2)
    r_max = max(round(Int, xsize/2), round(Int, ysize/2)) * sqrt(2.0)
    for ϕ = 0.0:0.01:2π, r = 1:1:r_max
        y = round(Int, ycenter + r*cos(ϕ))
        x = round(Int, xcenter + r*sin(ϕ))
        if 1<=y<=ysize && 1<=x<=xsize
            # vrad[y, x] = sqrt((vx[y, x]*sin(ϕ))^2 + (vy[y, x]*cos(ϕ))^2)
            dx = x - xcenter
            dy = y - ycenter
            vrad[y, x] = (dx*vx[y, x] + dy*vy[y, x]) / sqrt(dx^2 + dy^2)
        end
    end
    return vrad
end

"""
Create a 2D heatmap plot of a given scalar field.

$(SIGNATURES)

# Details

    - ax: PyPlot axis to plot on
    - field: 2D array containing the scalar field to be plotted
    - vmin: minimum value of the colorbar
    - vmax: maximum value of the colorbar
    - t: time step of the field [Ma]
    - cmap: color map to use
    - fontsize: font size of the labels
    - hide_labels: switch to hide heatmap axis labels

# Returns

    - im: PyPlot image object
"""
function plot_field(
    ax, field, vmin, vmax, t, cmap, fontsize=10, hide_labels=false)
    # im = ax.contourf(field, cmap=cmap, aspect="equal")
    im = ax.imshow(field, vmin=vmin, vmax=vmax, cmap=cmap, aspect="equal")
    # ax.set_aspect("equal")
    ax.locator_params(nbins=3)
    if hide_labels
        ax.set_xticklabels([])
        ax.set_yticklabels([])
    else
        ax.set_xlabel("x [km]", fontsize=fontsize)
        ax.set_ylabel("y [km]", fontsize=fontsize)
        ax.set_title("t = $(sprint(f3p, t)) Ma", fontsize=fontsize)
    end
    return im
end

"""
Get azimuth (polar coordinate angle) and radius of a set of markers for given 
center coordinate.

$(SIGNATURES)

# Details

    - xm: x-coordinates of the markers [m]
    - ym: y-coordinates of the markers [m]
    - xcenter: x-coordinate of the center of the radial average [m]
    - ycenter: y-coordinate of the center of the radial average [m]

# Returns

    - marker_azimuths: azimuths of the markers [deg]
    - marker_radii: azimuths of the markers [deg]
"""
function get_marker_azimuths_radii(xm, ym, xcenter, ycenter)
    atan2(y, x) = ifelse(-atan(y,x)+5π/2 >=2π, -atan(y,x)+π/2, -atan(y,x)+5π/2)
    Δy = ym .- ycenter
    Δx = xm .- xcenter
    return (atan2.(Δy, Δx) .* 180.0 ./ π,  sqrt.(Δx.^2 + Δy.^2)) 
end

"""
Find indices of markers within a given radius of a given center coordinate.

$(SIGNATURES)

# Details

    - xm: x-coordinates of the markers [m]
    - ym: y-coordinates of the markers [m]
    - rmin: minimum radius of the search [m]
    - rmax: maximum radius of the search [m]
    - xcenter: x-coordinate of the center of the radial average [m]
    - ycenter: y-coordinate of the center of the radial average [m]

# Returns

    - selection: indices of the markers within the search radius
"""
function get_markers_in_radius(xm, ym, rmin, rmax, xcenter, ycenter)
    selection = collect(1:size(xm, 1))
    R = get_marker_radii(xm, ym, selection, xcenter, ycenter)
    return findall(r -> rmin<=r<=rmax, R)
end

"""
Create set of plots from results stored in `.jld2` simulation output files.

$(SIGNATURES)

# Details

    - input_path: path to the directory containing the simulation output files

# Returns

    - nothing
"""
function generate_plots(input_path)
    files = filter!(
        x -> isfile(x) && endswith(x, ".jld2"), readdir(input_path, join=true))
    n_steps = length(files)
    file = jldopen(files[end], "r")
    Nx = file["Nx"]
    Ny = file["Ny"]
    Nx1 = file["Nx1"]
    Ny1 = file["Ny1"]
    x = file["x"]
    y = file["y"]
    xvx = file["xvx"]
    yvx = file["yvx"]
    xvy = file["xvy"]
    yvy = file["yvy"]
    xp = file["xp"]
    yp = file["yp"]
    xxm = file["xxm"]
    yym = file["yym"]
    rplanet = file["rplanet"]
    xcenter = file["xcenter"]
    ycenter = file["ycenter"]
    dx = file["dx"]
    dy = file["dy"]
    end_marknum = file["marknum"]
    phim0 = file["phim0"]
    ratio_al = file["ratio_al"]
    XWsolidm_init = file["XWsolidm_init"]
    close(file)
    radius_max = trunc(Int, min(rplanet/dx, rplanet/dy))
    radius_range = collect(1:1:radius_max)
    n_ticks = 4

    xlim_b, ylim_b = extrema.((x, y))
    xlim_vx, ylim_vx = extrema.((xvx, yvx))
    xlim_vy, ylim_vy = extrema.((xvy, yvy))
    xlim_p, ylim_p = extrema.((xp, yp))
    xlim_m, ylim_m = extrema.((xxm, yym))
    xticks_b, yticks_b = label_ticks(xlim_b, ylim_b, n=n_ticks)
    xticks_vx, yticks_vx = label_ticks(xlim_vx, ylim_vx, n=n_ticks)
    xticks_vy, yticks_vy = label_ticks(xlim_vy, ylim_vy, n=n_ticks)
    xticks_p, yticks_p = label_ticks(xp, yp, n=n_ticks)
    xticks_m, yticks_m = label_ticks(xxm, yym, n=n_ticks)

    RHO = Array{Float64}(undef, Ny1, Nx1, n_steps)
    ETA = Array{Float64}(undef, Ny, Nx, n_steps)
    vx = Array{Float64}(undef, Ny1, Nx1, n_steps)
    vy = Array{Float64}(undef, Ny1, Nx1, n_steps)
    vxp = Array{Float64}(undef, Ny1, Nx1, n_steps)
    vyp = Array{Float64}(undef, Ny1, Nx1, n_steps)
    PHI = Array{Float64}(undef, Ny1, Nx1, n_steps)
    EII = Array{Float64}(undef, Ny1, Nx1, n_steps)
    SII = Array{Float64}(undef, Ny1, Nx1, n_steps)
    gx = Array{Float64}(undef, Ny1, Nx1, n_steps)
    gy = Array{Float64}(undef, Ny1, Nx1, n_steps)
    tk2 = Array{Float64}(undef, Ny1, Nx1, n_steps)
    HS = Array{Float64}(undef, Ny1, Nx1, n_steps)
    HA = Array{Float64}(undef, Ny1, Nx1, n_steps)
    pr = Array{Float64}(undef, Ny1, Nx1, n_steps)
    pf = Array{Float64}(undef, Ny1, Nx1, n_steps)
    KX = Array{Float64}(undef, Ny1, Nx1, n_steps)
    qxD = Array{Float64}(undef, Ny1, Nx1, n_steps)
    qyD = Array{Float64}(undef, Ny1, Nx1, n_steps)
    RX = Array{Float64}(undef, Ny1, Nx1, n_steps)
    ETAPHI = Array{Float64}(undef, Ny1, Nx1, n_steps)
    DMP = Array{Float64}(undef, Ny1, Nx1, n_steps)
    DHP = Array{Float64}(undef, Ny1, Nx1, n_steps)
    XWS = Array{Float64}(undef, Ny1, Nx1, n_steps)
    APHI = Array{Float64}(undef, Ny1, Nx1, n_steps)
    timestep = Array{Float64}(undef, n_steps)
    dt = Array{Float64}(undef, n_steps)
    timesum = Array{Float64}(undef, n_steps)
    timesum_Ma = Array{Float64}(undef, n_steps)
    marknum = Array{Int}(undef, n_steps)
    max_T = Array{Float64}(undef, n_steps)

    xm = zeros(Float64, end_marknum, n_steps)
    ym = zeros(Float64, end_marknum, n_steps)
    tm = zeros(Int, end_marknum, n_steps)
    tkm = zeros(Float64, end_marknum, n_steps)
    phim = zeros(Float64, end_marknum, n_steps)
    XWsolidm0 = zeros(Float64, end_marknum, n_steps)
  
    @showprogress 1 "reading files..." for (i, f) in enumerate(files)
        jldopen(f, "r") do file
            RHO[:,:,i] = file["RHO"]
            ETA[:,:,i] = file["ETA"]
            vx[:, :, i] = file["vx"]
            vy[:, :, i] = file["vy"]
            vxp[:, :, i] = file["vxp"]
            vyp[:, :, i] = file["vyp"]
            PHI[:,:,i] = file["PHI"]
            EII[:,:,i] = file["EII"]
            SII[:,:,i] = file["SII"]
            gx[:,:,i] = file["gx"]
            gy[:,:,i] = file["gy"]
            tk2[:,:,i] = file["tk2"]
            HS[:,:,i] = file["HS"]
            HA[:,:,i] = file["HA"]
            pr[:, :, i] = file["pr"]
            pf[:, :, i] = file["pf"]
            KX[:,:,i] = file["KX"]
            qxD[:,:,i] = file["qxD"]
            qyD[:,:,i] = file["qyD"]
            RX[:,:,i] = file["RX"]
            ETAPHI[:,:,i] = file["ETAPHI"]
            DMP[:,:,i] = file["DMP"]
            DHP[:,:,i] = file["DHP"]
            XWS[:,:,i] = file["XWS"]
            APHI[:,:,i] = file["APHI"]
            timestep[i] = file["timestep"]
            dt[i] = file["dt"]
            timesum[i] = file["timesum"]
            timesum_Ma[i] = timesum[i] / (365.25 * 24 * 3600) * 1e-6
            max_T[i] = maximum(tk2[:,:,i])
            marknum[i] = file["marknum"]
            xm[1:marknum[i], i] = file["xm"]
            ym[1:marknum[i], i] = file["ym"]
            tm[1:marknum[i], i] = file["tm"]
            tkm[1:marknum[i], i] = file["tkm"]
            phim[1:marknum[i], i] = file["phim"]
            XWsolidm0[1:marknum[i], i] = file["XWsolidm0"]
        end
    end

    @info "plotting t vs max T"
    fig, ax = plt.subplots()
    ax.plot(timesum_Ma[2:end], max_T[2:end])
    ax.set_xlabel("time [Ma]")
    ax.set_ylabel("max T [K]")
    ax.set_title("Maximum temperature")
    fig.savefig(input_path*"/fig_t_maxT-phim0=$(phim0)_ratioAl26=$(ratio_al)_XWsolid0=$(XWsolidm_init[1]).pdf")
    plt.close()

    @info "plotting t vs Δρ"
    inner_r = round(Int, 0.1*length(radius_range))
    outer_r = round(Int, 0.9*length(radius_range))
    fig, ax = plt.subplots()
    Δρ = zero(timesum_Ma)
    for i in 1:1:n_steps
        Δρ[i] = abs(ravg(RHO[:, :, i], inner_r) - ravg(RHO[:, :, i], outer_r))
    end
    ax.plot(timesum_Ma[2:end], Δρ[2:end])
    ax.set_xlabel("time [Ma]")
    ax.set_ylabel("Δρ [kg/m³]")
    ax.set_title("Density contrast 0.1R - 0.9R") 
    fig.savefig(input_path*"/fig_t_deltaRho-phim0=$(phim0)_ratioAl26=$(ratio_al)_XWsolid0=$(XWsolidm_init[1]).pdf")
    plt.close()

    @info "plotting t vs T, Δρ"
    inner_r = round(Int, 0.1*length(radius_range))
    outer_r = round(Int, 0.9*length(radius_range))
    fig, ax1 = plt.subplots()
    color = "tab:blue"
    ax1.plot(timesum_Ma[2:end], Δρ[2:end], color=color)
    ax1.set_xlabel("time [Ma]")
    ax1.set_ylabel("Δρ [kg/m³]", color=color)
    ax1.set_title("Maximum temperature and density contrast 0.1R - 0.9R") 
    ax2 = ax1.twinx()
    color = "tab:red"
    ax2.plot(timesum_Ma[2:end], max_T[2:end], color=color)
    ax2.set_ylabel("max T [K]", color=color)
    fig.tight_layout()
    fig.savefig(input_path*"/fig_t_deltaRho-maxT-phim0=$(phim0)_ratioAl26=$(ratio_al)_XWsolid0=$(XWsolidm_init[1]).pdf")
    plt.close()

    @info "plotting t vs mean T"
    fig, ax = plt.subplots()
    mean_T = reshape(mean(tk2, dims=(1, 2)), n_steps)
    ax.plot(timesum_Ma[2:end], mean_T[2:end])
    ax.set_xlabel("time [Ma]")
    ax.set_ylabel("mean T [K]")
    ax.set_title("Mean temperature")
    fig.savefig(input_path*"/fig_t_meanT-phim0=$(phim0)_ratioAl26=$(ratio_al)_XWsolid0=$(XWsolidm_init[1]).pdf")
    plt.close()
# -----------------------------------------------------------------------------
    n_ages = 4
    idx_ages = trunc.(Int, range(20, n_steps, n_ages))
# -----------------------------------------------------------------------------
    @info "plotting mean T vs radius"
    fig, ax = plt.subplots()
    for idx in idx_ages
        age = timesum_Ma[idx]
        mean_T = ravg.(Ref(tk2[:, :, idx]), radius_range)
        ax.plot(mean_T, radius_range, label="$(sprint(f3p, age)) Ma")
    end
    ax.set_xlabel("mean T [K]")
    ax.set_ylabel("radius [km]")
    ax.set_title("Mean temperature")
    ax.legend()
    fig.savefig(input_path*"/fig_meanT_radius-phim0=$(phim0)_ratioAl26=$(ratio_al)_XWsolid0=$(XWsolidm_init[1]).pdf")
    plt.close()

    @info "plotting mean ρ vs radius"
    fig, ax = plt.subplots()
    for idx in idx_ages
        age = timesum_Ma[idx]
        mean_ρ = ravg.(Ref(RHO[:, :, idx]), radius_range)
        ax.plot(mean_ρ, radius_range, label="$(sprint(f3p, age)) Ma")
    end
    ax.set_xlabel("mean density [kg/m³]")
    ax.set_ylabel("radius [km]")
    ax.set_title("Mean density")
    ax.legend()
    fig.savefig(input_path*"/fig_meanRho_radius-phim0=$(phim0)_ratioAl26=$(ratio_al)_XWsolid0=$(XWsolidm_init[1]).pdf") 
    plt.close()

    @info "plotting radial Darcy velocity vs radius"
    fig, ax = plt.subplots()
    for idx in idx_ages
        age = timesum_Ma[idx]
        rad_qD = ravg.(Ref(vrad(qxD[:, :, idx], qyD[:, :, idx])), radius_range)
        ax.plot(rad_qD, radius_range, label="$(sprint(f3p, age)) Ma")
    end
    ax.set_xlabel("radial Darcy velocity [m/s]")
    ax.set_ylabel("radius [km]")
    ax.set_title("Mean radial Darcy velocity")
    ax.legend()
    fig.savefig(input_path*"/fig_radialDarcyvel_radius-phim0=$(phim0)_ratioAl26=$(ratio_al)_XWsolid0=$(XWsolidm_init[1]).pdf")
    plt.close()
# -----------------------------------------------------------------------------
    # cm = plt.get_cmap(:inferno)
    cm = plt.get_cmap(:plasma)
    t_idxs = round.(Int, LinRange(2, n_steps, 9))
# -----------------------------------------------------------------------------
    @info "plotting density panel"    
    fig, axs = plt.subplots(3, 3, figsize=(8, 8), constrained_layout=true)
    vmin, vmax = extrema(RHO[:, :, t_idxs])
    im = nothing
    for (ax, i) in zip(vcat(permutedims(axs)...), t_idxs)
        im = plot_field(ax, RHO[:, :, i], vmin, vmax, timesum_Ma[i], cm)
    end
    cbar = fig.colorbar(im, ax=axs, shrink=0.6)
    cbar.set_label("ρ [kg/m³]", fontsize=10)
    fig.suptitle("Density")
    fig.savefig(input_path*"/fig_rho_panel-phim0=$(phim0)_ratioAl26=$(ratio_al)_XWsolid0=$(XWsolidm_init[1]).pdf")
    plt.close()

    @info "plotting temperature panel"        
    fig, axs = plt.subplots(3, 3, figsize=(8, 8), constrained_layout=true)
    vmin, vmax = extrema(tk2[:, :, t_idxs])
    im = nothing
    for (ax, i) in zip(vcat(permutedims(axs)...), t_idxs)
        im = plot_field(ax, tk2[:, :, i], vmin, vmax, timesum_Ma[i], cm)
    end
    cbar = fig.colorbar(im, ax=axs, shrink=0.6)
    cbar.set_label("T [K]", fontsize=10)
    fig.suptitle("Temperature")
    fig.savefig(input_path*"/fig_T_panel-phim0=$(phim0)_ratioAl26=$(ratio_al)_XWsolid0=$(XWsolidm_init[1]).pdf")
    plt.close()

    @info "plotting porosity panel"    
    fig, axs = plt.subplots(3, 3, figsize=(8, 8), constrained_layout=true)
    vmin, vmax = extrema(PHI[:, :, t_idxs])
    im = nothing
    for (ax, i) in zip(vcat(permutedims(axs)...), t_idxs)
        im = plot_field(ax, PHI[:, :, i], vmin, vmax, timesum_Ma[i], cm)
    end
    cbar = fig.colorbar(im, ax=axs, shrink=0.6)
    cbar.set_label("ϕ", fontsize=10)
    fig.suptitle("Porosity")
    fig.savefig(input_path*"/fig_phi_panel-phim0=$(phim0)_ratioAl26=$(ratio_al)_XWsolid0=$(XWsolidm_init[1]).pdf")
    plt.close()

    @info "plotting wet silicate panel"        
    fig, axs = plt.subplots(3, 3, figsize=(8, 8), constrained_layout=true)
    vmin, vmax = extrema(XWS[:, :, t_idxs])
    im = nothing
    for (ax, i) in zip(vcat(permutedims(axs)...), t_idxs)
        im = plot_field(ax, XWS[:, :, i], vmin, vmax, timesum_Ma[i], cm)
    end
    cbar = fig.colorbar(im, ax=axs, shrink=0.6)
    cbar.set_label("XWˢ", fontsize=10)
    fig.suptitle("Wet silicate molar fraction in solid")
    fig.savefig(input_path*"/fig_XWS_panel-phim0=$(phim0)_ratioAl26=$(ratio_al)_XWsolid0=$(XWsolidm_init[1]).pdf")
    plt.close()

# -----------------------------------------------------------------------------
    levels = 250
    cm = plt.get_cmap(:jet, levels)
    marker_radius_margin = 500
    angle_step = 5
    t_begin = 1
    t_end = n_steps
    rs = [0.3, 0.6, 0.9] .* rplanet
    for r in rs
        rmin = r - marker_radius_margin
        rmax = r + marker_radius_margin
        azimuths, radii = get_marker_azimuths_radii(xm, ym, xcenter, ycenter)
        angles = range(0, step=angle_step, stop=359)
        selection = zero(angles)
        for (i, ϕ) in enumerate(angles)
            ϕ_markers = findall(
                (ϕ .<= azimuths[:, 1] .< ϕ+1) .& (rmin .< radii[:, 1] .<= rmax))
            if length(ϕ_markers) == 0
            ϕ_markers = findall(
                (ϕ-1 .<= azimuths[:, 1] .< ϕ+1) .& (rmin .< radii[:, 1] .<= rmax))
            end
            selection[i] = sample(rgen, ϕ_markers)
        end
        # X, Y = repeat(angles, inner=(1, size(timesum_Ma[t_begin:t_end], 1))), repeat(
        #     timesum_Ma[t_begin:t_end], inner=(1, size(angles, 1)))'

    # -------------------------------------------------------------------------
        @info "plotting marker porosity at r=$(r)m"
        fig = plt.figure(figsize=(12, 8), dpi=300)
        ax = fig.gca(projection="3d")
        im = ax.contour3D(
            angles, timesum_Ma[t_begin:t_end]', phim[selection, :]', levels, cmap=cm, alpha=0.5, antialiased=true)
        ax.set_ylim(reverse(ax.get_ylim()))
        ax.view_init(elev=30.0, azim=-50.0)
        ax.set_xlabel("azimuth [°]")
        ax.set_ylabel("time [Ma]")
        ax.set_zlabel("porosity")
        cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.2)
        cbar.set_label("porosity", fontsize=10)
        cbar.set_alpha(1.0)
        cbar.draw_all()
        fig.savefig(input_path*"/fig_marker_phi-r=$(r)_phim0=$(phim0)_ratioAl26=$(ratio_al)_XWsolid0=$(XWsolidm_init[1]).pdf")
        plt.close()

        @info "plotting marker temperature at r=$(r)m"
        fig = plt.figure(figsize=(12, 8), dpi=300)
        ax = fig.gca(projection="3d")
        im = ax.contour3D(
            angles, timesum_Ma[t_begin:t_end]', tkm[selection, :]', levels, cmap=cm, alpha=0.5, antialiased=true)
        ax.set_ylim(reverse(ax.get_ylim()))
        ax.view_init(elev=30.0, azim=-50.0)
        ax.set_xlabel("azimuth [°]")
        ax.set_ylabel("time [Ma]")
        ax.set_zlabel("temperature [K]")
        cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.2)
        cbar.set_label("temperature [K]", fontsize=10)
        cbar.set_alpha(1.0)
        cbar.draw_all()
        fig.savefig(input_path*"/fig_marker_T-r=$(r)_phim0=$(phim0)_ratioAl26=$(ratio_al)_XWsolid0=$(XWsolidm_init[1]).pdf")
        plt.close()

        @info "plotting marker radius at r=$(r)m"
        fig = plt.figure(figsize=(12, 8), dpi=300)
        ax = fig.gca(projection="3d")
        im = ax.contour3D(
            angles, timesum_Ma[t_begin:t_end]', radii[selection, :]', levels, cmap=cm, alpha=0.5, antialiased=true)
        ax.set_ylim(reverse(ax.get_ylim()))
        ax.view_init(elev=30.0, azim=-50.0)
        ax.set_xlabel("azimuth [°]")
        ax.set_ylabel("time [Ma]")
        ax.set_zlabel("radius [m]")
        cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.1)
        cbar.set_label("radius [m]", fontsize=10)
        cbar.set_alpha(1.0)
        cbar.draw_all()
        fig.savefig(input_path*"/fig_marker_r-r=$(r)_phim0=$(phim0)_ratioAl26=$(ratio_al)_XWsolid0=$(XWsolidm_init[1]).pdf")
        plt.close()

        @info "plotting marker XWˢ at r=$(r)m"
        fig = plt.figure(figsize=(12, 8), dpi=300)
        ax = fig.gca(projection="3d")
        im = ax.contour3D(
            angles, timesum_Ma[t_begin:t_end]', XWsolidm0[selection, :]', levels, cmap=cm, alpha=0.5, antialiased=true)
        ax.set_ylim(reverse(ax.get_ylim()))
        ax.view_init(elev=30.0, azim=-50.0)
        ax.set_xlabel("azimuth [°]")
        ax.set_ylabel("time [Ma]")
        ax.set_zlabel("XWˢ")
        cbar = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.1)
        cbar.set_label("XWˢ", fontsize=10)
        cbar.set_alpha(1.0)
        cbar.draw_all()
        fig.savefig(input_path*"/fig_marker_XWS-r=$(r)_phim0=$(phim0)_ratioAl26=$(ratio_al)_XWsolid0=$(XWsolidm_init[1]).pdf")
        plt.close()
    end
end



"""
Parse command line arguments and feed them to the main function.

$(SIGNATURES)

# Details:
    
    - nothing

# Returns

    - parsed_args: parsed command line arguments
"""
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table! s begin
        "input_path"
            help = "input path where simulation output is stored"
            required = true
    end
    return parse_args(s)
end

"""
Runs the plotting routine.

$(SIGNATURES)

# Details

    - nothing

# Returns

    - nothing 
"""
function main()
    parsed_args = parse_commandline()
    input_path = parsed_args["input_path"]
    if !isdir(input_path)
        throw(ArgumentError("input_path must be a valid directory"))
    end
    @info "reading from $input_path"
    generate_plots(input_path)
end

main()
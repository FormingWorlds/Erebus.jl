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
    ycenter = round(Int, ysize / 2)
    xcenter = round(Int, xsize / 2)
    r = min(r, ycenter, xcenter, ysize-ycenter, xsize-xcenter)
    circle = zeros(Bool, ysize, xsize)
    for ϕ = 0:0.01:2π
        y = round(Int, ycenter + r*cos(ϕ))
        x = round(Int, xcenter + r*sin(ϕ))
        if 1<=y<=ysize && 1<=x<=xsize
            @inbounds circle[y, x] = true
        end
    end
    @inbounds return mean(f[circle])
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
    ycenter = round(Int, ysize / 2)
    xcenter = round(Int, xsize / 2)
    r_max = max(round(Int, xsize/2), round(Int, ysize/2)) * sqrt(2.0)
    for ϕ = 0.0:0.01:2π, r = 1:1:r_max
        y = round(Int, ycenter + r*cos(ϕ))
        x = round(Int, xcenter + r*sin(ϕ))
        if 1<=y<=ysize && 1<=x<=xsize
            dx = x - xcenter
            dy = y - ycenter
            @inbounds vrad[y, x] = (
                dx*vx[y, x] + dy*vy[y, x]) / sqrt(dx^2 + dy^2)
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
    - t: time step of the field [Myr]
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
        ax.set_title("t = $(sprint(f3p, t)) Myr", fontsize=fontsize)
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
Return boolean mask indicating which grid points are within a given radius of a
circular planet.

$(SIGNATURES)

# Details

    - xsize: number of grid points in x-direction
    - ysize: number of grid points in y-direction
    - dx: grid spacing in x-direction [m]
    - dy: grid spacing in y-direction [m]
    - r: radius of the planet [m]
    - xcenter: x-coordinate of the center of the planet [m]
    - ycenter: y-coordinate of the center of the planet [m]

# Returns

    - mask: boolean mask indicating which grid points are within the planet
"""
function get_planet_mask(xsize, ysize, dx, dy, r, xcenter, ycenter)
    mask = zeros(Bool, ysize, xsize)
    for y = 1:ysize, x = 1:xsize
        Δx = x*dx - xcenter
        Δy = y*dy - ycenter
        mask[y, x] = sqrt(Δx^2 + Δy^2) <= r
    end
    return mask
end

"""
Create additional plots from results stored in `.jld2` simulation output files
 across six specific simulation runs.

$(SIGNATURES)

# Details

    - paths: vector of paths to the directories containing the simulation output files in the following order:
        - run "Iyer et al. (2012)"
        - run "Bland and Travis (2017)"
        - run "Travis et al. (2018)"
        - run "Gerya (2019)"
        - run "modified Travis et al. (2018)"
        - run "modified Gerya (2019)"

# Returns

    - nothing
"""
function generate_special_plots(paths)
    input_path = paths[1]
    titles = [
        "Iyer et al. (2012)",
        "Bland and Travis (2017)",
        "Travis et al. (2018)",
        "Gerya (2019)",
        "modified Travis et al. (2018)",
        "modified Gerya (2019)"
    ]
    files = filter!(
        x -> isfile(x) && endswith(x, ".jld2"), readdir(input_path, join=true))
    n_steps = max_n_steps = length(files)
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
    # end_marknum = file["marknum"]
    phim0 = file["phim0"]
    ratio_al = file["ratio_al"]
    XWsolidm_init = file["XWsolidm_init"]
    # rrcoef = file["reaction_rate_coeff_mode"]
    # mkrprop = file["marker_property_mode"]
    close(file)

    for input_path in paths[2:end]
        files = filter!(
        x -> isfile(x) && endswith(x, ".jld2"), readdir(input_path, join=true))
        n_steps = min(length(files), n_steps)
        max_n_steps = max(length(files), max_n_steps)
    end

    radius_max = trunc(Int, min(rplanet/dx, rplanet/dy))
    radius_range = collect(1:1:radius_max)
    n_ticks = 4

    planet_NxNy = zeros(Ny, Nx)
    planet_Nx1Ny1 = zeros(Ny1, Nx1)

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

    HR = Array{Float64}(undef, Ny1, Nx1, max_n_steps, 6)
    HS = Array{Float64}(undef, Ny1, Nx1, max_n_steps, 6)
    HA = Array{Float64}(undef, Ny1, Nx1, max_n_steps, 6)
    DHP = Array{Float64}(undef, Ny1, Nx1, max_n_steps, 6)
    timesum_Myr = Array{Float64}(undef, max_n_steps)

    for (ii, input_path) in enumerate(paths)
        files = filter!(
            x -> isfile(x) && endswith(x, ".jld2"), readdir(input_path, join=true))
        @showprogress 1 "reading files $ii" for (i, f) in enumerate(files)
            jldopen(f, "r") do file
                HR[:, :, i, ii] = file["HR"]
                HS[:,:,i, ii] = file["HS"]
                HA[:,:,i, ii] = file["HA"]
                DHP[:,:,i, ii] = file["DHP"]
                timesum_Myr[i] = file["timesum"] / (365.25 * 24 * 3600) * 1e-6
            end
            end
    end

# -----------------------------------------------------------------------------
    planet_NxNy = get_planet_mask(Nx, Ny, dx, dy, rplanet, xcenter, ycenter)
    planet_Nx1Ny1 = get_planet_mask(Nx1, Ny1, dx, dy, rplanet, xcenter, ycenter)
    inner_f = 0.1
    outer_f = 0.9
# 
# -----------------------------------------------------------------------------
    @info "plotting heating terms"    
    fig, axs = plt.subplots(
        2, 3, figsize=(12, 8), constrained_layout=true, sharex=true, sharey=true)
    HR_vmin, HR_vmax = extrema(HR)
    HS_vmin, HS_vmax = extrema(HS)
    HA_vmin, HA_vmax = extrema(HA)
    DHP_min, DHP_max = extrema(DHP)
    # vmin = min(HR_vmin, HS_vmin, HA_vmin, DHP_min)
    vmin =10e-22
    vmax = max(HR_vmax, HS_vmax, HA_vmax, DHP_max)
    HR_avg_t = reshape(mean(HR[planet_Nx1Ny1,:,:], dims=1), :, 6)
    HA_avg_t = reshape(mean(HA[planet_Nx1Ny1,:,:], dims=1), :, 6)
    HS_avg_t = reshape(mean(HS[planet_Nx1Ny1,:,:], dims=1), :, 6)
    DHP_avg_t = reshape(mean(DHP[planet_Nx1Ny1,:,:], dims=1), :, 6)
    handles = labels = nothing
    for (i, ax) in enumerate(vcat(permutedims(axs)...))
        l1, = ax.plot(
            timesum_Myr[2:n_steps], HR_avg_t[2:n_steps, i], label="Hᵣ (radioactive)")
        l2, = ax.plot(
            timesum_Myr[2:n_steps], HA_avg_t[2:n_steps, i], label="Hₐ (adiabatic)")
        l3, = ax.plot(
            timesum_Myr[2:n_steps], HS_avg_t[2:n_steps, i], label="Hₛ (shear)")
        l4, = ax.plot(
            timesum_Myr[2:n_steps], DHP_avg_t[2:n_steps, i], label="HL (reaction latent)")
        if i>3; ax.set_xlabel("time [Myr]"); end
        if i==1||i==4; ax.set_ylabel("heat production [Wm⁻³]"); end
        ax.set_ylim([vmin, vmax])
        ax.set_title(titles[i])
        ax.set_yscale("log")
        handles, labels = ax.get_legend_handles_labels()
    end
    legend = fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, -0.05), fancybox=true, shadow=false, ncol=4)
    supertitle = fig.suptitle("Heat production contributions")
    fig.savefig("fig_heat-production_phim0=$(phim0)_ratioAl26=$(ratio_al)_XWsolid0=$(XWsolidm_init[1]).pdf", bbox_extra_artists=[legend,supertitle], bbox_inches="tight")
    plt.close()
    return nothing
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
        "input_path_1"
            help = "input path where simulation output 1 is stored"
            required = true
        "input_path_2"
            help = "input path where simulation output 2 is stored"
            required = true
        "input_path_3"
            help = "input path where simulation output 3 is stored"
            required = true
        "input_path_4"
            help = "input path where simulation output 4 is stored"
            required = true
        "input_path_5"
            help = "input path where simulation output 5 is stored"
            required = true
        "input_path_6"
            help = "input path where simulation output 6 is stored"
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
    input_path_1 = parsed_args["input_path_1"]
    input_path_2 = parsed_args["input_path_2"]
    input_path_3 = parsed_args["input_path_3"]
    input_path_4 = parsed_args["input_path_4"]
    input_path_5 = parsed_args["input_path_5"]
    input_path_6 = parsed_args["input_path_6"]
    if !isdir(input_path_1)
        throw(ArgumentError("input_path must be a valid directory"))
    end
    @info "reading from $input_path_1, $input_path_2, $input_path_3, $input_path_4, $input_path_5, $input_path_6"
    generate_special_plots([input_path_1, input_path_2, input_path_3, input_path_4, input_path_5, input_path_6])
end

main()
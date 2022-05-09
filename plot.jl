using ArgParse
using DocStringExtensions   
using JLD2
using LinearAlgebra
using Plots
using Printf

theme(:dark)

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
    limits = extrema(grid)
    cbartickrange = range(limits..., n_ticks)
    cbarticklabels = map(t -> @sprintf("%.1E", t), cbartickrange)
    return limits, (cbartickrange, cbarticklabels)
end

"""
Plot results from `.jld2` simulation output files.

$(SIGNATURES)

# Details

    - input_path: path to the directory containing the simulation output files

# Returns

    - nothing
"""
function plot_results(input_path)
    files = filter!(
        x -> isfile(x) && endswith(x, ".jld2"), readdir(input_path, join=true))
    n_steps = length(files)
    file = jldopen(files[end], "r")
    
    Nx = file["Nx"]
    Ny = file["Ny"]
    Nx1 = file["Nx1"]
    Ny1 = file["Ny1"]
    marknum = file["marknum"]

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

    # xm = Array{Float64}(undef, marknum, n_steps)
    # ym = Array{Float64}(undef, marknum, n_steps)
    # tm = Array{Int}(undef, marknum, n_steps)

    for (i, f) in enumerate(files)
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
            # xm[:, i] = file["xm"]
            # ym[:, i] = file["ym"]
            # tm[:, i] = file["tm"]
        end
    end

    gmag = Array{Float64}(undef, Ny1, Nx1, n_steps)
    @views @. gmag[2:Ny, 2:Nx, :] = sqrt((0.5 * (gx[2:Ny, 2:Nx, :] + gx[2:Ny, 1:Nx-1, :]))^2 + (0.5 * (gy[2:Ny, 2:Nx, :] + gy[1:Ny-1, 2:Nx, :]))^2)

    # limits_pr, cbarticks_pr = colorbar_properties(pr)
    
    gr(fmt = :png)

    anim_1 = @animate for i in 1:n_steps
        @views RHO_ = RHO[:, :, i]
        @views ETA_ = ETA[:, :, i]
        @views vx_ = vx[:, :, i]
        @views vy_ = vy[:, :, i]
        @views PHI_ = PHI[:, :, i]
        @views EII_ = EII[2:Ny, 2:Nx, i]
        @views SII_ = SII[2:Ny, 2:Nx, i]
        @views gmag_ = gmag[2:Ny, 2:Nx, i]
        @views tk2_ = tk2[:, :, i]
        @views pr_ = pr[:, :, i]
        # @views xm_ = xm[1:100:end, i]
        # @views ym_ = ym[1:100:end, i]
        # @views tm_ = tm[1:100:end, i]
        A = heatmap(xp, yp, RHO_; title="RHO", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p)
        B = heatmap(x, y, ETA_; title="ETA", aspect_ratio=:equal, cmap=:jet, xticks=xticks_b, yticks=yticks_b, xlim=xlim_b, ylim=ylim_b)
        C = heatmap(xvx, yvx, vx_; title="vx", aspect_ratio=:equal, cmap=:jet, xticks=xticks_vx, yticks=yticks_vx, xlim=xlim_vx, ylim=ylim_vx)
        D = heatmap(xvy, yvy, vy_; title="vy", aspect_ratio=:equal, cmap=:jet,
        xticks=xticks_vy, yticks=yticks_vy, xlim=xlim_vy, ylim=ylim_vy)
        E = heatmap(
            xp, yp, log10.(PHI_ ./ (1 .- PHI_)); title="log₁₀(ϕ/(1-ϕ))", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p)
        F = heatmap(
            xp[2:Nx], yp[2:Nx], log10.(EII_); title="log₁₀(Eᴵᴵ)", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p)
        G = heatmap(
            xp[2:Nx], yp[2:Nx], log10.(SII_); title="log₁₀(Sᴵᴵ)", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p)
        H = heatmap(
            xp[2:Nx], yp[2:Nx], gmag_; title="gmag", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p)
        I = heatmap(
            xp, yp, tk2_; title="T", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p)
        layout = @layout [
            [a b c]
            [d e f]
            [g h i]
        ]
        plot(A, B, C, D, E, F, G, H, I; layout, size=(1920, 1920)) 
    end
    gif(anim_1, input_path*"/HydrologyPlanetesimals_1.mp4")        

    X, Y = meshgrid(xp[10:20:Nx1], yp[10:20:Ny1])

    anim_2 = @animate for i in 1:n_steps
        @views ETA_ = ETA[:, :, i]
        @views pr_ = pr[:, :, i]
        @views vxp_ = vxp[:, :, i]
        @views vyp_ = vyp[:, :, i]
        @views tk2_ = tk2[:, :, i]
        @views HS_ = HS[:, :, i]
        @views HA_ = HA[:, :, i]
        @views pr_ = pr[:, :, i]
        @views pf_ = pf[:, :, i]
        @views KX_ = KX[:, :, i]
        @views PHI_ = PHI[:, :, i]
        @views qxD_ = qxD[:, :, i]
        @views qyD_ = qyD[:, :, i]
        @views RX_ = RX[:, :, i]
        @views ETAPHI_ = ETAPHI[:, :, i]
        @views RHO_ = RHO[:, :, i]
        A = heatmap(x, y, log10.(ETA_); title="log₁₀(ETA)", aspect_ratio=:equal, cmap=:jet, xticks=xticks_b, yticks=yticks_b, xlim=xlim_b, ylim=ylim_b)
        quiver!(X, Y, quiver=(vxp_[10:20:Nx1, 10:20:Ny1], vyp_[10:20:Nx1, 10:20:Ny1]))
        B = heatmap(xp, yp, pr_; title="pr", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p)
        C = heatmap(xp, yp, vxp_; title="vxp", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p)
        D = heatmap(xp, yp, vyp_; title="vyp", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p)
        E = heatmap(xp, yp, HS_; title="HS", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p)
        F = heatmap(xp, yp, HA_; title="HA", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p)
        G = heatmap(xp, yp, pr_ .- pf_; title="Ptotal - Pfluid", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p)
        H = heatmap(xp, yp, log10.(KX_); title="log₁₀(KX)", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p)
        I = heatmap(xp, yp, tk2_; title="T", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p)
        J = heatmap(xp, yp, log10.(PHI_ ./ (1 .- PHI_)); title="log₁₀(ϕ/(1-ϕ))", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p)
        K = heatmap(xp, yp, pf_; title="Pfluid", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p)
        L = heatmap(xvx, yvx, qxD_; title="qxD", aspect_ratio=:equal, cmap=:jet, xticks=xticks_vx, yticks=yticks_vx, xlim=xlim_vx, ylim=ylim_vx)
        M = heatmap(xvy, yvy, qyD_; title="qyD", aspect_ratio=:equal, cmap=:jet, xticks=xticks_vy, yticks=yticks_vy, xlim=xlim_vy, ylim=ylim_vy)
        N = heatmap(xvx, yvx, log10.(RX_); title="log₁₀(RX)", aspect_ratio=:equal, cmap=:jet, xticks=xticks_vx, yticks=yticks_vx, xlim=xlim_vx, ylim=ylim_vx)
        O = heatmap(xp, yp, log10.(ETAPHI_); title="log₁₀(ETAPHI)", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p)
        P = heatmap(xp, yp, RHO_; title="RHO", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p)

        layout = @layout [
            [a b c d]
            [e f g h]
            [i j k l]
            [m n o p]
        ]
        plot(
            A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P; layout, size=(1920, 1920)) 
    end
    gif(anim_2, input_path*"/HydrologyPlanetesimals_2.mp4")   
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
            help = "input path where simulation data is stored"
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
    plot_results(input_path)
end

main()

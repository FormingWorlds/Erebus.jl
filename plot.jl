using ArgParse
using DocStringExtensions   
using JLD2
using LinearAlgebra
using Plots
using Printf
using ProgressMeter

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
    limits = get_extrema(grid)
    cbartickrange = range(limits..., n_ticks)
    cbarticklabels = map(t -> @sprintf("%.1E", t), cbartickrange)
    return limits, (cbartickrange, cbarticklabels)
end


function get_extrema(grid)
    a, b = extrema(grid)
    if a==b ||  (a==-Inf && b==Inf)
        return 0.0, 1.0
    else
        return a, b
    end
end

"""
Create movies of results from `.jld2` simulation output files.

$(SIGNATURES)

# Details

    - input_path: path to the directory containing the simulation output files

# Returns

    - nothing
"""
function movie_results(input_path)
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
    close(file)

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
    marknum = Array{Float64}(undef, n_steps)
    max_T = Array{Float64}(undef, n_steps)

    # xm = Array{Float64}(undef, marknum, n_steps)
    # ym = Array{Float64}(undef, marknum, n_steps)
    # tm = Array{Int}(undef, marknum, n_steps)

    @showprogress 1 "Reading files..." for (i, f) in enumerate(files)
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
            # dtm[i] = file["dtm"]
            dt[i] = file["dt"]
            timesum[i] = file["timesum"]
            timesum_Ma[i] = timesum[i] / (365.25 * 24 * 3600) * 1e-20
            marknum[i] = file["marknum"]
            max_T[i] = maximum(tk2[:,:,i])
            # xm[:, i] = file["xm"]
            # ym[:, i] = file["ym"]
            # tm[:, i] = file["tm"]
        end
    end
    min_RHO, max_RHO = get_extrema(RHO[:,:,3:end])
    min_ETA, max_ETA = get_extrema(ETA[:,:,3:end])
    min_vx, max_vx = get_extrema(vx[:,:,3:end])
    min_vy, max_vy = get_extrema(vy[:,:,3:end])
    min_vxp, max_vxp = get_extrema(vxp[:,:,3:end])
    min_vyp, max_vyp = get_extrema(vyp[:,:,3:end])
    min_PHI, max_PHI = get_extrema(PHI[:,:,3:end])
    min_tk2, max_tk2 = get_extrema(tk2[:,:,3:end])
    min_HS, max_HS = get_extrema(HS[:,:,3:end])
    min_HA, max_HA = get_extrema(HA[:,:,3:end])
    min_pr, max_pr = get_extrema(pr[:,:,3:end])
    min_pf, max_pf = get_extrema(pf[:,:,3:end])
    min_qxD, max_qxD = get_extrema(qxD[:,:,3:end])
    min_qyD, max_qyD = get_extrema(qyD[:,:,3:end])
    min_DMP, max_DMP = get_extrema(DMP[:,:,3:end])
    min_DHP, max_DHP = get_extrema(DHP[:,:,3:end])
    min_XWS, max_XWS = get_extrema(XWS[:,:,3:end])
    min_APHI, max_APHI = get_extrema(APHI[:,:,3:end])
    
    gmag = Array{Float64}(undef, Ny1, Nx1, n_steps)
    @views @. gmag[2:Ny, 2:Nx, :] = sqrt((0.5 * (gx[2:Ny, 2:Nx, :] + gx[2:Ny, 1:Nx-1, :]))^2 + (0.5 * (gy[2:Ny, 2:Nx, :] + gy[1:Ny-1, 2:Nx, :]))^2)

    log_PHI1mPHI = log10.(PHI ./ (1 .- PHI))
    log_EII = log10.(EII)
    log_SII = log10.(SII)
    log_ETA = log10.(ETA)
    prmpf = pr .- pf 
    log_KX = log10.(KX)
    log_RX = log10.(RX)
    log_ETAPHI = log10.(ETAPHI)
    # log_XWS = log10.(XWS)


    min_gmax, max_gmag = get_extrema(gmag[:,:,3:end])
    min_log_PHI1mPHI, max_log_PHI1mPHI = get_extrema(log_PHI1mPHI[:,:,3:end])
    min_log_EII, max_log_EII = get_extrema(log_EII[2:end-1,2:end-1,3:end])
    min_log_SII, max_log_SII = get_extrema(log_SII[2:end-1,2:end-1,3:end])
    min_log_ETA, max_log_ETA = get_extrema(log_ETA[:,:,3:end])
    min_prmpf, max_prmpf = get_extrema(prmpf[:,:,3:end])
    min_log_KX, max_log_KX = get_extrema(log_KX[:,:,3:end])
    min_log_RX, max_log_RX = get_extrema(log_RX[:,:,3:end])
    min_log_ETAPHI, max_log_ETAPHI = get_extrema(log_ETAPHI[:,:,3:end])
    # min_log_XWS, max_log_XWS = get_extrema(log_XWS[:,:,3:end])
    
    # limits_pr, cbarticks_pr = colorbar_properties(pr)
    
    gr(fmt = :png)

    invalid(x) = all(iszero, x) || all(isnan, x) || all(isinf, x)

    p_1 = Progress(n_steps, 1, "Assembling animation 1...")
    p_2 = Progress(n_steps, 1, "Assembling animation 2...")

    anim_1 = @animate for i in 1:n_steps
        @views RHO_ = RHO[:, :, i]
        @views ETA_ = ETA[:, :, i]
        @views vx_ = vx[:, :, i]
        @views vy_ = vy[:, :, i]
        @views log_PHI1mPHI_ = log_PHI1mPHI[:, :, i]
        @views log_EII_ = log_EII[2:Ny, 2:Nx, i]
        @views log_SII_ = log_SII[2:Ny, 2:Nx, i]
        @views gmag_ = gmag[2:Ny, 2:Nx, i]
        @views tk2_ = tk2[:, :, i]
        @views pr_ = pr[:, :, i]
        @views DMP_ = DMP[:, :, i]
        @views DHP_ = DHP[:, :, i]
        @views XWS_ = XWS[:, :, i]
        # @views log_XWS_ = log_XWS[:, :, i]
        # @views xm_ = xm[1:100:end, i]
        # @views ym_ = ym[1:100:end, i]
        # @views tm_ = tm[1:100:end, i]
        if invalid(RHO_)
            RHO_[2, 2] = 1e-20
        end
        A = heatmap(xp, yp, RHO_; title="RHO", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_RHO, max_RHO))
        if invalid(ETA_)
            ETA_[2, 2] = 1e-20
        end
        B = heatmap(x, y, ETA_; title="ETA", aspect_ratio=:equal, cmap=:jet, xticks=xticks_b, yticks=yticks_b, xlim=xlim_b, ylim=ylim_b, clim=(min_ETA, max_ETA))
        if invalid(vx_)
            vx_[2, 2] = 1e-20
        end
        C = heatmap(xvx, yvx, vx_; title="vx", aspect_ratio=:equal, cmap=:jet, xticks=xticks_vx, yticks=yticks_vx, xlim=xlim_vx, ylim=ylim_vx, clim=(min_vx, max_vx))
        if invalid(vy_)
            vy_[2, 2] = 1e-20
        end
        D = heatmap(xvy, yvy, vy_; title="vy", aspect_ratio=:equal, cmap=:jet,
        xticks=xticks_vy, yticks=yticks_vy, xlim=xlim_vy, ylim=ylim_vy, clim=(min_vy, max_vy))
        if invalid(log_PHI1mPHI_)
            log_PHI1mPHI_[2, 2] = 1e-20
        end
        E = heatmap(
            xp, yp, log_PHI1mPHI_; title="log₁₀(ϕ/(1-ϕ))", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_log_PHI1mPHI, max_log_PHI1mPHI))
        if invalid(log_EII_)
            log_EII_[2, 2] = 1e-20
        end
        F = heatmap(
            xp[2:Nx], yp[2:Ny], log_EII_; title="log₁₀(Eᴵᴵ)", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_log_EII, max_log_EII))
        if invalid(log_SII_)
            log_SII_[2, 2] = 1e-20
        end
        G = heatmap(
            xp[2:Nx], yp[2:Ny], log_SII_; title="log₁₀(Sᴵᴵ)", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_log_SII, max_log_SII))
        if invalid(gmag_)
            gmag_[2, 2] = 1e-20
        end
        H = heatmap(
            xp[2:Nx], yp[2:Nx], gmag_; title="gmag", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_gmax, max_gmag))
        if invalid(tk2_)
            tk2_[2, 2] = 1e-20
        end
        I = heatmap(
            xp, yp, tk2_; title="T", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_tk2, max_tk2))
        if invalid(DMP_)
            DMP_[2, 2] = 1e-20
        end
        J= heatmap(
            xp, yp, DMP_; title="DMP", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_DMP, max_DMP))
        if invalid(DHP_)
            DHP_[2, 2] = 1e-20
        end
        K = heatmap(
            xp, yp, DHP_; title="DHP", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_DHP, max_DHP))
        if invalid(XWS_)
            XWS_[2, 2] = 1e-20
        end
        L = heatmap(
            xp, yp, XWS_; title="XWS", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_XWS, max_XWS))
        layout = @layout [
            [a b c j]
            [d e f k]
            [g h i l]
        ]
        plot(A, B, C, D, E, F, G, H, I, J, K, L; layout, size=(2560, 1920),
        plot_title="t_Ma=$(timesum_Ma[i]), dt=$(dt[i]), step=$(timestep[i]), marknum=$(marknum[i]), maxT_K=$(max_T[i])")
        next!(p_1) 
    end
    gif(anim_1, input_path*"/HydrologyPlanetesimals_1.mp4")        

    # X, Y = meshgrid(xp[10:20:Nx1], yp[10:20:Ny1])

    anim_2 = @animate for i in 1:n_steps
        @views log_ETA_ = log_ETA[:, :, i]
        @views log_PHI1mPHI_ = log_PHI1mPHI[:, :, i]
        @views pr_ = pr[:, :, i]
        @views vxp_ = vxp[:, :, i]
        @views vyp_ = vyp[:, :, i]
        @views tk2_ = tk2[:, :, i]
        @views HS_ = HS[:, :, i]
        @views HA_ = HA[:, :, i]
        @views pr_ = pr[:, :, i]
        @views pf_ = pf[:, :, i]
        @views prmpf_ = prmpf[:, :, i]
        @views log_KX_ = log_KX[:, :, i]
        @views PHI_ = PHI[:, :, i]
        @views qxD_ = qxD[:, :, i]
        @views qyD_ = qyD[:, :, i]
        @views log_RX_ = log_RX[:, :, i]
        @views log_ETAPHI_ = log_ETAPHI[:, :, i]
        @views RHO_ = RHO[:, :, i]
        @views APHI_ = APHI[:, :, i]
        
        if invalid(log_ETA_)
            log_ETA_[2, 2] = 1e-20
        end
        A = heatmap(x, y, log_ETA_; title="log₁₀(ETA)", aspect_ratio=:equal, cmap=:jet, xticks=xticks_b, yticks=yticks_b, xlim=xlim_b, ylim=ylim_b, clim=(min_log_ETA, max_log_ETA))
        # quiver!(X, Y, quiver=(vxp_[10:20:Nx1, 10:20:Ny1], vyp_[10:20:Nx1, 10:20:Ny1]))
        if invalid(pr_)
            pr_[2, 2] = 1e-20
        end
        B = heatmap(xp, yp, pr_; title="pr", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_pr, max_pr))
        if invalid(vxp_)
            vxp_[2, 2] = 1e-20
        end
        C = heatmap(xp, yp, vxp_; title="vxp", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_vxp, max_vxp))
        if invalid(vyp_)
            vyp_[2, 2] = 1e-20
        end
        D = heatmap(xp, yp, vyp_; title="vyp", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_vyp, max_vyp))
        if invalid(HS_)
            HS_[2, 2] = 1e-20
        end
        E = heatmap(xp, yp, HS_; title="HS", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_HS, max_HS))
        if invalid(HA_)
            HA_[2, 2] = 1e-20
        end
        F = heatmap(xp, yp, HA_; title="HA", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_HA, max_HA))
        if invalid(prmpf_)
            prmpf_[2, 2] = 1e-20
        end
        G = heatmap(xp, yp, prmpf_; title="Ptotal - Pfluid", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_prmpf, max_prmpf))
        if invalid(log_KX_)
            log_KX_[2, 2] = 1e-20
        end
        H = heatmap(xp, yp, log_KX_; title="log₁₀(KX)", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_log_KX, max_log_KX))
        if invalid(tk2_)
            tk2_[2, 2] = 1e-20
        end
        I = heatmap(xp, yp, tk2_; title="T", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_tk2, max_tk2))
        if invalid(log_PHI1mPHI_)
            PHI1mPHI_[2, 2] = 1e-20
        end
        J = heatmap(xp, yp, log_PHI1mPHI_; title="log₁₀(ϕ/(1-ϕ))", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_log_PHI1mPHI, max_log_PHI1mPHI))
        if invalid(pf_)
            pf_[2, 2] = 1e-20
        end
        K = heatmap(xp, yp, pf_; title="Pfluid", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_pf, max_pf))
        if invalid(qxD_)
            qxD_[2, 2] = 1e-20
        end
        L = heatmap(xvx, yvx, qxD_; title="qxD", aspect_ratio=:equal, cmap=:jet, xticks=xticks_vx, yticks=yticks_vx, xlim=xlim_vx, ylim=ylim_vx, clim=(min_qxD, max_qxD))
        if invalid(qyD_)
            qyD_[2, 2] = 1e-20
        end
        M = heatmap(xvy, yvy, qyD_; title="qyD", aspect_ratio=:equal, cmap=:jet, xticks=xticks_vy, yticks=yticks_vy, xlim=xlim_vy, ylim=ylim_vy, clim=(min_qyD, max_qyD))
        if invalid(log_RX_)
            log_RX_[2, 2] = 1e-20
        end
        N = heatmap(xvx, yvx, log_RX_; title="log₁₀(RX)", aspect_ratio=:equal, cmap=:jet, xticks=xticks_vx, yticks=yticks_vx, xlim=xlim_vx, ylim=ylim_vx, clim=(min_log_RX, max_log_RX))
        if invalid(log_ETAPHI_)
            log_ETAPHI_[2, 2] = 1e-20
        end
        O = heatmap(xp, yp, log_ETAPHI_; title="log₁₀(ETAPHI)", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_log_ETAPHI, max_log_ETAPHI))
        if invalid(RHO_)
            RHO_[2, 2] = 1e-20
        end
        P = heatmap(xp, yp, RHO_; title="RHO", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_RHO, max_RHO))
        if invalid(APHI_)
            APHI_[2, 2] = 1e-20
        end
        Q = heatmap(
            xp, yp, APHI_; title="APHI", aspect_ratio=:equal, cmap=:jet, xticks=xticks_p, yticks=yticks_p, xlim=xlim_p, ylim=ylim_p, clim=(min_APHI, max_APHI))

        layout = @layout [
            [a b c d]
            [e f g h]
            [i j k l]
            [m n o p q]
        ]
        plot(
            A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q; layout, size=(2560, 1920),
            plot_title="t_Ma=$(timesum_Ma[i]), dt=$(dt[i]), step=$(timestep[i]), marknum=$(marknum[i]), maxT_K=$(max_T[i])")
        next!(p_2)
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
    @info "reading from $input_path"
    movie_results(input_path)
end

main()

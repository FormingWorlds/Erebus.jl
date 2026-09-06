#!/usr/bin/env julia
"""
Export 2D hydrothermal reaction benchmark time series and final 2D fields.
"""

using Erebus
using JLD2
using Printf
using JSON

function export_reaction_data()
    output_dir = length(ARGS) >= 1 ? ARGS[1] : "output_hydrothermal_reaction_on_32"
    files = sort(
        filter(
            f ->
                startswith(f, "output_") &&
                endswith(f, ".jld2") &&
                f != "output_00000.jld2",
            readdir(output_dir),
        ),
    )
    if isempty(files)
        error("No output checkpoint files found in $output_dir")
    end

    time_Ma = Float64[]
    mean_T = Float64[]
    max_T = Float64[]
    mean_XW = Float64[]
    max_XW = Float64[]
    mean_q = Float64[]
    water_solid = Float64[]
    water_fluid = Float64[]

    stride = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : max(1, div(length(files), 250))
    selected_indices = unique(vcat(collect(1:stride:length(files)), [length(files)]))
    selected_files = files[selected_indices]
    println("Processing $(length(selected_files)) checkpoints out of $(length(files)) (stride=$stride)...")

    # Read geometry once from first file
    local dx, dy, rplanet, xcenter, ycenter, xp, yp, mask
    jldopen(joinpath(output_dir, selected_files[1]), "r") do f
        dx = f["dx"]
        dy = f["dy"]
        rplanet = f["rplanet"]
        xcenter = f["xcenter"]
        ycenter = f["ycenter"]
        xp = f["xp"]
        yp = f["yp"]
        mask = [(x - xcenter)^2 + (y - ycenter)^2 <= rplanet^2 for y in yp, x in xp]
    end

    mask_count = sum(mask)
    MH2O = 0.01801528
    MD = 0.031548
    rho_f = 1000.0
    rho_s = 3300.0

    for (idx, file) in enumerate(selected_files)
        path = joinpath(output_dir, file)
        jldopen(path, "r") do f
            ts = f["timesum"]
            push!(time_Ma, ts / (365.25 * 86400 * 1e6))

            phi = f["PHI"]
            tk = f["tk2"]
            mean_T_val = sum(tk[mask]) / mask_count
            push!(mean_T, mean_T_val)
            push!(max_T, maximum(tk))

            XWS = haskey(f, "XWS") ? f["XWS"] : zeros(size(tk))
            mean_XW_val = sum(XWS[mask]) / mask_count
            push!(mean_XW, mean_XW_val)
            push!(max_XW, maximum(XWS))

            qx = f["qxD"]
            qy = f["qyD"]
            qmag = sqrt.(qx.^2 .+ qy.^2)
            push!(mean_q, sum(qmag[mask]) / mask_count)

            wf = sum(phi[mask]) * dx * dy * rho_f
            push!(water_fluid, wf)

            mass_frac = @. (MH2O * XWS) / (MD + MH2O * XWS)
            ws = sum(((1.0 .- phi) .* mass_frac)[mask]) * dx * dy * rho_s
            push!(water_solid, ws)
        end
    end

    # Extract final fields from last file
    last_path = joinpath(output_dir, files[end])
    local tk_last, pf_last, XWS_last, DQPF_last, DHP_last, x_coords, y_coords
    jldopen(last_path, "r") do f
        tk_last = f["tk2"]
        pf_last = f["pf"]
        XWS_last = haskey(f, "XWS") ? f["XWS"] : zeros(size(tk_last))
        DQPF_last = haskey(f, "DQPF") ? f["DQPF"] : zeros(size(tk_last))
        DHP_last = haskey(f, "DHP") ? f["DHP"] : zeros(size(tk_last))
        x_coords = collect(f["x"])
        y_coords = collect(f["y"])
    end

    out_dict = Dict(
        "time_Ma" => time_Ma,
        "mean_T" => mean_T,
        "max_T" => max_T,
        "mean_XW" => mean_XW,
        "max_XW" => max_XW,
        "mean_q" => mean_q,
        "water_solid" => water_solid,
        "water_fluid" => water_fluid,
        "x" => x_coords,
        "y" => y_coords,
        "rplanet" => rplanet,
        "xcenter" => xcenter,
        "ycenter" => ycenter,
        "tk" => tk_last,
        "pf" => pf_last,
        "XWS" => XWS_last,
        "DQPF" => DQPF_last,
        "DHP" => DHP_last
    )

    json_path = joinpath(output_dir, "reaction_plot_data.json")
    open(json_path, "w") do io
        JSON.print(io, out_dict)
    end
    println("Exported JSON data to $json_path")
end

export_reaction_data()

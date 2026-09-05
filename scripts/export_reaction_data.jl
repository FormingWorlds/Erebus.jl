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

    local last_data

    for (i, file) in enumerate(files)
        path = joinpath(output_dir, file)
        data = load_state(path)
        last_data = data
        
        push!(time_Ma, data["timesum"] / (365.25 * 86400 * 1e6))
        
        # global water balance approximation
        dx = data["dx"]
        dy = data["dy"]
        phi = data["PHI"]
        mask = [(x - data["xcenter"])^2 + (y - data["ycenter"])^2 <= data["rplanet"]^2 for y in data["yp"], x in data["xp"]]
        
        tk = data["tk2"]
        mean_T_val = sum(tk[mask]) / sum(mask)
        push!(mean_T, mean_T_val)
        
        # Max temperature
        push!(max_T, maximum(tk))
        
        XWS = haskey(data, "XWS") ? data["XWS"] : zeros(size(tk))
        # Mean XW
        mean_XW_val = sum(XWS[mask]) / sum(mask)
        push!(mean_XW, mean_XW_val)
        push!(max_XW, maximum(XWS))
        
        # Circulation proxy (mean Darcy magnitude)
        qx = data["qxD"]
        qy = data["qyD"]
        # Very rough proxy
        qmag = sqrt.(qx.^2 .+ qy.^2)
        push!(mean_q, sum(qmag[mask]) / sum(mask))
        
        # Fluid water mass
        rho_f = 1000.0 # nominal fluid density
        wf = sum(phi[mask]) * dx * dy * rho_f
        push!(water_fluid, wf)
        
        # Solid water mass (convert molar fraction to mass fraction)
        MH2O = 0.01801528
        MD = 0.031548
        mass_frac = @. (MH2O * XWS) / (MD + MH2O * XWS)
        rho_s = 3300.0
        ws = sum(((1.0 .- phi) .* mass_frac)[mask]) * dx * dy * rho_s
        push!(water_solid, ws)
    end

    # Extract final fields
    tk = last_data["tk2"]
    pf = last_data["pf"]
    XWS = haskey(last_data, "XWS") ? last_data["XWS"] : zeros(size(tk))
    DQPF = haskey(last_data, "DQPF") ? last_data["DQPF"] : zeros(size(tk))
    DHP = haskey(last_data, "DHP") ? last_data["DHP"] : zeros(size(tk))
    
    out_dict = Dict(
        "time_Ma" => time_Ma,
        "mean_T" => mean_T,
        "max_T" => max_T,
        "mean_XW" => mean_XW,
        "max_XW" => max_XW,
        "mean_q" => mean_q,
        "water_solid" => water_solid,
        "water_fluid" => water_fluid,
        "x" => collect(last_data["x"]),
        "y" => collect(last_data["y"]),
        "rplanet" => last_data["rplanet"],
        "xcenter" => last_data["xcenter"],
        "ycenter" => last_data["ycenter"],
        "tk" => tk,
        "pf" => pf,
        "XWS" => XWS,
        "DQPF" => DQPF,
        "DHP" => DHP
    )
    
    json_path = joinpath(output_dir, "reaction_plot_data.json")
    open(json_path, "w") do io
        JSON.print(io, out_dict)
    end
    println("Exported JSON data to $json_path")
end

export_reaction_data()

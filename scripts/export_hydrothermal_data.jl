#!/usr/bin/env julia
"""
Export 2D hydrothermal benchmark simulation fields to JSON for publication plotting.
"""

using Erebus
using JLD2
using Printf

function export_data()
    output_dir = length(ARGS) >= 1 ? ARGS[1] : "output_hydrothermal"
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
    latest_file = joinpath(output_dir, files[end])
    println("Exporting from: $latest_file")
    data = load_state(latest_file)

    Nx = data["Nx"]
    Ny = data["Ny"]
    Nx1 = data["Nx1"]
    Ny1 = data["Ny1"]
    x = collect(data["x"])
    y = collect(data["y"])
    tk = data["tk2"]
    pr = data["pr"]
    pf = data["pf"]
    qx = data["qxD"]
    qy = data["qyD"]
    rx = data["RX"]

    # Effective stress
    peff = pr .- pf

    # Convert 2D arrays to nested JSON format
    json_path = joinpath(output_dir, "benchmark_plot_data.json")
    open(json_path, "w") do io
        println(io, "{")
        println(io, "  \"timestep\": $(data["timestep"]),")
        println(io, "  \"timesum_Ma\": $(data["timesum"] / (365.25 * 86400 * 1e6)),")
        println(io, "  \"rplanet\": $(data["rplanet"]),")
        println(io, "  \"rcrust\": $(data["rcrust"]),")
        println(io, "  \"xcenter\": $(data["xcenter"]),")
        println(io, "  \"ycenter\": $(data["ycenter"]),")
        println(io, "  \"x\": [", join(x, ", "), "],")
        println(io, "  \"y\": [", join(y, ", "), "],")

        # 2D tk (Ny1 x Nx1)
        println(io, "  \"tk\": [")
        for i in 1:size(tk, 1)
            row_str = join(tk[i, :], ", ")
            comma = i < size(tk, 1) ? "," : ""
            println(io, "    [", row_str, "]", comma)
        end
        println(io, "  ],")

        # 2D pr (Ny1 x Nx1)
        println(io, "  \"pr\": [")
        for i in 1:size(pr, 1)
            row_str = join(pr[i, :], ", ")
            comma = i < size(pr, 1) ? "," : ""
            println(io, "    [", row_str, "]", comma)
        end
        println(io, "  ],")

        # 2D pf (Ny1 x Nx1)
        println(io, "  \"pf\": [")
        for i in 1:size(pf, 1)
            row_str = join(pf[i, :], ", ")
            comma = i < size(pf, 1) ? "," : ""
            println(io, "    [", row_str, "]", comma)
        end
        println(io, "  ],")

        # 2D qx (Ny1 x Nx1)
        println(io, "  \"qxD\": [")
        for i in 1:size(qx, 1)
            row_str = join(qx[i, :], ", ")
            comma = i < size(qx, 1) ? "," : ""
            println(io, "    [", row_str, "]", comma)
        end
        println(io, "  ],")

        # 2D qy (Ny1 x Nx1)
        println(io, "  \"qyD\": [")
        for i in 1:size(qy, 1)
            row_str = join(qy[i, :], ", ")
            comma = i < size(qy, 1) ? "," : ""
            println(io, "    [", row_str, "]", comma)
        end
        println(io, "  ],")

        # 2D peff (Ny1 x Nx1)
        println(io, "  \"peff\": [")
        for i in 1:size(peff, 1)
            row_str = join(peff[i, :], ", ")
            comma = i < size(peff, 1) ? "," : ""
            println(io, "    [", row_str, "]", comma)
        end
        println(io, "  ]")

        return println(io, "}")
    end
    return println("Exported JSON data to $json_path")
end

export_data()

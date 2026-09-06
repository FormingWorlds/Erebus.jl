#!/usr/bin/env julia
using JLD2
using Printf

function export_frames()
    input_dir = length(ARGS) >= 1 ? ARGS[1] : "output_hydrothermal_reaction_on_128"
    out_dir = length(ARGS) >= 2 ? ARGS[2] : "/tmp/erebus_movie_frames"
    n_frames = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 150

    mkpath(out_dir)

    files = sort(
        filter(
            f -> startswith(f, "output_") && endswith(f, ".jld2") && f != "output_00000.jld2",
            readdir(input_dir)
        )
    )
    N = length(files)
    if N == 0
        error("No checkpoints found in $input_dir")
    end

    step = max(1, div(N, n_frames))
    indices = unique(vcat(collect(1:step:N), [N]))
    selected = files[indices]

    println("Exporting $(length(selected)) frames from $N checkpoints...")

    # Geometry from first file
    # Geometry from first file
    jldopen(joinpath(input_dir, selected[1]), "r") do f
        xp = collect(f["xp"])
        yp = collect(f["yp"])
        rplanet = f["rplanet"]
        xcenter = f["xcenter"]
        ycenter = f["ycenter"]

        # Save geometry metadata
        open(joinpath(out_dir, "meta.bin"), "w") do io
            write(io, Int32(length(xp)))
            write(io, Int32(length(yp)))
            write(io, Float32(rplanet))
            write(io, Float32(xcenter))
            write(io, Float32(ycenter))
            write(io, Float32.(xp))
            write(io, Float32.(yp))
        end
    end

    for (k, file) in enumerate(selected)
        path = joinpath(input_dir, file)
        jldopen(path, "r") do f
            timesum = Float32(f["timesum"])
            tk = Float32.(f["tk2"])
            XWS = haskey(f, "XWS") ? Float32.(f["XWS"]) : zeros(Float32, size(tk))
            DQPF = haskey(f, "DQPF") ? Float32.(f["DQPF"]) : zeros(Float32, size(tk))
            DHP = haskey(f, "DHP") ? Float32.(f["DHP"]) : zeros(Float32, size(tk))

            frame_file = joinpath(out_dir, @sprintf("frame_%04d.bin", k))
            open(frame_file, "w") do io
                write(io, timesum)
                write(io, tk)
                write(io, XWS)
                write(io, DQPF)
                write(io, DHP)
            end
        end
    end
    println("Done exporting $(length(selected)) frames to $out_dir")
end

export_frames()

#!/usr/bin/env julia
# tools/generate_golden_manifest.jl
# Generate golden runs for all shipped configs and write tools/golden_manifest.json.

using Pkg
const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
Pkg.activate(ROOT_DIR; io=devnull)

@eval using Dates
@eval using Erebus
@eval using JSON
@eval using LinearAlgebra
@eval using SHA
@eval include(joinpath(normpath(joinpath(@__DIR__, "..")), "tools", "golden_run.jl"))

const SHIPPED_CONFIGS = [
    "configs/test_quick.toml",
    "configs/hydrothermal_benchmark.toml",
    "configs/core_formation_benchmark.toml",
    "configs/lunar_growth_tutorial.toml",
    "configs/magma_ocean_cooling_turb_on_32.toml",
]

const MANIFEST_PATH = joinpath(ROOT_DIR, "tools", "golden_manifest.json")

function get_golden_dir()
    return get(ENV, "EREBUS_GOLDEN_DIR", joinpath(homedir(), ".erebus", "golden"))
end

function generate_manifest()
    golden_dir = get_golden_dir()
    mkpath(golden_dir)

    commit_hash = readchomp(`git -C $ROOT_DIR rev-parse HEAD`)
    host = gethostname()
    julia_ver = string(VERSION)
    blas_vendor = string(LinearAlgebra.BLAS.get_config())

    manifest = Dict{String,Any}(
        "commit" => commit_hash,
        "hostname" => host,
        "julia_version" => julia_ver,
        "blas_config" => blas_vendor,
        "generated_at" => Dates.format(now(UTC), "yyyy-mm-ddTHH:MM:SSZ"),
        "files" => Dict{String,Any}(),
    )

    for cfg_rel in SHIPPED_CONFIGS
        cfg_path = joinpath(ROOT_DIR, cfg_rel)
        cfg_name = splitext(basename(cfg_rel))[1]
        out_jld2 = joinpath(golden_dir, "$(cfg_name).jld2")

        println("Generating golden baseline for $cfg_rel -> $out_jld2...")
        run_golden(cfg_path, out_jld2)

        file_bytes = read(out_jld2)
        hash_sha256 = bytes2hex(sha256(file_bytes))

        manifest["files"][cfg_name] = Dict{String,Any}(
            "config_path" => cfg_rel,
            "filename" => "$(cfg_name).jld2",
            "sha256" => hash_sha256,
            "size_bytes" => length(file_bytes),
        )
    end

    open(MANIFEST_PATH, "w") do io
        JSON.print(io, manifest, 4)
        return println(io)
    end
    println("Successfully generated golden manifest: $MANIFEST_PATH")
    return manifest
end

function main()
    generate_manifest()
    return exit(0)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

using Test
using Erebus
using Erebus.Config

include("test_helpers.jl")

function extract_reference_diagnostics(state, cfg)
    dx = cfg.grid.xsize / (cfg.grid.Nx - 1)
    dy = cfg.grid.ysize / (cfg.grid.Ny - 1)
    dV = dx * dy

    RHO = state.grids.RHO
    RHOCP = state.grids.RHOCP
    tk2 = state.grids.tk2

    total_mass = sum(RHO) * dV
    total_energy = sum(RHOCP .* tk2) * dV
    max_temperature = maximum(tk2)
    core_radius = Float64(state.rcore)
    degassed_mass = Float64(state.M_vent_total)

    return Dict(
        "total_mass" => total_mass,
        "total_energy" => total_energy,
        "max_temperature" => max_temperature,
        "core_radius" => core_radius,
        "degassed_mass" => degassed_mass,
    )
end

function load_reference_json(path::String)
    data = Dict{String,Dict{String,Float64}}()
    current_run = ""
    for line in eachline(path)
        trimmed = strip(line)
        if startswith(trimmed, "\"") && endswith(trimmed, "{")
            m = match(r"\"([^\"]+)\"\s*:\s*\{", trimmed)
            if m !== nothing
                current_run = m.captures[1]
                data[current_run] = Dict{String,Float64}()
            end
        elseif startswith(trimmed, "\"") && occursin(":", trimmed) && !isempty(current_run)
            m = match(r"\"([^\"]+)\"\s*:\s*([0-9.eE+-]+)", trimmed)
            if m !== nothing
                key = m.captures[1]
                val = parse(Float64, m.captures[2])
                data[current_run][key] = val
            end
        end
    end
    return data
end

@testset "Physical Reference Runs Regression Baselines" begin
    ref_path = joinpath(@__DIR__, "data", "reference_runs.json")
    @test isfile(ref_path)
    ref_data = load_reference_json(ref_path)
    @test haskey(ref_data, "test_quick")

    @testset "test_quick.toml (2 timesteps)" begin
        quick_cfg_path = joinpath(@__DIR__, "..", "configs", "test_quick.toml")
        cfg = load_config(quick_cfg_path)
        mktempdir() do outdir
            state = Erebus.simulation_loop(cfg; output_path=outdir)
            diag = extract_reference_diagnostics(state, cfg)
            expected = ref_data["test_quick"]

            @test isapprox(diag["total_mass"], expected["total_mass"]; rtol=1e-3)
            @test isapprox(diag["total_energy"], expected["total_energy"]; rtol=1e-3)
            @test isapprox(diag["max_temperature"], expected["max_temperature"]; rtol=1e-3)
            @test isapprox(
                diag["core_radius"], expected["core_radius"]; rtol=1e-3, atol=1e-6
            )
            @test isapprox(
                diag["degassed_mass"], expected["degassed_mass"]; rtol=1e-3, atol=1e-6
            )
        end
    end

    @testset "magma_ocean_cooling_turb_on_32.toml (3 timesteps)" begin
        mo_cfg_path = joinpath(
            @__DIR__, "..", "configs", "magma_ocean_cooling_turb_on_32.toml"
        )
        cfg = load_config(mo_cfg_path)
        cfg = override_config(cfg, Dict("time.n_steps" => 3))
        mktempdir() do outdir
            state = Erebus.simulation_loop(cfg; output_path=outdir)
            diag = extract_reference_diagnostics(state, cfg)
            expected = ref_data["magma_ocean_cooling_turb_on_32"]

            @test isapprox(diag["total_mass"], expected["total_mass"]; rtol=1e-3)
            @test isapprox(diag["total_energy"], expected["total_energy"]; rtol=1e-3)
            @test isapprox(diag["max_temperature"], expected["max_temperature"]; rtol=1e-3)
            @test isapprox(
                diag["core_radius"], expected["core_radius"]; rtol=1e-3, atol=1e-6
            )
            @test isapprox(
                diag["degassed_mass"], expected["degassed_mass"]; rtol=1e-3, atol=1e-6
            )
        end
    end

    @testset "core_formation_benchmark.toml (3 timesteps)" begin
        cf_cfg_path = joinpath(@__DIR__, "..", "configs", "core_formation_benchmark.toml")
        cfg = load_config(cf_cfg_path)
        cfg = override_config(cfg, Dict("time.n_steps" => 3))
        mktempdir() do outdir
            state = Erebus.simulation_loop(cfg; output_path=outdir)
            diag = extract_reference_diagnostics(state, cfg)
            expected = ref_data["core_formation_benchmark"]

            @test isapprox(diag["total_mass"], expected["total_mass"]; rtol=1e-3)
            @test isapprox(diag["total_energy"], expected["total_energy"]; rtol=1e-3)
            @test isapprox(diag["max_temperature"], expected["max_temperature"]; rtol=1e-3)
            @test isapprox(
                diag["core_radius"], expected["core_radius"]; rtol=1e-3, atol=1e-6
            )
            @test isapprox(
                diag["degassed_mass"], expected["degassed_mass"]; rtol=1e-3, atol=1e-6
            )
        end
    end
end

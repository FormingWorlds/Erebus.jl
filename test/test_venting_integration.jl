using Erebus
using ExtendableSparse
using JLD2
using LinearSolve
using Random
using StaticArrays
using Test

include("../src/test_constants.jl")
const rgen = MersenneTwister(seed)

@testset "Cold Surface Venting Integration" begin
    @testset "Runtime loop with venting inactive (baseline)" begin
        output_dir = mktempdir()
        try
            quick_toml = joinpath(@__DIR__, "..", "configs", "test_quick.toml")
            cfg = load_config(quick_toml)
            # 2 steps with venting inactive
            cfg_off = SimulationConfig(
                grid=cfg.grid,
                geometry=cfg.geometry,
                time=TimeConfig(
                    dt_initial=cfg.time.dt_initial,
                    dt_longest=cfg.time.dt_longest,
                    dtcoefdn=cfg.time.dtcoefdn,
                    dtcoefup=cfg.time.dtcoefup,
                    dtstep=cfg.time.dtstep,
                    dxymax=cfg.time.dxymax,
                    vpratio=cfg.time.vpratio,
                    DTmax=cfg.time.DTmax,
                    start_time=cfg.time.start_time,
                    endtime=cfg.time.endtime,
                    start_step=1,
                    n_steps=2,
                ),
                solver=cfg.solver,
                poroelasticity=cfg.poroelasticity,
                thermodynamics=cfg.thermodynamics,
                reaction=cfg.reaction,
                materials=cfg.materials,
                output=OutputConfig(output_dir=output_dir, savematstep=2),
                disk=cfg.disk,
                melting=cfg.melting,
                venting=VentingConfig(active=false),
            )
            Erebus.simulation_loop(cfg_off; output_path=output_dir)

            files = readdir(output_dir)
            @test "output_00000.jld2" in files
            @test "output_00002.jld2" in files

            data2 = load_state(joinpath(output_dir, "output_00002.jld2"))
            @test data2["timestep"] == 2
            @test haskey(data2, "M_vent_total")
            @test isapprox(data2["M_vent_total"], 0.0; atol=1e-12)
            @test !any(isnan, data2["tk2"])
            @test !any(isinf, data2["tk2"])
            @test all(data2["tk2"] .> 0.0)
        finally
            rm(output_dir; recursive=true, force=true)
        end
    end

    @testset "Runtime loop with cold surface venting active" begin
        output_dir = mktempdir()
        try
            quick_toml = joinpath(@__DIR__, "..", "configs", "test_quick.toml")
            cfg = load_config(quick_toml)
            # 2 steps with venting active, latent cooling enabled, darcy sink mode
            cfg_on = SimulationConfig(
                grid=cfg.grid,
                geometry=cfg.geometry,
                time=TimeConfig(
                    dt_initial=cfg.time.dt_initial,
                    dt_longest=cfg.time.dt_longest,
                    dtcoefdn=cfg.time.dtcoefdn,
                    dtcoefup=cfg.time.dtcoefup,
                    dtstep=cfg.time.dtstep,
                    dxymax=cfg.time.dxymax,
                    vpratio=cfg.time.vpratio,
                    DTmax=cfg.time.DTmax,
                    start_time=cfg.time.start_time,
                    endtime=cfg.time.endtime,
                    start_step=1,
                    n_steps=2,
                ),
                solver=cfg.solver,
                poroelasticity=cfg.poroelasticity,
                thermodynamics=cfg.thermodynamics,
                reaction=cfg.reaction,
                materials=cfg.materials,
                output=OutputConfig(output_dir=output_dir, savematstep=2),
                disk=DiskConfig(
                    t_dispersal_myr=0.01,
                    dt_dispersal_myr=0.002,
                    p_amb_disk=10.0,
                    p_amb_space=1.0e-4,
                    dispersal_active=true,
                ),
                melting=cfg.melting,
                venting=VentingConfig(
                    active=true,
                    mode=:darcy_sink,
                    k_vent=1.0e-11,
                    conductance_factor=1.0,
                    L_sublimation=2.83e6,
                    latent_cooling=true,
                ),
            )
            Erebus.simulation_loop(cfg_on; output_path=output_dir)

            files = readdir(output_dir)
            @test "output_00000.jld2" in files
            @test "output_00002.jld2" in files

            data2 = load_state(joinpath(output_dir, "output_00002.jld2"))
            @test data2["timestep"] == 2
            @test haskey(data2, "M_vent_total")
            @test data2["M_vent_total"] > 0.0
            @test haskey(data2, "P_amb")
            @test data2["P_amb"] > 0.0
            @test haskey(data2, "S_vent")
            @test size(data2["S_vent"]) == (cfg.grid.Ny + 1, cfg.grid.Nx + 1)
            @test all(data2["S_vent"] .>= 0.0)
            @test any(data2["S_vent"] .> 0.0)

            # Solution field health checks
            @test !any(isnan, data2["tk2"])
            @test !any(isinf, data2["tk2"])
            @test all(data2["tk2"] .> 0.0)
            @test !any(isnan, data2["pf"])
            @test !any(isinf, data2["pf"])
            @test !any(isnan, data2["pr"])
            @test !any(isinf, data2["pr"])

            # Marker porosity bounds verification
            @test haskey(data2, "phim")
            @test all(data2["phim"] .>= cfg.poroelasticity.phimin)
            @test all(data2["phim"] .<= cfg.poroelasticity.phimax)
        finally
            rm(output_dir; recursive=true, force=true)
        end
    end
end

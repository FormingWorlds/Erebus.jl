@testset "Integration" begin
    @testset "simulation_loop() execution and state persistence" begin
        output_dir = mktempdir()
        try
            # Execute simulation loop on test grid
            Erebus.simulation_loop(output_dir)

            # Verify output files exist
            files = readdir(output_dir)
            @test "output_00000.jld2" in files
            @test "output_00010.jld2" in files

            # Verify initial state file contents
            data0 = JLD2.load(joinpath(output_dir, "output_00000.jld2"))
            @test haskey(data0, "timestep")
            @test data0["timestep"] == 0
            @test haskey(data0, "timesum")
            @test haskey(data0, "marknum")
            @test data0["marknum"] > 0

            # Verify completed state file contents
            data10 = JLD2.load(joinpath(output_dir, "output_00010.jld2"))
            @test haskey(data10, "timestep")
            @test data10["timestep"] == 10
            @test data10["timesum"] > data0["timesum"]
            @test haskey(data10, "tk2")
            @test !any(isnan, data10["tk2"])
            @test !any(isinf, data10["tk2"])
            @test all(data10["tk2"] .> 0.0)
            @test haskey(data10, "tkm")
            @test all(data10["tkm"] .> 0.0)
        finally
            rm(output_dir, recursive=true, force=true)
        end
    end

    @testset "simulation_loop() execution with TOML config" begin
        output_dir = mktempdir()
        try
            quick_toml = joinpath(@__DIR__, "..", "configs", "test_quick.toml")
            cfg = load_config(quick_toml)
            Erebus.simulation_loop(cfg; output_path=output_dir)
            files = readdir(output_dir)
            @test "output_00000.jld2" in files
            @test "output_00002.jld2" in files
            data2 = JLD2.load(joinpath(output_dir, "output_00002.jld2"))
            @test data2["timestep"] == 2
            @test !any(isnan, data2["tk2"])
        finally
            rm(output_dir, recursive=true, force=true)
        end
    end

    @testset "2D hydrothermal circulation benchmark execution on 32x32 grid" begin
        output_dir = mktempdir()
        try
            bench_toml = joinpath(@__DIR__, "..", "configs", "hydrothermal_benchmark.toml")
            cfg = load_config(bench_toml)
            # Run 2 steps for integration testing
            cfg_test = SimulationConfig(
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
                materials=cfg.materials,
                output=OutputConfig(output_dir=output_dir, savematstep=2),
            )
            Erebus.simulation_loop(cfg_test; output_path=output_dir)

            files = readdir(output_dir)
            @test "output_00000.jld2" in files
            @test "output_00002.jld2" in files

            data2 = load_state(joinpath(output_dir, "output_00002.jld2"))
            @test data2["timestep"] == 2
            @test !any(isnan, data2["tk2"])
            @test !any(isinf, data2["tk2"])
            @test all(data2["tk2"] .> 0.0)

            # Check Darcy flux fields are finite
            @test !any(isnan, data2["qxD"])
            @test !any(isnan, data2["qyD"])
            @test !any(isinf, data2["qxD"])
            @test !any(isinf, data2["qyD"])

            # Check thermal conductivity fields are strictly positive on physical interior
            @test all(data2["KX"][2:(end - 1), 2:(end - 1)] .> 0.0)
            @test all(data2["KY"][2:(end - 1), 2:(end - 1)] .> 0.0)

            # Check fluid drag ratio fields are strictly positive on physical interior
            @test all(data2["RX"][2:(end - 1), 2:(end - 1)] .> 0.0)
            @test all(data2["RY"][2:(end - 1), 2:(end - 1)] .> 0.0)
        finally
            rm(output_dir, recursive=true, force=true)
        end
    end

    @testset "checkpoint restart continuity" begin
        dir_orig = mktempdir()
        dir_resume = mktempdir()
        try
            quick_toml = joinpath(@__DIR__, "..", "configs", "test_quick.toml")
            cfg = load_config(quick_toml)

            # 1. Run 3 steps continuously in dir_orig, saving every step
            cfg_orig = SimulationConfig(
                grid=cfg.grid,
                geometry=cfg.geometry,
                time=TimeConfig(
                    dt_initial=cfg.time.dt_initial,
                    dt_longest=cfg.time.dt_longest,
                    start_step=1,
                    n_steps=3,
                ),
                solver=cfg.solver,
                poroelasticity=cfg.poroelasticity,
                thermodynamics=cfg.thermodynamics,
                materials=cfg.materials,
                output=OutputConfig(output_dir=dir_orig, savematstep=1),
            )
            Erebus.simulation_loop(cfg_orig; output_path=dir_orig)

            # 2. Resume from step 2 checkpoint in dir_resume and run to step 3
            ckpt2_path = joinpath(dir_orig, "output_00002.jld2")
            @test isfile(ckpt2_path)

            cfg_resume = SimulationConfig(
                grid=cfg.grid,
                geometry=cfg.geometry,
                time=TimeConfig(
                    dt_initial=cfg.time.dt_initial,
                    dt_longest=cfg.time.dt_longest,
                    start_step=1,
                    n_steps=3,
                ),
                solver=cfg.solver,
                poroelasticity=cfg.poroelasticity,
                thermodynamics=cfg.thermodynamics,
                materials=cfg.materials,
                output=OutputConfig(
                    output_dir=dir_resume, savematstep=1, restart_from=ckpt2_path
                ),
            )
            Erebus.simulation_loop(
                cfg_resume; output_path=dir_resume, restart_from=ckpt2_path
            )

            @test isfile(joinpath(dir_resume, "output_00003.jld2"))
            data_orig = load_state(joinpath(dir_orig, "output_00003.jld2"))
            data_res = load_state(joinpath(dir_resume, "output_00003.jld2"))

            # Verify identical state progression across checkpoint boundary
            @test data_res["timestep"] == 3
            @test data_res["timesum"] ≈ data_orig["timesum"] rtol=1e-12
            @test isapprox(data_res["tk2"], data_orig["tk2"]; rtol=1e-8, atol=1e-8)
            @test isapprox(data_res["pr"], data_orig["pr"]; rtol=1e-8, atol=1e-8)
            @test isapprox(data_res["pf"], data_orig["pf"]; rtol=1e-8, atol=1e-8)
            @test isapprox(data_res["vx"], data_orig["vx"]; rtol=1e-8, atol=1e-14)
            @test isapprox(data_res["vy"], data_orig["vy"]; rtol=1e-8, atol=1e-14)
            @test isapprox(data_res["qxD"], data_orig["qxD"]; rtol=1e-8, atol=1e-14)
            @test isapprox(data_res["qyD"], data_orig["qyD"]; rtol=1e-8, atol=1e-14)
        finally
            rm(dir_orig, recursive=true, force=true)
            rm(dir_resume, recursive=true, force=true)
        end
    end
end

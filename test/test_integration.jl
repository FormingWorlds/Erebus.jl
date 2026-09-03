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
            Erebus.simulation_loop(cfg; output_path = output_dir)
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
end

using Test
using Erebus
using JLD2

@testset "Streaming Telemetry & Output Modes" begin
    @testset "Telemetry File Initialization & Streaming" begin
        mktempdir() do tmpdir
            io = init_telemetry(tmpdir, "test_telemetry.csv")
            @test isopen(io)

            stream_telemetry_row!(
                io,
                1,
                0.01,
                100.0,
                50000.0,
                15000.0,
                1200.0,
                800.0,
                0.25,
                0.10,
                1.5e15,
                1.2e15,
                2.0e14,
                1.0e13,
                0.40,
                0.05,
            )
            close(io)

            csv_path = joinpath(tmpdir, "test_telemetry.csv")
            @test isfile(csv_path)

            lines = readlines(csv_path)
            @test length(lines) == 2

            header_cols = split(lines[1], ",")
            @test length(header_cols) == 16
            @test header_cols[1] == "step"
            @test header_cols[2] == "time_Ma"
            @test header_cols[10] == "M_outgassed_total"
            @test header_cols[15] == "F_melt_mean"
            @test header_cols[16] == "dt_aphimax_max"

            row_cols = split(lines[2], ",")
            @test parse(Int, row_cols[1]) == 1
            @test isapprox(parse(Float64, row_cols[2]), 0.01; atol=1e-6)
            @test isapprox(parse(Float64, row_cols[4]), 50000.0; atol=1e-3)
            @test isapprox(parse(Float64, row_cols[10]), 1.5e15; rtol=1e-6)

            # Test append=true on existing file
            io_append = init_telemetry(tmpdir, "test_telemetry.csv"; append=true)
            stream_telemetry_row!(
                io_append,
                2,
                0.02,
                100.0,
                50000.0,
                15000.0,
                400.0,
                320.0,
                0.24,
                0.09,
                1.6e15,
                1.3e15,
                2.1e14,
                1.1e13,
                0.38,
                0.04,
            )
            close(io_append)
            lines_app = readlines(csv_path)
            @test length(lines_app) == 3  # Header + 2 data rows (no duplicate header)
            @test split(lines_app[3], ",")[1] == "2"
        end
    end

    @testset "Output Configuration Validation" begin
        cfg_def = default_config()
        @test cfg_def.output.mode == :snapshots
        @test cfg_def.output.telemetrystep == 1
        @test cfg_def.output.telemetry_file == "telemetry.csv"
        @test cfg_def.output.save_final == true

        # Rejection of invalid output modes
        @test_throws ArgumentError validate_config(
            override_config(cfg_def, Dict("output.mode" => "invalid_mode"))
        )

        # Rejection of invalid telemetry step (< 1)
        @test_throws ArgumentError validate_config(
            override_config(cfg_def, Dict("output.telemetrystep" => 0))
        )

        # Rejection of empty telemetry file
        @test_throws ArgumentError validate_config(
            override_config(cfg_def, Dict("output.telemetry_file" => ""))
        )
    end

    @testset "Simulation Output Mode Integration" begin
        mktempdir() do tmpdir
            # 1. Test :both mode (snapshots + telemetry)
            out_both = joinpath(tmpdir, "both")
            cfg_both = override_config(
                load_config("configs/test_quick.toml"),
                Dict(
                    "output.output_dir" => out_both,
                    "output.mode" => "both",
                    "output.savematstep" => 1,
                    "output.telemetrystep" => 1,
                ),
            )
            simulation_loop(cfg_both)

            telemetry_path = joinpath(out_both, "telemetry.csv")
            @test isfile(telemetry_path)
            lines = readlines(telemetry_path)
            @test length(lines) >= 3  # Header + at least 2 steps

            # Checkpoints must exist
            @test isfile(joinpath(out_both, "output_00000.jld2"))
            @test isfile(joinpath(out_both, "output_00001.jld2"))
            @test isfile(joinpath(out_both, "output_00002.jld2"))

            # 2. Test :telemetry mode with save_final=true
            out_telem = joinpath(tmpdir, "telem")
            cfg_telem = override_config(
                load_config("configs/test_quick.toml"),
                Dict(
                    "output.output_dir" => out_telem,
                    "output.mode" => "telemetry",
                    "output.savematstep" => 1,
                    "output.telemetrystep" => 1,
                    "output.save_final" => true,
                ),
            )
            simulation_loop(cfg_telem)

            @test isfile(joinpath(out_telem, "telemetry.csv"))
            # Step 0 and 1 should be skipped because mode is :telemetry
            @test !isfile(joinpath(out_telem, "output_00000.jld2"))
            @test !isfile(joinpath(out_telem, "output_00001.jld2"))
            # Final step 2 should be saved because save_final is true
            @test isfile(joinpath(out_telem, "output_00002.jld2"))
        end
    end
end

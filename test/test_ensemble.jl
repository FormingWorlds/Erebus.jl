using Test
using Erebus

@testset "Ensemble Parameter Sweeps" begin
    @testset "Configuration Override Mechanism" begin
        cfg_base = default_config()
        overrides = Dict(
            "thermodynamics.phim0" => 0.35,
            "geometry.rplanet" => 55000.0,
            "geometry.rcrust" => 52000.0,
            "solver.seed" => 999,
            "output.mode" => "telemetry",
        )
        cfg_mod = override_config(cfg_base, overrides)

        @test isapprox(cfg_mod.thermodynamics.phim0, 0.35; atol=1e-6)
        @test isapprox(cfg_mod.geometry.rplanet, 55000.0; atol=1e-3)
        @test isapprox(cfg_mod.geometry.rcrust, 52000.0; atol=1e-3)
        @test cfg_mod.solver.seed == 999
        @test cfg_mod.output.mode == :telemetry

        # Malformed key without dot separator
        @test_throws ArgumentError override_config(cfg_base, Dict("invalidkey" => 123))

        # Unknown section
        @test_throws ArgumentError override_config(
            cfg_base, Dict("unknown_section.foo" => 123)
        )
    end

    @testset "Ensemble Specification Validation" begin
        cfg_base = default_config()

        spec_valid = EnsembleSweepSpec(
            cfg_base;
            output_dir="test_ens",
            sampling_method=:lhs,
            num_samples=5,
            parameters=Dict("thermodynamics.phim0" => [0.1, 0.4]),
        )
        @test spec_valid.num_samples == 5
        @test spec_valid.sampling_method == :lhs

        # Invalid sampling method
        @test_throws ArgumentError EnsembleSweepSpec(cfg_base; sampling_method=:bogus)

        # Invalid num_samples (< 1)
        @test_throws ArgumentError EnsembleSweepSpec(cfg_base; num_samples=0)

        # Positional constructor validation
        @test_throws ArgumentError EnsembleSweepSpec(
            cfg_base, "test_ens", :invalid_method, 5, Dict{String,Any}(), 42
        )
        @test_throws ArgumentError EnsembleSweepSpec(
            cfg_base, "test_ens", :lhs, 0, Dict{String,Any}(), 42
        )
    end

    @testset "Parameter Sampling Algorithms" begin
        cfg_base = default_config()

        # 1. Grid sampling (Cartesian product)
        spec_grid = EnsembleSweepSpec(
            cfg_base;
            sampling_method=:grid,
            parameters=Dict(
                "thermodynamics.phim0" => [0.2, 0.35], "solver.tile_size" => [4, 8]
            ),
        )
        runs_grid = sample_parameters(spec_grid)
        @test length(runs_grid) == 4
        all_phim0 = [r[2]["thermodynamics.phim0"] for r in runs_grid]
        @test count(x -> isapprox(x, 0.2; atol=1e-6), all_phim0) == 2
        @test count(x -> isapprox(x, 0.35; atol=1e-6), all_phim0) == 2

        # 2. Latin Hypercube Sampling
        N_lhs = 8
        lo, hi = 0.1, 0.5
        spec_lhs = EnsembleSweepSpec(
            cfg_base;
            sampling_method=:lhs,
            num_samples=N_lhs,
            parameters=Dict("thermodynamics.phim0" => [lo, hi]),
            seed=42,
        )
        runs_lhs = sample_parameters(spec_lhs)
        @test length(runs_lhs) == N_lhs

        vals_lhs = [r[2]["thermodynamics.phim0"] for r in runs_lhs]
        @test all(v -> lo <= v <= hi, vals_lhs)

        # Verify stratification: exactly 1 point in each of the N bins
        bin_width = (hi - lo) / N_lhs
        bin_counts = zeros(Int, N_lhs)
        for v in vals_lhs
            bin_idx = clamp(floor(Int, (v - lo) / bin_width) + 1, 1, N_lhs)
            bin_counts[bin_idx] += 1
        end
        @test all(c -> c == 1, bin_counts)

        # 3. Random sampling
        spec_rand = EnsembleSweepSpec(
            cfg_base;
            sampling_method=:random,
            num_samples=6,
            parameters=Dict("geometry.rplanet" => [50000.0, 60000.0]),
        )
        runs_rand = sample_parameters(spec_rand)
        @test length(runs_rand) == 6
        vals_rand = [r[2]["geometry.rplanet"] for r in runs_rand]
        @test all(v -> 50000.0 <= v <= 60000.0, vals_rand)
    end

    @testset "Ensemble Execution & Catalog Output" begin
        mktempdir() do tmpdir
            base_cfg = load_config("configs/test_quick.toml")
            # Set to telemetry mode for compact storage
            base_cfg_telem = override_config(
                base_cfg,
                Dict(
                    "output.mode" => "telemetry",
                    "output.telemetrystep" => 1,
                    "time.n_steps" => 1,
                ),
            )

            ens_dir = joinpath(tmpdir, "sweep")
            spec = EnsembleSweepSpec(
                base_cfg_telem;
                output_dir=ens_dir,
                sampling_method=:lhs,
                num_samples=2,
                parameters=Dict("thermodynamics.phim0" => [0.15, 0.35]),
                seed=100,
            )

            catalog = run_ensemble(spec; max_workers=2, verbose=false)
            @test length(catalog) == 2
            @test catalog[1]["status"] == "success"
            @test catalog[2]["status"] == "success"
            @test haskey(catalog[1], "walltime_s")

            catalog_csv = joinpath(ens_dir, "catalog.csv")
            @test isfile(catalog_csv)
            lines = readlines(catalog_csv)
            @test length(lines) == 3  # Header + 2 run rows
            @test occursin("run_id", lines[1])
            @test occursin("status", lines[1])
            @test occursin("thermodynamics.phim0", lines[1])
        end
    end

    @testset "Catalog CSV Escaping" begin
        mktempdir() do tmpdir
            catalog = [
                Dict{String,Any}(
                    "run_id" => "run_0001",
                    "status" => "failed",
                    "walltime_s" => 1.23,
                    "error_msg" => "ArgumentError: parameter must be :tiled or :buffered, got :foo",
                    "description" => "Test \"quoted\" value\nwith newline",
                ),
            ]
            csv_path = joinpath(tmpdir, "catalog.csv")
            save_ensemble_catalog(catalog, csv_path)
            @test isfile(csv_path)
            content = read(csv_path, String)
            @test occursin(
                "\"ArgumentError: parameter must be :tiled or :buffered, got :foo\"",
                content,
            )
            @test occursin("\"Test \"\"quoted\"\" value\nwith newline\"", content)
        end
    end
end

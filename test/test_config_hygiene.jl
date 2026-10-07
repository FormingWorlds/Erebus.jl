using Test
using Erebus
using TOML
using StaticArrays

@testset "Configuration Hygiene and Shipped Configurations" begin
    @testset "Removed Dead Config Fields Rejection" begin
        # 1. MeltingConfig dead fields
        toml_melt = """
        [melting]
        latent_heat_mode = "apparent_cp"
        """
        @test_throws ArgumentError load_config(toml_melt)
        err_melt = try
            load_config(toml_melt)
        catch e
            e
        end
        @test occursin("latent_heat_mode", err_melt.msg)
        @test occursin("removed", err_melt.msg)

        # 2. VolatilesConfig dead fields
        dead_volatiles = [
            ("h2_active", "true"),
            ("h2_law", "\"gaillard2003\""),
            ("t_organic_devol", "550.0"),
            ("dt_organic_devol", "50.0"),
            ("organic_n_initial_ppm", "500.0"),
            ("melt_feo_wtpct", "10.0"),
            ("x_sio2", "0.56"),
            ("x_al2o3", "0.11"),
            ("x_tio2", "0.01"),
        ]
        for (fld, val_str) in dead_volatiles
            toml_vol = """
            [volatiles]
            $fld = $val_str
            """
            @test_throws ArgumentError load_config(toml_vol)
            err_vol = try
                load_config(toml_vol)
            catch e
                e
            end
            @test occursin(fld, err_vol.msg)
            @test occursin("removed", err_vol.msg)
        end

        # 3. RedoxConfig dead fields
        toml_rdx = """
        [redox]
        venting_redox = true
        """
        @test_throws ArgumentError load_config(toml_rdx)
        err_rdx = try
            load_config(toml_rdx)
        catch e
            e
        end
        @test occursin("venting_redox", err_rdx.msg)
        @test occursin("removed", err_rdx.msg)

        # 4. AtmosphereConfig dead fields
        dead_atm = [("kappa_vis_default", "1.0e-3"), ("b_diff_ref", "1.0e21")]
        for (fld, val_str) in dead_atm
            toml_atm = """
            [atmosphere]
            $fld = $val_str
            """
            @test_throws ArgumentError load_config(toml_atm)
            err_atm = try
                load_config(toml_atm)
            catch e
                e
            end
            @test occursin(fld, err_atm.msg)
            @test occursin("removed", err_atm.msg)
        end

        # 5. TelescopingConfig dead fields
        toml_tele = """
        [telescoping]
        target_radius = 1737000.0
        """
        @test_throws ArgumentError load_config(toml_tele)
        err_tele = try
            load_config(toml_tele)
        catch e
            e
        end
        @test occursin("target_radius", err_tele.msg)
        @test occursin("removed", err_tele.msg)

        # 6. MagmaOceanDegassingConfig removed field
        toml_degas = """
        [magma_degassing]
        crystallization_degassing = true
        """
        @test_throws ArgumentError load_config(toml_degas)
        err_degas = try
            load_config(toml_degas)
        catch e
            e
        end
        @test occursin("crystallization_degassing", err_degas.msg)
    end

    @testset "Strict Validation Bounds: TimeConfig" begin
        cfg_base = default_config()

        # dtcoefdn bounds: (0, 1)
        @test_throws ArgumentError validate_config(
            SimulationConfig(time=TimeConfig(dtcoefdn=0.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(time=TimeConfig(dtcoefdn=1.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(time=TimeConfig(dtcoefdn=-0.2))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(time=TimeConfig(dtcoefdn=NaN))
        )

        # dtcoefup bounds: [1, 10]
        @test_throws ArgumentError validate_config(
            SimulationConfig(time=TimeConfig(dtcoefup=0.9))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(time=TimeConfig(dtcoefup=12.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(time=TimeConfig(dtcoefup=NaN))
        )

        # dtstep bounds: >= 1
        @test_throws ArgumentError validate_config(
            SimulationConfig(time=TimeConfig(dtstep=0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(time=TimeConfig(dtstep=-5))
        )

        # dxymax bounds: (0, 1]
        @test_throws ArgumentError validate_config(
            SimulationConfig(time=TimeConfig(dxymax=0.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(time=TimeConfig(dxymax=1.5))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(time=TimeConfig(dxymax=NaN))
        )

        # DTmax bounds: > 0 and finite
        @test_throws ArgumentError validate_config(
            SimulationConfig(time=TimeConfig(DTmax=0.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(time=TimeConfig(DTmax=-10.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(time=TimeConfig(DTmax=Inf))
        )
    end

    @testset "Strict Validation Bounds: SolverConfig" begin
        # yerrmax bounds: (0, 1e6]
        @test_throws ArgumentError validate_config(
            SimulationConfig(solver=SolverConfig(yerrmax=0.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(solver=SolverConfig(yerrmax=1.0e8))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(solver=SolverConfig(yerrmax=NaN))
        )

        # dphimax bounds: (0, 1]
        @test_throws ArgumentError validate_config(
            SimulationConfig(solver=SolverConfig(dphimax=0.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(solver=SolverConfig(dphimax=100.01))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(solver=SolverConfig(dphimax=NaN))
        )

        # etawt bounds: [0, 1)
        @test_throws ArgumentError validate_config(
            SimulationConfig(solver=SolverConfig(etawt=-0.1))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(solver=SolverConfig(etawt=1.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(solver=SolverConfig(etawt=NaN))
        )

        # Valid boundary values succeed
        cfg_valid_edge = SimulationConfig(
            time=TimeConfig(dtcoefdn=0.01, dtcoefup=1.0, dxymax=1.0, dtstep=1, DTmax=5.0),
            solver=SolverConfig(yerrmax=1.0e6, dphimax=1.0, etawt=0.0),
        )
        @test validate_config(cfg_valid_edge) === nothing
        @test isapprox(cfg_valid_edge.time.dxymax, 1.0; atol=1e-12)
        @test isapprox(cfg_valid_edge.solver.dphimax, 1.0; atol=1e-12)
    end

    @testset "Ensemble Loader and Dispatch Error Handling" begin
        # Passing an ensemble file to load_config raises dedicated ArgumentError
        toml_ensemble = """
        [base]
        config = "test_quick.toml"

        [ensemble]
        num_samples = 2
        """
        @test_throws ArgumentError load_config(toml_ensemble)
        err = try
            load_config(toml_ensemble)
        catch e
            e
        end
        @test occursin("ensemble section", err.msg)
        @test occursin("load_ensemble_config", err.msg)

        # load_ensemble_config loads and parses successfully
        spec_path = normpath(
            joinpath(@__DIR__, "..", "configs", "test_ensemble_sweep.toml")
        )
        spec = load_ensemble_config(spec_path)
        @test spec isa EnsembleSweepSpec
        @test spec.num_samples == 3
        @test spec.sampling_method === :lhs
        @test spec.seed == 42
        @test isapprox(spec.base_config.solver.dphimax, 0.1; atol=1e-6)
    end

    @testset "Ensemble Parameter Sampling: Integer Handling and Seed Preservation" begin
        cfg_base = default_config()

        # 1. Latin Hypercube sampling with integer parameters
        spec_int_lhs = EnsembleSweepSpec(
            cfg_base;
            output_dir=mktempdir(),
            sampling_method=:lhs,
            num_samples=5,
            parameters=Dict(
                "grid.Nx" => [33, 65],
                "time.n_steps" => [5, 20],
                "thermodynamics.phim0" => [0.15, 0.35],
            ),
            seed=123,
        )
        runs_lhs = sample_parameters(spec_int_lhs)
        @test length(runs_lhs) == 5
        for (run_id, params, run_cfg) in runs_lhs
            @test params["grid.Nx"] isa Integer
            @test params["time.n_steps"] isa Integer
            @test params["thermodynamics.phim0"] isa Float64
            @test 33 <= run_cfg.grid.Nx <= 65
            @test 5 <= run_cfg.time.n_steps <= 20
            @test 0.15 <= run_cfg.thermodynamics.phim0 <= 0.35
            @test run_cfg.grid.Nx == params["grid.Nx"]
        end

        # 2. Random sampling with integer parameters
        spec_int_rand = EnsembleSweepSpec(
            cfg_base;
            output_dir=mktempdir(),
            sampling_method=:random,
            num_samples=5,
            parameters=Dict("grid.Ny" => [33, 65], "time.dtstep" => [50, 200]),
            seed=456,
        )
        runs_rand = sample_parameters(spec_int_rand)
        @test length(runs_rand) == 5
        for (run_id, params, run_cfg) in runs_rand
            @test params["grid.Ny"] isa Integer
            @test params["time.dtstep"] isa Integer
            @test 33 <= run_cfg.grid.Ny <= 65
            @test 50 <= run_cfg.time.dtstep <= 200
        end

        # 3. Swept seed preservation
        spec_seed = EnsembleSweepSpec(
            cfg_base;
            output_dir=mktempdir(),
            sampling_method=:grid,
            parameters=Dict("solver.seed" => [777, 888]),
            seed=42,
        )
        runs_seed = sample_parameters(spec_seed)
        @test length(runs_seed) == 2
        @test runs_seed[1][3].solver.seed == 777
        @test runs_seed[2][3].solver.seed == 888
        # 4. Empty parameters sampling
        spec_empty = EnsembleSweepSpec(
            cfg_base;
            output_dir=mktempdir(),
            sampling_method=:random,
            num_samples=3,
            parameters=Dict{String,Vector{Any}}(),
            seed=42,
        )
        runs_empty = sample_parameters(spec_empty)
        @test length(runs_empty) == 3
        @test isempty(runs_empty[1][2])

        # 5. Invalid sampling method throws ArgumentError in EnsembleSweepSpec constructor
        @test_throws ArgumentError EnsembleSweepSpec(
            cfg_base;
            output_dir=mktempdir(),
            sampling_method=:invalid_method,
            num_samples=2,
            parameters=Dict("grid.Nx" => [33, 65]),
            seed=42,
        )

        # 6. Integer parameter reflection helper
        @test Erebus._is_integer_parameter(cfg_base, "nonexistent", [1.0, 2.0]) === false
        @test Erebus._is_integer_parameter(cfg_base, "time.dt_initial", [1.0, 2.0]) ===
            false
        @test Erebus._is_integer_parameter(cfg_base, "time.n_steps", [1.0, 2.0]) === true
        @test Erebus._is_integer_parameter(cfg_base, "custom_key", [1, 2]) === true
    end

    @testset "Ensemble Config Relative Path Resolution" begin
        # Test loading from relative path in current directory
        spec_rel = load_ensemble_config("configs/test_ensemble_sweep.toml")
        @test spec_rel isa EnsembleSweepSpec
        @test spec_rel.base_config.time.n_steps == 2
    end

    @testset "CLI Bare Flag Parsing" begin
        parsed_bare = Erebus.parse_commandline(args=["--show_timer", "test_config.toml"])
        @test parsed_bare["show_timer"] === true
        @test parsed_bare["config_or_output"] == "test_config.toml"

        parsed_nobare = Erebus.parse_commandline(args=["test_config.toml"])
        @test parsed_nobare["show_timer"] === false
        @test parsed_nobare["config_or_output"] == "test_config.toml"
    end
end

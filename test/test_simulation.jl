@testset "Simulation" begin
    @testset "setup_dynamic_simulation_parameters(): initial state physical invariants" begin
        # Baseline default configuration
        (timestep, dt, timesum, marknum, hrsolidm, hrfluidm, YERRNOD) = Erebus.setup_dynamic_simulation_parameters()

        # Invariant Family 2: Positivity and Boundedness
        @test timestep >= 0
        @test dt > 0.0
        @test timesum >= 0.0
        @test marknum > 0

        # Radiogenic decay heat vector positivity
        @test length(hrsolidm) == 3
        @test all(hrsolidm .>= 0.0)
        @test length(hrfluidm) == 3
        @test all(hrfluidm .>= 0.0)

        # Plastic yielding error vector initialization
        @test length(YERRNOD) == nplast
        @test all(iszero, YERRNOD)

        # Explicit custom configuration propagation
        custom_time = Erebus.TimeConfig(start_step=7, dt_initial=50.0, start_time=3.5)
        cfg = Erebus.SimulationConfig(time=custom_time)
        (ts_c, dt_c, time_c, mark_c, _, _, _) = Erebus.setup_dynamic_simulation_parameters(
            cfg
        )
        @test ts_c == 7
        @test dt_c ≈ 50.0 * cfg.time.yearlength rtol=1e-12
        @test time_c ≈ 3.5 * cfg.time.yearlength rtol=1e-12
        @test mark_c == start_marknum
    end # testset "setup_dynamic_simulation_parameters()"

    @testset "s_to_Ma(): time conversion physical invariants and benchmarks" begin
        # 1. Exact astronomical benchmark: 1 Ma = 1e6 years
        @test Erebus.s_to_Ma(1e6 * yearlength) ≈ 1.0 rtol=1e-12
        @test Erebus.s_to_Ma(yearlength) ≈ 1.0e-6 rtol=1e-12

        # 2. Origin preservation
        @test iszero(Erebus.s_to_Ma(0.0))

        # 3. Linearity and homogeneity: f(a * x) = a * f(x)
        s_test = 3.15576e13 # ~1 Ma
        @test Erebus.s_to_Ma(2.5 * s_test) ≈ 2.5 * Erebus.s_to_Ma(s_test) rtol=1e-12
        @test Erebus.s_to_Ma(-s_test) ≈ -Erebus.s_to_Ma(s_test) rtol=1e-12

        # 4. Strict monotonicity
        @test Erebus.s_to_Ma(2.0e13) > Erebus.s_to_Ma(1.0e13)

        # 5. 3-class discrimination guards
        res = Erebus.s_to_Ma(1e6 * yearlength)
        @test res > 0.0 # sign
        @test abs(res - 1e6) > 10.0 # scale guard: not missing yearlength
        @test abs(res - (1e6 * yearlength)) > 10.0 # conversion guard: not identity
    end # testset "s_to_Ma()"

    @testset "Simulation Orchestration Error Contracts" begin
        # Nonexistent config file path
        @test_throws ArgumentError run_simulation("nonexistent_config_file_path.toml")

        # Config with nonexistent checkpoint restart
        bad_restart_cfg = SimulationConfig(
            output=OutputConfig(restart_from="nonexistent_checkpoint.jld2")
        )
        @test_throws ArgumentError run_simulation(bad_restart_cfg)

        # Dynamic coordinates support in setup_dynamic_simulation_parameters
        coords = Erebus.GridCoordinates(15, 15; xsize=100000.0, ysize=100000.0)
        (ts, dt, time, marknum, hrsolid, hrfluid, YERRNOD) = Erebus.setup_dynamic_simulation_parameters(;
            coords=coords
        )
        @test marknum == coords.start_marknum
        @test isapprox(dt, 1.0e11; rtol=1e-12)

        # s_to_Ma non-finite propagation
        @test isnan(Erebus.s_to_Ma(NaN))
        @test isinf(Erebus.s_to_Ma(Inf))
        @test isinf(Erebus.s_to_Ma(-Inf))
    end

    @testset "save_state() and load_state(): DT0 persistence and recovery" begin
        # 1. Direct round-trip unit test of DT0 in JLD2 checkpoint
        mktempdir() do tmpdir
            dummy_file = joinpath(tmpdir, "test_ckpt.jld2")
            test_dt0 = [i * 1.5 + j * 0.7 for i in 1:5, j in 1:5]
            JLD2.jldsave(dummy_file; DT0=test_dt0)
            loaded = load_state(dummy_file)
            @test haskey(loaded, "DT0")
            @test loaded["DT0"] == test_dt0

            # Absence test: file without DT0 does not have the key
            dummy_no_dt0 = joinpath(tmpdir, "test_no_dt0.jld2")
            JLD2.jldsave(dummy_no_dt0; timestep=1)
            loaded_no = load_state(dummy_no_dt0)
            @test !haskey(loaded_no, "DT0")
        end

        # 2. Simulation checkpoint emission and restart restoration
        mktempdir() do tmpdir
            quick_toml = joinpath(@__DIR__, "..", "configs", "test_quick.toml")
            cfg = load_config(quick_toml)
            cfg_run = SimulationConfig(
                time=TimeConfig(
                    n_steps=2,
                    dt_initial=cfg.time.dt_initial,
                    dt_longest=cfg.time.dt_longest,
                ),
                solver=cfg.solver,
                poroelasticity=cfg.poroelasticity,
                thermodynamics=cfg.thermodynamics,
                materials=cfg.materials,
                output=OutputConfig(output_dir=tmpdir, savematstep=1),
            )
            Erebus.simulation_loop(cfg_run; output_path=tmpdir)
            ckpt1_path = joinpath(tmpdir, "output_00001.jld2")
            @test isfile(ckpt1_path)
            state1 = load_state(ckpt1_path)
            @test haskey(state1, "DT0")
            @test size(state1["DT0"]) == (cfg_run.grid.Ny + 1, cfg_run.grid.Nx + 1)
            @test all(isfinite, state1["DT0"])
            @test any(!iszero, state1["DT0"])
            @test any(state1["DT0"] .> 0.0)

            # Test restart restoring DT0
            dir_res = mktempdir()
            try
                cfg_res = SimulationConfig(
                    time=TimeConfig(
                        start_step=1,
                        n_steps=2,
                        dt_initial=cfg.time.dt_initial,
                        dt_longest=cfg.time.dt_longest,
                    ),
                    solver=cfg.solver,
                    poroelasticity=cfg.poroelasticity,
                    thermodynamics=cfg.thermodynamics,
                    materials=cfg.materials,
                    output=OutputConfig(
                        output_dir=dir_res, savematstep=1, restart_from=ckpt1_path
                    ),
                )
                Erebus.simulation_loop(
                    cfg_res; output_path=dir_res, restart_from=ckpt1_path
                )
                ckpt2_res_path = joinpath(dir_res, "output_00002.jld2")
                @test isfile(ckpt2_res_path)
                state2_res = load_state(ckpt2_res_path)
                @test haskey(state2_res, "DT0")
                @test all(isfinite, state2_res["DT0"])
            finally
                rm(dir_res, recursive=true, force=true)
            end
        end
    end
end

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
end

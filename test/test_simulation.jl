using Test
using Erebus
using JLD2

@testset "Simulation" begin
    nplast = 100_000
    yearlength = Erebus.default_config().time.yearlength
    @testset "setup_dynamic_simulation_parameters(): initial state physical invariants" begin
        # Baseline default configuration
        (timestep, dt, timesum, marknum, hrsolidm, hrfluidm, YERRNOD) = Erebus.setup_dynamic_simulation_parameters()

        cfg_default = Erebus.default_config()
        @test timestep == 1
        @test dt ≈ cfg_default.time.dt_initial * cfg_default.time.yearlength rtol=1e-12
        @test timesum ≈ cfg_default.time.start_time * cfg_default.time.yearlength rtol=1e-12
        @test marknum == Erebus.start_marknum

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
        @test mark_c == Erebus.start_marknum
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

    @testset "End-to-end speciation and surface thermal coupling in simulation_loop" begin
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
                thermodynamics=ThermalConfig(surface_radiation=true),
                materials=cfg.materials,
                reaction=ReactionConfig(cfl_reaction=0.5, dphi_reaction_max=0.01),
                atmosphere=AtmosphereConfig(active=true),
                volatiles=VolatilesConfig(
                    active=true, speciation_active=true, graphite_saturation=true
                ),
                output=OutputConfig(output_dir=tmpdir, savematstep=1),
            )
            Erebus.simulation_loop(cfg_run; output_path=tmpdir)
            ckpt1_path = joinpath(tmpdir, "output_00001.jld2")
            @test isfile(ckpt1_path)
            state1 = load_state(ckpt1_path)
            @test haskey(state1, "tk1")
            @test all(isfinite, state1["tk1"])
            @test any(state1["tk1"] .> 0.0)
            @test haskey(state1, "atm_T_surf_eq")
            @test isfinite(state1["atm_T_surf_eq"])
            @test state1["atm_T_surf_eq"] > 0.0
            @test haskey(state1, "DQPF")
            @test all(isfinite, state1["DQPF"])

            ckpt2_path = joinpath(tmpdir, "output_00002.jld2")
            @test isfile(ckpt2_path)
            state2 = load_state(ckpt2_path)
            @test isfinite(state2["atm_T_surf_eq"])
            @test state2["atm_T_surf_eq"] > 0.0
        end
    end

    @testset "Simulation loop return value and state verification" begin
        mktempdir() do tmpdir
            quick_toml = joinpath(@__DIR__, "..", "configs", "test_quick.toml")
            cfg = load_config(quick_toml)
            cfg_run = SimulationConfig(
                grid=cfg.grid,
                time=cfg.time,
                poroelasticity=cfg.poroelasticity,
                output=OutputConfig(
                    output_dir=tmpdir, mode=:both, savematstep=2, telemetrystep=1
                ),
            )
            res = Erebus.simulation_loop(cfg_run; output_path=tmpdir)
            @test res isa NamedTuple
            @test haskey(res, :markers)
            @test haskey(res, :grids)
            @test haskey(res, :atm)
            @test haskey(res, :transfers)
            @test haskey(res, :timesum)
            @test haskey(res, :dt)
            @test haskey(res, :timestep)

            # Checkpoint written at step 2
            ckpt2_path = joinpath(tmpdir, "output_00002.jld2")
            @test isfile(ckpt2_path)
            state2 = load_state(ckpt2_path)
            @test res.markers.xm == state2["xm"]
            @test res.timestep == 2

            # Telemetry verification: timesum equals start_time + sum of accepted dt values
            telem_path = joinpath(tmpdir, "telemetry.csv")
            @test isfile(telem_path)
            lines = readlines(telem_path)
            @test length(lines) >= 3
            dt1_yr = parse(Float64, split(lines[2], ",")[3])
            dt2_yr = parse(Float64, split(lines[3], ",")[3])
            dt_sum_s = (dt1_yr + dt2_yr) * cfg.time.yearlength
            initial_timesum = cfg.time.start_time * cfg.time.yearlength
            @test res.timesum ≈ initial_timesum + dt_sum_s rtol=1e-12
        end
    end

    @testset "Degassing transfers and 2D water mass conservation" begin
        mktempdir() do tmpdir
            cfg_path = joinpath(
                @__DIR__, "..", "configs", "magma_ocean_cooling_turb_on_32.toml"
            )
            cfg = load_config(cfg_path)
            cfg_run = SimulationConfig(
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
                output=OutputConfig(output_dir=tmpdir, savematstep=2),
                disk=cfg.disk,
                melting=cfg.melting,
                volatiles=VolatilesConfig(initial_water_wtpct=2.0),
                magma_degassing=MagmaOceanDegassingConfig(active=true, mode=:dynamic_flux),
            )
            res = Erebus.simulation_loop(cfg_run; output_path=tmpdir)
            @test haskey(res, :transfers)
            @test res.transfers isa Vector{TransferRecord}
            degas_h = filter(
                r -> r.channel === :degassing && r.element === :H, res.transfers
            )
            @test !isempty(degas_h)
            sum_dM2_H = sum(r -> r.dM2, degas_h)
            @test !iszero(sum_dM2_H)

            coords = GridCoordinates(cfg.grid)
            Am = marker_area(coords)
            rho_s = cfg.materials.rhosolidm[1]
            H_ratio = 2.01588 / 18.01528

            init_data = load_state(joinpath(tmpdir, "output_00000.jld2"))
            init_XH2O = init_data["XH2Om"]
            init_tm = init_data["tm"]
            final_XH2O = res.markers.XH2Om

            delta_H_markers = 0.0
            for m in 1:length(init_XH2O)
                if init_tm[m] < 3
                    w_init = init_XH2O[m] * 0.01
                    w_final = final_XH2O[m] * 0.01
                    dw = w_init - w_final
                    if dw > 0.0
                        delta_H_markers += dw * (rho_s * Am) * H_ratio
                    end
                end
            end

            @test isapprox(sum_dM2_H, delta_H_markers; rtol=1e-8)
        end
    end
end

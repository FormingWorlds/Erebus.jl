using Test
using Random
using Dates
using Logging
using JLD2
using Erebus

@testset "SimulationState, Checkpoint Schema v2, and Restart Continuity" begin
    @testset "A10: override_config Unknown Field Validation" begin
        cfg = default_config()
        # Verify unknown section throws ArgumentError
        @test_throws ArgumentError Erebus.override_config(
            cfg, Dict("invalid_section.seed" => 123)
        )
        # Verify unknown field in valid section throws ArgumentError
        @test_throws ArgumentError Erebus.override_config(
            cfg, Dict("solver.does_not_exist" => 1)
        )
        @test_throws ArgumentError Erebus.override_config(
            cfg, Dict("thermodynamics.nonexistent_property" => 42.0)
        )
    end

    @testset "Global Timer Binding Elimination" begin
        src_dir = joinpath(@__DIR__, "..", "src")
        # Assert that no mutable global timer binding exists anywhere in src/
        to_matches = String[]
        for (root, _, files) in walkdir(src_dir)
            for f in files
                endswith(f, ".jl") || continue
                path = joinpath(root, f)
                for (ln, line) in enumerate(eachline(path))
                    if occursin(r"\bconst\s+to\b", line)
                        push!(to_matches, "$path:$ln: $line")
                    end
                end
            end
        end
        @test isempty(to_matches)
        @test !isdefined(Erebus, :to)
    end

    @testset "A3: CLI Config Section Preservation on Restart" begin
        cfg_base = Erebus.override_config(
            default_config(),
            Dict(
                "melting.active" => true,
                "volatiles.active" => true,
                "redox.active" => true,
                "accretion.active" => true,
            ),
        )
        @test cfg_base.melting.active == true
        @test cfg_base.volatiles.active == true
        @test cfg_base.redox.active == true
        @test cfg_base.accretion.active == true

        # Rebuild config with new restart_from path using internal helper
        dummy_restart_path = "/tmp/dummy_checkpoint.jld2"
        rebuilt = Erebus.rebuild_cli_restart_config(cfg_base, dummy_restart_path)

        # Assert every section matches cfg_base except output.restart_from
        for fn in fieldnames(SimulationConfig)
            if fn === :output
                @test getproperty(rebuilt, fn).restart_from == dummy_restart_path
                @test getproperty(rebuilt, fn).output_dir == cfg_base.output.output_dir
            else
                @test getproperty(rebuilt, fn) == getproperty(cfg_base, fn)
            end
        end
    end

    @testset "Checkpoint Schema v2 Version Validation" begin
        mktempdir() do tmpdir
            # Schema version 1 dummy checkpoint
            v1_path = joinpath(tmpdir, "v1_checkpoint.jld2")
            jldsave(v1_path; schema_version=1, timestep=1, dt=1.0, timesum=1.0)
            @test_throws Erebus.CheckpointError Erebus.load_simulation_state(v1_path)

            # Unversioned dummy checkpoint
            unversioned_path = joinpath(tmpdir, "unversioned_checkpoint.jld2")
            jldsave(unversioned_path; timestep=1, dt=1.0, timesum=1.0)
            @test_throws Erebus.CheckpointError Erebus.load_simulation_state(
                unversioned_path
            )
        end
    end

    @testset "Section 3 Restart Configuration Comparison" begin
        cfg_saved = default_config()

        # 1. Output differences and time progression differences are allowed
        cfg_allowed = Erebus.override_config(
            cfg_saved,
            Dict(
                "output.output_dir" => "/tmp/new_dir",
                "time.n_steps" => 10,
                "time.endtime" => 5.0e7,
                "time.start_step" => 3,
            ),
        )
        diffs_allowed = Erebus.compare_restart_configs(cfg_saved, cfg_allowed)
        @test isempty(diffs_allowed)

        # 2. Physics section differences are rejected unless force_restart_config is set
        cfg_physics_changed = Erebus.override_config(
            cfg_saved,
            Dict(
                "thermodynamics.hr_al" => !cfg_saved.thermodynamics.hr_al,
                "poroelasticity.betasolid" => 5.0e-11,
            ),
        )
        diffs_physics = Erebus.compare_restart_configs(cfg_saved, cfg_physics_changed)
        @test "thermodynamics.hr_al" in diffs_physics
        @test "poroelasticity.betasolid" in diffs_physics
        @test length(diffs_physics) == 2

        # 3. Solver seed difference is allowed on restart
        cfg_seed = Erebus.override_config(cfg_saved, Dict("solver.seed" => 99999))
        @test isempty(Erebus.compare_restart_configs(cfg_saved, cfg_seed))
    end

    @testset "Missing Checkpoint File Throws CheckpointError" begin
        @test_throws Erebus.CheckpointError Erebus.load_simulation_state(
            "/nonexistent/checkpoint_path_12345.jld2"
        )
    end

    @testset "SimulationState Property Delegation and propertynames" begin
        mktempdir() do tmpdir
            cfg = Erebus.override_config(
                load_config(joinpath(@__DIR__, "..", "configs", "test_quick.toml")),
                Dict(
                    "grid.Nx" => 17,
                    "grid.Ny" => 17,
                    "time.n_steps" => 1,
                    "output.output_dir" => tmpdir,
                ),
            )
            s = simulation_loop(cfg)
            @test hasproperty(s, :ETA)
            @test hasproperty(s, :xcenter)
            @test hasproperty(s, :grids)
            @test hasproperty(s, :markers)
            @test s.ETA === s.grids.ETA
            @test s.xcenter === s.accumulators.xcenter
            pnames = propertynames(s)
            @test :ETA in pnames
            @test :xcenter in pnames
            @test :grids in pnames
        end
    end

    @testset "Bitwise Restart Continuity (Straight vs Split Run)" begin
        mktempdir() do tmpdir
            straight_dir = joinpath(tmpdir, "straight")
            split_dir = joinpath(tmpdir, "split")

            cfg_base = Erebus.override_config(
                load_config(joinpath(@__DIR__, "..", "configs", "test_quick.toml")),
                Dict(
                    "grid.Nx" => 17,
                    "grid.Ny" => 17,
                    "time.n_steps" => 4,
                    "output.output_dir" => straight_dir,
                    "output.savematstep" => 1,
                    "output.mode" => :snapshots,
                ),
            )

            s_straight = simulation_loop(cfg_base)

            ckpt_step2 = joinpath(straight_dir, "output_00002.jld2")
            cfg_resume = Erebus.override_config(
                cfg_base,
                Dict(
                    "time.n_steps" => 4,
                    "output.output_dir" => split_dir,
                    "output.restart_from" => ckpt_step2,
                ),
            )
            s_resume = simulation_loop(cfg_resume)

            # Assert bitwise reproducibility across all grids
            for fn in fieldnames(GridArrays)
                v1 = getfield(s_straight.grids, fn)
                v2 = getfield(s_resume.grids, fn)
                @test isequal(v1, v2)
            end

            # Assert bitwise reproducibility across all marker arrays
            for fn in fieldnames(typeof(s_straight.markers.core))
                v1 = getfield(s_straight.markers.core, fn)
                v2 = getfield(s_resume.markers.core, fn)
                @test isequal(v1, v2)
            end

            # Assert bitwise reproducibility of RNG, accumulators, transfers, atm, timing
            @test s_straight.rng == s_resume.rng
            for fn in fieldnames(SimulationAccumulators)
                v1 = getfield(s_straight.accumulators, fn)
                v2 = getfield(s_resume.accumulators, fn)
                @test isequal(v1, v2)
            end
            @test s_straight.transfers == s_resume.transfers
            @test s_straight.atm == s_resume.atm
            @test s_straight.timestep == s_resume.timestep
            @test s_straight.timesum == s_resume.timesum
            @test s_straight.dt == s_resume.dt
        end
    end

    @testset "Timer Hierarchy Cleanliness on Restart" begin
        mktempdir() do tmpdir
            cfg_base = Erebus.override_config(
                load_config(joinpath(@__DIR__, "..", "configs", "test_quick.toml")),
                Dict(
                    "grid.Nx" => 17,
                    "grid.Ny" => 17,
                    "time.n_steps" => 2,
                    "output.output_dir" => tmpdir,
                    "output.savematstep" => 1,
                    "output.mode" => :snapshots,
                ),
            )
            simulation_loop(cfg_base)
            ckpt_path = joinpath(tmpdir, "output_00002.jld2")
            (state, _, _) = Erebus.load_simulation_state(ckpt_path)
            # timer_stack must be empty to avoid nesting resumed steps
            @test isempty(state.timer.timer_stack)
            @test state.timestep == 2
        end
    end

    @testset "Force Restart Rejects Grid Geometry Override" begin
        mktempdir() do tmpdir
            cfg_base = Erebus.override_config(
                load_config(joinpath(@__DIR__, "..", "configs", "test_quick.toml")),
                Dict(
                    "grid.Nx" => 17,
                    "grid.Ny" => 17,
                    "time.n_steps" => 2,
                    "output.output_dir" => tmpdir,
                    "output.savematstep" => 1,
                    "output.mode" => :snapshots,
                ),
            )
            simulation_loop(cfg_base)
            ckpt_path = joinpath(tmpdir, "output_00002.jld2")

            cfg_mismatch = Erebus.override_config(
                cfg_base, Dict("grid.Nx" => 21, "grid.Ny" => 21)
            )
            @test_throws Erebus.CheckpointError Erebus.load_simulation_state(
                ckpt_path; current_cfg=cfg_mismatch, force_restart_config=true
            )
        end
    end

    @testset "Bitwise Restart Continuity (Multi-Group, Off-Center, Atmosphere)" begin
        mktempdir() do tmpdir
            straight_dir = joinpath(tmpdir, "straight_multi")
            split_dir = joinpath(tmpdir, "split_multi")

            cfg_base = Erebus.override_config(
                load_config(joinpath(@__DIR__, "..", "configs", "test_quick.toml")),
                Dict(
                    "grid.Nx" => 17,
                    "grid.Ny" => 17,
                    "grid.xsize" => 70_000.0,
                    "grid.ysize" => 70_000.0,
                    "geometry.xcenter" => 30_000.0,
                    "geometry.ycenter" => 30_000.0,
                    "geometry.rplanet" => 20_000.0,
                    "geometry.rcrust" => 20_000.0,
                    "volatiles.active" => true,
                    "atmosphere.active" => true,
                    "metal_partition.active" => true,
                    "time.n_steps" => 4,
                    "output.output_dir" => straight_dir,
                    "output.savematstep" => 1,
                    "output.mode" => :snapshots,
                ),
            )

            s_straight = simulation_loop(cfg_base)

            ckpt_step2 = joinpath(straight_dir, "output_00002.jld2")
            cfg_resume = Erebus.override_config(
                cfg_base,
                Dict(
                    "time.n_steps" => 4,
                    "output.output_dir" => split_dir,
                    "output.restart_from" => ckpt_step2,
                ),
            )
            s_resume = simulation_loop(cfg_resume)

            # Assert bitwise reproducibility across all grids
            for fn in fieldnames(GridArrays)
                v1 = getfield(s_straight.grids, fn)
                v2 = getfield(s_resume.grids, fn)
                @test isequal(v1, v2)
            end

            # Assert bitwise reproducibility across core markers
            for fn in fieldnames(typeof(s_straight.markers.core))
                v1 = getfield(s_straight.markers.core, fn)
                v2 = getfield(s_resume.markers.core, fn)
                @test isequal(v1, v2)
            end

            # Assert bitwise reproducibility across all optional marker groups
            @test keys(s_straight.markers.groups) == keys(s_resume.markers.groups)
            for grp_name in keys(s_straight.markers.groups)
                g1 = s_straight.markers.groups[grp_name]
                g2 = s_resume.markers.groups[grp_name]
                for fn in fieldnames(typeof(g1))
                    @test isequal(getfield(g1, fn), getfield(g2, fn))
                end
            end

            # Assert bitwise reproducibility of RNG, accumulators, transfers, atm, timing
            @test s_straight.rng == s_resume.rng
            for fn in fieldnames(SimulationAccumulators)
                v1 = getfield(s_straight.accumulators, fn)
                v2 = getfield(s_resume.accumulators, fn)
                @test isequal(v1, v2)
            end
            @test s_straight.accumulators.xcenter ≈ 30_000.0
            @test s_straight.accumulators.ycenter ≈ 30_000.0
            @test s_resume.accumulators.xcenter ≈ 30_000.0
            @test s_resume.accumulators.ycenter ≈ 30_000.0

            @test s_straight.transfers == s_resume.transfers
            @test s_straight.atm == s_resume.atm
            @test s_straight.timestep == s_resume.timestep
            @test isequal(s_straight.timesum, s_resume.timesum)
            @test isequal(s_straight.dt, s_resume.dt)
        end
    end

    @testset "SimulationState and GridArrays Complete Dictionary and Indexing Interface" begin
        # 1. CheckpointError showerror
        io = IOBuffer()
        showerror(io, CheckpointError("test message"))
        err_str = String(take!(io))
        @test occursin("CheckpointError:", err_str)
        @test occursin("test message", err_str)

        mktempdir() do tmpdir
            cfg = Erebus.override_config(
                load_config(joinpath(@__DIR__, "..", "configs", "test_quick.toml")),
                Dict(
                    "grid.Nx" => 17,
                    "grid.Ny" => 17,
                    "time.n_steps" => 1,
                    "output.output_dir" => tmpdir,
                ),
            )
            s = simulation_loop(cfg)
            g = s.grids

            # 2. GridArrays keys, pairs, haskey, getindex, propertynames, copy
            @test :ETA in keys(g)
            @test haskey(g, :ETA)
            @test haskey(g, "ETA")
            @test !haskey(g, :nonexistent_field)
            @test !haskey(g, "nonexistent_field")
            @test g[:ETA] === g.ETA
            @test g["ETA"] === g.ETA
            @test :ETA in propertynames(g)

            g_copy = copy(g)
            @test g_copy.ETA == g.ETA
            @test g_copy.ETA !== g.ETA

            pairs_list = collect(pairs(g))
            @test length(pairs_list) == length(fieldnames(GridArrays))
            @test any(p -> p.first === :ETA && p.second === g.ETA, pairs_list)

            # 3. GridArrays copy with Q_metric present
            g_metric = GridArrays(
                (
                    fn === :Q_metric ? zeros(Float64, 18, 18) : getfield(g, fn) for
                    fn in fieldnames(GridArrays)
                )...,
            )
            g_metric_copy = copy(g_metric)
            @test g_metric_copy.Q_metric !== g_metric.Q_metric
            @test size(g_metric_copy.Q_metric) == (18, 18)

            # 4. SimulationAccumulators copy
            acc = s.accumulators
            acc_copy = copy(acc)
            @test acc_copy !== acc
            @test acc_copy.xcenter ≈ acc.xcenter
            @test acc_copy.ycenter ≈ acc.ycenter
            @test acc_copy.max_v_seg_prev ≈ acc.max_v_seg_prev

            # 5. SimulationState getproperty and error handling
            @test s.timestep == 1
            @test s.xcenter ≈ acc.xcenter
            @test s.ETA === g.ETA
            @test_throws ErrorException s.nonexistent_property

            # 6. SimulationState keys, haskey, pairs
            all_keys = keys(s)
            @test :grids in all_keys
            @test :markers in all_keys
            @test :ETA in all_keys
            @test :transfer_log in all_keys
            @test :S_vent in all_keys
            @test haskey(s, :grids)
            @test haskey(s, "grids")
            @test haskey(s, :ETA)
            @test haskey(s, "ETA")
            @test haskey(s, :transfer_log)
            @test haskey(s, "transfer_log")
            @test haskey(s, :S_vent)
            @test haskey(s, "S_vent")
            @test !haskey(s, :nonexistent_symbol)
            @test !haskey(s, "nonexistent_symbol")

            spairs = collect(pairs(s))
            @test any(p -> p.first === :transfer_log, spairs)
            @test any(p -> p.first === :S_vent, spairs)

            # 7. SimulationState getindex and KeyError
            @test s[:transfer_log] === s.transfers
            @test s["transfer_log"] === s.transfers
            @test s[:S_vent] === s.grids.S_vent_grid
            @test s["S_vent"] === s.grids.S_vent_grid
            @test s[:timestep] == s.timestep
            @test s["timestep"] == s.timestep
            @test s[:xcenter] ≈ s.accumulators.xcenter
            @test s["xcenter"] ≈ s.accumulators.xcenter
            @test s[:ETA] === s.grids.ETA
            @test s["ETA"] === s.grids.ETA
            @test_throws KeyError s[:nonexistent_key]
            @test_throws KeyError s["nonexistent_key"]

            # 8. SimulationState copy and NamedTuple
            nt = NamedTuple(s)
            @test nt.timestep == s.timestep
            @test nt.grids === s.grids
            @test nt.markers === s.markers

            scopy = copy(s)
            @test scopy.timestep == s.timestep
            @test scopy.dt ≈ s.dt
            @test scopy.timesum ≈ s.timesum
            @test scopy.grids.ETA == s.grids.ETA
            @test scopy.grids.ETA !== s.grids.ETA
        end
    end

    @testset "Checkpoint Schema v2 Validation and Error Branches" begin
        mktempdir() do tmpdir
            # 1. Corrupted checkpoint archive (JLD2 cannot load)
            corrupt_path = joinpath(tmpdir, "corrupt.jld2")
            write(corrupt_path, "not a valid jld2 file content")
            @test_throws Erebus.CheckpointError Erebus.load_simulation_state(corrupt_path)
            nonexistent_path = joinpath(tmpdir, "does_not_exist.jld2")
            @test_throws Erebus.CheckpointError Erebus.load_simulation_state(
                nonexistent_path
            )

            # 2. Checkpoint missing required keys ("cfg", "timestep", "dt", "timesum", "rng", "transfers")
            missing_req_path = joinpath(tmpdir, "missing_req.jld2")
            jldsave(missing_req_path; schema_version=2, timestep=1)
            @test_throws Erebus.CheckpointError Erebus.load_simulation_state(
                missing_req_path
            )

            # 3. Invalid 'cfg' payload
            invalid_cfg_path = joinpath(tmpdir, "invalid_cfg.jld2")
            jldsave(
                invalid_cfg_path;
                schema_version=2,
                cfg="not_a_SimulationConfig",
                timestep=1,
                dt=1.0,
                timesum=1.0,
                rng=MersenneTwister(42),
                transfers=TransferRecord[],
            )
            @test_throws Erebus.CheckpointError Erebus.load_simulation_state(
                invalid_cfg_path
            )

            # 4. Save state with explicit .jld2 path and verify readback
            cfg = Erebus.override_config(
                load_config(joinpath(@__DIR__, "..", "configs", "test_quick.toml")),
                Dict(
                    "grid.Nx" => 17,
                    "grid.Ny" => 17,
                    "time.n_steps" => 1,
                    "output.output_dir" => tmpdir,
                    "output.savematstep" => 1,
                ),
            )
            s = simulation_loop(cfg)
            ckpt_path = joinpath(tmpdir, "output_00001.jld2")
            explicit_path = joinpath(tmpdir, "explicit_save.jld2")
            saved_ret = Erebus.save_state(
                explicit_path, s, Erebus.GridCoordinates(cfg.grid), cfg
            )
            @test saved_ret == explicit_path
            (s_read, _, _) = Erebus.load_simulation_state(explicit_path)
            @test s_read.timestep == s.timestep
            @test isequal(s_read.grids.ETA, s.grids.ETA)

            # 5. Missing required grid array
            missing_eta_path = joinpath(tmpdir, "missing_eta.jld2")
            data = JLD2.load(ckpt_path)
            delete!(data, "ETA")
            JLD2.jldsave(missing_eta_path; (Symbol(k) => v for (k, v) in data)...)
            @test_throws Erebus.CheckpointError Erebus.load_simulation_state(
                missing_eta_path
            )

            # 6. Missing required core array xm
            missing_xm_path = joinpath(tmpdir, "missing_xm.jld2")
            data_core = JLD2.load(ckpt_path)
            delete!(data_core, "xm")
            JLD2.jldsave(missing_xm_path; (Symbol(k) => v for (k, v) in data_core)...)
            @test_throws Erebus.CheckpointError Erebus.load_simulation_state(
                missing_xm_path
            )

            # 7. Fallback loading without coords, accumulators, timer, atm
            data_fallback = JLD2.load(ckpt_path)
            delete!(data_fallback, "coords")
            delete!(data_fallback, "accumulators")
            delete!(data_fallback, "timer")
            delete!(data_fallback, "atm_state")
            fallback_path = joinpath(tmpdir, "fallback.jld2")
            JLD2.jldsave(fallback_path; (Symbol(k) => v for (k, v) in data_fallback)...)
            (s_fb, coords_fb, cfg_fb) = Erebus.load_simulation_state(fallback_path)
            @test s_fb.timestep == 1
            @test s_fb.accumulators.rplanet ≈ 50000.0

            # 8. Missing Q_metric in checkpoint loads as nothing
            missing_q_path = joinpath(tmpdir, "missing_q.jld2")
            data_q = JLD2.load(ckpt_path)
            delete!(data_q, "Q_metric")
            JLD2.jldsave(missing_q_path; (Symbol(k) => v for (k, v) in data_q)...)
            (s_no_q, _, _) = Erebus.load_simulation_state(missing_q_path)
            @test s_no_q.grids.Q_metric === nothing
            @test s_no_q.timestep == 1

            # 9. Config mismatch without force_restart_config throws CheckpointError
            cfg_override = Erebus.override_config(
                cfg, Dict("thermodynamics.hr_al" => !cfg.thermodynamics.hr_al)
            )
            @test_throws Erebus.CheckpointError Erebus.load_simulation_state(
                ckpt_path; current_cfg=cfg_override, force_restart_config=false
            )

            # 10. Force restart with permitted override triggering warning
            s_warn = @test_logs (:warn, r"thermodynamics\.hr_al") match_mode=:any begin
                (sw, _, _) = Erebus.load_simulation_state(
                    ckpt_path; current_cfg=cfg_override, force_restart_config=true
                )
                sw
            end
            @test s_warn.timestep == 1
            @test isequal(s_warn.dt, s.dt)
        end
    end

    @testset "Checkpoint Optional Marker Groups Loading Coverage" begin
        mktempdir() do tmpdir
            cfg_base = Erebus.override_config(
                load_config(joinpath(@__DIR__, "..", "configs", "test_quick.toml")),
                Dict(
                    "grid.Nx" => 17,
                    "grid.Ny" => 17,
                    "time.n_steps" => 1,
                    "output.output_dir" => tmpdir,
                    "output.savematstep" => 1,
                ),
            )
            simulation_loop(cfg_base)
            ckpt_path = joinpath(tmpdir, "output_00001.jld2")

            cfg_groups = Erebus.override_config(
                cfg_base,
                Dict(
                    "volatiles.active" => true,
                    "metal_partition.active" => true,
                    "phase_tracking.active" => true,
                    "redox.active" => true,
                    "volatile_mixture.active" => true,
                    "accretion.active" => true,
                ),
            )
            (s_loaded, _, _) = Erebus.load_simulation_state(
                ckpt_path; current_cfg=cfg_groups, force_restart_config=true
            )
            @test haskey(s_loaded.markers.groups, :redox)
            @test haskey(s_loaded.markers.groups, :hcnspo)
            @test haskey(s_loaded.markers.groups, :phase)
            @test haskey(s_loaded.markers.groups, :accretion)
        end
    end

    @testset "CLI and run_simulation Coverage" begin
        mktempdir() do tmpdir
            toml_path = joinpath(tmpdir, "test_cli.toml")
            cfg = Erebus.override_config(
                load_config(joinpath(@__DIR__, "..", "configs", "test_quick.toml")),
                Dict(
                    "grid.Nx" => 17,
                    "grid.Ny" => 17,
                    "time.n_steps" => 1,
                    "output.output_dir" => joinpath(tmpdir, "cli_out"),
                    "output.savematstep" => 1,
                ),
            )
            Erebus.save_config(toml_path, cfg)

            orig_logger = global_logger()
            orig_args = copy(ARGS)
            try
                # 1. run_simulation with toml file path
                s1 = run_simulation(toml_path)
                @test s1 isa SimulationState
                @test s1.timestep == 1

                # 2. run_simulation with directory path override
                out_sub = joinpath(tmpdir, "sub_out")
                s2 = run_simulation(cfg; output_path=out_sub)
                @test s2 isa SimulationState
                @test isdir(out_sub)

                # 3. run_simulation with preloaded config and restart_from
                ckpt = joinpath(tmpdir, "cli_out", "output_00001.jld2")
                cfg_step2 = Erebus.override_config(cfg, Dict("time.n_steps" => 2))
                s3 = run_simulation(cfg_step2; restart_from=ckpt, force_restart_config=true)
                @test s3 isa SimulationState
                @test s3.timestep == 2

                # 4. rebuild_cli_restart_config
                rebuilt = Erebus.rebuild_cli_restart_config(cfg, ckpt)
                @test rebuilt.output.restart_from == ckpt
                @test rebuilt.grid.Nx == cfg.grid.Nx

                # 5. run_simulation with CLI ARGS dispatch
                toml_step2 = joinpath(tmpdir, "test_cli_step2.toml")
                Erebus.save_config(toml_step2, cfg_step2)
                empty!(ARGS)
                push!(ARGS, toml_step2)
                push!(ARGS, "--restart", ckpt)
                push!(ARGS, "--force-restart-config")
                push!(ARGS, "--show_timer", "true")
                s_cli = run_simulation("")
                @test s_cli isa SimulationState
                @test s_cli.timestep == 2
            finally
                empty!(ARGS)
                append!(ARGS, orig_args)
                global_logger(orig_logger)
            end
        end
    end
end

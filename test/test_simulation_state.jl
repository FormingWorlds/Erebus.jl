using Test
using Random
using Dates
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
            @test s_straight.timesum == s_resume.timesum
            @test s_straight.dt == s_resume.dt
        end
    end
end

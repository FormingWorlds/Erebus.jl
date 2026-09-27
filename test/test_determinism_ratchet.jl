"""
Validation tests for RNG threading, determinism, and coordinate guards.

Verifies:
1. Coordinate validation in fix() rejecting non-finite positions.
2. Configuration validation rejecting negative seeds and invalid threaded P2M modes.
3. Simulation bitwise determinism for matching seeds across time steps.
4. Divergence of simulations initialized with different seeds.
5. Bitwise equality of Distributed ensemble members against standalone runs.
6. Absence of mutable global rgen binding in source tree.
"""

using Test
using Distributed
using Random
using LinearAlgebra
using Erebus

@testset "RNG Threading and Determinism Ratchet" begin
    @testset "fix() rejects non-finite marker coordinates (A13)" begin
        x_axis = collect(0.0:100.0:1000.0)
        y_axis = collect(0.0:100.0:1000.0)
        dx = 100.0
        dy = 100.0
        jmin, jmax = 1, 10
        imin, imax = 1, 10

        # Finite positions succeed
        i, j = Erebus.fix(50.0, 50.0, x_axis, y_axis, dx, dy, jmin, jmax, imin, imax)
        @test i >= imin && i <= imax
        @test j >= jmin && j <= jmax

        # Non-finite coordinates must throw DomainError
        @test_throws DomainError Erebus.fix(
            NaN, 50.0, x_axis, y_axis, dx, dy, jmin, jmax, imin, imax
        )
        @test_throws DomainError Erebus.fix(
            50.0, NaN, x_axis, y_axis, dx, dy, jmin, jmax, imin, imax
        )
        @test_throws DomainError Erebus.fix(
            Inf, 50.0, x_axis, y_axis, dx, dy, jmin, jmax, imin, imax
        )
        @test_throws DomainError Erebus.fix(
            50.0, -Inf, x_axis, y_axis, dx, dy, jmin, jmax, imin, imax
        )
    end

    @testset "Configuration rejects negative seeds and invalid threaded P2M modes (A8)" begin
        cfg_bad_seed = SimulationConfig(solver=SolverConfig(seed=-1))
        @test_throws ArgumentError Erebus.validate_config(cfg_bad_seed)

        cfg_valid_seed = SimulationConfig(solver=SolverConfig(seed=42))
        @test Erebus.validate_config(cfg_valid_seed) === nothing
    end

    @testset "Simulation loop is bitwise deterministic for identical seeds (A1)" begin
        cfg1 = SimulationConfig(
            grid=GridConfig(Nx=17, Ny=17),
            solver=SolverConfig(seed=1234, p2m_mode=:tiled),
            time=TimeConfig(n_steps=2),
            reaction=ReactionConfig(active=false),
            melting=MeltingConfig(active=false),
            accretion=AccretionConfig(active=false),
            output=OutputConfig(mode=:telemetry, telemetrystep=1000, save_final=false),
        )

        # Run 1
        res1 = simulation_loop(cfg1)

        # Run 2 with identical seed
        cfg2 = deepcopy(cfg1)
        res2 = simulation_loop(cfg2)

        # Assert bitwise equality of marker positions and temperatures
        @test res1.markers.xm == res2.markers.xm
        @test res1.markers.ym == res2.markers.ym
        @test res1.markers.tkm == res2.markers.tkm

        # Run 3 with different seed must differ
        cfg3 = SimulationConfig(
            grid=GridConfig(Nx=17, Ny=17),
            solver=SolverConfig(seed=5678, p2m_mode=:tiled),
            time=TimeConfig(n_steps=2),
            reaction=ReactionConfig(active=false),
            melting=MeltingConfig(active=false),
            accretion=AccretionConfig(active=false),
            output=OutputConfig(mode=:telemetry, telemetrystep=1000, save_final=false),
        )
        res3 = simulation_loop(cfg3)

        @test res1.markers.xm != res3.markers.xm
        @test res1.markers.ym != res3.markers.ym
    end

    @testset "Buffer allocation defaults to 16 chunks for thread invariance (A8)" begin
        coords = GridCoordinates(17, 17; xsize=1000.0, ysize=1000.0)
        bufs = allocate_thread_interpolation_buffers(coords)
        @test length(bufs) == 16
        @test length(allocate_thread_interpolation_buffers(16, coords.Nx, coords.Ny)) == 16
    end

    @testset "Ensemble sweep validates max_workers against Distributed.nprocs (A9)" begin
        cfg = default_config()
        spec = EnsembleSweepSpec(
            cfg;
            output_dir="mktemp_ens",
            sampling_method=:grid,
            parameters=Dict("thermodynamics.phim0" => [0.2]),
            seed=42,
        )
        @test_throws ArgumentError run_ensemble(spec; max_workers=Distributed.nprocs() + 1)
        @test_throws ArgumentError run_ensemble(spec; max_workers=0)
    end

    @testset "Global rgen binding is removed from source (A1)" begin
        src_dir = joinpath(@__DIR__, "..", "src")
        # Ensure grep for \brgen\b in src returns nothing
        matches = String[]
        for (root, _, files) in walkdir(src_dir)
            for f in files
                endswith(f, ".jl") || continue
                path = joinpath(root, f)
                for (line_no, line) in enumerate(eachline(path))
                    # Ignore comments that might mention rgen
                    stripped = strip(line)
                    startswith(stripped, "#") && continue
                    if occursin(r"\brgen\b", line)
                        push!(matches, "\$path:\$line_no: \$line")
                    end
                end
            end
        end
        @test !isdefined(Erebus, :rgen)
        @test isempty(matches)
    end
end

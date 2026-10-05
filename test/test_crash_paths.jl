using Test
using Erebus
using StaticArrays
using Random
using ArgParse

@testset "Crash Path Robustness" begin
    @testset "Direct solver configuration validation" begin
        cfg_pardiso = SimulationConfig(solver=SolverConfig(use_pardiso=true))
        @test_throws ArgumentError validate_config(cfg_pardiso)
    end

    @testset "Iterative and matrix-free experimental gate validation" begin
        # :iterative without experimental = true should throw ArgumentError
        cfg_iter = SimulationConfig(
            solver=SolverConfig(hydromech_solver=:iterative, experimental=false)
        )
        @test_throws ArgumentError validate_config(cfg_iter)

        # :matrix_free without experimental = true should throw ArgumentError
        cfg_mf = SimulationConfig(
            solver=SolverConfig(
                hydromech_solver=:matrix_free, darcy_elimination=true, experimental=false
            ),
        )
        @test_throws ArgumentError validate_config(cfg_mf)

        # With experimental = true, both should pass validation
        cfg_iter_exp = SimulationConfig(
            solver=SolverConfig(hydromech_solver=:iterative, experimental=true)
        )
        @test validate_config(cfg_iter_exp) === nothing

        cfg_mf_exp = SimulationConfig(
            solver=SolverConfig(
                hydromech_solver=:matrix_free, darcy_elimination=true, experimental=true
            ),
        )
        @test validate_config(cfg_mf_exp) === nothing
    end

    @testset "Command-line interface argument parsing" begin
        parsed_timer = Erebus.parse_commandline(args=["--show_timer", "custom_config.toml"])
        @test parsed_timer["show_timer"] === true
        @test parsed_timer["config_or_output"] == "custom_config.toml"

        parsed_default = Erebus.parse_commandline(args=["custom_config.toml"])
        @test parsed_default["show_timer"] === false
    end

    @testset "Pyrolysis density and mass balance" begin
        cfg_refr = RefractoryConfig(active=true, kinetics_active=true, T_pyro_min=200.0)
        marknum = 10
        tkm = fill(500.0, marknum)
        dt = 1000.0
        phim = fill(0.1, marknum)
        X_C = fill(0.02, marknum)
        X_N = fill(0.001, marknum)
        X_H = fill(0.005, marknum)
        tm = [1, 1, 2, 2, 3, 3, 1, 2, 1, 2]
        rhosolidm = SVector{3,Float64}([3000.0, 3200.0, 0.0])

        # Test that calling update_marker_pyrolysis! with rhosolidm does not throw TypeError
        res = Erebus.update_marker_pyrolysis!(
            tkm, dt, phim, X_C, X_N, X_H, cfg_refr; tm=tm, rhosolidm=rhosolidm
        )
        @test isapprox(res.total_dC_gas, 1.02238e-5, rtol=1e-3)
        @test isapprox(res.total_dN_gas, 1.27798e-6, rtol=1e-3)
        @test isapprox(res.total_dH_gas, 7.78865e-4, rtol=1e-3)
        @test haskey(res, :total_dCO_gas) || haskey(res, :total_dC_gas)

        # Edge case: zero carbon with hydrogen present under redox speciation
        redox_props = (
            nFe0_m=zeros(Float64, marknum),
            nFe2_m=fill(0.1, marknum),
            nFe3_m=zeros(Float64, marknum),
            deltaIW_m=fill(-2.0, marknum),
            nC_graphite_m=zeros(Float64, marknum),
            nCO_m=zeros(Float64, marknum),
            nCO2_m=zeros(Float64, marknum),
            nCH4_m=zeros(Float64, marknum),
        )
        X_C_zero = zeros(Float64, marknum)
        X_H_zero = fill(0.005, marknum)
        res_zero_c = Erebus.update_marker_pyrolysis!(
            tkm,
            dt,
            phim,
            X_C_zero,
            X_N,
            X_H_zero,
            cfg_refr;
            tm=tm,
            rhosolidm=rhosolidm,
            redox_props=redox_props,
            redox_cfg=RedoxConfig(active=true),
        )
        @test iszero(res_zero_c.total_dC_gas)
        @test isapprox(res_zero_c.total_dH_gas, 7.78865e-4, rtol=1e-3)
    end

    @testset "Marker array telescoping invariant guard and RedoxGroup extension" begin
        cfg = SimulationConfig(
            redox=RedoxConfig(active=true), volatiles=VolatilesConfig(active=true)
        )
        old_coords = GridCoordinates(cfg.grid)
        new_coords = Erebus.compute_telescoped_coordinates(old_coords)
        marknum = old_coords.Nxm * old_coords.Nym
        markers = init_marker_arrays(marknum, cfg, old_coords)

        # Invariant assertion should pass for initial markers
        Erebus.assert_marker_arrays_invariants(markers, marknum)

        # Telescope marker arrays
        new_marknum = Erebus.telescope_marker_arrays!(
            markers; old_coords=old_coords, new_coords=new_coords, cfg=cfg
        )
        @test new_marknum > marknum
        @test length(markers.core.xm) == new_marknum
        @test length(markers.core.w3d_m) == new_marknum
        @test length(markers.X_graphite_m) == new_marknum
        @test length(markers.nFe0_m) == new_marknum
        @test length(markers.deltaIW_m) == new_marknum
        @test all(m -> -6.0 <= markers.deltaIW_m[m] <= 6.0, 1:new_marknum)

        # Check invariant holds for new marknum
        Erebus.assert_marker_arrays_invariants(markers, new_marknum)

        # Check deliberate desynchronization throws DimensionMismatch
        push!(markers.core.xm, 0.0)
        @test_throws DimensionMismatch Erebus.assert_marker_arrays_invariants(
            markers, new_marknum
        )
        pop!(markers.core.xm)
    end

    @testset "Atmospheric oxygen excess handling" begin
        # Hydrogen depleted, excess oxygen beyond stoichiometric capacity
        # nO_max = 0.5 * nH + 2.0 * nC + 2.0 * nS
        # Create an inventory where elem.O > mO_max
        nH = 1.0e3
        nC = 1.0e4
        nS = 1.0e3
        mH = nH * 0.001008
        mC = nC * 0.012011
        mS = nS * 0.03206
        nO_max = 0.5 * nH + 2.0 * nC + 2.0 * nS
        mO_max = nO_max * 0.015999
        mO_excess = mO_max * 1.5 # 50% surplus oxygen

        elem = ElementInventory(mH, mC, 1.0e2, mS, mO_excess)
        T_surf = 300.0
        P_surf = 1.0e5

        # speciate_closed_system should return maximum oxidation state and excess_O without throwing ConvergenceError
        spec_res = Erebus.speciate_closed_system(elem, T_surf, P_surf; handle_excess_O=true)
        @test spec_res.log10_fO2 ≈ 20.0
        @test isapprox(spec_res.excess_O, mO_excess - mO_max, rtol=1e-3)

        # equilibrate_atmospheric_speciation! should transfer surplus oxygen to dO_buffer with positive credit
        atm = AtmosphereState(; elem=elem)
        Erebus.equilibrate_atmospheric_speciation!(atm, T_surf, P_surf)
        surplus_expected = mO_excess - mO_max
        @test isapprox(atm.dO_buffer, surplus_expected; rtol=1e-3)
        @test isapprox(atm.dO_buffer + atm.elem.O, mO_excess; rtol=1e-5)
    end

    @testset "Checkpoint restart preserves MetalGroup and Xfe_bulk" begin
        cfg_save = SimulationConfig(
            coreformation=CoreFormationConfig(percolation_active=false, Xfe_bulk=0.325),
            metal_partition=MetalPartitionConfig(active=false),
        )
        coords = GridCoordinates(cfg_save.grid)
        marknum = coords.Nxm * coords.Nym
        markers = init_marker_arrays(marknum, cfg_save, coords)
        Erebus.define_markers!(
            markers;
            coords=coords,
            rplanet_val=cfg_save.geometry.rplanet,
            Xfe_bulk_val=cfg_save.coreformation.Xfe_bulk,
        )

        function _make_test_grids(c)
            Nx, Ny, Nx1, Ny1 = c.Nx, c.Ny, c.Nx1, c.Ny1
            basic = (
                :ETA,
                :ETA0,
                :GGG,
                :EXY,
                :SXY,
                :SXY0,
                :wyx,
                :COH,
                :TEN,
                :FRI,
                :YNY,
                :ETA5,
                :ETA00,
                :YNY5,
                :YNY00,
                :YNY_inv_ETA,
                :DSXY,
                :DSY,
            )
            g_args = Any[]
            for fn in fieldnames(GridArrays)
                if fn === :Q_metric
                    push!(g_args, nothing)
                elseif fn in basic
                    push!(
                        g_args,
                        if fn in (:YNY, :YNY5, :YNY00)
                            zeros(Bool, Ny, Nx)
                        else
                            zeros(Float64, Ny, Nx)
                        end,
                    )
                else
                    push!(g_args, zeros(Float64, Ny1, Nx1))
                end
            end
            return GridArrays(g_args...)
        end

        grids = _make_test_grids(coords)
        timer = Erebus.TimerOutput()
        acc = SimulationAccumulators(
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            cfg_save.geometry.rplanet,
            0.0,
            0,
            0.0,
            0.0,
            coords.xcenter,
            coords.ycenter,
            0.0,
            nothing,
            nothing,
            nothing,
            nothing,
        )
        state = SimulationState(
            grids,
            markers,
            acc,
            Any[],
            nothing,
            Random.MersenneTwister(42),
            timer,
            1,
            100.0,
            100.0,
        )

        mktempdir() do tmpdir
            ckpt_path = joinpath(tmpdir, "checkpoint_test.jld2")
            Erebus.save_state(ckpt_path, state, coords, cfg_save)

            # Reload simulation state with a restart config where percolation_active = true
            cfg_restart = SimulationConfig(
                coreformation=CoreFormationConfig(percolation_active=true, Xfe_bulk=0.325),
                metal_partition=MetalPartitionConfig(active=false),
            )
            state_loaded, coords_loaded, cfg_loaded = Erebus.load_simulation_state(
                ckpt_path; cfg=cfg_restart, force_restart_config=true
            )

            @test haskey(state_loaded.markers.groups, :metal)
            @test cfg_loaded === cfg_restart
            rock_idx = findfirst(==(1), state_loaded.markers.core.tm)
            @test state_loaded.markers.Xfe_bulk[rock_idx] ≈ 0.325
            @test all(iszero, state_loaded.markers.Xfem)
        end
    end
end

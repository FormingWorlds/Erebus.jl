using Test
using Erebus
using JLD2: JLD2

@testset "Flagship Tutorial: Growth to Lunar Mass" begin
    config_path = normpath(
        joinpath(@__DIR__, "..", "configs", "lunar_growth_tutorial.toml")
    )

    # ---------------------------------------------------------------------
    # 1. Configuration Validation and Schema Compliance
    # ---------------------------------------------------------------------
    @testset "Tutorial Configuration Loading & Validation" begin
        @test isfile(config_path)
        cfg = load_config(config_path)

        # Multi-stage accretion parameters
        @test cfg.accretion.active
        @test cfg.accretion.mode === :multistage
        @test cfg.accretion.stage1_mode === :safronov
        @test cfg.accretion.stage2_mode === :pebble_auto
        @test cfg.accretion.stage3_mode === :safronov
        @test isapprox(cfg.accretion.M_initial, 1.753e18, rtol=1e-12)
        @test isapprox(cfg.accretion.M_target, 7.35e22, rtol=1e-12)
        @test isapprox(cfg.accretion.R_initial, 50000.0, rtol=1e-12)
        @test isapprox(cfg.accretion.R_target, 1737000.0, rtol=1e-12)
        @test isapprox(cfg.accretion.rho_bulk, 3348.0, rtol=1e-12)

        # Telescoping grid expansion
        @test cfg.telescoping.active
        @test cfg.telescoping.max_telescope_levels == 6
        @test isapprox(cfg.telescoping.r_threshold_fraction, 0.70, rtol=1e-12)
        @test isapprox(cfg.telescoping.target_radius, 1737000.0, rtol=1e-12)

        # Protoplanetary disk gas dispersal
        @test cfg.disk.enabled
        @test cfg.disk.dispersal_active
        @test isapprox(cfg.disk.t_dispersal_myr, 2.0, rtol=1e-12)
        @test isapprox(cfg.disk.dt_dispersal_myr, 0.10, rtol=1e-12)

        # Coupled atmosphere, core formation, and melting
        @test cfg.atmosphere.active
        @test cfg.atmosphere.mode === :guillot
        @test cfg.coreformation.percolation_active
        @test cfg.coreformation.settling_active
        @test isapprox(cfg.coreformation.Xfe_bulk, 0.015, rtol=1e-12)
        @test cfg.melting.active

        # Full validation contract
        @test validate_config(cfg) === nothing
    end

    # ---------------------------------------------------------------------
    # 2. Coupled Physical Transition Dynamics
    # ---------------------------------------------------------------------
    @testset "Physical Transition & Accretion Dispatch" begin
        cfg = load_config(config_path)
        sec_yr = 3.15576e7
        Omega_K_1au = compute_keplerian_frequency(1.495978707e11, 1.98847e30)

        # Initial seed at t = 0 (Stage 1 Safronov)
        rate_seed = compute_accretion_rate(
            0.0, cfg.accretion.M_initial, cfg.accretion.R_initial, cfg.accretion, cfg.disk
        )
        expected_seed = compute_safronov_accretion_rate(
            cfg.accretion.M_initial,
            cfg.accretion.R_initial,
            cfg.accretion.Sigma_pl_0,
            cfg.accretion.v_disp_kms * 1000.0,
            Omega_K_1au,
        )
        @test isapprox(rate_seed, expected_seed, rtol=1e-10)
        @test rate_seed > 0.0

        # Intermediate mass embryo during gas disk phase (Stage 2 Pebble Accretion)
        M_inter = 1.0e22
        R_inter = 890000.0
        t_early = 0.5 * 1.0e6 * sec_yr
        rate_peb_early = compute_accretion_rate(
            t_early, M_inter, R_inter, cfg.accretion, cfg.disk
        )
        @test rate_peb_early > rate_seed
        @test rate_peb_early > 0.0

        # Post-dispersal at t = 3.0 Ma: pebble accretion halted, evaluates to Stage 3
        t_late = 3.0 * 1.0e6 * sec_yr
        rate_late = compute_accretion_rate(
            t_late, M_inter, R_inter, cfg.accretion, cfg.disk
        )
        expected_late_stage3 = compute_safronov_accretion_rate(
            M_inter,
            R_inter,
            cfg.accretion.Sigma_pl_0,
            cfg.accretion.v_disp_kms * 1000.0,
            Omega_K_1au,
        )
        @test isapprox(rate_late, expected_late_stage3, rtol=1e-10)
        @test rate_late > 0.0
    end

    # ---------------------------------------------------------------------
    # 3. Demonstration Simulation Execution & Checkpoints
    # ---------------------------------------------------------------------
    @testset "Demonstration Simulation Run" begin
        mktempdir() do tmp_dir
            cfg_run = load_config(config_path)
            # Configure minimal 2-step execution for test verification
            cfg_demo = SimulationConfig(
                grid=GridConfig(Nx=33, Ny=33, xsize=140000.0, ysize=140000.0),
                geometry=cfg_run.geometry,
                time=TimeConfig(
                    dt_initial=10.0,
                    dt_longest=10.0,
                    start_time=0.0,
                    endtime=1.0e6,
                    start_step=1,
                    n_steps=2,
                ),
                solver=SolverConfig(titermax=2, nplast=2),
                poroelasticity=cfg_run.poroelasticity,
                thermodynamics=cfg_run.thermodynamics,
                materials=cfg_run.materials,
                accretion=cfg_run.accretion,
                telescoping=cfg_run.telescoping,
                disk=cfg_run.disk,
                melting=cfg_run.melting,
                atmosphere=cfg_run.atmosphere,
                volatiles=cfg_run.volatiles,
                retention=cfg_run.retention,
                coreformation=cfg_run.coreformation,
                output=OutputConfig(output_dir=tmp_dir, savematstep=1, visstep=1),
            )

            # Validate demo configuration
            @test validate_config(cfg_demo) === nothing

            # Execute simulation loop
            simulation_loop(cfg_demo; output_path=tmp_dir)

            chk1 = joinpath(tmp_dir, "output_00001.jld2")
            chk2 = joinpath(tmp_dir, "output_00002.jld2")
            @test isfile(chk1)
            @test isfile(chk2)

            data2 = JLD2.load(chk2)
            @test data2["telescope_level"] == 1
            @test data2["Nx"] == 65
            @test data2["Ny"] == 65
            @test isapprox(data2["xsize"], 280000.0, rtol=1e-12)
            @test data2["rplanet"] >= cfg_run.geometry.rplanet
            @test all(isfinite, data2["tk2"])
            @test all(isfinite, data2["pr"])
        end
    end
end

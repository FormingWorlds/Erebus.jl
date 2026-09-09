using Test
using Erebus
using Erebus.Config
using Erebus.Physics
using Erebus.Particles
using Erebus.Numerics
using Erebus.Simulation
using LinearAlgebra
using StaticArrays
using TOML
using JLD2

@testset "Metal-Silicate Volatile Partitioning & Core Segregation Transport" begin
    @testset "MetalPartitionConfig Schema & Defaults" begin
        cfg_def = MetalPartitionConfig()
        @test !cfg_def.active
        @test cfg_def.model_carbon === :grewal2019
        @test cfg_def.model_nitrogen === :grewal2019
        @test cfg_def.model_hydrogen === :clesi2018
        @test cfg_def.model_sulfur === :boujibar2014
        @test isapprox(cfg_def.D_H_const, 0.5; atol=1.0e-12)
        @test isapprox(cfg_def.D_C_const, 500.0; atol=1.0e-12)
        @test isapprox(cfg_def.D_N_const, 20.0; atol=1.0e-12)
        @test isapprox(cfg_def.D_S_const, 200.0; atol=1.0e-12)
        @test isapprox(cfg_def.equilibration_rate, 1.0; atol=1.0e-12)
        @test cfg_def.dynamic_sulfur_density
        @test isapprox(cfg_def.D_min, 1.0e-4; atol=1.0e-12)
        @test isapprox(cfg_def.D_max, 1.0e5; atol=1.0e-12)
        @test isapprox(cfg_def.initial_metal_h_ppm, 0.0; atol=1.0e-12)
        @test isapprox(cfg_def.initial_metal_c_ppm, 0.0; atol=1.0e-12)
        @test isapprox(cfg_def.initial_metal_n_ppm, 0.0; atol=1.0e-12)
        @test isapprox(cfg_def.initial_metal_s_ppm, 0.0; atol=1.0e-12)
        @test isapprox(cfg_def.core_radius_fraction, 0.5; atol=1.0e-12)
        @test isapprox(cfg_def.phi_core_threshold, 0.40; atol=1.0e-12)

        # SimulationConfig inclusion
        sim_cfg = SimulationConfig()
        @test sim_cfg.metal_partition isa MetalPartitionConfig
        @test !sim_cfg.metal_partition.active
    end

    @testset "MetalPartitionConfig Parameter Validation" begin
        # Active requires volatiles.active
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                volatiles=VolatilesConfig(; active=false),
                metal_partition=MetalPartitionConfig(; active=true),
            ),
        )

        # Valid active configuration
        valid_cfg = SimulationConfig(;
            volatiles=VolatilesConfig(; active=true),
            metal_partition=MetalPartitionConfig(; active=true),
        )
        @test validate_config(valid_cfg) === nothing

        # Invalid D_min
        @test_throws ArgumentError validate_config(
            SimulationConfig(; metal_partition=MetalPartitionConfig(; D_min=-1.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; metal_partition=MetalPartitionConfig(; D_min=0.0))
        )

        # Invalid D_max (< D_min)
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                metal_partition=MetalPartitionConfig(; D_min=10.0, D_max=5.0)
            ),
        )

        # Invalid equilibration_rate
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                metal_partition=MetalPartitionConfig(; equilibration_rate=-0.1)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                metal_partition=MetalPartitionConfig(; equilibration_rate=1.5)
            ),
        )

        # Invalid core thresholds
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                metal_partition=MetalPartitionConfig(; core_radius_fraction=-0.1)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                metal_partition=MetalPartitionConfig(; phi_core_threshold=1.2)
            ),
        )

        # Invalid constant partition coefficients
        @test_throws ArgumentError validate_config(
            SimulationConfig(; metal_partition=MetalPartitionConfig(; D_C_const=-5.0))
        )

        # Invalid models
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                metal_partition=MetalPartitionConfig(; model_carbon=:invalid_model)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                metal_partition=MetalPartitionConfig(; model_nitrogen=:invalid_model)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                metal_partition=MetalPartitionConfig(; model_hydrogen=:invalid_model)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                metal_partition=MetalPartitionConfig(; model_sulfur=:invalid_model)
            ),
        )
    end

    @testset "MetalPartitionConfig TOML Serialization Round-Trip" begin
        custom_cfg = SimulationConfig(;
            volatiles=VolatilesConfig(; active=true),
            metal_partition=MetalPartitionConfig(;
                active=true,
                model_carbon=:fischer2020,
                model_nitrogen=:grewal2019,
                model_hydrogen=:clesi2018,
                model_sulfur=:boujibar2014,
                D_C_const=450.0,
                equilibration_rate=0.75,
                core_radius_fraction=0.45,
                phi_core_threshold=0.35,
            ),
        )
        buf = IOBuffer()
        save_config(buf, custom_cfg)
        toml_str = String(take!(buf))
        @test occursin("[metal_partition]", toml_str)
        @test occursin("model_carbon = \"fischer2020\"", toml_str)

        loaded_cfg = load_config(toml_str)
        @test loaded_cfg.metal_partition.active
        @test loaded_cfg.metal_partition.model_carbon === :fischer2020
        @test isapprox(loaded_cfg.metal_partition.D_C_const, 450.0; atol=1.0e-12)
        @test isapprox(loaded_cfg.metal_partition.equilibration_rate, 0.75; atol=1.0e-12)
        @test isapprox(loaded_cfg.metal_partition.core_radius_fraction, 0.45; atol=1.0e-12)
        @test isapprox(loaded_cfg.metal_partition.phi_core_threshold, 0.35; atol=1.0e-12)
    end

    @testset "Partition Coefficient Parameterizations & Invariants" begin
        # 1. Type stability and baseline inference
        T_ref = 1800.0
        P_ref = 5.0e7 # 0.5 kbar planetesimal core
        dIW_ref = -2.0
        wS_fe = 0.0 # pure Fe metal
        wS_eutectic = 0.31 # Fe-FeS eutectic

        @test @inferred(
            compute_metal_silicate_partition_coefficient(:C, T_ref, P_ref, dIW_ref, wS_fe)
        ) isa Float64
        @test @inferred(
            compute_metal_silicate_partition_coefficient(:N, T_ref, P_ref, dIW_ref, wS_fe)
        ) isa Float64
        @test @inferred(
            compute_metal_silicate_partition_coefficient(:H, T_ref, P_ref, dIW_ref, wS_fe)
        ) isa Float64
        @test @inferred(
            compute_metal_silicate_partition_coefficient(:S, T_ref, P_ref, dIW_ref, wS_fe)
        ) isa Float64

        # 2. Domain guards
        @test_throws DomainError compute_metal_silicate_partition_coefficient(
            :C, 0.0, P_ref, dIW_ref, wS_fe
        )
        @test_throws DomainError compute_metal_silicate_partition_coefficient(
            :C, -100.0, P_ref, dIW_ref, wS_fe
        )
        @test_throws DomainError compute_metal_silicate_partition_coefficient(
            :C, T_ref, -1.0, dIW_ref, wS_fe
        )
        @test_throws DomainError compute_metal_silicate_partition_coefficient(
            :C, T_ref, P_ref, NaN, wS_fe
        )
        @test_throws DomainError compute_metal_silicate_partition_coefficient(
            :C, T_ref, P_ref, dIW_ref, -0.05
        )
        @test_throws DomainError compute_metal_silicate_partition_coefficient(
            :C, T_ref, P_ref, dIW_ref, 1.05
        )
        @test_throws DomainError compute_metal_silicate_partition_coefficient(
            :C, T_ref, P_ref, dIW_ref, wS_fe; D_min=-1.0
        )
        @test_throws DomainError compute_metal_silicate_partition_coefficient(
            :C, T_ref, P_ref, dIW_ref, wS_fe; D_min=10.0, D_max=5.0
        )
        @test_throws ArgumentError compute_metal_silicate_partition_coefficient(
            :UnknownSpecies, T_ref, P_ref, dIW_ref, wS_fe
        )
        @test_throws ArgumentError compute_metal_silicate_partition_coefficient(
            :C, T_ref, P_ref, dIW_ref, wS_fe; model=:nonexistent
        )

        # 3. Carbon partition suppression by dissolved sulfur (Grewal et al. 2019b)
        D_C_fe = compute_metal_silicate_partition_coefficient(
            :C, T_ref, P_ref, dIW_ref, wS_fe; model=:grewal2019
        )
        D_C_eutectic = compute_metal_silicate_partition_coefficient(
            :C, T_ref, P_ref, dIW_ref, wS_eutectic; model=:grewal2019
        )
        @test D_C_fe > 1000.0 # Strongly siderophile in pure metallic iron
        @test D_C_eutectic < 50.0 # Strongly suppressed by sulfur at Fe-FeS eutectic
        @test D_C_fe / D_C_eutectic > 50.0

        # 4. Nitrogen partition stability relative to carbon (Grewal et al. 2019a, 2019b)
        D_N_fe = compute_metal_silicate_partition_coefficient(
            :N, T_ref, P_ref, dIW_ref, wS_fe; model=:grewal2019
        )
        D_N_eutectic = compute_metal_silicate_partition_coefficient(
            :N, T_ref, P_ref, dIW_ref, wS_eutectic; model=:grewal2019
        )
        @test D_N_fe > 10.0
        @test D_N_eutectic > 10.0
        @test (D_N_fe / D_N_eutectic) < 5.0 # Nitrogen is much less sensitive to sulfur than carbon
        @test (D_C_fe / D_N_fe) > 20.0 # High C/N partition in pure iron
        @test (D_C_eutectic / D_N_eutectic) < 2.0 # C/N partition drops near unity in sulfur-rich metal

        # 5. Hydrogen partitioning (Clesi et al. 2018)
        D_H_val = compute_metal_silicate_partition_coefficient(
            :H, T_ref, P_ref, dIW_ref, wS_fe; model=:clesi2018
        )
        @test 0.1 <= D_H_val <= 2.0 # Moderately siderophile to slightly lithophile

        # 6. Sulfur partitioning (Boujibar et al. 2014)
        D_S_val = compute_metal_silicate_partition_coefficient(
            :S, T_ref, P_ref, dIW_ref, wS_fe; model=:boujibar2014
        )
        @test D_S_val > 50.0 # Highly chalcophile/siderophile

        # 7. Fischer et al. (2020) carbon model
        D_C_fischer = compute_metal_silicate_partition_coefficient(
            :C, T_ref, P_ref, dIW_ref, wS_fe; model=:fischer2020
        )
        @test D_C_fischer > 500.0
        @test isapprox(
            compute_metal_silicate_partition_coefficient(
                :C, T_ref, P_ref, dIW_ref, wS_eutectic; model=:fischer2020
            ),
            D_C_fischer;
            atol=1.0e-12,
        )

        # 8. Constant model and clamping bounds
        D_c_test = compute_metal_silicate_partition_coefficient(
            :C, T_ref, P_ref, dIW_ref, wS_fe; model=:constant, D_const=777.0
        )
        @test isapprox(D_c_test, 777.0; atol=1.0e-12)
        D_c_clamp_hi = compute_metal_silicate_partition_coefficient(
            :C, T_ref, P_ref, dIW_ref, wS_fe; model=:constant, D_const=1.0e8, D_max=1.0e4
        )
        @test isapprox(D_c_clamp_hi, 1.0e4; atol=1.0e-12)
        D_c_clamp_lo = compute_metal_silicate_partition_coefficient(
            :C, T_ref, P_ref, dIW_ref, wS_fe; model=:constant, D_const=1.0e-8, D_min=1.0e-3
        )
        @test isapprox(D_c_clamp_lo, 1.0e-3; atol=1.0e-12)

        # 9. Batch coefficients function
        cfg_part = MetalPartitionConfig()
        batch = @inferred(
            compute_metal_silicate_partition_coefficients(
                T_ref, P_ref, dIW_ref, wS_fe, cfg_part
            )
        )
        @test isapprox(batch.D_C, D_C_fe; atol=1.0e-12)
        @test isapprox(batch.D_N, D_N_fe; atol=1.0e-12)
        @test isapprox(batch.D_H, D_H_val; atol=1.0e-12)
        @test isapprox(batch.D_S, D_S_val; atol=1.0e-12)
    end

    @testset "Setup Marker Metal Volatile Properties" begin
        nmarks = 50
        (Xfe_H_m, Xfe_C_m, Xfe_N_m, Xfe_S_m) = setup_marker_metal_volatile_properties(
            nmarks;
            initial_h_ppm=12.0,
            initial_c_ppm=150.0,
            initial_n_ppm=25.0,
            initial_s_ppm=310000.0,
        )
        @test length(Xfe_H_m) == nmarks
        @test length(Xfe_C_m) == nmarks
        @test length(Xfe_N_m) == nmarks
        @test length(Xfe_S_m) == nmarks
        @test all(x -> isapprox(x, 12.0; atol=1.0e-12), Xfe_H_m)
        @test all(x -> isapprox(x, 150.0; atol=1.0e-12), Xfe_C_m)
        @test all(x -> isapprox(x, 25.0; atol=1.0e-12), Xfe_N_m)
        @test all(x -> isapprox(x, 310000.0; atol=1.0e-12), Xfe_S_m)
    end

    @testset "Equilibrate Metal-Silicate Volatiles Mass Conservation" begin
        # Set up a single marker with silicate melt and molten metal
        cfg = MetalPartitionConfig(;
            active=true,
            model_carbon=:grewal2019,
            model_nitrogen=:grewal2019,
            model_hydrogen=:clesi2018,
            model_sulfur=:boujibar2014,
            equilibration_rate=1.0,
        )

        nmarks = 1
        m = 1
        rho_sil = 3300.0
        rho_met = 7000.0
        F_fe = 0.80 # 80% molten metal
        F_melt = 0.40 # 40% molten silicate
        phi_fe_bulk = 0.20 # 20 vol% bulk metal
        phi_fe_liq = phi_fe_bulk * F_fe

        T_val = 1750.0
        P_val = 1.0e8
        fO2_val = -2.0

        # Initial concentrations
        init_H2O_wtpct = 1.5
        init_C_ppm = 120.0
        init_N_ppm = 35.0
        init_S_ppm = 800.0

        init_fe_H_ppm = 0.0
        init_fe_C_ppm = 0.0
        init_fe_N_ppm = 0.0
        init_fe_S_ppm = 0.0

        Xfe_bulk = [phi_fe_bulk]
        Xfem = [phi_fe_liq]
        XH2Om = [init_H2O_wtpct]
        XCm = [init_C_ppm]
        XNm = [init_N_ppm]
        XSm = [init_S_ppm]
        Xfe_H_m = [init_fe_H_ppm]
        Xfe_C_m = [init_fe_C_ppm]
        Xfe_N_m = [init_fe_N_ppm]
        Xfe_S_m = [init_fe_S_ppm]

        # Stoichiometric conversion for H: 1 wt% H2O = f_H ppmw H
        f_H = (2.0 * 1.00794 / 18.01528) * 1.0e4
        m_sil = (1.0 - phi_fe_bulk) * rho_sil
        m_met = phi_fe_bulk * F_fe * rho_met

        # Initial elemental volatile masses in marker
        M_H_init = m_sil * (init_H2O_wtpct * f_H) + m_met * init_fe_H_ppm
        M_C_init = m_sil * init_C_ppm + m_met * init_fe_C_ppm
        M_N_init = m_sil * init_N_ppm + m_met * init_fe_N_ppm
        M_S_init = m_sil * init_S_ppm + m_met * init_fe_S_ppm

        # Perform equilibration
        equilibrate_metal_silicate_volatiles!(
            m,
            F_fe,
            F_melt,
            T_val,
            P_val,
            fO2_val,
            Xfe_bulk,
            Xfem,
            XH2Om,
            XCm,
            XNm,
            XSm,
            Xfe_H_m,
            Xfe_C_m,
            Xfe_N_m,
            Xfe_S_m,
            cfg;
            rho_silicate=rho_sil,
            rho_metal=rho_met,
            equilibration_fraction=cfg.equilibration_rate,
        )

        # Final elemental volatile masses in marker
        M_H_final = m_sil * (XH2Om[1] * f_H) + m_met * Xfe_H_m[1]
        M_C_final = m_sil * XCm[1] + m_met * Xfe_C_m[1]
        M_N_final = m_sil * XNm[1] + m_met * Xfe_N_m[1]
        M_S_final = m_sil * XSm[1] + m_met * Xfe_S_m[1]

        # Verify exact mass conservation
        @test isapprox(M_H_final, M_H_init; rtol=1.0e-12)
        @test isapprox(M_C_final, M_C_init; rtol=1.0e-12)
        @test isapprox(M_N_final, M_N_init; rtol=1.0e-12)
        @test isapprox(M_S_final, M_S_init; rtol=1.0e-12)

        # Volatiles must have partitioned into metallic phase
        @test Xfe_C_m[1] > 0.0
        @test Xfe_N_m[1] > 0.0
        @test Xfe_S_m[1] > 0.0
        @test Xfe_H_m[1] > 0.0

        # Silicate concentrations must have dropped accordingly
        @test XCm[1] < init_C_ppm
        @test XNm[1] < init_N_ppm
        @test XSm[1] < init_S_ppm
        @test XH2Om[1] < init_H2O_wtpct

        # Test subsolidus no-op: when F_melt == 0, equilibration does not modify inventories
        C_c_prev = XCm[1]
        C_c_met_prev = Xfe_C_m[1]
        equilibrate_metal_silicate_volatiles!(
            m,
            F_fe,
            0.0, # 0% silicate melt
            T_val,
            P_val,
            fO2_val,
            Xfe_bulk,
            Xfem,
            XH2Om,
            XCm,
            XNm,
            XSm,
            Xfe_H_m,
            Xfe_C_m,
            Xfe_N_m,
            Xfe_S_m,
            cfg;
            rho_silicate=rho_sil,
            rho_metal=rho_met,
            equilibration_fraction=1.0,
        )
        @test isapprox(XCm[1], C_c_prev; atol=1.0e-12)
        @test isapprox(Xfe_C_m[1], C_c_met_prev; atol=1.0e-12)

        # Test partial kinetic equilibration (alpha_eq = 0.5)
        XH2Om[1] = init_H2O_wtpct
        XCm[1] = init_C_ppm
        XNm[1] = init_N_ppm
        XSm[1] = init_S_ppm
        Xfe_H_m[1] = 0.0
        Xfe_C_m[1] = 0.0
        Xfe_N_m[1] = 0.0
        Xfe_S_m[1] = 0.0

        equilibrate_metal_silicate_volatiles!(
            m,
            F_fe,
            F_melt,
            T_val,
            P_val,
            fO2_val,
            Xfe_bulk,
            Xfem,
            XH2Om,
            XCm,
            XNm,
            XSm,
            Xfe_H_m,
            Xfe_C_m,
            Xfe_N_m,
            Xfe_S_m,
            cfg;
            rho_silicate=rho_sil,
            rho_metal=rho_met,
            equilibration_fraction=0.5,
        )

        M_C_partial = m_sil * XCm[1] + m_met * Xfe_C_m[1]
        @test isapprox(M_C_partial, M_C_init; rtol=1.0e-12)
        @test 0.0 < Xfe_C_m[1] < (M_C_init / m_met)
    end

    @testset "Advective Volatile Segregation & Conservation in apply_metal_segregation!" begin
        cfg_sim = default_config()
        coords = GridCoordinates(cfg_sim.grid)
        Ny_val, Nx_val = coords.Ny, coords.Nx

        cfg_core = CoreFormationConfig(;
            percolation_active=true,
            settling_active=true,
            phi_crit_perc=0.05,
            phi_residual=0.01,
            sulfur_fraction=0.05,
            metal_density_mode=:sanloup2000,
        )
        cfg_part = MetalPartitionConfig(; active=true, dynamic_sulfur_density=true)

        nmarks = coords.start_marknum
        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0, Fm) = setup_marker_properties(
            nmarks, coords
        )
        Xfem, Xfem0, Xfe_bulk = setup_marker_metal_properties(nmarks)
        (Xfe_H_m, Xfe_C_m, Xfe_N_m, Xfe_S_m) = setup_marker_metal_volatile_properties(
            nmarks;
            initial_h_ppm=15.0,
            initial_c_ppm=400.0,
            initial_n_ppm=40.0,
            initial_s_ppm=50000.0,
        )

        rplanet = cfg_sim.geometry.rplanet
        xc = cfg_sim.geometry.xcenter
        yc = cfg_sim.geometry.ycenter

        for i in 1:nmarks
            frac = i / nmarks
            r_val = 0.6 * rplanet * sqrt(frac)
            theta = 2.0 * pi * frac
            xm[i] = xc + r_val * cos(theta)
            ym[i] = yc + r_val * sin(theta)
            tm[i] = 1 # rock marker
            tkm[i] = 1600.0 # hot
            Fm[i] = 0.50 # silicate melt
            Xfe_bulk[i] = 0.25
            Xfem[i] = 0.25
            Xfe_H_m[i] = 15.0
            Xfe_C_m[i] = 400.0
            Xfe_N_m[i] = 40.0
            Xfe_S_m[i] = 50000.0
        end

        rho_m = cfg_core.rho_metal
        init_total_metal = sum(Xfe_bulk .* rho_m)
        init_total_H = sum(Xfe_bulk .* rho_m .* (Xfe_H_m .* 1.0e-6))
        init_total_C = sum(Xfe_bulk .* rho_m .* (Xfe_C_m .* 1.0e-6))
        init_total_N = sum(Xfe_bulk .* rho_m .* (Xfe_N_m .* 1.0e-6))
        init_total_S = sum(Xfe_bulk .* rho_m .* (Xfe_S_m .* 1.0e-6))

        @test init_total_metal > 0.0
        @test init_total_C > 0.0

        gx = zeros(Float64, Ny_val, Nx_val)
        gy = zeros(Float64, Ny_val, Nx_val)
        for j in 1:Nx_val, i in 1:Ny_val
            dx = coords.x[j] - xc
            dy = coords.y[i] - yc
            rm = max(sqrt(dx^2 + dy^2), 1000.0)
            gx[i, j] = -9.81 * (dx / rm)
            gy[i, j] = -9.81 * (dy / rm)
        end
        ETA = fill(1.0e19, Ny_val - 1, Nx_val - 1)

        dt = 1.0e9

        seg_res = apply_metal_segregation!(
            xm,
            ym,
            tm,
            tkm,
            phim,
            Xfe_bulk,
            Xfem,
            nmarks,
            dt,
            cfg_core;
            coords=coords,
            xcenter=xc,
            ycenter=yc,
            rplanet=rplanet,
            gx=gx,
            gy=gy,
            rho_silicate=3300.0,
            eta_silicate=1.0e19,
            ETA=ETA,
            Fm=Fm,
            T_solidus_silicate=1400.0,
            T_liquidus_silicate=1800.0,
            cfg_partition=cfg_part,
            Xfe_H_m=Xfe_H_m,
            Xfe_C_m=Xfe_C_m,
            Xfe_N_m=Xfe_N_m,
            Xfe_S_m=Xfe_S_m,
        )

        @test seg_res.max_v_seg >= 0.0

        final_total_metal = sum(Xfe_bulk .* rho_m)
        final_total_H = sum(Xfe_bulk .* rho_m .* (Xfe_H_m .* 1.0e-6))
        final_total_C = sum(Xfe_bulk .* rho_m .* (Xfe_C_m .* 1.0e-6))
        final_total_N = sum(Xfe_bulk .* rho_m .* (Xfe_N_m .* 1.0e-6))
        final_total_S = sum(Xfe_bulk .* rho_m .* (Xfe_S_m .* 1.0e-6))

        @test isapprox(final_total_metal, init_total_metal; rtol=1.0e-11)
        @test isapprox(final_total_H, init_total_H; rtol=1.0e-11)
        @test isapprox(final_total_C, init_total_C; rtol=1.0e-11)
        @test isapprox(final_total_N, init_total_N; rtol=1.0e-11)
        @test isapprox(final_total_S, init_total_S; rtol=1.0e-11)

        @test all(x -> x >= 0.0, Xfe_H_m)
        @test all(x -> x >= 0.0, Xfe_C_m)
        @test all(x -> x >= 0.0, Xfe_N_m)
        @test all(x -> x >= 0.0, Xfe_S_m)
    end

    @testset "Compute Core Volatile Budgets & Geochemical Benchmark" begin
        nmarks = 100
        xc = 70000.0
        yc = 70000.0
        rplanet = 50000.0
        xm = fill(xc, nmarks)
        ym = fill(yc, nmarks)
        tm = fill(1, nmarks)

        for i in 1:50
            xm[i] = xc + 3000.0 * (i / 50)
            ym[i] = yc
        end
        for i in 51:100
            xm[i] = xc + 40000.0
            ym[i] = yc
        end

        Xfe_bulk = fill(0.50, nmarks)
        for i in 51:100
            Xfe_bulk[i] = 0.05
        end

        Xfe_H_m = fill(5.0, nmarks)
        Xfe_C_m = fill(350.0, nmarks)
        Xfe_N_m = fill(25.0, nmarks)
        Xfe_S_m = fill(20000.0, nmarks)

        budgets = @inferred(
            compute_core_volatile_budgets(
                xm,
                ym,
                tm,
                Xfe_bulk,
                Xfe_H_m,
                Xfe_C_m,
                Xfe_N_m,
                Xfe_S_m,
                nmarks;
                xcenter=xc,
                ycenter=yc,
                rplanet=rplanet,
                rho_metal=7000.0,
                core_radius_fraction=0.4,
                phi_core_threshold=0.40,
            )
        )

        @test budgets.M_core_metal > 0.0
        @test budgets.M_core_C > 0.0
        @test budgets.M_core_N > 0.0
        @test budgets.M_core_S > 0.0
        @test budgets.M_core_H > 0.0

        @test isapprox(budgets.w_core_C_ppm, 350.0; atol=1.0e-10)
        @test isapprox(budgets.w_core_N_ppm, 25.0; atol=1.0e-10)
        @test isapprox(budgets.w_core_H_ppm, 5.0; atol=1.0e-10)
        @test isapprox(budgets.w_core_S_ppm, 20000.0; atol=1.0e-10)
        @test isapprox(budgets.w_core_S_wtpct, 2.0; atol=1.0e-10)

        @test 20.0 <= budgets.w_core_C_ppm <= 5000.0
        @test 5.0 <= budgets.w_core_N_ppm <= 150.0
        @test 0.1 <= budgets.w_core_S_wtpct <= 10.0
        @test 0.1 <= budgets.w_core_H_ppm <= 50.0
    end

    @testset "Simulation Loop & Checkpoint Round-Trip with Metal Partitioning" begin
        output_dir = mktempdir()
        try
            cfg = SimulationConfig(;
                grid=GridConfig(; Nx=8, Ny=8, xsize=70000.0, ysize=70000.0),
                time=TimeConfig(; n_steps=2, dt_initial=1.0e8),
                thermodynamics=ThermalConfig(; hr_al=true, hr_fe=true),
                melting=MeltingConfig(; active=true),
                coreformation=CoreFormationConfig(;
                    percolation_active=true, settling_active=true, sulfur_fraction=0.10
                ),
                volatiles=VolatilesConfig(;
                    active=true,
                    initial_water_wtpct=0.5,
                    initial_carbon_ppm=100.0,
                    initial_nitrogen_ppm=20.0,
                    initial_sulfur_ppm=500.0,
                ),
                metal_partition=MetalPartitionConfig(;
                    active=true,
                    model_carbon=:grewal2019,
                    model_nitrogen=:grewal2019,
                    initial_metal_c_ppm=10.0,
                    initial_metal_n_ppm=2.0,
                    initial_metal_s_ppm=10000.0,
                ),
                output=OutputConfig(; output_dir=output_dir, savematstep=1),
            )

            @test run_simulation(cfg) === nothing

            ckpt1 = joinpath(output_dir, "output_00001.jld2")
            ckpt2 = joinpath(output_dir, "output_00002.jld2")
            @test isfile(ckpt1)
            @test isfile(ckpt2)

            data2 = load_state(ckpt2)
            @test haskey(data2, "Xfe_H_m")
            @test haskey(data2, "Xfe_C_m")
            @test haskey(data2, "Xfe_N_m")
            @test haskey(data2, "Xfe_S_m")
            @test haskey(data2, "core_budgets")
            @test data2["timestep"] == 2

            # Resume from step 1 for 1 step
            cfg_resume = SimulationConfig(;
                grid=cfg.grid,
                time=TimeConfig(; n_steps=2, start_step=2, dt_initial=1.0e8),
                thermodynamics=cfg.thermodynamics,
                melting=cfg.melting,
                coreformation=cfg.coreformation,
                volatiles=cfg.volatiles,
                metal_partition=cfg.metal_partition,
                output=OutputConfig(;
                    output_dir=output_dir, restart_from=ckpt1, savematstep=1
                ),
            )
            @test run_simulation(cfg_resume) === nothing
            @test isfile(ckpt2)
        finally
            rm(output_dir; recursive=true, force=true)
        end
    end
end

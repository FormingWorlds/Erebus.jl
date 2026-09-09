using Test
using Erebus
using StaticArrays
using TOML

@testset "Volatile Retention Floors & Venting Drainage Coupling" begin
    @testset "RetentionConfig Defaults & Type Stability" begin
        cfg = RetentionConfig()
        @test !cfg.active
        @test isapprox(cfg.h2o_retention_ppm, 50.0; atol=1.0e-12)
        @test isapprox(cfg.carbon_retention_ppm, 50.0; atol=1.0e-12)
        @test isapprox(cfg.nitrogen_retention_ppm, 5.0; atol=1.0e-12)
        @test isapprox(cfg.sulfur_retention_ppm, 100.0; atol=1.0e-12)
        @test isapprox(cfg.T_solidus_ref, 1400.0; atol=1.0e-12)
        @test isapprox(cfg.dT_retention, 200.0; atol=1.0e-12)
        @test cfg.retention_law === :nams_exponential
        @test cfg.venting_drainage_active
        @test isapprox(cfg.chi_vent, 1.0; atol=1.0e-12)

        sim_cfg = SimulationConfig()
        @test !sim_cfg.retention.active
        @test isapprox(sim_cfg.volatiles.initial_nitrogen_ppm, 50.0; atol=1.0e-12)
    end

    @testset "RetentionConfig Parameter Bounds Validation" begin
        # Valid active configuration
        valid_cfg = SimulationConfig(;
            retention=RetentionConfig(;
                active=true,
                h2o_retention_ppm=100.0,
                carbon_retention_ppm=75.0,
                nitrogen_retention_ppm=10.0,
                sulfur_retention_ppm=150.0,
                T_solidus_ref=1350.0,
                dT_retention=250.0,
                retention_law=:constant_floor,
                chi_vent=0.8,
            ),
        )
        @test validate_config(valid_cfg) === nothing

        # Invalid negative H2O retention floor
        @test_throws ArgumentError validate_config(
            SimulationConfig(; retention=RetentionConfig(; h2o_retention_ppm=-1.0))
        )
        # Invalid non-finite carbon retention floor
        @test_throws ArgumentError validate_config(
            SimulationConfig(; retention=RetentionConfig(; carbon_retention_ppm=NaN))
        )
        # Invalid negative nitrogen retention floor
        @test_throws ArgumentError validate_config(
            SimulationConfig(; retention=RetentionConfig(; nitrogen_retention_ppm=-5.0))
        )
        # Invalid negative sulfur retention floor
        @test_throws ArgumentError validate_config(
            SimulationConfig(; retention=RetentionConfig(; sulfur_retention_ppm=-10.0))
        )
        # Invalid non-positive solidus reference
        @test_throws ArgumentError validate_config(
            SimulationConfig(; retention=RetentionConfig(; T_solidus_ref=0.0))
        )
        # Invalid non-positive dT_retention
        @test_throws ArgumentError validate_config(
            SimulationConfig(; retention=RetentionConfig(; dT_retention=-50.0))
        )
        # Invalid retention law symbol
        @test_throws ArgumentError validate_config(
            SimulationConfig(; retention=RetentionConfig(; retention_law=:unphysical_law))
        )
        # Invalid chi_vent < 0
        @test_throws ArgumentError validate_config(
            SimulationConfig(; retention=RetentionConfig(; chi_vent=-0.1))
        )
        # Invalid chi_vent > 1
        @test_throws ArgumentError validate_config(
            SimulationConfig(; retention=RetentionConfig(; chi_vent=1.5))
        )
        # Invalid negative initial_nitrogen_ppm
        @test_throws ArgumentError validate_config(
            SimulationConfig(; volatiles=VolatilesConfig(; initial_nitrogen_ppm=-10.0))
        )
    end

    @testset "RetentionConfig TOML Serialization Round-Trip" begin
        orig_cfg = SimulationConfig(;
            retention=RetentionConfig(;
                active=true,
                h2o_retention_ppm=80.0,
                carbon_retention_ppm=60.0,
                nitrogen_retention_ppm=8.0,
                sulfur_retention_ppm=120.0,
                T_solidus_ref=1450.0,
                dT_retention=180.0,
                retention_law=:linear_melt_blend,
                venting_drainage_active=false,
                chi_vent=0.5,
            ),
            volatiles=VolatilesConfig(; initial_nitrogen_ppm=75.0),
        )

        toml_str = save_config(orig_cfg)
        loaded_cfg = load_config(toml_str)

        @test loaded_cfg.retention.active
        @test isapprox(loaded_cfg.retention.h2o_retention_ppm, 80.0; atol=1.0e-12)
        @test isapprox(loaded_cfg.retention.carbon_retention_ppm, 60.0; atol=1.0e-12)
        @test isapprox(loaded_cfg.retention.nitrogen_retention_ppm, 8.0; atol=1.0e-12)
        @test isapprox(loaded_cfg.retention.sulfur_retention_ppm, 120.0; atol=1.0e-12)
        @test isapprox(loaded_cfg.retention.T_solidus_ref, 1450.0; atol=1.0e-12)
        @test isapprox(loaded_cfg.retention.dT_retention, 180.0; atol=1.0e-12)
        @test loaded_cfg.retention.retention_law === :linear_melt_blend
        @test !loaded_cfg.retention.venting_drainage_active
        @test isapprox(loaded_cfg.retention.chi_vent, 0.5; atol=1.0e-12)
        @test isapprox(loaded_cfg.volatiles.initial_nitrogen_ppm, 75.0; atol=1.0e-12)
    end

    @testset "compute_volatile_retention_floor Functional Laws & Asymptotics" begin
        ret_exp = RetentionConfig(;
            active=true,
            h2o_retention_ppm=50.0,
            carbon_retention_ppm=50.0,
            nitrogen_retention_ppm=5.0,
            sulfur_retention_ppm=100.0,
            T_solidus_ref=1400.0,
            dT_retention=200.0,
            retention_law=:nams_exponential,
        )

        # Inactive returns 0.0
        ret_inact = RetentionConfig(; active=false)
        @test isapprox(
            compute_volatile_retention_floor(1300.0, :H2O, ret_inact), 0.0; atol=1.0e-12
        )
        @test isapprox(
            compute_volatile_retention_floor(1300.0, :C, ret_inact), 0.0; atol=1.0e-12
        )

        # Sub-solidus behavior: T <= T_solidus_ref returns exact base floor
        @test isapprox(
            compute_volatile_retention_floor(1000.0, :H2O, ret_exp), 50.0; atol=1.0e-12
        )
        @test isapprox(
            compute_volatile_retention_floor(1400.0, :H2O, ret_exp), 50.0; atol=1.0e-12
        )
        @test isapprox(
            compute_volatile_retention_floor(1400.0, :C, ret_exp), 50.0; atol=1.0e-12
        )
        @test isapprox(
            compute_volatile_retention_floor(1400.0, :N, ret_exp), 5.0; atol=1.0e-12
        )
        @test isapprox(
            compute_volatile_retention_floor(1400.0, :S, ret_exp), 100.0; atol=1.0e-12
        )

        # Super-solidus exponential decay
        c_1600 = compute_volatile_retention_floor(1600.0, :H2O, ret_exp)
        expected_1600 = 50.0 * exp(-(1600.0 - 1400.0) / 200.0) # 50 / e ≈ 18.39397
        @test isapprox(c_1600, expected_1600; atol=1.0e-10)
        @test c_1600 < 50.0
        @test c_1600 > 0.0

        c_1800 = compute_volatile_retention_floor(1800.0, :H2O, ret_exp)
        expected_1800 = 50.0 * exp(-(1800.0 - 1400.0) / 200.0) # 50 / e^2 ≈ 6.76676
        @test isapprox(c_1800, expected_1800; atol=1.0e-10)
        @test c_1800 < c_1600

        # Linear melt blend law
        ret_lin = RetentionConfig(;
            active=true, h2o_retention_ppm=50.0, retention_law=:linear_melt_blend
        )
        @test isapprox(
            compute_volatile_retention_floor(1500.0, :H2O, ret_lin; F_melt=0.0),
            50.0;
            atol=1.0e-12,
        )
        @test isapprox(
            compute_volatile_retention_floor(1500.0, :H2O, ret_lin; F_melt=0.4),
            30.0;
            atol=1.0e-12,
        )
        @test isapprox(
            compute_volatile_retention_floor(1500.0, :H2O, ret_lin; F_melt=1.0),
            0.0;
            atol=1.0e-12,
        )

        # Constant floor law
        ret_const = RetentionConfig(;
            active=true, h2o_retention_ppm=50.0, retention_law=:constant_floor
        )
        @test isapprox(
            compute_volatile_retention_floor(1000.0, :H2O, ret_const), 50.0; atol=1.0e-12
        )
        @test isapprox(
            compute_volatile_retention_floor(2500.0, :H2O, ret_const; F_melt=0.9),
            50.0;
            atol=1.0e-12,
        )

        # Helper functions match dispatcher
        @test isapprox(
            compute_h2o_retention_floor(1500.0, ret_exp),
            compute_volatile_retention_floor(1500.0, :H2O, ret_exp);
            atol=1.0e-12,
        )
        @test isapprox(
            compute_carbon_retention_floor(1500.0, ret_exp),
            compute_volatile_retention_floor(1500.0, :C, ret_exp);
            atol=1.0e-12,
        )
        @test isapprox(
            compute_nitrogen_retention_floor(1500.0, ret_exp),
            compute_volatile_retention_floor(1500.0, :N, ret_exp);
            atol=1.0e-12,
        )
        @test isapprox(
            compute_sulfur_retention_floor(1500.0, ret_exp),
            compute_volatile_retention_floor(1500.0, :S, ret_exp);
            atol=1.0e-12,
        )

        # Unknown species throws
        @test_throws ArgumentError compute_volatile_retention_floor(1500.0, :He, ret_exp)
    end

    @testset "compute_volatile_exsolution Decompression Vacuum Limit & Retention Floor" begin
        # Compare baseline (retention disabled) vs active retention
        # At near-vacuum surface pressure (P = 1.0 Pa), melt solubility S_H2O(P) ≈ 0.
        # Bulk rock water = 1.0 wt% = 0.01 mass fraction. Melt fraction F_m = 0.5.
        P_vac = 1.0
        T_magma = 1500.0
        F_melt = 0.5
        w_H2O_bulk = 0.01 # 1 wt% = 10,000 ppmw
        C_C_bulk = 500.0
        C_N_bulk = 50.0
        C_S_bulk = 1000.0
        delta_IW = -1.0

        # Baseline: without retention floor
        ex_base = compute_volatile_exsolution(
            F_melt,
            P_vac,
            T_magma,
            w_H2O_bulk,
            C_C_bulk,
            C_N_bulk,
            C_S_bulk,
            delta_IW;
            carbon_active=true,
            sulfur_active=true,
            retention_active=false,
        )

        # Near-vacuum: melt capacity is tiny (< 1 ppm), so dissolved volatiles drop to ~0
        @test ex_base.w_H2O_diss < 1.0e-5
        @test ex_base.C_C_diss_ppm < 0.1
        @test ex_base.C_N_diss_ppm < 0.1

        # With active retention floor (50 ppm H2O floor, 50 ppm C floor, 5 ppm N floor, 100 ppm S floor)
        ret_cfg = RetentionConfig(;
            active=true,
            h2o_retention_ppm=50.0,
            carbon_retention_ppm=50.0,
            nitrogen_retention_ppm=5.0,
            sulfur_retention_ppm=100.0,
            T_solidus_ref=1400.0,
            dT_retention=200.0,
            retention_law=:nams_exponential,
        )

        ex_ret = compute_volatile_exsolution(
            F_melt,
            P_vac,
            T_magma,
            w_H2O_bulk,
            C_C_bulk,
            C_N_bulk,
            C_S_bulk,
            delta_IW;
            carbon_active=true,
            sulfur_active=true,
            retention_active=true,
            retention_cfg=ret_cfg,
        )

        # Expected retention floor at 1500 K: 50 * exp(-100/200) = 50 * exp(-0.5) ≈ 30.3265 ppm
        floor_H2O_ppm = compute_h2o_retention_floor(T_magma, ret_cfg; F_melt=F_melt)
        floor_H2O_frac = floor_H2O_ppm * 1.0e-6
        floor_C_ppm = compute_carbon_retention_floor(T_magma, ret_cfg; F_melt=F_melt)
        floor_N_ppm = compute_nitrogen_retention_floor(T_magma, ret_cfg; F_melt=F_melt)
        floor_S_ppm = compute_sulfur_retention_floor(T_magma, ret_cfg; F_melt=F_melt)

        # Retained dissolved concentrations must not drop below the retention floor
        @test ex_ret.w_H2O_diss >= floor_H2O_frac
        @test ex_ret.C_C_diss_ppm >= floor_C_ppm
        @test ex_ret.C_N_diss_ppm >= floor_N_ppm
        @test ex_ret.C_S_diss_ppm >= floor_S_ppm

        # Retained dissolved water is significantly higher than in the unphysical baseline
        @test ex_ret.w_H2O_diss > ex_base.w_H2O_diss
        @test ex_ret.C_C_diss_ppm > ex_base.C_C_diss_ppm
        @test ex_ret.C_N_diss_ppm > ex_base.C_N_diss_ppm

        # Exact mass conservation
        @test isapprox(ex_ret.w_H2O_ex + ex_ret.w_H2O_diss, w_H2O_bulk; atol=1.0e-12)
        @test isapprox(
            ex_ret.w_C_ex + ex_ret.C_C_diss_ppm * 1.0e-6, C_C_bulk * 1.0e-6; atol=1.0e-12
        )
        @test isapprox(
            ex_ret.w_N_ex + ex_ret.C_N_diss_ppm * 1.0e-6, C_N_bulk * 1.0e-6; atol=1.0e-12
        )
        @test isapprox(
            ex_ret.w_S_ex + ex_ret.C_S_diss_ppm * 1.0e-6, C_S_bulk * 1.0e-6; atol=1.0e-12
        )
    end

    @testset "compute_volatile_exsolution Sub-Floor Under-Saturation Invariant" begin
        # When bulk volatile concentration is less than the retention floor,
        # zero exsolution should occur, and 100% of volatiles remain in the rock.
        ret_cfg = RetentionConfig(;
            active=true,
            h2o_retention_ppm=100.0, # 100 ppm floor = 1.0e-4 mass fraction
            carbon_retention_ppm=50.0,
            nitrogen_retention_ppm=10.0,
            sulfur_retention_ppm=100.0,
            retention_law=:constant_floor,
        )

        w_H2O_low = 5.0e-5 # 50 ppm < 100 ppm floor
        C_C_low = 30.0    # 30 ppm < 50 ppm floor
        C_N_low = 5.0     # 5 ppm < 10 ppm floor
        C_S_low = 50.0    # 50 ppm < 100 ppm floor

        ex = compute_volatile_exsolution(
            0.5,
            1.0,
            1400.0,
            w_H2O_low,
            C_C_low,
            C_N_low,
            C_S_low,
            -1.0;
            carbon_active=true,
            sulfur_active=true,
            retention_active=true,
            retention_cfg=ret_cfg,
        )

        # Zero exsolution
        @test isapprox(ex.w_H2O_ex, 0.0; atol=1.0e-14)
        @test isapprox(ex.w_C_ex, 0.0; atol=1.0e-14)
        @test isapprox(ex.w_N_ex, 0.0; atol=1.0e-14)
        @test isapprox(ex.w_S_ex, 0.0; atol=1.0e-14)
        @test isapprox(ex.w_total_ex, 0.0; atol=1.0e-14)

        # 100% retained in rock
        @test isapprox(ex.w_H2O_diss, w_H2O_low; atol=1.0e-14)
        @test isapprox(ex.C_C_diss_ppm, C_C_low; atol=1.0e-12)
        @test isapprox(ex.C_N_diss_ppm, C_N_low; atol=1.0e-12)
        @test isapprox(ex.C_S_diss_ppm, C_S_low; atol=1.0e-12)
    end

    @testset "update_single_marker_volatile_exsolution! Wiring & Porosity Coupling" begin
        cfg_vol = VolatilesConfig(; active=true, carbon_active=true, sulfur_active=true)
        ret_cfg = RetentionConfig(; active=true, h2o_retention_ppm=50.0)

        marknum = 1
        XH2Om = [1.0] # 1 wt%
        XCm = [500.0]
        XNm = [50.0]
        XSm = [1000.0]
        phim = [0.01]

        # Call with retention_cfg active
        dw = update_single_marker_volatile_exsolution!(
            1,
            0.5,
            1.0e6,
            1500.0,
            XH2Om,
            XCm,
            XNm,
            XSm,
            phim,
            cfg_vol;
            retention_cfg=ret_cfg,
        )

        @test dw > 0.0
        @test phim[1] > 0.01 # Exsolution increased porosity
        # Retained water in wt% must be >= floor in wt% (50 ppm = 0.005 wt%)
        @test XH2Om[1] >= 0.005
        @test XCm[1] >= 0.0
        @test XNm[1] >= 0.0
        @test XSm[1] >= 0.0
    end

    @testset "drain_vented_marker_volatiles! Depletion Kinetics & Retention Floor Protection" begin
        coords = GridCoordinates(GridConfig(; Nx=5, Ny=5, xsize=100_000.0, ysize=100_000.0))
        marknum = 10
        xm = fill(50_000.0, marknum)
        ym = fill(50_000.0, marknum)
        tm = fill(2, marknum) # rock markers
        tkm = fill(1400.0, marknum)
        XH2Om = fill(1.0, marknum) # 1.0 wt% = 10,000 ppm
        XCm = fill(500.0, marknum) # 500 ppm
        XNm = fill(50.0, marknum)  # 50 ppm
        XSm = fill(1000.0, marknum) # 1000 ppm

        ret_cfg = RetentionConfig(;
            active=true,
            h2o_retention_ppm=50.0,
            carbon_retention_ppm=50.0,
            nitrogen_retention_ppm=5.0,
            sulfur_retention_ppm=100.0,
            T_solidus_ref=1400.0,
            retention_law=:constant_floor,
            venting_drainage_active=true,
            chi_vent=1.0,
        )

        # Constant venting rate grid: S_vent = 1.0e-11 s^-1 across domain
        S_vent_grid = fill(1.0e-11, coords.Ny1, coords.Nx1)
        dt = 1.0e10 # 10^10 s ~ 317 years, s_vent * dt = 0.1

        # Calculate drainage for one step
        res = drain_vented_marker_volatiles!(
            xm,
            ym,
            tm,
            tkm,
            XH2Om,
            XCm,
            XNm,
            XSm,
            S_vent_grid,
            dt,
            marknum,
            ret_cfg;
            coords=coords,
            rhosolid=3000.0,
        )

        @test res.M_vent_H2O > 0.0
        @test res.M_vent_C > 0.0
        @test res.M_vent_N > 0.0
        @test res.M_vent_S > 0.0
        @test isapprox(
            res.M_vent_volatiles_total,
            res.M_vent_H2O + res.M_vent_C + res.M_vent_N + res.M_vent_S;
            atol=1.0e-12,
        )

        # Volatiles in markers decreased
        @test XH2Om[1] < 1.0
        @test XCm[1] < 500.0
        @test XNm[1] < 50.0
        @test XSm[1] < 1000.0

        # Now simulate prolonged extreme venting: 50 additional steps with large dt
        dt_extreme = 1.0e12 # s_vent * dt = 10.0 (extreme drainage)
        for step in 1:50
            drain_vented_marker_volatiles!(
                xm,
                ym,
                tm,
                tkm,
                XH2Om,
                XCm,
                XNm,
                XSm,
                S_vent_grid,
                dt_extreme,
                marknum,
                ret_cfg;
                coords=coords,
                rhosolid=3000.0,
            )
        end

        # Invariant: Even after extreme venting, markers NEVER drop below retention floor!
        # H2O floor = 50 ppm = 0.005 wt%
        @test isapprox(XH2Om[1], 0.005; atol=1.0e-6)
        @test XH2Om[1] >= 0.005
        # C floor = 50 ppm
        @test isapprox(XCm[1], 50.0; atol=1.0e-4)
        @test XCm[1] >= 50.0
        # N floor = 5 ppm
        @test isapprox(XNm[1], 5.0; atol=1.0e-4)
        @test XNm[1] >= 5.0
        # S floor = 100 ppm
        @test isapprox(XSm[1], 100.0; atol=1.0e-4)
        @test XSm[1] >= 100.0

        # Further venting should drain zero mass because all mobile volatiles are exhausted
        res_exhausted = drain_vented_marker_volatiles!(
            xm,
            ym,
            tm,
            tkm,
            XH2Om,
            XCm,
            XNm,
            XSm,
            S_vent_grid,
            dt_extreme,
            marknum,
            ret_cfg;
            coords=coords,
            rhosolid=3000.0,
        )
        @test isapprox(res_exhausted.M_vent_H2O, 0.0; atol=1.0e-12)
        @test isapprox(res_exhausted.M_vent_C, 0.0; atol=1.0e-12)
        @test isapprox(res_exhausted.M_vent_N, 0.0; atol=1.0e-12)
        @test isapprox(res_exhausted.M_vent_S, 0.0; atol=1.0e-12)
    end

    @testset "drain_vented_marker_volatiles! Sticky-Air Phase & Zero-Rate Invariance" begin
        coords = GridCoordinates(GridConfig(; Nx=5, Ny=5, xsize=100_000.0, ysize=100_000.0))
        marknum = 5
        xm = fill(50_000.0, marknum)
        ym = fill(50_000.0, marknum)
        tm_air = fill(3, marknum) # sticky air markers
        tkm = fill(1400.0, marknum)
        XH2Om = fill(1.0, marknum)
        XCm = fill(500.0, marknum)
        XNm = fill(50.0, marknum)
        XSm = fill(1000.0, marknum)

        ret_cfg = RetentionConfig(; active=true, venting_drainage_active=true)
        S_vent_grid = fill(1.0e-11, coords.Ny1, coords.Nx1)

        # Sticky air markers must not be drained
        res_air = drain_vented_marker_volatiles!(
            xm,
            ym,
            tm_air,
            tkm,
            XH2Om,
            XCm,
            XNm,
            XSm,
            S_vent_grid,
            1.0e10,
            marknum,
            ret_cfg;
            coords=coords,
        )
        @test isapprox(res_air.M_vent_volatiles_total, 0.0; atol=1.0e-14)
        @test isapprox(XH2Om[1], 1.0; atol=1.0e-14)

        # Zero venting rate grid
        tm_rock = fill(2, marknum)
        res_zero_s = drain_vented_marker_volatiles!(
            xm,
            ym,
            tm_rock,
            tkm,
            XH2Om,
            XCm,
            XNm,
            XSm,
            zeros(Float64, coords.Ny1, coords.Nx1),
            1.0e10,
            marknum,
            ret_cfg;
            coords=coords,
        )
        @test isapprox(res_zero_s.M_vent_volatiles_total, 0.0; atol=1.0e-14)
        @test isapprox(XH2Om[1], 1.0; atol=1.0e-14)
    end
end

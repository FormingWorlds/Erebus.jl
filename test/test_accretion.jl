using Test
using LinearAlgebra
using Erebus
using Erebus.Particles
using Erebus.Config

@testset "Accretion Engine Physics & Scaling Laws" begin
    # ---------------------------------------------------------------------
    # 1. Keplerian orbital kinematics
    # ---------------------------------------------------------------------
    @testset "Keplerian Kinematics" begin
        # Earth at 1 AU around 1 M_sun: Period ~ 1 year (3.15576e7 s)
        a_1au = 1.495978707e11
        M_sun = 1.98847e30
        omega_k = compute_keplerian_frequency(a_1au, M_sun)
        v_k = compute_keplerian_velocity(a_1au, M_sun)

        # Expected angular frequency ~ 2 pi / (1 yr in seconds)
        expected_omega = 2.0 * pi / 3.15576e7
        @test isapprox(omega_k, expected_omega, rtol=0.01)

        # Expected orbital velocity ~ 29.8 km/s
        @test isapprox(v_k, 29785.0, rtol=0.01)

        # Invariant: v_K == Omega_K * a
        @test isapprox(v_k, omega_k * a_1au, rtol=1e-12)

        # Domain error checks
        @test_throws DomainError compute_keplerian_frequency(-1.0, M_sun)
        @test_throws DomainError compute_keplerian_frequency(a_1au, -1.0)
        @test_throws DomainError compute_keplerian_velocity(0.0, M_sun)
        @test_throws DomainError compute_keplerian_velocity(a_1au, 0.0)
    end

    # ---------------------------------------------------------------------
    # 2. Sound speed and scale heights
    # ---------------------------------------------------------------------
    @testset "Sound Speed & Scale Heights" begin
        T_gas = 280.0
        c_s = compute_sound_speed(T_gas)
        # Expected sound speed for H2/He at 280 K ~ 1180 m/s
        @test c_s > 900.0
        @test c_s < 1400.0

        a_m = 1.495978707e11
        omega_k = compute_keplerian_frequency(a_m)
        H_gas = compute_gas_scale_height(c_s, omega_k)
        @test H_gas > 0.0
        # Aspect ratio H/r is typically ~ 0.03 - 0.07 in protoplanetary disks
        aspect_ratio = H_gas / a_m
        @test 0.02 < aspect_ratio < 0.10

        # Pebble scale height with settling
        St_large = 1.0
        alpha_turb = 1.0e-3
        H_peb_settled = compute_pebble_scale_height(H_gas, St_large, alpha_turb)
        @test H_peb_settled < H_gas
        @test isapprox(
            H_peb_settled / H_gas, sqrt(alpha_turb / (alpha_turb + St_large)), rtol=1e-10
        )

        # Limiting behavior: very small Stokes number => dust well-mixed with gas
        St_tiny = 1.0e-8
        H_peb_mixed = compute_pebble_scale_height(H_gas, St_tiny, alpha_turb)
        @test isapprox(H_peb_mixed, H_gas, rtol=1e-4)

        # Domain errors
        @test_throws DomainError compute_sound_speed(-10.0)
        @test_throws DomainError compute_gas_scale_height(-1.0, omega_k)
        @test_throws DomainError compute_pebble_scale_height(H_gas, -0.1, alpha_turb)
        @test_throws DomainError compute_pebble_scale_height(H_gas, 0.1, -1e-4)
    end

    # ---------------------------------------------------------------------
    # 3. Bondi, Hill, and Pebble Accretion Regimes
    # ---------------------------------------------------------------------
    @testset "Pebble Accretion Regimes" begin
        M_p = 1.0e21
        M_star = 1.98847e30
        a_m = 2.0 * 1.495978707e11
        T_gas = 150.0
        c_s = compute_sound_speed(T_gas)
        omega_k = compute_keplerian_frequency(a_m, M_star)
        Sigma_peb = compute_pebble_surface_density(2.0; Sigma_peb_0=50.0, p_peb=1.0)
        @test isapprox(Sigma_peb, 25.0, rtol=1e-10)

        R_B = compute_bondi_radius(M_p, c_s)
        R_H = compute_hill_radius(M_p, a_m, M_star)
        @test R_B > 0.0
        @test isapprox(R_H, 1.649714e8, rtol=1e-5)
        # Verify default stellar mass (M_SUN_KG)
        @test isapprox(compute_hill_radius(M_p, a_m), R_H, rtol=1e-12)
        @test_throws DomainError compute_hill_radius(-1.0, a_m, M_star)
        @test_throws DomainError compute_hill_radius(M_p, 0.0, M_star)
        @test_throws DomainError compute_hill_radius(M_p, a_m, 0.0)

        # Hill radius must scale as M^(1/3)
        R_H_8x = compute_hill_radius(8.0 * M_p, a_m, M_star)
        @test isapprox(R_H_8x, 2.0 * R_H, rtol=1e-10)

        # Accretion rate evaluations
        M_dot_bondi = compute_pebble_accretion_rate(
            M_p, M_star, a_m, Sigma_peb, 0.05, c_s, 1.0e-3; regime=:pebble_bondi
        )
        M_dot_hill = compute_pebble_accretion_rate(
            M_p, M_star, a_m, Sigma_peb, 0.05, c_s, 1.0e-3; regime=:pebble_hill
        )
        M_dot_auto = compute_pebble_accretion_rate(
            M_p, M_star, a_m, Sigma_peb, 0.05, c_s, 1.0e-3; regime=:pebble_auto
        )

        @test M_dot_bondi > 0.0
        @test M_dot_hill > 0.0
        @test M_dot_auto > 0.0
        # For M_p = 1e21 kg < M_trans, :pebble_auto selects Bondi regime
        @test isapprox(M_dot_auto, M_dot_bondi, rtol=1e-10)

        # For giant embryo mass M_p = 1e25 kg > M_trans, :pebble_auto selects Hill regime
        M_dot_auto_large = compute_pebble_accretion_rate(
            1.0e25, M_star, a_m, Sigma_peb, 0.05, c_s, 1.0e-3; regime=:pebble_auto
        )
        M_dot_hill_large = compute_pebble_accretion_rate(
            1.0e25, M_star, a_m, Sigma_peb, 0.05, c_s, 1.0e-3; regime=:pebble_hill
        )
        @test isapprox(M_dot_auto_large, M_dot_hill_large, rtol=1e-10)

        # Invalid regime throws ArgumentError
        @test_throws ArgumentError compute_pebble_accretion_rate(
            M_p, M_star, a_m, Sigma_peb, 0.05, c_s, 1.0e-3; regime=:invalid_mode
        )

        # Scaling: accretion rate strictly increases with target mass
        M_dot_larger = compute_pebble_accretion_rate(
            10.0 * M_p, M_star, a_m, Sigma_peb, 0.05, c_s, 1.0e-3; regime=:pebble_hill
        )
        @test M_dot_larger > M_dot_hill
    end

    # ---------------------------------------------------------------------
    # 4. Safronov Gravitational Focusing
    # ---------------------------------------------------------------------
    @testset "Safronov Gravitational Focusing" begin
        M = 5.0e19
        R = 35000.0
        Sigma_pl = 100.0
        v_disp = 100.0
        omega_k = 1.0e-7

        M_dot_saf = compute_safronov_accretion_rate(M, R, Sigma_pl, v_disp, omega_k)
        @test M_dot_saf > 0.0

        # Geometric cross-section limit: as v_disp -> very large, Theta -> 0, F_g -> 1
        v_disp_huge = 1.0e7
        M_dot_geom = compute_safronov_accretion_rate(M, R, Sigma_pl, v_disp_huge, omega_k)
        expected_geom = pi * (R^2) * Sigma_pl * omega_k
        @test isapprox(M_dot_geom, expected_geom, rtol=1e-4)

        # Gravitational focusing always increases collision cross-section
        @test M_dot_saf > M_dot_geom
    end

    # ---------------------------------------------------------------------
    # 5. Analytical Growth Modes & Timing Activation
    # ---------------------------------------------------------------------
    @testset "Growth Modes & Saturation" begin
        acc_cfg = AccretionConfig(
            active=true,
            mode=:constant_rate,
            dM_dt_constant=5.0e11,
            M_initial=1.0e17,
            M_target=1.0e18,
            R_initial=20000.0,
            R_target=40000.0,
            t_start_myr=0.1,
            t_duration_myr=1.0,
        )
        sec_yr = 3.15576e7
        t_before = 0.05 * 1.0e6 * sec_yr
        t_active = 0.5 * 1.0e6 * sec_yr
        t_after = 2.0 * 1.0e6 * sec_yr

        # Gated outside time window
        @test isapprox(
            compute_accretion_rate(t_before, 1.0e17, 20000.0, acc_cfg), 0.0, atol=1e-12
        )
        @test isapprox(
            compute_accretion_rate(t_after, 1.0e17, 20000.0, acc_cfg), 0.0, atol=1e-12
        )

        # Active within window
        rate_active = compute_accretion_rate(t_active, 1.0e17, 20000.0, acc_cfg)
        @test isapprox(rate_active, 5.0e11, rtol=1e-10)

        # Target saturation cutoff
        rate_saturated_m = compute_accretion_rate(t_active, 1.0e19, 20000.0, acc_cfg)
        rate_saturated_r = compute_accretion_rate(t_active, 1.0e17, 45000.0, acc_cfg)
        @test isapprox(rate_saturated_m, 0.0, atol=1e-12)
        @test isapprox(rate_saturated_r, 0.0, atol=1e-12)

        # Inactive config returns 0
        cfg_off = AccretionConfig(active=false)
        @test isapprox(
            compute_accretion_rate(t_active, 1.0e17, 20000.0, cfg_off), 0.0, atol=1e-12
        )

        # Linear radius growth mode
        acc_cfg_lin = AccretionConfig(
            active=true,
            mode=:linear_radius,
            dR_dt_constant=1.0e-3,
            rho_bulk=3000.0,
            t_start_myr=0.0,
            t_duration_myr=2.0,
        )
        R_curr = 25000.0
        rate_lin = compute_accretion_rate(t_active, 1.0e17, R_curr, acc_cfg_lin)
        expected_rate_lin = 4.0 * pi * (R_curr^2) * 3000.0 * 1.0e-3
        @test isapprox(rate_lin, expected_rate_lin, rtol=1e-10)
    end

    # ---------------------------------------------------------------------
    # 6. Exact 3D-to-2D Spherical Geometric Mapping
    # ---------------------------------------------------------------------
    @testset "Exact Spherical Shell Mapping" begin
        R0 = 20000.0
        rho = 3200.0
        delta_M = 5.0e17

        delta_R = compute_radius_increment(R0, delta_M, rho)
        @test delta_R > 0.0

        # Exact spherical volume conservation check:
        # V_new - V_old = 4/3 pi ( (R0 + dR)^3 - R0^3 ) == delta_M / rho
        V_old = (4.0 / 3.0) * pi * (R0^3)
        V_new = (4.0 / 3.0) * pi * ((R0 + delta_R)^3)
        delta_V = V_new - V_old
        expected_delta_V = delta_M / rho
        @test isapprox(delta_V, expected_delta_V, rtol=1e-11)

        # Zero mass increment produces zero radius increment
        @test isapprox(compute_radius_increment(R0, 0.0, rho), 0.0, atol=1e-14)

        # Numerical stability against catastrophic cancellation for small dM down to 1e-2 kg
        for dM_small in [1.0e6, 1.0e3, 1.0, 1.0e-2]
            dR_small = compute_radius_increment(R0, dM_small, rho)
            @test dR_small > 0.0
            dR_taylor = dM_small / (4.0 * pi * rho * (R0^2))
            @test isapprox(dR_small, dR_taylor, rtol=1e-4)
        end

        # Domain error on unphysical input
        @test_throws DomainError compute_radius_increment(-100.0, delta_M, rho)
        @test_throws DomainError compute_radius_increment(R0, -1.0, rho)
        @test_throws DomainError compute_radius_increment(R0, delta_M, 0.0)
    end

    # ---------------------------------------------------------------------
    # 7. Impact Heating Thermodynamics
    # ---------------------------------------------------------------------
    @testset "Impact Heating Thermodynamics" begin
        M = 2.0e19
        R = 30000.0
        cp = 1000.0
        h_imp = 0.6
        v_inf = 500.0

        u_acc, delta_T = compute_impact_heating(M, R; h_impact=h_imp, c_p=cp, v_inf=v_inf)
        @test u_acc > 0.0
        @test delta_T > 0.0

        # Invariant: u_acc = G M / R + 0.5 v_inf^2
        expected_u = (6.67430e-11 * M / R) + 0.5 * (v_inf^2)
        @test isapprox(u_acc, expected_u, rtol=1e-12)

        # Temperature rise is proportional to retention efficiency h_imp
        @test isapprox(delta_T, h_imp * expected_u / cp, rtol=1e-12)

        # Zero retention produces zero temperature rise
        _, delta_T_zero = compute_impact_heating(M, R; h_impact=0.0, c_p=cp, v_inf=v_inf)
        @test isapprox(delta_T_zero, 0.0, atol=1e-14)

        # Domain errors
        @test_throws DomainError compute_impact_heating(-1.0, R)
        @test_throws DomainError compute_impact_heating(M, 0.0)
        @test_throws DomainError compute_impact_heating(M, R; h_impact=-0.1)
        @test_throws DomainError compute_impact_heating(M, R; h_impact=1.2)
        @test_throws DomainError compute_impact_heating(M, R; c_p=0.0)
    end

    # ---------------------------------------------------------------------
    # 8. Snowline Volatile Coupling
    # ---------------------------------------------------------------------
    @testset "Snowline Volatile Coupling" begin
        # Outside water snowline (cold)
        T_cold = 120.0
        XW_cold, H2O_cold = evaluate_snowline_water_content(
            T_cold;
            T_snowline_cond=160.0,
            XW_wet=0.45,
            XW_dry=0.02,
            H2O_wet_wtpct=12.0,
            H2O_dry_wtpct=0.05,
        )
        @test isapprox(XW_cold, 0.45, rtol=1e-10)
        @test isapprox(H2O_cold, 12.0, rtol=1e-10)

        # Inside water snowline (warm)
        T_warm = 220.0
        XW_warm, H2O_warm = evaluate_snowline_water_content(
            T_warm;
            T_snowline_cond=160.0,
            XW_wet=0.45,
            XW_dry=0.02,
            H2O_wet_wtpct=12.0,
            H2O_dry_wtpct=0.05,
        )
        @test isapprox(XW_warm, 0.02, rtol=1e-10)
        @test isapprox(H2O_warm, 0.05, rtol=1e-10)
    end

    # ---------------------------------------------------------------------
    # 9. Dynamic Marker Boundary Expansion
    # ---------------------------------------------------------------------
    @testset "Dynamic Boundary Expansion" begin
        # Create small test grid of markers around center (50000, 50000)
        xc = 50000.0
        yc = 50000.0
        R_init = 10000.0
        delta_R = 5000.0
        R_new = R_init + delta_R

        radii = [5000.0, 8000.0, 12000.0, 14000.0, 18000.0, 22000.0]
        n_m = length(radii)
        xm = [xc + r for r in radii]
        ym = fill(yc, n_m)

        # Markers 1, 2 inside R_init (crust/core)
        # Markers 3, 4, 5, 6 outside R_init (sticky air, tm = 3)
        tm = [2, 2, 3, 3, 3, 3]
        tkm = fill(150.0, n_m)
        phim = fill(0.10, n_m)
        XWsolidm0 = fill(0.0, n_m)
        Xfe_bulk = fill(0.0, n_m)
        Xfem = fill(0.0, n_m)
        t_acc = fill(0.0, n_m)

        # Advance boundary by delta_R
        n_conv = advance_accretion_boundary!(
            R_init,
            delta_R,
            xm,
            ym,
            tm,
            tkm,
            phim,
            XWsolidm0,
            Xfe_bulk,
            Xfem;
            xcenter=xc,
            ycenter=yc,
            T_accreted=250.0,
            phi_accreted=0.30,
            XWsolid_accreted=0.40,
            Xfe_accreted=0.15,
            t_accreted=t_acc,
            current_time=1.0e12,
        )

        # Markers 3 and 4 (radii 12000 and 14000) are inside R_new (15000)
        # Markers 5 and 6 (radii 18000 and 22000) are outside R_new
        @test n_conv == 2

        # Converted markers must have updated properties
        @test tm[3] == 2
        @test tm[4] == 2
        @test isapprox(tkm[3], 250.0, rtol=1e-10)
        @test isapprox(phim[3], 0.30, rtol=1e-10)
        @test isapprox(XWsolidm0[3], 0.40, rtol=1e-10)
        @test isapprox(Xfe_bulk[3], 0.15, rtol=1e-10)
        @test isapprox(t_acc[3], 1.0e12, rtol=1e-10)

        # Unconverted air markers outside R_new remain tm = 3
        @test tm[5] == 3
        @test tm[6] == 3

        # Primordial markers inside R_init must not be overwritten
        @test tm[1] == 2
        @test isapprox(tkm[1], 150.0, rtol=1e-10)
        @test isapprox(phim[1], 0.10, rtol=1e-10)
        @test isapprox(t_acc[1], 0.0, atol=1e-14)
    end

    # ---------------------------------------------------------------------
    # 10. Radiogenic Onion-Shell Clock Inheritance
    # ---------------------------------------------------------------------
    @testset "Radiogenic Onion-Shell Inheritance" begin
        # Markers accreted at early epoch (0.5 Myr) vs late epoch (2.5 Myr)
        t1 = 0.5 * 1.0e6 * SEC_PER_YEAR
        t2 = 2.5 * 1.0e6 * SEC_PER_YEAR

        hr_early, _, _ = Erebus.calculate_radioactive_heating(true, false, t1)
        hr_late, _, _ = Erebus.calculate_radioactive_heating(true, false, t2)

        # Early interior experiences much higher radiogenic heating rate than late accreted shell
        @test hr_early[2] > hr_late[2]
        tau_al = Erebus.t_half_al / log(2.0)
        decay_factor = exp(-(t2 - t1) / tau_al)
        @test isapprox(hr_late[2] / hr_early[2], decay_factor, rtol=1e-5)
        @test decay_factor < 0.16
    end

    # ---------------------------------------------------------------------
    # 11. Marker Replenishment Accretion Array Synchronization
    # ---------------------------------------------------------------------
    @testset "Marker Replenishment Synchronization" begin
        # Setup population tracking with 25 markers in a small grid
        coords = GridCoordinates(GridConfig(; Nx=5, Ny=5, xsize=10000.0, ysize=10000.0))
        mdis, mnum = setup_marker_geometry_helpers(coords)
        marknum = 25
        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0, Fm) = setup_marker_properties(
            marknum, coords
        )
        (rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm, tkm_rhocptotalm, etafluidcur_inv_kphim, inv_gggtotalm, fricttotalm, cohestotalm, tenstotalm, rhofluidcur, alphasolidcur, alphafluidcur) = setup_marker_properties_helpers(
            marknum
        )

        for i in 1:marknum
            xm[i] = 1000.0 + (i - 1) * 300.0
            ym[i] = 1000.0 + (i - 1) * 300.0
            tm[i] = 1
            tkm[i] = 300.0
        end

        t_acc = fill(5.0e12, marknum)

        new_marknum = replenish_markers!(
            xm,
            ym,
            tm,
            tkm,
            phim,
            sxxm,
            sxym,
            etavpm,
            phinewm,
            pfm0,
            XWsolidm,
            XWsolidm0,
            rhototalm,
            rhocptotalm,
            etatotalm,
            hrtotalm,
            ktotalm,
            inv_gggtotalm,
            fricttotalm,
            cohestotalm,
            tenstotalm,
            rhofluidcur,
            alphasolidcur,
            alphafluidcur,
            tkm_rhocptotalm,
            etafluidcur_inv_kphim,
            mdis,
            mnum;
            Fm=Fm,
            randomized=false,
            coords=coords,
            t_accreted=t_acc,
        )

        @test new_marknum > marknum
        @test length(t_acc) == new_marknum
        @test length(t_acc) == length(xm)
    end

    @testset "Accreted Marker Constitutive Properties and Metal Melt Initialization" begin
        cfg = SimulationConfig(;
            grid=GridConfig(Nx=5, Ny=5, xsize=40_000.0, ysize=40_000.0),
            geometry=GeometryConfig(
                rplanet=8_000.0, rcrust=8_000.0, xcenter=20_000.0, ycenter=20_000.0
            ),
            accretion=AccretionConfig(
                active=true,
                mode=:constant_rate,
                R_initial=8_000.0,
                R_target=18_000.0,
                M_initial=1.0e18,
                M_target=1.0e20,
                dM_dt_constant=1.0e14,
                phi_accreted=0.25,
                Xfe_bulk_accreted=0.12,
            ),
            coreformation=CoreFormationConfig(
                percolation_active=true, Xfe_bulk=0.12, T_eutectic=1213.0, dT_metal=50.0
            ),
            redox=RedoxConfig(active=true),
        )
        coords = GridCoordinates(cfg.grid)
        markers = init_marker_arrays(2, cfg, coords)
        markers.core.xm[1] = 20_000.0 + 10_000.0
        markers.core.ym[1] = 20_000.0
        markers.core.tm[1] = 3
        markers.core.etavpm[1] = cfg.materials.etasolidm[3]
        markers.core.fricttotalm[1] = 0.0
        markers.core.alphasolidcur[1] = 0.0
        markers.core.alphafluidcur[1] = 0.0
        markers.core.sxxm[1] = 5.0e5
        markers.core.sxym[1] = -2.0e5
        markers.core.Fm[1] = 0.4

        markers.core.xm[2] = 20_000.0 + 16_000.0
        markers.core.ym[2] = 20_000.0
        markers.core.tm[2] = 3

        T_hot = 1238.0
        n_conv = advance_accretion_boundary!(
            8_000.0,
            4_000.0,
            markers;
            cfg=cfg,
            coords=coords,
            xcenter=20_000.0,
            ycenter=20_000.0,
            T_accreted=T_hot,
            phi_accreted=0.25,
            Xfe_accreted=0.12,
        )
        @test n_conv == 1
        @test markers.core.tm[1] == 2
        @test markers.core.tm[2] == 3
        @test isapprox(markers.core.etavpm[1], cfg.materials.etasolidm[2]; rtol=1e-12)
        @test isapprox(
            markers.core.inv_gggtotalm[1], inv(cfg.materials.gggsolidm[2]); rtol=1e-12
        )
        @test isapprox(
            markers.core.fricttotalm[1], cfg.materials.frictsolidm[2]; rtol=1e-12
        )
        @test isapprox(
            markers.core.cohestotalm[1], cfg.materials.cohessolidm[2]; rtol=1e-12
        )
        @test isapprox(markers.core.tenstotalm[1], cfg.materials.tenssolidm[2]; rtol=1e-12)
        @test isapprox(
            markers.core.alphasolidcur[1], cfg.materials.alphasolidm[2]; rtol=1e-12
        )
        @test isapprox(
            markers.core.alphafluidcur[1], cfg.materials.alphafluidm[2]; rtol=1e-12
        )
        @test isapprox(markers.core.rhofluidcur[1], cfg.materials.rhofluidm[2]; rtol=1e-12)
        expected_rho = cfg.materials.rhosolidm[2] * 0.75 + cfg.materials.rhofluidm[2] * 0.25
        expected_rhocp =
            cfg.materials.rhocpsolidm[2] * 0.75 + cfg.materials.rhocpfluidm[2] * 0.25
        expected_k = cfg.materials.ksolidm[2] * 0.75 + cfg.materials.kfluidm[2] * 0.25
        @test isapprox(markers.core.rhototalm[1], expected_rho; rtol=1e-12)
        @test isapprox(markers.core.rhocptotalm[1], expected_rhocp; rtol=1e-12)
        @test isapprox(markers.core.ktotalm[1], expected_k; rtol=1e-12)
        @test isapprox(
            markers.core.etatotalm[1],
            cfg.materials.etasolidm[2] * exp(-28.0 * 0.25);
            rtol=1e-12,
        )
        @test isapprox(markers.core.tkm_rhocptotalm[1], T_hot * expected_rhocp; rtol=1e-12)
        @test isapprox(
            markers.core.pfm0[1],
            max(cfg.geometry.psurface, cfg.disk.p_amb_disk);
            rtol=1e-12,
        )
        @test isapprox(markers.core.Fm[1], 0.0; atol=1e-15)
        @test isapprox(markers.core.sxxm[1], 0.0; atol=1e-15)
        @test isapprox(markers.core.sxym[1], 0.0; atol=1e-15)
        @test isapprox(markers.core.hrtotalm[1], 0.0; atol=1e-15)
        @test isapprox(
            markers.core.etafluidcur_inv_kphim[1],
            Erebus.ηᶠcur_inv_kᵠ(
                cfg.materials.kphim0[2],
                0.25,
                cfg.materials.etafluidm[2];
                phim0=cfg.thermodynamics.phim0,
            );
            rtol=1e-12,
        )

        expected_F_fe = compute_metal_melt_fraction(T_hot; T_eutectic=1213.0, dT_metal=50.0)
        @test isapprox(expected_F_fe, 0.5; rtol=1e-12)
        @test isapprox(markers.groups.metal.Xfem[1], 0.12 * expected_F_fe; rtol=1e-12)
        @test isapprox(markers.groups.metal.Xfem0[1], 0.12 * expected_F_fe; rtol=1e-12)

        # Fallback path when cfg === nothing
        n_conv_nocfg = advance_accretion_boundary!(
            12_000.0,
            6_000.0,
            markers;
            cfg=nothing,
            coords=coords,
            xcenter=20_000.0,
            ycenter=20_000.0,
            T_accreted=300.0,
            phi_accreted=0.20,
            Xfe_accreted=0.10,
        )
        @test n_conv_nocfg == 1
        @test markers.core.tm[2] == 2
        @test isapprox(markers.core.fricttotalm[2], Erebus.frictsolidm[2]; rtol=1e-12)
        @test isapprox(markers.core.alphasolidcur[2], Erebus.alphasolidm[2]; rtol=1e-12)
        @test isapprox(markers.groups.metal.Xfem[2], 0.0; atol=1e-15)
    end

    @testset "Chondritic Metal Mass-to-Volume Fraction Conversion" begin
        phi_2arg = chondritic_phi_fe(3500.0, 5450.0)
        phi_kw = chondritic_phi_fe(3500.0, 5450.0; w_fe=X_FE_REF_CHONDRITE)
        @test isapprox(phi_2arg, phi_kw; rtol=1e-15)
        @test isapprox(
            metal_volume_to_mass_fraction(phi_kw, 5450.0, 3500.0),
            X_FE_REF_CHONDRITE;
            rtol=1e-14,
        )
        @test isapprox(chondritic_phi_fe(3500.0, 5450.0; w_fe=0.0), 0.0; atol=1e-15)
        @test isapprox(chondritic_phi_fe(3500.0, 5450.0; w_fe=1.0), 1.0; rtol=1e-15)
        @test isapprox(chondritic_phi_fe(0.0, 5450.0; w_fe=0.2), 0.0; atol=1e-15)
        @test isapprox(
            chondritic_phi_fe(3500.0, 0.0; w_fe=0.2),
            chondritic_phi_fe(3500.0, 5450.0; w_fe=0.2);
            rtol=1e-15,
        )
    end

    @testset "Accretion Timestep Synchronization and Reconciliation" begin
        cfg = SimulationConfig(;
            grid=GridConfig(Nx=7, Ny=7, xsize=40_000.0, ysize=40_000.0),
            geometry=GeometryConfig(
                rplanet=8_000.0, rcrust=8_000.0, xcenter=20_000.0, ycenter=20_000.0
            ),
            time=TimeConfig(start_time=0.0),
            accretion=AccretionConfig(
                active=true,
                mode=:constant_rate,
                t_start_myr=0.0,
                R_initial=8_000.0,
                R_target=18_000.0,
                M_initial=1.0e18,
                M_target=1.0e20,
                dM_dt_constant=1.0e12,
                rho_bulk=3000.0,
            ),
            coreformation=CoreFormationConfig(percolation_active=true, Xfe_bulk=0.10),
        )
        mktempdir() do tmp_dir
            state, coords, ws, telemetry_io, _ = init_simulation(cfg; output_path=tmp_dir)
            telemetry_io !== nothing && close(telemetry_io)

            state.dt = 1000.0
            snapshot = snapshot_step_state(state, coords, ws)
            accrete!(state, coords, cfg)
            dM_trial = 1.0e12 * 1000.0
            dR_trial = compute_radius_increment(8_000.0, dM_trial, 3000.0)
            @test isapprox(state.accumulators.M_planet_val, 1.0e18 + dM_trial; rtol=1e-12)
            @test isapprox(state.accumulators.rplanet, 8_000.0 + dR_trial; rtol=1e-12)

            # Place a test marker between r(dt=250) and r(dt=1000) that was converted on trial
            dM_final = 1.0e12 * 250.0
            dR_final = compute_radius_increment(8_000.0, dM_final, 3000.0)
            r_mid = 8_000.0 + 0.5 * (dR_final + dR_trial)
            state.markers.core.xm[1] = 20_000.0 + r_mid
            state.markers.core.ym[1] = 20_000.0
            snapshot.markers.core.xm[1] = 20_000.0 + r_mid
            snapshot.markers.core.ym[1] = 20_000.0
            snapshot.markers.core.tm[1] = 3
            state.markers.core.tm[1] = 2

            # Shrink state.dt to 250.0 s and reconcile
            state.dt = 250.0
            reconcile_accretion_timestep!(state, snapshot, cfg)
            @test isapprox(state.accumulators.M_planet_val, 1.0e18 + dM_final; rtol=1e-12)
            @test isapprox(state.accumulators.M_accreted_total, dM_final; rtol=1e-12)
            @test isapprox(state.accumulators.rplanet, 8_000.0 + dR_final; rtol=1e-12)
            @test state.markers.core.tm[1] == 3

            # Verify zero-growth branch when dM_acc == 0
            state.dt = 0.0
            reconcile_accretion_timestep!(state, snapshot, cfg)
            @test isapprox(state.accumulators.M_planet_val, 1.0e18; rtol=1e-12)
            @test isapprox(state.accumulators.M_accreted_total, 0.0; atol=1e-15)
            @test isapprox(state.accumulators.rplanet, 8_000.0; rtol=1e-12)

            # Verify no-op early return when accretion is inactive
            state.accumulators.rplanet = 9_500.0
            reconcile_accretion_timestep!(state, snapshot, SimulationConfig())
            @test isapprox(state.accumulators.rplanet, 9_500.0; rtol=1e-12)
        end

        # Full simulation_loop integration where adaptive timestep shrinks state.dt
        mktempdir() do tmp_dir2
            cfg_loop = SimulationConfig(;
                grid=GridConfig(Nx=9, Ny=9, xsize=40_000.0, ysize=40_000.0),
                geometry=GeometryConfig(
                    rplanet=10_000.0, rcrust=10_000.0, xcenter=20_000.0, ycenter=20_000.0
                ),
                time=TimeConfig(
                    dt_initial=1.0e4,
                    dt_longest=1.0e4,
                    dxymax=1.0e-6,
                    start_time=0.0,
                    endtime=1.0e7,
                    start_step=1,
                    n_steps=1,
                ),
                solver=SolverConfig(max_plastic_iterations=2),
                output=OutputConfig(output_dir=tmp_dir2, savematstep=0, visstep=0),
                accretion=AccretionConfig(
                    active=true,
                    mode=:constant_rate,
                    t_start_myr=0.0,
                    R_initial=10_000.0,
                    R_target=18_000.0,
                    M_initial=1.0e18,
                    M_target=1.0e22,
                    dM_dt_constant=1.0e7,
                    rho_bulk=3000.0,
                ),
            )
            state_loop = simulation_loop(cfg_loop; output_path=tmp_dir2)
            dt_trial_sec = cfg_loop.time.dt_initial * cfg_loop.time.yearlength
            @test state_loop.timesum < dt_trial_sec
            expected_dM = cfg_loop.accretion.dM_dt_constant * state_loop.timesum
            expected_dR = compute_radius_increment(10_000.0, expected_dM, 3000.0)
            @test isapprox(
                state_loop.accumulators.M_accreted_total, expected_dM; rtol=1e-12
            )
            @test isapprox(
                state_loop.accumulators.M_planet_val, 1.0e18 + expected_dM; rtol=1e-12
            )
            @test isapprox(
                state_loop.accumulators.rplanet, 10_000.0 + expected_dR; rtol=1e-12
            )
        end
    end
end

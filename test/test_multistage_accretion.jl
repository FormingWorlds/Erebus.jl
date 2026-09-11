using Test
using Erebus

@testset "Multi-Stage Accretion Sequence Physics & Transitions" begin
    # Physical constants for testing
    M_sun = 1.98847e30
    AU = 1.495978707e11
    a_test = 2.5 * AU
    T_disk = 150.0
    c_s = compute_sound_speed(T_disk)
    v_K = compute_keplerian_velocity(a_test, M_sun)
    Omega_K = compute_keplerian_frequency(a_test, M_sun)
    St_test = 0.05

    # ---------------------------------------------------------------------
    # 1. Sub-Keplerian Gas Headwind Velocity
    # ---------------------------------------------------------------------
    @testset "Gas Headwind Velocity Physics" begin
        v_hw = compute_headwind_velocity(c_s, v_K)
        eta_expected = 1.5 * (c_s / v_K)^2
        v_hw_expected = eta_expected * v_K

        # Pinned against analytical reference value (42.133 m/s at 2.5 AU, 150 K)
        @test isapprox(v_hw, 42.13301, rtol=1e-4)
        @test isapprox(v_hw, v_hw_expected, rtol=1e-12)
        @test v_hw > 0.0
        @test v_hw < c_s

        # Error guards on non-positive or non-finite inputs
        @test_throws DomainError compute_headwind_velocity(-c_s, v_K)
        @test_throws DomainError compute_headwind_velocity(c_s, 0.0)
        @test_throws DomainError compute_headwind_velocity(NaN, v_K)
        @test_throws DomainError compute_headwind_velocity(c_s, Inf)
    end

    # ---------------------------------------------------------------------
    # 2. Pebble Accretion Onset Mass (Visser & Ormel 2016)
    # ---------------------------------------------------------------------
    @testset "Pebble Onset Mass (Visser & Ormel 2016)" begin
        M_onset = compute_pebble_onset_mass(M_sun, a_test, St_test, c_s; f_onset=1.0)
        v_hw = compute_headwind_velocity(c_s, v_K)
        M_onset_expected = (v_hw^3) * St_test / (Erebus.G_GRAV * Omega_K)

        # Pinned against independent reference calculation (1.1124e21 kg at 2.5 AU, 150 K, St=0.05)
        @test isapprox(M_onset, 1.1124e21, rtol=1e-3)
        @test isapprox(M_onset, M_onset_expected, rtol=1e-12)
        @test M_onset > 0.0

        # Physical order of magnitude at 2.5 AU, T=150 K, St=0.05:
        # M_onset ~ 1e21 kg (R ~ 445 km for rho ~ 3000 kg/m^3)
        R_onset = cbrt(3.0 * M_onset / (4.0 * pi * 3000.0))
        @test 20000.0 <= R_onset <= 600000.0
        @test 1.0e19 <= M_onset <= 1.0e23

        # Analytical scaling: linear in Stokes number
        M_onset_2St = compute_pebble_onset_mass(
            M_sun, a_test, 2.0 * St_test, c_s; f_onset=1.0
        )
        @test isapprox(M_onset_2St, 2.0 * M_onset, rtol=1e-12)

        # Analytical scaling: linear in calibration factor f_onset
        M_onset_2f = compute_pebble_onset_mass(M_sun, a_test, St_test, c_s; f_onset=2.5)
        @test isapprox(M_onset_2f, 2.5 * M_onset, rtol=1e-12)

        # Zero Stokes limit: M_onset -> 0
        M_onset_zero = compute_pebble_onset_mass(M_sun, a_test, 0.0, c_s; f_onset=1.0)
        @test isapprox(M_onset_zero, 0.0, atol=1e-30)
        @test M_onset_zero >= 0.0

        # Domain error guards
        @test_throws DomainError compute_pebble_onset_mass(0.0, a_test, St_test, c_s)
        @test_throws DomainError compute_pebble_onset_mass(M_sun, -a_test, St_test, c_s)
        @test_throws DomainError compute_pebble_onset_mass(M_sun, a_test, -0.01, c_s)
        @test_throws DomainError compute_pebble_onset_mass(M_sun, a_test, St_test, -c_s)
        @test_throws DomainError compute_pebble_onset_mass(
            M_sun, a_test, St_test, c_s; f_onset=-1.0
        )
        @test_throws DomainError compute_pebble_onset_mass(
            M_sun, a_test, St_test, c_s; f_onset=0.0
        )
    end

    # ---------------------------------------------------------------------
    # 3. Pebble Isolation Mass (Lambrechts et al. 2014)
    # ---------------------------------------------------------------------
    @testset "Pebble Isolation Mass Physics" begin
        M_iso = compute_pebble_isolation_mass(M_sun, a_test, c_s; f_iso=0.5)
        h_aspect = c_s / v_K
        M_iso_expected = 0.5 * M_sun * (h_aspect^3)

        # Pinned against independent reference value (5.7246e25 kg, approx 9.58 Earth masses)
        @test isapprox(M_iso, 5.7246e25, rtol=1e-3)
        @test isapprox(M_iso, M_iso_expected, rtol=1e-12)
        @test M_iso > 0.0

        # Physical order of magnitude: M_iso ~ 10 Earth masses ~ 1e25 to 1e26 kg
        @test M_iso > 1.0e24
        @test M_iso < 1.0e27

        # Scaling: linear in f_iso
        M_iso_double = compute_pebble_isolation_mass(M_sun, a_test, c_s; f_iso=1.0)
        @test isapprox(M_iso_double, 2.0 * M_iso, rtol=1e-12)

        # Hierarchy invariant: M_onset << M_iso
        M_onset = compute_pebble_onset_mass(M_sun, a_test, St_test, c_s)
        @test M_onset < 1.0e-3 * M_iso

        # Domain error guards
        @test_throws DomainError compute_pebble_isolation_mass(-M_sun, a_test, c_s)
        @test_throws DomainError compute_pebble_isolation_mass(M_sun, -a_test, c_s)
        @test_throws DomainError compute_pebble_isolation_mass(M_sun, a_test, 0.0)
        @test_throws DomainError compute_pebble_isolation_mass(
            M_sun, a_test, c_s; f_iso=-0.5
        )
    end

    # ---------------------------------------------------------------------
    # 4. Configuration Schema and Validation
    # ---------------------------------------------------------------------
    @testset "Multistage Configuration Validation" begin
        # Default construction
        cfg = default_config()
        @test cfg.accretion.mode in Set([
            :constant_rate,
            :linear_radius,
            :exponential,
            :safronov,
            :pebble_bondi,
            :pebble_hill,
            :pebble_auto,
            :multistage,
        ])
        @test cfg.accretion.f_onset > 0.0
        @test cfg.accretion.f_iso > 0.0
        @test cfg.accretion.transition_width > 0.0

        # Valid multistage configuration
        cfg_multi = SimulationConfig(
            accretion=AccretionConfig(
                active=true,
                mode=:multistage,
                stage1_mode=:safronov,
                stage2_mode=:pebble_auto,
                stage3_mode=:safronov,
                f_onset=1.0,
                f_iso=0.5,
                transition_smoothing=true,
                transition_width=0.10,
                M_target=1.0e28,
                R_target=5000000.0,
            ),
            grid=GridConfig(Nx=11, Ny=11, xsize=12000000.0, ysize=12000000.0),
            geometry=GeometryConfig(xcenter=6000000.0, ycenter=6000000.0),
        )
        @test validate_config(cfg_multi) === nothing

        # Valid disabling of M_iso (<= 0.0)
        cfg_dis1 = SimulationConfig(
            accretion=AccretionConfig(active=true, mode=:multistage, M_iso=0.0)
        )
        cfg_dis2 = SimulationConfig(
            accretion=AccretionConfig(active=true, mode=:multistage, M_iso=-1.0)
        )
        @test validate_config(cfg_dis1) === nothing
        @test validate_config(cfg_dis2) === nothing

        # Invalid stage1_mode
        cfg_bad_stage = SimulationConfig(
            accretion=AccretionConfig(
                active=true, mode=:multistage, stage1_mode=:invalid_mode
            ),
        )
        @test_throws ArgumentError validate_config(cfg_bad_stage)

        # Invalid f_onset
        cfg_bad_fonset = SimulationConfig(
            accretion=AccretionConfig(active=true, mode=:multistage, f_onset=-1.0)
        )
        @test_throws ArgumentError validate_config(cfg_bad_fonset)

        # Invalid transition_width (> 0.5)
        cfg_bad_width = SimulationConfig(
            accretion=AccretionConfig(active=true, mode=:multistage, transition_width=0.8)
        )
        @test_throws ArgumentError validate_config(cfg_bad_width)

        # Invalid M_onset >= M_iso when both are specified and M_iso > 0
        cfg_bad_hierarchy = SimulationConfig(
            accretion=AccretionConfig(
                active=true, mode=:multistage, M_onset=1.0e25, M_iso=1.0e24
            ),
        )
        @test_throws ArgumentError validate_config(cfg_bad_hierarchy)
    end

    # ---------------------------------------------------------------------
    # 5. Multi-Stage Regime Switching and Rate Calculations
    # ---------------------------------------------------------------------
    @testset "Sequential Accretion Regimes Dispatch" begin
        disk_cfg = DiskConfig(
            orbital_distance_au=2.5, stellar_mass_msun=1.0, t_ambient=150.0
        )
        sec_yr = 3.15576e7
        t_active = 0.5 * 1.0e6 * sec_yr

        # Sharp transition setup (smoothing = false) for unambiguous regime verification
        acc_sharp = AccretionConfig(
            active=true,
            mode=:multistage,
            stage1_mode=:safronov,
            stage2_mode=:pebble_auto,
            stage3_mode=:safronov,
            transition_smoothing=false,
            t_start_myr=0.0,
            t_duration_myr=2.0,
            M_initial=1.0e17,
            M_target=1.0e28,
            R_initial=20000.0,
            R_target=5000000.0,
            Sigma_pl_0=100.0,
            Sigma_peb_0=50.0,
            v_disp_kms=0.1,
            stokes_number=0.05,
        )

        M_onset = compute_pebble_onset_mass(M_sun, a_test, 0.05, c_s; f_onset=1.0)
        M_iso = compute_pebble_isolation_mass(M_sun, a_test, c_s; f_iso=0.5)

        # Case A: Low mass seed (Stage 1, Safronov)
        M_sub = 0.1 * M_onset
        R_sub = 25000.0
        rate_stage1 = compute_accretion_rate(t_active, M_sub, R_sub, acc_sharp, disk_cfg)
        expected_safronov_sub = compute_safronov_accretion_rate(
            M_sub, R_sub, acc_sharp.Sigma_pl_0, acc_sharp.v_disp_kms * 1000.0, Omega_K
        )
        @test isapprox(rate_stage1, expected_safronov_sub, rtol=1e-10)
        @test rate_stage1 > 0.0

        # Case B: Intermediate mass body (Stage 2, Pebble Accretion)
        M_mid = 10.0 * M_onset
        R_mid = 80000.0
        rate_stage2 = compute_accretion_rate(t_active, M_mid, R_mid, acc_sharp, disk_cfg)
        Sigma_peb = compute_pebble_surface_density(
            disk_cfg.orbital_distance_au;
            Sigma_peb_0=acc_sharp.Sigma_peb_0,
            p_peb=acc_sharp.p_peb,
        )
        expected_pebble_mid = compute_pebble_accretion_rate(
            M_mid,
            M_sun,
            a_test,
            Sigma_peb,
            acc_sharp.stokes_number,
            c_s,
            acc_sharp.alpha_turbulence;
            regime=:pebble_auto,
        )
        @test isapprox(rate_stage2, expected_pebble_mid, rtol=1e-10)
        @test rate_stage2 > 0.0

        # Case C: Giant embryo exceeding isolation mass (Stage 3, Safronov)
        M_super = 2.0 * M_iso
        R_super = 400000.0
        rate_stage3 = compute_accretion_rate(
            t_active, M_super, R_super, acc_sharp, disk_cfg
        )
        expected_safronov_super = compute_safronov_accretion_rate(
            M_super, R_super, acc_sharp.Sigma_pl_0, acc_sharp.v_disp_kms * 1000.0, Omega_K
        )
        @test isapprox(rate_stage3, expected_safronov_super, rtol=1e-10)
        @test rate_stage3 > 0.0

        # Custom explicit M_onset threshold override
        M_custom_onset = 5.0e19
        acc_custom = AccretionConfig(
            active=true,
            mode=:multistage,
            stage1_mode=:safronov,
            stage2_mode=:pebble_auto,
            stage3_mode=:safronov,
            M_onset=M_custom_onset,
            transition_smoothing=false,
            t_start_myr=0.0,
            t_duration_myr=2.0,
            M_target=1.0e28,
            R_target=5000000.0,
        )
        rate_below_custom = compute_accretion_rate(
            t_active, 0.5 * M_custom_onset, 20000.0, acc_custom, disk_cfg
        )
        rate_above_custom = compute_accretion_rate(
            t_active, 2.0 * M_custom_onset, 40000.0, acc_custom, disk_cfg
        )
        expected_below_custom = compute_safronov_accretion_rate(
            0.5 * M_custom_onset,
            20000.0,
            acc_custom.Sigma_pl_0,
            acc_custom.v_disp_kms * 1000.0,
            Omega_K,
        )
        expected_above_custom = compute_pebble_accretion_rate(
            2.0 * M_custom_onset,
            M_sun,
            a_test,
            Sigma_peb,
            acc_custom.stokes_number,
            c_s,
            acc_custom.alpha_turbulence;
            regime=:pebble_auto,
        )
        @test isapprox(rate_below_custom, expected_below_custom, rtol=1e-10)
        @test isapprox(rate_above_custom, expected_above_custom, rtol=1e-10)
    end

    # ---------------------------------------------------------------------
    # 6. Smoothstep Transition Rate Continuity
    # ---------------------------------------------------------------------
    @testset "Smoothstep Boundary Continuity" begin
        disk_cfg = DiskConfig(
            orbital_distance_au=2.5, stellar_mass_msun=1.0, t_ambient=150.0
        )
        sec_yr = 3.15576e7
        t_active = 0.5 * 1.0e6 * sec_yr

        acc_smooth = AccretionConfig(
            active=true,
            mode=:multistage,
            stage1_mode=:safronov,
            stage2_mode=:pebble_auto,
            stage3_mode=:safronov,
            transition_smoothing=true,
            transition_width=0.15,
            t_start_myr=0.0,
            t_duration_myr=2.0,
            M_initial=1.0e17,
            M_target=1.0e28,
            R_initial=20000.0,
            R_target=5000000.0,
        )

        M_onset = compute_pebble_onset_mass(M_sun, a_test, 0.05, c_s; f_onset=1.0)
        R_onset = cbrt(3.0 * M_onset / (4.0 * pi * 3000.0))

        # Test fine mass sampling across the onset boundary [-20%, +20%] with 101 points
        masses = [M_onset * (1.0 + eps) for eps in range(-0.20, 0.20, length=101)]
        rates = Float64[]
        for M in masses
            R = cbrt(3.0 * M / (4.0 * pi * 3000.0))
            push!(rates, compute_accretion_rate(t_active, M, R, acc_smooth, disk_cfg))
        end

        # Positivity check across transition
        @test all(r -> r > 0.0, rates)

        # Monotonicity check across onset transition window
        @test all(rates[i + 1] >= rates[i] for i in 1:(length(rates) - 1))

        # Continuity: adjacent relative difference must be bounded without step jumps
        for i in 1:(length(rates) - 1)
            rel_diff = abs(rates[i + 1] - rates[i]) / max(rates[i], rates[i + 1])
            @test rel_diff < 0.20
        end

        # Midpoint of transition window at M = M_onset should be an intermediate blend
        rate_mid = compute_accretion_rate(t_active, M_onset, R_onset, acc_smooth, disk_cfg)
        r_saf = compute_safronov_accretion_rate(
            M_onset, R_onset, acc_smooth.Sigma_pl_0, acc_smooth.v_disp_kms * 1000.0, Omega_K
        )
        Sigma_peb = compute_pebble_surface_density(
            disk_cfg.orbital_distance_au;
            Sigma_peb_0=acc_smooth.Sigma_peb_0,
            p_peb=acc_smooth.p_peb,
        )
        r_peb = compute_pebble_accretion_rate(
            M_onset,
            M_sun,
            a_test,
            Sigma_peb,
            acc_smooth.stokes_number,
            c_s,
            acc_smooth.alpha_turbulence;
            regime=:pebble_auto,
        )
        expected_blend = 0.5 * (r_saf + r_peb)
        @test isapprox(rate_mid, expected_blend, rtol=1e-3)
    end

    # ---------------------------------------------------------------------
    # 7. Edge Cases, Window Clamping, and Contract Verification
    # ---------------------------------------------------------------------
    @testset "Edge Cases and Domain Contracts" begin
        disk_cfg = DiskConfig(
            orbital_distance_au=2.5, stellar_mass_msun=1.0, t_ambient=150.0
        )
        sec_yr = 3.15576e7
        t_active = 0.5 * 1.0e6 * sec_yr

        acc_base = AccretionConfig(
            active=true,
            mode=:multistage,
            stage1_mode=:safronov,
            stage2_mode=:pebble_auto,
            stage3_mode=:safronov,
        )

        # Domain errors for non-physical states
        @test_throws DomainError compute_multistage_accretion_rate(
            t_active, -1.0, 20000.0, acc_base, disk_cfg
        )
        @test_throws DomainError compute_multistage_accretion_rate(
            t_active, 1.0e20, 0.0, acc_base, disk_cfg
        )
        @test_throws DomainError compute_multistage_accretion_rate(
            -1.0, 1.0e20, 20000.0, acc_base, disk_cfg
        )
        @test_throws DomainError compute_multistage_accretion_rate(
            t_active, NaN, 20000.0, acc_base, disk_cfg
        )

        # Error on unrecognized stage mode
        acc_bad_symbol = AccretionConfig(
            active=true, mode=:multistage, stage1_mode=:nonexistent_regime
        )
        @test_throws ArgumentError compute_multistage_accretion_rate(
            t_active, 1.0e18, 20000.0, acc_bad_symbol, disk_cfg
        )

        # M_iso <= 0 disables Stage 3: verify pebble accretion continues at very large mass
        acc_no_iso = AccretionConfig(
            active=true,
            mode=:multistage,
            stage1_mode=:safronov,
            stage2_mode=:pebble_auto,
            stage3_mode=:safronov,
            M_iso=0.0,
            M_target=1.0e30,
            R_target=1.0e8,
        )
        M_huge = 1.0e28
        R_huge = 5000000.0
        rate_huge = compute_accretion_rate(t_active, M_huge, R_huge, acc_no_iso, disk_cfg)
        Sigma_peb = compute_pebble_surface_density(
            disk_cfg.orbital_distance_au;
            Sigma_peb_0=acc_no_iso.Sigma_peb_0,
            p_peb=acc_no_iso.p_peb,
        )
        expected_pebble_huge = compute_pebble_accretion_rate(
            M_huge,
            M_sun,
            a_test,
            Sigma_peb,
            acc_no_iso.stokes_number,
            c_s,
            acc_no_iso.alpha_turbulence;
            regime=:pebble_auto,
        )
        @test isapprox(rate_huge, expected_pebble_huge, rtol=1e-10)
        @test rate_huge > 0.0

        # Close threshold smoothing: M_onset and M_iso within 5%
        # Window clamping ensures Stage 3 is reachable for M >= M_iso * (1 + w_eff)
        M_on_close = 1.0e20
        M_iso_close = 1.05e20
        acc_close = AccretionConfig(
            active=true,
            mode=:multistage,
            stage1_mode=:safronov,
            stage2_mode=:pebble_auto,
            stage3_mode=:safronov,
            M_onset=M_on_close,
            M_iso=M_iso_close,
            transition_smoothing=true,
            transition_width=0.10,
            M_target=1.0e30,
            R_target=1.0e8,
        )
        # For M far above M_iso, must evaluate strictly to stage 3
        rate_stage3_reach = compute_accretion_rate(
            t_active, 2.0 * M_iso_close, 400000.0, acc_close, disk_cfg
        )
        expected_stage3_reach = compute_safronov_accretion_rate(
            2.0 * M_iso_close,
            400000.0,
            acc_close.Sigma_pl_0,
            acc_close.v_disp_kms * 1000.0,
            Omega_K,
        )
        @test isapprox(rate_stage3_reach, expected_stage3_reach, rtol=1e-10)
        @test rate_stage3_reach > 0.0

        # Direct transition when M_iso <= M_onset (bypassing Stage 2)
        acc_direct = AccretionConfig(
            active=true,
            mode=:multistage,
            stage1_mode=:safronov,
            stage2_mode=:pebble_auto,
            stage3_mode=:safronov,
            M_onset=2.0e20,
            M_iso=1.0e20,
            transition_smoothing=false,
            M_target=1.0e30,
            R_target=1.0e8,
        )
        rate_direct_sub = compute_accretion_rate(
            t_active, 0.5e20, 20000.0, acc_direct, disk_cfg
        )
        rate_direct_super = compute_accretion_rate(
            t_active, 1.5e20, 40000.0, acc_direct, disk_cfg
        )
        exp_saf_sub = compute_safronov_accretion_rate(
            0.5e20, 20000.0, acc_direct.Sigma_pl_0, acc_direct.v_disp_kms * 1000.0, Omega_K
        )
        exp_saf_super = compute_safronov_accretion_rate(
            1.5e20, 40000.0, acc_direct.Sigma_pl_0, acc_direct.v_disp_kms * 1000.0, Omega_K
        )
        @test isapprox(rate_direct_sub, exp_saf_sub, rtol=1e-10)
        @test isapprox(rate_direct_super, exp_saf_super, rtol=1e-10)

        # TOML string parsing of "NaN" and unquoted nan
        toml_str = """
        [accretion]
        active = true
        mode = "multistage"
        M_onset = "NaN"
        M_iso = nan
        """
        parsed_cfg = parse_config_string(toml_str)
        @test isnan(parsed_cfg.accretion.M_onset)
        @test isnan(parsed_cfg.accretion.M_iso)
    end
end

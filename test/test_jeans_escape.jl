using Test
using TOML
using Erebus
using Erebus.Config
using Erebus.Physics

@testset "Atmospheric Jeans Kinetic Escape & Volatile Loss" begin
    @testset "Molecular Masses & Physical Constants" begin
        # Fundamental physical constants
        G_const = 6.67430e-11
        kB_const = 1.380649e-23
        NA_const = 6.02214076e23

        # Species molecular masses against standard atomic weights [kg]
        @test isapprox(MASS_H2O_KG, 18.01528e-3 / NA_const; rtol=1e-5)
        @test isapprox(MASS_N2_KG, 28.0134e-3 / NA_const; rtol=1e-5)
        @test isapprox(MASS_NH3_KG, 17.03052e-3 / NA_const; rtol=1e-5)
        @test isapprox(MASS_CO_KG, 28.0101e-3 / NA_const; rtol=1e-5)
        @test isapprox(MASS_CO2_KG, 44.0095e-3 / NA_const; rtol=1e-5)

        # Species lookup function
        @test get_species_molecular_mass(:H2O) == MASS_H2O_KG
        @test get_species_molecular_mass(:h2o) == MASS_H2O_KG
        @test get_species_molecular_mass(:N2) == MASS_N2_KG
        @test get_species_molecular_mass(:n2) == MASS_N2_KG
        @test get_species_molecular_mass(:NH3) == MASS_NH3_KG
        @test get_species_molecular_mass(:CO) == MASS_CO_KG
        @test get_species_molecular_mass(:CO2) == MASS_CO2_KG
        @test get_species_molecular_mass(:CH4) == MASS_CH4_KG
        @test get_species_molecular_mass(:H2) == MASS_H2_KG
        @test get_species_molecular_mass(:H2S) == MASS_H2S_KG
        @test get_species_molecular_mass(:S2) == MASS_S2_KG
        @test get_species_molecular_mass(:SO2) == MASS_SO2_KG

        @test_throws ArgumentError get_species_molecular_mass(:xenon)
        @test_throws ArgumentError get_species_molecular_mass(:argon)
    end

    @testset "compute_escape_velocity Invariants" begin
        # Earth escape velocity: M = 5.9722e24 kg, R = 6.371e6 m -> ~11.186 km/s
        v_esc_earth = compute_escape_velocity(5.9722e24, 6.371e6)
        @test isapprox(v_esc_earth, 11186.0; rtol=1e-3)

        # 50 km planetesimal (rho = 2500 kg/m^3 -> M = 4/3 * pi * R^3 * rho ≈ 1.309e18 kg)
        R_planetesimal = 50_000.0
        M_planetesimal = (4.0 / 3.0) * π * R_planetesimal^3 * 2500.0
        v_esc_planetesimal = compute_escape_velocity(M_planetesimal, R_planetesimal)
        # Expected: sqrt(2 * 6.6743e-11 * 1.309e18 / 50000) ≈ 59.08 m/s
        @test isapprox(v_esc_planetesimal, 59.08; rtol=1e-2)

        # Scaling: quadrupling mass doubles escape velocity
        @test isapprox(
            compute_escape_velocity(4.0 * M_planetesimal, R_planetesimal),
            2.0 * v_esc_planetesimal;
            rtol=1e-12,
        )

        # Scaling: quadrupling radius halves escape velocity
        @test isapprox(
            compute_escape_velocity(M_planetesimal, 4.0 * R_planetesimal),
            0.5 * v_esc_planetesimal;
            rtol=1e-12,
        )

        # DomainError guards
        @test_throws DomainError compute_escape_velocity(0.0, R_planetesimal)
        @test_throws DomainError compute_escape_velocity(-1.0e20, R_planetesimal)
        @test_throws DomainError compute_escape_velocity(M_planetesimal, 0.0)
        @test_throws DomainError compute_escape_velocity(M_planetesimal, -500.0)
        @test_throws DomainError compute_escape_velocity(NaN, R_planetesimal)
        @test_throws DomainError compute_escape_velocity(M_planetesimal, NaN)
    end

    @testset "compute_thermal_velocity Invariants" begin
        # Water at T = 200 K
        T_200 = 200.0
        v_th_h2o = compute_thermal_velocity(T_200, MASS_H2O_KG)
        # Expected: sqrt(2 * 1.380649e-23 * 200 / 2.991507e-26) ≈ 429.69 m/s
        @test isapprox(v_th_h2o, 429.69; rtol=1e-3)

        # Nitrogen at T = 200 K
        v_th_n2 = compute_thermal_velocity(T_200, MASS_N2_KG)
        # Heavier molecule has lower thermal speed
        @test v_th_n2 < v_th_h2o
        @test isapprox(v_th_n2, 344.56; rtol=1e-3)

        # Scaling: quadrupling temperature doubles thermal speed
        @test isapprox(
            compute_thermal_velocity(4.0 * T_200, MASS_H2O_KG), 2.0 * v_th_h2o; rtol=1e-12
        )

        # DomainError guards
        @test_throws DomainError compute_thermal_velocity(0.0, MASS_H2O_KG)
        @test_throws DomainError compute_thermal_velocity(-100.0, MASS_H2O_KG)
        @test_throws DomainError compute_thermal_velocity(T_200, 0.0)
        @test_throws DomainError compute_thermal_velocity(T_200, -1.0e-26)
        @test_throws DomainError compute_thermal_velocity(NaN, MASS_H2O_KG)
        @test_throws DomainError compute_thermal_velocity(T_200, NaN)
    end

    @testset "compute_jeans_parameter Invariants" begin
        R_50km = 50_000.0
        M_50km = (4.0 / 3.0) * π * R_50km^3 * 2500.0
        T_200 = 200.0

        lambda_h2o_small = compute_jeans_parameter(M_50km, R_50km, T_200, MASS_H2O_KG)
        # Expected: (v_esc / v_th)^2 = (59.08 / 429.69)^2 ≈ 0.0189 << 1
        @test isapprox(lambda_h2o_small, 0.0189; rtol=2e-2)

        # Exact identity check: lambda == (v_esc / v_th)^2
        v_esc = compute_escape_velocity(M_50km, R_50km)
        v_th = compute_thermal_velocity(T_200, MASS_H2O_KG)
        @test isapprox(lambda_h2o_small, (v_esc / v_th)^2; rtol=1e-12)

        # 500 km protoplanet (R 10x larger, M 1000x larger -> lambda 100x larger)
        R_500km = 500_000.0
        M_500km = (4.0 / 3.0) * π * R_500km^3 * 2500.0
        lambda_h2o_mid = compute_jeans_parameter(M_500km, R_500km, T_200, MASS_H2O_KG)
        @test isapprox(lambda_h2o_mid, 100.0 * lambda_h2o_small; rtol=1e-10)
        @test isapprox(lambda_h2o_mid, 1.89; rtol=2e-2)

        # DomainError guards
        @test_throws DomainError compute_jeans_parameter(0.0, R_50km, T_200, MASS_H2O_KG)
        @test_throws DomainError compute_jeans_parameter(M_50km, 0.0, T_200, MASS_H2O_KG)
        @test_throws DomainError compute_jeans_parameter(M_50km, R_50km, 0.0, MASS_H2O_KG)
        @test_throws DomainError compute_jeans_parameter(M_50km, R_50km, T_200, 0.0)
    end

    @testset "compute_jeans_escape_flux & Mass Loss Rate Invariants" begin
        T_200 = 200.0
        v_th = compute_thermal_velocity(T_200, MASS_H2O_KG)
        n_exo = 1.0e18 # m^-3

        # Effusion limit: lambda -> 0
        flux_effusion = compute_jeans_escape_flux(
            n_exo, T_200, MASS_H2O_KG, 0.0; hydrodynamic=false
        )
        # Expected: n * v_th / (2 * sqrt(pi))
        expected_effusion = n_exo * v_th / (2.0 * sqrt(π))
        @test isapprox(flux_effusion, expected_effusion; rtol=1e-12)

        # Gravitational suppression: lambda = 5.0
        flux_lambda5 = compute_jeans_escape_flux(
            n_exo, T_200, MASS_H2O_KG, 5.0; hydrodynamic=false
        )
        # Factor: (1 + 5) * exp(-5) ≈ 6 * 0.0067379 = 0.0404
        @test isapprox(flux_lambda5 / flux_effusion, 6.0 * exp(-5.0); rtol=1e-12)
        @test flux_lambda5 < flux_effusion

        # Strong gravitational suppression: lambda = 120 (Earth-like or massive)
        flux_suppressed = compute_jeans_escape_flux(n_exo, T_200, MASS_H2O_KG, 120.0)
        @test iszero(flux_suppressed)

        # Zero density gives zero flux
        @test iszero(compute_jeans_escape_flux(0.0, T_200, MASS_H2O_KG, 2.0))
        @test iszero(compute_jeans_escape_flux(-1.0e10, T_200, MASS_H2O_KG, 2.0))

        # Total planetary mass loss rate: 4 * pi * R^2 * flux * m
        R_50km = 50_000.0
        rho_exo = n_exo * MASS_H2O_KG
        M_50km = (4.0 / 3.0) * π * R_50km^3 * 2500.0
        loss_rate = compute_jeans_mass_loss_rate(
            M_50km, R_50km, T_200, MASS_H2O_KG, rho_exo; hydrodynamic=false
        )
        lambda_val = compute_jeans_parameter(M_50km, R_50km, T_200, MASS_H2O_KG)
        expected_loss =
            4.0 *
            π *
            R_50km^2 *
            flux_effusion *
            (1.0 + lambda_val) *
            exp(-lambda_val) *
            MASS_H2O_KG
        @test isapprox(loss_rate, expected_loss; rtol=1e-10)

        # Hydrodynamic blow-off regime when lambda < 1.5 and hydrodynamic=true
        c_s = sqrt(1.4 * 1.380649e-23 * T_200 / MASS_H2O_KG)
        flux_hydro = compute_jeans_escape_flux(
            n_exo, T_200, MASS_H2O_KG, 0.0; hydrodynamic=true
        )
        @test isapprox(flux_hydro, n_exo * c_s; rtol=1e-12)
        @test flux_hydro > flux_effusion

        loss_rate_hydro = compute_jeans_mass_loss_rate(
            M_50km, R_50km, T_200, MASS_H2O_KG, rho_exo; hydrodynamic=true
        )
        expected_loss_hydro = 4.0 * π * R_50km^2 * rho_exo * c_s
        @test isapprox(loss_rate_hydro, expected_loss_hydro; rtol=1e-10)

        # Zero density gives zero loss rate
        @test iszero(compute_jeans_mass_loss_rate(M_50km, R_50km, T_200, MASS_H2O_KG, 0.0))

        # Hermite spline transition continuity across lambda in [1.0, 2.0]
        flux_099 = compute_jeans_escape_flux(
            n_exo, T_200, MASS_H2O_KG, 0.999; hydrodynamic=true
        )
        flux_101 = compute_jeans_escape_flux(
            n_exo, T_200, MASS_H2O_KG, 1.001; hydrodynamic=true
        )
        @test isapprox(flux_099, flux_101; rtol=1e-2)

        flux_199 = compute_jeans_escape_flux(
            n_exo, T_200, MASS_H2O_KG, 1.999; hydrodynamic=true
        )
        flux_201 = compute_jeans_escape_flux(
            n_exo, T_200, MASS_H2O_KG, 2.001; hydrodynamic=true
        )
        @test isapprox(flux_199, flux_201; rtol=1e-2)

        # Monotonicity check across transition: flux decreases as gravity increases
        flux_10 = compute_jeans_escape_flux(
            n_exo, T_200, MASS_H2O_KG, 1.0; hydrodynamic=true
        )
        flux_15 = compute_jeans_escape_flux(
            n_exo, T_200, MASS_H2O_KG, 1.5; hydrodynamic=true
        )
        flux_20 = compute_jeans_escape_flux(
            n_exo, T_200, MASS_H2O_KG, 2.0; hydrodynamic=true
        )
        @test flux_10 > flux_15 > flux_20

        # Exact consistency between compute_jeans_escape_flux and compute_jeans_mass_loss_rate
        for test_R in [20_000.0, 50_000.0, 100_000.0, 500_000.0]
            m_loss = compute_jeans_mass_loss_rate(
                M_50km, test_R, T_200, MASS_H2O_KG, rho_exo; hydrodynamic=true
            )
            lam_test = compute_jeans_parameter(M_50km, test_R, T_200, MASS_H2O_KG)
            f_test = compute_jeans_escape_flux(
                n_exo, T_200, MASS_H2O_KG, lam_test; hydrodynamic=true
            )
            @test isapprox(m_loss, 4.0 * π * test_R^2 * f_test * MASS_H2O_KG; rtol=1e-12)
        end
    end

    @testset "Atmospheric Scale Height & Surface Pressure" begin
        R_50km = 50_000.0
        M_50km = (4.0 / 3.0) * π * R_50km^3 * 2500.0
        T_200 = 200.0

        # Scale height: H = k_B * T / (m * g)
        H_h2o = compute_atmospheric_scale_height(M_50km, R_50km, T_200, MASS_H2O_KG)
        # g = 6.6743e-11 * 1.309e18 / (50000)^2 ≈ 0.03495 m/s^2
        # H = 1.380649e-23 * 200 / (2.991507e-26 * 0.03495) ≈ 2.641e6 m ≈ 2641 km >> R_planet
        @test H_h2o > R_50km

        # Surface pressure from total atmospheric mass
        M_atm_test = 1.0e12 # kg
        P_surf = compute_surface_atmospheric_pressure(M_atm_test, M_50km, R_50km)
        @test P_surf > 0.0
        # Linearity: doubling atmospheric mass doubles surface pressure
        @test isapprox(
            compute_surface_atmospheric_pressure(2.0 * M_atm_test, M_50km, R_50km),
            2.0 * P_surf;
            rtol=1e-12,
        )
        # Zero atmospheric mass gives zero pressure
        @test iszero(compute_surface_atmospheric_pressure(0.0, M_50km, R_50km))

        # DomainError guards
        @test_throws DomainError compute_atmospheric_scale_height(
            0.0, R_50km, T_200, MASS_H2O_KG
        )
        @test_throws DomainError compute_surface_atmospheric_pressure(
            M_atm_test, 0.0, R_50km
        )
        @test_throws DomainError compute_surface_atmospheric_pressure(
            M_atm_test, M_50km, 0.0
        )
        @test_throws DomainError compute_surface_atmospheric_pressure(0.0, 0.0, R_50km)
        @test_throws DomainError compute_surface_atmospheric_pressure(0.0, M_50km, 0.0)
        @test_throws DomainError compute_surface_atmospheric_pressure(-1.0, M_50km, R_50km)
    end

    @testset "evolve_atmospheric_species_inventory Dynamics & Conservation" begin
        R_50km = 50_000.0
        M_50km = (4.0 / 3.0) * π * R_50km^3 * 2500.0
        T_200 = 200.0
        dt_1yr = 365.25 * 86400.0 # seconds in 1 year

        # Case 1: Pure decay (zero venting influx)
        M_init = 1.0e10 # kg
        res_decay = evolve_atmospheric_species_inventory(
            M_init, 0.0, dt_1yr, M_50km, R_50km, T_200, MASS_H2O_KG
        )
        @test res_decay.M_atm < M_init
        @test res_decay.M_escaped_step > 0.0
        # Exact mass conservation: M_atm + M_escaped == M_init
        @test isapprox(res_decay.M_atm + res_decay.M_escaped_step, M_init; rtol=1e-12)

        # Case 2: Steady venting influx into small body (pure effusion)
        vent_rate = 100.0 # kg / s
        res_vent = evolve_atmospheric_species_inventory(
            0.0, vent_rate, dt_1yr, M_50km, R_50km, T_200, MASS_H2O_KG; hydrodynamic=false
        )
        total_vented = vent_rate * dt_1yr
        @test isapprox(res_vent.M_atm + res_vent.M_escaped_step, total_vented; rtol=1e-12)
        # On a 50 km body, escape is fast, so most vented mass escapes
        @test res_vent.M_escaped_step > res_vent.M_atm

        # Case 2b: Steady venting influx into small body with hydrodynamic blowoff
        res_vent_hydro = evolve_atmospheric_species_inventory(
            0.0, vent_rate, dt_1yr, M_50km, R_50km, T_200, MASS_H2O_KG; hydrodynamic=true
        )
        @test isapprox(
            res_vent_hydro.M_atm + res_vent_hydro.M_escaped_step, total_vented; rtol=1e-12
        )
        @test res_vent_hydro.M_escaped_step > res_vent_hydro.M_atm
        c_s_h2o = sqrt(1.4 * 1.380649e-23 * T_200 / MASS_H2O_KG)
        k_hydro = c_s_h2o / R_50km
        @test isapprox(res_vent_hydro.M_atm, vent_rate / k_hydro; rtol=1e-5)

        # Case 3: Massive body retention (R = 5000 km, M = 1.3e24 kg)
        R_massive = 5_000_000.0
        M_massive = (4.0 / 3.0) * π * R_massive^3 * 3000.0
        res_massive = evolve_atmospheric_species_inventory(
            0.0, vent_rate, dt_1yr, M_massive, R_massive, T_200, MASS_H2O_KG
        )
        # Retained: almost all vented mass remains in the atmosphere
        @test isapprox(res_massive.M_atm, total_vented; rtol=1e-5)
        @test isapprox(res_massive.M_escaped_step, 0.0; atol=1e-4)
        @test isapprox(
            res_massive.M_atm + res_massive.M_escaped_step, total_vented; rtol=1e-12
        )

        # Zero dt leaves state unchanged
        res_zero_dt = evolve_atmospheric_species_inventory(
            M_init, vent_rate, 0.0, M_50km, R_50km, T_200, MASS_H2O_KG
        )
        @test isapprox(res_zero_dt.M_atm, M_init; rtol=1e-12)
        @test iszero(res_zero_dt.M_escaped_step)

        # Quantitative closed-form check in the small planetesimal regime (λ < 1)
        g_50 = 6.67430e-11 * M_50km / (R_50km^2)
        H_50 = (1.380649e-23 * T_200) / (MASS_H2O_KG * g_50)
        v_th_50 = sqrt(2.0 * 1.380649e-23 * T_200 / MASS_H2O_KG)
        lam_50 = compute_jeans_parameter(M_50km, R_50km, T_200, MASS_H2O_KG)
        k_expected = (v_th_50 / (2.0 * sqrt(π) * H_50)) * (1.0 + lam_50) * exp(-lam_50)
        @test isapprox(lam_50, 0.01893002; rtol=1e-5)
        @test isapprox(k_expected, 4.5880267e-5; rtol=1e-5)

        # Steady-state atmospheric mass under continuous venting: M_ss = M_dot / k_escape
        M_ss_expected = vent_rate / k_expected
        @test isapprox(res_vent.M_atm, M_ss_expected; rtol=1e-5)

        # R_exobase keyword argument support and scaling
        R_exo_high = 1.2 * R_50km
        res_exo = evolve_atmospheric_species_inventory(
            0.0, vent_rate, dt_1yr, M_50km, R_50km, T_200, MASS_H2O_KG; R_exobase=R_exo_high
        )
        @test res_exo.M_atm > 0.0
        @test isapprox(res_exo.M_atm + res_exo.M_escaped_step, total_vented; rtol=1e-12)
        @test res_exo.M_atm != res_vent.M_atm

        # Guard: R_exobase < R_planet throws DomainError
        @test_throws DomainError evolve_atmospheric_species_inventory(
            0.0,
            vent_rate,
            dt_1yr,
            M_50km,
            R_50km,
            T_200,
            MASS_H2O_KG;
            R_exobase=0.8 * R_50km,
        )

        # DomainError guards
        @test_throws DomainError evolve_atmospheric_species_inventory(
            -1.0, vent_rate, dt_1yr, M_50km, R_50km, T_200, MASS_H2O_KG
        )
        @test_throws DomainError evolve_atmospheric_species_inventory(
            M_init, -10.0, dt_1yr, M_50km, R_50km, T_200, MASS_H2O_KG
        )
        @test_throws DomainError evolve_atmospheric_species_inventory(
            M_init, vent_rate, -100.0, M_50km, R_50km, T_200, MASS_H2O_KG
        )
    end

    @testset "EscapeConfig Schema & Bounds Validation" begin
        esc_default = EscapeConfig()
        @test !esc_default.active
        @test esc_default.M_planet > 0.0
        @test esc_default.R_planet > 0.0
        @test esc_default.T_exobase > 0.0
        @test esc_default.R_exobase >= esc_default.R_planet

        # TOML serialization round-trip
        quick_toml = joinpath(@__DIR__, "..", "configs", "test_quick.toml")
        cfg = load_config(quick_toml)
        @test !cfg.escape.active

        # Modify escape section and verify round-trip
        toml_content = read(quick_toml, String)
        escape_block = """
        [geometry]
        rplanet = 60000.0

        [escape]
        active = true
        M_planet = 2.5e18
        R_planet = 60000.0
        T_exobase = 220.0
        R_exobase = 65000.0
        """
        full_toml = toml_content * "\n" * escape_block
        cfg_loaded = load_config(full_toml)
        @test cfg_loaded.escape.active
        @test isapprox(cfg_loaded.geometry.rplanet, 60000.0; rtol=1e-12)
        @test isapprox(cfg_loaded.escape.M_planet, 2.5e18; rtol=1e-12)
        @test isapprox(cfg_loaded.escape.R_planet, 60000.0; rtol=1e-12)
        @test isapprox(cfg_loaded.escape.T_exobase, 220.0; rtol=1e-12)
        @test isapprox(cfg_loaded.escape.R_exobase, 65000.0; rtol=1e-12)

        # Save and re-load round-trip
        io = IOBuffer()
        save_config(io, cfg_loaded)
        saved_str = String(take!(io))
        cfg_reloaded = load_config(saved_str)
        @test cfg_reloaded.escape.active == cfg_loaded.escape.active
        @test isapprox(
            cfg_reloaded.geometry.rplanet, cfg_loaded.geometry.rplanet; rtol=1e-12
        )
        @test isapprox(cfg_reloaded.escape.M_planet, cfg_loaded.escape.M_planet; rtol=1e-12)
        @test isapprox(cfg_reloaded.escape.R_planet, cfg_loaded.escape.R_planet; rtol=1e-12)
        @test isapprox(
            cfg_reloaded.escape.T_exobase, cfg_loaded.escape.T_exobase; rtol=1e-12
        )
        @test isapprox(
            cfg_reloaded.escape.R_exobase, cfg_loaded.escape.R_exobase; rtol=1e-12
        )

        # Validation error bounds
        cfg_bad_mass = SimulationConfig(; escape=EscapeConfig(; M_planet=-1.0))
        @test_throws ArgumentError validate_config(cfg_bad_mass)

        cfg_bad_radius = SimulationConfig(; escape=EscapeConfig(; R_planet=0.0))
        @test_throws ArgumentError validate_config(cfg_bad_radius)

        cfg_bad_T = SimulationConfig(; escape=EscapeConfig(; T_exobase=-50.0))
        @test_throws ArgumentError validate_config(cfg_bad_T)

        cfg_bad_R_exo = SimulationConfig(;
            escape=EscapeConfig(; R_planet=50000.0, R_exobase=40000.0)
        )
        @test_throws ArgumentError validate_config(cfg_bad_R_exo)

        cfg_mismatched_radius = SimulationConfig(;
            geometry=GeometryConfig(; rplanet=50000.0),
            escape=EscapeConfig(; active=true, R_planet=60000.0),
        )
        @test_throws ArgumentError validate_config(cfg_mismatched_radius)

        cfg_bad_gamma = SimulationConfig(; escape=EscapeConfig(; gamma=0.0))
        @test_throws ArgumentError validate_config(cfg_bad_gamma)

        cfg_bad_species = SimulationConfig(;
            escape=EscapeConfig(; species=:unsupported_gas)
        )
        @test_throws ArgumentError validate_config(cfg_bad_species)

        cfg_bad_species_list = SimulationConfig(;
            escape=EscapeConfig(; species_list=[:H2O, :invalid_gas])
        )
        @test_throws ArgumentError validate_config(cfg_bad_species_list)

        cfg_mismatched_venting_escape_multi = SimulationConfig(;
            venting=VentingConfig(; active=true, species=:CO2),
            escape=EscapeConfig(;
                active=true, multi_species=true, species_list=[:H2O, :CO]
            ),
        )
        @test_throws ArgumentError validate_config(cfg_mismatched_venting_escape_multi)

        cfg_mismatched_venting_escape_single = SimulationConfig(;
            venting=VentingConfig(; active=true, species=:CO2),
            escape=EscapeConfig(; active=true, multi_species=false, species=:H2O),
        )
        @test_throws ArgumentError validate_config(cfg_mismatched_venting_escape_single)
    end

    @testset "Simulation Loop Integration with Atmospheric Escape" begin
        output_dir = mktempdir()
        try
            quick_toml = joinpath(@__DIR__, "..", "configs", "test_quick.toml")
            cfg = load_config(quick_toml)

            # 2 steps with escape.active = true and venting.active = false (ungated escape test)
            cfg_escape = SimulationConfig(
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
                output=OutputConfig(output_dir=output_dir, savematstep=2),
                disk=cfg.disk,
                melting=cfg.melting,
                venting=VentingConfig(active=false),
                escape=EscapeConfig(
                    active=true,
                    M_planet=1.309e18,
                    R_planet=50_000.0,
                    T_exobase=200.0,
                    R_exobase=50_000.0,
                ),
            )
            Erebus.simulation_loop(cfg_escape; output_path=output_dir)

            files = readdir(output_dir)
            @test "output_00000.jld2" in files
            @test "output_00002.jld2" in files

            data2 = load_state(joinpath(output_dir, "output_00002.jld2"))
            @test data2["timestep"] == 2
            @test haskey(data2, "M_atm_total")
            @test haskey(data2, "M_escaped_total")
            @test isapprox(data2["M_atm_total"], 0.0; atol=1e-12)
            @test isapprox(data2["M_escaped_total"], 0.0; atol=1e-12)

            # Test restart/resume preserving atmospheric inventory fields
            cfg_resume = SimulationConfig(
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
                    start_step=2,
                    n_steps=3,
                ),
                solver=cfg.solver,
                poroelasticity=cfg.poroelasticity,
                thermodynamics=cfg.thermodynamics,
                reaction=cfg.reaction,
                materials=cfg.materials,
                output=OutputConfig(output_dir=output_dir, savematstep=1),
                disk=cfg.disk,
                melting=cfg.melting,
                venting=VentingConfig(active=false),
                escape=EscapeConfig(
                    active=true,
                    M_planet=1.309e18,
                    R_planet=50_000.0,
                    T_exobase=200.0,
                    R_exobase=50_000.0,
                ),
            )
            Erebus.simulation_loop(
                cfg_resume;
                output_path=output_dir,
                restart_from=joinpath(output_dir, "output_00002.jld2"),
            )
            data3 = load_state(joinpath(output_dir, "output_00003.jld2"))
            @test data3["timestep"] == 3
            @test haskey(data3, "M_atm_total")
            @test haskey(data3, "M_escaped_total")
            @test isapprox(data3["M_atm_total"], 0.0; atol=1e-12)
            @test isapprox(data3["M_escaped_total"], 0.0; atol=1e-12)
        finally
            rm(output_dir; recursive=true, force=true)
        end

        # Active coupled venting and atmospheric escape test
        output_dir_coupled = mktempdir()
        try
            quick_toml = joinpath(@__DIR__, "..", "configs", "test_quick.toml")
            cfg = load_config(quick_toml)
            cfg_coupled = SimulationConfig(
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
                output=OutputConfig(output_dir=output_dir_coupled, savematstep=2),
                disk=cfg.disk,
                melting=cfg.melting,
                venting=VentingConfig(active=true),
                escape=EscapeConfig(
                    active=true,
                    M_planet=1.309e18,
                    R_planet=50_000.0,
                    T_exobase=200.0,
                    R_exobase=50_000.0,
                ),
            )
            Erebus.simulation_loop(cfg_coupled; output_path=output_dir_coupled)
            data_c = load_state(joinpath(output_dir_coupled, "output_00002.jld2"))
            @test data_c["timestep"] == 2
            @test haskey(data_c, "M_atm_total")
            @test haskey(data_c, "M_escaped_total")
            @test haskey(data_c, "M_vent_total")
            @test data_c["M_vent_total"] >= 0.0
            @test data_c["M_atm_total"] >= 0.0
            @test data_c["M_escaped_total"] >= 0.0
            # Exact mass balance: 3D vented mass equals retained atmospheric mass plus escaped mass
            @test isapprox(
                data_c["M_atm_total"] + data_c["M_escaped_total"],
                data_c["M_vent_total"];
                rtol=1e-10,
                atol=1e-12,
            )
        finally
            rm(output_dir_coupled; recursive=true, force=true)
        end
    end
end

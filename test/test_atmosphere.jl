using Test
using Erebus
using StaticArrays
using LinearAlgebra
using Random

@testset "Coupled 1D Atmosphere, Disk Envelope, and Escape Physics" begin
    M_sun = 1.98847e30
    AU = 1.495978707e11
    a_test = 2.5 * AU
    T_disk = 150.0
    c_s = compute_sound_speed(T_disk)
    R_50km = 50_000.0
    M_50km = (4.0 / 3.0) * pi * (R_50km^3) * 3000.0
    g_50km = Erebus.G_GRAV * M_50km / (R_50km^2)

    # ---------------------------------------------------------------------
    # 1. Gravitational Capture Radius (Bondi vs Hill)
    # ---------------------------------------------------------------------
    @testset "Gravitational Capture Radius Physics" begin
        # Analytical Bondi and Hill radii
        R_Bondi = Erebus.G_GRAV * M_50km / (c_s^2)
        R_Hill = a_test * cbrt(M_50km / (3.0 * M_sun))
        R_cap = compute_gravitational_capture_radius(M_50km, M_sun, a_test, c_s)

        # Basic properties
        @test isapprox(R_cap, min(R_Bondi, R_Hill), rtol=1e-12)
        @test R_cap > 0.0
        @test R_cap < a_test

        # Regime 1: High sound speed -> Bondi-limited
        c_s_hot = 5000.0
        R_cap_hot = compute_gravitational_capture_radius(M_50km, M_sun, a_test, c_s_hot)
        R_Bondi_hot = Erebus.G_GRAV * M_50km / (c_s_hot^2)
        @test isapprox(R_cap_hot, R_Bondi_hot, rtol=1e-12)
        @test R_cap_hot < R_Hill

        # Regime 2: Low sound speed -> Hill-limited
        c_s_cold = 0.5
        R_cap_cold = compute_gravitational_capture_radius(M_50km, M_sun, a_test, c_s_cold)
        @test isapprox(R_cap_cold, R_Hill, rtol=1e-12)

        # Monotonicity with planet mass
        R_cap_2M = compute_gravitational_capture_radius(2.0 * M_50km, M_sun, a_test, c_s)
        @test R_cap_2M > R_cap

        # Error guards on non-positive or non-finite inputs
        @test_throws DomainError compute_gravitational_capture_radius(
            0.0, M_sun, a_test, c_s
        )
        @test_throws DomainError compute_gravitational_capture_radius(
            -M_50km, M_sun, a_test, c_s
        )
        @test_throws DomainError compute_gravitational_capture_radius(
            M_50km, 0.0, a_test, c_s
        )
        @test_throws DomainError compute_gravitational_capture_radius(
            M_50km, M_sun, 0.0, c_s
        )
        @test_throws DomainError compute_gravitational_capture_radius(
            M_50km, M_sun, a_test, 0.0
        )
        @test_throws DomainError compute_gravitational_capture_radius(
            NaN, M_sun, a_test, c_s
        )
    end

    # ---------------------------------------------------------------------
    # 2. Disk Gas Envelope Mass & Recycling Limit (Ormel et al. 2015)
    # ---------------------------------------------------------------------
    @testset "Disk Gas Envelope Mass and Recycling Limit" begin
        # For a massive embryo with R_cap > R_planet
        M_embryo = 1.0e23 # kg
        R_embryo = 1_500_000.0 # m (1500 km)
        c_s_test = 300.0 # m/s
        R_cap = compute_gravitational_capture_radius(M_embryo, M_sun, a_test, c_s_test)
        @test R_cap > R_embryo

        rho_disk = 1.0e-9 # kg/m^3
        f_rec = 0.10

        M_env = compute_disk_envelope_mass(
            M_embryo, R_embryo, R_cap, rho_disk, c_s_test; f_rec=f_rec
        )
        @test M_env > 0.0

        # Capped by Ormel et al. (2015) recycling limit
        M_rec_max = f_rec * (4.0 * pi / 3.0) * (R_cap^3) * rho_disk
        @test M_env <= M_rec_max * (1.0 + 1e-12)

        # Asymptotic limit: Zero disk gas density yields zero envelope mass
        @test iszero(compute_disk_envelope_mass(M_embryo, R_embryo, R_cap, 0.0, c_s_test))

        # Asymptotic limit: If planet radius equals or exceeds capture radius, envelope is zero
        @test iszero(compute_disk_envelope_mass(M_embryo, R_cap, R_cap, rho_disk, c_s_test))
        @test iszero(
            compute_disk_envelope_mass(M_embryo, 2.0 * R_cap, R_cap, rho_disk, c_s_test)
        )

        # Monotonicity with disk density
        M_env_2rho = compute_disk_envelope_mass(
            M_embryo, R_embryo, R_cap, 2.0 * rho_disk, c_s_test; f_rec=f_rec
        )
        @test M_env_2rho > M_env

        # Error guards
        @test_throws DomainError compute_disk_envelope_mass(
            -1.0, R_50km, R_cap, rho_disk, c_s
        )
        @test_throws DomainError compute_disk_envelope_mass(
            M_50km, -1.0, R_cap, rho_disk, c_s
        )
        @test_throws DomainError compute_disk_envelope_mass(
            M_50km, R_50km, -1.0, rho_disk, c_s
        )
        @test_throws DomainError compute_disk_envelope_mass(
            M_50km, R_50km, R_cap, -1.0, c_s
        )
        @test_throws DomainError compute_disk_envelope_mass(
            M_50km, R_50km, R_cap, rho_disk, 0.0
        )
        @test_throws DomainError compute_disk_envelope_mass(
            M_50km, R_50km, R_cap, rho_disk, c_s; f_rec=-0.1
        )
    end

    # ---------------------------------------------------------------------
    # 3. Atmospheric Optical Depth & Species Weighting
    # ---------------------------------------------------------------------
    @testset "Atmospheric Optical Depth & Species Weighting" begin
        opacities = Dict(:H2O => 1.0e-2, :CO2 => 1.0e-3, :CH4 => 2.0e-3)
        M_atm = Dict(:H2O => 1.0e14, :CO2 => 2.0e14, :CH4 => 0.5e14)
        area = 4.0 * pi * (R_50km^2)

        tau_LW = compute_atmospheric_optical_depth(M_atm, R_50km, opacities)
        tau_expected = (1.0e-2 * 1.0e14 + 1.0e-3 * 2.0e14 + 2.0e-3 * 0.5e14) / area

        # Assert precision match
        @test isapprox(tau_LW, tau_expected, rtol=1e-12)
        @test tau_LW > 0.0

        # Asymptotic limit: Zero atmospheric inventory yields zero optical depth
        M_zero = Dict(:H2O => 0.0, :CO2 => 0.0)
        @test iszero(compute_atmospheric_optical_depth(M_zero, R_50km, opacities))

        # Linear scaling with mass
        M_double = Dict(k => 2.0 * v for (k, v) in M_atm)
        @test isapprox(
            compute_atmospheric_optical_depth(M_double, R_50km, opacities),
            2.0 * tau_LW,
            rtol=1e-12,
        )

        # Inverse square scaling with radius: R -> 2R => area -> 4 area => tau -> tau / 4
        @test isapprox(
            compute_atmospheric_optical_depth(M_atm, 2.0 * R_50km, opacities),
            0.25 * tau_LW,
            rtol=1e-12,
        )

        # Fallback to kappa_default for unlisted species
        opacities_partial = Dict(:H2O => 1.0e-2)
        M_unlisted = Dict(:H2O => 1.0e14, :CO => 1.0e14)
        tau_unlisted = compute_atmospheric_optical_depth(
            M_unlisted, R_50km, opacities_partial; kappa_default=1.0e-3
        )
        tau_unlisted_expected = (1.0e-2 * 1.0e14 + 1.0e-3 * 1.0e14) / area
        @test isapprox(tau_unlisted, tau_unlisted_expected, rtol=1e-12)

        # Error guards
        @test_throws DomainError compute_atmospheric_optical_depth(
            Dict(:H2O => -1.0), R_50km, opacities
        )
        @test_throws DomainError compute_atmospheric_optical_depth(M_atm, 0.0, opacities)
        @test_throws DomainError compute_atmospheric_optical_depth(
            M_atm, -R_50km, opacities
        )
        @test_throws DomainError compute_atmospheric_optical_depth(
            M_atm, R_50km, Dict(:H2O => -0.01)
        )
    end

    # ---------------------------------------------------------------------
    # 4. Guillot (2010) Semi-Grey Radiative Equilibrium Profile
    # ---------------------------------------------------------------------
    @testset "Guillot (2010) Radiative Equilibrium Profile" begin
        T_int = 100.0
        T_irr = 300.0
        gamma = 0.10
        albedo = 0.20

        # Limit tau -> 0 (optically thin skin temperature)
        # Guillot (2010) Eq. 49 at tau = 0:
        # T^4(0) = (1/2) * T_int^4 + (1/2 + (√3/4) * gamma) * T_eqm^4
        T_surf_thin = compute_guillot_surface_temperature(
            0.0, T_int, T_irr; gamma=gamma, albedo=albedo
        )
        T_eqm_test = (0.25 * (1.0 - albedo) * T_irr^4)^0.25
        T4_thin_expected =
            0.5 * (T_int^4) + (0.5 + (sqrt(3.0) / 4.0) * gamma) * (T_eqm_test^4)
        @test isapprox(T_surf_thin, (T4_thin_expected)^0.25, rtol=1e-12)
        @test T_surf_thin > 0.0

        # Deep optical depth (greenhouse warming)
        tau_thick = 10.0
        T_surf_thick = compute_guillot_surface_temperature(
            tau_thick, T_int, T_irr; gamma=gamma, albedo=albedo
        )
        @test T_surf_thick > T_surf_thin

        # Monotonicity: higher optical depth increases surface temperature
        tau_dense = 50.0
        T_surf_dense = compute_guillot_surface_temperature(
            tau_dense, T_int, T_irr; gamma=gamma, albedo=albedo
        )
        @test T_surf_dense > T_surf_thick

        # Zero interior flux limit
        T_surf_noint = compute_guillot_surface_temperature(
            0.0, 0.0, T_irr; gamma=gamma, albedo=albedo
        )
        T4_noint_expected = (0.5 + (sqrt(3.0) / 4.0) * gamma) * (T_eqm_test^4)
        @test isapprox(T_surf_noint, (T4_noint_expected)^0.25, rtol=1e-12)

        # Explicit T_eqm argument consistency (verifies albedo is not double-counted)
        T_surf_with_eqm = compute_guillot_surface_temperature(
            0.0, T_int, T_irr; T_eqm=T_eqm_test, gamma=gamma, albedo=albedo
        )
        @test isapprox(T_surf_with_eqm, T_surf_thin, rtol=1e-12)

        # Albedo sensitivity: Higher albedo cools surface equilibrium temperature
        T_surf_alb10 = compute_guillot_surface_temperature(
            2.0, T_int, T_irr; gamma=gamma, albedo=0.10
        )
        T_surf_alb80 = compute_guillot_surface_temperature(
            2.0, T_int, T_irr; gamma=gamma, albedo=0.80
        )
        @test T_surf_alb80 < T_surf_alb10

        # Error guards
        @test_throws DomainError compute_guillot_surface_temperature(
            -1.0, T_int, T_irr; gamma=gamma, albedo=albedo
        )
        @test_throws DomainError compute_guillot_surface_temperature(
            0.0, -10.0, T_irr; gamma=gamma, albedo=albedo
        )
        @test_throws DomainError compute_guillot_surface_temperature(
            0.0, T_int, -10.0; gamma=gamma, albedo=albedo
        )
        @test_throws DomainError compute_guillot_surface_temperature(
            0.0, T_int, T_irr; gamma=0.0, albedo=albedo
        )
        @test_throws DomainError compute_guillot_surface_temperature(
            0.0, T_int, T_irr; gamma=gamma, albedo=1.5
        )
        @test_throws DomainError compute_guillot_surface_temperature(
            0.0, T_int, T_irr; gamma=gamma, albedo=-0.1
        )
    end

    # ---------------------------------------------------------------------
    # 5. Greenhouse-Attenuated Surface Heat Transfer Coefficient
    # ---------------------------------------------------------------------
    @testset "Greenhouse-Attenuated Radiation HTC" begin
        T_surf = 400.0
        T_amb = 200.0
        h_rad_bare = compute_radiation_htc(T_surf, T_amb; emissivity=0.9)

        # Asymptotic limit: Zero optical depth returns bare radiation HTC
        h_rad_0 = compute_effective_radiation_htc(T_surf, T_amb, 0.0; emissivity=0.9)
        @test isapprox(h_rad_0, h_rad_bare, rtol=1e-12)

        # Greenhouse blanket attenuation: h_eff = h_bare / (1 + 0.75 * tau)
        tau_test = 4.0
        h_rad_att = compute_effective_radiation_htc(T_surf, T_amb, tau_test; emissivity=0.9)
        h_rad_expected = h_rad_bare / (1.0 + 0.75 * tau_test)
        @test isapprox(h_rad_att, h_rad_expected, rtol=1e-12)
        @test h_rad_att < h_rad_bare
        @test h_rad_att > 0.0

        # Monotonicity: strictly decreasing with optical depth
        h_rad_thick = compute_effective_radiation_htc(T_surf, T_amb, 20.0; emissivity=0.9)
        @test h_rad_thick < h_rad_att

        # Module Physics submodule export consistency
        @test Erebus.Physics.compute_effective_radiation_htc ===
            Erebus.compute_effective_radiation_htc

        # Error guards
        @test_throws DomainError compute_effective_radiation_htc(
            T_surf, T_amb, -1.0; emissivity=0.9
        )
        @test_throws DomainError compute_effective_radiation_htc(
            -T_surf, T_amb, tau_test; emissivity=0.9
        )
        @test_throws DomainError compute_effective_radiation_htc(
            T_surf, 0.0, tau_test; emissivity=0.9
        )
        @test_throws DomainError compute_effective_radiation_htc(
            T_surf, T_amb, tau_test; emissivity=-0.1
        )
        @test_throws DomainError compute_effective_radiation_htc(
            T_surf, T_amb, tau_test; emissivity=1.5
        )
    end

    # ---------------------------------------------------------------------
    # 6. Hydrodynamic Boil-Off Loss Rate
    # ---------------------------------------------------------------------
    @testset "Hydrodynamic Boil-Off Loss Rate" begin
        tau_boil = 1.0e4 * Erebus.SEC_PER_YEAR
        M_env = 1.0e15 # kg
        M_target = 0.5e15 # kg

        # Excess envelope: rate = (M_env - M_target) / tau_boil
        rate = compute_boiloff_rate(M_env, M_target, tau_boil)
        @test isapprox(rate, (M_env - M_target) / tau_boil, rtol=1e-12)
        @test rate > 0.0

        # Under-filled or equal envelope: zero boil-off rate
        @test iszero(compute_boiloff_rate(M_target, M_target, tau_boil))
        @test iszero(compute_boiloff_rate(0.5 * M_target, M_target, tau_boil))

        # Linear scaling with excess mass
        rate_2x = compute_boiloff_rate(2.0 * M_env - M_target, M_target, tau_boil)
        @test isapprox(rate_2x, 2.0 * rate, rtol=1e-12)

        # Error guards
        @test_throws DomainError compute_boiloff_rate(-1.0, M_target, tau_boil)
        @test_throws DomainError compute_boiloff_rate(M_env, -1.0, tau_boil)
        @test_throws DomainError compute_boiloff_rate(M_env, M_target, 0.0)
        @test_throws DomainError compute_boiloff_rate(M_env, M_target, -tau_boil)
    end

    # ---------------------------------------------------------------------
    # 7. Zahnle-Kasting (1986) Crossover Mass and Drag Fractionation
    # ---------------------------------------------------------------------
    @testset "Zahnle-Kasting (1986) Hydrodynamic Crossover Escape" begin
        m_H2 = Erebus.MASS_H2_KG
        m_H2O = Erebus.MASS_H2O_KG
        m_CO2 = Erebus.MASS_CO2_KG
        T_exo = 400.0
        b_diff = 1.0e21 # m^-1 s^-1
        X_carrier = 0.90
        Phi_H2 = 1.0e18 # molecules / (m^2 s)

        # Crossover mass calculation
        m_c = compute_crossover_mass(m_H2, T_exo, Phi_H2, g_50km, X_carrier; b_diff=b_diff)
        m_c_expected =
            m_H2 + (Erebus.K_BOLTZMANN * T_exo * Phi_H2) / (b_diff * g_50km * X_carrier)
        @test isapprox(m_c, m_c_expected, rtol=1e-12)
        @test m_c > m_H2

        # Asymptotic limit: Zero carrier escape flux gives m_c = m_carrier
        m_c_zero = compute_crossover_mass(
            m_H2, T_exo, 0.0, g_50km, X_carrier; b_diff=b_diff
        )
        @test isapprox(m_c_zero, m_H2, rtol=1e-12)

        # Monotonicity with escape flux
        m_c_high = compute_crossover_mass(
            m_H2, T_exo, 2.0 * Phi_H2, g_50km, X_carrier; b_diff=b_diff
        )
        @test m_c_high > m_c

        # Drag fractionation factor x_j
        # Case A: Species is the carrier -> drag efficiency = 1.0
        @test isapprox(compute_crossover_drag_fraction(m_H2, m_c, m_H2), 1.0, rtol=1e-12)

        # Case B: Species lighter than carrier -> drag efficiency = 1.0
        @test isapprox(
            compute_crossover_drag_fraction(0.5 * m_H2, m_c, m_H2), 1.0, rtol=1e-12
        )

        # Case C: Species mass exactly at crossover mass -> drag efficiency = 0.0
        @test iszero(compute_crossover_drag_fraction(m_c, m_c, m_H2))

        # Case D: Species mass exceeds crossover mass -> drag efficiency = 0.0
        @test iszero(compute_crossover_drag_fraction(2.0 * m_c, m_c, m_H2))

        # Case E: Intermediate mass m_H2 < m_sp < m_c -> 0 < x_j < 1
        m_mid = 0.5 * (m_H2 + m_c)
        x_mid = compute_crossover_drag_fraction(m_mid, m_c, m_H2)
        @test isapprox(x_mid, 0.5, rtol=1e-12)

        # Error guards
        @test_throws DomainError compute_crossover_mass(
            0.0, T_exo, Phi_H2, g_50km, X_carrier
        )
        @test_throws DomainError compute_crossover_mass(
            m_H2, 0.0, Phi_H2, g_50km, X_carrier
        )
        @test_throws DomainError compute_crossover_mass(
            m_H2, T_exo, -1.0, g_50km, X_carrier
        )
        @test_throws DomainError compute_crossover_mass(m_H2, T_exo, Phi_H2, 0.0, X_carrier)
        @test_throws DomainError compute_crossover_mass(m_H2, T_exo, Phi_H2, g_50km, 0.0)
        @test_throws DomainError compute_crossover_mass(
            m_H2, T_exo, Phi_H2, g_50km, X_carrier; b_diff=0.0
        )
        @test_throws DomainError compute_crossover_drag_fraction(0.0, m_c, m_H2)
        @test_throws DomainError compute_crossover_drag_fraction(m_H2O, m_H2 - 1.0, m_H2)
        @test_throws DomainError compute_crossover_drag_fraction(m_H2O, m_c, 0.0)
    end

    # ---------------------------------------------------------------------
    # 8. Coupled Atmospheric Time Step & Exact Mass Conservation
    # ---------------------------------------------------------------------
    @testset "Coupled Atmosphere Time Step Conservation & Dynamics" begin
        dt_s = 100.0 * Erebus.SEC_PER_YEAR
        opacities = Dict(:H2O => 1.0e-2, :CO2 => 1.0e-3, :H2 => 1.0e-5)

        M_atm_init = Dict(:H2O => 1.0e14, :CO2 => 5.0e13, :H2 => 1.0e12)
        M_esc_init = Dict(:H2O => 0.0, :CO2 => 0.0, :H2 => 0.0)

        # Influxes from interior venting
        vent_rates = Dict(:H2O => 1.0e5, :CO2 => 2.0e4, :H2 => 1.0e3) # kg/s

        atm_state = AtmosphereState(
            deepcopy(M_atm_init),
            deepcopy(M_esc_init),
            0.0, # P_surf
            0.0, # T_surf_eq
            0.0, # tau_LW
            0.0, # M_env_bound
            0.0, # F_net_rad
            0.0, # h_rad_eff
        )

        cfg_atm = AtmosphereConfig(;
            active=true, mode=:guillot, opacities=opacities, crossover_active=true
        )

        # Advance one step
        evolve_coupled_atmosphere_step!(
            atm_state,
            vent_rates,
            dt_s,
            M_50km,
            R_50km,
            T_disk,
            cfg_atm;
            rho_disk=0.0,
            c_s=c_s,
            M_star=M_sun,
            a_orb=a_test,
        )

        # Invariant 1: Exact global mass conservation to machine precision
        for sp in keys(M_atm_init)
            M_in = M_atm_init[sp] + vent_rates[sp] * dt_s
            M_out = atm_state.M_atm[sp] + atm_state.M_escaped[sp]
            @test isapprox(M_out, M_in, rtol=1e-12)
        end

        # Discrimination test: lighter species (H2O, 18 g/mol) undergoes higher escaped fraction than heavier (CO2, 44 g/mol)
        frac_esc_h2o =
            atm_state.M_escaped[:H2O] / (M_atm_init[:H2O] + vent_rates[:H2O] * dt_s)
        frac_esc_co2 =
            atm_state.M_escaped[:CO2] / (M_atm_init[:CO2] + vent_rates[:CO2] * dt_s)
        @test frac_esc_h2o > frac_esc_co2
        @test frac_esc_h2o > 0.0

        # Invariant 2: Surface pressure matches analytical gravity relation
        M_tot = sum(values(atm_state.M_atm))
        P_expected = compute_surface_atmospheric_pressure(M_tot, M_50km, R_50km)
        @test isapprox(atm_state.P_surf, P_expected, rtol=1e-12)
        @test atm_state.P_surf > 0.0

        # Invariant 3: Optical depth is strictly positive
        @test atm_state.tau_LW > 0.0

        # Invariant 4: Surface temperature is strictly positive and elevated by greenhouse above skin temperature
        T_skin = compute_guillot_surface_temperature(
            0.0,
            100.0,
            T_disk;
            T_eqm=T_disk,
            gamma=cfg_atm.gamma_guillot,
            albedo=cfg_atm.albedo,
        )
        @test atm_state.T_surf_eq > 0.0
        @test atm_state.T_surf_eq > T_skin

        # Invariant 5: Effective radiation HTC is attenuated
        @test atm_state.h_rad_eff > 0.0
        h_bare = compute_radiation_htc(atm_state.T_surf_eq, T_disk)
        @test atm_state.h_rad_eff < h_bare

        # Invariant 6: Net radiative flux is computed and finite with distinct T_int
        @test isfinite(atm_state.F_net_rad)
        atm_tint = AtmosphereState(
            Dict(:H2O => 1.0e14), Dict(:H2O => 0.0), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        )
        evolve_coupled_atmosphere_step!(
            atm_tint, Dict(:H2O => 0.0), dt_s, M_50km, R_50km, T_disk, cfg_atm; T_int=300.0
        )
        F_net_expected = atm_tint.h_rad_eff * (300.0 - T_disk)
        @test isapprox(atm_tint.F_net_rad, F_net_expected, rtol=1e-12)
        @test atm_tint.F_net_rad > 0.0

        # Invariant 7: Disk gas envelope capture and growth for a massive embryo
        R_embryo = 1.5e6
        M_embryo = 1.0e23
        atm_disk = AtmosphereState(
            Dict(:H2 => 0.0, :H2O => 1.0e14),
            Dict(:H2 => 0.0, :H2O => 0.0),
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        )
        rho_disk_val = 1.0e-9
        evolve_coupled_atmosphere_step!(
            atm_disk,
            Dict(:H2 => 0.0, :H2O => 0.0),
            dt_s,
            M_embryo,
            R_embryo,
            T_disk,
            cfg_atm;
            rho_disk=rho_disk_val,
            c_s=c_s,
            M_star=M_sun,
            a_orb=a_test,
        )
        @test atm_disk.M_env_bound > 0.0
        @test atm_disk.M_atm[:H2] > 0.0
        # Bound envelope mass cannot exceed available atmospheric H2
        @test atm_disk.M_env_bound <= atm_disk.M_atm[:H2] * (1.0 + 1e-12)

        # Invariant 8: Envelope boil-off timescale sensitivity during disk dispersal
        atm_boil1 = AtmosphereState(
            Dict(:H2 => 1.0e15), Dict(:H2 => 0.0), 0.0, 0.0, 0.0, 1.0e15, 0.0, 0.0
        )
        cfg_boil_fast = AtmosphereConfig(active=true, tau_boil=100.0 * Erebus.SEC_PER_YEAR)
        evolve_coupled_atmosphere_step!(
            atm_boil1,
            Dict(:H2 => 0.0),
            10.0 * Erebus.SEC_PER_YEAR,
            M_embryo,
            R_embryo,
            T_disk,
            cfg_boil_fast;
            rho_disk=0.0,
            escape_active=false,
        )
        atm_boil2 = AtmosphereState(
            Dict(:H2 => 1.0e15), Dict(:H2 => 0.0), 0.0, 0.0, 0.0, 1.0e15, 0.0, 0.0
        )
        cfg_boil_slow = AtmosphereConfig(active=true, tau_boil=1000.0 * Erebus.SEC_PER_YEAR)
        evolve_coupled_atmosphere_step!(
            atm_boil2,
            Dict(:H2 => 0.0),
            10.0 * Erebus.SEC_PER_YEAR,
            M_embryo,
            R_embryo,
            T_disk,
            cfg_boil_slow;
            rho_disk=0.0,
            escape_active=false,
        )
        # Shorter tau_boil unbinds and boils off envelope mass faster
        @test atm_boil1.M_env_bound < atm_boil2.M_env_bound
        @test isapprox(atm_boil1.M_atm[:H2] + atm_boil1.M_escaped[:H2], 1.0e15; rtol=1e-12)
        @test isapprox(atm_boil2.M_atm[:H2] + atm_boil2.M_escaped[:H2], 1.0e15; rtol=1e-12)

        # Invariant 9: Multispecies escape when H2 is absent (resolving carrier deadlock)
        atm_noh2 = AtmosphereState(
            Dict(:H2O => 1.0e14, :CO2 => 1.0e14),
            Dict(:H2O => 0.0, :CO2 => 0.0),
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        )
        evolve_coupled_atmosphere_step!(
            atm_noh2, Dict(:H2O => 0.0, :CO2 => 0.0), dt_s, M_50km, R_50km, T_disk, cfg_atm
        )
        @test isapprox(atm_noh2.M_atm[:H2O] + atm_noh2.M_escaped[:H2O], 1.0e14, rtol=1e-12)
        @test isapprox(atm_noh2.M_atm[:CO2] + atm_noh2.M_escaped[:CO2], 1.0e14, rtol=1e-12)
        @test atm_noh2.M_escaped[:H2O] > 0.0
        @test atm_noh2.M_escaped[:CO2] > 0.0

        # Guard against NaN dt_s and invalid T_int
        atm_nan = AtmosphereState(
            Dict(:H2O => 1.0e14), Dict(:H2O => 0.0), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        )
        @test_throws DomainError evolve_coupled_atmosphere_step!(
            atm_nan, Dict(:H2O => 0.0), NaN, M_50km, R_50km, T_disk, cfg_atm
        )
        @test_throws DomainError evolve_coupled_atmosphere_step!(
            atm_nan, Dict(:H2O => 0.0), dt_s, M_50km, R_50km, T_disk, cfg_atm; T_int=-10.0
        )
        @test_throws DomainError evolve_coupled_atmosphere_step!(
            atm_nan, Dict(:H2O => 0.0), dt_s, M_50km, R_50km, T_disk, cfg_atm; T_int=NaN
        )

        # Invariant 10: Atmospheric modes and skin temperature floor
        cfg_iso = AtmosphereConfig(mode=:isothermal, T_skin_floor=60.0)
        evolve_coupled_atmosphere_step!(
            atm_noh2, Dict(:H2O => 0.0, :CO2 => 0.0), dt_s, M_50km, R_50km, T_disk, cfg_iso;
        )
        @test isapprox(atm_noh2.T_surf_eq, T_disk; rtol=1e-12)

        cfg_grey = AtmosphereConfig(mode=:grey)
        evolve_coupled_atmosphere_step!(
            atm_noh2,
            Dict(:H2O => 0.0, :CO2 => 0.0),
            dt_s,
            M_50km,
            R_50km,
            T_disk,
            cfg_grey;
        )
        T_grey_expected = T_disk * (1.0 + 0.75 * atm_noh2.tau_LW)^0.25
        @test isapprox(atm_noh2.T_surf_eq, T_grey_expected; rtol=1e-12)

        cfg_floor = AtmosphereConfig(mode=:isothermal, T_skin_floor=250.0)
        evolve_coupled_atmosphere_step!(
            atm_noh2,
            Dict(:H2O => 0.0, :CO2 => 0.0),
            dt_s,
            M_50km,
            R_50km,
            T_disk,
            cfg_floor;
        )
        @test isapprox(atm_noh2.T_surf_eq, 250.0; rtol=1e-12)

        # Error guards on non-positive / non-finite inputs
        @test_throws DomainError evolve_coupled_atmosphere_step!(
            atm_state, Dict(:H2O => 0.0), dt_s, 0.0, R_50km, T_disk, cfg_atm
        )
        @test_throws DomainError evolve_coupled_atmosphere_step!(
            atm_state, Dict(:H2O => 0.0), dt_s, M_50km, 0.0, T_disk, cfg_atm
        )
        @test_throws DomainError evolve_coupled_atmosphere_step!(
            atm_state, Dict(:H2O => 0.0), dt_s, M_50km, R_50km, 0.0, cfg_atm
        )
    end

    # ---------------------------------------------------------------------
    # 9. Configuration & TOML Schema Validation
    # ---------------------------------------------------------------------
    @testset "AtmosphereConfig Schema & Validation" begin
        # Default construction
        cfg_def = AtmosphereConfig()
        @test cfg_def.active == false
        @test cfg_def.mode == :guillot
        @test isapprox(cfg_def.albedo, 0.20; rtol=1e-12)
        @test isapprox(cfg_def.gamma_guillot, 0.10; rtol=1e-12)
        @test haskey(cfg_def.opacities, :H2O)

        # SimulationConfig inclusion
        sim_cfg = SimulationConfig(; atmosphere=AtmosphereConfig(; active=true))
        @test sim_cfg.atmosphere.active == true

        # Validation success
        @test validate_config(sim_cfg) === nothing

        # Validation failures
        # 1. Negative opacity
        bad_op = SimulationConfig(;
            atmosphere=AtmosphereConfig(; active=true, opacities=Dict(:H2O => -0.01))
        )
        @test_throws ArgumentError validate_config(bad_op)

        # 2. Invalid mode
        bad_mode = SimulationConfig(;
            atmosphere=AtmosphereConfig(; active=true, mode=:invalid_mode)
        )
        @test_throws ArgumentError validate_config(bad_mode)

        # 3. Invalid albedo >= 1
        bad_alb = SimulationConfig(; atmosphere=AtmosphereConfig(; active=true, albedo=1.5))
        @test_throws ArgumentError validate_config(bad_alb)

        # 4. Non-positive gamma
        bad_gam = SimulationConfig(;
            atmosphere=AtmosphereConfig(; active=true, gamma_guillot=0.0)
        )
        @test_throws ArgumentError validate_config(bad_gam)

        # 5. Non-positive boil-off timescale
        bad_boil = SimulationConfig(;
            atmosphere=AtmosphereConfig(; active=true, tau_boil=0.0)
        )
        @test_throws ArgumentError validate_config(bad_boil)

        # 6. Non-positive kappa_ir_default
        bad_kir = SimulationConfig(;
            atmosphere=AtmosphereConfig(; active=true, kappa_ir_default=0.0)
        )
        @test_throws ArgumentError validate_config(bad_kir)

        # 7. Non-positive kappa_vis_default
        bad_kvis = SimulationConfig(;
            atmosphere=AtmosphereConfig(; active=true, kappa_vis_default=-0.01)
        )
        @test_throws ArgumentError validate_config(bad_kvis)

        # 8. Invalid f_rec (must be in (0, 1])
        bad_frec_hi = SimulationConfig(;
            atmosphere=AtmosphereConfig(; active=true, f_rec=1.5)
        )
        @test_throws ArgumentError validate_config(bad_frec_hi)
        bad_frec_lo = SimulationConfig(;
            atmosphere=AtmosphereConfig(; active=true, f_rec=-0.1)
        )
        @test_throws ArgumentError validate_config(bad_frec_lo)

        # TOML Round-Trip Serialization
        toml_str = serialize_config(sim_cfg)
        @test occursin("[atmosphere]", toml_str)
        @test occursin("active = true", toml_str)

        cfg_parsed = parse_config_string(toml_str)
        @test cfg_parsed.atmosphere.active == true
        @test cfg_parsed.atmosphere.mode == :guillot
        @test isapprox(cfg_parsed.atmosphere.albedo, 0.20, rtol=1e-12)
    end

    # ---------------------------------------------------------------------
    # 10. Simulation Loop Integration with Coupled Atmosphere
    # ---------------------------------------------------------------------
    @testset "Simulation Loop Integration with Coupled Atmosphere" begin
        output_dir = mktempdir()
        try
            quick_toml = joinpath(@__DIR__, "..", "configs", "test_quick.toml")
            cfg = load_config(quick_toml)

            cfg_atm = SimulationConfig(
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
                    active=false, species_list=[:H2O, :H2, :CO2, :N2, :H2S]
                ),
                atmosphere=AtmosphereConfig(active=true, mode=:guillot),
            )

            Erebus.simulation_loop(cfg_atm; output_path=output_dir)

            files = readdir(output_dir)
            @test "output_00000.jld2" in files
            @test "output_00002.jld2" in files

            data2 = load_state(joinpath(output_dir, "output_00002.jld2"))
            @test data2["timestep"] == 2
            @test haskey(data2, "M_planet_val")
            @test data2["M_planet_val"] > 0.0
            @test haskey(data2, "atm_P_surf")
            @test haskey(data2, "atm_T_surf_eq")
            @test haskey(data2, "atm_tau_LW")
            @test haskey(data2, "atm_h_rad_eff")
            @test haskey(data2, "atm_F_net_rad")
            @test haskey(data2, "atm_M_atm")
            @test haskey(data2, "atm_M_escaped")
            @test data2["atm_P_surf"] >= 0.0
            @test data2["atm_T_surf_eq"] > 0.0
            @test data2["atm_tau_LW"] >= 0.0
            @test data2["atm_h_rad_eff"] > 0.0
            @test isfinite(data2["atm_F_net_rad"])

            # Test retention active with venting inactive (scoping crash regression guard)
            cfg_ret_guard = SimulationConfig(
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
                    n_steps=1,
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
                retention=RetentionConfig(active=true, venting_drainage_active=true),
                escape=EscapeConfig(
                    active=false, species_list=[:H2O, :H2, :CO2, :N2, :H2S]
                ),
                atmosphere=AtmosphereConfig(active=true, mode=:guillot),
            )
            Erebus.simulation_loop(cfg_ret_guard; output_path=output_dir)

            # Test restart/resume preserving atmospheric state
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
                output=OutputConfig(output_dir=output_dir, savematstep=3),
                disk=cfg.disk,
                melting=cfg.melting,
                venting=VentingConfig(active=false),
                escape=EscapeConfig(
                    active=false, species_list=[:H2O, :H2, :CO2, :N2, :H2S]
                ),
                atmosphere=AtmosphereConfig(active=true, mode=:guillot),
            )

            Erebus.simulation_loop(
                cfg_resume;
                output_path=output_dir,
                restart_from=joinpath(output_dir, "output_00002.jld2"),
            )

            files_after = readdir(output_dir)
            @test "output_00003.jld2" in files_after
            data3 = load_state(joinpath(output_dir, "output_00003.jld2"))
            @test data3["timestep"] == 3
            @test haskey(data3, "atm_P_surf")
            @test haskey(data3, "atm_T_surf_eq")
            @test data3["atm_T_surf_eq"] > 0.0
            @test haskey(data3, "atm_M_atm")
            @test haskey(data3, "atm_M_escaped")
            @test haskey(data3, "atm_M_env_bound")
            @test data3["atm_M_atm"] == data2["atm_M_atm"]
        finally
            rm(output_dir; recursive=true, force=true)
        end

        # Coupled Atmosphere Volatile Transfer: Darcy Porosity Venting + Retention Drainage
        output_dir_vent = mktempdir()
        try
            quick_toml = joinpath(@__DIR__, "..", "configs", "test_quick.toml")
            cfg = load_config(quick_toml)

            cfg_vent = SimulationConfig(
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
                reaction=ReactionConfig(active=false),
                materials=cfg.materials,
                output=OutputConfig(output_dir=output_dir_vent, savematstep=2),
                disk=cfg.disk,
                melting=cfg.melting,
                venting=VentingConfig(
                    active=true, mode=:darcy_sink, k_vent=1.0e-11, conductance_factor=1.0
                ),
                volatiles=VolatilesConfig(
                    active=true,
                    carbon_active=true,
                    sulfur_active=true,
                    initial_water_wtpct=1.0,
                    initial_carbon_ppm=500.0,
                    initial_nitrogen_ppm=50.0,
                    initial_sulfur_ppm=1000.0,
                ),
                retention=RetentionConfig(
                    active=true,
                    venting_drainage_active=true,
                    h2o_retention_ppm=50.0,
                    carbon_retention_ppm=50.0,
                    nitrogen_retention_ppm=5.0,
                    sulfur_retention_ppm=100.0,
                ),
                escape=EscapeConfig(
                    active=false, species_list=[:H2O, :H2, :CO2, :N2, :H2S]
                ),
                atmosphere=AtmosphereConfig(active=true, mode=:guillot),
            )

            Erebus.simulation_loop(cfg_vent; output_path=output_dir_vent)

            data2 = load_state(joinpath(output_dir_vent, "output_00002.jld2"))
            @test data2["timestep"] == 2
            @test haskey(data2, "M_vent_total")
            @test haskey(data2, "M_vent_H2O_total")
            @test data2["M_vent_total"] > 0.0
            @test data2["M_vent_H2O_total"] > 0.0
            # Both pore water and mineral water contribute additively to atmospheric H2O budget
            @test isapprox(
                data2["atm_M_atm"][:H2O],
                data2["M_vent_total"] + data2["M_vent_H2O_total"],
                rtol=1e-10,
            )
            # Mobile mineral volatiles enter respective atmospheric species
            @test data2["atm_M_atm"][:CO2] > 0.0
            @test data2["atm_M_atm"][:N2] > 0.0
            @test data2["atm_M_atm"][:H2S] > 0.0
            @test isapprox(
                data2["atm_M_atm"][:CO2],
                data2["M_vent_C_total"] * (44.0095 / 12.011),
                rtol=1e-10,
            )
            @test isapprox(data2["atm_M_atm"][:N2], data2["M_vent_N_total"], rtol=1e-10)
            @test isapprox(
                data2["atm_M_atm"][:H2S],
                data2["M_vent_S_total"] * (34.08 / 32.06),
                rtol=1e-10,
            )
        finally
            rm(output_dir_vent; recursive=true, force=true)
        end
    end

    # ---------------------------------------------------------------------
    # 11. Multispecies Escape: Binary Reduction vs IsoFATE
    # ---------------------------------------------------------------------
    @testset "Multispecies Escape: Binary Reduction vs IsoFATE" begin
        function isofate_binary(phi, x1, x2, m1, m2, T, g0, b12)
            kT = Erebus.K_BOLTZMANN * T
            H1 = kT / (m1 * g0)
            H2 = kT / (m2 * g0)
            phi_c = b12 * x1 * (m2 - m1) / H1
            if phi < phi_c
                return phi / m1, 0.0, phi_c
            end
            mbar = m1 * x1 + m2 * x2
            Phi1 = (x1 * phi + x1 * x2 * (m2 - m1) * b12 / H2) / mbar
            Phi2 = (x2 * phi + x1 * x2 * (m1 - m2) * b12 / H1) / mbar
            return Phi1, Phi2, phi_c
        end

        rng = MersenneTwister(101)
        for _ in 1:200
            m1 = Erebus.ATOMIC_MASS_UNIT * (1.0 + 3.0 * rand(rng))
            m2 = m1 * (1.5 + 38.5 * rand(rng))
            x2 = 0.01 + 0.89 * rand(rng)
            x1 = 1.0 - x2
            T = 200.0 + 1800.0 * rand(rng)
            g0 = 1.0 + 29.0 * rand(rng)
            b12 = (0.5 + 4.5 * rand(rng)) * 1e21 * (T / 1000.0)^0.75
            b = [Inf b12; b12 Inf]
            kT = Erebus.K_BOLTZMANN * T
            phi_c = b12 * x1 * (m2 - m1) * m1 * g0 / kT

            for frac in (0.2, 0.999999, 1.000001, 1.7, 30.0)
                phi = frac * phi_c
                Phi = solve_multispecies_escape_closure(phi, [x1, x2], [m1, m2], T, g0, b)
                P1, P2, _ = isofate_binary(phi, x1, x2, m1, m2, T, g0, b12)
                err = max(abs(Phi[1] - P1), abs(Phi[2] - P2)) / max(P1, 1e-300)
                @test err < 1e-11
                @test isapprox(m1 * Phi[1] + m2 * Phi[2], phi, rtol=1e-12)
            end
        end
    end

    # ---------------------------------------------------------------------
    # 12. Multispecies Escape: Ternary Reduction vs Gu & Chen (2023)
    # ---------------------------------------------------------------------
    @testset "Multispecies Escape: Ternary Reduction vs Gu & Chen (2023)" begin
        m = [1.008, 4.0026, 2.0141] .* Erebus.ATOMIC_MASS_UNIT
        T = 1000.0
        g0 = 10.0
        kT = Erebus.K_BOLTZMANN * T
        b12 = 1.04e20 * (T^0.732)
        b13 = 7.183e19 * (T^0.728)
        b23 = 5.087e19 * (T^0.728)
        b = [Inf b12 b13; b12 Inf b23; b13 b23 Inf]
        X3 = 1e-13
        X2 = 0.15
        X1 = 1.0 - X2 - X3
        X = [X1, X2, X3]
        a2, a3 = b13 / b12, b13 / b23
        Phid(bij, mi) = bij * (mi - m[1]) * g0 / kT
        phi_DL_He = Phid(b12, m[2])
        phi_DL_D = Phid(b13, m[3])
        phi_crit_He = m[1] * X1 * phi_DL_He
        phi_crit_D = m[1] * phi_DL_D / (1.0 + a2 * X2 / X1)

        # Supercritical: He escaping
        for frac in (1.5, 5.0, 50.0)
            phi = frac * phi_crit_He
            Phi, C, act = solve_multispecies_escape_closure(
                phi, X, m, T, g0, b; return_diag=true
            )
            f2, f3 = X2 / X1, X3 / X1
            Phi3_ref =
                f3 * (Phi[1] + a3 * Phi[2] + a2 * phi_DL_He * X2 - phi_DL_D) /
                (1.0 + a3 * f2)
            @test isapprox(Phi[3], Phi3_ref, rtol=1e-10)
            @test act == Set([1, 2, 3])
        end

        # Subcritical: He retained, D escaping
        for frac in (1.5, 3.0)
            phi = frac * phi_crit_D
            if phi >= phi_crit_He
                continue
            end
            Phi, C, act = solve_multispecies_escape_closure(
                phi, X, m, T, g0, b; return_diag=true
            )
            @test act == Set([1, 3])
            @test Phi[2] == 0.0
            Phi3_ref = X3 * (Phi[1] * (1.0 + a2 * X2 / X1) - phi_DL_D) / (X1 + a3 * X2)
            @test isapprox(Phi[3], Phi3_ref, rtol=1e-10)
        end
    end

    # ---------------------------------------------------------------------
    # 13. Multispecies Escape: Chassefiere (1996) Binary Partition
    # ---------------------------------------------------------------------
    @testset "Multispecies Escape: Chassefiere (1996) Binary Partition" begin
        m1 = 1.008 * Erebus.ATOMIC_MASS_UNIT
        m2 = 15.999 * Erebus.ATOMIC_MASS_UNIT
        x1, x2 = 0.667, 0.333
        T = 1000.0
        g0 = 9.8
        kT = Erebus.K_BOLTZMANN * T
        b12 = 4.8e21
        b = [Inf b12; b12 Inf]
        phi_star = compute_escape_activation_threshold([x1, x2], [m1, m2], T, g0, b)

        for frac in (1.2, 2.0, 10.0)
            phi = frac * phi_star
            Phi = solve_multispecies_escape_closure(phi, [x1, x2], [m1, m2], T, g0, b)
            diff_drift = Phi[1] / x1 - Phi[2] / x2
            expected_diff = b12 * (m2 - m1) * g0 / kT
            @test isapprox(diff_drift, expected_diff, rtol=1e-11)
            @test isapprox(m1 * Phi[1] + m2 * Phi[2], phi, rtol=1e-12)
        end
    end

    # ---------------------------------------------------------------------
    # 14. Multispecies Escape: Three-Species Limiting Flux (Z90 Eq. 42)
    # ---------------------------------------------------------------------
    @testset "Multispecies Escape: Three-Species Limiting Flux" begin
        m = [2.01588, 44.0095, 28.014] .* Erebus.ATOMIC_MASS_UNIT
        T, g0 = 400.0, 3.73
        kT = Erebus.K_BOLTZMANN * T
        b12 = 2.3e19 * (T^0.75)
        b13 = 2.65e19 * (T^0.75)
        b23 = 1.0e19 * (T^0.75)
        b = [Inf b12 b13; b12 Inf b23; b13 b23 Inf]
        f2, f3 = 1.0, 0.5
        X1 = 1.0 / (1.0 + f2 + f3)
        X = [X1, f2 * X1, f3 * X1]

        phi12 = b12 * (m[2] - m[1]) * g0 / kT / (1.0 + f2 + (b12 / b13) * f3)
        phi13 = b13 * (m[3] - m[1]) * g0 / kT / (1.0 + f3 + (b13 / b12) * f2)
        first_thresh = min(phi12, phi13)
        phi_thresh = m[1] * first_thresh

        eps_val = 1e-6
        Phi_lo = solve_multispecies_escape_closure(
            phi_thresh * (1.0 - eps_val), X, m, T, g0, b
        )
        Phi_hi = solve_multispecies_escape_closure(
            phi_thresh * (1.0 + eps_val), X, m, T, g0, b
        )

        idx = phi12 < phi13 ? 2 : 3
        other_idx = phi12 < phi13 ? 3 : 2
        @test iszero(Phi_lo[idx])
        @test Phi_hi[idx] > 0.0
        @test iszero(Phi_lo[other_idx])
    end

    # ---------------------------------------------------------------------
    # 15. Binary Diffusion Parameters and Scaling Rules
    # ---------------------------------------------------------------------
    @testset "Binary Diffusion Matrix and Scaling Properties" begin
        # Direct lookup of tabulated pairs at 1000 K
        b_h_he = get_binary_diffusion_parameter(:H, :He, 1000.0)
        @test isapprox(b_h_he, 1.6e22, rtol=1e-12)

        # Symmetry: b_ij == b_ji
        @test isapprox(
            get_binary_diffusion_parameter(:He, :H, 1000.0),
            get_binary_diffusion_parameter(:H, :He, 1000.0),
            rtol=1e-12,
        )

        # Temperature scaling exponent: b ∝ T^0.75
        b_2000 = get_binary_diffusion_parameter(:H, :He, 2000.0)
        @test isapprox(b_2000, b_h_he * (2000.0 / 1000.0)^0.75, rtol=1e-12)

        # Self-diffusion is Infinite on diagonal
        @test isinf(get_binary_diffusion_parameter(:H, :H, 1000.0))

        # Matrix assembly
        species_test = [:H, :He, :O]
        b_mat = assemble_binary_diffusion_matrix(species_test, 1000.0)
        @test size(b_mat) == (3, 3)
        @test isinf(b_mat[1, 1])
        @test isinf(b_mat[2, 2])
        @test isinf(b_mat[3, 3])
        @test isapprox(b_mat[1, 2], b_mat[2, 1], rtol=1e-12)
        @test isapprox(b_mat[1, 3], b_mat[3, 1], rtol=1e-12)

        # Error guards
        @test_throws DomainError get_binary_diffusion_parameter(:H, :He, 0.0)
        @test_throws DomainError get_binary_diffusion_parameter(:H, :He, -100.0)
        @test_throws DomainError get_binary_diffusion_parameter(:H, :He, NaN)
        @test_throws DomainError solve_multispecies_escape_closure(
            -1.0, [0.5, 0.5], [1.0, 2.0], 1000.0, 10.0, b_mat[1:2, 1:2]
        )
        @test_throws DomainError solve_multispecies_escape_closure(
            1.0, [0.5, 0.5], [1.0, 2.0], 0.0, 10.0, b_mat[1:2, 1:2]
        )
        @test_throws DomainError solve_multispecies_escape_closure(
            1.0, [0.5, 0.5], [1.0, 2.0], 1000.0, 0.0, b_mat[1:2, 1:2]
        )
        @test_throws DimensionMismatch solve_multispecies_escape_closure(
            1.0, [0.5, 0.5], [1.0], 1000.0, 10.0, b_mat[1:2, 1:2]
        )
    end

    # ---------------------------------------------------------------------
    # 16. Standard Atomic Weights and Noble Gas Species Consistency
    # ---------------------------------------------------------------------
    @testset "Standard Atomic Weights and Noble Gas Species Consistency" begin
        # Single source of truth verification
        @test SPECIES_AMU_ESCAPE === SPECIES_AMU
        @test haskey(SPECIES_AMU, :Ar)
        @test haskey(SPECIES_AMU, :Ne)

        # IUPAC standard atomic weights (Meija et al. 2016)
        @test isapprox(SPECIES_AMU[:Ar], 39.948, atol=1e-3)
        @test isapprox(SPECIES_AMU[:Ne], 20.1797, atol=1e-4)

        # Molecular mass lookup in physics module
        @test isapprox(
            get_species_molecular_mass(:Ar), 39.948 * ATOMIC_MASS_UNIT, rtol=1e-12
        )
        @test isapprox(
            get_species_molecular_mass(:ar), 39.948 * ATOMIC_MASS_UNIT, rtol=1e-12
        )
        @test isapprox(
            get_species_molecular_mass(:Ne), 20.1797 * ATOMIC_MASS_UNIT, rtol=1e-12
        )
        @test isapprox(
            get_species_molecular_mass(:ne), 20.1797 * ATOMIC_MASS_UNIT, rtol=1e-12
        )

        # Binary diffusion matrix with noble gases
        sp_list = [:H, :He, :Ne, :Ar]
        b_mat = assemble_binary_diffusion_matrix(sp_list, 1000.0)
        @test size(b_mat) == (4, 4)
        @test isapprox(b_mat[1, 4], 6.5e21, rtol=1e-12) # H-Ar
        @test isapprox(b_mat[2, 4], 4.4906e21, rtol=1e-12) # He-Ar
        @test isapprox(b_mat[3, 4], 1.6161e21, rtol=1e-12) # Ne-Ar

        # Evolve step with arbitrary volatile mixture containing Ar and Ne
        atm_state = AtmosphereState()
        atm_state.M_atm[:H] = 1.0e15
        atm_state.M_atm[:Ar] = 1.0e14
        atm_state.M_atm[:Ne] = 1.0e14
        cfg = AtmosphereConfig(crossover_active=true)
        evolve_coupled_atmosphere_step!(
            atm_state, Dict{Symbol,Float64}(), 3600.0, 1.0e21, 5.0e5, 300.0, cfg
        )
        # Verify mass conservation holds with noble gases
        m_tot_final = sum(values(atm_state.M_atm)) + sum(values(atm_state.M_escaped))
        @test isapprox(m_tot_final, 1.2e15, rtol=1e-12)
    end
end

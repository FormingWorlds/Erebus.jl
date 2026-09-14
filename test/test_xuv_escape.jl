using Test
using Erebus
using StaticArrays
using LinearAlgebra

include("test_helpers.jl")

@testset "Hydrodynamic Boil-Off & XUV Crossover Escape Physics" begin
    # Reference planet parameters
    M_sun = 1.98847e30 # kg
    AU = 1.495978707e11 # m
    d_au = 1.0 # 1 AU
    a_orb = d_au * AU
    R_p = 50_000.0 # 50 km
    M_p = (4.0 / 3.0) * π * (R_p^3) * 3000.0 # ~1.57e18 kg

    # -------------------------------------------------------------------------
    # 1. Stellar XUV Flux Evolution & Distance Scaling
    # -------------------------------------------------------------------------
    @testset "Stellar XUV Flux Evolution" begin
        t_sat = 1.0e8 # 100 Myr
        F_sat = 1.361 # W/m^2 at 1 AU

        # Early saturated regime (t <= t_sat)
        F_early = compute_stellar_xuv_flux(1.0e6, d_au; F_xuv_1au_sat=F_sat, t_sat_yr=t_sat)
        @test isapprox(F_early, F_sat; rtol=1e-12)

        F_at_sat = compute_stellar_xuv_flux(
            t_sat, d_au; F_xuv_1au_sat=F_sat, t_sat_yr=t_sat
        )
        @test isapprox(F_at_sat, F_sat; rtol=1e-12)

        # Decaying regime (t > t_sat, beta = 1.23)
        t_late = 1.0e9 # 1 Gyr (10 * t_sat)
        F_late = compute_stellar_xuv_flux(
            t_late, d_au; F_xuv_1au_sat=F_sat, t_sat_yr=t_sat, beta=1.23
        )
        F_expected = F_sat * ((t_late / t_sat)^(-1.23))
        @test isapprox(F_late, F_expected; rtol=1e-12)
        @test F_late < F_sat

        # Inverse square law with distance
        F_2au = compute_stellar_xuv_flux(1.0e6, 2.0; F_xuv_1au_sat=F_sat, t_sat_yr=t_sat)
        @test isapprox(F_2au, F_sat / 4.0; rtol=1e-12)

        # Domain errors
        @test_throws DomainError compute_stellar_xuv_flux(-100.0, d_au)
        @test_throws DomainError compute_stellar_xuv_flux(1.0e6, 0.0)
        @test_throws DomainError compute_stellar_xuv_flux(1.0e6, -1.0)
    end

    # -------------------------------------------------------------------------
    # 2. Roche Lobe / Tidal Reduction Factor
    # -------------------------------------------------------------------------
    @testset "Roche Lobe Tidal Correction" begin
        # For small planetesimal at 1 AU, R_Hill >> R_p -> K_tide ≈ 1.0
        K_tide = compute_roche_lobe_correction(M_p, M_sun, a_orb, R_p)
        @test isapprox(K_tide, 1.0; atol=1e-2)
        @test 0.0 < K_tide <= 1.0

        # Close-in orbit increases tidal gravity, reducing K_tide
        a_close = 0.02 * AU
        K_tide_close = compute_roche_lobe_correction(M_p, M_sun, a_close, R_p)
        @test K_tide_close < K_tide

        # Domain error guards
        @test_throws DomainError compute_roche_lobe_correction(0.0, M_sun, a_orb, R_p)
        @test_throws DomainError compute_roche_lobe_correction(M_p, 0.0, a_orb, R_p)
        @test_throws DomainError compute_roche_lobe_correction(M_p, M_sun, 0.0, R_p)
        @test_throws DomainError compute_roche_lobe_correction(M_p, M_sun, a_orb, 0.0)
        # Rp >= R_Hill raises DomainError
        @test_throws DomainError compute_roche_lobe_correction(M_p, M_sun, 0.004 * AU, R_p)
    end

    # -------------------------------------------------------------------------
    # 3. Energy-Limited Hydrodynamic Escape Rate & Base Flux
    # -------------------------------------------------------------------------
    @testset "Energy-Limited Escape Rate & Base Flux" begin
        F_xuv = 10.0 # W/m^2
        epsilon = 0.15 # 15% efficiency

        res = compute_energy_limited_escape_flux(
            M_p, R_p, F_xuv; epsilon=epsilon, R_xuv=R_p, K_tide=1.0
        )

        @test res.M_dot_xuv > 0.0
        @test res.phi_xuv > 0.0

        # Analytical formula check:
        # M_dot = (epsilon * π * R_xuv^2 * F_xuv) / (G * M_p / R_p * K_tide)
        M_dot_ana = (epsilon * π * (R_p^2) * F_xuv) / (Erebus.G_GRAV * M_p / R_p)
        @test isapprox(res.M_dot_xuv, M_dot_ana; rtol=1e-12)

        phi_ana = M_dot_ana / (4.0 * π * (R_p^2))
        @test isapprox(res.phi_xuv, phi_ana; rtol=1e-12)

        # Proportionality: doubling F_xuv doubles escape rate
        res_2x = compute_energy_limited_escape_flux(
            M_p, R_p, 2.0 * F_xuv; epsilon=epsilon, R_xuv=R_p, K_tide=1.0
        )
        @test isapprox(res_2x.M_dot_xuv, 2.0 * res.M_dot_xuv; rtol=1e-12)

        # Inverse scaling with planet mass
        res_2M = compute_energy_limited_escape_flux(
            2.0 * M_p, R_p, F_xuv; epsilon=epsilon, R_xuv=R_p, K_tide=1.0
        )
        @test isapprox(res_2M.M_dot_xuv, 0.5 * res.M_dot_xuv; rtol=1e-12)

        # Domain error guards
        @test_throws DomainError compute_energy_limited_escape_flux(0.0, R_p, F_xuv)
        @test_throws DomainError compute_energy_limited_escape_flux(M_p, 0.0, F_xuv)
        @test_throws DomainError compute_energy_limited_escape_flux(M_p, R_p, -1.0)
        @test_throws DomainError compute_energy_limited_escape_flux(
            M_p, R_p, F_xuv; epsilon=0.0
        )
    end

    # -------------------------------------------------------------------------
    # 4. Multi-Species Fractionation & Crossover Drag Under XUV Wind
    # -------------------------------------------------------------------------
    @testset "Multi-Species Crossover Drag & Mass Fractionation" begin
        # Mixed atmosphere: H2 (carrier) + CO2 + N2 + H2O
        species_list = [:H2, :H2O, :N2, :CO2]
        m_H2 = Erebus.MASS_H2_KG
        m_H2O = Erebus.MASS_H2O_KG
        m_N2 = Erebus.MASS_N2_KG
        m_CO2 = Erebus.MASS_CO2_KG
        m_vec = [m_H2, m_H2O, m_N2, m_CO2]

        X_vec = [0.70, 0.15, 0.10, 0.05]
        T_exo = 1000.0 # Exobase temperature [K]
        g_exo = Erebus.G_GRAV * M_p / (R_p^2)
        b_mat = assemble_binary_diffusion_matrix(species_list, T_exo)

        # 1. Low flux regime: only H2 escapes, heavier species retained
        phi_low = 1.0e-12 # kg / (m^2 s)
        Phi_low = solve_multispecies_escape_closure(
            phi_low, X_vec, m_vec, T_exo, g_exo, b_mat
        )
        @test Phi_low[1] > 0.0 # H2 escapes
        @test iszero(Phi_low[4]) # CO2 retained

        # 2. High flux regime (vigorous XUV-driven blow-off): all species entrain
        phi_high = 1.0e-6 # kg / (m^2 s)
        Phi_high = solve_multispecies_escape_closure(
            phi_high, X_vec, m_vec, T_exo, g_exo, b_mat
        )
        @test all(Phi_high .> 0.0) # All species dragged along

        # 3. Fractionation factors
        frac = compute_escape_fractionation_factors(Phi_high, X_vec, species_list)
        @test frac isa Dict{Tuple{Symbol,Symbol},Float64}
        # H2 is lighter than CO2, so its escape is favored: alpha(H2, CO2) >= 1.0
        @test frac[(:H2, :CO2)] >= 1.0
    end

    # -------------------------------------------------------------------------
    # 5. Coupled Atmosphere Evolution Step with XUV Escape
    # -------------------------------------------------------------------------
    @testset "Coupled Atmosphere Step with XUV Escape" begin
        cfg_atm = AtmosphereConfig(; active=true, mode=:guillot)
        cfg_esc = EscapeConfig(;
            active=true,
            hydrodynamic=true,
            xuv_driven=true,
            epsilon_xuv=0.20,
            F_xuv_1au_sat=2.0,
        )

        atm_state = AtmosphereState()
        # Initialize atmosphere with 10^12 kg H2 and 10^12 kg CO2
        atm_state.M_atm[:H2] = 1.0e12
        atm_state.M_atm[:CO2] = 1.0e12

        dt = 3600.0 * 24.0 * 365.0 # 1 year
        vent_rates = Dict{Symbol,Float64}(:H2 => 0.0, :CO2 => 0.0)

        # Advance step
        evolve_coupled_atmosphere_step!(
            atm_state,
            vent_rates,
            dt,
            M_p,
            R_p,
            200.0,
            cfg_atm;
            escape_cfg=cfg_esc,
            sim_time_s=1.0e6 * Erebus.SEC_PER_YEAR,
            a_orb=a_orb,
            escape_active=true,
        )

        # Mass was lost to space via escape
        @test atm_state.M_atm[:H2] < 1.0e12
        @test atm_state.M_escaped[:H2] > 0.0
        # Exact mass conservation: initial = current + escaped
        @test isapprox(atm_state.M_atm[:H2] + atm_state.M_escaped[:H2], 1.0e12; rtol=1e-10)
        @test isapprox(
            atm_state.M_atm[:CO2] + atm_state.M_escaped[:CO2], 1.0e12; rtol=1e-10
        )
    end
end

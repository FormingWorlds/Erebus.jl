using Test
using Erebus
using Erebus.Config
using Erebus.Physics

@testset "HCNS Volatile Solubility, Speciation, and Saturation Ceilings" begin
    @testset "Extended Water Solubility Laws" begin
        p_100bar = 1.0e7 # 100 bar = 10 MPa

        # Burnham & Dixon standard basalt law
        w_bd = compute_water_solubility_melt(p_100bar; law=:burnham_dixon, As=0.40)
        # As * sqrt(10 MPa) = 0.40 * sqrt(10) ≈ 1.2649 wt%
        @test isapprox(w_bd, 0.40 * sqrt(10.0); rtol=1e-5)

        # Sossi et al. (2023) peridotite melt law
        w_sossi = compute_water_solubility_melt(p_100bar; law=:sossi_peridotite)
        @test w_sossi > 0.0
        # Expected: (524.0 * sqrt(100 bar)) * 1e-4 = 5240 * 1e-4 = 0.524 wt%
        @test isapprox(w_sossi, 0.524; rtol=1e-5)

        # Dixon et al. (1995) basalt law
        w_dixon = compute_water_solubility_melt(p_100bar; law=:basalt_dixon)
        @test w_dixon > 0.0
        # Expected: (965.0 * sqrt(100 bar)) * 1e-4 = 0.965 wt%
        @test isapprox(w_dixon, 0.965; rtol=1e-5)

        # Newcombe et al. (2017) lunar basalt law
        w_newcombe = compute_water_solubility_melt(p_100bar; law=:newcombe_lunar)
        @test w_newcombe > 0.0
        # Expected: (683.0 * sqrt(100 bar)) * 1e-4 = 0.683 wt%
        @test isapprox(w_newcombe, 0.683; rtol=1e-5)

        # Zero and negative pressure limit
        @test iszero(compute_water_solubility_melt(0.0; law=:sossi_peridotite))
        @test iszero(compute_water_solubility_melt(-1.0e5; law=:basalt_dixon))

        # DomainError guards
        @test_throws DomainError compute_water_solubility_melt(NaN; law=:sossi_peridotite)
        @test_throws DomainError compute_water_solubility_melt(Inf; law=:basalt_dixon)
        @test_throws DomainError compute_water_solubility_melt(
            -1.0; As=-5.0, law=:burnham_dixon
        )
        @test_throws ArgumentError compute_water_solubility_melt(p_100bar; law=:nonexistent)
    end

    @testset "Molecular H2 Solubility Laws" begin
        p_10bar = 1.0e6 # 10 bar = 1 MPa = 1.0e6 Pa

        # Hirschmann et al. (2012)
        h2_hirsch = compute_h2_solubility_melt(p_10bar; law=:hirschmann2012)
        # 10^(1.10083602 + 0.52413928 * log10(10)) = 10^1.6249753 ≈ 42.167 ppmw
        expected_hirsch = 10.0^(1.10083602 + 0.52413928)
        @test isapprox(h2_hirsch, expected_hirsch; rtol=1e-5)

        # Gaillard et al. (2003)
        h2_gaillard = compute_h2_solubility_melt(p_10bar; law=:gaillard2003)
        # 0.163 * (10.0^1.252) ≈ 2.91197 ppmw
        expected_gaillard = 0.163 * (10.0^1.252)
        @test isapprox(h2_gaillard, expected_gaillard; rtol=1e-5)

        # Both positive and finite
        @test h2_hirsch > 0.0
        @test h2_gaillard > 0.0

        # Zero and negative pressure
        @test iszero(compute_h2_solubility_melt(0.0; law=:hirschmann2012))
        @test iszero(compute_h2_solubility_melt(-1.0e5; law=:gaillard2003))

        # Guards
        @test_throws DomainError compute_h2_solubility_melt(NaN; law=:hirschmann2012)
        @test_throws DomainError compute_h2_solubility_melt(Inf; law=:hirschmann2012)
        @test_throws ArgumentError compute_h2_solubility_melt(p_10bar; law=:invalid)
    end

    @testset "Dasgupta et al. (2022) Compositional Nitrogen Solubility" begin
        p_100bar = 1.0e7
        p_tot_100bar = 1.0e7
        T_1600K = 1600.0
        d_IW = -2.0

        n_dasgupta = compute_nitrogen_solubility_dasgupta(
            p_100bar, p_tot_100bar, T_1600K, d_IW; x_SiO2=0.50, x_Al2O3=0.15, x_TiO2=0.02
        )
        @test n_dasgupta.total_ppm > 0.0
        @test n_dasgupta.physical_ppm > 0.0
        @test n_dasgupta.chemical_ppm > 0.0
        @test isapprox(
            n_dasgupta.total_ppm,
            n_dasgupta.physical_ppm + n_dasgupta.chemical_ppm;
            rtol=1e-10,
        )

        # Reducing conditions increase nitride solubility
        n_more_reduced = compute_nitrogen_solubility_dasgupta(
            p_100bar, p_tot_100bar, T_1600K, -4.0; x_SiO2=0.50, x_Al2O3=0.15, x_TiO2=0.02
        )
        @test n_more_reduced.chemical_ppm > n_dasgupta.chemical_ppm
        @test isapprox(n_more_reduced.physical_ppm, n_dasgupta.physical_ppm; rtol=1e-10)

        # Zero pressure
        n_zero = compute_nitrogen_solubility_dasgupta(0.0, p_tot_100bar, T_1600K, d_IW)
        @test iszero(n_zero.total_ppm)
        @test iszero(n_zero.physical_ppm)
        @test iszero(n_zero.chemical_ppm)

        # DomainError guards
        @test_throws DomainError compute_nitrogen_solubility_dasgupta(
            p_100bar, p_tot_100bar, 0.0, d_IW
        )
        @test_throws DomainError compute_nitrogen_solubility_dasgupta(
            p_100bar, p_tot_100bar, -100.0, d_IW
        )
        @test_throws DomainError compute_nitrogen_solubility_dasgupta(
            NaN, p_tot_100bar, T_1600K, d_IW
        )
        @test_throws DomainError compute_nitrogen_solubility_dasgupta(
            p_100bar, p_tot_100bar, T_1600K, NaN
        )
        @test_throws DomainError compute_nitrogen_solubility_dasgupta(
            p_100bar, p_tot_100bar, T_1600K, d_IW; x_SiO2=-0.1
        )
    end

    @testset "Carbon Species Solubility Laws (CO, CH4, CO2)" begin
        p_co_Pa = 1.0e6   # 10 bar
        p_ch4_Pa = 1.0e6  # 10 bar
        p_co2_Pa = 1.0e6  # 10 bar
        p_tot_Pa = 5.0e6  # 50 bar
        T_1500K = 1500.0

        # CO: Armstrong et al. (2015) vs Yoshioka et al. (2019)
        co_arm = compute_co_solubility_melt(p_co_Pa, p_tot_Pa; law=:armstrong2015)
        co_yosh = compute_co_solubility_melt(p_co_Pa, p_tot_Pa; law=:yoshioka2019_morb)
        @test co_arm > 0.0
        @test co_yosh > 0.0
        # Expected Armstrong: log10 = -0.738 + 0.876 * log10(10) - 5.44e-5 * 50 = -0.738 + 0.876 - 0.00272 = 0.13528
        expected_co = 10.0^(-0.738 + 0.876 * 1.0 - 5.44e-5 * 50.0)
        @test isapprox(co_arm, expected_co; rtol=1e-5)

        # CH4: Ardia et al. (2013)
        ch4_ard = compute_ch4_solubility_melt(p_ch4_Pa, p_tot_Pa; law=:ardia2013)
        @test ch4_ard > 0.0
        # Expected Ardia: p_ch4_gpa * exp(4.93 - 1.93 * p_tot_gpa)
        expected_ch4 = (1.0e6 * 1.0e-9) * exp(4.93 - 1.93 * (5.0e6 * 1.0e-9))
        @test isapprox(ch4_ard, expected_ch4; rtol=1e-5)

        # CO2: Dixon et al. (1995)
        co2_dix = compute_co2_solubility_melt(p_co2_Pa, T_1500K; law=:dixon1995)
        @test co2_dix > 0.0
        # Expected Dixon: x = 3.8e-7 * 10 * exp(-23 * 9 / (83.15 * 1500))
        x_co2 = 3.8e-7 * 10.0 * exp(-23.0 * 9.0 / (83.15 * 1500.0))
        expected_co2 = 1.0e4 * (4400.0 * x_co2) / (36.6 - 44.0 * x_co2)
        @test isapprox(co2_dix, expected_co2; rtol=1e-5)

        # Composite carbon solubility
        c_comp = compute_carbon_solubility_melt(
            p_co_Pa, p_ch4_Pa, p_co2_Pa, p_tot_Pa, T_1500K
        )
        @test isapprox(
            c_comp.total_ppm, c_comp.co_ppm + c_comp.ch4_ppm + c_comp.co2_ppm; rtol=1e-10
        )
        @test isapprox(c_comp.co_ppm, co_arm; rtol=1e-10)
        @test isapprox(c_comp.ch4_ppm, ch4_ard; rtol=1e-10)
        @test isapprox(c_comp.co2_ppm, co2_dix; rtol=1e-10)

        # Graphite saturation capping in composite carbon solubility
        c_uncapped = compute_carbon_solubility_melt(
            1.0e9, 1.0e5, 1.0e9, p_tot_Pa, T_1500K; graphite_saturation=false
        )
        c_capped = compute_carbon_solubility_melt(
            1.0e9, 1.0e5, 1.0e9, p_tot_Pa, T_1500K; graphite_saturation=true, delta_IW=-1.0
        )
        @test c_capped.total_ppm < c_uncapped.total_ppm
        @test c_capped.co_ppm < c_uncapped.co_ppm
        @test c_capped.co2_ppm < c_uncapped.co2_ppm

        # Zero pressure limits
        @test iszero(compute_co_solubility_melt(0.0, p_tot_Pa))
        @test iszero(compute_ch4_solubility_melt(0.0, p_tot_Pa))
        @test iszero(compute_co2_solubility_melt(0.0, T_1500K))

        # DomainError guards
        @test_throws DomainError compute_co_solubility_melt(NaN, p_tot_Pa)
        @test_throws DomainError compute_ch4_solubility_melt(NaN, p_tot_Pa)
        @test_throws DomainError compute_co2_solubility_melt(p_co2_Pa, 0.0)
        @test_throws DomainError compute_co2_solubility_melt(p_co2_Pa, -100.0)
        @test_throws ArgumentError compute_co_solubility_melt(
            p_co_Pa, p_tot_Pa; law=:unknown
        )
    end

    @testset "Graphite Saturation Ceiling (French 1966 / Holloway 1992)" begin
        T_1500K = 1500.0
        log10_fO2_reduced = -15.0
        log10_fO2_oxidized = -5.0

        gr_red = compute_graphite_saturation_fugacity(T_1500K, log10_fO2_reduced)
        @test gr_red.f_CO_max_bar > 0.0
        @test gr_red.f_CO2_max_bar > 0.0
        # Expected:
        # log10(f_CO) = 5785 / 1500 + 4.545 + 0.5 * (-15) = 3.85667 + 4.545 - 7.5 = 0.90167
        # log10(f_CO2) = 20590 / 1500 - 0.043 + (-15) = 13.72667 - 0.043 - 15 = -1.31633
        @test isapprox(log10(gr_red.f_CO_max_bar), 5785.0 / 1500.0 + 4.545 - 7.5; rtol=1e-5)
        @test isapprox(
            log10(gr_red.f_CO2_max_bar), 20590.0 / 1500.0 - 0.043 - 15.0; rtol=1e-5
        )

        # Under oxidized conditions, CO2 ceiling rises, CO ceiling rises as fO2^0.5
        gr_ox = compute_graphite_saturation_fugacity(T_1500K, log10_fO2_oxidized)
        @test gr_ox.f_CO2_max_bar > gr_red.f_CO2_max_bar
        @test gr_ox.f_CO_max_bar > gr_red.f_CO_max_bar

        # Guards
        @test_throws DomainError compute_graphite_saturation_fugacity(0.0, -10.0)
        @test_throws DomainError compute_graphite_saturation_fugacity(-300.0, -10.0)
        @test_throws DomainError compute_graphite_saturation_fugacity(T_1500K, NaN)
        @test_throws DomainError compute_graphite_saturation_fugacity(T_1500K, Inf)
    end

    @testset "Sulfur Solubility Laws (Boulliung & Wood 2023, Gaillard 2022)" begin
        T_1500K = 1500.0
        p_S2_Pa = 1.0e5 # 1 bar
        d_IW = -1.0

        # Boulliung 2023 basalt
        s_boul = compute_sulfur_solubility_melt(
            p_S2_Pa, T_1500K, d_IW; law=:boulliung2023, sulfide_melt=:basalt
        )
        @test s_boul > 0.0

        # Different compositions have different capacities
        s_and = compute_sulfur_solubility_melt(
            p_S2_Pa, T_1500K, d_IW; law=:boulliung2023, sulfide_melt=:andesite
        )
        @test s_and > 0.0
        @test s_boul > s_and

        # Adding sulfate capacity
        s_with_sulfate = compute_sulfur_solubility_melt(
            p_S2_Pa,
            T_1500K,
            2.0;
            law=:boulliung2023,
            sulfide_melt=:basalt,
            include_sulfate=true,
        )
        s_no_sulfate = compute_sulfur_solubility_melt(
            p_S2_Pa,
            T_1500K,
            2.0;
            law=:boulliung2023,
            sulfide_melt=:basalt,
            include_sulfate=false,
        )
        @test s_with_sulfate > s_no_sulfate

        # Gaillard 2022 law
        s_gail = compute_sulfur_solubility_melt(
            p_S2_Pa, T_1500K, d_IW; law=:gaillard2022, x_FeO=10.0
        )
        @test s_gail > 0.0

        # SCSS capping in sulfur solubility
        s_uncapped = compute_sulfur_solubility_melt(
            1.0e7, T_1500K, -1.0; law=:boulliung2023, scss_active=false
        )
        s_capped = compute_sulfur_solubility_melt(
            1.0e7,
            T_1500K,
            -1.0;
            law=:boulliung2023,
            scss_active=true,
            p_total_Pa=1.0e7,
            x_FeO=10.0,
        )
        scss_cap = compute_scss(T_1500K, 1.0e7; x_FeO=10.0, law=:smythe2017)
        @test s_capped <= scss_cap
        @test isapprox(s_capped, min(s_uncapped, scss_cap); rtol=1e-10)

        # Zero pressure limit
        @test iszero(compute_sulfur_solubility_melt(0.0, T_1500K, d_IW))

        # Guards
        @test_throws DomainError compute_sulfur_solubility_melt(p_S2_Pa, 0.0, d_IW)
        @test_throws DomainError compute_sulfur_solubility_melt(p_S2_Pa, T_1500K, 100.0)
        @test_throws DomainError compute_sulfur_solubility_melt(p_S2_Pa, T_1500K, NaN)
        @test_throws DomainError compute_sulfur_solubility_melt(NaN, T_1500K, d_IW)
        @test_throws ArgumentError compute_sulfur_solubility_melt(
            p_S2_Pa, T_1500K, d_IW; law=:invalid
        )
    end

    @testset "SCSS Ceilings (Smythe et al. 2017 / O'Neill & Mavrogenes 2002)" begin
        T_1500K = 1500.0
        p_100bar = 1.0e7

        scss_val = compute_scss(T_1500K, p_100bar; x_FeO=10.0, law=:smythe2017)
        @test scss_val > 500.0
        @test scss_val < 3000.0
        # Expected: exp(7.50 - 4500 / 1500 + 0.90 * log(10) - 2.5e-4 * 100 / 1500)
        # = exp(7.50 - 3.0 + 2.07233 - 1.6667e-5) ≈ exp(6.5723) ≈ 715 ppmw
        expected_scss = exp(
            7.50 - 4500.0 / 1500.0 + 0.90 * log(10.0) - 2.5e-4 * (100.0 / 1500.0)
        )
        @test isapprox(scss_val, expected_scss; rtol=1e-5)

        # Higher FeO increases SCSS
        scss_high_fe = compute_scss(T_1500K, p_100bar; x_FeO=20.0, law=:smythe2017)
        @test scss_high_fe > scss_val

        # Higher temperature increases SCSS
        scss_hot = compute_scss(1700.0, p_100bar; x_FeO=10.0, law=:smythe2017)
        @test scss_hot > scss_val

        # Differentiation between Smythe et al. (2017) and O'Neill & Mavrogenes (2002)
        p_high = 1.0e8 # 1000 bar
        scss_smythe = compute_scss(T_1500K, p_high; x_FeO=10.0, law=:smythe2017)
        scss_oneill = compute_scss(T_1500K, p_high; x_FeO=10.0, law=:oneill2002)
        @test scss_oneill > scss_smythe
        @test isapprox(
            scss_oneill, exp(7.50 - 4500.0 / 1500.0 + 0.90 * log(10.0)); rtol=1e-5
        )

        # Guards
        @test_throws DomainError compute_scss(0.0, p_100bar)
        @test_throws DomainError compute_scss(-100.0, p_100bar)
        @test_throws DomainError compute_scss(T_1500K, NaN)
        @test_throws DomainError compute_scss(T_1500K, p_100bar; x_FeO=-5.0)
        @test_throws ArgumentError compute_scss(T_1500K, p_100bar; law=:unknown)
    end

    @testset "Gas Speciation Solver (solve_chnos_speciation)" begin
        p_tot = 1.0e7 # 100 bar = 10 MPa
        T_1500K = 1500.0

        # Reducing case: IW - 2
        spec_red = solve_chnos_speciation(
            p_tot, T_1500K, -2.0; z_H=0.80, z_C=0.15, z_N=0.03, z_S=0.02
        )
        @test spec_red.p_H2_Pa > 0.0
        @test spec_red.p_H2O_Pa > 0.0
        @test spec_red.p_CO_Pa > 0.0
        @test spec_red.p_CO2_Pa > 0.0
        @test spec_red.p_N2_Pa > 0.0
        @test spec_red.p_H2S_Pa > 0.0

        # Under reducing conditions: H2 > H2O, CO > CO2
        @test spec_red.p_H2_Pa > spec_red.p_H2O_Pa
        @test spec_red.p_CO_Pa > spec_red.p_CO2_Pa

        # Oxidizing case: IW + 2
        spec_ox = solve_chnos_speciation(
            p_tot, T_1500K, 2.0; z_H=0.80, z_C=0.15, z_N=0.03, z_S=0.02
        )
        @test spec_ox.p_H2O_Pa > spec_ox.p_H2_Pa
        @test spec_ox.p_CO2_Pa > spec_ox.p_CO_Pa
        @test spec_ox.p_SO2_Pa > spec_red.p_SO2_Pa

        # Elemental conservation checks
        pN_tot = 0.03 * p_tot
        pS_tot = 0.02 * p_tot
        # Nitrogen conservation: 2 * p_N2 + p_NH3 = pN_tot
        @test isapprox(2.0 * spec_red.p_N2_Pa + spec_red.p_NH3_Pa, pN_tot; rtol=1e-10)
        @test isapprox(2.0 * spec_ox.p_N2_Pa + spec_ox.p_NH3_Pa, pN_tot; rtol=1e-10)
        # Sulfur conservation: 2 * p_S2 + p_H2S + p_SO2 = pS_tot
        @test isapprox(
            2.0 * spec_red.p_S2_Pa + spec_red.p_H2S_Pa + spec_red.p_SO2_Pa,
            pS_tot;
            rtol=1e-10,
        )
        @test isapprox(
            2.0 * spec_ox.p_S2_Pa + spec_ox.p_H2S_Pa + spec_ox.p_SO2_Pa, pS_tot; rtol=1e-10
        )

        # Low temperature stability (no NaN, all non-negative and finite)
        spec_cold5 = solve_chnos_speciation(1.0e7, 5.0, 0.0)
        spec_cold50 = solve_chnos_speciation(1.0e7, 50.0, 0.0)
        for sp in (spec_cold5, spec_cold50)
            @test isfinite(sp.p_H2_Pa) && sp.p_H2_Pa >= 0.0
            @test isfinite(sp.p_H2O_Pa) && sp.p_H2O_Pa >= 0.0
            @test isfinite(sp.p_CO_Pa) && sp.p_CO_Pa >= 0.0
            @test isfinite(sp.p_CO2_Pa) && sp.p_CO2_Pa >= 0.0
            @test isfinite(sp.p_CH4_Pa) && sp.p_CH4_Pa >= 0.0
            @test isfinite(sp.p_N2_Pa) && sp.p_N2_Pa >= 0.0
            @test isfinite(sp.p_NH3_Pa) && sp.p_NH3_Pa >= 0.0
            @test isfinite(sp.p_H2S_Pa) && sp.p_H2S_Pa >= 0.0
            @test isfinite(sp.p_S2_Pa) && sp.p_S2_Pa >= 0.0
            @test isfinite(sp.p_SO2_Pa) && sp.p_SO2_Pa >= 0.0
        end

        # Zero pressure limit
        spec_zero = solve_chnos_speciation(0.0, T_1500K, 0.0)
        @test iszero(spec_zero.p_H2_Pa)
        @test iszero(spec_zero.p_H2O_Pa)
        @test iszero(spec_zero.p_CO_Pa)

        # DomainError guards
        @test_throws DomainError solve_chnos_speciation(NaN, T_1500K, 0.0)
        @test_throws DomainError solve_chnos_speciation(p_tot, 0.0, 0.0)
        @test_throws DomainError solve_chnos_speciation(p_tot, T_1500K, 100.0)
        @test_throws DomainError solve_chnos_speciation(
            p_tot, T_1500K, 0.0; z_H=0.0, z_C=0.0, z_N=0.0, z_S=0.0
        )
    end

    @testset "Species Molecular Masses and Constants" begin
        @test isapprox(MASS_H2_KG, 2.01588 / (6.02214076e23 * 1000.0); rtol=1e-4)
        @test isapprox(MASS_CH4_KG, 16.04246 / (6.02214076e23 * 1000.0); rtol=1e-4)
        @test isapprox(MASS_H2S_KG, 34.08088 / (6.02214076e23 * 1000.0); rtol=1e-4)
        @test isapprox(MASS_S2_KG, 64.130 / (6.02214076e23 * 1000.0); rtol=1e-4)
        @test isapprox(MASS_SO2_KG, 64.066 / (6.02214076e23 * 1000.0); rtol=1e-4)

        # get_species_molecular_mass lookups
        @test isapprox(get_species_molecular_mass(:H2), MASS_H2_KG; rtol=1e-12)
        @test isapprox(get_species_molecular_mass(:CH4), MASS_CH4_KG; rtol=1e-12)
        @test isapprox(get_species_molecular_mass(:H2S), MASS_H2S_KG; rtol=1e-12)
        @test isapprox(get_species_molecular_mass(:S2), MASS_S2_KG; rtol=1e-12)
        @test isapprox(get_species_molecular_mass(:SO2), MASS_SO2_KG; rtol=1e-12)
        @test isapprox(get_species_molecular_mass(:h2), MASS_H2_KG; rtol=1e-12)

        # Unknown species guard
        @test_throws ArgumentError get_species_molecular_mass(:X_UNKNOWN)
        @test_throws ArgumentError get_species_molecular_mass(:Argon)
    end

    @testset "VolatilesConfig Extended Parameters & Validation" begin
        # Default config
        cfg_def = VolatilesConfig()
        @test cfg_def.water_law === :burnham_dixon
        @test cfg_def.h2_active == false
        @test cfg_def.h2_law === :hirschmann2012
        @test cfg_def.carbon_active == false
        @test cfg_def.co_law === :armstrong2015
        @test cfg_def.ch4_law === :ardia2013
        @test cfg_def.co2_law === :dixon1995
        @test cfg_def.graphite_saturation == true
        @test cfg_def.sulfur_active == false
        @test cfg_def.sulfide_law === :boulliung2023
        @test cfg_def.sulfide_melt === :basalt
        @test cfg_def.include_sulfate == false
        @test cfg_def.scss_active == true
        @test cfg_def.scss_law === :smythe2017
        @test isapprox(cfg_def.melt_feo_wtpct, 10.0; rtol=1e-12)

        # Round-trip serialization
        sim_cfg = SimulationConfig(
            volatiles=VolatilesConfig(
                water_law=:sossi_peridotite,
                h2_active=true,
                h2_law=:gaillard2003,
                carbon_active=true,
                co_law=:yoshioka2019_morb,
                sulfur_active=true,
                sulfide_law=:gaillard2022,
                include_sulfate=true,
                melt_feo_wtpct=12.5,
            ),
        )
        toml_str = save_config(sim_cfg)
        loaded = load_config(toml_str)
        @test loaded.volatiles.water_law === :sossi_peridotite
        @test loaded.volatiles.h2_active == true
        @test loaded.volatiles.h2_law === :gaillard2003
        @test loaded.volatiles.carbon_active == true
        @test loaded.volatiles.co_law === :yoshioka2019_morb
        @test loaded.volatiles.sulfur_active == true
        @test loaded.volatiles.sulfide_law === :gaillard2022
        @test loaded.volatiles.include_sulfate == true
        @test isapprox(loaded.volatiles.melt_feo_wtpct, 12.5; rtol=1e-12)

        # Validation errors
        cfg_base = default_config()
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                grid=cfg_base.grid,
                geometry=cfg_base.geometry,
                time=cfg_base.time,
                solver=cfg_base.solver,
                poroelasticity=cfg_base.poroelasticity,
                thermodynamics=cfg_base.thermodynamics,
                reaction=cfg_base.reaction,
                materials=cfg_base.materials,
                output=cfg_base.output,
                disk=cfg_base.disk,
                melting=cfg_base.melting,
                venting=cfg_base.venting,
                volatiles=VolatilesConfig(water_law=:bad_law),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                grid=cfg_base.grid,
                geometry=cfg_base.geometry,
                time=cfg_base.time,
                solver=cfg_base.solver,
                poroelasticity=cfg_base.poroelasticity,
                thermodynamics=cfg_base.thermodynamics,
                reaction=cfg_base.reaction,
                materials=cfg_base.materials,
                output=cfg_base.output,
                disk=cfg_base.disk,
                melting=cfg_base.melting,
                venting=cfg_base.venting,
                volatiles=VolatilesConfig(melt_feo_wtpct=-1.0),
            ),
        )
    end
end

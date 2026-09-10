using Test
using LinearAlgebra
using Erebus
using Erebus.Config
using Erebus.Physics
using Erebus.Particles

@testset "HCNSPO Volatile Mixtures & Refractory Phases" begin
    # =========================================================================
    # 1. Configuration & Schema Validation
    # =========================================================================
    @testset "RefractoryConfig & VolatileMixtureConfig Schema" begin
        # Default construction
        refr_cfg = RefractoryConfig()
        @test 0.50 <= refr_cfg.f_refr_C <= 0.80   # Bergin et al. (2026) IOM / dust range
        @test 0.80 <= refr_cfg.f_refr_S <= 0.95   # Kama et al. (2019) FeS sulfide range
        @test 0.05 <= refr_cfg.f_refr_N <= 0.25   # Refractory organic N / nitrides
        @test 0.90 <= refr_cfg.f_refr_P <= 1.00   # Schreibersite / phosphates
        @test refr_cfg.T_pyrolysis_C >= 500.0
        @test refr_cfg.T_dehydrate_H >= 650.0

        mix_cfg = VolatileMixtureConfig()
        @test mix_cfg.X_ice_H2O > 0.50
        @test mix_cfg.T_cond_H2O >
            mix_cfg.T_cond_NH3 >
            mix_cfg.T_cond_CO2 >=
            mix_cfg.T_cond_H2S >
            mix_cfg.T_cond_CH4 >
            mix_cfg.T_cond_CO >
            mix_cfg.T_cond_N2
        @test isapprox(mix_cfg.T_eutectic_ammonia, 176.0, atol=1e-6)
        @test isapprox(mix_cfg.lambda_nh3_depression, (273.15 - 176.0) / 0.33, atol=1e-6)
        @test isapprox(mix_cfg.alpha_P, 0.0, atol=1e-6)
        @test isapprox(mix_cfg.P_ref, 1.0, atol=1e-6)
        @test_throws DomainError VolatileMixtureConfig(alpha_P=-0.05)
        @test_throws DomainError VolatileMixtureConfig(P_ref=0.0)

        # Custom valid construction
        custom_refr = RefractoryConfig(
            active=true,
            f_refr_C=0.65,
            f_refr_N=0.15,
            f_refr_S=0.88,
            f_refr_P=0.99,
            f_refr_H=0.08,
            T_pyrolysis_C=650.0,
            T_dehydrate_H=780.0,
        )
        @test custom_refr.active == true
        @test isapprox(custom_refr.f_refr_C, 0.65, atol=1e-6)
        @test isapprox(custom_refr.f_refr_S, 0.88, atol=1e-6)

        custom_mix = VolatileMixtureConfig(
            active=true,
            X_ice_H2O=0.80,
            X_ice_CO2=0.10,
            X_ice_CO=0.02,
            X_ice_CH4=0.02,
            X_ice_NH3=0.04,
            X_ice_N2=0.01,
            X_ice_H2S=0.01,
            X_ice_PH3=0.0,
            T_cond_H2O=165.0,
            T_cond_NH3=130.0,
            T_cond_CO2=78.0,
            T_cond_H2S=75.0,
            T_cond_CH4=42.0,
            T_cond_CO=24.0,
            T_cond_N2=17.0,
        )
        @test custom_mix.active == true
        @test isapprox(custom_mix.X_ice_NH3, 0.04, atol=1e-6)

        # Validation constraints: fractions must be in [0, 1]
        @test_throws DomainError RefractoryConfig(f_refr_C=-0.1)
        @test_throws DomainError RefractoryConfig(f_refr_C=1.1)
        @test_throws DomainError RefractoryConfig(f_refr_S=-0.05)
        @test_throws DomainError RefractoryConfig(f_refr_S=1.05)
        @test_throws DomainError RefractoryConfig(f_refr_N=-0.01)
        @test_throws DomainError RefractoryConfig(f_refr_N=1.01)
        @test_throws DomainError RefractoryConfig(f_refr_P=-0.01)
        @test_throws DomainError RefractoryConfig(f_refr_P=1.01)
        @test_throws DomainError RefractoryConfig(T_pyrolysis_C=-100.0)
        @test_throws DomainError RefractoryConfig(T_dehydrate_H=-50.0)

        @test_throws DomainError VolatileMixtureConfig(X_ice_H2O=-0.1)
        @test_throws DomainError VolatileMixtureConfig(X_ice_H2O=1.5)
        @test_throws DomainError VolatileMixtureConfig(X_ice_NH3=-0.01)
        @test_throws DomainError VolatileMixtureConfig(T_eutectic_ammonia=-10.0)
        @test_throws DomainError VolatileMixtureConfig(T_cond_H2O=-160.0)
        # Unphysical inversion of snowline order across species
        @test_throws ArgumentError VolatileMixtureConfig(T_cond_H2O=100.0, T_cond_CO2=150.0)
        @test_throws ArgumentError VolatileMixtureConfig(T_cond_H2S=140.0, T_cond_NH3=135.0)
        @test_throws ArgumentError VolatileMixtureConfig(T_cond_PH3=50.0, T_cond_CH4=45.0)
    end

    # =========================================================================
    # 2. Pore Fluid Freezing Point & Eutectic Thermodynamics
    # =========================================================================
    @testset "Ammonia-Water Freezing Point Depression" begin
        # Pure water limit (X_nh3 = 0.0): T_freeze == 273.15 K
        T_freeze_pure = compute_mixture_freezing_point(0.0)
        @test isapprox(T_freeze_pure, 273.15, atol=1e-6)

        # Dilute ammonia solution (X_nh3 = 0.05): continuous depression
        T_freeze_dilute = compute_mixture_freezing_point(0.05)
        expected_dilute = 273.15 - ((273.15 - 176.0) / 0.33) * 0.05
        @test isapprox(T_freeze_dilute, expected_dilute, atol=1e-2)
        @test T_freeze_dilute < T_freeze_pure

        # Continuous transition at eutectic composition without discontinuity
        T_below_eutectic = compute_mixture_freezing_point(0.329)
        T_at_eutectic = compute_mixture_freezing_point(0.330)
        T_above_eutectic = compute_mixture_freezing_point(0.331)
        @test abs(T_at_eutectic - T_below_eutectic) < 0.5
        @test isapprox(T_at_eutectic, 176.0, atol=1e-2)
        @test isapprox(T_above_eutectic, 176.0, atol=1e-6)

        # Monotonicity test across fine range X_nh3 in [0, 0.4]
        nh3_sweep = range(0.0, 0.40, length=100)
        T_sweep = [compute_mixture_freezing_point(x) for x in nh3_sweep]
        for k in 1:(length(T_sweep) - 1)
            @test T_sweep[k] >= T_sweep[k + 1] - 1e-12
            @test T_sweep[k] >= 176.0 - 1e-12
        end

        # Test passing custom VolatileMixtureConfig
        custom_mix = VolatileMixtureConfig(
            lambda_nh3_depression=200.0, lambda_solute_depression=60.0, T_freeze_floor=180.0
        )
        T_cfg = compute_mixture_freezing_point(0.10, 0.05; cfg=custom_mix)
        @test isapprox(T_cfg, 273.15 - 200.0 * 0.10 - 60.0 * 0.05, atol=1e-2)
        T_cfg_floor = compute_mixture_freezing_point(0.80; cfg=custom_mix)
        @test isapprox(T_cfg_floor, 180.0, atol=1e-6)

        # Domain errors
        @test_throws DomainError compute_mixture_freezing_point(-0.01)
        @test_throws DomainError compute_mixture_freezing_point(1.05)
        @test_throws DomainError compute_mixture_freezing_point(0.1, -0.05)
        @test_throws DomainError compute_mixture_freezing_point(NaN)
    end

    # =========================================================================
    # 3. Fluid Mixture Density & Viscosity
    # =========================================================================
    @testset "Multi-Component Fluid EOS & Viscosity" begin
        T_test = 300.0
        P_test = 1.0e6   # 1 MPa

        # Pure water density at 300 K, 1 MPa ~ 1000 kg/m^3
        rho_pure = compute_mixture_fluid_density(T_test, P_test, 0.0)
        @test 990.0 < rho_pure < 1010.0

        # Ammonia is less dense than liquid water (~680 kg/m^3 at 300 K vs 1000 kg/m^3)
        # Therefore X_nh3 = 0.10 must reduce fluid density
        rho_mix = compute_mixture_fluid_density(T_test, P_test, 0.10)
        @test rho_mix < rho_pure
        @test rho_mix > 900.0
        # Check expected reduction: rho_mix ~ rho_pure * (1 - 0.25 * 0.10) ~ 975 kg/m^3
        @test isapprox(rho_mix, rho_pure * (1.0 - 0.25 * 0.10), rtol=0.02)

        # Sub-freezing ice viscosity: T < T_freeze yields eta_ice (1.0e12 Pa s)
        T_subfreezing = 200.0
        eta_ice_pure = compute_mixture_fluid_viscosity(
            T_subfreezing, 0.0; T_melt=273.15, eta_ice=1.0e12
        )
        @test isapprox(eta_ice_pure, 1.0e12, rtol=1e-10)

        # For an ammonia-rich fluid (X_nh3 = 0.20), T_melt ~ 243 K.
        # At T = 250 K, pure water is frozen (T < 273.15 K -> eta = 1e12),
        # BUT ammonia-water is LIQUID (T > 243 K -> eta ~ 1e-3 Pa s)!
        eta_h2o_250k = compute_mixture_fluid_viscosity(
            250.0, 0.0; T_melt=273.15, eta_ice=1.0e12
        )
        eta_mix_250k = compute_mixture_fluid_viscosity(
            250.0, 0.20; T_melt=243.15, eta_ice=1.0e12
        )
        @test isapprox(eta_h2o_250k, 1.0e12, rtol=1e-10)
        @test eta_mix_250k < 0.10   # Liquid mobility enabled at 250 K!
        @test eta_mix_250k > 1.0e-4

        # Domain error contracts
        @test_throws DomainError compute_mixture_fluid_density(-10.0, P_test, 0.0)
        @test_throws DomainError compute_mixture_fluid_density(T_test, -1.0e5, 0.0)
        @test_throws DomainError compute_mixture_fluid_density(T_test, P_test, -0.1)
        @test_throws DomainError compute_mixture_fluid_viscosity(-5.0, 0.0)
        @test_throws DomainError compute_mixture_fluid_viscosity(T_test, 1.2)
    end

    # =========================================================================
    # 4. Multi-Snowline Disk Accretion Condensation
    # =========================================================================
    @testset "Multi-Snowline Disk Condensation Model" begin
        mix_cfg = VolatileMixtureConfig(active=true)
        refr_cfg = RefractoryConfig(
            active=true, f_refr_C=0.60, f_refr_S=0.89, f_refr_N=0.10, f_refr_P=0.98
        )
        P_midplane = 10.0   # 10 Pa

        # Case A: Warm inner disk (T = 300 K > T_snowline_H2O = 160 K)
        # Inside the water snowline: NO volatile ices condense!
        # Refractory carbon, sulfur, nitrogen, and phosphorus ARE delivered by refractory grains!
        state_inner = evaluate_disk_volatile_condensation(
            300.0, P_midplane, mix_cfg, refr_cfg
        )
        @test state_inner.condensed_H2O == false
        @test state_inner.condensed_NH3 == false
        @test state_inner.condensed_CO2 == false
        @test state_inner.condensed_CH4 == false
        @test state_inner.condensed_CO == false
        @test state_inner.condensed_N2 == false
        @test iszero(state_inner.X_ice_H2O)
        @test iszero(state_inner.X_ice_CO2)
        # Refractory fractions preserved
        @test isapprox(state_inner.f_refr_C, 0.60, atol=1e-10)
        @test isapprox(state_inner.f_refr_S, 0.89, atol=1e-10)
        @test isapprox(state_inner.f_refr_N, 0.10, atol=1e-10)
        @test isapprox(state_inner.f_refr_P, 0.98, atol=1e-10)

        # Case B: Intermediate disk (T = 120 K, outside H2O snowline, inside CO2/CH4 snowlines)
        # H2O and NH3 ices condense, but CO2, CH4, CO, N2 remain gas!
        state_inter = evaluate_disk_volatile_condensation(
            120.0, P_midplane, mix_cfg, refr_cfg
        )
        @test state_inter.condensed_H2O == true
        @test state_inter.condensed_NH3 == true
        @test state_inter.condensed_CO2 == false
        @test state_inter.condensed_CH4 == false
        @test state_inter.condensed_CO == false
        @test state_inter.condensed_N2 == false
        @test state_inter.X_ice_H2O > 0.0
        @test state_inter.X_ice_NH3 > 0.0
        @test iszero(state_inter.X_ice_CO2)
        @test iszero(state_inter.X_ice_CH4)

        # Case C: Cold intermediate disk (T = 60 K, outside CO2 and H2S snowlines)
        # H2O, NH3, CO2, and H2S condense!
        state_cold = evaluate_disk_volatile_condensation(
            60.0, P_midplane, mix_cfg, refr_cfg
        )
        @test state_cold.condensed_H2O == true
        @test state_cold.condensed_NH3 == true
        @test state_cold.condensed_CO2 == true
        @test state_cold.condensed_H2S == true
        @test state_cold.condensed_CH4 == false
        @test state_cold.condensed_CO == false
        @test state_cold.X_ice_CO2 > 0.0
        @test state_cold.X_ice_H2S > 0.0
        @test iszero(state_cold.X_ice_CH4)

        # Case D: Cryogenic outer disk (T = 10 K, beyond CO and N2 snowlines)
        # All volatile ices condense!
        state_cryo = evaluate_disk_volatile_condensation(
            10.0, P_midplane, mix_cfg, refr_cfg
        )
        @test state_cryo.condensed_H2O == true
        @test state_cryo.condensed_NH3 == true
        @test state_cryo.condensed_CO2 == true
        @test state_cryo.condensed_H2S == true
        @test state_cryo.condensed_CH4 == true
        @test state_cryo.condensed_CO == true
        @test state_cryo.condensed_N2 == true
        @test state_cryo.X_ice_CO > 0.0
        @test state_cryo.X_ice_N2 > 0.0
        @test state_cryo.X_ice_CH4 > 0.0

        # Case E: Extreme cryogenic disk (T = 5 K): all species condensed
        state_extreme_cryo = evaluate_disk_volatile_condensation(
            5.0, P_midplane, mix_cfg, refr_cfg
        )
        @test state_extreme_cryo.condensed_H2O && state_extreme_cryo.condensed_N2

        # Case F: High-temperature disk (T = 1500 K): zero volatile ices condensed
        state_hot = evaluate_disk_volatile_condensation(
            1500.0, P_midplane, mix_cfg, refr_cfg
        )
        @test !state_hot.condensed_H2O &&
            !state_hot.condensed_CO2 &&
            !state_hot.condensed_N2
        @test state_hot.f_refr_C > 0.0
        @test state_hot.f_refr_S > 0.0

        # Case G: Intermediate negative assertions (T = 30 K)
        # CO2 (75 K) and CH4 (45 K) are condensed, but CO (25 K) and N2 (18 K) are NOT condensed
        state_30k = evaluate_disk_volatile_condensation(30.0, P_midplane, mix_cfg, refr_cfg)
        @test state_30k.condensed_CO2 == true
        @test state_30k.condensed_CH4 == true
        @test state_30k.condensed_CO == false
        @test state_30k.condensed_N2 == false

        # Case H: Clausius-Clapeyron pressure shifting (alpha_P > 0)
        st_lo_P = evaluate_disk_volatile_condensation(
            161.0, 0.1, mix_cfg, refr_cfg; P_ref=1.0, alpha_P=0.05
        )
        st_hi_P = evaluate_disk_volatile_condensation(
            161.0, 10.0, mix_cfg, refr_cfg; P_ref=1.0, alpha_P=0.05
        )
        @test st_lo_P.condensed_H2O == false
        @test st_hi_P.condensed_H2O == true

        # Case I: Independent active gating
        st_mix_off = evaluate_disk_volatile_condensation(
            10.0, P_midplane, VolatileMixtureConfig(active=false), RefractoryConfig(active=true)
        )
        @test st_mix_off.condensed_H2O == false
        @test iszero(st_mix_off.X_ice_H2O)
        @test st_mix_off.f_refr_C > 0.0

        st_refr_off = evaluate_disk_volatile_condensation(
            10.0, P_midplane, VolatileMixtureConfig(active=true), RefractoryConfig(active=false)
        )
        @test st_refr_off.condensed_H2O == true
        @test st_refr_off.X_ice_H2O > 0.0
        @test iszero(st_refr_off.f_refr_C)

        # Case J: Pressure clamping at near-vacuum (p_factor strictly positive, floor = 1.6 K for H2O)
        st_vacuum_cold = evaluate_disk_volatile_condensation(
            1.0, 1.0e-12, mix_cfg, refr_cfg; P_ref=1.0, alpha_P=0.5
        )
        @test st_vacuum_cold.condensed_H2O == true
        st_vacuum_warm = evaluate_disk_volatile_condensation(
            10.0, 1.0e-12, mix_cfg, refr_cfg; P_ref=1.0, alpha_P=0.5
        )
        @test st_vacuum_warm.condensed_H2O == false

        # Monotonicity of total condensed ice fraction with decreasing disk temperature
        temps = [300.0, 150.0, 100.0, 70.0, 40.0, 20.0, 10.0]
        ice_totals = [
            let st = evaluate_disk_volatile_condensation(T, P_midplane, mix_cfg, refr_cfg)
                st.X_ice_H2O +
                st.X_ice_NH3 +
                st.X_ice_CO2 +
                st.X_ice_H2S +
                st.X_ice_CH4 +
                st.X_ice_CO +
                st.X_ice_N2
            end for T in temps
        ]
        for k in 1:(length(ice_totals) - 1)
            @test ice_totals[k] <= ice_totals[k + 1]
        end

        # Domain error contracts
        @test_throws DomainError evaluate_disk_volatile_condensation(
            -5.0, P_midplane, mix_cfg, refr_cfg
        )
        @test_throws DomainError evaluate_disk_volatile_condensation(
            100.0, -1.0, mix_cfg, refr_cfg
        )
        @test_throws DomainError evaluate_disk_volatile_condensation(
            100.0, 1.0, mix_cfg, refr_cfg; P_ref=-1.0
        )
        @test_throws DomainError evaluate_disk_volatile_condensation(
            100.0, 1.0, mix_cfg, refr_cfg; alpha_P=-0.1
        )
    end

    # =========================================================================
    # 5. Refractory Pyrolysis & Thermal Devolatilization
    # =========================================================================
    @testset "Refractory IOM Pyrolysis & Breakdown" begin
        refr_cfg = RefractoryConfig(active=true, T_pyrolysis_C=600.0, T_dehydrate_H=750.0)
        C_refr_init = 1000.0   # 1000 ppm refractory carbon
        N_refr_init = 100.0    # 100 ppm refractory nitrogen
        H_refr_init = 50.0     # 50 ppm refractory hydrogen

        # Low temperature (T = 400 K < T_pyrolysis): no pyrolysis occurs
        pyro_low = evaluate_refractory_pyrolysis(
            400.0, C_refr_init, N_refr_init, H_refr_init, refr_cfg
        )
        @test isapprox(pyro_low.C_graphite_residue, 0.0, atol=1e-10)
        @test isapprox(pyro_low.C_refr_remaining, C_refr_init, atol=1e-10)
        @test isapprox(pyro_low.C_devolatilized_gas, 0.0, atol=1e-10)
        @test isapprox(pyro_low.N_devolatilized_gas, 0.0, atol=1e-10)
        @test isapprox(pyro_low.H_dehydrated_gas, 0.0, atol=1e-10)
        @test isapprox(pyro_low.H_refr_remaining, H_refr_init, atol=1e-10)

        # Pyrolysis equality boundary (T = 600.0 K == T_pyrolysis_C): exact boundary
        pyro_boundary = evaluate_refractory_pyrolysis(
            600.0, C_refr_init, N_refr_init, H_refr_init, refr_cfg
        )
        @test isapprox(pyro_boundary.C_refr_remaining, C_refr_init, atol=1e-10)
        @test isapprox(pyro_boundary.C_graphite_residue, 0.0, atol=1e-10)
        @test isapprox(pyro_boundary.C_devolatilized_gas, 0.0, atol=1e-10)

        # High temperature (T = 800 K): partial pyrolysis and onset of dehydration
        pyro_mid = evaluate_refractory_pyrolysis(
            800.0, C_refr_init, N_refr_init, H_refr_init, refr_cfg
        )
        @test pyro_mid.C_refr_remaining < C_refr_init
        @test pyro_mid.C_graphite_residue > 0.0
        @test pyro_mid.C_devolatilized_gas > 0.0
        @test pyro_mid.N_devolatilized_gas > 0.0
        @test pyro_mid.H_dehydrated_gas > 0.0
        @test pyro_mid.H_refr_remaining < H_refr_init
        @test isapprox(
            H_refr_init, pyro_mid.H_refr_remaining + pyro_mid.H_dehydrated_gas, rtol=1e-12
        )

        # Pyrolysis cap saturation (T = 900 K and T = 1200 K): capped at 80%
        pyro_cap = evaluate_refractory_pyrolysis(
            1200.0, C_refr_init, N_refr_init, H_refr_init, refr_cfg
        )
        @test isapprox(pyro_cap.C_refr_remaining, 0.20 * C_refr_init, rtol=1e-10)
        @test isapprox(pyro_cap.C_graphite_residue, 0.60 * 0.80 * C_refr_init, rtol=1e-10)
        @test isapprox(pyro_cap.C_devolatilized_gas, 0.40 * 0.80 * C_refr_init, rtol=1e-10)
        @test isapprox(pyro_cap.H_refr_remaining, 0.0, atol=1e-10)

        # Exact mass conservation: C_init == C_remaining + C_graphite + C_devol
        @test isapprox(
            C_refr_init,
            pyro_mid.C_refr_remaining +
            pyro_mid.C_graphite_residue +
            pyro_mid.C_devolatilized_gas,
            rtol=1e-12,
        )
        @test isapprox(
            N_refr_init,
            pyro_mid.N_refr_remaining + pyro_mid.N_devolatilized_gas,
            rtol=1e-12,
        )

        # Domain error contracts
        @test_throws DomainError evaluate_refractory_pyrolysis(
            -10.0, C_refr_init, N_refr_init, refr_cfg
        )
        @test_throws DomainError evaluate_refractory_pyrolysis(
            500.0, -1.0, N_refr_init, refr_cfg
        )
        @test_throws DomainError evaluate_refractory_pyrolysis(
            500.0, C_refr_init, -10.0, refr_cfg
        )
        @test_throws DomainError evaluate_refractory_pyrolysis(
            500.0, C_refr_init, N_refr_init, -5.0, refr_cfg
        )
        @test_throws DomainError evaluate_refractory_pyrolysis(
            700.0, C_refr_init, N_refr_init; f_graphite=-0.1
        )
        @test_throws DomainError evaluate_refractory_pyrolysis(
            700.0, C_refr_init, N_refr_init; f_graphite=1.1
        )
        @test_throws DomainError evaluate_refractory_pyrolysis(
            700.0, C_refr_init, N_refr_init; DeltaT_pyro=0.0
        )
        @test_throws DomainError evaluate_refractory_pyrolysis(
            700.0, C_refr_init, N_refr_init; DeltaT_dehydrate=-1.0
        )
    end

    # =========================================================================
    # 6. Marker Allocation & Accretion Injection
    # =========================================================================
    @testset "Marker HCNSPO Setup & Accretion Conversion" begin
        marknum = 100
        mix_cfg = VolatileMixtureConfig(active=true)
        refr_cfg = RefractoryConfig(active=true)

        # Setup arrays
        hcnspo_props = setup_marker_hcnspo_properties(marknum, mix_cfg, refr_cfg)
        @test length(hcnspo_props.X_ice_H2O_m) == marknum
        @test length(hcnspo_props.X_ice_NH3_m) == marknum
        @test length(hcnspo_props.X_ice_CO2_m) == marknum
        @test length(hcnspo_props.X_refr_C_m) == marknum
        @test length(hcnspo_props.X_refr_S_m) == marknum
        @test length(hcnspo_props.X_refr_N_m) == marknum
        @test length(hcnspo_props.X_refr_P_m) == marknum
        @test length(hcnspo_props.X_refr_H_m) == marknum

        # Verify initial values are bounded and finite
        @test all(isfinite, hcnspo_props.X_ice_H2O_m)
        @test all(isfinite, hcnspo_props.X_refr_C_m)
        @test all(>=(0.0), hcnspo_props.X_refr_S_m)

        # Accretion injection test
        # Convert sticky air markers to rock inside expanded accretion radius with cold disk state (T = 50 K)
        xm = collect(range(60000.0, 80000.0, length=marknum))
        ym = fill(70000.0, marknum)
        tm = fill(3, marknum)   # start as sticky air
        tkm = fill(150.0, marknum)
        phim = fill(0.40, marknum)
        XWsolidm0 = zeros(Float64, marknum)

        R_curr = 5000.0
        delta_R = 3000.0
        disk_state = evaluate_disk_volatile_condensation(50.0, 10.0, mix_cfg, refr_cfg)

        n_conv = advance_accretion_boundary_hcnspo!(
            R_curr,
            delta_R,
            xm,
            ym,
            tm,
            tkm,
            phim,
            XWsolidm0,
            hcnspo_props,
            disk_state;
            xcenter=70000.0,
            ycenter=70000.0,
        )

        @test n_conv > 0
        # Converted markers must have tm == 2 (rock) and contain accreted ice & refractory inventories
        for m in 1:marknum
            r = sqrt((xm[m] - 70000.0)^2 + (ym[m] - 70000.0)^2)
            if r <= R_curr + delta_R
                @test tm[m] == 2
                @test hcnspo_props.X_ice_H2O_m[m] > 0.0
                @test hcnspo_props.X_ice_CO2_m[m] > 0.0
                @test hcnspo_props.X_refr_C_m[m] > 0.0
                @test hcnspo_props.X_refr_S_m[m] > 0.0
            end
        end

        # Direct test of production advance_accretion_boundary! with hcnspo_props & disk_state
        tm_direct = fill(3, marknum)
        tkm_direct = fill(150.0, marknum)
        phim_direct = fill(0.40, marknum)
        XWsolidm0_direct = zeros(Float64, marknum)
        hcnspo_direct = setup_marker_hcnspo_properties(marknum, mix_cfg, refr_cfg)

        n_conv_direct = advance_accretion_boundary!(
            R_curr,
            delta_R,
            xm,
            ym,
            tm_direct,
            tkm_direct,
            phim_direct,
            XWsolidm0_direct;
            xcenter=70000.0,
            ycenter=70000.0,
            T_accreted=150.0,
            phi_accreted=0.35,
            XWsolid_accreted=disk_state.X_ice_H2O,
            hcnspo_props=hcnspo_direct,
            disk_state=disk_state,
        )
        @test n_conv_direct == n_conv
        for m in 1:marknum
            r = sqrt((xm[m] - 70000.0)^2 + (ym[m] - 70000.0)^2)
            if r <= R_curr + delta_R
                @test tm_direct[m] == 2
                @test isapprox(hcnspo_direct.X_ice_H2O_m[m], disk_state.X_ice_H2O, atol=1e-10)
                @test isapprox(hcnspo_direct.X_refr_C_m[m], disk_state.f_refr_C, atol=1e-10)
            end
        end
    end

    # =========================================================================
    # 7. TOML Configuration Serialization & Round-Trip
    # =========================================================================
    @testset "Configuration Schema & TOML Round-Trip" begin
        refr_cfg = RefractoryConfig(
            active=true,
            f_refr_C=0.72,
            f_refr_N=0.15,
            f_refr_S=0.91,
            f_refr_P=0.95,
            f_refr_H=0.08,
            T_pyrolysis_C=620.0,
            T_dehydrate_H=760.0,
        )
        vol_cfg = VolatileMixtureConfig(
            active=true,
            X_ice_H2O=0.80,
            X_ice_CO2=0.10,
            X_ice_CH4=0.03,
            X_ice_NH3=0.04,
            lambda_nh3_depression=280.0,
            T_freeze_floor=178.0,
        )
        sim_cfg = SimulationConfig(refractory=refr_cfg, volatile_mixture=vol_cfg)
        validate_config(sim_cfg)

        mktempdir() do tmpdir
            toml_path = joinpath(tmpdir, "test_config.toml")
            save_config(toml_path, sim_cfg)
            @test isfile(toml_path)

            loaded_cfg = load_config(toml_path)
            validate_config(loaded_cfg)

            @test loaded_cfg.refractory.active == true
            @test isapprox(loaded_cfg.refractory.f_refr_C, 0.72)
            @test isapprox(loaded_cfg.refractory.f_refr_N, 0.15)
            @test isapprox(loaded_cfg.refractory.f_refr_S, 0.91)
            @test isapprox(loaded_cfg.refractory.f_refr_P, 0.95)
            @test isapprox(loaded_cfg.refractory.f_refr_H, 0.08)
            @test isapprox(loaded_cfg.refractory.T_pyrolysis_C, 620.0)
            @test isapprox(loaded_cfg.refractory.T_dehydrate_H, 760.0)

            @test loaded_cfg.volatile_mixture.active == true
            @test isapprox(loaded_cfg.volatile_mixture.X_ice_H2O, 0.80)
            @test isapprox(loaded_cfg.volatile_mixture.X_ice_CO2, 0.10)
            @test isapprox(loaded_cfg.volatile_mixture.X_ice_CH4, 0.03)
            @test isapprox(loaded_cfg.volatile_mixture.X_ice_NH3, 0.04)
            @test isapprox(loaded_cfg.volatile_mixture.lambda_nh3_depression, 280.0)
            @test isapprox(loaded_cfg.volatile_mixture.T_freeze_floor, 178.0)
        end
    end
end

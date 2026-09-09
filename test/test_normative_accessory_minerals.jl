using Erebus
using Test
using JLD2

@testset "Normative Accessory Minerals & Meteorite Diagnostics" begin
    @testset "PhaseTrackingConfig Schema & Validation" begin
        # Default config
        cfg = PhaseTrackingConfig()
        @test cfg.active == false
        @test cfg.T_eutectic ≈ 1213.0
        @test cfg.dT_transition ≈ 50.0
        @test cfg.bulk_P_ppm ≈ 1000.0
        @test cfg.schreibersite_ni_frac ≈ 0.25
        @test cfg.cohenite_carbide_max ≈ 0.0667
        @test cfg.nitride_mode === :roaldite
        @test cfg.track_regional_modes == true
        @test cfg.r_core_norm ≈ 0.5
        @test cfg.r_mantle_norm ≈ 0.85

        # Custom config
        cfg_custom = PhaseTrackingConfig(
            active=true,
            T_eutectic=1260.0,
            dT_transition=30.0,
            bulk_P_ppm=1500.0,
            schreibersite_ni_frac=0.30,
            cohenite_carbide_max=0.06,
            nitride_mode=:carlsbergite,
            track_regional_modes=true,
            r_core_norm=0.45,
            r_mantle_norm=0.80,
        )
        @test cfg_custom.active == true
        @test cfg_custom.T_eutectic ≈ 1260.0
        @test cfg_custom.nitride_mode === :carlsbergite

        # Validation passes for valid configs
        sim_cfg = SimulationConfig(phase_tracking=cfg_custom)
        @test validate_config(sim_cfg) === nothing

        # Validation errors on unphysical inputs
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                phase_tracking=PhaseTrackingConfig(active=true, T_eutectic=-10.0)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                phase_tracking=PhaseTrackingConfig(active=true, dT_transition=-5.0)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                phase_tracking=PhaseTrackingConfig(active=true, bulk_P_ppm=-100.0)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                phase_tracking=PhaseTrackingConfig(active=true, schreibersite_ni_frac=-0.1)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                phase_tracking=PhaseTrackingConfig(active=true, schreibersite_ni_frac=1.2)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                phase_tracking=PhaseTrackingConfig(active=true, cohenite_carbide_max=-0.01)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                phase_tracking=PhaseTrackingConfig(active=true, cohenite_carbide_max=0.0)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                phase_tracking=PhaseTrackingConfig(active=true, cohenite_carbide_max=1.5)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                phase_tracking=PhaseTrackingConfig(
                    active=true, nitride_mode=:invalid_nitride
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                phase_tracking=PhaseTrackingConfig(
                    active=true, r_core_norm=0.9, r_mantle_norm=0.5
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                phase_tracking=PhaseTrackingConfig(
                    active=true, r_core_norm=0.5, r_mantle_norm=0.5
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                phase_tracking=PhaseTrackingConfig(active=true, r_core_norm=-0.1)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                phase_tracking=PhaseTrackingConfig(active=true, r_core_norm=0.0)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                phase_tracking=PhaseTrackingConfig(active=true, r_core_norm=1.0)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                phase_tracking=PhaseTrackingConfig(active=true, r_mantle_norm=1.2)
            ),
        )
    end

    @testset "TOML Round-Trip Serialization" begin
        cfg_orig = SimulationConfig(
            phase_tracking=PhaseTrackingConfig(
                active=true,
                T_eutectic=1220.0,
                dT_transition=40.0,
                bulk_P_ppm=1200.0,
                schreibersite_ni_frac=0.28,
                cohenite_carbide_max=0.065,
                nitride_mode=:osbornite,
                track_regional_modes=true,
                r_core_norm=0.48,
                r_mantle_norm=0.82,
            ),
        )

        toml_str = serialize_config(cfg_orig)
        @test occursin("[phase_tracking]", toml_str)
        @test occursin("active = true", toml_str)
        @test occursin("T_eutectic = 1220.0", toml_str)
        @test occursin("nitride_mode = \"osbornite\"", toml_str)

        cfg_parsed = parse_config_string(toml_str)
        @test cfg_parsed.phase_tracking.active == true
        @test cfg_parsed.phase_tracking.T_eutectic ≈ 1220.0
        @test cfg_parsed.phase_tracking.dT_transition ≈ 40.0
        @test cfg_parsed.phase_tracking.bulk_P_ppm ≈ 1200.0
        @test cfg_parsed.phase_tracking.schreibersite_ni_frac ≈ 0.28
        @test cfg_parsed.phase_tracking.cohenite_carbide_max ≈ 0.065
        @test cfg_parsed.phase_tracking.nitride_mode === :osbornite
        @test cfg_parsed.phase_tracking.r_core_norm ≈ 0.48
        @test cfg_parsed.phase_tracking.r_mantle_norm ≈ 0.82
    end

    @testset "Stoichiometric Mineral Conversion" begin
        # 1. Troilite (FeS)
        # S molar mass: 32.065, Fe: 55.845, FeS: 87.910 => factor ~ 2.7416
        w_S = 0.05 # 5 wt% S in metallic alloy
        w_troilite, w_fe_consumed_S = compute_troilite_stoichiometry(w_S)
        @test w_troilite ≈ w_S * (87.910 / 32.065) atol=1e-5
        @test w_fe_consumed_S ≈ w_S * (55.845 / 32.065) atol=1e-5
        @test w_troilite ≈ w_S + w_fe_consumed_S atol=1e-12

        # 2. Schreibersite ((Fe,Ni)3P)
        # P molar mass: 30.97376, Ni frac: 0.25 => 3 * (0.75*55.845 + 0.25*58.6934) + 30.97376
        w_P = 0.002 # 2000 ppmw P
        w_schreib, w_met_consumed_P = compute_schreibersite_stoichiometry(w_P; ni_frac=0.25)
        M_metal_avg = 0.75 * 55.845 + 0.25 * 58.6934
        M_schreib_calc = 3.0 * M_metal_avg + 30.97376
        f_schreib = M_schreib_calc / 30.97376
        @test w_schreib ≈ w_P * f_schreib atol=1e-6
        @test w_schreib ≈ w_P + w_met_consumed_P atol=1e-12

        # 3. Cohenite ((Fe,Ni)3C) and Graphite (C)
        # Carbon below saturation: all in cohenite
        w_C_low = 0.001 # 1000 ppmw C
        w_cohenite, w_graphite, w_fe_consumed_C = compute_cohenite_graphite_stoichiometry(
            w_C_low; carbide_max=0.0667
        )
        f_cohenite = (3.0 * 55.845 + 12.011) / 12.011 # ~ 14.948
        @test w_cohenite ≈ w_C_low * f_cohenite atol=1e-5
        @test iszero(w_graphite)
        @test w_cohenite ≈ w_C_low + w_fe_consumed_C atol=1e-12

        # Carbon above saturation: cohenite capped at carbide_max * f_cohenite, excess in graphite
        w_C_high = 0.08 # 8 wt% C (exceeds stoichiometric cohenite)
        w_coh_high, w_gra_high, w_fe_high = compute_cohenite_graphite_stoichiometry(
            w_C_high; carbide_max=0.0667
        )
        @test w_coh_high ≈ 0.0667 * f_cohenite atol=1e-5
        @test w_gra_high ≈ (w_C_high - 0.0667) atol=1e-12
        @test w_fe_high ≈ 0.0667 * (f_cohenite - 1.0) atol=1e-5

        # Custom carbide_max (e.g. 0.05)
        w_coh_cust, w_gra_cust, w_fe_cust = compute_cohenite_graphite_stoichiometry(
            w_C_high; carbide_max=0.05
        )
        @test w_coh_cust ≈ 0.05 * f_cohenite atol=1e-5
        @test w_gra_cust ≈ (w_C_high - 0.05) atol=1e-12
        @test w_fe_cust ≈ 0.05 * (f_cohenite - 1.0) atol=1e-5

        # 4. Nitrides (Roaldite, Carlsbergite, Osbornite)
        w_N = 0.0001 # 100 ppmw N
        w_roaldite, _ = compute_nitride_stoichiometry(w_N; mode=:roaldite)
        f_roaldite = (4.0 * 55.845 + 14.007) / 14.007 # ~ 16.948
        @test w_roaldite ≈ w_N * f_roaldite atol=1e-5

        w_carlsbergite, _ = compute_nitride_stoichiometry(w_N; mode=:carlsbergite)
        f_carlsbergite = (51.996 + 14.007) / 14.007 # ~ 4.712
        @test w_carlsbergite ≈ w_N * f_carlsbergite atol=1e-5

        w_osbornite, _ = compute_nitride_stoichiometry(w_N; mode=:osbornite)
        f_osbornite = (47.867 + 14.007) / 14.007 # ~ 4.417
        @test w_osbornite ≈ w_N * f_osbornite atol=1e-5

        # Domain error on negative inputs
        @test_throws DomainError compute_troilite_stoichiometry(-0.01)
        @test_throws DomainError compute_schreibersite_stoichiometry(-0.01)
        @test_throws DomainError compute_cohenite_graphite_stoichiometry(-0.01)
        @test_throws DomainError compute_nitride_stoichiometry(-0.01)
    end

    @testset "Thermal Eutectic Melting Transition" begin
        cfg = PhaseTrackingConfig(
            active=true, T_eutectic=1213.0, dT_transition=50.0, bulk_P_ppm=1000.0
        )

        w_S = 0.04
        w_C = 0.002
        w_N = 0.0002
        w_P = 0.001

        # 1. Sub-eutectic state (T = 1000 K < T_eutectic): 100% solid mineral assemblage
        res_sub = compute_normative_mineral_assemblage(1000.0, w_S, w_C, w_N, w_P, cfg)
        @test res_sub.F_solid ≈ 1.0
        @test res_sub.F_liquid ≈ 0.0
        @test res_sub.w_liquid_alloy ≈ 0.0
        @test res_sub.w_troilite > 0.0
        @test res_sub.w_schreibersite > 0.0
        @test res_sub.w_cohenite > 0.0
        @test res_sub.w_metal_matrix > 0.0
        # Mass conservation: all phases sum to 1.0 of the metallic system
        total_sub = (
            res_sub.w_troilite +
            res_sub.w_schreibersite +
            res_sub.w_cohenite +
            res_sub.w_graphite +
            res_sub.w_nitride +
            res_sub.w_metal_matrix +
            res_sub.w_liquid_alloy
        )
        @test total_sub ≈ 1.0 atol=1e-12

        # 2. Super-eutectic fully molten state (T = 1400 K > T_eutectic + dT_transition)
        res_super = compute_normative_mineral_assemblage(1400.0, w_S, w_C, w_N, w_P, cfg)
        @test res_super.F_solid ≈ 0.0
        @test res_super.F_liquid ≈ 1.0
        @test res_super.w_troilite ≈ 0.0
        @test res_super.w_schreibersite ≈ 0.0
        @test res_super.w_cohenite ≈ 0.0
        @test res_super.w_graphite ≈ 0.0
        @test res_super.w_nitride ≈ 0.0
        @test res_super.w_metal_matrix ≈ 0.0
        @test res_super.w_liquid_alloy ≈ 1.0

        # 3. Partial melting transition state (T = 1238 K, midpoint: F_solid = 0.5)
        T_mid = 1213.0 + 25.0
        res_mid = compute_normative_mineral_assemblage(T_mid, w_S, w_C, w_N, w_P, cfg)
        @test res_mid.F_solid ≈ 0.5 atol=1e-6
        @test res_mid.F_liquid ≈ 0.5 atol=1e-6
        @test res_mid.w_liquid_alloy ≈ 0.5 atol=1e-6
        @test res_mid.w_troilite ≈ 0.5 * res_sub.w_troilite atol=1e-6
        @test res_mid.w_schreibersite ≈ 0.5 * res_sub.w_schreibersite atol=1e-6
        @test res_mid.w_cohenite ≈ 0.5 * res_sub.w_cohenite atol=1e-6

        total_mid = (
            res_mid.w_troilite +
            res_mid.w_schreibersite +
            res_mid.w_cohenite +
            res_mid.w_graphite +
            res_mid.w_nitride +
            res_mid.w_metal_matrix +
            res_mid.w_liquid_alloy
        )
        @test total_mid ≈ 1.0 atol=1e-12

        # High sulfur (w_S = 0.50): troilite capped at 1.0, assemblage strictly conserved
        res_high_S = compute_normative_mineral_assemblage(1100.0, 0.50, 0.0, 0.0, 0.0, cfg)
        @test res_high_S.w_troilite ≈ 1.0 atol=1e-12
        @test iszero(res_high_S.w_metal_matrix)
        @test (res_high_S.w_troilite + res_high_S.w_liquid_alloy) ≈ 1.0 atol=1e-12

        # Multiple high volatiles (w_S=0.30, w_C=0.05, w_P=0.01): total solid strictly <= 1.0
        res_multi_high = compute_normative_mineral_assemblage(
            1100.0, 0.30, 0.05, 0.001, 0.01, cfg
        )
        total_solid_multi = (
            res_multi_high.w_troilite +
            res_multi_high.w_schreibersite +
            res_multi_high.w_cohenite +
            res_multi_high.w_graphite +
            res_multi_high.w_nitride +
            res_multi_high.w_metal_matrix
        )
        @test total_solid_multi ≈ res_multi_high.F_solid atol=1e-12
        @test (total_solid_multi + res_multi_high.w_liquid_alloy) ≈ 1.0 atol=1e-12

        # Strict input domain checking
        @test_throws DomainError compute_normative_mineral_assemblage(
            -10.0, w_S, w_C, w_N, w_P, cfg
        )
        @test_throws DomainError compute_normative_mineral_assemblage(
            1000.0, -0.01, w_C, w_N, w_P, cfg
        )
        @test_throws DomainError compute_normative_mineral_assemblage(
            1000.0, w_S, -0.01, w_N, w_P, cfg
        )
        @test_throws DomainError compute_normative_mineral_assemblage(
            1000.0, w_S, w_C, -0.001, w_P, cfg
        )
        @test_throws DomainError compute_normative_mineral_assemblage(
            1000.0, w_S, w_C, w_N, -0.001, cfg
        )
        @test_throws DomainError compute_normative_mineral_assemblage(
            1000.0, w_S, w_C, w_N, w_P, PhaseTrackingConfig(dT_transition=0.0)
        )
        @test_throws DomainError compute_normative_mineral_assemblage(
            1000.0, w_S, w_C, w_N, w_P, PhaseTrackingConfig(dT_transition=-10.0)
        )
    end

    @testset "Marker Setup, Property Computation & Replenishment" begin
        cfg = PhaseTrackingConfig(active=true, T_eutectic=1213.0)
        marknum = 100
        arrays = setup_marker_phase_tracking_properties(marknum, cfg)
        @test length(arrays.Xmin_troilite_m) == marknum
        @test length(arrays.Xmin_schreibersite_m) == marknum
        @test length(arrays.Xmin_cohenite_m) == marknum
        @test length(arrays.Xmin_graphite_m) == marknum
        @test length(arrays.Xmin_nitride_m) == marknum
        @test length(arrays.Xmin_metal_matrix_m) == marknum

        # Inactive returns nothing
        cfg_off = PhaseTrackingConfig(active=false)
        arrays_off = setup_marker_phase_tracking_properties(marknum, cfg_off)
        @test arrays_off.Xmin_troilite_m === nothing
    end

    @testset "Regional Mineral Modes & Meteorite Classification" begin
        # 100 markers across a 50 km radius planetesimal
        marknum = 100
        xm = collect(range(70000.0 - 45000.0, 70000.0 + 45000.0, length=marknum))
        ym = fill(70000.0, marknum)
        tm = fill(1, marknum) # rock/core markers
        Xfe_bulk = fill(0.20, marknum)
        # Inner markers (r < 25 km): warm core/mantle (T = 1300 K)
        # Outer markers (r >= 25 km): cold mantle/crust (T = 800 K)
        tkm = zeros(marknum)
        for m in 1:marknum
            r = abs(xm[m] - 70000.0)
            tkm[m] = r < 25000.0 ? 1350.0 : 850.0
        end

        Xfe_S_m = fill(40000.0, marknum) # 4 wt% S
        Xfe_C_m = fill(2000.0, marknum)  # 2000 ppmw C
        Xfe_N_m = fill(100.0, marknum)   # 100 ppmw N

        cfg = PhaseTrackingConfig(
            active=true, T_eutectic=1213.0, r_core_norm=0.5, r_mantle_norm=0.85
        )

        modes = compute_regional_mineral_modes(
            xm,
            ym,
            tm,
            tkm,
            Xfe_bulk,
            Xfe_S_m,
            Xfe_C_m,
            Xfe_N_m,
            marknum;
            cfg=cfg,
            rplanet=50000.0,
            xcenter=70000.0,
            ycenter=70000.0,
        )

        @test modes.M_total_troilite >= 0.0
        @test modes.M_total_schreibersite >= 0.0
        @test modes.M_total_cohenite >= 0.0
        @test modes.M_crust_troilite > 0.0 # cold crust retains solid troilite
        @test modes.M_core_liquid_alloy > 0.0 # warm core has molten alloy
        @test modes.classification in
            (:IAB_winonaite_primitive, :magmatic_differentiated, :transitional)

        # Test marknum = 0 edge case
        empty_modes = compute_regional_mineral_modes(
            Float64[],
            Float64[],
            Int[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            0;
            cfg=cfg,
        )
        @test iszero(empty_modes.M_total_metal)
        @test empty_modes.classification === :transitional

        # Test invalid rho_metal throws DomainError
        @test_throws DomainError compute_regional_mineral_modes(
            xm,
            ym,
            tm,
            tkm,
            Xfe_bulk,
            Xfe_S_m,
            Xfe_C_m,
            Xfe_N_m,
            marknum;
            cfg=cfg,
            rho_metal=-100.0,
        )
        @test_throws DomainError compute_regional_mineral_modes(
            xm,
            ym,
            tm,
            tkm,
            Xfe_bulk,
            Xfe_S_m,
            Xfe_C_m,
            Xfe_N_m,
            marknum;
            cfg=cfg,
            rho_metal=0.0,
        )
    end
end

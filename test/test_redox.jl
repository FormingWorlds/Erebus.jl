using Test
using Erebus

@testset "Redox Buffer Coefficients and Equilibrium Physics" begin
    # Test conditions: magmatic temperature 1400 K, lithostatic pressure 1e8 Pa (1 kbar)
    T = 1400.0
    P_1bar = 1.0e5
    P_1kbar = 1.0e8

    # Primary pin: Iron-Wüstite must match existing compute_iron_wustite_fO2 identically
    lfo2_iw_engine = Erebus.log10_fo2_of_buffer(:IW, T, P_1bar)
    lfo2_iw_legacy = Erebus.compute_iron_wustite_fO2(T; delta_IW=0.0)
    @test isapprox(lfo2_iw_engine, lfo2_iw_legacy; atol=1e-12)

    # Frost (1991) Table 1 Iron-Wüstite calibration pin
    lfo2_iw_frost = Erebus.log10_fo2_of_buffer(:IW_Frost, T, P_1bar)
    @test isapprox(lfo2_iw_frost, -27489.0 / T + 6.702; atol=1e-6)
    # Quantified offset: Frost IW sits ~0.64 dex above Campbell IW at 1400 K
    @test isapprox(lfo2_iw_frost - lfo2_iw_engine, 0.6433; atol=0.01)

    # 1. Petrologic ordering guard at 1400 K and 1 bar (Frost 1991):
    # MH > NNO > QFM > WM > IW > QIF
    lfo2_mh = Erebus.log10_fo2_of_buffer(:MH, T, P_1bar)
    lfo2_nno = Erebus.log10_fo2_of_buffer(:NNO, T, P_1bar)
    lfo2_qfm = Erebus.log10_fo2_of_buffer(:QFM, T, P_1bar)
    lfo2_wm = Erebus.log10_fo2_of_buffer(:WM, T, P_1bar)
    lfo2_iw = Erebus.log10_fo2_of_buffer(:IW, T, P_1bar)
    lfo2_qif = Erebus.log10_fo2_of_buffer(:QIF, T, P_1bar)

    @test lfo2_mh > lfo2_nno
    @test lfo2_nno > lfo2_qfm
    @test lfo2_qfm > lfo2_wm
    @test lfo2_wm > lfo2_iw
    @test lfo2_iw > lfo2_qif

    # 2. 3-class discrimination guards on QFM at 1400 K (Frost 1991):
    # Expected: -25096.3 / 1400 + 8.735 = -9.1909
    # Scale guard: -10.0 < log10_fO2 < -8.0
    @test -10.0 < lfo2_qfm < -8.0
    # Sign guard: oxygen fugacity in bar is far below atmospheric (negative in log10)
    @test lfo2_qfm < 0.0
    # Exponent/prefactor guard: delta relative to IW is ~ +4.4 log units
    delta_qfm_iw = lfo2_qfm - lfo2_iw
    @test isapprox(delta_qfm_iw, 4.3852; atol=0.05)
    @test abs(delta_qfm_iw - 8.0) > 2.0

    # 3. Pressure sensitivity guard (C > 0 increases log10_fO2 with pressure):
    lfo2_qfm_highp = Erebus.log10_fo2_of_buffer(:QFM, T, P_1kbar)
    @test lfo2_qfm_highp > lfo2_qfm
    # Delta at 1 kbar (999 bar offset): C * (P_bar - 1) / T = 0.110 * 999 / 1400 ≈ 0.0785 log units
    @test isapprox(lfo2_qfm_highp - lfo2_qfm, 0.110 * 999.0 / 1400.0; atol=1e-4)

    # 4. Input error contracts
    @test_throws DomainError Erebus.log10_fo2_of_buffer(:IW, -100.0, P_1bar)
    @test_throws DomainError Erebus.log10_fo2_of_buffer(:IW, 0.0, P_1bar)
    @test_throws DomainError Erebus.log10_fo2_of_buffer(:IW, NaN, P_1bar)
    @test_throws DomainError Erebus.log10_fo2_of_buffer(:IW, T, -1.0)
    @test_throws DomainError Erebus.log10_fo2_of_buffer(:IW, T, NaN)
    @test_throws ArgumentError Erebus.log10_fo2_of_buffer(:UNKNOWN_BUFFER, T, P_1bar)
end

@testset "Graphite CCO Buffer Inversion Physics" begin
    # Test CCO buffer across pressures at T = 1300 K (French 1966)
    T = 1300.0
    P_10bar = 1.0e6   # 10 bar in Pa
    P_100bar = 1.0e7  # 100 bar in Pa

    lfo2_cco_10 = Erebus.log10_fo2_of_buffer(:CCO, T, P_10bar)
    lfo2_cco_100 = Erebus.log10_fo2_of_buffer(:CCO, T, P_100bar)

    # Primary pin: Inverting fO2 back into compute_graphite_saturation_fugacity
    # must recover the total gas pressure f_CO + f_CO2 = P_bar
    sat_10 = Erebus.compute_graphite_saturation_fugacity(T, lfo2_cco_10)
    p_tot_10 = sat_10.f_CO_max_bar + sat_10.f_CO2_max_bar
    @test isapprox(p_tot_10, 10.0; rtol=1e-6)

    sat_100 = Erebus.compute_graphite_saturation_fugacity(T, lfo2_cco_100)
    p_tot_100 = sat_100.f_CO_max_bar + sat_100.f_CO2_max_bar
    @test isapprox(p_tot_100, 100.0; rtol=1e-6)

    # Pressure dependence guard: higher pressure shifts equilibrium toward CO2, increasing fO2
    @test lfo2_cco_100 > lfo2_cco_10

    # Scale guard: at 1300 K and 10 bar, CCO sits at ΔIW ≈ -0.92
    lfo2_iw = Erebus.log10_fo2_of_buffer(:IW, T, P_10bar)
    delta_cco_iw = lfo2_cco_10 - lfo2_iw
    @test -1.5 < delta_cco_iw < -0.5
    @test isapprox(delta_cco_iw, -0.92; atol=0.1)

    # Low pressure stability guard: verify no catastrophic cancellation down to 1e-10 Pa
    P_low = 1.0e-10
    lfo2_cco_low = Erebus.log10_fo2_of_buffer(:CCO, T, P_low)
    @test isfinite(lfo2_cco_low)
    @test lfo2_cco_low < lfo2_cco_10

    # Error contracts
    @test_throws DomainError Erebus.log10_fo2_of_buffer(:CCO, 50.0, P_10bar) # T < 100 K
    @test_throws DomainError Erebus.log10_fo2_of_buffer(:CCO, T, 0.0)        # P <= 0
end

@testset "Bidirectional Buffer Translation and Invariants" begin
    T = 1500.0
    P = 5.0e7 # 500 bar

    # 1. Exact round-trip identity: convert(convert(v, A, B), B, A) == v
    for val in (-5.0, -2.5, -1.0, 0.0, 1.5, 3.0)
        rt = Erebus.convert_redox_buffer(
            Erebus.convert_redox_buffer(val, :IW, :QFM, T, P), :QFM, :IW, T, P
        )
        @test isapprox(rt, val; atol=1e-12)

        rt_nno = Erebus.convert_redox_buffer(
            Erebus.convert_redox_buffer(val, :QFM, :NNO, T, P), :NNO, :QFM, T, P
        )
        @test isapprox(rt_nno, val; atol=1e-12)
    end

    # 2. Closed translation triangle: IW -> QFM -> NNO -> IW
    val_init = -1.5
    val_qfm = Erebus.convert_redox_buffer(val_init, :IW, :QFM, T, P)
    val_nno = Erebus.convert_redox_buffer(val_qfm, :QFM, :NNO, T, P)
    val_final = Erebus.convert_redox_buffer(val_nno, :NNO, :IW, T, P)
    @test isapprox(val_final, val_init; atol=1e-12)

    # 3. Consistency with delta_buffer helpers
    lfo2 = Erebus.delta_buffer_to_log10_fo2(-1.0, :IW, T, P)
    delta_qfm = Erebus.log10_fo2_to_delta_buffer(lfo2, :QFM, T, P)
    delta_direct = Erebus.convert_redox_buffer(-1.0, :IW, :QFM, T, P)
    @test isapprox(delta_qfm, delta_direct; atol=1e-12)

    # 4. Error contracts on conversion
    @test_throws DomainError Erebus.convert_redox_buffer(NaN, :IW, :QFM, T, P)
    @test_throws DomainError Erebus.convert_redox_buffer(-1.0, :IW, :QFM, -10.0, P)
    @test_throws ArgumentError Erebus.convert_redox_buffer(-1.0, :INVALID, :QFM, T, P)
    @test_throws ArgumentError Erebus.convert_redox_buffer(-1.0, :IW, :INVALID, T, P)
end

@testset "Local Controlling Buffer Regime Selection" begin
    # Case 1: Metal-saturated (Fe0 present above tolerance) -> :IW
    buf1, reg1 = Erebus.local_controlling_buffer(0.05, 0.0, 0.0)
    @test buf1 === :IW
    @test reg1 === :metal_saturated

    # Case 2: Graphite-saturated (graphite present, no metal) -> :CCO
    buf2, reg2 = Erebus.local_controlling_buffer(0.0, 0.02, 0.0)
    @test buf2 === :CCO
    @test reg2 === :graphite_saturated

    # Case 3: Silicate melt with ferric iron -> :QFM display
    buf3, reg3 = Erebus.local_controlling_buffer(0.0, 0.0, 0.15)
    @test buf3 === :QFM
    @test reg3 === :silicate_melt

    # Case 4: Degassed pore / trace gas fallback -> :QFM
    buf4, reg4 = Erebus.local_controlling_buffer(0.0, 0.0, 0.0)
    @test buf4 === :QFM
    @test reg4 === :gas_ratio_fallback

    # Boundedness: tolerance threshold gating
    buf_tol, _ = Erebus.local_controlling_buffer(1e-7, 0.0, 0.0; tol=1e-6)
    @test buf_tol === :QFM

    # Domain error contracts
    @test_throws DomainError Erebus.local_controlling_buffer(-0.1, 0.0, 0.0)
    @test_throws DomainError Erebus.local_controlling_buffer(1.2, 0.0, 0.0)
    @test_throws DomainError Erebus.local_controlling_buffer(NaN, 0.0, 0.0)
    @test_throws DomainError Erebus.local_controlling_buffer(0.0, -0.05, 0.0)
    @test_throws DomainError Erebus.local_controlling_buffer(0.0, 0.0, 1.05)
    @test_throws DomainError Erebus.local_controlling_buffer(0.0, 0.0, 0.0; tol=-1e-4)
end

@testset "Evans 2012 Redox Budget Electron Accounting" begin
    # Reference: Evans (2012) Earth-Sci. Rev. 113, 11-32, DOI: 10.1016/j.earscirev.2012.03.003
    # Mantle reference state (M): Fe2+, C0, S2-, H+, O2-, P5+
    # Crust reference state (C): Fe3+, C4+, S6+, H+, O2-, P5+

    # 1. Pure component electron counts relative to Mantle reference state
    c_fe0 = Erebus.RedoxComponents(n_Fe0=1.0)
    rb_fe0 = Erebus.compute_redox_budget(c_fe0; reference=:mantle)
    @test isapprox(rb_fe0, -2.0; atol=1e-12)

    c_fe3 = Erebus.RedoxComponents(n_Fe3=1.0)
    rb_fe3 = Erebus.compute_redox_budget(c_fe3; reference=:mantle)
    @test isapprox(rb_fe3, 1.0; atol=1e-12)

    c_h2 = Erebus.RedoxComponents(n_H2=1.0)
    rb_h2 = Erebus.compute_redox_budget(c_h2; reference=:mantle)
    @test isapprox(rb_h2, -2.0; atol=1e-12)

    c_h2o = Erebus.RedoxComponents(n_H2O=1.0)
    rb_h2o = Erebus.compute_redox_budget(c_h2o; reference=:mantle)
    @test isapprox(rb_h2o, 0.0; atol=1e-12)

    c_co2 = Erebus.RedoxComponents(n_CO2=1.0)
    rb_co2 = Erebus.compute_redox_budget(c_co2; reference=:mantle)
    @test isapprox(rb_co2, 4.0; atol=1e-12)

    c_ch4 = Erebus.RedoxComponents(n_CH4=1.0)
    rb_ch4 = Erebus.compute_redox_budget(c_ch4; reference=:mantle)
    @test isapprox(rb_ch4, -4.0; atol=1e-12)

    # Specific redox budget: RB / mass_kg
    mass = 100.0 # kg
    rb_spec = Erebus.compute_specific_redox_budget(c_fe0, mass; reference=:mantle)
    @test isapprox(rb_spec, -0.02; atol=1e-12)

    # 2. Conservation Invariant: Serpentinization Reaction (both zero and non-zero backgrounds)
    c_before = Erebus.RedoxComponents(n_Fe2=3.0, n_H2O=1.0)
    c_after = Erebus.serpentinize_redox_budget(c_before, 1.0)
    rb_init = Erebus.compute_redox_budget(c_before; reference=:mantle)
    rb_final = Erebus.compute_redox_budget(c_after; reference=:mantle)
    @test isapprox(rb_init, rb_final; atol=1e-12)
    @test isapprox(c_after.n_Fe3, 2.0; atol=1e-12)
    @test isapprox(c_after.n_H2, 1.0; atol=1e-12)

    # Non-zero background: ensure conservation holds with coexisting Fe0 and Fe3+
    c_bg = Erebus.RedoxComponents(n_Fe0=10.0, n_Fe2=10.0, n_Fe3=2.0, n_H2O=5.0)
    c_bg_after = Erebus.serpentinize_redox_budget(c_bg, 2.0)
    rb_bg_init = Erebus.compute_redox_budget(c_bg; reference=:mantle)
    rb_bg_final = Erebus.compute_redox_budget(c_bg_after; reference=:mantle)
    @test isapprox(rb_bg_init, rb_bg_final; atol=1e-12)
    @test rb_bg_init != 0.0

    # Over-consumption error contract: cannot react more Fe2+ or H2O than available
    @test_throws DomainError Erebus.serpentinize_redox_budget(c_before, 5.0)
    @test_throws DomainError Erebus.serpentinize_redox_budget(c_before, -0.5)

    # 3. Conservation Invariant: Core Segregation
    c_bulk = Erebus.RedoxComponents(n_Fe0=50.0, n_Fe2=100.0, n_Fe3=5.0)
    c_mantle, c_core = Erebus.segregate_core_redox_budget(c_bulk, 0.8)
    rb_bulk = Erebus.compute_redox_budget(c_bulk; reference=:mantle)
    rb_mantle = Erebus.compute_redox_budget(c_mantle; reference=:mantle)
    rb_core = Erebus.compute_redox_budget(c_core; reference=:mantle)
    @test isapprox(rb_bulk, rb_mantle + rb_core; atol=1e-12)
    @test rb_core < 0.0
    @test rb_mantle > rb_bulk

    # Core segregation domain contract
    @test_throws DomainError Erebus.segregate_core_redox_budget(c_bulk, -0.1)
    @test_throws DomainError Erebus.segregate_core_redox_budget(c_bulk, 1.1)

    # 4. Conservation Invariant: Degassing and Gas Venting
    c_rock0 = Erebus.RedoxComponents(n_Fe2=100.0, n_H2=10.0, n_CO=5.0, n_CO2=5.0)
    c_vent = Erebus.RedoxComponents(n_H2=8.0, n_CO=3.0)
    c_rock1 = Erebus.vent_gas_redox_budget(c_rock0, c_vent)
    rb_rock0 = Erebus.compute_redox_budget(c_rock0; reference=:mantle)
    rb_rock1 = Erebus.compute_redox_budget(c_rock1; reference=:mantle)
    rb_vent = Erebus.compute_redox_budget(c_vent; reference=:mantle)
    @test isapprox(rb_rock0, rb_rock1 + rb_vent; atol=1e-12)
    @test rb_vent < 0.0
    @test rb_rock1 > rb_rock0

    # Over-venting domain contract: cannot vent more than available in rock
    c_over_vent = Erebus.RedoxComponents(n_H2=15.0)
    @test_throws DomainError Erebus.vent_gas_redox_budget(c_rock0, c_over_vent)

    # Non-volatile vent domain contract: cannot vent condensed phases
    c_nonvol_vent = Erebus.RedoxComponents(n_Fe0=1.0)
    @test_throws DomainError Erebus.vent_gas_redox_budget(c_rock0, c_nonvol_vent)

    # 5. Crust reference state conversion invariance and reaction conservation
    rb_fe0_crust = Erebus.compute_redox_budget(c_fe0; reference=:crust)
    @test isapprox(rb_fe0_crust, -3.0; atol=1e-12)

    # Crust reference conservation across serpentinization
    rb_crust_init = Erebus.compute_redox_budget(c_before; reference=:crust)
    rb_crust_final = Erebus.compute_redox_budget(c_after; reference=:crust)
    @test isapprox(rb_crust_init, rb_crust_final; atol=1e-12)

    # Crust reference conservation across core segregation
    rb_crust_bulk = Erebus.compute_redox_budget(c_bulk; reference=:crust)
    rb_crust_mantle = Erebus.compute_redox_budget(c_mantle; reference=:crust)
    rb_crust_core = Erebus.compute_redox_budget(c_core; reference=:crust)
    @test isapprox(rb_crust_bulk, rb_crust_mantle + rb_crust_core; atol=1e-12)

    # Crust reference conservation across venting
    rb_crust_rock0 = Erebus.compute_redox_budget(c_rock0; reference=:crust)
    rb_crust_rock1 = Erebus.compute_redox_budget(c_rock1; reference=:crust)
    rb_crust_vent = Erebus.compute_redox_budget(c_vent; reference=:crust)
    @test isapprox(rb_crust_rock0, rb_crust_rock1 + rb_crust_vent; atol=1e-12)

    # 6. Struct domain contracts
    @test_throws DomainError Erebus.RedoxComponents(n_Fe0=-1.0)
    @test_throws DomainError Erebus.RedoxComponents(n_H2=NaN)
    @test_throws ArgumentError Erebus.compute_redox_budget(c_fe0; reference=:INVALID_REF)
    @test_throws DomainError Erebus.compute_specific_redox_budget(
        c_fe0, -10.0; reference=:mantle
    )
    @test_throws DomainError Erebus.compute_specific_redox_budget(
        c_fe0, 0.0; reference=:mantle
    )
end

@testset "Marker Redox Components and Local Delta IW Dynamics" begin
    # 1. marker_redox_components adapter construction
    c_m = Erebus.marker_redox_components(1.5, 3.0, 0.2; n_H2O=10.0, n_CO2=2.0)
    @test isapprox(c_m.n_Fe0, 1.5; atol=1e-12)
    @test isapprox(c_m.n_Fe2, 3.0; atol=1e-12)
    @test isapprox(c_m.n_Fe3, 0.2; atol=1e-12)
    @test isapprox(c_m.n_H2O, 10.0; atol=1e-12)
    @test isapprox(c_m.n_CO2, 2.0; atol=1e-12)
    @test isapprox(c_m.n_H2, 0.0; atol=1e-12)

    # Negative inputs are clamped non-negative
    c_neg = Erebus.marker_redox_components(-1.0, 2.0, 0.0)
    @test isapprox(c_neg.n_Fe0, 0.0; atol=1e-12)

    # 2. local_delta_iw: Metal-saturated regime
    T = 1500.0
    P = 1.0e8
    c_metal = Erebus.marker_redox_components(1.0, 0.2, 0.0)
    diw_metal = Erebus.local_delta_iw(c_metal, T, P)
    @test diw_metal < 0.0
    # Thermodynamic target: x_FeO = 0.2 / 1.2 = 1/6; deltaIW = 2 * log10(1/6) ≈ -1.5563
    @test isapprox(diw_metal, -1.5563; atol=1e-3)

    # Upper clamp limit respected even if negative
    diw_clamped = Erebus.local_delta_iw(c_metal, T, P; deltaIW_max=-2.0)
    @test isapprox(diw_clamped, -2.0; atol=1e-6)

    # 3. local_delta_iw: Silicate ferric/ferrous buffer regime (metal-absent)
    # Reference neutral state with default initial_x_ferric=0.05 (Fe3+/Fe2+ = 0.05/0.95)
    c_neutral = Erebus.marker_redox_components(0.0, 0.95, 0.05)
    diw_neutral = Erebus.local_delta_iw(c_neutral, T, P)
    @test isapprox(diw_neutral, 0.0; atol=1e-6)

    # Custom initial_x_ferric reference ratio consistency
    diw_custom_ref = Erebus.local_delta_iw(
        Erebus.marker_redox_components(0.0, 0.90, 0.10), T, P; initial_x_ferric=0.10
    )
    @test isapprox(diw_custom_ref, 0.0; atol=1e-6)

    # Oxidized silicate: Fe3+ / Fe2+ = 0.20 / 0.80 = 0.25 (should buffer > 0)
    c_ox = Erebus.marker_redox_components(0.0, 0.80, 0.20)
    diw_ox = Erebus.local_delta_iw(c_ox, T, P)
    @test diw_ox > diw_neutral
    @test isapprox(diw_ox, 2.7068; atol=1e-3)

    # Deeply reduced silicate without Fe3+: clamps to deltaIW_min
    c_red = Erebus.marker_redox_components(0.0, 1.0, 0.0)
    diw_red = Erebus.local_delta_iw(c_red, T, P; deltaIW_min=-5.0)
    @test isapprox(diw_red, -5.0; atol=1e-12)

    # Completely oxidized silicate without Fe2+: clamps to deltaIW_max
    c_all_fe3 = Erebus.marker_redox_components(0.0, 0.0, 1.0)
    diw_all_fe3 = Erebus.local_delta_iw(c_all_fe3, T, P; deltaIW_max=5.0)
    @test isapprox(diw_all_fe3, 5.0; atol=1e-12)

    # 4. Continuity across metallic iron exhaustion (n_Fe0 -> 0)
    # Ensure smooth transition without discontinuous jumps
    fe0_sweep = [1.0e-2, 1.0e-3, 5.0e-4, 2.0e-4, 1.0e-4, 5.0e-5, 1.0e-5, 0.0]
    diw_prev = Erebus.local_delta_iw(
        Erebus.marker_redox_components(fe0_sweep[1], 0.80, 0.20), T, P
    )
    for fe0 in fe0_sweep[2:end]
        c_step = Erebus.marker_redox_components(fe0, 0.80, 0.20)
        diw_curr = Erebus.local_delta_iw(c_step, T, P)
        @test diw_curr >= diw_prev - 1e-6 # Monotonic oxidation as reducing metal depletes
        @test abs(diw_curr - diw_prev) < 2.0 # Smooth, continuous transition
        diw_prev = diw_curr
    end

    # 5. Dynamic update_marker_redox! response to core segregation and serpentinization
    cfg_rdx = RedoxConfig(;
        active=true, segregation_redox=true, serpentinization_redox=true
    )
    props = Erebus.setup_marker_redox_properties(2, cfg_rdx; initial_xfe_bulk=[0.3, 0.3])
    @test props.deltaIW_m[1] <= 0.0
    @test props.deltaIW_m[2] <= 0.0

    # Core segregation drains metal to 0 on marker 1, marker 2 stays metal-bearing
    tkm_test = [1200.0, 1200.0]
    pfm_test = [1.0e7, 1.0e7]
    Erebus.update_marker_redox!(
        props, tkm_test, pfm_test, cfg_rdx; Xfem=[0.0, 0.3], XWsolidm=[0.0, 0.0]
    )
    @test isapprox(props.nFe0_m[1], 0.0; atol=1e-12)
    @test props.nFe0_m[2] > 0.0
    @test props.deltaIW_m[1] > props.deltaIW_m[2]

    # Serpentinization on marker 1 oxidizes Fe2+ to Fe3+, raising deltaIW further
    diw_before_serp = props.deltaIW_m[1]
    Erebus.update_marker_redox!(
        props, tkm_test, pfm_test, cfg_rdx; Xfem=[0.0, 0.3], XWsolidm=[0.8, 0.0]
    )
    @test props.deltaIW_m[1] > diw_before_serp
    @test props.nFe3_m[1] > props.nFe3_m[2]

    # 6. local_delta_iw: Gas buffer fallback (Fe-free, H2O + H2)
    c_gas = Erebus.marker_redox_components(0.0, 0.0, 0.0; n_H2=1.0, n_H2O=10.0)
    diw_gas = Erebus.local_delta_iw(c_gas, T, P)
    @test isfinite(diw_gas)

    # 7. Error contracts
    @test_throws DomainError Erebus.local_delta_iw(c_neutral, -100.0, P)
    @test_throws DomainError Erebus.local_delta_iw(c_neutral, T, -1.0)
    @test_throws ArgumentError Erebus.local_delta_iw(
        c_neutral, T, P; deltaIW_min=2.0, deltaIW_max=-2.0
    )

    # 6. Complete End-to-End Differentiation Electron Conservation Test
    # Stage 0: Primordial bulk assemblage (Fe0 metal, silicates, water, organics)
    c_initial = Erebus.RedoxComponents(;
        n_Fe0=50.0,
        n_Fe2=100.0,
        n_Fe3=5.0,
        n_H2=0.0,
        n_H2O=20.0,
        n_C_graphite=10.0,
        n_CO=0.0,
        n_CO2=2.0,
        n_CH4=0.0,
    )
    rb_initial = Erebus.compute_redox_budget(c_initial; reference=:mantle)

    # Stage 1: Hydrothermal serpentinization consumes 10 mol H2O -> oxidizes Fe2+ to Fe3+ + produces H2
    c_serp = Erebus.serpentinize_redox_budget(c_initial, 10.0)
    rb_serp = Erebus.compute_redox_budget(c_serp; reference=:mantle)
    @test isapprox(rb_serp, rb_initial; atol=1e-12)

    # Stage 2: Core segregation (80% metal extracts to core)
    c_mantle, c_core = Erebus.segregate_core_redox_budget(c_serp, 0.80)
    rb_core = Erebus.compute_redox_budget(c_core; reference=:mantle)
    rb_mantle = Erebus.compute_redox_budget(c_mantle; reference=:mantle)
    @test isapprox(rb_mantle + rb_core, rb_initial; atol=1e-12)

    # Stage 3: Gas venting (degas all produced H2 and CO2)
    c_vent = Erebus.RedoxComponents(; n_H2=c_mantle.n_H2, n_CO2=c_mantle.n_CO2)
    c_crust_residue = Erebus.vent_gas_redox_budget(c_mantle, c_vent)
    rb_vent = Erebus.compute_redox_budget(c_vent; reference=:mantle)
    rb_crust_residue = Erebus.compute_redox_budget(c_crust_residue; reference=:mantle)

    # Total conservation check: Residue + Core + Vented Gas == Initial
    @test isapprox(rb_crust_residue + rb_core + rb_vent, rb_initial; atol=1e-12)

    # 7. setup_marker_redox_properties initialization
    cfg_inactive = Erebus.RedoxConfig(; active=false)
    p_inact = Erebus.setup_marker_redox_properties(10, cfg_inactive)
    @test p_inact.nFe0_m === nothing
    @test p_inact.deltaIW_m === nothing
    @test p_inact.nC_graphite_m === nothing
    @test p_inact.nCO_m === nothing
    @test p_inact.nCO2_m === nothing
    @test p_inact.nCH4_m === nothing

    cfg_act = Erebus.RedoxConfig(; active=true, initial_x_ferric=0.08)
    xfe_init = fill(0.20, 10)
    p_act = Erebus.setup_marker_redox_properties(10, cfg_act; initial_xfe_bulk=xfe_init)
    @test length(p_act.nFe0_m) == 10
    @test length(p_act.nFe2_m) == 10
    @test length(p_act.nFe3_m) == 10
    @test length(p_act.deltaIW_m) == 10
    @test length(p_act.nC_graphite_m) == 10
    @test length(p_act.nCO_m) == 10
    @test all(p_act.nFe0_m .> 0.0)
    @test all(p_act.nFe2_m .> 0.0)
    @test all(p_act.nFe3_m .> 0.0)
    @test all(p_act.deltaIW_m .< 0.0) # metal present -> deltaIW negative
    @test all(isapprox.(p_act.nC_graphite_m, 0.0; atol=1e-12))
    @test all(isapprox.(p_act.nCO_m, 0.0; atol=1e-12))

    # Domain error contract on negative marknum
    @test_throws DomainError Erebus.setup_marker_redox_properties(-5, cfg_act)
end

@testset "Organic Carbon Pyrolysis Electron Conservation (Evans 2012)" begin
    # Initial rock reservoir containing metal, iron oxides, and water
    c_init = Erebus.RedoxComponents(;
        n_Fe0=1.0,
        n_Fe2=5.0,
        n_Fe3=0.8,
        n_H2=0.0,
        n_H2O=2.0,
        n_C_graphite=0.0,
        n_CO=0.0,
        n_CO2=0.0,
        n_CH4=0.0,
    )
    rb_m_init = Erebus.compute_redox_budget(c_init; reference=:mantle)

    # 1. Pure graphite residue formation (nu = 0 in mantle reference)
    c_gr = Erebus.pyrolyze_redox_budget(c_init, 2.0, 2.0, 0.0, 0.0, 0.0)
    @test isapprox(c_gr.n_C_graphite, 2.0; atol=1e-12)
    @test isapprox(c_gr.n_CO, 0.0; atol=1e-12)
    rb_m_gr = Erebus.compute_redox_budget(c_gr; reference=:mantle)
    @test isapprox(rb_m_gr, rb_m_init; atol=1e-12)

    # 2. Disproportionation: 2 C -> CO2 (nu = +4) + CH4 (nu = -4), net carbon e- delta = 0
    c_disp = Erebus.pyrolyze_redox_budget(c_init, 2.0, 0.0, 0.0, 1.0, 1.0)
    @test isapprox(c_disp.n_CO2, 1.0; atol=1e-12)
    @test isapprox(c_disp.n_CH4, 1.0; atol=1e-12)
    rb_m_disp = Erebus.compute_redox_budget(c_disp; reference=:mantle)
    @test isapprox(rb_m_disp, rb_m_init; atol=1e-12)

    # 3. Carbon oxidation to CO: C -> CO (+2 e-) coupled to iron reduction (Fe3+ -> Fe2+)
    c_co = Erebus.pyrolyze_redox_budget(c_init, 0.5, 0.0, 0.5, 0.0, 0.0; auto_balance=true)
    @test isapprox(c_co.n_CO, 0.5; atol=1e-12)
    # 0.5 mol CO releases 1.0 e-. Fe3+ (0.8 mol) absorbs 0.8 e-, remaining 0.2 e- reduces 0.1 mol Fe2+ -> Fe0
    @test isapprox(c_co.n_Fe3, 0.0; atol=1e-12)
    @test isapprox(c_co.n_Fe0, 1.0 + 0.1; atol=1e-12)
    rb_m_co = Erebus.compute_redox_budget(c_co; reference=:mantle)
    @test isapprox(rb_m_co, rb_m_init; atol=1e-12)

    # 4. Mixed pyrolysis products: graphite (60%), CO (20%), CO2 (10%), CH4 (10%)
    c_mix = Erebus.pyrolyze_redox_budget(c_init, 2.0, 1.2, 0.4, 0.2, 0.2; auto_balance=true)
    @test isapprox(c_mix.n_C_graphite, 1.2; atol=1e-12)
    @test isapprox(c_mix.n_CO, 0.4; atol=1e-12)
    @test isapprox(c_mix.n_CO2, 0.2; atol=1e-12)
    @test isapprox(c_mix.n_CH4, 0.2; atol=1e-12)
    rb_m_mix = Erebus.compute_redox_budget(c_mix; reference=:mantle)
    @test isapprox(rb_m_mix, rb_m_init; atol=1e-12)

    # 5. Whole-body 4-stage differentiation electron conservation
    c_s1 = Erebus.serpentinize_redox_budget(c_init, 0.5)
    c_s2 = Erebus.pyrolyze_redox_budget(c_s1, 1.0, 0.6, 0.3, 0.1, 0.0; auto_balance=true)
    c_mantle, c_core = Erebus.segregate_core_redox_budget(c_s2, 0.90)
    c_gas_vent = Erebus.RedoxComponents(;
        n_H2=c_mantle.n_H2, n_CO=c_mantle.n_CO, n_CO2=c_mantle.n_CO2
    )
    c_final_rock = Erebus.vent_gas_redox_budget(c_mantle, c_gas_vent)

    rb_rock = Erebus.compute_redox_budget(c_final_rock; reference=:mantle)
    rb_core = Erebus.compute_redox_budget(c_core; reference=:mantle)
    rb_gas = Erebus.compute_redox_budget(c_gas_vent; reference=:mantle)
    @test isapprox(rb_rock + rb_core + rb_gas, rb_m_init; atol=1e-12)

    # 6. Gauge invariance: reference frame (:mantle vs :crust) yields identical reaction stoichiometry
    c_crust = Erebus.pyrolyze_redox_budget(
        c_init, 0.5, 0.0, 0.5, 0.0, 0.0; auto_balance=true, reference=:crust
    )
    @test isapprox(c_crust.n_Fe3, c_co.n_Fe3; atol=1e-12)
    @test isapprox(c_crust.n_Fe0, c_co.n_Fe0; atol=1e-12)

    # 7. Error contracts
    @test_throws DomainError Erebus.pyrolyze_redox_budget(c_init, -1.0, 0.0, 0.0, 0.0, 0.0)
    @test_throws DomainError Erebus.pyrolyze_redox_budget(c_init, 1.0, 0.5, 0.0, 0.0, 0.0)
    @test_throws DomainError Erebus.pyrolyze_redox_budget(
        Erebus.RedoxComponents(), 1.0, 0.0, 1.0, 0.0, 0.0; auto_balance=true
    )
    # Insufficient rock reductants during carbon reduction (CH4 production without rock reductants)
    c_all_fe3 = Erebus.RedoxComponents(; n_Fe3=1.0)
    @test_throws DomainError Erebus.pyrolyze_redox_budget(
        c_all_fe3, 1.0, 0.0, 0.0, 0.0, 1.0; auto_balance=true
    )
    @test_throws ArgumentError Erebus.pyrolyze_redox_budget(
        c_init, 1.0, 1.0, 0.0, 0.0, 0.0; reference=:unknown
    )
end

@testset "Graphite CCO Oxygen Fugacity Buffering and Transitions" begin
    T = 1300.0
    P_100bar = 1.0e7

    # 1. Graphite CCO buffer locking when metallic iron is depleted
    c_gr = Erebus.marker_redox_components(0.0, 1.0, 0.05; n_C_graphite=0.1)
    diw_gr = Erebus.local_delta_iw(
        c_gr, T, P_100bar; graphite_buffer_active=true, w_graphite_threshold=1.0e-6
    )

    # Theoretical CCO offset at 1300 K, 100 bar
    lfo2_cco = Erebus.log10_fo2_of_buffer(:CCO, T, P_100bar)
    lfo2_iw = Erebus.log10_fo2_of_buffer(:IW, T, P_100bar)
    diw_cco_expected = lfo2_cco - lfo2_iw
    @test isapprox(diw_gr, diw_cco_expected; atol=1e-3)
    @test diw_gr > 0.0
    @test 0.5 < diw_gr < 1.2

    # 2. Buffer bypass when graphite buffering is deactivated
    diw_no_gr = Erebus.local_delta_iw(c_gr, T, P_100bar; graphite_buffer_active=false)
    c_sil_only = Erebus.marker_redox_components(0.0, 1.0, 0.05)
    diw_sil_expected = Erebus.local_delta_iw(c_sil_only, T, P_100bar)
    @test isapprox(diw_no_gr, diw_sil_expected; atol=1e-6)

    # 3. Sensitivity to w_graphite_threshold at intermediate graphite contents
    c_int_gr = Erebus.marker_redox_components(0.0, 0.80, 0.20; n_C_graphite=2.0e-4)
    diw_high_thresh = Erebus.local_delta_iw(
        c_int_gr, T, P_100bar; graphite_buffer_active=true, w_graphite_threshold=1.0e-2
    )
    diw_low_thresh = Erebus.local_delta_iw(
        c_int_gr, T, P_100bar; graphite_buffer_active=true, w_graphite_threshold=1.0e-6
    )
    # Lower threshold reaches full CCO buffer while high threshold stays closer to oxidized silicate
    @test diw_low_thresh < diw_high_thresh
    @test isapprox(diw_low_thresh, diw_cco_expected; atol=1e-2)

    # 4. Smooth continuous transition across metal exhaustion: IW -> CCO
    fe0_sweep = [1.0e-2, 1.0e-3, 5.0e-4, 2.0e-4, 1.0e-4, 5.0e-5, 1.0e-5, 0.0]
    diw_prev = Erebus.local_delta_iw(
        Erebus.marker_redox_components(fe0_sweep[1], 1.0, 0.05; n_C_graphite=0.1),
        T,
        P_100bar;
        graphite_buffer_active=true,
    )
    for fe0 in fe0_sweep[2:end]
        c_step = Erebus.marker_redox_components(fe0, 1.0, 0.05; n_C_graphite=0.1)
        diw_curr = Erebus.local_delta_iw(c_step, T, P_100bar; graphite_buffer_active=true)
        @test diw_curr >= diw_prev - 1e-6
        @test abs(diw_curr - diw_prev) < 2.0
        diw_prev = diw_curr
    end
    @test isapprox(diw_prev, diw_cco_expected; atol=1e-3)

    # 5. Smooth continuous transition across graphite exhaustion: CCO -> QFM
    c_ox_sil = Erebus.marker_redox_components(0.0, 0.80, 0.20)
    diw_ox_sil = Erebus.local_delta_iw(c_ox_sil, T, P_100bar)
    gr_sweep = [1.0e-1, 1.0e-2, 1.0e-3, 5.0e-4, 1.0e-4, 5.0e-5, 1.0e-5, 0.0]
    diw_prev_gr = Erebus.local_delta_iw(
        Erebus.marker_redox_components(0.0, 0.80, 0.20; n_C_graphite=gr_sweep[1]),
        T,
        P_100bar;
        graphite_buffer_active=true,
    )
    for gr in gr_sweep[2:end]
        c_step = Erebus.marker_redox_components(0.0, 0.80, 0.20; n_C_graphite=gr)
        diw_curr = Erebus.local_delta_iw(c_step, T, P_100bar; graphite_buffer_active=true)
        @test diw_curr >= diw_prev_gr - 1e-6
        @test abs(diw_curr - diw_prev_gr) < 2.0
        diw_prev_gr = diw_curr
    end
    @test isapprox(diw_prev_gr, diw_ox_sil; atol=1e-3)
end

@testset "Multi-Timestep Pyrolysis and Redox Coupling Invariants" begin
    tkm = [1300.0]
    dt = 100.0
    phim = [0.1]
    X_refr_C_m = [0.5]
    X_refr_N_m = [0.001]
    X_refr_H_m = [0.002]

    cfg_ref = Erebus.RefractoryConfig(; active=true, kinetics_active=true, f_refr_C=0.6)
    cfg_rdx = Erebus.RedoxConfig(; active=true, pyrolysis_redox=true)

    props = Erebus.setup_marker_redox_properties(1, cfg_rdx; tkm=tkm, pfm=[1.0e7])
    fe3_init = props.nFe3_m[1]
    @test fe3_init > 0.0

    # Execute 5 consecutive timesteps
    for step in 1:5
        Erebus.update_marker_pyrolysis!(
            tkm,
            dt,
            phim,
            X_refr_C_m,
            X_refr_N_m,
            X_refr_H_m,
            cfg_ref;
            redox_props=props,
            redox_cfg=cfg_rdx,
        )
        Erebus.update_marker_redox!(
            props, tkm, [1.0e7], cfg_rdx; Xfem=nothing, XWsolidm=nothing
        )
    end

    # 1. Smelting reduction is preserved: Fe3+ stays depleted and does not revert
    @test isapprox(props.nFe3_m[1], 0.0; atol=1e-12)
    @test props.nFe0_m[1] > 0.0

    # 2. Multi-species carbon speciation: both CO and CO2 are formed
    @test props.nCO_m[1] > 0.0
    @test props.nCO2_m[1] > 0.0
    @test props.nC_graphite_m[1] > 0.0

    # 3. Oxygen starvation handling: total carbon is bounded and non-negative
    @test isfinite(props.nC_graphite_m[1])
    @test isfinite(props.nCO_m[1])
    @test isfinite(props.nCO2_m[1])
end

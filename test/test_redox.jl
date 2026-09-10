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
    # Test CCO buffer across pressures at T = 1300 K
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

    # Scale guard: at 1300 K and 10 bar, CCO sits near IW to IW+2
    lfo2_iw = Erebus.log10_fo2_of_buffer(:IW, T, P_10bar)
    delta_cco_iw = lfo2_cco_10 - lfo2_iw
    @test -2.0 < delta_cco_iw < 4.0
    @test lfo2_cco_10 < 0.0
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
end

@testset "Evans 2012 Redox Budget Electron Accounting" begin
    # Reference: Evans (2012) Earth-Sci. Rev. 113, 11-32, DOI: 10.1016/j.earscirev.2012.03.003
    # Mantle reference state (M): Fe2+, C0, S2-, H+, O2-, P5+
    # Crust reference state (C): Fe3+, C4+, S6+, H+, O2-, P5+

    # 1. Pure component electron counts relative to Mantle reference state
    # 1 mol Fe0 -> nu = -2 mol e-
    c_fe0 = Erebus.RedoxComponents(n_Fe0=1.0)
    rb_fe0 = Erebus.compute_redox_budget(c_fe0; reference=:mantle)
    @test isapprox(rb_fe0, -2.0; atol=1e-12)

    # 1 mol Fe3+ -> nu = +1 mol e-
    c_fe3 = Erebus.RedoxComponents(n_Fe3=1.0)
    rb_fe3 = Erebus.compute_redox_budget(c_fe3; reference=:mantle)
    @test isapprox(rb_fe3, 1.0; atol=1e-12)

    # 1 mol H2 -> nu = -2 mol e-
    c_h2 = Erebus.RedoxComponents(n_H2=1.0)
    rb_h2 = Erebus.compute_redox_budget(c_h2; reference=:mantle)
    @test isapprox(rb_h2, -2.0; atol=1e-12)

    # 1 mol H2O -> nu = 0.0
    c_h2o = Erebus.RedoxComponents(n_H2O=1.0)
    rb_h2o = Erebus.compute_redox_budget(c_h2o; reference=:mantle)
    @test isapprox(rb_h2o, 0.0; atol=1e-12)

    # 1 mol CO2 -> nu = +4 mol e-
    c_co2 = Erebus.RedoxComponents(n_CO2=1.0)
    rb_co2 = Erebus.compute_redox_budget(c_co2; reference=:mantle)
    @test isapprox(rb_co2, 4.0; atol=1e-12)

    # 1 mol CH4 -> nu = -4 mol e-
    c_ch4 = Erebus.RedoxComponents(n_CH4=1.0)
    rb_ch4 = Erebus.compute_redox_budget(c_ch4; reference=:mantle)
    @test isapprox(rb_ch4, -4.0; atol=1e-12)

    # Specific redox budget: RB / mass_kg
    mass = 100.0 # kg
    rb_spec = Erebus.compute_specific_redox_budget(c_fe0, mass; reference=:mantle)
    @test isapprox(rb_spec, -0.02; atol=1e-12)

    # 2. Conservation Invariant: Serpentinization Reaction
    # 3 FeO (Fe2+) + H2O -> Fe3O4 (1 Fe2+ + 2 Fe3+) + H2
    # Before: 3 mol Fe2+ (0) + 1 mol H2O (0) => RB = 0.0
    # After: 2 mol Fe3+ (+2) + 1 mol H2 (-2) => RB = 0.0
    c_before = Erebus.RedoxComponents(n_Fe2=3.0, n_H2O=1.0)
    c_after = Erebus.serpentinize_redox_budget(c_before, 1.0) # react 1 mol H2O
    rb_init = Erebus.compute_redox_budget(c_before; reference=:mantle)
    rb_final = Erebus.compute_redox_budget(c_after; reference=:mantle)
    @test isapprox(rb_init, rb_final; atol=1e-12)
    @test isapprox(c_after.n_Fe3, 2.0; atol=1e-12)
    @test isapprox(c_after.n_H2, 1.0; atol=1e-12)

    # 3. Conservation Invariant: Core Segregation
    # Metal segregation moves Fe0 to core, conserving total whole-body electrons
    c_bulk = Erebus.RedoxComponents(n_Fe0=50.0, n_Fe2=100.0, n_Fe3=5.0)
    c_mantle, c_core = Erebus.segregate_core_redox_budget(c_bulk, 0.8) # 80% metal to core
    rb_bulk = Erebus.compute_redox_budget(c_bulk; reference=:mantle)
    rb_mantle = Erebus.compute_redox_budget(c_mantle; reference=:mantle)
    rb_core = Erebus.compute_redox_budget(c_core; reference=:mantle)
    @test isapprox(rb_bulk, rb_mantle + rb_core; atol=1e-12)
    # Core is highly reduced (negative RB_M)
    @test rb_core < 0.0
    # Mantle is more oxidized after metal segregation
    @test rb_mantle > rb_bulk

    # 4. Conservation Invariant: Degassing and Gas Venting
    # Loss of reduced gas (H2, CO) leaves an oxidized residual rock
    c_rock0 = Erebus.RedoxComponents(n_Fe2=100.0, n_H2=10.0, n_CO=5.0, n_CO2=5.0)
    c_vent = Erebus.RedoxComponents(n_H2=8.0, n_CO=3.0)
    c_rock1 = Erebus.vent_gas_redox_budget(c_rock0, c_vent)
    rb_rock0 = Erebus.compute_redox_budget(c_rock0; reference=:mantle)
    rb_rock1 = Erebus.compute_redox_budget(c_rock1; reference=:mantle)
    rb_vent = Erebus.compute_redox_budget(c_vent; reference=:mantle)
    @test isapprox(rb_rock0, rb_rock1 + rb_vent; atol=1e-12)
    # Vented gas is reducing (rb_vent < 0), so rock becomes more oxidized (rb_rock1 > rb_rock0)
    @test rb_vent < 0.0
    @test rb_rock1 > rb_rock0

    # 5. Crust reference state conversion invariance
    # In crust reference state (Fe3+, C4+, S6+):
    # Fe0 has nu = -3, Fe2+ has nu = -1, Fe3+ has nu = 0
    rb_fe0_crust = Erebus.compute_redox_budget(c_fe0; reference=:crust)
    @test isapprox(rb_fe0_crust, -3.0; atol=1e-12)

    # 6. Error contracts
    @test_throws ArgumentError Erebus.compute_redox_budget(c_fe0; reference=:INVALID_REF)
    @test_throws DomainError Erebus.compute_specific_redox_budget(
        c_fe0, -10.0; reference=:mantle
    )
    @test_throws DomainError Erebus.compute_specific_redox_budget(
        c_fe0, 0.0; reference=:mantle
    )
end

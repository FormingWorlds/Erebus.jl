using Test
using Erebus

@testset "Redox and Metal Segregation Tests" begin
    @testset "Solid Metal Fe(0) Invariant During Thermal Cycling" begin
        # Setup markers with solid metallic iron at sub-eutectic temperature
        tkm = [800.0]
        pfm = [1.0e7]
        Xfe_bulk = [0.20]
        Xfem = [0.0] # Sub-eutectic: molten metal volume fraction is zero
        rho_s = 3000.0
        rho_metal = 7000.0
        rho_ratio = rho_metal / rho_s
        redox_cfg = RedoxConfig(; active=true, segregation_redox=true)
        redox_props = setup_marker_redox_properties(
            1, redox_cfg; initial_xfe_bulk=Xfe_bulk, rhosolid=rho_s, rho_metal=rho_metal
        )

        # Baseline: initialize marker redox state at 800 K
        update_marker_redox!(
            redox_props,
            tkm,
            pfm,
            redox_cfg;
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            rhosolid=rho_s,
            rho_metal=rho_metal,
        )

        n_fe0_expected = (0.20 * rho_ratio) / Erebus.M_Fe
        @test isapprox(redox_props.nFe0_m[1], n_fe0_expected; atol=1e-12)

        # Thermal cycling: heat above Fe-FeS eutectic (1300 K)
        tkm[1] = 1300.0
        Xfem[1] = 0.20 # Fully molten metal
        update_marker_redox!(
            redox_props,
            tkm,
            pfm,
            redox_cfg;
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            rhosolid=rho_s,
            rho_metal=rho_metal,
        )
        @test isapprox(redox_props.nFe0_m[1], n_fe0_expected; atol=1e-12)

        # Cool back below eutectic to 800 K: solid metal inventory must be preserved
        tkm[1] = 800.0
        Xfem[1] = 0.0 # Re-solidified metal
        update_marker_redox!(
            redox_props,
            tkm,
            pfm,
            redox_cfg;
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            rhosolid=rho_s,
            rho_metal=rho_metal,
        )
        @test isapprox(redox_props.nFe0_m[1], n_fe0_expected; atol=1e-12)

        # Metal segregation: bulk metal fraction drains from 0.20 to 0.05
        Xfe_bulk[1] = 0.05
        update_marker_redox!(
            redox_props,
            tkm,
            pfm,
            redox_cfg;
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            rhosolid=rho_s,
            rho_metal=rho_metal,
        )
        n_fe0_drained = (0.05 * rho_ratio) / Erebus.M_Fe
        @test isapprox(redox_props.nFe0_m[1], n_fe0_drained; atol=1e-12)
        @test redox_props.nFe0_m[1] < n_fe0_expected
    end

    @testset "IOM Smelted Iron Preservation and Density Scaling" begin
        # Setup marker at 1000 K where organic pyrolysis actively smelts iron oxide
        # 1000 K is above the 900 K methanation window, so CO/CO2 oxidize carbon
        # and reduce FeO to solid Fe(0) below eutectic (1213 K)
        tkm = [1000.0]
        dt = 1.0e6 * 365.25 * 86400.0
        phim = [0.10]
        X_refr_C = [0.03]
        X_refr_N = [0.001]
        X_refr_H = [0.002]
        refractory_cfg = RefractoryConfig(; active=true, kinetics_active=true)

        redox_cfg = RedoxConfig(; active=true, pyrolysis_redox=true, segregation_redox=true)
        Xfe_bulk = [0.05]
        Xfem = [0.0]
        rho_s = 3000.0
        rho_metal = 7000.0
        rho_ratio = rho_metal / rho_s
        redox_props = setup_marker_redox_properties(
            1, redox_cfg; initial_xfe_bulk=Xfe_bulk, rhosolid=rho_s, rho_metal=rho_metal
        )

        # Initialise iron inventories before smelting (Fe3 = 0 so electrons reduce Fe2 to Fe0)
        redox_props.nFe3_m[1] = 0.0
        redox_props.nFe2_m[1] = 2.0
        redox_props.nFe0_m[1] = (0.05 * rho_ratio) / Erebus.M_Fe

        Xfe_bulk_initial = Xfe_bulk[1]
        update_marker_pyrolysis!(
            tkm,
            dt,
            phim,
            X_refr_C,
            X_refr_N,
            X_refr_H,
            refractory_cfg;
            redox_props=redox_props,
            redox_cfg=redox_cfg,
            Xfem=Xfem,
            Xfe_bulk=Xfe_bulk,
            rhosolid=rho_s,
            rho_metal=rho_metal,
        )

        # Pyrolysis must produce smelted metallic iron and increase Xfe_bulk
        @test Xfe_bulk[1] > Xfe_bulk_initial
        @test isapprox(Xfem[1], 0.0; atol=1e-15)

        # Verify density scaling: dphi_fe0 = dw_fe0 * (rho_s / rho_metal)
        dphi_fe = Xfe_bulk[1] - Xfe_bulk_initial
        smelted_fe0 = redox_props.nFe0_m[1]
        expected_dphi_fe =
            (smelted_fe0 - (0.05 * rho_ratio) / Erebus.M_Fe) * Erebus.M_Fe * (rho_s / rho_metal)
        @test isapprox(dphi_fe, expected_dphi_fe; rtol=1e-12)
        @test isapprox(Xfe_bulk[1], Xfe_bulk_initial + expected_dphi_fe; rtol=1e-12)

        # Verify preservation across follow-up update_marker_redox! call
        update_marker_redox!(
            redox_props,
            tkm,
            [1.0e7],
            redox_cfg;
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            rhosolid=rho_s,
            rho_metal=rho_metal,
        )
        @test isapprox(redox_props.nFe0_m[1], smelted_fe0; atol=1e-12)
    end

    @testset "Metal-Silicate Volatile Partitioning Mass Conservation" begin
        # Setup marker with coexisting metal and silicate melt
        m = 1
        T_val = 1600.0 # Above silicate solidus and metal liquidus
        P_val = 1.0e8  # 100 MPa
        ΔIW = -2.0
        F_fe = 0.50
        F_melt = 0.20
        rho_sil = 3000.0
        rho_met = 7000.0

        Xfe_bulk = [0.15]
        Xfem = [0.15 * F_fe]
        XH2Om = [1.5] # wt%
        XC = [500.0] # ppmw
        XN = [50.0]  # ppmw
        XS = [1000.0] # ppmw

        init_fe_C = 50.0
        init_fe_S = 200.0
        Xfe_H = [10.0]
        Xfe_C = [init_fe_C]
        Xfe_N = [5.0]
        Xfe_S = [init_fe_S]

        cfg = MetalPartitionConfig(; active=true, equilibration_rate=0.80)

        # Total initial volatile masses per unit volume
        phi_fe = Xfe_bulk[1]
        phi_sil = 1.0 - phi_fe
        m_sil = phi_sil * rho_sil
        m_met = phi_fe * rho_met

        M_tot_C_init = m_sil * XC[1] + m_met * Xfe_C[1]
        M_tot_N_init = m_sil * XN[1] + m_met * Xfe_N[1]
        M_tot_S_init = m_sil * XS[1] + m_met * Xfe_S[1]

        equilibrate_metal_silicate_volatiles!(
            m,
            F_fe,
            F_melt,
            T_val,
            P_val,
            ΔIW,
            Xfe_bulk,
            Xfem,
            XH2Om,
            XC,
            XN,
            XS,
            Xfe_H,
            Xfe_C,
            Xfe_N,
            Xfe_S,
            cfg;
            rho_silicate=rho_sil,
            rho_metal=rho_met,
        )

        # Strict mass conservation checks across phase exchange
        M_tot_C_final = m_sil * XC[1] + m_met * Xfe_C[1]
        M_tot_N_final = m_sil * XN[1] + m_met * Xfe_N[1]
        M_tot_S_final = m_sil * XS[1] + m_met * Xfe_S[1]

        @test isapprox(M_tot_C_final, M_tot_C_init; rtol=1e-12)
        @test isapprox(M_tot_N_final, M_tot_N_init; rtol=1e-12)
        @test isapprox(M_tot_S_final, M_tot_S_init; rtol=1e-12)

        # Siderophile volatile uptake in metal phase
        @test Xfe_C[1] > init_fe_C
        @test Xfe_S[1] > init_fe_S
    end

    @testset "Partial Melt Denominator Stability at Incipient Melting" begin
        # Incipient melting: F_melt = 1.0e-5 (near zero)
        m = 1
        T_val = 1450.0
        P_val = 2.0e8
        ΔIW = -1.5
        F_fe = 0.80
        F_melt = 1.0e-5
        rho_sil = 3000.0
        rho_met = 7000.0

        Xfe_bulk = [0.10]
        Xfem = [0.08]
        XH2O = [0.5]
        XC = [200.0]
        XN = [30.0]
        XS = [500.0]

        Xfe_H = [5.0]
        Xfe_C = [20.0]
        Xfe_N = [2.0]
        Xfe_S = [100.0]

        cfg = MetalPartitionConfig(; active=true, equilibration_rate=1.0)

        phi_fe = Xfe_bulk[1]
        phi_sil = 1.0 - phi_fe
        m_sil = phi_sil * rho_sil
        m_met = phi_fe * rho_met
        M_tot_C_init = m_sil * XC[1] + m_met * Xfe_C[1]

        # Verify no singularity or NaN occurs at incipient melt
        equilibrate_metal_silicate_volatiles!(
            m,
            F_fe,
            F_melt,
            T_val,
            P_val,
            ΔIW,
            Xfe_bulk,
            Xfem,
            XH2O,
            XC,
            XN,
            XS,
            Xfe_H,
            Xfe_C,
            Xfe_N,
            Xfe_S,
            cfg;
            rho_silicate=rho_sil,
            rho_metal=rho_met,
        )

        @test isfinite(XC[1])
        @test isfinite(Xfe_C[1])
        @test XC[1] <= M_tot_C_init / m_sil
        @test Xfe_C[1] <= M_tot_C_init / m_met
        M_tot_C_final = m_sil * XC[1] + m_met * Xfe_C[1]
        @test isapprox(M_tot_C_final, M_tot_C_init; rtol=1e-12)
    end

    @testset "Core Membership Classification on Undifferentiated vs Differentiated Bodies" begin
        # 1. Undifferentiated cold body with uniform chondritic metal (Xfe_bulk = 0.05)
        # All markers are within rplanet, but none exceed phi_core_threshold (0.40)
        marknum = 100
        R_planet = 50000.0
        xc = 70000.0
        yc = 70000.0
        rho_metal = 7000.0

        xm = fill(xc, marknum)
        ym = range(yc - 0.8 * R_planet, yc + 0.8 * R_planet; length=marknum)
        tm = fill(1, marknum)
        Xfe_bulk_undiff = fill(0.05, marknum)
        Xfe_H_m = fill(10.0, marknum)
        Xfe_C_m = fill(100.0, marknum)
        Xfe_N_m = fill(20.0, marknum)
        Xfe_S_m = fill(5000.0, marknum)
        w3d_m = fill(2.0 * R_planet, marknum)

        budgets_undiff = compute_core_volatile_budgets(
            xm,
            collect(ym),
            tm,
            Xfe_bulk_undiff,
            Xfe_H_m,
            Xfe_C_m,
            Xfe_N_m,
            Xfe_S_m,
            marknum;
            xcenter=xc,
            ycenter=yc,
            rplanet=R_planet,
            rho_metal=rho_metal,
            core_radius_fraction=0.50,
            phi_core_threshold=0.40,
            w3d_m=w3d_m,
            V_marker=1.0e6,
        )

        # Core metal mass must be identically 0.0 on undifferentiated body
        @test isapprox(budgets_undiff.M_core_metal, 0.0; atol=1e-12)
        @test isapprox(budgets_undiff.M_core_C, 0.0; atol=1e-12)
        expected_undiff_metal = marknum * 0.05 * rho_metal * (1.0e6 * 2.0 * R_planet)
        @test isapprox(budgets_undiff.M_total_metal, expected_undiff_metal; rtol=1e-12)

        # 2. Differentiated body with segregated core in central 20 markers
        Xfe_bulk_diff = copy(Xfe_bulk_undiff)
        core_indices = 41:60
        Xfe_bulk_diff[core_indices] .= 0.80

        budgets_diff = compute_core_volatile_budgets(
            xm,
            collect(ym),
            tm,
            Xfe_bulk_diff,
            Xfe_H_m,
            Xfe_C_m,
            Xfe_N_m,
            Xfe_S_m,
            marknum;
            xcenter=xc,
            ycenter=yc,
            rplanet=R_planet,
            rho_metal=rho_metal,
            core_radius_fraction=0.50,
            phi_core_threshold=0.40,
            w3d_m=w3d_m,
            V_marker=1.0e6,
        )

        # Segregated central markers must be included in core budgets
        expected_diff_core_metal =
            length(core_indices) * 0.80 * rho_metal * (1.0e6 * 2.0 * R_planet)
        @test isapprox(budgets_diff.M_core_metal, expected_diff_core_metal; rtol=1e-12)
        expected_diff_core_C = expected_diff_core_metal * (100.0 * 1.0e-6)
        @test isapprox(budgets_diff.M_core_C, expected_diff_core_C; rtol=1e-12)
        @test isapprox(budgets_diff.w_core_C_ppm, 100.0; atol=1e-10)
    end

    @testset "Buffer Oxygen Scaling with Specific Marker Masses" begin
        redox_cfg = RedoxConfig(; active=true)
        redox_props = setup_marker_redox_properties(2, redox_cfg)
        redox_props.nFe3_m[1] = 0.50 # mol/kg
        redox_props.nFe3_m[2] = 0.50 # mol/kg
        redox_props.nFe2_m[1] = 1.00 # mol/kg
        redox_props.nFe2_m[2] = 1.00 # mol/kg

        marker_indices = [1, 2]
        weights = [0.5, 0.5]
        marker_masses = [1.0e6, 2.0e6]
        dO = 100.0 # kg of oxygen to be extracted (reducing Fe3O4 to FeO)

        dn_O_total = dO / Erebus.M_O

        apply_buffer_oxygen!(
            redox_props, marker_indices, weights, dO; marker_masses=marker_masses
        )

        # Verify specific moles update matches (dn_O_total * w_k) / m_k
        dn_O_spec_1 = (dn_O_total * 0.5) / 1.0e6
        dn_O_spec_2 = (dn_O_total * 0.5) / 2.0e6

        @test isapprox(redox_props.nFe3_m[1], 0.50 - 2.0 * dn_O_spec_1; atol=1e-12)
        @test isapprox(redox_props.nFe3_m[2], 0.50 - 2.0 * dn_O_spec_2; atol=1e-12)
        @test isapprox(redox_props.nFe2_m[1], 1.00 + 2.0 * dn_O_spec_1; atol=1e-12)
        @test isapprox(redox_props.nFe2_m[2], 1.00 + 2.0 * dn_O_spec_2; atol=1e-12)

        # Default marker_masses === nothing fallback (m_mass = 1.0 kg)
        fe3_before = redox_props.nFe3_m[1]
        fe2_before = redox_props.nFe2_m[1]
        dn_O_fallback = (1.0e-4 / Erebus.M_O) * 0.5 / 1.0
        apply_buffer_oxygen!(
            redox_props, marker_indices, weights, 1.0e-4; marker_masses=nothing
        )
        @test isapprox(redox_props.nFe3_m[1], fe3_before - 2.0 * dn_O_fallback; atol=1e-12)
        @test isapprox(redox_props.nFe2_m[1], fe2_before + 2.0 * dn_O_fallback; atol=1e-12)

        # Insufficient Fe3O4 to supply oxygen demand (dO > 0)
        @test_throws DomainError apply_buffer_oxygen!(
            redox_props, marker_indices, weights, 1.0e10; marker_masses=marker_masses
        )

        # Insufficient FeO to absorb oxygen (dO < 0)
        @test_throws DomainError apply_buffer_oxygen!(
            redox_props, marker_indices, weights, -1.0e10; marker_masses=marker_masses
        )

        # Marker mass dimension validation
        @test_throws ArgumentError apply_buffer_oxygen!(
            redox_props, marker_indices, weights, dO; marker_masses=[1.0e6]
        )
        @test_throws DomainError apply_buffer_oxygen!(
            redox_props, marker_indices, weights, dO; marker_masses=[-1.0, 1.0e6]
        )
        @test_throws ArgumentError apply_buffer_oxygen!(
            redox_props, [1, 1], [0.5, 0.5], dO; marker_masses=[1.0e6, 1.0e6]
        )
    end

    @testset "High-Temperature Smelting and Molten Metal Allocation" begin
        # Smelting above Fe-FeS eutectic (1213 K) increments both Xfe_bulk and Xfem
        tkm = [1300.0]
        dt = 1.0e6 * 365.25 * 86400.0
        phim = [0.10]
        X_refr_C = [0.03]
        X_refr_N = [0.001]
        X_refr_H = [0.002]
        refractory_cfg = RefractoryConfig(; active=true, kinetics_active=true)
        redox_cfg = RedoxConfig(; active=true, pyrolysis_redox=true, segregation_redox=true)

        Xfe_bulk = [0.05]
        Xfem = [0.02]
        rho_s = 3000.0
        rho_metal = 7000.0
        rho_ratio = rho_metal / rho_s
        redox_props = setup_marker_redox_properties(
            1, redox_cfg; initial_xfe_bulk=Xfe_bulk, rhosolid=rho_s, rho_metal=rho_metal
        )
        redox_props.nFe3_m[1] = 0.0
        redox_props.nFe2_m[1] = 2.0
        redox_props.nFe0_m[1] = (0.05 * rho_ratio) / Erebus.M_Fe

        Xfe_bulk_init = Xfe_bulk[1]
        Xfem_init = Xfem[1]

        # Case 1: tm and rhosolidm with valid phase index
        tm = [1]
        rhosolidm = [3000.0, 3200.0]
        update_marker_pyrolysis!(
            tkm,
            dt,
            phim,
            X_refr_C,
            X_refr_N,
            X_refr_H,
            refractory_cfg;
            redox_props=redox_props,
            redox_cfg=redox_cfg,
            Xfem=Xfem,
            Xfe_bulk=Xfe_bulk,
            rhosolid=rho_s,
            rho_metal=rho_metal,
            tm=tm,
            rhosolidm=rhosolidm,
        )

        F_fe = Erebus.compute_metal_melt_fraction(1300.0)
        dphi_fe = Xfe_bulk[1] - Xfe_bulk_init
        @test isapprox(dphi_fe, (redox_props.nFe0_m[1] - (0.05 * rho_ratio) / Erebus.M_Fe) * Erebus.M_Fe * (rho_s / rho_metal); rtol=1e-12)
        @test isapprox(Xfem[1], Xfem_init + dphi_fe * F_fe; rtol=1e-12)

        # Case 2: tm with out-of-bounds phase index fallback to rhosolidm[1]
        X_refr_C[1] = 0.03
        redox_props.nFe2_m[1] = 2.0
        Xfe_bulk_c1 = Xfe_bulk[1]
        Xfem_c1 = Xfem[1]
        n_fe0_c1 = redox_props.nFe0_m[1]
        tm_oob = [99]
        update_marker_pyrolysis!(
            tkm,
            dt,
            phim,
            X_refr_C,
            X_refr_N,
            X_refr_H,
            refractory_cfg;
            redox_props=redox_props,
            redox_cfg=redox_cfg,
            Xfem=Xfem,
            Xfe_bulk=Xfe_bulk,
            rhosolid=rho_s,
            rho_metal=rho_metal,
            tm=tm_oob,
            rhosolidm=rhosolidm,
        )
        dphi_fe_c2 = (redox_props.nFe0_m[1] - n_fe0_c1) * Erebus.M_Fe * (rhosolidm[1] / rho_metal)
        @test isapprox(Xfe_bulk[1], Xfe_bulk_c1 + dphi_fe_c2; rtol=1e-12)
        @test isapprox(Xfem[1], Xfem_c1 + dphi_fe_c2 * F_fe; rtol=1e-12)

        # Case 3: tm === nothing and rhosolidm !== nothing fallback
        X_refr_C[1] = 0.03
        redox_props.nFe2_m[1] = 2.0
        Xfe_bulk_c2 = Xfe_bulk[1]
        Xfem_c2 = Xfem[1]
        n_fe0_c2 = redox_props.nFe0_m[1]
        update_marker_pyrolysis!(
            tkm,
            dt,
            phim,
            X_refr_C,
            X_refr_N,
            X_refr_H,
            refractory_cfg;
            redox_props=redox_props,
            redox_cfg=redox_cfg,
            Xfem=Xfem,
            Xfe_bulk=Xfe_bulk,
            rhosolid=rho_s,
            rho_metal=rho_metal,
            tm=nothing,
            rhosolidm=rhosolidm,
        )
        dphi_fe_c3 = (redox_props.nFe0_m[1] - n_fe0_c2) * Erebus.M_Fe * (rhosolidm[1] / rho_metal)
        @test isapprox(Xfe_bulk[1], Xfe_bulk_c2 + dphi_fe_c3; rtol=1e-12)
        @test isapprox(Xfem[1], Xfem_c2 + dphi_fe_c3 * F_fe; rtol=1e-12)
    end
end

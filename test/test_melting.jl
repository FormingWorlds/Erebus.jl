using Test
using Erebus
using Erebus.Config
using Erebus.Physics
using Erebus.Particles
using StaticArrays
using TOML

@testset "Silicate Rock Melting & Magma Rheology" begin
    @testset "MeltingConfig Schema & Bounds Validation" begin
        cfg_def = MeltingConfig()
        @test cfg_def.active == false
        @test isapprox(cfg_def.T_solidus[1], 1400.0; rtol=1e-12)
        @test isapprox(cfg_def.T_liquidus[1], 1800.0; rtol=1e-12)
        @test isapprox(cfg_def.L_melt, 4.0e5; rtol=1e-12)
        @test isapprox(cfg_def.rho_melt, 2800.0; rtol=1e-12)
        @test isapprox(cfg_def.alpha_eta, 28.0; rtol=1e-12)
        @test isapprox(cfg_def.phi_crit, 0.4; rtol=1e-12)
        @test isapprox(cfg_def.eta_melt, 10.0; rtol=1e-12)
        @test iszero(cfg_def.dpdt_clapeyron)
        @test cfg_def.latent_heat_mode == :apparent_cp

        # Inclusion in SimulationConfig
        sim_cfg = default_config()
        @test sim_cfg.melting isa MeltingConfig
        @test validate_config(sim_cfg) === nothing

        # Validation bounds: invalid solidus/liquidus
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(;
                    active=true,
                    T_solidus=SVector{3,Float64}([1800.0, 1800.0, NaN]),
                    T_liquidus=SVector{3,Float64}([1400.0, 1400.0, NaN]),
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(;
                    active=true, T_solidus=SVector{3,Float64}([-100.0, 1400.0, NaN])
                ),
            ),
        )

        # Validation bounds: non-positive latent heat
        @test_throws ArgumentError validate_config(
            SimulationConfig(; melting=MeltingConfig(; active=true, L_melt=0.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; melting=MeltingConfig(; active=true, L_melt=-1.0e5))
        )

        # Validation bounds: unphysical melt density
        @test_throws ArgumentError validate_config(
            SimulationConfig(; melting=MeltingConfig(; active=true, rho_melt=0.0))
        )

        # Validation bounds: unphysical phi_crit
        @test_throws ArgumentError validate_config(
            SimulationConfig(; melting=MeltingConfig(; active=true, phi_crit=0.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; melting=MeltingConfig(; active=true, phi_crit=1.2))
        )

        # Validation bounds: unphysical alpha_eta
        @test_throws ArgumentError validate_config(
            SimulationConfig(; melting=MeltingConfig(; active=true, alpha_eta=-5.0))
        )

        # Validation bounds: unphysical Clapeyron slope
        @test_throws ArgumentError validate_config(
            SimulationConfig(; melting=MeltingConfig(; active=true, dpdt_clapeyron=-1.0e-7))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; melting=MeltingConfig(; active=true, dpdt_clapeyron=NaN))
        )

        # Validation bounds: invalid latent_heat_mode
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(; active=true, latent_heat_mode=:source_term)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(; active=true, latent_heat_mode=:invalid_mode)
            ),
        )

        # Validation bounds: soft_turbulence cannot be true when active=false
        @test_throws ArgumentError validate_config(
            SimulationConfig(; melting=MeltingConfig(; active=false, soft_turbulence=true))
        )
    end

    @testset "MeltingConfig TOML Serialization Round-Trip" begin
        custom_melting = MeltingConfig(;
            active=true,
            T_solidus=SVector{3,Float64}([1350.0, 1350.0, NaN]),
            T_liquidus=SVector{3,Float64}([1850.0, 1850.0, NaN]),
            L_melt=4.5e5,
            rho_melt=2750.0,
            alpha_eta=30.0,
            phi_crit=0.45,
            eta_melt=5.0,
            dpdt_clapeyron=1.3e-7,
            latent_heat_mode=:apparent_cp,
        )
        sim_cfg = SimulationConfig(; melting=custom_melting)
        toml_str = Erebus.Config.save_config(sim_cfg)
        loaded_cfg = Erebus.Config.load_config(toml_str)

        @test loaded_cfg.melting.active == true
        @test isapprox(loaded_cfg.melting.T_solidus[1], 1350.0; rtol=1e-12)
        @test isapprox(loaded_cfg.melting.T_liquidus[1], 1850.0; rtol=1e-12)
        @test isapprox(loaded_cfg.melting.L_melt, 4.5e5; rtol=1e-12)
        @test isapprox(loaded_cfg.melting.rho_melt, 2750.0; rtol=1e-12)
        @test isapprox(loaded_cfg.melting.alpha_eta, 30.0; rtol=1e-12)
        @test isapprox(loaded_cfg.melting.phi_crit, 0.45; rtol=1e-12)
        @test isapprox(loaded_cfg.melting.eta_melt, 5.0; rtol=1e-12)
        @test isapprox(loaded_cfg.melting.dpdt_clapeyron, 1.3e-7; rtol=1e-12)
        @test loaded_cfg.melting.latent_heat_mode == :apparent_cp
    end

    @testset "compute_melt_fraction: Analytical Limits & Monotonicity" begin
        T_sol = 1400.0
        T_liq = 1800.0

        # Sub-solidus: zero melt
        @test iszero(compute_melt_fraction(1200.0, 0.0, 1; T_sol=T_sol, T_liq=T_liq))
        @test iszero(compute_melt_fraction(T_sol, 0.0, 1; T_sol=T_sol, T_liq=T_liq))

        # Super-liquidus: complete melt
        @test isapprox(
            compute_melt_fraction(2000.0, 0.0, 1; T_sol=T_sol, T_liq=T_liq), 1.0; rtol=1e-12
        )
        @test isapprox(
            compute_melt_fraction(T_liq, 0.0, 1; T_sol=T_sol, T_liq=T_liq), 1.0; rtol=1e-12
        )

        # Midpoint: linear partition
        T_mid = 0.5 * (T_sol + T_liq)
        @test isapprox(
            compute_melt_fraction(T_mid, 0.0, 1; T_sol=T_sol, T_liq=T_liq), 0.5; rtol=1e-12
        )

        # Monotonicity test across interval
        T_vals = range(T_sol - 100.0, T_liq + 100.0; length=50)
        F_vals = [
            compute_melt_fraction(T, 0.0, 1; T_sol=T_sol, T_liq=T_liq) for T in T_vals
        ]
        @test issorted(F_vals)
        @test all(0.0 .<= F_vals .<= 1.0)

        # Sticky air (tm == 3): non-melting material
        @test iszero(compute_melt_fraction(2000.0, 0.0, 3; T_sol=T_sol, T_liq=T_liq))

        # Clapeyron pressure dependence: positive slope raises solidus and liquidus
        dpdt = 5.0e-8 # 50 K / GPa = 5e-8 K/Pa
        P_litho = 2.0e8 # 200 MPa
        F_P0 = compute_melt_fraction(1500.0, 0.0, 1; T_sol=T_sol, T_liq=T_liq, dpdt=dpdt)
        F_P = compute_melt_fraction(1500.0, P_litho, 1; T_sol=T_sol, T_liq=T_liq, dpdt=dpdt)
        @test F_P < F_P0
        T_s_eff = T_sol + dpdt * P_litho
        T_l_eff = T_liq + dpdt * P_litho
        expected_F_P = (1500.0 - T_s_eff) / (T_l_eff - T_s_eff)
        @test isapprox(F_P, expected_F_P; rtol=1e-12)

        # High pressure test: melting interval does not collapse at P > 400 MPa
        P_high = 1.0e9 # 1 GPa
        T_s_high = T_sol + dpdt * P_high # 1450 K
        T_l_high = T_liq + dpdt * P_high # 1850 K
        @test isapprox(T_l_high - T_s_high, T_liq - T_sol; rtol=1e-12)
        F_high = compute_melt_fraction(
            1650.0, P_high, 1; T_sol=T_sol, T_liq=T_liq, dpdt=dpdt
        )
        @test isapprox(F_high, (1650.0 - T_s_high) / (T_l_high - T_s_high); rtol=1e-12)

        # Discrimination guards: non-zero, non-unity, distinct from quadratic
        F_quarter = compute_melt_fraction(1500.0, 0.0, 1; T_sol=T_sol, T_liq=T_liq)
        @test isapprox(F_quarter, 0.25; rtol=1e-12)
        @test abs(F_quarter - 0.25^2) > 0.1 # Distinguishes linear from quadratic

        # Domain error contract: invalid and non-finite inputs
        @test_throws DomainError compute_melt_fraction(
            -10.0, 0.0, 1; T_sol=T_sol, T_liq=T_liq
        )
        @test_throws DomainError compute_melt_fraction(
            NaN, 0.0, 1; T_sol=T_sol, T_liq=T_liq
        )
        @test_throws DomainError compute_melt_fraction(
            1500.0, NaN, 1; T_sol=T_sol, T_liq=T_liq
        )
        @test_throws DomainError compute_melt_fraction(
            1500.0, 0.0, 1; T_sol=NaN, T_liq=T_liq
        )
        @test_throws DomainError compute_melt_fraction(
            1500.0, 0.0, 1; T_sol=T_sol, T_liq=Inf
        )
        @test_throws DomainError compute_melt_fraction(
            1500.0, 0.0, 1; T_sol=1800.0, T_liq=1400.0
        )
    end

    @testset "rhocp_apparent_silicate: Latent Heat & Conservation Integral" begin
        T_sol = 1400.0
        T_liq = 1800.0
        L_melt = 4.0e5
        rho_s = 3300.0
        rhocp_base = 3.3e6

        # Inactive mode: returns unaugmented heat capacity
        @test isapprox(
            rhocp_apparent_silicate(
                1600.0,
                0.0,
                rhocp_base,
                rho_s,
                1;
                T_sol=T_sol,
                T_liq=T_liq,
                L_melt=L_melt,
                active=false,
            ),
            rhocp_base;
            rtol=1e-12,
        )

        # Active mode below solidus: no apparent Cp change
        @test isapprox(
            rhocp_apparent_silicate(
                1200.0,
                0.0,
                rhocp_base,
                rho_s,
                1;
                T_sol=T_sol,
                T_liq=T_liq,
                L_melt=L_melt,
                active=true,
            ),
            rhocp_base;
            rtol=1e-12,
        )

        # Active mode above liquidus: no apparent Cp change
        @test isapprox(
            rhocp_apparent_silicate(
                1900.0,
                0.0,
                rhocp_base,
                rho_s,
                1;
                T_sol=T_sol,
                T_liq=T_liq,
                L_melt=L_melt,
                active=true,
            ),
            rhocp_base;
            rtol=1e-12,
        )

        # Active mode inside melting interval: effective Cp includes latent spike
        rhocp_eff = rhocp_apparent_silicate(
            1600.0,
            0.0,
            rhocp_base,
            rho_s,
            1;
            T_sol=T_sol,
            T_liq=T_liq,
            L_melt=L_melt,
            active=true,
        )
        expected_spike = rho_s * L_melt / (T_liq - T_sol)
        @test isapprox(rhocp_eff, rhocp_base + expected_spike; rtol=1e-12)
        @test rhocp_eff > rhocp_base

        # Energy conservation integral: ∫ (rhocp_eff - rhocp_base) dT == rho_s * L_melt
        N_pts = 1000
        T_grid = range(T_sol, T_liq; length=N_pts)
        dT = step(T_grid)
        excess_integral = sum(
            (
                rhocp_apparent_silicate(
                    T + 0.5 * dT,
                    0.0,
                    rhocp_base,
                    rho_s,
                    1;
                    T_sol=T_sol,
                    T_liq=T_liq,
                    L_melt=L_melt,
                    active=true,
                ) - rhocp_base
            ) * dT for T in range(T_sol, T_liq - dT; length=N_pts - 1)
        )
        @test isapprox(excess_integral, rho_s * L_melt; rtol=1e-6)

        # Sticky air: unaffected
        @test isapprox(
            rhocp_apparent_silicate(
                1600.0,
                0.0,
                3.0e6,
                1.0,
                3;
                T_sol=T_sol,
                T_liq=T_liq,
                L_melt=L_melt,
                active=true,
            ),
            3.0e6;
            rtol=1e-12,
        )

        # Clapeyron shift preserves latent heat buffering at high pressure
        dpdt = 5.0e-8
        P_high = 1.0e9 # 1 GPa
        rhocp_eff_high = rhocp_apparent_silicate(
            1650.0,
            P_high,
            rhocp_base,
            rho_s,
            1;
            T_sol=T_sol,
            T_liq=T_liq,
            L_melt=L_melt,
            active=true,
            dpdt=dpdt,
        )
        @test isapprox(rhocp_eff_high, rhocp_base + expected_spike; rtol=1e-12)

        # Domain error contract: negative and non-finite inputs
        @test_throws DomainError rhocp_apparent_silicate(
            -50.0,
            0.0,
            rhocp_base,
            rho_s,
            1;
            T_sol=T_sol,
            T_liq=T_liq,
            L_melt=L_melt,
            active=true,
        )
        @test_throws DomainError rhocp_apparent_silicate(
            NaN,
            0.0,
            rhocp_base,
            rho_s,
            1;
            T_sol=T_sol,
            T_liq=T_liq,
            L_melt=L_melt,
            active=true,
        )
        @test_throws DomainError rhocp_apparent_silicate(
            1500.0,
            NaN,
            rhocp_base,
            rho_s,
            1;
            T_sol=T_sol,
            T_liq=T_liq,
            L_melt=L_melt,
            active=true,
        )
        @test_throws DomainError rhocp_apparent_silicate(
            1500.0,
            0.0,
            rhocp_base,
            rho_s,
            1;
            T_sol=NaN,
            T_liq=T_liq,
            L_melt=L_melt,
            active=true,
        )
        @test_throws DomainError rhocp_apparent_silicate(
            1500.0,
            0.0,
            rhocp_base,
            rho_s,
            1;
            T_sol=T_sol,
            T_liq=Inf,
            L_melt=L_melt,
            active=true,
        )
        # Consistent invariant check: invalid solidus/liquidus throws for all tm and active states
        @test_throws DomainError rhocp_apparent_silicate(
            1500.0,
            0.0,
            rhocp_base,
            rho_s,
            3;
            T_sol=1800.0,
            T_liq=1400.0,
            L_melt=L_melt,
            active=true,
        )
        @test_throws DomainError rhocp_apparent_silicate(
            1500.0,
            0.0,
            rhocp_base,
            rho_s,
            1;
            T_sol=1800.0,
            T_liq=1400.0,
            L_melt=L_melt,
            active=false,
        )
    end

    @testset "compute_melt_weakened_viscosity: Rheological Transition" begin
        eta_s = 1.0e19
        alpha_eta = 28.0
        phi_crit = 0.4
        eta_melt = 10.0
        etamin = 1.0e12
        etamax = 1.0e23

        # Zero melt: returns unweakened solid viscosity
        @test isapprox(
            compute_melt_weakened_viscosity(
                eta_s,
                0.0,
                1;
                alpha_eta=alpha_eta,
                phi_crit=phi_crit,
                eta_melt=eta_melt,
                etamin=etamin,
                etamax=etamax,
            ),
            eta_s;
            rtol=1e-12,
        )

        # Sub-critical melt (F_m = 0.2): exponential weakening
        eta_02 = compute_melt_weakened_viscosity(
            eta_s,
            0.2,
            1;
            alpha_eta=alpha_eta,
            phi_crit=phi_crit,
            eta_melt=eta_melt,
            etamin=etamin,
            etamax=etamax,
        )
        expected_02 = eta_s * exp(-alpha_eta * 0.2)
        @test isapprox(eta_02, expected_02; rtol=1e-12)
        @test eta_02 < eta_s

        # Monotonicity: higher melt fraction produces lower viscosity
        eta_01 = compute_melt_weakened_viscosity(
            eta_s,
            0.1,
            1;
            alpha_eta=alpha_eta,
            phi_crit=phi_crit,
            eta_melt=eta_melt,
            etamin=etamin,
            etamax=etamax,
        )
        eta_03 = compute_melt_weakened_viscosity(
            eta_s,
            0.3,
            1;
            alpha_eta=alpha_eta,
            phi_crit=phi_crit,
            eta_melt=eta_melt,
            etamin=etamin,
            etamax=etamax,
        )
        @test eta_s > eta_01 > eta_02 > eta_03

        # Critical disaggregation (F_m >= phi_crit): drops to suspension/magma regime
        eta_crit = compute_melt_weakened_viscosity(
            eta_s,
            phi_crit,
            1;
            alpha_eta=alpha_eta,
            phi_crit=phi_crit,
            eta_melt=eta_melt,
            etamin=etamin,
            etamax=etamax,
        )
        eta_ocean = compute_melt_weakened_viscosity(
            eta_s,
            0.8,
            1;
            alpha_eta=alpha_eta,
            phi_crit=phi_crit,
            eta_melt=eta_melt,
            etamin=etamin,
            etamax=etamax,
        )
        @test eta_ocean <= eta_crit
        @test eta_ocean >= etamin # Clamped at etamin

        # Analytical value check in crystal suspension regime (unclamped by etamin)
        eta_at_crit = eta_s * exp(-alpha_eta * phi_crit)
        # At F_m = 0.7: midpoint between phi_crit = 0.4 and 1.0 (frac = 0.5)
        expected_07 = exp(0.5 * log(eta_at_crit) + 0.5 * log(eta_melt))
        eta_07_unclamped = compute_melt_weakened_viscosity(
            eta_s,
            0.7,
            1;
            alpha_eta=alpha_eta,
            phi_crit=phi_crit,
            eta_melt=eta_melt,
            etamin=1.0,
            etamax=etamax,
        )
        @test isapprox(eta_07_unclamped, expected_07; rtol=1e-12)

        # Pure melt limit at F_m = 1.0 reaches eta_melt
        eta_10_unclamped = compute_melt_weakened_viscosity(
            eta_s,
            1.0,
            1;
            alpha_eta=alpha_eta,
            phi_crit=phi_crit,
            eta_melt=eta_melt,
            etamin=1.0,
            etamax=etamax,
        )
        @test isapprox(eta_10_unclamped, eta_melt; rtol=1e-12)

        # Sticky air (tm == 3): unaffected, returns clamped eta_solid
        @test isapprox(
            compute_melt_weakened_viscosity(
                1.0e16,
                0.5,
                3;
                alpha_eta=alpha_eta,
                phi_crit=phi_crit,
                eta_melt=eta_melt,
                etamin=etamin,
                etamax=etamax,
            ),
            1.0e16;
            rtol=1e-12,
        )

        # 3-class discrimination guards
        @test abs(eta_02 - eta_s * exp(+alpha_eta * 0.2)) > 1e18 # Exponent sign guard
        @test eta_02 > 0.0                                       # Positivity
        @test etamin <= eta_02 <= etamax                         # Clamping bounds

        # Domain error contract: non-finite and non-positive inputs
        @test_throws DomainError compute_melt_weakened_viscosity(eta_s, NaN, 1)
        @test_throws DomainError compute_melt_weakened_viscosity(eta_s, Inf, 1)
        @test_throws DomainError compute_melt_weakened_viscosity(NaN, 0.5, 1)
        @test_throws DomainError compute_melt_weakened_viscosity(0.0, 0.5, 1)
        @test_throws DomainError compute_melt_weakened_viscosity(-1.0e19, 0.5, 1)
    end

    @testset "Marker Properties Setup with Fm" begin
        marknum = 100
        props = setup_marker_properties(marknum)
        @test length(props) == 13
        xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0, Fm =
            props

        @test length(Fm) == marknum
        @test eltype(Fm) == Float64
        @test all(iszero, Fm)
    end

    @testset "Marker Melting Coupling: compute_marker_properties! & update_marker_viscosity!" begin
        marknum = 3
        coords = GridCoordinates(GridConfig(; Nx=17, Ny=17, xsize=10_000.0, ysize=10_000.0))
        props = setup_marker_properties(marknum, coords)
        xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0, Fm =
            props
        helpers = setup_marker_properties_helpers(marknum)
        (
            rhototalm,
            rhocptotalm,
            etatotalm,
            hrtotalm,
            ktotalm,
            tkm_rhocptotalm,
            etafluidcur_inv_kphim,
            inv_gggtotalm,
            fricttotalm,
            cohestotalm,
            tenstotalm,
            rhofluidcur,
            alphasolidcur,
            alphafluidcur,
        ) = helpers

        tm .= 1
        phim .= 0.1
        # Marker 1: partially molten
        tkm[1] = 1600.0
        pfm0[1] = 2.0e7
        # Marker 2: sub-solidus
        tkm[2] = 1000.0
        pfm0[2] = 0.0
        # Marker 3: super-liquidus
        tkm[3] = 1900.0
        pfm0[3] = 0.0

        hrsolidm = SVector{3,Float64}([0.0, 0.0, 0.0])
        hrfluidm = SVector{3,Float64}([0.0, 0.0, 0.0])
        mode = 9

        # Execute compute_marker_properties! with melting active and Clapeyron slope
        dpdt_val = 1.0e-7 # 1e-7 K/Pa -> at 20 MPa, shift = 2 K
        for m in 1:marknum
            compute_marker_properties!(
                m,
                tm,
                tkm,
                rhototalm,
                rhocptotalm,
                etatotalm,
                hrtotalm,
                ktotalm,
                tkm_rhocptotalm,
                etafluidcur_inv_kphim,
                hrsolidm,
                hrfluidm,
                phim,
                XWsolidm0,
                mode,
                rhofluidcur;
                pm=pfm0,
                Fm=Fm,
                melting_active=true,
                T_solidus_val=SVector{3,Float64}([1400.0, 1400.0, NaN]),
                T_liquidus_val=SVector{3,Float64}([1800.0, 1800.0, NaN]),
                L_melt_val=4.0e5,
                rho_melt_val=2800.0,
                alpha_eta_val=28.0,
                phi_crit_val=0.4,
                eta_melt_val=10.0,
                dpdt_clapeyron_val=dpdt_val,
            )
        end

        # Marker 1: T_s = 1400 + 1e-7 * 2e7 = 1402 K, T_l = 1800 + 1e-7 * 2e7 = 1802 K
        # F_m = (1600 - 1402) / (1802 - 1402) = 198 / 400 = 0.495
        expected_Fm1 = (1600.0 - 1402.0) / (1802.0 - 1402.0)
        @test isapprox(Fm[1], expected_Fm1; rtol=1e-12)
        # Marker 2: Sub-solidus -> F_m = 0
        @test iszero(Fm[2])
        # Marker 3: Super-liquidus -> F_m = 1
        @test isapprox(Fm[3], 1.0; rtol=1e-12)

        # Apparent heat capacity on marker 1 includes latent heat
        expected_rhocp1 = rhocp_apparent_silicate(
            1600.0,
            2.0e7,
            Erebus.rhocpsolidm[1],
            Erebus.rhosolidm[1],
            1;
            T_sol=1400.0,
            T_liq=1800.0,
            L_melt=4.0e5,
            active=true,
            dpdt=dpdt_val,
        )
        expected_total_rhocp1 = total(
            expected_rhocp1, Erebus.compute_rhocpfluidm(1600.0, mode), phim[1]
        )
        @test isapprox(rhocptotalm[1], expected_total_rhocp1; rtol=1e-12)

        # Execute update_marker_viscosity! with melting active
        YNY = zeros(coords.Ny, coords.Nx)
        YNY_inv_ETA = zeros(coords.Ny, coords.Nx)
        for m in 1:marknum
            update_marker_viscosity!(
                m,
                xm,
                ym,
                tm,
                tkm,
                etatotalm,
                etavpm,
                YNY,
                YNY_inv_ETA;
                coords=coords,
                Fm=Fm,
                melting_active=true,
                alpha_eta_val=28.0,
                phi_crit_val=0.4,
                eta_melt_val=10.0,
            )
        end

        # Partially molten marker 1 viscosity is weakened relative to sub-solidus marker 2
        @test etatotalm[1] < etatotalm[2]
        expected_eta1 = compute_melt_weakened_viscosity(
            Erebus.etatotal_rocks(tkm[1], tm[1]),
            Fm[1],
            tm[1];
            alpha_eta=28.0,
            phi_crit=0.4,
            eta_melt=10.0,
        )
        @test isapprox(etatotalm[1], expected_eta1; rtol=1e-12)
    end
end

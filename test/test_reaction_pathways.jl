using Erebus
using Erebus.Config
using Erebus.Physics
using Erebus.Numerics
using Erebus.Geometry
using Erebus.Particles
using Test
using TOML

_mean(x) = sum(x) / length(x)

@testset "Reaction Pathways & Thermodynamic Coupling" begin
    @testset "ReactionConfig Schema & Bounds Validation" begin
        cfg_def = ReactionConfig()
        @test cfg_def.active == true
        @test cfg_def.hydration_active == true
        @test cfg_def.dehydration_active == true
        @test cfg_def.hydration_mode == 1
        @test cfg_def.dehydration_mode == 2
        @test isapprox(cfg_def.delta_H, 40_000.0; rtol=1e-12)
        @test isapprox(cfg_def.delta_S, 60.0; rtol=1e-12)
        @test isapprox(cfg_def.c_I, 543.0; rtol=1e-12)
        @test isapprox(cfg_def.alpha_relaxation, 0.5; rtol=1e-12)
        @test isapprox(cfg_def.p_cavitation, 1.0e7; rtol=1e-12)

        # SimulationConfig inclusion
        sim_cfg = default_config()
        @test sim_cfg.reaction isa ReactionConfig
        @test validate_config(sim_cfg) === nothing

        # Validation error cases
        @test_throws ArgumentError validate_config(
            SimulationConfig(; reaction=ReactionConfig(; hydration_mode=5))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; reaction=ReactionConfig(; dehydration_mode=0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; reaction=ReactionConfig(; delta_H=-1000.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; reaction=ReactionConfig(; delta_S=0.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; reaction=ReactionConfig(; alpha_relaxation=0.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; reaction=ReactionConfig(; alpha_relaxation=1.5))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; reaction=ReactionConfig(; pfcoeff=-0.1))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(; reaction=ReactionConfig(; p_cavitation=-1.0e6))
        )
    end

    @testset "ReactionConfig TOML Serialization Round-Trip" begin
        sim_cfg = default_config()
        custom_react = ReactionConfig(;
            hydration_active=false,
            dtreaction_hydration=2.0e10,
            dtreaction_dehydration=5.0e7,
            delta_H=42_000.0,
            delta_S=62.0,
            alpha_relaxation=0.75,
            p_cavitation=2.0e7,
        )
        sim_cfg_custom = SimulationConfig(; reaction=custom_react)
        toml_str = sprint(save_config, sim_cfg_custom)
        @test occursin("[reaction]", toml_str)
        @test occursin("hydration_active = false", toml_str)

        loaded_cfg = load_config(toml_str)
        @test loaded_cfg.reaction.hydration_active == false
        @test isapprox(loaded_cfg.reaction.dtreaction_hydration, 2.0e10; rtol=1e-12)
        @test isapprox(loaded_cfg.reaction.dtreaction_dehydration, 5.0e7; rtol=1e-12)
        @test isapprox(loaded_cfg.reaction.delta_H, 42_000.0; rtol=1e-12)
        @test isapprox(loaded_cfg.reaction.delta_S, 62.0; rtol=1e-12)
        @test isapprox(loaded_cfg.reaction.alpha_relaxation, 0.75; rtol=1e-12)
        @test isapprox(loaded_cfg.reaction.p_cavitation, 2.0e7; rtol=1e-12)
    end

    @testset "Equilibrium Direction and Continuous Phase Boundary" begin
        # Clausius-Clapeyron equilibrium temperature: Teq = (ΔH + Pf*ΔV) / ΔS
        # At Pf = 0: Teq = 40000 / 60 = 666.67 K
        # At Pf = 30 MPa (3.0e7 Pa): Teq = (40000 + 3e7 * 1.2867e-6) / 60 ≈ 670.64 K
        Pf = 3.0e7
        T_cold = 500.0
        T_hot = 750.0
        XW_init = 0.5

        # Cold regime: Gibbs free energy favors hydrated serpentine (forward reaction)
        dG_cold_eq = compute_gibbs_free_energy(
            T_cold, Pf, 1.0 - XW_init, XW_init, 1.0e11, 1.0e10
        )
        K_cold = compute_reaction_constant(T_cold, Pf, dG_cold_eq)
        XW_eq_cold = inv(K_cold + 1.0)
        @test XW_eq_cold > XW_init
        @test XW_eq_cold > 0.90

        # Hot regime: Gibbs free energy favors dry silicate + steam (reverse reaction)
        dG_hot_eq = compute_gibbs_free_energy(
            T_hot, Pf, 1.0 - XW_init, XW_init, 1.0e11, 1.0e10
        )
        K_hot = compute_reaction_constant(T_hot, Pf, dG_hot_eq)
        XW_eq_hot = inv(K_hot + 1.0)
        @test XW_eq_hot < XW_init
        @test XW_eq_hot < 0.35

        # Finite-rate relaxation at Δt < Δtr
        dG_cold_half = compute_gibbs_free_energy(
            T_cold, Pf, 1.0 - XW_init, XW_init, 5.0e9, 1.0e10
        )
        K_half = compute_reaction_constant(T_cold, Pf, dG_cold_half)
        XW_half = inv(K_half + 1.0)
        @test XW_init < XW_half < XW_eq_cold
    end

    @testset "Two-Way Kinetics Timescales with ReactionConfig" begin
        cfg = ReactionConfig(;
            c_I=543.0, b_I=2.5e-4, A_I=1.0e-11, To_B=293.0, Tscl_B=10.0, Sxo_B=2.0e-11
        )
        phi = 0.2

        # Mode 1 Gaussian (hydration) has minimum timescale (fastest rate) at T = c_I
        tr_opt = compute_Δtreaction(543.0, phi, 1; cfg=cfg)
        tr_sub = compute_Δtreaction(443.0, phi, 1; cfg=cfg)
        tr_sup = compute_Δtreaction(643.0, phi, 1; cfg=cfg)
        @test tr_opt < tr_sub
        @test tr_opt < tr_sup
        @test isapprox(tr_sub, tr_sup; rtol=1e-8)

        # Mode 2 Pseudo-Arrhenius (dehydration) decreases monotonically with temperature
        tr_dehyd_cool = compute_Δtreaction(650.0, phi, 2; cfg=cfg)
        tr_dehyd_hot = compute_Δtreaction(750.0, phi, 2; cfg=cfg)
        @test tr_dehyd_hot < tr_dehyd_cool
        @test tr_dehyd_hot > 0.0

        # Mode 9 Constant Timescale with explicit direction selection
        tr_const_hyd = compute_Δtreaction(500.0, phi, 9; cfg=cfg, is_hydration=true)
        tr_const_deh = compute_Δtreaction(700.0, phi, 9; cfg=cfg, is_hydration=false)
        @test isapprox(tr_const_hyd, cfg.dtreaction_hydration; rtol=1e-12)
        @test isapprox(tr_const_deh, cfg.dtreaction_dehydration; rtol=1e-12)
    end

    @testset "Exothermic Hydration: Physical Invariants & Fluid Suction" begin
        c17 = GridCoordinates(GridConfig(; Nx=17, Ny=17, xsize=10_000.0, ysize=10_000.0))
        marknum = 100
        props = setup_marker_properties(marknum, c17)
        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0) =
            props

        tm .= 1
        xm .= 5000.0
        ym .= 5000.0
        phim .= 0.20
        phinewm .= 0.20
        pfm0 .= 1.0e7

        # Initially anhydrous silicate
        XWsolidm0 .= 0.05
        XWsolidm .= 0.05

        tk2 = fill(500.0, c17.Ny1, c17.Nx1)
        pf = fill(1.0e7, c17.Ny1, c17.Nx1)

        DMP = zeros(Float64, c17.Ny1, c17.Nx1)
        DHP = zeros(Float64, c17.Ny1, c17.Nx1)
        DQPF = zeros(Float64, c17.Ny1, c17.Nx1)
        DMPSUM = zeros(Float64, c17.Ny1, c17.Nx1)
        DHPSUM = zeros(Float64, c17.Ny1, c17.Nx1)
        DQPFSUM = zeros(Float64, c17.Ny1, c17.Nx1)
        WTPSUM = zeros(Float64, c17.Ny1, c17.Nx1)

        dt = 1.0e8

        perform_thermochemical_reaction!(
            DMP,
            DHP,
            DMPSUM,
            DHPSUM,
            WTPSUM,
            pf,
            tk2,
            tm,
            xm,
            ym,
            XWsolidm0,
            XWsolidm,
            phim,
            phinewm,
            pfm0,
            marknum,
            dt,
            2,
            1;
            coords=c17,
            DQPF=DQPF,
            DQPFSUM=DQPFSUM,
        )

        # 1. Hydration progression: wet silicate fraction increases
        @test _mean(XWsolidm) > 0.05

        # 2. Porosity decreases
        @test _mean(phinewm) < 0.20

        # 3. Latent heat release: DHP > 0 (exothermic heating)
        max_dhp = maximum(DHP)
        @test max_dhp > 0.0
        @test minimum(DHP) >= -1e-15

        # 4. Fluid suction: DQPF < 0 (water consumed from pore space)
        min_dqpf = minimum(DQPF)
        @test min_dqpf < 0.0
        @test maximum(DQPF) <= 1e-15

        # 5. Closed-box exact water mass conservation
        MH2O_mol = 0.018
        MD_mol = 0.120
        delta_m_mineral_water = sum(
            begin
                rho_s0 =
                    (MD_mol + MH2O_mol * XWsolidm0[m]) /
                    (Erebus.VDˢ * (1.0 - XWsolidm0[m]) + Erebus.VWˢ * XWsolidm0[m])
                rho_f0 = 1000.0
                rho0 = rho_s0 * (1.0 - phim[m]) + rho_f0 * phim[m]

                rho_s1 =
                    (MD_mol + MH2O_mol * XWsolidm[m]) /
                    (Erebus.VDˢ * (1.0 - XWsolidm[m]) + Erebus.VWˢ * XWsolidm[m])
                rho_f1 = 1000.0
                rho1 = rho_s1 * (1.0 - phinewm[m]) + rho_f1 * phinewm[m]
                RV = rho1 / rho0

                mw1 =
                    (1.0 - phinewm[m]) * (1.0 / RV) * (XWsolidm[m] * MH2O_mol) /
                    (Erebus.VDˢ * (1.0 - XWsolidm[m]) + Erebus.VWˢ * XWsolidm[m])
                mw0 =
                    (1.0 - phim[m]) * (XWsolidm0[m] * MH2O_mol) /
                    (Erebus.VDˢ * (1.0 - XWsolidm0[m]) + Erebus.VWˢ * XWsolidm0[m])
                mw1 - mw0
            end for m in 1:marknum
        )
        delta_m_pore_water = sum(
            begin
                rho_s0 =
                    (MD_mol + MH2O_mol * XWsolidm0[m]) /
                    (Erebus.VDˢ * (1.0 - XWsolidm0[m]) + Erebus.VWˢ * XWsolidm0[m])
                rho_f0 = 1000.0
                rho0 = rho_s0 * (1.0 - phim[m]) + rho_f0 * phim[m]

                rho_s1 =
                    (MD_mol + MH2O_mol * XWsolidm[m]) /
                    (Erebus.VDˢ * (1.0 - XWsolidm[m]) + Erebus.VWˢ * XWsolidm[m])
                rho_f1 = 1000.0
                rho1 = rho_s1 * (1.0 - phinewm[m]) + rho_f1 * phinewm[m]
                RV = rho1 / rho0

                1000.0 * (phinewm[m] * (1.0 / RV) - phim[m])
            end for m in 1:marknum
        )
        @test isapprox(delta_m_mineral_water + delta_m_pore_water, 0.0; atol=1e-10)
    end

    @testset "Endothermic Dehydration: Physical Invariants & Pore Overpressure" begin
        c17 = GridCoordinates(GridConfig(; Nx=17, Ny=17, xsize=10_000.0, ysize=10_000.0))
        marknum = 100
        props = setup_marker_properties(marknum, c17)
        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0) =
            props

        tm .= 1
        xm .= 5000.0
        ym .= 5000.0
        phim .= 0.10
        phinewm .= 0.10
        pfm0 .= 3.0e7

        # Initially hydrous serpentine
        XWsolidm0 .= 0.95
        XWsolidm .= 0.95

        tk2 = fill(750.0, c17.Ny1, c17.Nx1)
        pf = fill(3.0e7, c17.Ny1, c17.Nx1)

        DMP = zeros(Float64, c17.Ny1, c17.Nx1)
        DHP = zeros(Float64, c17.Ny1, c17.Nx1)
        DQPF = zeros(Float64, c17.Ny1, c17.Nx1)
        DMPSUM = zeros(Float64, c17.Ny1, c17.Nx1)
        DHPSUM = zeros(Float64, c17.Ny1, c17.Nx1)
        DQPFSUM = zeros(Float64, c17.Ny1, c17.Nx1)
        WTPSUM = zeros(Float64, c17.Ny1, c17.Nx1)

        dt = 1.0e6

        perform_thermochemical_reaction!(
            DMP,
            DHP,
            DMPSUM,
            DHPSUM,
            WTPSUM,
            pf,
            tk2,
            tm,
            xm,
            ym,
            XWsolidm0,
            XWsolidm,
            phim,
            phinewm,
            pfm0,
            marknum,
            dt,
            2,
            1;
            coords=c17,
            DQPF=DQPF,
            DQPFSUM=DQPFSUM,
        )

        # 1. Dehydration progression: wet silicate fraction decreases
        @test _mean(XWsolidm) < 0.95

        # 2. Porosity increases
        @test _mean(phinewm) > 0.10

        # 3. Latent heat sink: DHP < 0 (endothermic cooling)
        min_dhp = minimum(DHP)
        @test min_dhp < 0.0
        @test maximum(DHP) <= 1e-15

        # 4. Fluid overpressure: DQPF > 0
        max_dqpf = maximum(DQPF)
        @test max_dqpf > 0.0
        @test minimum(DQPF) >= -1e-15

        # 5. Closed-box exact water mass conservation
        MH2O_mol = 0.018
        MD_mol = 0.120
        delta_m_mineral_water = sum(
            begin
                rho_s0 =
                    (MD_mol + MH2O_mol * XWsolidm0[m]) /
                    (Erebus.VDˢ * (1.0 - XWsolidm0[m]) + Erebus.VWˢ * XWsolidm0[m])
                rho_f0 = 1000.0
                rho0 = rho_s0 * (1.0 - phim[m]) + rho_f0 * phim[m]

                rho_s1 =
                    (MD_mol + MH2O_mol * XWsolidm[m]) /
                    (Erebus.VDˢ * (1.0 - XWsolidm[m]) + Erebus.VWˢ * XWsolidm[m])
                rho_f1 = 1000.0
                rho1 = rho_s1 * (1.0 - phinewm[m]) + rho_f1 * phinewm[m]
                RV = rho1 / rho0

                mw1 =
                    (1.0 - phinewm[m]) * (1.0 / RV) * (XWsolidm[m] * MH2O_mol) /
                    (Erebus.VDˢ * (1.0 - XWsolidm[m]) + Erebus.VWˢ * XWsolidm[m])
                mw0 =
                    (1.0 - phim[m]) * (XWsolidm0[m] * MH2O_mol) /
                    (Erebus.VDˢ * (1.0 - XWsolidm0[m]) + Erebus.VWˢ * XWsolidm0[m])
                mw1 - mw0
            end for m in 1:marknum
        )
        delta_m_pore_water = sum(
            begin
                rho_s0 =
                    (MD_mol + MH2O_mol * XWsolidm0[m]) /
                    (Erebus.VDˢ * (1.0 - XWsolidm0[m]) + Erebus.VWˢ * XWsolidm0[m])
                rho_f0 = 1000.0
                rho0 = rho_s0 * (1.0 - phim[m]) + rho_f0 * phim[m]

                rho_s1 =
                    (MD_mol + MH2O_mol * XWsolidm[m]) /
                    (Erebus.VDˢ * (1.0 - XWsolidm[m]) + Erebus.VWˢ * XWsolidm[m])
                rho_f1 = 1000.0
                rho1 = rho_s1 * (1.0 - phinewm[m]) + rho_f1 * phinewm[m]
                RV = rho1 / rho0

                1000.0 * (phinewm[m] * (1.0 / RV) - phim[m])
            end for m in 1:marknum
        )
        @test isapprox(delta_m_mineral_water + delta_m_pore_water, 0.0; atol=1e-10)
    end

    @testset "Dynamic Hydrofracture Coupling to Fluid Overpressure" begin
        kphi = 1.0e-15
        sigma_t = 1.0e7

        # Sub-critical effective pressure: Pt - Pf = 2.0e7 > -sigma_t -> no fracture
        keff_sub = compute_hydrofracture_permeability(
            kphi, 2.0e7, sigma_t; active=true, kappa_frac=1.0e3, gamma=1.0, kmax=1.0e-9
        )
        @test isapprox(keff_sub, kphi; rtol=1e-10)

        # Super-critical overpressure: Pt - Pf = -2.0e7 <= -sigma_t -> fracture active
        keff_frac = compute_hydrofracture_permeability(
            kphi, -2.0e7, sigma_t; active=true, kappa_frac=1.0e3, gamma=1.0, kmax=1.0e-9
        )
        @test keff_frac > kphi
        @test keff_frac <= 1.0e-9
    end

    @testset "Reaction Activation Switches and Picard Under-Relaxation" begin
        c17 = GridCoordinates(GridConfig(; Nx=17, Ny=17, xsize=10_000.0, ysize=10_000.0))
        marknum = 10
        props = setup_marker_properties(marknum, c17)
        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0) =
            props

        tm .= 1
        xm .= 5000.0
        ym .= 5000.0
        phim .= 0.20
        phinewm .= 0.20
        pfm0 .= 1.0e7
        XWsolidm0 .= 0.05
        XWsolidm .= 0.05

        tk2 = fill(500.0, c17.Ny1, c17.Nx1)
        pf = fill(1.0e7, c17.Ny1, c17.Nx1)

        DMP = zeros(Float64, c17.Ny1, c17.Nx1)
        DHP = zeros(Float64, c17.Ny1, c17.Nx1)
        DQPF = zeros(Float64, c17.Ny1, c17.Nx1)
        DMPSUM = zeros(Float64, c17.Ny1, c17.Nx1)
        DHPSUM = zeros(Float64, c17.Ny1, c17.Nx1)
        DQPFSUM = zeros(Float64, c17.Ny1, c17.Nx1)
        WTPSUM = zeros(Float64, c17.Ny1, c17.Nx1)

        # 1. Total deactivation switch
        cfg_inactive = ReactionConfig(; active=false)
        perform_thermochemical_reaction!(
            DMP,
            DHP,
            DMPSUM,
            DHPSUM,
            WTPSUM,
            pf,
            tk2,
            tm,
            xm,
            ym,
            XWsolidm0,
            XWsolidm,
            phim,
            phinewm,
            pfm0,
            marknum,
            1.0e8,
            2,
            1;
            coords=c17,
            DQPF=DQPF,
            DQPFSUM=DQPFSUM,
            cfg=cfg_inactive,
        )
        @test all(iszero, DMP)
        @test all(iszero, DHP)
        @test all(iszero, DQPF)

        # 2. Hydration deactivation switch
        cfg_no_hyd = ReactionConfig(; hydration_active=false)
        perform_thermochemical_reaction!(
            DMP,
            DHP,
            DMPSUM,
            DHPSUM,
            WTPSUM,
            pf,
            tk2,
            tm,
            xm,
            ym,
            XWsolidm0,
            XWsolidm,
            phim,
            phinewm,
            pfm0,
            marknum,
            1.0e8,
            2,
            1;
            coords=c17,
            DQPF=DQPF,
            DQPFSUM=DQPFSUM,
            cfg=cfg_no_hyd,
        )
        @test isapprox(_mean(XWsolidm), 0.05; atol=1e-12)

        # 3. Picard relaxation midpoint check
        XW_prev = 0.10
        XWsolidm .= XW_prev
        XWsolidm0 .= 0.05
        cfg_unrelaxed = ReactionConfig(; alpha_relaxation=1.0)
        perform_thermochemical_reaction!(
            DMP,
            DHP,
            DMPSUM,
            DHPSUM,
            WTPSUM,
            pf,
            tk2,
            tm,
            xm,
            ym,
            XWsolidm0,
            XWsolidm,
            phim,
            phinewm,
            pfm0,
            marknum,
            1.0e8,
            2,
            2;
            coords=c17,
            DQPF=DQPF,
            DQPFSUM=DQPFSUM,
            cfg=cfg_unrelaxed,
        )
        XW_unrelaxed = _mean(XWsolidm)

        # Run relaxed with alpha = 0.5 from same initial state
        XWsolidm .= XW_prev
        XWsolidm0 .= 0.05
        cfg_relax = ReactionConfig(; alpha_relaxation=0.5)
        perform_thermochemical_reaction!(
            DMP,
            DHP,
            DMPSUM,
            DHPSUM,
            WTPSUM,
            pf,
            tk2,
            tm,
            xm,
            ym,
            XWsolidm0,
            XWsolidm,
            phim,
            phinewm,
            pfm0,
            marknum,
            1.0e8,
            2,
            2;
            coords=c17,
            DQPF=DQPF,
            DQPFSUM=DQPFSUM,
            cfg=cfg_relax,
        )
        XW_relaxed = _mean(XWsolidm)
        @test isapprox(XW_relaxed, 0.5 * XW_unrelaxed + 0.5 * XW_prev; rtol=1e-10)

        # 4. Dehydration deactivation switch
        tk_hot = fill(750.0, c17.Ny1, c17.Nx1)
        XWsolidm0 .= 0.95
        XWsolidm .= 0.95
        cfg_no_deh = ReactionConfig(; dehydration_active=false)
        perform_thermochemical_reaction!(
            DMP,
            DHP,
            DMPSUM,
            DHPSUM,
            WTPSUM,
            pf,
            tk_hot,
            tm,
            xm,
            ym,
            XWsolidm0,
            XWsolidm,
            phim,
            phinewm,
            pfm0,
            marknum,
            1.0e8,
            2,
            1;
            coords=c17,
            DQPF=DQPF,
            DQPFSUM=DQPFSUM,
            cfg=cfg_no_deh,
        )
        @test isapprox(_mean(XWsolidm), 0.95; atol=1e-12)

        # 5. Timestep 1 backload branch (mutates XWsolidm0 and phim to match markers)
        XWsolidm .= 0.50
        XWsolidm0 .= 0.10
        phinewm .= 0.25
        phim .= 0.10
        perform_thermochemical_reaction!(
            DMP,
            DHP,
            DMPSUM,
            DHPSUM,
            WTPSUM,
            pf,
            tk2,
            tm,
            xm,
            ym,
            XWsolidm0,
            XWsolidm,
            phim,
            phinewm,
            pfm0,
            marknum,
            1.0e8,
            1,
            1;
            coords=c17,
            DQPF=DQPF,
            DQPFSUM=DQPFSUM,
        )
        @test isapprox(_mean(XWsolidm0), _mean(XWsolidm); rtol=1e-10)
        @test isapprox(_mean(phim), _mean(phinewm); rtol=1e-10)

        # 6. Cavitation floor clamping
        pf_neg = fill(-5.0e7, c17.Ny1, c17.Nx1)
        cfg_cav = ReactionConfig(; p_cavitation=1.0e7)
        perform_thermochemical_reaction!(
            DMP,
            DHP,
            DMPSUM,
            DHPSUM,
            WTPSUM,
            pf_neg,
            tk2,
            tm,
            xm,
            ym,
            XWsolidm0,
            XWsolidm,
            phim,
            phinewm,
            pfm0,
            marknum,
            1.0e8,
            2,
            1;
            coords=c17,
            DQPF=DQPF,
            DQPFSUM=DQPFSUM,
            cfg=cfg_cav,
        )
        @test all(pfm0 .>= -1.0e7)

        # 7. compute_thermochemical_iteration_outcome with non-default pferrmax
        dmp_dummy = zeros(5, 5)
        pf_a = fill(1.0e6, 5, 5)
        pf_b = fill(1.0e6 + 50.0, 5, 5)
        @test compute_thermochemical_iteration_outcome(
            dmp_dummy, pf_a, pf_b, 3; pferrmax=100.0
        )
        @test !compute_thermochemical_iteration_outcome(
            dmp_dummy, pf_a, pf_b, 3; pferrmax=10.0
        )

        # 8. DomainError guard on NaN and out-of-bounds marker wet fraction
        XW_nan = copy(XWsolidm0)
        XW_nan[1] = NaN
        @test_throws DomainError perform_thermochemical_reaction!(
            DMP,
            DHP,
            DMPSUM,
            DHPSUM,
            WTPSUM,
            pf,
            tk2,
            tm,
            xm,
            ym,
            XW_nan,
            XWsolidm,
            phim,
            phinewm,
            pfm0,
            marknum,
            1.0e8,
            2,
            1;
            coords=c17,
        )
        XW_neg = copy(XWsolidm0)
        XW_neg[1] = -0.1
        @test_throws DomainError perform_thermochemical_reaction!(
            DMP,
            DHP,
            DMPSUM,
            DHPSUM,
            WTPSUM,
            pf,
            tk2,
            tm,
            xm,
            ym,
            XW_neg,
            XWsolidm,
            phim,
            phinewm,
            pfm0,
            marknum,
            1.0e8,
            2,
            1;
            coords=c17,
        )
        XW_excess = copy(XWsolidm0)
        XW_excess[1] = 1.2
        @test_throws DomainError perform_thermochemical_reaction!(
            DMP,
            DHP,
            DMPSUM,
            DHPSUM,
            WTPSUM,
            pf,
            tk2,
            tm,
            xm,
            ym,
            XW_excess,
            XWsolidm,
            phim,
            phinewm,
            pfm0,
            marknum,
            1.0e8,
            2,
            1;
            coords=c17,
        )

        # 9. Large-extent hydration step with fast kinetics
        cfg_fast_hyd = ReactionConfig(; dtreaction_hydration=1.0e5, hydration_mode=9)
        XW_init_hyd = fill(0.05, marknum)
        XW_out_hyd = fill(0.05, marknum)
        phi_init_hyd = fill(0.50, marknum)
        phi_out_hyd = fill(0.50, marknum)
        pfm0_hyd = fill(3.0e7, marknum)
        perform_thermochemical_reaction!(
            DMP,
            DHP,
            DMPSUM,
            DHPSUM,
            WTPSUM,
            pf,
            tk2,
            tm,
            xm,
            ym,
            XW_init_hyd,
            XW_out_hyd,
            phi_init_hyd,
            phi_out_hyd,
            pfm0_hyd,
            marknum,
            1.0e7,
            2,
            1;
            coords=c17,
            DQPF=DQPF,
            DQPFSUM=DQPFSUM,
            cfg=cfg_fast_hyd,
        )
        @test _mean(XW_out_hyd) > 0.85
        @test _mean(phi_out_hyd) < 0.35
        @test any(DHP .> 0.0)
        @test any(DQPF .< 0.0)
    end
end

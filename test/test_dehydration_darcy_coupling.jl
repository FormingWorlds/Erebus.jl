using Test
using LinearAlgebra
using SparseArrays
using Random
using Erebus
using Erebus.Config
using Erebus.Numerics
using Erebus.Physics
using Erebus.Particles
using TOML

Random.seed!(42)

@testset "Dehydration-Darcy Fluid Overpressure & Venting Speciation" begin
    @testset "ReactionConfig Schema & Serialization" begin
        cfg_def = ReactionConfig()
        @test cfg_def.fluid_overpressure_coupling == true

        cfg_false = ReactionConfig(fluid_overpressure_coupling=false)
        @test cfg_false.fluid_overpressure_coupling == false

        sim_cfg = default_config()
        @test sim_cfg.reaction.fluid_overpressure_coupling == true
        @test validate_config(sim_cfg) === nothing

        # TOML round-trip
        toml_str = sprint(save_config, sim_cfg)
        @test occursin("fluid_overpressure_coupling = true", toml_str)
        mktemp() do toml_path, io
            save_config(toml_path, sim_cfg)
            cfg_loaded = load_config(toml_path)
            @test cfg_loaded.reaction.fluid_overpressure_coupling == true
        end

        sim_cfg_false = SimulationConfig(
            reaction=ReactionConfig(fluid_overpressure_coupling=false)
        )
        mktemp() do toml_path, io
            save_config(toml_path, sim_cfg_false)
            cfg_loaded_false = load_config(toml_path)
            @test cfg_loaded_false.reaction.fluid_overpressure_coupling == false
        end
    end

    @testset "Analytical Dehydration Overpressure Pulse (Zero-Permeability Limit)" begin
        coords = GridCoordinates(GridConfig(Nx=7, Ny=7, xsize=35000.0, ysize=35000.0))
        Ny, Nx = coords.Ny, coords.Nx
        Ny1, Nx1 = coords.Ny1, coords.Nx1

        eta_solid = 1.0e22
        eta_bulk = 1.0e25
        beta_drained = 1.0e-10
        beta_fluid = 1.0e-9
        phi_val = 0.1
        dt_val = 1.0e8

        ETA = fill(eta_solid, Ny, Nx)
        ETAP = fill(eta_solid, Ny1, Nx1)
        GGG = fill(1.0e10, Ny, Nx)
        GGGP = fill(1.0e10, Ny1, Nx1)
        SXY0 = zeros(Ny, Nx)
        SXX0 = zeros(Ny, Nx)
        RHOX = fill(3000.0, Ny1, Nx1)
        RHOY = fill(3000.0, Ny1, Nx1)
        RHOFX = fill(1000.0, Ny1, Nx1)
        RHOFY = fill(1000.0, Ny1, Nx1)
        # Impermeable domain: very high Darcy flow resistance
        RX = fill(1.0e20, Ny1, Nx1)
        RY = fill(1.0e20, Ny1, Nx1)
        PHI = fill(phi_val, Ny1, Nx1)
        ETAPHI = fill(eta_bulk, Ny1, Nx1)
        BETAPHI = fill(beta_drained, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = zeros(Ny1, Nx1)
        pr0 = zeros(Ny1, Nx1)
        pf0 = zeros(Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)

        i_int, j_int = 4, 4
        dqpf_val = 2.0e-12 # [1/s] fluid volume production rate

        # 1. Solve with fluid_overpressure_coupling = true
        DQPF_on = zeros(Ny1, Nx1)
        DQPF_on[i_int, j_int] = dqpf_val
        R_6_on = zeros(Nx1 * Ny1 * 6)

        L_6_on = assemble_hydromechanical_lse!(
            ETA,
            ETAP,
            GGG,
            GGGP,
            SXY0,
            SXX0,
            RHOX,
            RHOY,
            RHOFX,
            RHOFY,
            RX,
            RY,
            ETAPHI,
            BETAPHI,
            PHI,
            gx,
            gy,
            pr0,
            pf0,
            DMP,
            dt_val,
            R_6_on;
            coords=coords,
            DQPF=DQPF_on,
            fluid_overpressure_coupling=true,
        )
        S_6_on = L_6_on \ R_6_on
        pr_6_on = zeros(Ny1, Nx1)
        pf_6_on = zeros(Ny1, Nx1)
        vx_6_on = zeros(Ny1, Nx1)
        vy_6_on = zeros(Ny1, Nx1)
        qx_6_on = zeros(Ny1, Nx1)
        qy_6_on = zeros(Ny1, Nx1)
        process_hydromechanical_solution!(
            S_6_on, vx_6_on, vy_6_on, pr_6_on, qx_6_on, qy_6_on, pf_6_on; coords=coords
        )

        # 2. Solve with fluid_overpressure_coupling = false
        R_6_off = zeros(Nx1 * Ny1 * 6)
        L_6_off = assemble_hydromechanical_lse!(
            ETA,
            ETAP,
            GGG,
            GGGP,
            SXY0,
            SXX0,
            RHOX,
            RHOY,
            RHOFX,
            RHOFY,
            RX,
            RY,
            ETAPHI,
            BETAPHI,
            PHI,
            gx,
            gy,
            pr0,
            pf0,
            DMP,
            dt_val,
            R_6_off;
            coords=coords,
            DQPF=DQPF_on,
            fluid_overpressure_coupling=false,
        )
        S_6_off = L_6_off \ R_6_off
        pr_6_off = zeros(Ny1, Nx1)
        pf_6_off = zeros(Ny1, Nx1)
        vx_6_off = zeros(Ny1, Nx1)
        vy_6_off = zeros(Ny1, Nx1)
        qx_6_off = zeros(Ny1, Nx1)
        qy_6_off = zeros(Ny1, Nx1)
        process_hydromechanical_solution!(
            S_6_off,
            vx_6_off,
            vy_6_off,
            pr_6_off,
            qx_6_off,
            qy_6_off,
            pf_6_off;
            coords=coords,
        )

        # Overpressure pulse generation
        dp_fluid = pf_6_on[i_int, j_int] - pf_6_off[i_int, j_int]
        dp_total = pr_6_on[i_int, j_int] - pr_6_off[i_int, j_int]
        @test dp_fluid > 0.0
        @test dp_total > 0.0
        @test dp_fluid > dp_total
        # In low-permeability medium, fluid overpressure reaches order of MPa
        @test 1.0e6 < dp_fluid < 1.0e7

        # Invariant: Negative DQPF (hydration) produces symmetric negative pressure pulse
        DQPF_neg = zeros(Ny1, Nx1)
        DQPF_neg[i_int, j_int] = -dqpf_val
        R_6_neg = zeros(Nx1 * Ny1 * 6)
        L_6_neg = assemble_hydromechanical_lse!(
            ETA,
            ETAP,
            GGG,
            GGGP,
            SXY0,
            SXX0,
            RHOX,
            RHOY,
            RHOFX,
            RHOFY,
            RX,
            RY,
            ETAPHI,
            BETAPHI,
            PHI,
            gx,
            gy,
            pr0,
            pf0,
            DMP,
            dt_val,
            R_6_neg;
            coords=coords,
            DQPF=DQPF_neg,
            fluid_overpressure_coupling=true,
        )
        S_6_neg = L_6_neg \ R_6_neg
        pr_6_neg = zeros(Ny1, Nx1)
        pf_6_neg = zeros(Ny1, Nx1)
        vx_6_neg = zeros(Ny1, Nx1)
        vy_6_neg = zeros(Ny1, Nx1)
        qx_6_neg = zeros(Ny1, Nx1)
        qy_6_neg = zeros(Ny1, Nx1)
        process_hydromechanical_solution!(
            S_6_neg,
            vx_6_neg,
            vy_6_neg,
            pr_6_neg,
            qx_6_neg,
            qy_6_neg,
            pf_6_neg;
            coords=coords,
        )
        dp_fluid_neg = pf_6_neg[i_int, j_int] - pf_6_off[i_int, j_int]
        @test dp_fluid_neg < 0.0
        @test isapprox(dp_fluid_neg, -dp_fluid; rtol=1e-6)

        # Invariant: Turning coupling off suppresses overpressure response (matches background psurface)
        @test isapprox(maximum(pf_6_off), 1000.0; rtol=1e-5)
        @test isapprox(maximum(pr_6_off), 1000.0; rtol=1e-5)

        # 3. Solve 4-variable condensed system with fluid_overpressure_coupling = false
        R_4_off = zeros(Nx1 * Ny1 * 4)
        L_4_off = assemble_hydromechanical_4var_lse!(
            ETA,
            ETAP,
            GGG,
            GGGP,
            SXY0,
            SXX0,
            RHOX,
            RHOY,
            RHOFX,
            RHOFY,
            RX,
            RY,
            ETAPHI,
            BETAPHI,
            PHI,
            gx,
            gy,
            pr0,
            pf0,
            DMP,
            dt_val,
            R_4_off;
            coords=coords,
            DQPF=DQPF_on,
            fluid_overpressure_coupling=false,
        )
        S_4_off = L_4_off \ R_4_off
        pr_4_off = zeros(Ny1, Nx1)
        pf_4_off = zeros(Ny1, Nx1)
        vx_4_off = zeros(Ny1, Nx1)
        vy_4_off = zeros(Ny1, Nx1)
        process_hydromechanical_4var_solution!(
            S_4_off, vx_4_off, vy_4_off, pr_4_off, pf_4_off; coords=coords
        )
        @test isapprox(maximum(pf_4_off), 1000.0; rtol=1e-5)
        @test isapprox(maximum(pr_4_off), 1000.0; rtol=1e-5)
    end

    @testset "4-Variable Condensed vs 6-Variable Solution Equivalence Under DQPF" begin
        coords = GridCoordinates(GridConfig(Nx=7, Ny=7, xsize=35000.0, ysize=35000.0))
        Ny, Nx = coords.Ny, coords.Nx
        Ny1, Nx1 = coords.Ny1, coords.Nx1

        ETA = fill(1.0e21, Ny, Nx)
        ETAP = fill(1.0e21, Ny1, Nx1)
        GGG = fill(1.0e10, Ny, Nx)
        GGGP = fill(1.0e10, Ny1, Nx1)
        SXY0 = zeros(Ny, Nx)
        SXX0 = zeros(Ny, Nx)
        RHOX = fill(3000.0, Ny1, Nx1)
        RHOY = fill(3000.0, Ny1, Nx1)
        RHOFX = fill(1000.0, Ny1, Nx1)
        RHOFY = fill(1000.0, Ny1, Nx1)
        # Moderate permeability
        RX = fill(1.0e-3 / 1.0e-14, Ny1, Nx1)
        RY = fill(1.0e-3 / 1.0e-14, Ny1, Nx1)
        PHI = fill(0.1, Ny1, Nx1)
        ETAPHI = fill(1.0e24, Ny1, Nx1)
        BETAPHI = fill(1.0e-10, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = zeros(Ny1, Nx1)
        pr0 = zeros(Ny1, Nx1)
        pf0 = zeros(Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)
        dt_val = 1.0e7

        i_int, j_int = 4, 4
        dqpf_val = 1.5e-12
        DQPF = zeros(Ny1, Nx1)
        DQPF[i_int, j_int] = dqpf_val

        # Solve 6-variable system
        R6 = zeros(Nx1 * Ny1 * 6)
        L6 = assemble_hydromechanical_lse!(
            ETA,
            ETAP,
            GGG,
            GGGP,
            SXY0,
            SXX0,
            RHOX,
            RHOY,
            RHOFX,
            RHOFY,
            RX,
            RY,
            ETAPHI,
            BETAPHI,
            PHI,
            gx,
            gy,
            pr0,
            pf0,
            DMP,
            dt_val,
            R6;
            coords=coords,
            DQPF=DQPF,
            fluid_overpressure_coupling=true,
        )
        S6 = L6 \ R6
        pr_6 = zeros(Ny1, Nx1)
        pf_6 = zeros(Ny1, Nx1)
        vx_6 = zeros(Ny1, Nx1)
        vy_6 = zeros(Ny1, Nx1)
        qx_6 = zeros(Ny1, Nx1)
        qy_6 = zeros(Ny1, Nx1)
        process_hydromechanical_solution!(
            S6, vx_6, vy_6, pr_6, qx_6, qy_6, pf_6; coords=coords
        )

        # Solve 4-variable condensed system
        R4 = zeros(Nx1 * Ny1 * 4)
        L4 = assemble_hydromechanical_4var_lse!(
            ETA,
            ETAP,
            GGG,
            GGGP,
            SXY0,
            SXX0,
            RHOX,
            RHOY,
            RHOFX,
            RHOFY,
            RX,
            RY,
            ETAPHI,
            BETAPHI,
            PHI,
            gx,
            gy,
            pr0,
            pf0,
            DMP,
            dt_val,
            R4;
            coords=coords,
            DQPF=DQPF,
            fluid_overpressure_coupling=true,
        )
        S4 = L4 \ R4
        pr_4 = zeros(Ny1, Nx1)
        pf_4 = zeros(Ny1, Nx1)
        vx_4 = zeros(Ny1, Nx1)
        vy_4 = zeros(Ny1, Nx1)
        process_hydromechanical_4var_solution!(S4, vx_4, vy_4, pr_4, pf_4; coords=coords)

        # Fluid pressure equivalence
        max_pf_diff = maximum(abs, pf_4 .- pf_6)
        rel_pf_err = max_pf_diff / maximum(abs, pf_6)
        @test rel_pf_err < 1.0e-8
        @test pf_4[i_int, j_int] > 0.0

        # Total pressure equivalence
        max_pr_diff = maximum(abs, pr_4 .- pr_6)
        rel_pr_err = max_pr_diff / max(maximum(abs, pr_6), 1.0)
        @test rel_pr_err < 1.0e-8
    end

    @testset "Outward Darcy Filtration Driven by Dehydration Overpressure" begin
        coords = GridCoordinates(GridConfig(Nx=7, Ny=7, xsize=35000.0, ysize=35000.0))
        Ny, Nx = coords.Ny, coords.Nx
        Ny1, Nx1 = coords.Ny1, coords.Nx1

        ETA = fill(1.0e21, Ny, Nx)
        ETAP = fill(1.0e21, Ny1, Nx1)
        GGG = fill(1.0e10, Ny, Nx)
        GGGP = fill(1.0e10, Ny1, Nx1)
        SXY0 = zeros(Ny, Nx)
        SXX0 = zeros(Ny, Nx)
        RHOX = fill(3000.0, Ny1, Nx1)
        RHOY = fill(3000.0, Ny1, Nx1)
        RHOFX = fill(1000.0, Ny1, Nx1)
        RHOFY = fill(1000.0, Ny1, Nx1)
        RX = fill(1.0e-3 / 1.0e-12, Ny1, Nx1)
        RY = fill(1.0e-3 / 1.0e-12, Ny1, Nx1)
        PHI = fill(0.1, Ny1, Nx1)
        ETAPHI = fill(1.0e24, Ny1, Nx1)
        BETAPHI = fill(1.0e-10, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = zeros(Ny1, Nx1)
        pr0 = zeros(Ny1, Nx1)
        pf0 = zeros(Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)
        dt_val = 1.0e7

        i_int, j_int = 4, 4
        dqpf_val = 2.0e-12
        DQPF = zeros(Ny1, Nx1)
        DQPF[i_int, j_int] = dqpf_val

        R4 = zeros(Nx1 * Ny1 * 4)
        L4 = assemble_hydromechanical_4var_lse!(
            ETA,
            ETAP,
            GGG,
            GGGP,
            SXY0,
            SXX0,
            RHOX,
            RHOY,
            RHOFX,
            RHOFY,
            RX,
            RY,
            ETAPHI,
            BETAPHI,
            PHI,
            gx,
            gy,
            pr0,
            pf0,
            DMP,
            dt_val,
            R4;
            coords=coords,
            DQPF=DQPF,
            fluid_overpressure_coupling=true,
        )
        S4 = L4 \ R4
        pr_4 = zeros(Ny1, Nx1)
        pf_4 = zeros(Ny1, Nx1)
        vx_4 = zeros(Ny1, Nx1)
        vy_4 = zeros(Ny1, Nx1)
        process_hydromechanical_4var_solution!(S4, vx_4, vy_4, pr_4, pf_4; coords=coords)

        qxD = zeros(Ny1, Nx1)
        qyD = zeros(Ny1, Nx1)
        reconstruct_darcy_fluxes!(qxD, qyD, pf_4, RHOFX, RHOFY, RX, RY, gx, gy, coords)

        # Darcy flux directed outward from dehydration source cell (i_int=4, j_int=4)
        @test qxD[i_int, j_int] > 0.0
        @test qxD[i_int, j_int - 1] < 0.0
        @test qyD[i_int, j_int] > 0.0
        @test qyD[i_int - 1, j_int] < 0.0

        # Divergence of Darcy flux at the dehydration cell is positive
        div_qD =
            (qxD[i_int, j_int] - qxD[i_int, j_int - 1]) / coords.dx +
            (qyD[i_int, j_int] - qyD[i_int - 1, j_int]) / coords.dy
        @test div_qD > 0.9 * dqpf_val
        @test div_qD <= dqpf_val * 1.05
    end

    @testset "Unified Surface Venting Water Mass Balance & Speciation" begin
        delta_m_vent_3d = 5.0e14
        vented_vols = (
            M_vent_H2O=3.0e14,
            M_vent_C=2.0e13,
            M_vent_N=5.0e12,
            M_vent_S=1.0e13,
            M_vent_volatiles_total=3.35e14,
        )
        L_3D_equiv = 1.2
        dt_val = 1000.0
        P_surf = 1.0e5
        T_surf = 1200.0
        rplanet = 5.0e4
        marknum = 10
        tm = fill(1, marknum)
        xm = fill(0.9 * rplanet, marknum)
        ym = zeros(marknum)
        redox_props = (deltaIW_m=fill(-2.0, marknum),)

        # 1. Speciation inactive: additive water mass conservation (pore + mineral)
        cfg_nospec = SimulationConfig(
            volatiles=VolatilesConfig(speciation_active=false),
            venting=VentingConfig(active=true),
            retention=RetentionConfig(active=true, venting_drainage_active=true),
        )

        rates_add = compute_surface_venting_rates(
            cfg_nospec,
            delta_m_vent_3d,
            vented_vols,
            L_3D_equiv,
            dt_val,
            P_surf,
            T_surf,
            redox_props,
            marknum,
            tm,
            xm,
            ym,
            rplanet,
        )
        expected_H2O_step = delta_m_vent_3d + vented_vols.M_vent_H2O * L_3D_equiv
        @test rates_add[:H2O] ≈ expected_H2O_step / dt_val
        @test rates_add[:CO2] ≈
            (vented_vols.M_vent_C * L_3D_equiv * (44.0095 / 12.011)) / dt_val
        @test rates_add[:N2] ≈ (vented_vols.M_vent_N * L_3D_equiv) / dt_val
        @test rates_add[:H2S] ≈
            (vented_vols.M_vent_S * L_3D_equiv * (34.08 / 32.06)) / dt_val

        # 2. Pore water only (mineral drainage inactive)
        cfg_pore_only = SimulationConfig(
            volatiles=VolatilesConfig(speciation_active=false),
            venting=VentingConfig(active=true),
            retention=RetentionConfig(active=false, venting_drainage_active=false),
        )

        rates_pore = compute_surface_venting_rates(
            cfg_pore_only,
            delta_m_vent_3d,
            vented_vols,
            L_3D_equiv,
            dt_val,
            P_surf,
            T_surf,
            redox_props,
            marknum,
            tm,
            xm,
            ym,
            rplanet,
        )
        @test rates_pore[:H2O] ≈ delta_m_vent_3d / dt_val
        @test rates_pore[:CO2] ≈ 0.0
        @test rates_pore[:N2] ≈ 0.0
        @test rates_pore[:H2S] ≈ 0.0

        # 3. Mineral drainage only (Darcy pore venting inactive)
        cfg_drain_only = SimulationConfig(
            volatiles=VolatilesConfig(speciation_active=false),
            venting=VentingConfig(active=false),
            retention=RetentionConfig(active=true, venting_drainage_active=true),
        )

        rates_drain = compute_surface_venting_rates(
            cfg_drain_only,
            delta_m_vent_3d,
            vented_vols,
            L_3D_equiv,
            dt_val,
            P_surf,
            T_surf,
            redox_props,
            marknum,
            tm,
            xm,
            ym,
            rplanet,
        )
        @test rates_drain[:H2O] ≈ (vented_vols.M_vent_H2O * L_3D_equiv) / dt_val
        @test rates_drain[:CO2] ≈
            (vented_vols.M_vent_C * L_3D_equiv * (44.0095 / 12.011)) / dt_val

        # 3b. Non-speciation venting with non-water pore species (:CO2)
        cfg_co2_vent = SimulationConfig(
            volatiles=VolatilesConfig(speciation_active=false),
            venting=VentingConfig(active=true, species=:CO2),
            retention=RetentionConfig(active=true, venting_drainage_active=true),
        )
        rates_co2 = compute_surface_venting_rates(
            cfg_co2_vent,
            delta_m_vent_3d,
            vented_vols,
            L_3D_equiv,
            dt_val,
            P_surf,
            T_surf,
            redox_props,
            marknum,
            tm,
            xm,
            ym,
            rplanet,
        )
        # Mineral water routes to :H2O, pore fluid routes to :CO2 (plus mineral C)
        @test rates_co2[:H2O] ≈ (vented_vols.M_vent_H2O * L_3D_equiv) / dt_val
        @test rates_co2[:CO2] ≈
            (delta_m_vent_3d + vented_vols.M_vent_C * L_3D_equiv * (44.0095 / 12.011)) /
              dt_val

        # 4. Speciation active: thermodynamic speciation of full additive inventory
        cfg_spec = SimulationConfig(
            volatiles=VolatilesConfig(speciation_active=true, graphite_saturation=false),
            venting=VentingConfig(active=true),
            retention=RetentionConfig(active=true, venting_drainage_active=true),
        )

        rates_spec = compute_surface_venting_rates(
            cfg_spec,
            delta_m_vent_3d,
            vented_vols,
            L_3D_equiv,
            dt_val,
            P_surf,
            T_surf,
            redox_props,
            marknum,
            tm,
            xm,
            ym,
            rplanet,
        )

        mH2O = rates_spec[:H2O] * dt_val
        mH2 = rates_spec[:H2] * dt_val
        mCO = rates_spec[:CO] * dt_val
        mCO2 = rates_spec[:CO2] * dt_val
        mCH4 = rates_spec[:CH4] * dt_val
        mN2 = rates_spec[:N2] * dt_val
        mNH3 = rates_spec[:NH3] * dt_val
        mH2S = rates_spec[:H2S] * dt_val
        mS2 = rates_spec[:S2] * dt_val
        mSO2 = rates_spec[:SO2] * dt_val

        nH_in = 2.0 * expected_H2O_step / 18.01528e-3
        nC_in = (vented_vols.M_vent_C * L_3D_equiv) / 12.011e-3
        nN_in = (vented_vols.M_vent_N * L_3D_equiv) / 14.007e-3
        nS_in = (vented_vols.M_vent_S * L_3D_equiv) / 32.06e-3

        nH_out =
            2.0 * (mH2 / 2.01588e-3) +
            2.0 * (mH2O / 18.01528e-3) +
            4.0 * (mCH4 / 16.0425e-3) +
            3.0 * (mNH3 / 17.0305e-3) +
            2.0 * (mH2S / 34.0809e-3)
        nC_out = mCO / 28.0101e-3 + mCO2 / 44.0095e-3 + mCH4 / 16.0425e-3
        nN_out = 2.0 * (mN2 / 28.0134e-3) + mNH3 / 17.0305e-3
        nS_out = mH2S / 34.0809e-3 + 2.0 * (mS2 / 64.12e-3) + mSO2 / 64.066e-3

        @test isapprox(nH_out, nH_in; rtol=1.0e-10)
        @test isapprox(nC_out, nC_in; rtol=1.0e-10)
        @test isapprox(nN_out, nN_in; rtol=1.0e-10)
        @test isapprox(nS_out, nS_in; rtol=1.0e-10)

        # 5. Early return guards
        rates_dt0 = compute_surface_venting_rates(
            cfg_spec,
            delta_m_vent_3d,
            vented_vols,
            L_3D_equiv,
            0.0,
            P_surf,
            T_surf,
            redox_props,
            marknum,
            tm,
            xm,
            ym,
            rplanet,
        )
        @test all(iszero, values(rates_dt0))

        rates_no_vent = compute_surface_venting_rates(
            SimulationConfig(
                venting=VentingConfig(active=false), retention=RetentionConfig(active=false)
            ),
            delta_m_vent_3d,
            nothing,
            L_3D_equiv,
            dt_val,
            P_surf,
            T_surf,
            redox_props,
            marknum,
            tm,
            xm,
            ym,
            rplanet,
        )
        @test all(iszero, values(rates_no_vent))

        rates_zero_p0 = compute_surface_venting_rates(
            cfg_spec,
            0.0,
            nothing,
            L_3D_equiv,
            dt_val,
            0.0,
            0.0,
            redox_props,
            marknum,
            tm,
            xm,
            ym,
            rplanet,
        )
        @test all(iszero, values(rates_zero_p0))

        # 6. Escape and venting configuration validations
        cfg_bad_esc = SimulationConfig(
            escape=EscapeConfig(active=true, multi_species=true, species_list=[:H2O, :CO2]),
            volatiles=VolatilesConfig(speciation_active=true),
        )
        @test_throws ArgumentError validate_config(cfg_bad_esc)

        cfg_bad_esc_single = SimulationConfig(
            escape=EscapeConfig(active=true, multi_species=false, species=:H2O),
            volatiles=VolatilesConfig(speciation_active=true),
        )
        @test_throws ArgumentError validate_config(cfg_bad_esc_single)

        cfg_bad_vent_spec = SimulationConfig(
            venting=VentingConfig(active=true, species=:CO2),
            volatiles=VolatilesConfig(speciation_active=true),
        )
        @test_throws ArgumentError validate_config(cfg_bad_vent_spec)

        cfg_bad_drain_esc = SimulationConfig(
            escape=EscapeConfig(active=true, multi_species=true, species_list=[:H2O, :CO2]),
            retention=RetentionConfig(active=true, venting_drainage_active=true),
            volatiles=VolatilesConfig(active=true, speciation_active=false),
        )
        @test_throws ArgumentError validate_config(cfg_bad_drain_esc)

        # 7. Surface-mean delta_IW with center coordinates
        xc_test = 70_000.0
        yc_test = 70_000.0
        rp_test = 50_000.0
        xm_c = [70_000.0, 115_000.0]
        ym_c = [70_000.0, 70_000.0]
        tm_c = Int32[1, 1]
        redox_c = (deltaIW_m=[-4.0, +1.5],)
        diw_surf = compute_surface_mean_delta_iw(
            redox_c, tm_c, xm_c, ym_c, 2, xc_test, yc_test, rp_test, 0.0
        )
        @test diw_surf ≈ 1.5
    end

    @testset "Dynamic Equilibrium Gas Speciation Under Reducing Surface Venting" begin
        m_H2O = 1.0e12
        m_C = 5.0e11
        m_N = 5.0e10
        m_S = 1.0e11
        p_surf = 1.0e5
        T_surf = 1200.0

        # 1. Without graphite saturation: 100% of all elements stay in the gas phase
        gas_spec = speciate_vented_volatiles(
            m_H2O, m_C, m_N, m_S, p_surf, T_surf, -3.0; graphite_saturation=false
        )

        @test gas_spec[:H2] > 0.0
        @test gas_spec[:CO] > 0.0
        @test gas_spec[:H2] > gas_spec[:H2O]
        @test gas_spec[:CO] > gas_spec[:CO2]

        nH_in = 2.0 * m_H2O / 18.01528e-3
        nC_in = m_C / 12.011e-3
        nN_in = m_N / 14.007e-3
        nS_in = m_S / 32.06e-3

        nH_out =
            2.0 * (gas_spec[:H2] / 2.01588e-3) +
            2.0 * (gas_spec[:H2O] / 18.01528e-3) +
            4.0 * (gas_spec[:CH4] / 16.0425e-3) +
            3.0 * (gas_spec[:NH3] / 17.0305e-3) +
            2.0 * (gas_spec[:H2S] / 34.0809e-3)

        nC_out =
            gas_spec[:CO] / 28.0101e-3 +
            gas_spec[:CO2] / 44.0095e-3 +
            gas_spec[:CH4] / 16.0425e-3

        nN_out = 2.0 * (gas_spec[:N2] / 28.0134e-3) + gas_spec[:NH3] / 17.0305e-3

        nS_out =
            gas_spec[:H2S] / 34.0809e-3 +
            2.0 * (gas_spec[:S2] / 64.12e-3) +
            gas_spec[:SO2] / 64.066e-3

        @test isapprox(nH_out, nH_in; rtol=1.0e-10)
        @test isapprox(nC_out, nC_in; rtol=1.0e-10)
        @test isapprox(nN_out, nN_in; rtol=1.0e-10)
        @test isapprox(nS_out, nS_in; rtol=1.0e-10)

        # 2. With graphite saturation (default): solid carbon precipitates under reducing conditions
        sat_spec = speciate_vented_volatiles(
            m_H2O, m_C, m_N, m_S, p_surf, T_surf, -3.0; graphite_saturation=true
        )
        @test sat_spec[:H2] > 0.0
        @test sat_spec[:CO] > 0.0
        nC_sat =
            sat_spec[:CO] / 28.0101e-3 +
            sat_spec[:CO2] / 44.0095e-3 +
            sat_spec[:CH4] / 16.0425e-3
        @test nC_sat < nC_in
    end
end

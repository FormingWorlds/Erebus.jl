using Test
using Random
using Erebus
using Erebus.Config
using Erebus.Physics
using Erebus.Particles
using Erebus.Numerics
using StaticArrays
using TOML

Random.seed!(42)

@testset "Two-Phase Silicate Melt Segregation & Magma Ascent" begin
    @testset "MagmaTransportConfig Schema & Validation" begin
        cfg_def = MagmaTransportConfig()
        @test cfg_def.active == false
        @test isapprox(cfg_def.k_melt_ref, 1.0e-11; rtol=1e-12)
        @test isapprox(cfg_def.perm_exponent, 3.0; rtol=1e-12)
        @test isapprox(cfg_def.phi0, 0.10; rtol=1e-12)
        @test isapprox(cfg_def.phi_residual, 0.01; rtol=1e-12)
        @test isapprox(cfg_def.phi_crit, 0.40; rtol=1e-12)
        @test isapprox(cfg_def.phi_pack, 1.0; rtol=1e-12)
        @test isapprox(cfg_def.eta_melt, 10.0; rtol=1e-12)
        @test isapprox(cfg_def.r_grain, 1.0e-3; rtol=1e-12)
        @test isapprox(cfg_def.hindered_exponent, 2.0; rtol=1e-12)
        @test isapprox(cfg_def.F_perc_end, 0.35; rtol=1e-12)
        @test isapprox(cfg_def.F_settle_start, 0.45; rtol=1e-12)
        @test isapprox(cfg_def.cfl_melt, 0.5; rtol=1e-12)
        @test cfg_def.max_subcycles == 2000
        @test cfg_def.segregation_heating == true
        @test cfg_def.latent_crystallization == true
        @test cfg_def.exsolution_active == true
        @test cfg_def.track_depletion == true
        @test cfg_def.compaction_active == false
        @test isapprox(cfg_def.bulk_viscosity_ratio, 1.0; rtol=1e-12)
        @test isapprox(cfg_def.min_bulk_porosity, 0.005; rtol=1e-12)
        @test isapprox(cfg_def.compaction_length_min, 100.0; rtol=1e-12)
        @test isapprox(cfg_def.compaction_length_max, 50000.0; rtol=1e-12)
        @test cfg_def.ponding_active == false
        @test cfg_def.eruption_active == false
        @test isapprox(cfg_def.tensile_strength, 1.0e7; rtol=1e-12)

        # Integration in SimulationConfig
        sim_cfg = default_config()
        @test sim_cfg.magma_transport isa MagmaTransportConfig
        @test validate_config(sim_cfg) === nothing

        # Validation dependency: magma_transport requires active melting
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(; active=false),
                magma_transport=MagmaTransportConfig(; active=true),
            ),
        )

        m_active = MeltingConfig(; active=true)

        # Validation bounds: negative reference permeability
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, k_melt_ref=-1.0e-11),
            ),
        )

        # Validation bounds: non-positive permeability exponent
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, perm_exponent=0.0),
            ),
        )

        # Validation bounds: non-positive reference porosity
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, phi0=0.0),
            ),
        )

        # Validation bounds: unphysical residual porosity
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, phi_residual=-0.01),
            ),
        )

        # Validation bounds: unphysical critical melt fraction
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, phi_crit=0.0),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, phi_crit=1.0),
            ),
        )

        # Validation bounds: unphysical maximum packing fraction
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, phi_pack=1.05),
            ),
        )

        # Validation bounds: regime hierarchy ordering violation (phi_crit < F_perc_end)
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(;
                    active=true, phi_crit=0.30, F_perc_end=0.35, F_settle_start=0.45
                ),
            ),
        )

        # Validation bounds: regime hierarchy ordering violation (phi_pack < F_settle_start)
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(;
                    active=true, phi_pack=0.40, F_settle_start=0.45
                ),
            ),
        )

        # Validation bounds: regime hierarchy ordering violation (phi_residual >= F_perc_end)
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(;
                    active=true, phi_residual=0.36, F_perc_end=0.35
                ),
            ),
        )

        # Validation bounds: non-positive melt viscosity
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, eta_melt=0.0),
            ),
        )

        # Validation bounds: non-positive grain size
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, r_grain=-1.0e-3),
            ),
        )

        # Validation bounds: negative hindered settling exponent
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, hindered_exponent=-0.5),
            ),
        )

        # Validation bounds: unphysical CFL parameter
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, cfl_melt=0.0),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, cfl_melt=1.5),
            ),
        )

        # Validation bounds: invalid maximum subcycles
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, max_subcycles=0),
            ),
        )

        # Validation bounds: invalid bulk viscosity ratio
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(;
                    active=true, bulk_viscosity_ratio=-1.0
                ),
            ),
        )

        # Validation bounds: invalid min bulk porosity
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, min_bulk_porosity=0.0),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, min_bulk_porosity=0.45),
            ),
        )

        # Validation bounds: invalid compaction length bounds
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(;
                    active=true, compaction_length_min=50000.0, compaction_length_max=100.0
                ),
            ),
        )

        # Validation bounds: invalid tensile strength
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(;
                    active=true, tensile_strength=-1.0e6
                ),
            ),
        )
    end

    @testset "MagmaTransportConfig TOML Round-Trip" begin
        custom_magma = MagmaTransportConfig(;
            active=true,
            k_melt_ref=2.5e-11,
            perm_exponent=2.8,
            phi0=0.04,
            phi_residual=0.008,
            phi_crit=0.38,
            phi_pack=0.65,
            eta_melt=8.0,
            r_grain=1.5e-3,
            hindered_exponent=2.3,
            F_perc_end=0.32,
            F_settle_start=0.42,
            cfl_melt=0.4,
            max_subcycles=80,
            segregation_heating=true,
            latent_crystallization=true,
            exsolution_active=true,
            track_depletion=true,
            compaction_active=true,
            bulk_viscosity_ratio=1.5,
            min_bulk_porosity=0.008,
            compaction_length_min=200.0,
            compaction_length_max=40000.0,
            ponding_active=true,
            eruption_active=true,
            tensile_strength=1.5e7,
        )
        sim_cfg = SimulationConfig(;
            melting=MeltingConfig(; active=true), magma_transport=custom_magma
        )
        toml_str = Erebus.Config.save_config(sim_cfg)
        loaded_cfg = Erebus.Config.load_config(toml_str)

        @test loaded_cfg.magma_transport.active == true
        @test isapprox(loaded_cfg.magma_transport.k_melt_ref, 2.5e-11; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.perm_exponent, 2.8; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.phi0, 0.04; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.phi_residual, 0.008; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.phi_crit, 0.38; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.phi_pack, 0.65; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.eta_melt, 8.0; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.r_grain, 1.5e-3; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.hindered_exponent, 2.3; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.F_perc_end, 0.32; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.F_settle_start, 0.42; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.cfl_melt, 0.4; rtol=1e-12)
        @test loaded_cfg.magma_transport.max_subcycles == 80
        @test loaded_cfg.magma_transport.segregation_heating == true
        @test loaded_cfg.magma_transport.latent_crystallization == true
        @test loaded_cfg.magma_transport.exsolution_active == true
        @test loaded_cfg.magma_transport.track_depletion == true
        @test loaded_cfg.magma_transport.compaction_active == true
        @test isapprox(loaded_cfg.magma_transport.bulk_viscosity_ratio, 1.5; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.min_bulk_porosity, 0.008; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.compaction_length_min, 200.0; rtol=1e-12)
        @test isapprox(
            loaded_cfg.magma_transport.compaction_length_max, 40000.0; rtol=1e-12
        )
        @test loaded_cfg.magma_transport.ponding_active == true
        @test loaded_cfg.magma_transport.eruption_active == true
        @test isapprox(loaded_cfg.magma_transport.tensile_strength, 1.5e7; rtol=1e-12)
    end

    @testset "compaction_viscosity: McKenzie (1984) Bulk Viscosity" begin
        eta_s = 1.0e19
        # Above regularizer threshold: zeta = bulk_ratio * eta_s / F_m
        zeta1 = compaction_viscosity(eta_s, 0.10; bulk_ratio=1.0, phi_min=0.005)
        @test isapprox(zeta1, 1.0e20; rtol=1e-12)

        # Monotonicity with melt fraction
        zeta2 = compaction_viscosity(eta_s, 0.20; bulk_ratio=1.0, phi_min=0.005)
        @test zeta2 < zeta1
        @test isapprox(zeta2, 5.0e19; rtol=1e-12)

        # Scaling with bulk ratio
        zeta_scaled = compaction_viscosity(eta_s, 0.10; bulk_ratio=2.0, phi_min=0.005)
        @test isapprox(zeta_scaled, 2.0 * zeta1; rtol=1e-12)

        # Below regularizer threshold: clamped to phi_min
        zeta_floor = compaction_viscosity(eta_s, 0.001; bulk_ratio=1.0, phi_min=0.005)
        @test isapprox(zeta_floor, 1.0e19 / 0.005; rtol=1e-12)

        # Domain errors
        @test_throws DomainError compaction_viscosity(-1.0e19, 0.10)
        @test_throws DomainError compaction_viscosity(eta_s, -0.05)
        @test_throws DomainError compaction_viscosity(eta_s, 1.5)
        @test_throws DomainError compaction_viscosity(eta_s, 0.10; bulk_ratio=-1.0)
        @test_throws DomainError compaction_viscosity(eta_s, 0.10; phi_min=0.0)
    end

    @testset "compaction_length: Matrix Compaction Length Scale" begin
        eta_s = 1.0e19
        eta_m = 10.0
        k_m = 1.0e-11
        F_m = 0.10

        delta_c = compaction_length(eta_s, eta_m, k_m, F_m; bulk_ratio=1.0, phi_min=0.005)
        zeta = 1.0e19 / 0.10
        expected_delta = sqrt((zeta + (4.0 / 3.0) * eta_s) * k_m / eta_m)
        @test isapprox(delta_c, expected_delta; rtol=1e-12)
        @test 100.0 <= delta_c <= 50000.0

        # Clamping bounds
        delta_clamped_min = compaction_length(
            eta_s, eta_m, k_m, F_m; delta_min=20000.0, delta_max=50000.0
        )
        @test isapprox(delta_clamped_min, 20000.0; rtol=1e-12)

        delta_clamped_max = compaction_length(
            eta_s, eta_m, k_m, F_m; delta_min=100.0, delta_max=500.0
        )
        @test isapprox(delta_clamped_max, 500.0; rtol=1e-12)

        # Domain errors
        @test_throws DomainError compaction_length(-1.0, eta_m, k_m, F_m)
        @test_throws DomainError compaction_length(eta_s, -1.0, k_m, F_m)
        @test_throws DomainError compaction_length(eta_s, eta_m, -1.0, F_m)
        @test_throws DomainError compaction_length(
            eta_s, eta_m, k_m, F_m; delta_min=1000.0, delta_max=100.0
        )
    end

    @testset "compaction_pressure: Dynamic Compaction Stress" begin
        eta_s = 1.0e19
        F_m = 0.10

        # Compacting matrix: div_v < 0 => P_comp > 0
        div_v_comp = -1.0e-14
        P_comp = compaction_pressure(div_v_comp, eta_s, F_m)
        @test P_comp > 0.0
        @test isapprox(P_comp, 1.0e6; rtol=1e-12)

        # Dilating matrix: div_v > 0 => P_comp < 0
        div_v_dil = 1.0e-14
        P_dil = compaction_pressure(div_v_dil, eta_s, F_m)
        @test P_dil < 0.0
        @test isapprox(P_dil, -1.0e6; rtol=1e-12)

        # Domain errors
        @test_throws DomainError compaction_pressure(NaN, eta_s, F_m)
        @test_throws DomainError compaction_pressure(div_v_comp, -eta_s, F_m)
    end

    @testset "silicate_melt_permeability: McKenzie (1984) Power Law" begin
        k0 = 1.0e-11
        phi0 = 0.05
        n = 3.0
        phi_res = 0.005

        # Sub-residual melt fraction yields zero permeability
        @test iszero(
            silicate_melt_permeability(0.0; k0=k0, phi0=phi0, n=n, phi_residual=phi_res)
        )
        @test iszero(
            silicate_melt_permeability(0.003; k0=k0, phi0=phi0, n=n, phi_residual=phi_res)
        )
        @test iszero(
            silicate_melt_permeability(phi_res; k0=k0, phi0=phi0, n=n, phi_residual=phi_res)
        )

        # Reference porosity phi0 yields reference permeability k0 identically
        k_ref = silicate_melt_permeability(
            phi0; k0=k0, phi0=phi0, n=n, phi_residual=phi_res
        )
        @test isapprox(k_ref, k0; rtol=1e-12)

        # Power law scaling verification: double effective porosity yields 2^n times permeability
        k_double = silicate_melt_permeability(
            phi_res + 2.0 * (phi0 - phi_res); k0=k0, phi0=phi0, n=n, phi_residual=phi_res
        )
        @test isapprox(k_double, k0 * (2.0^n); rtol=1e-12)

        # Monotonicity test
        phi_vals = range(phi_res + 0.001, 0.35; length=30)
        k_vals = [
            silicate_melt_permeability(p; k0=k0, phi0=phi0, n=n, phi_residual=phi_res) for
            p in phi_vals
        ]
        @test issorted(k_vals)
        @test all(k_vals .> 0.0)

        # Domain error guards
        @test_throws DomainError silicate_melt_permeability(-0.01)
        @test_throws DomainError silicate_melt_permeability(NaN)
        @test_throws DomainError silicate_melt_permeability(0.05; k0=-1.0e-11)
        @test_throws DomainError silicate_melt_permeability(0.05; phi0=0.0)
        @test_throws DomainError silicate_melt_permeability(
            0.05; phi0=0.005, phi_residual=0.005
        )
        @test_throws DomainError silicate_melt_permeability(0.05; n=0.0)
        @test_throws DomainError silicate_melt_permeability(0.05; phi_residual=-0.01)
    end

    @testset "silicate_melt_segregation_velocity: Regime Blending & Limits" begin
        k0 = 1.0e-11
        phi0 = 0.05
        n_perm = 3.0
        phi_res = 0.005
        eta_liq = 10.0
        r_gr = 1.0e-3
        n_hind = 2.0
        drho = 300.0
        g = 0.2
        F_p_end = 0.35
        F_s_start = 0.45

        # Sub-residual melt fraction yields zero segregation velocity
        v_sub = silicate_melt_segregation_velocity(
            0.002,
            drho,
            g,
            eta_liq;
            k_melt_ref=k0,
            phi0=phi0,
            perm_exponent=n_perm,
            phi_residual=phi_res,
            r_grain=r_gr,
            hindered_exponent=n_hind,
            F_perc_end=F_p_end,
            F_settle_start=F_s_start,
        )
        @test iszero(v_sub)

        # Zero gravity yields zero velocity
        v_g0 = silicate_melt_segregation_velocity(
            0.1,
            drho,
            0.0,
            eta_liq;
            k_melt_ref=k0,
            phi0=phi0,
            perm_exponent=n_perm,
            phi_residual=phi_res,
            r_grain=r_gr,
            hindered_exponent=n_hind,
            F_perc_end=F_p_end,
            F_settle_start=F_s_start,
        )
        @test iszero(v_g0)

        # Negative density contrast yields zero buoyant velocity
        v_neg = silicate_melt_segregation_velocity(
            0.1,
            -10.0,
            g,
            eta_liq;
            k_melt_ref=k0,
            phi0=phi0,
            perm_exponent=n_perm,
            phi_residual=phi_res,
            r_grain=r_gr,
            hindered_exponent=n_hind,
            F_perc_end=F_p_end,
            F_settle_start=F_s_start,
        )
        @test iszero(v_neg)

        # Low melt fraction (Darcy percolation regime: F_m = 0.1 <= F_p_end)
        v_darcy = silicate_melt_segregation_velocity(
            0.1,
            drho,
            g,
            eta_liq;
            k_melt_ref=k0,
            phi0=phi0,
            perm_exponent=n_perm,
            phi_residual=phi_res,
            r_grain=r_gr,
            hindered_exponent=n_hind,
            F_perc_end=F_p_end,
            F_settle_start=F_s_start,
        )
        k_expected = silicate_melt_permeability(
            0.1; k0=k0, phi0=phi0, n=n_perm, phi_residual=phi_res
        )
        v_darcy_expected = (k_expected / (eta_liq * 0.1)) * drho * g
        @test isapprox(v_darcy, v_darcy_expected; rtol=1e-12)

        # High melt fraction (Stokes crystal suspension regime: F_m = 0.6 >= F_s_start)
        # Richardson-Zaki hindered settling: v_susp = v_stokes * F_m^n
        v_stokes = silicate_melt_segregation_velocity(
            0.6,
            drho,
            g,
            eta_liq;
            k_melt_ref=k0,
            phi0=phi0,
            perm_exponent=n_perm,
            phi_residual=phi_res,
            r_grain=r_gr,
            hindered_exponent=n_hind,
            F_perc_end=F_p_end,
            F_settle_start=F_s_start,
        )
        v_stokes_base = (2.0 / 9.0) * (r_gr^2) * drho * g / eta_liq
        v_stokes_expected = v_stokes_base * (0.6^n_hind)
        @test isapprox(v_stokes, v_stokes_expected; rtol=1e-12)

        # Pure melt limit (F_m = 1.0): unhindered Stokes velocity
        v_pure = silicate_melt_segregation_velocity(
            1.0,
            drho,
            g,
            eta_liq;
            k_melt_ref=k0,
            phi0=phi0,
            perm_exponent=n_perm,
            phi_residual=phi_res,
            r_grain=r_gr,
            hindered_exponent=n_hind,
            F_perc_end=F_p_end,
            F_settle_start=F_s_start,
        )
        @test isapprox(v_pure, v_stokes_base; rtol=1e-12)

        # Transition regime (F_perc_end < F_m < F_settle_start): smooth Hermite interpolation
        F_mid = 0.5 * (F_p_end + F_s_start)
        v_trans = silicate_melt_segregation_velocity(
            F_mid,
            drho,
            g,
            eta_liq;
            k_melt_ref=k0,
            phi0=phi0,
            perm_exponent=n_perm,
            phi_residual=phi_res,
            r_grain=r_gr,
            hindered_exponent=n_hind,
            F_perc_end=F_p_end,
            F_settle_start=F_s_start,
        )
        v_p_end = silicate_melt_segregation_velocity(
            F_p_end,
            drho,
            g,
            eta_liq;
            k_melt_ref=k0,
            phi0=phi0,
            perm_exponent=n_perm,
            phi_residual=phi_res,
            r_grain=r_gr,
            hindered_exponent=n_hind,
            F_perc_end=F_p_end,
            F_settle_start=F_s_start,
        )
        v_s_start = silicate_melt_segregation_velocity(
            F_s_start,
            drho,
            g,
            eta_liq;
            k_melt_ref=k0,
            phi0=phi0,
            perm_exponent=n_perm,
            phi_residual=phi_res,
            r_grain=r_gr,
            hindered_exponent=n_hind,
            F_perc_end=F_p_end,
            F_settle_start=F_s_start,
        )
        @test v_trans >= min(v_p_end, v_s_start)
        @test v_trans <= max(v_p_end, v_s_start)

        # Domain error guards
        @test_throws DomainError silicate_melt_segregation_velocity(-0.1, drho, g, eta_liq)
        @test_throws DomainError silicate_melt_segregation_velocity(1.1, drho, g, eta_liq)
        @test_throws DomainError silicate_melt_segregation_velocity(0.2, NaN, g, eta_liq)
        @test_throws DomainError silicate_melt_segregation_velocity(
            0.2, drho, -0.1, eta_liq
        )
        @test_throws DomainError silicate_melt_segregation_velocity(0.2, drho, g, 0.0)
    end

    @testset "silicate_melt_dissipation_heating: Gravitational Energy Release" begin
        v_seg = 1.0e-5
        drho = 300.0
        g = 0.5
        F_m = 0.2
        psi = silicate_melt_dissipation_heating(F_m, drho, g, v_seg)
        @test isapprox(psi, drho * g * F_m * v_seg; rtol=1e-12)

        # Zero dissipation when stationary or neutral buoyancy
        @test iszero(silicate_melt_dissipation_heating(0.0, drho, g, v_seg))
        @test iszero(silicate_melt_dissipation_heating(F_m, 0.0, g, v_seg))
        @test iszero(silicate_melt_dissipation_heating(F_m, drho, 0.0, v_seg))
        @test iszero(silicate_melt_dissipation_heating(F_m, drho, g, 0.0))

        # Domain error guards
        @test_throws DomainError silicate_melt_dissipation_heating(-0.1, drho, g, v_seg)
        @test_throws DomainError silicate_melt_dissipation_heating(1.1, drho, g, v_seg)
        @test_throws DomainError silicate_melt_dissipation_heating(F_m, NaN, g, v_seg)
        @test_throws DomainError silicate_melt_dissipation_heating(F_m, drho, NaN, v_seg)
        @test_throws DomainError silicate_melt_dissipation_heating(F_m, drho, g, NaN)
    end

    @testset "setup_marker_magma_properties Allocation" begin
        marknum = 100
        magma_props = setup_marker_magma_properties(marknum)
        @test length(magma_props) == 1
        F_extract_m = magma_props[1]
        @test length(F_extract_m) == marknum
        @test all(iszero, F_extract_m)
        @test eltype(F_extract_m) === Float64
    end

    @testset "apply_silicate_melt_segregation! Conservation & Ascent" begin
        # 16x16 grid setup over 80 km x 80 km domain
        Nx = 16
        Ny = 16
        xsize = 80000.0
        ysize = 80000.0
        dx = xsize / (Nx - 1)
        dy = ysize / (Ny - 1)
        coords = GridCoordinates(Nx, Ny; xsize=xsize, ysize=ysize)
        xcenter = coords.xcenter
        ycenter = coords.ycenter
        rplanet = 30000.0
        g_surf = 0.2

        # Populate markers: 4 markers per cell in planetary interior
        marknum = 4 * Nx * Ny
        xm = zeros(Float64, marknum)
        ym = zeros(Float64, marknum)
        tm = fill(1, marknum)
        tkm = fill(1600.0, marknum)
        Fm = zeros(Float64, marknum)
        F_extract_m = zeros(Float64, marknum)

        idx = 1
        for i in 1:Ny, j in 1:Nx
            xc = (j - 1) * dx
            yc = (i - 1) * dy
            for sx in (-0.25, 0.25), sy in (-0.25, 0.25)
                xm[idx] = xc + sx * dx
                ym[idx] = yc + sy * dy
                rmark = distance(xm[idx], ym[idx], xcenter, ycenter)
                if rmark <= rplanet
                    tm[idx] = 1 # Rock mantle
                    # Partially molten plume between 10 km and 18 km radius
                    if 10000.0 <= rmark <= 18000.0
                        Fm[idx] = 0.20
                        tkm[idx] = 1600.0
                    elseif rmark < 6000.0
                        # Super-liquidus magma ocean core: Fm = 1.0 > phi_pack
                        Fm[idx] = 1.0
                        tkm[idx] = 2000.0
                    else
                        Fm[idx] = 0.0
                        # Exterior cold mantle below solidus (1400 K)
                        tkm[idx] = 1200.0
                    end
                else
                    tm[idx] = 3 # Sticky air
                    tkm[idx] = 200.0
                    Fm[idx] = 0.0
                end
                idx += 1
            end
        end

        # Test with phi_pack = 0.60 < 1.0: markers with Fm = 1.0 must not crash
        cfg_magma = MagmaTransportConfig(;
            active=true,
            k_melt_ref=1.0e-10,
            phi0=0.05,
            phi_residual=0.005,
            phi_crit=0.38,
            F_perc_end=0.32,
            F_settle_start=0.42,
            phi_pack=0.60,
            eta_melt=1.0,
            cfl_melt=0.4,
            max_subcycles=20,
            segregation_heating=true,
            latent_crystallization=true,
            track_depletion=true,
        )

        Q_seg_grid = zeros(Float64, Ny, Nx)
        Q_lat_grid = zeros(Float64, Ny, Nx)
        dt = 5.0e10 # Transport timestep (~1500 years)

        M_melt_initial = sum(Fm)
        @test M_melt_initial > 0.0

        # Run silicate melt segregation
        res = apply_silicate_melt_segregation!(
            xm,
            ym,
            tm,
            tkm,
            Fm,
            marknum,
            dt,
            cfg_magma;
            coords=coords,
            xcenter=xcenter,
            ycenter=ycenter,
            rplanet=rplanet,
            g_surf=g_surf,
            Q_seg_grid=Q_seg_grid,
            Q_lat_grid=Q_lat_grid,
            rho_silicate=3300.0,
            rho_melt=2800.0,
            T_solidus_silicate=1400.0,
            T_liquidus_silicate=1800.0,
            L_melt=4.0e5,
            F_extract_m=F_extract_m,
        )

        # Machine-precision mass conservation: remaining melt + crystallized melt
        M_melt_final = sum(Fm)
        M_conserved = M_melt_final + res.total_crystallized_mass
        relative_mass_drift = abs(M_conserved - M_melt_initial) / M_melt_initial
        @test relative_mass_drift < 1.0e-12

        # Radial upward / outward melt segregation
        @test res.max_v_seg > 0.0
        @test res.n_subcycles >= 1

        # Energy dissipation heating occurred and populated Q_seg_grid
        @test res.total_dissipation_energy > 0.0
        @test maximum(Q_seg_grid) > 0.0

        # Latent heat release matches crystallized mass to machine precision
        @test res.total_crystallized_mass > 0.0
        @test maximum(Q_lat_grid) > 0.0
        dV_cell = dx * dy
        E_lat_grid = sum(Q_lat_grid) * dt * dV_cell
        # Total crystallized mass on markers: each unit of Fm represents dV_cell / 4 mass equivalent
        E_lat_expected = res.total_crystallized_mass * (dV_cell / 4) * 2800.0 * 4.0e5
        @test maximum(Q_seg_grid) > 0.0

        # Depletion tracking accumulated on donor markers
        @test maximum(F_extract_m) > 0.0

        # Inactive simulation no-op
        cfg_inactive = MagmaTransportConfig(; active=false)
        res_inactive = apply_silicate_melt_segregation!(
            xm, ym, tm, tkm, Fm, marknum, dt, cfg_inactive; coords=coords
        )
        @test iszero(res_inactive.max_v_seg)
        @test res_inactive.n_subcycles == 0

        # Neutral and negative buoyancy tests (Finding 2)
        res_neutral = apply_silicate_melt_segregation!(
            xm,
            ym,
            tm,
            tkm,
            Fm,
            marknum,
            dt,
            cfg_magma;
            coords=coords,
            xcenter=xcenter,
            ycenter=ycenter,
            rplanet=rplanet,
            g_surf=g_surf,
            rho_silicate=2800.0,
            rho_melt=2800.0,
        )
        @test iszero(res_neutral.max_v_seg)
        @test res_neutral.n_subcycles == 0

        res_negative = apply_silicate_melt_segregation!(
            xm,
            ym,
            tm,
            tkm,
            Fm,
            marknum,
            dt,
            cfg_magma;
            coords=coords,
            xcenter=xcenter,
            ycenter=ycenter,
            rplanet=rplanet,
            g_surf=g_surf,
            rho_silicate=2700.0,
            rho_melt=2800.0,
        )
        @test iszero(res_negative.max_v_seg)
        @test res_negative.n_subcycles == 0
    end

    @testset "McKenzie (1984) Dynamic Compaction Pressure Coupling" begin
        Nx = 16
        Ny = 16
        xsize = 40000.0
        ysize = 40000.0
        coords = GridCoordinates(Nx, Ny; xsize=xsize, ysize=ysize)
        dx = coords.dx
        dy = coords.dy

        marknum = Nx * Ny * 4
        xm = zeros(Float64, marknum)
        ym = zeros(Float64, marknum)
        tm = ones(Int, marknum)
        tkm = fill(1600.0, marknum)
        Fm = fill(0.10, marknum)

        idx = 1
        for j in 1:Nx, i in 1:Ny
            for kx in (0.25, 0.75), ky in (0.25, 0.75)
                xm[idx] = (j - 1 + kx) * dx
                ym[idx] = (i - 1 + ky) * dy
                idx += 1
            end
        end

        div_v_mat = fill(-1.0e-14, Ny, Nx)
        ETA_mat = fill(1.0e19, Ny, Nx)

        cfg_comp = MagmaTransportConfig(;
            active=true,
            compaction_active=true,
            bulk_viscosity_ratio=1.0,
            min_bulk_porosity=0.005,
            compaction_length_min=100.0,
            compaction_length_max=50000.0,
        )

        dt = 1.0e9
        res = apply_silicate_melt_segregation!(
            xm,
            ym,
            tm,
            tkm,
            Fm,
            marknum,
            dt,
            cfg_comp;
            coords=coords,
            xcenter=20000.0,
            ycenter=20000.0,
            rplanet=50000.0,
            g_surf=0.1,
            div_v=div_v_mat,
            ETA=ETA_mat,
            rho_silicate=3000.0,
            rho_melt=2800.0,
        )

        @test res.max_compaction_pressure > 0.0
        @test isapprox(res.max_compaction_pressure, 1.0e6; rtol=1e-3)
        @test res.mean_compaction_length > 0.0
        @test 100.0 <= res.mean_compaction_length <= 50000.0
        @test res.n_subcycles >= 1

        # Non-zero compaction pressure gradient coupling verification
        div_v_grad = zeros(Float64, Ny, Nx)
        for i in 1:Ny
            div_v_grad[i, :] .= -1.0e-14 * (1.0 + 2.0 * (i - 1) / (Ny - 1))
        end
        Fm_on = copy(Fm)
        res_on = apply_silicate_melt_segregation!(
            copy(xm),
            copy(ym),
            copy(tm),
            copy(tkm),
            Fm_on,
            marknum,
            dt,
            cfg_comp;
            coords=coords,
            xcenter=20000.0,
            ycenter=20000.0,
            rplanet=50000.0,
            g_surf=0.1,
            div_v=div_v_grad,
            ETA=ETA_mat,
            rho_silicate=3000.0,
            rho_melt=2800.0,
        )

        cfg_no_comp = MagmaTransportConfig(; active=true, compaction_active=false)
        Fm_off = copy(Fm)
        res_off = apply_silicate_melt_segregation!(
            copy(xm),
            copy(ym),
            copy(tm),
            copy(tkm),
            Fm_off,
            marknum,
            dt,
            cfg_no_comp;
            coords=coords,
            xcenter=20000.0,
            ycenter=20000.0,
            rplanet=50000.0,
            g_surf=0.1,
            ETA=ETA_mat,
            rho_silicate=3000.0,
            rho_melt=2800.0,
        )

        # Dynamic compaction pressure gradient actively modifies melt distribution
        diff_fm = maximum(abs.(Fm_on .- Fm_off))
        @test diff_fm > 0.0

        # Production path: evaluate divergence from staggered vx/vy velocity field
        Nx1 = Nx + 1
        Ny1 = Ny + 1
        vx_prod = zeros(Float64, Ny1, Nx1)
        vy_prod = zeros(Float64, Ny1, Nx1)
        C_rate = -2.0e-14
        for j in 1:Nx1, i in 1:Ny1
            vx_prod[i, j] = C_rate * (j - 1) * dx
        end
        Fm_prod = copy(Fm)
        res_prod = apply_silicate_melt_segregation!(
            copy(xm),
            copy(ym),
            copy(tm),
            copy(tkm),
            Fm_prod,
            marknum,
            dt,
            cfg_comp;
            coords=coords,
            xcenter=20000.0,
            ycenter=20000.0,
            rplanet=50000.0,
            g_surf=0.1,
            vx=vx_prod,
            vy=vy_prod,
            ETA=ETA_mat,
            rho_silicate=3000.0,
            rho_melt=2800.0,
        )
        @test res_prod.max_compaction_pressure > 0.0
        @test isapprox(res_prod.max_compaction_pressure, 2.0e6; rtol=1e-2)
    end

    @testset "Scott & Stevenson (1984) 1D Compaction Column Benchmark" begin
        Nx = 3
        Ny = 40
        xsize = 6000.0
        ysize = 20000.0
        dx = xsize / Nx
        dy = ysize / Ny
        coords = GridCoordinates(Nx, Ny; xsize=xsize, ysize=ysize)

        n_per_cell = 4
        marknum = Nx * Ny * n_per_cell
        xm = zeros(Float64, marknum)
        ym = zeros(Float64, marknum)
        tm = ones(Int, marknum)
        tkm = fill(1600.0, marknum)
        phi_0 = 0.08
        Fm = fill(phi_0, marknum)

        idx = 1
        for j in 1:Nx, i in 1:Ny
            for kx in (0.25, 0.75), ky in (0.25, 0.75)
                xm[idx] = (j - 1 + kx) * dx
                ym[idx] = (i - 1 + ky) * dy
                idx += 1
            end
        end

        eta_s = 1.0e19
        eta_m = 10.0
        k0 = 1.0e-11
        drho = 200.0
        g = 0.5

        delta_c = compaction_length(eta_s, eta_m, k0, phi_0; bulk_ratio=1.0, phi_min=0.005)

        cfg_comp = MagmaTransportConfig(;
            active=true,
            compaction_active=true,
            bulk_viscosity_ratio=1.0,
            min_bulk_porosity=0.005,
            k_melt_ref=k0,
            eta_melt=eta_m,
            phi0=phi_0,
            perm_exponent=3.0,
            phi_residual=0.005,
        )

        gy_mat = fill(-g, Ny, Nx)
        gx_mat = zeros(Float64, Ny, Nx)
        ETA_mat = fill(eta_s, Ny, Nx)

        mass_init = sum(Fm)
        dt = 5.0e7
        res = apply_silicate_melt_segregation!(
            xm,
            ym,
            tm,
            tkm,
            Fm,
            marknum,
            dt,
            cfg_comp;
            coords=coords,
            xcenter=3000.0,
            ycenter=0.0,
            rplanet=25000.0,
            g_surf=g,
            gx=gx_mat,
            gy=gy_mat,
            ETA=ETA_mat,
            rho_silicate=3000.0,
            rho_melt=3000.0 - drho,
        )

        @test res.n_subcycles >= 1
        @test res.max_v_seg > 0.0
        @test res.mean_compaction_length > 0.0
        @test isapprox(res.mean_compaction_length, delta_c; rtol=0.05)
        @test 10000.0 <= res.mean_compaction_length <= 13000.0
        @test isapprox(sum(Fm), mass_init; rtol=1e-12)

        # 1D compaction column melt redistribution: melt ascends into upper column
        M_upper = sum(Fm[m] for m in 1:marknum if ym[m] >= ysize / 2.0)
        M_lower = sum(Fm[m] for m in 1:marknum if ym[m] < ysize / 2.0)
        @test M_upper > M_lower
    end

    @testset "Solitary Porosity Wave Speed Benchmark" begin
        eta_s = 1.0e19
        eta_m = 10.0
        k0 = 1.0e-11
        drho = 250.0
        g = 1.0
        n_perm = 3.0

        phi_low = 0.05
        phi_high = 0.10

        v_low = silicate_melt_segregation_velocity(
            phi_low,
            drho,
            g,
            eta_m;
            k_melt_ref=k0,
            phi0=0.10,
            perm_exponent=n_perm,
            phi_residual=0.01,
        )
        v_high = silicate_melt_segregation_velocity(
            phi_high,
            drho,
            g,
            eta_m;
            k_melt_ref=k0,
            phi0=0.10,
            perm_exponent=n_perm,
            phi_residual=0.01,
        )

        phi_eff_ratio = (phi_high - 0.01) / (phi_low - 0.01)
        expected_ratio = (phi_eff_ratio^n_perm) * (phi_low / phi_high)
        actual_ratio = v_high / v_low
        @test isapprox(actual_ratio, expected_ratio; rtol=1e-12)
        @test v_high > v_low
    end

    @testset "Decompression Volatile Exsolution on Ascending Melt" begin
        Nx = 8
        Ny = 8
        xsize = 20000.0
        ysize = 20000.0
        dx = xsize / Nx
        dy = ysize / Ny
        coords = GridCoordinates(Nx, Ny; xsize=xsize, ysize=ysize)

        marknum = Nx * Ny * 4
        xm = zeros(Float64, marknum)
        ym = zeros(Float64, marknum)
        tm = ones(Int, marknum)
        tkm = fill(1500.0, marknum)
        Fm = fill(0.12, marknum)
        XH2Om = fill(2.0, marknum)
        phim = zeros(Float64, marknum)

        idx = 1
        for j in 1:Nx, i in 1:Ny
            for _ in 1:4
                xm[idx] = (j - 0.5) * dx + (rand() - 0.5) * 0.4 * dx
                ym[idx] = (i - 0.5) * dy + (rand() - 0.5) * 0.4 * dy
                idx += 1
            end
        end

        cfg_magma = MagmaTransportConfig(; active=true, exsolution_active=true)
        cfg_vol = VolatilesConfig(; active=true)

        pr_mat = zeros(Float64, Ny, Nx)
        for i in 1:Ny
            pr_mat[i, :] .= 5.0e6 + (i - 1) / (Ny - 1) * 9.5e7
        end

        res = apply_silicate_melt_segregation!(
            xm,
            ym,
            tm,
            tkm,
            Fm,
            marknum,
            1.0e8,
            cfg_magma;
            coords=coords,
            xcenter=10000.0,
            ycenter=10000.0,
            rplanet=18000.0,
            g_surf=0.2,
            pr=pr_mat,
            XH2Om=XH2Om,
            phim=phim,
            cfg_volatiles=cfg_vol,
        )

        @test res.total_exsolved_volatiles > 0.0
        @test maximum(phim) > 0.0

        # Exact mass conservation: volatile loss equals pore volume gain via rho_s / rho_f
        w_lost = (2.0 * marknum - sum(XH2Om)) / 100.0
        phi_gained = sum(phim)
        @test w_lost > 0.0
        @test isapprox(phi_gained, w_lost * (3000.0 / 1000.0); rtol=1e-10)

        # Standalone exsolve_magma_volatiles! verification
        XH2O_single = [2.0]
        phi_single = [0.0]
        dw_single = exsolve_magma_volatiles!(
            [10000.0],
            [10000.0],
            [1],
            [1500.0],
            [0.12],
            1,
            XH2O_single,
            nothing,
            nothing,
            nothing,
            phi_single,
            cfg_vol;
            pr=fill(5.0e6, 8, 8),
            dx_val=2500.0,
            dy_val=2500.0,
            Nx_val=8,
            Ny_val=8,
            rho_silicate=3000.0,
            g_surf=0.2,
        )
        @test dw_single > 0.0
        @test isapprox(phi_single[1], dw_single * (3000.0 / 1000.0); rtol=1e-10)
    end

    @testset "Crustal Magma Sill Ponding & Overpressure Hydrofracture Eruption" begin
        Nx = 8
        Ny = 8
        xsize = 20000.0
        ysize = 20000.0
        dx = xsize / Nx
        dy = ysize / Ny
        coords = GridCoordinates(Nx, Ny; xsize=xsize, ysize=ysize)

        marknum = Nx * Ny * 4
        xm = zeros(Float64, marknum)
        ym = zeros(Float64, marknum)
        tm = ones(Int, marknum)
        tkm = zeros(Float64, marknum)
        Fm = zeros(Float64, marknum)

        idx = 1
        for j in 1:Nx, i in 1:Ny
            for _ in 1:4
                xm[idx] = (j - 0.5) * dx + (rand() - 0.5) * 0.4 * dx
                ym[idx] = (i - 0.5) * dy + (rand() - 0.5) * 0.4 * dy
                if i == 1
                    tkm[idx] = 1200.0
                    Fm[idx] = 0.0
                else
                    tkm[idx] = 1600.0
                    Fm[idx] = 0.15
                end
                idx += 1
            end
        end

        cfg_ponding = MagmaTransportConfig(;
            active=true,
            ponding_active=true,
            eruption_active=false,
            latent_crystallization=false,
            tensile_strength=1.0e7,
        )

        Fm_pond = copy(Fm)
        res_pond = apply_silicate_melt_segregation!(
            xm,
            ym,
            tm,
            tkm,
            Fm_pond,
            marknum,
            1.0e8,
            cfg_ponding;
            coords=coords,
            xcenter=10000.0,
            ycenter=20000.0,
            rplanet=25000.0,
            g_surf=0.2,
            T_solidus_silicate=1400.0,
        )

        lid_melt_pond = 0.0
        for m in 1:marknum
            if ym[m] < dy
                lid_melt_pond += Fm_pond[m]
            end
        end
        @test iszero(lid_melt_pond)

        # Paired runs isolating the tensile strength threshold trigger:
        # Case A: P_comp < tensile_strength => overpressure insufficient, melt blocked
        div_v_erupt = fill(-5.0e-13, Ny, Nx)
        cfg_sub_threshold = MagmaTransportConfig(;
            active=true,
            compaction_active=true,
            ponding_active=true,
            eruption_active=true,
            latent_crystallization=false,
            tensile_strength=1.0e9,
        )
        Fm_blocked = copy(Fm)
        res_blocked = apply_silicate_melt_segregation!(
            copy(xm),
            copy(ym),
            copy(tm),
            copy(tkm),
            Fm_blocked,
            marknum,
            1.0e8,
            cfg_sub_threshold;
            coords=coords,
            xcenter=10000.0,
            ycenter=20000.0,
            rplanet=25000.0,
            g_surf=0.2,
            T_solidus_silicate=1400.0,
            div_v=div_v_erupt,
        )
        lid_melt_blocked = sum(Fm_blocked[m] for m in 1:marknum if ym[m] < dy)
        @test iszero(lid_melt_blocked)

        # Case B: P_comp > tensile_strength => tensile hydrofracture breaches lid
        cfg_erupt = MagmaTransportConfig(;
            active=true,
            compaction_active=true,
            ponding_active=true,
            eruption_active=true,
            latent_crystallization=false,
            tensile_strength=1.0e6,
        )
        Fm_erupt = copy(Fm)
        res_erupt = apply_silicate_melt_segregation!(
            copy(xm),
            copy(ym),
            copy(tm),
            copy(tkm),
            Fm_erupt,
            marknum,
            1.0e8,
            cfg_erupt;
            coords=coords,
            xcenter=10000.0,
            ycenter=20000.0,
            rplanet=25000.0,
            g_surf=0.2,
            T_solidus_silicate=1400.0,
            div_v=div_v_erupt,
        )

        lid_melt_erupt = sum(Fm_erupt[m] for m in 1:marknum if ym[m] < dy)
        @test res_erupt.max_compaction_pressure > cfg_erupt.tensile_strength
        @test lid_melt_erupt > 0.0

        # Validate config checks
        cfg_invalid = SimulationConfig(
            magma_transport=MagmaTransportConfig(
                active=true, eruption_active=true, compaction_active=false
            ),
            melting=MeltingConfig(active=true),
        )
        @test_throws ArgumentError validate_config(cfg_invalid)
    end
end

using Test
using Erebus
using Erebus.Config
using Erebus.Physics
using Erebus.Particles
using Erebus.Geometry
using Erebus: start_hrsolidm, start_hrfluidm
using StaticArrays
using TOML
using JLD2

@testset "Iron Core Formation Physics & Configuration" begin
    @testset "CoreFormationConfig Schema & Defaults" begin
        cfg_def = CoreFormationConfig()
        @test cfg_def.percolation_active == false
        @test cfg_def.settling_active == false
        @test isapprox(cfg_def.rho_metal, 7200.0; rtol=1e-12)
        @test isapprox(cfg_def.eta_metal, 1.0e-2; rtol=1e-12)
        @test isapprox(cfg_def.k_metal, 40.0; rtol=1e-12)
        @test isapprox(cfg_def.rhocp_metal, 4.0e6; rtol=1e-12)
        @test isapprox(cfg_def.Xfe_bulk, 0.20; rtol=1e-12)
        @test isapprox(cfg_def.phi_pack, 0.65; rtol=1e-12)
        @test isapprox(cfg_def.T_eutectic, 1213.0; rtol=1e-12)
        @test isapprox(cfg_def.dT_metal, 50.0; rtol=1e-12)
        @test isapprox(cfg_def.k_metal_ref, 1.0e-9; rtol=1e-12)
        @test isapprox(cfg_def.perm_exponent, 3.0; rtol=1e-12)
        @test isapprox(cfg_def.phi_crit_perc, 0.05; rtol=1e-12)
        @test isapprox(cfg_def.phi_residual, 0.02; rtol=1e-12)
        @test cfg_def.droplet_size_mode == :weber_mean
        @test isapprox(cfg_def.droplet_diameter_fixed, 5.0e-3; rtol=1e-12)
        @test isapprox(cfg_def.sigma_metal_silicate, 1.0; rtol=1e-12)
        @test isapprox(cfg_def.We_crit, 10.0; rtol=1e-12)
        @test isapprox(cfg_def.hindered_exponent, 4.5; rtol=1e-12)
        @test cfg_def.hadamard_rybczynski == false
        @test isapprox(cfg_def.F_settle_start, 0.40; rtol=1e-12)
        @test isapprox(cfg_def.F_perc_end, 0.50; rtol=1e-12)
        @test cfg_def.segregation_heating == true
        @test isapprox(cfg_def.cfl_settling, 0.5; rtol=1e-12)
        @test cfg_def.max_subcycles == 2000

        # Default configuration integration
        sim_cfg = default_config()
        @test sim_cfg.coreformation isa CoreFormationConfig
        @test validate_config(sim_cfg) === nothing
    end

    @testset "CoreFormationConfig Bounds Validation" begin
        # Metal density must exceed silicate density (rhosolidm[1] = 3300) and be finite
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true,
                    rho_metal=3000.0,
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true,
                    rho_metal=Inf,
                ),
            ),
        )

        # Non-positive viscosity or thermal properties
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true,
                    eta_metal=-1.0e-2,
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    settling_active=true,
                    k_metal=0.0,
                ),
                melting=MeltingConfig(; active=true),
            ),
        )

        # Percolation ordering: phi_residual <= phi_crit_perc < phi_pack
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true,
                    phi_residual=0.10,
                    phi_crit_perc=0.05,
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true,
                    phi_crit_perc=0.70,
                    phi_pack=0.65,
                ),
            ),
        )

        # Regime handover interval: F_settle_start <= F_perc_end
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    settling_active=true,
                    F_settle_start=0.60,
                    F_perc_end=0.40,
                ),
                melting=MeltingConfig(; active=true),
            ),
        )

        # Droplet mode symbol
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    settling_active=true,
                    droplet_size_mode=:invalid_mode,
                ),
                melting=MeltingConfig(; active=true),
            ),
        )

        # Eutectic temperature must be below silicate solidus for percolation
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true,
                    T_eutectic=1500.0,
                ),
            ),
        )

        # Settling requires melting active
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    settling_active=true,
                ),
                melting=MeltingConfig(; active=false),
            ),
        )

        # CFL and subcycle limits
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true,
                    cfl_settling=0.0,
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true,
                    max_subcycles=0,
                ),
            ),
        )
    end

    @testset "CoreFormationConfig TOML Serialization Roundtrip" begin
        cfg_orig = SimulationConfig(;
            coreformation=CoreFormationConfig(;
                percolation_active=true,
                settling_active=true,
                rho_metal=7800.0,
                eta_metal=5.0e-3,
                phi_pack=0.70,
                T_eutectic=1230.0,
                droplet_size_mode=:fixed,
                droplet_diameter_fixed=2.5e-3,
            ),
            melting=MeltingConfig(; active=true),
        )
        toml_str = save_config(cfg_orig)
        @test occursin("[coreformation]", toml_str)
        @test occursin("percolation_active = true", toml_str)
        @test occursin("settling_active = true", toml_str)

        cfg_loaded = load_config(toml_str)
        @test cfg_loaded.coreformation.percolation_active == true
        @test cfg_loaded.coreformation.settling_active == true
        @test isapprox(cfg_loaded.coreformation.rho_metal, 7800.0; rtol=1e-12)
        @test isapprox(cfg_loaded.coreformation.eta_metal, 5.0e-3; rtol=1e-12)
        @test isapprox(cfg_loaded.coreformation.phi_pack, 0.70; rtol=1e-12)
        @test isapprox(cfg_loaded.coreformation.T_eutectic, 1230.0; rtol=1e-12)
        @test cfg_loaded.coreformation.droplet_size_mode == :fixed
        @test isapprox(cfg_loaded.coreformation.droplet_diameter_fixed, 2.5e-3; rtol=1e-12)
    end

    @testset "Liquid Metal Permeability Physics" begin
        # 1. Error contracts: DomainError on invalid porosity or unphysical parameters
        @test_throws DomainError metal_permeability(-0.01)
        @test_throws DomainError metal_permeability(1.0)
        @test_throws DomainError metal_permeability(0.1; k_metal_ref=-1.0e-9)
        @test_throws DomainError metal_permeability(0.1; phi_crit_perc=-0.01)
        @test_throws DomainError metal_permeability(0.1; phi_crit_perc=1.0)
        @test_throws DomainError metal_permeability(0.1; perm_exponent=-1.0)

        # 2. Sub-threshold and threshold limit: exactly zero permeability
        @test iszero(metal_permeability(0.0; k_metal_ref=1.0e-9, phi_crit_perc=0.05))
        @test iszero(metal_permeability(0.03; k_metal_ref=1.0e-9, phi_crit_perc=0.05))
        @test iszero(metal_permeability(0.05; k_metal_ref=1.0e-9, phi_crit_perc=0.05))

        # 3. Analytical value above threshold (Kozeny-Carman formulation)
        phi_m = 0.15
        phi_crit = 0.05
        k_ref = 1.0e-9
        phi0 = 0.1
        phi_mob = phi_m - phi_crit # 0.10
        # Expected: k_ref * (phi_mob / phi0)^3 * ((1 - phi_mob)/(1 - phi0))^(-2) = 1e-9
        k_val = metal_permeability(phi_m; k_metal_ref=k_ref, phi_crit_perc=phi_crit, phi0=phi0)
        @test isapprox(k_val, 1.0e-9; rtol=1e-10)

        # 4. Strict monotonicity above threshold
        k_low = metal_permeability(0.10; k_metal_ref=k_ref, phi_crit_perc=phi_crit, phi0=phi0)
        k_high = metal_permeability(0.20; k_metal_ref=k_ref, phi_crit_perc=phi_crit, phi0=phi0)
        @test k_low > 0.0
        @test k_high > k_low

        # 5. Prefactor scaling (k_ref doubles permeability)
        k_doubled = metal_permeability(phi_m; k_metal_ref=2.0e-9, phi_crit_perc=phi_crit, phi0=phi0)
        @test isapprox(k_doubled, 2.0 * k_val; rtol=1e-10)
    end

    @testset "Stokes Settling Velocity Physics" begin
        # Baseline parameters
        r_d = 5.0e-3      # 5 mm
        drho = 4000.0     # kg/m^3
        g_acc = 0.05      # m/s^2 (planetesimal gravity)
        eta = 10.0        # Pa s (magma suspension)

        # 1. Analytical check: v = (2/9) * drho * g * r^2 / eta
        # v = (2/9) * 4000 * 0.05 * (2.5e-5) / 10 = (2/9) * 5e-4 = 1.1111111111e-4 m/s
        v_analytical = (2.0 / 9.0) * drho * g_acc * (r_d^2) / eta
        v_calc = stokes_settling_velocity(r_d, drho, g_acc, eta)
        @test isapprox(v_calc, v_analytical; rtol=1e-12)

        # 2. 3-Class Discrimination Guards:
        # Exponent guard: r^2 quadratic scaling
        v_2r = stokes_settling_velocity(2.0 * r_d, drho, g_acc, eta)
        @test isapprox(v_2r, 4.0 * v_analytical; rtol=1e-10)
        @test abs(v_2r - 2.0 * v_analytical) > 1.0e-5

        # Sign guard: positive downward settling
        @test v_calc > 0.0

        # Scale guard: mm/s to cm/s order of magnitude for magma droplets
        @test 1.0e-5 < v_calc < 1.0e-2

        # 3. Degenerate zero limits
        @test iszero(stokes_settling_velocity(0.0, drho, g_acc, eta))
        @test iszero(stokes_settling_velocity(r_d, 0.0, g_acc, eta))
        @test iszero(stokes_settling_velocity(r_d, drho, 0.0, eta))

        # 4. Fluid droplet Hadamard-Rybczynski factor (1.5x)
        v_hr = stokes_settling_velocity(r_d, drho, g_acc, eta; hadamard_rybczynski=true)
        @test isapprox(v_hr, 1.5 * v_analytical; rtol=1e-10)

        # 5. Error contracts
        @test_throws DomainError stokes_settling_velocity(-1.0e-3, drho, g_acc, eta)
        @test_throws DomainError stokes_settling_velocity(r_d, drho, g_acc, 0.0)
        @test_throws DomainError stokes_settling_velocity(r_d, drho, g_acc, -10.0)
        @test_throws DomainError stokes_settling_velocity(NaN, drho, g_acc, eta)
    end

    @testset "Weber Droplet Breakup Equilibrium" begin
        rho_sil = 2800.0
        v_rel = 0.01      # 1 cm/s
        sigma = 1.0       # N/m
        we_crit = 10.0

        # 1. Analytical value: d = We_crit * sigma / (rho_sil * v_rel^2)
        d_analytical = we_crit * sigma / (rho_sil * (v_rel^2))
        d_calc = weber_equilibrium_diameter(rho_sil, v_rel, sigma; We_crit=we_crit)
        @test isapprox(d_calc, d_analytical; rtol=1e-12)

        # 2. Invariant: Weber number at calculated diameter equals We_crit
        we_recalc = rho_sil * (v_rel^2) * d_calc / sigma
        @test isapprox(we_recalc, we_crit; rtol=1e-10)

        # 3. Inverse quadratic scaling on relative velocity
        d_2v = weber_equilibrium_diameter(rho_sil, 2.0 * v_rel, sigma; We_crit=we_crit)
        @test isapprox(d_2v, 0.25 * d_analytical; rtol=1e-10)
        @test d_2v < d_calc

        # 4. Error contracts
        @test_throws DomainError weber_equilibrium_diameter(-2800.0, v_rel, sigma)
        @test_throws DomainError weber_equilibrium_diameter(rho_sil, 0.0, sigma)
        @test_throws DomainError weber_equilibrium_diameter(rho_sil, -0.01, sigma)
        @test_throws DomainError weber_equilibrium_diameter(rho_sil, v_rel, -1.0)
        @test_throws DomainError weber_equilibrium_diameter(rho_sil, v_rel, sigma; We_crit=0.0)
    end

    @testset "Richardson-Zaki Hindered Settling" begin
        # 1. Boundary limits: 1.0 at zero metal, 0.0 at/above packing
        @test isapprox(richardson_zaki_hindrance(0.0; hindered_exponent=4.5, phi_pack=0.65), 1.0; rtol=1e-12)
        @test iszero(richardson_zaki_hindrance(0.65; hindered_exponent=4.5, phi_pack=0.65))
        @test iszero(richardson_zaki_hindrance(0.80; hindered_exponent=4.5, phi_pack=0.65))

        # 2. Analytical value at mid-packing
        phi_m = 0.325 # phi_m / phi_pack = 0.5
        h_expected = 0.5^4.5
        h_calc = richardson_zaki_hindrance(phi_m; hindered_exponent=4.5, phi_pack=0.65)
        @test isapprox(h_calc, h_expected; rtol=1e-10)

        # 3. Monotonic decrease on [0, phi_pack]
        h1 = richardson_zaki_hindrance(0.10; hindered_exponent=4.5, phi_pack=0.65)
        h2 = richardson_zaki_hindrance(0.20; hindered_exponent=4.5, phi_pack=0.65)
        @test h1 > h2 > 0.0

        # 4. Error contracts
        @test_throws DomainError richardson_zaki_hindrance(-0.01)
        @test_throws DomainError richardson_zaki_hindrance(1.05)
        @test_throws DomainError richardson_zaki_hindrance(0.1; hindered_exponent=-1.0)
        @test_throws DomainError richardson_zaki_hindrance(0.1; phi_pack=0.0)
        @test_throws DomainError richardson_zaki_hindrance(0.1; phi_pack=1.5)
    end

    @testset "Eutectic Metal Melt Fraction" begin
        T_eut = 1213.0
        dT_m = 50.0

        # 1. Sub-eutectic: zero melt
        @test iszero(compute_metal_melt_fraction(1100.0; T_eutectic=T_eut, dT_metal=dT_m))
        @test iszero(compute_metal_melt_fraction(1213.0; T_eutectic=T_eut, dT_metal=dT_m))

        # 2. Super-eutectic: full melt
        @test isapprox(compute_metal_melt_fraction(1263.0; T_eutectic=T_eut, dT_metal=dT_m), 1.0; rtol=1e-12)
        @test isapprox(compute_metal_melt_fraction(1400.0; T_eutectic=T_eut, dT_metal=dT_m), 1.0; rtol=1e-12)

        # 3. Mid-window linear ramp
        f_mid = compute_metal_melt_fraction(1238.0; T_eutectic=T_eut, dT_metal=dT_m)
        @test isapprox(f_mid, 0.5; rtol=1e-10)

        # 4. Monotonic increase
        f_low = compute_metal_melt_fraction(1225.0; T_eutectic=T_eut, dT_metal=dT_m)
        f_high = compute_metal_melt_fraction(1250.0; T_eutectic=T_eut, dT_metal=dT_m)
        @test 0.0 < f_low < f_high < 1.0

        # 5. Error contracts
        @test_throws DomainError compute_metal_melt_fraction(-10.0)
        @test_throws DomainError compute_metal_melt_fraction(1200.0; T_eutectic=0.0)
        @test_throws DomainError compute_metal_melt_fraction(1200.0; dT_metal=0.0)
        @test_throws DomainError compute_metal_melt_fraction(NaN)
    end

    @testset "Unified Metal Segregation Velocity & Regimes" begin
        phi_m = 0.15
        drho = 3900.0
        g_acc = 0.05
        eta_s = 10.0

        # 1. Inactive case: returns zero
        v_none = metal_segregation_velocity(
            phi_m, 0.20, drho, g_acc, eta_s;
            percolation_active=false, settling_active=false,
        )
        @test iszero(v_none)

        # 2. Percolation only case:
        v_perc_solid = metal_segregation_velocity(
            phi_m, 0.0, drho, g_acc, eta_s;
            percolation_active=true, settling_active=false,
            phi_crit_perc=0.05, F_perc_end=0.50,
        )
        @test v_perc_solid > 0.0

        # Shuts off at or above disaggregation when settling is inactive
        v_perc_high_melt = metal_segregation_velocity(
            phi_m, 0.55, drho, g_acc, eta_s;
            percolation_active=true, settling_active=false,
            phi_crit_perc=0.05, F_perc_end=0.50,
        )
        @test iszero(v_perc_high_melt)

        # 3. Settling only case:
        v_settle_sub = metal_segregation_velocity(
            phi_m, 0.20, drho, g_acc, eta_s;
            percolation_active=false, settling_active=true,
            F_settle_start=0.40,
        )
        @test iszero(v_settle_sub)

        v_settle_high = metal_segregation_velocity(
            phi_m, 0.60, drho, g_acc, eta_s;
            percolation_active=false, settling_active=true,
            F_settle_start=0.40,
        )
        @test v_settle_high > 0.0

        # 4. Hybrid transition continuity across [0.40, 0.50]
        v_hyb_solid = metal_segregation_velocity(
            phi_m, 0.39, drho, g_acc, eta_s;
            percolation_active=true, settling_active=true,
            F_settle_start=0.40, F_perc_end=0.50,
        )
        v_hyb_mid = metal_segregation_velocity(
            phi_m, 0.45, drho, g_acc, eta_s;
            percolation_active=true, settling_active=true,
            F_settle_start=0.40, F_perc_end=0.50,
        )
        v_hyb_melt = metal_segregation_velocity(
            phi_m, 0.51, drho, g_acc, eta_s;
            percolation_active=true, settling_active=true,
            F_settle_start=0.40, F_perc_end=0.50,
        )
        @test isapprox(v_hyb_solid, v_perc_solid; rtol=1e-10)
        @test isapprox(v_hyb_melt, v_settle_high; rtol=1e-10)
        @test v_hyb_mid > 0.0

        # 5. Residual metal fraction sensitivity
        v_res_low = metal_segregation_velocity(
            phi_m, 0.0, drho, g_acc, eta_s;
            percolation_active=true, settling_active=false,
            phi_crit_perc=0.05, phi_residual=0.01,
        )
        v_res_high = metal_segregation_velocity(
            phi_m, 0.0, drho, g_acc, eta_s;
            percolation_active=true, settling_active=false,
            phi_crit_perc=0.05, phi_residual=0.04,
        )
        @test v_res_low > v_res_high > 0.0

        # Trapped below residual
        v_trapped = metal_segregation_velocity(
            0.03, 0.0, drho, g_acc, eta_s;
            percolation_active=true, settling_active=false,
            phi_crit_perc=0.05, phi_residual=0.04,
        )
        @test iszero(v_trapped)

        # 6. Error contracts
        @test_throws DomainError metal_segregation_velocity(
            -0.01, 0.2, drho, g_acc, eta_s;
            percolation_active=true, settling_active=true,
        )
        @test_throws DomainError metal_segregation_velocity(
            phi_m, 0.2, drho, g_acc, -1.0;
            percolation_active=true, settling_active=true,
        )
        @test_throws DomainError metal_segregation_velocity(
            phi_m, 0.2, drho, g_acc, eta_s;
            percolation_active=true, phi_residual=-0.01,
        )
        @test_throws DomainError metal_segregation_velocity(
            phi_m, 0.2, drho, g_acc, eta_s;
            percolation_active=true, phi_crit_perc=0.05, phi_residual=0.06,
        )
    end

    @testset "Segregation Dissipation Energetics" begin
        phi_m = 0.10
        drho = 3900.0
        g_acc = 0.05
        v_seg = 1.0e-4

        # 1. Analytical value: Q = phi_m * drho * g * v_seg
        Q_analytical = phi_m * drho * g_acc * v_seg
        Q_calc = segregation_dissipation_heating(phi_m, drho, g_acc, v_seg)
        @test isapprox(Q_calc, Q_analytical; rtol=1e-12)

        # 2. Positivity and non-negativity
        @test Q_calc > 0.0

        # 3. Degenerate zero limits
        @test iszero(segregation_dissipation_heating(0.0, drho, g_acc, v_seg))
        @test iszero(segregation_dissipation_heating(phi_m, 0.0, g_acc, v_seg))
        @test iszero(segregation_dissipation_heating(phi_m, drho, 0.0, v_seg))
        @test iszero(segregation_dissipation_heating(phi_m, drho, g_acc, 0.0))

        # 4. Error contracts
        @test_throws DomainError segregation_dissipation_heating(-0.01, drho, g_acc, v_seg)
        @test_throws DomainError segregation_dissipation_heating(1.05, drho, g_acc, v_seg)
        @test_throws DomainError segregation_dissipation_heating(phi_m, -1.0, g_acc, v_seg)
        @test_throws DomainError segregation_dissipation_heating(phi_m, drho, -0.05, v_seg)
        @test_throws DomainError segregation_dissipation_heating(phi_m, drho, g_acc, -1.0e-4)
    end

    @testset "Metal Property Blending Mixtures" begin
        rho_s = 3300.0
        rho_m = 7200.0
        k_s = 3.0
        k_m = 40.0
        rhocp_s = 3.3e6
        rhocp_m = 4.0e6

        # 1. Density mixture limits
        @test isapprox(metal_blended_density(rho_s, rho_m, 0.0), rho_s; rtol=1e-12)
        @test isapprox(metal_blended_density(rho_s, rho_m, 1.0), rho_m; rtol=1e-12)
        @test isapprox(metal_blended_density(rho_s, rho_m, 0.5), 0.5 * (rho_s + rho_m); rtol=1e-12)
        @test metal_blended_density(rho_s, rho_m, 0.2) < metal_blended_density(rho_s, rho_m, 0.4)

        # 2. Thermal conductivity mixture limits
        @test isapprox(metal_blended_conductivity(k_s, k_m, 0.0), k_s; rtol=1e-12)
        @test isapprox(metal_blended_conductivity(k_s, k_m, 1.0), k_m; rtol=1e-12)
        @test isapprox(metal_blended_conductivity(k_s, k_m, 0.5; mode=:arithmetic), 21.5; rtol=1e-12)
        k_geom = metal_blended_conductivity(k_s, k_m, 0.5; mode=:geometric)
        @test isapprox(k_geom, sqrt(k_s * k_m); rtol=1e-12)

        # 3. Heat capacity mixture limits
        @test isapprox(metal_blended_heat_capacity(rhocp_s, rhocp_m, 0.0), rhocp_s; rtol=1e-12)
        @test isapprox(metal_blended_heat_capacity(rhocp_s, rhocp_m, 1.0), rhocp_m; rtol=1e-12)
        @test isapprox(metal_blended_heat_capacity(rhocp_s, rhocp_m, 0.5), 0.5 * (rhocp_s + rhocp_m); rtol=1e-12)

        # 4. Error contracts
        @test_throws DomainError metal_blended_density(-3300.0, rho_m, 0.5)
        @test_throws DomainError metal_blended_density(rho_s, 0.0, 0.5)
        @test_throws DomainError metal_blended_density(rho_s, rho_m, -0.1)
        @test_throws DomainError metal_blended_conductivity(0.0, k_m, 0.5)
        @test_throws DomainError metal_blended_heat_capacity(rhocp_s, -1.0, 0.5)
        @test_throws ArgumentError metal_blended_conductivity(k_s, k_m, 0.5; mode=:unsupported)
    end

    @testset "Suspension Rouse Number Diagnostic" begin
        # Planetesimal parameters: v_st ~ 1e-4 m/s, u_conv ~ 1e-5 m/s -> R > 1
        v_settle = 1.0e-4
        u_conv = 1.0e-5
        r_rouse = suspension_rouse_number(v_settle, u_conv)
        @test isapprox(r_rouse, 10.0; rtol=1e-12)
        @test r_rouse > 1.0

        # Vigorous magma ocean: v_st ~ 1e-4 m/s, u_conv ~ 1e-2 m/s -> R < 1
        r_vigorous = suspension_rouse_number(v_settle, 1.0e-2)
        @test isapprox(r_vigorous, 0.01; rtol=1e-12)
        @test r_vigorous < 1.0

        # Zero settling limit
        @test iszero(suspension_rouse_number(0.0, u_conv))

        # Error contracts
        @test_throws DomainError suspension_rouse_number(-1.0e-4, u_conv)
        @test_throws DomainError suspension_rouse_number(v_settle, 0.0)
        @test_throws DomainError suspension_rouse_number(v_settle, -1.0e-3)
    end

    @testset "Marker Metal Property Setup & Allocation" begin
        marknum = 300
        (Xfem, Xfem0, Xfe_bulk) = setup_marker_metal_properties(marknum)
        @test length(Xfem) == marknum
        @test length(Xfem0) == marknum
        @test length(Xfe_bulk) == marknum
        @test eltype(Xfem) == Float64
        @test eltype(Xfe_bulk) == Float64
        @test all(iszero, Xfem)
        @test all(iszero, Xfe_bulk)

        # Extended tuple unpack from setup_marker_properties
        extended_props = setup_marker_properties(marknum; include_metal=true)
        @test length(extended_props) == 16
        @test length(extended_props[14]) == marknum # Xfem
        @test length(extended_props[15]) == marknum # Xfem0
        @test length(extended_props[16]) == marknum # Xfe_bulk
    end

    @testset "define_markers! Metal Tracking & Zoning" begin
        coords = GridCoordinates(GridConfig(; Nx=20, Ny=20, xsize=140000.0, ysize=140000.0))
        marknum = coords.Nxm * coords.Nym
        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0, Fm) = setup_marker_properties(
            marknum, coords
        )
        (rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm, tkm_rhocptotalm, etafluidcur_inv_kphim, inv_gggtotalm, fricttotalm, cohestotalm, tenstotalm, rhofluidcur, alphasolidcur, alphafluidcur) = setup_marker_properties_helpers(
            marknum
        )
        (Xfem, Xfem0, Xfe_bulk) = setup_marker_metal_properties(marknum)

        # Initialize markers with metal tracking
        define_markers!(
            xm, ym, tm, phim, etavpm, rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm,
            tkm, inv_gggtotalm, fricttotalm, cohestotalm, tenstotalm, rhofluidcur, alphasolidcur, alphafluidcur,
            XWsolidm0;
            randomized=false, coords=coords,
            Xfe_bulk=Xfe_bulk, Xfem=Xfem, Xfem0=Xfem0,
            Xfe_bulk_val=0.22, T_eutectic_val=1213.0, dT_metal_val=50.0,
        )

        for m in 1:marknum
            rmark = sqrt((xm[m] - 70000.0)^2 + (ym[m] - 70000.0)^2)
            if rmark < 50000.0 # planet interior
                @test isapprox(Xfe_bulk[m], 0.22; rtol=1e-12)
                # At initial ambient temperature (~273-300 K), metal is solid (sub-eutectic)
                @test iszero(Xfem[m])
                @test iszero(Xfem0[m])
            else # sticky space air
                @test iszero(Xfe_bulk[m])
                @test iszero(Xfem[m])
            end
        end
    end

    @testset "compute_marker_properties! Metal Blending & Melting" begin
        marknum = 10
        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0, Fm) = setup_marker_properties(marknum)
        (rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm, tkm_rhocptotalm, etafluidcur_inv_kphim, inv_gggtotalm, fricttotalm, cohestotalm, tenstotalm, rhofluidcur, alphasolidcur, alphafluidcur) = setup_marker_properties_helpers(marknum)
        (Xfem, Xfem0, Xfe_bulk) = setup_marker_metal_properties(marknum)

        # Setup planet rock marker (tm = 1) with bulk metal
        m = 1
        tm[m] = 1
        XWsolidm0[m] = 0.0
        phim[m] = 1.0e-4
        Xfe_bulk[m] = 0.25

        # 1. Below eutectic: metal is solid, no liquid metal
        tkm[m] = 1100.0
        compute_marker_properties!(
            m, tm, tkm, rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm, tkm_rhocptotalm,
            etafluidcur_inv_kphim, start_hrsolidm, start_hrfluidm, phim, XWsolidm0, 9, rhofluidcur;
            Xfe_bulk=Xfe_bulk, Xfem=Xfem, coreformation_active=true,
            T_eutectic_val=1213.0, dT_metal_val=50.0, rho_metal_val=7200.0, k_metal_val=40.0, rhocp_metal_val=4.0e6,
        )
        @test iszero(Xfem[m])
        # Density should equal pure rock density (3300)
        @test isapprox(rhototalm[m], 3300.0; rtol=1e-4)

        # 2. Above eutectic: metal melts and blends density and conductivity
        tkm[m] = 1300.0 # above T_eutectic + dT_metal (1263 K) -> full melt fraction 1.0
        compute_marker_properties!(
            m, tm, tkm, rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm, tkm_rhocptotalm,
            etafluidcur_inv_kphim, start_hrsolidm, start_hrfluidm, phim, XWsolidm0, 9, rhofluidcur;
            Xfe_bulk=Xfe_bulk, Xfem=Xfem, coreformation_active=true,
            T_eutectic_val=1213.0, dT_metal_val=50.0, rho_metal_val=7200.0, k_metal_val=40.0, rhocp_metal_val=4.0e6,
        )
        @test isapprox(Xfem[m], 0.25; rtol=1e-12)
        # Expected blended density: 0.75 * 3300 + 0.25 * 7200 = 2475 + 1800 = 4275 kg/m^3
        @test isapprox(rhototalm[m], 4275.0; rtol=1e-4)
        @test rhototalm[m] > 3300.0
        # Conductivity should increase above solid rock value (~3.0) toward metal (40.0)
        @test ktotalm[m] > 3.0
    end

    @testset "Xfem0 Timestep Advancing" begin
        marknum = 20
        (Xfem, Xfem0, Xfe_bulk) = setup_marker_metal_properties(marknum)
        fill!(Xfem0, 0.0)
        fill!(Xfem, 0.15)
        @test all(iszero, Xfem0)
        # Advance previous timestep molten metal fraction
        Xfem0 .= Xfem
        @test all(x -> isapprox(x, 0.15; rtol=1e-12), Xfem0)
        @test all(iszero, Xfem0 .- Xfem)
    end

    @testset "replenish_markers! Metal Inventory Conservation" begin
        # Setup population tracking with 25 markers in a small grid
        coords = GridCoordinates(GridConfig(; Nx=5, Ny=5, xsize=10000.0, ysize=10000.0))
        mdis, mnum = setup_marker_geometry_helpers(coords)
        marknum = 25
        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0, Fm) = setup_marker_properties(marknum, coords)
        (rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm, tkm_rhocptotalm, etafluidcur_inv_kphim, inv_gggtotalm, fricttotalm, cohestotalm, tenstotalm, rhofluidcur, alphasolidcur, alphafluidcur) = setup_marker_properties_helpers(marknum)
        (Xfem, Xfem0, Xfe_bulk) = setup_marker_metal_properties(marknum)

        # Place markers at positions with non-uniform distinct metal values
        for i in 1:marknum
            xm[i] = 1000.0 + (i - 1) * 300.0
            ym[i] = 1000.0 + (i - 1) * 300.0
            tm[i] = 1
            tkm[i] = 300.0
            Xfe_bulk[i] = 0.10 + 0.005 * i
            Xfem[i] = 0.02 + 0.001 * i
            Xfem0[i] = 0.02 + 0.001 * i
        end

        marknum_new = replenish_markers!(
            xm, ym, tm, tkm, phim, sxxm, sxym, etavpm, phinewm, pfm0, XWsolidm, XWsolidm0,
            rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm, inv_gggtotalm, fricttotalm,
            cohestotalm, tenstotalm, rhofluidcur, alphasolidcur, alphafluidcur, tkm_rhocptotalm,
            etafluidcur_inv_kphim, mdis, mnum;
            Fm=Fm, randomized=false, coords=coords,
            Xfem=Xfem, Xfem0=Xfem0, Xfe_bulk=Xfe_bulk,
        )

        @test marknum_new > marknum
        @test length(Xfem) == marknum_new
        @test length(Xfem0) == marknum_new
        @test length(Xfe_bulk) == marknum_new
        # Re-seeded markers must strictly inherit parent values (not an average)
        original_bulk_set = Set(Xfe_bulk[1:marknum])
        original_melt_set = Set(Xfem[1:marknum])
        for k in (marknum + 1):marknum_new
            @test Xfe_bulk[k] in original_bulk_set
            @test Xfem[k] in original_melt_set
            @test isapprox(Xfem0[k], Xfem[k]; rtol=1e-12)
        end
    end

    @testset "Checkpoint Serialization Roundtrip with Metal Fields" begin
        marknum = 50
        (Xfem, Xfem0, Xfe_bulk) = setup_marker_metal_properties(marknum)
        fill!(Xfe_bulk, 0.24)
        fill!(Xfem, 0.12)
        fill!(Xfem0, 0.12)

        ckpt_file = joinpath(tempdir(), "test_core_ckpt_$(rand(1000:9999)).jld2")
        try
            JLD2.jldsave(ckpt_file; Xfem=Xfem, Xfe_bulk=Xfe_bulk, Xfem0=Xfem0)
            loaded = load_state(ckpt_file)
            @test haskey(loaded, "Xfem")
            @test haskey(loaded, "Xfe_bulk")
            @test isapprox(loaded["Xfe_bulk"][1], 0.24; rtol=1e-12)
            @test isapprox(loaded["Xfem"][1], 0.12; rtol=1e-12)
        finally
            rm(ckpt_file; force=true)
        end

        # Test simulation_loop execution and checkpoint persistence with core formation
        output_dir = mktempdir()
        try
            quick_toml = joinpath(@__DIR__, "..", "configs", "test_quick.toml")
            cfg = load_config(quick_toml)
            cfg_core = SimulationConfig(
                grid=cfg.grid,
                time=TimeConfig(n_steps=1, dt_initial=cfg.time.dt_initial),
                solver=cfg.solver,
                geometry=cfg.geometry,
                materials=cfg.materials,
                thermodynamics=cfg.thermodynamics,
                reaction=cfg.reaction,
                melting=cfg.melting,
                coreformation=CoreFormationConfig(
                    percolation_active=true,
                    settling_active=true,
                    Xfe_bulk=0.22,
                ),
                output=OutputConfig(savematstep=1, output_dir=output_dir),
            )
            Erebus.simulation_loop(cfg_core; output_path=output_dir)
            out_file0 = joinpath(output_dir, "output_00000.jld2")
            out_file1 = joinpath(output_dir, "output_00001.jld2")
            @test isfile(out_file0)
            @test isfile(out_file1)
            loaded0 = load_state(out_file0)
            @test haskey(loaded0, "Xfem")
            @test haskey(loaded0, "Xfe_bulk")
            @test haskey(loaded0, "Xfem0")
            loaded1 = load_state(out_file1)
            @test haskey(loaded1, "Xfem")
            @test haskey(loaded1, "Xfe_bulk")
            @test haskey(loaded1, "Xfem0")

            # Verify backward compatibility: run without core formation does not save Xfem
            output_dir_nocore = mktempdir()
            try
                cfg_nocore = SimulationConfig(
                    grid=cfg.grid,
                    time=TimeConfig(n_steps=0, dt_initial=cfg.time.dt_initial),
                    solver=cfg.solver,
                    geometry=cfg.geometry,
                    materials=cfg.materials,
                    thermodynamics=cfg.thermodynamics,
                    reaction=cfg.reaction,
                    melting=cfg.melting,
                    coreformation=CoreFormationConfig(
                        percolation_active=false,
                        settling_active=false,
                    ),
                    output=OutputConfig(savematstep=1, output_dir=output_dir_nocore),
                )
                Erebus.simulation_loop(cfg_nocore; output_path=output_dir_nocore)
                out_nocore = joinpath(output_dir_nocore, "output_00000.jld2")
                @test isfile(out_nocore)
                loaded_nocore = load_state(out_nocore)
                @test !haskey(loaded_nocore, "Xfem")
                @test !haskey(loaded_nocore, "Xfe_bulk")
            finally
                rm(output_dir_nocore, recursive=true, force=true)
            end
        finally
            rm(output_dir, recursive=true, force=true)
        end
    end
end

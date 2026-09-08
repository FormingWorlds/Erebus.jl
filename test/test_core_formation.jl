using Test
using Erebus
using Erebus.Config
using Erebus.Physics
using Erebus.Particles
using Erebus.Geometry
using Erebus: start_hrsolidm, start_hrfluidm
using Erebus.Numerics: assemble_thermal_lse!
using StaticArrays
using TOML
using JLD2

@testset "Iron Core Formation Physics & Configuration" begin
    @testset "CoreFormationConfig Schema & Defaults" begin
        cfg_def = CoreFormationConfig()
        @test cfg_def.percolation_active == false
        @test cfg_def.settling_active == false
        @test isapprox(cfg_def.sulfur_fraction, 0.31; rtol=1e-12)
        @test cfg_def.metal_density_mode == :sanloup2000
        @test isapprox(cfg_def.rho_metal, 5450.0; rtol=1e-12)
        @test isapprox(cfg_def.rho_metal_solid, 5700.0; rtol=1e-12)
        @test isapprox(cfg_def.L_metal, 2.7e5; rtol=1e-12)
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
        @test isapprox(cfg_def.phi0, 0.1; rtol=1e-12)
        @test cfg_def.droplet_size_mode == :capillary_mean
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
        # Sulfur fraction bounds: must be in [0, 0.40] and finite
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true, sulfur_fraction=-0.05
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true, sulfur_fraction=0.45
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true, sulfur_fraction=NaN
                ),
            ),
        )

        # Metal density mode validation: must be :sanloup2000, :morard2014, or :constant
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true, metal_density_mode=:invalid_density_mode
                ),
            ),
        )

        # Metal density must exceed silicate density (rhosolidm[1] = 3300) and be finite
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true, rho_metal=3000.0
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(; percolation_active=true, rho_metal=Inf)
            ),
        )

        # Solid metal density must be finite and strictly > rho_metal
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true, rho_metal=5450.0, rho_metal_solid=5450.0
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true, rho_metal=5450.0, rho_metal_solid=5400.0
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true, rho_metal_solid=Inf
                ),
            ),
        )

        # L_metal must be non-negative and finite
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(; percolation_active=true, L_metal=-1.0e5)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(; percolation_active=true, L_metal=NaN)
            ),
        )

        # Non-positive viscosity or thermal properties
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true, eta_metal=-1.0e-2
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(; settling_active=true, k_metal=0.0),
                melting=MeltingConfig(; active=true),
            ),
        )

        # Percolation ordering: phi_residual <= phi_crit_perc < phi_pack
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true, phi_residual=0.10, phi_crit_perc=0.05
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true, phi_crit_perc=0.70, phi_pack=0.65
                ),
            ),
        )

        # Regime handover interval: F_settle_start <= F_perc_end
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    settling_active=true, F_settle_start=0.60, F_perc_end=0.40
                ),
                melting=MeltingConfig(; active=true),
            ),
        )

        # Xfe_bulk cannot exceed phi_pack
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true, Xfe_bulk=0.75, phi_pack=0.65
                ),
            ),
        )

        # phi0 must be in (0, 1)
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(; percolation_active=true, phi0=-0.05)
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(; percolation_active=true, phi0=1.05)
            ),
        )

        # Droplet mode symbol
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    settling_active=true, droplet_size_mode=:invalid_mode
                ),
                melting=MeltingConfig(; active=true),
            ),
        )

        # Eutectic temperature must be below silicate solidus for percolation
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true, T_eutectic=1500.0
                ),
            ),
        )

        # Settling requires melting active
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(; settling_active=true),
                melting=MeltingConfig(; active=false),
            ),
        )

        # CFL and subcycle limits
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true, cfl_settling=0.0
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                coreformation=CoreFormationConfig(;
                    percolation_active=true, max_subcycles=0
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
                rho_metal_solid=8200.0,
                L_metal=2.5e5,
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
        @test isapprox(cfg_loaded.coreformation.rho_metal_solid, 8200.0; rtol=1e-12)
        @test isapprox(cfg_loaded.coreformation.L_metal, 2.5e5; rtol=1e-12)
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
        k_val = metal_permeability(
            phi_m; k_metal_ref=k_ref, phi_crit_perc=phi_crit, phi0=phi0
        )
        @test isapprox(k_val, 1.0e-9; rtol=1e-10)

        # 4. Strict monotonicity above threshold
        k_low = metal_permeability(
            0.10; k_metal_ref=k_ref, phi_crit_perc=phi_crit, phi0=phi0
        )
        k_high = metal_permeability(
            0.20; k_metal_ref=k_ref, phi_crit_perc=phi_crit, phi0=phi0
        )
        @test k_low > 0.0
        @test k_high > k_low

        # 5. Prefactor scaling (k_ref doubles permeability)
        k_doubled = metal_permeability(
            phi_m; k_metal_ref=2.0e-9, phi_crit_perc=phi_crit, phi0=phi0
        )
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

        # 4. Fluid droplet Hadamard-Rybczynski factor
        v_hr_inviscid = stokes_settling_velocity(
            r_d, drho, g_acc, eta; hadamard_rybczynski=true
        )
        @test isapprox(v_hr_inviscid, 1.5 * v_analytical; rtol=1e-10)
        # Viscosity ratio: (3*eta_s + 3*eta_m) / (2*eta_s + 3*eta_m) = (30 + 3) / (20 + 3) = 33/23
        v_hr_ratio = stokes_settling_velocity(
            r_d, drho, g_acc, eta; hadamard_rybczynski=true, eta_metal=1.0
        )
        @test isapprox(v_hr_ratio, (33.0 / 23.0) * v_analytical; rtol=1e-10)

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
        @test_throws DomainError weber_equilibrium_diameter(
            rho_sil, v_rel, sigma; We_crit=0.0
        )
    end

    @testset "Richardson-Zaki Hindered Settling" begin
        # 1. Boundary limits: 1.0 at zero metal, 0.0 at/above packing
        @test isapprox(
            richardson_zaki_hindrance(0.0; hindered_exponent=4.5, phi_pack=0.65),
            1.0;
            rtol=1e-12,
        )
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
        @test isapprox(
            compute_metal_melt_fraction(1263.0; T_eutectic=T_eut, dT_metal=dT_m),
            1.0;
            rtol=1e-12,
        )
        @test isapprox(
            compute_metal_melt_fraction(1400.0; T_eutectic=T_eut, dT_metal=dT_m),
            1.0;
            rtol=1e-12,
        )

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
            phi_m, 0.20, drho, g_acc, eta_s; percolation_active=false, settling_active=false
        )
        @test iszero(v_none)

        # 2. Percolation only case:
        v_perc_solid = metal_segregation_velocity(
            phi_m,
            0.0,
            drho,
            g_acc,
            eta_s;
            percolation_active=true,
            settling_active=false,
            phi_crit_perc=0.05,
            F_perc_end=0.50,
        )
        @test v_perc_solid > 0.0

        # Shuts off at or above disaggregation when settling is inactive
        v_perc_high_melt = metal_segregation_velocity(
            phi_m,
            0.55,
            drho,
            g_acc,
            eta_s;
            percolation_active=true,
            settling_active=false,
            phi_crit_perc=0.05,
            F_perc_end=0.50,
        )
        @test iszero(v_perc_high_melt)

        # 3. Settling only case:
        v_settle_sub = metal_segregation_velocity(
            phi_m,
            0.20,
            drho,
            g_acc,
            eta_s;
            percolation_active=false,
            settling_active=true,
            F_settle_start=0.40,
        )
        @test iszero(v_settle_sub)

        v_settle_high = metal_segregation_velocity(
            phi_m,
            0.60,
            drho,
            g_acc,
            eta_s;
            percolation_active=false,
            settling_active=true,
            F_settle_start=0.40,
        )
        @test v_settle_high > 0.0

        # 4. Hybrid transition continuity across [0.40, 0.50]
        v_hyb_solid = metal_segregation_velocity(
            phi_m,
            0.39,
            drho,
            g_acc,
            eta_s;
            percolation_active=true,
            settling_active=true,
            F_settle_start=0.40,
            F_perc_end=0.50,
        )
        v_hyb_mid = metal_segregation_velocity(
            phi_m,
            0.45,
            drho,
            g_acc,
            eta_s;
            percolation_active=true,
            settling_active=true,
            F_settle_start=0.40,
            F_perc_end=0.50,
        )
        v_hyb_melt = metal_segregation_velocity(
            phi_m,
            0.51,
            drho,
            g_acc,
            eta_s;
            percolation_active=true,
            settling_active=true,
            F_settle_start=0.40,
            F_perc_end=0.50,
        )
        @test isapprox(v_hyb_solid, v_perc_solid; rtol=1e-10)
        @test isapprox(v_hyb_melt, v_settle_high; rtol=1e-10)
        @test v_hyb_mid > 0.0

        # 5. Residual metal fraction sensitivity
        v_res_low = metal_segregation_velocity(
            phi_m,
            0.0,
            drho,
            g_acc,
            eta_s;
            percolation_active=true,
            settling_active=false,
            phi_crit_perc=0.05,
            phi_residual=0.01,
        )
        v_res_high = metal_segregation_velocity(
            phi_m,
            0.0,
            drho,
            g_acc,
            eta_s;
            percolation_active=true,
            settling_active=false,
            phi_crit_perc=0.05,
            phi_residual=0.04,
        )
        @test v_res_low > v_res_high > 0.0

        # Trapped below residual
        v_trapped = metal_segregation_velocity(
            0.03,
            0.0,
            drho,
            g_acc,
            eta_s;
            percolation_active=true,
            settling_active=false,
            phi_crit_perc=0.05,
            phi_residual=0.04,
        )
        @test iszero(v_trapped)

        # 6. Error contracts
        @test_throws DomainError metal_segregation_velocity(
            -0.01, 0.2, drho, g_acc, eta_s; percolation_active=true, settling_active=true
        )
        @test_throws DomainError metal_segregation_velocity(
            phi_m, 0.2, drho, g_acc, -1.0; percolation_active=true, settling_active=true
        )
        @test_throws DomainError metal_segregation_velocity(
            phi_m, 0.2, drho, g_acc, eta_s; percolation_active=true, phi_residual=-0.01
        )
        @test_throws DomainError metal_segregation_velocity(
            phi_m,
            0.2,
            drho,
            g_acc,
            eta_s;
            percolation_active=true,
            phi_crit_perc=0.05,
            phi_residual=0.06,
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
        @test_throws DomainError segregation_dissipation_heating(
            phi_m, drho, g_acc, -1.0e-4
        )
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
        @test isapprox(
            metal_blended_density(rho_s, rho_m, 0.5), 0.5 * (rho_s + rho_m); rtol=1e-12
        )
        @test metal_blended_density(rho_s, rho_m, 0.2) <
            metal_blended_density(rho_s, rho_m, 0.4)

        # 2. Thermal conductivity mixture limits
        @test isapprox(metal_blended_conductivity(k_s, k_m, 0.0), k_s; rtol=1e-12)
        @test isapprox(metal_blended_conductivity(k_s, k_m, 1.0), k_m; rtol=1e-12)
        @test isapprox(
            metal_blended_conductivity(k_s, k_m, 0.5; mode=:arithmetic), 21.5; rtol=1e-12
        )
        k_geom = metal_blended_conductivity(k_s, k_m, 0.5; mode=:geometric)
        @test isapprox(k_geom, sqrt(k_s * k_m); rtol=1e-12)

        # 3. Heat capacity mixture limits
        @test isapprox(
            metal_blended_heat_capacity(rhocp_s, rhocp_m, 0.0), rhocp_s; rtol=1e-12
        )
        @test isapprox(
            metal_blended_heat_capacity(rhocp_s, rhocp_m, 1.0), rhocp_m; rtol=1e-12
        )
        @test isapprox(
            metal_blended_heat_capacity(rhocp_s, rhocp_m, 0.5),
            0.5 * (rhocp_s + rhocp_m);
            rtol=1e-12,
        )

        # 4. Error contracts
        @test_throws DomainError metal_blended_density(-3300.0, rho_m, 0.5)
        @test_throws DomainError metal_blended_density(rho_s, 0.0, 0.5)
        @test_throws DomainError metal_blended_density(rho_s, rho_m, -0.1)
        @test_throws DomainError metal_blended_conductivity(0.0, k_m, 0.5)
        @test_throws DomainError metal_blended_heat_capacity(rhocp_s, -1.0, 0.5)
        @test_throws ArgumentError metal_blended_conductivity(
            k_s, k_m, 0.5; mode=:unsupported
        )
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
            xm,
            ym,
            tm,
            phim,
            etavpm,
            rhototalm,
            rhocptotalm,
            etatotalm,
            hrtotalm,
            ktotalm,
            tkm,
            inv_gggtotalm,
            fricttotalm,
            cohestotalm,
            tenstotalm,
            rhofluidcur,
            alphasolidcur,
            alphafluidcur,
            XWsolidm0;
            randomized=false,
            coords=coords,
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            Xfem0=Xfem0,
            Xfe_bulk_val=0.22,
            T_eutectic_val=1213.0,
            dT_metal_val=50.0,
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
        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0, Fm) = setup_marker_properties(
            marknum
        )
        (rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm, tkm_rhocptotalm, etafluidcur_inv_kphim, inv_gggtotalm, fricttotalm, cohestotalm, tenstotalm, rhofluidcur, alphasolidcur, alphafluidcur) = setup_marker_properties_helpers(
            marknum
        )
        (Xfem, Xfem0, Xfe_bulk) = setup_marker_metal_properties(marknum)

        # Setup planet rock marker (tm = 1) with bulk metal
        m = 1
        tm[m] = 1
        XWsolidm0[m] = 0.0
        phim[m] = 1.0e-4
        Xfe_bulk[m] = 0.25

        # 1. Below eutectic: metal is solid, carrying solid mass, no liquid metal
        tkm[m] = 1100.0
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
            start_hrsolidm,
            start_hrfluidm,
            phim,
            XWsolidm0,
            9,
            rhofluidcur;
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            coreformation_active=true,
            metal_density_mode_val=:constant,
            T_eutectic_val=1213.0,
            dT_metal_val=50.0,
            rho_metal_val=7200.0,
            rho_metal_solid_val=7800.0,
            L_metal_val=2.7e5,
            k_metal_val=40.0,
            rhocp_metal_val=4.0e6,
        )
        @test iszero(Xfem[m])
        # Density should equal blended solid metal + rock density (0.75 * 3300 + 0.25 * 7800 = 4425 kg/m^3)
        @test isapprox(rhototalm[m], 4425.0; rtol=1e-3)
        @test rhototalm[m] > 3300.0
        # Volumetric heat capacity is unbuffered below melting
        @test isapprox(rhocptotalm[m], 0.75 * 3.3e6 + 0.25 * 4.0e6; rtol=1e-2)

        # 2. In melting range: partial melt with latent heat buffering
        tkm[m] = 1238.0 # 50% through 50 K melting interval
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
            start_hrsolidm,
            start_hrfluidm,
            phim,
            XWsolidm0,
            9,
            rhofluidcur;
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            coreformation_active=true,
            metal_density_mode_val=:constant,
            T_eutectic_val=1213.0,
            dT_metal_val=50.0,
            rho_metal_val=7200.0,
            rho_metal_solid_val=7800.0,
            L_metal_val=2.7e5,
            k_metal_val=40.0,
            rhocp_metal_val=4.0e6,
        )
        @test isapprox(Xfem[m], 0.125; rtol=1e-4)
        # Blended density: 0.75 * 3300 + 0.25 * 7500 = 4350 kg/m^3
        @test isapprox(rhototalm[m], 4350.0; rtol=1e-3)
        # Apparent heat capacity buffered by metal latent heat:
        # 0.75 * 3.3e6 + 0.25 * (4.0e6 + 7800.0 * 2.7e5 / 50.0) = 2.475e6 + 1.153e7 = 1.4005e7 J/(m^3 K)
        @test isapprox(
            rhocptotalm[m], 0.75 * 3.3e6 + 0.25 * (4.0e6 + 7800.0 * 2.7e5 / 50.0); rtol=1e-2
        )

        # 3. Above eutectic: fully molten metal blends liquid density and conductivity
        tkm[m] = 1300.0 # above T_eutectic + dT_metal (1263 K) -> full melt fraction 1.0
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
            start_hrsolidm,
            start_hrfluidm,
            phim,
            XWsolidm0,
            9,
            rhofluidcur;
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            coreformation_active=true,
            metal_density_mode_val=:constant,
            T_eutectic_val=1213.0,
            dT_metal_val=50.0,
            rho_metal_val=7200.0,
            rho_metal_solid_val=7800.0,
            L_metal_val=2.7e5,
            k_metal_val=40.0,
            rhocp_metal_val=4.0e6,
        )
        @test isapprox(Xfem[m], 0.25; rtol=1e-12)
        # Expected blended density: 0.75 * 3300 + 0.25 * 7200 = 2475 + 1800 = 4275 kg/m^3
        @test isapprox(rhototalm[m], 4275.0; rtol=1e-4)
        @test rhototalm[m] > 3300.0
        # Conductivity should increase above solid rock value (~3.0) toward metal (40.0)
        @test ktotalm[m] > 3.0
        # Heat capacity returns to unbuffered mixture
        @test isapprox(rhocptotalm[m], 0.75 * 3.3e6 + 0.25 * 4.0e6; rtol=1e-2)

        # 4. Variable EOS liquid metal density (Sanloup 2000 and Morard 2014)
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
            start_hrsolidm,
            start_hrfluidm,
            phim,
            XWsolidm0,
            9,
            rhofluidcur;
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            coreformation_active=true,
            metal_density_mode_val=:sanloup2000,
            sulfur_fraction_val=0.31,
            T_eutectic_val=1213.0,
            dT_metal_val=50.0,
        )
        rho_liq_sanloup = compute_liquid_metal_density(0.31; T=1300.0, law=:sanloup2000)
        @test isapprox(rhototalm[m], 0.75 * 3300.0 + 0.25 * rho_liq_sanloup; rtol=1e-4)

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
            start_hrsolidm,
            start_hrfluidm,
            phim,
            XWsolidm0,
            9,
            rhofluidcur;
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            coreformation_active=true,
            metal_density_mode_val=:morard2014,
            sulfur_fraction_val=0.31,
            T_eutectic_val=1213.0,
            dT_metal_val=50.0,
        )
        rho_liq_morard = compute_liquid_metal_density(0.31; T=1300.0, law=:morard2014)
        @test isapprox(rhototalm[m], 0.75 * 3300.0 + 0.25 * rho_liq_morard; rtol=1e-4)

        # 5. Radiogenic heat deposition of 60Fe into metallic iron phase
        hrmetal_vec = @SVector [8.0e-5, 8.0e-5, 0.0]
        hrtotalm_base = start_hrsolidm[tm[m]]
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
            start_hrsolidm,
            start_hrfluidm,
            phim,
            XWsolidm0,
            9,
            rhofluidcur;
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            coreformation_active=true,
            hrmetalm=hrmetal_vec,
            T_eutectic_val=1213.0,
            dT_metal_val=50.0,
        )
        phi_fe_val = Xfe_bulk[m]
        expected_hr = (1.0 - phi_fe_val) * hrtotalm_base + phi_fe_val * hrmetal_vec[tm[m]]
        @test isapprox(hrtotalm[m], expected_hr; rtol=1e-10)
        @test hrtotalm[m] > hrtotalm_base

        # 6. Radiogenic heat deposition of 60Fe into crust markers (tm=2)
        tm[m] = 2 # crust rock
        hrtotalm_base_crust = start_hrsolidm[tm[m]]
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
            start_hrsolidm,
            start_hrfluidm,
            phim,
            XWsolidm0,
            9,
            rhofluidcur;
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            coreformation_active=true,
            hrmetalm=hrmetal_vec,
            T_eutectic_val=1213.0,
            dT_metal_val=50.0,
        )
        expected_hr_crust =
            (1.0 - phi_fe_val) * hrtotalm_base_crust + phi_fe_val * hrmetal_vec[tm[m]]
        @test isapprox(hrtotalm[m], expected_hr_crust; rtol=1e-10)
        @test hrtotalm[m] > 0.0
        tm[m] = 1 # restore
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
        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0, Fm) = setup_marker_properties(
            marknum, coords
        )
        (rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm, tkm_rhocptotalm, etafluidcur_inv_kphim, inv_gggtotalm, fricttotalm, cohestotalm, tenstotalm, rhofluidcur, alphasolidcur, alphafluidcur) = setup_marker_properties_helpers(
            marknum
        )
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
            xm,
            ym,
            tm,
            tkm,
            phim,
            sxxm,
            sxym,
            etavpm,
            phinewm,
            pfm0,
            XWsolidm,
            XWsolidm0,
            rhototalm,
            rhocptotalm,
            etatotalm,
            hrtotalm,
            ktotalm,
            inv_gggtotalm,
            fricttotalm,
            cohestotalm,
            tenstotalm,
            rhofluidcur,
            alphasolidcur,
            alphafluidcur,
            tkm_rhocptotalm,
            etafluidcur_inv_kphim,
            mdis,
            mnum;
            Fm=Fm,
            randomized=false,
            coords=coords,
            Xfem=Xfem,
            Xfem0=Xfem0,
            Xfe_bulk=Xfe_bulk,
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
                    percolation_active=true, settling_active=true, Xfe_bulk=0.22
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
                        percolation_active=false, settling_active=false
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

    @testset "Sub-Cycled Drift-Flux Segregation Solver" begin
        # Setup grid and markers for segregation tests
        Nx_test = 16
        Ny_test = 16
        coords = GridCoordinates(
            GridConfig(; Nx=Nx_test, Ny=Ny_test, xsize=140000.0, ysize=140000.0)
        )
        marknum = coords.Nxm * coords.Nym
        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0, Fm) = setup_marker_properties(
            marknum, coords
        )
        (rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm, tkm_rhocptotalm, etafluidcur_inv_kphim, inv_gggtotalm, fricttotalm, cohestotalm, tenstotalm, rhofluidcur, alphasolidcur, alphafluidcur) = setup_marker_properties_helpers(
            marknum
        )
        (Xfem, Xfem0, Xfe_bulk) = setup_marker_metal_properties(marknum)

        define_markers!(
            xm,
            ym,
            tm,
            phim,
            etavpm,
            rhototalm,
            rhocptotalm,
            etatotalm,
            hrtotalm,
            ktotalm,
            tkm,
            inv_gggtotalm,
            fricttotalm,
            cohestotalm,
            tenstotalm,
            rhofluidcur,
            alphasolidcur,
            alphafluidcur,
            XWsolidm0;
            randomized=false,
            coords=coords,
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            Xfem0=Xfem0,
            Xfe_bulk_val=0.20,
            T_eutectic_val=1213.0,
            dT_metal_val=50.0,
        )

        @testset "Inactive and Early Guard Conditions" begin
            cfg_inactive = CoreFormationConfig(
                percolation_active=false, settling_active=false
            )
            res_inactive = apply_metal_segregation!(
                xm,
                ym,
                tm,
                tkm,
                phim,
                Xfe_bulk,
                Xfem,
                marknum,
                1.0e10,
                cfg_inactive;
                coords=coords,
            )
            @test iszero(res_inactive.max_v_seg)
            @test iszero(res_inactive.n_subcycles)
            @test iszero(res_inactive.total_dissipation_energy)

            # Zero dt guard
            cfg_active = CoreFormationConfig(percolation_active=true, settling_active=false)
            res_zerodt = apply_metal_segregation!(
                xm,
                ym,
                tm,
                tkm,
                phim,
                Xfe_bulk,
                Xfem,
                marknum,
                0.0,
                cfg_active;
                coords=coords,
            )
            @test iszero(res_zerodt.max_v_seg)
            @test iszero(res_zerodt.n_subcycles)
        end

        @testset "Frozen State Invariance (Below Eutectic)" begin
            # Set temperature below eutectic (1000 K < 1213 K)
            fill!(tkm, 1000.0)
            fill!(Xfem, 0.0)
            initial_sum = sum(Xfe_bulk)

            cfg_perc = CoreFormationConfig(percolation_active=true, settling_active=false)
            res_cold = apply_metal_segregation!(
                xm,
                ym,
                tm,
                tkm,
                phim,
                Xfe_bulk,
                Xfem,
                marknum,
                1.0e10,
                cfg_perc;
                coords=coords,
            )
            @test iszero(res_cold.max_v_seg)
            @test iszero(res_cold.n_subcycles)
            @test isapprox(sum(Xfe_bulk), initial_sum; atol=1e-14)
        end

        @testset "Machine-Precision Mass Conservation under Percolation" begin
            # Heat interior above eutectic to melt iron
            for m in 1:marknum
                if tm[m] < 3 &&
                    distance(xm[m], ym[m], coords.xcenter, coords.ycenter) <= 50000.0
                    tkm[m] = 1300.0
                    Xfem[m] = Xfe_bulk[m]
                    phim[m] = 0.02 # solid silicate matrix
                end
            end

            initial_sum = sum(Xfe_bulk)
            @test initial_sum > 0.0

            cfg_perc = CoreFormationConfig(
                percolation_active=true,
                settling_active=false,
                phi_crit_perc=0.05,
                phi_residual=0.02,
                cfl_settling=0.5,
            )

            dt_step = 1.0e10 # 317 years
            res = apply_metal_segregation!(
                xm,
                ym,
                tm,
                tkm,
                phim,
                Xfe_bulk,
                Xfem,
                marknum,
                dt_step,
                cfg_perc;
                coords=coords,
                rplanet=50000.0,
            )

            final_sum = sum(Xfe_bulk)
            # Mass conservation must hold to machine precision (10^-12)
            rel_err = abs(final_sum - initial_sum) / initial_sum
            @test rel_err < 1.0e-12
            @test res.max_v_seg > 0.0
            @test res.n_subcycles >= 1

            # Physical segregation check: central markers gain iron, outer markers lose iron
            r_core = 15000.0
            r_outer = 45000.0
            core_fe_sum = 0.0
            core_count = 0
            outer_fe_sum = 0.0
            outer_count = 0
            for m in 1:marknum
                if tm[m] < 3
                    rm = distance(xm[m], ym[m], coords.xcenter, coords.ycenter)
                    if rm <= r_core
                        core_fe_sum += Xfe_bulk[m]
                        core_count += 1
                    elseif rm >= 35000.0 && rm <= r_outer
                        outer_fe_sum += Xfe_bulk[m]
                        outer_count += 1
                    end
                end
            end
            avg_core_fe = core_fe_sum / core_count
            avg_outer_fe = outer_fe_sum / outer_count
            @test avg_core_fe >= 0.20
            @test avg_outer_fe <= 0.20
            @test avg_core_fe > avg_outer_fe
            @test all(x -> 0.0 <= x <= 1.0, Xfe_bulk)
        end

        @testset "CFL Subcycling and Energetics" begin
            # Large timestep triggers multiple subcycles
            dt_large = 1.0e11 # ~3168 years
            Q_seg = zeros(Float64, coords.Ny1, coords.Nx1)

            cfg_perc = CoreFormationConfig(
                percolation_active=true,
                settling_active=false,
                cfl_settling=0.25,
                segregation_heating=true,
            )

            initial_sum = sum(Xfe_bulk)
            res_sub = apply_metal_segregation!(
                xm,
                ym,
                tm,
                tkm,
                phim,
                Xfe_bulk,
                Xfem,
                marknum,
                dt_large,
                cfg_perc;
                coords=coords,
                rplanet=50000.0,
                Q_seg_grid=Q_seg,
            )

            @test res_sub.n_subcycles > 1
            @test res_sub.dt_sub <= dt_large
            @test res_sub.total_dissipation_energy > 0.0
            # Dissipation grid must have positive energy inside planet
            @test maximum(Q_seg) > 0.0

            # Mass conservation across all subcycles
            final_sum = sum(Xfe_bulk)
            @test abs(final_sum - initial_sum) / initial_sum < 1.0e-12
        end

        @testset "Magma Ocean Stokes Settling Mode" begin
            # High silicate melt fraction (Fm = 0.70) triggers settling
            for m in 1:marknum
                if tm[m] < 3 &&
                    distance(xm[m], ym[m], coords.xcenter, coords.ycenter) <= 50000.0
                    phim[m] = 0.70
                    tkm[m] = 1700.0
                    Xfem[m] = Xfe_bulk[m]
                end
            end

            cfg_settle = CoreFormationConfig(
                percolation_active=false,
                settling_active=true,
                F_settle_start=0.40,
                droplet_size_mode=:fixed,
                droplet_diameter_fixed=1.0e-2,
            )

            initial_sum = sum(Xfe_bulk)
            res_settle = apply_metal_segregation!(
                xm,
                ym,
                tm,
                tkm,
                phim,
                Xfe_bulk,
                Xfem,
                marknum,
                1.0e9,
                cfg_settle;
                coords=coords,
                rplanet=50000.0,
            )

            @test res_settle.max_v_seg > 0.0
            @test abs(sum(Xfe_bulk) - initial_sum) / initial_sum < 1.0e-12
            @test all(x -> 0.0 <= x <= 1.0, Xfe_bulk)
        end

        @testset "Multi-Step Sequential Evolution and Packing Bounds" begin
            # Reset markers
            fill!(Xfe_bulk, 0.0)
            for m in 1:marknum
                if tm[m] < 3 &&
                    distance(xm[m], ym[m], coords.xcenter, coords.ycenter) <= 50000.0
                    Xfe_bulk[m] = 0.25
                    tkm[m] = 1400.0
                    Xfem[m] = Xfe_bulk[m]
                    phim[m] = 0.55 # hybrid percolation + settling
                end
            end

            cfg_multi = CoreFormationConfig(
                percolation_active=true,
                settling_active=true,
                phi_pack=0.65,
                cfl_settling=0.5,
            )

            base_sum = sum(Xfe_bulk)
            for step in 1:5
                apply_metal_segregation!(
                    xm,
                    ym,
                    tm,
                    tkm,
                    phim,
                    Xfe_bulk,
                    Xfem,
                    marknum,
                    5.0e9,
                    cfg_multi;
                    coords=coords,
                    rplanet=50000.0,
                )
                curr_sum = sum(Xfe_bulk)
                # Exact conservation at every single timestep
                @test abs(curr_sum - base_sum) / base_sum < 1.0e-12
                # Packing fraction invariant
                @test all(x -> x <= cfg_multi.phi_pack + 1.0e-10, Xfe_bulk)
                @test all(x -> x >= 0.0, Xfe_bulk)
            end
        end

        @testset "Non-Uniform Marker Initialization and Packing Invariant" begin
            # Alternating high/low metal markers in shared cells (Sonnet finding 1 repro)
            fill!(Xfe_bulk, 0.0)
            for m in 1:marknum
                if tm[m] < 3 &&
                    distance(xm[m], ym[m], coords.xcenter, coords.ycenter) <= 50000.0
                    Xfe_bulk[m] = isodd(m) ? 0.60 : 0.02
                    tkm[m] = 1350.0
                    Xfem[m] = Xfe_bulk[m]
                    phim[m] = 0.05
                end
            end

            cfg_hetero = CoreFormationConfig(
                percolation_active=true,
                settling_active=false,
                phi_pack=0.65,
                cfl_settling=0.5,
            )

            initial_hetero_sum = sum(Xfe_bulk)
            for step in 1:20
                apply_metal_segregation!(
                    xm,
                    ym,
                    tm,
                    tkm,
                    phim,
                    Xfe_bulk,
                    Xfem,
                    marknum,
                    5.0e9,
                    cfg_hetero;
                    coords=coords,
                    rplanet=50000.0,
                )
                curr_sum = sum(Xfe_bulk)
                # Exact conservation
                @test abs(curr_sum - initial_hetero_sum) / initial_hetero_sum < 1.0e-12
                # Strict marker-level packing fraction ceiling: NO marker can exceed phi_pack
                @test all(x -> x <= cfg_hetero.phi_pack + 1.0e-10, Xfe_bulk)
                @test all(x -> x >= 0.0, Xfe_bulk)
                @test maximum(Xfe_bulk) <= cfg_hetero.phi_pack + 1.0e-10
            end
        end

        @testset "Input Guard Validation on Metal Packing Exceedance" begin
            # Test that apply_metal_segregation! rejects invalid marker metal fractions
            m_int = findfirst(
                m ->
                    tm[m] < 3 &&
                    distance(xm[m], ym[m], coords.xcenter, coords.ycenter) <= 50000.0,
                1:marknum,
            )
            @test m_int isa Int
            @test 1 <= m_int <= marknum
            fill!(Xfe_bulk, 0.20)
            Xfe_bulk[m_int] = 0.80 # exceeds phi_pack = 0.65
            cfg_guard = CoreFormationConfig(percolation_active=true, phi_pack=0.65)
            @test_throws DomainError apply_metal_segregation!(
                xm,
                ym,
                tm,
                tkm,
                phim,
                Xfe_bulk,
                Xfem,
                marknum,
                1.0e9,
                cfg_guard;
                coords=coords,
                rplanet=50000.0,
            )

            Xfe_bulk[m_int] = -0.05 # negative
            @test_throws DomainError apply_metal_segregation!(
                xm,
                ym,
                tm,
                tkm,
                phim,
                Xfe_bulk,
                Xfem,
                marknum,
                1.0e9,
                cfg_guard;
                coords=coords,
                rplanet=50000.0,
            )
            Xfe_bulk[m_int] = 0.20 # restore
        end

        @testset "Variable Metal Density Mode Segregation Coupling" begin
            # Verify segregation velocity responds to sulfur fraction and EOS law
            # via dynamic liquid metal density in apply_metal_segregation!
            for m in 1:marknum
                if tm[m] < 3 &&
                    distance(xm[m], ym[m], coords.xcenter, coords.ycenter) <= 50000.0
                    tkm[m] = 1350.0
                    Xfem[m] = Xfe_bulk[m]
                    phim[m] = 0.02
                end
            end

            cfg_low_s = CoreFormationConfig(;
                percolation_active=true,
                settling_active=false,
                metal_density_mode=:sanloup2000,
                sulfur_fraction=0.10,
                phi_crit_perc=0.05,
            )
            cfg_high_s = CoreFormationConfig(;
                percolation_active=true,
                settling_active=false,
                metal_density_mode=:sanloup2000,
                sulfur_fraction=0.35,
                phi_crit_perc=0.05,
            )

            res_low_s = apply_metal_segregation!(
                copy(xm),
                copy(ym),
                tm,
                tkm,
                phim,
                copy(Xfe_bulk),
                copy(Xfem),
                marknum,
                1.0e10,
                cfg_low_s;
                coords=coords,
                rplanet=50000.0,
            )
            res_high_s = apply_metal_segregation!(
                copy(xm),
                copy(ym),
                tm,
                tkm,
                phim,
                copy(Xfe_bulk),
                copy(Xfem),
                marknum,
                1.0e10,
                cfg_high_s;
                coords=coords,
                rplanet=50000.0,
            )
            @test res_low_s.max_v_seg > res_high_s.max_v_seg
            @test res_low_s.max_v_seg > 0.0
            @test res_high_s.max_v_seg > 0.0
        end
    end

    @testset "Coupled Core Formation Integration & Benchmarks" begin
        @testset "Benchmark Configuration Loading & Schema Validation" begin
            bench_toml = joinpath(
                @__DIR__, "..", "configs", "core_formation_benchmark.toml"
            )
            @test isfile(bench_toml)

            cfg = load_config(bench_toml)
            @test cfg.coreformation isa CoreFormationConfig
            @test cfg.coreformation.percolation_active == true
            @test cfg.coreformation.settling_active == true
            @test isapprox(cfg.coreformation.sulfur_fraction, 0.31; rtol=1e-12)
            @test cfg.coreformation.metal_density_mode == :sanloup2000
            @test isapprox(cfg.coreformation.rho_metal, 5450.0; rtol=1e-12)
            @test isapprox(cfg.coreformation.rho_metal_solid, 5700.0; rtol=1e-12)
            @test isapprox(cfg.coreformation.eta_metal, 1.0e-2; rtol=1e-12)
            @test isapprox(cfg.coreformation.k_metal, 40.0; rtol=1e-12)
            @test isapprox(cfg.coreformation.rhocp_metal, 4.0e6; rtol=1e-12)
            @test isapprox(cfg.coreformation.Xfe_bulk, 0.20; rtol=1e-12)
            @test isapprox(cfg.coreformation.phi_pack, 0.65; rtol=1e-12)
            @test isapprox(cfg.coreformation.T_eutectic, 1213.0; rtol=1e-12)
            @test isapprox(cfg.coreformation.dT_metal, 50.0; rtol=1e-12)
            @test isapprox(cfg.coreformation.k_metal_ref, 1.0e-9; rtol=1e-12)
            @test isapprox(cfg.coreformation.perm_exponent, 3.0; rtol=1e-12)
            @test isapprox(cfg.coreformation.phi_crit_perc, 0.05; rtol=1e-12)
            @test isapprox(cfg.coreformation.phi_residual, 0.02; rtol=1e-12)
            @test isapprox(cfg.coreformation.phi0, 0.1; rtol=1e-12)
            @test cfg.coreformation.droplet_size_mode == :capillary_mean
            @test isapprox(cfg.coreformation.droplet_diameter_fixed, 5.0e-3; rtol=1e-12)
            @test isapprox(cfg.coreformation.sigma_metal_silicate, 1.0; rtol=1e-12)
            @test isapprox(cfg.coreformation.We_crit, 10.0; rtol=1e-12)
            @test isapprox(cfg.coreformation.hindered_exponent, 4.5; rtol=1e-12)
            @test cfg.coreformation.hadamard_rybczynski == false
            @test isapprox(cfg.coreformation.F_settle_start, 0.40; rtol=1e-12)
            @test isapprox(cfg.coreformation.F_perc_end, 0.50; rtol=1e-12)
            @test cfg.coreformation.segregation_heating == true
            @test isapprox(cfg.coreformation.cfl_settling, 0.5; rtol=1e-12)
            @test cfg.coreformation.max_subcycles == 2000

            # Direct validation check
            @test validate_config(cfg) === nothing
        end

        @testset "Full Multi-Step Coupled Simulation with Iron Segregation" begin
            output_dir = mktempdir()
            try
                bench_toml = joinpath(
                    @__DIR__, "..", "configs", "core_formation_benchmark.toml"
                )
                cfg_base = load_config(bench_toml)
                # Short 3-step run for integration test
                cfg = SimulationConfig(
                    grid=cfg_base.grid,
                    geometry=cfg_base.geometry,
                    time=TimeConfig(
                        n_steps=3,
                        dt_initial=cfg_base.time.dt_initial,
                        dt_longest=cfg_base.time.dt_longest,
                    ),
                    solver=cfg_base.solver,
                    poroelasticity=cfg_base.poroelasticity,
                    thermodynamics=cfg_base.thermodynamics,
                    reaction=cfg_base.reaction,
                    melting=cfg_base.melting,
                    coreformation=cfg_base.coreformation,
                    materials=cfg_base.materials,
                    output=OutputConfig(savematstep=1, output_dir=output_dir),
                )

                Erebus.simulation_loop(cfg; output_path=output_dir)

                for step in 0:3
                    fpath = joinpath(output_dir, "output_$(lpad(step, 5, '0')).jld2")
                    @test isfile(fpath)
                end

                initial_state = load_state(joinpath(output_dir, "output_00000.jld2"))
                final_state = load_state(joinpath(output_dir, "output_00003.jld2"))

                @test haskey(initial_state, "Xfe_bulk")
                @test haskey(initial_state, "Xfem")
                @test haskey(initial_state, "Xfem0")
                @test haskey(final_state, "Xfe_bulk")
                @test haskey(final_state, "Xfem")
                @test haskey(final_state, "Xfem0")

                Xfe_init = initial_state["Xfe_bulk"]
                Xfe_final = final_state["Xfe_bulk"]

                @test length(Xfe_init) == 16384
                @test length(Xfe_final) == length(Xfe_init)

                # Boundedness invariants
                @test all(x -> 0.0 <= x <= cfg.coreformation.phi_pack + 1.0e-10, Xfe_final)
                @test all(x -> 0.0 <= x <= 1.0, final_state["Xfem"])

                # Conservation of initial state metal inventory
                sum_init = sum(Xfe_init)
                sum_final = sum(Xfe_final)
                @test sum_init > 0.0
                @test sum_final > 0.0
                rel_diff = abs(sum_final - sum_init) / sum_init
                @test rel_diff < 1.0e-10

                # Physical core segregation assertions:
                # 1. Interior was preheated to 1350 K, exceeding Fe-FeS eutectic (1213 K)
                @test any(x -> x > 1213.0, final_state["tkm"])
                # 2. Molten metal fraction is non-zero
                @test any(x -> x > 0.0, final_state["Xfem"])
                # 3. Metal segregated inward into core
                @test maximum(final_state["Xfe_bulk"]) > 0.20
            finally
                rm(output_dir, recursive=true, force=true)
            end
        end

        @testset "Coupled Dissipation Heating Thermal Impact" begin
            # Verify that Q_seg_grid is positive and bounded when segregation occurs
            Nx_t = 16
            Ny_t = 16
            coords_t = GridCoordinates(
                GridConfig(; Nx=Nx_t, Ny=Ny_t, xsize=140000.0, ysize=140000.0)
            )
            marknum_t = coords_t.Nxm * coords_t.Nym
            (xm_t, ym_t, tm_t, tkm_t, sxxm_t, sxym_t, etavpm_t, phim_t, phinewm_t, pfm0_t, XWsolidm_t, XWsolidm0_t, Fm_t) = setup_marker_properties(
                marknum_t, coords_t
            )
            (rhototalm_t, rhocptotalm_t, etatotalm_t, hrtotalm_t, ktotalm_t, tkm_rhocptotalm_t, etafluidcur_inv_kphim_t, inv_gggtotalm_t, fricttotalm_t, cohestotalm_t, tenstotalm_t, rhofluidcur_t, alphasolidcur_t, alphafluidcur_t) = setup_marker_properties_helpers(
                marknum_t
            )
            (Xfem_t, Xfem0_t, Xfe_bulk_t) = setup_marker_metal_properties(marknum_t)

            define_markers!(
                xm_t,
                ym_t,
                tm_t,
                phim_t,
                etavpm_t,
                rhototalm_t,
                rhocptotalm_t,
                etatotalm_t,
                hrtotalm_t,
                ktotalm_t,
                tkm_t,
                inv_gggtotalm_t,
                fricttotalm_t,
                cohestotalm_t,
                tenstotalm_t,
                rhofluidcur_t,
                alphasolidcur_t,
                alphafluidcur_t,
                XWsolidm0_t;
                randomized=false,
                coords=coords_t,
                Xfe_bulk=Xfe_bulk_t,
                Xfem=Xfem_t,
                Xfem0=Xfem0_t,
                Xfe_bulk_val=0.20,
                T_eutectic_val=1213.0,
                dT_metal_val=50.0,
            )

            # Warm interior to melt metal
            for m in 1:marknum_t
                if tm_t[m] < 3 &&
                    distance(xm_t[m], ym_t[m], coords_t.xcenter, coords_t.ycenter) <= 50000.0
                    tkm_t[m] = 1350.0
                    Xfem_t[m] = 1.0
                    phim_t[m] = 0.05
                end
            end

            Q_seg_grid = zeros(Float64, coords_t.Ny1, coords_t.Nx1)
            cfg_heat = CoreFormationConfig(
                percolation_active=true,
                settling_active=false,
                segregation_heating=true,
                cfl_settling=0.5,
            )

            res = apply_metal_segregation!(
                xm_t,
                ym_t,
                tm_t,
                tkm_t,
                phim_t,
                Xfe_bulk_t,
                Xfem_t,
                marknum_t,
                1.0e10,
                cfg_heat;
                coords=coords_t,
                rplanet=50000.0,
                Q_seg_grid=Q_seg_grid,
            )

            @test res.total_dissipation_energy > 0.0
            @test maximum(Q_seg_grid) > 0.0
            @test minimum(Q_seg_grid) >= 0.0
            @test all(isfinite, Q_seg_grid)
        end

        @testset "Dynamic Convective Velocity and Rouse Number Coupling" begin
            # Flow field: convective velocity ~ 1e-4 m/s, settling velocity ~ 1e-5 m/s -> R ~ 0.1
            u_conv = 1.0e-4
            v_settle = 1.0e-5
            r_susp = suspension_rouse_number(v_settle, u_conv)
            @test isapprox(r_susp, 0.1; rtol=1e-12)
            @test r_susp < 1.0

            # Quiescent mantle: convective velocity ~ 1e-6 m/s, settling velocity ~ 1e-4 m/s -> R ~ 100
            u_slow = 1.0e-6
            v_fast = 1.0e-4
            r_settle = suspension_rouse_number(v_fast, u_slow)
            @test isapprox(r_settle, 100.0; rtol=1e-12)
            @test r_settle > 1.0
            @test r_settle > r_susp
        end

        @testset "Weber Droplet Size Modes Discrimination" begin
            cfg_cap = CoreFormationConfig(
                percolation_active=false,
                settling_active=true,
                droplet_size_mode=:capillary_mean,
                sigma_metal_silicate=1.0,
                We_crit=10.0,
            )
            cfg_bond = CoreFormationConfig(
                percolation_active=false,
                settling_active=true,
                droplet_size_mode=:bond_mean,
                sigma_metal_silicate=1.0,
                We_crit=10.0,
            )
            cfg_mean = CoreFormationConfig(
                percolation_active=false,
                settling_active=true,
                droplet_size_mode=:weber_mean,
                sigma_metal_silicate=1.0,
                We_crit=10.0,
            )
            cfg_turb = CoreFormationConfig(
                percolation_active=false,
                settling_active=true,
                droplet_size_mode=:weber_turbulent,
                sigma_metal_silicate=1.0,
                We_crit=10.0,
            )

            @test cfg_cap.droplet_size_mode === :capillary_mean
            @test cfg_bond.droplet_size_mode === :bond_mean
            @test cfg_mean.droplet_size_mode === :weber_mean

            # Test gravity-capillary balance d = sqrt(We * sigma / (drho * g))
            drho_val = 3900.0
            g_low = 0.05
            g_high = 0.20
            d_mean_low = sqrt(
                cfg_cap.We_crit * cfg_cap.sigma_metal_silicate / (drho_val * g_low)
            )
            d_mean_high = sqrt(
                cfg_cap.We_crit * cfg_cap.sigma_metal_silicate / (drho_val * g_high)
            )
            @test d_mean_low > d_mean_high
            # 4x higher gravity -> 2x smaller droplet
            @test isapprox(d_mean_low / d_mean_high, 2.0; rtol=1e-10)

            # Capillary and Bond modes compute identical diameter to Weber mean
            d_bond_low = sqrt(
                cfg_bond.We_crit * cfg_bond.sigma_metal_silicate / (drho_val * g_low)
            )
            @test isapprox(d_bond_low, d_mean_low; rtol=1e-12)

            # Weber turbulent relative velocity balance: d = We * sigma / (rho * v^2)
            v_fast = 1.0e-2
            v_slow = 1.0e-3
            d_turb_fast = weber_equilibrium_diameter(3300.0, v_fast, 1.0; We_crit=10.0)
            d_turb_slow = weber_equilibrium_diameter(3300.0, v_slow, 1.0; We_crit=10.0)
            @test d_turb_fast < d_turb_slow
            # 10x faster relative velocity -> 100x smaller droplet
            @test isapprox(d_turb_slow / d_turb_fast, 100.0; rtol=1e-10)
        end

        @testset "Thermal RHS Assembly with Segregation Heating Source" begin
            # Verify assemble_thermal_lse! accepts Q_seg without mutating HR
            coords_local = GridCoordinates(
                GridConfig(; Nx=16, Ny=16, xsize=140000.0, ysize=140000.0)
            )
            Ny1, Nx1 = coords_local.Ny1, coords_local.Nx1
            tk1 = fill(1400.0, Ny1, Nx1)
            RHOCP = fill(3.0e6, Ny1, Nx1)
            KX = fill(3.0, Ny1, Nx1)
            KY = fill(3.0, Ny1, Nx1)
            HR = fill(1.0e-7, Ny1, Nx1)
            HA = zeros(Ny1, Nx1)
            HS = zeros(Ny1, Nx1)
            DHP = zeros(Ny1, Nx1)
            RT1 = zeros(Ny1 * Nx1)
            RT2 = zeros(Ny1 * Nx1)
            Q_seg = fill(5.0e-7, Ny1, Nx1)

            HR_initial = copy(HR)
            LT1 = assemble_thermal_lse!(
                tk1, RHOCP, KX, KY, HR, HA, HS, DHP, RT1, 1.0e9; coords=coords_local
            )
            # Assemble with Q_seg passed directly
            LT2 = assemble_thermal_lse!(
                tk1,
                RHOCP,
                KX,
                KY,
                HR,
                HA,
                HS,
                DHP,
                RT2,
                1.0e9;
                coords=coords_local,
                Q_seg=Q_seg,
            )

            # HR must not be modified in place
            @test HR == HR_initial
            # Interior points must differ exactly by Q_seg
            for j in 2:(Nx1 - 1), i in 2:(Ny1 - 1)
                gk = (j - 1) * Ny1 + i
                @test isapprox(RT2[gk] - RT1[gk], Q_seg[i, j]; atol=1e-14)
            end
        end

        @testset "Segregation Regime Silicate Melt vs Porosity Discrimination" begin
            # Verify that pore water porosity (phim) does NOT trigger Stokes settling when Fm is zero
            coords_discr = GridCoordinates(
                GridConfig(; Nx=9, Ny=9, xsize=70000.0, ysize=70000.0)
            )
            marknum_d = 25
            props_d = setup_marker_properties(marknum_d, coords_discr)
            xm_d = props_d[1]
            ym_d = props_d[2]
            tm_d = props_d[3]
            tkm_d = props_d[4]
            phim_d = props_d[8]
            Fm_d = props_d[13]

            (Xfem_d, Xfem0_d, Xfe_bulk_d) = setup_marker_metal_properties(marknum_d)
            fill!(xm_d, coords_discr.xcenter)
            fill!(ym_d, coords_discr.ycenter + 10000.0) # safely inside planet radius (10 km < 35 km)
            fill!(tm_d, 1) # rock
            fill!(tkm_d, 1250.0) # molten metal (T > 1213 K)
            fill!(Xfe_bulk_d, 0.20)
            fill!(Xfem_d, 0.20)
            fill!(phim_d, 0.35) # high pore water porosity
            fill!(Fm_d, 0.0)    # zero silicate melt

            cfg_settle_only = CoreFormationConfig(
                percolation_active=false, settling_active=true, F_settle_start=0.40
            )

            # Case A: High porosity, but zero silicate melt fraction -> zero settling velocity
            res_nosilicate = apply_metal_segregation!(
                xm_d,
                ym_d,
                tm_d,
                tkm_d,
                phim_d,
                Xfe_bulk_d,
                Xfem_d,
                marknum_d,
                1.0e8,
                cfg_settle_only;
                coords=coords_discr,
                rplanet=35000.0,
                Fm=Fm_d,
            )
            @test iszero(res_nosilicate.max_v_seg)

            # Case B: High silicate melt fraction (Fm = 0.60 > F_settle_start) -> active settling
            fill!(Fm_d, 0.60)
            res_withsilicate = apply_metal_segregation!(
                xm_d,
                ym_d,
                tm_d,
                tkm_d,
                phim_d,
                Xfe_bulk_d,
                Xfem_d,
                marknum_d,
                1.0e8,
                cfg_settle_only;
                coords=coords_discr,
                rplanet=35000.0,
                Fm=Fm_d,
            )
            @test res_withsilicate.max_v_seg > 0.0
        end

        @testset "Adaptive Timestep Segregation CFL Constraint" begin
            coords_dt = GridCoordinates(
                GridConfig(; Nx=17, Ny=17, xsize=70000.0, ysize=70000.0)
            )
            dx = coords_dt.dx
            dy = coords_dt.dy
            vx = zeros(18, 17)
            vy = zeros(17, 18)
            vxf = zeros(17, 18)
            vyf = zeros(18, 17)

            # When max_v_seg > 0, dt must satisfy segregation subcycle ceiling
            v_fast = 1.0e-3 # m/s
            cfl_settle = 0.5
            max_sub = 100
            dt_seg_limit = max_sub * cfl_settle * min(dx, dy) / v_fast

            dt_test = compute_adaptive_timestep(
                vx,
                vy,
                vxf,
                vyf,
                1.0e12, # requested large dt
                0.0;    # aphimax
                coords=coords_dt,
                dt_longest_val=1.0e15,
                max_v_seg=v_fast,
                max_subcycles=max_sub,
                cfl_settling=cfl_settle,
            )
            @test dt_test <= dt_seg_limit + 1.0e-6

            # When max_v_seg == 0, segregation CFL does not constrain dt
            dt_unconstrained = compute_adaptive_timestep(
                vx,
                vy,
                vxf,
                vyf,
                1.0e12,
                0.0;
                coords=coords_dt,
                dt_longest_val=1.0e15,
                max_v_seg=0.0,
            )
            @test dt_unconstrained > dt_seg_limit
        end

        @testset "Thermochemical Iteration Invariance for Metal Segregation" begin
            # Setup a planet domain where markers have molten metal in high silicate melt
            coords_titer = GridCoordinates(
                GridConfig(; Nx=16, Ny=16, xsize=140000.0, ysize=140000.0)
            )
            marknum_t = coords_titer.Nxm * coords_titer.Nym
            props_t = setup_marker_properties(marknum_t, coords_titer)
            xm_t = props_t[1]
            ym_t = props_t[2]
            tm_t = props_t[3]
            tkm_t = props_t[4]
            phim_t = props_t[8]
            Fm_t = props_t[13]
            (Xfem_t, Xfem0_t, Xfe_bulk_t) = setup_marker_metal_properties(marknum_t)

            define_markers!(
                xm_t,
                ym_t,
                tm_t,
                phim_t,
                zeros(marknum_t),
                zeros(marknum_t),
                zeros(marknum_t),
                zeros(marknum_t),
                zeros(marknum_t),
                zeros(marknum_t),
                tkm_t,
                zeros(marknum_t),
                zeros(marknum_t),
                zeros(marknum_t),
                zeros(marknum_t),
                zeros(marknum_t),
                zeros(marknum_t),
                zeros(marknum_t),
                zeros(marknum_t);
                randomized=false,
                coords=coords_titer,
                Xfe_bulk=Xfe_bulk_t,
                Xfem=Xfem_t,
                Xfem0=Xfem0_t,
                Xfe_bulk_val=0.20,
                T_eutectic_val=1213.0,
                dT_metal_val=50.0,
            )

            for m in 1:marknum_t
                if tm_t[m] < 3 &&
                    distance(xm_t[m], ym_t[m], coords_titer.xcenter, coords_titer.ycenter) <=
                   50000.0
                    tkm_t[m] = 1600.0
                    Fm_t[m] = 0.60
                    Xfem_t[m] = Xfe_bulk_t[m]
                end
            end

            cfg_titer = CoreFormationConfig(
                percolation_active=true, settling_active=true, cfl_settling=0.5
            )
            dt_step = 1.0e9

            # Single iteration pass
            Xfe_pass1 = copy(Xfe_bulk_t)
            Xfem_pass1 = copy(Xfem_t)
            res_pass1 = apply_metal_segregation!(
                xm_t,
                ym_t,
                tm_t,
                tkm_t,
                phim_t,
                Xfe_pass1,
                Xfem_pass1,
                marknum_t,
                dt_step,
                cfg_titer;
                coords=coords_titer,
                rplanet=50000.0,
                Fm=Fm_t,
            )

            # Multiple passes without snapshot/restore (compounding error)
            Xfe_unrestored = copy(Xfe_bulk_t)
            Xfem_unrestored = copy(Xfem_t)
            for _ in 1:3
                apply_metal_segregation!(
                    xm_t,
                    ym_t,
                    tm_t,
                    tkm_t,
                    phim_t,
                    Xfe_unrestored,
                    Xfem_unrestored,
                    marknum_t,
                    dt_step,
                    cfg_titer;
                    coords=coords_titer,
                    rplanet=50000.0,
                    Fm=Fm_t,
                )
            end

            # Multiple passes WITH snapshot/restore (simulation_loop invariant)
            Xfe_restored = copy(Xfe_bulk_t)
            Xfem_restored = copy(Xfem_t)
            Xfe_snap = copy(Xfe_bulk_t)
            Xfem_snap = copy(Xfem_t)
            for _ in 1:3
                copyto!(Xfe_restored, Xfe_snap)
                copyto!(Xfem_restored, Xfem_snap)
                apply_metal_segregation!(
                    xm_t,
                    ym_t,
                    tm_t,
                    tkm_t,
                    phim_t,
                    Xfe_restored,
                    Xfem_restored,
                    marknum_t,
                    dt_step,
                    cfg_titer;
                    coords=coords_titer,
                    rplanet=50000.0,
                    Fm=Fm_t,
                )
            end

            @test res_pass1.max_v_seg > 0.0
            # Unrestored passes compound and artificially over-transport metal mass
            @test sum(abs.(Xfe_unrestored .- Xfe_bulk_t)) >
                sum(abs.(Xfe_pass1 .- Xfe_bulk_t))
            # Restored snapshot matches single pass to exact floating-point precision
            @test isapprox(Xfe_restored, Xfe_pass1; atol=1e-14)
            @test isapprox(sum(Xfe_restored), sum(Xfe_bulk_t); rtol=1e-12)
            # Verify Xfem is updated consistently with newly segregated bulk metal
            @test isapprox(
                Xfe_pass1 .* compute_metal_melt_fraction.(
                    tkm_t; T_eutectic=cfg_titer.T_eutectic, dT_metal=cfg_titer.dT_metal
                ),
                Xfem_pass1;
                atol=1e-14,
            )
        end

        @testset "Snapshot Buffer Resizing Under Marker Replenishment" begin
            n_initial = 100
            n_grown = 120
            Xfe_bulk_m = fill(0.20, n_initial)
            Xfem_m = fill(0.10, n_initial)
            Xfe_snap = zeros(Float64, n_initial)
            Xfem_snap = zeros(Float64, n_initial)

            # Initial snapshot
            copyto!(Xfe_snap, Xfe_bulk_m)
            copyto!(Xfem_snap, Xfem_m)

            # Replenishment grows arrays
            append!(Xfe_bulk_m, fill(0.15, n_grown - n_initial))
            append!(Xfem_m, fill(0.05, n_grown - n_initial))

            # Resizing logic from simulation_loop step start
            if length(Xfe_snap) != length(Xfe_bulk_m)
                resize!(Xfe_snap, length(Xfe_bulk_m))
            end
            if length(Xfem_snap) != length(Xfem_m)
                resize!(Xfem_snap, length(Xfem_m))
            end
            copyto!(Xfe_snap, Xfe_bulk_m)
            copyto!(Xfem_snap, Xfem_m)
            @test length(Xfe_snap) == n_grown
            @test length(Xfem_snap) == n_grown
            @test isapprox(Xfe_snap[end], 0.15; atol=1e-14)

            # Restoration inside titer loop
            if length(Xfe_bulk_m) != length(Xfe_snap)
                resize!(Xfe_bulk_m, length(Xfe_snap))
            end
            copyto!(Xfe_bulk_m, Xfe_snap)
            @test length(Xfe_bulk_m) == n_grown
        end

        @testset "60Fe Radiogenic Heating with Inactive Core Formation" begin
            # Verify that when hr_fe is active but core formation is inactive,
            # Xfe_bulk is allocated and markers receive 60Fe decay power
            output_dir = mktempdir()
            try
                cfg_fe_only = SimulationConfig(;
                    grid=GridConfig(; Nx=8, Ny=8, xsize=70000.0, ysize=70000.0),
                    time=TimeConfig(; n_steps=2, dt_initial=1.0e8),
                    thermodynamics=ThermalConfig(; hr_al=false, hr_fe=true),
                    coreformation=CoreFormationConfig(;
                        percolation_active=false, settling_active=false
                    ),
                    output=OutputConfig(; output_dir=output_dir, savematstep=1),
                )
                @test run_simulation(cfg_fe_only) === nothing
                final_state = load_state(joinpath(output_dir, "output_00002.jld2"))
                @test haskey(final_state, "Xfe_bulk")
                @test maximum(final_state["tk2"]) > 150.0
                @test final_state["timestep"] == 2
            finally
                rm(output_dir; recursive=true, force=true)
            end
        end
    end
end

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
using Logging

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

        # 5. Exact Darcy law independence above phi_crit_perc and trap fraction sensitivity
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
        km_exact = metal_permeability(phi_m; phi_crit_perc=0.05)
        v_darcy_exact =
            (km_exact / (phi_m * CoreFormationConfig().eta_metal)) * drho * g_acc
        @test isapprox(v_res_low, v_darcy_exact; rtol=1e-12)
        @test isapprox(v_res_high, v_darcy_exact; rtol=1e-12)

        # Trap fraction controls irreducible retention and vanishes across settling transition
        @test isapprox(
            effective_metal_trap_fraction(0.0, 0.04; F_settle_start=0.40, F_perc_end=0.50),
            0.04;
            rtol=1e-12,
        )
        @test iszero(
            effective_metal_trap_fraction(0.60, 0.04; F_settle_start=0.40, F_perc_end=0.50)
        )

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

        r_floor_seg = 0.5e-12 * min(coords.dx, coords.dy)
        weighted_metal_sum(Xfe_arr) = sum(
            Xfe_arr[m] *
            2.0 *
            max(distance(xm[m], ym[m], coords.xcenter, coords.ycenter), r_floor_seg) for
            m in 1:marknum if tm[m] < 3 &&
                distance(xm[m], ym[m], coords.xcenter, coords.ycenter) <= 50000.0
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
            initial_sum = weighted_metal_sum(Xfe_bulk)

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
            @test isapprox(weighted_metal_sum(Xfe_bulk), initial_sum; atol=1e-14)
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

            initial_sum = weighted_metal_sum(Xfe_bulk)
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

            final_sum = weighted_metal_sum(Xfe_bulk)
            # 3D shell volume-weighted mass conservation holds to machine precision
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

            initial_sum = weighted_metal_sum(Xfe_bulk)
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
            final_sum = weighted_metal_sum(Xfe_bulk)
            @test abs(final_sum - initial_sum) / initial_sum < 1.0e-12
        end

        @testset "Magma Ocean Stokes Settling Mode" begin
            # High silicate melt fraction (Fm = 0.70) triggers settling
            Fm_high = zeros(Float64, marknum)
            for m in 1:marknum
                if tm[m] < 3 &&
                    distance(xm[m], ym[m], coords.xcenter, coords.ycenter) <= 50000.0
                    phim[m] = 0.05
                    Fm_high[m] = 0.70
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

            initial_sum = weighted_metal_sum(Xfe_bulk)
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
                Fm=Fm_high,
            )

            @test 0.0 < res_settle.max_v_seg < 1.0
            @test abs(weighted_metal_sum(Xfe_bulk) - initial_sum) / initial_sum < 1.0e-12
            @test all(x -> 0.0 <= x <= 1.0, Xfe_bulk)
        end

        @testset "Multi-Step Sequential Evolution and Packing Bounds" begin
            # Reset markers
            fill!(Xfe_bulk, 0.0)
            Fm_multi = zeros(Float64, marknum)
            for m in 1:marknum
                if tm[m] < 3 &&
                    distance(xm[m], ym[m], coords.xcenter, coords.ycenter) <= 50000.0
                    Xfe_bulk[m] = 0.25
                    tkm[m] = 1620.0
                    Xfem[m] = Xfe_bulk[m]
                    phim[m] = 0.05
                    Fm_multi[m] = 0.55 # hybrid percolation + settling
                end
            end

            cfg_multi = CoreFormationConfig(
                percolation_active=true,
                settling_active=true,
                phi_pack=0.65,
                cfl_settling=0.5,
            )

            base_sum = weighted_metal_sum(Xfe_bulk)
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
                    Fm=Fm_multi,
                )
                curr_sum = weighted_metal_sum(Xfe_bulk)
                # Exact conservation at every single timestep
                @test abs(curr_sum - base_sum) / base_sum < 1.0e-12
                # Packing fraction invariant
                @test all(x -> x <= cfg_multi.phi_pack + 1.0e-10, Xfe_bulk)
                @test all(x -> x >= 0.0, Xfe_bulk)
            end
        end

        @testset "Non-Uniform Marker Initialization and Packing Invariant" begin
            # Alternating high/low metal markers in shared cells
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

            initial_hetero_sum = weighted_metal_sum(Xfe_bulk)
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
                curr_sum = weighted_metal_sum(Xfe_bulk)
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

        @testset "3D Spherical Volume Conservation vs 2D Slab Mode" begin
            fill!(Xfe_bulk, 0.0)
            w3d_vec = zeros(Float64, marknum)
            for m in 1:marknum
                rm = distance(xm[m], ym[m], coords.xcenter, coords.ycenter)
                w3d_vec[m] = 2.0 * max(rm, r_floor_seg)
                if tm[m] < 3 && rm <= 50000.0
                    Xfe_bulk[m] = 0.20
                    tkm[m] = 1400.0
                    Xfem[m] = 0.20
                    phim[m] = 0.02
                end
            end
            cfg_mode = CoreFormationConfig(;
                percolation_active=true, settling_active=false, phi_pack=0.90
            )
            Xfe_3d = copy(Xfe_bulk)
            Xfem_3d = copy(Xfem)
            init_3d = weighted_metal_sum(Xfe_3d)
            init_2d_unweighted = sum(Xfe_3d)
            apply_metal_segregation!(
                xm,
                ym,
                tm,
                tkm,
                phim,
                Xfe_3d,
                Xfem_3d,
                marknum,
                2.0e10,
                cfg_mode;
                coords=coords,
                rplanet=50000.0,
                spherical_metric=true,
                w3d_m=w3d_vec,
            )
            @test isapprox(weighted_metal_sum(Xfe_3d), init_3d; rtol=1.0e-12)
            # Inward radial transport into smaller spherical shells increases 2D cross-sectional area fraction
            @test sum(Xfe_3d) > init_2d_unweighted

            Xfe_2d = copy(Xfe_bulk)
            Xfem_2d = copy(Xfem)
            apply_metal_segregation!(
                xm,
                ym,
                tm,
                tkm,
                phim,
                Xfe_2d,
                Xfem_2d,
                marknum,
                2.0e10,
                cfg_mode;
                coords=coords,
                rplanet=50000.0,
                spherical_metric=false,
            )
            @test isapprox(sum(Xfe_2d), init_2d_unweighted; rtol=1.0e-12)
        end

        @testset "Staggered Gravity Symmetry and P-Node Dissipation Indexing" begin
            coords_sym = GridCoordinates(5, 5; xsize=40000.0, ysize=40000.0)
            xc_sym = 20000.0
            yc_sym = 20000.0
            dx_sym = coords_sym.dx
            dy_sym = coords_sym.dy
            n_sym = 16 * 4
            xm_sym = zeros(Float64, n_sym)
            ym_sym = zeros(Float64, n_sym)
            tm_sym = ones(Int, n_sym)
            tkm_sym = fill(1400.0, n_sym)
            phim_sym = fill(0.02, n_sym)
            Xfe_sym = fill(0.20, n_sym)
            Xfem_sym = fill(0.20, n_sym)
            idx_s = 1
            for j in 1:4, i in 1:4
                for kx in (0.25, 0.75), ky in (0.25, 0.75)
                    xm_sym[idx_s] = (j - 1 + kx) * dx_sym
                    ym_sym[idx_s] = (i - 1 + ky) * dy_sym
                    idx_s += 1
                end
            end
            gx_stag = zeros(Float64, 6, 6)
            gy_stag = zeros(Float64, 6, 6)
            for j in 1:5, i in 1:6
                rx = coords_sym.xvx[j] - xc_sym
                ry = coords_sym.yvx[i] - yc_sym
                rm = hypot(rx, ry)
                gx_stag[i, j] = rm > 0.0 ? -1.0 * (rx / rm) : 0.0
            end
            for j in 1:6, i in 1:5
                rx = coords_sym.xvy[j] - xc_sym
                ry = coords_sym.yvy[i] - yc_sym
                rm = hypot(rx, ry)
                gy_stag[i, j] = rm > 0.0 ? -1.0 * (ry / rm) : 0.0
            end
            Q_seg_sym = zeros(Float64, 6, 6)
            cfg_sym = CoreFormationConfig(;
                percolation_active=true, settling_active=false, segregation_heating=true
            )
            apply_metal_segregation!(
                xm_sym,
                ym_sym,
                tm_sym,
                tkm_sym,
                phim_sym,
                Xfe_sym,
                Xfem_sym,
                n_sym,
                5.0e9,
                cfg_sym;
                coords=coords_sym,
                xcenter=xc_sym,
                ycenter=yc_sym,
                rplanet=25000.0,
                gx=gx_stag,
                gy=gy_stag,
                Q_seg_grid=Q_seg_sym,
            )
            @test all(iszero, Q_seg_sym[1, :])
            @test all(iszero, Q_seg_sym[:, 1])
            @test all(iszero, Q_seg_sym[6, :])
            @test all(iszero, Q_seg_sym[:, 6])
            @test isapprox(Q_seg_sym[2, 2], Q_seg_sym[5, 5]; rtol=1.0e-12)
            @test isapprox(Q_seg_sym[2, 5], Q_seg_sym[5, 2]; rtol=1.0e-12)
            @test isapprox(Q_seg_sym[3, 3], Q_seg_sym[4, 4]; rtol=1.0e-12)
            @test isapprox(Q_seg_sym[3, 4], Q_seg_sym[4, 3]; rtol=1.0e-12)
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
                # The transfer conserves the sum at fixed marker lengths; marker motion
                # changes these out-of-plane lengths, so the current sum is close only.
                w_init = initial_state["w3d_m"]
                sum_init = sum(Xfe_init .* w_init)
                @test sum_init > 0.0
                @test abs(sum(Xfe_final .* w_init) - sum_init) / sum_init < 1.0e-6
                sum_final = sum(Xfe_final .* final_state["w3d_m"])
                @test abs(sum_final - sum_init) / sum_init < 1.0e-4
                @test maximum(abs, Xfe_final .- Xfe_init) > 0.1

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
                    Xfem_t[m] = Xfe_bulk_t[m]
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

        @testset "Segregation heat follows the moved metal" begin
            function metal_body(; X0, rplanet=50000.0)
                coords_e = GridCoordinates(
                    GridConfig(; Nx=33, Ny=49, xsize=140000.0, ysize=140000.0)
                )
                nm = coords_e.Nxm * coords_e.Nym
                xm_e, ym_e = zeros(nm), zeros(nm)
                k = 0
                for j in 1:(coords_e.Nxm), i in 1:(coords_e.Nym)
                    k += 1
                    xm_e[k] = (j - 0.5) * coords_e.xsize / coords_e.Nxm
                    ym_e[k] = (i - 0.5) * coords_e.ysize / coords_e.Nym
                end
                r_e = hypot.(xm_e .- coords_e.xcenter, ym_e .- coords_e.ycenter)
                tm_e = [r <= rplanet ? 1 : 3 for r in r_e]
                Xfe_e = [r <= rplanet ? X0 : 0.0 for r in r_e]
                return (;
                    coords=coords_e,
                    nm,
                    xm=xm_e,
                    ym=ym_e,
                    r=r_e,
                    tm=tm_e,
                    tk=[r <= rplanet ? 1350.0 : 200.0 for r in r_e],
                    phi=[r <= rplanet ? 0.05 : 1.0 for r in r_e],
                    Xfe=Xfe_e,
                    Xfem=[r <= rplanet ? X0 : 0.0 for r in r_e],
                    Q=zeros(coords_e.Ny1, coords_e.Nx1),
                )
            end
            cfg_e = CoreFormationConfig(;
                percolation_active=true, settling_active=false, segregation_heating=true
            )
            run_body!(b, cfg, dt; spherical) = apply_metal_segregation!(
                b.xm,
                b.ym,
                b.tm,
                b.tk,
                b.phi,
                b.Xfe,
                b.Xfem,
                b.nm,
                dt,
                cfg;
                coords=b.coords,
                rplanet=50000.0,
                Q_seg_grid=b.Q,
                spherical_metric=spherical,
            )

            # A packed body cannot move metal, so it releases no heat.
            for spherical in (false, true)
                packed = metal_body(; X0=cfg_e.phi_pack)
                res_packed = run_body!(packed, cfg_e, 1.0e9; spherical=spherical)
                free = metal_body(; X0=0.20)
                res_free = run_body!(free, cfg_e, 1.0e9; spherical=spherical)
                @test iszero(res_packed.max_v_seg)
                @test iszero(res_packed.max_v_char)
                @test res_packed.n_subcycles == 0
                @test !iszero(res_free.max_v_seg)
                @test !iszero(maximum(free.Q))
                @test iszero(maximum(packed.Q))
            end

            # Over a long time the heat equals the potential energy released.
            for spherical in (false, true)
                body = metal_body(; X0=0.20)
                weight = spherical ? 2.0 .* body.r : ones(body.nm)
                Am = body.coords.dxm * body.coords.dym
                drho = cfg_e.rho_metal - 3000.0
                potential(X) = sum(
                    X[m] * weight[m] * Am * drho * 0.1 * body.r[m]^2 / (2.0 * 50000.0)
                    for m in 1:(body.nm) if body.tm[m] < 3
                )
                pe_start = potential(body.Xfe)
                heat = 0.0
                dt_e = 1.0e10
                xp = body.coords.xp
                yp = body.coords.yp
                for _ in 1:200
                    fill!(body.Q, 0.0)
                    run_body!(body, cfg_e, dt_e; spherical=spherical)
                    for j in 2:(body.coords.Nx1 - 1), i in 2:(body.coords.Ny1 - 1)
                        r_node = hypot(
                            xp[j] - body.coords.xcenter, yp[i] - body.coords.ycenter
                        )
                        w_node = spherical ? 2.0 * r_node : 1.0
                        heat +=
                            body.Q[i, j] * w_node * body.coords.dx * body.coords.dy * dt_e
                    end
                end
                release = pe_start - potential(body.Xfe)
                @test !iszero(release)
                @test 0.9 < heat / release < 1.15
            end
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
            ym_d[1:12] .= coords_discr.ycenter + 15000.0
            ym_d[13:marknum_d] .= coords_discr.ycenter + 15000.0 - coords_discr.dy
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
            @test 0.0 < res_withsilicate.max_v_seg < 1.0
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
            r_floor_t = 0.5e-12 * min(coords_titer.dx, coords_titer.dy)
            weighted_sum_t(arr) = sum(
                arr[m] *
                2.0 *
                max(
                    distance(xm_t[m], ym_t[m], coords_titer.xcenter, coords_titer.ycenter),
                    r_floor_t,
                ) for m in 1:marknum_t if tm_t[m] < 3 &&
                    distance(
                        xm_t[m], ym_t[m], coords_titer.xcenter, coords_titer.ycenter
                    ) <= 50000.0
            )
            @test isapprox(Xfe_restored, Xfe_pass1; atol=1e-14)
            @test isapprox(
                weighted_sum_t(Xfe_restored), weighted_sum_t(Xfe_bulk_t); rtol=1e-12
            )
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
                res_fe = run_simulation(cfg_fe_only)
                @test res_fe isa SimulationState
                final_state = load_state(joinpath(output_dir, "output_00002.jld2"))
                @test haskey(final_state, "Xfe_bulk")
                @test maximum(final_state["tk2"]) > 150.0
                @test final_state["timestep"] == 2
            finally
                rm(output_dir; recursive=true, force=true)
            end
        end

        @testset "Silicate Suspension Viscosity and Melt-Blended Density" begin
            eta_s = 1.0e19
            alpha_e = 28.0
            phi_c = 0.40
            eta_m = 10.0

            # Exponential melt weakening below phi_crit
            for F in (0.0, 0.20, 0.40)
                expected = eta_s * exp(-alpha_e * F)
                got = compute_silicate_suspension_viscosity(
                    eta_s, F; alpha_eta=alpha_e, phi_crit=phi_c, eta_melt=eta_m
                )
                @test isapprox(got, expected; rtol=1e-12)
            end

            # Continuity across phi_crit and monotonic smoothing across transition band
            v_minus = compute_silicate_suspension_viscosity(
                eta_s, phi_c - 1.0e-7; alpha_eta=alpha_e, phi_crit=phi_c, eta_melt=eta_m
            )
            v_at = compute_silicate_suspension_viscosity(
                eta_s, phi_c; alpha_eta=alpha_e, phi_crit=phi_c, eta_melt=eta_m
            )
            v_plus = compute_silicate_suspension_viscosity(
                eta_s, phi_c + 1.0e-7; alpha_eta=alpha_e, phi_crit=phi_c, eta_melt=eta_m
            )
            v_mid = compute_silicate_suspension_viscosity(
                eta_s, 0.45; alpha_eta=alpha_e, phi_crit=phi_c, eta_melt=eta_m
            )
            v_narrow = compute_silicate_suspension_viscosity(
                eta_s,
                0.45;
                alpha_eta=alpha_e,
                phi_crit=phi_c,
                eta_melt=eta_m,
                F_blend_width=0.05,
            )
            @test abs(log10(v_plus) - log10(v_minus)) < 1.0e-5
            @test isapprox(v_at, eta_s * exp(-alpha_e * phi_c); rtol=1e-12)
            @test v_at > v_mid > compute_silicate_suspension_viscosity(eta_s, 0.50)
            @test isapprox(
                v_narrow, eta_m * (1.0 - (1.0 - 0.45) / (1.0 - phi_c))^(-2.5); rtol=1e-10
            )

            # Soft-solid min(log_roscoe, log_crit) cap and etamin/etamax bounds
            v_soft_at = compute_silicate_suspension_viscosity(
                1.0e3, phi_c; alpha_eta=alpha_e, phi_crit=phi_c, eta_melt=eta_m
            )
            v_soft_mid = compute_silicate_suspension_viscosity(
                1.0e3, 0.45; alpha_eta=alpha_e, phi_crit=phi_c, eta_melt=eta_m
            )
            @test isapprox(v_soft_at, v_soft_mid; rtol=1e-12)
            @test isapprox(
                compute_silicate_suspension_viscosity(eta_s, 1.0; etamin=500.0),
                500.0;
                rtol=1e-12,
            )
            @test isapprox(
                compute_silicate_suspension_viscosity(eta_s, 0.0; etamax=1.0e15),
                1.0e15;
                rtol=1e-12,
            )

            # Exact Roscoe (1952) / Solomatov (2015) branch above blend window
            for F in (0.50, 0.60, 0.80, 1.0)
                roscoe_val = eta_m * (1.0 - (1.0 - F) / (1.0 - phi_c))^(-2.5)
                got = compute_silicate_suspension_viscosity(
                    eta_s, F; alpha_eta=alpha_e, phi_crit=phi_c, eta_melt=eta_m
                )
                @test isapprox(got, roscoe_val; rtol=1e-10)
            end

            # Input domain contracts
            @test_throws DomainError compute_silicate_suspension_viscosity(eta_s, -0.01)
            @test_throws DomainError compute_silicate_suspension_viscosity(eta_s, 1.01)
            @test_throws DomainError compute_silicate_suspension_viscosity(eta_s, NaN)
            @test_throws DomainError compute_silicate_suspension_viscosity(0.0, 0.2)
            @test_throws DomainError compute_silicate_suspension_viscosity(
                eta_s, 0.2; phi_crit=1.0
            )
            @test_throws DomainError compute_silicate_suspension_viscosity(
                eta_s, 0.2; alpha_eta=-1.0
            )
            @test_throws DomainError compute_silicate_suspension_viscosity(
                eta_s, 0.2; eta_melt=0.0
            )
            @test_throws DomainError compute_silicate_suspension_viscosity(
                eta_s, 0.2; etamin=0.0
            )
            @test_throws DomainError compute_silicate_suspension_viscosity(
                eta_s, 0.2; etamin=20.0, etamax=10.0
            )
            @test_throws DomainError compute_silicate_suspension_viscosity(
                eta_s, 0.2; F_blend_width=0.0
            )

            # Exact parity between public metal_segregation_velocity and internal cell regime helper
            cfg_parity = CoreFormationConfig(;
                percolation_active=true,
                settling_active=true,
                droplet_size_mode=:fixed,
                droplet_diameter_fixed=5.0e-3,
            )
            for F_chk in (0.0, 0.45, 0.80)
                eta_chk = compute_silicate_suspension_viscosity(eta_s, F_chk)
                v0_chk = stokes_settling_velocity(
                    cfg_parity.droplet_diameter_fixed / 2.0,
                    2200.0,
                    0.08,
                    eta_chk;
                    hadamard_rybczynski=cfg_parity.hadamard_rybczynski,
                    eta_metal=cfg_parity.eta_metal,
                )
                v_pub = metal_segregation_velocity(
                    0.20,
                    F_chk,
                    2200.0,
                    0.08,
                    eta_chk;
                    percolation_active=true,
                    settling_active=true,
                    k_metal_ref=cfg_parity.k_metal_ref,
                    eta_metal=cfg_parity.eta_metal,
                    phi_crit_perc=cfg_parity.phi_crit_perc,
                    phi_residual=cfg_parity.phi_residual,
                    phi0=cfg_parity.phi0,
                    perm_exponent=cfg_parity.perm_exponent,
                    r_drop=cfg_parity.droplet_diameter_fixed / 2.0,
                    hindered_exponent=cfg_parity.hindered_exponent,
                    phi_pack=cfg_parity.phi_pack,
                    hadamard_rybczynski=cfg_parity.hadamard_rybczynski,
                    F_settle_start=cfg_parity.F_settle_start,
                    F_perc_end=cfg_parity.F_perc_end,
                )
                sp_int = Erebus._metal_cell_regime_speeds(
                    0.20, F_chk, 2200.0, 0.08, v0_chk, cfg_parity.phi_residual, cfg_parity
                )
                @test isapprox(sp_int.v_seg, v_pub; rtol=1e-12)
            end
            @test_throws DomainError metal_segregation_velocity(
                0.10, -0.1, 2000.0, 0.1, 1.0e18
            )

            # Lighter silicate melt increases metal-silicate density contrast and settling speed proportionally
            coords_rho = GridCoordinates(
                GridConfig(; Nx=7, Ny=7, xsize=70000.0, ysize=70000.0)
            )
            n_rho = coords_rho.start_marknum
            xm_r, ym_r, tm_r = Float64[], Float64[], Int[]
            tkm_r, phim_r, Xfe_r, Xfem_r, Fm_r = Float64[],
            Float64[], Float64[], Float64[],
            Float64[]
            for j in 1:(coords_rho.Nxm), i in 1:(coords_rho.Nym)
                x_v = (j - 0.5) * coords_rho.dxm
                y_v = (i - 0.5) * coords_rho.dym
                in_p = hypot(x_v - coords_rho.xcenter, y_v - coords_rho.ycenter) <= 30000.0
                push!(xm_r, x_v)
                push!(ym_r, y_v)
                push!(tm_r, in_p ? 1 : 3)
                push!(tkm_r, in_p ? 1750.0 : 200.0)
                push!(phim_r, in_p ? 0.05 : 1.0)
                push!(Xfe_r, in_p ? 0.20 : 0.0)
                push!(Xfem_r, in_p ? 0.20 : 0.0)
                push!(Fm_r, in_p ? 0.80 : 0.0)
            end
            cfg_rho = CoreFormationConfig(;
                percolation_active=false,
                settling_active=true,
                droplet_size_mode=:fixed,
                droplet_diameter_fixed=5.0e-3,
                metal_density_mode=:constant,
                rho_metal=5450.0,
            )
            res_dense_melt = apply_metal_segregation!(
                copy(xm_r),
                copy(ym_r),
                copy(tm_r),
                copy(tkm_r),
                copy(phim_r),
                copy(Xfe_r),
                copy(Xfem_r),
                n_rho,
                1.0e5,
                cfg_rho;
                coords=coords_rho,
                rplanet=30000.0,
                Fm=Fm_r,
                rho_silicate=3300.0,
                rho_silicate_melt=2700.0,
            )
            res_heavy_melt = apply_metal_segregation!(
                copy(xm_r),
                copy(ym_r),
                copy(tm_r),
                copy(tkm_r),
                copy(phim_r),
                copy(Xfe_r),
                copy(Xfem_r),
                n_rho,
                1.0e5,
                cfg_rho;
                coords=coords_rho,
                rplanet=30000.0,
                Fm=Fm_r,
                rho_silicate=3300.0,
                rho_silicate_melt=3600.0,
            )
            drho_dense_exp = 5450.0 - (0.20 * 3300.0 + 0.80 * 2700.0)
            drho_heavy_exp = 5450.0 - (0.20 * 3300.0 + 0.80 * 3600.0)
            @test isfinite(res_dense_melt.max_v_seg) &&
                1.0e-12 < res_dense_melt.max_v_seg < 1.0
            @test isapprox(
                res_dense_melt.max_v_seg / res_heavy_melt.max_v_seg,
                drho_dense_exp / drho_heavy_exp;
                rtol=0.05,
            )

            # Explicit ETA keyword grid viscosity overrides scalar eta_silicate
            cfg_eta_kw = CoreFormationConfig(;
                percolation_active=false,
                settling_active=true,
                droplet_size_mode=:fixed,
                droplet_diameter_fixed=5.0e-3,
                F_settle_start=0.0,
                F_perc_end=0.2,
            )
            res_eta_low = apply_metal_segregation!(
                copy(xm_r),
                copy(ym_r),
                copy(tm_r),
                copy(tkm_r),
                copy(phim_r),
                copy(Xfe_r),
                copy(Xfem_r),
                n_rho,
                1.0e5,
                cfg_eta_kw;
                coords=coords_rho,
                rplanet=30000.0,
                Fm=fill(0.10, n_rho),
                ETA=fill(1.0e6, coords_rho.Ny, coords_rho.Nx),
            )
            res_eta_high = apply_metal_segregation!(
                copy(xm_r),
                copy(ym_r),
                copy(tm_r),
                copy(tkm_r),
                copy(phim_r),
                copy(Xfe_r),
                copy(Xfem_r),
                n_rho,
                1.0e5,
                cfg_eta_kw;
                coords=coords_rho,
                rplanet=30000.0,
                Fm=fill(0.10, n_rho),
                ETA=fill(1.0e10, coords_rho.Ny, coords_rho.Nx),
            )
            @test isfinite(res_eta_low.max_v_seg) && 1.0e-15 < res_eta_low.max_v_seg < 1.0
            @test res_eta_low.max_v_seg > 100.0 * res_eta_high.max_v_seg
        end

        @testset "Effective Metal Trap Fraction and Flux Retention" begin
            phi_res = 0.04
            @test isapprox(
                effective_metal_trap_fraction(
                    0.0, phi_res; F_settle_start=0.40, F_perc_end=0.50
                ),
                phi_res;
                rtol=1e-12,
            )
            @test isapprox(
                effective_metal_trap_fraction(
                    0.45, phi_res; F_settle_start=0.40, F_perc_end=0.50
                ),
                0.5 * phi_res;
                rtol=1e-12,
            )
            @test iszero(
                effective_metal_trap_fraction(
                    0.50, phi_res; F_settle_start=0.40, F_perc_end=0.50
                ),
            )
            @test isapprox(
                effective_metal_trap_fraction(
                    0.80, phi_res; settling_active=false, percolation_active=true
                ),
                phi_res;
                rtol=1e-12,
            )
            @test iszero(
                effective_metal_trap_fraction(
                    0.80, phi_res; settling_active=false, percolation_active=false
                ),
            )
            @test_throws DomainError effective_metal_trap_fraction(-0.1, phi_res)
            @test_throws DomainError effective_metal_trap_fraction(0.2, 1.1)
            @test_throws DomainError effective_metal_trap_fraction(
                0.2, phi_res; F_settle_start=0.60, F_perc_end=0.50
            )

            # Column flux retention with validated config (phi_residual = 0.04 <= phi_crit_perc = 0.05)
            coords_tr = GridCoordinates(
                GridConfig(; Nx=15, Ny=15, xsize=70000.0, ysize=70000.0)
            )
            n_tr = coords_tr.start_marknum
            xm_t2, ym_t2, tm_t2, phim_t2 = Float64[], Float64[], Int[], Float64[]
            Xfe_solid = Float64[]
            for j in 1:(coords_tr.Nxm), i in 1:(coords_tr.Nym)
                x_v = (j - 0.5) * coords_tr.dxm
                y_v = (i - 0.5) * coords_tr.dym
                in_p = hypot(x_v - coords_tr.xcenter, y_v - coords_tr.ycenter) <= 30000.0
                push!(xm_t2, x_v)
                push!(ym_t2, y_v)
                push!(tm_t2, in_p ? 1 : 3)
                push!(phim_t2, in_p ? 0.05 : 1.0)
                push!(Xfe_solid, in_p ? 0.08 : 0.0)
            end
            Xfe_mush = copy(Xfe_solid)
            Xfe_liq = copy(Xfe_solid)
            cfg_tr = CoreFormationConfig(;
                percolation_active=true,
                settling_active=true,
                droplet_size_mode=:fixed,
                droplet_diameter_fixed=5.0e-3,
                phi_crit_perc=0.05,
                phi_residual=0.04,
                k_metal_ref=1.0e-8,
                F_settle_start=0.40,
                F_perc_end=0.50,
            )
            @test Erebus.validate_config(
                SimulationConfig(;
                    melting=MeltingConfig(; active=true), coreformation=cfg_tr
                ),
            ) === nothing
            res_cold_trap = apply_metal_segregation!(
                copy(xm_t2),
                copy(ym_t2),
                copy(tm_t2),
                fill(1350.0, n_tr),
                copy(phim_t2),
                Xfe_solid,
                copy(Xfe_solid),
                n_tr,
                5.0e9,
                cfg_tr;
                coords=coords_tr,
                rplanet=30000.0,
                Fm=fill(0.0, n_tr),
                eta_silicate=1.0e6,
            )
            res_mush_trap = apply_metal_segregation!(
                copy(xm_t2),
                copy(ym_t2),
                copy(tm_t2),
                fill(1580.0, n_tr),
                copy(phim_t2),
                Xfe_mush,
                copy(Xfe_mush),
                n_tr,
                5.0e9,
                cfg_tr;
                coords=coords_tr,
                rplanet=30000.0,
                Fm=fill(0.45, n_tr),
                eta_silicate=1.0e5,
            )
            res_magma_drain = apply_metal_segregation!(
                copy(xm_t2),
                copy(ym_t2),
                copy(tm_t2),
                fill(1750.0, n_tr),
                copy(phim_t2),
                Xfe_liq,
                copy(Xfe_liq),
                n_tr,
                5.0e9,
                cfg_tr;
                coords=coords_tr,
                rplanet=30000.0,
                Fm=fill(0.75, n_tr),
                eta_silicate=1.0e5,
            )
            outer_idx = [
                m for m in 1:n_tr if tm_t2[m] < 3 &&
                    22000.0 <=
                    hypot(xm_t2[m] - coords_tr.xcenter, ym_t2[m] - coords_tr.ycenter) <=
                    29000.0
            ]
            inner_idx = [
                m for m in 1:n_tr if tm_t2[m] < 3 &&
                    hypot(xm_t2[m] - coords_tr.xcenter, ym_t2[m] - coords_tr.ycenter) <=
                    10000.0
            ]
            rim_idx = [
                m for m in 1:n_tr if tm_t2[m] < 3 &&
                    29000.0 <
                    hypot(xm_t2[m] - coords_tr.xcenter, ym_t2[m] - coords_tr.ycenter) <=
                    30000.0
            ]
            w_t2 = [
                2.0 * hypot(xm_t2[m] - coords_tr.xcenter, ym_t2[m] - coords_tr.ycenter) for
                m in 1:n_tr
            ]
            planet_mask_t2 = tm_t2 .< 3
            w_inv_0 = sum(0.08 * w_t2[m] for m in 1:n_tr if planet_mask_t2[m])
            mean_outer_solid = sum(Xfe_solid[outer_idx]) / length(outer_idx)
            mean_outer_mush = sum(Xfe_mush[outer_idx]) / length(outer_idx)
            mean_outer_liq = sum(Xfe_liq[outer_idx]) / length(outer_idx)
            mean_inner_liq = sum(Xfe_liq[inner_idx]) / length(inner_idx)
            @test isfinite(res_cold_trap.max_v_seg) &&
                1.0e-12 < res_cold_trap.max_v_seg < 1.0
            @test isfinite(res_mush_trap.max_v_seg) &&
                1.0e-12 < res_mush_trap.max_v_seg < 1.0
            @test isfinite(res_magma_drain.max_v_seg) &&
                1.0e-12 < res_magma_drain.max_v_seg < 10.0
            @test 0.045 <= mean_outer_solid <= 0.068
            @test mean_outer_solid < 0.075
            @test 0.019 <= mean_outer_mush < mean_outer_solid - 0.005
            @test mean_outer_liq < 0.015
            @test mean_inner_liq > 0.085
            # Rim markers in partially filled boundary cells must drain rather than stall
            @test !isempty(rim_idx)
            @test all(m -> 0.0 <= Xfe_liq[m] < 1.0e-6, rim_idx)
            # In-planet solid-regime cells never drain below phi_residual
            solid_cell_means = [
                sum(
                    Xfe_solid[m] * w_t2[m] for m in 1:n_tr if planet_mask_t2[m] &&
                        clamp(Int(floor(xm_t2[m] / coords_tr.dx)) + 1, 1, coords_tr.Nx) ==
                        jc &&
                        clamp(Int(floor(ym_t2[m] / coords_tr.dy)) + 1, 1, coords_tr.Ny) ==
                        ic
                ) / sum(
                    w_t2[m] for m in 1:n_tr if planet_mask_t2[m] &&
                        clamp(Int(floor(xm_t2[m] / coords_tr.dx)) + 1, 1, coords_tr.Nx) ==
                        jc &&
                        clamp(Int(floor(ym_t2[m] / coords_tr.dy)) + 1, 1, coords_tr.Ny) ==
                        ic
                ) for jc in 1:(coords_tr.Nx) for ic in 1:(coords_tr.Ny) if any(
                    m ->
                        planet_mask_t2[m] &&
                        clamp(Int(floor(xm_t2[m] / coords_tr.dx)) + 1, 1, coords_tr.Nx) ==
                        jc &&
                        clamp(Int(floor(ym_t2[m] / coords_tr.dy)) + 1, 1, coords_tr.Ny) ==
                        ic,
                    1:n_tr,
                )
            ]
            @test all(v -> v >= cfg_tr.phi_residual - 1.0e-12, solid_cell_means)
            for X_chk in (Xfe_solid, Xfe_mush, Xfe_liq)
                @test isapprox(
                    sum(X_chk[m] * w_t2[m] for m in 1:n_tr if planet_mask_t2[m]),
                    w_inv_0;
                    rtol=1e-12,
                )
            end
            @test res_cold_trap.n_subcycles < cfg_tr.max_subcycles

            # Steady-state percolation drainage stops cleanly at phi_residual (trap floor)
            coords_ss = GridCoordinates(
                GridConfig(; Nx=5, Ny=9, xsize=25000.0, ysize=45000.0)
            )
            xm_ss, ym_ss, tm_ss = Float64[], Float64[], Int[]
            tkm_ss, phim_ss, Xfe_ss, Xfem_ss, Fm_ss = Float64[],
            Float64[], Float64[], Float64[],
            Float64[]
            for j in 1:(coords_ss.Nxm), i in 1:(coords_ss.Nym)
                x_v = (j - 0.5) * coords_ss.dxm
                y_v = (i - 0.5) * coords_ss.dym
                i_c = clamp(Int(floor(y_v / coords_ss.dy)) + 1, 1, coords_ss.Ny)
                push!(xm_ss, x_v)
                push!(ym_ss, y_v)
                push!(tm_ss, 1)
                push!(tkm_ss, 1350.0)
                push!(phim_ss, 0.05)
                push!(Fm_ss, 0.0)
                x0 = i_c <= 5 ? 0.08 : 0.02
                push!(Xfe_ss, x0)
                push!(Xfem_ss, x0)
            end
            n_ss = length(xm_ss)
            cfg_ss = CoreFormationConfig(;
                percolation_active=true,
                settling_active=false,
                phi_crit_perc=0.01,
                phi_residual=0.03,
                k_metal_ref=1.0e-7,
            )
            res_ss = apply_metal_segregation!(
                xm_ss,
                ym_ss,
                tm_ss,
                tkm_ss,
                phim_ss,
                Xfe_ss,
                Xfem_ss,
                n_ss,
                1.0e11,
                cfg_ss;
                coords=coords_ss,
                xcenter=coords_ss.xsize / 2,
                ycenter=coords_ss.ysize,
                rplanet=1.0e6,
                gx=zeros(coords_ss.Ny1, coords_ss.Nx1),
                gy=fill(0.2, coords_ss.Ny1, coords_ss.Nx1),
                Fm=Fm_ss,
                spherical_metric=false,
                sensible_heat_exchange=false,
            )
            ss_donor_means = [
                sum(
                    Xfe_ss[m] for
                    m in 1:n_ss if Int(floor(xm_ss[m] / coords_ss.dx)) + 1 == jc &&
                        Int(floor(ym_ss[m] / coords_ss.dy)) + 1 == ic
                ) / count(
                    m ->
                        Int(floor(xm_ss[m] / coords_ss.dx)) + 1 == jc &&
                        Int(floor(ym_ss[m] / coords_ss.dy)) + 1 == ic,
                    1:n_ss,
                ) for ic in 1:5, jc in 1:(coords_ss.Nx - 1)
            ]
            @test res_ss.n_subcycles < cfg_ss.max_subcycles
            @test all(v -> isapprox(v, 0.03; atol=1e-4), ss_donor_means)
        end

        @testset "Partial Metal Melting Solid-Metal Invariant" begin
            cfg_g = GridConfig(; Nx=5, Ny=11, xsize=25000.0, ysize=55000.0)
            coords_pm = GridCoordinates(cfg_g)
            dx_p, dy_p = coords_pm.dx, coords_pm.dy
            dxm_p, dym_p = coords_pm.dxm, coords_pm.dym
            # Fill every basic cell uniformly with sub-grid markers
            xm_p, ym_p, tm_p = Float64[], Float64[], Int[]
            tkm_p, phim_p, Xfe_p, Xfem_p, Fm_p = Float64[],
            Float64[], Float64[], Float64[],
            Float64[]
            for j in 1:(coords_pm.Nxm), i in 1:(coords_pm.Nym)
                x_val = (j - 0.5) * dxm_p
                y_val = (i - 0.5) * dym_p
                push!(xm_p, x_val)
                push!(ym_p, y_val)
                push!(tm_p, 1)
                push!(tkm_p, 1238.0) # mid eutectic
                push!(phim_p, 0.05)
                push!(Fm_p, 0.0)
                i_cell = clamp(Int(floor(y_val / dy_p)) + 1, 1, coords_pm.Ny)
                if i_cell <= 7
                    push!(Xfe_p, 0.30)
                    # Alternate markers: half completely solid (0.0), half 100% molten (0.30) -> cell mean molten 0.15
                    push!(Xfem_p, iseven(i + j) ? 0.30 : 0.0)
                else
                    push!(Xfe_p, 0.05)
                    push!(Xfem_p, 0.05)
                end
            end
            n_pm = length(xm_p)
            gx_p = zeros(coords_pm.Ny1, coords_pm.Nx1)
            gy_p = fill(0.2, coords_pm.Ny1, coords_pm.Nx1)
            cfg_pm = CoreFormationConfig(;
                percolation_active=true,
                settling_active=false,
                phi_crit_perc=0.02,
                phi_residual=0.0,
                k_metal_ref=1.0e-7,
            )
            solid_markers_mask = [
                i <= n_pm && iszero(Xfem_p[i]) && Xfe_p[i] > 0.0 for i in 1:n_pm
            ]
            donor_mask = [Int(floor(ym_p[m] / dy_p)) + 1 <= 7 for m in 1:n_pm]
            recv_mask = [8 <= Int(floor(ym_p[m] / dy_p)) + 1 <= 10 for m in 1:n_pm]
            solid_before = Xfe_p[solid_markers_mask]
            donor_sum_before = sum(Xfe_p[donor_mask])
            recv_sum_before = sum(Xfe_p[recv_mask])
            recv_mean_before = recv_sum_before / count(recv_mask)
            total_before = sum(Xfe_p)

            res_pm = apply_metal_segregation!(
                xm_p,
                ym_p,
                tm_p,
                tkm_p,
                phim_p,
                Xfe_p,
                Xfem_p,
                n_pm,
                1.0e10,
                cfg_pm;
                coords=coords_pm,
                xcenter=coords_pm.xsize / 2,
                ycenter=coords_pm.ysize,
                rplanet=1.0e6,
                gx=gx_p,
                gy=gy_p,
                Fm=Fm_p,
                spherical_metric=false,
                sensible_heat_exchange=false,
            )

            donor_cell_means = [
                sum(
                    Xfe_p[m] for m in 1:n_pm if Int(floor(xm_p[m] / dx_p)) + 1 == jc &&
                        Int(floor(ym_p[m] / dy_p)) + 1 == ic
                ) / count(
                    m ->
                        Int(floor(xm_p[m] / dx_p)) + 1 == jc &&
                        Int(floor(ym_p[m] / dy_p)) + 1 == ic,
                    1:n_pm,
                ) for ic in 1:7, jc in 1:(coords_pm.Nx - 1)
            ]
            donor_sum_after = sum(Xfe_p[donor_mask])
            recv_sum_after = sum(Xfe_p[recv_mask])
            recv_mean_after = recv_sum_after / count(recv_mask)
            @test isfinite(res_pm.max_v_seg) && 1.0e-12 < res_pm.max_v_seg < 1.0
            @test 0.0 <= res_pm.raw_conservation_rel_error < 1.0e-12
            @test all(x -> x >= 0.15 - 1.0e-12, donor_cell_means)
            @test any(x -> x < 0.25, donor_cell_means)
            @test recv_mean_after > recv_mean_before + 0.01
            @test isapprox(
                donor_sum_before - donor_sum_after,
                recv_sum_after - recv_sum_before;
                rtol=1e-11,
            )
            @test all(
                i -> isapprox(Xfe_p[solid_markers_mask][i], solid_before[i]; atol=1e-12),
                eachindex(solid_before),
            )
            @test isapprox(sum(Xfe_p), total_before; rtol=1e-12)

            # With phi_residual > 0, donor equilibrium floor is m_solid + phi_residual * W
            Xfe_pt_res = [Int(floor(ym_p[m] / dy_p)) + 1 <= 7 ? 0.30 : 0.02 for m in 1:n_pm]
            Xfem_pt_res = [
                Int(floor(ym_p[m] / dy_p)) + 1 <= 7 ? (iseven(m) ? 0.30 : 0.0) : 0.02 for
                m in 1:n_pm
            ]
            cfg_pm_res = CoreFormationConfig(;
                percolation_active=true,
                settling_active=false,
                phi_crit_perc=0.01,
                phi_residual=0.02,
                k_metal_ref=1.0e-7,
            )
            res_pm_res = apply_metal_segregation!(
                copy(xm_p),
                copy(ym_p),
                copy(tm_p),
                copy(tkm_p),
                copy(phim_p),
                Xfe_pt_res,
                Xfem_pt_res,
                n_pm,
                1.0e11,
                cfg_pm_res;
                coords=coords_pm,
                xcenter=coords_pm.xsize / 2,
                ycenter=coords_pm.ysize,
                rplanet=1.0e6,
                gx=gx_p,
                gy=gy_p,
                Fm=Fm_p,
                spherical_metric=false,
                sensible_heat_exchange=false,
            )
            donor_res_means = [
                sum(
                    Xfe_pt_res[m] for m in 1:n_pm if Int(floor(xm_p[m] / dx_p)) + 1 == jc &&
                        Int(floor(ym_p[m] / dy_p)) + 1 == ic
                ) / count(
                    m ->
                        Int(floor(xm_p[m] / dx_p)) + 1 == jc &&
                        Int(floor(ym_p[m] / dy_p)) + 1 == ic,
                    1:n_pm,
                ) for ic in 1:5, jc in 1:(coords_pm.Nx - 1)
            ]
            @test res_pm_res.n_subcycles < cfg_pm_res.max_subcycles
            @test all(v -> isapprox(v, 0.17; atol=1e-4), donor_res_means)
        end

        @testset "Godunov Hindered-Settling Monotonicity Above Flux Peak" begin
            @testset "Godunov Flux and Wave Speed Unit Contracts" begin
                v0_u = 1.0e-3
                m_u = 4.5
                pp_u = 0.65
                q_fn(phi) = v0_u * phi * (1.0 - phi / pp_u)^m_u
                phi_star = pp_u / (m_u + 1.0)

                # Rarefaction / expanding shock pd <= pr -> minimum on [pd, pr]
                @test isapprox(
                    Erebus._godunov_settling_flux(0.10, 0.25, v0_u, m_u, pp_u),
                    min(q_fn(0.10), q_fn(0.25));
                    rtol=1e-12,
                )
                # Sonic point straddle pd > phi* > pr -> flux peak q(phi*), distinct from q(pd)
                f_peak = Erebus._godunov_settling_flux(0.30, 0.05, v0_u, m_u, pp_u)
                @test isapprox(f_peak, q_fn(phi_star); rtol=1e-12)
                @test abs(f_peak - q_fn(0.30)) > 0.25 * f_peak
                # Boundary & invalid mobility cases
                @test iszero(Erebus._godunov_settling_flux(pp_u, pp_u, v0_u, m_u, pp_u))
                @test iszero(Erebus._godunov_settling_flux(0.30, 0.05, 0.0, m_u, pp_u))
                # Wave speed across [0.0, 0.30] captures peak derivative |q'(0)| = v0_u
                @test isapprox(
                    Erebus._godunov_settling_wave_speed(0.30, 0.0, v0_u, m_u, pp_u),
                    v0_u;
                    rtol=1e-12,
                )
                # Inflection point straddle [0.12, 0.35] matches fine sampling max |q'|
                ws_inf = Erebus._godunov_settling_wave_speed(0.12, 0.35, v0_u, m_u, pp_u)
                ws_grid = maximum(
                    Erebus._settling_flux_derivative_abs(p, v0_u, m_u, pp_u) for
                    p in range(0.12, 0.35; length=1001)
                )
                @test isapprox(ws_inf, ws_grid; rtol=1e-4)
            end

            cfg_g = GridConfig(; Nx=5, Ny=17, xsize=25000.0, ysize=85000.0)
            coords_gd = GridCoordinates(cfg_g)
            dx_g, dy_g = coords_gd.dx, coords_gd.dy
            dxm_g, dym_g = coords_gd.dxm, coords_gd.dym
            xm_g, ym_g, tm_g = Float64[], Float64[], Int[]
            tkm_g, phim_g, Xfe_g, Xfem_g, Fm_g = Float64[],
            Float64[], Float64[], Float64[],
            Float64[]
            for j in 1:(coords_gd.Nxm), i in 1:(coords_gd.Nym)
                x_val = (j - 0.5) * dxm_g
                y_val = (i - 0.5) * dym_g
                i_cell = clamp(Int(floor(y_val / dy_g)) + 1, 1, coords_gd.Ny)
                push!(xm_g, x_val)
                push!(ym_g, y_val)
                push!(tm_g, 1)
                push!(tkm_g, 1750.0)
                push!(phim_g, 0.05)
                push!(Fm_g, 0.80)
                # Step initial condition at phi0 = 0.30 > phi* = 0.65/5.5 = 0.1182
                x0 = i_cell <= 7 ? 0.30 : 0.0
                push!(Xfe_g, x0)
                push!(Xfem_g, x0)
            end
            n_gd = length(xm_g)
            gx_g = zeros(coords_gd.Ny1, coords_gd.Nx1)
            gy_g = fill(0.1, coords_gd.Ny1, coords_gd.Nx1)
            cfg_gd = CoreFormationConfig(;
                percolation_active=false,
                settling_active=true,
                droplet_size_mode=:fixed,
                droplet_diameter_fixed=5.0e-3,
                hindered_exponent=4.5,
                phi_pack=0.65,
                cfl_settling=0.45,
            )

            col_j = 3
            cell_profile() = [
                sum(
                    Xfe_g[m] for
                    m in 1:n_gd if Int(floor(xm_g[m] / dx_g)) + 1 == col_j &&
                        Int(floor(ym_g[m] / dy_g)) + 1 == ic
                ) / count(
                    m ->
                        Int(floor(xm_g[m] / dx_g)) + 1 == col_j &&
                        Int(floor(ym_g[m] / dy_g)) + 1 == ic,
                    1:n_gd,
                ) for ic in 1:(coords_gd.Ny - 1)
            ]
            prof_init = cell_profile()
            tv_init = sum(abs, diff(prof_init))

            res_gd = apply_metal_segregation!(
                xm_g,
                ym_g,
                tm_g,
                tkm_g,
                phim_g,
                Xfe_g,
                Xfem_g,
                n_gd,
                2.0e10,
                cfg_gd;
                coords=coords_gd,
                xcenter=coords_gd.xsize / 2,
                ycenter=coords_gd.ysize,
                rplanet=1.0e6,
                gx=gx_g,
                gy=gy_g,
                Fm=Fm_g,
                spherical_metric=false,
                sensible_heat_exchange=false,
            )
            prof_final = cell_profile()
            tv_final = sum(abs, diff(prof_final))

            @test 5 <= res_gd.n_subcycles <= cfg_gd.max_subcycles
            @test isfinite(res_gd.max_v_char) && 1.0e-12 < res_gd.max_v_char < 1.0
            @test 0.0 <= res_gd.raw_conservation_rel_error < 1.0e-12
            @test tv_final <= max(tv_init, cfg_gd.phi_pack) + 1.0e-10
            @test maximum(prof_final[1:7]) <= 0.30 + 1.0e-12
            @test tv_init <= 0.35
            @test all(x -> 0.0 <= x <= cfg_gd.phi_pack + 1.0e-10, Xfe_g)

            # Uniform above-peak column stays strictly monotone non-decreasing under Godunov
            Xfe_uni = fill(0.30, n_gd)
            res_uni = apply_metal_segregation!(
                copy(xm_g),
                copy(ym_g),
                copy(tm_g),
                copy(tkm_g),
                copy(phim_g),
                Xfe_uni,
                copy(Xfe_uni),
                n_gd,
                5.0e9,
                cfg_gd;
                coords=coords_gd,
                xcenter=coords_gd.xsize / 2,
                ycenter=coords_gd.ysize,
                rplanet=1.0e6,
                gx=gx_g,
                gy=gy_g,
                Fm=Fm_g,
                spherical_metric=false,
                sensible_heat_exchange=false,
            )
            prof_uni = [
                sum(
                    Xfe_uni[m] for m in 1:n_gd if Int(floor(xm_g[m] / dx_g)) + 1 == col_j &&
                        Int(floor(ym_g[m] / dy_g)) + 1 == ic
                ) / count(
                    m ->
                        Int(floor(xm_g[m] / dx_g)) + 1 == col_j &&
                        Int(floor(ym_g[m] / dy_g)) + 1 == ic,
                    1:n_gd,
                ) for ic in 1:(coords_gd.Ny - 1)
            ]
            @test all(d -> d >= -1.0e-12, diff(prof_uni))
            @test sum(abs, diff(prof_uni)) <= cfg_gd.phi_pack + 1.0e-10
            @test 0.0 <= res_uni.raw_conservation_rel_error < 1.0e-12
        end

        @testset "Characteristic Wave-Speed CFL and Return Contract" begin
            expected_fields = (
                :max_v_seg,
                :max_v_char,
                :n_subcycles,
                :dt_sub,
                :total_dissipation_energy,
                :total_sensible_energy,
                :raw_conservation_rel_error,
            )

            coords_cfl = GridCoordinates(
                GridConfig(; Nx=5, Ny=9, xsize=25000.0, ysize=45000.0)
            )
            dx_c, dy_c = coords_cfl.dx, coords_cfl.dy
            dxm_c, dym_c = coords_cfl.dxm, coords_cfl.dym
            xm_c, ym_c, tm_c = Float64[], Float64[], Int[]
            tkm_c, phim_c, Xfe_c, Xfem_c, Fm_c = Float64[],
            Float64[], Float64[], Float64[],
            Float64[]
            for j in 1:(coords_cfl.Nxm), i in 1:(coords_cfl.Nym)
                x_val = (j - 0.5) * dxm_c
                y_val = (i - 0.5) * dym_c
                i_cell = clamp(Int(floor(y_val / dy_c)) + 1, 1, coords_cfl.Ny)
                push!(xm_c, x_val)
                push!(ym_c, y_val)
                push!(tm_c, 1)
                push!(tkm_c, 1350.0)
                push!(phim_c, 0.05)
                push!(Fm_c, 0.0)
                x0 = i_cell == 4 ? 0.20 : 0.0
                push!(Xfe_c, x0)
                push!(Xfem_c, x0)
            end
            n_c = length(xm_c)
            gx_c = zeros(coords_cfl.Ny1, coords_cfl.Nx1)
            gy_c = fill(0.1, coords_cfl.Ny1, coords_cfl.Nx1)
            cfg_c = CoreFormationConfig(;
                percolation_active=true,
                settling_active=false,
                phi_crit_perc=0.05,
                phi_residual=0.0,
                perm_exponent=3.0,
                k_metal_ref=1.0e-11,
                cfl_settling=0.5,
                max_subcycles=2000,
            )

            # Central finite-difference verification of percolation characteristic wave speed
            q_perc_fn(p) =
                metal_permeability(
                    p;
                    k_metal_ref=cfg_c.k_metal_ref,
                    phi_crit_perc=cfg_c.phi_crit_perc,
                    perm_exponent=cfg_c.perm_exponent,
                    phi0=cfg_c.phi0,
                ) / cfg_c.eta_metal *
                2450.0 *
                0.1
            eps_fd = 1.0e-7
            dq_fd = (q_perc_fn(0.20 + eps_fd) - q_perc_fn(0.20 - eps_fd)) / (2.0 * eps_fd)
            v_perc_020 = q_perc_fn(0.20) / 0.20
            @test isapprox(
                Erebus._percolation_wave_speed(0.20, v_perc_020, cfg_c), dq_fd; rtol=1e-6
            )

            # Inactive path check
            res_off = apply_metal_segregation!(
                copy(xm_c),
                copy(ym_c),
                copy(tm_c),
                copy(tkm_c),
                copy(phim_c),
                copy(Xfe_c),
                copy(Xfem_c),
                n_c,
                1.0e4,
                CoreFormationConfig(; percolation_active=false, settling_active=false);
                coords=coords_cfl,
            )
            @test keys(res_off) == expected_fields

            # Finite sub-unity exponent wave speeds near packing and percolation threshold
            @test isapprox(
                Erebus._settling_flux_derivative_abs(0.20, 1.0e-3, 0.0, 0.65),
                1.0e-3;
                rtol=1e-12,
            )
            ws_sub = Erebus._godunov_settling_wave_speed(0.20, 0.65, 1.0e-3, 0.5, 0.65)
            ws_cell_sub = Erebus._settling_flux_derivative_abs(
                0.65 * (1.0 - 1.0e-10), 1.0e-3, 0.5, 0.65
            )
            @test isfinite(ws_sub) && isapprox(ws_sub, ws_cell_sub; rtol=1e-12)
            cfg_sub_p = CoreFormationConfig(;
                percolation_active=true,
                settling_active=false,
                phi_crit_perc=0.05,
                phi_residual=0.0,
                perm_exponent=0.5,
                k_metal_ref=1.0e-11,
            )
            @test isfinite(
                Erebus._percolation_wave_speed(0.05 + 1.0e-11, 1.0e-10, cfg_sub_p)
            ) && Erebus._percolation_wave_speed(0.05 + 1.0e-11, 1.0e-10, cfg_sub_p) <
                  1.0e-3

            # Active percolation wave speed ratio check and single-substep CFL formula
            res_char = apply_metal_segregation!(
                copy(xm_c),
                copy(ym_c),
                copy(tm_c),
                copy(tkm_c),
                copy(phim_c),
                copy(Xfe_c),
                copy(Xfem_c),
                n_c,
                1.0e5,
                cfg_c;
                coords=coords_cfl,
                xcenter=coords_cfl.xsize / 2,
                ycenter=coords_cfl.ysize,
                rplanet=1.0e6,
                gx=gx_c,
                gy=gy_c,
                Fm=Fm_c,
                spherical_metric=false,
                sensible_heat_exchange=false,
            )
            @test keys(res_char) == expected_fields
            ratio_char = res_char.max_v_char / res_char.max_v_seg
            expected_ratio = 0.20 * (3.0 / 0.15 + 2.0 / 0.85)
            @test isapprox(ratio_char, expected_ratio; rtol=0.02)
            @test res_char.raw_conservation_rel_error < 1.0e-12

            # First substep obeys dt_sub == cfl_settling * dr_min / max_v_char
            cfg_one = CoreFormationConfig(;
                percolation_active=true,
                settling_active=false,
                phi_crit_perc=0.05,
                phi_residual=0.0,
                perm_exponent=3.0,
                k_metal_ref=1.0e-7,
                cfl_settling=0.5,
                max_subcycles=1,
            )
            res_one = @test_logs (:warn, r"max_subcycles") apply_metal_segregation!(
                copy(xm_c),
                copy(ym_c),
                copy(tm_c),
                copy(tkm_c),
                copy(phim_c),
                copy(Xfe_c),
                copy(Xfem_c),
                n_c,
                1.0e10,
                cfg_one;
                coords=coords_cfl,
                xcenter=coords_cfl.xsize / 2,
                ycenter=coords_cfl.ysize,
                rplanet=1.0e6,
                gx=gx_c,
                gy=gy_c,
                Fm=Fm_c,
                spherical_metric=false,
                sensible_heat_exchange=false,
            )
            @test isapprox(
                res_one.dt_sub,
                cfg_one.cfl_settling * min(dx_c, dy_c) / res_one.max_v_char;
                rtol=1e-12,
            )

            # Safe throttle when max_subcycles reached
            cfg_throttle = CoreFormationConfig(;
                percolation_active=true,
                settling_active=false,
                phi_crit_perc=0.05,
                phi_residual=0.0,
                perm_exponent=3.0,
                k_metal_ref=1.0e-7,
                cfl_settling=0.5,
                max_subcycles=2,
            )
            res_th = @test_logs (:warn, r"max_subcycles") apply_metal_segregation!(
                copy(xm_c),
                copy(ym_c),
                copy(tm_c),
                copy(tkm_c),
                copy(phim_c),
                copy(Xfe_c),
                copy(Xfem_c),
                n_c,
                1.0e10,
                cfg_throttle;
                coords=coords_cfl,
                xcenter=coords_cfl.xsize / 2,
                ycenter=coords_cfl.ysize,
                rplanet=1.0e6,
                gx=gx_c,
                gy=gy_c,
                Fm=Fm_c,
                spherical_metric=false,
                sensible_heat_exchange=false,
            )
            @test res_th.n_subcycles == 2
            @test 0.0 < res_th.dt_sub * res_th.n_subcycles < 0.01 * 1.0e10

            # Mush/settling donor over empty magma receiver must stay CFL-bounded on single substep
            for F_donor_cfl in (0.45, 0.60)
                Xfe_mc = Float64[]
                Xfem_mc = Float64[]
                Fm_mc = Float64[]
                top_mask_mc = Bool[]
                for j in 1:(coords_cfl.Nxm), i in 1:(coords_cfl.Nym)
                    y_val = (i - 0.5) * dym_c
                    i_cell = clamp(Int(floor(y_val / dy_c)) + 1, 1, coords_cfl.Ny)
                    in_top = i_cell <= 4
                    push!(top_mask_mc, in_top)
                    push!(Xfe_mc, in_top ? 0.25 : 0.0)
                    push!(Xfem_mc, in_top ? 0.25 : 0.0)
                    push!(Fm_mc, in_top ? F_donor_cfl : 1.0)
                end
                m_top_0 = sum(Xfe_mc[top_mask_mc])
                cfg_mc = CoreFormationConfig(;
                    percolation_active=false,
                    settling_active=true,
                    droplet_size_mode=:fixed,
                    droplet_diameter_fixed=5.0e-3,
                    cfl_settling=0.5,
                    max_subcycles=1,
                )
                res_mc = @test_logs (:warn, r"max_subcycles") apply_metal_segregation!(
                    copy(xm_c),
                    copy(ym_c),
                    copy(tm_c),
                    fill(1750.0, n_c),
                    copy(phim_c),
                    Xfe_mc,
                    Xfem_mc,
                    n_c,
                    1.0e15,
                    cfg_mc;
                    coords=coords_cfl,
                    xcenter=coords_cfl.xsize / 2,
                    ycenter=coords_cfl.ysize,
                    rplanet=1.0e6,
                    gx=gx_c,
                    gy=gy_c,
                    Fm=Fm_mc,
                    eta_silicate=1.0e5,
                    spherical_metric=false,
                    sensible_heat_exchange=false,
                )
                frac_drain = (m_top_0 - sum(Xfe_mc[top_mask_mc])) / m_top_0
                @test 0.0 < frac_drain <= cfg_mc.cfl_settling + 1.0e-10
                @test isfinite(res_mc.max_v_char) && 1.0e-6 < res_mc.max_v_char < 10.0
            end

            # Packed core receiver must not inject its own percolation wave speed into donor CFL
            Xfe_pk = Float64[]
            Fm_pk = zeros(n_c)
            for j in 1:(coords_cfl.Nxm), i in 1:(coords_cfl.Nym)
                y_val = (i - 0.5) * dym_c
                i_cell = clamp(Int(floor(y_val / dy_c)) + 1, 1, coords_cfl.Ny)
                x0 = i_cell <= 4 ? 0.08 : (i_cell == 5 ? 0.64 : 0.65)
                push!(Xfe_pk, x0)
            end
            cfg_pk = CoreFormationConfig(;
                percolation_active=true,
                settling_active=false,
                phi_crit_perc=0.05,
                phi_residual=0.02,
                perm_exponent=3.0,
                k_metal_ref=1.0e-9,
                cfl_settling=0.5,
                max_subcycles=2000,
            )
            res_pk = apply_metal_segregation!(
                copy(xm_c),
                copy(ym_c),
                copy(tm_c),
                copy(tkm_c),
                copy(phim_c),
                Xfe_pk,
                copy(Xfe_pk),
                n_c,
                1.0e11,
                cfg_pk;
                coords=coords_cfl,
                xcenter=coords_cfl.xsize / 2,
                ycenter=coords_cfl.ysize,
                rplanet=1.0e6,
                gx=gx_c,
                gy=gy_c,
                Fm=Fm_pk,
                spherical_metric=false,
                sensible_heat_exchange=false,
                dissipation_heating=false,
            )
            @test isfinite(res_pk.max_v_char) && res_pk.max_v_char < 0.025
            @test res_pk.n_subcycles < 600
            @test isapprox(res_pk.dt_sub * res_pk.n_subcycles, 1.0e11; rtol=1e-9)
            @test iszero(res_pk.total_dissipation_energy)
        end

        @testset "Conservative Two-Phase Counterflow Sensible-Heat Exchange" begin
            cfg_g = GridConfig(; Nx=5, Ny=11, xsize=25000.0, ysize=55000.0)
            coords_sh = GridCoordinates(cfg_g)
            dx_s, dy_s = coords_sh.dx, coords_sh.dy
            dxm_s, dym_s = coords_sh.dxm, coords_sh.dym
            xm_s, ym_s, tm_s = Float64[], Float64[], Int[]
            tkm_s, phim_s, Xfe_s, Xfem_s, Fm_s = Float64[],
            Float64[], Float64[], Float64[],
            Float64[]
            for j in 1:(coords_sh.Nxm), i in 1:(coords_sh.Nym)
                x_val = (j - 0.5) * dxm_s
                y_val = (i - 0.5) * dym_s
                i_cell = clamp(Int(floor(y_val / dy_s)) + 1, 1, coords_sh.Ny)
                push!(xm_s, x_val)
                push!(ym_s, y_val)
                push!(tm_s, 1)
                # Top donor warm (1500 K), bottom receiver cool (1000 K)
                push!(tkm_s, i_cell <= 5 ? 1500.0 : 1000.0)
                push!(phim_s, 0.05)
                push!(Fm_s, 0.0)
                x0 = i_cell <= 5 ? 0.30 : 0.08
                push!(Xfe_s, x0)
                push!(Xfem_s, x0)
            end
            n_s = length(xm_s)
            gx_s = zeros(coords_sh.Ny1, coords_sh.Nx1)
            gy_s = fill(0.1, coords_sh.Ny1, coords_sh.Nx1)
            Q_sens_grid = zeros(coords_sh.Ny1, coords_sh.Nx1)

            cfg_sh = CoreFormationConfig(;
                percolation_active=true,
                settling_active=false,
                phi_crit_perc=0.05,
                phi_residual=0.0,
                k_metal_ref=1.0e-7,
                rhocp_metal=5.0e6,
                cfl_settling=0.5,
                segregation_heating=false,
            )
            rhocp_sil = 3.0e6

            # Record pre-segregation cell enthalpy sum and T_init on basic cells
            dV = dx_s * dy_s
            n_cy = coords_sh.Ny - 1
            n_cx = coords_sh.Nx - 1
            cell_phi0 = zeros(n_cy, n_cx)
            cell_T0 = zeros(n_cy, n_cx)
            cell_cnt = zeros(Int, n_cy, n_cx)
            for m in 1:n_s
                jc = clamp(Int(floor(xm_s[m] / dx_s)) + 1, 1, n_cx)
                ic = clamp(Int(floor(ym_s[m] / dy_s)) + 1, 1, n_cy)
                cell_phi0[ic, jc] += Xfe_s[m]
                cell_T0[ic, jc] += tkm_s[m]
                cell_cnt[ic, jc] += 1
            end
            cell_phi0 ./= cell_cnt
            cell_T0 ./= cell_cnt
            C_init = (1.0 .- cell_phi0) .* rhocp_sil .+ cell_phi0 .* cfg_sh.rhocp_metal
            H_init = sum(C_init .* cell_T0) * dV

            dt_sh = 1.0e9
            res_sh = apply_metal_segregation!(
                xm_s,
                ym_s,
                tm_s,
                tkm_s,
                phim_s,
                Xfe_s,
                Xfem_s,
                n_s,
                dt_sh,
                cfg_sh;
                coords=coords_sh,
                xcenter=coords_sh.xsize / 2,
                ycenter=coords_sh.ysize,
                rplanet=1.0e6,
                gx=gx_s,
                gy=gy_s,
                Fm=Fm_s,
                Q_seg_grid=Q_sens_grid,
                rhocp_silicate=rhocp_sil,
                sensible_heat_exchange=true,
                dissipation_heating=false,
                spherical_metric=false,
            )

            cell_phi1 = zeros(n_cy, n_cx)
            for m in 1:n_s
                jc = clamp(Int(floor(xm_s[m] / dx_s)) + 1, 1, n_cx)
                ic = clamp(Int(floor(ym_s[m] / dy_s)) + 1, 1, n_cy)
                cell_phi1[ic, jc] += Xfe_s[m]
            end
            cell_phi1 ./= cell_cnt
            C_final = (1.0 .- cell_phi1) .* rhocp_sil .+ cell_phi1 .* cfg_sh.rhocp_metal

            Q_cell_int = zeros(n_cy, n_cx)
            for j in 1:n_cx, i in 1:n_cy
                Q_cell_int[i, j] = Q_sens_grid[i + 1, j + 1]
            end
            T_pred = cell_T0 .+ (Q_cell_int .* dt_sh) ./ C_init
            H_final = sum(C_final .* T_pred) * dV
            grid_sens_sum = sum(Q_sens_grid) * dV * dt_sh
            expected_cap_shift = -sum((C_final .- C_init) .* T_pred) * dV

            @test 5 <= res_sh.n_subcycles <= cfg_sh.max_subcycles
            @test isapprox(grid_sens_sum, expected_cap_shift; rtol=1e-10)
            @test isapprox(res_sh.total_sensible_energy, expected_cap_shift; rtol=1e-10)
            @test isapprox(H_final, H_init; rtol=1e-10)
            @test all(T -> 1000.0 - 1.0e-10 <= T <= 1500.0 + 1.0e-10, T_pred)

            # Horizontal E-W gravity counterflow exchange conserves enthalpy cleanly
            cfg_ew = GridConfig(; Nx=11, Ny=5, xsize=55000.0, ysize=25000.0)
            coords_ew = GridCoordinates(cfg_ew)
            n_ew = coords_ew.start_marknum
            xm_ew, ym_ew, tm_ew, tkm_ew, phim_ew, Xfe_ew = Float64[],
            Float64[], Int[], Float64[], Float64[],
            Float64[]
            for j in 1:(coords_ew.Nxm), i in 1:(coords_ew.Nym)
                x_v = (j - 0.5) * coords_ew.dxm
                y_v = (i - 0.5) * coords_ew.dym
                jc = clamp(Int(floor(x_v / coords_ew.dx)) + 1, 1, coords_ew.Nx)
                push!(xm_ew, x_v)
                push!(ym_ew, y_v)
                push!(tm_ew, 1)
                push!(tkm_ew, jc <= 5 ? 1500.0 : 1000.0)
                push!(phim_ew, 0.05)
                push!(Xfe_ew, jc <= 5 ? 0.30 : 0.08)
            end
            n_cy_ew, n_cx_ew = coords_ew.Ny - 1, coords_ew.Nx - 1
            phi0_ew = zeros(n_cy_ew, n_cx_ew)
            T0_ew = zeros(n_cy_ew, n_cx_ew)
            cnt_ew = zeros(Int, n_cy_ew, n_cx_ew)
            for m in 1:n_ew
                jc = clamp(Int(floor(xm_ew[m] / coords_ew.dx)) + 1, 1, n_cx_ew)
                ic = clamp(Int(floor(ym_ew[m] / coords_ew.dy)) + 1, 1, n_cy_ew)
                phi0_ew[ic, jc] += Xfe_ew[m]
                T0_ew[ic, jc] += tkm_ew[m]
                cnt_ew[ic, jc] += 1
            end
            phi0_ew ./= cnt_ew
            T0_ew ./= cnt_ew
            C0_ew = (1.0 .- phi0_ew) .* rhocp_sil .+ phi0_ew .* cfg_sh.rhocp_metal
            dV_ew = coords_ew.dx * coords_ew.dy
            Q_ew = zeros(coords_ew.Ny1, coords_ew.Nx1)
            res_ew = apply_metal_segregation!(
                xm_ew,
                ym_ew,
                tm_ew,
                tkm_ew,
                phim_ew,
                Xfe_ew,
                copy(Xfe_ew),
                n_ew,
                dt_sh,
                cfg_sh;
                coords=coords_ew,
                xcenter=coords_ew.xsize,
                ycenter=coords_ew.ysize / 2,
                rplanet=1.0e6,
                gx=fill(0.1, coords_ew.Ny1, coords_ew.Nx1),
                gy=zeros(coords_ew.Ny1, coords_ew.Nx1),
                Fm=zeros(n_ew),
                Q_seg_grid=Q_ew,
                rhocp_silicate=rhocp_sil,
                sensible_heat_exchange=true,
                dissipation_heating=false,
                spherical_metric=false,
            )
            phi1_ew = zeros(n_cy_ew, n_cx_ew)
            for m in 1:n_ew
                jc = clamp(Int(floor(xm_ew[m] / coords_ew.dx)) + 1, 1, n_cx_ew)
                ic = clamp(Int(floor(ym_ew[m] / coords_ew.dy)) + 1, 1, n_cy_ew)
                phi1_ew[ic, jc] += Xfe_ew[m]
            end
            phi1_ew ./= cnt_ew
            C1_ew = (1.0 .- phi1_ew) .* rhocp_sil .+ phi1_ew .* cfg_sh.rhocp_metal
            Q_cell_ew = [Q_ew[i + 1, j + 1] for i in 1:n_cy_ew, j in 1:n_cx_ew]
            T_pred_ew = T0_ew .+ (Q_cell_ew .* dt_sh) ./ C0_ew
            @test isapprox(
                sum(C1_ew .* T_pred_ew) * dV_ew, sum(C0_ew .* T0_ew) * dV_ew; rtol=1e-10
            )
            @test Q_cell_ew[2, 5] < -1.0e-3 && Q_cell_ew[2, 6] > 1.0e-3
            @test all(T -> 1000.0 - 1.0e-10 <= T <= 1500.0 + 1.0e-10, T_pred_ew)
            @test iszero(res_ew.total_dissipation_energy)

            # Partial-cell w_c/w_full normalization and spherical metric conservation
            coords_part = GridCoordinates(
                GridConfig(; Nx=9, Ny=9, xsize=90000.0, ysize=90000.0)
            )
            n_part = coords_part.start_marknum
            xm_pt, ym_pt, tm_pt, tkm_pt, phim_pt, Xfe_pt = Float64[],
            Float64[], Int[], Float64[], Float64[],
            Float64[]
            rplanet_pt = 36000.0
            for j in 1:(coords_part.Nxm), i in 1:(coords_part.Nym)
                x_v = (j - 0.5) * coords_part.dxm
                y_v = (i - 0.5) * coords_part.dym
                rm = hypot(x_v - coords_part.xcenter, y_v - coords_part.ycenter)
                in_p = rm <= rplanet_pt
                push!(xm_pt, x_v)
                push!(ym_pt, y_v)
                push!(tm_pt, in_p ? 1 : 3)
                push!(tkm_pt, rm <= 18000.0 ? 1500.0 : 1000.0)
                push!(phim_pt, in_p ? 0.05 : 1.0)
                push!(Xfe_pt, !in_p ? 0.0 : (rm > 18000.0 ? 0.25 : 0.05))
            end
            Q_part = zeros(coords_part.Ny1, coords_part.Nx1)
            res_part = apply_metal_segregation!(
                xm_pt,
                ym_pt,
                tm_pt,
                tkm_pt,
                phim_pt,
                Xfe_pt,
                copy(Xfe_pt),
                n_part,
                5.0e8,
                cfg_sh;
                coords=coords_part,
                rplanet=rplanet_pt,
                Fm=zeros(n_part),
                Q_seg_grid=Q_part,
                rhocp_silicate=rhocp_sil,
                sensible_heat_exchange=true,
                dissipation_heating=false,
                spherical_metric=true,
            )
            dV_pt = coords_part.dx * coords_part.dy
            @test 1 <= res_part.n_subcycles <= cfg_sh.max_subcycles
            @test isapprox(
                res_part.total_sensible_energy, sum(Q_part) * dV_pt * 5.0e8; rtol=1e-10
            )

            # Independent w_c/w_full check: half-filled rock+air column scales Q_seg_grid by 0.5 vs full column
            coords_half = GridCoordinates(
                GridConfig(; Nx=5, Ny=5, xsize=20000.0, ysize=20000.0)
            )
            xm_f, ym_f, tm_f, tkm_f, phim_f, Xfe_f = Float64[],
            Float64[], Int[], Float64[], Float64[],
            Float64[]
            xm_h, ym_h, tm_h, tkm_h, phim_h, Xfe_h = Float64[],
            Float64[], Int[], Float64[], Float64[],
            Float64[]
            for j in 1:(coords_half.Nxm), i in 1:(coords_half.Nym)
                x_v = (j - 0.5) * coords_half.dxm
                y_v = (i - 0.5) * coords_half.dym
                i_c = clamp(Int(floor(y_v / coords_half.dy)) + 1, 1, coords_half.Ny)
                t0 = i_c <= 2 ? 1500.0 : 1000.0
                x0 = i_c <= 2 ? 0.25 : 0.05
                push!(xm_f, x_v)
                push!(ym_f, y_v)
                push!(tm_f, 1)
                push!(tkm_f, t0)
                push!(phim_f, 0.05)
                push!(Xfe_f, x0)
                is_rock_half = isodd(i + j)
                push!(xm_h, x_v)
                push!(ym_h, y_v)
                push!(tm_h, is_rock_half ? 1 : 3)
                push!(tkm_h, is_rock_half ? t0 : 200.0)
                push!(phim_h, is_rock_half ? 0.05 : 1.0)
                push!(Xfe_h, is_rock_half ? x0 : 0.0)
            end
            gx_h = zeros(coords_half.Ny1, coords_half.Nx1)
            gy_h = fill(0.1, coords_half.Ny1, coords_half.Nx1)
            Q_full_col = zeros(coords_half.Ny1, coords_half.Nx1)
            Q_half_col = zeros(coords_half.Ny1, coords_half.Nx1)
            res_full_col = apply_metal_segregation!(
                xm_f,
                ym_f,
                tm_f,
                tkm_f,
                phim_f,
                Xfe_f,
                copy(Xfe_f),
                length(xm_f),
                5.0e8,
                cfg_sh;
                coords=coords_half,
                xcenter=coords_half.xsize / 2,
                ycenter=coords_half.ysize,
                rplanet=1.0e6,
                gx=gx_h,
                gy=gy_h,
                Fm=zeros(length(xm_f)),
                Q_seg_grid=Q_full_col,
                rhocp_silicate=rhocp_sil,
                sensible_heat_exchange=true,
                dissipation_heating=false,
                spherical_metric=false,
            )
            res_half_col = apply_metal_segregation!(
                xm_h,
                ym_h,
                tm_h,
                tkm_h,
                phim_h,
                Xfe_h,
                copy(Xfe_h),
                length(xm_h),
                5.0e8,
                cfg_sh;
                coords=coords_half,
                xcenter=coords_half.xsize / 2,
                ycenter=coords_half.ysize,
                rplanet=1.0e6,
                gx=gx_h,
                gy=gy_h,
                Fm=zeros(length(xm_h)),
                Q_seg_grid=Q_half_col,
                rhocp_silicate=rhocp_sil,
                sensible_heat_exchange=true,
                dissipation_heating=false,
                spherical_metric=false,
            )
            @test isapprox(sum(Q_half_col), 0.5 * sum(Q_full_col); rtol=1e-10)
            @test isapprox(
                res_half_col.total_sensible_energy,
                0.5 * res_full_col.total_sensible_energy;
                rtol=1e-10,
            )

            # Integrated solve_thermal_energy! check with segregation_heating = false
            bench_toml = joinpath(
                @__DIR__, "..", "configs", "core_formation_benchmark.toml"
            )
            cfg_b = load_config(bench_toml)
            out_int = mktempdir()
            try
                cf_sens_only = CoreFormationConfig(;
                    percolation_active=true,
                    settling_active=false,
                    phi_crit_perc=0.05,
                    phi_residual=0.02,
                    k_metal_ref=1.0e-9,
                    segregation_heating=false,
                )
                cfg_int = SimulationConfig(;
                    grid=GridConfig(; Nx=17, Ny=17, xsize=140000.0, ysize=140000.0),
                    geometry=cfg_b.geometry,
                    time=TimeConfig(; n_steps=2, dt_initial=5.0e6, dt_longest=5.0e6),
                    solver=cfg_b.solver,
                    poroelasticity=cfg_b.poroelasticity,
                    thermodynamics=cfg_b.thermodynamics,
                    melting=cfg_b.melting,
                    coreformation=cf_sens_only,
                    materials=cfg_b.materials,
                    output=OutputConfig(; output_dir=out_int, savematstep=100),
                )
                st_i, coords_i, ws_i, tel_i, _ = Erebus.init_simulation(
                    cfg_int; output_path=out_int
                )
                tel_i isa IO && close(tel_i)
                st_i.dt = cfg_int.time.dt_initial
                for m in 1:length(st_i.markers)
                    if st_i.markers.core.tm[m] < 3
                        rm = hypot(
                            st_i.markers.core.xm[m] - coords_i.xcenter,
                            st_i.markers.core.ym[m] - coords_i.ycenter,
                        )
                        st_i.markers.core.tkm[m] = rm <= 20000.0 ? 1240.0 : 1350.0
                    end
                end
                Erebus.prepare_step_ambient!(st_i, coords_i, cfg_int, ws_i)
                Erebus.radiogenic_heating!(st_i, coords_i, cfg_int)
                Erebus.snapshot_radiogenic_step_start!(ws_i, st_i, cfg_int)
                Erebus.interpolate_markers_to_grid!(
                    st_i,
                    coords_i,
                    cfg_int;
                    p2m_workspace=ws_i.p2m,
                    thread_buffers=ws_i.thread_buffers,
                    interp_arrays=ws_i.interp_arrays,
                    rock_sums=ws_i.rock_sums,
                )
                Erebus.solve_gravity!(
                    st_i, coords_i, cfg_int; F_grav=ws_i.F_grav, RP=ws_i.RP, SP=ws_i.SP
                )
                Erebus.apply_surface_radiation!(st_i, coords_i, cfg_int)
                Erebus.snapshot_step_start_inventories!(ws_i, st_i, cfg_int)
                Erebus.update_step_start_pressures!(st_i)
                with_logger(NullLogger()) do
                    return Erebus.solve_stokes_darcy!(
                        st_i, coords_i, cfg_int, ws_i; titer=1, dt_step_initial=st_i.dt
                    )
                end
                res_i = Erebus.solve_thermal_energy!(
                    st_i, coords_i, cfg_int, ws_i; titer=1, DHP_pyro=nothing
                )
                @test all(isfinite, st_i.grids.Q_seg_grid) &&
                    maximum(abs, st_i.grids.Q_seg_grid) > 1.0e-8
                @test 0.0 < st_i.accumulators.max_v_seg_prev < 1.0
                @test isfinite(res_i.dt_next)

                cf_off = CoreFormationConfig(;
                    percolation_active=true,
                    settling_active=false,
                    phi_crit_perc=0.99,
                    phi_residual=0.02,
                    k_metal_ref=1.0e-9,
                    segregation_heating=false,
                )
                cfg_off = SimulationConfig(;
                    grid=cfg_int.grid,
                    geometry=cfg_int.geometry,
                    time=cfg_int.time,
                    solver=cfg_int.solver,
                    poroelasticity=cfg_int.poroelasticity,
                    thermodynamics=cfg_int.thermodynamics,
                    melting=cfg_int.melting,
                    coreformation=cf_off,
                    materials=cfg_int.materials,
                    output=cfg_int.output,
                )
                st_n, coords_n, ws_n, tel_n, _ = Erebus.init_simulation(
                    cfg_off; output_path=out_int
                )
                tel_n isa IO && close(tel_n)
                st_n.dt = cfg_off.time.dt_initial
                for m in 1:length(st_n.markers)
                    if st_n.markers.core.tm[m] < 3
                        rm = hypot(
                            st_n.markers.core.xm[m] - coords_n.xcenter,
                            st_n.markers.core.ym[m] - coords_n.ycenter,
                        )
                        st_n.markers.core.tkm[m] = rm <= 20000.0 ? 1240.0 : 1350.0
                    end
                end
                Erebus.prepare_step_ambient!(st_n, coords_n, cfg_off, ws_n)
                Erebus.radiogenic_heating!(st_n, coords_n, cfg_off)
                Erebus.snapshot_radiogenic_step_start!(ws_n, st_n, cfg_off)
                Erebus.interpolate_markers_to_grid!(
                    st_n,
                    coords_n,
                    cfg_off;
                    p2m_workspace=ws_n.p2m,
                    thread_buffers=ws_n.thread_buffers,
                    interp_arrays=ws_n.interp_arrays,
                    rock_sums=ws_n.rock_sums,
                )
                Erebus.solve_gravity!(
                    st_n, coords_n, cfg_off; F_grav=ws_n.F_grav, RP=ws_n.RP, SP=ws_n.SP
                )
                Erebus.apply_surface_radiation!(st_n, coords_n, cfg_off)
                Erebus.snapshot_step_start_inventories!(ws_n, st_n, cfg_off)
                Erebus.update_step_start_pressures!(st_n)
                with_logger(NullLogger()) do
                    return Erebus.solve_stokes_darcy!(
                        st_n, coords_n, cfg_off, ws_n; titer=1, dt_step_initial=st_n.dt
                    )
                end
                Erebus.solve_thermal_energy!(
                    st_n, coords_n, cfg_off, ws_n; titer=1, DHP_pyro=nothing
                )
                @test 0.1 < maximum(abs, st_i.grids.tk2 .- st_n.grids.tk2) < 5.0

                # With segregation_heating=true, total deposited Q_seg_grid adds dissipation on top of sensible exchange
                cf_diss_on = CoreFormationConfig(;
                    percolation_active=true,
                    settling_active=false,
                    phi_crit_perc=0.05,
                    phi_residual=0.02,
                    k_metal_ref=1.0e-9,
                    segregation_heating=true,
                )
                cfg_diss = SimulationConfig(;
                    grid=cfg_int.grid,
                    geometry=cfg_int.geometry,
                    time=cfg_int.time,
                    solver=cfg_int.solver,
                    poroelasticity=cfg_int.poroelasticity,
                    thermodynamics=cfg_int.thermodynamics,
                    melting=cfg_int.melting,
                    coreformation=cf_diss_on,
                    materials=cfg_int.materials,
                    output=OutputConfig(; output_dir=out_int, savematstep=100),
                )
                st_d, coords_d, ws_d, tel_d, _ = Erebus.init_simulation(
                    cfg_diss; output_path=out_int
                )
                tel_d isa IO && close(tel_d)
                st_d.dt = cfg_diss.time.dt_initial
                for m in 1:length(st_d.markers)
                    if st_d.markers.core.tm[m] < 3
                        rm = hypot(
                            st_d.markers.core.xm[m] - coords_d.xcenter,
                            st_d.markers.core.ym[m] - coords_d.ycenter,
                        )
                        st_d.markers.core.tkm[m] = rm <= 20000.0 ? 1240.0 : 1350.0
                    end
                end
                Erebus.prepare_step_ambient!(st_d, coords_d, cfg_diss, ws_d)
                Erebus.radiogenic_heating!(st_d, coords_d, cfg_diss)
                Erebus.snapshot_radiogenic_step_start!(ws_d, st_d, cfg_diss)
                Erebus.interpolate_markers_to_grid!(
                    st_d,
                    coords_d,
                    cfg_diss;
                    p2m_workspace=ws_d.p2m,
                    thread_buffers=ws_d.thread_buffers,
                    interp_arrays=ws_d.interp_arrays,
                    rock_sums=ws_d.rock_sums,
                )
                Erebus.solve_gravity!(
                    st_d, coords_d, cfg_diss; F_grav=ws_d.F_grav, RP=ws_d.RP, SP=ws_d.SP
                )
                Erebus.apply_surface_radiation!(st_d, coords_d, cfg_diss)
                Erebus.snapshot_step_start_inventories!(ws_d, st_d, cfg_diss)
                Erebus.update_step_start_pressures!(st_d)
                with_logger(NullLogger()) do
                    return Erebus.solve_stokes_darcy!(
                        st_d, coords_d, cfg_diss, ws_d; titer=1, dt_step_initial=st_d.dt
                    )
                end
                Erebus.solve_thermal_energy!(
                    st_d, coords_d, cfg_diss, ws_d; titer=1, DHP_pyro=nothing
                )
                @test sum(st_d.grids.Q_seg_grid) > sum(st_i.grids.Q_seg_grid) + 1.0e-10

                # Magma-ocean settling step uses physical suspension viscosity rather than Stokes solver.etamin
                cf_set_only = CoreFormationConfig(;
                    percolation_active=false,
                    settling_active=true,
                    droplet_size_mode=:fixed,
                    droplet_diameter_fixed=5.0e-3,
                    segregation_heating=true,
                )
                cfg_set = SimulationConfig(;
                    grid=cfg_int.grid,
                    geometry=cfg_int.geometry,
                    time=TimeConfig(; n_steps=1, dt_initial=1.0e9, dt_longest=1.0e9),
                    solver=SolverConfig(;
                        etamin=1.0e14,
                        etamax=cfg_b.solver.etamax,
                        hydromech_solver=cfg_b.solver.hydromech_solver,
                    ),
                    poroelasticity=cfg_b.poroelasticity,
                    thermodynamics=cfg_b.thermodynamics,
                    melting=cfg_b.melting,
                    coreformation=cf_set_only,
                    materials=cfg_b.materials,
                    output=OutputConfig(; output_dir=out_int, savematstep=100),
                )
                st_s, coords_s, ws_s, tel_s, _ = Erebus.init_simulation(
                    cfg_set; output_path=out_int
                )
                tel_s isa IO && close(tel_s)
                st_s.dt = cfg_set.time.dt_initial
                for m in 1:length(st_s.markers)
                    if st_s.markers.core.tm[m] < 3
                        st_s.markers.core.tkm[m] = 1750.0
                    end
                end
                Erebus.prepare_step_ambient!(st_s, coords_s, cfg_set, ws_s)
                Erebus.radiogenic_heating!(st_s, coords_s, cfg_set)
                Erebus.snapshot_radiogenic_step_start!(ws_s, st_s, cfg_set)
                Erebus.interpolate_markers_to_grid!(
                    st_s,
                    coords_s,
                    cfg_set;
                    p2m_workspace=ws_s.p2m,
                    thread_buffers=ws_s.thread_buffers,
                    interp_arrays=ws_s.interp_arrays,
                    rock_sums=ws_s.rock_sums,
                )
                Erebus.solve_gravity!(
                    st_s, coords_s, cfg_set; F_grav=ws_s.F_grav, RP=ws_s.RP, SP=ws_s.SP
                )
                Erebus.apply_surface_radiation!(st_s, coords_s, cfg_set)
                Erebus.snapshot_step_start_inventories!(ws_s, st_s, cfg_set)
                Erebus.update_step_start_pressures!(st_s)
                with_logger(NullLogger()) do
                    return Erebus.solve_stokes_darcy!(
                        st_s, coords_s, cfg_set, ws_s; titer=1, dt_step_initial=st_s.dt
                    )
                end
                Erebus.solve_thermal_energy!(
                    st_s, coords_s, cfg_set, ws_s; titer=1, DHP_pyro=nothing
                )
                @test 0.22 < maximum(st_s.markers.Xfe_bulk) <= cf_set_only.phi_pack
                @test 1.0e-7 < st_s.accumulators.max_v_seg_prev < 1.0

                st_m, coords_m, ws_m, tel_m, _ = Erebus.init_simulation(
                    cfg_set; output_path=out_int
                )
                tel_m isa IO && close(tel_m)
                st_m.dt = cfg_set.time.dt_initial
                for m in 1:length(st_m.markers)
                    if st_m.markers.core.tm[m] < 3
                        st_m.markers.core.tkm[m] = 1570.0
                    end
                end
                Erebus.prepare_step_ambient!(st_m, coords_m, cfg_set, ws_m)
                Erebus.radiogenic_heating!(st_m, coords_m, cfg_set)
                Erebus.snapshot_radiogenic_step_start!(ws_m, st_m, cfg_set)
                Erebus.interpolate_markers_to_grid!(
                    st_m,
                    coords_m,
                    cfg_set;
                    p2m_workspace=ws_m.p2m,
                    thread_buffers=ws_m.thread_buffers,
                    interp_arrays=ws_m.interp_arrays,
                    rock_sums=ws_m.rock_sums,
                )
                Erebus.solve_gravity!(
                    st_m, coords_m, cfg_set; F_grav=ws_m.F_grav, RP=ws_m.RP, SP=ws_m.SP
                )
                Erebus.apply_surface_radiation!(st_m, coords_m, cfg_set)
                Erebus.snapshot_step_start_inventories!(ws_m, st_m, cfg_set)
                Erebus.update_step_start_pressures!(st_m)
                with_logger(NullLogger()) do
                    return Erebus.solve_stokes_darcy!(
                        st_m, coords_m, cfg_set, ws_m; titer=1, dt_step_initial=st_m.dt
                    )
                end
                Erebus.solve_thermal_energy!(
                    st_m, coords_m, cfg_set, ws_m; titer=1, DHP_pyro=nothing
                )
                @test 1.0e-18 < st_m.accumulators.max_v_seg_prev < 2.0e-17
            finally
                rm(out_int; recursive=true, force=true)
            end
        end

        @testset "compute_core_radius and Simulation rcore Diagnostics" begin
            @test iszero(compute_core_radius(0.0, 5450.0, 0.65))
            r_expected = 15000.0
            M_exact = (4.0 / 3.0) * π * r_expected^3 * 5450.0 * 0.65
            @test isapprox(
                compute_core_radius(M_exact, 5450.0, 0.65), r_expected; rtol=1e-12
            )
            @test_throws DomainError compute_core_radius(-1.0, 5450.0, 0.65)
            @test_throws DomainError compute_core_radius(M_exact, 0.0, 0.65)
            @test_throws DomainError compute_core_radius(M_exact, 5450.0, 0.0)
            @test_throws DomainError compute_core_radius(NaN, 5450.0, 0.65)

            # Analytic packed spherical core pin test via compute_core_volatile_budgets
            coords_pin = GridCoordinates(
                GridConfig(; Nx=41, Ny=41, xsize=140000.0, ysize=140000.0)
            )
            n_pin = coords_pin.start_marknum
            xm_pin, ym_pin, tm_pin, Xfe_pin, w3d_pin = Float64[],
            Float64[], Int[], Float64[],
            Float64[]
            r_core_known = 20000.0
            for j in 1:(coords_pin.Nxm), i in 1:(coords_pin.Nym)
                x_v = (j - 0.5) * coords_pin.dxm
                y_v = (i - 0.5) * coords_pin.dym
                rm = hypot(x_v - coords_pin.xcenter, y_v - coords_pin.ycenter)
                push!(xm_pin, x_v)
                push!(ym_pin, y_v)
                push!(tm_pin, rm <= 50000.0 ? 1 : 3)
                push!(Xfe_pin, rm <= r_core_known ? 0.65 : (rm <= 50000.0 ? 0.10 : 0.0))
                push!(w3d_pin, 2.0 * rm)
            end
            cb_pin = compute_core_volatile_budgets(
                xm_pin,
                ym_pin,
                tm_pin,
                Xfe_pin,
                nothing,
                nothing,
                nothing,
                nothing,
                n_pin;
                coords=coords_pin,
                w3d_m=w3d_pin,
                rplanet=50000.0,
                rho_metal=5450.0,
                phi_core_threshold=0.40,
            )
            rcore_pin = cbrt(3.0 * cb_pin.V_core_total / (4.0 * π))
            @test isapprox(rcore_pin, r_core_known; rtol=0.05)
            @test isapprox(
                rcore_pin,
                compute_core_radius(cb_pin.M_core_metal, 5450.0, 0.65);
                rtol=1e-12,
            )

            # Full simulation loop coreformation run without metal_partition sets non-zero rcore; pure redox keeps rcore = 0
            bench_toml = joinpath(
                @__DIR__, "..", "configs", "core_formation_benchmark.toml"
            )
            cfg_b = load_config(bench_toml)
            out1 = mktempdir()
            out2 = mktempdir()
            try
                cfg_cf = SimulationConfig(;
                    grid=GridConfig(; Nx=17, Ny=17, xsize=140000.0, ysize=140000.0),
                    geometry=cfg_b.geometry,
                    time=TimeConfig(; n_steps=3, dt_initial=5000.0, dt_longest=5000.0),
                    solver=cfg_b.solver,
                    poroelasticity=cfg_b.poroelasticity,
                    thermodynamics=cfg_b.thermodynamics,
                    melting=cfg_b.melting,
                    coreformation=cfg_b.coreformation,
                    materials=cfg_b.materials,
                    output=OutputConfig(; output_dir=out1, savematstep=1),
                )
                state_cf = run_simulation(cfg_cf)
                @test 0.0 < state_cf.accumulators.rcore < cfg_cf.geometry.rplanet

                cfg_redox_only = SimulationConfig(;
                    grid=GridConfig(; Nx=9, Ny=9, xsize=140000.0, ysize=140000.0),
                    geometry=cfg_b.geometry,
                    time=TimeConfig(; n_steps=2, dt_initial=1.0e4, dt_longest=1.0e4),
                    coreformation=CoreFormationConfig(;
                        percolation_active=false, settling_active=false
                    ),
                    redox=RedoxConfig(; active=true),
                    materials=cfg_b.materials,
                    output=OutputConfig(; output_dir=out2, savematstep=1),
                )
                state_rx = run_simulation(cfg_redox_only)
                @test iszero(state_rx.accumulators.rcore)
            finally
                rm(out1; recursive=true, force=true)
                rm(out2; recursive=true, force=true)
            end
        end
    end
end

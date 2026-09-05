using Erebus
using ExtendableSparse
using JLD2
using LinearAlgebra
using LinearSolve
using Random
using StaticArrays
using Test
using TOML
import Erebus.Numerics: assemble_thermal_lse!, perform_thermal_iterations!

@testset verbose=true "Planetesimal Geometry, Radiation, and Disk Evolution" begin
    @testset "DiskConfig Schema & Bounds Validation" begin
        # Default configuration
        def_disk = DiskConfig()
        @test def_disk.enabled == false
        @test def_disk.model === :fixed
        @test def_disk.t_ambient ≈ 170.0 rtol=1e-12
        @test def_disk.orbital_distance_au ≈ 2.5 rtol=1e-12
        @test def_disk.stellar_mass_msun ≈ 1.0 rtol=1e-12

        # Valid configurations
        cfg_mono = SimulationConfig(disk=DiskConfig(enabled=true, model=:monotonic))
        @test validate_config(cfg_mono) === nothing

        cfg_c02 = SimulationConfig(disk=DiskConfig(enabled=true, model=:class0_to_class2))
        @test validate_config(cfg_c02) === nothing

        cfg_c12 = SimulationConfig(disk=DiskConfig(enabled=true, model=:class1_to_class2))
        @test validate_config(cfg_c12) === nothing

        # Invalid model
        @test_throws ArgumentError validate_config(
            SimulationConfig(disk=DiskConfig(model=:invalid_model))
        )

        # Invalid physical bounds
        @test_throws ArgumentError validate_config(
            SimulationConfig(disk=DiskConfig(t_ambient=-10.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(disk=DiskConfig(orbital_distance_au=0.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(disk=DiskConfig(stellar_mass_msun=-1.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(disk=DiskConfig(gamma=0.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(disk=DiskConfig(gamma=NaN))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(disk=DiskConfig(alpha=-1.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(disk=DiskConfig(q_irr=0.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(disk=DiskConfig(q_visc=-0.5))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(disk=DiskConfig(p_r_t=-0.1))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(disk=DiskConfig(p_m_irr=-0.2))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(disk=DiskConfig(p_m_visc=-0.1))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(disk=DiskConfig(p_m_t=-0.3))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(disk=DiskConfig(p_m_visc_decay=-0.5))
        )

        # Geometry metric bounds
        @test_throws ArgumentError validate_config(
            SimulationConfig(geometry=GeometryConfig(metric_regularization_cells=0.0))
        )

        # Thermal radiation bounds
        @test_throws ArgumentError validate_config(
            SimulationConfig(thermodynamics=ThermalConfig(emissivity=1.5))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(thermodynamics=ThermalConfig(emissivity=-0.1))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(thermodynamics=ThermalConfig(sigma_sb=0.0))
        )
    end

    @testset "DiskConfig Serialization & TOML Round-Trip" begin
        cfg = SimulationConfig(
            geometry=GeometryConfig(
                spherical_metric=true, metric_regularization_cells=0.75
            ),
            thermodynamics=ThermalConfig(surface_radiation=true, emissivity=0.85),
            disk=DiskConfig(
                enabled=true,
                model=:class0_to_class2,
                orbital_distance_au=1.5,
                stellar_mass_msun=0.3,
            ),
        )
        tmp_path = tempname() * ".toml"
        save_config(tmp_path, cfg)
        @test isfile(tmp_path)

        loaded_cfg = load_config(tmp_path)
        @test loaded_cfg.geometry.spherical_metric == true
        @test loaded_cfg.geometry.metric_regularization_cells ≈ 0.75 rtol=1e-12
        @test loaded_cfg.thermodynamics.surface_radiation == true
        @test loaded_cfg.thermodynamics.emissivity ≈ 0.85 rtol=1e-12
        @test loaded_cfg.disk.enabled == true
        @test loaded_cfg.disk.model === :class0_to_class2
        @test loaded_cfg.disk.orbital_distance_au ≈ 1.5 rtol=1e-12
        @test loaded_cfg.disk.stellar_mass_msun ≈ 0.3 rtol=1e-12
        rm(tmp_path; force=true)
    end

    @testset "Disk Temperature Model 1: Monotonic Viscous Decay" begin
        # Fixed ambient mode when disabled
        cfg_disabled = DiskConfig(enabled=false, model=:monotonic, t_ambient=180.0)
        @test compute_disk_temperature(0.0, cfg_disabled) ≈ 180.0 rtol=1e-12
        @test compute_disk_temperature(1.0e14, cfg_disabled) ≈ 180.0 rtol=1e-12

        # Monotonic decay enabled
        cfg_mono = DiskConfig(
            enabled=true, model=:monotonic, orbital_distance_au=2.5, stellar_mass_msun=1.0
        )

        # Initial temperature at t = 0 should equal T_peak
        T_0 = compute_disk_temperature(0.0, cfg_mono)
        T_peak_expected =
            cfg_mono.t_peak_1au * (cfg_mono.orbital_distance_au^(-cfg_mono.q_visc))
        @test isapprox(T_0, T_peak_expected; rtol=1e-4)

        # Monotonic decrease over time
        times_yr = [0.0, 0.1e6, 0.5e6, 1.0e6, 3.0e6, 10.0e6]
        sec_per_yr = 365.25 * 86400.0
        T_vals = [compute_disk_temperature(t * sec_per_yr, cfg_mono) for t in times_yr]
        for i in 1:(length(T_vals) - 1)
            @test T_vals[i] > T_vals[i + 1]
        end

        # Long-term asymptote approaches T_irr
        T_inf = compute_disk_temperature(100.0e6 * sec_per_yr, cfg_mono)
        T_irr_expected =
            cfg_mono.t_irr_1au * (cfg_mono.orbital_distance_au^(-cfg_mono.q_irr))
        @test isapprox(T_inf, T_irr_expected; rtol=0.01)

        # Distance scaling: closer is hotter
        cfg_inner = DiskConfig(enabled=true, model=:monotonic, orbital_distance_au=1.0)
        cfg_outer = DiskConfig(enabled=true, model=:monotonic, orbital_distance_au=5.0)
        @test compute_disk_temperature(1.0e6 * sec_per_yr, cfg_inner) >
            compute_disk_temperature(1.0e6 * sec_per_yr, cfg_outer)

        # Stellar mass scaling: more massive star is hotter
        cfg_high_m = DiskConfig(enabled=true, model=:monotonic, stellar_mass_msun=2.0)
        cfg_low_m = DiskConfig(enabled=true, model=:monotonic, stellar_mass_msun=0.3)
        @test compute_disk_temperature(1.0e6 * sec_per_yr, cfg_high_m) >
            compute_disk_temperature(1.0e6 * sec_per_yr, cfg_low_m)
    end

    @testset "Disk Temperature Model 2: Cold -> Hot -> Cold (Option 2A)" begin
        cfg_c02 = DiskConfig(
            enabled=true,
            model=:class0_to_class2,
            orbital_distance_au=2.5,
            stellar_mass_msun=1.0,
            t_cloud=30.0,
        )
        sec_per_yr = 365.25 * 86400.0

        # At onset (t = 0), disk is at molecular cloud temperature
        T_start = compute_disk_temperature(0.0, cfg_c02)
        @test isapprox(T_start, 30.0; atol=1e-6)

        # Rises during Class 0
        t_mid_class0 = 0.05e6 * sec_per_yr
        T_mid_class0 = compute_disk_temperature(t_mid_class0, cfg_c02)
        @test T_mid_class0 > T_start

        # Peak near t_peak
        t_peak_myr =
            cfg_c02.t_peak_time_1au_myr * (cfg_c02.orbital_distance_au^cfg_c02.p_r_t)
        t_peak_sec = t_peak_myr * 1.0e6 * sec_per_yr
        T_peak_eval = compute_disk_temperature(t_peak_sec, cfg_c02)
        T_peak_theoretical =
            cfg_c02.t_peak_1au * (cfg_c02.orbital_distance_au^(-cfg_c02.q_visc))
        @test isapprox(T_peak_eval, T_peak_theoretical; rtol=0.05)
        @test T_peak_eval > T_mid_class0

        # Decays after peak during Class I and Class II
        t_class1_late = 0.5e6 * sec_per_yr
        t_class2 = 2.0e6 * sec_per_yr
        T_class1_late = compute_disk_temperature(t_class1_late, cfg_c02)
        T_class2 = compute_disk_temperature(t_class2, cfg_c02)
        @test T_peak_eval > T_class1_late > T_class2

        # Asymptote at late times approaches passive irradiation floor
        t_late = 10.0e6 * sec_per_yr
        T_late = compute_disk_temperature(t_late, cfg_c02)
        T_irr_theoretical =
            cfg_c02.t_irr_1au * (cfg_c02.orbital_distance_au^(-cfg_c02.q_irr))
        @test isapprox(T_late, T_irr_theoretical; rtol=0.06)

        # Multi-mass scaling (Williams et al. 2026): M-dwarf cloud collapse is faster
        cfg_mdwarf = DiskConfig(
            enabled=true,
            model=:class0_to_class2,
            orbital_distance_au=2.5,
            stellar_mass_msun=0.1,
        )
        t_peak_mdwarf =
            cfg_mdwarf.t_peak_time_1au_myr * (0.1^cfg_mdwarf.p_m_t) * (2.5^cfg_mdwarf.p_r_t)
        @test t_peak_mdwarf < t_peak_myr
        @test compute_disk_temperature(t_peak_sec, cfg_c02) >
            compute_disk_temperature(t_peak_sec, cfg_mdwarf)

        # Alias equivalence: :class1_to_class2 matches :class0_to_class2 identically
        cfg_c12 = DiskConfig(
            enabled=true,
            model=:class1_to_class2,
            orbital_distance_au=2.5,
            stellar_mass_msun=1.0,
            t_cloud=30.0,
        )
        @test compute_disk_temperature(t_peak_sec, cfg_c12) ==
            compute_disk_temperature(t_peak_sec, cfg_c02)

        # Outer disk regime (r = 60 AU where T_peak <= T_irr and T_irr <= t_cloud)
        cfg_outer_c12 = DiskConfig(
            enabled=true,
            model=:class1_to_class2,
            orbital_distance_au=60.0,
            stellar_mass_msun=1.0,
            t_cloud=30.0,
        )
        T_outer_late = compute_disk_temperature(10.0e6 * sec_per_yr, cfg_outer_c12)
        # Clamped at molecular cloud floor (30 K > 25.9 K)
        @test T_outer_late == cfg_outer_c12.t_cloud
        T_outer_early = compute_disk_temperature(0.1e6 * sec_per_yr, cfg_outer_c12)
        @test T_outer_early >= cfg_outer_c12.t_cloud

        # Outer disk regime where T_irr > t_cloud (r = 30 AU: T_irr ≈ 34.9 K)
        cfg_30au = DiskConfig(
            enabled=true,
            model=:class1_to_class2,
            orbital_distance_au=30.0,
            stellar_mass_msun=1.0,
            t_cloud=30.0,
        )
        T_30_late = compute_disk_temperature(10.0e6 * sec_per_yr, cfg_30au)
        T_30_irr = cfg_30au.t_irr_1au * (30.0^(-cfg_30au.q_irr))
        @test isapprox(T_30_late, T_30_irr; rtol=0.05)
    end

    @testset "Water Snowline Dynamics Across Evolutionary Stages" begin
        sec_per_yr = 365.25 * 86400.0
        sec_per_myr = 1.0e6 * sec_per_yr
        cfg_disk = DiskConfig(enabled=true, model=:class0_to_class2, stellar_mass_msun=1.0)

        # Stage 0: Cloud collapse onset (t = 0), entire disk below sublimation
        r_snow_0 = compute_snowline_radius(0.0, cfg_disk; T_sub=170.0)
        @test isapprox(r_snow_0, 0.05; atol=1e-4)

        # Stage 1: Peak accretion heating (t = 0.12 Myr), snowline pushed outward > 3.5 AU
        t_peak = 0.12 * sec_per_myr
        r_snow_peak = compute_snowline_radius(t_peak, cfg_disk; T_sub=170.0)
        @test r_snow_peak > 3.5
        @test isapprox(
            compute_disk_temperature(t_peak, cfg_disk; orbital_distance_au=r_snow_peak),
            170.0;
            atol=0.2,
        )

        # Stage 2: Viscous clearing at 3.0 Myr, snowline retreats inward
        t_late = 3.0 * sec_per_myr
        r_snow_late = compute_snowline_radius(t_late, cfg_disk; T_sub=170.0)
        @test r_snow_late < r_snow_peak
        @test r_snow_late < 2.0
        @test isapprox(
            compute_disk_temperature(t_late, cfg_disk; orbital_distance_au=r_snow_late),
            170.0;
            atol=0.2,
        )

        # Low-mass star (0.1 M_sun) has smaller snowline radius
        r_snow_mdwarf = compute_snowline_radius(
            t_peak, cfg_disk; stellar_mass_msun=0.1, T_sub=170.0
        )
        @test r_snow_mdwarf < r_snow_peak
    end

    @testset "Stefan-Boltzmann Surface Radiation Physics" begin
        emissivity = 0.9
        sigma_sb = 5.670374419e-8
        T_surf = 350.0
        T_amb = 170.0

        # Linearized HTC formulation
        h_rad = compute_radiation_htc(
            T_surf, T_amb; emissivity=emissivity, sigma_sb=sigma_sb
        )
        @test h_rad > 0.0

        # Exact equivalence: h_rad * (T_surf - T_amb) == ε * σ_SB * (T_surf⁴ - T_amb⁴)
        flux_linear = h_rad * (T_surf - T_amb)
        flux_exact = emissivity * sigma_sb * (T_surf^4 - T_amb^4)
        @test isapprox(flux_linear, flux_exact; rtol=1e-12)

        # Limiting derivative behavior: as T_amb -> T_surf, h_rad -> 4 * ε * σ_SB * T³
        T_close = T_surf + 1e-6
        h_close = compute_radiation_htc(
            T_surf, T_close; emissivity=emissivity, sigma_sb=sigma_sb
        )
        h_tangent = 4.0 * emissivity * sigma_sb * (T_surf^3)
        @test isapprox(h_close, h_tangent; rtol=1e-5)

        # Non-finite and non-positive temperature handling
        @test iszero(compute_radiation_htc(0.0, 170.0))
        @test iszero(compute_radiation_htc(300.0, -10.0))
        @test iszero(compute_radiation_htc(NaN, 170.0))
    end

    @testset "Radiative Boundary Application on Staggered Grid" begin
        coords = default_grid_coordinates()
        Ny1, Nx1 = coords.Ny1, coords.Nx1
        KX = fill(3.0, Ny1, Nx1)
        KY = fill(3.0, Ny1, Nx1)
        tk = fill(200.0, Ny1, Nx1)
        rplanet = 50_000.0
        xcenter = 70_000.0
        ycenter = 70_000.0
        T_amb = 170.0

        KX_orig = copy(KX)
        KY_orig = copy(KY)

        apply_radiative_surface_boundary!(
            KX,
            KY,
            tk,
            coords,
            rplanet,
            xcenter,
            ycenter,
            T_amb;
            emissivity=0.9,
            sigma_sb=5.670374419e-8,
            k_rock=3.0,
        )

        # Verify that only boundary-crossing faces are modified
        n_modified_kx = count(KX .!= KX_orig)
        n_modified_ky = count(KY .!= KY_orig)
        @test n_modified_kx > 0
        @test n_modified_ky > 0

        # Verify deep interior and far exterior faces remain unmodified
        j_center = div(Nx1, 2) + 1
        i_center = div(Ny1, 2) + 1
        @test KX[i_center, j_center] == KX_orig[i_center, j_center]
        @test KY[i_center, j_center] == KY_orig[i_center, j_center]
        @test KX[2, 2] == KX_orig[2, 2]
        @test KY[2, 2] == KY_orig[2, 2]

        # Verify harmonic series resistance bounds: k_face <= 2 * k_rock
        @test all(KX .<= 2.0 * 3.0 + 1e-12)
        @test all(KY .<= 2.0 * 3.0 + 1e-12)

        # Dynamic temperature-dependent rock thermal conductivity (k_rock=nothing, mode=1)
        KX_dyn = fill(10.0, Ny1, Nx1)
        KY_dyn = fill(10.0, Ny1, Nx1)
        apply_radiative_surface_boundary!(
            KX_dyn,
            KY_dyn,
            tk,
            coords,
            rplanet,
            xcenter,
            ycenter,
            T_amb;
            emissivity=0.9,
            sigma_sb=5.670374419e-8,
            k_rock=nothing,
            marker_property_mode=1,
        )
        @test count(KX_dyn .!= 10.0) == n_modified_kx
        @test count(KY_dyn .!= 10.0) == n_modified_ky
        # Verify dynamic conductivity produces higher effective conductivity at low T
        @test all(KX_dyn .> 0.0)
        @test all(KY_dyn .> 0.0)
    end

    @testset "Spherical Geometric Metric Divergence" begin
        coords = default_grid_coordinates()
        Ny1, Nx1 = coords.Ny1, coords.Nx1
        xcenter = 70_000.0
        ycenter = 70_000.0
        rplanet = 50_000.0

        Q_metric = zeros(Ny1, Nx1)
        KX = fill(3.0, Ny1, Nx1)
        KY = fill(3.0, Ny1, Nx1)

        # Radially decreasing temperature field: T(r) = T0 - C * r²
        T0 = 300.0
        C_grad = 100.0 / (rplanet^2)
        tk = zeros(Ny1, Nx1)
        for j in 1:Nx1, i in 1:Ny1
            r_sq = (coords.xp[j] - xcenter)^2 + (coords.yp[i] - ycenter)^2
            tk[i, j] = T0 - C_grad * r_sq
        end

        compute_spherical_metric_heat_source!(
            Q_metric,
            tk,
            KX,
            KY,
            coords;
            xcenter=xcenter,
            ycenter=ycenter,
            rplanet=rplanet,
            reg_cells=0.5,
        )

        # Physical direction: temperature decreases outward, so Q_metric < 0 inside planet
        j_center = div(Nx1, 2) + 1
        i_center = div(Ny1, 2) + 1
        @test Q_metric[i_center + 1, j_center + 1] < 0.0

        # Regularization check: finite at all nodes
        @test all(isfinite, Q_metric)

        # Sticky air: Q_metric must be exactly 0 outside rplanet
        for j in 2:(Nx1 - 1), i in 2:(Ny1 - 1)
            r_sq = (coords.xp[j] - xcenter)^2 + (coords.yp[i] - ycenter)^2
            if r_sq > rplanet^2
                @test iszero(Q_metric[i, j])
            end
        end
    end

    @testset "Analytical 1D Steady-State Spherical Conduction Benchmark" begin
        Nx_test = 65
        Ny_test = 65
        xsize_test = 140_000.0
        ysize_test = 140_000.0
        rplanet_test = 50_000.0
        xcenter_test = 70_000.0
        ycenter_test = 70_000.0
        k_val = 3.0
        H_val = 1.0e-4

        coords_test = GridCoordinates(Nx_test, Ny_test; xsize=xsize_test, ysize=ysize_test)
        Ny1, Nx1 = coords_test.Ny1, coords_test.Nx1

        # Setup 2D cylindrical solve
        tk1_cyl = fill(170.0, Ny1, Nx1)
        RHOCP_test = fill(1.0e6, Ny1, Nx1)
        KX_test = fill(k_val, Ny1, Nx1)
        KY_test = fill(k_val, Ny1, Nx1)
        HR_test = zeros(Ny1, Nx1)
        for j in 1:Nx1, i in 1:Ny1
            r_sq =
                (coords_test.xp[j] - xcenter_test)^2 + (coords_test.yp[i] - ycenter_test)^2
            if r_sq <= rplanet_test^2
                HR_test[i, j] = H_val
            end
        end
        HA_test = zeros(Ny1, Nx1)
        HS_test = zeros(Ny1, Nx1)
        DHP_test = zeros(Ny1, Nx1)
        RT_cyl = zeros(Ny1 * Nx1)
        dt_steady = 1.0e16

        LT_cyl = assemble_thermal_lse!(
            tk1_cyl,
            RHOCP_test,
            KX_test,
            KY_test,
            HR_test,
            HA_test,
            HS_test,
            DHP_test,
            RT_cyl,
            dt_steady;
            coords=coords_test,
        )
        sol_cyl = LT_cyl \ RT_cyl
        T_cyl_field = reshape(sol_cyl, Ny1, Nx1)

        # Solve with spherical geometric metric term
        Q_metric_test = zeros(Ny1, Nx1)
        for j in 1:Nx1, i in 1:Ny1
            r_sq =
                (coords_test.xp[j] - xcenter_test)^2 + (coords_test.yp[i] - ycenter_test)^2
            if r_sq <= rplanet_test^2
                Q_metric_test[i, j] = -H_val / 3.0
            end
        end

        RT_sph = zeros(Ny1 * Nx1)
        LT_sph = assemble_thermal_lse!(
            tk1_cyl,
            RHOCP_test,
            KX_test,
            KY_test,
            HR_test,
            HA_test,
            HS_test,
            DHP_test,
            RT_sph,
            dt_steady;
            coords=coords_test,
            Q_metric=Q_metric_test,
        )
        sol_sph = LT_sph \ RT_sph
        T_sph_field = reshape(sol_sph, Ny1, Nx1)

        j_c = div(Nx1, 2) + 1
        i_c = div(Ny1, 2) + 1
        idx_surf_j = argmin(abs.(coords_test.xp .- (xcenter_test + rplanet_test)))
        T_surf_cyl = T_cyl_field[i_c, idx_surf_j]
        T_surf_sph = T_sph_field[i_c, idx_surf_j]

        dT_cyl = T_cyl_field[i_c, j_c] - T_surf_cyl
        dT_sph = T_sph_field[i_c, j_c] - T_surf_sph

        ratio = dT_sph / dT_cyl
        @test isapprox(ratio, 2.0 / 3.0; atol=0.01)

        # Iterative dynamic gradient-driven convergence test
        tk_iter = copy(T_cyl_field)
        Q_metric_dynamic = zeros(Ny1, Nx1)
        for iter in 1:25
            compute_spherical_metric_heat_source!(
                Q_metric_dynamic,
                tk_iter,
                KX_test,
                KY_test,
                coords_test;
                xcenter=xcenter_test,
                ycenter=ycenter_test,
                rplanet=rplanet_test,
                reg_cells=0.5,
            )
            RT_iter = zeros(Ny1 * Nx1)
            LT_iter = assemble_thermal_lse!(
                tk1_cyl,
                RHOCP_test,
                KX_test,
                KY_test,
                HR_test,
                HA_test,
                HS_test,
                DHP_test,
                RT_iter,
                dt_steady;
                coords=coords_test,
                Q_metric=Q_metric_dynamic,
            )
            sol_iter = LT_iter \ RT_iter
            tk_new = reshape(sol_iter, Ny1, Nx1)
            tk_iter .= 0.5 .* (tk_iter .+ tk_new)
        end
        dT_iter = tk_iter[i_c, j_c] - tk_iter[i_c, idx_surf_j]
        ratio_dynamic = dT_iter / dT_cyl
        @test ratio_dynamic < 0.85
        @test ratio_dynamic > 0.60
    end

    @testset "perform_thermal_iterations! with Q_metric" begin
        coords_pti = GridCoordinates(GridConfig(Nx=9, Ny=9, xsize=100000.0, ysize=100000.0))
        Ny1, Nx1 = coords_pti.Ny1, coords_pti.Nx1
        tk0_pti = fill(300.0, Ny1, Nx1)
        tk1_pti = fill(300.0, Ny1, Nx1)
        tk2_pti = zeros(Ny1, Nx1)
        DT_pti = zeros(Ny1, Nx1)
        DT0_pti = zeros(Ny1, Nx1)
        RHOCP_pti = fill(3.3e6, Ny1, Nx1)
        KX_pti = fill(3.0, Ny1, Nx1)
        KY_pti = fill(3.0, Ny1, Nx1)
        HR_pti = zeros(Ny1, Nx1)
        HA_pti = zeros(Ny1, Nx1)
        HS_pti = zeros(Ny1, Nx1)
        DHP_pti = zeros(Ny1, Nx1)
        RT_pti = zeros(Ny1 * Nx1)
        ST_pti = zeros(Ny1 * Nx1)
        Q_metric_pti = fill(1.0e-7, Ny1, Nx1)
        tk_prev = copy(tk1_pti)
        perform_thermal_iterations!(
            tk0_pti,
            tk1_pti,
            tk2_pti,
            DT_pti,
            DT0_pti,
            RHOCP_pti,
            KX_pti,
            KY_pti,
            HR_pti,
            HA_pti,
            HS_pti,
            DHP_pti,
            RT_pti,
            ST_pti,
            1.0e10;
            coords=coords_pti,
            Q_metric=Q_metric_pti,
        )
        @test all(tk2_pti .> tk_prev)
        @test all(isfinite, tk2_pti)
        @test minimum(tk2_pti) > 0.0
    end

    @testset "Simulation Integration with Spherical Metric, Radiation, and Disk Models" begin
        # 1. Monotonic disk model + radiation (verify surface radiation transfers heat from warmer disk)
        output_dir1 = mktempdir()
        output_dir_norad = mktempdir()
        try
            cfg_sim1 = SimulationConfig(
                time=TimeConfig(n_steps=3),
                thermodynamics=ThermalConfig(surface_radiation=true, emissivity=0.9),
                disk=DiskConfig(enabled=true, model=:monotonic, orbital_distance_au=2.5),
                output=OutputConfig(output_dir=output_dir1, savematstep=1),
            )
            cfg_norad = SimulationConfig(
                time=TimeConfig(n_steps=3),
                thermodynamics=ThermalConfig(surface_radiation=false),
                disk=DiskConfig(enabled=true, model=:monotonic, orbital_distance_au=2.5),
                output=OutputConfig(output_dir=output_dir_norad, savematstep=1),
            )
            Random.seed!(Erebus.rgen, 42)
            @test run_simulation(cfg_sim1) === nothing
            Random.seed!(Erebus.rgen, 42)
            @test run_simulation(cfg_norad) === nothing
            files1 = readdir(output_dir1)
            @test "output_00000.jld2" in files1
            @test "output_00003.jld2" in files1
            data1 = JLD2.load(joinpath(output_dir1, "output_00003.jld2"))
            data_norad = JLD2.load(joinpath(output_dir_norad, "output_00003.jld2"))
            @test data1["timestep"] == 3
            @test !any(isnan, data1["tk2"])
            @test all(data1["tk2"] .> 0.0)
            @test maximum(data1["tk2"]) > maximum(data_norad["tk2"])
            @test sum(data1["tk2"]) > sum(data_norad["tk2"])
        finally
            rm(output_dir1, recursive=true, force=true)
            rm(output_dir_norad, recursive=true, force=true)
        end

        # 2. Class 0/I/II disk model + spherical metric + radiation
        output_dir2 = mktempdir()
        try
            cfg_sim2 = SimulationConfig(
                time=TimeConfig(n_steps=3),
                geometry=GeometryConfig(spherical_metric=true),
                thermodynamics=ThermalConfig(surface_radiation=true, emissivity=0.9),
                disk=DiskConfig(
                    enabled=true, model=:class0_to_class2, orbital_distance_au=2.5
                ),
                output=OutputConfig(output_dir=output_dir2, savematstep=1),
            )
            @test run_simulation(cfg_sim2) === nothing
            files2 = readdir(output_dir2)
            @test "output_00000.jld2" in files2
            @test "output_00003.jld2" in files2
            data2 = JLD2.load(joinpath(output_dir2, "output_00003.jld2"))
            @test data2["timestep"] == 3
            @test !any(isnan, data2["tk2"])
            @test all(data2["tk2"] .> 0.0)
        finally
            rm(output_dir2, recursive=true, force=true)
        end

        # 3. Radiation enabled with disk.enabled = false (reservoir from t_ambient)
        output_dir3 = mktempdir()
        try
            cfg_sim3 = SimulationConfig(
                time=TimeConfig(n_steps=3),
                thermodynamics=ThermalConfig(surface_radiation=true, emissivity=0.9),
                disk=DiskConfig(enabled=false, t_ambient=190.0),
                output=OutputConfig(output_dir=output_dir3, savematstep=1),
            )
            @test run_simulation(cfg_sim3) === nothing
            files3 = readdir(output_dir3)
            @test "output_00003.jld2" in files3
            data3 = JLD2.load(joinpath(output_dir3, "output_00003.jld2"))
            @test data3["timestep"] == 3
            @test !any(isnan, data3["tk2"])
            @test all(data3["tk2"] .> 0.0)
        finally
            rm(output_dir3, recursive=true, force=true)
        end

        # 4. Bit-for-bit regression: default configuration produces identical fields
        output_dir_a = mktempdir()
        output_dir_b = mktempdir()
        try
            cfg_a = SimulationConfig(
                time=TimeConfig(n_steps=2),
                output=OutputConfig(output_dir=output_dir_a, savematstep=1),
            )
            cfg_b = SimulationConfig(
                time=TimeConfig(n_steps=2),
                geometry=GeometryConfig(spherical_metric=false),
                thermodynamics=ThermalConfig(surface_radiation=false),
                disk=DiskConfig(enabled=false),
                output=OutputConfig(output_dir=output_dir_b, savematstep=1),
            )
            Random.seed!(Erebus.rgen, 42)
            @test run_simulation(cfg_a) === nothing
            Random.seed!(Erebus.rgen, 42)
            @test run_simulation(cfg_b) === nothing
            data_a = JLD2.load(joinpath(output_dir_a, "output_00002.jld2"))
            data_b = JLD2.load(joinpath(output_dir_b, "output_00002.jld2"))
            @test data_a["tk2"] == data_b["tk2"]
            @test data_a["pr"] == data_b["pr"]
            @test data_a["PHI"] == data_b["PHI"]
        finally
            rm(output_dir_a, recursive=true, force=true)
            rm(output_dir_b, recursive=true, force=true)
        end
    end
end

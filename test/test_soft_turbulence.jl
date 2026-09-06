using Test
using Erebus
using Erebus.Config
using Erebus.Physics
using Erebus.Particles
using Erebus.Geometry
using StaticArrays
using TOML
using JLD2

@testset "Sub-Grid Soft Turbulence & Regularized Conductivity" begin
    @testset "Config Schema & Validation" begin
        # Default configuration has soft turbulence disabled
        melt_def = MeltingConfig()
        @test melt_def.soft_turbulence == false
        @test isapprox(melt_def.turb_exponent, 1.0 / 3.0; rtol=1e-12)
        @test isapprox(melt_def.eta_fluid_silicate, 100.0; rtol=1e-12)
        @test isapprox(melt_def.F_turb_start, 0.30; rtol=1e-12)
        @test isapprox(melt_def.F_turb_end, 0.50; rtol=1e-12)
        @test isapprox(melt_def.dT_turb_min, 10.0; rtol=1e-12)
        @test isapprox(melt_def.T_surface_ref, 300.0; rtol=1e-12)
        @test isapprox(melt_def.k_turb_cutoff, 1.0e6; rtol=1e-12)
        @test isapprox(melt_def.k_turb_floor, 1.0e-3; rtol=1e-12)

        # Valid active soft turbulence configuration
        cfg_valid = SimulationConfig(;
            melting=MeltingConfig(;
                active=true,
                soft_turbulence=true,
                turb_exponent=1.0 / 3.0,
                eta_fluid_silicate=100.0,
                F_turb_start=0.30,
                F_turb_end=0.50,
                dT_turb_min=10.0,
                T_surface_ref=300.0,
                k_turb_cutoff=1.0e6,
                k_turb_floor=1.0e-3,
            ),
        )
        @test validate_config(cfg_valid) === nothing

        # Invalid: eta_fluid_silicate <= 0 or non-finite
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(;
                    active=true, soft_turbulence=true, eta_fluid_silicate=0.0
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(;
                    active=true, soft_turbulence=true, eta_fluid_silicate=-10.0
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(;
                    active=true, soft_turbulence=true, eta_fluid_silicate=NaN
                ),
            ),
        )

        # Invalid: F_turb_start < 0 or F_turb_start >= F_turb_end or F_turb_end > 1
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(;
                    active=true, soft_turbulence=true, F_turb_start=-0.1
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(;
                    active=true, soft_turbulence=true, F_turb_start=0.5, F_turb_end=0.4
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(; active=true, soft_turbulence=true, F_turb_end=1.2)
            ),
        )

        # Invalid: turb_exponent <= 0 or non-finite
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(;
                    active=true, soft_turbulence=true, turb_exponent=0.0
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(;
                    active=true, soft_turbulence=true, turb_exponent=-0.5
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(;
                    active=true, soft_turbulence=true, turb_exponent=NaN
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(;
                    active=true, soft_turbulence=true, turb_exponent=Inf
                ),
            ),
        )

        # Invalid: dT_turb_min <= 0 or non-finite
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(; active=true, soft_turbulence=true, dT_turb_min=0.0)
            ),
        )

        # Invalid: T_surface_ref <= 0 or non-finite
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(;
                    active=true, soft_turbulence=true, T_surface_ref=-50.0
                ),
            ),
        )

        # Invalid: k_turb_cutoff <= k_turb_floor
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(;
                    active=true, soft_turbulence=true, k_turb_cutoff=1.0, k_turb_floor=10.0
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(;
                    active=true, soft_turbulence=true, k_turb_floor=-1.0
                ),
            ),
        )
    end

    @testset "TOML Round-Trip Serialization" begin
        cfg_custom = SimulationConfig(;
            melting=MeltingConfig(;
                active=true,
                soft_turbulence=true,
                turb_exponent=0.5,
                eta_fluid_silicate=50.0,
                F_turb_start=0.25,
                F_turb_end=0.55,
                dT_turb_min=15.0,
                T_surface_ref=280.0,
                k_turb_cutoff=5.0e5,
                k_turb_floor=5.0e-4,
            ),
        )
        toml_str = save_config(cfg_custom)
        cfg_reloaded = load_config(toml_str)
        @test cfg_reloaded.melting.soft_turbulence == true
        @test isapprox(cfg_reloaded.melting.turb_exponent, 0.5; rtol=1e-12)
        @test isapprox(cfg_reloaded.melting.eta_fluid_silicate, 50.0; rtol=1e-12)
        @test isapprox(cfg_reloaded.melting.F_turb_start, 0.25; rtol=1e-12)
        @test isapprox(cfg_reloaded.melting.F_turb_end, 0.55; rtol=1e-12)
        @test isapprox(cfg_reloaded.melting.dT_turb_min, 15.0; rtol=1e-12)
        @test isapprox(cfg_reloaded.melting.T_surface_ref, 280.0; rtol=1e-12)
        @test isapprox(cfg_reloaded.melting.k_turb_cutoff, 5.0e5; rtol=1e-12)
        @test isapprox(cfg_reloaded.melting.k_turb_floor, 5.0e-4; rtol=1e-12)
    end

    @testset "Physics: regularized_soft_turbulence_conductivity" begin
        k_cond = 3.0
        eta_num = 1.0e12
        eta_fluid = 100.0
        # Default exponent 1/3: k_turb = 3.0 * (1e12 / 100)^(1/3)
        k_turb_expected = 3.0 * (1.0e12 / 100.0)^(1.0 / 3.0)
        # Exponent 0.5: k_turb_05 = 3.0 * sqrt(1e12 / 100) = 3.0e5 W/(m K)
        k_turb_05 = 3.0 * sqrt(1.0e12 / 100.0)
        @test isapprox(k_turb_05, 3.0e5; rtol=1e-12)

        # Exponent scaling check: 1/3 vs 1/2
        k_full_third = regularized_soft_turbulence_conductivity(
            k_cond,
            eta_num,
            eta_fluid,
            0.60,
            1500.0,
            300.0;
            turb_exponent=1.0 / 3.0,
            F_start=0.30,
            F_end=0.50,
            dT_min=10.0,
        )
        @test isapprox(k_full_third, k_turb_expected; rtol=1e-12)

        k_full_half = regularized_soft_turbulence_conductivity(
            k_cond,
            eta_num,
            eta_fluid,
            0.60,
            1500.0,
            300.0;
            turb_exponent=0.5,
            F_start=0.30,
            F_end=0.50,
            dT_min=10.0,
        )
        @test isapprox(k_full_half, k_turb_05; rtol=1e-12)

        # Domain errors for non-finite, negative, and invalid parameters
        @test_throws DomainError regularized_soft_turbulence_conductivity(
            0.0, eta_num, eta_fluid, 0.5, 1500.0, 300.0
        )
        @test_throws DomainError regularized_soft_turbulence_conductivity(
            -3.0, eta_num, eta_fluid, 0.5, 1500.0, 300.0
        )
        @test_throws DomainError regularized_soft_turbulence_conductivity(
            NaN, eta_num, eta_fluid, 0.5, 1500.0, 300.0
        )
        @test_throws DomainError regularized_soft_turbulence_conductivity(
            k_cond, 0.0, eta_fluid, 0.5, 1500.0, 300.0
        )
        @test_throws DomainError regularized_soft_turbulence_conductivity(
            k_cond, -1.0, eta_fluid, 0.5, 1500.0, 300.0
        )
        @test_throws DomainError regularized_soft_turbulence_conductivity(
            k_cond, NaN, eta_fluid, 0.5, 1500.0, 300.0
        )
        @test_throws DomainError regularized_soft_turbulence_conductivity(
            k_cond, eta_num, 0.0, 0.5, 1500.0, 300.0
        )
        @test_throws DomainError regularized_soft_turbulence_conductivity(
            k_cond, eta_num, -10.0, 0.5, 1500.0, 300.0
        )
        @test_throws DomainError regularized_soft_turbulence_conductivity(
            k_cond, eta_num, NaN, 0.5, 1500.0, 300.0
        )
        @test_throws DomainError regularized_soft_turbulence_conductivity(
            k_cond, eta_num, eta_fluid, NaN, 1500.0, 300.0
        )
        @test_throws DomainError regularized_soft_turbulence_conductivity(
            k_cond, eta_num, eta_fluid, 0.5, NaN, 300.0
        )
        @test_throws DomainError regularized_soft_turbulence_conductivity(
            k_cond, eta_num, eta_fluid, 0.5, 1500.0, NaN
        )
        @test_throws DomainError regularized_soft_turbulence_conductivity(
            k_cond, eta_num, eta_fluid, 0.5, 1500.0, 300.0; turb_exponent=0.0
        )
        @test_throws DomainError regularized_soft_turbulence_conductivity(
            k_cond, eta_num, eta_fluid, 0.5, 1500.0, 300.0; turb_exponent=-0.5
        )
        @test_throws DomainError regularized_soft_turbulence_conductivity(
            k_cond, eta_num, eta_fluid, 0.5, 1500.0, 300.0; turb_exponent=NaN
        )
        @test_throws DomainError regularized_soft_turbulence_conductivity(
            k_cond, eta_num, eta_fluid, 0.5, 1500.0, 300.0; F_start=0.5, F_end=0.3
        )

        # Invariant: k_eff >= k_cond everywhere, even when eta_num < eta_fluid
        k_low_visc = regularized_soft_turbulence_conductivity(
            k_cond, 10.0, 100.0, 0.40, 1500.0, 300.0
        )
        @test k_low_visc >= k_cond

        # Asymptotic limit 1: Below F_turb_start (solid / sub-threshold rock)
        # Should return k_cond exactly
        k_sub = regularized_soft_turbulence_conductivity(
            k_cond, eta_num, eta_fluid, 0.20, 1500.0, 300.0; F_start=0.30, F_end=0.50
        )
        @test isapprox(k_sub, k_cond; rtol=1e-12)

        # Asymptotic limit 2: Zero temperature contrast (isothermal)
        # Should return k_cond exactly
        k_iso = regularized_soft_turbulence_conductivity(
            k_cond,
            eta_num,
            eta_fluid,
            0.80,
            300.0,
            300.0;
            F_start=0.30,
            F_end=0.50,
            dT_min=10.0,
        )
        # Asymptotic limit 2b: Partial thermal contrast (0 < dT < dT_min)
        # Weight w_T = (dT / dT_min)^2
        k_partial_dt = regularized_soft_turbulence_conductivity(
            k_cond,
            eta_num,
            eta_fluid,
            1.0,
            305.0,
            300.0;
            turb_exponent=1.0 / 3.0,
            F_start=0.30,
            F_end=0.50,
            dT_min=10.0,
        )
        # w_F = 1.0, w_T = (5/10)^2 = 0.25, w = 0.25
        k_expected_partial = 10.0^(0.75 * log10(k_cond) + 0.25 * log10(k_turb_expected))
        @test isapprox(k_partial_dt, k_expected_partial; rtol=1e-10)

        # Test at dT = 2 K (w_T = 0.04)
        k_partial_dt2 = regularized_soft_turbulence_conductivity(
            k_cond,
            eta_num,
            eta_fluid,
            1.0,
            302.0,
            300.0;
            turb_exponent=1.0 / 3.0,
            F_start=0.30,
            F_end=0.50,
            dT_min=10.0,
        )
        k_expected_partial2 = 10.0^(0.96 * log10(k_cond) + 0.04 * log10(k_turb_expected))
        @test isapprox(k_partial_dt2, k_expected_partial2; rtol=1e-10)

        # Asymptotic limit 3: Fully molten (F_m >= F_turb_end) and large dT
        # Weight w = 1.0, so should return k_turb_expected exactly
        k_full = regularized_soft_turbulence_conductivity(
            k_cond,
            eta_num,
            eta_fluid,
            0.60,
            1500.0,
            300.0;
            turb_exponent=1.0 / 3.0,
            F_start=0.30,
            F_end=0.50,
            dT_min=10.0,
        )
        @test isapprox(k_full, k_turb_expected; rtol=1e-12)

        # Smoothness and Monotonicity across transition window [0.30, 0.50]
        F_vals = 0.28:0.02:0.52
        k_vals = [
            regularized_soft_turbulence_conductivity(
                k_cond,
                eta_num,
                eta_fluid,
                F,
                1500.0,
                300.0;
                turb_exponent=1.0 / 3.0,
                F_start=0.30,
                F_end=0.50,
                dT_min=10.0,
            ) for F in F_vals
        ]

        # Strictly monotonic non-decreasing
        for idx in 1:(length(k_vals) - 1)
            @test k_vals[idx + 1] >= k_vals[idx]
        end

        # Value at midpoint F = 0.40 (smoothstep xi = 0.5, w = 0.5)
        # Logarithmic midpoint: log10(k_mid) = 0.5 * log10(3) + 0.5 * log10(k_turb_expected)
        k_mid = regularized_soft_turbulence_conductivity(
            k_cond,
            eta_num,
            eta_fluid,
            0.40,
            1500.0,
            300.0;
            turb_exponent=1.0 / 3.0,
            F_start=0.30,
            F_end=0.50,
            dT_min=10.0,
        )
        @test isapprox(k_mid, sqrt(k_cond * k_turb_expected); rtol=1e-10)

        # Clamping at k_cutoff
        # For liquid iron: eta_fluid = 0.01, k_turb = 3 * (1e12 / 0.01)^0.5 = 3e7 > 1e6
        k_clamped = regularized_soft_turbulence_conductivity(
            k_cond, 1.0e12, 0.01, 1.0, 1500.0, 300.0; turb_exponent=0.5, k_cutoff=1.0e6
        )
        @test isapprox(k_clamped, 1.0e6; rtol=1e-12)

        # Clamping at k_floor
        k_floored = regularized_soft_turbulence_conductivity(
            1.0e-5, 1.0, 1.0, 0.0, 1500.0, 300.0; k_floor=1.0e-3
        )
        @test isapprox(k_floored, 1.0e-3; rtol=1e-12)
    end

    @testset "Particle Integration: compute_marker_properties!" begin
        marknum = 4
        tm = zeros(Int64, marknum)
        tkm = zeros(Float64, marknum)
        rhototalm = zeros(Float64, marknum)
        rhocptotalm = zeros(Float64, marknum)
        etatotalm = zeros(Float64, marknum)
        hrtotalm = zeros(Float64, marknum)
        ktotalm_off = zeros(Float64, marknum)
        ktotalm_on = zeros(Float64, marknum)
        tkm_rhocptotalm = zeros(Float64, marknum)
        etafluidcur_inv_kphim = zeros(Float64, marknum)
        phim = zeros(Float64, marknum)
        XWsolidm0 = zeros(Float64, marknum)
        pfm0 = zeros(Float64, marknum)
        Fm = zeros(Float64, marknum)
        rhofluidcur = zeros(Float64, marknum)

        tm[1] = 1 # rock, partially molten (~0.45)
        tkm[1] = 1580.0
        phim[1] = 0.05

        tm[2] = 1 # rock, cold sub-solidus
        tkm[2] = 1000.0
        phim[2] = 0.05

        tm[3] = 1 # rock, fully molten (1900 K, F_m = 1.0)
        tkm[3] = 1900.0
        phim[3] = 0.05

        tm[4] = 3 # sticky air
        tkm[4] = 300.0
        phim[4] = 0.99

        hrsolidm = SVector{3,Float64}([0.0, 0.0, 0.0])
        hrfluidm = SVector{3,Float64}([0.0, 0.0, 0.0])
        mode = 9

        # Case 1: soft_turbulence = false
        for m in 1:marknum
            compute_marker_properties!(
                m,
                tm,
                tkm,
                rhototalm,
                rhocptotalm,
                etatotalm,
                hrtotalm,
                ktotalm_off,
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
                soft_turbulence=false,
                T_solidus_val=SVector{3,Float64}([1400.0, 1400.0, NaN]),
                T_liquidus_val=SVector{3,Float64}([1800.0, 1800.0, NaN]),
            )
        end

        # Case 2: soft_turbulence = true
        for m in 1:marknum
            compute_marker_properties!(
                m,
                tm,
                tkm,
                rhototalm,
                rhocptotalm,
                etatotalm,
                hrtotalm,
                ktotalm_on,
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
                soft_turbulence=true,
                turb_exponent_val=1.0 / 3.0,
                T_solidus_val=SVector{3,Float64}([1400.0, 1400.0, NaN]),
                T_liquidus_val=SVector{3,Float64}([1800.0, 1800.0, NaN]),
                eta_fluid_silicate_val=100.0,
                F_turb_start_val=0.30,
                F_turb_end_val=0.50,
                dT_turb_min_val=10.0,
                T_surface_ref_val=300.0,
                k_turb_cutoff_val=1.0e6,
                k_turb_floor_val=1.0e-3,
            )
        end

        # Marker 2: Sub-solidus -> conductivity identical in both cases
        @test isapprox(ktotalm_on[2], ktotalm_off[2]; rtol=1e-12)
        @test ktotalm_on[2] < 10.0

        # Marker 3: Super-liquidus -> convective conductivity enhanced by orders of magnitude!
        @test ktotalm_off[3] < 10.0
        @test ktotalm_on[3] > 5.0e3
        @test ktotalm_on[3] / ktotalm_off[3] > 1.5e3
        @test ktotalm_on[3] <= 1.0e6

        # Marker 1: Partial melt (~0.45) -> smoothly intermediate between base and fully turbulent
        @test ktotalm_on[1] > ktotalm_off[1]
        @test ktotalm_on[1] < ktotalm_on[3]

        # Marker 4: Sticky air -> unaffected
        @test isapprox(ktotalm_on[4], ktotalm_off[4]; rtol=1e-12)

        # Marker Monotonicity across melt range [0.25, 0.625] with real marker path
        # Verifies no dip occurs at melting onset (F_m ~ 0.30 - 0.35)
        T_sweep = 1500.0:20.0:1650.0
        k_sweep = Float64[]
        for T_cur in T_sweep
            tkm[1] = T_cur
            compute_marker_properties!(
                1,
                tm,
                tkm,
                rhototalm,
                rhocptotalm,
                etatotalm,
                hrtotalm,
                ktotalm_on,
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
                soft_turbulence=true,
                turb_exponent_val=1.0 / 3.0,
                T_solidus_val=SVector{3,Float64}([1400.0, 1400.0, NaN]),
                T_liquidus_val=SVector{3,Float64}([1800.0, 1800.0, NaN]),
            )
            push!(k_sweep, ktotalm_on[1])
        end
        @test all(k_sweep .>= ktotalm_off[2])
        for idx in 1:(length(k_sweep) - 1)
            @test k_sweep[idx + 1] >= k_sweep[idx]
        end
    end

    @testset "Grid Interpolation: KX & KY receive enhanced conductivity" begin
        coords = default_grid_coordinates()
        Nx1 = coords.Nx1
        Ny1 = coords.Ny1
        KX = zeros(Float64, Ny1, Nx1)
        KY = zeros(Float64, Ny1, Nx1)
        KXSUM = zeros(Float64, Ny1, Nx1)
        KYSUM = zeros(Float64, Ny1, Nx1)
        WTXSUM = zeros(Float64, Ny1, Nx1)
        WTY_SUM = zeros(Float64, Ny1, Nx1)
        RHOXSUM = zeros(Float64, Ny1, Nx1)
        RHOFXSUM = zeros(Float64, Ny1, Nx1)
        PHIXSUM = zeros(Float64, Ny1, Nx1)
        RXSUM = zeros(Float64, Ny1, Nx1)
        RHOYSUM = zeros(Float64, Ny1, Nx1)
        RHOFYSUM = zeros(Float64, Ny1, Nx1)
        PHIYSUM = zeros(Float64, Ny1, Nx1)
        RY_SUM = zeros(Float64, Ny1, Nx1)

        # Place 4 markers in the center cell with high conductivity
        xm = [coords.xp[3], coords.xp[3], coords.xp[3], coords.xp[3]]
        ym = [coords.yp[3], coords.yp[3], coords.yp[3], coords.yp[3]]
        ktotalm = [2.5e5, 2.5e5, 2.5e5, 2.5e5]
        rhototalm = [3000.0, 3000.0, 3000.0, 3000.0]
        rhofluidcur = [1000.0, 1000.0, 1000.0, 1000.0]
        phim = [0.1, 0.1, 0.1, 0.1]
        etafluidcur_inv_kphim = [1.0, 1.0, 1.0, 1.0]

        for m in 1:4
            marker_to_vx_nodes!(
                m,
                xm[m],
                ym[m],
                rhototalm,
                rhofluidcur,
                ktotalm,
                phim,
                etafluidcur_inv_kphim,
                RHOXSUM,
                RHOFXSUM,
                KXSUM,
                PHIXSUM,
                RXSUM,
                WTXSUM;
                coords=coords,
            )
            marker_to_vy_nodes!(
                m,
                xm[m],
                ym[m],
                rhototalm,
                rhofluidcur,
                ktotalm,
                phim,
                etafluidcur_inv_kphim,
                RHOYSUM,
                RHOFYSUM,
                KYSUM,
                PHIYSUM,
                RY_SUM,
                WTY_SUM;
                coords=coords,
            )
        end
        for j in 1:Nx1, i in 1:Ny1
            if WTXSUM[i, j] > 0.0
                KX[i, j] = KXSUM[i, j] / WTXSUM[i, j]
            end
            if WTY_SUM[i, j] > 0.0
                KY[i, j] = KYSUM[i, j] / WTY_SUM[i, j]
            end
        end

        # Center node has enhanced conductivity
        @test KX[3, 3] > 1.0e4
        @test KY[3, 3] > 1.0e4
    end

    @testset "Mini-Simulation Execution with Soft Turbulence" begin
        output_dir = mktempdir()
        try
            cfg_test = SimulationConfig(;
                grid=GridConfig(; Nx=17, Ny=17, xsize=140000.0, ysize=140000.0),
                time=TimeConfig(; dt_initial=100.0, dt_longest=100.0, n_steps=2),
                materials=MaterialConfig(;
                    tkm0=SVector{3,Float64}([1650.0, 1650.0, 300.0])
                ),
                melting=MeltingConfig(;
                    active=true,
                    soft_turbulence=true,
                    turb_exponent=1.0 / 3.0,
                    eta_fluid_silicate=100.0,
                    F_turb_start=0.30,
                    F_turb_end=0.50,
                    T_surface_ref=300.0,
                ),
                output=OutputConfig(; output_dir="output_soft_turb_test", savematstep=2),
            )
            @test validate_config(cfg_test) === nothing
            Erebus.simulation_loop(cfg_test; output_path=output_dir)
            files = readdir(output_dir)
            @test "output_00000.jld2" in files
            @test "output_00002.jld2" in files
            data2 = JLD2.load(joinpath(output_dir, "output_00002.jld2"))
            @test data2["timestep"] == 2
            @test !any(isnan, data2["tk2"])
            @test !any(isnan, data2["pr"])
            # Verify that hot silicate markers have enhanced turbulent conductivity
            @test any(data2["ktotalm"] .> 1.0e3)
            # Verify that air markers remain at normal conductivities
            air_idx = findall(data2["tm"] .== 3)
            @test !isempty(air_idx)
            @test all(data2["ktotalm"][air_idx] .<= 3000.0)
        finally
            rm(output_dir; recursive=true, force=true)
        end
    end
end

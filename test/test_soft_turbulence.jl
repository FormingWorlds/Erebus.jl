using Test
using Erebus
using Erebus.Config
using Erebus.Physics
using Erebus.Particles
using Erebus.Geometry
using StaticArrays
using TOML

@testset "Sub-Grid Soft Turbulence & Regularized Conductivity" begin
    @testset "Config Schema & Validation" begin
        # Default configuration has soft turbulence disabled
        melt_def = MeltingConfig()
        @test melt_def.soft_turbulence == false
        @test isapprox(melt_def.eta_fluid_silicate, 100.0; rtol=1e-12)
        @test isapprox(melt_def.F_turb_start, 0.30; rtol=1e-12)
        @test isapprox(melt_def.F_turb_end, 0.50; rtol=1e-12)
        @test isapprox(melt_def.F_turb_crit, 0.40; rtol=1e-12)
        @test isapprox(melt_def.dT_turb_min, 10.0; rtol=1e-12)
        @test isapprox(melt_def.T_surface_ref, 300.0; rtol=1e-12)
        @test isapprox(melt_def.k_turb_cutoff, 1.0e6; rtol=1e-12)
        @test isapprox(melt_def.k_turb_floor, 1.0e-3; rtol=1e-12)

        # Valid active soft turbulence configuration
        cfg_valid = SimulationConfig(;
            melting=MeltingConfig(;
                active=true,
                soft_turbulence=true,
                eta_fluid_silicate=100.0,
                F_turb_start=0.30,
                F_turb_end=0.50,
                F_turb_crit=0.40,
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
                melting=MeltingConfig(;
                    active=true, soft_turbulence=true, F_turb_end=1.2
                ),
            ),
        )

        # Invalid: F_turb_crit outside [F_turb_start, F_turb_end]
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(;
                    active=true, soft_turbulence=true, F_turb_start=0.3, F_turb_end=0.5, F_turb_crit=0.25
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(;
                    active=true, soft_turbulence=true, F_turb_start=0.3, F_turb_end=0.5, F_turb_crit=0.55
                ),
            ),
        )

        # Invalid: dT_turb_min <= 0 or non-finite
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(;
                    active=true, soft_turbulence=true, dT_turb_min=0.0
                ),
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
                eta_fluid_silicate=50.0,
                F_turb_start=0.25,
                F_turb_end=0.55,
                F_turb_crit=0.38,
                dT_turb_min=15.0,
                T_surface_ref=280.0,
                k_turb_cutoff=5.0e5,
                k_turb_floor=5.0e-4,
            ),
        )
        toml_str = save_config(cfg_custom)
        cfg_reloaded = load_config(toml_str)
        @test cfg_reloaded.melting.soft_turbulence == true
        @test isapprox(cfg_reloaded.melting.eta_fluid_silicate, 50.0; rtol=1e-12)
        @test isapprox(cfg_reloaded.melting.F_turb_start, 0.25; rtol=1e-12)
        @test isapprox(cfg_reloaded.melting.F_turb_end, 0.55; rtol=1e-12)
        @test isapprox(cfg_reloaded.melting.F_turb_crit, 0.38; rtol=1e-12)
        @test isapprox(cfg_reloaded.melting.dT_turb_min, 15.0; rtol=1e-12)
        @test isapprox(cfg_reloaded.melting.T_surface_ref, 280.0; rtol=1e-12)
        @test isapprox(cfg_reloaded.melting.k_turb_cutoff, 5.0e5; rtol=1e-12)
        @test isapprox(cfg_reloaded.melting.k_turb_floor, 5.0e-4; rtol=1e-12)
    end

    @testset "Physics: regularized_soft_turbulence_conductivity" begin
        k_cond = 3.0
        eta_num = 1.0e12
        eta_fluid = 100.0
        # Turbulent conductivity target: k_turb = 3.0 * sqrt(1e12 / 100) = 3.0 * 1e5 = 3.0e5 W/(m K)
        k_turb_expected = 3.0 * sqrt(1.0e12 / 100.0)
        @test isapprox(k_turb_expected, 3.0e5; rtol=1e-12)

        # Domain errors
        @test_throws DomainError regularized_soft_turbulence_conductivity(
            0.0, eta_num, eta_fluid, 0.5, 1500.0, 300.0
        )
        @test_throws DomainError regularized_soft_turbulence_conductivity(
            -3.0, eta_num, eta_fluid, 0.5, 1500.0, 300.0
        )
        @test_throws DomainError regularized_soft_turbulence_conductivity(
            k_cond, 0.0, eta_fluid, 0.5, 1500.0, 300.0
        )
        @test_throws DomainError regularized_soft_turbulence_conductivity(
            k_cond, eta_num, -10.0, 0.5, 1500.0, 300.0
        )

        # Asymptotic limit 1: Below F_turb_start (solid / sub-threshold rock)
        # Should return k_cond exactly
        k_sub = regularized_soft_turbulence_conductivity(
            k_cond, eta_num, eta_fluid, 0.20, 1500.0, 300.0;
            F_start=0.30, F_end=0.50
        )
        @test isapprox(k_sub, k_cond; rtol=1e-12)

        # Asymptotic limit 2: Zero temperature contrast (isothermal)
        # Should return k_cond exactly
        k_iso = regularized_soft_turbulence_conductivity(
            k_cond, eta_num, eta_fluid, 0.80, 300.0, 300.0;
            F_start=0.30, F_end=0.50, dT_min=10.0
        )
        @test isapprox(k_iso, k_cond; rtol=1e-12)

        # Asymptotic limit 3: Fully molten (F_m >= F_turb_end) and large dT
        # Weight w = 1.0, so should return k_turb_expected exactly
        k_full = regularized_soft_turbulence_conductivity(
            k_cond, eta_num, eta_fluid, 0.60, 1500.0, 300.0;
            F_start=0.30, F_end=0.50, dT_min=10.0
        )
        @test isapprox(k_full, k_turb_expected; rtol=1e-12)

        # Smoothness and Monotonicity across transition window [0.30, 0.50]
        F_vals = 0.28:0.02:0.52
        k_vals = [
            regularized_soft_turbulence_conductivity(
                k_cond, eta_num, eta_fluid, F, 1500.0, 300.0;
                F_start=0.30, F_end=0.50, dT_min=10.0
            ) for F in F_vals
        ]

        # Strictly monotonic non-decreasing
        for idx in 1:(length(k_vals) - 1)
            @test k_vals[idx + 1] >= k_vals[idx]
        end

        # Value at midpoint F = 0.40 (smoothstep xi = 0.5, w = 0.5)
        # Logarithmic midpoint: log10(k_mid) = 0.5 * log10(3) + 0.5 * log10(3e5) = log10(sqrt(3 * 3e5))
        k_mid = regularized_soft_turbulence_conductivity(
            k_cond, eta_num, eta_fluid, 0.40, 1500.0, 300.0;
            F_start=0.30, F_end=0.50, dT_min=10.0
        )
        @test isapprox(k_mid, sqrt(k_cond * k_turb_expected); rtol=1e-10)

        # Clamping at k_cutoff
        # For liquid iron: eta_fluid = 0.01, k_turb = 3 * sqrt(1e12 / 0.01) = 3e7 > 1e6
        k_clamped = regularized_soft_turbulence_conductivity(
            k_cond, 1.0e12, 0.01, 1.0, 1500.0, 300.0;
            k_cutoff=1.0e6
        )
        @test isapprox(k_clamped, 1.0e6; rtol=1e-12)

        # Clamping at k_floor
        k_floored = regularized_soft_turbulence_conductivity(
            1.0e-5, 1.0, 1.0, 0.0, 1500.0, 300.0;
            k_floor=1.0e-3
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
                m, tm, tkm, rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm_off,
                tkm_rhocptotalm, etafluidcur_inv_kphim, hrsolidm, hrfluidm, phim,
                XWsolidm0, mode, rhofluidcur;
                pm=pfm0, Fm=Fm, melting_active=true, soft_turbulence=false,
                T_solidus_val=SVector{3,Float64}([1400.0, 1400.0, NaN]),
                T_liquidus_val=SVector{3,Float64}([1800.0, 1800.0, NaN]),
            )
        end

        # Case 2: soft_turbulence = true
        for m in 1:marknum
            compute_marker_properties!(
                m, tm, tkm, rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm_on,
                tkm_rhocptotalm, etafluidcur_inv_kphim, hrsolidm, hrfluidm, phim,
                XWsolidm0, mode, rhofluidcur;
                pm=pfm0, Fm=Fm, melting_active=true, soft_turbulence=true,
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
        @test ktotalm_on[3] > 1.0e4
        @test ktotalm_on[3] <= 1.0e6

        # Marker 1: Partial melt (~0.45) -> smoothly intermediate between base and fully turbulent
        @test ktotalm_on[1] > ktotalm_off[1]
        @test ktotalm_on[1] < ktotalm_on[3]

        # Marker 4: Sticky air -> unaffected
        @test isapprox(ktotalm_on[4], ktotalm_off[4]; rtol=1e-12)
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
                m, xm[m], ym[m], rhototalm, rhofluidcur, ktotalm, phim,
                etafluidcur_inv_kphim, RHOXSUM, RHOFXSUM, KXSUM, PHIXSUM,
                RXSUM, WTXSUM; coords=coords,
            )
            marker_to_vy_nodes!(
                m, xm[m], ym[m], rhototalm, rhofluidcur, ktotalm, phim,
                etafluidcur_inv_kphim, RHOYSUM, RHOFYSUM, KYSUM, PHIYSUM,
                RY_SUM, WTY_SUM; coords=coords,
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
        cfg_test = SimulationConfig(;
            time=TimeConfig(; dt_initial=100.0, n_steps=2),
            melting=MeltingConfig(;
                active=true,
                soft_turbulence=true,
                eta_fluid_silicate=100.0,
                F_turb_start=0.30,
                F_turb_end=0.50,
                T_surface_ref=300.0,
            ),
        )
        @test validate_config(cfg_test) === nothing
        (ts_c, dt_c, time_c, mark_c, _, _, _) = Erebus.setup_dynamic_simulation_parameters(cfg_test)
        @test ts_c >= 0
        @test dt_c > 0.0
    end
end

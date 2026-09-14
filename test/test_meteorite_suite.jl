using Test
using Erebus

@testset "Synthetic Meteorite Suite Physics and Classification" begin
    # 1. Two-point cooling rate analytical pins
    t1 = 0.0
    t2 = 3.15576e13 # exactly 1 Myr
    T1 = 800.0
    T2 = 700.0

    cr_myr = Erebus.compute_cooling_rate(t1, T1, t2, T2; units=:K_per_Myr)
    @test isapprox(cr_myr, 100.0; atol=1e-10)

    cr_yr = Erebus.compute_cooling_rate(t1, T1, t2, T2; units=:K_per_yr)
    @test isapprox(cr_yr, 1.0e-4; atol=1e-12)
    @test isapprox(cr_myr, cr_yr * 1.0e6; atol=1e-10)

    # Error contracts
    @test_throws DomainError Erebus.compute_cooling_rate(t2, T1, t1, T2) # t2 <= t1
    @test_throws ArgumentError Erebus.compute_cooling_rate(t1, T1, t2, T2; units=:INVALID)

    # 2. Time-series cooling rate at closure temperature (773.15 K)
    sec_myr = 3.15576e13
    times = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0] .* sec_myr
    # Temperature trajectory: heating to 1000 K at 2 Myr, cooling to 500 K at 5 Myr
    # Between 3 Myr (850 K) and 4 Myr (700 K), passes through 773.15 K
    # Cooling rate in this interval: (850 - 700) / 1 Myr = 150 K/Myr
    temps = [300.0, 700.0, 1000.0, 850.0, 700.0, 500.0]

    cr_series = Erebus.compute_cooling_rate(times, temps, 773.15; units=:K_per_Myr)
    @test isapprox(cr_series, 150.0; atol=1e-10)

    # Case: body never reached closure temperature (e.g. cold body, peak T = 600 K)
    temps_cold = [300.0, 450.0, 600.0, 500.0, 400.0, 300.0]
    @test isnan(Erebus.compute_cooling_rate(times, temps_cold, 773.15))

    # Case: body still hotter than closure temperature at simulation end
    temps_hot = [300.0, 700.0, 1100.0, 1050.0, 950.0, 850.0]
    @test isnan(Erebus.compute_cooling_rate(times, temps_hot, 773.15))

    # Error contracts for series
    @test_throws DimensionMismatch Erebus.compute_cooling_rate([0.0, 1.0], [500.0])
    @test_throws ArgumentError Erebus.compute_cooling_rate([1.0], [500.0])
    @test_throws DomainError Erebus.compute_cooling_rate(times, temps, -10.0)

    # 3. Analytical conductive cooling rate scaling
    R_50km = 50_000.0
    R_100km = 100_000.0
    d_center_50 = 50_000.0
    d_center_100 = 100_000.0

    cr_cond_50 = Erebus.compute_conductive_cooling_rate(R_50km, d_center_50)
    cr_cond_100 = Erebus.compute_conductive_cooling_rate(R_100km, d_center_100)

    # R^2 scaling guard: a 100 km body cools 4x slower than a 50 km body
    @test isapprox(cr_cond_50 / cr_cond_100, 4.0; atol=1e-3)

    # Depth scaling guard: near-surface shell cools faster than core
    d_shallow_50 = 5_000.0
    cr_shallow_50 = Erebus.compute_conductive_cooling_rate(R_50km, d_shallow_50)
    @test cr_shallow_50 > cr_cond_50

    # Diffusivity scaling guard: 2x kappa -> 2x cooling rate
    cr_2x_kap = Erebus.compute_conductive_cooling_rate(R_50km, d_center_50; kappa=2.0e-6)
    @test isapprox(cr_2x_kap / cr_cond_50, 2.0; atol=1e-10)

    # Error contracts
    @test_throws DomainError Erebus.compute_conductive_cooling_rate(-50_000.0, 10_000.0)
    @test_throws DomainError Erebus.compute_conductive_cooling_rate(50_000.0, 60_000.0)
    @test_throws DomainError Erebus.compute_conductive_cooling_rate(50_000.0, -100.0)
    @test_throws DomainError Erebus.compute_conductive_cooling_rate(
        50_000.0, 10_000.0; kappa=-1e-6
    )
    @test_throws DomainError Erebus.compute_conductive_cooling_rate(
        50_000.0, 10_000.0; delta_T=-100.0
    )
    @test_throws ArgumentError Erebus.compute_conductive_cooling_rate(
        50_000.0, 10_000.0; units=:INVALID
    )

    # 4. Metallographic cooling rate proxy (Yang & Goldstein 2006)
    proxy_fast = Erebus.compute_metallographic_cooling_rate_proxy(500.0) # 500 K/Myr
    proxy_slow = Erebus.compute_metallographic_cooling_rate_proxy(1.0)   # 1 K/Myr

    @test isapprox(proxy_fast.cooling_rate_metallographic_K_per_Myr, 500.0; atol=1e-12)
    @test isapprox(proxy_slow.cooling_rate_metallographic_K_per_Myr, 1.0; atol=1e-12)
    # Slow cooling yields higher central Ni due to longer diffusion / kamacite growth
    @test proxy_slow.taenite_central_Ni_wtpct > proxy_fast.taenite_central_Ni_wtpct
    @test 10.0 <= proxy_fast.taenite_central_Ni_wtpct <= 45.0
    @test 10.0 <= proxy_slow.taenite_central_Ni_wtpct <= 45.0

    # Pinned calibration points: 1 K/Myr -> 25 wt% Ni, 10 K/Myr -> 20 wt% Ni, 100 K/Myr -> 15 wt% Ni
    @test isapprox(proxy_slow.taenite_central_Ni_wtpct, 25.0; atol=1e-5)
    proxy_10 = Erebus.compute_metallographic_cooling_rate_proxy(10.0)
    @test isapprox(proxy_10.taenite_central_Ni_wtpct, 20.0; atol=1e-5)
    proxy_100 = Erebus.compute_metallographic_cooling_rate_proxy(100.0)
    @test isapprox(proxy_100.taenite_central_Ni_wtpct, 15.0; atol=1e-5)

    # NaN / non-positive cooling rate handling
    proxy_nan = Erebus.compute_metallographic_cooling_rate_proxy(NaN)
    @test isnan(proxy_nan.cooling_rate_metallographic_K_per_Myr)
    @test isnan(proxy_nan.taenite_central_Ni_wtpct)
    proxy_zero = Erebus.compute_metallographic_cooling_rate_proxy(0.0)
    @test isnan(proxy_zero.cooling_rate_metallographic_K_per_Myr)

    # 5. Petrologic type classification (Van Schmus & Wood 1967; Huss et al. 2006)
    # Achondrites
    @test Erebus.classify_petrologic_type(1600.0, 0.70) === :achondrite
    @test Erebus.classify_petrologic_type(1400.0, 0.25) === :primitive_achondrite

    # Chondrites: Types 1-7
    @test Erebus.classify_petrologic_type(350.0, 0.0; X_water=0.08) === Symbol("Type 1")
    @test Erebus.classify_petrologic_type(480.0, 0.0; X_water=0.03) === Symbol("Type 2")
    @test Erebus.classify_petrologic_type(700.0, 0.0; X_water=0.001) === Symbol("Type 3")
    @test Erebus.classify_petrologic_type(920.0, 0.0) === Symbol("Type 4")
    @test Erebus.classify_petrologic_type(1020.0, 0.0) === Symbol("Type 5")
    @test Erebus.classify_petrologic_type(1150.0, 0.0) === Symbol("Type 6")
    @test Erebus.classify_petrologic_type(1250.0, 0.05) === Symbol("Type 7")

    # Error contracts
    @test_throws DomainError Erebus.classify_petrologic_type(-10.0, 0.0)
    @test_throws DomainError Erebus.classify_petrologic_type(500.0, -0.1)
    @test_throws DomainError Erebus.classify_petrologic_type(500.0, 1.2)
    @test_throws DomainError Erebus.classify_petrologic_type(500.0, 0.0; X_water=-0.05)
    @test_throws DomainError Erebus.classify_petrologic_type(500.0, 0.0; X_water=1.5)

    # 6. Full synthetic meteorite suite generation
    R_p = 50_000.0
    r_arr = [45_000.0, 35_000.0, 20_000.0, 5_000.0]
    T_peak_arr = [400.0, 750.0, 1100.0, 1600.0]
    F_melt_arr = [0.0, 0.0, 0.0, 0.85]
    cr_arr = [100.0, 50.0, 25.0, 5.0]
    x_w_arr = [0.07, 0.002, 0.0, 0.0]
    x_c_arr = [0.03, 0.02, 0.005, 0.0]

    suite = Erebus.generate_synthetic_meteorite_suite(
        r_arr, T_peak_arr, F_melt_arr, cr_arr, R_p; X_water_m=x_w_arr, X_refr_C_m=x_c_arr
    )

    @test length(suite) == 4
    # Sample 1: Shallow shell (depth = 5 km), aqueously altered Type 1
    @test isapprox(suite[1].depth_m, 5_000.0; atol=1e-12)
    @test isapprox(suite[1].radius_m, 45_000.0; atol=1e-12)
    @test suite[1].petrologic_type === Symbol("Type 1")
    @test isapprox(suite[1].cooling_rate_773K_K_per_Myr, 100.0; atol=1e-12)
    @test isapprox(suite[1].X_water_final, 0.07; atol=1e-12)

    # Sample 2: Intermediate shell (depth = 15 km), pristine Type 3
    @test isapprox(suite[2].depth_m, 15_000.0; atol=1e-12)
    @test suite[2].petrologic_type === Symbol("Type 3")

    # Sample 3: Deep mantle (depth = 30 km), metamorphosed Type 6
    @test isapprox(suite[3].depth_m, 30_000.0; atol=1e-12)
    @test suite[3].petrologic_type === Symbol("Type 6")

    # Sample 4: Core / deep interior (depth = 45 km), differentiated Achondrite
    @test isapprox(suite[4].depth_m, 45_000.0; atol=1e-12)
    @test suite[4].petrologic_type === :achondrite
    @test isapprox(suite[4].F_melt_peak, 0.85; atol=1e-12)

    # Error contracts
    @test_throws DimensionMismatch Erebus.generate_synthetic_meteorite_suite(
        r_arr[1:3], T_peak_arr, F_melt_arr, cr_arr, R_p
    )
    @test_throws DomainError Erebus.generate_synthetic_meteorite_suite(
        r_arr, T_peak_arr, F_melt_arr, cr_arr, -50_000.0
    )
end

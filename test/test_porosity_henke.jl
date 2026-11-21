using Erebus
using Test

# Independent reference for the cold pressing and hot sintering laws of Henke et al.
# (2012), coded from the raw constants of the paper and of the i2elvis code
# (markmart.c), in BigFloat. These helpers do not call the Erebus functions.

const _HK_SECONDS_PER_YEAR = 3.15576e7

"""
Reference grain radius `R` and contact radius `a` [m] for porosity `phi` after
sintering from `phi0`, by bisection of the volume balance of Henke et al. (2012)
Eqs. 34 to 37 with 8 contacts, in BigFloat.
"""
function _hk_grain_root(phi, phi0, R0; iterations::Int=140)
    R0b = big(R0)
    X = (1 - big(phi0)) / (1 - big(phi))
    X = max(X, big(1) / 3 + big(1.0e-12))
    c = cbrt(X)
    # Eq. 36 with Eqs. 34, 35 and 37: R0^3 = R^3 - Z/8 (R - R0 c)^2 (2 R + R0 c)
    balance(R) = R^3 - (R - R0b * c)^2 * (2R + R0b * c) - R0b^3
    lo, hi = R0b, 2R0b * c
    for _ in 1:iterations
        mid = (lo + hi) / 2
        if balance(mid) < 0
            lo = mid
        else
            hi = mid
        end
    end
    R = (lo + hi) / 2
    a = sqrt(max(R^2 - R0b^2 * c^2, big(0)))
    return R, a
end

"""
Reference porosity decrease rate [1/s] (Henke et al. 2012 Eqs. 39 to 44, with the
negative exponent and the i2elvis units: P in bar, R in cm, R_gas in cal/(K mol)).
"""
function _hk_rate(phi, phi0, P_bar, T; R0=1.0e-6, A=4.0e-5, E=85000.0, n=1.5)
    iszero(phi) && return big(0)
    R, a = _hk_grain_root(phi, phi0, R0)
    sigma = 2 * sqrt(big(3)) * (R^2 - a^2) * big(P_bar) / (big(pi) * a^2)
    F = phi <= 0.1 ? 10 * big(phi) : big(1)
    return F * big(A) * (1 - big(phi)) * sigma^big(n) / (100R)^3 *
           exp(-big(E) / (big(1.9872) * big(T)))
end

"""
Fine RK4 integration of dphi/dt = -rate(phi) in BigFloat from `phi_start` over `t_end`.
"""
function _hk_rk4(phi_start, phi0, P_bar, T, t_end; nsteps::Int=400)
    h = big(t_end) / nsteps
    phi = big(phi_start)
    f(p) = -_hk_rate(p, phi0, P_bar, T)
    for _ in 1:nsteps
        k1 = f(phi)
        k2 = f(phi + h * k1 / 2)
        k3 = f(phi + h * k2 / 2)
        k4 = f(phi + h * k3)
        phi += h * (k1 + 2k2 + 2k3 + k4) / 6
    end
    return phi
end

"""
Time [s] to sinter from `phi0` down to `phi0 - delta`, from the quadrature of
dt = dphi / rate with the substitution delta' = delta u^2 (smooth integrand).
"""
function _hk_onset_time(delta, phi0, P_bar, T; nodes::Int=400)
    d = big(delta)
    h = big(1) / nodes
    integrand(u) = iszero(u) ? big(0) : 2d * u / _hk_rate(phi0 - d * u^2, phi0, P_bar, T)
    s = integrand(big(0)) + integrand(big(1))
    for k in 1:(nodes - 1)
        s += (isodd(k) ? 4 : 2) * integrand(k * h)
    end
    return s * h / 3
end

# Output of tools/replay_i2elvis_porosity.jl: the porosity block of i2elvis markmart.c
# (commit e2c8487), compiled with cc and run once per state. Columns: phi, T [K],
# P [Pa], phi_init, dt [s], phi after 1 step, porosity change of the step.
const _HK_I2ELVIS_REPLAY = (
    (0.8, 500.0, 50000.0, 0.8, 1.0e10, 0.461274653673172, 0.0),
    (0.6, 650.0, 1.0e6, 0.75, 1.0e10, 0.4202621281147003, 0.0),
    (0.84, 699.5, 1000.0, 0.95, 1.0e10, 0.8399999737739563, 0.0),
    (0.38, 900.0, 100000.0, 0.4, 1.9e7, 0.3796263635158539, 0.000373630870498235),
    (0.08, 1000.0, 500000.0, 0.4, 200000.0, 0.07992299646139145, 7.700542985652849e-5),
    (0.3, 1100.0, 2.0e6, 0.5, 10000.0, 0.2931424379348755, 0.006857580976981374),
    (0.84, 800.0, 1000.0, 0.95, 1.8e14, 0.8391586542129517, 0.0008413204775132069),
    (0.4, 700.0, 200000.0, 0.42, 6.0e12, 0.3996112644672394, 0.0003887338021696076),
)

@testset "Porosity laws of Henke et al. (2012)" begin
    setprecision(BigFloat, 128) do
        @testset "Cold pressing law phi_H(P), Henke Eq. 15 and L16 Eq. 10" begin
            # Pinned by hand: 0.42 + 0.46 = 0.88 at P = 0 and 0.42 + 0.46 / 2 at P = P0
            @test isapprox(Erebus.phi_cold(0.0), 0.88; rtol=1e-14)
            @test isapprox(Erebus.phi_cold(0.13), 0.65; rtol=1e-14)
            # 1.3 bar = 10 P0: 0.42 + 0.46 / (10^1.72 + 1) = 0.428601...
            phi_13 = Erebus.phi_cold(1.3)
            @test isapprox(phi_13, 0.4286013; atol=1e-7)
            # Exponent guard: 1.7 in place of 1.72 moves phi_H(1.3 bar) by 4e-4
            @test !isapprox(phi_13, 0.42 + 0.46 / (exp10(1.7) + 1.0); atol=3.0e-4)
            # Scale guard: 1.3 Pa = 1.3e-5 bar is still uncompacted powder
            @test isapprox(Erebus.phi_cold(1.3e-5), 0.88; atol=1e-6)
            # Limits and monotonic decrease with pressure
            @test isapprox(Erebus.phi_cold(1.0e6), 0.42; atol=1e-9)
            Ps = [0.0, 0.01, 0.13, 1.0, 10.0, 100.0]
            vals = Erebus.phi_cold.(Ps)
            @test all(diff(vals) .< 0.0)
            @test all(0.42 .<= vals .<= 0.88)
            @test_throws DomainError Erebus.phi_cold(-1.0e-3)
            @test_throws DomainError Erebus.phi_cold(NaN)
        end

        @testset "Starting porosity phi_0 = min(phi_init, phi_H(P), phi_cap)" begin
            # P = 0: raw 0.88 is capped at 0.85
            @test isapprox(Erebus.porosity_phi0(0.0, 0.95), 0.85; rtol=1e-14)
            @test isapprox(Erebus.porosity_phi0(0.0, 0.95; phi_cap=0.9), 0.88; rtol=1e-14)
            # phi_init below the pressing curve wins
            @test isapprox(Erebus.porosity_phi0(0.13, 0.30), 0.30; rtol=1e-14)
            # pressing curve below phi_init wins
            @test isapprox(Erebus.porosity_phi0(0.13, 0.75), 0.65; rtol=1e-14)
            # Edge case: phi_init = 0 gives a dense rock at every pressure
            @test iszero(Erebus.porosity_phi0(50.0, 0.0))
            @test_throws DomainError Erebus.porosity_phi0(1.0, 1.2)
            @test_throws DomainError Erebus.porosity_phi0(1.0, -0.1)
            @test_throws DomainError Erebus.porosity_phi0(1.0, 0.5; phi_cap=1.5)
            @test_throws DomainError Erebus.porosity_phi0(-1.0, 0.5)
        end

        @testset "Grain geometry root, Henke Eqs. 34 to 37 (8 contacts)" begin
            R0 = 1.0e-6
            for (phi, phi0) in ((0.40, 0.42), (0.30, 0.65), (0.10, 0.42), (0.419999, 0.42))
                R, a = Erebus.grain_geometry(phi, phi0, R0)
                R_ref, a_ref = _hk_grain_root(phi, phi0, R0)
                @test isapprox(R, Float64(R_ref); rtol=1e-12)
                @test isapprox(a, Float64(a_ref); rtol=1e-8)
                # Volume conservation, Eq. 36 with the cap volume of Eqs. 34 and 35
                h = R - sqrt(R^2 - a^2)
                V_cap = pi / 6 * h^2 * (3R - h)
                @test isapprox(4pi / 3 * R^3 - 8 * V_cap, 4pi / 3 * R0^3; rtol=1e-12)
                # The physical root: R in [R0, 2 R0 X^(1/3)], a < R
                X = (1 - phi0) / (1 - phi)
                @test R0 <= R <= 2R0 * cbrt(X)
                @test 0.0 < a < R
            end
            # No sintering: X = 1 gives R = R0 and a = 0
            R1, a1 = Erebus.grain_geometry(0.42, 0.42, R0)
            @test isapprox(R1, R0; rtol=1e-14)
            @test isapprox(a1, 0.0; atol=1e-14 * R0)
            # Spurious root guard: the second root of the cubic at X -> 1 is 2.73 R0
            @test R1 < 1.01R0
            # Contact-overlap limit X <= 1/3: geometry frozen at X = 1/3
            Rlim, alim = Erebus.grain_geometry(0.0, 0.80, R0)
            Rfreeze, afreeze = Erebus.grain_geometry(0.40, 0.80, R0)
            @test isapprox(Rlim, Rfreeze; rtol=1e-6)
            @test isapprox(Rlim, 2R0 * cbrt(1 / 3); rtol=1e-5)
            @test isapprox(alim, afreeze; rtol=1e-6)
            @test_throws DomainError Erebus.grain_geometry(0.5, 0.4, R0)
            @test_throws DomainError Erebus.grain_geometry(0.3, 0.4, 0.0)
            @test_throws DomainError Erebus.grain_geometry(-0.1, 0.4, R0)
            @test_throws DomainError Erebus.grain_geometry(0.3, 1.0, R0)
        end

        @testset "Contact stress, Henke Eqs. 40 and 41, by hand at 1 state" begin
            R0 = 1.0e-6
            P = 1.0
            R, a = Erebus.grain_geometry(0.40, 0.42, R0)
            sigma = Erebus.contact_stress(P, R, a)
            # By hand: C_av = 2 sqrt(3) (R^2 - a^2), sigma = C_av P / (pi a^2)
            R_ref, a_ref = _hk_grain_root(0.40, 0.42, R0)
            sigma_ref =
                2 * sqrt(3.0) * Float64(R_ref^2 - a_ref^2) * P / (pi * Float64(a_ref^2))
            @test isapprox(sigma, sigma_ref; rtol=1e-8)
            # sigma >> P early in sintering; 6 contacts would change C_av by 4/3 only in
            # the geometry, so the scale guard is on the stress amplification itself
            @test 10.0 < sigma / P < 100.0
            # Linear in P, zero at P = 0
            @test isapprox(Erebus.contact_stress(2P, R, a), 2sigma; rtol=1e-14)
            @test iszero(Erebus.contact_stress(0.0, R, a))
            @test_throws DomainError Erebus.contact_stress(-1.0, R, a)
            @test_throws DomainError Erebus.contact_stress(P, R, 0.0)
            @test_throws DomainError Erebus.contact_stress(P, R, 2R)
        end

        @testset "Explicit sintering rate, Henke Eqs. 39, 43, 44" begin
            R0 = 1.0e-6
            for (phi, phi0, P, T) in ((0.38, 0.40, 1.0, 900.0), (0.08, 0.40, 5.0, 1000.0))
                r = Erebus.porosity_decrease_rate(phi, phi0, P, T)
                r_ref = Float64(_hk_rate(phi, phi0, P, T))
                @test isapprox(r, r_ref; rtol=1e-8)
                @test 0.0 < r < 1.0
            end
            # Scale guard: the 0.01 decrease time at 900 K, 1 bar, phi 0.38 is ~16 yr
            r900 = Erebus.porosity_decrease_rate(0.38, 0.40, 1.0, 900.0)
            @test 10.0 < 0.01 / r900 / _HK_SECONDS_PER_YEAR < 30.0
            # Arrhenius: hotter is faster, by exp(E/R_gas (1/900 - 1/1000)) = 116
            r1000 = Erebus.porosity_decrease_rate(0.38, 0.40, 1.0, 1000.0)
            @test isapprox(
                r1000 / r900, exp(85000 / 1.9872 * (1 / 900 - 1 / 1000)); rtol=1e-10
            )
            # Eq. 43: below phi = 0.1 the rate carries 10 phi; continuous at 0.1
            R, a = Erebus.grain_geometry(0.08, 0.40, R0)
            s = Erebus.contact_stress(5.0, R, a)
            r_low = Erebus.sintering_rate(0.08, s, R, 1000.0)
            r_nofactor =
                4.0e-5 * (1 - 0.08) * s^1.5 / (100R)^3 * exp(-85000 / (1.9872 * 1000))
            @test isapprox(r_low, 0.8 * r_nofactor; rtol=1e-12)
            @test isapprox(
                Erebus.sintering_rate(0.1, s, R, 1000.0),
                Erebus.sintering_rate(nextfloat(0.1), s, R, 1000.0);
                rtol=1e-12,
            )
            # Edge cases: no stress, no porosity, no pressure give no sintering
            @test iszero(Erebus.sintering_rate(0.3, 0.0, R, 1000.0))
            @test iszero(Erebus.porosity_decrease_rate(0.0, 0.40, 5.0, 1000.0))
            @test iszero(Erebus.porosity_decrease_rate(0.3, 0.40, 0.0, 1000.0))
            # Onset: unsintered grains have a = 0 and an unbounded rate
            @test isinf(Erebus.porosity_decrease_rate(0.40, 0.40, 5.0, 1000.0))
            @test_throws DomainError Erebus.sintering_rate(0.3, s, R, 0.0)
            @test_throws DomainError Erebus.sintering_rate(0.3, -1.0, R, 900.0)
            @test_throws DomainError Erebus.sintering_rate(1.2, s, R, 900.0)
            @test_throws DomainError Erebus.sintering_rate(0.3, s, -R, 900.0)
        end

        @testset "Sub-stepped backward Euler converges with order 1 to the fine integration" begin
            phi0, phi_start, P, T = 0.40, 0.38, 1.0, 900.0
            t_end = 200.0 * _HK_SECONDS_PER_YEAR
            ref = _hk_rk4(phi_start, phi0, P, T, t_end; nsteps=800)
            ref2 = _hk_rk4(phi_start, phi0, P, T, t_end; nsteps=400)
            # The reference is converged (RK4 change 2e-12) far below the errors here
            @test abs(ref - ref2) < big(1.0e-11)
            errs = Float64[]
            for dphi_max in (0.004, 0.002, 0.001)
                phi, nsub, capped = Erebus.sinter_porosity_step(
                    phi_start, phi0, P, T, t_end; dphi_max=dphi_max
                )
                @test !capped
                @test 2 <= nsub <= 200
                push!(errs, abs(phi - Float64(ref)))
            end
            # Total decrease over 200 yr is ~0.05: several substeps at dphi_max
            @test 0.02 < phi_start - Float64(ref) < 0.1
            # Order 1: the error halves with dphi_max (ratio 2, not 4 or 1)
            @test 1.6 < errs[1] / errs[2] < 2.5
            @test 1.6 < errs[2] / errs[3] < 2.5
            # Default 0.01 substeps: error 2.9e-3, on the same order-1 line
            phi_def, _, _ = Erebus.sinter_porosity_step(phi_start, phi0, P, T, t_end)
            @test abs(phi_def - Float64(ref)) < 4.0e-3
            @test 1.8 < abs(phi_def - Float64(ref)) / errs[1] < 2.8
            # Porosity decreases, stays positive
            @test 0.0 < phi_def < phi_start
        end

        @testset "Onset asymptote phi_0 - phi = C t^(2/5)" begin
            phi0, P, T = 0.40, 5.0, 1000.0
            R0, A, E, n = 1.0e-6, 4.0e-5, 85000.0, 1.5
            # C = ((n+1) K)^(1/(n+1)), K = A (1-phi0) (3 sqrt(3) (1-phi0) P / pi)^n
            #     (100 R0)^-3 exp(-E / (R_gas T)); see the derivation on the page
            K =
                A * (1 - phi0) * (3 * sqrt(3.0) * (1 - phi0) * P / pi)^n / (100R0)^3 *
                exp(-E / (1.9872 * T))
            C = ((n + 1) * K)^(1 / (n + 1))
            devs = Float64[]
            for delta in (1.0e-5, 1.0e-4, 1.0e-3)
                t = _hk_onset_time(delta, phi0, P, T)
                push!(devs, abs(delta / (C * Float64(t)^(2 / 5)) - 1))
            end
            # Relative deviation of the asymptote is O(delta): x10 in delta gives ~x10
            @test devs[1] < 1.0e-4
            @test devs[3] < 1.0e-2
            @test 5.0 < devs[2] / devs[1] < 20.0
            # Exponent guard: t^(1/2) or t^(1/3) miss the quadrature by > 10 %
            t3 = Float64(_hk_onset_time(1.0e-3, phi0, P, T))
            t4 = Float64(_hk_onset_time(1.0e-4, phi0, P, T))
            @test isapprox(log(1.0e-3 / 1.0e-4) / log(t3 / t4), 0.4; rtol=1e-2)
            # The Erebus step from phi_0 with fine substeps follows the asymptote
            t2 = Float64(_hk_onset_time(1.0e-2, phi0, P, T))
            phi_a, _, _ = Erebus.sinter_porosity_step(phi0, phi0, P, T, t2; dphi_max=1.0e-4)
            phi_b, _, _ = Erebus.sinter_porosity_step(phi0, phi0, P, T, t2; dphi_max=1.0e-3)
            @test isapprox(phi0 - phi_a, 1.0e-2; rtol=2.0e-2)
            @test abs(phi0 - phi_a - 1.0e-2) < abs(phi0 - phi_b - 1.0e-2)
            # One backward Euler step from phi_0 has the (2/5)^(2/5) onset error
            phi_1, _, _ = Erebus.sinter_porosity_step(phi0, phi0, P, T, t2; dphi_max=1.0)
            @test isapprox((phi0 - phi_1) / 1.0e-2, (2 / 5)^(2 / 5); rtol=5.0e-2)
        end

        @testset "Step contracts: monotone, bounded, short cuts" begin
            # P = 0 at phi = phi_0: no contact stress, no change
            phi, nsub, capped = Erebus.sinter_porosity_step(0.40, 0.40, 0.0, 1000.0, 1.0e12)
            @test isapprox(phi, 0.40; rtol=0, atol=0)
            @test nsub == 0 && !capped
            # phi = 0 stays 0
            phi0_, _, _ = Erebus.sinter_porosity_step(0.0, 0.40, 5.0, 1000.0, 1.0e12)
            @test iszero(phi0_)
            # Long, hot step: phi stays in [0, phi^n] and the substep cap is reported
            phi_l, nsub_l, capped_l = Erebus.sinter_porosity_step(
                0.40, 0.40, 50.0, 1200.0, 1.0e15; max_substeps=5
            )
            @test 0.0 <= phi_l < 0.40
            @test nsub_l == 5 && capped_l
            # Low porosity: Eq. 43 slows the decrease, phi never reaches below 0
            phi_s, _, _ = Erebus.sinter_porosity_step(0.05, 0.40, 50.0, 1200.0, 1.0e14)
            @test 0.0 < phi_s < 0.05
            @test_throws DomainError Erebus.sinter_porosity_step(0.5, 0.4, 1.0, 900.0, 1.0)
            @test_throws DomainError Erebus.sinter_porosity_step(0.3, 0.4, 1.0, 900.0, -1.0)
        end

        @testset "Porous conductivity, L16 Eq. 12 and markmart.c (3 branches, floor)" begin
            k = 4.3
            # Branch 1 (phi < 0.2): k exp(-phi / 0.08)
            @test isapprox(Erebus.k_porous(k, 0.1), k * exp(-1.25); rtol=1e-14)
            # Branch 3 (phi > 0.4): k exp(-1.2 - phi / 0.167), sign of a as in L16
            @test isapprox(Erebus.k_porous(k, 0.6), k * exp(-1.2 - 0.6 / 0.167); rtol=1e-14)
            # Branch 2 (0.2 <= phi <= 0.4): (k1^4 + k2^4)^(1/4), one-sided values
            k1(p) = k * exp(-p / 0.08)
            k2(p) = k * exp(-1.2 - p / 0.167)
            @test isapprox(
                Erebus.k_porous(k, 0.2), (k1(0.2)^4 + k2(0.2)^4)^0.25; rtol=1e-14
            )
            @test isapprox(Erebus.k_porous(k, prevfloat(0.2)), k1(0.2); rtol=1e-12)
            @test isapprox(
                Erebus.k_porous(k, 0.4), (k1(0.4)^4 + k2(0.4)^4)^0.25; rtol=1e-14
            )
            @test isapprox(Erebus.k_porous(k, nextfloat(0.4)), k2(0.4); rtol=1e-12)
            # The reference is not continuous at 0.2: the jump is +26 %
            @test 1.2 < Erebus.k_porous(k, 0.2) / Erebus.k_porous(k, prevfloat(0.2)) < 1.3
            # Pinned by hand: 4.3 exp(-1.25) and 4.3 exp(-1.2 - 0.6 / 0.167)
            @test isapprox(Erebus.k_porous(k, 0.1), 1.231971; atol=1e-6)
            @test isapprox(Erebus.k_porous(k, 0.6), 0.0356431; atol=1e-7)
            # Sign guard: a = +1.2 would give k(0.6) = 0.091 k instead of 0.0083 k
            @test 0.005 * k < Erebus.k_porous(k, 0.6) < 0.01 * k
            # Dense limit and monotone decrease
            @test isapprox(Erebus.k_porous(k, 0.0), k; rtol=1e-14)
            ps = 0.0:0.05:0.85
            @test all(
                diff([Erebus.k_porous(k, p) for p in ps if !isapprox(p, 0.2; atol=1e-12)]) .< 0.0,
            )
            # Floor 1e-3 W/(m K) with a low solid conductivity
            @test isapprox(Erebus.k_porous(0.01, 0.85), 1.0e-3; rtol=1e-14)
            @test_throws DomainError Erebus.k_porous(-1.0, 0.3)
            @test_throws DomainError Erebus.k_porous(k, 1.1)
        end

        @testset "Replay of the i2elvis porosity block (1 explicit step)" begin
            nsinter = 0
            for (phi, T, P, phi_init, dt, phi_c, change_c) in _HK_I2ELVIS_REPLAY
                # i2elvis stores phi and T as float; P [Pa] -> bar as markmart.c
                phi32 = Float64(Float32(phi))
                T32 = Float64(Float32(T))
                P_bar = P * 1.0e-5
                if T32 < 700.0
                    phi_e = min(phi32, Erebus.phi_cold(P_bar), 0.85)
                    @test isapprox(phi_e, phi_c; rtol=1e-6)
                    @test iszero(change_c)
                else
                    nsinter += 1
                    phi_0 = Erebus.porosity_phi0(P_bar, phi_init)
                    change_e = Erebus.porosity_decrease_rate(phi32, phi_0, P_bar, T32) * dt
                    @test isapprox(change_e, change_c; rtol=1e-6)
                    @test isapprox(phi32 - change_e, phi_c; rtol=1e-6)
                    # 1 step changes phi by 1e-6 to 3 % of phi: the comparison is not trivial
                    @test 1.0e-6 < change_c / phi32 < 0.03
                end
            end
            @test nsinter == 5
        end
    end
end

"""
Return the error thrown by `f()`, or `nothing`.
"""
function _hk_error(f)
    try
        f()
    catch e
        return e
    end
    return nothing
end

@testset "Configuration section [porosity]" begin
    @testset "Defaults: model none and the reference constants" begin
        cfg = Erebus.default_config()
        @test cfg.porosity.model === :none
        @test isapprox(cfg.porosity.phi_cap, 0.85; rtol=1e-14)
        @test isapprox(cfg.porosity.grain_radius, 1.0e-6; rtol=1e-14)
        @test isapprox(cfg.porosity.T_sinter, 700.0; rtol=1e-14)
        @test isapprox(cfg.porosity.A, 4.0e-5; rtol=1e-14)
        @test isapprox(cfg.porosity.E_act, 85000.0; rtol=1e-14)
        @test isapprox(cfg.porosity.n, 1.5; rtol=1e-14)
        @test haskey(Erebus.config_to_dict(cfg), "porosity")
    end

    @testset "henke2012 parses, round-trips, and needs the coupled processes off" begin
        toml = """
        [porosity]
        model = "henke2012"
        phi_init = 0.5
        [reaction]
        active = false
        """
        cfg = Erebus.load_config(toml)
        @test cfg.porosity.model === :henke2012
        @test isapprox(cfg.porosity.phi_init, 0.5; rtol=1e-14)
        cfg_rt = Erebus.parse_config_string(Erebus.serialize_config(cfg))
        @test cfg_rt.porosity.model === :henke2012
        @test isapprox(cfg_rt.porosity.phi_init, 0.5; rtol=1e-14)
        # The reaction switch is on by default, so henke2012 alone is rejected
        err = _hk_error(() -> Erebus.load_config("[porosity]\nmodel = \"henke2012\"\n"))
        @test err isa ArgumentError
        @test occursin("henke2012", err.msg)
        base = Erebus.default_config()
        por = Erebus.PorosityConfig(; model=:henke2012)
        react_off = Erebus.ReactionConfig(; active=false)
        ok = Erebus.SimulationConfig(; porosity=por, reaction=react_off)
        @test Erebus.validate_config(ok) === nothing
        for (field, value) in (
            (:venting, Erebus.VentingConfig(; active=true)),
            (:hydrothermal, Erebus.HydrothermalConfig(; active=true)),
            (:accretion, Erebus.AccretionConfig(; active=true)),
            (:telescoping, Erebus.TelescopingConfig(; active=true)),
        )
            kw = Dict{Symbol,Any}(:porosity => por, :reaction => react_off, field => value)
            e = _hk_error(() -> Erebus.validate_config(Erebus.SimulationConfig(; kw...)))
            @test e isa ArgumentError
            @test occursin("henke2012", e.msg)
        end
        # The same switches are allowed with model none
        @test Erebus.validate_config(base) === nothing
    end

    @testset "Field contracts of [porosity]" begin
        react_off = Erebus.ReactionConfig(; active=false)
        for kw in (
            (; model=:henke2013),
            (; phi_init=1.2),
            (; phi_init=-0.1),
            (; phi_cap=1.0),
            (; phi_cap=0.0),
            (; grain_radius=0.0),
            (; T_sinter=0.0),
            (; A=-1.0),
            (; E_act=-85000.0),
            (; n=0.0),
        )
            por = Erebus.PorosityConfig(; kw...)
            e = _hk_error(
                () -> Erebus.validate_config(
                    Erebus.SimulationConfig(; porosity=por, reaction=react_off)
                ),
            )
            @test e isa ArgumentError
            @test occursin("porosity", e.msg)
        end
    end
end

"""
Small henke2012 configuration (17 x 17 grid, dry rock, coupled processes off).
"""
function _hk_cfg(overrides::Pair...)
    d = Dict{String,Any}(
        "grid.Nx" => 17,
        "grid.Ny" => 17,
        "porosity.model" => :henke2012,
        "porosity.phi_init" => 0.75,
        "reaction.active" => false,
        "venting.active" => false,
        "materials.XWsolidm_init" => [0.0, 0.0, NaN],
        "output.mode" => :telemetry,
        "output.save_final" => false,
        "output.output_dir" => tempname(),
        "solver.seed" => 7,
    )
    for (k, v) in overrides
        d[k] = v
    end
    return Erebus.override_config(Erebus.default_config(), d)
end

"""
Initialize a simulation state for `cfg` and close the telemetry stream.
"""
function _hk_state(cfg)
    st, coords, ws, tel, _ = Erebus.init_simulation(cfg)
    tel !== nothing && close(tel)
    return st, coords, ws
end

@testset "Dry porosity marker group and step (model henke2012)" begin
    @testset "Group exists only with henke2012 and draws no random numbers" begin
        coords = Erebus.GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        rng_a = Erebus.Random.MersenneTwister(3)
        rng_b = Erebus.Random.MersenneTwister(3)
        m_none = Erebus.init_marker_arrays(
            64, _hk_cfg("porosity.model" => :none), coords; rng=rng_a
        )
        m_h = Erebus.init_marker_arrays(64, _hk_cfg(), coords; rng=rng_b)
        @test !haskey(m_none.groups, :porosity)
        @test !haskey(m_none, :phi_d)
        @test haskey(m_h.groups, :porosity)
        @test length(m_h.phi_d) == 64
        @test length(m_h.tkm_max) == 64 && length(m_h.dphi_sinter) == 64
        @test isapprox(rand(rng_a), rand(rng_b); rtol=0, atol=0)
        # The full setup draws the same marker positions and temperatures
        st_n, _, _ = _hk_state(_hk_cfg("porosity.model" => :none))
        st_h, _, _ = _hk_state(_hk_cfg())
        @test isapprox(st_n.markers.xm, st_h.markers.xm; rtol=0, atol=0)
        @test isapprox(st_n.markers.tkm, st_h.markers.tkm; rtol=0, atol=0)
        @test isapprox(rand(st_n.rng), rand(st_h.rng); rtol=0, atol=0)
    end

    @testset "Setup pressure: uniform-density limit of dP/dr = -rho g" begin
        G = Erebus.G_GRAV
        for (mode, c_grav) in ((:enclosed_mass, 2 / 3), (:poisson2d, 1.0))
            # phi_init = 0.3 < 0.42 <= phi_H(P): phi_0 = 0.3 at every depth
            cfg = _hk_cfg("porosity.phi_init" => 0.3, "geometry.gravity_mode" => mode)
            r, P = Erebus.dry_hydrostatic_profile(cfg)
            R = cfg.geometry.rplanet
            rho = 3300.0 * (1 - 0.3)
            P_an = c_grav * pi * G * rho^2 .* (R^2 .- r .^ 2)
            @test isapprox(P[1], P_an[1]; rtol=1e-6)
            @test maximum(abs.(P .- P_an)) < 1.0e-6 * P_an[1]
            # P(R_p) = 0 and dP/dr <= 0
            @test isapprox(P[end], 0.0; atol=1e-9 * P[1])
            @test all(diff(P) .<= 0.0)
            # Scale guard: ~18.6 bar (sphere) or ~28 bar (cylinder) at the centre of 50 km
            @test 1.0e6 < P[1] < 4.0e6
        end
        # With phi_init = 0.75 the deep rock is cold pressed: denser, higher pressure
        r, P = Erebus.dry_hydrostatic_profile(_hk_cfg())
        rho_top = 3300.0 * 0.25
        P_top_only = (2 / 3) * pi * G * rho_top^2 * r[end]^2
        @test P[1] > 2.0 * P_top_only
    end

    @testset "Setup: water porosity at phimin, dry porosity cold pressed" begin
        cfg = _hk_cfg()
        st, coords, _ = _hk_state(cfg)
        mk = st.markers
        rock = findall(<(3), mk.tm)
        air = findall(==(3), mk.tm)
        @test all(isapprox.(mk.phim[rock], cfg.poroelasticity.phimin; rtol=1e-14))
        @test all(iszero, mk.phi_d[air])
        @test all(iszero, mk.dphi_sinter)
        @test isapprox(mk.tkm_max[rock], mk.tkm[rock]; rtol=1e-14)
        r_prof, P_prof = Erebus.dry_hydrostatic_profile(cfg)
        xc, yc = cfg.geometry.xcenter, cfg.geometry.ycenter
        for m in rock[1:7:end]
            rm = hypot(mk.xm[m] - xc, mk.ym[m] - yc)
            Pm = Erebus._interp_profile(r_prof, P_prof, rm)
            @test isapprox(mk.phi_d[m], Erebus.porosity_phi0(Pm * 1.0e-5, 0.75); rtol=1e-12)
        end
        # Deep rock is pressed to near 0.42, the outer rock keeps 0.75
        @test 0.42 < minimum(mk.phi_d[rock]) < 0.45
        @test isapprox(maximum(mk.phi_d[rock]), 0.75; rtol=1e-12)
    end

    @testset "Marker update: cold pressing, threshold, hysteresis, P = 0" begin
        pc = _hk_cfg().porosity
        dt = 1.0e3 * 3.15576e7
        # Cold: phi_H(0.13 bar) = 0.65, no sintering below 700 K
        u = Erebus.dry_porosity_update(0.75, 0.13, 650.0, dt, pc)
        @test isapprox(u.phi, 0.65; rtol=1e-14)
        @test iszero(u.dphi_sinter)
        # Threshold: 699.9 K does not sinter, 700 K does (phi_0 = phi_H(1 bar) = 0.433)
        below = Erebus.dry_porosity_update(0.4, 1.0, 699.9, dt, pc)
        at = Erebus.dry_porosity_update(0.4, 1.0, 700.0, dt, pc)
        @test isapprox(below.phi, 0.4; rtol=1e-14)
        @test at.phi < 0.4
        @test isapprox(at.dphi_sinter, 0.4 - at.phi; rtol=1e-12)
        # Hysteresis: 650 K -> 900 K -> 600 K, porosity does not recover
        s1 = Erebus.dry_porosity_update(0.75, 2.0, 650.0, dt, pc)
        s2 = Erebus.dry_porosity_update(s1.phi, 2.0, 900.0, dt, pc)
        s3 = Erebus.dry_porosity_update(s2.phi, 2.0, 600.0, dt, pc)
        @test s2.phi < s1.phi
        @test isapprox(s3.phi, s2.phi; rtol=1e-14)
        @test iszero(s3.dphi_sinter)
        @test isapprox(s2.dphi_sinter, s1.phi - s2.phi; rtol=1e-12)
        # P = 0 at phi = phi_0: no contact stress, no sintering even when hot
        z = Erebus.dry_porosity_update(0.75, 0.0, 1200.0, dt, pc)
        @test isapprox(z.phi, 0.75; rtol=1e-14)
        @test iszero(z.dphi_sinter)
        @test !z.capped
        @test_throws DomainError Erebus.dry_porosity_update(0.5, -1.0, 900.0, dt, pc)
    end

    @testset "Step in vent_and_degas!: total pressure grids.pr, APHI skipped" begin
        cfg = _hk_cfg()
        st, coords, _ = _hk_state(cfg)
        mk = st.markers
        rock = findall(<(3), mk.tm)
        air = findall(==(3), mk.tm)
        st.timestep = 2
        st.dt = 1.0e6
        fill!(st.grids.DT, 0.0)
        fill!(st.grids.APHI, 1.0e-9)
        P_pa = 2.0e5
        fill!(st.grids.pr, cfg.geometry.psurface + P_pa)
        fill!(st.grids.pf, 0.0)
        mk.tkm[rock] .= 1000.0
        phi_before = copy(mk.phi_d)
        Erebus.vent_and_degas!(st, coords, cfg)
        for m in rock
            exp_u = Erebus.dry_porosity_update(
                phi_before[m], 2.0, 1000.0, 1.0e6, cfg.porosity
            )
            @test isapprox(mk.phi_d[m], exp_u.phi; rtol=1e-12)
        end
        @test all(mk.phi_d[rock] .< phi_before[rock])
        @test all(isapprox.(mk.tkm_max[rock], 1000.0; rtol=1e-14))
        @test all(mk.dphi_sinter[rock] .> 0.0)
        # Water porosity stays at phimin: the compaction rate APHI is not applied
        @test all(isapprox.(mk.phim[rock], cfg.poroelasticity.phimin; rtol=1e-14))
        @test all(iszero, mk.phi_d[air])
    end

    @testset "Rock properties use (1 - phi_d) and k_porous at every call site" begin
        # compute_marker_properties!: density, rho cp, conductivity, both radiogenic sites
        phi_d = [0.6]
        arrs() = (zeros(1), zeros(1), zeros(1), zeros(1), zeros(1), zeros(1), zeros(1))
        hrs = [1.0e-6, 2.0e-6, 0.0]
        hrf = [5.0e-7, 5.0e-7, 0.0]
        for hrmetal in (nothing, [4.0e-6, 4.0e-6, 0.0])
            rho, rhocp, eta, hr, k, tkrc, eik = arrs()
            Erebus.compute_marker_properties!(
                1,
                [2],
                [300.0],
                rho,
                rhocp,
                eta,
                hr,
                k,
                tkrc,
                eik,
                hrs,
                hrf,
                [1.0e-4],
                [0.0],
                1;
                compute_hr=true,
                hrmetalm=hrmetal,
                phi_dry=phi_d,
            )
            @test isapprox(rho[1], 3300.0 * 0.4; rtol=1e-12)
            @test isapprox(rhocp[1], 3.3e6 * 0.4; rtol=1e-12)
            @test isapprox(
                k[1], Erebus.k_porous(Erebus.compute_ksolidm(300.0, 1), 0.6); rtol=1e-12
            )
            # No pore-fluid term: the dense dry value times (1 - phi_d)
            rho_w, rhocp_w, eta_w, hr_w, k_w, tkrc_w, eik_w = arrs()
            Erebus.compute_marker_properties!(
                1,
                [2],
                [300.0],
                rho_w,
                rhocp_w,
                eta_w,
                hr_w,
                k_w,
                tkrc_w,
                eik_w,
                hrs,
                hrf,
                [1.0e-4],
                [0.0],
                1;
                compute_hr=true,
                hrmetalm=hrmetal,
                phi_dry=[0.0],
            )
            @test isapprox(hr[1], 0.4 * hr_w[1]; rtol=1e-12)
            @test isapprox(rho[1], 0.4 * rho_w[1]; rtol=1e-12)
        end

        # The 3 particle-to-mesh paths give the same marker properties
        for mode in (:tiled, :buffered)
            cfg = _hk_cfg("solver.p2m_mode" => mode)
            st, coords, ws = _hk_state(cfg)
            Erebus.interpolate_markers_to_grid!(
                st,
                coords,
                cfg;
                thread_buffers=ws.thread_buffers,
                interp_arrays=ws.interp_arrays,
            )
            mk = st.markers
            rock = findall(<(3), mk.tm)
            @test isapprox(mk.rhototalm[rock], 3300.0 .* (1 .- mk.phi_d[rock]); rtol=1e-12)
            k_exp = [
                Erebus.k_porous(
                    Erebus.compute_ksolidm(mk.tkm[m], Erebus.marker_property_mode),
                    mk.phi_d[m],
                ) for m in rock
            ]
            @test isapprox(mk.ktotalm[rock], k_exp; rtol=1e-12)
        end

        # Sub-grid diffusion relaxes with the dry rho cp and k_porous
        coords = Erebus.GridCoordinates(9, 9; xsize=10_000.0, ysize=10_000.0)
        Ny1, Nx1 = coords.Ny1, coords.Nx1
        for phid in (0.0, 0.6)
            tkm = [400.0]
            tk1 = fill(300.0, Ny1, Nx1)
            Erebus.apply_subgrid_temperature_diffusion!(
                [5000.0],
                [5000.0],
                [2],
                tkm,
                [1.0e-4],
                tk1,
                zeros(Ny1, Nx1),
                zeros(Ny1, Nx1),
                zeros(Ny1, Nx1),
                1.0e9,
                1,
                1;
                coords=coords,
                dsubgridt=1.0,
                phi_dry=[phid],
            )
            k_d = Erebus.k_porous(Erebus.compute_ksolidm(400.0, 1), phid)
            rc_d = 3.3e6 * (1 - phid)
            fac = exp(-k_d * 1.0e9 / rc_d * 2 * (1 / coords.dx^2 + 1 / coords.dy^2))
            @test isapprox(tkm[1], 300.0 + 100.0 * fac; rtol=1e-12)
        end

        # Radiogenic heating scales with (1 - phi_d) and has no pore-water term
        cfg = _hk_cfg()
        st, coords, _ = _hk_state(cfg)
        rock = findall(<(3), st.markers.tm)
        st.markers.phi_d[rock] .= 0.0
        Erebus.radiogenic_heating!(st, coords, cfg)
        hr_dense = copy(st.markers.hrtotalm[rock])
        st.markers.phi_d[rock] .= 0.6
        Erebus.radiogenic_heating!(st, coords, cfg)
        @test isapprox(st.markers.hrtotalm[rock], 0.4 .* hr_dense; rtol=1e-12)
        @test all(hr_dense .> 0.0)

        # Surface radiation face conductivity uses k_porous at the surface porosity
        coords = Erebus.GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        KX0 = ones(coords.Ny1, coords.Nx1)
        KY0 = ones(coords.Ny1, coords.Nx1)
        tk = fill(300.0, coords.Ny1, coords.Nx1)
        KXd, KYd = copy(KX0), copy(KY0)
        Erebus.apply_radiative_surface_boundary!(
            KXd, KYd, tk, coords, 50_000.0, 70_000.0, 70_000.0, 170.0; phi_dry_surface=0.75
        )
        KXs, KYs = copy(KX0), copy(KY0)
        Erebus.apply_radiative_surface_boundary!(
            KXs, KYs, tk, coords, 50_000.0, 70_000.0, 70_000.0, 170.0
        )
        face = findfirst(!isapprox(1.0), KXd)
        h = Erebus.compute_effective_radiation_htc(300.0, 170.0, 0.0)
        kb = Erebus.k_porous(Erebus.compute_ksolidm(300.0, 1), 0.75)
        @test isapprox(KXd[face], 2kb * h * coords.dx / (2kb + h * coords.dx); rtol=1e-12)
        @test KXd[face] < KXs[face]
    end

    @testset "Marker lifecycle: replenish and checkpoint keep the group" begin
        cfg = _hk_cfg()
        coords = Erebus.GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        markers = Erebus.init_marker_arrays(10, cfg, coords)
        markers.phi_d .= range(0.1, 0.55; length=10)
        markers.tkm_max .= 900.0
        markers.dphi_sinter .= 0.02
        mdis = fill(1.0, coords.Ny, coords.Nx)
        mnum = zeros(Int, coords.Ny, coords.Nx)
        n_new = Erebus.replenish_markers!(markers, mdis, mnum; coords=coords, cfg=cfg)
        @test 10 < n_new < 100_000
        @test length(markers.phi_d) == n_new
        @test length(markers.tkm_max) == n_new && length(markers.dphi_sinter) == n_new
        # New markers copy a source marker: values stay in the source set
        @test all(0.1 - 1e-12 .<= markers.phi_d .<= 0.55 + 1e-12)
        @test all(isapprox.(markers.tkm_max, 900.0; rtol=1e-14))
        @test all(isapprox.(markers.dphi_sinter, 0.02; rtol=1e-14))
        # Checkpoint restore rebuilds the group from the saved arrays
        st, coords2, _ = _hk_state(cfg)
        st.markers.phi_d .= 0.33
        path = joinpath(mktempdir(), "ckpt.jld2")
        Erebus.save_state(path, st, coords2, cfg)
        data = Erebus.JLD2.load(path)
        restored = Erebus.restore_checkpoint_markers(data, st.markers.core, cfg)
        @test haskey(restored.groups, :porosity)
        @test all(isapprox.(restored.phi_d, 0.33; rtol=1e-14))
        # A checkpoint of a model-none run cannot restart a henke2012 run silently
        st_n, coords_n, _ = _hk_state(_hk_cfg("porosity.model" => :none))
        path_n = joinpath(mktempdir(), "ckpt_none.jld2")
        Erebus.save_state(path_n, st_n, coords_n, _hk_cfg("porosity.model" => :none))
        data_n = Erebus.JLD2.load(path_n)
        @test_throws Erebus.CheckpointError Erebus.restore_checkpoint_markers(
            data_n, st_n.markers.core, cfg
        )
        restored_n = Erebus.restore_checkpoint_markers(
            data_n, st_n.markers.core, _hk_cfg("porosity.model" => :none)
        )
        @test !haskey(restored_n.groups, :porosity)
    end

    @testset "Two-step run: Stokes pressure drives sintering of hot rock" begin
        cfg = _hk_cfg(
            "time.n_steps" => 2,
            "time.dt_initial" => 1.0,
            "time.dt_longest" => 1.0,
            "materials.tkm0" => [1000.0, 1000.0, 170.0],
            "thermodynamics.surface_radiation" => true,
        )
        st0, _, _ = _hk_state(cfg)
        phi_start = copy(st0.markers.phi_d)
        st = Erebus.with_logger(Erebus.NullLogger()) do
            Erebus.simulation_loop(cfg; output_path=mktempdir())
        end
        mk = st.markers
        rock = findall(<(3), mk.tm)
        @test st.timestep == 2
        @test all(isfinite, mk.phi_d) && all(isfinite, mk.tkm)
        # Porosity stays in [0, phi_init] and never rises; water porosity stays at phimin
        @test all(0.0 .<= mk.phi_d[rock] .<= 0.75)
        @test all(isapprox.(mk.phim[rock], cfg.poroelasticity.phimin; rtol=1e-14))
        # The interior sinters under the Stokes pressure: phi_d falls, dphi_sinter records it
        @test 0.0 < maximum(mk.dphi_sinter[rock]) < 0.5
        @test minimum(mk.phi_d[rock]) < minimum(phi_start[rock])
        @test all(mk.tkm_max[rock] .>= 700.0)
    end
end

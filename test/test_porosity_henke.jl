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
# ("Effect of cold pressing" to the [0, 1] clamp, the anchors the script checks), compiled
# with cc. Columns: phi, T [K], P [Pa], phi_init, dt [s], phi after 1 step, change.
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

@testset "Porosity laws of Henke et al. (2012)" begin
    setprecision(BigFloat, 128) do
        @testset "Cold pressing law phi_H(P), Henke Eq. 15 and L16 Eq. 10" begin
            # Pinned by hand: 0.42 + 0.46 = 0.88 at P = 0 and 0.42 + 0.46 / 2 at P = P0
            @test isapprox(Erebus.phi_cold(0.0), 0.88; rtol=1e-14)
            @test isapprox(Erebus.phi_cold(0.13), 0.65; rtol=1e-14)
            # 1.3 bar = 10 P0: 0.42 + 0.46 / (10^1.72 + 1) = 0.42860123
            phi_13 = Erebus.phi_cold(1.3)
            @test isapprox(phi_13, 0.4286012; atol=5e-8)
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
            # The cap has the domain (0, 1) of the configuration check
            @test_throws DomainError Erebus.porosity_phi0(1.0, 0.5; phi_cap=1.0)
            @test_throws DomainError Erebus.porosity_phi0(1.0, 0.5; phi_cap=0.0)
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
            # By hand: the BigFloat root gives R / R0 = 1.000128687, (a / R0)^2 = 0.0226049357;
            # sigma / P = 2 sqrt(3) (1.000257391 - 0.0226049357) / (pi 0.0226049357) = 47.6894122
            @test isapprox(sigma, 47.6894122; atol=1e-6)
            # Prefactor guards: sqrt(3) in place of 2 sqrt(3), or no pi, miss by > 50 %
            @test !isapprox(sigma, 47.6894122 / 2; rtol=0.2)
            @test !isapprox(sigma, 47.6894122 * pi; rtol=0.2)
            # Sign guard: C_av with R^2 + a^2 gives sigma / P = 49.894, 4.6 % higher
            @test !isapprox(sigma, 49.894; rtol=0.02)
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
            # Stress exponent guard: n = 1 in place of 1.5 changes the rate by s^0.5 = 3.75
            r_n1 = Erebus.sintering_rate(0.08, s, R, 1000.0; n=1.0)
            @test isapprox(r_low / r_n1, sqrt(s); rtol=1e-12)
            @test 3.0 < sqrt(s) < 5.0
            @test_throws DomainError Erebus.porosity_decrease_rate(0.3, 0.40, -1.0, 1000.0)
            @test_throws DomainError Erebus.porosity_decrease_rate(0.3, 0.40, 1.0, 0.0)
            @test_throws DomainError Erebus.porosity_decrease_rate(0.5, 0.40, 1.0, 900.0)
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
            # dt = 0 short cut: no substep, phi unchanged even at the onset (rate Inf)
            phi_z, nsub_z, capped_z = Erebus.sinter_porosity_step(
                0.40, 0.40, 5.0, 1000.0, 0.0
            )
            @test isapprox(phi_z, 0.40; rtol=0, atol=0)
            @test nsub_z == 0 && !capped_z
            @test_throws DomainError Erebus.sinter_porosity_step(0.5, 0.4, 1.0, 900.0, 1.0)
            @test_throws DomainError Erebus.sinter_porosity_step(0.3, 0.4, 1.0, 900.0, -1.0)
            @test_throws ArgumentError Erebus.sinter_porosity_step(
                0.3, 0.4, 1.0, 900.0, 1.0; max_substeps=0
            )
            @test_throws DomainError Erebus.sinter_porosity_step(
                0.3, 0.4, 1.0, 900.0, 1.0; dphi_max=0.0
            )
            @test_throws DomainError Erebus.sinter_porosity_step(
                0.3, 0.4, 1.0, 900.0, 1.0; dphi_max=-0.01
            )
            # A rate that overflows to Inf below the onset shrinks the substep to 0
            @test_throws DomainError Erebus.sinter_porosity_step(
                0.5, 0.6, 1.0, 950.0, 3.15576e10; A=1.0e300
            )
            # dphi_max below DPHI_SINTERED cannot step over the rounding of 1 - phi near
            # X = 1 (rate Inf on a few ULP below phi_0): the contract rejects it up front
            for (p, dm) in ((0.1, 5.0e-17), (0.42, 1.0e-16), (0.3, prevfloat(1.0e-12)))
                e = _hk_error(
                    () -> Erebus.sinter_porosity_step(p, p, 1.0, 900.0, 1.0e10; dphi_max=dm)
                )
                @test e isa DomainError
                @test occursin("dphi_max", e.msg)
                @test !occursin("overflow", e.msg)
            end
            # dphi_max must be finite; a BigFloat converts to Float64 first
            for dm in (Inf, big"1e400")
                @test_throws DomainError Erebus.sinter_porosity_step(
                    0.3, 0.4, 1.0, 900.0, 1.0e10; dphi_max=dm
                )
            end
            phi_b, nsub_b, cap_b = Erebus.sinter_porosity_step(
                0.3, 0.4, 1.0, 900.0, 1.0e10; dphi_max=big(0.01)
            )
            phi_f, nsub_f, cap_f = Erebus.sinter_porosity_step(
                0.3, 0.4, 1.0, 900.0, 1.0e10; dphi_max=0.01
            )
            @test isapprox(phi_b, phi_f; rtol=0, atol=0)
            @test nsub_b === nsub_f && cap_b === cap_f
            # Edge case: dphi_max = DPHI_SINTERED steps over the rounding and sinters
            phi_m, nsub_m, _ = Erebus.sinter_porosity_step(
                0.42, 0.42, 1.0, 900.0, 1.0e10; dphi_max=Erebus.DPHI_SINTERED
            )
            @test 0.0 < phi_m < 0.42
            @test 1 <= nsub_m <= 1000
            # A rate far above 1/h: 1 substep closes the pores to the smallest positive
            # float, not to phi_old 2^-200 (the bisection runs to adjacent floats)
            phi_c, nsub_c, capped_c = Erebus.sinter_porosity_step(
                0.005, 0.40, 2.0, 1000.0, 1.0e100; A=1.0e250
            )
            @test 0.0 <= phi_c <= nextfloat(0.0)
            @test nsub_c == 1 && !capped_c
            # h = 0 returns phi_old exactly, also at phi_0 where the rate is Inf
            r01(p) = Erebus.porosity_decrease_rate(p, 0.1, 1.0, 900.0)
            @test isapprox(
                Erebus._backward_euler_porosity(0.1, 0.0, r01), 0.1; rtol=0, atol=0
            )
            @test isapprox(
                Erebus._backward_euler_porosity(0.05, 0.0, r01), 0.05; rtol=0, atol=0
            )
        end

        @testset "An overflowing rate throws within 5 s, not a subnormal stall" begin
            # With dphi_max < dphi < 1.8 dphi_max the shrink factor is above 0.5, so the
            # smallest subnormal substep times it rounds back to itself
            overflow_cases = (
                (0.015, 0.6, 1.0, 950.0, 3.15576e10, (; A=1.0e300)),
                (0.015, 0.015, 1.0, 950.0, 3.15576e10, (; A=1.0e305)),
                (0.75, 0.75, 2.0, 1000.0, 3.0e15, (; A=1.0e300, dphi_max=0.5)),
                (0.015, 0.6, 1.0, 950.0, 3.15576e10, (; R0=1.0e-110)),
            )
            for (phi_n, phi_0, P_bar, T, dt, kw) in overflow_cases
                t_wall = @elapsed e = _hk_error(
                    () -> Erebus.sinter_porosity_step(phi_n, phi_0, P_bar, T, dt; kw...)
                )
                @test e isa DomainError
                @test occursin("overflow", e.msg)
                @test isapprox(e.val.phi_n, phi_n; rtol=0, atol=0)
                @test t_wall < 5.0
            end
            # The default constants never reach the guard: 1 kyr at 1 bar, 950 K sinters
            phi, nsub, capped = Erebus.sinter_porosity_step(
                0.015, 0.6, 1.0, 950.0, 3.15576e10
            )
            @test 0.0 < phi < 0.015
            @test nsub >= 1 && !capped
        end

        @testset "Porous conductivity, L16 Eq. 12 and markmart.c (3 branches, floor)" begin
            k = 4.3
            # Branch 1 (phi < 0.2): k exp(-phi / 0.08)
            @test isapprox(Erebus.k_porous(k, 0.1), k * exp(-1.25); rtol=1e-14)
            # Branch 3 (phi > 0.4): k exp(-1.2 - phi / 0.167), sign of a as in L16
            @test isapprox(Erebus.k_porous(k, 0.6), k * exp(-1.2 - 0.6 / 0.167); rtol=1e-14)
            # Branch 2 (0.2 <= phi <= 0.4), by hand at phi = 0.3: k1 = 4.3 exp(-3.75) =
            # 0.1011263, k2 = 4.3 exp(-1.2 - 1.796407) = 0.2148549, (k1^4 + k2^4)^(1/4)
            k_03 = Erebus.k_porous(k, 0.3)
            @test isapprox(k_03, 0.2174439; atol=1e-7)
            # Exponent guards: the power-2 mean gives 0.2374641 (+9 %), max(k1, k2) -1.2 %
            @test !isapprox(k_03, 0.2374641; rtol=0.05)
            @test !isapprox(k_03, 0.2148549; rtol=5.0e-3)
            # Branch ends by hand: 0.4441049 at 0.2 and 0.1181631 at 0.4; one-sided
            # values k1(0.2) = 0.3529655 and k2(0.4) = 0.1180562
            @test isapprox(Erebus.k_porous(k, 0.2), 0.4441049; atol=1e-7)
            @test isapprox(Erebus.k_porous(k, prevfloat(0.2)), 0.3529655; atol=1e-7)
            @test isapprox(Erebus.k_porous(k, 0.4), 0.1181631; atol=1e-7)
            @test isapprox(Erebus.k_porous(k, nextfloat(0.4)), 0.1180562; atol=1e-7)
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
Copy of the keyword-built configuration struct `x` with the fields in `kw` replaced.
"""
function _hk_with(x; kw...)
    names = fieldnames(typeof(x))
    vals = NamedTuple{names}(Tuple(getfield(x, f) for f in names))
    return typeof(x)(; merge(vals, (; kw...))...)
end

"""
Largest residual of dP/dr + rho(P) g over max(rho g) for a setup profile `(r, P)`.
`w = 2` takes the gravity G M(r) / r^2 of a sphere with the trapezoid mass integral,
`w = 1` the 2-D Poisson gravity (8/3) pi G m(r) / r with m(r) the integral of rho r'.
"""
function _hk_hydro_residual(r, P, phi_init, w)
    G = Erebus.G_GRAV
    rho = [3300.0 * (1 - Erebus.porosity_phi0(p * 1.0e-5, phi_init)) for p in P]
    m = zeros(length(r))
    for i in 2:length(r)
        m[i] =
            m[i - 1] + 0.5 * (rho[i - 1] * r[i - 1]^w + rho[i] * r[i]^w) * (r[i] - r[i - 1])
    end
    c = w == 2 ? 4pi * G : (8 / 3) * pi * G
    g = [i == 1 ? 0.0 : c * m[i] / r[i]^w for i in eachindex(r)]
    rhog = rho .* g
    res = [
        (P[i + 1] - P[i]) / (r[i + 1] - r[i]) + 0.5 * (rhog[i] + rhog[i + 1]) for
        i in 1:(length(r) - 1)
    ]
    return maximum(abs.(res)) / maximum(rhog)
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
        [materials]
        XWsolidm_init = [0.0, 0.0, nan]
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
        dry = Erebus.MaterialConfig(; XWsolidm_init=[0.0, 0.0, NaN])
        ok = Erebus.SimulationConfig(; porosity=por, reaction=react_off, materials=dry)
        @test Erebus.validate_config(ok) === nothing
        for (field, value) in (
            (:venting, Erebus.VentingConfig(; active=true)),
            (:hydrothermal, Erebus.HydrothermalConfig(; active=true)),
            (:accretion, Erebus.AccretionConfig(; active=true)),
            (:telescoping, Erebus.TelescopingConfig(; active=true)),
        )
            kw = Dict{Symbol,Any}(
                :porosity => por, :reaction => react_off, :materials => dry, field => value
            )
            e = _hk_error(() -> Erebus.validate_config(Erebus.SimulationConfig(; kw...)))
            @test e isa ArgumentError
            @test occursin("henke2012", e.msg)
        end
        # The same switches are allowed with model none
        @test Erebus.validate_config(base) === nothing
    end

    @testset "henke2012 needs dry rock and no volatile, magma or iron process" begin
        por = Erebus.PorosityConfig(; model=:henke2012)
        react_off = Erebus.ReactionConfig(; active=false)
        dry = Erebus.MaterialConfig(; XWsolidm_init=[0.0, 0.0, NaN])
        check(kw) = _hk_error(
            () -> Erebus._validate_porosity_config(
                Erebus.SimulationConfig(;
                    porosity=por, reaction=react_off, materials=dry, kw...
                ),
            ),
        )
        @test check((;)) === nothing
        # Edge case: refractory budgets without pyrolysis kinetics leave phim alone
        @test check((; refractory=Erebus.RefractoryConfig(; active=true))) === nothing
        # Kinetics switch of an inactive refractory section: pyrolysis does not run
        @test check((;
            refractory=Erebus.RefractoryConfig(; active=false, kinetics_active=true)
        )) === nothing
        # Each excluded process, with the name the message must give
        for (name, kw) in (
            ("XWsolidm_init", (; materials=Erebus.MaterialConfig())),
            (
                "XWsolidm_init",
                (; materials=Erebus.MaterialConfig(; XWsolidm_init=[0.1, 0.0, NaN])),
            ),
            (
                "XWsolidm_init",
                (; materials=Erebus.MaterialConfig(; XWsolidm_init=[0.0, 0.1, NaN])),
            ),
            ("volatiles", (; volatiles=Erebus.VolatilesConfig(; active=true))),
            (
                "refractory",
                (; refractory=Erebus.RefractoryConfig(; active=true, kinetics_active=true)),
            ),
            (
                "magma_transport",
                (; magma_transport=Erebus.MagmaTransportConfig(; active=true)),
            ),
            (
                "magma_degassing",
                (; magma_degassing=Erebus.MagmaOceanDegassingConfig(; active=true)),
            ),
            (
                "percolation_active",
                (; coreformation=Erebus.CoreFormationConfig(; percolation_active=true)),
            ),
            (
                "settling_active",
                (; coreformation=Erebus.CoreFormationConfig(; settling_active=true)),
            ),
            ("hr_fe", (; thermodynamics=Erebus.ThermalConfig(; hr_fe=true))),
            (
                "metal_partition",
                (; metal_partition=Erebus.MetalPartitionConfig(; active=true)),
            ),
            (
                "phase_tracking",
                (; phase_tracking=Erebus.PhaseTrackingConfig(; active=true)),
            ),
            (
                "volatile_mixture",
                (; volatile_mixture=Erebus.VolatileMixtureConfig(; active=true)),
            ),
            ("redox", (; redox=Erebus.RedoxConfig(; active=true))),
        )
            e = check(kw)
            @test e isa ArgumentError
            @test occursin("henke2012", e.msg) && occursin(name, e.msg)
        end
        # Model none keeps the hydrous default and every process
        none = Erebus.SimulationConfig(;
            thermodynamics=Erebus.ThermalConfig(; hr_fe=true),
            volatiles=Erebus.VolatilesConfig(; active=true),
        )
        @test Erebus._validate_porosity_config(none) === nothing
        @test isapprox(none.materials.XWsolidm_init[1], 0.5; rtol=1e-14)
    end

    @testset "init_simulation checks [porosity] for a keyword-built configuration" begin
        out = Erebus.OutputConfig(;
            mode=:telemetry, save_final=false, output_dir=tempname()
        )
        dry = Erebus.MaterialConfig(; XWsolidm_init=[0.0, 0.0, NaN])
        # henke2012 with the default reaction switch on
        cfg_r = Erebus.SimulationConfig(;
            porosity=Erebus.PorosityConfig(; model=:henke2012), materials=dry, output=out
        )
        e = _hk_error(() -> Erebus.init_simulation(cfg_r; output_path=mktempdir()))
        @test e isa ArgumentError
        @test occursin("reaction", e.msg)
        # An unknown model symbol does not run as model none
        cfg_u = Erebus.SimulationConfig(;
            porosity=Erebus.PorosityConfig(; model=:henke_2012),
            reaction=Erebus.ReactionConfig(; active=false),
            output=out,
        )
        e_u = _hk_error(() -> Erebus.init_simulation(cfg_u; output_path=mktempdir()))
        @test e_u isa ArgumentError
        @test occursin("porosity.model", e_u.msg)
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

    @testset "The marker step checks the [porosity] fields of a direct PorosityConfig" begin
        coords = Erebus.GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        pr = fill(1.0e5 + 2.0e5, coords.Ny1, coords.Nx1)
        function step!(por)
            pg = Erebus.PorosityGroup([0.5], [0.75], [0.0], [0.0])
            e = _hk_error(
                () -> Erebus.advance_marker_dry_porosity!(
                    pg,
                    [70_000.0],
                    [70_000.0],
                    [1],
                    [900.0],
                    pr,
                    1.0e10;
                    coords=coords,
                    psurface=1.0e5,
                    porosity=por,
                ),
            )
            return e, pg
        end
        # Control: valid constants sinter the hot marker below 0.5
        e_ok, pg_ok = step!(Erebus.PorosityConfig(; model=:henke2012))
        @test e_ok === nothing
        @test 0.0 < pg_ok.phi_d[1] < 0.49
        # NaN T_sinter would skip sintering silently (T >= NaN is false)
        for kw in (
            (; T_sinter=NaN),
            (; T_sinter=-10.0),
            (; T_sinter=0.0),
            (; grain_radius=0.0),
            (; grain_radius=NaN),
            (; A=-1.0),
            (; A=Inf),
            (; E_act=-85000.0),
            (; n=0.0),
            (; phi_cap=1.0),
            (; phi_init=1.2),
            (; model=:henke2013),
        )
            e, pg = step!(Erebus.PorosityConfig(; model=:henke2012, kw...))
            @test e isa ArgumentError
            @test occursin("porosity.$(only(keys(kw)))", e.msg)
            # The check runs before any marker changes
            @test isapprox(pg.phi_d[1], 0.5; rtol=0, atol=0)
        end
    end

    @testset "henke2012 needs poroelasticity.phimin <= 1e-3" begin
        cfg_of(phimin) = Erebus.SimulationConfig(;
            porosity=Erebus.PorosityConfig(; model=:henke2012),
            reaction=Erebus.ReactionConfig(; active=false),
            materials=Erebus.MaterialConfig(; XWsolidm_init=[0.0, 0.0, NaN]),
            poroelasticity=Erebus.PoroelasticConfig(; phimin=phimin),
        )
        @test Erebus.validate_config(cfg_of(1.0e-4)) === nothing
        @test Erebus.validate_config(cfg_of(1.0e-3)) === nothing
        for bad in (nextfloat(1.0e-3), 2.0e-3, 0.3)
            e = _hk_error(() -> Erebus.validate_config(cfg_of(bad)))
            @test e isa ArgumentError
            @test occursin("poroelasticity.phimin", e.msg)
        end
        # Model none keeps any phimin of the poroelastic section
        none_03 = Erebus.SimulationConfig(;
            poroelasticity=Erebus.PoroelasticConfig(; phimin=0.3)
        )
        @test Erebus.validate_config(none_03) === nothing
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
        @test length(m_h.phi0) == 64
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
        # The 2-D Poisson gravity of the Stokes solve has the source (8/3) pi G rho, so a
        # uniform cylinder has the central pressure of a uniform sphere
        for (mode, c_grav) in ((:enclosed_mass, 2 / 3), (:poisson2d, 2 / 3))
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
            # Scale guard: ~18.6 bar at the centre of 50 km
            @test 1.0e6 < P[1] < 4.0e6
        end
        # phi_cap = 0.3 below phi_init = 0.75 and phi_H(P) >= 0.42: phi_0 = 0.3 everywhere,
        # the same uniform limit; without the cap phi_0 >= 0.42 lowers P[1] by about 30 %
        for mode in (:enclosed_mass, :poisson2d)
            cfg_c = _hk_cfg("porosity.phi_cap" => 0.3, "geometry.gravity_mode" => mode)
            r_c, P_c = Erebus.dry_hydrostatic_profile(cfg_c)
            R = cfg_c.geometry.rplanet
            P_an = (2 / 3) * pi * G * (3300.0 * 0.7)^2 * R^2
            @test isapprox(P_c[1], P_an; rtol=1e-6)
            @test isapprox(P_c[end], 0.0; atol=1e-9 * P_an)
        end
        # With phi_init = 0.75 the deep rock is cold pressed: denser, higher pressure
        r, P = Erebus.dry_hydrostatic_profile(_hk_cfg())
        rho_top = 3300.0 * 0.25
        P_top_only = (2 / 3) * pi * G * rho_top^2 * r[end]^2
        @test P[1] > 2.0 * P_top_only
        # The profile solves dP/dr = -rho(P) g with rho(P) = 3300 (1 - phi_0(P)): the
        # residual with an independent trapezoid mass integral is 1e-4 of max(rho g) (at
        # the centre cell); 2 or 5 fixed-point iterations leave 0.30 or 3.4e-3
        rho = [3300.0 * (1 - Erebus.porosity_phi0(p * 1.0e-5, 0.75)) for p in P]
        mass = zeros(length(r))
        for i in 2:length(r)
            mass[i] =
                mass[i - 1] +
                2pi * (rho[i - 1] * r[i - 1]^2 + rho[i] * r[i]^2) * (r[i] - r[i - 1])
        end
        g = [i == 1 ? 0.0 : G * mass[i] / r[i]^2 for i in eachindex(r)]
        rhog = rho .* g
        res = [
            (P[i + 1] - P[i]) / (r[i + 1] - r[i]) + 0.5 * (rhog[i] + rhog[i + 1]) for
            i in 1:(length(r) - 1)
        ]
        @test maximum(abs.(res)) < 1.0e-3 * maximum(rhog)
        # Scale guard: the deep rock is pressed to rho = 3300 (1 - 0.42 ...) > 1900 kg/m^3
        @test 1900.0 < rho[1] < 3300.0 * 0.58
        # Depth where P reaches P_0 = 0.13 bar, for R_p = 50 km and 20 km (known limit)
        for (Rp, depth_ref) in ((50_000.0, 564.0), (20_000.0, 1667.0))
            cfg_d = _hk_cfg(
                "geometry.rplanet" => Rp,
                "geometry.rcrust" => Rp,
                "grid.xsize" => 2.8 * Rp,
                "grid.ysize" => 2.8 * Rp,
                "geometry.xcenter" => 1.4 * Rp,
                "geometry.ycenter" => 1.4 * Rp,
            )
            rd, Pd = Erebus.dry_hydrostatic_profile(cfg_d)
            k = findfirst(<(0.13e5), Pd)
            r0 =
                rd[k - 1] + (Pd[k - 1] - 0.13e5) / (Pd[k - 1] - Pd[k]) * (rd[k] - rd[k - 1])
            @test isapprox(Rp - r0, depth_ref; rtol=2.0e-3)
        end
        @test_throws ArgumentError Erebus.dry_hydrostatic_profile(_hk_cfg(); nr=1)
        # The 2-D Poisson profile of compressible powder (phi_init = 0.75) solves the
        # cylinder balance, not the sphere balance
        cfg_p = _hk_cfg("geometry.gravity_mode" => :poisson2d)
        r_p, P_p = Erebus.dry_hydrostatic_profile(cfg_p)
        @test _hk_hydro_residual(r_p, P_p, 0.75, 1) < 1.0e-6
        # Guard: the sphere balance leaves a residual of 4.7e-3 on this profile
        @test 1.0e-3 < _hk_hydro_residual(r_p, P_p, 0.75, 2) < 1.0e-2
        @test 1.0e6 < P_p[1] < 1.5e6
        # A gravity mode that is neither :enclosed_mass nor :poisson2d is an error
        cfg_g = _hk_with(
            _hk_cfg(); geometry=_hk_with(_hk_cfg().geometry; gravity_mode=:invalid)
        )
        e_g = _hk_error(() -> Erebus.dry_hydrostatic_profile(cfg_g))
        @test e_g isa ArgumentError
        @test occursin("gravity_mode", e_g.msg)
    end

    @testset "Setup pressure equals the Stokes pressure of the first step" begin
        # Uniform phi_0 = 0.3; centre pressure of the setup over pr - psurface after 1 step
        for (mode, metric) in ((:poisson2d, false), (:enclosed_mass, true))
            cfg = _hk_cfg(
                "porosity.phi_init" => 0.3,
                "geometry.gravity_mode" => mode,
                "geometry.spherical_metric" => metric,
                "time.n_steps" => 1,
                "time.dt_initial" => 1.0,
                "time.dt_longest" => 1.0,
            )
            _, P = Erebus.dry_hydrostatic_profile(cfg)
            st = Erebus.with_logger(Erebus.NullLogger()) do
                Erebus.simulation_loop(cfg; output_path=mktempdir())
            end
            coords = Erebus.GridCoordinates(cfg.grid)
            i, j, w = Erebus.fix_weights(
                cfg.geometry.xcenter,
                cfg.geometry.ycenter,
                coords.xp,
                coords.yp,
                coords.dx,
                coords.dy,
                coords.jmin_p,
                coords.jmax_p,
                coords.imin_p,
                coords.imax_p,
            )
            P_stokes =
                Erebus.dot4(Erebus.grid_vector(i, j, st.grids.pr), w) -
                cfg.geometry.psurface
            # Setup and Stokes centre pressures agree to 5 %
            @test isapprox(P[1] / P_stokes, 1.0; rtol=0.05)
            # Gravity guard: the cylinder gravity 4 pi G would give a ratio of 1.5
            @test !isapprox(P[1] / P_stokes, 1.5; rtol=0.2)
        end
    end

    @testset "Setup: water porosity at phimin, dry porosity cold pressed" begin
        cfg = _hk_cfg()
        st, coords, ws = _hk_state(cfg)
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
        # Adiabatic heating weights the solid term with 1 - PHI of the grid, PHI = phimin:
        # at the centre 1 - phi_d = 0.58, so the term is 1.72 times its dry value
        Erebus.interpolate_markers_to_grid!(
            st,
            coords,
            cfg;
            thread_buffers=ws.thread_buffers,
            interp_arrays=ws.interp_arrays,
        )
        ic = argmin(abs.(coords.xp .- xc))
        jc = argmin(abs.(coords.yp .- yc))
        dense = st.grids.RHO[jc, ic] / 3300.0
        @test isapprox(st.grids.PHI[jc, ic], cfg.poroelasticity.phimin; rtol=1e-9)
        @test isapprox(dense, 0.58; rtol=2e-3)
        @test isapprox((1 - st.grids.PHI[jc, ic]) / dense, 1.72; rtol=5e-3)
        # The porosity before sintering starts at the cold pressed porosity
        @test isapprox(mk.phi0[rock], mk.phi_d[rock]; rtol=1e-14)
        @test all(iszero, mk.phi0[air])
        # A cap below phi_init bounds the setup porosity: outer rock at 0.5, not 0.75
        st_c, _, _ = _hk_state(_hk_cfg("porosity.phi_cap" => 0.5))
        rock_c = findall(<(3), st_c.markers.tm)
        @test isapprox(maximum(st_c.markers.phi_d[rock_c]), 0.5; rtol=1e-12)
        @test isapprox(maximum(st_c.markers.phi0[rock_c]), 0.5; rtol=1e-12)
    end

    @testset "Setup: melting pressure pfm0 is psurface plus the dry profile" begin
        for mode in (:enclosed_mass, :poisson2d)
            cfg = _hk_cfg(
                "geometry.gravity_mode" => mode,
                "geometry.spherical_metric" => (mode === :enclosed_mass),
            )
            st, _, _ = _hk_state(cfg)
            mk = st.markers
            rock = findall(<(3), mk.tm)
            r_prof, P_prof = Erebus.dry_hydrostatic_profile(cfg)
            xc, yc = cfg.geometry.xcenter, cfg.geometry.ycenter
            ps = cfg.geometry.psurface
            P_exp = [
                ps +
                Erebus._interp_profile(r_prof, P_prof, hypot(mk.xm[m] - xc, mk.ym[m] - yc))
                for m in rock
            ]
            @test isapprox(mk.pfm0[rock], P_exp; rtol=1e-12)
            # Scale guard: about 12.5 bar at the centre, not the 1.56 or 2.34 times larger
            # value of the water-porosity estimate
            @test 1.0e6 < maximum(mk.pfm0[rock]) - ps < 1.4e6
            @test all(iszero, mk.pfm0[findall(==(3), mk.tm)])
        end
    end

    @testset "Melting pressure pfm0 follows the Stokes total pressure after step 1" begin
        # Compressible pores (production beta values) make pf differ from pr after step 1
        for n_steps in (2, 3)
            cfg = _hk_cfg(
                "time.n_steps" => n_steps,
                "melting.active" => true,
                "time.dt_initial" => 100.0,
                "time.dt_longest" => 100.0,
                "poroelasticity.betasolid" => 2.5e-11,
                "poroelasticity.betafluid" => 4.0e-10,
            )
            _, coords, _ = _hk_state(cfg)
            st = Erebus.with_logger(Erebus.NullLogger()) do
                Erebus.simulation_loop(cfg; output_path=mktempdir())
            end
            mk = st.markers
            rock = findall(<(3), mk.tm)
            xc, yc = cfg.geometry.xcenter, cfg.geometry.ycenter
            m = rock[argmin([hypot(mk.xm[k] - xc, mk.ym[k] - yc) for k in rock])]
            i, j, w = Erebus.fix_weights(
                mk.xm[m],
                mk.ym[m],
                coords.xp,
                coords.yp,
                coords.dx,
                coords.dy,
                coords.jmin_p,
                coords.jmax_p,
                coords.imin_p,
                coords.imax_p,
            )
            pr_c = Erebus.dot4(Erebus.grid_vector(i, j, st.grids.pr), w)
            pf_c = Erebus.dot4(Erebus.grid_vector(i, j, st.grids.pf), w)
            @test st.timestep == n_steps
            @test isapprox(mk.pfm0[m], pr_c; rtol=1e-6)
            # Guard: the fluid pressure is 0.5 % to 1 % lower here, so pf would fail
            @test 1.0e-3 < (pr_c - pf_c) / pr_c < 5.0e-2
            # Scale guard: about 12.5 bar of overburden at the centre
            @test 1.0e6 < mk.pfm0[m] - cfg.geometry.psurface < 1.4e6
        end
    end

    @testset "A direct solve_stokes_darcy! call sets pfm0 from the total pressure" begin
        # The model step overwrites pfm0 after this call, so call it alone as a modular
        # step at step 2 with compressible pores, where pf differs from pr
        cfg = _hk_cfg(
            "poroelasticity.betasolid" => 2.5e-11, "poroelasticity.betafluid" => 4.0e-10
        )
        st, coords, ws = _hk_state(cfg)
        st.timestep = 2
        st.dt = 100.0 * _HK_SECONDS_PER_YEAR
        Erebus.prepare_step_ambient!(st, coords, cfg, ws)
        Erebus.radiogenic_heating!(st, coords, cfg)
        Erebus.interpolate_markers_to_grid!(
            st,
            coords,
            cfg;
            p2m_workspace=ws.p2m,
            thread_buffers=ws.thread_buffers,
            interp_arrays=ws.interp_arrays,
        )
        Erebus.solve_gravity!(st, coords, cfg; F_grav=ws.F_grav, RP=ws.RP, SP=ws.SP)
        Erebus.update_step_start_pressures!(st)
        Erebus.solve_stokes_darcy!(st, coords, cfg, ws; titer=1, dt_step_initial=st.dt)
        mk = st.markers
        rock = findall(<(3), mk.tm)
        xc, yc = cfg.geometry.xcenter, cfg.geometry.ycenter
        m = rock[argmin([hypot(mk.xm[k] - xc, mk.ym[k] - yc) for k in rock])]
        i, j, w = Erebus.fix_weights(
            mk.xm[m],
            mk.ym[m],
            coords.xp,
            coords.yp,
            coords.dx,
            coords.dy,
            coords.jmin_p,
            coords.jmax_p,
            coords.imin_p,
            coords.imax_p,
        )
        pr_c = Erebus.dot4(Erebus.grid_vector(i, j, st.grids.pr), w)
        pf_c = Erebus.dot4(Erebus.grid_vector(i, j, st.grids.pf), w)
        @test isapprox(mk.pfm0[m], pr_c; rtol=1e-6)
        # Guard: the fluid pressure is about 60 % lower here, so pf would fail
        @test 0.1 < (pr_c - pf_c) / pr_c < 0.9
        # Scale guard: about 12 bar of overburden at the centre
        @test 1.0e6 < mk.pfm0[m] - cfg.geometry.psurface < 1.4e6
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
        # Cold pressing and sintering in 1 step: dphi_sinter counts only the sintering
        # below phi_0 = phi_H(1 bar) = 0.4334, not the pressing from 0.75
        cs = Erebus.dry_porosity_update(0.75, 1.0, 900.0, dt, pc)
        @test isapprox(cs.phi0, Erebus.phi_cold(1.0); rtol=1e-14)
        @test isapprox(cs.dphi_sinter, cs.phi0 - cs.phi; rtol=1e-12)
        @test 0.05 < cs.dphi_sinter < 0.2
        # The cap of the configuration bounds phi_0: phi_cap 0.5 < phi_init 0.75 at P = 0
        pc_cap = _hk_cfg("porosity.phi_cap" => 0.5).porosity
        cap = Erebus.dry_porosity_update(0.75, 0.0, 300.0, dt, pc_cap)
        @test isapprox(cap.phi, 0.5; rtol=1e-14)
        @test isapprox(cap.phi0, 0.5; rtol=1e-14)
        @test_throws DomainError Erebus.dry_porosity_update(0.5, -1.0, 900.0, dt, pc)
        @test_throws DomainError Erebus.dry_porosity_update(0.5, 1.0, 0.0, dt, pc)
        @test_throws DomainError Erebus.dry_porosity_update(0.5, 1.0, -10.0, dt, pc)
        @test_throws DomainError Erebus.dry_porosity_update(1.1, 1.0, 900.0, dt, pc)
        @test_throws DomainError Erebus.dry_porosity_update(-0.1, 1.0, 900.0, dt, pc)
        @test_throws DomainError Erebus.dry_porosity_update(
            0.5, 1.0, 900.0, dt, pc; phi0_prev=1.5
        )
        # Below T_sinter no later check sees phi0_prev, so the bound itself must reject it
        @test_throws DomainError Erebus.dry_porosity_update(
            0.5, 1.0, 500.0, dt, pc; phi0_prev=-0.1
        )
    end

    @testset "phi_0 remembers the highest pressure the marker has felt" begin
        pc = _hk_cfg().porosity
        # Cold pressing at 1 bar and 500 K: phi_d = phi_0 = phi_H(1 bar) = 0.4334
        u1 = Erebus.dry_porosity_update(0.75, 1.0, 500.0, 1.0e10, pc)
        @test isapprox(u1.phi, 0.4334; atol=1e-4)
        @test isapprox(u1.phi0, u1.phi; rtol=1e-14)
        # Then 0.2 bar and 900 K for 100 yr: a pressure drop never raises phi_0
        dt = 100.0 * _HK_SECONDS_PER_YEAR
        u2 = Erebus.dry_porosity_update(u1.phi, 0.2, 900.0, dt, pc; phi0_prev=u1.phi0)
        @test isapprox(u2.phi0, u1.phi0; rtol=1e-14)
        @test u2.phi0 < Erebus.porosity_phi0(0.2, 0.75)
        # Same sintering as a marker that started at that phi_0 (X = 1 at the onset)
        phi_ref, _, _ = Erebus.sinter_porosity_step(u1.phi, u1.phi, 0.2, 900.0, dt)
        @test isapprox(u2.dphi_sinter, u1.phi - phi_ref; rtol=1e-12)
        @test isapprox(u2.phi, phi_ref; rtol=1e-14)
        # Guard: phi_0 of the current pressure would treat the pressed rock as sintered
        # and give a 79 times smaller decrease
        u_now = Erebus.dry_porosity_update(u1.phi, 0.2, 900.0, dt, pc)
        @test 50.0 < u2.dphi_sinter / u_now.dphi_sinter < 120.0
        # Before sintering, a pressure rise lowers phi_0 to phi_0(P) and cold presses
        u3 = Erebus.dry_porosity_update(u1.phi, 50.0, 500.0, dt, pc; phi0_prev=u1.phi0)
        @test isapprox(u3.phi0, Erebus.porosity_phi0(50.0, 0.75); rtol=1e-14)
        @test isapprox(u3.phi, u3.phi0; rtol=1e-14)
        @test iszero(u3.dphi_sinter)
    end

    @testset "phi_0 freezes at the onset of sintering (Henke Eq. 37)" begin
        pc = _hk_cfg().porosity
        dt = 5000.0 * _HK_SECONDS_PER_YEAR
        # Sintering at 0.15 bar and 950 K from phi_init 0.75: phi_0 = phi_H(0.15 bar)
        a1 = Erebus.dry_porosity_update(0.75, 0.15, 950.0, dt, pc)
        @test isapprox(a1.phi0, Erebus.phi_cold(0.15); rtol=1e-14)
        @test 0.05 < a1.dphi_sinter < 0.2
        # Then 0.25 bar and 850 K: phi_0(0.25 bar) = 0.533 lies between phi and phi_0
        p0_25 = Erebus.porosity_phi0(0.25, 0.75)
        @test a1.phi < p0_25 < a1.phi0
        a2 = Erebus.dry_porosity_update(
            a1.phi, 0.25, 850.0, dt, pc; phi0_prev=a1.phi0, sintered=true
        )
        # The step keeps phi_0 of 0.15 bar: same decrease as sintering from that phi_0
        phi_ref, _, _ = Erebus.sinter_porosity_step(a1.phi, a1.phi0, 0.25, 850.0, dt)
        @test isapprox(a2.phi0, a1.phi0; rtol=1e-14)
        @test isapprox(a2.dphi_sinter, a1.phi - phi_ref; rtol=1e-12)
        @test isapprox(a2.phi, phi_ref; rtol=1e-14)
        # Guard: phi_0 lowered to 0.533 shrinks the contact faces, a 9 times larger step
        phi_low, _, _ = Erebus.sinter_porosity_step(a1.phi, p0_25, 0.25, 850.0, dt)
        @test 5.0 < (a1.phi - phi_low) / a2.dphi_sinter < 15.0
        # Sintering at 0.05 bar, then 5 bar: no cold pressing of the sintered rock, and
        # the cumulative decrease by sintering stays phi_0 - phi_d < phi_0
        dtb = 5.0e4 * _HK_SECONDS_PER_YEAR
        b1 = Erebus.dry_porosity_update(0.75, 0.05, 950.0, dtb, pc)
        b2 = Erebus.dry_porosity_update(
            b1.phi,
            5.0,
            950.0,
            dtb,
            pc;
            phi0_prev=b1.phi0,
            sintered=b1.dphi_sinter > Erebus.DPHI_SINTERED,
        )
        @test isapprox(b2.phi0, b1.phi0; rtol=1e-14)
        @test isapprox(b1.dphi_sinter + b2.dphi_sinter, b2.phi0 - b2.phi; atol=1e-12)
        @test 0.0 < b2.phi < b1.phi
        @test b1.dphi_sinter + b2.dphi_sinter < b2.phi0
        # Edge case: cold and at 50 bar, the sintered rock keeps phi and phi_0
        c = Erebus.dry_porosity_update(
            a1.phi, 50.0, 500.0, dt, pc; phi0_prev=a1.phi0, sintered=true
        )
        @test isapprox(c.phi, a1.phi; rtol=1e-14)
        @test isapprox(c.phi0, a1.phi0; rtol=1e-14)
        @test iszero(c.dphi_sinter)
        # Without the flag the marker counts as unsintered and is cold pressed at 50 bar
        c0 = Erebus.dry_porosity_update(a1.phi, 50.0, 500.0, dt, pc; phi0_prev=a1.phi0)
        @test isapprox(c0.phi, Erebus.porosity_phi0(50.0, 0.75); rtol=1e-14)
        # A sintered marker has phi_d <= phi_0 and needs a pressure >= 0
        @test_throws DomainError Erebus.dry_porosity_update(
            0.6, 1.0, 900.0, dt, pc; phi0_prev=0.5, sintered=true
        )
        @test_throws DomainError Erebus.dry_porosity_update(
            0.4, -1.0, 500.0, dt, pc; phi0_prev=0.5, sintered=true
        )
        # Edge case phi_d = phi_0: the contract accepts a sintered marker at its phi_0;
        # it keeps phi when cold and sinters when hot
        e0 = Erebus.dry_porosity_update(
            a1.phi0, 50.0, 500.0, dt, pc; phi0_prev=a1.phi0, sintered=true
        )
        @test isapprox(e0.phi, a1.phi0; rtol=1e-14)
        @test isapprox(e0.phi0, a1.phi0; rtol=1e-14)
        e1 = Erebus.dry_porosity_update(
            a1.phi0, 0.15, 950.0, dt, pc; phi0_prev=a1.phi0, sintered=true
        )
        @test isapprox(e1.phi, a1.phi; rtol=1e-12)
        @test isapprox(e1.dphi_sinter, a1.dphi_sinter; rtol=1e-12)
    end

    @testset "dry_porosity_update checks dt and phi0_prev below and above T_sinter" begin
        pc = _hk_cfg().porosity
        # 500 K takes the cold branch, 900 K the sintering branch
        for T in (500.0, 900.0)
            for bad_dt in (-10.0, NaN, Inf)
                @test_throws DomainError Erebus.dry_porosity_update(0.5, 1.0, T, bad_dt, pc)
            end
            # A sintered marker has phi_0 < 1: 1 means no pressure history
            e = _hk_error(
                () -> Erebus.dry_porosity_update(
                    0.5, 1.0, T, 1.0e10, pc; phi0_prev=1.0, sintered=true
                ),
            )
            @test e isa DomainError
            @test occursin("phi0_prev", e.msg)
            # A BigFloat phi0_prev below 1 that rounds to 1.0 in Float64 is rejected too
            e_big = _hk_error(
                () -> Erebus.dry_porosity_update(
                    0.5, 1.0, T, 1.0e10, pc; phi0_prev=1 - big(1.0e-30), sintered=true
                ),
            )
            @test e_big isa DomainError
            @test occursin("phi0_prev", e_big.msg)
            # dt = 0 is a valid step: cold pressing only, no sintering
            z = Erebus.dry_porosity_update(0.75, 1.0, T, 0.0, pc)
            @test isapprox(z.phi, Erebus.phi_cold(1.0); rtol=1e-14)
            @test iszero(z.dphi_sinter)
        end
    end

    @testset "dry_porosity_update uses T_sinter, grain_radius, E_act and n of the config" begin
        pc0 = _hk_cfg().porosity
        step(pc) = Erebus.dry_porosity_update(
            0.5, 2.0, 750.0, 1.0e10, pc; phi0_prev=0.5, sintered=true
        )
        d0 = step(pc0).dphi_sinter
        # Scale guard: 750 K and 2 bar sinter 1.4e-3 in 1e10 s with the default constants
        @test 1.0e-3 < d0 < 2.0e-3
        # T_sinter = 760 K stops sintering at 750 K; T_sinter = 750 K (T >= T_sinter) does not
        u_hi = step(_hk_with(pc0; T_sinter=760.0))
        @test iszero(u_hi.dphi_sinter)
        @test isapprox(u_hi.phi, 0.5; rtol=0, atol=0)
        @test isapprox(step(_hk_with(pc0; T_sinter=750.0)).dphi_sinter, d0; rtol=0, atol=0)
        # Larger grains and a higher activation energy sinter slower; a higher stress
        # exponent sinters faster, since the contact stress exceeds 1 bar
        for (kw, step_kw, slower) in (
            ((; grain_radius=1.0e-5), (; R0=1.0e-5), true),
            ((; E_act=90000.0), (; E_act=90000.0), true),
            ((; n=2.0), (; n=2.0), false),
        )
            u = step(_hk_with(pc0; kw...))
            phi_ref, _, _ = Erebus.sinter_porosity_step(
                0.5, 0.5, 2.0, 750.0, 1.0e10; step_kw...
            )
            @test isapprox(u.phi, phi_ref; rtol=0, atol=0)
            @test slower ? (u.dphi_sinter < 0.5 * d0) : (u.dphi_sinter > 2.0 * d0)
        end
    end

    @testset "A porosity decrease below DPHI_SINTERED does not count as sintering" begin
        pc = _hk_cfg().porosity
        # The threshold lies far above the rounding of 1 - phi near X = 1 (a few 1e-16)
        @test isapprox(Erebus.DPHI_SINTERED, 1.0e-12; rtol=0, atol=0)
        @test 1.0e3 * eps(1.0) < Erebus.DPHI_SINTERED < 1.0e-9
        # A pushed marker (phi_d = 0.6 < phi_0 = 0.75) at 1e-12 bar for 1 s: h r(0.6)
        # is 1.5e-27 ULP, so the backward Euler root rounds to 0.6 itself
        u = Erebus.dry_porosity_update(0.6, 1.0e-12, 700.0, 1.0, pc)
        @test isapprox(u.phi, 0.6; rtol=0, atol=0)
        @test iszero(u.dphi_sinter)
        # The marker is still unsintered: a later 50 bar, 500 K step cold presses it
        u2 = Erebus.dry_porosity_update(
            u.phi,
            50.0,
            500.0,
            1.0,
            pc;
            phi0_prev=u.phi0,
            sintered=u.dphi_sinter > Erebus.DPHI_SINTERED,
        )
        @test isapprox(u2.phi, Erebus.porosity_phi0(50.0, 0.75); rtol=1e-14)
        @test u2.phi < 0.45
        # The root never over-sinters: phi_old - phi <= h r(phi), and the next float
        # below phi does (the root is bracketed to 1 ULP)
        for (phi_old, h, P_bar) in (
            (0.6, 1.0, 1.0e-12),
            (prevfloat(0.6), 1.0, 1.0e-12),
            (0.7, 1.0e9, 1.0),
            (0.3, 1.0e11, 1.0),
        )
            rate(p) = Erebus.porosity_decrease_rate(p, 0.75, P_bar, 950.0)
            phi = Erebus._backward_euler_porosity(phi_old, h, rate)
            @test phi_old - phi <= h * rate(phi)
            @test phi_old - prevfloat(phi) > h * rate(prevfloat(phi))
        end
        # Edge case at the onset phi_d = phi_0 = 0.75: the exact backward Euler decrease
        # is 5.69e-16 = 5.1 ULP (256-bit bisection of the reference rate), far below
        # DPHI_SINTERED, so it does not count as sintering
        v = Erebus.dry_porosity_update(0.75, 1.0e-12, 700.0, 1.0, pc)
        @test 3.0 <= v.dphi_sinter / eps(0.75) <= 6.0
        @test isapprox(v.dphi_sinter, v.phi0 - v.phi; rtol=0, atol=0)
        # phi_0 = 0.1 and 0.42 at 1e-30 bar: the rounding of 1 - phi near X = 1 lowers
        # phi by 14 and 2 ULP (exact decrease 3e-26); both stay below the threshold
        for p0 in (0.1, 0.42)
            pc0 = _hk_cfg("porosity.phi_init" => p0).porosity
            w = Erebus.dry_porosity_update(p0, 1.0e-30, 700.0, 1.0, pc0)
            @test 0.0 <= w.dphi_sinter <= 20 * eps(p0)
            @test w.dphi_sinter <= Erebus.DPHI_SINTERED
        end
        # The marker step leaves the flag unset: a later 50 bar, 500 K step cold presses
        # the marker of the onset case from 0.75 to phi_H(50 bar) = 0.42
        cfg = _hk_cfg()
        coords = Erebus.GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        mk = Erebus.init_marker_arrays(0, cfg, coords)
        Erebus.push_marker!(mk; xm=70_000.0, ym=70_000.0, tm=1, tkm=700.0, phi_d=0.75)
        ps = cfg.geometry.psurface
        for (P_pa, T) in ((1.0e-7, 700.0), (50.0e5, 500.0))
            mk.tkm[1] = T
            Erebus.advance_marker_dry_porosity!(
                mk.groups.porosity,
                mk.xm,
                mk.ym,
                mk.tm,
                mk.tkm,
                fill(ps + P_pa, coords.Ny1, coords.Nx1),
                1.0;
                coords=coords,
                psurface=ps,
                porosity=cfg.porosity,
            )
        end
        @test 0.0 < mk.dphi_sinter[1] <= 10 * eps(0.75)
        @test isapprox(mk.phi0[1], Erebus.porosity_phi0(50.0, 0.75); rtol=1e-14)
        @test isapprox(mk.phi_d[1], mk.phi0[1]; rtol=1e-14)
        # Guard: a frozen phi_0 would stay at 0.75
        @test mk.phi0[1] < 0.43
        # The invariant dphi_sinter = phi0 - phi_d then holds to within the threshold
        mismatch = mk.dphi_sinter[1] - (mk.phi0[1] - mk.phi_d[1])
        @test abs(mismatch) <= Erebus.DPHI_SINTERED
    end

    @testset "Marker sintered just above DPHI_SINTERED keeps phi_0 at 50 bar" begin
        cfg = _hk_cfg()
        coords = Erebus.GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        mk = Erebus.init_marker_arrays(0, cfg, coords)
        Erebus.push_marker!(mk; xm=70_000.0, ym=70_000.0, tm=1, tkm=700.0, phi_d=0.75)
        ps = cfg.geometry.psurface
        for (P_pa, T, dt) in ((1.0e-3, 700.0, 1.0e4), (50.0e5, 500.0, 1.0))
            mk.tkm[1] = T
            Erebus.advance_marker_dry_porosity!(
                mk.groups.porosity,
                mk.xm,
                mk.ym,
                mk.tm,
                mk.tkm,
                fill(ps + P_pa, coords.Ny1, coords.Nx1),
                dt;
                coords=coords,
                psurface=ps,
                porosity=cfg.porosity,
            )
        end
        # Guard: the first step sinters about 5.7e-12, past the threshold, far below 1e-6
        @test Erebus.DPHI_SINTERED < mk.dphi_sinter[1] < 1.0e-9
        @test isapprox(mk.phi0[1], 0.75; rtol=0, atol=0)
        @test 0.74 < mk.phi_d[1] < 0.75
    end

    @testset "Ten-step run keeps dphi_sinter = phi0 - phi_d on every rock marker" begin
        # Densification raises g and the Stokes pressure under sintered rock each step
        cfg = _hk_cfg(
            "time.n_steps" => 10,
            "materials.tkm0" => [900.0, 900.0, 170.0],
            "time.dt_initial" => 1000.0,
            "time.dt_longest" => 5000.0,
        )
        st = Erebus.with_logger(Erebus.NullLogger()) do
            Erebus.simulation_loop(cfg; output_path=mktempdir())
        end
        mk = st.markers
        rock = findall(<(3), mk.tm)
        @test st.timestep == 10
        mismatch = mk.dphi_sinter[rock] .- (mk.phi0[rock] .- mk.phi_d[rock])
        @test maximum(abs.(mismatch)) < 1.0e-12
        # The run sinters: up to about 0.4 of porosity, never more than phi_0
        @test 0.1 < maximum(mk.dphi_sinter[rock]) < 0.75
        @test all(mk.dphi_sinter[rock] .<= mk.phi0[rock])
        @test all(0.0 .<= mk.phi_d[rock] .<= mk.phi0[rock])
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
        # A dilatant compaction rate (APHI < 0) would raise phim by 1 % in this step
        fill!(st.grids.APHI, -1.0e-8)
        P_pa = 2.0e5
        fill!(st.grids.pr, cfg.geometry.psurface + P_pa)
        fill!(st.grids.pf, 0.0)
        mk.tkm[rock] .= 1000.0
        phi_before = copy(mk.phi_d)
        phi0_before = copy(mk.phi0)
        Erebus.vent_and_degas!(st, coords, cfg)
        for m in rock
            exp_u = Erebus.dry_porosity_update(
                phi_before[m], 2.0, 1000.0, 1.0e6, cfg.porosity; phi0_prev=phi0_before[m]
            )
            @test isapprox(mk.phi_d[m], exp_u.phi; rtol=1e-12)
        end
        @test all(mk.phi_d[rock] .< phi_before[rock])
        @test all(isapprox.(mk.tkm_max[rock], 1000.0; rtol=1e-14))
        @test all(mk.dphi_sinter[rock] .> 0.0)
        # phi0 keeps the deeper setup pressure where it exceeds 2 bar
        phi0_2bar = Erebus.porosity_phi0(2.0, 0.75)
        @test all(isapprox.(mk.phi0[rock], min.(phi0_before[rock], phi0_2bar); rtol=1e-14))
        @test any(mk.phi0[rock] .< phi0_2bar - 1.0e-3)
        # Water porosity stays at phimin: the compaction rate APHI is not applied
        @test all(isapprox.(mk.phim[rock], cfg.poroelasticity.phimin; rtol=1e-14))
        @test all(iszero, mk.phi_d[air])
        # A second step adds its sintering to dphi_sinter (no cold pressing at fixed P)
        d1 = copy(mk.dphi_sinter)
        phi1 = copy(mk.phi_d)
        Erebus.vent_and_degas!(st, coords, cfg)
        @test all(mk.phi_d[rock] .< phi1[rock])
        @test all(
            isapprox.(
                mk.dphi_sinter[rock], d1[rock] .+ (phi1[rock] .- mk.phi_d[rock]); rtol=1e-12
            ),
        )
    end

    @testset "Substep cap: counted per marker and reported by vent_and_degas!" begin
        # An extreme creep prefactor and time step need more than 1000 substeps at every
        # start porosity from 0.42 to 0.75 (A = 1e250 with dt = 1e30 s caps only 0.75)
        A_x, dt_x = 1.0e250, 1.0e100
        _, nsub, capped = Erebus.sinter_porosity_step(0.75, 0.75, 2.0, 1000.0, dt_x; A=A_x)
        @test capped
        @test nsub == 1000
        # Edge case: the default constants over 1 kyr stay far below the cap
        _, nsub_ref, capped_ref = Erebus.sinter_porosity_step(
            0.75, 0.75, 2.0, 1000.0, 1.0e3 * _HK_SECONDS_PER_YEAR
        )
        @test !capped_ref
        @test 1 <= nsub_ref < 500
        cfg = _hk_cfg("porosity.A" => A_x)
        st, coords, _ = _hk_state(cfg)
        mk = st.markers
        # 5 hot rock markers sinter; the others stay below T_sinter
        hot = findall(<(3), mk.tm)[1:300:1500]
        @test all(mk.tkm .< 700.0)
        mk.tkm[hot] .= 1000.0
        pr = fill(cfg.geometry.psurface + 2.0e5, size(st.grids.pr))
        ncap = Erebus.advance_marker_dry_porosity!(
            mk.groups.porosity,
            mk.xm,
            mk.ym,
            mk.tm,
            mk.tkm,
            pr,
            dt_x;
            coords=coords,
            psurface=cfg.geometry.psurface,
            porosity=cfg.porosity,
        )
        @test ncap == length(hot)
        # The rock sinters to the smallest positive float, not 0: for phi <= 0.1 the rate
        # is proportional to phi, so the backward Euler root stays above 0
        @test all(p -> 0.0 <= p <= nextfloat(0.0), mk.phi_d[hot])
        # The step of the model warns with the number of capped markers
        st2, coords2, _ = _hk_state(cfg)
        hot2 = findall(<(3), st2.markers.tm)[1:300:1500]
        st2.markers.tkm[hot2] .= 1000.0
        st2.timestep = 2
        st2.dt = dt_x
        fill!(st2.grids.DT, 0.0)
        fill!(st2.grids.pr, cfg.geometry.psurface + 2.0e5)
        @test_logs (:warn, "Dry porosity: 5 markers reached the sintering substep cap") match_mode =
            :any Erebus.vent_and_degas!(st2, coords2, cfg)
        @test all(p -> 0.0 <= p <= nextfloat(0.0), st2.markers.phi_d[hot2])
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

        # The 3 particle-to-mesh paths give the same marker properties: tiled, buffered
        # without thread buffers (serial) and with them (threaded when nthreads > 1)
        for (mode, use_tb) in ((:tiled, true), (:buffered, false), (:buffered, true))
            cfg = _hk_cfg("solver.p2m_mode" => mode)
            st, coords, ws = _hk_state(cfg)
            mk = st.markers
            fill!(mk.rhototalm, 0.0)
            fill!(mk.ktotalm, 0.0)
            Erebus.interpolate_markers_to_grid!(
                st,
                coords,
                cfg;
                thread_buffers=use_tb ? ws.thread_buffers : nothing,
                interp_arrays=ws.interp_arrays,
            )
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
        # The step wrapper diffuse_and_update_markers! passes phi_d (dsubgridt = 1); at
        # step 1 the marker keeps the relaxed difference to the uniform 300 K grid
        cfg_s = _hk_cfg("solver.dsubgridt" => 1.0)
        st, coords, ws = _hk_state(cfg_s)
        mk = st.markers
        rock = findall(<(3), mk.tm)
        mk.phi_d[rock] .= 0.6
        mk.tkm[rock] .= 400.0
        fill!(st.grids.tk1, 300.0)
        fill!(st.grids.DT, 0.0)
        st.timestep = 1
        st.dt = 1.0e13
        Erebus.diffuse_and_update_markers!(st, coords, cfg_s, ws)
        k_d = Erebus.k_porous(
            Erebus.compute_ksolidm(400.0, Erebus.marker_property_mode), 0.6
        )
        x_d = k_d * 1.0e13 / (3.3e6 * 0.4) * 2 * (1 / coords.dx^2 + 1 / coords.dy^2)
        @test isapprox(
            mk.tkm[rock], fill(300.0 + 100.0 * exp(-x_d), length(rock)); rtol=1e-12
        )
        # Scale guard: the dry rock barely relaxes; wet rock at phimin gives exp(-x) ~ 0.6
        @test 0.9 < exp(-x_d) < 1.0

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
        face_y = findfirst(!isapprox(1.0), KYd)
        @test isapprox(KYd[face_y], 2kb * h * coords.dy / (2kb + h * coords.dy); rtol=1e-12)
        @test KYd[face_y] < KYs[face_y]
        # The step wrapper apply_surface_radiation! uses phi_0(P = 0) = 0.75 on both faces
        cfg_r = _hk_cfg("thermodynamics.surface_radiation" => true)
        st, coords, _ = _hk_state(cfg_r)
        fill!(st.grids.KX, 1.0)
        fill!(st.grids.KY, 1.0)
        fill!(st.grids.tk1, 300.0)
        Erebus.apply_surface_radiation!(st, coords, cfg_r)
        T_amb = Erebus.compute_ambient_conditions(st.timesum, cfg_r.disk)[1]
        h_r = Erebus.compute_effective_radiation_htc(
            300.0,
            T_amb,
            0.0;
            emissivity=cfg_r.thermodynamics.emissivity,
            sigma_sb=cfg_r.thermodynamics.sigma_sb,
        )
        kb_r = Erebus.k_porous(
            Erebus.compute_ksolidm(300.0, Erebus.marker_property_mode), 0.75
        )
        fx = findfirst(!isapprox(1.0), st.grids.KX)
        fy = findfirst(!isapprox(1.0), st.grids.KY)
        @test isapprox(
            st.grids.KX[fx], 2kb_r * h_r * coords.dx / (2kb_r + h_r * coords.dx); rtol=1e-12
        )
        @test isapprox(
            st.grids.KY[fy], 2kb_r * h_r * coords.dy / (2kb_r + h_r * coords.dy); rtol=1e-12
        )
        # phi_cap = 0.5 below phi_init = 0.75: the faces use phi_0(P = 0) = 0.5
        cfg_rc = _hk_cfg(
            "thermodynamics.surface_radiation" => true, "porosity.phi_cap" => 0.5
        )
        st, coords, _ = _hk_state(cfg_rc)
        fill!(st.grids.KX, 1.0)
        fill!(st.grids.KY, 1.0)
        fill!(st.grids.tk1, 300.0)
        Erebus.apply_surface_radiation!(st, coords, cfg_rc)
        kb_c = Erebus.k_porous(
            Erebus.compute_ksolidm(300.0, Erebus.marker_property_mode), 0.5
        )
        @test isapprox(
            st.grids.KX[fx], 2kb_c * h_r * coords.dx / (2kb_c + h_r * coords.dx); rtol=1e-12
        )
        @test isapprox(
            st.grids.KY[fy], 2kb_c * h_r * coords.dy / (2kb_c + h_r * coords.dy); rtol=1e-12
        )
        # Guard: the less porous powder conducts better, so the face value rises
        @test st.grids.KX[fx] > 1.01 * (2kb_r * h_r * coords.dx / (2kb_r + h_r * coords.dx))
    end

    @testset "Marker lifecycle: replenish and checkpoint keep the group" begin
        cfg = _hk_cfg()
        coords = Erebus.GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        markers = Erebus.init_marker_arrays(10, cfg, coords)
        # Spread sources; every group value is a distinct function of the core temperature
        markers.xm .= range(5_000.0, 135_000.0; length=10)
        markers.ym .= range(135_000.0, 5_000.0; length=10)
        markers.tkm .= range(300.0, 390.0; length=10)
        markers.phi_d .= 0.1 .+ (markers.tkm .- 300.0) ./ 200.0
        markers.phi0 .= markers.phi_d .+ 0.05
        markers.tkm_max .= markers.tkm .+ 100.0
        markers.dphi_sinter .= (markers.tkm .- 300.0) ./ 1.0e4
        mdis = zeros(coords.Nym, coords.Nxm)
        mnum = zeros(Int, coords.Nym, coords.Nxm)
        n_new = Erebus.replenish_markers!(markers, mdis, mnum; coords=coords, cfg=cfg)
        @test 10 < n_new < 100_000
        @test length(markers.phi_d) == n_new && length(markers.phi0) == n_new
        @test length(markers.tkm_max) == n_new && length(markers.dphi_sinter) == n_new
        # New markers copy the group values of the same source marker as the core values
        @test all(0.1 - 1e-12 .<= markers.phi_d .<= 0.55 + 1e-12)
        @test isapprox(markers.phi_d, 0.1 .+ (markers.tkm .- 300.0) ./ 200.0; rtol=1e-14)
        @test isapprox(markers.phi0, markers.phi_d .+ 0.05; rtol=1e-14)
        @test isapprox(markers.tkm_max, markers.tkm .+ 100.0; rtol=1e-14)
        @test isapprox(markers.dphi_sinter, (markers.tkm .- 300.0) ./ 1.0e4; atol=1e-16)
        @test length(unique(markers.tkm[11:end])) == 10
        # Checkpoint round trip keeps every non-uniform group array
        st, coords2, _ = _hk_state(cfg)
        nm = length(st.markers)
        st.markers.phi_d .= range(0.2, 0.7; length=nm)
        st.markers.phi0 .= st.markers.phi_d .+ 0.01
        st.markers.tkm_max .= range(650.0, 950.0; length=nm)
        st.markers.dphi_sinter .= range(0.0, 0.05; length=nm)
        saved = deepcopy(st.markers.groups.porosity)
        path = joinpath(mktempdir(), "ckpt.jld2")
        Erebus.save_state(path, st, coords2, cfg)
        data = Erebus.JLD2.load(path)
        restored = Erebus.restore_checkpoint_markers(data, st.markers.core, cfg)
        @test haskey(restored.groups, :porosity)
        for fn in fieldnames(Erebus.PorosityGroup)
            @test isapprox(
                getfield(restored.groups.porosity, fn), getfield(saved, fn); rtol=0, atol=0
            )
        end
        # A checkpoint without one of the group arrays is a CheckpointError
        for fn in fieldnames(Erebus.PorosityGroup)
            partial = copy(data)
            delete!(partial, string(fn))
            e = _hk_error(
                () -> Erebus.restore_checkpoint_markers(partial, st.markers.core, cfg)
            )
            @test e isa Erebus.CheckpointError
            @test occursin(string(fn), e.msg)
        end
        # The optional-group helpers of the setup and the checkpoint carry the group
        opt = Erebus._extract_optional_marker_arrays(st.markers, cfg, nm, false)
        @test opt.porosity_props === st.markers.groups.porosity
        rebuilt = Erebus._reconstruct_checkpoint_marker_arrays(st.markers.core, opt, nm)
        @test rebuilt.groups.porosity === st.markers.groups.porosity
        pg = st.markers.groups.porosity
        d_raw = Erebus._collect_checkpoint_marker_dict(nothing; porosity_props=pg)
        d_mk = Erebus._collect_checkpoint_marker_dict(st.markers)
        for fn in fieldnames(Erebus.PorosityGroup)
            @test d_raw[fn] === getfield(pg, fn)
            @test d_mk[fn] === getfield(pg, fn)
        end
        # Without the group (porosity_props === nothing) no array is written or rebuilt
        d_none = Erebus._collect_checkpoint_marker_dict(nothing)
        @test !any(haskey(d_none, fn) for fn in fieldnames(Erebus.PorosityGroup))
        opt_none = merge(opt, (; porosity_props=nothing))
        opt_old = Base.structdiff(opt, NamedTuple{(:porosity_props,)})
        for o in (opt_none, opt_old)
            rb = Erebus._reconstruct_checkpoint_marker_arrays(st.markers.core, o, nm)
            @test !haskey(rb.groups, :porosity)
            @test length(rb) == nm
        end
        # The positional save_state writes the group it is given
        m_pos = only(
            m for m in methods(Erebus.save_state) if length(Base.method_argnames(m)) > 20
        )
        dir_pos = mktempdir() * "/"
        lead = Dict(
            :output_path => dir_pos,
            :timestep => 3,
            :dt => 1.0,
            :timesum => 0.0,
            :marknum => nm,
        )
        args = [
            if haskey(lead, a)
                lead[a]
            elseif hasfield(Erebus.GridArrays, a)
                getfield(st.grids, a)
            else
                getfield(st.markers.core, a)
            end for a in Base.method_argnames(m_pos)[2:end]
        ]
        Erebus.save_state(args...; coords=coords2, cfg=cfg, porosity_props=pg)
        data_pos = Erebus.JLD2.load(joinpath(dir_pos, "output_00003.jld2"))
        for fn in fieldnames(Erebus.PorosityGroup)
            @test isapprox(data_pos[string(fn)], getfield(pg, fn); rtol=0, atol=0)
        end
        # A porosity entry that is not a real vector of length marknum is a CheckpointError
        # for each of the 4 fields; the message names the dry porosity array
        for fn in fieldnames(Erebus.PorosityGroup),
            bad in (nothing, fill("x", nm), zeros(3), zeros(nm + 1), zeros(nm, 2))

            partial = copy(data)
            partial[string(fn)] = bad
            e = _hk_error(
                () -> Erebus.restore_checkpoint_markers(partial, st.markers.core, cfg)
            )
            @test e isa Erebus.CheckpointError
            @test occursin("dry porosity array $fn", e.msg)
        end
        # A restart checks [porosity] of the configuration it runs with
        bad_cfg = _hk_with(cfg; reaction=Erebus.ReactionConfig(; active=true))
        e_l = _hk_error(
            () -> Erebus.load_simulation_state(
                path; current_cfg=bad_cfg, force_restart_config=true
            ),
        )
        @test e_l isa ArgumentError
        @test occursin("reaction", e_l.msg)
        path_bad = joinpath(mktempdir(), "ckpt_bad.jld2")
        Erebus.save_state(path_bad, st, coords2, bad_cfg)
        @test_throws ArgumentError Erebus.load_simulation_state(path_bad)
        loaded, _, _ = Erebus.load_simulation_state(path)
        @test isapprox(loaded.markers.phi0, st.markers.phi0; rtol=0, atol=0)
        # push_marker! and resize! give phi0 = 1 (no pressure history), not 0
        mk_p = Erebus.init_marker_arrays(0, cfg, coords)
        Erebus.push_marker!(mk_p; xm=70_000.0, ym=70_000.0, tm=1, tkm=500.0, phi_d=0.6)
        @test isapprox(mk_p.phi0[1], 1.0; rtol=1e-14)
        @test iszero(mk_p.dphi_sinter[1]) && iszero(mk_p.tkm_max[1])
        pr0 = fill(cfg.geometry.psurface, coords.Ny1, coords.Nx1)
        Erebus.advance_marker_dry_porosity!(
            mk_p.groups.porosity,
            mk_p.xm,
            mk_p.ym,
            mk_p.tm,
            mk_p.tkm,
            pr0,
            1.0;
            coords=coords,
            psurface=cfg.geometry.psurface,
            porosity=cfg.porosity,
        )
        # At P = 0 and 500 K the pushed marker keeps phi_d = 0.6 (phi0 = 0 would give 0)
        @test isapprox(mk_p.phi_d[1], 0.6; rtol=1e-14)
        @test isapprox(mk_p.phi0[1], 0.75; rtol=1e-14)
        resize!(mk_p, 4)
        @test isapprox(mk_p.phi0[2:4], ones(3); rtol=1e-14)
        # resize! keeps the phi_0 of the existing marker
        @test isapprox(mk_p.phi0[1], 0.75; rtol=1e-14)
        @test isapprox(mk_p.phi_d[1], 0.6; rtol=1e-14)
        @test all(iszero, mk_p.phi_d[2:4]) && all(iszero, mk_p.dphi_sinter[2:4])
        # init_marker_arrays gives phi0 = 1 too (define_markers! outside init_simulation)
        mk_i = Erebus.init_marker_arrays(5, cfg, coords)
        @test all(p -> isapprox(p, 1.0; rtol=0, atol=0), mk_i.phi0)
        @test all(iszero, mk_i.phi_d) && all(iszero, mk_i.dphi_sinter)
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
        # porosity.model cannot change on restart, also with force_restart_config = true:
        # henke2012 -> none would drop phi_d, none -> henke2012 has no phi_d to load
        for (ckpt, cfg_to) in ((path, _hk_cfg("porosity.model" => :none)), (path_n, cfg))
            e_m = _hk_error(
                () -> Erebus.with_logger(Erebus.NullLogger()) do
                    Erebus.load_simulation_state(
                        ckpt; current_cfg=cfg_to, force_restart_config=true
                    )
                end,
            )
            @test e_m isa Erebus.CheckpointError
            @test occursin("porosity.model cannot change on restart", e_m.msg)
        end
        # Control: the same checkpoint and model restart with a forced other change
        loaded_f, _, _ = Erebus.with_logger(Erebus.NullLogger()) do
            Erebus.load_simulation_state(
                path;
                current_cfg=_hk_cfg("porosity.phi_init" => 0.7),
                force_restart_config=true,
            )
        end
        @test isapprox(loaded_f.markers.phi_d, st.markers.phi_d; rtol=0, atol=0)
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

    @testset "Without pore water the water porosity limits no time step and carries no heat" begin
        # 200 km body at 900 K: sintering shifts pr - pr0 by tens of Pa per step, which
        # the water porosity phimin would turn into a time step limit of seconds
        Rp = 200_000.0
        body = (
            "grid.xsize" => 2.8 * Rp,
            "grid.ysize" => 2.8 * Rp,
            "geometry.rplanet" => Rp,
            "geometry.rcrust" => Rp,
            "geometry.xcenter" => 1.4 * Rp,
            "geometry.ycenter" => 1.4 * Rp,
            "materials.tkm0" => [900.0, 900.0, 170.0],
            "time.n_steps" => 4,
            "time.dt_initial" => 1000.0,
            "time.dt_longest" => 5000.0,
        )
        runs = map((:henke2012, :none)) do model
            dir = mktempdir()
            st = Erebus.with_logger(Erebus.NullLogger()) do
                Erebus.simulation_loop(
                    _hk_cfg(body..., "porosity.model" => model); output_path=dir
                )
            end
            rows = split.(readlines(joinpath(dir, "telemetry.csv")), ",")
            (; st, head=rows[1], last=parse.(Float64, rows[end][1:16]))
        end
        st_h, st_n = runs[1].st, runs[2].st
        yr = Erebus.SECONDS_PER_YEAR
        @test st_h.timestep == 4 && st_n.timestep == 4
        # The step of henke2012 stays within a factor 3 of model none (1358 and 719 yr)
        @test 1 / 3 < st_h.dt / st_n.dt < 3
        @test 100.0 * yr < st_h.dt <= 5000.0 * yr
        # Shear heating stays below 1e-6 W/m^3, far below the 1.7e-3 of a 800 s step
        @test maximum(abs, st_h.grids.HS) < 1.0e-6
        @test maximum(abs, st_h.grids.HS) < 0.1 * maximum(abs, st_h.grids.HR)
        # The pore-water heat is 0 without pore water; with it it is 0.6 HR here
        @test all(iszero, st_h.grids.Q_fluid_grid)
        @test maximum(abs, st_n.grids.Q_fluid_grid) > 0.1 * maximum(abs, st_n.grids.HR)
        # Telemetry: phi_max and phi_mean are the rock phi_d, dt_aphimax_max is 0
        head = runs[1].head
        col(name) = findfirst(==(name), head)
        rock = findall(<(3), st_h.markers.tm)
        phi_d = st_h.markers.phi_d[rock]
        row_h = runs[1].last
        @test isapprox(row_h[col("phi_max")], maximum(phi_d); rtol=1e-14)
        @test isapprox(row_h[col("phi_mean")], sum(phi_d) / length(phi_d); rtol=1e-14)
        # Scale guard: the rock is cold pressed below 0.43, far above phimin = 1e-4
        @test 0.1 < row_h[col("phi_mean")] < row_h[col("phi_max")] < 0.43
        @test iszero(row_h[col("dt_aphimax_max")])
        # Model none keeps the water porosity of the grid and its compaction diagnostic
        row_n = runs[2].last
        PHI = st_n.grids.PHI
        @test isapprox(row_n[col("phi_max")], maximum(PHI); rtol=1e-14)
        @test isapprox(row_n[col("phi_mean")], sum(PHI) / length(PHI); rtol=1e-14)
        @test 1.0e-3 < row_n[col("dt_aphimax_max")] < 0.1
    end

    @testset "Restart from a step-4 checkpoint repeats an 8-step sintering run exactly" begin
        cfg8 = _hk_cfg(
            "materials.tkm0" => [900.0, 900.0, 170.0],
            "time.n_steps" => 8,
            "time.dt_initial" => 1000.0,
            "time.dt_longest" => 5000.0,
            "output.mode" => :snapshots,
            "output.savematstep" => 4,
        )
        cfg4 = _hk_with(cfg8; time=_hk_with(cfg8.time; n_steps=4))
        quiet(f) = Erebus.with_logger(f, Erebus.NullLogger())
        st_full = quiet(() -> Erebus.simulation_loop(cfg8; output_path=mktempdir()))
        dir1 = mktempdir()
        st_half = quiet(() -> Erebus.simulation_loop(cfg4; output_path=dir1))
        st_res = quiet(
            () -> Erebus.simulation_loop(
                cfg8;
                output_path=mktempdir(),
                restart_from=joinpath(dir1, "output_00004.jld2"),
            ),
        )
        @test st_half.timestep == 4
        @test st_full.timestep == 8 && st_res.timestep == 8
        # Guard: the rock sinters between step 4 and step 8, so the restart has work to do
        rock = findall(<(3), st_full.markers.tm)
        gain = st_full.markers.dphi_sinter[rock] .- st_half.markers.dphi_sinter[rock]
        @test 1.0e-3 < maximum(gain) < 0.75
        # The checkpoint keeps the full state and the RNG, so the 2 runs agree bit for bit
        for fn in (:phi_d, :phi0, :tkm_max, :dphi_sinter, :tkm)
            @test isapprox(
                getproperty(st_res.markers, fn),
                getproperty(st_full.markers, fn);
                rtol=0,
                atol=0,
            )
        end
        for fn in (:tk2, :pr)
            @test isapprox(
                getfield(st_res.grids, fn), getfield(st_full.grids, fn); rtol=0, atol=0
            )
        end
        @test isapprox(st_res.timesum, st_full.timesum; rtol=0, atol=0)
    end
end

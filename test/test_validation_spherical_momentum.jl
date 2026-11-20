using Logging: NullLogger, with_logger
using SHA
using Test
using TOML
using Erebus

isdefined(@__MODULE__, :SphericalMomentumValidation) ||
    include(joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_momentum.jl"))

include(joinpath(@__DIR__, "..", "tools", "check_validation.jl"))
isdefined(@__MODULE__, :SphericalMomentumExport) ||
    include(joinpath(@__DIR__, "..", "tools", "export_validation_spherical_momentum.jl"))

const SMV = SphericalMomentumValidation
const SMVE = SphericalMomentumExport
const QB = Rational{BigInt}

"""
    MomentumInterruptingBench

Benchmark stand-in whose every property access throws `InterruptException`, as a Ctrl-C
during the benchmark does.
"""
struct MomentumInterruptingBench end
Base.getproperty(::MomentumInterruptingBench, ::Symbol) = throw(InterruptException())

"""
    MomentumCrashingBench

Benchmark stand-in whose every property access throws an `ErrorException`, as a solver
failure during the benchmark does.
"""
struct MomentumCrashingBench end
function Base.getproperty(::MomentumCrashingBench, ::Symbol)
    return error("solver failure inside the benchmark")
end

"""
    write_momentum_toml(path, dict)

Write `dict` as TOML to `path` and return `path`.

Parameters
----------
path : String
    Output file path.
dict : AbstractDict
    Table to write.

Returns
-------
String
    `path`.
"""
function write_momentum_toml(path::String, dict::AbstractDict)
    open(path, "w") do io
        return TOML.print(io, dict)
    end
    return path
end

"""
    capture_momentum_stdio(f)

Call `f()` with standard output and standard error redirected to files and return the
result with the captured text.

Parameters
----------
f : Function
    Function of no arguments.

Returns
-------
NamedTuple
    `value`: the return value of `f`; `out`, `err`: the captured standard output and
    standard error.
"""
function capture_momentum_stdio(f)
    return mktempdir() do dir
        out_path, err_path = joinpath(dir, "out.txt"), joinpath(dir, "err.txt")
        value = open(out_path, "w") do out
            open(err_path, "w") do err
                return redirect_stdio(f; stdout=out, stderr=err)
            end
        end
        return (; value, out=read(out_path, String), err=read(err_path, String))
    end
end

"""
    rpoly_mul(a, b)

Return the product of 2 polynomials in one variable with exact rational coefficients; entry
`k + 1` holds the coefficient of the power `k`.

Parameters
----------
a, b : Vector{Rational{BigInt}}
    Coefficients.

Returns
-------
Vector{Rational{BigInt}}
    Coefficients of `a b`.
"""
function rpoly_mul(a::Vector{QB}, b::Vector{QB})
    c = zeros(QB, length(a) + length(b) - 1)
    for i in eachindex(a), j in eachindex(b)
        c[i + j - 1] += a[i] * b[j]
    end
    return c
end

"""
    rpoly_pow(a, n)

Return the `n`-th power of the polynomial `a` (see `rpoly_mul`).

Parameters
----------
a : Vector{Rational{BigInt}}
    Coefficients.
n : Int
    Power, at least 0.

Returns
-------
Vector{Rational{BigInt}}
    Coefficients of `a^n`.
"""
function rpoly_pow(a::Vector{QB}, n::Int)
    c = QB[1]
    for _ in 1:n
        c = rpoly_mul(c, a)
    end
    return c
end

"""
    rpoly_add(a, b)

Return the sum of 2 polynomials (see `rpoly_mul`).

Parameters
----------
a, b : Vector{Rational{BigInt}}
    Coefficients.

Returns
-------
Vector{Rational{BigInt}}
    Coefficients of `a + b`.
"""
function rpoly_add(a::Vector{QB}, b::Vector{QB})
    n = max(length(a), length(b))
    return [(k <= length(a) ? a[k] : QB(0)) + (k <= length(b) ? b[k] : QB(0)) for k in 1:n]
end

"""
    rpoly_der(a)

Return the derivative of the polynomial `a` (see `rpoly_mul`).

Parameters
----------
a : Vector{Rational{BigInt}}
    Coefficients.

Returns
-------
Vector{Rational{BigInt}}
    Coefficients of `da/dx`.
"""
function rpoly_der(a::Vector{QB})
    return length(a) == 1 ? QB[0] : [k * a[k + 1] for k in 1:(length(a) - 1)]
end

"""
    rpoly_shift(a, m)

Return the polynomial `a` times `x^m` for `m >= 0` (see `rpoly_mul`).

Parameters
----------
a : Vector{Rational{BigInt}}
    Coefficients.
m : Int
    Power of the factor, at least 0.

Returns
-------
Vector{Rational{BigInt}}
    Coefficients of `x^m a`.
"""
rpoly_shift(a::Vector{QB}, m::Int) = vcat(zeros(QB, m), a)

"""
    rpoly_int01(a)

Return the exact integral of the polynomial `a` over `[0, 1]`.

Parameters
----------
a : Vector{Rational{BigInt}}
    Coefficients.

Returns
-------
Rational{BigInt}
    `int_0^1 a(x) dx`.
"""
rpoly_int01(a::Vector{QB}) = sum(a[k] // (k) for k in eachindex(a))

"""
    rpoly_iszero(a)

Return whether every coefficient of the polynomial `a` is exactly 0.

Parameters
----------
a : Vector{Rational{BigInt}}
    Coefficients.

Returns
-------
Bool
"""
rpoly_iszero(a::Vector{QB}) = all(iszero, a)

"""
    rho_poly_from_s(c)

Return the coefficients in `rho` of the polynomial `sum_k c[k + 1] rho^(2 k)` given by its
coefficients `c` in `s = rho^2`.

Parameters
----------
c : Vector{Rational{BigInt}}
    Coefficients in `s`.

Returns
-------
Vector{Rational{BigInt}}
    Coefficients in `rho`.
"""
function rho_poly_from_s(c::Vector{QB})
    out = zeros(QB, 2 * length(c) - 1)
    for k in eachindex(c)
        out[2k - 1] = c[k]
    end
    return out
end

"""
    radial_residual(omega, force, k)

Return the coefficients of `(rho^(k + 1) Omega')' + rho^k f` in `rho`. The on form
(`k = 3`) is the stress divergence `(1 / r^3) (r^3 tau)'` and the off form (`k = 2`)
`(1 / r^2) (r^2 tau)'` of `tau = eta r Omega'`, both times `rho^k`: the angular component
of the momentum equation, Eq. (3) of Hernlund and Tackley (2008), for an azimuthal flow
with `d = 2` (sphere) and `d = 1` (cylinder, the 2-D form) and the azimuthal body force
`f` added as in their Eq. (14).

Parameters
----------
omega : Vector{Rational{BigInt}}
    `Omega(rho)` in units of `F0 R / eta`, coefficients in `rho`.
force : Vector{Rational{BigInt}}
    Azimuthal force `f(rho) / F0`, coefficients in `rho`.
k : Int
    3 for `spherical_flow = true`, 2 for `false`.

Returns
-------
Vector{Rational{BigInt}}
    Residual coefficients; all 0 for an exact solution.
"""
function radial_residual(omega::Vector{QB}, force::Vector{QB}, k::Int)
    flux = rpoly_shift(rpoly_der(omega), k + 1)
    return rpoly_add(rpoly_der(flux), rpoly_shift(force, k))
end

const ONE_MINUS_RHO2 = QB[1, 0, -1]

@testset "Spherical Momentum Validation" begin
    @testset "case A: closed forms solve the radial momentum equation with and without the option" begin
        # The closed forms written out: Omega in units of F0 R / eta, f in units of F0.
        w2 = ONE_MINUS_RHO2
        om_on = rpoly_mul(rpoly_pow(w2, 4), QB[1, 0, -6]) .// 100
        om_off = rpoly_mul(rpoly_pow(w2, 4), QB[3, 0, -13]) .// 200
        f_a = rpoly_shift(rpoly_mul(rpoly_pow(w2, 2), QB[1, 0, -32 // 5, 0, 39 // 5]), 1)
        @test rpoly_iszero(radial_residual(om_on, f_a, 3))
        @test rpoly_iszero(radial_residual(om_off, f_a, 2))
        # The same force does not solve the other form: the metric term is O(1).
        @test !rpoly_iszero(radial_residual(om_off, f_a, 3))
        @test !rpoly_iszero(radial_residual(om_on, f_a, 2))
        # C4: the moments int rho^2 f and int rho^3 f vanish, so Omega' = 0 at rho = 1 in
        # both forms and the flow is 0 outside R.
        @test iszero(rpoly_int01(rpoly_shift(f_a, 2)))
        @test iszero(rpoly_int01(rpoly_shift(f_a, 3)))
        @test !iszero(rpoly_int01(rpoly_shift(f_a, 1)))
        # The module holds these closed forms, in s = rho^2.
        @test rho_poly_from_s(SMV.omega_coefficients(:A; spherical_flow=true)) == om_on
        @test rho_poly_from_s(SMV.omega_coefficients(:A; spherical_flow=false)) == om_off
        @test rpoly_shift(rho_poly_from_s(SMV.force_coefficients(:A)), 1) == f_a
        # Centre rotation: 1/100 with the option and 3/200 without; ratio 2/3.
        @test SMV.omega_coefficients(:A; spherical_flow=true)[1] == 1 // 100
        @test SMV.omega_coefficients(:A; spherical_flow=false)[1] == 3 // 200
        # C5: Omega and Omega' vanish at rho = 1 in both forms.
        for om in (om_on, om_off)
            @test iszero(sum(om))
            @test iszero(sum(rpoly_der(om)))
        end
    end

    @testset "case A: Float64 evaluation against finite differences of the closed form" begin
        h = 1.0e-3
        d4(g, x) = (g(x - 2h) - 8g(x - h) + 8g(x + h) - g(x + 2h)) / (12h)
        for sf in (true, false)
            k = sf ? 3 : 2
            om(rho) = SMV.omega_hat(:A, rho; spherical_flow=sf)
            for rho in (0.1, 0.33, 0.5, 0.71, 0.9)
                # (1 / rho^k) d/drho (rho^(k + 1) Omega') + f / F0 = 0, 4th-order steps.
                flux(x) = x^(k + 1) * d4(om, x)
                div = d4(flux, rho) / rho^k
                f = SMV.force_hat(:A, rho)
                @test isapprox(div + f, 0.0; atol=1.0e-6 * abs(f) + 1.0e-9)
            end
            # Outside R the flow is 0, and on the surface it is 0 to round-off.
            @test iszero(SMV.omega_hat(:A, 1.2; spherical_flow=sf))
            @test isapprox(
                SMV.omega_hat(:A, 1.0 - 1.0e-12; spherical_flow=sf), 0.0; atol=1.0e-15
            )
        end
        @test isapprox(SMV.omega_hat(:A, 0.0; spherical_flow=true), 0.01; rtol=1.0e-15)
        @test isapprox(SMV.omega_hat(:A, 0.0; spherical_flow=false), 0.015; rtol=1.0e-15)
        # The velocity is u_theta e_theta, u_theta = rho Omega: (vx, vy) = Omega (-Y, X).
        X, Y = 0.3, -0.4
        om = SMV.omega_hat(:A, 0.5; spherical_flow=true)
        vx, vy = SMV.exact_velocity(:A, X, Y; spherical_flow=true)
        @test isapprox(vx, -Y * om; rtol=1.0e-14)
        @test isapprox(vy, X * om; rtol=1.0e-14)
        fx, fy = SMV.exact_force(:A, X, Y)
        @test isapprox(-fx * Y + fy * X, 0.5 * SMV.force_hat(:A, 0.5); rtol=1.0e-14)
        # The force is azimuthal: no radial part.
        @test isapprox(fx * X + fy * Y, 0.0; atol=1.0e-16)
    end

    @testset "case A0: Omega(0) = 0, its closed form and its force" begin
        w2 = ONE_MINUS_RHO2
        om_a0 = rpoly_shift(rpoly_pow(w2, 4), 2) .// -10
        f_a0 = rpoly_shift(rpoly_mul(rpoly_pow(w2, 2), QB[1, 0, -46 // 5, 0, 13]), 1)
        @test rpoly_iszero(radial_residual(om_a0, f_a0, 3))
        @test !rpoly_iszero(radial_residual(om_a0, f_a0, 2))
        @test iszero(rpoly_int01(rpoly_shift(f_a0, 3)))
        @test iszero(om_a0[1])
        @test iszero(sum(om_a0)) && iszero(sum(rpoly_der(om_a0)))
        @test rho_poly_from_s(SMV.omega_coefficients(:A0; spherical_flow=true)) == om_a0
        @test rpoly_shift(rho_poly_from_s(SMV.force_coefficients(:A0)), 1) == f_a0
        # The centre rotates as rho^2: 0 at the centre and negative inside.
        @test iszero(SMV.omega_hat(:A0, 0.0; spherical_flow=true))
        @test -0.1 < SMV.omega_hat(:A0, 0.4; spherical_flow=true) < 0.0
        # A0 has no closed form without the option.
        @test_throws ArgumentError SMV.omega_hat(:A0, 0.5; spherical_flow=false)
        @test_throws ArgumentError SMV.exact_velocity(:A0, 0.1, 0.2; spherical_flow=false)
    end

    @testset "case B: spherical divergence 0, polynomial force finite at the centre" begin
        # rho^2 div_sph u = rho^2 (d_x ux + d_y uy) + X ux + Y uy, exactly 0.
        @test SMV.spherical_divergence_is_zero()
        @test SMV.spherical_divergence_is_zero(3)
        @test_throws ArgumentError SMV.spherical_divergence_is_zero(1)
        # The force is a polynomial plus q / rho^2 with q of degree >= 3: rho^2 does not
        # divide q, but q / rho^2 <= C rho, so the force is finite and 0 at the centre.
        @test !isempty(SMV.CASE_B_EXACT.fx_rem) && !isempty(SMV.CASE_B_EXACT.fy_rem)
        @test minimum(sum(k) for k in keys(SMV.CASE_B_EXACT.fx_rem)) == 3
        fx0, fy0 = SMV.exact_force(:B, 0.0, 0.0)
        @test iszero(fx0) && iszero(fy0)
        for eps in (1.0e-3, 1.0e-6, 1.0e-9), t in (0.3, 1.1, 2.5)
            f = SMV.exact_force(:B, eps * cos(t), eps * sin(t))
            @test all(isfinite, f)
            @test hypot(f...) <= 100 * eps
        end
        # C5: u = 0 for r >= R; just inside the circle to round-off, where the polynomial
        # is evaluated, and on rational points of the circle exactly.
        r_in = 1.0 - 1.0e-12
        for t in range(0.0, 2pi; length=13)
            @assert (r_in * cos(t))^2 + (r_in * sin(t))^2 < 1
            vx, vy = SMV.exact_velocity(:B, r_in * cos(t), r_in * sin(t))
            @test hypot(vx, vy) <= 1.0e-13
        end
        @test SMV.velocity_vanishes_on_circle(:B)
        @test SMV.velocity_vanishes_on_circle(:A)
        @test SMV.velocity_vanishes_on_circle(:A0)
        @test isapprox(
            SMV.omega_hat(:A0, 1.0 - 1.0e-12; spherical_flow=true), 0.0; atol=1.0e-15
        )
        @test all(iszero, SMV.exact_velocity(:B, 0.9, 0.6))
        @test all(iszero, SMV.exact_force(:B, 0.9, 0.6))
        # u_r != 0 and tau_33 = 2 u_r / r: O(1) inside.
        X, Y = 0.4, 0.1
        vx, vy = SMV.exact_velocity(:B, X, Y)
        ur = (vx * X + vy * Y) / hypot(X, Y)
        @test isapprox(SMV.exact_tau33(:B, X, Y), 2 * ur / hypot(X, Y); rtol=1.0e-12)
        @test 0.01 < abs(ur) < 1.0
        @test iszero(SMV.exact_tau33(:A, X, Y))
        @test_throws ArgumentError SMV.exact_velocity(:B, 0.1, 0.2; spherical_flow=false)
    end

    @testset "case B: force from finite differences of the exact velocity (second evaluation)" begin
        # tau_ij = d_i u_j + d_j u_i, tau_33 = 2 (X ux + Y uy) / rho^2, and the force
        # f_i = -(d_j tau_ij + (tau_ij X_j - tau_33 X_i) / rho^2), 4th-order differences.
        u(X, Y) = SMV.exact_velocity(:B, X, Y)
        d4(g, x, h) = (g(x - 2h) - 8g(x - h) + 8g(x + h) - g(x + 2h)) / (12h)
        h = 1.0e-3
        function tau(X, Y)
            dxu = (d4(x -> u(x, Y)[1], X, h), d4(x -> u(x, Y)[2], X, h))
            dyu = (d4(y -> u(X, y)[1], Y, h), d4(y -> u(X, y)[2], Y, h))
            vx, vy = u(X, Y)
            return (
                2 * dxu[1], 2 * dyu[2], dyu[1] + dxu[2], 2 * (X * vx + Y * vy) / (X^2 + Y^2)
            )
        end
        fmax = 0.0
        for (X, Y) in ((0.3, 0.2), (-0.5, 0.35), (0.12, -0.7), (-0.25, -0.25), (0.6, 0.05))
            t = tau(X, Y)
            r2 = X^2 + Y^2
            dx_txx = d4(x -> tau(x, Y)[1], X, h)
            dy_txy = d4(y -> tau(X, y)[3], Y, h)
            dx_txy = d4(x -> tau(x, Y)[3], X, h)
            dy_tyy = d4(y -> tau(X, y)[2], Y, h)
            fx = -(dx_txx + dy_txy + (t[1] * X + t[3] * Y - t[4] * X) / r2)
            fy = -(dx_txy + dy_tyy + (t[3] * X + t[2] * Y - t[4] * Y) / r2)
            fx_m, fy_m = SMV.exact_force(:B, X, Y)
            scale = max(abs(fx_m), abs(fy_m))
            fmax = max(fmax, scale)
            @test isapprox(fx, fx_m; atol=1.0e-6 * scale)
            @test isapprox(fy, fy_m; atol=1.0e-6 * scale)
        end
        # Scale guard: the force is O(10) in units of eta U0 / R^2.
        @test 1.0 < fmax < 1.0e3
    end

    @testset "case B force from the spherical annulus equations of Hernlund and Tackley (2008)" begin
        # Eqs. (2), (3) with d = 2 and Eqs. (4) to (7) for div u = 0, eta = 1:
        # tau_rr = 2 d_r v_r, tau_thth = 2 v_r / r, tau_phiphi = 2 (d_phi v_phi + v_r) / r,
        # tau_rphi = d_phi v_r / r + r d_r(v_phi / r).
        function polar_u(r, p)
            vx, vy = SMV.exact_velocity(:B, r * cos(p), r * sin(p))
            return (vx * cos(p) + vy * sin(p), -vx * sin(p) + vy * cos(p))
        end
        h = 1.0e-3
        d4(g, x) = (g(x - 2h) - 8g(x - h) + 8g(x + h) - g(x + 2h)) / (12h)
        function stresses(r, p)
            vr, vp = polar_u(r, p)
            dr_vr = d4(x -> polar_u(x, p)[1], r)
            dp_vr = d4(y -> polar_u(r, y)[1], p)
            dp_vp = d4(y -> polar_u(r, y)[2], p)
            dr_om = d4(x -> polar_u(x, p)[2] / x, r)
            return (
                rr=2 * dr_vr,
                pp=2 * (dp_vp / r + vr / r),
                tt=2 * vr / r,
                rp=dp_vr / r + r * dr_om,
            )
        end
        for (r, p) in ((0.4, 0.3), (0.62, 2.1), (0.25, -1.2), (0.8, 4.0))
            t = stresses(r, p)
            d_r =
                d4(x -> x^2 * stresses(x, p).rr, r) / r^2 +
                d4(y -> stresses(r, y).rp, p) / r - (t.pp + t.tt) / r
            d_p =
                d4(x -> x^2 * stresses(x, p).rp, r) / r^2 +
                d4(y -> stresses(r, y).pp, p) / r +
                t.rp / r
            # The force balances the stress divergence: f = -D, in Cartesian components.
            fx = -(d_r * cos(p) - d_p * sin(p))
            fy = -(d_r * sin(p) + d_p * cos(p))
            fx_m, fy_m = SMV.exact_force(:B, r * cos(p), r * sin(p))
            scale = max(abs(fx_m), abs(fy_m))
            @test isapprox(fx, fx_m; atol=1.0e-6 * scale)
            @test isapprox(fy, fy_m; atol=1.0e-6 * scale)
        end
        # The cylinder form (d = 1: no tau_thth, r instead of r^2) gives another force.
        r, p = 0.4, 0.3
        t = stresses(r, p)
        d_r1 =
            d4(x -> x * stresses(x, p).rr, r) / r + d4(y -> stresses(r, y).rp, p) / r -
            t.pp / r
        fx_m, fy_m = SMV.exact_force(:B, r * cos(p), r * sin(p))
        f_r = fx_m * cos(p) + fy_m * sin(p)
        @test !isapprox(-d_r1, f_r; rtol=0.05)
    end

    # Case C in units of U0, U0 / R and eta U0 / R: u_r = -10 rho (1 - rho^2)^4 from the
    # potential (1 - rho^2)^5, P_f - p_s = lam (1 - rho^2)^5 and Z / eta = zr.
    lam_c = SMV.MU_F * SMV.R_PLANET^2 / (SMV.K_PERM * SMV.ETA_ROCK)
    zr_c = (1 - SMV.PHI) / SMV.ETA_RATIO_C
    ur_c(rho) = -10 * rho * (1 - rho^2)^4
    points_c = ((0.3, 0.2), (-0.5, 0.35), (0.12, -0.7), (-0.25, -0.25), (0.6, 0.05))

    @testset "case C: divergence, solid and fluid continuity and Darcy of the compaction flow" begin
        h = 1.0e-3
        d4(g, x) = (g(x - 2h) - 8g(x - h) + 8g(x + h) - g(x + 2h)) / (12h)
        u(X, Y) = SMV.exact_velocity(:C, X, Y)
        p(X, Y) = SMV.exact_pressures(:C, X, Y)
        @test isapprox(lam_c, 0.25; rtol=1.0e-12)
        @test isapprox(zr_c, 0.8; rtol=1.0e-12)
        for (X, Y) in points_c
            rho = hypot(X, Y)
            vx, vy = u(X, Y)
            @test isapprox(vx, ur_c(rho) * X / rho; rtol=1.0e-12)
            @test isapprox(vy, ur_c(rho) * Y / rho; rtol=1.0e-12)
            # div3 = d_x ux + d_y uy + (X ux + Y uy) / rho^2 by finite differences.
            div_fd =
                d4(x -> u(x, Y)[1], X) + d4(y -> u(X, y)[2], Y) + (X * vx + Y * vy) / rho^2
            div = SMV.exact_div3(:C, X, Y)
            @test isapprox(div, div_fd; rtol=1.0e-8)
            pt, pf = p(X, Y)
            # Solid continuity: div3 u + (P_t - P_f) / Z = 0.
            @test isapprox(div + (pt - pf) / zr_c, 0.0; atol=1.0e-12 * abs(div))
            # Darcy with rho_f = 0: q = -(k / mu) grad P_f = -grad(P_f) / lam = -u.
            qx(x, y) = -d4(s -> p(s, y)[2], x) / lam_c
            qy(x, y) = -d4(s -> p(x, s)[2], y) / lam_c
            @test isapprox(qx(X, Y), -vx; rtol=1.0e-8)
            @test isapprox(qy(X, Y), -vy; rtol=1.0e-8)
            # Fluid continuity: div3 q - (P_t - P_f) / Z = 0, from second differences.
            div_q =
                d4(x -> qx(x, Y), X) +
                d4(y -> qy(X, y), Y) +
                (X * qx(X, Y) + Y * qy(X, Y)) / rho^2
            @test isapprox(div_q - (pt - pf) / zr_c, 0.0; atol=1.0e-6 * abs(div))
        end
        # Centre: isotropic compaction with div3 = -30 and e_33 = -10, tau_33 = 0, P_t - p_s
        # = lam + 30 zr = 24.25 and the force 0.
        @test isapprox(SMV.exact_div3(:C, 0.0, 0.0), -30.0; rtol=1.0e-14)
        @test isapprox(SMV.exact_tau33(:C, 0.0, 0.0), 0.0; atol=1.0e-14)
        @test isapprox(collect(SMV.exact_force(:C, 0.0, 0.0)), [0.0, 0.0]; atol=1.0e-14)
        @test isapprox(
            collect(SMV.exact_pressures(:C, 0.0, 0.0)), [24.25, 0.25]; rtol=1.0e-14
        )
        @test isapprox(SMV.exact_tau33(:C, 1.0e-9, 0.0), 0.0; atol=1.0e-14)
        # Outside R every field is 0, and just inside it vanishes to round-off.
        for f in (
            SMV.exact_velocity(:C, 0.9, 0.6),
            SMV.exact_force(:C, 0.9, 0.6),
            SMV.exact_pressures(:C, 0.9, 0.6),
        )
            @test all(iszero, f)
        end
        @test iszero(SMV.exact_tau33(:C, 0.9, 0.6)) && iszero(SMV.exact_div3(:C, 0.9, 0.6))
        r_in = 1.0 - 1.0e-12
        for t in range(0.0, 2pi; length=13)
            X, Y = r_in * cos(t), r_in * sin(t)
            @test hypot(SMV.exact_velocity(:C, X, Y)...) <= 1.0e-12
            @test maximum(abs, SMV.exact_pressures(:C, X, Y)) <= 1.0e-10
            @test hypot(SMV.exact_force(:C, X, Y)...) <= 1.0e-9
        end
        @test SMV.velocity_vanishes_on_circle(:C)
        # P_t > 0 everywhere with p_s = 1e7 Pa; with p_s = 100 Pa the outer shell, where the
        # rock dilates, would have P_t < 0 (tension).
        sc = SMV.stress_scale(:C)
        pt_min = minimum(
            SMV.exact_pressures(:C, rho, 0.0)[1] for rho in range(0.0, 1.0; length=2001)
        )
        @test -5.0 < pt_min < 0.0
        @test 0.0 < SMV.P_AMB_C + sc * pt_min < SMV.P_AMB_C
        @test SMV.P_AMB + sc * pt_min < -1.0e6
        @test isapprox(SMV.P_AMB_C, 1.0e7; rtol=1.0e-15)
        # P_f - p_s >= 0: the fluid pressure peaks at the centre.
        @test all(
            SMV.exact_pressures(:C, rho, 0.0)[2] >= 0.0 for
            rho in range(0.0, 1.0; length=101)
        )
    end

    @testset "case C force against the radial momentum balance of a 3-D sphere" begin
        # sigma_rr = -P_t + 2 (u_r' - div3 / 3), sigma_thth = sigma_phph = -P_t +
        # 2 (u_r / r - div3 / 3), div3 = u_r' + 2 u_r / r and P_t = P_f - zr div3, eta = 1;
        # f_r = -(d_r sigma_rr + 2 (sigma_rr - sigma_thth) / r), BigFloat differences in r.
        fr_max, worst, no_trace_gap = 0.0, 0.0, 0.0
        setprecision(BigFloat, 256) do
            h = BigFloat(10)^-25
            urb(r) = -10 * r * (1 - r^2)^4
            function sig(r, trace)
                du = (urb(r + h) - urb(r - h)) / (2h)
                div = du + 2 * urb(r) / r
                pt = BigFloat(lam_c) * (1 - r^2)^5 - BigFloat(zr_c) * div
                t = trace ? div / 3 : zero(div)
                return (-pt + 2 * (du - t), -pt + 2 * (urb(r) / r - t))
            end
            function f_r(r, trace)
                srr, stt = sig(r, trace)
                dsrr = (sig(r + h, trace)[1] - sig(r - h, trace)[1]) / (2h)
                return Float64(-(dsrr + 2 * (srr - stt) / r))
            end
            for r in 0.05:0.05:0.95
                fr = f_r(BigFloat(r), true)
                fx, fy = SMV.exact_force(:C, r, 0.0)
                gx, gy = SMV.exact_force(:C, r / sqrt(2), r / sqrt(2))
                worst = max(
                    worst,
                    abs(fx - fr),
                    abs(fy),
                    abs((gx + gy) / sqrt(2) - fr),
                    abs(gx - gy),
                )
                fr_max = max(fr_max, abs(fr))
                no_trace_gap = max(no_trace_gap, abs(f_r(BigFloat(r), false) - fx))
            end
        end
        @test worst <= 1.0e-6 * fr_max
        # Scale guard: max |f_r| is O(100) in units of eta U0 / R^2.
        @test 10.0 < fr_max < 1.0e3
        # Without the trace, tau_ij = 2 eps_ij, the balance needs a force that differs by
        # more than 10 % of the largest force.
        @test no_trace_gap > 0.1 * fr_max
    end

    @testset "case C second invariants and dissipation against the 3-D strain rate" begin
        h = 1.0e-3
        d4(g, x) = (g(x - 2h) - 8g(x - h) + 8g(x + h) - g(x + 2h)) / (12h)
        for (X, Y) in points_c
            rho = hypot(X, Y)
            du = d4(ur_c, rho)
            div = du + 2 * ur_c(rho) / rho
            # Radial flow: e'_rr = u_r' - div3 / 3, e'_thth = e'_phph = u_r / r - div3 / 3.
            err, ett = du - div / 3, ur_c(rho) / rho - div / 3
            eii = sqrt((err^2 + 2 * ett^2) / 2)
            @test isapprox(SMV.exact_eii(:C, X, Y), eii; rtol=1.0e-8)
            @test isapprox(
                SMV.exact_eii(:C, X, Y), 80 / sqrt(3) * rho^2 * (1 - rho^2)^3; rtol=1.0e-12
            )
            # S_II = 2 eta E_II in the viscous limit, eta = 1.
            @test isapprox(SMV.exact_sii(:C, X, Y), 2 * eii; rtol=1.0e-8)
            # Dissipation: tau : eps = 4 E_II^2, compaction (P_t - P_f)^2 / Z, Darcy
            # (mu / k) |q|^2 = lam u_r^2.
            pt, pf = SMV.exact_pressures(:C, X, Y)
            hs = 4 * eii^2 + (pt - pf)^2 / zr_c + lam_c * ur_c(rho)^2
            @test isapprox(SMV.exact_hs(:C, X, Y), hs; rtol=1.0e-8)
        end
        # The centre compacts isotropically: E_II = S_II = 0 there. Without the trace the
        # out-of-plane part 0.75 e_33^2 would give E_II(0) = sqrt(0.75) 10.
        @test isapprox(SMV.exact_eii(:C, 0.0, 0.0), 0.0; atol=1.0e-12)
        @test isapprox(SMV.exact_sii(:C, 0.0, 0.0), 0.0; atol=1.0e-12)
        @test isapprox(SMV.exact_hs(:C, 0.0, 0.0), zr_c * 900; rtol=1.0e-12)
        @test 1.0 < SMV.exact_eii(:C, 0.5, 0.0) < 10.0
        @test all(iszero, (SMV.exact_eii(:C, 0.9, 0.6), SMV.exact_hs(:C, 0.9, 0.6)))
    end

    @testset "scales, guards and the case list" begin
        @test isapprox(
            SMV.velocity_scale(:A), SMV.F0 * SMV.R_PLANET^2 / SMV.ETA_ROCK; rtol=1.0e-15
        )
        @test isapprox(SMV.velocity_scale(:B), SMV.U0; rtol=1.0e-15)
        @test isapprox(SMV.force_scale(:A0), SMV.F0; rtol=1.0e-15)
        @test isapprox(
            SMV.force_scale(:B), SMV.ETA_ROCK * SMV.U0 / SMV.R_PLANET^2; rtol=1.0e-15
        )
        @test isapprox(SMV.stress_scale(:A), SMV.F0 * SMV.R_PLANET; rtol=1.0e-15)
        @test isapprox(
            SMV.stress_scale(:B), SMV.ETA_ROCK * SMV.U0 / SMV.R_PLANET; rtol=1.0e-15
        )
        # Velocity scale 2.5e-9 m/s for A: the flow stays far below the time step limit.
        @test 1.0e-9 < SMV.velocity_scale(:A) < 1.0e-8
        for bad in (:D, :a, :B0)
            @test_throws ArgumentError SMV.velocity_scale(bad)
            @test_throws ArgumentError SMV.exact_velocity(bad, 0.1, 0.1)
            @test_throws ArgumentError SMV.omega_coefficients(bad; spherical_flow=true)
        end
        @test_throws ArgumentError SMV.force_coefficients(:B)
        @test_throws ArgumentError SMV.omega_hat(:B, 0.5; spherical_flow=true)
        # Case C: a radial compaction flow, with the option on only.
        @test SMV.CASES == (:A, :A0, :B, :C)
        @test isapprox(SMV.velocity_scale(:C), SMV.U0; rtol=1.0e-15)
        @test isapprox(SMV.stress_scale(:C), SMV.stress_scale(:B); rtol=1.0e-15)
        @test_throws ArgumentError SMV.force_coefficients(:C)
        @test_throws ArgumentError SMV.omega_coefficients(:C; spherical_flow=true)
        @test_throws ArgumentError SMV.exact_velocity(:C, 0.1, 0.2; spherical_flow=false)
        @test_throws ArgumentError SMV.make_config(33; spherical_flow=false, case=:C)
        for f in (
            SMV.exact_pressures, SMV.exact_div3, SMV.exact_eii, SMV.exact_sii, SMV.exact_hs
        )
            @test_throws ArgumentError f(:B, 0.1, 0.2)
            @test_throws ArgumentError f(:D, 0.1, 0.2)
            @test_throws DomainError f(:C, NaN, 0.2)
            @test_throws DomainError f(:C, 0.1, Inf)
        end
        @test_throws DomainError SMV.exact_force(:C, Inf, 0.0)
        # The force vanishes on r = R only for a power of at least 4 in the potential.
        @test_throws ArgumentError SMV._case_c_polynomials(3, QB(1, 4), QB(4, 5))
        grid33 = SMV.grid_of(
            Erebus.GridCoordinates(33, 33; xsize=SMV.L_BOX, ysize=SMV.L_BOX)
        )
        @test_throws ArgumentError SMV.exact_pnode_fields(grid33, :B)
        @test_throws ArgumentError SMV.pnode_errors(
            ones(grid33.Ny1, grid33.Nx1), zeros(grid33.Ny1, grid33.Nx1), grid33
        )
        @test_throws DomainError SMV.anchor_deviation(
            ones(3, 3), ones(3, 3), (; Ny1=3, Nx1=3), 0.0
        )
        @test isapprox(SMV.ambient_pressure(:C), 1.0e7; rtol=1.0e-15)
        @test isapprox(SMV.ambient_pressure(:B), SMV.P_AMB; rtol=1.0e-15)
        @test isapprox(SMV.eta_zeta_ratio(:C), 1.0; rtol=1.0e-15)
        @test isapprox(SMV.eta_zeta_ratio(:A0), 1.0e-6; rtol=1.0e-15)
        @test_throws ArgumentError SMV.ambient_pressure(:D)
        @test_throws ArgumentError SMV.eta_zeta_ratio(:D)
        for bad in (0.0, -1.0, NaN, Inf)
            @test_throws DomainError SMV.velocity_scale(:A; r_planet=bad)
            @test_throws DomainError SMV.velocity_scale(:A; eta=bad)
            @test_throws DomainError SMV.force_scale(:B; eta=bad)
        end
        for bad in (NaN, Inf, -Inf)
            @test_throws DomainError SMV.velocity_scale(:A; f0=bad)
            @test_throws DomainError SMV.velocity_scale(:B; u0=bad)
            @test_throws DomainError SMV.exact_velocity(:A, bad, 0.0)
            @test_throws DomainError SMV.exact_force(:B, 0.0, bad)
            @test_throws DomainError SMV.omega_hat(:A, bad; spherical_flow=true)
        end
        @test_throws DomainError SMV.omega_hat(:A, -0.1; spherical_flow=true)
        @test_throws DomainError SMV.force_hat(:A, -0.1)
    end

    @testset "exact on/off L2-norm ratio of the case A profile" begin
        # sqrt(int rho u_on^2 / int rho u_off^2) over [0, 1], u = rho Omega.
        w2 = ONE_MINUS_RHO2
        u_on = rpoly_shift(rpoly_mul(rpoly_pow(w2, 4), QB[1, 0, -6]) .// 100, 1)
        u_off = rpoly_shift(rpoly_mul(rpoly_pow(w2, 4), QB[3, 0, -13]) .// 200, 1)
        num = rpoly_int01(rpoly_shift(rpoly_mul(u_on, u_on), 1))
        den = rpoly_int01(rpoly_shift(rpoly_mul(u_off, u_off), 1))
        @test SMV.profile_ratio_squared() == num / den
        ratio = SMV.profile_norm_ratio()
        @test isapprox(ratio^2, Float64(num / den); rtol=1.0e-15)
        # Midpoint quadrature of the Float64 closed forms gives the same value.
        n = 20_000
        xs = ((1:n) .- 0.5) ./ n
        q(sf) = sum(x^3 * SMV.omega_hat(:A, x; spherical_flow=sf)^2 for x in xs) / n
        @test isapprox(ratio, sqrt(q(true) / q(false)); rtol=1.0e-8)
        # The option slows the flow: the profile norm drops to 0.853 of the 2-D value.
        @test 0.8 < ratio < 0.9
        @test isapprox(ratio, 0.853; atol=5.0e-4)
        @test !isapprox(ratio, 2 / 3; rtol=0.05)
    end

    @testset "case set-up: configuration and input contracts" begin
        cfg = SMV.make_config(SMV.NX_REFERENCE; spherical_flow=true)
        @test cfg.geometry.spherical_flow
        @test cfg.geometry.spherical_metric
        @test cfg.solver.darcy_elimination
        @test isapprox(cfg.geometry.psurface, SMV.P_AMB; rtol=1.0e-15)
        # eta / zeta = phi / etaphikoef = 1e-6: a nearly incompressible matrix.
        @test isapprox(SMV.PHI / cfg.solver.etaphikoef, 1.0e-6; rtol=1.0e-12)
        @test isapprox(cfg.time.dt_initial, SMV.DT_YR; rtol=1.0e-15)
        off = SMV.make_config(33; spherical_flow=false, darcy_elimination=false)
        @test !off.geometry.spherical_flow
        @test !off.solver.darcy_elimination
        @test_throws ArgumentError SMV.make_config(SMV.NX_MIN - 2; spherical_flow=true)
        for n_even in (18, 34, 64)
            @test_throws ArgumentError SMV.make_config(n_even; spherical_flow=true)
        end
        @test_throws ArgumentError SMV.run_case(SMV.NX_MIN - 2; case=:A)
        @test_throws ArgumentError SMV.run_case(33; case=:D)
        @test_throws ArgumentError SMV.run_case(33; case=:C, spherical_flow=false)
        # Case C: eta / zeta = 1 through etaphikoef = phi, and p_s = 1e7 Pa.
        cfg_c = SMV.make_config(33; spherical_flow=true, case=:C)
        @test isapprox(SMV.PHI / cfg_c.solver.etaphikoef, 1.0; rtol=1.0e-12)
        @test isapprox(cfg_c.geometry.psurface, 1.0e7; rtol=1.0e-15)
        @test isapprox(cfg_c.disk.p_amb_disk, 1.0e7; rtol=1.0e-15)
        @test isapprox(cfg_c.materials.kphim0[1], 1.0e-13; rtol=1.0e-15)
        @test isapprox(cfg_c.materials.etafluidm[1], 1.0e-3; rtol=1.0e-15)
        @test isapprox(cfg_c.thermodynamics.phim0, 0.2; rtol=1.0e-15)
        @test !cfg_c.poroelasticity.hydrofracture && !cfg_c.venting.active
        @test isapprox(
            cfg.solver.etaphikoef,
            SMV.make_config(65; spherical_flow=true, case=:B).solver.etaphikoef;
            rtol=1.0e-15,
        )
        @test_throws ArgumentError SMV.run_case(33; case=:B, spherical_flow=false)
        @test_throws ArgumentError SMV.run_case(33; case=:A0, spherical_flow=false)
        @test_throws ArgumentError SMV.evaluate_spherical_momentum_benchmark(;
            nx_list=(17,)
        )
        @test_throws ArgumentError SMV.evaluate_spherical_momentum_benchmark(;
            nx_list=(33.0,)
        )
    end

    @testset "set_fields!: uniform matrix, body force and C2 before the solve (case $case)" for case in
                                                                                                SMV.CASES

        cfg = SMV.make_config(33; spherical_flow=true, case=case)
        mktempdir() do out
            state, coords, _, tio, _ = SMV._init_case(cfg, out, case)
            tio === nothing || close(tio)
            g = state.grids
            @test state.timestep == 1
            @test isapprox(state.dt, SMV.DT_YR * cfg.time.yearlength; rtol=1.0e-15)
            @test all(x -> isapprox(x, SMV.ETA_ROCK; rtol=1.0e-15), g.ETA)
            # C2: G dt >= 1e6 eta on the basic and the P nodes.
            @test maximum(g.ETA ./ (g.GGG .* state.dt)) <= 1.0e-6
            @test maximum(g.ETAP ./ (g.GGGP .* state.dt)) <= 1.0e-6
            @test all(iszero, g.S330) && all(iszero, g.SXX0) && all(iszero, g.SXY0)
            # rho g = f at every vx and vy node, 0 outside R.
            fs = SMV.force_scale(case)
            grid = SMV.grid_of(coords)
            n_out = 0
            for j in 2:(coords.Nx1 - 2), i in 2:(coords.Ny1 - 1)
                X = (coords.xvx[j] - SMV.X_CENTRE) / SMV.R_PLANET
                Y = (coords.yvx[i] - SMV.Y_CENTRE) / SMV.R_PLANET
                fx = fs * SMV.exact_force(case, X, Y)[1]
                @test isapprox(g.RHOX[i, j] * g.gx[i, j], fx; atol=1.0e-12 * fs)
                n_out += X^2 + Y^2 >= 1 ? 1 : 0
            end
            @test 100 < n_out < coords.Nx1 * coords.Ny1
            j, i = 20, 12
            X = (coords.xvy[j] - SMV.X_CENTRE) / SMV.R_PLANET
            Y = (coords.yvy[i] - SMV.Y_CENTRE) / SMV.R_PLANET
            @test isapprox(
                g.RHOY[i, j] * g.gy[i, j], fs * SMV.exact_force(case, X, Y)[2]; rtol=1.0e-14
            )
            @test isapprox(grid.dx, SMV.L_BOX / 32; rtol=1.0e-15)
            # C7: no fluid body force; the pressures start at the ambient pressure.
            @test maximum(abs, g.RHOFX) <= 1.0e-30 && maximum(abs, g.RHOFY) <= 1.0e-30
            p_s = SMV.ambient_pressure(case)
            @test all(x -> isapprox(x, p_s; rtol=1.0e-15), g.pr)
            @test all(x -> isapprox(x, p_s; rtol=1.0e-15), g.pf)
            @test isapprox(state.accumulators.P_amb, p_s; rtol=1.0e-15)
        end
    end

    @testset "error measures on exact fields" begin
        coords = Erebus.GridCoordinates(33, 33; xsize=SMV.L_BOX, ysize=SMV.L_BOX)
        grid = SMV.grid_of(coords)
        for (case, sf) in ((:A, true), (:A, false), (:A0, true), (:B, true))
            v = SMV.exact_fields(grid, case; spherical_flow=sf)
            e = SMV.velocity_errors(v.vx, v.vy, grid, case; spherical_flow=sf)
            @test isapprox(e.l2, 0.0; atol=1.0e-15)
            @test isapprox(e.linf, 0.0; atol=1.0e-15)
            # A field 1 % too large has the relative error 0.01 everywhere.
            e1 = SMV.velocity_errors(
                1.01 .* v.vx, 1.01 .* v.vy, grid, case; spherical_flow=sf
            )
            @test isapprox(e1.l2, 0.01; rtol=1.0e-10)
            @test isapprox(e1.l2_inside, 0.01; rtol=1.0e-10)
            @test isapprox(e1.linf, 0.01; rtol=1.0e-10)
            @test isapprox(e1.l2_outside, 0.0; atol=1.0e-15)
        end
        # A flow outside R counts in the full-box and the outside error only.
        v = SMV.exact_fields(grid, :B)
        vx = copy(v.vx)
        vx[2, 2] = 0.1 * maximum(abs, v.vx)
        e = SMV.velocity_errors(vx, v.vy, grid, :B)
        @test isapprox(e.l2_inside, 0.0; atol=1.0e-15)
        @test isapprox(e.l2^2, e.l2_inside^2 + e.l2_outside^2; rtol=1.0e-12)
        @test 0.0 < e.l2_outside < 1.0
        # The fit of Omega(rho) with degree 5 in rho^2 recovers Omega(0) of an exact field.
        for (case, sf) in ((:A, true), (:A, false), (:A0, true))
            v = SMV.exact_fields(grid, case; spherical_flow=sf)
            om0 = SMV.omega_hat(case, 0.0; spherical_flow=sf)
            fit = SMV.centre_fit_rotation(v.vx, v.vy, grid; scale=SMV.velocity_scale(case))
            @test isapprox(fit, om0; atol=1.0e-12)
            # The basic-node vorticity of the exact field is within (dx / R)^2 of Omega(0).
            node = SMV.centre_node_rotation(
                v.vx, v.vy, grid; scale=SMV.velocity_scale(case)
            )
            @test isapprox(node, om0; atol=0.1 * (grid.dx / SMV.R_PLANET)^2)
        end
        # A degree-2 fit has a bias: the case A profile has degree 5 in rho^2.
        v = SMV.exact_fields(grid, :A)
        fit2 = SMV.centre_fit_rotation(
            v.vx, v.vy, grid; scale=SMV.velocity_scale(:A), degree=2
        )
        @test !isapprox(fit2, 0.01; rtol=1.0e-3)
        # Each axis carries half of Omega(0): the vy row alone gives half the rotation.
        vs = SMV.velocity_scale(:A)
        @test isapprox(
            SMV.centre_fit_rotation(zero(v.vx), v.vy, grid; scale=vs), 0.005; atol=1.0e-12
        )
        @test isapprox(
            SMV.centre_fit_rotation(v.vx, zero(v.vy), grid; scale=vs), 0.005; atol=1.0e-12
        )
        node_vy = SMV.centre_node_rotation(zero(v.vx), v.vy, grid; scale=vs)
        node_all = SMV.centre_node_rotation(v.vx, v.vy, grid; scale=vs)
        @test isapprox(node_vy, 0.5 * node_all; rtol=1.0e-12)
        @test_throws ArgumentError SMV.centre_fit_rotation(
            v.vx, v.vy, grid; rho_min=0.8, rho_max=0.85
        )
        @test_throws ArgumentError SMV.centre_fit_rotation(
            v.vx, v.vy, grid; rho_min=0.5, rho_max=0.4
        )
        # The window must lie in [0, 1], finite, and the degree must be non-negative.
        for (lo, hi) in ((0.0, 0.85), (0.18, 1.3), (0.18, Inf), (NaN, 0.85))
            @test_throws ArgumentError SMV.centre_fit_rotation(
                v.vx, v.vy, grid; rho_min=lo, rho_max=hi
            )
        end
        @test_throws ArgumentError SMV.centre_fit_rotation(v.vx, v.vy, grid; degree=-1)
        @test_throws ArgumentError SMV.centre_fit_rotation(v.vx, v.vy, grid; degree=-3)
        # The 2 axes sample the same radii: 16 nodes on 4 radii cannot fit 6 coefficients.
        @test_throws ArgumentError SMV.centre_fit_rotation(
            v.vx, v.vy, grid; rho_min=0.55, rho_max=0.85
        )
        coarse = SMV.grid_of(
            Erebus.GridCoordinates(17, 17; xsize=SMV.L_BOX, ysize=SMV.L_BOX)
        )
        vc = SMV.exact_fields(coarse, :A)
        @test_throws ArgumentError SMV.centre_fit_rotation(vc.vx, vc.vy, coarse; scale=vs)
        # The velocity unit must be finite and positive.
        for s in (0.0, -1.0, Inf, NaN)
            @test_throws DomainError SMV.centre_fit_rotation(v.vx, v.vy, grid; scale=s)
            @test_throws DomainError SMV.centre_node_rotation(v.vx, v.vy, grid; scale=s)
        end
        # The smallest grid of the case resolves the default window with no bias.
        nmin = SMV.NX_MIN
        gmin = SMV.grid_of(
            Erebus.GridCoordinates(nmin, nmin; xsize=SMV.L_BOX, ysize=SMV.L_BOX)
        )
        vmin = SMV.exact_fields(gmin, :A)
        @test isapprox(
            SMV.centre_fit_rotation(vmin.vx, vmin.vy, gmin; scale=vs), 0.01; atol=1.0e-12
        )
        gbelow = SMV.grid_of(
            Erebus.GridCoordinates(nmin - 2, nmin - 2; xsize=SMV.L_BOX, ysize=SMV.L_BOX)
        )
        vbelow = SMV.exact_fields(gbelow, :A)
        @test_throws ArgumentError SMV.centre_fit_rotation(
            vbelow.vx, vbelow.vy, gbelow; scale=vs
        )
        # The centre of an even grid lies on a P node: no basic node at the centre.
        even = SMV.grid_of(Erebus.GridCoordinates(34, 34; xsize=SMV.L_BOX, ysize=SMV.L_BOX))
        ve = SMV.exact_fields(even, :A)
        @test_throws ArgumentError SMV.centre_node_rotation(ve.vx, ve.vy, even)
        # S33 of the exact flow has the error 0; an S33 with half the prefactor 1/2.
        s33 = SMV.exact_s33(grid, :B)
        @test isapprox(SMV.s33_error(s33, grid), 0.0; atol=1.0e-15)
        @test isapprox(SMV.s33_error(0.5 .* s33, grid), 0.5; rtol=1.0e-12)
        @test 1.0e4 < maximum(abs, s33) < 1.0e6
        # The pressure deviation skips the anchor ring and reads the interior P nodes.
        pr = fill(SMV.P_AMB, grid.Ny1, grid.Nx1)
        pr[2, 2] += 1.0e5
        @test isapprox(SMV.pressure_deviation(pr, grid), 0.0; atol=1.0e-15)
        pr[17, 17] -= 250.0
        @test isapprox(SMV.pressure_deviation(pr, grid), 250.0; rtol=1.0e-12)
        # The profile norm counts the nodes inside R.
        va = SMV.exact_fields(grid, :A)
        pa = SMV.profile_norm(va.vx, va.vy, grid)
        @test isapprox(SMV.profile_norm(2 .* va.vx, 2 .* va.vy, grid), 2 * pa; rtol=1.0e-14)
        # A velocity outside R leaves it unchanged.
        vx_out = copy(va.vx)
        vx_out[2, 2] = 10 * maximum(abs, va.vx)
        @test isapprox(SMV.profile_norm(vx_out, va.vy, grid), pa; rtol=1.0e-14)
        # The vx and vy parts add in quadrature, and each carries half of the flow.
        px = SMV.profile_norm(va.vx, zero(va.vy), grid)
        py = SMV.profile_norm(zero(va.vx), va.vy, grid)
        @test isapprox(px^2 + py^2, pa^2; rtol=1.0e-14)
        @test isapprox(px, py; rtol=1.0e-12)
        # Case C fields at the P nodes: the exact fields have the error 0, a field 1 % too
        # large the error 0.01, and a value outside R counts only in the outside error.
        ex = SMV.exact_pnode_fields(grid, :C)
        for f in (:s33, :pt, :pf, :eii, :sii, :hs)
            e0 = SMV.pnode_errors(getfield(ex, f), getfield(ex, f), grid)
            @test isapprox(e0.l2_inside, 0.0; atol=1.0e-15) &&
                isapprox(e0.out, 0.0; atol=1.0e-15)
            e1 = SMV.pnode_errors(1.01 .* getfield(ex, f), getfield(ex, f), grid)
            @test isapprox(e1.l2_inside, 0.01; rtol=1.0e-10)
        end
        p_s = SMV.P_AMB_C
        pt = p_s .+ ex.pt
        pt[2, 2] += 1.0e-3 * maximum(abs, ex.pt)
        e_out = SMV.pnode_errors(pt, ex.pt, grid; offset=p_s)
        @test isapprox(e_out.l2_inside, 0.0; atol=1.0e-12)
        @test isapprox(e_out.out, 1.0e-3; rtol=1.0e-6)
        # The heating peaks at the centre with 720 eta U0^2 / R^2, about 3e-5 W/m^3.
        @test 1.0e-5 < maximum(ex.hs) < 1.0e-4
        # The anchor rings read the P nodes next to the box edge only.
        pf = fill(p_s, grid.Ny1, grid.Nx1)
        @test isapprox(SMV.anchor_deviation(pf, pf, grid, p_s), 0.0; atol=1.0e-15)
        pf[17, 17] += 1.0e6
        @test isapprox(SMV.anchor_deviation(pf, pf, grid, p_s), 0.0; atol=1.0e-15)
        pf[2, 17] += 1.0e3
        @test isapprox(SMV.anchor_deviation(pf, pf, grid, p_s), 1.0e-4; rtol=1.0e-9)
    end

    @testset "a rigid rotation has a discrete out-of-plane strain rate next to the centre" begin
        # Omega(0) != 0 is a rigid rotation near the centre: no strain in the continuum, but
        # the discrete eps_33 of the continuity and momentum rows and of the stress update
        # is +-0.011 Omega at the 8 P nodes at r = sqrt(2.5) dx = 1.58 dx on every grid.
        function rotation_eps33(n)
            coords = Erebus.GridCoordinates(n, n; xsize=SMV.L_BOX, ysize=SMV.L_BOX)
            Ny1, Nx1 = coords.Ny1, coords.Nx1
            om, eta = 1.0e-12, SMV.ETA_ROCK
            P() = zeros(Ny1, Nx1)
            B() = zeros(Ny1 - 1, Nx1 - 1)
            g = (;
                vx=[-om * (coords.yvx[i] - SMV.Y_CENTRE) for i in 1:Ny1, j in 1:Nx1],
                vy=[om * (coords.xvy[j] - SMV.X_CENTRE) for i in 1:Ny1, j in 1:Nx1],
                ETAP=fill(eta, Ny1, Nx1),
                GGGP=fill(1.0e40, Ny1, Nx1),
                S330=P(),
                S33=P(),
                DS33=P(),
                EXX=P(),
                EXY=B(),
                SXX=P(),
                SXY=B(),
                EII=P(),
                SII=P(),
            )
            Erebus._spherical_stress_update!(
                g, coords, 1.0; xcenter=SMV.X_CENTRE, ycenter=SMV.Y_CENTRE
            )
            # S33 = (2 eta / 3) (2 eps_33 - eps_xx - eps_yy) with eps_xx = eps_yy = 0.
            near, r_near, far = 0.0, 0.0, 0.0
            inner = Tuple{Float64,Float64}[]
            e33_all = g.S33 ./ (4 / 3 * eta * om)
            for j in 2:(Nx1 - 1), i in 2:(Ny1 - 1)
                e33 = abs(e33_all[i, j])
                dxp, dyp = coords.xp[j] - SMV.X_CENTRE, coords.yp[i] - SMV.Y_CENTRE
                r = hypot(dxp, dyp)
                e33 > near && ((near, r_near) = (e33, r / coords.dx))
                r > 0.5 * SMV.R_PLANET && (far = max(far, e33))
                isapprox(r / coords.dx, sqrt(2.5); rtol=1.0e-9) &&
                    push!(inner, (atan(dyp, dxp), e33_all[i, j]))
            end
            return (; near, r_near, far, inner=sort!(inner), e33_all)
        end
        e33, e65 = rotation_eps33(33), rotation_eps33(65)
        @test isapprox(e33.near, 0.011; rtol=0.01)
        @test isapprox(e65.near, e33.near; rtol=1.0e-10)
        @test isapprox(e33.r_near, sqrt(2.5); rtol=1.0e-12)
        @test isapprox(e65.r_near, sqrt(2.5); rtol=1.0e-12)
        # The 8 nodes at r = sqrt(2.5) dx carry +-0.011, alternating around the centre
        # as sin(4 theta) does, so they sum to 0.
        vals = last.(e33.inner)
        @test length(vals) == 8
        @test all(isapprox(abs(v), e33.near; rtol=1.0e-10) for v in vals)
        @test all(vals[k] * vals[mod1(k + 1, 8)] < 0 for k in 1:8)
        @test isapprox(sum(vals), 0.0; atol=1.0e-12 * e33.near)
        # Away from the centre it falls at second order.
        @test isapprox(log2(e33.far / e65.far), 2.0; atol=0.35)
        # The rotation is linear, so eps_33 depends only on the offset from the centre in
        # cells: every interior P node of Nx = 33, the edge ring included, matches Nx = 65.
        a33 = e33.e33_all
        n33 = size(a33, 1)
        a65_mid = e65.e33_all[(2:(n33 - 1)) .+ 16, (2:(n33 - 1)) .+ 16]
        @test isapprox(a33[2:(n33 - 1), 2:(n33 - 1)], a65_mid; rtol=1.0e-10)
        @test maximum(abs, a33[2, 2:(n33 - 1)]) > 1.0e-3 * e33.near
    end

    # Runs at Nx = 33 and 65 only. The Nx = 129 and order metrics take their recorded
    # values, which `validate_record_dict` checks against their reference and tolerance.
    bench_runs = SMV.evaluate_spherical_momentum_benchmark(; nx_list=(33, 65))
    runs = bench_runs.runs
    measured = bench_runs.metrics
    a0 = Dict(n => SMV.run_case(n; case=:A0) for n in (33, 65))
    record_path = joinpath(
        @__DIR__, "..", "benchmarks", "records", "spherical_momentum.toml"
    )
    raw_record_all = TOML.parsefile(record_path)
    recorded = Dict(
        String(m["id"]) => Float64(m["value"]) for m in raw_record_all["metrics"]
    )

    @testset "conditions C1, C2, C3 and C6 after the solve ($(key[2]), on = $(key[3]), $(key[4]), Nx = $(key[1]))" for key in
                                                                                                                       sort!(
        collect(keys(runs)); by=string
    )
        r = runs[key]
        # C6: one plastic pass, converged, and the solve keeps the time step.
        @test r.converged
        @test r.n_plastic_passes == 1
        @test r.dt_dev <= 1.0e-12
        # C1: eta / zeta = 1e-6 (case C: 1) on every interior P node.
        ratio = SMV.eta_zeta_ratio(key[2])
        @test isapprox(r.eta_zeta, ratio; rtol=1.0e-12)
        @test isapprox(r.eta_zeta_min, ratio; rtol=1.0e-12)
        # C2: G dt >= 1e6 eta on the basic and the P nodes.
        @test r.elastic <= 1.0e-6
        @test r.elastic_p <= 1.0e-6
        # C3: the solve keeps the body force.
        @test iszero(r.g_dev)
        @test all(isfinite, r.vx) && all(isfinite, r.vy) && all(isfinite, r.pr)
        @test all(isfinite, r.S33) && all(isfinite, r.pf)
        # No fluid body force after the solve, and the anchor rings of P_t and P_f hold
        # the ambient pressure.
        @test r.rhof_max <= 1.0e-30
        @test r.anchor_dev <= 1.0e-9
    end

    @testset "conditions hold for the A0 runs (Nx = $n)" for n in (33, 65)
        r = a0[n]
        @test r.converged && r.n_plastic_passes == 1
        @test r.dt_dev <= 1.0e-12
        @test r.eta_zeta <= 1.0e-6 * (1 + 1.0e-12)
        @test r.elastic <= 1.0e-6 && r.elastic_p <= 1.0e-6
        @test iszero(r.g_dev)
        @test all(isfinite, r.vx) && all(isfinite, r.vy) && all(isfinite, r.pr)
    end

    @testset "case A velocity converges at second order (spherical_flow = $sf)" for sf in (
        true, false
    )
        tag = sf ? "a_on" : "a_off"
        l33, l65 = runs[(33, :A, sf, :base)].l2, runs[(65, :A, sf, :base)].l2
        @test l33 <= 0.053
        @test l65 <= (sf ? 0.0122 : 0.012)
        @test isapprox(log2(l33 / l65), 2.0; atol=0.35)
        @test isapprox(measured["l2_$(tag)_nx65"], l65; rtol=1.0e-14)
        # The error lies inside R: the flow outside is 0 in both forms.
        r = runs[(65, :A, sf, :base)]
        @test r.l2_outside < r.l2_inside
        # Scale guard: the largest velocity is a few 1e-12 m/s, as the closed form.
        @test 1.0e-12 < r.v_max < 1.0e-11
    end

    @testset "case A L2 error at Nx = 65 with the option (pinned)" begin
        l65 = runs[(65, :A, true, :base)].l2
        @test isapprox(l65, 8.6724558541e-3; rtol=1.0e-6)
        @test 0.005 < l65 < 0.012
    end

    @testset "case A centre rotation: fit at second order, basic node slow with the option" begin
        om_on = SMV.omega_hat(:A, 0.0; spherical_flow=true)
        om_off = SMV.omega_hat(:A, 0.0; spherical_flow=false)
        fit(n, sf) = runs[(n, :A, sf, :base)].centre_fit
        node(n, sf) = runs[(n, :A, sf, :base)].centre_node
        for (sf, om) in ((true, om_on), (false, om_off))
            e33, e65 = fit(33, sf) / om - 1, fit(65, sf) / om - 1
            # The fit lies above Omega(0) and falls at order 2.06 in both settings.
            @test 0.0 < e65 < e33 < 0.05
            @test isapprox(log2(e33 / e65), 2.0; atol=0.35)
            @test isapprox(log2(e33 / e65), 2.06; atol=0.01)
        end
        # The basic-node value with the option falls at order 0.45 from Nx = 33 to 65 and
        # at order 0.74 from 65 to 129 in the record; without the option it falls at order 2.
        n33, n65 = abs(node(33, true) / om_on - 1), abs(node(65, true) / om_on - 1)
        @test 0.2 < log2(n33 / n65) < 1.0
        @test isapprox(log2(n33 / n65), 0.45; atol=0.01)
        rec_node_order = log2(
            recorded["centre_node_a_on_nx65"] / recorded["centre_node_a_on_nx129"]
        )
        @test isapprox(rec_node_order, 0.74; atol=0.01)
        @test n65 <= 0.01
        d33, d65 = abs(node(33, false) / om_off - 1), abs(node(65, false) / om_off - 1)
        @test 1.65 < log2(d33 / d65) < 2.5
        # The ratio of the fitted centre values is 2/3 within 0.08 %.
        ratio = measured["centre_fit_ratio_nx65"]
        @test isapprox(ratio, 2 / 3; rtol=0.0008)
        @test isapprox(ratio, fit(65, true) / fit(65, false); rtol=1.0e-14)
        @test !isapprox(ratio, 1.0; rtol=0.1)
    end

    @testset "case A0: Omega(0) = 0 converges at second order in L2, Linf and pressure" begin
        r33, r65 = a0[33], a0[65]
        @test isapprox(log2(r33.l2 / r65.l2), 2.0; atol=0.35)
        @test isapprox(log2(r33.linf / r65.linf), 2.0; atol=0.35)
        @test 1.65 < log2(r33.dp_max / r65.dp_max) < 2.6
        # Case A with Omega(0) != 0 keeps a pressure error of about 1100 Pa at the centre.
        a33, a65 = runs[(33, :A, true, :base)].dp_max, runs[(65, :A, true, :base)].dp_max
        @test 0.0 < log2(a33 / a65) < 0.5
        @test 900.0 < a65 < 1300.0
        @test isapprox(a65, 1100.0; atol=20.0)
        @test r65.dp_max < 0.15 * a65
        # Without the option the same pressure error falls with the grid.
        o33, o65 = runs[(33, :A, false, :base)].dp_max, runs[(65, :A, false, :base)].dp_max
        @test 1.65 < log2(o33 / o65) < 2.5
    end

    @testset "case B velocity and S33 converge at second order" begin
        b33, b65 = runs[(33, :B, true, :base)], runs[(65, :B, true, :base)]
        @test b65.l2 <= 0.0138
        @test b65.l2_inside <= 0.0138
        @test isapprox(log2(b33.l2 / b65.l2), 2.0; atol=0.35)
        @test isapprox(log2(b33.l2_inside / b65.l2_inside), 2.0; atol=0.35)
        # grids.S33 from the stress update against tau_33 = 2 eta u_r / r.
        @test b65.s33_l2 <= 0.0152
        @test isapprox(log2(b33.s33_l2 / b65.s33_l2), 2.0; atol=0.35)
        exact = SMV.exact_s33(b65.grid, :B)
        i, j = argmax(abs.(exact)).I
        @test sign(b65.S33[i, j]) == sign(exact[i, j])
        @test isapprox(b65.S33[i, j], exact[i, j]; rtol=0.05)
        # Scale guard: tau_33 is O(eta U0 / R) = 2e6 Pa times 0.11.
        @test 1.0e5 < maximum(abs, b65.S33) < 1.0e6
    end

    @testset "case B at Nx = 65: pinned errors and the S33 projection" begin
        b65 = runs[(65, :B, true, :base)]
        @test isapprox(b65.l2, 9.828284022e-3; rtol=1.0e-5)
        @test isapprox(b65.s33_l2, 1.0838697880e-2; rtol=1.0e-5)
        # The projection of S33 on tau_33 inside R moves linearly with a wrong prefactor.
        exact = SMV.exact_s33(b65.grid, :B)
        num = den = 0.0
        for j in 2:(b65.grid.Nx1 - 1), i in 2:(b65.grid.Ny1 - 1)
            X, Y = SMV._nd(b65.grid.xp[j], b65.grid.yp[i])
            X^2 + Y^2 < 1 || continue
            num += b65.S33[i, j] * exact[i, j]
            den += exact[i, j]^2
        end
        @test isapprox(num / den, 0.99846; rtol=5.0e-4)
        @test !isapprox(num / den, 1.0; rtol=1.0e-3)
    end

    @testset "S33 outside R, on the ghost rows and finite (case B, Nx = $n)" for (
        n, bound
    ) in (
        (33, 2.5e-3), (65, 5.4e-4)
    )
        r = runs[(n, :B, true, :base)]
        S = r.S33
        Ny1, Nx1 = size(S)
        peak = maximum(abs, SMV.exact_s33(r.grid, :B))
        @test all(isfinite, S)
        # tau_33 = 0 for r >= R; the solved S33 there is a grid error that falls with Nx.
        s_out = maximum(
            abs(S[i, j]) for j in 2:(Nx1 - 1), i in 2:(Ny1 - 1) if
            sum(abs2, SMV._nd(r.grid.xp[j], r.grid.yp[i])) >= 1
        )
        @test s_out <= bound * peak
        # The ghost rows and columns copy their interior neighbours.
        @test isequal(S[1, 2:(Nx1 - 1)], S[2, 2:(Nx1 - 1)])
        @test isequal(S[Ny1, 2:(Nx1 - 1)], S[Ny1 - 1, 2:(Nx1 - 1)])
        @test isequal(S[:, 1], S[:, 2])
        @test isequal(S[:, Nx1], S[:, Nx1 - 1])
        @test !iszero(S[Ny1 - 1, 2:(Nx1 - 1)]) && !iszero(S[2:(Ny1 - 1), 2])
    end

    @testset "case C conditions C1 to C8 after the solve (Nx = $n)" for n in (33, 65)
        r = runs[(n, :C, true, :base)]
        # C1: eta / zeta = 1 on every interior P node; C2: viscous limit.
        @test isapprox(r.eta_zeta, 1.0; rtol=1.0e-12)
        @test isapprox(r.eta_zeta_min, 1.0; rtol=1.0e-12)
        @test r.elastic <= 1.0e-6 && r.elastic_p <= 1.0e-6
        # C3: the solve keeps g; C6: the time step and one plastic pass without yield.
        @test iszero(r.g_dev)
        @test r.dt_dev <= 1.0e-12
        @test r.converged && r.n_plastic_passes == 1
        # C7: no fluid body force; C8: the anchor rings hold p_s = 1e7 Pa.
        @test r.rhof_max <= 1.0e-30
        @test r.anchor_dev <= 1.0e-9
        # P_t > 0 everywhere in the solved field, as in the exact one.
        @test 0.0 < minimum(r.pr[2:(end - 1), 2:(end - 1)]) < SMV.P_AMB_C
        @test maximum(r.pr[2:(end - 1), 2:(end - 1)]) > 4 * SMV.P_AMB_C
    end

    @testset "case C velocity, stresses, pressures and heating converge at second order" begin
        c33, c65 = runs[(33, :C, true, :base)], runs[(65, :C, true, :base)]
        for f in (:l2, :l2_inside, :s33_l2, :pt_l2, :pf_l2, :eii_l2, :sii_l2, :hs_l2)
            e33, e65 = getfield(c33, f), getfield(c65, f)
            @test isapprox(log2(e33 / e65), 2.0; atol=0.35)
            @test e65 <= 0.01
        end
        # The flow is 0 outside R: the full-box error is close to the error inside.
        @test c65.l2_outside < 0.2 * c65.l2_inside
        # S_II = 2 eta E_II in the viscous limit, so their errors agree closely.
        @test isapprox(c65.sii_l2, c65.eii_l2; rtol=0.05)
        # Sign and scale guards: S33 has the sign of tau_33 where it peaks, and tau_33
        # peaks at r = R / 2 with -5.6 eta U0 / R = -1.1e7 Pa.
        exact = SMV.exact_s33(c65.grid, :C)
        i, j = argmax(abs.(exact)).I
        @test sign(c65.S33[i, j]) == sign(exact[i, j])
        @test isapprox(c65.S33[i, j], exact[i, j]; rtol=0.05)
        @test isapprox(SMV.exact_tau33(:C, 0.5, 0.0), -45 / 8; rtol=1.0e-12)
        @test 5.0e6 < maximum(abs, c65.S33) < 2.0e7
        @test 1.0e-10 < maximum(abs, c65.vx) < 1.0e-8
    end

    @testset "case C: the S33 error next to the centre falls in share with the grid" begin
        # Share of the squared S33 error inside R on the 32 P nodes with r <= 3 dx.
        function near_share(r)
            g = r.grid
            ex = SMV.exact_s33(g, :C)
            near = total = 0.0
            n_near = 0
            for j in 2:(g.Nx1 - 1), i in 2:(g.Ny1 - 1)
                X, Y = SMV._nd(g.xp[j], g.yp[i])
                X^2 + Y^2 < 1 || continue
                d2 = (r.S33[i, j] - ex[i, j])^2
                total += d2
                if hypot(X, Y) * SMV.R_PLANET <= 3 * g.dx
                    near += d2
                    n_near += 1
                end
            end
            return near / total, n_near
        end
        s33, n33 = near_share(runs[(33, :C, true, :base)])
        s65, n65 = near_share(runs[(65, :C, true, :base)])
        @test n33 == 32 && n65 == 32
        @test isapprox(s33, 0.57; atol=0.01)
        @test isapprox(s65, 0.22; atol=0.01)
    end

    @testset "case C at Nx = 65: pinned velocity error" begin
        c65 = runs[(65, :C, true, :base)]
        @test isapprox(c65.l2_inside, 4.3581222816e-3; rtol=1.0e-6)
        # The flow without the trace in the momentum rows has the error 0.10 at Nx = 65.
        @test c65.l2_inside < 0.05 * 0.10
        @test isapprox(measured["l2_inside_c_nx65"], c65.l2_inside; rtol=1.0e-14)
    end

    @testset "case C outside R: S33, P_t and P_f stay at the grid error (Nx = $n)" for (
        n, bounds
    ) in (
        (33, (1.5e-3, 2.1e-4, 1.1e-3)), (65, (3.7e-4, 6.3e-5, 3.0e-4))
    )
        # Largest error for r >= R over the largest exact value inside, about 1.4 times
        # the measured value.
        r = runs[(n, :C, true, :base)]
        @test r.s33_out <= bounds[1]
        @test r.pt_out <= bounds[2]
        @test r.pf_out <= bounds[3]
        @test all(isfinite, (r.s33_out, r.pt_out, r.pf_out))
    end

    @testset "case C outside R: the errors fall with the grid" begin
        c33, c65 = runs[(33, :C, true, :base)], runs[(65, :C, true, :base)]
        for f in (:s33_out, :pt_out, :pf_out)
            @test getfield(c65, f) < getfield(c33, f) / 3
            @test getfield(c65, f) < 1.0e-3
        end
    end

    @testset "on/off L2-norm ratio of the case A profile against its exact value" begin
        exact = SMV.profile_norm_ratio()
        @test isapprox(exact, sqrt(8 / 11); rtol=1.0e-15)
        r65 = measured["profile_ratio_nx65"]
        @test isapprox(r65, exact; rtol=0.00075)
        @test !isapprox(r65, 1.0; rtol=0.05)
        @test isapprox(
            r65,
            runs[(65, :A, true, :base)].norm / runs[(65, :A, false, :base)].norm;
            rtol=1.0e-14,
        )
    end

    @testset "6-variable path equals the 4-variable path (case $case)" for case in
                                                                           (:A, :B, :C)
        tag = lowercase(String(case))
        six, four = runs[(33, case, true, :elim_off)], runs[(33, case, true, :base)]
        @test measured["elim_$(tag)_nx33"] <= 1.0e-7
        @test maximum(abs, six.vx .- four.vx) / maximum(abs, four.vx) <= 1.0e-7
        @test maximum(abs, six.pr .- four.pr) / SMV.stress_scale(case) <= 1.0e-7
        @test maximum(abs, six.pf .- four.pf) / SMV.stress_scale(case) <= 1.0e-7
        other = case == :A ? :B : :A
        @test_throws ArgumentError SMV.elimination_difference(
            six, runs[(33, other, true, :base)]
        )
    end

    recorded_only = filter(kv -> !haskey(measured, kv.first), recorded)
    bench = (; runs, metrics=merge(measured, recorded_only))

    @testset "metrics: 37 from the runs, 33 recorded values for Nx = 129 and the orders" begin
        @test length(measured) == 37
        @test length(recorded_only) == 33
        @test all(endswith(k, "129") for k in keys(recorded_only))
        @test Set(keys(bench.metrics)) == Set(SMV.EXPECTED_METRIC_IDS)
    end

    @testset "recorded orders from Nx = 65 to 129 match their recorded errors" begin
        order_ids = filter(k -> occursin("_order_", k), SMV.EXPECTED_METRIC_IDS)
        @test length(order_ids) == 15
        for oid in order_ids
            base_id = replace(oid, "_order" => "", "_65_129" => "")
            e65, e129 = recorded[base_id * "_nx65"], recorded[base_id * "_nx129"]
            @test isapprox(recorded[oid], log2(e65 / e129); rtol=1.0e-12)
        end
    end

    @testset "metrics_from_runs with every grid of NX_LIST" begin
        om0(sf) = SMV.omega_hat(:A, 0.0; spherical_flow=sf)
        fake(case, sf, e) = (;
            case,
            l2=e,
            l2_inside=0.9e,
            s33_l2=1.1e,
            pt_l2=1.2e,
            pf_l2=1.3e,
            eii_l2=1.4e,
            sii_l2=1.5e,
            hs_l2=1.6e,
            centre_fit=case == :A ? om0(sf) * (1 + e) : NaN,
            centre_node=case == :A ? om0(sf) * (1 - e / 2) : NaN,
            norm=sf ? 0.85 : 1.0,
            vx=fill(1.0, 2, 2),
            vy=fill(-2.0, 2, 2),
            pr=fill(SMV.P_AMB, 2, 2),
            pf=fill(SMV.P_AMB, 2, 2),
        )
        fake_runs = Dict{Tuple,Any}()
        for (n, e) in ((33, 0.04), (65, 0.01), (129, 0.0025))
            fake_runs[(n, :A, true, :base)] = fake(:A, true, e)
            fake_runs[(n, :A, false, :base)] = fake(:A, false, 1.5e)
            fake_runs[(n, :B, true, :base)] = fake(:B, true, 2e)
            fake_runs[(n, :C, true, :base)] = fake(:C, true, 3e)
        end
        for case in (:A, :B, :C)
            fake_runs[(33, case, true, :elim_off)] = fake(case, true, 0.04)
        end
        all_metrics = SMV.metrics_from_runs(fake_runs)
        @test Set(keys(all_metrics)) == Set(SMV.EXPECTED_METRIC_IDS)
        # The errors quarter from Nx = 65 to 129, so every order is 2.
        for mid in SMV.EXPECTED_METRIC_IDS
            occursin("_order_", mid) || continue
            @test isapprox(all_metrics[mid], 2.0; rtol=1.0e-12)
        end
        @test isapprox(all_metrics["centre_fit_a_on_nx65"], 0.01; rtol=1.0e-12)
        @test isapprox(all_metrics["centre_node_a_on_nx65"], 0.005; rtol=1.0e-12)
        @test isapprox(
            all_metrics["centre_fit_ratio_nx129"], (2 / 3) * 1.0025 / 1.00375; rtol=1.0e-12
        )
        @test isapprox(all_metrics["profile_ratio_nx65"], 0.85; rtol=1.0e-14)
        @test isapprox(all_metrics["elim_b_nx33"], 0.0; atol=1.0e-15)
        @test isapprox(all_metrics["s33_b_nx65"], 0.022; rtol=1.0e-14)
        @test isapprox(all_metrics["pf_c_nx65"], 0.039; rtol=1.0e-14)
        @test isapprox(all_metrics["hs_c_nx129"], 0.012; rtol=1.0e-14)
        @test isapprox(all_metrics["l2_inside_c_nx33"], 0.108; rtol=1.0e-14)
        @test isapprox(all_metrics["elim_c_nx33"], 0.0; atol=1.0e-15)
        delete!(fake_runs, (129, :A, true, :base))
        delete!(fake_runs, (33, :B, true, :elim_off))
        partial = SMV.metrics_from_runs(fake_runs)
        @test !haskey(partial, "l2_order_a_on_65_129")
        @test !haskey(partial, "centre_fit_ratio_nx129")
        @test !haskey(partial, "elim_b_nx33")
        @test haskey(partial, "l2_order_a_off_65_129")
    end

    @testset "elimination difference: the largest of the 4 field differences" begin
        base_fields(case) = (;
            case,
            vx=fill(2.0, 3, 3),
            vy=fill(-4.0, 3, 3),
            pr=fill(SMV.P_AMB, 3, 3),
            pf=fill(SMV.P_AMB, 3, 3),
        )
        four = base_fields(:B)
        @test iszero(SMV.elimination_difference(base_fields(:B), four))
        sigma = SMV.stress_scale(:B)
        # P_f is scaled by the stress unit as P_t, not by P_f - p_s, which is 0 in A, B.
        for (field, delta, expected) in (
            (:vx, 0.3, 0.15),
            (:vy, -0.4, 0.1),
            (:pr, 0.2 * sigma, 0.2),
            (:pf, 0.3 * sigma, 0.3),
        )
            six = base_fields(:B)
            getfield(six, field)[2, 3] += delta
            @test isapprox(SMV.elimination_difference(six, four), expected; rtol=1.0e-12)
        end
        # The pressure scale of case A is F0 R = 5e6 Pa.
        six = base_fields(:A)
        six.pr[1, 1] += 5.0e6
        @test isapprox(SMV.elimination_difference(six, base_fields(:A)), 1.0; rtol=1.0e-12)
    end

    @testset "the off run fails the on bands of case A by more than twice the band" begin
        off = runs[(33, :A, false, :base)]
        as_on = SMV.velocity_errors(off.vx, off.vy, off.grid, :A; spherical_flow=true)
        om_on = SMV.omega_hat(:A, 0.0; spherical_flow=true)
        shifted = Dict(
            "l2_a_on_nx33" => as_on.l2,
            "centre_fit_a_on_nx33" => abs(off.centre_fit / om_on - 1),
            "centre_node_a_on_nx33" => abs(off.centre_node / om_on - 1),
        )
        for (mid, v) in shifted
            band = SMV.PHYSICS_REGRESSION_BANDS[mid].atol
            @test abs(v - recorded[mid]) > 2 * band
            moved = copy(bench.metrics)
            moved[mid] = v
            ok, errs = SMV.check_physics_regression_bands(moved, record_path)
            @test !ok
            @test any(occursin(mid, e) for e in errs)
        end
        # Against the solution with the option the off run has the L2 error 0.63.
        @test 0.5 < as_on.l2 < 0.8
        @test isapprox(as_on.l2, 0.63; atol=0.01)
        @test 0.4 < abs(off.centre_fit / om_on - 1) < 0.6
    end

    @testset "Validation record schema and drift guards" begin
        @test isfile(record_path)
        raw_record = TOML.parsefile(record_path)
        valid_ok, valid_errors = SMV.validate_record_dict(
            raw_record; require_doc_pages=true
        )
        @test valid_ok
        @test isempty(valid_errors)
        rec_ids = [String(m["id"]) for m in raw_record["metrics"]]
        @test rec_ids == collect(SMV.EXPECTED_METRIC_IDS)
        @test all(String(m["units"]) == "dimensionless" for m in raw_record["metrics"])
        @test Set(keys(SMV.PHYSICS_REGRESSION_BANDS)) == Set(SMV.EXPECTED_METRIC_IDS)
        missing_key = deepcopy(raw_record)
        delete!(missing_key["config"], "eta_zeta_ratio")
        ok_missing, errs_missing = SMV.validate_record_dict(
            missing_key; require_doc_pages=false, check_doc_sync=false
        )
        @test !ok_missing
        @test any(occursin("eta_zeta_ratio", e) for e in errs_missing)
        missing_c = deepcopy(raw_record)
        delete!(missing_c["config"], "eta_zeta_ratio_c")
        ok_missing_c, errs_missing_c = SMV.validate_record_dict(
            missing_c; require_doc_pages=false, check_doc_sync=false
        )
        @test !ok_missing_c
        @test any(occursin("eta_zeta_ratio_c", e) for e in errs_missing_c)
        # The 37 measured metrics lie within their bands around the recorded values.
        for (mid, value) in measured
            band = SMV.PHYSICS_REGRESSION_BANDS[mid]
            @test isapprox(value, recorded[mid]; atol=band.atol, rtol=band.rtol)
        end
        ok_full, errs_full = SMV.check_physics_regression_bands(bench.metrics, record_path)
        @test ok_full
        @test isempty(errs_full)
        target = "l2_b_nx65"
        moved = copy(bench.metrics)
        moved[target] += 2 * SMV.PHYSICS_REGRESSION_BANDS[target].atol
        ok_moved, errs_moved = SMV.check_physics_regression_bands(moved, record_path)
        @test !ok_moved
        @test occursin(target, only(errs_moved))
        missing_metric = copy(bench.metrics)
        delete!(missing_metric, target)
        ok_gone, errs_gone = SMV.check_physics_regression_bands(missing_metric, record_path)
        @test !ok_gone
        @test any(occursin("Measured metrics missing '$target'", e) for e in errs_gone)
        no_metric = deepcopy(raw_record)
        filter!(m -> m["id"] != target, no_metric["metrics"])
        mktempdir() do tmpdir
            path = write_momentum_toml(joinpath(tmpdir, "short.toml"), no_metric)
            ok_short, errs_short = SMV.check_physics_regression_bands(bench.metrics, path)
            @test !ok_short
            @test any(occursin("record missing metric '$target'", e) for e in errs_short)
        end
        @test_throws DomainError SMV.check_physics_regression_bands(bench.metrics, "")
        @test !first(SMV.check_physics_regression_bands(bench.metrics, record_path * ".x"))
    end

    @testset "Record matches the configuration and the reference values" begin
        raw_record = TOML.parsefile(record_path)
        cfg = SMV.make_config(SMV.NX_REFERENCE; spherical_flow=true)
        sha = SMV.canonical_config_sha256(cfg)
        @test sha == raw_record["metadata"]["config_sha256"]
        @test length(sha) == 64
        @test SMV.canonical_config_sha256(
            SMV.make_config(SMV.NX_REFERENCE; spherical_flow=false)
        ) != sha
        text = SMV.canonical_config_text(cfg)
        for name in SMV.HASHED_CASE_CONSTANTS
            @test occursin("case.$(name)=", text)
        end
        @test length(SMV.HASHED_CASE_CONSTANTS) == 34
        @test all(
            n in SMV.HASHED_CASE_CONSTANTS for n in (:ETA_RATIO_C, :P_AMB_C, :C_POWER)
        )
        @test bytes2hex(SHA.sha256(text)) == sha
        summary = SMV.expected_config_summary(cfg)
        rec_cfg = raw_record["config"]
        @test Set(keys(summary)) == Set(keys(rec_cfg))
        for (key, expected) in summary
            if expected isa AbstractFloat
                @test isapprox(Float64(rec_cfg[key]), expected; rtol=1.0e-12, atol=0.0)
            else
                @test rec_cfg[key] == expected
            end
        end
        @test isapprox(summary["eta_zeta_ratio"], 1.0e-6; rtol=1.0e-12)
        @test isapprox(summary["eta_zeta_ratio_c"], 1.0; rtol=1.0e-12)
        @test isapprox(summary["P_amb_c_Pa"], 1.0e7; rtol=1.0e-15)
        @test summary["c_power"] == 5
        refs = SMV.expected_metric_references()
        @test Set(keys(refs)) == Set(SMV.EXPECTED_METRIC_IDS)
        @test isapprox(refs["centre_fit_ratio_nx65"], 2 / 3; rtol=1.0e-15)
        @test isapprox(refs["profile_ratio_nx129"], sqrt(8 / 11); rtol=1.0e-15)
        @test isapprox(refs["s33_order_b_65_129"], 2.0; rtol=1.0e-15)
        for m in raw_record["metrics"]
            @test isapprox(
                Float64(m["reference"]), refs[String(m["id"])]; rtol=1.0e-12, atol=1.0e-12
            )
        end
    end

    @testset "Reproducibility check, CLI table and export round trip" begin
        raw_record = TOML.parsefile(record_path)
        ok_repro, errs_repro = SMV.check_validation_reproducibility(bench, record_path)
        @test ok_repro
        @test isempty(errs_repro)
        @test run_momentum_validation_check(record_path; bench=bench)
        @test run_validation_cli(["--check", record_path]; momentum_bench=bench) == 0
        @test !run_momentum_validation_check(record_path * ".missing"; bench=bench)
        @test_throws DomainError run_momentum_validation_check(""; bench=bench)
        @test validation_case("spherical_momentum").bench_keyword == :momentum_bench

        mktempdir() do tmpdir
            empty_bench = (; runs=Dict{Tuple,Any}(), metrics=Dict{String,Float64}())
            no_arg = capture_momentum_stdio(
                () -> run_validation_cli(
                    ["--check"];
                    bench=empty_bench,
                    darcy_bench=empty_bench,
                    compaction_bench=empty_bench,
                    diffusion_bench=empty_bench,
                    momentum_bench=bench,
                ),
            )
            @test no_arg.value == 1
            @test occursin("spherical_momentum.toml", no_arg.out)
            @test occursin("[reproducibility]", no_arg.err)
            @test last(VALIDATION_CASES).id == "spherical_momentum"

            exported = joinpath(tmpdir, "exported.toml")
            commit = String(raw_record["metadata"]["src_commit"])
            @test SMVE.run_export_cli(
                [exported]; src_commit=commit, bench=bench, check_clean_tree=false
            ) == 0
            raw_exported = TOML.parsefile(exported)
            ok_valid, errs_valid = SMV.validate_record_dict(
                raw_exported; require_doc_pages=false, check_doc_sync=false
            )
            @test ok_valid
            @test isempty(errs_valid)
            @test run_momentum_validation_check(exported; bench=bench)
            @test [String(m["id"]) for m in raw_exported["metrics"]] == collect(SMV.EXPECTED_METRIC_IDS)

            # A moved value is refused without the flag and written with it.
            target = "l2_a_on_nx65"
            moved_metrics = copy(bench.metrics)
            moved_metrics[target] += 2 * SMV.PHYSICS_REGRESSION_BANDS[target].atol
            moved_bench = (; runs=bench.runs, metrics=moved_metrics)
            @test SMVE.run_export_cli(
                [exported]; src_commit=commit, bench=moved_bench, check_clean_tree=false
            ) == 1
            @test SMVE.run_export_cli(
                ["--accept-physics-change", exported];
                src_commit=commit,
                bench=moved_bench,
                check_clean_tree=false,
            ) == 0
            @test SMVE.run_export_cli(
                ["--bad-flag"]; src_commit=commit, bench=bench, check_clean_tree=false
            ) == 2
            @test SMVE.run_export_cli(
                [exported];
                src_commit=commit,
                bench=bench,
                check_clean_tree=true,
                dirty_status=" M src/config.jl",
            ) == 1
            @test_throws DomainError SMVE.run_export_cli(
                [exported]; src_commit="abc", bench=bench, check_clean_tree=false
            )
            # A value that fails its own tolerance is never written.
            failing = copy(bench.metrics)
            failing["elim_a_nx33"] = 1.0
            @test SMVE.run_export_cli(
                ["--accept-physics-change", exported];
                src_commit=commit,
                bench=(; runs=bench.runs, metrics=failing),
                check_clean_tree=false,
            ) == 1
            cases = (
                (0.01, 0.0, "le_atol"), (0.01, 1.0, "le_atol"), (1.0, 1.001, "isapprox")
            )
            @test isequal(
                [SMVE.metric_status(v, ref, 0.05, 0.0, c) for (v, ref, c) in cases],
                ["PASS", "FAIL", "PASS"],
            )
            @test isequal(SMVE.metric_status(1.0, 1.0, 0.0, 0.0, "other"), "FAIL")
            @test occursin(r"^[0-9a-f]{40}$", SMVE.last_src_commit())
            # Outside a git repository git log fails.
            no_repo = mkpath(joinpath(tmpdir, "no_repo"))
            withenv("GIT_DIR" => joinpath(no_repo, "absent.git")) do
                @test_throws ProcessFailedException redirect_stdio(
                    () -> SMVE.last_src_commit(no_repo); stderr=devnull
                )
            end
            @test first(SMVE.check_clean_source_tree()) isa Bool

            # A changed configuration value and hash are rejected.
            changed_cfg = deepcopy(raw_record)
            changed_cfg["config"]["eta_Pa_s"] *= 2
            changed_cfg["metadata"]["config_sha256"] = repeat("0", 64)
            cfg_path = write_momentum_toml(joinpath(tmpdir, "changed.toml"), changed_cfg)
            ok_cfg, errs_cfg = SMV.check_validation_reproducibility(bench, cfg_path)
            @test !ok_cfg
            @test any(occursin("config.eta_Pa_s", e) for e in errs_cfg)
            @test any(occursin("config_sha256", e) for e in errs_cfg)
            no_key = deepcopy(raw_record)
            delete!(no_key["config"], "F0_N_m3")
            @test SMV._config_mismatches(
                no_key, SMV.make_config(65; spherical_flow=true)
            ) == ["Recorded config missing key 'F0_N_m3'"]
            extra = deepcopy(raw_record)
            extra["config"]["extra_key"] = 1.0
            @test SMV._config_mismatches(extra, SMV.make_config(65; spherical_flow=true)) ==
                ["Recorded config has unexpected key 'extra_key'"]
            odd = deepcopy(raw_record)
            iratio = findfirst(m -> m["id"] == "centre_fit_ratio_nx65", odd["metrics"])
            odd["metrics"][iratio]["reference"] = 0.5
            odd["metrics"][iratio]["value"] = 0.5
            odd["metrics"][end]["id"] = "not_a_metric"
            odd_path = write_momentum_toml(joinpath(tmpdir, "odd.toml"), odd)
            ok_odd, errs_odd = SMV.check_validation_reproducibility(bench, odd_path)
            @test !ok_odd
            @test any(occursin("Recorded metric IDs", e) for e in errs_odd)
            @test any(occursin("does not match recomputed reference", e) for e in errs_odd)
            @test any(occursin("Benchmark output missing metric", e) for e in errs_odd)
            failing_bench = (; runs=bench.runs, metrics=failing)
            this_julia = deepcopy(raw_record)
            this_julia["metadata"]["julia_version"] = string(VERSION)
            this_julia_path = write_momentum_toml(
                joinpath(tmpdir, "this_julia.toml"), this_julia
            )
            ok_fail, errs_fail = SMV.check_validation_reproducibility(
                failing_bench, this_julia_path
            )
            @test !ok_fail
            @test any(occursin("fails tolerance", e) for e in errs_fail)
            @test any(occursin("differs from recorded value", e) for e in errs_fail)
            other_julia = deepcopy(raw_record)
            other_julia["metadata"]["julia_version"] = "1.0.0"
            other_path = write_momentum_toml(joinpath(tmpdir, "julia.toml"), other_julia)
            @test first(SMV.check_validation_reproducibility(bench, other_path))
            @test first(
                SMV.check_validation_reproducibility(
                    bench, record_path; check_recorded_values=false
                ),
            )
            # A recorded value moved beyond REPLAY_RTOL but inside every tolerance fails
            # only the replay, and only for the Julia minor version of the record.
            tampered = deepcopy(raw_exported)
            tmid = findfirst(m -> m["id"] == "l2_a_on_nx65", tampered["metrics"])
            tampered["metrics"][tmid]["value"] += 1.0e-5
            tampered_path = write_momentum_toml(joinpath(tmpdir, "tampered.toml"), tampered)
            ok_rep, errs_rep = SMV.check_validation_reproducibility(bench, tampered_path)
            @test !ok_rep
            @test occursin("differs from recorded value", only(errs_rep))
            @test first(
                SMV.check_validation_reproducibility(
                    bench, tampered_path; check_recorded_values=false
                ),
            )
            tampered["metadata"]["julia_version"] = "1.0.0"
            tampered_other = write_momentum_toml(
                joinpath(tmpdir, "tampered_other.toml"), tampered
            )
            @test first(SMV.check_validation_reproducibility(bench, tampered_other))
        end
    end

    @testset "checks stop on an interrupt and report other errors with a stack trace" begin
        commit = String(raw_record_all["metadata"]["src_commit"])
        @test_throws InterruptException run_momentum_validation_check(
            record_path; bench=MomentumInterruptingBench()
        )
        crash = capture_momentum_stdio(
            () -> run_momentum_validation_check(record_path; bench=MomentumCrashingBench())
        )
        @test crash.value === false
        @test occursin("solver failure inside the benchmark", crash.err)
        @test occursin("Stacktrace", crash.err)
        mktempdir() do tmpdir
            target = joinpath(tmpdir, "record.toml")
            @test_throws InterruptException SMVE.run_export_cli(
                [target];
                src_commit=commit,
                bench=MomentumInterruptingBench(),
                check_clean_tree=false,
            )
            failed = capture_momentum_stdio(
                () -> SMVE.run_export_cli(
                    [target];
                    src_commit=commit,
                    bench=MomentumCrashingBench(),
                    check_clean_tree=false,
                ),
            )
            @test failed.value == 1
            @test occursin("Stacktrace", failed.err)
            @test !isfile(target)
            not_repo = mkpath(joinpath(tmpdir, "not_a_repo"))
            withenv("GIT_DIR" => joinpath(not_repo, "absent.git")) do
                for cli in (SMVE.run_export_cli, SMVE.export_main)
                    no_git = capture_momentum_stdio(
                        () -> cli(
                            [target];
                            src_commit=commit,
                            bench=bench,
                            check_clean_tree=true,
                            repo_root=not_repo,
                        ),
                    )
                    @test no_git.value == 1
                    @test occursin("status --porcelain", no_git.err)
                    @test !isfile(target)
                end
            end
        end
    end

    @testset "export command line: one record path, unknown flags and an empty path" begin
        commit = String(raw_record_all["metadata"]["src_commit"])
        mktempdir() do tmpdir
            cd(tmpdir) do
                code = SMVE.run_export_cli(
                    ["-h"]; src_commit=commit, bench=bench, check_clean_tree=false
                )
                @test code == 2
                @test !isfile(joinpath(tmpdir, "-h"))
            end
            p1, p2 = joinpath(tmpdir, "first.toml"), joinpath(tmpdir, "second.toml")
            @test SMVE.run_export_cli(
                [p1, p2]; src_commit=commit, bench=bench, check_clean_tree=false
            ) == 2
            @test !isfile(p1) && !isfile(p2)
            empty_run = capture_momentum_stdio(
                () -> SMVE.export_main(
                    [""]; src_commit=commit, bench=bench, check_clean_tree=false
                ),
            )
            @test empty_run.value == 1
            @test occursin("Record path must be non-empty", empty_run.err)
            @test SMVE.export_main(
                [p2]; src_commit=commit, bench=bench, check_clean_tree=false
            ) == 0
            @test isfile(p2)
            @test SMVE.momentum_validation_module() === SMV
            @test_throws ErrorException SMVE.momentum_validation_module((Module(:Empty),))
        end
    end

    @testset "malformed records are reported, not thrown" begin
        mktempdir() do tmpdir
            broken = joinpath(tmpdir, "broken.toml")
            write(broken, "[metadata\nid = ")
            parsed = capture_momentum_stdio(
                () -> run_momentum_validation_check(broken; bench=bench)
            )
            @test parsed.value === false
            @test occursin("[schema] Cannot parse", parsed.err)
            for check in (
                () -> SMV.check_validation_reproducibility(bench, broken),
                () -> SMV.check_physics_regression_bands(bench.metrics, broken),
            )
                ok_parse, errs_parse = check()
                @test !ok_parse
                @test occursin("Cannot parse", only(errs_parse))
            end
            @test_throws DomainError SMV.check_validation_reproducibility(bench, "")
            @test !first(SMV.check_validation_reproducibility(bench, broken * ".x"))
            schema_bad = deepcopy(raw_record_all)
            delete!(schema_bad, "config")
            schema_path = write_momentum_toml(joinpath(tmpdir, "schema.toml"), schema_bad)
            ok_schema, errs_schema = SMV.check_validation_reproducibility(
                bench, schema_path
            )
            @test !ok_schema
            @test any(occursin("config", e) for e in errs_schema)
        end
    end

    @testset "figure data export: CSV columns and the exported files" begin
        script_host = Module(:MomentumFigureData)
        Core.eval(script_host, :(include(path) = Base.include($script_host, path)))
        Base.include(
            script_host,
            joinpath(@__DIR__, "..", "scripts", "export_validation_momentum_data.jl"),
        )
        write_columns = getfield(script_host, :write_columns)
        export_all = getfield(script_host, :export_all)
        export_data_main = getfield(script_host, :export_data_main)
        mktempdir() do tmpdir
            path = joinpath(tmpdir, "columns.csv")
            write(path, "kept\n")
            for (header, columns) in (
                (["a"], Vector{Float64}[]),
                (["a", "b", "c"], [[1.0, 2.0], [3.0, 4.0]]),
                (["a", "b"], [[1.0, 2.0], [3.0, 4.0, 5.0]]),
            )
                @test_throws ArgumentError Base.invokelatest(
                    write_columns, path, header, columns
                )
            end
            @test read(path, String) == "kept\n"
            outdir = joinpath(tmpdir, "data")
            @test Base.invokelatest(export_all, outdir; runs=runs) == outdir
            @test sort(readdir(outdir)) == sort([
                "convergence.csv",
                "convergence_c.csv",
                "profile_c_s33_65.csv",
                "profile_c_u_65.csv",
                "profile_off_65.csv",
                "profile_on_65.csv",
                "reference_curves.csv",
                "reference_curves_c.csv",
            ])
            conv_c = readlines(joinpath(outdir, "convergence_c.csv"))
            @test conv_c[1] == "nx,l2_c,l2_inside_c,s33_c,pt_c,pf_c,eii_c,sii_c,hs_c"
            @test length(conv_c) == 1 + length(SMV.NX_LIST)
            row65_c = parse.(Float64, split(conv_c[3], ","))
            @test isapprox(row65_c[3], recorded["l2_inside_c_nx65"]; rtol=1.0e-15)
            @test isapprox(row65_c[9], recorded["hs_c_nx65"]; rtol=1.0e-15)
            curves_c = readlines(joinpath(outdir, "reference_curves_c.csv"))
            @test curves_c[1] == "rho,u_r,tau33"
            ref_c = [parse.(Float64, split(l, ",")) for l in curves_c[2:end]]
            @test all(isapprox(v, 0.0; atol=1.0e-15) for v in ref_c[1])
            @test all(isapprox(v, 0.0; atol=1.0e-12) for v in ref_c[end][2:end])
            # u_r = -10 rho (1 - rho^2)^4 is most negative at rho = 1/3, -2.6.
            @test isapprox(minimum(v[2] for v in ref_c), -10 / 3 * (8 / 9)^4; rtol=1.0e-4)
            for (file, col) in (("profile_c_u_65.csv", 2), ("profile_c_s33_65.csv", 3))
                lines = readlines(joinpath(outdir, file))
                @test lines[1] == (col == 2 ? "rho,u_r" : "rho,tau33")
                vals_c = [parse.(Float64, split(l, ",")) for l in lines[2:end]]
                # The row and the column of nodes half a cell off the axes, 4 half-axes
                # with r < R: 22 vx or vy nodes and 23 P nodes per half-axis at Nx = 65.
                @test length(vals_c) == 4 * (col == 2 ? 22 : 23)
                peak = maximum(abs(v[col]) for v in ref_c)
                exact_c(rho) =
                    if col == 2
                        SMV.exact_velocity(:C, rho, 0.0)[1]
                    else
                        SMV.exact_tau33(:C, rho, 0.0)
                    end
                @test all(abs(v[2] - exact_c(v[1])) < 0.05 * peak for v in vals_c)
                @test all(0.0 < v[1] < 1.0 for v in vals_c)
            end
            conv = readlines(joinpath(outdir, "convergence.csv"))
            @test conv[1] ==
                "nx,l2_a_on,l2_a_off,l2_b,s33_b,centre_fit_a_on,centre_node_a_on"
            @test length(conv) == 1 + length(SMV.NX_LIST)
            row65 = parse.(Float64, split(conv[3], ","))
            @test isapprox(row65[2], recorded["l2_a_on_nx65"]; rtol=1.0e-15)
            @test isapprox(row65[7], recorded["centre_node_a_on_nx65"]; rtol=1.0e-15)
            curves = readlines(joinpath(outdir, "reference_curves.csv"))
            @test curves[1] == "rho,u_on,u_off"
            first_row = parse.(Float64, split(curves[2], ","))
            last_row = parse.(Float64, split(curves[end], ","))
            @test all(isapprox(v, 0.0; atol=1.0e-15) for v in first_row)
            @test isapprox(last_row[1], 1.0; rtol=1.0e-15)
            @test all(isapprox(v, 0.0; atol=1.0e-15) for v in last_row[2:end])
            prof = readlines(joinpath(outdir, "profile_on_65.csv"))
            @test prof[1] == "rho,u"
            vals = [parse.(Float64, split(l, ",")) for l in prof[2:end]]
            # 4 half-axes with the nodes inside R: 4 x 23 at Nx = 65.
            @test length(vals) == 4 * 23
            umax = maximum(abs(v[2]) for v in vals)
            exact(rho) = rho * SMV.omega_hat(:A, rho; spherical_flow=true)
            @test all(abs(v[2] - exact(v[1])) < 0.05 * umax for v in vals)
            @test 1.0e-3 < umax < 2.0e-3
            cd(tmpdir) do
                for bad in (String[], ["a", "b"], [""], ["  "], ["-h"], ["--out"])
                    @test_throws ArgumentError Base.invokelatest(
                        export_data_main, bad; nx=SMV.NX_MIN
                    )
                end
            end
            @test !ispath(joinpath(tmpdir, "-h"))
            # Without runs the script runs the case itself, here on the coarsest grid.
            small = joinpath(tmpdir, "small")
            @test Base.invokelatest(export_data_main, [small]; nx=SMV.NX_MIN) == small
            off_small = readlines(joinpath(small, "profile_off_$(SMV.NX_MIN).csv"))
            # 4 half-axes with the nodes inside R: 4 x 10 at Nx = 29.
            @test SMV.NX_MIN == 29
            @test length(off_small) - 1 == 4 * 10
            u_off = [parse.(Float64, split(l, ","))[2] for l in off_small[2:end]]
            @test 0.0 < maximum(u_off) < 3.0e-3
        end
    end
end

"""
    SphericalMomentumValidation

Steady Stokes flow of a uniform viscous matrix driven by a body force inside a sphere,
against exact solutions of the momentum equation of the spherical slice.

With `geometry.spherical_flow = true` the momentum rows solve
`(1 / w) d_j (w tau_ij) - tau_33 d_i w / w - d_i P + f_i = 0` with `w = 2 r` and the 3-D
deviatoric stress `tau_ij`, and with `false` the 2-D rows `d_j tau2D_ij - d_i P + f_i = 0`
with the 2-D deviatoric stress `tau2D_xx = -tau2D_yy = S_xx = eta (eps_xx - eps_yy)` and
`tau2D_xy = 2 eta eps_xy`. The case sets the body force `f = rho g` in `gx` and `gy` and
compares the solved velocity with 3 exact flows that vanish for `r >= R_PLANET`:

- case `:A`, azimuthal shear `u = r Omega(r) e_theta` with a closed form for both settings;
- case `:A0`, azimuthal shear with `Omega(0) = 0`, for the option on;
- case `:B`, a manufactured flow with `u_r != 0` and `tau_33 = O(1)`, for the option on;
- case `:C`, a radial compaction flow `u = U0 R grad (1 - rho^2)^5` of a compressible
  matrix with `eta / zeta = 1` and the Darcy flux `q = -u`, for the option on. For a radial
  flow the slice equations are those of a 3-D sphere, with the volumetric trace
  `-2 eta div3 / 3` in `tau_ij` and `tau_33` of order 1.

Cases `:A`, `:A0` and `:B` have `eta / zeta = ETA_RATIO`; case `:C` has `ETA_RATIO_C` and
the ambient pressure `P_AMB_C`. The case builds the state with `init_simulation`, sets the grid fields directly and calls
`solve_stokes_darcy!` once on timestep 1.
"""
module SphericalMomentumValidation

using Erebus
using Erebus.Config
using Logging
using SHA
using StaticArrays
using TOML

include(joinpath(@__DIR__, "..", "tools", "validation_records.jl"))

const R_PLANET = 50_000.0
const L_BOX = 140_000.0
const X_CENTRE = 70_000.0
const Y_CENTRE = 70_000.0
const ETA_ROCK = 1.0e20
const ETA_RATIO = 1.0e-6
const ETA_MIN = 1.0e12
const ETA_MAX = 1.0e30
const PHI = 0.2
const K_PERM = 1.0e-13
const MU_F = 1.0e-3
const RHO_SOLID = 3000.0
const RHO_FLUID = 0.0
const RHOCP = 3.0e6
const P_AMB = 100.0
const T_SURF = 200.0
const SHEAR_MODULUS = 1.0e35
const STRENGTH_FILL = 1.0e20
const DT_YR = 1.0e3
const F0 = 100.0
const U0 = 1.0e-9
const B_POWER = 6
const C_POWER = 5
const ETA_RATIO_C = 1.0
const P_AMB_C = 1.0e7
const FIT_RHO_MIN = 0.18
const FIT_RHO_MAX = 0.85
const FIT_DEGREE = 5
const NX_LIST = (33, 65, 129)
const NX_REFERENCE = 65
const NX_MIN = 29
const CASES = (:A, :A0, :B, :C)

const QB = Rational{BigInt}

"""
    _spoly_mul(a, b)

Return the product of 2 polynomials with exact rational coefficients; entry `k + 1` holds
the coefficient of the power `k`.

Parameters
----------
a, b : Vector{Rational{BigInt}}
    Coefficients.

Returns
-------
Vector{Rational{BigInt}}
    Coefficients of `a b`.
"""
function _spoly_mul(a::Vector{QB}, b::Vector{QB})
    c = zeros(QB, length(a) + length(b) - 1)
    for i in eachindex(a), j in eachindex(b)
        c[i + j - 1] += a[i] * b[j]
    end
    return c
end

"""
    _spoly_pow(a, n)

Return the `n`-th power of a polynomial (see `_spoly_mul`).

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
function _spoly_pow(a::Vector{QB}, n::Int)
    c = QB[1]
    for _ in 1:n
        c = _spoly_mul(c, a)
    end
    return c
end

# Closed forms in s = rho^2: Omega in units of F0 R / eta, and the azimuthal force divided
# by rho, f_theta / (F0 rho). The test checks them against the radial momentum equation.
const OMEGA_A_ON = _spoly_mul(_spoly_pow(QB[1, -1], 4), QB[1, -6]) .// 100
const OMEGA_A_OFF = _spoly_mul(_spoly_pow(QB[1, -1], 4), QB[3, -13]) .// 200
const FORCE_A = _spoly_mul(_spoly_pow(QB[1, -1], 2), QB[1, -32 // 5, 39 // 5])
const OMEGA_A0 = _spoly_mul(_spoly_pow(QB[1, -1], 4), QB[0, 1]) .// -10
const FORCE_A0 = _spoly_mul(_spoly_pow(QB[1, -1], 2), QB[1, -46 // 5, 13])

"""
    _check_case(case; spherical_flow=true)

Check that `case` is a case of the module and that it has an exact solution for the
setting `spherical_flow`.

Parameters
----------
case : Symbol
    `:A`, `:A0`, `:B` or `:C`.
spherical_flow : Bool, optional
    Value of `geometry.spherical_flow` (default `true`).

Returns
-------
Nothing

Raises
------
ArgumentError
    If `case` is not in `CASES`, or `case` is `:A0`, `:B` or `:C` and `spherical_flow` is
    `false` (only case `:A` has an exact solution without the option).
"""
function _check_case(case::Symbol; spherical_flow::Bool=true)
    case in CASES || throw(ArgumentError("case must be one of $CASES, got :$case"))
    (spherical_flow || case == :A) || throw(
        ArgumentError("case :$case has an exact solution only with spherical_flow = true"),
    )
    return nothing
end

"""
    eta_zeta_ratio(case)

Return the ratio `eta / zeta` of the shear to the compaction viscosity of `case`:
`ETA_RATIO_C` for case `:C` and `ETA_RATIO` for the others.

Parameters
----------
case : Symbol
    `:A`, `:A0`, `:B` or `:C`.

Returns
-------
Float64

Raises
------
ArgumentError
    If `case` is not in `CASES`.
"""
function eta_zeta_ratio(case::Symbol)
    _check_case(case)
    return case == :C ? ETA_RATIO_C : ETA_RATIO
end

"""
    ambient_pressure(case)

Return the ambient pressure `p_s` [Pa] of `case`, held on the anchor rings of `P_t` and
`P_f`: `P_AMB_C` for case `:C`, where it keeps `P_t > 0`, and `P_AMB` for the others.

Parameters
----------
case : Symbol
    `:A`, `:A0`, `:B` or `:C`.

Returns
-------
Float64
    Pressure [Pa].

Raises
------
ArgumentError
    If `case` is not in `CASES`.
"""
function ambient_pressure(case::Symbol)
    _check_case(case)
    return case == :C ? P_AMB_C : P_AMB
end

"""
    _check_point(X, Y)

Check a nondimensional position `(X, Y) = (x - xc, y - yc) / R_PLANET`.

Parameters
----------
X, Y : Real
    Nondimensional coordinates relative to the centre.

Returns
-------
Nothing

Raises
------
DomainError
    If `X` or `Y` is not finite.
"""
function _check_point(X::Real, Y::Real)
    isfinite(X) || throw(DomainError(X, "X must be finite"))
    isfinite(Y) || throw(DomainError(Y, "Y must be finite"))
    return nothing
end

"""
    _check_rho(rho)

Check a nondimensional radius `rho = r / R_PLANET`.

Parameters
----------
rho : Real
    Nondimensional radius.

Returns
-------
Nothing

Raises
------
DomainError
    If `rho` is not finite or negative.
"""
function _check_rho(rho::Real)
    (isfinite(rho) && rho >= 0) ||
        throw(DomainError(rho, "rho must be finite and non-negative"))
    return nothing
end

"""
    omega_coefficients(case; spherical_flow)

Return the exact angular velocity `Omega` of an azimuthal case in units of `F0 R / eta`, as
coefficients in `s = (r / R)^2` (entry `k + 1` holds `s^k`). Case `:A` with the option:
`(1 - s)^4 (1 - 6 s) / 100`; without: `(1 - s)^4 (3 - 13 s) / 200`; case `:A0`:
`-s (1 - s)^4 / 10`.

Parameters
----------
case : Symbol
    `:A` or `:A0`.
spherical_flow : Bool
    Value of `geometry.spherical_flow`.

Returns
-------
Vector{Rational{BigInt}}
    Coefficients in `s`.

Raises
------
ArgumentError
    If `case` is not `:A` or `:A0`, or `case = :A0` with `spherical_flow = false`.
"""
function omega_coefficients(case::Symbol; spherical_flow::Bool)
    _check_case(case; spherical_flow=spherical_flow)
    case in (:B, :C) && throw(ArgumentError("case :$case is not azimuthal"))
    case == :A0 && return copy(OMEGA_A0)
    return copy(spherical_flow ? OMEGA_A_ON : OMEGA_A_OFF)
end

"""
    force_coefficients(case)

Return the azimuthal force of an azimuthal case divided by `F0 r / R`, as coefficients in
`s = (r / R)^2`. Case `:A`: `(1 - s)^2 (1 - 32 s / 5 + 39 s^2 / 5)`, the same force for
both settings; case `:A0`: `(1 - s)^2 (1 - 46 s / 5 + 13 s^2)`.

Parameters
----------
case : Symbol
    `:A` or `:A0`.

Returns
-------
Vector{Rational{BigInt}}
    Coefficients in `s`.

Raises
------
ArgumentError
    If `case` is not `:A` or `:A0`.
"""
function force_coefficients(case::Symbol)
    _check_case(case)
    case in (:B, :C) && throw(ArgumentError("case :$case is not azimuthal"))
    return copy(case == :A ? FORCE_A : FORCE_A0)
end

"""
    _seval(c, s)

Return the polynomial with rational coefficients `c` at `s` in Float64.

Parameters
----------
c : Vector{Rational{BigInt}}
    Coefficients.
s : Float64
    Argument.

Returns
-------
Float64
"""
function _seval(c::Vector{QB}, s::Float64)
    v = 0.0
    for k in reverse(eachindex(c))
        v = v * s + Float64(c[k])
    end
    return v
end

"""
    omega_hat(case, rho; spherical_flow)

Return the exact angular velocity `Omega` of an azimuthal case in units of `F0 R / eta` at
`rho = r / R_PLANET`, 0 for `rho >= 1` (see `omega_coefficients`).

Parameters
----------
case : Symbol
    `:A` or `:A0`.
rho : Real
    Nondimensional radius, finite and non-negative.
spherical_flow : Bool
    Value of `geometry.spherical_flow`.

Returns
-------
Float64

Raises
------
ArgumentError
    From `omega_coefficients`.
DomainError
    If `rho` is not finite or negative.
"""
function omega_hat(case::Symbol, rho::Real; spherical_flow::Bool)
    c = omega_coefficients(case; spherical_flow=spherical_flow)
    _check_rho(rho)
    rho >= 1 && return 0.0
    return _seval(c, Float64(rho)^2)
end

"""
    force_hat(case, rho)

Return the azimuthal force of an azimuthal case in units of `F0` at `rho = r / R_PLANET`,
0 for `rho >= 1` (see `force_coefficients`).

Parameters
----------
case : Symbol
    `:A` or `:A0`.
rho : Real
    Nondimensional radius, finite and non-negative.

Returns
-------
Float64

Raises
------
ArgumentError
    From `force_coefficients`.
DomainError
    If `rho` is not finite or negative.
"""
function force_hat(case::Symbol, rho::Real)
    c = force_coefficients(case)
    _check_rho(rho)
    rho >= 1 && return 0.0
    return Float64(rho) * _seval(c, Float64(rho)^2)
end

# Polynomials in (X, Y) with exact rational coefficients: (power of X, power of Y) => value.
const XYPoly = Dict{Tuple{Int,Int},QB}

"""
    _xy_add(a, b)

Return the sum of 2 polynomials in `(X, Y)` without zero terms.

Parameters
----------
a, b : Dict{Tuple{Int,Int},Rational{BigInt}}
    Terms by the powers of `X` and `Y`.

Returns
-------
Dict{Tuple{Int,Int},Rational{BigInt}}
"""
function _xy_add(a::XYPoly, b::XYPoly)
    c = copy(a)
    for (k, v) in b
        c[k] = get(c, k, QB(0)) + v
    end
    return filter!(kv -> !iszero(kv.second), c)
end

"""
    _xy_scale(a, s)

Return the polynomial `a` in `(X, Y)` times the number `s`.

Parameters
----------
a : Dict{Tuple{Int,Int},Rational{BigInt}}
    Terms.
s : Real
    Exact factor (integer or rational).

Returns
-------
Dict{Tuple{Int,Int},Rational{BigInt}}
"""
function _xy_scale(a::XYPoly, s)
    return filter!(kv -> !iszero(kv.second), XYPoly(k => QB(s) * v for (k, v) in a))
end

"""
    _xy_mul(a, b)

Return the product of 2 polynomials in `(X, Y)`.

Parameters
----------
a, b : Dict{Tuple{Int,Int},Rational{BigInt}}
    Terms.

Returns
-------
Dict{Tuple{Int,Int},Rational{BigInt}}
"""
function _xy_mul(a::XYPoly, b::XYPoly)
    c = XYPoly()
    for (ka, va) in a, (kb, vb) in b
        k = (ka[1] + kb[1], ka[2] + kb[2])
        c[k] = get(c, k, QB(0)) + va * vb
    end
    return filter!(kv -> !iszero(kv.second), c)
end

"""
    _xy_dx(a)

Return the derivative of a polynomial in `(X, Y)` with respect to `X`.

Parameters
----------
a : Dict{Tuple{Int,Int},Rational{BigInt}}
    Terms.

Returns
-------
Dict{Tuple{Int,Int},Rational{BigInt}}
"""
_xy_dx(a::XYPoly) = XYPoly((k[1] - 1, k[2]) => k[1] * v for (k, v) in a if k[1] > 0)

"""
    _xy_dy(a)

Return the derivative of a polynomial in `(X, Y)` with respect to `Y`.

Parameters
----------
a : Dict{Tuple{Int,Int},Rational{BigInt}}
    Terms.

Returns
-------
Dict{Tuple{Int,Int},Rational{BigInt}}
"""
_xy_dy(a::XYPoly) = XYPoly((k[1], k[2] - 1) => k[2] * v for (k, v) in a if k[2] > 0)

"""
    _xy_divrem_rho2(a)

Return the quotient and the remainder of a polynomial in `(X, Y)` divided by `X^2 + Y^2`.
Each term `X^p Y^q` with `p >= 2` is rewritten as
`X^(p - 2) Y^q (X^2 + Y^2) - X^(p - 2) Y^(q + 2)` until every remaining term has `p <= 1`;
the remainder is 0 exactly when `X^2 + Y^2` divides the polynomial.

Parameters
----------
a : Dict{Tuple{Int,Int},Rational{BigInt}}
    Terms.

Returns
-------
Tuple
    `(quotient, remainder)`, each a `Dict{Tuple{Int,Int},Rational{BigInt}}`; every term of
    the remainder has the power of `X` at most 1.
"""
function _xy_divrem_rho2(a::XYPoly)
    rem = copy(a)
    quot = XYPoly()
    while true
        high = [k for (k, v) in rem if k[1] >= 2 && !iszero(v)]
        isempty(high) && break
        k = maximum(high)
        v = rem[k]
        delete!(rem, k)
        q = (k[1] - 2, k[2])
        quot[q] = get(quot, q, QB(0)) + v
        lower = (k[1] - 2, k[2] + 2)
        rem[lower] = get(rem, lower, QB(0)) - v
    end
    return (filter!(kv -> !iszero(kv.second), quot), filter!(kv -> !iszero(kv.second), rem))
end

"""
    _xy_div_rho2(a)

Return the exact quotient of a polynomial in `(X, Y)` by `X^2 + Y^2` (see
`_xy_divrem_rho2`).

Parameters
----------
a : Dict{Tuple{Int,Int},Rational{BigInt}}
    Terms.

Returns
-------
Dict{Tuple{Int,Int},Rational{BigInt}}
    The quotient.

Raises
------
ArgumentError
    If the division leaves a remainder.
"""
function _xy_div_rho2(a::XYPoly)
    quot, rem = _xy_divrem_rho2(a)
    isempty(rem) || throw(ArgumentError("the polynomial is not divisible by X^2 + Y^2"))
    return quot
end

"""
    _xy_terms(a)

Return the terms of a polynomial in `(X, Y)` as Float64 for fast evaluation.

Parameters
----------
a : Dict{Tuple{Int,Int},Rational{BigInt}}
    Terms.

Returns
-------
Vector{Tuple{Int,Int,Float64}}
    `(power of X, power of Y, coefficient)`, sorted by the powers.
"""
_xy_terms(a::XYPoly) = sort!([(k[1], k[2], Float64(v)) for (k, v) in a])

"""
    _xy_eval(terms, X, Y)

Return the value of a polynomial in `(X, Y)` from its terms (see `_xy_terms`).

Parameters
----------
terms : Vector{Tuple{Int,Int,Float64}}
    Terms.
X, Y : Float64
    Arguments.

Returns
-------
Float64
"""
function _xy_eval(terms::Vector{Tuple{Int,Int,Float64}}, X::Float64, Y::Float64)
    v = 0.0
    for (p, q, c) in terms
        v += c * X^p * Y^q
    end
    return v
end

"""
    _xy_eval_exact(a, X, Y)

Return the exact value of a polynomial in `(X, Y)` at rational arguments.

Parameters
----------
a : Dict{Tuple{Int,Int},Rational{BigInt}}
    Terms.
X, Y : Rational{BigInt}
    Arguments.

Returns
-------
Rational{BigInt}
"""
function _xy_eval_exact(a::XYPoly, X::QB, Y::QB)
    return sum((v * X^k[1] * Y^k[2] for (k, v) in a); init=QB(0))
end

const XY_X = XYPoly((1, 0) => QB(1))
const XY_Y = XYPoly((0, 1) => QB(1))
const XY_RHO2 = XYPoly((2, 0) => QB(1), (0, 2) => QB(1))

"""
    _xy_from_s(c)

Return the polynomial `sum_k c[k + 1] (X^2 + Y^2)^k` in `(X, Y)`.

Parameters
----------
c : Vector{Rational{BigInt}}
    Coefficients in `s = X^2 + Y^2`.

Returns
-------
Dict{Tuple{Int,Int},Rational{BigInt}}
"""
function _xy_from_s(c::Vector{QB})
    p = XYPoly()
    pow = XYPoly((0, 0) => QB(1))
    for v in c
        p = _xy_add(p, _xy_scale(pow, v))
        pow = _xy_mul(pow, XY_RHO2)
    end
    return p
end

"""
    _stress_divergence(ux, uy)

Return the stress divergence of the documented momentum equation with the option on, for
an incompressible velocity with unit viscosity and uniform pressure:
`D_i = d_j tau_ij + (tau_ij X_j - tau_33 X_i) / rho^2` with `tau_ij = d_i u_j + d_j u_i`
and `tau_33 = 2 (X ux + Y uy) / rho^2`. `tau_33` is an exact polynomial; each `D_i` is a
polynomial plus a remainder divided by `rho^2` (`_xy_divrem_rho2`).

Parameters
----------
ux, uy : Dict{Tuple{Int,Int},Rational{BigInt}}
    Velocity components.

Returns
-------
NamedTuple
    `dx`, `dy`: the polynomial parts of the divergence; `dx_rem`, `dy_rem`: the remainders
    whose quotient by `rho^2` completes it; `tau33`: `tau_33`.

Raises
------
ArgumentError
    If `rho^2` does not divide `X ux + Y uy`.
"""
function _stress_divergence(ux::XYPoly, uy::XYPoly)
    txx = _xy_scale(_xy_dx(ux), 2)
    tyy = _xy_scale(_xy_dy(uy), 2)
    txy = _xy_add(_xy_dy(ux), _xy_dx(uy))
    radial = _xy_add(_xy_mul(ux, XY_X), _xy_mul(uy, XY_Y))
    t33 = _xy_scale(_xy_div_rho2(radial), 2)
    hx = _xy_add(
        _xy_add(_xy_mul(txx, XY_X), _xy_mul(txy, XY_Y)), _xy_scale(_xy_mul(t33, XY_X), -1)
    )
    hy = _xy_add(
        _xy_add(_xy_mul(txy, XY_X), _xy_mul(tyy, XY_Y)), _xy_scale(_xy_mul(t33, XY_Y), -1)
    )
    qx, rx = _xy_divrem_rho2(hx)
    qy, ry = _xy_divrem_rho2(hy)
    dx = _xy_add(_xy_add(_xy_dx(txx), _xy_dy(txy)), qx)
    dy = _xy_add(_xy_add(_xy_dx(txy), _xy_dy(tyy)), qy)
    return (; dx, dy, dx_rem=rx, dy_rem=ry, tau33=t33)
end

"""
    _min_degree(a)

Return the lowest total degree of the terms of a polynomial in `(X, Y)`, `typemax(Int)`
for the zero polynomial.

Parameters
----------
a : Dict{Tuple{Int,Int},Rational{BigInt}}
    Terms.

Returns
-------
Int
"""
_min_degree(a::XYPoly) = isempty(a) ? typemax(Int) : minimum(sum(k) for k in keys(a))

"""
    _case_b_polynomials(power)

Return the exact polynomials of case `:B`: the potential `psi = X Y (1 - rho^2)^power`, the
velocity `u` with `r u = curl(r^3 psi e_z)`, that is
`ux = 3 Y psi + rho^2 d_Y psi`, `uy = -(3 X psi + rho^2 d_X psi)`, the force `f = -D(u)`
(`_stress_divergence`) as a polynomial plus a remainder over `rho^2`, and `tau_33`, in
units of `U0`, `eta U0 / R^2` and `eta U0 / R`. Every remainder term has a total degree of
at least 3, so the remainder over `rho^2` is at most a constant times `rho` and the force
is continuous with the value 0 at the centre.

Parameters
----------
power : Int
    Power of `(1 - rho^2)`, at least 2.

Returns
-------
NamedTuple
    `ux`, `uy`, `fx`, `fy`, `fx_rem`, `fy_rem`, `tau33` as exact polynomials.

Raises
------
ArgumentError
    If `power < 2`, or a remainder term has a total degree below 3.
"""
function _case_b_polynomials(power::Int)
    power >= 2 || throw(ArgumentError("the power must be at least 2, got $power"))
    w = _xy_add(XYPoly((0, 0) => QB(1)), _xy_scale(XY_RHO2, -1))
    wp = XYPoly((0, 0) => QB(1))
    for _ in 1:power
        wp = _xy_mul(wp, w)
    end
    phi = _xy_mul(_xy_mul(XY_X, XY_Y), wp)
    ux = _xy_add(_xy_scale(_xy_mul(XY_Y, phi), 3), _xy_mul(XY_RHO2, _xy_dy(phi)))
    uy = _xy_scale(
        _xy_add(_xy_scale(_xy_mul(XY_X, phi), 3), _xy_mul(XY_RHO2, _xy_dx(phi))), -1
    )
    d = _stress_divergence(ux, uy)
    min(_min_degree(d.dx_rem), _min_degree(d.dy_rem)) >= 3 ||
        throw(ArgumentError("the force remainder over rho^2 does not vanish at the centre"))
    return (;
        ux,
        uy,
        fx=_xy_scale(d.dx, -1),
        fy=_xy_scale(d.dy, -1),
        fx_rem=_xy_scale(d.dx_rem, -1),
        fy_rem=_xy_scale(d.dy_rem, -1),
        tau33=d.tau33,
    )
end

const CASE_B_EXACT = _case_b_polynomials(B_POWER)
const CASE_B_TERMS = (
    ux=_xy_terms(CASE_B_EXACT.ux),
    uy=_xy_terms(CASE_B_EXACT.uy),
    fx=_xy_terms(CASE_B_EXACT.fx),
    fy=_xy_terms(CASE_B_EXACT.fy),
    fx_rem=_xy_terms(CASE_B_EXACT.fx_rem),
    fy_rem=_xy_terms(CASE_B_EXACT.fy_rem),
    tau33=_xy_terms(CASE_B_EXACT.tau33),
)

"""
    _case_c_polynomials(power, lam, zr)

Return the exact polynomials of case `:C`, the radial compaction flow `u = grad Phi` with
`Phi = (1 - rho^2)^power`, in units of `U0` (velocity), `U0 / R` (strain rate),
`eta U0 / R` (stress and pressure), `eta U0 / R^2` (force) and `eta U0^2 / R^2` (heating).
The fields follow the continuum equations of the rows with `spherical_flow = true`:

- `e33 = (X ux + Y uy) / rho^2`, an exact polynomial (`_xy_div_rho2`), and
  `div3 = d_X ux + d_Y uy + e33`;
- `tau_ij = d_i u_j + d_j u_i - (2 / 3) delta_ij div3` and `tau33 = 2 (e33 - div3 / 3)`;
- `pf = lam Phi`, so that the Darcy flux `-(k / mu) grad P_f` is `-u`, and
  `pt = pf - zr div3` from the solid continuity `div3 + (P_t - P_f) / Z = 0`;
- the force `f = -(d_j tau_ij + (tau_ij X_j - tau33 X_i) / rho^2 - d_i pt)`, where `rho^2`
  divides the hoop numerator exactly for a radial flow;
- `eii2 = EXX^2 + EXY^2 + (3 / 4) e33d^2` with `EXX = (d_X ux - d_Y uy) / 2`,
  `EXY = (d_Y ux + d_X uy) / 2` and `e33d = (2 e33 - d_X ux - d_Y uy) / 3`, the square of
  the second invariant as the stress update defines it;
- `hs = 4 eii2 + zr div3^2 + lam |u|^2`, the deviatoric, compaction and Darcy heating of
  `compute_shear_heating!` with `S33 = tau33`.

Parameters
----------
power : Int
    Power of `(1 - rho^2)` in `Phi`, at least 4, so that the force is 0 on `rho = 1`.
lam : Rational{BigInt}
    `mu R^2 / (k eta)`, the fluid pressure unit `mu U0 R / k` over `eta U0 / R`.
zr : Rational{BigInt}
    `Z / eta = (1 - phi) zeta / eta`.

Returns
-------
NamedTuple
    `ux`, `uy`, `e33`, `div3`, `tau33`, `pf`, `pt`, `fx`, `fy`, `eii2`, `hs`, each a
    `Dict{Tuple{Int,Int},Rational{BigInt}}`.

Raises
------
ArgumentError
    If `power < 4`, or a division by `rho^2` leaves a remainder.
"""
function _case_c_polynomials(power::Int, lam::QB, zr::QB)
    power >= 4 || throw(ArgumentError("the power must be at least 4, got $power"))
    w = _xy_add(XYPoly((0, 0) => QB(1)), _xy_scale(XY_RHO2, -1))
    phi = XYPoly((0, 0) => QB(1))
    for _ in 1:power
        phi = _xy_mul(phi, w)
    end
    ux, uy = _xy_dx(phi), _xy_dy(phi)
    dxux, dyuy = _xy_dx(ux), _xy_dy(uy)
    e33 = _xy_div_rho2(_xy_add(_xy_mul(ux, XY_X), _xy_mul(uy, XY_Y)))
    div3 = _xy_add(_xy_add(dxux, dyuy), e33)
    trace = _xy_scale(div3, -2 // 3)
    txx = _xy_add(_xy_scale(dxux, 2), trace)
    tyy = _xy_add(_xy_scale(dyuy, 2), trace)
    txy = _xy_add(_xy_dy(ux), _xy_dx(uy))
    tau33 = _xy_add(_xy_scale(e33, 2), trace)
    pf = _xy_scale(phi, lam)
    pt = _xy_add(pf, _xy_scale(div3, -zr))
    hx = _xy_add(
        _xy_add(_xy_mul(txx, XY_X), _xy_mul(txy, XY_Y)), _xy_scale(_xy_mul(tau33, XY_X), -1)
    )
    hy = _xy_add(
        _xy_add(_xy_mul(txy, XY_X), _xy_mul(tyy, XY_Y)), _xy_scale(_xy_mul(tau33, XY_Y), -1)
    )
    dx = _xy_add(_xy_add(_xy_dx(txx), _xy_dy(txy)), _xy_div_rho2(hx))
    dy = _xy_add(_xy_add(_xy_dx(txy), _xy_dy(tyy)), _xy_div_rho2(hy))
    fx = _xy_scale(_xy_add(dx, _xy_scale(_xy_dx(pt), -1)), -1)
    fy = _xy_scale(_xy_add(dy, _xy_scale(_xy_dy(pt), -1)), -1)
    exx = _xy_scale(_xy_add(dxux, _xy_scale(dyuy, -1)), 1 // 2)
    exy = _xy_scale(txy, 1 // 2)
    e33d = _xy_scale(_xy_add(_xy_scale(e33, 2), _xy_scale(_xy_add(dxux, dyuy), -1)), 1 // 3)
    eii2 = _xy_add(
        _xy_add(_xy_mul(exx, exx), _xy_mul(exy, exy)),
        _xy_scale(_xy_mul(e33d, e33d), 3 // 4),
    )
    u2 = _xy_add(_xy_mul(ux, ux), _xy_mul(uy, uy))
    hs = _xy_add(
        _xy_add(_xy_scale(eii2, 4), _xy_scale(_xy_mul(div3, div3), zr)), _xy_scale(u2, lam)
    )
    return (; ux, uy, e33, div3, tau33, pf, pt, fx, fy, eii2, hs)
end

const CASE_C_EXACT = _case_c_polynomials(
    C_POWER,
    rationalize(BigInt, MU_F * R_PLANET^2 / (K_PERM * ETA_ROCK)),
    rationalize(BigInt, (1 - PHI) / ETA_RATIO_C),
)
const CASE_C_TERMS = NamedTuple{keys(CASE_C_EXACT)}(map(_xy_terms, values(CASE_C_EXACT)))

"""
    spherical_divergence_is_zero(power=B_POWER)

Return whether the case `:B` velocity has the spherical divergence 0 exactly:
`rho^2 (d_X ux + d_Y uy) + X ux + Y uy = 0` as a polynomial with rational coefficients.

Parameters
----------
power : Int, optional
    Power of `(1 - rho^2)` in the potential, at least 2 (default `B_POWER`).

Returns
-------
Bool

Raises
------
ArgumentError
    From `_case_b_polynomials` if `power < 2`.
"""
function spherical_divergence_is_zero(power::Int=B_POWER)
    p = power == B_POWER ? CASE_B_EXACT : _case_b_polynomials(power)
    div2 = _xy_add(_xy_dx(p.ux), _xy_dy(p.uy))
    total = _xy_add(
        _xy_mul(XY_RHO2, div2), _xy_add(_xy_mul(XY_X, p.ux), _xy_mul(XY_Y, p.uy))
    )
    return isempty(total)
end

"""
    velocity_vanishes_on_circle(case)

Return whether the exact velocity polynomial of `case` (with the option on) vanishes
exactly at the rational points `(3/5, 4/5)`, `(-5/13, 12/13)`, `(8/17, -15/17)` and
`(1, 0)` of the unit circle, so that `u = 0` on `r = R` (condition C5).

Parameters
----------
case : Symbol
    `:A`, `:A0`, `:B` or `:C`.

Returns
-------
Bool

Raises
------
ArgumentError
    If `case` is not in `CASES`.
"""
function velocity_vanishes_on_circle(case::Symbol)
    _check_case(case)
    ux, uy = if case == :B
        (CASE_B_EXACT.ux, CASE_B_EXACT.uy)
    elseif case == :C
        (CASE_C_EXACT.ux, CASE_C_EXACT.uy)
    else
        om = _xy_from_s(omega_coefficients(case; spherical_flow=true))
        (_xy_scale(_xy_mul(XY_Y, om), -1), _xy_mul(XY_X, om))
    end
    points = (
        (QB(3, 5), QB(4, 5)),
        (QB(-5, 13), QB(12, 13)),
        (QB(8, 17), QB(-15, 17)),
        (QB(1), QB(0)),
    )
    return all(
        iszero(_xy_eval_exact(ux, X, Y)) && iszero(_xy_eval_exact(uy, X, Y)) for
        (X, Y) in points
    )
end

"""
    exact_velocity(case, X, Y; spherical_flow=true)

Return the exact velocity of `case` at `(X, Y) = (x - xc, y - yc) / R_PLANET` in units of
`velocity_scale(case)`, `(0, 0)` for `X^2 + Y^2 >= 1`. The azimuthal cases have
`(vx, vy) = Omega (-Y, X)` with `Omega = omega_hat`.

Parameters
----------
case : Symbol
    `:A`, `:A0`, `:B` or `:C`.
X, Y : Real
    Nondimensional coordinates, finite.
spherical_flow : Bool, optional
    Value of `geometry.spherical_flow` (default `true`); only case `:A` has a solution for
    `false`.

Returns
-------
Tuple{Float64,Float64}
    `(vx, vy)`, dimensionless.

Raises
------
ArgumentError
    From `_check_case`.
DomainError
    If `X` or `Y` is not finite.
"""
function exact_velocity(case::Symbol, X::Real, Y::Real; spherical_flow::Bool=true)
    _check_case(case; spherical_flow=spherical_flow)
    _check_point(X, Y)
    x, y = Float64(X), Float64(Y)
    s = x^2 + y^2
    s >= 1 && return (0.0, 0.0)
    if case in (:B, :C)
        t = case == :B ? CASE_B_TERMS : CASE_C_TERMS
        return (_xy_eval(t.ux, x, y), _xy_eval(t.uy, x, y))
    end
    om = _seval(omega_coefficients(case; spherical_flow=spherical_flow), s)
    return (-y * om, x * om)
end

"""
    exact_force(case, X, Y; spherical_flow=true)

Return the body force of `case` at `(X, Y)` in units of `force_scale(case)`, `(0, 0)` for
`X^2 + Y^2 >= 1`. The azimuthal cases have `(fx, fy) = g(s) (-Y, X)` with
`g = force_coefficients`; case `:B` has `-D(u)` of `_case_b_polynomials`, a polynomial
plus a remainder over `rho^2` that tends to 0 at the centre, where the force is `(0, 0)`;
case `:C` has the polynomial force of `_case_c_polynomials`, `(0, 0)` at the centre.

Parameters
----------
case : Symbol
    `:A`, `:A0`, `:B` or `:C`.
X, Y : Real
    Nondimensional coordinates, finite.
spherical_flow : Bool, optional
    Value of `geometry.spherical_flow` (default `true`). Case `:A` has the same force for
    both settings.

Returns
-------
Tuple{Float64,Float64}
    `(fx, fy)`, dimensionless.

Raises
------
ArgumentError
    From `_check_case`.
DomainError
    If `X` or `Y` is not finite.
"""
function exact_force(case::Symbol, X::Real, Y::Real; spherical_flow::Bool=true)
    _check_case(case; spherical_flow=spherical_flow)
    _check_point(X, Y)
    x, y = Float64(X), Float64(Y)
    s = x^2 + y^2
    s >= 1 && return (0.0, 0.0)
    if case == :B
        # The remainder terms have degree >= 3, so their quotient by s tends to 0 at s = 0.
        fx = _xy_eval(CASE_B_TERMS.fx, x, y)
        fy = _xy_eval(CASE_B_TERMS.fy, x, y)
        s > 0 || return (fx, fy)
        return (
            fx + _xy_eval(CASE_B_TERMS.fx_rem, x, y) / s,
            fy + _xy_eval(CASE_B_TERMS.fy_rem, x, y) / s,
        )
    end
    case == :C && return (_xy_eval(CASE_C_TERMS.fx, x, y), _xy_eval(CASE_C_TERMS.fy, x, y))
    g = _seval(force_coefficients(case), s)
    return (-y * g, x * g)
end

"""
    exact_tau33(case, X, Y)

Return the exact out-of-plane deviatoric stress `tau_33 = 2 eta (u_r / r - div3 / 3)` of
`case` with the option on at `(X, Y)`, in units of `stress_scale(case) = eta U0 / R` for
cases `:B` (with `div3 = 0`) and `:C`; 0 for the azimuthal cases and for `X^2 + Y^2 >= 1`.
Case `:C` evaluates an exact polynomial, finite at the centre.

Parameters
----------
case : Symbol
    `:A`, `:A0`, `:B` or `:C`.
X, Y : Real
    Nondimensional coordinates, finite.

Returns
-------
Float64

Raises
------
ArgumentError
    From `_check_case`.
DomainError
    If `X` or `Y` is not finite.
"""
function exact_tau33(case::Symbol, X::Real, Y::Real)
    _check_case(case)
    _check_point(X, Y)
    x, y = Float64(X), Float64(Y)
    (case in (:B, :C) && x^2 + y^2 < 1) || return 0.0
    return _xy_eval((case == :B ? CASE_B_TERMS : CASE_C_TERMS).tau33, x, y)
end

"""
    _check_case_c(case)

Check that `case` is case `:C`, the only case with a nonzero divergence and pressures.

Parameters
----------
case : Symbol
    Case to check.

Returns
-------
Nothing

Raises
------
ArgumentError
    If `case` is not in `CASES` or not `:C`.
"""
function _check_case_c(case::Symbol)
    _check_case(case)
    case == :C || throw(ArgumentError("only case :C has this field, got :$case"))
    return nothing
end

"""
    _case_c_eval(case, field, X, Y)

Evaluate the case `:C` polynomial `field` of `_case_c_polynomials` at `(X, Y)`, 0 for
`X^2 + Y^2 >= 1`.

Parameters
----------
case : Symbol
    Must be `:C`.
field : Symbol
    Name of the polynomial in `CASE_C_TERMS`.
X, Y : Real
    Nondimensional coordinates, finite.

Returns
-------
Float64

Raises
------
ArgumentError
    If `case` is not in `CASES` or not `:C`.
DomainError
    If `X` or `Y` is not finite.
"""
function _case_c_eval(case::Symbol, field::Symbol, X::Real, Y::Real)
    _check_case_c(case)
    _check_point(X, Y)
    x, y = Float64(X), Float64(Y)
    x^2 + y^2 >= 1 && return 0.0
    return _xy_eval(getfield(CASE_C_TERMS, field), x, y)
end

"""
    exact_div3(case, X, Y)

Return the exact 3-D divergence `div3 = d_x ux + d_y uy + u_r / r` of case `:C` at `(X, Y)`
in units of `U0 / R`; `-30` at the centre and 0 for `X^2 + Y^2 >= 1`.

Parameters
----------
case : Symbol
    `:C`.
X, Y : Real
    Nondimensional coordinates, finite.

Returns
-------
Float64

Raises
------
ArgumentError
    If `case` is not `:C`.
DomainError
    If `X` or `Y` is not finite.
"""
exact_div3(case::Symbol, X::Real, Y::Real) = _case_c_eval(case, :div3, X, Y)

"""
    exact_pressures(case, X, Y)

Return the exact total and fluid pressure of case `:C` at `(X, Y)` relative to the ambient
pressure, `(P_t - p_s, P_f - p_s)`, in units of `stress_scale(:C)`; `(0, 0)` for
`X^2 + Y^2 >= 1`.

Parameters
----------
case : Symbol
    `:C`.
X, Y : Real
    Nondimensional coordinates, finite.

Returns
-------
Tuple{Float64,Float64}

Raises
------
ArgumentError
    If `case` is not `:C`.
DomainError
    If `X` or `Y` is not finite.
"""
function exact_pressures(case::Symbol, X::Real, Y::Real)
    return (_case_c_eval(case, :pt, X, Y), _case_c_eval(case, :pf, X, Y))
end

"""
    exact_eii(case, X, Y)

Return the exact second invariant `E_II = sqrt(EXX^2 + EXY^2 + 0.75 e33d^2)` of the
strain rate of case `:C` at `(X, Y)`, as the stress update defines it, in units of
`U0 / R`; 0 at the centre and for `X^2 + Y^2 >= 1`.

Parameters
----------
case : Symbol
    `:C`.
X, Y : Real
    Nondimensional coordinates, finite.

Returns
-------
Float64

Raises
------
ArgumentError
    If `case` is not `:C`.
DomainError
    If `X` or `Y` is not finite.
"""
exact_eii(case::Symbol, X::Real, Y::Real) = sqrt(max(_case_c_eval(case, :eii2, X, Y), 0.0))

"""
    exact_sii(case, X, Y)

Return the exact second invariant `S_II = 2 eta E_II` of the deviatoric stress of case `:C`
in the viscous limit at `(X, Y)`, in units of `stress_scale(:C)`.

Parameters
----------
case : Symbol
    `:C`.
X, Y : Real
    Nondimensional coordinates, finite.

Returns
-------
Float64

Raises
------
ArgumentError
    If `case` is not `:C`.
DomainError
    If `X` or `Y` is not finite.
"""
exact_sii(case::Symbol, X::Real, Y::Real) = 2 * exact_eii(case, X, Y)

"""
    exact_hs(case, X, Y)

Return the exact heating of `compute_shear_heating!` for case `:C` at `(X, Y)`: the
deviatoric part `4 eta E_II^2`, the compaction part `(P_t - P_f)^2 / Z` and the Darcy part
`(mu / k) |q|^2`, in units of `eta U0^2 / R^2`; 0 for `X^2 + Y^2 >= 1`.

Parameters
----------
case : Symbol
    `:C`.
X, Y : Real
    Nondimensional coordinates, finite.

Returns
-------
Float64

Raises
------
ArgumentError
    If `case` is not `:C`.
DomainError
    If `X` or `Y` is not finite.
"""
exact_hs(case::Symbol, X::Real, Y::Real) = _case_c_eval(case, :hs, X, Y)

"""
    _check_positive(name, v)

Check that a scale parameter is finite and positive.

Parameters
----------
name : String
    Name of the parameter for the message.
v : Real
    Value.

Returns
-------
Nothing

Raises
------
DomainError
    If `v` is not finite or not positive.
"""
function _check_positive(name::String, v::Real)
    (isfinite(v) && v > 0) || throw(DomainError(v, "$name must be finite and positive"))
    return nothing
end

"""
    velocity_scale(case; f0=F0, u0=U0, r_planet=R_PLANET, eta=ETA_ROCK)

Return the velocity unit of `exact_velocity`: `f0 r_planet^2 / eta` for the azimuthal
cases and `u0` for cases `:B` and `:C`.

Parameters
----------
case : Symbol
    `:A`, `:A0`, `:B` or `:C`.
f0 : Real, optional
    Force scale [N/m^3], finite (default `F0`).
u0 : Real, optional
    Velocity scale of cases `:B` and `:C` [m/s], finite (default `U0`).
r_planet : Real, optional
    Radius [m], finite and positive (default `R_PLANET`).
eta : Real, optional
    Shear viscosity [Pa s], finite and positive (default `ETA_ROCK`).

Returns
-------
Float64
    Velocity unit [m/s].

Raises
------
ArgumentError
    If `case` is not in `CASES`.
DomainError
    If `f0` or `u0` is not finite, or `r_planet` or `eta` is not finite and positive.
"""
function velocity_scale(
    case::Symbol; f0::Real=F0, u0::Real=U0, r_planet::Real=R_PLANET, eta::Real=ETA_ROCK
)
    _check_case(case)
    isfinite(f0) || throw(DomainError(f0, "F0 must be finite"))
    isfinite(u0) || throw(DomainError(u0, "U0 must be finite"))
    _check_positive("r_planet", r_planet)
    _check_positive("eta", eta)
    return case in (:B, :C) ? Float64(u0) : Float64(f0) * Float64(r_planet)^2 / Float64(eta)
end

"""
    force_scale(case; f0=F0, u0=U0, r_planet=R_PLANET, eta=ETA_ROCK)

Return the force unit of `exact_force`: `f0` for the azimuthal cases and
`eta u0 / r_planet^2` for cases `:B` and `:C`.

Parameters
----------
case : Symbol
    `:A`, `:A0`, `:B` or `:C`.
f0, u0, r_planet, eta : Real, optional
    As for `velocity_scale`.

Returns
-------
Float64
    Force unit [N/m^3].

Raises
------
ArgumentError
    If `case` is not in `CASES`.
DomainError
    As for `velocity_scale`.
"""
function force_scale(
    case::Symbol; f0::Real=F0, u0::Real=U0, r_planet::Real=R_PLANET, eta::Real=ETA_ROCK
)
    v = velocity_scale(case; f0=f0, u0=u0, r_planet=r_planet, eta=eta)
    return Float64(eta) * v / Float64(r_planet)^2
end

"""
    stress_scale(case; f0=F0, u0=U0, r_planet=R_PLANET, eta=ETA_ROCK)

Return the stress unit of the case, `eta velocity_scale / r_planet`: `f0 r_planet` for the
azimuthal cases and `eta u0 / r_planet` for cases `:B` and `:C`.

Parameters
----------
case : Symbol
    `:A`, `:A0`, `:B` or `:C`.
f0, u0, r_planet, eta : Real, optional
    As for `velocity_scale`.

Returns
-------
Float64
    Stress unit [Pa].

Raises
------
ArgumentError
    If `case` is not in `CASES`.
DomainError
    As for `velocity_scale`.
"""
function stress_scale(
    case::Symbol; f0::Real=F0, u0::Real=U0, r_planet::Real=R_PLANET, eta::Real=ETA_ROCK
)
    v = velocity_scale(case; f0=f0, u0=u0, r_planet=r_planet, eta=eta)
    return Float64(eta) * v / Float64(r_planet)
end

"""
    profile_ratio_squared()

Return the exact square of the on/off ratio of the L2 norms of the case `:A` profile,
`int_0^1 rho u_on^2 drho / int_0^1 rho u_off^2 drho` with `u = rho Omega`. With
`s = rho^2` each integral is `(1 / 2) int_0^1 s Omega(s)^2 ds`.

Returns
-------
Rational{BigInt}
"""
function profile_ratio_squared()
    # Entry k of the coefficients holds s^(k - 1), whose integral over [0, 1] is 1 / k.
    integral(c) =
        sum(v // k for (k, v) in enumerate(_spoly_mul(QB[0, 1], _spoly_mul(c, c))))
    return integral(OMEGA_A_ON) / integral(OMEGA_A_OFF)
end

"""
    profile_norm_ratio()

Return the exact on/off ratio of the L2 norms of the case `:A` profile,
`sqrt(profile_ratio_squared())`, rounded to Float64.

Returns
-------
Float64
"""
profile_norm_ratio() = Float64(sqrt(BigFloat(profile_ratio_squared())))

"""
    make_config(N; spherical_flow, darcy_elimination=true, case=:A)

Return the `SimulationConfig` of the case on an `N` x `N` grid: uniform viscosity
`ETA_ROCK`, `eta / zeta = eta_zeta_ratio(case)` through
`solver.etaphikoef = PHI / eta_zeta_ratio(case)`, the time step `DT_YR`, the surface
pressure `ambient_pressure(case)` and no venting. Cases `:A`, `:A0` and `:B` share one
configuration.

Parameters
----------
N : Int
    Basic grid resolution in x and y, odd so that the centre lies on a basic node, and at
    least `NX_MIN`, the smallest grid with `FIT_DEGREE + 2` distinct radii in the window
    of `centre_fit_rotation`.
spherical_flow : Bool
    Value of `geometry.spherical_flow`.
darcy_elimination : Bool, optional
    Value of `solver.darcy_elimination`: `true` solves the 4-variable system, `false` the
    6-variable system (default `true`).
case : Symbol, optional
    `:A`, `:A0`, `:B` or `:C` (default `:A`).

Returns
-------
SimulationConfig

Raises
------
ArgumentError
    If `N < NX_MIN` or `N` is even, `case` is not in `CASES`, or `case` has no exact
    solution for `spherical_flow`.
"""
function make_config(
    N::Int; spherical_flow::Bool, darcy_elimination::Bool=true, case::Symbol=:A
)
    _check_case(case; spherical_flow=spherical_flow)
    N >= NX_MIN ||
        throw(ArgumentError("N must be at least $NX_MIN to resolve the centre fit, got $N"))
    isodd(N) || throw(
        ArgumentError("N must be odd so that the centre lies on a basic node, got $N")
    )
    eta3 = SVector{3,Float64}(ETA_ROCK, ETA_ROCK, ETA_ROCK)
    mu3 = SVector{3,Float64}(MU_F, MU_F, MU_F)
    return SimulationConfig(;
        grid=GridConfig(; xsize=L_BOX, ysize=L_BOX, Nx=N, Ny=N),
        geometry=GeometryConfig(;
            rplanet=R_PLANET,
            rcrust=R_PLANET,
            xcenter=X_CENTRE,
            ycenter=Y_CENTRE,
            psurface=ambient_pressure(case),
            spherical_metric=true,
            spherical_flow=spherical_flow,
            gravity_mode=:enclosed_mass,
        ),
        time=TimeConfig(; dt_initial=DT_YR, dt_longest=1.0e12, n_steps=1),
        solver=SolverConfig(;
            darcy_elimination=darcy_elimination,
            etamax=ETA_MAX,
            etamin=ETA_MIN,
            etaphikoef=PHI / eta_zeta_ratio(case),
        ),
        materials=MaterialConfig(;
            etasolidm=eta3,
            etasolidmm=eta3,
            etafluidm=mu3,
            etafluidmm=mu3,
            kphim0=SVector{3,Float64}(K_PERM, K_PERM, 1.0e-30),
            tkm0=SVector{3,Float64}(T_SURF, T_SURF, T_SURF),
        ),
        thermodynamics=ThermalConfig(;
            phim0=PHI,
            tmfluidphase=100.0,
            fluid_viscosity_mode=:constant,
            fluid_viscosity_eta0=MU_F,
            thermal_buoyancy=false,
        ),
        poroelasticity=PoroelasticConfig(;
            hydrofracture=false, betasolid=0.0, betafluid=0.0
        ),
        venting=VentingConfig(;
            active=false, mode=:darcy_sink, species=:H2O, k_vent=1.0e-6, ice_sealing=false
        ),
        reaction=ReactionConfig(; active=false, fluid_overpressure_coupling=true),
        disk=DiskConfig(;
            enabled=true, model=:fixed, t_ambient=T_SURF, p_amb_disk=ambient_pressure(case)
        ),
        output=OutputConfig(; save_final=false, savematstep=100, visstep=100),
    )
end

"""
    grid_of(coords)

Return the node coordinates that the error measures read, copied from `coords`.

Parameters
----------
coords : GridCoordinates
    Staggered grid coordinates.

Returns
-------
NamedTuple
    `Nx1`, `Ny1`, `dx`, `dy` [m] and the vectors `xvx`, `yvx`, `xvy`, `yvy`, `xp`, `yp` [m].
"""
function grid_of(coords)
    return (;
        Nx1=coords.Nx1,
        Ny1=coords.Ny1,
        dx=coords.dx,
        dy=coords.dy,
        xvx=collect(coords.xvx),
        yvx=collect(coords.yvx),
        xvy=collect(coords.xvy),
        yvy=collect(coords.yvy),
        xp=collect(coords.xp),
        yp=collect(coords.yp),
    )
end

"""
    _nd(x, y)

Return the nondimensional position `((x - X_CENTRE) / R_PLANET, (y - Y_CENTRE) / R_PLANET)`.

Parameters
----------
x, y : Float64
    Position [m].

Returns
-------
Tuple{Float64,Float64}
"""
_nd(x, y) = ((x - X_CENTRE) / R_PLANET, (y - Y_CENTRE) / R_PLANET)

"""
    set_fields!(state, coords, case)

Set the grid fields of the case: uniform viscosity `ETA_ROCK` (basic nodes `ETA`, `ETA0`
and P nodes `ETAP`), shear modulus `SHEAR_MODULUS`, strengths `STRENGTH_FILL` without
friction, porosity `PHI`, no compressibility, no stress history (`SXX0`, `SXY0`, `S330`),
density `RHO_SOLID` and fluid density `RHO_FLUID`, pressures `ambient_pressure(case)`, the
face resistance `MU_F / K_PERM`, and the body force of `case` as `gx = fx / RHOX` and
`gy = fy / RHOY` at the vx and vy nodes.

Parameters
----------
state : SimulationState
    State from `init_simulation`; its grid fields are overwritten.
coords : GridCoordinates
    Staggered grid coordinates of `state`.
case : Symbol
    `:A`, `:A0`, `:B` or `:C`. Case `:A` has the same force with and without the option.

Returns
-------
Nothing

Raises
------
ArgumentError
    If `case` is not in `CASES`.
"""
function set_fields!(state, coords, case::Symbol)
    _check_case(case)
    g = state.grids
    for f in (g.ETA, g.ETA0, g.ETAP)
        fill!(f, ETA_ROCK)
    end
    fill!(g.FRI, 0.0)
    for f in (g.COH, g.TEN)
        fill!(f, STRENGTH_FILL)
    end
    for f in (g.GGG, g.GGGP)
        fill!(f, SHEAR_MODULUS)
    end
    fill!(g.tk1, T_SURF)
    fill!(g.tk2, T_SURF)
    for f in (g.PHI, g.PHIX, g.PHIY)
        fill!(f, PHI)
    end
    for f in (g.BETAPHI, g.DQPF, g.DMP, g.SXX0, g.SXY0, g.S330, g.gx, g.gy)
        fill!(f, 0.0)
    end
    for f in (g.RHO, g.RHOX, g.RHOY)
        fill!(f, RHO_SOLID)
    end
    fill!(g.RHOCP, RHOCP)
    fill!(g.RHOFX, RHO_FLUID)
    fill!(g.RHOFY, RHO_FLUID)
    for f in (g.pr, g.pr0, g.pf, g.pf0)
        fill!(f, ambient_pressure(case))
    end
    fill!(g.RX, MU_F / K_PERM)
    fill!(g.RY, MU_F / K_PERM)
    fs = force_scale(case)
    for j in 1:(coords.Nx1), i in 1:(coords.Ny1)
        g.gx[i, j] =
            fs * exact_force(case, _nd(coords.xvx[j], coords.yvx[i])...)[1] / g.RHOX[i, j]
        g.gy[i, j] =
            fs * exact_force(case, _nd(coords.xvy[j], coords.yvy[i])...)[2] / g.RHOY[i, j]
    end
    return nothing
end

"""
    _init_case(cfg, out, case)

Build the state of the case with `init_simulation`, set timestep 1, the ambient pressure
`ambient_pressure(case)`, the time step of `cfg` in seconds and the grid fields of `set_fields!`.

Parameters
----------
cfg : SimulationConfig
    Configuration from `make_config`.
out : String
    Output folder of `init_simulation`.
case : Symbol
    `:A`, `:A0`, `:B` or `:C`.

Returns
-------
Tuple
    The result of `init_simulation`: `(state, coords, ws, tio, rest)`.

Raises
------
ArgumentError
    If `case` is not in `CASES`.
"""
function _init_case(cfg, out, case::Symbol)
    _check_case(case)
    state, coords, ws, tio, rest = with_logger(NullLogger()) do
        return Erebus.init_simulation(cfg; output_path=out)
    end
    state.timestep = 1
    state.accumulators.P_amb = ambient_pressure(case)
    state.dt = cfg.time.dt_initial * cfg.time.yearlength
    set_fields!(state, coords, case)
    return (state, coords, ws, tio, rest)
end

"""
    _vx_nodes(grid)

Return the interior vx nodes `(i, j)` of a grid: rows `2:Ny1-1` and columns `2:Nx1-2`,
without the ghost row and the wall columns.

Parameters
----------
grid : NamedTuple
    Grid from `grid_of`.

Returns
-------
Iterator of Tuple{Int,Int}
"""
_vx_nodes(grid) = ((i, j) for j in 2:(grid.Nx1 - 2) for i in 2:(grid.Ny1 - 1))

"""
    _vy_nodes(grid)

Return the interior vy nodes `(i, j)` of a grid: rows `2:Ny1-2` and columns `2:Nx1-1`.

Parameters
----------
grid : NamedTuple
    Grid from `grid_of`.

Returns
-------
Iterator of Tuple{Int,Int}
"""
_vy_nodes(grid) = ((i, j) for j in 2:(grid.Nx1 - 1) for i in 2:(grid.Ny1 - 2))

"""
    exact_fields(grid, case; spherical_flow=true)

Return the exact velocity of `case` [m/s] at the vx and vy nodes of `grid`.

Parameters
----------
grid : NamedTuple
    Grid from `grid_of`.
case : Symbol
    `:A`, `:A0`, `:B` or `:C`.
spherical_flow : Bool, optional
    Value of `geometry.spherical_flow` (default `true`).

Returns
-------
NamedTuple
    `vx`, `vy`: arrays of size `(Ny1, Nx1)` [m/s].

Raises
------
ArgumentError
    From `_check_case`.
"""
function exact_fields(grid, case::Symbol; spherical_flow::Bool=true)
    _check_case(case; spherical_flow=spherical_flow)
    v = velocity_scale(case)
    vx = [
        v * exact_velocity(case, _nd(grid.xvx[j], grid.yvx[i])...; spherical_flow)[1] for
        i in 1:(grid.Ny1), j in 1:(grid.Nx1)
    ]
    vy = [
        v * exact_velocity(case, _nd(grid.xvy[j], grid.yvy[i])...; spherical_flow)[2] for
        i in 1:(grid.Ny1), j in 1:(grid.Nx1)
    ]
    return (; vx, vy)
end

"""
    exact_s33(grid, case)

Return the exact out-of-plane deviatoric stress `tau_33` of `case` [Pa] at the P nodes of
`grid` (see `exact_tau33`).

Parameters
----------
grid : NamedTuple
    Grid from `grid_of`.
case : Symbol
    `:A`, `:A0`, `:B` or `:C`.

Returns
-------
Matrix{Float64}
    Size `(Ny1, Nx1)` [Pa].

Raises
------
ArgumentError
    From `_check_case`.
"""
function exact_s33(grid, case::Symbol)
    _check_case(case)
    s = stress_scale(case)
    return [
        s * exact_tau33(case, _nd(grid.xp[j], grid.yp[i])...) for
        i in 1:(grid.Ny1), j in 1:(grid.Nx1)
    ]
end

"""
    velocity_errors(vx, vy, grid, case; spherical_flow=true)

Return the relative errors of a velocity field against the exact flow of `case` over the
interior vx and vy nodes with equal weights:
`l2 = sqrt(sum (v - v_exact)^2 / sum v_exact^2)` over all nodes, `l2_inside` and
`l2_outside` with the numerator restricted to the nodes with `r < R` and `r >= R`
(`l2^2 = l2_inside^2 + l2_outside^2`), and `linf = max |v - v_exact| / max |v_exact|`.

Parameters
----------
vx, vy : AbstractMatrix{Float64}
    Velocity components [m/s] of size `(Ny1, Nx1)`.
grid : NamedTuple
    Grid from `grid_of`.
case : Symbol
    `:A`, `:A0`, `:B` or `:C`.
spherical_flow : Bool, optional
    Setting of the exact flow (default `true`).

Returns
-------
NamedTuple
    `l2`, `l2_inside`, `l2_outside`, `linf`, dimensionless.

Raises
------
ArgumentError
    From `_check_case`.
"""
function velocity_errors(vx, vy, grid, case::Symbol; spherical_flow::Bool=true)
    ex = exact_fields(grid, case; spherical_flow=spherical_flow)
    se_in = se_out = sr = emax = vmax = 0.0
    for (v, e, xs, ys, nodes) in (
        (vx, ex.vx, grid.xvx, grid.yvx, _vx_nodes(grid)),
        (vy, ex.vy, grid.xvy, grid.yvy, _vy_nodes(grid)),
    )
        for (i, j) in nodes
            d2 = (v[i, j] - e[i, j])^2
            X, Y = _nd(xs[j], ys[i])
            X^2 + Y^2 < 1 ? (se_in += d2) : (se_out += d2)
            sr += e[i, j]^2
            emax = max(emax, abs(v[i, j] - e[i, j]))
            vmax = max(vmax, abs(e[i, j]))
        end
    end
    return (;
        l2=sqrt((se_in + se_out) / sr),
        l2_inside=sqrt(se_in / sr),
        l2_outside=sqrt(se_out / sr),
        linf=emax / vmax,
    )
end

"""
    centre_node_rotation(vx, vy, grid; scale=velocity_scale(:A))

Return the rotation rate at the basic node on the centre, half the vorticity
`(d_x vy - d_y vx) / 2` from the 4 velocity nodes around it, in units of
`scale / R_PLANET`.

Parameters
----------
vx, vy : AbstractMatrix{Float64}
    Velocity components [m/s].
grid : NamedTuple
    Grid from `grid_of`; the centre must lie on a basic node (odd `N`).
scale : Float64, optional
    Velocity unit [m/s] (default `velocity_scale(:A)`), finite and positive.

Returns
-------
Float64
    Rotation rate, dimensionless.

Raises
------
DomainError
    If `scale` is not finite or not positive.
ArgumentError
    If no basic node lies on the centre.
"""
function centre_node_rotation(vx, vy, grid; scale::Float64=velocity_scale(:A))
    _check_positive("scale", scale)
    bx = argmin(abs.(grid.xvx .- X_CENTRE))
    ay = argmin(abs.(grid.yvy .- Y_CENTRE))
    a = searchsortedlast(grid.yvx, Y_CENTRE)
    b = searchsortedlast(grid.xvy, X_CENTRE)
    aligned =
        abs(grid.xvx[bx] - X_CENTRE) < 1.0e-6 * grid.dx &&
        abs(grid.yvy[ay] - Y_CENTRE) < 1.0e-6 * grid.dy &&
        1 <= a < length(grid.yvx) &&
        1 <= b < length(grid.xvy)
    aligned || throw(ArgumentError("no basic node lies on the centre"))
    om = ((vy[ay, b + 1] - vy[ay, b]) / grid.dx - (vx[a + 1, bx] - vx[a, bx]) / grid.dy) / 2
    return om * R_PLANET / scale
end

"""
    centre_fit_rotation(vx, vy, grid; scale=velocity_scale(:A), rho_min=FIT_RHO_MIN, rho_max=FIT_RHO_MAX, degree=FIT_DEGREE)

Return `Omega(0)` from a least-squares fit of `Omega(rho) = sum_{k=0..degree} c_k rho^(2 k)`
to the velocity nodes on the 2 axes through the centre with `rho_min <= rho <= rho_max`:
`Omega = vy / (x - xc)` on the vy nodes of the row `y = yc` and `-vx / (y - yc)` on the vx
nodes of the column `x = xc`, where the velocity is purely azimuthal for an azimuthal flow.
The 2 axes sample the same radii, so the fit needs at least `degree + 2` distinct radii in
the range. The exact `Omega` of the azimuthal cases has degree 5 in `rho^2`, so the fit
with the default degree has no bias.

Parameters
----------
vx, vy : AbstractMatrix{Float64}
    Velocity components [m/s].
grid : NamedTuple
    Grid from `grid_of`; the centre must lie on a basic node (odd `N`).
scale : Float64, optional
    Velocity unit [m/s] (default `velocity_scale(:A)`), finite and positive.
rho_min, rho_max : Float64, optional
    Radius range of the nodes, finite with `0 < rho_min < rho_max <= 1` (defaults
    `FIT_RHO_MIN`, `FIT_RHO_MAX`).
degree : Int, optional
    Degree of the fit in `rho^2`, non-negative (default `FIT_DEGREE`).

Returns
-------
Float64
    `Omega(0)` in units of `scale / R_PLANET`.

Raises
------
DomainError
    If `scale` is not finite or not positive.
ArgumentError
    If the range is not finite with `0 < rho_min < rho_max <= 1`, `degree < 0`, no basic
    node lies on the centre, or fewer than `degree + 2` distinct radii lie in the range.
"""
function centre_fit_rotation(
    vx,
    vy,
    grid;
    scale::Float64=velocity_scale(:A),
    rho_min::Float64=FIT_RHO_MIN,
    rho_max::Float64=FIT_RHO_MAX,
    degree::Int=FIT_DEGREE,
)
    _check_positive("scale", scale)
    (isfinite(rho_min) && isfinite(rho_max) && 0 < rho_min < rho_max <= 1) || throw(
        ArgumentError("need finite 0 < rho_min < rho_max <= 1, got $rho_min, $rho_max")
    )
    degree >= 0 || throw(ArgumentError("degree must be non-negative, got $degree"))
    bx = argmin(abs.(grid.xvx .- X_CENTRE))
    ay = argmin(abs.(grid.yvy .- Y_CENTRE))
    (
        abs(grid.xvx[bx] - X_CENTRE) < 1.0e-6 * grid.dx &&
        abs(grid.yvy[ay] - Y_CENTRE) < 1.0e-6 * grid.dy
    ) || throw(ArgumentError("no basic node lies on the centre"))
    rhos = Float64[]
    oms = Float64[]
    for j in 2:(grid.Nx1 - 1)
        x = grid.xvy[j] - X_CENTRE
        rho = abs(x) / R_PLANET
        if rho_min <= rho <= rho_max
            push!(rhos, rho)
            push!(oms, vy[ay, j] / x)
        end
    end
    for i in 2:(grid.Ny1 - 1)
        y = grid.yvx[i] - Y_CENTRE
        rho = abs(y) / R_PLANET
        if rho_min <= rho <= rho_max
            push!(rhos, rho)
            push!(oms, -vx[i, bx] / y)
        end
    end
    n_radii = length(unique(round.(rhos; sigdigits=12)))
    n_radii >= degree + 2 || throw(
        ArgumentError("$n_radii distinct radii lie in the fit range, need $(degree + 2)"),
    )
    A = [rho^(2k) for rho in rhos, k in 0:degree]
    return (A \ oms)[1] * R_PLANET / scale
end

"""
    s33_error(S33, grid; case=:B)

Return the relative L2 error of the out-of-plane stress `S33` [Pa] against the exact
`tau_33` of `case` over the interior P nodes with `r < R`:
`sqrt(sum (S33 - tau_33)^2 / sum tau_33^2)`.

Parameters
----------
S33 : AbstractMatrix{Float64}
    `grids.S33` after the solve [Pa], size `(Ny1, Nx1)`.
grid : NamedTuple
    Grid from `grid_of`.
case : Symbol, optional
    Case with a nonzero `tau_33`, `:B` or `:C` (default `:B`).

Returns
-------
Float64
    Relative error, dimensionless.

Raises
------
ArgumentError
    If `case` is not `:B` or `:C` (the azimuthal cases have `tau_33 = 0`).
"""
function s33_error(S33, grid; case::Symbol=:B)
    case in (:B, :C) || throw(ArgumentError("only cases :B and :C have a nonzero tau_33"))
    return pnode_errors(S33, exact_s33(grid, case), grid).l2_inside
end

"""
    pnode_errors(field, exact, grid; offset=0.0)

Return the errors of a field at the P nodes against its exact values over the interior P
nodes: `l2_inside = sqrt(sum (field - offset - exact)^2 / sum exact^2)` over the nodes with
`r < R`, and `out = max |field - offset - exact| / max |exact|` with the numerator over the
nodes with `r >= R` and the denominator over the nodes with `r < R`.

Parameters
----------
field : AbstractMatrix{Float64}
    Solved field at the P nodes, size `(Ny1, Nx1)`.
exact : AbstractMatrix{Float64}
    Exact field relative to `offset`, same size and units.
grid : NamedTuple
    Grid from `grid_of`.
offset : Float64, optional
    Value subtracted from `field`, for example the ambient pressure (default `0.0`).

Returns
-------
NamedTuple
    `l2_inside`, `out`, dimensionless.

Raises
------
ArgumentError
    If `exact` is 0 on every node with `r < R`.
"""
function pnode_errors(field, exact, grid; offset::Float64=0.0)
    se = sr = emax_out = vmax = 0.0
    for j in 2:(grid.Nx1 - 1), i in 2:(grid.Ny1 - 1)
        X, Y = _nd(grid.xp[j], grid.yp[i])
        d = field[i, j] - offset - exact[i, j]
        if X^2 + Y^2 < 1
            se += d^2
            sr += exact[i, j]^2
            vmax = max(vmax, abs(exact[i, j]))
        else
            emax_out = max(emax_out, abs(d))
        end
    end
    sr > 0 || throw(ArgumentError("the exact field is 0 inside R"))
    return (; l2_inside=sqrt(se / sr), out=emax_out / vmax)
end

"""
    exact_pnode_fields(grid, case)

Return the exact fields of case `:C` at the P nodes of `grid`: `tau_33` [Pa], `P_t - p_s`
and `P_f - p_s` [Pa], `E_II` [1/s], `S_II` [Pa] and the heating `H_S` [W/m^3].

Parameters
----------
grid : NamedTuple
    Grid from `grid_of`.
case : Symbol
    `:C`.

Returns
-------
NamedTuple
    `s33`, `pt`, `pf`, `eii`, `sii`, `hs`: matrices of size `(Ny1, Nx1)`.

Raises
------
ArgumentError
    If `case` is not `:C`.
"""
function exact_pnode_fields(grid, case::Symbol)
    _check_case_c(case)
    sc = stress_scale(case)
    rate = velocity_scale(case) / R_PLANET
    at(f) =
        [f(case, _nd(grid.xp[j], grid.yp[i])...) for i in 1:(grid.Ny1), j in 1:(grid.Nx1)]
    p = at(exact_pressures)
    return (;
        s33=sc .* at(exact_tau33),
        pt=sc .* first.(p),
        pf=sc .* last.(p),
        eii=rate .* at(exact_eii),
        sii=sc .* at(exact_sii),
        hs=(sc * rate) .* at(exact_hs),
    )
end

"""
    anchor_deviation(pr, pf, grid, p_s)

Return the largest relative deviation of `P_t` and `P_f` from the ambient pressure on the
anchor rings, the P nodes `i in (2, Ny1 - 1)` or `j in (2, Nx1 - 1)` where the solver sets
both pressures to `geometry.psurface`.

Parameters
----------
pr, pf : AbstractMatrix{Float64}
    Total and fluid pressure [Pa].
grid : NamedTuple
    Grid from `grid_of`.
p_s : Float64
    Ambient pressure [Pa], finite and positive.

Returns
-------
Float64
    `max |P - p_s| / p_s` over both rings, dimensionless.

Raises
------
DomainError
    If `p_s` is not finite or not positive.
"""
function anchor_deviation(pr, pf, grid, p_s::Float64)
    _check_positive("p_s", p_s)
    Ny1, Nx1 = grid.Ny1, grid.Nx1
    dev = 0.0
    for j in 2:(Nx1 - 1), i in 2:(Ny1 - 1)
        (i in (2, Ny1 - 1) || j in (2, Nx1 - 1)) || continue
        dev = max(dev, abs(pr[i, j] - p_s), abs(pf[i, j] - p_s))
    end
    return dev / p_s
end

"""
    _case_c_heating!(g)

Compute the heating `g.HS` of `compute_shear_heating!` from the solved fields, with
`S33 = g.S33` as `solve_thermal_energy!` passes it with `spherical_flow = true`.

Parameters
----------
g : Grids
    Grid fields after the solve; `g.HS` is overwritten.

Returns
-------
Nothing
"""
function _case_c_heating!(g)
    Erebus.compute_shear_heating!(
        g.HS,
        g.ETA,
        g.SXY,
        g.ETAP,
        g.SXX,
        g.RX,
        g.RY,
        g.qxD,
        g.qyD,
        g.PHI,
        g.ETAPHI,
        g.pr,
        g.pf;
        S33=g.S33,
    )
    return nothing
end

"""
    pressure_deviation(pr, grid)

Return the largest `abs(pr - P_AMB)` [Pa] over the P nodes `3:Ny1-2` x `3:Nx1-2`, inside
the anchor ring of the total pressure. The exact total pressure is uniform.

Parameters
----------
pr : AbstractMatrix{Float64}
    Total pressure [Pa].
grid : NamedTuple
    Grid from `grid_of`.

Returns
-------
Float64
    Largest deviation [Pa].
"""
function pressure_deviation(pr, grid)
    return maximum(abs(pr[i, j] - P_AMB) for j in 3:(grid.Nx1 - 2), i in 3:(grid.Ny1 - 2))
end

"""
    profile_norm(vx, vy, grid)

Return `sqrt(sum v^2)` [m/s] over the interior vx and vy nodes with `r < R`. Each node
stands for the same area, so the ratio of 2 such norms approximates the ratio of the
continuum L2 norms `sqrt(int rho u^2 drho)`.

Parameters
----------
vx, vy : AbstractMatrix{Float64}
    Velocity components [m/s].
grid : NamedTuple
    Grid from `grid_of`.

Returns
-------
Float64
    Discrete norm [m/s].
"""
function profile_norm(vx, vy, grid)
    total = 0.0
    for (v, xs, ys, nodes) in (
        (vx, grid.xvx, grid.yvx, _vx_nodes(grid)), (vy, grid.xvy, grid.yvy, _vy_nodes(grid))
    )
        for (i, j) in nodes
            X, Y = _nd(xs[j], ys[i])
            X^2 + Y^2 < 1 && (total += v[i, j]^2)
        end
    end
    return sqrt(total)
end

"""
    _condition_ratios(g, dt)

Return the condition ratios of the case after a solve: the largest and the smallest
`ETAP / ETAPHI` (`eta / zeta`, condition C1) on the interior P nodes, and the largest and
the smallest `ETA / (GGG dt)` on the basic nodes and `ETAP / (GGGP dt)` on the interior P
nodes (condition C2). A negative smallest ratio flags a negative shear modulus.

Parameters
----------
g : Grids
    Grid fields after the solve.
dt : Float64
    Time step [s].

Returns
-------
NamedTuple
    `eta_zeta`, `eta_zeta_min`, `elastic`, `elastic_min`, `elastic_p`, `elastic_p_min`,
    dimensionless.
"""
function _condition_ratios(g, dt)
    inner = (2:(size(g.ETAP, 1) - 1), 2:(size(g.ETAP, 2) - 1))
    ratio = g.ETAP[inner...] ./ g.ETAPHI[inner...]
    basic = g.ETA ./ (g.GGG .* dt)
    pnode = g.ETAP[inner...] ./ (g.GGGP[inner...] .* dt)
    return (;
        eta_zeta=maximum(ratio),
        eta_zeta_min=minimum(ratio),
        elastic=maximum(basic),
        elastic_min=minimum(basic),
        elastic_p=maximum(pnode),
        elastic_p_min=minimum(pnode),
    )
end

"""
    run_case(N; case, spherical_flow=true, darcy_elimination=true)

Run one solve of the case and compare the velocity with the exact flow of `case`; for case
`:C` also compare `S33`, `P_t`, `P_f`, `EII`, `SII` and the heating `HS` of
`compute_shear_heating!` with their exact fields.

Parameters
----------
N : Int
    Basic grid resolution in x and y, odd and at least `NX_MIN`.
case : Symbol
    `:A`, `:A0`, `:B` or `:C`.
spherical_flow : Bool, optional
    Value of `geometry.spherical_flow` (default `true`); `false` only for case `:A`.
darcy_elimination : Bool, optional
    Value of `solver.darcy_elimination` (default `true`).

Returns
-------
NamedTuple
    `case`, `spherical_flow`, `N`; `converged` (`plastic_converged`), `n_plastic_passes`,
    `dt` [s] and `dt_dev` (`abs(state.dt / dt - 1)` after the solve, condition C6);
    `n_yielding`: number of basic nodes that yield in the last plastic pass (`YNY5`, C6);
    the ratios of `_condition_ratios` (C1 and C2); `g_dev`: largest change of `gx` and
    `gy` by the solve [m/s^2] (C3); `rhof_max`: largest `abs(RHOFX)` and `abs(RHOFY)` after
    the solve [kg/m^3]; `anchor_dev` (`anchor_deviation`); copies of `vx`, `vy` [m/s],
    `qxD`, `qyD` (Darcy flux at the vx and vy nodes) [m/s], `pr`, `pf` [Pa], `S33` and
    its increment `DS33` [Pa];
    `grid` (`grid_of`); the fields of `velocity_errors`; `centre_node` and `centre_fit`:
    rotation rate at the centre from `centre_node_rotation` and `centre_fit_rotation` in
    units of `F0 R / eta` (`NaN` for cases `:B` and `:C`); `s33_l2` (`s33_error`, `NaN`
    for the azimuthal cases); `pt_l2`, `pf_l2`, `eii_l2`, `sii_l2`, `hs_l2`: L2 errors for
    `r < R` of `pnode_errors`, and `s33_out`, `pt_out`, `pf_out`: its errors for
    `r >= R` (`NaN` except for case `:C`); `dp_max` (`pressure_deviation`) [Pa];
    `v_max`: largest velocity component [m/s]; `norm` (`profile_norm`) [m/s].

Raises
------
ArgumentError
    If `make_config` rejects `N`, `case` is not in `CASES`, or `case` has no exact
    solution for `spherical_flow`.
"""
function run_case(
    N::Int; case::Symbol, spherical_flow::Bool=true, darcy_elimination::Bool=true
)
    _check_case(case; spherical_flow=spherical_flow)
    cfg = make_config(
        N; spherical_flow=spherical_flow, darcy_elimination=darcy_elimination, case=case
    )
    return mktempdir() do out
        state, coords, ws, tio, _ = _init_case(cfg, out, case)
        tio === nothing || close(tio)
        g = state.grids
        dt = state.dt
        gx_in, gy_in = copy(g.gx), copy(g.gy)
        res = with_logger(NullLogger()) do
            return Erebus.solve_stokes_darcy!(
                state, coords, cfg, ws; titer=1, dt_step_initial=dt
            )
        end
        grid = grid_of(coords)
        vx, vy = copy(g.vx), copy(g.vy)
        errs = velocity_errors(vx, vy, grid, case; spherical_flow=spherical_flow)
        azimuthal = case in (:A, :A0)
        v = velocity_scale(case)
        c = case == :C ? _case_c_errors(g, grid) : _no_case_c_errors()
        return (;
            case,
            spherical_flow,
            N,
            converged=res.plastic_converged,
            n_plastic_passes=res.n_plastic_passes,
            dt,
            dt_dev=abs(state.dt / dt - 1),
            n_yielding=count(!iszero, g.YNY5),
            _condition_ratios(g, dt)...,
            g_dev=max(maximum(abs, g.gx .- gx_in), maximum(abs, g.gy .- gy_in)),
            rhof_max=max(maximum(abs, g.RHOFX), maximum(abs, g.RHOFY)),
            anchor_dev=anchor_deviation(g.pr, g.pf, grid, ambient_pressure(case)),
            vx,
            vy,
            qxD=copy(g.qxD),
            qyD=copy(g.qyD),
            pr=copy(g.pr),
            pf=copy(g.pf),
            S33=copy(g.S33),
            DS33=copy(g.DS33),
            grid,
            errs...,
            centre_node=azimuthal ? centre_node_rotation(vx, vy, grid; scale=v) : NaN,
            centre_fit=azimuthal ? centre_fit_rotation(vx, vy, grid; scale=v) : NaN,
            s33_l2=azimuthal ? NaN : s33_error(g.S33, grid; case=case),
            c...,
            dp_max=pressure_deviation(g.pr, grid),
            v_max=max(maximum(abs, vx), maximum(abs, vy)),
            norm=profile_norm(vx, vy, grid),
        )
    end
end

"""
    _case_c_errors(g, grid)

Return the errors of case `:C` at the P nodes after the solve (`pnode_errors`): `P_t` and
`P_f` against their exact values relative to `P_AMB_C`, `grids.EII` and `grids.SII` of the
stress update, and the heating `grids.HS` of `_case_c_heating!` (which overwrites
`g.HS`).

Parameters
----------
g : Grids
    Grid fields after the solve of case `:C`.
grid : NamedTuple
    Grid from `grid_of`.

Returns
-------
NamedTuple
    `pt_l2`, `pf_l2`, `eii_l2`, `sii_l2`, `hs_l2`, `s33_out`, `pt_out`, `pf_out`,
    dimensionless.
"""
function _case_c_errors(g, grid)
    ex = exact_pnode_fields(grid, :C)
    _case_c_heating!(g)
    pt = pnode_errors(g.pr, ex.pt, grid; offset=P_AMB_C)
    pf = pnode_errors(g.pf, ex.pf, grid; offset=P_AMB_C)
    return (;
        pt_l2=pt.l2_inside,
        pf_l2=pf.l2_inside,
        eii_l2=pnode_errors(g.EII, ex.eii, grid).l2_inside,
        sii_l2=pnode_errors(g.SII, ex.sii, grid).l2_inside,
        hs_l2=pnode_errors(g.HS, ex.hs, grid).l2_inside,
        s33_out=pnode_errors(g.S33, ex.s33, grid).out,
        pt_out=pt.out,
        pf_out=pf.out,
    )
end

"""
    _no_case_c_errors()

Return the fields of `_case_c_errors` set to `NaN`, for the cases other than `:C`.

Returns
-------
NamedTuple
"""
function _no_case_c_errors()
    return (;
        pt_l2=NaN,
        pf_l2=NaN,
        eii_l2=NaN,
        sii_l2=NaN,
        hs_l2=NaN,
        s33_out=NaN,
        pt_out=NaN,
        pf_out=NaN,
    )
end

"""
    elimination_difference(six, four)

Return the largest difference between the 6-variable and the 4-variable solve of one case:
`max(|vx6 - vx4| / max|vx4|, |vy6 - vy4| / max|vy4|, |pr6 - pr4| / stress_scale(case),
|pf6 - pf4| / stress_scale(case), |qxD6 - qxD4| / velocity_scale(case),
|qyD6 - qyD4| / velocity_scale(case))` over all nodes. The exact pressures of cases `:A`
and `:B` are uniform and their Darcy fluxes 0, so the pressure differences are scaled by
the stress unit and the flux differences by the velocity unit of the case.

Parameters
----------
six, four : NamedTuple
    `run_case` results with `darcy_elimination` off and on, or any `NamedTuple` with
    `case`, `vx`, `vy`, `pr`, `pf`, `qxD` and `qyD`.

Returns
-------
Float64
    Largest relative difference, dimensionless.

Raises
------
ArgumentError
    If the 2 results belong to different cases.
"""
function elimination_difference(six, four)
    six.case == four.case || throw(ArgumentError("the 2 runs must be of the same case"))
    rel(f, scale) = maximum(abs, getfield(six, f) .- getfield(four, f)) / scale
    return max(
        rel(:vx, maximum(abs, four.vx)),
        rel(:vy, maximum(abs, four.vy)),
        rel(:pr, stress_scale(four.case)),
        rel(:pf, stress_scale(four.case)),
        rel(:qxD, velocity_scale(four.case)),
        rel(:qyD, velocity_scale(four.case)),
    )
end

"""
    _order(coarse, fine)

Return the convergence order `log2(coarse / fine)` of an error that falls from `coarse`
to `fine` when the grid spacing halves.

Parameters
----------
coarse, fine : Float64
    Errors on the coarser and the finer grid, positive.

Returns
-------
Float64
"""
_order(coarse, fine) = log2(coarse / fine)

# Record rows of case C: metric name and the field of the `run_case` result.
const CASE_C_ROWS = (
    ("l2", :l2),
    ("l2_inside", :l2_inside),
    ("s33", :s33_l2),
    ("pt", :pt_l2),
    ("pf", :pf_l2),
    ("eii", :eii_l2),
    ("sii", :sii_l2),
    ("hs", :hs_l2),
)

"""
    metrics_from_runs(runs)

Return the validation metrics that the `run_case` results in `runs` determine.

The errors exist for each grid of `NX_LIST` present in `runs`. An order needs the
`NX_REFERENCE` and the finest grid of `NX_LIST`, a ratio the grid in both settings of case
`:A`, and an elimination row the run `(minimum(NX_LIST), case, true, :elim_off)` with its
base run. The centre errors are relative to the exact `Omega(0)`.

Parameters
----------
runs : AbstractDict
    Map from `(N, case, spherical_flow, variant)` to the `run_case` result, with `variant`
    `:base` or `:elim_off`.

Returns
-------
Dict{String,Float64}
    Metric values by id (see `EXPECTED_METRIC_IDS`); an id whose runs are missing is
    absent.
"""
function metrics_from_runs(runs::AbstractDict)::Dict{String,Float64}
    n_fine = maximum(NX_LIST)
    n_coarse = minimum(NX_LIST)
    metrics = Dict{String,Float64}()
    base(n, case, sf) = get(runs, (n, case, sf, :base), nothing)
    om0(sf) = omega_hat(:A, 0.0; spherical_flow=sf)
    for (sf, tag) in ((true, "a_on"), (false, "a_off"))
        for n in NX_LIST
            r = base(n, :A, sf)
            r === nothing && continue
            metrics["l2_$(tag)_nx$(n)"] = r.l2
            metrics["centre_fit_$(tag)_nx$(n)"] = abs(r.centre_fit / om0(sf) - 1)
            sf && (metrics["centre_node_$(tag)_nx$(n)"] = abs(r.centre_node / om0(sf) - 1))
        end
        mid, fine = base(NX_REFERENCE, :A, sf), base(n_fine, :A, sf)
        if mid !== nothing && fine !== nothing
            metrics["l2_order_$(tag)_$(NX_REFERENCE)_$(n_fine)"] = _order(mid.l2, fine.l2)
            metrics["centre_fit_order_$(tag)_$(NX_REFERENCE)_$(n_fine)"] = _order(
                abs(mid.centre_fit / om0(sf) - 1), abs(fine.centre_fit / om0(sf) - 1)
            )
        end
    end
    for n in NX_LIST
        r = base(n, :B, true)
        r === nothing && continue
        metrics["l2_b_nx$(n)"] = r.l2
        metrics["l2_inside_b_nx$(n)"] = r.l2_inside
        metrics["s33_b_nx$(n)"] = r.s33_l2
    end
    mid, fine = base(NX_REFERENCE, :B, true), base(n_fine, :B, true)
    if mid !== nothing && fine !== nothing
        for (name, f) in (("l2", :l2), ("l2_inside", :l2_inside), ("s33", :s33_l2))
            metrics["$(name)_order_b_$(NX_REFERENCE)_$(n_fine)"] = _order(
                getfield(mid, f), getfield(fine, f)
            )
        end
    end
    for n in (NX_REFERENCE, n_fine)
        on, off = base(n, :A, true), base(n, :A, false)
        if on !== nothing && off !== nothing
            metrics["centre_fit_ratio_nx$(n)"] = on.centre_fit / off.centre_fit
            metrics["profile_ratio_nx$(n)"] = on.norm / off.norm
        end
    end
    for n in NX_LIST
        r = base(n, :C, true)
        r === nothing && continue
        for (name, f) in CASE_C_ROWS
            metrics["$(name)_c_nx$(n)"] = getfield(r, f)
        end
    end
    mid, fine = base(NX_REFERENCE, :C, true), base(n_fine, :C, true)
    if mid !== nothing && fine !== nothing
        for (name, f) in CASE_C_ROWS
            metrics["$(name)_order_c_$(NX_REFERENCE)_$(n_fine)"] = _order(
                getfield(mid, f), getfield(fine, f)
            )
        end
    end
    for (case, tag) in ((:A, "a"), (:B, "b"), (:C, "c"))
        four = base(n_coarse, case, true)
        six = get(runs, (n_coarse, case, true, :elim_off), nothing)
        if four !== nothing && six !== nothing
            metrics["elim_$(tag)_nx$(n_coarse)"] = elimination_difference(six, four)
        end
    end
    return metrics
end

"""
    evaluate_spherical_momentum_benchmark(; nx_list=NX_LIST)

Run case `:A` with `spherical_flow` on and off and cases `:B` and `:C` with it on at every
grid in `nx_list`, and at the coarsest grid of `NX_LIST` (when in `nx_list`) cases `:A`,
`:B` and `:C` with `darcy_elimination = false`, and return the runs and the validation
metrics they determine.

Parameters
----------
nx_list : Tuple{Vararg{Int}}, optional
    Grids to run, a subset of `NX_LIST` (default `NX_LIST`). With the default the result
    holds every metric of `EXPECTED_METRIC_IDS`.

Returns
-------
NamedTuple
    `runs`: `Dict` from `(N, case, spherical_flow, variant)` to the `run_case` result;
    `metrics`: `Dict{String,Float64}` from `metrics_from_runs`.

Raises
------
ArgumentError
    If a grid in `nx_list` is not in `NX_LIST`.
"""
function evaluate_spherical_momentum_benchmark(; nx_list=NX_LIST)
    all(n -> n isa Integer && n in NX_LIST, nx_list) ||
        throw(ArgumentError("nx_list must be a subset of $NX_LIST, got $nx_list"))
    runs = Dict{Tuple,Any}()
    for n in nx_list
        for (case, sf) in ((:A, true), (:A, false), (:B, true), (:C, true))
            runs[(n, case, sf, :base)] = run_case(n; case=case, spherical_flow=sf)
            n == maximum(NX_LIST) && GC.gc()
        end
        if n == minimum(NX_LIST)
            for case in (:A, :B, :C)
                runs[(n, case, true, :elim_off)] = run_case(
                    n; case=case, darcy_elimination=false
                )
            end
        end
    end
    return (; runs, metrics=metrics_from_runs(runs))
end

const HASHED_CASE_CONSTANTS = (
    :L_BOX,
    :R_PLANET,
    :X_CENTRE,
    :Y_CENTRE,
    :NX_REFERENCE,
    :NX_LIST,
    :ETA_ROCK,
    :ETA_RATIO,
    :ETA_MIN,
    :ETA_MAX,
    :PHI,
    :K_PERM,
    :MU_F,
    :RHO_SOLID,
    :RHO_FLUID,
    :RHOCP,
    :P_AMB,
    :T_SURF,
    :SHEAR_MODULUS,
    :STRENGTH_FILL,
    :DT_YR,
    :F0,
    :U0,
    :B_POWER,
    :C_POWER,
    :ETA_RATIO_C,
    :P_AMB_C,
    :FIT_RHO_MIN,
    :FIT_RHO_MAX,
    :FIT_DEGREE,
    :OMEGA_A_ON,
    :OMEGA_A_OFF,
    :FORCE_A,
    :CASES,
)

const BENCHMARK_CONFIG_SECTIONS = (
    :grid,
    :geometry,
    :time,
    :solver,
    :materials,
    :thermodynamics,
    :poroelasticity,
    :venting,
    :reaction,
    :disk,
)

const EXPECTED_METRIC_IDS = (
    "l2_a_on_nx33",
    "l2_a_on_nx65",
    "l2_a_on_nx129",
    "l2_a_off_nx33",
    "l2_a_off_nx65",
    "l2_a_off_nx129",
    "l2_b_nx33",
    "l2_b_nx65",
    "l2_b_nx129",
    "l2_inside_b_nx33",
    "l2_inside_b_nx65",
    "l2_inside_b_nx129",
    "s33_b_nx33",
    "s33_b_nx65",
    "s33_b_nx129",
    "centre_fit_a_on_nx33",
    "centre_fit_a_on_nx65",
    "centre_fit_a_on_nx129",
    "centre_fit_a_off_nx33",
    "centre_fit_a_off_nx65",
    "centre_fit_a_off_nx129",
    "centre_node_a_on_nx33",
    "centre_node_a_on_nx65",
    "centre_node_a_on_nx129",
    "l2_c_nx33",
    "l2_c_nx65",
    "l2_c_nx129",
    "l2_inside_c_nx33",
    "l2_inside_c_nx65",
    "l2_inside_c_nx129",
    "s33_c_nx33",
    "s33_c_nx65",
    "s33_c_nx129",
    "pt_c_nx33",
    "pt_c_nx65",
    "pt_c_nx129",
    "pf_c_nx33",
    "pf_c_nx65",
    "pf_c_nx129",
    "eii_c_nx33",
    "eii_c_nx65",
    "eii_c_nx129",
    "sii_c_nx33",
    "sii_c_nx65",
    "sii_c_nx129",
    "hs_c_nx33",
    "hs_c_nx65",
    "hs_c_nx129",
    "l2_order_a_on_65_129",
    "l2_order_a_off_65_129",
    "l2_order_b_65_129",
    "l2_inside_order_b_65_129",
    "s33_order_b_65_129",
    "centre_fit_order_a_on_65_129",
    "centre_fit_order_a_off_65_129",
    "l2_order_c_65_129",
    "l2_inside_order_c_65_129",
    "s33_order_c_65_129",
    "pt_order_c_65_129",
    "pf_order_c_65_129",
    "eii_order_c_65_129",
    "sii_order_c_65_129",
    "hs_order_c_65_129",
    "centre_fit_ratio_nx65",
    "centre_fit_ratio_nx129",
    "profile_ratio_nx65",
    "profile_ratio_nx129",
    "elim_a_nx33",
    "elim_b_nx33",
    "elim_c_nx33",
)

# Regression band of an error row by its grid: about 5 to 20 % of the recorded value. The
# P_t rows of case C have errors 4 times smaller than the others and half the band.
const ERROR_BAND_BY_GRID = (nx33=2.0e-3, nx65=5.0e-4, nx129=1.5e-4)
const ERROR_BAND_PT_C = (nx33=1.0e-3, nx65=2.5e-4, nx129=7.5e-5)

const PHYSICS_REGRESSION_BANDS = Dict{
    String,NamedTuple{(:atol, :rtol),Tuple{Float64,Float64}}
}(
    (
        mid => (;
            atol=getfield(
                startswith(mid, "pt_c_") ? ERROR_BAND_PT_C : ERROR_BAND_BY_GRID,
                Symbol(last(split(mid, "_"))),
            ),
            rtol=0.0,
        ) for mid in EXPECTED_METRIC_IDS if
        !occursin("_order_", mid) && !occursin("_ratio_", mid) && !startswith(mid, "elim")
    )...,
    (
        mid => (; atol=0.05, rtol=0.0) for
        mid in EXPECTED_METRIC_IDS if occursin("_order_", mid)
    )...,
    "centre_fit_ratio_nx65" => (; atol=2.0e-4, rtol=0.0),
    "centre_fit_ratio_nx129" => (; atol=1.0e-4, rtol=0.0),
    "profile_ratio_nx65" => (; atol=2.0e-4, rtol=0.0),
    "profile_ratio_nx129" => (; atol=1.0e-4, rtol=0.0),
    "elim_a_nx33" => (; atol=1.0e-7, rtol=0.0),
    "elim_b_nx33" => (; atol=1.0e-7, rtol=0.0),
    "elim_c_nx33" => (; atol=1.0e-7, rtol=0.0),
)

const REPLAY_RTOL = 1.0e-6
const REPLAY_ATOL_BOUNDS = Dict{String,Float64}(
    "elim_a_nx33" => 1.0e-7, "elim_b_nx33" => 1.0e-7, "elim_c_nx33" => 1.0e-7
)

"""
    expected_metric_references()

Return the reference value of each validation metric: exact error `0.0`, second-order
convergence `2.0`, the exact centre rotation ratio `2/3` of case `:A`
(`omega_hat(:A, 0; true) / omega_hat(:A, 0; false)`), the exact profile norm ratio
`profile_norm_ratio()`, and `0.0` for the bounds of the elimination rows.

Returns
-------
Dict{String,Float64}
    Mapping from metric id to its reference value.
"""
function expected_metric_references()::Dict{String,Float64}
    refs = Dict{String,Float64}()
    centre_ratio = Float64(OMEGA_A_ON[1] / OMEGA_A_OFF[1])
    for mid in EXPECTED_METRIC_IDS
        refs[mid] = if occursin("_order_", mid)
            2.0
        elseif startswith(mid, "centre_fit_ratio")
            centre_ratio
        elseif startswith(mid, "profile_ratio")
            profile_norm_ratio()
        else
            0.0
        end
    end
    return refs
end

"""
    canonical_config_text(cfg; cfg_c=make_config(cfg.grid.Nx; spherical_flow=true, case=:C))

Return the text that `canonical_config_sha256` hashes: for each section in
`BENCHMARK_CONFIG_SECTIONS` in that order, one line `section.field=value` per field of the
section of `cfg` with the fields sorted by name, then the same lines of `cfg_c` with the
prefix `case_c.`, followed by one line `case.NAME=value` for each module constant in
`HASHED_CASE_CONSTANTS` in that order (see `_canonical_field_repr`).

Parameters
----------
cfg : Union{SimulationConfig,NamedTuple}
    Simulation configuration (or section container) of the case, from `make_config`.
cfg_c : Union{SimulationConfig,NamedTuple}, optional
    Configuration of case `:C`, which sets its own ambient pressure and viscosity ratio
    (default: `make_config` of case `:C` at the grid of `cfg`).

Returns
-------
String
    The canonical text.
"""
function canonical_config_text(
    cfg; cfg_c=make_config(cfg.grid.Nx; spherical_flow=true, case=:C)
)::String
    io = IOBuffer()
    for (prefix, c) in (("", cfg), ("case_c.", cfg_c)), sec_sym in BENCHMARK_CONFIG_SECTIONS
        sec = getproperty(c, sec_sym)
        for f_sym in sort(collect(fieldnames(typeof(sec))))
            val = getfield(sec, f_sym)
            println(io, prefix, sec_sym, ".", f_sym, "=", _canonical_field_repr(val))
        end
    end
    for name in HASHED_CASE_CONSTANTS
        println(io, "case.", name, "=", _canonical_field_repr(getfield(@__MODULE__, name)))
    end
    return String(take!(io))
end

"""
    canonical_config_sha256(cfg; cfg_c=make_config(cfg.grid.Nx; spherical_flow=true, case=:C))

Return the SHA-256 digest of `canonical_config_text(cfg; cfg_c)`, which covers the
configuration sections of the case and of case `:C`, the grid list `NX_LIST` and every
constant that `set_fields!`, the closed forms and the error measures use
(`HASHED_CASE_CONSTANTS`).

Parameters
----------
cfg : Union{SimulationConfig,NamedTuple}
    Simulation configuration (or section container) of the case, from `make_config`.
cfg_c : Union{SimulationConfig,NamedTuple}, optional
    Configuration of case `:C` (default: `make_config` of case `:C` at the grid of `cfg`).

Returns
-------
String
    64-character lowercase hexadecimal SHA-256 digest.
"""
function canonical_config_sha256(
    cfg; cfg_c=make_config(cfg.grid.Nx; spherical_flow=true, case=:C)
)::String
    return bytes2hex(sha256(canonical_config_text(cfg; cfg_c=cfg_c)))
end

"""
    expected_config_summary(cfg)

Return the configuration summary stored in the `[config]` section of the validation record.

Parameters
----------
cfg : SimulationConfig
    Configuration of the case with `spherical_flow = true`, from `make_config`. The keys
    `eta_zeta_ratio_c` and `P_amb_c_Pa` come from the configuration of case `:C` at the
    same grid.

Returns
-------
Dict{String,Any}
    Physical and numerical parameters of the case, keyed as in `REQUIRED_CONFIG_KEYS`.
"""
function expected_config_summary(cfg::SimulationConfig)::Dict{String,Any}
    cfg_c = make_config(cfg.grid.Nx; spherical_flow=true, case=:C)
    return Dict{String,Any}(
        "Nx" => collect(NX_LIST),
        "R_planet_m" => cfg.geometry.rplanet,
        "L_box_m" => cfg.grid.xsize,
        "eta_Pa_s" => cfg.materials.etasolidm[1],
        "eta_zeta_ratio" => cfg.thermodynamics.phim0 / cfg.solver.etaphikoef,
        "phi" => cfg.thermodynamics.phim0,
        "k_m2" => cfg.materials.kphim0[1],
        "mu_f_Pa_s" => cfg.materials.etafluidm[1],
        "rho_kg_m3" => RHO_SOLID,
        "rho_f_kg_m3" => RHO_FLUID,
        "shear_modulus_Pa" => SHEAR_MODULUS,
        "F0_N_m3" => F0,
        "U0_m_s" => U0,
        "b_power" => B_POWER,
        "c_power" => C_POWER,
        "eta_zeta_ratio_c" => cfg_c.thermodynamics.phim0 / cfg_c.solver.etaphikoef,
        "P_amb_c_Pa" => cfg_c.geometry.psurface,
        "dt_yr" => cfg.time.dt_initial,
        "P_amb_Pa" => cfg.geometry.psurface,
        "fit_rho_min" => FIT_RHO_MIN,
        "fit_rho_max" => FIT_RHO_MAX,
        "fit_degree" => FIT_DEGREE,
        "spherical_metric" => cfg.geometry.spherical_metric,
        "darcy_elimination" => cfg.solver.darcy_elimination,
    )
end

"""
    _parse_record(record_path)

Parse the validation record `record_path`.

Parameters
----------
record_path : String
    Path to an existing TOML validation record file.

Returns
-------
Union{Dict{String,Any},String}
    The parsed record, or the message of the parser error when the file is not valid TOML.
"""
function _parse_record(record_path::String)
    return try
        TOML.parsefile(record_path)
    catch err
        err isa TOML.ParserError || rethrow()
        "Cannot parse validation record $record_path: $(sprint(showerror, err))"
    end
end

"""
    check_physics_regression_bands(metrics, record_path)

Check that each metric lies within its regression band (`PHYSICS_REGRESSION_BANDS`) around
the recorded value in `record_path`. A file that is not valid TOML returns the message of
the parser error.

Parameters
----------
metrics : Dict{String,Float64}
    Measured metric dictionary.
record_path : String
    Path to the validation record TOML file.

Returns
-------
Tuple{Bool,Vector{String}}
    Pass flag and list of regression-band failure messages.

Raises
------
DomainError
    If `isempty(record_path)`.
"""
function check_physics_regression_bands(
    metrics::Dict{String,Float64}, record_path::String
)::Tuple{Bool,Vector{String}}
    isempty(record_path) && throw(DomainError(record_path, "Record path must be non-empty"))
    isfile(record_path) || return (false, ["Validation record not found: $record_path"])
    raw = _parse_record(record_path)
    raw isa String && return (false, [raw])
    recorded_by_id = Dict{String,Float64}()
    for m in get(raw, "metrics", Any[])
        if m isa AbstractDict &&
            get(m, "id", nothing) isa AbstractString &&
            get(m, "value", nothing) isa Real
            recorded_by_id[String(m["id"])] = Float64(m["value"])
        end
    end
    errors = String[]
    for mid in sort!(collect(keys(PHYSICS_REGRESSION_BANDS)))
        if !haskey(recorded_by_id, mid)
            push!(
                errors,
                "Validation record missing metric '$mid' for physics regression band",
            )
        elseif !haskey(metrics, mid)
            push!(errors, "Measured metrics missing '$mid' for physics regression band")
        else
            band = PHYSICS_REGRESSION_BANDS[mid]
            if !isapprox(metrics[mid], recorded_by_id[mid]; atol=band.atol, rtol=band.rtol)
                push!(
                    errors,
                    "Metric '$mid' = $(metrics[mid]) outside physics regression band around recorded $(recorded_by_id[mid]) (atol=$(band.atol), rtol=$(band.rtol))",
                )
            end
        end
    end
    return (isempty(errors), errors)
end

"""
    _config_mismatches(raw, cfg)

Return the messages for the `[config]` entries of a record that are missing, differ from
`expected_config_summary(cfg)` or are not part of it; floats compare with `rtol = 1e-12`.

Parameters
----------
raw : Dict{String,Any}
    Parsed validation record.
cfg : SimulationConfig
    Configuration of the case, from `make_config`.

Returns
-------
Vector{String}
    One message per missing, differing or unexpected key.
"""
function _config_mismatches(raw, cfg)
    errors = String[]
    rec_cfg = get(raw, "config", Dict{String,Any}())
    exp_cfg = expected_config_summary(cfg)
    for k in sort!(collect(setdiff(keys(rec_cfg), keys(exp_cfg))))
        push!(errors, "Recorded config has unexpected key '$k'")
    end
    for (k, v_exp) in exp_cfg
        if !haskey(rec_cfg, k)
            push!(errors, "Recorded config missing key '$k'")
            continue
        end
        v_rec = rec_cfg[k]
        matches = if v_rec isa AbstractFloat || v_exp isa AbstractFloat
            (v_rec isa Real && v_exp isa Real) &&
                isapprox(Float64(v_rec), Float64(v_exp); rtol=1.0e-12, atol=0.0)
        else
            v_rec == v_exp
        end
        matches ||
            push!(errors, "Recorded config.$k ($v_rec) does not match expected ($v_exp)")
    end
    return errors
end

"""
    _metric_mismatches(bench, raw_metrics, run_replay)

Return the messages for the recorded metrics that the measured metrics of `bench` fail:
a missing measurement, a reference other than `expected_metric_references`, a measured
value outside the tolerance of the record and, with `run_replay`, a measured value that
differs from the recorded value beyond `REPLAY_RTOL` and `REPLAY_ATOL_BOUNDS`.

Parameters
----------
bench : NamedTuple
    Result of `evaluate_spherical_momentum_benchmark`.
raw_metrics : Vector
    The `[[metrics]]` tables of a record that passed `validate_record_dict`.
run_replay : Bool
    Whether to compare the measured with the recorded values.

Returns
-------
Vector{String}
    Diagnostic messages.
"""
function _metric_mismatches(bench, raw_metrics, run_replay::Bool)
    errors = String[]
    exp_refs = expected_metric_references()
    for m in raw_metrics
        mid = String(m["id"])
        if !haskey(bench.metrics, mid)
            push!(errors, "Benchmark output missing metric '$mid'")
            continue
        end
        measured = bench.metrics[mid]
        recorded = Float64(m["value"])
        ref = Float64(m["reference"])
        atol = Float64(m["atol"])
        rtol = Float64(m["rtol"])
        comp = String(m["comparison"])
        if haskey(exp_refs, mid) &&
            !isapprox(ref, exp_refs[mid]; rtol=1.0e-12, atol=1.0e-12)
            push!(
                errors,
                "Recorded reference for '$mid' ($ref) does not match recomputed reference ($(exp_refs[mid]))",
            )
        end
        tol_ok = if comp == "isapprox"
            isapprox(measured, ref; atol=atol, rtol=rtol)
        else
            isfinite(measured) && 0.0 <= measured <= atol
        end
        tol_ok || push!(
            errors,
            "Measured metric '$mid' = $measured fails tolerance ($comp, ref=$ref, atol=$atol, rtol=$rtol)",
        )
        if run_replay
            replay_atol = get(REPLAY_ATOL_BOUNDS, mid, 1.0e-12)
            if !isapprox(measured, recorded; rtol=REPLAY_RTOL, atol=replay_atol)
                push!(
                    errors,
                    "Measured metric '$mid' = $measured differs from recorded value $recorded (rtol=$REPLAY_RTOL, atol=$replay_atol)",
                )
            end
        end
    end
    return errors
end

"""
    check_validation_reproducibility(bench, record_path; check_recorded_values=true)

Compare freshly evaluated metrics, reference values, regression bands and the configuration
hash and summary with a committed validation record. The replay of recorded values
(`REPLAY_RTOL`, `REPLAY_ATOL_BOUNDS`) runs when `check_recorded_values` is `true` and the
Julia major.minor version matches `metadata.julia_version`; the regression bands always run.
A record that fails `validate_record_dict` returns its schema errors only, and a file that is
not valid TOML the message of the parser error.

Parameters
----------
bench : NamedTuple
    Result of `evaluate_spherical_momentum_benchmark`.
record_path : String
    Path to the TOML validation record file.
check_recorded_values : Bool, optional
    Whether to replay the recorded values (default `true`).

Returns
-------
Tuple{Bool,Vector{String}}
    Pass flag and diagnostic error messages.

Raises
------
DomainError
    If `isempty(record_path)`.
"""
function check_validation_reproducibility(
    bench, record_path::String; check_recorded_values::Bool=true
)::Tuple{Bool,Vector{String}}
    isempty(record_path) && throw(DomainError(record_path, "Record path must be non-empty"))
    isfile(record_path) || return (false, ["Validation record not found: $record_path"])
    raw = _parse_record(record_path)
    raw isa String && return (false, [raw])
    ok_schema, errs = validate_record_dict(
        raw; require_doc_pages=false, check_doc_sync=false
    )
    # The checks below read the sections and metric fields that a schema error leaves invalid.
    ok_schema || return (false, errs)
    errors = copy(errs)
    cfg = make_config(NX_REFERENCE; spherical_flow=true)
    meta_sec = raw["metadata"]
    expected_sha = canonical_config_sha256(cfg)
    recorded_sha = string(get(meta_sec, "config_sha256", ""))
    if recorded_sha != expected_sha
        push!(
            errors,
            "Recorded config_sha256 ($recorded_sha) does not match canonical config SHA-256 ($expected_sha)",
        )
    end
    run_replay =
        check_recorded_values &&
        _same_julia_major_minor(string(get(meta_sec, "julia_version", "")))
    append!(errors, _config_mismatches(raw, cfg))
    raw_metrics = raw["metrics"]
    recorded_ids = Set{String}(String(m["id"]) for m in raw_metrics)
    if recorded_ids != Set{String}(EXPECTED_METRIC_IDS)
        push!(
            errors,
            "Recorded metric IDs ($(sort!(collect(recorded_ids)))) do not match expected metric IDs ($(sort!(collect(EXPECTED_METRIC_IDS))))",
        )
    end
    append!(errors, _metric_mismatches(bench, raw_metrics, run_replay))
    ok_bands, band_errs = check_physics_regression_bands(bench.metrics, record_path)
    ok_bands || append!(errors, band_errs)
    return (isempty(errors), errors)
end

end # module

"""
    SphericalPressureDiffusionValidation

Transient pore pressure diffusion in a rigid porous sphere whose surface vents to a fixed
pressure.

The matrix is rigid (shear and compaction viscosity `ETA_ROCK`, no elastic softening), so
the total pressure row gives `Pt - pr0 = kbw (Pf - pf0)` and the fluid pressure row reduces
to `S dPf/dt = (k / mu) lap(Pf)` with the storage coefficient
`S = phi bf + (1 - phi) bphi bs / (bphi + bs)`. The fluid pressure starts at
`P_AMB + P0` inside `R_PLANET` and the vent holds `Pf = P_AMB` at the surface cells. The
exact solution after `m` backward-Euler steps is a Fourier series in `x = r / R` for the
sphere (`geometry.spherical_flow = true`) and a Fourier-Bessel series for the disk
(`false`). The case builds the state with `init_simulation`, sets the grid fields directly
and calls `solve_stokes_darcy!` for `N_STEPS` steps on timestep 2, so the poroelastic
compressibilities are on.
"""
module SphericalPressureDiffusionValidation

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
const MU_F = 1.0e-3
const K_PERM = 1.0e-13
const PHI = 0.2
const BETA_PHI = 1.0e-10
const BETA_S = 1.0e-10
const BETA_S_ALT = 1.0e-11
const BETA_F = 4.5e-10
const P0 = 1.0e4
const P_AMB = 100.0
const T_SURF = 200.0
const K_VENT = 1.0e-11
const ETA_ROCK = 1.0e28
const ETA_MIN = 1.0e12
const SHEAR_MODULUS = 1.0e35
const STRENGTH_FILL = 1.0e20
const RX_AIR = 1.0e30
const RHO_SOLID = 3000.0
const RHOCP = 3.0e6
const RHO_FLUID = 1000.0
const TAU_END = 0.1
const N_STEPS = 20
const DTAU = TAU_END / N_STEPS
const M_RECORD = (5, 20)
const INNER_FRACTION = 0.9
const TAIL_BOUND = 1.0e-14
const SPHERE_MODES = 4000
const DISK_MODES = 200
const NX_LIST = (33, 65, 129)
const NX_REFERENCE = 65
const NX_MIN = 17

"""
    _j0j1(x)

Return the Bessel functions of the first kind `J0(x)` and `J1(x)` from their integral forms
`J0(x) = (1 / pi) int_0^pi cos(x sin t) dt` and `J1(x) = (1 / pi) int_0^pi cos(t - x sin t)
dt`, by the midpoint rule on `ceil(abs(x)) + 64` intervals. The integrands are smooth and
periodic, so the rule converges faster than any power of the interval width.

Parameters
----------
x : Real
    Argument, dimensionless, finite.

Returns
-------
Tuple{Float64,Float64}
    `(J0(x), J1(x))`.

Raises
------
DomainError
    If `x` is not finite.
"""
function _j0j1(x::Real)
    isfinite(x) || throw(DomainError(x, "the argument must be finite"))
    xf = Float64(x)
    n = ceil(Int, abs(xf)) + 64
    s0 = s1 = 0.0
    for k in 0:(n - 1)
        t = pi * (k + 0.5) / n
        s0 += cos(xf * sin(t))
        s1 += cos(t - xf * sin(t))
    end
    return (s0 / n, s1 / n)
end

"""
    _j0_zeros(n)

Return the first `n` positive zeros of `J0`, each by Newton iteration from
`(k - 1/4) pi` with the derivative `J0' = -J1`.

Parameters
----------
n : Int
    Number of zeros, at least 1.

Returns
-------
Vector{Float64}
    Zeros `j_1 < j_2 < ... < j_n`.

Raises
------
ArgumentError
    If `n < 1`.
ErrorException
    If a Newton iteration has not converged after 50 steps.
"""
function _j0_zeros(n::Int)
    n >= 1 || throw(ArgumentError("at least one zero is needed, got $n"))
    zeros_j0 = Vector{Float64}(undef, n)
    for k in 1:n
        z = (k - 0.25) * pi
        converged = false
        for _ in 1:50
            j0, j1 = _j0j1(z)
            step = j0 / j1
            z += step
            if abs(step) <= 4 * eps(z)
                converged = true
                break
            end
        end
        converged || error("the Newton iteration for zero $k of J0 did not converge")
        zeros_j0[k] = z
    end
    return zeros_j0
end

const J0_ZEROS = Tuple(_j0_zeros(DISK_MODES))
const J1_AT_ZEROS = Tuple(_j0j1(z)[2] for z in J0_ZEROS)

"""
    _check_series_args(x, d, m, dtau)

Check the arguments of the exact pressure series.

Parameters
----------
x : Real
    Nondimensional radius `r / R_PLANET`.
d : Real
    Flow dimension.
m : Real
    Number of time steps.
dtau : Real
    Nondimensional time step.

Returns
-------
Nothing

Raises
------
DomainError
    If `x` is not finite or lies outside `[0, 1]`, `d` is a `Bool` or neither 2 nor 3, `m`
    is a `Bool` or not an integer of at least 1, or `dtau` is not finite and positive.
"""
function _check_series_args(x::Real, d::Real, m::Real, dtau::Real)
    (isfinite(x) && 0 <= x <= 1) || throw(DomainError(x, "x must lie within [0, 1]"))
    (!(d isa Bool) && d in (2, 3)) ||
        throw(DomainError(d, "the flow dimension must be 2 or 3"))
    (m isa Integer && !(m isa Bool) && m >= 1) ||
        throw(DomainError(m, "the number of steps must be an integer of at least 1"))
    (isfinite(dtau) && dtau > 0) ||
        throw(DomainError(dtau, "the time step must be finite and positive"))
    return nothing
end

"""
    _mode_factor(lam, m, dtau, be)

Return the time factor of the mode with eigenvalue `lam` after `m` steps of length `dtau`:
`(1 + lam dtau)^(-m)` for backward Euler and `exp(-lam m dtau)` for the continuum.

Parameters
----------
lam : Float64
    Eigenvalue of the mode, dimensionless.
m : Int
    Number of time steps.
dtau : Float64
    Nondimensional time step.
be : Bool
    `true` for the backward-Euler factor, `false` for the continuum factor.

Returns
-------
Float64
    Time factor within `(0, 1]`.
"""
_mode_factor(lam, m, dtau, be) = be ? (1 + lam * dtau)^(-m) : exp(-lam * m * dtau)

"""
    _series_sum(x, d, m, dtau, be; tail=TAIL_BOUND)

Return the nondimensional pressure `p / P0` of the series at `x = r / R_PLANET`. The sphere
(`d = 3`) sums `2 (-1)^(n + 1) sin(n pi x) / (n pi x) f_n` with `lam_n = (n pi)^2` over at
most `SPHERE_MODES` modes; the disk (`d = 2`) sums `2 J0(j_n x) / (j_n J1(j_n)) f_n` with
`lam_n = j_n^2` over the `DISK_MODES` zeros of `J0_ZEROS`. The sum stops at the first mode
whose time factor `f_n` (`_mode_factor`) is below `tail`; with `tail = 0` it runs over the
whole table.

Parameters
----------
x : Real
    Nondimensional radius within `[0, 1]`.
d : Int
    Flow dimension, 2 or 3.
m : Int
    Number of time steps, at least 1.
dtau : Real
    Nondimensional time step, positive.
be : Bool
    `true` for the backward-Euler factors, `false` for the continuum.
tail : Float64, optional
    Stop bound on the time factor (default `TAIL_BOUND`).

Returns
-------
Float64
    `p / P0`.

Raises
------
DomainError
    From `_check_series_args`.
ErrorException
    If `tail > 0` and the table ends before a time factor falls below `tail`.
"""
function _series_sum(x::Real, d, m, dtau::Real, be::Bool; tail::Float64=TAIL_BOUND)
    _check_series_args(x, d, m, dtau)
    xf, mi, dtauf = Float64(x), Int(m), Float64(dtau)
    total = 0.0
    n_max = d == 3 ? SPHERE_MODES : DISK_MODES
    for n in 1:n_max
        z = d == 3 ? n * pi : J0_ZEROS[n]
        f = _mode_factor(z^2, mi, dtauf, be)
        f < tail && return total
        if d == 3
            shape = iszero(xf) ? 1.0 : sin(z * xf) / (z * xf)
            total += 2 * (isodd(n) ? 1 : -1) * shape * f
        else
            total += 2 * _j0j1(z * xf)[1] / (z * J1_AT_ZEROS[n]) * f
        end
    end
    tail > 0 && error(
        "the series for d = $d, m = $m needs more than $n_max modes to reach the tail bound $tail",
    )
    return total
end

"""
    n_modes(d, m, dtau, be)

Return the number of modes that `exact_p` sums: the modes before the first time factor
below `TAIL_BOUND`.

Parameters
----------
d : Int
    Flow dimension, 2 or 3.
m : Int
    Number of time steps, at least 1.
dtau : Real
    Nondimensional time step, positive.
be : Bool
    `true` for the backward-Euler factors, `false` for the continuum.

Returns
-------
Int
    Number of modes.

Raises
------
DomainError
    From `_check_series_args`.
ErrorException
    If the table ends before a time factor falls below `TAIL_BOUND`.
"""
function n_modes(d, m, dtau, be::Bool)
    _check_series_args(0.0, d, m, dtau)
    mi, dtauf = Int(m), Float64(dtau)
    n_max = d == 3 ? SPHERE_MODES : DISK_MODES
    for n in 1:n_max
        z = d == 3 ? n * pi : J0_ZEROS[n]
        _mode_factor(z^2, mi, dtauf, be) < TAIL_BOUND && return n - 1
    end
    return error("the series for d = $d, m = $m needs more than $n_max modes")
end

"""
    exact_p(x, d, m; dtau=DTAU, be=true)

Return the exact nondimensional fluid overpressure `(Pf - P_AMB) / P0` at `x = r / R` after
`m` time steps of the nondimensional length `dtau = D dt / R^2`, `D = k / (mu S)`, for the
initial value 1 inside and 0 on the surface (`_series_sum`). With `be = true` each mode
decays by the backward-Euler factor `(1 + lam dtau)^(-m)`, as the solver steps; with
`false` by the continuum factor `exp(-lam m dtau)`. At `x = 1` every mode vanishes and the
result is `0.0`, after the same argument and table checks.

Parameters
----------
x : Real
    Nondimensional radius within `[0, 1]`.
d : Int
    Flow dimension, 3 for the sphere (`geometry.spherical_flow = true`) and 2 for the disk.
m : Int
    Number of time steps, at least 1.
dtau : Real, optional
    Nondimensional time step, positive (default `DTAU`).
be : Bool, optional
    Backward-Euler (`true`, default) or continuum (`false`) time factors.

Returns
-------
Float64
    `(Pf - P_AMB) / P0`, dimensionless.

Raises
------
DomainError
    If `x` lies outside `[0, 1]`, `d` is a `Bool` or neither 2 nor 3, `m` is a `Bool` or
    not an integer of at least 1, or `dtau` is not finite and positive.
ErrorException
    If the series needs more modes than its table holds. At `dtau = DTAU` with backward
    Euler the disk table (`DISK_MODES`) needs `m >= 5` and the sphere table
    (`SPHERE_MODES`) `m >= 3`.
"""
function exact_p(x::Real, d, m; dtau::Real=DTAU, be::Bool=true)
    _check_series_args(x, d, m, dtau)
    if isone(x)
        n_modes(d, m, dtau, be)
        return 0.0
    end
    return _series_sum(x, d, m, dtau, be)
end

"""
    exact_p_cells(xs, d, m; dtau=DTAU)

Return the backward-Euler `exact_p` at each radius of `xs`, with one series evaluation per
distinct radius: each radius is rounded to 12 decimals and the series is evaluated at the
rounded radius, so radii that agree to 12 decimals share one value in any order. The cells
of the grid repeat their radii 8 times by symmetry, which makes the disk series affordable
on the finest grid. The arguments and every radius are checked before any evaluation, also
when `xs` is empty.

Parameters
----------
xs : AbstractVector{<:Real}
    Nondimensional radii within `[0, 1]`.
d : Int
    Flow dimension, 2 or 3.
m : Int
    Number of time steps, at least 1 (at `DTAU`, at least 5 for the disk and 3 for the
    sphere, see `exact_p`).
dtau : Real, optional
    Nondimensional time step, positive (default `DTAU`).

Returns
-------
Vector{Float64}
    `exact_p` at each rounded radius of `xs`.

Raises
------
DomainError
    If `d`, `m` or `dtau` fails `_check_series_args`, or a radius of `xs` lies outside
    `[0, 1]` or is not finite.
ErrorException
    If the series needs more modes than its table holds (`n_modes`).
"""
function exact_p_cells(xs::AbstractVector{<:Real}, d, m; dtau::Real=DTAU)
    _check_series_args(0.0, d, m, dtau)
    n_modes(d, m, dtau, true)
    for x in xs
        _check_series_args(x, d, m, dtau)
    end
    cache = Dict{Float64,Float64}()
    return map(xs) do x
        xr = round(Float64(x); digits=12) + 0.0
        return get!(() -> exact_p(xr, d, m; dtau=dtau), cache, xr)
    end
end

"""
    _check_compressibilities(beta_phi, beta_s)

Check the pore and solid compressibilities of the closed forms.

Parameters
----------
beta_phi : Real
    Pore compressibility [1/Pa].
beta_s : Real
    Solid compressibility [1/Pa].

Returns
-------
Nothing

Raises
------
DomainError
    If either value is negative or not finite, or both are 0.
"""
function _check_compressibilities(beta_phi::Real, beta_s::Real)
    for b in (beta_phi, beta_s)
        (isfinite(b) && b >= 0) ||
            throw(DomainError(b, "a compressibility must be finite and non-negative"))
    end
    beta_phi + beta_s > 0 ||
        throw(DomainError(beta_phi + beta_s, "beta_phi + beta_s must be positive"))
    return nothing
end

"""
    _check_porosity(phi)

Check the porosity of the closed forms.

Parameters
----------
phi : Real
    Porosity, dimensionless.

Returns
-------
Nothing

Raises
------
DomainError
    If `phi` is not finite or lies outside `(0, 1)`.
"""
function _check_porosity(phi::Real)
    (isfinite(phi) && 0 < phi < 1) ||
        throw(DomainError(phi, "the porosity must be finite and lie within (0, 1)"))
    return nothing
end

"""
    storage_coefficient(beta_phi, beta_s; phi=PHI, beta_f=BETA_F)

Return the storage coefficient `S = phi bf + (1 - phi) bphi bs / (bphi + bs)` of the fluid
pressure row of a rigid matrix. It equals `bd kbw (1 / ksk - kbw)` with the drained
compressibility `bd = (bphi + bs) / (1 - phi)`, `kbw = 1 - bs / bd` and
`1 / ksk = 1 + phi (bf - bs) / (bd - bs)`, after the total pressure row
`Pt - pr0 = kbw (Pf - pf0)` is put into the fluid pressure row.

Parameters
----------
beta_phi : Real
    Pore compressibility [1/Pa], non-negative.
beta_s : Real
    Solid compressibility [1/Pa], non-negative, `beta_phi + beta_s > 0`.
phi : Real, optional
    Porosity, within `(0, 1)` (default `PHI`).
beta_f : Real, optional
    Fluid compressibility [1/Pa], at least `beta_s` (default `BETA_F`).

Returns
-------
Float64
    Storage coefficient [1/Pa].

Raises
------
DomainError
    If a compressibility is negative or not finite, `beta_phi + beta_s = 0`, `phi` is not
    finite or lies outside `(0, 1)`, or `beta_f` is not finite or below `beta_s` (where the
    solver clamps the Skempton coefficient to 1, condition C5).
"""
function storage_coefficient(beta_phi::Real, beta_s::Real; phi=PHI, beta_f=BETA_F)
    _check_compressibilities(beta_phi, beta_s)
    _check_porosity(phi)
    (isfinite(beta_f) && beta_f >= beta_s) ||
        throw(DomainError(beta_f, "beta_f must be finite and at least beta_s = $beta_s"))
    bphi, bs, bf, ph = Float64.((beta_phi, beta_s, beta_f, phi))
    return ph * bf + (1 - ph) * bphi * bs / (bphi + bs)
end

"""
    biot_willis_closed_form(beta_phi, beta_s; phi=PHI)

Return the Biot-Willis coefficient `kbw = 1 - bs / bd = (bphi + phi bs) / (bphi + bs)`
from the case constants, without the solver helpers.

Parameters
----------
beta_phi : Real
    Pore compressibility [1/Pa], non-negative.
beta_s : Real
    Solid compressibility [1/Pa], non-negative, `beta_phi + beta_s > 0`.
phi : Real, optional
    Porosity, within `(0, 1)` (default `PHI`).

Returns
-------
Float64
    `kbw`, dimensionless, within `[phi, 1]`.

Raises
------
DomainError
    If a compressibility is negative or not finite, `beta_phi + beta_s = 0`, or `phi` is
    not finite or lies outside `(0, 1)`.
"""
function biot_willis_closed_form(beta_phi::Real, beta_s::Real; phi=PHI)
    _check_compressibilities(beta_phi, beta_s)
    _check_porosity(phi)
    bphi, bs, ph = Float64.((beta_phi, beta_s, phi))
    return (bphi + ph * bs) / (bphi + bs)
end

"""
    time_step_s(beta_s)

Return the time step `dt = DTAU R^2 mu S / k` [s] that reaches `tau = TAU_END` in `N_STEPS`
steps with the storage coefficient `S = storage_coefficient(BETA_PHI, beta_s)`.

Parameters
----------
beta_s : Real
    Solid compressibility [1/Pa], within `[0, BETA_F]`.

Returns
-------
Float64
    Time step [s].

Raises
------
DomainError
    If `beta_s` is not finite or lies outside `[0, BETA_F]`, with `beta_s` as the value.
"""
function time_step_s(beta_s::Real)
    (isfinite(beta_s) && 0 <= beta_s <= BETA_F) || throw(
        DomainError(beta_s, "beta_s must be finite and lie within [0, BETA_F = $BETA_F]"),
    )
    return DTAU * R_PLANET^2 * MU_F * storage_coefficient(BETA_PHI, beta_s) / K_PERM
end

"""
    make_config(N; spherical_flow, darcy_elimination=true, beta_s=BETA_S, beta_f=BETA_F)

Return the `SimulationConfig` of the case on an `N` x `N` grid, with the time step
`time_step_s(beta_s)` in years of `Erebus.SECONDS_PER_YEAR`.

Parameters
----------
N : Int
    Basic grid resolution in x and y, odd and at least `NX_MIN`, so that the centre lies on
    a basic node (see `centre_radius`).
spherical_flow : Bool
    Value of `geometry.spherical_flow`.
darcy_elimination : Bool, optional
    Value of `solver.darcy_elimination`: `true` solves the 4-variable system, `false` the
    6-variable system with the Darcy fluxes as unknowns (default `true`).
beta_s : Float64, optional
    Solid compressibility [1/Pa], within `(0, BETA_F]` (default `BETA_S`).
beta_f : Float64, optional
    Fluid compressibility `poroelasticity.betafluid` [1/Pa], finite and at least `beta_s`
    (default `BETA_F`). The time step does not depend on it.

Returns
-------
SimulationConfig

Raises
------
ArgumentError
    If `N < NX_MIN`, `N` is even, `beta_s` lies outside `(0, BETA_F]` or `beta_f` is not
    finite or below `beta_s`.
"""
function make_config(
    N::Int;
    spherical_flow::Bool,
    darcy_elimination::Bool=true,
    beta_s::Float64=BETA_S,
    beta_f::Float64=BETA_F,
)
    N >= NX_MIN ||
        throw(ArgumentError("N must be at least $NX_MIN to resolve the sphere, got $N"))
    isodd(N) || throw(
        ArgumentError("N must be odd so that the centre lies on a basic node, got $N")
    )
    0 < beta_s <= BETA_F ||
        throw(ArgumentError("beta_s must lie within (0, BETA_F], got $beta_s"))
    (isfinite(beta_f) && beta_f >= beta_s) ||
        throw(ArgumentError("beta_f must be finite and at least beta_s, got $beta_f"))
    eta3 = SVector{3,Float64}(ETA_ROCK, ETA_ROCK, ETA_ROCK)
    mu3 = SVector{3,Float64}(MU_F, MU_F, MU_F)
    return SimulationConfig(;
        grid=GridConfig(; xsize=L_BOX, ysize=L_BOX, Nx=N, Ny=N),
        geometry=GeometryConfig(;
            rplanet=R_PLANET,
            rcrust=R_PLANET,
            xcenter=X_CENTRE,
            ycenter=Y_CENTRE,
            psurface=P_AMB,
            spherical_metric=true,
            spherical_flow=spherical_flow,
            gravity_mode=:enclosed_mass,
        ),
        time=TimeConfig(;
            dt_initial=time_step_s(beta_s) / Erebus.SECONDS_PER_YEAR,
            dt_longest=1.0e12,
            n_steps=1,
        ),
        solver=SolverConfig(;
            darcy_elimination=darcy_elimination, etamax=1.0e30, etamin=ETA_MIN
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
            hydrofracture=false, betasolid=beta_s, betafluid=beta_f
        ),
        venting=VentingConfig(;
            active=true, mode=:darcy_sink, species=:H2O, k_vent=K_VENT, ice_sealing=false
        ),
        reaction=ReactionConfig(; active=false, fluid_overpressure_coupling=true),
        disk=DiskConfig(; enabled=true, model=:fixed, t_ambient=T_SURF, p_amb_disk=P_AMB),
        output=OutputConfig(; save_final=false, savematstep=100, visstep=100),
    )
end

"""
    _inside(coords, i, j)

Return whether the P node `(i, j)` lies inside the sphere, `r <= R_PLANET`.

Parameters
----------
coords : GridCoordinates
    Grid coordinates of the case.
i, j : Int
    Row and column index of the P node.

Returns
-------
Bool
    `true` when the node distance from the centre is at most `R_PLANET`.
"""
_inside(coords, i, j) = hypot(coords.xp[j] - X_CENTRE, coords.yp[i] - Y_CENTRE) <= R_PLANET

"""
    _n_outside_neighbours(coords, i, j)

Return the number of the 4 neighbours `(i +- 1, j)`, `(i, j +- 1)` of the P node `(i, j)`
that lie outside the sphere; a rock cell with at least one is a surface (vent) cell.

Parameters
----------
coords : GridCoordinates
    Grid coordinates of the case.
i, j : Int
    Row and column index of the P node, with `2 <= i <= coords.Ny` and
    `2 <= j <= coords.Nx`.

Returns
-------
Int
    Number of outside neighbours, from 0 to 4.
"""
function _n_outside_neighbours(coords, i, j)
    return count(
        !_inside(coords, a, b) for
        (a, b) in ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1))
    )
end

"""
    _surface_rock(coords, i, j)

Return whether the P node `(i, j)` is a surface rock cell: it lies inside the sphere and at
least one of its 4 neighbours lies outside. These are the cells whose faces to the outside
carry the vent.

Parameters
----------
coords : GridCoordinates
    Grid coordinates of the case.
i, j : Int
    Row and column index of the P node, with `2 <= i <= coords.Ny` and
    `2 <= j <= coords.Nx`.

Returns
-------
Bool
    `true` for a surface rock cell.
"""
function _surface_rock(coords, i, j)
    return _inside(coords, i, j) && _n_outside_neighbours(coords, i, j) > 0
end

"""
    set_fields!(state, coords)

Set the grid fields of the rigid porous sphere: no gravity, uniform porosity `PHI` and pore
compressibility `BETA_PHI`, total and fluid pressure `P_AMB + P0` inside the sphere and
`P_AMB` outside, impermeable faces outside the sphere and the face resistance
`MU_F / K_PERM` between 2 inside nodes. Every node has the viscosity `ETA_ROCK` (basic
nodes, `ETA`, `ETA0`, and P nodes, `ETAP`) and the shear modulus `SHEAR_MODULUS`, so the
viscoelastic viscosity `eta G dt / (G dt + eta)` equals `eta` to a relative `1e-6` or
better; the strengths `COH`, `TEN` are `STRENGTH_FILL`, with no friction.

Parameters
----------
state : SimulationState
    State from `init_simulation`; its grid fields are overwritten.
coords : GridCoordinates
    Staggered grid coordinates of `state`.

Returns
-------
Nothing
"""
function set_fields!(state, coords)
    g = state.grids
    Ny1, Nx1 = coords.Ny1, coords.Nx1
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
    fill!(g.BETAPHI, BETA_PHI)
    for f in (g.DQPF, g.DMP, g.gx, g.gy, g.SXX0, g.SXY0)
        fill!(f, 0.0)
    end
    for f in (g.RHO, g.RHOX, g.RHOY)
        fill!(f, RHO_SOLID)
    end
    fill!(g.RHOCP, RHOCP)
    fill!(g.RHOFX, RHO_FLUID)
    fill!(g.RHOFY, RHO_FLUID)
    for j in 1:Nx1, i in 1:Ny1
        inside = _inside(coords, i, j)
        p = inside ? P_AMB + P0 : P_AMB
        g.pr[i, j] = g.pr0[i, j] = g.pf[i, j] = g.pf0[i, j] = p
        g.RX[i, j] =
            (j < Nx1 && inside && _inside(coords, i, j + 1)) ? MU_F / K_PERM : RX_AIR
        g.RY[i, j] =
            (i < Ny1 && inside && _inside(coords, i + 1, j)) ? MU_F / K_PERM : RX_AIR
    end
    return nothing
end

"""
    _init_case(cfg, out)

Build the state of the case with `init_simulation`, set timestep 2 (so that the solve
keeps the compressibilities), the ambient pressure `P_AMB`, the time step of `cfg` in
seconds and the grid fields of `set_fields!`.

Parameters
----------
cfg : SimulationConfig
    Configuration from `make_config`.
out : String
    Output folder of `init_simulation`.

Returns
-------
Tuple
    The result of `init_simulation`: `(state, coords, ws, tio, rest)`.
"""
function _init_case(cfg, out)
    state, coords, ws, tio, rest = with_logger(NullLogger()) do
        return Erebus.init_simulation(cfg; output_path=out)
    end
    state.timestep = 2
    state.accumulators.P_amb = P_AMB
    state.dt = cfg.time.dt_initial * cfg.time.yearlength
    set_fields!(state, coords)
    return (state, coords, ws, tio, rest)
end

"""
    drain_cap_ratio(N; spherical_flow, beta_s=BETA_S, k_vent=K_VENT)

Return the largest vent drain of one time step relative to the fluid that a vent cell
holds, `(sum of C_face over the vent faces of the cell) P0 dt / (phi - phimin)`, over the
vent cells of the initial state. `C_face` is the face conductance of
`apply_venting_surface_boundary!`, read from its vent rate at the initial overpressure `P0`
with a time step short enough that its drain cap cannot act. The cap of the solve acts on
a cell when the ratio exceeds 1; the case keeps it at most 0.1.

Parameters
----------
N : Int
    Basic grid resolution in x and y, odd and at least `NX_MIN`.
spherical_flow : Bool
    Value of `geometry.spherical_flow` (face weights of the conductance).
beta_s : Float64, optional
    Solid compressibility [1/Pa]; sets the time step (default `BETA_S`).
k_vent : Float64, optional
    Vent permeability [m^2] (default `K_VENT`).

Returns
-------
NamedTuple
    `ratio`: largest drain ratio; `n_vent`: number of vent cells; `max_faces`: largest
    number of vent faces of one cell.

Raises
------
ArgumentError
    If `make_config` rejects the arguments.
"""
function drain_cap_ratio(
    N::Int; spherical_flow::Bool, beta_s::Float64=BETA_S, k_vent::Float64=K_VENT
)
    cfg = make_config(N; spherical_flow=spherical_flow, beta_s=beta_s)
    coords = Erebus.GridCoordinates(cfg)
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    pf = [_inside(coords, i, j) ? P_AMB + P0 : P_AMB for i in 1:Ny1, j in 1:Nx1]
    s_vent = zeros(Ny1, Nx1)
    Erebus.apply_venting_surface_boundary!(
        nothing,
        nothing,
        fill(T_SURF, Ny1, Nx1),
        coords,
        R_PLANET,
        X_CENTRE,
        Y_CENTRE,
        P_AMB;
        species=cfg.venting.species,
        k_vent=k_vent,
        conductance_factor=cfg.venting.conductance_factor,
        mode=cfg.venting.mode,
        pf=pf,
        PHI=fill(PHI, Ny1, Nx1),
        phimin=cfg.poroelasticity.phimin,
        dt=1.0e-12,
        eta_fluid_surf=cfg.materials.etafluidmm[2],
        L_sub=cfg.venting.L_sublimation,
        S_vent_out=s_vent,
        spherical_flow=spherical_flow,
    )
    dt = cfg.time.dt_initial * cfg.time.yearlength
    phi_avail = PHI - cfg.poroelasticity.phimin
    ratio = maximum(s_vent) * dt / phi_avail
    vent = [(i, j) for j in 2:(Nx1 - 1), i in 2:(Ny1 - 1) if s_vent[i, j] > 0.0]
    max_faces = maximum(_n_outside_neighbours(coords, i, j) for (i, j) in vent)
    return (; ratio, n_vent=length(vent), max_faces)
end

"""
    centre_radius(N)

Return the nondimensional radius `r / R_PLANET` of the P nodes nearest the centre on an
`N` x `N` grid: for an odd `N` the centre lies on a basic node, so the 4 nearest P nodes
sit at `dx / sqrt(2)` with `dx = L_BOX / (N - 1)`. For an even `N` a P node sits on the
centre, which the case does not use.

Parameters
----------
N : Int
    Basic grid resolution in x and y, odd and at least `NX_MIN`.

Returns
-------
Float64
    `r / R_PLANET` of the centre P nodes.

Raises
------
ArgumentError
    If `N < NX_MIN` or `N` is even.
"""
function centre_radius(N::Int)
    N >= NX_MIN || throw(ArgumentError("N must be at least $NX_MIN, got $N"))
    isodd(N) || throw(ArgumentError("N must be odd, got $N"))
    return L_BOX / (sqrt(2.0) * (N - 1) * R_PLANET)
end

"""
    _rock_cells(coords)

Return the interior P nodes `(i, j)` inside the sphere with their radius `r / R_PLANET`.

Parameters
----------
coords : GridCoordinates
    Staggered grid coordinates.

Returns
-------
Vector{Tuple{Int,Int,Float64}}
    `(i, j, x)` for every rock P node with `2 <= i <= coords.Ny`, `2 <= j <= coords.Nx`.
"""
function _rock_cells(coords)
    cells = Tuple{Int,Int,Float64}[]
    for j in 2:(coords.Nx1 - 1), i in 2:(coords.Ny1 - 1)
        r = hypot(coords.xp[j] - X_CENTRE, coords.yp[i] - Y_CENTRE)
        r <= R_PLANET && push!(cells, (i, j, r / R_PLANET))
    end
    return cells
end

"""
    _cell_weight(coords, i, j, spherical_flow)

Return the volume weight of the continuity row of P node `(i, j)`: `2 w_P` [m] from
`flow_metric_face_weights` with `spherical_flow = true`, else 1 (dimensionless). Times
`dx dy` it gives the cell volume [m^3] of the sphere or the cell area [m^2] of the disk.

Parameters
----------
coords : GridCoordinates
    Staggered grid coordinates.
i, j : Int
    Row and column index of the P node.
spherical_flow : Bool
    Value of `geometry.spherical_flow`.

Returns
-------
Float64
    Weight: `2 w_P` [m] with `spherical_flow = true`, `1.0` (dimensionless) with `false`.
"""
function _cell_weight(coords, i, j, spherical_flow)
    spherical_flow || return 1.0
    w = Erebus.flow_metric_face_weights(i, j, coords; xcenter=X_CENTRE, ycenter=Y_CENTRE)
    return 2.0 * w[5]
end

"""
    _c4_ratios(g, coords, dt, beta_s, timestep)

Return the ratios of the rigid-matrix condition C4 after a solve: the largest ratio of
the compaction rate coefficient `1 / (ETAPHI (1 - PHI))` to the storage rate coefficient
`bd / dt` on the interior P nodes, with the drained compressibility
`bd = (BETAPHI + bs) / (1 - PHI)` of the solved fields and the solid compressibility
`bs` of the solve (0 on timestep 1, where `BETAPHI` is also 0, so the ratio is `Inf`);
`ETA / (GGG dt)` on the basic nodes and `ETAP / (GGGP dt)` on the interior P nodes; and
the largest `1 / ETA`.

Parameters
----------
g : Grids
    Grid fields after the solve.
coords : GridCoordinates
    Staggered grid coordinates.
dt : Float64
    Time step [s].
beta_s : Float64
    Solid compressibility of the configuration [1/Pa].
timestep : Int
    Timestep of the solve; the solver uses no solid compressibility on timestep 1.

Returns
-------
NamedTuple
    `c4_compaction`, `c4_elastic`, `c4_elastic_p` (dimensionless) and `eta_inv_max`
    [1/(Pa s)].
"""
function _c4_ratios(g, coords, dt, beta_s, timestep)
    inner = (2:(coords.Ny1 - 1), 2:(coords.Nx1 - 1))
    bs = timestep == 1 ? 0.0 : beta_s
    phi = g.PHI[inner...]
    storage = @. (g.BETAPHI[inner...] + bs) / (1 - phi) / dt
    c4_compaction = maximum(@. inv(g.ETAPHI[inner...] * (1 - phi)) / storage)
    c4_elastic = maximum(g.ETA ./ (g.GGG .* dt))
    c4_elastic_p = maximum(g.ETAP[inner...] ./ (g.GGGP[inner...] .* dt))
    return (; c4_compaction, c4_elastic, c4_elastic_p, eta_inv_max=maximum(inv, g.ETA))
end

"""
    _centre_values(g, coords, d, m)

Return the computed and exact nondimensional fluid overpressure at the P node nearest the
centre after `m` steps, and the total overpressure there.

Parameters
----------
g : Grids
    Grid fields after step `m`.
coords : GridCoordinates
    Staggered grid coordinates.
d : Int
    Flow dimension, 3 for the sphere and 2 for the disk.
m : Int
    Number of steps taken.

Returns
-------
NamedTuple
    `p`: `(Pf - P_AMB) / P0`; `pt`: `(Pt - P_AMB) / P0`; `p_exact`: `exact_p`; `x`:
    `r / R_PLANET` of the node.
"""
function _centre_values(g, coords, d, m)
    cells = _rock_cells(coords)
    i, j, x = cells[argmin([c[3] for c in cells])]
    return (;
        p=(g.pf[i, j] - P_AMB) / P0,
        pt=(g.pr[i, j] - P_AMB) / P0,
        p_exact=exact_p(x, d, m),
        x,
    )
end

"""
    _l2_inner(g, coords, d, m)

Return the relative L2 error of the nondimensional fluid overpressure against the
backward-Euler `exact_p` after `m` steps over the rock P nodes with
`x = r / R_PLANET <= INNER_FRACTION`, with the weight `x` for the sphere and 1 for the disk:
`sqrt(sum w (p - p_exact)^2 / sum w p_exact^2)`.

Parameters
----------
g : Grids
    Grid fields after step `m`.
coords : GridCoordinates
    Staggered grid coordinates.
d : Int
    Flow dimension, 3 for the sphere and 2 for the disk.
m : Int
    Number of steps taken.

Returns
-------
Float64
    Relative L2 error, dimensionless.

Raises
------
ArgumentError
    If no rock P node lies within `INNER_FRACTION * R_PLANET`.
"""
function _l2_inner(g, coords, d, m)
    cells = filter(c -> c[3] <= INNER_FRACTION, _rock_cells(coords))
    isempty(cells) &&
        throw(ArgumentError("no P node lies within $(INNER_FRACTION) R of the centre"))
    exact = exact_p_cells([c[3] for c in cells], d, m)
    sq_err = sq_ref = 0.0
    for ((i, j, x), e) in zip(cells, exact)
        w = d == 3 ? x : 1.0
        sq_err += w * ((g.pf[i, j] - P_AMB) / P0 - e)^2
        sq_ref += w * e^2
    end
    return sqrt(sq_err / sq_ref)
end

"""
    _final_metrics(g, coords, spherical_flow, dt, beta_s)

Return the checks of step `N_STEPS`, before `pr0` and `pf0` take the solved pressures.

The total pressure check is the largest `abs((Pt - Pt_init) - kbw (Pf - Pf_init)) / P0`
over the rock P nodes, with `Pt_init = Pf_init = P_AMB + P0` and the closed-form
`biot_willis_closed_form`. The vented volume rate sums `S_vent_grid` and the storage
volume rate sums `S (Pf - pf0) / dt` with `S = storage_coefficient(BETA_PHI, beta_s)`, both
with the cell weights of the continuity rows.

Parameters
----------
g : Grids
    Grid fields after the last solve, with `pf0` from the step before.
coords : GridCoordinates
    Staggered grid coordinates.
spherical_flow : Bool
    Value of `geometry.spherical_flow`.
dt : Float64
    Time step [s].
beta_s : Float64
    Solid compressibility [1/Pa].

Returns
-------
NamedTuple
    `pt_kbw`; `vent_dev`: largest `abs(Pf - P_AMB) / P0` over the surface rock cells,
    `Inf` when one of them does not vent (`S_vent_grid <= 0`); `vent_flux`,
    `storage_flux`: vented and storage volume rate ([m^3/s] with `spherical_flow = true`;
    [m^2/s] with `false`); `vent_closure`: `abs(vent_flux / -storage_flux - 1)`; `n_rock`.
"""
function _final_metrics(g, coords, spherical_flow, dt, beta_s)
    kbw = biot_willis_closed_form(BETA_PHI, beta_s)
    s = storage_coefficient(BETA_PHI, beta_s)
    p_init = P_AMB + P0
    pt_kbw = vent_flux = storage_flux = 0.0
    vent_dev = -Inf
    vent_closed = false
    cells = _rock_cells(coords)
    for (i, j, _) in cells
        dev = abs((g.pr[i, j] - p_init) - kbw * (g.pf[i, j] - p_init)) / P0
        pt_kbw = max(pt_kbw, dev)
        if _surface_rock(coords, i, j)
            vent_closed |= !(g.S_vent_grid[i, j] > 0.0)
            vent_dev = max(vent_dev, abs(g.pf[i, j] - P_AMB) / P0)
        end
        vol = _cell_weight(coords, i, j, spherical_flow) * coords.dx * coords.dy
        vent_flux += g.S_vent_grid[i, j] * vol
        storage_flux += s * (g.pf[i, j] - g.pf0[i, j]) / dt * vol
    end
    return (;
        pt_kbw,
        vent_dev=(vent_closed || !isfinite(vent_dev)) ? Inf : vent_dev,
        vent_flux,
        storage_flux,
        vent_closure=abs(vent_flux / -storage_flux - 1),
        n_rock=length(cells),
    )
end

"""
    run_case(N; spherical_flow, darcy_elimination=true, beta_s=BETA_S)

Run the `N_STEPS` steps of the case and compare the fluid pressure with the backward-Euler
series. After each step `pr0` and `pf0` take the solved pressures; the time step stays
`time_step_s(beta_s)`.

Parameters
----------
N : Int
    Basic grid resolution in x and y, odd and at least `NX_MIN`.
spherical_flow : Bool
    Value of `geometry.spherical_flow`.
darcy_elimination : Bool, optional
    Value of `solver.darcy_elimination` (default `true`).
beta_s : Float64, optional
    Solid compressibility [1/Pa] (default `BETA_S`).

Returns
-------
NamedTuple
    `converged`: every solve returned `plastic_converged`; `timestep`: timestep of the
    solves; `dt`: time step [s]; `dt_dev`: largest `abs(state.dt / dt - 1)` after a solve;
    `xp`, `yp`: P node coordinates [m]; `pf`, `pr`, `qxD`, `qyD`: copies of the fields after
    the last step [Pa, m/s]; `pf_m5`: copy of `pf` after step 5 [Pa]; `v_max`: largest
    solid speed component over all steps [m/s]; the fields of `_c4_ratios` after the
    first solve; `p_centre_m<m>`,
    `p_exact_m<m>`, `centre_err_m<m>` (`p / p_exact - 1`) for `m` in `M_RECORD`;
    `pt_centre_m20`, `x_centre`, `l2_m20` (`_l2_inner`); and the fields of
    `_final_metrics`.

Raises
------
ArgumentError
    If `make_config` rejects the arguments.
"""
function run_case(
    N::Int; spherical_flow::Bool, darcy_elimination::Bool=true, beta_s::Float64=BETA_S
)
    cfg = make_config(
        N; spherical_flow=spherical_flow, darcy_elimination=darcy_elimination, beta_s=beta_s
    )
    d = spherical_flow ? 3 : 2
    return mktempdir() do out
        state, coords, ws, tio, _ = _init_case(cfg, out)
        tio === nothing || close(tio)
        g = state.grids
        dt = state.dt
        converged = true
        dt_dev = v_max = 0.0
        c4 = final = nothing
        centres = Dict{Int,Any}()
        pf_m5 = copy(g.pf)
        for m in 1:N_STEPS
            res = with_logger(NullLogger()) do
                return Erebus.solve_stokes_darcy!(
                    state, coords, cfg, ws; titer=1, dt_step_initial=dt
                )
            end
            converged &= res.plastic_converged
            dt_dev = max(dt_dev, abs(state.dt / dt - 1))
            v_max = max(v_max, maximum(abs, g.vx), maximum(abs, g.vy))
            m == 1 && (c4 = _c4_ratios(g, coords, dt, beta_s, state.timestep))
            m in M_RECORD && (centres[m] = _centre_values(g, coords, d, m))
            m == first(M_RECORD) && (pf_m5 .= g.pf)
            m == N_STEPS && (final = _final_metrics(g, coords, spherical_flow, dt, beta_s))
            g.pr0 .= g.pr
            g.pf0 .= g.pf
        end
        c5, c20 = centres[5], centres[20]
        return (;
            converged,
            timestep=state.timestep,
            dt,
            dt_dev,
            xp=copy(coords.xp),
            yp=copy(coords.yp),
            pf=copy(g.pf),
            pf_m5,
            pr=copy(g.pr),
            qxD=copy(g.qxD),
            qyD=copy(g.qyD),
            v_max,
            c4...,
            p_centre_m5=c5.p,
            p_exact_m5=c5.p_exact,
            centre_err_m5=c5.p / c5.p_exact - 1,
            p_centre_m20=c20.p,
            p_exact_m20=c20.p_exact,
            centre_err_m20=c20.p / c20.p_exact - 1,
            pt_centre_m20=c20.pt,
            x_centre=c20.x,
            l2_m20=_l2_inner(g, coords, d, N_STEPS),
            final...,
        )
    end
end

"""
    storage_defect_run(N; spherical_flow, factor)

Run the first `first(M_RECORD)` steps of the case with the storage coefficient `S` of the
fluid pressure row scaled by `factor` at the time step of the case. Only the fluid
compressibility moves, `beta_f = BETA_F + (factor - 1) S0 / PHI` with
`S0 = storage_coefficient(BETA_PHI, BETA_S)`, so the drained compressibility, the
Biot-Willis coefficient and the time step `time_step_s(BETA_S)` keep their values. The
series keeps `DTAU`, so the centre error after the steps measures the storage defect.

Parameters
----------
N : Int
    Basic grid resolution in x and y, odd and at least `NX_MIN`.
spherical_flow : Bool
    Value of `geometry.spherical_flow`.
factor : Real
    Ratio of the storage coefficient of the run to `S0`, finite and positive.

Returns
-------
NamedTuple
    `converged`: every solve returned `plastic_converged`; `centre_err_m5`:
    `p / p_exact - 1` at the P node nearest the centre after the steps; `pt_kbw`: the total
    pressure check of `_final_metrics` with the Biot-Willis coefficient of the case;
    `storage`: `S` of the run [1/Pa]; `beta_f` [1/Pa]; `dt`: time step [s].

Raises
------
DomainError
    If `factor` is not finite or not positive.
ArgumentError
    If `make_config` rejects `N` or `beta_f`; `beta_f < BETA_S` for
    `factor < 1 - PHI (BETA_F - BETA_S) / S0`.
"""
function storage_defect_run(N::Int; spherical_flow::Bool, factor::Real)
    (isfinite(factor) && factor > 0) ||
        throw(DomainError(factor, "the storage factor must be finite and positive"))
    s0 = storage_coefficient(BETA_PHI, BETA_S)
    beta_f = BETA_F + (Float64(factor) - 1) * s0 / PHI
    cfg = make_config(N; spherical_flow=spherical_flow, beta_f=beta_f)
    m = first(M_RECORD)
    return mktempdir() do out
        state, coords, ws, tio, _ = _init_case(cfg, out)
        tio === nothing || close(tio)
        g = state.grids
        dt = state.dt
        converged = true
        for _ in 1:m
            res = with_logger(NullLogger()) do
                return Erebus.solve_stokes_darcy!(
                    state, coords, cfg, ws; titer=1, dt_step_initial=dt
                )
            end
            converged &= res.plastic_converged
            g.pr0 .= g.pr
            g.pf0 .= g.pf
        end
        c = _centre_values(g, coords, spherical_flow ? 3 : 2, m)
        pt_kbw = _final_metrics(g, coords, spherical_flow, dt, BETA_S).pt_kbw
        return (;
            converged,
            centre_err_m5=c.p / c.p_exact - 1,
            pt_kbw,
            storage=storage_coefficient(BETA_PHI, BETA_S; beta_f=beta_f),
            beta_f,
            dt,
        )
    end
end

"""
    _rock_mask(run)

Return the mask of the P nodes inside the sphere for the coordinates of a `run_case` result.

Parameters
----------
run : NamedTuple
    Result of `run_case` with `xp` and `yp`.

Returns
-------
BitMatrix
    `true` at the rock P nodes.
"""
function _rock_mask(run)
    return [hypot(x - X_CENTRE, y - Y_CENTRE) <= R_PLANET for y in run.yp, x in run.xp]
end

"""
    elimination_difference(six, four)

Return the largest difference between the 6-variable and the 4-variable solve after the
last step: `max(|pf6 - pf4| / P0, |pr6 - pr4| / P0, |qxD6 - qxD4| / max|qxD4|,
|qyD6 - qyD4| / max|qyD4|)` over all nodes.

Parameters
----------
six, four : NamedTuple
    `run_case` results with `darcy_elimination` off and on.

Returns
-------
Float64
    Largest relative difference, dimensionless.
"""
function elimination_difference(six, four)
    rel(f, scale) = maximum(abs, getfield(six, f) .- getfield(four, f)) / scale
    return max(
        rel(:pf, P0),
        rel(:pr, P0),
        rel(:qxD, maximum(abs, four.qxD)),
        rel(:qyD, maximum(abs, four.qyD)),
    )
end

"""
    s_invariance(base, alt)

Return the largest difference of the nondimensional fluid overpressure `(Pf - P_AMB) / P0`
over the rock P nodes between 2 runs with different solid compressibility after the last
step. Both runs reach the same `tau`, so the difference vanishes when the solver diffuses
with the storage coefficient `storage_coefficient`.

Parameters
----------
base, alt : NamedTuple
    `run_case` results on the same grid.

Returns
-------
Float64
    Largest difference, dimensionless.
"""
function s_invariance(base, alt)
    mask = _rock_mask(base)
    return maximum(abs, (alt.pf .- base.pf)[mask]) / P0
end

const HASHED_CASE_CONSTANTS = (
    :L_BOX,
    :R_PLANET,
    :X_CENTRE,
    :Y_CENTRE,
    :NX_REFERENCE,
    :NX_LIST,
    :INNER_FRACTION,
    :MU_F,
    :K_PERM,
    :PHI,
    :BETA_PHI,
    :BETA_S,
    :BETA_S_ALT,
    :BETA_F,
    :P0,
    :P_AMB,
    :T_SURF,
    :K_VENT,
    :ETA_ROCK,
    :ETA_MIN,
    :SHEAR_MODULUS,
    :STRENGTH_FILL,
    :RX_AIR,
    :RHO_SOLID,
    :RHOCP,
    :RHO_FLUID,
    :TAU_END,
    :N_STEPS,
    :M_RECORD,
    :TAIL_BOUND,
    :SPHERE_MODES,
    :DISK_MODES,
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
    "centre_err_m5_off_nx33",
    "centre_err_m5_off_nx65",
    "centre_err_m5_off_nx129",
    "centre_err_m5_on_nx33",
    "centre_err_m5_on_nx65",
    "centre_err_m5_on_nx129",
    "centre_err_m20_off_nx33",
    "centre_err_m20_off_nx65",
    "centre_err_m20_off_nx129",
    "centre_err_m20_on_nx33",
    "centre_err_m20_on_nx65",
    "centre_err_m20_on_nx129",
    "l2_m20_off_nx33",
    "l2_m20_off_nx65",
    "l2_m20_off_nx129",
    "l2_m20_on_nx33",
    "l2_m20_on_nx65",
    "l2_m20_on_nx129",
    "centre_order_m20_off_65_129",
    "centre_order_m20_on_65_129",
    "l2_order_m20_off_65_129",
    "l2_order_m20_on_65_129",
    "centre_ratio_m20_nx65",
    "centre_ratio_m20_nx129",
    "pt_kbw_off_nx33",
    "pt_kbw_off_nx65",
    "pt_kbw_off_nx129",
    "pt_kbw_on_nx33",
    "pt_kbw_on_nx65",
    "pt_kbw_on_nx129",
    "s_invariance_off_nx33",
    "s_invariance_on_nx33",
    "elim_off_nx33",
    "elim_on_nx33",
    "vent_dev_off_nx65",
    "vent_dev_on_nx65",
    "vent_closure_off_nx65",
    "vent_closure_on_nx65",
    "v_max_off_nx65",
    "v_max_on_nx65",
)

const PHYSICS_REGRESSION_BANDS = Dict{
    String,NamedTuple{(:atol, :rtol),Tuple{Float64,Float64}}
}(
    (
        mid => (; atol=9.0e-5, rtol=0.0) for
        mid in EXPECTED_METRIC_IDS if startswith(mid, "centre_err_m5")
    )...,
    (
        mid => (; atol=0.002, rtol=0.0) for mid in EXPECTED_METRIC_IDS if
        startswith(mid, "centre_err_m20") || startswith(mid, "l2_m20")
    )...,
    (
        mid => (; atol=0.05, rtol=0.0) for
        mid in EXPECTED_METRIC_IDS if occursin("_order_", mid)
    )...,
    (
        mid => (; atol=0.002, rtol=0.0) for
        mid in EXPECTED_METRIC_IDS if startswith(mid, "centre_ratio")
    )...,
    (
        mid => (; atol=1.0e-9, rtol=0.0) for
        mid in EXPECTED_METRIC_IDS if startswith(mid, "pt_kbw")
    )...,
    "s_invariance_off_nx33" => (; atol=1.0e-9, rtol=0.0),
    "s_invariance_on_nx33" => (; atol=1.0e-9, rtol=0.0),
    "elim_off_nx33" => (; atol=1.0e-10, rtol=0.0),
    "elim_on_nx33" => (; atol=1.0e-10, rtol=0.0),
    "vent_dev_off_nx65" => (; atol=1.0e-4, rtol=0.0),
    "vent_dev_on_nx65" => (; atol=1.0e-4, rtol=0.0),
    "vent_closure_off_nx65" => (; atol=1.0e-8, rtol=0.0),
    "vent_closure_on_nx65" => (; atol=1.0e-8, rtol=0.0),
    "v_max_off_nx65" => (; atol=1.0e-15, rtol=0.0),
    "v_max_on_nx65" => (; atol=1.0e-15, rtol=0.0),
)

const REPLAY_RTOL = 1.0e-6
const REPLAY_ATOL_BOUNDS = Dict{String,Float64}(
    (mid => 1.0e-11 for mid in EXPECTED_METRIC_IDS if startswith(mid, "pt_kbw"))...,
    "s_invariance_off_nx33" => 1.0e-11,
    "s_invariance_on_nx33" => 1.0e-11,
    "elim_off_nx33" => 1.0e-11,
    "elim_on_nx33" => 1.0e-11,
    "vent_closure_off_nx65" => 1.0e-10,
    "vent_closure_on_nx65" => 1.0e-10,
    "v_max_off_nx65" => 1.0e-16,
    "v_max_on_nx65" => 1.0e-16,
)

"""
    metrics_from_runs(runs)

Return the validation metrics that the `run_case` results in `runs` determine.

The errors and the total pressure check exist for each grid of `NX_LIST` present in `runs`.
The orders of a setting need both the `NX_REFERENCE` and the finest grid of `NX_LIST`; the
vent and rigid-matrix rows need the `NX_REFERENCE` grid. A centre ratio needs the grid in
both settings. The solid compressibility and the elimination rows need the runs
`(33, spherical_flow, :beta_alt)` and `(33, spherical_flow, :elim_off)` with the base run.

Parameters
----------
runs : AbstractDict
    Map from `(N, spherical_flow, variant)` to the `run_case` result, with `variant` one of
    `:base`, `:elim_off` and `:beta_alt`.

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
    base(n, sf) = get(runs, (n, sf, :base), nothing)
    for (sf, tag) in ((false, "off"), (true, "on"))
        for n in NX_LIST
            r = base(n, sf)
            r === nothing && continue
            metrics["centre_err_m5_$(tag)_nx$(n)"] = abs(r.centre_err_m5)
            metrics["centre_err_m20_$(tag)_nx$(n)"] = abs(r.centre_err_m20)
            metrics["l2_m20_$(tag)_nx$(n)"] = r.l2_m20
            metrics["pt_kbw_$(tag)_nx$(n)"] = r.pt_kbw
        end
        mid, fine = base(NX_REFERENCE, sf), base(n_fine, sf)
        if mid !== nothing && fine !== nothing
            metrics["centre_order_m20_$(tag)_$(NX_REFERENCE)_$(n_fine)"] = log2(
                abs(mid.centre_err_m20) / abs(fine.centre_err_m20)
            )
            metrics["l2_order_m20_$(tag)_$(NX_REFERENCE)_$(n_fine)"] = log2(
                mid.l2_m20 / fine.l2_m20
            )
        end
        if mid !== nothing
            metrics["vent_dev_$(tag)_nx$(NX_REFERENCE)"] = mid.vent_dev
            metrics["vent_closure_$(tag)_nx$(NX_REFERENCE)"] = mid.vent_closure
            metrics["v_max_$(tag)_nx$(NX_REFERENCE)"] = mid.v_max
        end
        coarse = base(n_coarse, sf)
        if coarse !== nothing && haskey(runs, (n_coarse, sf, :beta_alt))
            metrics["s_invariance_$(tag)_nx$(n_coarse)"] = s_invariance(
                coarse, runs[(n_coarse, sf, :beta_alt)]
            )
        end
        if coarse !== nothing && haskey(runs, (n_coarse, sf, :elim_off))
            metrics["elim_$(tag)_nx$(n_coarse)"] = elimination_difference(
                runs[(n_coarse, sf, :elim_off)], coarse
            )
        end
    end
    for n in (NX_REFERENCE, n_fine)
        on, off = base(n, true), base(n, false)
        if on !== nothing && off !== nothing
            metrics["centre_ratio_m20_nx$(n)"] = on.p_centre_m20 / off.p_centre_m20
        end
    end
    return metrics
end

"""
    evaluate_spherical_pressure_diffusion_benchmark(; nx_list=NX_LIST)

Run the case at every grid in `nx_list` with `spherical_flow` off and on, and at the
coarsest grid of `NX_LIST` (when in `nx_list`) also with `darcy_elimination = false` and
with `beta_s = BETA_S_ALT`, and return the runs and the validation metrics they determine.

Parameters
----------
nx_list : Tuple{Vararg{Int}}, optional
    Grids to run, a subset of `NX_LIST` (default `NX_LIST`). With the default the result
    holds all 40 metrics.

Returns
-------
NamedTuple
    `runs`: `Dict` from `(N, spherical_flow, variant)` to the `run_case` result;
    `metrics`: `Dict{String,Float64}` from metric id (see `EXPECTED_METRIC_IDS`) to value,
    from `metrics_from_runs`.

Raises
------
ArgumentError
    If a grid in `nx_list` is not in `NX_LIST`.
"""
function evaluate_spherical_pressure_diffusion_benchmark(; nx_list=NX_LIST)
    all(n -> n isa Integer && n in NX_LIST, nx_list) ||
        throw(ArgumentError("nx_list must be a subset of $NX_LIST, got $nx_list"))
    runs = Dict{Tuple,Any}()
    for n in nx_list, sf in (false, true)
        runs[(n, sf, :base)] = run_case(n; spherical_flow=sf)
        if n == minimum(NX_LIST)
            runs[(n, sf, :elim_off)] = run_case(
                n; spherical_flow=sf, darcy_elimination=false
            )
            runs[(n, sf, :beta_alt)] = run_case(n; spherical_flow=sf, beta_s=BETA_S_ALT)
        end
        n == maximum(NX_LIST) && GC.gc()
    end
    return (; runs, metrics=metrics_from_runs(runs))
end

"""
    expected_metric_references()

Return the reference value of each validation metric: exact error `0.0`, first-order
convergence `1.0`, the centre ratio `exact_p(x_c, 3, 20) / exact_p(x_c, 2, 20)` at the
centre radius `x_c = centre_radius(N)` of its grid, and `0.0` for the bounds of the total
pressure, solid compressibility, elimination, vent and solid velocity checks.

Returns
-------
Dict{String,Float64}
    Mapping from metric id to its reference value.
"""
function expected_metric_references()::Dict{String,Float64}
    refs = Dict{String,Float64}()
    for mid in EXPECTED_METRIC_IDS
        refs[mid] = if occursin("_order_", mid)
            1.0
        elseif startswith(mid, "centre_ratio")
            x_c = centre_radius(parse(Int, last(split(mid, "_nx"))))
            exact_p(x_c, 3, N_STEPS) / exact_p(x_c, 2, N_STEPS)
        else
            0.0
        end
    end
    return refs
end

"""
    canonical_config_text(cfg)

Return the text that `canonical_config_sha256` hashes: for each section in
`BENCHMARK_CONFIG_SECTIONS` in that order, one line `section.field=value` per field of the
section with the fields sorted by name, followed by one line `case.NAME=value` for each
module constant in `HASHED_CASE_CONSTANTS` in that order (see `_canonical_field_repr`).

Parameters
----------
cfg : Union{SimulationConfig,NamedTuple}
    Simulation configuration (or section container) of the case, from `make_config`.

Returns
-------
String
    The canonical text.
"""
function canonical_config_text(cfg)::String
    io = IOBuffer()
    for sec_sym in BENCHMARK_CONFIG_SECTIONS
        sec = getproperty(cfg, sec_sym)
        for f_sym in sort(collect(fieldnames(typeof(sec))))
            val = getfield(sec, f_sym)
            println(io, sec_sym, ".", f_sym, "=", _canonical_field_repr(val))
        end
    end
    for name in HASHED_CASE_CONSTANTS
        println(io, "case.", name, "=", _canonical_field_repr(getfield(@__MODULE__, name)))
    end
    return String(take!(io))
end

"""
    canonical_config_sha256(cfg)

Return the SHA-256 digest of `canonical_config_text(cfg)`, which covers the configuration
sections of the case, the grid list `NX_LIST` and every constant that `set_fields!`, the
series and the error measures use (`HASHED_CASE_CONSTANTS`).

Parameters
----------
cfg : Union{SimulationConfig,NamedTuple}
    Simulation configuration (or section container) of the case, from `make_config`.

Returns
-------
String
    64-character lowercase hexadecimal SHA-256 digest.
"""
function canonical_config_sha256(cfg)::String
    return bytes2hex(sha256(canonical_config_text(cfg)))
end

"""
    expected_config_summary(cfg)

Return the configuration summary stored in the `[config]` section of the validation record.

Parameters
----------
cfg : SimulationConfig
    Configuration of the case with `spherical_flow = true`, from `make_config`.

Returns
-------
Dict{String,Any}
    Physical and numerical parameters of the case, keyed as in `REQUIRED_CONFIG_KEYS`.
"""
function expected_config_summary(cfg::SimulationConfig)::Dict{String,Any}
    pe = cfg.poroelasticity
    return Dict{String,Any}(
        "Nx" => collect(NX_LIST),
        "R_planet_m" => cfg.geometry.rplanet,
        "L_box_m" => cfg.grid.xsize,
        "mu_f_Pa_s" => cfg.materials.etafluidm[1],
        "k_m2" => cfg.materials.kphim0[1],
        "phi" => cfg.thermodynamics.phim0,
        "beta_phi_per_Pa" => BETA_PHI,
        "beta_s_per_Pa" => pe.betasolid,
        "beta_s_alt_per_Pa" => BETA_S_ALT,
        "beta_f_per_Pa" => pe.betafluid,
        "storage_per_Pa" => storage_coefficient(BETA_PHI, pe.betasolid),
        "kbw" => biot_willis_closed_form(BETA_PHI, pe.betasolid),
        "P0_Pa" => P0,
        "P_amb_Pa" => cfg.disk.p_amb_disk,
        "k_vent_m2" => cfg.venting.k_vent,
        "T_surf_K" => cfg.disk.t_ambient,
        "eta_rock_Pa_s" => cfg.materials.etasolidm[1],
        "shear_modulus_Pa" => SHEAR_MODULUS,
        "tau_end" => TAU_END,
        "n_steps" => N_STEPS,
        "dt_s" => cfg.time.dt_initial * cfg.time.yearlength,
        "inner_fraction" => INNER_FRACTION,
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
    Result of `evaluate_spherical_pressure_diffusion_benchmark`.
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
    Result of `evaluate_spherical_pressure_diffusion_benchmark`.
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

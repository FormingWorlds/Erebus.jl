"""
    SphericalCompactionValidation

Compaction of a porous sphere with a fixed effective pressure at its surface.

The rock has the shear viscosity `eta` and the compaction viscosity `zeta = ETAPHI`, with
`eta / zeta = ETA_RATIO` small. The total pressure `Pt` is then uniform inside, and the
effective pressure `Pe = Pt - Pf` obeys `delta^2 lap(Pe) = Pe` with the compaction length
`delta^2 = k zeta (1 - phi) / mu`. The surface keeps `Pe = PE_R`: the total pressure is the
anchor `psurface = P_AMB + PE_R` and the vent cells hold `Pf = P_vent = P_AMB`. The exact
solution is `Pe = PE_R R sinh(r / delta) / (r sinh(R / delta))` for the sphere
(`geometry.spherical_flow = true`) and `PE_R I0(r / delta) / I0(R / delta)` for the disk
(`false`). The solid moves against the fluid, `v_s = -q = -(k / mu) Pe'(r) r_hat`.
The case builds the state with `init_simulation`, sets the grid fields directly and calls
`solve_stokes_darcy!` once on timestep 1, with the poroelastic compressibilities off.
"""
module SphericalCompactionValidation

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
const DELTA = R_PLANET / 2
const ZETA = DELTA^2 * MU_F / (K_PERM * (1 - PHI))
const ETA_RATIO = 1.0e-4
const ETA_AIR = 1.0e12
const PE_R = 1.0e6
const P_AMB = 100.0
const T_SURF = 200.0
const K_VENT = 1.0e-6
const DT_YR = 1.0e-10
const RX_AIR = 1.0e30
const RHO_SOLID = 3000.0
const RHOCP = 3.0e6
const RHO_FLUID = 1000.0
const STRENGTH_FILL = 1.0e20
const INNER_FRACTION = 0.9
const SERIES_SWITCH = 1.0e-2
const BESSEL_X_MAX = 700.0
const NX_LIST = (33, 65, 129)
const NX_REFERENCE = 65
const NX_MIN = 17

"""
    _bessel_series(x, order)

Return the modified Bessel function `I_order(x)` for `order = 0` or `1` from its power
series `sum_m (x / 2)^(2m + order) / (m! (m + order)!)`, summed until a term no longer
changes the sum.

Parameters
----------
x : Real
    Argument, dimensionless, within `[0, BESSEL_X_MAX]`.
order : Int
    Order of the function, 0 or 1.

Returns
-------
Float64
    `I_order(x)`.

Raises
------
DomainError
    If `x` is not finite or lies outside `[0, BESSEL_X_MAX]`, or `order` is neither 0
    nor 1.
ErrorException
    If the series has not converged after 10 000 terms.
"""
function _bessel_series(x::Real, order::Int)
    (isfinite(x) && 0 <= x <= BESSEL_X_MAX) ||
        throw(DomainError(x, "the argument must lie within [0, $BESSEL_X_MAX]"))
    order in (0, 1) || throw(DomainError(order, "the order must be 0 or 1"))
    # The sum runs in Float64 for every real input, so a Float32 term cannot overflow.
    xf = Float64(x)
    q = (xf / 2)^2
    term = order == 0 ? 1.0 : xf / 2
    total = 0.0
    for m in 0:9_999
        total += term
        term <= eps(total) && return total
        term *= q / ((m + 1) * (m + 1 + order))
    end
    return error("the Bessel series did not converge for x = $x")
end

"""
    _i0(x)

Return the modified Bessel function of the first kind of order 0.

Parameters
----------
x : Real
    Argument, dimensionless, within `[0, BESSEL_X_MAX]`.

Returns
-------
Float64
    `I0(x)`.

Raises
------
DomainError
    If `x` is not finite or lies outside `[0, BESSEL_X_MAX]`.
"""
_i0(x::Real) = _bessel_series(x, 0)

"""
    _i1(x)

Return the modified Bessel function of the first kind of order 1, the derivative of `I0`.

Parameters
----------
x : Real
    Argument, dimensionless, within `[0, BESSEL_X_MAX]`.

Returns
-------
Float64
    `I1(x)`.

Raises
------
DomainError
    If `x` is not finite or lies outside `[0, BESSEL_X_MAX]`.
"""
_i1(x::Real) = _bessel_series(x, 1)

const CENTRE_RATIO_EXACT = 2 * _i0(R_PLANET / DELTA) / sinh(R_PLANET / DELTA)

"""
    _check_radius(r, d)

Check the arguments of the exact profile.

Parameters
----------
r : Real
    Distance from the centre [m].
d : Real
    Flow dimension.

Returns
-------
Nothing

Raises
------
DomainError
    If `r` is not finite or lies outside `[0, R_PLANET]`, or `d` is neither 2 nor 3.
"""
function _check_radius(r::Real, d::Real)
    (isfinite(r) && 0 <= r <= R_PLANET) ||
        throw(DomainError(r, "r must lie within [0, R_PLANET]"))
    d in (2, 3) || throw(DomainError(d, "the flow dimension must be 2 or 3"))
    return nothing
end

"""
    exact_pe(r, d)

Return the exact effective pressure `Pe = Pt - Pf` at distance `r` from the centre in the
limit `eta / zeta -> 0`, with `Pe = PE_R` at `r = R_PLANET`.

The sphere (`d = 3`) has `Pe = PE_R R sinh(x) / (r sinh(R / delta))` with `x = r / delta`;
for `x < SERIES_SWITCH` the factor `sinh(x) / x` is the series `1 + x^2 / 6 + x^4 / 120`.
The disk (`d = 2`) has `Pe = PE_R I0(x) / I0(R / delta)`.

Parameters
----------
r : Real
    Distance from the centre [m], within `[0, R_PLANET]`.
d : Real
    Flow dimension, 3 for the sphere (`geometry.spherical_flow = true`) and 2 for the disk.

Returns
-------
Float64
    Effective pressure [Pa].

Raises
------
DomainError
    If `r` lies outside `[0, R_PLANET]` or `d` is neither 2 nor 3.
"""
function exact_pe(r::Real, d::Real)
    _check_radius(r, d)
    x = Float64(r) / DELTA
    xr = R_PLANET / DELTA
    d == 2 && return PE_R * _i0(x) / _i0(xr)
    sinhc = x < SERIES_SWITCH ? 1 + x^2 / 6 + x^4 / 120 : sinh(x) / x
    return PE_R * xr * sinhc / sinh(xr)
end

"""
    exact_pe_derivative(r, d)

Return the radial derivative `Pe'(r)` of `exact_pe`.

The sphere has `Pe' = PE_R R (x cosh(x) - sinh(x)) / (r^2 sinh(R / delta))` with
`x = r / delta`; for `x < SERIES_SWITCH` the value is the series
`Pe(0) r / (3 delta^2) (1 + x^2 / 10)`, which avoids the cancellation of the closed form.
The disk has `Pe' = PE_R I1(x) / (delta I0(R / delta))`.

Parameters
----------
r : Real
    Distance from the centre [m], within `[0, R_PLANET]`.
d : Real
    Flow dimension, 3 for the sphere and 2 for the disk.

Returns
-------
Float64
    Derivative [Pa/m], zero at the centre.

Raises
------
DomainError
    If `r` lies outside `[0, R_PLANET]` or `d` is neither 2 nor 3.
"""
function exact_pe_derivative(r::Real, d::Real)
    _check_radius(r, d)
    rf = Float64(r)
    x = rf / DELTA
    xr = R_PLANET / DELTA
    d == 2 && return PE_R * _i1(x) / (DELTA * _i0(xr))
    pe_centre = PE_R * xr / sinh(xr)
    x < SERIES_SWITCH && return pe_centre * rf / (3 * DELTA^2) * (1 + x^2 / 10)
    return PE_R * R_PLANET * (x * cosh(x) - sinh(x)) / (rf^2 * sinh(xr))
end

"""
    make_config(N; spherical_flow, darcy_elimination=true, eta_ratio=ETA_RATIO, dt_yr=DT_YR)

Return the `SimulationConfig` of the case on an `N` x `N` grid.

The compaction viscosity is `ETAPHI = etaphikoef ETAP / phi`, so `etaphikoef = PHI /
eta_ratio` gives `ETAPHI = ZETA` for the rock viscosity `eta_ratio ZETA` that
`set_fields!` writes.

Parameters
----------
N : Int
    Basic grid resolution in x and y, at least `NX_MIN`.
spherical_flow : Bool
    Value of `geometry.spherical_flow`.
darcy_elimination : Bool, optional
    Value of `solver.darcy_elimination`: `true` solves the 4-variable system, `false` the
    6-variable system with the Darcy fluxes as unknowns (default `true`).
eta_ratio : Float64, optional
    Ratio of shear to compaction viscosity of the rock, within `(0, 1)` (default
    `ETA_RATIO`).
dt_yr : Float64, optional
    Time step [yr], positive (default `DT_YR`).

Returns
-------
SimulationConfig

Raises
------
ArgumentError
    If `N < NX_MIN`, `eta_ratio` lies outside `(0, 1)` or `dt_yr <= 0`.
"""
function make_config(
    N::Int;
    spherical_flow::Bool,
    darcy_elimination::Bool=true,
    eta_ratio::Float64=ETA_RATIO,
    dt_yr::Float64=DT_YR,
)
    N >= NX_MIN ||
        throw(ArgumentError("N must be at least $NX_MIN to resolve the sphere, got $N"))
    0 < eta_ratio < 1 ||
        throw(ArgumentError("eta_ratio must lie within (0, 1), got $eta_ratio"))
    dt_yr > 0 || throw(ArgumentError("dt_yr must be positive, got $dt_yr"))
    eta3 = SVector{3,Float64}(STRENGTH_FILL, STRENGTH_FILL, STRENGTH_FILL)
    mu3 = SVector{3,Float64}(MU_F, MU_F, MU_F)
    return SimulationConfig(;
        grid=GridConfig(; xsize=L_BOX, ysize=L_BOX, Nx=N, Ny=N),
        geometry=GeometryConfig(;
            rplanet=R_PLANET,
            rcrust=R_PLANET,
            xcenter=X_CENTRE,
            ycenter=Y_CENTRE,
            psurface=P_AMB + PE_R,
            spherical_metric=true,
            spherical_flow=spherical_flow,
            gravity_mode=:enclosed_mass,
        ),
        time=TimeConfig(; dt_initial=dt_yr, dt_longest=1.0e12, n_steps=1),
        solver=SolverConfig(;
            darcy_elimination=darcy_elimination,
            etamax=1.0e30,
            etamin=ETA_AIR,
            etaphikoef=PHI / eta_ratio,
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
    _surface_rock(coords, i, j)

Return whether the P node `(i, j)` is a surface rock cell: it lies inside the sphere and at
least one of its 4 neighbours `(i +- 1, j)`, `(i, j +- 1)` lies outside. These are the cells
whose faces to the outside carry the vent.

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
    _inside(coords, i, j) || return false
    return !(
        _inside(coords, i - 1, j) &&
        _inside(coords, i + 1, j) &&
        _inside(coords, i, j - 1) &&
        _inside(coords, i, j + 1)
    )
end

"""
    set_fields!(state, coords; eta_ratio=ETA_RATIO)

Set the grid fields of the compacting sphere: no gravity, uniform porosity, total and fluid
pressure `P_AMB + PE_R`, impermeable faces outside the sphere and the face resistance
`MU_F / K_PERM` between 2 inside nodes. Every basic node next to an inside P node has the
rock viscosity `eta_ratio ZETA`, so that each inside P node has `ETAP = eta_ratio ZETA` and
`ETAPHI = ZETA`; the other basic nodes have `ETA_AIR`. The shear moduli `GGG`, `GGGP` and
the strengths `COH`, `TEN` are `STRENGTH_FILL`, with no friction: the time step enters the
momentum rows through the viscoelastic viscosity `eta G dt / (G dt + eta)`, which differs from
`eta` by the relative `eta / (G dt + eta)`, and no node reaches its yield stress.

Parameters
----------
state : SimulationState
    State from `init_simulation`; its grid fields are overwritten.
coords : GridCoordinates
    Staggered grid coordinates of `state`.
eta_ratio : Float64, optional
    Ratio of shear to compaction viscosity of the rock (default `ETA_RATIO`).

Returns
-------
Nothing
"""
function set_fields!(state, coords; eta_ratio::Float64=ETA_RATIO)
    g = state.grids
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    eta_rock = eta_ratio * ZETA
    for j in axes(g.ETA, 2), i in axes(g.ETA, 1)
        rock = any(_inside(coords, i + di, j + dj) for di in 0:1, dj in 0:1)
        g.ETA[i, j] = rock ? eta_rock : ETA_AIR
        g.ETA0[i, j] = g.ETA[i, j]
    end
    fill!(g.ETAP, ETA_AIR)
    fill!(g.FRI, 0.0)
    for f in (g.COH, g.TEN, g.GGG, g.GGGP)
        fill!(f, STRENGTH_FILL)
    end
    fill!(g.tk1, T_SURF)
    fill!(g.tk2, T_SURF)
    for f in (g.PHI, g.PHIX, g.PHIY)
        fill!(f, PHI)
    end
    for f in (g.BETAPHI, g.DQPF, g.DMP, g.gx, g.gy, g.SXX0, g.SXY0)
        fill!(f, 0.0)
    end
    for f in (g.RHO, g.RHOX, g.RHOY)
        fill!(f, RHO_SOLID)
    end
    fill!(g.RHOCP, RHOCP)
    fill!(g.RHOFX, RHO_FLUID)
    fill!(g.RHOFY, RHO_FLUID)
    for f in (g.pr, g.pr0, g.pf, g.pf0)
        fill!(f, P_AMB + PE_R)
    end
    for j in 1:Nx1, i in 1:Ny1
        inside = _inside(coords, i, j)
        g.RX[i, j] =
            (j < Nx1 && inside && _inside(coords, i, j + 1)) ? MU_F / K_PERM : RX_AIR
        g.RY[i, j] =
            (i < Ny1 && inside && _inside(coords, i + 1, j)) ? MU_F / K_PERM : RX_AIR
    end
    return nothing
end

"""
    run_case(N; spherical_flow, darcy_elimination=true, eta_ratio=ETA_RATIO, dt_yr=DT_YR)

Solve the case once and compare the effective pressure and the solid velocity with the
exact solution.

Parameters
----------
N : Int
    Basic grid resolution in x and y, at least `NX_MIN`.
spherical_flow : Bool
    Value of `geometry.spherical_flow`.
darcy_elimination : Bool, optional
    Value of `solver.darcy_elimination` (default `true`).
eta_ratio : Float64, optional
    Ratio of shear to compaction viscosity of the rock (default `ETA_RATIO`).
dt_yr : Float64, optional
    Time step [yr] (default `DT_YR`).

Returns
-------
NamedTuple
    `converged`: `plastic_converged` flag returned by `solve_stokes_darcy!`;
    `xp`, `yp`: P node coordinates [m]; and the fields of `_metrics`.

Raises
------
ArgumentError
    If `make_config` rejects the arguments.
"""
function run_case(
    N::Int;
    spherical_flow::Bool,
    darcy_elimination::Bool=true,
    eta_ratio::Float64=ETA_RATIO,
    dt_yr::Float64=DT_YR,
)
    cfg = make_config(
        N;
        spherical_flow=spherical_flow,
        darcy_elimination=darcy_elimination,
        eta_ratio=eta_ratio,
        dt_yr=dt_yr,
    )
    return mktempdir() do out
        state, coords, ws, tio, _ = with_logger(NullLogger()) do
            return Erebus.init_simulation(cfg; output_path=out)
        end
        tio === nothing || close(tio)
        state.timestep = 1
        state.accumulators.P_amb = P_AMB
        set_fields!(state, coords; eta_ratio=eta_ratio)
        res = with_logger(NullLogger()) do
            return Erebus.solve_stokes_darcy!(
                state, coords, cfg, ws; titer=1, dt_step_initial=state.dt
            )
        end
        return (;
            converged=res.plastic_converged,
            xp=copy(coords.xp),
            yp=copy(coords.yp),
            _metrics(state, coords, spherical_flow)...,
        )
    end
end

"""
    _vs_l2_inner(g, coords, d)

Return the relative L2 vector error of the solid velocity against the exact velocity
`-(K_PERM / MU_F) Pe'(r) r_vec / r` over the faces with `0 < r <= INNER_FRACTION * R_PLANET`.

The `vx[i, j]` face lies at `(coords.xvx[j], coords.yvx[i])` and the `vy[i, j]` face at
`(coords.xvy[j], coords.yvy[i])`. The result is `sqrt(sum(err^2) / sum(exact^2))` over all
selected `vx` and `vy` faces.

Parameters
----------
g : Grids
    Solved grid fields with `vx` and `vy` [m/s].
coords : GridCoordinates
    Staggered grid coordinates.
d : Int
    Flow dimension, 3 for the sphere and 2 for the disk.

Returns
-------
Float64
    Relative L2 error, dimensionless.

Raises
------
ArgumentError
    If no face is selected within `INNER_FRACTION * R_PLANET` of the centre.
"""
function _vs_l2_inner(g, coords, d::Int)
    sq_err = sq_ref = 0.0
    n_faces = 0
    r_max = INNER_FRACTION * R_PLANET
    scale = -K_PERM / MU_F
    for j in 1:(coords.Nx), i in 2:(coords.Ny)
        dxf = coords.xvx[j] - X_CENTRE
        rf = hypot(dxf, coords.yvx[i] - Y_CENTRE)
        (0 < rf <= r_max) || continue
        exact = scale * exact_pe_derivative(rf, d) * dxf / rf
        sq_err += (g.vx[i, j] - exact)^2
        sq_ref += exact^2
        n_faces += 1
    end
    for j in 2:(coords.Nx), i in 1:(coords.Ny)
        dyf = coords.yvy[i] - Y_CENTRE
        rf = hypot(coords.xvy[j] - X_CENTRE, dyf)
        (0 < rf <= r_max) || continue
        exact = scale * exact_pe_derivative(rf, d) * dyf / rf
        sq_err += (g.vy[i, j] - exact)^2
        sq_ref += exact^2
        n_faces += 1
    end
    n_faces > 0 || throw(
        ArgumentError(
            "no velocity face is selected within $(INNER_FRACTION) R of the centre"
        ),
    )
    return sqrt(sq_err / sq_ref)
end

"""
    _metrics(state, coords, spherical_flow)

Return the error measures and the checks of the case conditions of a solved case.

The P nodes with `r <= R_PLANET` (the rock) enter, except the 2 outer rows and columns of
the grid. The compaction rate of a rock cell is `Pe / (ETAPHI (1 - PHI))` [1/s]; the vented
and the compaction volume rate sum the vent rate `S_vent_grid` and the compaction rate with
the cell weights of the continuity rows.

Parameters
----------
state : SimulationState
    Solved state; only `state.grids` is read.
coords : GridCoordinates
    Staggered grid coordinates.
spherical_flow : Bool
    Value of `geometry.spherical_flow`; sets the flow dimension and the cell weights.

Returns
-------
NamedTuple
    `pe`, `pf`, `vx`, `vy`, `qxD`, `qyD`: copies of `pr - pf` [Pa], `pf` [Pa], the solid
    velocity [m/s] and the Darcy flux [m/s];
    `pe_centre`, `pe_centre_exact`: computed and exact `Pe` at the P node nearest the centre
    [Pa]; `centre_err`: `pe_centre / pe_centre_exact - 1`;
    `l2_inner`: relative L2 error of `Pe` over the rock nodes with `r <= 0.9 R`;
    `vs_l2_inner`: solid velocity error from `_vs_l2_inner`;
    `pt_spread`: `(max(Pt) - min(Pt)) / PE_R` over the rock nodes;
    `pf_vent_dev`: largest `abs(Pf - P_AMB) / PE_R` over the surface rock cells
    (`_surface_rock`), `Inf` when one of them does not vent (`S_vent_grid <= 0`) or there is
    none;
    `vent_flux`, `comp_flux`: vented and compaction volume rate ([m^3/s] with
    `spherical_flow = true`; [m^2/s], per metre of depth, with `false`);
    `vent_closure`: `abs(vent_flux / comp_flux - 1)`;
    `n_rock`: number of rock nodes; `n_bad_etaphi`: number of rock nodes with `ETAPHI`
    not equal to `ZETA` within `rtol = 1e-12`.

Raises
------
ArgumentError
    If no P node lies inside the planet or none within `INNER_FRACTION * R_PLANET`, or if
    `_vs_l2_inner` finds no face.
"""
function _metrics(state, coords, spherical_flow)
    g = state.grids
    d = spherical_flow ? 3 : 2
    sq_err = sq_ref = 0.0
    vent_flux = comp_flux = 0.0
    n_rock = n_inner = n_bad_etaphi = 0
    pt_min, pt_max, pf_vent_dev = Inf, -Inf, -Inf
    vent_closed = false
    r_min, pe_c, pe_c_exact = Inf, NaN, NaN
    for j in 3:(coords.Nx1 - 2), i in 3:(coords.Ny1 - 2)
        r = hypot(coords.xp[j] - X_CENTRE, coords.yp[i] - Y_CENTRE)
        r <= R_PLANET || continue
        n_rock += 1
        pe = g.pr[i, j] - g.pf[i, j]
        pe_exact = exact_pe(r, d)
        if r <= INNER_FRACTION * R_PLANET
            sq_err += (pe - pe_exact)^2
            sq_ref += pe_exact^2
            n_inner += 1
        end
        if r < r_min
            r_min, pe_c, pe_c_exact = r, pe, pe_exact
        end
        pt_min = min(pt_min, g.pr[i, j])
        pt_max = max(pt_max, g.pr[i, j])
        isapprox(g.ETAPHI[i, j], ZETA; rtol=1.0e-12) || (n_bad_etaphi += 1)
        if _surface_rock(coords, i, j)
            vent_closed |= !(g.S_vent_grid[i, j] > 0.0)
            pf_vent_dev = max(pf_vent_dev, abs(g.pf[i, j] - P_AMB) / PE_R)
        end
        w = if spherical_flow
            2.0 * Erebus.flow_metric_face_weights(
                i, j, coords; xcenter=X_CENTRE, ycenter=Y_CENTRE
            )[5]
        else
            1.0
        end
        vent_flux += g.S_vent_grid[i, j] * w * coords.dx * coords.dy
        comp_flux += pe / (g.ETAPHI[i, j] * (1 - g.PHI[i, j])) * w * coords.dx * coords.dy
    end
    n_rock > 0 || throw(ArgumentError("no P node lies inside the planet"))
    n_inner > 0 ||
        throw(ArgumentError("no P node lies within $(INNER_FRACTION) R of the centre"))
    return (;
        pe=g.pr .- g.pf,
        pf=copy(g.pf),
        vx=copy(g.vx),
        vy=copy(g.vy),
        qxD=copy(g.qxD),
        qyD=copy(g.qyD),
        pe_centre=pe_c,
        pe_centre_exact=pe_c_exact,
        centre_err=pe_c / pe_c_exact - 1,
        l2_inner=sqrt(sq_err / sq_ref),
        vs_l2_inner=_vs_l2_inner(g, coords, d),
        pt_spread=(pt_max - pt_min) / PE_R,
        pf_vent_dev=(vent_closed || !isfinite(pf_vent_dev)) ? Inf : pf_vent_dev,
        vent_flux,
        comp_flux,
        vent_closure=abs(vent_flux / comp_flux - 1),
        n_rock,
        n_bad_etaphi,
    )
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
    :DELTA,
    :ZETA,
    :ETA_RATIO,
    :ETA_AIR,
    :PE_R,
    :P_AMB,
    :T_SURF,
    :K_VENT,
    :DT_YR,
    :RX_AIR,
    :RHO_SOLID,
    :RHOCP,
    :RHO_FLUID,
    :STRENGTH_FILL,
    :SERIES_SWITCH,
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
    "centre_err_off_nx33",
    "centre_err_off_nx65",
    "centre_err_off_nx129",
    "centre_err_on_nx33",
    "centre_err_on_nx65",
    "centre_err_on_nx129",
    "profile_l2_inner_off_nx33",
    "profile_l2_inner_off_nx65",
    "profile_l2_inner_off_nx129",
    "profile_l2_inner_on_nx33",
    "profile_l2_inner_on_nx65",
    "profile_l2_inner_on_nx129",
    "vs_l2_inner_off_nx33",
    "vs_l2_inner_off_nx65",
    "vs_l2_inner_off_nx129",
    "vs_l2_inner_on_nx33",
    "vs_l2_inner_on_nx65",
    "vs_l2_inner_on_nx129",
    "centre_order_off_65_129",
    "centre_order_on_65_129",
    "profile_order_off_65_129",
    "profile_order_on_65_129",
    "vs_order_off_65_129",
    "vs_order_on_65_129",
    "centre_ratio_on_off_nx65",
    "centre_ratio_on_off_nx129",
    "vent_closure_off_nx65",
    "vent_closure_on_nx65",
    "pf_vent_dev_off_nx65",
    "pf_vent_dev_on_nx65",
)

const PHYSICS_REGRESSION_BANDS = Dict{
    String,NamedTuple{(:atol, :rtol),Tuple{Float64,Float64}}
}(
    (
        mid => (; atol=0.002, rtol=0.0) for
        mid in EXPECTED_METRIC_IDS if startswith(mid, "centre_err") ||
            startswith(mid, "profile_l2") ||
            startswith(mid, "vs_l2")
    )...,
    (
        mid => (; atol=0.05, rtol=0.0) for
        mid in EXPECTED_METRIC_IDS if occursin("_order_", mid)
    )...,
    "centre_ratio_on_off_nx65" => (; atol=0.002, rtol=0.0),
    "centre_ratio_on_off_nx129" => (; atol=0.002, rtol=0.0),
    "vent_closure_off_nx65" => (; atol=1.0e-8, rtol=0.0),
    "vent_closure_on_nx65" => (; atol=1.0e-8, rtol=0.0),
    "pf_vent_dev_off_nx65" => (; atol=1.0e-6, rtol=0.0),
    "pf_vent_dev_on_nx65" => (; atol=1.0e-6, rtol=0.0),
)

const REPLAY_RTOL = 1.0e-6
const REPLAY_ATOL_BOUNDS = Dict{String,Float64}(
    "vent_closure_off_nx65" => 1.0e-10,
    "vent_closure_on_nx65" => 1.0e-10,
    "pf_vent_dev_off_nx65" => 1.0e-8,
    "pf_vent_dev_on_nx65" => 1.0e-8,
)

"""
    metrics_from_runs(runs)

Return the validation metrics that the `run_case` results in `runs` determine.

The centre, profile and velocity errors exist for each grid of `NX_LIST` present in `runs`.
The orders of a setting need both the `NX_REFERENCE` and the finest grid of `NX_LIST`, and
its vent closure and vent pressure need the `NX_REFERENCE` grid. A centre ratio needs the
grid in both settings, at `NX_REFERENCE` and at the finest grid.

Parameters
----------
runs : AbstractDict
    Map from `(N, spherical_flow)` to the `run_case` result.

Returns
-------
Dict{String,Float64}
    Metric values by id (see `EXPECTED_METRIC_IDS`); an id whose runs are missing is
    absent.
"""
function metrics_from_runs(runs::AbstractDict)::Dict{String,Float64}
    n_fine = maximum(NX_LIST)
    metrics = Dict{String,Float64}()
    names = (
        ("centre", "centre_err", :centre_err),
        ("profile", "profile_l2_inner", :l2_inner),
        ("vs", "vs_l2_inner", :vs_l2_inner),
    )
    for (sf, tag) in ((false, "off"), (true, "on"))
        for n in NX_LIST, (_, id, field) in names
            haskey(runs, (n, sf)) || continue
            metrics["$(id)_$(tag)_nx$(n)"] = getfield(runs[(n, sf)], field)
        end
        if haskey(runs, (NX_REFERENCE, sf)) && haskey(runs, (n_fine, sf))
            for (short, _, field) in names
                metrics["$(short)_order_$(tag)_$(NX_REFERENCE)_$(n_fine)"] = log2(
                    getfield(runs[(NX_REFERENCE, sf)], field) /
                    getfield(runs[(n_fine, sf)], field),
                )
            end
        end
        if haskey(runs, (NX_REFERENCE, sf))
            m = runs[(NX_REFERENCE, sf)]
            metrics["vent_closure_$(tag)_nx$(NX_REFERENCE)"] = m.vent_closure
            metrics["pf_vent_dev_$(tag)_nx$(NX_REFERENCE)"] = m.pf_vent_dev
        end
    end
    for n in (NX_REFERENCE, n_fine)
        if haskey(runs, (n, true)) && haskey(runs, (n, false))
            metrics["centre_ratio_on_off_nx$(n)"] =
                runs[(n, true)].pe_centre / runs[(n, false)].pe_centre
        end
    end
    return metrics
end

"""
    evaluate_spherical_compaction_benchmark(; nx_list=NX_LIST)

Run the case at every grid in `nx_list` with `spherical_flow` off and on and return the runs
and the validation metrics they determine.

Parameters
----------
nx_list : Tuple{Vararg{Int}}, optional
    Grids to run, a subset of `NX_LIST` (default `NX_LIST`). With the default the result
    holds all 30 metrics.

Returns
-------
NamedTuple
    `runs`: `Dict` from `(N, spherical_flow)` to the `run_case` result;
    `metrics`: `Dict{String,Float64}` from metric id (see `EXPECTED_METRIC_IDS`) to value,
    from `metrics_from_runs`.

Raises
------
ArgumentError
    If a grid in `nx_list` is not in `NX_LIST`.
"""
function evaluate_spherical_compaction_benchmark(; nx_list=NX_LIST)
    all(n -> n isa Integer && n in NX_LIST, nx_list) ||
        throw(ArgumentError("nx_list must be a subset of $NX_LIST, got $nx_list"))
    runs = Dict{Tuple,Any}()
    for n in nx_list, sf in (false, true)
        runs[(n, sf)] = run_case(n; spherical_flow=sf)
        n == maximum(NX_LIST) && GC.gc()
    end
    return (; runs, metrics=metrics_from_runs(runs))
end

"""
    expected_metric_references()

Return the reference value of each validation metric: exact centre, profile and velocity
error `0.0`, first-order convergence `1.0`, centre ratio `CENTRE_RATIO_EXACT` and `0.0` for
the vent closure and vent pressure bounds.

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
            CENTRE_RATIO_EXACT
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
sections of the case, the grid list `NX_LIST` and every constant that `set_fields!` and the
error measures use (`HASHED_CASE_CONSTANTS`).

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
    return Dict{String,Any}(
        "Nx" => collect(NX_LIST),
        "R_planet_m" => cfg.geometry.rplanet,
        "L_box_m" => cfg.grid.xsize,
        "mu_f_Pa_s" => cfg.materials.etafluidm[1],
        "k_m2" => cfg.materials.kphim0[1],
        "phi" => cfg.thermodynamics.phim0,
        "delta_m" => DELTA,
        "zeta_Pa_s" => ZETA,
        "eta_ratio" => cfg.thermodynamics.phim0 / cfg.solver.etaphikoef,
        "eta_air_Pa_s" => cfg.solver.etamin,
        "Pe_R_Pa" => cfg.geometry.psurface - cfg.disk.p_amb_disk,
        "k_vent_m2" => cfg.venting.k_vent,
        "T_surf_K" => cfg.disk.t_ambient,
        "P_amb_Pa" => cfg.disk.p_amb_disk,
        "dt_yr" => cfg.time.dt_initial,
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
    Result of `evaluate_spherical_compaction_benchmark`.
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
    meta_sec = get(raw, "metadata", Dict{String,Any}())
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

    rec_cfg = get(raw, "config", Dict{String,Any}())
    for (k, v_exp) in expected_config_summary(cfg)
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

    raw_metrics = get(raw, "metrics", Any[])
    recorded_ids = Set{String}(
        String(m["id"]) for m in raw_metrics if m isa AbstractDict && haskey(m, "id")
    )
    if recorded_ids != Set{String}(EXPECTED_METRIC_IDS)
        push!(
            errors,
            "Recorded metric IDs ($(sort!(collect(recorded_ids)))) do not match expected metric IDs ($(sort!(collect(EXPECTED_METRIC_IDS))))",
        )
    end

    exp_refs = expected_metric_references()
    for m in raw_metrics
        !(m isa AbstractDict && haskey(m, "id")) && continue
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
        elseif comp == "le_atol"
            isfinite(measured) && 0.0 <= measured <= atol
        else
            false
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

    ok_bands, band_errs = check_physics_regression_bands(bench.metrics, record_path)
    ok_bands || append!(errors, band_errs)
    return (ok_schema && isempty(errors), errors)
end

end # module

"""
    SphericalDarcyFlowValidation

Steady radial Darcy flow from a uniform fluid source in a rigid porous sphere.

The exact excess fluid pressure for a source `S` [1/s], permeability `k` and fluid
viscosity `mu`, with the surface held at the vent pressure (Dirichlet limit of the
venting condition), is `p - P_vent = mu S (R^2 - r^2) / (2 d k)`, with `d = 3` for the
sphere (`geometry.spherical_flow = true`) and `d = 2` for the disk (`false`).
The ice saturation pressure at `T_SURF` is below `P_AMB`, so `P_vent = P_AMB`.
The exact Darcy flux is `q = S r_vec / d`, with `r_vec` pointing from the centre.
The case builds the state with `init_simulation`, sets the grid fields directly and calls
`solve_stokes_darcy!` once on timestep 1, where the poroelastic compressibilities are off.

A second case, the permeability jump, gives the core `r < A_JUMP` the permeability `K_CORE`
and keeps `K_PERM` outside. The flux `q = S r_vec / d` is unchanged, and the excess pressure
is `mu S (R^2 - a^2) / (2 d K_PERM) + mu S (a^2 - r^2) / (2 d K_CORE)` for `r < a`, with
`a = A_JUMP`, and the uniform-permeability form for `a <= r <= R`.
"""
module SphericalDarcyFlowValidation

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
const K_CORE = K_PERM / 100
const A_JUMP = R_PLANET / 2
const S_SOURCE = 1.0e-11
const T_SURF = 200.0
const P_AMB = 100.0
const PHI = 0.2
const ETA_SOLID = 1.0e28
const RX_AIR = 1.0e30
const K_VENT = 1.0e-6
const PF_INIT = 1.0e5
const RHO_SOLID = 3000.0
const RHOCP = 3.0e6
const RHO_FLUID = 1000.0
const STRENGTH_FILL = 1.0e20
const DT_YR = 1.0e-10
const INNER_FRACTION = 0.9
const JUMP_EXCLUDE_CELLS = 2.0
const NX_LIST = (33, 65, 129)
const NX_REFERENCE = 65
const NX_MIN = 17

"""
    make_config(N; spherical_flow, darcy_elimination=true)

Return the `SimulationConfig` of the case on an `N` x `N` grid.

Parameters
----------
N : Int
    Basic grid resolution in x and y, at least `NX_MIN`.
spherical_flow : Bool
    Value of `geometry.spherical_flow`.
darcy_elimination : Bool, optional
    Value of `solver.darcy_elimination`: `true` solves the 4-variable system, `false` the
    6-variable system with the Darcy fluxes as unknowns (default `true`).

Returns
-------
SimulationConfig

Raises
------
ArgumentError
    If `N < NX_MIN`.
"""
function make_config(N::Int; spherical_flow::Bool, darcy_elimination::Bool=true)
    N >= NX_MIN ||
        throw(ArgumentError("N must be at least $NX_MIN to resolve the sphere, got $N"))
    eta3 = SVector{3,Float64}(ETA_SOLID, ETA_SOLID, ETA_SOLID)
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
        time=TimeConfig(; dt_initial=DT_YR, dt_longest=1.0e9, n_steps=1),
        solver=SolverConfig(;
            darcy_elimination=darcy_elimination, etamax=1.0e30, etamin=1.0e12
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
        poroelasticity=PoroelasticConfig(; hydrofracture=false),
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
    exact_excess_pressure(r, d; jump=false)

Return the exact excess fluid pressure `p - P_AMB` of the steady radial flow at distance `r`
from the centre, with the surface `r = R_PLANET` at the ambient pressure.

For uniform permeability `K_PERM` the value is `mu S (R^2 - r^2) / (2 d K_PERM)`. With
`jump = true` the core `r < A_JUMP` has the permeability `K_CORE` and the value there is
`mu S (R^2 - a^2) / (2 d K_PERM) + mu S (a^2 - r^2) / (2 d K_CORE)`, with `a = A_JUMP`.
The pressure is continuous at `a`; the flux `S r / d` is the same in both settings.

Parameters
----------
r : Real
    Distance from the centre [m], within `0 <= r <= R_PLANET`.
d : Real
    Flow dimension, 3 for the sphere (`geometry.spherical_flow = true`) and 2 for the disk.
jump : Bool, optional
    Whether the core has the permeability `K_CORE` (default `false`).

Returns
-------
Float64
    Excess pressure [Pa].

Raises
------
DomainError
    If `r` lies outside `[0, R_PLANET]` or `d` is neither 2 nor 3.
"""
function exact_excess_pressure(r::Real, d::Real; jump::Bool=false)
    0 <= r <= R_PLANET || throw(DomainError(r, "r must lie within [0, R_PLANET]"))
    d in (2, 3) || throw(DomainError(d, "the flow dimension must be 2 or 3"))
    if jump && r < A_JUMP
        return MU_F * S_SOURCE * (R_PLANET^2 - A_JUMP^2) / (2d * K_PERM) +
               MU_F * S_SOURCE * (A_JUMP^2 - r^2) / (2d * K_CORE)
    end
    return MU_F * S_SOURCE * (R_PLANET^2 - r^2) / (2d * K_PERM)
end

"""
    _core(coords, i, j)

Return whether the P node `(i, j)` lies in the core, `r < A_JUMP`.

Parameters
----------
coords : GridCoordinates
    Grid coordinates of the case.
i, j : Int
    Row and column index of the P node.

Returns
-------
Bool
    `true` when the node distance from the centre is below `A_JUMP`.
"""
_core(coords, i, j) = hypot(coords.xp[j] - X_CENTRE, coords.yp[i] - Y_CENTRE) < A_JUMP

"""
    set_fields!(state, coords; jump=false)

Set the grid fields of a rigid sphere with a uniform fluid source: impermeable faces outside
the sphere, no gravity, uniform porosity and fluid pressure. The node permeability is
`K_PERM`, or `K_CORE` for `r < A_JUMP` when `jump` is `true`. The resistance `mu / k` of a
face between 2 inside nodes is the mean of the 2 node values.

Parameters
----------
state : SimulationState
    State from `init_simulation`; its grid fields are overwritten.
coords : GridCoordinates
    Staggered grid coordinates of `state`.
jump : Bool, optional
    Whether the core `r < A_JUMP` has the permeability `K_CORE` (default `false`).

Returns
-------
Nothing
"""
function set_fields!(state, coords; jump::Bool=false)
    g = state.grids
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    for f in (g.ETA, g.ETA0, g.ETAP)
        fill!(f, ETA_SOLID)
    end
    fill!(g.FRI, 0.0)
    for f in (g.COH, g.TEN, g.GGG, g.GGGP)
        fill!(f, STRENGTH_FILL)
    end
    fill!(g.tk1, T_SURF)
    fill!(g.tk2, T_SURF)
    for f in (g.PHI, g.PHIX, g.PHIY)
        fill!(f, PHI)
    end
    for f in (g.RHO, g.RHOX, g.RHOY)
        fill!(f, RHO_SOLID)
    end
    fill!(g.RHOCP, RHOCP)
    fill!(g.RHOFX, RHO_FLUID)
    fill!(g.RHOFY, RHO_FLUID)
    fill!(g.pr, P_AMB)
    fill!(g.pr0, P_AMB)
    fill!(g.pf, PF_INIT)
    fill!(g.pf0, PF_INIT)
    fill!(g.gx, 0.0)
    fill!(g.gy, 0.0)
    node_k(i, j) = jump && _core(coords, i, j) ? K_CORE : K_PERM
    face_rx(i1, j1, i2, j2) = 0.5 * (MU_F / node_k(i1, j1) + MU_F / node_k(i2, j2))
    for j in 1:Nx1, i in 1:Ny1
        inside = _inside(coords, i, j)
        g.RX[i, j] = if (j < Nx1 && inside && _inside(coords, i, j + 1))
            face_rx(i, j, i, j + 1)
        else
            RX_AIR
        end
        g.RY[i, j] = if (i < Ny1 && inside && _inside(coords, i + 1, j))
            face_rx(i, j, i + 1, j)
        else
            RX_AIR
        end
        g.DQPF[i, j] = inside ? S_SOURCE : 0.0
    end
    return nothing
end

"""
    run_case(N; spherical_flow, darcy_elimination=true, jump=false)

Solve the case once and compare the fluid pressure and the Darcy flux with the exact
solution.

Parameters
----------
N : Int
    Basic grid resolution in x and y, at least `NX_MIN`.
spherical_flow : Bool
    Value of `geometry.spherical_flow`.
darcy_elimination : Bool, optional
    Value of `solver.darcy_elimination` (default `true`).
jump : Bool, optional
    Whether the core `r < A_JUMP` has the permeability `K_CORE` (default `false`).

Returns
-------
NamedTuple
    `converged`: `plastic_converged` flag returned by `solve_stokes_darcy!`;
    `pf`: fluid pressure grid [Pa];
    `qxD`, `qyD`: Darcy flux grids [m/s];
    `l2_inner`: relative L2 error of `p - P_AMB` over P nodes with `r <= 0.9 R`;
    `flux_l2_inner`: relative L2 vector error of the Darcy flux over the faces with
    `r <= 0.9 R`, without the faces near `A_JUMP` when `jump` is `true` (see
    `_flux_l2_inner`);
    `dp_centre`: computed `p - P_AMB` at the P node nearest the centre;
    `dp_centre_exact`: exact value at that node;
    `vent_flux`, `source_flux`: vented and produced fluid volume rate, summed with the
    cell weights of the continuity rows ([m^3/s] with `spherical_flow = true`; [m^2/s],
    that is [m^3/s] per metre of depth, with `false`);
    `max_solid_speed`: largest solid velocity component [m/s].

Raises
------
ArgumentError
    If `N < NX_MIN`.
"""
function run_case(
    N::Int; spherical_flow::Bool, darcy_elimination::Bool=true, jump::Bool=false
)
    cfg = make_config(N; spherical_flow=spherical_flow, darcy_elimination=darcy_elimination)
    return mktempdir() do out
        state, coords, ws, tio, _ = with_logger(NullLogger()) do
            return Erebus.init_simulation(cfg; output_path=out)
        end
        tio === nothing || close(tio)
        state.timestep = 1
        state.accumulators.P_amb = P_AMB
        set_fields!(state, coords; jump=jump)
        res = with_logger(NullLogger()) do
            return Erebus.solve_stokes_darcy!(
                state, coords, cfg, ws; titer=1, dt_step_initial=state.dt
            )
        end
        return (;
            converged=res.plastic_converged,
            _metrics(state, coords, spherical_flow; jump=jump)...,
        )
    end
end

"""
    _flux_l2_inner(g, coords, d; jump=false)

Return the relative L2 vector error of the Darcy flux against the exact flux
`S r_vec / d` over the faces with `r <= INNER_FRACTION * R_PLANET`. With `jump` the faces with
`|r - A_JUMP| <= JUMP_EXCLUDE_CELLS * dx` are left out.

At each face the error is the stored component minus the same component of the exact flux,
with `r_vec` pointing from the centre to the face. The `qxD[i, j]` face lies at
`(coords.xvx[j], coords.yvx[i])` and the `qyD[i, j]` face at `(coords.xvy[j], coords.yvy[i])`.
The result is `sqrt(sum(err^2) / sum(exact^2))` over all selected `qxD` and `qyD` faces.

Parameters
----------
g : Grids
    Solved grid fields with `qxD` and `qyD` [m/s].
coords : GridCoordinates
    Staggered grid coordinates.
d : Float64
    Flow dimension, 3 for the sphere and 2 for the disk.
jump : Bool, optional
    Whether to leave out the faces within `JUMP_EXCLUDE_CELLS * dx` of `r = A_JUMP`
    (default `false`).

Returns
-------
Float64
    Relative L2 error, dimensionless.

Raises
------
ArgumentError
    If no face is selected within `INNER_FRACTION * R_PLANET` of the centre, or the exact
    flux is zero on all selected faces.
"""
function _flux_l2_inner(g, coords, d::Float64; jump::Bool=false)
    sq_err = sq_ref = 0.0
    n_faces = 0
    r_max = INNER_FRACTION * R_PLANET
    r_skip = JUMP_EXCLUDE_CELLS * coords.dx
    for j in 1:(coords.Nx), i in 2:(coords.Ny)
        dxf = coords.xvx[j] - X_CENTRE
        rf = hypot(dxf, coords.yvx[i] - Y_CENTRE)
        rf <= r_max || continue
        jump && abs(rf - A_JUMP) <= r_skip && continue
        exact = S_SOURCE * dxf / d
        sq_err += (g.qxD[i, j] - exact)^2
        sq_ref += exact^2
        n_faces += 1
    end
    for j in 2:(coords.Nx), i in 1:(coords.Ny)
        dyf = coords.yvy[i] - Y_CENTRE
        rf = hypot(coords.xvy[j] - X_CENTRE, dyf)
        rf <= r_max || continue
        jump && abs(rf - A_JUMP) <= r_skip && continue
        exact = S_SOURCE * dyf / d
        sq_err += (g.qyD[i, j] - exact)^2
        sq_ref += exact^2
        n_faces += 1
    end
    n_faces > 0 || throw(
        ArgumentError("no flux face is selected within $(INNER_FRACTION) R of the centre"),
    )
    sq_ref > 0.0 || throw(ArgumentError("the exact flux is zero on every selected face"))
    return sqrt(sq_err / sq_ref)
end

"""
    _metrics(state, coords, spherical_flow; jump=false)

Return the error measures of a solved case: the pressure and flux errors against the exact
solution, the centre pressure, the vented and produced rates and the largest solid speed.

Parameters
----------
state : SimulationState
    Solved state; only `state.grids` is read.
coords : GridCoordinates
    Staggered grid coordinates.
spherical_flow : Bool
    Value of `geometry.spherical_flow`; sets the flow dimension and the cell weights.
jump : Bool, optional
    Whether the case has the permeability jump (default `false`); selects the exact
    profile and the faces of `_flux_l2_inner`.

Returns
-------
NamedTuple
    The fields of `run_case` from `pf` onward.

Raises
------
ArgumentError
    If no P node lies inside the planet or none within `INNER_FRACTION * R_PLANET`, or if
    `_flux_l2_inner` finds no face.
"""
function _metrics(state, coords, spherical_flow; jump::Bool=false)
    g = state.grids
    d = spherical_flow ? 3.0 : 2.0
    sq_err = sq_ref = 0.0
    vent_flux = source_flux = 0.0
    n_inside = n_inner = 0
    r_min, dp_c, dp_c_exact = Inf, NaN, NaN
    for j in 3:(coords.Nx1 - 2), i in 3:(coords.Ny1 - 2)
        r = hypot(coords.xp[j] - X_CENTRE, coords.yp[i] - Y_CENTRE)
        r <= R_PLANET || continue
        n_inside += 1
        dp = g.pf[i, j] - P_AMB
        dp_exact = exact_excess_pressure(r, d; jump=jump)
        if r <= INNER_FRACTION * R_PLANET
            sq_err += (dp - dp_exact)^2
            sq_ref += dp_exact^2
            n_inner += 1
        end
        if r < r_min
            r_min, dp_c, dp_c_exact = r, dp, dp_exact
        end
        w = if spherical_flow
            2.0 * Erebus.flow_metric_face_weights(
                i, j, coords; xcenter=X_CENTRE, ycenter=Y_CENTRE
            )[5]
        else
            1.0
        end
        vent_flux += g.S_vent_grid[i, j] * w * coords.dx * coords.dy
        source_flux += g.DQPF[i, j] * w * coords.dx * coords.dy
    end
    n_inside > 0 || throw(ArgumentError("no P node lies inside the planet"))
    n_inner > 0 ||
        throw(ArgumentError("no P node lies within $(INNER_FRACTION) R of the centre"))
    return (;
        pf=copy(g.pf),
        qxD=copy(g.qxD),
        qyD=copy(g.qyD),
        l2_inner=sqrt(sq_err / sq_ref),
        flux_l2_inner=_flux_l2_inner(g, coords, d; jump=jump),
        dp_centre=dp_c,
        dp_centre_exact=dp_c_exact,
        vent_flux,
        source_flux,
        max_solid_speed=max(maximum(abs, g.vx), maximum(abs, g.vy)),
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
    :S_SOURCE,
    :MU_F,
    :K_PERM,
    :K_CORE,
    :A_JUMP,
    :JUMP_EXCLUDE_CELLS,
    :PHI,
    :ETA_SOLID,
    :RX_AIR,
    :PF_INIT,
    :P_AMB,
    :T_SURF,
    :RHO_SOLID,
    :RHOCP,
    :RHO_FLUID,
    :STRENGTH_FILL,
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
    "profile_l2_inner_off_nx33",
    "profile_l2_inner_off_nx65",
    "profile_l2_inner_off_nx129",
    "profile_l2_inner_on_nx33",
    "profile_l2_inner_on_nx65",
    "profile_l2_inner_on_nx129",
    "profile_order_off_65_129",
    "profile_order_on_65_129",
    "centre_ratio_on_off_nx65",
    "flux_l2_inner_off_nx33",
    "flux_l2_inner_off_nx65",
    "flux_l2_inner_off_nx129",
    "flux_l2_inner_on_nx33",
    "flux_l2_inner_on_nx65",
    "flux_l2_inner_on_nx129",
    "vent_closure_off_nx65",
    "vent_closure_on_nx65",
    "solid_speed_max_nx65",
    "profile_l2_inner_jump_off_nx33",
    "profile_l2_inner_jump_off_nx65",
    "profile_l2_inner_jump_off_nx129",
    "profile_l2_inner_jump_on_nx33",
    "profile_l2_inner_jump_on_nx65",
    "profile_l2_inner_jump_on_nx129",
    "centre_ratio_jump_on_off_nx65",
    "flux_l2_inner_jump_off_nx65",
    "flux_l2_inner_jump_off_nx129",
    "flux_l2_inner_jump_on_nx65",
    "flux_l2_inner_jump_on_nx129",
)

const PHYSICS_REGRESSION_BANDS = Dict{
    String,NamedTuple{(:atol, :rtol),Tuple{Float64,Float64}}
}(
    "profile_l2_inner_off_nx33" => (; atol=0.005, rtol=0.0),
    "profile_l2_inner_off_nx65" => (; atol=0.005, rtol=0.0),
    "profile_l2_inner_off_nx129" => (; atol=0.005, rtol=0.0),
    "profile_l2_inner_on_nx33" => (; atol=0.005, rtol=0.0),
    "profile_l2_inner_on_nx65" => (; atol=0.005, rtol=0.0),
    "profile_l2_inner_on_nx129" => (; atol=0.005, rtol=0.0),
    "profile_order_off_65_129" => (; atol=0.05, rtol=0.0),
    "profile_order_on_65_129" => (; atol=0.05, rtol=0.0),
    "flux_l2_inner_off_nx33" => (; atol=0.008, rtol=0.0),
    "flux_l2_inner_off_nx65" => (; atol=0.003, rtol=0.0),
    "flux_l2_inner_off_nx129" => (; atol=0.0015, rtol=0.0),
    "flux_l2_inner_on_nx33" => (; atol=0.008, rtol=0.0),
    "flux_l2_inner_on_nx65" => (; atol=0.003, rtol=0.0),
    "flux_l2_inner_on_nx129" => (; atol=0.0015, rtol=0.0),
    "centre_ratio_on_off_nx65" => (; atol=0.002, rtol=0.0),
    "vent_closure_off_nx65" => (; atol=1.0e-8, rtol=0.0),
    "vent_closure_on_nx65" => (; atol=1.0e-8, rtol=0.0),
    "solid_speed_max_nx65" => (; atol=1.0e-12, rtol=0.0),
    "profile_l2_inner_jump_off_nx33" => (; atol=0.005, rtol=0.0),
    "profile_l2_inner_jump_off_nx65" => (; atol=0.005, rtol=0.0),
    "profile_l2_inner_jump_off_nx129" => (; atol=0.005, rtol=0.0),
    "profile_l2_inner_jump_on_nx33" => (; atol=0.005, rtol=0.0),
    "profile_l2_inner_jump_on_nx65" => (; atol=0.005, rtol=0.0),
    "profile_l2_inner_jump_on_nx129" => (; atol=0.005, rtol=0.0),
    "flux_l2_inner_jump_off_nx65" => (; atol=0.003, rtol=0.0),
    "flux_l2_inner_jump_off_nx129" => (; atol=0.0015, rtol=0.0),
    "flux_l2_inner_jump_on_nx65" => (; atol=0.003, rtol=0.0),
    "flux_l2_inner_jump_on_nx129" => (; atol=0.0015, rtol=0.0),
    "centre_ratio_jump_on_off_nx65" => (; atol=0.002, rtol=0.0),
)

const REPLAY_RTOL = 1.0e-6
const REPLAY_ATOL_BOUNDS = Dict{String,Float64}(
    "vent_closure_off_nx65" => 1.0e-10,
    "vent_closure_on_nx65" => 1.0e-10,
    "solid_speed_max_nx65" => 1.0e-14,
)

"""
    metrics_from_runs(runs)

Return the validation metrics that the `run_case` results in `runs` determine.

The profile and flux errors exist for each grid of `NX_LIST` present in `runs`. The convergence order
of a setting needs both the `NX_REFERENCE` and the finest grid of `NX_LIST`, and its closure
needs the `NX_REFERENCE` grid. The centre ratio and the solid speed need the `NX_REFERENCE`
grid of both settings. The permeability-jump runs have the keys `(N, spherical_flow, :jump)`
and give the ids with `jump_` after the error name, in the same way, except that the flux
error exists for the grids `NX_REFERENCE` and finer only and there is no order, closure or
solid speed. The jump error depends on how the grid cuts the circle `r = A_JUMP` (with the option on
0.025 at `Nx = 63` and 0.034 at `Nx = 65`), so an order between 2 grids has no meaning there.

Parameters
----------
runs : AbstractDict
    Map from `(N, spherical_flow)` and `(N, spherical_flow, :jump)` to the `run_case`
    result.

Returns
-------
Dict{String,Float64}
    Metric values by id (see `EXPECTED_METRIC_IDS`); an id whose runs are missing is
    absent.
"""
function metrics_from_runs(runs::AbstractDict)::Dict{String,Float64}
    n_fine = maximum(NX_LIST)
    metrics = Dict{String,Float64}()
    for (sf, tag) in ((false, "off"), (true, "on"))
        for n in NX_LIST
            haskey(runs, (n, sf)) || continue
            metrics["profile_l2_inner_$(tag)_nx$(n)"] = runs[(n, sf)].l2_inner
            metrics["flux_l2_inner_$(tag)_nx$(n)"] = runs[(n, sf)].flux_l2_inner
        end
        if haskey(runs, (NX_REFERENCE, sf)) && haskey(runs, (n_fine, sf))
            metrics["profile_order_$(tag)_$(NX_REFERENCE)_$(n_fine)"] = log2(
                runs[(NX_REFERENCE, sf)].l2_inner / runs[(n_fine, sf)].l2_inner
            )
        end
        if haskey(runs, (NX_REFERENCE, sf))
            m = runs[(NX_REFERENCE, sf)]
            metrics["vent_closure_$(tag)_nx$(NX_REFERENCE)"] = abs(
                m.vent_flux / m.source_flux - 1.0
            )
        end
        for n in NX_LIST
            haskey(runs, (n, sf, :jump)) || continue
            metrics["profile_l2_inner_jump_$(tag)_nx$(n)"] = runs[(n, sf, :jump)].l2_inner
            if n >= NX_REFERENCE
                metrics["flux_l2_inner_jump_$(tag)_nx$(n)"] =
                    runs[(n, sf, :jump)].flux_l2_inner
            end
        end
    end
    if haskey(runs, (NX_REFERENCE, true)) && haskey(runs, (NX_REFERENCE, false))
        on, off = runs[(NX_REFERENCE, true)], runs[(NX_REFERENCE, false)]
        metrics["centre_ratio_on_off_nx$(NX_REFERENCE)"] = on.dp_centre / off.dp_centre
        metrics["solid_speed_max_nx$(NX_REFERENCE)"] = max(
            on.max_solid_speed, off.max_solid_speed
        )
    end
    if haskey(runs, (NX_REFERENCE, true, :jump)) &&
        haskey(runs, (NX_REFERENCE, false, :jump))
        on, off = runs[(NX_REFERENCE, true, :jump)], runs[(NX_REFERENCE, false, :jump)]
        metrics["centre_ratio_jump_on_off_nx$(NX_REFERENCE)"] = on.dp_centre / off.dp_centre
    end
    return metrics
end

"""
    evaluate_spherical_darcy_flow_benchmark(; nx_list=NX_LIST)

Run the uniform-permeability case and the permeability-jump case at every grid in `nx_list`
with `spherical_flow` off and on and return the runs and the validation metrics they
determine.

Parameters
----------
nx_list : Tuple{Vararg{Int}}, optional
    Grids to run, a subset of `NX_LIST` (default `NX_LIST`). With the default the result
    holds all 29 metrics.

Returns
-------
NamedTuple
    `runs`: `Dict` from `(N, spherical_flow)` and `(N, spherical_flow, :jump)` to the
    `run_case` result;
    `metrics`: `Dict{String,Float64}` from metric id (see `EXPECTED_METRIC_IDS`) to value,
    from `metrics_from_runs`.

Raises
------
ArgumentError
    If a grid in `nx_list` is not in `NX_LIST`.
"""
function evaluate_spherical_darcy_flow_benchmark(; nx_list=NX_LIST)
    all(n -> n isa Integer && n in NX_LIST, nx_list) ||
        throw(ArgumentError("nx_list must be a subset of $NX_LIST, got $nx_list"))
    runs = Dict{Tuple,Any}()
    for n in nx_list, sf in (false, true)
        runs[(n, sf)] = run_case(n; spherical_flow=sf)
        runs[(n, sf, :jump)] = run_case(n; spherical_flow=sf, jump=true)
        n == maximum(NX_LIST) && GC.gc()
    end
    return (; runs, metrics=metrics_from_runs(runs))
end

"""
    expected_metric_references()

Return the reference value of each validation metric: exact profile and flux error `0.0`,
first-order convergence `1.0`, centre ratio `2/3` (also with the permeability jump), and
`0.0` for the closure and solid-speed bounds.

Returns
-------
Dict{String,Float64}
    Mapping from metric id to its reference value.
"""
function expected_metric_references()::Dict{String,Float64}
    refs = Dict{String,Float64}()
    for mid in EXPECTED_METRIC_IDS
        refs[mid] = startswith(mid, "profile_order") ? 1.0 : 0.0
    end

    refs["centre_ratio_on_off_nx65"] = 2.0 / 3.0
    refs["centre_ratio_jump_on_off_nx65"] = 2.0 / 3.0
    return refs
end

"""
    _canonical_field_repr(val)

Return a fixed text form of a config value for the config hash: Bool, integer, float,
Symbol, String and `nothing` in fixed forms, tuples and vectors element by element,
dictionaries with their pairs sorted by key, and any other value by `string(val)`.

Parameters
----------
val : Any
    Value of a config field or of a hashed module constant.

Returns
-------
String
    Text form of `val`.

Raises
------
InexactError
    If `val` is an integer outside the `Int64` range.
"""
function _canonical_field_repr(val)::String
    if val isa Bool
        return val ? "true" : "false"
    elseif val isa Integer
        return string(Int64(val))
    elseif val isa AbstractFloat
        return repr(Float64(val))
    elseif val isa Symbol
        return ":" * String(val)
    elseif val isa AbstractString
        return repr(String(val))
    elseif val === nothing
        return "nothing"
    elseif val isa Tuple || val isa AbstractVector
        return "[" * join((_canonical_field_repr(x) for x in val), ", ") * "]"
    elseif val isa AbstractDict
        isempty(val) && return "{}"
        sorted_pairs = sort!(collect(pairs(val)); by=p -> string(first(p)))
        inner = join(
            [
                "$(_canonical_field_repr(k)) => $(_canonical_field_repr(v))" for
                (k, v) in sorted_pairs
            ],
            ", ",
        )
        return "{$inner}"
    end
    return string(val)
end

"""
    canonical_config_text(cfg)

Return the text that `canonical_config_sha256` hashes: for each section in
`BENCHMARK_CONFIG_SECTIONS` in that order, one line `section.field=value` per field of the
section with the fields sorted by name, followed by one line `case.NAME=value` for each
module constant in `HASHED_CASE_CONSTANTS` in that order. Booleans, integers, floats
(`repr` of the `Float64`), symbols, strings, `nothing`, tuples and vectors have a fixed text
form. The sections not in `BENCHMARK_CONFIG_SECTIONS` and the module constants not in
`HASHED_CASE_CONSTANTS` do not enter the text.

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
        "k_core_m2" => K_CORE,
        "a_jump_m" => A_JUMP,
        "S_source_per_s" => S_SOURCE,
        "phi" => cfg.thermodynamics.phim0,
        "eta_solid_Pa_s" => cfg.materials.etasolidm[1],
        "k_vent_m2" => cfg.venting.k_vent,
        "T_surf_K" => cfg.disk.t_ambient,
        "P_amb_Pa" => cfg.geometry.psurface,
        "inner_fraction" => INNER_FRACTION,
        "spherical_metric" => cfg.geometry.spherical_metric,
        "darcy_elimination" => cfg.solver.darcy_elimination,
    )
end

"""
    check_physics_regression_bands(metrics, record_path)

Check that each metric lies within its regression band (`PHYSICS_REGRESSION_BANDS`) around
the recorded value in `record_path`.

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
    raw = TOML.parsefile(record_path)
    recorded_by_id = Dict{String,Float64}()
    for m in get(raw, "metrics", Any[])
        if m isa AbstractDict && haskey(m, "id") && get(m, "value", nothing) isa Real
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
    _same_julia_major_minor(recorded_ver, current_ver=VERSION)

Return whether the Julia version of a record and the running version share major and minor.

Parameters
----------
recorded_ver : AbstractString
    Julia version stored in the record, for example `"1.12.1"`.
current_ver : VersionNumber, optional
    Running Julia version.

Returns
-------
Bool
    `true` when major and minor agree; `false` when they differ or `recorded_ver` does not
    start with a version.
"""
function _same_julia_major_minor(
    recorded_ver::AbstractString, current_ver::VersionNumber=VERSION
)::Bool
    m = match(r"^v?(\d+)\.(\d+)(?:\.|\b|$)", strip(recorded_ver))
    m === nothing && return false
    return parse(Int, m.captures[1]) == current_ver.major &&
           parse(Int, m.captures[2]) == current_ver.minor
end

"""
    check_validation_reproducibility(bench, record_path; check_recorded_values=true)

Compare freshly evaluated metrics, reference values, regression bands and the configuration
hash and summary with a committed validation record. The replay of recorded values
(`REPLAY_RTOL`, `REPLAY_ATOL_BOUNDS`) runs when `check_recorded_values` is `true` and the
Julia major.minor version matches `metadata.julia_version`; the regression bands always run.

Parameters
----------
bench : NamedTuple
    Result of `evaluate_spherical_darcy_flow_benchmark`.
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
    raw = TOML.parsefile(record_path)
    ok_schema, errs = validate_record_dict(
        raw; require_doc_pages=false, check_doc_sync=false
    )
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

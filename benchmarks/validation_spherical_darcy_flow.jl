"""
    SphericalDarcyFlowValidation

Steady radial Darcy flow from a uniform fluid source in a rigid porous sphere.

The exact excess fluid pressure for a source `S` [1/s], permeability `k` and fluid
viscosity `mu`, with the surface held at the vent pressure (Dirichlet limit of the
venting condition), is `p - P_vent = mu S (R^2 - r^2) / (2 d k)`, with `d = 3` for the
sphere (`geometry.spherical_flow = true`) and `d = 2` for the disk (`false`).
The case builds the state with `init_simulation`, sets the grid fields directly and calls
`solve_stokes_darcy!` once on timestep 1, where the poroelastic compressibilities are off.
"""
module SphericalDarcyFlowValidation

using Erebus
using Erebus.Config
using Logging
using StaticArrays

const R_PLANET = 50_000.0
const L_BOX = 140_000.0
const X_CENTRE = 70_000.0
const Y_CENTRE = 70_000.0
const MU_F = 1.0e-3
const K_PERM = 1.0e-13
const S_SOURCE = 1.0e-11
const T_SURF = 200.0
const P_AMB = 100.0
const PHI = 0.2
const ETA_SOLID = 1.0e28
const RX_AIR = 1.0e30
const K_VENT = 1.0e-6
const PF_INIT = 1.0e5
const DT_YR = 1.0e-10
const INNER_FRACTION = 0.9

"""
    make_config(N; spherical_flow)

Return the `SimulationConfig` of the case on an `N` x `N` grid.

Parameters
----------
N : Int
    Basic grid resolution in x and y.
spherical_flow : Bool
    Value of `geometry.spherical_flow`.

Returns
-------
SimulationConfig
"""
function make_config(N::Int; spherical_flow::Bool)
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
        solver=SolverConfig(; darcy_elimination=true, etamax=1.0e30, etamin=1.0e12),
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

_inside(coords, i, j) = hypot(coords.xp[j] - X_CENTRE, coords.yp[i] - Y_CENTRE) <= R_PLANET

"""
    set_fields!(state, coords)

Set the grid fields of a rigid, uniformly permeable sphere with a uniform fluid source:
impermeable faces outside the sphere, no gravity, uniform porosity and fluid pressure.
"""
function set_fields!(state, coords)
    g = state.grids
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    for f in (g.ETA, g.ETA0, g.ETAP)
        fill!(f, ETA_SOLID)
    end
    fill!(g.FRI, 0.0)
    for f in (g.COH, g.TEN, g.GGG, g.GGGP)
        fill!(f, 1.0e20)
    end
    fill!(g.tk1, T_SURF)
    fill!(g.tk2, T_SURF)
    for f in (g.PHI, g.PHIX, g.PHIY)
        fill!(f, PHI)
    end
    for f in (g.RHO, g.RHOX, g.RHOY)
        fill!(f, 3000.0)
    end
    fill!(g.RHOCP, 3.0e6)
    fill!(g.RHOFX, 1000.0)
    fill!(g.RHOFY, 1000.0)
    fill!(g.pr, P_AMB)
    fill!(g.pr0, P_AMB)
    fill!(g.pf, PF_INIT)
    fill!(g.pf0, PF_INIT)
    fill!(g.gx, 0.0)
    fill!(g.gy, 0.0)
    rx_rock = MU_F / K_PERM
    for j in 1:Nx1, i in 1:Ny1
        inside = _inside(coords, i, j)
        g.RX[i, j] = (j < Nx1 && inside && _inside(coords, i, j + 1)) ? rx_rock : RX_AIR
        g.RY[i, j] = (i < Ny1 && inside && _inside(coords, i + 1, j)) ? rx_rock : RX_AIR
        g.DQPF[i, j] = inside ? S_SOURCE : 0.0
    end
    return nothing
end

"""
    run_case(N; spherical_flow)

Solve the case once and compare the fluid pressure with the exact solution.

Returns
-------
NamedTuple
    `l2_inner`: relative L2 error of `p - P_vent` over P nodes with `r <= 0.9 R`;
    `dp_centre`: computed `p - P_vent` at the P node nearest the centre;
    `dp_centre_exact`: exact value at that node;
    `vent_flux`, `source_flux`: vented and produced fluid volume rate [m^3/s per slice
    weight], summed with the cell weights of the continuity rows;
    `max_solid_speed`: largest solid velocity component [m/s].
"""
function run_case(N::Int; spherical_flow::Bool)
    cfg = make_config(N; spherical_flow=spherical_flow)
    out = mktempdir()
    state, coords, ws, tio, _ = with_logger(NullLogger()) do
        Erebus.init_simulation(cfg; output_path=out)
    end
    tio === nothing || close(tio)
    state.timestep = 1
    state.accumulators.P_amb = P_AMB
    set_fields!(state, coords)
    with_logger(NullLogger()) do
        Erebus.solve_stokes_darcy!(state, coords, cfg, ws; titer=1, dt_step_initial=state.dt)
    end
    P_vent = Erebus.compute_venting_pressure(
        T_SURF, P_AMB; species=:H2O, L_sub=cfg.venting.L_sublimation
    )
    return _metrics(state, coords, P_vent, spherical_flow)
end

function _metrics(state, coords, P_vent, spherical_flow)
    g = state.grids
    d = spherical_flow ? 3.0 : 2.0
    sq_err = sq_ref = 0.0
    vent_flux = source_flux = 0.0
    r_min, dp_c, dp_c_exact = Inf, NaN, NaN
    for j in 3:(coords.Nx1 - 2), i in 3:(coords.Ny1 - 2)
        r = hypot(coords.xp[j] - X_CENTRE, coords.yp[i] - Y_CENTRE)
        r <= R_PLANET || continue
        dp = g.pf[i, j] - P_vent
        dp_exact = MU_F * S_SOURCE * (R_PLANET^2 - r^2) / (2d * K_PERM)
        if r <= INNER_FRACTION * R_PLANET
            sq_err += (dp - dp_exact)^2
            sq_ref += dp_exact^2
        end
        if r < r_min
            r_min, dp_c, dp_c_exact = r, dp, dp_exact
        end
        w = if spherical_flow
            2.0 *
            Erebus.flow_metric_face_weights(i, j, coords; xcenter=X_CENTRE, ycenter=Y_CENTRE)[5]
        else
            1.0
        end
        vent_flux += g.S_vent_grid[i, j] * w * coords.dx * coords.dy
        source_flux += g.DQPF[i, j] * w * coords.dx * coords.dy
    end
    return (;
        l2_inner=sqrt(sq_err / sq_ref),
        dp_centre=dp_c,
        dp_centre_exact=dp_c_exact,
        vent_flux,
        source_flux,
        max_solid_speed=max(maximum(abs, g.vx), maximum(abs, g.vy)),
    )
end

end # module

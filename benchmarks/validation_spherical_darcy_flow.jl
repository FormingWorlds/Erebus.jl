"""
    SphericalDarcyFlowValidation

Steady radial Darcy flow from a uniform fluid source in a rigid porous sphere.

The exact excess fluid pressure for a source `S` [1/s], permeability `k` and fluid
viscosity `mu`, with the surface held at the vent pressure (Dirichlet limit of the
venting condition), is `p - P_vent = mu S (R^2 - r^2) / (2 d k)`, with `d = 3` for the
sphere (`geometry.spherical_flow = true`) and `d = 2` for the disk (`false`).
The ice saturation pressure at `T_SURF` is below `P_AMB`, so `P_vent = P_AMB`.
The exact radial Darcy flux is `q_r = S r / d`.
The case builds the state with `init_simulation`, sets the grid fields directly and calls
`solve_stokes_darcy!` once on timestep 1, where the poroelastic compressibilities are off.
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
    run_case(N; spherical_flow, darcy_elimination=true)

Solve the case once and compare the fluid pressure and the radial Darcy flux with the exact
solution.

Parameters
----------
N : Int
    Basic grid resolution in x and y, at least `NX_MIN`.
spherical_flow : Bool
    Value of `geometry.spherical_flow`.
darcy_elimination : Bool, optional
    Value of `solver.darcy_elimination` (default `true`).

Returns
-------
NamedTuple
    `converged`: `plastic_converged` flag returned by `solve_stokes_darcy!`;
    `pf`: fluid pressure grid [Pa];
    `l2_inner`: relative L2 error of `p - P_AMB` over P nodes with `r <= 0.9 R`;
    `flux_l2_inner`: relative L2 error of the radial Darcy flux over faces with
    `r <= 0.9 R`;
    `dp_centre`: computed `p - P_AMB` at the P node nearest the centre;
    `dp_centre_exact`: exact value at that node;
    `vent_flux`, `source_flux`: vented and produced fluid volume rate [m^3/s per slice
    weight], summed with the cell weights of the continuity rows;
    `max_solid_speed`: largest solid velocity component [m/s].

Raises
------
ArgumentError
    If `N < NX_MIN`.
"""
function run_case(N::Int; spherical_flow::Bool, darcy_elimination::Bool=true)
    cfg = make_config(N; spherical_flow=spherical_flow, darcy_elimination=darcy_elimination)
    out = mktempdir()
    state, coords, ws, tio, _ = with_logger(NullLogger()) do
        return Erebus.init_simulation(cfg; output_path=out)
    end
    tio === nothing || close(tio)
    state.timestep = 1
    state.accumulators.P_amb = P_AMB
    set_fields!(state, coords)
    res = with_logger(NullLogger()) do
        return Erebus.solve_stokes_darcy!(
            state, coords, cfg, ws; titer=1, dt_step_initial=state.dt
        )
    end
    return (; converged=res.plastic_converged, _metrics(state, coords, spherical_flow)...)
end

"""
    _radial_flux_l2_inner(g, coords, d)

Return the relative L2 error of the radial Darcy flux against the exact `S r / d` over the
faces with `r <= INNER_FRACTION * R_PLANET`.

At each `qxD` face the radial flux is `qxD cos(theta) + qy sin(theta)` with `qy` the mean of
the 4 surrounding `qyD` faces, and likewise at each `qyD` face with `qx` the mean of the 4
surrounding `qxD` faces. The `qxD[i, j]` face lies at `(coords.xvx[j], coords.yvx[i])` and the
`qyD[i, j]` face at `(coords.xvy[j], coords.yvy[i])`.

Parameters
----------
g : Grids
    Solved grid fields with `qxD` and `qyD` [m/s].
coords : GridCoordinates
    Staggered grid coordinates.
d : Float64
    Flow dimension, 3 for the sphere and 2 for the disk.

Returns
-------
Float64
    Relative L2 error, dimensionless.
"""
function _radial_flux_l2_inner(g, coords, d::Float64)
    sq_err = sq_ref = 0.0
    r_max = INNER_FRACTION * R_PLANET
    # x faces: i in 2:Ny, j in 2:(Nx - 1); the 4 qyD neighbours are (i - 1:i, j:(j + 1)).
    for j in 2:(coords.Nx - 1), i in 2:coords.Ny
        dxf = coords.xvx[j] - X_CENTRE
        dyf = coords.yvx[i] - Y_CENTRE
        r = hypot(dxf, dyf)
        (0.0 < r <= r_max) || continue
        qy = 0.25 * (g.qyD[i - 1, j] + g.qyD[i - 1, j + 1] + g.qyD[i, j] + g.qyD[i, j + 1])
        q_r = (g.qxD[i, j] * dxf + qy * dyf) / r
        sq_err += (q_r - S_SOURCE * r / d)^2
        sq_ref += (S_SOURCE * r / d)^2
    end
    # y faces: i in 2:(Ny - 1), j in 2:Nx; the 4 qxD neighbours are (i:(i + 1), j - 1:j).
    for j in 2:coords.Nx, i in 2:(coords.Ny - 1)
        dxf = coords.xvy[j] - X_CENTRE
        dyf = coords.yvy[i] - Y_CENTRE
        r = hypot(dxf, dyf)
        (0.0 < r <= r_max) || continue
        qx = 0.25 * (g.qxD[i, j - 1] + g.qxD[i, j] + g.qxD[i + 1, j - 1] + g.qxD[i + 1, j])
        q_r = (qx * dxf + g.qyD[i, j] * dyf) / r
        sq_err += (q_r - S_SOURCE * r / d)^2
        sq_ref += (S_SOURCE * r / d)^2
    end
    return sqrt(sq_err / sq_ref)
end

function _metrics(state, coords, spherical_flow)
    g = state.grids
    d = spherical_flow ? 3.0 : 2.0
    sq_err = sq_ref = 0.0
    vent_flux = source_flux = 0.0
    r_min, dp_c, dp_c_exact = Inf, NaN, NaN
    for j in 3:(coords.Nx1 - 2), i in 3:(coords.Ny1 - 2)
        r = hypot(coords.xp[j] - X_CENTRE, coords.yp[i] - Y_CENTRE)
        r <= R_PLANET || continue
        dp = g.pf[i, j] - P_AMB
        dp_exact = MU_F * S_SOURCE * (R_PLANET^2 - r^2) / (2d * K_PERM)
        if r <= INNER_FRACTION * R_PLANET
            sq_err += (dp - dp_exact)^2
            sq_ref += dp_exact^2
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
    return (;
        pf=copy(g.pf),
        l2_inner=sqrt(sq_err / sq_ref),
        flux_l2_inner=_radial_flux_l2_inner(g, coords, d),
        dp_centre=dp_c,
        dp_centre_exact=dp_c_exact,
        vent_flux,
        source_flux,
        max_solid_speed=max(maximum(abs, g.vx), maximum(abs, g.vy)),
    )
end

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
    "flux_l2_inner_off_nx33" => (; atol=0.005, rtol=0.0),
    "flux_l2_inner_off_nx65" => (; atol=0.002, rtol=0.0),
    "flux_l2_inner_off_nx129" => (; atol=0.001, rtol=0.0),
    "flux_l2_inner_on_nx33" => (; atol=0.005, rtol=0.0),
    "flux_l2_inner_on_nx65" => (; atol=0.002, rtol=0.0),
    "flux_l2_inner_on_nx129" => (; atol=0.001, rtol=0.0),
    "centre_ratio_on_off_nx65" => (; atol=0.002, rtol=0.0),
    "vent_closure_off_nx65" => (; atol=1.0e-8, rtol=0.0),
    "vent_closure_on_nx65" => (; atol=1.0e-8, rtol=0.0),
    "solid_speed_max_nx65" => (; atol=1.0e-12, rtol=0.0),
)

const REPLAY_RTOL = 1.0e-6
const REPLAY_ATOL_BOUNDS = Dict{String,Float64}(
    "vent_closure_off_nx65" => 1.0e-10,
    "vent_closure_on_nx65" => 1.0e-10,
    "solid_speed_max_nx65" => 1.0e-14,
)

"""
    evaluate_spherical_darcy_flow_benchmark()

Run the case at every grid in `NX_LIST` with `spherical_flow` off and on and return the 18
validation metrics.

Returns
-------
NamedTuple
    `runs`: `Dict` from `(N, spherical_flow)` to the `run_case` result;
    `metrics`: `Dict{String,Float64}` from metric id (see `EXPECTED_METRIC_IDS`) to value.
"""
function evaluate_spherical_darcy_flow_benchmark()
    runs = Dict(
        (n, sf) => run_case(n; spherical_flow=sf) for n in NX_LIST, sf in (false, true)
    )
    metrics = Dict{String,Float64}()
    for (sf, tag) in ((false, "off"), (true, "on"))
        for n in NX_LIST
            metrics["profile_l2_inner_$(tag)_nx$(n)"] = runs[(n, sf)].l2_inner
            metrics["flux_l2_inner_$(tag)_nx$(n)"] = runs[(n, sf)].flux_l2_inner
        end
        metrics["profile_order_$(tag)_65_129"] = log2(
            runs[(65, sf)].l2_inner / runs[(129, sf)].l2_inner
        )
        m = runs[(NX_REFERENCE, sf)]
        metrics["vent_closure_$(tag)_nx65"] = abs(m.vent_flux / m.source_flux - 1.0)
    end
    metrics["centre_ratio_on_off_nx65"] =
        runs[(NX_REFERENCE, true)].dp_centre / runs[(NX_REFERENCE, false)].dp_centre
    metrics["solid_speed_max_nx65"] = max(
        runs[(NX_REFERENCE, false)].max_solid_speed,
        runs[(NX_REFERENCE, true)].max_solid_speed,
    )
    return (; runs, metrics)
end

"""
    expected_metric_references()

Return the reference value of each validation metric: exact profile and flux error `0.0`,
first-order convergence `1.0`, centre ratio `2/3`, and `0.0` for the closure and
solid-speed bounds.

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
    return refs
end

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
    end
    return string(val)
end

"""
    canonical_config_sha256(cfg)

Return the SHA-256 digest of a text that lists, for each section in
`BENCHMARK_CONFIG_SECTIONS` in that order, one line `section.field=value` per field of the
section with the fields sorted by name, followed by the 5 lines `case.NAME=value` for the
module constants `S_SOURCE`, `L_BOX`, `INNER_FRACTION`, `RX_AIR` and `PF_INIT`. Booleans,
integers, floats (`repr` of the `Float64`), symbols, strings, `nothing` and vectors have a
fixed text form. The sections not in `BENCHMARK_CONFIG_SECTIONS` and the other module
constants do not enter the digest.

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
    io = IOBuffer()
    for sec_sym in BENCHMARK_CONFIG_SECTIONS
        sec = getproperty(cfg, sec_sym)
        for f_sym in sort(collect(fieldnames(typeof(sec))))
            val = getfield(sec, f_sym)
            println(io, sec_sym, ".", f_sym, "=", _canonical_field_repr(val))
        end
    end
    println(io, "case.S_SOURCE=", _canonical_field_repr(S_SOURCE))
    println(io, "case.L_BOX=", _canonical_field_repr(L_BOX))
    println(io, "case.INNER_FRACTION=", _canonical_field_repr(INNER_FRACTION))
    println(io, "case.RX_AIR=", _canonical_field_repr(RX_AIR))
    println(io, "case.PF_INIT=", _canonical_field_repr(PF_INIT))
    return bytes2hex(sha256(take!(io)))
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

module MeltingShellValidation

# Melting shell case: the porous shell of a planetesimal (dry porosity laws of Henke et
# al. 2012, 26Al heating) with silicate melting and latent heat, conduction only. The
# melt has the density of the solid (no melt buoyancy) and soft turbulence is off.

using Erebus
using StaticArrays
using Logging
using SHA

isdefined(@__MODULE__, :PorousShellValidation) ||
    include(joinpath(@__DIR__, "validation_porous_shell.jl"))

const PSV = PorousShellValidation
# The record checks of tools/check_validation.jl read `validate_record_dict` of this module.
using .PorousShellValidation: validate_record_dict

const L_MELT = 4.0e5
const T_SOL = 1400.0
const T_LIQ = 1800.0

"""
    melt_fraction_ref(T; T_sol=T_SOL, T_liq=T_LIQ) -> Float64

Melt fraction of rock at temperature `T` [K]: linear between `T_sol` and `T_liq`, clamped
to [0, 1] (no pressure shift).
"""
function melt_fraction_ref(T; T_sol=T_SOL, T_liq=T_LIQ)
    return clamp((T - T_sol) / (T_liq - T_sol), 0.0, 1.0)
end

"""
    rhocp_app_ref(T, phi; L_melt=L_MELT, T_sol=T_SOL, T_liq=T_LIQ) -> Float64

Apparent volumetric heat capacity [J/(m^3 K)] of rock with dry porosity `phi` at `T` [K]:
`(1 - phi) (RHOCP_SOL + RHO_SOL L_melt / (T_liq - T_sol))` for `T_sol < T < T_liq`, the
sensible term alone outside.
"""
function rhocp_app_ref(T, phi; L_melt=L_MELT, T_sol=T_SOL, T_liq=T_LIQ)
    c = PSV.RHOCP_SOL + (T_sol < T < T_liq ? PSV.RHO_SOL * L_melt / (T_liq - T_sol) : 0.0)
    return (1.0 - phi) * c
end

"""
    enthalpy_s_ref(T; L_melt=L_MELT, T_sol=T_SOL, T_liq=T_LIQ, T0=PSV.T_SPACE) -> Float64

Enthalpy per unit volume of dense rock [J/m^3] from `T0` to `T`: sensible heat plus
`RHO_SOL L_melt F(T)`.
"""
function enthalpy_s_ref(T; L_melt=L_MELT, T_sol=T_SOL, T_liq=T_LIQ, T0=PSV.T_SPACE)
    F = melt_fraction_ref(T; T_sol=T_sol, T_liq=T_liq)
    return PSV.RHOCP_SOL * (T - T0) + PSV.RHO_SOL * L_melt * F
end

"""
    make_melting_shell_config(; Nx, phi_init, t_form_yr, n_steps, dt_yr,
                              L_melt=L_MELT, T_sol=T_SOL, T_liq=T_LIQ, kwargs...)

Configuration of the porous shell case (`make_porous_shell_config`, which takes the other
keywords) with `melting.active`, latent heat `L_melt` [J/kg] over the rock melting
interval `T_sol` to `T_liq` [K], the melt density equal to the rock density and soft
turbulence off.

# Throws

    - `ArgumentError`: if `T_sol >= T_liq`, or the errors of `make_porous_shell_config`
"""
function make_melting_shell_config(;
    L_melt::Float64=L_MELT, T_sol::Float64=T_SOL, T_liq::Float64=T_LIQ, kwargs...
)
    T_sol < T_liq || throw(ArgumentError("T_sol must be below T_liq, got $T_sol, $T_liq"))
    melting = MeltingConfig(;
        active=true,
        soft_turbulence=false,
        L_melt=L_melt,
        T_solidus=SVector{3,Float64}([T_sol, T_sol, NaN]),
        T_liquidus=SVector{3,Float64}([T_liq, T_liq, NaN]),
        rho_melt=PSV.RHO_SOL,
    )
    return PSV.make_porous_shell_config(; melting=melting, kwargs...)
end

"""
    melt_radius(xm, ym, T_max; xc, yc, R_p, dr, T_sol=T_SOL, nsector=16) -> Float64

Radius [m] of the molten centre from the peak temperatures `T_max` [K] of rock markers at
`(xm, ym)`. The markers with `r < R_p` about `(xc, yc)` fall into `nsector` equal angular
sectors and centre-anchored radial bins `[(b - 1) dr, min(b dr, R_p))`. In each sector the
scan runs outward over the non-empty bins while the bin median of `T_max` is above
`T_sol`; the sector radius is the outer edge of the last such bin (0 if the first non-empty
bin is not molten). The result is the median of the sector radii.

# Throws

    - `ArgumentError`: if `xm`, `ym` and `T_max` differ in length, or no marker is inside `R_p`
    - `DomainError`: if `dr` or `R_p` is not positive and finite, or `nsector < 1`
"""
function melt_radius(xm, ym, T_max; xc, yc, R_p, dr, T_sol=T_SOL, nsector::Int=16)
    length(xm) == length(ym) == length(T_max) ||
        throw(ArgumentError("xm, ym and T_max must have the same length"))
    (isfinite(dr) && dr > 0.0) || throw(DomainError(dr, "dr must be positive and finite"))
    (isfinite(R_p) && R_p > 0.0) ||
        throw(DomainError(R_p, "R_p must be positive and finite"))
    nsector >= 1 || throw(DomainError(nsector, "nsector must be at least 1"))
    nb = ceil(Int, R_p / dr)
    bins = [[Float64[] for _ in 1:nb] for _ in 1:nsector]
    for m in eachindex(xm)
        dx, dy = xm[m] - xc, ym[m] - yc
        r = hypot(dx, dy)
        r < R_p || continue
        s = min(floor(Int, mod(atan(dy, dx), 2pi) / (2pi / nsector)) + 1, nsector)
        push!(bins[s][min(floor(Int, r / dr) + 1, nb)], T_max[m])
    end
    any(sec -> any(!isempty, sec), bins) || throw(ArgumentError("no marker is inside R_p"))
    radii = Float64[]
    for sec in bins
        any(!isempty, sec) || continue
        rs = 0.0
        for b in 1:nb
            isempty(sec[b]) && continue
            PSV._median(sec[b]) > T_sol || break
            rs = min(b * dr, R_p)
        end
        push!(radii, rs)
    end
    return PSV._median(radii)
end

"""
    melt_radius_1d(T_max, R; T_sol=T_SOL) -> Float64

Radius [m] of the molten centre of a 1-D profile of peak temperatures `T_max` [K] on `Nr`
equal cells of a sphere of radius `R` [m]: the outer face of the last cell of the run of
cells from the centre with `T_max > T_sol`; 0 if the centre cell is not molten, `R` if
every cell is.

# Throws

    - `ArgumentError`: if `T_max` is empty
    - `DomainError`: if `R` is not positive and finite
"""
function melt_radius_1d(T_max::AbstractVector{Float64}, R::Float64; T_sol=T_SOL)
    isempty(T_max) && throw(ArgumentError("T_max must not be empty"))
    (isfinite(R) && R > 0.0) || throw(DomainError(R, "R must be positive and finite"))
    i = findfirst(<=(T_sol), T_max)
    return i === nothing ? R : (i - 1) * R / length(T_max)
end

"""
    reference_melt(; L_melt=L_MELT, T_sol=T_SOL, T_liq=T_LIQ) -> NamedTuple

The `melt` keyword of `solve_porous_sphere` for latent heat `L_melt` [J/kg] over the
melting interval `T_sol` to `T_liq` [K].
"""
reference_melt(; L_melt=L_MELT, T_sol=T_SOL, T_liq=T_LIQ) = (; L_melt, T_sol, T_liq)

"""
    run_melting_shell_case(; Nx, phi_init, t_form_yr, n_steps, dt_yr, L_melt=L_MELT,
                           seed=PSV.SEED_DEFAULT, Nr=2000, reference=nothing,
                           run_loop=simulation_loop) -> NamedTuple

Run the Erebus melting shell case with latent heat `L_melt` [J/kg] and compare it with the
1-D reference at the nominal latent heat `L_MELT`, with the same time step and number of
steps. The rock markers fall into centre-anchored radial bins of half a cell, `dx / 2`.
The centre temperature is the median of the marker peak temperatures `tkm_max` in the
first bin, against the 1-D peak temperature of the cell at the bin centre. The melt radius
is `melt_radius` with bins of `dx / 2`, against `melt_radius_1d`. The relative L2 error
of the binned peak temperature rise runs over the bins whose 1-D peak temperature is at
least `TMAX_L2_T_MIN`. The energy residual is `thermal_step_energy_residual` of the last
step over the free nodes.

# Arguments

    - `Nx`, `phi_init`, `t_form_yr` [yr], `n_steps`, `dt_yr` [yr]: as
      `make_melting_shell_config`
    - `L_melt`: latent heat of the Erebus run [J/kg]; the 1-D reference keeps `L_MELT`
    - `seed`: seed of the marker positions
    - `Nr`: cells of the 1-D reference
    - `reference`: a `solve_porous_sphere` result of the same case with `reference_melt()`
      and `stop_cold = false`, or `nothing` to compute it here
    - `run_loop`: `run_loop(cfg; output_path)` runs Erebus and returns the final state

# Returns

`(; dx, T_centre, T_centre_ref, T_centre_error, r_melt, r_melt_ref, r_melt_error_cells,
tmax_rel_l2, energy_lse_residual, enthalpy_residual_ref, T_peak, r, tmax, tmax_ref,
reference)`: temperatures [K], radii [m], the radius error in cells of `dx`, the L2 error
and residuals (dimensionless), the largest marker peak temperature [K], the binned
profiles [m, K] and the full 1-D result.

# Throws

    - `DomainError`, `ArgumentError`: for the arguments that `make_melting_shell_config`
      rejects
    - `ErrorException`: from `PSV.check_run_clock` and `PSV.check_fixed_surface_counts`
"""
function run_melting_shell_case(;
    Nx::Int,
    phi_init::Float64,
    t_form_yr::Float64,
    n_steps::Int,
    dt_yr::Float64,
    L_melt::Float64=L_MELT,
    seed::Int=PSV.SEED_DEFAULT,
    Nr::Int=2000,
    reference=nothing,
    run_loop=simulation_loop,
)
    cfg = make_melting_shell_config(; Nx, phi_init, t_form_yr, n_steps, dt_yr, L_melt, seed)
    st = mktempdir() do tmpdir
        redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                with_logger(NullLogger()) do
                    return run_loop(cfg; output_path=tmpdir)
                end
            end
        end
    end
    PSV.check_run_clock(st, cfg)
    PSV.check_fixed_surface_counts(st, cfg)
    ref = if reference === nothing
        PSV.solve_porous_sphere(;
            R=PSV.R_PLANET,
            phi_init,
            t_form_yr,
            t_end_yr=n_steps * dt_yr,
            dt_yr,
            Nr,
            stop_cold=false,
            melt=reference_melt(),
        )
    else
        reference
    end
    mk = st.markers
    pg = mk.groups.porosity
    xm, ym, tm = mk.core.xm, mk.core.ym, mk.core.tm
    dx = PSV.L_BOX / (Nx - 1)
    xc = yc = PSV.L_BOX / 2.0
    rock = findall(t -> t == 1 || t == 2, tm)
    r = hypot.(xm[rock] .- xc, ym[rock] .- yc)
    rc, tmax = PSV.radial_bin_medians(r, pg.tkm_max[rock], PSV.R_PLANET, dx / 2)
    isapprox(rc[1], dx / 4; rtol=1.0e-12) || error("the centre bin holds no rock marker")
    nr = length(ref.T_max)
    cell = [clamp(floor(Int, x / (PSV.R_PLANET / nr)) + 1, 1, nr) for x in rc]
    tmax_ref = ref.T_max[cell]
    hot = tmax_ref .>= TMAX_L2_T_MIN
    rF = melt_radius(
        xm[rock], ym[rock], pg.tkm_max[rock]; xc, yc, R_p=PSV.R_PLANET, dr=dx / 2
    )
    rF_ref = melt_radius_1d(ref.T_max, PSV.R_PLANET)
    return (;
        dx,
        T_centre=tmax[1],
        T_centre_ref=tmax_ref[1],
        T_centre_error=tmax[1] - tmax_ref[1],
        r_melt=rF,
        r_melt_ref=rF_ref,
        r_melt_error_cells=(rF - rF_ref) / dx,
        tmax_rel_l2=sqrt(
            sum(abs2, tmax[hot] .- tmax_ref[hot]) / sum(abs2, tmax_ref[hot] .- PSV.T_SPACE)
        ),
        energy_lse_residual=PSV.thermal_step_energy_residual(
            st.grids, st.dt, dx, dx; fixed=st.grids.fixed_tk
        ),
        enthalpy_residual_ref=ref.enthalpy_residual,
        T_peak=maximum(pg.tkm_max[rock]),
        r=rc,
        tmax,
        tmax_ref,
        reference=ref,
    )
end

"""
    lag_shift(; dt_yr=DT_YR_CASE, divisor=LAG_DT_DIVISOR, nx=NX_REFERENCE, Nr=NR_REFERENCE,
              coarse=nothing) -> NamedTuple

Shift of the 1-D reference of the recorded case from the time lag of the apparent heat
capacity: the reference at the case step `dt_yr` [yr] minus the reference at
`dt_yr / divisor`, over the same run time. The centre temperature is the peak temperature
of the cell at the centre of the first radial bin of `dx / 2` on an `nx` grid. `coarse`
is the reference at `dt_yr` with `Nr` cells, or `nothing` to compute it here.

# Returns

`(; dT_centre, dr_melt, coarse, fine)`: the shifts of the centre temperature [K] and of
the melt radius [m], and the 2 `solve_porous_sphere` results.

# Throws

    - `DomainError`: if `divisor < 1`
"""
function lag_shift(;
    dt_yr::Float64=DT_YR_CASE,
    divisor::Int=LAG_DT_DIVISOR,
    nx::Int=NX_REFERENCE,
    Nr::Int=NR_REFERENCE,
    coarse=nothing,
)
    divisor >= 1 || throw(DomainError(divisor, "divisor must be at least 1"))
    coarse === nothing && (coarse = reference_solution(; dt_yr, Nr))
    fine = reference_solution(; dt_yr=dt_yr / divisor, Nr)
    a, b = reference_values(coarse, nx), reference_values(fine, nx)
    return (;
        dT_centre=a.centre_tmax - b.centre_tmax,
        dr_melt=a.melt_radius - b.melt_radius,
        coarse,
        fine,
    )
end

# Case of the records: formation 1.15 Myr after the first solids, 140 steps of 25 kyr.
# The 1-D centre peaks near 1661 K, inside the melting interval.
const CASE_ID = "melting_shell_t1.15_phi0.25"
const RECORD_ID = "melting_shell"
const SUITE_RECORD_ID = "melting_shell_suite"
const PHI_INIT_CASE = 0.25
const T_FORM_YR_CASE = 1.15e6
const N_STEPS_CASE = 140
const DT_YR_CASE = 2.5e4
const NX_REFERENCE = 65
const NX_SUITE = 129
const NR_REFERENCE = 2000
const TMAX_L2_T_MIN = 1100.0
const CENTRE_TMAX_ATOL = ((65, 20.0), (129, 10.0))
const TMAX_REL_L2_LIMIT = ((65, 0.04), (129, 0.025))
const MELT_RADIUS_ATOL_CELLS = 1.5
const ENERGY_RESIDUAL_LIMIT = 1.0e-9
const ENTHALPY_RESIDUAL_LIMIT = 0.01
const LAG_DT_DIVISOR = 10
const L_MUTANT_FACTOR = 1.5

"""
    make_case_config(; Nx=NX_REFERENCE, seed=PSV.SEED_DEFAULT, L_melt=L_MELT)

Return the `SimulationConfig` of the recorded case on an `Nx` x `Nx` grid with latent
heat `L_melt` [J/kg].

# Throws

    - `DomainError`: if `Nx < 5`
"""
function make_case_config(;
    Nx::Int=NX_REFERENCE, seed::Int=PSV.SEED_DEFAULT, L_melt::Float64=L_MELT
)
    return make_melting_shell_config(;
        Nx,
        phi_init=PHI_INIT_CASE,
        t_form_yr=T_FORM_YR_CASE,
        n_steps=N_STEPS_CASE,
        dt_yr=DT_YR_CASE,
        L_melt,
        seed,
    )
end

"""
    reference_solution(; dt_yr=DT_YR_CASE, Nr=NR_REFERENCE) -> NamedTuple

Return the 1-D reference of the recorded case: `solve_porous_sphere` with the nominal
latent heat, `Nr` cells, the step `dt_yr` [yr] over the run time of the case
(`N_STEPS_CASE` steps of `DT_YR_CASE`) and no stop once the rock cools.
"""
function reference_solution(; dt_yr::Float64=DT_YR_CASE, Nr::Int=NR_REFERENCE)
    return PSV.solve_porous_sphere(;
        R=PSV.R_PLANET,
        phi_init=PHI_INIT_CASE,
        t_form_yr=T_FORM_YR_CASE,
        t_end_yr=N_STEPS_CASE * DT_YR_CASE,
        dt_yr,
        Nr,
        stop_cold=false,
        melt=reference_melt(),
    )
end

"""
    reference_values(ref, nx) -> NamedTuple

Return `(; centre_tmax, melt_radius)` of a 1-D result `ref`: the peak temperature [K] of
the cell at the centre `dx / 4` of the first radial bin of an `nx` grid, and
`melt_radius_1d` [m].

# Throws

    - `DomainError`: if `nx < 2`
"""
function reference_values(ref, nx::Int)
    nx >= 2 || throw(DomainError(nx, "the grid resolution must be at least 2"))
    nr = length(ref.T_max)
    r_bin = PSV.L_BOX / (nx - 1) / 4
    cell = clamp(floor(Int, r_bin / (PSV.R_PLANET / nr)) + 1, 1, nr)
    return (;
        centre_tmax=ref.T_max[cell], melt_radius=melt_radius_1d(ref.T_max, PSV.R_PLANET)
    )
end

# Value of `nx` in the table `tab` of (grid, value) pairs.
function _grid_value(tab, nx::Int, what::String)
    for (n, v) in tab
        n == nx && return v
    end
    return throw(ArgumentError("no $what for the grid Nx = $nx"))
end

"""
    case_metric_ids(nx) -> Vector{String}

The metric ids of a run on an `nx` grid against the 1-D reference: the centre
temperature `centre_tmax_nx<nx>`, the melt radius `melt_radius_nx<nx>` and the relative L2
error of the hot bins `tmax_rel_l2_hot_nx<nx>`.
"""
function case_metric_ids(nx::Int)
    return ["centre_tmax_nx$nx", "melt_radius_nx$nx", "tmax_rel_l2_hot_nx$nx"]
end

const EXPECTED_METRIC_IDS = (
    "centre_tmax_nx65",
    "melt_radius_nx65",
    "tmax_rel_l2_hot_nx65",
    "energy_lse_residual_nx65",
    "enthalpy_residual_1d",
)

"""
    metric_nx(id) -> Int

The grid of a metric id: the number after `_nx`, or 0 for a metric of the 1-D reference
(`_1d`).

# Throws

    - `ArgumentError`: if `id` names no grid
"""
function metric_nx(id::AbstractString)::Int
    endswith(id, "_1d") && return 0
    m = match(r"_nx(\d+)$", id)
    m === nothing && throw(ArgumentError("metric id '$id' names no grid"))
    return parse(Int, m.captures[1])
end

"""
    metric_units(id) -> String

The units of a metric: `"K"` for a centre temperature, `"m"` for a melt radius and
`"dimensionless"` for an L2 error or a residual.

# Throws

    - `ArgumentError`: if `id` has no known prefix
"""
function metric_units(id::AbstractString)::String
    startswith(id, "centre_tmax_") && return "K"
    startswith(id, "melt_radius_") && return "m"
    any(
        p -> startswith(id, p),
        ("tmax_rel_l2_hot_", "energy_lse_residual_", "enthalpy_residual_"),
    ) && return "dimensionless"
    return throw(ArgumentError("unknown metric id '$id'"))
end

"""
    metric_limit(id) -> NamedTuple

The pass limit `(; atol, rtol, comparison)` of a metric, fixed before any record: the
centre temperature within `CENTRE_TMAX_ATOL` of the 1-D value (20 K at `Nx = 65`, 10 K at
`Nx = 129`), the melt radius within `MELT_RADIUS_ATOL_CELLS` cells (`isapprox`), the L2
error of the hot bins at most `TMAX_REL_L2_LIMIT` (0.04, 0.025), the energy residual at
most `ENERGY_RESIDUAL_LIMIT` and the 1-D enthalpy residual at most
`ENTHALPY_RESIDUAL_LIMIT` (`le_atol`). Every limit has `rtol = 0`.

# Throws

    - `ArgumentError`: if `id` has no known prefix or its grid has no limit
"""
function metric_limit(id::AbstractString)
    T = NamedTuple{(:atol, :rtol, :comparison),Tuple{Float64,Float64,String}}
    if startswith(id, "centre_tmax_")
        return T((_grid_value(CENTRE_TMAX_ATOL, metric_nx(id), id), 0.0, "isapprox"))
    elseif startswith(id, "melt_radius_")
        dx = PSV.L_BOX / (metric_nx(id) - 1)
        return T((MELT_RADIUS_ATOL_CELLS * dx, 0.0, "isapprox"))
    elseif startswith(id, "tmax_rel_l2_hot_")
        return T((_grid_value(TMAX_REL_L2_LIMIT, metric_nx(id), id), 0.0, "le_atol"))
    elseif startswith(id, "energy_lse_residual_")
        return T((ENERGY_RESIDUAL_LIMIT, 0.0, "le_atol"))
    elseif startswith(id, "enthalpy_residual_")
        return T((ENTHALPY_RESIDUAL_LIMIT, 0.0, "le_atol"))
    end
    return throw(ArgumentError("unknown metric id '$id'"))
end

"""
    metric_passes(id, value, reference) -> Bool

Whether `value` passes the limit of `metric_limit(id)` around `reference`.
"""
function metric_passes(id::AbstractString, value::Float64, reference::Float64)::Bool
    lim = metric_limit(id)
    lim.comparison == "isapprox" &&
        return isapprox(value, reference; atol=lim.atol, rtol=lim.rtol)
    return iszero(reference) && isfinite(value) && 0.0 <= value <= lim.atol
end

"""
    metrics_from_run(run, nx) -> Dict{String,Float64}

The metrics of `case_metric_ids(nx)` of a `run_melting_shell_case` result.
"""
function metrics_from_run(run, nx::Int)::Dict{String,Float64}
    a, b, c = case_metric_ids(nx)
    return Dict(a => run.T_centre, b => run.r_melt, c => run.tmax_rel_l2)
end

"""
    evaluate_melting_shell_benchmark(; seed=PSV.SEED_DEFAULT) -> NamedTuple

Run the recorded case at `NX_REFERENCE` against one 1-D reference.

# Returns

`(; runs, metrics, reference)`: `Dict` from `Nx` to the `run_melting_shell_case` result,
the metrics of `EXPECTED_METRIC_IDS` and the 1-D result.
"""
function evaluate_melting_shell_benchmark(; seed::Int=PSV.SEED_DEFAULT)
    ref = reference_solution()
    run = run_melting_shell_case(;
        Nx=NX_REFERENCE,
        phi_init=PHI_INIT_CASE,
        t_form_yr=T_FORM_YR_CASE,
        n_steps=N_STEPS_CASE,
        dt_yr=DT_YR_CASE,
        seed,
        Nr=NR_REFERENCE,
        reference=ref,
    )
    metrics = metrics_from_run(run, NX_REFERENCE)
    metrics["energy_lse_residual_nx$(NX_REFERENCE)"] = run.energy_lse_residual
    metrics["enthalpy_residual_1d"] = ref.enthalpy_residual
    return (; runs=Dict{Int,Any}(NX_REFERENCE => run), metrics, reference=ref)
end

"""
    expected_metric_references(; ref=reference_solution()) -> Dict{String,Float64}

The reference value of each metric of `EXPECTED_METRIC_IDS`: the 1-D centre temperature
and melt radius of `reference_values(ref, NX_REFERENCE)`, and 0 for the others.
"""
function expected_metric_references(; ref=reference_solution())::Dict{String,Float64}
    v = reference_values(ref, NX_REFERENCE)
    refs = Dict{String,Float64}(id => 0.0 for id in EXPECTED_METRIC_IDS)
    refs["centre_tmax_nx$(NX_REFERENCE)"] = v.centre_tmax
    refs["melt_radius_nx$(NX_REFERENCE)"] = v.melt_radius
    return refs
end

"""
    expected_metric_limits() -> Dict{String,NamedTuple}

The `metric_limit` of each metric of `EXPECTED_METRIC_IDS`.
"""
expected_metric_limits() = Dict(id => metric_limit(id) for id in EXPECTED_METRIC_IDS)

# Bands around the recorded value that a new record may not leave without
# --accept-physics-change. The melt radius moves in bins of dx / 2: its band is one bin.
# The centre and L2 bands are 3 times their spread over the marker seeds 42, 1 and 2.
const SEED_BAND_CENTRE_TMAX = 0.382
const SEED_BAND_TMAX_L2 = 5.24e-4
const PHYSICS_REGRESSION_BANDS = Dict{
    String,NamedTuple{(:atol, :rtol),Tuple{Float64,Float64}}
}(
    "centre_tmax_nx65" => (; atol=SEED_BAND_CENTRE_TMAX, rtol=0.0),
    "melt_radius_nx65" => (; atol=375.0, rtol=0.0),
    "tmax_rel_l2_hot_nx65" => (; atol=SEED_BAND_TMAX_L2, rtol=0.0),
)

const REPLAY_RTOL = 1.0e-6
# One marker across a bin edge moves the melt radius by one bin of 375 m.
const REPLAY_ATOL_BOUNDS = Dict{String,Float64}(
    "melt_radius_nx65" => 375.0, "energy_lse_residual_nx65" => 1.0e-9
)

"""
    check_reference_literals_match_config(cfg=make_case_config())
        -> Tuple{Bool,Vector{String}}

Compare the constants of the 1-D reference with the configuration and the Erebus
constants: the checks of `PSV.check_reference_literals_match_config`, then `L_MELT`,
`T_SOL`, `T_LIQ` with the rock entries of `cfg.melting`, a zero Clapeyron slope, the melt
density with `RHO_SOL`, and `RHO_SOL` with the dry rock density `Erebus.MD / Erebus.VDˢ`
that the latent term of a marker uses.

# Returns

Pass flag and list of mismatch descriptions.
"""
function check_reference_literals_match_config(cfg::SimulationConfig=make_case_config())
    _, errors = PSV.check_reference_literals_match_config(cfg)
    m = cfg.melting
    checks = (
        ("L_MELT", L_MELT, m.L_melt),
        ("T_SOL", T_SOL, m.T_solidus[1]),
        ("T_LIQ", T_LIQ, m.T_liquidus[1]),
        ("dpdt_clapeyron", 0.0, m.dpdt_clapeyron),
        ("rho_melt", PSV.RHO_SOL, m.rho_melt),
        ("RHO_SOL (dry rock density)", PSV.RHO_SOL, Erebus.MD / Erebus.VDˢ),
    )
    for (name, lit, val) in checks
        isapprox(lit, val; rtol=1.0e-14, atol=0.0) || push!(
            errors, "Reference literal $name ($lit) does not match config value ($val)"
        )
    end
    return (isempty(errors), errors)
end

const HASHED_CASE_CONSTANTS = (
    :L_MELT,
    :T_SOL,
    :T_LIQ,
    :PHI_INIT_CASE,
    :T_FORM_YR_CASE,
    :N_STEPS_CASE,
    :DT_YR_CASE,
    :NX_REFERENCE,
    :NR_REFERENCE,
    :TMAX_L2_T_MIN,
    :CENTRE_TMAX_ATOL,
    :TMAX_REL_L2_LIMIT,
    :MELT_RADIUS_ATOL_CELLS,
    :ENERGY_RESIDUAL_LIMIT,
    :ENTHALPY_RESIDUAL_LIMIT,
)

# Value of a hashed constant: of this module, else of the porous shell module.
function _constant(name::Symbol)
    return isdefined(@__MODULE__, name) ? getfield(@__MODULE__, name) : getfield(PSV, name)
end

# Lines `prefix.section.field=value` of the hashed sections of `cfg`.
function _print_config_sections(io, cfg, prefix::String)
    for sec_sym in PSV.BENCHMARK_CONFIG_SECTIONS
        sec = getproperty(cfg, sec_sym)
        for f_sym in sort(collect(fieldnames(typeof(sec))))
            val = PSV._canonical_field_repr(getfield(sec, f_sym))
            println(io, prefix, sec_sym, ".", f_sym, "=", val)
        end
    end
    return nothing
end

"""
    canonical_config_text(cfg) -> String

The text that `canonical_config_sha256` hashes: one line `section.field=value` per field
of the sections of `PSV.BENCHMARK_CONFIG_SECTIONS`, fields sorted by name, then one line
`case.NAME=value` per constant of `HASHED_CASE_CONSTANTS` and of
`PSV.REFERENCE_LAW_CONSTANTS`.
"""
function canonical_config_text(cfg)::String
    io = IOBuffer()
    _print_config_sections(io, cfg, "")
    for name in (HASHED_CASE_CONSTANTS..., PSV.REFERENCE_LAW_CONSTANTS...)
        println(io, "case.", name, "=", PSV._canonical_field_repr(_constant(name)))
    end
    return String(take!(io))
end

"""
    canonical_config_sha256(cfg) -> String

The SHA-256 digest (64 lowercase hex characters) of `canonical_config_text(cfg)`.
"""
canonical_config_sha256(cfg)::String = bytes2hex(sha256(canonical_config_text(cfg)))

"""
    expected_config_summary(cfg) -> Dict{String,Any}

The `[config]` section of the record of the case `cfg`, keyed as in
`REQUIRED_CONFIG_KEYS["melting_shell"]`.
"""
function expected_config_summary(cfg::SimulationConfig)::Dict{String,Any}
    return Dict{String,Any}(
        "Nx" => cfg.grid.Nx,
        "Nr_reference" => NR_REFERENCE,
        "n_steps" => cfg.time.n_steps,
        "dt_yr" => cfg.time.dt_initial,
        "t_form_yr" => cfg.time.start_time,
        "phi_init" => cfg.porosity.phi_init,
        "porosity_model" => String(cfg.porosity.model),
        "R_planet_m" => cfg.geometry.rplanet,
        "L_box_m" => cfg.grid.xsize,
        "ratio_al" => cfg.thermodynamics.ratio_al,
        "L_melt_J_kg" => cfg.melting.L_melt,
        "T_solidus_K" => cfg.melting.T_solidus[1],
        "T_liquidus_K" => cfg.melting.T_liquidus[1],
        "rho_melt_kg_m3" => cfg.melting.rho_melt,
        "rhocp_rock_J_m3K" => cfg.materials.rhocpsolidm[1],
        "k_rock_W_mK" => cfg.materials.ksolidm[1],
        "T_space_K" => cfg.disk.t_ambient,
        "surface_temperature" => String(cfg.thermodynamics.surface_temperature),
        "tmax_l2_T_min_K" => TMAX_L2_T_MIN,
        "melt_radius_atol_cells" => MELT_RADIUS_ATOL_CELLS,
        "seed" => cfg.solver.seed,
        "spherical_metric" => cfg.geometry.spherical_metric,
    )
end

# Check of the recorded `[config]` section against `summary`; errors appended to `errors`.
function _check_config_summary!(errors::Vector{String}, rec_cfg, summary)
    PSV._push_extra_config_keys!(errors, rec_cfg, summary)
    for (k, v_exp) in summary
        if !haskey(rec_cfg, k)
            push!(errors, "Recorded config missing key '$k'")
        elseif !PSV._summary_value_matches(rec_cfg[k], v_exp)
            push!(
                errors,
                "Recorded config.$k ($(rec_cfg[k])) does not match expected ($v_exp)",
            )
        end
    end
    return errors
end

# Check of the recorded tolerance of metric `mid` against `metric_limit`.
function _check_recorded_limit!(errors::Vector{String}, mid::String, m)
    lim = metric_limit(mid)
    atol, rtol, comp = Float64(m["atol"]), Float64(m["rtol"]), String(m["comparison"])
    same =
        comp == lim.comparison &&
        isapprox(atol, lim.atol; rtol=1.0e-12, atol=0.0) &&
        isapprox(rtol, lim.rtol; rtol=1.0e-12, atol=0.0)
    same || push!(
        errors,
        "Recorded tolerance of '$mid' ($comp, atol=$atol, rtol=$rtol) differs from the limit of the code ($(lim.comparison), atol=$(lim.atol), rtol=$(lim.rtol))",
    )
    return errors
end

"""
    check_physics_regression_bands(metrics, record_path) -> Tuple{Bool,Vector{String}}

Check that each metric of `PHYSICS_REGRESSION_BANDS` lies within its band around the value
recorded in `record_path`.

# Throws

    - `DomainError`: if `record_path` is empty
"""
function check_physics_regression_bands(metrics::Dict{String,Float64}, record_path::String)
    isempty(record_path) && throw(DomainError(record_path, "Record path must be non-empty"))
    isfile(record_path) || return (false, ["Validation record not found: $record_path"])
    raw = PSV._parse_record(record_path)
    raw isa String && return (false, [raw])
    recorded = Dict{String,Float64}()
    for m in get(raw, "metrics", Any[])
        if m isa AbstractDict &&
            get(m, "id", nothing) isa AbstractString &&
            get(m, "value", nothing) isa Real
            recorded[String(m["id"])] = Float64(m["value"])
        end
    end
    errors = String[]
    for mid in sort!(collect(keys(PHYSICS_REGRESSION_BANDS)))
        band = PHYSICS_REGRESSION_BANDS[mid]
        if !haskey(recorded, mid)
            push!(
                errors,
                "Validation record missing metric '$mid' for physics regression band",
            )
        elseif !haskey(metrics, mid)
            push!(errors, "Measured metrics missing '$mid' for physics regression band")
        elseif !isapprox(metrics[mid], recorded[mid]; atol=band.atol, rtol=band.rtol)
            push!(
                errors,
                "Metric '$mid' = $(metrics[mid]) outside physics regression band around recorded $(recorded[mid]) (atol=$(band.atol), rtol=$(band.rtol))",
            )
        end
    end
    return (isempty(errors), errors)
end

"""
    check_validation_reproducibility(bench, record_path; check_recorded_values=true)
        -> Tuple{Bool,Vector{String}}

Compare a fresh `evaluate_melting_shell_benchmark` result `bench` with the record
`record_path`: the schema, the reference literals, the configuration hash and summary, the
metric ids, the `nx` (`metric_nx`), case and units of each metric, its recorded tolerance
against `metric_limit`, its reference against `expected_metric_references` of the 1-D
result of `bench`, that the measured value passes, the regression bands and, when
`check_recorded_values` holds and the Julia major.minor version matches, the replay of the
recorded values.

# Throws

    - `DomainError`: if `record_path` is empty
"""
function check_validation_reproducibility(
    bench, record_path::String; check_recorded_values::Bool=true
)
    isempty(record_path) && throw(DomainError(record_path, "Record path must be non-empty"))
    isfile(record_path) || return (false, ["Validation record not found: $record_path"])
    raw = PSV._parse_record(record_path)
    raw isa String && return (false, [raw])
    ok_schema, errs = PSV.validate_record_dict(
        raw; require_doc_pages=false, check_doc_sync=false
    )
    ok_schema || return (false, errs)
    errors = String[]
    cfg = make_case_config()
    append!(errors, last(check_reference_literals_match_config(cfg)))
    meta = raw["metadata"]
    String(meta["id"]) == RECORD_ID ||
        push!(errors, "Record id '$(meta["id"])' is not '$RECORD_ID'")
    sha = canonical_config_sha256(cfg)
    String(meta["config_sha256"]) == sha || push!(
        errors,
        "Recorded config_sha256 ($(meta["config_sha256"])) does not match canonical config SHA-256 ($sha)",
    )
    replay =
        check_recorded_values &&
        PSV._same_julia_major_minor(string(get(meta, "julia_version", "")))
    _check_config_summary!(errors, raw["config"], expected_config_summary(cfg))
    ids = [String(m["id"]) for m in raw["metrics"]]
    ids == collect(EXPECTED_METRIC_IDS) || push!(
        errors,
        "Recorded metric IDs $ids do not match expected metric IDs $(collect(EXPECTED_METRIC_IDS))",
    )
    refs = expected_metric_references(; ref=bench.reference)
    for m in raw["metrics"]
        mid = String(m["id"])
        mid in EXPECTED_METRIC_IDS || continue
        Int(m["nx"]) == metric_nx(mid) ||
            push!(errors, "Metric '$mid' has nx $(m["nx"]), not $(metric_nx(mid))")
        String(m["case"]) == CASE_ID ||
            push!(errors, "Metric '$mid' has case '$(m["case"])', not '$CASE_ID'")
        String(m["units"]) == metric_units(mid) || push!(
            errors,
            "Metric '$mid' has units '$(m["units"])', not '$(metric_units(mid))'",
        )
        _check_recorded_limit!(errors, mid, m)
        ref = Float64(m["reference"])
        isapprox(ref, refs[mid]; rtol=1.0e-12, atol=1.0e-12) || push!(
            errors,
            "Recorded reference for '$mid' ($ref) does not match recomputed reference ($(refs[mid]))",
        )
        if !haskey(bench.metrics, mid)
            push!(errors, "Benchmark output missing metric '$mid'")
            continue
        end
        measured = bench.metrics[mid]
        metric_passes(mid, measured, refs[mid]) || push!(
            errors,
            "Measured metric '$mid' = $measured fails its limit around $(refs[mid])",
        )
        if replay
            recorded = Float64(m["value"])
            atol = get(REPLAY_ATOL_BOUNDS, mid, 1.0e-12)
            isapprox(measured, recorded; rtol=REPLAY_RTOL, atol) || push!(
                errors,
                "Measured metric '$mid' = $measured differs from recorded value $recorded (rtol=$REPLAY_RTOL, atol=$atol)",
            )
        end
    end
    append!(errors, last(check_physics_regression_bands(bench.metrics, record_path)))
    return (isempty(errors), errors)
end

# Suite record: the case at Nx 129, too long for CI, checked without a rerun.
# 1-D references of the run at Nx 129, from `reference_values(reference_solution(), 129)`.
const SUITE_REFERENCE_VALUES = (
    ("centre_tmax_nx129", 1661.3024326477178), ("melt_radius_nx129", 13130.0)
)
const SUITE_RESULT_HEADER = (
    "tag",
    "Nx",
    "L_melt",
    "dx",
    "T_centre",
    "T_centre_ref",
    "r_melt",
    "r_melt_ref",
    "tmax_rel_l2_hot",
    "energy_res",
    "enthalpy_res_1d",
    "wall_s",
)
const SUITE_HASHED_CONSTANTS = (
    :NX_SUITE, :SUITE_REFERENCE_VALUES, HASHED_CASE_CONSTANTS...
)

"""
    suite_metric_ids() -> Vector{String}

The metric ids of the suite record: `case_metric_ids(NX_SUITE)`.
"""
suite_metric_ids() = case_metric_ids(NX_SUITE)

"""
    suite_reference_value(id) -> Float64

The pinned 1-D reference of the suite metric `id`: its entry of `SUITE_REFERENCE_VALUES`,
or 0 for the L2 error.

# Throws

    - `ArgumentError`: if `id` is not a suite metric
"""
function suite_reference_value(id::AbstractString)::Float64
    id in suite_metric_ids() || throw(ArgumentError("'$id' is not a suite metric"))
    for (k, v) in SUITE_REFERENCE_VALUES
        k == id && return v
    end
    return 0.0
end

"""
    suite_config_text() -> String

The text that `suite_config_sha256` hashes: the configuration sections of the case at
`NX_SUITE` as `suite.section.field=value`, then `case.NAME=value` per constant of
`SUITE_HASHED_CONSTANTS` and of `PSV.REFERENCE_LAW_CONSTANTS`.
"""
function suite_config_text()::String
    io = IOBuffer()
    _print_config_sections(io, make_case_config(; Nx=NX_SUITE), "suite.")
    for name in (SUITE_HASHED_CONSTANTS..., PSV.REFERENCE_LAW_CONSTANTS...)
        println(io, "case.", name, "=", PSV._canonical_field_repr(_constant(name)))
    end
    return String(take!(io))
end

"""
    suite_config_sha256() -> String

The SHA-256 digest of `suite_config_text()`.
"""
suite_config_sha256()::String = bytes2hex(sha256(suite_config_text()))

"""
    expected_suite_config_summary() -> Dict{String,Any}

The `[config]` section of the suite record: `expected_config_summary` of the case at
`NX_SUITE`.
"""
expected_suite_config_summary() = expected_config_summary(make_case_config(; Nx=NX_SUITE))

"""
    check_suite_record_static(path; check_doc_sync=true) -> Tuple{Bool,Vector{String}}

Check the suite record `path` without a rerun: the schema, the record id, the
configuration hash (`suite_config_sha256`), the `[config]` summary, the metric ids, the
case, grid and units of each metric, its tolerance against `metric_limit`, its reference
against the pinned `suite_reference_value` and that its value passes.

# Throws

    - `DomainError`: if `path` is empty
"""
function check_suite_record_static(path::String; check_doc_sync::Bool=true)
    isempty(path) && throw(DomainError(path, "Record path must be non-empty"))
    isfile(path) || return (false, ["Suite record not found: $path"])
    raw = PSV._parse_record(path)
    raw isa String && return (false, [raw])
    ok_schema, errs = PSV.validate_record_dict(raw; require_doc_pages=true, check_doc_sync)
    ok_schema || return (false, errs)
    errors = String[]
    meta = raw["metadata"]
    String(meta["id"]) == SUITE_RECORD_ID ||
        push!(errors, "Record id '$(meta["id"])' is not '$SUITE_RECORD_ID'")
    sha = suite_config_sha256()
    String(meta["config_sha256"]) == sha || push!(
        errors,
        "Recorded config_sha256 ($(meta["config_sha256"])) does not match the case ($sha)",
    )
    _check_config_summary!(errors, raw["config"], expected_suite_config_summary())
    ids = [String(m["id"]) for m in raw["metrics"]]
    ids == suite_metric_ids() ||
        push!(errors, "Recorded metric ids $ids do not match $(suite_metric_ids())")
    for m in raw["metrics"]
        mid = String(m["id"])
        mid in suite_metric_ids() || continue
        String(m["case"]) == CASE_ID ||
            push!(errors, "Metric '$mid' has case '$(m["case"])', not '$CASE_ID'")
        Int(m["nx"]) == NX_SUITE ||
            push!(errors, "Metric '$mid' has nx $(m["nx"]), not $NX_SUITE")
        String(m["units"]) == metric_units(mid) || push!(
            errors,
            "Metric '$mid' has units '$(m["units"])', not '$(metric_units(mid))'",
        )
        _check_recorded_limit!(errors, mid, m)
        value, ref, pinned = Float64(m["value"]),
        Float64(m["reference"]),
        suite_reference_value(mid)
        isapprox(ref, pinned; rtol=1.0e-12, atol=0.0) || push!(
            errors,
            "Recorded reference of '$mid' ($ref) differs from the pinned 1-D value ($pinned)",
        )
        metric_passes(mid, value, pinned) ||
            push!(errors, "Metric '$mid' = $value fails its limit around $pinned")
    end
    return (isempty(errors), errors)
end

end # module MeltingShellValidation

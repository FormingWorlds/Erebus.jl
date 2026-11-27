module MagmaOceanConductivityValidation

# Effective thermal conductivity of convecting melt in a porous planetesimal: dry porosity
# laws of Henke et al. (2012), 26Al heating, silicate melting with latent heat, and
# regularized soft-turbulence conductivity enhancement above melt fraction 0.30.

using Erebus
using StaticArrays
using Logging
using SHA

if isdefined(parentmodule(@__MODULE__), :MeltingShellValidation)
    const MSV = getfield(parentmodule(@__MODULE__), :MeltingShellValidation)
else
    include(joinpath(@__DIR__, "validation_melting_shell.jl"))
    const MSV = MeltingShellValidation
end
const PSV = MSV.PSV
const validate_record_dict = PSV.validate_record_dict

const L_MELT = MSV.L_MELT
const T_SOL = MSV.T_SOL
const T_LIQ = MSV.T_LIQ

const ETA_SOLID_REF = 1.0e19
const ALPHA_ETA_REF = 28.0
const PHI_CRIT_REF = 0.4
const ETA_MELT_REF = 10.0
const ETA_MIN_REF = 1.0e12
const ETA_MAX_REF = 1.0e23
const ETA_FLUID_SIL_REF = 100.0
const F_TURB_START_REF = 0.30
const F_TURB_END_REF = 0.50
const TURB_EXP_REF = 1.0 / 3.0
const DT_TURB_MIN_REF = 10.0
const T_SURF_REF = 300.0
const K_TURB_FLOOR_REF = 1.0e-3
const K_TURB_CUTOFF_REF = 1.0e6
const T_MUSH = T_SOL + F_TURB_START_REF * (T_LIQ - T_SOL)

const REFERENCE_TURB_CONSTANTS = (
    :ETA_SOLID_REF,
    :ALPHA_ETA_REF,
    :PHI_CRIT_REF,
    :ETA_MELT_REF,
    :ETA_MIN_REF,
    :ETA_MAX_REF,
    :ETA_FLUID_SIL_REF,
    :F_TURB_START_REF,
    :F_TURB_END_REF,
    :TURB_EXP_REF,
    :DT_TURB_MIN_REF,
    :T_SURF_REF,
    :K_TURB_FLOOR_REF,
    :K_TURB_CUTOFF_REF,
)

"""
    reference_soft_turb(; enabled=true, T_sol=T_SOL, T_liq=T_LIQ,
                        turb_exponent=TURB_EXP_REF, eta_fluid_sil=ETA_FLUID_SIL_REF,
                        face_mode=:harmonic) -> NamedTuple

Parameters of the 1-D soft-turbulence reference `PSV.k_soft_turb_ref` and
`PSV.solve_porous_sphere`.
"""
function reference_soft_turb(;
    enabled::Bool=true,
    T_sol::Float64=T_SOL,
    T_liq::Float64=T_LIQ,
    eta_solid::Float64=ETA_SOLID_REF,
    alpha_eta::Float64=ALPHA_ETA_REF,
    phi_crit::Float64=PHI_CRIT_REF,
    eta_melt::Float64=ETA_MELT_REF,
    eta_min::Float64=ETA_MIN_REF,
    eta_max::Float64=ETA_MAX_REF,
    eta_fluid_sil::Float64=ETA_FLUID_SIL_REF,
    F_start::Float64=F_TURB_START_REF,
    F_end::Float64=F_TURB_END_REF,
    turb_exponent::Float64=TURB_EXP_REF,
    dT_min::Float64=DT_TURB_MIN_REF,
    T_surf::Float64=T_SURF_REF,
    k_floor::Float64=K_TURB_FLOOR_REF,
    k_cutoff::Float64=K_TURB_CUTOFF_REF,
    face_mode::Symbol=:harmonic,
)
    return (;
        enabled,
        T_sol,
        T_liq,
        eta_solid,
        alpha_eta,
        phi_crit,
        eta_melt,
        eta_min,
        eta_max,
        eta_fluid_sil,
        F_start,
        F_end,
        turb_exponent,
        dT_min,
        T_surf,
        k_floor,
        k_cutoff,
        face_mode,
    )
end

"""
    make_magma_ocean_config(; Nx, phi_init, t_form_yr, n_steps, dt_yr,
                            soft_turbulence=true, turb_exponent=TURB_EXP_REF,
                            eta_fluid_silicate=ETA_FLUID_SIL_REF, L_melt=L_MELT,
                            T_sol=T_SOL, T_liq=T_LIQ, kwargs...) -> SimulationConfig

Configuration of the porous shell case with silicate melting, latent heat, melt density
equal to the solid rock density, and regularized soft-turbulence conductivity.

# Throws

    - `ArgumentError`: if `T_sol >= T_liq`, or the errors of `make_porous_shell_config`
    - `DomainError`: for out-of-range arguments in `make_porous_shell_config`
"""
function make_magma_ocean_config(;
    soft_turbulence::Bool=true,
    turb_exponent::Float64=TURB_EXP_REF,
    eta_fluid_silicate::Float64=ETA_FLUID_SIL_REF,
    L_melt::Float64=L_MELT,
    T_sol::Float64=T_SOL,
    T_liq::Float64=T_LIQ,
    kwargs...,
)
    T_sol < T_liq || throw(ArgumentError("T_sol must be below T_liq, got $T_sol, $T_liq"))
    melting = MeltingConfig(;
        active=true,
        soft_turbulence=soft_turbulence,
        turb_exponent=turb_exponent,
        eta_fluid_silicate=eta_fluid_silicate,
        L_melt=L_melt,
        T_solidus=SVector{3,Float64}([T_sol, T_sol, NaN]),
        T_liquidus=SVector{3,Float64}([T_liq, T_liq, NaN]),
        rho_melt=PSV.RHO_SOL,
    )
    return PSV.make_porous_shell_config(; melting=melting, kwargs...)
end

"""
    mush_radius(xm, ym, T_max; xc, yc, R_p, dr, T_mush=T_MUSH, nsector=16) -> Float64

Radius [m] of the convectively enhanced mush core where marker peak temperatures `T_max`
[K] exceed `T_mush` (`F_max > F_turb_start`), measured across `nsector` sectors in bins of
`dr` via `MeltingShellValidation.melt_radius`.
"""
function mush_radius(
    xm, ym, T_max; xc, yc, R_p, dr, T_mush::Float64=T_MUSH, nsector::Int=16
)
    return MSV.melt_radius(xm, ym, T_max; xc, yc, R_p, dr, T_sol=T_mush, nsector)
end

"""
    mush_radius_1d(T_max, R; T_mush=T_MUSH) -> Float64

Radius [m] of the convectively enhanced mush core (`T_max > T_mush`) of a 1-D profile on a
sphere of radius `R` [m].
"""
function mush_radius_1d(T_max::AbstractVector{Float64}, R::Float64; T_mush::Float64=T_MUSH)
    return MSV.melt_radius_1d(T_max, R; T_sol=T_mush)
end

"""
    plateau_drop(rc, tmax; r_outer=0.5 * PSV.R_PLANET) -> Float64

Temperature drop `tmax[1] - tmax[idx]` [K] across the inner convecting core from the
innermost radial bin `rc[1]` to the bin closest to `r_outer` [m].

# Throws

    - `ArgumentError`: if `rc` and `tmax` differ in length, have fewer than 2 entries, or
      `rc` is not strictly increasing
    - `DomainError`: if `r_outer` is not finite or lies outside `(rc[1], rc[end]]`
"""
function plateau_drop(
    rc::AbstractVector{Float64},
    tmax::AbstractVector{Float64};
    r_outer::Float64=0.5 * PSV.R_PLANET,
)::Float64
    length(rc) == length(tmax) ||
        throw(ArgumentError("rc and tmax must have the same length"))
    length(rc) >= 2 || throw(ArgumentError("rc must contain at least 2 radial bins"))
    issorted(rc; lt=(<=)) && all(diff(rc) .> 0.0) ||
        throw(ArgumentError("rc must be strictly increasing"))
    (isfinite(r_outer) && rc[1] < r_outer <= rc[end]) ||
        throw(DomainError(r_outer, "r_outer must be finite and in ($(rc[1]), $(rc[end])]"))
    idx = argmin(abs.(rc .- r_outer))
    return tmax[1] - tmax[idx]
end

"""
    hot_core_metric_step_check(; Nx=41, r_core=12_000.0, T_core=1850.0, dt_yr=2500.0,
                               Nr=2000) -> NamedTuple

Single-step comparison on an unrelaxed sintered (`phi_d = 0`) molten core (`T = T_core`
inside `r_core`, linear ramp to `T_SPACE` at `R_PLANET`) over one step `dt_yr` [yr]
without radiogenic heating, quantifying the step-1 explicit `Q_metric` cooling against the
2-D solve with `Q_metric = 0` and the 1-D implicit spherical step.

# Returns

`(; k_core_max, min_dT_metric, T_centre_2d_sph, T_centre_2d_no_metric, T_centre_1d)`.
"""
function hot_core_metric_step_check(;
    Nx::Int=41,
    r_core::Float64=12_000.0,
    T_core::Float64=1850.0,
    dt_yr::Float64=2500.0,
    Nr::Int=2000,
)
    R = PSV.R_PLANET
    L = PSV.L_BOX
    dt_s = dt_yr * PSV.YEARLENGTH_S
    T_prof(r) =
        if r <= r_core
            T_core
        elseif r >= R
            PSV.T_SPACE
        else
            T_core + (PSV.T_SPACE - T_core) * (r - r_core) / (R - r_core)
        end
    cfg = make_magma_ocean_config(; Nx, phi_init=0.0, t_form_yr=0.85e6, n_steps=1, dt_yr)
    st, coords, ws, tel, _ = Erebus.init_simulation(cfg; output_path=mktempdir())
    tel !== nothing && close(tel)
    xc = yc = L / 2.0
    for m in eachindex(st.markers.core.tm)
        if st.markers.core.tm[m] < 3
            rm = hypot(st.markers.core.xm[m] - xc, st.markers.core.ym[m] - yc)
            st.markers.core.tkm[m] = T_prof(rm)
            st.markers.groups.porosity.phi_d[m] = 0.0
        end
    end
    Erebus.interpolate_markers_to_grid!(
        st,
        coords,
        cfg;
        thread_buffers=ws.thread_buffers,
        interp_arrays=ws.interp_arrays,
        rock_sums=ws.rock_sums,
    )
    g = st.grids
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    Z = zeros(Ny1, Nx1)
    Q = zeros(Ny1, Nx1)
    RT_sph = zeros(Ny1 * Nx1)
    RT_no = zeros(Ny1 * Nx1)
    Erebus.compute_spherical_metric_heat_source!(
        Q,
        g.tk1,
        ws.rock_sums.KXu,
        ws.rock_sums.KYu,
        coords;
        xcenter=xc,
        ycenter=yc,
        rplanet=R,
        reg_cells=cfg.geometry.metric_regularization_cells,
    )
    LT_sph = Erebus.assemble_thermal_lse!(
        g.tk1,
        g.RHOCP,
        g.KX,
        g.KY,
        Z,
        Z,
        Z,
        Z,
        RT_sph,
        dt_s;
        coords=coords,
        Q_metric=Q,
        fixed=g.fixed_tk,
    )
    tk_sph = reshape(LT_sph.cscmatrix \ RT_sph, Ny1, Nx1)
    LT_no = Erebus.assemble_thermal_lse!(
        g.tk1,
        g.RHOCP,
        g.KX,
        g.KY,
        Z,
        Z,
        Z,
        Z,
        RT_no,
        dt_s;
        coords=coords,
        Q_metric=Z,
        fixed=g.fixed_tk,
    )
    tk_no = reshape(LT_no.cscmatrix \ RT_no, Ny1, Nx1)
    ic = (Nx + 1) ÷ 2
    min_dT_qm = minimum([
        Q[i, j] * dt_s / g.RHOCP[i, j] for
        i in 2:(Ny1 - 1), j in 2:(Nx1 - 1) if !g.fixed_tk[i, j]
    ])

    dr = R / Nr
    rf = [(i - 1) * dr for i in 1:(Nr + 1)]
    rc = [(i - 0.5) * dr for i in 1:Nr]
    A = 4pi .* rf .^ 2
    V = [(4pi / 3) * (rf[i + 1]^3 - rf[i]^3) for i in 1:Nr]
    T0v = [T_prof(r) for r in rc]
    st_nt = reference_soft_turb()
    kv = [PSV.k_soft_turb_ref(PSV.K_SOL, t; soft_turb=st_nt) for t in T0v]
    Cv = [PSV.rhocp_app_ref(t, 0.0; L_melt=L_MELT, T_sol=T_SOL, T_liq=T_LIQ) for t in T0v]
    cond = zeros(Nr + 1)
    for i in 2:Nr
        cond[i] = A[i] / (0.5 * dr / kv[i - 1] + 0.5 * dr / kv[i])
    end
    cond[Nr + 1] = A[Nr + 1] * kv[Nr] / (0.5 * dr)
    d = [Cv[i] * V[i] / dt_s + cond[i] + cond[i + 1] for i in 1:Nr]
    lower = [-cond[i] for i in 2:Nr]
    rhs = [Cv[i] * V[i] / dt_s * T0v[i] for i in 1:Nr]
    rhs[Nr] += cond[Nr + 1] * PSV.T_SPACE
    T1v = PSV._thomas(lower, d, lower, rhs)

    return (;
        k_core_max=maximum(ws.rock_sums.KXu),
        min_dT_metric=min_dT_qm,
        T_centre_2d_sph=tk_sph[ic, ic],
        T_centre_2d_no_metric=tk_no[ic, ic],
        T_centre_1d=T1v[1],
    )
end

# Recorded case: R_p = 20 km, phi_init = 0.25, formation at 0.85 Myr, 140 steps of 10 kyr.
const CASE_ID = "magma_ocean_t0.85_phi0.25"
const RECORD_ID = "magma_ocean_conductivity"
const SUITE_RECORD_ID = "magma_ocean_conductivity_suite"
const PHI_INIT_CASE = 0.25
const T_FORM_YR_CASE = 0.85e6
const N_STEPS_CASE = 140
const DT_YR_CASE = 1.0e4
const NX_REFERENCE = 65
const NX_SUITE = 129
const NR_REFERENCE = 2000
const TMAX_L2_T_MIN = 1100.0
const CENTRE_TMAX_ATOL = ((65, 8.0), (129, 5.0))
const PLATEAU_DROP_LIMIT = ((65, 8.0), (129, 6.0))
const TMAX_REL_L2_LIMIT = ((65, 0.04), (129, 0.025))
const RADIUS_ATOL_CELLS = 1.5
const ENERGY_RESIDUAL_LIMIT = 1.0e-9
const ENTHALPY_RESIDUAL_LIMIT = 0.01
const LAG_DT_DIVISOR = 10
const TURB_EXP_LOW_MUTANT = 0.2
const TURB_EXP_HIGH_MUTANT = 0.5

"""
    make_case_config(; Nx=NX_REFERENCE, seed=PSV.SEED_DEFAULT, soft_turbulence=true,
                     turb_exponent=TURB_EXP_REF) -> SimulationConfig

Return the `SimulationConfig` of the recorded magma-ocean conductivity case on an `Nx` x `Nx` grid.
"""
function make_case_config(;
    Nx::Int=NX_REFERENCE,
    seed::Int=PSV.SEED_DEFAULT,
    soft_turbulence::Bool=true,
    turb_exponent::Float64=TURB_EXP_REF,
)
    return make_magma_ocean_config(;
        Nx,
        phi_init=PHI_INIT_CASE,
        t_form_yr=T_FORM_YR_CASE,
        n_steps=N_STEPS_CASE,
        dt_yr=DT_YR_CASE,
        soft_turbulence,
        turb_exponent,
        seed,
    )
end

"""
    reference_solution(; dt_yr=DT_YR_CASE, Nr=NR_REFERENCE, face_mode=:harmonic) -> NamedTuple

Return the 1-D spherical reference of the recorded case with nominal latent heat and
soft-turbulence conductivity over `N_STEPS_CASE * DT_YR_CASE`.
"""
function reference_solution(;
    dt_yr::Float64=DT_YR_CASE, Nr::Int=NR_REFERENCE, face_mode::Symbol=:harmonic
)
    return PSV.solve_porous_sphere(;
        R=PSV.R_PLANET,
        phi_init=PHI_INIT_CASE,
        t_form_yr=T_FORM_YR_CASE,
        t_end_yr=N_STEPS_CASE * DT_YR_CASE,
        dt_yr,
        Nr,
        stop_cold=false,
        melt=MSV.reference_melt(),
        soft_turb=reference_soft_turb(; face_mode),
    )
end

"""
    reference_values(ref, nx; r_outer=0.5 * PSV.R_PLANET) -> NamedTuple

Return `(; centre_tmax, plateau_drop, mush_radius, melt_radius)` of a 1-D result `ref`
sampled at the radial bin centres of an `nx` grid (`dx / 4` at the centre and the bin
containing `r_outer` for the plateau drop).

# Throws

    - `DomainError`: if `nx < 2`
"""
function reference_values(ref, nx::Int; r_outer::Float64=0.5 * PSV.R_PLANET)
    nx >= 2 || throw(DomainError(nx, "the grid resolution must be at least 2"))
    nr = length(ref.T_max)
    dx = PSV.L_BOX / (nx - 1)
    dr_bin = dx / 2
    nb = ceil(Int, PSV.R_PLANET / dr_bin)
    r_bin1 = 0.5 * min(dr_bin, PSV.R_PLANET)
    b_half = min(floor(Int, r_outer / dr_bin) + 1, nb)
    r_half = 0.5 * ((b_half - 1) * dr_bin + min(b_half * dr_bin, PSV.R_PLANET))
    c1 = clamp(floor(Int, r_bin1 / (PSV.R_PLANET / nr)) + 1, 1, nr)
    ch = clamp(floor(Int, r_half / (PSV.R_PLANET / nr)) + 1, 1, nr)
    return (;
        centre_tmax=ref.T_max[c1],
        plateau_drop=ref.T_max[c1] - ref.T_max[ch],
        mush_radius=mush_radius_1d(ref.T_max, PSV.R_PLANET),
        melt_radius=MSV.melt_radius_1d(ref.T_max, PSV.R_PLANET),
    )
end

"""
    run_magma_ocean_case(; Nx, phi_init, t_form_yr, n_steps, dt_yr,
                         soft_turbulence=true, turb_exponent=TURB_EXP_REF,
                         eta_fluid_silicate=ETA_FLUID_SIL_REF, seed=PSV.SEED_DEFAULT,
                         Nr=NR_REFERENCE, reference=nothing,
                         run_loop=simulation_loop) -> NamedTuple

Run the 2-D `Erebus.jl` case and compare it against the nominal 1-D spherical reference at
the same time step and run duration.

# Returns

`(; dx, T_centre, T_centre_ref, T_centre_error, plateau_drop, plateau_drop_ref, r_mush,
r_mush_ref, r_mush_error_cells, r_melt, r_melt_ref, r_melt_error_cells, tmax_rel_l2,
energy_lse_residual, enthalpy_residual_ref, T_peak, r, tmax, tmax_ref, reference)`.
"""
function run_magma_ocean_case(;
    Nx::Int,
    phi_init::Float64,
    t_form_yr::Float64,
    n_steps::Int,
    dt_yr::Float64,
    soft_turbulence::Bool=true,
    turb_exponent::Float64=TURB_EXP_REF,
    eta_fluid_silicate::Float64=ETA_FLUID_SIL_REF,
    seed::Int=PSV.SEED_DEFAULT,
    Nr::Int=NR_REFERENCE,
    reference=nothing,
    run_loop=simulation_loop,
)
    cfg = make_magma_ocean_config(;
        Nx,
        phi_init,
        t_form_yr,
        n_steps,
        dt_yr,
        soft_turbulence,
        turb_exponent,
        eta_fluid_silicate,
        seed,
    )
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
            melt=MSV.reference_melt(),
            soft_turb=reference_soft_turb(),
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
    pdrop = plateau_drop(rc, tmax)
    rv = reference_values(ref, Nx)
    r_mush_val = mush_radius(
        xm[rock], ym[rock], pg.tkm_max[rock]; xc, yc, R_p=PSV.R_PLANET, dr=dx / 2
    )
    r_melt_val = MSV.melt_radius(
        xm[rock], ym[rock], pg.tkm_max[rock]; xc, yc, R_p=PSV.R_PLANET, dr=dx / 2
    )
    return (;
        dx,
        T_centre=tmax[1],
        T_centre_ref=rv.centre_tmax,
        T_centre_error=tmax[1] - rv.centre_tmax,
        plateau_drop=pdrop,
        plateau_drop_ref=rv.plateau_drop,
        r_mush=r_mush_val,
        r_mush_ref=rv.mush_radius,
        r_mush_error_cells=(r_mush_val - rv.mush_radius) / dx,
        r_melt=r_melt_val,
        r_melt_ref=rv.melt_radius,
        r_melt_error_cells=(r_melt_val - rv.melt_radius) / dx,
        tmax_rel_l2=sqrt(
            sum(abs2, tmax[hot] .- tmax_ref[hot]) / sum(abs2, tmax_ref[hot] .- PSV.T_SPACE)
        ),
        energy_lse_residual=PSV.thermal_step_energy_residual(
            st.grids, st.dt, dx, dx; fixed=st.grids.fixed_tk
        ),
        enthalpy_residual_ref=ref.enthalpy_residual,
        T_peak=maximum(pg.tkm_max[rock]),
        T_min_grid=minimum(st.grids.tk2[2:(end - 1), 2:(end - 1)]),
        r=rc,
        tmax,
        tmax_ref,
        reference=ref,
    )
end

"""
    lag_shift(; dt_yr=DT_YR_CASE, divisor=LAG_DT_DIVISOR, nx=NX_REFERENCE, Nr=NR_REFERENCE,
              coarse=nothing) -> NamedTuple

Shift of the 1-D reference of the recorded case between `dt_yr` and `dt_yr / divisor`
from the step-start evaluation of apparent heat capacity and soft-turbulence conductivity.

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
        dr_mush=a.mush_radius - b.mush_radius,
        dr_melt=a.melt_radius - b.melt_radius,
        coarse,
        fine,
    )
end

"""
    face_averaging_shift(; nx=NX_REFERENCE, Nr=NR_REFERENCE, reference=nothing) -> NamedTuple

Difference in the 1-D reference between arithmetic and harmonic cell-face conductivity
averaging (`arithmetic - harmonic`).
"""
function face_averaging_shift(;
    nx::Int=NX_REFERENCE, Nr::Int=NR_REFERENCE, reference=nothing
)
    harm = reference === nothing ? reference_solution(; Nr) : reference
    arith = reference_solution(; Nr, face_mode=:arithmetic)
    vh, va = reference_values(harm, nx), reference_values(arith, nx)
    return (;
        dT_centre=va.centre_tmax - vh.centre_tmax,
        dr_mush=va.mush_radius - vh.mush_radius,
        dr_melt=va.melt_radius - vh.melt_radius,
        harmonic=harm,
        arithmetic=arith,
    )
end

"""
    case_metric_ids(nx) -> Vector{String}

Metric ids of a run on an `nx` grid: centre peak temperature, inner plateau drop, mush
radius, melt radius, and relative L2 error of the hot bins.
"""
function case_metric_ids(nx::Int)
    return [
        "centre_tmax_nx$nx",
        "plateau_drop_nx$nx",
        "mush_radius_nx$nx",
        "melt_radius_nx$nx",
        "tmax_rel_l2_hot_nx$nx",
    ]
end

const EXPECTED_METRIC_IDS = (
    "centre_tmax_nx65",
    "plateau_drop_nx65",
    "mush_radius_nx65",
    "melt_radius_nx65",
    "tmax_rel_l2_hot_nx65",
    "energy_lse_residual_nx65",
    "enthalpy_residual_1d",
)

"""
    metric_nx(id) -> Int

Grid resolution parsed from a metric id (`0` for `_1d`).
"""
metric_nx(id::AbstractString)::Int = MSV.metric_nx(id)

"""
    metric_units(id) -> String

Units of a metric id (`"K"`, `"m"`, or `"dimensionless"`).

# Throws

    - `ArgumentError`: if `id` has no known prefix
"""
function metric_units(id::AbstractString)::String
    (startswith(id, "centre_tmax_") || startswith(id, "plateau_drop_")) && return "K"
    (startswith(id, "mush_radius_") || startswith(id, "melt_radius_")) && return "m"
    any(
        p -> startswith(id, p),
        ("tmax_rel_l2_hot_", "energy_lse_residual_", "enthalpy_residual_"),
    ) && return "dimensionless"
    return throw(ArgumentError("unknown metric id '$id'"))
end

"""
    metric_limit(id) -> NamedTuple

Pass limit `(; atol, rtol, comparison)` of a metric id, fixed before any record.

# Throws

    - `ArgumentError`: if `id` has no known prefix or its grid has no limit
"""
function metric_limit(id::AbstractString)
    T = NamedTuple{(:atol, :rtol, :comparison),Tuple{Float64,Float64,String}}
    if startswith(id, "centre_tmax_")
        return T((MSV._grid_value(CENTRE_TMAX_ATOL, metric_nx(id), id), 0.0, "isapprox"))
    elseif startswith(id, "plateau_drop_")
        return T((MSV._grid_value(PLATEAU_DROP_LIMIT, metric_nx(id), id), 0.0, "le_atol"))
    elseif startswith(id, "mush_radius_") || startswith(id, "melt_radius_")
        dx = PSV.L_BOX / (metric_nx(id) - 1)
        return T((RADIUS_ATOL_CELLS * dx, 0.0, "isapprox"))
    elseif startswith(id, "tmax_rel_l2_hot_")
        return T((MSV._grid_value(TMAX_REL_L2_LIMIT, metric_nx(id), id), 0.0, "le_atol"))
    elseif startswith(id, "energy_lse_residual_")
        return T((ENERGY_RESIDUAL_LIMIT, 0.0, "le_atol"))
    elseif startswith(id, "enthalpy_residual_")
        return T((ENTHALPY_RESIDUAL_LIMIT, 0.0, "le_atol"))
    end
    return throw(ArgumentError("unknown metric id '$id'"))
end

"""
    metric_passes(id, value, reference) -> Bool

Whether `value` satisfies `metric_limit(id)` around `reference`.
"""
function metric_passes(id::AbstractString, value::Float64, reference::Float64)::Bool
    lim = metric_limit(id)
    lim.comparison == "isapprox" &&
        return isapprox(value, reference; atol=lim.atol, rtol=lim.rtol)
    return iszero(reference) && isfinite(value) && 0.0 <= value <= lim.atol
end

"""
    metrics_from_run(run, nx) -> Dict{String,Float64}

Extract the 5 `case_metric_ids(nx)` values from a `run_magma_ocean_case` result.
"""
function metrics_from_run(run, nx::Int)::Dict{String,Float64}
    c, p, mu, me, l2 = case_metric_ids(nx)
    return Dict(
        c => run.T_centre,
        p => run.plateau_drop,
        mu => run.r_mush,
        me => run.r_melt,
        l2 => run.tmax_rel_l2,
    )
end

"""
    evaluate_magma_ocean_benchmark(; seed=PSV.SEED_DEFAULT) -> NamedTuple

Run the recorded case at `NX_REFERENCE` against the 1-D reference.
"""
function evaluate_magma_ocean_benchmark(; seed::Int=PSV.SEED_DEFAULT)
    ref = reference_solution()
    run = run_magma_ocean_case(;
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

Reference values of `EXPECTED_METRIC_IDS` (`0.0` for `le_atol` metrics).
"""
function expected_metric_references(; ref=reference_solution())::Dict{String,Float64}
    v = reference_values(ref, NX_REFERENCE)
    refs = Dict{String,Float64}(id => 0.0 for id in EXPECTED_METRIC_IDS)
    refs["centre_tmax_nx$(NX_REFERENCE)"] = v.centre_tmax
    refs["mush_radius_nx$(NX_REFERENCE)"] = v.mush_radius
    refs["melt_radius_nx$(NX_REFERENCE)"] = v.melt_radius
    return refs
end

"""
    expected_metric_limits() -> Dict{String,NamedTuple}

Pass limits of `EXPECTED_METRIC_IDS`.
"""
expected_metric_limits() = Dict(id => metric_limit(id) for id in EXPECTED_METRIC_IDS)

# Regression bands: 3x spread across marker seeds 42, 1 and 2 at Nx = 65; 1 bin (375 m) on radii.
const SEED_BAND_CENTRE_TMAX = 0.00906
const SEED_BAND_PLATEAU_DROP = 0.0237
const SEED_BAND_TMAX_L2 = 8.56e-4
const PHYSICS_REGRESSION_BANDS = Dict{
    String,NamedTuple{(:atol, :rtol),Tuple{Float64,Float64}}
}(
    "centre_tmax_nx65" => (; atol=SEED_BAND_CENTRE_TMAX, rtol=0.0),
    "plateau_drop_nx65" => (; atol=SEED_BAND_PLATEAU_DROP, rtol=0.0),
    "mush_radius_nx65" => (; atol=375.0, rtol=0.0),
    "melt_radius_nx65" => (; atol=375.0, rtol=0.0),
    "tmax_rel_l2_hot_nx65" => (; atol=SEED_BAND_TMAX_L2, rtol=0.0),
)

const REPLAY_RTOL = 1.0e-6
const REPLAY_ATOL_BOUNDS = Dict{String,Float64}(
    "mush_radius_nx65" => 375.0,
    "melt_radius_nx65" => 375.0,
    "energy_lse_residual_nx65" => 1.0e-9,
)

"""
    check_reference_literals_match_config(cfg=make_case_config()) -> Tuple{Bool,Vector{String}}

Check that all constants of the 1-D porous, melting and soft-turbulence reference laws
match the corresponding fields of `cfg`.
"""
function check_reference_literals_match_config(cfg::SimulationConfig=make_case_config())
    _, errors = MSV.check_reference_literals_match_config(cfg)
    m = cfg.melting
    mat = cfg.materials
    sol = cfg.solver
    th = cfg.thermodynamics
    checks = (
        ("ETA_SOLID_REF (etasolidmm[1])", ETA_SOLID_REF, mat.etasolidmm[1]),
        ("ETA_SOLID_REF (etasolidmm[2])", ETA_SOLID_REF, mat.etasolidmm[2]),
        ("ETA_SOLID_REF (etasolidm[1])", ETA_SOLID_REF, mat.etasolidm[1]),
        ("ETA_SOLID_REF (etasolidm[2])", ETA_SOLID_REF, mat.etasolidm[2]),
        ("T_SOL (tmsolidphase)", T_SOL, th.tmsolidphase),
        ("T_SOL (T_solidus[2])", T_SOL, m.T_solidus[2]),
        ("T_LIQ (T_liquidus[2])", T_LIQ, m.T_liquidus[2]),
        ("ETA_MIN_REF", ETA_MIN_REF, sol.etamin),
        ("ETA_MAX_REF", ETA_MAX_REF, sol.etamax),
        ("ALPHA_ETA_REF", ALPHA_ETA_REF, m.alpha_eta),
        ("PHI_CRIT_REF", PHI_CRIT_REF, m.phi_crit),
        ("ETA_MELT_REF", ETA_MELT_REF, m.eta_melt),
        ("ETA_FLUID_SIL_REF", ETA_FLUID_SIL_REF, m.eta_fluid_silicate),
        ("F_TURB_START_REF", F_TURB_START_REF, m.F_turb_start),
        ("F_TURB_END_REF", F_TURB_END_REF, m.F_turb_end),
        ("TURB_EXP_REF", TURB_EXP_REF, m.turb_exponent),
        ("DT_TURB_MIN_REF", DT_TURB_MIN_REF, m.dT_turb_min),
        ("T_SURF_REF", T_SURF_REF, m.T_surface_ref),
        ("K_TURB_FLOOR_REF", K_TURB_FLOOR_REF, m.k_turb_floor),
        ("K_TURB_CUTOFF_REF", K_TURB_CUTOFF_REF, m.k_turb_cutoff),
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
    :T_MUSH,
    :PHI_INIT_CASE,
    :T_FORM_YR_CASE,
    :N_STEPS_CASE,
    :DT_YR_CASE,
    :NX_REFERENCE,
    :NR_REFERENCE,
    :TMAX_L2_T_MIN,
    :CENTRE_TMAX_ATOL,
    :PLATEAU_DROP_LIMIT,
    :TMAX_REL_L2_LIMIT,
    :RADIUS_ATOL_CELLS,
    :ENERGY_RESIDUAL_LIMIT,
    :ENTHALPY_RESIDUAL_LIMIT,
    REFERENCE_TURB_CONSTANTS...,
)

"""
    _constant(name) -> Any

Look up constant `name` in `MagmaOceanConductivityValidation` or `PorousShellValidation`.
"""
function _constant(name::Symbol)
    return isdefined(@__MODULE__, name) ? getfield(@__MODULE__, name) : getfield(PSV, name)
end

"""
    canonical_config_text(cfg) -> String

Canonical text hashed by `canonical_config_sha256`.
"""
function canonical_config_text(cfg)::String
    io = IOBuffer()
    MSV._print_config_sections(io, cfg, "")
    for name in (HASHED_CASE_CONSTANTS..., PSV.REFERENCE_LAW_CONSTANTS...)
        println(io, "case.", name, "=", PSV._canonical_field_repr(_constant(name)))
    end
    return String(take!(io))
end

"""
    canonical_config_sha256(cfg) -> String

SHA-256 digest of `canonical_config_text(cfg)`.
"""
canonical_config_sha256(cfg)::String = bytes2hex(sha256(canonical_config_text(cfg)))

"""
    expected_config_summary(cfg) -> Dict{String,Any}

Expected `[config]` section of the record for `cfg`.
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
        "soft_turbulence" => cfg.melting.soft_turbulence,
        "turb_exponent" => cfg.melting.turb_exponent,
        "eta_fluid_silicate_Pa_s" => cfg.melting.eta_fluid_silicate,
        "F_turb_start" => cfg.melting.F_turb_start,
        "F_turb_end" => cfg.melting.F_turb_end,
        "etamin_Pa_s" => cfg.solver.etamin,
        "rhocp_rock_J_m3K" => cfg.materials.rhocpsolidm[1],
        "k_rock_W_mK" => cfg.materials.ksolidm[1],
        "T_space_K" => cfg.disk.t_ambient,
        "surface_temperature" => String(cfg.thermodynamics.surface_temperature),
        "tmax_l2_T_min_K" => TMAX_L2_T_MIN,
        "radius_atol_cells" => RADIUS_ATOL_CELLS,
        "seed" => cfg.solver.seed,
        "spherical_metric" => cfg.geometry.spherical_metric,
    )
end

"""
    _check_recorded_limit!(errors, mid, m) -> Vector{String}

Append an error message to `errors` when the recorded tolerance of `mid` in `m` differs from `metric_limit(mid)`.
"""
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

Check that each metric in `PHYSICS_REGRESSION_BANDS` lies within its band around the
recorded value in `record_path`.

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

Compare a fresh `evaluate_magma_ocean_benchmark` result `bench` against `record_path`.

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
    MSV._check_config_summary!(errors, raw["config"], expected_config_summary(cfg))
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

# Suite record at Nx = 129, checked statically in CI.
const SUITE_REFERENCE_VALUES = (
    ("centre_tmax_nx129", 1580.6387196054766),
    ("mush_radius_nx129", 15640.0),
    ("melt_radius_nx129", 16440.0),
)
const SUITE_RESULT_HEADER = (
    "tag",
    "Nx",
    "soft_turbulence",
    "turb_exponent",
    "dx",
    "T_centre",
    "T_centre_ref",
    "plateau_drop",
    "plateau_drop_ref",
    "r_mush",
    "r_mush_ref",
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

Metric ids of the `Nx = 129` suite record.
"""
suite_metric_ids() = case_metric_ids(NX_SUITE)

"""
    suite_reference_value(id) -> Float64

Pinned 1-D reference of suite metric `id` (`0.0` for `le_atol` metrics).
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

Canonical text hashed by `suite_config_sha256`.
"""
function suite_config_text()::String
    io = IOBuffer()
    MSV._print_config_sections(io, make_case_config(; Nx=NX_SUITE), "suite.")
    for name in (SUITE_HASHED_CONSTANTS..., PSV.REFERENCE_LAW_CONSTANTS...)
        println(io, "case.", name, "=", PSV._canonical_field_repr(_constant(name)))
    end
    return String(take!(io))
end

"""
    suite_config_sha256() -> String

SHA-256 digest of `suite_config_text()`.
"""
suite_config_sha256()::String = bytes2hex(sha256(suite_config_text()))

"""
    expected_suite_config_summary() -> Dict{String,Any}

Expected `[config]` section of the `Nx = 129` suite record.
"""
expected_suite_config_summary() = expected_config_summary(make_case_config(; Nx=NX_SUITE))

"""
    check_suite_record_static(path; check_doc_sync=true) -> Tuple{Bool,Vector{String}}

Static validation of the `Nx = 129` suite record without rerunning the simulation.

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
    MSV._check_config_summary!(errors, raw["config"], expected_suite_config_summary())
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
        value = Float64(m["value"])
        ref = Float64(m["reference"])
        pinned = suite_reference_value(mid)
        isapprox(ref, pinned; rtol=1.0e-12, atol=0.0) || push!(
            errors,
            "Recorded reference of '$mid' ($ref) differs from the pinned 1-D value ($pinned)",
        )
        metric_passes(mid, value, pinned) ||
            push!(errors, "Metric '$mid' = $value fails its limit around $pinned")
    end
    return (isempty(errors), errors)
end

end # module MagmaOceanConductivityValidation

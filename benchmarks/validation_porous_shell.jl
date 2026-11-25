module PorousShellValidation

# One-dimensional spherical reference for the porous shell of a planetesimal with the
# dry porosity laws of Henke et al. (2012): conduction with the porous conductivity of
# Lichtenberg et al. (2016a, Eq. 12), 26Al heating, cold pressing and hot sintering.
# The laws are coded here from the raw constants; nothing calls the Erebus porosity code.

using Erebus
using StaticArrays
using Logging
using SHA
using TOML

include(joinpath(@__DIR__, "..", "tools", "validation_records.jl"))

const YEARLENGTH_S = 31_557_600.0
const T_HALF_AL_S = 717_000.0 * 31_540_000.0
const TAU_AL_S = T_HALF_AL_S / log(2.0)
const F_AL = 1.9e23
const E_AL = 5.0470e-13
const RATIO_AL_L16 = 5.85e-5
const G_GRAV = 6.67430e-11

const RHO_SOL = 3300.0
const RHO_MET = 5450.0
const X_FE = 0.182
const RHO_BULK = 1.0 / (X_FE / RHO_MET + (1.0 - X_FE) / RHO_SOL)
const RHOCP_SOL = 3.3e6
const K_SOL = 3.0
const K_FLOOR = 1.0e-3
const T_SPACE = 290.0

const PHI_CAP = 0.85
const T_SINTER = 700.0
const GRAIN_R0 = 1.0e-6
const CREEP_A = 4.0e-5
const E_ACT_CAL = 85000.0
const R_GAS_CAL = 1.9872
const CREEP_N = 1.5
const DPHI_SHELL = 1.0e-3
const DPHI_SINTERED = 1.0e-12

"""
Cold-pressed porosity of Henke et al. (2012) Eq. 15 at pressure `P_bar` [bar].
"""
phi_cold(P_bar::Float64) = 0.42 + 0.46 / ((P_bar / 0.13)^1.72 + 1.0)

"""
Porous conductivity [W/(m K)] of L16 Eq. 12: 3 branches in `phi`, floor `K_FLOOR`.
"""
function k_porous_ref(k::Float64, phi::Float64)
    k1 = k * exp(-phi / 0.08)
    k2 = k * exp(-1.2 - phi / 0.167)
    kp = phi < 0.2 ? k1 : (phi <= 0.4 ? (k1^4 + k2^4)^0.25 : k2)
    return max(kp, K_FLOOR)
end

"""
Grain radius `R` and contact radius `a` [m] for sintering from `phi0` to `phi`: the root
in [R0, 2 R0 X^(1/3)] of R^3 - 3 R0 X^(1/3) R^2 + R0^3 (X + 1) = 0, X = (1 - phi0) /
(1 - phi), with X at least 1/3 + 1e-12 (the contact-overlap limit). In y = R / R0 the
cubic falls and is convex on the bracket (f(1) >= 0, f(2c) < 0), so Newton from y = 1
rises monotonically to the root.
"""
function grain_geometry_ref(phi::Float64, phi0::Float64)
    X = max((1.0 - phi0) / (1.0 - phi), 1.0 / 3.0 + 1.0e-12)
    c = cbrt(X)
    y = 1.0
    for _ in 1:100
        f = y^3 - 3c * y^2 + c^3 + 1.0
        f <= 0.0 && break
        step = f / (3y * (2c - y))
        y += step
        step <= 1.0e-16 * y && break
    end
    R = GRAIN_R0 * min(y, 2c)
    a = sqrt(max(R^2 - GRAIN_R0^2 * c^2, 0.0))
    return R, a
end

"""
Porosity decrease rate [1/s] of Henke et al. (2012) Eqs. 39 to 44 (P in bar, R in cm,
R_gas in cal/(K mol)); 0 below `T_SINTER`, for P = 0 or phi = 0; Inf at phi = phi0.
"""
function sinter_rate_ref(phi::Float64, phi0::Float64, P_bar::Float64, T::Float64)
    (T < T_SINTER || P_bar <= 0.0 || phi <= 0.0) && return 0.0
    R, a = grain_geometry_ref(phi, phi0)
    a > 0.0 || return Inf
    sigma = 2.0 * sqrt(3.0) * (R^2 - a^2) * P_bar / (pi * a^2)
    F = phi <= 0.1 ? 10.0 * phi : 1.0
    return F * CREEP_A * (1.0 - phi) * sigma^CREEP_N / (100.0 * R)^3 *
           exp(-E_ACT_CAL / (R_GAS_CAL * T))
end

# Gauss-Legendre nodes and weights on [-1, 1] (Newton on the Legendre polynomial).
function _gauss_legendre(n::Int)
    x, w = zeros(n), zeros(n)
    for i in 1:n
        z = cos(pi * (i - 0.25) / (n + 0.5))
        dp = 0.0
        for _ in 1:100
            p1, p2 = 1.0, 0.0
            for j in 1:n
                p1, p2 = ((2j - 1) * z * p1 - (j - 1) * p2) / j, p1
            end
            dp = n * (z * p1 - p2) / (z^2 - 1)
            dz = p1 / dp
            z -= dz
            abs(dz) < 1.0e-15 && break
        end
        x[i], w[i] = z, 2 / ((1 - z^2) * dp^2)
    end
    return x, w
end

const GL_X, GL_W = _gauss_legendre(12)
const GL_PANELS = 6

"""
Time [s] to sinter from `phi_hi` down to `phi_lo` at fixed `P_bar` and `T`: the integral
of 1 / r over [phi_lo, phi_hi] in u = ln(phi), by composite Gauss-Legendre quadrature
(1 / r = 0 where the rate is infinite, at the onset).
"""
function sinter_time_ref(phi_lo, phi_hi, phi0, P_bar, T)
    ua, ub = log(phi_lo), log(phi_hi)
    hp = (ub - ua) / GL_PANELS
    tsum = 0.0
    for p in 1:GL_PANELS
        c = ua + (p - 0.5) * hp
        for (x, w) in zip(GL_X, GL_W)
            phi = exp(c + 0.5hp * x)
            r = sinter_rate_ref(phi, phi0, P_bar, T)
            tsum += 0.5hp * w * (isinf(r) ? 0.0 : phi / r)
        end
    end
    return tsum
end

"""
Sintered porosity after `dt` [s] at fixed `P_bar` and `T`, from the exact solution of
d phi / dt = -r(phi): the phi with sinter_time_ref(phi, phi_n) = dt, by bisection in
ln(phi) on [1e-7, phi_n]. Porosity below 1e-7 does not change.
"""
function sinter_step_ref(
    phi::Float64, phi0::Float64, P_bar::Float64, T::Float64, dt::Float64
)
    (phi < 1.0e-7 || sinter_rate_ref(phi, phi0, P_bar, T) * dt <= 0.0) && return phi
    sinter_time_ref(1.0e-7, phi, phi0, P_bar, T) <= dt && return 1.0e-7
    lo, hi = log(1.0e-7), log(phi)
    for _ in 1:60
        mid = 0.5 * (lo + hi)
        sinter_time_ref(exp(mid), phi, phi0, P_bar, T) > dt ? (lo = mid) : (hi = mid)
    end
    return exp(0.5 * (lo + hi))
end

"""
Hydrostatic pressure [Pa] at the cell centres of a sphere of radius `R` with density
`RHO_SOL (1 - phi)` per cell and gravity from the enclosed mass, P(R) = 0.
"""
function hydrostatic_pressure(phi::Vector{Float64}, R::Float64)
    Nr = length(phi)
    dr = R / Nr
    rf = [(i - 1) * dr for i in 1:(Nr + 1)]
    rc = [(i - 0.5) * dr for i in 1:Nr]
    m_face = zeros(Nr + 1)
    for i in 1:Nr
        m_face[i + 1] =
            m_face[i] + RHO_SOL * (1.0 - phi[i]) * (4pi / 3) * (rf[i + 1]^3 - rf[i]^3)
    end
    # Gravity and density at the centre of each cell; the mass inside the centre radius
    g_c = [
        G_GRAV * (m_face[i] + RHO_SOL * (1.0 - phi[i]) * (4pi / 3) * (rc[i]^3 - rf[i]^3)) /
        rc[i]^2 for i in 1:Nr
    ]
    P = zeros(Nr)
    P[Nr] = RHO_SOL * (1.0 - phi[Nr]) * g_c[Nr] * 0.5 * dr
    for i in (Nr - 1):-1:1
        P[i] =
            P[i + 1] +
            0.5 *
            (
                RHO_SOL * (1.0 - phi[i]) * g_c[i] +
                RHO_SOL * (1.0 - phi[i + 1]) * g_c[i + 1]
            ) *
            dr
    end
    return P
end

"""
Step-averaged 26Al heating [W/kg bulk] over [t, t + dt] [s] from the CAI ratio `ratio_al`.
"""
function q_al_step(t::Float64, dt::Float64, ratio_al::Float64)
    q0 = F_AL * ratio_al * E_AL / TAU_AL_S * exp(-t / TAU_AL_S)
    x = dt / TAU_AL_S
    return x > 0.0 ? q0 * (-expm1(-x)) / x : q0
end

"""
    solve_porous_sphere(; R, phi_init, t_form_yr, t_end_yr, dt_yr, Nr=2000,
                        ratio_al=RATIO_AL_L16, T0=T_SPACE, stop_cold=true) -> NamedTuple

One-dimensional spherical finite-volume run: backward Euler conduction with k, rho cp and
heat from the porosity at the start of the step, Dirichlet `T_SPACE` at `R`, then the dry
porosity step at the new temperature and the pressure of the old porosity (the order of
an Erebus step). At the start the porosity is cold pressed on the hydrostatic profile.
The run stops at `t_end_yr` after formation, or with `stop_cold` once the peak has passed
and every cell is below `T_SINTER` (no further sintering).

Returns `(; r, phi, phi0, dphi_sinter, T, T_max, d_shell, d_phi04, t_peak_yr, T_peak,
t_stop_yr, energy_residual, n_steps)`: profiles at the cell centres, the depth [m] of the
unsintered shell (`dphi_sinter < 1e-3`) and of the region with `phi > 0.4` (outer face of
the outermost sintered cell), the time [yr after formation] and value of the peak
temperature, and the largest relative energy residual of a step.
"""
function solve_porous_sphere(;
    R::Float64,
    phi_init::Float64,
    t_form_yr::Float64,
    t_end_yr::Float64,
    dt_yr::Float64,
    Nr::Int=2000,
    ratio_al::Float64=RATIO_AL_L16,
    T0::Float64=T_SPACE,
    stop_cold::Bool=true,
)
    dr = R / Nr
    rf = [(i - 1) * dr for i in 1:(Nr + 1)]
    rc = [(i - 0.5) * dr for i in 1:Nr]
    A = 4pi .* rf .^ 2
    V = [(4pi / 3) * (rf[i + 1]^3 - rf[i]^3) for i in 1:Nr]
    dt = dt_yr * YEARLENGTH_S
    # Cold pressing on the hydrostatic profile of the cold-pressed porosity itself
    phi = fill(min(phi_init, PHI_CAP), Nr)
    for _ in 1:100
        P = hydrostatic_pressure(phi, R)
        new = [min(phi_init, phi_cold(p * 1.0e-5), PHI_CAP) for p in P]
        maximum(abs.(new .- phi)) < 1.0e-14 && (phi=new; break)
        phi = new
    end
    phi0 = copy(phi)
    dphi = zeros(Nr)
    T = fill(T0, Nr)
    T_max = copy(T)
    t = t_form_yr * YEARLENGTH_S
    t_peak, T_peak, eres, nstep = 0.0, T0, 0.0, 0
    nmax = ceil(Int, t_end_yr / dt_yr)
    for n in 1:nmax
        P = hydrostatic_pressure(phi, R)
        k = [k_porous_ref(K_SOL, p) for p in phi]
        C = RHOCP_SOL .* (1.0 .- phi)
        H = RHO_BULK .* (1.0 .- phi) .* q_al_step(t, dt, ratio_al)
        # Face conductances: harmonic mean inside, half cell to the Dirichlet surface
        cond = zeros(Nr + 1)
        for i in 2:Nr
            cond[i] = A[i] / (0.5 * dr / k[i - 1] + 0.5 * dr / k[i])
        end
        cond[Nr + 1] = A[Nr + 1] * k[Nr] / (0.5 * dr)
        d = [C[i] * V[i] / dt + cond[i] + cond[i + 1] for i in 1:Nr]
        lower = [-cond[i] for i in 2:Nr]
        rhs = [C[i] * V[i] / dt * T[i] + H[i] * V[i] for i in 1:Nr]
        rhs[Nr] += cond[Nr + 1] * T_SPACE
        Tn = _thomas(lower, d, lower, rhs)
        # Energy: stored change = heat - flux through the surface, per step
        stored = sum(C[i] * V[i] * (Tn[i] - T[i]) for i in 1:Nr)
        source = dt * (sum(H .* V) - cond[Nr + 1] * (Tn[Nr] - T_SPACE))
        eres = max(eres, abs(stored - source) / max(abs(source), abs(stored), eps()))
        for i in 1:Nr
            Pb = P[i] * 1.0e-5
            if dphi[i] <= DPHI_SINTERED
                phi0[i] = min(phi0[i], phi_cold(Pb), PHI_CAP)
                phi[i] = min(phi[i], phi0[i])
            end
            pnew = sinter_step_ref(phi[i], phi0[i], Pb, Tn[i], dt)
            dphi[i] += phi[i] - pnew
            phi[i] = pnew
        end
        T = Tn
        T_max .= max.(T_max, T)
        t += dt
        nstep = n
        Tm = maximum(T)
        Tm > T_peak && ((T_peak, t_peak) = (Tm, t))
        (stop_cold && t > t_peak + dt && Tm < T_SINTER && T_peak >= T_SINTER) && break
    end
    return (;
        r=rc,
        phi,
        phi0,
        dphi_sinter=dphi,
        T,
        T_max,
        d_shell=_shell_depth_1d(dphi .>= DPHI_SHELL, R, dr),
        d_phi04=_shell_depth_1d(phi .<= 0.4, R, dr),
        t_peak_yr=t_peak / YEARLENGTH_S - t_form_yr,
        T_peak,
        t_stop_yr=t / YEARLENGTH_S - t_form_yr,
        energy_residual=eres,
        n_steps=nstep,
    )
end

# Depth below R of the outer face of the outermost cell where `mask` holds; R if none.
function _shell_depth_1d(mask::AbstractVector{Bool}, R::Float64, dr::Float64)
    i = findlast(mask)
    return i === nothing ? R : R - i * dr
end

# Tridiagonal solve (Thomas algorithm): sub `a`, diagonal `b`, super `c`, rhs `d`.
function _thomas(a, b, c, d)
    n = length(b)
    cp = zeros(n)
    dp = zeros(n)
    cp[1] = n > 1 ? c[1] / b[1] : 0.0
    dp[1] = d[1] / b[1]
    for i in 2:n
        m = b[i] - a[i - 1] * cp[i - 1]
        cp[i] = i < n ? c[i] / m : 0.0
        dp[i] = (d[i] - a[i - 1] * dp[i - 1]) / m
    end
    x = zeros(n)
    x[n] = dp[n]
    for i in (n - 1):-1:1
        x[i] = dp[i] - cp[i] * x[i + 1]
    end
    return x
end

"""
Backward-Euler eigenmode series of a sphere of radius `R` with uniform heating, Dirichlet
`T_surf` at `R`, uniform start `T_surf`: temperature at `r` after `n` steps of `dt` [s],
diffusivity `kappa` [m^2/s], step heating `H[k]` [W/m^3] and heat capacity `rhocp`.
"""
function eigenmode_T(
    r::Float64; R, T_surf, kappa, rhocp, H::Vector{Float64}, dt, n_modes::Int=4000
)
    s = 0.0
    for m in 1:n_modes
        lam = (m * pi / R)^2
        b = 2R * (-1)^(m + 1) / (m * pi)
        a = 0.0
        for h in H
            a = (a + dt * h / rhocp * b) / (1.0 + dt * kappa * lam)
        end
        s += a * (r > 0 ? sin(m * pi * r / R) / r : m * pi / R)
    end
    return T_surf + s
end

const R_PLANET = 20_000.0
const L_BOX = 2.4 * R_PLANET
const RHOCP_AIR = 3.0e10
const PHIM0 = 1.0e-4
const SEED_DEFAULT = 42

"""
    make_porous_shell_config(; Nx, phi_init, t_form_yr, n_steps, dt_yr,
                             k_air=nothing, seed=42) -> SimulationConfig

Erebus configuration of the porous shell case: a dry planetesimal of radius `R_PLANET`
in a box of `L_BOX`, spherical metric, 26Al heating from the CAI ratio `RATIO_AL_L16`,
`porosity.model = :henke2012`, a fixed time step, no melting, reactions or pore water.
The sticky air stays at `T_SPACE` (fixed disk temperature, heat capacity `RHOCP_AIR`).
`k_air = nothing` gives the air the conductivity of the cold-pressed surface rock, so
that the face between air and rock conducts as the rock; pass 3000.0 for the value that
a configuration file gives.
"""
function make_porous_shell_config(;
    Nx::Int,
    phi_init::Float64,
    t_form_yr::Float64,
    n_steps::Int,
    dt_yr::Float64,
    k_air::Union{Nothing,Float64}=nothing,
    seed::Int=SEED_DEFAULT,
)::SimulationConfig
    Nx >= 5 || throw(DomainError(Nx, "Grid resolution Nx must be >= 5"))
    n_steps >= 1 || throw(DomainError(n_steps, "Number of steps n_steps must be >= 1"))
    (isfinite(dt_yr) && dt_yr > 0.0) ||
        throw(DomainError(dt_yr, "Timestep dt_yr must be positive and finite"))
    (isfinite(t_form_yr) && t_form_yr >= 0.0) ||
        throw(DomainError(t_form_yr, "Formation time t_form_yr must be >= 0 and finite"))
    0.0 <= phi_init <= 1.0 ||
        throw(DomainError(phi_init, "Initial porosity phi_init must be in [0, 1]"))
    k_a = if k_air === nothing
        k_porous_ref(K_SOL, min(phi_init, phi_cold(0.0), PHI_CAP))
    else
        k_air
    end
    (isfinite(k_a) && k_a > 0.0) ||
        throw(DomainError(k_a, "Air conductivity k_air must be positive and finite"))
    return SimulationConfig(;
        grid=GridConfig(; xsize=L_BOX, ysize=L_BOX, Nx=Nx, Ny=Nx),
        geometry=GeometryConfig(;
            rplanet=R_PLANET,
            rcrust=R_PLANET,
            xcenter=L_BOX / 2.0,
            ycenter=L_BOX / 2.0,
            spherical_metric=true,
            metric_regularization_cells=0.5,
            gravity_mode=:enclosed_mass,
        ),
        time=TimeConfig(;
            start_time=t_form_yr,
            endtime=t_form_yr + 3.0e7,
            dt_initial=dt_yr,
            dt_longest=dt_yr,
            DTmax=500.0,
            n_steps=n_steps,
        ),
        solver=SolverConfig(; seed=seed, hydromech_solver=:direct),
        thermodynamics=ThermalConfig(;
            hr_al=true,
            hr_fe=false,
            ratio_al=RATIO_AL_L16,
            phim0=PHIM0,
            thermal_buoyancy=false,
            surface_radiation=false,
        ),
        materials=MaterialConfig(;
            ksolidm=SVector{3,Float64}([K_SOL, K_SOL, k_a]),
            rhocpsolidm=SVector{3,Float64}([RHOCP_SOL, RHOCP_SOL, RHOCP_AIR]),
            rhocpfluidm=SVector{3,Float64}([1.0e6, 1.0e6, RHOCP_AIR]),
            alphasolidm=SVector{3,Float64}([0.0, 0.0, 0.0]),
            alphafluidm=SVector{3,Float64}([0.0, 0.0, 0.0]),
            tkm0=SVector{3,Float64}([T_SPACE, T_SPACE, T_SPACE]),
            XWsolidm_init=SVector{3,Float64}([0.0, 0.0, NaN]),
        ),
        poroelasticity=PoroelasticConfig(; phimin=PHIM0, hydrofracture=false),
        porosity=PorosityConfig(; model=:henke2012, phi_init=phi_init),
        output=OutputConfig(;
            savematstep=n_steps + 10, visstep=n_steps + 10, save_final=false
        ),
        reaction=ReactionConfig(; active=false),
        melting=MeltingConfig(; active=false, soft_turbulence=false),
        coreformation=CoreFormationConfig(;
            percolation_active=false, settling_active=false
        ),
        metal_partition=MetalPartitionConfig(; active=false),
        venting=VentingConfig(; active=false),
        volatiles=VolatilesConfig(; active=false),
        retention=RetentionConfig(; active=false),
        escape=EscapeConfig(; active=false),
        atmosphere=AtmosphereConfig(; active=false),
        hydrothermal=HydrothermalConfig(; active=false),
        accretion=AccretionConfig(; active=false),
        telescoping=TelescopingConfig(; active=false),
        refractory=RefractoryConfig(; active=false),
        volatile_mixture=VolatileMixtureConfig(; active=false),
        magma_transport=MagmaTransportConfig(; active=false),
        redox=RedoxConfig(; active=false),
        magma_degassing=MagmaOceanDegassingConfig(; active=false),
        disk=DiskConfig(;
            enabled=true, model=:fixed, t_ambient=T_SPACE, dispersal_active=false
        ),
    )
end

# Median of a non-empty vector (sorted copy).
function _median(v::AbstractVector{Float64})
    s = sort(v)
    n = length(s)
    return isodd(n) ? s[(n + 1) ÷ 2] : 0.5 * (s[n ÷ 2] + s[n ÷ 2 + 1])
end

"""
    radial_bin_medians(r, values, R, dr) -> (centres, medians)

Medians of `values` in the radial bins `[(b - 1) dr, min(b dr, R))` below `R`, for
`b = 1, ..., ceil(R / dr)`: the outermost bin is partial when `R / dr` is not a whole
number. The centre of a bin is the midpoint of its part below `R`. Values at `r >= R` and
`r < 0`, and empty bins, are left out; a value at `0 <= r < R` must be finite.

Parameters
----------
r : AbstractVector{Float64}
    Radius of each value [m].
values : AbstractVector{Float64}
    Values to bin, one per radius.
R : Float64
    Outer radius [m].
dr : Float64
    Bin width [m].

Returns
-------
Tuple{Vector{Float64},Vector{Float64}}
    Bin centres [m] and the median of `values` in each non-empty bin.

Raises
------
ArgumentError
    If `r` and `values` differ in length.
DomainError
    If `R` or `dr` is not positive and finite, or a value at `0 <= r < R` is not finite.
"""
function radial_bin_medians(
    r::AbstractVector{Float64}, values::AbstractVector{Float64}, R::Float64, dr::Float64
)
    length(r) == length(values) || throw(
        ArgumentError("r and values differ in length: $(length(r)) and $(length(values))"),
    )
    (isfinite(R) && R > 0.0) ||
        throw(DomainError(R, "the outer radius R must be positive and finite"))
    (isfinite(dr) && dr > 0.0) ||
        throw(DomainError(dr, "the bin width dr must be positive and finite"))
    nb = ceil(Int, R / dr)
    bins = [Float64[] for _ in 1:nb]
    for (ri, v) in zip(r, values)
        0.0 <= ri < R || continue
        isfinite(v) || throw(DomainError(v, "the value at r = $ri m must be finite"))
        b = min(floor(Int, ri / dr) + 1, nb)
        push!(bins[b], v)
    end
    keep = findall(!isempty, bins)
    centres = [0.5 * ((b - 1) * dr + min(b * dr, R)) for b in keep]
    return centres, [_median(bins[b]) for b in keep]
end

"""
    thermal_step_energy_residual(g, dt, dx, dy; fixed=nothing)

Return the relative residual of the energy identity of the last thermal step held in the
grid fields `g`. Over the interior P nodes the stored heat is
`sum(RHOCP (tk2 - tk1)) dx dy` and the heat from the sources is `dt sum(Q) dx dy`, with
`Q = HR + HA + HS + DHP + Q_metric + Q_lat_grid + Q_seg_grid + Q_fluid_grid`, the terms of
the right-hand side of `assemble_thermal_lse!`. The boundary nodes copy their neighbour
(zero flux), so the conduction terms cancel in the sum and the 2 heats agree. The scale is
`sum(abs(RHOCP (tk2 - tk1))) dx dy + dt sum(abs(Q)) dx dy`: an air node with a large
`RHOCP` changes the stored heat of a node by more than the whole planet, so a scale from
the sums alone would cancel.

With a `fixed` mask (`thermodynamics.surface_temperature = :fixed`) the sums run over the
free nodes only, and the heat that enters the free nodes through the faces to fixed nodes,
`dt sum(K (tk2_fixed - tk2) / h^2) dx dy` with the face conductivity `KX` or `KY` of
`g`, counts as a source. Its absolute value enters the scale.

Parameters
----------
g : Grids or NamedTuple
    Grid fields after a thermal solve: `tk1`, `tk2` [K], `RHOCP` [J/(m^3 K)], `HR`, `HA`,
    `HS`, `DHP`, `Q_metric` (or `nothing`), `Q_lat_grid`, `Q_seg_grid`, `Q_fluid_grid`
    [W/m^3]; with `fixed` also `KX`, `KY` [W/(m K)].
dt : Float64
    Time step of the thermal solve [s].
dx, dy : Float64
    Grid spacing [m].
fixed : Union{Nothing,AbstractMatrix{Bool}}
    Fixed-temperature P nodes, or `nothing` for all interior nodes.

Returns
-------
Float64
    `abs(stored - source) / scale`, dimensionless.

Raises
------
ArgumentError
    If the scale is not positive and finite (no temperature change and no heat).
"""
function thermal_step_energy_residual(
    g, dt::Float64, dx::Float64, dy::Float64; fixed=nothing
)
    Ny1, Nx1 = size(g.tk1)
    stored = source = stored_abs = source_abs = 0.0
    isfixed(i, j) = fixed !== nothing && fixed[i, j]
    for j in 2:(Nx1 - 1), i in 2:(Ny1 - 1)
        isfixed(i, j) && continue
        de = g.RHOCP[i, j] * (g.tk2[i, j] - g.tk1[i, j])
        q =
            g.HR[i, j] +
            g.HA[i, j] +
            g.HS[i, j] +
            g.DHP[i, j] +
            (g.Q_metric === nothing ? 0.0 : g.Q_metric[i, j]) +
            g.Q_lat_grid[i, j] +
            g.Q_seg_grid[i, j] +
            g.Q_fluid_grid[i, j]
        if fixed !== nothing
            T = g.tk2[i, j]
            isfixed(i, j - 1) && (q += g.KX[i, j - 1] * (g.tk2[i, j - 1] - T) / dx^2)
            isfixed(i, j + 1) && (q += g.KX[i, j] * (g.tk2[i, j + 1] - T) / dx^2)
            isfixed(i - 1, j) && (q += g.KY[i - 1, j] * (g.tk2[i - 1, j] - T) / dy^2)
            isfixed(i + 1, j) && (q += g.KY[i, j] * (g.tk2[i + 1, j] - T) / dy^2)
        end
        stored += de
        source += q
        stored_abs += abs(de)
        source_abs += abs(q)
    end
    scale = (stored_abs + dt * source_abs) * dx * dy
    (isfinite(scale) && scale > 0.0) ||
        throw(ArgumentError("the energy scale must be positive and finite, got $scale"))
    return abs(stored - dt * source) * dx * dy / scale
end

"""
    check_run_clock(st, cfg)

Check that a finished run took `cfg.time.n_steps` steps of the fixed step
`cfg.time.dt_initial`: the time step of the last step equals the case step, and the final
model time `st.timesum` equals `start_time + n_steps dt`, both in seconds and within a
relative `1e-12`.

Parameters
----------
st : SimulationState or NamedTuple
    Final state of `simulation_loop`; `dt` [s] and `timesum` [s] are read.
cfg : SimulationConfig
    Configuration of the run, from `make_porous_shell_config`.

Returns
-------
Nothing

Raises
------
ErrorException
    If the last time step or the final time differs from the case clock.
"""
function check_run_clock(st, cfg::SimulationConfig)
    yl = cfg.time.yearlength
    dt_s = cfg.time.dt_initial * yl
    isapprox(st.dt, dt_s; rtol=1.0e-12) ||
        error("the last time step is $(st.dt) s, not the case step $dt_s s")
    t_end_s = cfg.time.start_time * yl + cfg.time.n_steps * dt_s
    isapprox(st.timesum, t_end_s; rtol=1.0e-12) ||
        error("the final time is $(st.timesum) s, not start_time + n_steps dt = $t_end_s s")
    return nothing
end

"""
    run_porous_shell_case(; Nx, phi_init, t_form_yr, n_steps, dt_yr, k_air=nothing,
                          seed=42, Nr=2000, reference=nothing) -> NamedTuple

Run the Erebus case and the 1-D reference with the same time step and number of steps,
and compare them. The shell depth of Erebus is `Erebus.porous_shell_depth` with radial
bins of half a cell. The profiles of the maximum temperature and of the porosity are
medians of the rock markers in radial bins of half a cell, compared with the 1-D cell at
the bin centre. The energy residual is `thermal_step_energy_residual` of the last step.

`reference` is a `solve_porous_sphere` result of the same case (`stop_cold = false`),
or `nothing` to compute it here with `Nr` cells.

Returns `(; dx, d_shell, d_shell_q25, d_shell_q75, d_shell_ref, d_shell_error_cells,
tmax_rel_l2, phi_rms, T_peak, T_peak_ref, energy_lse_residual, r, tmax, tmax_ref, phi,
phi_ref, reference)`: depths [m], the relative L2 norm of `T_max - T_SPACE`, the rms
porosity difference [-], peak temperatures [K], the relative energy residual of the last
thermal step [-], the binned profiles and the full 1-D result.

Raises `DomainError` for the arguments that `make_porous_shell_config` rejects, and an
`ErrorException` from `check_run_clock` when the run did not take `n_steps` steps of
`dt_yr`.
"""
function run_porous_shell_case(;
    Nx::Int,
    phi_init::Float64,
    t_form_yr::Float64,
    n_steps::Int,
    dt_yr::Float64,
    k_air::Union{Nothing,Float64}=nothing,
    seed::Int=SEED_DEFAULT,
    Nr::Int=2000,
    reference=nothing,
)
    cfg = make_porous_shell_config(; Nx, phi_init, t_form_yr, n_steps, dt_yr, k_air, seed)
    st = mktempdir() do tmpdir
        redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                with_logger(NullLogger()) do
                    return simulation_loop(cfg; output_path=tmpdir)
                end
            end
        end
    end
    check_run_clock(st, cfg)
    ref = if reference === nothing
        solve_porous_sphere(;
            R=R_PLANET,
            phi_init,
            t_form_yr,
            t_end_yr=n_steps * dt_yr,
            dt_yr,
            Nr,
            stop_cold=false,
        )
    else
        reference
    end
    mk = st.markers
    pg = mk.groups.porosity
    xm, ym, tm = mk.core.xm, mk.core.ym, mk.core.tm
    dx = L_BOX / (Nx - 1)
    xc = yc = L_BOX / 2.0
    sh = Erebus.porous_shell_depth(
        xm, ym, tm, pg; xcenter=xc, ycenter=yc, rplanet=R_PLANET, dr=dx / 2
    )
    rock = findall(t -> t == 1 || t == 2, tm)
    r = hypot.(xm[rock] .- xc, ym[rock] .- yc)
    rc, tmax = radial_bin_medians(r, pg.tkm_max[rock], R_PLANET, dx / 2)
    _, phi = radial_bin_medians(r, pg.phi_d[rock], R_PLANET, dx / 2)
    cell = [clamp(floor(Int, x / (R_PLANET / Nr)) + 1, 1, Nr) for x in rc]
    tmax_ref = ref.T_max[cell]
    phi_ref = ref.phi[cell]
    return (;
        dx,
        d_shell=sh.median,
        d_shell_q25=sh.q25,
        d_shell_q75=sh.q75,
        d_shell_ref=ref.d_shell,
        d_shell_error_cells=(sh.median - ref.d_shell) / dx,
        tmax_rel_l2=sqrt(sum(abs2, tmax .- tmax_ref) / sum(abs2, tmax_ref .- T_SPACE)),
        phi_rms=sqrt(sum(abs2, phi .- phi_ref) / length(phi)),
        T_peak=maximum(pg.tkm_max[rock]),
        T_peak_ref=ref.T_peak,
        energy_lse_residual=thermal_step_energy_residual(st.grids, st.dt, dx, dx),
        r=rc,
        tmax,
        tmax_ref,
        phi,
        phi_ref,
        reference=ref,
    )
end

const CASE_ID = "porous_shell_t2.5_phi0.25"
const PHI_INIT_CASE = 0.25
const T_FORM_YR_CASE = 2.5e6
const N_STEPS_CASE = 80
const DT_YR_CASE = 5.0e4
const NX_LIST = (65,)
const NX_REFERENCE = 65
const NR_REFERENCE = 2000
const SHELL_DEPTH_REL_BOUND = 0.15
const TMAX_REL_L2_LIMIT = 0.10
const ENERGY_RESIDUAL_LIMIT = 1.0e-9

# Constants of the laws and of the 1-D reference, hashed in both records.
const REFERENCE_LAW_CONSTANTS = (
    :YEARLENGTH_S,
    :T_HALF_AL_S,
    :F_AL,
    :E_AL,
    :RATIO_AL_L16,
    :G_GRAV,
    :RHO_SOL,
    :RHO_MET,
    :X_FE,
    :RHOCP_SOL,
    :K_SOL,
    :K_FLOOR,
    :T_SPACE,
    :PHI_CAP,
    :T_SINTER,
    :GRAIN_R0,
    :CREEP_A,
    :E_ACT_CAL,
    :R_GAS_CAL,
    :CREEP_N,
    :DPHI_SHELL,
    :DPHI_SINTERED,
)

const HASHED_CASE_CONSTANTS = (
    :R_PLANET,
    :L_BOX,
    :RHOCP_AIR,
    :PHIM0,
    :NX_LIST,
    :NX_REFERENCE,
    :NR_REFERENCE,
    :PHI_INIT_CASE,
    :T_FORM_YR_CASE,
    :N_STEPS_CASE,
    :DT_YR_CASE,
    :SHELL_DEPTH_REL_BOUND,
    :TMAX_REL_L2_LIMIT,
    :ENERGY_RESIDUAL_LIMIT,
    REFERENCE_LAW_CONSTANTS...,
)

const BENCHMARK_CONFIG_SECTIONS = (
    :grid,
    :geometry,
    :time,
    :solver,
    :thermodynamics,
    :materials,
    :poroelasticity,
    :porosity,
    :reaction,
    :melting,
    :coreformation,
    :metal_partition,
    :venting,
    :magma_transport,
    :disk,
)

const EXPECTED_METRIC_IDS = (
    "shell_depth_nx65", "tmax_rel_l2_nx65", "energy_lse_residual_nx65"
)

# A shell depth changes only by whole radial bins of dx / 2: its band is one bin. The
# temperature band is 3 times its spread over the marker seeds 42, 1 and 2. The energy
# residual (round-off) has no band; its pass limit applies.
const PHYSICS_REGRESSION_BANDS = Dict{
    String,NamedTuple{(:atol, :rtol),Tuple{Float64,Float64}}
}(
    "shell_depth_nx65" => (; atol=375.0, rtol=0.0),
    "tmax_rel_l2_nx65" => (; atol=1.29e-3, rtol=1.0e-12),
)

const REPLAY_RTOL = 1.0e-6
# The shell depth moves in whole radial bins (375 m at Nx 65): one marker across a bin edge
# on another platform is round-off, so its replay allows one bin.
const REPLAY_ATOL_BOUNDS = Dict{String,Float64}(
    "shell_depth_nx65" => 375.0, "energy_lse_residual_nx65" => 1.0e-9
)

"""
    make_case_config(; Nx=NX_REFERENCE, seed=SEED_DEFAULT)

Return the `SimulationConfig` of the recorded case on an `Nx` x `Nx` grid.

Parameters
----------
Nx : Int, optional
    Basic grid resolution in x and y, at least 5 (default `NX_REFERENCE`).
seed : Int, optional
    Seed of the marker positions (default `SEED_DEFAULT`).

Returns
-------
SimulationConfig
    `make_porous_shell_config` with `PHI_INIT_CASE`, `T_FORM_YR_CASE`, `N_STEPS_CASE`,
    `DT_YR_CASE` and matched air (`k_air = nothing`).

Raises
------
DomainError
    If `Nx < 5`.
"""
function make_case_config(; Nx::Int=NX_REFERENCE, seed::Int=SEED_DEFAULT)
    return make_porous_shell_config(;
        Nx,
        phi_init=PHI_INIT_CASE,
        t_form_yr=T_FORM_YR_CASE,
        n_steps=N_STEPS_CASE,
        dt_yr=DT_YR_CASE,
        k_air=nothing,
        seed,
    )
end

"""
    reference_solution()

Return the 1-D reference of the recorded case: `solve_porous_sphere` with `NR_REFERENCE`
cells, the case time step and number of steps, and no stop once the rock cools.

Returns
-------
NamedTuple
    The `solve_porous_sphere` result.
"""
function reference_solution()
    return solve_porous_sphere(;
        R=R_PLANET,
        phi_init=PHI_INIT_CASE,
        t_form_yr=T_FORM_YR_CASE,
        t_end_yr=N_STEPS_CASE * DT_YR_CASE,
        dt_yr=DT_YR_CASE,
        Nr=NR_REFERENCE,
        stop_cold=false,
    )
end

"""
    shell_depth_atol(nx, d_shell_ref)

Return the bound of the shell depth error on an `nx` x `nx` grid: one cell `dx = L_BOX /
(nx - 1)` or `SHELL_DEPTH_REL_BOUND` of the reference depth, whichever is larger.

Parameters
----------
nx : Int
    Basic grid resolution, at least 2.
d_shell_ref : Float64
    Shell depth of the 1-D reference [m], at least 0.

Returns
-------
Float64
    Absolute tolerance [m].

Raises
------
DomainError
    If `nx < 2` or `d_shell_ref` is negative or not finite.
"""
function shell_depth_atol(nx::Int, d_shell_ref::Float64)
    nx >= 2 || throw(DomainError(nx, "the grid resolution must be at least 2"))
    (isfinite(d_shell_ref) && d_shell_ref >= 0.0) ||
        throw(DomainError(d_shell_ref, "the reference shell depth must be finite and >= 0"))
    return max(L_BOX / (nx - 1), SHELL_DEPTH_REL_BOUND * d_shell_ref)
end

"""
    metrics_from_runs(runs)

Return the validation metrics that the `run_porous_shell_case` results in `runs` determine.

The shell depth exists for each grid of `NX_LIST` present in `runs`. The maximum
temperature error and the energy residual need the `NX_REFERENCE` grid.

Parameters
----------
runs : AbstractDict
    Map from `Nx` to the `run_porous_shell_case` result.

Returns
-------
Dict{String,Float64}
    Metric values by id (see `EXPECTED_METRIC_IDS`); an id whose runs are missing is
    absent.
"""
function metrics_from_runs(runs::AbstractDict)::Dict{String,Float64}
    metrics = Dict{String,Float64}()
    for n in NX_LIST
        haskey(runs, n) && (metrics["shell_depth_nx$(n)"] = runs[n].d_shell)
    end
    if haskey(runs, NX_REFERENCE)
        m = runs[NX_REFERENCE]
        metrics["tmax_rel_l2_nx$(NX_REFERENCE)"] = m.tmax_rel_l2
        metrics["energy_lse_residual_nx$(NX_REFERENCE)"] = m.energy_lse_residual
    end
    return metrics
end

"""
    evaluate_porous_shell_benchmark(; seed=SEED_DEFAULT, nx_list=NX_LIST)

Run the case at every grid in `nx_list` against one 1-D reference and return the runs and
the validation metrics they determine.

Parameters
----------
seed : Int, optional
    Seed of the marker positions (default `SEED_DEFAULT`).
nx_list : Tuple{Vararg{Int}}, optional
    Grids to run, a subset of `NX_LIST` (default `NX_LIST`). With the default the result
    holds all 3 metrics.

Returns
-------
NamedTuple
    `runs`: `Dict` from `Nx` to the `run_porous_shell_case` result;
    `metrics`: `Dict{String,Float64}` from metric id (see `EXPECTED_METRIC_IDS`) to value,
    from `metrics_from_runs`.

Raises
------
ArgumentError
    If a grid in `nx_list` is not in `NX_LIST`.
"""
function evaluate_porous_shell_benchmark(; seed::Int=SEED_DEFAULT, nx_list=NX_LIST)
    all(n -> n isa Integer && n in NX_LIST, nx_list) ||
        throw(ArgumentError("nx_list must be a subset of $NX_LIST, got $nx_list"))
    ref = reference_solution()
    runs = Dict{Int,Any}()
    for n in nx_list
        runs[n] = run_porous_shell_case(;
            Nx=n,
            phi_init=PHI_INIT_CASE,
            t_form_yr=T_FORM_YR_CASE,
            n_steps=N_STEPS_CASE,
            dt_yr=DT_YR_CASE,
            seed,
            Nr=NR_REFERENCE,
            reference=ref,
        )
        GC.gc()
    end
    return (; runs, metrics=metrics_from_runs(runs))
end

"""
    expected_metric_references(; d_shell_ref=reference_solution().d_shell)

Return the reference value of each validation metric: the 1-D shell depth for the shell
depths and `0.0` for the maximum temperature error and the energy residual.

Parameters
----------
d_shell_ref : Float64, optional
    Shell depth of the 1-D reference [m] (default: computed by `reference_solution`).

Returns
-------
Dict{String,Float64}
    Mapping from metric id to its reference value.
"""
function expected_metric_references(;
    d_shell_ref::Float64=reference_solution().d_shell
)::Dict{String,Float64}
    refs = Dict{String,Float64}()
    for mid in EXPECTED_METRIC_IDS
        refs[mid] = startswith(mid, "shell_depth") ? d_shell_ref : 0.0
    end
    return refs
end

"""
    expected_metric_limits(; d_shell_ref=reference_solution().d_shell)

Return the pass limit of each validation metric, fixed before any run: the shell depth
within `shell_depth_atol(nx, d_shell_ref)` of the 1-D depth (`isapprox`), the maximum
temperature error at most `TMAX_REL_L2_LIMIT` and the energy residual at most
`ENERGY_RESIDUAL_LIMIT` (`le_atol`). Every limit has `rtol = 0`.

Parameters
----------
d_shell_ref : Float64, optional
    Shell depth of the 1-D reference [m] (default: computed by `reference_solution`).

Returns
-------
Dict{String,NamedTuple}
    Mapping from metric id to `(; atol, rtol, comparison)`.

Raises
------
DomainError
    If `d_shell_ref` is negative or not finite (from `shell_depth_atol`).
"""
function expected_metric_limits(; d_shell_ref::Float64=reference_solution().d_shell)
    limits = Dict{
        String,NamedTuple{(:atol, :rtol, :comparison),Tuple{Float64,Float64,String}}
    }()
    for n in NX_LIST
        limits["shell_depth_nx$(n)"] = (;
            atol=shell_depth_atol(n, d_shell_ref), rtol=0.0, comparison="isapprox"
        )
    end
    limits["tmax_rel_l2_nx$(NX_REFERENCE)"] = (;
        atol=TMAX_REL_L2_LIMIT, rtol=0.0, comparison="le_atol"
    )
    limits["energy_lse_residual_nx$(NX_REFERENCE)"] = (;
        atol=ENERGY_RESIDUAL_LIMIT, rtol=0.0, comparison="le_atol"
    )
    return limits
end

"""
    metric_units(id)

Return the units of a metric of the porous shell records, as the exporters write them:
`"m"` for a shell depth and `"dimensionless"` for a maximum temperature error or an energy
residual.

Parameters
----------
id : AbstractString
    Metric id, `shell_depth_<run>`, `tmax_rel_l2_<run>` or `energy_lse_residual_<run>`.

Returns
-------
String
    `"m"` or `"dimensionless"`.

Raises
------
ArgumentError
    If `id` has none of the 3 prefixes.
"""
function metric_units(id::AbstractString)::String
    startswith(id, "shell_depth_") && return "m"
    (startswith(id, "tmax_rel_l2_") || startswith(id, "energy_lse_residual_")) &&
        return "dimensionless"
    return throw(ArgumentError("unknown metric id '$id'"))
end

"""
    bench_shell_reference(bench)

Return the 1-D shell depth of the runs in `bench`, the reference of the shell depth
metrics, or compute it with `reference_solution` when `bench` holds no run.

Parameters
----------
bench : NamedTuple
    Result of `evaluate_porous_shell_benchmark`, or a tuple with `metrics` only.

Returns
-------
Float64
    Shell depth of the 1-D reference [m].

Raises
------
ArgumentError
    If the runs of `bench` disagree on the 1-D shell depth.
"""
function bench_shell_reference(bench)::Float64
    runs = hasproperty(bench, :runs) ? bench.runs : Dict{Int,Any}()
    depths = [r.d_shell_ref for r in values(runs) if hasproperty(r, :d_shell_ref)]
    isempty(depths) && return reference_solution().d_shell
    all(d -> isapprox(d, first(depths); rtol=1.0e-14, atol=0.0), depths) ||
        throw(ArgumentError("the runs disagree on the 1-D shell depth: $depths"))
    return first(depths)
end

"""
    check_reference_literals_match_config(cfg=make_case_config(); x_fe_ref=Erebus.X_FE_REF_CHONDRITE, g_grav=Erebus.G_GRAV)

Compare the constants of the 1-D reference with the configuration and the Erebus
constants they stand for: `YEARLENGTH_S` (`cfg.time.yearlength`), `T_HALF_AL_S`, `F_AL`,
`E_AL`, `RATIO_AL_L16` (`cfg.thermodynamics`), `RHO_SOL`, `RHOCP_SOL`, `K_SOL`
(`cfg.materials`, rock), `RHO_MET` (`cfg.coreformation.rho_metal`), `X_FE` (`x_fe_ref`),
`G_GRAV` (`g_grav`), `T_SPACE` (`cfg.disk.t_ambient`) and `PHI_CAP`, `T_SINTER`,
`GRAIN_R0`, `CREEP_A`, `E_ACT_CAL`, `CREEP_N` (`cfg.porosity`).

Parameters
----------
cfg : SimulationConfig, optional
    Configuration to compare with (default `make_case_config()`).
x_fe_ref : Float64, optional
    Iron mass fraction of the chondritic rock (default `Erebus.X_FE_REF_CHONDRITE`).
g_grav : Float64, optional
    Gravitational constant [m^3/(kg s^2)] (default `Erebus.G_GRAV`).

Returns
-------
Tuple{Bool,Vector{String}}
    Pass flag and list of mismatch descriptions.
"""
function check_reference_literals_match_config(
    cfg::SimulationConfig=make_case_config();
    x_fe_ref::Float64=Erebus.X_FE_REF_CHONDRITE,
    g_grav::Float64=Erebus.G_GRAV,
)::Tuple{Bool,Vector{String}}
    checks = (
        ("YEARLENGTH_S", YEARLENGTH_S, cfg.time.yearlength),
        ("T_HALF_AL_S", T_HALF_AL_S, cfg.thermodynamics.t_half_al),
        ("F_AL", F_AL, cfg.thermodynamics.f_al),
        ("E_AL", E_AL, cfg.thermodynamics.E_al),
        ("RATIO_AL_L16", RATIO_AL_L16, cfg.thermodynamics.ratio_al),
        ("RHO_SOL", RHO_SOL, cfg.materials.rhosolidm[1]),
        ("RHOCP_SOL", RHOCP_SOL, cfg.materials.rhocpsolidm[1]),
        ("K_SOL", K_SOL, cfg.materials.ksolidm[1]),
        ("RHO_MET", RHO_MET, cfg.coreformation.rho_metal),
        ("X_FE", X_FE, x_fe_ref),
        ("G_GRAV", G_GRAV, g_grav),
        ("T_SPACE", T_SPACE, cfg.disk.t_ambient),
        ("PHI_CAP", PHI_CAP, cfg.porosity.phi_cap),
        ("T_SINTER", T_SINTER, cfg.porosity.T_sinter),
        ("GRAIN_R0", GRAIN_R0, cfg.porosity.grain_radius),
        ("CREEP_A", CREEP_A, cfg.porosity.A),
        ("E_ACT_CAL", E_ACT_CAL, cfg.porosity.E_act),
        ("CREEP_N", CREEP_N, cfg.porosity.n),
    )
    errors = String[]
    for (name, lit_val, cfg_val) in checks
        if !isapprox(lit_val, cfg_val; rtol=1.0e-14, atol=0.0)
            push!(
                errors,
                "Reference literal $name ($lit_val) does not match config value ($cfg_val)",
            )
        end
    end
    return (isempty(errors), errors)
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
    Simulation configuration (or section container) of the case, from `make_case_config`.

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
sections of the case, the grid list `NX_LIST`, the case run and every constant of the 1-D
reference (`HASHED_CASE_CONSTANTS`).

Parameters
----------
cfg : Union{SimulationConfig,NamedTuple}
    Simulation configuration (or section container) of the case, from `make_case_config`.

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
    Configuration of the case on the `NX_REFERENCE` grid, from `make_case_config`.

Returns
-------
Dict{String,Any}
    Physical and numerical parameters of the case, keyed as in `REQUIRED_CONFIG_KEYS`.
"""
function expected_config_summary(cfg::SimulationConfig)::Dict{String,Any}
    return Dict{String,Any}(
        "Nx" => collect(NX_LIST),
        "Nr_reference" => NR_REFERENCE,
        "n_steps" => cfg.time.n_steps,
        "dt_yr" => cfg.time.dt_initial,
        "t_form_yr" => cfg.time.start_time,
        "phi_init" => cfg.porosity.phi_init,
        "porosity_model" => String(cfg.porosity.model),
        "R_planet_m" => cfg.geometry.rplanet,
        "L_box_m" => cfg.grid.xsize,
        "ratio_al" => cfg.thermodynamics.ratio_al,
        "k_rock_W_mK" => cfg.materials.ksolidm[1],
        "k_air_W_mK" => cfg.materials.ksolidm[3],
        "rhocp_rock_J_m3K" => cfg.materials.rhocpsolidm[1],
        "rhocp_air_J_m3K" => cfg.materials.rhocpsolidm[3],
        "T_space_K" => cfg.disk.t_ambient,
        "shell_dphi" => DPHI_SHELL,
        "shell_depth_rel_bound" => SHELL_DEPTH_REL_BOUND,
        "seed" => cfg.solver.seed,
        "spherical_metric" => cfg.geometry.spherical_metric,
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
    _push_extra_config_keys!(errors, rec_cfg, summary)

Append to `errors` one message per key of the recorded `[config]` section that the
expected summary does not have, in sorted order.

Parameters
----------
errors : Vector{String}
    Error messages, changed in place.
rec_cfg : AbstractDict
    Recorded `[config]` section.
summary : AbstractDict
    Expected configuration summary.

Returns
-------
Vector{String}
    `errors`.
"""
function _push_extra_config_keys!(
    errors::Vector{String}, rec_cfg::AbstractDict, summary::AbstractDict
)
    for k in sort!([String(k) for k in keys(rec_cfg) if !haskey(summary, k)])
        push!(
            errors, "Recorded config has key '$k' that the expected summary does not have"
        )
    end
    return errors
end

"""
    _push_units_error!(errors, mid, m)

Append to `errors` a message when the recorded units of the metric `mid` differ from
`metric_units(mid)`.

Parameters
----------
errors : Vector{String}
    Error messages, changed in place.
mid : String
    Metric id with a prefix that `metric_units` knows.
m : AbstractDict
    Recorded metric table with a string `units`.

Returns
-------
Vector{String}
    `errors`.

Raises
------
ArgumentError
    From `metric_units`, if `mid` has no known prefix.
"""
function _push_units_error!(errors::Vector{String}, mid::String, m::AbstractDict)
    units, expected = String(m["units"]), metric_units(mid)
    units == expected || push!(errors, "Metric '$mid' has units '$units', not '$expected'")
    return errors
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
    check_validation_reproducibility(bench, record_path; check_recorded_values=true, d_shell_ref=nothing)

Compare freshly evaluated metrics, reference values, regression bands, the reference
literals and the configuration hash and summary with a committed validation record. A
`[config]` key outside `expected_config_summary` is an error, and every metric must have
`nx = NX_REFERENCE`, `case = CASE_ID` and, for an id of `EXPECTED_METRIC_IDS`, the units
of `metric_units`. The recorded `atol`, `rtol` and `comparison` of each metric must
equal the limit of the code (`expected_metric_limits`), so that a widened tolerance in
the record fails. A recorded metric without a reference in `expected_metric_references`
or without a limit in `expected_metric_limits` is an error. The replay of recorded values (`REPLAY_RTOL`,
`REPLAY_ATOL_BOUNDS`) runs when `check_recorded_values` is `true` and the Julia
major.minor version matches
`metadata.julia_version`; the regression bands always run. A record that fails
`validate_record_dict` returns its schema errors only, and a file that is not valid TOML
the message of the parser error.

Parameters
----------
bench : NamedTuple
    Result of `evaluate_porous_shell_benchmark`.
record_path : String
    Path to the TOML validation record file.
check_recorded_values : Bool, optional
    Whether to replay the recorded values (default `true`).
d_shell_ref : Union{Float64,Nothing}, optional
    Shell depth of the 1-D reference [m], or `nothing` to take it from the runs of
    `bench` with `bench_shell_reference` (default).

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
    bench,
    record_path::String;
    check_recorded_values::Bool=true,
    d_shell_ref::Union{Float64,Nothing}=nothing,
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

    cfg = make_case_config()
    ok_lit, lit_errs = check_reference_literals_match_config(cfg)
    ok_lit || append!(errors, lit_errs)
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
    summary = expected_config_summary(cfg)
    _push_extra_config_keys!(errors, rec_cfg, summary)
    for (k, v_exp) in summary
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

    exp_refs = nothing
    exp_limits = nothing
    for m in raw_metrics
        !(m isa AbstractDict && haskey(m, "id")) && continue
        mid = String(m["id"])
        # The exporter writes every metric at NX_REFERENCE, of CASE_ID, in its units
        Int(m["nx"]) == NX_REFERENCE ||
            push!(errors, "Metric '$mid' has nx $(m["nx"]), not $NX_REFERENCE")
        String(m["case"]) == CASE_ID ||
            push!(errors, "Metric '$mid' has case '$(m["case"])', not '$CASE_ID'")
        mid in EXPECTED_METRIC_IDS && _push_units_error!(errors, mid, m)
        if !haskey(bench.metrics, mid)
            push!(errors, "Benchmark output missing metric '$mid'")
            continue
        end
        if exp_refs === nothing
            d_ref = d_shell_ref === nothing ? bench_shell_reference(bench) : d_shell_ref
            exp_refs = expected_metric_references(; d_shell_ref=d_ref)
            exp_limits = expected_metric_limits(; d_shell_ref=d_ref)
        end
        measured = bench.metrics[mid]
        recorded = Float64(m["value"])
        ref = Float64(m["reference"])
        atol = Float64(m["atol"])
        rtol = Float64(m["rtol"])
        comp = String(m["comparison"])
        # A metric that the code does not define has no reference or limit to check
        # its recorded values against, so it is an error.
        if !haskey(exp_refs, mid)
            push!(errors, "Recorded metric '$mid' has no recomputed reference in the code")
        elseif !isapprox(ref, exp_refs[mid]; rtol=1.0e-12, atol=1.0e-12)
            push!(
                errors,
                "Recorded reference for '$mid' ($ref) does not match recomputed reference ($(exp_refs[mid]))",
            )
        end
        if !haskey(exp_limits, mid)
            push!(
                errors,
                "Recorded metric '$mid' has no limit in expected_metric_limits of the code",
            )
            continue
        end
        lim = exp_limits[mid]
        same_limit =
            comp == lim.comparison &&
            isapprox(atol, lim.atol; rtol=1.0e-12, atol=0.0) &&
            isapprox(rtol, lim.rtol; rtol=1.0e-12, atol=0.0)
        same_limit || push!(
            errors,
            "Recorded tolerance of '$mid' ($comp, atol=$atol, rtol=$rtol) differs from the limit of the code ($(lim.comparison), atol=$(lim.atol), rtol=$(lim.rtol))",
        )
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

# Full runs of the suite: (tag, Nx, t_form [Myr], phi_init, steps of SUITE_DT_YR to the
# shell freeze, air conductivity [W/(m K)] or `nothing` for the matched air).
const SUITE_CASES = (
    ("s65_25_025", 65, 2.5, 0.25, 750, nothing),
    ("s129_25_025", 129, 2.5, 0.25, 750, nothing),
    ("s129_225_025", 129, 2.25, 0.25, 700, nothing),
    ("s129_175_010", 129, 1.75, 0.1, 490, nothing),
    ("s129_25_050", 129, 2.5, 0.5, 1580, nothing),
    ("s129_25_025_k", 129, 2.5, 0.25, 750, 3000.0),
)
# Shell depth [m] of the 1-D reference of each case of SUITE_CASES, by tag:
# solve_porous_sphere with NR_REFERENCE cells and the steps of the case, no stop once the
# rock cools.
const SUITE_REFERENCE_DEPTHS = (
    ("s65_25_025", 4560.0),
    ("s129_25_025", 4560.0),
    ("s129_225_025", 2340.0),
    ("s129_175_010", 2170.0),
    ("s129_25_050", 2440.0),
    ("s129_25_025_k", 4560.0),
)
const SUITE_DT_YR = 1.0e4
const SUITE_TMAX_NX_MIN = 65
const SUITE_RECORD_ID = "porous_shell_suite"
const SUITE_RESULT_HEADER = (
    "tag",
    "Nx",
    "t_form_myr",
    "phi_init",
    "n_steps",
    "k_air",
    "dx",
    "d_shell",
    "q25",
    "q75",
    "d_ref",
    "err_cells",
    "tmax_rel_l2",
    "phi_rms",
    "T_peak",
    "T_peak_ref",
    "energy_res",
    "wall_s",
)
const SUITE_HASHED_CONSTANTS = (
    :SUITE_REFERENCE_DEPTHS,
    :SUITE_DT_YR,
    :SUITE_TMAX_NX_MIN,
    :R_PLANET,
    :L_BOX,
    :RHOCP_AIR,
    :PHIM0,
    :NR_REFERENCE,
    :SHELL_DEPTH_REL_BOUND,
    :TMAX_REL_L2_LIMIT,
    REFERENCE_LAW_CONSTANTS...,
)

"""
    suite_case_config(case)

Return the `SimulationConfig` of one entry of `SUITE_CASES`.

Parameters
----------
case : Tuple
    Entry `(tag, Nx, t_form_myr, phi_init, n_steps, k_air)` of `SUITE_CASES`; `t_form_myr`
    in Myr after the first solids, `k_air` in W/(m K) or `nothing`.

Returns
-------
SimulationConfig
    `make_porous_shell_config` with the step `SUITE_DT_YR`.

Raises
------
DomainError
    For the arguments that `make_porous_shell_config` rejects.
"""
function suite_case_config(case::Tuple)
    _, nx, tf_myr, phi, n, k_air = case
    return make_porous_shell_config(;
        Nx=nx, phi_init=phi, t_form_yr=tf_myr * 1.0e6, n_steps=n, dt_yr=SUITE_DT_YR, k_air
    )
end

"""
    suite_reference_depth(tag)

Return the pinned shell depth of the 1-D reference of the suite case `tag`, from
`SUITE_REFERENCE_DEPTHS`.

Parameters
----------
tag : AbstractString
    Tag of a case of `SUITE_CASES`.

Returns
-------
Float64
    Shell depth of the 1-D reference [m].

Raises
------
ArgumentError
    If `tag` has no entry in `SUITE_REFERENCE_DEPTHS`.
"""
function suite_reference_depth(tag::AbstractString)::Float64
    for (t, d) in SUITE_REFERENCE_DEPTHS
        t == tag && return d
    end
    return throw(ArgumentError("suite case '$tag' has no pinned 1-D shell depth"))
end

"""
    suite_metric_ids()

Return the metric ids of the suite record in record order: for each case of `SUITE_CASES`
the shell depth `shell_depth_<tag>`, then `tmax_rel_l2_<tag>` when `Nx >= SUITE_TMAX_NX_MIN`.

Returns
-------
Vector{String}
    Metric ids.
"""
function suite_metric_ids()::Vector{String}
    ids = String[]
    for (tag, nx) in SUITE_CASES
        push!(ids, "shell_depth_$(tag)")
        nx >= SUITE_TMAX_NX_MIN && push!(ids, "tmax_rel_l2_$(tag)")
    end
    return ids
end

"""
    suite_metric_limit(id, nx, d_shell_ref)

Return the pass limit of a suite metric: the shell depth within
`shell_depth_atol(nx, d_shell_ref)`, that is `max(dx, 0.15 d_shell_ref)`, of the 1-D depth
(`isapprox`), and the maximum temperature error at most `TMAX_REL_L2_LIMIT` (`le_atol`).

Parameters
----------
id : AbstractString
    Metric id, `shell_depth_<tag>` or `tmax_rel_l2_<tag>`.
nx : Int
    Grid resolution of the run.
d_shell_ref : Float64
    Shell depth of the 1-D reference of the run [m].

Returns
-------
NamedTuple
    `(; atol, rtol, comparison)`.

Raises
------
ArgumentError
    If `id` is neither a shell depth nor a maximum temperature error.
DomainError
    From `shell_depth_atol`.
"""
function suite_metric_limit(id::AbstractString, nx::Int, d_shell_ref::Float64)
    startswith(id, "shell_depth_") &&
        return (; atol=shell_depth_atol(nx, d_shell_ref), rtol=0.0, comparison="isapprox")
    startswith(id, "tmax_rel_l2_") &&
        return (; atol=TMAX_REL_L2_LIMIT, rtol=0.0, comparison="le_atol")
    return throw(ArgumentError("unknown suite metric id '$id'"))
end

"""
    suite_config_text()

Return the text that `suite_config_sha256` hashes: per case of `SUITE_CASES` one line with
the case tuple and one line `suite.<tag>.section.field=value` per field of the sections in
`BENCHMARK_CONFIG_SECTIONS` of `suite_case_config`, then one line `case.NAME=value` per
constant of `SUITE_HASHED_CONSTANTS`.

Returns
-------
String
    The canonical text of the case set.
"""
function suite_config_text()::String
    io = IOBuffer()
    for case in SUITE_CASES
        tag = case[1]
        println(io, "suite.", tag, "=", _canonical_field_repr(collect(Any, case)))
        cfg = suite_case_config(case)
        for sec_sym in BENCHMARK_CONFIG_SECTIONS
            sec = getproperty(cfg, sec_sym)
            for f_sym in sort(collect(fieldnames(typeof(sec))))
                val = _canonical_field_repr(getfield(sec, f_sym))
                println(io, "suite.", tag, ".", sec_sym, ".", f_sym, "=", val)
            end
        end
    end
    for name in SUITE_HASHED_CONSTANTS
        println(io, "case.", name, "=", _canonical_field_repr(getfield(@__MODULE__, name)))
    end
    return String(take!(io))
end

"""
    suite_config_sha256()

Return the SHA-256 digest of `suite_config_text()`, the configuration hash of the case set.

Returns
-------
String
    64-character lowercase hexadecimal SHA-256 digest.
"""
suite_config_sha256()::String = bytes2hex(sha256(suite_config_text()))

"""
    expected_suite_config_summary()

Return the configuration summary stored in the `[config]` section of the suite record:
the case list as parallel arrays and the shared settings.

Returns
-------
Dict{String,Any}
    Keyed as in `REQUIRED_CONFIG_KEYS["porous_shell_suite"]`; `k_air_W_mK` holds the air
    conductivity of each run, `air_matched` whether it is the matched air.
"""
function expected_suite_config_summary()::Dict{String,Any}
    cfgs = [suite_case_config(c) for c in SUITE_CASES]
    return Dict{String,Any}(
        "tags" => [String(c[1]) for c in SUITE_CASES],
        "Nx" => [c[2] for c in SUITE_CASES],
        "t_form_myr" => [Float64(c[3]) for c in SUITE_CASES],
        "phi_init" => [Float64(c[4]) for c in SUITE_CASES],
        "n_steps" => [c[5] for c in SUITE_CASES],
        "air_matched" => [c[6] === nothing for c in SUITE_CASES],
        "k_air_W_mK" => [cfg.materials.ksolidm[3] for cfg in cfgs],
        "dt_yr" => SUITE_DT_YR,
        "Nr_reference" => NR_REFERENCE,
        "R_planet_m" => R_PLANET,
        "L_box_m" => L_BOX,
        "seed" => SEED_DEFAULT,
        "porosity_model" => "henke2012",
        "shell_dphi" => DPHI_SHELL,
        "shell_depth_rel_bound" => SHELL_DEPTH_REL_BOUND,
        "tmax_rel_l2_limit" => TMAX_REL_L2_LIMIT,
        "tmax_nx_min" => SUITE_TMAX_NX_MIN,
    )
end

# Whether a recorded config value equals the expected one: floats (also in arrays) within
# a relative 1e-12, other values exactly.
function _summary_value_matches(v_rec, v_exp)::Bool
    if v_exp isa AbstractVector
        return v_rec isa AbstractVector &&
               length(v_rec) == length(v_exp) &&
               all(_summary_value_matches(a, b) for (a, b) in zip(v_rec, v_exp))
    elseif v_exp isa AbstractFloat
        return v_rec isa Real &&
               !(v_rec isa Bool) &&
               isapprox(Float64(v_rec), v_exp; rtol=1.0e-12, atol=0.0)
    end
    return typeof(v_rec) == typeof(v_exp) && v_rec == v_exp
end

"""
    check_suite_record_static(path; check_doc_sync=true)

Check the suite record `path` without a rerun: the schema of `validate_record_dict`, the
record id, the configuration hash of the case set (`suite_config_sha256`), the `[config]`
summary with no key outside it, the case list against `SUITE_CASES`, the metric ids
(`suite_metric_ids`), the case, grid and units (`metric_units`) of the run that each
metric id names, the tolerance of each metric against `suite_metric_limit` and that each
value passes it. The reference of a shell depth must equal the pinned 1-D depth of
its case (`suite_reference_depth`), which also sets its rule, and the shell depth must be
finite and in `[0, R_PLANET]`. A file that is not valid TOML returns the parser error.

Parameters
----------
path : String
    Path to the suite record TOML file.
check_doc_sync : Bool, optional
    Whether the documentation page must hold the table of the record (default `true`).

Returns
-------
Tuple{Bool,Vector{String}}
    Pass flag and diagnostic error messages.

Raises
------
DomainError
    If `isempty(path)`.
"""
function check_suite_record_static(
    path::String; check_doc_sync::Bool=true
)::Tuple{Bool,Vector{String}}
    isempty(path) && throw(DomainError(path, "Record path must be non-empty"))
    isfile(path) || return (false, ["Suite record not found: $path"])
    raw = _parse_record(path)
    raw isa String && return (false, [raw])
    ok_schema, errs = validate_record_dict(raw; require_doc_pages=true, check_doc_sync)
    ok_schema || return (false, errs)
    errors = String[]
    meta = raw["metadata"]
    rec_id = String(meta["id"])
    rec_id == SUITE_RECORD_ID ||
        push!(errors, "Record id '$rec_id' is not '$SUITE_RECORD_ID'")
    expected_sha = suite_config_sha256()
    String(meta["config_sha256"]) == expected_sha || push!(
        errors,
        "Recorded config_sha256 ($(meta["config_sha256"])) does not match the case set ($expected_sha)",
    )
    rec_cfg = raw["config"]
    summary = expected_suite_config_summary()
    _push_extra_config_keys!(errors, rec_cfg, summary)
    for (k, v_exp) in summary
        if !haskey(rec_cfg, k)
            push!(errors, "Recorded config missing key '$k'")
        elseif !_summary_value_matches(rec_cfg[k], v_exp)
            push!(
                errors,
                "Recorded config.$k ($(rec_cfg[k])) does not match expected ($v_exp)",
            )
        end
    end
    metrics = raw["metrics"]
    ids = [String(m["id"]) for m in metrics]
    ids == suite_metric_ids() || push!(
        errors,
        "Recorded metric ids $ids do not match the case list $(suite_metric_ids())",
    )
    nx_of = Dict(String(c[1]) => c[2] for c in SUITE_CASES)
    cases = Set(String(m["case"]) for m in metrics)
    cases == Set(keys(nx_of)) || push!(
        errors,
        "Recorded cases $(sort!(collect(cases))) do not match SUITE_CASES $(sort!(collect(keys(nx_of))))",
    )
    # The run of a metric comes from its id; an id of no run fails the id check above
    tag_of = Dict(
        "$(prefix)_$(c[1])" => String(c[1]) for c in SUITE_CASES for
        prefix in ("shell_depth", "tmax_rel_l2")
    )
    for m in metrics
        mid = String(m["id"])
        haskey(tag_of, mid) || continue
        tag = tag_of[mid]
        String(m["case"]) == tag ||
            push!(errors, "Metric '$mid' has case '$(m["case"])', not '$tag'")
        nx = nx_of[tag]
        Int(m["nx"]) == nx || push!(errors, "Metric '$mid' has nx $(m["nx"]), not $nx")
        _push_units_error!(errors, mid, m)
        value, ref = Float64(m["value"]), Float64(m["reference"])
        if startswith(mid, "shell_depth_")
            (isfinite(value) && 0.0 <= value <= R_PLANET) || push!(
                errors,
                "Shell depth '$mid' = $value is not finite or outside [0, $R_PLANET]",
            )
            pinned = suite_reference_depth(tag)
            isapprox(ref, pinned; rtol=1.0e-12, atol=0.0) || push!(
                errors,
                "Recorded reference of '$mid' ($ref) differs from the pinned 1-D depth ($pinned) of case '$tag'",
            )
            # The rule below uses the pinned depth, so a changed reference cannot move it
            ref = pinned
        end
        lim = try
            suite_metric_limit(mid, nx, ref)
        catch err
            err isa Union{ArgumentError,DomainError} || rethrow()
            push!(errors, "Metric '$mid' has no limit: $(sprint(showerror, err))")
            continue
        end
        atol, rtol, comp = Float64(m["atol"]), Float64(m["rtol"]), String(m["comparison"])
        same_limit =
            comp == lim.comparison &&
            isapprox(atol, lim.atol; rtol=1.0e-12, atol=0.0) &&
            isapprox(rtol, lim.rtol; rtol=1.0e-12, atol=0.0)
        same_limit || push!(
            errors,
            "Recorded tolerance of '$mid' ($comp, atol=$atol, rtol=$rtol) differs from the rule ($(lim.comparison), atol=$(lim.atol), rtol=$(lim.rtol))",
        )
        passes = if lim.comparison == "isapprox"
            0.0 < ref < R_PLANET && isapprox(value, ref; atol=lim.atol, rtol=lim.rtol)
        else
            iszero(ref) && isfinite(value) && 0.0 <= value <= lim.atol
        end
        passes || push!(errors, "Metric '$mid' = $value fails the rule around $ref")
    end
    return (isempty(errors), errors)
end

end # module

module PorousShellValidation

# One-dimensional spherical reference for the porous shell of a planetesimal with the
# dry porosity laws of Henke et al. (2012): conduction with the porous conductivity of
# Lichtenberg et al. (2016a, Eq. 12), 26Al heating, cold pressing and hot sintering.
# The laws are coded here from the raw constants; nothing calls the Erebus porosity code.

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
                        ratio_al=RATIO_AL_L16, T0=T_SPACE) -> NamedTuple

One-dimensional spherical finite-volume run: backward Euler conduction with k, rho cp and
heat from the porosity at the start of the step, Dirichlet `T_SPACE` at `R`, then the dry
porosity step at the new temperature and the pressure of the old porosity (the order of
an Erebus step). At the start the porosity is cold pressed on the hydrostatic profile.
The run stops at `t_end_yr` after formation, or once the peak has passed and every cell
is below `T_SINTER` (no further sintering).

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
        (t > t_peak + dt && Tm < T_SINTER && T_peak >= T_SINTER) && break
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

end # module

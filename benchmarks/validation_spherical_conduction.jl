using Erebus
using StaticArrays
using LinearAlgebra
using Printf
using Logging
using TOML
using SHA

const R_PLANET = 50_000.0
const L_BOX = 140_000.0
const T_SURF = 170.0
const K_COND = 3.0
const PHIM0 = 1.0e-4
const N_STEPS_TRANSIENT = 20
const DT_YR = 500_000.0
const START_TIME_YR = 2.0e6

"""
    make_spherical_conduction_config(; kwargs...) -> SimulationConfig

Construct a pure-conduction spherical planetesimal configuration with all non-conductive physics disabled.

Parameters
----------
Nx : Int, optional
    Number of grid intervals in x and y (default 33).
n_steps : Int, optional
    Number of timesteps (default 20).
start_time : Float64, optional
    Simulation start time after CAI in years (default `2.0e6`).
dt_yr : Float64, optional
    Timestep duration in years (default `500_000.0`).
spherical_metric : Bool, optional
    Enable the 2.5D spherical radial conduction metric (default `true`).
reg_cells : Float64, optional
    Radial metric regularization scale in grid cells (default `0.5`).
disk_enabled : Bool, optional
    Enable sticky-air ambient temperature clamping (default `true`).
k_rock : Float64, optional
    Silicate thermal conductivity in W/(m*K) (default `3.0`).
k_air : Float64, optional
    Sticky-air thermal conductivity in W/(m*K) (default `3.0`).
rhocp_air : Float64, optional
    Sticky-air volumetric heat capacity in J/(m^3*K) (default `3.0e10`).
dsubgridt : Float64, optional
    Subgrid thermal diffusion coefficient (default `0.0`).

Returns
-------
SimulationConfig
    Configured pure-conduction `SimulationConfig`.
"""
function make_spherical_conduction_config(;
    Nx::Int=33,
    n_steps::Int=N_STEPS_TRANSIENT,
    start_time::Float64=START_TIME_YR,
    dt_yr::Float64=DT_YR,
    spherical_metric::Bool=true,
    reg_cells::Float64=0.5,
    disk_enabled::Bool=true,
    k_rock::Float64=K_COND,
    k_air::Float64=K_COND,
    rhocp_air::Float64=3.0e10,
    dsubgridt::Float64=0.0,
)::SimulationConfig
    out_cfg = OutputConfig(;
        savematstep=n_steps + 10, visstep=n_steps + 10, save_final=false
    )
    return SimulationConfig(;
        grid=GridConfig(; xsize=L_BOX, ysize=L_BOX, Nx=Nx, Ny=Nx),
        geometry=GeometryConfig(;
            rplanet=R_PLANET,
            rcrust=R_PLANET,
            xcenter=L_BOX / 2.0,
            ycenter=L_BOX / 2.0,
            spherical_metric=spherical_metric,
            metric_regularization_cells=reg_cells,
        ),
        time=TimeConfig(;
            start_time=start_time,
            endtime=30.0e6,
            dt_initial=dt_yr,
            dt_longest=dt_yr,
            DTmax=500.0,
            n_steps=n_steps,
        ),
        solver=SolverConfig(; dsubgridt=dsubgridt),
        thermodynamics=ThermalConfig(;
            hr_al=true,
            hr_fe=false,
            phim0=PHIM0,
            thermal_buoyancy=false,
            surface_radiation=false,
        ),
        materials=MaterialConfig(;
            ksolidm=SVector{3,Float64}([k_rock, k_rock, k_air]),
            kfluidm=SVector{3,Float64}([k_rock, k_rock, k_air]),
            rhocpsolidm=SVector{3,Float64}([3.3e6, 3.3e6, rhocp_air]),
            rhocpfluidm=SVector{3,Float64}([1.0e6, 1.0e6, rhocp_air]),
            kphim0=SVector{3,Float64}([1.0e-30, 1.0e-30, 1.0e-30]),
            etasolidm=SVector{3,Float64}([1.0e25, 1.0e25, 1.0e22]),
            etasolidmm=SVector{3,Float64}([1.0e25, 1.0e25, 1.0e22]),
            alphasolidm=SVector{3,Float64}([0.0, 0.0, 0.0]),
            alphafluidm=SVector{3,Float64}([0.0, 0.0, 0.0]),
        ),
        poroelasticity=PoroelasticConfig(;
            phimin=PHIM0, phimax=0.9999, hydrofracture=false, betasolid=0.0, betafluid=0.0
        ),
        output=out_cfg,
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
            enabled=disk_enabled, model=:fixed, t_ambient=T_SURF, dispersal_active=false
        ),
    )
end

"""
    canonical_config_sha256(cfg::SimulationConfig) -> String

Compute a deterministic 64-character lowercase hex SHA-256 digest of all fields in `SimulationConfig`.

Parameters
----------
cfg : SimulationConfig
    Simulation configuration to hash.

Returns
-------
String
    64-character lowercase hexadecimal SHA-256 digest.
"""
function canonical_config_sha256(cfg::SimulationConfig)::String
    io = IOBuffer()
    for sec_sym in sort(collect(fieldnames(SimulationConfig)))
        sec = getfield(cfg, sec_sym)
        for f_sym in sort(collect(fieldnames(typeof(sec))))
            val = getfield(sec, f_sym)
            println(io, sec_sym, ".", f_sym, "=", repr(val))
        end
    end
    return bytes2hex(sha256(take!(io)))
end

"""
    expected_config_summary(cfg::SimulationConfig) -> Dict{String,Any}

Return the canonical configuration summary dictionary stored in the `[config]` section of the validation record.

Parameters
----------
cfg : SimulationConfig
    Fine-grid validation configuration.

Returns
-------
Dict{String,Any}
    Dictionary of physical and numerical configuration parameters.
"""
function expected_config_summary(cfg::SimulationConfig)::Dict{String,Any}
    return Dict{String,Any}(
        "Nx_coarse" => 33,
        "Nx_fine" => cfg.grid.Nx,
        "n_steps" => cfg.time.n_steps,
        "dt_yr" => cfg.time.dt_initial,
        "R_planet_m" => cfg.geometry.rplanet,
        "k_W_mK" => cfg.materials.ksolidm[1],
        "rho_kg_m3" => cfg.materials.rhosolidm[1],
        "cp_J_kgK" => cfg.materials.rhocpsolidm[1] / cfg.materials.rhosolidm[1],
        "T_surf_K" => cfg.disk.t_ambient,
        "spherical_metric" => cfg.geometry.spherical_metric,
        "reg_cells" => cfg.geometry.metric_regularization_cells,
        "disk_enabled" => cfg.disk.enabled,
    )
end

const _REF_CFG_65 = make_spherical_conduction_config(; Nx=65)
const X_FE = 0.182
const RHO_SOL = _REF_CFG_65.materials.rhosolidm[1]
const RHO_MET = _REF_CFG_65.coreformation.rho_metal
const RHO_BULK = 1.0 / (X_FE / RHO_MET + (1.0 - X_FE) / RHO_SOL)
const RHOCP_M =
    (1.0 - PHIM0) * _REF_CFG_65.materials.rhocpsolidm[1] +
    PHIM0 * _REF_CFG_65.materials.rhocpfluidm[1]
const K_EFF =
    (1.0 - PHIM0) * _REF_CFG_65.materials.ksolidm[1] +
    PHIM0 * _REF_CFG_65.materials.kfluidm[1]
const KAPPA = K_EFF / RHOCP_M
const T_HALF_AL = _REF_CFG_65.thermodynamics.t_half_al
const TAU_AL = T_HALF_AL / log(2.0)
const LAMBDA_AL = 1.0 / TAU_AL
const T0_S = _REF_CFG_65.time.start_time * _REF_CFG_65.time.yearlength
const DT_S = _REF_CFG_65.time.dt_initial * _REF_CFG_65.time.yearlength

# Compute the initial radiogenic power density at t = 2 Myr independently of src/.
const Q_AL_0 =
    _REF_CFG_65.thermodynamics.f_al *
    _REF_CFG_65.thermodynamics.ratio_al *
    _REF_CFG_65.thermodynamics.E_al *
    exp(-T0_S / TAU_AL) / TAU_AL
const H0_2MYR = (1.0 - PHIM0) * RHO_BULK * Q_AL_0

"""
    analytical_T_BE(r::Real; R_b=R_PLANET, n_steps=N_STEPS_TRANSIENT, dt_s=DT_S, n_modes=500, include_step_average=true) -> Float64

Evaluate the exact time-discrete backward-Euler spherical eigenmode temperature at radius `r` with step-averaged radiogenic heating.

Parameters
----------
r : Real
    Radial coordinate in metres (`r >= 0`).
R_b : Real, optional
    Spherical boundary radius in metres (default `R_PLANET`).
n_steps : Int, optional
    Number of backward-Euler timesteps (default `N_STEPS_TRANSIENT`).
dt_s : Real, optional
    Timestep duration in seconds (default `DT_S`).
n_modes : Int, optional
    Number of spherical Dirichlet eigenmodes (default 500).
include_step_average : Bool, optional
    Include the exact step-integral decay factor `(1 - exp(-dt/tau))/(dt/tau)` (default `true`).

Returns
-------
Float64
    Temperature in K at radius `r`.

Raises
------
DomainError
    If `r < 0`, `R_b <= 0`, `dt_s <= 0`, `n_steps < 1`, or `n_modes < 1`.
"""
function analytical_T_BE(
    r::Real;
    R_b::Real=R_PLANET,
    n_steps::Int=N_STEPS_TRANSIENT,
    dt_s::Real=DT_S,
    n_modes::Int=500,
    include_step_average::Bool=true,
)::Float64
    r < 0.0 && throw(DomainError(r, "Radius r must be non-negative"))
    R_b <= 0.0 && throw(DomainError(R_b, "Boundary radius R_b must be positive"))
    dt_s <= 0.0 && throw(DomainError(dt_s, "Timestep dt_s must be positive"))
    n_steps < 1 && throw(DomainError(n_steps, "Number of steps n_steps must be >= 1"))
    n_modes < 1 && throw(DomainError(n_modes, "Number of eigenmodes n_modes must be >= 1"))
    r >= R_b && return Float64(T_SURF)

    # Compute the exact step-integral decay factor independently of src/.
    x_dt = Float64(dt_s) * LAMBDA_AL
    step_avg = include_step_average ? (1.0 - exp(-x_dt)) / x_dt : 1.0

    q = exp(-LAMBDA_AL * dt_s)
    S0 = (H0_2MYR * step_avg / RHOCP_M) * dt_s
    total = 0.0
    for m in 1:n_modes
        mu_m = KAPPA * (m * pi / R_b)^2
        p = 1.0 / (1.0 + mu_m * dt_s)
        I_m = S0 * p * (p^n_steps - q^n_steps) / (p - q)
        phi_m = if r == 0.0
            (-1.0)^(m + 1) * 2.0
        else
            (-1.0)^(m + 1) * (2.0 * R_b / (m * pi * r)) * sin(m * pi * r / R_b)
        end
        total += phi_m * I_m
    end
    return T_SURF + total
end

"""
    analytical_T_continuous(r::Real; R_b=R_PLANET, t_tot_s=N_STEPS_TRANSIENT * DT_S, n_modes=500) -> Float64

Evaluate the continuous-time spherical eigenmode solution at radius `r` and total elapsed time `t_tot_s`.

Parameters
----------
r : Real
    Radial coordinate in metres (`r >= 0`).
R_b : Real, optional
    Spherical boundary radius in metres (default `R_PLANET`).
t_tot_s : Real, optional
    Total integration time in seconds (default `N_STEPS_TRANSIENT * DT_S`).
n_modes : Int, optional
    Number of spherical Dirichlet eigenmodes (default 500).

Returns
-------
Float64
    Continuous-time analytical temperature in K at radius `r`.

Raises
------
DomainError
    If `r < 0`, `R_b <= 0`, `t_tot_s <= 0`, or `n_modes < 1`.
"""
function analytical_T_continuous(
    r::Real; R_b::Real=R_PLANET, t_tot_s::Real=N_STEPS_TRANSIENT * DT_S, n_modes::Int=500
)::Float64
    r < 0.0 && throw(DomainError(r, "Radius r must be non-negative"))
    R_b <= 0.0 && throw(DomainError(R_b, "Boundary radius R_b must be positive"))
    t_tot_s <= 0.0 && throw(DomainError(t_tot_s, "Total time t_tot_s must be positive"))
    n_modes < 1 && throw(DomainError(n_modes, "Number of eigenmodes n_modes must be >= 1"))
    r >= R_b && return Float64(T_SURF)

    S_c = H0_2MYR / RHOCP_M
    total = 0.0
    for m in 1:n_modes
        mu_m = KAPPA * (m * pi / R_b)^2
        I_m = S_c * (exp(-LAMBDA_AL * t_tot_s) - exp(-mu_m * t_tot_s)) / (mu_m - LAMBDA_AL)
        phi_m = if r == 0.0
            (-1.0)^(m + 1) * 2.0
        else
            (-1.0)^(m + 1) * (2.0 * R_b / (m * pi * r)) * sin(m * pi * r / R_b)
        end
        total += phi_m * I_m
    end
    return T_SURF + total
end

"""
    solve_1d_spherical_fv(Nr::Int=500; n_steps=N_STEPS_TRANSIENT, dt_s=DT_S, include_step_average=true) -> Tuple{Vector{Float64},Vector{Float64}}

Solve the 1D spherical radial heat diffusion equation using a conservative finite-volume tridiagonal discretization.

Parameters
----------
Nr : Int, optional
    Number of radial finite-volume shells (default 500).
n_steps : Int, optional
    Number of backward-Euler timesteps (default `N_STEPS_TRANSIENT`).
dt_s : Real, optional
    Timestep duration in seconds (default `DT_S`).
include_step_average : Bool, optional
    Include the exact step-integral decay factor `(1 - exp(-dt/tau))/(dt/tau)` (default `true`).

Returns
-------
Tuple{Vector{Float64},Vector{Float64}}
    Cell-centre radii `r_c` (m) and final temperatures `T` (K).

Raises
------
DomainError
    If `Nr < 2`, `dt_s <= 0`, or `n_steps < 1`.
"""
function solve_1d_spherical_fv(
    Nr::Int=500;
    n_steps::Int=N_STEPS_TRANSIENT,
    dt_s::Real=DT_S,
    include_step_average::Bool=true,
)::Tuple{Vector{Float64},Vector{Float64}}
    Nr < 2 && throw(DomainError(Nr, "Number of radial shells Nr must be >= 2"))
    dt_s <= 0.0 && throw(DomainError(dt_s, "Timestep dt_s must be positive"))
    n_steps < 1 && throw(DomainError(n_steps, "Number of steps n_steps must be >= 1"))

    dr = R_PLANET / Nr
    r_f = [i * dr for i in 0:Nr]
    r_c = [(i - 0.5) * dr for i in 1:Nr]
    V = [(4.0 * pi / 3.0) * (r_f[i + 1]^3 - r_f[i]^3) for i in 1:Nr]
    A = [4.0 * pi * r_f[i]^2 for i in 1:(Nr + 1)]

    dl = zeros(Float64, Nr - 1)
    d = zeros(Float64, Nr)
    du = zeros(Float64, Nr - 1)

    for i in 1:Nr
        cond_in = i == 1 ? 0.0 : A[i] * K_EFF / dr
        cond_out = i == Nr ? A[Nr + 1] * K_EFF / (0.5 * dr) : A[i + 1] * K_EFF / dr
        d[i] = RHOCP_M * V[i] / dt_s + cond_in + cond_out
        if i > 1
            dl[i - 1] = -cond_in
        end
        if i < Nr
            du[i] = -cond_out
        end
    end
    M = Tridiagonal(dl, d, du)

    x_dt = Float64(dt_s) * LAMBDA_AL
    step_avg = include_step_average ? (1.0 - exp(-x_dt)) / x_dt : 1.0

    T = fill(T_SURF, Nr)
    rhs = zeros(Float64, Nr)
    cond_R = A[Nr + 1] * K_EFF / (0.5 * dr)
    for k in 0:(n_steps - 1)
        H_k = H0_2MYR * exp(-k * LAMBDA_AL * dt_s) * step_avg
        for i in 1:Nr
            rhs[i] = (RHOCP_M * V[i] / dt_s) * T[i] + H_k * V[i]
        end
        rhs[Nr] += cond_R * T_SURF
        T = M \ rhs
    end
    return r_c, T
end

"""
    run_steady_spherical_conduction(Nx::Int; spherical_metric=true, metric_scale=1.0, clamp_exterior=true, perturb_stencil=false)

Solve the steady-state spherical Poisson conduction problem using `assemble_thermal_lse!` and `compute_spherical_metric_heat_source!`.
"""
function run_steady_spherical_conduction(
    Nx::Int;
    spherical_metric::Bool=true,
    metric_scale::Float64=1.0,
    clamp_exterior::Bool=true,
    perturb_stencil::Bool=false,
)
    coords = GridCoordinates(Nx, Nx; xsize=L_BOX, ysize=L_BOX)
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    dx = coords.dx
    dy = coords.dy
    cx = L_BOX / 2.0
    cy = L_BOX / 2.0

    H_val = 1.0e-4
    dt_steady = clamp_exterior ? 1.0 : 1.0e14

    KX = fill(K_COND, Ny1, Nx1)
    KY = fill(K_COND, Ny1, Nx1)
    T0 = fill(T_SURF, Ny1, Nx1)
    RHOCP = zeros(Float64, Ny1, Nx1)
    HR = zeros(Float64, Ny1, Nx1)
    HA = zeros(Float64, Ny1, Nx1)
    HS = zeros(Float64, Ny1, Nx1)
    DHP = zeros(Float64, Ny1, Nx1)
    Q_metric = zeros(Float64, Ny1, Nx1)

    for j in 1:Nx1, i in 1:Ny1
        r = hypot(coords.xp[j] - cx, coords.yp[i] - cy)
        if r <= R_PLANET
            HR[i, j] = H_val
            RHOCP[i, j] = clamp_exterior ? 0.0 : 1.0e6
        else
            HR[i, j] = 0.0
            RHOCP[i, j] = clamp_exterior ? 1.0e4 : 1.0e6
        end
    end

    RT_2d = zeros(Float64, Ny1 * Nx1)
    LT_2d = Erebus.assemble_thermal_lse!(
        T0, RHOCP, KX, KY, HR, HA, HS, DHP, RT_2d, dt_steady; coords=coords
    )
    tk_2d = reshape(LT_2d \ RT_2d, Ny1, Nx1)

    tk_iter = copy(tk_2d)
    RT_sph = zeros(Float64, Ny1 * Nx1)
    omega = 0.5
    n_iter = spherical_metric ? 50 : 1
    for _ in 1:n_iter
        if spherical_metric
            Erebus.compute_spherical_metric_heat_source!(
                Q_metric,
                tk_iter,
                KX,
                KY,
                coords;
                xcenter=cx,
                ycenter=cy,
                rplanet=R_PLANET,
                reg_cells=0.5,
            )
            if metric_scale != 1.0
                Q_metric .*= metric_scale
            end
        else
            fill!(Q_metric, 0.0)
        end
        LT_sph = Erebus.assemble_thermal_lse!(
            T0,
            RHOCP,
            KX,
            KY,
            HR,
            HA,
            HS,
            DHP,
            RT_sph,
            dt_steady;
            coords=coords,
            Q_metric=Q_metric,
        )
        tk_new = reshape(LT_sph \ RT_sph, Ny1, Nx1)
        tk_iter .= omega .* tk_new .+ (1.0 - omega) .* tk_iter
    end

    RHOCP_energy = fill(RHOCP_M, Ny1, Nx1)
    dt_energy = 1.0e11
    RT_e = zeros(Float64, Ny1 * Nx1)
    LT_e = Erebus.assemble_thermal_lse!(
        T0,
        RHOCP_energy,
        KX,
        KY,
        HR,
        HA,
        HS,
        DHP,
        RT_e,
        dt_energy;
        coords=coords,
        Q_metric=Q_metric,
    )
    if perturb_stencil
        mid_idx = ((Nx1 ÷ 2) - 1) * Ny1 + (Ny1 ÷ 2)
        LT_e[mid_idx, mid_idx] *= 1.01
    end
    tk_e = reshape(LT_e \ RT_e, Ny1, Nx1)

    dE_grid = 0.0
    W_src = 0.0
    for j in 2:(Nx1 - 1), i in 2:(Ny1 - 1)
        dE_grid += RHOCP_energy[i, j] * (tk_e[i, j] - T0[i, j]) * dx * dy
        W_src +=
            dt_energy *
            (HR[i, j] + HA[i, j] + HS[i, j] + DHP[i, j] + Q_metric[i, j]) *
            dx *
            dy
    end
    energy_res = abs(dE_grid - W_src) / abs(dE_grid)

    r_min = Inf
    dT_center_3d = 0.0
    dT_center_2d = 0.0
    sq_err_3d = 0.0
    sq_ref_3d = 0.0
    for j in 2:(Nx1 - 1), i in 2:(Ny1 - 1)
        r = hypot(coords.xp[j] - cx, coords.yp[i] - cy)
        if r <= R_PLANET
            dT_num = tk_iter[i, j] - T_SURF
            dT_2d_node = tk_2d[i, j] - T_SURF
            dT_ana = H_val * (R_PLANET^2 - r^2) / (6.0 * K_COND)
            sq_err_3d += (dT_num - dT_ana)^2
            sq_ref_3d += dT_ana^2
            if r < r_min
                r_min = r
                dT_center_3d = dT_num
                dT_center_2d = dT_2d_node
            end
        end
    end
    ratio_3d_2d = dT_center_3d / dT_center_2d
    l2_rel = sqrt(sq_err_3d / sq_ref_3d)

    return (; Nx, dx, r_min, dT_center_3d, dT_center_2d, ratio_3d_2d, l2_rel, energy_res)
end

"""
    run_transient_spherical_conduction(Nx::Int; kwargs...)

Run the transient 2D spherical radiogenic conduction simulation through `simulation_loop` and compare against `analytical_T_BE`.
"""
function run_transient_spherical_conduction(
    Nx::Int;
    n_steps::Int=N_STEPS_TRANSIENT,
    dt_yr::Float64=DT_YR,
    spherical_metric::Bool=true,
    reg_cells::Float64=0.5,
    k_rock::Float64=K_COND,
    k_air::Float64=K_COND,
    disk_enabled::Bool=true,
    rhocp_air::Float64=3.0e10,
    dsubgridt::Float64=0.0,
)
    cfg = make_spherical_conduction_config(;
        Nx=Nx,
        n_steps=n_steps,
        dt_yr=dt_yr,
        spherical_metric=spherical_metric,
        reg_cells=reg_cells,
        disk_enabled=disk_enabled,
        k_rock=k_rock,
        k_air=k_air,
        rhocp_air=rhocp_air,
        dsubgridt=dsubgridt,
    )
    st = mktempdir() do tmpdir
        redirect_stdout(devnull) do
            redirect_stderr(devnull) do
                with_logger(NullLogger()) do
                    return simulation_loop(cfg; output_path=tmpdir)
                end
            end
        end
    end

    coords = GridCoordinates(cfg.grid)
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    cx = L_BOX / 2.0
    cy = L_BOX / 2.0
    dt_s_local = dt_yr * cfg.time.yearlength

    r_min = Inf
    T_center = 0.0
    sq_err = 0.0
    sq_ref = 0.0
    T_max_air_shell = -Inf
    for j in 1:Nx1, i in 1:Ny1
        r = hypot(coords.xp[j] - cx, coords.yp[i] - cy)
        if r <= R_PLANET
            T_num = st.grids.tk2[i, j]
            T_ref = analytical_T_BE(r; n_steps=n_steps, dt_s=dt_s_local)
            sq_err += ((T_num - T_SURF) - (T_ref - T_SURF))^2
            sq_ref += (T_ref - T_SURF)^2
            if r < r_min
                r_min = r
                T_center = T_num
            end
        elseif r >= R_PLANET + coords.dx
            T_max_air_shell = max(T_max_air_shell, st.grids.tk2[i, j])
        end
    end

    dT_center = T_center - T_SURF
    dT_be_0 = analytical_T_BE(0.0; n_steps=n_steps, dt_s=dt_s_local) - T_SURF
    dT_be_rmin = analytical_T_BE(r_min; n_steps=n_steps, dt_s=dt_s_local) - T_SURF
    center_err_K = dT_center - dT_be_rmin
    rel_err_be_rmin = abs(center_err_K) / dT_be_rmin
    l2_rel = sqrt(sq_err / sq_ref)

    return (;
        Nx,
        r_min,
        T_center,
        dT_center,
        dT_be_0,
        dT_be_rmin,
        center_err_K,
        rel_err_be_rmin,
        l2_rel,
        T_max_air_shell,
    )
end

"""
    evaluate_spherical_conduction_benchmark()

Run the steady and transient spherical conduction benchmarks and return all 10 validation metrics.
"""
function evaluate_spherical_conduction_benchmark()
    s33 = run_steady_spherical_conduction(33)
    s65 = run_steady_spherical_conduction(65)
    t33 = run_transient_spherical_conduction(33)
    t65 = run_transient_spherical_conduction(65)

    r_fv, T_fv = solve_1d_spherical_fv(500)
    T_be_0 = analytical_T_BE(0.0; n_modes=500)
    dT_be_0 = T_be_0 - T_SURF
    dT_fv_r1 = T_fv[1] - T_SURF
    be_vs_1dfv_rel_diff = abs(dT_be_0 - dT_fv_r1) / dT_be_0

    dT_500 = analytical_T_BE(0.0; n_modes=500) - T_SURF
    dT_1000 = analytical_T_BE(0.0; n_modes=1000) - T_SURF
    mode_trunc_500_vs_1000 = abs(dT_500 - dT_1000) / dT_1000

    metrics = Dict{String,Float64}(
        "steady_center_dT_ratio_3d_2d_nx33" => s33.ratio_3d_2d,
        "steady_center_dT_ratio_3d_2d_nx65" => s65.ratio_3d_2d,
        "steady_profile_l2_rel_error_nx33" => s33.l2_rel,
        "steady_profile_l2_rel_error_nx65" => s65.l2_rel,
        "steady_energy_conservation_residual_nx65" => s65.energy_res,
        "transient_be_vs_1dfv_rel_diff" => be_vs_1dfv_rel_diff,
        "transient_center_dT_nx65" => t65.dT_center,
        "transient_profile_l2_rel_error_nx65" => t65.l2_rel,
        "transient_profile_l2_convergence_ratio_65_33" => t65.l2_rel / t33.l2_rel,
        "transient_max_air_shell_temperature_nx65" => t65.T_max_air_shell,
    )

    return (; s33, s65, t33, t65, be_vs_1dfv_rel_diff, mode_trunc_500_vs_1000, metrics)
end

"""
    check_validation_reproducibility(bench, record_path::String; check_recorded_values::Bool=(Threads.nthreads() == 1)) -> Tuple{Bool,Vector{String}}

Compare freshly evaluated benchmark metrics and configuration hashes against a committed TOML validation record.
"""
function check_validation_reproducibility(
    bench, record_path::String; check_recorded_values::Bool=(Threads.nthreads() == 1)
)::Tuple{Bool,Vector{String}}
    if !isfile(record_path)
        return (false, ["Validation record not found: $record_path"])
    end
    raw = TOML.parsefile(record_path)
    ok_schema, errs = validate_record_dict(raw; require_doc_pages=false)
    errors = copy(errs)

    cfg65 = make_spherical_conduction_config(; Nx=65)
    expected_sha = canonical_config_sha256(cfg65)
    recorded_sha = String(get(get(raw, "metadata", Dict()), "config_sha256", ""))
    if recorded_sha != expected_sha
        push!(
            errors,
            "Recorded config_sha256 ($recorded_sha) does not match canonical config SHA-256 ($expected_sha)",
        )
    end

    expected_cfg = expected_config_summary(cfg65)
    rec_cfg = get(raw, "config", Dict{String,Any}())
    for (k, v_exp) in expected_cfg
        if !haskey(rec_cfg, k) || rec_cfg[k] != v_exp
            push!(
                errors,
                "Recorded config.$k ($(get(rec_cfg, k, nothing))) does not match expected ($v_exp)",
            )
        end
    end

    for m in get(raw, "metrics", Any[])
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

        tol_ok = if comp == "isapprox"
            isapprox(measured, ref; atol=atol, rtol=rtol)
        elseif comp == "le_atol"
            isfinite(measured) && measured <= atol
        else
            false
        end
        if !tol_ok
            push!(
                errors,
                "Measured metric '$mid' = $measured fails tolerance ($comp, ref=$ref, atol=$atol, rtol=$rtol)",
            )
        end
        if check_recorded_values &&
            !isapprox(measured, recorded; rtol=1.0e-10, atol=1.0e-12)
            push!(
                errors,
                "Measured metric '$mid' = $measured differs from recorded value $recorded",
            )
        end
    end

    return (ok_schema && isempty(errors), errors)
end

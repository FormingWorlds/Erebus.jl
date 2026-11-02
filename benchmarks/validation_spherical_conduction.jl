module SphericalConductionValidation

using Erebus
using StaticArrays
using LinearAlgebra
using Logging
using TOML
using SHA

include(joinpath(@__DIR__, "..", "tools", "validation_records.jl"))

const R_PLANET = 50_000.0
const L_BOX = 140_000.0
const T_SURF = 170.0
const K_COND = 3.0
const PHIM0 = 1.0e-4
const REG_CELLS_DEFAULT = 0.5
const H_STEADY = 1.0e-4
const N_STEPS_TRANSIENT = 20
const DT_YR = 500_000.0
const START_TIME_YR = 2.0e6
const YEARLENGTH_S = 31_557_600.0

const X_FE = 0.182
const RHO_SOL = 3300.0
const RHO_MET = 5450.0
const RHOCP_SOL = 3.3e6
const RHOCP_FLU = 1.0e6
const RHOCP_AIR = 3.0e10
const T_HALF_AL = 717_000.0 * 31_540_000.0
const F_AL = 1.9e23
const RATIO_AL = 5.0e-5
const E_AL = 5.0470e-13

const RHO_BULK = 1.0 / (X_FE / RHO_MET + (1.0 - X_FE) / RHO_SOL)
const RHOCP_M = (1.0 - PHIM0) * RHOCP_SOL + PHIM0 * RHOCP_FLU
const K_EFF = (1.0 - PHIM0) * K_COND + PHIM0 * K_COND
const KAPPA = K_EFF / RHOCP_M
const TAU_AL = T_HALF_AL / log(2.0)
const LAMBDA_AL = 1.0 / TAU_AL
const T0_S = START_TIME_YR * YEARLENGTH_S
const DT_S = DT_YR * YEARLENGTH_S
const Q_AL_0 = F_AL * RATIO_AL * E_AL * exp(-T0_S / TAU_AL) / TAU_AL
const H0_2MYR = (1.0 - PHIM0) * RHO_BULK * Q_AL_0

const BENCHMARK_CONFIG_SECTIONS = (
    :grid, :geometry, :time, :thermodynamics, :materials, :poroelasticity, :disk
)

const EXPECTED_METRIC_IDS = (
    "steady_center_dT_ratio_3d_2d_nx33",
    "steady_center_dT_ratio_3d_2d_nx65",
    "steady_profile_l2_rel_error_nx33",
    "steady_profile_l2_rel_error_nx65",
    "steady_energy_conservation_residual_nx65",
    "transient_be_vs_1dfv_rel_diff",
    "transient_center_dT_nx65",
    "transient_profile_l2_rel_error_nx65",
    "transient_profile_l2_convergence_ratio_65_33",
    "transient_max_air_shell_temperature_nx65",
)

const PHYSICS_REGRESSION_BANDS = Dict{
    String,NamedTuple{(:atol, :rtol),Tuple{Float64,Float64}}
}(
    "steady_center_dT_ratio_3d_2d_nx33" => (; atol=0.01, rtol=0.0),
    "steady_center_dT_ratio_3d_2d_nx65" => (; atol=0.01, rtol=0.0),
    "steady_profile_l2_rel_error_nx33" => (; atol=0.005, rtol=0.0),
    "steady_profile_l2_rel_error_nx65" => (; atol=0.005, rtol=0.0),
    "steady_energy_conservation_residual_nx65" => (; atol=1.0e-11, rtol=0.0),
    "transient_center_dT_nx65" => (; atol=0.0, rtol=0.01),
    "transient_profile_l2_rel_error_nx65" => (; atol=0.005, rtol=0.0),
    "transient_profile_l2_convergence_ratio_65_33" => (; atol=0.02, rtol=0.0),
    "transient_max_air_shell_temperature_nx65" => (; atol=0.001, rtol=0.0),
)

const REPLAY_RTOL = 1.0e-6
const REPLAY_NOISE_ATOLS = Dict{String,Float64}(
    "steady_energy_conservation_residual_nx65" => 1.0e-12,
    "transient_be_vs_1dfv_rel_diff" => 5.0e-7,
)

"""
    make_spherical_conduction_config(; Nx=33, n_steps=N_STEPS_TRANSIENT, start_time=START_TIME_YR, dt_yr=DT_YR, spherical_metric=true, reg_cells=REG_CELLS_DEFAULT, disk_enabled=true, k_rock=K_COND, k_air=K_COND, rhocp_air=RHOCP_AIR, dsubgridt=0.0, f_al=F_AL) -> SimulationConfig

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
f_al : Float64, optional
    Chondritic Al abundance in atoms/kg (default `1.9e23`).

Returns
-------
SimulationConfig
    Configured pure-conduction `SimulationConfig`.

Raises
------
DomainError
    If `Nx < 5`, `n_steps < 1`, `dt_yr <= 0`, `reg_cells < 0`, `k_rock <= 0`, `k_air <= 0`, `rhocp_air <= 0`, or `f_al <= 0`.
"""
function make_spherical_conduction_config(;
    Nx::Int=33,
    n_steps::Int=N_STEPS_TRANSIENT,
    start_time::Float64=START_TIME_YR,
    dt_yr::Float64=DT_YR,
    spherical_metric::Bool=true,
    reg_cells::Float64=REG_CELLS_DEFAULT,
    disk_enabled::Bool=true,
    k_rock::Float64=K_COND,
    k_air::Float64=K_COND,
    rhocp_air::Float64=RHOCP_AIR,
    dsubgridt::Float64=0.0,
    f_al::Float64=F_AL,
)::SimulationConfig
    Nx < 5 && throw(DomainError(Nx, "Grid resolution Nx must be >= 5"))
    n_steps < 1 && throw(DomainError(n_steps, "Number of steps n_steps must be >= 1"))
    dt_yr <= 0.0 && throw(DomainError(dt_yr, "Timestep dt_yr must be positive"))
    reg_cells < 0.0 &&
        throw(DomainError(reg_cells, "Regularization cells reg_cells must be >= 0"))
    k_rock <= 0.0 && throw(DomainError(k_rock, "Rock conductivity k_rock must be positive"))
    k_air <= 0.0 && throw(DomainError(k_air, "Air conductivity k_air must be positive"))
    rhocp_air <= 0.0 &&
        throw(DomainError(rhocp_air, "Air heat capacity rhocp_air must be positive"))
    f_al <= 0.0 && throw(DomainError(f_al, "Aluminium abundance f_al must be positive"))

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
            f_al=f_al,
            phim0=PHIM0,
            thermal_buoyancy=false,
            surface_radiation=false,
        ),
        materials=MaterialConfig(;
            ksolidm=SVector{3,Float64}([k_rock, k_rock, k_air]),
            kfluidm=SVector{3,Float64}([k_rock, k_rock, k_air]),
            rhocpsolidm=SVector{3,Float64}([RHOCP_SOL, RHOCP_SOL, rhocp_air]),
            rhocpfluidm=SVector{3,Float64}([RHOCP_FLU, RHOCP_FLU, rhocp_air]),
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
    check_reference_literals_match_config(cfg::SimulationConfig=make_spherical_conduction_config(; Nx=65)) -> Tuple{Bool,Vector{String}}

Compare every literal physical constant used by the analytical and 1D finite-volume reference
solutions against the corresponding `SimulationConfig` fields so that a changed default fails
with an explicit diagnostic message.

Parameters
----------
cfg : SimulationConfig, optional
    Simulation configuration to compare against the benchmark literal constants.

Returns
-------
Tuple{Bool,Vector{String}}
    Boolean pass flag and list of mismatch descriptions.
"""
function check_reference_literals_match_config(
    cfg::SimulationConfig=make_spherical_conduction_config(; Nx=65)
)::Tuple{Bool,Vector{String}}
    checks = (
        ("R_PLANET", R_PLANET, cfg.geometry.rplanet),
        ("L_BOX", L_BOX, cfg.grid.xsize),
        ("T_SURF", T_SURF, cfg.disk.t_ambient),
        ("K_COND", K_COND, cfg.materials.ksolidm[1]),
        ("PHIM0", PHIM0, cfg.thermodynamics.phim0),
        ("REG_CELLS_DEFAULT", REG_CELLS_DEFAULT, cfg.geometry.metric_regularization_cells),
        ("N_STEPS_TRANSIENT", Float64(N_STEPS_TRANSIENT), Float64(cfg.time.n_steps)),
        ("DT_YR", DT_YR, cfg.time.dt_initial),
        ("START_TIME_YR", START_TIME_YR, cfg.time.start_time),
        ("YEARLENGTH_S", YEARLENGTH_S, cfg.time.yearlength),
        ("RHO_SOL", RHO_SOL, cfg.materials.rhosolidm[1]),
        ("RHO_MET", RHO_MET, cfg.coreformation.rho_metal),
        ("RHOCP_SOL", RHOCP_SOL, cfg.materials.rhocpsolidm[1]),
        ("RHOCP_FLU", RHOCP_FLU, cfg.materials.rhocpfluidm[1]),
        ("RHOCP_AIR", RHOCP_AIR, cfg.materials.rhocpsolidm[3]),
        ("T_HALF_AL", T_HALF_AL, cfg.thermodynamics.t_half_al),
        ("F_AL", F_AL, cfg.thermodynamics.f_al),
        ("RATIO_AL", RATIO_AL, cfg.thermodynamics.ratio_al),
        ("E_AL", E_AL, cfg.thermodynamics.E_al),
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

function _canonical_field_repr(val)::String
    if val isa AbstractDict
        isempty(val) && return repr(val)
        sorted_pairs = sort!(collect(pairs(val)); by=p -> string(first(p)))
        inner = join(
            ["$(repr(k)) => $(_canonical_field_repr(v))" for (k, v) in sorted_pairs], ", "
        )
        return "$(typeof(val))($inner)"
    end
    return repr(val)
end

"""
    canonical_config_sha256(cfg; H_steady=H_STEADY, L_box=L_BOX, reg_cells=REG_CELLS_DEFAULT, k_cond=K_COND) -> String

Compute a deterministic 64-character lowercase hex SHA-256 digest of the configuration sections
read by the spherical conduction benchmark (`:grid`, `:geometry`, `:time`, `:thermodynamics`,
`:materials`, `:poroelasticity`, `:disk`) together with the steady-state parameters `H_steady`,
`L_box`, `reg_cells`, and `k_cond`. Unrelated configuration sections are excluded so that adding
a field to an unused physics module does not invalidate the benchmark hash.

Parameters
----------
cfg : Union{SimulationConfig,NamedTuple}
    Simulation configuration (or section container) holding the 7 benchmark sections.
H_steady : Float64, optional
    Steady-state volumetric heat source in W/m^3 (default `H_STEADY`).
L_box : Float64, optional
    Cartesian box side length in m (default `L_BOX`).
reg_cells : Float64, optional
    Radial metric regularization width in grid cells (default `REG_CELLS_DEFAULT`).
k_cond : Float64, optional
    Steady-state thermal conductivity in W/(m*K) (default `K_COND`).

Returns
-------
String
    64-character lowercase hexadecimal SHA-256 digest.

Raises
------
DomainError
    If `H_steady <= 0`, `L_box <= 0`, `reg_cells < 0`, or `k_cond <= 0`.
"""
function canonical_config_sha256(
    cfg;
    H_steady::Float64=H_STEADY,
    L_box::Float64=L_BOX,
    reg_cells::Float64=REG_CELLS_DEFAULT,
    k_cond::Float64=K_COND,
)::String
    H_steady <= 0.0 &&
        throw(DomainError(H_steady, "Steady heat source H_steady must be positive"))
    L_box <= 0.0 && throw(DomainError(L_box, "Box size L_box must be positive"))
    reg_cells < 0.0 &&
        throw(DomainError(reg_cells, "Regularization cells reg_cells must be >= 0"))
    k_cond <= 0.0 && throw(DomainError(k_cond, "Thermal conductivity k_cond must be > 0"))

    io = IOBuffer()
    for sec_sym in BENCHMARK_CONFIG_SECTIONS
        sec = getproperty(cfg, sec_sym)
        for f_sym in sort(collect(fieldnames(typeof(sec))))
            val = getfield(sec, f_sym)
            println(io, sec_sym, ".", f_sym, "=", _canonical_field_repr(val))
        end
    end
    println(io, "steady.H_steady_W_m3=", _canonical_field_repr(H_steady))
    println(io, "steady.K_COND=", _canonical_field_repr(k_cond))
    println(io, "steady.L_BOX=", _canonical_field_repr(L_box))
    println(io, "steady.reg_cells=", _canonical_field_repr(reg_cells))
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
        "Nx_fine" => cfg.grid.Nx,
        "n_steps" => cfg.time.n_steps,
        "dt_yr" => cfg.time.dt_initial,
        "R_planet_m" => cfg.geometry.rplanet,
        "L_box_m" => cfg.grid.xsize,
        "k_W_mK" => cfg.materials.ksolidm[1],
        "rho_kg_m3" => cfg.materials.rhosolidm[1],
        "cp_J_kgK" => cfg.materials.rhocpsolidm[1] / cfg.materials.rhosolidm[1],
        "T_surf_K" => cfg.disk.t_ambient,
        "H_steady_W_m3" => H_STEADY,
        "spherical_metric" => cfg.geometry.spherical_metric,
        "reg_cells" => cfg.geometry.metric_regularization_cells,
        "disk_enabled" => cfg.disk.enabled,
    )
end

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

    # Exact integral of exponential decay over a single backward-Euler step [t_n, t_n + dt].
    x_dt = Float64(dt_s) * LAMBDA_AL
    step_avg = include_step_average ? (1.0 - exp(-x_dt)) / x_dt : 1.0

    q = exp(-LAMBDA_AL * dt_s)
    S0 = (H0_2MYR * step_avg / RHOCP_M) * dt_s
    total = 0.0
    for m in 1:n_modes
        mu_m = KAPPA * (m * pi / R_b)^2
        p = 1.0 / (1.0 + mu_m * dt_s)
        I_m = if abs(p - q) < 1.0e-12 * max(p, q)
            S0 * p * n_steps * (0.5 * (p + q))^(n_steps - 1)
        else
            S0 * p * (p^n_steps - q^n_steps) / (p - q)
        end
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
        I_m = if abs(mu_m - LAMBDA_AL) < 1.0e-12 * max(mu_m, LAMBDA_AL)
            S_c * t_tot_s * exp(-0.5 * (mu_m + LAMBDA_AL) * t_tot_s)
        else
            S_c * (exp(-LAMBDA_AL * t_tot_s) - exp(-mu_m * t_tot_s)) / (mu_m - LAMBDA_AL)
        end
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
    run_steady_spherical_conduction(Nx::Int; spherical_metric=true, metric_scale=1.0, clamp_exterior=true, perturb_stencil=false) -> NamedTuple

Solve the steady-state spherical Poisson conduction problem using `assemble_thermal_lse!` and `compute_spherical_metric_heat_source!`.

Parameters
----------
Nx : Int
    Number of grid intervals in x and y (`Nx >= 5`).
spherical_metric : Bool, optional
    Include the 2.5D spherical radial conduction metric source term (default `true`).
metric_scale : Float64, optional
    Multiplicative factor applied to `Q_metric` (default `1.0`).
clamp_exterior : Bool, optional
    Clamp exterior sticky-air nodes to `T_SURF` via high thermal inertia (default `true`).
perturb_stencil : Bool, optional
    Perturb one interior diagonal entry of the energy-check matrix by `1%` (default `false`).

Returns
-------
NamedTuple
    Steady-state centre temperature rises, 3D/2D ratio, relative L2 error, and energy residual.

Raises
------
DomainError
    If `Nx < 5`.
"""
function run_steady_spherical_conduction(
    Nx::Int;
    spherical_metric::Bool=true,
    metric_scale::Float64=1.0,
    clamp_exterior::Bool=true,
    perturb_stencil::Bool=false,
)
    Nx < 5 && throw(DomainError(Nx, "Grid resolution Nx must be >= 5"))
    coords = GridCoordinates(Nx, Nx; xsize=L_BOX, ysize=L_BOX)
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    dx = coords.dx
    dy = coords.dy
    cx = L_BOX / 2.0
    cy = L_BOX / 2.0

    H_val = H_STEADY
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
                reg_cells=REG_CELLS_DEFAULT,
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
    run_transient_spherical_conduction(Nx::Int; n_steps=N_STEPS_TRANSIENT, dt_yr=DT_YR, spherical_metric=true, reg_cells=REG_CELLS_DEFAULT, k_rock=K_COND, k_air=K_COND, disk_enabled=true, rhocp_air=RHOCP_AIR, dsubgridt=0.0, f_al=F_AL) -> NamedTuple

Run the transient 2D spherical radiogenic conduction simulation through `simulation_loop` and compare against `analytical_T_BE`.

Parameters
----------
Nx : Int
    Number of grid intervals in x and y (`Nx >= 5`).
n_steps : Int, optional
    Number of timesteps (default `N_STEPS_TRANSIENT`).
dt_yr : Float64, optional
    Timestep duration in years (default `DT_YR`).
spherical_metric : Bool, optional
    Enable the 2.5D spherical radial conduction metric (default `true`).
reg_cells : Float64, optional
    Radial metric regularization scale in grid cells (default `REG_CELLS_DEFAULT`).
k_rock : Float64, optional
    Silicate thermal conductivity in W/(m*K) (default `K_COND`).
k_air : Float64, optional
    Sticky-air thermal conductivity in W/(m*K) (default `K_COND`).
disk_enabled : Bool, optional
    Enable sticky-air ambient temperature clamping (default `true`).
rhocp_air : Float64, optional
    Sticky-air volumetric heat capacity in J/(m^3*K) (default `RHOCP_AIR`).
dsubgridt : Float64, optional
    Subgrid thermal diffusion coefficient (default `0.0`).
f_al : Float64, optional
    Simulated chondritic Al abundance in atoms/kg (default `F_AL`).

Returns
-------
NamedTuple
    Transient centre temperature rise, analytical reference values, relative L2 error, and exterior air temperature.

Raises
------
DomainError
    If `Nx < 5`, `n_steps < 1`, `dt_yr <= 0`, `reg_cells < 0`, `k_rock <= 0`, `k_air <= 0`, `rhocp_air <= 0`, or `f_al <= 0`.
"""
function run_transient_spherical_conduction(
    Nx::Int;
    n_steps::Int=N_STEPS_TRANSIENT,
    dt_yr::Float64=DT_YR,
    spherical_metric::Bool=true,
    reg_cells::Float64=REG_CELLS_DEFAULT,
    k_rock::Float64=K_COND,
    k_air::Float64=K_COND,
    disk_enabled::Bool=true,
    rhocp_air::Float64=RHOCP_AIR,
    dsubgridt::Float64=0.0,
    f_al::Float64=F_AL,
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
        f_al=f_al,
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
    dt_s_local = dt_yr * YEARLENGTH_S

    r_min = Inf
    T_center = 0.0
    sq_err = 0.0
    sq_ref = 0.0
    T_max_air_shell = -Inf
    ref_cache = Dict{Float64,Float64}()
    for j in 1:Nx1, i in 1:Ny1
        r = hypot(coords.xp[j] - cx, coords.yp[i] - cy)
        if r <= R_PLANET
            T_num = st.grids.tk2[i, j]
            T_ref = get!(ref_cache, r) do
                return analytical_T_BE(r; n_steps=n_steps, dt_s=dt_s_local)
            end
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
    dT_be_rmin = ref_cache[r_min] - T_SURF
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
    evaluate_spherical_conduction_benchmark() -> NamedTuple

Run the steady and transient spherical conduction benchmarks and return all 10 validation metrics.

Returns
-------
NamedTuple
    Steady and transient benchmark outputs, 1D finite-volume comparison at `r_fv[1]`, mode-truncation residual, and `metrics` dictionary.
"""
function evaluate_spherical_conduction_benchmark()
    s33 = run_steady_spherical_conduction(33)
    s65 = run_steady_spherical_conduction(65)
    t33 = run_transient_spherical_conduction(33)
    t65 = run_transient_spherical_conduction(65)

    r_fv, T_fv = solve_1d_spherical_fv(500)
    dT_be_r1 = analytical_T_BE(r_fv[1]; n_modes=500) - T_SURF
    dT_fv_r1 = T_fv[1] - T_SURF
    be_vs_1dfv_rel_diff = abs(dT_be_r1 - dT_fv_r1) / dT_be_r1

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
    expected_metric_references(bench) -> Dict{String,Float64}

Recompute the expected `reference` value for each of the 10 validation metrics from the
analytical reference functions and constants (`0.0` for `le_atol` upper-bound metrics).

Parameters
----------
bench : NamedTuple
    Evaluated benchmark summary returned by `evaluate_spherical_conduction_benchmark`.

Returns
-------
Dict{String,Float64}
    Mapping from metric `id` to its exact reference value.
"""
function expected_metric_references(bench)::Dict{String,Float64}
    dT_ref_rmin =
        analytical_T_BE(bench.t65.r_min; n_steps=N_STEPS_TRANSIENT, dt_s=DT_S) - T_SURF
    return Dict{String,Float64}(
        "steady_center_dT_ratio_3d_2d_nx33" => 2.0 / 3.0,
        "steady_center_dT_ratio_3d_2d_nx65" => 2.0 / 3.0,
        "steady_profile_l2_rel_error_nx33" => 0.0,
        "steady_profile_l2_rel_error_nx65" => 0.0,
        "steady_energy_conservation_residual_nx65" => 0.0,
        "transient_be_vs_1dfv_rel_diff" => 0.0,
        "transient_center_dT_nx65" => dT_ref_rmin,
        "transient_profile_l2_rel_error_nx65" => 0.0,
        "transient_profile_l2_convergence_ratio_65_33" => 0.0,
        "transient_max_air_shell_temperature_nx65" => T_SURF,
    )
end

"""
    check_physics_regression_bands(metrics::Dict{String,Float64}, record_path::String) -> Tuple{Bool,Vector{String}}

Check that each steady and transient physics metric lies within its machine-independent regression
band (`PHYSICS_REGRESSION_BANDS`) around the recorded value in `record_path`.

Parameters
----------
metrics : Dict{String,Float64}
    Measured benchmark metric dictionary.
record_path : String
    Path to the validation record TOML file.

Returns
-------
Tuple{Bool,Vector{String}}
    Boolean pass flag and list of regression-band failure messages.

Raises
------
DomainError
    If `isempty(record_path)`.
"""
function check_physics_regression_bands(
    metrics::Dict{String,Float64}, record_path::String
)::Tuple{Bool,Vector{String}}
    isempty(record_path) && throw(DomainError(record_path, "Record path must be non-empty"))
    if !isfile(record_path)
        return (false, ["Validation record not found: $record_path"])
    end
    raw = TOML.parsefile(record_path)
    recorded_by_id = Dict{String,Float64}()
    for m in get(raw, "metrics", Any[])
        if m isa AbstractDict &&
            haskey(m, "id") &&
            haskey(m, "value") &&
            m["value"] isa Real
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
            continue
        end
        if !haskey(metrics, mid)
            push!(errors, "Measured metrics missing '$mid' for physics regression band")
            continue
        end
        band = PHYSICS_REGRESSION_BANDS[mid]
        measured = metrics[mid]
        recorded = recorded_by_id[mid]
        if !isapprox(measured, recorded; atol=band.atol, rtol=band.rtol)
            push!(
                errors,
                "Metric '$mid' = $measured outside physics regression band around recorded $recorded (atol=$(band.atol), rtol=$(band.rtol))",
            )
        end
    end
    return (isempty(errors), errors)
end

"""
    check_validation_reproducibility(bench, record_path::String; check_recorded_values::Bool=true) -> Tuple{Bool,Vector{String}}

Compare freshly evaluated benchmark metrics, reference values, physics regression bands, and
section-scoped configuration hashes against a committed TOML validation record.

Parameters
----------
bench : NamedTuple
    Evaluated benchmark output from `evaluate_spherical_conduction_benchmark`.
record_path : String
    Path to the TOML validation record file.
check_recorded_values : Bool, optional
    Whether to check replay tolerance (`rtol = 1e-6`, `atol` from `REPLAY_NOISE_ATOLS`) against recorded values (default `true`).

Returns
-------
Tuple{Bool,Vector{String}}
    Boolean pass flag and vector of diagnostic error messages.

Raises
------
DomainError
    If `isempty(record_path)`.
"""
function check_validation_reproducibility(
    bench, record_path::String; check_recorded_values::Bool=true
)::Tuple{Bool,Vector{String}}
    isempty(record_path) && throw(DomainError(record_path, "Record path must be non-empty"))
    if !isfile(record_path)
        return (false, ["Validation record not found: $record_path"])
    end
    raw = TOML.parsefile(record_path)
    ok_schema, errs = validate_record_dict(raw; require_doc_pages=false)
    errors = copy(errs)

    cfg65 = make_spherical_conduction_config(; Nx=65)
    ok_lit, lit_errs = check_reference_literals_match_config(cfg65)
    if !ok_lit
        append!(errors, lit_errs)
    end

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
        if !haskey(rec_cfg, k)
            push!(errors, "Recorded config missing key '$k'")
            continue
        end
        v_rec = rec_cfg[k]
        matches = if v_rec isa AbstractFloat || v_exp isa AbstractFloat
            (v_rec isa Real && v_exp isa Real) &&
                isapprox(Float64(v_rec), Float64(v_exp); rtol=1.0e-12, atol=1.0e-12)
        else
            v_rec == v_exp
        end
        if !matches
            push!(errors, "Recorded config.$k ($v_rec) does not match expected ($v_exp)")
        end
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

    exp_refs = expected_metric_references(bench)
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
        if !tol_ok
            push!(
                errors,
                "Measured metric '$mid' = $measured fails tolerance ($comp, ref=$ref, atol=$atol, rtol=$rtol)",
            )
        end
        if check_recorded_values
            replay_atol = get(REPLAY_NOISE_ATOLS, mid, 1.0e-12)
            if !isapprox(measured, recorded; rtol=REPLAY_RTOL, atol=replay_atol)
                push!(
                    errors,
                    "Measured metric '$mid' = $measured differs from recorded value $recorded (rtol=$REPLAY_RTOL, atol=$replay_atol)",
                )
            end
        end
    end

    ok_bands, band_errs = check_physics_regression_bands(bench.metrics, record_path)
    if !ok_bands
        append!(errors, band_errs)
    end

    return (ok_schema && isempty(errors), errors)
end

end # module SphericalConductionValidation

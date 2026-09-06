# Benchmark generator for Erebus.jl soft turbulence model and planetesimal cooling
using JSON
using LinearAlgebra
using StaticArrays

push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))
using Erebus
using Erebus.Physics
using Erebus.Config

println("=== Erebus.jl Soft Turbulence Benchmark Generation ===")

# -----------------------------------------------------------------------------
# Part 1: Analytical and Regularization Comparison Curves
# -----------------------------------------------------------------------------
println("Generating analytical conductivity curves...")

Fm_vals = range(0.0, 1.0; length=500)
eta_num = 1.0e12   # Pa*s
k_cond = 3.0       # W/(m*K)
dT_val = 50.0      # K
T_surf = 300.0     # K
T_val = T_surf + dT_val

eta_fluids = [10.0, 100.0, 1000.0]

function k_i2elvis(Fm, eta_num, eta_fl, k0, dT)
    if Fm >= 0.40 && dT > 10.0
        return k0 * sqrt(eta_num / eta_fl)
    else
        return k0
    end
end

k_i2elvis_curves = Dict{String,Vector{Float64}}()
k_erebus_curves = Dict{String,Vector{Float64}}()
dk_dFm_erebus = Dict{String,Vector{Float64}}()

for eta_fl in eta_fluids
    k_i2 = [k_i2elvis(f, eta_num, eta_fl, k_cond, dT_val) for f in Fm_vals]
    k_eb = [
        regularized_soft_turbulence_conductivity(
            k_cond,
            eta_num,
            eta_fl,
            f,
            T_val,
            T_surf;
            F_start=0.30,
            F_end=0.50,
            dT_min=10.0,
            k_cutoff=1.0e6,
            k_floor=1.0e-3,
        ) for f in Fm_vals
    ]

    log_k = log10.(k_eb)
    dlogk = zeros(length(Fm_vals))
    dF = Fm_vals[2] - Fm_vals[1]
    for i in 2:(length(Fm_vals) - 1)
        dlogk[i] = (log_k[i + 1] - log_k[i - 1]) / (2.0 * dF)
    end
    dlogk[1] = (log_k[2] - log_k[1]) / dF
    dlogk[end] = (log_k[end] - log_k[end - 1]) / dF

    key = string(Int(eta_fl))
    k_i2elvis_curves[key] = k_i2
    k_erebus_curves[key] = k_eb
    dk_dFm_erebus[key] = dlogk
end

dT_vals = range(0.0, 30.0; length=200)
w_T_vals = [clamp(dt / 10.0, 0.0, 1.0) for dt in dT_vals]

# -----------------------------------------------------------------------------
# Part 2: 1D Implicit Planetesimal Magma Ocean Solidification Benchmark
# -----------------------------------------------------------------------------
println("Running planetesimal cooling simulations (Soft Turbulence ON vs OFF)...")

R_planet = 50_000.0   # 50 km radius
r_core = 30_000.0     # 30 km initial molten core radius
T_core_init = 1850.0  # K (above liquidus 1800 K, Fm = 1.0)
T_surf_init = 300.0   # K
T_solidus = 1400.0    # K
T_liquidus = 1800.0   # K
L_m = 4.0e5           # J/kg
rho_s = 3300.0        # kg/m^3
cp_s = 1000.0         # J/(kg*K)
rho_cp_base = rho_s * cp_s # 3.3e6 J/(m^3*K)

function run_planetesimal_cooling(
    Nr; soft_turb::Bool=true, t_total_yr=50_000.0, dt_yr=50.0, movie_snapshots::Bool=false
)
    dr = R_planet / Nr
    r = [(i - 0.5) * dr for i in 1:Nr] # cell centers
    r_face = [i * dr for i in 0:Nr]    # cell faces

    T = zeros(Nr)
    for i in 1:Nr
        if r[i] <= r_core
            T[i] = T_core_init
        else
            xi = (r[i] - r_core) / (R_planet - r_core)
            T[i] = T_core_init - xi * (T_core_init - T_surf_init)
        end
    end

    dt_sec = dt_yr * 365.25 * 86400.0
    n_steps = Int(round(t_total_yr / dt_yr))

    times_yr = Float64[]
    T_core_hist = Float64[]
    melt_radius_hist = Float64[]
    q_surf_hist = Float64[]

    snapshot_times = if movie_snapshots
        collect(0.0:500.0:t_total_yr)
    else
        [0.0, 5000.0, 15000.0, 30000.0, 50000.0]
    end

    snapshot_T = Dict{String,Vector{Float64}}()
    snapshot_Fm = Dict{String,Vector{Float64}}()
    snapshot_k = Dict{String,Vector{Float64}}()

    push!(times_yr, 0.0)
    push!(T_core_hist, T[1])

    Fm = [clamp((T[i] - T_solidus) / (T_liquidus - T_solidus), 0.0, 1.0) for i in 1:Nr]
    r_melt = 0.0
    for i in 1:Nr
        if Fm[i] >= 0.40
            r_melt = r[i]
        end
    end
    push!(melt_radius_hist, r_melt)
    push!(q_surf_hist, k_cond * (T[Nr] - T_surf_init) / (dr / 2.0))

    snapshot_T["0.0"] = copy(T)
    snapshot_Fm["0.0"] = copy(Fm)

    k_init = [
        if soft_turb
            regularized_soft_turbulence_conductivity(
                k_cond,
                eta_num,
                100.0,
                Fm[i],
                T[i],
                T_surf_init;
                F_start=0.30,
                F_end=0.50,
                dT_min=10.0,
            )
        else
            k_cond
        end for i in 1:Nr
    ]
    snapshot_k["0.0"] = copy(k_init)

    snap_idx = 2
    k_cell = zeros(Nr)
    k_face = zeros(Nr + 1)

    # Tridiagonal vectors: A = Tridiagonal(dl, d, du)
    dl = zeros(Nr - 1)
    d = zeros(Nr)
    du = zeros(Nr - 1)
    rhs = zeros(Nr)

    for step in 1:n_steps
        current_time_yr = step * dt_yr

        # 1. Update melt fraction and cell conductivities
        for i in 1:Nr
            Fm[i] = clamp((T[i] - T_solidus) / (T_liquidus - T_solidus), 0.0, 1.0)
            if soft_turb
                k_cell[i] = regularized_soft_turbulence_conductivity(
                    k_cond,
                    eta_num,
                    100.0,
                    Fm[i],
                    T[i],
                    T_surf_init;
                    F_start=0.30,
                    F_end=0.50,
                    dT_min=10.0,
                    k_cutoff=1.0e6,
                    k_floor=1.0e-3,
                )
            else
                k_cell[i] = k_cond
            end
        end

        # 2. Harmonic mean on internal faces
        k_face[1] = k_cell[1]
        for i in 2:Nr
            k_face[i] =
                2.0 * k_cell[i - 1] * k_cell[i] / (k_cell[i - 1] + k_cell[i] + 1e-12)
        end
        k_face[Nr + 1] = k_cond # surface boundary

        # 3. Assemble implicit finite volume linear system:
        # C_i * (T_i^{n+1} - T_i^n) = D_{i-1/2}*(T_{i-1}^{n+1} - T_i^{n+1}) + D_{i+1/2}*(T_{i+1}^{n+1} - T_i^{n+1})
        for i in 1:Nr
            # Apparent heat capacity buffering in melting mush
            rho_cp_eff = rho_cp_base
            if T[i] >= T_solidus && T[i] <= T_liquidus
                rho_cp_eff += rho_s * L_m / (T_liquidus - T_solidus)
            end

            C_i = r[i]^2 * dr * rho_cp_eff / dt_sec
            D_left = (i == 1) ? 0.0 : r_face[i]^2 * k_face[i] / dr
            D_right = if (i == Nr)
                r_face[i + 1]^2 * k_face[i + 1] / (dr / 2.0)
            else
                r_face[i + 1]^2 * k_face[i + 1] / dr
            end

            d[i] = C_i + D_left + D_right
            rhs[i] = C_i * T[i]

            if i == Nr
                rhs[i] += D_right * T_surf_init
            end

            if i > 1
                dl[i - 1] = -D_left
            end
            if i < Nr
                du[i] = -D_right
            end
        end

        # 4. Direct tridiagonal solve O(N)
        A = Tridiagonal(dl, d, du)
        T .= A \ rhs

        # 5. Track diagnostics
        push!(times_yr, current_time_yr)
        push!(T_core_hist, T[1])

        r_melt = 0.0
        for i in 1:Nr
            if Fm[i] >= 0.40
                r_melt = r[i]
            end
        end
        push!(melt_radius_hist, r_melt)
        q_surf = k_face[Nr + 1] * (T[Nr] - T_surf_init) / (dr / 2.0)
        push!(q_surf_hist, q_surf)

        if snap_idx <= length(snapshot_times) && current_time_yr >= snapshot_times[snap_idx] - 1e-6
            t_snap = snapshot_times[snap_idx]
            key = string(round(t_snap; digits=1))
            snapshot_T[key] = copy(T)
            snapshot_Fm[key] = copy(Fm)
            snapshot_k[key] = copy(k_cell)
            snap_idx += 1
        end
    end

    return Dict(
        "r" => collect(r),
        "times_yr" => times_yr,
        "T_core_hist" => T_core_hist,
        "melt_radius_hist" => melt_radius_hist,
        "q_surf_hist" => q_surf_hist,
        "snapshot_T" => snapshot_T,
        "snapshot_Fm" => snapshot_Fm,
        "snapshot_k" => snapshot_k,
    )
end

println("Running Case A: Soft Turbulence OFF (Nr=128)...")
res_off = run_planetesimal_cooling(128; soft_turb=false, t_total_yr=50_000.0, dt_yr=50.0)

println("Running Case B: Soft Turbulence ON (Nr=128, movie frames)...")
res_on = run_planetesimal_cooling(
    128; soft_turb=true, t_total_yr=50_000.0, dt_yr=50.0, movie_snapshots=true
)

println("Running Resolution Convergence Tests: Nr in [32, 64, 128, 256]...")
res_grid32 = run_planetesimal_cooling(32; soft_turb=true, t_total_yr=50_000.0, dt_yr=50.0)
res_grid64 = run_planetesimal_cooling(64; soft_turb=true, t_total_yr=50_000.0, dt_yr=50.0)
res_grid128 = res_on
res_grid256 = run_planetesimal_cooling(256; soft_turb=true, t_total_yr=50_000.0, dt_yr=50.0)

output_dir = normpath(joinpath(@__DIR__, "..", "output_files"))
mkpath(output_dir)
output_path = joinpath(output_dir, "soft_turbulence_benchmark_data.json")

all_data = Dict(
    "Fm_vals" => collect(Fm_vals),
    "k_i2elvis_curves" => k_i2elvis_curves,
    "k_erebus_curves" => k_erebus_curves,
    "dk_dFm_erebus" => dk_dFm_erebus,
    "dT_vals" => collect(dT_vals),
    "w_T_vals" => w_T_vals,
    "res_off" => res_off,
    "res_on" => res_on,
    "res_grid32" => res_grid32,
    "res_grid64" => res_grid64,
    "res_grid128" => res_grid128,
    "res_grid256" => res_grid256,
)

open(output_path, "w") do io
    return JSON.print(io, all_data)
end

println("Benchmark data saved successfully to $output_path")

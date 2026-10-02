# Benchmark exporter for Erebus.jl radiogenic isotope decay power
using JSON
using LinearAlgebra

push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))
using Erebus

println("=== Erebus.jl Radiogenic Decay Benchmark Export ===")

# Literature constants and anchors
# 26Al: Russell et al. (1996), Lichtenberg et al. (2019)
const T_HALF_AL_YR = 717_000.0
const T_HALF_AL_S = T_HALF_AL_YR * 31_540_000.0
const TAU_AL_S = T_HALF_AL_S / log(2.0)
const E_AL_J = 5.0470e-13
const F_AL = 1.9e23               # 27Al atoms per kg bulk rock
const RATIO_AL = 5.0e-5           # Canonical initial 26Al/27Al

# 60Fe: Tang & Dauphas (2012)
const T_HALF_FE_YR = 2_620_000.0
const T_HALF_FE_S = T_HALF_FE_YR * 31_540_000.0
const TAU_FE_S = T_HALF_FE_S / log(2.0)
const E_FE_J = 4.34e-13
const F_FE = 1.957e24             # 56Fe atoms per kg bulk rock (18.2 wt% Fe in CI)
const RATIO_FE = 1.15e-8          # Canonical initial 60Fe/56Fe

# Simulation time domain: 0 to 10 Myr
const N_POINTS = 101
t_Myr = collect(range(0.0, 10.0; length=N_POINTS))
t_sec = t_Myr .* (1.0e6 * 31_540_000.0)

Q_al = [Erebus.Q_radiogenic(F_AL, RATIO_AL, E_AL_J, TAU_AL_S, ts) for ts in t_sec]
Q_fe = [Erebus.Q_radiogenic(F_FE, RATIO_FE, E_FE_J, TAU_FE_S, ts) for ts in t_sec]

# Linear regression in log space to fit decay mean lifetime tau:
# ln(Q(t)) = ln(Q0) - t / tau
# Design matrix A = [ones(N) -t]
function fit_half_life(times_s, powers)
    N = length(times_s)
    log_powers = log.(powers)
    A = [ones(N) times_s]
    coeffs = A \ log_powers  # coeffs[1] = ln(Q0), coeffs[2] = -1/tau
    tau_fit_s = -1.0 / coeffs[2]
    t_half_fit_yr = (tau_fit_s * log(2.0)) / 31_540_000.0
    return t_half_fit_yr, tau_fit_s
end

t_half_al_fit_yr, tau_al_fit_s = fit_half_life(t_sec, Q_al)
t_half_fe_fit_yr, tau_fe_fit_s = fit_half_life(t_sec, Q_fe)

err_al = abs(t_half_al_fit_yr - T_HALF_AL_YR) / T_HALF_AL_YR
err_fe = abs(t_half_fe_fit_yr - T_HALF_FE_YR) / T_HALF_FE_YR

println(
    "26Al half-life literature: $(T_HALF_AL_YR / 1e6) Myr, fitted: $(t_half_al_fit_yr / 1e6) Myr, error: $(err_al)",
)
println(
    "60Fe half-life literature: $(T_HALF_FE_YR / 1e6) Myr, fitted: $(t_half_fe_fit_yr / 1e6) Myr, error: $(err_fe)",
)

pass_al = err_al < 0.01
pass_fe = err_fe < 0.01
all_passed = pass_al && pass_fe

if !all_passed
    error("Radiogenic decay benchmark failed: fitted half-life error exceeds 1%")
end

output_dir = normpath(joinpath(@__DIR__, "..", "output_files"))
mkpath(output_dir)
output_path = joinpath(output_dir, "radiogenic_decay_benchmark_data.json")

data = Dict(
    "t_Myr" => t_Myr,
    "Q_al_W_kg" => Q_al,
    "Q_fe_W_kg" => Q_fe,
    "Q0_al" => Q_al[1],
    "Q0_fe" => Q_fe[1],
    "t_half_al_Myr" => T_HALF_AL_YR / 1.0e6,
    "t_half_fe_Myr" => T_HALF_FE_YR / 1.0e6,
    "fitted_t_half_al_Myr" => t_half_al_fit_yr / 1.0e6,
    "fitted_t_half_fe_Myr" => t_half_fe_fit_yr / 1.0e6,
    "err_al" => err_al,
    "err_fe" => err_fe,
    "ratio_al" => RATIO_AL,
    "ratio_fe" => RATIO_FE,
    "passed" => all_passed,
)

open(output_path, "w") do f
    return JSON.print(f, data, 2)
end

println("Exported benchmark data to: $(output_path)")

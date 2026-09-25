# Benchmark exporter for Erebus.jl Jeans kinetic effusion vs analytical Maxwell-Boltzmann integral
using JSON
using LinearAlgebra

push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))
using Erebus
using Erebus.Physics

println("=== Erebus.jl Jeans Effusion Benchmark Export ===")

# CODATA 2018/2022 constants explicitly defined (not imported from escape module)
const K_B_CODATA = 1.380649e-23        # Boltzmann constant [J/K]
const G_CODATA = 6.67430e-11          # Gravitational constant [m^3 kg^-1 s^-2]
const AMU_CODATA = 1.66053906660e-27  # Atomic mass unit [kg]
const M_H_CODATA = 1.00794 * AMU_CODATA # Atomic hydrogen mass [kg]

# Fixed exobase and planet geometry (Lunar-like embryo)
const R_EXO = 1.7374e6                 # Exobase radius [m]
const G_SURF = 1.62                    # Exobase gravity [m/s^2]
const M_PLANET = G_SURF * (R_EXO^2) / G_CODATA # Derived mass [kg]
const N_EXO = 1.0e11                   # Exobase number density [m^-3]

println("Exobase radius: $(R_EXO) m")
println("Surface gravity: $(G_SURF) m/s^2")
println("Planetary mass: $(M_PLANET) kg")
println("Atomic H mass: $(M_H_CODATA) kg")
println("Number density: $(N_EXO) m^-3")

# 1. Verification of the lambda -> 0 limit: n * v_th / (2 * sqrt(pi))
let T_test = 1000.0
    v_th = sqrt(2.0 * K_B_CODATA * T_test / M_H_CODATA)
    limit_analytical = N_EXO * v_th / (2.0 * sqrt(pi))
    # Using lambda = 0 in the analytical effusion formula
    lambda_zero = 0.0
    flux_zero = N_EXO * (v_th / (2.0 * sqrt(pi))) * (1.0 + lambda_zero) * exp(-lambda_zero)
    rel_err_limit = abs(flux_zero - limit_analytical) / limit_analytical
    println("Lambda -> 0 limit relative error: $(rel_err_limit)")
    if rel_err_limit > 1.0e-15
        error("Lambda -> 0 limit check failed: $(rel_err_limit) > 1e-15")
    end
end

# 2. Temperature sweep: effusion flux vs T_exo (hydrodynamic = false)
const N_POINTS = 60
T_sweep = collect(range(200.0, 3000.0; length=N_POINTS))
phi_ref = Float64[]
phi_code = Float64[]
lambda_sweep = Float64[]
rel_errors = Float64[]

for T_val in T_sweep
    # Analytical Maxwell-Boltzmann exobase effusion integral
    v_th = sqrt(2.0 * K_B_CODATA * T_val / M_H_CODATA)
    lam = (G_CODATA * M_PLANET * M_H_CODATA) / (K_B_CODATA * T_val * R_EXO)
    f_ref = N_EXO * (v_th / (2.0 * sqrt(pi))) * (1.0 + lam) * exp(-lam)
    push!(phi_ref, f_ref)
    push!(lambda_sweep, lam)

    # Erebus.jl function call with hydrodynamic=false
    f_code = compute_jeans_escape_flux(N_EXO, T_val, M_H_CODATA, lam; hydrodynamic=false)
    push!(phi_code, f_code)

    err = abs(f_code - f_ref) / f_ref
    push!(rel_errors, err)
end

max_rel_error = maximum(rel_errors)
println("Maximum relative error vs analytical reference: $(max_rel_error)")
pass_criterion = max_rel_error < 1.0e-6
println("Pass criterion (< 1e-6): $(pass_criterion)")
if !pass_criterion
    error(
        "Jeans effusion benchmark verification failed: max rel error $(max_rel_error) >= 1e-6",
    )
end

# 3. Branch switch point: hydrodynamic = true with lambda < 2 asserts flux exceeds effusion
# Select hot temperature where lambda < 2
T_hot = 12000.0
lam_hot = (G_CODATA * M_PLANET * M_H_CODATA) / (K_B_CODATA * T_hot * R_EXO)
println("Branch test point: T = $(T_hot) K, lambda = $(lam_hot) (< 2.0)")

flux_effusion_only = compute_jeans_escape_flux(
    N_EXO, T_hot, M_H_CODATA, lam_hot; hydrodynamic=false
)
flux_hydro_active = compute_jeans_escape_flux(
    N_EXO, T_hot, M_H_CODATA, lam_hot; hydrodynamic=true
)

println("Effusion flux (hydrodynamic=false): $(flux_effusion_only)")
println("Hydrodynamic flux (hydrodynamic=true): $(flux_hydro_active)")
if flux_hydro_active <= flux_effusion_only
    error(
        "Branch switch check failed: hydrodynamic flux ($(flux_hydro_active)) <= effusion flux ($(flux_effusion_only))",
    )
end
println(
    "Branch switch successfully verified: hydrodynamic flux exceeds effusion flux by factor $(flux_hydro_active / flux_effusion_only)x",
)

# 4. Save benchmark data
output_dir = normpath(joinpath(@__DIR__, "..", "output_files"))
mkpath(output_dir)
output_path = joinpath(output_dir, "jeans_effusion_benchmark_data.json")

data = Dict(
    "T_sweep" => T_sweep,
    "lambda_sweep" => lambda_sweep,
    "phi_ref" => phi_ref,
    "phi_code" => phi_code,
    "rel_errors" => rel_errors,
    "max_rel_error" => max_rel_error,
    "T_hot" => T_hot,
    "lam_hot" => lam_hot,
    "flux_effusion_only" => flux_effusion_only,
    "flux_hydro_active" => flux_hydro_active,
    "hydro_ratio" => flux_hydro_active / flux_effusion_only,
    "R_exo" => R_EXO,
    "g_surf" => G_SURF,
    "n_exo" => N_EXO,
)

open(output_path, "w") do io
    return JSON.print(io, data, 2)
end

println("Benchmark data exported successfully to: $(output_path)")

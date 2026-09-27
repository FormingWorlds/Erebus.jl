# Benchmark data exporter for hydrofracture regularisation ramp and Darcy under-relaxation
using JSON

push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))
using Erebus

println("=== Erebus.jl Hydrofracture Ramp Benchmark Export ===")

# 1. Sweep normalized overpressure x = (-P_eff - sigma_t) / sigma_t
x_sweep = collect(range(-0.1, 0.3; length=401))

# Compare kink law (delta=0) with C^1 ramp laws (delta > 0)
delta_values = [0.0, 0.02, 0.05, 0.10]
ramp_curves = Dict{String,Any}()

for delta in delta_values
    s_vals = [hydrofracture_overpressure_ramp(x, delta) for x in x_sweep]
    ds_vals = Float64[]
    for x in x_sweep
        if delta == 0.0
            push!(ds_vals, x > 0.0 ? 1.0 : 0.0)
        else
            if x <= 0.0
                push!(ds_vals, 0.0)
            elseif x < delta
                push!(ds_vals, x / delta)
            else
                push!(ds_vals, 1.0)
            end
        end
    end
    f_vals = [1.0 + 1000.0 * s for s in s_vals]
    ramp_curves[string(delta)] = Dict("s" => s_vals, "ds" => ds_vals, "f" => f_vals)
end

# 2. Darcy under-relaxation convergence histories
# r^(k) - r_new = (1 - theta)^k (r^(0) - r_new)
k_iters = collect(0:20)
theta_values = [0.1, 0.3, 0.5, 1.0]
relaxation_curves = Dict{String,Any}()

for theta in theta_values
    res_vals = [(1.0 - theta)^k for k in k_iters]
    relaxation_curves[string(theta)] = res_vals
end

# Assemble benchmark payload
payload = Dict(
    "x_sweep" => x_sweep,
    "delta_values" => delta_values,
    "ramp_curves" => ramp_curves,
    "k_iters" => k_iters,
    "theta_values" => theta_values,
    "relaxation_curves" => relaxation_curves,
)

out_dir = normpath(joinpath(@__DIR__, "..", "output_files"))
mkpath(out_dir)
out_file = joinpath(out_dir, "hydrofracture_ramp_benchmark_data.json")

open(out_file, "w") do io
    return JSON.print(io, payload, 2)
end

println("Successfully exported benchmark data to $out_file")

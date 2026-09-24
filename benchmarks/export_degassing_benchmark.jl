# Benchmark exporter for Erebus.jl magma ocean volatile degassing
using JSON
using LinearAlgebra

push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))
using Erebus
using Erebus.Physics
using Erebus.Config

println("=== Erebus.jl Magma Ocean Degassing Benchmark Export ===")

const R_PLANET = 50_000.0
const P_SURF = 1.0e5       # 1 bar surface atmospheric pressure [Pa]
const T_MELT_REF = 1500.0   # Reference melt temperature [K]
const DELTA_IW = 0.0        # Neutral oxygen fugacity
const WATER_AS = 0.40       # Burnham/Dixon water solubility coefficient [wt% / MPa^0.5]

# Thermodynamic gas speciation at surface ambient conditions
spec = solve_chnos_speciation(P_SURF, T_MELT_REF, DELTA_IW)
p_H2O_Pa = spec.p_H2O_Pa
p_H2O_MPa = p_H2O_Pa * 1.0e-6

println("Surface pressure: $(P_SURF) Pa")
println("Water partial pressure: $(p_H2O_Pa) Pa ($(p_H2O_MPa) MPa)")
println("Water solubility coefficient As: $(WATER_AS) wt%/MPa^0.5")

# 1. Analytical reference curve from Burnham (1979) / Dixon et al. (1995) law:
# w_sat(F) = F * As * sqrt(p_H2O_MPa)
const N_REF = 100
F_ref = collect(range(0.1, 1.0; length=N_REF))
retained_ref = [F * WATER_AS * sqrt(p_H2O_MPa) for F in F_ref]

# 2. Simulation points evaluated via degas_magma_ocean_markers!
cfg = MagmaOceanDegassingConfig(;
    active=true,
    mode=:dynamic_flux,
    degas_depth_fraction=0.90,
    F_melt_threshold=0.40,
    water_As=WATER_AS,
    efficiency=1.0,
)

const N_SIM = 10
F_sim = collect(range(0.1, 1.0; length=N_SIM))
retained_sim = Float64[]
rel_errors = Float64[]

for F_val in F_sim
    xm = [0.0]
    ym = [0.95 * R_PLANET]
    tm = [2]
    tkm = [T_MELT_REF]
    Fm = [F_val]
    Fm_old = [F_val]
    XH2Om = [2.0] # 2 wt% (supersaturated)
    XCm = [0.0]
    XNm = [0.0]
    XSm = [0.0]

    degas_magma_ocean_markers!(
        xm, ym, tm, tkm, Fm, Fm_old, XH2Om, XCm, XNm, XSm,
        1, 1000.0, P_SURF, R_PLANET, cfg, T_MELT_REF;
        marker_volume=1.0,
    )

    w_ret = XH2Om[1]
    push!(retained_sim, w_ret)
    w_expected = F_val * WATER_AS * sqrt(p_H2O_MPa)
    err = abs(w_ret - w_expected) / w_expected
    push!(rel_errors, err)
end

max_rel_error = maximum(rel_errors)
println("Maximum relative error vs analytical law: $(max_rel_error)")
pass_criterion = max_rel_error < 1.0e-6
println("Pass criterion (< 1e-6): $(pass_criterion)")
if !pass_criterion
    error("Degassing benchmark verification failed: max rel error $(max_rel_error) >= 1e-6")
end

output_dir = normpath(joinpath(@__DIR__, "..", "output_files"))
mkpath(output_dir)
output_path = joinpath(output_dir, "degassing_benchmark_data.json")

data = Dict(
    "F_ref" => F_ref,
    "retained_ref_wtpct" => retained_ref,
    "F_sim" => F_sim,
    "retained_sim_wtpct" => retained_sim,
    "rel_errors" => rel_errors,
    "max_rel_error" => max_rel_error,
    "p_H2O_MPa" => p_H2O_MPa,
    "water_As" => WATER_AS,
    "T_melt_ref" => T_MELT_REF,
    "passed" => pass_criterion,
)

open(output_path, "w") do f
    JSON.print(f, data, 2)
end
println("Exported benchmark data to: $(output_path)")

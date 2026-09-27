# Benchmark exporter for Erebus.jl 2D analytical thermal slab diffusion
using JSON
using LinearAlgebra

push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))
using Erebus

println("=== Erebus.jl 2D Thermal Slab Conduction Benchmark Export ===")

const LX = 100_000.0   # 100 km box width [m]
const LY = 100_000.0   # 100 km box height [m]
const K_VAL = 3.0      # Thermal conductivity [W/(m K)]
const RHOCP_VAL = 3.0e6 # Volumetric heat capacity [J/(m^3 K)]
const KAPPA = K_VAL / RHOCP_VAL # Thermal diffusivity [m^2/s]
const T0 = 300.0       # Background temperature [K]
const DELTA_T = 50.0   # Temperature perturbation amplitude [K]

# Fundamental decay rate and characteristic diffusion timescale
const LAMBDA = KAPPA * π^2 * (1.0 / LX^2 + 1.0 / LY^2)
const TAU_DIFF = 1.0 / LAMBDA # ~5.066e14 s (~16.05 Ma)
const YEARLENGTH = 3.15576e7  # 1 year in seconds

function main()
    println("Domain size: $(LX / 1000.0) km x $(LY / 1000.0) km")
    println("Diffusivity kappa: $(KAPPA) m^2/s")
    println("Diffusion timescale tau: $(TAU_DIFF / YEARLENGTH / 1.0e6) Ma")

    # -----------------------------------------------------------------------------
    # 1. Multi-epoch centerline temperature profiles (N = 33)
    # -----------------------------------------------------------------------------
    N_GRID = 33
    coords_main = Erebus.GridCoordinates(N_GRID, N_GRID; xsize=LX, ysize=LY)
    Nx1, Ny1 = coords_main.Nx1, coords_main.Ny1

    # Target evaluation times: t = 0, 0.1 tau, 0.25 tau, 0.5 tau (~0, 1.6, 4.0, 8.0 Ma)
    eval_fractions = [0.0, 0.1, 0.25, 0.5]
    eval_times = eval_fractions .* TAU_DIFF
    times_Myr = eval_times ./ (YEARLENGTH * 1.0e6)

    # Setup initial temperature field
    tk_sim = zeros(Ny1, Nx1)
    for j in 1:Nx1, i in 1:Ny1
        tk_sim[i, j] =
            T0 + DELTA_T * cos(π * coords_main.xp[j] / LX) * cos(π * coords_main.yp[i] / LY)
    end

    RHOCP = fill(RHOCP_VAL, Ny1, Nx1)
    KX = fill(K_VAL, N_GRID, Nx1)
    KY = fill(K_VAL, Ny1, N_GRID)
    HR = zeros(Ny1, Nx1)
    HA = zeros(Ny1, Nx1)
    HS = zeros(Ny1, Nx1)
    DHP = zeros(Ny1, Nx1)
    RT = zeros(Ny1 * Nx1)

    # Centerline indices (closest to y = LY / 2)
    # Internal nodes run j = 2:(Nx1-1), i = 2:(Ny1-1)
    y_mid = LY / 2.0
    mid_i = argmin(abs.(coords_main.yp .- y_mid))
    x_centerline_km = coords_main.xp[2:(Nx1 - 1)] ./ 1000.0

    profiles_num = Vector{Vector{Float64}}()
    profiles_ana = Vector{Vector{Float64}}()
    max_profile_errors = Float64[]

    current_time = 0.0
    # Record t = 0
    ana_t0 = [
        T0 +
        DELTA_T * cos(π * coords_main.xp[j] / LX) * cos(π * coords_main.yp[mid_i] / LY) for
        j in 2:(Nx1 - 1)
    ]
    num_t0 = tk_sim[mid_i, 2:(Nx1 - 1)]
    push!(profiles_ana, ana_t0)
    push!(profiles_num, copy(num_t0))
    push!(max_profile_errors, maximum(abs.(num_t0 .- ana_t0)) / DELTA_T)

    # Step to subsequent times
    dt_step = 1.0e11 # ~3,170 years per step
    for t_target in eval_times[2:end]
        n_substeps = round(Int, (t_target - current_time) / dt_step)
        dt_actual = (t_target - current_time) / n_substeps
        for _ in 1:n_substeps
            LT = Erebus.assemble_thermal_lse!(
                tk_sim, RHOCP, KX, KY, HR, HA, HS, DHP, RT, dt_actual; coords=coords_main
            )
            sol = LT \ RT
            tk_sim .= reshape(sol, Ny1, Nx1)
        end
        current_time = t_target

        decay_factor = exp(-LAMBDA * current_time)
        ana_prof = [
            T0 +
            DELTA_T *
            cos(π * coords_main.xp[j] / LX) *
            cos(π * coords_main.yp[mid_i] / LY) *
            decay_factor for j in 2:(Nx1 - 1)
        ]
        num_prof = copy(tk_sim[mid_i, 2:(Nx1 - 1)])
        err = maximum(abs.(num_prof .- ana_prof)) / DELTA_T

        push!(profiles_ana, ana_prof)
        push!(profiles_num, num_prof)
        push!(max_profile_errors, err)
        println(
            "Time $(round(current_time / YEARLENGTH / 1.0e6; digits=2)) Ma: max rel error = $err",
        )
    end

    # -----------------------------------------------------------------------------
    # 2. Grid convergence study at t = 0.1 tau (~1.6 Ma)
    # -----------------------------------------------------------------------------
    resolutions = [17, 33, 65]
    dx_values_km = Float64[]
    l2_errors = Float64[]
    linf_errors = Float64[]
    energy_drifts = Float64[]

    t_test = 0.1 * TAU_DIFF
    decay_test = exp(-LAMBDA * t_test)

    for N in resolutions
        coords_res = Erebus.GridCoordinates(N, N; xsize=LX, ysize=LY)
        Nx1_res, Ny1_res = coords_res.Nx1, coords_res.Ny1

        tk_res = zeros(Ny1_res, Nx1_res)
        for j in 1:Nx1_res, i in 1:Ny1_res
            tk_res[i, j] =
                T0 +
                DELTA_T * cos(π * coords_res.xp[j] / LX) * cos(π * coords_res.yp[i] / LY)
        end

        RHOCP_res = fill(RHOCP_VAL, Ny1_res, Nx1_res)
        KX_res = fill(K_VAL, N, Nx1_res)
        KY_res = fill(K_VAL, Ny1_res, N)
        HR_res = zeros(Ny1_res, Nx1_res)
        HA_res = zeros(Ny1_res, Nx1_res)
        HS_res = zeros(Ny1_res, Nx1_res)
        DHP_res = zeros(Ny1_res, Nx1_res)
        RT_res = zeros(Ny1_res * Nx1_res)

        E_init =
            sum(tk_res[2:(Ny1_res - 1), 2:(Nx1_res - 1)]) *
            RHOCP_VAL *
            coords_res.dx *
            coords_res.dy

        # 10 steps to reach t_test
        n_steps_conv = 10
        dt_conv = t_test / n_steps_conv
        for _ in 1:n_steps_conv
            LT = Erebus.assemble_thermal_lse!(
                tk_res,
                RHOCP_res,
                KX_res,
                KY_res,
                HR_res,
                HA_res,
                HS_res,
                DHP_res,
                RT_res,
                dt_conv;
                coords=coords_res,
            )
            sol = LT \ RT_res
            tk_res .= reshape(sol, Ny1_res, Nx1_res)
        end

        # Analytical solution
        diffs = Float64[]
        for j in 2:(Nx1_res - 1), i in 2:(Ny1_res - 1)
            ana_val =
                T0 +
                DELTA_T *
                cos(π * coords_res.xp[j] / LX) *
                cos(π * coords_res.yp[i] / LY) *
                decay_test
            push!(diffs, abs(tk_res[i, j] - ana_val))
        end

        linf = maximum(diffs) / DELTA_T
        l2 = sqrt(sum(diffs .^ 2) / length(diffs)) / DELTA_T
        E_final =
            sum(tk_res[2:(Ny1_res - 1), 2:(Nx1_res - 1)]) *
            RHOCP_VAL *
            coords_res.dx *
            coords_res.dy
        drift = abs(E_final - E_init) / E_init

        push!(dx_values_km, coords_res.dx / 1000.0)
        push!(linf_errors, linf)
        push!(l2_errors, l2)
        push!(energy_drifts, drift)
        println(
            "Grid N=$N (dx=$(coords_res.dx/1000.0) km): Linf=$linf, L2=$l2, energy drift=$drift",
        )
    end

    # Verification criteria
    pass_criterion =
        maximum(max_profile_errors) < 1.0e-3 && maximum(energy_drifts) < 1.0e-12
    println("Pass criterion (< 1e-3 error, < 1e-12 energy drift): $pass_criterion")
    if !pass_criterion
        error("Thermal slab benchmark verification failed")
    end

    # -----------------------------------------------------------------------------
    # 3. Export to JSON
    # -----------------------------------------------------------------------------
    output_dir = normpath(joinpath(@__DIR__, "..", "output_files"))
    mkpath(output_dir)
    output_path = joinpath(output_dir, "thermal_slab_benchmark_data.json")

    data = Dict(
        "times_Myr" => times_Myr,
        "x_centerline_km" => x_centerline_km,
        "profiles_analytical" => profiles_ana,
        "profiles_numerical" => profiles_num,
        "max_profile_errors" => max_profile_errors,
        "resolutions" => resolutions,
        "dx_km" => dx_values_km,
        "l2_errors" => l2_errors,
        "linf_errors" => linf_errors,
        "energy_drifts" => energy_drifts,
        "LX_km" => LX / 1000.0,
        "LY_km" => LY / 1000.0,
        "T0_K" => T0,
        "delta_T_K" => DELTA_T,
        "tau_diff_Myr" => TAU_DIFF / (YEARLENGTH * 1.0e6),
        "passed" => pass_criterion,
    )

    open(output_path, "w") do f
        return JSON.print(f, data, 2)
    end
    return println("Exported benchmark data to: $(output_path)")
end

main()

# Benchmark exporter for Erebus.jl gravity for differentiated bodies
using JSON
using LinearAlgebra

push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))
using Erebus

println("=== Erebus.jl Two-Layer Differentiated Body Gravity Benchmark Export ===")

const R_PLANET = 50_000.0   # Planet radius [m]
const R_CORE = 25_000.0     # Core radius [m]
const RHO_CORE = 7000.0     # Core density [kg/m^3]
const RHO_MANTLE = 3000.0   # Mantle density [kg/m^3]
const RHO_AIR = 1.0         # Sticky-air density [kg/m^3]
const CONST_G = 6.6743e-11  # Gravitational constant [m^3/(kg s^2)]
const DOMAIN_SIZE = 140_000.0 # Computational domain size [m]

function analytical_gravity_3d(r::Float64)
    if r <= 0.0
        return 0.0
    elseif r <= R_CORE
        return (4.0 / 3.0) * π * CONST_G * RHO_CORE * r
    elseif r <= R_PLANET
        excess = (RHO_CORE - RHO_MANTLE) * (R_CORE^3) / (r^2)
        return (4.0 / 3.0) * π * CONST_G * (RHO_MANTLE * r + excess)
    else
        m_tot =
            (4.0 / 3.0) *
            π *
            (RHO_MANTLE * (R_PLANET^3) + (RHO_CORE - RHO_MANTLE) * (R_CORE^3))
        return CONST_G * m_tot / (r^2)
    end
end

function main()
    Nx = 65
    Ny = 65
    coords = Erebus.GridCoordinates(Nx, Ny; xsize=DOMAIN_SIZE, ysize=DOMAIN_SIZE)
    xc = coords.xcenter
    yc = coords.ycenter

    Nxm = (Nx - 1) * 4
    Nym = (Ny - 1) * 4
    dxm = DOMAIN_SIZE / Nxm
    dym = DOMAIN_SIZE / Nym
    Am = dxm * dym

    xm = Float64[]
    ym = Float64[]
    rhototalm = Float64[]
    tm = Int[]

    for j in 1:Nxm, i in 1:Nym
        x = (j - 0.5) * dxm
        y = (i - 0.5) * dym
        r = hypot(x - xc, y - yc)
        push!(xm, x)
        push!(ym, y)
        if r <= R_CORE
            push!(rhototalm, RHO_CORE)
            push!(tm, 1)
        elseif r <= R_PLANET
            push!(rhototalm, RHO_MANTLE)
            push!(tm, 2)
        else
            push!(rhototalm, RHO_AIR)
            push!(tm, 3)
        end
    end

    # 1. Compute enclosed-mass formulation
    gx_enc = zeros(Float64, coords.Ny1, coords.Nx1)
    gy_enc = zeros(Float64, coords.Ny1, coords.Nx1)
    r_bins, g_bins = Erebus.compute_gravity_enclosed_mass!(
        gx_enc,
        gy_enc;
        xm=xm,
        ym=ym,
        rhototalm=rhototalm,
        tm=tm,
        coords=coords,
        gravity_nr_factor=4,
        rplanet=R_PLANET,
    )

    # 2. Compute 2D Poisson solution
    RHO_p = zeros(Float64, coords.Ny1, coords.Nx1)
    for j in 1:coords.Nx1, i in 1:coords.Ny1
        r = hypot(coords.xp[j] - xc, coords.yp[i] - yc)
        if r <= R_CORE
            RHO_p[i, j] = RHO_CORE
        elseif r <= R_PLANET
            RHO_p[i, j] = RHO_MANTLE
        end
    end

    SP = zeros(Float64, coords.Nx1 * coords.Ny1)
    RP = zeros(Float64, coords.Nx1 * coords.Ny1)
    FI = zeros(Float64, coords.Ny1, coords.Nx1)
    gx_p2d = zeros(Float64, coords.Ny1, coords.Nx1)
    gy_p2d = zeros(Float64, coords.Ny1, coords.Nx1)
    Erebus.compute_gravity_solution!(SP, RP, RHO_p, FI, gx_p2d, gy_p2d; coords=coords)

    # 3. Sample profiles along horizontal centerline from xc to xc + 68 km
    r_eval_km = collect(range(0.5, 68.0; length=140))
    r_eval_m = r_eval_km .* 1000.0

    g_ana_vals = [analytical_gravity_3d(r) for r in r_eval_m]
    i_mid = argmin(abs.(coords.yvx .- yc))

    g_p2d_vals = Float64[]
    g_enc_grid_vals = Float64[]

    for r_target in r_eval_m
        j_left = findlast(coords.xvx .<= (xc + r_target))
        if j_left === nothing || j_left >= coords.Nx1
            push!(g_p2d_vals, 0.0)
            push!(g_enc_grid_vals, 0.0)
            continue
        end
        j_right = j_left + 1
        xl = coords.xvx[j_left]
        xr = coords.xvx[j_right]
        frac = (xc + r_target - xl) / (xr - xl)

        g_p2d = abs(gx_p2d[i_mid, j_left] * (1.0 - frac) + gx_p2d[i_mid, j_right] * frac)
        push!(g_p2d_vals, g_p2d)

        g_enc = abs(gx_enc[i_mid, j_left] * (1.0 - frac) + gx_enc[i_mid, j_right] * frac)
        push!(g_enc_grid_vals, g_enc)
    end

    # Verification checks
    g_surf_ana = analytical_gravity_3d(R_PLANET)
    g_surf_enc = g_bins[end]
    err_surf = abs(g_surf_enc - g_surf_ana) / g_surf_ana
    println("Surface gravity: analytical = $(g_surf_ana), enclosed_mass = $(g_surf_enc)")
    println("Surface relative error: $(err_surf)")
    if err_surf > 0.01
        error("Enclosed mass surface gravity error exceeds 1% tolerance: $(err_surf)")
    end

    # 4. Export to JSON
    output_dir = normpath(joinpath(@__DIR__, "..", "output_files"))
    mkpath(output_dir)
    output_path = joinpath(output_dir, "gravity_two_layer_benchmark_data.json")

    data = Dict(
        "r_eval_km" => r_eval_km,
        "g_analytical" => g_ana_vals,
        "g_poisson2d" => g_p2d_vals,
        "g_enclosed_grid" => g_enc_grid_vals,
        "r_bins_km" => r_bins ./ 1000.0,
        "g_bins" => g_bins,
        "R_planet_km" => R_PLANET / 1000.0,
        "R_core_km" => R_CORE / 1000.0,
        "rho_core" => RHO_CORE,
        "rho_mantle" => RHO_MANTLE,
        "g_surface_analytical" => g_surf_ana,
        "g_surface_enclosed" => g_surf_enc,
        "g_surface_error" => err_surf,
        "passed" => (err_surf <= 0.01),
    )

    open(output_path, "w") do f
        return JSON.print(f, data, 2)
    end
    println("Exported benchmark data to: $(output_path)")
    return nothing
end

main()

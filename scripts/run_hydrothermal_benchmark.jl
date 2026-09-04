#!/usr/bin/env julia
"""
2D Hydrothermal Circulation Planetesimal Benchmark Runner.

Executes the coupled Stokes-Darcy hydrothermal circulation benchmark
at 32x32 resolution (Nx=33, Ny=33, dx=dy=4375.0 m).
Couples Darcy thermal buoyancy, Arrhenius fluid viscosity, and dynamic hydrofracturing.

Usage:
    julia --project=. scripts/run_hydrothermal_benchmark.jl [config_path]
"""

using Erebus
using Erebus.Physics: compute_fluid_viscosity, compute_hydrofracture_permeability
using JLD2
using LinearAlgebra
using Printf

function main()
    config_path = if length(ARGS) >= 1
        ARGS[1]
    else
        joinpath(@__DIR__, "..", "configs", "hydrothermal_benchmark.toml")
    end
    println("="^70)
    println("  Erebus.jl 2D Hydrothermal Circulation Benchmark (32x32)")
    println("="^70)
    println("Loading configuration: $config_path")

    cfg = load_config(config_path)
    output_dir = cfg.output.output_dir
    println("Output directory: $output_dir")
    println("Domain size: $(cfg.grid.xsize) m x $(cfg.grid.ysize) m")
    println(
        "Grid resolution: $(cfg.grid.Nx) x $(cfg.grid.Ny) nodes (32x32 cells, dx=$(cfg.grid.xsize/(cfg.grid.Nx-1)) m)",
    )
    println(
        "Planetesimal: R_planet=$(cfg.geometry.rplanet) m, R_crust=$(cfg.geometry.rcrust) m"
    )
    println("Thermal buoyancy: $(cfg.thermodynamics.thermal_buoyancy)")
    println("Fluid viscosity mode: $(cfg.thermodynamics.fluid_viscosity_mode)")
    println("Dynamic hydrofracturing: $(cfg.poroelasticity.hydrofracture)")
    println(
        "Poroelastic compressibility: beta_s=$(cfg.poroelasticity.betasolid) Pa⁻¹, beta_f=$(cfg.poroelasticity.betafluid) Pa⁻¹",
    )
    println("-"^70)

    # Execute simulation
    println("Launching simulation loop...")
    run_simulation(config_path)
    println("Simulation execution completed.")
    println("-"^70)

    # Post-process results
    println("Diagnostic Analysis of Output Checkpoints:")
    files = sort(
        filter(f -> startswith(f, "output_") && endswith(f, ".jld2"), readdir(output_dir))
    )

    @printf(
        "%-10s %-12s %-12s %-14s %-14s %-12s\n",
        "Timestep",
        "Time [Ma]",
        "T_max [K]",
        "max(|qD|) [m/s]",
        "max(k_eff) [m²]",
        "P_eff_min [MPa]"
    )
    @printf(
        "%-10s %-12s %-12s %-14s %-14s %-12s\n",
        "-"^8,
        "-"^10,
        "-"^10,
        "-"^12,
        "-"^12,
        "-"^10
    )

    for f in files
        fpath = joinpath(output_dir, f)
        data = load_state(fpath)
        step = data["timestep"]
        if step == 0
            continue
        end
        t_Ma = data["timesum"] / (365.25 * 86400 * 1e6)
        tk = data["tk2"]
        t_max = maximum(tk)

        # Darcy flux magnitude interpolated to cell centers (P nodes)
        qx = data["qxD"]
        qy = data["qyD"]
        Ny1, Nx1 = size(tk)
        qxp = zeros(Ny1, Nx1)
        qyp = zeros(Ny1, Nx1)
        for i in 1:Ny1, j in 2:Nx1
            qxp[i, j] = 0.5 * (qx[i, j] + qx[i, j - 1])
        end
        for i in 2:Ny1, j in 1:Nx1
            qyp[i, j] = 0.5 * (qy[i, j] + qy[i - 1, j])
        end
        q_mag = sqrt.(qxp .^ 2 .+ qyp .^ 2)
        max_q = maximum(q_mag)

        # Effective permeability from RX = eta_f / k_eff evaluated with local eta_f(T)
        rx = data["RX"]
        keff_vals = Float64[]
        for i in 1:Ny1, j in 2:(Nx1 - 1)
            rx_val = rx[i, j]
            if rx_val > 0.0
                T_vx = 0.5 * (tk[i, j] + tk[i, j + 1])
                eta_f = compute_fluid_viscosity(
                    T_vx,
                    1;
                    mode=Symbol(cfg.thermodynamics.fluid_viscosity_mode),
                    Ea=cfg.thermodynamics.fluid_viscosity_Ea,
                    T0=cfg.thermodynamics.fluid_viscosity_T0,
                    eta0=cfg.thermodynamics.fluid_viscosity_eta0,
                    tmfluidphase=cfg.thermodynamics.tmfluidphase,
                )
                push!(keff_vals, eta_f / rx_val)
            end
        end
        max_k = isempty(keff_vals) ? 0.0 : maximum(keff_vals)

        # Effective stress (Pt - Pf)
        pt = data["pr"]
        pf = data["pf"]
        peff = pt .- pf
        min_peff_MPa = minimum(peff) * 1e-6

        @printf(
            "%-10d %-12.4f %-12.2f %-14.4e %-14.4e %-12.3f\n",
            step,
            t_Ma,
            t_max,
            max_q,
            max_k,
            min_peff_MPa
        )
    end
    return println("="^70)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

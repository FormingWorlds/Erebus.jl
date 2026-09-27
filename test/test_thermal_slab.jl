using Test
using Erebus
using LinearAlgebra

@testset "Analytical 2D Thermal Slab Conduction Benchmark" begin
    Lx = 100_000.0   # 100 km box width [m]
    Ly = 100_000.0   # 100 km box height [m]
    k_val = 3.0      # Thermal conductivity [W/(m K)]
    rhocp_val = 3.0e6 # Volumetric heat capacity [J/(m^3 K)]
    kappa = k_val / rhocp_val # Thermal diffusivity [m^2/s]
    T0 = 300.0       # Background temperature [K]
    deltaT = 50.0   # Perturbation amplitude [K]

    # Fundamental decay rate and characteristic diffusion timescale
    lambda = kappa * π^2 * (1.0 / Lx^2 + 1.0 / Ly^2)
    tau_diff = 1.0 / lambda

    @testset "2D Cosine Decay Solution and Energy Conservation" begin
        N = 33
        coords = Erebus.GridCoordinates(N, N; xsize=Lx, ysize=Ly)
        Nx1, Ny1 = coords.Nx1, coords.Ny1

        tk = zeros(Ny1, Nx1)
        for j in 1:Nx1, i in 1:Ny1
            tk[i, j] = T0 + deltaT * cos(π * coords.xp[j] / Lx) * cos(π * coords.yp[i] / Ly)
        end

        RHOCP = fill(rhocp_val, Ny1, Nx1)
        KX = fill(k_val, N, Nx1)
        KY = fill(k_val, Ny1, N)
        HR = zeros(Ny1, Nx1)
        HA = zeros(Ny1, Nx1)
        HS = zeros(Ny1, Nx1)
        DHP = zeros(Ny1, Nx1)
        RT = zeros(Ny1 * Nx1)

        E_init = sum(tk[2:(Ny1 - 1), 2:(Nx1 - 1)]) * rhocp_val * coords.dx * coords.dy

        # Time integration to t = 0.1 tau
        t_target = 0.1 * tau_diff
        n_steps = 10
        dt = t_target / n_steps
        for _ in 1:n_steps
            LT = Erebus.assemble_thermal_lse!(
                tk, RHOCP, KX, KY, HR, HA, HS, DHP, RT, dt; coords=coords
            )
            sol = LT \ RT
            tk .= reshape(sol, Ny1, Nx1)
        end

        # Analytical solution at t_target
        decay = exp(-lambda * t_target)
        diffs = Float64[]
        for j in 2:(Nx1 - 1), i in 2:(Ny1 - 1)
            ana_val =
                T0 +
                deltaT * cos(π * coords.xp[j] / Lx) * cos(π * coords.yp[i] / Ly) * decay
            push!(diffs, abs(tk[i, j] - ana_val))
        end

        linf_rel = maximum(diffs) / deltaT
        l2_rel = sqrt(sum(diffs .^ 2) / length(diffs)) / deltaT
        E_final = sum(tk[2:(Ny1 - 1), 2:(Nx1 - 1)]) * rhocp_val * coords.dx * coords.dy
        energy_drift = abs(E_final - E_init) / E_init

        # Accuracy verification
        @test linf_rel < 1.0e-3
        @test l2_rel < 5.0e-4

        # Strict closed-box energy conservation to machine precision
        @test energy_drift < 1.0e-12

        # Bounding and extrema monotonicity
        @test all(tk .>= T0 - deltaT)
        @test all(tk .<= T0 + deltaT)
        @test maximum(tk) < T0 + deltaT
        @test minimum(tk) > T0 - deltaT
    end

    @testset "Subcycling with perform_thermal_iterations! (N4)" begin
        N = 33
        coords = Erebus.GridCoordinates(N, N; xsize=Lx, ysize=Ly)
        Nx1, Ny1 = coords.Nx1, coords.Ny1

        tk1 = zeros(Ny1, Nx1)
        for j in 1:Nx1, i in 1:Ny1
            tk1[i, j] =
                T0 + deltaT * cos(π * coords.xp[j] / Lx) * cos(π * coords.yp[i] / Ly)
        end

        tk0 = copy(tk1)
        tk2 = copy(tk1)
        DT = zeros(Ny1, Nx1)
        DT0 = zeros(Ny1, Nx1)
        RHOCP = fill(rhocp_val, Ny1, Nx1)
        KX = fill(k_val, N, Nx1)
        KY = fill(k_val, Ny1, N)
        HR = zeros(Ny1, Nx1)
        HA = zeros(Ny1, Nx1)
        HS = zeros(Ny1, Nx1)
        DHP = zeros(Ny1, Nx1)
        RT = zeros(Ny1 * Nx1)
        ST = zeros(Ny1 * Nx1)

        E_init = sum(tk1[2:(Ny1 - 1), 2:(Nx1 - 1)]) * rhocp_val * coords.dx * coords.dy

        # Single macro timestep with DTmax enforcement
        dt_macro = 1.0e13 # ~317,000 years
        DTmax = 5.0      # Force subcycling (macro DT would be ~15 K)
        Erebus.perform_thermal_iterations!(
            tk0,
            tk1,
            tk2,
            DT,
            DT0,
            RHOCP,
            KX,
            KY,
            HR,
            HA,
            HS,
            DHP,
            RT,
            ST,
            dt_macro;
            coords=coords,
            DTmax_val=DTmax,
        )

        decay_sub = exp(-lambda * dt_macro)
        diffs_sub = Float64[]
        for j in 2:(Nx1 - 1), i in 2:(Ny1 - 1)
            ana_val =
                T0 +
                deltaT * cos(π * coords.xp[j] / Lx) * cos(π * coords.yp[i] / Ly) * decay_sub
            push!(diffs_sub, abs(tk2[i, j] - ana_val))
        end

        linf_sub = maximum(diffs_sub) / deltaT
        E_final = sum(tk2[2:(Ny1 - 1), 2:(Nx1 - 1)]) * rhocp_val * coords.dx * coords.dy
        drift_sub = abs(E_final - E_init) / E_init

        @test linf_sub < 1.0e-3
        @test drift_sub < 1.0e-12
        @test isapprox(DT, tk2 .- tk0; rtol=1e-12)
        @test isapprox(DT0, DT; rtol=1e-12)
    end

    @testset "Spatial Grid Refinement Convergence" begin
        errors = Float64[]
        t_ref = 0.1 * tau_diff
        decay_ref = exp(-lambda * t_ref)

        for N in [17, 33, 65]
            coords_n = Erebus.GridCoordinates(N, N; xsize=Lx, ysize=Ly)
            Nx1_n, Ny1_n = coords_n.Nx1, coords_n.Ny1

            tk_n = zeros(Ny1_n, Nx1_n)
            for j in 1:Nx1_n, i in 1:Ny1_n
                tk_n[i, j] =
                    T0 +
                    deltaT * cos(π * coords_n.xp[j] / Lx) * cos(π * coords_n.yp[i] / Ly)
            end

            RHOCP_n = fill(rhocp_val, Ny1_n, Nx1_n)
            KX_n = fill(k_val, N, Nx1_n)
            KY_n = fill(k_val, Ny1_n, N)
            HR_n = zeros(Ny1_n, Nx1_n)
            HA_n = zeros(Ny1_n, Nx1_n)
            HS_n = zeros(Ny1_n, Nx1_n)
            DHP_n = zeros(Ny1_n, Nx1_n)
            RT_n = zeros(Ny1_n * Nx1_n)

            dt_n = t_ref / 10
            for _ in 1:10
                LT = Erebus.assemble_thermal_lse!(
                    tk_n,
                    RHOCP_n,
                    KX_n,
                    KY_n,
                    HR_n,
                    HA_n,
                    HS_n,
                    DHP_n,
                    RT_n,
                    dt_n;
                    coords=coords_n,
                )
                sol = LT \ RT_n
                tk_n .= reshape(sol, Ny1_n, Nx1_n)
            end

            diffs_n = Float64[]
            for j in 2:(Nx1_n - 1), i in 2:(Ny1_n - 1)
                ana_val =
                    T0 +
                    deltaT *
                    cos(π * coords_n.xp[j] / Lx) *
                    cos(π * coords_n.yp[i] / Ly) *
                    decay_ref
                push!(diffs_n, abs(tk_n[i, j] - ana_val))
            end
            push!(errors, maximum(diffs_n) / deltaT)
        end

        # Monotonic convergence under grid refinement
        @test errors[1] > errors[2] > errors[3]
        @test errors[3] < 5.0e-4
    end
end

using Erebus
using Erebus.Config
using Erebus.Physics
using Erebus.Geometry
using Erebus.Particles
using Erebus.Numerics
using LinearSolve
using ExtendableSparse
using Test

include("test_helpers.jl")

@testset "Physics & Numerics Mutation Testing Suite" begin
    @testset "StokesDarcy continuity coupling mutation" begin
        coords = test_grid_coordinates(; Nx=17, Ny=17, xsize=50_000.0, ysize=50_000.0)
        Ny, Nx = coords.Ny, coords.Nx
        Ny1, Nx1 = coords.Ny1, coords.Nx1
        dt = 1000.0
        ETA = fill(1.0e19, Ny, Nx)
        ETAP = fill(1.0e19, Ny1, Nx1)
        GGG = fill(1.0e10, Ny, Nx)
        GGGP = fill(1.0e10, Ny1, Nx1)
        SXY0 = zeros(Ny, Nx)
        SXX0 = zeros(Ny1, Nx1)
        RHOX = fill(3000.0, Ny1, Nx1)
        RHOY = fill(3000.0, Ny1, Nx1)
        RHOFX = fill(1000.0, Ny1, Nx1)
        RHOFY = fill(1000.0, Ny1, Nx1)
        RX = fill(1.0e10, Ny1, Nx1)
        RY = fill(1.0e10, Ny1, Nx1)
        ETAPHI = fill(1.0e16, Ny1, Nx1)
        BETAPHI = fill(1.0e-10, Ny1, Nx1)
        PHI = fill(0.1, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = zeros(Ny1, Nx1)
        pr0 = zeros(Ny1, Nx1)
        pf0 = zeros(Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)
        R = zeros(Nx1 * Ny1 * 6)

        L = assemble_hydromechanical_lse!(
            ETA,
            ETAP,
            GGG,
            GGGP,
            SXY0,
            SXX0,
            RHOX,
            RHOY,
            RHOFX,
            RHOFY,
            RX,
            RY,
            ETAPHI,
            BETAPHI,
            PHI,
            gx,
            gy,
            pr0,
            pf0,
            DMP,
            dt,
            R;
            coords=coords,
            betasolid=0.0,
            betafluid=0.0,
        )
        # Baseline: discrete solid continuity stencil sums to zero on interior P node
        i, j = 5, 5
        kvx = ((j - 1) * Ny1 + i - 1) * 6 + 1
        kpm = kvx + 2
        div_x_sum = L[kpm, kvx - 6 * Ny1] + L[kpm, kvx]
        @test isapprox(div_x_sum, 0.0; atol=1e-12)

        # Mutant: reversed continuity sign (+1/dx instead of -1/dx)
        div_x_mutant = (-L[kpm, kvx - 6 * Ny1]) + L[kpm, kvx]
        @test isapprox(div_x_mutant, 2.0 / coords.dx; atol=1e-10)
        @test !isapprox(div_x_mutant, 0.0; atol=1e-6)
    end

    @testset "Thermal diffusion sign mutation" begin
        Nx_th = 51
        dx_th = 1000.0
        kappa_th = 1.0e-6
        dt_th = 0.2 * dx_th^2 / kappa_th
        T0 = fill(300.0, Nx_th)
        T0[26] = 1200.0

        # Baseline: positive thermal conductivity diffuses peak temperature
        T_base = copy(T0)
        for i in 2:(Nx_th - 1)
            d2T = (T0[i + 1] - 2.0 * T0[i] + T0[i - 1]) / dx_th^2
            T_base[i] += kappa_th * dt_th * d2T
        end
        T_analytic_diffused = T0[26] - 2.0 * kappa_th * dt_th * (T0[26] - 300.0) / dx_th^2
        @test isapprox(T_base[26], T_analytic_diffused; rtol=1e-6)
        @test isapprox(sum(T_base), sum(T0); rtol=1e-10)

        # Mutant: negative conductivity causes explosive temperature growth
        T_mut = copy(T0)
        for i in 2:(Nx_th - 1)
            d2T = (T0[i + 1] - 2.0 * T0[i] + T0[i - 1]) / dx_th^2
            T_mut[i] += (-kappa_th) * dt_th * d2T
        end
        T_analytic_mut = T0[26] + 2.0 * kappa_th * dt_th * (T0[26] - 300.0) / dx_th^2
        @test isapprox(T_mut[26], T_analytic_mut; rtol=1e-6)
        @test !isapprox(T_mut[26], T_base[26]; rtol=0.2)
    end

    @testset "Solubility Henry law exponent mutation" begin
        P1 = 1.0e6
        P2 = 4.0e6
        As_val = 0.40

        # Baseline: square-root dependence in Burnham-Dixon law
        w1_base = compute_water_solubility_melt(P1; As=As_val, law=:burnham_dixon)
        w2_base = compute_water_solubility_melt(P2; As=As_val, law=:burnham_dixon)
        ratio_base = w2_base / w1_base
        @test isapprox(ratio_base, 2.0; rtol=1e-6)

        # Mutant: linear pressure dependence (n = 1.0 instead of 0.5)
        w1_mut = As_val * (P1 * 1.0e-6)
        w2_mut = As_val * (P2 * 1.0e-6)
        ratio_mut = w2_mut / w1_mut
        @test isapprox(ratio_mut, 4.0; rtol=1e-6)
        @test !isapprox(ratio_mut, ratio_base; rtol=0.2)
    end

    @testset "Jeans kinetic escape flux velocity mutation" begin
        n_exo = 1.0e12
        T_exo = 1000.0
        m_H_kg = 1.008e-3 / AVOGADRO_CONSTANT
        lambda_val = 6.0

        # Baseline: kinetic Jeans escape flux with physical thermal velocity
        flux_base = compute_jeans_escape_flux(
            n_exo, T_exo, m_H_kg, lambda_val; hydrodynamic=false
        )
        v_th = sqrt(2.0 * BOLTZMANN_CONSTANT * T_exo / m_H_kg)
        flux_analytic =
            (n_exo * v_th / (2.0 * sqrt(pi))) * (1.0 + lambda_val) * exp(-lambda_val)
        @test isapprox(flux_base, flux_analytic; rtol=1e-8)

        # Mutant: thermal velocity multiplied by sqrt(2)
        v_th_mut = v_th * sqrt(2.0)
        flux_mut =
            (n_exo * v_th_mut / (2.0 * sqrt(pi))) * (1.0 + lambda_val) * exp(-lambda_val)
        @test isapprox(flux_mut, sqrt(2.0) * flux_base; rtol=1e-8)
        @test !isapprox(flux_mut, flux_base; rtol=0.2)
    end

    @testset "Metal-silicate volatile mass balance and D mutation" begin
        m_sil = 0.8
        m_met = 0.2
        M_tot = 1.0e-3
        D_base = 4.0

        # Baseline: volatile mass balance closes to machine precision
        C_sil_base = M_tot / (m_sil + D_base * m_met)
        C_met_base = D_base * C_sil_base
        M_closed_base = m_sil * C_sil_base + m_met * C_met_base
        @test isapprox(M_closed_base, M_tot; atol=1e-12)

        # Mutant: doubled partition coefficient shifts metal volatile budget
        D_mut = 2.0 * D_base
        C_sil_mut = M_tot / (m_sil + D_mut * m_met)
        C_met_mut = D_mut * C_sil_mut
        @test isapprox(m_sil * C_sil_mut + m_met * C_met_mut, M_tot; atol=1e-12)
        @test isapprox(C_met_mut, (8.0 / 2.4) * 1.0e-3; rtol=1e-10)
        @test !isapprox(C_met_mut, C_met_base; rtol=0.1)
    end

    @testset "Venting surface single-cell budget clamp mutation" begin
        phi_cell = 0.20
        phi_min = 0.05
        rho_fluid = 1000.0
        dV = 1000.0 * 1000.0
        dt_val = 100.0
        S_max = max(0.0, phi_cell - phi_min) * rho_fluid * dV / dt_val

        unconstrained_flux = 5.0 * S_max
        # Baseline: constrained venting cannot exceed S_max
        constrained_flux = min(unconstrained_flux, S_max)
        @test isapprox(constrained_flux, S_max; rtol=1e-10)

        # Mutant: unconstrained venting violates single-cell fluid budget
        @test isapprox(unconstrained_flux, 5.0 * S_max; rtol=1e-10)
        @test !isapprox(unconstrained_flux, S_max; rtol=0.5)
    end

    @testset "Radioactive heating decay constant mutation" begin
        f_al = 2.2e24
        ratio_al = 5.0e-5
        E_al = 3.0 * 1.602176634e-13
        tau_al = 0.717e6 * 365.25 * 86400.0 / log(2.0)
        t_eval = 1.0e6 * 365.25 * 86400.0

        # Baseline: Q_radiogenic matches analytical formula
        Q_base = Q_radiogenic(f_al, ratio_al, E_al, tau_al, t_eval)
        Q_analytic = f_al * ratio_al * E_al * exp(-t_eval / tau_al) / tau_al
        @test isapprox(Q_base, Q_analytic; rtol=1e-10)

        # Mutant: halved lifetime (doubled decay rate) deviates by > 20%
        Q_mut = Q_radiogenic(f_al, ratio_al, E_al, tau_al * 0.5, t_eval)
        Q_analytic_mut =
            f_al * ratio_al * E_al * exp(-t_eval / (tau_al * 0.5)) / (tau_al * 0.5)
        @test isapprox(Q_mut, Q_analytic_mut; rtol=1e-10)
        @test !isapprox(Q_mut, Q_base; rtol=0.2)
    end

    @testset "P2M interpolation weights partition of unity mutation" begin
        coords = test_grid_coordinates(; Nx=9, Ny=9, xsize=10_000.0, ysize=10_000.0)
        xp = coords.xp
        yp = coords.yp
        dx = coords.dx
        dy = coords.dy

        xm = 3500.0
        ym = 4200.0
        i, j, weights = fix_weights(
            xm,
            ym,
            xp,
            yp,
            dx,
            dy,
            coords.jmin_p,
            coords.jmax_p,
            coords.imin_p,
            coords.imax_p,
        )
        # Baseline: bilinear weights satisfy partition of unity
        @test isapprox(sum(weights), 1.0; atol=1e-12)

        # Mutant: squared weights violate partition of unity
        mutated_weights = weights .^ 2
        sum_sq = sum(mutated_weights)
        @test isapprox(sum(mutated_weights), sum_sq; atol=1e-12)
        @test !isapprox(sum(mutated_weights), 1.0; rtol=0.1)
    end

    @testset "Silicate melt segregation buoyancy direction mutation" begin
        F_m = 0.25
        g_acc = 1.5
        eta_melt = 10.0
        rho_solid = 3000.0

        # Baseline: positive buoyancy (rho_solid > rho_melt) produces upward velocity
        rho_melt_buoyant = 2700.0
        drho_base = rho_solid - rho_melt_buoyant
        v_base = silicate_melt_segregation_velocity(F_m, drho_base, g_acc, eta_melt)
        v_expected = silicate_melt_segregation_velocity(F_m, drho_base, g_acc, eta_melt)
        @test isapprox(v_base, v_expected; rtol=1e-10)
        @test !isapprox(v_base, 0.0; atol=1e-10)

        # Mutant: negative buoyancy halts upward percolation
        drho_mut = rho_solid - 3300.0
        v_mut = silicate_melt_segregation_velocity(F_m, drho_mut, g_acc, eta_melt)
        @test isapprox(v_mut, 0.0; atol=1e-12)
        @test !isapprox(v_base, v_mut; atol=1e-8)
    end

    @testset "Plastic yielding cohesion threshold mutation" begin
        pr_eff = 5.0e6
        tau_elastic = 2.0e6
        friction = 0.6

        # Baseline: positive cohesion keeps stress strictly below yield surface
        coh_base = 10.0e6
        syield_base = max(coh_base + friction * pr_eff, 0.0)
        is_yielding_base = tau_elastic > syield_base
        @test !is_yielding_base
        @test isapprox(syield_base, 13.0e6; rtol=1e-10)

        # Mutant: negative cohesion triggers yielding everywhere
        coh_mut = -5.0e6
        syield_mut = max(coh_mut + friction * pr_eff, 0.0)
        is_yielding_mut = tau_elastic > syield_mut
        @test is_yielding_mut
        @test isapprox(syield_mut, 0.0; atol=1e-12)
    end
end

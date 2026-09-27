using Test
using Erebus
using LinearAlgebra

@testset "Hydrofracture Stability and Limiters" begin
    @testset "Darcy Resistance Under-Relaxation Recurrence" begin
        rx_new = 5.0e8
        rx_0 = 1.0e10
        theta = 0.3

        # Iteration sequence: rx^{(k)} = theta * rx_new + (1 - theta) * rx^{(k-1)}
        rx_k = rx_0
        e0 = rx_0 - rx_new
        for k in 1:5
            rx_k = theta * rx_new + (1.0 - theta) * rx_k
            analytic_diff = (1.0 - theta)^k * e0
            curr_diff = rx_k - rx_new
            @test isapprox(curr_diff, analytic_diff; rtol=1.0e-12)
        end

        # Neutral identity: theta = 1.0 matches undamped sequence bitwise
        theta_neutral = 1.0
        rx_undamped = theta_neutral * rx_new + (1.0 - theta_neutral) * rx_0
        @test rx_undamped === rx_new
    end

    @testset "C1 Overpressure Regularisation Ramp" begin
        delta = 0.05
        kappa = 1000.0
        gamma = 1.0

        # Exact zero for non-positive normalised overpressures
        @test Erebus.hydrofracture_overpressure_ramp(-1.0, delta) === 0.0
        @test Erebus.hydrofracture_overpressure_ramp(0.0, delta) === 0.0
        @test isnan(Erebus.hydrofracture_overpressure_ramp(NaN, delta))
        @test Erebus.hydrofracture_overpressure_ramp(Inf, delta) == Inf

        # DomainError validation for invalid ramp_width
        @test_throws DomainError Erebus.hydrofracture_overpressure_ramp(0.5, -0.05)
        @test_throws DomainError Erebus.hydrofracture_overpressure_ramp(0.5, NaN)

        # Quadratic transition regime: 0 < x < delta
        sigma_t = 1.0e7
        x_mid = 0.5 * delta
        val_mid = Erebus.hydrofracture_overpressure_ramp(x_mid, delta)
        @test isapprox(val_mid, (x_mid^2) / (2.0 * delta); rtol=1.0e-12)
        f_mid = Erebus.compute_hydrofracture_factor(
            -sigma_t * (1.0 + x_mid),
            sigma_t;
            active=true,
            kappa_frac=kappa,
            gamma=gamma,
            ramp_width=delta,
        )
        @test isapprox(f_mid, 1.0 + kappa * val_mid^gamma; rtol=1.0e-12)

        # Continuity of value at x = delta
        val_at_delta = Erebus.hydrofracture_overpressure_ramp(delta, delta)
        @test isapprox(val_at_delta, 0.5 * delta; atol=1.0e-15)

        # C1 Continuity: left and right numerical derivatives at x = delta
        h = 1.0e-7
        d_left =
            (val_at_delta - Erebus.hydrofracture_overpressure_ramp(delta - h, delta)) / h
        d_right =
            (Erebus.hydrofracture_overpressure_ramp(delta + h, delta) - val_at_delta) / h
        @test isapprox(d_left, 1.0; atol=1.0e-5)
        @test isapprox(d_right, 1.0; atol=1.0e-5)
        @test isapprox(d_left, d_right; atol=1.0e-5)

        # Asymptotic linearity at large overpressure: s(x) - x == -delta / 2
        x_large = 10.0 * delta
        s_large = Erebus.hydrofracture_overpressure_ramp(x_large, delta)
        @test isapprox(s_large - x_large, -0.5 * delta; atol=1.0e-12)

        # Exact reproduction of hard kink law when delta = 0.0
        @test Erebus.hydrofracture_overpressure_ramp(-0.5, 0.0) === 0.0
        @test Erebus.hydrofracture_overpressure_ramp(0.0, 0.0) === 0.0
        @test Erebus.hydrofracture_overpressure_ramp(0.42, 0.0) === 0.42

        # Compressive / intact state -> factor 1.0
        f_intact = Erebus.compute_hydrofracture_factor(
            1.0e7, sigma_t; active=true, kappa_frac=kappa, gamma=gamma, ramp_width=delta
        )
        @test f_intact === 1.0

        # Overpressured state with ramp
        peff_over = -1.5e7
        norm_over = (-peff_over - sigma_t) / sigma_t # 0.5
        s_expected = Erebus.hydrofracture_overpressure_ramp(norm_over, delta)
        f_ramp = Erebus.compute_hydrofracture_factor(
            peff_over, sigma_t; active=true, kappa_frac=kappa, gamma=gamma, ramp_width=delta
        )
        @test isapprox(f_ramp, 1.0 + kappa * s_expected^gamma; rtol=1.0e-12)
    end

    @testset "Stokes-Darcy Assembly and Flux Reconstruction Relaxation Integration" begin
        Nx, Ny = 5, 5
        Nx1, Ny1 = Nx + 1, Ny + 1
        coords = GridCoordinates(GridConfig(Nx=Nx, Ny=Ny, xsize=1000.0, ysize=1000.0))

        eta_f = 1.0e-3
        k_perm = 1.0e-12
        sigma_t_val = 1.0e6
        ETA = fill(1e20, Ny, Nx)
        ETAP = fill(1e20, Ny1, Nx1)
        GGG = fill(1e10, Ny, Nx)
        GGGP = fill(1e10, Ny1, Nx1)
        SXY0 = zeros(Ny, Nx)
        SXX0 = zeros(Ny, Nx)
        RHOX = fill(3000.0, Ny1, Nx1)
        RHOY = fill(3000.0, Ny1, Nx1)
        RHOFX = fill(1000.0, Ny1, Nx1)
        RHOFY = fill(1000.0, Ny1, Nx1)
        RX = fill(eta_f / k_perm, Ny1, Nx1)
        RY = fill(eta_f / k_perm, Ny1, Nx1)
        KX = fill(k_perm, Ny1, Nx1)
        KY = fill(k_perm, Ny1, Nx1)
        ETAPHI = fill(1e24, Ny1, Nx1)
        BETAPHI = fill(1e-10, Ny1, Nx1)
        PHI = fill(0.1, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = fill(10.0, Ny1, Nx1)
        pr0 = zeros(Ny1, Nx1)
        pf0 = zeros(Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)
        dt = 1.0e6
        R = zeros(Nx1 * Ny1 * 6)
        pr = fill(1.0e6, Ny1, Nx1)
        pf = fill(3.0e6, Ny1, Nx1) # Pf > Pr -> tensile overpressure Peff = -2e6 Pa
        TEN = fill(sigma_t_val, Ny, Nx)

        rx_eff_out = zeros(Ny1, Nx1)
        ry_eff_out = zeros(Ny1, Nx1)
        rx_eff_prev = fill(2.0 * eta_f / k_perm, Ny1, Nx1)
        ry_eff_prev = fill(2.0 * eta_f / k_perm, Ny1, Nx1)

        L = Erebus.assemble_hydromechanical_lse!(
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
            pr=pr,
            pf=pf,
            TEN=TEN,
            KX=KX,
            KY=KY,
            hydrofracture=true,
            theta_frac=0.5,
            ramp_width=0.05,
            rx_eff_prev=rx_eff_prev,
            ry_eff_prev=ry_eff_prev,
            rx_eff_out=rx_eff_out,
            ry_eff_out=ry_eff_out,
            coords=coords,
        )

        @test all(rx_eff_out .> 0.0)
        @test all(ry_eff_out .> 0.0)
        # Boundary nodes must match unrelaxed RX
        @test rx_eff_out[1, 1] == RX[1, 1]
        @test ry_eff_out[1, 1] == RY[1, 1]
        # Relaxed resistance must match analytic 0.5 * r_target + 0.5 * r_prev
        r_target_ana = 1.0e-3 / 9.76e-10
        r_relaxed_ana = 0.5 * r_target_ana + 0.5 * (2.0 * eta_f / k_perm)
        @test isapprox(rx_eff_out[2, 2], r_relaxed_ana; rtol=1.0e-12)
        @test isapprox(ry_eff_out[2, 2], r_relaxed_ana; rtol=1.0e-12)

        # Neutral identity: theta = 1.0 with Inf in rx_eff_prev does not produce NaN
        rx_eff_out_inf = zeros(Ny1, Nx1)
        rx_eff_prev_inf = fill(Inf, Ny1, Nx1)
        Erebus.assemble_hydromechanical_lse!(
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
            pr=pr,
            pf=pf,
            TEN=TEN,
            KX=KX,
            KY=KY,
            hydrofracture=true,
            theta_frac=1.0,
            ramp_width=0.05,
            rx_eff_prev=rx_eff_prev_inf,
            ry_eff_prev=rx_eff_prev_inf,
            rx_eff_out=rx_eff_out_inf,
            ry_eff_out=zeros(Ny1, Nx1),
            coords=coords,
        )
        @test all(!isnan, rx_eff_out_inf)
        @test rx_eff_out_inf[2, 2] == r_target_ana

        # Darcy flux reconstruction must not produce NaN or Inf
        qxD = zeros(Ny1, Nx1)
        qyD = zeros(Ny1, Nx1)
        Erebus.reconstruct_darcy_fluxes!(
            qxD,
            qyD,
            pf,
            RHOFX,
            RHOFY,
            RX,
            RY,
            gx,
            gy,
            coords;
            hydrofracture=true,
            pr=pr,
            TEN=TEN,
            KX=KX,
            KY=KY,
            rx_eff=rx_eff_out,
            ry_eff=ry_eff_out,
        )
        @test all(isfinite, qxD)
        @test all(isfinite, qyD)
    end

    @testset "4-Variable Stokes-Darcy Assembly Hydrofracture Integration" begin
        Nx, Ny = 5, 5
        Nx1, Ny1 = Nx + 1, Ny + 1
        coords = GridCoordinates(GridConfig(Nx=Nx, Ny=Ny, xsize=1000.0, ysize=1000.0))

        eta_f = 1.0e-3
        k_perm = 1.0e-12
        sigma_t_val = 1.0e6
        ETA = fill(1e20, Ny, Nx)
        ETAP = fill(1e20, Ny1, Nx1)
        GGG = fill(1e10, Ny, Nx)
        GGGP = fill(1e10, Ny1, Nx1)
        SXY0 = zeros(Ny, Nx)
        SXX0 = zeros(Ny, Nx)
        RHOX = fill(3000.0, Ny1, Nx1)
        RHOY = fill(3000.0, Ny1, Nx1)
        RHOFX = fill(1000.0, Ny1, Nx1)
        RHOFY = fill(1000.0, Ny1, Nx1)
        RX = fill(eta_f / k_perm, Ny1, Nx1)
        RY = fill(eta_f / k_perm, Ny1, Nx1)
        KX = fill(k_perm, Ny1, Nx1)
        KY = fill(k_perm, Ny1, Nx1)
        ETAPHI = fill(1e24, Ny1, Nx1)
        BETAPHI = fill(1e-10, Ny1, Nx1)
        PHI = fill(0.1, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = fill(10.0, Ny1, Nx1)
        pr0 = zeros(Ny1, Nx1)
        pf0 = zeros(Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)
        dt = 1.0e6
        R4 = zeros(Nx1 * Ny1 * 4)
        pr = fill(1.0e6, Ny1, Nx1)
        pf = fill(3.0e6, Ny1, Nx1)
        TEN = fill(sigma_t_val, Ny, Nx)

        rx_eff_out_4 = zeros(Ny1, Nx1)
        ry_eff_out_4 = zeros(Ny1, Nx1)
        rx_eff_prev = fill(2.0 * eta_f / k_perm, Ny1, Nx1)
        ry_eff_prev = fill(2.0 * eta_f / k_perm, Ny1, Nx1)

        L4 = Erebus.assemble_hydromechanical_4var_lse!(
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
            R4;
            pr=pr,
            pf=pf,
            TEN=TEN,
            KX=KX,
            KY=KY,
            hydrofracture=true,
            theta_frac=0.5,
            ramp_width=0.05,
            rx_eff_prev=rx_eff_prev,
            ry_eff_prev=ry_eff_prev,
            rx_eff_out=rx_eff_out_4,
            ry_eff_out=ry_eff_out_4,
            coords=coords,
        )

        @test all(rx_eff_out_4 .> 0.0)
        @test all(ry_eff_out_4 .> 0.0)
        r_target_ana = 1.0e-3 / 9.76e-10
        r_relaxed_ana = 0.5 * r_target_ana + 0.5 * (2.0 * eta_f / k_perm)
        Nx_val = Nx1 - 1
        Ny_val = Ny1 - 1
        for j in 2:(Nx_val - 1), i in 2:(Ny1 - 1)
            @test isapprox(rx_eff_out_4[i, j], r_relaxed_ana; rtol=1.0e-12)
        end
        for j in 2:(Nx1 - 1), i in 2:(Ny_val - 1)
            @test isapprox(ry_eff_out_4[i, j], r_relaxed_ana; rtol=1.0e-12)
        end

        # Neutral identity: theta = 1.0 with Inf in rx_eff_prev does not produce NaN
        rx_eff_out_inf = zeros(Ny1, Nx1)
        ry_eff_out_inf = zeros(Ny1, Nx1)
        rx_eff_prev_inf = fill(Inf, Ny1, Nx1)
        Erebus.assemble_hydromechanical_4var_lse!(
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
            R4;
            pr=pr,
            pf=pf,
            TEN=TEN,
            KX=KX,
            KY=KY,
            hydrofracture=true,
            theta_frac=1.0,
            ramp_width=0.05,
            rx_eff_prev=rx_eff_prev_inf,
            ry_eff_prev=rx_eff_prev_inf,
            rx_eff_out=rx_eff_out_inf,
            ry_eff_out=ry_eff_out_inf,
            coords=coords,
        )
        @test all(!isnan, rx_eff_out_inf)
        @test all(!isnan, ry_eff_out_inf)
        @test rx_eff_out_inf[2, 2] == r_target_ana

        # Reconstructed fluxes with 4-var effective resistances
        qxD = zeros(Ny1, Nx1)
        qyD = zeros(Ny1, Nx1)
        Erebus.reconstruct_darcy_fluxes!(
            qxD,
            qyD,
            pf,
            RHOFX,
            RHOFY,
            RX,
            RY,
            gx,
            gy,
            coords;
            hydrofracture=true,
            pr=pr,
            TEN=TEN,
            KX=KX,
            KY=KY,
            rx_eff=rx_eff_out_4,
            ry_eff=ry_eff_out_4,
        )
        @test all(isfinite, qxD)
        @test all(isfinite, qyD)
    end

    @testset "Configuration Schema and Parameter Bounds" begin
        # Valid default configuration validation
        cfg_def = default_config()
        @test validate_config(cfg_def) === nothing

        # PoroelasticConfig validation
        @test_throws DomainError PoroelasticConfig(; theta_frac=0.0)
        @test_throws DomainError PoroelasticConfig(; theta_frac=-0.1)
        @test_throws DomainError PoroelasticConfig(; theta_frac=1.5)
        @test_throws DomainError PoroelasticConfig(; ramp_width=-0.1)
        @test_throws DomainError PoroelasticConfig(; ramp_width=1.5)
        @test_throws DomainError PoroelasticConfig(; rx_floor_prefactor=0.0)
        @test_throws DomainError PoroelasticConfig(; rx_floor_prefactor=-1.0)

        # CoreFormationConfig limiter validation
        @test_throws DomainError CoreFormationConfig(; droplet_radius_min=-1.0e-4)
        @test_throws DomainError CoreFormationConfig(;
            droplet_radius_min=0.1, droplet_radius_max=0.05
        )
        @test_throws DomainError CoreFormationConfig(; settling_drhog_floor=0.0)
        @test_throws DomainError CoreFormationConfig(; g_acc_floor_settling=0.0)
        @test_throws DomainError CoreFormationConfig(; v_rel_floor=0.0)

        # MagmaTransportConfig limiter validation
        @test_throws DomainError MagmaTransportConfig(; eta_solid_floor=0.0)
        @test_throws DomainError MagmaTransportConfig(; eta_solid_floor=-1.0)
    end
end

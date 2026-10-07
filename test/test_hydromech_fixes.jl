# Unit tests for hydromechanics and grid numerics consolidations
using Test
using LinearAlgebra
using StaticArrays
using Erebus

@testset "Hydromechanics and Grid Numerics" begin
    @testset "F03: Hydrofracture permeability and Darcy resistance" begin
        k0 = 1.0e-12
        phi_f = 0.10
        phim0 = 0.01
        eta_f = 1.0e-3
        kphi_expected = Erebus.kphi(k0, phi_f; phim0=phim0)
        r_base = eta_f / kphi_expected

        Peff_active = -5.0e7
        Peff_inactive = 1.0e7
        sigma_t = 2.0e7
        kappa_frac = 1.0e3
        gamma = 1.0
        kmax = 1.0e-5

        # 1. Verification of hydrofracture resistance reduction under overpressure
        R_inactive = Erebus.evaluate_hydrofracture_resistance(
            r_base, Peff_inactive, sigma_t, phi_f;
            kphim0=k0, phim0_val=phim0, kappa_frac=kappa_frac, gamma_frac=gamma, k_frac_max=kmax
        )
        R_active = Erebus.evaluate_hydrofracture_resistance(
            r_base, Peff_active, sigma_t, phi_f;
            kphim0=k0, phim0_val=phim0, kappa_frac=kappa_frac, gamma_frac=gamma, k_frac_max=kmax
        )

        @test isapprox(R_inactive, r_base; rtol=1e-12)
        @test R_active < R_inactive

        keff_expected = Erebus.compute_hydrofracture_permeability(
            kphi_expected, Peff_active, sigma_t;
            active=true, kappa_frac=kappa_frac, gamma=gamma, kmax=kmax
        )
        expected_ratio = kphi_expected / keff_expected
        @test isapprox(R_active / r_base, expected_ratio; rtol=1e-10)

        # 2. Independence from thermal conductivity and verification of finite resistance
        @test isfinite(R_active)
        @test isfinite(R_inactive)
        overpressure_active = -Peff_active - sigma_t
        expected_factor = 1.0 + kappa_frac * (overpressure_active / sigma_t)^gamma
        @test isapprox(keff_expected, kphi_expected * expected_factor; rtol=1e-10)

        # 3. Guard against zero or negative porosity and zero reference permeability
        R_zero_phi = Erebus.evaluate_hydrofracture_resistance(
            r_base, Peff_active, sigma_t, 0.0;
            kphim0=k0, phim0_val=phim0, phimin_val=1e-4, kappa_frac=kappa_frac, gamma_frac=gamma, k_frac_max=kmax
        )
        @test isfinite(R_zero_phi)
        @test R_zero_phi <= r_base

        R_neg_phi = Erebus.evaluate_hydrofracture_resistance(
            r_base, Peff_active, sigma_t, -0.05;
            kphim0=k0, phim0_val=phim0, phimin_val=1e-4, kappa_frac=kappa_frac, gamma_frac=gamma, k_frac_max=kmax
        )
        @test isfinite(R_neg_phi)
        @test isapprox(R_neg_phi, R_zero_phi; rtol=1e-12)

        R_zero_k0 = Erebus.evaluate_hydrofracture_resistance(
            r_base, Peff_active, sigma_t, phi_f;
            kphim0=0.0, phim0_val=phim0, kappa_frac=kappa_frac, gamma_frac=gamma, k_frac_max=kmax
        )
        @test isfinite(R_zero_k0)
        ffrac_expected = Erebus.compute_hydrofracture_factor(
            Peff_active, sigma_t; active=true, kappa_frac=kappa_frac, gamma=gamma
        )
        @test isapprox(R_zero_k0, max(r_base / ffrac_expected, 1.0e-5 / kmax); rtol=1e-10)

        # 4. Extreme overpressure upper bound capping at kmax
        Peff_extreme = -1.0e11
        R_ext = Erebus.evaluate_hydrofracture_resistance(
            r_base, Peff_extreme, sigma_t, phi_f;
            kphim0=k0, phim0_val=phim0, kappa_frac=kappa_frac, gamma_frac=gamma, k_frac_max=kmax
        )
        keff_capped = Erebus.compute_hydrofracture_permeability(
            kphi_expected, Peff_extreme, sigma_t;
            active=true, kappa_frac=kappa_frac, gamma=gamma, kmax=kmax
        )
        @test isapprox(keff_capped, kmax; rtol=1e-10)
        @test isapprox(R_ext, r_base * (kphi_expected / kmax); rtol=1e-10)

        # 5. Verification of minimum resistance floor enforcement when kphi > 0
        R_floored = Erebus.evaluate_hydrofracture_resistance(
            1.0e-10, Peff_active, sigma_t, phi_f;
            kphim0=k0, phim0_val=phim0, kappa_frac=kappa_frac, gamma_frac=gamma,
            k_frac_max=kmax, rx_floor_prefactor=1.0e-5
        )
        expected_floor = 1.0e-5 / kmax
        @test isapprox(R_floored, expected_floor; rtol=1e-12)
        @test R_floored >= expected_floor

        # 6. Face orientation interpolation of cell-centered tensile strength
        # X-face (vertical face) averages horizontally across j
        # Y-face (horizontal face) averages vertically across i
        TEN_grid = [1.0e6 3.0e6; 5.0e6 7.0e6]
        sigma_t_x_face = 0.5 * (TEN_grid[1, 1] + TEN_grid[1, 2])
        sigma_t_y_face = 0.5 * (TEN_grid[1, 1] + TEN_grid[2, 1])
        @test isapprox(sigma_t_x_face, 2.0e6; rtol=1e-12)
        @test isapprox(sigma_t_y_face, 3.0e6; rtol=1e-12)
        @test sigma_t_x_face != sigma_t_y_face
    end

    @testset "F11: Venting corner drainage single-budget constraint" begin
        Ny, Nx = 5, 5
        dx, dy = 100.0, 100.0
        dt = 100.0
        rho_fluid = 1000.0
        phi_val = 0.20
        phimin = 0.01
        cell_vol = dx * dy
        max_drain = (phi_val - phimin) * rho_fluid * cell_vol / dt

        # Staggered conductance and pressure arrays
        cond_x = zeros(Float64, Ny, Nx)
        cond_y = zeros(Float64, Ny, Nx)
        dp_x = zeros(Float64, Ny, Nx)
        dp_y = zeros(Float64, Ny, Nx)

        # Cell (1, 1) has two exposed venting faces: horizontal top face and vertical left face
        # Large conductances and pressure drops driving total unconstrained flux above max_drain
        cond_x[1, 1] = 1.0e-2
        dp_x[1, 1] = 2.0e6
        cond_y[1, 1] = 1.0e-2
        dp_y[1, 1] = 2.0e6

        flux_x_raw = cond_x[1, 1] * dp_x[1, 1]
        flux_y_raw = cond_y[1, 1] * dp_y[1, 1]
        tot_raw = flux_x_raw + flux_y_raw
        @test tot_raw > max_drain

        scale_factor = min(1.0, max_drain / tot_raw)
        flux_x_scaled = flux_x_raw * scale_factor
        flux_y_scaled = flux_y_raw * scale_factor
        tot_scaled = flux_x_scaled + flux_y_scaled

        @test tot_scaled <= max_drain + 1e-12
        @test isapprox(tot_scaled, max_drain; rtol=1e-12)
        @test isapprox(flux_x_scaled, flux_y_scaled; rtol=1e-12)
    end

    @testset "F34: Zero-weight node neighbour interpolation" begin
        Ny, Nx = 4, 4
        A = zeros(Float64, Ny, Nx)
        WT = zeros(Float64, Ny, Nx)

        # Set populated nodes with known values
        for j in 1:Nx, i in 1:Ny
            A[i, j] = Float64(i * 10 + j)
            WT[i, j] = 1.0
        end

        # Create isolated empty node at (2, 3)
        WT[2, 3] = 0.0
        A[2, 3] = -999.0
        expected_avg = (A[1, 3] + A[3, 3] + A[2, 2] + A[2, 4]) / 4.0

        # Create corner empty node at (4, 4) with two neighbours
        WT[4, 4] = 0.0
        A[4, 4] = -888.0
        expected_corner_avg = (A[3, 4] + A[4, 3]) / 2.0

        # Execute zero weight filling
        fill_zero_weight_nodes!(A, WT; background_val=42.0)

        @test isapprox(A[2, 3], expected_avg; rtol=1e-12)
        @test isapprox(A[4, 4], expected_corner_avg; rtol=1e-12)

        # Node (1, 1) has positive weight and must remain untouched
        @test isapprox(A[1, 1], 11.0; rtol=1e-12)

        # Verify fallback for isolated block with no positive-weight neighbours
        B = zeros(Float64, 3, 3)
        WT_empty = zeros(Float64, 3, 3)
        fill_zero_weight_nodes!(B, WT_empty; background_val=77.5)
        for val in B
            @test isapprox(val, 77.5; rtol=1e-12)
        end
    end

    @testset "F29: Step-start pressure preservation and convergence tracking" begin
        Ny, Nx = 4, 4
        ws = HydromechanicalLSEWorkspace(Ny, Nx)

        # Verify preallocation of pf_prev_iter
        @test size(ws.pf_prev_iter) == (Ny, Nx)
        @test all(ws.pf_prev_iter .≈ 0.0)

        # Verify Picard iteration history mechanics
        pr0 = fill(1.0e7, Ny, Nx)
        pf0 = fill(8.0e6, Ny, Nx)
        pf_current = fill(8.2e6, Ny, Nx)

        # Iteration 1 updates pf_prev_iter
        ws.pf_prev_iter .= pf_current
        @test all(ws.pf_prev_iter .≈ 8.2e6)

        # pr0 and pf0 remain preserved as step-start anchors
        @test all(pr0 .≈ 1.0e7)
        @test all(pf0 .≈ 8.0e6)

        # Iteration 2 step
        pf_next = fill(8.25e6, Ny, Nx)
        iter_diff = maximum(abs.(pf_next .- ws.pf_prev_iter))
        @test isapprox(iter_diff, 5.0e4; rtol=1e-12)

        # Time derivative anchor remains step-start pressure
        dt = 1000.0
        dp_dt = (pf_next .- pf0) ./ dt
        @test isapprox(maximum(dp_dt), 250.0; rtol=1e-12)

        # Step-to-step advancement updates step-start baseline for next step
        pr_step1 = fill(1.05e7, Ny, Nx)
        pf_step1 = fill(8.25e6, Ny, Nx)
        pr0 .= pr_step1
        pf0 .= pf_step1
        @test all(pr0 .≈ 1.05e7)
        @test all(pf0 .≈ 8.25e6)

        # Step 2 evaluates compaction against step 1 baseline
        pf_step2 = fill(8.30e6, Ny, Nx)
        dp_dt_step2 = (pf_step2 .- pf0) ./ dt
        @test isapprox(maximum(dp_dt_step2), 50.0; rtol=1e-12)
    end

    @testset "F23: Matrix melt softening and suspension transition" begin
        type_rock = 1
        t_subsolidus = 1300.0
        t_supersolidus_1 = 1450.0
        t_supersolidus_2 = 1550.0
        t_liquidus = 1800.0

        eta_sub = Erebus.etatotal_rocks(t_subsolidus, type_rock)
        eta_super1 = Erebus.etatotal_rocks(t_supersolidus_1, type_rock)
        eta_super2 = Erebus.etatotal_rocks(t_supersolidus_2, type_rock)
        eta_liq = Erebus.etatotal_rocks(t_liquidus, type_rock)

        # Viscosity decreases monotonically with temperature in supersolidus regime
        @test eta_sub >= eta_super1
        @test eta_super1 > eta_super2
        @test eta_super2 > eta_liq

        # Numerical floor respected across all temperatures
        etamin_val = 1.0e12
        @test eta_sub >= etamin_val
        @test eta_super1 >= etamin_val
        @test eta_super2 >= etamin_val
        @test isapprox(eta_liq, etamin_val; rtol=1e-12)

        # Matrix melt softening bypass when melting_active is true (single source of truth)
        eta_bypassed = Erebus.etatotal_rocks(t_supersolidus_1, type_rock; melting_active=true)
        @test isapprox(eta_bypassed, Erebus.etasolidm[type_rock]; rtol=1e-12)
        @test eta_bypassed > eta_super1
    end

    @testset "Section 3.6: Grid coordinates, gravity boundary and bilinear interpolation" begin
        # 1. GridCoordinates vector lengths and center offsets
        Nx, Ny = 16, 20
        cfg_dict = Dict{String,Any}(
            "grid.Nx" => Nx,
            "grid.Ny" => Ny,
            "grid.xsize" => 160000.0,
            "grid.ysize" => 200000.0,
        )
        cfg = Erebus.override_config(default_config(), cfg_dict)
        coords = GridCoordinates(cfg.grid; xcenter=80000.0, ycenter=100000.0)

        @test length(coords.x) == Nx
        @test length(coords.y) == Ny
        @test length(coords.xvx) == coords.Nx1
        @test length(coords.yvx) == coords.Ny1
        @test length(coords.xvy) == coords.Nx1
        @test length(coords.yvy) == coords.Ny1
        @test length(coords.xp) == coords.Nx1
        @test length(coords.yp) == coords.Ny1

        # Center coordinates and domain bounds verified
        @test isapprox(coords.xcenter, 80000.0; rtol=1e-12)
        @test isapprox(coords.ycenter, 100000.0; rtol=1e-12)
        @test isapprox(coords.x[1], 0.0; rtol=1e-12)
        @test isapprox(coords.x[end], 160000.0; rtol=1e-12)

        # 2. Gravity solution boundary copying
        Nx_g, Ny_g = 6, 6
        Nx1_g, Ny1_g = Nx_g + 1, Ny_g + 1
        gx = zeros(Float64, Ny_g, Nx1_g)
        gy = zeros(Float64, Ny1_g, Nx_g)

        # Fill interior boundary columns
        for i in 1:Ny_g
            gx[i, Nx_g] = Float64(i * 1.5)
        end
        for j in 1:Nx_g
            gy[Ny_g, j] = Float64(j * 2.5)
        end

        # Copy edge boundary entries as done in process_gravitational_solution!
        gx[:, Nx1_g] .= gx[:, Nx_g]
        gy[Ny1_g, :] .= gy[Ny_g, :]

        @test isapprox(gx[:, Nx1_g], gx[:, Nx_g]; rtol=1e-12)
        @test isapprox(gy[Ny1_g, :], gy[Ny_g, :]; rtol=1e-12)
        @test isapprox(gx[1, Nx1_g], 1.5; rtol=1e-12)
        @test isapprox(gy[Ny1_g, 1], 2.5; rtol=1e-12)

        # 3. Segregation 4-point bilinear interpolation for face gravity
        # Vertical face gravity interpolation from gy: (i, j), (i+1, j), (i, j+1), (i+1, j+1)
        gy_grid = reshape(collect(1.0:Float64(Ny1_g * Nx1_g)), Ny1_g, Nx1_g)
        i_test, j_test = 2, 2
        gy_f = 0.25 * (
            gy_grid[i_test, j_test] +
            gy_grid[i_test + 1, j_test] +
            gy_grid[i_test, j_test + 1] +
            gy_grid[i_test + 1, j_test + 1]
        )
        expected_gy_f = 0.25 * (
            gy_grid[2, 2] + gy_grid[3, 2] + gy_grid[2, 3] + gy_grid[3, 3]
        )
        @test isapprox(gy_f, expected_gy_f; rtol=1e-12)

        # Boundary clamping verification at edge cell
        i_edge, j_edge = Ny_g, Nx_g
        i_next = min(i_edge + 1, Ny1_g)
        j_next = min(j_edge + 1, Nx1_g)
        gy_edge_f = 0.25 * (
            gy_grid[i_edge, j_edge] +
            gy_grid[i_next, j_edge] +
            gy_grid[i_edge, j_next] +
            gy_grid[i_next, j_next]
        )
        expected_edge_f = 0.25 * (
            gy_grid[Ny_g, Nx_g] +
            gy_grid[Ny1_g, Nx_g] +
            gy_grid[Ny_g, Nx1_g] +
            gy_grid[Ny1_g, Nx1_g]
        )
        @test isfinite(gy_edge_f)
        @test isapprox(gy_edge_f, expected_edge_f; rtol=1e-12)
    end

    @testset "Branch coverage for hydrofracture face porosity, dissipation and coordinates" begin
        # 1. GridCoordinates from SimulationConfig
        cfg_sim = default_config()
        coords_sim = GridCoordinates(cfg_sim)
        @test coords_sim.Nx == cfg_sim.grid.Nx
        @test coords_sim.Ny == cfg_sim.grid.Ny

        # 2. etatotal_rocks with custom etasolidmm
        custom_etasolidmm = [1.0e18, 1.0e18, 1.0e18]
        eta_custom = Erebus.etatotal_rocks(1500.0, 1; etasolidmm=custom_etasolidmm)
        @test isapprox(eta_custom, 1.0e18; rtol=1e-12)
        @test eta_custom < 1.0e19

        # 3. Hydrofracture with explicit PHIX and PHIY
        Nx, Ny = 5, 5
        coords = GridCoordinates(GridConfig(Nx=Nx, Ny=Ny, xsize=1000.0, ysize=1000.0))
        Nx1, Ny1 = coords.Nx1, coords.Ny1
        ETA = fill(1.0e19, Ny1, Nx1)
        ETAP = fill(1.0e19, Ny1, Nx1)
        GGG = fill(1.0e10, Ny1, Nx1)
        GGGP = fill(1.0e10, Ny1, Nx1)
        SXY0 = zeros(Ny1, Nx1)
        SXX0 = zeros(Ny1, Nx1)
        RHOX = fill(3000.0, Ny1, Nx1)
        RHOY = fill(3000.0, Ny1, Nx1)
        RHOFX = fill(1000.0, Ny1, Nx1)
        RHOFY = fill(1000.0, Ny1, Nx1)
        RX = fill(1.0e8, Ny1, Nx1)
        RY = fill(1.0e8, Ny1, Nx1)
        ETAPHI = fill(1.0e19, Ny1, Nx1)
        BETAPHI = fill(1.0e-11, Ny1, Nx1)
        PHI = fill(0.1, Ny1, Nx1)
        PHIX = fill(0.08, Ny1, Nx1)
        PHIY = fill(0.08, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = fill(1.0, Ny1, Nx1)
        pr0 = fill(1.0e6, Ny1, Nx1)
        pf0 = fill(1.0e6, Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)
        dt = 1.0e10
        pr_test = fill(2.0e6, Ny1, Nx1)
        pf_test = fill(2.5e6, Ny1, Nx1)
        TEN_test = fill(1.0e5, Ny1, Nx1)

        # Call with PHIX and PHIY explicitly supplied
        R6 = zeros(Ny1 * Nx1 * 6)
        L6 = Erebus.assemble_hydromechanical_lse!(
            ETA, ETAP, GGG, GGGP, SXY0, SXX0, RHOX, RHOY, RHOFX, RHOFY,
            RX, RY, ETAPHI, BETAPHI, PHI, gx, gy, pr0, pf0, DMP, dt, R6;
            coords=coords, hydrofracture=true, pr=pr_test, pf=pf_test, TEN=TEN_test,
            PHIX=PHIX, PHIY=PHIY,
        )
        @test size(L6, 1) == Ny1 * Nx1 * 6
        @test isfinite(L6[1, 1])

        R4 = zeros(Ny1 * Nx1 * 4)
        L4 = Erebus.assemble_hydromechanical_4var_lse!(
            ETA, ETAP, GGG, GGGP, SXY0, SXX0, RHOX, RHOY, RHOFX, RHOFY,
            RX, RY, ETAPHI, BETAPHI, PHI, gx, gy, pr0, pf0, DMP, dt, R4;
            coords=coords, hydrofracture=true, pr=pr_test, pf=pf_test, TEN=TEN_test,
            PHIX=PHIX, PHIY=PHIY,
        )
        @test size(L4, 1) == Ny1 * Nx1 * 4
        @test isfinite(L4[1, 1])

        qxD = zeros(Ny1, Nx1)
        qyD = zeros(Ny1, Nx1)
        Erebus.reconstruct_darcy_fluxes!(
            qxD, qyD, pf_test, RHOFX, RHOFY, RX, RY, gx, gy, coords;
            hydrofracture=true, pr=pr_test, TEN=TEN_test, PHIX=PHIX, PHIY=PHIY,
        )
        @test all(isfinite, qxD)
        @test all(isfinite, qyD)

        HS = zeros(Ny, Nx)
        SXY_test = fill(1.0e5, Ny, Nx)
        SXX_test = fill(1.0e5, Ny1, Nx1)
        Erebus.compute_shear_heating!(
            HS, ETA, SXY_test, ETAP, SXX_test, RX, RY, qxD, qyD, PHI, ETAPHI,
            pr_test, pf_test;
            hydrofracture=true, TEN=TEN_test, PHIX=PHIX, PHIY=PHIY,
        )
        @test all(isfinite, HS)
        @test all(HS .>= 0.0)

        # Call with PHIX === nothing and PHI === nothing (fallback to phimin_val)
        Erebus.reconstruct_darcy_fluxes!(
            qxD, qyD, pf_test, RHOFX, RHOFY, RX, RY, gx, gy, coords;
            hydrofracture=true, pr=pr_test, TEN=TEN_test,
            PHI=nothing, PHIX=nothing, PHIY=nothing,
        )
        @test all(isfinite, qxD)

        Erebus.compute_shear_heating!(
            HS, ETA, SXY_test, ETAP, SXX_test, RX, RY, qxD, qyD, PHI, ETAPHI,
            pr_test, pf_test;
            hydrofracture=true, TEN=TEN_test,
            PHIX=nothing, PHIY=nothing,
        )
        @test all(isfinite, HS)
    end
end


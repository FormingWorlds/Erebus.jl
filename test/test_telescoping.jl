using Test
using LinearAlgebra
using Erebus
using Erebus.Config
using Erebus.Particles
using Erebus.Numerics
using JLD2

@testset "Telescoping Domain Mechanics & Conservation Invariants" begin
    # ---------------------------------------------------------------------
    # 1. TelescopingConfig Schema, Defaults, and Bounds Validation
    # ---------------------------------------------------------------------
    @testset "TelescopingConfig Schema & Defaults" begin
        cfg_def = TelescopingConfig()
        @test cfg_def.active == false
        @test isapprox(cfg_def.r_threshold_fraction, 0.70; rtol=1e-12)
        @test cfg_def.max_telescope_levels == 10
        @test isapprox(cfg_def.target_radius, 1_737_000.0; rtol=1e-12)
        @test cfg_def.buffer_markers_per_cell == 4

        # Custom config
        cfg_custom = TelescopingConfig(;
            active=true,
            r_threshold_fraction=0.75,
            max_telescope_levels=4,
            target_radius=1_500_000.0,
            buffer_markers_per_cell=4,
        )
        @test cfg_custom.active == true
        @test isapprox(cfg_custom.r_threshold_fraction, 0.75; rtol=1e-12)
        @test cfg_custom.max_telescope_levels == 4

        # Bounds validation in SimulationConfig
        sim_cfg_bad_frac = SimulationConfig(;
            telescoping=TelescopingConfig(; active=true, r_threshold_fraction=1.5)
        )
        @test_throws ArgumentError validate_config(sim_cfg_bad_frac)

        sim_cfg_bad_frac_low = SimulationConfig(;
            telescoping=TelescopingConfig(; active=true, r_threshold_fraction=0.0)
        )
        @test_throws ArgumentError validate_config(sim_cfg_bad_frac_low)

        sim_cfg_bad_levels = SimulationConfig(;
            telescoping=TelescopingConfig(; active=true, max_telescope_levels=0)
        )
        @test_throws ArgumentError validate_config(sim_cfg_bad_levels)

        sim_cfg_bad_radius = SimulationConfig(;
            telescoping=TelescopingConfig(; active=true, target_radius=-1000.0)
        )
        @test_throws ArgumentError validate_config(sim_cfg_bad_radius)

        sim_cfg_bad_buffer = SimulationConfig(;
            telescoping=TelescopingConfig(; active=true, buffer_markers_per_cell=0)
        )
        @test_throws ArgumentError validate_config(sim_cfg_bad_buffer)

        # Even grid dimensions rejected when telescoping is active
        sim_cfg_even_nx = SimulationConfig(;
            grid=GridConfig(Nx=20, Ny=17), telescoping=TelescopingConfig(; active=true)
        )
        @test_throws ArgumentError validate_config(sim_cfg_even_nx)

        sim_cfg_even_ny = SimulationConfig(;
            grid=GridConfig(Nx=17, Ny=20), telescoping=TelescopingConfig(; active=true)
        )
        @test_throws ArgumentError validate_config(sim_cfg_even_ny)
    end

    # ---------------------------------------------------------------------
    # 2. Trigger Evaluation (should_telescope_domain)
    # ---------------------------------------------------------------------
    @testset "Telescoping Trigger Logic" begin
        coords = GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        cfg_tele = TelescopingConfig(; active=true, r_threshold_fraction=0.70)

        # Half domain is 70_000.0 m. Threshold is 0.70 * 70_000 = 49_000.0 m
        @test !should_telescope_domain(40_000.0, coords, cfg_tele; level=0)
        @test !should_telescope_domain(49_000.0, coords, cfg_tele; level=0)
        @test should_telescope_domain(49_001.0, coords, cfg_tele; level=0)
        @test should_telescope_domain(60_000.0, coords, cfg_tele; level=0)

        # Inactive configuration
        cfg_inactive = TelescopingConfig(; active=false, r_threshold_fraction=0.70)
        @test !should_telescope_domain(60_000.0, coords, cfg_inactive; level=0)

        # Maximum level reached
        cfg_max = TelescopingConfig(;
            active=true, r_threshold_fraction=0.70, max_telescope_levels=2
        )
        @test should_telescope_domain(60_000.0, coords, cfg_max; level=1)
        @test !should_telescope_domain(60_000.0, coords, cfg_max; level=2)

        # Domain errors
        @test_throws DomainError should_telescope_domain(-100.0, coords, cfg_tele; level=0)
        @test_throws DomainError should_telescope_domain(NaN, coords, cfg_tele; level=0)
        @test_throws DomainError should_telescope_domain(Inf, coords, cfg_tele; level=0)
        @test_throws DomainError should_telescope_domain(-Inf, coords, cfg_tele; level=0)
    end

    # ---------------------------------------------------------------------
    # 3. Coordinate Doubling & Cell Size Invariant (compute_telescoped_coordinates)
    # ---------------------------------------------------------------------
    @testset "Coordinate Doubling & Spacing Invariants" begin
        c17 = GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        c33 = compute_telescoped_coordinates(c17)

        # Node counts: N_new = 2(N_old - 1) + 1
        @test c33.Nx == 33
        @test c33.Ny == 33
        @test c33.Nx1 == 34
        @test c33.Ny1 == 34

        # Physical dimensions: doubled
        @test isapprox(c33.xsize, 280_000.0; rtol=1e-12)
        @test isapprox(c33.ysize, 280_000.0; rtol=1e-12)

        # Spacing invariant: dx and dy remain exactly constant
        @test isapprox(c33.dx, c17.dx; rtol=1e-12)
        @test isapprox(c33.dy, c17.dy; rtol=1e-12)

        # Center coordinates: doubled
        @test isapprox(c33.xcenter, 2.0 * c17.xcenter; rtol=1e-12)
        @test isapprox(c33.ycenter, 2.0 * c17.ycenter; rtol=1e-12)

        # Center shift invariant
        shift_x = c33.xcenter - c17.xcenter
        shift_y = c33.ycenter - c17.ycenter
        @test isapprox(shift_x, c17.xsize / 2.0; rtol=1e-12)
        @test isapprox(shift_y, c17.ysize / 2.0; rtol=1e-12)

        # Multi-level telescoping: 17 -> 33 -> 65
        c65 = compute_telescoped_coordinates(c33)
        @test c65.Nx == 65
        @test c65.Ny == 65
        @test isapprox(c65.xsize, 560_000.0; rtol=1e-12)
        @test isapprox(c65.dx, c17.dx; rtol=1e-12)

        # Even grid dimensions rejected for staggered grid alignment
        c_even = GridCoordinates(20, 20; xsize=140_000.0, ysize=140_000.0)
        @test_throws ArgumentError compute_telescoped_coordinates(c_even)
    end

    # ---------------------------------------------------------------------
    # 4. Exact Staggered Grid Remapping (remap_staggered_grid_array)
    # ---------------------------------------------------------------------
    @testset "Staggered Grid Remapping & Exact Centered Alignment" begin
        # 1. Basic node array (17x17 -> 33x33)
        arr_old = rand(17, 17)
        arr_new = remap_staggered_grid_array(arr_old, (33, 33); background_val=-999.0)
        @test size(arr_new) == (33, 33)

        # Offset: ioff = (33 - 17) / 2 = 8, joff = 8
        ioff = 8
        joff = 8
        @test arr_new[(ioff + 1):(ioff + 17), (joff + 1):(joff + 17)] == arr_old
        # Verify outer border cells have background value
        @test all(isapprox.(arr_new[1:ioff, :], -999.0))
        @test all(isapprox.(arr_new[(ioff + 18):33, :], -999.0))
        @test all(isapprox.(arr_new[:, 1:joff], -999.0))
        @test all(isapprox.(arr_new[:, (joff + 18):33], -999.0))

        # Asymmetric / odd-difference dimensions rejected
        @test_throws ArgumentError remap_staggered_grid_array(arr_old, (32, 33))
        @test_throws ArgumentError remap_staggered_grid_array(arr_old, (33, 32))

        # 2. P-node array (18x18 -> 34x34)
        p_old = fill(1500.0, 18, 18)
        p_old[9:10, 9:10] .= 2000.0
        p_new = remap_staggered_grid_array(p_old, (34, 34); background_val=200.0)
        @test size(p_new) == (34, 34)
        ioff_p = (34 - 18) ÷ 2
        joff_p = (34 - 18) ÷ 2
        @test p_new[(ioff_p + 1):(ioff_p + 18), (joff_p + 1):(joff_p + 18)] == p_old
        @test isapprox(p_new[1, 1], 200.0)
        @test isapprox(p_new[34, 34], 200.0)
    end

    # ---------------------------------------------------------------------
    # 5. Marker Radial Distance & Center Invariance (telescope_marker_arrays!)
    # ---------------------------------------------------------------------
    @testset "Marker Radial Distance & Center Invariance" begin
        c17 = GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        c33 = compute_telescoped_coordinates(c17)

        N_m = 50
        # Place markers in planetesimal sphere around old center
        r_seed = 30_000.0
        theta = range(0, 2pi, length=N_m)
        xm = [c17.xcenter + 0.8 * r_seed * cos(th) for th in theta]
        ym = [c17.ycenter + 0.8 * r_seed * sin(th) for th in theta]
        tm = fill(2, N_m) # solid rock
        tkm = fill(1200.0, N_m)
        phim = fill(0.15, N_m)
        sxxm = zeros(N_m)
        sxym = zeros(N_m)
        etavpm = fill(1.0e20, N_m)
        phinewm = fill(0.15, N_m)
        pfm0 = zeros(N_m)
        XWsolidm = fill(0.1, N_m)
        XWsolidm0 = fill(0.1, N_m)
        Fm = zeros(N_m)
        rhototalm = fill(3200.0, N_m)
        rhocptotalm = fill(3.2e6, N_m)
        etatotalm = fill(1.0e20, N_m)
        hrtotalm = fill(1.0e-7, N_m)
        ktotalm = fill(3.0, N_m)
        inv_gggtotalm = fill(1.0e-10, N_m)
        fricttotalm = fill(0.6, N_m)
        cohestotalm = fill(1.0e7, N_m)
        tenstotalm = fill(1.0e7, N_m)
        rhofluidcur = fill(1000.0, N_m)
        alphasolidcur = fill(3.0e-5, N_m)
        alphafluidcur = fill(2.0e-4, N_m)
        tkm_rhocptotalm = fill(6.0e8, N_m)
        etafluidcur_inv_kphim = fill(1.0e10, N_m)
        Xfem = fill(0.05, N_m)
        Xfem0 = fill(0.05, N_m)
        Xfe_bulk = fill(0.20, N_m)
        XH2Om = fill(2.0, N_m)
        XCm = fill(500.0, N_m)
        XNm = fill(50.0, N_m)
        XSm = fill(1000.0, N_m)
        t_acc = fill(1.0e12, N_m)

        # Compute pre-telescoping radial distances
        r_old = [sqrt((xm[m] - c17.xcenter)^2 + (ym[m] - c17.ycenter)^2) for m in 1:N_m]

        new_marknum = telescope_marker_arrays!(
            xm,
            ym,
            tm,
            tkm,
            sxxm,
            sxym,
            etavpm,
            phim,
            phinewm,
            pfm0,
            XWsolidm,
            XWsolidm0,
            Fm,
            rhototalm,
            rhocptotalm,
            etatotalm,
            hrtotalm,
            ktotalm,
            inv_gggtotalm,
            fricttotalm,
            cohestotalm,
            tenstotalm,
            rhofluidcur,
            alphasolidcur,
            alphafluidcur,
            tkm_rhocptotalm,
            etafluidcur_inv_kphim;
            old_coords=c17,
            new_coords=c33,
            T_ambient=250.0,
            phi_ambient=0.35,
            Xfem=Xfem,
            Xfem0=Xfem0,
            Xfe_bulk=Xfe_bulk,
            XH2Om=XH2Om,
            XCm=XCm,
            XNm=XNm,
            XSm=XSm,
            t_accreted=t_acc,
        )

        @test new_marknum > N_m
        @test length(xm) == new_marknum
        @test length(t_acc) == new_marknum
        @test length(Xfem) == new_marknum

        # INVARIANT: Radial distance from planetesimal center is preserved to machine precision
        r_new = [sqrt((xm[m] - c33.xcenter)^2 + (ym[m] - c33.ycenter)^2) for m in 1:N_m]
        for m in 1:N_m
            @test isapprox(r_new[m], r_old[m]; rtol=1e-12)
        end

        # INVARIANT: New outer markers are sticky air (tm = 3) at ambient temperature
        for m in (N_m + 1):new_marknum
            @test tm[m] == 3
            @test isapprox(tkm[m], 250.0)
            @test isapprox(phim[m], 0.35)
            @test isapprox(Xfe_bulk[m], 0.0; atol=1e-15)
            @test isapprox(XH2Om[m], 0.0; atol=1e-15)
        end

        # Verify exact marker injection for arbitrary and non-square buffer marker counts (no silent shortfall)
        n_outer_cells = ((c33.Nx - 1) * (c33.Ny - 1)) - ((c17.Nx - 1) * (c17.Ny - 1)) # (32*32) - (16*16) = 768
        for n_buf in [1, 4, 5, 6, 7, 8, 9, 10]
            xm_t = [c17.xcenter]
            ym_t = [c17.ycenter]
            tm_t = [2]
            tkm_t = [300.0]
            sxxm_t = [0.0]
            sxym_t = [0.0]
            etavpm_t = [1.0e18]
            phim_t = [0.1]
            phinewm_t = [0.1]
            pfm0_t = [0.0]
            XWsolidm_t = [0.0]
            XWsolidm0_t = [0.0]
            Fm_t = [0.0]
            rhototalm_t = [3000.0]
            rhocptotalm_t = [3.0e6]
            etatotalm_t = [1.0e18]
            hrtotalm_t = [0.0]
            ktotalm_t = [3.0]
            inv_gggtotalm_t = [1.0e-10]
            fricttotalm_t = [0.6]
            cohestotalm_t = [1.0e7]
            tenstotalm_t = [6.0e6]
            rhofluidcur_t = [1000.0]
            alphasolidcur_t = [3.0e-5]
            alphafluidcur_t = [1.0e-4]
            tkm_rhocptotalm_t = [9.0e8]
            etafluidcur_inv_kphim_t = [1.0e10]

            tot_m = telescope_marker_arrays!(
                xm_t,
                ym_t,
                tm_t,
                tkm_t,
                sxxm_t,
                sxym_t,
                etavpm_t,
                phim_t,
                phinewm_t,
                pfm0_t,
                XWsolidm_t,
                XWsolidm0_t,
                Fm_t,
                rhototalm_t,
                rhocptotalm_t,
                etatotalm_t,
                hrtotalm_t,
                ktotalm_t,
                inv_gggtotalm_t,
                fricttotalm_t,
                cohestotalm_t,
                tenstotalm_t,
                rhofluidcur_t,
                alphasolidcur_t,
                alphafluidcur_t,
                tkm_rhocptotalm_t,
                etafluidcur_inv_kphim_t;
                old_coords=c17,
                new_coords=c33,
                buffer_markers_per_cell=n_buf,
            )
            @test tot_m == 1 + n_outer_cells * n_buf
            @test length(xm_t) == tot_m
        end
    end

    # ---------------------------------------------------------------------
    # 6. Physical Conservation Invariants
    # ---------------------------------------------------------------------
    @testset "Planetesimal Mass & Thermal Energy Conservation Invariants" begin
        c17 = GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        c33 = compute_telescoped_coordinates(c17)

        N_m = 100
        r_seed = 35_000.0
        xm = [c17.xcenter + r_seed * (rand() - 0.5) for _ in 1:N_m]
        ym = [c17.ycenter + r_seed * (rand() - 0.5) for _ in 1:N_m]
        tm = fill(2, N_m)
        tkm = [1000.0 + 500.0 * rand() for _ in 1:N_m]
        rhototalm = fill(3000.0, N_m)
        rhocptotalm = fill(3.0e6, N_m)
        phim = fill(0.2, N_m)
        Xfe_bulk = fill(0.15, N_m)
        XH2Om = fill(1.5, N_m)
        t_acc = fill(2.0e12, N_m)

        # Baseline physical integrals for solid planetesimal
        solid_count_old = count(m -> tm[m] in (1, 2), 1:N_m)
        thermal_energy_old = sum(rhocptotalm[m] * tkm[m] for m in 1:N_m if tm[m] in (1, 2))
        iron_inventory_old = sum(Xfe_bulk[m] for m in 1:N_m if tm[m] in (1, 2))
        water_inventory_old = sum(XH2Om[m] for m in 1:N_m if tm[m] in (1, 2))

        # Dummy arrays
        sxxm = zeros(N_m)
        sxym = zeros(N_m)
        etavpm = fill(1e20, N_m)
        phinewm = fill(0.2, N_m)
        pfm0 = zeros(N_m)
        XWsolidm = zeros(N_m)
        XWsolidm0 = zeros(N_m)
        Fm = zeros(N_m)
        etatotalm = fill(1e20, N_m)
        hrtotalm = zeros(N_m)
        ktotalm = fill(2.0, N_m)
        inv_gggtotalm = fill(1e-10, N_m)
        fricttotalm = fill(0.6, N_m)
        cohestotalm = fill(1e7, N_m)
        tenstotalm = fill(1e7, N_m)
        rhofluidcur = fill(1000.0, N_m)
        alphasolidcur = fill(3e-5, N_m)
        alphafluidcur = fill(2e-4, N_m)
        tkm_rhocptotalm = fill(6e8, N_m)
        etafluidcur_inv_kphim = fill(1e10, N_m)
        Xfem = fill(0.05, N_m)
        Xfem0 = fill(0.05, N_m)
        XCm = fill(500.0, N_m)
        XNm = fill(50.0, N_m)
        XSm = fill(1000.0, N_m)

        new_marknum = telescope_marker_arrays!(
            xm,
            ym,
            tm,
            tkm,
            sxxm,
            sxym,
            etavpm,
            phim,
            phinewm,
            pfm0,
            XWsolidm,
            XWsolidm0,
            Fm,
            rhototalm,
            rhocptotalm,
            etatotalm,
            hrtotalm,
            ktotalm,
            inv_gggtotalm,
            fricttotalm,
            cohestotalm,
            tenstotalm,
            rhofluidcur,
            alphasolidcur,
            alphafluidcur,
            tkm_rhocptotalm,
            etafluidcur_inv_kphim;
            old_coords=c17,
            new_coords=c33,
            T_ambient=250.0,
            phi_ambient=0.35,
            Xfem=Xfem,
            Xfem0=Xfem0,
            Xfe_bulk=Xfe_bulk,
            XH2Om=XH2Om,
            XCm=XCm,
            XNm=XNm,
            XSm=XSm,
            t_accreted=t_acc,
        )

        solid_count_new = count(m -> tm[m] in (1, 2), 1:new_marknum)
        thermal_energy_new = sum(
            rhocptotalm[m] * tkm[m] for m in 1:new_marknum if tm[m] in (1, 2)
        )
        iron_inventory_new = sum(Xfe_bulk[m] for m in 1:new_marknum if tm[m] in (1, 2))
        water_inventory_new = sum(XH2Om[m] for m in 1:new_marknum if tm[m] in (1, 2))

        # CONSERVATION LAWS: exactly zero change in solid planetesimal contents
        @test solid_count_new == solid_count_old
        @test isapprox(thermal_energy_new, thermal_energy_old; rtol=1e-14)
        @test isapprox(iron_inventory_new, iron_inventory_old; rtol=1e-14)
        @test isapprox(water_inventory_new, water_inventory_old; rtol=1e-14)

        # STRUCTURAL BUFFER INVARIANTS:
        # 1. Total marker count matches old plus outer cell injection
        n_outer_cells = ((c33.Nx - 1) * (c33.Ny - 1)) - ((c17.Nx - 1) * (c17.Ny - 1))
        @test new_marknum == N_m + n_outer_cells * 4

        # 2. Every newly added buffer marker lies inside [0, xsize_new] x [0, ysize_new]
        #    and outside the shifted inner rectangle [shift_x, shift_x + xsize_old]
        shift_x = c17.xsize / 2.0
        shift_y = c17.ysize / 2.0
        for m in (N_m + 1):new_marknum
            @test 0.0 <= xm[m] <= c33.xsize
            @test 0.0 <= ym[m] <= c33.ysize
            in_inner =
                (shift_x < xm[m] < shift_x + c17.xsize) &&
                (shift_y < ym[m] < shift_y + c17.ysize)
            @test !in_inner
            @test tm[m] == 3
            @test isapprox(Xfe_bulk[m], 0.0; atol=1e-15)
            @test isapprox(XH2Om[m], 0.0; atol=1e-15)
            @test isapprox(XCm[m], 0.0; atol=1e-15)
            @test isapprox(XNm[m], 0.0; atol=1e-15)
            @test isapprox(XSm[m], 0.0; atol=1e-15)
        end
    end

    # ---------------------------------------------------------------------
    # 7. Gravitational Poisson Re-factorization on Telescoped Grid
    # ---------------------------------------------------------------------
    @testset "Gravitational Poisson Re-factorization" begin
        c17 = GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        c33 = compute_telescoped_coordinates(c17)

        # Setup and factorize for 33x33 telescoped grid
        RP, SP = setup_gravitational_lse(c33)
        RHO = fill(3000.0, c33.Ny1, c33.Nx1)
        LP = assemble_gravitational_lse!(RHO, RP; coords=c33)
        F_grav = lu(LP.cscmatrix)

        FI = zeros(c33.Ny1, c33.Nx1)
        gx = zeros(c33.Ny1, c33.Nx1)
        gy = zeros(c33.Ny1, c33.Nx1)

        # Solve Poisson equation with new factorized operator
        Erebus.compute_gravity_solution!(SP, RP, RHO, FI, gx, gy; coords=c33)

        @test all(isfinite, FI)
        @test all(isfinite, gx)
        @test all(isfinite, gy)
        # Potential must be negative in interior (attractive gravity)
        @test minimum(FI) < 0.0
        # Dirichlet boundary: FI must be zero on the outer boundary
        @test isapprox(FI[1, 1], 0.0; atol=1e-15)
        @test isapprox(FI[c33.Ny1, c33.Nx1], 0.0; atol=1e-15)
    end

    # ---------------------------------------------------------------------
    # 8. Multi-Step Continuous Accretion Simulation with Domain Doubling
    # ---------------------------------------------------------------------
    @testset "Multi-Step Continuous Accretion & Domain Doubling Integration" begin
        output_dir = mktempdir()
        try
            sim_cfg = SimulationConfig(;
                grid=GridConfig(Nx=17, Ny=17, xsize=140_000.0, ysize=140_000.0),
                geometry=GeometryConfig(rplanet=45_000.0, rcrust=45_000.0),
                time=TimeConfig(
                    dt_initial=100.0,
                    dt_longest=100.0,
                    dtcoefdn=0.5,
                    dtcoefup=1.0,
                    dtstep=10,
                    dxymax=0.05,
                    vpratio=0.33,
                    DTmax=20.0,
                    yearlength=3.15e7,
                    start_time=0.0,
                    endtime=1.0e6,
                    start_step=1,
                    n_steps=2,
                ),
                solver=SolverConfig(titermax=2, nplast=1),
                output=OutputConfig(output_dir=output_dir, savematstep=1, visstep=0),
                accretion=AccretionConfig(
                    active=true,
                    mode=:constant_rate,
                    R_initial=45_000.0,
                    M_initial=1.0e21,
                    R_target=60_000.0,
                    M_target=2.0e21,
                    dR_dt_constant=100.0,
                    dM_dt_constant=1.0e15,
                ),
                telescoping=TelescopingConfig(
                    active=true, r_threshold_fraction=0.70, max_telescope_levels=2
                ),
            )

            Erebus.simulation_loop(sim_cfg; output_path=output_dir)

            # Check that output checkpoints exist
            step1_file = joinpath(output_dir, "output_00001.jld2")
            step2_file = joinpath(output_dir, "output_00002.jld2")
            @test isfile(step1_file)
            @test isfile(step2_file)

            data1 = JLD2.load(step1_file)
            @test data1["telescope_level"] == 1
            @test data1["Nx"] == 33
            @test data1["Ny"] == 33
            @test isapprox(data1["xsize"], 280_000.0; rtol=1e-12)
            @test isapprox(data1["rplanet"], 60_000.0; rtol=1e-6)

            data2 = JLD2.load(step2_file)
            @test data2["telescope_level"] == 1
            @test data2["Nx"] == 33
            @test isapprox(data2["rplanet"], 60_000.0; rtol=1e-6)
            @test all(isfinite, data2["tk2"])
            @test all(isfinite, data2["pr"])

            # Test restart from telescoped checkpoint
            restart_output_dir = mktempdir()
            try
                restart_cfg = SimulationConfig(;
                    grid=GridConfig(Nx=17, Ny=17, xsize=140_000.0, ysize=140_000.0),
                    geometry=GeometryConfig(rplanet=45_000.0, rcrust=45_000.0),
                    time=TimeConfig(
                        dt_initial=100.0,
                        dt_longest=100.0,
                        dtcoefdn=0.5,
                        dtcoefup=1.0,
                        dtstep=10,
                        dxymax=0.05,
                        vpratio=0.33,
                        DTmax=20.0,
                        yearlength=3.15e7,
                        start_time=0.0,
                        endtime=1.0e6,
                        start_step=1,
                        n_steps=2,
                    ),
                    solver=SolverConfig(titermax=2, nplast=1),
                    output=OutputConfig(
                        output_dir=restart_output_dir, savematstep=1, visstep=0
                    ),
                    accretion=AccretionConfig(
                        active=true,
                        mode=:constant_rate,
                        R_initial=45_000.0,
                        M_initial=1.0e21,
                        R_target=60_000.0,
                        M_target=2.0e21,
                        dR_dt_constant=100.0,
                        dM_dt_constant=1.0e15,
                    ),
                    telescoping=TelescopingConfig(
                        active=true, r_threshold_fraction=0.70, max_telescope_levels=2
                    ),
                )
                Erebus.simulation_loop(
                    restart_cfg; output_path=restart_output_dir, restart_from=step1_file
                )
                restarted_file = joinpath(restart_output_dir, "output_00002.jld2")
                @test isfile(restarted_file)
                restarted_data = JLD2.load(restarted_file)
                @test restarted_data["telescope_level"] == 1
                @test restarted_data["Nx"] == 33
            finally
                rm(restart_output_dir; recursive=true, force=true)
            end
        finally
            rm(output_dir; recursive=true, force=true)
        end
    end

    # ---------------------------------------------------------------------
    # 9. HCNSPO Marker Array Replenishment Under Domain Doubling
    # ---------------------------------------------------------------------
    @testset "telescope_marker_arrays!() with HCNSPO marker replenishment" begin
        c17 = GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        c33 = GridCoordinates(33, 33; xsize=280_000.0, ysize=280_000.0)
        N_m = 1000
        xm = rand(N_m) .* 100_000.0 .+ 20_000.0
        ym = rand(N_m) .* 100_000.0 .+ 20_000.0
        tm = fill(1, N_m)
        tkm = fill(300.0, N_m)
        sxxm = zeros(N_m)
        sxym = zeros(N_m)
        etavpm = fill(1.0e22, N_m)
        phim = fill(0.1, N_m)
        phinewm = fill(0.1, N_m)
        pfm0 = zeros(N_m)
        XWsolidm = fill(0.1, N_m)
        XWsolidm0 = fill(0.1, N_m)
        Fm = zeros(N_m)
        rhototalm = fill(3200.0, N_m)
        rhocptotalm = fill(3.2e6, N_m)
        etatotalm = fill(1.0e20, N_m)
        hrtotalm = fill(1.0e-7, N_m)
        ktotalm = fill(3.0, N_m)
        inv_gggtotalm = fill(1.0e-10, N_m)
        fricttotalm = fill(0.6, N_m)
        cohestotalm = fill(1.0e7, N_m)
        tenstotalm = fill(1.0e7, N_m)
        rhofluidcur = fill(1000.0, N_m)
        alphasolidcur = fill(3.0e-5, N_m)
        alphafluidcur = fill(2.0e-4, N_m)
        tkm_rhocptotalm = fill(6.0e8, N_m)
        etafluidcur_inv_kphim = fill(1.0e10, N_m)

        # Setup HCNSPO marker properties
        cfg_vm = VolatileMixtureConfig()
        cfg_refr = RefractoryConfig()
        hcnspo = setup_marker_hcnspo_properties(N_m, cfg_vm, cfg_refr)

        h2o_init = copy(hcnspo.X_ice_H2O_m)

        new_marknum = telescope_marker_arrays!(
            xm,
            ym,
            tm,
            tkm,
            sxxm,
            sxym,
            etavpm,
            phim,
            phinewm,
            pfm0,
            XWsolidm,
            XWsolidm0,
            Fm,
            rhototalm,
            rhocptotalm,
            etatotalm,
            hrtotalm,
            ktotalm,
            inv_gggtotalm,
            fricttotalm,
            cohestotalm,
            tenstotalm,
            rhofluidcur,
            alphasolidcur,
            alphafluidcur,
            tkm_rhocptotalm,
            etafluidcur_inv_kphim;
            old_coords=c17,
            new_coords=c33,
            hcnspo_props=hcnspo,
        )

        @test new_marknum > N_m
        for arr in values(hcnspo)
            @test length(arr) == new_marknum
            # Buffer markers are initialized to zero
            @test all(iszero, arr[(N_m + 1):new_marknum])
        end
        # Existing markers preserve values
        @test isapprox(hcnspo.X_ice_H2O_m[1:N_m], h2o_init; rtol=1e-12)
    end
end

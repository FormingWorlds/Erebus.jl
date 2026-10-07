using Erebus
using Erebus:
    setup_staggered_grid_properties,
    setup_staggered_grid_properties_helpers,
    setup_gravitational_lse,
    setup_hydromechanical_lse,
    setup_thermal_lse,
    setup_interpolated_properties,
    allocate_thread_interpolation_buffers,
    validate_config,
    default_config,
    default_grid_coordinates
using Random
using Test

@testset "Grid Single Source of Truth" begin
    @testset "Compile-time grid constants deleted" begin
        # Verify compiled grid dimension constants do not exist at module level
        @test !isdefined(Erebus, :Nx)
        @test !isdefined(Erebus, :Ny)
        @test !isdefined(Erebus, :Nx1)
        @test !isdefined(Erebus, :Ny1)
        @test !isdefined(Erebus, :xsize)
        @test !isdefined(Erebus, :ysize)
        @test !isdefined(Erebus, :xcenter)
        @test !isdefined(Erebus, :ycenter)
        @test !isdefined(Erebus, :dx)
        @test !isdefined(Erebus, :dy)
    end

    @testset "Poroelasticity production struct defaults and validation" begin
        cfg_def = default_config()
        # Verify default compressibility fields match production calibration
        @test isapprox(cfg_def.poroelasticity.betasolid, 2.5e-11; rtol=1e-12)
        @test isapprox(cfg_def.poroelasticity.betafluid, 4.0e-10; rtol=1e-12)

        # Warning triggered when both compressibilities are zero
        cfg_zero = SimulationConfig(
            poroelasticity=PoroelasticConfig(
                betasolid=0.0,
                betafluid=0.0,
                phimin=1.0e-4,
                phimax=0.9999,
                hydrofracture=false,
                kappa_frac=1.0e3,
                gamma_frac=1.0,
                k_frac_max=1.0e-9,
                theta_frac=1.0,
                ramp_width=0.0,
                rx_floor_prefactor=1.0e-5,
            ),
        )
        @test_logs (:warn, r"Both poroelasticity.betasolid and betafluid are 0.0") validate_config(
            cfg_zero
        )
    end

    @testset "65x65 grid allocates 66x66 arrays across structures" begin
        cfg65 = SimulationConfig(
            grid=GridConfig(Nx=65, Ny=65, xsize=140_000.0, ysize=140_000.0),
            geometry=GeometryConfig(
                rplanet=50_000.0, rcrust=50_000.0, xcenter=70_000.0, ycenter=70_000.0
            ),
        )
        coords65 = GridCoordinates(cfg65.grid)
        @test coords65.Nx == 65
        @test coords65.Nx1 == 66

        props = setup_staggered_grid_properties(coords65)
        ETA = props[1]
        RHOX = props[12]
        RHO = props[30]
        @test size(ETA) == (65, 65)
        @test size(RHOX) == (66, 66)
        @test size(RHO) == (66, 66)

        helpers = setup_staggered_grid_properties_helpers(coords65)
        ETA5 = helpers[1]
        EII = helpers[8]
        @test size(ETA5) == (65, 65)
        @test size(EII) == (66, 66)

        RP, SP = setup_gravitational_lse(coords65)
        @test length(RP) == 66 * 66
        @test length(SP) == 66 * 66

        R, S = setup_hydromechanical_lse(coords65)
        @test length(R) == 66 * 66 * 6
        @test length(S) == 66 * 66 * 6

        RT, ST = setup_thermal_lse(coords65)
        @test length(RT) == 66 * 66
        @test length(ST) == 66 * 66

        interp = setup_interpolated_properties(coords65)
        @test size(interp[1]) == (65, 65)
        @test size(interp[9]) == (66, 66)

        bufs = allocate_thread_interpolation_buffers(coords65)
        @test length(bufs) == 16
        @test size(bufs[1].WTSUM) == (65, 65)
    end

    @testset "Coords-required signature contracts" begin
        # Zero-argument invocation must throw MethodError once defaults are removed
        @test_throws MethodError setup_staggered_grid_properties()
        @test_throws MethodError setup_staggered_grid_properties_helpers()
        @test_throws MethodError setup_gravitational_lse()
        @test_throws MethodError setup_hydromechanical_lse()
        @test_throws MethodError setup_thermal_lse()
        @test_throws MethodError setup_interpolated_properties()
    end

    @testset "Default grid coordinates and module constants consistency" begin
        def_coords = default_grid_coordinates()
        cfg_grid = GridConfig()
        @test def_coords.Nx == cfg_grid.Nx
        @test def_coords.Ny == cfg_grid.Ny
        @test isapprox(def_coords.xsize, cfg_grid.xsize; atol=1e-12)
        @test isapprox(def_coords.ysize, cfg_grid.ysize; atol=1e-12)
        @test isapprox(Erebus.betasolid, PoroelasticConfig().betasolid; atol=1e-18)
        @test isapprox(Erebus.betafluid, PoroelasticConfig().betafluid; atol=1e-18)
    end

    @testset "DHP grid dimension validation in update_marker_pyrolysis!" begin
        refr_cfg = RefractoryConfig(active=true, kinetics_active=true)
        coords_10 = GridCoordinates(9, 9; xsize=10000.0, ysize=10000.0)
        # Mismatched DHP matrix size: 8x8 instead of 10x10
        bad_dhp = zeros(Float64, 8, 8)
        good_dhp = zeros(Float64, 10, 10)
        tkm_test = [500.0]
        phim_test = [0.01]
        c_test = [1000.0]
        n_test = [100.0]
        h_test = [50.0]
        xm_test = [5000.0]
        ym_test = [5000.0]
        @test_throws DimensionMismatch Erebus.update_marker_pyrolysis!(
            tkm_test,
            1000.0,
            phim_test,
            c_test,
            n_test,
            h_test,
            refr_cfg;
            xm=xm_test,
            ym=ym_test,
            coords=coords_10,
            DHP=bad_dhp,
        )
        res = Erebus.update_marker_pyrolysis!(
            tkm_test,
            1000.0,
            phim_test,
            c_test,
            n_test,
            h_test,
            refr_cfg;
            xm=xm_test,
            ym=ym_test,
            coords=coords_10,
            DHP=good_dhp,
        )
        @test isapprox(res.total_dC_gas, 0.0511192; rtol=1e-4)
    end

    @testset "Volatile budget area calculation fallback" begin
        # 100 synthetic markers in a disc of radius 50 km
        nmarks = 100
        rp = 50000.0
        xc = 70000.0
        yc = 70000.0
        xm = fill(xc + 10000.0, nmarks)
        ym = fill(yc, nmarks)
        tm = fill(1, nmarks)
        tkm = fill(1500.0, nmarks)
        Xfe = fill(0.5, nmarks)
        coords_custom = GridCoordinates(33, 33; xsize=140000.0, ysize=140000.0)

        # Passing coords uses marker_area(coords)
        budgets_with_coords = Erebus.compute_core_volatile_budgets(
            xm,
            ym,
            tm,
            Xfe,
            nothing,
            nothing,
            nothing,
            nothing,
            nmarks;
            coords=coords_custom,
            xcenter=xc,
            ycenter=yc,
            rplanet=rp,
        )
        # Omitting coords uses empirical pi * rp^2 / N_planet
        budgets_no_coords = Erebus.compute_core_volatile_budgets(
            xm,
            ym,
            tm,
            Xfe,
            nothing,
            nothing,
            nothing,
            nothing,
            nmarks;
            xcenter=xc,
            ycenter=yc,
            rplanet=rp,
        )
        expected_empirical_area = (pi * rp^2) / nmarks
        expected_coords_area = coords_custom.dxm * coords_custom.dym
        ratio = budgets_with_coords.M_core_metal / budgets_no_coords.M_core_metal
        @test isapprox(ratio, expected_coords_area / expected_empirical_area; rtol=1e-10)
        expected_total_m_core = (pi * rp^2) * (2.0 * 10000.0) * 7000.0 * 0.5
        @test isapprox(budgets_no_coords.M_core_metal, expected_total_m_core; rtol=1e-10)
    end
end

using Test
using Erebus
using StaticArrays
using TOML

include("test_helpers.jl")

@testset "Config" begin
    @testset "default_config() matches baseline constants" begin
        cfg = default_config()
        @test cfg isa SimulationConfig
        @test cfg.grid.Nx == 33
        @test cfg.grid.Ny == 33
        @test cfg.grid.xsize ≈ 140000.0 rtol=1e-12
        @test cfg.grid.ysize ≈ 140000.0 rtol=1e-12
        @test cfg.geometry.rplanet ≈ 50000.0 rtol=1e-12
        @test cfg.geometry.rcrust ≈ 50000.0 rtol=1e-12
        # Poroelastic baseline default matches constants.jl test baseline (0.0)
        @test iszero(cfg.poroelasticity.betasolid)
        @test iszero(cfg.poroelasticity.betafluid)
        @test cfg.poroelasticity.phimin ≈ 1.0e-4 rtol=1e-12
        @test cfg.poroelasticity.phimax ≈ 0.9999 rtol=1e-12
        @test cfg.time.dt_longest ≈ 1.0e11 / cfg.time.yearlength
        @test cfg.time.dt_initial ≈ 1.0e11 / cfg.time.yearlength
        @test cfg.time.start_time ≈ 2.25e6 rtol=1e-12
        @test cfg.time.endtime ≈ 15.0e6 rtol=1e-12
        @test cfg.time.n_steps == 10
        @test cfg.solver.use_pardiso == false
        @test cfg.solver.dphimax ≈ 0.1

        # Material arrays must match constants.jl element-by-element
        @test cfg.materials.rhosolidm ≈ SVector{3,Float64}([3300.0, 3300.0, 1.0])
        @test cfg.materials.rhofluidm ≈ SVector{3,Float64}([1000.0, 1000.0, 1.0])
        @test cfg.materials.etasolidm ≈ SVector{3,Float64}([1.0e+19, 1.0e+19, 1.0e+16])
        @test cfg.materials.etasolidmm ≈ SVector{3,Float64}([1.0e+19, 1.0e+19, 1.0e+16])
        @test cfg.materials.etafluidm ≈ SVector{3,Float64}([1.0e+12, 1.0e+12, 1.0e-03])
        @test cfg.materials.etafluidmm ≈ SVector{3,Float64}([1.0e-03, 1.0e-03, 1.0e-03])
        @test cfg.materials.rhocpsolidm ≈ SVector{3,Float64}([3.3e+06, 3.3e+06, 3.0e+06])
        @test cfg.materials.rhocpfluidm ≈ SVector{3,Float64}([1.0e+06, 1.0e+06, 3.0e+06])
        @test cfg.materials.alphasolidm ≈ SVector{3,Float64}([3.0e-05, 3.0e-05, 0.0])
        @test cfg.materials.alphafluidm ≈ SVector{3,Float64}([5.0e-05, 5.0e-05, 0.0])
        @test cfg.materials.ksolidm ≈ SVector{3,Float64}([3.0, 3.0, 3000.0])
        @test cfg.materials.kfluidm ≈ SVector{3,Float64}([50.0, 50.0, 3000.0])
        @test cfg.materials.gggsolidm ≈ SVector{3,Float64}([1.0e+10, 1.0e+10, 1.0e+10])
        @test cfg.materials.frictsolidm ≈ SVector{3,Float64}([0.6, 0.6, 0.0])
        @test cfg.materials.cohessolidm ≈ SVector{3,Float64}([1.0e+08, 1.0e+08, 1.0e+08])
        @test cfg.materials.tenssolidm ≈ SVector{3,Float64}([6.0e+07, 6.0e+07, 6.0e+07])
        @test cfg.materials.kphim0 ≈ SVector{3,Float64}([1.0e-13, 1.0e-13, 1.0e-17])
        @test cfg.materials.tkm0 ≈ SVector{3,Float64}([170.0, 170.0, 170.0])
    end

    @testset "load_config() from file" begin
        # Default TOML file provides calibrated production configuration
        default_toml = joinpath(@__DIR__, "..", "configs", "default.toml")
        cfg_def = load_config(default_toml)
        @test cfg_def.grid.Nx == 33
        @test cfg_def.grid.xsize ≈ 140000.0 rtol=1e-12
        @test cfg_def.poroelasticity.betasolid ≈ 2.5e-11 rtol=1e-12
        @test cfg_def.poroelasticity.betafluid ≈ 4.0e-10 rtol=1e-12
        @test cfg_def.materials.ksolidm ≈ SVector{3,Float64}([3.0, 3.0, 3000.0])

        # Quick test TOML file
        quick_toml = joinpath(@__DIR__, "..", "configs", "test_quick.toml")
        cfg_q = load_config(quick_toml)
        @test cfg_q.time.n_steps == 2
        @test cfg_q.output.output_dir == "output_test"
        @test cfg_q.output.savematstep == 2
        # Verify inherited defaults for omitted sections
        @test cfg_q.grid.Nx == 33
        @test cfg_q.solver.max_plastic_iterations == 10000
        @test cfg_q.solver.max_dt_reductions == 5

        # Missing file error check
        @test_throws SystemError load_config("nonexistent_path_to_config.toml")
    end

    @testset "all shipped configs have dphimax ≈ 0.1 (N3)" begin
        configs_dir = joinpath(@__DIR__, "..", "configs")
        toml_files = filter(
            f -> endswith(f, ".toml") && f != "test_ensemble_sweep.toml",
            readdir(configs_dir),
        )
        @test length(toml_files) >= 15
        for f in toml_files
            path = joinpath(configs_dir, f)
            cfg = load_config(path)
            @test cfg.solver.dphimax ≈ 0.1
        end
    end

    @testset "load_config() from string with partial overlays" begin
        toml_str = """
        [poroelasticity]
        betasolid = 1.0e-10

        [time]
        n_steps = 5
        """
        cfg = load_config(toml_str)
        @test cfg.grid.Nx == 33
        @test cfg.grid.Ny == 33
        @test cfg.poroelasticity.betasolid ≈ 1.0e-10 rtol=1e-12
        @test iszero(cfg.poroelasticity.betafluid) # baseline default preserved
        @test cfg.time.n_steps == 5
        @test cfg.geometry.rplanet ≈ 50000.0 rtol=1e-12 # default preserved
    end

    @testset "validate_config() schema and physical bounds enforcement" begin
        # Valid baseline
        @test validate_config(default_config()) === nothing

        # Dynamic grid resolutions are valid
        @test validate_config(
            SimulationConfig(grid=GridConfig(Nx=15, Ny=15, xsize=140000.0, ysize=140000.0))
        ) === nothing
        @test validate_config(
            SimulationConfig(grid=GridConfig(Nx=65, Ny=65, xsize=140000.0, ysize=140000.0))
        ) === nothing

        # Invalid grid constraints
        @reject_config grid=GridConfig(Nx=2, Ny=33, xsize=140000.0, ysize=140000.0)
        @reject_config grid=GridConfig(Nx=33, Ny=1, xsize=140000.0, ysize=140000.0)
        @reject_config grid=GridConfig(Nx=33, Ny=33, xsize=-50000.0, ysize=140000.0)
        @reject_config grid=GridConfig(Nx=33, Ny=33, xsize=140000.0, ysize=0.0)

        # Planet exceeding domain boundary (centered)
        @reject_config(
            grid=GridConfig(Nx=33, Ny=33, xsize=80000.0, ysize=80000.0),
            geometry=GeometryConfig(
                rplanet=50000.0,
                rcrust=50000.0,
                xcenter=40000.0,
                ycenter=40000.0,
                psurface=1e3,
            ),
        )

        # Planet exceeding domain boundary (off-center placement)
        @reject_config(
            grid=GridConfig(Nx=33, Ny=33, xsize=80000.0, ysize=80000.0),
            geometry=GeometryConfig(
                rplanet=30000.0,
                rcrust=30000.0,
                xcenter=5000.0,
                ycenter=40000.0,
                psurface=1e3,
            ),
        )

        # Crust radius exceeding planet radius
        @reject_config(
            geometry=GeometryConfig(
                rplanet=40000.0,
                rcrust=50000.0,
                xcenter=70000.0,
                ycenter=70000.0,
                psurface=1e3,
            ),
        )

        # Invalid poroelastic parameters
        @reject_config poroelasticity=PoroelasticConfig(
            betasolid=-1.0e-11, betafluid=4e-10, phimin=1e-4, phimax=0.9999
        )
        @reject_config poroelasticity=PoroelasticConfig(
            betasolid=2.5e-11, betafluid=-4e-10, phimin=1e-4, phimax=0.9999
        )
        @reject_config poroelasticity=PoroelasticConfig(
            betasolid=2.5e-11, betafluid=4e-10, phimin=0.9, phimax=0.1
        )
        @reject_config poroelasticity=PoroelasticConfig(
            betasolid=2.5e-11, betafluid=4e-10, phimin=-0.1, phimax=0.9
        )
        @reject_config poroelasticity=PoroelasticConfig(kappa_frac=-1.0)
        @reject_config poroelasticity=PoroelasticConfig(gamma_frac=0.0)
        @reject_config poroelasticity=PoroelasticConfig(k_frac_max=-1.0e-9)

        # Invalid time parameters
        @reject_config(
            time=TimeConfig(
                dt_initial=-1.0,
                dt_longest=100.0,
                dtcoefdn=0.5,
                dtcoefup=1.2,
                dtstep=200,
                dxymax=0.05,
                vpratio=0.33,
                DTmax=20.0,
                yearlength=3.15e7,
                start_time=0.0,
                endtime=1000.0,
                start_step=1,
                n_steps=10,
            ),
        )
        @reject_config(
            time=TimeConfig(
                dt_initial=200.0,
                dt_longest=100.0,
                dtcoefdn=0.5,
                dtcoefup=1.2,
                dtstep=200,
                dxymax=0.05,
                vpratio=0.33,
                DTmax=20.0,
                yearlength=3.15e7,
                start_time=0.0,
                endtime=1000.0,
                start_step=1,
                n_steps=10,
            ),
        )
        @reject_config(
            time=TimeConfig(
                dt_initial=100.0,
                dt_longest=100.0,
                dtcoefdn=0.5,
                dtcoefup=1.2,
                dtstep=200,
                dxymax=0.05,
                vpratio=0.33,
                DTmax=20.0,
                yearlength=3.15e7,
                start_time=0.0,
                endtime=1000.0,
                start_step=1,
                n_steps=0,
            ),
        )
        @reject_config time=TimeConfig(start_time=-1.0)
        @reject_config time=TimeConfig(start_time=10.0, endtime=5.0)

        # Invalid solver parameters: max_dt_reductions must be >= 1
        @reject_config solver=SolverConfig(max_dt_reductions=0)

        # Invalid output parameters: savematstep and visstep must be >= 1
        @reject_config output=OutputConfig(savematstep=0, visstep=1)
        @reject_config output=OutputConfig(savematstep=10, visstep=0)
        @reject_config output=OutputConfig(telemetrystep=0)
        @reject_config output=OutputConfig(mode=:invalid_mode)
        @reject_config output=OutputConfig(telemetry_file="")

        # Invalid solver parameters
        @reject_config solver=SolverConfig(p2m_mode=:unsupported)
        @reject_config solver=SolverConfig(tile_size=1)
        @reject_config solver=SolverConfig(hydromech_solver=:unknown_solver)
        @reject_config solver=SolverConfig(krylov_method=:unknown_method)
        @reject_config solver=SolverConfig(krylov_rtol=-1.0e-5)
        @reject_config solver=SolverConfig(krylov_atol=-1.0e-5)
        @reject_config solver=SolverConfig(krylov_maxiter=0)
        @reject_config solver=SolverConfig(krylov_restart=0)
        @reject_config solver=SolverConfig(preconditioner=:unknown_prec)

        # Invalid thermodynamics
        @reject_config thermodynamics=ThermalConfig(ratio_al=-0.1)
        @reject_config thermodynamics=ThermalConfig(
            tmfluidphase=1500.0, tmsolidphase=1400.0
        )
    end

    @testset "typo protection and unknown key rejection" begin
        # Unknown section
        bad_section_toml = """
        [unknown_section]
        foo = 123
        """
        @test_throws ArgumentError load_config(bad_section_toml)

        # Unknown key inside valid section
        bad_key_toml = """
        [poroelasticity]
        betasolidd = 1.0e-11
        """
        @test_throws ArgumentError load_config(bad_key_toml)

        # SVector length mismatch
        bad_svector_toml = """
        [materials]
        rhosolidm = [3300.0, 3300.0]
        """
        @test_throws ArgumentError load_config(bad_svector_toml)
    end

    @testset "save_config() roundtrip serialization" begin
        cfg_orig = SimulationConfig(
            poroelasticity=PoroelasticConfig(
                betasolid=3.0e-11, betafluid=5.0e-10, phimin=2e-4, phimax=0.99
            ),
            time=TimeConfig(
                dt_initial=1000.0,
                dt_longest=1000.0,
                dtcoefdn=0.5,
                dtcoefup=1.2,
                dtstep=200,
                dxymax=0.05,
                vpratio=0.33,
                DTmax=20.0,
                yearlength=3.15e7,
                start_time=0.0,
                endtime=1.0e6,
                start_step=1,
                n_steps=4,
            ),
            output=OutputConfig(output_dir="test_roundtrip", savematstep=2, visstep=1),
        )

        tmp_path = tempname() * ".toml"
        try
            save_config(tmp_path, cfg_orig)
            @test isfile(tmp_path)
            cfg_loaded = load_config(tmp_path)

            @test cfg_loaded.grid.Nx == cfg_orig.grid.Nx
            @test cfg_loaded.grid.Ny == cfg_orig.grid.Ny
            @test cfg_loaded.grid.xsize ≈ cfg_orig.grid.xsize
            @test cfg_loaded.grid.ysize ≈ cfg_orig.grid.ysize
            @test cfg_loaded.poroelasticity.betasolid ≈ cfg_orig.poroelasticity.betasolid
            @test cfg_loaded.poroelasticity.betafluid ≈ cfg_orig.poroelasticity.betafluid
            @test cfg_loaded.poroelasticity.phimin ≈ cfg_orig.poroelasticity.phimin
            @test cfg_loaded.poroelasticity.phimax ≈ cfg_orig.poroelasticity.phimax
            @test cfg_loaded.time.n_steps == cfg_orig.time.n_steps
            @test cfg_loaded.output.output_dir == cfg_orig.output.output_dir
            @test cfg_loaded.materials.rhosolidm ≈ cfg_orig.materials.rhosolidm
        finally
            rm(tmp_path, force=true)
        end
    end

    @testset "materials and thermodynamics validation" begin
        # Zero or negative shear modulus
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                materials=MaterialConfig(gggsolidm=SVector{3,Float64}([0.0, 1e10, 1e10]))
            ),
        )
        # Negative conductivity
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                materials=MaterialConfig(ksolidm=SVector{3,Float64}([-3.0, 3.0, 3000.0]))
            ),
        )
        # Invalid radiogenic half-life
        @reject_config thermodynamics=ThermalConfig(t_half_al=-1.0)
        # Material modification away from compiled constants
        @reject_config materials=MaterialConfig(
            rhosolidm=SVector{3,Float64}([4000.0, 3300.0, 1.0])
        )
        # Radiogenic heating calculation keyword arguments and toggling
        hr_sol_on, _, _ = Erebus.calculate_radioactive_heating(true, false, 0.0)
        @test hr_sol_on[1] > 0.0
        hr_sol_off, _, _ = Erebus.calculate_radioactive_heating(false, false, 0.0)
        @test all(iszero, hr_sol_off)
    end

    @testset "output restart_from validation" begin
        # Non-existent checkpoint file
        @reject_config output=OutputConfig(restart_from="nonexistent_checkpoint.jld2")
        # Non-.jld2 extension
        tmp_txt = tempname() * ".txt"
        touch(tmp_txt)
        try
            @reject_config output=OutputConfig(restart_from=tmp_txt)
        finally
            rm(tmp_txt, force=true)
        end
    end

    @testset "Unphysical Disk Configurations" begin
        @reject_config disk=DiskConfig(orbital_distance_au=-1.0)
        @reject_config disk=DiskConfig(stellar_mass_msun=-0.5)
        @reject_config disk=DiskConfig(model=:invalid_disk_model)
        @reject_config disk=DiskConfig(t_visc_0_myr=-0.1)
        @reject_config disk=DiskConfig(gamma=-1.0)
        @reject_config disk=DiskConfig(orbital_distance_au=NaN)
        @reject_config disk=DiskConfig(orbital_distance_au=Inf)
        @reject_config disk=DiskConfig(stellar_mass_msun=NaN)
        @reject_config disk=DiskConfig(t_visc_0_myr=Inf)
    end

    @testset "Unphysical Thermodynamics Configurations" begin
        @reject_config thermodynamics=ThermalConfig(Lᶠ=-100.0)
        @reject_config thermodynamics=ThermalConfig(ratio_al=1.5)
        @reject_config thermodynamics=ThermalConfig(ratio_fe=-0.05)
        @reject_config thermodynamics=ThermalConfig(emissivity=-0.1)
        @reject_config thermodynamics=ThermalConfig(emissivity=1.2)
        @reject_config thermodynamics=ThermalConfig(sigma_sb=-1.0)
        @reject_config thermodynamics=ThermalConfig(sigma_sb=NaN)
        @reject_config thermodynamics=ThermalConfig(sigma_sb=Inf)
        @reject_config thermodynamics=ThermalConfig(fluid_viscosity_mode=:invalid_mode)
        @reject_config thermodynamics=ThermalConfig(fluid_viscosity_T0=-10.0)
        @reject_config thermodynamics=ThermalConfig(fluid_viscosity_T0=NaN)
        @reject_config thermodynamics=ThermalConfig(fluid_viscosity_eta0=-1.0e-3)
        @reject_config thermodynamics=ThermalConfig(fluid_viscosity_eta0=Inf)
    end

    @testset "Unphysical Solver Configurations" begin
        @reject_config solver=SolverConfig(etamin=-1.0)
        @reject_config solver=SolverConfig(etamin=10.0, etamax=1.0)
        @reject_config solver=SolverConfig(etaphikoef=-0.1)
        @reject_config solver=SolverConfig(max_plastic_iterations=0)
    end

    @testset "Hydrothermal Configurations" begin
        cfg_default = HydrothermalConfig()
        @test cfg_default.active == false
        @test cfg_default.phi_start ≈ 0.30

        @reject_config hydrothermal=HydrothermalConfig(
            active=true, phi_start=0.8, phi_end=0.2
        )
        @reject_config hydrothermal=HydrothermalConfig(active=true, Ra_m_crit=-1.0)
        @reject_config hydrothermal=HydrothermalConfig(
            active=true, k_floor=10.0, k_cutoff=1.0
        )
        @reject_config hydrothermal=HydrothermalConfig(active=true, H_layer_min=-10.0)
        @reject_config hydrothermal=HydrothermalConfig(
            active=true, H_layer_min=20000.0, H_layer=10000.0
        )
    end

    @testset "Accretion Configurations" begin
        cfg_default = AccretionConfig()
        @test cfg_default.active == false
        @test cfg_default.mode === :pebble_hill
        @test cfg_default.M_initial ≈ 1.0e17
        @test cfg_default.R_initial ≈ 20000.0
        @test cfg_default.rho_bulk ≈ 3000.0
        @test cfg_default.M_target ≈ 1.0e20
        @test cfg_default.R_target ≈ 50000.0
        @test cfg_default.h_impact ≈ 0.5
        @test cfg_default.phi_accreted ≈ 0.35
        @test cfg_default.Xfe_bulk_accreted ≈ 0.10
        @test cfg_default.dM_dt_constant ≈ 1.5e6
        @test cfg_default.dR_dt_constant ≈ 5.0e-10

        # Valid config passes validation
        @test validate_config(SimulationConfig(accretion=AccretionConfig(active=true))) ===
            nothing

        # Invalid modes
        @reject_config accretion=AccretionConfig(active=true, mode=:invalid_mode)

        # Invalid mass and radius bounds
        @reject_config accretion=AccretionConfig(active=true, M_initial=-1.0)
        @reject_config accretion=AccretionConfig(
            active=true, M_target=1.0e16, M_initial=1.0e17
        )
        @reject_config accretion=AccretionConfig(active=true, R_initial=-1000.0)
        @reject_config accretion=AccretionConfig(
            active=true, R_target=10000.0, R_initial=20000.0
        )
        @reject_config accretion=AccretionConfig(active=true, rho_bulk=-1000.0)

        # Target radius exceeding domain boundary
        @reject_config accretion=AccretionConfig(active=true, R_target=80000.0)

        # Invalid timing and growth parameters
        @reject_config accretion=AccretionConfig(active=true, t_start_myr=-0.5)
        @reject_config accretion=AccretionConfig(active=true, t_duration_myr=0.0)
        @reject_config accretion=AccretionConfig(active=true, dM_dt_constant=0.0)
        @reject_config accretion=AccretionConfig(active=true, dR_dt_constant=-1.0e-5)
        @reject_config accretion=AccretionConfig(active=true, tau_growth_myr=0.0)

        # Invalid physical fractions
        @reject_config accretion=AccretionConfig(active=true, h_impact=-0.1)
        @reject_config accretion=AccretionConfig(active=true, h_impact=1.5)
        @reject_config accretion=AccretionConfig(active=true, phi_accreted=-0.1)
        @reject_config accretion=AccretionConfig(active=true, phi_accreted=1.2)
        @reject_config accretion=AccretionConfig(active=true, Xfe_bulk_accreted=-0.05)
        @reject_config accretion=AccretionConfig(active=true, Xfe_bulk_accreted=1.05)

        # Invalid volatile abundances
        @reject_config accretion=AccretionConfig(active=true, XC_accreted_ppm=-10.0)
        @reject_config accretion=AccretionConfig(active=true, XN_accreted_ppm=-1.0)
        @reject_config accretion=AccretionConfig(active=true, XS_accreted_ppm=-50.0)

        # Invalid pebble/turbulence parameters
        @reject_config accretion=AccretionConfig(active=true, stokes_number=0.0)
        @reject_config accretion=AccretionConfig(active=true, alpha_turbulence=-1e-3)

        # TOML deserialization roundtrip
        toml_overlay = """
        [accretion]
        active = true
        mode = "safronov"
        M_initial = 2.0e17
        R_initial = 25000.0
        h_impact = 0.8
        Sigma_pl_0 = 150.0
        """
        cfg_parsed = load_config(toml_overlay)
        @test cfg_parsed.accretion.active == true
        @test cfg_parsed.accretion.mode === :safronov
        @test cfg_parsed.accretion.M_initial ≈ 2.0e17
        @test cfg_parsed.accretion.R_initial ≈ 25000.0
        @test cfg_parsed.accretion.h_impact ≈ 0.8
        @test cfg_parsed.accretion.Sigma_pl_0 ≈ 150.0
    end

    @testset "@reject_config and @unpack_coords helper validation" begin
        # 1. Test @reject_config with both kwarg and raw expression forms
        @reject_config grid=GridConfig(Nx=1)
        @reject_config SimulationConfig(grid=GridConfig(Ny=1))

        # 2. Test @unpack_coords with concrete GridCoordinates
        test_coords = GridCoordinates(33, 33; xsize=100000.0, ysize=100000.0)
        @unpack_coords test_coords dx dy Nx Ny
        @test isapprox(dx_val, test_coords.dx; atol=1e-12)
        @test isapprox(dy_val, test_coords.dy; atol=1e-12)
        @test Nx_val == 33
        @test Ny_val == 33

        # 3. Test @unpack_coords with nothing (falls back to caller variables)
        dx, dy = 1000.0, 2000.0
        no_coords = nothing
        @unpack_coords no_coords dx dy
        @test isapprox(dx_val, 1000.0; atol=1e-12)
        @test isapprox(dy_val, 2000.0; atol=1e-12)
    end

    @testset "RedoxConfig & RefractoryConfig Validation and TOML Roundtrip" begin
        # 1. RedoxConfig constructor validation
        @test_throws ArgumentError RedoxConfig(; reference=:INVALID)
        @test_throws ArgumentError RedoxConfig(; deltaIW_min=3.0, deltaIW_max=1.0)
        @test_throws DomainError RedoxConfig(; initial_x_ferric=-0.1)
        @test_throws DomainError RedoxConfig(; initial_x_ferric=1.1)
        @test_throws DomainError RedoxConfig(; w_graphite_threshold=-1.0)
        @test_throws DomainError RedoxConfig(; w_graphite_threshold=1.5)

        # 2. RedoxConfig & RefractoryConfig constructor domain bounds
        @test_throws DomainError RedoxConfig(; active=true, initial_x_ferric=-0.05)
        @test_throws DomainError RefractoryConfig(;
            active=true, kinetics_active=true, A_C=0.0
        )
        @test_throws DomainError RefractoryConfig(;
            active=true, kinetics_active=true, Ea_C=-1.0
        )
        @test_throws DomainError RefractoryConfig(;
            active=true, kinetics_active=true, dh_pyro_C=-100.0
        )
        @test_throws DomainError RefractoryConfig(; active=true, T_pyro_min=-50.0)

        # 4. TOML load/save roundtrip for redox and refractory
        toml_overlay = """
        [redox]
        active = true
        reference = "crust"
        serpentinization_redox = true
        segregation_redox = true
        venting_redox = true
        deltaIW_min = -5.0
        deltaIW_max = 5.0
        initial_x_ferric = 0.08
        pyrolysis_redox = false
        graphite_buffer_active = false
        w_graphite_threshold = 2.0e-5

        [refractory]
        active = true
        kinetics_active = true
        A_C = 2.0e14
        Ea_C = 2.1e5
        dh_pyro_C = 4.5e5
        T_pyro_min = 320.0
        """
        cfg_parsed = load_config(toml_overlay)
        validate_config(cfg_parsed)
        @test cfg_parsed.redox.active == true
        @test cfg_parsed.redox.reference === :crust
        @test isapprox(cfg_parsed.redox.deltaIW_min, -5.0)
        @test isapprox(cfg_parsed.redox.deltaIW_max, 5.0)
        @test isapprox(cfg_parsed.redox.initial_x_ferric, 0.08)
        @test cfg_parsed.redox.pyrolysis_redox == false
        @test cfg_parsed.redox.graphite_buffer_active == false
        @test isapprox(cfg_parsed.redox.w_graphite_threshold, 2.0e-5)
        @test cfg_parsed.refractory.active == true
        @test cfg_parsed.refractory.kinetics_active == true
        @test isapprox(cfg_parsed.refractory.A_C, 2.0e14)
        @test isapprox(cfg_parsed.refractory.Ea_C, 2.1e5)
        @test isapprox(cfg_parsed.refractory.dh_pyro_C, 4.5e5)
        @test isapprox(cfg_parsed.refractory.T_pyro_min, 320.0)

        # Roundtrip via save_config
        dict_cfg = config_to_dict(cfg_parsed)
        @test haskey(dict_cfg, "redox")
        @test haskey(dict_cfg, "refractory")
        @test dict_cfg["redox"]["reference"] == "crust"
        @test dict_cfg["redox"]["pyrolysis_redox"] == false
        @test dict_cfg["redox"]["graphite_buffer_active"] == false
        @test isapprox(dict_cfg["redox"]["w_graphite_threshold"], 2.0e-5; atol=1e-12)
        @test isapprox(dict_cfg["refractory"]["T_pyro_min"], 320.0; atol=1e-12)
    end

    @testset "MagmaOceanDegassingConfig & EscapeConfig XUV validation and TOML roundtrip" begin
        # 1. MagmaOceanDegassingConfig constructor validation
        cfg_mo = MagmaOceanDegassingConfig(;
            active=true,
            mode=:equilibrium,
            F_melt_threshold=0.5,
            degas_depth_fraction=0.85,
            efficiency=0.9,
        )
        @test cfg_mo.active == true
        @test cfg_mo.mode === :equilibrium
        @test isapprox(cfg_mo.F_melt_threshold, 0.5)
        @test isapprox(cfg_mo.degas_depth_fraction, 0.85)
        @test isapprox(cfg_mo.efficiency, 0.9)

        # Domain errors for MagmaOceanDegassingConfig
        @test_throws ArgumentError MagmaOceanDegassingConfig(; mode=:invalid_mode)
        @test_throws DomainError MagmaOceanDegassingConfig(; F_melt_threshold=-0.1)
        @test_throws DomainError MagmaOceanDegassingConfig(; F_melt_threshold=1.1)
        @test_throws DomainError MagmaOceanDegassingConfig(; degas_depth_fraction=-0.05)
        @test_throws DomainError MagmaOceanDegassingConfig(; degas_depth_fraction=1.05)
        @test_throws DomainError MagmaOceanDegassingConfig(; efficiency=0.0)
        @test_throws DomainError MagmaOceanDegassingConfig(; efficiency=1.5)

        # 2. EscapeConfig XUV domain errors
        @test_throws DomainError EscapeConfig(; epsilon_xuv=0.0)
        @test_throws DomainError EscapeConfig(; epsilon_xuv=-0.1)
        @test_throws DomainError EscapeConfig(; F_xuv_1au_sat=-1.0)
        @test_throws DomainError EscapeConfig(; t_sat_yr=0.0)
        @test_throws DomainError EscapeConfig(; beta_xuv=-0.5)
        @test_throws DomainError EscapeConfig(; r_xuv_ratio=0.8)

        # 3. TOML loading and roundtrip
        toml_overlay = """
        [magma_degassing]
        active = true
        mode = "equilibrium"
        F_melt_threshold = 0.35
        degas_depth_fraction = 0.92
        water_As = 0.45
        redox_coupled = true
        efficiency = 0.85

        [escape]
        active = true
        xuv_driven = true
        epsilon_xuv = 0.20
        F_xuv_1au_sat = 2.5
        t_sat_yr = 5.0e7
        beta_xuv = 1.15
        r_xuv_ratio = 1.05
        tidal_correction = true
        """
        cfg_loaded = load_config(toml_overlay)
        validate_config(cfg_loaded)
        @test cfg_loaded.magma_degassing.active == true
        @test cfg_loaded.magma_degassing.mode === :equilibrium
        @test isapprox(cfg_loaded.magma_degassing.F_melt_threshold, 0.35)
        @test isapprox(cfg_loaded.magma_degassing.degas_depth_fraction, 0.92)
        @test isapprox(cfg_loaded.magma_degassing.water_As, 0.45)
        @test isapprox(cfg_loaded.magma_degassing.efficiency, 0.85)

        @test cfg_loaded.escape.active == true
        @test cfg_loaded.escape.xuv_driven == true
        @test isapprox(cfg_loaded.escape.epsilon_xuv, 0.20)
        @test isapprox(cfg_loaded.escape.F_xuv_1au_sat, 2.5)
        @test isapprox(cfg_loaded.escape.t_sat_yr, 5.0e7)
        @test isapprox(cfg_loaded.escape.beta_xuv, 1.15)
        @test isapprox(cfg_loaded.escape.r_xuv_ratio, 1.05)
        @test cfg_loaded.escape.tidal_correction == true

        # Roundtrip via dict
        dict_repr = config_to_dict(cfg_loaded)
        @test haskey(dict_repr, "magma_degassing")
        @test haskey(dict_repr, "escape")
        @test dict_repr["magma_degassing"]["mode"] == "equilibrium"
        @test isapprox(dict_repr["escape"]["epsilon_xuv"], 0.20)
    end

    @testset "PR 1b: MagmaOceanDegassingConfig deprecation and rejection" begin
        # 1. Direct constructor with crystallization_degassing throws ArgumentError
        @test_throws ArgumentError MagmaOceanDegassingConfig(;
            crystallization_degassing=true
        )

        # 2. TOML parser rejects crystallization_degassing key with descriptive error naming it
        bad_degas_toml = """
        [magma_degassing]
        crystallization_degassing = true
        """
        err = try
            load_config(bad_degas_toml)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin(
            "crystallization_degassing has been removed; saturation is evaluated in the melt frame",
            sprint(showerror, err),
        )
    end

    @testset "SolverConfig Validation, Cross-Validation, and TOML Roundtrip" begin
        # 1. Invalid preconditioner rejection
        @reject_config solver=SolverConfig(preconditioner=:block_jacobi)
        @reject_config solver=SolverConfig(preconditioner=:unknown_prec)

        # 2. Matrix-free cross-validation constraints
        @reject_config solver=SolverConfig(
            hydromech_solver=:matrix_free, darcy_elimination=false
        )
        @reject_config solver=SolverConfig(
            hydromech_solver=:matrix_free, darcy_elimination=true
        ) poroelasticity=PoroelasticConfig(hydrofracture=true)
        @reject_config solver=SolverConfig(
            hydromech_solver=:matrix_free, darcy_elimination=true
        ) venting=VentingConfig(active=true)

        # 3. Valid configurations pass validation
        cfg_iter = load_config("""
        [solver]
        hydromech_solver = "iterative"
        krylov_method = "fgmres"
        krylov_rtol = 1.0e-7
        krylov_atol = 1.0e-11
        krylov_maxiter = 150
        krylov_restart = 30
        darcy_elimination = true
        preconditioner = "diagonal"
        """)
        validate_config(cfg_iter)
        @test cfg_iter.solver.hydromech_solver === :iterative
        @test cfg_iter.solver.krylov_method === :fgmres
        @test isapprox(cfg_iter.solver.krylov_rtol, 1.0e-7)
        @test isapprox(cfg_iter.solver.krylov_atol, 1.0e-11)
        @test cfg_iter.solver.krylov_maxiter == 150
        @test cfg_iter.solver.krylov_restart == 30
        @test cfg_iter.solver.darcy_elimination == true
        @test cfg_iter.solver.preconditioner === :diagonal

        cfg_mf = load_config("""
        [solver]
        hydromech_solver = "matrix_free"
        krylov_method = "gmres"
        darcy_elimination = true
        preconditioner = "block_schur"
        """)
        validate_config(cfg_mf)
        @test cfg_mf.solver.hydromech_solver === :matrix_free
        @test cfg_mf.solver.darcy_elimination == true
        @test cfg_mf.solver.preconditioner === :block_schur

        # 4. TOML load and roundtrip
        toml_solver = """
        [solver]
        hydromech_solver = "iterative"
        krylov_method = "bicgstab"
        krylov_rtol = 1.0e-8
        krylov_atol = 1.0e-12
        krylov_maxiter = 300
        krylov_restart = 40
        darcy_elimination = true
        preconditioner = "diagonal"
        """
        cfg_loaded = load_config(toml_solver)
        validate_config(cfg_loaded)
        @test cfg_loaded.solver.hydromech_solver === :iterative
        @test cfg_loaded.solver.krylov_method === :bicgstab
        @test isapprox(cfg_loaded.solver.krylov_rtol, 1.0e-8)
        @test isapprox(cfg_loaded.solver.krylov_atol, 1.0e-12)
        @test cfg_loaded.solver.krylov_maxiter == 300
        @test cfg_loaded.solver.krylov_restart == 40
        @test cfg_loaded.solver.darcy_elimination == true
        @test cfg_loaded.solver.preconditioner === :diagonal

        dict_repr = config_to_dict(cfg_loaded)
        @test haskey(dict_repr, "solver")
        @test dict_repr["solver"]["hydromech_solver"] == "iterative"
        @test dict_repr["solver"]["preconditioner"] == "diagonal"

        # 5. Multigrid configuration and validation
        cfg_mg = load_config("""
        [solver]
        hydromech_solver = "iterative"
        krylov_method = "fgmres"
        preconditioner = "multigrid"
        mg_levels = 5
        mg_pre_smooth = 3
        mg_post_smooth = 3
        mg_smoother = "redblack_gauss_seidel"
        mg_omega = 0.75
        """)
        validate_config(cfg_mg)
        @test cfg_mg.solver.preconditioner === :multigrid
        @test cfg_mg.solver.mg_levels == 5
        @test cfg_mg.solver.mg_pre_smooth == 3
        @test cfg_mg.solver.mg_post_smooth == 3
        @test cfg_mg.solver.mg_smoother === :redblack_gauss_seidel
        @test isapprox(cfg_mg.solver.mg_omega, 0.75)

        @reject_config solver=SolverConfig(mg_levels=0)
        @reject_config solver=SolverConfig(mg_pre_smooth=0)
        @reject_config solver=SolverConfig(mg_post_smooth=0)
        @reject_config solver=SolverConfig(mg_omega=0.0)
        @reject_config solver=SolverConfig(mg_omega=1.5)
        @reject_config solver=SolverConfig(mg_smoother=:invalid_smoother)
    end

    @testset "tools/check_config_schema.jl" begin
        cmd = `$(Base.julia_cmd()) --project=$(normpath(joinpath(@__DIR__, ".."))) $(joinpath(@__DIR__, "..", "tools", "check_config_schema.jl")) --check`
        out = read(cmd, String)
        @test occursin("Schema verification passed", out)
        @test occursin("485 configuration fields", out)
    end
end

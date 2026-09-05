using Test
using Erebus
using StaticArrays
using TOML

@testset "Config" begin
    @testset "default_config() matches baseline constants" begin
        cfg = default_config()
        @test cfg isa SimulationConfig
        @test cfg.grid.Nx == 33
        @test cfg.grid.Ny == 33
        @test cfg.grid.xsize == 140000.0
        @test cfg.grid.ysize == 140000.0
        @test cfg.geometry.rplanet == 50000.0
        @test cfg.geometry.rcrust == 50000.0
        # Poroelastic baseline default matches constants.jl test baseline (0.0)
        @test cfg.poroelasticity.betasolid == 0.0
        @test cfg.poroelasticity.betafluid == 0.0
        @test cfg.poroelasticity.phimin == 1.0e-4
        @test cfg.poroelasticity.phimax == 0.9999
        @test cfg.time.dt_longest ≈ 1.0e11 / cfg.time.yearlength
        @test cfg.time.dt_initial ≈ 1.0e11 / cfg.time.yearlength
        @test cfg.time.start_time == 2.25e6
        @test cfg.time.endtime == 15.0e6
        @test cfg.time.n_steps == 10
        @test cfg.solver.use_pardiso == false

        # Material arrays must match constants.jl element-by-element
        @test cfg.materials.rhosolidm == SVector{3,Float64}([3300.0, 3300.0, 1.0])
        @test cfg.materials.rhofluidm == SVector{3,Float64}([1000.0, 1000.0, 1.0])
        @test cfg.materials.etasolidm == SVector{3,Float64}([1.0e+19, 1.0e+19, 1.0e+16])
        @test cfg.materials.etasolidmm == SVector{3,Float64}([1.0e+19, 1.0e+19, 1.0e+16])
        @test cfg.materials.etafluidm == SVector{3,Float64}([1.0e+12, 1.0e+12, 1.0e-03])
        @test cfg.materials.etafluidmm == SVector{3,Float64}([1.0e-03, 1.0e-03, 1.0e-03])
        @test cfg.materials.rhocpsolidm == SVector{3,Float64}([3.3e+06, 3.3e+06, 3.0e+06])
        @test cfg.materials.rhocpfluidm == SVector{3,Float64}([1.0e+06, 1.0e+06, 3.0e+06])
        @test cfg.materials.alphasolidm == SVector{3,Float64}([3.0e-05, 3.0e-05, 0.0])
        @test cfg.materials.alphafluidm == SVector{3,Float64}([5.0e-05, 5.0e-05, 0.0])
        @test cfg.materials.ksolidm == SVector{3,Float64}([3.0, 3.0, 3000.0])
        @test cfg.materials.kfluidm == SVector{3,Float64}([50.0, 50.0, 3000.0])
        @test cfg.materials.gggsolidm == SVector{3,Float64}([1.0e+10, 1.0e+10, 1.0e+10])
        @test cfg.materials.frictsolidm == SVector{3,Float64}([0.6, 0.6, 0.0])
        @test cfg.materials.cohessolidm == SVector{3,Float64}([1.0e+08, 1.0e+08, 1.0e+08])
        @test cfg.materials.tenssolidm == SVector{3,Float64}([6.0e+07, 6.0e+07, 6.0e+07])
        @test cfg.materials.kphim0 == SVector{3,Float64}([1.0e-13, 1.0e-13, 1.0e-17])
        @test cfg.materials.tkm0 == SVector{3,Float64}([170.0, 170.0, 170.0])
    end

    @testset "load_config() from file" begin
        # Default TOML file provides calibrated production configuration
        default_toml = joinpath(@__DIR__, "..", "configs", "default.toml")
        cfg_def = load_config(default_toml)
        @test cfg_def.grid.Nx == 33
        @test cfg_def.grid.xsize == 140000.0
        @test cfg_def.poroelasticity.betasolid == 2.5e-11
        @test cfg_def.poroelasticity.betafluid == 4.0e-10
        @test cfg_def.materials.ksolidm == SVector{3,Float64}([3.0, 3.0, 3000.0])

        # Quick test TOML file
        quick_toml = joinpath(@__DIR__, "..", "configs", "test_quick.toml")
        cfg_q = load_config(quick_toml)
        @test cfg_q.time.n_steps == 2
        @test cfg_q.output.output_dir == "output_test"
        @test cfg_q.output.savematstep == 2
        # Verify inherited defaults for omitted sections
        @test cfg_q.grid.Nx == 33
        @test cfg_q.geometry.rplanet == 50000.0
        @test cfg_q.solver.titermax == 10000

        # Missing file error check
        @test_throws SystemError load_config("nonexistent_path_to_config.toml")
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
        @test cfg.poroelasticity.betasolid == 1.0e-10
        @test cfg.poroelasticity.betafluid == 0.0 # baseline default preserved
        @test cfg.time.n_steps == 5
        @test cfg.geometry.rplanet == 50000.0 # default preserved
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
        @test_throws ArgumentError validate_config(
            SimulationConfig(grid=GridConfig(Nx=2, Ny=33, xsize=140000.0, ysize=140000.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(grid=GridConfig(Nx=33, Ny=1, xsize=140000.0, ysize=140000.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(grid=GridConfig(Nx=33, Ny=33, xsize=-50000.0, ysize=140000.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(grid=GridConfig(Nx=33, Ny=33, xsize=140000.0, ysize=0.0))
        )

        # Planet exceeding domain boundary (centered)
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                grid=GridConfig(Nx=33, Ny=33, xsize=80000.0, ysize=80000.0),
                geometry=GeometryConfig(
                    rplanet=50000.0,
                    rcrust=50000.0,
                    xcenter=40000.0,
                    ycenter=40000.0,
                    psurface=1e3,
                ),
            ),
        )

        # Planet exceeding domain boundary (off-center placement)
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                grid=GridConfig(Nx=33, Ny=33, xsize=80000.0, ysize=80000.0),
                geometry=GeometryConfig(
                    rplanet=30000.0,
                    rcrust=30000.0,
                    xcenter=5000.0,
                    ycenter=40000.0,
                    psurface=1e3,
                ),
            ),
        )

        # Crust radius exceeding planet radius
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                geometry=GeometryConfig(
                    rplanet=40000.0,
                    rcrust=50000.0,
                    xcenter=70000.0,
                    ycenter=70000.0,
                    psurface=1e3,
                ),
            ),
        )

        # Invalid poroelastic parameters
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                poroelasticity=PoroelasticConfig(
                    betasolid=-1.0e-11, betafluid=4e-10, phimin=1e-4, phimax=0.9999
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                poroelasticity=PoroelasticConfig(
                    betasolid=2.5e-11, betafluid=-4e-10, phimin=1e-4, phimax=0.9999
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                poroelasticity=PoroelasticConfig(
                    betasolid=2.5e-11, betafluid=4e-10, phimin=0.9, phimax=0.1
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                poroelasticity=PoroelasticConfig(
                    betasolid=2.5e-11, betafluid=4e-10, phimin=-0.1, phimax=0.9
                ),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(poroelasticity=PoroelasticConfig(kappa_frac=-1.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(poroelasticity=PoroelasticConfig(gamma_frac=0.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(poroelasticity=PoroelasticConfig(k_frac_max=-1.0e-9))
        )

        # Invalid time parameters
        @test_throws ArgumentError validate_config(
            SimulationConfig(
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
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
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
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
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
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(time=TimeConfig(start_time=-1.0))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(time=TimeConfig(start_time=10.0, endtime=5.0))
        )

        # Invalid solver parameters: titermax must be <= nplast to prevent plastic array overflow
        @test_throws ArgumentError validate_config(
            SimulationConfig(solver=SolverConfig(titermax=200_000, nplast=100_000))
        )

        # Invalid output parameters: savematstep and visstep must be >= 1
        @test_throws ArgumentError validate_config(
            SimulationConfig(output=OutputConfig(savematstep=0, visstep=1))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(output=OutputConfig(savematstep=10, visstep=0))
        )

        # Invalid thermodynamics
        @test_throws ArgumentError validate_config(
            SimulationConfig(thermodynamics=ThermalConfig(ratio_al=-0.1))
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                thermodynamics=ThermalConfig(tmfluidphase=1500.0, tmsolidphase=1400.0)
            ),
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
            @test cfg_loaded.grid.xsize == cfg_orig.grid.xsize
            @test cfg_loaded.grid.ysize == cfg_orig.grid.ysize
            @test cfg_loaded.poroelasticity.betasolid == cfg_orig.poroelasticity.betasolid
            @test cfg_loaded.poroelasticity.betafluid == cfg_orig.poroelasticity.betafluid
            @test cfg_loaded.poroelasticity.phimin == cfg_orig.poroelasticity.phimin
            @test cfg_loaded.poroelasticity.phimax == cfg_orig.poroelasticity.phimax
            @test cfg_loaded.time.n_steps == cfg_orig.time.n_steps
            @test cfg_loaded.output.output_dir == cfg_orig.output.output_dir
            @test cfg_loaded.materials.rhosolidm == cfg_orig.materials.rhosolidm
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
        @test_throws ArgumentError validate_config(
            SimulationConfig(thermodynamics=ThermalConfig(t_half_al=-1.0))
        )
        # Material modification away from compiled constants
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                materials=MaterialConfig(
                    rhosolidm=SVector{3,Float64}([4000.0, 3300.0, 1.0])
                ),
            ),
        )
        # Radiogenic heating calculation keyword arguments and toggling
        hr_sol_on, _ = Erebus.calculate_radioactive_heating(true, false, 0.0)
        @test hr_sol_on[1] > 0.0
        hr_sol_off, _ = Erebus.calculate_radioactive_heating(false, false, 0.0)
        @test all(hr_sol_off .== 0.0)
    end

    @testset "output restart_from validation" begin
        # Non-existent checkpoint file
        @test_throws ArgumentError validate_config(
            SimulationConfig(
                output=OutputConfig(restart_from="nonexistent_checkpoint.jld2")
            ),
        )
        # Non-.jld2 extension
        tmp_txt = tempname() * ".txt"
        touch(tmp_txt)
        try
            @test_throws ArgumentError validate_config(
                SimulationConfig(output=OutputConfig(restart_from=tmp_txt))
            )
        finally
            rm(tmp_txt, force=true)
        end
    end
end

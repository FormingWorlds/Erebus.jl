using Erebus
using Erebus.Numerics: assemble_gravitational_lse!
using Test

@testset "Dynamic Grid Resolution" begin
    @testset "GridCoordinates construction across resolutions" begin
        # 17x17 grid
        c17 = GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        @test c17.Nx == 17
        @test c17.Ny == 17
        @test c17.Nx1 == 18
        @test c17.Ny1 == 18
        @test length(c17.x) == 17
        @test length(c17.y) == 17
        @test length(c17.xp) == 18
        @test length(c17.yp) == 18
        @test length(c17.xvx) == 18
        @test length(c17.yvx) == 18
        @test length(c17.xvy) == 18
        @test length(c17.yvy) == 18
        @test c17.dx ≈ 140_000.0 / 16
        @test c17.dy ≈ 140_000.0 / 16
        @test c17.Nxm == 16 * 4
        @test c17.Nym == 16 * 4
        @test c17.start_marknum == (16 * 4)^2
        @test length(c17.xxm) == 64
        @test length(c17.yym) == 64

        # 33x33 grid
        c33 = default_grid_coordinates()
        @test c33.Nx == 33
        @test c33.Ny == 33
        @test c33.Nx1 == 34
        @test c33.Ny1 == 34
        @test length(c33.x) == 33
        @test length(c33.y) == 33
        @test length(c33.xp) == 34
        @test length(c33.yp) == 34
        @test c33.start_marknum == 16384

        # 65x65 grid
        c65 = GridCoordinates(65, 65; xsize=140_000.0, ysize=140_000.0)
        @test c65.Nx == 65
        @test c65.Ny == 65
        @test c65.Nx1 == 66
        @test c65.Ny1 == 66
        @test length(c65.x) == 65
        @test length(c65.xp) == 66
        @test c65.start_marknum == (64 * 4)^2

        # Construction from GridConfig
        cfg_grid = GridConfig(Nx=25, Ny=25, xsize=100_000.0, ysize=100_000.0)
        c25 = GridCoordinates(cfg_grid)
        @test c25.Nx == 25
        @test c25.Ny == 25
        @test c25.xsize ≈ 100_000.0 rtol=1e-12
        @test c25.dx ≈ 100_000.0 / 24 rtol=1e-12
    end

    @testset "Dynamic staggered grid allocation" begin
        c17 = GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        props = Erebus.setup_staggered_grid_properties(c17)
        ETA = props[1]
        RHOX = props[12]
        RHO = props[30]
        @test size(ETA) == (17, 17)
        @test size(RHOX) == (18, 18)
        @test size(RHO) == (18, 18)

        helpers = Erebus.setup_staggered_grid_properties_helpers(c17)
        ETA5 = helpers[1]
        EII = helpers[8]
        @test size(ETA5) == (17, 17)
        @test size(EII) == (18, 18)
    end

    @testset "Dynamic linear system setups and gravitational solve" begin
        c17 = GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        RP, SP = Erebus.setup_gravitational_lse(c17)
        @test length(RP) == 18 * 18
        @test length(SP) == 18 * 18

        R, S = Erebus.setup_hydromechanical_lse(c17)
        @test length(R) == 18 * 18 * 6
        @test length(S) == 18 * 18 * 6

        RT, ST = Erebus.setup_thermal_lse(c17)
        @test length(RT) == 18 * 18
        @test length(ST) == 18 * 18

        RHO = fill(3000.0, 18, 18)
        FI = zeros(18, 18)
        gx = zeros(18, 18)
        gy = zeros(18, 18)
        Erebus.compute_gravity_solution!(SP, RP, RHO, FI, gx, gy; coords=c17)
        @test all(isfinite, FI)
        @test all(isfinite, gx)
        @test all(isfinite, gy)
        @test any(!iszero, gx)
        @test any(!iszero, gy)
    end

    @testset "GridCoordinates degenerate bounds validation" begin
        @test_throws ArgumentError GridCoordinates(2, 17)
        @test_throws ArgumentError GridCoordinates(17, 1)
        @test_throws ArgumentError GridCoordinates(17, 17; xsize=0.0)
        @test_throws ArgumentError GridCoordinates(17, 17; ysize=-10.0)
        @test_throws ArgumentError GridCoordinates(17, 17; Nxmc=0)
        @test_throws ArgumentError GridCoordinates(17, 17; Nymc=-1)
    end

    @testset "End-to-end simulation on non-33x33 grid (17x17)" begin
        mktempdir() do tmpdir
            cfg17 = SimulationConfig(
                grid=GridConfig(Nx=17, Ny=17, xsize=140_000.0, ysize=140_000.0),
                geometry=GeometryConfig(
                    rplanet=50_000.0,
                    rcrust=50_000.0,
                    xcenter=70_000.0,
                    ycenter=70_000.0,
                    psurface=1.0e3,
                ),
                time=TimeConfig(
                    n_steps=2, start_step=1, dt_initial=1.0e10, dt_longest=1.0e11
                ),
                output=OutputConfig(savematstep=1, output_dir=tmpdir),
            )
            Erebus.simulation_loop(cfg17; output_path=tmpdir)
            ckpt_path = joinpath(tmpdir, "output_00002.jld2")
            @test isfile(ckpt_path)
            ckpt = Erebus.load_state(ckpt_path)
            @test ckpt["timestep"] == 2
            @test ckpt["timesum"] > 0.0
            @test ckpt["Nx"] == 17
            @test ckpt["Ny"] == 17

            # Checkpoint dimension mismatch detection
            cfg33 = SimulationConfig(
                grid=GridConfig(Nx=33, Ny=33, xsize=140_000.0, ysize=140_000.0),
                time=TimeConfig(n_steps=3, start_step=1),
                output=OutputConfig(savematstep=1, output_dir=tmpdir),
            )
            @test_throws DimensionMismatch Erebus.simulation_loop(
                cfg33;
                output_path=tmpdir,
                restart_from=joinpath(tmpdir, "output_00002.jld2"),
            )

            # Checkpoint domain size mismatch detection (same Nx, Ny but different xsize, ysize)
            cfg_wrong_size = SimulationConfig(
                grid=GridConfig(Nx=17, Ny=17, xsize=200_000.0, ysize=200_000.0),
                geometry=GeometryConfig(
                    rplanet=50_000.0,
                    rcrust=50_000.0,
                    xcenter=100_000.0,
                    ycenter=100_000.0,
                    psurface=1.0e3,
                ),
                time=TimeConfig(n_steps=3, start_step=1),
                output=OutputConfig(savematstep=1, output_dir=tmpdir),
            )
            @test_throws DimensionMismatch Erebus.simulation_loop(
                cfg_wrong_size;
                output_path=tmpdir,
                restart_from=joinpath(tmpdir, "output_00002.jld2"),
            )
        end
    end

    @testset "Non-square grid coordinates length consistency (Nx != Ny)" begin
        for (nx, ny) in [(33, 17), (17, 33), (65, 33)]
            c = GridCoordinates(nx, ny; xsize=140_000.0, ysize=100_000.0)
            @test length(c.xvx) == c.Nx1
            @test length(c.yvx) == c.Ny1
            @test length(c.xvy) == c.Nx1
            @test length(c.yvy) == c.Ny1
            @test length(c.xp) == c.Nx1
            @test length(c.yp) == c.Ny1
        end
    end

    @testset "Gravitational potential boundary condition on non-square domain" begin
        c = GridCoordinates(33, 17; xsize=140_000.0, ysize=70_000.0)
        RHO = fill(3000.0, c.Ny1, c.Nx1)
        RP = zeros(c.Ny1, c.Nx1)
        LP = assemble_gravitational_lse!(RHO, RP; coords=c)
        total_dofs = c.Nx1 * c.Ny1
        @test size(LP) == (total_dofs, total_dofs)
        # Dirichlet boundary: corner and boundary DOFs must have unit diagonal
        @test LP[1, 1] ≈ 1.0 rtol=1e-12
        @test LP[total_dofs, total_dofs] ≈ 1.0 rtol=1e-12
        @test iszero(RP[1, 1])
    end
end

@testset "Geometry" begin
    @testset "setup_staggered_grid_geometry()" begin
        # verification, from HTM-planetary.m line 38ff
        # Basic nodes
        x_ver=0:dx:xsize # Horizontal coordinates of basic grid points, m
        y_ver=0:dy:ysize # Vertical coordinates of basic grid points, m
        # Vx-Nodes
        xvx_ver=0:dx:(xsize + dy) # Horizontal coordinates of vx grid points, m
        yvx_ver=(-dy / 2):dy:(ysize + dy / 2) # Vertical coordinates of vx grid points, m
        # Vy-nodes
        xvy_ver=(-dx / 2):dx:(xsize + dx / 2) # Horizontal coordinates of vy grid points, m
        yvy_ver=0:dy:(ysize + dy) # Vertical coordinates of vy grid points, m
        # P-Nodes
        xp_ver=(-dx / 2):dx:(xsize + dx / 2) # Horizontal coordinates of P grid points, m
        yp_ver=(-dy / 2):dy:(ysize + dy / 2) # Vertical coordinates of P grid points, m
        # test
        @test x == x_ver
        @test y == y_ver
        @test xvx == xvx_ver
        @test yvx == yvx_ver
        @test xvy == xvy_ver
        @test yvy == yvy_ver
        @test xp == xp_ver
        @test yp == yp_ver
        @test length(x) == Nx
        @test length(y) == Ny
        @test length(xvx) == Nx1
        @test length(yvx) == Ny1
        @test length(xvy) == Nx1
        @test length(yvy) == Ny1
        @test length(xp) == Nx1
        @test length(yp) == Ny1
    end # testset "setup_staggered_grid_geometry()"

    @testset "setup_staggered_grid_properties()" begin
        # set up staggered grid properties
        (ETA, ETA0, GGG, EXY, SXY, SXY0, wyx, COH, TEN, FRI, YNY, RHOX, RHOFX, KX, PHIX, vx, vxf, RX, qxD, gx, RHOY, RHOFY, KY, PHIY, vy, vyf, RY, qyD, gy, RHO, RHOCP, ALPHA, ALPHAF, HR, HA, HS, ETAP, GGGP, EXX, SXX, SXX0, tk1, tk2, DT, DT0, vxp, vyp, vxpf, vypf, pr, pf, ps, pr0, pf0, ps0, ETAPHI, BETTAPHI, PHI, APHI, FI, DMP, DHP, XWS) = Erebus.setup_staggered_grid_properties()
        # verification, from HTM-planetary.m line 51ff, i2visHTM_hydration.m line 691, 692, 1587
        # Basic nodes
        ETA_ver = zeros(Ny, Nx) # Viscoplastic Viscosity, Pa*s
        ETA0_ver = zeros(Ny, Nx) # Viscous Viscosity, Pa*s
        GGG_ver = zeros(Ny, Nx) # Shear modulus, Pa
        EXY_ver = zeros(Ny, Nx) # EPSILONxy, 1/s
        SXY_ver = zeros(Ny, Nx) # SIGMAxy, 1/s
        SXY0_ver = zeros(Ny, Nx) # SIGMA0xy, 1/s
        wyx_ver = zeros(Ny, Nx) # Rotation rate, 1/s
        COH_ver = zeros(Ny, Nx) # Compressive strength, Pa
        TEN_ver = zeros(Ny, Nx) # Tensile strength, Pa
        FRI_ver = zeros(Ny, Nx) # Friction
        YNY_ver = zeros(Ny, Nx) # Plastic yielding mark, 1=yes,0=no
        # Vx-Nodes
        RHOX_ver = zeros(Ny1, Nx1) # Density, kg/m^3
        RHOFX_ver = zeros(Ny1, Nx1) # Fluid Density, kg/m^3
        KX_ver = zeros(Ny1, Nx1) # Thermal conductivity, W/m/K
        PHIX_ver = zeros(Ny1, Nx1) # Porosity
        vx_ver = zeros(Ny1, Nx1) # Solid vx-velocity m/s
        vxf_ver = zeros(Ny1, Nx1) # Fluid vx-velocity m/s
        RX_ver = zeros(Ny1, Nx1) # ETAfluid/Kphi ratio , m^2
        qxD_ver = zeros(Ny1, Nx1) # qx-Darcy flux m/s
        gx_ver = zeros(Ny1, Nx1) # gx-gravity, m/s^2
        # Vy-Nodes
        RHOY_ver = zeros(Ny1, Nx1) # Density, kg/m^3
        RHOFY_ver = zeros(Ny1, Nx1) # Fluid Density, kg/m^3
        KY_ver = zeros(Ny1, Nx1) # Thermal conductivity, W/m/K
        PHIY_ver = zeros(Ny1, Nx1) # Porosity
        vy_ver = zeros(Ny1, Nx1) # Solid vy-velocity m/s
        vyf_ver = zeros(Ny1, Nx1) # Fluid vy-velocity m/s
        RY_ver = zeros(Ny1, Nx1) # ETAfluid/Kphi ratio , m^2
        qyD_ver = zeros(Ny1, Nx1) # qy-Darcy flux m/s
        gy_ver = zeros(Ny1, Nx1) # gy-gravity, m/s^2
        # P-nodes
        RHO_ver = zeros(Ny1, Nx1) # Density, kg/m^3
        RHOCP_ver = zeros(Ny1, Nx1) # Volumetric heat capacity, J/m^3/K
        ALPHA_ver = zeros(Ny1, Nx1) # Thermal expansion, J/m^3/K
        ALPHAF_ver = zeros(Ny1, Nx1) # Fluid Thermal expansion, J/m^3/K
        HR_ver = zeros(Ny1, Nx1) # Radioactive heating, W/m^3
        HA_ver = zeros(Ny1, Nx1) # Adiabatic heating, W/m^3
        HS_ver = zeros(Ny1, Nx1) # Shear heating, W/m^3
        ETAP_ver = zeros(Ny1, Nx1) # Viscosity, Pa*s
        GGGP_ver = zeros(Ny1, Nx1) # Shear modulus, Pa
        # RMK: EXX, SXX (but oddly not SXX0) are first defined with
        # size (Ny, Nx) in the verification code (lines 93-95), but later
        # redefined with size (Ny1, Nx1) (lines 1158, 1160). Possibly an 
        # oversight; we assume the latter size from the beginning.
        EXX_ver = zeros(Ny1, Nx1) # EPSILONxx, 1/s
        SXX_ver = zeros(Ny1, Nx1) # SIGMA'xx, 1/s
        SXX0_ver = zeros(Ny1, Nx1) # SIGMA0'xx, 1/s
        tk1_ver = zeros(Ny1, Nx1) # Old temperature, K
        tk2_ver = zeros(Ny1, Nx1) # New temperature, K
        DT_ver = zeros(Ny1, Nx1) # temperature difference, K
        DT0_ver = zeros(Ny1, Nx1) # previous temperature difference, K
        vxp_ver = zeros(Ny1, Nx1) # Solid Vx in pressure nodes, m/s
        vyp_ver = zeros(Ny1, Nx1) # Solid Vy in pressure nodes, m/s
        vxpf_ver = zeros(Ny1, Nx1) # Fluid Vx in pressure nodes, m/s
        vypf_ver = zeros(Ny1, Nx1) # Fluid Vy in pressure nodes, m/s
        pr_ver = zeros(Ny1, Nx1) # Total Pressure, Pa
        pf_ver = zeros(Ny1, Nx1) # Fluid Pressure, Pa
        ps_ver = zeros(Ny1, Nx1) # Solid Pressure, Pa
        pr0_ver = zeros(Ny1, Nx1) # Old Total Pressure, Pa
        pf0_ver = zeros(Ny1, Nx1) # Old Fluid Pressure, Pa
        ps0_ver = zeros(Ny1, Nx1) # Old Solid Pressure, Pa
        ETAPHI_ver = zeros(Ny1, Nx1) # Bulk Viscosity, Pa*s
        BETTAPHI_ver = zeros(Ny1, Nx1) # Bulk compresibility, Pa*s
        PHI_ver = zeros(Ny1, Nx1) # porosity
        APHI_ver = zeros(Ny1, Nx1) # Dln((1-PHI)/PHI)/Dt
        FI_ver = zeros(Ny1, Nx1) # Gravity potential, J/kg
        DMP_ver = zeros(Ny1, Nx1)
        DHP_ver = zeros(Ny1, Nx1)
        XWS_ver = zeros(Ny1, Nx1)
        # test
        @test ETA == ETA_ver
        @test ETA0 == ETA0_ver
        @test GGG == GGG_ver
        @test EXY == EXY_ver
        @test SXY == SXY_ver
        @test SXY0 == SXY0_ver
        @test wyx == wyx_ver
        @test COH == COH_ver
        @test TEN == TEN_ver
        @test FRI == FRI_ver
        @test YNY == YNY_ver
        @test RHOX == RHOX_ver
        @test RHOFX == RHOFX_ver
        @test KX == KX_ver
        @test PHIX == PHIX_ver
        @test vx == vx_ver
        @test vxf == vxf_ver
        @test RX == RX_ver
        @test qxD == qxD_ver
        @test gx == gx_ver
        @test RHOY == RHOY_ver
        @test RHOFY == RHOFY_ver
        @test KY == KY_ver
        @test PHIY == PHIY_ver
        @test vy == vy_ver
        @test vyf == vyf_ver
        @test RY == RY_ver
        @test qyD == qyD_ver
        @test gy == gy_ver
        @test RHO == RHO_ver
        @test RHOCP == RHOCP_ver
        @test ALPHA == ALPHA_ver
        @test ALPHAF == ALPHAF_ver
        @test HR == HR_ver
        @test HA == HA_ver
        @test HS == HS_ver
        @test ETAP == ETAP_ver
        @test GGGP == GGGP_ver
        @test EXX == EXX_ver
        @test SXX == SXX_ver
        @test SXX0 == SXX0_ver
        @test tk1 == tk1_ver
        @test tk2 == tk2_ver
        @test DT == DT_ver
        @test DT0 == DT0_ver
        @test vxp == vxp_ver
        @test vyp == vyp_ver
        @test vxpf == vxpf_ver
        @test vypf == vypf_ver
        @test pr == pr_ver
        @test pf == pf_ver
        @test ps == ps_ver
        @test pr0 == pr0_ver
        @test pf0 == pf0_ver
        @test ps0 == ps0_ver
        @test ETAPHI == ETAPHI_ver
        @test BETTAPHI == BETTAPHI_ver
        @test PHI == PHI_ver
        @test APHI == APHI_ver
        @test FI == FI_ver
        @test DMP == DMP_ver
        @test DHP == DHP_ver
        @test XWS == XWS_ver
    end # testset "setup_staggered_grid_properties()"

    @testset "setup_staggered_grid_properties_helpers()" begin
        # setup staggered grid properties helpers
        (ETA5, ETA00, YNY5, YNY00, YNY_inv_ETA, DSXY, DSY, EII, SII, DSXX, tk0) = Erebus.setup_staggered_grid_properties_helpers()
        # test
        @test ETA5 == zeros(Float64, Ny, Nx)
        @test ETA00 == zeros(Float64, Ny, Nx)
        @test YNY5 == zeros(Bool, Ny, Nx)
        @test YNY00 == zeros(Bool, Ny, Nx)
        @test YNY_inv_ETA == zeros(Float64, Ny, Nx)
        @test DSXY == zeros(Float64, Ny, Nx)
        @test DSY == zeros(Float64, Ny, Nx)
        @test EII == zeros(Float64, Ny1, Nx1)
        @test SII == zeros(Float64, Ny1, Nx1)
        @test DSXX == zeros(Float64, Ny1, Nx1)
        @test tk0 == zeros(Float64, Ny1, Nx1)
    end # testset "setup_staggered_grid_properties_helpers()"

    @testset "grid_vector()" begin
        grid = rand(rgen, 10, 10)
        @test Erebus.grid_vector(1, 1, grid) ==
            @SVector [grid[1, 1], grid[2, 1], grid[1, 2], grid[2, 2]]
    end # testset "grid_vector()"

    @testset "grid_average()" begin
        grid = rand(rgen, 10, 10)
        @test Erebus.grid_average(1, 1, grid) ==
            0.25 * (grid[1, 1] + grid[2, 1] + grid[1, 2] + grid[2, 2])
    end # testset "grid_average()"

    @testset "apply_insulating_boundary_conditions!()" begin
        max_size = 10
        for j in 3:1:max_size, i in 3:1:max_size
            t = rand(rgen, i, j)
            Erebus.apply_insulating_boundary_conditions!(t)
            @test t[1, 2:(j - 1)] == t[2, 2:(j - 1)]
            @test t[i, 2:(j - 1)] == t[i - 1, 2:(j - 1)]
            @test t[:, 1] == t[:, 2]
            @test t[:, j] == t[:, j - 1]
        end
    end # testset "apply_insulating_boundary_conditions!()"
end

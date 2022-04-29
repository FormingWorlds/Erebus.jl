using ExtendableSparse
using HydrologyPlanetesimals
using Parameters
using StaticArrays
using Test

@testset verbose=true "HydrologyPlanetesimals.jl" begin

    @testset "setup_dynamic_simulation_parameters()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        # set up dynamic simulation parameters
        timestep,
        dt,
        timesum,
        marknum,
        hrsolidm,
        hrfluidm,
        YERRNOD = HydrologyPlanetesimals.setup_dynamic_simulation_parameters(sp)
        # verification & test
        @test timestep == sp.start_step
        @test dt == sp.dtelastic
        @test timesum == sp.start_time
        @test marknum == sp.start_marknum
        @test hrsolidm == sp.start_hrsolidm
        @test hrfluidm == sp.start_hrfluidm
        @test YERRNOD == zeros(sp.nplast)
    end # testset "setup_dynamic_simulation_parameters()"

    @testset "setup_staggered_grid_geometry()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        dx, dy = sp.dx, sp.dy
        xsize, ysize = sp.xsize, sp.ysize
        # set up staggered grid geometry
        x = sp.x
        y = sp.y
        xvx = sp.xvx
        yvx = sp.xvy
        xvy = sp.xvy
        yvy = sp.yvy
        xp = sp.xp
        yp = sp.yp
        # verification, from madcph.m line 38ff
        # Basic nodes
        x_ver=0:dx:xsize; # Horizontal coordinates of basic grid points, m
        y_ver=0:dy:ysize; # Vertical coordinates of basic grid points, m
        # Vx-Nodes
        xvx_ver=0:dx:xsize+dy; # Horizontal coordinates of vx grid points, m
        yvx_ver=-dy/2:dy:ysize+dy/2; # Vertical coordinates of vx grid points, m
        # Vy-nodes
        xvy_ver=-dx/2:dx:xsize+dx/2; # Horizontal coordinates of vy grid points, m
        yvy_ver=0:dy:ysize+dy; # Vertical coordinates of vy grid points, m
        # P-Nodes
        xp_ver=-dx/2:dx:xsize+dx/2; # Horizontal coordinates of P grid points, m
        yp_ver=-dy/2:dy:ysize+dy/2; # Vertical coordinates of P grid points, m
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
        sp = HydrologyPlanetesimals.StaticParameters() 
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        # set up staggered grid properties
        (
            ETA,
            ETA0,
            GGG,
            EXY,
            SXY,
            SXY0,
            wyx,
            COH,
            TEN,
            FRI,
            YNY,
            RHOX,
            RHOFX,
            KX,
            PHIX,
            vx,
            vxf,
            RX,
            qxD,
            gx,
            RHOY,
            RHOFY,
            KY,
            PHIY,
            vy,
            vyf,
            RY,
            qyD,
            gy,
            RHO,
            RHOCP,
            ALPHA,
            ALPHAF,
            HR,
            HA,
            HS,
            ETAP,
            GGGP,
            EXX,
            SXX,
            SXX0,
            tk1,
            tk2,
            vxp,
            vyp,
            vxpf,
            vypf,
            pr,
            pf,
            ps,
            pr0,
            pf0,
            ps0,
            ETAPHI,
            BETTAPHI,
            PHI,
            APHI,
            FI
        ) = HydrologyPlanetesimals.setup_staggered_grid_properties(sp)
        # verification, from madcph.m line 51ff
        # Basic nodes
        ETA_ver = zeros(Ny,Nx) # Viscoplastic Viscosity, Pa*s
        ETA0_ver = zeros(Ny,Nx) # Viscous Viscosity, Pa*s
        GGG_ver = zeros(Ny,Nx) # Shear modulus, Pa
        EXY_ver = zeros(Ny,Nx) # EPSILONxy, 1/s
        SXY_ver = zeros(Ny,Nx) # SIGMAxy, 1/s
        SXY0_ver = zeros(Ny,Nx) # SIGMA0xy, 1/s
        wyx_ver = zeros(Ny,Nx) # Rotation rate, 1/s
        COH_ver = zeros(Ny,Nx) # Compressive strength, Pa
        TEN_ver = zeros(Ny,Nx) # Tensile strength, Pa
        FRI_ver = zeros(Ny,Nx) # Friction
        YNY_ver = zeros(Ny,Nx) # Plastic yielding mark, 1=yes,0=no
        # Vx-Nodes
        RHOX_ver = zeros(Ny1,Nx1) # Density, kg/m^3
        RHOFX_ver = zeros(Ny1,Nx1) # Fluid Density, kg/m^3
        KX_ver = zeros(Ny1,Nx1) # Thermal conductivity, W/m/K
        PHIX_ver = zeros(Ny1,Nx1) # Porosity
        vx_ver = zeros(Ny1,Nx1) # Solid vx-velocity m/s
        vxf_ver = zeros(Ny1,Nx1) # Fluid vx-velocity m/s
        RX_ver = zeros(Ny1,Nx1) # ETAfluid/Kphi ratio , m^2
        qxD_ver = zeros(Ny1,Nx1) # qx-Darcy flux m/s
        gx_ver = zeros(Ny1,Nx1) # gx-gravity, m/s^2
        # Vy-Nodes
        RHOY_ver = zeros(Ny1,Nx1) # Density, kg/m^3
        RHOFY_ver = zeros(Ny1,Nx1) # Fluid Density, kg/m^3
        KY_ver = zeros(Ny1,Nx1) # Thermal conductivity, W/m/K
        PHIY_ver = zeros(Ny1,Nx1) # Porosity
        vy_ver = zeros(Ny1,Nx1) # Solid vy-velocity m/s
        vyf_ver = zeros(Ny1,Nx1) # Fluid vy-velocity m/s
        RY_ver = zeros(Ny1,Nx1) # ETAfluid/Kphi ratio , m^2
        qyD_ver = zeros(Ny1,Nx1) # qy-Darcy flux m/s
        gy_ver = zeros(Ny1,Nx1) # gy-gravity, m/s^2
        # P-nodes
        RHO_ver = zeros(Ny1,Nx1) # Density, kg/m^3
        RHOCP_ver = zeros(Ny1,Nx1) # Volumetric heat capacity, J/m^3/K
        ALPHA_ver = zeros(Ny1,Nx1) # Thermal expansion, J/m^3/K
        ALPHAF_ver = zeros(Ny1,Nx1) # Fluid Thermal expansion, J/m^3/K
        HR_ver = zeros(Ny1,Nx1) # Radioactive heating, W/m^3
        HA_ver = zeros(Ny1,Nx1) # Adiabatic heating, W/m^3
        HS_ver = zeros(Ny1,Nx1) # Shear heating, W/m^3
        ETAP_ver = zeros(Ny1,Nx1) # Viscosity, Pa*s
        GGGP_ver = zeros(Ny1,Nx1) # Shear modulus, Pa
        # RMK: EXX, SXX (but oddly not SXX0) are first defined with
        # size (Ny, Nx) in the verification code (lines 93-95), but later
        # redefined with size (Ny1, Nx1) (lines 1158, 1160). Possibly an 
        # oversight; we assume the latter size from the beginning.
        EXX_ver = zeros(Ny1,Nx1) # EPSILONxx, 1/s
        SXX_ver = zeros(Ny1,Nx1) # SIGMA'xx, 1/s
        SXX0_ver = zeros(Ny1,Nx1) # SIGMA0'xx, 1/s
        tk1_ver = zeros(Ny1,Nx1) # Old temperature, K
        tk2_ver = zeros(Ny1,Nx1) # New temperature, K
        vxp_ver = zeros(Ny1,Nx1) # Solid Vx in pressure nodes, m/s
        vyp_ver = zeros(Ny1,Nx1) # Solid Vy in pressure nodes, m/s
        vxpf_ver = zeros(Ny1,Nx1) # Fluid Vx in pressure nodes, m/s
        vypf_ver = zeros(Ny1,Nx1) # Fluid Vy in pressure nodes, m/s
        pr_ver = zeros(Ny1,Nx1) # Total Pressure, Pa
        pf_ver = zeros(Ny1,Nx1) # Fluid Pressure, Pa
        ps_ver = zeros(Ny1,Nx1) # Solid Pressure, Pa
        pr0_ver = zeros(Ny1,Nx1) # Old Total Pressure, Pa
        pf0_ver = zeros(Ny1,Nx1) # Old Fluid Pressure, Pa
        ps0_ver = zeros(Ny1,Nx1) # Old Solid Pressure, Pa
        ETAPHI_ver = zeros(Ny1,Nx1) # Bulk Viscosity, Pa*s
        BETTAPHI_ver = zeros(Ny1,Nx1) # Bulk compresibility, Pa*s
        PHI_ver = zeros(Ny1,Nx1) # porosity
        APHI_ver = zeros(Ny1,Nx1) # Dln((1-PHI)/PHI)/Dt
        FI_ver = zeros(Ny1,Nx1) # Gravity potential, J/kg
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
    end # testset "setup_staggered_grid_properties()"

    @testset "setup_staggered_grid_properties_helpers()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        # setup staggered grid properties helpers
        (
            ETA5,
            ETA00,
            YNY5,
            YNY00,
            YNY_inv_ETA,
            DSXY,
            ETAcomp,
            SXYcomp,
            dRHOXdx,
            dRHOXdy,
            dRHOYdx,
            dRHOYdy,
            ETAPcomp,
            SXXcomp,
            SYYcomp,
            EII,
            SII,
            DSXX,
            tk0
        ) = HydrologyPlanetesimals.setup_staggered_grid_properties_helpers(sp)
        # test
        @test ETA5 == zeros(Float64, Ny, Nx)
        @test ETA00 == zeros(Float64, Ny, Nx)
        @test YNY5 == zeros(Bool, Ny, Nx)
        @test YNY00 == zeros(Bool, Ny, Nx)
        @test YNY_inv_ETA == zeros(Float64, Ny, Nx)
        @test DSXY == zeros(Float64, Ny, Nx)
        @test ETAcomp == zeros(Float64, Ny, Nx)
        @test SXYcomp == zeros(Float64, Ny, Nx)
        @test dRHOXdx == zeros(Float64, Ny1, Nx1)
        @test dRHOXdy == zeros(Float64, Ny1, Nx1)
        @test dRHOYdx == zeros(Float64, Ny1, Nx1)
        @test dRHOYdy == zeros(Float64, Ny1, Nx1)
        @test ETAPcomp == zeros(Float64, Ny1, Nx1)
        @test SXXcomp == zeros(Float64, Ny1, Nx1)
        @test SYYcomp == zeros(Float64, Ny1, Nx1)
        @test EII == zeros(Float64, Ny1, Nx1)
        @test SII == zeros(Float64, Ny1, Nx1)
        @test DSXX == zeros(Float64, Ny1, Nx1)
        @test tk0 == zeros(Float64, Ny1, Nx1)
    end # testset "setup_staggered_grid_properties_helpers()"

    @testset "setup_interpolated_properties()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        # setup interpolated grid properties
        (
            ETA0SUM,
            ETASUM,
            GGGSUM,
            SXYSUM,
            COHSUM,
            TENSUM,
            FRISUM,
            WTSUM,
            RHOXSUM,
            RHOFXSUM,
            KXSUM,
            PHIXSUM,
            RXSUM,
            WTXSUM,
            RHOYSUM,
            RHOFYSUM,
            KYSUM,
            PHIYSUM,
            RYSUM,
            WTYSUM,
            RHOSUM,
            RHOCPSUM,
            ALPHASUM,
            ALPHAFSUM,
            HRSUM,
            GGGPSUM,
            SXXSUM,
            TKSUM,
            PHISUM,
            WTPSUM
        ) = HydrologyPlanetesimals.setup_interpolated_properties(sp)
        # verification and test
        @test ETA0SUM == zeros(Float64, Ny, Nx, Base.Threads.nthreads())
        @test ETASUM == zeros(Float64, Ny, Nx, Base.Threads.nthreads())
        @test GGGSUM == zeros(Float64, Ny, Nx, Base.Threads.nthreads())
        @test SXYSUM == zeros(Float64, Ny, Nx, Base.Threads.nthreads())
        @test COHSUM == zeros(Float64, Ny, Nx, Base.Threads.nthreads())
        @test TENSUM == zeros(Float64, Ny, Nx, Base.Threads.nthreads())
        @test FRISUM == zeros(Float64, Ny, Nx, Base.Threads.nthreads())
        @test WTSUM == zeros(Float64, Ny, Nx, Base.Threads.nthreads())
        @test RHOXSUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test RHOFXSUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test KXSUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test PHIXSUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test RXSUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test WTXSUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test RHOYSUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test RHOFYSUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test KYSUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test PHIYSUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test RYSUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test WTYSUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test RHOSUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test RHOCPSUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test ALPHASUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test ALPHAFSUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test HRSUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test GGGPSUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test SXXSUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test TKSUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test PHISUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
        @test WTPSUM == zeros(Float64, Ny1, Nx1, Base.Threads.nthreads())
    end # testset "setup_interpolated_properties()"

    @testset "setup_marker_properties()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        xsize, ysize = sp.xsize, sp.ysize
        Nxmc, Nymc = sp.Nxmc, sp.Nymc
        Nxm, Nym = sp.Nxm, sp.Nym
        dxm, dym = sp.dxm, sp.dym
        marknum = sp.start_marknum
        # setup marker properties
        (
            xm,
            ym,
            tm,
            tkm,
            sxxm,
            sxym,
            etavpm,
            phim
        ) = HydrologyPlanetesimals.setup_marker_properties(sp)
        # verification, from madcph.m line 115ff
        Nxmc_ver = 4; # Number of markers per cell in horizontal direction
        Nymc_ver = 4; # Number of markers per cell in vertical direction
        Nxm_ver = (Nx-1)*Nxmc; # Marker grid resolution in horizontal direction
        Nym_ver = (Ny-1)*Nymc; # Marker grid resolution in vertical direction
        dxm_ver = xsize/Nxm; # Marker grid step in horizontal direction,m
        dym_ver = ysize/Nym; # Marker grid step in vertical direction,m
        marknum_ver = Nxm*Nym; # Number of markers
        xm_ver = zeros(marknum); # Horizontal coordinates, m
        ym_ver = zeros(marknum); # Vertical coordinates, m
        tm_ver = zeros(Int, marknum); # Material type
        tkm_ver = zeros(marknum); # Marker temperature, K
        sxxm_ver = zeros(marknum); # SIGMA'xx, Pa
        sxym_ver = zeros(marknum); # SIGMAxy, Pa
        etavpm_ver = zeros(marknum); # Visco-plastic viscosity, Pa
        phim_ver = zeros(marknum); # Marker porosity
        # test
        @test Nxmc == Nxmc_ver
        @test Nymc == Nymc_ver
        @test Nxm == Nxm_ver
        @test Nym == Nym_ver
        @test dxm == dxm_ver
        @test dym == dym_ver
        @test marknum == marknum_ver
        @test xm == xm_ver
        @test ym == ym_ver
        @test tm == tm_ver
        @test tkm == tkm_ver
        @test sxxm == sxxm_ver
        @test sxym == sxym_ver
        @test etavpm == etavpm_ver
        @test phim == phim_ver
    end # testset "setup_marker_properties()"

    @testset "setup_marker_properties_helpers()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        marknum = sp.start_marknum
        # setup marker properties
        (
            rhototalm,
            rhocptotalm,
            etasolidcur,
            etafluidcur,
            etatotalm,
            hrtotalm,
            ktotalm,
            tkm_rhocptotalm,
            etafluidcur_inv_kphim,
            inv_gggtotalm,
            fricttotalm,
            cohestotalm,
            tenstotalm,
            rhofluidcur,
            alphasolidcur,
            alphafluidcur
        ) = HydrologyPlanetesimals.setup_marker_properties_helpers(sp)
        # test
        @test rhocptotalm == zeros(Float64, marknum)
        @test rhocptotalm == zeros(Float64, marknum)
        @test etasolidcur == zeros(Float64, marknum)
        @test etafluidcur == zeros(Float64, marknum)
        @test etatotalm == zeros(Float64, marknum)
        @test hrtotalm == zeros(Float64, marknum)
        @test ktotalm == zeros(Float64, marknum)
        @test tkm_rhocptotalm == zeros(Float64, marknum)
        @test etafluidcur_inv_kphim == zeros(Float64, marknum)
        @test inv_gggtotalm == zeros(Float64, marknum)
        @test fricttotalm == zeros(Float64, marknum)
        @test cohestotalm == zeros(Float64, marknum)
        @test tenstotalm == zeros(Float64, marknum)
        @test rhofluidcur == zeros(Float64, marknum)
        @test alphasolidcur == zeros(Float64, marknum)
        @test alphafluidcur == zeros(Float64, marknum)
    end # testset "setup_marker_properties_helpers()"
    
    @testset "define_markers!() & compute_marker_properties!()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nxm, Nym = sp.Nxm, sp.Nym
        dxm, dym = sp.dxm, sp.dym
        xsize, ysize = sp.xsize, sp.ysize
        rplanet, rcrust = sp.rplanet, sp.rcrust
        etasolidm = sp.etasolidm
        psurface = sp.psurface
        phim0, phimin = sp.phim0, sp.phimin
        marknum = sp.start_marknum
        start_hrsolidm, start_hrfluidm = sp.start_hrsolidm, sp.start_hrfluidm
        hrsolidm, hrfluidm = sp.start_hrsolidm, sp.start_hrfluidm
        gggsolidm = sp.gggsolidm
        frictsolidm = sp.frictsolidm
        cohessolidm = sp.cohessolidm
        tenssolidm = sp.tenssolidm
        alphasolidm = sp.alphasolidm
        alphafluidm = sp.alphafluidm
        rhosolidm = sp.rhosolidm
        rhofluidm = sp.rhofluidm
        rhocpsolidm = sp.rhocpsolidm
        rhocpfluidm = sp.rhocpfluidm
        tmsilicate = sp.tmsilicate
        tmiron = sp.tmiron
        etamin = sp.etamin
        etasolidm = sp.etasolidm
        etasolidmm = sp.etasolidmm
        etafluidm = sp.etafluidm
        etafluidmm = sp.etafluidmm
        ksolidm = sp.ksolidm
        kfluidm = sp.kfluidm
        kphim0 = sp.kphim0
        phim0 = sp.phim0
        (
            xm,
            ym,
            tm,
            tkm,
            sxxm,
            sxym,
            etavpm,
            phim
        ) = HydrologyPlanetesimals.setup_marker_properties(sp)
        (
            rhototalm,
            rhocptotalm,
            etasolidcur,
            etafluidcur,
            etatotalm,
            hrtotalm,
            ktotalm,
            tkm_rhocptotalm,
            etafluidcur_inv_kphim,
            inv_gggtotalm,
            fricttotalm,
            cohestotalm,
            tenstotalm,
            rhofluidcur,
            alphasolidcur,
            alphafluidcur
        ) = HydrologyPlanetesimals.setup_marker_properties_helpers(sp)
        xm_ver = zeros(marknum)
        ym_ver = zeros(marknum)
        tm_ver = zeros(Int, marknum)
        tkm_ver = zeros(marknum)
        phim_ver = zeros(marknum)
        etavpm_ver = zeros(marknum)
        # define markers
        HydrologyPlanetesimals.define_markers!(
            xm,
            ym,
            tm,
            phim,
            etavpm,
            rhototalm,
            rhocptotalm,
            etatotalm,
            hrtotalm,
            ktotalm,
            etafluidcur,
            tkm,
            inv_gggtotalm,
            fricttotalm,
            cohestotalm,
            tenstotalm,
            rhofluidcur,
            alphasolidcur,
            alphafluidcur,
            sp,
            randomized=false
        )
        for m=1:1:marknum
            HydrologyPlanetesimals.compute_marker_properties!(
                m,
                tm,
                tkm,
                rhototalm,
                rhocptotalm,
                etasolidcur,
                etafluidcur,
                etatotalm,
                hrtotalm,
                ktotalm,
                tkm_rhocptotalm,
                etafluidcur_inv_kphim,
                hrsolidm,
                hrfluidm,
                phim,
                sp
            )
        end
        # verification, from madcph.m, line 180ff
        m=1; # Marker counter
        for jm=1:1:Nxm
            for im=1:1:Nym
                # Define marker coordinates
                xm_ver[m]=dxm/2+(jm-1)*dxm # +(rand-0.5)*dxm;
                ym_ver[m]=dym/2+(im-1)*dym # +(rand-0.5)*dym;
                # Marker properties
                rmark=((xm_ver[m]-xsize/2)^2+(ym_ver[m]-ysize/2)^2)^0.5;
                if(rmark<rplanet)
                    # Planet
                    tm_ver[m]=1; # mantle
                    if(rmark>rcrust) 
                        tm_ver[m]=2; # crust
                    end
                    tkm_ver[m]=300; # Temperature
                    phim_ver[m]=phim0 # *(1+1.0*(rand-0.5)); # Porosity
                    etavpm_ver[m]=etasolidm[tm_ver[m]];#*exp(-28*phim_ver[m]); % Matrix viscosity
                else
                    # Sticky space (to have internal free surface)
                    tm_ver[m]=3; # Material type
                    tkm_ver[m]=273; # Temperature
                    phim_ver[m]=phimin; # Porosity
                    etavpm_ver[m]=etasolidm[tm_ver[m]]; # Matrix viscosity
                end
                # Update marker counter
                m=m+1;
            end
        end
        # test
        for m=1:1:marknum
            # define_markers()!
            @test xm[m] == xm_ver[m]
            @test ym[m] == ym_ver[m]
            @test tm[m] == tm_ver[m]
            @test tkm[m] == tkm_ver[m]
            @test phim[m] == phim_ver[m]
            @test etavpm[m] == etavpm_ver[m]
            # compute_marker_properties()!
            # for air type markers
            if tm[m] == 3
                @test rhototalm[m] == rhosolidm[tm[m]]
                @test rhocptotalm[m] == rhocpsolidm[tm[m]]
                @test etatotalm[m] == etasolidm[tm[m]]
                @test hrtotalm[m] == start_hrsolidm[tm[m]]
                @test ktotalm[m] == ksolidm[tm[m]]           
                @test etafluidcur[m] == etafluidm[tm[m]]
            # for rock type markers
            elseif tm[m] < 3
                @test rhototalm[m] == HydrologyPlanetesimals.total(
                    rhosolidm[tm[m]], rhofluidm[tm[m]], phim[m])
                @test rhocptotalm[m] == HydrologyPlanetesimals.total(
                    rhocpsolidm[tm[m]], rhocpfluidm[tm[m]], phim[m])
                @test etasolidcur[m] == ifelse(
                    tkm[m]>tmsilicate, etasolidmm[tm[m]], etasolidm[tm[m]])
                @test etafluidcur[m] == ifelse(
                    tkm[m]>tmiron, etafluidmm[tm[m]], etafluidm[tm[m]])
                @test etatotalm[m] == max(
                    etamin, etasolidcur[m], etafluidcur[m])
                @test hrtotalm[m] == HydrologyPlanetesimals.total(
                    hrsolidm[tm[m]], hrfluidm[tm[m]], phim[m])
                @test ktotalm[m] == HydrologyPlanetesimals.ktotal(
                    ksolidm[tm[m]], kfluidm[tm[m]], phim[m])
            end
            # for all markers
            @test inv_gggtotalm[m] == inv(gggsolidm[tm[m]])
            @test fricttotalm[m] == frictsolidm[tm[m]]
            @test cohestotalm[m] == cohessolidm[tm[m]]
            @test tenstotalm[m] == tenssolidm[tm[m]]
            @test rhofluidcur[m] == rhofluidm[tm[m]]
            # @test alphasolidcur[m] == alphasolidm[tm[m]]
            @test alphafluidcur[m] == alphafluidm[tm[m]]
            @test tkm_rhocptotalm[m] == tkm[m] * rhocptotalm[m]
            @test etafluidcur_inv_kphim[m] == (
                etafluidcur[m] / HydrologyPlanetesimals.kphi(
                    kphim0[tm[m]], phim0, phim[m])
            )
        end
    end # testset "define_markers!() & compute_marker_properties!()"

    @testset "update_marker_viscosity!()" begin
        sp = HydrologyPlanetesimals.StaticParameters(Nxmc=1, Nymc=1)
        Nx, Ny = sp.Nx, sp.Ny
        dx, dy = sp.dx, sp.dy
        x, y = sp.x, sp.y
        etasolidm, etasolidmm = sp.etasolidm, sp.etasolidmm
        etafluidm, etafluidmm = sp.etafluidm, sp.etafluidmm
        tmsilicate, tmiron = sp.tmsilicate, sp.tmiron
        marknum = sp.start_marknum
        xm = rand(-dx:0.1:x[end]+dx, marknum)
        ym = rand(-dy:0.1:y[end]+dy, marknum)
        tm = rand(1:3, marknum)
        tkm = rand(sp.tmiron-100:0.1:sp.tmsilicate+100, marknum)
        etatotalm = [etasolidm[tm[m]] for m in 1:1:marknum]
        ETA = rand(Ny, Nx)
        YNY = rand(Bool, Ny, Nx)
        YNY_inv_ETA = YNY ./ ETA
        etavpm = zeros(marknum)
        etavpm_ver = zero(etavpm)
        # update marker Viscosity
        for m=1:1:marknum
            HydrologyPlanetesimals.update_marker_viscosity!(
                m, xm, ym, tm, tkm, etatotalm, etavpm, YNY, YNY_inv_ETA, sp)
        end
        # verification, from madcph.m, line 1321ff
        for m=1:1:marknum
            # Interpolation viscosity from basic nodes
            # Define i;j indexes for the upper left node
            j=trunc(Int, (xm[m]-x[1])/dx)+1
            i=trunc(Int, (ym[m]-y[1])/dy)+1
            if j<1
                j=1
            elseif j>Nx-1
                j=Nx-1
            end
            if i<1
                i=1
            elseif i>Ny-1
                i=Ny-1
            end
            # Compute distances
            dxmj=xm[m]-x[j]
            dymi=ym[m]-y[i]
            # Compute weights
            wtmij=(1-dxmj/dx)*(1-dymi/dy)
            wtmi1j=(1-dxmj/dx)*(dymi/dy);    
            wtmij1=(dxmj/dx)*(1-dymi/dy)
            wtmi1j1=(dxmj/dx)*(dymi/dy)
            # Matrix viscosity
            if tm[m]<3
                # Rocks
                etasolidcur_ver=etasolidm[tm[m]]
                if tkm[m]>tmsilicate
                    etasolidcur_ver=etasolidmm[tm[m]]
                end
                etatotalm_ver=etasolidcur_ver;#*exp(-28*phim[m])
            else
                # Sticky air
                etatotalm_ver=etasolidm[tm[m]]
            end
            if YNY[i,j]>0 || YNY[i+1,j]>0 || YNY[i,j+1]>0 || YNY[i+1,j+1]>0
        #         etavpm[m]=ETA[i,j]*wtmij+ETA[i+1,j]*wtmi1j+...
        #                 ETA[i,j+1]*wtmij1+ETA[i+1,j+1]*wtmi1j1
        #         etavpm[m]=1/(1/ETA[i,j]*wtmij+1/ETA[i+1,j]*wtmi1j+...
        #                 1/ETA[i,j+1]*wtmij1+1/ETA[i+1,j+1]*wtmi1j1)
                etavpm_ver[m]=1/(YNY[i,j]/ETA[i,j]*wtmij+YNY[i+1,j]/ETA[i+1,j]*wtmi1j+ YNY[i,j+1]/ETA[i,j+1]*wtmij1+YNY[i+1,j+1]/ETA[i+1,j+1]*wtmi1j1)
                if etavpm_ver[m]>=etatotalm_ver
                    etavpm_ver[m]=etatotalm_ver
                end
            else
                etavpm_ver[m]=etatotalm_ver
            end
        end
        # test
        for m=1:1:marknum
            @test etavpm[m] ≈ etavpm_ver[m] atol=1e-6
        end
    end # testset "update_marker_viscosity!()"

    @testset "distance()" begin
        @test HydrologyPlanetesimals.distance(0, 0, 0, 0) == 0
        @test HydrologyPlanetesimals.distance(1, 0, 0, 0) == 1
        @test HydrologyPlanetesimals.distance(0, 1, 0, 0) == 1
        @test HydrologyPlanetesimals.distance(0, 0, 1, 0) == 1
        @test HydrologyPlanetesimals.distance(0, 0, 0, 1) == 1
        @test HydrologyPlanetesimals.distance(0, 0, 1, 1) ≈ sqrt(2)
        @test HydrologyPlanetesimals.distance(1, 1, 0, 0) ≈ sqrt(2)
        @test HydrologyPlanetesimals.distance(-1, -1, 1, 1) ≈ sqrt(8)
        # random tests
        num_samples = 100
        x = rand(0:0.001:1000, 2, num_samples)
        y = rand(0:0.001:1000, 2, num_samples)
        for i in 1:num_samples
            @test HydrologyPlanetesimals.distance(
                x[:, i]..., y[:, i]...) ≈ sqrt(sum((x[:, i] - y[:, i]).^2))
        end
    end # testset "distance()"

    @testset "total()" begin
        @test HydrologyPlanetesimals.total(0, 0, 0) == 0
        @test HydrologyPlanetesimals.total(1, 0, 0) == 1
        @test HydrologyPlanetesimals.total(0, 1, 0) == 0
        @test HydrologyPlanetesimals.total(0, 0, 1) == 0
        @test HydrologyPlanetesimals.total(1, 2, 0.5) == 1.5
    end # testset "total()"

    @testset "grid_vector()" begin
        grid = rand(10, 10)
        @test HydrologyPlanetesimals.grid_vector(1, 1, grid) == @SVector [
            grid[1, 1], grid[2, 1], grid[1, 2], grid[2, 2]
        ]
    end # testset "grid_vector()"

    @testset "ktotal()" begin
        # verification, from madcph.m, line 1761
        ktotalm_ver(ksolidm, kfluid, phim) = (
            ksolidm*kfluid/2+(
                (ksolidm*(3*phim-2)+kfluid*(1-3*phim))^2
                )/16
            )^0.5 - (ksolidm*(3*phim-2)+ kfluid*(1-3*phim))/4
        @test HydrologyPlanetesimals.ktotal(1., 2., 3.) == ktotalm_ver(
            1., 2., 3.)
    end # testset "ktotal()"

    @testset "kphi()" begin
        # verification, from madcph.m, line 333
        kphim(kphim0, phim, phim0)=kphim0*(phim/phim0)^3/((1-phim)/(1-phim0))^2
        @test HydrologyPlanetesimals.kphi(1., 2., 3.) == kphim(1., 2., 3.)
    end # testset "kphi()"

    @testset "etatotal_rocks()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        etasolidm = sp.etasolidm
        etasolidmm = sp.etasolidmm
        etafluidm = sp.etafluidm
        etafluidmm = sp.etafluidmm
        etamin = sp.etamin
        tmin = min(sp.tmiron, sp.tmsilicate) - 10
        tmid = min(sp.tmiron, sp.tmsilicate) + 0.5*abs(sp.tmiron-sp.tmsilicate)
        tmax = max(sp.tmiron, sp.tmsilicate) + 10
        for type in 1:2
            @test HydrologyPlanetesimals.etatotal_rocks(tmin, type, sp) == max(
                etamin, etasolidm[type], etafluidm[type])
            if sp.tmiron <= sp.tmsilicate
                @test HydrologyPlanetesimals.etatotal_rocks(tmid, type, sp) ==
                max(etamin, etasolidm[type], etafluidmm[type])
            else
                @test HydrologyPlanetesimals.etatotal_rocks(tmid, type, sp) ==
                max(etamin, etasolidmm[type], etafluidm[type])
            end
            @test HydrologyPlanetesimals.etatotal_rocks(tmax, type, sp) == max(
                etamin, etasolidmm[type], etafluidmm[type])
        end
    end # testset "etatotal_rocks()"

    @testset "Q_radiogenic()" begin
        # verification, from madcph.m, line 276
        Q(f, ratio, E, tau, timesum)=f*ratio*E*exp(-timesum/tau)/tau
        @test HydrologyPlanetesimals.Q_radiogenic(1., 2., 3., 4., 5.) == Q(
            1., 2., 3., 4., 5.)
        @test HydrologyPlanetesimals.Q_radiogenic(1., 2., 3., 4., 0.) == Q(
            1., 2., 3., 4., 0.)
    end # testset "Q_radiogenic()"

    @testset "calculate_radioactive_heating()" begin
        hr_al = false
        hr_fe = false
        v = @SVector [0., 0., 0.]
        sp = HydrologyPlanetesimals.StaticParameters(hr_al=hr_al, hr_fe=hr_fe)
        @test HydrologyPlanetesimals.calculate_radioactive_heating(
            1000., sp) == (v, v)

        hr_al = true
        hr_fe = false        
        sp2 = HydrologyPlanetesimals.StaticParameters(hr_al=hr_al, hr_fe=hr_fe)
        Q_al = HydrologyPlanetesimals.Q_radiogenic(
            sp2.f_al, sp2.ratio_al, sp2.E_al, sp2.tau_al, 1000.)
        u = @SVector [Q_al * sp2.rhosolidm[1], Q_al * sp2.rhosolidm[2], 0.]
        @test HydrologyPlanetesimals.calculate_radioactive_heating(
            1000., sp2) == (u, v)

        hr_al = false
        hr_fe = true
        sp3 = HydrologyPlanetesimals.StaticParameters(hr_al=hr_al, hr_fe=hr_fe)
        Q_fe = HydrologyPlanetesimals.Q_radiogenic(
            sp3.f_fe, sp3.ratio_fe, sp3.E_fe, sp3.tau_fe, 1000.)
        w = @SVector [Q_fe * sp2.rhofluidm[1], 0., 0.]
        @test HydrologyPlanetesimals.calculate_radioactive_heating(
            1000., sp3) == (v, w)
    end # testset "calculate_radioactive_heating()"

    @testset "fix_weights() elementary" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        dx, dy = sp.dx, sp.dy
        xsize, ysize = sp.xsize, sp.ysize
        # verification, from madcph.m, line 38ff
        # Basic nodes
        x=0:dx:xsize
        y=0:dy:ysize
        # Vx-Nodes
        xvx=0:dx:xsize+dy
        yvx=-dy/2:dy:ysize+dy/2
        # Vy-nodes
        xvy=-dx/2:dx:xsize+dx/2
        yvy=0:dy:ysize+dy
        # P-Nodes
        xp=-dx/2:dx:xsize+dx/2
        yp=-dy/2:dy:ysize+dy/2

        @testset "basic nodes" begin
        # verification, from madcph.m, line 373ff
        jmin, jmax = sp.jmin_basic, sp.jmax_basic
        imin, imax = sp.imin_basic, sp.imax_basic
        function fix_basic(xm, ym, x_axis, y_axis, dx, dy)
            j=trunc(Int, (xm-x_axis[1])/dx)+1;
            i=trunc(Int, (ym-y_axis[1])/dy)+1;
            if j<1
                j=1
            elseif j>Nx-1 
                j=Nx-1
            end
            if i<1 
                i=1
            elseif i>Ny-1
                i=Ny-1
            end
            dxmj=xm-x[j];
            dymi=ym-y[i];
            wtmij=(1-dxmj/dx)*(1-dymi/dy);
            wtmi1j=(1-dxmj/dx)*(dymi/dy);    
            wtmij1=(dxmj/dx)*(1-dymi/dy);
            wtmi1j1=(dxmj/dx)*(dymi/dy);
            return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
        end
        # top left
        xm = -x[1]
        ym = -y[1]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, x, y, dx, dy, jmin, jmax, imin, imax) == fix_basic(
            xm, ym, x, y, dx, dy)
        # bottom left
        xm = x[1]
        ym = y[end]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, x, y, dx, dy, jmin, jmax, imin, imax) == fix_basic(
            xm, ym, x, y, dx, dy)
        # top right
        xm = x[end]
        ym = y[1]
        j=trunc(Int, (xm-x[1])/dx)+1;
        i=trunc(Int, (ym-y[1])/dy)+1;
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, x, y, dx, dy, jmin, jmax, imin, imax) == fix_basic(
            xm, ym, x, y, dx, dy)
        # bottom right
        xm = x[end]
        ym = y[end]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, x, y, dx, dy, jmin, jmax, imin, imax) == fix_basic(
            xm, ym, x, y, dx, dy)
        end # testset "basic nodes"

        @testset "Vx nodes" begin
        # verification, from madcph.m, line 434ff
        jmin, jmax = sp.jmin_vx, sp.jmax_vx
        imin, imax = sp.imin_vx, sp.imax_vx
        function fix_vx(xm, ym, x_axis, y_axis, dx, dy)
            j=trunc(Int, (xm-x_axis[1])/dx)+1;
            i=trunc(Int, (ym-y_axis[1])/dy)+1;
            if j<1
                j=1
            elseif j>Nx-1 
                j=Nx-1
            end
            if i<1 
                i=1
            elseif i>Ny
                i=Ny
            end
            dxmj=xm-xvx[j];
            dymi=ym-yvx[i];
            wtmij=(1-dxmj/dx)*(1-dymi/dy);
            wtmi1j=(1-dxmj/dx)*(dymi/dy);    
            wtmij1=(dxmj/dx)*(1-dymi/dy);
            wtmi1j1=(dxmj/dx)*(dymi/dy);
            return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
        end
        # top left
        xm = -xvx[1]
        ym = -yvx[1]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xvx, yvx, dx, dy, jmin, jmax, imin, imax) == fix_vx(
            xm, ym, xvx, yvx, dx, dy)
        # bottom left
        xm = xvx[1]
        ym = yvx[end]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xvx, yvx, dx, dy, jmin, jmax, imin, imax) == fix_vx(
            xm, ym, xvx, yvx, dx, dy)
        # top right
        xm = xvx[end]
        ym = yvx[1]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xvx, yvx, dx, dy, jmin, jmax, imin, imax) == fix_vx(
            xm, ym, xvx, yvx, dx, dy)
        # bottom right
        xm = xvx[end]
        ym = yvx[end]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xvx, yvx, dx, dy, jmin, jmax, imin, imax) == fix_vx(
            xm, ym, xvx, yvx, dx, dy)
        end # testset "Vx nodes"

        @testset "Vy nodes" begin
        # verification, from madcph.m, line 484ff
        jmin, jmax = sp.jmin_vy, sp.jmax_vy
        imin, imax = sp.imin_vy, sp.imax_vy
        function fix_vy(xm, ym, x_axis, y_axis, dx, dy)
            j=trunc(Int, (xm-x_axis[1])/dx)+1;
            i=trunc(Int, (ym-y_axis[1])/dy)+1;
            if j<1
                j=1
            elseif j>Nx 
                j=Nx
            end
            if i<1 
                i=1
            elseif i>Ny-1
                i=Ny-1
            end
            dxmj=xm-xvy[j];
            dymi=ym-yvy[i];
            wtmij=(1-dxmj/dx)*(1-dymi/dy);
            wtmi1j=(1-dxmj/dx)*(dymi/dy);    
            wtmij1=(dxmj/dx)*(1-dymi/dy);
            wtmi1j1=(dxmj/dx)*(dymi/dy);
            return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
        end
        # top left
        xm = -xvy[1]
        ym = -yvy[1]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xvy, yvy, dx, dy, jmin, jmax, imin, imax) == fix_vy(
            xm, ym, xvy, yvy, dx, dy)
        # bottom left
        xm = xvy[1]
        ym = yvy[end]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xvy, yvy, dx, dy, jmin, jmax, imin, imax) == fix_vy(
            xm, ym, xvy, yvy, dx, dy)
        # top right
        xm = xvy[end]
        ym = yvy[1]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xvy, yvy, dx, dy, jmin, jmax, imin, imax) == fix_vy(
            xm, ym, xvy, yvy, dx, dy)
        # bottom right
        xm = xvy[end]
        ym = yvy[end]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xvy, yvy, dx, dy, jmin, jmax, imin, imax) == fix_vy(
            xm, ym, xvy, yvy, dx, dy)
        end # testset "Vy nodes"
    
        @testset "P nodes" begin
        # verification, from madcph.m, line 538ff
        jmin, jmax = sp.jmin_p, sp.jmax_p
        imin, imax = sp.imin_p, sp.imax_p
        function fix_p(xm, ym, x_axis, y_axis, dx, dy)
            j=trunc(Int, (xm-x_axis[1])/dx)+1;
            i=trunc(Int, (ym-y_axis[1])/dy)+1;
            if j<1
                j=1
            elseif j>Nx 
                j=Nx
            end
            if i<1 
                i=1
            elseif i>Ny
                i=Ny
            end
            dxmj=xm-xp[j];
            dymi=ym-yp[i];
            wtmij=(1-dxmj/dx)*(1-dymi/dy);
            wtmi1j=(1-dxmj/dx)*(dymi/dy);    
            wtmij1=(dxmj/dx)*(1-dymi/dy);
            wtmi1j1=(dxmj/dx)*(dymi/dy);
            return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
        end
        # top left
        xm = -xp[1]
        ym = -yp[1]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xp, yp, dx, dy, jmin, jmax, imin, imax) == fix_p(
            xm, ym, xp, yp, dx, dy)
        # bottom left
        xm = xp[1]
        ym = yp[end]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xp, yp, dx, dy, jmin, jmax, imin, imax) == fix_p(
            xm, ym, xp, yp, dx, dy)
        # top right
        xm = xp[end]
        ym = yp[1]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xp, yp, dx, dy, jmin, jmax, imin, imax) == fix_p(
            xm, ym, xp, yp, dx, dy)
        # bottom right
        xm = xp[end]
        ym = yp[end]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xp, yp, dx, dy, jmin, jmax, imin, imax) == fix_p(
            xm, ym, xp, yp, dx, dy)
        end # testset "P nodes"    
    end # testset "fix_weights() elementary"

    @testset "fix_weights() advanced" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        dx, dy = sp.dx, sp.dy
        xsize, ysize = sp.xsize, sp.ysize
        # verification, from madcph.m, line 38ff
        # basic nodes
        x=0:dx:xsize
        y=0:dy:ysize
        # Vx nodes
        xvx=0:dx:xsize+dy
        yvx=-dy/2:dy:ysize+dy/2
        # Vy nodes
        xvy=-dx/2:dx:xsize+dx/2
        yvy=0:dy:ysize+dy
        # P Nodes
        xp=-dx/2:dx:xsize+dx/2
        yp=-dy/2:dy:ysize+dy/2
        # simulating markers
        num_markers = 10_000
        @testset "basic nodes" begin
            # verification, from madcph.m, line 373ff
            jmin, jmax = sp.jmin_basic, sp.jmax_basic
            imin, imax = sp.imin_basic, sp.imax_basic
            function fix_basic(xm, ym, x_axis, y_axis, dx, dy)
                j=trunc(Int, (xm-x_axis[1])/dx)+1;
                i=trunc(Int, (ym-y_axis[1])/dy)+1;
                if j<1
                    j=1
                elseif j>Nx-1 
                    j=Nx-1
                end
                if i<1 
                    i=1
                elseif i>Ny-1
                    i=Ny-1
                end
                dxmj=xm-x[j];
                dymi=ym-y[i];
                wtmij=(1-dxmj/dx)*(1-dymi/dy);
                wtmi1j=(1-dxmj/dx)*(dymi/dy);    
                wtmij1=(dxmj/dx)*(1-dymi/dy);
                wtmi1j1=(dxmj/dx)*(dymi/dy);
                return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
            end
            # simulating markers
            xm = rand(-x[1]:0.1:x[end]+dx, num_markers)
            ym = rand(-y[1]:0.1:y[end]+dy, num_markers)
            for m in 1:num_markers
                i, j, weights = HydrologyPlanetesimals.fix_weights(
                    xm[m],
                    ym[m],
                    x,
                    y,
                    dx,
                    dy,
                    jmin,
                    jmax,
                    imin,
                    imax
                )
                i_ver, j_ver, weights_ver = fix_basic(
                    xm[m], ym[m], x, y, dx, dy)
                @test i == i_ver
                @test j == j_ver
                @test weights == weights_ver
            end
        end # testset "basic nodes"
        
        @testset "Vx nodes" begin
            # verification, from madcph.m, line 434ff
            jmin, jmax = sp.jmin_vx, sp.jmax_vx
            imin, imax = sp.imin_vx, sp.imax_vx
            function fix_vx(xm, ym, x_axis, y_axis, dx, dy)
                j=trunc(Int, (xm-x_axis[1])/dx)+1;
                i=trunc(Int, (ym-y_axis[1])/dy)+1;
                if j<1
                    j=1
                elseif j>Nx-1 
                    j=Nx-1
                end
                if i<1 
                    i=1
                elseif i>Ny
                    i=Ny
                end
                dxmj=xm-xvx[j];
                dymi=ym-yvx[i];
                wtmij=(1-dxmj/dx)*(1-dymi/dy);
                wtmi1j=(1-dxmj/dx)*(dymi/dy);    
                wtmij1=(dxmj/dx)*(1-dymi/dy);
                wtmi1j1=(dxmj/dx)*(dymi/dy);
                return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
            end
            # simulating markers
            xm = rand(-xvx[1]:0.1:xvx[end]+dx, num_markers)
            ym = rand(-yvx[1]:0.1:yvx[end]+dy, num_markers)
            for m in 1:num_markers
                i, j, weights = HydrologyPlanetesimals.fix_weights(
                    xm[m],
                    ym[m],
                    xvx,
                    yvx,
                    dx,
                    dy,
                    jmin,
                    jmax,
                    imin,
                    imax
                )
                i_ver, j_ver, weights_ver = fix_vx(
                    xm[m], ym[m], xvx, yvx, dx, dy)
                @debug "fix_weights Vx" i i_ver j j_ver weights weights_ver
                @test i == i_ver
                @test j == j_ver
                @test weights == weights_ver
            end
        end # testset "Vx nodes"

        @testset "Vy nodes" begin
            # verification, from madcph.m, line 484ff
            jmin, jmax = sp.jmin_vy, sp.jmax_vy
            imin, imax = sp.imin_vy, sp.imax_vy
            function fix_vy(xm, ym, x_axis, y_axis, dx, dy)
                j=trunc(Int, (xm-x_axis[1])/dx)+1;
                i=trunc(Int, (ym-y_axis[1])/dy)+1;
                if j<1
                    j=1
                elseif j>Nx 
                    j=Nx
                end
                if i<1 
                    i=1
                elseif i>Ny-1
                    i=Ny-1
                end
                dxmj=xm-xvy[j];
                dymi=ym-yvy[i];
                wtmij=(1-dxmj/dx)*(1-dymi/dy);
                wtmi1j=(1-dxmj/dx)*(dymi/dy);    
                wtmij1=(dxmj/dx)*(1-dymi/dy);
                wtmi1j1=(dxmj/dx)*(dymi/dy);
                return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
            end
            # simulating markers
            xm = rand(-xvy[1]:0.1:xvy[end]+dx, num_markers)
            ym = rand(-yvy[1]:0.1:yvy[end]+dy, num_markers)
            for m in 1:num_markers
                i, j, weights = HydrologyPlanetesimals.fix_weights(
                    xm[m],
                    ym[m],
                    xvy,
                    yvy,
                    dx,
                    dy,
                    jmin,
                    jmax,
                    imin,
                    imax
                )
                i_ver, j_ver, weights_ver = fix_vy(
                    xm[m], ym[m], xvy, yvy, dx, dy)
                @test i == i_ver
                @test j == j_ver
                @test weights == weights_ver
            end
        end # testset "Vy nodes"
    
        @testset "P nodes" begin
            # verification, from madcph.m, line 538ff
            jmin, jmax = sp.jmin_p, sp.jmax_p
            imin, imax = sp.imin_p, sp.imax_p
            function fix_p(xm, ym, x_axis, y_axis, dx, dy)
                j=trunc(Int, (xm-x_axis[1])/dx)+1;
                i=trunc(Int, (ym-y_axis[1])/dy)+1;
                if j<1
                    j=1
                elseif j>Nx 
                    j=Nx
                end
                if i<1 
                    i=1
                elseif i>Ny
                    i=Ny
                end
                dxmj=xm-xp[j];
                dymi=ym-yp[i];
                wtmij=(1-dxmj/dx)*(1-dymi/dy);
                wtmi1j=(1-dxmj/dx)*(dymi/dy);    
                wtmij1=(dxmj/dx)*(1-dymi/dy);
                wtmi1j1=(dxmj/dx)*(dymi/dy);
                return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
            end
            # simulating markers
            xm = rand(-xp[1]:0.1:xp[end]+dx, num_markers)
            ym = rand(-yp[1]:0.1:yp[end]+dy, num_markers)
            for m in 1:num_markers
                i, j, weights = HydrologyPlanetesimals.fix_weights(
                    xm[m],
                    ym[m],
                    xp,
                    yp,
                    dx,
                    dy,
                    jmin,
                    jmax,
                    imin,
                    imax
                )
                i_ver, j_ver, weights_ver = fix_p(
                    xm[m], ym[m], xp, yp, dx, dy)
                @test i == i_ver
                @test j == j_ver
                @test weights == weights_ver
            end
        end # testset "P nodes"    
    end # testset "fix_weights() advanced"

    @testset "reduce_add_3darray!()" begin
        A = rand(100, 100, 10)
        A_ver = copy(A)
        # reduce-sum A along third axis
        HydrologyPlanetesimals.reduce_add_3darray!(A)
        # verification
        A_ver = reduce(+, A_ver, dims=3)
        # test
        @test A[:, :, 1] ≈ A_ver atol=1e-6
    end # testset "reduce_add_3darray!"
   
    @testset "interpolate_add_to_grid!()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        dx, dy = sp.dx, sp.dy
        xsize, ysize = sp.xsize, sp.ysize
        jmin, jmax = sp.jmin_basic, sp.jmax_basic
        imin, imax = sp.imin_basic, sp.imax_basic
        x=0:dx:xsize
        y=0:dy:ysize
        # simulate markers
        num_markers = 10_000
        xm = rand(-x[1]:0.1:x[end]+dx, num_markers)
        ym = rand(-y[1]:0.1:y[end]+dy, num_markers)
        property = rand(num_markers)
        # sample interpolation array
        for m=1:1:num_markers
            grid = zeros(Ny, Nx, Base.Threads.nthreads())
            i, j, weights = HydrologyPlanetesimals.fix_weights(
                xm[m], ym[m], x, y, dx, dy, jmin, jmax, imin, imax)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, property[m], grid)
            @test grid[i, j, Base.Threads.threadid()] ==
                property[m] * weights[1] 
            @test grid[i+1, j, Base.Threads.threadid()] ==
                property[m] * weights[2]
            @test grid[i, j+1, Base.Threads.threadid()] ==
                property[m] * weights[3]
            @test grid[i+1, j+1, Base.Threads.threadid()] ==
                property[m] * weights[4]
        end
    end # testset "interpolate_add_to_grid!()"

    @testset "interpolate_to_marker!()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        dx, dy = sp.dx, sp.dy
        xsize, ysize = sp.xsize, sp.ysize
        jmin, jmax = sp.jmin_basic, sp.jmax_basic
        imin, imax = sp.imin_basic, sp.imax_basic
        x=0:dx:xsize
        y=0:dy:ysize
        num_markers = 10_000
        xm = rand(-x[1]:0.1:x[end]+dx, num_markers)
        ym = rand(-y[1]:0.1:y[end]+dy, num_markers)
        property = zeros(num_markers)
        grid = rand(Ny, Nx)
        # interpolate to markers & test
        for m=1:1:num_markers
            i, j, weights = HydrologyPlanetesimals.fix_weights(
                xm[m], ym[m], x, y, dx, dy, jmin, jmax, imin, imax)
            HydrologyPlanetesimals.interpolate_to_marker!(
                m, i, j, weights, property, grid)
            @test property[m] ≈ (
                grid[i, j] * weights[1]
                + grid[i+1, j] * weights[2]
                + grid[i, j+1] * weights[3]
                + grid[i+1, j+1] * weights[4]
            ) atol=1e-6
        end
    end # testset "interpolate_to_marker!()"

    @testset "interpolate_add_to_marker!()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        dx, dy = sp.dx, sp.dy
        xsize, ysize = sp.xsize, sp.ysize
        jmin, jmax = sp.jmin_basic, sp.jmax_basic
        imin, imax = sp.imin_basic, sp.imax_basic
        x=0:dx:xsize
        y=0:dy:ysize
        num_markers = 10_000
        xm = rand(-x[1]:0.1:x[end]+dx, num_markers)
        ym = rand(-y[1]:0.1:y[end]+dy, num_markers)
        property = rand(num_markers)
        property_ver = copy(property)
        grid = rand(Ny, Nx)
        # interpolate to markers, verification, and test
        for m=1:1:num_markers
            i, j, weights = HydrologyPlanetesimals.fix_weights(
                xm[m], ym[m], x, y, dx, dy, jmin, jmax, imin, imax)
            HydrologyPlanetesimals.interpolate_add_to_marker!(
                m, i, j, weights, property, grid)
            @test property[m] ≈ property_ver[m] + (
                grid[i, j] * weights[1]
                + grid[i+1, j] * weights[2]
                + grid[i, j+1] * weights[3]
                + grid[i+1, j+1] * weights[4]
            ) atol=1e-6
        end
    end # testset "interpolate_add_to_marker!()"


    @testset "marker_to_basic_nodes!()" begin
        sp = HydrologyPlanetesimals.StaticParameters(Nxmc=1, Nymc=1)
        x, y = sp.x, sp.y
        dx, dy = sp.dx, sp.dy
        jmin_basic, jmax_basic = sp.jmin_basic, sp.jmax_basic
        imin_basic, imax_basic = sp.imin_basic, sp.imax_basic
        marknum = sp.start_marknum
        (
            xm,
            ym,
            _,
            _,
            _,
            sxym,
            etavpm,
            _
        ) = HydrologyPlanetesimals.setup_marker_properties(
            sp, randomized=true)
        (
            rhototalm,
            rhocptotalm,
            _,
            _,
            etatotalm,
            _,
            _,
            _,
            _,
            inv_gggtotalm,
            fricttotalm,
            cohestotalm,
            tenstotalm,
            _,
            _,
            _
        ) = HydrologyPlanetesimals.setup_marker_properties_helpers(
            sp, randomized=true)
        (
            ETA0SUM,
            ETASUM,
            GGGSUM,
            SXYSUM,
            COHSUM,
            TENSUM,
            FRISUM,
            WTSUM,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _
        ) = HydrologyPlanetesimals.setup_interpolated_properties(sp)
        ETA0SUM_ver = zero(ETA0SUM)
        ETASUM_ver = zero(ETASUM)
        GGGSUM_ver = zero(GGGSUM)
        SXYSUM_ver = zero(SXYSUM)
        COHSUM_ver = zero(COHSUM)
        TENSUM_ver = zero(TENSUM)
        FRISUM_ver = zero(FRISUM)
        WTSUM_ver = zero(WTSUM)
        # interpolate markers to basic nodes
        for m=1:1:marknum
            HydrologyPlanetesimals.marker_to_basic_nodes!(
                m,
                xm[m],
                ym[m],
                etatotalm,
                etavpm,
                inv_gggtotalm,
                sxym,
                cohestotalm,
                tenstotalm,
                fricttotalm,
                ETA0SUM,
                ETASUM,
                GGGSUM,
                SXYSUM,
                COHSUM,
                TENSUM,
                FRISUM,
                WTSUM,
                sp
            ) 
        end
        # verification
        for m=1:1:marknum
            i, j, weights = HydrologyPlanetesimals.fix_weights(
                xm[m],
                ym[m],
                x,
                y,
                dx,
                dy,
                jmin_basic,
                jmax_basic,
                imin_basic,
                imax_basic
            )
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, etatotalm[m], ETA0SUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, etavpm[m], ETASUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, inv_gggtotalm[m], GGGSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, sxym[m], SXYSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, cohestotalm[m], COHSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, tenstotalm[m], TENSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, fricttotalm[m], FRISUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, 1.0, WTSUM_ver)
        end
        # test
        @test ETA0SUM == ETA0SUM_ver
        @test ETASUM == ETASUM_ver
        @test GGGSUM == GGGSUM_ver
        @test SXYSUM == SXYSUM_ver
        @test COHSUM == COHSUM_ver
        @test TENSUM == TENSUM_ver
        @test FRISUM == FRISUM_ver
        @test WTSUM == WTSUM_ver
    end # testset "marker_to_basic_nodes!()"

    @testset "marker_to_vx_nodes!()" begin
        sp = HydrologyPlanetesimals.StaticParameters(Nxmc=1, Nymc=1)
        xvx, yvx = sp.xvx, sp.yvx
        dx, dy = sp.dx, sp.dy
        jmin_vx, jmax_vx = sp.jmin_vx, sp.jmax_vx
        imin_vx, imax_vx = sp.imin_vx, sp.imax_vx
        marknum = sp.start_marknum
        (
            xm,
            ym,
            _,
            _,
            _,
            _,
            _,
            phim
        ) = HydrologyPlanetesimals.setup_marker_properties(
            sp, randomized=true)
        (
            rhototalm,
            _,
            _,
            _,
            etatotalm,
            _,
            ktotalm,
            _,
            etafluidcur_inv_kphim,
            _,
            _,
            _,
            _,
            rhofluidcur,
            _,
            _
        ) = HydrologyPlanetesimals.setup_marker_properties_helpers(
            sp, randomized=true)
        (
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            RHOXSUM,
            RHOFXSUM,
            KXSUM,
            PHIXSUM,
            RXSUM,
            WTXSUM,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _
        ) = HydrologyPlanetesimals.setup_interpolated_properties(sp)
        RHOXSUM_ver = zero(RHOXSUM)
        RHOFXSUM_ver = zero(RHOFXSUM)
        KXSUM_ver = zero(KXSUM)
        PHIXSUM_ver = zero(PHIXSUM)
        RXSUM_ver = zero(RXSUM)
        WTXSUM_ver = zero(WTXSUM)
        # interpolate markers to Vx nodes
        for m=1:1:marknum
            HydrologyPlanetesimals.marker_to_vx_nodes!(
                m,
                xm[m],
                ym[m],
                rhototalm,
                rhofluidcur,
                ktotalm,
                phim,
                etafluidcur_inv_kphim,
                RHOXSUM,
                RHOFXSUM,
                KXSUM,
                PHIXSUM,
                RXSUM,
                WTXSUM,
                sp
            )
        end
        # verification
        for m=1:1:marknum
            i, j, weights = HydrologyPlanetesimals.fix_weights(
                xm[m],
                ym[m],
                xvx,
                yvx,
                dx,
                dy,
                jmin_vx,
                jmax_vx,
                imin_vx,
                imax_vx
            )
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, rhototalm[m], RHOXSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, rhofluidcur[m], RHOFXSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, ktotalm[m], KXSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, phim[m], PHIXSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, etafluidcur_inv_kphim[m], RXSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, 1.0, WTXSUM_ver)
        end
        # test
        @test RHOXSUM == RHOXSUM_ver
        @test RHOFXSUM == RHOFXSUM_ver
        @test KXSUM == KXSUM_ver
        @test PHIXSUM == PHIXSUM_ver
        @test RXSUM == RXSUM_ver
        @test WTXSUM == WTXSUM_ver
    end # testset "marker_to_vx_nodes!()"

    @testset "marker_to_vy_nodes!()" begin
        sp = HydrologyPlanetesimals.StaticParameters(Nxmc=1, Nymc=1)
        xvy, yvy = sp.xvy, sp.yvy
        dx, dy = sp.dx, sp.dy
        jmin_vy, jmax_vy = sp.jmin_vy, sp.jmax_vy
        imin_vy, imax_vy = sp.imin_vy, sp.imax_vy
        marknum = sp.start_marknum
        (
            xm,
            ym,
            _,
            _,
            _,
            _,
            _,
            phim
        ) = HydrologyPlanetesimals.setup_marker_properties(
            sp, randomized=true)
        (
            rhototalm,
            _,
            _,
            _,
            etatotalm,
            _,
            ktotalm,
            _,
            etafluidcur_inv_kphim,
            _,
            _,
            _,
            _,
            rhofluidcur,
            _,
            _
        ) = HydrologyPlanetesimals.setup_marker_properties_helpers(
            sp, randomized=true)
        (
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            RHOYSUM,
            RHOFYSUM,
            KYSUM,
            PHIYSUM,
            RYSUM,
            WTYSUM,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _
        ) = HydrologyPlanetesimals.setup_interpolated_properties(sp)
        RHOYSUM_ver = zero(RHOYSUM)
        RHOFYSUM_ver = zero(RHOFYSUM)
        KYSUM_ver = zero(KYSUM)
        PHIYSUM_ver = zero(PHIYSUM)
        RYSUM_ver = zero(RYSUM)
        WTYSUM_ver = zero(WTYSUM)
        # interpolate markers to Vy nodes
        for m=1:1:marknum
            HydrologyPlanetesimals.marker_to_vy_nodes!(
                m,
                xm[m],
                ym[m],
                rhototalm,
                rhofluidcur,
                ktotalm,
                phim,
                etafluidcur_inv_kphim,
                RHOYSUM,
                RHOFYSUM,
                KYSUM,
                PHIYSUM,
                RYSUM,
                WTYSUM,
                sp
            )
        end
        # verification
        for m=1:1:marknum
            i, j, weights = HydrologyPlanetesimals.fix_weights(
                xm[m],
                ym[m],
                xvy,
                yvy,
                dx,
                dy,
                jmin_vy,
                jmax_vy,
                imin_vy,
                imax_vy
            )
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, rhototalm[m], RHOYSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, rhofluidcur[m], RHOFYSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, ktotalm[m], KYSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, phim[m], PHIYSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, etafluidcur_inv_kphim[m], RYSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, 1.0, WTYSUM_ver)
        end
        # test
        @test RHOYSUM == RHOYSUM_ver
        @test RHOFYSUM == RHOFYSUM_ver
        @test KYSUM == KYSUM_ver
        @test PHIYSUM == PHIYSUM_ver
        @test RYSUM == RYSUM_ver
        @test WTYSUM == WTYSUM_ver
    end # testset "marker_to_vy_nodes!()"

    @testset "marker_to_P_nodes!()" begin
        sp = HydrologyPlanetesimals.StaticParameters(Nxmc=1, Nymc=1)
        xp, yp = sp.xp, sp.yp
        dx, dy = sp.dx, sp.dy
        jmin_p, jmax_p = sp.jmin_p, sp.jmax_p
        imin_p, imax_p = sp.imin_p, sp.imax_p
        marknum = sp.start_marknum
        (
            xm,
            ym,
            _,
            _,
            sxxm,
            _,
            _,
            phim
        ) = HydrologyPlanetesimals.setup_marker_properties(
            sp, randomized=true)
        (
            rhototalm,
            rhocptotalm,
            _,
            _,
            _,
            hrtotalm,
            _,
            tkm_rhocptotalm,
            _,
            inv_gggtotalm,
            _,
            _,
            _,
            _,
            alphasolidcur,
            alphafluidcur
        ) = HydrologyPlanetesimals.setup_marker_properties_helpers(
            sp, randomized=true)
        (
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            _,
            GGGPSUM,
            SXXSUM,
            RHOSUM,
            RHOCPSUM,
            ALPHASUM,
            ALPHAFSUM,
            HRSUM,
            PHISUM,
            TKSUM,
            WTPSUM
        ) = HydrologyPlanetesimals.setup_interpolated_properties(sp)
        GGGPSUM_ver = zero(GGGPSUM)
        SXXSUM_ver = zero(SXXSUM)
        RHOSUM_ver = zero(RHOSUM)
        RHOCPSUM_ver = zero(RHOCPSUM)
        ALPHASUM_ver = zero(ALPHASUM)
        ALPHAFSUM_ver = zero(ALPHAFSUM)
        HRSUM_ver = zero(HRSUM)
        PHISUM_ver = zero(PHISUM)
        TKSUM_ver = zero(TKSUM)
        WTPSUM_ver = zero(WTPSUM)
        # interpolate markers to P nodes
        for m=1:1:marknum
            HydrologyPlanetesimals.marker_to_p_nodes!(
                m,
                xm[m],
                ym[m],
                inv_gggtotalm,
                sxxm,
                rhototalm,
                rhocptotalm,
                alphasolidcur,
                alphafluidcur,
                hrtotalm,
                phim,
                tkm_rhocptotalm,
                GGGPSUM,
                SXXSUM,
                RHOSUM,
                RHOCPSUM,
                ALPHASUM,
                ALPHAFSUM,
                HRSUM,
                PHISUM,
                TKSUM,
                WTPSUM,
                sp
            )
        end
        # verification
        for m=1:1:marknum
            i, j, weights = HydrologyPlanetesimals.fix_weights(
                xm[m],
                ym[m],
                xp,
                yp,
                dx,
                dy,
                jmin_p,
                jmax_p,
                imin_p,
                imax_p
            )
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, inv_gggtotalm[m], GGGPSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, sxxm[m], SXXSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, rhototalm[m], RHOSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, rhocptotalm[m], RHOCPSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, alphafluidcur[m], ALPHAFSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, hrtotalm[m], HRSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, phim[m], PHISUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, tkm_rhocptotalm[m], TKSUM_ver)
            HydrologyPlanetesimals.interpolate_add_to_grid!(
                i, j, weights, 1.0, WTPSUM_ver)
        end
        # test

        @test WTPSUM == WTPSUM_ver
    end # testset "marker_to_p_nodes!()"



    @testset "compute node properties: basic, Vx, Vy, P" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        dx, dy = sp.dx, sp.dy
        xsize, ysize = sp.xsize, sp.ysize
        # verification, from madcph.m, line 38ff
        # basic nodes
        x=0:dx:xsize
        y=0:dy:ysize
        # Vx nodes
        xvx=0:dx:xsize+dy
        yvx=-dy/2:dy:ysize+dy/2
        # Vy nodes
        xvy=-dx/2:dx:xsize+dx/2
        yvy=0:dy:ysize+dy
        # P Nodes
        xp=-dx/2:dx:xsize+dx/2
        yp=-dy/2:dy:ysize+dy/2
        # simulating markers
        num_markers = 10_000
        @testset "compute_basic_node_properties!()" begin    
            jmin, jmax = sp.jmin_basic, sp.jmax_basic
            imin, imax = sp.imin_basic, sp.imax_basic
            ETA0SUM = zeros(Ny, Nx, Base.Threads.nthreads())
            ETASUM = zeros(Ny, Nx, Base.Threads.nthreads())
            GGGSUM = zeros(Ny, Nx, Base.Threads.nthreads())
            SXYSUM = zeros(Ny, Nx, Base.Threads.nthreads())
            COHSUM = zeros(Ny, Nx, Base.Threads.nthreads())
            TENSUM = zeros(Ny, Nx, Base.Threads.nthreads())
            FRISUM = zeros(Ny, Nx, Base.Threads.nthreads())
            WTSUM = zeros(Ny, Nx, Base.Threads.nthreads())
            ETA0 = zeros(Float64, Ny, Nx)
            ETA = zeros(Float64, Ny, Nx)
            GGG = zeros(Float64, Ny, Nx)
            SXY0 = zeros(Float64, Ny, Nx)
            COH = zeros(Float64, Ny, Nx)
            TEN = zeros(Float64, Ny, Nx)
            FRI = zeros(Float64, Ny, Nx)
            YNY = zeros(Bool, Ny, Nx)
            ETA0SUM_ver = zeros(Ny, Nx)
            ETASUM_ver = zeros(Ny, Nx)
            GGGSUM_ver = zeros(Ny, Nx)
            SXYSUM_ver = zeros(Ny, Nx)
            COHSUM_ver = zeros(Ny, Nx)
            TENSUM_ver = zeros(Ny, Nx)
            FRISUM_ver = zeros(Ny, Nx)
            WTSUM_ver = zeros(Ny, Nx)
            ETA0_ver = zeros(Float64, Ny, Nx)
            ETA_ver = zeros(Float64, Ny, Nx)
            GGG_ver = zeros(Float64, Ny, Nx)
            SXY0_ver = zeros(Float64, Ny, Nx)
            COH_ver = zeros(Float64, Ny, Nx)
            TEN_ver = zeros(Float64, Ny, Nx)
            FRI_ver = zeros(Float64, Ny, Nx)
            YNY_ver = zeros(Bool, Ny, Nx)
            # simulate markers
            xm = rand(-x[1]:0.1:x[end]+dx, num_markers)
            ym = rand(-y[1]:0.1:y[end]+dy, num_markers)
            property = rand(7, num_markers)
            # calculate grid properties
            Base.Threads.@threads for m=1:1:num_markers
                i, j, weights = HydrologyPlanetesimals.fix_weights(
                    xm[m], ym[m], x, y, dx, dy, jmin, jmax, imin, imax)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[1, m], ETA0SUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[2, m], ETASUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, inv(property[3, m]), GGGSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[4, m], SXYSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[5, m], COHSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[6, m], TENSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[7, m], FRISUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, 1.0, WTSUM)
            end
            HydrologyPlanetesimals.compute_basic_node_properties!(
                ETA0SUM,
                ETASUM,
                GGGSUM,
                SXYSUM,
                COHSUM,
                TENSUM,
                FRISUM,
                WTSUM,
                ETA0,
                ETA,
                GGG,
                SXY0,
                COH,
                TEN,
                FRI,
                YNY
            )
            # verification properties, from madcph.m, lines 373ff, 606ff
            for m=1:1:num_markers
                j=trunc(Int, (xm[m]-x[1])/dx)+1
                i=trunc(Int, (ym[m]-y[1])/dy)+1
                if j<1
                    j=1
                elseif j>Nx-1
                    j=Nx-1
                end
                if i<1
                    i=1
                elseif i>Ny-1 
                    i=Ny-1
                end
                # Compute distances
                dxmj=xm[m]-x[j]
                dymi=ym[m]-y[i]
                # Compute weights
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy);    
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Update properties
                # i;j Node
                ETA0SUM_ver[i,j]=ETA0SUM_ver[i,j]+property[1,m]*wtmij
                ETASUM_ver[i,j]=ETASUM_ver[i,j]+property[2,m]*wtmij
                GGGSUM_ver[i,j]=GGGSUM_ver[i,j]+1/property[3,m]*wtmij
                SXYSUM_ver[i,j]=SXYSUM_ver[i,j]+property[4,m]*wtmij
                COHSUM_ver[i,j]=COHSUM_ver[i,j]+property[5,m]*wtmij
                TENSUM_ver[i,j]=TENSUM_ver[i,j]+property[6,m]*wtmij
                FRISUM_ver[i,j]=FRISUM_ver[i,j]+property[7,m]*wtmij
                WTSUM_ver[i,j]=WTSUM_ver[i,j]+wtmij
                # i+1;j Node
                ETA0SUM_ver[i+1,j]=ETA0SUM_ver[i+1,j]+property[1,m]*wtmi1j
                ETASUM_ver[i+1,j]=ETASUM_ver[i+1,j]+property[2,m]*wtmi1j
                GGGSUM_ver[i+1,j]=GGGSUM_ver[i+1,j]+1/property[3,m]*wtmi1j
                SXYSUM_ver[i+1,j]=SXYSUM_ver[i+1,j]+property[4,m]*wtmi1j
                COHSUM_ver[i+1,j]=COHSUM_ver[i+1,j]+property[5,m]*wtmi1j
                TENSUM_ver[i+1,j]=TENSUM_ver[i+1,j]+property[6,m]*wtmi1j
                FRISUM_ver[i+1,j]=FRISUM_ver[i+1,j]+property[7,m]*wtmi1j
                WTSUM_ver[i+1,j]=WTSUM_ver[i+1,j]+wtmi1j
                # i;j+1 Node
                ETA0SUM_ver[i,j+1]=ETA0SUM_ver[i,j+1]+property[1,m]*wtmij1
                ETASUM_ver[i,j+1]=ETASUM_ver[i,j+1]+property[2,m]*wtmij1
                GGGSUM_ver[i,j+1]=GGGSUM_ver[i,j+1]+1/property[3,m]*wtmij1
                SXYSUM_ver[i,j+1]=SXYSUM_ver[i,j+1]+property[4,m]*wtmij1
                COHSUM_ver[i,j+1]=COHSUM_ver[i,j+1]+property[5,m]*wtmij1
                TENSUM_ver[i,j+1]=TENSUM_ver[i,j+1]+property[6,m]*wtmij1
                FRISUM_ver[i,j+1]=FRISUM_ver[i,j+1]+property[7,m]*wtmij1
                WTSUM_ver[i,j+1]=WTSUM_ver[i,j+1]+wtmij1
                # i+1;j+1 Node
                ETA0SUM_ver[i+1,j+1]=ETA0SUM_ver[i+1,j+1]+property[1,m]*wtmi1j1
                ETASUM_ver[i+1,j+1]=ETASUM_ver[i+1,j+1]+property[2,m]*wtmi1j1
                GGGSUM_ver[i+1,j+1]=GGGSUM_ver[i+1,j+1]+1/property[3,m]*wtmi1j1
                SXYSUM_ver[i+1,j+1]=SXYSUM_ver[i+1,j+1]+property[4,m]*wtmi1j1
                COHSUM_ver[i+1,j+1]=COHSUM_ver[i+1,j+1]+property[5,m]*wtmi1j1
                TENSUM_ver[i+1,j+1]=TENSUM_ver[i+1,j+1]+property[6,m]*wtmi1j1
                FRISUM_ver[i+1,j+1]=FRISUM_ver[i+1,j+1]+property[7,m]*wtmi1j1
                WTSUM_ver[i+1,j+1]=WTSUM_ver[i+1,j+1]+wtmi1j1
            end
            for j=1:1:Nx
                for i=1:1:Ny
                    if WTSUM_ver[i,j]>0 
                        ETA0_ver[i,j]=ETA0SUM_ver[i,j]/WTSUM_ver[i,j]
                        ETA_ver[i,j]=ETASUM_ver[i,j]/WTSUM_ver[i,j]
                        if(ETA_ver[i,j]<ETA0_ver[i,j])
                            YNY_ver[i,j]=1
                        end
                        GGG_ver[i,j]=1/(GGGSUM_ver[i,j]/WTSUM_ver[i,j])
                        SXY0_ver[i,j]=SXYSUM_ver[i,j]/WTSUM_ver[i,j]
                        COH_ver[i,j]=COHSUM_ver[i,j]/WTSUM_ver[i,j]
                        TEN_ver[i,j]=TENSUM_ver[i,j]/WTSUM_ver[i,j]
                        FRI_ver[i,j]=FRISUM_ver[i,j]/WTSUM_ver[i,j]
                    end
                end
            end 
            # test
            for j=1:1:Nx, i=1:1:Ny
                @test ETA0[i, j] ≈ ETA0_ver[i, j] atol=1e-6
                @test ETA[i, j] ≈ ETA_ver[i, j] atol=1e-6
                @test GGG[i, j] ≈ GGG_ver[i, j] atol=1e-6
                @test SXY0[i, j] ≈ SXY0_ver[i, j] atol=1e-6
                @test COH[i, j] ≈ COH_ver[i, j] atol=1e-6
                @test TEN[i, j] ≈ TEN_ver[i, j] atol=1e-6
                @test FRI[i, j] ≈ FRI_ver[i, j] atol=1e-6
                @test YNY[i, j] ≈ YNY_ver[i, j] atol=1e-6
            end
        end # testset "compute_basic_node_properties!()"

        @testset "compute_vx_node_properties!()" begin
            jmin, jmax = sp.jmin_vx, sp.jmax_vx
            imin, imax = sp.imin_vx, sp.imax_vx 
            RHOXSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            RHOFXSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            KXSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            PHIXSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            RXSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            WTXSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            RHOX = zeros(Float64, Ny1, Nx1)
            RHOFX = zeros(Float64, Ny1, Nx1)
            KX = zeros(Float64, Ny1, Nx1)
            PHIX = zeros(Float64, Ny1, Nx1)
            RX = zeros(Float64, Ny1, Nx1)
            RHOXSUM_ver = zeros(Ny1, Nx1)
            RHOFXSUM_ver = zeros(Ny1, Nx1)
            KXSUM_ver = zeros(Ny1, Nx1)
            PHIXSUM_ver = zeros(Ny1, Nx1)
            RXSUM_ver = zeros(Ny1, Nx1)
            WTXSUM_ver = zeros(Ny1, Nx1)
            RHOX_ver = zeros(Float64, Ny1, Nx1)
            RHOFX_ver = zeros(Float64, Ny1, Nx1)
            KX_ver = zeros(Float64, Ny1, Nx1)
            PHIX_ver = zeros(Float64, Ny1, Nx1)
            RX_ver = zeros(Float64, Ny1, Nx1)
            # simulate markers
            xm = rand(-xvx[1]:0.1:xvx[end]+dx, num_markers)
            ym = rand(-yvx[1]:0.1:yvx[end]+dy, num_markers)
            property = rand(5, num_markers)
            # calculate grid properties
            for m=1:1:num_markers
                i, j, weights = HydrologyPlanetesimals.fix_weights(
                    xm[m], ym[m], xvx, yvx, dx, dy, jmin, jmax, imin, imax)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[1, m], RHOXSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[2, m], RHOFXSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[3, m], KXSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[4, m], PHIXSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[5, m], RXSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, 1.0, WTXSUM)
            end
            HydrologyPlanetesimals.compute_vx_node_properties!(
                RHOXSUM,
                RHOFXSUM,
                KXSUM,
                PHIXSUM,
                RXSUM,
                WTXSUM,
                RHOX,
                RHOFX,
                KX,
                PHIX,
                RX
            )
            # verification properties, from madcph.m, lines 434ff, 624ff
            for m=1:1:num_markers
                j=trunc(Int, (xm[m]-xvx[1])/dx)+1
                i=trunc(Int, (ym[m]-yvx[1])/dy)+1
                if j<1
                    j=1
                elseif j>Nx-1
                    j=Nx-1
                end
                if(i<1)
                    i=1
                elseif i>Ny
                    i=Ny
                end
                # Compute distances
                dxmj=xm[m]-xvx[j]
                dymi=ym[m]-yvx[i]
                # Compute weights
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy);    
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Update properties
                # i;j Node
                RHOXSUM_ver[i,j]=RHOXSUM_ver[i,j]+property[1, m]*wtmij
                RHOFXSUM_ver[i,j]=RHOFXSUM_ver[i,j]+property[2, m]*wtmij
                KXSUM_ver[i,j]=KXSUM_ver[i,j]+property[3, m]*wtmij
                PHIXSUM_ver[i,j]=PHIXSUM_ver[i,j]+property[4, m]*wtmij
                RXSUM_ver[i,j]=RXSUM_ver[i,j]+property[5, m]*wtmij
                WTXSUM_ver[i,j]=WTXSUM_ver[i,j]+wtmij
                # i+1;j Node
                RHOXSUM_ver[i+1,j]=RHOXSUM_ver[i+1,j]+property[1, m]*wtmi1j
                RHOFXSUM_ver[i+1,j]=RHOFXSUM_ver[i+1,j]+property[2, m]*wtmi1j
                KXSUM_ver[i+1,j]=KXSUM_ver[i+1,j]+property[3, m]*wtmi1j
                PHIXSUM_ver[i+1,j]=PHIXSUM_ver[i+1,j]+property[4, m]*wtmi1j
                RXSUM_ver[i+1,j]=RXSUM_ver[i+1,j]+property[5, m]*wtmi1j
                WTXSUM_ver[i+1,j]=WTXSUM_ver[i+1,j]+wtmi1j
                # i;j+1 Node
                RHOXSUM_ver[i,j+1]=RHOXSUM_ver[i,j+1]+property[1, m]*wtmij1
                RHOFXSUM_ver[i,j+1]=RHOFXSUM_ver[i,j+1]+property[2, m]*wtmij1
                KXSUM_ver[i,j+1]=KXSUM_ver[i,j+1]+property[3, m]*wtmij1
                PHIXSUM_ver[i,j+1]=PHIXSUM_ver[i,j+1]+property[4, m]*wtmij1
                RXSUM_ver[i,j+1]=RXSUM_ver[i,j+1]+property[5, m]*wtmij1
                WTXSUM_ver[i,j+1]=WTXSUM_ver[i,j+1]+wtmij1
                # i+1;j+1 Node
                RHOXSUM_ver[i+1,j+1]=RHOXSUM_ver[i+1,j+1]+property[1, m]*wtmi1j1
                RHOFXSUM_ver[i+1,j+1]=RHOFXSUM_ver[i+1,j+1]+property[2, m]*wtmi1j1
                KXSUM_ver[i+1,j+1]=KXSUM_ver[i+1,j+1]+property[3, m]*wtmi1j1
                PHIXSUM_ver[i+1,j+1]=PHIXSUM_ver[i+1,j+1]+property[4, m]*wtmi1j1
                RXSUM_ver[i+1,j+1]=RXSUM_ver[i+1,j+1]+property[5, m]*wtmi1j1
                WTXSUM_ver[i+1,j+1]=WTXSUM_ver[i+1,j+1]+wtmi1j1
            end
            for j=1:1:Nx1
                for i=1:1:Ny1
                    if(WTXSUM_ver[i,j]>0)
                        RHOX_ver[i,j]=RHOXSUM_ver[i,j]/WTXSUM_ver[i,j]
                        RHOFX_ver[i,j]=RHOFXSUM_ver[i,j]/WTXSUM_ver[i,j]
                        KX_ver[i,j]=KXSUM_ver[i,j]/WTXSUM_ver[i,j]
                        PHIX_ver[i,j]=PHIXSUM_ver[i,j]/WTXSUM_ver[i,j]
                        RX_ver[i,j]=RXSUM_ver[i,j]/WTXSUM_ver[i,j]
                    end
                end
            end
            # test
            for j=1:1:Nx1, i=1:1:Ny1
                @test RHOX[i, j] ≈ RHOX_ver[i, j] atol=1e-6
                @test RHOFX[i, j] ≈ RHOFX_ver[i, j] atol=1e-6
                @test KX[i, j] ≈ KX_ver[i, j] atol=1e-6
                @test PHIX[i, j] ≈ PHIX_ver[i, j] atol=1e-6
                @test RX[i, j] ≈ RX_ver[i, j] atol=1e-6
            end
        end # testset "compute_vx_node_properties!()"

        @testset "compute_vy_node_properties!()" begin
            jmin, jmax = sp.jmin_vy, sp.jmax_vy
            imin, imax = sp.imin_vy, sp.imax_vy
            RHOYSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            RHOFYSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            KYSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            PHIYSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            RYSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            WTYSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            RHOY = zeros(Float64, Ny1, Nx1)
            RHOFY = zeros(Float64, Ny1, Nx1)
            KY = zeros(Float64, Ny1, Nx1)
            PHIY = zeros(Float64, Ny1, Nx1)
            RY = zeros(Float64, Ny1, Nx1)
            RHOYSUM_ver = zeros(Ny1, Nx1)
            RHOFYSUM_ver = zeros(Ny1, Nx1)
            KYSUM_ver = zeros(Ny1, Nx1)
            PHIYSUM_ver = zeros(Ny1, Nx1)
            RYSUM_ver = zeros(Ny1, Nx1)
            WTYSUM_ver = zeros(Ny1, Nx1)
            RHOY_ver = zeros(Float64, Ny1, Nx1)
            RHOFY_ver = zeros(Float64, Ny1, Nx1)
            KY_ver = zeros(Float64, Ny1, Nx1)
            PHIY_ver = zeros(Float64, Ny1, Nx1)
            RY_ver = zeros(Float64, Ny1, Nx1)
            # simulate markers
            xm = rand(-xvy[1]:0.1:xvy[end]+dx, num_markers)
            ym = rand(-yvy[1]:0.1:yvy[end]+dy, num_markers)
            property = rand(5, num_markers)
            # calculate grid properties
            for m=1:1:num_markers
                i, j, weights = HydrologyPlanetesimals.fix_weights(
                    xm[m], ym[m], xvy, yvy, dx, dy, jmin, jmax, imin, imax)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[1, m], RHOYSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[2, m], RHOFYSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[3, m], KYSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[4, m], PHIYSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[5, m], RYSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, 1.0, WTYSUM)
            end
            HydrologyPlanetesimals.compute_vy_node_properties!(
                RHOYSUM,
                RHOFYSUM,
                KYSUM,
                PHIYSUM,
                RYSUM,
                WTYSUM,
                RHOY,
                RHOFY,
                KY,
                PHIY,
                RY
            )
            # verification properties, from madcph.m, lines 486ff, 636ff
            for m=1:1:num_markers
                j=trunc(Int, (xm[m]-xvy[1])/dx)+1
                i=trunc(Int, (ym[m]-yvy[1])/dy)+1
                if j<1 
                    j=1
                elseif j>Nx 
                    j=Nx
                end
                if i<1 
                    i=1
                elseif i>Ny-1 
                    i=Ny-1
                end
                # Compute distances
                dxmj=xm[m]-xvy[j]
                dymi=ym[m]-yvy[i]
                # Compute weights
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy);    
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Update properties
                # i;j Node
                RHOYSUM_ver[i,j]=RHOYSUM_ver[i,j]+property[1, m]*wtmij
                RHOFYSUM_ver[i,j]=RHOFYSUM_ver[i,j]+property[2, m]*wtmij
                KYSUM_ver[i,j]=KYSUM_ver[i,j]+property[3, m]*wtmij
                PHIYSUM_ver[i,j]=PHIYSUM_ver[i,j]+property[4, m]*wtmij
                RYSUM_ver[i,j]=RYSUM_ver[i,j]+property[5, m]*wtmij
                WTYSUM_ver[i,j]=WTYSUM_ver[i,j]+wtmij
                # i+1;j Node
                RHOYSUM_ver[i+1,j]=RHOYSUM_ver[i+1,j]+property[1, m]*wtmi1j
                RHOFYSUM_ver[i+1,j]=RHOFYSUM_ver[i+1,j]+property[2, m]*wtmi1j
                KYSUM_ver[i+1,j]=KYSUM_ver[i+1,j]+property[3, m]*wtmi1j
                PHIYSUM_ver[i+1,j]=PHIYSUM_ver[i+1,j]+property[4, m]*wtmi1j
                RYSUM_ver[i+1,j]=RYSUM_ver[i+1,j]+property[5, m]*wtmi1j
                WTYSUM_ver[i+1,j]=WTYSUM_ver[i+1,j]+wtmi1j
                # i;j+1 Node
                RHOYSUM_ver[i,j+1]=RHOYSUM_ver[i,j+1]+property[1, m]*wtmij1
                RHOFYSUM_ver[i,j+1]=RHOFYSUM_ver[i,j+1]+property[2, m]*wtmij1
                KYSUM_ver[i,j+1]=KYSUM_ver[i,j+1]+property[3, m]*wtmij1
                PHIYSUM_ver[i,j+1]=PHIYSUM_ver[i,j+1]+property[4, m]*wtmij1
                RYSUM_ver[i,j+1]=RYSUM_ver[i,j+1]+property[5, m]*wtmij1
                WTYSUM_ver[i,j+1]=WTYSUM_ver[i,j+1]+wtmij1
                # i+1;j+1 Node
                RHOYSUM_ver[i+1,j+1]=RHOYSUM_ver[i+1,j+1]+property[1, m]*wtmi1j1
                RHOFYSUM_ver[i+1,j+1]=RHOFYSUM_ver[i+1,j+1]+property[2, m]*wtmi1j1
                KYSUM_ver[i+1,j+1]=KYSUM_ver[i+1,j+1]+property[3, m]*wtmi1j1
                PHIYSUM_ver[i+1,j+1]=PHIYSUM_ver[i+1,j+1]+property[4, m]*wtmi1j1
                RYSUM_ver[i+1,j+1]=RYSUM_ver[i+1,j+1]+property[5, m]*wtmi1j1
                WTYSUM_ver[i+1,j+1]=WTYSUM_ver[i+1,j+1]+wtmi1j1
            end
            for j=1:1:Nx1
                for i=1:1:Ny1
                    if WTYSUM_ver[i,j]>0 
                        RHOY_ver[i,j]=RHOYSUM_ver[i,j]/WTYSUM_ver[i,j]
                        RHOFY_ver[i,j]=RHOFYSUM_ver[i,j]/WTYSUM_ver[i,j]
                        KY_ver[i,j]=KYSUM_ver[i,j]/WTYSUM_ver[i,j]
                        PHIY_ver[i,j]=PHIYSUM_ver[i,j]/WTYSUM_ver[i,j]
                        RY_ver[i,j]=RYSUM_ver[i,j]/WTYSUM_ver[i,j]
                    end
                end
            end
            #test
            for j=1:1:Nx1, i=1:1:Ny1
                @test RHOY[i,j] ≈ RHOY_ver[i,j] atol=1e-6
                @test RHOFY[i,j] ≈ RHOFY_ver[i,j] atol=1e-6
                @test KY[i,j] ≈ KY_ver[i,j] atol=1e-6
                @test PHIY[i,j] ≈ PHIY_ver[i,j] atol=1e-6
                @test RY[i,j] ≈ RY_ver[i,j] atol=1e-6
            end
        end # testset "compute_vy_node_properties!()"

        @testset "compute_p_node_properties!()" begin
            jmin, jmax = sp.jmin_p, sp.jmax_p
            imin, imax = sp.imin_p, sp.imax_p
            RHOSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            RHOCPSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            ALPHASUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            ALPHAFSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            HRSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            GGGPSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            SXXSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            TKSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            PHISUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            WTPSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            RHO = zeros(Float64, Ny1, Nx1)
            RHOCP = zeros(Float64, Ny1, Nx1)
            ALPHA = zeros(Float64, Ny1, Nx1)
            ALPHAF = zeros(Float64, Ny1, Nx1)
            HR = zeros(Float64, Ny1, Nx1)
            GGGP = zeros(Float64, Ny1, Nx1)
            SXX0 = zeros(Float64, Ny1, Nx1)
            tk1 = zeros(Float64, Ny1, Nx1)
            PHI = zeros(Float64, Ny1, Nx1)
            BETTAPHI = zeros(Float64, Ny1, Nx1)
            RHOSUM_ver = zeros(Ny1, Nx1)
            RHOCPSUM_ver = zeros(Ny1, Nx1)
            ALPHASUM_ver = zeros(Ny1, Nx1)
            ALPHAFSUM_ver = zeros(Ny1, Nx1)
            HRSUM_ver = zeros(Ny1, Nx1)
            GGGPSUM_ver = zeros(Ny1, Nx1)
            SXXSUM_ver = zeros(Ny1, Nx1)
            TKSUM_ver = zeros(Ny1, Nx1)
            PHISUM_ver = zeros(Ny1, Nx1)
            WTPSUM_ver = zeros(Ny1, Nx1)
            RHO_ver = zeros(Float64, Ny1, Nx1)
            RHOCP_ver = zeros(Float64, Ny1, Nx1)
            ALPHA_ver = zeros(Float64, Ny1, Nx1)
            ALPHAF_ver = zeros(Float64, Ny1, Nx1)
            HR_ver = zeros(Float64, Ny1, Nx1)
            GGGP_ver = zeros(Float64, Ny1, Nx1)
            SXX0_ver = zeros(Float64, Ny1, Nx1)
            tk1_ver = zeros(Float64, Ny1, Nx1)
            PHI_ver = zeros(Float64, Ny1, Nx1)
            BETTAPHI_ver = zeros(Float64, Ny1, Nx1)
            # simulate markers
            xm = rand(-xp[1]:0.1:xp[end]+dx, num_markers)
            ym = rand(-yp[1]:0.1:yp[end]+dy, num_markers)
            property = rand(9, num_markers)
            # calculate grid properties
            for m=1:1:num_markers
                i, j, weights = HydrologyPlanetesimals.fix_weights(
                    xm[m], ym[m], xp, yp, dx, dy, jmin, jmax, imin, imax)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[1, m], RHOSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[2, m], RHOCPSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[3, m], ALPHASUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[4, m], ALPHAFSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[5, m], HRSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, inv(property[6, m]), GGGPSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[7, m], SXXSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[2, m] * property[8, m], TKSUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, property[9, m], PHISUM)
                HydrologyPlanetesimals.interpolate_add_to_grid!(
                    i, j, weights, 1.0, WTPSUM)
            end
            HydrologyPlanetesimals.compute_p_node_properties!(
                RHOSUM,
                RHOCPSUM,
                ALPHASUM,
                ALPHAFSUM,
                HRSUM,
                GGGPSUM,
                SXXSUM,
                TKSUM,
                PHISUM,
                WTPSUM,
                RHO,
                RHOCP,
                ALPHA,
                ALPHAF,
                HR,
                GGGP,
                SXX0,
                tk1,
                PHI,
                BETTAPHI
            )
            # verification properties, from madcph.m, lines 538ff, 648ff
            for m=1:1:num_markers
                j=trunc(Int, (xm[m]-xp[1])/dx)+1
                i=trunc(Int, (ym[m]-yp[1])/dy)+1
                if j<1 
                    j=1
                elseif j>Nx
                    j=Nx
                end
                if i<1 
                    i=1
                elseif i>Ny 
                    i=Ny
                end
                # Compute distances
                dxmj=xm[m]-xp[j]
                dymi=ym[m]-yp[i]
                # Compute weights
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy);    
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Update properties
                # i;j Node
                GGGPSUM_ver[i,j]=GGGPSUM_ver[i,j]+1/property[6, m]*wtmij
                SXXSUM_ver[i,j]=SXXSUM_ver[i,j]+property[7, m]*wtmij
                RHOSUM_ver[i,j]=RHOSUM_ver[i,j]+property[1, m]*wtmij
                RHOCPSUM_ver[i,j]=RHOCPSUM_ver[i,j]+property[2, m]*wtmij
                ALPHASUM_ver[i,j]=ALPHASUM_ver[i,j]+property[3, m]*wtmij
                ALPHAFSUM_ver[i,j]=ALPHAFSUM_ver[i,j]+property[4, m]*wtmij
                HRSUM_ver[i,j]=HRSUM_ver[i,j]+property[5, m]*wtmij
                TKSUM_ver[i,j]=TKSUM_ver[i,j]+property[8, m]*
                    property[2, m]*wtmij
                PHISUM_ver[i,j]=PHISUM_ver[i,j]+property[9, m]*wtmij
                WTPSUM_ver[i,j]=WTPSUM_ver[i,j]+wtmij
                # i+1;j Node
                GGGPSUM_ver[i+1,j]=GGGPSUM_ver[i+1,j]+1/property[6, m]*wtmi1j
                SXXSUM_ver[i+1,j]=SXXSUM_ver[i+1,j]+property[7, m]*wtmi1j
                RHOSUM_ver[i+1,j]=RHOSUM_ver[i+1,j]+property[1, m]*wtmi1j
                RHOCPSUM_ver[i+1,j]=RHOCPSUM_ver[i+1,j]+property[2, m]*wtmi1j
                ALPHASUM_ver[i+1,j]=ALPHASUM_ver[i+1,j]+property[3, m]*wtmi1j
                ALPHAFSUM_ver[i+1,j]=ALPHAFSUM_ver[i+1,j]+property[4, m]*wtmi1j
                HRSUM_ver[i+1,j]=HRSUM_ver[i+1,j]+property[5, m]*wtmi1j
                TKSUM_ver[i+1,j]=TKSUM_ver[i+1,j]+property[8, m]*
                    property[2, m]*wtmi1j
                PHISUM_ver[i+1,j]=PHISUM_ver[i+1,j]+property[9, m]*wtmi1j
                WTPSUM_ver[i+1,j]=WTPSUM_ver[i+1,j]+wtmi1j
                # i;j+1 Node
                GGGPSUM_ver[i,j+1]=GGGPSUM_ver[i,j+1]+1/property[6, m]*wtmij1
                SXXSUM_ver[i,j+1]=SXXSUM_ver[i,j+1]+property[7, m]*wtmij1
                RHOSUM_ver[i,j+1]=RHOSUM_ver[i,j+1]+property[1, m]*wtmij1
                RHOCPSUM_ver[i,j+1]=RHOCPSUM_ver[i,j+1]+property[2, m]*wtmij1
                ALPHASUM_ver[i,j+1]=ALPHASUM_ver[i,j+1]+property[3, m]*wtmij1
                ALPHAFSUM_ver[i,j+1]=ALPHAFSUM_ver[i,j+1]+property[4, m]*wtmij1
                HRSUM_ver[i,j+1]=HRSUM_ver[i,j+1]+property[5, m]*wtmij1
                TKSUM_ver[i,j+1]=TKSUM_ver[i,j+1]+property[8, m]*
                    property[2, m]*wtmij1
                PHISUM_ver[i,j+1]=PHISUM_ver[i,j+1]+property[9, m]*wtmij1
                WTPSUM_ver[i,j+1]=WTPSUM_ver[i,j+1]+wtmij1
                # i+1;j+1 Node
                GGGPSUM_ver[i+1,j+1]=GGGPSUM_ver[i+1,j+1]+1/property[6, m]*
                    wtmi1j1
                SXXSUM_ver[i+1,j+1]=SXXSUM_ver[i+1,j+1]+property[7, m]*wtmi1j1
                RHOSUM_ver[i+1,j+1]=RHOSUM_ver[i+1,j+1]+property[1, m]*wtmi1j1
                RHOCPSUM_ver[i+1,j+1]=RHOCPSUM_ver[i+1,j+1]+
                    property[2, m]*wtmi1j1
                ALPHASUM_ver[i+1,j+1]=ALPHASUM_ver[i+1,j+1]+
                    property[3, m]*wtmi1j1
                ALPHAFSUM_ver[i+1,j+1]=ALPHAFSUM_ver[i+1,j+1]+
                    property[4, m]*wtmi1j1
                HRSUM_ver[i+1,j+1]=HRSUM_ver[i+1,j+1]+property[5, m]*wtmi1j1
                TKSUM_ver[i+1,j+1]=TKSUM_ver[i+1,j+1]+property[8, m]*
                    property[2, m]*wtmi1j1
                PHISUM_ver[i+1,j+1]=PHISUM_ver[i+1,j+1]+property[9, m]*wtmi1j1
                WTPSUM_ver[i+1,j+1]=WTPSUM_ver[i+1,j+1]+wtmi1j1
            end
            for j=1:1:Nx1
                for i=1:1:Ny1
                    if WTPSUM_ver[i,j]>0
                        GGGP_ver[i,j]=1/(GGGPSUM_ver[i,j]/WTPSUM_ver[i,j])
                        SXX0_ver[i,j]=SXXSUM_ver[i,j]/WTPSUM_ver[i,j]
                        RHO_ver[i,j]=RHOSUM_ver[i,j]/WTPSUM_ver[i,j]
                        RHOCP_ver[i,j]=RHOCPSUM_ver[i,j]/WTPSUM_ver[i,j]
                        ALPHA_ver[i,j]=ALPHASUM_ver[i,j]/WTPSUM_ver[i,j]
                        ALPHAF_ver[i,j]=ALPHAFSUM_ver[i,j]/WTPSUM_ver[i,j]
                        HR_ver[i,j]=HRSUM_ver[i,j]/WTPSUM_ver[i,j]
                        PHI_ver[i,j]=PHISUM_ver[i,j]/WTPSUM_ver[i,j]
                        BETTAPHI_ver[i,j]=1/GGGP_ver[i,j]*PHI_ver[i,j]
                        tk1_ver[i,j]=TKSUM_ver[i,j]/RHOCPSUM_ver[i,j]
                    end
                end
            end
            # test
            for j=1:1:Nx, i=1:1:Ny
                @test RHO[i, j] ≈ RHO_ver[i, j] atol=1e-6
                @test RHOCP[i, j] ≈ RHOCP_ver[i, j] atol=1e-6
                @test ALPHA[i, j] ≈ ALPHA_ver[i, j] atol=1e-6
                @test ALPHAF[i, j] ≈ ALPHAF_ver[i, j] atol=1e-6
                @test HR[i, j] ≈ HR_ver[i, j] atol=1e-6
                @test GGGP[i, j] ≈ GGGP_ver[i, j] atol=1e-6
                @test SXX0[i, j] ≈ SXX0_ver[i, j] atol=1e-6
                @test tk1[i, j] ≈ tk1_ver[i, j] atol=1e-6
                @test PHI[i, j] ≈ PHI_ver[i, j] atol=1e-6
                @test BETTAPHI[i, j] ≈ BETTAPHI_ver[i, j] atol=1e-6
            end
        end # testset "compute_p_node_properties!()"
    end # testset "compute node properties" 

    @testset "apply_insulating_boundary_conditions!()" begin
        max_size = 10
        for j=3:1:max_size, i=3:1:max_size
            t = rand(i, j)
            HydrologyPlanetesimals.apply_insulating_boundary_conditions!(t)
            @test t[1, 2:j-1] == t[2, 2:j-1]
            @test t[i, 2:j-1] == t[i-1, 2:j-1]
            @test t[:, 1] == t[:, 2]
            @test t[:, j] == t[:, j-1]
        end
    end # testset "apply_insulating_boundary_conditions!()"

    @testset "compute_gravity_solution!()" begin
        xsize = 35_000.0
        ysize = 35_000.0
        rplanet = 12_500.0
        rcrust = 12_000.0
        Nx = 35
        Ny = 35
        sp = HydrologyPlanetesimals.StaticParameters(
            xsize=xsize,
            ysize=ysize,
            rplanet=rplanet,
            rcrust=rcrust,
            Nx=Nx,
            Ny=Ny
        )
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        dx, dy = sp.dx, sp.dy
        xsize, ysize = sp.xsize, sp.ysize
        G = sp.G
        # P nodes
        xp=-dx/2:dx:xsize+dx/2
        yp=-dy/2:dy:ysize+dy/2
        # LP = ExtendableSparseMatrix(Nx1*Ny1, Nx1*Ny1)
        SP = zeros(Float64, Nx1*Ny1)
        RP = zeros(Float64, Nx1*Ny1)
        gx = zeros(Float64, Ny1, Nx1)
        gy = zeros(Float64, Ny1, Nx1)
        LP_ver = zeros(Nx1*Ny1, Nx1*Ny1)
        RP_ver = zeros(Float64, Nx1*Ny1)
        FI_ver = zeros(Float64, Ny1, Nx1)
        gx_ver = zeros(Float64, Ny1, Nx1)
        gy_ver = zeros(Float64, Ny1, Nx1)
        # simulate density field RHO
        RHO = rand(1:0.1:7000, Ny1, Nx1)
        # compute gravity solution
        HydrologyPlanetesimals.compute_gravity_solution!(
            SP,
            RP,
            RHO,
            xp,
            yp,
            gx,
            gy,
            sp
        )
        # verification, from madcph.m, lines 680ff
        for j=1:1:Nx1
            for i=1:1:Ny1
                # Define global index in algebraic space
                gk=(j-1)*Ny1+i
                # Distance from the model centre
                rnode=((xp[j]-xsize/2)^2+(yp[i]-ysize/2)^2)^0.5
                # External points
                if rnode>xsize/2 || i==1 || i==Ny1 || j==1 || j==Nx1
                    # Boundary Condition
                    # PHI=0
                    LP_ver[gk,gk]=1; # Left part
                    RP_ver[gk]=0; # Right part
                else
                    # Internal points: Temperature eq.
                    # d2PHI/dx^2+d2PHI/dy^2=2/3*4*G*pi*RHO
                    #          PHI2
                    #           |
                    #           |
                    #  PHI1----PHI3----PHI5
                    #           |
                    #           |
                    #          PHI4
                    #
                    # Density gradients
                    dRHOdx=(RHO[i,j+1]-RHO[i,j-1])/2/dx
                    dRHOdy=(RHO[i+1,j]-RHO[i-1,j])/2/dy
                    # Left part
                    LP_ver[gk,gk-Ny1]=1/dx^2; # PHI1
                    LP_ver[gk,gk-1]=1/dy^2; # PHI2
                    LP_ver[gk,gk]=-2/dx^2-2/dy^2; # PHI3
                    LP_ver[gk,gk+1]=1/dy^2; # PHI4
                    LP_ver[gk,gk+Ny1]=1/dx^2; # PHI5
                    # Right part
                    RP_ver[gk]=2/3*4*G*pi*RHO[i,j]
                end
            end
        end
        # Solving matrixes
        SP_ver=LP_ver\RP_ver # Obtaining algebraic vector of solutions SP[]
        # Reload solutions SP[] to geometrical array PHI[]
        # Going through all grid points
        for j=1:1:Nx1
            for i=1:1:Ny1
                # Compute global index
                gk=(j-1)*Ny1+i
                # Reload solution
                FI_ver[i,j]=SP_ver[gk]
            end
        end
        # Compute gravity acceleration
        # gx
        for j=1:1:Nx
            for i=1:1:Ny1
                # gx=-dPHI/dx
                gx_ver[i,j]=-(FI_ver[i,j+1]-FI_ver[i,j])/dx
            end
        end
        # gy
        for j=1:1:Nx1
            for i=1:1:Ny
                # gy=-dPHI/dy
                gy_ver[i,j]=-(FI_ver[i+1,j]-FI_ver[i,j])/dy
            end
        end
        # test
        for j=1:1:Nx, i=1:1:Ny
            @test gx[i, j] ≈ gx_ver[i, j] rtol=1e-6
            @test gy[i, j] ≈ gy_ver[i, j] rtol=1e-6
        end
    end # testset "compute_gravity_solution!()"

    @testset "recompute_bulk_viscosity!()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1 
        etaphikoef = sp.etaphikoef
        ETAP = zeros(Ny1, Nx1)
        ETAPHI = zeros(Ny1, Nx1)
        ETAP_ver = zeros(Ny1, Nx1)
        ETAPHI_ver = zeros(Ny1, Nx1)
        # simulate data
        ETA = rand(Ny, Nx)
        PHI = rand(Ny1, Nx1)
        # compute bulk viscosity
        HydrologyPlanetesimals.recompute_bulk_viscosity!(
            ETA,
            ETAP,
            ETAPHI,
            PHI,
            etaphikoef
        )
        # verification, from madcph.m, lines 771ff
        for i=2:1:Ny
            for j=2:1:Nx
                ETAP_ver[i,j]=1/((1/ETA[i-1,j-1]+1/ETA[i,j-1]+1/ETA[i-1,j]+1/ETA[i,j])/4)
                ETAPHI_ver[i,j]=etaphikoef*ETAP_ver[i,j]/PHI[i,j]
            end
        end       
        # test
        for j=1:1:Nx, i=1:1:Ny
            @test ETAP[i, j] ≈ ETAP_ver[i, j] rtol=1e-6
            @test ETAPHI[i, j] ≈ ETAPHI_ver[i, j] rtol=1e-6
        end
    end # testset "recompute_bulk_viscosity!()"

    @testset "get_viscosities_stresses_density_gradients()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        dx, dy, dt = sp.dx, sp.dy, sp.dtelastic
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        # simulate data
        ETA = rand(Ny, Nx)
        ETAP = rand(Ny1, Nx1)
        GGG = rand(Ny, Nx)
        GGGP = rand(Ny1, Nx1)
        SXY0 = rand(Ny, Nx)
        SXX0 = rand(Ny1, Nx1)
        RHOX = rand(Ny1, Nx1)
        RHOY = rand(Ny1, Nx1)
        ETAcomp = zeros(Ny, Nx)
        ETAPcomp = zeros(Ny1, Nx1)
        SXYcomp = zeros(Ny, Nx)
        SXXcomp = zeros(Ny1, Nx1)
        SYYcomp = zeros(Ny1, Nx1)
        dRHOXdx = zeros(Ny1, Nx1)
        dRHOXdy = zeros(Ny1, Nx1)
        dRHOYdx = zeros(Ny1, Nx1)
        dRHOYdy = zeros(Ny1, Nx1)
        # compute viscosities, stresses, density gradients
        HydrologyPlanetesimals.get_viscosities_stresses_density_gradients!(
            ETA,
            ETAP,
            GGG,
            GGGP,
            SXY0,
            SXX0,
            RHOX,
            RHOY,
            dx,
            dy,
            dt,
            Nx,
            Ny,
            Nx1,
            Ny1,
            ETAcomp,
            ETAPcomp,
            SXYcomp,
            SXXcomp,
            SYYcomp,
            dRHOXdx,
            dRHOXdy,
            dRHOYdx,
            dRHOYdy
        )
        # verification, from madcph.m, lines 832ff, 905ff
        for j=1:1:Nx, i=1:1:Ny
            # x-Stokes
            if i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1
                # pass: external points
            else
                # x-Stokes internal points
                # Computational viscosity
                ETA1=ETA[i-1,j]*GGG[i-1,j]*dt/(GGG[i-1,j]*dt+ETA[i-1,j])
                ETA2=ETA[i,j]*GGG[i,j]*dt/(GGG[i,j]*dt+ETA[i,j])
                ETAP1=ETAP[i,j]*GGGP[i,j]*dt/(GGGP[i,j]*dt+ETAP[i,j])
                ETAP2=ETAP[i,j+1]*GGGP[i,j+1]*dt/(GGGP[i,j+1]*dt+ETAP[i,j+1])
                # Old stresses
                SXY1=SXY0[i-1,j]*ETA[i-1,j]/(GGG[i-1,j]*dt+ETA[i-1,j])
                SXY2=SXY0[i,j]*ETA[i,j]/(GGG[i,j]*dt+ETA[i,j])
                SXX1=SXX0[i,j]*ETAP[i,j]/(GGGP[i,j]*dt+ETAP[i,j])
                SXX2=SXX0[i,j+1]*ETAP[i,j+1]/(GGGP[i,j+1]*dt+ETAP[i,j+1])
                # Density gradients
                dRHOdx=(RHOX[i,j+1]-RHOX[i,j-1])/2/dx
                dRHOdy=(RHOX[i+1,j]-RHOX[i-1,j])/2/dy
                # test
                @test ETAcomp[i-1, j] ≈ ETA1 rtol=1e-6
                @test ETAcomp[i, j] ≈ ETA2 rtol=1e-6
                @test ETAPcomp[i, j] ≈ ETAP1 rtol=1e-6
                @test ETAPcomp[i, j+1] ≈ ETAP2 rtol=1e-6
                @test SXYcomp[i-1, j] ≈ SXY1 rtol=1e-6
                @test SXYcomp[i, j] ≈ SXY2 rtol=1e-6
                @test SXXcomp[i, j] ≈ SXX1 rtol=1e-6
                @test SXXcomp[i, j+1] ≈ SXX2 rtol=1e-6
                @test dRHOXdx[i, j] ≈ dRHOdx rtol=1e-6
                @test dRHOXdy[i, j] ≈ dRHOdy rtol=1e-6        
            end
            # y-Stokes
            if j==1 || j==Nx1 || i==1 || i==Ny || i==Ny1
                # pass: external points
            else
                # Computational viscosity
                ETA1=ETA[i,j-1]*GGG[i,j-1]*dt/(GGG[i,j-1]*dt+ETA[i,j-1])
                ETA2=ETA[i,j]*GGG[i,j]*dt/(GGG[i,j]*dt+ETA[i,j])
                ETAP1=ETAP[i,j]*GGGP[i,j]*dt/(GGGP[i,j]*dt+ETAP[i,j])
                ETAP2=ETAP[i+1,j]*GGGP[i+1,j]*dt/(GGGP[i+1,j]*dt+ETAP[i+1,j])
                # Old stresses
                SXY1=SXY0[i,j-1]*ETA[i,j-1]/(GGG[i,j-1]*dt+ETA[i,j-1])
                SXY2=SXY0[i,j]*ETA[i,j]/(GGG[i,j]*dt+ETA[i,j])
                SYY1=-SXX0[i,j]*ETAP[i,j]/(GGGP[i,j]*dt+ETAP[i,j])
                SYY2=-SXX0[i+1,j]*ETAP[i+1,j]/(GGGP[i+1,j]*dt+ETAP[i+1,j])
                # Density gradients
                dRHOdx=(RHOY[i,j+1]-RHOY[i,j-1])/2/dx
                dRHOdy=(RHOY[i+1,j]-RHOY[i-1,j])/2/dy
                # test
                @test ETAcomp[i, j-1] ≈ ETA1 rtol=1e-6
                @test ETAcomp[i, j] ≈ ETA2 rtol=1e-6
                @test ETAPcomp[i, j] ≈ ETAP1 rtol=1e-6
                @test ETAPcomp[i+1, j] ≈ ETAP2 rtol=1e-6
                @test SXYcomp[i, j-1] ≈ SXY1 rtol=1e-6
                @test SXYcomp[i, j] ≈ SXY2 rtol=1e-6
                @test SYYcomp[i, j] ≈ SYY1 rtol=1e-6
                @test SYYcomp[i+1, j] ≈ SYY2 rtol=1e-6
                @test dRHOYdx[i, j] ≈ dRHOdx rtol=1e-6
                @test dRHOYdy[i, j] ≈ dRHOdy rtol=1e-6     
            end       
        end
    end # testset "get_viscosities_stresses_density_gradients()"

    @testset "setup_hydromechanical_lse()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx1, Ny1 = sp.Nx1, sp.Ny1        
        # setup hydromechanical LSE
        L, R, S = HydrologyPlanetesimals.setup_hydromechanical_lse(sp)
        # test
        @test typeof(L) == ExtendableSparseMatrix{Float64, Int64}
        @test size(L) == (Nx1*Ny1*6, Nx1*Ny1*6)
        @test typeof(R) == Vector{Float64}
        @test size(R) == (Nx1*Ny1*6,)
        @test typeof(S) == Vector{Float64}
        @test size(S) == (Nx1*Ny1*6,)
    end

    @testset "setup_thermal_lse()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx1, Ny1 = sp.Nx1, sp.Ny1        
        # setup thermal LSE
        LP, RP, SP = HydrologyPlanetesimals.setup_thermal_lse(sp)
        # test
        @test typeof(LP) == ExtendableSparseMatrix{Float64, Int64}
        @test size(LP) == (Nx1*Ny1, Nx1*Ny1)
        @test typeof(RP) == Vector{Float64}
        @test size(RP) == (Nx1*Ny1,)
        @test typeof(SP) == Vector{Float64}
        @test size(SP) == (Nx1*Ny1,)
    end
    
    @testset "setup_gravitational_lse()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx1, Ny1 = sp.Nx1, sp.Ny1        
        # setup gravitational LSE
        LT, RT, ST = HydrologyPlanetesimals.setup_gravitational_lse(sp)
        # test
        @test typeof(LT) == ExtendableSparseMatrix{Float64, Int64}
        @test size(LT) == (Nx1*Ny1, Nx1*Ny1)
        @test typeof(RT) == Vector{Float64}
        @test size(RT) == (Nx1*Ny1,)
        @test typeof(ST) == Vector{Float64}
        @test size(ST) == (Nx1*Ny1,)
    end
    
    @testset "assemble_hydromechanical_lse()" begin
        xsize = 35_000.0
        ysize = 35_000.0
        rplanet = 12_500.0
        rcrust = 12_000.0
        Nx = 35
        Ny = 35
        sp = HydrologyPlanetesimals.StaticParameters(
            xsize=xsize,
            ysize=ysize,
            rplanet=rplanet,
            rcrust=rcrust,
            Nx=Nx,
            Ny=Ny
        )
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        dx, dy = sp.dx, sp.dy
        xsize, ysize = sp.xsize, sp.ysize
        vxleft = sp.vxleft
        vxright = sp.vxright
        vytop = sp.vytop
        vybottom = sp.vybottom
        bctop = sp.bctop
        bcbottom = sp.bcbottom
        bcleft = sp.bcleft
        bcright = sp.bcright
        bcftop = sp.bcftop
        bcfbottom = sp.bcfbottom
        bcfleft = sp.bcfleft
        bcfright = sp.bcfright
        pscale  = sp.pscale
        psurface = sp.psurface
        etaphikoef = sp.etaphikoef
        dt = sp.dtelastic
        # simulate data
        ETA = rand(Ny, Nx)
        ETAP = rand(Ny1, Nx1)
        GGG = rand(Ny, Nx)
        GGGP = rand(Ny1, Nx1)
        SXY0 = rand(Ny, Nx)
        SXX0 = rand(Ny, Nx)
        RHOX = rand(Ny1, Nx1)
        RHOY = rand(Ny1, Nx1)
        RHOFX = rand(Ny1, Nx1)
        RHOFY = rand(Ny1, Nx1)
        RX = rand(Ny1, Nx1)
        RY = rand(Ny1, Nx1)
        ETAPHI = rand(Ny1, Nx1)
        BETTAPHI = rand(Ny1, Nx1)
        PHI = rand(Ny1, Nx1)
        gx = rand(Ny1, Nx1)
        gy = rand(Ny1, Nx1) 
        pr0 = rand(Ny1, Nx1)
        pf0 = rand(Ny1, Nx1)
        ETAcomp = zeros(Ny, Nx)
        ETAPcomp = zeros(Ny1, Nx1)
        SXYcomp = zeros(Ny, Nx)
        SXXcomp = zeros(Ny, Nx)
        SYYcomp = zeros(Ny, Nx)
        dRHOXdx = zeros(Ny1, Nx1)
        dRHOXdy = zeros(Ny1, Nx1)
        dRHOYdx = zeros(Ny1, Nx1)
        dRHOYdy = zeros(Ny1, Nx1)
        # LSE
        L = ExtendableSparseMatrix(Nx1*Ny1*6, Nx1*Ny1*6)
        R = zeros(Nx1*Ny1*6)
        L_ver = zeros(Nx1*Ny1*6, Nx1*Ny1*6)
        R_ver = zeros(Nx1*Ny1*6)
        # assemble hydromechanical LSE
        HydrologyPlanetesimals.assemble_hydromechanical_lse!(
            ETAcomp,
            ETAPcomp,
            SXYcomp,
            SXXcomp,
            SYYcomp,
            dRHOXdx,
            dRHOXdy,
            dRHOYdx,
            dRHOYdy,
            RHOX,
            RHOY,
            RHOFX,
            RHOFY,
            RX,
            RY,
            ETAPHI,
            BETTAPHI,
            PHI,
            gx,
            gy,
            pr0,
            pf0,
            dt,
            L,
            R,
            sp
        )
        # verification
        # verification, from madcph.m, lines 779ff
        # Hydro-Mechanical Solution
        # Composing global matrixes L_ver[], R_ver[] for Stokes & continuity equations
        for j=1:1:Nx1
            for i=1:1:Ny1
                # Define global indexes in algebraic space
                kvx=((j-1)*Ny1+i-1)*6+1; # Vx solid
                kvy=kvx+1; # Vy solid
                kpm=kvx+2; # Ptotal
                kqx=kvx+3; # qx Darcy
                kqy=kvx+4; # qy Darcy
                kpf=kvx+5; # P fluid
                
                # Vx equation External points
                if i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1
                    # Boundary Condition 
                    # Ghost unknowns 1*Vx=0
                    if j==Nx1
                        L_ver[kvx,kvx]=1; # Left part
                        R_ver[kvx]=0; # Right part
                    end
                    # Left Boundary
                    if j==1
                        L_ver[kvx,kvx]=1; # Left part
                        R_ver[kvx]=vxleft; # Right part
                    end
                    # Right Boundary
                    if j==Nx 
                        L_ver[kvx,kvx]=1; # Left part
                        R_ver[kvx]=vxright; # Right part
                    end
                    # Top boundary
                    if i==1 && j>1 && j<Nx
                        L_ver[kvx,kvx]=1; # Left part
                        L_ver[kvx,kvx+6]=bctop; # Left part
                        R_ver[kvx]=0; # Right part
                    end
                    # Top boundary
                    if i==Ny1 && j>1 && j<Nx
                        L_ver[kvx,kvx]=1; # Left part
                        L_ver[kvx,kvx-6]=bcbottom; # Left part
                        R_ver[kvx]=0; # Right part
                    end
                else
                # Internal points: x-Stokes eq.
                #            Vx2
                #             |
                #        Vy1  |  Vy3
                #             |
                #     Vx1-P1-Vx3-P2-Vx5
                #             |
                #        Vy2  |  Vy4
                #             |
                #            Vx4
                #
                # Computational viscosity
                ETA1=ETA[i-1,j]*GGG[i-1,j]*dt/(GGG[i-1,j]*dt+ETA[i-1,j])
                ETA2=ETA[i,j]*GGG[i,j]*dt/(GGG[i,j]*dt+ETA[i,j])
                ETAP1=ETAP[i,j]*GGGP[i,j]*dt/(GGGP[i,j]*dt+ETAP[i,j])
                ETAP2=ETAP[i,j+1]*GGGP[i,j+1]*dt/(GGGP[i,j+1]*dt+ETAP[i,j+1])
                # Old stresses
                SXY1=SXY0[i-1,j]*ETA[i-1,j]/(GGG[i-1,j]*dt+ETA[i-1,j])
                SXY2=SXY0[i,j]*ETA[i,j]/(GGG[i,j]*dt+ETA[i,j])
                SXX1=SXX0[i,j]*ETAP[i,j]/(GGGP[i,j]*dt+ETAP[i,j])
                SXX2=SXX0[i,j+1]*ETAP[i,j+1]/(GGGP[i,j+1]*dt+ETAP[i,j+1])
                # Density gradients
                dRHOdx=(RHOX[i,j+1]-RHOX[i,j-1])/2/dx
                dRHOdy=(RHOX[i+1,j]-RHOX[i-1,j])/2/dy
                # Left part
                L_ver[kvx,kvx-Ny1*6]=ETAP1/dx^2; # Vx1
                L_ver[kvx,kvx-6]=ETA1/dy^2; # Vx2
                L_ver[kvx,kvx]=-(ETAP1+ETAP2)/dx^2-  (ETA1+ETA2)/dy^2-  dRHOdx*gx[i,j]*dt; # Vx3
                L_ver[kvx,kvx+6]=ETA2/dy^2; # Vx4
                L_ver[kvx,kvx+Ny1*6]=ETAP2/dx^2; # Vx5
                L_ver[kvx,kvy]=ETAP1/dx/dy-ETA2/dx/dy-dRHOdy*gx[i,j]*dt/4;  # Vy2
                L_ver[kvx,kvy+Ny1*6]=-ETAP2/dx/dy+ETA2/dx/dy-dRHOdy*gx[i,j]*dt/4;  # Vy4
                L_ver[kvx,kvy-6]=-ETAP1/dx/dy+ETA1/dx/dy-dRHOdy*gx[i,j]*dt/4;  # Vy1
                L_ver[kvx,kvy+Ny1*6-6]=ETAP2/dx/dy-ETA1/dx/dy-dRHOdy*gx[i,j]*dt/4;  # Vy3
                L_ver[kvx,kpm]=pscale/dx; # P1
                L_ver[kvx,kpm+Ny1*6]=-pscale/dx; # P2
                # Right part
                R_ver[kvx]=-RHOX[i,j]*gx[i,j]-(SXY2-SXY1)/dy-(SXX2-SXX1)/dx
                end
                
                # Vy equation External points
                if j==1 || j==Nx1 || i==1 || i==Ny || i==Ny1
                    # Boundary Condition
                    # Ghost unknowns 1*Vx=0
                    if i==Ny1
                        L_ver[kvy,kvy]=1; # Left part
                        R_ver[kvy]=0; # Right part
                    end
                    # Top boundary
                    if i==1
                        L_ver[kvy,kvy]=1; # Left part
                        R_ver[kvy]=vytop; # Right part
                    end
                    # Bottom boundary
                    if i==Ny
                        L_ver[kvy,kvy]=1; # Left part
                        R_ver[kvy]=vybottom; # Right part
                    end
                    # Left boundary
                    if j==1 && i>1 && i<Ny
                        L_ver[kvy,kvy]=1; # Left part
                        L_ver[kvy,kvy+6*Ny1]=bcleft; # Left part
                        R_ver[kvy]=0; # Right part
                    end
                    # Right boundary
                    if j==Nx1 && i>1 && i<Ny
                        L_ver[kvy,kvy]=1; # Left part
                        L_ver[kvy,kvy-6*Ny1]=bcright; # Left part
                        R_ver[kvy]=0; # Right part
                    end
                else
                # Internal points: y-Stokes eq.
                #            Vy2
                #             |
                #         Vx1 P1 Vx3
                #             |
                #     Vy1----Vy3----Vy5
                #             |
                #         Vx2 P2 Vx4
                #             |
                #            Vy4
                #
                # Computational viscosity
                ETA1=ETA[i,j-1]*GGG[i,j-1]*dt/(GGG[i,j-1]*dt+ETA[i,j-1])
                ETA2=ETA[i,j]*GGG[i,j]*dt/(GGG[i,j]*dt+ETA[i,j])
                ETAP1=ETAP[i,j]*GGGP[i,j]*dt/(GGGP[i,j]*dt+ETAP[i,j])
                ETAP2=ETAP[i+1,j]*GGGP[i+1,j]*dt/(GGGP[i+1,j]*dt+ETAP[i+1,j])
                # Old stresses
                SXY1=SXY0[i,j-1]*ETA[i,j-1]/(GGG[i,j-1]*dt+ETA[i,j-1])
                SXY2=SXY0[i,j]*ETA[i,j]/(GGG[i,j]*dt+ETA[i,j])
                SYY1=-SXX0[i,j]*ETAP[i,j]/(GGGP[i,j]*dt+ETAP[i,j])
                SYY2=-SXX0[i+1,j]*ETAP[i+1,j]/(GGGP[i+1,j]*dt+ETAP[i+1,j])
                # Density gradients
                dRHOdx=(RHOY[i,j+1]-RHOY[i,j-1])/2/dx
                dRHOdy=(RHOY[i+1,j]-RHOY[i-1,j])/2/dy
                # Left part
                L_ver[kvy,kvy-Ny1*6]=ETA1/dx^2; # Vy1
                L_ver[kvy,kvy-6]=ETAP1/dy^2; # Vy2
                L_ver[kvy,kvy]=-(ETAP1+ETAP2)/dy^2-  (ETA1+ETA2)/dx^2-  dRHOdy*gy[i,j]*dt; # Vy3
                L_ver[kvy,kvy+6]=ETAP2/dy^2; # Vy4
                L_ver[kvy,kvy+Ny1*6]=ETA2/dx^2; # Vy5
                L_ver[kvy,kvx]=ETAP1/dx/dy-ETA2/dx/dy-dRHOdx*gy[i,j]*dt/4; #Vx3
                L_ver[kvy,kvx+6]=-ETAP2/dx/dy+ETA2/dx/dy-dRHOdx*gy[i,j]*dt/4; #Vx4
                L_ver[kvy,kvx-Ny1*6]=-ETAP1/dx/dy+ETA1/dx/dy-dRHOdx*gy[i,j]*dt/4; #Vx1
                L_ver[kvy,kvx+6-Ny1*6]=ETAP2/dx/dy-ETA1/dx/dy-dRHOdx*gy[i,j]*dt/4; #Vx2
                L_ver[kvy,kpm]=pscale/dy; # P1
                L_ver[kvy,kpm+6]=-pscale/dy; # P2
                
                # Right part
                R_ver[kvy]=-RHOY[i,j]*gy[i,j]-(SXY2-SXY1)/dx-(SYY2-SYY1)/dy
                end
                
                # P equation External points
                if i==1 || j==1 || i==Ny1 || j==Nx1
                    # Boundary Condition
                    # 1*P=0
                    L_ver[kpm,kpm]=1; # Left part
                    R_ver[kpm]=0; # Right part
                else
                # Internal points: continuity eq.
                # dVx/dx+dVy/dy=0
                #            Vy1
                #             |
                #        Vx1--P--Vx2
                #             |
                #            Vy2
                #
                # Left part
                L_ver[kpm,kvx-Ny1*6]=-1/dx; # Vx1
                L_ver[kpm,kvx]=1/dx; # Vx2
                L_ver[kpm,kvy-6]=-1/dy; # Vy1
                L_ver[kpm,kvy]=1/dy; # Vy2
                L_ver[kpm,kpm]= pscale/(1-PHI[i,j])*(1/ETAPHI[i,j]+BETTAPHI[i,j]/dt); # Ptotal
                L_ver[kpm,kpf]=-pscale/(1-PHI[i,j])*(1/ETAPHI[i,j]+BETTAPHI[i,j]/dt); # Pfluid
                # Right part
                R_ver[kpm]=(pr0[i,j]-pf0[i,j])/(1-PHI[i,j])*BETTAPHI[i,j]/dt
                end

                # qxDarcy equation External points
                if i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1
                    # Boundary Condition
                    # 1*qx=0
                    L_ver[kqx,kqx]=1; # Left part
                    R_ver[kqx]=0; # Right part
                    # Top boundary
                    if i==1 && j>1 && j<Nx
                        L_ver[kqx,kqx+6]=bcftop; # Left part
                    end
                    # Bottom boundary
                    if i==Ny1 && j>1 && j<Nx
                        L_ver[kqx,kqx-6]=bcfbottom; # Left part
                    end
                else
                # Internal points: x-Darcy eq.
                # Rx*qxDarcy+dP/dx=RHOfluid*gx
                #     P1-qxD-P2
                # Left part
                L_ver[kqx,kqx]=RX[i,j]; # qxD
                L_ver[kqx,kpf]=-pscale/dx; # P1
                L_ver[kqx,kpf+Ny1*6]=pscale/dx; # P2
                # Right part
                R_ver[kqx]=RHOFX[i,j]*gx[i,j]
                end
                
                # qyDarcy equation External points
                if j==1 || j==Nx1 || i==1 || i==Ny || i==Ny1
                    # Boundary Condition
                    # 1*Vy=0
                    L_ver[kqy,kqy]=1; # Left part
                    R_ver[kqy]=0; # Right part
                    # Left boundary
                    if j==1 && i>1 && i<Ny
                        L_ver[kqy,kqy+6*Ny1]=bcfleft; # Left part
                    end
                    # Right boundary
                    if j==Nx1 && i>1 && i<Ny
                        L_ver[kqy,kqy-6*Ny1]=bcfright; # Left part
                    end
                else
                # Internal points: y-Stokes eq.
                # Internal points: x-Darcy eq.
                # Rx*qxDarcy+dP/dx=RHOfluid*gx
                #      P1
                #      |
                #     qxD
                #      |
                #      P2
                # Left part
                L_ver[kqy,kqy]=RY[i,j]; # qxD
                L_ver[kqy,kpf]=-pscale/dy; # P1
                L_ver[kqy,kpf+6]=pscale/dy; # P
                # Right part
                R_ver[kqy]=RHOFY[i,j]*gy[i,j]
                end
                
                # Pfluid equation External points
                if i==1 || j==1 || i==Ny1 || j==Nx1 || (i==2 && j==2)
                    # Boundary Condition
                    # 1*Pfluid=0
                    L_ver[kpf,kpf]=1; # Left part
                    R_ver[kpf]=0; # Right part
                    # Real BC
                    if i==2 && j==2
                        L_ver[kpf,kpf]=1*pscale; #Left part
                        R_ver[kpf]=psurface; # Right part
                    end
                else
                # Internal points: continuity eq.
                # dqxD/dx+dqyD/dy-(Ptotal-Pfluid)/ETHAphi=0
                #            qyD1
                #              |
                #        qxD1--P--qxD2
                #              |
                #            qyD2
                #
                # Left part
                L_ver[kpf,kqx-Ny1*6]=-1/dx; # qxD1
                L_ver[kpf,kqx]=1/dx; # qxD2
                L_ver[kpf,kqy-6]=-1/dy; # qyD1
                L_ver[kpf,kqy]=1/dy; # qyD2
                L_ver[kpf,kpm]=-pscale/(1-PHI[i,j])*(1/ETAPHI[i,j]+BETTAPHI[i,j]/dt); # Ptotal
                L_ver[kpf,kpf]= pscale/(1-PHI[i,j])*(1/ETAPHI[i,j]+BETTAPHI[i,j]/dt); # Pfluid
                # Right part
                R_ver[kpf]=-(pr0[i,j]-pf0[i,j])/(1-PHI[i,j])*BETTAPHI[i,j]/dt
                end
            end
        end
        # test
        for j=1:1:Nx1*6, i=1:1:Ny1*6
            @test L[i, j] ≈ L_ver[i, j] rtol=1e-6
            @test R[i] ≈ R_ver[i] rtol=1e-6
        end
    end # testset "assemble_hydromechanical_lse()"

    @testset "process_hydromechanical_solution!()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        pscale = sp.pscale
        # simulate data
        S = rand(Nx1*Ny1*6)
        vx = zeros(Ny1, Nx1)
        vy = zeros(Ny1, Nx1)
        pr = zeros(Ny1, Nx1)
        qxD = zeros(Ny1, Nx1)
        qyD = zeros(Ny1, Nx1)
        pf = zeros(Ny1, Nx1)
        vx_ver = zeros(Ny1, Nx1)
        vy_ver = zeros(Ny1, Nx1)
        pr_ver = zeros(Ny1, Nx1)
        qxD_ver = zeros(Ny1, Nx1)
        qyD_ver = zeros(Ny1, Nx1)
        pf_ver = zeros(Ny1, Nx1)
        # process solution
        HydrologyPlanetesimals.process_hydromechanical_solution!(
            S,
            vx,
            vy,
            pr,
            qxD,
            qyD,
            pf,
            pscale,
            Nx1,
            Ny1
        )
        # verification, from madcph.m, line 1058ff
        for j=1:1:Nx1
            for i=1:1:Ny1
                # Define global indexes in algebraic space
                kvx=((j-1)*Ny1+i-1)*6+1; # Vx solid
                kvy=kvx+1; # Vy solid
                kpm=kvx+2; # Ptotal
                kqx=kvx+3; # qx Darcy
                kqy=kvx+4; # qy Darcy
                kpf=kvx+5; # P fluid
                # Reload solution
                vx_ver[i,j]=S[kvx]
                vy_ver[i,j]=S[kvy]
                pr_ver[i,j]=S[kpm]*pscale
                qxD_ver[i,j]=S[kqx]
                qyD_ver[i,j]=S[kqy]
                pf_ver[i,j]=S[kpf]*pscale
            end
        end
        # test
        for j=1:1:Nx1, i=1:1:Ny1
            @test vx[i, j] ≈ vx_ver[i, j] rtol=1e-6
            @test vy[i, j] ≈ vy_ver[i, j] rtol=1e-6
            @test pr[i, j] ≈ pr_ver[i, j] rtol=1e-6
            @test qxD[i, j] ≈ qxD_ver[i, j] rtol=1e-6
            @test qyD[i, j] ≈ qyD_ver[i, j] rtol=1e-6
            @test pf[i, j] ≈ pf_ver[i, j] rtol=1e-6
        end
    end # testset "process_hydromechanical_solution!()"

    @testset "compute_Aϕ!()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        dt = sp.dtelastic
        # simulate data
        APHI = rand(Ny1, Nx1)
        APHI_ver = rand(Ny1, Nx1)
        ETAPHI = rand(Ny1, Nx1)
        BETTAPHI = rand(Ny1, Nx1)
        PHI = rand(Ny1, Nx1)
        pr = rand(Ny1, Nx1)
        pf = rand(Ny1, Nx1)
        pr0 = rand(Ny1, Nx1)
        pf0 = rand(Ny1, Nx1)
        # compute Aϕ
        aphimax = HydrologyPlanetesimals.compute_Aϕ!(
            APHI,
            ETAPHI,
            BETTAPHI,
            PHI,
            pr,
            pf,
            pr0,
            pf0,
            dt,
            sp
        )
        # verification, from madcph.m, line 1078ff
        APHI_ver = zeros(Ny1, Nx1)
        aphimax_ver=0
        for j=2:1:Nx
            for i=2:1:Ny
                APHI_ver[i,j]=((pr[i,j]-pf[i,j])/ETAPHI[i,j]+  ((pr[i,j]-pr0[i,j])-(pf[i,j]-pf0[i,j]))/dt*BETTAPHI[i,j])/(1-PHI[i,j])/PHI[i,j]
                aphimax_ver=max(aphimax_ver,abs(APHI_ver[i,j]))
            end
        end
        # test
        for j=2:1:Nx, i=2:1:Ny
            @test APHI[i, j] ≈ APHI_ver[i, j] rtol=1e-6
        end
        @test aphimax ≈ aphimax_ver rtol=1e-6
    end # testset "compute_Aϕ!()"

    @testset "compute_fluid_velocity!()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        bcftop, bcfbottom = sp.bcftop, sp.bcfbottom
        bcfleft, bcfright = sp.bcfleft, sp.bcfright
        # simulate data
        PHIX = rand(Ny1, Nx1)
        PHIY = rand(Ny1, Nx1)
        qxD = rand(Ny1, Nx1)
        qyD = rand(Ny1, Nx1)
        vx = rand(Ny1, Nx1)
        vy = rand(Ny1, Nx1)
        vxf = zeros(Ny1, Nx1)
        vyf = zeros(Ny1, Nx1)
        vxf_ver = zeros(Ny1, Nx1)
        vyf_ver = zeros(Ny1, Nx1)
        # compute fluid velocities
        HydrologyPlanetesimals.compute_fluid_velocities!(
            PHIX,
            PHIY,
            qxD,
            qyD,
            vx,
            vy,
            vxf,
            vyf,
            sp
        )
        # verification, from madcph.m line 1090ff
        for j=1:1:Nx
            for i=2:1:Ny
                vxf_ver[i,j]=qxD[i,j]/PHIX[i,j]
            end
        end
        # Apply BC
        # Top
        vxf_ver[1,:]= -bcftop*vxf_ver[2,:];    
        # Bottom
        vxf_ver[Ny1,:]= -bcfbottom*vxf_ver[Ny,:];    
        # Vy fluid
        for j=2:1:Nx
            for i=1:1:Ny
                vyf_ver[i,j]=qyD[i,j]/PHIY[i,j]
            end
        end
        # Apply BC
        # Left
        vyf_ver[:,1]= -bcfleft*vyf_ver[:,2];    
        # Right
        vyf_ver[:, Nx1]= -bcfright*vyf_ver[:, Nx];     
        # Add solid velocity
        # vxf0=vxf; vxf=vxf+vx
        vxf_ver.=vxf_ver.+vx
        # vyf0=vyf; vyf=vyf+vy
        vyf_ver.=vyf_ver.+vy
        # test
        for j=1:1:Nx1, i=1:1:Ny1
            @test vxf[i, j] ≈ vxf_ver[i, j] rtol=1e-6
            @test vyf[i, j] ≈ vyf_ver[i, j] rtol=1e-6
        end
    end # testset "compute_fluid_velocity!()"

    @testset "compute_displacement_timestep()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        dt = sp.dtelastic
        dx, dy = sp.dx, sp.dy
        dxymax, dphimax = sp.dxymax, sp.dphimax        
        # simulate data
        aphimax = rand()
        vx = rand(Ny1, Nx1)
        vy = rand(Ny1, Nx1)
        vxf = rand(Ny1, Nx1)
        vyf = rand(Ny1, Nx1)
        # compute displacement timestep
        dtm = HydrologyPlanetesimals.compute_displacement_timestep(
            vx,
            vy,
            vxf,
            vyf,
            dt,
            dx,
            dy,
            dxymax,
            aphimax,
            dphimax
        )
        # verification, from madcph.m, line 1117ff
        dtm_ver=dt
        maxvx=maximum(abs.(vx))
        maxvy=maximum(abs.(vy))
        if dtm_ver*maxvx>dxymax*dx
            dtm_ver=dxymax*dx/maxvx
        end
        if dtm_ver*maxvy>dxymax*dy
            dtm_ver=dxymax*dy/maxvy
        end
        # Fluid velocity
        maxvxf=maximum(abs.(vxf))
        maxvyf=maximum(abs.(vyf))
        if dtm_ver*maxvxf>dxymax*dx
            dtm_ver=dxymax*dx/maxvxf
        end
        if dtm_ver*maxvyf>dxymax*dy
            dtm_ver=dxymax*dy/maxvyf
        end
        # Porosity change
        if aphimax*dtm_ver>dphimax
            dtm_ver=dphimax/aphimax
        end
        # test
        @test dtm ≈ dtm_ver rtol=1e-6
    end # testset "compute_displacement_timestep()"

    @testset "compute_stress_strainrate!()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        dtm = sp.dtelastic
        dx, dy = sp.dx, sp.dy
        # simulate data
        vx = rand(Ny1, Nx1)
        vy = rand(Ny1, Nx1)
        ETA = rand(Ny, Nx)
        GGG = rand(Ny, Nx)
        ETAP = rand(Ny1, Nx1)
        GGGP = rand(Ny1, Nx1)
        SXX0 = rand(Ny1, Nx1)
        SXY0 = rand(Ny, Nx)
        EXY = rand(Ny, Nx)
        SXY = rand(Ny, Nx)
        DSXY = rand(Ny, Nx)
        EXX = rand(Ny1, Nx1)
        SXX = rand(Ny1, Nx1)
        DSXX = rand(Ny1, Nx1)
        EII = zeros(Ny1, Nx1)
        SII = zeros(Ny1, Nx1)
        # compute stress, strainrate
        HydrologyPlanetesimals.compute_stress_strainrate!(
            vx,
            vy,
            ETA,
            GGG,
            ETAP,
            GGGP,
            SXX0,
            SXY0,
            EXX,
            EXY,
            SXX,
            SXY,
            DSXX,
            DSXY,
            EII,
            SII,
            dtm,
            sp
        )
        # verification, from madcph.m, line 1144ff
        EXY_ver = zeros(Ny, Nx); # Strain rate EPSILONxy, 1/s
        SXY_ver = zeros(Ny, Nx); # Stress SIGMAxy, Pa
        DSXY_ver = zeros(Ny, Nx); # Stress change SIGMAxy, Pa
        for j=1:1:Nx
            for i=1:1:Ny
                # EXY;SXY; DSXY
                EXY_ver[i,j]=0.5*((vx[i+1,j]-vx[i,j])/dy+(vy[i,j+1]-vy[i,j])/dx)
                SXY_ver[i,j]=2*ETA[i,j]*EXY_ver[i,j]*GGG[i,j]*dtm/(GGG[i,j]*dtm+ETA[i,j])+SXY0[i,j]*ETA[i,j]/(GGG[i,j]*dtm+ETA[i,j])
                DSXY_ver[i,j]=SXY_ver[i,j]-SXY0[i,j]
            end
        end
        # Compute EPSILONxx; SIGMA'xx in pressure nodes
        EXX_ver = zeros(Ny1, Nx1); # Strain rate EPSILONxx, 1/s
        EII_ver = zeros(Ny1, Nx1); # Second strain rate invariant, 1/s
        SXX_ver = zeros(Ny1, Nx1); # Stress SIGMA'xx, Pa
        SII_ver = zeros(Ny1, Nx1); # Second stress invariant, Pa
        DSXX_ver = zeros(Ny1, Nx1); # Stress change SIGMA'xx, Pa
        DIVV_ver = zeros(Ny1, Nx1); # div[v]
        for j=2:1:Nx
            for i=2:1:Ny
                # DIVV
                DIVV_ver[i,j]=(vx[i,j]-vx[i,j-1])/dx+(vy[i,j]-vy[i-1,j])/dy
                # EXX
                EXX_ver[i,j]=((vx[i,j]-vx[i,j-1])/dx-(vy[i,j]-vy[i-1,j])/dy)/2
                # SXX
                SXX_ver[i,j]=2*ETAP[i,j]*EXX_ver[i,j]*GGGP[i,j]*dtm/(GGGP[i,j]*dtm+ETAP[i,j])+SXX0[i,j]*ETAP[i,j]/(GGGP[i,j]*dtm+ETAP[i,j])
                DSXX_ver[i,j]=SXX_ver[i,j]-SXX0[i,j]
                # EII
                EII_ver[i,j]=(EXX_ver[i,j]^2+((EXY_ver[i,j]+EXY_ver[i-1,j]+EXY_ver[i,j-1]+EXY_ver[i-1,j-1])/4)^2)^0.5
                # SII
                SII_ver[i,j]=(SXX_ver[i,j]^2+((SXY_ver[i,j]+SXY_ver[i-1,j]+SXY_ver[i,j-1]+SXY_ver[i-1,j-1])/4)^2)^0.5
            end
        end
        # test
        for j=1:1:Nx, i=1:1:Ny
            @test EXY[i,j] ≈ EXY_ver[i,j] rtol=1e-6
            @test SXY[i,j] ≈ SXY_ver[i,j] rtol=1e-6
            @test DSXY[i,j] ≈ DSXY_ver[i,j] rtol=1e-6
        end
        for j=2:1:Nx, i=2:1:Ny
            @test EXX[i,j] ≈ EXX_ver[i,j] rtol=1e-6
            @test SXX[i,j] ≈ SXX_ver[i,j] rtol=1e-6
            # @test DSXX[i,j] ≈ DSXX_ver[i,j] rtol=1e-6
            @test EII[i,j] ≈ EII_ver[i,j] rtol=1e-6
            @test SII[i,j] ≈ SII_ver[i,j] rtol=1e-6
        end
    end # testset "compute_stress_strainrate!()"

    @testset "symmetrize_p_node_observables!()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        # simulate data
        SXX = rand(Ny1, Nx1)
        APHI = rand(Ny1, Nx1)
        PHI = rand(Ny1, Nx1)
        pr = rand(Ny1, Nx1)
        pf = rand(Ny1, Nx1)
        ps = zeros(Ny1, Nx1)
        SXX_ver = copy(SXX)
        APHI_ver = copy(APHI)
        PHI_ver = copy(PHI)
        pr_ver = copy(pr)
        pf_ver = copy(pf)
        ps_ver = zeros(Ny1, Nx1)
        # symmetrize p node variables
        HydrologyPlanetesimals.symmetrize_p_node_observables!(
            SXX,
            APHI,
            PHI,
            pr,
            pf,
            ps,
            Nx,
            Ny,
            Nx1,
            Ny1
        )
        # verification, from madcph.m, line 1196ff
        # Apply Symmetry to Pressure nodes
        # External P-nodes: symmetry
        # Top
        SXX_ver[1,2:Nx]=SXX_ver[2,2:Nx]
        APHI_ver[1,2:Nx]=APHI_ver[2,2:Nx];    
        PHI_ver[1,2:Nx]=PHI_ver[2,2:Nx];    
        pr_ver[1,2:Nx]=pr_ver[2,2:Nx];    
        pf_ver[1,2:Nx]=pf_ver[2,2:Nx];    
        # Bottom
        SXX_ver[Ny1,2:Nx]=SXX_ver[Ny,2:Nx]
        APHI_ver[Ny1,2:Nx]=APHI_ver[Ny,2:Nx];    
        PHI_ver[Ny1,2:Nx]=PHI_ver[Ny,2:Nx];    
        pr_ver[Ny1,2:Nx]=pr_ver[Ny,2:Nx];    
        pf_ver[Ny1,2:Nx]=pf_ver[Ny,2:Nx];    
        # Left
        SXX_ver[:,1]=SXX_ver[:,2]
        APHI_ver[:,1]=APHI_ver[:,2];    
        PHI_ver[:,1]=PHI_ver[:,2];    
        pr_ver[:,1]=pr_ver[:,2];    
        pf_ver[:,1]=pf_ver[:,2];    
        # Right
        SXX_ver[:, Nx1]=SXX_ver[:, Nx]
        APHI_ver[:, Nx1]=APHI_ver[:, Nx];    
        PHI_ver[:, Nx1]=PHI_ver[:, Nx];    
        pr_ver[:, Nx1]=pr_ver[:, Nx];    
        pf_ver[:, Nx1]=pf_ver[:, Nx]; 
        # Compute solid pressure
        ps_ver=(pr_ver .- pf_ver.*PHI_ver)./(1 .- PHI_ver)
        # test
        @test SXX[1, 2:Nx] == SXX_ver[1, 2:Nx]
        @test APHI[1, 2:Nx] == APHI_ver[1, 2:Nx]
        @test PHI[1, 2:Nx] == PHI_ver[1, 2:Nx]
        @test pr[1, 2:Nx] == pr_ver[1, 2:Nx]
        @test pf[1, 2:Nx] == pf_ver[1, 2:Nx]
        @test SXX[Ny1, 2:Nx] == SXX_ver[Ny1, 2:Nx]
        @test APHI[Ny1, 2:Nx] == APHI_ver[Ny1, 2:Nx]
        @test PHI[Ny1, 2:Nx] == PHI_ver[Ny1, 2:Nx]
        @test pr[Ny1, 2:Nx] == pr_ver[Ny1, 2:Nx]
        @test pf[Ny1, 2:Nx] == pf_ver[Ny1, 2:Nx]
        @test SXX[:, 1] == SXX_ver[:, 1]
        @test APHI[:, 1] == APHI_ver[:, 1]
        @test PHI[:, 1] == PHI_ver[:, 1]
        @test pr[:, 1] == pr_ver[:, 1]
        @test pf[:, 1] == pf_ver[:, 1]
        @test SXX[:, Nx1] == SXX_ver[:, Nx1]
        @test APHI[:, Nx1] == APHI_ver[:, Nx1]
        @test PHI[:, Nx1] == PHI_ver[:, Nx1]
        @test pr[:, Nx1] == pr_ver[:, Nx1]
        @test pf[:, Nx1] == pf_ver[:, Nx1]
        @test ps == ps_ver
    end # testset "symmetrize_p_node_observables!()"

    @testset "positive_max()" begin
        # simulate data
        A = rand(-100:0.1:100, 1000, 1000)
        B = rand(-100:0.1:100, 1000, 1000)
        R = zeros(1000, 1000)
        # compute positive max
        HydrologyPlanetesimals.positive_max!(A, B, R)
        # test
        for i in eachindex(R)
            @test R[i] == max(A[i], B[i], 0.0)
        end
    end # testset "positive_max()"

    @testset "compute_nodal_adjustment!()" begin
        sp = HydrologyPlanetesimals.StaticParameters(
            etamin=0.2,
            etamax=0.4,
            etawt=0.5,
            nplast=1
        )
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        dt = sp.dtelastic
        nplast = sp.nplast
        yerrmax = sp.yerrmax
        etamin, etamax, etawt = sp.etamin, sp.etamax, sp.etawt
        iplast = 1
        # simulate data
        ETA = rand(Ny, Nx)
        ETA0 = rand(Ny, Nx)
        ETA5 = zeros(Ny, Nx)
        GGG = rand(Ny, Nx)
        SXX = rand(Ny1, Nx1)
        SXY = rand(Ny, Nx)
        pr = rand(Ny1, Nx1)
        pf = rand(Ny1, Nx1)
        COH = rand(Ny, Nx)
        TEN = rand(Ny, Nx)
        FRI = rand(Ny, Nx)
        SIIB = zeros(Ny, Nx)
        siiel = zeros(Ny, Nx)
        prB = zeros(Ny, Nx)
        pfB = zeros(Ny, Nx)
        syieldc = zeros(Ny, Nx)
        syieldt = zeros(Ny, Nx)
        syield = zeros(Ny, Nx)
        etapl = zeros(Ny, Nx)
        YNY = zeros(Bool, Ny, Nx)
        YNY5 = zeros(Bool, Ny, Nx)
        DSY = rand(Ny, Nx)
        DDD = rand(Ny, Nx)
        YNPL = zeros(Bool, Ny, Nx)
        YERRNOD = zeros(nplast)
        YERRNOD_ver = zeros(nplast)
        # compute nodal adjustment
        complete = HydrologyPlanetesimals.compute_nodal_adjustment!(
            ETA,
            ETA0,
            ETA5,
            GGG,
            SXX,
            SXY,
            pr,
            pf,
            COH,
            TEN,
            FRI,
            SIIB,
            siiel,
            prB,
            pfB,
            syieldc,
            syieldt,
            syield,
            etapl,
            YNY,
            YNY5,
            YERRNOD,
            DSY,
            YNPL,
            DDD,
            dt,
            iplast,
            sp
        )
        # verification, from madcph.m, line 1232ff
        ETA5_ver=copy(ETA0)
        YNY5_ver = zeros(Ny, Nx)
        DSY_ver = zeros(Ny, Nx)
        ynpl=0
        ddd=0
        for i=1:1:Ny
            for j=1:1:Nx
                # Compute second stress invariant
                SIIB_ver=(SXY[i,j]^2+((SXX[i,j]+SXX[i+1,j]+SXX[i,j+1]+SXX[i+1,j+1])/4)^2)^0.5
                # Compute second invariant for a purely elastic stress build-up
                siiel_ver=SIIB_ver*(GGG[i,j]*dt+ETA[i,j])/ETA[i,j]
                # Compute total & fluid pressure
                prB_ver=(pr[i,j]+pr[i+1,j]+pr[i,j+1]+pr[i+1,j+1])/4
                pfB_ver=(pf[i,j]+pf[i+1,j]+pf[i,j+1]+pf[i+1,j+1])/4
                # Compute yielding stress
                syieldc_ver=COH[i,j]+FRI[i,j]*(prB_ver-pfB_ver); # Confined fracture
                syieldt_ver=TEN[i,j]+(prB_ver-pfB_ver); # Tensile fracture
                syield_ver=max(syieldt_ver,syieldc_ver,0); # Non-negative strength requirement
                # Update error for old yielding nodes
                ynn=0
                if(YNY[i,j]>0)
                    ynn=1
                    DSY_ver[i,j]=SIIB_ver-syield_ver
                    ddd=ddd+DSY_ver[i,j]^2
                    ynpl=ynpl+1
                end
                # Correcting viscosity for yielding
                if syield_ver<siiel_ver
                    # New viscosity for the basic node
                    etapl_ver=dt*GGG[i,j]*syield_ver/(siiel_ver-syield_ver)
                    if etapl_ver<ETA0[i,j]
                        # Recompute nodal visocity
                        ETA5_ver[i,j]=etapl_ver^(1-etawt)*ETA[i,j]^etawt
                        # Mark yielding nodes
                        YNY5_ver[i,j]=1
                        # Apply viscosity cutoff values
                        if ETA5_ver[i,j]<etamin
                            ETA5_ver[i,j]=etamin
                        elseif ETA5_ver[i,j]>etamax
                            ETA5_ver[i,j]=etamax
                        end
                        # Update Error for new yielding nodes
                        if ynn==0
                            DSY_ver[i,j]=SIIB_ver-syield_ver
                            ddd=ddd+DSY_ver[i,j]^2
                            ynpl=ynpl+1
                        end
                    else
                        ETA5_ver[i,j]=ETA0[i,j]
                        YNY5_ver[i,j]=0
                    end
                else
                    ETA5_ver[i,j]=ETA0[i,j]
                    YNY5_ver[i,j]=0
                end
            end
        end
        # Compute yielding error for markers
        if(ynpl>0)
            YERRNOD_ver[iplast]=(ddd/ynpl)^0.5
        end
        # test
        for j=1:1:Nx, i=1:1:Ny
            @test ETA5[i, j] ≈ ETA5_ver[i, j] rtol=1e-6
            @test YNY5[i, j] == YNY5_ver[i, j]
        end
        @test YERRNOD[iplast] ≈ YERRNOD_ver[iplast] rtol=1e-6
        @test complete == (ynpl==0 || iplast==nplast || YERRNOD[iplast]<yerrmax)
    end # testset "compute_nodal_adjustment!()"

    @testset "finalize_plastic_iteration_pass!()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        dt = sp.dtelastic
        dtkoef = sp.dtkoef
        dtstep = sp.dtstep
        iplast = 1
        # simulate data
        ETA = zeros(Ny, Nx)
        ETA5 = rand(Ny, Nx)
        ETA00 = rand(Ny, Nx)
        YNY = zeros(Bool, Ny, Nx)
        YNY5 = rand(Bool, Ny, Nx)
        YNY00 = rand(Bool, Ny, Nx)
        YNY_inv_ETA = zeros(Ny, Nx)
        ETA_ver = zeros(Ny, Nx)
        YNY_ver = zeros(Bool, Ny, Nx)
        YNY_inv_ETA_ver = zeros(Ny, Nx)
        dt_ver = sp.dtelastic
        # finalize_plastic_iteration_pass
        dt = HydrologyPlanetesimals.finalize_plastic_iteration_pass!(
            ETA,
            ETA5,
            ETA00,
            YNY,
            YNY5,
            YNY00,
            YNY_inv_ETA,
            dt,
            dtkoef,
            dtstep,
            iplast
        )
        # verification, from madcph.m, line 1301ff
        if trunc(Int, iplast/dtstep)*dtstep==iplast
            # Decrease timestep
            dt_ver=dt_ver/dtkoef
            # Reset old viscoplastic viscosity
            ETA_ver=ETA00
            YNY_ver=YNY00
        else
            # Use new viscoplastic viscosity
            ETA_ver=ETA5
            YNY_ver=YNY5
        end
        YNY_inv_ETA_ver.=YNY_ver./ETA_ver
        # test
        @test dt == dt_ver
        @test ETA == ETA_ver
        @test YNY == YNY_ver
        @test YNY_inv_ETA == YNY_inv_ETA_ver
    end # testset "finalize_plastic_iteration_pass!()"

    @testset "apply_subgrid_stress_diffusion!()" begin
        sp = HydrologyPlanetesimals.StaticParameters(
            Nxmc=1, Nymc=1, dsubgrids=0.5, dtelastic=0.5)
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        marknum = sp.start_marknum
        x, y = sp.x, sp.y
        xp, yp = sp.xp, sp.yp
        dx, dy = sp.dx, sp.dy
        dtm = sp.dtelastic
        dsubgrids = sp.dsubgrids
        etam = @SVector ones(3)
        # simulate markers
        xm = rand(-dx:0.1:x[end]+dx, marknum)
        ym = rand(-dy:0.1:y[end]+dy, marknum)
        tm = rand(1:3, marknum)
        gggm = rand(3)
        inv_gggtotalm = inv.([gggm[tm[m]] for m in 1:marknum])
        sxxm = rand(marknum)
        sxym = rand(marknum)
        SXX0 = rand(Ny1, Nx1)
        SXY0 = rand(Ny, Nx)
        DSXX = rand(Ny1, Nx1)
        DSXY = rand(Ny, Nx)
        SXXSUM = rand(Ny1, Nx1, Base.Threads.nthreads())
        SXYSUM = rand(Ny, Nx, Base.Threads.nthreads())
        WTPSUM = rand(Ny1, Nx1, Base.Threads.nthreads())
        WTSUM = rand(Ny, Nx, Base.Threads.nthreads())
        sxxm_ver = copy(sxxm)
        sxym_ver = copy(sxym)
        DSXX_ver = copy(DSXX)
        DSXY_ver = copy(DSXY)
        # apply subgrid stress diffusion
        HydrologyPlanetesimals.apply_subgrid_stress_diffusion!(
            xm,
            ym,
            tm,
            inv_gggtotalm,
            sxxm,
            sxym,
            SXX0,
            SXY0,
            DSXX,
            DSXY,
            SXXSUM,
            SXYSUM,
            WTPSUM,
            WTSUM,
            dtm,
            marknum,
            sp
        )
        # verification, from madcph.m, line 1374ff
        # Apply subgrid stress diffusion to markers
        if(dsubgrids>0)
        SXYSUM_ver = zeros(Ny, Nx)
        WTSUM_ver = zeros(Ny, Nx)
        SXXSUM_ver = zeros(Ny1, Nx1)
        WTPSUM_ver = zeros(Ny1, Nx1)
        for m=1:1:marknum
            # SIGMA'xx
            # Define i;j indexes for the upper left node
            j=trunc(Int, (xm[m]-xp[1])/dx)+1
            i=trunc(Int, (ym[m]-yp[1])/dy)+1
            if j<1
                j=1
            elseif j>Nx
                j=Nx
            end
            if i<1
                i=1
            elseif i>Ny
                i=Ny
            end
            # Compute distances
            dxmj=xm[m]-xp[j]
            dymi=ym[m]-yp[i]
            # Compute weights
            wtmij=(1-dxmj/dx)*(1-dymi/dy)
            wtmi1j=(1-dxmj/dx)*(dymi/dy);    
            wtmij1=(dxmj/dx)*(1-dymi/dy)
            wtmi1j1=(dxmj/dx)*(dymi/dy)
            # Compute marker-node SIGMA'xx difference
            dsxxm0=sxxm_ver[m]-(SXX0[i,j]*wtmij+SXX0[i+1,j]*wtmi1j+ SXX0[i,j+1]*wtmij1+SXX0[i+1,j+1]*wtmi1j1)
            # Relax stress difference
            dsxxm1=dsxxm0*exp(-dsubgrids*dtm/(etam[tm[m]]/gggm[tm[m]]))
            # Correct marker stress
            ddsxxm_ver=dsxxm1-dsxxm0
            sxxm_ver[m]=sxxm_ver[m]+ddsxxm_ver
            # Update subgrid diffusion on nodes
            # i;j Node
            SXXSUM_ver[i,j]=SXXSUM_ver[i,j]+ddsxxm_ver*wtmij
            WTPSUM_ver[i,j]=WTPSUM_ver[i,j]+wtmij
            # i+1;j Node
            SXXSUM_ver[i+1,j]=SXXSUM_ver[i+1,j]+ddsxxm_ver*wtmi1j
            WTPSUM_ver[i+1,j]=WTPSUM_ver[i+1,j]+wtmi1j
            # i;j+1 Node
            SXXSUM_ver[i,j+1]=SXXSUM_ver[i,j+1]+ddsxxm_ver*wtmij1
            WTPSUM_ver[i,j+1]=WTPSUM_ver[i,j+1]+wtmij1
            # i+1;j+1 Node
            SXXSUM_ver[i+1,j+1]=SXXSUM_ver[i+1,j+1]+ddsxxm_ver*wtmi1j1
            WTPSUM_ver[i+1,j+1]=WTPSUM_ver[i+1,j+1]+wtmi1j1
            # SIGMAxy
            # Define i;j indexes for the upper left node
            j=trunc(Int, (xm[m]-x[1])/dx)+1
            i=trunc(Int, (ym[m]-y[1])/dy)+1
            if j<1
                j=1
            elseif j>Nx-1
                j=Nx-1
            end
            if i<1
                i=1
            elseif i>Ny-1
                i=Ny-1
            end
            # Compute distances
            dxmj=xm[m]-x[j]
            dymi=ym[m]-y[i]
            # Compute weights
            wtmij=(1-dxmj/dx)*(1-dymi/dy)
            wtmi1j=(1-dxmj/dx)*(dymi/dy);    
            wtmij1=(dxmj/dx)*(1-dymi/dy)
            wtmi1j1=(dxmj/dx)*(dymi/dy)
            # Compute marker-node SIGMAxy difference
            dsxym0=sxym_ver[m]-(SXY0[i,j]*wtmij+SXY0[i+1,j]*wtmi1j+  SXY0[i,j+1]*wtmij1+SXY0[i+1,j+1]*wtmi1j1)
            # Relax stress difference
            dsxym1=dsxym0*exp(-dsubgrids*dtm/(etam[tm[m]]/gggm[tm[m]]))
            # Correct marker stress
            ddsxym_ver=dsxym1-dsxym0
            sxym_ver[m]=sxym_ver[m]+ddsxym_ver
            # Update subgrid diffusion on nodes
            # i;j Node
            SXYSUM_ver[i,j]=SXYSUM_ver[i,j]+ddsxym_ver*wtmij
            WTSUM_ver[i,j]=WTSUM_ver[i,j]+wtmij
            # i+1;j Node
            SXYSUM_ver[i+1,j]=SXYSUM_ver[i+1,j]+ddsxym_ver*wtmi1j
            WTSUM_ver[i+1,j]=WTSUM_ver[i+1,j]+wtmi1j
            # i;j+1 Node
            SXYSUM_ver[i,j+1]=SXYSUM_ver[i,j+1]+ddsxym_ver*wtmij1
            WTSUM_ver[i,j+1]=WTSUM_ver[i,j+1]+wtmij1
            # i+1;j+1 Node
            SXYSUM_ver[i+1,j+1]=SXYSUM_ver[i+1,j+1]+ddsxym_ver*wtmi1j1
            WTSUM_ver[i+1,j+1]=WTSUM_ver[i+1,j+1]+wtmi1j1
        end
        # Compute DSXXsubgrid_ver
        DSXXsubgrid_ver = zeros(Ny1, Nx1)
        # P-nodes
        for j=2:1:Nx
            for i=2:1:Ny
                if(WTPSUM_ver[i,j]>0)
                    DSXXsubgrid_ver[i,j]=SXXSUM_ver[i,j]/WTPSUM_ver[i,j]
                end
            end
        end
        # Correct DSXX_ver
        DSXX_ver=DSXX_ver-DSXXsubgrid_ver
        # Compute DSXYsubgrid_ver
        DSXYsubgrid_ver = zeros(Ny, Nx)
        # Basic nodes
        for j=1:1:Nx
            for i=1:1:Ny
                if(WTSUM_ver[i,j]>0)
                    DSXYsubgrid_ver[i,j]=SXYSUM_ver[i,j]/WTSUM_ver[i,j]
                end
            end
        end
        # Correct DSXY_ver
        DSXY_ver=DSXY_ver-DSXYsubgrid_ver
        end
    # test
    @test sxxm ≈ sxxm_ver atol=1e-6
    @test sxym ≈ sxym_ver atol=1e-6
    @test DSXX ≈ DSXX_ver atol=1e-6
    @test DSXY ≈ DSXY_ver atol=1e-6
    end # testset "apply_subgrid_stress_diffusion!()"

    @testset "update_marker_stress!()" begin
        sp = HydrologyPlanetesimals.StaticParameters(Nxmc=1, Nymc=1)
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        marknum = sp.start_marknum
        x, y = sp.x, sp.y
        xp, yp = sp.xp, sp.yp
        dx, dy = sp.dx, sp.dy
        # simulate markers
        xm = rand(-dx:0.1:x[end]+dx, marknum)
        ym = rand(-dy:0.1:y[end]+dy, marknum)
        sxxm = rand(marknum)
        sxym = rand(marknum)
        DSXX = rand(Ny1, Nx1)
        DSXY = rand(Ny, Nx)
        sxxm_ver = copy(sxxm)
        sxym_ver = copy(sxym)
        # update marker stress
        HydrologyPlanetesimals.update_marker_stress!(
            xm, ym, sxxm, sxym, DSXX, DSXY, marknum, sp)
        # verification, from madcph.m, line 1495ff
        for m=1:1:marknum
            # SIGMA'xx
            # Define i;j indexes for the upper left node
            j=trunc(Int, (xm[m]-xp[1])/dx)+1
            i=trunc(Int, (ym[m]-yp[1])/dy)+1
            if j<2
                j=2
            elseif j>Nx-1
                j=Nx-1
            end
            if i<2
                i=2
            elseif i>Ny-1
                i=Ny-1
            end
            # Compute distances
            dxmj=xm[m]-xp[j]
            dymi=ym[m]-yp[i]
            # Compute weights
            wtmij=(1-dxmj/dx)*(1-dymi/dy)
            wtmi1j=(1-dxmj/dx)*(dymi/dy);    
            wtmij1=(dxmj/dx)*(1-dymi/dy)
            wtmi1j1=(dxmj/dx)*(dymi/dy)
            # Update marker by SIGMA'xx change 
            sxxm_ver[m]=sxxm_ver[m]+(DSXX[i,j]*wtmij+DSXX[i+1,j]*wtmi1j+ DSXX[i,j+1]*wtmij1+DSXX[i+1,j+1]*wtmi1j1)
        
            # SIGMAxy
            # Define i;j indexes for the upper left node
            j=trunc(Int, (xm[m]-x[1])/dx)+1
            i=trunc(Int, (ym[m]-y[1])/dy)+1
            if j<1 
                j=1
            elseif j>Nx-1 
                j=Nx-1
            end
            if i<1
                i=1
            elseif i>Ny-1
                i=Ny-1
            end
            # Compute distances
            dxmj=xm[m]-x[j]
            dymi=ym[m]-y[i]
            # Compute weights
            wtmij=(1-dxmj/dx)*(1-dymi/dy)
            wtmi1j=(1-dxmj/dx)*(dymi/dy);    
            wtmij1=(dxmj/dx)*(1-dymi/dy)
            wtmi1j1=(dxmj/dx)*(dymi/dy)
            # Update marker by SIGMA'xx change 
            sxym_ver[m]=sxym_ver[m]+(DSXY[i,j]*wtmij+DSXY[i+1,j]*wtmi1j+ DSXY[i,j+1]*wtmij1+DSXY[i+1,j+1]*wtmi1j1)
        end
        # test
        @test sxxm ≈ sxxm_ver atol=1e-6
        @test sxym ≈ sxym_ver atol=1e-6
    end # testset "update_marker_stress!()"

    @testset "compute_shear_heating!()" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        HS = zeros(Ny1, Nx1)
        SXYEXY = zeros(Ny, Nx)
        ETA = rand(Ny, Nx)
        SXY = rand(Ny, Nx)
        ETAP = rand(Ny1, Nx1)
        SXX = rand(Ny1, Nx1)
        RX = rand(Ny1, Nx1)
        RY = rand(Ny1, Nx1)
        qxD = rand(Ny1, Nx1)
        qyD = rand(Ny1, Nx1)
        PHI = rand(Ny1, Nx1)
        ETAPHI = rand(Ny1, Nx1)
        pr = rand(Ny1, Nx1)
        pf = rand(Ny1, Nx1)
        # compute shear compute shear heating
        HydrologyPlanetesimals.compute_shear_heating!(
            HS,
            SXYEXY,
            ETA,
            SXY,
            ETAP,
            SXX,
            RX,
            RY,
            qxD,
            qyD,
            PHI,
            ETAPHI,
            pr,
            pf,
            sp
        )
        # verification, from madcph.m, line 1551ff
        HS_ver = zeros(Ny1, Nx1); # Adiabatic heating, W/m^3
        for j=2:1:Nx
            for i=2:1:Ny
                # Average SXY*EXY
                SXYEXY_ver=(SXY[i,j]^2/ETA[i,j]+SXY[i-1,j]^2/ETA[i-1,j]+ SXY[i,j-1]^2/ETA[i,j-1]+SXY[i-1,j-1]^2/ETA[i-1,j-1])/4
                # HS
                HS_ver[i,j]=SXX[i,j]^2/ETAP[i,j]+SXYEXY_ver+ (pr[i,j]-pf[i,j])^2/(1-PHI[i,j])/ETAPHI[i,j]+ (RX[i,j-1]*qxD[i,j-1]^2+RX[i,j]*qxD[i,j]^2)/2+ (RY[i-1,j]*qyD[i-1,j]^2+RY[i,j]*qyD[i,j]^2)/2
            end
        end
        # test
        @test HS ≈ HS_ver atol=1e-6
    end # testset "compute_shear_heating!()"

    @testset "compute_adiabatic_heating!()" begin
        sp = HydrologyPlanetesimals.StaticParameters(Nx=20, Ny=20)
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        dx, dy = sp.dx, sp.dy
        HA = zeros(Ny1, Nx1)
        tk1 = rand(Ny1, Nx1)
        ALPHA = rand(Ny1, Nx1)
        ALPHAF = rand(Ny1, Nx1)
        PHI = rand(Ny1, Nx1)
        vx = rand(Ny1, Nx1)
        vy = rand(Ny1, Nx1)
        vxf = rand(Ny1, Nx1)
        vyf = rand(Ny1, Nx1)
        ps = rand(Ny1, Nx1)
        pf = rand(Ny1, Nx1)
        HA_ver = zero(HA)
        # compute adiabatic heating
        HydrologyPlanetesimals.compute_adiabatic_heating!(
            HA, tk1, ALPHA, ALPHAF, PHI, vx, vy, vxf, vyf, ps, pf, sp)
        # verification, from madcph.m, line 1573ff
        for j=2:1:Nx
            for i=2:1:Ny
                # HA
                # Indirect calculation of dpdt
                # Average vy; vx; vxf; vyf
                VXP=(vx[i,j]+vx[i,j-1])/2
                VYP=(vy[i,j]+vy[i-1,j])/2
                VXFP=(vxf[i,j]+vxf[i,j-1])/2
                VYFP=(vyf[i,j]+vyf[i-1,j])/2
                # Evaluate DPsolid/Dt with upwind differences
                if VXP<0
                    dpsdx=(ps[i,j]-ps[i,j-1])/dx
                else
                    dpsdx=(ps[i,j+1]-ps[i,j])/dx
                end
                if VYP<0
                    dpsdy=(ps[i,j]-ps[i-1,j])/dy
                else
                    dpsdy=(ps[i+1,j]-ps[i,j])/dy
                end
                dpsdt=VXP*dpsdx+VYP*dpsdy
                # Evaluate DPfluid/Dt with upwind differences
                if VXFP>0
                    dpfdx=(pf[i,j]-pf[i,j-1])/dx
                else
                    dpfdx=(pf[i,j+1]-pf[i,j])/dx
                end
                if VYFP>0
                    dpfdy=(pf[i,j]-pf[i-1,j])/dy
                else
                    dpfdy=(pf[i+1,j]-pf[i,j])/dy
                end
                dpfdt=VXFP*dpsdx+VYFP*dpsdy
        #         # Direct calculation of dpdt
        #         dpsdt=(ps[i,j]-ps0[i,j])/dt
        #         dpfdt=(pf[i,j]-pf0[i,j])/dt
                # HA
                HA_ver[i,j]=(1-PHI[i,j])*tk1[i,j]*ALPHA[i,j]*dpsdt+ PHI[i,j]*tk1[i,j]*ALPHAF[i,j]*dpfdt
            end
        end
        # test
        @test HA ≈ HA_ver atol=1e-6
    end # testset "compute_adiabatic_heating!()"

    @testset "perform_thermal_iterations!()" begin
        sp = HydrologyPlanetesimals.StaticParameters(Nx=4, Ny=4)
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        dx, dy = sp.dx, sp.dy
        DTmax = sp.DTmax
        dtm = 0.0001
        tk0 = rand(Ny1, Nx1)
        tk1 = rand(Ny1, Nx1)
        tk2 = rand(Ny1, Nx1)
        RHOCP = rand(Ny1, Nx1)
        KX = rand(Ny1, Nx1)
        KY = rand(Ny1, Nx1)
        HR = rand(Ny1, Nx1)
        HA = rand(Ny1, Nx1)
        HS = rand(Ny1, Nx1)
        DT = zeros(Ny1, Nx1)
        DT0 = zeros(Ny1, Nx1)
        DT_ver = zeros(Ny1, Nx1)
        DT0_ver = zeros(Ny1, Nx1)
        tk0_ver = copy(tk0)
        tk1_ver = copy(tk1)
        tk2_ver = copy(tk2)
        # perform thermal iterations
        HydrologyPlanetesimals.perform_thermal_iterations!(
            tk0, tk1, tk2, DT, DT0, RHOCP, KX, KY, HR, HA, HS, dtm, sp)
        # verification, from madcph.m, line 1618ff
        LT = zeros(Ny1*Nx1, Ny1*Nx1)
        RT = zeros(Ny1*Nx1)
        ST = zeros(Ny1*Nx1)
        tk0_ver=tk1_ver
        dtt=dtm
        dttsum=0
        titer=1
        while dttsum<dtm
        # Composing global matrixes LT[], RT[]
        # Going through all points of the 2D grid &
        # composing respective equations
        for j=1:1:Nx1
            for i=1:1:Ny1
                # Define global index in algebraic space
                gk=(j-1)*Ny1+i
                # External points
                if i==1 || i==Ny1 || j==1 || j==Nx1
                    # Boundary Condition
                    # Top BC: T=273
                    if i==1 && j>1 && j<Nx1
                        LT[gk,gk]=1; # Left part
                        LT[gk,gk+1]=-1; # Left part
                        RT[gk]=0; # Right part
                    end
                    # Bottom BC: T=1500
                    if i==Ny1 && j>1 && j<Nx1
                        LT[gk,gk]=1; # Left part
                        LT[gk,gk-1]=-1; # Left part
                        RT[gk]=0; # Right part
                    end
                    # Left BC: dT/dx=0
                    if j==1
                        LT[gk,gk]=1; # Left part
                        LT[gk,gk+Ny1]=-1; # Left part
                        RT[gk]=0; # Right part
                    end
                    # Right BC: dT/dx=0
                    if j==Nx1
                        LT[gk,gk]=1; # Left part
                        LT[gk,gk-Ny1]=-1; # Left part
                        RT[gk]=0; # Right part
                    end
                else
                # Internal points: Temperature eq.
                # RHO*CP*dT/dt=-dqx/dx-dqy/dy+Hr+Hs+Ha
                #          Tdt2
                #           |
                #          Ky1
                #           |
                #Tdt1-Kx1-T03;Tdt3-Kx2-Tdt5
                #           |
                #          Ky2
                #           |
                #          Tdt4
                #
                # Left part
                Kx1=KX[i,j-1]; 
                Kx2=KX[i,j]; 
                Ky1=KY[i-1,j]; 
                Ky2=KY[i,j]; 
                LT[gk,gk-Ny1]=-Kx1/dx^2; # T1
                LT[gk,gk-1]=-Ky1/dy^2; # FI2
                LT[gk,gk]=RHOCP[i,j]/dtt+(Kx1+Kx2)/dx^2+(Ky1+Ky2)/dy^2; # FI3
                LT[gk,gk+1]=-Ky2/dy^2; # FI4
                LT[gk,gk+Ny1]=-Kx2/dx^2; # FI5
                # Right part
                RT[gk]=RHOCP[i,j]/dtt*tk1[i,j]+HR[i,j]+HA[i,j]+HS[i,j]
                end
            end
        end
        # Solving matrixes
        ST=LT\RT; # Obtaining algebraic vector of solutions ST[]
        # Reload solutions ST[] to geometrical array Tdt[]
        # Going through all grid points
        for j=1:1:Nx1
            for i=1:1:Ny1
                # Compute global index
                gk=(j-1)*Ny1+i
                # Reload solution
                tk2_ver[i,j]=ST[gk]
            end
        end
        # Compute DT
        DT_ver=tk2_ver-tk1_ver
        titer
        dtt
        if titer==1
            # Apply thermal timestepping condition
            maxDTcurrent=maximum(abs, DT_ver)
            if maxDTcurrent>DTmax 
                dtt=dtt/maxDTcurrent*DTmax
            else
                dttsum=dttsum+dtt; # Update dttsum
            end
        else
            dttsum=dttsum+dtt; # Update dttsum
            # Adjust timestep
            if dtt>dtm-dttsum
                dtt=dtm-dttsum
            end
        end
        titer=titer+1; # Update iteration counter
        end
        # Compute/save overall temperature changes
        DT_ver=tk2_ver-tk0_ver
        DT0_ver=DT_ver
        # test
        @test DT ≈ DT_ver atol=1e-6
        @test DT0 ≈ DT0_ver atol=1e-6
    end # testset "perform_thermal_iterations!()"

    @testset "apply_subgrid_temperature_diffusion!()" begin
        sp = HydrologyPlanetesimals.StaticParameters(Nx=5,Ny=5,
            Nxmc=1, Nymc=1, dsubgridt=0.5)
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        x, y = sp.x, sp.y
        xp, yp = sp.xp, sp.yp
        dx, dy = sp.dx, sp.dy
        rhocpsolidm, rhocpfluidm = sp.rhocpsolidm, sp.rhocpfluidm
        ksolidm, kfluidm = sp.ksolidm, sp.kfluidm
        marknum = sp.start_marknum
        dtm = sp.dtelastic
        dsubgridt = sp.dsubgridt
        # simulate markers
        xm = rand(-dx:0.1:x[end]+dx, marknum)
        ym = rand(-dy:0.1:y[end]+dy, marknum)
        tm = rand(1:3, marknum)
        tkm = rand(marknum)
        phim = rand(marknum)
        tk1 = rand(Ny1, Nx1)
        DT = rand(Ny1, Nx1)
        TKSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
        RHOCPSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
        tkm_ver = copy(tkm)
        DT_ver = copy(DT)
        # apply subgrid stress diffusion
        HydrologyPlanetesimals.apply_subgrid_temperature_diffusion!(
            xm, ym, tm, tkm, phim, tk1, DT, TKSUM, RHOCPSUM, dtm, marknum, sp)
        # verification, from madcph.m, line 1731ff
        # Apply subgrid stress diffusion to markers
        if dsubgridt>0
            TKSUM_ver = zeros(Ny1, Nx1)
            RHOCPSUM_ver = zeros(Ny1, Nx1)
            for m=1:1:marknum
                # Define i;j indexes for the upper left node
                j=trunc(Int, (xm[m]-xp[1])/dx)+1
                i=trunc(Int, (ym[m]-yp[1])/dy)+1
                if j<1
                    j=1
                elseif j>Nx
                    j=Nx
                end
                if i<1
                    i=1
                elseif i>Ny
                    i=Ny
                end
                # Compute distances
                dxmj=xm[m]-xp[j]
                dymi=ym[m]-yp[i]
                # Compute weights
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy);    
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Compute marker-node T difference
                dtkm0=tkm_ver[m]-(tk1[i,j]*wtmij+tk1[i+1,j]*wtmi1j+ tk1[i,j+1]*wtmij1+tk1[i+1,j+1]*wtmi1j1)
                # Compute marker parameters
                if tm[m]<3
                    # Rocks
                    rhocptotalm=rhocpsolidm[tm[m]]*(1-phim[m])+rhocpfluidm[tm[m]]*phim[m]
                    ktotalm=(ksolidm[tm[m]]*kfluidm[tm[m]]/2+((ksolidm[tm[m]]*(3*phim[m]-2)+ kfluidm[tm[m]]*(1-3*phim[m]))^2)/16)^0.5-(ksolidm[tm[m]]*(3*phim[m]-2)+ kfluidm[tm[m]]*(1-3*phim[m]))/4
                else
                    # Sticky air
                    rhocptotalm=rhocpsolidm[tm[m]]
                    ktotalm=ksolidm[tm[m]]
                end    # Relax temperature difference
                dtkm1=dtkm0*exp(-dsubgridt*ktotalm*dtm/rhocptotalm*(2/dx^2+2/dy^2))
                # Correct marker temperature
                ddtkm=dtkm1-dtkm0
                # @info "ver" m dtkm0 ddtkm
                tkm_ver[m]=tkm_ver[m]+ddtkm
                # Update subgrid diffusion on nodes
                # i;j Node
                TKSUM_ver[i,j]=TKSUM_ver[i,j]+ddtkm*rhocptotalm*wtmij
                RHOCPSUM_ver[i,j]=RHOCPSUM_ver[i,j]+rhocptotalm*wtmij
                # i+1;j Node
                TKSUM_ver[i+1,j]=TKSUM_ver[i+1,j]+ddtkm*rhocptotalm*wtmi1j
                RHOCPSUM_ver[i+1,j]=RHOCPSUM_ver[i+1,j]+rhocptotalm*wtmi1j
                # i;j+1 Node
                TKSUM_ver[i,j+1]=TKSUM_ver[i,j+1]+ddtkm*rhocptotalm*wtmij1
                RHOCPSUM_ver[i,j+1]=RHOCPSUM_ver[i,j+1]+rhocptotalm*wtmij1
                # i+1;j+1 Node
                TKSUM_ver[i+1,j+1]=TKSUM_ver[i+1,j+1]+ddtkm*rhocptotalm*wtmi1j1
                RHOCPSUM_ver[i+1,j+1]=RHOCPSUM_ver[i+1,j+1]+rhocptotalm*wtmi1j1
            end
            # Compute DTsubgrid
            DTsubgrid = zeros(Ny1, Nx1)
            # P-nodes
            for j=1:1:Nx1
                for i=1:1:Ny1
                    if RHOCPSUM_ver[i,j]>0
                        DTsubgrid[i,j]=TKSUM_ver[i,j]/RHOCPSUM_ver[i,j]
                    end
                end
            end
            # Correct DT
            DT_ver=DT_ver-DTsubgrid
        end
        # test
        @test tkm ≈ tkm_ver atol=1e-6
        @test DT ≈ DT_ver atol=1e-6
    end # testset "apply_subgrid_temperature_diffusion!()"
end


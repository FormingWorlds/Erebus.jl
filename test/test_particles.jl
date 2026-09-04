@testset "Particles" begin
    @testset "setup_interpolated_properties()" begin
        # setup interpolated grid properties
        (ETA0SUM, ETASUM, GGGSUM, SXYSUM, COHSUM, TENSUM, FRISUM, WTSUM, RHOXSUM, RHOFXSUM, KXSUM, PHIXSUM, RXSUM, WTXSUM, RHOYSUM, RHOFYSUM, KYSUM, PHIYSUM, RYSUM, WTYSUM, RHOSUM, RHOCPSUM, ALPHASUM, ALPHAFSUM, HRSUM, GGGPSUM, SXXSUM, TKSUM, PHISUM, DMPSUM, DHPSUM, XWSSUM, WTPSUM) = Erebus.setup_interpolated_properties()
        # verification and test
        @test ETA0SUM == zeros(Float64, Ny, Nx)
        @test ETASUM == zeros(Float64, Ny, Nx)
        @test GGGSUM == zeros(Float64, Ny, Nx)
        @test SXYSUM == zeros(Float64, Ny, Nx)
        @test COHSUM == zeros(Float64, Ny, Nx)
        @test TENSUM == zeros(Float64, Ny, Nx)
        @test FRISUM == zeros(Float64, Ny, Nx)
        @test WTSUM == zeros(Float64, Ny, Nx)
        @test RHOXSUM == zeros(Float64, Ny1, Nx1)
        @test RHOFXSUM == zeros(Float64, Ny1, Nx1)
        @test KXSUM == zeros(Float64, Ny1, Nx1)
        @test PHIXSUM == zeros(Float64, Ny1, Nx1)
        @test RXSUM == zeros(Float64, Ny1, Nx1)
        @test WTXSUM == zeros(Float64, Ny1, Nx1)
        @test RHOYSUM == zeros(Float64, Ny1, Nx1)
        @test RHOFYSUM == zeros(Float64, Ny1, Nx1)
        @test KYSUM == zeros(Float64, Ny1, Nx1)
        @test PHIYSUM == zeros(Float64, Ny1, Nx1)
        @test RYSUM == zeros(Float64, Ny1, Nx1)
        @test WTYSUM == zeros(Float64, Ny1, Nx1)
        @test RHOSUM == zeros(Float64, Ny1, Nx1)
        @test RHOCPSUM == zeros(Float64, Ny1, Nx1)
        @test ALPHASUM == zeros(Float64, Ny1, Nx1)
        @test ALPHAFSUM == zeros(Float64, Ny1, Nx1)
        @test HRSUM == zeros(Float64, Ny1, Nx1)
        @test GGGPSUM == zeros(Float64, Ny1, Nx1)
        @test SXXSUM == zeros(Float64, Ny1, Nx1)
        @test TKSUM == zeros(Float64, Ny1, Nx1)
        @test PHISUM == zeros(Float64, Ny1, Nx1)
        @test DMPSUM == zeros(Float64, Ny1, Nx1)
        @test DHPSUM == zeros(Float64, Ny1, Nx1)
        @test XWSSUM == zeros(Float64, Ny1, Nx1)
        @test WTPSUM == zeros(Float64, Ny1, Nx1)
    end # testset "setup_interpolated_properties()"

    @testset "reset_interpolated_properties!()" begin
        (ETA0SUM, ETASUM, GGGSUM, SXYSUM, COHSUM, TENSUM, FRISUM, WTSUM, RHOXSUM, RHOFXSUM, KXSUM, PHIXSUM, RXSUM, WTXSUM, RHOYSUM, RHOFYSUM, KYSUM, PHIYSUM, RYSUM, WTYSUM, RHOSUM, RHOCPSUM, ALPHASUM, ALPHAFSUM, HRSUM, GGGPSUM, SXXSUM, TKSUM, PHISUM, WTPSUM) = Erebus.setup_interpolated_properties()
        ETA0SUM[1, 1] = one(1.0)
        ETASUM[1, 1] = one(1.0)
        GGGSUM[1, 1] = one(1.0)
        SXYSUM[1, 1] = one(1.0)
        COHSUM[1, 1] = one(1.0)
        TENSUM[1, 1] = one(1.0)
        FRISUM[1, 1] = one(1.0)
        WTSUM[1, 1] = one(1.0)
        RHOXSUM[1, 1] = one(1.0)
        RHOFXSUM[1, 1] = one(1.0)
        KXSUM[1, 1] = one(1.0)
        PHIXSUM[1, 1] = one(1.0)
        RXSUM[1, 1] = one(1.0)
        WTXSUM[1, 1] = one(1.0)
        RHOYSUM[1, 1] = one(1.0)
        RHOFYSUM[1, 1] = one(1.0)
        KYSUM[1, 1] = one(1.0)
        PHIYSUM[1, 1] = one(1.0)
        RYSUM[1, 1] = one(1.0)
        WTYSUM[1, 1] = one(1.0)
        RHOSUM[1, 1] = one(1.0)
        RHOCPSUM[1, 1] = one(1.0)
        ALPHASUM[1, 1] = one(1.0)
        ALPHAFSUM[1, 1] = one(1.0)
        HRSUM[1, 1] = one(1.0)
        GGGPSUM[1, 1] = one(1.0)
        SXXSUM[1, 1] = one(1.0)
        TKSUM[1, 1] = one(1.0)
        PHISUM[1, 1] = one(1.0)
        WTPSUM[1, 1] = one(1.0)
        Erebus.reset_interpolated_properties!(
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
            WTPSUM,
        )
        @test ETA0SUM == zeros(Float64, Ny, Nx)
        @test ETASUM == zeros(Float64, Ny, Nx)
        @test GGGSUM == zeros(Float64, Ny, Nx)
        @test SXYSUM == zeros(Float64, Ny, Nx)
        @test COHSUM == zeros(Float64, Ny, Nx)
        @test TENSUM == zeros(Float64, Ny, Nx)
        @test FRISUM == zeros(Float64, Ny, Nx)
        @test WTSUM == zeros(Float64, Ny, Nx)
        @test RHOXSUM == zeros(Float64, Ny1, Nx1)
        @test RHOFXSUM == zeros(Float64, Ny1, Nx1)
        @test KXSUM == zeros(Float64, Ny1, Nx1)
        @test PHIXSUM == zeros(Float64, Ny1, Nx1)
        @test RXSUM == zeros(Float64, Ny1, Nx1)
        @test WTXSUM == zeros(Float64, Ny1, Nx1)
        @test RHOYSUM == zeros(Float64, Ny1, Nx1)
        @test RHOFYSUM == zeros(Float64, Ny1, Nx1)
        @test KYSUM == zeros(Float64, Ny1, Nx1)
        @test PHIYSUM == zeros(Float64, Ny1, Nx1)
        @test RYSUM == zeros(Float64, Ny1, Nx1)
        @test WTYSUM == zeros(Float64, Ny1, Nx1)
        @test RHOSUM == zeros(Float64, Ny1, Nx1)
        @test RHOCPSUM == zeros(Float64, Ny1, Nx1)
        @test ALPHASUM == zeros(Float64, Ny1, Nx1)
        @test ALPHAFSUM == zeros(Float64, Ny1, Nx1)
        @test HRSUM == zeros(Float64, Ny1, Nx1)
        @test GGGPSUM == zeros(Float64, Ny1, Nx1)
        @test SXXSUM == zeros(Float64, Ny1, Nx1)
        @test TKSUM == zeros(Float64, Ny1, Nx1)
        @test PHISUM == zeros(Float64, Ny1, Nx1)
        @test WTPSUM == zeros(Float64, Ny1, Nx1)
    end # testset "reset_interpolated_properties!()"

    @testset "reset_interpolated_properties!()" begin
        DMPSUM = rand(rgen, Ny1, Nx1)
        DHPSUM = rand(rgen, Ny1, Nx1)
        WTPSUM = rand(rgen, Ny1, Nx1)
        Erebus.reset_thermochemical_properties!(DMPSUM, DHPSUM, WTPSUM)
        @test DMPSUM == zeros(Float64, Ny1, Nx1)
        @test DHPSUM == zeros(Float64, Ny1, Nx1)
        @test WTPSUM == zeros(Float64, Ny1, Nx1)
    end # testset "reset_thermochemical_properties!()"

    @testset "setup_marker_properties()" begin
        marknum = start_marknum
        # setup marker properties
        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0) = Erebus.setup_marker_properties(
            marknum
        )
        # verification, from HTM-planetary.m, line 115ff
        Nxmc_ver = 4 # Number of markers per cell in horizontal direction
        Nymc_ver = 4 # Number of markers per cell in vertical direction
        Nxm_ver = (Nx-1)*Nxmc # Marker grid resolution in horizontal direction
        Nym_ver = (Ny-1)*Nymc # Marker grid resolution in vertical direction
        dxm_ver = xsize/Nxm # Marker grid step in horizontal direction,m
        dym_ver = ysize/Nym # Marker grid step in vertical direction,m
        marknum_ver = Nxm*Nym # Number of markers
        xm_ver = zeros(marknum) # Horizontal coordinates, m
        ym_ver = zeros(marknum) # Vertical coordinates, m
        tm_ver = zeros(Int, marknum) # Material type
        tkm_ver = zeros(marknum) # Marker temperature, K
        sxxm_ver = zeros(marknum) # SIGMA'xx, Pa
        sxym_ver = zeros(marknum) # SIGMAxy, Pa
        etavpm_ver = zeros(marknum) # Visco-plastic viscosity, Pa
        phim_ver = zeros(marknum) # Marker porosity
        phinewm_ver = zeros(marknum) # reacted Marker porosity
        pfm0_ver = zeros(marknum) # previous marker fluid pressure
        XWsolidm_ver = zeros(marknum) # marker melt molar fraction
        XWsolidm0_ver = zeros(marknum) # previous marker melt molar fraction
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
        @test phinewm == phinewm_ver
        @test pfm0 == pfm0_ver
        @test XWsolidm == XWsolidm_ver
        @test XWsolidm0 == XWsolidm0_ver
    end # testset "setup_marker_properties()"

    @testset "setup_marker_properties_helpers()" begin
        marknum = start_marknum
        # setup marker properties
        (rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm, tkm_rhocptotalm, etafluidcur_inv_kphim, inv_gggtotalm, fricttotalm, cohestotalm, tenstotalm, rhofluidcur, alphasolidcur, alphafluidcur) = Erebus.setup_marker_properties_helpers(
            marknum
        )
        # test
        @test rhocptotalm == zeros(Float64, marknum)
        @test rhocptotalm == zeros(Float64, marknum)
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

    @testset "setup_marker_geometry_helpers()" begin
        # setup marker geometry helpers
        (
            mdis,
            mnum,
            # mtyp,
            # mpor
        ) = Erebus.setup_marker_geometry_helpers()
        # test
        @test typeof(mdis) == Matrix{Float64}
        @test size(mdis) == (Nym, Nxm)
        @test rand(rgen, mdis) == mdis_init
        @test typeof(mnum) == Matrix{Int}
        @test size(mnum) == (Nym, Nxm)
    end # testset " setup_marker_geometry_helpers()"

    @testset "define_markers!() & compute_marker_properties!()" begin
        marknum = start_marknum
        mode = marker_property_mode
        hrsolidm, hrfluidm = start_hrsolidm, start_hrfluidm
        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0) = Erebus.setup_marker_properties(
            marknum
        )
        (rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm, tkm_rhocptotalm, etafluidcur_inv_kphim, inv_gggtotalm, fricttotalm, cohestotalm, tenstotalm, rhofluidcur, alphasolidcur, alphafluidcur) = Erebus.setup_marker_properties_helpers(
            marknum
        )
        xm_ver = zeros(marknum)
        ym_ver = zeros(marknum)
        tm_ver = zeros(Int, marknum)
        tkm_ver = zeros(marknum)
        phim_ver = zeros(marknum)
        etavpm_ver = zeros(marknum)
        rhototalm_ver = zeros(marknum)
        rhocptotalm_ver = zeros(marknum)
        hrtotalm_ver = zeros(marknum)
        ktotalm_ver = zeros(marknum)
        gggtotalm_ver = zeros(marknum)
        fricttotalm_ver = zeros(marknum)
        cohestotalm_ver = zeros(marknum)
        tenstotalm_ver = zeros(marknum)
        etasolidcur_ver = zeros(marknum)
        etafluidcur_ver = zeros(marknum)
        kphim_ver = zeros(marknum)
        rhofluidcur_ver = zeros(marknum)
        etatotalm_ver = zeros(marknum)
        XWsolidm0_ver = zeros(marknum)
        MH2O = MH₂O
        VDsolid = VDˢ
        VWsolid = VWˢ
        rhoH2Ofluid = ρH₂Oᶠ
        rhoH2Ofluidice = ρH₂Oᶠⁱ
        # define markers
        Erebus.define_markers!(
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
            tkm,
            inv_gggtotalm,
            fricttotalm,
            cohestotalm,
            tenstotalm,
            rhofluidcur,
            alphasolidcur,
            alphafluidcur,
            XWsolidm0,
            randomized=false,
        )
        for m in 1:1:marknum
            Erebus.compute_marker_properties!(
                m,
                tm,
                tkm,
                rhototalm,
                rhocptotalm,
                etatotalm,
                hrtotalm,
                ktotalm,
                tkm_rhocptotalm,
                etafluidcur_inv_kphim,
                hrsolidm,
                hrfluidm,
                phim,
                XWsolidm0,
                mode,
            )
        end
        # verification, from HTM-planetary.m, line 180ff;
        # HTM-hydration.m, line 154ff
        m=1 # Marker counter
        for jm in 1:1:Nxm
            for im in 1:1:Nym
                # Define marker coordinates
                xm_ver[m]=dxm/2+(jm-1)*dxm # +(rand-0.5)*dxm;
                ym_ver[m]=dym/2+(im-1)*dym # +(rand-0.5)*dym;
                # Marker properties
                rmark=((xm_ver[m]-xsize/2)^2+(ym_ver[m]-ysize/2)^2)^0.5
                if (rmark<rplanet)
                    # Planet
                    tm_ver[m]=1 # mantle
                    if (rmark>rcrust)
                        tm_ver[m]=2 # crust
                    end
                    tkm_ver[m]=tkm0[tm_ver[m]] # Temperature
                    phim_ver[m]=phim0 # *(1+1.0*(rand-0.5)); # Porosity
                    etavpm_ver[m]=etasolidm[tm_ver[m]]#*exp(-28*phim_ver[m]); % Matrix viscosity
                    XWsolidm0_ver[m] = XWsolidm_init[tm_ver[m]]
                else
                    # Sticky space (to have internal free surface)
                    tm_ver[m]=3 # Material type
                    tkm_ver[m]=tkm0[tm_ver[m]] # Temperature
                    phim_ver[m]=phimin # Porosity
                    etavpm_ver[m]=etasolidm[tm_ver[m]] # Matrix viscosity
                end
                # Update marker counter
                m=m+1
            end
        end
        # verification, from HTM-planetary.m, line 327ff
        # HTM-hydration.m, line 269ff
        for m in 1:1:marknum
            # Compute marker parameters
            if tm[m]<3
                rhocpfluidm_ver = Erebus.compute_rhocpfluidm(tkm[m], mode)
                ksolidm_ver = Erebus.compute_ksolidm(tkm[m], mode)
                kfluidm_ver = Erebus.compute_kfluidm(tkm[m], mode)
                # rocks
                # Compute density of solid and fluid
                XDsolidm0=1-XWsolidm0_ver[m]
                rhosolidm0=(MD+MH2O*XWsolidm0_ver[m])/(
                    VWsolid*XWsolidm0_ver[m]+VDsolid*XDsolidm0
                )
                if tkm[m] > tmfluidphase
                    rhofluidm0=rhoH2Ofluid
                else
                    rhofluidm0=rhoH2Ofluidice
                end
                kphim_ver[m] = (
                    kphim0[tm_ver[m]]*(phim_ver[m]/phim0)^3/((1-phim_ver[m])/(1-phim0))^2
                ) #Permeability
                rhototalm_ver[m] = (rhosolidm0*(1-phim_ver[m]) + rhofluidm0*phim_ver[m])
                rhocptotalm_ver[m] = (
                    rhocpsolidm[tm_ver[m]]*(1-phim_ver[m]) + rhocpfluidm_ver*phim_ver[m]
                )
                etasolidcur_ver[m] = etasolidm[tm_ver[m]]
                if tkm_ver[m]>tmsolidphase
                    etasolidcur_ver[m] = etasolidmm[tm_ver[m]]
                end
                hrtotalm_ver[m] =
                    hrsolidm[tm_ver[m]]*(1-phim_ver[m])+hrfluidm[tm_ver[m]]*phim_ver[m]
                ktotalm_ver[m] =
                    (
                        ksolidm_ver*kfluidm_ver/2+(
                            (ksolidm_ver*(3*phim_ver[m]-2) + kfluidm_ver*(1-3*phim_ver[m]))^2
                        )/16
                    )^0.5-(ksolidm_ver*(3*phim_ver[m]-2) + kfluidm_ver*(1-3*phim_ver[m]))/4
                gggtotalm_ver[m] = gggsolidm[tm_ver[m]]
                fricttotalm_ver[m] = frictsolidm[tm_ver[m]]
                cohestotalm_ver[m] = cohessolidm[tm_ver[m]]
                tenstotalm_ver[m] = tenssolidm[tm_ver[m]]
                etafluidcur_ver[m] = etafluidm[tm_ver[m]]
                rhofluidcur_ver[m] = rhofluidm[tm_ver[m]]
                if tkm_ver[m]>tmfluidphase
                    etafluidcur_ver[m] = etafluidmm[tm_ver[m]]
                end
                etatotalm_ver[m] = max(etamin, etafluidcur_ver[m], etasolidcur_ver[m])#*exp(-28*phim_ver[m])));
            else
                # Sticky air
                kphim_ver[m] =
                    kphim0[tm_ver[m]]*(phim_ver[m]/phim0)^3/((1-phim_ver[m])/(1-phim0))^2 #Permeability
                rhototalm_ver[m] = rhosolidm[tm_ver[m]]
                rhocptotalm_ver[m] = rhocpsolidm[tm_ver[m]]
                etatotalm_ver[m] = etasolidm[tm_ver[m]]
                hrtotalm_ver[m] = hrsolidm[tm_ver[m]]
                ktotalm_ver[m] = ksolidm[tm_ver[m]]
                gggtotalm_ver[m] = gggsolidm[tm_ver[m]]
                fricttotalm_ver[m] = frictsolidm[tm_ver[m]]
                cohestotalm_ver[m] = cohessolidm[tm_ver[m]]
                tenstotalm_ver[m] = tenssolidm[tm_ver[m]]
                rhofluidcur_ver[m] = rhofluidm[tm_ver[m]]
                etafluidcur_ver[m] = etafluidm[tm_ver[m]]
            end
        end
        # test
        @test xm == xm_ver
        @test ym == ym_ver
        @test tm == tm_ver
        @test tkm == tkm_ver
        @test phim == phim_ver
        @test etavpm == etavpm_ver
        @test rhototalm == rhototalm_ver
        @test rhocptotalm == rhocptotalm_ver
        @test hrtotalm == hrtotalm_ver
        @test ktotalm == ktotalm_ver
        @test inv_gggtotalm == inv.(gggtotalm_ver)
        @test fricttotalm == fricttotalm_ver
        @test cohestotalm == cohestotalm_ver
        @test tenstotalm == tenstotalm_ver
        @test rhofluidcur == rhofluidcur_ver
        @test etatotalm == etatotalm_ver
        @test tkm_rhocptotalm ≈ tkm_ver .* rhocptotalm_ver
        @test etafluidcur_inv_kphim ≈ etafluidcur_ver ./ kphim_ver

        # Test compute_marker_properties! with rhofluidcur integration
        rhof_test = zeros(3)
        tm_test = [1, 2, 3] # core (cold ice), crust (hot water), sticky air
        tk_test = [200.0, 373.0, 300.0]
        rhotot_test = zeros(3)
        rhocptot_test = zeros(3)
        etatot_test = zeros(3)
        hrtot_test = zeros(3)
        ktot_test = zeros(3)
        tkm_rhocptot_test = zeros(3)
        eta_inv_k_test = zeros(3)
        phi_test = [0.1, 0.1, 1.0]
        xw_test = [0.0, 0.0, 0.0]

        for idx in 1:3
            Erebus.compute_marker_properties!(
                idx,
                tm_test,
                tk_test,
                rhotot_test,
                rhocptot_test,
                etatot_test,
                hrtot_test,
                ktot_test,
                tkm_rhocptot_test,
                eta_inv_k_test,
                hrsolidm,
                hrfluidm,
                phi_test,
                xw_test,
                mode,
                rhof_test;
                thermal_buoyancy=true,
                alphafluid=2.0e-4,
                tmfluidphase_val=273.0,
            )
        end
        # Sub-freezing marker (T = 200 K <= 273 K) -> ice density
        @test rhof_test[1] == 917.0
        # Hot rock marker (T = 373 K > 273 K, ΔT = 100 K) -> thermal expansion
        expected_hot_rhof = 1000.0 * (1.0 - 2.0e-4 * 100.0) # 980.0 kg/m³
        @test rhof_test[2] ≈ expected_hot_rhof rtol=1e-12
        # Sticky air marker -> sticky air fluid density
        @test rhof_test[3] == 1.0

        # With thermal_buoyancy = false
        for idx in 1:3
            Erebus.compute_marker_properties!(
                idx,
                tm_test,
                tk_test,
                rhotot_test,
                rhocptot_test,
                etatot_test,
                hrtot_test,
                ktot_test,
                tkm_rhocptot_test,
                eta_inv_k_test,
                hrsolidm,
                hrfluidm,
                phi_test,
                xw_test,
                mode,
                rhof_test;
                thermal_buoyancy=false,
                alphafluid=2.0e-4,
                tmfluidphase_val=273.0,
            )
        end
        @test rhof_test[1] == 917.0
        @test rhof_test[2] == 1000.0
        @test rhof_test[3] == 1.0

        # Fluid viscosity tests with compute_marker_properties!
        for idx in 1:3
            Erebus.compute_marker_properties!(
                idx,
                tm_test,
                tk_test,
                rhotot_test,
                rhocptot_test,
                etatot_test,
                hrtot_test,
                ktot_test,
                tkm_rhocptot_test,
                eta_inv_k_test,
                hrsolidm,
                hrfluidm,
                phi_test,
                xw_test,
                mode,
                rhof_test;
                fluid_viscosity_mode=:arrhenius,
                fluid_viscosity_Ea=15.0e3,
                fluid_viscosity_T0=293.15,
                fluid_viscosity_eta0=1.0e-3,
                tmfluidphase_val=273.0,
            )
        end
        # Marker 1 (T = 200 K): frozen, eta_fluid = 1e12
        expected_kphi_1 = Erebus.kphi(kphim0[tm_test[1]], phi_test[1])
        @test eta_inv_k_test[1] ≈ 1.0e12 / expected_kphi_1 rtol=1e-12

        # Marker 2 (T = 373 K): hot fluid with Arrhenius viscosity decrease
        expected_eta_2 = Erebus.compute_fluid_viscosity(
            373.0, 2; mode=:arrhenius, Ea=15.0e3, T0=293.15, eta0=1.0e-3, tmfluidphase=273.0
        )
        expected_kphi_2 = Erebus.kphi(kphim0[tm_test[2]], phi_test[2])
        @test eta_inv_k_test[2] ≈ expected_eta_2 / expected_kphi_2 rtol=1e-12
        @test expected_eta_2 < 1.0e-3

        # Mode :constant
        for idx in 1:3
            Erebus.compute_marker_properties!(
                idx,
                tm_test,
                tk_test,
                rhotot_test,
                rhocptot_test,
                etatot_test,
                hrtot_test,
                ktot_test,
                tkm_rhocptot_test,
                eta_inv_k_test,
                hrsolidm,
                hrfluidm,
                phi_test,
                xw_test,
                mode,
                rhof_test;
                fluid_viscosity_mode=:constant,
                fluid_viscosity_eta0=1.0e-3,
                tmfluidphase_val=273.0,
            )
        end
        @test eta_inv_k_test[2] ≈ 1.0e-3 / expected_kphi_2 rtol=1e-12
    end # testset "define_markers!() & compute_marker_properties!()"

    @testset "update_marker_viscosity!()" begin
        marknum = start_marknum
        xm = rand(rgen, (-dx):0.1:(x[end] + dx), marknum)
        ym = rand(rgen, (-dy):0.1:(y[end] + dy), marknum)
        tm = rand(rgen, 1:3, marknum)
        tkm = rand(rgen, (tmfluidphase - 100):0.1:(tmsolidphase + 100), marknum)
        etatotalm = [etasolidm[tm[m]] for m in 1:1:marknum]
        ETA = rand(rgen, Ny, Nx)
        YNY = rand(rgen, Bool, Ny, Nx)
        YNY_inv_ETA = YNY ./ ETA
        etavpm = zeros(marknum)
        etavpm_ver = zero(etavpm)
        # update marker Viscosity
        for m in 1:1:marknum
            Erebus.update_marker_viscosity!(
                m, xm, ym, tm, tkm, etatotalm, etavpm, YNY, YNY_inv_ETA
            )
        end
        # verification, from HTM-planetary.m, line 1321ff
        for m in 1:1:marknum
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
            wtmi1j=(1-dxmj/dx)*(dymi/dy)
            wtmij1=(dxmj/dx)*(1-dymi/dy)
            wtmi1j1=(dxmj/dx)*(dymi/dy)
            # Matrix viscosity
            if tm[m]<3
                # Rocks
                etasolidcur_ver=etasolidm[tm[m]]
                if tkm[m]>tmsolidphase
                    etasolidcur_ver=etasolidmm[tm[m]]
                end
                etatotalm_ver=etasolidcur_ver#*exp(-28*phim[m])
            else
                # Sticky air
                etatotalm_ver=etasolidm[tm[m]]
            end
            if YNY[i, j]>0 || YNY[i + 1, j]>0 || YNY[i, j + 1]>0 || YNY[i + 1, j + 1]>0
                #         etavpm[m]=ETA[i,j]*wtmij+ETA[i+1,j]*wtmi1j+...
                #                 ETA[i,j+1]*wtmij1+ETA[i+1,j+1]*wtmi1j1
                #         etavpm[m]=1/(1/ETA[i,j]*wtmij+1/ETA[i+1,j]*wtmi1j+...
                #                 1/ETA[i,j+1]*wtmij1+1/ETA[i+1,j+1]*wtmi1j1)
                etavpm_ver[m]=1/(
                    YNY[i, j]/ETA[i, j]*wtmij+YNY[i + 1, j]/ETA[i + 1, j]*wtmi1j+YNY[
                        i, j + 1
                    ]/ETA[i, j + 1]*wtmij1+YNY[i + 1, j + 1]/ETA[i + 1, j + 1]*wtmi1j1
                )
                if etavpm_ver[m]>=etatotalm_ver
                    etavpm_ver[m]=etatotalm_ver
                end
            else
                etavpm_ver[m]=etatotalm_ver
            end
        end
        # test
        @test etavpm ≈ etavpm_ver rtol=1e-9
    end # testset "update_marker_viscosity!()"

    @testset "fix()" begin
        @testset "basic nodes" begin
            for yy in (-dy):(dy / 3):(ysize + dy), xx in (-dx):(dx / 3):(xsize + dx)
                i, j = Erebus.fix(
                    xx, yy, x, y, dx, dy, jmin_basic, jmax_basic, imin_basic, imax_basic
                )
                # verification, from HTM-planetary.m, line 395ff
                j_ver=trunc(Int, (xx-x[1])/dx)+1
                i_ver=trunc(Int, (yy-y[1])/dy)+1
                if j_ver<1
                    j_ver=1
                elseif j_ver>Nx-1
                    j_ver=Nx-1
                end
                if i_ver<1
                    i_ver=1
                elseif i_ver>Ny-1
                    i_ver=Ny-1
                end
                @test i == i_ver
                @test j == j_ver
            end
        end
        @testset "Vx nodes" begin
            for yy in (-dy):(dy / 3):(ysize + dy), xx in (-dx):(dx / 3):(xsize + dx)
                i, j = Erebus.fix(
                    xx, yy, xvx, yvx, dx, dy, jmin_vx, jmax_vx, imin_vx, imax_vx
                )
                # verification, from HTM-planetary.m, line 455ff
                j_ver=trunc(Int, (xx-xvx[1])/dx)+1
                i_ver=trunc(Int, (yy-yvx[1])/dy)+1
                if j_ver<1
                    j_ver=1
                elseif j_ver>Nx-1
                    j_ver=Nx-1
                end
                if i_ver<1
                    i_ver=1
                elseif i_ver>Ny
                    i_ver=Ny
                end
                @test i == i_ver
                @test j == j_ver
            end
        end
        @testset "Vy nodes" begin
            for yy in (-dy):(dy / 3):(ysize + dy), xx in (-dx):(dx / 3):(xsize + dx)
                i, j = Erebus.fix(
                    xx, yy, xvy, yvy, dx, dy, jmin_vy, jmax_vy, imin_vy, imax_vy
                )
                # verification, from HTM-planetary.m, line 506ff
                j_ver=trunc(Int, (xx-xvy[1])/dx)+1
                i_ver=trunc(Int, (yy-yvy[1])/dy)+1
                if j_ver<1
                    j_ver=1
                elseif j_ver>Nx
                    j_ver=Nx
                end
                if i_ver<1
                    i_ver=1
                elseif i_ver>Ny-1
                    i_ver=Ny-1
                end
                @test i == i_ver
                @test j == j_ver
            end
        end
        @testset "P nodes" begin
            for yy in (-dy):(dy / 3):(ysize + dy), xx in (-dx):(dx / 3):(xsize + dx)
                i, j = Erebus.fix(xx, yy, xp, yp, dx, dy, jmin_p, jmax_p, imin_p, imax_p)
                # verification, from HTM-planetary.m, line 558ff
                j_ver=trunc(Int, (xx-xp[1])/dx)+1
                i_ver=trunc(Int, (yy-yp[1])/dy)+1
                if j_ver<1
                    j_ver=1
                elseif j_ver>Nx
                    j_ver=Nx
                end
                if i_ver<1
                    i_ver=1
                elseif i_ver>Ny
                    i_ver=Ny
                end
                @test i == i_ver
                @test j == j_ver
            end
        end
    end # testset "fix"

    @testset "fix_distances()" begin
        @testset "basic nodes" begin
            for yy in (-dy):(dy / 3):(ysize + dy), xx in (-dx):(dx / 3):(xsize + dx)
                i, j, dxmj, dymi = Erebus.fix_distances(
                    xx, yy, x, y, dx, dy, jmin_basic, jmax_basic, imin_basic, imax_basic
                )
                # verification, from HTM-planetary.m, line 395ff
                j_ver=trunc(Int, (xx-x[1])/dx)+1
                i_ver=trunc(Int, (yy-y[1])/dy)+1
                if j_ver<1
                    j_ver=1
                elseif j_ver>Nx-1
                    j_ver=Nx-1
                end
                if i_ver<1
                    i_ver=1
                elseif i_ver>Ny-1
                    i_ver=Ny-1
                end
                dxmj_ver=xx-x[j_ver]
                dymi_ver=yy-y[i_ver]
                @test dxmj == dxmj_ver
                @test dymi == dymi_ver
            end
        end
        @testset "Vx nodes" begin
            for yy in (-dy):(dy / 3):(ysize + dy), xx in (-dx):(dx / 3):(xsize + dx)
                i, j, dxmj, dymi = Erebus.fix_distances(
                    xx, yy, xvx, yvx, dx, dy, jmin_vx, jmax_vx, imin_vx, imax_vx
                )
                # verification, from HTM-planetary.m, line 455ff
                j_ver=trunc(Int, (xx-xvx[1])/dx)+1
                i_ver=trunc(Int, (yy-yvx[1])/dy)+1
                if j_ver<1
                    j_ver=1
                elseif j_ver>Nx-1
                    j_ver=Nx-1
                end
                if i_ver<1
                    i_ver=1
                elseif i_ver>Ny
                    i_ver=Ny
                end
                dxmj_ver=xx-xvx[j_ver]
                dymi_ver=yy-yvx[i_ver]
                @test dxmj == dxmj_ver
                @test dymi == dymi_ver
            end
        end
        @testset "Vy nodes" begin
            for yy in (-dy):(dy / 3):(ysize + dy), xx in (-dx):(dx / 3):(xsize + dx)
                i, j, dxmj, dymi = Erebus.fix_distances(
                    xx, yy, xvy, yvy, dx, dy, jmin_vy, jmax_vy, imin_vy, imax_vy
                )
                # verification, from HTM-planetary.m, line 506ff
                j_ver=trunc(Int, (xx-xvy[1])/dx)+1
                i_ver=trunc(Int, (yy-yvy[1])/dy)+1
                if j_ver<1
                    j_ver=1
                elseif j_ver>Nx
                    j_ver=Nx
                end
                if i_ver<1
                    i_ver=1
                elseif i_ver>Ny-1
                    i_ver=Ny-1
                end
                dxmj_ver=xx-xvy[j_ver]
                dymi_ver=yy-yvy[i_ver]
                @test dxmj == dxmj_ver
                @test dymi == dymi_ver
            end
        end
        @testset "P nodes" begin
            for yy in (-dy):(dy / 3):(ysize + dy), xx in (-dx):(dx / 3):(xsize + dx)
                i, j, dxmj, dymi = Erebus.fix_distances(
                    xx, yy, xp, yp, dx, dy, jmin_p, jmax_p, imin_p, imax_p
                )
                # verification, from HTM-planetary.m, line 558ff
                j_ver=trunc(Int, (xx-xp[1])/dx)+1
                i_ver=trunc(Int, (yy-yp[1])/dy)+1
                if j_ver<1
                    j_ver=1
                elseif j_ver>Nx
                    j_ver=Nx
                end
                if i_ver<1
                    i_ver=1
                elseif i_ver>Ny
                    i_ver=Ny
                end
                dxmj_ver=xx-xp[j_ver]
                dymi_ver=yy-yp[i_ver]
                @test dxmj == dxmj_ver
                @test dymi == dymi_ver
            end
        end
    end # testset "fix_distances"

    @testset "fix_weights() elementary" begin
        @testset "basic nodes" begin
            # verification, from HTM-planetary.m, line 373ff
            jmin, jmax = jmin_basic, jmax_basic
            imin, imax = imin_basic, imax_basic
            function fix_basic(xm, ym, x_axis, y_axis, dx, dy)
                j=trunc(Int, (xm-x_axis[1])/dx)+1
                i=trunc(Int, (ym-y_axis[1])/dy)+1
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
                dxmj=xm-x[j]
                dymi=ym-y[i]
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
            end
            # top left
            xm = -x[1]
            ym = -y[1]
            @test Erebus.fix_weights(xm, ym, x, y, dx, dy, jmin, jmax, imin, imax) ==
                fix_basic(xm, ym, x, y, dx, dy)
            # bottom left
            xm = x[1]
            ym = y[end]
            @test Erebus.fix_weights(xm, ym, x, y, dx, dy, jmin, jmax, imin, imax) ==
                fix_basic(xm, ym, x, y, dx, dy)
            # top right
            xm = x[end]
            ym = y[1]
            j=trunc(Int, (xm-x[1])/dx)+1
            i=trunc(Int, (ym-y[1])/dy)+1
            @test Erebus.fix_weights(xm, ym, x, y, dx, dy, jmin, jmax, imin, imax) ==
                fix_basic(xm, ym, x, y, dx, dy)
            # bottom right
            xm = x[end]
            ym = y[end]
            @test Erebus.fix_weights(xm, ym, x, y, dx, dy, jmin, jmax, imin, imax) ==
                fix_basic(xm, ym, x, y, dx, dy)
        end # testset "basic nodes"

        @testset "Vx nodes" begin
            # verification, from HTM-planetary.m, line 434ff
            jmin, jmax = jmin_vx, jmax_vx
            imin, imax = imin_vx, imax_vx
            function fix_vx(xm, ym, x_axis, y_axis, dx, dy)
                j=trunc(Int, (xm-x_axis[1])/dx)+1
                i=trunc(Int, (ym-y_axis[1])/dy)+1
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
                dxmj=xm-xvx[j]
                dymi=ym-yvx[i]
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
            end
            # top left
            xm = -xvx[1]
            ym = -yvx[1]
            @test Erebus.fix_weights(xm, ym, xvx, yvx, dx, dy, jmin, jmax, imin, imax) ==
                fix_vx(xm, ym, xvx, yvx, dx, dy)
            # bottom left
            xm = xvx[1]
            ym = yvx[end]
            @test Erebus.fix_weights(xm, ym, xvx, yvx, dx, dy, jmin, jmax, imin, imax) ==
                fix_vx(xm, ym, xvx, yvx, dx, dy)
            # top right
            xm = xvx[end]
            ym = yvx[1]
            @test Erebus.fix_weights(xm, ym, xvx, yvx, dx, dy, jmin, jmax, imin, imax) ==
                fix_vx(xm, ym, xvx, yvx, dx, dy)
            # bottom right
            xm = xvx[end]
            ym = yvx[end]
            @test Erebus.fix_weights(xm, ym, xvx, yvx, dx, dy, jmin, jmax, imin, imax) ==
                fix_vx(xm, ym, xvx, yvx, dx, dy)
        end # testset "Vx nodes"

        @testset "Vy nodes" begin
            # verification, from HTM-planetary.m, line 484ff
            jmin, jmax = jmin_vy, jmax_vy
            imin, imax = imin_vy, imax_vy
            function fix_vy(xm, ym, x_axis, y_axis, dx, dy)
                j=trunc(Int, (xm-x_axis[1])/dx)+1
                i=trunc(Int, (ym-y_axis[1])/dy)+1
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
                dxmj=xm-xvy[j]
                dymi=ym-yvy[i]
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
            end
            # top left
            xm = -xvy[1]
            ym = -yvy[1]
            @test Erebus.fix_weights(xm, ym, xvy, yvy, dx, dy, jmin, jmax, imin, imax) ==
                fix_vy(xm, ym, xvy, yvy, dx, dy)
            # bottom left
            xm = xvy[1]
            ym = yvy[end]
            @test Erebus.fix_weights(xm, ym, xvy, yvy, dx, dy, jmin, jmax, imin, imax) ==
                fix_vy(xm, ym, xvy, yvy, dx, dy)
            # top right
            xm = xvy[end]
            ym = yvy[1]
            @test Erebus.fix_weights(xm, ym, xvy, yvy, dx, dy, jmin, jmax, imin, imax) ==
                fix_vy(xm, ym, xvy, yvy, dx, dy)
            # bottom right
            xm = xvy[end]
            ym = yvy[end]
            @test Erebus.fix_weights(xm, ym, xvy, yvy, dx, dy, jmin, jmax, imin, imax) ==
                fix_vy(xm, ym, xvy, yvy, dx, dy)
        end # testset "Vy nodes"

        @testset "P nodes" begin
            # verification, from HTM-planetary.m, line 538ff
            jmin, jmax = jmin_p, jmax_p
            imin, imax = imin_p, imax_p
            function fix_p(xm, ym, x_axis, y_axis, dx, dy)
                j=trunc(Int, (xm-x_axis[1])/dx)+1
                i=trunc(Int, (ym-y_axis[1])/dy)+1
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
                dxmj=xm-xp[j]
                dymi=ym-yp[i]
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
            end
            # top left
            xm = -xp[1]
            ym = -yp[1]
            @test Erebus.fix_weights(xm, ym, xp, yp, dx, dy, jmin, jmax, imin, imax) ==
                fix_p(xm, ym, xp, yp, dx, dy)
            # bottom left
            xm = xp[1]
            ym = yp[end]
            @test Erebus.fix_weights(xm, ym, xp, yp, dx, dy, jmin, jmax, imin, imax) ==
                fix_p(xm, ym, xp, yp, dx, dy)
            # top right
            xm = xp[end]
            ym = yp[1]
            @test Erebus.fix_weights(xm, ym, xp, yp, dx, dy, jmin, jmax, imin, imax) ==
                fix_p(xm, ym, xp, yp, dx, dy)
            # bottom right
            xm = xp[end]
            ym = yp[end]
            @test Erebus.fix_weights(xm, ym, xp, yp, dx, dy, jmin, jmax, imin, imax) ==
                fix_p(xm, ym, xp, yp, dx, dy)
        end # testset "P nodes"    
    end # testset "fix_weights() elementary"

    @testset "fix_weights() advanced" begin
        # simulating markers
        marknum = 10_000
        @testset "basic nodes" begin
            # verification, from HTM-planetary.m, line 373ff
            jmin, jmax = jmin_basic, jmax_basic
            imin, imax = imin_basic, imax_basic
            function fix_basic(xm, ym, x_axis, y_axis, dx, dy)
                j=trunc(Int, (xm-x_axis[1])/dx)+1
                i=trunc(Int, (ym-y_axis[1])/dy)+1
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
                dxmj=xm-x[j]
                dymi=ym-y[i]
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
            end
            # simulating markers
            xm = rand(rgen, (-x[1]):0.1:(x[end] + dx), marknum)
            ym = rand(rgen, (-y[1]):0.1:(y[end] + dy), marknum)
            for m in 1:marknum
                i, j, weights = Erebus.fix_weights(
                    xm[m], ym[m], x, y, dx, dy, jmin, jmax, imin, imax
                )
                i_ver, j_ver, weights_ver = fix_basic(xm[m], ym[m], x, y, dx, dy)
                @test i == i_ver
                @test j == j_ver
                @test weights ≈ weights_ver rtol=1e-9
            end
        end # testset "basic nodes"
        @testset "Vx nodes" begin
            # verification, from HTM-planetary.m, line 434ff
            jmin, jmax = jmin_vx, jmax_vx
            imin, imax = imin_vx, imax_vx
            function fix_vx(xm, ym, x_axis, y_axis, dx, dy)
                j=trunc(Int, (xm-x_axis[1])/dx)+1
                i=trunc(Int, (ym-y_axis[1])/dy)+1
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
                dxmj=xm-xvx[j]
                dymi=ym-yvx[i]
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
            end
            # simulating markers
            xm = rand(rgen, (-xvx[1]):0.1:(xvx[end] + dx), marknum)
            ym = rand(rgen, (-yvx[1]):0.1:(yvx[end] + dy), marknum)
            for m in 1:marknum
                i, j, weights = Erebus.fix_weights(
                    xm[m], ym[m], xvx, yvx, dx, dy, jmin, jmax, imin, imax
                )
                i_ver, j_ver, weights_ver = fix_vx(xm[m], ym[m], xvx, yvx, dx, dy)
                @debug "fix_weights Vx" i i_ver j j_ver weights weights_ver
                @test i == i_ver
                @test j == j_ver
                @test weights ≈ weights_ver rtol=1e-9
            end
        end # testset "Vx nodes"
        @testset "Vy nodes" begin
            # verification, from HTM-planetary.m, line 484ff
            jmin, jmax = jmin_vy, jmax_vy
            imin, imax = imin_vy, imax_vy
            function fix_vy(xm, ym, x_axis, y_axis, dx, dy)
                j=trunc(Int, (xm-x_axis[1])/dx)+1
                i=trunc(Int, (ym-y_axis[1])/dy)+1
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
                dxmj=xm-xvy[j]
                dymi=ym-yvy[i]
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
            end
            # simulating markers
            xm = rand(rgen, (-xvy[1]):0.1:(xvy[end] + dx), marknum)
            ym = rand(rgen, (-yvy[1]):0.1:(yvy[end] + dy), marknum)
            for m in 1:marknum
                i, j, weights = Erebus.fix_weights(
                    xm[m], ym[m], xvy, yvy, dx, dy, jmin, jmax, imin, imax
                )
                i_ver, j_ver, weights_ver = fix_vy(xm[m], ym[m], xvy, yvy, dx, dy)
                @test i == i_ver
                @test j == j_ver
                @test weights ≈ weights_ver rtol=1e-9
            end
        end # testset "Vy nodes"
        @testset "P nodes" begin
            # verification, from HTM-planetary.m, line 538ff
            jmin, jmax = jmin_p, jmax_p
            imin, imax = imin_p, imax_p
            function fix_p(xm, ym, x_axis, y_axis, dx, dy)
                j=trunc(Int, (xm-x_axis[1])/dx)+1
                i=trunc(Int, (ym-y_axis[1])/dy)+1
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
                dxmj=xm-xp[j]
                dymi=ym-yp[i]
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
            end
            # simulating markers
            xm = rand(rgen, (-xp[1]):0.1:(xp[end] + dx), marknum)
            ym = rand(rgen, (-yp[1]):0.1:(yp[end] + dy), marknum)
            for m in 1:marknum
                i, j, weights = Erebus.fix_weights(
                    xm[m], ym[m], xp, yp, dx, dy, jmin, jmax, imin, imax
                )
                i_ver, j_ver, weights_ver = fix_p(xm[m], ym[m], xp, yp, dx, dy)
                @test i == i_ver
                @test j == j_ver
                @test weights ≈ weights_ver rtol=1e-9
            end
        end # testset "P nodes"    
    end # testset "fix_weights() advanced"

    @testset "reduce_add_3darray!()" begin
        A = rand(rgen, 100, 100, 10)
        A_ver = deepcopy(A)
        # reduce-sum A along third axis
        Erebus.reduce_add_3darray!(A)
        # verification
        A_ver = reduce(+, A_ver, dims=3)
        # test
        @test A[:, :, 1] ≈ A_ver rtol=1e-9
    end # testset "reduce_add_3darray!"

    @testset "interpolate_add_to_grid!()" begin
        # simulate markers
        marknum = start_marknum + 10_000
        xm = zeros(marknum)
        ym = zeros(marknum)
        for jm in 1:1:Nxm, im in 1:1:Nym
            # calculate marker counter
            m = (jm-1) * Nym + im
            # define marker coordinates
            xm[m] = dxm/2 + (jm-1) * dxm
            ym[m] = dym/2 + (im-1) * dym
        end
        xm[(start_marknum + 1):end] = rand(rgen, (-dx):0.1:(xsize + dx), 10_000)
        ym[(start_marknum + 1):end] = rand(rgen, (-dy):0.1:(ysize + dy), 10_000)
        property = rand(rgen, marknum)
        @testset "basic nodes" begin
            # sample interpolation array
            gridsum = zeros(Ny, Nx)
            weightsum = zeros(Ny, Nx)
            weightsum_2 = zeros(Ny, Nx)
            gridsum_ver = zeros(Ny, Nx)
            weightsum_ver = zeros(Ny, Nx)
            for m in 1:1:marknum
                i, j, weights = Erebus.fix_weights(
                    xm[m],
                    ym[m],
                    x,
                    y,
                    dx,
                    dy,
                    jmin_basic,
                    jmax_basic,
                    imin_basic,
                    imax_basic,
                )
                Erebus.interpolate_add_to_grid!(i, j, weights, property[m], gridsum)
                Erebus.interpolate_add_to_grid!(i, j, weights, one(1.0), weightsum)
                # verification
                j_ver=trunc(Int, (xm[m]-x[1])/dx)+1
                i_ver=trunc(Int, (ym[m]-y[1])/dy)+1
                if j_ver<1
                    j_ver=1
                elseif j_ver>Nx-1
                    j_ver=Nx-1
                end
                if i_ver<1
                    i_ver=1
                elseif i_ver>Ny-1
                    i_ver=Ny-1
                end
                # Compute distances
                dxmj=xm[m]-x[j_ver]
                dymi=ym[m]-y[i_ver]
                # Compute weights
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                gridsum_ver[i_ver, j_ver] += property[m]*wtmij
                weightsum_ver[i_ver, j_ver] += wtmij
                gridsum_ver[i_ver + 1, j_ver] += property[m]*wtmi1j
                weightsum_ver[i_ver + 1, j_ver] += wtmi1j
                gridsum_ver[i_ver, j_ver + 1] += property[m]*wtmij1
                weightsum_ver[i_ver, j_ver + 1] += wtmij1
                gridsum_ver[i_ver + 1, j_ver + 1] += property[m]*wtmi1j1
                weightsum_ver[i_ver + 1, j_ver + 1] += wtmi1j1
                if weightsum[i, j]!=weightsum_ver[i, j]
                    @show m i j i_ver j_ver weights wtmij wtmi1j wtmij1 wtmi1j1 weightsum[
                        i, j
                    ] weightsum_2[i, j] weightsum_ver[i, j] weightsum[i + 1, j] weightsum_ver[
                        i + 1, j
                    ] weightsum[i, j + 1] weightsum_ver[i, j + 1] weightsum[i + 1, j + 1] weightsum_ver[
                        i + 1, j + 1
                    ] gridsum[i, j] gridsum_ver[i, j] gridsum[i + 1, j] gridsum_ver[
                        i + 1, j
                    ] gridsum[i, j + 1] gridsum_ver[i, j + 1] gridsum[i + 1, j + 1] gridsum_ver[
                        i + 1, j + 1
                    ]
                end
            end
            # Test
            for i in 1:1:Nx, j in 1:1:Ny
                @test gridsum[i, j] ≈ gridsum_ver[i, j] rtol=1e-9
                @test weightsum[i, j] ≈ weightsum_ver[i, j] rtol=1e-9
            end
        end
        @testset "Vx nodes" begin
            # sample interpolation array
            gridsum = zeros(Ny1, Nx1)
            weightsum = zeros(Ny1, Nx1)
            gridsum_ver = zeros(Ny1, Nx1)
            weightsum_ver = zeros(Ny1, Nx1)
            for m in 1:1:marknum
                i, j, weights = Erebus.fix_weights(
                    xm[m], ym[m], xvx, yvx, dx, dy, jmin_vx, jmax_vx, imin_vx, imax_vx
                )
                Erebus.interpolate_add_to_grid!(i, j, weights, property[m], gridsum)
                Erebus.interpolate_add_to_grid!(i, j, weights, one(1.0), weightsum)
                # verification
                j_ver=trunc(Int, (xm[m]-xvx[1])/dx)+1
                i_ver=trunc(Int, (ym[m]-yvx[1])/dy)+1
                if j_ver<1
                    j_ver=1
                elseif j_ver>Nx-1
                    j_ver=Nx-1
                end
                if i_ver<1
                    i_ver=1
                elseif i_ver>Ny
                    i_ver=Ny
                end
                # Compute distances
                dxmj=xm[m]-xvx[j]
                dymi=ym[m]-yvx[i]
                # Compute weights
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Update properties
                gridsum_ver[i_ver, j_ver] += property[m]*wtmij
                weightsum_ver[i_ver, j_ver] += wtmij
                gridsum_ver[i_ver + 1, j_ver] += property[m]*wtmi1j
                weightsum_ver[i_ver + 1, j_ver] += wtmi1j
                gridsum_ver[i_ver, j_ver + 1] += property[m]*wtmij1
                weightsum_ver[i_ver, j_ver + 1] += wtmij1
                gridsum_ver[i_ver + 1, j_ver + 1] += property[m]*wtmi1j1
                weightsum_ver[i_ver + 1, j_ver + 1] += wtmi1j1
            end
            # Test
            for i in 1:1:Nx1, j in 1:1:Ny1
                @test gridsum[i, j] ≈ gridsum_ver[i, j] rtol=1e-9
                @test weightsum[i, j] ≈ weightsum_ver[i, j] rtol=1e-9
            end
        end
        @testset "Vy nodes" begin
            # sample interpolation array
            gridsum = zeros(Ny1, Nx1)
            weightsum = zeros(Ny1, Nx1)
            gridsum_ver = zeros(Ny1, Nx1)
            weightsum_ver = zeros(Ny1, Nx1)
            for m in 1:1:marknum
                i, j, weights = Erebus.fix_weights(
                    xm[m], ym[m], xvy, yvy, dx, dy, jmin_vy, jmax_vy, imin_vy, imax_vy
                )
                Erebus.interpolate_add_to_grid!(i, j, weights, property[m], gridsum)
                Erebus.interpolate_add_to_grid!(i, j, weights, one(1.0), weightsum)
                # verification
                j_ver=trunc(Int, (xm[m]-xvy[1])/dx)+1
                i_ver=trunc(Int, (ym[m]-yvy[1])/dy)+1
                if j_ver<1
                    j_ver=1
                elseif j_ver>Nx
                    j_ver=Nx
                end
                if i_ver<1
                    i_ver=1
                elseif i_ver>Ny-1
                    i_ver=Ny-1
                end
                # Compute distances
                dxmj=xm[m]-xvy[j_ver]
                dymi=ym[m]-yvy[i_ver]
                # Compute weights
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Update properties
                gridsum_ver[i_ver, j_ver] += property[m]*wtmij
                weightsum_ver[i_ver, j_ver] += wtmij
                gridsum_ver[i_ver + 1, j_ver] += property[m]*wtmi1j
                weightsum_ver[i_ver + 1, j_ver] += wtmi1j
                gridsum_ver[i_ver, j_ver + 1] += property[m]*wtmij1
                weightsum_ver[i_ver, j_ver + 1] += wtmij1
                gridsum_ver[i_ver + 1, j_ver + 1] += property[m]*wtmi1j1
                weightsum_ver[i_ver + 1, j_ver + 1] += wtmi1j1
            end
            # Test
            for i in 1:1:Nx1, j in 1:1:Ny1
                @test gridsum[i, j] ≈ gridsum_ver[i, j] rtol=1e-9
                @test weightsum[i, j] ≈ weightsum_ver[i, j] rtol=1e-9
            end
        end
        @testset "P nodes" begin
            # sample interpolation array
            gridsum = zeros(Ny1, Nx1)
            weightsum = zeros(Ny1, Nx1)
            gridsum_ver = zeros(Ny1, Nx1)
            weightsum_ver = zeros(Ny1, Nx1)
            for m in 1:1:marknum
                i, j, weights = Erebus.fix_weights(
                    xm[m], ym[m], xp, yp, dx, dy, jmin_p, jmax_p, imin_p, imax_p
                )
                Erebus.interpolate_add_to_grid!(i, j, weights, property[m], gridsum)
                Erebus.interpolate_add_to_grid!(i, j, weights, one(1.0), weightsum)
                # verification
                j_ver=trunc(Int, (xm[m]-xp[1])/dx)+1
                i_ver=trunc(Int, (ym[m]-yp[1])/dy)+1
                if j_ver<1
                    j_ver=1
                elseif j_ver>Nx
                    j_ver=Nx
                end
                if i_ver<1
                    i_ver=1
                elseif i_ver>Ny
                    i_ver=Ny
                end
                # Compute distances
                dxmj=xm[m]-xp[j_ver]
                dymi=ym[m]-yp[i_ver]
                # Compute weights
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Update properties
                gridsum_ver[i_ver, j_ver] += property[m]*wtmij
                weightsum_ver[i_ver, j_ver] += wtmij
                gridsum_ver[i_ver + 1, j_ver] += property[m]*wtmi1j
                weightsum_ver[i_ver + 1, j_ver] += wtmi1j
                gridsum_ver[i_ver, j_ver + 1] += property[m]*wtmij1
                weightsum_ver[i_ver, j_ver + 1] += wtmij1
                gridsum_ver[i_ver + 1, j_ver + 1] += property[m]*wtmi1j1
                weightsum_ver[i_ver + 1, j_ver + 1] += wtmi1j1
            end
            # Test
            for i in 1:1:Nx1, j in 1:1:Ny1
                @test gridsum[i, j] ≈ gridsum_ver[i, j] rtol=1e-9
                @test weightsum[i, j] ≈ weightsum_ver[i, j] rtol=1e-9
            end
        end
    end # testset "interpolate_add_to_grid

    @testset "interpolate_to_marker!()" begin
        jmin, jmax = jmin_basic, jmax_basic
        imin, imax = imin_basic, imax_basic
        marknum = 10_000
        xm = rand(rgen, (-x[1]):0.1:(x[end] + dx), marknum)
        ym = rand(rgen, (-y[1]):0.1:(y[end] + dy), marknum)
        property = zeros(marknum)
        grid = rand(rgen, Ny, Nx)
        # interpolate to markers & test
        for m in 1:1:marknum
            i, j, weights = Erebus.fix_weights(
                xm[m], ym[m], x, y, dx, dy, jmin, jmax, imin, imax
            )
            Erebus.interpolate_to_marker!(m, i, j, weights, property, grid)
            @test property[m] ≈ (
                grid[i, j] * weights[1] +
                grid[i + 1, j] * weights[2] +
                grid[i, j + 1] * weights[3] +
                grid[i + 1, j + 1] * weights[4]
            ) rtol=1e-9
        end
    end # testset "interpolate_to_marker!()"

    @testset "interpolate_add_to_marker!()" begin
        jmin, jmax = jmin_basic, jmax_basic
        imin, imax = imin_basic, imax_basic
        marknum = 10_000
        xm = rand(rgen, (-x[1]):0.1:(x[end] + dx), marknum)
        ym = rand(rgen, (-y[1]):0.1:(y[end] + dy), marknum)
        property = rand(rgen, marknum)
        property_ver = deepcopy(property)
        grid = rand(rgen, Ny, Nx)
        # interpolate to markers, verification, and test
        for m in 1:1:marknum
            i, j, weights = Erebus.fix_weights(
                xm[m], ym[m], x, y, dx, dy, jmin, jmax, imin, imax
            )
            Erebus.interpolate_add_to_marker!(m, i, j, weights, property, grid)
            @test property[m] ≈
                property_ver[m] + (
                grid[i, j] * weights[1] +
                grid[i + 1, j] * weights[2] +
                grid[i, j + 1] * weights[3] +
                grid[i + 1, j + 1] * weights[4]
            ) rtol=1e-9
        end
    end # testset "interpolate_add_to_marker!()"

    @testset "marker_to_basic_nodes!()" begin
        marknum = start_marknum
        (xm, ym, _, _, _, sxym, etavpm, _) = Erebus.setup_marker_properties(
            marknum, randomized=true
        )
        (rhototalm, rhocptotalm, etatotalm, _, _, _, _, inv_gggtotalm, fricttotalm, cohestotalm, tenstotalm, _, _, _) = Erebus.setup_marker_properties_helpers(
            marknum, randomized=true
        )
        (ETA0SUM, ETASUM, GGGSUM, SXYSUM, COHSUM, TENSUM, FRISUM, WTSUM, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _) = Erebus.setup_interpolated_properties()
        ETA0SUM_ver = zero(ETA0SUM)
        ETASUM_ver = zero(ETASUM)
        GGGSUM_ver = zero(GGGSUM)
        SXYSUM_ver = zero(SXYSUM)
        COHSUM_ver = zero(COHSUM)
        TENSUM_ver = zero(TENSUM)
        FRISUM_ver = zero(FRISUM)
        WTSUM_ver = zero(WTSUM)
        # interpolate markers to basic nodes
        for m in 1:1:marknum
            Erebus.marker_to_basic_nodes!(
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
            )
        end
        # verification
        for m in 1:1:marknum
            i, j, weights = Erebus.fix_weights(
                xm[m], ym[m], x, y, dx, dy, jmin_basic, jmax_basic, imin_basic, imax_basic
            )
            Erebus.interpolate_add_to_grid!(i, j, weights, etatotalm[m], ETA0SUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, etavpm[m], ETASUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, inv_gggtotalm[m], GGGSUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, sxym[m], SXYSUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, cohestotalm[m], COHSUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, tenstotalm[m], TENSUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, fricttotalm[m], FRISUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, 1.0, WTSUM_ver)
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
        marknum = start_marknum
        (xm, ym, _, _, _, _, _, phim) = Erebus.setup_marker_properties(
            marknum, randomized=true
        )
        (rhototalm, _, etatotalm, _, ktotalm, _, etafluidcur_inv_kphim, _, _, _, _, rhofluidcur, _, _) = Erebus.setup_marker_properties_helpers(
            marknum, randomized=true
        )
        (_, _, _, _, _, _, _, _, RHOXSUM, RHOFXSUM, KXSUM, PHIXSUM, RXSUM, WTXSUM, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _) = Erebus.setup_interpolated_properties()
        RHOXSUM_ver = zero(RHOXSUM)
        RHOFXSUM_ver = zero(RHOFXSUM)
        KXSUM_ver = zero(KXSUM)
        PHIXSUM_ver = zero(PHIXSUM)
        RXSUM_ver = zero(RXSUM)
        WTXSUM_ver = zero(WTXSUM)
        # interpolate markers to Vx nodes
        for m in 1:1:marknum
            Erebus.marker_to_vx_nodes!(
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
            )
        end
        # verification
        for m in 1:1:marknum
            i, j, weights = Erebus.fix_weights(
                xm[m], ym[m], xvx, yvx, dx, dy, jmin_vx, jmax_vx, imin_vx, imax_vx
            )
            Erebus.interpolate_add_to_grid!(i, j, weights, rhototalm[m], RHOXSUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, rhofluidcur[m], RHOFXSUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, ktotalm[m], KXSUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, phim[m], PHIXSUM_ver)
            Erebus.interpolate_add_to_grid!(
                i, j, weights, etafluidcur_inv_kphim[m], RXSUM_ver
            )
            Erebus.interpolate_add_to_grid!(i, j, weights, 1.0, WTXSUM_ver)
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
        marknum = start_marknum
        (xm, ym, _, _, _, _, _, phim) = Erebus.setup_marker_properties(
            marknum, randomized=true
        )
        (rhototalm, _, etatotalm, _, ktotalm, _, etafluidcur_inv_kphim, _, _, _, _, rhofluidcur, _, _) = Erebus.setup_marker_properties_helpers(
            marknum, randomized=true
        )
        (_, _, _, _, _, _, _, _, _, _, _, _, _, _, RHOYSUM, RHOFYSUM, KYSUM, PHIYSUM, RYSUM, WTYSUM, _, _, _, _, _, _, _, _, _, _) = Erebus.setup_interpolated_properties()
        RHOYSUM_ver = zero(RHOYSUM)
        RHOFYSUM_ver = zero(RHOFYSUM)
        KYSUM_ver = zero(KYSUM)
        PHIYSUM_ver = zero(PHIYSUM)
        RYSUM_ver = zero(RYSUM)
        WTYSUM_ver = zero(WTYSUM)
        # interpolate markers to Vy nodes
        for m in 1:1:marknum
            Erebus.marker_to_vy_nodes!(
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
            )
        end
        # verification
        for m in 1:1:marknum
            i, j, weights = Erebus.fix_weights(
                xm[m], ym[m], xvy, yvy, dx, dy, jmin_vy, jmax_vy, imin_vy, imax_vy
            )
            Erebus.interpolate_add_to_grid!(i, j, weights, rhototalm[m], RHOYSUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, rhofluidcur[m], RHOFYSUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, ktotalm[m], KYSUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, phim[m], PHIYSUM_ver)
            Erebus.interpolate_add_to_grid!(
                i, j, weights, etafluidcur_inv_kphim[m], RYSUM_ver
            )
            Erebus.interpolate_add_to_grid!(i, j, weights, 1.0, WTYSUM_ver)
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
        marknum = start_marknum
        (xm, ym, _, _, sxxm, _, _, phim) = Erebus.setup_marker_properties(
            marknum, randomized=true
        )
        (rhototalm, rhocptotalm, _, hrtotalm, _, tkm_rhocptotalm, _, inv_gggtotalm, _, _, _, _, alphasolidcur, alphafluidcur) = Erebus.setup_marker_properties_helpers(
            marknum, randomized=true
        )
        (_, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, _, GGGPSUM, SXXSUM, RHOSUM, RHOCPSUM, ALPHASUM, ALPHAFSUM, HRSUM, PHISUM, TKSUM, WTPSUM) = Erebus.setup_interpolated_properties()
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
        for m in 1:1:marknum
            Erebus.marker_to_p_nodes!(
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
            )
        end
        # verification
        for m in 1:1:marknum
            i, j, weights = Erebus.fix_weights(
                xm[m], ym[m], xp, yp, dx, dy, jmin_p, jmax_p, imin_p, imax_p
            )
            Erebus.interpolate_add_to_grid!(i, j, weights, inv_gggtotalm[m], GGGPSUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, sxxm[m], SXXSUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, rhototalm[m], RHOSUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, rhocptotalm[m], RHOCPSUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, alphafluidcur[m], ALPHAFSUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, hrtotalm[m], HRSUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, phim[m], PHISUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, tkm_rhocptotalm[m], TKSUM_ver)
            Erebus.interpolate_add_to_grid!(i, j, weights, 1.0, WTPSUM_ver)
        end
        # test
        @test WTPSUM == WTPSUM_ver
    end # testset "marker_to_p_nodes!()"

    @testset "compute node properties: basic, Vx, Vy, P, thermodynamic" begin
        # simulating markers
        marknum = 10_000
        @testset "compute_basic_node_properties!()" begin
            jmin, jmax = jmin_basic, jmax_basic
            imin, imax = imin_basic, imax_basic
            ETA0SUM = zeros(Ny, Nx)
            ETASUM = zeros(Ny, Nx)
            GGGSUM = zeros(Ny, Nx)
            SXYSUM = zeros(Ny, Nx)
            COHSUM = zeros(Ny, Nx)
            TENSUM = zeros(Ny, Nx)
            FRISUM = zeros(Ny, Nx)
            WTSUM = zeros(Ny, Nx)
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
            xm = rand(rgen, (-x[1]):0.1:(x[end] + dx), marknum)
            ym = rand(rgen, (-y[1]):0.1:(y[end] + dy), marknum)
            property = rand(rgen, 7, marknum)*1e6
            # calculate grid properties
            for m in 1:1:marknum
                i, j, weights = Erebus.fix_weights(
                    xm[m], ym[m], x, y, dx, dy, jmin, jmax, imin, imax
                )
                Erebus.interpolate_add_to_grid!(i, j, weights, property[1, m], ETA0SUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, property[2, m], ETASUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, inv(property[3, m]), GGGSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, property[4, m], SXYSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, property[5, m], COHSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, property[6, m], TENSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, property[7, m], FRISUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, 1.0, WTSUM)
            end
            Erebus.compute_basic_node_properties!(
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
                YNY,
            )
            # verification properties, from HTM-planetary.m, line 373ff, 606ff
            for m in 1:1:marknum
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
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Update properties
                # i;j Node
                ETA0SUM_ver[i, j]=ETA0SUM_ver[i, j]+property[1, m]*wtmij
                ETASUM_ver[i, j]=ETASUM_ver[i, j]+property[2, m]*wtmij
                GGGSUM_ver[i, j]=GGGSUM_ver[i, j]+1/property[3, m]*wtmij
                SXYSUM_ver[i, j]=SXYSUM_ver[i, j]+property[4, m]*wtmij
                COHSUM_ver[i, j]=COHSUM_ver[i, j]+property[5, m]*wtmij
                TENSUM_ver[i, j]=TENSUM_ver[i, j]+property[6, m]*wtmij
                FRISUM_ver[i, j]=FRISUM_ver[i, j]+property[7, m]*wtmij
                WTSUM_ver[i, j]=WTSUM_ver[i, j]+wtmij
                # i+1;j Node
                ETA0SUM_ver[i + 1, j]=ETA0SUM_ver[i + 1, j]+property[1, m]*wtmi1j
                ETASUM_ver[i + 1, j]=ETASUM_ver[i + 1, j]+property[2, m]*wtmi1j
                GGGSUM_ver[i + 1, j]=GGGSUM_ver[i + 1, j]+1/property[3, m]*wtmi1j
                SXYSUM_ver[i + 1, j]=SXYSUM_ver[i + 1, j]+property[4, m]*wtmi1j
                COHSUM_ver[i + 1, j]=COHSUM_ver[i + 1, j]+property[5, m]*wtmi1j
                TENSUM_ver[i + 1, j]=TENSUM_ver[i + 1, j]+property[6, m]*wtmi1j
                FRISUM_ver[i + 1, j]=FRISUM_ver[i + 1, j]+property[7, m]*wtmi1j
                WTSUM_ver[i + 1, j]=WTSUM_ver[i + 1, j]+wtmi1j
                # i;j+1 Node
                ETA0SUM_ver[i, j + 1]=ETA0SUM_ver[i, j + 1]+property[1, m]*wtmij1
                ETASUM_ver[i, j + 1]=ETASUM_ver[i, j + 1]+property[2, m]*wtmij1
                GGGSUM_ver[i, j + 1]=GGGSUM_ver[i, j + 1]+1/property[3, m]*wtmij1
                SXYSUM_ver[i, j + 1]=SXYSUM_ver[i, j + 1]+property[4, m]*wtmij1
                COHSUM_ver[i, j + 1]=COHSUM_ver[i, j + 1]+property[5, m]*wtmij1
                TENSUM_ver[i, j + 1]=TENSUM_ver[i, j + 1]+property[6, m]*wtmij1
                FRISUM_ver[i, j + 1]=FRISUM_ver[i, j + 1]+property[7, m]*wtmij1
                WTSUM_ver[i, j + 1]=WTSUM_ver[i, j + 1]+wtmij1
                # i+1;j+1 Node
                ETA0SUM_ver[i + 1, j + 1]=ETA0SUM_ver[i + 1, j + 1]+property[1, m]*wtmi1j1
                ETASUM_ver[i + 1, j + 1]=ETASUM_ver[i + 1, j + 1]+property[2, m]*wtmi1j1
                GGGSUM_ver[i + 1, j + 1]=GGGSUM_ver[i + 1, j + 1]+1/property[3, m]*wtmi1j1
                SXYSUM_ver[i + 1, j + 1]=SXYSUM_ver[i + 1, j + 1]+property[4, m]*wtmi1j1
                COHSUM_ver[i + 1, j + 1]=COHSUM_ver[i + 1, j + 1]+property[5, m]*wtmi1j1
                TENSUM_ver[i + 1, j + 1]=TENSUM_ver[i + 1, j + 1]+property[6, m]*wtmi1j1
                FRISUM_ver[i + 1, j + 1]=FRISUM_ver[i + 1, j + 1]+property[7, m]*wtmi1j1
                WTSUM_ver[i + 1, j + 1]=WTSUM_ver[i + 1, j + 1]+wtmi1j1
            end
            for j in 1:1:Nx
                for i in 1:1:Ny
                    if WTSUM_ver[i, j]>0
                        ETA0_ver[i, j]=ETA0SUM_ver[i, j]/WTSUM_ver[i, j]
                        ETA_ver[i, j]=ETASUM_ver[i, j]/WTSUM_ver[i, j]
                        if (ETA_ver[i, j]<ETA0_ver[i, j])
                            YNY_ver[i, j]=1
                        end
                        GGG_ver[i, j]=1/(GGGSUM_ver[i, j]/WTSUM_ver[i, j])
                        SXY0_ver[i, j]=SXYSUM_ver[i, j]/WTSUM_ver[i, j]
                        COH_ver[i, j]=COHSUM_ver[i, j]/WTSUM_ver[i, j]
                        TEN_ver[i, j]=TENSUM_ver[i, j]/WTSUM_ver[i, j]
                        FRI_ver[i, j]=FRISUM_ver[i, j]/WTSUM_ver[i, j]
                    end
                end
            end
            # test
            for j in 1:1:Nx, i in 1:1:Ny
                @test ETA0[i, j] ≈ ETA0_ver[i, j] rtol=1e-9
                @test ETA[i, j] ≈ ETA_ver[i, j] rtol=1e-9
                @test GGG[i, j] ≈ GGG_ver[i, j] rtol=1e-9
                @test SXY0[i, j] ≈ SXY0_ver[i, j] rtol=1e-9
                @test COH[i, j] ≈ COH_ver[i, j] rtol=1e-9
                @test TEN[i, j] ≈ TEN_ver[i, j] rtol=1e-9
                @test FRI[i, j] ≈ FRI_ver[i, j] rtol=1e-9
                @test YNY[i, j] ≈ YNY_ver[i, j] rtol=1e-9
            end
        end # testset "compute_basic_node_properties!()"

        @testset "compute_vx_node_properties!()" begin
            jmin, jmax = jmin_vx, jmax_vx
            imin, imax = imin_vx, imax_vx
            RHOXSUM = zeros(Ny1, Nx1)
            RHOFXSUM = zeros(Ny1, Nx1)
            KXSUM = zeros(Ny1, Nx1)
            PHIXSUM = zeros(Ny1, Nx1)
            RXSUM = zeros(Ny1, Nx1)
            WTXSUM = zeros(Ny1, Nx1)
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
            xm = rand(rgen, (-xvx[1]):0.1:(xvx[end] + dx), marknum)
            ym = rand(rgen, (-yvx[1]):0.1:(yvx[end] + dy), marknum)
            property = rand(rgen, 5, marknum)*1e6
            # calculate grid properties
            for m in 1:1:marknum
                i, j, weights = Erebus.fix_weights(
                    xm[m], ym[m], xvx, yvx, dx, dy, jmin, jmax, imin, imax
                )
                Erebus.interpolate_add_to_grid!(i, j, weights, property[1, m], RHOXSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, property[2, m], RHOFXSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, property[3, m], KXSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, property[4, m], PHIXSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, property[5, m], RXSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, 1.0, WTXSUM)
            end
            Erebus.compute_vx_node_properties!(
                RHOXSUM, RHOFXSUM, KXSUM, PHIXSUM, RXSUM, WTXSUM, RHOX, RHOFX, KX, PHIX, RX
            )
            # verification properties, from HTM-planetary.m, line 434ff, 624ff
            for m in 1:1:marknum
                j=trunc(Int, (xm[m]-xvx[1])/dx)+1
                i=trunc(Int, (ym[m]-yvx[1])/dy)+1
                if j<1
                    j=1
                elseif j>Nx-1
                    j=Nx-1
                end
                if (i<1)
                    i=1
                elseif i>Ny
                    i=Ny
                end
                # Compute distances
                dxmj=xm[m]-xvx[j]
                dymi=ym[m]-yvx[i]
                # Compute weights
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Update properties
                # i;j Node
                RHOXSUM_ver[i, j]=RHOXSUM_ver[i, j]+property[1, m]*wtmij
                RHOFXSUM_ver[i, j]=RHOFXSUM_ver[i, j]+property[2, m]*wtmij
                KXSUM_ver[i, j]=KXSUM_ver[i, j]+property[3, m]*wtmij
                PHIXSUM_ver[i, j]=PHIXSUM_ver[i, j]+property[4, m]*wtmij
                RXSUM_ver[i, j]=RXSUM_ver[i, j]+property[5, m]*wtmij
                WTXSUM_ver[i, j]=WTXSUM_ver[i, j]+wtmij
                # i+1;j Node
                RHOXSUM_ver[i + 1, j]=RHOXSUM_ver[i + 1, j]+property[1, m]*wtmi1j
                RHOFXSUM_ver[i + 1, j]=RHOFXSUM_ver[i + 1, j]+property[2, m]*wtmi1j
                KXSUM_ver[i + 1, j]=KXSUM_ver[i + 1, j]+property[3, m]*wtmi1j
                PHIXSUM_ver[i + 1, j]=PHIXSUM_ver[i + 1, j]+property[4, m]*wtmi1j
                RXSUM_ver[i + 1, j]=RXSUM_ver[i + 1, j]+property[5, m]*wtmi1j
                WTXSUM_ver[i + 1, j]=WTXSUM_ver[i + 1, j]+wtmi1j
                # i;j+1 Node
                RHOXSUM_ver[i, j + 1]=RHOXSUM_ver[i, j + 1]+property[1, m]*wtmij1
                RHOFXSUM_ver[i, j + 1]=RHOFXSUM_ver[i, j + 1]+property[2, m]*wtmij1
                KXSUM_ver[i, j + 1]=KXSUM_ver[i, j + 1]+property[3, m]*wtmij1
                PHIXSUM_ver[i, j + 1]=PHIXSUM_ver[i, j + 1]+property[4, m]*wtmij1
                RXSUM_ver[i, j + 1]=RXSUM_ver[i, j + 1]+property[5, m]*wtmij1
                WTXSUM_ver[i, j + 1]=WTXSUM_ver[i, j + 1]+wtmij1
                # i+1;j+1 Node
                RHOXSUM_ver[i + 1, j + 1]=RHOXSUM_ver[i + 1, j + 1]+property[1, m]*wtmi1j1
                RHOFXSUM_ver[i + 1, j + 1]=RHOFXSUM_ver[i + 1, j + 1]+property[2, m]*wtmi1j1
                KXSUM_ver[i + 1, j + 1]=KXSUM_ver[i + 1, j + 1]+property[3, m]*wtmi1j1
                PHIXSUM_ver[i + 1, j + 1]=PHIXSUM_ver[i + 1, j + 1]+property[4, m]*wtmi1j1
                RXSUM_ver[i + 1, j + 1]=RXSUM_ver[i + 1, j + 1]+property[5, m]*wtmi1j1
                WTXSUM_ver[i + 1, j + 1]=WTXSUM_ver[i + 1, j + 1]+wtmi1j1
            end
            for j in 1:1:Nx1
                for i in 1:1:Ny1
                    if (WTXSUM_ver[i, j]>0)
                        RHOX_ver[i, j]=RHOXSUM_ver[i, j]/WTXSUM_ver[i, j]
                        RHOFX_ver[i, j]=RHOFXSUM_ver[i, j]/WTXSUM_ver[i, j]
                        KX_ver[i, j]=KXSUM_ver[i, j]/WTXSUM_ver[i, j]
                        PHIX_ver[i, j]=PHIXSUM_ver[i, j]/WTXSUM_ver[i, j]
                        RX_ver[i, j]=RXSUM_ver[i, j]/WTXSUM_ver[i, j]
                    end
                end
            end
            # test
            for j in 1:1:Nx1, i in 1:1:Ny1
                @test RHOX[i, j] ≈ RHOX_ver[i, j] rtol=1e-9
                @test RHOFX[i, j] ≈ RHOFX_ver[i, j] rtol=1e-9
                @test KX[i, j] ≈ KX_ver[i, j] rtol=1e-9
                @test PHIX[i, j] ≈ PHIX_ver[i, j] rtol=1e-9
                @test RX[i, j] ≈ RX_ver[i, j] rtol=1e-9
            end
        end # testset "compute_vx_node_properties!()"

        @testset "compute_vy_node_properties!()" begin
            jmin, jmax = jmin_vy, jmax_vy
            imin, imax = imin_vy, imax_vy
            RHOYSUM = zeros(Ny1, Nx1)
            RHOFYSUM = zeros(Ny1, Nx1)
            KYSUM = zeros(Ny1, Nx1)
            PHIYSUM = zeros(Ny1, Nx1)
            RYSUM = zeros(Ny1, Nx1)
            WTYSUM = zeros(Ny1, Nx1)
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
            xm = rand(rgen, (-xvy[1]):0.1:(xvy[end] + dx), marknum)
            ym = rand(rgen, (-yvy[1]):0.1:(yvy[end] + dy), marknum)
            property = rand(rgen, 5, marknum)*1e6
            # calculate grid properties
            for m in 1:1:marknum
                i, j, weights = Erebus.fix_weights(
                    xm[m], ym[m], xvy, yvy, dx, dy, jmin, jmax, imin, imax
                )
                Erebus.interpolate_add_to_grid!(i, j, weights, property[1, m], RHOYSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, property[2, m], RHOFYSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, property[3, m], KYSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, property[4, m], PHIYSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, property[5, m], RYSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, 1.0, WTYSUM)
            end
            Erebus.compute_vy_node_properties!(
                RHOYSUM, RHOFYSUM, KYSUM, PHIYSUM, RYSUM, WTYSUM, RHOY, RHOFY, KY, PHIY, RY
            )
            # verification properties, from HTM-planetary.m, line 486ff, 636ff
            for m in 1:1:marknum
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
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Update properties
                # i;j Node
                RHOYSUM_ver[i, j]=RHOYSUM_ver[i, j]+property[1, m]*wtmij
                RHOFYSUM_ver[i, j]=RHOFYSUM_ver[i, j]+property[2, m]*wtmij
                KYSUM_ver[i, j]=KYSUM_ver[i, j]+property[3, m]*wtmij
                PHIYSUM_ver[i, j]=PHIYSUM_ver[i, j]+property[4, m]*wtmij
                RYSUM_ver[i, j]=RYSUM_ver[i, j]+property[5, m]*wtmij
                WTYSUM_ver[i, j]=WTYSUM_ver[i, j]+wtmij
                # i+1;j Node
                RHOYSUM_ver[i + 1, j]=RHOYSUM_ver[i + 1, j]+property[1, m]*wtmi1j
                RHOFYSUM_ver[i + 1, j]=RHOFYSUM_ver[i + 1, j]+property[2, m]*wtmi1j
                KYSUM_ver[i + 1, j]=KYSUM_ver[i + 1, j]+property[3, m]*wtmi1j
                PHIYSUM_ver[i + 1, j]=PHIYSUM_ver[i + 1, j]+property[4, m]*wtmi1j
                RYSUM_ver[i + 1, j]=RYSUM_ver[i + 1, j]+property[5, m]*wtmi1j
                WTYSUM_ver[i + 1, j]=WTYSUM_ver[i + 1, j]+wtmi1j
                # i;j+1 Node
                RHOYSUM_ver[i, j + 1]=RHOYSUM_ver[i, j + 1]+property[1, m]*wtmij1
                RHOFYSUM_ver[i, j + 1]=RHOFYSUM_ver[i, j + 1]+property[2, m]*wtmij1
                KYSUM_ver[i, j + 1]=KYSUM_ver[i, j + 1]+property[3, m]*wtmij1
                PHIYSUM_ver[i, j + 1]=PHIYSUM_ver[i, j + 1]+property[4, m]*wtmij1
                RYSUM_ver[i, j + 1]=RYSUM_ver[i, j + 1]+property[5, m]*wtmij1
                WTYSUM_ver[i, j + 1]=WTYSUM_ver[i, j + 1]+wtmij1
                # i+1;j+1 Node
                RHOYSUM_ver[i + 1, j + 1]=RHOYSUM_ver[i + 1, j + 1]+property[1, m]*wtmi1j1
                RHOFYSUM_ver[i + 1, j + 1]=RHOFYSUM_ver[i + 1, j + 1]+property[2, m]*wtmi1j1
                KYSUM_ver[i + 1, j + 1]=KYSUM_ver[i + 1, j + 1]+property[3, m]*wtmi1j1
                PHIYSUM_ver[i + 1, j + 1]=PHIYSUM_ver[i + 1, j + 1]+property[4, m]*wtmi1j1
                RYSUM_ver[i + 1, j + 1]=RYSUM_ver[i + 1, j + 1]+property[5, m]*wtmi1j1
                WTYSUM_ver[i + 1, j + 1]=WTYSUM_ver[i + 1, j + 1]+wtmi1j1
            end
            for j in 1:1:Nx1
                for i in 1:1:Ny1
                    if WTYSUM_ver[i, j]>0
                        RHOY_ver[i, j]=RHOYSUM_ver[i, j]/WTYSUM_ver[i, j]
                        RHOFY_ver[i, j]=RHOFYSUM_ver[i, j]/WTYSUM_ver[i, j]
                        KY_ver[i, j]=KYSUM_ver[i, j]/WTYSUM_ver[i, j]
                        PHIY_ver[i, j]=PHIYSUM_ver[i, j]/WTYSUM_ver[i, j]
                        RY_ver[i, j]=RYSUM_ver[i, j]/WTYSUM_ver[i, j]
                    end
                end
            end
            #test
            for j in 1:1:Nx1, i in 1:1:Ny1
                @test RHOY[i, j] ≈ RHOY_ver[i, j] rtol=1e-9
                @test RHOFY[i, j] ≈ RHOFY_ver[i, j] rtol=1e-9
                @test KY[i, j] ≈ KY_ver[i, j] rtol=1e-9
                @test PHIY[i, j] ≈ PHIY_ver[i, j] rtol=1e-9
                @test RY[i, j] ≈ RY_ver[i, j] rtol=1e-9
            end
        end # testset "compute_vy_node_properties!()"

        @testset "compute_p_node_properties!()" begin
            jmin, jmax = jmin_p, jmax_p
            imin, imax = imin_p, imax_p
            RHOSUM = zeros(Ny1, Nx1)
            RHOCPSUM = zeros(Ny1, Nx1)
            ALPHASUM = zeros(Ny1, Nx1)
            ALPHAFSUM = zeros(Ny1, Nx1)
            HRSUM = zeros(Ny1, Nx1)
            GGGPSUM = zeros(Ny1, Nx1)
            SXXSUM = zeros(Ny1, Nx1)
            TKSUM = zeros(Ny1, Nx1)
            PHISUM = zeros(Ny1, Nx1)
            WTPSUM = zeros(Ny1, Nx1)
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
            xm = rand(rgen, (-xp[1]):0.1:(xp[end] + dx), marknum)
            ym = rand(rgen, (-yp[1]):0.1:(yp[end] + dy), marknum)
            property = rand(rgen, 9, marknum)*1e6
            # calculate grid properties
            for m in 1:1:marknum
                i, j, weights = Erebus.fix_weights(
                    xm[m], ym[m], xp, yp, dx, dy, jmin, jmax, imin, imax
                )
                Erebus.interpolate_add_to_grid!(i, j, weights, property[1, m], RHOSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, property[2, m], RHOCPSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, property[3, m], ALPHASUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, property[4, m], ALPHAFSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, property[5, m], HRSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, inv(property[6, m]), GGGPSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, property[7, m], SXXSUM)
                Erebus.interpolate_add_to_grid!(
                    i, j, weights, property[2, m] * property[8, m], TKSUM
                )
                Erebus.interpolate_add_to_grid!(i, j, weights, property[9, m], PHISUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, 1.0, WTPSUM)
            end
            Erebus.compute_p_node_properties!(
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
                BETTAPHI,
            )
            # verification properties, from HTM-planetary.m, line 538ff, 648ff
            for m in 1:1:marknum
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
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Update properties
                # i;j Node
                GGGPSUM_ver[i, j]=GGGPSUM_ver[i, j]+1/property[6, m]*wtmij
                SXXSUM_ver[i, j]=SXXSUM_ver[i, j]+property[7, m]*wtmij
                RHOSUM_ver[i, j]=RHOSUM_ver[i, j]+property[1, m]*wtmij
                RHOCPSUM_ver[i, j]=RHOCPSUM_ver[i, j]+property[2, m]*wtmij
                ALPHASUM_ver[i, j]=ALPHASUM_ver[i, j]+property[3, m]*wtmij
                ALPHAFSUM_ver[i, j]=ALPHAFSUM_ver[i, j]+property[4, m]*wtmij
                HRSUM_ver[i, j]=HRSUM_ver[i, j]+property[5, m]*wtmij
                TKSUM_ver[i, j]=TKSUM_ver[i, j]+property[8, m] * property[2, m] * wtmij
                PHISUM_ver[i, j]=PHISUM_ver[i, j]+property[9, m]*wtmij
                WTPSUM_ver[i, j]=WTPSUM_ver[i, j]+wtmij
                # i+1;j Node
                GGGPSUM_ver[i + 1, j]=GGGPSUM_ver[i + 1, j]+1/property[6, m]*wtmi1j
                SXXSUM_ver[i + 1, j]=SXXSUM_ver[i + 1, j]+property[7, m]*wtmi1j
                RHOSUM_ver[i + 1, j]=RHOSUM_ver[i + 1, j]+property[1, m]*wtmi1j
                RHOCPSUM_ver[i + 1, j]=RHOCPSUM_ver[i + 1, j]+property[2, m]*wtmi1j
                ALPHASUM_ver[i + 1, j]=ALPHASUM_ver[i + 1, j]+property[3, m]*wtmi1j
                ALPHAFSUM_ver[i + 1, j]=ALPHAFSUM_ver[i + 1, j]+property[4, m]*wtmi1j
                HRSUM_ver[i + 1, j]=HRSUM_ver[i + 1, j]+property[5, m]*wtmi1j
                TKSUM_ver[i + 1, j]=TKSUM_ver[i + 1, j]+property[8, m] *
                                                        property[2, m] *
                                                        wtmi1j
                PHISUM_ver[i + 1, j]=PHISUM_ver[i + 1, j]+property[9, m]*wtmi1j
                WTPSUM_ver[i + 1, j]=WTPSUM_ver[i + 1, j]+wtmi1j
                # i;j+1 Node
                GGGPSUM_ver[i, j + 1]=GGGPSUM_ver[i, j + 1]+1/property[6, m]*wtmij1
                SXXSUM_ver[i, j + 1]=SXXSUM_ver[i, j + 1]+property[7, m]*wtmij1
                RHOSUM_ver[i, j + 1]=RHOSUM_ver[i, j + 1]+property[1, m]*wtmij1
                RHOCPSUM_ver[i, j + 1]=RHOCPSUM_ver[i, j + 1]+property[2, m]*wtmij1
                ALPHASUM_ver[i, j + 1]=ALPHASUM_ver[i, j + 1]+property[3, m]*wtmij1
                ALPHAFSUM_ver[i, j + 1]=ALPHAFSUM_ver[i, j + 1]+property[4, m]*wtmij1
                HRSUM_ver[i, j + 1]=HRSUM_ver[i, j + 1]+property[5, m]*wtmij1
                TKSUM_ver[i, j + 1]=TKSUM_ver[i, j + 1]+property[8, m] *
                                                        property[2, m] *
                                                        wtmij1
                PHISUM_ver[i, j + 1]=PHISUM_ver[i, j + 1]+property[9, m]*wtmij1
                WTPSUM_ver[i, j + 1]=WTPSUM_ver[i, j + 1]+wtmij1
                # i+1;j+1 Node
                GGGPSUM_ver[i + 1, j + 1]=GGGPSUM_ver[i + 1, j + 1]+1/property[6, m] *
                                                                    wtmi1j1
                SXXSUM_ver[i + 1, j + 1]=SXXSUM_ver[i + 1, j + 1]+property[7, m]*wtmi1j1
                RHOSUM_ver[i + 1, j + 1]=RHOSUM_ver[i + 1, j + 1]+property[1, m]*wtmi1j1
                RHOCPSUM_ver[i + 1, j + 1]=RHOCPSUM_ver[i + 1, j + 1] +
                                           property[2, m]*wtmi1j1
                ALPHASUM_ver[i + 1, j + 1]=ALPHASUM_ver[i + 1, j + 1] +
                                           property[3, m]*wtmi1j1
                ALPHAFSUM_ver[i + 1, j + 1]=ALPHAFSUM_ver[i + 1, j + 1] +
                                            property[4, m]*wtmi1j1
                HRSUM_ver[i + 1, j + 1]=HRSUM_ver[i + 1, j + 1]+property[5, m]*wtmi1j1
                TKSUM_ver[i + 1, j + 1]=TKSUM_ver[i + 1, j + 1]+property[8, m] *
                                                                property[2, m] *
                                                                wtmi1j1
                PHISUM_ver[i + 1, j + 1]=PHISUM_ver[i + 1, j + 1]+property[9, m]*wtmi1j1
                WTPSUM_ver[i + 1, j + 1]=WTPSUM_ver[i + 1, j + 1]+wtmi1j1
            end
            for j in 1:1:Nx1
                for i in 1:1:Ny1
                    if WTPSUM_ver[i, j]>0
                        GGGP_ver[i, j]=1/(GGGPSUM_ver[i, j]/WTPSUM_ver[i, j])
                        SXX0_ver[i, j]=SXXSUM_ver[i, j]/WTPSUM_ver[i, j]
                        RHO_ver[i, j]=RHOSUM_ver[i, j]/WTPSUM_ver[i, j]
                        RHOCP_ver[i, j]=RHOCPSUM_ver[i, j]/WTPSUM_ver[i, j]
                        ALPHA_ver[i, j]=ALPHASUM_ver[i, j]/WTPSUM_ver[i, j]
                        ALPHAF_ver[i, j]=ALPHAFSUM_ver[i, j]/WTPSUM_ver[i, j]
                        HR_ver[i, j]=HRSUM_ver[i, j]/WTPSUM_ver[i, j]
                        PHI_ver[i, j]=PHISUM_ver[i, j]/WTPSUM_ver[i, j]
                        BETTAPHI_ver[i, j]=1/GGGP_ver[i, j]*PHI_ver[i, j]
                        tk1_ver[i, j]=TKSUM_ver[i, j]/RHOCPSUM_ver[i, j]
                    end
                end
            end
            # test
            for j in 1:1:Nx1, i in 1:1:Ny1
                @test RHO[i, j] ≈ RHO_ver[i, j] rtol=1e-9
                @test RHOCP[i, j] ≈ RHOCP_ver[i, j] rtol=1e-9
                @test ALPHA[i, j] ≈ ALPHA_ver[i, j] rtol=1e-9
                @test ALPHAF[i, j] ≈ ALPHAF_ver[i, j] rtol=1e-9
                @test HR[i, j] ≈ HR_ver[i, j] rtol=1e-9
                @test GGGP[i, j] ≈ GGGP_ver[i, j] rtol=1e-9
                @test SXX0[i, j] ≈ SXX0_ver[i, j] rtol=1e-9
                @test tk1[i, j] ≈ tk1_ver[i, j] rtol=1e-9
                @test PHI[i, j] ≈ PHI_ver[i, j] rtol=1e-9
                @test BETTAPHI[i, j] ≈ BETTAPHI_ver[i, j] rtol=1e-9
            end
        end # testset "compute_p_node_properties!()"

        @testset "compute_thermodynamic_xfer!()" begin
            jmin, jmax = jmin_p, jmax_p
            imin, imax = imin_p, imax_p
            DMPSUM = zeros(Ny1, Nx1)
            DHPSUM = zeros(Ny1, Nx1)
            WTPSUM = zeros(Ny1, Nx1)
            DMP = rand(rgen, Ny1, Nx1)
            DHP = rand(rgen, Ny1, Nx1)
            DMPSUM_ver = zeros(Ny1, Nx1)
            DHPSUM_ver = zeros(Ny1, Nx1)
            WTPSUM_ver = zeros(Ny1, Nx1)
            DMP_ver = copy(DMP)
            DHP_ver = copy(DHP)
            # simulate markers
            xm = rand(rgen, (-xp[1]):0.1:(xp[end] + dx), marknum)
            ym = rand(rgen, (-yp[1]):0.1:(yp[end] + dy), marknum)
            property = rand(rgen, 2, marknum)*1e3
            # calculate grid properties
            for m in 1:1:marknum
                i, j, weights = Erebus.fix_weights(
                    xm[m], ym[m], xp, yp, dx, dy, jmin, jmax, imin, imax
                )
                Erebus.interpolate_add_to_grid!(i, j, weights, property[1, m], DMPSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, property[2, m], DHPSUM)
                Erebus.interpolate_add_to_grid!(i, j, weights, one(1.0), WTPSUM)
            end
            Erebus.compute_thermodynamic_xfer!(DMPSUM, DHPSUM, WTPSUM, DMP, DHP)
            # verification properties, from i2visHTM_hydration.m, line 566f, 669ff,  
            for m in 1:1:marknum
                DMm = property[1, m]
                DHm = property[2, m]
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
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Interpolation to pressure nodes 
                # Update subgrid diffusion on nodes
                # i,j Node
                DMPSUM_ver[i, j]=DMPSUM_ver[i, j]+DMm*wtmij
                DHPSUM_ver[i, j]=DHPSUM_ver[i, j]+DHm*wtmij
                WTPSUM_ver[i, j]=WTPSUM_ver[i, j]+wtmij
                # i+1,j Node
                DMPSUM_ver[i + 1, j]=DMPSUM_ver[i + 1, j]+DMm*wtmi1j
                DHPSUM_ver[i + 1, j]=DHPSUM_ver[i + 1, j]+DHm*wtmi1j
                WTPSUM_ver[i + 1, j]=WTPSUM_ver[i + 1, j]+wtmi1j
                # i,j+1 Node
                DMPSUM_ver[i, j + 1]=DMPSUM_ver[i, j + 1]+DMm*wtmij1
                DHPSUM_ver[i, j + 1]=DHPSUM_ver[i, j + 1]+DHm*wtmij1
                WTPSUM_ver[i, j + 1]=WTPSUM_ver[i, j + 1]+wtmij1
                # i+1,j+1 Node
                DMPSUM_ver[i + 1, j + 1]=DMPSUM_ver[i + 1, j + 1]+DMm*wtmi1j1
                DHPSUM_ver[i + 1, j + 1]=DHPSUM_ver[i + 1, j + 1]+DHm*wtmi1j1
                WTPSUM_ver[i + 1, j + 1]=WTPSUM_ver[i + 1, j + 1]+wtmi1j1
            end
            # P-nodes
            DMP_ver=zeros(Ny1, Nx1)
            DHP_ver=zeros(Ny1, Nx1)
            for j in 1:1:Nx1
                for i in 1:1:Ny1
                    if WTPSUM_ver[i, j]>0
                        DMP_ver[i, j]=DMPSUM[i, j]/WTPSUM_ver[i, j]
                        DHP_ver[i, j]=DHPSUM[i, j]/WTPSUM_ver[i, j]
                    end
                end
            end
            # test
            for j in 1:1:Nx1, i in 1:1:Ny1
                @test DMP[i, j] ≈ DMP_ver[i, j] rtol=1e-9
                @test DHP[i, j] ≈ DHP_ver[i, j] rtol=1e-9
            end
        end # testset "compute_thermodynamic_xfer!()"

        @testset "molarfraction_marker_to_p_nodes! & compute_thermodynamic_xfer!()" begin
            jmin, jmax = jmin_p, jmax_p
            imin, imax = imin_p, imax_p
            XWSSUM = zeros(Ny1, Nx1)
            WTPSUM = zeros(Ny1, Nx1)
            XWS = rand(rgen, Ny1, Nx1)
            XWSSUM_ver = zeros(Ny1, Nx1)
            WTPSUM_ver = zeros(Ny1, Nx1)
            XWS_ver = copy(XWS)
            # simulate markers
            xm = rand(rgen, (-xp[1]):0.1:(xp[end] + dx), marknum)
            ym = rand(rgen, (-yp[1]):0.1:(yp[end] + dy), marknum)
            XWsolidm0 = rand(rgen, marknum)
            # calculate grid properties
            Erebus.update_p_nodes_melt_composition!(
                xm, ym, XWsolidm0, XWS, XWSSUM, WTPSUM, marknum
            )
            # verification properties, from i2visHTM_hydration.m, line 1547ff
            for m in 1:1:marknum
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
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Update fluid composition on nodes
                # i,j Node
                XWSSUM_ver[i, j]=XWSSUM_ver[i, j]+XWsolidm0[m]*wtmij
                WTPSUM_ver[i, j]=WTPSUM_ver[i, j]+wtmij
                # i+1,j Node
                XWSSUM_ver[i + 1, j]=XWSSUM_ver[i + 1, j]+XWsolidm0[m]*wtmi1j
                WTPSUM_ver[i + 1, j]=WTPSUM_ver[i + 1, j]+wtmi1j
                # i,j+1 Node
                XWSSUM_ver[i, j + 1]=XWSSUM_ver[i, j + 1]+XWsolidm0[m]*wtmij1
                WTPSUM_ver[i, j + 1]=WTPSUM_ver[i, j + 1]+wtmij1
                # i+1,j+1 Node
                XWSSUM_ver[i + 1, j + 1]=XWSSUM_ver[i + 1, j + 1]+XWsolidm0[m]*wtmi1j1
                WTPSUM_ver[i + 1, j + 1]=WTPSUM_ver[i + 1, j + 1]+wtmi1j1
            end
            # P-nodes
            XWS_ver=zeros(Ny1, Nx1)
            for j in 1:1:Nx1
                for i in 1:1:Ny1
                    if WTPSUM_ver[i, j]>0
                        XWS_ver[i, j]=XWSSUM_ver[i, j]/WTPSUM_ver[i, j]
                    end
                end
            end
            # test
            for j in 1:1:Nx1, i in 1:1:Ny1
                @test XWS[i, j] ≈ XWS_ver[i, j] rtol=1e-9
            end
        end # testset "compute_molarfraction!()"
    end # testset "compute node properties" 

    @testset "apply_subgrid_stress_diffusion!()" begin
        marknum = start_marknum
        dtm = dt_longest
        etam = @SVector ones(3)
        # simulate markers
        xm = rand(rgen, (-dx):0.1:(x[end] + dx), marknum)
        ym = rand(rgen, (-dy):0.1:(y[end] + dy), marknum)
        tm = rand(rgen, 1:3, marknum)
        gggm = rand(rgen, 3)
        inv_gggtotalm = inv.([gggm[tm[m]] for m in 1:marknum])
        sxxm = rand(rgen, marknum)
        sxym = rand(rgen, marknum)
        SXX0 = rand(rgen, Ny1, Nx1)
        SXY0 = rand(rgen, Ny, Nx)
        DSXX = rand(rgen, Ny1, Nx1)
        DSXY = rand(rgen, Ny, Nx)
        SXXSUM = rand(rgen, Ny1, Nx1, Base.Threads.nthreads())
        SXYSUM = rand(rgen, Ny, Nx, Base.Threads.nthreads())
        WTPSUM = rand(rgen, Ny1, Nx1, Base.Threads.nthreads())
        WTSUM = rand(rgen, Ny, Nx, Base.Threads.nthreads())
        sxxm_ver = deepcopy(sxxm)
        sxym_ver = deepcopy(sxym)
        DSXX_ver = deepcopy(DSXX)
        DSXY_ver = deepcopy(DSXY)
        # apply subgrid stress diffusion
        Erebus.apply_subgrid_stress_diffusion!(
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
        )
        # verification, from HTM-planetary.m, line 1374ff
        # Apply subgrid stress diffusion to markers
        if dsubgrids>0
            SXYSUM_ver = zeros(Ny, Nx)
            WTSUM_ver = zeros(Ny, Nx)
            SXXSUM_ver = zeros(Ny1, Nx1)
            WTPSUM_ver = zeros(Ny1, Nx1)
            for m in 1:1:marknum
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
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Compute marker-node SIGMA'xx difference
                dsxxm0=sxxm_ver[m]-(
                    SXX0[i, j]*wtmij+SXX0[i + 1, j]*wtmi1j+SXX0[i, j + 1]*wtmij1+SXX0[
                        i + 1, j + 1
                    ]*wtmi1j1
                )
                # Relax stress difference
                dsxxm1=dsxxm0*exp(-dsubgrids*dtm/(etam[tm[m]]/gggm[tm[m]]))
                # Correct marker stress
                ddsxxm_ver=dsxxm1-dsxxm0
                sxxm_ver[m]=sxxm_ver[m]+ddsxxm_ver
                # Update subgrid diffusion on nodes
                # i;j Node
                SXXSUM_ver[i, j]=SXXSUM_ver[i, j]+ddsxxm_ver*wtmij
                WTPSUM_ver[i, j]=WTPSUM_ver[i, j]+wtmij
                # i+1;j Node
                SXXSUM_ver[i + 1, j]=SXXSUM_ver[i + 1, j]+ddsxxm_ver*wtmi1j
                WTPSUM_ver[i + 1, j]=WTPSUM_ver[i + 1, j]+wtmi1j
                # i;j+1 Node
                SXXSUM_ver[i, j + 1]=SXXSUM_ver[i, j + 1]+ddsxxm_ver*wtmij1
                WTPSUM_ver[i, j + 1]=WTPSUM_ver[i, j + 1]+wtmij1
                # i+1;j+1 Node
                SXXSUM_ver[i + 1, j + 1]=SXXSUM_ver[i + 1, j + 1]+ddsxxm_ver*wtmi1j1
                WTPSUM_ver[i + 1, j + 1]=WTPSUM_ver[i + 1, j + 1]+wtmi1j1
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
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Compute marker-node SIGMAxy difference
                dsxym0=sxym_ver[m]-(
                    SXY0[i, j]*wtmij+SXY0[i + 1, j]*wtmi1j+SXY0[i, j + 1]*wtmij1+SXY0[
                        i + 1, j + 1
                    ]*wtmi1j1
                )
                # Relax stress difference
                dsxym1=dsxym0*exp(-dsubgrids*dtm/(etam[tm[m]]/gggm[tm[m]]))
                # Correct marker stress
                ddsxym_ver=dsxym1-dsxym0
                sxym_ver[m]=sxym_ver[m]+ddsxym_ver
                # Update subgrid diffusion on nodes
                # i;j Node
                SXYSUM_ver[i, j]=SXYSUM_ver[i, j]+ddsxym_ver*wtmij
                WTSUM_ver[i, j]=WTSUM_ver[i, j]+wtmij
                # i+1;j Node
                SXYSUM_ver[i + 1, j]=SXYSUM_ver[i + 1, j]+ddsxym_ver*wtmi1j
                WTSUM_ver[i + 1, j]=WTSUM_ver[i + 1, j]+wtmi1j
                # i;j+1 Node
                SXYSUM_ver[i, j + 1]=SXYSUM_ver[i, j + 1]+ddsxym_ver*wtmij1
                WTSUM_ver[i, j + 1]=WTSUM_ver[i, j + 1]+wtmij1
                # i+1;j+1 Node
                SXYSUM_ver[i + 1, j + 1]=SXYSUM_ver[i + 1, j + 1]+ddsxym_ver*wtmi1j1
                WTSUM_ver[i + 1, j + 1]=WTSUM_ver[i + 1, j + 1]+wtmi1j1
            end
            # Compute DSXXsubgrid_ver
            DSXXsubgrid_ver = zeros(Ny1, Nx1)
            # P-nodes
            for j in 2:1:Nx
                for i in 2:1:Ny
                    if (WTPSUM_ver[i, j]>0)
                        DSXXsubgrid_ver[i, j]=SXXSUM_ver[i, j]/WTPSUM_ver[i, j]
                    end
                end
            end
            # Correct DSXX_ver
            DSXX_ver=DSXX_ver-DSXXsubgrid_ver
            # Compute DSXYsubgrid_ver
            DSXYsubgrid_ver = zeros(Ny, Nx)
            # Basic nodes
            for j in 1:1:Nx
                for i in 1:1:Ny
                    if (WTSUM_ver[i, j]>0)
                        DSXYsubgrid_ver[i, j]=SXYSUM_ver[i, j]/WTSUM_ver[i, j]
                    end
                end
            end
            # Correct DSXY_ver
            DSXY_ver=DSXY_ver-DSXYsubgrid_ver
        end
        # test
        @test sxxm ≈ sxxm_ver rtol=1e-9
        @test sxym ≈ sxym_ver rtol=1e-9
        @test DSXX ≈ DSXX_ver rtol=1e-9
        @test DSXY ≈ DSXY_ver rtol=1e-9
    end # testset "apply_subgrid_stress_diffusion!()"

    @testset "update_marker_stress!()" begin
        marknum = start_marknum
        # simulate markers
        xm = rand(rgen, (-dx):0.1:(x[end] + dx), marknum)
        ym = rand(rgen, (-dy):0.1:(y[end] + dy), marknum)
        sxxm = rand(rgen, marknum)
        sxym = rand(rgen, marknum)
        DSXX = rand(rgen, Ny1, Nx1)
        DSXY = rand(rgen, Ny, Nx)
        sxxm_ver = deepcopy(sxxm)
        sxym_ver = deepcopy(sxym)
        # update marker stress
        Erebus.update_marker_stress!(xm, ym, sxxm, sxym, DSXX, DSXY, marknum)
        # verification, from HTM-planetary.m, line 1495ff
        for m in 1:1:marknum
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
            wtmi1j=(1-dxmj/dx)*(dymi/dy)
            wtmij1=(dxmj/dx)*(1-dymi/dy)
            wtmi1j1=(dxmj/dx)*(dymi/dy)
            # Update marker by SIGMA'xx change 
            sxxm_ver[m]=sxxm_ver[m]+(
                DSXX[i, j]*wtmij+DSXX[i + 1, j]*wtmi1j+DSXX[i, j + 1]*wtmij1+DSXX[
                    i + 1, j + 1
                ]*wtmi1j1
            )

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
            wtmi1j=(1-dxmj/dx)*(dymi/dy)
            wtmij1=(dxmj/dx)*(1-dymi/dy)
            wtmi1j1=(dxmj/dx)*(dymi/dy)
            # Update marker by SIGMA'xx change 
            sxym_ver[m]=sxym_ver[m]+(
                DSXY[i, j]*wtmij+DSXY[i + 1, j]*wtmi1j+DSXY[i, j + 1]*wtmij1+DSXY[
                    i + 1, j + 1
                ]*wtmi1j1
            )
        end
        # test
        @test sxxm ≈ sxxm_ver rtol=1e-9
        @test sxym ≈ sxym_ver rtol=1e-9
    end # testset "update_marker_stress!()"

    @testset "apply_subgrid_temperature_diffusion!()" begin
        marknum = start_marknum
        dtm = dt_longest
        mode = 2
        # simulate markers
        xm = rand(rgen, (-dx):0.1:(x[end] + dx), marknum)
        ym = rand(rgen, (-dy):0.1:(y[end] + dy), marknum)
        tm = rand(rgen, 1:3, marknum)
        tkm = rand(rgen, marknum)
        phim = rand(rgen, marknum)
        tk1 = rand(rgen, Ny1, Nx1)
        DT = rand(rgen, Ny1, Nx1)
        TKSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
        RHOCPSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
        tkm_ver = deepcopy(tkm)
        DT_ver = deepcopy(DT)
        # apply subgrid stress diffusion
        Erebus.apply_subgrid_temperature_diffusion!(
            xm, ym, tm, tkm, phim, tk1, DT, TKSUM, RHOCPSUM, dtm, marknum, mode
        )
        # verification, from HTM-planetary.m, line 1731ff
        # Apply subgrid stress diffusion to markers
        if dsubgridt>0
            TKSUM_ver = zeros(Ny1, Nx1)
            RHOCPSUM_ver = zeros(Ny1, Nx1)
            for m in 1:1:marknum
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
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Compute marker-node T difference
                dtkm0=tkm_ver[m]-(
                    tk1[i, j]*wtmij+tk1[i + 1, j]*wtmi1j+tk1[i, j + 1]*wtmij1+tk1[
                        i + 1, j + 1
                    ]*wtmi1j1
                )
                # Compute marker parameters
                if tm[m]<3
                    # Rocks
                    rhocptotalm=rhocpsolidm[tm[m]]*(1-phim[m])+rhocpfluidm[tm[m]]*phim[m]
                    ktotalm=(
                        ksolidm[tm[m]]*kfluidm[tm[m]]/2+(
                            (ksolidm[tm[m]]*(3*phim[m]-2) + kfluidm[tm[m]]*(1-3*phim[m]))^2
                        )/16
                    )^0.5-(ksolidm[tm[m]]*(3*phim[m]-2) + kfluidm[tm[m]]*(1-3*phim[m]))/4
                else
                    # Sticky air
                    rhocptotalm=rhocpsolidm[tm[m]]
                    ktotalm=ksolidm[tm[m]]
                end    # Relax temperature difference
                dtkm1=dtkm0*exp(-dsubgridt*ktotalm*dtm/rhocptotalm*(2/dx^2+2/dy^2))
                # Correct marker temperature
                ddtkm=dtkm1-dtkm0
                tkm_ver[m]=tkm_ver[m]+ddtkm
                # Update subgrid diffusion on nodes
                # i;j Node
                TKSUM_ver[i, j]=TKSUM_ver[i, j]+ddtkm*rhocptotalm*wtmij
                RHOCPSUM_ver[i, j]=RHOCPSUM_ver[i, j]+rhocptotalm*wtmij
                # i+1;j Node
                TKSUM_ver[i + 1, j]=TKSUM_ver[i + 1, j]+ddtkm*rhocptotalm*wtmi1j
                RHOCPSUM_ver[i + 1, j]=RHOCPSUM_ver[i + 1, j]+rhocptotalm*wtmi1j
                # i;j+1 Node
                TKSUM_ver[i, j + 1]=TKSUM_ver[i, j + 1]+ddtkm*rhocptotalm*wtmij1
                RHOCPSUM_ver[i, j + 1]=RHOCPSUM_ver[i, j + 1]+rhocptotalm*wtmij1
                # i+1;j+1 Node
                TKSUM_ver[i + 1, j + 1]=TKSUM_ver[i + 1, j + 1]+ddtkm*rhocptotalm*wtmi1j1
                RHOCPSUM_ver[i + 1, j + 1]=RHOCPSUM_ver[i + 1, j + 1]+rhocptotalm*wtmi1j1
            end
            # Compute DTsubgrid
            DTsubgrid = zeros(Ny1, Nx1)
            # P-nodes
            for j in 1:1:Nx1
                for i in 1:1:Ny1
                    if RHOCPSUM_ver[i, j]>0
                        DTsubgrid[i, j]=TKSUM_ver[i, j]/RHOCPSUM_ver[i, j]
                    end
                end
            end
            # Correct DT
            DT_ver=DT_ver-DTsubgrid
        end
        # test
        @test tkm ≈ tkm_ver rtol=1e-9
        @test DT ≈ DT_ver rtol=1e-9
    end # testset "apply_subgrid_temperature_diffusion!()"

    @testset "update_marker_temperature!()" begin
        marknum = start_marknum
        # simulate markers
        xm = rand(rgen, (-dx):0.1:(x[end] + dx), marknum)
        ym = rand(rgen, (-dy):0.1:(y[end] + dy), marknum)
        DT = rand(rgen, Ny1, Nx1)
        tk2 = rand(rgen, Ny1, Nx1)
        for timestep in 1:2
            tkm = rand(rgen, marknum)
            tkm_ver = deepcopy(tkm)
            # update marker temperature
            Erebus.update_marker_temperature!(xm, ym, tkm, DT, tk2, timestep, marknum)
            # verification, from HTM-planetary.m, line 1805ff 
            for m in 1:1:marknum
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
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Update properties
                tkm_ver[m]=tkm_ver[m]+DT[i, j]*wtmij+DT[i + 1, j]*wtmi1j+DT[i, j + 1]*wtmij1+DT[
                    i + 1, j + 1
                ]*wtmi1j1
                # Interpolate tk2 at 1st timestep
                if timestep==1
                    tkm_ver[m]=tk2[i, j]*wtmij+tk2[i + 1, j]*wtmi1j+tk2[i, j + 1]*wtmij1+tk2[
                        i + 1, j + 1
                    ]*wtmi1j1
                end
            end
            # test
            @test tkm ≈ tkm_ver rtol=1e-9
        end # for timestep = 1:2
    end # testset "update_marker_temperature!()"

    @testset "update_marker_porosity!()" begin
        marknum = start_marknum
        dtm = dt_longest
        # simulate markers
        xm = rand(rgen, (-dx):0.1:(x[end] + dx), marknum)
        ym = rand(rgen, (-dy):0.1:(y[end] + dy), marknum)
        tm = rand(rgen, 1:3, marknum)
        APHI = rand(rgen, Ny1, Nx1)
        phim = rand(rgen, marknum)
        phim_ver = deepcopy(phim)
        # update marker porosity
        Erebus.update_marker_porosity!(xm, ym, tm, phim, APHI, dtm, marknum)
        # verification, from HTM-planetary.m, line 1805ff 
        for m in 1:1:marknum
            if tm[m]<3
                # Interpolate APHI
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
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Compute Dln[(1-PHI)/PHI]/Dt
                aphim=APHI[i, j]*wtmij+APHI[i + 1, j]*wtmi1j+APHI[i, j + 1]*wtmij1+APHI[
                    i + 1, j + 1
                ]*wtmi1j1
                # Change Porosity
                phim_ver[m]=phim_ver[m]/((1-phim_ver[m])*exp(aphim*dtm)+phim_ver[m])
                if phim_ver[m]<phimin
                    phim_ver[m]=phimin
                elseif phim_ver[m]>phimax
                    phim_ver[m]=phimax
                end
            end
        end
        # test
        @test phim ≈ phim_ver rtol=1e-9
    end # testset "update_marker_porosity!()"

    @testset "add_vrk4()" begin
        vrk4 = @SVector zeros(4)
        v = rand(rgen)
        for rk in 1:4
            vrk4 = Erebus.add_vrk4(vrk4, v, rk)
        end
        @test vrk4 == @SVector [v, v, v, v]
    end # testset "add_vrk4()"

    @testset "compute_velocities!()" begin
        vx = rand(rgen, Ny1, Nx1)
        vy = rand(rgen, Ny1, Nx1)
        vxf = rand(rgen, Ny1, Nx1)
        vyf = rand(rgen, Ny1, Nx1)
        vxp = zeros(Ny1, Nx1)
        vyp = zeros(Ny1, Nx1)
        vxpf = zeros(Ny1, Nx1)
        vypf = zeros(Ny1, Nx1)
        vxp_ver = zero(vxp)
        vyp_ver = zero(vyp)
        vxpf_ver = zero(vxpf)
        vypf_ver = zero(vypf)
        # compute velocities
        Erebus.compute_velocities!(vx, vy, vxf, vyf, vxp, vyp, vxpf, vypf)
        # verification, from HTM-planetary.m, line 1879ff:
        # Compute fluid velocity in pressure nodes
        # vxpf
        for j in 2:1:Nx
            for i in 2:1:Ny
                vxpf_ver[i, j]=(vxf[i, j]+vxf[i, j - 1])/2
            end
        end
        # Apply BC
        # Top
        @. vxpf_ver[1, 2:(Nx - 1)]=-bcftop*vxpf_ver[2, 2:(Nx - 1)]
        # Bottom
        @. vxpf_ver[Ny1, 2:(Nx - 1)]=-bcfbottom*vxpf_ver[Ny, 2:(Nx - 1)]
        # Left
        @. vxpf_ver[:, 1]=2*vxleft-vxpf_ver[:, 2]
        # Right
        @. vxpf_ver[:, Nx1]=2*vxright-vxpf_ver[:, Nx]
        # vypf
        for j in 2:1:Nx
            for i in 2:1:Ny
                vypf_ver[i, j]=(vyf[i, j]+vyf[i - 1, j])/2
            end
        end
        # Apply BC
        # Left
        @. vypf_ver[2:(Ny - 1), 1]=-bcfleft*vypf_ver[2:(Ny - 1), 2]
        # Right
        @. vypf_ver[2:(Ny - 1), Nx1]=-bcfright*vypf_ver[2:(Ny - 1), Nx] # Free slip    
        # Top
        @. vypf_ver[1, :]=2*vytop-vypf_ver[2, :]
        # Bottom
        @. vypf_ver[Ny1, :]=2*vybottom-vypf_ver[Ny, :]
        # Compute velocity in pressure nodes
        # vx
        for j in 2:1:Nx
            for i in 2:1:Ny
                vxp_ver[i, j]=(vx[i, j]+vx[i, j - 1])/2
            end
        end
        # Apply BC
        # Top
        @. vxp_ver[1, 2:(Nx - 1)]=-bctop*vxp_ver[2, 2:(Nx - 1)]
        # Bottom
        @. vxp_ver[Ny1, 2:(Nx - 1)]=-bcbottom*vxp_ver[Ny, 2:(Nx - 1)]
        # Left
        @. vxp_ver[:, 1]=2*vxleft-vxp_ver[:, 2]
        # Right
        @. vxp_ver[:, Nx1]=2*vxright-vxp_ver[:, Nx]
        # vy
        for j in 2:1:Nx
            for i in 2:1:Ny
                vyp_ver[i, j]=(vy[i, j]+vy[i - 1, j])/2
            end
        end
        # Apply BC
        # Left
        @. vyp_ver[2:(Ny - 1), 1]=-bcleft*vyp_ver[2:(Ny - 1), 2]
        # Right
        @. vyp_ver[2:(Ny - 1), Nx1]=-bcright*vyp_ver[2:(Ny - 1), Nx] # Free slip    
        # Top
        @. vyp_ver[1, :]=2*vytop-vyp_ver[2, :]
        # Bottom
        @. vyp_ver[Ny1, :]=2*vybottom-vyp_ver[Ny, :]
        # test
        @test vxp ≈ vxp_ver rtol=1e-9
        @test vyp ≈ vyp_ver rtol=1e-9
        @test vxpf ≈ vxpf_ver rtol=1e-9
        @test vypf ≈ vypf_ver rtol=1e-9
    end # testset "compute_velocities!()"

    @testset "compute_rotation_rate!()" begin
        vx = rand(rgen, Ny1, Nx1)
        vy = rand(rgen, Ny1, Nx1)
        wyx = zeros(Ny, Nx)
        wyx_ver = zero(wyx)
        # compute rotation rate
        Erebus.compute_rotation_rate!(vx, vy, wyx)
        # verification, from HTM-planetary.m, line 1942ff:
        # Compute rotation rate wyx=1/2[dVy/dx-dVx/dy] for basic nodes
        for i in 1:1:Ny
            for j in 1:1:Nx
                wyx_ver[i, j]=0.5*((vy[i, j + 1]-vy[i, j])/dx-(vx[i + 1, j]-vx[i, j])/dy)
            end
        end
        # test
        @test wyx ≈ wyx_ver rtol=1e-9
    end # testset "compute_rotation_rate!()"

    @testset "move_markers_rk4!()" begin
        dtm = 0.9
        marknum = start_marknum + 10_000
        mode = marker_property_mode
        xm = zeros(marknum)
        ym = zeros(marknum)
        for jm in 1:1:Nxm, im in 1:1:Nym
            # calculate marker counter
            m = (jm-1) * Nym + im
            # define marker coordinates
            xm[m] = dxm/2 + (jm-1) * dxm
            ym[m] = dym/2 + (im-1) * dym
        end
        tm = rand(rgen, 1:3, marknum)
        tkm = rand(rgen, 263:0.1:310, marknum)
        phim = rand(phimin:1e-4:phimax, marknum)
        sxym = rand(rgen, marknum)
        sxxm = rand(rgen, marknum)
        vx = rand(rgen, Ny1, Nx1)
        vy = rand(rgen, Ny1, Nx1)
        vxf = rand(rgen, Ny1, Nx1)
        vyf = rand(rgen, Ny1, Nx1)
        tk2 = rand(rgen, 263:0.1:310, Ny1, Nx1)
        wyx = rand(rgen, Ny, Nx)
        xm_ver = deepcopy(xm)
        ym_ver = deepcopy(ym)
        tkm_ver = deepcopy(tkm)
        sxym_ver = deepcopy(sxym)
        sxxm_ver = deepcopy(sxxm)
        # move markers with RK4
        Erebus.move_markers_rk4!(
            xm,
            ym,
            tm,
            tkm,
            phim,
            sxym,
            sxxm,
            vx,
            vy,
            vxf,
            vyf,
            wyx,
            tk2,
            marknum,
            dtm,
            mode,
        )
        # verification, from HTM-planetary.m, line 1947ff:
        # Move markers with 4th order Runge-Kutta
        vxm = zeros(4, 1)
        vym = zeros(4, 1)
        for m in 1:1:marknum
            # Interpolate solid temperature for the initial marker location
            # Define i;j indexes for the upper left node
            j=trunc(Int, (xm_ver[m]-xp[1])/dx)+1
            i=trunc(Int, (ym_ver[m]-yp[1])/dy)+1
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
            dxmj=xm_ver[m]-xp[j]
            dymi=ym_ver[m]-yp[i]
            # Compute weights
            wtmij=(1-dxmj/dx)*(1-dymi/dy)
            wtmi1j=(1-dxmj/dx)*(dymi/dy)
            wtmij1=(dxmj/dx)*(1-dymi/dy)
            wtmi1j1=(dxmj/dx)*(dymi/dy)
            # Compute Tsolid
            tksm0=tk2[i, j]*wtmij+tk2[i + 1, j]*wtmi1j+tk2[i, j + 1]*wtmij1+tk2[
                i + 1, j + 1
            ]*wtmi1j1
            # Interpolate local rotation rate
            # Define i;j indexes for the upper left node
            j=trunc(Int, (xm_ver[m]-x[1])/dx)+1
            i=trunc(Int, (ym_ver[m]-y[1])/dy)+1
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
            dxmj=xm_ver[m]-x[j]
            dymi=ym_ver[m]-y[i]
            # Compute weights
            wtmij=(1-dxmj/dx)*(1-dymi/dy)
            wtmi1j=(1-dxmj/dx)*(dymi/dy)
            wtmij1=(dxmj/dx)*(1-dymi/dy)
            wtmi1j1=(dxmj/dx)*(dymi/dy)
            # Compute vx velocity
            omegam=wyx[i, j]*wtmij+wyx[i + 1, j]*wtmi1j+wyx[i, j + 1]*wtmij1+wyx[
                i + 1, j + 1
            ]*wtmi1j1
            # Analytical stress rotation using SIGMA"xx=-SIGMA"yy
            THETA=dtm*omegam # Incremental rotation angle()
            sxxmnew=sxxm_ver[m]*cos(THETA)^2-sxxm_ver[m]*sin(THETA)^2-sxym_ver[m]*sin(
                2*THETA
            )
            sxymnew=sxxm_ver[m]*sin(2*THETA)+sxym_ver[m]*cos(2*THETA)
            sxxm_ver[m]=sxxmnew
            sxym_ver[m]=sxymnew

            # Save initial marker coordinates
            xA=xm_ver[m]
            yA=ym_ver[m]
            for rk in 1:1:4
                # Interpolate vx
                # Define i;j indexes for the upper left node
                j=trunc(Int, (xm_ver[m]-xvx[1])/dx)+1
                i=trunc(Int, (ym_ver[m]-yvx[1])/dy)+1
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
                # Compute distances
                dxmj=xm_ver[m]-xvx[j]
                dymi=ym_ver[m]-yvx[i]
                # Compute weights
                # Compute vx velocity for the top & bottom of the cell()
                vxm13=vx[i, j]*(1-dxmj/dx)+vx[i, j + 1]*dxmj/dx
                vxm24=vx[i + 1, j]*(1-dxmj/dx)+vx[i + 1, j + 1]*dxmj/dx
                # Compute correction
                if dxmj/dx>=0.5
                    if j<Nx-1
                        vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(
                            vx[i, j]-2*vx[i, j + 1]+vx[i, j + 2]
                        )
                        vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(
                            vx[i + 1, j]-2*vx[i + 1, j + 1]+vx[i + 1, j + 2]
                        )
                    end
                else
                    if j>1
                        vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(
                            vx[i, j - 1]-2*vx[i, j]+vx[i, j + 1]
                        )
                        vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(
                            vx[i + 1, j - 1]-2*vx[i + 1, j]+vx[i + 1, j + 1]
                        )
                    end
                end
                # Compute vx
                vxm[rk]=(1-dymi/dy)*vxm13+(dymi/dy)*vxm24

                # Interpolate vy
                # Define i;j indexes for the upper left node
                j=trunc(Int, (xm_ver[m]-xvy[1])/dx)+1
                i=trunc(Int, (ym_ver[m]-yvy[1])/dy)+1
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
                dxmj=xm_ver[m]-xvy[j]
                dymi=ym_ver[m]-yvy[i]
                # Compute weights
                # Compute vy velocity for the left & right of the cell()
                vym12=vy[i, j]*(1-dymi/dy)+vy[i + 1, j]*dymi/dy
                vym34=vy[i, j + 1]*(1-dymi/dy)+vy[i + 1, j + 1]*dymi/dy
                # Compute correction
                if dymi/dy>=0.5
                    if i<Ny-1
                        vym12=vym12+1/2*((dymi/dy-0.5)^2)*(
                            vy[i, j]-2*vy[i + 1, j]+vy[i + 2, j]
                        )
                        vym34=vym34+1/2*((dymi/dy-0.5)^2)*(
                            vy[i, j + 1]-2*vy[i + 1, j + 1]+vy[i + 2, j + 1]
                        )
                    end
                else
                    if i>1
                        vym12=vym12+1/2*((dymi/dy-0.5)^2)*(
                            vy[i - 1, j]-2*vy[i, j]+vy[i + 1, j]
                        )
                        vym34=vym34+1/2*((dymi/dy-0.5)^2)*(
                            vy[i - 1, j + 1]-2*vy[i, j + 1]+vy[i + 1, j + 1]
                        )
                    end
                end
                # Compute vy
                vym[rk]=(1-dxmj/dx)*vym12+(dxmj/dx)*vym34

                # Change coordinates to obtain B;C;D points
                if rk==1 || rk==2
                    xm_ver[m]=xA+dtm/2*vxm[rk]
                    ym_ver[m]=yA+dtm/2*vym[rk]
                elseif rk==3
                    xm_ver[m]=xA+dtm*vxm[rk]
                    ym_ver[m]=yA+dtm*vym[rk]
                end
            end
            # Restore initial coordinates
            xm_ver[m]=xA
            ym_ver[m]=yA
            # Compute effective velocity
            vxmeff=1/6*(vxm[1]+2*vxm[2]+2*vxm[3]+vxm[4])
            vymeff=1/6*(vym[1]+2*vym[2]+2*vym[3]+vym[4])
            # Move markers
            xm_ver[m]=xm_ver[m]+dtm*vxmeff
            ym_ver[m]=ym_ver[m]+dtm*vymeff

            # Backtracing markers with fluid velocity
            xcur=xm_ver[m]
            ycur=ym_ver[m]
            xA=xcur
            yA=ycur
            for rk in 1:1:4
                # Interpolate vx
                # Define i;j indexes for the upper left node
                j=trunc(Int, (xcur-xvx[1])/dx)+1
                i=trunc(Int, (ycur-yvx[1])/dy)+1
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
                # Compute distances
                dxmj=xcur-xvx[j]
                dymi=ycur-yvx[i]
                # Compute weights
                # Compute vx velocity for the top & bottom of the cell()
                vxm13=vxf[i, j]*(1-dxmj/dx)+vxf[i, j + 1]*dxmj/dx
                vxm24=vxf[i + 1, j]*(1-dxmj/dx)+vxf[i + 1, j + 1]*dxmj/dx
                # Compute correction
                if dxmj/dx>=0.5
                    if j<Nx-1
                        vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(
                            vxf[i, j]-2*vxf[i, j + 1]+vxf[i, j + 2]
                        )
                        vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(
                            vxf[i + 1, j]-2*vxf[i + 1, j + 1]+vxf[i + 1, j + 2]
                        )
                    end
                else
                    if j>1
                        vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(
                            vxf[i, j - 1]-2*vxf[i, j]+vxf[i, j + 1]
                        )
                        vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(
                            vxf[i + 1, j - 1]-2*vxf[i + 1, j]+vxf[i + 1, j + 1]
                        )
                    end
                end
                # Compute vx
                vxm[rk]=(1-dymi/dy)*vxm13+(dymi/dy)*vxm24

                # Interpolate vy
                # Define i;j indexes for the upper left node
                j=trunc(Int, (xcur-xvy[1])/dx)+1
                i=trunc(Int, (ycur-yvy[1])/dy)+1
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
                dxmj=xcur-xvy[j]
                dymi=ycur-yvy[i]
                # Compute weights
                # Compute vy velocity for the left & right of the cell()
                vym12=vyf[i, j]*(1-dymi/dy)+vyf[i + 1, j]*dymi/dy
                vym34=vyf[i, j + 1]*(1-dymi/dy)+vyf[i + 1, j + 1]*dymi/dy
                # Compute correction
                if dymi/dy>=0.5
                    if i<Ny-1
                        vym12=vym12+1/2*((dymi/dy-0.5)^2)*(
                            vyf[i, j]-2*vyf[i + 1, j]+vyf[i + 2, j]
                        )
                        vym34=vym34+1/2*((dymi/dy-0.5)^2)*(
                            vyf[i, j + 1]-2*vyf[i + 1, j + 1]+vyf[i + 2, j + 1]
                        )
                    end
                else
                    if i>1
                        vym12=vym12+1/2*((dymi/dy-0.5)^2)*(
                            vyf[i - 1, j]-2*vyf[i, j]+vyf[i + 1, j]
                        )
                        vym34=vym34+1/2*((dymi/dy-0.5)^2)*(
                            vyf[i - 1, j + 1]-2*vyf[i, j + 1]+vyf[i + 1, j + 1]
                        )
                    end
                end
                # Compute vy
                vym[rk]=(1-dxmj/dx)*vym12+(dxmj/dx)*vym34

                # Change coordinates to obtain B;C;D points
                if rk==1 || rk==2
                    xcur=xA-dtm/2*vxm[rk]
                    ycur=yA-dtm/2*vym[rk]
                elseif rk==3
                    xcur=xA-dtm*vxm[rk]
                    ycur=yA-dtm*vym[rk]
                end
            end
            # Compute effective velocity
            vxmeff=1/6*(vxm[1]+2*vxm[2]+2*vxm[3]+vxm[4])
            vymeff=1/6*(vym[1]+2*vym[2]+2*vym[3]+vym[4])
            # Trace the node backward
            xcur=xA-dtm*vxmeff
            ycur=yA-dtm*vymeff
            # Interpolate fluid temperature
            # Define i;j indexes for the upper left node
            j=trunc(Int, (xcur-xp[1])/dx)+1
            i=trunc(Int, (ycur-yp[1])/dy)+1
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
            dxmj=xcur-xp[j]
            dymi=ycur-yp[i]
            # Compute weights
            wtmij=(1-dxmj/dx)*(1-dymi/dy)
            wtmi1j=(1-dxmj/dx)*(dymi/dy)
            wtmij1=(dxmj/dx)*(1-dymi/dy)
            wtmi1j1=(dxmj/dx)*(dymi/dy)
            # Compute nodal Tfluid
            tkfm0=tk2[i, j]*wtmij+tk2[i + 1, j]*wtmi1j+tk2[i, j + 1]*wtmij1+tk2[
                i + 1, j + 1
            ]*wtmi1j1
            # Compute Tfluid-Tsolid for the marker
            dtkfsm=tkfm0-tksm0
            # Correct marker temperature
            tkm_ver[m]=(
                (1-phim[m])*tkm_ver[m]*rhocpsolidm[tm[m]] +
                phim[m]*(tkm_ver[m]+dtkfsm)*rhocpfluidm[tm[m]]
            ) / ((1-phim[m])*rhocpsolidm[tm[m]]+phim[m]*rhocpfluidm[tm[m]])
        end  # end of marker loop
        # test

        @testset "xm" begin
            for m in 1:1:marknum
                @test xm[m] ≈ xm_ver[m] rtol=1e-12
            end
        end
        @testset "ym" begin
            for m in 1:1:marknum
                @test ym[m] ≈ ym_ver[m] rtol=1e-12
            end
        end
        @testset "tkm" begin
            for m in 1:1:marknum
                @test tkm[m] ≈ tkm_ver[m] rtol=1e-4
            end
        end
        @testset "sxym" begin
            for m in 1:1:marknum
                @test sxym[m] == sxym_ver[m]
            end
        end
        @testset "sxxm" begin
            for m in 1:1:marknum
                @test sxxm[m] == sxxm_ver[m]
            end
        end
    end # testset "move_markers_rk4!()"

    @testset "backtrace_pressures_rk4!()" begin
        vx = rand(rgen, Ny1, Nx1) * 2e-9 .- 1e-9
        vy = rand(rgen, Ny1, Nx1) * 2e-9 .- 1e-9
        vxf = rand(rgen, Ny1, Nx1) * 2e-9 .- 1e-9
        vyf = rand(rgen, Ny1, Nx1) * 2e-9 .- 1e-9
        pr = rand(rgen, Ny1, Nx1) * 1e6
        # pr = ones(Ny1, Nx1) * 1e6
        pf = rand(rgen, Ny1, Nx1) * 1e6
        # pf = ones(Ny1, Nx1) * 1e6
        ps = rand(rgen, Ny1, Nx1) * 1e6
        # ps = ones(Ny1, Nx1) * 1e6
        pr0 = rand(rgen, Ny1, Nx1) * 1e6
        pf0 = rand(rgen, Ny1, Nx1) * 1e6
        ps0 = rand(rgen, Ny1, Nx1) * 1e6
        pr_ver = deepcopy(pr)
        pf_ver = deepcopy(pf)
        ps_ver = deepcopy(ps)
        pr0_ver = deepcopy(pr0)
        pf0_ver = deepcopy(pf0)
        ps0_ver = deepcopy(ps0)
        dtm = 0.9
        # verification, from madcph.p, line 2363ff
        Erebus.backtrace_pressures_rk4!(pr, pr0, ps, ps0, pf, pf0, vx, vy, vxf, vyf, dtm)
        vxm = zeros(4, 1)
        vym = zeros(4, 1)
        pr0_ver=deepcopy(pr_ver)
        ps0_ver=deepcopy(ps_ver)
        for jj in 2:1:Nx
            for ii in 2:1:Ny
                # Save initial nodal coordinates
                xcur=xp[jj]
                ycur=yp[ii]
                xA=xcur
                yA=ycur
                for rk in 1:1:4
                    # Interpolate vx
                    # Define i;j indexes for the upper left node
                    j=trunc(Int, (xcur-xvx[1])/dx)+1
                    i=trunc(Int, (ycur-yvx[1])/dy)+1
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
                    # Compute distances
                    dxmj=xcur-xvx[j]
                    dymi=ycur-yvx[i]
                    # Compute weights
                    # Compute vx velocity for the top & bottom of the cell()
                    vxm13=vx[i, j]*(1-dxmj/dx)+vx[i, j + 1]*dxmj/dx
                    vxm24=vx[i + 1, j]*(1-dxmj/dx)+vx[i + 1, j + 1]*dxmj/dx
                    # Compute correction
                    if dxmj/dx>=0.5
                        if j<Nx-1
                            vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(
                                vx[i, j]-2*vx[i, j + 1]+vx[i, j + 2]
                            )
                            vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(
                                vx[i + 1, j]-2*vx[i + 1, j + 1]+vx[i + 1, j + 2]
                            )
                        end
                    else
                        if j>1
                            vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(
                                vx[i, j - 1]-2*vx[i, j]+vx[i, j + 1]
                            )
                            vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(
                                vx[i + 1, j - 1]-2*vx[i + 1, j]+vx[i + 1, j + 1]
                            )
                        end
                    end
                    # Compute vx
                    vxm[rk]=(1-dymi/dy)*vxm13+(dymi/dy)*vxm24

                    # Interpolate vy
                    # Define i;j indexes for the upper left node
                    j=trunc(Int, (xcur-xvy[1])/dx)+1
                    i=trunc(Int, (ycur-yvy[1])/dy)+1
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
                    dxmj=xcur-xvy[j]
                    dymi=ycur-yvy[i]
                    # Compute weights
                    # Compute vy velocity for the left & right of the cell()
                    vym12=vy[i, j]*(1-dymi/dy)+vy[i + 1, j]*dymi/dy
                    vym34=vy[i, j + 1]*(1-dymi/dy)+vy[i + 1, j + 1]*dymi/dy
                    # Compute correction
                    if dymi/dy>=0.5
                        if i<Ny-1
                            vym12=vym12+1/2*((dymi/dy-0.5)^2)*(
                                vy[i, j]-2*vy[i + 1, j]+vy[i + 2, j]
                            )
                            vym34=vym34+1/2*((dymi/dy-0.5)^2)*(
                                vy[i, j + 1]-2*vy[i + 1, j + 1]+vy[i + 2, j + 1]
                            )
                        end
                    else
                        if i>1
                            vym12=vym12+1/2*((dymi/dy-0.5)^2)*(
                                vy[i - 1, j]-2*vy[i, j]+vy[i + 1, j]
                            )
                            vym34=vym34+1/2*((dymi/dy-0.5)^2)*(
                                vy[i - 1, j + 1]-2*vy[i, j + 1]+vy[i + 1, j + 1]
                            )
                        end
                    end
                    # Compute vy
                    vym[rk]=(1-dxmj/dx)*vym12+(dxmj/dx)*vym34

                    # Change coordinates to obtain B;C;D points
                    if rk==1 || rk==2
                        xcur=xA-dtm/2*vxm[rk]
                        ycur=yA-dtm/2*vym[rk]
                    elseif rk==3
                        xcur=xA-dtm*vxm[rk]
                        ycur=yA-dtm*vym[rk]
                    end
                end
                # Compute effective velocity
                vxmeff=1/6*(vxm[1]+2*vxm[2]+2*vxm[3]+vxm[4])
                vymeff=1/6*(vym[1]+2*vym[2]+2*vym[3]+vym[4])
                # Trace the node backward
                xcur=xA-dtm*vxmeff
                ycur=yA-dtm*vymeff
                # Interpolate nodal property
                # Define i;j indexes for the upper left node
                j=trunc(Int, (xcur-xp[1])/dx)+1
                i=trunc(Int, (ycur-yp[1])/dy)+1
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
                dxmj=xcur-xp[j]
                dymi=ycur-yp[i]
                # Compute weights
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Compute nodal total pressure
                pr0_ver[ii, jj]=pr_ver[i, j]*wtmij+pr_ver[i + 1, j]*wtmi1j+pr_ver[i, j + 1]*wtmij1+pr_ver[
                    i + 1, j + 1
                ]*wtmi1j1
                ps0_ver[ii, jj]=ps_ver[i, j]*wtmij+ps_ver[i + 1, j]*wtmi1j+ps_ver[i, j + 1]*wtmij1+ps_ver[
                    i + 1, j + 1
                ]*wtmi1j1
            end
        end

        # Backtracing Pressure nodes: Pfluid
        pf0_ver=deepcopy(pf_ver)
        for jj in 2:1:Nx
            for ii in 2:1:Ny
                # Save initial nodal coordinates
                xcur=xp[jj]
                ycur=yp[ii]
                xA=xcur
                yA=ycur
                for rk in 1:1:4
                    # Interpolate vx
                    # Define i;j indexes for the upper left node
                    j=trunc(Int, (xcur-xvx[1])/dx)+1
                    i=trunc(Int, (ycur-yvx[1])/dy)+1
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
                    # Compute distances
                    dxmj=xcur-xvx[j]
                    dymi=ycur-yvx[i]
                    # Compute weights
                    # Compute vx velocity for the top & bottom of the cell()
                    vxm13=vxf[i, j]*(1-dxmj/dx)+vxf[i, j + 1]*dxmj/dx
                    vxm24=vxf[i + 1, j]*(1-dxmj/dx)+vxf[i + 1, j + 1]*dxmj/dx
                    # Compute correction
                    if dxmj/dx>=0.5
                        if j<Nx-1
                            vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(
                                vxf[i, j]-2*vxf[i, j + 1]+vxf[i, j + 2]
                            )
                            vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(
                                vxf[i + 1, j]-2*vxf[i + 1, j + 1]+vxf[i + 1, j + 2]
                            )
                        end
                    else
                        if j>1
                            vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(
                                vxf[i, j - 1]-2*vxf[i, j]+vxf[i, j + 1]
                            )
                            vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(
                                vxf[i + 1, j - 1]-2*vxf[i + 1, j]+vxf[i + 1, j + 1]
                            )
                        end
                    end
                    # Compute vx
                    vxm[rk]=(1-dymi/dy)*vxm13+(dymi/dy)*vxm24

                    # Interpolate vy
                    # Define i;j indexes for the upper left node
                    j=trunc(Int, (xcur-xvy[1])/dx)+1
                    i=trunc(Int, (ycur-yvy[1])/dy)+1
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
                    dxmj=xcur-xvy[j]
                    dymi=ycur-yvy[i]
                    # Compute weights
                    # Compute vy velocity for the left & right of the cell()
                    vym12=vyf[i, j]*(1-dymi/dy)+vyf[i + 1, j]*dymi/dy
                    vym34=vyf[i, j + 1]*(1-dymi/dy)+vyf[i + 1, j + 1]*dymi/dy
                    # Compute correction
                    if dymi/dy>=0.5
                        if i<Ny-1
                            vym12=vym12+1/2*((dymi/dy-0.5)^2)*(
                                vyf[i, j]-2*vyf[i + 1, j]+vyf[i + 2, j]
                            )
                            vym34=vym34+1/2*((dymi/dy-0.5)^2)*(
                                vyf[i, j + 1]-2*vyf[i + 1, j + 1]+vyf[i + 2, j + 1]
                            )
                        end
                    else
                        if i>1
                            vym12=vym12+1/2*((dymi/dy-0.5)^2)*(
                                vyf[i - 1, j]-2*vyf[i, j]+vyf[i + 1, j]
                            )
                            vym34=vym34+1/2*((dymi/dy-0.5)^2)*(
                                vyf[i - 1, j + 1]-2*vyf[i, j + 1]+vyf[i + 1, j + 1]
                            )
                        end
                    end
                    # Compute vy
                    vym[rk]=(1-dxmj/dx)*vym12+(dxmj/dx)*vym34

                    # Change coordinates to obtain B;C;D points
                    if rk==1 || rk==2
                        xcur=xA-dtm/2*vxm[rk]
                        ycur=yA-dtm/2*vym[rk]
                    elseif rk==3
                        xcur=xA-dtm*vxm[rk]
                        ycur=yA-dtm*vym[rk]
                    end
                end

                # Compute effective velocity
                vxmeff=1/6*(vxm[1]+2*vxm[2]+2*vxm[3]+vxm[4])
                vymeff=1/6*(vym[1]+2*vym[2]+2*vym[3]+vym[4])
                # Trace the node backward
                xcur=xA-dtm*vxmeff
                ycur=yA-dtm*vymeff
                # Interpolate nodal property
                # SIGMA'xx; P
                # Define i;j indexes for the upper left node
                j=trunc(Int, (xcur-xp[1])/dx)+1
                i=trunc(Int, (ycur-yp[1])/dy)+1
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
                dxmj=xcur-xp[j]
                dymi=ycur-yp[i]
                # Compute weights
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy)
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Compute nodal total pressure
                pf0_ver[ii, jj]=pf_ver[i, j]*wtmij+pf_ver[i + 1, j]*wtmi1j+pf_ver[i, j + 1]*wtmij1+pf_ver[
                    i + 1, j + 1
                ]*wtmi1j1
            end
        end
        # test
        for jj in 2:1:Nx, ii in 2:1:Ny
            @test pr0[ii, jj] ≈ pr0_ver[ii, jj] rtol=1e-9
            @test ps0[ii, jj] ≈ ps0_ver[ii, jj] rtol=1e-9
            @test pf0[ii, jj] ≈ pf0_ver[ii, jj] rtol=1e-9
        end
    end # testset "backtrace_pressures_rk4!()"

    @testset "replenish_markers!()" begin
        marknum = start_marknum
        xm = rand(rgen, (-dx):0.1:(x[end] / 2), marknum)
        ym = rand(rgen, (-dy):0.1:(y[end] / 2), marknum)
        tm = rand(rgen, 1:3, marknum)
        tkm = rand(rgen, marknum)
        sxxm = rand(rgen, marknum)
        sxym = rand(rgen, marknum)
        phim = rand(rgen, marknum)
        etavpm = rand(rgen, marknum)
        phinewm = rand(rgen, marknum)
        pfm0 = rand(rgen, marknum)
        XWsolidm = rand(rgen, marknum)
        XWsolidm0 = rand(rgen, marknum)
        rhototalm = rand(rgen, marknum)
        rhocptotalm = rand(rgen, marknum)
        etatotalm = rand(rgen, marknum)
        hrtotalm = rand(rgen, marknum)
        ktotalm = rand(rgen, marknum)
        inv_gggtotalm = rand(rgen, marknum)
        fricttotalm = rand(rgen, marknum)
        cohestotalm = rand(rgen, marknum)
        tenstotalm = rand(rgen, marknum)
        rhofluidcur = rand(rgen, marknum)
        alphasolidcur = rand(rgen, marknum)
        alphafluidcur = rand(rgen, marknum)
        tkm_rhocptotalm = rand(rgen, marknum)
        etafluidcur_inv_kphim = rand(rgen, marknum)
        mdis, mnum = Erebus.setup_marker_geometry_helpers()
        xm_ver = deepcopy(xm)
        ym_ver = deepcopy(ym)
        tm_ver = deepcopy(tm)
        tkm_ver = deepcopy(tkm)
        sxxm_ver = deepcopy(sxxm)
        sxym_ver = deepcopy(sxym)
        phim_ver = deepcopy(phim)
        etavpm_ver = deepcopy(etavpm)
        phinewm_ver = deepcopy(phinewm)
        pfm0_ver = deepcopy(pfm0)
        XWsolidm_ver = deepcopy(XWsolidm)
        XWsolidm0_ver = deepcopy(XWsolidm0)
        # replenish_markers
        marknum_new = Erebus.replenish_markers!(
            xm,
            ym,
            tm,
            tkm,
            phim,
            sxxm,
            sxym,
            etavpm,
            phinewm,
            pfm0,
            XWsolidm,
            XWsolidm0,
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
            etafluidcur_inv_kphim,
            mdis,
            mnum;
            randomized=false,
        )
        # verification, from HTM-planetary.m, line 2491ff
        # Add markers to empty areas
        marknumold=marknum
        mdis_ver=1e30 .* ones(Nym, Nxm)
        mnum_ver = zeros(Int, Nym, Nxm)
        mtyp = zeros(Int, Nym, Nxm)
        mpor = zeros(Nym, Nxm)
        for m in 1:1:marknum
            # Check markers with the nearest nodes
            # Define i;j indexes for the upper left node
            j=trunc(Int, (xm_ver[m]-xxm[1])/dxm)+1
            i=trunc(Int, (ym_ver[m]-yym[1])/dym)+1
            if j<1
                j=1
            elseif j>Nxm-1
                j=Nxm-1
            end
            if i<1
                i=1
            elseif i>Nym-1
                i=Nym-1
            end
            # Check nodes
            # i;j Node
            # Compute distance
            dxmj=xm_ver[m]-xxm[j]
            dymi=ym_ver[m]-yym[i]
            dismij=(dxmj^2+dymi^2)^0.5
            if dismij<mdis_ver[i, j]
                mdis_ver[i, j]=dismij
                mnum_ver[i, j]=m
                mtyp[i, j]=tm_ver[m]
                mpor[i, j]=phim[m]
            end
            # i+1;j Node
            # Compute distance
            dxmj=xm_ver[m]-xxm[j]
            dymi=ym_ver[m]-yym[i + 1]
            dismi1j=(dxmj^2+dymi^2)^0.5
            if dismi1j<mdis_ver[i + 1, j]
                mdis_ver[i + 1, j]=dismi1j
                mnum_ver[i + 1, j]=m
                mtyp[i + 1, j]=tm_ver[m]
                mpor[i + 1, j]=phim[m]
            end
            # i;j+1 Node
            # Compute distance
            dxmj=xm_ver[m]-xxm[j + 1]
            dymi=ym_ver[m]-yym[i]
            dismij1=(dxmj^2+dymi^2)^0.5
            if dismij1<mdis_ver[i, j + 1]
                mdis_ver[i, j + 1]=dismij1
                mnum_ver[i, j + 1]=m
                mtyp[i, j + 1]=tm_ver[m]
                mpor[i, j + 1]=phim[m]
            end
            # i+1;j+1 Node
            # Compute distance
            dxmj=xm_ver[m]-xxm[j + 1]
            dymi=ym_ver[m]-yym[i + 1]
            dismi1j1=(dxmj^2+dymi^2)^0.5
            if dismi1j1<mdis_ver[i + 1, j + 1]
                mdis_ver[i + 1, j + 1]=dismi1j1
                mnum_ver[i + 1, j + 1]=m
                mtyp[i + 1, j + 1]=tm_ver[m]
                mpor[i + 1, j + 1]=phim[m]
            end
        end
        dii=5*Nxmc
        djj=5*Nymc
        for j in 1:1:Nxm
            for i in 1:1:Nym
                if mnum_ver[i, j]==0
                    # Serch surrounding nodes
                    for jj in (j - djj):1:(j + djj)
                        for ii in (i - dii):1:(i + dii)
                            if ii>=1 && ii<=Nym && jj>=1 && jj<=Nxm && mnum_ver[ii, jj]>0
                                # Compute distance
                                m=mnum_ver[ii, jj]
                                dxmj=xm_ver[m]-xxm[j]
                                dymi=ym_ver[m]-yym[i]
                                dismij=(dxmj^2+dymi^2)^0.5
                                if dismij<mdis_ver[i, j]
                                    mdis_ver[i, j]=dismij
                                    mnum_ver[i, j]=-m
                                    mtyp[i, j]=-tm_ver[m]
                                    mpor[i, j]=phim[m]
                                end
                            end
                        end
                    end
                    # Add New marker
                    if mnum_ver[i, j]<0
                        # Add marker number
                        marknum=marknum+1
                        # Assign marker coordinates
                        push!(xm_ver, xxm[j])#+(rand(rgen)-0.5)*dxm)
                        push!(ym_ver, yym[i])#+(rand(rgen)-0.5)*dym)
                        # Copy marker properties
                        m=-mnum_ver[i, j]
                        push!(tm_ver, tm[m]) # Material type()
                        push!(tkm_ver, tkm[m]) # Temperature
                        push!(phim_ver, phim[m]) # Porosity
                        push!(sxxm_ver, sxxm[m]) # SIGMA'xx, Pa
                        push!(sxym_ver, sxym[m]) # SIGMAxy, Pa
                        push!(etavpm_ver, etavpm[m]) # Visco-plastic viscosity, Pa
                        push!(phinewm_ver, phinewm[m]) # New porosity
                        push!(pfm0_ver, pfm0[m]) # pfm0
                        push!(XWsolidm_ver, XWsolidm[m]) # XMsolid
                        push!(XWsolidm0_ver, XWsolidm0[m]) # XMsolidnew
                    end
                end
            end
        end
        marknum_new_ver=marknum
        # test
        @test length(xm) == length(xm_ver)
        @test length(ym) == length(ym_ver)
        @test marknum_new == marknum_new_ver
        @test tm == tm_ver
        @test tkm == tkm_ver
        @test phim == phim_ver
        @test sxxm == sxxm_ver
        @test sxym == sxym_ver
        @test etavpm == etavpm_ver
        @test phinewm == phinewm_ver
        @test pfm0 == pfm0_ver
        @test XWsolidm == XWsolidm_ver
        @test XWsolidm0 == XWsolidm0_ver
    end # testset "replenish_markers!()"
end

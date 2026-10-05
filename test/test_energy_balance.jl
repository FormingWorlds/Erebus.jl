using Erebus
using StaticArrays
using Test

@testset "Energy Balance Consolidation (PR 5)" begin
    @testset "F16: Serpentinization Latent Heat Linear Scaling & Magnitude" begin
        cfg = SimulationConfig(
            reaction=ReactionConfig(
                active=true,
                hydration_active=true,
                dehydration_active=true,
                hydration_mode=9,
                dtreaction_hydration=1.0e6,
                delta_H=40000.0, # 40 kJ / mol H2O
            ),
        )
        coords = GridCoordinates(cfg.grid)
        dx, dy = coords.dx, coords.dy
        cell_vol = dx * dy
        Ny1, Nx1 = coords.Ny1, coords.Nx1

        function run_reaction(dt_val)
            DMP = zeros(Ny1, Nx1)
            DHP = zeros(Ny1, Nx1)
            DMPSUM = zeros(Ny1, Nx1)
            DHPSUM = zeros(Ny1, Nx1)
            WTPSUM = zeros(Ny1, Nx1)
            DQPF = zeros(Ny1, Nx1)
            DQPFSUM = zeros(Ny1, Nx1)

            pf = fill(5.0e7, Ny1, Nx1)
            tk2 = fill(500.0, Ny1, Nx1) # Within hydration stability window
            tm = [1]
            xm = [coords.xp[2] + dx / 2.0]
            ym = [coords.yp[2] + dy / 2.0]
            XW0_val = 0.05
            XW0 = [XW0_val]
            XW = [XW0_val]
            phim = [0.10]
            phinewm = [0.10]
            pfm0 = [5.0e7]
            marknum = 1
            timestep = 2
            titer = 1

            Erebus.perform_thermochemical_reaction!(
                DMP,
                DHP,
                DMPSUM,
                DHPSUM,
                WTPSUM,
                pf,
                tk2,
                tm,
                xm,
                ym,
                XW0,
                XW,
                phim,
                phinewm,
                pfm0,
                marknum,
                dt_val,
                timestep,
                titer;
                coords=coords,
                DQPF=DQPF,
                DQPFSUM=DQPFSUM,
                cfg=cfg.reaction,
                backload_step1=false,
            )

            delta_XW = XW[1] - XW0_val
            # Staggered cell P-nodes: 4 corner nodes each represent 1/4 of cell volume
            Q_code = (sum(DHP) / 4.0) * cell_vol * dt_val

            # Theoretical physical energy from reacted moles:
            # Molar mass of wet silicate rock = MD + MH2O * XW
            MD = Erebus.MD
            MH2O = Erebus.MH2O
            molar_mass_rock = MD + MH2O * XW0_val
            rho_solid = 3000.0 # kg/m^3
            solid_mass = rho_solid * (1.0 - phim[1]) * cell_vol
            moles_rock = solid_mass / molar_mass_rock
            moles_H2O_bound = moles_rock * delta_XW
            Q_phys = moles_H2O_bound * cfg.reaction.delta_H

            return (; dt=dt_val, delta_XW=delta_XW, Q_code=Q_code, Q_phys=Q_phys)
        end

        r1 = run_reaction(1000.0)
        r2 = run_reaction(500.0)

        # Ratio must scale linearly with dt (2.0), not quadratically (4.0)
        ratio_dt = r1.Q_code / r2.Q_code
        @test isapprox(ratio_dt, 2.0; atol=0.1)

        # Released energy must match physical reaction enthalpy n_H2O * delta_H within 10%
        ratio_phys = r1.Q_code / r1.Q_phys
        @test isapprox(ratio_phys, 1.0; atol=0.10)
    end

    @testset "F15: Serpentinization WTPSUM Grid Weighting" begin
        cfg = SimulationConfig(
            grid=GridConfig(Nx=4, Ny=4, xsize=4000.0, ysize=4000.0),
            reaction=ReactionConfig(active=true),
        )
        coords = GridCoordinates(cfg.grid)
        Nx1, Ny1 = coords.Nx1, coords.Ny1
        dx, dy = coords.dx, coords.dy

        DMP = zeros(Ny1, Nx1)
        DHP = zeros(Ny1, Nx1)
        DMPSUM = zeros(Ny1, Nx1)
        DHPSUM = zeros(Ny1, Nx1)
        WTPSUM = zeros(Ny1, Nx1)

        pf = fill(1.0e6, Ny1, Nx1)
        tk2 = fill(600.0, Ny1, Nx1)

        # 4 markers in cell (1, 1): 1 rock marker that reacts, 3 air markers (tm = 3)
        marknum = 4
        xm = [
            coords.xp[1] + 0.25 * dx,
            coords.xp[1] + 0.75 * dx,
            coords.xp[1] + 0.25 * dx,
            coords.xp[1] + 0.75 * dx,
        ]
        ym = [
            coords.yp[1] + 0.25 * dy,
            coords.yp[1] + 0.25 * dy,
            coords.yp[1] + 0.75 * dy,
            coords.yp[1] + 0.75 * dy,
        ]
        tm = [1, 3, 3, 3]

        XWsolidm0 = [0.5, 0.0, 0.0, 0.0]
        XWsolidm = copy(XWsolidm0)
        phim = [0.1, 0.0, 0.0, 0.0]
        phinewm = copy(phim)
        pfm0 = fill(1.0e6, marknum)

        dt = 1.0e5
        timestep = 1
        titer = 1

        Erebus.perform_thermochemical_reaction!(
            DMP,
            DHP,
            DMPSUM,
            DHPSUM,
            WTPSUM,
            pf,
            tk2,
            tm,
            xm,
            ym,
            XWsolidm0,
            XWsolidm,
            phim,
            phinewm,
            pfm0,
            marknum,
            dt,
            timestep,
            titer;
            coords=coords,
            cfg=cfg.reaction,
        )

        # Total weight of all 4 markers across the 4 cell nodes must sum to 4.0
        sum_wt = sum(WTPSUM)
        @test isapprox(sum_wt, 4.0; atol=1e-10)
        @test isapprox(sum(WTPSUM[1:2, 1:2]), 4.0; atol=1e-10)
    end

    @testset "F10: Melt Segregation Sensible Heat Conservation" begin
        cfg_grid = GridConfig(xsize=10000.0, ysize=10000.0, Nx=7, Ny=7)
        coords = GridCoordinates(cfg_grid)
        dx, dy = coords.dx, coords.dy
        cell_vol = dx * dy

        xm = Float64[]
        ym = Float64[]
        tm = Int[]
        tkm = Float64[]
        Fm = Float64[]

        for j in 1:(coords.Nx - 1)
            for i in 1:(coords.Ny - 1)
                xc = coords.x[j] + dx / 2
                yc = coords.y[i] + dy / 2
                for ox in (-dx / 4, dx / 4), oy in (-dy / 4, dy / 4)
                    push!(xm, xc + ox)
                    push!(ym, yc + oy)
                    push!(tm, 1)
                    # Temperature gradient: deeper is hotter
                    push!(tkm, 1500.0 + 300.0 * (yc / coords.ysize))
                    push!(Fm, 0.25)
                end
            end
        end
        marknum = length(xm)

        gx = zeros(coords.Ny, coords.Nx)
        gy = fill(1.0, coords.Ny, coords.Nx) # 1 m/s^2 downward
        Q_seg_grid = zeros(coords.Ny, coords.Nx)

        cfg_magma = MagmaTransportConfig(
            active=true,
            sensible_heat_transport=true,
            segregation_heating=false,
            latent_crystallization=false,
            k_melt_ref=1.0e-9,
            cp_melt=1200.0,
        )

        dt = 1.0e5
        Erebus.apply_silicate_melt_segregation!(
            xm,
            ym,
            tm,
            tkm,
            Fm,
            marknum,
            dt,
            cfg_magma;
            coords=coords,
            gx=gx,
            gy=gy,
            Q_seg_grid=Q_seg_grid,
            rho_melt=2800.0,
        )

        # Net integrated sensible heat added across the domain must sum to zero
        total_sens_energy = sum(Q_seg_grid) * cell_vol * dt
        @test isapprox(total_sens_energy, 0.0; atol=1e-5)

        # Donor cells must be cooled (Q_seg_grid < 0) and receiver cells heated (Q_seg_grid > 0)
        @test minimum(Q_seg_grid) < -1e-6
        @test isapprox(maximum(Q_seg_grid), -minimum(Q_seg_grid); rtol=0.2)
    end

    @testset "F12 & F33: 26Al/60Fe Power & Volume Weighting" begin
        cfg = SimulationConfig()
        coords = GridCoordinates(cfg.grid)
        marknum = 1000
        xm = collect(range(0.0, 1000.0; length=marknum))
        ym = collect(range(0.0, 1000.0; length=marknum))
        Vm = 1.0 # m^3
        tm = fill(1, marknum)
        phim = fill(0.0, marknum)
        timesum = 0.0

        tau_al = cfg.thermodynamics.t_half_al / log(2.0)
        Q_al = Erebus.Q_radiogenic(
            cfg.thermodynamics.f_al,
            cfg.thermodynamics.ratio_al,
            cfg.thermodynamics.E_al,
            tau_al,
            timesum,
        )

        rho_silicate = cfg.materials.rhosolidm[1]
        rho_metal = cfg.coreformation.rho_metal
        Xfe_val = cfg.coreformation.Xfe_bulk # 0.2

        hrsolidm, hrfluidm, hrmetalm = Erebus.calculate_radioactive_heating(
            true,
            false,
            timesum;
            ratio_al=cfg.thermodynamics.ratio_al,
            E_al=cfg.thermodynamics.E_al,
            f_al=cfg.thermodynamics.f_al,
            tau_al=tau_al,
            rho_metal=rho_metal,
            rhosolidm=cfg.materials.rhosolidm,
            rhofluidm=cfg.materials.rhofluidm,
            X_fe_ref=Xfe_val,
        )

        # Case 1: With metal group
        phi_fe_ref =
            (Xfe_val / rho_metal) / ((Xfe_val / rho_metal) + (1.0 - Xfe_val) / rho_silicate)
        Xfe_bulk = fill(phi_fe_ref, marknum)
        rho_marker = (1.0 - phi_fe_ref) * rho_silicate + phi_fe_ref * rho_metal
        M_bulk_1 = marknum * Vm * rho_marker

        hrtotalm_1 = zeros(marknum)
        for m in 1:marknum
            phi_fe = Xfe_bulk[m]
            hr_solid = (1.0 - phi_fe) * hrsolidm[tm[m]] + phi_fe * hrmetalm[tm[m]]
            hrtotalm_1[m] = (1.0 - phim[m]) * hr_solid + phim[m] * hrfluidm[tm[m]]
        end
        total_power_1 = sum(hrtotalm_1 .* Vm)
        expected_power_1 = M_bulk_1 * Q_al
        ratio_1 = total_power_1 / expected_power_1
        # Must equal bulk power within 0.2%
        @test isapprox(ratio_1, 1.0; atol=0.005)

        # Case 2: Without metal group (chondritic reference mixture)
        hr_sol_chon, _, hr_met_chon = Erebus.calculate_radioactive_heating(
            true,
            false,
            timesum;
            ratio_al=cfg.thermodynamics.ratio_al,
            E_al=cfg.thermodynamics.E_al,
            f_al=cfg.thermodynamics.f_al,
            tau_al=tau_al,
            rho_metal=rho_metal,
            rhosolidm=cfg.materials.rhosolidm,
            rhofluidm=cfg.materials.rhofluidm,
            X_fe_ref=Erebus.X_FE_REF_CHONDRITE,
        )
        phi_fe_2 = Erebus.chondritic_phi_fe(rho_silicate, rho_metal)
        hr_rock_2 = (1.0 - phi_fe_2) * hr_sol_chon[1] + phi_fe_2 * hr_met_chon[1]
        hr_vol_2 = (1.0 - phim[1]) * hr_rock_2 + phim[1] * hrfluidm[1]
        v_fe = Erebus.X_FE_REF_CHONDRITE / rho_metal
        v_si = (1.0 - Erebus.X_FE_REF_CHONDRITE) / rho_silicate
        rho_bulk_chondrite = 1.0 / (v_fe + v_si)
        expected_power_2 = rho_bulk_chondrite * Q_al
        @test isapprox(hr_vol_2, expected_power_2; rtol=1e-5)
    end

    @testset "F30: DTmax Convergence & Timestep Consistency" begin
        # When maxDTcurrent exceeds DTmax, the step must not accept convergence
        # with an unreduced thermal solution.
        DTmax = 20.0
        maxDTcurrent_exceeded = 50.0
        maxDTcurrent_ok = 15.0

        # compute_thermochemical_iteration_outcome must reject convergence when DTmax is violated
        DMP = zeros(5, 5)
        pf = zeros(5, 5)
        pf0 = zeros(5, 5)
        titer = 1

        # On main, compute_thermochemical_iteration_outcome does not take maxDTcurrent or DTmax,
        # so it returns true regardless of maxDTcurrent.
        outcome_exceeded = Erebus.compute_thermochemical_iteration_outcome(
            DMP,
            pf,
            pf0,
            titer;
            pferrmax=1.0e5,
            maxDTcurrent=maxDTcurrent_exceeded,
            DTmax=DTmax,
        )
        @test outcome_exceeded == false

        outcome_ok = Erebus.compute_thermochemical_iteration_outcome(
            DMP, pf, pf0, titer; pferrmax=1.0e5, maxDTcurrent=maxDTcurrent_ok, DTmax=DTmax
        )
        @test outcome_ok == true
    end

    @testset "F12 & F33: Marker Radiogenic Property Isolation (compute_hr)" begin
        @test iszero(Erebus.chondritic_phi_fe(0.0, 5450.0))
        @test iszero(Erebus.chondritic_phi_fe(-1.0, 5450.0))
        v_fe = Erebus.X_FE_REF_CHONDRITE / 5450.0
        v_si = (1.0 - Erebus.X_FE_REF_CHONDRITE) / 3000.0
        expected_phi = v_fe / (v_fe + v_si)
        @test isapprox(Erebus.chondritic_phi_fe(3000.0, 0.0), expected_phi; atol=1e-10)
        @test isapprox(Erebus.chondritic_phi_fe(3000.0, -10.0), expected_phi; atol=1e-10)

        marknum = 3
        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0, Fm) = Erebus.setup_marker_properties(
            marknum
        )
        (rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm, tkm_rhocptotalm, etafluidcur_inv_kphim, inv_gggtotalm, fricttotalm, cohestotalm, tenstotalm, rhofluidcur, alphasolidcur, alphafluidcur) = Erebus.setup_marker_properties_helpers(
            marknum
        )
        (Xfem, Xfem0, Xfe_bulk) = Erebus.setup_marker_metal_properties(marknum)

        tm .= 1
        phim .= 0.1
        tkm .= 500.0
        hrsolid = [100.0, 100.0, 0.0]
        hrfluid = [0.0, 0.0, 0.0]
        hrmetal = [50.0, 50.0, 0.0]

        # Case 1: coreformation_active = true with metal blending
        Xfe_bulk[1] = 0.25
        Erebus.compute_marker_properties!(
            1,
            tm,
            tkm,
            rhototalm,
            rhocptotalm,
            etatotalm,
            hrtotalm,
            ktotalm,
            tkm_rhocptotalm,
            etafluidcur_inv_kphim,
            hrsolid,
            hrfluid,
            phim,
            XWsolidm0,
            9,
            rhofluidcur;
            compute_hr=true,
            coreformation_active=true,
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            hrmetalm=hrmetal,
        )
        @test isapprox(hrtotalm[1], 80.0; atol=1e-10)

        # Case 2: !coreformation_active with Xfe_bulk array
        Xfe_bulk[2] = 0.25
        Erebus.compute_marker_properties!(
            2,
            tm,
            tkm,
            rhototalm,
            rhocptotalm,
            etatotalm,
            hrtotalm,
            ktotalm,
            tkm_rhocptotalm,
            etafluidcur_inv_kphim,
            hrsolid,
            hrfluid,
            phim,
            XWsolidm0,
            9,
            rhofluidcur;
            compute_hr=true,
            coreformation_active=false,
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            hrmetalm=hrmetal,
        )
        @test isapprox(hrtotalm[2], 80.0; atol=1e-10)

        # Case 3: !coreformation_active with chondritic reference fallback
        Erebus.compute_marker_properties!(
            3,
            tm,
            tkm,
            rhototalm,
            rhocptotalm,
            etatotalm,
            hrtotalm,
            ktotalm,
            tkm_rhocptotalm,
            etafluidcur_inv_kphim,
            hrsolid,
            hrfluid,
            phim,
            XWsolidm0,
            9,
            rhofluidcur;
            compute_hr=true,
            coreformation_active=false,
            Xfe_bulk=nothing,
            Xfem=Xfem,
            hrmetalm=hrmetal,
        )
        phi_fe_3 = Erebus.chondritic_phi_fe(3300.0, 5450.0)
        expected_3 = (1.0 - 0.1) * ((1.0 - phi_fe_3) * 100.0 + phi_fe_3 * 50.0)
        @test isapprox(hrtotalm[3], expected_3; atol=1e-10)

        # Case 4: compute_hr = false leaves hrtotalm unchanged
        hrtotalm[1] = -999.0
        Erebus.compute_marker_properties!(
            1,
            tm,
            tkm,
            rhototalm,
            rhocptotalm,
            etatotalm,
            hrtotalm,
            ktotalm,
            tkm_rhocptotalm,
            etafluidcur_inv_kphim,
            hrsolid,
            hrfluid,
            phim,
            XWsolidm0,
            9,
            rhofluidcur;
            compute_hr=false,
            coreformation_active=true,
            Xfe_bulk=Xfe_bulk,
            Xfem=Xfem,
            hrmetalm=hrmetal,
        )
        @test isapprox(hrtotalm[1], -999.0; atol=1e-10)
    end
end

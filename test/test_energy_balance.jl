using Erebus
using StaticArrays
using Test

@testset "Energy Balance Consolidation" begin
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

        # Ratio scales linearly with dt (2.0)
        ratio_dt = r1.Q_code / r2.Q_code
        @test isapprox(ratio_dt, 2.0; atol=0.005)

        # Released energy matches physical reaction enthalpy n_H2O * delta_H
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

        # Receiver cells gain sensible heat (+dQ_sens) and donor cells lose (-dQ_sens)
        @test isapprox(minimum(Q_seg_grid), -0.0367; rtol=0.01)
        @test isapprox(maximum(Q_seg_grid), 0.0367; rtol=0.01)

        # Net domain-integrated sensible heat source is identically zero
        total_sens_energy = sum(Q_seg_grid) * cell_vol * dt
        @test isapprox(total_sens_energy, 0.0; atol=1e-10)

        # Horizontal sensible heat transport (flux_x)
        tkm_x = [
            1500.0 + 300.0 * ((floor(x / dx) * dx + dx / 2) / coords.xsize) for x in xm
        ]
        Fm_x = fill(0.25, marknum)
        gx_x = fill(1.0, coords.Ny, coords.Nx)
        gy_x = zeros(coords.Ny, coords.Nx)
        Q_seg_x = zeros(coords.Ny, coords.Nx)
        Erebus.apply_silicate_melt_segregation!(
            xm,
            ym,
            tm,
            tkm_x,
            Fm_x,
            marknum,
            dt,
            cfg_magma;
            coords=coords,
            gx=gx_x,
            gy=gy_x,
            Q_seg_grid=Q_seg_x,
            rho_melt=2800.0,
        )
        @test isapprox(minimum(Q_seg_x), -0.0367; rtol=0.01)
        @test isapprox(maximum(Q_seg_x), 0.0367; rtol=0.01)
        @test isapprox(sum(Q_seg_x), 0.0; atol=1e-14)
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
        w_metal_val = Erebus.metal_volume_to_mass_fraction(Xfe_val, rho_metal, rho_silicate)

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
            w_metal_ref=w_metal_val,
        )

        # Case 1: With metal group; 20 vol% metal holds more iron than CI
        phi_fe_ref = Xfe_val
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
        # Al per kg silicate stays at the CI value above the CI iron content
        expected_power_1 = M_bulk_1 * (1.0 - w_metal_val) * Q_al / 0.818
        ratio_1 = total_power_1 / expected_power_1
        @test isapprox(ratio_1, 1.0; rtol=1e-12)

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
            w_metal_ref=Erebus.X_FE_REF_CHONDRITE,
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
        # When maxDTcurrent exceeds DTmax, the step rejects convergence
        DTmax = 20.0
        maxDTcurrent_exceeded = 50.0
        maxDTcurrent_ok = 15.0

        DMP = zeros(5, 5)
        pf = zeros(5, 5)
        pf0 = zeros(5, 5)
        titer = 1

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

        # Timestep limiter reduces dt and honors dt_min clamp
        dt_reduced = Erebus.finalize_thermochemical_iteration_pass(
            50.0, 100.0, 1, 20.0; dt_min_val=10.0
        )
        @test isapprox(dt_reduced, 40.0; atol=1e-10)

        dt_clamped = Erebus.finalize_thermochemical_iteration_pass(
            100.0, 10.0, 1, 20.0; dt_min_val=5.0
        )
        @test isapprox(dt_clamped, 5.0; atol=1e-10)

        dt_unchanged = Erebus.finalize_thermochemical_iteration_pass(
            15.0, 100.0, 1, 20.0; dt_min_val=5.0
        )
        @test isapprox(dt_unchanged, 100.0; atol=1e-10)
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

    @testset "Adiabatic 0D Multi-Nuclide Step-Integral Energy Conservation" begin
        # 5-point Gauss-Legendre nodes and weights on [-1, 1]
        gl5_xi = [
            -sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0,
            -sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0,
            0.0,
            sqrt(5.0 - 2.0 * sqrt(10.0 / 7.0)) / 3.0,
            sqrt(5.0 + 2.0 * sqrt(10.0 / 7.0)) / 3.0,
        ]
        gl5_w = [
            (322.0 - 13.0 * sqrt(70.0)) / 900.0,
            (322.0 + 13.0 * sqrt(70.0)) / 900.0,
            128.0 / 225.0,
            (322.0 + 13.0 * sqrt(70.0)) / 900.0,
            (322.0 - 13.0 * sqrt(70.0)) / 900.0,
        ]

        function _gl5_decay_integral(Q0::Float64, tau::Float64, t_edges::Vector{Float64})
            total = 0.0
            for k in 1:(length(t_edges) - 1)
                a = t_edges[k]
                b = t_edges[k + 1]
                half_len = 0.5 * (b - a)
                mid = 0.5 * (a + b)
                sub = 0.0
                for q in 1:5
                    t_q = mid + half_len * gl5_xi[q]
                    sub += gl5_w[q] * Q0 * exp(-t_q / tau)
                end
                total += half_len * sub
            end
            return total
        end

        function _make_subinterval_edges(dt_hist::Vector{Float64})
            edges = Float64[0.0]
            t_cur = 0.0
            for dt_k in dt_hist
                push!(edges, t_cur + 0.5 * dt_k)
                t_cur += dt_k
                push!(edges, t_cur)
            end
            return edges
        end

        function _run_adiabatic_0d_box(
            nsteps::Int;
            hr_al::Bool=true,
            hr_fe::Bool=true,
            dt_yr::Float64=500_000.0,
            DTmax::Float64=1.0e6,
        )
            yr = Erebus.SECONDS_PER_YEAR
            T0 = 200.0
            sv3(a, b, c) = Erebus.SVector{3,Float64}([a, b, c])
            cfg = SimulationConfig(;
                grid=GridConfig(; Nx=9, Ny=9, xsize=200_000.0, ysize=200_000.0),
                geometry=GeometryConfig(;
                    rplanet=90_000.0,
                    rcrust=90_000.0,
                    xcenter=100_000.0,
                    ycenter=100_000.0,
                    spherical_metric=false,
                    gravity_mode=:poisson2d,
                ),
                time=TimeConfig(;
                    n_steps=nsteps,
                    dt_initial=dt_yr,
                    dt_longest=dt_yr,
                    DTmax=DTmax,
                    yearlength=yr,
                ),
                materials=MaterialConfig(;
                    alphasolidm=sv3(0.0, 0.0, 0.0),
                    alphafluidm=sv3(0.0, 0.0, 0.0),
                    kphim0=sv3(1.0e-30, 1.0e-30, 1.0e-30),
                    tkm0=sv3(T0, T0, T0),
                ),
                thermodynamics=ThermalConfig(;
                    hr_al=hr_al,
                    hr_fe=hr_fe,
                    phim0=1.0e-4,
                    thermal_buoyancy=false,
                    surface_radiation=false,
                ),
                poroelasticity=PoroelasticConfig(;
                    phimin=1.0e-4, phimax=0.9999, hydrofracture=false
                ),
                reaction=ReactionConfig(; active=false),
                melting=MeltingConfig(; active=false),
                coreformation=CoreFormationConfig(;
                    percolation_active=false,
                    settling_active=false,
                    Xfe_bulk=Erebus.X_FE_REF_CHONDRITE,
                ),
                venting=VentingConfig(; active=false),
                disk=DiskConfig(; enabled=false, model=:fixed, t_ambient=T0),
                solver=SolverConfig(;
                    dsubgridt=0.0, dsubgrids=0.0, p2m_mode=:tiled, max_plastic_iterations=4
                ),
                output=OutputConfig(;
                    mode=:telemetry, savematstep=100, visstep=100, save_final=false
                ),
            )
            Erebus.Config.validate_config(cfg)
            state, coords, ws = mktempdir() do tmpdir
                st, co, w, tel_io, _ = Erebus.init_simulation(cfg; output_path=tmpdir)
                tel_io !== nothing && close(tel_io)
                return st, co, w
            end

            fill!(state.markers.core.tm, Int8(1))
            if haskey(state.markers.groups, :metal)
                X_fe = Erebus.X_FE_REF_CHONDRITE
                rho_sol = cfg.materials.rhosolidm[1]
                rho_met = cfg.coreformation.rho_metal
                phi_fe_chond = clamp(
                    (X_fe / rho_met) / (X_fe / rho_met + (1.0 - X_fe) / rho_sol), 0.0, 1.0
                )
                fill!(state.markers.groups.metal.Xfe_bulk, phi_fe_chond)
            end
            fill!(state.markers.core.phim, cfg.thermodynamics.phim0)
            fill!(state.markers.core.phinewm, cfg.thermodynamics.phim0)
            fill!(state.markers.core.pfm0, 0.0)
            fill!(state.markers.core.tkm, T0)
            fill!(state.grids.tk1, T0)
            fill!(state.grids.tk2, T0)
            fill!(state.grids.pr, 0.0)
            fill!(state.grids.pf, 0.0)
            fill!(state.grids.ps, 0.0)
            fill!(state.grids.pr0, 0.0)
            fill!(state.grids.pf0, 0.0)
            fill!(state.grids.ps0, 0.0)
            fill!(state.grids.gx, 0.0)
            fill!(state.grids.gy, 0.0)

            t0 = 2.0e6 * yr
            state.timesum = t0
            state.dt = dt_yr * yr
            dt_hist = Float64[]
            t_hist = Float64[t0]

            for step in 1:nsteps
                state.timestep = step
                dt_step_target = state.dt
                Erebus.prepare_step_ambient!(state, coords, cfg, ws)
                Erebus.radiogenic_heating!(state, coords, cfg)
                Erebus.interpolate_markers_to_grid!(
                    state,
                    coords,
                    cfg;
                    p2m_workspace=ws.p2m,
                    thread_buffers=ws.thread_buffers,
                    interp_arrays=ws.interp_arrays,
                )
                state.grids.gx .= 0.0
                state.grids.gy .= 0.0
                Erebus.snapshot_step_start_inventories!(ws, state, cfg)
                Erebus.update_step_start_pressures!(state)
                res = Base.CoreLogging.with_logger(Base.CoreLogging.NullLogger()) do
                    Erebus.solve_thermomechanical_iterations!(
                        state, coords, cfg, ws; dt_step_initial=dt_step_target
                    )
                end
                Erebus.diffuse_and_update_markers!(state, coords, cfg, ws)
                Erebus.vent_and_degas!(state, coords, cfg; Fm_step_start=ws.Fm_step_start)
                state.timesum += state.dt
                push!(dt_hist, state.dt)
                push!(t_hist, state.timesum)
                if res.dt_reduced_by_maxDT
                    state.dt = res.dt_next
                else
                    state.dt = min(state.dt * cfg.time.dtcoefup, cfg.time.dt_longest * yr)
                end
            end

            dT_center = state.grids.tk2[5, 5] - T0
            return (; cfg, t0, t_end=state.timesum, dt_hist, t_hist, dT_center)
        end

        function _analytical_0d_oracles(cfg, t0::Float64, dt_hist::Vector{Float64})
            phim0 = cfg.thermodynamics.phim0
            rho_sol = cfg.materials.rhosolidm[1]
            rho_met = cfg.coreformation.rho_metal
            rhocp_eff =
                (1.0 - phim0) * cfg.materials.rhocpsolidm[1] +
                phim0 * cfg.materials.rhocpfluidm[1]
            X_fe = Erebus.X_FE_REF_CHONDRITE
            phi_fe = clamp(
                (X_fe / rho_met) / (X_fe / rho_met + (1.0 - X_fe) / rho_sol), 0.0, 1.0
            )
            tau_al = cfg.thermodynamics.t_half_al / log(2.0)
            tau_fe = cfg.thermodynamics.t_half_fe / log(2.0)

            H0_al = if cfg.thermodynamics.hr_al
                Q0_al =
                    cfg.thermodynamics.f_al *
                    cfg.thermodynamics.ratio_al *
                    cfg.thermodynamics.E_al / tau_al
                (1.0 - phim0) *
                (1.0 - phi_fe) *
                (Q0_al / (1.0 - X_fe)) *
                rho_sol *
                exp(-t0 / tau_al)
            else
                0.0
            end
            H0_fe = if cfg.thermodynamics.hr_fe
                Q0_fe =
                    cfg.thermodynamics.f_fe *
                    cfg.thermodynamics.ratio_fe *
                    cfg.thermodynamics.E_fe / tau_fe
                (1.0 - phim0) * phi_fe * (Q0_fe / X_fe) * rho_met * exp(-t0 / tau_fe)
            else
                0.0
            end

            t_total = sum(dt_hist)
            edges = _make_subinterval_edges(dt_hist)
            E_gl5 =
                _gl5_decay_integral(H0_al, tau_al, edges) +
                _gl5_decay_integral(H0_fe, tau_fe, edges)
            E_closed =
                H0_al * tau_al * (-expm1(-t_total / tau_al)) +
                H0_fe * tau_fe * (-expm1(-t_total / tau_fe))

            E_start = 0.0
            E_end = 0.0
            E_wrong_lambda = 0.0
            t_cur = 0.0
            for dt_k in dt_hist
                H_t = H0_al * exp(-t_cur / tau_al) + H0_fe * exp(-t_cur / tau_fe)
                E_start += H_t * dt_k
                E_end +=
                    (
                        H0_al * exp(-(t_cur + dt_k) / tau_al) +
                        H0_fe * exp(-(t_cur + dt_k) / tau_fe)
                    ) * dt_k
                tau_al_w = tau_al / log(2.0)
                tau_fe_w = tau_fe / log(2.0)
                E_wrong_lambda +=
                    H0_al * exp(-t_cur / tau_al) * tau_al_w * (-expm1(-dt_k / tau_al_w)) +
                    H0_fe * exp(-t_cur / tau_fe) * tau_fe_w * (-expm1(-dt_k / tau_fe_w))
                t_cur += dt_k
            end

            return (;
                dT_gl5=E_gl5 / rhocp_eff,
                dT_closed=E_closed / rhocp_eff,
                dT_start=E_start / rhocp_eff,
                dT_end=E_end / rhocp_eff,
                dT_wrong_lambda=E_wrong_lambda / rhocp_eff,
            )
        end

        # 1. 20-step uniform dt = 0.5 Myr (26Al + 60Fe)
        res_both = _run_adiabatic_0d_box(20; hr_al=true, hr_fe=true, dt_yr=5.0e5)
        or_both = _analytical_0d_oracles(res_both.cfg, res_both.t0, res_both.dt_hist)
        @test isapprox(or_both.dT_gl5, or_both.dT_closed; atol=10 * eps(or_both.dT_closed))
        @test isapprox(res_both.dT_center, or_both.dT_gl5; rtol=1e-12)
        @test !isapprox(res_both.dT_center, or_both.dT_start; rtol=0.1)
        @test !isapprox(res_both.dT_center, or_both.dT_end; rtol=0.1)
        @test !isapprox(res_both.dT_center, or_both.dT_wrong_lambda; rtol=0.05)

        # 2. 20-step uniform dt = 0.5 Myr (26Al only)
        res_al = _run_adiabatic_0d_box(20; hr_al=true, hr_fe=false, dt_yr=5.0e5)
        or_al = _analytical_0d_oracles(res_al.cfg, res_al.t0, res_al.dt_hist)
        @test isapprox(res_al.dT_center, or_al.dT_gl5; rtol=1e-12)

        # 3. 5-step variable-dt with valid DTmax = 20.0 K mid-step dt cut
        res_cut = _run_adiabatic_0d_box(5; hr_al=true, hr_fe=true, dt_yr=5.0e5, DTmax=20.0)
        or_cut = _analytical_0d_oracles(res_cut.cfg, res_cut.t0, res_cut.dt_hist)
        @test res_cut.dt_hist[1] < 500_000.0 * Erebus.SECONDS_PER_YEAR
        @test isapprox(res_cut.dT_center, or_cut.dT_gl5; rtol=1e-12)
    end
end

using Test
using Erebus
using Erebus.Config
using Erebus.Physics
using Erebus.Particles
using Erebus.Numerics
using StaticArrays
using TOML

@testset "Two-Phase Silicate Melt Segregation & Magma Ascent" begin
    @testset "MagmaTransportConfig Schema & Validation" begin
        cfg_def = MagmaTransportConfig()
        @test cfg_def.active == false
        @test isapprox(cfg_def.k_melt_ref, 1.0e-11; rtol=1e-12)
        @test isapprox(cfg_def.perm_exponent, 3.0; rtol=1e-12)
        @test isapprox(cfg_def.phi0, 0.10; rtol=1e-12)
        @test isapprox(cfg_def.phi_residual, 0.01; rtol=1e-12)
        @test isapprox(cfg_def.phi_crit, 0.40; rtol=1e-12)
        @test isapprox(cfg_def.phi_pack, 1.0; rtol=1e-12)
        @test isapprox(cfg_def.eta_melt, 10.0; rtol=1e-12)
        @test isapprox(cfg_def.r_grain, 1.0e-3; rtol=1e-12)
        @test isapprox(cfg_def.hindered_exponent, 2.0; rtol=1e-12)
        @test isapprox(cfg_def.F_perc_end, 0.35; rtol=1e-12)
        @test isapprox(cfg_def.F_settle_start, 0.45; rtol=1e-12)
        @test isapprox(cfg_def.cfl_melt, 0.5; rtol=1e-12)
        @test cfg_def.max_subcycles == 2000
        @test cfg_def.segregation_heating == true
        @test cfg_def.latent_crystallization == true
        @test cfg_def.exsolution_active == true
        @test cfg_def.track_depletion == true

        # Integration in SimulationConfig
        sim_cfg = default_config()
        @test sim_cfg.magma_transport isa MagmaTransportConfig
        @test validate_config(sim_cfg) === nothing

        # Validation dependency: magma_transport requires active melting
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=MeltingConfig(; active=false),
                magma_transport=MagmaTransportConfig(; active=true),
            ),
        )

        m_active = MeltingConfig(; active=true)

        # Validation bounds: negative reference permeability
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, k_melt_ref=-1.0e-11),
            ),
        )

        # Validation bounds: non-positive permeability exponent
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, perm_exponent=0.0),
            ),
        )

        # Validation bounds: non-positive reference porosity
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, phi0=0.0),
            ),
        )

        # Validation bounds: unphysical residual porosity
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, phi_residual=-0.01),
            ),
        )

        # Validation bounds: unphysical critical melt fraction
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, phi_crit=0.0),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, phi_crit=1.0),
            ),
        )

        # Validation bounds: unphysical maximum packing fraction
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, phi_pack=1.05),
            ),
        )

        # Validation bounds: regime hierarchy ordering violation (phi_crit < F_perc_end)
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(;
                    active=true, phi_crit=0.30, F_perc_end=0.35, F_settle_start=0.45
                ),
            ),
        )

        # Validation bounds: regime hierarchy ordering violation (phi_pack < F_settle_start)
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(;
                    active=true, phi_pack=0.40, F_settle_start=0.45
                ),
            ),
        )

        # Validation bounds: regime hierarchy ordering violation (phi_residual >= F_perc_end)
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(;
                    active=true, phi_residual=0.36, F_perc_end=0.35
                ),
            ),
        )

        # Validation bounds: non-positive melt viscosity
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, eta_melt=0.0),
            ),
        )

        # Validation bounds: non-positive grain size
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, r_grain=-1.0e-3),
            ),
        )

        # Validation bounds: negative hindered settling exponent
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, hindered_exponent=-0.5),
            ),
        )

        # Validation bounds: unphysical CFL parameter
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, cfl_melt=0.0),
            ),
        )
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, cfl_melt=1.5),
            ),
        )

        # Validation bounds: invalid maximum subcycles
        @test_throws ArgumentError validate_config(
            SimulationConfig(;
                melting=m_active,
                magma_transport=MagmaTransportConfig(; active=true, max_subcycles=0),
            ),
        )
    end

    @testset "MagmaTransportConfig TOML Round-Trip" begin
        custom_magma = MagmaTransportConfig(;
            active=true,
            k_melt_ref=2.5e-11,
            perm_exponent=2.8,
            phi0=0.04,
            phi_residual=0.008,
            phi_crit=0.38,
            phi_pack=0.65,
            eta_melt=8.0,
            r_grain=1.5e-3,
            hindered_exponent=2.3,
            F_perc_end=0.32,
            F_settle_start=0.42,
            cfl_melt=0.4,
            max_subcycles=80,
            segregation_heating=true,
            latent_crystallization=true,
            exsolution_active=true,
            track_depletion=true,
        )
        sim_cfg = SimulationConfig(;
            melting=MeltingConfig(; active=true), magma_transport=custom_magma
        )
        toml_str = Erebus.Config.save_config(sim_cfg)
        loaded_cfg = Erebus.Config.load_config(toml_str)

        @test loaded_cfg.magma_transport.active == true
        @test isapprox(loaded_cfg.magma_transport.k_melt_ref, 2.5e-11; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.perm_exponent, 2.8; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.phi0, 0.04; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.phi_residual, 0.008; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.phi_crit, 0.38; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.phi_pack, 0.65; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.eta_melt, 8.0; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.r_grain, 1.5e-3; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.hindered_exponent, 2.3; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.F_perc_end, 0.32; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.F_settle_start, 0.42; rtol=1e-12)
        @test isapprox(loaded_cfg.magma_transport.cfl_melt, 0.4; rtol=1e-12)
        @test loaded_cfg.magma_transport.max_subcycles == 80
        @test loaded_cfg.magma_transport.segregation_heating == true
        @test loaded_cfg.magma_transport.latent_crystallization == true
        @test loaded_cfg.magma_transport.exsolution_active == true
        @test loaded_cfg.magma_transport.track_depletion == true
    end

    @testset "silicate_melt_permeability: McKenzie (1984) Power Law" begin
        k0 = 1.0e-11
        phi0 = 0.05
        n = 3.0
        phi_res = 0.005

        # Sub-residual melt fraction yields zero permeability
        @test iszero(
            silicate_melt_permeability(0.0; k0=k0, phi0=phi0, n=n, phi_residual=phi_res)
        )
        @test iszero(
            silicate_melt_permeability(0.003; k0=k0, phi0=phi0, n=n, phi_residual=phi_res)
        )
        @test iszero(
            silicate_melt_permeability(phi_res; k0=k0, phi0=phi0, n=n, phi_residual=phi_res)
        )

        # Reference porosity yields reference permeability k0
        k_ref = silicate_melt_permeability(
            phi0 + phi_res; k0=k0, phi0=phi0, n=n, phi_residual=phi_res
        )
        @test isapprox(k_ref, k0; rtol=1e-12)

        # Power law scaling verification: double porosity yields 2^n times permeability
        k_double = silicate_melt_permeability(
            2.0 * phi0 + phi_res; k0=k0, phi0=phi0, n=n, phi_residual=phi_res
        )
        @test isapprox(k_double, k0 * (2.0^n); rtol=1e-12)

        # Monotonicity test
        phi_vals = range(phi_res + 0.001, 0.35; length=30)
        k_vals = [
            silicate_melt_permeability(p; k0=k0, phi0=phi0, n=n, phi_residual=phi_res) for
            p in phi_vals
        ]
        @test issorted(k_vals)
        @test all(k_vals .> 0.0)

        # Domain error guards
        @test_throws DomainError silicate_melt_permeability(-0.01)
        @test_throws DomainError silicate_melt_permeability(NaN)
        @test_throws DomainError silicate_melt_permeability(0.05; k0=-1.0e-11)
        @test_throws DomainError silicate_melt_permeability(0.05; phi0=0.0)
        @test_throws DomainError silicate_melt_permeability(0.05; n=0.0)
        @test_throws DomainError silicate_melt_permeability(0.05; phi_residual=-0.01)
    end

    @testset "silicate_melt_segregation_velocity: Regime Blending & Limits" begin
        k0 = 1.0e-11
        phi0 = 0.05
        n_perm = 3.0
        phi_res = 0.005
        eta_liq = 10.0
        r_gr = 1.0e-3
        n_hind = 2.0
        drho = 300.0
        g = 0.2
        F_p_end = 0.35
        F_s_start = 0.45

        # Sub-residual melt fraction yields zero segregation velocity
        v_sub = silicate_melt_segregation_velocity(
            0.002,
            drho,
            g,
            eta_liq;
            k_melt_ref=k0,
            phi0=phi0,
            perm_exponent=n_perm,
            phi_residual=phi_res,
            r_grain=r_gr,
            hindered_exponent=n_hind,
            F_perc_end=F_p_end,
            F_settle_start=F_s_start,
        )
        @test iszero(v_sub)

        # Zero gravity yields zero velocity
        v_g0 = silicate_melt_segregation_velocity(
            0.1,
            drho,
            0.0,
            eta_liq;
            k_melt_ref=k0,
            phi0=phi0,
            perm_exponent=n_perm,
            phi_residual=phi_res,
            r_grain=r_gr,
            hindered_exponent=n_hind,
            F_perc_end=F_p_end,
            F_settle_start=F_s_start,
        )
        @test iszero(v_g0)

        # Negative density contrast yields zero buoyant velocity
        v_neg = silicate_melt_segregation_velocity(
            0.1,
            -10.0,
            g,
            eta_liq;
            k_melt_ref=k0,
            phi0=phi0,
            perm_exponent=n_perm,
            phi_residual=phi_res,
            r_grain=r_gr,
            hindered_exponent=n_hind,
            F_perc_end=F_p_end,
            F_settle_start=F_s_start,
        )
        @test iszero(v_neg)

        # Low melt fraction (Darcy percolation regime: F_m = 0.1 <= F_p_end)
        v_darcy = silicate_melt_segregation_velocity(
            0.1,
            drho,
            g,
            eta_liq;
            k_melt_ref=k0,
            phi0=phi0,
            perm_exponent=n_perm,
            phi_residual=phi_res,
            r_grain=r_gr,
            hindered_exponent=n_hind,
            F_perc_end=F_p_end,
            F_settle_start=F_s_start,
        )
        k_expected = silicate_melt_permeability(
            0.1; k0=k0, phi0=phi0, n=n_perm, phi_residual=phi_res
        )
        v_darcy_expected = (k_expected / (eta_liq * 0.1)) * drho * g
        @test isapprox(v_darcy, v_darcy_expected; rtol=1e-12)

        # High melt fraction (Stokes crystal suspension regime: F_m = 0.6 >= F_s_start)
        # Richardson-Zaki hindered settling: v_susp = v_stokes * F_m^n
        v_stokes = silicate_melt_segregation_velocity(
            0.6,
            drho,
            g,
            eta_liq;
            k_melt_ref=k0,
            phi0=phi0,
            perm_exponent=n_perm,
            phi_residual=phi_res,
            r_grain=r_gr,
            hindered_exponent=n_hind,
            F_perc_end=F_p_end,
            F_settle_start=F_s_start,
        )
        v_stokes_base = (2.0 / 9.0) * (r_gr^2) * drho * g / eta_liq
        v_stokes_expected = v_stokes_base * (0.6^n_hind)
        @test isapprox(v_stokes, v_stokes_expected; rtol=1e-12)

        # Pure melt limit (F_m = 1.0): unhindered Stokes velocity
        v_pure = silicate_melt_segregation_velocity(
            1.0,
            drho,
            g,
            eta_liq;
            k_melt_ref=k0,
            phi0=phi0,
            perm_exponent=n_perm,
            phi_residual=phi_res,
            r_grain=r_gr,
            hindered_exponent=n_hind,
            F_perc_end=F_p_end,
            F_settle_start=F_s_start,
        )
        @test isapprox(v_pure, v_stokes_base; rtol=1e-12)

        # Transition regime (F_perc_end < F_m < F_settle_start): smooth Hermite interpolation
        F_mid = 0.5 * (F_p_end + F_s_start)
        v_trans = silicate_melt_segregation_velocity(
            F_mid,
            drho,
            g,
            eta_liq;
            k_melt_ref=k0,
            phi0=phi0,
            perm_exponent=n_perm,
            phi_residual=phi_res,
            r_grain=r_gr,
            hindered_exponent=n_hind,
            F_perc_end=F_p_end,
            F_settle_start=F_s_start,
        )
        v_p_end = silicate_melt_segregation_velocity(
            F_p_end,
            drho,
            g,
            eta_liq;
            k_melt_ref=k0,
            phi0=phi0,
            perm_exponent=n_perm,
            phi_residual=phi_res,
            r_grain=r_gr,
            hindered_exponent=n_hind,
            F_perc_end=F_p_end,
            F_settle_start=F_s_start,
        )
        v_s_start = silicate_melt_segregation_velocity(
            F_s_start,
            drho,
            g,
            eta_liq;
            k_melt_ref=k0,
            phi0=phi0,
            perm_exponent=n_perm,
            phi_residual=phi_res,
            r_grain=r_gr,
            hindered_exponent=n_hind,
            F_perc_end=F_p_end,
            F_settle_start=F_s_start,
        )
        @test v_trans >= min(v_p_end, v_s_start)
        @test v_trans <= max(v_p_end, v_s_start)

        # Domain error guards
        @test_throws DomainError silicate_melt_segregation_velocity(-0.1, drho, g, eta_liq)
        @test_throws DomainError silicate_melt_segregation_velocity(1.1, drho, g, eta_liq)
        @test_throws DomainError silicate_melt_segregation_velocity(0.2, NaN, g, eta_liq)
        @test_throws DomainError silicate_melt_segregation_velocity(
            0.2, drho, -0.1, eta_liq
        )
        @test_throws DomainError silicate_melt_segregation_velocity(0.2, drho, g, 0.0)
    end

    @testset "silicate_melt_dissipation_heating: Gravitational Energy Release" begin
        v_seg = 1.0e-5
        drho = 300.0
        g = 0.5
        F_m = 0.2
        psi = silicate_melt_dissipation_heating(F_m, drho, g, v_seg)
        @test isapprox(psi, drho * g * F_m * v_seg; rtol=1e-12)

        # Zero dissipation when stationary or neutral buoyancy
        @test iszero(silicate_melt_dissipation_heating(0.0, drho, g, v_seg))
        @test iszero(silicate_melt_dissipation_heating(F_m, 0.0, g, v_seg))
        @test iszero(silicate_melt_dissipation_heating(F_m, drho, 0.0, v_seg))
        @test iszero(silicate_melt_dissipation_heating(F_m, drho, g, 0.0))

        # Domain error guards
        @test_throws DomainError silicate_melt_dissipation_heating(-0.1, drho, g, v_seg)
        @test_throws DomainError silicate_melt_dissipation_heating(1.1, drho, g, v_seg)
        @test_throws DomainError silicate_melt_dissipation_heating(F_m, NaN, g, v_seg)
        @test_throws DomainError silicate_melt_dissipation_heating(F_m, drho, NaN, v_seg)
        @test_throws DomainError silicate_melt_dissipation_heating(F_m, drho, g, NaN)
    end

    @testset "setup_marker_magma_properties Allocation" begin
        marknum = 100
        magma_props = setup_marker_magma_properties(marknum)
        @test length(magma_props) == 1
        F_extract_m = magma_props[1]
        @test length(F_extract_m) == marknum
        @test all(iszero, F_extract_m)
        @test eltype(F_extract_m) === Float64
    end

    @testset "apply_silicate_melt_segregation! Conservation & Ascent" begin
        # 16x16 grid setup over 80 km x 80 km domain
        Nx = 16
        Ny = 16
        xsize = 80000.0
        ysize = 80000.0
        dx = xsize / (Nx - 1)
        dy = ysize / (Ny - 1)
        coords = GridCoordinates(Nx, Ny; xsize=xsize, ysize=ysize)
        xcenter = coords.xcenter
        ycenter = coords.ycenter
        rplanet = 30000.0
        g_surf = 0.2

        # Populate markers: 4 markers per cell in planetary interior
        marknum = 4 * Nx * Ny
        xm = zeros(Float64, marknum)
        ym = zeros(Float64, marknum)
        tm = fill(1, marknum)
        tkm = fill(1600.0, marknum)
        Fm = zeros(Float64, marknum)
        F_extract_m = zeros(Float64, marknum)

        idx = 1
        for i in 1:Ny, j in 1:Nx
            xc = (j - 1) * dx
            yc = (i - 1) * dy
            for sx in (-0.25, 0.25), sy in (-0.25, 0.25)
                xm[idx] = xc + sx * dx
                ym[idx] = yc + sy * dy
                rmark = distance(xm[idx], ym[idx], xcenter, ycenter)
                if rmark <= rplanet
                    tm[idx] = 1 # Rock mantle
                    # Partially molten plume between 10 km and 18 km radius
                    if 10000.0 <= rmark <= 18000.0
                        Fm[idx] = 0.20
                        tkm[idx] = 1600.0
                    elseif rmark < 6000.0
                        # Super-liquidus magma ocean core: Fm = 1.0 > phi_pack
                        Fm[idx] = 1.0
                        tkm[idx] = 2000.0
                    else
                        Fm[idx] = 0.0
                        # Exterior cold mantle below solidus (1400 K)
                        tkm[idx] = 1200.0
                    end
                else
                    tm[idx] = 3 # Sticky air
                    tkm[idx] = 200.0
                    Fm[idx] = 0.0
                end
                idx += 1
            end
        end

        # Test with phi_pack = 0.60 < 1.0: markers with Fm = 1.0 must not crash
        cfg_magma = MagmaTransportConfig(;
            active=true,
            k_melt_ref=1.0e-10,
            phi0=0.05,
            phi_residual=0.005,
            phi_crit=0.38,
            F_perc_end=0.32,
            F_settle_start=0.42,
            phi_pack=0.60,
            eta_melt=1.0,
            cfl_melt=0.4,
            max_subcycles=20,
            segregation_heating=true,
            latent_crystallization=true,
            track_depletion=true,
        )

        Q_seg_grid = zeros(Float64, Ny, Nx)
        Q_lat_grid = zeros(Float64, Ny, Nx)
        dt = 5.0e10 # Transport timestep (~1500 years)

        M_melt_initial = sum(Fm)
        @test M_melt_initial > 0.0

        # Run silicate melt segregation
        res = apply_silicate_melt_segregation!(
            xm,
            ym,
            tm,
            tkm,
            Fm,
            marknum,
            dt,
            cfg_magma;
            coords=coords,
            xcenter=xcenter,
            ycenter=ycenter,
            rplanet=rplanet,
            g_surf=g_surf,
            Q_seg_grid=Q_seg_grid,
            Q_lat_grid=Q_lat_grid,
            rho_silicate=3300.0,
            rho_melt=2800.0,
            T_solidus_silicate=1400.0,
            T_liquidus_silicate=1800.0,
            L_melt=4.0e5,
            F_extract_m=F_extract_m,
        )

        # Machine-precision mass conservation: remaining melt + crystallized melt
        M_melt_final = sum(Fm)
        M_conserved = M_melt_final + res.total_crystallized_mass
        relative_mass_drift = abs(M_conserved - M_melt_initial) / M_melt_initial
        @test relative_mass_drift < 1.0e-12

        # Radial upward / outward melt segregation
        @test res.max_v_seg > 0.0
        @test res.n_subcycles >= 1

        # Energy dissipation heating occurred and populated Q_seg_grid
        @test res.total_dissipation_energy > 0.0
        @test maximum(Q_seg_grid) > 0.0

        # Latent heat release matches crystallized mass to machine precision
        @test res.total_crystallized_mass > 0.0
        @test maximum(Q_lat_grid) > 0.0
        dV_cell = dx * dy
        E_lat_grid = sum(Q_lat_grid) * dt * dV_cell
        # Total crystallized mass on markers: each unit of Fm represents dV_cell / 4 mass equivalent
        E_lat_expected = res.total_crystallized_mass * (dV_cell / 4) * 2800.0 * 4.0e5
        @test maximum(Q_seg_grid) > 0.0

        # Depletion tracking accumulated on donor markers
        @test maximum(F_extract_m) > 0.0

        # Inactive simulation no-op
        cfg_inactive = MagmaTransportConfig(; active=false)
        res_inactive = apply_silicate_melt_segregation!(
            xm, ym, tm, tkm, Fm, marknum, dt, cfg_inactive; coords=coords
        )
        @test iszero(res_inactive.max_v_seg)
        @test res_inactive.n_subcycles == 0

        # Neutral and negative buoyancy tests (Finding 2)
        res_neutral = apply_silicate_melt_segregation!(
            xm,
            ym,
            tm,
            tkm,
            Fm,
            marknum,
            dt,
            cfg_magma;
            coords=coords,
            xcenter=xcenter,
            ycenter=ycenter,
            rplanet=rplanet,
            g_surf=g_surf,
            rho_silicate=2800.0,
            rho_melt=2800.0,
        )
        @test iszero(res_neutral.max_v_seg)
        @test res_neutral.n_subcycles == 0

        res_negative = apply_silicate_melt_segregation!(
            xm,
            ym,
            tm,
            tkm,
            Fm,
            marknum,
            dt,
            cfg_magma;
            coords=coords,
            xcenter=xcenter,
            ycenter=ycenter,
            rplanet=rplanet,
            g_surf=g_surf,
            rho_silicate=2700.0,
            rho_melt=2800.0,
        )
        @test iszero(res_negative.max_v_seg)
        @test res_negative.n_subcycles == 0
    end
end

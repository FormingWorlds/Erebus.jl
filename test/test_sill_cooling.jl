using Test
using LinearAlgebra
using SparseArrays
using Random
using Erebus
using Erebus.Config
using Erebus.Numerics
using Erebus.Physics
using Erebus.Particles
using Erebus.Geometry
using ExtendableSparse

# Standard erf implementation (Abramowitz & Stegun 7.1.26, max error < 1.5e-7)
function _sill_erf(x::Real)
    p = 0.3275911
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    sign_x = sign(x)
    abs_x = abs(x)
    t = 1.0 / (1.0 + p * abs_x)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-abs_x * abs_x)
    return sign_x * y
end

# Jaeger (1957) analytical solution for conductive cooling of an intrusive sheet / sill
function jaeger_sill_temperature(y::Real, t::Real, b::Real, T0::Real, Tc::Real, kappa::Real)
    if t <= 0.0
        return abs(y) <= b ? T0 : Tc
    end
    diff_term = 2.0 * sqrt(kappa * t)
    term1 = _sill_erf((b - y) / diff_term)
    term2 = _sill_erf((b + y) / diff_term)
    return Tc + 0.5 * (T0 - Tc) * (term1 + term2)
end

@testset "Crustal Sill Cooling & Magma-Hydrothermal Coupling" begin
    @testset "Isothermal Zero Sensible Heating Discrimination" begin
        coords = GridCoordinates(11, 11; xsize=10000.0, ysize=10000.0)
        marknum = 200
        rng = MersenneTwister(42)

        # Markers distributed deterministically with small jitter
        xm = [
            clamp(
                coords.dx * (0.5 + mod(m - 1, 9) + 0.1 * (rand(rng) - 0.5)), 100.0, 9900.0
            ) for m in 1:marknum
        ]
        ym = [
            clamp(
                coords.dy * (0.5 + div(m - 1, 9) + 0.1 * (rand(rng) - 0.5)), 100.0, 9900.0
            ) for m in 1:marknum
        ]
        tm = ones(Int32, marknum)
        tkm = fill(1500.0, marknum) # Exact isothermal domain
        Fm = [clamp(0.10 + 0.05 * sin(m), 0.0, 1.0) for m in 1:marknum]

        cfg_magma = MagmaTransportConfig(;
            active=true,
            sensible_heat_transport=true,
            cp_melt=1200.0,
            segregation_heating=false,
            latent_crystallization=false,
            compaction_active=false,
            phi_residual=0.01,
        )

        Q_seg_grid = zeros(Float64, coords.Ny1, coords.Nx1)
        dt_val = 1.0e6

        res = apply_silicate_melt_segregation!(
            xm,
            ym,
            tm,
            tkm,
            Fm,
            marknum,
            dt_val,
            cfg_magma;
            coords=coords,
            xcenter=5000.0,
            ycenter=5000.0,
            rplanet=5000.0,
            g_surf=0.5,
            Q_seg_grid=Q_seg_grid,
        )

        @test res.max_v_seg > 0.0
        # In an isothermal domain, temperature advection must be identically zero
        @test isapprox(maximum(abs.(Q_seg_grid)), 0.0; atol=1e-12)
    end

    @testset "Sensible Heat Advection in Thermal Gradient" begin
        coords = GridCoordinates(11, 11; xsize=10000.0, ysize=10000.0)
        Ny, Nx = coords.Ny, coords.Nx
        Ny1, Nx1 = coords.Ny1, coords.Nx1
        rng = MersenneTwister(123)

        marknum = 400
        xm = [
            clamp(
                coords.dx * (0.5 + mod(m - 1, 9) + 0.1 * (rand(rng) - 0.5)), 100.0, 9900.0
            ) for m in 1:marknum
        ]
        ym = [
            clamp(
                coords.dy * (0.5 + div(m - 1, 9) + 0.1 * (rand(rng) - 0.5)), 100.0, 9900.0
            ) for m in 1:marknum
        ]
        tm = ones(Int32, marknum)
        # Hot bottom, cooler top
        tkm = [1600.0 - 600.0 * (ym[m] / 10000.0) for m in 1:marknum]
        Fm = [clamp(0.20 - 0.10 * (ym[m] / 10000.0), 0.0, 1.0) for m in 1:marknum]

        cfg_magma = MagmaTransportConfig(;
            active=true,
            sensible_heat_transport=true,
            cp_melt=1200.0,
            segregation_heating=false,
            latent_crystallization=false,
            compaction_active=false,
            phi_residual=0.01,
            max_subcycles=50,
        )

        ws = MagmaSegregationWorkspace(Ny, Nx)
        Q_seg_grid = zeros(Ny1, Nx1)
        dt_val = 1.0e6

        res = apply_silicate_melt_segregation!(
            xm,
            ym,
            tm,
            tkm,
            Fm,
            marknum,
            dt_val,
            cfg_magma;
            coords=coords,
            xcenter=5000.0,
            ycenter=5000.0,
            rplanet=5000.0,
            g_surf=0.5,
            Q_seg_grid=Q_seg_grid,
            workspace=ws,
        )

        @test res.max_v_seg > 0.0
        # Hot rising melt deposits positive sensible heat in cooler crust
        @test maximum(Q_seg_grid) > 0.0
        @test res.total_sensible_energy > 0.0
        @test isfinite(res.total_sensible_energy)
        @test res.dt_sub > 0.0

        # All deposited heat resides on interior nodes (never on ghost rim)
        @test isapprox(sum(abs.(Q_seg_grid[1, :])), 0.0; atol=1e-12)
        @test isapprox(sum(abs.(Q_seg_grid[Ny1, :])), 0.0; atol=1e-12)
        @test isapprox(sum(abs.(Q_seg_grid[:, 1])), 0.0; atol=1e-12)
        @test isapprox(sum(abs.(Q_seg_grid[:, Nx1])), 0.0; atol=1e-12)
    end

    @testset "Crustal Sill Solidification & Latent Heat Kinetics" begin
        coords = GridCoordinates(8, 8; xsize=8000.0, ysize=8000.0)
        Ny, Nx = coords.Ny, coords.Nx
        Ny1, Nx1 = coords.Ny1, coords.Nx1

        marknum = 640
        rng = MersenneTwister(999)
        xm = [
            clamp(
                coords.dx * (0.5 + mod(m - 1, 7) + 0.1 * (rand(rng) - 0.5)), 100.0, 7900.0
            ) for m in 1:marknum
        ]
        ym = [
            clamp(
                coords.dy * (0.5 + div(m - 1, 7) + 0.1 * (rand(rng) - 0.5)), 100.0, 7900.0
            ) for m in 1:marknum
        ]
        tm = ones(Int32, marknum)
        tkm = zeros(marknum)
        Fm_inst = zeros(marknum)

        for m in 1:marknum
            if ym[m] < 2000.0
                tkm[m] = 800.0
                Fm_inst[m] = 0.0
            elseif ym[m] <= 4000.0
                tkm[m] = 1200.0
                Fm_inst[m] = 0.20
            else
                tkm[m] = 1500.0
                Fm_inst[m] = 0.25
            end
        end

        Fm_kin = copy(Fm_inst)
        initial_melt_mass = sum(Fm_inst)

        cfg_inst = MagmaTransportConfig(;
            active=true,
            latent_crystallization=true,
            sill_cooling_active=true,
            crystallization_timescale=0.0,
            sensible_heat_transport=false,
            segregation_heating=false,
            compaction_active=false,
            phi_residual=0.01,
        )

        ws = MagmaSegregationWorkspace(Ny, Nx)
        Q_lat_inst = zeros(Ny1, Nx1)
        dt_val = 1.0e6

        res_inst = apply_silicate_melt_segregation!(
            xm,
            ym,
            tm,
            tkm,
            Fm_inst,
            marknum,
            dt_val,
            cfg_inst;
            coords=coords,
            xcenter=4000.0,
            ycenter=4000.0,
            rplanet=4000.0,
            g_surf=0.5,
            Q_lat_grid=Q_lat_inst,
            T_solidus_silicate=1400.0,
            T_liquidus_silicate=1800.0,
            L_melt=4.0e5,
            workspace=ws,
        )

        @test res_inst.total_crystallized_mass > 0.0
        @test maximum(Q_lat_inst) > 0.0
        @test isfinite(res_inst.total_crystallized_mass)

        # Kinetic freezing with finite timescale tau = 1.0e7 s
        cfg_kin = MagmaTransportConfig(;
            active=true,
            latent_crystallization=true,
            sill_cooling_active=true,
            crystallization_timescale=1.0e7,
            sensible_heat_transport=false,
            segregation_heating=false,
            compaction_active=false,
            phi_residual=0.01,
        )

        Q_lat_kin = zeros(Ny1, Nx1)
        res_kin = apply_silicate_melt_segregation!(
            xm,
            ym,
            tm,
            tkm,
            Fm_kin,
            marknum,
            dt_val,
            cfg_kin;
            coords=coords,
            xcenter=4000.0,
            ycenter=4000.0,
            rplanet=4000.0,
            g_surf=0.5,
            Q_lat_grid=Q_lat_kin,
            T_solidus_silicate=1400.0,
            T_liquidus_silicate=1800.0,
            L_melt=4.0e5,
            workspace=ws,
        )

        @test sum(Fm_kin) < initial_melt_mass
        @test sum(Fm_kin) > 0.0
        @test res_kin.total_crystallized_mass < res_inst.total_crystallized_mass
        @test maximum(Q_lat_kin) < maximum(Q_lat_inst)
    end

    @testset "Hydrothermal Convection Enhancement Above Emplaced Sill" begin
        cfg_hydro = HydrothermalConfig(;
            active=true,
            phi_start=0.20,
            phi_end=0.60,
            kphi_ref=1.0e-12,
            H_layer=5000.0,
            dT_min=5.0,
            gravity=0.5,
        )

        k_cond = 2.5 # W/(m K)
        phi_crust = 0.35
        tm_crust = 1

        # Above cooling sill: hot hydrothermal fluid zone (T = 450 K)
        T_sill_plume = 450.0
        k_eff_plume = apply_hydrothermal_convection_closure(
            k_cond, T_sill_plume, phi_crust, tm_crust; cfg=cfg_hydro, H_eff=5000.0
        )

        # Sub-critical / cold crust: below surface reference (T = 270 K <= 273.15 K)
        T_cold = 270.0
        k_eff_cold = apply_hydrothermal_convection_closure(
            k_cond, T_cold, phi_crust, tm_crust; cfg=cfg_hydro, H_eff=5000.0
        )

        Nu_plume = k_eff_plume / k_cond
        Nu_cold = k_eff_cold / k_cond
        @test Nu_plume > 1.5
        @test isapprox(Nu_cold, 1.0; atol=1e-12)
        @test k_eff_plume > k_eff_cold

        # Sill-bounded layer thickness scaling: thicker permeable zone yields stronger convection
        k_eff_thin = apply_hydrothermal_convection_closure(
            k_cond, T_sill_plume, phi_crust, tm_crust; cfg=cfg_hydro, H_eff=2000.0
        )
        @test k_eff_plume > k_eff_thin > k_cond
    end

    @testset "Jaeger (1957) 1D Analytical Sill Cooling Benchmark" begin
        # Sheet sill: half-thickness b = 50 m, T0 = 1400 K, Tc = 400 K
        b = 50.0
        T0 = 1400.0
        Tc = 400.0
        rho = 2800.0
        cp = 1000.0
        k_therm = 2.8
        kappa = k_therm / (rho * cp) # 1.0e-6 m^2/s

        # 1D finite-difference diffusion grid across y in [-200, 200] m
        L = 400.0
        N = 201
        coords = GridCoordinates(3, N; xsize=40.0, ysize=L)
        Ny1 = coords.Ny1
        Nx1 = coords.Nx1

        tk1 = fill(Tc, Ny1, Nx1)
        for i in 2:(Ny1 - 1)
            yj = -200.0 + (i - 2) * coords.dy
            if abs(yj) < b
                for j in 1:Nx1
                    tk1[i, j] = T0
                end
            elseif isapprox(abs(yj), b; atol=0.25 * coords.dy)
                for j in 1:Nx1
                    tk1[i, j] = 0.5 * (T0 + Tc)
                end
            end
        end

        RHOCP = fill(rho * cp, Ny1, Nx1)
        KX = fill(k_therm, Ny1, coords.Nx)
        KY = fill(k_therm, coords.Ny, Nx1)
        HR = zeros(Ny1, Nx1)
        HA = zeros(Ny1, Nx1)
        HS = zeros(Ny1, Nx1)
        DHP = zeros(Ny1, Nx1)
        RT = zeros(Ny1 * Nx1)

        t_total = 5.0e8 # ~15.8 yr (diffusion length 2*sqrt(kappa*t) ≈ 44.7 m, resolved across 22.4 cells)
        dt = 5.0e6
        steps = Int(round(t_total / dt))
        LT = ExtendableSparseMatrix(Ny1 * Nx1, Ny1 * Nx1)

        # Drive the actual Erebus production thermal solver
        for _ in 1:steps
            assemble_thermal_lse!(
                tk1, RHOCP, KX, KY, HR, HA, HS, DHP, RT, dt; coords=coords, LT=LT
            )
            tk2_vec = LT \ RT
            tk1 .= reshape(tk2_vec, Ny1, Nx1)
        end

        y_interior = [-200.0 + (i - 2) * coords.dy for i in 2:(Ny1 - 1)]
        T_num = tk1[2:(Ny1 - 1), 2]
        T_ana = [
            jaeger_sill_temperature(yj, t_total, b, T0, Tc, kappa) for yj in y_interior
        ]

        l2_err = norm(T_num - T_ana) / norm(T_ana)
        max_err = maximum(abs.(T_num - T_ana))

        # L2 relative error must be under 0.1%, max absolute error under 1.0 K
        @test l2_err < 0.001
        @test max_err < 1.0

        # Centerline temperature verification
        T_center_ana = Tc + (T0 - Tc) * _sill_erf(b / (2.0 * sqrt(kappa * t_total)))
        center_idx = (length(y_interior) + 1) ÷ 2
        @test isapprox(T_num[center_idx], T_center_ana; atol=1.0)
    end

    @testset "Stefan Latent Heat Buffering Benchmark" begin
        T_sol = 1300.0
        T_liq = 1500.0
        L_melt = 4.0e5 # J/kg
        cp = 1000.0
        rho = 2800.0
        rhocp = rho * cp
        Ste = cp * (T_liq - T_sol) / L_melt # 0.5

        # In mushy interval (1400 K), rhocp_apparent_silicate must provide (1 + 1/Ste) buffering
        T_mush = 1400.0
        rhocp_buffered = rhocp_apparent_silicate(
            T_mush, 0.0, rhocp, rho, 1; T_sol=T_sol, T_liq=T_liq, L_melt=L_melt, active=true
        )
        buffering_factor = rhocp_buffered / rhocp
        @test isapprox(buffering_factor, 1.0 + 1.0 / Ste; rtol=1e-12)

        # Outside mushy interval: identical to baseline sensible heat capacity
        rhocp_subsolidus = rhocp_apparent_silicate(
            1200.0, 0.0, rhocp, rho, 1; T_sol=T_sol, T_liq=T_liq, L_melt=L_melt, active=true
        )
        @test isapprox(rhocp_subsolidus, rhocp; rtol=1e-12)

        rhocp_superliquidus = rhocp_apparent_silicate(
            1600.0, 0.0, rhocp, rho, 1; T_sol=T_sol, T_liq=T_liq, L_melt=L_melt, active=true
        )
        @test isapprox(rhocp_superliquidus, rhocp; rtol=1e-12)
    end

    @testset "Contact Metamorphic Dehydration Aureole Benchmark" begin
        T_country = 400.0
        T_sill = 1400.0

        # Thermodynamic reaction equilibrium threshold
        cfg_react = ReactionConfig(; delta_H=2.4e5, delta_S=370.0)
        T_dehydration_thermo = cfg_react.delta_H / cfg_react.delta_S # ~648.6 K
        T_contact = 0.5 * (T_sill + T_country)
        @test T_contact > T_dehydration_thermo

        # Solve hydromechanical overpressure pulse driven by dehydration DQPF
        coords = GridCoordinates(7, 7; xsize=14000.0, ysize=14000.0)
        Ny, Nx = coords.Ny, coords.Nx
        Ny1, Nx1 = coords.Ny1, coords.Nx1

        ETA = fill(1.0e21, Ny, Nx)
        ETAP = fill(1.0e21, Ny1, Nx1)
        GGG = fill(1.0e10, Ny, Nx)
        GGGP = fill(1.0e10, Ny1, Nx1)
        SXY0 = zeros(Ny, Nx)
        SXX0 = zeros(Ny, Nx)
        RHOX = fill(2800.0, Ny1, Nx1)
        RHOY = fill(2800.0, Ny1, Nx1)
        RHOFX = fill(1000.0, Ny1, Nx1)
        RHOFY = fill(1000.0, Ny1, Nx1)
        RX = fill(1.0e18, Ny1, Nx1)
        RY = fill(1.0e18, Ny1, Nx1)
        PHI = fill(0.05, Ny1, Nx1)
        ETAPHI = fill(1.0e24, Ny1, Nx1)
        BETAPHI = fill(1.0e-10, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = zeros(Ny1, Nx1)
        pr0 = zeros(Ny1, Nx1)
        pf0 = zeros(Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)
        dt_val = 1.0e7

        # Dehydration source term in contact aureole (cells adjacent to sill)
        DQPF = zeros(Ny1, Nx1)
        DQPF[3, 4] = 1.0e-12 # fluid mass production rate [1/s]
        DQPF[4, 4] = 1.0e-12

        R_6 = zeros(Nx1 * Ny1 * 6)
        L_6 = assemble_hydromechanical_lse!(
            ETA,
            ETAP,
            GGG,
            GGGP,
            SXY0,
            SXX0,
            RHOX,
            RHOY,
            RHOFX,
            RHOFY,
            RX,
            RY,
            ETAPHI,
            BETAPHI,
            PHI,
            gx,
            gy,
            pr0,
            pf0,
            DMP,
            dt_val,
            R_6;
            coords=coords,
            DQPF=DQPF,
            fluid_overpressure_coupling=true,
        )

        S_6 = L_6 \ R_6
        pr_6 = zeros(Ny1, Nx1)
        pf_6 = zeros(Ny1, Nx1)
        vx_6 = zeros(Ny1, Nx1)
        vy_6 = zeros(Ny1, Nx1)
        qx_6 = zeros(Ny1, Nx1)
        qy_6 = zeros(Ny1, Nx1)
        process_hydromechanical_solution!(
            S_6, vx_6, vy_6, pr_6, qx_6, qy_6, pf_6; coords=coords
        )

        # Dehydration creates significant fluid overpressure Pf > Pr in contact aureole
        delta_P = pf_6[4, 4] - pr_6[4, 4]
        @test delta_P > 5.0e4 # > 50 kPa overpressure
        @test maximum(pf_6) > maximum(pr_6)
        @test isfinite(delta_P)
    end
end

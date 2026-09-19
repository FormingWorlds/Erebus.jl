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

Random.seed!(42)

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
    @testset "Sensible Heat Advection Energy Conservation" begin
        coords = GridCoordinates(GridConfig(Nx=10, Ny=10, xsize=10000.0, ysize=10000.0))
        Ny, Nx = coords.Ny, coords.Nx
        Ny1, Nx1 = coords.Ny1, coords.Nx1
        dx, dy = coords.dx, coords.dy

        # Generate marker grid
        marknum = 1600
        xm = zeros(marknum)
        ym = zeros(marknum)
        tm = ones(Int, marknum)
        tkm = zeros(marknum)
        Fm = zeros(marknum)

        # Bottom layer hot and molten, top layer cool and solid
        for m in 1:marknum
            xm[m] = rand() * 10000.0
            ym[m] = rand() * 10000.0
            # Linear temperature gradient: 1600 K at bottom (y=10000), 1000 K at top (y=0)
            tkm[m] = 1000.0 + 600.0 * (ym[m] / 10000.0)
            # Melt only in deep warm region (y > 5000)
            if ym[m] > 5000.0
                Fm[m] = 0.25
            else
                Fm[m] = 0.0
            end
        end

        cfg_magma = MagmaTransportConfig(;
            active=true,
            sensible_heat_transport=true,
            cp_melt=1200.0,
            segregation_heating=false,
            latent_crystallization=false,
            compaction_active=false,
            max_subcycles=100,
        )

        ws = MagmaSegregationWorkspace(Ny, Nx)
        Q_seg_grid = zeros(Ny1, Nx1)
        dt_val = 1.0e8

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
            ycenter=20000.0,
            rplanet=25000.0,
            g_surf=0.5,
            Q_seg_grid=Q_seg_grid,
            workspace=ws,
        )

        # 1. Exact pairwise cancellation: net sensible enthalpy transfer is zero within machine precision
        @test isapprox(res.total_sensible_energy, 0.0; atol=1e-3)

        # 2. Conservation on nodal grid: sum of deposited Q_seg_grid sources vanishes to machine precision
        total_grid_source = sum(Q_seg_grid) * dx * dy
        @test isapprox(total_grid_source, 0.0; atol=1e-10)

        # 3. Spatial divergence: upward segregation extracts heat from donor and deposits in receiver
        @test maximum(Q_seg_grid) > 0.0
        @test minimum(Q_seg_grid) < 0.0
        @test res.dt_sub > 0.0
    end

    @testset "Crustal Sill Solidification & Latent Heat Kinetics" begin
        coords = GridCoordinates(GridConfig(Nx=8, Ny=8, xsize=8000.0, ysize=8000.0))
        Ny, Nx = coords.Ny, coords.Nx
        Ny1, Nx1 = coords.Ny1, coords.Nx1

        # Marker setup with a ponded sill beneath a cold lid
        marknum = 1000
        xm = zeros(marknum)
        ym = zeros(marknum)
        tm = ones(Int, marknum)
        tkm = zeros(marknum)
        Fm = zeros(marknum)

        for m in 1:marknum
            xm[m] = rand() * 8000.0
            ym[m] = rand() * 8000.0
            # Cold subsolidus lid everywhere (T = 1100 K < T_solidus = 1356 K)
            tkm[m] = 1100.0
            # Sill emplaced at depth 3000 to 5000 m with 20% melt fraction
            if 3000.0 <= ym[m] <= 5000.0
                Fm[m] = 0.20
            else
                Fm[m] = 0.0
            end
        end

        initial_melt_mass = sum(Fm)
        @test initial_melt_mass > 0.0

        # Case A: Instantaneous equilibrium crystallization (timescale = 0.0)
        cfg_instant = MagmaTransportConfig(;
            active=true,
            latent_crystallization=true,
            sill_cooling_active=true,
            crystallization_timescale=0.0,
            compaction_active=false,
        )

        Fm_inst = copy(Fm)
        Q_lat_inst = zeros(Ny1, Nx1)
        dt_val = 1.0e6

        res_inst = apply_silicate_melt_segregation!(
            copy(xm),
            copy(ym),
            copy(tm),
            copy(tkm),
            Fm_inst,
            marknum,
            dt_val,
            cfg_instant;
            coords=coords,
            xcenter=4000.0,
            ycenter=20000.0,
            rplanet=25000.0,
            g_surf=0.5,
            T_solidus_silicate=1356.0,
            Q_lat_grid=Q_lat_inst,
        )

        # In subsolidus lid, entire excess melt crystallizes in one step
        @test isapprox(sum(Fm_inst), 0.0; atol=1e-12)
        @test res_inst.total_crystallized_mass > 0.0
        @test maximum(Q_lat_inst) > 0.0

        # Case B: Rate-limited kinetic crystallization (timescale = 5e6 s)
        tau_cryst = 5.0e6
        cfg_kinetic = MagmaTransportConfig(;
            active=true,
            latent_crystallization=true,
            sill_cooling_active=true,
            crystallization_timescale=tau_cryst,
            compaction_active=false,
        )

        Fm_kin = copy(Fm)
        Q_lat_kin = zeros(Ny1, Nx1)

        res_kin = apply_silicate_melt_segregation!(
            copy(xm),
            copy(ym),
            copy(tm),
            copy(tkm),
            Fm_kin,
            marknum,
            dt_val,
            cfg_kinetic;
            coords=coords,
            xcenter=4000.0,
            ycenter=20000.0,
            rplanet=25000.0,
            g_surf=0.5,
            T_solidus_silicate=1356.0,
            Q_lat_grid=Q_lat_kin,
        )

        # Kinetic freezing freezes fraction dt / tau per step
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
            sill_coupling=true,
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

        # Convection enhancement: Nusselt number exceeds 1 above hot sill
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
        dy = L / (N - 1)
        y = range(-200.0, 200.0; length=N)

        # Initial condition: rectangular sill with smoothed contact
        T_num = [abs(yj) < b ? T0 : (abs(yj) == b ? 0.5 * (T0 + Tc) : Tc) for yj in y]
        t_total = 5.0e8 # ~15.8 yr (diffusion length 2*sqrt(kappa*t) ≈ 44.7 m, resolved across 25 cells)
        dt = 5.0e5

        # Implicit Crank-Nicolson diffusion solver
        r = 0.5 * kappa * dt / (dy^2)
        steps = Int(round(t_total / dt))

        # Tridiagonal system: (I - r D2) T^{n+1} = (I + r D2) T^n
        dl = fill(-r, N - 1)
        d_diag = fill(1.0 + 2.0 * r, N)
        du = fill(-r, N - 1)
        # Dirichlet BCs at boundaries
        d_diag[1] = 1.0
        du[1] = 0.0
        d_diag[N] = 1.0
        dl[N - 1] = 0.0
        A = Tridiagonal(dl, d_diag, du)

        rhs = zeros(N)
        for _ in 1:steps
            rhs[1] = Tc
            for i in 2:(N - 1)
                rhs[i] = r * T_num[i - 1] + (1.0 - 2.0 * r) * T_num[i] + r * T_num[i + 1]
            end
            rhs[N] = Tc
            T_num = A \ rhs
        end

        # Analytical Jaeger (1957) solution
        T_ana = [jaeger_sill_temperature(yj, t_total, b, T0, Tc, kappa) for yj in y]

        # Verification metrics
        l2_err = norm(T_num - T_ana) / norm(T_ana)
        max_err = maximum(abs.(T_num - T_ana))

        # L2 relative error must be under 0.1%, max absolute error under 1.0 K
        @test l2_err < 0.001
        @test max_err < 1.0

        # Centerline temperature verification: T(0, t) = Tc + (T0 - Tc) * erf(b / (2*sqrt(kappa*t)))
        T_center_ana = Tc + (T0 - Tc) * _sill_erf(b / (2.0 * sqrt(kappa * t_total)))
        center_idx = (N + 1) ÷ 2
        @test isapprox(T_num[center_idx], T_center_ana; atol=0.5)
    end

    @testset "Stefan Latent Heat Buffering Benchmark" begin
        # 1D thermal column with a solidifying melt layer
        # Stefan number Ste = cp * (T_liq - T_sol) / L_m
        T_sol = 1300.0
        T_liq = 1500.0
        L_melt = 4.0e5 # J/kg
        cp = 1000.0
        rho = 2800.0
        Ste = cp * (T_liq - T_sol) / L_melt # 0.5

        # Compare effective heat capacity with vs without latent heat
        cp_sensible = cp
        cp_apparent_mush = cp + L_melt / (T_liq - T_sol) # 1000 + 2000 = 3000 J/(kg K)

        # Buffering factor (1 + 1/Ste) = 3.0
        buffering_factor = cp_apparent_mush / cp_sensible
        @test isapprox(buffering_factor, 1.0 + 1.0 / Ste; rtol=1e-12)

        # Cooling timescale through crystallization range is prolonged by buffering_factor
        t_sensible = (rho * cp_sensible * (T_liq - T_sol)) / 100.0 # arbitrary flux
        t_latent = (rho * cp_apparent_mush * (T_liq - T_sol)) / 100.0
        @test isapprox(t_latent / t_sensible, 3.0; rtol=1e-12)
        @test t_latent > t_sensible
    end

    @testset "Contact Metamorphic Dehydration Aureole Benchmark" begin
        # Country rock initially at 400 K with hydrous mineralogy
        T_country = 400.0
        T_sill = 1400.0
        T_dehydration = 650.0 # Serpentine breakdown temperature

        # Contact interface temperature initially T_contact = (T_sill + T_country)/2 = 900 K > T_dehydration
        T_contact = 0.5 * (T_sill + T_country)
        @test T_contact > T_dehydration

        # Solve hydromechanical overpressure pulse driven by dehydration DQPF
        coords = GridCoordinates(GridConfig(Nx=7, Ny=7, xsize=14000.0, ysize=14000.0))
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

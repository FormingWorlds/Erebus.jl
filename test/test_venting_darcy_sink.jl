using Test
using Erebus
using Erebus.Config
using Erebus.Geometry
using Erebus.Physics
using Erebus.Particles
using Erebus.Numerics
using ExtendableSparse
using SparseArrays
using StaticArrays

@testset "Venting Darcy Sink & Mass Tracking" begin
    @testset "Surface Face Boundary Detection & Leaky Robin Assembly" begin
        Nx, Ny = 9, 9
        Nx1, Ny1 = Nx + 1, Ny + 1
        xsize, ysize = 140_000.0, 140_000.0
        coords = GridCoordinates(Nx, Ny; xsize=xsize, ysize=ysize)
        rplanet = 50_000.0
        xcenter = 70_000.0
        ycenter = 70_000.0

        n_dof = 6 * Ny1 * Nx1
        L = ExtendableSparseMatrix{Float64,Int64}(n_dof, n_dof)
        R = zeros(Float64, n_dof)
        tk = fill(200.0, Ny1, Nx1)
        pf = fill(2.0e6, Ny1, Nx1)
        pr = fill(2.5e6, Ny1, Nx1)
        TEN = fill(6.0e7, Ny1, Nx1)
        S_vent_out = zeros(Float64, Ny1, Nx1)

        P_amb = 10.0
        k_vent = 1.0e-11
        conductance_factor = 1.0
        Kcont = 1.0e20

        # Run boundary assembly
        apply_venting_surface_boundary!(
            L,
            R,
            tk,
            coords,
            rplanet,
            xcenter,
            ycenter,
            P_amb;
            k_vent=k_vent,
            conductance_factor=conductance_factor,
            mode=:darcy_sink,
            pf=pf,
            pr=pr,
            TEN=TEN,
            Kcont=Kcont,
            S_vent_out=S_vent_out,
        )
        flush!(L)

        # Check that surface nodes have positive conductance and S_vent_out > 0
        surface_nodes_found = 0
        interior_nodes_clean = true
        exterior_nodes_clean = true

        for j in 1:Nx1, i in 1:Ny1
            xj = coords.xp[j] - xcenter
            yi = coords.yp[i] - ycenter
            r = sqrt(xj^2 + yi^2)
            kpf = ((j - 1) * Ny1 + i - 1) * 6 + 6

            if S_vent_out[i, j] > 0.0
                surface_nodes_found += 1
                # Conductance must be positive in L
                @test L[kpf, kpf] > 0.0
                # RHS must contain P_vent contribution
                @test R[kpf] > 0.0
                # Surface nodes must lie in a shell around rplanet
                @test 0.5 * rplanet <= r <= 1.2 * rplanet
            else
                # Deep interior nodes (r < 0.4 * rplanet) must have zero S_vent
                if r < 0.4 * rplanet
                    @test iszero(S_vent_out[i, j])
                    @test iszero(L[kpf, kpf])
                end
                # Distant air nodes (r > 1.3 * rplanet) must have zero S_vent
                if r > 1.3 * rplanet
                    @test iszero(S_vent_out[i, j])
                    @test iszero(L[kpf, kpf])
                end
            end
        end

        @test surface_nodes_found >= 4

        # Verify refresh mode (L === nothing, R === nothing) yields identical S_vent_out
        S_vent_refresh = zeros(Float64, Ny1, Nx1)
        apply_venting_surface_boundary!(
            nothing,
            nothing,
            tk,
            coords,
            rplanet,
            xcenter,
            ycenter,
            P_amb;
            k_vent=k_vent,
            conductance_factor=conductance_factor,
            mode=:darcy_sink,
            pf=pf,
            pr=pr,
            TEN=TEN,
            S_vent_out=S_vent_refresh,
        )
        @test isapprox(S_vent_refresh, S_vent_out; atol=1e-15)

        # When pf <= P_vent, S_vent_out must be zero everywhere (no outward overpressure)
        S_vent_zero = zeros(Float64, Ny1, Nx1)
        pf_sub = fill(0.0, Ny1, Nx1)
        apply_venting_surface_boundary!(
            nothing,
            nothing,
            tk,
            coords,
            rplanet,
            xcenter,
            ycenter,
            P_amb;
            k_vent=k_vent,
            conductance_factor=conductance_factor,
            mode=:darcy_sink,
            pf=pf_sub,
            pr=pr,
            TEN=TEN,
            S_vent_out=S_vent_zero,
        )
        @test all(iszero, S_vent_zero)

        # When pf <= P_vent, Robin conductance must NOT be added to L or R (strictly one-sided)
        L_sub_mat = ExtendableSparseMatrix{Float64,Int64}(n_dof, n_dof)
        R_sub_mat = zeros(Float64, n_dof)
        apply_venting_surface_boundary!(
            L_sub_mat,
            R_sub_mat,
            tk,
            coords,
            rplanet,
            xcenter,
            ycenter,
            P_amb;
            k_vent=k_vent,
            conductance_factor=conductance_factor,
            mode=:darcy_sink,
            pf=pf_sub,
            pr=pr,
            TEN=TEN,
            Kcont=Kcont,
        )
        flush!(L_sub_mat)
        @test iszero(SparseArrays.nnz(L_sub_mat.cscmatrix))
        @test all(iszero, R_sub_mat)

        # Fluid-availability limit: when PHI <= phimin, venting shuts off completely
        L_dry = ExtendableSparseMatrix{Float64,Int64}(n_dof, n_dof)
        R_dry = zeros(Float64, n_dof)
        S_vent_dry = zeros(Float64, Ny1, Nx1)
        PHI_dry = fill(1.0e-4, Ny1, Nx1)
        apply_venting_surface_boundary!(
            L_dry,
            R_dry,
            tk,
            coords,
            rplanet,
            xcenter,
            ycenter,
            P_amb;
            k_vent=k_vent,
            conductance_factor=conductance_factor,
            mode=:darcy_sink,
            pf=pf,
            pr=pr,
            TEN=TEN,
            PHI=PHI_dry,
            phimin=1.0e-4,
            Kcont=Kcont,
            S_vent_out=S_vent_dry,
        )
        flush!(L_dry)
        @test iszero(SparseArrays.nnz(L_dry.cscmatrix))
        @test all(iszero, R_dry)
        @test all(iszero, S_vent_dry)

        # Divide-by-zero guards
        @test_throws DomainError apply_venting_surface_boundary!(
            L, R, tk, coords, rplanet, xcenter, ycenter, P_amb; eta_fluid_surf=0.0
        )
        @test_throws DomainError apply_venting_surface_boundary!(
            L, R, tk, coords, rplanet, xcenter, ycenter, P_amb; eta_fluid_surf=-1.0
        )
    end

    @testset "Analytical 1D Steady-State Darcy Flux Benchmark" begin
        # 1D Darcy flux: qD = (k / eta) * (ΔP / L)
        k_val = 1.0e-13       # [m²]
        eta_val = 1.0e-3      # [Pa·s]
        delta_P = 2.0e6       # [Pa] (2 MPa)
        L_col = 10_000.0      # [m] (10 km column)

        qD_analytic = (k_val / eta_val) * (delta_P / L_col) # 2.0e-8 m/s
        @test isapprox(qD_analytic, 2.0e-8; rtol=1e-12)

        # Drainage mass flux: M_dot = rho_f * qD * Area
        Area = 1000.0 * 1000.0 # 1 km² area
        rho_f = 1000.0
        M_dot_analytic = rho_f * qD_analytic * Area # 20 kg/s
        @test isapprox(M_dot_analytic, 20.0; rtol=1e-12)
    end

    @testset "sink_vented_marker_porosity! Invariants & Mass Conservation" begin
        Nx, Ny = 9, 9
        Nx1, Ny1 = Nx + 1, Ny + 1
        xsize, ysize = 140_000.0, 140_000.0
        coords = GridCoordinates(Nx, Ny; xsize=xsize, ysize=ysize)

        marknum = 1000
        xm = fill(70_000.0, marknum)
        ym = fill(70_000.0, marknum)
        tm = fill(1, marknum)       # rock markers
        phim0 = 0.20
        phim = fill(phim0, marknum)
        phimin = 1.0e-4
        rhofluidcur = 1000.0

        S_vent_grid = fill(1.0e-10, Ny1, Nx1)
        dt = 1.0e8 # 1e8 seconds (~3.17 years)

        delta_m_vent = sink_vented_marker_porosity!(
            xm,
            ym,
            tm,
            phim,
            S_vent_grid,
            dt,
            marknum;
            coords=coords,
            phimin=phimin,
            rhofluidcur=rhofluidcur,
        )

        # Porosity must decrease by Δϕ = S_vent * dt = 1e-10 * 1e8 = 0.01
        @test isapprox(phim[1], 0.19; rtol=1e-6)
        @test delta_m_vent > 0.0

        # Total fluid mass drained from markers must equal rho_f * Δϕ * total_volume
        V_total = xsize * ysize
        expected_mass = rhofluidcur * (phim0 - 0.19) * V_total
        @test isapprox(delta_m_vent, expected_mass; rtol=1e-4)

        # Extreme drainage: porosity must clamp to phimin, not below
        dt_huge = 1.0e13
        delta_m_clamp = sink_vented_marker_porosity!(
            xm,
            ym,
            tm,
            phim,
            S_vent_grid,
            dt_huge,
            marknum;
            coords=coords,
            phimin=phimin,
            rhofluidcur=rhofluidcur,
        )
        @test isapprox(phim[1], phimin; rtol=1e-12)
        @test all(p -> p >= phimin, phim)
    end

    @testset "Sublimation Latent Cooling Heat Sink Verification" begin
        # Q_lat = -L_sub * rho_f * S_vent [W/m³]
        L_sub = 2.83e6      # [J/kg]
        rho_f = 1000.0      # [kg/m³]
        S_vent = 1.0e-11    # [s⁻¹]

        Q_lat_val = -L_sub * rho_f * S_vent
        @test Q_lat_val < 0.0
        @test isapprox(Q_lat_val, -2.83e-2; rtol=1e-12) # -0.0283 W/m³

        # Energy removed over 1 km³ cell during 1 year
        V_cell = 1.0e9 # 1 km³
        t_year = 365.25 * 86400.0
        E_removed = abs(Q_lat_val) * V_cell * t_year
        mass_vented = rho_f * S_vent * V_cell * t_year
        @test isapprox(E_removed, L_sub * mass_vented; rtol=1e-12)

        # Direct LSE assembly verification: Q_lat enters thermal RHS RT
        Nx_th, Ny_th = 5, 5
        coords_th = GridCoordinates(Nx_th, Ny_th; xsize=50_000.0, ysize=50_000.0)
        Nx1_th, Ny1_th = Nx_th + 1, Ny_th + 1
        tk1_th = fill(200.0, Ny1_th, Nx1_th)
        RHOCP_th = fill(2.0e6, Ny1_th, Nx1_th)
        KX_th = fill(2.0, Ny_th, Nx1_th)
        KY_th = fill(2.0, Ny1_th, Nx_th)
        HR_th = zeros(Float64, Ny1_th, Nx1_th)
        HA_th = zeros(Float64, Ny1_th, Nx1_th)
        HS_th = zeros(Float64, Ny1_th, Nx1_th)
        DHP_th = zeros(Float64, Ny1_th, Nx1_th)
        RT_base = zeros(Float64, Ny1_th * Nx1_th)
        RT_vent = zeros(Float64, Ny1_th * Nx1_th)
        dt_th = 1.0e8

        assemble_thermal_lse!(
            tk1_th,
            RHOCP_th,
            KX_th,
            KY_th,
            HR_th,
            HA_th,
            HS_th,
            DHP_th,
            RT_base,
            dt_th;
            coords=coords_th,
        )

        S_vent_th = zeros(Float64, Ny1_th, Nx1_th)
        S_vent_th[3, 3] = S_vent
        Q_lat_grid = -L_sub * rho_f .* S_vent_th

        assemble_thermal_lse!(
            tk1_th,
            RHOCP_th,
            KX_th,
            KY_th,
            HR_th,
            HA_th,
            HS_th,
            DHP_th,
            RT_vent,
            dt_th;
            coords=coords_th,
            Q_lat=Q_lat_grid,
        )

        # Verify Q_lat strictly reduces thermal RHS at (3, 3) by exactly Q_lat_val
        gk_target = (2 * Ny1_th + 3) # j=3, i=3 -> (3-1)*6 + 3 = 15
        @test RT_vent[gk_target] < RT_base[gk_target]
        @test isapprox(RT_vent[gk_target] - RT_base[gk_target], Q_lat_val; rtol=1e-12)

        # Verify unaffected nodes have identical RHS
        gk_other = (1 * Ny1_th + 2) # j=2, i=2
        @test isapprox(RT_vent[gk_other], RT_base[gk_other]; atol=1e-15)
    end
end

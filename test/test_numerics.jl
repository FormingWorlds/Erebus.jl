@testset "Numerics" begin
    @testset "assemble_gravitational_lse!(): operator structure and discrete Poisson invariants" begin
        RP = zeros(Float64, Nx1 * Ny1)
        RHO = zeros(Float64, Ny1, Nx1)
        LP = Erebus.assemble_gravitational_lse!(RHO, RP)

        # 1. Zero-density invariant: when rho = 0, RHS must be zero everywhere
        @test all(iszero, RP)

        # 2. Boundary condition rows: on boundaries, LP[gk, gk] = 1.0 and off-diagonals are 0.0
        xc_val = xcenter
        yc_val = ycenter
        r_limit = min(xc_val, yc_val)
        boundary_nodes = 0
        interior_nodes = 0
        for j in 1:Nx1, i in 1:Ny1
            gk = (j - 1) * Ny1 + i
            if Erebus.is_gravitational_boundary(
                i, j, Ny1, Nx1, xp, yp, xc_val, yc_val, r_limit
            )
                boundary_nodes += 1
                @test isapprox(LP[gk, gk], 1.0; rtol=1e-12)
                @test isapprox(RP[gk], 0.0; atol=1e-12)
            else
                interior_nodes += 1
                # 3. Interior diagonal is negative: -2*(1/dx^2 + 1/dy^2)
                expected_diag = -2.0 * (inv(dx^2) + inv(dy^2))
                @test isapprox(LP[gk, gk], expected_diag; rtol=1e-12)
                # 4. Sum of coefficients on interior rows is zero (harmonic null space for constant potential)
                row_sum =
                    LP[gk, gk - Ny1] +
                    LP[gk, gk - 1] +
                    LP[gk, gk] +
                    LP[gk, gk + 1] +
                    LP[gk, gk + Ny1]
                @test isapprox(row_sum, 0.0; atol=1e-12)
            end
        end
        @test boundary_nodes > 0
        @test interior_nodes > 0

        # 5. Symmetry of interior Laplacian couplings
        for j in 2:(Nx1 - 1), i in 2:(Ny1 - 1)
            gk = (j - 1) * Ny1 + i
            if !Erebus.is_gravitational_boundary(
                i, j, Ny1, Nx1, xp, yp, xc_val, yc_val, r_limit
            )
                # Check vertical coupling if neighbor is also interior
                if !Erebus.is_gravitational_boundary(
                    i + 1, j, Ny1, Nx1, xp, yp, xc_val, yc_val, r_limit
                )
                    gk_next = gk + 1
                    @test isapprox(LP[gk, gk_next], LP[gk_next, gk]; rtol=1e-12)
                end
                # Check horizontal coupling if neighbor is also interior
                if !Erebus.is_gravitational_boundary(
                    i, j + 1, Ny1, Nx1, xp, yp, xc_val, yc_val, r_limit
                )
                    gk_next = gk + Ny1
                    @test isapprox(LP[gk, gk_next], LP[gk_next, gk]; rtol=1e-12)
                end
            end
        end

        # 6. RHS linearity with density
        RHO1 = fill(3000.0, Ny1, Nx1)
        RP1 = zeros(Float64, Nx1 * Ny1)
        Erebus.assemble_gravitational_rhs!(RHO1, RP1)
        RHO2 = fill(6000.0, Ny1, Nx1)
        RP2 = zeros(Float64, Nx1 * Ny1)
        Erebus.assemble_gravitational_rhs!(RHO2, RP2)
        @test isapprox(RP2, 2.0 .* RP1; rtol=1e-12)

        # 7. Physical prefactor verification against independently derived constant (8πG/3)
        RP_unit = zeros(Float64, Nx1 * Ny1)
        Erebus.assemble_gravitational_rhs!(fill(1.0, Ny1, Nx1), RP_unit)
        gk_mid = (div(Nx1, 2) - 1) * Ny1 + div(Ny1, 2)
        # Independently derived 8πG/3 with G = 6.672e-11 m^3/(kg s^2)
        expected_prefactor = 5.58952164926696e-10
        @test isapprox(RP_unit[gk_mid], expected_prefactor; rtol=1e-10)
        @test abs(RP_unit[gk_mid] - (4.0 * π * G)) > 1e-10 # prefactor discrimination (not 4πG)
        @test abs(RP_unit[gk_mid] - (4.0 / 3.0 * π * G)) > 1e-10 # prefactor discrimination (not 4πG/3)
        @test RP_unit[gk_mid] > 0.0
    end

    @testset "process_gravitational_solution!(): gradient operator invariants" begin
        FI = zeros(Float64, Ny1, Nx1)
        gx = zeros(Float64, Ny1, Nx1)
        gy = zeros(Float64, Ny1, Nx1)
        SP = zeros(Float64, Nx1 * Ny1)

        # 1. Linear potential along x: Φ(x) = -g0 * x => gx = g0, gy = 0
        g0 = 0.05
        for j in 1:Nx1, i in 1:Ny1
            FI[i, j] = -g0 * xp[j]
        end
        SP .= vec(FI)
        Erebus.process_gravitational_solution!(SP, FI, gx, gy)
        @test all(isapprox.(gx[:, 1:Nx], g0; rtol=1e-12))
        @test all(isapprox.(gy[1:Ny, :], 0.0; atol=1e-12))

        # 2. Linear potential along y: Φ(y) = -g0 * y => gx = 0, gy = g0
        for j in 1:Nx1, i in 1:Ny1
            FI[i, j] = -g0 * yp[i]
        end
        SP .= vec(FI)
        Erebus.process_gravitational_solution!(SP, FI, gx, gy)
        @test all(isapprox.(gx[:, 1:Nx], 0.0; atol=1e-12))
        @test all(isapprox.(gy[1:Ny, :], g0; rtol=1e-12))

        # 3. Gauge invariance: uniform shift in potential does not alter gravitational acceleration
        SP_shifted = SP .+ 1.2345e5
        gx_shift = zeros(Float64, Ny1, Nx1)
        gy_shift = zeros(Float64, Ny1, Nx1)
        Erebus.process_gravitational_solution!(SP_shifted, FI, gx_shift, gy_shift)
        @test isapprox(gx_shift, gx; rtol=1e-12)
        @test isapprox(gy_shift, gy; rtol=1e-12)
    end

    @testset "compute_gravity_solution!(): planetary radial benchmark and symmetry" begin
        SP = zeros(Float64, Nx1 * Ny1)
        RP = zeros(Float64, Nx1 * Ny1)
        FI = zeros(Float64, Ny1, Nx1)
        gx = zeros(Float64, Ny1, Nx1)
        gy = zeros(Float64, Ny1, Nx1)

        # Construct a symmetric uniform planetesimal core centered in domain
        RHO = zeros(Float64, Ny1, Nx1)
        r_planet = 40000.0
        rho_planet = 3300.0
        r_limit = min(xcenter, ycenter)
        for j in 1:Nx1, i in 1:Ny1
            r_dist = sqrt((xp[j] - xcenter)^2 + (yp[i] - ycenter)^2)
            if r_dist <= r_planet
                RHO[i, j] = rho_planet
            end
        end

        Erebus.compute_gravity_solution!(SP, RP, RHO, FI, gx, gy)

        # 1. Attractor direction invariant:
        # Gravity must point inward towards center (xcenter, ycenter) inside the active domain
        for j in 1:Nx, i in 1:Ny1
            x_mid = 0.5 * (xp[j] + xp[j + 1])
            r_node = hypot(x_mid - xcenter, yp[i] - ycenter)
            if r_node < (r_limit - 2.0 * dx)
                if x_mid > xcenter + dx
                    @test gx[i, j] < 0.0 # points left toward center
                elseif x_mid < xcenter - dx
                    @test gx[i, j] > 0.0 # points right toward center
                end
            end
        end
        for j in 1:Nx1, i in 1:Ny
            y_mid = 0.5 * (yp[i] + yp[i + 1])
            r_node = hypot(xp[j] - xcenter, y_mid - ycenter)
            if r_node < (r_limit - 2.0 * dy)
                if y_mid > ycenter + dy
                    @test gy[i, j] < 0.0 # points up toward center
                elseif y_mid < ycenter - dy
                    @test gy[i, j] > 0.0 # points down toward center
                end
            end
        end

        # 2. Quadrant reflection anti-symmetry: gx(x_c + d) = -gx(x_c - d)
        j_center = div(Nx1 + 1, 2)
        i_center = div(Ny1 + 1, 2)
        for dj in 1:(div(Nx1, 2) - 2)
            j_right = j_center + dj
            j_left = j_center - dj
            if j_right <= Nx && j_left >= 1
                @test isapprox(gx[i_center, j_right], -gx[i_center, j_left]; rtol=1e-6)
            end
        end

        # 3. Scale and magnitude discrimination guard
        g_max = maximum(hypot.(gx[1:Ny, 1:Nx], gy[1:Ny, 1:Nx]))
        @test 0.005 < g_max < 0.1 # order of magnitude for 40 km body at 3300 kg/m^3
        @test g_max > 0.0 # sign
    end

    @testset "recompute_bulk_viscosity!(): harmonic averaging and compaction scaling" begin
        ETAP = zeros(Ny1, Nx1)
        ETAPHI = zeros(Ny1, Nx1)
        PHI = fill(0.15, Ny1, Nx1)

        # 1. Constant field invariance: harmonic mean of uniform ETA is identical to ETA
        eta0 = 1.0e18
        ETA_const = fill(eta0, Ny, Nx)
        Erebus.recompute_bulk_viscosity!(ETA_const, ETAP, ETAPHI, PHI, etaphikoef)
        for j in 2:Nx, i in 2:Ny
            @test isapprox(ETAP[i, j], eta0; rtol=1e-12)
            expected_etaphi = etaphikoef * eta0 / PHI[i, j]
            @test isapprox(ETAPHI[i, j], expected_etaphi; rtol=1e-12)
        end

        # 2. Extremum bounding: min(ETA_4) <= ETAP[i, j] <= max(ETA_4)
        ETA_var = 1.0e17 .+ rand(rgen, Ny, Nx) * 1.0e18
        Erebus.recompute_bulk_viscosity!(ETA_var, ETAP, ETAPHI, PHI, etaphikoef)
        for j in 2:Nx, i in 2:Ny
            eta_min = min(
                ETA_var[i - 1, j - 1], ETA_var[i, j - 1], ETA_var[i - 1, j], ETA_var[i, j]
            )
            eta_max = max(
                ETA_var[i - 1, j - 1], ETA_var[i, j - 1], ETA_var[i - 1, j], ETA_var[i, j]
            )
            @test ETAP[i, j] >= eta_min - 1e-6
            @test ETAP[i, j] <= eta_max + 1e-6
        end

        # 3. Compaction resistance porosity inverse scaling: eta_phi ~ 1/phi
        PHI_low = fill(0.05, Ny1, Nx1)
        PHI_high = fill(0.20, Ny1, Nx1)
        ETAPHI_low = zeros(Ny1, Nx1)
        ETAPHI_high = zeros(Ny1, Nx1)
        Erebus.recompute_bulk_viscosity!(ETA_const, ETAP, ETAPHI_low, PHI_low, etaphikoef)
        Erebus.recompute_bulk_viscosity!(ETA_const, ETAP, ETAPHI_high, PHI_high, etaphikoef)
        for j in 2:Nx, i in 2:Ny
            @test ETAPHI_low[i, j] > ETAPHI_high[i, j] # strict monotonicity
            @test isapprox(ETAPHI_low[i, j] / ETAPHI_high[i, j], 0.20 / 0.05; rtol=1e-12)
        end
    end

    @testset "get_viscosities_stresses_density_gradients!(): Maxwell viscoelastic asymptotes" begin
        dt_test = 1.0e10
        ETAcomp = zeros(Ny, Nx)
        ETAPcomp = zeros(Ny1, Nx1)
        SXYcomp = zeros(Ny, Nx)
        SXXcomp = zeros(Ny1, Nx1)
        SYYcomp = zeros(Ny1, Nx1)
        dRHOXdx = zeros(Ny1, Nx1)
        dRHOXdy = zeros(Ny1, Nx1)
        dRHOYdx = zeros(Ny1, Nx1)
        dRHOYdy = zeros(Ny1, Nx1)

        # 1. Viscous limit: G*dt >> ETA => ETAcomp -> ETA, SXXcomp -> 0
        ETA_fluid = fill(1.0e12, Ny, Nx)
        ETAP_fluid = fill(1.0e12, Ny1, Nx1)
        GGG_high = fill(1.0e11, Ny, Nx) # G*dt = 1e21 >> 1e12
        GGGP_high = fill(1.0e11, Ny1, Nx1)
        SXY0 = fill(1.0e6, Ny, Nx)
        SXX0 = fill(1.0e6, Ny1, Nx1)
        RHOX = fill(3000.0, Ny1, Nx1)
        RHOY = fill(3000.0, Ny1, Nx1)

        Erebus.get_viscosities_stresses_density_gradients!(
            ETA_fluid,
            ETAP_fluid,
            GGG_high,
            GGGP_high,
            SXY0,
            SXX0,
            RHOX,
            RHOY,
            dt_test,
            ETAcomp,
            ETAPcomp,
            SXYcomp,
            SXXcomp,
            SYYcomp,
            dRHOXdx,
            dRHOXdy,
            dRHOYdx,
            dRHOYdy,
        )
        for j in 2:Nx, i in 2:Ny
            @test isapprox(ETAcomp[i, j], ETA_fluid[i, j]; rtol=1e-6)
            @test isapprox(ETAPcomp[i, j], ETAP_fluid[i, j]; rtol=1e-6)
            @test isapprox(SXYcomp[i, j], 0.0; atol=1e-2)
            @test isapprox(SXXcomp[i, j], 0.0; atol=1e-2)
            @test isapprox(SYYcomp[i, j], 0.0; atol=1e-2)
        end

        # 2. Elastic limit: ETA >> G*dt => ETAcomp -> G*dt, SXXcomp -> SXX0
        ETA_solid = fill(1.0e24, Ny, Nx)
        ETAP_solid = fill(1.0e24, Ny1, Nx1)
        GGG_low = fill(1.0e9, Ny, Nx) # G*dt = 1e19 << 1e24
        GGGP_low = fill(1.0e9, Ny1, Nx1)
        Erebus.get_viscosities_stresses_density_gradients!(
            ETA_solid,
            ETAP_solid,
            GGG_low,
            GGGP_low,
            SXY0,
            SXX0,
            RHOX,
            RHOY,
            dt_test,
            ETAcomp,
            ETAPcomp,
            SXYcomp,
            SXXcomp,
            SYYcomp,
            dRHOXdx,
            dRHOXdy,
            dRHOYdx,
            dRHOYdy,
        )
        for j in 2:Nx, i in 2:Ny
            @test isapprox(ETAcomp[i, j], GGG_low[i, j] * dt_test; rtol=1e-4)
            @test isapprox(ETAPcomp[i, j], GGGP_low[i, j] * dt_test; rtol=1e-4)
            @test isapprox(SXYcomp[i, j], SXY0[i, j]; rtol=1e-4)
            @test isapprox(SXXcomp[i, j], SXX0[i, j]; rtol=1e-4)
            @test isapprox(SYYcomp[i, j], -SXX0[i, j]; rtol=1e-4)
        end

        # 3. Density gradient antisymmetry on centered symmetric profile
        RHOX_sym = zeros(Ny1, Nx1)
        RHOY_sym = zeros(Ny1, Nx1)
        for j in 1:Nx1, i in 1:Ny1
            RHOX_sym[i, j] = 3000.0 + 300.0 * cos(π * (xp[j] - xcenter) / xsize)
            RHOY_sym[i, j] = 3000.0 + 300.0 * cos(π * (yp[i] - ycenter) / ysize)
        end
        Erebus.get_viscosities_stresses_density_gradients!(
            ETA_fluid,
            ETAP_fluid,
            GGG_high,
            GGGP_high,
            SXY0,
            SXX0,
            RHOX_sym,
            RHOY_sym,
            dt_test,
            ETAcomp,
            ETAPcomp,
            SXYcomp,
            SXXcomp,
            SYYcomp,
            dRHOXdx,
            dRHOXdy,
            dRHOYdx,
            dRHOYdy,
        )
        j_c = div(Nx1, 2) # 17 for Nx1 = 34
        i_c = div(Ny1, 2) # 17 for Ny1 = 34
        for k in 0:(j_c - 2)
            @test dRHOXdx[i_c, j_c - k] > 0.0
            @test dRHOXdx[i_c, j_c + 1 + k] < 0.0
            @test isapprox(dRHOXdx[i_c, j_c - k], -dRHOXdx[i_c, j_c + 1 + k]; rtol=1e-12)
        end
        for k in 0:(i_c - 2)
            @test dRHOYdy[i_c - k, j_c] > 0.0
            @test dRHOYdy[i_c + 1 + k, j_c] < 0.0
            @test isapprox(dRHOYdy[i_c - k, j_c], -dRHOYdy[i_c + 1 + k, j_c]; rtol=1e-12)
        end
    end

    @testset "setup_*_lse() constructors: dimensions and dynamic coordinates" begin
        # 1. Default grid sizes
        R_h, S_h = Erebus.setup_hydromechanical_lse()
        @test size(R_h) == (Nx1 * Ny1 * 6,)
        @test size(S_h) == (Nx1 * Ny1 * 6,)
        @test eltype(R_h) === Float64
        @test eltype(S_h) === Float64

        RP_t, SP_t = Erebus.setup_thermal_lse()
        @test size(RP_t) == (Nx1 * Ny1,)
        @test size(SP_t) == (Nx1 * Ny1,)

        RT_g, ST_g = Erebus.setup_gravitational_lse()
        @test size(RT_g) == (Nx1 * Ny1,)
        @test size(ST_g) == (Nx1 * Ny1,)

        # 2. Dynamic grid coordinates
        dyn_coords = Erebus.GridCoordinates(15, 15; xsize=100000.0, ysize=100000.0)
        R_dyn, S_dyn = Erebus.setup_hydromechanical_lse(dyn_coords)
        @test size(R_dyn) == (dyn_coords.Nx1 * dyn_coords.Ny1 * 6,)
        @test size(S_dyn) == (dyn_coords.Nx1 * dyn_coords.Ny1 * 6,)
    end

    @testset "assemble_hydromechanical_lse!(): Stokes null space, fluid conservation, compaction balance" begin
        dt = dt_longest
        ETA = fill(1.0e18, Ny, Nx)
        ETAP = fill(1.0e18, Ny1, Nx1)
        GGG = fill(1.0e10, Ny, Nx)
        GGGP = fill(1.0e10, Ny1, Nx1)
        SXY0 = zeros(Ny, Nx)
        SXX0 = zeros(Ny1, Nx1)
        RHOX = fill(3300.0, Ny1, Nx1)
        RHOY = fill(3300.0, Ny1, Nx1)
        RHOFX = fill(1000.0, Ny1, Nx1)
        RHOFY = fill(1000.0, Ny1, Nx1)
        RX = fill(1.0e10, Ny1, Nx1)
        RY = fill(1.0e10, Ny1, Nx1)
        ETAPHI = fill(1.0e16, Ny1, Nx1)
        BETTAPHI = fill(1.0e-10, Ny1, Nx1)
        PHI = fill(0.1, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = zeros(Ny1, Nx1)
        pr0 = zeros(Ny1, Nx1)
        pf0 = zeros(Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)
        R = zeros(Nx1 * Ny1 * 6)

        L = Erebus.assemble_hydromechanical_lse!(
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
            BETTAPHI,
            PHI,
            gx,
            gy,
            pr0,
            pf0,
            DMP,
            dt,
            R;
            betasolid=0.0,
            betafluid=0.0,
        )

        # 1. Rigid Body Translation Invariant: sum of solid velocity stencil on interior Vx rows is 0
        for j in 3:(Nx - 2), i in 3:(Ny - 2)
            kvx = ((j - 1) * Ny1 + i - 1) * 6 + 1
            vx_sum =
                L[kvx, kvx - 6 * Ny1] +
                L[kvx, kvx - 6] +
                L[kvx, kvx] +
                L[kvx, kvx + 6] +
                L[kvx, kvx + 6 * Ny1]
            @test isapprox(vx_sum, 0.0; atol=1e-12)
        end

        # 2. Rigid Body Translation Invariant: sum of solid velocity stencil on interior Vy rows is 0
        for j in 3:(Nx - 2), i in 3:(Ny - 2)
            kvy = ((j - 1) * Ny1 + i - 1) * 6 + 2
            vy_sum =
                L[kvy, kvy - 6 * Ny1] +
                L[kvy, kvy - 6] +
                L[kvy, kvy] +
                L[kvy, kvy + 6] +
                L[kvy, kvy + 6 * Ny1]
            @test isapprox(vy_sum, 0.0; atol=1e-12)
        end

        # 3. Discrete Fluid Mass Conservation: Darcy divergence coefficients sum to 0 in row kpf
        for j in 3:(Nx - 2), i in 3:(Ny - 2)
            kvx = ((j - 1) * Ny1 + i - 1) * 6 + 1
            kpf = kvx + 5
            kqx = kvx + 3
            kqy = kvx + 4
            qx_sum = L[kpf, kqx - Ny1 * 6] + L[kpf, kqx]
            @test isapprox(qx_sum, 0.0; atol=1e-12)
            qy_sum = L[kpf, kqy - 6] + L[kpf, kqy]
            @test isapprox(qy_sum, 0.0; atol=1e-12)
        end

        # 4. Two-Phase Volume Balance: solid compaction rate equals fluid expulsion rate
        for j in 3:(Nx - 2), i in 3:(Ny - 2)
            kvx = ((j - 1) * Ny1 + i - 1) * 6 + 1
            kpm = kvx + 2
            kpf = kvx + 5
            @test isapprox(L[kpm, kpm] + L[kpm, kpf], 0.0; atol=1e-12)
            @test isapprox(L[kpf, kpm] + L[kpf, kpf], 0.0; atol=1e-12)
            @test isapprox(L[kpm, kpm], -L[kpf, kpm]; rtol=1e-12)
        end

        # 5. Boundary Condition Rows: ghost and external rows have 1.0 on diagonal
        for j in 1:Nx1, i in 1:Ny1
            kvx = ((j - 1) * Ny1 + i - 1) * 6 + 1
            if i == 1 || i == Ny1 || j == 1 || j == Nx || j == Nx1
                @test isapprox(L[kvx, kvx], 1.0; rtol=1e-12)
            end
        end
    end

    @testset "process_hydromechanical_solution!(): round-trip injection preservation" begin
        vx_in = zeros(Ny1, Nx1)
        vy_in = zeros(Ny1, Nx1)
        pr_in = zeros(Ny1, Nx1)
        qxD_in = zeros(Ny1, Nx1)
        qyD_in = zeros(Ny1, Nx1)
        pf_in = zeros(Ny1, Nx1)
        for j in 1:Nx1, i in 1:Ny1
            vx_in[i, j] = 1.0e-9 * (i + 2 * j)
            vy_in[i, j] = 2.0e-9 * (3 * i - j)
            pr_in[i, j] = 1.0e6 * (i * j)
            qxD_in[i, j] = 5.0e-10 * (i^2 + j)
            qyD_in[i, j] = 7.0e-10 * (j^2 - i)
            pf_in[i, j] = 0.5e6 * (i + j^2)
        end

        # Pack into algebraic vector S with Kcont pressure scaling
        S = zeros(Float64, Nx1 * Ny1 * 6)
        for j in 1:Nx1, i in 1:Ny1
            kvx = ((j - 1) * Ny1 + i - 1) * 6 + 1
            S[kvx] = vx_in[i, j]
            S[kvx + 1] = vy_in[i, j]
            S[kvx + 2] = pr_in[i, j] / Kcont
            S[kvx + 3] = qxD_in[i, j]
            S[kvx + 4] = qyD_in[i, j]
            S[kvx + 5] = pf_in[i, j] / Kcont
        end

        vx_out = zeros(Ny1, Nx1)
        vy_out = zeros(Ny1, Nx1)
        pr_out = zeros(Ny1, Nx1)
        qxD_out = zeros(Ny1, Nx1)
        qyD_out = zeros(Ny1, Nx1)
        pf_out = zeros(Ny1, Nx1)

        Erebus.process_hydromechanical_solution!(
            S, vx_out, vy_out, pr_out, qxD_out, qyD_out, pf_out
        )

        # 1. Round-trip bit-for-bit extraction equality
        @test isapprox(vx_out, vx_in; rtol=1e-12)
        @test isapprox(vy_out, vy_in; rtol=1e-12)
        @test isapprox(pr_out, pr_in; rtol=1e-12)
        @test isapprox(qxD_out, qxD_in; rtol=1e-12)
        @test isapprox(qyD_out, qyD_in; rtol=1e-12)
        @test isapprox(pf_out, pf_in; rtol=1e-12)

        # 2. Component isolation: pure fluid pressure perturbation does not bleed into velocities
        S_pf_only = zeros(Float64, Nx1 * Ny1 * 6)
        for j in 1:Nx1, i in 1:Ny1
            kpf = ((j - 1) * Ny1 + i - 1) * 6 + 6
            S_pf_only[kpf] = 1.0e6 / Kcont
        end
        Erebus.process_hydromechanical_solution!(
            S_pf_only, vx_out, vy_out, pr_out, qxD_out, qyD_out, pf_out
        )
        @test all(iszero, vx_out)
        @test all(iszero, vy_out)
        @test all(iszero, pr_out)
        @test all(iszero, qxD_out)
        @test all(iszero, qyD_out)
        @test all(isapprox.(pf_out, 1.0e6; rtol=1e-12))
    end

    @testset "compute_Aϕ!(): compaction rate equilibrium and sign invariants" begin
        dt = dt_longest
        APHI = zeros(Ny1, Nx1)
        ETAPHI = fill(1.0e16, Ny1, Nx1)
        BETTAPHI = fill(1.0e-10, Ny1, Nx1)
        PHI = fill(0.15, Ny1, Nx1)

        # 1. Compaction Equilibrium: when pr == pf and steady state, Aphi = 0
        pr_eq = fill(5.0e6, Ny1, Nx1)
        pf_eq = fill(5.0e6, Ny1, Nx1)
        aphimax = Erebus.compute_Aϕ!(
            APHI, ETAPHI, BETTAPHI, PHI, pr_eq, pf_eq, pr_eq, pf_eq, dt
        )
        @test isapprox(aphimax, 0.0; atol=1e-12)
        @test all(isapprox.(APHI[2:Ny, 2:Nx], 0.0; atol=1e-12))

        # 2. Overpressure: pr > pf drives positive compaction rate Aphi > 0
        pr_high = fill(10.0e6, Ny1, Nx1)
        pf_low = fill(2.0e6, Ny1, Nx1)
        aphimax = Erebus.compute_Aϕ!(
            APHI, ETAPHI, BETTAPHI, PHI, pr_high, pf_low, pr_high, pf_low, dt
        )
        @test aphimax > 0.0
        @test all(APHI[2:Ny, 2:Nx] .> 0.0)

        # 3. Underpressure: pr < pf drives dilation rate Aphi < 0
        aphimax = Erebus.compute_Aϕ!(
            APHI, ETAPHI, BETTAPHI, PHI, pf_low, pr_high, pf_low, pr_high, dt
        )
        @test all(APHI[2:Ny, 2:Nx] .< 0.0)
    end

    @testset "compute_fluid_velocities!(): two-phase relative velocity and Galilean invariance" begin
        PHIX = fill(0.2, Ny1, Nx1)
        PHIY = fill(0.2, Ny1, Nx1)
        qxD = zeros(Ny1, Nx1)
        qyD = zeros(Ny1, Nx1)
        vx = zeros(Ny1, Nx1)
        vy = zeros(Ny1, Nx1)
        vxf = zeros(Ny1, Nx1)
        vyf = zeros(Ny1, Nx1)

        # 1. Synchronous limit: zero Darcy flux implies v_fluid == v_solid
        vx_val = 1.5e-9
        vy_val = -2.5e-9
        vx .= vx_val
        vy .= vy_val
        Erebus.compute_fluid_velocities!(PHIX, PHIY, qxD, qyD, vx, vy, vxf, vyf)
        @test isapprox(vxf[1:Ny1, 1:Nx], vx[1:Ny1, 1:Nx]; rtol=1e-12)
        @test isapprox(vyf[1:Ny, 1:Nx1], vy[1:Ny, 1:Nx1]; rtol=1e-12)

        # 2. Pure Darcy seepage with stationary solid matrix: v_fluid = qD / phi
        vx .= 0.0
        vy .= 0.0
        qxD .= 2.0e-10
        qyD .= 4.0e-10
        Erebus.compute_fluid_velocities!(PHIX, PHIY, qxD, qyD, vx, vy, vxf, vyf)
        for j in 1:Nx, i in 2:Ny
            @test isapprox(vxf[i, j], qxD[i, j] / PHIX[i, j]; rtol=1e-12)
        end
        for j in 2:Nx, i in 1:Ny
            @test isapprox(vyf[i, j], qyD[i, j] / PHIY[i, j]; rtol=1e-12)
        end

        # 3. Galilean invariance: translating both solid and fluid by V0 shifts fluid velocity by V0
        V0 = 5.0e-9
        vx_shifted = vx .+ V0
        vy_shifted = vy .+ V0
        vxf_shifted = zeros(Ny1, Nx1)
        vyf_shifted = zeros(Ny1, Nx1)
        Erebus.compute_fluid_velocities!(
            PHIX, PHIY, qxD, qyD, vx_shifted, vy_shifted, vxf_shifted, vyf_shifted
        )
        @test isapprox(vxf_shifted[1:Ny1, 1:Nx], vxf[1:Ny1, 1:Nx] .+ V0; rtol=1e-12)
        @test isapprox(vyf_shifted[1:Ny, 1:Nx1], vyf[1:Ny, 1:Nx1] .+ V0; rtol=1e-12)
    end

    @testset "compute_displacement_timestep(): CFL and porosity-rate constraints" begin
        dt = dt_longest

        # 1. Zero velocity and zero porosity rate: timestep preserved exactly
        dt_zero = Erebus.compute_displacement_timestep(
            zeros(Ny1, Nx1), zeros(Ny1, Nx1), zeros(Ny1, Nx1), zeros(Ny1, Nx1), dt, 0.0
        )
        @test dt_zero ≈ dt rtol=1e-12

        # 2. Solid velocity CFL constraint: dt * max(|vx|) <= dxymax * dx
        U_fast = 100.0 * (dxymax * dx / dt)
        vx_fast = fill(U_fast, Ny1, Nx1)
        dt_vx = Erebus.compute_displacement_timestep(
            vx_fast, zeros(Ny1, Nx1), zeros(Ny1, Nx1), zeros(Ny1, Nx1), dt, 0.0
        )
        @test dt_vx * U_fast ≈ dxymax * dx rtol=1e-12
        @test dt_vx < dt

        # 3. Fluid velocity CFL constraint: dt * max(|vxf|) <= dxymax * dx
        V_fast = 200.0 * (dxymax * dx / dt)
        vxf_fast = fill(V_fast, Ny1, Nx1)
        dt_vxf = Erebus.compute_displacement_timestep(
            zeros(Ny1, Nx1), zeros(Ny1, Nx1), vxf_fast, zeros(Ny1, Nx1), dt, 0.0
        )
        @test dt_vxf * V_fast ≈ dxymax * dx rtol=1e-12
        @test dt_vxf < dt

        # 4. Porosity rate constraint: dt * aphimax <= dphimax
        aphi_fast = 50.0 * (dphimax / dt)
        dt_phi = Erebus.compute_displacement_timestep(
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            dt,
            aphi_fast,
        )
        @test dt_phi * aphi_fast ≈ dphimax rtol=1e-12
        @test dt_phi < dt

        # 5. Monotonicity and positivity under combined random loads
        aphimax = rand(rgen)
        vx = rand(rgen, Ny1, Nx1)
        vy = rand(rgen, Ny1, Nx1)
        vxf = rand(rgen, Ny1, Nx1)
        vyf = rand(rgen, Ny1, Nx1)
        dtm = Erebus.compute_displacement_timestep(vx, vy, vxf, vyf, dt, aphimax)
        @test 0.0 < dtm <= dt
        @test dtm * maximum(abs, vx) <= dxymax * dx + 1e-12
        @test dtm * maximum(abs, vy) <= dxymax * dy + 1e-12
        @test dtm * maximum(abs, vxf) <= dxymax * dx + 1e-12
        @test dtm * maximum(abs, vyf) <= dxymax * dy + 1e-12
        @test dtm * aphimax <= dphimax + 1e-12
    end # testset "compute_displacement_timestep()"

    @testset "compute_stress_strainrate!(): continuum rotation and shear invariants" begin
        dtm = dt_longest
        vx = zeros(Ny1, Nx1)
        vy = zeros(Ny1, Nx1)
        ETA = fill(1.0e18, Ny, Nx)
        GGG = fill(1.0e10, Ny, Nx)
        ETAP = fill(1.0e18, Ny1, Nx1)
        GGGP = fill(1.0e10, Ny1, Nx1)
        SXX0 = zeros(Ny1, Nx1)
        SXY0 = zeros(Ny, Nx)
        EXX = zeros(Ny1, Nx1)
        EXY = zeros(Ny, Nx)
        SXX = zeros(Ny1, Nx1)
        SXY = zeros(Ny, Nx)
        DSXX = zeros(Ny1, Nx1)
        DSXY = zeros(Ny, Nx)
        EII = zeros(Ny1, Nx1)
        SII = zeros(Ny1, Nx1)

        # 1. Rigid Body Rotation Invariant: v = (-omega*y, omega*x) produces zero strain rate and zero stress
        omega = 1.0e-14
        for j in 1:Nx1, i in 1:Ny1
            vx[i, j] = -omega * (yp[i] - ycenter)
            vy[i, j] = omega * (xp[j] - xcenter)
        end
        Erebus.compute_stress_strainrate!(
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
        )
        @test all(isapprox.(EXY[2:Ny, 2:Nx], 0.0; atol=1e-18))
        @test all(isapprox.(EXX[2:Ny, 2:Nx], 0.0; atol=1e-18))
        @test all(isapprox.(EII[2:Ny, 2:Nx], 0.0; atol=1e-18))
        @test all(isapprox.(SII[2:Ny, 2:Nx], 0.0; atol=1e-4)) # stress zero within float tolerance

        # 2. Simple Shear Benchmark: v = (gamma_dot * y, 0) => EXY = 0.5 * gamma_dot, EXX = 0
        gamma_dot = 1.0e-12
        for j in 1:Nx1, i in 1:Ny1
            vx[i, j] = gamma_dot * yp[i]
            vy[i, j] = 0.0
        end
        Erebus.compute_stress_strainrate!(
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
        )
        for j in 2:Nx, i in 2:Ny
            @test isapprox(EXY[i, j], 0.5 * gamma_dot; rtol=1e-6)
            @test isapprox(EXX[i, j], 0.0; atol=1e-18)
        end

        # 3. Second Invariant Non-negativity
        @test all(EII .>= 0.0)
        @test all(SII .>= 0.0)
    end

    @testset "symmetrize_p_node_observables!(): Neumann symmetry and idempotence" begin
        SXX = 1.0e6 .* rand(rgen, Ny1, Nx1)
        APHI = 1.0e-12 .* rand(rgen, Ny1, Nx1)
        PHI = 0.1 .+ 0.1 .* rand(rgen, Ny1, Nx1)
        pr = 1.0e6 .* rand(rgen, Ny1, Nx1)
        pf = 0.5e6 .* rand(rgen, Ny1, Nx1)
        ps = zeros(Ny1, Nx1)

        # Record interior values before symmetrization
        interior_pr_before = copy(pr[2:Ny, 2:Nx])
        interior_pf_before = copy(pf[2:Ny, 2:Nx])

        Erebus.symmetrize_p_node_observables!(SXX, APHI, PHI, pr, pf, ps)

        # 1. Interior preservation: internal values must be unmodified
        @test isapprox(pr[2:Ny, 2:Nx], interior_pr_before; rtol=1e-12)
        @test isapprox(pf[2:Ny, 2:Nx], interior_pf_before; rtol=1e-12)

        # 2. Neumann zero-gradient boundary conditions on ghost nodes for all 5 observables
        for arr in (SXX, APHI, PHI, pr, pf)
            @test isapprox(arr[1, 2:Nx], arr[2, 2:Nx]; rtol=1e-12) # top
            @test isapprox(arr[Ny1, 2:Nx], arr[Ny, 2:Nx]; rtol=1e-12) # bottom
            @test isapprox(arr[:, 1], arr[:, 2]; rtol=1e-12) # left
            @test isapprox(arr[:, Nx1], arr[:, Nx]; rtol=1e-12) # right
        end

        # 3. Solid pressure formula exact verification
        ps_expected = (pr .- pf .* PHI) ./ (1.0 .- PHI)
        @test isapprox(ps, ps_expected; rtol=1e-12)

        # 4. Idempotence: applying symmetrization a second time yields identical results
        pr_pass1 = copy(pr)
        Erebus.symmetrize_p_node_observables!(SXX, APHI, PHI, pr, pf, ps)
        @test isapprox(pr, pr_pass1; rtol=1e-12)
    end

    @testset "positive_max(): non-negativity and upper bound axioms" begin
        A = rand(rgen, -100:0.1:100, 50, 50)
        B = rand(rgen, -100:0.1:100, 50, 50)
        R = zeros(50, 50)
        Erebus.positive_max!(A, B, R)

        # Axiom 1: Non-negativity everywhere
        @test all(r -> r >= 0.0, R)

        # Axiom 2: Upper bound (R >= A and R >= B everywhere)
        @test all(R .>= A)
        @test all(R .>= B)

        # Axiom 3: Negative saturation: strictly negative inputs produce exact zero
        A_neg = fill(-10.0, 10, 10)
        B_neg = fill(-25.0, 10, 10)
        R_neg = fill(99.0, 10, 10)
        Erebus.positive_max!(A_neg, B_neg, R_neg)
        @test all(iszero, R_neg)

        # Axiom 4: Positive maximum selection
        A_pos = [2.0 5.0; 8.0 1.0]
        B_pos = [3.0 1.0; 6.0 4.0]
        R_pos = zeros(2, 2)
        Erebus.positive_max!(A_pos, B_pos, R_pos)
        @test R_pos ≈ [3.0 5.0; 8.0 4.0] rtol=1e-12

        # Axiom 5: Commutativity / symmetry: positive_max(A, B) == positive_max(B, A)
        R_sym = zeros(50, 50)
        Erebus.positive_max!(B, A, R_sym)
        @test R ≈ R_sym rtol=1e-12
    end # testset "positive_max()"

    @testset "compute_nodal_adjustment!(): Drucker-Prager plasticity and pore pressure trigger" begin
        dt = dt_longest
        iplast = 1
        ETA = fill(1.0e20, Ny, Nx)
        ETA0 = fill(1.0e20, Ny, Nx)
        ETA5 = zeros(Ny, Nx)
        GGG = fill(1.0e10, Ny, Nx)
        SXX = fill(1.0e7, Ny1, Nx1)
        SXY = fill(1.0e7, Ny, Nx)
        COH = fill(1.0e7, Ny, Nx)
        TEN = fill(6.0e7, Ny, Nx)
        FRI = fill(0.6, Ny, Nx) # friction coeff
        YNY = zeros(Bool, Ny, Nx)
        YNY5 = zeros(Bool, Ny, Nx)
        DSY = zeros(Ny, Nx)
        YERRNOD = zeros(nplast)

        # 1. Stable Regime: very high confining pressure (pr = 100 MPa, pf = 0) suppresses yielding
        pr_high = fill(1.0e8, Ny1, Nx1)
        pf_zero = fill(0.0, Ny1, Nx1)
        complete = Erebus.compute_nodal_adjustment!(
            ETA,
            ETA0,
            ETA5,
            GGG,
            SXX,
            SXY,
            pr_high,
            pf_zero,
            COH,
            TEN,
            FRI,
            YNY,
            YNY5,
            YERRNOD,
            DSY,
            dt,
            iplast,
        )
        @test complete # no yield, iteration is complete
        @test all(.!YNY5[2:(Ny - 1), 2:(Nx - 1)]) # no yielding in interior
        @test all(
            isapprox.(
                ETA5[2:(Ny - 1), 2:(Nx - 1)], ETA0[2:(Ny - 1), 2:(Nx - 1)]; rtol=1e-12
            ),
        ) # viscosity unchanged

        # 2. Plastic Yielding Regime: elevated pore pressure (pf = 99.9 MPa) reduces effective pressure
        # P_eff = pr - pf = 0.1 MPa => yield strength sigma_y drops, triggering failure
        pf_elevated = fill(9.99e7, Ny1, Nx1)
        complete = Erebus.compute_nodal_adjustment!(
            ETA,
            ETA0,
            ETA5,
            GGG,
            SXX,
            SXY,
            pr_high,
            pf_elevated,
            COH,
            TEN,
            FRI,
            YNY,
            YNY5,
            YERRNOD,
            DSY,
            dt,
            iplast,
        )
        @test !complete # plastic yielding occurred, iteration not complete
        @test any(YNY5) # yielding flags triggered
        @test any(ETA5 .< ETA0) # apparent viscosity reduced to bring stress to envelope
        @test all(ETA5 .> 0.0) # apparent viscosity remains positive
    end

    @testset "finalize_plastic_iteration_pass!(): convergence vs stalled recovery" begin
        dt_init = dt_longest
        ETA = zeros(Ny, Nx)
        ETA5 = fill(1.0e17, Ny, Nx)
        ETA00 = fill(1.0e19, Ny, Nx)
        YNY = zeros(Bool, Ny, Nx)
        YNY5 = fill(true, Ny, Nx)
        YNY00 = fill(false, Ny, Nx)
        YNY_inv_ETA = zeros(Ny, Nx)

        # 1. Normal convergence pass (iplast % dtstep != 0)
        iplast_normal = 1
        dt_out = Erebus.finalize_plastic_iteration_pass!(
            ETA, ETA5, ETA00, YNY, YNY5, YNY00, YNY_inv_ETA, dt_init, iplast_normal
        )
        @test isapprox(dt_out, dt_init; rtol=1e-12) # dt preserved
        @test isapprox(ETA, ETA5; rtol=1e-12) # new viscosity adopted
        @test all(YNY .== YNY5)
        @test isapprox(YNY_inv_ETA, YNY5 ./ ETA5; rtol=1e-12)

        # 2. Stalled divergence recovery pass (iplast % dtstep == 0)
        iplast_stalled = dtstep
        dt_out_stalled = Erebus.finalize_plastic_iteration_pass!(
            ETA, ETA5, ETA00, YNY, YNY5, YNY00, YNY_inv_ETA, dt_init, iplast_stalled
        )
        @test isapprox(dt_out_stalled, dt_init * dtcoefdn; rtol=1e-12) # dt decelerated
        @test dt_out_stalled < dt_init
        @test isapprox(ETA, ETA00; rtol=1e-12) # viscosity rolled back
        @test all(YNY .== YNY00)
    end

    @testset "finalize_thermochemical_iteration_pass(): thermal relaxation step control" begin
        # 1. First iteration (titer == 1) with sub-threshold temperature change: dt unchanged
        dt_sub = Erebus.finalize_thermochemical_iteration_pass(DTmax * 0.5, dt_longest, 1)
        @test dt_sub ≈ dt_longest rtol=1e-12

        # 2. First iteration (titer == 1) with excessive temperature change: dt scaled proportionally
        maxDT_large = DTmax * 4.0
        dt_cut = Erebus.finalize_thermochemical_iteration_pass(maxDT_large, dt_longest, 1)
        @test dt_cut ≈ dt_longest * (DTmax / maxDT_large) rtol=1e-12
        @test dt_cut < dt_longest

        # 3. Subsequent iterations (titer > 1): dt is preserved regardless of DT magnitude
        dt_iter2 = Erebus.finalize_thermochemical_iteration_pass(maxDT_large, dt_longest, 2)
        dt_iter3 = Erebus.finalize_thermochemical_iteration_pass(maxDT_large, dt_longest, 3)
        @test dt_iter2 ≈ dt_longest rtol=1e-12
        @test dt_iter3 ≈ dt_longest rtol=1e-12

        # 4. Positivity and monotonicity: larger excess DT produces strictly smaller dt in titer 1
        dt_cut_larger = Erebus.finalize_thermochemical_iteration_pass(
            DTmax * 8.0, dt_longest, 1
        )
        @test dt_cut_larger < dt_cut
        @test dt_cut_larger > 0.0
    end # testset "finalize_thermochemical_iteration_pass()"

    @testset "compute_thermochemical_iteration_outcome: convergence decision boundaries" begin
        # 1. Converged case: small pressure error and past titer threshold (titer > 2)
        pf_converged = rand(rgen, Ny1, Nx1)
        pf0 = copy(pf_converged) # zero pressure error
        DMP_active = fill(1.0e-5, Ny1, Nx1)
        @test Erebus.compute_thermochemical_iteration_outcome(
            DMP_active, pf_converged, pf0, 3
        ) == true

        # 2. Converged case: small pressure error and zero mass reaction (DMP <= 0) even at titer 1
        DMP_zero = zeros(Ny1, Nx1)
        @test Erebus.compute_thermochemical_iteration_outcome(
            DMP_zero, pf_converged, pf0, 1
        ) == true

        # 3. Non-converged case: excessive pressure error (pferrcur >= pferrmax)
        pf_diverged = copy(pf0)
        pf_diverged[2, 2] += pferrmax * 2.0
        @test Erebus.compute_thermochemical_iteration_outcome(
            DMP_active, pf_diverged, pf0, 3
        ) == false
        @test Erebus.compute_thermochemical_iteration_outcome(
            DMP_zero, pf_diverged, pf0, 1
        ) == false

        # 4. Non-converged case: active reaction (DMP > 0) at early iterations (titer <= 2)
        @test Erebus.compute_thermochemical_iteration_outcome(
            DMP_active, pf_converged, pf0, 1
        ) == false
        @test Erebus.compute_thermochemical_iteration_outcome(
            DMP_active, pf_converged, pf0, 2
        ) == false
    end # testset "compute_thermochemical_iteration_outcome"

    @testset "assemble_thermal_lse!(): operator symmetry, diagonal dominance, isothermal state" begin
        dt = 1.0e11
        tk1 = fill(300.0, Ny1, Nx1)
        RHOCP = fill(3.3e6, Ny1, Nx1)
        KX = fill(3.0, Ny1, Nx1)
        KY = fill(3.0, Ny1, Nx1)
        HR = zeros(Ny1, Nx1)
        HA = zeros(Ny1, Nx1)
        HS = zeros(Ny1, Nx1)
        DHP = zeros(Ny1, Nx1)
        RT = zeros(Ny1 * Nx1)

        LT = Erebus.assemble_thermal_lse!(tk1, RHOCP, KX, KY, HR, HA, HS, DHP, RT, dt)

        # 1. Symmetry of interior thermal conductance with spatially varying conductivities
        KX_var = zeros(Ny1, Nx1)
        KY_var = zeros(Ny1, Nx1)
        for j in 1:Nx1, i in 1:Ny1
            KX_var[i, j] = 2.0 + 1.0 * sin(xp[j] / xsize * π)
            KY_var[i, j] = 2.0 + 1.0 * cos(yp[i] / ysize * π)
        end
        LT_var = Erebus.assemble_thermal_lse!(
            tk1, RHOCP, KX_var, KY_var, HR, HA, HS, DHP, zeros(Ny1 * Nx1), dt
        )
        for j in 2:(Nx1 - 2), i in 2:(Ny1 - 2)
            gk = (j - 1) * Ny1 + i
            gk_right = gk + Ny1
            gk_down = gk + 1
            @test isapprox(LT_var[gk, gk_right], LT_var[gk_right, gk]; rtol=1e-12)
            @test isapprox(LT_var[gk, gk_down], LT_var[gk_down, gk]; rtol=1e-12)
        end

        # 2. Strict diagonal dominance & M-matrix property on interior rows
        for j in 2:Nx, i in 2:Ny
            gk = (j - 1) * Ny1 + i
            diag_val = LT[gk, gk]
            off_diag_sum =
                abs(LT[gk, gk - Ny1]) +
                abs(LT[gk, gk - 1]) +
                abs(LT[gk, gk + 1]) +
                abs(LT[gk, gk + Ny1])
            @test diag_val > off_diag_sum # strictly dominant
            @test isapprox(diag_val - off_diag_sum, RHOCP[i, j] / dt; rtol=1e-12)
        end

        # 3. Isothermal Stationary Invariant: uniform T0 with zero heat sources satisfies LT * T0 == RT
        T_vec = fill(300.0, Nx1 * Ny1)
        LT_dense = collect(LT)
        L_dot_T = LT_dense * T_vec
        for j in 2:Nx, i in 2:Ny
            gk = (j - 1) * Ny1 + i
            @test isapprox(L_dot_T[gk], RT[gk]; rtol=1e-12)
        end

        # 4. Linear superposition of metric heat source
        Q_metric = fill(1.5e-4, Ny1, Nx1)
        RT_metric = zeros(Ny1 * Nx1)
        LT_m = Erebus.assemble_thermal_lse!(
            tk1, RHOCP, KX, KY, HR, HA, HS, DHP, RT_metric, dt; Q_metric=Q_metric
        )
        for j in 2:Nx, i in 2:Ny
            gk = (j - 1) * Ny1 + i
            @test isapprox(RT_metric[gk], RT[gk] + Q_metric[i, j]; rtol=1e-12)
        end
    end

    @testset "poroelastic hydromechanical coupling" begin
        Ny, Nx = Erebus.Ny, Erebus.Nx
        Ny1, Nx1 = Erebus.Ny1, Erebus.Nx1
        dt = 10.0

        ETA = fill(1e22, Ny, Nx)
        ETAP = fill(1e22, Ny1, Nx1)
        GGG = fill(1e10, Ny, Nx)
        GGGP = fill(1e10, Ny1, Nx1)
        SXY0 = zeros(Ny, Nx)
        SXX0 = zeros(Ny, Nx)
        RHOX = fill(3000.0, Ny1, Nx1)
        RHOY = fill(3000.0, Ny1, Nx1)
        RHOFX = fill(1000.0, Ny1, Nx1)
        RHOFY = fill(1000.0, Ny1, Nx1)
        RX = fill(1e-3 / 1e-13, Ny1, Nx1)
        RY = fill(1e-3 / 1e-13, Ny1, Nx1)
        PHI = fill(0.1, Ny1, Nx1)
        ETAPHI = fill(1e22 / 0.1, Ny1, Nx1)
        BETAPHI = fill(0.1 / 1e10, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = zeros(Ny1, Nx1)
        pr0 = zeros(Ny1, Nx1)
        pf0 = fill(1e6, Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)
        R = zeros(Nx1*Ny1*6)

        betasolid = 2.5e-11
        betafluid = 4.0e-10
        L = Erebus.assemble_hydromechanical_lse!(
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
            dt,
            R;
            betasolid=betasolid,
            betafluid=betafluid,
        )

        # Independent theoretical coefficients (calculated without calling solver functions)
        bd_th = (0.1 / 1e10 + betasolid) / (1.0 - 0.1)
        kbw_th = 1.0 - betasolid / bd_th
        ksk_th = (bd_th - betasolid) / ((bd_th - betasolid) + 0.1 * (betafluid - betasolid))
        C_expected = -Kcont * (inv(1e22 / 0.1) / (1.0 - 0.1) + bd_th * kbw_th / dt)
        D_pm_expected = Kcont * (inv(1e22 / 0.1) / (1.0 - 0.1) + bd_th / dt)
        D_pf_expected =
            Kcont * (inv(1e22 / 0.1) / (1.0 - 0.1) + bd_th * kbw_th / ksk_th / dt)

        # Verify each matrix block independently against theoretical values
        for j in 4:(Nx - 2), i in 4:(Ny - 2)
            kvx = ((j-1)*Ny1 + i-1) * 6 + 1
            kpm = kvx + 2
            kpf = kvx + 5
            @test L[kpm, kpf] ≈ C_expected rtol=1e-10
            @test L[kpf, kpm] ≈ C_expected rtol=1e-10
            @test L[kpm, kpm] ≈ D_pm_expected rtol=1e-10
            @test L[kpf, kpf] ≈ D_pf_expected rtol=1e-10
        end

        # Verify solution solvability and finiteness
        prob = LinearProblem(L, R)
        sol = solve(prob, UMFPACKFactorization())
        @test !any(isnan, sol.u)
        @test !any(isinf, sol.u)
    end # testset "poroelastic hydromechanical coupling"

    @testset "Terzaghi 1D consolidation numerical simulation verification" begin
        Ny, Nx = Erebus.Ny, Erebus.Nx
        Ny1, Nx1 = Erebus.Ny1, Erebus.Nx1
        dy = Erebus.dy

        # Height between draining boundary anchors i=2 and i=Ny
        H = (Ny - 2) * dy
        k_perm = 1e-13
        eta_f = 1e-3
        betasolid = 2.5e-11
        betafluid = 4.0e-10
        G_p = 1e10
        phi_0 = 0.1
        beta_phi = phi_0 / G_p

        bd = Erebus.compute_drained_compressibility(beta_phi, phi_0, betasolid)
        kbw = Erebus.compute_biot_willis_coefficient(bd, betasolid)
        ksk = Erebus.compute_skempton_coefficient(bd, phi_0, betasolid, betafluid)
        S = bd * kbw / ksk
        c_v = k_perm / (eta_f * S)

        u0 = 1.0e6
        dt = 2.0e7 # Timestep [s] scaled with column height (H^2/cv)

        ETA = fill(1e25, Ny, Nx)
        ETAP = fill(1e25, Ny1, Nx1)
        GGG = fill(1e10, Ny, Nx)
        GGGP = fill(1e10, Ny1, Nx1)
        SXY0 = zeros(Ny, Nx)
        SXX0 = zeros(Ny, Nx)
        RHOX = fill(3000.0, Ny1, Nx1)
        RHOY = fill(3000.0, Ny1, Nx1)
        RHOFX = fill(1000.0, Ny1, Nx1)
        RHOFY = fill(1000.0, Ny1, Nx1)
        RX = fill(eta_f / k_perm, Ny1, Nx1)
        RY = fill(eta_f / k_perm, Ny1, Nx1)
        PHI = fill(phi_0, Ny1, Nx1)
        ETAPHI = fill(1e25, Ny1, Nx1)
        BETAPHI = fill(beta_phi, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = zeros(Ny1, Nx1)
        pr0 = zeros(Ny1, Nx1)
        pf0 = fill(u0, Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)
        R = zeros(Nx1*Ny1*6)

        vx = zeros(Ny1, Nx1)
        vy = zeros(Ny1, Nx1)
        pr = zeros(Ny1, Nx1)
        qxD = zeros(Ny1, Nx1)
        qyD = zeros(Ny1, Nx1)
        pf = zeros(Ny1, Nx1)

        nsteps = 3
        t_total = nsteps * dt

        # Step Erebus numerical solver forward in time
        for step in 1:nsteps
            R .= 0.0
            L = Erebus.assemble_hydromechanical_lse!(
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
                dt,
                R;
                betasolid=betasolid,
                betafluid=betafluid,
            )
            prob = LinearProblem(L, R)
            sol = solve(prob, UMFPACKFactorization())
            Erebus.process_hydromechanical_solution!(sol.u, vx, vy, pr, qxD, qyD, pf)
            pf0 .= pf
        end

        # Analytical 1D consolidation Fourier series solution
        function analytical_2drain(y, t, H_col, c_coeff, p0; nterms=100)
            val = 0.0
            for m in 0:nterms
                M = (2m + 1) * pi
                val +=
                    (4.0 * p0 / M) * sin(M * y / H_col) * exp(-M^2 * c_coeff * t / H_col^2)
            end
            return val
        end

        # Verify numerical solution against analytical solution across column depth
        # Account for psurface boundary condition at draining anchors
        mid_col = div(Nx1, 2)
        for i in 3:(Ny - 1)
            y = (i - 2) * dy
            u_ana = analytical_2drain(y, t_total, H, c_v, u0 - psurface) + psurface
            u_num = pf[i, mid_col]
            rel_err = abs(u_num - u_ana) / u0
            @test rel_err < 0.035
        end

        # Physical Invariant 1: Peak excess pore pressure decays monotonically over consolidation
        @test maximum(pf[:, mid_col]) < u0

        # Physical Invariant 2: Positivity and bounded pressure profile inside column
        @test all(pf[3:(Ny - 1), mid_col] .>= psurface - 1e-6)
        @test all(isfinite, pf)

        # Physical Invariant 3: Spatial symmetry of pressure profile across draining boundaries
        top_idx = 3
        bot_idx = Ny - 1
        @test isapprox(pf[top_idx, mid_col], pf[bot_idx, mid_col]; rtol=0.05)
    end # testset "Terzaghi 1D consolidation numerical simulation verification"

    @testset "assemble_hydromechanical_lse!() dynamic hydrofracture integration" begin
        eta_f = 1.0e-3
        k_perm = 1.0e-12
        k_frac_max = 1.0e-9
        kappa_frac = 100.0
        gamma_frac = 1.0
        sigma_t_val = 5.0e6

        ETA = fill(1e22, Ny, Nx)
        ETAP = fill(1e22, Ny1, Nx1)
        GGG = fill(1e10, Ny, Nx)
        GGGP = fill(1e10, Ny1, Nx1)
        SXY0 = zeros(Ny, Nx)
        SXX0 = zeros(Ny, Nx)
        RHOX = fill(3000.0, Ny1, Nx1)
        RHOY = fill(3000.0, Ny1, Nx1)
        RHOFX = fill(1000.0, Ny1, Nx1)
        RHOFY = fill(1000.0, Ny1, Nx1)
        RX = fill(eta_f / k_perm, Ny1, Nx1)
        RY = fill(eta_f / k_perm, Ny1, Nx1)
        KX = fill(k_perm, Ny1, Nx1)
        KY = fill(k_perm, Ny1, Nx1)
        PHI = fill(0.1, Ny1, Nx1)
        ETAPHI = fill(1e25, Ny1, Nx1)
        BETAPHI = fill(1e-10, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = zeros(Ny1, Nx1)
        pr0 = zeros(Ny1, Nx1)
        pf0 = zeros(Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)
        dt = 1e8
        R = zeros(Nx1*Ny1*6)

        pr = fill(1.0e7, Ny1, Nx1)
        TEN = fill(sigma_t_val, Ny, Nx)

        L_base = Erebus.assemble_hydromechanical_lse!(
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
            dt,
            R;
            hydrofracture=false,
        )

        pf_intact = fill(5.0e6, Ny1, Nx1)
        L_intact = Erebus.assemble_hydromechanical_lse!(
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
            dt,
            R;
            hydrofracture=true,
            pr=pr,
            pf=pf_intact,
            TEN=TEN,
            KX=KX,
            KY=KY,
            kappa_frac=kappa_frac,
            gamma_frac=gamma_frac,
            k_frac_max=k_frac_max,
        )
        i_test, j_test = 5, 5
        kqx = ((j_test - 1)*Ny1 + i_test - 1) * 6 + 4
        kqy = kqx + 1
        @test L_intact[kqx, kqx] ≈ L_base[kqx, kqx]
        @test L_intact[kqy, kqy] ≈ L_base[kqy, kqy]

        pf_frac = fill(2.0e7, Ny1, Nx1)
        L_frac = Erebus.assemble_hydromechanical_lse!(
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
            dt,
            R;
            hydrofracture=true,
            pr=pr,
            pf=pf_frac,
            TEN=TEN,
            KX=KX,
            KY=KY,
            kappa_frac=kappa_frac,
            gamma_frac=gamma_frac,
            k_frac_max=k_frac_max,
        )
        expected_keff = k_perm * 101.0
        expected_rx = RX[i_test, j_test] * (k_perm / expected_keff)
        @test L_frac[kqx, kqx] < L_base[kqx, kqx]
        @test L_frac[kqx, kqx] ≈ expected_rx rtol=1e-12
        @test L_frac[kqy, kqy] ≈ expected_rx rtol=1e-12

        pf_extreme = fill(2.0e8, Ny1, Nx1)
        L_extreme = Erebus.assemble_hydromechanical_lse!(
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
            dt,
            R;
            hydrofracture=true,
            pr=pr,
            pf=pf_extreme,
            TEN=TEN,
            KX=KX,
            KY=KY,
            kappa_frac=kappa_frac,
            gamma_frac=gamma_frac,
            k_frac_max=k_frac_max,
        )
        expected_ceiling_rx = eta_f / k_frac_max
        @test L_extreme[kqx, kqx] ≈ expected_ceiling_rx rtol=1e-12
        @test L_extreme[kqy, kqy] ≈ expected_ceiling_rx rtol=1e-12

        HS_intact = zeros(Ny, Nx)
        HS_frac = zeros(Ny, Nx)
        SXY_test = fill(1e5, Ny, Nx)
        SXX_test = fill(1e5, Ny1, Nx1)
        qxD_test = fill(1e-6, Ny1, Nx1)
        qyD_test = fill(1e-6, Ny1, Nx1)
        Erebus.compute_shear_heating!(
            HS_intact,
            ETA,
            SXY_test,
            ETAP,
            SXX_test,
            RX,
            RY,
            qxD_test,
            qyD_test,
            PHI,
            ETAPHI,
            pr,
            pf_intact;
            hydrofracture=true,
            TEN=TEN,
            KX=KX,
            KY=KY,
            kappa_frac=kappa_frac,
            gamma_frac=gamma_frac,
            k_frac_max=k_frac_max,
        )
        Erebus.compute_shear_heating!(
            HS_frac,
            ETA,
            SXY_test,
            ETAP,
            SXX_test,
            RX,
            RY,
            qxD_test,
            qyD_test,
            PHI,
            ETAPHI,
            pr,
            pf_frac;
            hydrofracture=true,
            TEN=TEN,
            KX=KX,
            KY=KY,
            kappa_frac=kappa_frac,
            gamma_frac=gamma_frac,
            k_frac_max=k_frac_max,
        )
        @test HS_frac[5, 5] < HS_intact[5, 5]
    end # testset "assemble_hydromechanical_lse!() dynamic hydrofracture integration"

    @testset "assemble_hydromechanical_lse!() Darcy thermal buoyancy integration" begin
        g_accel = 10.0
        rho0_f = 1000.0
        alpha_f = 2.0e-4
        delta_T = 50.0
        eta_f = 1.0e-3
        k_perm = 1.0e-13

        ETA = fill(1e22, Ny, Nx)
        ETAP = fill(1e22, Ny1, Nx1)
        GGG = fill(1e10, Ny, Nx)
        GGGP = fill(1e10, Ny1, Nx1)
        SXY0 = zeros(Ny, Nx)
        SXX0 = zeros(Ny, Nx)
        RHOX = fill(3000.0, Ny1, Nx1)
        RHOY = fill(3000.0, Ny1, Nx1)
        PHI = fill(0.1, Ny1, Nx1)
        ETAPHI = fill(1e25, Ny1, Nx1)
        BETAPHI = fill(1e-10, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = fill(g_accel, Ny1, Nx1)
        pr0 = zeros(Ny1, Nx1)
        pf0 = zeros(Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)
        dt = 1e8
        R_const = zeros(Nx1*Ny1*6)
        R_buoy = zeros(Nx1*Ny1*6)

        RX = fill(eta_f / k_perm, Ny1, Nx1)
        RY = fill(eta_f / k_perm, Ny1, Nx1)

        RHOFX_const = fill(rho0_f, Ny1, Nx1)
        RHOFY_const = fill(rho0_f, Ny1, Nx1)
        Erebus.assemble_hydromechanical_lse!(
            ETA,
            ETAP,
            GGG,
            GGGP,
            SXY0,
            SXX0,
            RHOX,
            RHOY,
            RHOFX_const,
            RHOFY_const,
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
            dt,
            R_const,
        )

        rho_f_hot = rho0_f * (1.0 - alpha_f * delta_T)
        RHOFX_buoy = copy(RHOFX_const)
        RHOFY_buoy = copy(RHOFY_const)
        RHOFY_buoy[Ny - 2, 5] = rho_f_hot

        Erebus.assemble_hydromechanical_lse!(
            ETA,
            ETAP,
            GGG,
            GGGP,
            SXY0,
            SXX0,
            RHOX,
            RHOY,
            RHOFX_buoy,
            RHOFY_buoy,
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
            dt,
            R_buoy,
        )

        kqy_test = ((5 - 1)*Ny1 + (Ny-2) - 1) * 6 + 5
        delta_R = R_const[kqy_test] - R_buoy[kqy_test]
        expected_delta_R = (rho0_f - rho_f_hot) * g_accel
        @test delta_R ≈ expected_delta_R rtol=1e-12
        @test delta_R ≈ (rho0_f * alpha_f * delta_T * g_accel) rtol=1e-12
    end # testset "assemble_hydromechanical_lse!() Darcy thermal buoyancy integration"
end

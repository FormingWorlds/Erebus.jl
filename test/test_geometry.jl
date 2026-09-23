@testset "Geometry" begin
    @testset "setup_staggered_grid_geometry(): metric monotonicity and staggered topology" begin
        # 1. Grid spacing positivity and scale
        @test dx > 0.0
        @test dy > 0.0
        @test isapprox(dx, xsize / (Nx - 1); rtol=1e-12)
        @test isapprox(dy, ysize / (Ny - 1); rtol=1e-12)

        # 2. Strict monotonicity and uniform spacing
        for j in 1:(Nx - 1)
            @test x[j + 1] > x[j]
            @test isapprox(x[j + 1] - x[j], dx; rtol=1e-12)
        end
        for i in 1:(Ny - 1)
            @test y[i + 1] > y[i]
            @test isapprox(y[i + 1] - y[i], dy; rtol=1e-12)
        end

        # 3. Domain coverage: basic grid spans [0, xsize] x [0, ysize]
        @test isapprox(x[1], 0.0; atol=1e-12)
        @test isapprox(x[end], xsize; rtol=1e-12)
        @test isapprox(y[1], 0.0; atol=1e-12)
        @test isapprox(y[end], ysize; rtol=1e-12)

        # 4. Staggered grid topology and mid-cell placement
        # Pressure nodes are centered between basic nodes for interior cells
        for j in 2:Nx
            @test isapprox(xp[j], 0.5 * (x[j] + x[j - 1]); rtol=1e-12)
        end
        for i in 2:Ny
            @test isapprox(yp[i], 0.5 * (y[i] + y[i - 1]); rtol=1e-12)
        end

        # Ghost node extensions: pressure grid has half-cell padding beyond physical domain
        @test isapprox(xp[1], -0.5 * dx; rtol=1e-12)
        @test isapprox(xp[end], xsize + 0.5 * dx; rtol=1e-12)
        @test isapprox(yp[1], -0.5 * dy; rtol=1e-12)
        @test isapprox(yp[end], ysize + 0.5 * dy; rtol=1e-12)

        # 5. Velocity node offsets
        # Vx nodes: horizontally aligned with x, vertically staggered by -dy/2
        @test isapprox(xvx[1], 0.0; atol=1e-12)
        @test isapprox(yvx[1], -0.5 * dy; rtol=1e-12)
        @test isapprox(yvx[end], ysize + 0.5 * dy; rtol=1e-12)

        # Vy nodes: horizontally staggered by -dx/2, vertically aligned with y
        @test isapprox(xvy[1], -0.5 * dx; atol=1e-12)
        @test isapprox(xvy[end], xsize + 0.5 * dx; rtol=1e-12)
        @test isapprox(yvy[1], 0.0; atol=1e-12)

        # 6. Grid dimension consistency
        @test length(x) == Nx
        @test length(y) == Ny
        @test length(xvx) == Nx1
        @test length(yvx) == Ny1
        @test length(xvy) == Nx1
        @test length(yvy) == Ny1
        @test length(xp) == Nx1
        @test length(yp) == Ny1
        @test Nx1 == Nx + 1
        @test Ny1 == Ny + 1

        # 7. Discrimination guard: staggered node is not coincident with basic node
        @test abs(xp[1] - x[1]) >= 0.49 * dx
        @test abs(yp[1] - y[1]) >= 0.49 * dy
    end # testset "setup_staggered_grid_geometry()"

    @testset "setup_staggered_grid_properties(): tensor dimension conformance" begin
        props = Erebus.setup_staggered_grid_properties()
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
            DT,
            DT0,
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
            FI,
            DMP,
            DHP,
            XWS,
        ) = props

        # 1. Basic node tensor dimensions: (Ny, Nx)
        basic_tensors = [ETA, ETA0, GGG, EXY, SXY, SXY0, wyx, COH, TEN, FRI, YNY]
        for t in basic_tensors
            @test size(t) == (Ny, Nx)
            @test all(iszero, t)
        end

        # 2. Velocity and Pressure staggered node tensor dimensions: (Ny1, Nx1)
        staggered_tensors = [
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
            DT,
            DT0,
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
            FI,
            DMP,
            DHP,
            XWS,
        ]
        for t in staggered_tensors
            @test size(t) == (Ny1, Nx1)
            @test all(iszero, t)
        end
    end # testset "setup_staggered_grid_properties()"

    @testset "setup_staggered_grid_properties(randomized=true): broadcasting and field variance" begin
        props_rand = Erebus.setup_staggered_grid_properties(; randomized=true)
        SXX_rand = props_rand[40]
        SXX0_rand = props_rand[41]
        FI_rand = props_rand[60]

        # 1. Non-constancy and finite values across randomized fields
        @test all(isfinite, SXX_rand)
        @test all(isfinite, SXX0_rand)
        @test all(isfinite, FI_rand)

        # 2. SXX / SXX0 broadcast bounds [-1e3, 1e3]
        @test -1.0e3 <= minimum(SXX_rand) < 0.0 < maximum(SXX_rand) <= 1.0e3
        @test -1.0e3 <= minimum(SXX0_rand) < 0.0 < maximum(SXX0_rand) <= 1.0e3

        # 3. Discrimination guard: FI must be non-constant in [-1e2, 1e2]
        # (catches regression of the 2e2.=1e2 bug which produced a constant field of 100.0)
        @test minimum(FI_rand) < 0.0 < maximum(FI_rand)
        @test minimum(FI_rand) != maximum(FI_rand)
        @test abs(maximum(FI_rand) - minimum(FI_rand)) > 10.0
    end

    @testset "setup_staggered_grid_properties_helpers(): helper array dimensions" begin
        helpers = Erebus.setup_staggered_grid_properties_helpers()
        (ETA5, ETA00, YNY5, YNY00, YNY_inv_ETA, DSXY, DSY, EII, SII, DSXX, tk0) = helpers

        # Basic helpers (Ny, Nx)
        for h in [ETA5, ETA00, DSXY, DSY]
            @test size(h) == (Ny, Nx)
            @test eltype(h) == Float64
            @test all(iszero, h)
        end
        for h in [YNY5, YNY00]
            @test size(h) == (Ny, Nx)
            @test eltype(h) == Bool
            @test all(!, h)
        end

        # Staggered helpers (Ny1, Nx1)
        for h in [EII, SII, DSXX, tk0]
            @test size(h) == (Ny1, Nx1)
            @test eltype(h) == Float64
            @test all(iszero, h)
        end
    end # testset "setup_staggered_grid_properties_helpers()"

    @testset "grid_vector(): 4-point cell stencil ordering" begin
        grid = rand(rgen, 8, 8)

        # 1. Stencil ordering: [top-left (i, j), bottom-left (i+1, j), top-right (i, j+1), bottom-right (i+1, j+1)]
        v11 = Erebus.grid_vector(1, 1, grid)
        @test v11 ≈ @SVector [grid[1, 1], grid[2, 1], grid[1, 2], grid[2, 2]]
        @test length(v11) == 4

        # 2. Interior cell ordering consistency
        v34 = Erebus.grid_vector(3, 4, grid)
        @test v34[1] ≈ grid[3, 4]
        @test v34[2] ≈ grid[4, 4]
        @test v34[3] ≈ grid[3, 5]
        @test v34[4] ≈ grid[4, 5]

        # 3. Positivity preservation
        grid_pos = rand(rgen, 5, 5) .+ 1.0
        v_pos = Erebus.grid_vector(2, 2, grid_pos)
        @test all(v_pos .> 0.0)
    end # testset "grid_vector()"

    @testset "grid_average(): partition of unity and convex hull" begin
        # 1. Constant field partition of unity: average of constant c is c identically
        const_val = 42.5
        grid_const = fill(const_val, 6, 6)
        @test isapprox(Erebus.grid_average(1, 1, grid_const), const_val; rtol=1e-12)
        @test isapprox(Erebus.grid_average(3, 3, grid_const), const_val; rtol=1e-12)

        # 2. Convex hull invariant: min(cell) <= average <= max(cell)
        grid_rand = rand(rgen, 6, 6)
        for j in 1:5, i in 1:5
            cell_vals = [
                grid_rand[i, j],
                grid_rand[i + 1, j],
                grid_rand[i, j + 1],
                grid_rand[i + 1, j + 1],
            ]
            avg = Erebus.grid_average(i, j, grid_rand)
            @test minimum(cell_vals) <= avg <= maximum(cell_vals)
        end

        # 3. Linearity: avg(a*G1 + b*G2) == a*avg(G1) + b*avg(G2)
        g1 = rand(rgen, 6, 6)
        g2 = rand(rgen, 6, 6)
        a, b = 2.5, -1.8
        comb = a .* g1 .+ b .* g2
        avg_comb = Erebus.grid_average(2, 2, comb)
        expected_comb =
            a * Erebus.grid_average(2, 2, g1) + b * Erebus.grid_average(2, 2, g2)
        @test isapprox(avg_comb, expected_comb; rtol=1e-12)

        # 4. Discrete symmetry under transpose: for a symmetric matrix, diagonal cell average is symmetric
        sym_mat = [1.0 2.0; 2.0 3.0]
        @test isapprox(
            Erebus.grid_average(1, 1, sym_mat), 0.25 * (1.0 + 2.0 + 2.0 + 3.0); rtol=1e-12
        )
    end # testset "grid_average()"

    @testset "apply_insulating_boundary_conditions!(): zero Neumann flux and idempotency" begin
        # 1. Zero Neumann flux: normal derivative across boundary is zero (ghost cell == adjacent interior cell)
        ny_test, nx_test = 8, 8
        t_field = rand(rgen, ny_test, nx_test)
        t_orig = copy(t_field)

        Erebus.apply_insulating_boundary_conditions!(t_field)

        # Top and bottom zero flux: dT/dy = 0
        @test t_field[1, 2:(nx_test - 1)] ≈ t_field[2, 2:(nx_test - 1)]
        @test t_field[ny_test, 2:(nx_test - 1)] ≈ t_field[ny_test - 1, 2:(nx_test - 1)]

        # Left and right zero flux: dT/dx = 0
        @test t_field[:, 1] ≈ t_field[:, 2]
        @test t_field[:, nx_test] ≈ t_field[:, nx_test - 1]

        # 2. Interior preservation: physical interior (3:(ny-2), 3:(nx-2)) must not be modified
        @test t_field[3:(ny_test - 2), 3:(nx_test - 2)] ≈
            t_orig[3:(ny_test - 2), 3:(nx_test - 2)]

        # 3. Idempotency: BC(BC(T)) == BC(T)
        t_double = copy(t_field)
        Erebus.apply_insulating_boundary_conditions!(t_double)
        @test t_double ≈ t_field
    end # testset "apply_insulating_boundary_conditions!()"

    @testset "Marker out-of-plane weight" begin
        # 1. Sanity anchor: uniform-density markers on the shipped seeding lattice inside a disk of radius R.
        # Sum of rho * marker_area * L(r) equals (4/3)*pi*R^3*rho to lattice tolerance 2*dxm/R.
        coords = GridCoordinates(65, 65; xsize=140_000.0, ysize=140_000.0, Nxmc=4, Nymc=4)
        R_planet = 50_000.0
        rho_ref = 3000.0
        xc, yc = coords.xcenter, coords.ycenter
        dxm, dym = coords.dxm, coords.dym
        tol_lattice = 2.0 * dxm / R_planet

        sum_3d_disk = 0.0
        for jm in 1:coords.Nxm, im in 1:coords.Nym
            xm = dxm / 2.0 + (jm - 1) * dxm
            ym = dym / 2.0 + (im - 1) * dym
            r = sqrt((xm - xc)^2 + (ym - yc)^2)
            if r <= R_planet
                sum_3d_disk +=
                    rho_ref *
                    marker_area(coords) *
                    marker_out_of_plane_length(xm, ym, xc, yc)
            end
        end
        exact_3d_disk = (4.0 / 3.0) * pi * R_planet^3 * rho_ref
        @test isapprox(sum_3d_disk, exact_3d_disk; rtol=tol_lattice)

        # 2. Shell 0.9R..R: sum equals (4/3)*pi*(R^3 - (0.9R)^3)*rho to lattice tolerance,
        # and the ratio to 2D shell mass equals (4/3)*(R^3 - r^3)/(R^2 - r^2) = 1.90175... * R (~1.902R to 1%).
        # This excludes both uniform (4/3)R (-29.9%) and uniform 2R (+5.2%) errors.
        r_inner = 0.9 * R_planet
        sum_3d_shell = 0.0
        sum_2d_shell = 0.0
        for jm in 1:coords.Nxm, im in 1:coords.Nym
            xm = dxm / 2.0 + (jm - 1) * dxm
            ym = dym / 2.0 + (im - 1) * dym
            r = sqrt((xm - xc)^2 + (ym - yc)^2)
            if r_inner <= r <= R_planet
                sum_3d_shell +=
                    rho_ref *
                    marker_area(coords) *
                    marker_out_of_plane_length(xm, ym, xc, yc)
                sum_2d_shell += rho_ref * marker_area(coords)
            end
        end
        exact_3d_shell = (4.0 / 3.0) * pi * (R_planet^3 - r_inner^3) * rho_ref
        @test isapprox(sum_3d_shell, exact_3d_shell; rtol=tol_lattice)

        ratio_shell = sum_3d_shell / sum_2d_shell
        expected_ratio = (4.0 / 3.0) * (R_planet^3 - r_inner^3) / (R_planet^2 - r_inner^2)
        @test isapprox(ratio_shell, expected_ratio; rtol=0.01)
        @test isapprox(ratio_shell / R_planet, 1.902; rtol=0.01)

        # 3. Ring at r = 0.5R: the 3D mass equals the 2D mass times R to 1e-12 (exact).
        # A uniform 2R error gives 2R, while a uniform (4/3)R error gives (4/3)R.
        n_ring = 100
        r_ring = 0.5 * R_planet
        theta_ring = range(0.0, 2.0 * pi; length=n_ring + 1)[1:n_ring]
        ring_3d_mass = 0.0
        ring_2d_mass = 0.0
        for th in theta_ring
            x_ring = xc + r_ring * cos(th)
            y_ring = yc + r_ring * sin(th)
            ring_3d_mass +=
                rho_ref *
                marker_area(coords) *
                marker_out_of_plane_length(x_ring, y_ring, xc, yc)
            ring_2d_mass += rho_ref * marker_area(coords)
        end
        @test isapprox(ring_3d_mass, ring_2d_mass * R_planet; rtol=1e-12)

        # 4. Physical edge cases and non-finite input validation
        # Marker at the center has out-of-plane length 0
        @test isapprox(marker_out_of_plane_length(xc, yc, xc, yc), 0.0; atol=1e-12)

        # Non-finite coordinates throw DomainError
        @test_throws DomainError marker_out_of_plane_length(NaN, yc, xc, yc)
        @test_throws DomainError marker_out_of_plane_length(xc, Inf, xc, yc)
        @test_throws DomainError marker_out_of_plane_length(xc, yc, -Inf, yc)
        @test_throws DomainError marker_out_of_plane_length(xc, yc, xc, NaN)

        # Non-square grid marker area: dx*dy/(Nxmc*Nymc)
        coords_nonsquare = GridCoordinates(
            33, 65; xsize=100_000.0, ysize=200_000.0, Nxmc=3, Nymc=5
        )
        expected_nonsquare_area =
            (coords_nonsquare.dx * coords_nonsquare.dy) /
            (coords_nonsquare.Nxmc * coords_nonsquare.Nymc)
        @test isapprox(marker_area(coords_nonsquare), expected_nonsquare_area; rtol=1e-12)
    end
end

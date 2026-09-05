@testset "Particles" begin
    @testset "setup_interpolated_properties(): dimensions and zero-initialization" begin
        props = Erebus.setup_interpolated_properties()
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
            DMPSUM,
            DHPSUM,
            XWSSUM,
            WTPSUM,
        ) = props

        # 1. Basic grid accumulation tensors: (Ny, Nx)
        for arr in [ETA0SUM, ETASUM, GGGSUM, SXYSUM, COHSUM, TENSUM, FRISUM, WTSUM]
            @test size(arr) == (Ny, Nx)
            @test all(iszero, arr)
        end

        # 2. Staggered grid accumulation tensors: (Ny1, Nx1)
        staggered_accum = [
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
            DMPSUM,
            DHPSUM,
            XWSSUM,
            WTPSUM,
        ]
        for arr in staggered_accum
            @test size(arr) == (Ny1, Nx1)
            @test all(iszero, arr)
        end
    end # testset "setup_interpolated_properties()"

    @testset "reset_interpolated_properties!(): zeroing guarantees" begin
        props = Erebus.setup_interpolated_properties()
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
            DMPSUM,
            DHPSUM,
            XWSSUM,
            WTPSUM,
        ) = props

        # Fill with non-zero values
        for arr in [
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
        ]
            fill!(arr, 42.0)
        end
        @test any(!iszero, ETA0SUM)
        @test any(!iszero, WTPSUM)

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

        for arr in [
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
        ]
            @test all(iszero, arr)
        end
    end # testset "reset_interpolated_properties!()"

    @testset "reset_thermochemical_properties!(): thermochemical zeroing" begin
        DMPSUM = fill(1.0e-5, Ny1, Nx1)
        DHPSUM = fill(2.0e-3, Ny1, Nx1)
        WTPSUM = fill(4.0, Ny1, Nx1)

        Erebus.reset_thermochemical_properties!(DMPSUM, DHPSUM, WTPSUM)
        @test all(iszero, DMPSUM)
        @test all(iszero, DHPSUM)
        @test all(iszero, WTPSUM)
    end # testset "reset_thermochemical_properties!()"

    @testset "setup_marker_properties() & helpers: vector allocation and types" begin
        n_markers = 500
        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0) = Erebus.setup_marker_properties(
            n_markers
        )

        for v in [xm, ym, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0]
            @test length(v) == n_markers
            @test eltype(v) == Float64
            @test all(iszero, v)
        end
        @test length(tm) == n_markers
        @test eltype(tm) == Int
        @test all(iszero, tm)

        helpers = Erebus.setup_marker_properties_helpers(n_markers)
        for h in helpers
            @test length(h) == n_markers
            @test eltype(h) == Float64
            @test all(iszero, h)
        end

        mdis, mnum = Erebus.setup_marker_geometry_helpers()
        @test size(mdis) == (Nym, Nxm)
        @test size(mnum) == (Nym, Nxm)
        @test eltype(mdis) == Float64
        @test eltype(mnum) == Int
    end # testset "setup_marker_properties()"

    @testset "define_markers!() & compute_marker_properties!(): planetary zoning and bounds" begin
        marknum = start_marknum
        (xm, ym, tm, tkm, sxxm, sxym, etavpm, phim, phinewm, pfm0, XWsolidm, XWsolidm0) = Erebus.setup_marker_properties(
            marknum
        )
        (rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm, tkm_rhocptotalm, etafluidcur_inv_kphim, inv_gggtotalm, fricttotalm, cohestotalm, tenstotalm, rhofluidcur, alphasolidcur, alphafluidcur) = Erebus.setup_marker_properties_helpers(
            marknum
        )

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
            XWsolidm0;
            randomized=false,
            rcrust_val=40_000.0,
            rplanet_val=50_000.0,
        )

        # 1. Planetary spatial zoning invariant: 3-zone core, crust, sticky-air discrimination
        @test count(==(1), tm) > 0
        @test count(==(2), tm) > 0
        @test count(==(3), tm) > 0

        cx, cy = xcenter, ycenter
        r_crust_test = 40_000.0
        r_planet_test = 50_000.0
        for m in 1:marknum
            r = Erebus.distance(xm[m], ym[m], cx, cy)
            if r < r_planet_test
                if r > r_crust_test
                    @test tm[m] == 2  # Crust
                else
                    @test tm[m] == 1  # Core / mantle
                end
            else
                @test tm[m] == 3  # Sticky air / space
            end
        end

        # 2. Marker property computation
        hrsolidm, hrfluidm = start_hrsolidm, start_hrfluidm
        for m in 1:marknum
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
                marker_property_mode,
            )
        end

        # 3. Physical property bounds across all markers
        for m in 1:marknum
            # Density positivity
            @test rhototalm[m] > 0.0
            @test 1.0 <= rhototalm[m] <= 5000.0

            # Heat capacity positivity
            @test rhocptotalm[m] > 0.0

            # Thermal conductivity positivity
            @test ktotalm[m] > 0.0

            # Viscosity floor
            @test etatotalm[m] >= etamin

            # Porosity bounds
            @test phimin <= phim[m] <= phimax

            # Sticky air has zero radiogenic heating
            if tm[m] == 3
                @test isapprox(hrtotalm[m], 0.0; atol=1e-12)
            else
                @test hrtotalm[m] >= 0.0
            end
        end

        # 4. Thermal buoyancy and fluid density integration
        rhof_test = zeros(3)
        tm_test = [1, 2, 3]
        tk_test = [200.0, 373.0, 300.0]
        rhotot_test = zeros(3)
        rhocptot_test = zeros(3)
        etatot_test = zeros(3)
        hrtot_test = zeros(3)
        ktot_test = zeros(3)
        tkm_rhocptot_test = zeros(3)
        eta_inv_k_test = zeros(3)
        phi_test = [0.1, 0.1, phimin]
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
                marker_property_mode,
                rhof_test;
                thermal_buoyancy=true,
                alphafluid=2.0e-4,
                tmfluidphase_val=273.0,
            )
        end

        # Sub-freezing marker (T = 200 K) -> frozen ice density
        @test isapprox(rhof_test[1], 917.0; rtol=1e-12)
        # Hot rock marker (T = 373 K) -> thermal expansion (density drops below 1000)
        @test rhof_test[2] < 1000.0
        @test isapprox(rhof_test[2], 1000.0 * (1.0 - 2.0e-4 * 100.0); rtol=1e-12)
        # Sticky air marker -> density is 1.0
        @test isapprox(rhof_test[3], 1.0; rtol=1e-12)
    end # testset "define_markers!() & compute_marker_properties!()"

    @testset "fix(), fix_distances(), and fix_weights(): bilinear interpolation axioms" begin
        # 1. Coordinate clamping outside domain
        # Far left and above domain: must clamp to (imin, jmin)
        i_cl, j_cl = Erebus.fix(
            -10.0 * dx,
            -10.0 * dy,
            x,
            y,
            dx,
            dy,
            jmin_basic,
            jmax_basic,
            imin_basic,
            imax_basic,
        )
        @test j_cl == jmin_basic
        @test i_cl == imin_basic

        # Far right and below domain: must clamp to (imax, jmax) across all staggered grids
        # Basic grid clamp
        i_cr, j_cr = Erebus.fix(
            xsize + 10.0 * dx,
            ysize + 10.0 * dy,
            x,
            y,
            dx,
            dy,
            jmin_basic,
            jmax_basic,
            imin_basic,
            imax_basic,
        )
        @test j_cr == jmax_basic
        @test i_cr == imax_basic

        # Vx grid clamp
        i_vx, j_vx = Erebus.fix(
            xsize + 10.0 * dx,
            ysize + 10.0 * dy,
            xvx,
            yvx,
            dx,
            dy,
            jmin_vx,
            jmax_vx,
            imin_vx,
            imax_vx,
        )
        @test j_vx == jmax_vx
        @test i_vx == imax_vx

        # Vy grid clamp
        i_vy, j_vy = Erebus.fix(
            xsize + 10.0 * dx,
            ysize + 10.0 * dy,
            xvy,
            yvy,
            dx,
            dy,
            jmin_vy,
            jmax_vy,
            imin_vy,
            imax_vy,
        )
        @test j_vy == jmax_vy
        @test i_vy == imax_vy

        # P grid clamp
        i_p, j_p = Erebus.fix(
            xsize + 10.0 * dx,
            ysize + 10.0 * dy,
            xp,
            yp,
            dx,
            dy,
            jmin_p,
            jmax_p,
            imin_p,
            imax_p,
        )
        @test j_p == jmax_p
        @test i_p == imax_p

        # 2. fix_distances(): bilinear distance decomposition axioms
        cx_dist, cy_dist = 0.35 * xsize, 0.65 * ysize
        i_d, j_d, dxmj, dymi = Erebus.fix_distances(
            cx_dist, cy_dist, x, y, dx, dy, jmin_basic, jmax_basic, imin_basic, imax_basic
        )
        @test isapprox(dxmj, cx_dist - x[j_d]; rtol=1e-12)
        @test isapprox(dymi, cy_dist - y[i_d]; rtol=1e-12)
        @test 0.0 <= dxmj <= dx
        @test 0.0 <= dymi <= dy

        # 3. Partition of unity: sum of weights equals 1.0 across all grid types
        coords_test = [
            (0.23 * xsize, 0.45 * ysize),
            (0.78 * xsize, 0.12 * ysize),
            (0.50 * xsize, 0.50 * ysize),
        ]
        for (cx, cy) in coords_test
            # Basic grid
            _, _, w_b = Erebus.fix_weights(
                cx, cy, x, y, dx, dy, jmin_basic, jmax_basic, imin_basic, imax_basic
            )
            @test isapprox(sum(w_b), 1.0; rtol=1e-12)
            @test all(0.0 .<= w_b .<= 1.0)

            # Pressure grid
            _, _, w_p = Erebus.fix_weights(
                cx, cy, xp, yp, dx, dy, jmin_p, jmax_p, imin_p, imax_p
            )
            @test isapprox(sum(w_p), 1.0; rtol=1e-12)
            @test all(0.0 .<= w_p .<= 1.0)

            # Vx grid
            _, _, w_vx = Erebus.fix_weights(
                cx, cy, xvx, yvx, dx, dy, jmin_vx, jmax_vx, imin_vx, imax_vx
            )
            @test isapprox(sum(w_vx), 1.0; rtol=1e-12)
            @test all(0.0 .<= w_vx .<= 1.0)

            # Vy grid
            _, _, w_vy = Erebus.fix_weights(
                cx, cy, xvy, yvy, dx, dy, jmin_vy, jmax_vy, imin_vy, imax_vy
            )
            @test isapprox(sum(w_vy), 1.0; rtol=1e-12)
            @test all(0.0 .<= w_vy .<= 1.0)
        end

        # 3. Exact nodal reproduction: marker placed exactly at node (x[j], y[i])
        i_node, j_node = 3, 4
        _, _, w_exact = Erebus.fix_weights(
            x[j_node],
            y[i_node],
            x,
            y,
            dx,
            dy,
            jmin_basic,
            jmax_basic,
            imin_basic,
            imax_basic,
        )
        @test isapprox(w_exact[1], 1.0; atol=1e-12)  # Top-left node weight is 1.0
        @test isapprox(w_exact[2], 0.0; atol=1e-12)
        @test isapprox(w_exact[3], 0.0; atol=1e-12)
        @test isapprox(w_exact[4], 0.0; atol=1e-12)

        # 4. Cell center symmetry: marker at center of cell has equal 0.25 weights
        cx_mid = x[j_node] + 0.5 * dx
        cy_mid = y[i_node] + 0.5 * dy
        _, _, w_mid = Erebus.fix_weights(
            cx_mid, cy_mid, x, y, dx, dy, jmin_basic, jmax_basic, imin_basic, imax_basic
        )
        for k in 1:4
            @test isapprox(w_mid[k], 0.25; rtol=1e-12)
        end

        # 5. Cell edge linearity: marker on horizontal edge has 2 non-zero weights summing to 1
        cx_edge = x[j_node] + 0.3 * dx
        cy_edge = y[i_node]
        _, _, w_edge = Erebus.fix_weights(
            cx_edge, cy_edge, x, y, dx, dy, jmin_basic, jmax_basic, imin_basic, imax_basic
        )
        @test isapprox(w_edge[1], 0.7; rtol=1e-12)
        @test isapprox(w_edge[3], 0.3; rtol=1e-12)
        @test isapprox(w_edge[2], 0.0; atol=1e-12)
        @test isapprox(w_edge[4], 0.0; atol=1e-12)
    end # testset "fix(), fix_distances(), and fix_weights()"

    @testset "reduce_add_3darray!() & interpolate_add_to_grid!(): conservation" begin
        # 1. 3D array reduction sum conservation
        ny_t, nx_t, nth = 5, 5, 4
        src_3d = rand(rgen, ny_t, nx_t, nth)
        expected_sum = sum(src_3d)
        Erebus.reduce_add_3darray!(src_3d)
        reduced = src_3d[:, :, 1]
        @test isapprox(sum(reduced), expected_sum; rtol=1e-12)
        @test size(reduced) == (ny_t, nx_t)

        # 2. Interpolate add to grid conservation: partition of unity ensures total added equals input
        grid_target = zeros(6, 6)
        weights = @SVector [0.2, 0.3, 0.1, 0.4]  # sum == 1.0
        val_add = 125.0
        Erebus.interpolate_add_to_grid!(2, 2, weights, val_add, grid_target)
        @test isapprox(sum(grid_target), val_add; rtol=1e-12)
        @test grid_target[2, 2] ≈ 0.2 * val_add
        @test grid_target[3, 2] ≈ 0.3 * val_add
        @test grid_target[2, 3] ≈ 0.1 * val_add
        @test grid_target[3, 3] ≈ 0.4 * val_add
    end # testset "reduce_add_3darray!()"

    @testset "interpolate_to_marker!(): convex hull and constant reproduction" begin
        grid_vals = rand(rgen, 6, 6)
        weights = @SVector [0.15, 0.25, 0.35, 0.25]
        marker_prop = zeros(1)

        Erebus.interpolate_to_marker!(1, 2, 3, weights, marker_prop, grid_vals)

        # Convex hull: interpolated value must lie within extrema of the 4 nodes
        cell_vals = [grid_vals[2, 3], grid_vals[3, 3], grid_vals[2, 4], grid_vals[3, 4]]
        @test minimum(cell_vals) <= marker_prop[1] <= maximum(cell_vals)

        # Constant field reproduction
        c_val = 137.0
        grid_const = fill(c_val, 6, 6)
        Erebus.interpolate_to_marker!(1, 2, 3, weights, marker_prop, grid_const)
        @test isapprox(marker_prop[1], c_val; rtol=1e-12)

        # interpolate_add_to_marker! accumulates correctly
        Erebus.interpolate_add_to_marker!(1, 2, 3, weights, marker_prop, grid_const)
        @test isapprox(marker_prop[1], 2.0 * c_val; rtol=1e-12)
    end # testset "interpolate_to_marker!()"

    @testset "marker_to_*_nodes!(): weight and property accumulation" begin
        marknum = 100
        xm = rand(rgen, (0.1 * xsize):0.01:(0.9 * xsize), marknum)
        ym = rand(rgen, (0.1 * ysize):0.01:(0.9 * ysize), marknum)
        tm = fill(1, marknum)  # All rocks

        # Basic nodes accumulation
        ETA0SUM = zeros(Ny, Nx)
        ETASUM = zeros(Ny, Nx)
        GGGSUM = zeros(Ny, Nx)
        SXYSUM = zeros(Ny, Nx)
        COHSUM = zeros(Ny, Nx)
        TENSUM = zeros(Ny, Nx)
        FRISUM = zeros(Ny, Nx)
        WTSUM = zeros(Ny, Nx)
        etatotalm = fill(1.0e20, marknum)
        etavpm = fill(1.0e20, marknum)
        inv_gggtotalm = fill(1.0e-10, marknum)
        sxym = zeros(marknum)
        cohestotalm = fill(1.0e7, marknum)
        tenstotalm = fill(1.0e7, marknum)
        fricttotalm = fill(0.6, marknum)

        for m in 1:marknum
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

        # Total weights accumulated must equal number of interior markers
        @test isapprox(sum(WTSUM), Float64(marknum); rtol=1e-10)
        # Viscosity average on nodes where markers exist must equal constant value
        for j in 1:Nx, i in 1:Ny
            if WTSUM[i, j] > 0.0
                @test isapprox(ETASUM[i, j] / WTSUM[i, j], 1.0e20; rtol=1e-12)
            end
        end

        # P nodes accumulation
        GGGPSUM = zeros(Ny1, Nx1)
        SXXSUM = zeros(Ny1, Nx1)
        RHOSUM = zeros(Ny1, Nx1)
        RHOCPSUM = zeros(Ny1, Nx1)
        ALPHASUM = zeros(Ny1, Nx1)
        ALPHAFSUM = zeros(Ny1, Nx1)
        HRSUM = zeros(Ny1, Nx1)
        PHISUM = zeros(Ny1, Nx1)
        TKSUM = zeros(Ny1, Nx1)
        WTPSUM = zeros(Ny1, Nx1)
        rhototalm = fill(3200.0, marknum)
        rhocptotalm = fill(4.0e6, marknum)
        alphasolidcur = fill(3.0e-5, marknum)
        alphafluidcur = fill(2.0e-4, marknum)
        hrtotalm = fill(1.0e-7, marknum)
        phim = fill(0.1, marknum)
        tkm_rhocptotalm = fill(300.0 * 4.0e6, marknum)
        sxxm = zeros(marknum)

        for m in 1:marknum
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
        @test isapprox(sum(WTPSUM), Float64(marknum); rtol=1e-10)
        for j in 1:Nx1, i in 1:Ny1
            if WTPSUM[i, j] > 0.0
                @test isapprox(RHOSUM[i, j] / WTPSUM[i, j], 3200.0; rtol=1e-12)
            end
        end
    end # testset "marker_to_*_nodes!()"

    @testset "compute_velocities!(): staggered averaging and boundary conditions" begin
        vx = rand(rgen, Ny1, Nx1)
        vy = rand(rgen, Ny1, Nx1)
        vxf = rand(rgen, Ny1, Nx1)
        vyf = rand(rgen, Ny1, Nx1)
        vxp = zeros(Ny1, Nx1)
        vyp = zeros(Ny1, Nx1)
        vxpf = zeros(Ny1, Nx1)
        vypf = zeros(Ny1, Nx1)

        Erebus.compute_velocities!(vx, vy, vxf, vyf, vxp, vyp, vxpf, vypf)

        # 1. Interior midpoint averaging: vxp[i, j] = 0.5 * (vx[i, j] + vx[i, j-1])
        for j in 2:Nx, i in 2:Ny
            @test isapprox(vxp[i, j], 0.5 * (vx[i, j] + vx[i, j - 1]); rtol=1e-12)
            @test isapprox(vyp[i, j], 0.5 * (vy[i, j] + vy[i - 1, j]); rtol=1e-12)
            @test isapprox(vxpf[i, j], 0.5 * (vxf[i, j] + vxf[i, j - 1]); rtol=1e-12)
            @test isapprox(vypf[i, j], 0.5 * (vyf[i, j] + vyf[i - 1, j]); rtol=1e-12)
        end

        # 2. Free-slip boundary condition relations:
        # Top boundary vxp[1, :] = -bctop * vxp[2, :]
        for j in 2:(Nx - 1)
            @test isapprox(vxp[1, j], -bctop * vxp[2, j]; rtol=1e-12)
            @test isapprox(vxpf[1, j], -bcftop * vxpf[2, j]; rtol=1e-12)
        end
    end # testset "compute_velocities!()"

    @testset "compute_rotation_rate!(): rigid body and irrotational flow" begin
        wyx = zeros(Ny, Nx)

        # 1. Rigid body rotation: vx = -Omega * y, vy = Omega * x
        # Rotation rate wyx = 0.5 * (dvy/dx - dvx/dy) = 0.5 * (Omega - (-Omega)) = Omega
        Omega = 2.5e-5
        vx_rot = zeros(Ny1, Nx1)
        vy_rot = zeros(Ny1, Nx1)
        for j in 1:Nx1, i in 1:Ny1
            vx_rot[i, j] = -Omega * yvx[i]
            vy_rot[i, j] = Omega * xvy[j]
        end

        Erebus.compute_rotation_rate!(vx_rot, vy_rot, wyx)
        for j in 1:Nx, i in 1:Ny
            @test isapprox(wyx[i, j], Omega; rtol=1e-10)
        end

        # 2. Irrotational flow: vx = c * x, vy = -c * y (pure shear / extension)
        # dvy/dx = 0, dvx/dy = 0 -> wyx = 0 identically
        vx_irrot = zeros(Ny1, Nx1)
        vy_irrot = zeros(Ny1, Nx1)
        rate = 1.0e-6
        for j in 1:Nx1, i in 1:Ny1
            vx_irrot[i, j] = rate * xvx[j]
            vy_irrot[i, j] = -rate * yvy[i]
        end

        Erebus.compute_rotation_rate!(vx_irrot, vy_irrot, wyx)
        for j in 1:Nx, i in 1:Ny
            @test isapprox(wyx[i, j], 0.0; atol=1e-12)
        end
    end # testset "compute_rotation_rate!()"

    @testset "move_markers_rk4!(): advective translation and zero-velocity invariance" begin
        marknum = 20
        xm_init = fill(xsize / 2.0, marknum)
        ym_init = fill(ysize / 2.0, marknum)
        xm = copy(xm_init)
        ym = copy(ym_init)
        tm = fill(1, marknum)
        tkm = fill(300.0, marknum)
        phim = fill(0.1, marknum)
        sxym = zeros(marknum)
        sxxm = zeros(marknum)
        tk2 = fill(300.0, Ny1, Nx1)
        wyx = zeros(Ny, Nx)
        dt_adv = 100.0

        # 1. Zero velocity field: markers must remain stationary
        Erebus.move_markers_rk4!(
            xm,
            ym,
            tm,
            tkm,
            phim,
            sxym,
            sxxm,
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            wyx,
            tk2,
            marknum,
            dt_adv,
            1,
        )
        @test isapprox(xm, xm_init; rtol=1e-12)
        @test isapprox(ym, ym_init; rtol=1e-12)

        # 2. Uniform velocity field: exact displacement dx = U * dt, dy = V * dt
        U_val = 1.0e-6
        V_val = -5.0e-7
        vx_uniform = fill(U_val, Ny1, Nx1)
        vy_uniform = fill(V_val, Ny1, Nx1)

        Erebus.move_markers_rk4!(
            xm,
            ym,
            tm,
            tkm,
            phim,
            sxym,
            sxxm,
            vx_uniform,
            vy_uniform,
            vx_uniform,
            vy_uniform,
            wyx,
            tk2,
            marknum,
            dt_adv,
            1,
        )
        expected_xm = xm_init .+ U_val * dt_adv
        expected_ym = ym_init .+ V_val * dt_adv
        @test isapprox(xm, expected_xm; rtol=1e-10)
        @test isapprox(ym, expected_ym; rtol=1e-10)

        # 3. Analytical hyperbolic trajectory: exact closed-form solution discriminating RK4 from explicit Euler
        # For vx = γ * y and vy = γ * x with symmetric initial coordinates x0 = y0,
        # the exact continuous trajectory is x(t) = y(t) = x0 * exp(γ * t).
        # Multi-stage RK4 integrates this to 4th-order accuracy O((γ*dt)⁵),
        # while explicit Euler produces only x0 * (1 + γ * dt), missing higher-order terms by ~3.5 m.
        γ_rate = 1.0e-4
        dt_rk4 = 100.0
        vx_grad = [γ_rate * yvx[i] for i in 1:Ny1, j in 1:Nx1]
        vy_grad = [γ_rate * xvy[j] for i in 1:Ny1, j in 1:Nx1]

        x0_sym = 70_000.0
        y0_sym = 70_000.0
        xm_sym = [x0_sym]
        ym_sym = [y0_sym]
        tm_sym = [1]
        tkm_sym = [300.0]
        phim_sym = [0.1]
        sxym_sym = [0.0]
        sxxm_sym = [0.0]

        Erebus.move_markers_rk4!(
            xm_sym,
            ym_sym,
            tm_sym,
            tkm_sym,
            phim_sym,
            sxym_sym,
            sxxm_sym,
            vx_grad,
            vy_grad,
            vx_grad,
            vy_grad,
            wyx,
            tk2,
            1,
            dt_rk4,
            1,
        )

        dx_disp = xm_sym[1] - x0_sym
        dy_disp = ym_sym[1] - y0_sym
        exact_disp = x0_sym * (exp(γ_rate * dt_rk4) - 1.0)
        euler_disp = x0_sym * γ_rate * dt_rk4

        # RK4 matches exact analytical solution to within 1 cm (order of 58 nm error)
        @test isapprox(dx_disp, exact_disp; atol=1e-2)
        @test isapprox(dy_disp, exact_disp; atol=1e-2)

        # Discrimination guard: RK4 displacement strictly exceeds explicit Euler step by > 3.0 m
        @test abs(dx_disp - euler_disp) > 3.0
        @test abs(dy_disp - euler_disp) > 3.0
        @test dx_disp > euler_disp
        @test dy_disp > euler_disp
    end # testset "move_markers_rk4!()"

    @testset "backtrace_pressures_rk4!(): zero-velocity pressure invariance" begin
        pr = rand(rgen, Ny1, Nx1) .* 1.0e7
        ps = rand(rgen, Ny1, Nx1) .* 1.0e7
        pf = rand(rgen, Ny1, Nx1) .* 1.0e7
        pr0 = zeros(Ny1, Nx1)
        ps0 = zeros(Ny1, Nx1)
        pf0 = zeros(Ny1, Nx1)
        dt_val = 10.0

        # With zero velocity, backtracing pressure nodes should reproduce current pressure
        Erebus.backtrace_pressures_rk4!(
            pr,
            pr0,
            ps,
            ps0,
            pf,
            pf0,
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            dt_val,
        )
        for j in 2:Nx, i in 2:Ny
            @test isapprox(pr0[i, j], pr[i, j]; rtol=1e-10)
            @test isapprox(ps0[i, j], ps[i, j]; rtol=1e-10)
            @test isapprox(pf0[i, j], pf[i, j]; rtol=1e-10)
        end
    end # testset "backtrace_pressures_rk4!()"

    @testset "replenish_markers!(): population recovery and domain bounds" begin
        marknum = 100
        # Concentrate markers in one quadrant to simulate a depleted void elsewhere
        xm = rand(rgen, 0.0:0.1:(xsize / 2.0), marknum)
        ym = rand(rgen, 0.0:0.1:(ysize / 2.0), marknum)
        tm = fill(1, marknum)
        tkm = fill(300.0, marknum)
        sxxm = zeros(marknum)
        sxym = zeros(marknum)
        etavpm = fill(1.0e21, marknum)
        phim = fill(0.1, marknum)
        phinewm = copy(phim)
        pfm0 = fill(1.0e7, marknum)
        XWsolidm = fill(0.5, marknum)
        XWsolidm0 = copy(XWsolidm)
        (rhototalm, rhocptotalm, etatotalm, hrtotalm, ktotalm, tkm_rhocptotalm, etafluidcur_inv_kphim, inv_gggtotalm, fricttotalm, cohestotalm, tenstotalm, rhofluidcur, alphasolidcur, alphafluidcur) = Erebus.setup_marker_properties_helpers(
            marknum
        )
        mdis, mnum = Erebus.setup_marker_geometry_helpers()

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

        # 1. Population recovery: new markers added to depleted cells
        @test marknum_new >= marknum
        @test length(xm) == marknum_new
        @test length(ym) == marknum_new

        # 2. Newly added marker coordinates are strictly within physical domain
        for m in 1:marknum_new
            @test -dx <= xm[m] <= xsize + dx
            @test -dy <= ym[m] <= ysize + dy
            @test tm[m] in (1, 2, 3)
            @test isfinite(tkm[m])
            @test isfinite(phim[m])
        end
    end # testset "replenish_markers!()"

    @testset "update_marker_porosity!(): compaction, dilation, and clamp bounds" begin
        xm_p = [0.5 * xsize, 0.5 * xsize, 0.5 * xsize, 0.5 * xsize]
        ym_p = [0.5 * ysize, 0.5 * ysize, 0.5 * ysize, 0.5 * ysize]
        tm_p = [1, 1, 1, 3]
        phi_0 = 0.20
        phim_p = [phi_0, phi_0, phi_0, phi_0]
        dt_p = 1.0e4

        # 1. Compaction: positive volumetric strain rate APHI > 0 strictly decreases porosity
        APHI_compact = fill(1.0e-5, Ny1, Nx1)
        phim_compact = copy(phim_p)
        Erebus.update_marker_porosity!(
            xm_p, ym_p, tm_p, phim_compact, APHI_compact, dt_p, 4
        )
        @test phim_compact[1] < phi_0
        expected_compact = phi_0 / ((1.0 - phi_0) * exp(1.0e-5 * dt_p) + phi_0)
        @test isapprox(phim_compact[1], expected_compact; rtol=1e-12)
        @test phim_compact[4] == phi_0  # Sticky air marker remains unchanged

        # 2. Dilation: negative volumetric strain rate APHI < 0 strictly increases porosity
        APHI_dilate = fill(-1.0e-5, Ny1, Nx1)
        phim_dilate = copy(phim_p)
        Erebus.update_marker_porosity!(xm_p, ym_p, tm_p, phim_dilate, APHI_dilate, dt_p, 4)
        @test phim_dilate[1] > phi_0
        expected_dilate = phi_0 / ((1.0 - phi_0) * exp(-1.0e-5 * dt_p) + phi_0)
        @test isapprox(phim_dilate[1], expected_dilate; rtol=1e-12)

        # 3. Discrimination guard: sign flip test
        # If exp(-aphim*dt) was computed instead of exp(aphim*dt), compaction and dilation would invert
        @test phim_compact[1] < phi_0 < phim_dilate[1]
        @test abs(phim_compact[1] - phim_dilate[1]) > 1.0e-3

        # 4. Extreme clamping bounds: phimin and phimax
        APHI_extreme_compact = fill(1.0, Ny1, Nx1)
        phim_min_test = [0.20]
        Erebus.update_marker_porosity!(
            [0.5 * xsize],
            [0.5 * ysize],
            [1],
            phim_min_test,
            APHI_extreme_compact,
            1.0e3,
            1;
            phimin=0.01,
            phimax=0.50,
        )
        @test phim_min_test[1] ≈ 0.01

        APHI_extreme_dilate = fill(-1.0, Ny1, Nx1)
        phim_max_test = [0.20]
        Erebus.update_marker_porosity!(
            [0.5 * xsize],
            [0.5 * ysize],
            [1],
            phim_max_test,
            APHI_extreme_dilate,
            1.0e3,
            1;
            phimin=0.01,
            phimax=0.50,
        )
        @test phim_max_test[1] ≈ 0.50
    end # testset "update_marker_porosity!()"

    @testset "update_marker_temperature!(): timestep initialization and incremental heating" begin
        xm_t = [0.5 * xsize]
        ym_t = [0.5 * ysize]
        tkm_t = [0.0]
        tk2_t = fill(350.0, Ny1, Nx1)
        DT_t = fill(5.0, Ny1, Nx1)

        # 1. Timestep 1: initializes tkm directly from tk2 grid
        Erebus.update_marker_temperature!(xm_t, ym_t, tkm_t, DT_t, tk2_t, 1, 1)
        @test isapprox(tkm_t[1], 350.0; rtol=1e-12)

        # 2. Timestep > 1: increments tkm by interpolated DT
        Erebus.update_marker_temperature!(xm_t, ym_t, tkm_t, DT_t, tk2_t, 2, 1)
        @test isapprox(tkm_t[1], 355.0; rtol=1e-12)

        # Repeated step
        Erebus.update_marker_temperature!(xm_t, ym_t, tkm_t, DT_t, tk2_t, 3, 1)
        @test isapprox(tkm_t[1], 360.0; rtol=1e-12)
    end # testset "update_marker_temperature!()"
end

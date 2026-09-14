using Test
using Erebus
using Random
import Erebus:
    setup_interpolated_properties,
    allocate_thread_interpolation_buffers,
    reduce_thread_buffers!,
    marker_to_basic_nodes!,
    marker_to_vx_nodes!,
    marker_to_vy_nodes!,
    marker_to_p_nodes!

@testset "Cell-Tiled P2M Scatter" begin
    @testset "Workspace Allocation & Sizing" begin
        coords = GridCoordinates(33, 33; xsize=140000.0, ysize=140000.0)
        marknum = 1000
        ws = P2MTiledWorkspace(coords, marknum, 4)

        @test ws.tile_size == 4
        @test ws.Ncx == 32
        @test ws.Ncy == 32
        @test ws.Ntx == 8
        @test ws.Nty == 8
        @test ws.Ntiles == 64
        @test length(ws.tile_markers) == marknum
        @test length(ws.tiles_by_color) == 4
        total_colored_tiles = sum(length(tiles) for tiles in ws.tiles_by_color)
        @test total_colored_tiles == 64

        # Verify tile_size guard (< 2 rejected)
        @test_throws ArgumentError P2MTiledWorkspace(coords, marknum, 1)
    end

    @testset "Marker Binning & Spatial Partitioning" begin
        coords = GridCoordinates(33, 33; xsize=140000.0, ysize=140000.0)
        marknum = 2500
        rng = MersenneTwister(12345)
        xm = coords.x[1] .+ rand(rng, marknum) .* coords.xsize
        ym = coords.y[1] .+ rand(rng, marknum) .* coords.ysize

        ws = P2MTiledWorkspace(coords, marknum, 4)
        bin_markers_into_tiles!(ws, xm, ym, coords, marknum)

        @test ws.tile_offsets[1] == 1
        @test ws.tile_offsets[ws.Ntiles + 1] == marknum + 1

        # Verify all markers appear once
        sorted_markers = sort(ws.tile_markers[1:marknum])
        @test sorted_markers == collect(1:marknum)

        # Verify each marker's cell coordinates match its tile
        all_markers_in_tile_bounds = true
        for t in 1:ws.Ntiles
            lo = ws.tile_offsets[t]
            hi = ws.tile_offsets[t + 1] - 1
            for idx in lo:hi
                m = ws.tile_markers[idx]
                if ws.marker_tiles[m] != t
                    all_markers_in_tile_bounds = false
                end
            end
        end
        @test all_markers_in_tile_bounds
    end

    @testset "Workspace Resizing & Compatibility" begin
        coords1 = GridCoordinates(33, 33; xsize=140000.0, ysize=140000.0)
        ws = P2MTiledWorkspace(coords1, 1000, 4)

        # Resizing when marker count grows
        ws_expanded = ensure_workspace_compatible(ws, coords1, 2500, 4)
        @test length(ws_expanded.tile_markers) >= 2500
        @test ws_expanded.Ntiles == 64

        # Reallocation when grid size changes
        coords2 = GridCoordinates(65, 65; xsize=280000.0, ysize=280000.0)
        ws_newgrid = ensure_workspace_compatible(ws_expanded, coords2, 2500, 4)
        @test ws_newgrid.Ncx == 64
        @test ws_newgrid.Ncy == 64
        @test ws_newgrid.Ntx == 16
        @test ws_newgrid.Nty == 16
        @test ws_newgrid.Ntiles == 256
    end

    @testset "Equivalence between Tiled and Buffered Scatter" begin
        coords = GridCoordinates(17, 17; xsize=140000.0, ysize=140000.0)
        marknum = 1200
        rng = MersenneTwister(999)
        xm = coords.x[1] .+ rand(rng, marknum) .* coords.xsize
        ym = coords.y[1] .+ rand(rng, marknum) .* coords.ysize

        # Synthetic properties for testing
        etatotalm = 1.0e18 .+ rand(rng, marknum) .* 1.0e19
        etavpm = 1.0e18 .+ rand(rng, marknum) .* 1.0e19
        inv_gggtotalm = fill(1.0e-10, marknum)
        sxym = rand(rng, marknum) .* 1.0e6
        cohestotalm = fill(1.0e8, marknum)
        tenstotalm = fill(6.0e7, marknum)
        fricttotalm = fill(0.6, marknum)
        rhototalm = 3300.0 .+ rand(rng, marknum) .* 100.0
        rhofluidcur = fill(1000.0, marknum)
        ktotalm = fill(3.0, marknum)
        phim = 0.05 .+ rand(rng, marknum) .* 0.15
        etafluidcur_inv_kphim = fill(1.0e10, marknum)
        sxxm = rand(rng, marknum) .* 1.0e6
        rhocptotalm = fill(3.3e6, marknum)
        alphasolidcur = fill(3.0e-5, marknum)
        alphafluidcur = fill(5.0e-5, marknum)
        hrtotalm = rand(rng, marknum) .* 1.0e-6
        tkm_rhocptotalm = fill(5.0e8, marknum)

        # Setup master grids for buffered execution
        (ETA0_buf, ETA_buf, GGG_buf, SXY_buf, COH_buf, TEN_buf, FRI_buf, WT_buf, RHOX_buf, RHOFX_buf, KX_buf, PHIX_buf, RX_buf, WTX_buf, RHOY_buf, RHOFY_buf, KY_buf, PHIY_buf, RY_buf, WTY_buf, RHO_buf, RHOCP_buf, ALPHA_buf, ALPHAF_buf, HR_buf, GGGP_buf, SXX_buf, TK_buf, PHI_buf, DMP_buf, DHP_buf, XWS_buf, WTP_buf) = setup_interpolated_properties(
            coords
        )

        # Setup master grids for tiled execution
        (ETA0_til, ETA_til, GGG_til, SXY_til, COH_til, TEN_til, FRI_til, WT_til, RHOX_til, RHOFX_til, KX_til, PHIX_til, RX_til, WTX_til, RHOY_til, RHOFY_til, KY_til, PHIY_til, RY_til, WTY_til, RHO_til, RHOCP_til, ALPHA_til, ALPHAF_til, HR_til, GGGP_til, SXX_til, TK_til, PHI_til, DMP_til, DHP_til, XWS_til, WTP_til) = setup_interpolated_properties(
            coords
        )

        # 1. Buffered execution using thread buffers
        buf = allocate_thread_interpolation_buffers(1, coords)[1]
        for m in 1:marknum
            marker_to_basic_nodes!(
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
                buf.ETA0SUM,
                buf.ETASUM,
                buf.GGGSUM,
                buf.SXYSUM,
                buf.COHSUM,
                buf.TENSUM,
                buf.FRISUM,
                buf.WTSUM;
                coords=coords,
            )
            marker_to_vx_nodes!(
                m,
                xm[m],
                ym[m],
                rhototalm,
                rhofluidcur,
                ktotalm,
                phim,
                etafluidcur_inv_kphim,
                buf.RHOXSUM,
                buf.RHOFXSUM,
                buf.KXSUM,
                buf.PHIXSUM,
                buf.RXSUM,
                buf.WTXSUM;
                coords=coords,
            )
            marker_to_vy_nodes!(
                m,
                xm[m],
                ym[m],
                rhototalm,
                rhofluidcur,
                ktotalm,
                phim,
                etafluidcur_inv_kphim,
                buf.RHOYSUM,
                buf.RHOFYSUM,
                buf.KYSUM,
                buf.PHIYSUM,
                buf.RYSUM,
                buf.WTYSUM;
                coords=coords,
            )
            marker_to_p_nodes!(
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
                buf.GGGPSUM,
                buf.SXXSUM,
                buf.RHOSUM,
                buf.RHOCPSUM,
                buf.ALPHASUM,
                buf.ALPHAFSUM,
                buf.HRSUM,
                buf.PHISUM,
                buf.TKSUM,
                buf.WTPSUM;
                coords=coords,
            )
        end
        reduce_thread_buffers!(
            ETA0_buf,
            ETA_buf,
            GGG_buf,
            SXY_buf,
            COH_buf,
            TEN_buf,
            FRI_buf,
            WT_buf,
            RHOX_buf,
            RHOFX_buf,
            KX_buf,
            PHIX_buf,
            RX_buf,
            WTX_buf,
            RHOY_buf,
            RHOFY_buf,
            KY_buf,
            PHIY_buf,
            RY_buf,
            WTY_buf,
            RHO_buf,
            RHOCP_buf,
            ALPHA_buf,
            ALPHAF_buf,
            HR_buf,
            GGGP_buf,
            SXX_buf,
            TK_buf,
            PHI_buf,
            WTP_buf,
            [buf],
        )

        # 2. Tiled execution
        ws = P2MTiledWorkspace(coords, marknum, 4)
        bin_markers_into_tiles!(ws, xm, ym, coords, marknum)
        for color in 1:4
            tiles = ws.tiles_by_color[color]
            Threads.@threads :dynamic for t in tiles
                lo = ws.tile_offsets[t]
                hi = ws.tile_offsets[t + 1] - 1
                lo > hi && continue
                for idx in lo:hi
                    m = ws.tile_markers[idx]
                    scatter_marker_to_master_grids!(
                        m,
                        xm[m],
                        ym[m],
                        coords,
                        etatotalm,
                        etavpm,
                        inv_gggtotalm,
                        sxym,
                        cohestotalm,
                        tenstotalm,
                        fricttotalm,
                        ETA0_til,
                        ETA_til,
                        GGG_til,
                        SXY_til,
                        COH_til,
                        TEN_til,
                        FRI_til,
                        WT_til,
                        rhototalm,
                        rhofluidcur,
                        ktotalm,
                        phim,
                        etafluidcur_inv_kphim,
                        RHOX_til,
                        RHOFX_til,
                        KX_til,
                        PHIX_til,
                        RX_til,
                        WTX_til,
                        RHOY_til,
                        RHOFY_til,
                        KY_til,
                        PHIY_til,
                        RY_til,
                        WTY_til,
                        sxxm,
                        rhocptotalm,
                        alphasolidcur,
                        alphafluidcur,
                        hrtotalm,
                        tkm_rhocptotalm,
                        GGGP_til,
                        SXX_til,
                        RHO_til,
                        RHOCP_til,
                        ALPHA_til,
                        ALPHAF_til,
                        HR_til,
                        PHI_til,
                        TK_til,
                        WTP_til,
                    )
                end
            end
        end

        # Assert equivalence across grids within floating-point tolerance
        @test isapprox(ETA0_buf, ETA0_til; rtol=1e-12, atol=1e-12)
        @test isapprox(ETA_buf, ETA_til; rtol=1e-12, atol=1e-12)
        @test isapprox(WT_buf, WT_til; rtol=1e-12, atol=1e-12)
        @test isapprox(RHOX_buf, RHOX_til; rtol=1e-12, atol=1e-12)
        @test isapprox(WTX_buf, WTX_til; rtol=1e-12, atol=1e-12)
        @test isapprox(RHOY_buf, RHOY_til; rtol=1e-12, atol=1e-12)
        @test isapprox(WTY_buf, WTY_til; rtol=1e-12, atol=1e-12)
        @test isapprox(RHO_buf, RHO_til; rtol=1e-12, atol=1e-12)
        @test isapprox(WTP_buf, WTP_til; rtol=1e-12, atol=1e-12)
    end

    @testset "Deterministic Concurrent Execution & Race-Freedom" begin
        coords = GridCoordinates(25, 25; xsize=100000.0, ysize=100000.0)
        marknum = 2000
        rng = MersenneTwister(42)
        xm = coords.x[1] .+ rand(rng, marknum) .* coords.xsize
        ym = coords.y[1] .+ rand(rng, marknum) .* coords.ysize
        rhototalm = 3200.0 .+ rand(rng, marknum) .* 200.0
        run_scatter = function ()
            (ETA0, ETA, GGG, SXY, COH, TEN, FRI, WT, RHOX, RHOFX, KX, PHIX, RX, WTX, RHOY, RHOFY, KY, PHIY, RY, WTY, RHO, RHOCP, ALPHA, ALPHAF, HR, GGGP, SXX, TK, PHI, DMP, DHP, XWS, WTP) = setup_interpolated_properties(
                coords
            )

            dummy_vec = zeros(Float64, marknum)
            ws = P2MTiledWorkspace(coords, marknum, 4)
            bin_markers_into_tiles!(ws, xm, ym, coords, marknum)
            for color in 1:4
                tiles = ws.tiles_by_color[color]
                Threads.@threads :dynamic for t in tiles
                    lo = ws.tile_offsets[t]
                    hi = ws.tile_offsets[t + 1] - 1
                    lo > hi && continue
                    for idx in lo:hi
                        m = ws.tile_markers[idx]
                        scatter_marker_to_master_grids!(
                            m,
                            xm[m],
                            ym[m],
                            coords,
                            dummy_vec,
                            dummy_vec,
                            dummy_vec,
                            dummy_vec,
                            dummy_vec,
                            dummy_vec,
                            dummy_vec,
                            ETA0,
                            ETA,
                            GGG,
                            SXY,
                            COH,
                            TEN,
                            FRI,
                            WT,
                            rhototalm,
                            dummy_vec,
                            dummy_vec,
                            dummy_vec,
                            dummy_vec,
                            RHOX,
                            RHOFX,
                            KX,
                            PHIX,
                            RX,
                            WTX,
                            RHOY,
                            RHOFY,
                            KY,
                            PHIY,
                            RY,
                            WTY,
                            dummy_vec,
                            dummy_vec,
                            dummy_vec,
                            dummy_vec,
                            dummy_vec,
                            dummy_vec,
                            GGGP,
                            SXX,
                            RHO,
                            RHOCP,
                            ALPHA,
                            ALPHAF,
                            HR,
                            PHI,
                            TK,
                            WTP,
                        )
                    end
                end
            end
            return RHO, WTP
        end

        # Run twice with multithreading
        RHO_1, WTP_1 = run_scatter()
        RHO_2, WTP_2 = run_scatter()

        # Deterministic bitwise equivalence across concurrent runs
        @test RHO_1 == RHO_2
        @test WTP_1 == WTP_2
        @test sum(WTP_1) > 0.0
    end
end

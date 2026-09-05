using Erebus
using Erebus.Particles
using Test
using Random
using StaticArrays

@testset "Multi-Threading & Adaptive Timestepping" begin
    @testset "ThreadInterpolationBuffers allocation and reset" begin
        coords = GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        nthreads = 4
        buffers = allocate_thread_interpolation_buffers(nthreads, coords)
        @test length(buffers) == 4
        for b in buffers
            @test size(b.ETA0SUM) == (17, 17)
            @test size(b.RHOXSUM) == (18, 18)
            @test size(b.RHOYSUM) == (18, 18)
            @test size(b.RHOSUM) == (18, 18)
            @test all(iszero, b.ETA0SUM)
        end

        # Dirty buffers and reset
        buffers[1].ETA0SUM .= 1.0
        buffers[2].RHOSUM .= 2.5
        reset_thread_buffers!(buffers)
        @test all(iszero, buffers[1].ETA0SUM)
        @test all(iszero, buffers[2].RHOSUM)
    end

    @testset "ThreadInterpolationBuffers reduction equivalence" begin
        coords = GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        buffers = allocate_thread_interpolation_buffers(2, coords)
        (ETA0SUM, ETASUM, GGGSUM, SXYSUM, COHSUM, TENSUM, FRISUM, WTSUM, RHOXSUM, RHOFXSUM, KXSUM, PHIXSUM, RXSUM, WTXSUM, RHOYSUM, RHOFYSUM, KYSUM, PHIYSUM, RYSUM, WTYSUM, RHOSUM, RHOCPSUM, ALPHASUM, ALPHAFSUM, HRSUM, GGGPSUM, SXXSUM, TKSUM, PHISUM, DMPSUM, DHPSUM, XWSSUM, WTPSUM) = setup_interpolated_properties(
            coords
        )

        # Populate buffer 1 and buffer 2
        buffers[1].ETA0SUM .= 1.5
        buffers[2].ETA0SUM .= 2.5
        buffers[1].RHOSUM .= 100.0
        buffers[2].RHOSUM .= 200.0

        reduce_thread_buffers!(
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
            buffers,
        )

        @test all(isapprox.(ETA0SUM, 4.0; rtol=1e-12))
        @test all(isapprox.(RHOSUM, 300.0; rtol=1e-12))
    end

    @testset "Marker interpolation: serial vs threaded buffer reduction" begin
        coords = GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        marknum = 200
        rng = MersenneTwister(42)

        # Setup serial target arrays
        (ETA0SUM_s, ETASUM_s, GGGSUM_s, SXYSUM_s, COHSUM_s, TENSUM_s, FRISUM_s, WTSUM_s, RHOXSUM_s, RHOFXSUM_s, KXSUM_s, PHIXSUM_s, RXSUM_s, WTXSUM_s, RHOYSUM_s, RHOFYSUM_s, KYSUM_s, PHIYSUM_s, RYSUM_s, WTYSUM_s, RHOSUM_s, RHOCPSUM_s, ALPHASUM_s, ALPHAFSUM_s, HRSUM_s, GGGPSUM_s, SXXSUM_s, TKSUM_s, PHISUM_s, DMPSUM_s, DHPSUM_s, XWSSUM_s, WTPSUM_s) = setup_interpolated_properties(
            coords
        )

        # Setup threaded target arrays and 2 thread buffers
        (ETA0SUM_t, ETASUM_t, GGGSUM_t, SXYSUM_t, COHSUM_t, TENSUM_t, FRISUM_t, WTSUM_t, RHOXSUM_t, RHOFXSUM_t, KXSUM_t, PHIXSUM_t, RXSUM_t, WTXSUM_t, RHOYSUM_t, RHOFYSUM_t, KYSUM_t, PHIYSUM_t, RYSUM_t, WTYSUM_t, RHOSUM_t, RHOCPSUM_t, ALPHASUM_t, ALPHAFSUM_t, HRSUM_t, GGGPSUM_t, SXXSUM_t, TKSUM_t, PHISUM_t, DMPSUM_t, DHPSUM_t, XWSSUM_t, WTPSUM_t) = setup_interpolated_properties(
            coords
        )

        buffers = allocate_thread_interpolation_buffers(2, coords)

        # Synthesize random marker coordinates inside domain
        xm = rand(rng, marknum) .* (coords.xsize - 2.0*coords.dx) .+ coords.dx
        ym = rand(rng, marknum) .* (coords.ysize - 2.0*coords.dy) .+ coords.dy
        etatotalm = fill(1e18, marknum)
        etavpm = fill(1e18, marknum)
        inv_gggtotalm = fill(1e-10, marknum)
        sxym = rand(rng, marknum) .* 1e6
        sxxm = rand(rng, marknum) .* 1e6
        cohestotalm = fill(1e7, marknum)
        tenstotalm = fill(1e6, marknum)
        fricttotalm = fill(0.6, marknum)
        rhototalm = fill(3000.0, marknum)
        rhofluidcur = fill(1000.0, marknum)
        ktotalm = fill(2.5, marknum)
        phim = fill(0.1, marknum)
        etafluidcur_inv_kphim = fill(1e12, marknum)
        rhocptotalm = fill(3e6, marknum)
        alphasolidcur = fill(3e-5, marknum)
        alphafluidcur = fill(1e-4, marknum)
        hrtotalm = fill(1e-7, marknum)
        tkm_rhocptotalm = fill(9e8, marknum)

        # Serial accumulation
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
                ETA0SUM_s,
                ETASUM_s,
                GGGSUM_s,
                SXYSUM_s,
                COHSUM_s,
                TENSUM_s,
                FRISUM_s,
                WTSUM_s;
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
                RHOXSUM_s,
                RHOFXSUM_s,
                KXSUM_s,
                PHIXSUM_s,
                RXSUM_s,
                WTXSUM_s;
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
                RHOYSUM_s,
                RHOFYSUM_s,
                KYSUM_s,
                PHIYSUM_s,
                RYSUM_s,
                WTYSUM_s;
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
                GGGPSUM_s,
                SXXSUM_s,
                RHOSUM_s,
                RHOCPSUM_s,
                ALPHASUM_s,
                ALPHAFSUM_s,
                HRSUM_s,
                PHISUM_s,
                TKSUM_s,
                WTPSUM_s;
                coords=coords,
            )
        end

        # Two-threaded partition accumulation
        reset_thread_buffers!(buffers)
        mid = marknum ÷ 2
        for m in 1:mid
            buf = buffers[1]
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
        for m in (mid + 1):marknum
            buf = buffers[2]
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
            ETA0SUM_t,
            ETASUM_t,
            GGGSUM_t,
            SXYSUM_t,
            COHSUM_t,
            TENSUM_t,
            FRISUM_t,
            WTSUM_t,
            RHOXSUM_t,
            RHOFXSUM_t,
            KXSUM_t,
            PHIXSUM_t,
            RXSUM_t,
            WTXSUM_t,
            RHOYSUM_t,
            RHOFYSUM_t,
            KYSUM_t,
            PHIYSUM_t,
            RYSUM_t,
            WTYSUM_t,
            RHOSUM_t,
            RHOCPSUM_t,
            ALPHASUM_t,
            ALPHAFSUM_t,
            HRSUM_t,
            GGGPSUM_t,
            SXXSUM_t,
            TKSUM_t,
            PHISUM_t,
            WTPSUM_t,
            buffers,
        )

        # Verify equivalence
        @test ETA0SUM_s ≈ ETA0SUM_t
        @test ETASUM_s ≈ ETASUM_t
        @test GGGSUM_s ≈ GGGSUM_t
        @test SXYSUM_s ≈ SXYSUM_t
        @test WTSUM_s ≈ WTSUM_t
        @test RHOXSUM_s ≈ RHOXSUM_t
        @test RHOYSUM_s ≈ RHOYSUM_t
        @test RHOSUM_s ≈ RHOSUM_t
        @test RHOCPSUM_s ≈ RHOCPSUM_t
        @test SXXSUM_s ≈ SXXSUM_t
        @test WTPSUM_s ≈ WTPSUM_t
    end

    @testset "Multi-criterion adaptive timestepping" begin
        coords = GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        vx = zeros(18, 18)
        vy = zeros(18, 18)
        vxf = zeros(18, 18)
        vyf = zeros(18, 18)

        # Zero velocity -> dt should not be limited by CFL
        dt0 = 1e8
        dt_cfl = compute_adaptive_timestep(
            vx, vy, vxf, vyf, dt0, 0.0; coords=coords, dt_longest_val=1e10
        )
        @test dt_cfl ≈ dt0

        # Fast velocity triggers CFL limiter: dt <= dxymax * dx / max(vx)
        vx .= 1e-3 # 1 mm/s = very fast in geodynamics
        dt_cfl_fast = compute_adaptive_timestep(
            vx, vy, vxf, vyf, dt0, 0.0; coords=coords, dt_longest_val=1e10
        )
        @test dt_cfl_fast < dt0
        @test dt_cfl_fast <= 0.1 * coords.dx / 1e-3 + 1e-6

        # Fast compaction triggers aphimax limiter
        vx .= 0.0
        aphimax = 1e-4
        dt_compaction = compute_adaptive_timestep(
            vx,
            vy,
            vxf,
            vyf,
            dt0,
            aphimax;
            coords=coords,
            dphimax_val=0.01,
            dt_longest_val=1e10,
        )
        @test dt_compaction <= 0.01 / aphimax + 1e-6

        # Temperature variation limiter (maxDTcurrent > DTmax)
        dt_thermal = compute_adaptive_timestep(
            vx,
            vy,
            vxf,
            vyf,
            dt0,
            0.0;
            coords=coords,
            maxDTcurrent=50.0,
            DTmax_val=10.0,
            dt_longest_val=1e10,
        )
        @test dt_thermal <= dt0 * (10.0 / 50.0)

        # Minimum dt clamping
        dt_clamped = compute_adaptive_timestep(
            vx, vy, vxf, vyf, 0.01, 0.0; coords=coords, dt_min=5.0, dt_longest_val=1e10
        )
        @test dt_clamped ≈ 5.0 rtol=1e-12

        # Custom dx_val / dy_val override
        vx[2, 2] = 1.0e-3
        dt_default_dx = compute_adaptive_timestep(
            vx, vy, vxf, vyf, dt0, 0.0; coords=coords, dt_longest_val=1e10
        )
        dt_custom_dx = compute_adaptive_timestep(
            vx,
            vy,
            vxf,
            vyf,
            dt0,
            0.0;
            coords=coords,
            dx_val=1000.0,
            dy_val=1000.0,
            dt_longest_val=1e10,
        )
        @test dt_custom_dx < dt_default_dx
        @test dt_custom_dx <= 0.1 * 1000.0 / 1.0e-3 + 1e-6
    end

    @testset "Live Threads.@threads execution with threadid() bounds safety" begin
        coords = GridCoordinates(17, 17; xsize=140_000.0, ysize=140_000.0)
        buffers = allocate_thread_interpolation_buffers(coords)
        @test length(buffers) >= max(Threads.nthreads(), Threads.maxthreadid())

        marknum = 500
        rng = MersenneTwister(123)
        xm = rand(rng, marknum) .* (coords.xsize - 2.0*coords.dx) .+ coords.dx
        ym = rand(rng, marknum) .* (coords.ysize - 2.0*coords.dy) .+ coords.dy
        etatotalm = fill(1e18, marknum)
        etavpm = fill(1e18, marknum)
        inv_gggtotalm = fill(1e-10, marknum)
        sxym = rand(rng, marknum) .* 1e6
        cohestotalm = fill(1e7, marknum)
        tenstotalm = fill(1e6, marknum)
        fricttotalm = fill(0.6, marknum)

        reset_thread_buffers!(buffers)
        # Verify static chunked threading execution with buffer reduction
        nchunks = length(buffers)
        chunk_size = cld(marknum, nchunks)
        Threads.@threads :static for c in 1:nchunks
            m_start = (c - 1) * chunk_size + 1
            m_end = min(c * chunk_size, marknum)
            buf = buffers[c]
            for m in m_start:m_end
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
            end
        end

        # Test move_markers_rk4! with multi-threading
        tm = fill(1, marknum)
        tkm = fill(300.0, marknum)
        phim = fill(0.1, marknum)
        sxxm = zeros(marknum)
        vx = fill(1.0e-9, 18, 18)
        vy = fill(-1.0e-9, 18, 18)
        vxf = copy(vx)
        vyf = copy(vy)
        wyx = zeros(17, 17)
        tk2 = fill(300.0, 18, 18)
        xm_orig = copy(xm)
        ym_orig = copy(ym)
        Erebus.Particles.move_markers_rk4!(
            xm,
            ym,
            tm,
            tkm,
            phim,
            sxym,
            sxxm,
            vx,
            vy,
            vxf,
            vyf,
            wyx,
            tk2,
            marknum,
            1.0e8,
            9;
            coords=coords,
        )
        @test all(xm .!= xm_orig)
        @test all(ym .!= ym_orig)
    end

    @testset "Non-compounding adaptive timestepping across repeated iterations" begin
        # Simulate repeated Picard iterations with constant maxDTcurrent > DTmax
        dt = 1.0e11
        dt_ref = dt
        vx = zeros(18, 18)
        vy = zeros(18, 18)
        vxf = zeros(18, 18)
        vyf = zeros(18, 18)
        aphimax = 0.0
        maxDTcurrent = 100.0
        DTmax = 20.0

        # Iteration 1
        dt1 = compute_adaptive_timestep(
            vx,
            vy,
            vxf,
            vyf,
            dt,
            aphimax;
            dt_ref=dt_ref,
            maxDTcurrent=maxDTcurrent,
            DTmax_val=DTmax,
        )
        # Next 10 iterations within the same timestep: dt must remain stable and not collapse
        dt_cur = dt1
        for _ in 1:10
            dt_cur = compute_adaptive_timestep(
                vx,
                vy,
                vxf,
                vyf,
                dt_cur,
                aphimax;
                dt_ref=dt_ref,
                maxDTcurrent=maxDTcurrent,
                DTmax_val=DTmax,
            )
        end
        @test dt_cur ≈ dt1
        @test dt_cur > 1.0
        @test dt_cur ≈ dt_ref * (DTmax / maxDTcurrent)
    end
end

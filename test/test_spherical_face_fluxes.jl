using Test
using Erebus
using Erebus.Config
using ExtendableSparse
using SparseArrays

# Face multipliers of P cell (i, j), as in the continuity rows.
function _sff_m(i, j, coords, xc, yc)
    return Erebus.flow_metric_face_weights(i, j, coords; xcenter=xc, ycenter=yc)
end

# Open surface faces of each rock node: sum of mF and face count, written from the grid.
function _sff_face_sums(coords, rplanet, xc, yc)
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    msum = zeros(Ny1, Nx1)
    nf = zeros(Int, Ny1, Nx1)
    rock(i, j) = hypot(coords.xp[j] - xc, coords.yp[i] - yc) <= rplanet
    inner(i, j) = 2 <= i <= Ny1 - 1 && 2 <= j <= Nx1 - 1
    side = Dict(:W => 1, :E => 2, :N => 3, :S => 4)
    seen = Set{Symbol}()
    for j in 1:Nx1, i in 1:Ny1
        rock(i, j) && inner(i, j) || continue
        m = _sff_m(i, j, coords, xc, yc)
        for (s, di, dj) in ((:W, 0, -1), (:E, 0, 1), (:N, -1, 0), (:S, 1, 0))
            rock(i + di, j + dj) && continue
            msum[i, j] += m[side[s]]
            nf[i, j] += 1
            push!(seen, s)
        end
    end
    @assert seen == Set((:W, :E, :N, :S))
    return msum, nf
end

function _sff_vent(coords, rplanet, xc, yc; kw...)
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    n = 6 * Ny1 * Nx1
    L = ExtendableSparseMatrix{Float64,Int64}(n, n)
    R = zeros(n)
    S = zeros(Ny1, Nx1)
    Erebus.apply_venting_surface_boundary!(
        L,
        R,
        fill(250.0, Ny1, Nx1),
        coords,
        rplanet,
        xc,
        yc,
        10.0;
        mode=:darcy_sink,
        pf=fill(2.0e6, Ny1, Nx1),
        pr=fill(2.5e6, Ny1, Nx1),
        TEN=fill(6.0e7, Ny1, Nx1),
        Kcont=1.0e20,
        S_vent_out=S,
        kw...,
    )
    flush!(L)
    return sparse(L), R, S
end

@testset "Spherical face fluxes" begin
    @testset "venting face weights" begin
        coords = Erebus.GridCoordinates(21, 21; xsize=140_000.0, ysize=140_000.0)
        rp, xc, yc = 50_000.0, 71_300.0, 68_900.0
        msum, nf = _sff_face_sums(coords, rp, xc, yc)
        @test count(>(0), nf) > 20
        @test any(==(2), nf)
        L0, R0, S0 = _sff_vent(coords, rp, xc, yc)
        L1, R1, S1 = _sff_vent(coords, rp, xc, yc; spherical_flow=true)
        k(i, j) = ((j - 1) * coords.Ny1 + i - 1) * 6 + 6
        for j in 1:coords.Nx1, i in 1:coords.Ny1
            if nf[i, j] == 0
                @test iszero(S0[i, j])
                @test iszero(S1[i, j])
                continue
            end
            f = msum[i, j] / nf[i, j]
            @test S1[i, j] ≈ f * S0[i, j] rtol = 1e-12
            @test L1[k(i, j), k(i, j)] ≈ f * L0[k(i, j), k(i, j)] rtol = 1e-12
            @test R1[k(i, j)] ≈ f * R0[k(i, j)] rtol = 1e-12
        end
        @test S0[nf .> 0] != S1[nf .> 0]
        # Option off unchanged: the keyword false is the default path.
        Lf, Rf, Sf = _sff_vent(coords, rp, xc, yc; spherical_flow=false)
        @test isequal(Lf, L0) && isequal(Rf, R0) && isequal(Sf, S0)
    end

    @testset "venting drain limit uses the weighted rate" begin
        coords = Erebus.GridCoordinates(21, 21; xsize=140_000.0, ysize=140_000.0)
        rp, xc, yc = 50_000.0, 71_300.0, 68_900.0
        msum, nf = _sff_face_sums(coords, rp, xc, yc)
        _, _, S0 = _sff_vent(coords, rp, xc, yc)
        _, _, S1 = _sff_vent(coords, rp, xc, yc; spherical_flow=true)
        # A limit between the unweighted and the weighted rate of each node.
        dt, phimin = 1.0e10, 1.0e-4
        smax = 0.5 .* (S0 .+ S1)
        PHI = phimin .+ smax .* dt
        _, _, S0l = _sff_vent(coords, rp, xc, yc; PHI, phimin, dt)
        _, _, S1l = _sff_vent(coords, rp, xc, yc; PHI, phimin, dt, spherical_flow=true)
        # Outer faces have mF > 1, so the weighted rate is the larger one.
        lo = nf .> 0
        @test all(S1[lo] .> S0[lo])
        @test S0l[lo] ≈ S0[lo] rtol = 1e-12
        @test S1l[lo] ≈ smax[lo] rtol = 1e-12
    end

    @testset "solid divergence of the compaction pressure" begin
        div!(D, vx, vy, coords, on, xc, yc) =
            Erebus._solid_divergence_cell!(D, vx, vy, coords, on, xc, yc)
        function fields(N, f)
            coords = Erebus.GridCoordinates(N, N; xsize=100_000.0, ysize=96_000.0)
            xc, yc = 51_300.0, 47_900.0
            Ny1, Nx1 = coords.Ny1, coords.Nx1
            vx = [f(coords.xvx[j] - xc, coords.yvx[i] - yc)[1] for i in 1:Ny1, j in 1:Nx1]
            vy = [f(coords.xvy[j] - xc, coords.yvy[i] - yc)[2] for i in 1:Ny1, j in 1:Nx1]
            return coords, xc, yc, vx, vy
        end
        a = 3.0e-14
        coords, xc, yc, vx, vy = fields(17, (X, Y) -> (a * X, a * Y))
        D = zeros(coords.Ny, coords.Nx)
        div!(D, vx, vy, coords, true, xc, yc)
        @test all(d -> isapprox(d, 3a; rtol=1e-10), D)
        div!(D, vx, vy, coords, false, xc, yc)
        @test all(d -> isapprox(d, 2a; rtol=1e-10), D)

        # Radial v = C r / r^3 has zero 3-D divergence and 2-D divergence -C / r^3.
        C = 1.0e-3
        err = Float64[]
        for N in (32, 64, 128)
            coords, xc, yc, vx, vy = fields(
                N, (X, Y) -> (C * X / hypot(X, Y)^3, C * Y / hypot(X, Y)^3)
            )
            D1 = zeros(coords.Ny, coords.Nx)
            D0 = zeros(coords.Ny, coords.Nx)
            div!(D1, vx, vy, coords, true, xc, yc)
            div!(D0, vx, vy, coords, false, xc, yc)
            r = [
                hypot(coords.xp[j + 1] - xc, coords.yp[i + 1] - yc) for
                i in 1:coords.Ny, j in 1:coords.Nx
            ]
            band = (r .> 15_000.0) .& (r .< 40_000.0)
            ref = sqrt(sum(abs2, C ./ r[band] .^ 3) / count(band))
            push!(err, sqrt(sum(abs2, D1[band]) / count(band)) / ref)
            e0 = sqrt(sum(abs2, D0[band] .+ C ./ r[band] .^ 3) / count(band)) / ref
            @test e0 < 0.05
        end
        @test err[3] < 1e-2
        @test 3.4 < err[1] / err[2] < 4.6
        @test 3.4 < err[2] / err[3] < 4.6

        # The option needs grid-size velocities.
        @test_throws ArgumentError div!(
            zeros(coords.Ny, coords.Nx), vx[1:(end - 1), :], vy, coords, true, xc, yc
        )
    end

    @testset "segregation step uses the weighted divergence" begin
        coords = Erebus.GridCoordinates(17, 17; xsize=100_000.0, ysize=100_000.0)
        xc, yc = 51_300.0, 48_700.0
        a = 3.0e-14
        Ny1, Nx1 = coords.Ny1, coords.Nx1
        vx = [a * (coords.xvx[j] - xc) for i in 1:Ny1, j in 1:Nx1]
        vy = [a * (coords.yvy[i] - yc) for i in 1:Ny1, j in 1:Nx1]
        cfg = MagmaTransportConfig(; active=true, compaction_active=true)
        function run(on)
            ws = Erebus.MagmaSegregationWorkspace(coords.Ny, coords.Nx)
            Erebus.apply_silicate_melt_segregation!(
                [xc],
                [yc],
                [1],
                [300.0],
                [0.0],
                1,
                1.0e10,
                cfg;
                coords,
                xcenter=xc,
                ycenter=yc,
                rplanet=40_000.0,
                vx,
                vy,
                spherical_flow=on,
                workspace=ws,
            )
            return copy(ws.div_v)
        end
        @test all(d -> isapprox(d, 3a; rtol=1e-10), run(true))
        @test all(d -> isapprox(d, 2a; rtol=1e-10), run(false))
        @test_throws ArgumentError Erebus.apply_silicate_melt_segregation!(
            [xc],
            [yc],
            [1],
            [300.0],
            [0.0],
            1,
            1.0e10,
            cfg;
            coords,
            xcenter=xc,
            ycenter=yc,
            div_v=zeros(coords.Ny, coords.Nx),
            spherical_flow=true,
        )
    end

    @testset "run with venting and spherical_flow" begin
        cfg = override_config(
            load_config(joinpath(@__DIR__, "..", "configs", "test_quick.toml")),
            Dict{String,Any}(
                "geometry.spherical_metric" => true,
                "geometry.spherical_flow" => true,
                "venting.active" => true,
                "time.n_steps" => 2,
                "output.savematstep" => 2,
            ),
        )
        @test Erebus._validate_spherical_flow(cfg) === nothing
        out = mktempdir()
        simulation_loop(cfg; output_path=out)
        st = Erebus.load_simulation_state(joinpath(out, "output_00002.jld2"); cfg=cfg)[1]
        @test st.timestep == 2
        S = st.grids.S_vent_grid
        @test all(isfinite, S)
        @test any(>(0), S)
    end
end

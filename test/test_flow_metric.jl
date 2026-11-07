using Erebus
using Random
using Test

# Spherical divergence on the staggered grid: D = sum over faces of m_F * u_F with the face multipliers of
# flow_metric_face_weights, for a P cell (i, j) with faces at xvx[j-1], xvx[j], yvy[i-1], yvy[i].
function sph_div(coords, i, j, ux, uy; xcenter, ycenter)
    mW, mE, mN, mS, _ = Erebus.flow_metric_face_weights(
        i, j, coords; xcenter=xcenter, ycenter=ycenter
    )
    return (mE * ux[i, j] - mW * ux[i, j - 1]) / coords.dx +
           (mS * uy[i, j] - mN * uy[i - 1, j]) / coords.dy
end

# Face values of a field f(X, Y) -> (ux, uy), X and Y relative to the centre.
function face_field(coords, f; xcenter, ycenter)
    ux = [f(coords.xvx[j] - xcenter, coords.yvx[i] - ycenter)[1] for i in 1:coords.Ny1, j in 1:coords.Nx1]
    uy = [f(coords.xvy[j] - xcenter, coords.yvy[i] - ycenter)[2] for i in 1:coords.Ny1, j in 1:coords.Nx1]
    return ux, uy
end

# Non-square grid (dx != dy, Nx != Ny) of size 1.0 x 1.3.
grid(N) = Erebus.GridCoordinates(N, round(Int, 1.3 * (N - 1)) + 1; xsize=1.0, ysize=1.3)

@testset verbose = true "Spherical flow metric face weights" begin
    c33 = grid(33)
    @test !isapprox(c33.dx, c33.dy; rtol=1e-3)
    centres = (
        (0.47, 0.61),                      # off grid
        (c33.xp[17], c33.yp[20]),          # on a P node
        (c33.xvx[16], c33.yp[22]),         # on a Vx face midpoint
    )

    @testset "identity: u = (X, Y) gives D = 3 (flux form)" begin
        for (xc, yc) in centres
            ux, uy = face_field(c33, (X, Y) -> (X, Y); xcenter=xc, ycenter=yc)
            err = maximum(
                abs(sph_div(c33, i, j, ux, uy; xcenter=xc, ycenter=yc) - 3.0) for
                j in 3:(c33.Nx - 1), i in 3:(c33.Ny - 1)
            )
            @test isapprox(err, 0.0; atol=1e-12)
        end
    end

    @testset "radial field u = rhat r^2: div3D = 4r, second order" begin
        xc, yc = 0.47, 0.61
        f(X, Y) = (r = hypot(X, Y); (X * r, Y * r))
        errs = map((33, 65, 129)) do N
            c = grid(N)
            ux, uy = face_field(c, f; xcenter=xc, ycenter=yc)
            e = 0.0
            for j in 3:(c.Nx - 1), i in 3:(c.Ny - 1)
                r = hypot(c.xp[j] - xc, c.yp[i] - yc)
                0.2 < r < 0.45 || continue
                e = max(e, abs(sph_div(c, i, j, ux, uy; xcenter=xc, ycenter=yc) - 4r) / (4r))
            end
            e
        end
        @test errs[1] < 5e-3
        @test isapprox(errs[1] / errs[2], 4.0; atol=0.3)
        @test isapprox(errs[2] / errs[3], 4.0; atol=0.3)
    end

    @testset "shared faces carry one weight, so the flux form telescopes" begin
        xc, yc = 0.47, 0.61
        wface(i, j) = begin
            mW, mE, mN, mS, wP = Erebus.flow_metric_face_weights(
                i, j, c33; xcenter=xc, ycenter=yc
            )
            (mW * wP, mE * wP, mN * wP, mS * wP, wP)
        end
        dmax = 0.0
        for j in 3:(c33.Nx - 2), i in 3:(c33.Ny - 1)
            dmax = max(dmax, abs(wface(i, j)[2] - wface(i, j + 1)[1]))
        end
        for j in 3:(c33.Nx - 1), i in 3:(c33.Ny - 2)
            dmax = max(dmax, abs(wface(i, j)[4] - wface(i + 1, j)[3]))
        end
        @test isapprox(dmax, 0.0; atol=1e-14)

        # Sum of wP D dx dy over the block equals the boundary flux of w u.
        rng = MersenneTwister(7)
        ux = randn(rng, c33.Ny1, c33.Nx1)
        uy = randn(rng, c33.Ny1, c33.Nx1)
        js, is = 3:(c33.Nx - 1), 3:(c33.Ny - 1)
        vol = sum(
            wface(i, j)[5] * sph_div(c33, i, j, ux, uy; xcenter=xc, ycenter=yc) * c33.dx * c33.dy
            for j in js, i in is
        )
        w(x, y) = hypot(x - xc, y - yc)
        bnd =
            sum(
                (w(c33.xvx[last(js)], c33.yp[i]) * ux[i, last(js)] -
                 w(c33.xvx[first(js) - 1], c33.yp[i]) * ux[i, first(js) - 1]) * c33.dy for i in is
            ) + sum(
                (w(c33.xp[j], c33.yvy[last(is)]) * uy[last(is), j] -
                 w(c33.xp[j], c33.yvy[first(is) - 1]) * uy[first(is) - 1, j]) * c33.dx for j in js
            )
        scale = sum(abs, ux) * c33.dy + sum(abs, uy) * c33.dx
        @test isapprox(vol, bnd; atol=1e-12 * scale)
    end

    @testset "positivity of wP and the guard" begin
        rng = MersenneTwister(1)
        m = Inf
        for _ in 1:3000
            xc, yc = rand(rng) * c33.xsize, rand(rng) * c33.ysize
            for j in 3:(c33.Nx - 1), i in 3:(c33.Ny - 1)
                wP = Erebus.flow_metric_face_weights(i, j, c33; xcenter=xc, ycenter=yc)[5]
                m = min(m, wP / min(c33.dx, c33.dy))
            end
        end
        @test isapprox(m, 0.3349; rtol=1e-3)
        # Inconsistent faces (east face west of the west face) give wP < 0.
        @test_throws DomainError Erebus.flow_metric_face_weights(1.0, -1.0, 0.0, 0.5, 0.5, 0.5, 0.1, 0.1)
        @test_throws DomainError Erebus.flow_metric_face_weights(NaN, 0.1, 0.05, 0.0, 0.1, 0.05, 0.1, 0.1)
    end
end

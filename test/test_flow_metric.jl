using Erebus
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
    ux = [
        f(coords.xvx[j] - xcenter, coords.yvx[i] - ycenter)[1] for
        i in 1:coords.Ny1, j in 1:coords.Nx1
    ]
    uy = [
        f(coords.xvy[j] - xcenter, coords.yvy[i] - ycenter)[2] for
        i in 1:coords.Ny1, j in 1:coords.Nx1
    ]
    return ux, uy
end

# Non-square grid: dy / dx = 1.7 / 1.3, so a swap of dx and dy changes every result.
grid(N) = Erebus.GridCoordinates(N, round(Int, 1.3 * (N - 1)) + 1; xsize=1.0, ysize=1.7)

@testset verbose = true "Spherical flow metric face weights" begin
    c33 = grid(33)
    @test isapprox(c33.dy / c33.dx, 1.7 / 1.3; rtol=0.01)
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
            # The plain 2-D divergence of the same field is 2: the metric adds u_r / r = 1.
            err2 = maximum(
                abs(
                    (ux[i, j] - ux[i, j - 1]) / c33.dx +
                    (uy[i, j] - uy[i - 1, j]) / c33.dy - 2.0,
                ) for j in 3:(c33.Nx - 1), i in 3:(c33.Ny - 1)
            )
            @test isapprox(err2, 0.0; atol=1e-12)
        end
    end

    @testset "radial field u = rhat r^2: div3D = 4r, second order" begin
        xc, yc = 0.47, 0.61
        f(X, Y) = (r=hypot(X, Y); (X * r, Y * r))
        errs = map((33, 65, 129)) do N
            c = grid(N)
            ux, uy = face_field(c, f; xcenter=xc, ycenter=yc)
            e = 0.0
            for j in 3:(c.Nx - 1), i in 3:(c.Ny - 1)
                r = hypot(c.xp[j] - xc, c.yp[i] - yc)
                0.2 < r < 0.45 || continue
                e = max(
                    e, abs(sph_div(c, i, j, ux, uy; xcenter=xc, ycenter=yc) - 4r) / (4r)
                )
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
        # A cell-centred weight instead of the face weight breaks the shared face.
        @test isapprox(wface(5, 6)[2], hypot(c33.xvx[6] - xc, c33.yp[5] - yc); rtol=1e-14)
    end

    @testset "rigid rotation u = (-Y, X): div3D = 0, second order away from the centre (N 65 to 257)" begin
        xc, yc = 0.47, 0.61
        errs = map((65, 129, 257)) do N
            c = grid(N)
            ux, uy = face_field(c, (X, Y) -> (-Y, X); xcenter=xc, ycenter=yc)
            maximum(
                abs(sph_div(c, i, j, ux, uy; xcenter=xc, ycenter=yc)) for
                j in 3:(c.Nx - 1), i in 3:(c.Ny - 1) if
                0.2 < hypot(c.xp[j] - xc, c.yp[i] - yc) < 0.45
            )
        end
        @test errs[1] < 5e-4
        @test isapprox(errs[1] / errs[2], 4.0; atol=0.3)
        @test isapprox(errs[2] / errs[3], 4.0; atol=0.3)
    end

    @testset "lower bound of wP and the guard" begin
        # The minimum (dx + dy) / 6 is reached with the centre on a P node.
        lo = (c33.dx + c33.dy) / 6
        wnode = Erebus.flow_metric_face_weights(
            20, 17, c33; xcenter=c33.xp[17], ycenter=c33.yp[20]
        )[5]
        @test isapprox(wnode, lo; rtol=1e-12)
        m = Inf
        for a in 0:0.1:1, b in 0:0.1:1
            xc, yc = c33.xp[17] + a * c33.dx, c33.yp[20] + b * c33.dy
            for j in 3:(c33.Nx - 1), i in 3:(c33.Ny - 1)
                m = min(
                    m, Erebus.flow_metric_face_weights(i, j, c33; xcenter=xc, ycenter=yc)[5]
                )
            end
        end
        @test isapprox(m, lo; rtol=1e-12)
        ok = (-0.05, 0.05, 0.0, -0.05, 0.05, 0.0, 0.1, 0.1)
        @test all(isfinite, Erebus.flow_metric_face_weights(ok...))
        bad(k, v) = ntuple(n -> n == k ? v : ok[n], 8)
        @test_throws DomainError Erebus.flow_metric_face_weights(bad(1, NaN)...)
        @test_throws DomainError Erebus.flow_metric_face_weights(bad(6, Inf)...)
        @test_throws DomainError Erebus.flow_metric_face_weights(bad(7, 0.0)...)
        @test_throws DomainError Erebus.flow_metric_face_weights(bad(8, -0.1)...)
        @test_throws DomainError Erebus.flow_metric_face_weights(bad(2, -0.05)...)
        @test_throws DomainError Erebus.flow_metric_face_weights(bad(5, -0.05)...)
        @test Erebus.Physics.flow_metric_face_weights === Erebus.flow_metric_face_weights
    end
end

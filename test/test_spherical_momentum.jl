using Erebus
using LinearAlgebra
using SparseArrays
using Test

# Material fields for one assembly; viscous limit (G dt >> eta), no gravity, no history.
function _smo_fields(coords, eta_b, eta_p)
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    z = zeros(Ny1, Nx1)
    return (
        ETA=eta_b,
        ETAP=eta_p,
        GGG=fill(1.0e40, Ny1, Nx1),
        GGGP=fill(1.0e40, Ny1, Nx1),
        SXY0=copy(z),
        SXX0=copy(z),
        RHOX=fill(3000.0, Ny1, Nx1),
        RHOY=fill(3000.0, Ny1, Nx1),
        RHOFX=fill(1000.0, Ny1, Nx1),
        RHOFY=fill(1000.0, Ny1, Nx1),
        RX=fill(1.0e8, Ny1, Nx1),
        RY=fill(1.0e8, Ny1, Nx1),
        ETAPHI=fill(1.0e19, Ny1, Nx1),
        BETAPHI=fill(1.0e-11, Ny1, Nx1),
        PHI=fill(0.1, Ny1, Nx1),
        gx=copy(z),
        gy=copy(z),
        pr0=copy(z),
        pf0=copy(z),
        DMP=copy(z),
    )
end

function _smo_assemble(nvar, coords, f; kwargs...)
    args = (
        f.ETA,
        f.ETAP,
        f.GGG,
        f.GGGP,
        f.SXY0,
        f.SXX0,
        f.RHOX,
        f.RHOY,
        f.RHOFX,
        f.RHOFY,
        f.RX,
        f.RY,
        f.ETAPHI,
        f.BETAPHI,
        f.PHI,
        f.gx,
        f.gy,
        f.pr0,
        f.pf0,
        f.DMP,
        1.0e10,
    )
    R = zeros(coords.Ny1 * coords.Nx1 * nvar)
    L = if nvar == 6
        assemble_hydromechanical_lse!(args..., R; coords=coords, kwargs...)
    else
        assemble_hydromechanical_4var_lse!(args..., R; coords=coords, kwargs...)
    end
    return sparse(L), R
end

function _smo_grid(N; aspect=1.17)
    Ny = round(Int, N * aspect)
    coords = GridCoordinates(
        N, Ny; xsize=140000.0, ysize=140000.0 * aspect * (N - 1) / (Ny - 1)
    )
    xc = 0.5 * coords.xsize + 0.37 * coords.dx
    yc = 0.5 * coords.ysize - 0.21 * coords.dy
    return coords, xc, yc
end

# Continuum pseudo-3-D momentum (1/w) d_j(w tau_ij) - tau_33 d_i w / w for w = 2r, from
# closed-form velocity gradients and central differences of the stress (step h).
function _smo_continuum(vel, grad, eta, x, y, xc, yc; h=1.0)
    function tau(x, y)
        X, Y = x - xc, y - yc
        r2 = X^2 + Y^2
        (vx, vy) = vel(x, y)
        (dxvx, dyvx, dxvy, dyvy) = grad(x, y)
        e33 = (X * vx + Y * vy) / r2
        d3 = (dxvx + dyvy + e33) / 3
        e = 2 * eta(x, y)
        return (e * (dxvx - d3), e * (dyvy - d3), e * 0.5 * (dyvx + dxvy), e * (e33 - d3))
    end
    X, Y = x - xc, y - yc
    r2 = X^2 + Y^2
    t = tau(x, y)
    dtxx = (tau(x + h, y)[1] - tau(x - h, y)[1]) / (2h)
    dtxy_y = (tau(x, y + h)[3] - tau(x, y - h)[3]) / (2h)
    dtyy = (tau(x, y + h)[2] - tau(x, y - h)[2]) / (2h)
    dtxy_x = (tau(x + h, y)[3] - tau(x - h, y)[3]) / (2h)
    fx = dtxx + dtxy_y + (t[1] * X + t[3] * Y) / r2 - t[4] * X / r2
    fy = dtxy_x + dtyy + (t[3] * X + t[2] * Y) / r2 - t[4] * Y / r2
    return fx, fy
end

# Relative rms of (L u)_momentum - f over the band rmin < r < rmax, for both rows.
function _smo_error(nvar, N, vel, grad, eta; on=true, rmin=15000.0, rmax=50000.0)
    coords, xc, yc = _smo_grid(N)
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    eta_b = [eta(coords.xvx[j], coords.yvy[i]) for i in 1:Ny1, j in 1:Nx1]
    eta_p = [eta(coords.xp[j], coords.yp[i]) for i in 1:Ny1, j in 1:Nx1]
    f = _smo_fields(coords, eta_b, eta_p)
    L, _ = _smo_assemble(nvar, coords, f; spherical_flow=on, xcenter=xc, ycenter=yc)
    u = zeros(size(L, 1))
    for j in 1:Nx1, i in 1:Ny1
        k = ((j - 1) * Ny1 + i - 1) * nvar
        u[k + 1] = vel(coords.xvx[j], coords.yvx[i])[1]
        u[k + 2] = vel(coords.xvy[j], coords.yvy[i])[2]
    end
    Lu = L * u
    num, den = 0.0, 0.0
    for j in 2:(Nx1 - 2), i in 2:(Ny1 - 2)
        k = ((j - 1) * Ny1 + i - 1) * nvar
        for (row, x, y, c) in (
            (k + 1, coords.xvx[j], coords.yvx[i], 1),
            (k + 2, coords.xvy[j], coords.yvy[i], 2),
        )
            rmin < hypot(x - xc, y - yc) < rmax || continue
            fc = _smo_continuum(vel, grad, eta, x, y, xc, yc)
            num += (Lu[row] - fc[c])^2
            den += fc[c]^2
        end
    end
    return sqrt(num / den), coords, xc, yc
end

@testset "Spherical momentum operator" begin
    xc0, yc0 = 0.5 * 140000.0, 0.5 * 140000.0 * 1.17
    L0 = 40000.0
    @testset "Radial flow v = C r_hat / r^2 has a uniform total pressure" begin
        # The sphere has zero stress divergence for this flow; scale by a 1/r^4 term.
        for nvar in (6, 4)
            errs = Float64[]
            for N in (33, 65, 129)
                coords, xc, yc = _smo_grid(N)
                C = 1.0e-9 * L0^3
                vel(x, y) = (
                    C * (x - xc) / hypot(x - xc, y - yc)^3,
                    C * (y - yc) / hypot(x - xc, y - yc)^3,
                )
                function grad(x, y)
                    X, Y = x - xc, y - yc
                    r = hypot(X, Y)
                    return (
                        C * (r^2 - 3X^2) / r^5,
                        -3C * X * Y / r^5,
                        -3C * X * Y / r^5,
                        C * (r^2 - 3Y^2) / r^5,
                    )
                end
                eta(x, y) = 1.0e19
                # Scale: rms of |d_x tau_xx| alone, since the full residual is 0 here.
                Lr = _smo_assemble(
                    nvar,
                    coords,
                    _smo_fields(
                        coords,
                        fill(1.0e19, coords.Ny1, coords.Nx1),
                        fill(1.0e19, coords.Ny1, coords.Nx1),
                    );
                    spherical_flow=true,
                    xcenter=xc,
                    ycenter=yc,
                )[1]
                u = zeros(size(Lr, 1))
                for j in 1:coords.Nx1, i in 1:coords.Ny1
                    k = ((j - 1) * coords.Ny1 + i - 1) * nvar
                    u[k + 1] = vel(coords.xvx[j], coords.yvx[i])[1]
                    u[k + 2] = vel(coords.xvy[j], coords.yvy[i])[2]
                end
                Lu = Lr * u
                num, den = 0.0, 0.0
                for j in 2:(coords.Nx1 - 2), i in 2:(coords.Ny1 - 2)
                    k = ((j - 1) * coords.Ny1 + i - 1) * nvar
                    for (row, x, y) in (
                        (k + 1, coords.xvx[j], coords.yvx[i]),
                        (k + 2, coords.xvy[j], coords.yvy[i]),
                    )
                        r = hypot(x - xc, y - yc)
                        15000.0 < r < 50000.0 || continue
                        num += Lu[row]^2
                        den += (12 * 1.0e19 * C / r^4)^2
                    end
                end
                push!(errs, sqrt(num / den))
            end
            for k in 1:2
                @test 3.6 < errs[k] / errs[k + 1] < 4.4
            end
            @test errs[3] < 1.0e-2
        end
    end

    @testset "Smooth non-radial flow with variable viscosity" begin
        U = 1.0e-9
        for nvar in (6, 4)
            errs = Float64[]
            local coords, xc, yc
            for N in (33, 65, 129)
                coords, xc, yc = _smo_grid(N)
                vel(x, y) = (
                    U * (1 + (x - xc) / L0 + (x - xc) * (y - yc) / L0^2),
                    U * ((y - yc)^2 / L0^2 - (x - xc) / L0),
                )
                grad(x, y) = (
                    U * (1 / L0 + (y - yc) / L0^2),
                    U * (x - xc) / L0^2,
                    -U / L0,
                    U * 2 * (y - yc) / L0^2,
                )
                eta(x, y) = 1.0e19 * (1 + 0.3 * (x - xc) / L0 + 0.2 * (y - yc) / L0)
                push!(errs, _smo_error(nvar, N, vel, grad, eta)[1])
            end
            for k in 1:2
                @test 3.6 < errs[k] / errs[k + 1] < 4.4
            end
            @test errs[3] < 1.0e-2
        end
    end

    @testset "Option off keeps the 2-D momentum" begin
        coords, xc, yc = _smo_grid(65)
        U = 1.0e-9
        vel(x, y) = (U * (1 + (x - xc) / L0), U * (-(x - xc) / L0))
        grad(x, y) = (U / L0, 0.0, -U / L0, 0.0)
        eta(x, y) = 1.0e19
        # The pseudo-3-D operator differs from the 2-D one by O(1) for a flow through the centre.
        e_off = _smo_error(6, 65, vel, grad, eta; on=false)[1]
        e_on = _smo_error(6, 65, vel, grad, eta; on=true)[1]
        @test e_off > 10 * e_on
        @test isapprox(e_on, 0.0; atol=0.05)
    end
end

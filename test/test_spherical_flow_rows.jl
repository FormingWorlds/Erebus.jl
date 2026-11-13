using Erebus
using ExtendableSparse
using LinearAlgebra
using SparseArrays
using Test

# Material fields for one assembly on a grid with dx != dy; RX, RY vary per face.
function _sfr_fields(coords)
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    RX = [1.0e8 * (1.0 + 0.3 * sin(0.7 * i + 1.3 * j)) for i in 1:Ny1, j in 1:Nx1]
    RY = [1.0e8 * (1.0 + 0.3 * cos(0.9 * i + 0.4 * j)) for i in 1:Ny1, j in 1:Nx1]
    return (
        ETA=fill(1.0e19, Ny1, Nx1),
        ETAP=fill(1.0e19, Ny1, Nx1),
        GGG=fill(1.0e10, Ny1, Nx1),
        GGGP=fill(1.0e10, Ny1, Nx1),
        SXY0=zeros(Ny1, Nx1),
        SXX0=zeros(Ny1, Nx1),
        RHOX=fill(3000.0, Ny1, Nx1),
        RHOY=fill(3000.0, Ny1, Nx1),
        RHOFX=fill(1000.0, Ny1, Nx1),
        RHOFY=fill(1000.0, Ny1, Nx1),
        RX=RX,
        RY=RY,
        ETAPHI=fill(1.0e19, Ny1, Nx1),
        BETAPHI=fill(1.0e-11, Ny1, Nx1),
        PHI=fill(0.1, Ny1, Nx1),
        gx=[0.4 + 0.01 * j for i in 1:Ny1, j in 1:Nx1],
        gy=[1.0 + 0.02 * i for i in 1:Ny1, j in 1:Nx1],
        pr0=fill(1.0e6, Ny1, Nx1),
        pf0=fill(1.0e6, Ny1, Nx1),
        DMP=zeros(Ny1, Nx1),
    )
end

function _sfr_assemble(nvar, coords, f; kwargs...)
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

# Face sum D_k(x) = sum over face columns of L[k, c] * (x[c] - x[k]) (Pf, 4-var), or
# sum L[k, c] * x[c] over the face velocity columns (divergence rows).
_sfr_rowsum(L, k, cols, x) = sum(L[k, c] * x[c] for c in cols)
_sfr_diffsum(L, k, cols, x) = sum(L[k, c] * (x[c] - x[k]) for c in cols)

function _sfr_grid(N; aspect=1.17)
    Ny = round(Int, N * aspect)
    coords = GridCoordinates(
        N, Ny; xsize=140000.0, ysize=140000.0 * aspect * (N - 1) / (Ny - 1)
    )
    xc = 0.5 * coords.xsize + 0.37 * coords.dx
    yc = 0.5 * coords.ysize - 0.21 * coords.dy
    return coords, xc, yc
end

@testset "Spherical flow weights in the continuity rows" begin
    coords, xc, yc = _sfr_grid(17)
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    dx, dy = coords.dx, coords.dy
    f = _sfr_fields(coords)
    on = (spherical_flow=true, xcenter=xc, ycenter=yc)
    off = (xcenter=xc, ycenter=yc)
    @test !isapprox(dx, dy; rtol=1.0e-3)

    # Continuity cells: near the centre and far from it; one anchor row; one ghost row.
    i_c = findmin(abs.(coords.yp[3:(Ny1 - 2)] .- yc))[2] + 2
    j_c = findmin(abs.(coords.xp[3:(Nx1 - 2)] .- xc))[2] + 2
    cells = ((i_c, j_c), (4, 4))
    @test hypot(coords.xp[4] - xc, coords.yp[4] - yc) > 5 * dx

    @testset "6-variable rows" begin
        L0, R0 = _sfr_assemble(6, coords, f; off...)
        L1, R1 = _sfr_assemble(6, coords, f; on...)
        for (i, j) in cells
            kvx = ((j - 1) * Ny1 + i - 1) * 6 + 1
            kvy, kpm, kqx, kqy, kpf = kvx + 1, kvx + 2, kvx + 3, kvx + 4, kvx + 5
            mW, mE, mN, mS, _ = Erebus.flow_metric_face_weights(
                i, j, coords; xcenter=xc, ycenter=yc
            )
            # Option off: unweighted divergence stencil.
            @test isapprox(L0[kpm, kvx], 1.0 / dx; rtol=1.0e-14)
            @test isapprox(L0[kpf, kqy - 6], -1.0 / dy; rtol=1.0e-14)
            for (row, cols) in (
                (kpm, (kvx - 6 * Ny1, kvx, kvy - 6, kvy)),
                (kpf, (kqx - 6 * Ny1, kqx, kqy - 6, kqy)),
            )
                for (c, m) in zip(cols, (mW, mE, mN, mS))
                    @test isapprox(L1[row, c], m * L0[row, c]; rtol=1.0e-14)
                end
                for c in (kpm, kpf)
                    @test isapprox(L1[row, c], L0[row, c]; rtol=1.0e-14)
                end
            end
            @test isapprox(R1[kpm], R0[kpm]; rtol=1.0e-14, atol=0.0)
            @test isapprox(R1[kpf], R0[kpf]; rtol=1.0e-14, atol=0.0)
        end
        for (i, j) in ((2, 6), (1, 6))
            kpm = ((j - 1) * Ny1 + i - 1) * 6 + 3
            @test isequal(L1[kpm, :], L0[kpm, :])
            @test isequal(L1[kpm + 3, :], L0[kpm + 3, :])
        end
        # Momentum rows are unchanged.
        kvx = ((j_c - 1) * Ny1 + i_c - 1) * 6 + 1
        @test isequal(L1[kvx, :], L0[kvx, :])
        @test isequal(L1[kvx + 1, :], L0[kvx + 1, :])
    end

    @testset "4-variable rows and the gravity part of the right-hand side" begin
        rx0 = zeros(Ny1, Nx1)
        ry0 = zeros(Ny1, Nx1)
        L0, R0 = _sfr_assemble(4, coords, f; off..., rx_eff_out=rx0, ry_eff_out=ry0)
        L1, R1 = _sfr_assemble(4, coords, f; on...)
        for (i, j) in cells
            kvx = ((j - 1) * Ny1 + i - 1) * 4 + 1
            kvy, kpm, kpf = kvx + 1, kvx + 2, kvx + 3
            mW, mE, mN, mS, _ = Erebus.flow_metric_face_weights(
                i, j, coords; xcenter=xc, ycenter=yc
            )
            ms = (mW, mE, mN, mS)
            for (c, m) in zip((kvx - 4 * Ny1, kvx, kvy - 4, kvy), ms)
                @test isapprox(L1[kpm, c], m * L0[kpm, c]; rtol=1.0e-14)
            end
            @test isapprox(L1[kpm, kpm], L0[kpm, kpm]; rtol=1.0e-14)
            fcols = (kpf - 4 * Ny1, kpf + 4 * Ny1, kpf - 4, kpf + 4)
            for (c, m) in zip(fcols, ms)
                @test isapprox(L1[kpf, c], m * L0[kpf, c]; rtol=1.0e-14)
            end
            dsum = sum((m - 1.0) * L0[kpf, c] for (c, m) in zip(fcols, ms))
            @test isapprox(L1[kpf, kpf], L0[kpf, kpf] - dsum; rtol=1.0e-12)
            @test isapprox(L1[kpf, kpm], L0[kpf, kpm]; rtol=1.0e-14)
            # Gravity term of each face flux, from the resistances the assembly used.
            gW = f.RHOFX[i, j - 1] * f.gx[i, j - 1] / (dx * rx0[i, j - 1])
            gE = -f.RHOFX[i, j] * f.gx[i, j] / (dx * rx0[i, j])
            gN = f.RHOFY[i - 1, j] * f.gy[i - 1, j] / (dy * ry0[i - 1, j])
            gS = -f.RHOFY[i, j] * f.gy[i, j] / (dy * ry0[i, j])
            dR = (mW - 1) * gW + (mE - 1) * gE + (mN - 1) * gN + (mS - 1) * gS
            @test !isapprox(dR, 0.0; atol=1.0e-6 * abs(gW))
            @test isapprox(R1[kpf] - R0[kpf], dR; rtol=1.0e-9)
            @test isapprox(R1[kpm], R0[kpm]; rtol=1.0e-14, atol=0.0)
        end
        i, j = 2, 6
        kpf = ((j - 1) * Ny1 + i - 1) * 4 + 4
        @test isequal(L1[kpf, :], L0[kpf, :])
    end

    @testset "Option off leaves the system unchanged" begin
        for nvar in (6, 4)
            La, Ra = _sfr_assemble(nvar, coords, f; off...)
            Lb, Rb = _sfr_assemble(nvar, coords, f; spherical_flow=false, off...)
            @test isequal(La, Lb)
            @test isequal(Ra, Rb)
        end
        @test_throws ArgumentError _sfr_assemble(4, coords, f; spherical_flow=true)
        @test_throws ArgumentError _sfr_assemble(
            6, coords, f; spherical_flow=true, xcenter=xc
        )
        @test_throws ArgumentError _sfr_assemble(
            6, coords, f; spherical_flow=true, ycenter=yc
        )
    end

    @testset "6-variable and 4-variable solutions agree with the option on" begin
        L6, R6 = _sfr_assemble(6, coords, f; on...)
        L4, R4 = _sfr_assemble(4, coords, f; on...)
        S6 = L6 \ R6
        S4 = L4 \ R4
        v6 = [zeros(Ny1, Nx1) for _ in 1:6]
        v4 = [zeros(Ny1, Nx1) for _ in 1:4]
        process_hydromechanical_solution!(S6, v6...; coords=coords)
        process_hydromechanical_4var_solution!(S4, v4...; coords=coords)
        @test norm(v4[1] - v6[1]) / norm(v6[1]) < 1.0e-10
        @test norm(v4[2] - v6[2]) / norm(v6[2]) < 1.0e-10
        @test norm(v4[3] - v6[3]) / norm(v6[3]) < 1.0e-8
        @test norm(v4[4] - v6[6]) / norm(v6[6]) < 1.0e-8
        S6off =
            sparse(_sfr_assemble(6, coords, f; off...)[1]) \
            _sfr_assemble(6, coords, f; off...)[2]
        @test !isapprox(S6off, S6; rtol=1.0e-6)
    end
end

@testset "Spherical divergence of radial fields" begin
    coords, xc, yc = _sfr_grid(33)
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    f = _sfr_fields(coords)
    f = merge(f, (RX=fill(1.0e8, Ny1, Nx1), RY=fill(1.0e8, Ny1, Nx1)))
    L0, _ = _sfr_assemble(4, coords, f; xcenter=xc, ycenter=yc)
    L1, _ = _sfr_assemble(4, coords, f; spherical_flow=true, xcenter=xc, ycenter=yc)
    L60, _ = _sfr_assemble(6, coords, f; xcenter=xc, ycenter=yc)
    L61, _ = _sfr_assemble(6, coords, f; spherical_flow=true, xcenter=xc, ycenter=yc)

    # x: Pf = r^2 on P nodes (4-var); v = (X, Y) on velocity and Darcy nodes.
    x4 = zeros(Ny1 * Nx1 * 4)
    x6 = zeros(Ny1 * Nx1 * 6)
    for j in 1:Nx1, i in 1:Ny1
        k4 = ((j - 1) * Ny1 + i - 1) * 4
        k6 = ((j - 1) * Ny1 + i - 1) * 6
        X, Y = coords.xvx[j] - xc, coords.yvy[i] - yc
        x4[k4 + 1], x4[k4 + 2] = X, Y
        x4[k4 + 4] = (coords.xp[j] - xc)^2 + (coords.yp[i] - yc)^2
        x6[k6 + 1], x6[k6 + 2], x6[k6 + 4], x6[k6 + 5] = X, Y, X, Y
    end
    ratios = Float64[]
    for j in 3:(Nx1 - 2), i in 3:(Ny1 - 2)
        k4 = ((j - 1) * Ny1 + i - 1) * 4
        k6 = ((j - 1) * Ny1 + i - 1) * 6
        vc4 = (k4 + 1 - 4 * Ny1, k4 + 1, k4 + 2 - 4, k4 + 2)
        pc4 = (k4 + 4 - 4 * Ny1, k4 + 4 + 4 * Ny1, k4 + 4 - 4, k4 + 4 + 4)
        vc6 = (k6 + 1 - 6 * Ny1, k6 + 1, k6 + 2 - 6, k6 + 2)
        qc6 = (k6 + 4 - 6 * Ny1, k6 + 4, k6 + 5 - 6, k6 + 5)
        push!(ratios, _sfr_rowsum(L1, k4 + 3, vc4, x4) / _sfr_rowsum(L0, k4 + 3, vc4, x4))
        push!(ratios, _sfr_diffsum(L1, k4 + 4, pc4, x4) / _sfr_diffsum(L0, k4 + 4, pc4, x4))
        push!(ratios, _sfr_rowsum(L61, k6 + 3, vc6, x6) / _sfr_rowsum(L60, k6 + 3, vc6, x6))
        push!(ratios, _sfr_rowsum(L61, k6 + 6, qc6, x6) / _sfr_rowsum(L60, k6 + 6, qc6, x6))
    end
    # Divergence 3 (sphere) against 2 (disk) for (X, Y); Laplacian 6a against 4a for a r^2.
    @test length(ratios) == 4 * (Nx1 - 4) * (Ny1 - 4)
    @test all(r -> isapprox(r, 1.5; rtol=1.0e-9), ratios)
end

@testset "Second order for a non-polynomial radial field" begin
    # Pf = r^4: 3-D Laplacian 20 r^2; face resistances from k = k0 (1 + r^2 / R^2),
    # 3-D divergence of k grad(r^2): k0 (6 + 10 r^2 / R^2).
    errs4 = Float64[]
    errsk = Float64[]
    for N in (33, 65, 129)
        coords, xc, yc = _sfr_grid(N)
        Ny1, Nx1 = coords.Ny1, coords.Nx1
        Rk = 40000.0
        k(x, y) = 1.0 + ((x - xc)^2 + (y - yc)^2) / Rk^2
        RXk = [1.0e8 / k(coords.xvx[j], coords.yp[i]) for i in 1:Ny1, j in 1:Nx1]
        RYk = [1.0e8 / k(coords.xp[j], coords.yvy[i]) for i in 1:Ny1, j in 1:Nx1]
        f = _sfr_fields(coords)
        fu = merge(f, (RX=fill(1.0e8, Ny1, Nx1), RY=fill(1.0e8, Ny1, Nx1)))
        fk = merge(f, (RX=RXk, RY=RYk))
        Lu, _ = _sfr_assemble(4, coords, fu; spherical_flow=true, xcenter=xc, ycenter=yc)
        Lk, _ = _sfr_assemble(4, coords, fk; spherical_flow=true, xcenter=xc, ycenter=yc)
        x2 = zeros(Ny1 * Nx1 * 4)
        x4 = zeros(Ny1 * Nx1 * 4)
        for j in 1:Nx1, i in 1:Ny1
            r2 = (coords.xp[j] - xc)^2 + (coords.yp[i] - yc)^2
            x2[((j - 1) * Ny1 + i - 1) * 4 + 4] = r2
            x4[((j - 1) * Ny1 + i - 1) * 4 + 4] = r2^2
        end
        e4 = 0.0
        ek = 0.0
        for j in 3:(Nx1 - 2), i in 3:(Ny1 - 2)
            r2 = (coords.xp[j] - xc)^2 + (coords.yp[i] - yc)^2
            sqrt(r2) > 15000.0 || continue
            k4 = ((j - 1) * Ny1 + i - 1) * 4 + 4
            pc = (k4 - 4 * Ny1, k4 + 4 * Ny1, k4 - 4, k4 + 4)
            # Scale: the uniform-k operator gives 6 for r^2 exactly.
            s = _sfr_diffsum(Lu, k4, pc, x2) / 6
            e4 = max(e4, abs(_sfr_diffsum(Lu, k4, pc, x4) / s - 20 * r2) / (20 * r2))
            ek = max(
                ek,
                abs(_sfr_diffsum(Lk, k4, pc, x2) / s - (6 + 10 * r2 / Rk^2)) /
                (6 + 10 * r2 / Rk^2),
            )
        end
        push!(errs4, e4)
        push!(errsk, ek)
    end
    for e in (errs4, errsk), k in 1:2
        @test 3.6 < e[k] / e[k + 1] < 4.4
    end
    # Measured at N 129 (r > 15 km): 1.93e-3 and 2.36e-4; ratios 3.95, 4.01 and 3.93, 4.01.
    @test errs4[3] < 2.5e-3
    @test errsk[3] < 3.0e-4
end

@testset "spherical_flow refusals" begin
    base = (geometry=GeometryConfig(spherical_flow=true),)
    @test Erebus._validate_spherical_flow(SimulationConfig(; base...)) === nothing
    bad = (
        (geometry=GeometryConfig(spherical_flow=true, spherical_metric=false),),
        (base..., solver=SolverConfig(hydromech_solver=:iterative, experimental=true)),
        (
            base...,
            solver=SolverConfig(
                hydromech_solver=:matrix_free, experimental=true, darcy_elimination=true
            ),
        ),
        (base..., mpi=MPIConfig(enable=true)),
        (base..., venting=VentingConfig(active=true)),
    )
    for kw in bad
        cfg = SimulationConfig(; kw...)
        @test_throws ArgumentError Erebus._validate_spherical_flow(cfg)
        @test_throws ArgumentError validate_config(cfg)
    end
    @test_throws ArgumentError init_simulation(SimulationConfig(; bad[end]...))
end

@testset "solve_stokes_darcy! passes spherical_flow to the assembly" begin
    base = load_config("configs/test_quick.toml")
    res = Dict{Tuple{Bool,Bool},Any}()
    for elim in (false, true), sf in (false, true)
        cfg = override_config(
            base,
            Dict{String,Any}(
                "geometry.spherical_metric" => true,
                "geometry.spherical_flow" => sf,
                "solver.darcy_elimination" => elim,
                "time.n_steps" => 1,
                "output.savematstep" => 1,
            ),
        )
        out = mktempdir()
        simulation_loop(cfg; output_path=out)
        res[(elim, sf)] = load_state(joinpath(out, "output_00001.jld2"))
    end
    for elim in (false, true)
        on, off = res[(elim, true)], res[(elim, false)]
        @test all(isfinite, on["pf"])
        @test !isapprox(on["pf"], off["pf"]; rtol=1.0e-6)
        @test !isapprox(on["vx"], off["vx"]; rtol=1.0e-6)
    end
    # Both solver forms solve the same weighted system.
    @test isapprox(res[(true, true)]["pf"], res[(false, true)]["pf"]; rtol=1.0e-6)
end

using Test
using Erebus
using Logging

isdefined(@__MODULE__, :PorousShellValidation) ||
    include(joinpath(@__DIR__, "..", "benchmarks", "validation_porous_shell.jl"))

# Bound on |marker T - grid T| for rock markers of the porous shell case at Nx 65 after
# 20 steps [K]: 119 K with dsubgridt 0, 32 K with the default 1 (probe at b61f1363).
const DRIFT_BOUND_K = 60.0

"""
    _subgrid_cell(; n_rock=6, n_air=6)

Rock and sticky-air markers that share the four P nodes around the box centre, with a
uniform node temperature of 300 K. Rock markers start at 400 K to 450 K, air at 290 K.
"""
function _subgrid_cell(; n_rock::Int=6, n_air::Int=6)
    coords = GridCoordinates(GridConfig(; Nx=10, Ny=10, xsize=1.0e5, ysize=1.0e5))
    n = n_rock + n_air
    xc, yc = 0.5 * coords.xsize, 0.5 * coords.ysize
    xm = [xc + 400.0 * cos(2pi * k / n) for k in 1:n]
    ym = [yc + 300.0 * sin(2pi * k / n) for k in 1:n]
    tm = vcat(fill(1, n_rock), fill(3, n_air))
    tkm = vcat(collect(range(400.0, 450.0; length=n_rock)), fill(290.0, n_air))
    tk1 = fill(300.0, coords.Ny1, coords.Nx1)
    return (; coords, n, xm, ym, tm, tkm, tk1)
end

"""
    _rock_heat_change(c, rhocp_rock; kw...)

Run one subgrid step and the grid-to-marker increment of its node counterpart (conduction
increment zero) and return the change of the rock heat content, sum of rhocp dT over rock
markers, with the sum of the absolute marker corrections as its scale.
"""
function _rock_heat_change(c, rhocp_rock; dt=3.0e12, phi_dry=zeros(c.n), kw...)
    tkm = copy(c.tkm)
    DT = zeros(size(c.tk1))
    TKSUM = zeros(size(c.tk1))
    RHOCPSUM = zeros(size(c.tk1))
    Erebus.apply_subgrid_temperature_diffusion!(
        c.xm,
        c.ym,
        c.tm,
        tkm,
        zeros(c.n),
        c.tk1,
        DT,
        TKSUM,
        RHOCPSUM,
        dt,
        c.n,
        1;
        coords=c.coords,
        dsubgridt=1.0,
        phi_dry=phi_dry,
        kw...,
    )
    tk_sub = copy(tkm)
    Erebus.update_marker_temperature!(c.xm, c.ym, tkm, DT, c.tk1, 2, c.n; coords=c.coords)
    rock = findall(<(3), c.tm)
    dE = sum(rhocp_rock[m] * (tkm[m] - c.tkm[m]) for m in rock)
    scale = sum(abs(rhocp_rock[m] * (tk_sub[m] - c.tkm[m])) for m in rock)
    return (; dE, scale, tkm, tk_sub)
end

@testset "Subgrid temperature diffusion" begin
    @testset "default relaxation factor is 1" begin
        # The keyword default of the function is the same factor: it relaxes.
        c = _subgrid_cell(; n_air=0)
        tkm = copy(c.tkm)
        Erebus.apply_subgrid_temperature_diffusion!(
            c.xm,
            c.ym,
            c.tm,
            tkm,
            zeros(c.n),
            c.tk1,
            zeros(size(c.tk1)),
            zeros(size(c.tk1)),
            zeros(size(c.tk1)),
            3.0e12,
            c.n,
            1;
            coords=c.coords,
            phi_dry=zeros(c.n),
        )
        @test !isapprox(tkm, c.tkm; rtol=1e-6)
        tkm1 = copy(c.tkm)
        Erebus.apply_subgrid_temperature_diffusion!(
            c.xm,
            c.ym,
            c.tm,
            tkm1,
            zeros(c.n),
            c.tk1,
            zeros(size(c.tk1)),
            zeros(size(c.tk1)),
            zeros(size(c.tk1)),
            3.0e12,
            c.n,
            1;
            coords=c.coords,
            phi_dry=zeros(c.n),
            dsubgridt=1.0,
        )
        @test isapprox(tkm, tkm1; rtol=1e-14)
        @test isapprox(SolverConfig().dsubgridt, 1.0; rtol=1e-12)
        # Scale guard: a factor of 10 or 0.1 is not the default.
        @test !isapprox(SolverConfig().dsubgridt, 10.0; rtol=1e-3)
        @test !isapprox(SolverConfig().dsubgridt, 0.1; rtol=1e-3)
    end

    @testset "rock heat is kept when sticky air shares the nodes" begin
        c = _subgrid_cell()
        rhocp_rock = fill(Erebus.rhocpsolidm[1], c.n)
        r = _rock_heat_change(c, rhocp_rock)
        # The rock markers relax toward 300 K: the scale is the size of the corrections.
        @test isapprox(
            r.scale,
            sum(Erebus.rhocpsolidm[1] .* abs.(r.tk_sub[1:6] .- c.tkm[1:6]));
            rtol=1e-12,
        )
        @test !isapprox(r.tk_sub[1:6], c.tkm[1:6]; rtol=1e-6)
        @test abs(r.dE) <= 1.0e-12 * r.scale
        # Sticky-air markers keep their temperature in the subgrid step.
        air = findall(>=(3), c.tm)
        @test isapprox(r.tk_sub[air], c.tkm[air]; rtol=0.0, atol=0.0)
        # Bounds: the relaxed rock markers stay between the node and their start values.
        rock = findall(<(3), c.tm)
        @test all(300.0 .<= r.tk_sub[rock] .<= c.tkm[rock])
    end

    @testset "relaxation uses the marker heat capacity and conductivity" begin
        c = _subgrid_cell(; n_air=0)
        dt = 3.0e12
        # Different values per marker, so the heat capacity does not cancel in the node mean.
        rhocp_m = [6.6e6, 4.5e6, 3.3e6, 8.1e6, 3.1e6, 5.2e6]
        k_m = [1.5, 0.4, 3.0, 2.2, 0.9, 5.0]
        r = _rock_heat_change(c, rhocp_m; dt, rhocptotalm=rhocp_m, ktotalm=k_m)
        coords = c.coords
        x = dt * (2.0 * (inv(coords.dx^2) + inv(coords.dy^2)))
        # Expected marker correction before the node counterpart: dT0 (exp(-k x / rhocp) - 1).
        dT0 = c.tkm .- 300.0
        expected = c.tkm .+ dT0 .* (exp.(-k_m .* x ./ rhocp_m) .- 1.0)
        @test isapprox(r.tk_sub, expected; rtol=1e-12)
        # Scale guards: the constant-property relaxation gives a different value.
        k0, c0 = Erebus.ksolidm[1], Erebus.rhocpsolidm[1]
        wrong = c.tkm .+ dT0 .* (exp.(-k0 * x / c0) .- 1.0)
        @test !isapprox(r.tk_sub, wrong; rtol=1e-3)
        @test abs(r.dE) <= 1.0e-12 * r.scale
    end

    @testset "zero marker property falls back to the constitutive value" begin
        c = _subgrid_cell(; n_air=0)
        z = zeros(c.n)
        a = _rock_heat_change(c, fill(Erebus.rhocpsolidm[1], c.n))
        b = _rock_heat_change(c, fill(Erebus.rhocpsolidm[1], c.n); rhocptotalm=z, ktotalm=z)
        @test isapprox(b.tk_sub, a.tk_sub; rtol=0.0, atol=0.0)
        @test all(isfinite, b.tk_sub)
    end

    @testset "zero time step and zero heat capacity change nothing" begin
        c = _subgrid_cell(; n_air=0)
        phi1 = vcat(1.0, zeros(c.n - 1))
        z = _rock_heat_change(c, fill(Erebus.rhocpsolidm[1], c.n); dt=0.0, phi_dry=phi1)
        @test all(isfinite, z.tk_sub)
        @test isapprox(z.tk_sub, c.tkm; rtol=0.0, atol=0.0)
        # A marker of porosity 1 holds no heat: it takes the grid temperature, no NaN.
        tkm = copy(c.tkm)
        Erebus.apply_subgrid_temperature_diffusion!(
            c.xm,
            c.ym,
            c.tm,
            tkm,
            zeros(c.n),
            c.tk1,
            zeros(size(c.tk1)),
            zeros(size(c.tk1)),
            zeros(size(c.tk1)),
            3.0e12,
            c.n,
            1;
            coords=c.coords,
            dsubgridt=1.0,
            phi_dry=vcat(1.0, zeros(c.n - 1)),
        )
        @test all(isfinite, tkm)
        @test isapprox(tkm[1], 300.0; rtol=1e-12)
        @test !isapprox(tkm[2:end], c.tkm[2:end]; rtol=1e-6)
    end

    @testset "error contract" begin
        c = _subgrid_cell()
        args(dt) = (
            c.xm,
            c.ym,
            c.tm,
            copy(c.tkm),
            zeros(c.n),
            c.tk1,
            zeros(size(c.tk1)),
            zeros(size(c.tk1)),
            zeros(size(c.tk1)),
            dt,
            c.n,
            1,
        )
        f = Erebus.apply_subgrid_temperature_diffusion!
        @test_throws DomainError f(args(1.0)...; coords=c.coords, dsubgridt=-1.0)
        @test_throws DomainError f(args(1.0)...; coords=c.coords, dsubgridt=NaN)
        @test_throws DomainError f(args(-1.0)...; coords=c.coords, dsubgridt=1.0)
        @test_throws DomainError f(args(Inf)...; coords=c.coords, dsubgridt=1.0)
    end

    @testset "the step passes the marker heat capacity and conductivity" begin
        cfg = override_config(
            default_config(), Dict{String,Any}("geometry.spherical_metric" => false)
        )
        state, coords, ws, _, _ = init_simulation(cfg)
        core = state.markers.core
        rocks = findall(<(3), core.tm)
        core.rhocptotalm[rocks] .= 2.0 * Erebus.rhocpsolidm[1]
        core.ktotalm[rocks] .= 0.25 * Erebus.ksolidm[1]
        core.tkm[rocks] .= 400.0
        state.grids.tk1 .= 300.0
        state.grids.tk2 .= 300.0
        state.grids.DT .= 0.0
        state.timestep = 2
        args = (core.xm, core.ym, core.tm)
        n = length(state.markers)
        function direct(; kw...)
            tkm = copy(core.tkm)
            g = state.grids
            Erebus.apply_subgrid_temperature_diffusion!(
                args...,
                tkm,
                core.phim,
                g.tk1,
                zeros(size(g.DT)),
                zeros(size(g.DT)),
                zeros(size(g.DT)),
                state.dt,
                n,
                1;
                coords=coords,
                dsubgridt=cfg.solver.dsubgridt,
                phi_dry=Erebus.dry_porosity(state.markers),
                kw...,
            )
            return tkm[rocks]
        end
        with_marker = direct(; rhocptotalm=core.rhocptotalm, ktotalm=core.ktotalm)
        without = direct()
        diffuse_and_update_markers!(state, coords, cfg, ws)
        @test isapprox(core.tkm[rocks], with_marker; rtol=1e-12)
        # Discrimination: the constitutive values give other temperatures.
        @test !isapprox(with_marker, without; rtol=1e-6)
    end

    @testset "the relaxation takes no rock marker below its cold floor" begin
        # Two rock markers in one cell between cold (290 K) and hot (1800 K) nodes, at a
        # large Fourier number: the counterpart of the hot-side marker must not push the
        # cold-side marker below 290 K.
        coords = GridCoordinates(GridConfig(; Nx=51, Ny=51, xsize=5.0e4, ysize=5.0e4))
        j0 = 25
        xL = coords.xp[j0]
        xm = [xL + 0.02 * coords.dx, xL + 0.48 * coords.dx]
        ym = fill(coords.yp[j0] + 0.5 * coords.dy, 2)
        tm = [1, 1]
        tkm = [290.0, 600.0]
        tk1 = fill(290.0, coords.Ny1, coords.Nx1)
        tk1[:, (j0 + 1):end] .= 1800.0
        rhocp = fill(3.3e6, 2)
        k = fill(3.0, 2)
        DT = zeros(size(tk1))
        TKSUM = zeros(size(tk1))
        RHOCPSUM = zeros(size(tk1))
        tkm_pre = copy(tkm)
        Erebus.apply_subgrid_temperature_diffusion!(
            xm,
            ym,
            tm,
            tkm,
            zeros(2),
            tk1,
            DT,
            TKSUM,
            RHOCPSUM,
            1.0e12,
            2,
            1;
            coords=coords,
            dsubgridt=1.0,
            rhocptotalm=rhocp,
            ktotalm=k,
        )
        Erebus.limit_relaxed_marker_temperatures!(
            xm,
            ym,
            tm,
            tkm,
            rhocp,
            DT,
            zeros(size(tk1)),
            tk1,
            tk1,
            TKSUM,
            RHOCPSUM,
            2;
            coords=coords,
            tkm_floor=tkm_pre,
        )
        Erebus.update_marker_temperature!(xm, ym, tkm, DT, tk1, 2, 2; coords=coords)
        @test all(isfinite, tkm)
        @test minimum(tkm) >= 290.0 - 1.0e-9
        # The relaxation only moves heat between the two markers: the rejected cooling of
        # the cold one goes back to the other, so the heat of the pair is kept.
        dE = sum(rhocp .* (tkm .- tkm_pre))
        @test abs(dE) <= 1.0e-12 * sum(rhocp .* tkm_pre)
        @test tkm[2] >= tkm[1]
    end

    @testset "each term of the cold floor can set the bound" begin
        coords = GridCoordinates(GridConfig(; Nx=10, Ny=10, xsize=1.0e5, ysize=1.0e5))
        sz = (coords.Ny1, coords.Nx1)
        xm, ym = [0.5 * coords.xsize + 10.0], [0.5 * coords.ysize + 10.0]
        function limited(tkm0, floor, tk1v, tk2v, dTc)
            tkm = [tkm0]
            DT = fill(-30.0, sz)
            Erebus.limit_relaxed_marker_temperatures!(
                xm,
                ym,
                [1],
                tkm,
                [3.3e6],
                DT,
                fill(dTc, sz),
                fill(tk1v, sz),
                fill(tk2v, sz),
                zeros(sz),
                zeros(sz),
                1;
                coords=coords,
                tkm_floor=[floor],
            )
            Erebus.update_marker_temperature!(
                xm, ym, tkm, DT, fill(tk2v, sz), 2, 1; coords=coords
            )
            return tkm[1]
        end
        # The start of the marker is the lowest term: 260 - 30 stops at 250.
        @test isapprox(limited(260.0, 250.0, 300.0, 300.0, 0.0), 250.0; rtol=1e-12)
        # The start node temperature is the lowest term: 280 - 30 stops at 270.
        @test isapprox(limited(280.0, 400.0, 270.0, 300.0, 0.0), 270.0; rtol=1e-12)
        # Conduction cooling stops at the end node temperature, 280.
        @test isapprox(limited(300.0, 250.0, 260.0, 280.0, -30.0), 280.0; rtol=1e-12)
        # Subgrid warming above the floor is kept under conduction cooling: 250 - 30 stops at 250.
        @test isapprox(limited(250.0, 200.0, 300.0, 280.0, -30.0), 250.0; rtol=1e-12)
        # Above the floor, the subgrid part is not limited: 260 - 30 = 230.
        @test isapprox(limited(260.0, 100.0, 100.0, 300.0, 0.0), 230.0; rtol=1e-12)
    end

    @testset "conduction acts on top of the cold floor" begin
        # Two markers between cold and hot nodes: the subgrid part stops at the floor, and the
        # conduction increment dTc of the nodes then acts in full on both markers.
        coords = GridCoordinates(GridConfig(; Nx=51, Ny=51, xsize=5.0e4, ysize=5.0e4))
        j0 = 25
        xm = [coords.xp[j0] + 0.02 * coords.dx, coords.xp[j0] + 0.48 * coords.dx]
        ym = fill(coords.yp[j0] + 0.5 * coords.dy, 2)
        rhocp, k = fill(3.3e6, 2), fill(3.0, 2)
        function pair(dTc)
            tkm = [290.0, 600.0]
            tkm_pre = copy(tkm)
            tk1 = fill(290.0, coords.Ny1, coords.Nx1)
            tk1[:, (j0 + 1):end] .= 1800.0
            DTc = fill(dTc, size(tk1))
            tk2 = tk1 .+ DTc
            DT = copy(DTc)
            TKSUM, RHOCPSUM = zeros(size(tk1)), zeros(size(tk1))
            Erebus.apply_subgrid_temperature_diffusion!(
                xm,
                ym,
                [1, 1],
                tkm,
                zeros(2),
                tk1,
                DT,
                TKSUM,
                RHOCPSUM,
                1.0e12,
                2,
                1;
                coords=coords,
                dsubgridt=1.0,
                rhocptotalm=rhocp,
                ktotalm=k,
            )
            Erebus.limit_relaxed_marker_temperatures!(
                xm,
                ym,
                [1, 1],
                tkm,
                rhocp,
                DT,
                DTc,
                tk1,
                tk2,
                TKSUM,
                RHOCPSUM,
                2;
                coords=coords,
                tkm_floor=tkm_pre,
            )
            Erebus.update_marker_temperature!(xm, ym, tkm, DT, tk2, 2, 2; coords=coords)
            return tkm
        end
        for (dTc, expected) in (
            (0.0, [290.0, 600.0]),
            (50.0, [340.0, 650.0]),
            (150.0, [440.0, 750.0]),
            (-20.0, [270.0, 580.0]),
        )
            t = pair(dTc)
            @test isapprox(t, expected; rtol=1e-12)
            # Heat of the pair: only the conduction increment of both markers is added.
            @test isapprox(sum(t) - 890.0, 2.0 * dTc; atol=1e-9)
        end
    end

    @testset "redistributed heat counts toward the cold floor" begin
        # Heat that a limited marker takes from its neighbours must not move them below their
        # own floor plus conduction, and must not create heat.
        coords = GridCoordinates(GridConfig(; Nx=51, Ny=51, xsize=5.0e4, ysize=5.0e4))
        j0 = 25
        function cell(tkm, k, offsets, dTc, dt)
            n = length(tkm)
            xm = coords.xp[j0] .+ offsets .* coords.dx
            ym = fill(coords.yp[j0] + 0.5 * coords.dy, n)
            rhocp = fill(3.3e6, n)
            tkm = copy(tkm)
            tkm_pre = copy(tkm)
            tk1 = fill(290.0, coords.Ny1, coords.Nx1)
            tk1[:, (j0 + 1):end] .= 1800.0
            DTc = fill(dTc, size(tk1))
            tk2 = tk1 .+ DTc
            DT = copy(DTc)
            TKSUM, RHOCPSUM = zeros(size(tk1)), zeros(size(tk1))
            Erebus.apply_subgrid_temperature_diffusion!(
                xm,
                ym,
                ones(Int, n),
                tkm,
                zeros(n),
                tk1,
                DT,
                TKSUM,
                RHOCPSUM,
                dt,
                n,
                1;
                coords=coords,
                dsubgridt=1.0,
                rhocptotalm=rhocp,
                ktotalm=k,
            )
            stats = Erebus.limit_relaxed_marker_temperatures!(
                xm,
                ym,
                ones(Int, n),
                tkm,
                rhocp,
                DT,
                DTc,
                tk1,
                tk2,
                TKSUM,
                RHOCPSUM,
                n;
                coords=coords,
                tkm_floor=tkm_pre,
            )
            Erebus.update_marker_temperature!(xm, ym, tkm, DT, tk2, 2, n; coords=coords)
            return tkm, stats
        end
        # A marker colder than its nodes gives the subgrid heat back down to its own floor.
        t, stats = cell([290.0, 200.0], [0.5, 5.0], [0.02, 0.48], 0.0, 1.0e11)
        @test isapprox(t, [290.0, 200.0]; rtol=1e-10)
        @test isapprox(stats.heat_lost, 0.0; atol=1.0)
        # Three markers: every marker keeps its conduction warming, and only that is added.
        for dTc in (50.0, 150.0)
            t, stats = cell(
                fill(290.0, 3), [0.5, 2.0, 5.0], [0.02, 0.31, 0.48], dTc, 1.0e12
            )
            @test minimum(t) >= 290.0 + dTc - 1e-9
            @test isapprox(sum(t), 3 * (290.0 + dTc); rtol=1e-12)
            @test isapprox(stats.heat_lost, 0.0; atol=1.0)
        end
    end

    @testset "the floor stage moves only subgrid heat" begin
        # Solve cooling on one side of the cell: the solve heat is kept in full, and a marker
        # far above its floor is not held up by rounding of the redistributed heat.
        coords = GridCoordinates(GridConfig(; Nx=21, Ny=21, xsize=2.0e4, ysize=2.0e4))
        I0, J0 = 10, 10
        function cell(tkm0, k, rhocp, ox, oy, tk1_4, dtc_4, dt)
            n = length(tkm0)
            xm = coords.xp[J0] .+ ox .* coords.dx
            ym = coords.yp[I0] .+ oy .* coords.dy
            tk1 = fill(tk1_4[1], coords.Ny1, coords.Nx1)
            DTc = zeros(coords.Ny1, coords.Nx1)
            for (q, (di, dj)) in enumerate(((0, 0), (1, 0), (0, 1), (1, 1)))
                tk1[I0 + di, J0 + dj] = tk1_4[q]
                DTc[I0 + di, J0 + dj] = dtc_4[q]
            end
            tk2 = tk1 .+ DTc
            tkm = copy(tkm0)
            DT = copy(DTc)
            TKSUM, RHOCPSUM = zeros(size(tk1)), zeros(size(tk1))
            Erebus.apply_subgrid_temperature_diffusion!(
                xm,
                ym,
                ones(Int, n),
                tkm,
                zeros(n),
                tk1,
                DT,
                TKSUM,
                RHOCPSUM,
                dt,
                n,
                1;
                coords=coords,
                dsubgridt=1.0,
                rhocptotalm=rhocp,
                ktotalm=k,
            )
            stats = Erebus.limit_relaxed_marker_temperatures!(
                xm,
                ym,
                ones(Int, n),
                tkm,
                rhocp,
                DT,
                DTc,
                tk1,
                tk2,
                TKSUM,
                RHOCPSUM,
                n;
                coords=coords,
                tkm_floor=copy(tkm0),
            )
            Erebus.update_marker_temperature!(xm, ym, tkm, DT, tk2, 2, n; coords=coords)
            return tkm, stats
        end
        t, stats = cell(
            [300.0, 302.0],
            [0.01, 5.0],
            [3.3e6, 6.6e6],
            [0.02, 0.98],
            [0.5, 0.5],
            [400.0, 400.0, 300.0, 300.0],
            [-150.0, -150.0, 0.0, 0.0],
            1.0e10,
        )
        @test isapprox(t, [250.0, 250.5]; rtol=1e-12)
        @test isapprox(stats.heat_lost, 0.0; atol=1.0)
        rhocp = [1.127898595768074e6, 1.1278189357731966e6, 1.4882315746412517e6]
        tkm0 = [310.0558627179961, 161.6159635491205, 301.8305152622005]
        t, stats = cell(
            tkm0,
            [0.34638470927542087, 0.03583300936468784, 0.007741122651759497],
            rhocp,
            [0.8497201644381142, 0.0050136960617534814, 0.817269692998967],
            [0.08403734228921009, 0.7646781607754731, 0.3598581006977638],
            [1215.9853824485313, 1215.9853824485313, 301.8305152622005, 301.8305152622005],
            zeros(4),
            4.934862247478527e12,
        )
        # No solve increment: the heat of the cell does not change.
        @test abs(sum(rhocp .* (t .- tkm0))) <= 1.0e-12 * sum(rhocp .* tkm0)
        @test isapprox(stats.heat_lost, 0.0; atol=1.0)
    end

    @testset "the cold floor needs both inputs of the right size" begin
        coords = GridCoordinates(GridConfig(; Nx=10, Ny=10, xsize=1.0e5, ysize=1.0e5))
        sz = (coords.Ny1, coords.Nx1)
        f(; kw...) = Erebus.limit_marker_temperature_increment!(
            [5.0e4],
            [5.0e4],
            [1],
            [300.0],
            [3.3e6],
            zeros(sz),
            zeros(sz),
            zeros(sz),
            zeros(sz),
            1;
            coords=coords,
            kw...,
        )
        @test_throws ArgumentError f(; tkm_floor=[300.0])
        @test_throws ArgumentError f(; tk_floor=zeros(sz))
        @test_throws ArgumentError f(; tkm_floor=Float64[], tk_floor=zeros(sz))
        @test_throws ArgumentError f(; tkm_floor=[300.0], tk_floor=zeros(2, 2))
        @test_throws ArgumentError f(;
            tkm_floor=[300.0], tk_floor=zeros(sz), DT_conduction=zeros(sz)
        )
        @test_throws ArgumentError Erebus.limit_relaxed_marker_temperatures!(
            [5.0e4],
            [5.0e4],
            [1],
            [300.0],
            [3.3e6],
            zeros(sz),
            zeros(2, 2),
            zeros(sz),
            zeros(sz),
            zeros(sz),
            zeros(sz),
            1;
            coords=coords,
            tkm_floor=[300.0],
        )
    end

    @testset "the step keeps rock markers above the cold floor" begin
        cfg = override_config(
            default_config(), Dict{String,Any}("geometry.spherical_metric" => false)
        )
        state, coords, ws, _, _ = init_simulation(cfg)
        core = state.markers.core
        g = state.grids
        rocks = findall(<(3), core.tm)
        core.rhocptotalm[rocks] .= 3.3e6
        core.ktotalm[rocks] .= 3.0
        core.tkm[rocks] .= 290.0
        # Hot and cold node columns alternate, so many cells hold a steep gradient.
        g.tk1 .= 290.0
        g.tk1[:, 2:2:end] .= 1800.0
        g.tk2 .= g.tk1
        g.DT .= 0.0
        g.DT0 .= 0.0
        state.dt = 1.0e13
        state.timestep = 2
        tkm_pre = copy(core.tkm)
        diffuse_and_update_markers!(state, coords, cfg, ws)
        Erebus.update_marker_temperature!(
            core.xm, core.ym, core.tkm, g.DT, g.tk2, 2, length(state.markers); coords=coords
        )
        @test all(isfinite, core.tkm[rocks])
        @test minimum(core.tkm[rocks]) >= 290.0 - 1.0e-9
        # The relaxation acted: some rock markers warmed toward the hot nodes.
        @test !isapprox(core.tkm[rocks], tkm_pre[rocks]; rtol=1e-6)
    end

    @testset "the step helper equals the explicit calls" begin
        cfg = override_config(
            default_config(), Dict{String,Any}("geometry.spherical_metric" => false)
        )
        state, coords, ws, _, _ = init_simulation(cfg)
        core = state.markers.core
        g = state.grids
        rocks = findall(<(3), core.tm)
        core.rhocptotalm[rocks] .= 3.3e6
        core.ktotalm[rocks] .= range(0.5, 5.0; length=length(rocks))
        core.tkm[rocks] .= 290.0
        # Markers start below all nodes, so the floor is the marker start (290 K), not a node.
        g.tk1 .= 500.0
        g.tk1[:, 2:2:end] .= 1800.0
        g.tk2 .= g.tk1 .+ 5.0
        g.DT .= 5.0
        state.dt = 1.0e13
        state.timestep = 2
        n = length(state.markers)
        tkm = copy(core.tkm)
        DT = copy(g.DT)
        TK, RC = zeros(size(DT)), zeros(size(DT))
        Erebus.apply_subgrid_temperature_diffusion!(
            core.xm,
            core.ym,
            core.tm,
            tkm,
            core.phim,
            g.tk1,
            DT,
            TK,
            RC,
            state.dt,
            n,
            Erebus.marker_property_mode;
            coords=coords,
            dsubgridt=1.0,
            phi_dry=Erebus.dry_porosity(state.markers),
            rhocptotalm=core.rhocptotalm,
            ktotalm=core.ktotalm,
        )
        Erebus.limit_relaxed_marker_temperatures!(
            core.xm,
            core.ym,
            core.tm,
            tkm,
            core.rhocptotalm,
            DT,
            copy(g.DT),
            g.tk1,
            g.tk2,
            TK,
            RC,
            n;
            coords=coords,
            tkm_floor=copy(core.tkm),
        )
        diffuse_and_update_markers!(state, coords, cfg, ws)
        @test isapprox(core.tkm[rocks], tkm[rocks]; rtol=1e-12)
        @test isapprox(g.DT, DT; rtol=1e-12)
        # A zero time step leaves the marker temperatures unchanged.
        tkm0 = copy(core.tkm)
        state.dt = 0.0
        g.DT .= 0.0
        diffuse_and_update_markers!(state, coords, cfg, ws)
        @test isapprox(core.tkm[rocks], tkm0[rocks]; rtol=0.0, atol=0.0)
    end

    @testset "rock marker temperature follows the grid at the air boundary" begin
        P = PorousShellValidation
        Nx = 65
        cfg = P.make_porous_shell_config(;
            Nx, phi_init=0.25, t_form_yr=2.5e6, n_steps=20, dt_yr=5.0e4
        )
        st = mktempdir() do d
            redirect_stdout(devnull) do
                with_logger(NullLogger()) do
                    simulation_loop(cfg; output_path=d)
                end
            end
        end
        mk = st.markers
        co = GridCoordinates(cfg.grid)
        xc = yc = P.L_BOX / 2
        drift = Float64[]
        for m in eachindex(mk.core.tm)
            mk.core.tm[m] < 3 || continue
            hypot(mk.core.xm[m] - xc, mk.core.ym[m] - yc) <= P.R_PLANET || continue
            i, j, w = Erebus.fix_weights(
                mk.core.xm[m],
                mk.core.ym[m],
                co.xp,
                co.yp,
                co.dx,
                co.dy,
                co.jmin_p,
                co.jmax_p,
                co.imin_p,
                co.imax_p,
            )
            Tg = Erebus.dot4(Erebus.grid_vector(i, j, st.grids.tk2), w)
            push!(drift, mk.core.tkm[m] - Tg)
        end
        @test !isempty(drift)
        @test all(isfinite, drift)
        @test maximum(abs, drift) <= DRIFT_BOUND_K
        # Bound: the heat exchange of the relaxation takes no rock marker below the
        # space temperature, the lowest temperature of the case.
        rock = findall(<(3), mk.core.tm)
        @test minimum(mk.core.tkm[rock]) >= P.T_SPACE - 1.0e-9
    end
end

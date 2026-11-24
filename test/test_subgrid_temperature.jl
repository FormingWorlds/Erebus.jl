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
function _rock_heat_change(c, rhocp_rock; dt=3.0e12, kw...)
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
        phi_dry=zeros(c.n),
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
        @test isapprox(SolverConfig().dsubgridt, 1.0; rtol=1e-12)
        # Scale guard: a factor of 10 or 0.1 is not the default.
        @test !isapprox(SolverConfig().dsubgridt, 10.0; rtol=1e-3)
        @test !isapprox(SolverConfig().dsubgridt, 0.1; rtol=1e-3)
    end

    @testset "rock heat is kept when sticky air shares the nodes" begin
        c = _subgrid_cell()
        rhocp_rock = fill(Erebus.rhocpsolidm[1], c.n)
        r = _rock_heat_change(c, rhocp_rock)
        # The rock markers relax toward 300 K, so the scale is not zero.
        @test r.scale > 1.0e6
        @test abs(r.dE) <= 1.0e-12 * r.scale
        # Sticky-air markers keep their temperature in the subgrid step.
        air = findall(>=(3), c.tm)
        @test r.tk_sub[air] == c.tkm[air]
        # Bounds: the relaxed rock markers stay between the node and their start values.
        rock = findall(<(3), c.tm)
        @test all(300.0 .<= r.tk_sub[rock] .<= c.tkm[rock])
    end

    @testset "relaxation uses the marker heat capacity and conductivity" begin
        c = _subgrid_cell(; n_air=0)
        dt = 3.0e12
        rhocp_m = fill(2.0 * Erebus.rhocpsolidm[1], c.n)
        k_m = fill(0.5 * Erebus.ksolidm[1], c.n)
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
        @test b.tk_sub == a.tk_sub
        @test all(isfinite, b.tk_sub)
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
    end
end

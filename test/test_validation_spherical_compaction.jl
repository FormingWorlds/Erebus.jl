using SHA
using Test
using TOML
using Erebus

isdefined(@__MODULE__, :SphericalCompactionValidation) ||
    include(joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_compaction.jl"))
include(joinpath(@__DIR__, "..", "tools", "check_validation.jl"))
isdefined(@__MODULE__, :SphericalCompactionExport) ||
    include(joinpath(@__DIR__, "..", "tools", "export_validation_spherical_compaction.jl"))

const SCV = SphericalCompactionValidation

"""
    InterruptingBench

Benchmark stand-in whose every property access throws `InterruptException`, as a Ctrl-C
during the benchmark does.
"""
struct InterruptingBench end
Base.getproperty(::InterruptingBench, ::Symbol) = throw(InterruptException())

"""
    CrashingBench

Benchmark stand-in whose every property access throws an `ErrorException`, as a solver
failure during the benchmark does.
"""
struct CrashingBench end
Base.getproperty(::CrashingBench, ::Symbol) = error("solver failure inside the benchmark")

"""
    write_compaction_toml(path, dict)

Write `dict` as TOML to `path` and return `path`.

Parameters
----------
path : String
    Output file path.
dict : AbstractDict
    Table to write.

Returns
-------
String
    `path`.
"""
function write_compaction_toml(path::String, dict::AbstractDict)
    open(path, "w") do io
        return TOML.print(io, dict)
    end
    return path
end

"""
    capture_stdio(f)

Call `f()` with standard output and standard error redirected to files and return the
result with the captured text.

Parameters
----------
f : Function
    Function of no arguments.

Returns
-------
NamedTuple
    `value`: the return value of `f`; `out`, `err`: the captured standard output and
    standard error.
"""
function capture_stdio(f)
    return mktempdir() do dir
        out_path, err_path = joinpath(dir, "out.txt"), joinpath(dir, "err.txt")
        value = open(out_path, "w") do out
            open(err_path, "w") do err
                return redirect_stdio(f; stdout=out, stderr=err)
            end
        end
        return (; value, out=read(out_path, String), err=read(err_path, String))
    end
end

"""
    surface_rock_cells(coords)

Return the P cells `(i, j)` inside the sphere with at least one of their 4 neighbours
outside it, over the rows and columns that `_metrics` reads.

Parameters
----------
coords : NamedTuple
    Coordinates from `compaction_coords`.

Returns
-------
Vector{Tuple{Int,Int}}
    Indices of the surface rock cells.
"""
function surface_rock_cells(coords)
    inside(i, j) =
        hypot(coords.xp[j] - SCV.X_CENTRE, coords.yp[i] - SCV.Y_CENTRE) <= SCV.R_PLANET
    cells = Tuple{Int,Int}[]
    for j in 3:(coords.Nx1 - 2), i in 3:(coords.Ny1 - 2)
        inside(i, j) || continue
        neighbours = ((i - 1, j), (i + 1, j), (i, j - 1), (i, j + 1))
        any(!inside(a, b) for (a, b) in neighbours) && push!(cells, (i, j))
    end
    return cells
end

"""
    compaction_coords(N; shift=0.0)

Return a `NamedTuple` with the staggered coordinates that the error measures read, laid out
as `GridCoordinates` does for an `N` x `N` grid of the box `SCV.L_BOX`, with every coordinate
moved by `shift` [m] in x and y.

Parameters
----------
N : Int
    Basic grid resolution in x and y.
shift : Float64, optional
    Offset of all coordinates [m] (default `0.0`).

Returns
-------
NamedTuple
    `Nx`, `Ny`, `Nx1`, `Ny1`, `dx`, `dy`, `xvx`, `yvx`, `xvy`, `yvy`, `xp` and `yp`.
"""
function compaction_coords(N::Int; shift::Float64=0.0)
    d = SCV.L_BOX / (N - 1)
    node(lo, hi) = collect(range(lo, hi; length=N + 1)) .+ shift
    return (;
        Nx=N,
        Ny=N,
        Nx1=N + 1,
        Ny1=N + 1,
        dx=d,
        dy=d,
        xvx=node(0.0, SCV.L_BOX + d),
        yvx=node(-d / 2, SCV.L_BOX + d / 2),
        xvy=node(-d / 2, SCV.L_BOX + d / 2),
        yvy=node(0.0, SCV.L_BOX + d),
        xp=node(-d / 2, SCV.L_BOX + d / 2),
        yp=node(-d / 2, SCV.L_BOX + d / 2),
    )
end

"""
    exact_solid_velocity_grids(coords, d, scale)

Return the grids `vx` and `vy` as a `NamedTuple` holding `scale` times the exact solid
velocity `-(k / mu) Pe'(r) r_hat` [m/s] at every face of `coords` with `r <= R_PLANET`, and
zero outside.

Parameters
----------
coords : NamedTuple
    Coordinates from `compaction_coords`.
d : Int
    Flow dimension, 3 for the sphere and 2 for the disk.
scale : Float64
    Factor applied to the exact velocity.

Returns
-------
NamedTuple
    `vx` and `vy`, each of size `(coords.Ny1, coords.Nx1)`.
"""
function exact_solid_velocity_grids(coords, d::Int, scale::Float64)
    function comp(x, y, c)
        r = hypot(x, y)
        (r <= SCV.R_PLANET && r > 0) || return 0.0
        return -scale * SCV.K_PERM / SCV.MU_F * SCV.exact_pe_derivative(r, d) * c / r
    end
    vx = [
        comp(x - SCV.X_CENTRE, y - SCV.Y_CENTRE, x - SCV.X_CENTRE) for
        y in coords.yvx, x in coords.xvx
    ]
    vy = [
        comp(x - SCV.X_CENTRE, y - SCV.Y_CENTRE, y - SCV.Y_CENTRE) for
        y in coords.yvy, x in coords.xvy
    ]
    return (; vx, vy)
end

"""
    trapezoid_bessel(x, n)

Return the modified Bessel function `I_n(x)` for `n = 0` or `1` from the integral form
`(1 / pi) int_0^pi exp(x cos t) cos(n t) dt`, with the trapezoidal rule on 400 intervals.
The integrand is smooth and periodic, so the rule converges faster than any power.

Parameters
----------
x : Float64
    Argument, dimensionless.
n : Int
    Order, 0 or 1.

Returns
-------
Float64
    `I_n(x)`.
"""
function trapezoid_bessel(x::Float64, n::Int)
    m = 400
    h = pi / m
    f(t) = exp(x * cos(t)) * cos(n * t)
    s = 0.5 * (f(0.0) + f(pi))
    for k in 1:(m - 1)
        s += f(k * h)
    end
    return s * h / pi
end

@testset "Spherical Compaction Validation" begin
    @testset "modified Bessel functions against the integral form" begin
        @test isapprox(SCV._i0(2.0), 2.2795853023360673; rtol=1.0e-14)
        @test isapprox(SCV._i0(0.0), 1.0; rtol=1.0e-15)
        @test isapprox(SCV._i1(0.0), 0.0; atol=1.0e-300)
        for x in (0.5, 2.0, 5.0, 20.0)
            @test isapprox(SCV._i0(x), trapezoid_bessel(x, 0); rtol=1.0e-12)
            @test isapprox(SCV._i1(x), trapezoid_bessel(x, 1); rtol=1.0e-12)
            # I0' = I1: a wrong series index or factor breaks this.
            h = 1.0e-5 * max(x, 1.0)
            fd = (SCV._i0(x + h) - SCV._i0(x - h)) / (2h)
            @test isapprox(fd, SCV._i1(x); rtol=1.0e-8)
        end
        for f in (SCV._i0, SCV._i1)
            @test_throws DomainError f(-1.0)
            @test_throws DomainError f(NaN)
            @test_throws DomainError f(Inf)
            @test_throws DomainError f(SCV.BESSEL_X_MAX * 1.01)
        end
        @test_throws DomainError SCV._bessel_series(2.0, -1)
        @test_throws DomainError SCV._bessel_series(2.0, 2)
    end

    @testset "Bessel functions and exact profile sum in Float64 for any real input" begin
        # Float32 and Float16 arguments that are exact in Float64 give the Float64 value.
        for x in (100.0f0, 700.0f0, Float16(20.0), 2, BigFloat(2.0))
            for f in (SCV._i0, SCV._i1)
                v = f(x)
                @test v isa Float64
                @test isapprox(v, f(Float64(x)); rtol=1.0e-15)
            end
        end
        @test isapprox(SCV._i1(100.0f0), trapezoid_bessel(100.0, 1); rtol=1.0e-12)
        r_big = BigFloat(1.0e4)
        for d in (2, 3)
            @test SCV.exact_pe(r_big, d) isa Float64
            @test SCV.exact_pe_derivative(r_big, d) isa Float64
            @test isapprox(
                SCV.exact_pe(Float32(1.0e4), d), SCV.exact_pe(1.0e4, d); rtol=1.0e-15
            )
        end
    end

    @testset "exact effective pressure: surface value, centre limit, Helmholtz equation" begin
        R, delta = SCV.R_PLANET, SCV.DELTA
        @test isapprox(delta, R / 2; rtol=1.0e-15)
        for d in (2, 3)
            @test isapprox(SCV.exact_pe(R, d), SCV.PE_R; rtol=1.0e-14)
            # Pe grows from the centre to the surface.
            rs = range(0.0, R; length=101)
            @test all(diff([SCV.exact_pe(r, d) for r in rs]) .> 0)
            # delta^2 (Pe'' + (d - 1) Pe' / r) = Pe in d dimensions.
            for r in (0.2R, 0.5R, 0.8R)
                h = 1.0e-3 * R
                p0, pp, pm = SCV.exact_pe(r, d),
                SCV.exact_pe(r + h, d),
                SCV.exact_pe(r - h, d)
                lap = (pp - 2p0 + pm) / h^2 + (d - 1) / r * (pp - pm) / (2h)
                @test isapprox(delta^2 * lap, p0; rtol=1.0e-5)
            end
        end
        # Centre values: Pe_R (R / delta) / sinh(R / delta) and Pe_R / I0(R / delta).
        pc3 = SCV.exact_pe(0.0, 3)
        pc2 = SCV.exact_pe(0.0, 2)
        @test isapprox(pc3, SCV.PE_R * 2 / sinh(2.0); rtol=1.0e-14)
        @test isapprox(pc2, SCV.PE_R / 2.2795853023360673; rtol=1.0e-14)
        @test isapprox(SCV.exact_pe(1.0e-6 * R, 3), pc3; rtol=1.0e-12)
        @test isapprox(
            SCV.exact_pe(SCV.SERIES_SWITCH * delta * (1 - 1.0e-9), 3), pc3; rtol=1.0e-4
        )
        # Sphere to disk centre ratio 2 I0(2) / sinh(2), neither 1 nor the Darcy value 2/3.
        @test isapprox(pc3 / pc2, SCV.CENTRE_RATIO_EXACT; rtol=1.0e-14)
        @test isapprox(SCV.CENTRE_RATIO_EXACT, 1.257057; rtol=1.0e-6)
        @test !isapprox(SCV.CENTRE_RATIO_EXACT, 1.0; rtol=0.1)
        # Scale guard: the centre holds about half of the surface value.
        @test 0.4 * SCV.PE_R < pc2 < pc3 < 0.6 * SCV.PE_R
        @test_throws DomainError SCV.exact_pe(-1.0, 3)
        @test_throws DomainError SCV.exact_pe(R * 1.01, 3)
        @test_throws DomainError SCV.exact_pe(0.5R, 1)
        @test_throws DomainError SCV.exact_pe(NaN, 2)
    end

    @testset "exact derivative of the effective pressure and its small-r series" begin
        R, delta = SCV.R_PLANET, SCV.DELTA
        for d in (2, 3)
            @test isapprox(SCV.exact_pe_derivative(0.0, d), 0.0; atol=1.0e-300)
            for r in (1.0e-3R, 0.1R, 0.5R, 0.9R, R)
                h = 1.0e-4 * R
                lo = max(r - h, 0.0)
                hi = min(r + h, R)
                fd = (SCV.exact_pe(hi, d) - SCV.exact_pe(lo, d)) / (hi - lo)
                @test isapprox(SCV.exact_pe_derivative(r, d), fd; rtol=1.0e-4)
            end
        end
        # The series and the closed form agree on both sides of the switch.
        r_sw = SCV.SERIES_SWITCH * delta
        below = SCV.exact_pe_derivative(r_sw * (1 - 1.0e-12), 3)
        above = SCV.exact_pe_derivative(r_sw * (1 + 1.0e-12), 3)
        @test isapprox(below, above; rtol=1.0e-9)
        # Near the centre Pe' = Pe_c r / (3 delta^2) for the sphere, Pe_c r / (2 delta^2) for the disk.
        r_small = 1.0e-6 * delta
        for d in (2, 3)
            slope = SCV.exact_pe(0.0, d) * r_small / (d * delta^2)
            @test isapprox(SCV.exact_pe_derivative(r_small, d), slope; rtol=1.0e-9)
        end
        @test_throws DomainError SCV.exact_pe_derivative(-1.0, 3)
        @test_throws DomainError SCV.exact_pe_derivative(R * 1.01, 2)
        @test_throws DomainError SCV.exact_pe_derivative(0.5R, 4)
    end

    @testset "case set-up: compaction viscosity, vent pressure and input contracts" begin
        @test isapprox(SCV.ZETA, 7.8125e18; rtol=1.0e-14)
        # delta^2 = k zeta (1 - phi) / mu holds with the case zeta.
        @test isapprox(
            SCV.K_PERM * SCV.ZETA * (1 - SCV.PHI) / SCV.MU_F, SCV.DELTA^2; rtol=1.0e-14
        )
        cfg = SCV.make_config(SCV.NX_REFERENCE; spherical_flow=true)
        @test isapprox(cfg.solver.etaphikoef, SCV.PHI / SCV.ETA_RATIO; rtol=1.0e-14)
        @test isapprox(cfg.geometry.psurface - SCV.P_AMB, SCV.PE_R; rtol=1.0e-14)
        p_vent = Erebus.compute_venting_pressure(
            SCV.T_SURF, SCV.P_AMB; species=:H2O, L_sub=cfg.venting.L_sublimation
        )
        @test isapprox(p_vent, SCV.P_AMB; rtol=1.0e-12)
        @test Erebus.compute_water_vapor_pressure(SCV.T_SURF) < SCV.P_AMB
        @test_throws ArgumentError SCV.make_config(SCV.NX_MIN - 1; spherical_flow=true)
        @test_throws ArgumentError SCV.make_config(33; spherical_flow=true, eta_ratio=0.0)
        @test_throws ArgumentError SCV.make_config(33; spherical_flow=true, eta_ratio=1.0)
        @test_throws ArgumentError SCV.make_config(33; spherical_flow=true, dt_yr=-1.0)
        @test_throws ArgumentError SCV.run_case(SCV.NX_MIN - 1; spherical_flow=false)
        @test_throws ArgumentError SCV.evaluate_spherical_compaction_benchmark(;
            nx_list=(17,)
        )
        @test_throws ArgumentError SCV.evaluate_spherical_compaction_benchmark(;
            nx_list=(33.0,)
        )
    end

    record_path = joinpath(
        @__DIR__, "..", "benchmarks", "records", "spherical_compaction.toml"
    )
    raw_record_all = TOML.parsefile(record_path)
    recorded = Dict(
        String(m["id"]) => Float64(m["value"]) for m in raw_record_all["metrics"]
    )
    # Runs at Nx = 33 and 65 only. The Nx = 129 and order metrics take their recorded values,
    # which `validate_record_dict` checks against their reference and tolerance.
    bench_runs = SCV.evaluate_spherical_compaction_benchmark(; nx_list=(33, 65))
    runs = bench_runs.runs
    measured = bench_runs.metrics
    n_fine = maximum(SCV.NX_LIST)
    recorded_only = filter(
        kv -> !haskey(measured, kv.first) && endswith(kv.first, string(n_fine)), recorded
    )
    bench = (; runs, metrics=merge(measured, recorded_only))

    @testset "metrics: 17 from the runs, 13 recorded values for Nx = 129 and the orders" begin
        @test length(measured) == 17
        @test length(recorded_only) == 13
        @test isempty(intersect(keys(measured), keys(recorded_only)))
        @test Set(keys(bench.metrics)) == Set(SCV.EXPECTED_METRIC_IDS)
    end

    @testset "solver converges and the fields are finite (spherical_flow = $sf)" for sf in (
        false, true
    )
        for n in (33, 65)
            r = runs[(n, sf)]
            @test r.converged
            for field in (:pe, :pf, :vx, :vy, :qxD, :qyD)
                @test all(isfinite, getfield(r, field))
            end
        end
    end

    @testset "ETAPHI equals zeta on every rock P cell (spherical_flow = $sf)" for sf in (
        false, true
    )
        for n in (33, 65)
            @test runs[(n, sf)].n_bad_etaphi == 0
            @test runs[(n, sf)].n_rock == count(
                hypot(x - SCV.X_CENTRE, y - SCV.Y_CENTRE) <= SCV.R_PLANET for
                y in runs[(n, sf)].yp[3:(end - 2)], x in runs[(n, sf)].xp[3:(end - 2)]
            )
        end
    end

    @testset "total pressure stays uniform over the rock (spherical_flow = $sf)" for sf in (
        false, true
    )
        # Measured 8.0e-5 and 8.7e-5 of Pe_R at Nx = 33 and 65.
        @test runs[(33, sf)].pt_spread < 3.0e-4
        @test runs[(65, sf)].pt_spread < 3.0e-4
    end

    @testset "vent cells hold P_vent and the vented volume equals the compaction volume (spherical_flow = $sf)" for sf in
                                                                                                                    (
        false, true
    )
        for n in (33, 65)
            m = runs[(n, sf)]
            @test m.pf_vent_dev < 1.0e-6
            @test m.vent_closure < 1.0e-8
            @test isapprox(m.vent_flux, m.comp_flux; rtol=1.0e-8)
        end
    end

    @testset "centre effective pressure against the exact profile (spherical_flow = $sf)" for sf in
                                                                                              (
        false, true
    )
        c33 = runs[(33, sf)].centre_err
        c65 = runs[(65, sf)].centre_err
        # The staircase surface lies inside R, so the centre value is too high.
        @test 0.0 < c65 < 0.035
        @test 0.0 < c33 < 0.11
        @test isapprox(log2(c33 / c65), 1.0; atol=0.35)
    end

    @testset "effective pressure profile and solid velocity (spherical_flow = $sf)" for sf in
                                                                                        (
        false, true
    )
        a, b = runs[(33, sf)], runs[(65, sf)]
        @test isapprox(b.l2_inner, 0.0; atol=0.045)
        @test isapprox(b.vs_l2_inner, 0.0; atol=0.045)
        @test isapprox(log2(a.l2_inner / b.l2_inner), 1.0; atol=0.35)
        @test isapprox(log2(a.vs_l2_inner / b.vs_l2_inner), 1.0; atol=0.35)
    end

    @testset "sphere to disk centre ratio is 2 I0(2) / sinh(2)" begin
        on, off = runs[(65, true)], runs[(65, false)]
        ratio = on.pe_centre / off.pe_centre
        @test isapprox(ratio, SCV.CENTRE_RATIO_EXACT; rtol=0.012)
        @test !isapprox(ratio, 1.0; rtol=0.1)
        @test isapprox(measured["centre_ratio_on_off_nx65"], ratio; rtol=1.0e-14)
    end

    @testset "model error: eta / zeta = 1e-3 moves the centre error by less than 10 %" begin
        r3 = SCV.run_case(65; spherical_flow=true, eta_ratio=1.0e-3)
        r4 = runs[(65, true)]
        @test r3.converged
        @test abs(r3.centre_err - r4.centre_err) < 0.1 * r4.centre_err
        # The total pressure spread grows with eta / zeta.
        @test r3.pt_spread < 3.0e-3
        @test r3.pt_spread > 3 * r4.pt_spread
    end

    @testset "a long time step reaches the vent drain cap and breaks the closure (spherical_flow = $sf)" for (
        sf, factor
    ) in (
        (true, 114.2), (false, 108.9)
    )
        long = SCV.run_case(33; spherical_flow=sf, dt_yr=1.0e-2)
        @test long.converged
        @test !isapprox(long.vent_closure, 0.0; atol=1.0e-3)
        @test runs[(33, sf)].vent_closure < 1.0e-8
        # The matrix holds the cap from the pressure before the solve; S_vent_grid is
        # recomputed from the solved pressure without it, so it exceeds the compaction rate.
        @test isapprox(long.vent_flux / long.comp_flux, factor; rtol=0.01)
    end

    @testset "6-variable path equals the 4-variable path (spherical_flow = $sf)" for sf in (
        false, true
    )
        four = runs[(33, sf)]
        six = SCV.run_case(33; spherical_flow=sf, darcy_elimination=false)
        rel_diff(field) =
            maximum(abs, getfield(six, field) .- getfield(four, field)) /
            maximum(abs, getfield(four, field))
        @test six.converged
        @test rel_diff(:pf) < 1.0e-10
        @test rel_diff(:pe) < 1.0e-8
        # Every face, also in the sticky air, where eta_air / zeta = 1.3e-7 turns the
        # round-off of the total pressure into velocity: measured 3.7e-7 at Nx = 33.
        @test rel_diff(:vx) < 2.0e-5
        @test rel_diff(:vy) < 2.0e-5
        # The faces with r <= R: measured 6.7e-8 at Nx = 33.
        coords = compaction_coords(33)
        @test isapprox(coords.xp, four.xp; rtol=1.0e-14)
        rock(xs, ys) =
            [hypot(x - SCV.X_CENTRE, y - SCV.Y_CENTRE) <= SCV.R_PLANET for y in ys, x in xs]
        rock_diff(field, mask) =
            maximum(abs, (getfield(six, field) .- getfield(four, field))[mask]) /
            maximum(abs, getfield(four, field))
        @test rock_diff(:vx, rock(coords.xvx, coords.yvx)) < 2.0e-6
        @test rock_diff(:vy, rock(coords.xvy, coords.yvy)) < 2.0e-6
        @test rel_diff(:qxD) < 1.0e-10
        @test rel_diff(:qyD) < 1.0e-10
        @test isapprox(six.vs_l2_inner, four.vs_l2_inner; rtol=1.0e-6)
    end

    @testset "solid velocity error is 0 for the exact field and scales with the departure" begin
        coords = compaction_coords(41)
        for d in (2, 3)
            err(scale) =
                SCV._vs_l2_inner(exact_solid_velocity_grids(coords, d, scale), coords, d)
            @test isapprox(err(1.0), 0.0; atol=1.0e-12)
            @test isapprox(err(1.1), 0.1; atol=1.0e-12)
            @test isapprox(err(-1.0), 2.0; atol=1.0e-12)
        end
    end

    @testset "solid velocity error counts each component and only faces within 0.9 R" begin
        coords = compaction_coords(41)
        rx(i, j) = hypot(coords.xvx[j] - SCV.X_CENTRE, coords.yvx[i] - SCV.Y_CENTRE)
        ry(i, j) = hypot(coords.xvy[j] - SCV.X_CENTRE, coords.yvy[i] - SCV.Y_CENTRE)
        r_in = SCV.INNER_FRACTION * SCV.R_PLANET
        for d in (2, 3)
            exact = exact_solid_velocity_grids(coords, d, 1.0)
            # The grid is symmetric under x <-> y, so each component holds half of the
            # norm: flipping one component gives sqrt(4 / 2).
            @test isapprox(
                SCV._vs_l2_inner((; vx=(-exact.vx), vy=exact.vy), coords, d),
                sqrt(2.0);
                rtol=1.0e-12,
            )
            @test isapprox(
                SCV._vs_l2_inner((; vx=exact.vx, vy=(-exact.vy)), coords, d),
                sqrt(2.0);
                rtol=1.0e-12,
            )
            # Faces with 0.9 R < r do not enter; the faces within 0.9 R do.
            outer_x = [
                r_in < rx(i, j) ? 5.0 : 1.0 for
                i in axes(exact.vx, 1), j in axes(exact.vx, 2)
            ]
            outer_y = [
                r_in < ry(i, j) ? 5.0 : 1.0 for
                i in axes(exact.vy, 1), j in axes(exact.vy, 2)
            ]
            @test any(
                r_in < rx(i, j) <= SCV.R_PLANET for
                i in axes(exact.vx, 1), j in axes(exact.vx, 2)
            )
            @test isapprox(
                SCV._vs_l2_inner(
                    (; vx=exact.vx .* outer_x, vy=exact.vy .* outer_y), coords, d
                ),
                0.0;
                atol=1.0e-12,
            )
            inner_x = [
                rx(i, j) <= r_in ? 5.0 : 1.0 for
                i in axes(exact.vx, 1), j in axes(exact.vx, 2)
            ]
            @test isapprox(
                SCV._vs_l2_inner((; vx=exact.vx .* inner_x, vy=exact.vy), coords, d),
                sqrt(16 / 2);
                rtol=1.0e-12,
            )
        end
    end

    @testset "vent pressure deviation covers every surface rock cell" begin
        coords = compaction_coords(41)
        cells = surface_rock_cells(coords)
        n = coords.Ny1
        grids(pf, s_vent) = (;
            pr=fill(SCV.P_AMB + SCV.PE_R, n, n),
            pf,
            ETAPHI=fill(SCV.ZETA, n, n),
            PHI=fill(SCV.PHI, n, n),
            S_vent_grid=s_vent,
            vx=zeros(n, n),
            vy=zeros(n, n),
            qxD=zeros(n, n),
            qyD=zeros(n, n),
        )
        pf = fill(SCV.P_AMB + SCV.PE_R, n, n)
        s_vent = zeros(n, n)
        for (i, j) in cells
            pf[i, j] = SCV.P_AMB
            s_vent[i, j] = 1.0e-12
        end
        @test length(cells) > 50
        @test sort([
            (i, j) for j in 3:(coords.Nx1 - 2), i in 3:(coords.Ny1 - 2) if
            SCV._surface_rock(coords, i, j)
        ]) == sort(cells)
        m_all = SCV._metrics((; grids=grids(pf, s_vent)), coords, false)
        @test isapprox(m_all.pf_vent_dev, 0.0; atol=1.0e-15)
        @test isequal(m_all.pf_vent_dev, 0.0)
        inside(i, j) =
            hypot(coords.xp[j] - SCV.X_CENTRE, coords.yp[i] - SCV.Y_CENTRE) <= SCV.R_PLANET
        # For each side, a surface cell whose only outside neighbour lies on that side:
        # the cells next to the poles have y-facing exterior faces only.
        sides = ((-1, 0), (1, 0), (0, -1), (0, 1))
        only_outside(side) = filter(cells) do (i, j)
            outside = [!inside(i + di, j + dj) for (di, dj) in sides]
            return count(outside) == 1 && !inside(i + side[1], j + side[2])
        end
        # The grid is symmetric under reflection and x <-> y: each side has as many cells.
        n_side = [length(only_outside(side)) for side in sides]
        @test all(==(n_side[1]), n_side)
        @test 1 <= n_side[1] < length(cells) / 4
        for side in sides, (i0, j0) in only_outside(side)
            # A surface cell 0.3 Pe_R above P_vent sets the deviation.
            pf_high = copy(pf)
            pf_high[i0, j0] = SCV.P_AMB + 0.3 * SCV.PE_R
            m_high = SCV._metrics((; grids=grids(pf_high, s_vent)), coords, false)
            @test isapprox(m_high.pf_vent_dev, 0.3; rtol=1.0e-12)
            # A surface cell that does not vent makes the deviation infinite, whatever its Pf.
            s_closed = copy(s_vent)
            s_closed[i0, j0] = 0.0
            m_closed = SCV._metrics((; grids=grids(pf, s_closed)), coords, false)
            @test isequal(m_closed.pf_vent_dev, Inf)
        end
    end

    @testset "metrics reject an empty mask" begin
        near = compaction_coords(41)
        far = compaction_coords(41; shift=10 * SCV.R_PLANET)
        @test_throws ArgumentError SCV._vs_l2_inner(
            exact_solid_velocity_grids(far, 3, 1.0), far, 3
        )
        @test isfinite(SCV._vs_l2_inner(exact_solid_velocity_grids(near, 3, 1.0), near, 3))
        zero_grid = zeros(far.Ny1, far.Nx1)
        grids = (;
            pr=zero_grid,
            pf=zero_grid,
            ETAPHI=fill(SCV.ZETA, far.Ny1, far.Nx1),
            PHI=fill(SCV.PHI, far.Ny1, far.Nx1),
            S_vent_grid=zero_grid,
            vx=zero_grid,
            vy=zero_grid,
        )
        @test_throws ArgumentError SCV._metrics((; grids), far, false)
    end

    @testset "metrics_from_runs with every grid of NX_LIST" begin
        fake(e) = (;
            centre_err=e,
            l2_inner=2e,
            vs_l2_inner=3e,
            pe_centre=1.0 + e,
            vent_closure=1.0e-12,
            pf_vent_dev=1.0e-9,
        )
        fake_runs = Dict{Tuple,Any}()
        for (n, e) in ((33, 0.04), (65, 0.02), (129, 0.01)), sf in (false, true)
            fake_runs[(n, sf)] = fake(sf ? e : 1.5e)
        end
        all_metrics = SCV.metrics_from_runs(fake_runs)
        @test Set(keys(all_metrics)) == Set(SCV.EXPECTED_METRIC_IDS)
        # The error halves from Nx = 65 to 129, so every order is 1.
        for tag in ("on", "off"), name in ("centre", "profile", "vs")
            @test isapprox(all_metrics["$(name)_order_$(tag)_65_129"], 1.0; rtol=1.0e-12)
        end
        @test isapprox(all_metrics["centre_ratio_on_off_nx129"], 1.01 / 1.015; rtol=1.0e-14)
        delete!(fake_runs, (129, true))
        partial = SCV.metrics_from_runs(fake_runs)
        @test !haskey(partial, "centre_order_on_65_129")
        @test !haskey(partial, "centre_ratio_on_off_nx129")
        @test haskey(partial, "centre_order_off_65_129")
    end

    @testset "shared record helpers: canonical text and Julia version match" begin
        @test isequal(SCV._canonical_field_repr([1.0, 2.5]), "[1.0, 2.5]")
        @test isequal(
            SCV._canonical_field_repr(Dict(:b => 2.0, :a => 1)), "{:a => 1, :b => 2.0}"
        )
        @test SCV._canonical_field_repr(Dict()) == "{}"
        @test SCV._canonical_field_repr((true, :x, "s", nothing)) ==
            "[true, :x, \"s\", nothing]"
        @test SCV._canonical_field_repr(v"1.2.3") == "1.2.3"
        @test SCV._same_julia_major_minor("1.12.1", v"1.12.7")
        @test SCV._same_julia_major_minor("v1.12", v"1.12.0")
        @test !SCV._same_julia_major_minor("1.11.0", v"1.12.0")
        @test !SCV._same_julia_major_minor("julia", v"1.12.0")
        # The conduction and Darcy modules use the same shared definitions.
        conduction = validation_module(validation_case("spherical_conduction"))
        darcy = validation_module(validation_case("spherical_darcy_flow"))
        @test isequal(conduction._canonical_field_repr((1, 2.0)), "[1, 2.0]")
        @test darcy._same_julia_major_minor("1.12.0", v"1.12.3")
        @test validation_case("not_a_record") === nothing
    end

    @testset "validation modules load into the given host module" begin
        host = Module(:ValidationHost)
        case = validation_case("spherical_darcy_flow")
        @test !isdefined(host, case.module_name)
        loaded = validation_module(case; host=host)
        @test isdefined(host, case.module_name)
        @test loaded === getfield(host, case.module_name)
        # A second call returns the same module and loads nothing new.
        @test validation_module(case; host=host) === loaded
        @test validation_module(case) !== loaded
        # The export tool reads the validation module of the module that includes it.
        @test SphericalCompactionExport.compaction_validation_module() === SCV
        @test SphericalCompactionExport.compaction_validation_module((@__MODULE__,)) === SCV
        @test_throws ErrorException SphericalCompactionExport.compaction_validation_module((
            Module(:Empty),
        ))
        Base.include(
            host,
            joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_compaction.jl"),
        )
        Base.include(
            host,
            joinpath(@__DIR__, "..", "tools", "export_validation_spherical_compaction.jl"),
        )
        host_scv = getfield(host, :SphericalCompactionValidation)
        host_export = getfield(host, :SphericalCompactionExport)
        @test Base.invokelatest(host_export.compaction_validation_module) === host_scv
        @test host_scv !== SCV
    end

    @testset "Validation record schema and drift guards" begin
        @test isfile(record_path)
        raw_record = TOML.parsefile(record_path)
        valid_ok, valid_errors = SCV.validate_record_dict(
            raw_record; require_doc_pages=true
        )
        @test valid_ok
        @test isempty(valid_errors)

        rec_ids = [String(m["id"]) for m in raw_record["metrics"]]
        @test rec_ids == collect(SCV.EXPECTED_METRIC_IDS)
        @test Set(keys(SCV.PHYSICS_REGRESSION_BANDS)) == Set(SCV.EXPECTED_METRIC_IDS)

        missing_key = deepcopy(raw_record)
        delete!(missing_key["config"], "zeta_Pa_s")
        ok_missing, errs_missing = SCV.validate_record_dict(
            missing_key; require_doc_pages=false, check_doc_sync=false
        )
        @test !ok_missing
        @test any(occursin("zeta_Pa_s", e) for e in errs_missing)

        for (mid, value) in bench.metrics
            band = SCV.PHYSICS_REGRESSION_BANDS[mid]
            @test isapprox(value, recorded[mid]; atol=band.atol, rtol=band.rtol)
        end
        ok_full, errs_full = SCV.check_physics_regression_bands(bench.metrics, record_path)
        @test ok_full
        @test isempty(errs_full)

        target = "vs_l2_inner_on_nx65"
        moved = copy(bench.metrics)
        moved[target] += 2 * SCV.PHYSICS_REGRESSION_BANDS[target].atol
        ok_moved, errs_moved = SCV.check_physics_regression_bands(moved, record_path)
        @test !ok_moved
        @test length(errs_moved) == 1
        @test occursin(target, only(errs_moved))
        missing_metric = copy(bench.metrics)
        delete!(missing_metric, target)
        ok_gone, errs_gone = SCV.check_physics_regression_bands(missing_metric, record_path)
        @test !ok_gone
        @test any(occursin("Measured metrics missing '$target'", e) for e in errs_gone)
        @test_throws DomainError SCV.check_physics_regression_bands(bench.metrics, "")
        @test !first(SCV.check_physics_regression_bands(bench.metrics, record_path * ".x"))
    end

    @testset "Record matches the configuration and the reference values" begin
        raw_record = TOML.parsefile(record_path)
        cfg = SCV.make_config(SCV.NX_REFERENCE; spherical_flow=true)
        sha = SCV.canonical_config_sha256(cfg)
        @test sha == raw_record["metadata"]["config_sha256"]
        @test length(sha) == 64
        @test SCV.canonical_config_sha256(
            SCV.make_config(SCV.NX_REFERENCE; spherical_flow=false)
        ) != sha
        @test SCV.canonical_config_sha256(
            SCV.make_config(SCV.NX_REFERENCE; spherical_flow=true, eta_ratio=1.0e-3)
        ) != sha

        text = SCV.canonical_config_text(cfg)
        used = (
            :L_BOX,
            :R_PLANET,
            :X_CENTRE,
            :Y_CENTRE,
            :NX_REFERENCE,
            :NX_LIST,
            :INNER_FRACTION,
            :MU_F,
            :K_PERM,
            :PHI,
            :DELTA,
            :ZETA,
            :ETA_RATIO,
            :ETA_AIR,
            :PE_R,
            :P_AMB,
            :T_SURF,
            :K_VENT,
            :DT_YR,
            :RX_AIR,
            :RHO_SOLID,
            :RHOCP,
            :RHO_FLUID,
            :STRENGTH_FILL,
            :SERIES_SWITCH,
        )
        @test Set(used) == Set(SCV.HASHED_CASE_CONSTANTS)
        for name in used
            @test occursin("case.$(name)=", text)
        end
        @test bytes2hex(SHA.sha256(text)) == sha

        summary = SCV.expected_config_summary(cfg)
        rec_cfg = raw_record["config"]
        @test Set(keys(summary)) == Set(keys(rec_cfg))
        for (key, expected) in summary
            if expected isa AbstractFloat
                @test isapprox(Float64(rec_cfg[key]), expected; rtol=1.0e-12, atol=0.0)
            else
                @test rec_cfg[key] == expected
            end
        end

        refs = SCV.expected_metric_references()
        @test Set(keys(refs)) == Set(SCV.EXPECTED_METRIC_IDS)
        @test isapprox(
            refs["centre_ratio_on_off_nx129"], SCV.CENTRE_RATIO_EXACT; rtol=1.0e-15
        )
        for m in raw_record["metrics"]
            @test isapprox(
                Float64(m["reference"]), refs[String(m["id"])]; rtol=1.0e-12, atol=1.0e-12
            )
        end
    end

    @testset "Reproducibility check, CLI table and export round trip" begin
        raw_record = TOML.parsefile(record_path)
        ok_repro, errs_repro = SCV.check_validation_reproducibility(bench, record_path)
        @test ok_repro
        @test isempty(errs_repro)
        @test run_compaction_validation_check(record_path; bench=bench)
        @test run_validation_cli(["--check", record_path]; compaction_bench=bench) == 0
        @test run_validation_cli(["--bad-flag"]) == 2
        @test !run_compaction_validation_check(record_path * ".missing"; bench=bench)
        @test_throws DomainError run_compaction_validation_check(""; bench=bench)

        mktempdir() do tmpdir
            # A record with an id outside the table is a schema error, not an exception.
            unknown = deepcopy(raw_record)
            unknown["metadata"]["id"] = "not_a_record"
            unknown_path = write_compaction_toml(joinpath(tmpdir, "unknown.toml"), unknown)
            @test run_validation_cli([unknown_path]; compaction_bench=bench) == 1
            @test run_validation_cli([joinpath(tmpdir, "absent.toml")]) == 1
            broken_path = joinpath(tmpdir, "broken.toml")
            write(broken_path, "[metadata\nid = ")
            @test run_validation_cli([broken_path]) == 1
            schema_bad = deepcopy(raw_record)
            delete!(schema_bad, "config")
            schema_path = write_compaction_toml(joinpath(tmpdir, "schema.toml"), schema_bad)
            @test !run_compaction_validation_check(schema_path; bench=bench)

            # The no-argument run checks every record of the table in order; the conduction
            # and Darcy benches here are empty, so only the compaction record passes.
            empty_bench = (; runs=Dict{Tuple,Any}(), metrics=Dict{String,Float64}())
            out_path = joinpath(tmpdir, "cli_out.txt")
            err_path = joinpath(tmpdir, "cli_err.txt")
            code = open(out_path, "w") do out
                open(err_path, "w") do err
                    redirect_stdio(; stdout=out, stderr=err) do
                        run_validation_cli(
                            ["--check"];
                            bench=empty_bench,
                            darcy_bench=empty_bench,
                            compaction_bench=bench,
                        )
                    end
                end
            end
            @test code == 1
            out_text = read(out_path, String)
            err_text = read(err_path, String)
            @test occursin("spherical_compaction.toml", out_text)
            @test occursin("[reproducibility]", err_text)
            @test [c.id for c in VALIDATION_CASES] == ["spherical_conduction", "spherical_darcy_flow", "spherical_compaction"]

            # Export round trip.
            exported = joinpath(tmpdir, "exported_spherical_compaction.toml")
            commit = String(raw_record["metadata"]["src_commit"])
            @test SphericalCompactionExport.run_export_cli(
                [exported]; src_commit=commit, bench=bench, check_clean_tree=false
            ) == 0
            raw_exported = TOML.parsefile(exported)
            ok_valid, errs_valid = SCV.validate_record_dict(
                raw_exported; require_doc_pages=false, check_doc_sync=false
            )
            @test ok_valid
            @test isempty(errs_valid)
            @test run_compaction_validation_check(exported; bench=bench)
            @test [String(m["id"]) for m in raw_exported["metrics"]] == [String(m["id"]) for m in raw_record["metrics"]]

            # A moved value is refused without the flag and written with it.
            target = "centre_err_on_nx65"
            band = SCV.PHYSICS_REGRESSION_BANDS[target]
            moved_metrics = copy(bench.metrics)
            moved_metrics[target] += 2 * band.atol
            moved_bench = (; runs=bench.runs, metrics=moved_metrics)
            @test SphericalCompactionExport.run_export_cli(
                [exported]; src_commit=commit, bench=moved_bench, check_clean_tree=false
            ) == 1
            @test SphericalCompactionExport.run_export_cli(
                ["--accept-physics-change", exported];
                src_commit=commit,
                bench=moved_bench,
                check_clean_tree=false,
            ) == 0
            @test SphericalCompactionExport.run_export_cli(
                ["--bad-flag"]; src_commit=commit, bench=bench, check_clean_tree=false
            ) == 2
            @test SphericalCompactionExport.run_export_cli(
                [exported];
                src_commit=commit,
                bench=bench,
                check_clean_tree=true,
                dirty_status=" M src/config.jl",
            ) == 1
            @test_throws DomainError SphericalCompactionExport.run_export_cli(
                [exported]; src_commit="abc", bench=bench, check_clean_tree=false
            )
            # A value that fails its own tolerance is never written.
            failing = copy(bench.metrics)
            failing["vent_closure_on_nx65"] = 1.0
            @test SphericalCompactionExport.run_export_cli(
                ["--accept-physics-change", exported];
                src_commit=commit,
                bench=(; runs=bench.runs, metrics=failing),
                check_clean_tree=false,
            ) == 1

            status = SphericalCompactionExport.metric_status
            cases = (
                (0.01, 0.0, "le_atol"), (0.01, 1.0, "le_atol"), (1.0, 1.001, "isapprox")
            )
            @test isequal(
                [status(v, ref, 0.05, 0.0, c) for (v, ref, c) in cases],
                ["PASS", "FAIL", "PASS"],
            )
            @test isequal(status(1.0, 1.0, 0.0, 0.0, "other"), "FAIL")
            @test occursin(r"^[0-9a-f]{40}$", SphericalCompactionExport.last_src_commit())
            @test first(SphericalCompactionExport.check_clean_source_tree()) isa Bool

            # A changed configuration value and hash are rejected.
            changed_cfg = deepcopy(raw_record)
            changed_cfg["config"]["k_m2"] *= 2
            changed_cfg["metadata"]["config_sha256"] = repeat("0", 64)
            cfg_path = write_compaction_toml(joinpath(tmpdir, "changed.toml"), changed_cfg)
            ok_cfg, errs_cfg = SCV.check_validation_reproducibility(bench, cfg_path)
            @test !ok_cfg
            @test any(occursin("config.k_m2", e) for e in errs_cfg)
            @test any(occursin("config_sha256", e) for e in errs_cfg)
        end
    end

    @testset "checks stop on an interrupt and report other errors with a stack trace" begin
        commit = String(raw_record_all["metadata"]["src_commit"])
        empty_bench = (; runs=Dict{Tuple,Any}(), metrics=Dict{String,Float64}())
        @test_throws InterruptException run_validation_cli(
            ["--check"];
            bench=InterruptingBench(),
            darcy_bench=empty_bench,
            compaction_bench=bench,
        )
        @test_throws InterruptException run_compaction_validation_check(
            record_path; bench=InterruptingBench()
        )
        crash = capture_stdio(
            () -> run_compaction_validation_check(record_path; bench=CrashingBench())
        )
        @test crash.value === false
        @test occursin("solver failure inside the benchmark", crash.err)
        @test occursin("Stacktrace", crash.err)
        mktempdir() do tmpdir
            target = joinpath(tmpdir, "record.toml")
            @test_throws InterruptException SphericalCompactionExport.run_export_cli(
                [target];
                src_commit=commit,
                bench=InterruptingBench(),
                check_clean_tree=false,
            )
            failed = capture_stdio(
                () -> SphericalCompactionExport.run_export_cli(
                    [target];
                    src_commit=commit,
                    bench=CrashingBench(),
                    check_clean_tree=false,
                ),
            )
            @test failed.value == 1
            @test occursin("Stacktrace", failed.err)
            @test !isfile(target)
            # A clean-tree check whose `git status` fails is reported and gives exit code 1.
            not_repo = mkpath(joinpath(tmpdir, "not_a_repo"))
            withenv("GIT_DIR" => joinpath(not_repo, "absent.git")) do
                for cli in (
                    SphericalCompactionExport.run_export_cli,
                    SphericalCompactionExport.export_main,
                )
                    no_git = capture_stdio(
                        () -> cli(
                            [target];
                            src_commit=commit,
                            bench=bench,
                            check_clean_tree=true,
                            repo_root=not_repo,
                        ),
                    )
                    @test no_git.value == 1
                    @test occursin("failed process", no_git.err)
                    @test occursin("status --porcelain", no_git.err)
                    @test !isfile(target)
                end
            end
        end
    end

    @testset "command lines: one record path, unknown flags and an empty path" begin
        commit = String(raw_record_all["metadata"]["src_commit"])
        mktempdir() do tmpdir
            missing_path = joinpath(tmpdir, "missing.toml")
            two = capture_stdio(
                () ->
                    run_validation_cli([missing_path, record_path]; compaction_bench=bench),
            )
            @test two.value == 2
            @test occursin("Usage", two.err)
            @test run_validation_cli(["-h"]) == 2
            @test run_validation_cli(["--check", "-x", record_path]) == 2
            cd(tmpdir) do
                code = SphericalCompactionExport.run_export_cli(
                    ["-h"]; src_commit=commit, bench=bench, check_clean_tree=false
                )
                @test code == 2
                @test !isfile(joinpath(tmpdir, "-h"))
            end
            p1, p2 = joinpath(tmpdir, "first.toml"), joinpath(tmpdir, "second.toml")
            @test SphericalCompactionExport.run_export_cli(
                [p1, p2]; src_commit=commit, bench=bench, check_clean_tree=false
            ) == 2
            @test !isfile(p1) && !isfile(p2)
            # The script entry point reports an empty path and exits with 1.
            empty_run = capture_stdio(
                () -> SphericalCompactionExport.export_main(
                    [""]; src_commit=commit, bench=bench, check_clean_tree=false
                ),
            )
            @test empty_run.value == 1
            @test occursin("Record path must be non-empty", empty_run.err)
            @test SphericalCompactionExport.export_main(
                [p2]; src_commit=commit, bench=bench, check_clean_tree=false
            ) == 0
            @test isfile(p2)
        end
    end

    @testset "malformed records are reported, not thrown" begin
        mktempdir() do tmpdir
            broken = joinpath(tmpdir, "broken.toml")
            write(broken, "[metadata\nid = ")
            parsed = capture_stdio(
                () -> run_compaction_validation_check(broken; bench=bench)
            )
            @test parsed.value === false
            @test occursin("[schema] Cannot parse", parsed.err)
            bad_meta = joinpath(tmpdir, "bad_meta.toml")
            write(bad_meta, "metadata = 1\n[config]\n[[metrics]]\n")
            ok_meta, errs_meta = SCV.check_validation_reproducibility(bench, bad_meta)
            @test !ok_meta
            @test any(occursin("metadata", e) for e in errs_meta)
            bad_value = deepcopy(raw_record_all)
            bad_value["metrics"][1]["value"] = "not a number"
            bad_value_path = write_compaction_toml(
                joinpath(tmpdir, "bad_value.toml"), bad_value
            )
            ok_value, errs_value = SCV.check_validation_reproducibility(
                bench, bad_value_path
            )
            @test !ok_value
            @test any(occursin("numeric", e) for e in errs_value)
            # A record that is not TOML fails the direct checks with a message.
            for check in (
                () -> SCV.check_validation_reproducibility(bench, broken),
                () -> SCV.check_physics_regression_bands(bench.metrics, broken),
            )
                ok_parse, errs_parse = check()
                @test !ok_parse
                @test occursin("Cannot parse", only(errs_parse))
            end
            # Valid TOML with a number or a Bool where the schema has a string (or an integer
            # for `nx`) gives schema errors in every check, not a MethodError.
            wrong_types = (
                ("metadata", "src_commit", 123),
                ("metadata", "config_sha256", 123),
                ("metadata", "id", 123),
                ("metadata", "title", true),
                ("metadata", "page", 7),
                ("metric", "id", 123),
                ("metric", "case", 1),
                ("metric", "description", false),
                ("metric", "comparison", 1),
                ("metric", "units", 1),
                ("metric", "status", true),
                ("metric", "nx", "33"),
            )
            for (section, key, value) in wrong_types
                rec = deepcopy(raw_record_all)
                table = section == "metadata" ? rec["metadata"] : rec["metrics"][1]
                table[key] = value
                path = write_compaction_toml(joinpath(tmpdir, "$(section)_$key.toml"), rec)
                ok_dict, errs_dict = SCV.validate_record_dict(TOML.parsefile(path))
                @test !ok_dict
                @test any(occursin("'$key' must be", e) for e in errs_dict)
                ok_rep, errs_rep = SCV.check_validation_reproducibility(bench, path)
                @test !ok_rep
                @test errs_rep == errs_dict
                checked = capture_stdio(
                    () -> run_compaction_validation_check(path; bench=bench)
                )
                @test checked.value === false
                @test occursin("[schema]", checked.err)
                cli = capture_stdio(
                    () -> run_validation_cli([path]; compaction_bench=bench)
                )
                @test cli.value == 1
                @test occursin("[schema]", cli.err)
            end
            # A missing metric key is reported once, without the type checks of that metric.
            no_units = deepcopy(raw_record_all)
            delete!(no_units["metrics"][1], "units")
            ok_units, errs_units = SCV.validate_record_dict(no_units)
            @test !ok_units
            @test errs_units == ["Metric #1 missing key 'units'"]
        end
    end

    @testset "figure data export: CSV columns and the exported files" begin
        script_host = Module(:CompactionFigureData)
        # A bare `Module` has no `include`; the script calls it at top level.
        Core.eval(script_host, :(include(path) = Base.include($script_host, path)))
        Base.include(
            script_host,
            joinpath(@__DIR__, "..", "scripts", "export_validation_compaction_data.jl"),
        )
        write_columns = getfield(script_host, :write_columns)
        export_all = getfield(script_host, :export_all)
        mktempdir() do tmpdir
            path = joinpath(tmpdir, "columns.csv")
            write(path, "kept\n")
            for (header, columns) in (
                (["a"], Vector{Float64}[]),
                (["a", "b", "c"], [[1.0, 2.0], [3.0, 4.0]]),
                (["a", "b"], [[1.0, 2.0], [3.0, 4.0, 5.0]]),
            )
                @test_throws ArgumentError Base.invokelatest(
                    write_columns, path, header, columns
                )
            end
            # A rejected call leaves the file as it was.
            @test read(path, String) == "kept\n"
            Base.invokelatest(write_columns, path, ["a", "b"], [[1.0, 2.0], [3, 4]])
            @test readlines(path) == ["a,b", "1.0,3.0", "2.0,4.0"]

            outdir = joinpath(tmpdir, "data")
            @test Base.invokelatest(export_all, outdir) == outdir
            files = sort(readdir(outdir))
            @test files == sort([
                "convergence.csv",
                "profile_off_65.csv",
                "profile_on_65.csv",
                "reference_curves.csv",
            ])
            conv = readlines(joinpath(outdir, "convergence.csv"))
            @test conv[1] == "nx,l2_off,l2_on,centre_off,centre_on"
            @test length(conv) == 1 + length(SCV.NX_LIST)
            # The Nx = 65 row carries the recorded errors.
            row65 = parse.(Float64, split(conv[3], ","))
            @test isapprox(row65[2], recorded["profile_l2_inner_off_nx65"]; rtol=1.0e-15)
            @test isapprox(row65[5], recorded["centre_err_on_nx65"]; rtol=1.0e-15)
            # The exact curves end at Pe_R on the surface for both settings.
            last_curve = parse.(
                Float64,
                split(readlines(joinpath(outdir, "reference_curves.csv"))[end], ","),
            )
            @test isapprox(last_curve[1], SCV.R_PLANET; rtol=1.0e-15)
            @test isapprox(last_curve[2], SCV.PE_R; rtol=1.0e-12)
            @test isapprox(last_curve[3], SCV.PE_R; rtol=1.0e-12)
            # The node profile of the sphere follows the exact curve within the profile bound.
            prof = readlines(joinpath(outdir, "profile_on_65.csv"))
            @test prof[1] == "r_m,pe_Pa,pe_exact_Pa"
            vals = [parse.(Float64, split(l, ",")) for l in prof[2:end]]
            @test length(vals) == runs[(65, true)].n_rock
            @test all(abs(v[2] / v[3] - 1) < 0.1 for v in vals)
        end
    end
end

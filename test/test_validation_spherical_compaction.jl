using SHA
using Test
using TOML
using Erebus

include(joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_compaction.jl"))
include(joinpath(@__DIR__, "..", "tools", "check_validation.jl"))
include(joinpath(@__DIR__, "..", "tools", "export_validation_spherical_compaction.jl"))

const SCV = SphericalCompactionValidation

"""
    write_compaction_toml(path, dict)

Write `dict` as TOML to `path` and return `path`.
"""
function write_compaction_toml(path::String, dict::AbstractDict)
    open(path, "w") do io
        return TOML.print(io, dict)
    end
    return path
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
        @test runs[(33, sf)].converged
        @test runs[(65, sf)].converged
        @test all(isfinite, runs[(65, sf)].pe) && all(isfinite, runs[(65, sf)].vx)
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
            @test m.pf_vent_dev < 1.0e-5
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

    @testset "a long time step reaches the vent drain cap and breaks the closure" begin
        long = SCV.run_case(33; spherical_flow=true, dt_yr=1.0e-2)
        @test long.converged
        @test !isapprox(long.vent_closure, 0.0; atol=1.0e-3)
        @test runs[(33, true)].vent_closure < 1.0e-8
    end

    @testset "6-variable path equals the 4-variable path (spherical_flow = $sf)" for sf in (
        false, true
    )
        four = runs[(33, sf)]
        six = SCV.run_case(33; spherical_flow=sf, darcy_elimination=false)
        @test six.converged
        @test maximum(abs, six.pf .- four.pf) / maximum(abs, four.pf) < 1.0e-10
        @test maximum(abs, six.pe .- four.pe) / maximum(abs, four.pe) < 1.0e-8
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
end

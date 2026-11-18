using SHA
using Test
using TOML
using Erebus

include(joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_darcy_flow.jl"))
include(joinpath(@__DIR__, "..", "tools", "check_validation.jl"))
include(joinpath(@__DIR__, "..", "tools", "export_validation_spherical_darcy_flow.jl"))

const SDF = SphericalDarcyFlowValidation

"""
    write_toml(path, dict)

Write `dict` as TOML to `path` and return `path`.
"""
function write_toml(path::String, dict::AbstractDict)
    open(path, "w") do io
        return TOML.print(io, dict)
    end
    return path
end

"""
    synthetic_coords(N; shift=0.0)

Return a `NamedTuple` with the staggered coordinates that the flux and error measures read,
laid out as `GridCoordinates` does for an `N` x `N` grid of the box `SDF.L_BOX`, with every
coordinate moved by `shift` [m] in x and y.

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
function synthetic_coords(N::Int; shift::Float64=0.0)
    d = SDF.L_BOX / (N - 1)
    node(lo, hi) = collect(range(lo, hi; length=N + 1)) .+ shift
    return (;
        Nx=N,
        Ny=N,
        Nx1=N + 1,
        Ny1=N + 1,
        dx=d,
        dy=d,
        xvx=node(0.0, SDF.L_BOX + d),
        yvx=node(-d / 2, SDF.L_BOX + d / 2),
        xvy=node(-d / 2, SDF.L_BOX + d / 2),
        yvy=node(0.0, SDF.L_BOX + d),
        xp=node(-d / 2, SDF.L_BOX + d / 2),
        yp=node(-d / 2, SDF.L_BOX + d / 2),
    )
end

"""
    exact_flux_grids(coords, d, scale)

Return the grids `qxD` and `qyD` as a `NamedTuple` holding `scale` times the exact flux
`S r_vec / d` [m/s] at every face of `coords`, with `d` the flow dimension.

Parameters
----------
coords : NamedTuple
    Coordinates from `synthetic_coords`.
d : Float64
    Flow dimension, 3 for the sphere and 2 for the disk.
scale : Float64
    Factor applied to the exact flux.

Returns
-------
NamedTuple
    `qxD` and `qyD`, each of size `(coords.Ny1, coords.Nx1)`.
"""
function exact_flux_grids(coords, d::Float64, scale::Float64)
    qxD = [
        scale * SDF.S_SOURCE * (x - SDF.X_CENTRE) / d for _ in coords.yvx, x in coords.xvx
    ]
    qyD = [
        scale * SDF.S_SOURCE * (y - SDF.Y_CENTRE) / d for y in coords.yvy, _ in coords.xvy
    ]
    return (; qxD, qyD)
end

@testset "Spherical Darcy Flow Validation" begin
    record_path = joinpath(
        @__DIR__, "..", "benchmarks", "records", "spherical_darcy_flow.toml"
    )
    raw_record_all = TOML.parsefile(record_path)
    recorded = Dict(
        String(m["id"]) => Float64(m["value"]) for m in raw_record_all["metrics"]
    )
    # Runs at Nx = 33 and 65 only. The Nx = 129 and order metrics take their recorded values,
    # which `validate_record_dict` checks against their reference and tolerance.
    bench_runs = SDF.evaluate_spherical_darcy_flow_benchmark(; nx_list=(33, 65))
    runs = bench_runs.runs
    measured = bench_runs.metrics
    n_fine = maximum(SDF.NX_LIST)
    recorded_only = filter(
        kv ->
            !haskey(measured, kv.first) &&
            (endswith(kv.first, "_nx$(n_fine)") || startswith(kv.first, "profile_order")),
        recorded,
    )
    bench = (; runs, metrics=merge(measured, recorded_only))

    @testset "metrics: 19 from the runs, 10 recorded values for Nx = 129 and the orders" begin
        @test length(measured) == 19
        @test length(recorded_only) == 10
        @test isempty(intersect(keys(measured), keys(recorded_only)))
        @test Set(keys(bench.metrics)) == Set(SDF.EXPECTED_METRIC_IDS)
    end

    @testset "surface pressure is the ambient pressure" begin
        cfg = SDF.make_config(SDF.NX_REFERENCE; spherical_flow=true)
        p_vent = Erebus.compute_venting_pressure(
            SDF.T_SURF, SDF.P_AMB; species=:H2O, L_sub=cfg.venting.L_sublimation
        )
        @test isapprox(p_vent, SDF.P_AMB; rtol=1.0e-14, atol=0.0)
        @test Erebus.compute_water_vapor_pressure(SDF.T_SURF) < SDF.P_AMB
    end

    @testset "solver converges (spherical_flow = $sf)" for sf in (false, true)
        @test runs[(33, sf)].converged
        @test runs[(65, sf)].converged
        @test all(isfinite, runs[(65, sf)].pf) &&
            all(isfinite, runs[(65, sf)].qxD) &&
            all(isfinite, runs[(65, sf)].qyD)
    end

    @testset "grid below the minimum is rejected" begin
        @test_throws ArgumentError SDF.make_config(SDF.NX_MIN - 1; spherical_flow=true)
        @test_throws ArgumentError SDF.run_case(SDF.NX_MIN - 1; spherical_flow=false)
        @test SDF.make_config(SDF.NX_MIN; spherical_flow=true).grid.Nx == SDF.NX_MIN
        @test_throws ArgumentError SDF.evaluate_spherical_darcy_flow_benchmark(;
            nx_list=(17,)
        )
        @test_throws ArgumentError SDF.evaluate_spherical_darcy_flow_benchmark(;
            nx_list=(33.0,)
        )
    end

    @testset "excess pressure profile against the exact solution (spherical_flow = $sf)" for sf in
                                                                                             (
        false, true
    )
        e33 = runs[(33, sf)].l2_inner
        e65 = runs[(65, sf)].l2_inner
        @test isapprox(e65, 0.0; atol=0.08)
        @test isapprox(log2(e33 / e65), 1.0; atol=0.3)
    end

    @testset "Darcy flux vector against S r_vec / d (spherical_flow = $sf)" for sf in (
        false, true
    )
        f33 = runs[(33, sf)].flux_l2_inner
        f65 = runs[(65, sf)].flux_l2_inner
        @test isapprox(f33, 0.0; atol=0.21)
        @test isapprox(f65, 0.0; atol=0.041)
        # The error falls faster than first order between the 2 grids.
        @test f65 < f33 / 2
    end

    @testset "flux vector error is 0 for the exact flux and scales with the departure" begin
        coords = synthetic_coords(41)
        for d in (2.0, 3.0)
            err(scale) = SDF._flux_l2_inner(exact_flux_grids(coords, d, scale), coords, d)
            @test isapprox(err(1.0), 0.0; atol=1.0e-12)
            @test isapprox(err(1.1), 0.1; atol=1.0e-12)
            @test isapprox(err(-1.0), 2.0; atol=1.0e-12)
        end
    end

    @testset "metrics reject an empty mask" begin
        near = synthetic_coords(41)
        far = synthetic_coords(41; shift=10 * SDF.R_PLANET)
        flux = exact_flux_grids(far, 3.0, 1.0)
        @test_throws ArgumentError SDF._flux_l2_inner(flux, far, 3.0)
        @test isfinite(SDF._flux_l2_inner(exact_flux_grids(near, 3.0, 1.0), near, 3.0))
        zero_grid = zeros(far.Ny1, far.Nx1)
        grids = (;
            pf=zero_grid,
            qxD=flux.qxD,
            qyD=flux.qyD,
            S_vent_grid=zero_grid,
            DQPF=zero_grid,
            vx=zero_grid,
            vy=zero_grid,
        )
        @test_throws ArgumentError SDF._metrics((; grids), far, false)
    end

    @testset "vented flux equals the discrete source, matrix stays rigid (spherical_flow = $sf)" for sf in
                                                                                                     (
        false, true
    )
        m = runs[(65, sf)]
        @test isapprox(m.vent_flux, m.source_flux; rtol=1.0e-8)
        @test isapprox(m.max_solid_speed, 0.0; atol=1.0e-12)
    end

    @testset "centre excess pressure on/off is 2/3" begin
        on, off = runs[(65, true)], runs[(65, false)]
        @test isapprox(on.dp_centre / off.dp_centre, 2 / 3; atol=0.01)
        # The staircase surface lowers both centre values by the same fraction.
        @test isapprox(
            on.dp_centre / on.dp_centre_exact,
            off.dp_centre / off.dp_centre_exact;
            rtol=0.01,
        )
    end

    @testset "permeability jump: exact profile and solver against the piecewise solution" begin
        a, R = SDF.A_JUMP, SDF.R_PLANET
        for d in (2.0, 3.0)
            pref = SDF.MU_F * SDF.S_SOURCE / (2d)
            # The excess pressure vanishes at the surface and is continuous across the jump.
            @test isapprox(SDF.exact_excess_pressure(R, d; jump=true), 0.0; atol=1.0e-9)
            p_out = SDF.exact_excess_pressure(a, d; jump=true)
            @test isapprox(p_out, pref * (R^2 - a^2) / SDF.K_PERM; rtol=1.0e-12)
            @test isapprox(
                SDF.exact_excess_pressure(a * (1 - 1.0e-9), d; jump=true),
                p_out;
                rtol=1.0e-6,
            )
            # The centre value adds the core drop mu S a^2 / (2 d K_CORE).
            @test isapprox(
                SDF.exact_excess_pressure(0.0, d; jump=true),
                p_out + pref * a^2 / SDF.K_CORE;
                rtol=1.0e-12,
            )
            # Without the jump the profile is the uniform one and the exponent is not 1.
            @test isapprox(
                SDF.exact_excess_pressure(0.0, d), pref * R^2 / SDF.K_PERM; rtol=1.0e-12
            )
            @test SDF.exact_excess_pressure(0.0, d; jump=true) >
                10 * SDF.exact_excess_pressure(0.0, d)
        end
        @test SDF.K_CORE ≈ SDF.K_PERM / 100 rtol = 1.0e-12
        @test SDF.A_JUMP ≈ SDF.R_PLANET / 2 rtol = 1.0e-12

        for sf in (false, true)
            r33, r65 = runs[(33, sf, :jump)], runs[(65, sf, :jump)]
            @test r33.converged && r65.converged
            @test all(isfinite, r65.pf) && all(isfinite, r65.qxD) && all(isfinite, r65.qyD)
            # Measured 0.034 at Nx = 65 and 0.138 at Nx = 33; the error depends on how the
            # grid cuts the circle r = A_JUMP, so the test asserts limits and no order.
            @test isapprox(r65.l2_inner, 0.0; atol=0.06)
            @test r65.l2_inner < r33.l2_inner
            # Measured 0.024; the faces near the jump are excluded.
            @test isfinite(r65.flux_l2_inner) && r65.flux_l2_inner < 0.04
            # Measured 0.983: the centre value is below the exact one at Nx = 65.
            @test isapprox(r65.dp_centre / r65.dp_centre_exact, 1.0; atol=0.03)
            @test r65.dp_centre < r65.dp_centre_exact
            @test isapprox(r65.vent_flux, r65.source_flux; rtol=1.0e-8)
        end
        # The grid cut of r = A_JUMP sets the error: the finer grid Nx = 65 has a larger
        # error than Nx = 63 (measured 0.034 and 0.025), so the record keeps no order.
        r63 = SDF.run_case(63; spherical_flow=true, jump=true)
        @test r63.l2_inner < runs[(65, true, :jump)].l2_inner
        # The centre value scales with the flow dimension d as 1 / d.
        on, off = runs[(65, true, :jump)], runs[(65, false, :jump)]
        @test isapprox(on.dp_centre / off.dp_centre, 2 / 3; atol=0.005)
        @test_throws DomainError SDF.exact_excess_pressure(-1.0, 3.0; jump=true)
        @test_throws DomainError SDF.exact_excess_pressure(SDF.R_PLANET * 1.01, 3.0)
        @test_throws DomainError SDF.exact_excess_pressure(0.0, 1.0)
    end

    @testset "jump flux error leaves out the faces within 2 cells of the jump" begin
        coords = synthetic_coords(41)
        near = [
            abs(hypot(x - SDF.X_CENTRE, y - SDF.Y_CENTRE) - SDF.A_JUMP) <=
            SDF.JUMP_EXCLUDE_CELLS * coords.dx for y in coords.yvx, x in coords.xvx
        ]
        near_y = [
            abs(hypot(x - SDF.X_CENTRE, y - SDF.Y_CENTRE) - SDF.A_JUMP) <=
            SDF.JUMP_EXCLUDE_CELLS * coords.dx for y in coords.yvy, x in coords.xvy
        ]
        @test any(near) && any(near_y)
        # Each component on its own: both loops must leave out the faces near the jump.
        for (field, mask) in ((:qxD, near), (:qyD, near_y))
            bad = exact_flux_grids(coords, 3.0, 1.0)
            getfield(bad, field)[mask] .*= 3.0
            @test !isapprox(SDF._flux_l2_inner(bad, coords, 3.0), 0.0; atol=0.1)
            @test isapprox(
                SDF._flux_l2_inner(bad, coords, 3.0; jump=true), 0.0; atol=1.0e-12
            )
            # A face outside the band still counts with the jump.
            far = exact_flux_grids(coords, 3.0, 1.0)
            getfield(far, field)[.!mask] .*= 3.0
            @test !isapprox(SDF._flux_l2_inner(far, coords, 3.0; jump=true), 0.0; atol=0.1)
        end
    end

    @testset "6-variable path equals the 4-variable path (spherical_flow = $sf, jump = $jump)" for sf in
                                                                                                   (
            false, true
        ),
        jump in (false, true)

        four = jump ? runs[(33, sf, :jump)] : runs[(33, sf)]
        six = SDF.run_case(33; spherical_flow=sf, darcy_elimination=false, jump=jump)
        @test six.converged
        @test maximum(abs, six.pf .- four.pf) / maximum(abs, four.pf) < 1.0e-8
        @test isapprox(six.flux_l2_inner, four.flux_l2_inner; rtol=1.0e-6)
        @test maximum(abs, six.qxD .- four.qxD) / maximum(abs, four.qxD) < 1.0e-8
        @test maximum(abs, six.qyD .- four.qyD) / maximum(abs, four.qyD) < 1.0e-8
    end

    @testset "Validation record schema and drift guards" begin
        @test isfile(record_path)
        raw_record = TOML.parsefile(record_path)
        valid_ok, valid_errors = SDF.validate_record_dict(
            raw_record; require_doc_pages=true
        )
        @test valid_ok
        @test isempty(valid_errors)

        rec_ids = [String(m["id"]) for m in raw_record["metrics"]]
        @test rec_ids == collect(SDF.EXPECTED_METRIC_IDS)
        @test Set(keys(SDF.PHYSICS_REGRESSION_BANDS)) == Set(SDF.EXPECTED_METRIC_IDS)
        @test Set(keys(bench.metrics)) == Set(SDF.EXPECTED_METRIC_IDS)

        # The record needs the config keys of its own id, and an unknown id is an error.
        missing_key = deepcopy(raw_record)
        delete!(missing_key["config"], "k_vent_m2")
        ok_missing, errs_missing = SDF.validate_record_dict(
            missing_key; require_doc_pages=false, check_doc_sync=false
        )
        @test !ok_missing
        @test any(occursin("k_vent_m2", e) for e in errs_missing)
        unknown_id = deepcopy(raw_record)
        unknown_id["metadata"]["id"] = "not_a_record"
        ok_unknown, errs_unknown = SDF.validate_record_dict(
            unknown_id; require_doc_pages=false, check_doc_sync=false
        )
        @test !ok_unknown
        @test any(occursin("Unknown record id", e) for e in errs_unknown)

        # Every metric stays within the regression band of the record: measured at Nx = 33
        # and 65, the recorded value for the others.
        for (mid, value) in bench.metrics
            band = SDF.PHYSICS_REGRESSION_BANDS[mid]
            @test isapprox(value, recorded[mid]; atol=band.atol, rtol=band.rtol)
        end
        ok_full, errs_full = SDF.check_physics_regression_bands(bench.metrics, record_path)
        @test ok_full
        @test isempty(errs_full)

        # A value moved by more than its band is the only one reported, with its metric id.
        target = "profile_l2_inner_on_nx65"
        moved = copy(bench.metrics)
        moved[target] += 2 * SDF.PHYSICS_REGRESSION_BANDS[target].atol
        @test Set(keys(moved)) == Set(SDF.EXPECTED_METRIC_IDS)
        ok_moved, errs_moved = SDF.check_physics_regression_bands(moved, record_path)
        @test !ok_moved
        @test length(errs_moved) == 1
        @test occursin(target, only(errs_moved))
        @test occursin("outside physics regression band", only(errs_moved))
    end

    @testset "Record matches the configuration and the reference values" begin
        raw_record = TOML.parsefile(record_path)
        cfg = SDF.make_config(SDF.NX_REFERENCE; spherical_flow=true)

        sha = SDF.canonical_config_sha256(cfg)
        @test sha == raw_record["metadata"]["config_sha256"]
        @test length(sha) == 64
        @test SDF.canonical_config_sha256(cfg) == sha
        @test SDF.canonical_config_sha256(
            SDF.make_config(SDF.NX_REFERENCE; spherical_flow=false)
        ) != sha
        @test SDF.canonical_config_sha256(SDF.make_config(33; spherical_flow=true)) != sha

        # The hashed text lists the geometry, the grid list and the constants of `set_fields!`.
        text = SDF.canonical_config_text(cfg)
        used = (
            :L_BOX,
            :R_PLANET,
            :X_CENTRE,
            :Y_CENTRE,
            :NX_REFERENCE,
            :NX_LIST,
            :INNER_FRACTION,
            :ETA_SOLID,
            :RX_AIR,
            :PF_INIT,
            :PHI,
            :K_PERM,
            :K_CORE,
            :A_JUMP,
            :JUMP_EXCLUDE_CELLS,
            :MU_F,
            :S_SOURCE,
            :P_AMB,
            :T_SURF,
            :RHO_SOLID,
            :RHOCP,
            :RHO_FLUID,
            :STRENGTH_FILL,
        )
        @test Set(used) == Set(SDF.HASHED_CASE_CONSTANTS)
        for name in used
            @test occursin("case.$(name)=", text)
        end
        @test occursin("case.NX_LIST=[33, 65, 129]\n", text)
        @test bytes2hex(SHA.sha256(text)) == sha

        summary = SDF.expected_config_summary(cfg)
        rec_cfg = raw_record["config"]
        @test Set(keys(summary)) == Set(keys(rec_cfg))
        for (key, expected) in summary
            if expected isa AbstractFloat
                @test isapprox(Float64(rec_cfg[key]), expected; rtol=1.0e-12, atol=0.0)
            else
                @test rec_cfg[key] == expected
            end
        end

        refs = SDF.expected_metric_references()
        @test Set(keys(refs)) == Set(String(m["id"]) for m in raw_record["metrics"])
        for m in raw_record["metrics"]
            @test isapprox(
                Float64(m["reference"]), refs[String(m["id"])]; rtol=1.0e-12, atol=1.0e-12
            )
        end
    end

    @testset "Reproducibility check and export round trip" begin
        raw_record = TOML.parsefile(record_path)
        ok_repro, errs_repro = SDF.check_validation_reproducibility(bench, record_path)
        @test ok_repro
        @test isempty(errs_repro)
        @test run_darcy_validation_check(record_path; bench=bench)
        @test run_validation_cli(["--check", record_path]; darcy_bench=bench) == 0

        status = SphericalDarcyFlowExport.metric_status
        # (value, reference) pairs: only a finite value within atol of a zero reference passes.
        cases = ((0.01, 0.0), (0.01, 1.0), (0.01, NaN), (NaN, 0.0), (0.06, 0.0))
        statuses = [status(v, ref, 0.05, 0.0, "le_atol") for (v, ref) in cases]
        @test statuses == ["PASS", "FAIL", "FAIL", "FAIL", "FAIL"]

        mktempdir() do tmpdir
            exported = joinpath(tmpdir, "exported_spherical_darcy_flow.toml")
            @test SphericalDarcyFlowExport.run_export_cli(
                [exported];
                src_commit=String(raw_record["metadata"]["src_commit"]),
                bench=bench,
                check_clean_tree=false,
            ) == 0
            @test isfile(exported)
            raw_exported = TOML.parsefile(exported)
            ok_valid, errs_valid = SDF.validate_record_dict(
                raw_exported; require_doc_pages=false, check_doc_sync=false
            )
            @test ok_valid
            @test isempty(errs_valid)
            @test run_darcy_validation_check(exported; bench=bench)
            @test [String(m["id"]) for m in raw_exported["metrics"]] == [String(m["id"]) for m in raw_record["metrics"]]

            # A value moved beyond its band is rejected with its metric id.
            target = "flux_l2_inner_on_nx65"
            band = SDF.PHYSICS_REGRESSION_BANDS[target]
            moved_value = deepcopy(raw_record)
            entry = only(filter(m -> m["id"] == target, moved_value["metrics"]))
            entry["value"] += 2 * band.atol
            moved_path = write_toml(joinpath(tmpdir, "moved_value.toml"), moved_value)
            ok_moved, errs_moved = SDF.check_validation_reproducibility(bench, moved_path)
            @test !ok_moved
            @test any(occursin(target, e) for e in errs_moved)
            @test !run_darcy_validation_check(moved_path; bench=bench)

            # A changed configuration value is rejected with its key.
            changed_cfg = deepcopy(raw_record)
            changed_cfg["config"]["k_m2"] *= 2
            cfg_path = write_toml(joinpath(tmpdir, "changed_config.toml"), changed_cfg)
            ok_cfg, errs_cfg = SDF.check_validation_reproducibility(bench, cfg_path)
            @test !ok_cfg
            @test any(occursin("config.k_m2", e) for e in errs_cfg)
            @test !run_darcy_validation_check(cfg_path; bench=bench)

            # A changed configuration hash is rejected.
            changed_sha = deepcopy(raw_record)
            changed_sha["metadata"]["config_sha256"] = repeat("0", 64)
            sha_path = write_toml(joinpath(tmpdir, "changed_sha.toml"), changed_sha)
            ok_sha, errs_sha = SDF.check_validation_reproducibility(bench, sha_path)
            @test !ok_sha
            @test any(occursin("config_sha256", e) for e in errs_sha)
        end
    end
end

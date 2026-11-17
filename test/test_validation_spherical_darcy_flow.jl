using Test
using TOML
using Erebus

include(joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_darcy_flow.jl"))

const SDF = SphericalDarcyFlowValidation

@testset "Spherical Darcy Flow Validation" begin
    runs = Dict(
        (n, sf) => SDF.run_case(n; spherical_flow=sf) for n in (33, 65), sf in (false, true)
    )

    @testset "excess pressure profile against the exact solution (spherical_flow = $sf)" for sf in
                                                                                             (
        false, true
    )
        e33 = runs[(33, sf)].l2_inner
        e65 = runs[(65, sf)].l2_inner
        @test isapprox(e65, 0.0; atol=0.08)
        @test isapprox(log2(e33 / e65), 1.0; atol=0.3)
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

    @testset "Validation record schema and drift guards" begin
        record_path = joinpath(
            @__DIR__, "..", "benchmarks", "records", "spherical_darcy_flow.toml"
        )
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

        # Rerun values at Nx = 33 and 65 stay within the regression band of the record.
        measured = Dict{String,Float64}()
        for (sf, tag) in ((false, "off"), (true, "on"))
            for n in (33, 65)
                measured["profile_l2_inner_$(tag)_nx$(n)"] = runs[(n, sf)].l2_inner
            end
            m65 = runs[(65, sf)]
            measured["vent_closure_$(tag)_nx65"] = abs(m65.vent_flux / m65.source_flux - 1)
        end
        measured["centre_ratio_on_off_nx65"] =
            runs[(65, true)].dp_centre / runs[(65, false)].dp_centre
        measured["solid_speed_max_nx65"] = max(
            runs[(65, false)].max_solid_speed, runs[(65, true)].max_solid_speed
        )
        recorded = Dict(
            String(m["id"]) => Float64(m["value"]) for m in raw_record["metrics"]
        )
        for (mid, value) in measured
            band = SDF.PHYSICS_REGRESSION_BANDS[mid]
            @test isapprox(value, recorded[mid]; atol=band.atol, rtol=band.rtol)
        end

        # A value moved by more than its band is reported with its metric id.
        moved = copy(measured)
        moved["profile_l2_inner_on_nx65"] +=
            2 * SDF.PHYSICS_REGRESSION_BANDS["profile_l2_inner_on_nx65"].atol
        ok_moved, errs_moved = SDF.check_physics_regression_bands(moved, record_path)
        @test !ok_moved
        @test any(occursin("profile_l2_inner_on_nx65", e) for e in errs_moved)
    end
end

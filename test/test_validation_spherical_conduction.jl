using Test
using TOML
using Erebus

if !isdefined(@__MODULE__, :validate_record_dict)
    include(joinpath(@__DIR__, "..", "tools", "validation_records.jl"))
end
if !isdefined(@__MODULE__, :run_doc_numbers_cli)
    include(joinpath(@__DIR__, "..", "tools", "check_doc_numbers.jl"))
end
if !isdefined(@__MODULE__, :evaluate_spherical_conduction_benchmark)
    include(joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_conduction.jl"))
end
if !isdefined(@__MODULE__, :run_validation_check)
    include(joinpath(@__DIR__, "..", "tools", "check_validation.jl"))
end
if !isdefined(@__MODULE__, :write_spherical_conduction_record)
    include(joinpath(@__DIR__, "..", "tools", "export_validation_spherical_conduction.jl"))
end

@testset "Spherical Conduction Validation" begin
    record_path = joinpath(
        @__DIR__, "..", "benchmarks", "records", "spherical_conduction.toml"
    )

    @testset "Validation record schema and drift guards" begin
        @test isfile(record_path)
        @test check_validation_records()

        raw_record = TOML.parsefile(record_path)
        valid_ok, valid_errors = validate_record_dict(raw_record; require_doc_pages=false)
        @test valid_ok
        @test isempty(valid_errors)

        rendered = render_validation_markdown_block(raw_record)
        @test occursin("<!-- VALIDATION_RECORD_START: spherical_conduction -->", rendered)
        @test occursin("<!-- VALIDATION_RECORD_END: spherical_conduction -->", rendered)

        mktempdir() do tmpdir
            rel_page = "validation/thermal_conduction.md"
            full_page = joinpath(tmpdir, rel_page)
            mkpath(dirname(full_page))
            doc_record = deepcopy(raw_record)
            doc_record["metadata"]["page"] = rel_page
            write(full_page, "# Spherical Conduction\n\n" * rendered * "\n")

            ok_sync, errs_sync = validate_record_dict(
                doc_record; docs_dir=tmpdir, require_doc_pages=true
            )
            @test ok_sync
            @test isempty(errs_sync)

            # Mutating one metric within schema tolerance must still fail markdown synchronization.
            drifted_record = deepcopy(doc_record)
            drifted_record["metrics"][1]["value"] += 0.001
            ok_drift, errs_drift = validate_record_dict(
                drifted_record; docs_dir=tmpdir, require_doc_pages=true
            )
            @test !ok_drift
            @test any(occursin("out of sync", e) for e in errs_drift)

            # Removing the markdown markers from the page must fail synchronization.
            write(full_page, "# Spherical Conduction\n\nNo fenced block here.\n")
            ok_missing, errs_missing = validate_record_dict(
                doc_record; docs_dir=tmpdir, require_doc_pages=true
            )
            @test !ok_missing
            @test any(occursin("out of sync", e) for e in errs_missing)
        end

        # Invalid commit and config hashes and failing metric statuses must be rejected.
        bad_commit = deepcopy(raw_record)
        bad_commit["metadata"]["git_commit"] = "not-a-40-char-sha"
        ok_commit, _ = validate_record_dict(bad_commit; require_doc_pages=false)
        @test !ok_commit

        bad_sha = deepcopy(raw_record)
        bad_sha["metadata"]["config_sha256"] = "short"
        ok_sha, _ = validate_record_dict(bad_sha; require_doc_pages=false)
        @test !ok_sha

        bad_status = deepcopy(raw_record)
        bad_status["metrics"][1]["status"] = "FAIL"
        ok_status, _ = validate_record_dict(bad_status; require_doc_pages=false)
        @test !ok_status

        bad_value = deepcopy(raw_record)
        bad_value["metrics"][1]["value"] = bad_value["metrics"][1]["reference"] + 10.0
        ok_val, _ = validate_record_dict(bad_value; require_doc_pages=false)
        @test !ok_val

        bad_comp = deepcopy(raw_record)
        bad_comp["metrics"][1]["comparison"] = "unknown_mode"
        ok_comp, _ = validate_record_dict(bad_comp; require_doc_pages=false)
        @test !ok_comp

        # Schema structure and required-key rejection branches.
        @test !validate_record_dict(Dict{String,Any}())[1]
        @test !validate_record_dict(
            Dict{String,Any}("metadata" => "bad", "config" => Dict(), "metrics" => Any[])
        )[1]
        missing_meta = deepcopy(raw_record)
        delete!(missing_meta["metadata"], "id")
        @test !validate_record_dict(missing_meta)[1]
        bad_cfg_tbl = deepcopy(raw_record)
        bad_cfg_tbl["config"] = "bad"
        @test !validate_record_dict(bad_cfg_tbl)[1]
        missing_cfg_key = deepcopy(raw_record)
        delete!(missing_cfg_key["config"], "Nx_fine")
        @test !validate_record_dict(missing_cfg_key)[1]
        empty_metrics = deepcopy(raw_record)
        empty_metrics["metrics"] = Any[]
        @test !validate_record_dict(empty_metrics)[1]
        nontbl_metric = deepcopy(raw_record)
        nontbl_metric["metrics"] = Any["not_a_table"]
        @test !validate_record_dict(nontbl_metric)[1]
        missing_mkey = deepcopy(raw_record)
        delete!(missing_mkey["metrics"][1], "atol")
        @test !validate_record_dict(missing_mkey)[1]
        dup_mid = deepcopy(raw_record)
        dup_mid["metrics"][2]["id"] = dup_mid["metrics"][1]["id"]
        @test !validate_record_dict(dup_mid)[1]
        @test !validate_record_dict(
            raw_record;
            docs_dir=joinpath(@__DIR__, "nonexistent_docs"),
            require_doc_pages=true,
        )[1]

        # Canonical config SHA-256 must match the committed record and detect changes.
        base_cfg = make_spherical_conduction_config(; Nx=65)
        @test canonical_config_sha256(base_cfg) == raw_record["metadata"]["config_sha256"]
        altered_cfg = make_spherical_conduction_config(; Nx=65, k_rock=3.5)
        @test canonical_config_sha256(altered_cfg) !=
            raw_record["metadata"]["config_sha256"]
    end

    @testset "Analytical and 1D finite-volume reference solutions" begin
        @test_throws DomainError analytical_T_BE(-1.0)
        @test_throws DomainError analytical_T_BE(0.0; R_b=0.0)
        @test_throws DomainError analytical_T_BE(0.0; dt_s=0.0)
        @test_throws DomainError analytical_T_BE(0.0; n_modes=0)
        @test_throws DomainError analytical_T_continuous(-1.0)
        @test_throws DomainError analytical_T_continuous(0.0; t_tot_s=0.0)
        @test_throws DomainError solve_1d_spherical_fv(1)
        @test_throws DomainError solve_1d_spherical_fv(10; dt_s=0.0)

        # Isothermal surface boundary and radial monotonicity from centre to surface.
        T_center = analytical_T_BE(0.0; n_modes=500)
        T_mid = analytical_T_BE(0.5 * R_PLANET; n_modes=500)
        T_bound = analytical_T_BE(R_PLANET; n_modes=500)
        T_cont_center = analytical_T_continuous(0.0; n_modes=500)
        T_cont_mid = analytical_T_continuous(0.5 * R_PLANET; n_modes=500)
        @test isapprox(T_bound, T_SURF; atol=1.0e-10)
        @test isapprox(analytical_T_continuous(R_PLANET; n_modes=500), T_SURF; atol=1.0e-10)
        @test T_bound < T_mid < T_center
        @test T_bound < T_cont_mid < T_cont_center

        # Fine-timestep backward-Euler series converges to the continuous solution as dt -> 0.
        T_be_fine = analytical_T_BE(
            0.0; n_steps=2000, dt_s=DT_S / 100.0, n_modes=500, include_step_average=true
        )
        @test isapprox(T_be_fine - T_SURF, T_cont_center - T_SURF; rtol=2.0e-3)

        # Removable resonance singularities in continuous and discrete eigenmode series.
        R_res_cont = pi * sqrt(KAPPA / LAMBDA_AL)
        R_res_be = pi * sqrt((1.0 / exp(-LAMBDA_AL * DT_S) - 1.0) / (KAPPA * DT_S))
        @test isfinite(analytical_T_continuous(0.0; R_b=R_res_cont, n_modes=10))
        @test isfinite(analytical_T_BE(0.0; R_b=R_res_be, n_modes=10))

        # Mode truncation convergence between 500 and 1000 eigenmodes.
        dT_500 = T_center - T_SURF
        dT_1000 = analytical_T_BE(0.0; n_modes=1000) - T_SURF
        mode_rel_diff = abs(dT_500 - dT_1000) / dT_1000
        @test mode_rel_diff <= 1.0e-7

        # Cross-verification against independent 1D spherical finite-volume reference.
        r_fv, T_fv = solve_1d_spherical_fv(500)
        dT_be_0 = T_center - T_SURF
        dT_fv_r1 = T_fv[1] - T_SURF
        rel_diff_be_1dfv = abs(dT_be_0 - dT_fv_r1) / dT_be_0
        @test rel_diff_be_1dfv <= 1.0e-5

        # Sensitivity test: omitting step averaging breaks reference agreement by >= 10,000x.
        T_be_no_avg = analytical_T_BE(0.0; n_modes=500, include_step_average=false)
        dT_be_no_avg = T_be_no_avg - T_SURF
        broken_rel_diff = abs(dT_be_no_avg - dT_fv_r1) / dT_fv_r1
        min_broken_ref_diff = 0.10
        @test broken_rel_diff >= min_broken_ref_diff
    end

    bench = evaluate_spherical_conduction_benchmark()

    @testset "Steady spherical Poisson conduction" begin
        s33 = bench.s33
        s65 = bench.s65

        @test isapprox(s33.ratio_3d_2d, 2.0 / 3.0; atol=0.01, rtol=0.015)
        @test isapprox(s65.ratio_3d_2d, 2.0 / 3.0; atol=0.01, rtol=0.015)
        @test s33.l2_rel <= 0.15
        @test s65.l2_rel <= 0.08
        @test s65.l2_rel < s33.l2_rel
        @test s65.energy_res <= 1.0e-11

        # Canary sensitivity test: no-op configuration reproduces baseline and passes guards.
        s33_noop = run_steady_spherical_conduction(
            33;
            spherical_metric=true,
            metric_scale=1.0,
            clamp_exterior=true,
            perturb_stencil=false,
        )
        @test isapprox(s33_noop.ratio_3d_2d, s33.ratio_3d_2d; rtol=1.0e-12)
        @test isapprox(s33_noop.l2_rel, s33.l2_rel; rtol=1.0e-12)

        # Sensitivity tests: each broken physics assumption fails its guard by >= 2x.
        min_ratio_break_65 = 2.0 * 0.01
        min_l2_break_65 = 2.0 * 0.08
        min_energy_break_65 = 2.0 * 1.0e-11

        s65_no_metric = run_steady_spherical_conduction(65; spherical_metric=false)
        @test abs(s65_no_metric.ratio_3d_2d - 2.0 / 3.0) >= min_ratio_break_65
        @test s65_no_metric.l2_rel >= min_l2_break_65

        s65_unclamped = run_steady_spherical_conduction(65; clamp_exterior=false)
        @test s65_unclamped.l2_rel >= min_l2_break_65

        s65_perturbed = run_steady_spherical_conduction(65; perturb_stencil=true)
        @test s65_perturbed.energy_res >= min_energy_break_65

        s65_doubled_source = run_steady_spherical_conduction(65; metric_scale=2.0)
        @test abs(s65_doubled_source.ratio_3d_2d - 2.0 / 3.0) >= min_ratio_break_65
        @test s65_doubled_source.l2_rel >= min_l2_break_65
    end

    @testset "Transient spherical radiogenic conduction and reproducibility" begin
        t33 = bench.t33
        t65 = bench.t65

        # Coarse and fine centre temperature rise and L2 bounds against step-averaged reference.
        @test abs(t33.center_err_K) <= 60.0
        @test t33.l2_rel <= 0.25
        @test isapprox(t65.dT_center, t65.dT_be_rmin; atol=25.0, rtol=0.05)

        # Centre error magnitude |sim - T_BE| in K falls by a factor >= 2 from Nx=33 to Nx=65.
        center_err_reduction = abs(t33.center_err_K) / abs(t65.center_err_K)
        min_reduction_factor = 2.0
        @test center_err_reduction >= min_reduction_factor

        # Interior profile L2 error, mesh convergence ratio, and sticky-air isothermal clamp.
        conv_ratio_65_33 = t65.l2_rel / t33.l2_rel
        @test t65.l2_rel <= 0.10
        @test conv_ratio_65_33 <= 0.60
        @test isapprox(t65.T_max_air_shell, T_SURF; atol=0.01)

        # Canary sensitivity test: no-op configuration reproduces baseline.
        t33_noop = run_transient_spherical_conduction(33; dsubgridt=0.0)
        @test isapprox(t33_noop.dT_center, t33.dT_center; rtol=1.0e-12)
        @test isapprox(t33_noop.l2_rel, t33.l2_rel; rtol=1.0e-12)

        # Sensitivity test 1: disabling spherical metric breaks centre, L2, and convergence by >= 2x.
        t33_no_metric = run_transient_spherical_conduction(33; spherical_metric=false)
        t65_no_metric = run_transient_spherical_conduction(65; spherical_metric=false)
        min_center_rel_break = 2.0 * 0.05
        min_l2_break = 2.0 * 0.10
        min_conv_break = 2.0 * 0.60
        @test abs(t65_no_metric.dT_center - t65_no_metric.dT_be_rmin) /
              t65_no_metric.dT_be_rmin >= min_center_rel_break
        @test t65_no_metric.l2_rel >= min_l2_break
        @test (t65_no_metric.l2_rel / t33_no_metric.l2_rel) >= min_conv_break

        # Sensitivity test 2: disabling sticky-air reset breaks air-shell clamp by >= 2x.
        t65_no_reset = run_transient_spherical_conduction(65; n_steps=5, disk_enabled=false)
        min_air_clamp_break = 2.0 * 0.01
        @test abs(t65_no_reset.T_max_air_shell - T_SURF) >= min_air_clamp_break

        # Additional physics sensitivity tests at Nx=33 (conductivity, air heat capacity, regularization).
        t33_zero_k = run_transient_spherical_conduction(33; k_rock=1.0e-6, k_air=1.0e-6)
        @test abs(t33_zero_k.dT_center - t33_zero_k.dT_be_rmin) / t33_zero_k.dT_be_rmin >=
            min_center_rel_break
        @test t33_zero_k.l2_rel >= min_l2_break

        t33_low_air_cp = run_transient_spherical_conduction(33; rhocp_air=3.3e6)
        @test abs(t33_low_air_cp.T_max_air_shell - T_SURF) >= min_air_clamp_break
        @test t33_low_air_cp.l2_rel >= min_l2_break

        t33_overreg = run_transient_spherical_conduction(33; reg_cells=15.0)
        @test abs(t33_overreg.dT_center - t33_overreg.dT_be_rmin) /
              t33_overreg.dT_be_rmin >= min_center_rel_break

        # Full benchmark evaluation and record reproducibility check.
        ok_repro, repro_errs = check_validation_reproducibility(bench, record_path)
        @test ok_repro
        @test isempty(repro_errs)
        @test run_validation_check(record_path; bench=bench)

        # Mutating a recorded metric or config hash in the record must fail reproducibility.
        mktempdir() do tmpdir
            exported_path = joinpath(tmpdir, "exported_spherical_conduction.toml")
            write_spherical_conduction_record(exported_path; bench=bench)
            @test isfile(exported_path)
            @test run_validation_check(exported_path; bench=bench)

            @test !run_validation_check(joinpath(tmpdir, "missing.toml"); bench=bench)
            @test !check_validation_records(joinpath(tmpdir, "missing_dir"))
            empty_records_dir = joinpath(tmpdir, "empty_records")
            mkpath(empty_records_dir)
            @test !check_validation_records(empty_records_dir)

            # Schema-invalid record file in directory and CLI entry points.
            invalid_schema_dir = joinpath(tmpdir, "invalid_records")
            mkpath(invalid_schema_dir)
            invalid_schema_path = joinpath(invalid_schema_dir, "bad.toml")
            write(invalid_schema_path, "[metadata]\nid = \"bad\"\n")
            @test !check_validation_records(invalid_schema_dir)
            @test !run_validation_check(invalid_schema_path; bench=bench)
            @test run_validation_cli(["--check", record_path]; bench=bench) == 0
            @test run_validation_cli(["--unknown-flag"]; bench=bench) == 2
            @test run_validation_cli([invalid_schema_path]; bench=bench) == 1
            @test run_doc_numbers_cli(String[]) == 0
            @test run_doc_numbers_cli([joinpath(tmpdir, "missing_map.toml")]) == 1
            @test run_doc_numbers_cli(String[]; records_dir=empty_records_dir) == 1

            # Drift within schema tolerance must be caught by recorded-value reproducibility.
            tampered_path = joinpath(tmpdir, "spherical_conduction.toml")
            raw_rec = TOML.parsefile(record_path)
            raw_rec["metrics"][1]["value"] += 0.001
            open(tampered_path, "w") do io
                TOML.print(io, raw_rec)
            end
            ok_tampered, errs_tampered = check_validation_reproducibility(
                bench, tampered_path; check_recorded_values=true
            )
            @test !ok_tampered
            @test any(occursin("differs from recorded value", e) for e in errs_tampered)
            @test !run_validation_check(tampered_path; bench=bench)

            bad_cfg_path = joinpath(tmpdir, "spherical_conduction_cfg.toml")
            raw_cfg_rec = TOML.parsefile(record_path)
            raw_cfg_rec["metadata"]["config_sha256"] = repeat("0", 64)
            open(bad_cfg_path, "w") do io
                TOML.print(io, raw_cfg_rec)
            end
            ok_bad_cfg, errs_bad_cfg = check_validation_reproducibility(bench, bad_cfg_path)
            @test !ok_bad_cfg
            @test !isempty(errs_bad_cfg)
            @test !run_validation_check(bad_cfg_path; bench=bench)
        end
    end
end

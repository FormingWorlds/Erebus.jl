using Test
using TOML
using Erebus

include(joinpath(@__DIR__, "..", "benchmarks", "validation_spherical_conduction.jl"))
include(joinpath(@__DIR__, "..", "tools", "check_validation.jl"))
include(joinpath(@__DIR__, "..", "tools", "export_validation_spherical_conduction.jl"))

@testset "Spherical Conduction Validation" begin
    record_path = joinpath(
        @__DIR__, "..", "benchmarks", "records", "spherical_conduction.toml"
    )

    @testset "Validation record schema and drift guards" begin
        @test isfile(record_path)
        @test SphericalConductionValidation.check_validation_records()

        raw_record = TOML.parsefile(record_path)
        valid_ok, valid_errors = SphericalConductionValidation.validate_record_dict(
            raw_record; require_doc_pages=false
        )
        @test valid_ok
        @test isempty(valid_errors)

        # Exact metric ID set in canonical order.
        rec_ids = [String(m["id"]) for m in raw_record["metrics"]]
        @test rec_ids == collect(SphericalConductionValidation.EXPECTED_METRIC_IDS)

        rendered = SphericalConductionValidation.render_validation_markdown_block(
            raw_record
        )
        @test occursin("<!-- VALIDATION_RECORD_START: spherical_conduction -->", rendered)
        @test occursin("<!-- VALIDATION_RECORD_END: spherical_conduction -->", rendered)
        @test occursin("Source commit", rendered)

        mktempdir() do tmpdir
            rel_page = "validation/thermal_conduction.md"
            full_page = joinpath(tmpdir, rel_page)
            mkpath(dirname(full_page))
            doc_record = deepcopy(raw_record)
            doc_record["metadata"]["page"] = rel_page
            write(full_page, "# Spherical Conduction\n\n" * rendered * "\n")

            ok_sync, errs_sync = SphericalConductionValidation.validate_record_dict(
                doc_record; docs_dir=tmpdir, require_doc_pages=true
            )
            @test ok_sync
            @test isempty(errs_sync)

            # Mutating one metric within schema tolerance must still fail markdown synchronization.
            drifted_record = deepcopy(doc_record)
            drifted_record["metrics"][1]["value"] += 0.001
            ok_drift, errs_drift = SphericalConductionValidation.validate_record_dict(
                drifted_record; docs_dir=tmpdir, require_doc_pages=true
            )
            @test !ok_drift
            @test any(occursin("out of sync", e) for e in errs_drift)

            # Removing the markdown markers from the page must fail synchronization.
            write(full_page, "# Spherical Conduction\n\nNo fenced block here.\n")
            ok_missing, errs_missing = SphericalConductionValidation.validate_record_dict(
                doc_record; docs_dir=tmpdir, require_doc_pages=true
            )
            @test !ok_missing
            @test any(occursin("out of sync", e) for e in errs_missing)
        end

        # Invalid src_commit and config hashes, failing statuses, and non-zero le_atol references must fail.
        bad_commit = deepcopy(raw_record)
        bad_commit["metadata"]["src_commit"] = "not-a-40-char-sha"
        ok_commit, _ = SphericalConductionValidation.validate_record_dict(
            bad_commit; require_doc_pages=false
        )
        @test !ok_commit

        bad_sha = deepcopy(raw_record)
        bad_sha["metadata"]["config_sha256"] = "short"
        ok_sha, _ = SphericalConductionValidation.validate_record_dict(
            bad_sha; require_doc_pages=false
        )
        @test !ok_sha

        bad_status = deepcopy(raw_record)
        bad_status["metrics"][1]["status"] = "FAIL"
        ok_status, _ = SphericalConductionValidation.validate_record_dict(
            bad_status; require_doc_pages=false
        )
        @test !ok_status

        bad_value = deepcopy(raw_record)
        bad_value["metrics"][1]["value"] = bad_value["metrics"][1]["reference"] + 10.0
        ok_val, _ = SphericalConductionValidation.validate_record_dict(
            bad_value; require_doc_pages=false
        )
        @test !ok_val

        bad_comp = deepcopy(raw_record)
        bad_comp["metrics"][1]["comparison"] = "unknown_mode"
        ok_comp, _ = SphericalConductionValidation.validate_record_dict(
            bad_comp; require_doc_pages=false
        )
        @test !ok_comp

        bad_le_ref = deepcopy(raw_record)
        bad_le_ref["metrics"][3]["reference"] = 0.5
        ok_le_ref, errs_le_ref = SphericalConductionValidation.validate_record_dict(
            bad_le_ref; require_doc_pages=false
        )
        @test !ok_le_ref
        @test any(occursin("must have reference = 0.0", e) for e in errs_le_ref)

        # Schema structure and required-key rejection branches.
        @test !SphericalConductionValidation.validate_record_dict(Dict{String,Any}())[1]
        @test !SphericalConductionValidation.validate_record_dict(
            Dict{String,Any}("metadata" => "bad", "config" => Dict(), "metrics" => Any[])
        )[1]
        missing_meta = deepcopy(raw_record)
        delete!(missing_meta["metadata"], "id")
        @test !SphericalConductionValidation.validate_record_dict(missing_meta)[1]
        bad_cfg_tbl = deepcopy(raw_record)
        bad_cfg_tbl["config"] = "bad"
        @test !SphericalConductionValidation.validate_record_dict(bad_cfg_tbl)[1]
        missing_cfg_key = deepcopy(raw_record)
        delete!(missing_cfg_key["config"], "Nx_fine")
        @test !SphericalConductionValidation.validate_record_dict(missing_cfg_key)[1]
        empty_metrics = deepcopy(raw_record)
        empty_metrics["metrics"] = Any[]
        @test !SphericalConductionValidation.validate_record_dict(empty_metrics)[1]
        nontbl_metric = deepcopy(raw_record)
        nontbl_metric["metrics"] = Any["not_a_table"]
        @test !SphericalConductionValidation.validate_record_dict(nontbl_metric)[1]
        missing_mkey = deepcopy(raw_record)
        delete!(missing_mkey["metrics"][1], "atol")
        @test !SphericalConductionValidation.validate_record_dict(missing_mkey)[1]
        dup_mid = deepcopy(raw_record)
        dup_mid["metrics"][2]["id"] = dup_mid["metrics"][1]["id"]
        @test !SphericalConductionValidation.validate_record_dict(dup_mid)[1]
        @test !SphericalConductionValidation.validate_record_dict(
            raw_record;
            docs_dir=joinpath(@__DIR__, "nonexistent_docs"),
            require_doc_pages=true,
        )[1]

        # Canonical config SHA-256: hashes only benchmark physics sections and steady inputs.
        base_cfg = SphericalConductionValidation.make_spherical_conduction_config(; Nx=65)
        base_sha = SphericalConductionValidation.canonical_config_sha256(base_cfg)
        @test base_sha == raw_record["metadata"]["config_sha256"]

        altered_cfg = SphericalConductionValidation.make_spherical_conduction_config(;
            Nx=65, k_rock=3.5
        )
        @test SphericalConductionValidation.canonical_config_sha256(altered_cfg) != base_sha
        @test SphericalConductionValidation.canonical_config_sha256(
            base_cfg; H_steady=2.0e-4
        ) != base_sha
        @test SphericalConductionValidation.canonical_config_sha256(
            base_cfg; L_box=150_000.0
        ) != base_sha
        @test SphericalConductionValidation.canonical_config_sha256(
            base_cfg; reg_cells=1.0
        ) != base_sha
        @test SphericalConductionValidation.canonical_config_sha256(base_cfg; k_cond=4.0) !=
            base_sha

        # Unrelated sections or extra fields in unrelated sections must not change the hash.
        unrelated_cfg = SimulationConfig(;
            grid=base_cfg.grid,
            geometry=base_cfg.geometry,
            time=base_cfg.time,
            solver=SolverConfig(; seed=99999, dsubgridt=0.5),
            thermodynamics=base_cfg.thermodynamics,
            materials=base_cfg.materials,
            poroelasticity=base_cfg.poroelasticity,
            output=OutputConfig(; savematstep=1, visstep=1, save_final=true),
            reaction=ReactionConfig(; active=true),
            melting=MeltingConfig(; active=true),
            coreformation=CoreFormationConfig(; percolation_active=true),
            metal_partition=MetalPartitionConfig(; active=true),
            venting=VentingConfig(; active=true),
            volatiles=VolatilesConfig(; active=true),
            retention=RetentionConfig(; active=true),
            escape=EscapeConfig(; active=true),
            atmosphere=AtmosphereConfig(; active=true),
            hydrothermal=HydrothermalConfig(; active=true),
            accretion=AccretionConfig(; active=true),
            telescoping=TelescopingConfig(; active=true),
            refractory=RefractoryConfig(; active=true),
            volatile_mixture=VolatileMixtureConfig(; active=true),
            magma_transport=MagmaTransportConfig(; active=true),
            redox=RedoxConfig(; active=true),
            magma_degassing=MagmaOceanDegassingConfig(; active=true),
            disk=base_cfg.disk,
        )
        @test SphericalConductionValidation.canonical_config_sha256(unrelated_cfg) ==
            base_sha

        mock_extended_cfg = (;
            grid=base_cfg.grid,
            geometry=base_cfg.geometry,
            time=base_cfg.time,
            thermodynamics=base_cfg.thermodynamics,
            materials=base_cfg.materials,
            poroelasticity=base_cfg.poroelasticity,
            disk=base_cfg.disk,
            unrelated_new_section=(; extra_new_field=42, another_field="ignored"),
        )
        @test SphericalConductionValidation.canonical_config_sha256(mock_extended_cfg) ==
            base_sha

        # Reference literal constants must match SimulationConfig defaults, and fail when a default changes.
        ok_lit, errs_lit = SphericalConductionValidation.check_reference_literals_match_config(
            base_cfg
        )
        @test ok_lit
        @test isempty(errs_lit)

        changed_default_cfg = SphericalConductionValidation.make_spherical_conduction_config(;
            Nx=65, f_al=2.0e23, k_rock=3.5
        )
        ok_bad_lit, errs_bad_lit = SphericalConductionValidation.check_reference_literals_match_config(
            changed_default_cfg
        )
        @test !ok_bad_lit
        @test any(occursin("F_AL", e) for e in errs_bad_lit)
        @test any(occursin("K_COND", e) for e in errs_bad_lit)
    end

    @testset "Analytical and 1D finite-volume reference solutions" begin
        # Input domain contracts.
        @test_throws DomainError SphericalConductionValidation.analytical_T_BE(-1.0)
        @test_throws DomainError SphericalConductionValidation.analytical_T_BE(0.0; R_b=0.0)
        @test_throws DomainError SphericalConductionValidation.analytical_T_BE(
            0.0; dt_s=0.0
        )
        @test_throws DomainError SphericalConductionValidation.analytical_T_BE(
            0.0; n_steps=0
        )
        @test_throws DomainError SphericalConductionValidation.analytical_T_BE(
            0.0; n_modes=0
        )
        @test_throws DomainError SphericalConductionValidation.analytical_T_continuous(-1.0)
        @test_throws DomainError SphericalConductionValidation.analytical_T_continuous(
            0.0; R_b=-10.0
        )
        @test_throws DomainError SphericalConductionValidation.analytical_T_continuous(
            0.0; t_tot_s=0.0
        )
        @test_throws DomainError SphericalConductionValidation.analytical_T_continuous(
            0.0; n_modes=0
        )
        @test_throws DomainError SphericalConductionValidation.solve_1d_spherical_fv(1)
        @test_throws DomainError SphericalConductionValidation.solve_1d_spherical_fv(
            10; dt_s=-1.0
        )
        @test_throws DomainError SphericalConductionValidation.solve_1d_spherical_fv(
            10; n_steps=0
        )
        @test_throws DomainError SphericalConductionValidation.make_spherical_conduction_config(;
            Nx=3
        )
        @test_throws DomainError SphericalConductionValidation.run_steady_spherical_conduction(
            3
        )
        @test_throws DomainError SphericalConductionValidation.run_transient_spherical_conduction(
            3
        )

        # Boundary condition at r = R_planet and radial monotonicity for r in [0, R_planet].
        T_surf_be = SphericalConductionValidation.analytical_T_BE(
            SphericalConductionValidation.R_PLANET
        )
        T_surf_cont = SphericalConductionValidation.analytical_T_continuous(
            SphericalConductionValidation.R_PLANET
        )
        @test isapprox(T_surf_be, SphericalConductionValidation.T_SURF; atol=1.0e-12)
        @test isapprox(T_surf_cont, SphericalConductionValidation.T_SURF; atol=1.0e-12)

        radii = [0.0, 10_000.0, 25_000.0, 40_000.0, 50_000.0]
        T_profile = [SphericalConductionValidation.analytical_T_BE(r) for r in radii]
        T_cont_profile = [
            SphericalConductionValidation.analytical_T_continuous(r) for r in radii
        ]
        @test all(diff(T_profile) .< 0.0)
        @test all(diff(T_cont_profile) .< 0.0)
        @test isapprox(T_cont_profile[1], T_profile[1]; rtol=0.05)

        # Removable resonance limits (m = 1: mu_1 = LAMBDA_AL and p = q) compared against R * (1 + 1e-6).
        R_res_cont =
            pi * sqrt(
                SphericalConductionValidation.KAPPA /
                SphericalConductionValidation.LAMBDA_AL,
            )
        T_cont_res = SphericalConductionValidation.analytical_T_continuous(
            0.0; R_b=R_res_cont, n_modes=10
        )
        T_cont_nbr = SphericalConductionValidation.analytical_T_continuous(
            0.0; R_b=R_res_cont * (1.0 + 1.0e-6), n_modes=10
        )
        @test isfinite(T_cont_res)
        @test isapprox(T_cont_res, T_cont_nbr; rtol=1.0e-4)

        q_be = exp(
            -SphericalConductionValidation.LAMBDA_AL * SphericalConductionValidation.DT_S
        )
        R_res_be =
            pi * sqrt(
                SphericalConductionValidation.KAPPA * SphericalConductionValidation.DT_S /
                (1.0 / q_be - 1.0),
            )
        T_be_res = SphericalConductionValidation.analytical_T_BE(
            0.0; R_b=R_res_be, n_modes=10
        )
        T_be_nbr = SphericalConductionValidation.analytical_T_BE(
            0.0; R_b=R_res_be * (1.0 + 1.0e-6), n_modes=10
        )
        @test isfinite(T_be_res)
        @test isapprox(T_be_res, T_be_nbr; rtol=1.0e-4)

        # Benchmark case must stay well away from the p = q resonance cliff across all 500 modes.
        min_pq_rel = minimum(1:500) do m
            mu_m =
                SphericalConductionValidation.KAPPA *
                (m * pi / SphericalConductionValidation.R_PLANET)^2
            p_m = 1.0 / (1.0 + mu_m * SphericalConductionValidation.DT_S)
            return abs(p_m - q_be) / q_be
        end
        min_pq_separation_tol = 1.0e-6
        @test min_pq_rel >= min_pq_separation_tol

        # Mode truncation convergence: 500 modes vs 1000 modes <= 1e-7.
        dT_500 =
            SphericalConductionValidation.analytical_T_BE(0.0; n_modes=500) -
            SphericalConductionValidation.T_SURF
        dT_1000 =
            SphericalConductionValidation.analytical_T_BE(0.0; n_modes=1000) -
            SphericalConductionValidation.T_SURF
        mode_trunc_rel = abs(dT_500 - dT_1000) / dT_1000
        @test mode_trunc_rel <= 1.0e-7

        # Cross-verification of analytical_T_BE(r_fv[1]) against 1D spherical finite-volume solver at r_fv[1].
        r_fv, T_fv = SphericalConductionValidation.solve_1d_spherical_fv(500)
        dT_be_r1 =
            SphericalConductionValidation.analytical_T_BE(r_fv[1]; n_modes=500) -
            SphericalConductionValidation.T_SURF
        dT_fv_r1 = T_fv[1] - SphericalConductionValidation.T_SURF
        be_vs_1dfv_rel = abs(dT_be_r1 - dT_fv_r1) / dT_be_r1
        @test be_vs_1dfv_rel <= 1.0e-5

        # Reference sensitivity test (benchmark-only): omitting step-average factor breaks agreement by >= 2x.
        dT_be_no_step_avg =
            SphericalConductionValidation.analytical_T_BE(
                r_fv[1]; n_modes=500, include_step_average=false
            ) - SphericalConductionValidation.T_SURF
        broken_ref_rel = abs(dT_be_no_step_avg - dT_fv_r1) / dT_fv_r1
        min_broken_ref_diff = 0.10
        @test broken_ref_rel >= min_broken_ref_diff
    end

    bench = SphericalConductionValidation.evaluate_spherical_conduction_benchmark()

    @testset "Steady spherical Poisson conduction" begin
        s33 = bench.s33
        s65 = bench.s65

        @test isapprox(s33.ratio_3d_2d, 2.0 / 3.0; atol=0.01, rtol=0.015)
        @test isapprox(s65.ratio_3d_2d, 2.0 / 3.0; atol=0.01, rtol=0.015)
        @test s33.l2_rel <= 0.15
        @test s65.l2_rel <= 0.08
        @test s65.l2_rel < s33.l2_rel
        @test s65.energy_res <= 1.0e-11

        # Determinism repeat check: repeating baseline steady call reproduces metrics identically.
        s33_repeat = SphericalConductionValidation.run_steady_spherical_conduction(
            33;
            spherical_metric=true,
            metric_scale=1.0,
            clamp_exterior=true,
            perturb_stencil=false,
        )
        @test isapprox(s33_repeat.ratio_3d_2d, s33.ratio_3d_2d; rtol=1.0e-12)
        @test isapprox(s33_repeat.l2_rel, s33.l2_rel; rtol=1.0e-12)

        # Sensitivity test on solver physics path: disabling spherical_metric breaks ratio and L2 guards by >= 2x.
        min_ratio_break_65 = 2.0 * 0.01
        min_l2_break_65 = 2.0 * 0.08
        min_energy_break_65 = 2.0 * 1.0e-11

        s65_no_metric = SphericalConductionValidation.run_steady_spherical_conduction(
            65; spherical_metric=false
        )
        @test abs(s65_no_metric.ratio_3d_2d - 2.0 / 3.0) >= min_ratio_break_65
        @test s65_no_metric.l2_rel >= min_l2_break_65

        # Benchmark-only sensitivity tests: exterior BC clamping, stencil diagonal perturbation, and source scaling.
        s65_unclamped = SphericalConductionValidation.run_steady_spherical_conduction(
            65; clamp_exterior=false
        )
        @test s65_unclamped.l2_rel >= min_l2_break_65

        s65_perturbed = SphericalConductionValidation.run_steady_spherical_conduction(
            65; perturb_stencil=true
        )
        @test s65_perturbed.energy_res >= min_energy_break_65

        s65_doubled_source = SphericalConductionValidation.run_steady_spherical_conduction(
            65; metric_scale=2.0
        )
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
        @test isapprox(t65.T_max_air_shell, SphericalConductionValidation.T_SURF; atol=0.01)

        # Determinism repeat check: repeating baseline transient call reproduces metrics identically.
        t33_repeat = SphericalConductionValidation.run_transient_spherical_conduction(
            33; dsubgridt=0.0
        )
        @test isapprox(t33_repeat.dT_center, t33.dT_center; rtol=1.0e-12)
        @test isapprox(t33_repeat.l2_rel, t33.l2_rel; rtol=1.0e-12)

        # Sensitivity test 1 (simulation physics path): scaling simulated f_al by 1.02 with reference unchanged
        # must fail the machine-independent physics regression band (rtol=0.01 on dT, atol=0.005 on L2) by >= 2x.
        t65_f_al_102 = SphericalConductionValidation.run_transient_spherical_conduction(
            65; f_al=1.02 * SphericalConductionValidation.F_AL
        )
        metrics_f_al_102 = copy(bench.metrics)
        metrics_f_al_102["transient_center_dT_nx65"] = t65_f_al_102.dT_center
        metrics_f_al_102["transient_profile_l2_rel_error_nx65"] = t65_f_al_102.l2_rel
        ok_band_102, errs_band_102 = SphericalConductionValidation.check_physics_regression_bands(
            metrics_f_al_102, record_path
        )
        @test !ok_band_102
        @test any(occursin("transient_center_dT_nx65", e) for e in errs_band_102)
        @test any(occursin("transient_profile_l2_rel_error_nx65", e) for e in errs_band_102)
        dT_band_tol = 0.01 * t65.dT_center
        dT_102_margin = abs(t65_f_al_102.dT_center - t65.dT_center) / dT_band_tol
        l2_102_margin = abs(t65_f_al_102.l2_rel - t65.l2_rel) / 0.005
        min_dT_margin_factor = 2.0 - 1.0e-9
        min_l2_margin_factor = 2.0
        @test dT_102_margin >= min_dT_margin_factor
        @test l2_102_margin >= min_l2_margin_factor

        # Sensitivity test 2 (simulation physics path): disabling spherical metric breaks centre, L2, and convergence by >= 2x.
        t33_no_metric = SphericalConductionValidation.run_transient_spherical_conduction(
            33; spherical_metric=false
        )
        t65_no_metric = SphericalConductionValidation.run_transient_spherical_conduction(
            65; spherical_metric=false
        )
        min_center_rel_break = 2.0 * 0.05
        min_l2_break = 2.0 * 0.10
        min_conv_break = 2.0 * 0.60
        @test abs(t65_no_metric.dT_center - t65_no_metric.dT_be_rmin) /
              t65_no_metric.dT_be_rmin >= min_center_rel_break
        @test t65_no_metric.l2_rel >= min_l2_break
        @test (t65_no_metric.l2_rel / t33_no_metric.l2_rel) >= min_conv_break

        # Sensitivity test 3 (simulation physics path): disabling sticky-air reset breaks air-shell clamp by >= 2x.
        t65_no_reset = SphericalConductionValidation.run_transient_spherical_conduction(
            65; n_steps=5, disk_enabled=false
        )
        min_air_clamp_break = 2.0 * 0.01
        @test abs(t65_no_reset.T_max_air_shell - SphericalConductionValidation.T_SURF) >=
            min_air_clamp_break

        # Benchmark-only parameter sensitivity tests at Nx=33 (conductivity, air heat capacity, regularization).
        t33_zero_k = SphericalConductionValidation.run_transient_spherical_conduction(
            33; k_rock=1.0e-6, k_air=1.0e-6
        )
        @test abs(t33_zero_k.dT_center - t33_zero_k.dT_be_rmin) / t33_zero_k.dT_be_rmin >=
            min_center_rel_break
        @test t33_zero_k.l2_rel >= min_l2_break

        t33_low_air_cp = SphericalConductionValidation.run_transient_spherical_conduction(
            33; rhocp_air=3.3e6
        )
        @test abs(t33_low_air_cp.T_max_air_shell - SphericalConductionValidation.T_SURF) >=
            min_air_clamp_break
        @test t33_low_air_cp.l2_rel >= min_l2_break

        t33_overreg = SphericalConductionValidation.run_transient_spherical_conduction(
            33; reg_cells=15.0
        )
        @test abs(t33_overreg.dT_center - t33_overreg.dT_be_rmin) /
              t33_overreg.dT_be_rmin >= min_center_rel_break

        # Recompute every recorded reference from the reference functions and compare against the TOML record.
        raw_rec_check = TOML.parsefile(record_path)
        exp_refs = SphericalConductionValidation.expected_metric_references(bench)
        for m in raw_rec_check["metrics"]
            mid = String(m["id"])
            @test isapprox(
                Float64(m["reference"]), exp_refs[mid]; rtol=1.0e-12, atol=1.0e-12
            )
        end

        # Full benchmark evaluation, physics regression bands, and record reproducibility check.
        ok_bands, band_errs = SphericalConductionValidation.check_physics_regression_bands(
            bench.metrics, record_path
        )
        @test ok_bands
        @test isempty(band_errs)

        ok_repro, repro_errs = SphericalConductionValidation.check_validation_reproducibility(
            bench, record_path
        )
        @test ok_repro
        @test isempty(repro_errs)
        @test run_validation_check(record_path; bench=bench)

        # Mutating a recorded metric, reference, metric set, or config hash in the record must fail reproducibility.
        mktempdir() do tmpdir
            exported_path = joinpath(tmpdir, "exported_spherical_conduction.toml")
            src_commit_str = String(raw_rec_check["metadata"]["src_commit"])
            write_spherical_conduction_record(
                exported_path; src_commit=src_commit_str, bench=bench
            )
            @test isfile(exported_path)
            @test run_validation_check(exported_path; bench=bench)

            @test !run_validation_check(joinpath(tmpdir, "missing.toml"); bench=bench)
            @test !SphericalConductionValidation.check_validation_records(
                joinpath(tmpdir, "missing_dir")
            )
            empty_records_dir = joinpath(tmpdir, "empty_records")
            mkpath(empty_records_dir)
            @test !SphericalConductionValidation.check_validation_records(empty_records_dir)

            # Schema-invalid record file in directory and CLI entry points.
            invalid_schema_dir = joinpath(tmpdir, "invalid_records")
            mkpath(invalid_schema_dir)
            invalid_schema_path = joinpath(invalid_schema_dir, "bad.toml")
            write(invalid_schema_path, "[metadata]\nid = \"bad\"\n")
            @test !SphericalConductionValidation.check_validation_records(
                invalid_schema_dir
            )
            @test !run_validation_check(invalid_schema_path; bench=bench)
            @test run_validation_cli(["--check", record_path]; bench=bench) == 0
            @test run_validation_cli(["--unknown-flag"]; bench=bench) == 2
            @test run_validation_cli([invalid_schema_path]; bench=bench) == 1

            # Drift within schema tolerance must be caught by recorded-value reproducibility.
            tampered_path = joinpath(tmpdir, "spherical_conduction.toml")
            raw_rec = TOML.parsefile(record_path)
            raw_rec["metrics"][1]["value"] += 0.001
            open(tampered_path, "w") do io
                TOML.print(io, raw_rec)
            end
            ok_tampered, errs_tampered = SphericalConductionValidation.check_validation_reproducibility(
                bench, tampered_path; check_recorded_values=true
            )
            @test !ok_tampered
            @test any(occursin("differs from recorded value", e) for e in errs_tampered)
            @test !run_validation_check(tampered_path; bench=bench)

            # Removing a metric from the record must fail the exact metric ID set check.
            dropped_metric_path = joinpath(tmpdir, "dropped_metric.toml")
            raw_drop = TOML.parsefile(record_path)
            pop!(raw_drop["metrics"])
            open(dropped_metric_path, "w") do io
                TOML.print(io, raw_drop)
            end
            ok_drop, errs_drop = SphericalConductionValidation.check_validation_reproducibility(
                bench, dropped_metric_path
            )
            @test !ok_drop
            @test any(occursin("metric IDs", e) for e in errs_drop)

            bad_cfg_path = joinpath(tmpdir, "spherical_conduction_cfg.toml")
            raw_cfg_rec = TOML.parsefile(record_path)
            raw_cfg_rec["metadata"]["config_sha256"] = repeat("0", 64)
            open(bad_cfg_path, "w") do io
                TOML.print(io, raw_cfg_rec)
            end
            ok_bad_cfg, errs_bad_cfg = SphericalConductionValidation.check_validation_reproducibility(
                bench, bad_cfg_path
            )
            @test !ok_bad_cfg
            @test !isempty(errs_bad_cfg)
            @test !run_validation_check(bad_cfg_path; bench=bench)
        end
    end
end

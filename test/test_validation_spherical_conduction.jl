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

        bad_le_neg = deepcopy(raw_record)
        bad_le_neg["metrics"][3]["value"] = -0.01
        ok_le_neg, errs_le_neg = SphericalConductionValidation.validate_record_dict(
            bad_le_neg; require_doc_pages=false
        )
        @test !ok_le_neg
        @test any(occursin("fails comparison 'le_atol'", e) for e in errs_le_neg)

        bad_type_rec = deepcopy(raw_record)
        bad_type_rec["metrics"][1]["value"] = "not_a_number"
        ok_bad_type, errs_bad_type = SphericalConductionValidation.validate_record_dict(
            bad_type_rec; require_doc_pages=false
        )
        @test !ok_bad_type
        @test any(occursin("numeric", e) for e in errs_bad_type)

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

        # Canonical config SHA-256: hashes the 9 benchmark physics sections and steady inputs,
        # without depending on container/SVector type printing.
        @test SphericalConductionValidation.BENCHMARK_CONFIG_SECTIONS == (
            :grid,
            :geometry,
            :time,
            :solver,
            :thermodynamics,
            :materials,
            :poroelasticity,
            :coreformation,
            :disk,
        )
        @test SphericalConductionValidation._canonical_field_repr(
            SVector{3,Float64}(1.0, 2.5, 3.0)
        ) == "[1.0, 2.5, 3.0]"
        @test SphericalConductionValidation._canonical_field_repr(
            Dict(:b => 2.0, :a => 1)
        ) == "{:a => 1, :b => 2.0}"

        base_cfg = SphericalConductionValidation.make_spherical_conduction_config(; Nx=65)
        base_sha = SphericalConductionValidation.canonical_config_sha256(base_cfg)
        @test base_sha == raw_record["metadata"]["config_sha256"]

        altered_cfg = SphericalConductionValidation.make_spherical_conduction_config(;
            Nx=65, k_rock=3.5
        )
        @test SphericalConductionValidation.canonical_config_sha256(altered_cfg) != base_sha
        seed_cfg = SphericalConductionValidation.make_spherical_conduction_config(;
            Nx=65, seed=99999
        )
        @test SphericalConductionValidation.canonical_config_sha256(seed_cfg) != base_sha
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

        # Unrelated sections (such as volatiles, reaction, output) or extra fields in unrelated sections must not change the hash.
        unrelated_cfg = SimulationConfig(;
            grid=base_cfg.grid,
            geometry=base_cfg.geometry,
            time=base_cfg.time,
            solver=base_cfg.solver,
            thermodynamics=base_cfg.thermodynamics,
            materials=base_cfg.materials,
            poroelasticity=base_cfg.poroelasticity,
            output=OutputConfig(; savematstep=1, visstep=1, save_final=true),
            reaction=ReactionConfig(; active=true),
            melting=MeltingConfig(; active=true),
            coreformation=base_cfg.coreformation,
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
            solver=base_cfg.solver,
            thermodynamics=base_cfg.thermodynamics,
            materials=base_cfg.materials,
            poroelasticity=base_cfg.poroelasticity,
            coreformation=base_cfg.coreformation,
            disk=base_cfg.disk,
            unrelated_new_section=(; extra_new_field=42, another_field="ignored"),
        )
        @test SphericalConductionValidation.canonical_config_sha256(mock_extended_cfg) ==
            base_sha

        # Reference literal constants must match non-builder SimulationConfig defaults and X_FE_REF_CHONDRITE.
        ok_lit, errs_lit = SphericalConductionValidation.check_reference_literals_match_config(
            base_cfg
        )
        @test ok_lit
        @test isempty(errs_lit)

        changed_default_cfg = SimulationConfig(;
            thermodynamics=ThermalConfig(; ratio_al=6.0e-5),
            coreformation=CoreFormationConfig(; rho_metal=6000.0),
        )
        ok_bad_lit, errs_bad_lit = SphericalConductionValidation.check_reference_literals_match_config(
            changed_default_cfg; x_fe_ref=0.19
        )
        @test !ok_bad_lit
        @test length(errs_bad_lit) == 3
        @test any(occursin("RATIO_AL", e) for e in errs_bad_lit)
        @test any(occursin("RHO_MET", e) for e in errs_bad_lit)
        @test any(occursin("X_FE", e) for e in errs_bad_lit)
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
        @test_throws DomainError SphericalConductionValidation.make_spherical_conduction_config(;
            n_steps=0
        )
        @test_throws DomainError SphericalConductionValidation.make_spherical_conduction_config(;
            dt_yr=-1.0
        )
        @test_throws DomainError SphericalConductionValidation.make_spherical_conduction_config(;
            reg_cells=-0.1
        )
        @test_throws DomainError SphericalConductionValidation.make_spherical_conduction_config(;
            k_rock=0.0
        )
        @test_throws DomainError SphericalConductionValidation.make_spherical_conduction_config(;
            k_air=0.0
        )
        @test_throws DomainError SphericalConductionValidation.make_spherical_conduction_config(;
            rhocp_air=0.0
        )
        @test_throws DomainError SphericalConductionValidation.make_spherical_conduction_config(;
            f_al=0.0
        )
        cfg_dom = SphericalConductionValidation.make_spherical_conduction_config(; Nx=33)
        @test_throws DomainError SphericalConductionValidation.canonical_config_sha256(
            cfg_dom; H_steady=0.0
        )
        @test_throws DomainError SphericalConductionValidation.canonical_config_sha256(
            cfg_dom; L_box=0.0
        )
        @test_throws DomainError SphericalConductionValidation.canonical_config_sha256(
            cfg_dom; reg_cells=-1.0
        )
        @test_throws DomainError SphericalConductionValidation.canonical_config_sha256(
            cfg_dom; k_cond=0.0
        )
        @test_throws DomainError SphericalConductionValidation.check_physics_regression_bands(
            Dict{String,Float64}(), ""
        )
        @test_throws DomainError SphericalConductionValidation.check_validation_reproducibility(
            nothing, ""
        )
        @test_throws DomainError run_validation_check("")
        @test_throws DomainError write_spherical_conduction_record("")
        @test_throws DomainError write_spherical_conduction_record(
            "bad.toml"; src_commit="not_hex_sha"
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
        # Real guard tolerance on ratio_3d_2d is atol=0.01 + rtol=0.015*(2/3) = 0.02, so 2x is 0.04.
        min_ratio_break_65 = 2.0 * (0.01 + 0.015 * (2.0 / 3.0))
        min_l2_break_65 = 2.0 * 0.08
        min_energy_break_65 = 2.0 * 1.0e-11

        s65_no_metric = SphericalConductionValidation.run_steady_spherical_conduction(
            65; spherical_metric=false
        )
        @test abs(s65_no_metric.ratio_3d_2d - 2.0 / 3.0) >= min_ratio_break_65
        @test s65_no_metric.l2_rel >= min_l2_break_65

        # Benchmark-only sensitivity tests: exterior BC clamping (Nx=65), stencil diagonal perturbation, and source scaling.
        s65_unclamped = SphericalConductionValidation.run_steady_spherical_conduction(
            65; clamp_exterior=false
        )
        @test s65_unclamped.l2_rel >= min_l2_break_65

        s33_perturbed = SphericalConductionValidation.run_steady_spherical_conduction(
            33; perturb_stencil=true
        )
        @test s33_perturbed.energy_res >= min_energy_break_65

        s33_doubled_source = SphericalConductionValidation.run_steady_spherical_conduction(
            33; metric_scale=2.0
        )
        @test abs(s33_doubled_source.ratio_3d_2d - 2.0 / 3.0) >= min_ratio_break_65
        @test s33_doubled_source.l2_rel >= 2.0 * 0.15
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

        # Determinism repeat check: repeating transient call reproduces metrics identically.
        t11_rep1 = SphericalConductionValidation.run_transient_spherical_conduction(
            11; n_steps=1, dsubgridt=0.0
        )
        t11_rep2 = SphericalConductionValidation.run_transient_spherical_conduction(
            11; n_steps=1, dsubgridt=0.0
        )
        @test isapprox(t11_rep2.dT_center, t11_rep1.dT_center; rtol=1.0e-12)
        @test isapprox(t11_rep2.l2_rel, t11_rep1.l2_rel; rtol=1.0e-12)

        # Sensitivity test 1 (simulation physics path): scaling simulated f_al by 1.02 with reference unchanged
        # must fail the physics regression band (rtol=0.01 on dT, atol=0.005 on L2).
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
        l2_102_margin = abs(t65_f_al_102.l2_rel - t65.l2_rel) / 0.005
        min_l2_margin_factor = 2.0
        @test l2_102_margin >= min_l2_margin_factor

        # Sensitivity test 2 (simulation physics path): disabling spherical metric at Nx=33 and Nx=65
        # breaks L2 (2.05x of 0.10), convergence l2(65)/l2(33) (2.73x of 0.60), and centre dT (3.53x of rtol=0.05,
        # 1.78x of combined atol=25 K + rtol=0.05*ref = 50.35 K).
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
        t33_no_reset = SphericalConductionValidation.run_transient_spherical_conduction(
            33; n_steps=14, disk_enabled=false
        )
        min_air_clamp_break = 2.0 * 0.01
        @test abs(t33_no_reset.T_max_air_shell - SphericalConductionValidation.T_SURF) >=
            min_air_clamp_break

        # Additional simulation parameter sensitivity tests (conductivity, air heat capacity, regularization).
        t17_zero_k = SphericalConductionValidation.run_transient_spherical_conduction(
            17; k_rock=1.0e-6, k_air=1.0e-6
        )
        @test abs(t17_zero_k.dT_center - t17_zero_k.dT_be_rmin) / t17_zero_k.dT_be_rmin >=
            min_center_rel_break
        @test t17_zero_k.l2_rel >= min_l2_break

        t17_low_air_cp = SphericalConductionValidation.run_transient_spherical_conduction(
            17; rhocp_air=3.3e6
        )
        @test abs(t17_low_air_cp.T_max_air_shell - SphericalConductionValidation.T_SURF) >=
            min_air_clamp_break

        t21_overreg = SphericalConductionValidation.run_transient_spherical_conduction(
            21; reg_cells=15.0
        )
        @test abs(t21_overreg.dT_center - t21_overreg.dT_be_rmin) /
              t21_overreg.dT_be_rmin >= min_center_rel_break

        # Alternative marker RNG seed (seed = 2) must stay inside every physics regression band.
        t33_seed2 = SphericalConductionValidation.run_transient_spherical_conduction(
            33; seed=2
        )
        t65_seed2 = SphericalConductionValidation.run_transient_spherical_conduction(
            65; seed=2
        )
        metrics_seed2 = copy(bench.metrics)
        metrics_seed2["transient_center_dT_nx65"] = t65_seed2.dT_center
        metrics_seed2["transient_profile_l2_rel_error_nx65"] = t65_seed2.l2_rel
        metrics_seed2["transient_profile_l2_convergence_ratio_65_33"] =
            t65_seed2.l2_rel / t33_seed2.l2_rel
        metrics_seed2["transient_max_air_shell_temperature_nx65"] =
            t65_seed2.T_max_air_shell
        ok_seed2, errs_seed2 = SphericalConductionValidation.check_physics_regression_bands(
            metrics_seed2, record_path
        )
        @test ok_seed2
        @test isempty(errs_seed2)

        # Recompute every recorded reference from the reference functions and compare against the TOML record.
        raw_rec_check = TOML.parsefile(record_path)
        @test raw_rec_check["config"]["Nx_coarse"] ==
            SphericalConductionValidation.NX_COARSE
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

        # Mutating a recorded metric, reference, metric set, config hash, Julia version, or export guard must behave as specified.
        mktempdir() do tmpdir
            exported_path = joinpath(tmpdir, "exported_spherical_conduction.toml")
            src_commit_str = String(raw_rec_check["metadata"]["src_commit"])
            write_spherical_conduction_record(
                exported_path; src_commit=src_commit_str, bench=bench
            )
            @test isfile(exported_path)
            @test run_validation_check(exported_path; bench=bench)

            # Export guard (T6): refuses to overwrite when a metric leaves its physics band around the existing record
            # unless --accept-physics-change is given, and refuses when validate_record_dict fails.
            bench_band_shift = (; bench..., metrics=copy(bench.metrics))
            bench_band_shift.metrics["transient_center_dT_nx65"] =
                bench.metrics["transient_center_dT_nx65"] * 1.02
            @test_throws ErrorException write_spherical_conduction_record(
                exported_path;
                src_commit=src_commit_str,
                bench=bench_band_shift,
                accept_physics_change=false,
            )
            @test run_export_cli(
                [exported_path];
                src_commit=src_commit_str,
                bench=bench_band_shift,
                check_clean_tree=false,
            ) == 1
            rec_after_refusal = TOML.parsefile(exported_path)
            @test isapprox(
                Float64(rec_after_refusal["metrics"][7]["value"]),
                bench.metrics["transient_center_dT_nx65"];
                rtol=1.0e-14,
            )

            @test run_export_cli(
                ["--accept-physics-change", exported_path];
                src_commit=src_commit_str,
                bench=bench_band_shift,
                check_clean_tree=false,
            ) == 0
            rec_after_accept = TOML.parsefile(exported_path)
            @test isapprox(
                Float64(rec_after_accept["metrics"][7]["value"]),
                bench_band_shift.metrics["transient_center_dT_nx65"];
                rtol=1.0e-14,
            )

            # Failing metric (outside reference tolerance) must be refused even with --accept-physics-change.
            bench_fail_ref = (; bench..., metrics=copy(bench.metrics))
            bench_fail_ref.metrics["steady_center_dT_ratio_3d_2d_nx33"] = 1.0
            @test run_export_cli(
                ["--accept-physics-change", exported_path];
                src_commit=src_commit_str,
                bench=bench_fail_ref,
                check_clean_tree=false,
            ) == 1
            @test run_export_cli(
                ["--unknown-flag"];
                src_commit=src_commit_str,
                bench=bench,
                check_clean_tree=false,
            ) == 2

            # Restore baseline exported_path for remaining checks.
            write_spherical_conduction_record(
                exported_path;
                src_commit=src_commit_str,
                bench=bench,
                accept_physics_change=true,
            )

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

            # Drift within schema tolerance must be caught by recorded-value reproducibility on same Julia major.minor.
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

            # On a different Julia major.minor, replay compare is skipped while physics bands still run.
            diff_minor_path = joinpath(tmpdir, "diff_julia_minor.toml")
            raw_diff_minor = TOML.parsefile(record_path)
            raw_diff_minor["metadata"]["julia_version"] = "1.99.0"
            raw_diff_minor["metrics"][1]["value"] += 0.001
            open(diff_minor_path, "w") do io
                TOML.print(io, raw_diff_minor)
            end
            ok_diff_minor, errs_diff_minor = SphericalConductionValidation.check_validation_reproducibility(
                bench, diff_minor_path; check_recorded_values=true
            )
            @test ok_diff_minor
            @test isempty(errs_diff_minor)

            raw_diff_minor["metrics"][1]["value"] += 0.02
            open(diff_minor_path, "w") do io
                TOML.print(io, raw_diff_minor)
            end
            ok_diff_minor_fail, errs_diff_minor_fail = SphericalConductionValidation.check_validation_reproducibility(
                bench, diff_minor_path; check_recorded_values=true
            )
            @test !ok_diff_minor_fail
            @test any(
                occursin("outside physics regression band", e) for e in errs_diff_minor_fail
            )

            # Changing Nx_coarse in the record must fail reproducibility (L-g).
            bad_nx_coarse_path = joinpath(tmpdir, "bad_nx_coarse.toml")
            raw_nx_coarse = TOML.parsefile(record_path)
            raw_nx_coarse["config"]["Nx_coarse"] = 17
            open(bad_nx_coarse_path, "w") do io
                TOML.print(io, raw_nx_coarse)
            end
            ok_nx_coarse, errs_nx_coarse = SphericalConductionValidation.check_validation_reproducibility(
                bench, bad_nx_coarse_path
            )
            @test !ok_nx_coarse
            @test any(occursin("config.Nx_coarse", e) for e in errs_nx_coarse)

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

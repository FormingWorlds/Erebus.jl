using TOML
using Printf

const DEFAULT_RECORDS_DIR = normpath(joinpath(@__DIR__, "..", "benchmarks", "records"))
const VALIDATION_DOCS_DIR = normpath(joinpath(@__DIR__, "..", "docs", "src"))

const REQUIRED_METADATA_KEYS = (
    "id",
    "title",
    "src_commit",
    "timestamp_utc",
    "julia_version",
    "config_sha256",
    "reference_type",
    "reference_citation",
)

# The CI record of the melting shell and its suite record at Nx 129 share the keys.
const MELTING_SHELL_CONFIG_KEYS = (
    "Nx",
    "Nr_reference",
    "n_steps",
    "dt_yr",
    "t_form_yr",
    "phi_init",
    "porosity_model",
    "R_planet_m",
    "L_box_m",
    "ratio_al",
    "L_melt_J_kg",
    "T_solidus_K",
    "T_liquidus_K",
    "rho_melt_kg_m3",
    "rhocp_rock_J_m3K",
    "k_rock_W_mK",
    "T_space_K",
    "surface_temperature",
    "tmax_l2_T_min_K",
    "melt_radius_atol_cells",
    "seed",
    "spherical_metric",
)

const REQUIRED_CONFIG_KEYS = Dict{String,NTuple{N,String} where N}(
    "spherical_conduction" => (
        "Nx_coarse",
        "Nx_fine",
        "n_steps",
        "dt_yr",
        "R_planet_m",
        "L_box_m",
        "k_W_mK",
        "rho_kg_m3",
        "cp_J_kgK",
        "T_surf_K",
        "H_steady_W_m3",
        "spherical_metric",
        "reg_cells",
        "disk_enabled",
    ),
    "spherical_compaction" => (
        "Nx",
        "R_planet_m",
        "L_box_m",
        "mu_f_Pa_s",
        "k_m2",
        "phi",
        "delta_m",
        "zeta_Pa_s",
        "eta_ratio",
        "eta_air_Pa_s",
        "Pe_R_Pa",
        "k_vent_m2",
        "T_surf_K",
        "P_amb_Pa",
        "dt_yr",
        "inner_fraction",
        "spherical_metric",
        "darcy_elimination",
    ),
    "spherical_darcy_flow" => (
        "Nx",
        "R_planet_m",
        "L_box_m",
        "mu_f_Pa_s",
        "k_m2",
        "k_core_m2",
        "a_jump_m",
        "S_source_per_s",
        "phi",
        "eta_solid_Pa_s",
        "k_vent_m2",
        "T_surf_K",
        "P_amb_Pa",
        "inner_fraction",
        "spherical_metric",
        "darcy_elimination",
    ),
    "spherical_pressure_diffusion" => (
        "Nx",
        "R_planet_m",
        "L_box_m",
        "mu_f_Pa_s",
        "k_m2",
        "phi",
        "beta_phi_per_Pa",
        "beta_s_per_Pa",
        "beta_s_alt_per_Pa",
        "beta_f_per_Pa",
        "storage_per_Pa",
        "kbw",
        "P0_Pa",
        "P_amb_Pa",
        "k_vent_m2",
        "T_surf_K",
        "eta_rock_Pa_s",
        "shear_modulus_Pa",
        "tau_end",
        "n_steps",
        "dt_s",
        "inner_fraction",
        "spherical_metric",
        "darcy_elimination",
    ),
    "spherical_momentum" => (
        "Nx",
        "R_planet_m",
        "L_box_m",
        "eta_Pa_s",
        "eta_zeta_ratio",
        "phi",
        "k_m2",
        "mu_f_Pa_s",
        "rho_kg_m3",
        "rho_f_kg_m3",
        "shear_modulus_Pa",
        "F0_N_m3",
        "U0_m_s",
        "b_power",
        "c_power",
        "eta_zeta_ratio_c",
        "P_amb_c_Pa",
        "dt_yr",
        "P_amb_Pa",
        "fit_rho_min",
        "fit_rho_max",
        "fit_degree",
        "spherical_metric",
        "darcy_elimination",
    ),
    "porous_shell" => (
        "Nx",
        "Nr_reference",
        "n_steps",
        "dt_yr",
        "t_form_yr",
        "phi_init",
        "porosity_model",
        "R_planet_m",
        "L_box_m",
        "ratio_al",
        "k_rock_W_mK",
        "k_air_W_mK",
        "rhocp_rock_J_m3K",
        "rhocp_air_J_m3K",
        "T_space_K",
        "shell_dphi",
        "shell_depth_rel_bound",
        "seed",
        "spherical_metric",
    ),
    "porous_shell_suite" => (
        "tags",
        "Nx",
        "t_form_myr",
        "phi_init",
        "n_steps",
        "air_matched",
        "k_air_W_mK",
        "dt_yr",
        "Nr_reference",
        "R_planet_m",
        "L_box_m",
        "seed",
        "porosity_model",
        "shell_dphi",
        "shell_depth_rel_bound",
        "tmax_rel_l2_limit",
        "tmax_nx_min",
    ),
    "melting_shell" => MELTING_SHELL_CONFIG_KEYS,
    "melting_shell_suite" => MELTING_SHELL_CONFIG_KEYS,
)

const REQUIRED_METRIC_KEYS = (
    "id",
    "case",
    "description",
    "nx",
    "value",
    "reference",
    "atol",
    "rtol",
    "comparison",
    "units",
    "status",
)

"""
    render_validation_markdown_block(data::Dict{String,Any}) -> String

Render a deterministic Markdown summary table for a validation record.

Parameters
----------
data : Dict{String,Any}
    Parsed validation record dictionary containing `metadata` and `metrics`.

Returns
-------
String
    Fenced Markdown block bounded by `VALIDATION_RECORD_START` and `VALIDATION_RECORD_END`.
"""
function render_validation_markdown_block(data::Dict{String,Any})::String
    meta = data["metadata"]
    rec_id = String(meta["id"])
    commit_short = first(String(meta["src_commit"]), 8)
    cfg_short = first(String(meta["config_sha256"]), 12)
    lines = String[
        "```@raw html",
        "<!-- VALIDATION_RECORD_START: $(rec_id) -->",
        "```",
        "",
        "**Validation record**: `benchmarks/records/$(rec_id).toml` (Source commit `$(commit_short)`, config SHA-256 `$(cfg_short)`)",
        "",
        "| Case | Metric | Grid `Nx` | Measured | Reference | Tolerance | Units | Status |",
        "| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: |",
    ]
    for m in data["metrics"]
        case_str = String(m["case"])
        desc_str = String(m["description"])
        nx_str = Int(m["nx"]) > 0 ? string(Int(m["nx"])) : "1D"
        val_str = @sprintf("%.4e", Float64(m["value"]))
        ref_str = @sprintf("%.4e", Float64(m["reference"]))
        comp = String(m["comparison"])
        tol_str = if comp == "le_atol"
            @sprintf("<= %.2e", Float64(m["atol"]))
        else
            @sprintf("atol=%.2e, rtol=%.2e", Float64(m["atol"]), Float64(m["rtol"]))
        end
        unit_str = isempty(String(m["units"])) ? "dimensionless" : String(m["units"])
        stat_str = String(m["status"])
        push!(
            lines,
            "| `$(case_str)` | $(desc_str) | $(nx_str) | `$(val_str)` | `$(ref_str)` | `$(tol_str)` | $(unit_str) | **$(stat_str)** |",
        )
    end
    append!(lines, ["", "```@raw html", "<!-- VALIDATION_RECORD_END: $(rec_id) -->", "```"])
    return join(lines, "\n")
end

"""
    update_record_doc_block!(data::Dict{String,Any}; docs_dir::String=VALIDATION_DOCS_DIR) -> Bool

Replace the record table on the documentation page of a validation record.

Parameters
----------
data : Dict{String,Any}
    Parsed TOML dictionary of a validation record.
docs_dir : String, optional
    Root directory of documentation markdown files.

Returns
-------
Bool
    `true` if the page holds the rendered table afterwards, `false` if the page or its markers
    do not exist.
"""
function update_record_doc_block!(
    data::Dict{String,Any}; docs_dir::String=VALIDATION_DOCS_DIR
)::Bool
    rec_id = String(data["metadata"]["id"])
    doc_path = joinpath(docs_dir, String(data["metadata"]["page"]))
    isfile(doc_path) || return false
    doc_text = read(doc_path, String)
    start_mark = "```@raw html\n<!-- VALIDATION_RECORD_START: $(rec_id) -->"
    end_mark = "<!-- VALIDATION_RECORD_END: $(rec_id) -->\n```"
    first_at = findfirst(start_mark, doc_text)
    last_at = findfirst(end_mark, doc_text)
    (first_at === nothing || last_at === nothing) && return false
    block = render_validation_markdown_block(data)
    write(
        doc_path,
        doc_text[1:prevind(doc_text, first(first_at))] *
        block *
        doc_text[nextind(doc_text, last(last_at)):end],
    )
    return true
end

"""
    validate_record_dict(data::Dict{String,Any}; docs_dir::String=VALIDATION_DOCS_DIR, require_doc_pages::Bool=false, check_doc_sync::Bool=true) -> Tuple{Bool,Vector{String}}

Validate the schema, metric tolerances, and optional documentation block of a validation record.
Metrics using `comparison = "le_atol"` represent upper-bound magnitude checks and must store
`reference = 0.0`. The required `[config]` keys come from `REQUIRED_CONFIG_KEYS` under
`metadata.id`; a record id without an entry is an error. The `REQUIRED_METADATA_KEYS`,
`metadata.page` and the metric keys `id`, `case`, `description`, `comparison`, `units` and
`status` must hold strings, and the metric key `nx` an integer; another TOML type is a
schema error.

Parameters
----------
data : Dict{String,Any}
    Parsed TOML dictionary of a validation record.
docs_dir : String, optional
    Root directory of documentation markdown files.
require_doc_pages : Bool, optional
    Whether a missing target documentation page is treated as an error.
check_doc_sync : Bool, optional
    Whether to compare the page with the record. A candidate record that is not yet the
    committed one passes `false`.

Returns
-------
Tuple{Bool,Vector{String}}
    Boolean pass flag and vector of diagnostic error messages.
"""
function validate_record_dict(
    data::Dict{String,Any};
    docs_dir::String=VALIDATION_DOCS_DIR,
    require_doc_pages::Bool=false,
    check_doc_sync::Bool=true,
)::Tuple{Bool,Vector{String}}
    errors = String[]
    for section in ("metadata", "config", "metrics")
        if !haskey(data, section)
            push!(errors, "Missing top-level section '$section'")
        end
    end
    !isempty(errors) && return (false, errors)

    meta = data["metadata"]
    if !(meta isa AbstractDict)
        push!(errors, "Section 'metadata' must be a table")
        return (false, errors)
    end
    meta_keys_ok = true
    for k in REQUIRED_METADATA_KEYS
        if !haskey(meta, k) || isempty(strip(string(meta[k])))
            push!(errors, "Metadata missing non-empty key '$k'")
            meta_keys_ok = false
        elseif !(meta[k] isa AbstractString)
            push!(errors, "Metadata key '$k' must be a string")
            meta_keys_ok = false
        end
    end
    if haskey(meta, "page") && !(meta["page"] isa AbstractString)
        push!(errors, "Metadata key 'page' must be a string")
    end
    src_commit = get(meta, "src_commit", nothing)
    if src_commit isa AbstractString && !occursin(r"^[0-9a-f]{40}$", src_commit)
        push!(errors, "metadata.src_commit must be a 40-character lowercase hex SHA-1")
    end
    config_sha = get(meta, "config_sha256", nothing)
    if config_sha isa AbstractString && !occursin(r"^[0-9a-f]{64}$", config_sha)
        push!(errors, "metadata.config_sha256 must be a 64-character lowercase hex SHA-256")
    end

    cfg = data["config"]
    if !(cfg isa AbstractDict)
        push!(errors, "Section 'config' must be a table")
        return (false, errors)
    end
    rec_id = get(meta, "id", nothing) isa AbstractString ? String(meta["id"]) : ""
    if haskey(REQUIRED_CONFIG_KEYS, rec_id)
        for k in REQUIRED_CONFIG_KEYS[rec_id]
            if !haskey(cfg, k)
                push!(errors, "Config section missing required key '$k'")
            end
        end
    else
        push!(
            errors,
            "Unknown record id '$rec_id' has no required config keys in REQUIRED_CONFIG_KEYS",
        )
    end

    metrics = data["metrics"]
    if !(metrics isa Vector) || isempty(metrics)
        push!(errors, "Record must contain at least one [[metrics]] entry")
        return (false, errors)
    end

    seen_ids = Set{String}()
    all_metric_keys_ok = true
    for (idx, m) in enumerate(metrics)
        if !(m isa AbstractDict)
            push!(errors, "Metric #$idx must be a table")
            all_metric_keys_ok = false
            continue
        end
        missing_key = false
        for k in REQUIRED_METRIC_KEYS
            if !haskey(m, k)
                push!(errors, "Metric #$idx missing key '$k'")
                missing_key = true
            end
        end
        if missing_key
            all_metric_keys_ok = false
            continue
        end
        bad_type = false
        for k in ("id", "case", "description", "comparison", "units", "status")
            if !(m[k] isa AbstractString)
                push!(errors, "Metric #$idx key '$k' must be a string")
                bad_type = true
            end
        end
        if !(m["nx"] isa Integer) || m["nx"] isa Bool
            push!(errors, "Metric #$idx key 'nx' must be an integer")
            bad_type = true
        end
        if bad_type
            all_metric_keys_ok = false
            continue
        end

        mid = String(m["id"])
        if mid in seen_ids
            push!(errors, "Duplicate metric id '$mid'")
        end
        push!(seen_ids, mid)

        if !(
            m["value"] isa Real &&
            m["reference"] isa Real &&
            m["atol"] isa Real &&
            m["rtol"] isa Real
        )
            push!(
                errors, "Metric '$mid' must have numeric value, reference, atol, and rtol"
            )
            continue
        end
        val = Float64(m["value"])
        ref = Float64(m["reference"])
        atol = Float64(m["atol"])
        rtol = Float64(m["rtol"])
        comp = string(m["comparison"])
        status = string(m["status"])

        if status != "PASS"
            push!(errors, "Metric '$mid' has status '$status' (expected 'PASS')")
        end
        passed = if comp == "isapprox"
            isapprox(val, ref; atol=atol, rtol=rtol)
        elseif comp == "le_atol"
            if ref != 0.0
                push!(
                    errors,
                    "Metric '$mid' with comparison 'le_atol' must have reference = 0.0 (got $ref)",
                )
            end
            isfinite(val) && 0.0 <= val <= atol
        else
            push!(errors, "Metric '$mid' has unsupported comparison '$comp'")
            false
        end
        if !passed
            push!(
                errors,
                "Metric '$mid' value $val fails comparison '$comp' against reference $ref (atol=$atol, rtol=$rtol)",
            )
        end
    end

    if check_doc_sync && isempty(errors) && haskey(meta, "page")
        doc_path = joinpath(docs_dir, String(meta["page"]))
        if isfile(doc_path)
            doc_text = read(doc_path, String)
            expected_block = render_validation_markdown_block(data)
            if !occursin(expected_block, doc_text)
                push!(
                    errors,
                    "Documentation page '$(meta["page"])' is out of sync with validation record '$(meta["id"])'",
                )
            end
        elseif require_doc_pages
            push!(errors, "Target documentation page not found: $doc_path")
        end
    end

    return (isempty(errors), errors)
end

"""
    check_validation_records(records_dir::String=DEFAULT_RECORDS_DIR; docs_dir::String=VALIDATION_DOCS_DIR, require_doc_pages::Bool=false) -> Bool

Validate all TOML validation records in `records_dir`.

Parameters
----------
records_dir : String, optional
    Directory containing validation record `.toml` files.
docs_dir : String, optional
    Documentation source directory for optional page synchronization checks.
require_doc_pages : Bool, optional
    Whether missing documentation pages are treated as errors.

Returns
-------
Bool
    `true` when at least one record exists and all records pass validation.
"""
function check_validation_records(
    records_dir::String=DEFAULT_RECORDS_DIR;
    docs_dir::String=VALIDATION_DOCS_DIR,
    require_doc_pages::Bool=false,
)::Bool
    if !isdir(records_dir)
        println(stderr, "Validation records directory not found: $records_dir")
        return false
    end
    record_files = sort(filter(f -> endswith(f, ".toml"), readdir(records_dir)))
    if isempty(record_files)
        println(stderr, "No validation record TOML files found in $records_dir")
        return false
    end

    all_ok = true
    for fname in record_files
        fpath = joinpath(records_dir, fname)
        data = TOML.parsefile(fpath)
        ok, errs = validate_record_dict(
            data; docs_dir=docs_dir, require_doc_pages=require_doc_pages
        )
        if !ok
            all_ok = false
            for e in errs
                println(stderr, "[$fname] $e")
            end
        end
    end
    return all_ok
end

"""
    _canonical_field_repr(val)

Return a fixed text form of a config value for the config hash: Bool, integer, float,
Symbol, String and `nothing` in fixed forms, tuples and vectors element by element,
dictionaries with their pairs sorted by key, and any other value by `string(val)`.

Parameters
----------
val : Any
    Value of a config field or of a hashed module constant.

Returns
-------
String
    Text form of `val`.

Raises
------
InexactError
    If `val` is an integer outside the `Int64` range.
"""
function _canonical_field_repr(val)::String
    if val isa Bool
        return val ? "true" : "false"
    elseif val isa Integer
        return string(Int64(val))
    elseif val isa AbstractFloat
        return repr(Float64(val))
    elseif val isa Symbol
        return ":" * String(val)
    elseif val isa AbstractString
        return repr(String(val))
    elseif val === nothing
        return "nothing"
    elseif val isa Tuple || val isa AbstractVector
        return "[" * join((_canonical_field_repr(x) for x in val), ", ") * "]"
    elseif val isa AbstractDict
        isempty(val) && return "{}"
        sorted_pairs = sort!(collect(pairs(val)); by=p -> string(first(p)))
        inner = join(
            [
                "$(_canonical_field_repr(k)) => $(_canonical_field_repr(v))" for
                (k, v) in sorted_pairs
            ],
            ", ",
        )
        return "{$inner}"
    end
    return string(val)
end

"""
    _same_julia_major_minor(recorded_ver, current_ver=VERSION)

Return whether the Julia version of a record and the running version share major and minor.

Parameters
----------
recorded_ver : AbstractString
    Julia version stored in the record, for example `"1.12.1"`.
current_ver : VersionNumber, optional
    Running Julia version.

Returns
-------
Bool
    `true` when major and minor agree; `false` when they differ or `recorded_ver` does not
    start with a version.
"""
function _same_julia_major_minor(
    recorded_ver::AbstractString, current_ver::VersionNumber=VERSION
)::Bool
    m = match(r"^v?(\d+)\.(\d+)(?:\.|\b|$)", strip(recorded_ver))
    m === nothing && return false
    return parse(Int, m.captures[1]) == current_ver.major &&
           parse(Int, m.captures[2]) == current_ver.minor
end

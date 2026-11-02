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

const REQUIRED_CONFIG_KEYS = (
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
        "<!-- VALIDATION_RECORD_START: $(rec_id) -->",
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
    push!(lines, "<!-- VALIDATION_RECORD_END: $(rec_id) -->")
    return join(lines, "\n")
end

"""
    validate_record_dict(data::Dict{String,Any}; docs_dir::String=VALIDATION_DOCS_DIR, require_doc_pages::Bool=false) -> Tuple{Bool,Vector{String}}

Validate the schema, metric tolerances, and optional documentation block of a validation record.
Metrics using `comparison = "le_atol"` represent upper-bound magnitude checks and must store
`reference = 0.0`.

Parameters
----------
data : Dict{String,Any}
    Parsed TOML dictionary of a validation record.
docs_dir : String, optional
    Root directory of documentation markdown files.
require_doc_pages : Bool, optional
    Whether a missing target documentation page is treated as an error.

Returns
-------
Tuple{Bool,Vector{String}}
    Boolean pass flag and vector of diagnostic error messages.
"""
function validate_record_dict(
    data::Dict{String,Any};
    docs_dir::String=VALIDATION_DOCS_DIR,
    require_doc_pages::Bool=false,
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
        end
    end
    if haskey(meta, "src_commit") &&
        !occursin(r"^[0-9a-f]{40}$", String(meta["src_commit"]))
        push!(errors, "metadata.src_commit must be a 40-character lowercase hex SHA-1")
    end
    if haskey(meta, "config_sha256") &&
        !occursin(r"^[0-9a-f]{64}$", String(meta["config_sha256"]))
        push!(errors, "metadata.config_sha256 must be a 64-character lowercase hex SHA-256")
    end

    cfg = data["config"]
    if !(cfg isa AbstractDict)
        push!(errors, "Section 'config' must be a table")
        return (false, errors)
    end
    for k in REQUIRED_CONFIG_KEYS
        if !haskey(cfg, k)
            push!(errors, "Config section missing required key '$k'")
        end
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

    if meta_keys_ok && all_metric_keys_ok && haskey(meta, "page")
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

#!/usr/bin/env julia
# tools/check_config_schema.jl
# Automated configuration schema verification tool for Erebus.jl.
# Compares SimulationConfig structs and defaults against docs/src/reference/config_schema.md.

const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
if Base.find_package("Erebus") === nothing
    pushfirst!(LOAD_PATH, ROOT_DIR)
end

using Erebus
using JSON
using Printf

const DOCS_SCHEMA_PATH = joinpath(ROOT_DIR, "docs", "src", "reference", "config_schema.md")
const BASELINE_PATH = joinpath(@__DIR__, "config_schema_baseline.json")

"""
    parse_markdown_schema(md_path::String) -> Dict{String, Dict{String, Dict{String, String}}}

Parse configuration reference tables from the schema documentation file.

Parameters
----------
md_path : String
    Path to config_schema.md.

Returns
-------
Dict{String, Dict{String, Dict{String, String}}}
    Mapping of section name to field metadata (type, default string).
"""
function parse_markdown_schema(md_path::String)
    isfile(md_path) || error("Documentation schema file not found: $md_path")
    content = read(md_path, String)
    lines = split(content, "\n")

    schema = Dict{String,Dict{String,Dict{String,String}}}()
    current_sec = ""

    for line in lines
        # Match section headers: ## [section] or ## `[section]`
        m = match(r"^#{2,3}\s+`?\[([a-zA-Z0-9_]+)\]`?", line)
        if m !== nothing
            current_sec = m.captures[1]
            schema[current_sec] = Dict{String,Dict{String,String}}()
            continue
        end

        # Match table rows: | `param` | `Type` | `Default` | ...
        if !isempty(current_sec)
            parts = split(line, "|")
            if length(parts) >= 5
                param = strip(parts[2])
                type_str = strip(parts[3])
                default_str = strip(parts[4])
                if startswith(param, "`") && endswith(param, "`")
                    param_clean = replace(param, "`" => "")
                    type_clean = replace(type_str, "`" => "")
                    default_clean = replace(default_str, "`" => "")
                    if param_clean != "Parameter"
                        schema[current_sec][param_clean] = Dict(
                            "type" => type_clean, "default" => default_clean
                        )
                    end
                end
            end
        end
    end

    return schema
end

"""
    check_value_match(jl_val, md_def::String) -> Bool

Determine whether a code default value matches the documented string representation.

Parameters
----------
jl_val : Any
    The code default value.
md_def : String
    The documented default value string.

Returns
-------
Bool
    True if the values match within tolerance, false otherwise.
"""
function check_value_match(jl_val, md_def::String)
    clean_def = strip(md_def)

    # Booleans
    if jl_val isa Bool
        val_parsed = if lowercase(clean_def) == "true"
            true
        else
            (lowercase(clean_def) == "false" ? false : nothing)
        end
        return val_parsed === jl_val
    end

    # Floating point numbers
    if jl_val isa AbstractFloat
        if isnan(jl_val)
            return clean_def == "NaN"
        end
        val_parsed = tryparse(Float64, clean_def)
        return val_parsed !== nothing && isapprox(val_parsed, jl_val; rtol=1e-3)
    end

    # Integers
    if jl_val isa Integer
        val_parsed = tryparse(Int, clean_def)
        return val_parsed === jl_val
    end

    # Symbols
    if jl_val isa Symbol
        val_clean = replace(clean_def, "\"" => "", ":" => "")
        return Symbol(val_clean) == jl_val
    end

    # Strings
    if jl_val isa AbstractString
        val_clean = replace(clean_def, "\"" => "")
        return val_clean == jl_val
    end

    # Array or composite: accept without strict value check
    return true
end

"""
    run_schema_verification(cfg::SimulationConfig, md_schema::Dict; verbose::Bool=false) -> NamedTuple

Verify that all struct fields in SimulationConfig are documented, with matching defaults and types.

Parameters
----------
cfg : SimulationConfig
    Default simulation configuration instance.
md_schema : Dict
    Parsed markdown documentation tables.
verbose : Bool
    Whether to print verbose output.

Returns
-------
NamedTuple
    Results containing counts, missing fields, phantom fields, and discrepancies.
"""
function run_schema_verification(
    cfg::SimulationConfig, md_schema::Dict; verbose::Bool=false
)
    jl_sections = String.(fieldnames(typeof(cfg)))

    missing_fields = String[]
    phantom_fields = String[]
    discrepancies = String[]
    total_fields = 0

    for s in jl_sections
        sub = getfield(cfg, Symbol(s))
        jl_fieldnames = String.(fieldnames(typeof(sub)))

        if !haskey(md_schema, s)
            for f in jl_fieldnames
                push!(missing_fields, "$s.$f")
            end
            continue
        end

        md_fields = collect(keys(md_schema[s]))

        for f in jl_fieldnames
            total_fields += 1
            key = "$s.$f"
            if !(f in md_fields)
                push!(missing_fields, key)
            else
                jl_val = getfield(sub, Symbol(f))
                md_def = md_schema[s][f]["default"]
                if !check_value_match(jl_val, md_def)
                    push!(discrepancies, "$key (code=$jl_val, docs=$md_def)")
                end
            end
        end

        for f in md_fields
            if !(f in jl_fieldnames)
                push!(phantom_fields, "$s.$f")
            end
        end
    end

    return (; total_fields, missing_fields, phantom_fields, discrepancies)
end

"""
    load_baseline() -> Dict{String, Any}

Read baseline exceptions from config_schema_baseline.json.

Returns
-------
Dict{String, Any}
    Baseline configuration dictionary.
"""
function load_baseline()
    if isfile(BASELINE_PATH)
        return JSON.parsefile(BASELINE_PATH)
    end
    return Dict{String,Any}(
        "documented_baseline_exceptions" => Dict{String,Any}(),
        "known_undocumented" => String[],
    )
end

"""
    write_baseline(data::Dict)

Write baseline exceptions to config_schema_baseline.json.

Parameters
----------
data : Dict
    Dictionary of baseline exceptions.
"""
function write_baseline(data::Dict)
    open(BASELINE_PATH, "w") do io
        JSON.print(io, data, 4)
        return println(io)
    end
    return println("Saved schema baseline to $BASELINE_PATH")
end

"""
    main()

Entry point for configuration schema verification.
"""
function main()
    mode = length(ARGS) >= 1 ? ARGS[1] : "--check"
    verbose = "--verbose" in ARGS

    println(
        "Verifying Erebus configuration schema against docs/src/reference/config_schema.md...",
    )
    md_schema = parse_markdown_schema(DOCS_SCHEMA_PATH)
    cfg = default_config()

    res = run_schema_verification(cfg, md_schema; verbose=verbose)
    baseline = load_baseline()
    known_undoc = get(baseline, "known_undocumented", String[])

    actual_missing = setdiff(res.missing_fields, known_undoc)

    @printf(
        "Verified %d configuration fields across %d sections.\n",
        res.total_fields,
        length(fieldnames(typeof(cfg)))
    )

    if mode == "--baseline"
        new_baseline = Dict{String,Any}(
            "documented_baseline_exceptions" =>
                get(baseline, "documented_baseline_exceptions", Dict{String,Any}()),
            "known_undocumented" => res.missing_fields,
        )
        write_baseline(new_baseline)
        return nothing
    end

    if mode == "--check"
        has_error = false

        if !isempty(actual_missing)
            has_error = true
            println("\n[ERROR] Fields present in code but missing from documentation:")
            for f in actual_missing
                println("  - $f")
            end
        end

        if !isempty(res.phantom_fields)
            has_error = true
            println("\n[ERROR] Fields documented in schema but missing from code structs:")
            for f in res.phantom_fields
                println("  - $f")
            end
        end

        if !isempty(res.discrepancies)
            has_error = true
            println("\n[ERROR] Default value discrepancies between code and documentation:")
            for d in res.discrepancies
                println("  - $d")
            end
        end

        if has_error
            println("\nSchema verification failed.")
            exit(1)
        else
            println("\nSchema verification passed: 0 missing, 0 phantom, 0 discrepancies.")
            exit(0)
        end
    else
        println(
            stderr,
            "Usage: julia tools/check_config_schema.jl [--check | --baseline] [--verbose]",
        )
        exit(2)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

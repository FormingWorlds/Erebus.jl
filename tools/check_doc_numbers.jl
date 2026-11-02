#!/usr/bin/env julia
# Validate documented numbers against matching test assertions.

using TOML

include(joinpath(@__DIR__, "validation_records.jl"))

const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_MAP_PATH = joinpath(
    ROOT_DIR, "docs", "src", "reference", "benchmarks", "number_map.toml"
)
const TEST_DIR = joinpath(ROOT_DIR, "test")
const DOCS_DIR = joinpath(ROOT_DIR, "docs", "src")

"""
    check_doc_numbers(map_path::String=DEFAULT_MAP_PATH; docs_dir::String=DOCS_DIR, test_dir::String=TEST_DIR) -> Tuple{Bool,Vector{String}}

Verify each documented number entry against test files and documentation pages.

Parameters
----------
map_path : String, optional
    Path to the TOML mapping file.
docs_dir : String, optional
    Root directory of documentation markdown files.
test_dir : String, optional
    Root directory of test files.

Returns
-------
Tuple{Bool,Vector{String}}
    Status flag and list of validation error messages.
"""
function check_doc_numbers(
    map_path::String=DEFAULT_MAP_PATH; docs_dir::String=DOCS_DIR, test_dir::String=TEST_DIR
)
    if !isfile(map_path)
        return false, ["Number map file not found: $map_path"]
    end

    data = TOML.parsefile(map_path)
    entries = get(data, "entry", Any[])
    if isempty(entries)
        return true, String[]
    end

    errors = String[]
    for (i, entry) in enumerate(entries)
        page = get(entry, "page", nothing)
        number = get(entry, "number", nothing)
        test_name = get(entry, "test_name", nothing)

        if page === nothing || number === nothing || test_name === nothing
            push!(
                errors,
                "Entry $i is missing required keys (page, number, test_name): $entry",
            )
            continue
        end

        num_str = string(number)

        page_path = joinpath(docs_dir, page)
        if !isfile(page_path)
            push!(errors, "Documentation page not found: $page_path (entry $i)")
        else
            page_content = read(page_path, String)
            if !occursin(num_str, page_content)
                push!(
                    errors,
                    "Documentation page $page does not contain number '$num_str' (entry $i)",
                )
            end
        end

        test_path = joinpath(test_dir, test_name)
        if !isfile(test_path)
            push!(errors, "Test file not found: $test_path (entry $i)")
        else
            test_content = read(test_path, String)
            if !occursin(num_str, test_content)
                push!(
                    errors,
                    "Test file $test_name does not contain number '$num_str' (entry $i)",
                )
            end
        end
    end

    return isempty(errors), errors
end

"""
    run_doc_numbers_cli(args::Vector{String}=ARGS; records_dir::String=DEFAULT_RECORDS_DIR) -> Int

Run the CLI check for documented numbers and validation records.

Parameters
----------
args : Vector{String}, optional
    Command-line arguments (optional custom path to `number_map.toml`).
records_dir : String, optional
    Directory containing validation record `.toml` files.

Returns
-------
Int
    Exit code (`0` on success, `1` on validation failure).
"""
function run_doc_numbers_cli(
    args::Vector{String}=ARGS; records_dir::String=DEFAULT_RECORDS_DIR
)::Int
    map_path = length(args) >= 1 ? args[1] : DEFAULT_MAP_PATH
    valid, errors = check_doc_numbers(map_path)
    if !valid
        println(stderr, "check_doc_numbers failed with $(length(errors)) error(s):")
        for err in errors
            println(stderr, "  - ", err)
        end
        return 1
    end
    if !check_validation_records(records_dir; require_doc_pages=true)
        println(stderr, "check_doc_numbers: validation records check failed.")
        return 1
    end
    println("check_doc_numbers: all entries and validation records verified successfully.")
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(run_doc_numbers_cli(ARGS))
end

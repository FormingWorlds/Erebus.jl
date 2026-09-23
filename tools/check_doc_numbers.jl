#!/usr/bin/env julia
# Validate documented numbers against matching test assertions.

using TOML

const ROOT_DIR = normpath(joinpath(@__DIR__, ".."))
const DEFAULT_MAP_PATH = joinpath(ROOT_DIR, "docs", "src", "validation", "number_map.toml")
const TEST_DIR = joinpath(ROOT_DIR, "test")
const DOCS_DIR = joinpath(ROOT_DIR, "docs", "src")

"""
    check_doc_numbers(map_path::String=DEFAULT_MAP_PATH)

Verify each documented number entry against test files and documentation pages.

# Parameters
- `map_path::String`: Path to the TOML mapping file.

# Returns
- `Tuple{Bool, Vector{String}}`: Status flag and list of validation error messages.
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

function main()
    map_path = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_MAP_PATH
    valid, errors = check_doc_numbers(map_path)
    if !valid
        println(stderr, "check_doc_numbers failed with $(length(errors)) error(s):")
        for err in errors
            println(stderr, "  - ", err)
        end
        exit(1)
    end
    println("check_doc_numbers: all entries verified successfully.")
    return exit(0)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

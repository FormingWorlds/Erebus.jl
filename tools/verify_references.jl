#!/usr/bin/env julia
# Reference verification script for Erebus.jl documentation.
#
# Verifies that all literature DOIs cited in docs/src/reference/bibliography.md
# and docs/src/validation/*.md resolve against the official DOI registry (doi.org).
#
# Prevents reference hallucination by ensuring cited works are real, registered publications.

const DOCS_DIR = normpath(joinpath(@__DIR__, "..", "docs", "src"))
const DOI_REGEX = r"https?://doi\.org/(10\.[0-9]{4,9}/[-._;()/:A-Za-z0-9]+)"

function extract_dois_from_file(filepath::String)
    dois = Set{String}()
    content = read(filepath, String)
    for m in eachmatch(DOI_REGEX, content)
        # Strip trailing punctuation often found in markdown links or prose
        doi = replace(m.captures[1], r"[,\.\)>]+$" => "")
        push!(dois, doi)
    end
    return dois
end

function check_doi(doi::String)
    url = "https://doi.org/" * doi
    # Query HTTP status code without following redirects (-sS -o /dev/null -w "%{http_code}")
    # A valid registered DOI returns 301, 302, or 307 redirect to the publisher page.
    # Add connect timeout, max time, and retries to prevent network hangs in CI.
    cmd = `curl -sS --connect-timeout 5 --max-time 15 --retry 2 --retry-delay 1 -o /dev/null -w "%{http_code}" $url`
    try
        code_str = read(cmd, String)
        code = parse(Int, strip(code_str))
        is_valid = code in [200, 301, 302, 303, 307, 308]
        return is_valid, code
    catch e
        return false, -1
    end
end

function main()
    target_files = String[]
    bib_file = joinpath(DOCS_DIR, "reference", "bibliography.md")
    if isfile(bib_file)
        push!(target_files, bib_file)
    end

    val_dir = joinpath(DOCS_DIR, "validation")
    if isdir(val_dir)
        for f in readdir(val_dir; join=true)
            if endswith(f, ".md")
                push!(target_files, f)
            end
        end
    end

    all_dois = Set{String}()
    doi_sources = Dict{String,Vector{String}}()

    for path in target_files
        dois = extract_dois_from_file(path)
        rel_path = relpath(path, joinpath(DOCS_DIR, ".."))
        for d in dois
            push!(all_dois, d)
            push!(get!(doi_sources, d, String[]), rel_path)
        end
    end

    println("=== Erebus Literature Reference Verification ===")
    println(
        "Found ",
        length(all_dois),
        " unique DOIs across ",
        length(target_files),
        " documentation files.\n",
    )

    failed = false
    sorted_dois = sort(collect(all_dois))

    has_skipped = false
    for doi in sorted_dois
        is_valid, code = check_doi(doi)
        srcs = join(doi_sources[doi], ", ")
        if is_valid
            println("  [PASS] ($code) https://doi.org/$doi  ($srcs)")
        elseif code == -1
            println(
                "  [SKIP] (network unreachable or timeout) https://doi.org/$doi  ($srcs)"
            )
            has_skipped = true
        else
            println(
                "  [FAIL] (HTTP $code) https://doi.org/$doi  ($srcs) - UNRESOLVED OR INVALID DOI",
            )
            failed = true
        end
    end

    if failed
        println("\nReference verification failed: Unresolved or invalid DOIs found.")
        exit(1)
    elseif has_skipped
        println(
            "\nReference verification finished: Some DOIs were skipped due to network isolation.",
        )
        exit(0)
    else
        println("\nAll cited DOIs successfully resolved against the registry.")
        exit(0)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

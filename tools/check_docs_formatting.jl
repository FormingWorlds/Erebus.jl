#!/usr/bin/env julia
# Documentation formatting validator for Erebus.jl.
#
# Checks Markdown source and built HTML for formatting defects:
# 1. List items starting with math ($...), which causes Documenter.jl / Markdown
#    stdlib to parse the formula as a block display equation rather than inline text.
# 2. Built HTML files containing math-container blocks inside list items (<li>).

const DOCS_SRC_DIR = normpath(joinpath(@__DIR__, "..", "docs", "src"))
const DOCS_BUILD_DIR = normpath(joinpath(@__DIR__, "..", "docs", "build"))

const LIST_MATH_REGEX = r"^\s*([-*+]|\d+\.)\s*\$"

function check_markdown_files(src_dir::String)
    violations = Tuple{String,Int,String}[]
    for (root, _, files) in walkdir(src_dir)
        for file in files
            endswith(file, ".md") || continue
            filepath = joinpath(root, file)
            rel_path = relpath(filepath, joinpath(src_dir, ".."))
            lines = readlines(filepath)
            for (idx, line) in enumerate(lines)
                if occursin(LIST_MATH_REGEX, line)
                    push!(violations, (rel_path, idx, strip(line)))
                end
            end
        end
    end
    return violations
end

function check_built_html(build_dir::String)
    violations = Tuple{String,Int}[]
    isdir(build_dir) || return violations
    for (root, _, files) in walkdir(build_dir)
        for file in files
            endswith(file, ".html") || continue
            filepath = joinpath(root, file)
            rel_path = relpath(filepath, build_dir)
            content = read(filepath, String)
            matches = collect(eachmatch(r"<li>\s*<p class=\"math-container\">", content))
            if !isempty(matches)
                push!(violations, (rel_path, length(matches)))
            end
        end
    end
    return violations
end

function main()
    println("=== Documentation Formatting Lint ===")
    md_violations = check_markdown_files(DOCS_SRC_DIR)
    html_violations = check_built_html(DOCS_BUILD_DIR)

    has_error = false

    if !isempty(md_violations)
        has_error = true
        println(
            "\n[FAIL] Found $(length(md_violations)) list item(s) starting with math in Markdown source:",
        )
        for (file, line_num, line_text) in md_violations
            println("  $file:$line_num: $line_text")
        end
        println(
            "\nGuidance: Prefix the list item with a descriptive word or variable name so Documenter.jl parses math inline.",
        )
    else
        println("[PASS] Markdown source: No list items start with math.")
    end

    if !isempty(html_violations)
        has_error = true
        println(
            "\n[FAIL] Found $(length(html_violations)) HTML file(s) with block math inside list items:",
        )
        for (file, count) in html_violations
            println("  $file: $count broken list item(s)")
        end
    else
        println("[PASS] Built HTML: No block math containers inside list items.")
    end

    if has_error
        exit(1)
    else
        println("\nAll documentation formatting checks passed.")
        exit(0)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

content = read("tools/check_test_quality.jl", String)

count_by_rule_old = """
function count_by_rule(violations::Vector{Violation})
    counts = Dict{String,Int}("float_equality" => 0, "weak_assert" => 0, "min_asserts" => 0)
    for v in violations
        rule_str = string(v.rule)
        counts[rule_str] = get(counts, rule_str, 0) + 1
    end
    return counts
end
"""

count_by_rule_new = """
function count_by_rule(violations::Vector{Violation})
    counts = Dict{String,Int}()
    for v in violations
        file_basename = basename(v.file)
        key = string(file_basename, ":", v.rule)
        counts[key] = get(counts, key, 0) + 1
    end
    return counts
end
"""
content = replace(content, count_by_rule_old => count_by_rule_new)

main_regression_old = """
        if has_regression
            println("\\nRegressions detected:")
            for v in violations
                rule_str = string(v.rule)
                base_count = get(baseline, rule_str, 0)
                if counts[rule_str] > base_count
                    rel_path = relpath(v.file, normpath(joinpath(TEST_DIR, "..")))
                    println("  ", rel_path, ":", v.line, " [", v.rule, "] ", v.message)
                end
            end
            exit(1)
"""

main_regression_new = """
        if has_regression
            println("\\nRegressions detected:")
            # To avoid printing every violation in a file that regressed, 
            # we just print the file-level regression summary, or all violations in that file.
            for v in violations
                file_basename = basename(v.file)
                key = string(file_basename, ":", v.rule)
                base_count = get(baseline, key, 0)
                if counts[key] > base_count
                    rel_path = relpath(v.file, normpath(joinpath(TEST_DIR, "..")))
                    println("  ", rel_path, ":", v.line, " [", v.rule, "] ", v.message)
                end
            end
            exit(1)
"""
content = replace(content, main_regression_old => main_regression_new)

write("tools/check_test_quality.jl", content)

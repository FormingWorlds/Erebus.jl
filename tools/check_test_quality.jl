#!/usr/bin/env julia
# AST-based test-quality linter for the Erebus.jl test suite.
#
# Enforces the test quality standards defined in .github/.claude/rules/erebus-tests.md:
# 1. No float equality comparisons with `==` (use `≈` or `isapprox`).
# 2. No standalone weak assertions (`!== nothing`, `length > 0`, etc.).
# 3. Leaf testsets must contain at least 2 assertions (no single-assert tests).
#
# Modes:
#   --baseline   Write current violation counts to tools/test_quality_baseline.json
#   --check      CI mode. Compare violation counts to baseline; fail if any count increases.

const TEST_DIR = normpath(joinpath(@__DIR__, "..", "test"))
const BASELINE_PATH = normpath(joinpath(@__DIR__, "test_quality_baseline.json"))

struct Violation
    file::String
    line::Int
    rule::Symbol
    message::String
end

function parse_all_expressions(content::String)
    exprs = Tuple{Any,Int}[]
    pos = 1
    len = length(content)
    while pos <= len
        line_num = count(c -> c == '\n', SubString(content, 1, pos)) + 1
        ex, next_pos = Meta.parse(content, pos)
        ex === nothing && break
        push!(exprs, (ex, line_num))
        pos = next_pos
    end
    return exprs
end

function is_float_literal(x)
    return isa(x, AbstractFloat)
end

function contains_float_equality(node)
    if Meta.isexpr(node, :call) && length(node.args) == 3
        op = node.args[1]
        arg1 = node.args[2]
        arg2 = node.args[3]
        if (op === :(==) || op === :.==) &&
            (is_float_literal(arg1) || is_float_literal(arg2))
            return true
        end
    elseif Meta.isexpr(node, :comparison)
        for i in 2:2:length(node.args)
            if node.args[i] === :(==) || node.args[i] === :.==
                left = node.args[i - 1]
                right = node.args[i + 1]
                if is_float_literal(left) || is_float_literal(right)
                    return true
                end
            end
        end
    end
    if isa(node, Expr)
        for child in node.args
            if contains_float_equality(child)
                return true
            end
        end
    end
    return false
end

function check_float_equality(ex, file::String, line::Int, violations::Vector{Violation})
    # Look for @test with float equality
    if Meta.isexpr(ex, :macrocall) && length(ex.args) >= 3
        macroname = ex.args[1]
        if macroname === Symbol("@test")
            test_arg = ex.args[3]
            if contains_float_equality(test_arg)
                push!(
                    violations,
                    Violation(
                        file,
                        line,
                        :float_equality,
                        "Float equality comparison with `==` or `.==` (use `≈` or `isapprox`)",
                    ),
                )
            end
        end
    end
end

function check_weak_asserts(ex, file::String, line::Int, violations::Vector{Violation})
    if Meta.isexpr(ex, :macrocall) && length(ex.args) >= 3
        macroname = ex.args[1]
        if macroname === Symbol("@test")
            test_arg = ex.args[3]
            # Check for: x !== nothing, x != nothing
            if Meta.isexpr(test_arg, :call) && length(test_arg.args) == 3
                op = test_arg.args[1]
                arg1 = test_arg.args[2]
                arg2 = test_arg.args[3]
                if (op === :(!==) || op === :(!=)) &&
                    (arg1 === :nothing || arg2 === :nothing)
                    push!(
                        violations,
                        Violation(
                            file,
                            line,
                            :weak_assert,
                            "Weak assertion testing `!== nothing` or `!= nothing`",
                        ),
                    )
                end
                # Check for: length(x) > 0 or length(x) >= 1
                if (
                    op === :(>) &&
                    Meta.isexpr(arg1, :call) &&
                    arg1.args[1] === :length &&
                    arg2 == 0
                ) || (
                    op === :(>=) &&
                    Meta.isexpr(arg1, :call) &&
                    arg1.args[1] === :length &&
                    arg2 == 1
                )
                    push!(
                        violations,
                        Violation(
                            file,
                            line,
                            :weak_assert,
                            "Weak assertion testing `length(x) > 0`",
                        ),
                    )
                end
            end
        end
    end
end

function collect_assertions_in_testset(block_ex)
    assert_count = 0
    has_sub_testsets = false

    function walk_inner(node)
        if Meta.isexpr(node, :macrocall) && length(node.args) >= 1
            macroname = node.args[1]
            if macroname === Symbol("@test") ||
                macroname === Symbol("@test_throws") ||
                macroname === Symbol("@test_broken")
                assert_count += 1
            elseif macroname === Symbol("@testset")
                has_sub_testsets = true
                return nothing
            end
        end
        if isa(node, Expr)
            for child in node.args
                walk_inner(child)
            end
        end
    end

    walk_inner(block_ex)
    return assert_count, has_sub_testsets
end

function check_testsets(ex, file::String, line::Int, violations::Vector{Violation})
    if Meta.isexpr(ex, :macrocall) && length(ex.args) >= 3
        macroname = ex.args[1]
        if macroname === Symbol("@testset")
            for arg in ex.args[3:end]
                if Meta.isexpr(arg, :block)
                    assert_count, has_sub_testsets = collect_assertions_in_testset(arg)
                    if !has_sub_testsets && assert_count < 2
                        push!(
                            violations,
                            Violation(
                                file,
                                line,
                                :min_asserts,
                                "Testset contains fewer than 2 assertions (contains $assert_count assertion(s); minimum 2 required)",
                            ),
                        )
                    end
                end
            end
        end
    end
end

function walk_ast(ex, file::String, line::Int, violations::Vector{Violation})
    check_float_equality(ex, file, line, violations)
    check_weak_asserts(ex, file, line, violations)
    check_testsets(ex, file, line, violations)

    if isa(ex, Expr)
        current_line = line
        for child in ex.args
            if isa(child, LineNumberNode)
                current_line = child.line
            elseif isa(child, Expr)
                walk_ast(child, file, current_line, violations)
            end
        end
    end
end

function lint_file(filepath::String)
    violations = Violation[]
    content = read(filepath, String)
    exprs = parse_all_expressions(content)
    for (ex, line) in exprs
        walk_ast(ex, filepath, line, violations)
    end
    return violations
end

function lint_all_tests()
    all_violations = Violation[]
    for (root, _, files) in walkdir(TEST_DIR)
        for file in files
            if endswith(file, ".jl") && file != "runtests.jl"
                filepath = joinpath(root, file)
                v = lint_file(filepath)
                append!(all_violations, v)
            end
        end
    end
    return all_violations
end

function count_by_rule(violations::Vector{Violation})
    counts = Dict{String,Int}("float_equality" => 0, "weak_assert" => 0, "min_asserts" => 0)
    for v in violations
        rule_str = string(v.rule)
        counts[rule_str] = get(counts, rule_str, 0) + 1
    end
    return counts
end

function write_baseline(path::String, counts::Dict{String,Int})
    open(path, "w") do io
        println(io, "{")
        sorted_keys = sort(collect(keys(counts)))
        for (idx, k) in enumerate(sorted_keys)
            comma = idx < length(sorted_keys) ? "," : ""
            println(io, "    \"", k, "\": ", counts[k], comma)
        end
        return println(io, "}")
    end
end

function read_baseline(path::String)
    counts = Dict{String,Int}()
    for line in eachline(path)
        m = match(r"\"([^\"]+)\"\s*:\s*(\d+)", line)
        if m !== nothing
            counts[m.captures[1]] = parse(Int, m.captures[2])
        end
    end
    return counts
end

function main()
    mode = "--check"
    if length(ARGS) >= 1
        mode = ARGS[1]
    end

    violations = lint_all_tests()
    counts = count_by_rule(violations)

    if mode == "--baseline"
        write_baseline(BASELINE_PATH, counts)
        println("Saved baseline to ", BASELINE_PATH, ":")
        for (rule, c) in sort(collect(counts))
            println("  ", rpad(rule, 20), ": ", c)
        end
        exit(0)
    elseif mode == "--check"
        if !isfile(BASELINE_PATH)
            println(stderr, "Error: Baseline file not found at ", BASELINE_PATH)
            println(stderr, "Run `julia tools/check_test_quality.jl --baseline` first.")
            exit(2)
        end

        baseline = read_baseline(BASELINE_PATH)
        has_regression = false

        println("=== Test Quality Lint Report ===")
        for (rule, current_count) in sort(collect(counts))
            base_count = get(baseline, rule, 0)
            status = current_count <= base_count ? "PASS" : "FAIL (REGRESSION)"
            if current_count > base_count
                has_regression = true
            end
            println(
                rpad(rule, 20),
                ": current = ",
                rpad(current_count, 4),
                " baseline = ",
                rpad(base_count, 4),
                " [",
                status,
                "]",
            )
        end

        if has_regression
            println("\nRegressions detected:")
            for v in violations
                rule_str = string(v.rule)
                base_count = get(baseline, rule_str, 0)
                if counts[rule_str] > base_count
                    rel_path = relpath(v.file, normpath(joinpath(TEST_DIR, "..")))
                    println("  ", rel_path, ":", v.line, " [", v.rule, "] ", v.message)
                end
            end
            exit(1)
        else
            println("\nAll test quality rules passed against baseline.")
            exit(0)
        end
    else
        println(stderr, "Usage: julia tools/check_test_quality.jl [--baseline | --check]")
        exit(2)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

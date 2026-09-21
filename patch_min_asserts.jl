content = read("tools/check_test_quality.jl", String)
old_func = """
function collect_assertions_in_testset(block_ex)
    assert_count = 0
    has_sub_testsets = false

    function walk_inner(node)
        if Meta.isexpr(node, :macrocall) && length(node.args) >= 1
            macroname = node.args[1]
            if macroname === Symbol("@test") ||
                macroname === Symbol("@test_throws") ||
                macroname === Symbol("@test_broken") ||
                macroname === Symbol("@reject_config")
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
"""
new_func = """
function collect_assertions_in_testset(block_ex)
    assert_count = 0
    throws_count = 0
    has_sub_testsets = false

    function walk_inner(node)
        if Meta.isexpr(node, :macrocall) && length(node.args) >= 1
            macroname = node.args[1]
            if macroname === Symbol("@test_throws")
                assert_count += 1
                throws_count += 1
            elseif macroname === Symbol("@test") ||
                macroname === Symbol("@test_broken") ||
                macroname === Symbol("@reject_config")
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
    return assert_count, has_sub_testsets, throws_count
end
"""
content = replace(content, old_func => new_func)

old_check = """
                    assert_count, has_sub_testsets = collect_assertions_in_testset(arg)
                    if (assert_count == 1 && !has_sub_testsets) || (!has_sub_testsets && assert_count == 0)
"""
new_check = """
                    assert_count, has_sub_testsets, throws_count = collect_assertions_in_testset(arg)
                    if (assert_count == 1 && throws_count == 0 && !has_sub_testsets) || (!has_sub_testsets && assert_count == 0)
"""
content = replace(content, old_check => new_check)
write("tools/check_test_quality.jl", content)

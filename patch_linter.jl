content = read("tools/check_test_quality.jl", String)

# 1. Update weak assert for typeof and eltype
content = replace(
    content,
    ":typeof\n                    push!(" => ":typeof || arg1.args[1] === :eltype\n                    push!(",
)
content = replace(
    content,
    "Weak assertion testing `typeof(x) == Type`" => "Weak assertion testing `typeof(x) == Type` or `eltype(x) == Type`",
)

# 2. Add bare positivity to weak assert check
bare_positivity_code = """
                # Check for: typeof(x) == Type or typeof(x) === Type
"""
bare_positivity_replacement = """
                # Check for bare positivity: @test x > 0, @test x >= 0, @test x < 0, etc.
                if (op === :(>) || op === :(>=) || op === :(<) || op === :(<=)) &&
                    (arg2 == 0 || arg2 == 0.0 || arg1 == 0 || arg1 == 0.0)
                    push!(
                        violations,
                        Violation(
                            file,
                            line,
                            :weak_assert,
                            "Weak assertion testing bare positivity (e.g. `x > 0`)",
                        ),
                    )
                end
                
                # Check for: typeof(x) == Type or typeof(x) === Type
"""
content = replace(content, bare_positivity_code => bare_positivity_replacement)

# 3. Update min_asserts to not skip parent testsets
min_asserts_code = """
                    assert_count, has_sub_testsets = collect_assertions_in_testset(arg)
                    if !has_sub_testsets && assert_count < 2
"""
min_asserts_replacement = """
                    assert_count, has_sub_testsets = collect_assertions_in_testset(arg)
                    if assert_count == 1 || (!has_sub_testsets && assert_count == 0)
"""
content = replace(content, min_asserts_code => min_asserts_replacement)

# Update return inside collect_assertions_in_testset
collect_assert_code = """
            elseif macroname === Symbol("@testset")
                has_sub_testsets = true
                return nothing
            end
"""
collect_assert_replacement = """
            elseif macroname === Symbol("@testset")
                has_sub_testsets = true
                # Do not return here. Stop recursing into THIS testset, but continue processing siblings.
                # However, walk_inner automatically skips the rest of this node because we return from this block.
                # Actually, wait. walk_inner is called on the arguments.
                return nothing
            end
"""
# Actually, the original just does `return nothing`. That stops processing `node` (the @testset macrocall).
# So it doesn't process its children. But `walk_inner` was called in a loop over `node.args` in the parent.
# Wait, `collect_assertions_in_testset` takes `block_ex`, which is the block INSIDE the parent `@testset`.
# It iterates over `block_ex.args` (lines inside the block).
# If a line is a `@testset`, `walk_inner(node)` returns `nothing`, which is correct: it skips recursing into the child testset.
# So `collect_assertions_in_testset` is already correct. We just need to change the condition `!has_sub_testsets && assert_count < 2`.

# 4. Update float equality heuristic
contains_float_code = """
function contains_float_equality(node)
    if Meta.isexpr(node, :call) && length(node.args) == 3
        op = node.args[1]
        arg1 = node.args[2]
        arg2 = node.args[3]
        if (op === :(==) || op === :.==) &&
            (contains_float_literal(arg1) || contains_float_literal(arg2))
            return true
        end
"""
contains_float_replacement = """
function is_suspected_float_var(x)
    if isa(x, Symbol)
        s = lowercase(string(x))
        return contains(s, "zero") || contains(s, "heat") || contains(s, "temp") || 
               contains(s, "press") || contains(s, "mass") || contains(s, "dens") ||
               contains(s, "vol") || contains(s, "flux") || contains(s, "energy") || 
               contains(s, "frac") || contains(s, "time") || endswith(s, "val")
    elseif Meta.isexpr(x, :ref)
        return true
    end
    return false
end

function contains_float_equality(node)
    if Meta.isexpr(node, :call) && length(node.args) == 3
        op = node.args[1]
        arg1 = node.args[2]
        arg2 = node.args[3]
        if (op === :(==) || op === :.==)
            has_float_literal = contains_float_literal(arg1) || contains_float_literal(arg2)
            has_float_var = is_suspected_float_var(arg1) || is_suspected_float_var(arg2)
            # If there's a literal float, flag it.
            # If both are variables/refs and one matches the float heuristic, flag it.
            if has_float_literal || (has_float_var && !isa(arg1, String) && !isa(arg2, String) && !isa(arg1, Integer) && !isa(arg2, Integer))
                return true
            end
        end
"""
content = replace(content, contains_float_code => contains_float_replacement)

# Also update the comparison loop
comp_code = """
    elseif Meta.isexpr(node, :comparison)
        for i in 2:2:length(node.args)
            if node.args[i] === :(==) || node.args[i] === :.==
                left = node.args[i - 1]
                right = node.args[i + 1]
                if contains_float_literal(left) || contains_float_literal(right)
                    return true
                end
"""
comp_replacement = """
    elseif Meta.isexpr(node, :comparison)
        for i in 2:2:length(node.args)
            if node.args[i] === :(==) || node.args[i] === :.==
                left = node.args[i - 1]
                right = node.args[i + 1]
                has_float_literal = contains_float_literal(left) || contains_float_literal(right)
                has_float_var = is_suspected_float_var(left) || is_suspected_float_var(right)
                if has_float_literal || (has_float_var && !isa(left, String) && !isa(right, String) && !isa(left, Integer) && !isa(right, Integer))
                    return true
                end
"""
content = replace(content, comp_code => comp_replacement)

write("tools/check_test_quality.jl", content)

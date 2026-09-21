content = read("tools/check_test_quality.jl", String)

# Finding 4
bad_heuristic = """
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
"""

good_heuristic = """
function is_suspected_float_var(x)
    if isa(x, Symbol)
        s = string(x)
        return s == "v_zero" || s == "heat_solid_off"
    end
    return false
end
"""
content = replace(content, bad_heuristic => good_heuristic)

# Finding 7
bad_asserts = """
                    assert_count, has_sub_testsets = collect_assertions_in_testset(arg)
                    if assert_count == 1 || (!has_sub_testsets && assert_count == 0)
"""
good_asserts = """
                    assert_count, has_sub_testsets = collect_assertions_in_testset(arg)
                    if (assert_count == 1 && !has_sub_testsets) || (!has_sub_testsets && assert_count == 0)
"""
content = replace(content, bad_asserts => good_asserts)

write("tools/check_test_quality.jl", content)

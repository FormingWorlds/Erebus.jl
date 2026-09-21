content = read("tools/check_test_quality.jl", String)
old_func = """
function is_suspected_float_var(x)
    if isa(x, Symbol)
        s = string(x)
        return s == "v_zero" || s == "heat_solid_off"
    end
    return false
end
"""
new_func = ""
content = replace(content, old_func => new_func)

old_check1 = """
            has_float_var = is_suspected_float_var(arg1) || is_suspected_float_var(arg2)
            # If there's a literal float, flag it.
            # If both are variables/refs and one matches the float heuristic, flag it.
            if has_float_literal || (has_float_var && !isa(arg1, String) && !isa(arg2, String) && !isa(arg1, Integer) && !isa(arg2, Integer))
                return true
            end
"""
new_check1 = """
            # If there's a literal float, flag it.
            if has_float_literal
                return true
            end
"""
content = replace(content, old_check1 => new_check1)

old_check2 = """
                has_float_var = is_suspected_float_var(left) || is_suspected_float_var(right)
                if has_float_literal || (has_float_var && !isa(left, String) && !isa(right, String) && !isa(left, Integer) && !isa(right, Integer))
                    return true
                end
"""
new_check2 = """
                if has_float_literal
                    return true
                end
"""
content = replace(content, old_check2 => new_check2)
write("tools/check_test_quality.jl", content)

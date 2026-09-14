using Test
using Erebus

"""
Assert that `validate_config` rejects an invalid configuration with `ArgumentError`.
"""
macro reject_config(args...)
    if length(args) == 1
        arg = args[1]
        if Meta.isexpr(arg, :(=)) || Meta.isexpr(arg, :kw)
            cfg_expr = :(SimulationConfig(; $(args...)))
        else
            cfg_expr = arg
        end
    else
        cfg_expr = :(SimulationConfig(; $(args...)))
    end
    return quote
        @test_throws ArgumentError validate_config($(esc(cfg_expr)))
    end
end

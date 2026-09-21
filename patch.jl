using Meta

function is_safe_equality_operand(x)
    isa(x, Integer) && return true
    isa(x, String) && return true
    isa(x, Symbol) && return false
    isa(x, QuoteNode) && return true
    isa(x, Bool) && return true
    x === :nothing && return true
    x === :missing && return true
    if Meta.isexpr(x, :call)
        fn = x.args[1]
        if fn in (
            :count,
            :length,
            :size,
            :sizeof,
            :firstindex,
            :lastindex,
            :ndims,
            :axes,
            :typeof,
            :eltype,
            :String,
            :Symbol,
        )
            return true
        end
    end
    return false
end

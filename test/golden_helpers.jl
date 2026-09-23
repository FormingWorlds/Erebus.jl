# test/golden_helpers.jl
# Golden-run comparison helper for exact bitwise verification.

module GoldenHelpers

export compare_golden

"""
    compare_golden(a, b; prefix="")

Compare two simulation states bitwise without numerical tolerance.

Recursively compares dictionaries, named tuples, arrays, and scalar fields.
Throws an `ErrorException` naming the field and first differing index on mismatch.

# Parameters
- `a`: First state container (`Dict` or `NamedTuple`).
- `b`: Second state container (`Dict` or `NamedTuple`).
- `prefix::String`: Field path prefix for nested diagnostics.

# Returns
- `Bool`: `true` if all fields match bitwise.

# Raises
- `ErrorException`: If keys, sizes, or array elements mismatch.
"""
function compare_golden(a, b; prefix::String="")
    raw_keys_a = collect(keys(a))
    raw_keys_b = collect(keys(b))
    str_keys_a = sort(string.(raw_keys_a))
    str_keys_b = sort(string.(raw_keys_b))

    if str_keys_a != str_keys_b
        missing_in_b = setdiff(str_keys_a, str_keys_b)
        missing_in_a = setdiff(str_keys_b, str_keys_a)
        error(
            "Key set mismatch at '$(prefix)':\n" *
            "  Missing in b: $(missing_in_b)\n" *
            "  Missing in a: $(missing_in_a)",
        )
    end

    for k in str_keys_a
        field_name = isempty(prefix) ? k : "$(prefix).$(k)"
        val_a = _get_field(a, k)
        val_b = _get_field(b, k)

        _compare_field(val_a, val_b, field_name)
    end

    return true
end

function _get_field(container::AbstractDict, k::AbstractString)
    if haskey(container, k)
        return container[k]
    elseif haskey(container, Symbol(k))
        return container[Symbol(k)]
    end
    return error("Key '$k' not found in Dict")
end

function _get_field(container::AbstractDict, k::Symbol)
    return _get_field(container, string(k))
end

function _get_field(container::NamedTuple, k::AbstractString)
    sym = Symbol(k)
    if haskey(container, sym)
        return getfield(container, sym)
    end
    return error("Key '$k' not found in NamedTuple")
end

function _get_field(container::NamedTuple, k::Symbol)
    return _get_field(container, string(k))
end

function _compare_field(
    val_a::Union{AbstractDict,NamedTuple},
    val_b::Union{AbstractDict,NamedTuple},
    field_name::String,
)
    return compare_golden(val_a, val_b; prefix=field_name)
end

function _compare_field(val_a::AbstractArray, val_b::AbstractArray, field_name::String)
    if eltype(val_a) !== eltype(val_b)
        error(
            "Type mismatch for field '$(field_name)': eltype(a) = $(eltype(val_a)), eltype(b) = $(eltype(val_b))",
        )
    end
    if size(val_a) != size(val_b)
        error(
            "Size mismatch for field '$(field_name)': size(a) = $(size(val_a)), size(b) = $(size(val_b))",
        )
    end

    for idx in eachindex(val_a, val_b)
        elem_a = val_a[idx]
        elem_b = val_b[idx]
        if elem_a isa Union{AbstractDict,NamedTuple} &&
            elem_b isa Union{AbstractDict,NamedTuple}
            compare_golden(elem_a, elem_b; prefix="$(field_name)[$(idx)]")
        elseif !_bitwise_equal(elem_a, elem_b)
            error(
                "Bitwise mismatch in field '$(field_name)' at index $(idx):\n" *
                "  a[$(idx)] = $(elem_a)\n" *
                "  b[$(idx)] = $(elem_b)",
            )
        end
    end
end

function _compare_field(val_a, val_b, field_name::String)
    if isstructtype(typeof(val_a)) && !isbitstype(typeof(val_a)) && typeof(val_a) === typeof(val_b)
        for fname in fieldnames(typeof(val_a))
            _compare_field(getfield(val_a, fname), getfield(val_b, fname), "$(field_name).$(fname)")
        end
        return nothing
    end
    if !_bitwise_equal(val_a, val_b)
        error(
            "Bitwise mismatch in field '$(field_name)':\n" *
            "  a = $(val_a)\n" *
            "  b = $(val_b)",
        )
    end
end

function _bitwise_equal(x::AbstractFloat, y::AbstractFloat)
    typeof(x) === typeof(y) || return false
    if isnan(x) && isnan(y)
        return true
    end
    return isequal(x, y)
end

function _bitwise_equal(x, y)
    typeof(x) === typeof(y) || return false
    return isequal(x, y)
end

end # module GoldenHelpers

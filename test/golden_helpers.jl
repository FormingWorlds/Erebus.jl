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
    keys_a = sort(collect(keys(a)))
    keys_b = sort(collect(keys(b)))

    if keys_a != keys_b
        missing_in_b = setdiff(keys_a, keys_b)
        missing_in_a = setdiff(keys_b, keys_a)
        error(
            "Key set mismatch at '$(prefix)':\n" *
            "  Missing in b: $(missing_in_b)\n" *
            "  Missing in a: $(missing_in_a)",
        )
    end

    for k in keys_a
        field_name = isempty(prefix) ? string(k) : "$(prefix).$(k)"
        val_a = _get_field(a, k)
        val_b = _get_field(b, k)

        _compare_field(val_a, val_b, field_name)
    end

    return true
end

function _get_field(container::AbstractDict, k)
    return container[k]
end

function _get_field(container::NamedTuple, k::Symbol)
    return getfield(container, k)
end

function _get_field(container::NamedTuple, k::AbstractString)
    return getfield(container, Symbol(k))
end

function _compare_field(val_a::Union{AbstractDict,NamedTuple}, val_b::Union{AbstractDict,NamedTuple}, field_name::String)
    compare_golden(val_a, val_b; prefix=field_name)
end

function _compare_field(val_a::AbstractArray, val_b::AbstractArray, field_name::String)
    if size(val_a) != size(val_b)
        error("Size mismatch for field '$(field_name)': size(a) = $(size(val_a)), size(b) = $(size(val_b))")
    end

    for idx in eachindex(val_a, val_b)
        elem_a = val_a[idx]
        elem_b = val_b[idx]
        if !_bitwise_equal(elem_a, elem_b)
            error(
                "Bitwise mismatch in field '$(field_name)' at index $(idx):\n" *
                "  a[$(idx)] = $(elem_a)\n" *
                "  b[$(idx)] = $(elem_b)",
            )
        end
    end
end

function _compare_field(val_a, val_b, field_name::String)
    if !_bitwise_equal(val_a, val_b)
        error(
            "Bitwise mismatch in field '$(field_name)':\n" *
            "  a = $(val_a)\n" *
            "  b = $(val_b)",
        )
    end
end

function _bitwise_equal(x::AbstractFloat, y::AbstractFloat)
    if isnan(x) && isnan(y)
        return true
    end
    return x == y
end

function _bitwise_equal(x, y)
    return x == y
end

end # module GoldenHelpers

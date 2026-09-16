module ErebusMetalExt

using Erebus
using Metal
using KernelAbstractions

function Erebus.to_device(::Metal.MetalBackend, a::AbstractArray{T}) where {T}
    T === Float64 && throw(
        ArgumentError(
            "Apple Metal does not support 64-bit float (Float64) hardware operations. Convert arrays to Float32 before transferring to Metal device.",
        ),
    )
    return Metal.MtlArray(a)
end

function Erebus.to_device(::Type{<:Metal.MtlArray}, a::AbstractArray{T}) where {T}
    T === Float64 && throw(
        ArgumentError(
            "Apple Metal does not support 64-bit float (Float64) hardware operations. Convert arrays to Float32 before transferring to Metal device.",
        ),
    )
    return Metal.MtlArray(a)
end

end

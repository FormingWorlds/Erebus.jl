module ErebusCUDAExt

using Erebus
using CUDA
using KernelAbstractions

function Erebus.to_device(::CUDA.CUDABackend, a::AbstractArray{T}) where {T}
    return CUDA.CuArray(a)
end

function Erebus.to_device(::Type{<:CUDA.CuArray}, a::AbstractArray{T}) where {T}
    return CUDA.CuArray(a)
end

end

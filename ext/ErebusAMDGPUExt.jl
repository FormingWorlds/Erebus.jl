module ErebusAMDGPUExt

using Erebus
using AMDGPU
using KernelAbstractions

function Erebus.to_device(::AMDGPU.ROCBackend, a::AbstractArray{T}) where {T}
    return AMDGPU.ROCArray(a)
end

function Erebus.to_device(::Type{<:AMDGPU.ROCArray}, a::AbstractArray{T}) where {T}
    return AMDGPU.ROCArray(a)
end

end

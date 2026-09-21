import re

with open("src/physics/magma_ocean_degassing.jl", "r") as f:
    text = f.read()

text = text.replace('using .Erebus: Erebus, SPECIES_AMU\n', '')

helper_code = """
function partial_pressures_to_masses(p_dict::Dict{Symbol,Float64}, P_total::Float64, col_coeff::Float64, amu_dict::Dict{Symbol,Float64})
    mu_bar = sum(p_dict[k] * amu_dict[k] for k in keys(p_dict)) / P_total
    return Dict{Symbol,Float64}(k => (p_dict[k] * amu_dict[k] / mu_bar) * col_coeff for k in keys(p_dict))
end
"""
if "function partial_pressures_to_masses" not in text:
    text = text.replace("using DocStringExtensions\n", "using DocStringExtensions\n" + helper_code)

with open("src/physics/magma_ocean_degassing.jl", "w") as f:
    f.write(text)

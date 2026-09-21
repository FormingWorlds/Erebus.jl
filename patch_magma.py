import re

with open("src/physics/magma_ocean_degassing.jl", "r") as f:
    text = f.read()

# Remove local mu_dict definition
text = re.sub(r'mu_dict = Dict\{Symbol,Float64\}\(\s*:H2 => 2\.01588,[\s\S]*?:SO2 => 64\.066,\s*\)\n', '', text)

# We also need to add `partial_pressures_to_masses` at the top of the file, or somewhere.
helper_code = """
function partial_pressures_to_masses(p_dict::Dict{Symbol,Float64}, P_total::Float64, col_coeff::Float64, amu_dict::Dict{Symbol,Float64})
    mu_bar = sum(p_dict[k] * amu_dict[k] for k in keys(p_dict)) / P_total
    return Dict{Symbol,Float64}(k => (p_dict[k] * amu_dict[k] / mu_bar) * col_coeff for k in keys(p_dict))
end
"""
# insert after module docstring or imports
text = text.replace('using .Erebus: Erebus\n', 'using .Erebus: Erebus, SPECIES_AMU\n' + helper_code)

# Replace instances with helper
text = re.sub(r'mu_bar = sum\(p_dict\[k\] \* mu_dict\[k\] for k in keys\(p_dict\)\) / p_surf_pure\n\s*m_atm_dict = Dict\{Symbol,Float64\}\(k => \(p_dict\[k\] \* mu_dict\[k\] / mu_bar\) \* col_coeff for k in keys\(p_dict\)\)', 'm_atm_dict = partial_pressures_to_masses(p_dict, p_surf_pure, col_coeff, SPECIES_AMU)', text)

text = re.sub(r'mu_bar = sum\(p_dict\[k\] \* mu_dict\[k\] for k in keys\(p_dict\)\) / P_trial\n\s*m_atm_dict = Dict\{Symbol,Float64\}\(k => \(p_dict\[k\] \* mu_dict\[k\] / mu_bar\) \* col_coeff for k in keys\(p_dict\)\)', 'm_atm_dict = partial_pressures_to_masses(p_dict, P_trial, col_coeff, SPECIES_AMU)', text)

with open("src/physics/magma_ocean_degassing.jl", "w") as f:
    f.write(text)

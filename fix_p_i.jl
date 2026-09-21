content = read("src/physics/magma_ocean_degassing.jl", String)
content = replace(
    content,
    "final_p_i = Dict{Symbol,Float64}(k => v / col_coeff for (k, v) in final_atm_i)" => "final_p_i = copy(best_res.p_dict)",
)

content = replace(
    content,
    """
        m_atm_S = (
            get(final_atm_i, :H2S, 0.0) * (32.060 / 34.08088) +
            get(final_atm_i, :SO2, 0.0) * (32.060 / 64.066) +
            get(final_atm_i, :S2, 0.0) * 1.0
        )
        final_p_i = copy(best_res.p_dict)
        final_P_surf = sum(values(final_p_i))
""" => """
           m_atm_S = (
               get(final_atm_i, :H2S, 0.0) * (32.060 / 34.08088) +
               get(final_atm_i, :SO2, 0.0) * (32.060 / 64.066) +
               get(final_atm_i, :S2, 0.0) * 1.0
           )
           for (sp, p) in final_p_i
               final_p_i[sp] = p * scale_global
           end
           final_P_surf = sum(values(final_p_i))
   """,
)
write("src/physics/magma_ocean_degassing.jl", content)

content = read("src/simulation/loop.jl", String)
content = replace(
    content,
    "vent_rates[:CO2]  = get(vent_rates, :CO2, 0.0)  + m_C_step / dt" => "vent_rates[:CO2]  = get(vent_rates, :CO2, 0.0)  + (m_C_step * (44.0095 / 12.011)) / dt",
)
content = replace(
    content,
    "vent_rates[:H2S]  = get(vent_rates, :H2S, 0.0)  + m_S_step / dt" => "vent_rates[:H2S]  = get(vent_rates, :H2S, 0.0)  + (m_S_step * (34.08 / 32.06)) / dt",
)
write("src/simulation/loop.jl", content)

content = read("test/test_atmosphere.jl", String)
content = replace(
    content,
    "rates_co2 = get(sim.vent_rates, :CO2, 0.0) * sim.time_props.dt_val" => "rates_co2 = get(sim.vent_rates, :CO2, 0.0) * sim.time_props.dt_val / (44.0095 / 12.011)",
)
content = replace(
    content,
    "rates_h2s = get(sim.vent_rates, :H2S, 0.0) * sim.time_props.dt_val" => "rates_h2s = get(sim.vent_rates, :H2S, 0.0) * sim.time_props.dt_val / (34.08 / 32.06)",
)
write("test/test_atmosphere.jl", content)

content2 = read("test/test_dehydration_darcy_coupling.jl", String)
content2 = replace(
    content2,
    "expected_rates_add[:CO2] += vent_C / dt" => "expected_rates_add[:CO2] += vent_C * (44.0095 / 12.011) / dt",
)
content2 = replace(
    content2,
    "expected_rates_add[:H2S] += vent_S / dt" => "expected_rates_add[:H2S] += vent_S * (34.08 / 32.06) / dt",
)
content2 = replace(
    content2,
    "expected_rates_drain[:CO2] = vent_C / dt" => "expected_rates_drain[:CO2] = vent_C * (44.0095 / 12.011) / dt",
)
content2 = replace(
    content2,
    "expected_rates_drain[:H2S] = vent_S / dt" => "expected_rates_drain[:H2S] = vent_S * (34.08 / 32.06) / dt",
)
content2 = replace(
    content2,
    "@test isapprox(rates_co2, vent_C, rtol=1e-5)" => "@test isapprox(rates_co2, vent_C * (44.0095 / 12.011), rtol=1e-5)",
)
content2 = replace(
    content2,
    "@test isapprox(rates_h2s, vent_S, rtol=1e-5)" => "@test isapprox(rates_h2s, vent_S * (34.08 / 32.06), rtol=1e-5)",
)
write("test/test_dehydration_darcy_coupling.jl", content2)

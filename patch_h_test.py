import re
with open("test/test_core_volatile_partitioning.jl", "r") as f:
    text = f.read()

text = text.replace(
    "M_C_whole_init = m_sil_total * init_C_ppm + m_met * init_fe_C_ppm",
    "M_H_whole_init = m_sil_total * (init_H2O_wtpct * f_H) + m_met * init_fe_H_ppm\n        M_C_whole_init = m_sil_total * init_C_ppm + m_met * init_fe_C_ppm"
)
text = text.replace(
    "M_C_whole_final = m_sil_total * XCm[1] + m_met * Xfe_C_m[1]",
    "M_H_whole_final = m_sil_total * (XH2Om[1] * f_H) + m_met * Xfe_H_m[1]\n        M_C_whole_final = m_sil_total * XCm[1] + m_met * Xfe_C_m[1]"
)
text = text.replace(
    "@test isapprox(M_C_whole_init, M_C_whole_final; rtol=1e-10)",
    "@test isapprox(M_H_whole_init, M_H_whole_final; rtol=1e-10)\n        @test isapprox(M_C_whole_init, M_C_whole_final; rtol=1e-10)"
)

with open("test/test_core_volatile_partitioning.jl", "w") as f:
    f.write(text)

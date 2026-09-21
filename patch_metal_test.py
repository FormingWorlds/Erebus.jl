import re

with open("test/test_core_volatile_partitioning.jl", "r") as f:
    text = f.read()

# Replace m_sil to m_sil_melt
text = text.replace("m_sil = (1.0 - phi_fe_bulk) * F_melt * rho_sil", "m_sil_melt = (1.0 - phi_fe_bulk) * F_melt * rho_sil")
text = text.replace("m_sil * init_C_ppm", "m_sil_melt * init_C_ppm")
text = text.replace("m_sil * init_N_ppm", "m_sil_melt * init_N_ppm")
text = text.replace("m_sil * init_S_ppm", "m_sil_melt * init_S_ppm")

text = text.replace("M_C_final = m_sil * XCm[m] + m_met * Xfe_C_m[m]", "M_C_final = m_sil_melt * XCm[m] + m_met * Xfe_C_m[m]")
text = text.replace("M_N_final = m_sil * XNm[m] + m_met * Xfe_N_m[m]", "M_N_final = m_sil_melt * XNm[m] + m_met * Xfe_N_m[m]")
text = text.replace("M_S_final = m_sil * XSm[m] + m_met * Xfe_S_m[m]", "M_S_final = m_sil_melt * XSm[m] + m_met * Xfe_S_m[m]")

# We need to add whole-marker conservation test!
whole_marker_test = """
        # Whole-marker mass conservation
        m_sil_total = (1.0 - phi_fe_bulk) * rho_sil
        M_C_whole_init = m_sil_total * init_C_ppm + m_met * init_fe_C_ppm
        M_N_whole_init = m_sil_total * init_N_ppm + m_met * init_fe_N_ppm
        M_S_whole_init = m_sil_total * init_S_ppm + m_met * init_fe_S_ppm
        
        M_C_whole_final = m_sil_total * XCm[m] + m_met * Xfe_C_m[m]
        M_N_whole_final = m_sil_total * XNm[m] + m_met * Xfe_N_m[m]
        M_S_whole_final = m_sil_total * XSm[m] + m_met * Xfe_S_m[m]
        
        @test isapprox(M_C_whole_init, M_C_whole_final; rtol=1e-10)
        @test isapprox(M_N_whole_init, M_N_whole_final; rtol=1e-10)
        @test isapprox(M_S_whole_init, M_S_whole_final; rtol=1e-10)
"""
text = text.replace("@test isapprox(M_S_init, M_S_final; rtol=1e-10)", "@test isapprox(M_S_init, M_S_final; rtol=1e-10)\n" + whole_marker_test)

with open("test/test_core_volatile_partitioning.jl", "w") as f:
    f.write(text)

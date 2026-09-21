import re

with open("test/test_core_volatile_partitioning.jl", "r") as f:
    text = f.read()

whole_marker_test = """
        # Whole-marker mass conservation
        m_sil_total = (1.0 - phi_fe_bulk) * rho_sil
        M_C_whole_init = m_sil_total * init_C_ppm + m_met * init_fe_C_ppm
        M_N_whole_init = m_sil_total * init_N_ppm + m_met * init_fe_N_ppm
        M_S_whole_init = m_sil_total * init_S_ppm + m_met * init_fe_S_ppm
        
        M_C_whole_final = m_sil_total * XCm[1] + m_met * Xfe_C_m[1]
        M_N_whole_final = m_sil_total * XNm[1] + m_met * Xfe_N_m[1]
        M_S_whole_final = m_sil_total * XSm[1] + m_met * Xfe_S_m[1]
        
        @test isapprox(M_C_whole_init, M_C_whole_final; rtol=1e-10)
        @test isapprox(M_N_whole_init, M_N_whole_final; rtol=1e-10)
        @test isapprox(M_S_whole_init, M_S_whole_final; rtol=1e-10)
"""
text = text.replace("@test isapprox(M_S_final, M_S_init; rtol=1.0e-12)", "@test isapprox(M_S_final, M_S_init; rtol=1.0e-12)\n" + whole_marker_test)

with open("test/test_core_volatile_partitioning.jl", "w") as f:
    f.write(text)

import re

with open("test/test_core_volatile_partitioning.jl", "r") as f:
    text = f.read()

# I will delete the lines that do M_H_final, M_C_final using m_sil_melt and replace them 
# with my whole-marker logic. 
# Wait, the whole-marker test is already there. So I can just delete the old ones.

text = re.sub(r'\n\s*M_H_final = m_sil_melt \+ m_met \* Xfe_H_m\[1\]\n', '\n', text)
text = re.sub(r'\s*M_H_final = m_sil_melt \* \(XH2Om\[1\] \* f_H\) \+ m_met \* Xfe_H_m\[1\]\n', '\n', text)
text = re.sub(r'\s*M_C_final = m_sil_melt \* XCm\[1\] \+ m_met \* Xfe_C_m\[1\]\n', '\n', text)
text = re.sub(r'\s*M_N_final = m_sil_melt \* XNm\[1\] \+ m_met \* Xfe_N_m\[1\]\n', '\n', text)
text = re.sub(r'\s*M_S_final = m_sil_melt \* XSm\[1\] \+ m_met \* Xfe_S_m\[1\]\n', '\n', text)

text = re.sub(r'\s*@test isapprox\(M_H_final, M_H_init; rtol=1\.0e-12\)\n', '\n', text)
text = re.sub(r'\s*@test isapprox\(M_C_final, M_C_init; rtol=1\.0e-12\)\n', '\n', text)
text = re.sub(r'\s*@test isapprox\(M_N_final, M_N_init; rtol=1\.0e-12\)\n', '\n', text)
text = re.sub(r'\s*@test isapprox\(M_S_final, M_S_init; rtol=1\.0e-12\)\n', '\n', text)

# For partial equilibration test at line 485:
# M_C_partial = m_sil_melt * XCm[1] + m_met * Xfe_C_m[1]
# @test isapprox(M_C_partial, M_C_init; rtol=1.0e-12)
text = re.sub(r'\s*M_C_partial = m_sil_melt \* XCm\[1\] \+ m_met \* Xfe_C_m\[1\]\n', '\n        M_C_partial = m_sil_total * XCm[1] + m_met * Xfe_C_m[1]\n', text)
text = re.sub(r'\s*@test isapprox\(M_C_partial, M_C_init; rtol=1\.0e-12\)\n', '\n        @test isapprox(M_C_partial, M_C_whole_init; rtol=1.0e-12)\n', text)


with open("test/test_core_volatile_partitioning.jl", "w") as f:
    f.write(text)

import re

with open("test/test_magma_degassing.jl", "r") as f:
    text = f.read()

# I want to replace the circular assertions:
# @test isapprox(sol_red.M_melt_N + sol_red.M_atm_i[:N2] + get(sol_red.M_atm_i, :NH3, 0.0) * (14.007/17.03052), M_tot_N, rtol=1e-5)
# with an independent check using p_i.
new_check = """
        # Non-circular mass conservation check via atmospheric column
        col_coeff_test = (4 * pi * Rp^2) / g
        mu_bar_red = sum(sol_red.p_i[k] * Erebus.SPECIES_AMU[k] for k in keys(sol_red.p_i)) / sol_red.P_surf
        M_atm_N_calc_red = (get(sol_red.p_i, :N2, 0.0) * Erebus.SPECIES_AMU[:N2] / mu_bar_red * col_coeff_test) + 
                           (get(sol_red.p_i, :NH3, 0.0) * Erebus.SPECIES_AMU[:NH3] / mu_bar_red * col_coeff_test) * (14.007/17.03052)
        @test isapprox(sol_red.M_melt_N + M_atm_N_calc_red, M_tot_N, rtol=1e-4)

        mu_bar_ox = sum(sol_ox.p_i[k] * Erebus.SPECIES_AMU[k] for k in keys(sol_ox.p_i)) / sol_ox.P_surf
        M_atm_N_calc_ox = (get(sol_ox.p_i, :N2, 0.0) * Erebus.SPECIES_AMU[:N2] / mu_bar_ox * col_coeff_test) + 
                          (get(sol_ox.p_i, :NH3, 0.0) * Erebus.SPECIES_AMU[:NH3] / mu_bar_ox * col_coeff_test) * (14.007/17.03052)
        @test isapprox(sol_ox.M_melt_N + M_atm_N_calc_ox, M_tot_N, rtol=1e-4)
"""
text = text.replace(
    "@test isapprox(sol_red.M_melt_N + sol_red.M_atm_i[:N2] + get(sol_red.M_atm_i, :NH3, 0.0) * (14.007/17.03052), M_tot_N, rtol=1e-5)",
    new_check
)
text = re.sub(r'@test isapprox\(sol_ox.M_melt_N \+ sol_ox.M_atm_i\[:N2\].*?\n', '', text)

with open("test/test_magma_degassing.jl", "w") as f:
    f.write(text)

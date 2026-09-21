import re

with open("test/test_magma_degassing.jl", "r") as f:
    text = f.read()

text = text.replace(
    "# Under reducing conditions, the atmosphere is H2-dominated, lowering mu_bar.\n        # This increases the relative mass fraction of heavy gases (like N2) in the atmosphere,\n        # drawing more Nitrogen out of the melt.\n        @test sol_red.M_melt_N < sol_ox.M_melt_N",
    "# Under reducing conditions, Nitrogen dissolves chemically as nitride (N3-), which vastly\n        # increases its solubility in the melt compared to physical dissolution under oxidizing conditions.\n        @test sol_red.M_melt_N > sol_ox.M_melt_N"
)

# Replace circular conservation assertions
# From: @test isapprox(sol_red.M_melt_N + sol_red.M_atm_i[:N2] + get(sol_red.M_atm_i, :NH3, 0.0) * (14.007/17.03052), M_tot_N, rtol=1e-5)
# I will just write a snippet to compute it independently.

with open("test/test_magma_degassing.jl", "w") as f:
    f.write(text)

import re

with open("src/physics/magma_ocean_degassing.jl", "r") as f:
    text = f.read()

text = text.replace("for outer_iter in 1:8", "for outer_iter in 1:100")
text = text.replace("z_H = 0.5 * (z_H + z_H_new)", "z_H = 0.2 * z_H_new + 0.8 * z_H")
text = text.replace("z_C = 0.5 * (z_C + z_C_new)", "z_C = 0.2 * z_C_new + 0.8 * z_C")
text = text.replace("z_N = 0.5 * (z_N + z_N_new)", "z_N = 0.2 * z_N_new + 0.8 * z_N")
text = text.replace("z_S = 0.5 * (z_S + z_S_new)", "z_S = 0.2 * z_S_new + 0.8 * z_S")

with open("src/physics/magma_ocean_degassing.jl", "w") as f:
    f.write(text)

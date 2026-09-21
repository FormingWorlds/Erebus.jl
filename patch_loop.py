import re

with open("src/simulation/loop.jl", "r") as f:
    text = f.read()

# Replace hardcoded metric
text = text.replace("(4.0 / 3.0) * rplanet_val", "compute_l3d_metric(rplanet_val)")

with open("src/simulation/loop.jl", "w") as f:
    f.write(text)

with open("src/Erebus.jl", "r") as f:
    erebus = f.read()

erebus = erebus.replace("compute_magma_ocean_core_cooling_timescale,", "compute_magma_ocean_core_cooling_timescale,\n    compute_l3d_metric,")

with open("src/Erebus.jl", "w") as f:
    f.write(erebus)

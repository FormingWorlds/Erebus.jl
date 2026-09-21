import re

with open("test/test_simulation.jl", "r") as f:
    text = f.read()

text = text.replace("    yearlength = 365.25 * 24.0 * 3600.0\n", "    yearlength = Erebus.setup_simulation_config().time.yearlength\n")

with open("test/test_simulation.jl", "w") as f:
    f.write(text)

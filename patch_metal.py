import re

with open("src/physics/metal_partitioning.jl", "r") as f:
    text = f.read()

text = text.replace("XCm[m] = max(0.0, C_sil + dC_sil)", "XCm[m] = max(0.0, C_sil + dC_sil * F_melt_val)")
text = text.replace("XNm[m] = max(0.0, C_sil + dC_sil)", "XNm[m] = max(0.0, C_sil + dC_sil * F_melt_val)")
text = text.replace("XSm[m] = max(0.0, C_sil + dC_sil)", "XSm[m] = max(0.0, C_sil + dC_sil * F_melt_val)")
text = text.replace("XH2Om[m] = max(0.0, C_sil_H + dC_sil)", "XH2Om[m] = max(0.0, C_sil_H + dC_sil * F_melt_val)")

with open("src/physics/metal_partitioning.jl", "w") as f:
    f.write(text)

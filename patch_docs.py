import re

def fix_math(filepath):
    with open(filepath, "r") as f:
        text = f.read()

    # jeans_escape.md
    text = text.replace(
        "surface area ratio $L_{\\text{3D}} = A_{\\text{3D}} / P_{\\text{2D}} = \\frac{4\\pi R_{\\text{planet}}^2}{2\\pi R_{\\text{planet}}} = (4/3) R_{\\text{planet}}$",
        "volumetric ratio $L_{\\text{3D}} = V_{\\text{sphere}} / A_{\\text{disk}} = \\frac{(4/3)\\pi R_{\\text{planet}}^3}{\\pi R_{\\text{planet}}^2} = (4/3) R_{\\text{planet}}$"
    )

    # degassing_and_venting.md
    text = text.replace(
        "this 2D surface boundary flux is scaled to 3D by the ratio of spherical surface area to 2D circular boundary perimeter: $L_{\\text{3D}} = A_{\\text{sphere}} / P_{\\text{2D}} = \\frac{4 \\pi R_{\\text{planet}}^2}{2 \\pi R_{\\text{planet}}} = (4/3) R_{\\text{planet}}$",
        "this 2D volume inventory is scaled to 3D by the ratio of spherical volume to 2D disk area: $L_{\\text{3D}} = V_{\\text{sphere}} / A_{\\text{disk}} = \\frac{(4/3) \\pi R_{\\text{planet}}^3}{\\pi R_{\\text{planet}}^2} = (4/3) R_{\\text{planet}}$"
    )

    with open(filepath, "w") as f:
        f.write(text)

fix_math("docs/src/validation/jeans_escape.md")
fix_math("docs/src/explanations/degassing_and_venting.md")

try:
    fix_math("docs/src/validation/coupled_atmosphere.md")
except:
    pass


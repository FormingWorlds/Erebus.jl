#!/usr/bin/env python3
"""
Generate SVG volatile reservoir architecture diagram for Erebus.jl documentation.

Renders the planetary volatile reservoirs (Core, Mantle/Interior, Magma Ocean,
Coupled Atmosphere, Disk Envelope, and Space Escape) and their exchange pathways.
"""

import os

# Neutral scientific palette
PALETTE = {
    "slate_blue": "#1F4E79",
    "crimson": "#C0392B",
    "green": "#27AE60",
    "amber": "#D97706",
    "gold": "#D4AC0D",
    "purple": "#6C3483",
    "teal": "#16A085",
    "charcoal": "#1A1A1A",
    "dark_gray": "#222222",
    "grid_gray": "#E0E0E0",
    "off_white": "#F5F5F5",
    "white": "#FFFFFF",
}

ASSETS_DIR = os.path.join(os.path.dirname(__file__), "..", "docs", "src", "assets")


def generate_svg():
    os.makedirs(ASSETS_DIR, exist_ok=True)
    svg_path = os.path.join(ASSETS_DIR, "volatile_reservoirs.svg")

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1000 620" width="100%" height="100%">
  <defs>
    <style>
      .title {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; font-size: 19px; font-weight: 700; fill: {PALETTE["charcoal"]}; }}
      .subtitle {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; font-size: 13px; font-weight: 400; fill: {PALETTE["dark_gray"]}; }}
      .box-title {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; font-size: 15px; font-weight: 700; }}
      .box-text {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; font-size: 12px; font-weight: 400; }}
      .arrow-text {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; font-size: 11px; font-weight: 600; fill: {PALETTE["charcoal"]}; }}
      .arrow-sub {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; font-size: 10px; font-weight: 400; fill: {PALETTE["dark_gray"]}; }}
    </style>
    <marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 0 1.5 L 8 5 L 0 8.5 z" fill="{PALETTE["charcoal"]}" />
    </marker>
    <marker id="arrow-cobalt" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 0 1.5 L 8 5 L 0 8.5 z" fill="{PALETTE["slate_blue"]}" />
    </marker>
    <marker id="arrow-amber" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 0 1.5 L 8 5 L 0 8.5 z" fill="{PALETTE["amber"]}" />
    </marker>
    <marker id="arrow-magma" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 0 1.5 L 8 5 L 0 8.5 z" fill="{PALETTE["crimson"]}" />
    </marker>
  </defs>

  <!-- Background Canvas -->
  <rect width="1000" height="620" fill="#FFFFFF" rx="8" />
  <rect x="10" y="10" width="980" height="600" fill="none" stroke="{PALETTE["grid_gray"]}" stroke-width="1" rx="8" />

  <!-- Diagram Title -->
  <text x="40" y="45" class="title">Erebus.jl Volatile Reservoirs and Mass Exchange Architecture</text>
  <text x="40" y="68" class="subtitle">Closed-system elemental tracking, high-temperature speciation, redox buffering, and hydrodynamic escape</text>

  <!-- 1. Disk Gas Envelope Reservoir (Top Left) -->
  <g transform="translate(50, 100)">
    <rect width="250" height="135" rx="6" fill="#F4F7FA" stroke="{PALETTE["slate_blue"]}" stroke-width="1.8" />
    <text x="18" y="28" class="box-title" fill="{PALETTE["slate_blue"]}">Protoplanetary Disk Gas</text>
    <text x="18" y="48" class="box-text" fill="{PALETTE["dark_gray"]}">• Nebular H₂-He background</text>
    <text x="18" y="66" class="box-text" fill="{PALETTE["dark_gray"]}">• Bound envelope M_env_bound</text>
    <text x="18" y="84" class="box-text" fill="{PALETTE["dark_gray"]}">• Gravitational capture (R_cap)</text>
    <text x="18" y="102" class="box-text" fill="{PALETTE["dark_gray"]}">• Photoevaporative boil-off (τ_boil)</text>
    <text x="18" y="120" class="box-text" fill="{PALETTE["dark_gray"]}">• Disk dispersal clearance</text>
  </g>

  <!-- 2. Coupled Atmosphere Reservoir (Center Top) -->
  <g transform="translate(370, 100)">
    <rect width="270" height="230" rx="6" fill="#EDF3F9" stroke="{PALETTE["slate_blue"]}" stroke-width="2.5" />
    <text x="18" y="30" class="box-title" fill="{PALETTE["slate_blue"]}">Coupled 1D Atmosphere</text>
    <text x="18" y="52" class="box-text" fill="{PALETTE["charcoal"]}">AtmosphereState (Elemental &amp; Speciated)</text>
    
    <rect x="14" y="62" width="242" height="60" rx="4" fill="#FFFFFF" stroke="{PALETTE["grid_gray"]}" stroke-width="1" />
    <text x="22" y="80" class="box-text" fill="{PALETTE["slate_blue"]}" font-weight="700">ElementInventory (elem):</text>
    <text x="22" y="98" class="box-text" fill="{PALETTE["dark_gray"]}">H, C, N, S, O [kg] (Conserved)</text>
    <text x="22" y="114" class="box-text" fill="{PALETTE["dark_gray"]}">Machine-precision closure &lt; 10⁻¹²</text>

    <rect x="14" y="130" width="242" height="52" rx="4" fill="#FFFFFF" stroke="{PALETTE["grid_gray"]}" stroke-width="1" />
    <text x="22" y="148" class="box-text" fill="{PALETTE["slate_blue"]}" font-weight="700">SpeciesInventory (species):</text>
    <text x="22" y="166" class="box-text" fill="{PALETTE["dark_gray"]}">10 gases: H₂O, H₂, CO₂, CO, CH₄,</text>
    <text x="22" y="178" class="box-text" fill="{PALETTE["dark_gray"]}">N₂, NH₃, H₂S, S₂, SO₂ [kg]</text>

    <text x="18" y="202" class="box-text" fill="{PALETTE["dark_gray"]}">• Equilibrium log₁₀ fO₂ &amp; dO_buffer</text>
    <text x="18" y="220" class="box-text" fill="{PALETTE["dark_gray"]}">• Radiative Guillot P_surf &amp; T_surf_eq</text>
  </g>

  <!-- 3. Space / Escaped Reservoir (Top Right) -->
  <g transform="translate(710, 100)">
    <rect width="240" height="155" rx="6" fill="#FFF6F0" stroke="{PALETTE["amber"]}" stroke-width="1.8" />
    <text x="18" y="28" class="box-title" fill="{PALETTE["amber"]}">Space / Escaped Reservoir</text>
    <text x="18" y="48" class="box-text" fill="{PALETTE["dark_gray"]}">ElementInventory (escaped)</text>
    <text x="18" y="68" class="box-text" fill="{PALETTE["dark_gray"]}">• Transonic blow-off (λ &lt; 1)</text>
    <text x="18" y="86" class="box-text" fill="{PALETTE["dark_gray"]}">• Kinetic Jeans effusion (λ ≥ 2)</text>
    <text x="18" y="104" class="box-text" fill="{PALETTE["dark_gray"]}">• XUV energy-limited loss</text>
    <text x="18" y="122" class="box-text" fill="{PALETTE["dark_gray"]}">• Multispecies crossover diffusion</text>
    <text x="18" y="140" class="box-text" fill="{PALETTE["dark_gray"]}">• Heavy trace drag / fractionation</text>
  </g>

  <!-- 4. Magma Ocean Reservoir (Middle Left) -->
  <g transform="translate(50, 390)">
    <rect width="250" height="165" rx="6" fill="#FDF3F2" stroke="{PALETTE["crimson"]}" stroke-width="1.8" />
    <text x="18" y="28" class="box-title" fill="{PALETTE["crimson"]}">Surface Magma Ocean</text>
    <text x="18" y="48" class="box-text" fill="{PALETTE["dark_gray"]}">Lagrangian Silicate Melt Markers</text>
    <text x="18" y="68" class="box-text" fill="{PALETTE["dark_gray"]}">• Melt fraction F_melt ≥ F_thresh</text>
    <text x="18" y="86" class="box-text" fill="{PALETTE["dark_gray"]}">• HCNS solubility equilibrium laws</text>
    <text x="18" y="104" class="box-text" fill="{PALETTE["dark_gray"]}">• Dixon (1995) H₂O &amp; CO₂</text>
    <text x="18" y="122" class="box-text" fill="{PALETTE["dark_gray"]}">• Libourel (2003) N₂ solubility</text>
    <text x="18" y="140" class="box-text" fill="{PALETTE["dark_gray"]}">• O'Neill (2002) sulfide/sulfate S</text>
    <text x="18" y="156" class="box-text" fill="{PALETTE["dark_gray"]}">• Signed resorption transfer (V4)</text>
  </g>

  <!-- 5. Mantle & Crustal Interior (Middle Center/Right) -->
  <g transform="translate(370, 390)">
    <rect width="270" height="165" rx="6" fill="#F9F5F8" stroke="{PALETTE["purple"]}" stroke-width="1.8" />
    <text x="18" y="28" class="box-title" fill="{PALETTE["purple"]}">Porous Mantle &amp; Crust</text>
    <text x="18" y="48" class="box-text" fill="{PALETTE["dark_gray"]}">Two-Phase Stokes-Darcy Matrix</text>
    <text x="18" y="68" class="box-text" fill="{PALETTE["dark_gray"]}">• Serpentine &amp; clay mineral water</text>
    <text x="18" y="86" class="box-text" fill="{PALETTE["dark_gray"]}">• Porous Darcy fluid transport (φ)</text>
    <text x="18" y="104" class="box-text" fill="{PALETTE["dark_gray"]}">• Hydrothermal convection cells</text>
    <text x="18" y="122" class="box-text" fill="{PALETTE["dark_gray"]}">• Retention drainage ceilings</text>
    <text x="18" y="140" class="box-text" fill="{PALETTE["dark_gray"]}">• Redox buffer (3 FeO + ½ O₂ ↔ Fe₃O₄)</text>
    <text x="18" y="156" class="box-text" fill="{PALETTE["dark_gray"]}">• Marker oxygen exchange ΔO_buffer</text>
  </g>

  <!-- 6. Metallic Core Reservoir (Bottom Right) -->
  <g transform="translate(710, 390)">
    <rect width="240" height="165" rx="6" fill="#F5F4F6" stroke="{PALETTE["dark_gray"]}" stroke-width="1.8" />
    <text x="18" y="28" class="box-title" fill="{PALETTE["dark_gray"]}">Segregated Metallic Core</text>
    <text x="18" y="48" class="box-text" fill="{PALETTE["dark_gray"]}">Siderophile Metal Alloy (Fe-Ni)</text>
    <text x="18" y="68" class="box-text" fill="{PALETTE["dark_gray"]}">• Fe-FeS eutectic melting (1213 K)</text>
    <text x="18" y="86" class="box-text" fill="{PALETTE["dark_gray"]}">• Darcy percolation &amp; Stokes settling</text>
    <text x="18" y="104" class="box-text" fill="{PALETTE["dark_gray"]}">• Siderophile volatile sequestration</text>
    <text x="18" y="122" class="box-text" fill="{PALETTE["dark_gray"]}">• Core partitioning: S, C, H, N</text>
    <text x="18" y="140" class="box-text" fill="{PALETTE["dark_gray"]}">• Troilite &amp; accessory crystallization</text>
    <text x="18" y="156" class="box-text" fill="{PALETTE["dark_gray"]}">• Deep planetary isolation</text>
  </g>

  <!-- Flow Arrows & Labels -->

  <!-- Disk <-> Atmosphere -->
  <path d="M 300 155 L 360 155" fill="none" stroke="{PALETTE["slate_blue"]}" stroke-width="2" marker-end="url(#arrow-cobalt)" />
  <path d="M 360 175 L 300 175" fill="none" stroke="{PALETTE["slate_blue"]}" stroke-width="2" marker-end="url(#arrow-cobalt)" />
  <text x="310" y="145" class="arrow-text" fill="{PALETTE["slate_blue"]}">Capture</text>
  <text x="312" y="195" class="arrow-text" fill="{PALETTE["slate_blue"]}">Boil-Off</text>

  <!-- Atmosphere -> Space Escape -->
  <path d="M 640 165 L 700 165" fill="none" stroke="{PALETTE["amber"]}" stroke-width="2.5" marker-end="url(#arrow-amber)" />
  <text x="646" y="150" class="arrow-text" fill="{PALETTE["amber"]}">Escape Loss</text>
  <text x="646" y="185" class="arrow-sub">debit elem</text>
  <text x="646" y="197" class="arrow-sub">credit escaped</text>

  <!-- Magma Ocean <-> Atmosphere -->
  <path d="M 230 390 L 380 335" fill="none" stroke="{PALETTE["crimson"]}" stroke-width="2" marker-end="url(#arrow-magma)" />
  <path d="M 390 345 L 245 402" fill="none" stroke="{PALETTE["crimson"]}" stroke-width="2" marker-end="url(#arrow-magma)" />
  <text x="235" y="345" class="arrow-text" fill="{PALETTE["crimson"]}">Equilibrium Degassing</text>
  <text x="245" y="360" class="arrow-sub">Signed Resorption (dM &lt; 0)</text>

  <!-- Mantle Venting -> Atmosphere -->
  <path d="M 425 390 L 425 338" fill="none" stroke="{PALETTE["purple"]}" stroke-width="2.5" marker-end="url(#arrow)" />
  <text x="432" y="358" class="arrow-text" fill="{PALETTE["purple"]}">Hydrothermal Venting</text>
  <text x="432" y="373" class="arrow-sub">Darcy Pore &amp; Drainage</text>

  <!-- Mantle <-> Atmosphere Redox Buffer Exchange -->
  <path d="M 575 390 L 575 338" fill="none" stroke="{PALETTE["charcoal"]}" stroke-width="1.8" stroke-dasharray="4,3" marker-end="url(#arrow)" />
  <text x="582" y="358" class="arrow-text">ΔO_buffer</text>
  <text x="582" y="373" class="arrow-sub">3 FeO + ½ O₂ ↔ Fe₃O₄</text>

  <!-- Mantle -> Core Differentiation -->
  <path d="M 640 470 L 700 470" fill="none" stroke="{PALETTE["dark_gray"]}" stroke-width="2.2" marker-end="url(#arrow)" />
  <text x="645" y="455" class="arrow-text">Metal Partitioning</text>
  <text x="645" y="490" class="arrow-sub">Core Percolation</text>
</svg>
"""

    with open(svg_path, "w", encoding="utf-8") as f:
        f.write(svg)
    print(f"Generated volatile reservoir SVG at: {svg_path}")


if __name__ == "__main__":
    generate_svg()

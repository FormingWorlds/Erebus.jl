#!/usr/bin/env python3
"""
Generate SVG architecture and execution flowchart for Erebus.jl documentation.

Renders the multi-physics execution pipeline across initialization, environment
coupling, core staggered solvers, post-solve transport, and persistence.
"""

import os

# Neutral scientific palette
PALETTE = {
    "slate_blue": "#1F4E79",
    "steel_blue": "#2B6CB0",
    "crimson": "#C0392B",
    "amber": "#D97706",
    "emerald": "#27AE60",
    "purple": "#6C3483",
    "charcoal": "#1A1A1A",
    "dark_slate": "#334155",
    "muted_gray": "#64748B",
    "border_gray": "#D0D7DE",
    "card_bg": "#F8FAFC",
    "header_bg": "#F1F5F9",
    "white": "#FFFFFF",
}

ASSETS_DIR = os.path.join(os.path.dirname(__file__), "..", "docs", "src", "assets")


def generate_svg():
    os.makedirs(ASSETS_DIR, exist_ok=True)
    svg_path = os.path.join(ASSETS_DIR, "erebus_architecture_flowchart.svg")

    svg = f"""<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1240 1050" width="100%" height="100%">
  <defs>
    <style>
      .title-text {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; font-weight: 700; fill: {PALETTE["charcoal"]}; font-size: 21px; }}
      .subtitle-text {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; font-weight: 400; fill: {PALETTE["dark_slate"]}; font-size: 13px; }}
      .section-title {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; font-weight: 700; font-size: 12px; letter-spacing: 0.5px; }}
      .tag-text {{ font-family: ui-monospace, Menlo, Consolas, "Liberation Mono", monospace; font-weight: 600; font-size: 9.5px; letter-spacing: 0.5px; }}
      .box-title {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; font-weight: 600; font-size: 13px; }}
      .box-heading {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; font-weight: 600; font-size: 11px; }}
      .box-body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif; font-weight: 400; fill: {PALETTE["dark_slate"]}; font-size: 11px; }}
      .box-code {{ font-family: ui-monospace, Menlo, Consolas, "Liberation Mono", monospace; font-weight: 600; font-size: 10.5px; }}
      .arrow-text {{ font-family: ui-monospace, Menlo, Consolas, "Liberation Mono", monospace; font-weight: 600; font-size: 10px; fill: {PALETTE["slate_blue"]}; }}
    </style>

    <!-- Arrow Markers -->
    <marker id="arrow-slate" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 0 1.5 L 8 5 L 0 8.5 z" fill="{PALETTE["slate_blue"]}" />
    </marker>
    <marker id="arrow-steel" viewBox="0 0 10 10" refX="8" refY="5" markerWidth="6" markerHeight="6" orient="auto-start-reverse">
      <path d="M 0 1.5 L 8 5 L 0 8.5 z" fill="{PALETTE["steel_blue"]}" />
    </marker>
  </defs>

  <!-- Background Canvas -->
  <rect width="1240" height="1050" fill="{PALETTE["white"]}" rx="8" />
  <rect x="10" y="10" width="1220" height="1030" fill="none" stroke="{PALETTE["border_gray"]}" stroke-width="1" rx="8" />

  <!-- Header Section -->
  <g transform="translate(60, 42)">
    <text class="title-text" y="0">Erebus.jl Architecture and Execution Flow</text>
    <text class="subtitle-text" y="22">Two-Dimensional Coupled Thermo-Hydro-Mechanical-Chemical Planetesimal Evolution</text>
    <line x1="0" y1="32" x2="1120" y2="32" stroke="{PALETTE["border_gray"]}" stroke-width="1" />
  </g>

  <!-- SECTION 1: INITIALIZATION -->
  <g transform="translate(60, 95)">
    <rect width="1120" height="136" rx="8" fill="{PALETTE["card_bg"]}" stroke="{PALETTE["border_gray"]}" stroke-width="1.2" />
    <path d="M 0 8 Q 0 0 8 0 L 1112 0 Q 1120 0 1120 8 L 1120 28 L 0 28 Z" fill="{PALETTE["header_bg"]}" />
    <line x1="0" y1="28" x2="1120" y2="28" stroke="{PALETTE["border_gray"]}" stroke-width="1" />
    <text class="section-title" x="16" y="19" fill="{PALETTE["slate_blue"]}">1. SIMULATION INITIALIZATION AND DOMAIN SETUP</text>
    
    <rect x="940" y="6" width="165" height="16" rx="4" fill="{PALETTE["slate_blue"]}" fill-opacity="0.12" />
    <text class="tag-text" x="1022" y="18" text-anchor="middle" fill="{PALETTE["slate_blue"]}">STAGE: PRE-FLIGHT</text>

    <!-- Card 1: Configuration -->
    <g transform="translate(16, 38)">
      <rect width="348" height="86" rx="6" fill="{PALETTE["white"]}" stroke="{PALETTE["border_gray"]}" stroke-width="1" />
      <text class="box-title" x="12" y="20" fill="{PALETTE["slate_blue"]}">TOML Configuration Parsing</text>
      <text class="box-body" x="12" y="38">• Validates <tspan class="box-code" fill="{PALETTE["charcoal"]}">SimulationConfig</tspan> schema via <tspan class="box-code" fill="{PALETTE["charcoal"]}">load_config</tspan></text>
      <text class="box-body" x="12" y="55">• Sets material regimes, thermodynamics, and physical bounds</text>
      <text class="box-body" x="12" y="72">• Direct sparse linear solver setup (<tspan class="box-code" fill="{PALETTE["charcoal"]}">UMFPACK</tspan>)</text>
    </g>

    <!-- Card 2: Mesh & Fields -->
    <g transform="translate(386, 38)">
      <rect width="348" height="86" rx="6" fill="{PALETTE["white"]}" stroke="{PALETTE["border_gray"]}" stroke-width="1" />
      <text class="box-title" x="12" y="20" fill="{PALETTE["slate_blue"]}">Staggered Eulerian Grids</text>
      <text class="box-body" x="12" y="38">• Basic nodes <tspan class="box-code" fill="{PALETTE["charcoal"]}">(Nx, Ny)</tspan>: shear stress <tspan class="box-code" fill="{PALETTE["charcoal"]}">σxy</tspan>, viscosity <tspan class="box-code" fill="{PALETTE["charcoal"]}">η</tspan></text>
      <text class="box-body" x="12" y="55">• Cell centers <tspan class="box-code" fill="{PALETTE["charcoal"]}">(Nx-1, Ny-1)</tspan>: normal stress, <tspan class="box-code" fill="{PALETTE["charcoal"]}">Pt, Pf, T, φ</tspan></text>
      <text class="box-body" x="12" y="72">• Staggered faces: velocities <tspan class="box-code" fill="{PALETTE["charcoal"]}">vx, vy</tspan> and Darcy flux <tspan class="box-code" fill="{PALETTE["charcoal"]}">qx, qy</tspan></text>
    </g>

    <!-- Card 3: Markers -->
    <g transform="translate(756, 38)">
      <rect width="348" height="86" rx="6" fill="{PALETTE["white"]}" stroke="{PALETTE["border_gray"]}" stroke-width="1" />
      <text class="box-title" x="12" y="20" fill="{PALETTE["slate_blue"]}">Marker-in-Cell Initialization</text>
      <text class="box-body" x="12" y="38">• Particles placed in rock, core, and sticky-air zones</text>
      <text class="box-body" x="12" y="55">• Tracks phase <tspan class="box-code" fill="{PALETTE["charcoal"]}">tm</tspan>, melt fraction <tspan class="box-code" fill="{PALETTE["charcoal"]}">Fm</tspan>, metal fraction <tspan class="box-code" fill="{PALETTE["charcoal"]}">Xfe</tspan></text>
      <text class="box-body" x="12" y="72">• Initializes volatile mass fractions: <tspan class="box-code" fill="{PALETTE["charcoal"]}">H, C, N, S, O</tspan></text>
    </g>
  </g>

  <!-- Connector: 1 -> 2 -->
  <g transform="translate(620, 231)">
    <line x1="0" y1="0" x2="0" y2="24" stroke="{PALETTE["slate_blue"]}" stroke-width="1.8" marker-end="url(#arrow-slate)" />
    <text class="arrow-text" x="10" y="15">initialize state</text>
  </g>

  <!-- SECTION 2: ENVIRONMENT & HEATING -->
  <g transform="translate(60, 258)">
    <rect width="1120" height="126" rx="8" fill="{PALETTE["card_bg"]}" stroke="{PALETTE["border_gray"]}" stroke-width="1.2" />
    <path d="M 0 8 Q 0 0 8 0 L 1112 0 Q 1120 0 1120 8 L 1120 28 L 0 28 Z" fill="{PALETTE["header_bg"]}" />
    <line x1="0" y1="28" x2="1120" y2="28" stroke="{PALETTE["border_gray"]}" stroke-width="1" />
    <text class="section-title" x="16" y="19" fill="{PALETTE["steel_blue"]}">2. ENVIRONMENT AND HEATING SOURCE COUPLING</text>
    
    <rect x="915" y="6" width="190" height="16" rx="4" fill="{PALETTE["steel_blue"]}" fill-opacity="0.12" />
    <text class="tag-text" x="1010" y="18" text-anchor="middle" fill="{PALETTE["steel_blue"]}">TIME EVOLUTION STAGE</text>

    <!-- Card 1: Disk -->
    <g transform="translate(16, 36)">
      <rect width="348" height="78" rx="6" fill="{PALETTE["white"]}" stroke="{PALETTE["border_gray"]}" stroke-width="1" />
      <text class="box-title" x="12" y="19" fill="{PALETTE["steel_blue"]}">Protoplanetary Disk and Envelope</text>
      <text class="box-body" x="12" y="36">• Dispersal sigmoid: <tspan class="box-code" fill="{PALETTE["charcoal"]}">w_disp(t) → Tamb(t), Pamb(t)</tspan></text>
      <text class="box-body" x="12" y="52">• Flared disk irradiation and viscous heating models</text>
      <text class="box-body" x="12" y="68">• Gas envelope capture and radiative surface balance</text>
    </g>

    <!-- Card 2: Sources -->
    <g transform="translate(386, 36)">
      <rect width="348" height="78" rx="6" fill="{PALETTE["white"]}" stroke="{PALETTE["border_gray"]}" stroke-width="1" />
      <text class="box-title" x="12" y="19" fill="{PALETTE["steel_blue"]}">Radiogenic and Accretion Sources</text>
      <text class="box-body" x="12" y="36">• Radiogenic decay: <tspan class="box-code" fill="{PALETTE["charcoal"]}">radiogenic_heating!</tspan> (²⁶Al, ⁶⁰Fe)</text>
      <text class="box-body" x="12" y="52">• Mass accretion: <tspan class="box-code" fill="{PALETTE["charcoal"]}">accrete!</tspan> and impact kinetic energy</text>
      <text class="box-body" x="12" y="68">• Potential energy dissipation heat: <tspan class="box-code" fill="{PALETTE["charcoal"]}">Qseg</tspan></text>
    </g>

    <!-- Card 3: Boundary & Gravity -->
    <g transform="translate(756, 36)">
      <rect width="348" height="78" rx="6" fill="{PALETTE["white"]}" stroke="{PALETTE["border_gray"]}" stroke-width="1" />
      <text class="box-title" x="12" y="19" fill="{PALETTE["steel_blue"]}">Boundaries and Gravity Field</text>
      <text class="box-body" x="12" y="36">• Sticky-air domain: rigid lid exterior boundary</text>
      <text class="box-body" x="12" y="52">• Radiative Robin boundary condition at surface</text>
      <text class="box-body" x="12" y="68">• Gravity profile: 2D Poisson solve or 3D enclosed mass</text>
    </g>
  </g>

  <!-- Connector: 2 -> 3 -->
  <g transform="translate(620, 384)">
    <line x1="0" y1="0" x2="0" y2="24" stroke="{PALETTE["slate_blue"]}" stroke-width="1.8" marker-end="url(#arrow-slate)" />
    <text class="arrow-text" x="10" y="15">assemble and solve</text>
  </g>

  <!-- SECTION 3: COUPLED MULTI-PHYSICS SOLVERS -->
  <g transform="translate(60, 411)">
    <rect width="1120" height="324" rx="8" fill="{PALETTE["card_bg"]}" stroke="{PALETTE["border_gray"]}" stroke-width="1.2" />
    <path d="M 0 8 Q 0 0 8 0 L 1112 0 Q 1120 0 1120 8 L 1120 28 L 0 28 Z" fill="{PALETTE["header_bg"]}" />
    <line x1="0" y1="28" x2="1120" y2="28" stroke="{PALETTE["border_gray"]}" stroke-width="1" />
    <text class="section-title" x="16" y="19" fill="{PALETTE["charcoal"]}">3. COUPLED MULTI-PHYSICS CORE SOLVERS</text>
    
    <rect x="915" y="6" width="190" height="16" rx="4" fill="{PALETTE["charcoal"]}" fill-opacity="0.1" />
    <text class="tag-text" x="1010" y="18" text-anchor="middle" fill="{PALETTE["charcoal"]}">MONOLITHIC AND SPLIT SOLVERS</text>

    <!-- Column 1: Thermochemical & Melting -->
    <g transform="translate(16, 38)">
      <rect width="260" height="272" rx="6" fill="{PALETTE["white"]}" stroke="{PALETTE["crimson"]}" stroke-width="1.2" />
      <path d="M 0 6 Q 0 0 6 0 L 254 0 Q 260 0 260 6 L 260 24 L 0 24 Z" fill="{PALETTE["crimson"]}" fill-opacity="0.12" />
      <text class="box-title" x="10" y="17" fill="{PALETTE["crimson"]}">Thermochemical and Melting</text>
      
      <text class="box-heading" x="10" y="42" fill="{PALETTE["crimson"]}">Phase State and Melting:</text>
      <text class="box-body" x="10" y="58">• Solidus/liquidus: <tspan class="box-code" fill="{PALETTE["charcoal"]}">Ts(P), Tl(P)</tspan></text>
      <text class="box-body" x="10" y="74">• Melt fraction: <tspan class="box-code" fill="{PALETTE["charcoal"]}">Fm ∈ [0, 1]</tspan></text>
      <text class="box-body" x="10" y="90">• Latent heat buffer: <tspan class="box-code" fill="{PALETTE["charcoal"]}">ρ·cp,eff</tspan></text>
      
      <text class="box-heading" x="10" y="114" fill="{PALETTE["crimson"]}">Convective Closures:</text>
      <text class="box-body" x="10" y="130">• Soft turbulence: <tspan class="box-code" fill="{PALETTE["charcoal"]}">Nu ~ Ra^(1/3)</tspan></text>
      <text class="box-body" x="10" y="146">• Hydrothermal: <tspan class="box-code" fill="{PALETTE["charcoal"]}">Rayleigh-Darcy Nu_eff</tspan></text>
      <text class="box-body" x="10" y="162">• Regularized conductivity: <tspan class="box-code" fill="{PALETTE["charcoal"]}">keff</tspan></text>
      
      <text class="box-heading" x="10" y="186" fill="{PALETTE["crimson"]}">Rheological Softening:</text>
      <text class="box-body" x="10" y="202">• Exponential matrix softening</text>
      <text class="box-body" x="10" y="218">• Suspension breakdown (<tspan class="box-code" fill="{PALETTE["charcoal"]}">Fm &gt; 0.4</tspan>)</text>
      <text class="box-body" x="10" y="234">• Viscosity bounds: <tspan class="box-code" fill="{PALETTE["charcoal"]}">[ηmin, ηmax]</tspan></text>
      <text class="box-body" x="10" y="254">• Clay dehydration pore fluid</text>
    </g>

    <!-- Column 2: Stokes-Darcy Hydromechanics -->
    <g transform="translate(292, 38)">
      <rect width="260" height="272" rx="6" fill="{PALETTE["white"]}" stroke="{PALETTE["slate_blue"]}" stroke-width="1.2" />
      <path d="M 0 6 Q 0 0 6 0 L 254 0 Q 260 0 260 6 L 260 24 L 0 24 Z" fill="{PALETTE["slate_blue"]}" fill-opacity="0.12" />
      <text class="box-title" x="10" y="17" fill="{PALETTE["slate_blue"]}">Stokes-Darcy Hydromechanics</text>
      
      <text class="box-heading" x="10" y="42" fill="{PALETTE["slate_blue"]}">Coupled Flow Equations:</text>
      <text class="box-body" x="10" y="58">• Solid matrix: <tspan class="box-code" fill="{PALETTE["charcoal"]}">-∇Pt + ∇·σ' + ρt g = 0</tspan></text>
      <text class="box-body" x="10" y="74">• Darcy flux: <tspan class="box-code" fill="{PALETTE["charcoal"]}">q = -(k/ηf)(∇Pf - ρf g)</tspan></text>
      <text class="box-body" x="10" y="90">• Fluid conservation and continuity</text>

      <text class="box-heading" x="10" y="114" fill="{PALETTE["slate_blue"]}">Poroelastic Mechanics:</text>
      <text class="box-body" x="10" y="130">• Matrix and fluid compressibility</text>
      <text class="box-body" x="10" y="146">• Biot-Willis coefficient: <tspan class="box-code" fill="{PALETTE["charcoal"]}">KBW</tspan></text>
      <text class="box-body" x="10" y="162">• Dynamic compaction and dilation</text>

      <text class="box-heading" x="10" y="186" fill="{PALETTE["slate_blue"]}">Terzaghi Stress and Failure:</text>
      <text class="box-body" x="10" y="202">• Effective stress: <tspan class="box-code" fill="{PALETTE["charcoal"]}">Peff = Pt - Pf</tspan></text>
      <text class="box-body" x="10" y="218">• Mohr-Coulomb shear yielding</text>
      <text class="box-body" x="10" y="234">• Hydrofracture permeability ramp</text>
      <text class="box-body" x="10" y="254">• Sparse LSE solve: <tspan class="box-code" fill="{PALETTE["charcoal"]}">Ku = f</tspan></text>
    </g>

    <!-- Column 3: Metal & Magma Segregation -->
    <g transform="translate(568, 38)">
      <rect width="260" height="272" rx="6" fill="{PALETTE["white"]}" stroke="{PALETTE["amber"]}" stroke-width="1.2" />
      <path d="M 0 6 Q 0 0 6 0 L 254 0 Q 260 0 260 6 L 260 24 L 0 24 Z" fill="{PALETTE["amber"]}" fill-opacity="0.12" />
      <text class="box-title" x="10" y="17" fill="{PALETTE["amber"]}">Metal and Magma Segregation</text>
      
      <text class="box-heading" x="10" y="42" fill="{PALETTE["amber"]}">Dual-Regime Metal Drift-Flux:</text>
      <text class="box-body" x="10" y="58">• Darcy percolation: <tspan class="box-code" fill="{PALETTE["charcoal"]}">Fm ≤ 0.40</tspan></text>
      <text class="box-body" x="10" y="74">• Stokes settling: <tspan class="box-code" fill="{PALETTE["charcoal"]}">Fm ≥ 0.50</tspan></text>
      <text class="box-body" x="10" y="90">• Cubic Hermite smooth handover</text>

      <text class="box-heading" x="10" y="114" fill="{PALETTE["amber"]}">Microphysics Corrections:</text>
      <text class="box-body" x="10" y="130">• Richardson-Zaki hindered drag</text>
      <text class="box-body" x="10" y="146">• Weber droplet size equilibrium</text>
      <text class="box-body" x="10" y="162">• Dynamic sulfur density: <tspan class="box-code" fill="{PALETTE["charcoal"]}">ρ(wS)</tspan></text>

      <text class="box-heading" x="10" y="186" fill="{PALETTE["amber"]}">Silicate Magma Ascent:</text>
      <text class="box-body" x="10" y="202">• Buoyant melt percolation</text>
      <text class="box-body" x="10" y="218">• Subsolidus crystallization</text>
      <text class="box-body" x="10" y="234">• Partitioning: <tspan class="box-code" fill="{PALETTE["charcoal"]}">D(H, C, N, S)</tspan></text>
      <text class="box-body" x="10" y="254">• Conservative dissipation: <tspan class="box-code" fill="{PALETTE["charcoal"]}">Qseg</tspan></text>
    </g>

    <!-- Column 4: Thermal Energy Solve -->
    <g transform="translate(844, 38)">
      <rect width="260" height="272" rx="6" fill="{PALETTE["white"]}" stroke="{PALETTE["emerald"]}" stroke-width="1.2" />
      <path d="M 0 6 Q 0 0 6 0 L 254 0 Q 260 0 260 6 L 260 24 L 0 24 Z" fill="{PALETTE["emerald"]}" fill-opacity="0.12" />
      <text class="box-title" x="10" y="17" fill="{PALETTE["emerald"]}">Thermal Energy Solve</text>
      
      <text class="box-heading" x="10" y="42" fill="{PALETTE["emerald"]}">Energy Conservation:</text>
      <text class="box-body" x="10" y="58">• Capacity: <tspan class="box-code" fill="{PALETTE["charcoal"]}">ρ·cp,eff · ∂T/∂t</tspan></text>
      <text class="box-body" x="10" y="74">• Heat conduction: <tspan class="box-code" fill="{PALETTE["charcoal"]}">∇·(keff ∇T)</tspan></text>
      <text class="box-body" x="10" y="90">• Fluid heat advection: <tspan class="box-code" fill="{PALETTE["charcoal"]}">q · ∇T</tspan></text>

      <text class="box-heading" x="10" y="114" fill="{PALETTE["emerald"]}">Coupled Source Terms:</text>
      <text class="box-body" x="10" y="130">• Radiogenic volumetric: <tspan class="box-code" fill="{PALETTE["charcoal"]}">Qrad</tspan></text>
      <text class="box-body" x="10" y="146">• Dissipation: <tspan class="box-code" fill="{PALETTE["charcoal"]}">Qseg + Qvisc</tspan></text>
      <text class="box-body" x="10" y="162">• Latent phase heat: <tspan class="box-code" fill="{PALETTE["charcoal"]}">Qlat</tspan></text>

      <text class="box-heading" x="10" y="186" fill="{PALETTE["emerald"]}">Implicit Discretization:</text>
      <text class="box-body" x="10" y="202">• Radiative Robin boundary</text>
      <text class="box-body" x="10" y="218">• Harmonic interface conductivity</text>
      <text class="box-body" x="10" y="234">• Sparse thermal LSE assembly</text>
      <text class="box-body" x="10" y="254">• Updates temperature: <tspan class="box-code" fill="{PALETTE["charcoal"]}">T(n+1)</tspan></text>
    </g>
  </g>

  <!-- Connector: 3 -> 4 -->
  <g transform="translate(620, 735)">
    <line x1="0" y1="0" x2="0" y2="24" stroke="{PALETTE["slate_blue"]}" stroke-width="1.8" marker-end="url(#arrow-slate)" />
    <text class="arrow-text" x="10" y="15">transport and degassing</text>
  </g>

  <!-- SECTION 4: TRANSPORT & EVOLUTION -->
  <g transform="translate(60, 762)">
    <rect width="1120" height="136" rx="8" fill="{PALETTE["card_bg"]}" stroke="{PALETTE["border_gray"]}" stroke-width="1.2" />
    <path d="M 0 8 Q 0 0 8 0 L 1112 0 Q 1120 0 1120 8 L 1120 28 L 0 28 Z" fill="{PALETTE["header_bg"]}" />
    <line x1="0" y1="28" x2="1120" y2="28" stroke="{PALETTE["border_gray"]}" stroke-width="1" />
    <text class="section-title" x="16" y="19" fill="{PALETTE["purple"]}">4. POST-SOLVE TRANSPORT, DEGASSING AND DOMAIN ADVANCEMENT</text>
    
    <rect x="915" y="6" width="190" height="16" rx="4" fill="{PALETTE["purple"]}" fill-opacity="0.12" />
    <text class="tag-text" x="1010" y="18" text-anchor="middle" fill="{PALETTE["purple"]}">POST-SOLVE STAGE</text>

    <!-- Card 1: Volatiles & Atmosphere -->
    <g transform="translate(16, 38)">
      <rect width="348" height="86" rx="6" fill="{PALETTE["white"]}" stroke="{PALETTE["border_gray"]}" stroke-width="1" />
      <text class="box-title" x="12" y="20" fill="{PALETTE["purple"]}">Volatiles, Venting and Atmosphere</text>
      <text class="box-body" x="12" y="38">• Venting and degassing: <tspan class="box-code" fill="{PALETTE["charcoal"]}">vent_and_degas!</tspan></text>
      <text class="box-body" x="12" y="55">• Atmosphere evolution: <tspan class="box-code" fill="{PALETTE["charcoal"]}">evolve_atmosphere!</tspan></text>
      <text class="box-body" x="12" y="72">• Hydrodynamic boil-off and kinetic Jeans escape</text>
    </g>

    <!-- Card 2: Marker Advection -->
    <g transform="translate(386, 38)">
      <rect width="348" height="86" rx="6" fill="{PALETTE["white"]}" stroke="{PALETTE["border_gray"]}" stroke-width="1" />
      <text class="box-title" x="12" y="20" fill="{PALETTE["purple"]}">Marker Advection and Replenishment</text>
      <text class="box-body" x="12" y="38">• RK4 particle advection: <tspan class="box-code" fill="{PALETTE["charcoal"]}">advect_markers!</tspan></text>
      <text class="box-body" x="12" y="55">• Step-start pressure backtracking: <tspan class="box-code" fill="{PALETTE["charcoal"]}">Pt, Pf</tspan></text>
      <text class="box-body" x="12" y="72">• Cell replenishment: <tspan class="box-code" fill="{PALETTE["charcoal"]}">replenish!</tspan></text>
    </g>

    <!-- Card 3: Timestep & Telescoping -->
    <g transform="translate(756, 38)">
      <rect width="348" height="86" rx="6" fill="{PALETTE["white"]}" stroke="{PALETTE["border_gray"]}" stroke-width="1" />
      <text class="box-title" x="12" y="20" fill="{PALETTE["purple"]}">Adaptive Timestep and Telescoping</text>
      <text class="box-body" x="12" y="38">• Courant CFL timestep control on matrix and Darcy flow</text>
      <text class="box-body" x="12" y="55">• Metal segregation subcycling: <tspan class="box-code" fill="{PALETTE["charcoal"]}">Nsub ≤ Nmax</tspan></text>
      <text class="box-body" x="12" y="72">• Telescoping doubling: <tspan class="box-code" fill="{PALETTE["charcoal"]}">R(t) &gt; 0.70·(xsize/2)</tspan> at constant dx</text>
    </g>
  </g>

  <!-- Connector: 4 -> 5 -->
  <g transform="translate(620, 898)">
    <line x1="0" y1="0" x2="0" y2="24" stroke="{PALETTE["slate_blue"]}" stroke-width="1.8" marker-end="url(#arrow-slate)" />
  </g>

  <!-- SECTION 5: CHECKPOINTING & RECURRENCE -->
  <g transform="translate(60, 925)">
    <rect width="1120" height="66" rx="8" fill="{PALETTE["card_bg"]}" stroke="{PALETTE["slate_blue"]}" stroke-width="1.2" />
    
    <text class="box-title" x="20" y="26" fill="{PALETTE["slate_blue"]}">Time Advancement:</text>
    <text class="box-code" x="165" y="26" fill="{PALETTE["charcoal"]}">t ← t + dt,  step ← step + 1</text>
    
    <text class="box-title" x="420" y="26" fill="{PALETTE["slate_blue"]}">Data Persistence:</text>
    <text class="box-code" x="555" y="26" fill="{PALETTE["charcoal"]}">JLD2 Checkpoint Snapshots and Telemetry Streaming</text>

    <text class="box-body" x="20" y="50">• Terminates when <tspan class="box-code" fill="{PALETTE["charcoal"]}">t ≥ t_end</tspan> or maximum steps reached</text>
    <text class="box-body" x="420" y="50">• Zero numerical diffusion on moving boundaries via Lagrangian markers</text>
  </g>

  <!-- Feedback Loop Arrow: 5 -> 2 -->
  <g>
    <path d="M 60 958 L 30 958 L 30 321 L 54 321" fill="none" stroke="{PALETTE["slate_blue"]}" stroke-width="1.8" stroke-dasharray="6,4" marker-end="url(#arrow-slate)" />
    <text class="arrow-text" transform="translate(22, 640) rotate(-90)" text-anchor="middle">next time step (dt)</text>
  </g>

</svg>
"""

    with open(svg_path, "w", encoding="utf-8") as f:
        f.write(svg)
    print(f"Generated architecture SVG at: {svg_path}")


if __name__ == "__main__":
    generate_svg()

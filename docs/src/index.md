# Erebus.jl

```@raw html
<p align="center">
  <b>Two-Dimensional Thermo-Hydro-Mechanical Evolution of Early Planetesimals</b>
</p>
```

**Erebus.jl** is a high-performance numerical simulation code developed in Julia for modeling the coupled thermo-hydro-mechanical evolution of porous planetesimals in the early Solar System. It solves the combined equations of solid Stokes matrix deformation, Darcy pore fluid percolation, poroelastic compressibility, and radiogenic heat transport on a staggered finite-difference grid coupled with a Marker-in-Cell (MIC) advection framework.

---

## Key Physical Capabilities

- **Coupled Stokes-Darcy Flow**:
  Simultaneous solution of viscous/visco-elasto-plastic solid matrix flow and Darcy fluid transport through a permeable, deformable porous medium.

- **Poroelastic Compressibility**:
  Coupled solid matrix compressibility ($\beta_s$) and fluid compressibility ($\beta_f$) with dynamic drained compressibility ($\beta_d$), Biot-Willis coefficient ($K_{\text{BW}}$), and Skempton coefficient ($B$).

- **Terzaghi Effective Stress & Plasticity**:
  Mohr-Coulomb yielding criterion and tensile cut-off formulated in terms of Terzaghi effective stress, capturing compaction, faulting, and pore fluid overpressure regimes.

- **Radiogenic Decay Kinetics**:
  Time-dependent volumetric heating from short-lived radionuclides ($^{26}\text{Al}$ and $^{60}\text{Fe}$), driving water ice melting, hydrothermal circulation, and dehydration reactions.

- **Marker-in-Cell Advection**:
  Conservative transport of composition, temperature, melt fraction, and porosity on moving lagrangian markers interpolated onto the Eulerian staggered grid.

- **Reproducible TOML Configuration**:
  Declarative simulation parameters structured across grid, geometry, timestepping, solver controls, poroelasticity, thermodynamics, materials, and output storage.

---

## Documentation Structure

The documentation follows the Diataxis framework, separated into four distinct quadrants:

```@raw html
<div style="display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin-top: 20px;">
  <div style="border: 1px solid #e0e0e0; border-radius: 8px; padding: 16px;">
    <h3><a href="tutorials/quickstart/">Tutorials</a></h3>
    <p>Learning-oriented guides taking you through a complete simulation run and analytical benchmark verification.</p>
  </div>
  <div style="border: 1px solid #e0e0e0; border-radius: 8px; padding: 16px;">
    <h3><a href="howto/installation/">How-To Guides</a></h3>
    <p>Task-oriented instructions for configuring simulation setups, exploring parameter spaces, and post-processing output files.</p>
  </div>
  <div style="border: 1px solid #e0e0e0; border-radius: 8px; padding: 16px;">
    <h3><a href="explanations/model_overview/">Explanations</a></h3>
    <p>In-depth theoretical background on the governing Stokes-Darcy equations, poroelastic constitutive laws, and numerical algorithms.</p>
  </div>
  <div style="border: 1px solid #e0e0e0; border-radius: 8px; padding: 16px;">
    <h3><a href="reference/config_schema/">Reference</a></h3>
    <p>Comprehensive technical specifications for all TOML configuration options, submodule APIs, and literature citations.</p>
  </div>
</div>
```

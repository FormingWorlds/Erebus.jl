# Scientific Context: Early Planetesimal Evolution

`Erebus.jl` models the internal physical and geochemical evolution of icy and rocky planetesimals during the first tens of millions of years of Solar System history.

---

## The Physical Scenario

Planetesimals formed within the solar protoplanetary disk through the gravitational collapse of pebble swarms or streaming instabilities. Bodies that accreted within the first 2 to 3 million years after the formation of Calcium-Aluminum-rich Inclusions (CAIs) incorporated significant quantities of short-lived radioactive isotopes, predominantly $^{26}\text{Al}$ (half-life $t_{1/2} \approx 0.717\text{ Ma}$) and $^{60}\text{Fe}$ ($t_{1/2} \approx 2.62\text{ Ma}$).

As these radionuclides decayed, they generated intense internal volumetric heating:

1. **Ice Melting**:
   Initial temperatures of $\approx 170\text{ K}$ rose rapidly. Once internal temperatures reached $273\text{ K}$, primordial water ice melted, absorbing latent heat ($L^f \approx 333\text{ kJ/kg}$) and generating mobile pore water within the porous silicate matrix.

2. **Hydrothermal Circulation**:
   Liquid water percolated through the interconnected pore network under Darcy flow, driven by thermal buoyancy and compaction-induced fluid pressure gradients.

3. **Clay Dehydration & Gas Generation**:
   At elevated temperatures ($T \approx 500\text{--}900\text{ K}$), hydrous phyllosilicates (such as serpentine) underwent thermal dehydration reactions, releasing bound hydroxyl groups as free supercritical fluid.

4. **Matrix Compaction and Pore Overpressure**:
   Viscous, elastic, and plastic deformation of the silicate matrix reduced porosity ($\phi$). When the compaction rate exceeded the rate of fluid escape via Darcy percolation, pore fluid pressure ($P_f$) rose to match or exceed the lithostatic solid pressure ($P_t$).

5. **Failure & Hydrofracture**:
   When the Terzaghi effective stress became tensile ($\sigma_{\text{eff}} \le -\sigma_{\text{tensile}}$), the silicate matrix experienced hydrofracturing, creating high-permeability pathways that vented fluids to the planetesimal surface.

---

## Thermo-Hydro-Mechanical Coupling

`Erebus.jl` resolves these interacting regimes through a fully coupled numerical framework:
- **Thermal Solver**: Computes conduction, radioactive decay, latent heat buffering, and fluid advective heat transport.
- **Hydromechanical Solver**: Solves coupled Stokes solid deformation, Darcy fluid flux, and poroelastic volume changes simultaneously in a monolithic linear system.
- **Marker-in-Cell (MIC)**: Advects material phases, temperature, composition, and porosity without numerical diffusion across moving boundaries.

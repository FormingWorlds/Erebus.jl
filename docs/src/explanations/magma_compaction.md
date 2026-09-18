# Buoyant Silicate Melt Percolation and Matrix Compaction

This section documents the physical theory, governing continuum equations, and numerical implementation of buoyant silicate melt migration through a compacting solid rock matrix in `Erebus.jl`.

---

## Physical Background and Motivation

Silicate partial melting in planetesimals and planetary embryos produces buoyant liquid magma.
Because liquid silicate melt is less dense than the crystalline silicate matrix ($\rho_m < \rho_s$), positive buoyancy drives the liquid upward toward the planetary surface.
As melt percolates through permeable crystalline grain networks, the solid rock matrix must compact to accommodate the departure of melt.
Conversely, regions where melt accumulates must expand or generate dynamic compaction overpressure.

McKenzie (1984) formulated the fundamental two-phase continuum equations governing the flow of melt through a deformable, viscous solid matrix.
In contrast to passive Darcy percolation through rigid porous media, matrix compaction couples the liquid pressure field directly to matrix deformation.
This coupling produces dynamic compaction pressure gradients, determines the characteristic compaction length scale, and governs the propagation of non-linear solitary porosity waves.

---

## Governing Continuum Equations

### Two-Phase Conservation Laws

Let $\phi_m$ be the volume fraction of liquid silicate melt (porosity), $\mathbf{v}_m$ the liquid melt velocity, and $\mathbf{v}_s$ the solid matrix velocity.
Assuming constant solid and liquid phase densities $\rho_s$ and $\rho_m$, conservation of mass for each phase yields:

$$\frac{\partial (1 - \phi_m)}{\partial t} + \nabla \cdot \left( (1 - \phi_m) \mathbf{v}_s \right) = -\Gamma_m$$

$$\frac{\partial \phi_m}{\partial t} + \nabla \cdot (\phi_m \mathbf{v}_m) = \Gamma_m$$

where $\Gamma_m$ is the net volumetric melting rate.
Summing both mass conservation equations gives the total mixture continuity condition:

$$\nabla \cdot \left( (1 - \phi_m) \mathbf{v}_s + \phi_m \mathbf{v}_m \right) = 0$$

### Darcy-McKenzie Momentum Balance

Melt filtration relative to the deforming solid matrix satisfies Darcy's law:

$$\phi_m (\mathbf{v}_m - \mathbf{v}_s) = -\frac{k_m}{\eta_m} (\nabla P_m - \rho_m \mathbf{g})$$

where:
- Permeability $k_m$ is the effective silicate melt permeability [$\text{m}^2$].
- Dynamic viscosity $\eta_m$ is the liquid silicate viscosity [$\text{Pa}\cdot\text{s}$].
- Melt pressure $P_m$ is the liquid melt pressure [$\text{Pa}$].
- Gravitational acceleration $\mathbf{g}$ is the gravity vector [$\text{m}/\text{s}^2$].

In a deformable compacting medium, the liquid melt pressure $P_m$ differs from the mean solid matrix pressure $P_s$ by the dynamic compaction pressure $P_{\text{comp}}$:

$$P_s - P_m = P_{\text{comp}} = -\zeta_m (\nabla \cdot \mathbf{v}_s)$$

where $\zeta_m$ is the effective bulk viscosity of the solid rock matrix [$\text{Pa}\cdot\text{s}$].

Under lithostatic equilibrium ($\nabla P_s \approx \rho_s \mathbf{g}$), the effective driving force for melt migration is:

$$\mathbf{f}_{\text{drive}} = (\rho_s - \rho_m) \mathbf{g} - \nabla P_{\text{comp}}$$

When matrix compaction resistance is negligible ($\zeta_m \to 0$), melt migration reduces to pure buoyant percolation.
When melt converges beneath an impermeable cold lid, dynamic compaction overpressure builds up, opposing buoyancy and halting upward ascent.

---

## Constitutive Relationships

### Matrix Bulk Viscosity ($\zeta_m$)

The resistance of the porous crystalline matrix to volumetric compaction or dilation scales inversely with melt fraction (McKenzie 1984; Scott & Stevenson 1986):

$$\zeta_m = \xi_{\text{bulk}} \frac{\eta_s}{\max(\phi_m, \phi_{\min})}$$

where:
- Shear viscosity $\eta_s$ is the shear viscosity of the solid rock matrix [$\text{Pa}\cdot\text{s}$].
- Viscosity ratio $\xi_{\text{bulk}}$ is the bulk-to-shear viscosity ratio (`bulk_viscosity_ratio`, default 1.0).
- Regularization threshold $\phi_{\min}$ is a numerical threshold (`min_bulk_porosity`, default 0.005) that prevents infinite bulk viscosity in subsolidus rock.

### Compaction Length Scale ($\delta_c$)

The compaction length $\delta_c$ defines the characteristic spatial scale over which matrix compaction balances Darcy resistance:

$$\delta_c = \sqrt{\frac{\zeta_m + \frac{4}{3}\eta_s}{\eta_m} k_m}$$

In `Erebus.jl`, the compaction length is computed locally and clamped within user-configurable bounds (`compaction_length_min` and `compaction_length_max`):

$$\delta_c = \text{clamp}\left( \sqrt{\frac{\zeta_m + \frac{4}{3}\eta_s}{\eta_m} k_m}, \, \delta_{c,\min}, \, \delta_{c,\max} \right)$$

Physical scales for planetesimal silicate interiors range from tens of meters in crystal-rich mush to tens of kilometers in partially molten mantles.

### Silicate Melt Permeability ($k_m$)

Permeable channel flow follows the McKenzie (1984) power-law formulation with residual melt retention threshold:

$$k_m(\phi_m) = \begin{cases}
0, & \phi_m \le \phi_{\text{residual}} \\
k_0 \left( \frac{\phi_m - \phi_{\text{residual}}}{\phi_0 - \phi_{\text{residual}}} \right)^n, & \phi_{\text{residual}} < \phi_m \le \phi_{\text{crit}} \\
k_{\text{crit}} + k_0 \left( \frac{\phi_m - \phi_{\text{crit}}}{\phi_0 - \phi_{\text{residual}}} \right), & \phi_m > \phi_{\text{crit}}
\end{cases}$$

where $k_0$ is the reference permeability at reference porosity $\phi_0$, $n$ is the power-law exponent ($n \approx 2\text{ to }3$), $\phi_{\text{residual}}$ is the residual melt fraction retained in pores by surface tension, and $\phi_{\text{crit}} \approx 0.40$ is the critical disaggregation threshold.

---

## Crustal Ponding and Volcanic Eruptions

### Subsolidus Sill Ponding

When buoyant melt ascends into cold crust where local temperature falls below the silicate solidus ($T < T_{\text{solidus}}$), permeability collapses and matrix viscosity increases by orders of magnitude.
If `ponding_active = true`, radial melt flux at subsolidus boundaries is blocked.
Silicate melt accumulates beneath cold, subsolidus lithospheric layers, simulating the formation of sub-crustal magma sills.

### Overpressure Eruption (`eruption_active`)

When compaction overpressure in a ponded melt accumulation exceeds rock tensile strength $\sigma_t$ (`tensile_strength`, default $10\text{ MPa}$):

$$\max(P_{\text{comp}}) > \sigma_t$$

When this criterion is met, tensile hydrofracture conduits penetrate the subsolidus crustal lid, permitting magma ascent through the lid layer.

---

## Decompression Volatile Exsolution

As hydrous silicate melt ascends toward shallower depths, lithostatic pressure decreases.
Volatile species dissolved in the melt (water, carbon monoxide, carbon dioxide, methane, nitrogen, and sulfur) exceed their equilibrium solubility limits.
When `exsolution_active = true`, decompression drives multi-species volatile exsolution into pore gas or supercritical fluid:

$$\Delta w_v = w_v(P_{\text{deep}}) - w_{v,\text{eq}}(P_{\text{shallow}}, T)$$

Exsolved volatiles enter the marker pore fluid fraction $\phi_m$, increasing fluid pressure and driving cold surface venting or atmospheric accumulation.

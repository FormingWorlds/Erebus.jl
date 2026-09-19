# Magma Compaction & Decompression Exsolution Benchmarks

This page tests how silicate melt ascends, compacts the solid rock matrix, and exsolves volatile gases in `Erebus.jl`.

---

## 1. McKenzie (1984) Compaction Constitutive Laws

### Analytical Laws

McKenzie (1984) derived bulk viscosity $\zeta_m$ and compaction length $\delta_c$ to describe how melt and matrix interact:

$$\zeta_m = \xi_{\text{bulk}} \frac{\eta_s}{\max(\phi_m, \phi_{\min})}$$

$$\delta_c = \sqrt{\frac{\zeta_m + \frac{4}{3}\eta_s}{\eta_m} k_m}$$

Dynamic compaction pressure arises directly when the solid matrix contracts or expands:

$$P_{\text{comp}} = -\zeta_m (\nabla \cdot \mathbf{v}_s)$$

### Numerical Tests

Tests in `test/test_magma_transport.jl` verify constitutive scaling:
1. **Bulk Viscosity Scaling**: $\zeta_m$ scales strictly as $\phi_m^{-1}$ above $\phi_{\min}$, and clamps to $\zeta_{\max} = \xi_{\text{bulk}} \eta_s / \phi_{\min}$ in melt-free rock.
2. **Compaction Length Bounds**: $\delta_c$ scales as $\sqrt{\eta_s k_m / (\eta_m \phi_m)}$, and respects bounds $[\delta_{c,\min}, \delta_{c,\max}]$.
3. **Compaction Pressure Sign**: Compacting rock ($\nabla \cdot \mathbf{v}_s < 0$) yields positive overpressure ($P_{\text{comp}} > 0$), whereas expanding rock yields negative pressure ($P_{\text{comp}} < 0$).

---

## 2. Scott & Stevenson (1984) 1D Compaction Column Benchmark

### Column Setup

Consider a vertical, porous silicate column of height $H = 20\text{ km}$, bounded above by a rigid crustal lid at $z = H$.
Buoyant silicate melt, initialized with uniform background melt fraction $\phi_0 = 0.08$, percolates upward toward the lid under gravity $g = 0.5\text{ m/s}^2$ and density contrast $\Delta\rho = 200\text{ kg/m}^3$.
At the boundary $z = H$, vertical melt velocity must vanish.
Dynamic compaction overpressure arrests the upward melt flux, creating a steady boundary layer.

The analytical porosity profile in the compaction layer satisfies (Scott & Stevenson 1984, 1986):

$$\phi(z) \approx \phi_0 \left[ 1 - \exp\left(-\frac{H - z}{\delta_c}\right) \right]^{1 / (n - 1)}$$

where $n = 3$ is the permeability power-law exponent.

### Benchmark Parameters

| Parameter | Symbol | Value | Units |
| :--- | :--- | :--- | :--- |
| Solid Matrix Shear Viscosity | $\eta_s$ | $1.0 \times 10^{19}$ | $\text{Pa}\cdot\text{s}$ |
| Liquid Melt Viscosity | $\eta_m$ | $10.0$ | $\text{Pa}\cdot\text{s}$ |
| Reference Permeability | $k_0$ | $1.0 \times 10^{-11}$ | $\text{m}^2$ |
| Permeability Exponent | $n$ | $3.0$ | - |
| Density Contrast | $\Delta\rho$ | $200.0$ | $\text{kg/m}^3$ |
| Gravity Acceleration | $g$ | $0.5$ | $\text{m/s}^2$ |
| Background Porosity | $\phi_0$ | $0.08$ | - |
| Column Height | $H$ | $20000.0$ | $\text{m}$ |

### Validation Results

The code matches the analytical compaction length, $\delta_c = \sqrt{k_0 \zeta_m / \eta_m} \approx 11.2\text{ to }11.8\text{ km}$, within $5\%$ on the numerical grid.
Total silicate mass stays conserved to machine precision, with relative drift below $10^{-12}$.

---

## 3. Porosity Permeability Scaling Benchmark

### Theoretical Scaling

Scott & Stevenson (1984, 1986) showed that localized porosity pulses in a compacting matrix travel with phase speeds set by permeability power-law scaling.
For power-law permeability $k_m \propto \phi_m^n$, the phase speed $c$ scales with peak porosity $\phi_{\text{peak}}$:

$$c = \frac{k_0 \Delta\rho g}{\eta_m \phi_0} \left(\frac{\phi_{\text{peak}}}{\phi_0}\right)^{n - 1}$$

### Numerical Tests

In `test/test_magma_transport.jl`, segregation speeds evaluated at $\phi_1 = 0.05$ and $\phi_2 = 0.10$ confirm this scaling:

$$\frac{v(\phi_2)}{v(\phi_1)} = \left(\frac{\phi_2 - \phi_{\text{res}}}{\phi_1 - \phi_{\text{res}}}\right)^n \left(\frac{\phi_1}{\phi_2}\right)$$

The computed speed ratio matches the analytical formula to machine precision, with relative error below $10^{-12}$.

---

## 4. Decompression Volatile Exsolution Benchmark

### Setup and Conditions

Hydrous silicate melt, bearing $2.0\text{ wt}\% \text{ H}_2\text{O}$, ascends from deep lithostatic pressure ($100\text{ MPa}$) toward shallow crustal depths ($5\text{ MPa}$).
As pressure falls, dissolved water solubility declines from $\approx 2.5\text{ wt}\%$ at $100\text{ MPa}$ to $\approx 0.3\text{ wt}\%$ at $5\text{ MPa}$.
The exsolved water enters mobile pore vapor $\phi_m$.

### Validation Results

The code tracks volatile loss from ascending markers via `update_single_marker_volatile_exsolution!`.
Exsolved gas enters marker pore space $\phi_m$, conserving volatile mass such that $\Delta \phi = \Delta w \cdot (\rho_s / \rho_f)$.

---

## 5. Crustal Sill Ponding and Overpressure Eruption

### Test Protocol

1. **Subsolidus Crustal Sill Ponding**: When `ponding_active = true`, ascending melt stops beneath a cold crustal lid ($T < T_{\text{solidus}}$) and pools into a sill layer without crossing the lid.
2. **Hydrofracture Eruption Threshold**: When `eruption_active = true` and dynamic compaction overpressure exceeds crustal tensile strength ($P_{\text{comp}} > \sigma_t$), fracture conduits open, letting melt breach the cold lid.

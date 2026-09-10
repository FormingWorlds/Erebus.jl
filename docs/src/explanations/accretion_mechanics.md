# Planetesimal Accretion Mechanics and Boundary Evolution

This section explains the physical principles, orbital dynamics, thermodynamic formulations, and numerical implementation of planetesimal accretion in `Erebus.jl`.

---

## Planetesimal Accretion Regimes

Planetesimals grow in circumstellar disks by capturing pebbles from the gas disk or by colliding with other planetesimals in a collision swarm:

### 1. Pebble Accretion

Pebble accretion occurs when small aerodynamic particles (pebbles with Stokes numbers $\mathrm{St} \sim 10^{-3} - 1.0$) drift inward through the circumstellar gas disk. When a pebble enters the gravitational sphere of influence of a planetesimal, aerodynamic drag dissipates its orbital kinetic energy, causing it to spiral onto the planetesimal:

- **Bondi Pebble Accretion ($M < M_{\mathrm{trans}}$):** For small bodies, gas pressure gradients and sub-Keplerian headwind dominate the relative approach velocity $\Delta v$. The Bondi radius $R_B = G M / c_s^2$ defines the scale where gravitational attraction matches gas thermal energy. The capture radius is bounded by the Hill radius, $r_{\mathrm{acc},B} = \min(R_H, 2 \sqrt{(\mathrm{St}/\Omega_K)(GM/\Delta v)})$.
- **Hill Pebble Accretion ($M \ge M_{\mathrm{trans}}$):** For larger bodies, the body's Hill sphere $R_H = a(M / 3M_\star)^{1/3}$ exceeds the Bondi sphere. Stellar tidal shear governs the encounter velocity $v_H = \Omega_K R_H$, and the capture radius expands to $r_{\mathrm{acc},H} = R_H \mathrm{St}^{1/3}$.
- **Regime Transition:** The transition mass $M_{\mathrm{trans}} = \sqrt{1/3} \Delta v^3 / (G \Omega_K)$ marks where Bondi and Hill regimes transition (Lambrechts & Johansen, 2012). In `Erebus.jl`, the `:pebble_auto` mode selects the Bondi regime for $M < M_{\mathrm{trans}}$ and the Hill regime for $M \ge M_{\mathrm{trans}}$.

### 2. Safronov Gravitational Focusing

In gas-poor environments or after disk clearing, planetesimal swarms grow through mutual inelastic collisions. The Safronov number $\Theta = v_{\mathrm{esc}}^2 / (2 \sigma_v^2)$ measures the ratio of escape velocity to velocity dispersion. Gravitational focusing expands the geometric cross-section by the factor $F_g = 1 + 2\Theta$, driving rapid runaway or oligarchic growth when velocity dispersion remains low.

---

## Expanding Boundary Mechanics on Marker-in-Cell Grids

`Erebus.jl` discretizes continuous multiphase rock and fluid mechanics using a staggered Eulerian grid coupled with Lagrangian markers:

1. **Sticky-Air Representation:** Space surrounding the planetesimal is represented as "sticky air" markers ($tm = 3$). These markers possess low viscosity and zero radiogenic heat production, enforcing a traction-free internal free surface at the planetesimal-space interface.
2. **Radial Boundary Advance:** When mass increments $\Delta M = \dot{M} \, \Delta t$ occur, the planetesimal radius expands by:
   $$\Delta R = R \left( \left( 1 + \frac{3\,\Delta M}{4\pi\,\rho_{\mathrm{bulk}}\,R^3} \right)^{1/3} - 1 \right)$$
   All sticky-air markers inside the expanded shell $r \le R + \Delta R$ convert to solid crust markers ($tm = 2$).
3. **Property Initialization:** Converted markers receive the accreted shell porosity $\phi_{\mathrm{accreted}}$, bulk iron fraction $X_{\mathrm{fe,bulk}}$, volatile fractions ($\mathrm{H}_2\mathrm{O}, \mathrm{C}, \mathrm{N}, \mathrm{S}$), and impact-heated temperature $T_{\mathrm{accreted}}$.
4. **Time Stamp Inheritance:** Each converted marker receives a timestamp $t_{\mathrm{acc}} = t$. This preserves the accretion epoch for post-processing and diagnostic analysis.

---

## Impact Heating Thermodynamics

Accreting projectiles deliver kinetic energy to the planetesimal surface:

$$u_{\mathrm{acc}} = \frac{G M}{R} + \frac{1}{2} v_\infty^2$$

A fraction $h_{\mathrm{impact}} \in [0, 1]$ is retained as internal heat, while the remaining fraction radiates to space:

$$\Delta T_{\mathrm{impact}} = \frac{h_{\mathrm{impact}} \, u_{\mathrm{acc}}}{c_p}$$

This heating establishes a warm surface shell during rapid accretion phases, driving early near-surface volatile exsolution and sintering porous regolith.

---

## Snowline Migration and Volatile Budget

As circumstellar gas cools or planetesimals migrate, the local temperature at semi-major axis $a$ determines whether volatile species condense as ices:

- **Outside Water Snowline ($T_{\mathrm{disk}} \le 160\text{ K}$):** Converted markers receive hydrated rock ($XW = 0.40$) and high bulk water content ($10\text{ wt\%}$), representing carbonaceous chondrite precursors.
- **Inside Water Snowline ($T_{\mathrm{disk}} > 160\text{ K}$):** Converted markers receive dry anhydrous rock ($XW = 0.0, 0.1\text{ wt\% H}_2\mathrm{O}$), representing enstatite or ordinary chondrite precursors.

---

## Radiogenic Onion-Shell Thermal Structure

Short-lived $^{26}\mathrm{Al}$ decays with mean lifetime $\tau = 1.035\text{ Myr}$ ($t_{1/2} \approx 0.717\text{ Myr}$):

$$Q(t) = Q_0 \, e^{-t / \tau}$$

Because early-accreted interior material experiences heating when $^{26}\mathrm{Al}$ concentrations are maximal ($t < 1\text{ Myr}$), the planetesimal center undergoes extensive heating, melting, and core segregation. Outer crustal shells accreted at later epochs ($t > 1.5\text{ Myr}$) inherit significantly decayed radionuclide inventories, remaining cold, unmelted, and volatile-rich. This mechanism naturally produces the concentric metamorphic zoning observed in chondritic meteorites.

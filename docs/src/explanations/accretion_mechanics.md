# Planetesimal Accretion Mechanics and Boundary Evolution

This section explains the physical principles, orbital dynamics, thermodynamic formulations, and numerical implementation of planetesimal accretion in `Erebus.jl`.

For numerical benchmark comparisons, scaling figures, and verification test suites, see [Planetesimal Accretion & Impact Heating](../validation/planetesimal_accretion.md).

---

## Orbital Kinematics and Disk Structure

At heliocentric semi-major axis $a$ around a central star of mass $M_\star$, the Keplerian orbital frequency and circular velocity are:

$$\Omega_K = \sqrt{\frac{G M_\star}{a^3}}, \quad v_K = \Omega_K \, a = \sqrt{\frac{G M_\star}{a}}$$

Circumstellar gas at temperature $T_{\mathrm{gas}}$, adiabatic index $\gamma = 1.0$ (isothermal disk gas), and mean molecular weight $\mu = 2.34$ has isothermal sound speed and pressure scale height:

$$c_s = \sqrt{\frac{\gamma \, k_B \, T_{\mathrm{gas}}}{\mu \, m_p}}, \quad H_g = \frac{c_s}{\Omega_K}$$

Gas drag causes solid pebbles with dimensionless Stokes number $\mathrm{St}$ to settle toward the disk midplane against turbulent diffusion parameterized by $\alpha_{\mathrm{turb}}$ (Youdin and Lithwick, 2007):

$$H_{\mathrm{peb}} = H_g \sqrt{\frac{\alpha_{\mathrm{turb}}}{\alpha_{\mathrm{turb}} + \mathrm{St}}}$$

In the strongly coupled limit ($\mathrm{St} \ll \alpha_{\mathrm{turb}}$), pebbles stay well-mixed across the gas column ($H_{\mathrm{peb}} \to H_g$). For decoupled pebbles ($\mathrm{St} \gg \alpha_{\mathrm{turb}}$), pebbles settle into a thin midplane layer ($H_{\mathrm{peb}} \ll H_g$).

The pebble surface density follows a radial power law relative to reference radius $a_0 = 1\text{ AU}$:

$$\Sigma_{\mathrm{peb}}(a) = \Sigma_{\mathrm{peb},0} \left(\frac{a}{a_0}\right)^{-p_{\mathrm{peb}}}$$

The corresponding midplane pebble volume density is $\rho_{\mathrm{peb}} = \Sigma_{\mathrm{peb}} / (\sqrt{2\pi} H_{\mathrm{peb}})$.

---

## Planetesimal Accretion Regimes

Planetesimals grow by capturing pebbles from the circumstellar disk or through mutual collisions in planetesimal swarms:

### 1. Pebble Accretion

Aerodynamic drag dissipates the kinetic energy of small particles ($\mathrm{St} \sim 10^{-3} - 1.0$) inside the gravitational sphere of influence of the planetesimal, causing pebbles to spiral onto the surface (Ormel and Klahr, 2010; Lambrechts and Johansen, 2012).

The characteristic interaction lengths are the Bondi radius $R_B$ and Hill radius $R_H$:

$$R_B = \frac{G M}{c_s^2}, \quad R_H = a \left(\frac{M}{3 M_\star}\right)^{1/3}$$

The transition between the Bondi and Hill regimes occurs at the transition mass:

$$M_{\mathrm{trans}} = \sqrt{\frac{1}{3}} \, \frac{\Delta v^3}{G \Omega_K}$$

where the sub-Keplerian gas headwind velocity is $\Delta v = \eta v_K \approx 1.5 (c_s / v_K)^2 v_K$.

#### Bondi Regime ($M < M_{\mathrm{trans}}$)

For small planetesimals, gas headwind dominates the approach velocity. The Bondi capture radius is bounded by the Hill sphere:

$$r_{\mathrm{acc},B} = \min\left(R_H, \, 2 \sqrt{\frac{\mathrm{St}}{\Omega_K} \frac{G M}{\Delta v}}\right)$$

The mass accretion rate evaluates from either 2D sheet accretion or 3D spherical capture depending on whether $r_{\mathrm{acc},B} \gtrless H_{\mathrm{peb}}$:

$$\dot{M}_{B,\mathrm{2D}} = 2 \, r_{\mathrm{acc},B} \, \Sigma_{\mathrm{peb}} \, \Delta v, \quad \dot{M}_{B,\mathrm{3D}} = \pi \, r_{\mathrm{acc},B}^2 \, \rho_{\mathrm{peb}} \, \Delta v$$

#### Hill Regime ($M \ge M_{\mathrm{trans}}$)

For larger bodies, stellar tidal shear governs the encounter velocity $v_H = \Omega_K R_H$, and the capture radius expands:

$$r_{\mathrm{acc},H} = R_H \, \mathrm{St}^{1/3}$$

The corresponding 2D and 3D accretion rates are:

$$\dot{M}_{H,\mathrm{2D}} = 2 \, r_{\mathrm{acc},H} \, \Sigma_{\mathrm{peb}} \, (\Omega_K R_H), \quad \dot{M}_{H,\mathrm{3D}} = \pi \, r_{\mathrm{acc},H}^2 \, \rho_{\mathrm{peb}} \, (\Omega_K R_H)$$

In `:pebble_auto` mode, the solver evaluates $M_{\mathrm{trans}}$ dynamically and switches between the Bondi and Hill formulations.

### 2. Safronov Gravitational Focusing

In gas-depleted disks or planetesimal swarms with surface density $\Sigma_{\mathrm{pl}}$ and velocity dispersion $\sigma_v$, gravitational focusing enhances the geometric cross-section (Safronov, 1972):

$$\Theta = \frac{v_{\mathrm{esc}}^2}{2 \sigma_v^2} = \frac{G M}{R \, \sigma_v^2}, \quad F_g = 1 + 2\Theta$$

$$\dot{M}_{\mathrm{saf}} = \pi \, R^2 \, \Sigma_{\mathrm{pl}} \, \Omega_K \, F_g$$

When $\sigma_v \gg v_{\mathrm{esc}}$, $\Theta \to 0$ and $F_g \to 1$ (geometric limit). When $\sigma_v \ll v_{\mathrm{esc}}$, $\Theta \gg 1$ and $F_g \approx 2\Theta$, producing runaway gravitational focusing.

### 3. Analytical Growth Modes

`Erebus.jl` also provides three analytical growth modes:
- **Constant mass rate (`:constant_rate`):** $\dot{M} = \dot{M}_0$.
- **Linear radius expansion (`:linear_radius`):** $\dot{R} = \dot{R}_0 \implies \dot{M} = 4\pi R^2 \rho_{\mathrm{bulk}} \dot{R}_0$.
- **Exponential growth (`:exponential`):** $\dot{M} = M / \tau_{\mathrm{growth}}$.

Accretion occurs within the time window $t \in [t_{\mathrm{start}}, t_{\mathrm{start}} + t_{\mathrm{duration}}]$, and shuts off when $M \ge M_{\mathrm{target}}$ or $R \ge R_{\mathrm{target}}$.

---

## Expanding Boundary Mechanics on Marker Grids

`Erebus.jl` tracks the expanding planetesimal radius on a staggered Eulerian grid coupled with Lagrangian markers:

1. **Sticky-Air Representation:** Space surrounding the body is populated with low-viscosity sticky-air markers ($tm = 3$).
2. **Exact 3D Spherical Volume Mapping:** A mass increment $\Delta M = \dot{M} \, \Delta t$ expands the planetesimal radius by:
   $$\Delta R = \left( R^3 + \frac{3\,\Delta M}{4\pi\,\rho_{\mathrm{bulk}}} \right)^{1/3} - R$$
   This formula conserves 3D spherical volume $\Delta V = \Delta M / \rho_{\mathrm{bulk}}$ without accumulating first-order truncation error. All sticky-air markers in the shell $r \le R + \Delta R$ convert to solid crust markers ($tm = 2$).
3. **Property Initialization:** Converted markers receive accreted shell porosity $\phi_{\mathrm{accreted}}$, bulk iron fraction $X_{\mathrm{fe,bulk}}$, volatile mass fractions ($\mathrm{H}_2\mathrm{O}, \mathrm{C}, \mathrm{N}, \mathrm{S}$), and impact temperature $T_{\mathrm{accreted}}$.
4. **Accretion Epoch Tracking:** Converted markers store their accretion timestamp $t_{\mathrm{acc}} = t$.

---

## Impact Heating Thermodynamics

Projectiles strike the surface with specific kinetic energy:

$$u_{\mathrm{acc}} = \frac{G M}{R} + \frac{1}{2} v_\infty^2$$

A fraction $h_{\mathrm{impact}} \in [0, 1]$ is retained as internal heat, while the remainder radiates to space:

$$\Delta T_{\mathrm{impact}} = \frac{h_{\mathrm{impact}} \, u_{\mathrm{acc}}}{c_p}$$

$$T_{\mathrm{accreted}} = T_{\mathrm{ambient}} + \Delta T_{\mathrm{impact}}$$

This heating establishes a warm outer boundary layer during rapid accretion stages.

---

## Disk Snowline Volatile Coupling

Disk temperature at heliocentric distance $a$ determines the volatile budget of newly accreted material:

- **Outside Water Snowline ($T_{\mathrm{disk}} \le 160\text{ K}$):** Converted markers receive hydrated rock ($XW = 0.40$) with $10\text{ wt\%}$ bulk $\mathrm{H}_2\mathrm{O}$, representing carbonaceous chondrite precursors.
- **Inside Water Snowline ($T_{\mathrm{disk}} > 160\text{ K}$):** Converted markers receive dry anhydrous rock ($XW = 0.0$, $0.1\text{ wt\% H}_2\mathrm{O}$), representing ordinary and enstatite chondrite precursors.

---

## Radiogenic Onion-Shell Thermal Structure

Short-lived $^{26}\mathrm{Al}$ decays with mean lifetime $\tau = 1.035\text{ Myr}$ ($t_{1/2} \approx 0.717\text{ Myr}$):

$$Q(t) = Q_0 \, e^{-t / \tau}$$

Because early-accreted interior material experiences peak radiogenic heating ($t < 1\text{ Myr}$), the deep interior undergoes melting and core segregation. Outer crustal shells accreted at later epochs ($t > 1.5\text{ Myr}$) inherit decayed radionuclide inventories and remain cold, unmelted, and volatile-rich, reproducing the concentric metamorphic zoning observed in chondritic parent bodies.

---

## Multi-Stage Accretion Transitions: Onset and Isolation Mass

In protoplanetary disks, growth from small planetesimals to planetary embryos spans three distinct physical regimes separated by characteristic mass thresholds:

### 1. Sub-Keplerian Gas Headwind
Gas pressure support reduces the azimuthal velocity of disk gas below the circular Keplerian velocity $v_K = \sqrt{G M_\star / a}$:

$$v_{\text{gas}} = v_K (1 - \eta)$$

$$\eta = -\frac{1}{2} \left(\frac{c_s}{v_K}\right)^2 \frac{\partial \ln P}{\partial \ln r} \approx 1.5 \left(\frac{c_s}{v_K}\right)^2$$

Keplerian solid bodies encounter this sub-Keplerian gas as a continuous headwind with relative velocity:

$$v_{\text{hw}} = \eta v_K = 1.5 \frac{c_s^2}{v_K}$$

### 2. Pebble Accretion Onset Mass ($M_{\text{onset}}$)
Pebble accretion operates in the settling regime when a pebble entering the gravitational capture radius settles onto the planetesimal within one stopping time $t_s = \tau_s / \Omega_K$, rather than being swept past by the headwind. Equating the headwind Bondi radius $R_B = G M / v_{\text{hw}}^2$ to the drift deflection distance $v_{\text{hw}} t_s$ yields the settling onset threshold (Visser & Ormel 2016):

$$M_{\text{onset}} = f_{\text{onset}} \frac{v_{\text{hw}}^3 \tau_s}{G \Omega_K}$$

where $f_{\text{onset}} = 1.0$ is the hydrodynamic calibration factor, $\tau_s$ is the dimensionless Stokes number, and $\Omega_K = \sqrt{G M_\star / a^3}$ is the orbital frequency. For $M < M_{\text{onset}}$, gravity is insufficient to overcome gas drag, and accretion proceeds via pairwise planetesimal collisions (Safronov 1972).

### 3. Pebble Isolation Mass ($M_{\text{iso}}$)
When the growing protoplanet perturbs the surrounding gas disk, spiral density waves carve a partial gap, creating an exterior pressure bump where $\partial P / \partial r > 0$. Inward-drifting pebbles become trapped at this pressure maximum, cutting off further pebble accretion onto the central embryo (Lambrechts et al. 2014; Bitsch et al. 2018):

$$M_{\text{iso}} = f_{\text{iso}} M_\star \left(\frac{H_g}{a}\right)^3 = f_{\text{iso}} M_\star \left(\frac{c_s}{v_K}\right)^3$$

where $H_g = c_s / \Omega_K$ is the disk scale height, $M_\star$ is stellar mass, and $f_{\text{iso}} \approx 0.5$ matches the standard hydrodynamic threshold of $\approx 20 M_\oplus (H_g / 0.05 a)^3$. Beyond $M_{\text{iso}}$, planet growth transitions to giant impact collisions between isolated embryos.

---

## References

- Bitsch, B., Morbidelli, A., Johansen, A., Lega, E., Lambrechts, M., & Crida, A. (2018). Pebble-isolation mass: Scaling law and implications for the formation of super-Earths and gas giants. *Astronomy & Astrophysics*, 612, A30. [https://doi.org/10.1051/0004-6361/201731931](https://doi.org/10.1051/0004-6361/201731931)
- Lambrechts, M., & Johansen, A. (2012). Rapid growth of gas-giant cores by pebble accretion. *Astronomy & Astrophysics*, 544, A32. [https://doi.org/10.1051/0004-6361/201219127](https://doi.org/10.1051/0004-6361/201219127)
- Lambrechts, M., Johansen, A., & Morbidelli, A. (2014). Separating gas-giant and ice-giant planets by halting pebble accretion. *Astronomy & Astrophysics*, 572, A35. [https://doi.org/10.1051/0004-6361/201423814](https://doi.org/10.1051/0004-6361/201423814)
- Ormel, C. W., & Klahr, H. H. (2010). The effect of gas drag on the growth of protoplanets. *Astronomy & Astrophysics*, 520, A43. [https://doi.org/10.1051/0004-6361/201014903](https://doi.org/10.1051/0004-6361/201014903)
- Safronov, V. S. (1972). *Evolution of the Protoplanetary Cloud and Formation of the Earth and Planets*. NASA TT F-677.
- Visser, R. G., & Ormel, C. W. (2016). On the growth of pebble-accreting planetesimals. *Astronomy & Astrophysics*, 586, A66. [https://doi.org/10.1051/0004-6361/201527361](https://doi.org/10.1051/0004-6361/201527361)

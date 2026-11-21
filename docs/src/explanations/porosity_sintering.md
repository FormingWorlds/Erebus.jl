# Porosity: cold pressing and hot sintering

This page explains the optional dry porosity model of `Erebus.jl`.
The model follows the laws of Henke et al. (2012) as Lichtenberg et al. (2016) used them.
Porosity $\phi$ is the volume fraction of pores in a rock, without units.
A young planetesimal forms from a loose dust powder with a high porosity.
Two processes lower this porosity: cold pressing and hot sintering.
The porosity in turn sets the density, the heat capacity and the thermal conductivity of the rock.

---

## What does `model = "none"` mean?

The short answer: with `model = "none"` the dry porosity model is off, and `Erebus.jl` runs as without it.

The section `[porosity]` of the configuration file selects the model.
The default `model = "none"` creates no extra marker arrays and draws no extra random numbers.
Rock porosity then evolves by the poroelastic compaction of the other sections, as before.
`model = "henke2012"` switches the dry porosity model on.
Each rock marker then carries 3 more values: the dry porosity $\phi_d$, its highest temperature, and the porosity it lost by sintering.

With `model = "henke2012"` the pores hold no fluid.
The water porosity of a rock marker stays at its floor `phimin`, and the compaction rate of the Stokes-Darcy solve does not change it.
The processes that need pore water are therefore not allowed together with this model.
The configuration check rejects `henke2012` when reaction, venting, hydrothermal convection, accretion or telescoping is on.

---

## How does pressure compact a cold powder?

The short answer: the weight of the overlying rock packs the grains closer, and the porosity drops to a value set by the pressure alone.

This process is cold pressing.
It does not depend on temperature.
Henke et al. (2012) fitted laboratory data of pressed powders with their Eq. 15.
Written for the porosity, the fit reads

$$\phi_H(P) = 0.42 + \frac{0.46}{(P / P_0)^{1.72} + 1}, \qquad P_0 = 0.13\ \text{bar}.$$

Here $P$ is the highest pressure the powder has felt, in bar.
Lichtenberg et al. (2016) use the same law as their Eq. 10.
The law gives $\phi_H = 0.88$ at $P = 0$ and $\phi_H = 0.65$ at $P = P_0$.
At $P = 10\,P_0 = 1.3$ bar it gives $\phi_H = 0.4286013$, close to the limit 0.42 of a dense random packing.

The porosity of a rock before sintering is the smallest of 3 values.
These are the initial porosity $\phi_{init}$ (`phi_init`), the pressing law $\phi_H(P)$, and an upper cap $\phi_{cap} = 0.85$ (`phi_cap`):

$$\phi_0 = \min\left(\phi_{init},\ \phi_H(P),\ \phi_{cap}\right).$$

The cap of 0.85 comes from the i2elvis code that Lichtenberg et al. (2016) ran.
Each time step first sets $\phi_d \leftarrow \min(\phi_d, \phi_0)$ at every temperature.
The porosity therefore never rises when the pressure drops.

---

## How does heat close the pores?

The short answer: above about 700 K the grains creep at their contacts, the contact faces grow, and the pores shrink.

This process is hot sintering.
`Erebus.jl` applies it at marker temperatures $T \ge T_{sinter}$, with $T_{sinter} = 700$ K (`T_sinter`), as Lichtenberg et al. (2016) did.

### The grain geometry

Henke et al. (2012) describe the rock as equal spheres of initial radius $R_0$ in a regular packing.
Each sphere touches 8 neighbours.
Sintering flattens each contact into a circular face of radius $a$, in m.
The grain keeps its volume, so its radius $R$, in m, grows above $R_0$.
With the ratio $X = (1 - \phi_0)/(1 - \phi)$ of relative densities, their Eqs. 34 to 37 give

$$R^3 - 3 R^2 R_0 X^{1/3} + R_0^3 X + R_0^3 = 0, \qquad a^2 = R^2 - R_0^2 X^{2/3}.$$

The cubic has its physical root between $R_0$ and $2 R_0 X^{1/3}$.
Before sintering $X = 1$, so $R = R_0$ and $a = 0$: the grains touch at points.
`Erebus.jl` finds this root by Newton iteration inside a bisection bracket.

### The stress on the contact faces

The applied pressure acts on the small contact faces only.
The stress on a face, $\sigma$ in bar, is therefore much larger than $P$ early in sintering (Henke et al., 2012, Eqs. 40 and 41):

$$\sigma = \frac{C_{av}\, P}{\pi a^2}, \qquad C_{av} = 2\sqrt{3}\,(R^2 - a^2).$$

Here $C_{av}$, in $\text{m}^2$, is the mean cross-section of the cell of one grain.

### The sintering rate

The porosity falls at a rate set by power-law creep of olivine.
Henke et al. (2012) took the constants of this law from laboratory work on olivine spheres (their Sect. 6.2).
`Erebus.jl` uses

$$\frac{d\phi}{dt} = -F(\phi)\, A\, (1 - \phi)\, \frac{\sigma^{n}}{(100\,R)^3}\, \exp\left(-\frac{E_{act}}{R_{gas} T}\right).$$

The symbols and their units are these.
The prefactor is $A = 4 \times 10^{-5}$ $\text{cm}^3\,\text{bar}^{-n}\,\text{s}^{-1}$ (`A`).
The stress exponent is $n = 1.5$ (`n`).
The activation energy is $E_{act} = 85000$ cal/mol (`E_act`).
The gas constant is $R_{gas} = 1.9872$ cal/(K mol).
The factor $100\,R$ converts the grain radius from m to cm.
The factor $F(\phi) = 10\,\phi$ for $\phi \le 0.1$ and $F = 1$ above (Henke et al., 2012, Eq. 43).
It slows the last stage, when single pores close.

The printed equations carry sign errors.
In Eq. 44 of Henke et al. (2012) the exponent reads $+E_{act}/(T R_{gas})$.
In Eq. 13 of Lichtenberg et al. (2016) the rate $\partial\phi/\partial t$ has no minus sign, so porosity would grow.
`Erebus.jl` uses the form above: a negative exponent, and porosity that only falls.
The i2elvis code uses the same signs.

### The grain radius $R_0$

`Erebus.jl` uses $R_0 = 10^{-6}$ m (`grain_radius`), the value in the i2elvis code.
Henke et al. (2012) estimate matrix particle radii of about 0.1 to 0.2 µm for chondrites.
The sintering rate scales as $R^{-3}$ through the creep law, so $R_0$ matters.

---

## How does `Erebus.jl` step the sintering law in time?

The short answer: with a backward Euler step, split into substeps that each lower $\phi$ by at most 0.01.

The rate is unbounded at the start of sintering, because $a = 0$ there.
An explicit step would then remove all porosity at once.
`Erebus.jl` solves instead, for each substep of length $h$ in s,

$$\phi^{k+1} = \phi^{k} - h\, r(\phi^{k+1}),$$

where $r \ge 0$ is the rate magnitude from the law above.
A bisection finds $\phi^{k+1}$ in $\left(0, \phi^k\right]$.
The substep length adapts so that each substep lowers $\phi$ by at most 0.01.
At most 1000 substeps run per marker and time step; the last one then takes the remaining time, and the run reports the count.
The step returns $\phi^n$ unchanged when $P = 0$ or $\phi^n = 0$.

Near the start the porosity drop grows as a power of time.
For a small drop $\Delta\phi = \phi_0 - \phi$, the law gives

$$\Delta\phi(t) = C\, t^{2/5}, \qquad C = \left[(n+1)\,K\right]^{1/(n+1)}.$$

Here $K = F(\phi_0)\, A\, (1-\phi_0)\, \left(3\sqrt{3}\,(1-\phi_0)\,P/\pi\right)^{n} (100 R_0)^{-3} \exp\left(-E_{act}/(R_{gas} T)\right)$, with $P$ in bar and $R_0$ in m.

!!! details "Derivation of the onset law"
    Let $\epsilon = 1 - X \approx \Delta\phi / (1 - \phi_0)$ be small.
    Write $R = R_0 (1 + \delta)$ and expand the cubic to second order.
    The terms give $-3\delta + \epsilon^2/3 = 0$, so $\delta = \epsilon^2 / 9$.
    Then $a^2 = R^2 - R_0^2 X^{2/3} \approx \tfrac{2}{3} R_0^2 \epsilon$ and $C_{av} \approx 2\sqrt{3} R_0^2$.
    The stress becomes $\sigma \approx 3\sqrt{3}\,(1-\phi_0)\,P / (\pi\, \Delta\phi)$.
    The rate is then $d\Delta\phi/dt = K\, \Delta\phi^{-n}$, with $K$ as above.
    Integration from $\Delta\phi = 0$ at $t = 0$ gives $\Delta\phi^{n+1} = (n+1) K t$.
    For $n = 1.5$ the exponent of time is $1/(n+1) = 2/5$.
    One backward Euler step of length $t$ from $\phi_0$ solves $\Delta\phi^{n+1} = K t$ instead.
    It therefore gives $\left(2/5\right)^{2/5} \approx 0.69$ of the exact drop, for any step length.
    This is why the step is split into substeps.

---

## How does porosity change the heat flow?

The short answer: pores lower the density, the heat capacity, the heat production and, most strongly, the thermal conductivity.

The pores of the dry rock are empty.
For a rock marker with dry porosity $\phi_d$ and solid values $\rho_s$, $\left(\rho c_p\right)_s$ and $H_s$, the bulk values are

$$\rho = \rho_s (1 - \phi_d), \qquad \rho c_p = (\rho c_p)_s (1 - \phi_d), \qquad H = H_s (1 - \phi_d).$$

Here $\rho$ is the density in $\text{kg/m}^3$, $\rho c_p$ the volumetric heat capacity in $\text{J}/(\text{m}^3\,\text{K})$, and $H$ the radiogenic heating in $\text{W/m}^3$.
Lichtenberg et al. (2016) use the same density law as their Eq. 11.

The thermal conductivity $k$, in W/(m K), follows the 3-branch form of Lichtenberg et al. (2016, Eq. 12):

$$k_{por} = \begin{cases} k_1 & \phi < 0.2 \\ \left(k_1^4 + k_2^4\right)^{1/4} & 0.2 \le \phi \le 0.4 \\ k_2 & \phi > 0.4 \end{cases} \qquad k_1 = k\, e^{-\phi/0.08}, \quad k_2 = k\, e^{-1.2 - \phi/0.167}.$$

Here $k$ is the conductivity of the dense solid.
A floor of $10^{-3}$ W/(m K) applies.
Henke et al. (2012) fit $k_1$ and $k_2$ to meteorite and powder data (their Eqs. 28 and 29).
Their Eq. 30 joins the two fits into 1 smooth function for all porosities.
The 3-branch form instead jumps at $\phi = 0.2$, by about 26 %, and very slightly at $\phi = 0.4$.
The sign of the constant $-1.2$ follows Lichtenberg et al. (2016) and the i2elvis code.

These values replace the pore-fluid mixture of the other models at every place that sets rock properties.
These places are the marker-to-grid mapping, the sub-grid temperature diffusion, the radiogenic heating, and the radiative surface boundary.

---

## Where does the pressure of a marker come from?

The short answer: from a hydrostatic profile at the start, and from the total pressure of the Stokes solve after that.

At the start the body has no Stokes pressure yet.
`Erebus.jl` integrates the hydrostatic balance $dP/dr = -\rho(r)\, g(r)$ inward from $P = 0$ at the surface.
The density $\rho = \rho_s (1 - \phi_0(P))$ depends on the pressure, so density and pressure are iterated to a fixed point.
The gravity $g(r)$ is that of the enclosed mass, for a sphere or for a cylinder, as `gravity_mode` sets.
The first cold pressing then uses this profile.

In each time step the marker pressure is the total pressure of the Stokes solve, interpolated to the marker.
The surface anchor `psurface` is subtracted, negative values are set to 0, and Pa are converted to bar.
The fluid pressure is not used.
The step runs where the marker temperature becomes the new temperature of the step.

---

## Known limits

The materials of `Erebus.jl` are fixed constants.
The dense rock density is therefore 3300 $\text{kg/m}^3$ for dry rock (`XWsolidm_init = 0`), not the 3500 $\text{kg/m}^3$ of Lichtenberg et al. (2016).
The temperature and pressure terms of the density in the i2elvis code are not included.

For $X \le 1/3$ the contact faces of the geometry overlap and the cubic has no root.
This happens only for $\phi_0 > 2/3$ once $\phi \le 3\phi_0 - 2$.
`Erebus.jl` then holds the geometry at $X = 1/3 + 10^{-12}$.
The rate still falls to 0 at $\phi = 0$ through the factor $F(\phi)$.

Molten rock keeps its dry porosity; the model does not reset it on melting.

The step needs the Stokes pressure.
A run that skips the Stokes solve keeps the porosity of the first cold pressing.

---

## Which tests check this?

No validation case runs this model yet.
The unit tests in `test/test_porosity_henke.jl` check each law and each place in the code that uses it.

!!! details "Test details"
    Testsets of `Porosity laws of Henke et al. (2012)`:
    the cold pressing law at 0, 0.13 and 1.3 bar;
    the rule for $\phi_0$ and its cap;
    the grain root against a BigFloat bisection and volume conservation to 1e-12;
    the contact stress by hand;
    the explicit rate against a BigFloat reference coded from the raw constants;
    order 1 convergence of the sub-stepped step to a fine BigFloat Runge-Kutta integration;
    the onset law $C t^{2/5}$ against a BigFloat quadrature;
    the 3 branches and the floor of the conductivity.
    The testset `Replay of the i2elvis porosity block` compares 1 explicit step with the i2elvis code.
    It uses 8 marker states and a relative tolerance of 1e-6.
    `tools/replay_i2elvis_porosity.jl` compiles that block from a local i2elvis clone and prints the recorded table.
    Testsets of `Dry porosity marker group and step` check the marker arrays and the setup pressure against the uniform-density limit.
    They also check the step in `vent_and_degas!`, the rock properties at every call site, and the marker lifecycle.
    A 2-step run checks that the Stokes pressure sinters hot rock.

    Command: `julia --project -e 'using Test; include("test/test_porosity_henke.jl")'`

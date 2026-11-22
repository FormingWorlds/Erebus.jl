# Porosity: cold pressing and hot sintering

This page explains the optional dry porosity model of `Erebus.jl`.
The model follows the laws of Henke et al. (2012) as Lichtenberg et al. (2016a) used them.
Porosity $\phi$ is the volume fraction of pores in a rock, without units.
A young planetesimal forms from a loose dust powder with a high porosity.
Two processes lower this porosity: cold pressing and hot sintering.
The porosity in turn sets the density, the heat capacity and the thermal conductivity of the rock.

---

## What does `model = "none"` mean?

The short answer: with `model = "none"` the dry porosity model is off, and `Erebus.jl` runs as without it.

The section `[porosity]` of the configuration file selects the model.
The default `model = "none"` creates no extra marker arrays and draws no extra random numbers.
Rock porosity then evolves by the poroelastic compaction of the other sections.
`model = "henke2012"` switches the dry porosity model on.
Each rock marker then carries 4 more values.
These are the dry porosity $\phi_d$, the porosity before sintering $\phi_0$, the highest temperature of the marker in K, and the porosity it lost by sintering.

With `model = "henke2012"` the pores hold no fluid.
The water porosity of a rock marker stays at its floor `phimin`.
`phimin` is the smallest water porosity of the section `[poroelasticity]`, $10^{-4}$ by default.
The compaction rate of the Stokes-Darcy solve does not change the water porosity.
It therefore does not limit the time step either.
The Darcy velocity of this pore water moves no marker.
The time step therefore takes the solid velocity and the solid displacement limit `time.dxymax` instead.
The heat that pore water carries through the grid is 0.
The shear heating has no compaction and Darcy terms, and the adiabatic heating has no pore-water term.
The fluid pressure moves with the rock from one step to the next.
It does not decide when the outer iteration of a time step has converged.
The configuration check rejects a `phimin` of 0 or less, or above $10^{-3}$, with this model.
The adiabatic heating and the compaction term of the Stokes-Darcy solve still use the water porosity of rock.
There the bound keeps the solid fraction $1 - \phi$ within 0.1 % of 1.
The configuration check also bounds the permeability of this water, as the next section explains.
The columns `phi_max` and `phi_mean` of `telemetry.csv` give the maximum and the mean of $\phi_d$ over the rock markers, not the water porosity.
The processes that need pore water, melt or metal are therefore not allowed together with this model.
The configuration check rejects `henke2012` when reaction, venting, hydrothermal convection, accretion or telescoping is on.
It also rejects volatiles, pyrolysis kinetics, magma transport, magma degassing, core formation, iron heating (`hr_fe`) and metal partitioning.
It rejects the tracking of accessory metal phases (`phase_tracking`) and the volatile ice mixture (`volatile_mixture`) too.
It rejects the redox model (`redox`) as well, because it tracks metallic iron and gas species.
It rejects hydrofracture (`poroelasticity.hydrofracture`), because no pore fluid pressure can open a crack.
The rock must be dry: `materials.XWsolidm_init = [0.0, 0.0, nan]`.

---

## Why does the configuration bound the permeability of rock?

The short answer: the Stokes-Darcy solve still lets the water porosity of rock drain, and a small permeability keeps this drainage from moving the rock.

Every pore-water term that heats rock, moves a marker or limits the time step is off with `model = "henke2012"`.
One coupling stays: the compaction term of the solid momentum and continuity equations.
Through it the solid compacts where the water of the porosity `phimin` drains.
The water drains faster for a larger Darcy mobility $M$, in $\text{m}^2\,\text{Pa}^{-1}\,\text{s}^{-1}$.
$M$ is the permeability $k$ at the water porosity $\phi_{min}$, divided by the smallest water viscosity $\eta_f$ of rock:

$$M = \frac{k(\phi_{min})}{\eta_f}, \qquad k(\phi) = k_0 \left(\frac{\phi}{\phi_{ref}}\right)^3 \left(\frac{1 - \phi}{1 - \phi_{ref}}\right)^{-2}.$$

Here $\phi_{min}$ is `poroelasticity.phimin`, $k_0$ is `materials.kphim0` of mantle and crust, in $\text{m}^2$, and $\phi_{ref}$ is `thermodynamics.phim0`.
$\eta_f$ is in Pa s.
The Arrhenius law of the water viscosity has a floor of $10^{-5}$ Pa s, so this floor is $\eta_f$ by default.
The default configuration has $k(\phi_{min}) = 8.0 \times 10^{-24}$ $\text{m}^2$ and $M$ = 8.0e-19 $\text{m}^2\,\text{Pa}^{-1}\,\text{s}^{-1}$.
With `phimin` = 1e-3, the largest allowed value, it has $M$ = 8.0e-16 $\text{m}^2\,\text{Pa}^{-1}\,\text{s}^{-1}$.

We measured the effect of $M$ with runs of 4 time steps at `phimin` = 1e-3.
Each pair of runs differs only in $M$, by a factor of $10^6$.
The temperature and $\phi_d$ change little.
Over the whole range of $M$ in the 200 km and 50 km bodies, the temperature changes by at most $5.1 \times 10^{-5}$ of its own change over the 4 steps.
$\phi_d$ changes by at most $3.4 \times 10^{-4}$ of its change.
The solid velocity changes most.
At small $M$ its change grows in proportion to $M$.
At large $M$ the water drains completely, and more $M$ changes nothing.
A 200 km body at 900 K then moves 1.2 % faster than at $M = 0$.

The configuration check therefore rejects $M$ above 1e-15 $\text{m}^2\,\text{Pa}^{-1}\,\text{s}^{-1}$ with `henke2012`.
At $\eta_f = 10^{-5}$ Pa s this bound is $k(\phi_{min}) \le 10^{-20}$ $\text{m}^2$.
The table gives the change between $k(\phi_{min}) = 10^{-20}$ $\text{m}^2$ and $10^{-26}$ $\text{m}^2$.
The velocity column is relative to the largest solid speed.
The temperature and $\phi_d$ columns are relative to their largest change over the 4 steps.
The water viscosity of these runs is $1.6 \times 10^{-5}$ Pa s at 900 K and $1.1 \times 10^{-5}$ Pa s at 1100 K.

| Body radius | Temperature | Solid velocity | Temperature change | $\phi_d$ change |
|:------------|:------------|:---------------|:-------------------|:-----------------|
| 200 km | 900 K | $1.6 \times 10^{-5}$ | $2.4 \times 10^{-8}$ | $1.0 \times 10^{-6}$ |
| 50 km | 900 K | $1.4 \times 10^{-4}$ | $1.6 \times 10^{-8}$ | $2.6 \times 10^{-6}$ |
| 20 km | 900 K | $4.4 \times 10^{-4}$ | $1.9 \times 10^{-6}$ | $2.4 \times 10^{-4}$ |
| 200 km | 1100 K mantle, 300 K crust | $1.7 \times 10^{-5}$ | $7.9 \times 10^{-9}$ | $6.8 \times 10^{-7}$ |
| 20 km | 1100 K | $3.7 \times 10^{-4}$ | $8.3 \times 10^{-7}$ | $3.5 \times 10^{-5}$ |
| 200 km, yielding, `:enclosed_mass` | 1100 K mantle, 300 K crust | $5.4 \times 10^{-5}$ | $3.8 \times 10^{-9}$ | $3.5 \times 10^{-6}$ |
| 200 km, yielding, `:poisson2d` | 1100 K mantle, 300 K crust | $3.6 \times 10^{-5}$ | $3.1 \times 10^{-9}$ | $1.7 \times 10^{-6}$ |

The 2 yielding runs have a friction coefficient of 0.005 and a cohesion of $10^4$ Pa, and 38 and 30 basic nodes yield in both runs of each pair.

The 20 km body at 900 K reaches a velocity change of $10^{-3}$ at $k(\phi_{min}) = 3 \times 10^{-20}$ $\text{m}^2$, where $M = 1.9 \times 10^{-15}$ $\text{m}^2\,\text{Pa}^{-1}\,\text{s}^{-1}$.
Smaller bodies respond more strongly, so the bound is measured for radii of 20 km and more.

---

## How does pressure compact a cold powder?

The short answer: the weight of the overlying rock packs the grains closer, and the porosity drops to a value set by the pressure alone.

This process is cold pressing.
It does not depend on temperature.
Güttler et al. (2009, in Henke et al., 2012) pressed fine dust powders in the laboratory and fitted the result.
Henke et al. (2012) use this fit as their Eq. 15.
Written for the porosity, the fit reads

$$\phi_H(P) = 0.42 + \frac{0.46}{(P / P_0)^{1.72} + 1}, \qquad P_0 = 0.13\ \text{bar}.$$

Here $P$ is the highest pressure the powder has felt, in bar.
The model computes pressures in Pa and multiplies them by $10^{-5}$ bar/Pa for this law.
Lichtenberg et al. (2016a) use the same law as their Eq. 10.
The law gives $\phi_H = 0.88$ at $P = 0$ and $\phi_H = 0.65$ at $P = P_0$.
At $P = 10\,P_0 = 1.3$ bar it gives $\phi_H = 0.4286012$, close to the limit 0.42 of a dense random packing.

The porosity of a rock before sintering is the smallest of 3 values.
These are the initial porosity $\phi_{init}$ (`phi_init`), the pressing law $\phi_H(P)$, and an upper cap $\phi_{cap} = 0.85$ (`phi_cap`):

$$\phi_0 = \min\left(\phi_{init},\ \phi_H(P),\ \phi_{cap}\right).$$

The cap of 0.85 comes from the i2elvis code that Lichtenberg et al. (2016a) ran.
Each marker stores its $\phi_0$.
Until the marker starts to sinter, each time step replaces it by the smaller of the stored value and $\phi_0$ at the present pressure.
The stored $\phi_0$ therefore belongs to the highest pressure the marker has felt, because $\phi_H$ falls with pressure.
The step then sets $\phi_d \leftarrow \min(\phi_d, \phi_0)$ at every temperature.
The porosity therefore never rises when the pressure drops.

Once a marker has sintered, its $\phi_0$ stays fixed and cold pressing stops for it.
Henke et al. (2012, Eq. 37) define $\phi_0$ through the relative density of the initial packing of the grains.
In that packing the grains touch at points, before any sintering.
A later rise of pressure therefore does not change $\phi_0$ and does not shrink the contact faces of sintered grains.
The i2elvis code differs here: it computes $\phi_0$ from the current pressure in every time step.

A marker counts as sintered once its cumulative porosity loss by sintering exceeds $10^{-12}$ (`DPHI_SINTERED`).
A smaller loss is numerical, not sintering.
Near the onset of sintering, the rounding of $1 - \phi$ in floating point lowers $\phi$ by a few $10^{-16}$.
A marker with such a small loss is still cold pressed when the pressure rises.

For a marker that starts to sinter at $\phi_d = \phi_0$, the porosity lost by sintering equals $\phi_0 - \phi_d$ to within $10^{-12}$.
Every marker of the setup starts this way.
A marker added later can lie below its $\phi_0$ when it starts to sinter.
For such a marker, $\phi_0 - \phi_d$ exceeds the porosity lost by sintering by the gap $\phi_0 - \phi_d$ at the onset of sintering, after any cold pressing.

---

## How does heat close the pores?

The short answer: above about 700 K the grains creep at their contacts, the contact faces grow, and the pores shrink.

This process is hot sintering.
`Erebus.jl` applies it at marker temperatures $T \ge T_{sinter}$, with $T_{sinter} = 700$ K (`T_sinter`), as Lichtenberg et al. (2016a) did.

### The grain geometry

Henke et al. (2012) describe the rock as equal spheres of initial radius $R_0$, in m, in a regular packing.
Each sphere touches 8 neighbours.
Sintering flattens each contact into a circular face of radius $a$, in m.
The grain keeps its volume, so its radius $R$, in m, grows above $R_0$.
With the ratio $X = (1 - \phi_0)/(1 - \phi)$ of relative densities, their Eqs. 34 to 37 give

$$R^3 - 3 R^2 R_0 X^{1/3} + R_0^3 X + R_0^3 = 0, \qquad a^2 = R^2 - R_0^2 X^{2/3}.$$

The relation for $a$ is the geometry of Kakar and Chaklader (1967, in Henke et al., 2012).
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
Henke et al. (2012) took the constants of this law from the experiments on small olivine spheres of Schwenn and Goetze (1978, in Henke et al., 2012).
`Erebus.jl` uses

$$\frac{d\phi}{dt} = -F(\phi)\, A\, (1 - \phi)\, \frac{\sigma^{n}}{(100\,R)^3}\, \exp\left(-\frac{E_{act}}{R_{gas} T}\right).$$

The symbols and their units are these.
The time is $t$, in s.
The prefactor is $A = 4 \times 10^{-5}$ $\text{cm}^3\,\text{bar}^{-n}\,\text{s}^{-1}$ (`A`).
The stress exponent is $n = 1.5$ (`n`).
The activation energy is $E_{act} = 85000$ cal/mol (`E_act`).
The gas constant is $R_{gas} = 1.9872$ cal/(K mol).
The temperature $T$ is in K.
The factor $100\,R$ converts the grain radius from m to cm.
The factor $F(\phi) = 10\,\phi$ for $\phi \le 0.1$ and $F = 1$ above (Henke et al., 2012, Eq. 43).
It slows the last stage, when single pores close.

The printed equations carry sign errors.
In Eq. 44 of Henke et al. (2012) the exponent reads $+E_{act}/(T R_{gas})$.
In Eq. 13 of Lichtenberg et al. (2016a) the rate $\partial\phi/\partial t$ has no minus sign, so porosity would grow.
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
Let $k$ count the substeps within 1 time step of the model.
`Erebus.jl` solves, for each substep of length $h$ in s,

$$\phi^{k+1} = \phi^{k} - h\, r(\phi^{k+1}),$$

where $r \ge 0$ is the rate magnitude from the law above, in 1/s.
A bisection finds $\phi^{k+1}$ in $\left(0, \phi^k\right]$.
It halves the bracket until its 2 ends are adjacent floating-point numbers.
The substep length adapts so that each substep lowers $\phi$ by at most 0.01.
This limit (`dphi_max`) must be at least $10^{-12}$.
A smaller limit could not step over the rounding of $1 - \phi$ near $\phi_0$.
At most 1000 substeps run per marker and time step.
The last one then takes the remaining time, and the run reports the number of markers that reached this cap.
Let $\phi_{start}$ be the porosity at the start of the time step.
The step returns $\phi_{start}$ unchanged when $P = 0$, $\phi_{start} = 0$ or the time step is 0.

Near the start the porosity drop grows as a power of time.
Let $t$ be the time since sintering started, in s.
For a small drop $\Delta\phi = \phi_0 - \phi$, the law gives

$$\Delta\phi(t) = C\, t^{1/(n+1)}, \qquad C = \left[(n+1)\,K\right]^{1/(n+1)}.$$

For $n = 1.5$ the drop grows as $t^{2/5}$.
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
    It therefore gives $\left(1/(n+1)\right)^{1/(n+1)} = \left(2/5\right)^{2/5} \approx 0.69$ of the exact drop, for any step length.
    This is why the step is split into substeps.

---

## How does porosity change the heat flow?

The short answer: pores lower the density, the heat capacity, the heat production and, most strongly, the thermal conductivity.

The pores of the dry rock are empty.
For a rock marker with dry porosity $\phi_d$ and solid values $\rho_s$, $\left(\rho c_p\right)_s$ and $H_s$, the bulk values are

$$\rho = \rho_s (1 - \phi_d), \qquad \rho c_p = (\rho c_p)_s (1 - \phi_d), \qquad H = H_s (1 - \phi_d).$$

Here $\rho$ is the density in $\text{kg/m}^3$, $\rho c_p$ the volumetric heat capacity in $\text{J}/(\text{m}^3\,\text{K})$, and $H$ the radiogenic heating in $\text{W/m}^3$.
Lichtenberg et al. (2016a) use the same density law as their Eq. 11.

The thermal conductivity of the porous rock, $k_{por}$ in W/(m K), follows the 3-branch form of Lichtenberg et al. (2016a, Eq. 12):

$$k_{por} = \begin{cases} k_1 & \phi < 0.2 \\ \left(k_1^4 + k_2^4\right)^{1/4} & 0.2 \le \phi \le 0.4 \\ k_2 & \phi > 0.4 \end{cases} \qquad k_1 = k\, e^{-\phi/0.08}, \quad k_2 = k\, e^{-1.2 - \phi/0.167}.$$

Here $k$ is the conductivity of the dense solid, also in W/(m K).
A floor of $10^{-3}$ W/(m K) applies.
Henke et al. (2012) fit $k_1$ and $k_2$ to meteorite and powder data (their Eqs. 28 and 29).
Their Eq. 30 joins the two fits into 1 smooth function for all porosities.
The 3-branch form instead jumps at $\phi = 0.2$, by about 26 %, and very slightly at $\phi = 0.4$.
The sign of the constant $-1.2$ follows Lichtenberg et al. (2016a) and the i2elvis code.

These values replace the pore-fluid mixture of the other models for the density, the heat capacity, the conductivity and the radiogenic heating of rock.
The code sets these 4 values in the marker-to-grid mapping, the sub-grid temperature diffusion and the radiogenic heating.
The radiative surface boundary uses $k_{por}$ too, at the porosity $\phi_0(P = 0)$ of unpressed powder.

---

## Where does the pressure of a marker come from?

The short answer: from a hydrostatic profile at the start, and from the total pressure of the Stokes solve after that.

At the start the body has no Stokes pressure yet.
Let $r$ be the distance from the centre of the body, in m, and $R_p$ the radius of the body, in m.
`Erebus.jl` integrates the hydrostatic balance $dP/dr = -\rho(r)\, g(r)$ inward from $P = 0$ at $r = R_p$.
Here $P$ is in Pa and the gravity $g(r)$ is in $\text{m/s}^2$.
The density $\rho = \rho_s (1 - \phi_0(P))$ depends on the pressure, so density and pressure are iterated to a fixed point.
The porosity law takes the pressure in bar, $10^{-5}$ times the pressure in Pa.

The gravity is that of the gravity mode of the Stokes solve, set by `gravity_mode`.
For `:enclosed_mass` it is the gravity of a sphere, $g = G M(r)/r^2$.
Here $M(r)$ is the mass inside radius $r$, in kg.
$G = 6.674 \times 10^{-11}$ $\text{m}^3\,\text{kg}^{-1}\,\text{s}^{-2}$ is the gravitational constant.
For `:poisson2d` the Stokes solve takes $g$ from the 2-D Poisson equation $\nabla^2 \Phi = (8/3)\pi G \rho$.
Here $\Phi$ is the gravitational potential, in J/kg.
The profile uses the same source term, so $g = (8/3)\pi G\, m(r)/r$.
Here $m(r) = \int_0^r \rho\, r'\, dr'$, in kg/m, integrates over the radius $r'$, in m, from the centre to $r$.
The mass per unit length of a cylinder of radius $r$ is $2\pi\, m(r)$.
A body of uniform density then has the central pressure $(2/3)\pi G \rho^2 R_p^2$ in both modes.
At the first step the centre pressure of the profile and of the Stokes solve agree to within 5 %.
The first cold pressing uses this profile.

In each time step the marker pressure is the total pressure of the Stokes solve, interpolated to the marker.
The surface anchor `psurface` is subtracted, negative values are set to 0, and Pa are converted to bar.
The fluid pressure is not used.
The melting law also takes the total pressure of the Stokes solve at the marker, from the second step on.
The step runs where the marker temperature becomes the new temperature of the step.

The plastic yield stress of the rock needs a pore pressure as well.
Its friction term grows with the effective pressure, the total pressure minus the pore pressure.
The Stokes solve fixes both pressures to `psurface` at the edge of the domain, so they are absolute pressures.
The dry pores hold only the gas at the surface, so their pressure is `psurface`.
The effective pressure of the yield stress is therefore the total pressure minus `psurface`.
The porosity step takes the same pressure.
The fluid pressure of the water porosity `phimin` comes close to the total pressure.
With it, the effective pressure would almost vanish, and rock of low cohesion would yield at a small shear stress.

---

## What does the model leave out?

The short answer: it leaves out some density terms, the shrinking of the body, the dry porosity in the adiabatic heating and the Stokes-Darcy solve, the compaction work, a reset on melting and the sintered surface porosity.
Its grid must also resolve the surface layer.

The materials of `Erebus.jl` are fixed constants.
The dense rock density is therefore 3300 $\text{kg/m}^3$ for dry rock (`XWsolidm_init = 0`), not the 3500 $\text{kg/m}^3$ of Lichtenberg et al. (2016a).
The temperature and pressure terms of the density in the i2elvis code are not included.

The body does not shrink when its pores close.
Sintering lowers $\phi_d$, while the grid and the radius $R_p$ of the body stay fixed.
The rock density $\rho_s (1 - \phi_d)$ then rises, where $\rho_s$ is the dense rock density in $\text{kg/m}^3$.
The rock mass of the body therefore grows.
The i2elvis code does the same: its density function also multiplies the dense rock density by $1 - \phi$.
A 6-step run of a 50 km body at 900 K shows the size of the gain.
With `:enclosed_mass` gravity the rock mass of the sphere rises by 37 % from the setup.
With `:poisson2d` gravity the rock mass per unit length of the 2-D slice rises by 46 %.
Henke et al. (2012, Sect. 7.1) instead keep the mass of each shell fixed and compute its radius from the porosity.
Their body shrinks: in their porous model the radius falls by about 20 % (Henke et al., 2012, Sect. 7.2.2).
In `Erebus.jl` the gravity and the pressure rise with the mass, not with the shrinking radius.
For a body of uniform density, let the density rise by a factor $f$, without units.
At fixed radius the centre pressure rises by $f^2$.
For a sphere of fixed mass, as with `:enclosed_mass` gravity, it rises by $f^{4/3}$, because the radius shrinks as $f^{-1/3}$.
For the 2-D slice of `:poisson2d` gravity the mass per unit length stays fixed instead.
Its radius shrinks as $f^{-1/2}$, and the centre pressure rises by $f$.
The radiogenic heating and the heat capacity per volume both scale with $1 - \phi_d$.
The heating rate of the rock in K/s therefore does not change.

The adiabatic heating does not use the dry porosity $\phi_d$.
It weights the solid term with $1 - \phi$, where $\phi$ is the water porosity of the grid.
In this model that porosity stays at `phimin`, so the term counts the rock as dense.
At the centre of the test body of 50 km, at the start, $1 - \phi_d = 0.58$.
The adiabatic heating there is therefore 1.72 times the value with $1 - \phi_d$.

The Stokes-Darcy solve does not use the dry porosity either.
It keeps the water porosity `phimin` in rock, with its bulk viscosity and its Darcy flow.
Their compaction and Darcy terms do not enter the shear heating, and the pore-water term does not enter the adiabatic heating.
The heat that this pore water would carry is set to 0, and its compaction rate does not limit the time step.
In a 4-step test run of a 200 km body at 900 K the whole shear heating stays below 1e-6 $\text{W/m}^3$.
Only the compaction term of the solid equations still feels this water, within the bound on its mobility.

The energy equation has no term for the work of compaction.
Henke et al. (2012, Eqs. 20 to 23) add 2 heat sources that compaction causes.
These are the compressional work and the gravitational energy that a shrinking body releases.
The heat sources of `Erebus.jl` are radiogenic, adiabatic, shear and latent heating, and none of them contains these 2 terms.

For $X \le 1/3$ the contact faces of the geometry overlap and the cubic has no root.
This happens only for $\phi_0 \ge 2/3$ once $\phi \le 3\phi_0 - 2$.
`Erebus.jl` then holds the geometry at $X = 1/3 + 10^{-12}$.
The rate still falls to 0 at $\phi = 0$ through the factor $F(\phi)$.

Molten rock keeps its dry porosity.
The model does not reset it on melting.

The step needs the Stokes pressure.
A run that skips the Stokes solve keeps the porosity of the first cold pressing.

The radiative surface boundary takes the conductivity of unpressed powder.
It uses $k_{por}$ at $\phi_0(P = 0) = \min(\phi_{init}, 0.88, \phi_{cap})$, not the porosity of the rock markers at the surface.
A sintered or cold pressed surface layer therefore has a lower boundary conductivity than its markers.

The Stokes pressure lives on the grid, with cells of width $\Delta x$, in m.
Near the surface the step interpolates the pressure between nodes inside the body and nodes in the sticky air.
The top marker layer then feels a pressure of the order of the pressure at depth $\Delta x$.
If $\Delta x$ exceeds the depth where $P$ reaches $P_0 = 0.13$ bar, the first step cold presses this layer below the profile value.
In Pa, $P_0 = 1.3 \times 10^4$ Pa.
A grid resolves the uncompacted surface layer only if

$$\Delta x < \frac{P_0}{\rho\, g},$$

with $P_0$ in Pa, $\rho$ the density in $\text{kg/m}^3$ and $g$ the gravity at the surface in $\text{m/s}^2$.
For $\phi_{init} = 0.75$ the surface density is $3300 \times 0.25 = 825$ $\text{kg/m}^3$.
With `:enclosed_mass` gravity, $P$ reaches $P_0$ at a depth of 1667 m for $R_p = 20$ km and 564 m for $R_p = 50$ km.
The test grids of 17 by 17 nodes have $\Delta x = 8.75$ km for $R_p = 50$ km and do not resolve this layer.

---

## Which tests check this?

No validation case runs this model yet.
The unit tests in `test/test_porosity_henke.jl` check each law and each place in the code that uses it.

!!! details "Test details"
    Testsets of `Porosity laws of Henke et al. (2012)`:
    the cold pressing law at 0, 0.13 and 1.3 bar;
    the rule for $\phi_0$ and its cap;
    the grain root against a BigFloat bisection and volume conservation to 1e-12;
    the contact stress pinned by hand, with prefactor and sign guards;
    the explicit rate against a BigFloat reference coded from the raw constants;
    order 1 convergence of the sub-stepped step to a fine BigFloat Runge-Kutta integration;
    the onset law $C t^{2/5}$ against a BigFloat quadrature;
    the 3 branches and the floor of the conductivity, pinned by hand;
    a rate that overflows throws within 5 s instead of shrinking the substep forever;
    a substep limit below 1e-12 throws;
    1 substep at a rate far above $1/h$ ends at the smallest positive float.
    The testset `Replay of the i2elvis porosity block` compares 1 explicit step with the i2elvis code.
    It uses 8 marker states and a relative tolerance of 1e-6.
    `tools/replay_i2elvis_porosity.jl` compiles that block from a local i2elvis clone and prints the recorded table.
    Testsets of `Configuration section [porosity]` check the processes that `henke2012` excludes, also at `init_simulation`.
    They check that a `phimin` of 0 or less or above 1e-3, and hydrofracture, are rejected with `henke2012`.
    A Darcy mobility of rock above 1e-15 $\text{m}^2\,\text{Pa}^{-1}\,\text{s}^{-1}$ is rejected, also from the crust permeability, a constant water viscosity or the viscosity below the freezing point; the default 8.0e-19 and 8.0e-16 at `phimin` = 1e-3 pass.
    They check that the marker step and `dry_porosity_update` reject a field of a directly built `PorosityConfig` outside its domain.
    The marker step also rejects a bad `dt` or `psurface` before it changes a marker, also when no marker is rock.
    It checks `dt` and `psurface` in Float64, so a BigFloat of 1e400 is rejected and a BigFloat of 1e5 gives the Float64 result.
    The yield stress, `dry_porosity_update` and the cold pressing law check their pressure the same way.
    Testsets of `Dry porosity marker group and step` check the marker arrays and the setup pressure against the uniform-density limit.
    They check the residual of the hydrostatic balance and the setup pressure against the Stokes pressure of the first step, to 5 %.
    They check the setup pressure `pfm0` of the melting law and the 2-D Poisson profile of compressible powder.
    A 3-step run with melting checks that `pfm0` follows the total pressure, not the fluid pressure, to 1e-6.
    A direct call of `solve_stokes_darcy!` at step 2 checks the same, to 1e-6.
    They check the memory of $\phi_0$ and its freeze at the onset of sintering.
    A decrease below 1e-12 does not count as sintering, so a later cold pressing still lowers $\phi_0$.
    A 10-step run checks that the sintering decrease equals $\phi_0 - \phi_d$ to 1e-12 on every rock marker.
    They check the step in `vent_and_degas!` over 2 steps, the substep cap, the rock properties at every call site, and the marker lifecycle.
    A 2-step run checks that the Stokes pressure sinters hot rock.
    The marker step uses the `T_sinter`, `grain_radius`, `E_act` and `n` of the configuration.
    The setup profile and the radiative surface faces use `phi_cap` where it lies below `phi_init`.
    A 4-step run of a 200 km body at 900 K keeps the time step within a factor 3 of model none and the shear heating below 1e-6 $\text{W/m}^3$.
    The same run checks that the pore-water heat is 0 and that the telemetry porosity columns give $\phi_d$.
    The step stays within a factor 3 also with a water porosity `phim0` of 1e-4 and a permeability `kphim0` of 1e-12 $\text{m}^2$.
    This and the next runs lie above the mobility bound; a test helper sets them past the configuration check and repeats the step functions of `simulation_loop`, bit for bit at a valid configuration.
    With `phimin` = 1e-3 and `kphim0` of 1e-9 or 1e-12 $\text{m}^2$, a 4-step run keeps the time step, the shear and adiabatic heating stay below 1e-6 $\text{W/m}^3$, and each step takes 1 outer iteration.
    No cell fractures and no fracture flag flips where the fictitious fluid pressure exceeds the total pressure plus the tensile strength.
    Without pore water the fluid pressure is carried with the solid velocity, and `time.dxymax_fluid` does not change the step of a melting run.
    Fluid and solid pressures 1 MPa per column off the solution leave the heating, the temperature and the convergence of a dry thermal step unchanged; with pore water they change them.
    Active yielding gives the same yielding nodes and step at `psurface` of 1e3 Pa and 1e8 Pa, in both gravity modes.
    The telemetry porosity counts crust markers and is 0 without rock markers.
    The plastic yield stress takes the pore pressure `psurface`: a 3-step run of a 200 km body with a cohesion of 1e4 Pa keeps the step of the default cohesion.
    A 6-step run at 900 K checks the rise of the rock mass, 37 % (3-D) and 46 % (2-D).
    A restart that changes `porosity.model` is refused, with `force_restart_config` false and true.
    A restart from the checkpoint of step 4 repeats an 8-step run bit for bit.

    Command: `julia --project -t 2 -e 'using Test; include("test/test_porosity_henke.jl")'`
    With 2 threads the test also runs the threaded marker-to-grid mapping.

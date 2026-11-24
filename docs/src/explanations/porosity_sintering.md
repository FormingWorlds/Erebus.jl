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
The plastic iteration does not weigh it either.
The configuration check rejects a `phimin` of 0 or less, or above $10^{-3}$, with this model.
The adiabatic heating still uses the water porosity of rock.
There the bound keeps the solid fraction $1 - \phi$ within 0.1 % of 1.
The Stokes-Darcy solve holds no water at all, as the next section explains.
The columns `phi_max` and `phi_mean` of `telemetry.csv` give the maximum and the mean of $\phi_d$ over the rock markers, not the water porosity.
The processes that need pore water, melt or metal are therefore not allowed together with this model.
The configuration check rejects `henke2012` when reaction, venting, hydrothermal convection, accretion or telescoping is on.
It also rejects volatiles, pyrolysis kinetics, magma transport, magma degassing, core formation, iron heating (`hr_fe`) and metal partitioning.
It rejects the tracking of accessory metal phases (`phase_tracking`) and the volatile ice mixture (`volatile_mixture`) too.
It rejects the redox model (`redox`) as well, because it tracks metallic iron and gas species.
It rejects hydrofracture (`poroelasticity.hydrofracture`), because no pore fluid pressure can open a crack.
The rock must be dry: `materials.XWsolidm_init = [0.0, 0.0, nan]`.

---

## How does the Stokes-Darcy solve treat rock without pore water?

The short answer: with `model = "henke2012"` the solve holds no water, so no value of the water can move the rock.

The Stokes-Darcy solve has 2 pressures at each cell centre, the total pressure $P_t$ and the fluid pressure $P_f$, both in Pa.
With pore water, the solid continuity equation has a compaction term, and a second equation moves the water by Darcy flow.
Without pore water every cell centre is dry, in rock and in the sticky air around the body.
The pores of the sticky air hold no water either: their water porosity is also `phimin`.
A dry cell has no compaction term and no fluid pressure in its solid continuity equation.
Its drained compressibility is the solid compressibility $\beta_s$ (`poroelasticity.betasolid`), in $\text{Pa}^{-1}$, and its Biot-Willis coefficient is 0.
Its fluid equation is $P_f = P_t$, and no Darcy flux crosses its faces.
The solid continuity equation of a dry cell is then

$$\nabla \cdot \mathbf{v}_s + \beta_s \frac{P_t - P_t^{0}}{\Delta t} = 0.$$

Here $\mathbf{v}_s$ is the solid velocity in m/s, $P_t^{0}$ the total pressure of the previous step in Pa, and $\Delta t$ the time step in s.
The first step takes $\beta_s = 0$, so the rock is incompressible there.
The permeability `materials.kphim0`, the water viscosity and the water compressibility `poroelasticity.betafluid` therefore change no velocity, pressure or temperature.
The water porosity `phimin` enters only the adiabatic heating, through the solid fraction $1 - \phi$, and through it the temperature.
A restart with another `phimin` or `kphim0` moves the rock as a restart with the values of the checkpoint.
The markers still compute the Darcy resistance $\eta_f / k(\phi_{min})$, where $\eta_f$ is the water viscosity in Pa s and $k$ the permeability in $\text{m}^2$.
Its largest value for the default mantle is the ice viscosity $10^{12}$ Pa s over $k = 8.0 \times 10^{-24}$ $\text{m}^2$.
The configuration check rejects a `kphim0` that is not positive, and a `phimin` or `kphim0` so small that this resistance overflows.
It also rejects the matrix-free solver (`solver.hydromech_solver = "matrix_free"`), whose operator has no dry cells, and the iterative solver, which no test runs with dry cells.

The tests of the section *Which tests check this?* run a 200 km body for 4 steps with the rock reference permeability `kphim0` from $10^{-30}$ to $10^{-18}$ $\text{m}^2$ and `betafluid` from 0 to $10^{-5}$ $\text{Pa}^{-1}$.
The solid velocity, the temperature and $\phi_d$ change by less than 1e-10 of their size.

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
The fluid pressure of the dry solve equals the total pressure.
With it, the effective pressure would vanish, and rock of low cohesion would yield at a small shear stress.

---

## What does the model leave out?

The short answer: it leaves out some density terms, the shrinking of the body, the dry porosity in the adiabatic heating, the compaction work, a reset on melting and the sintered surface porosity.
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

## How deep is the porous shell?

The short answer: the shell ends where the peak temperature of the rock reached the sintering temperature of 700 K, and `porous_shell_depth` measures its depth on the rock markers.

The porous shell is the outer layer of a planetesimal where the rock never sintered.
Its markers keep the porosity of cold pressing, because they did not stay hot enough for long enough.
The function `porous_shell_depth` measures the depth of this layer on the rock markers.

A marker counts as sintered when sintering has lowered its porosity by at least $10^{-3}$, that is, when `dphi_sinter` $\geq 10^{-3}$.
Cold pressing does not count, because it changes `phi_d` but not `dphi_sinter`.
The function divides the slice into 16 equal sectors around the centre of the planetesimal.
It divides each sector into radial bins of width $\Delta r$, from the surface inward.
In each sector, the first bin whose median marker is sintered ends the shell.
The depth of the sector is the depth of the outer edge of that bin.
The result is the median of the 16 sector depths, with the 25 % and 75 % quartiles.

- Rock markers above the radius $R_p$ do not count.
  The pressure there is 0, so they never sinter.
- An empty bin is skipped.
- A sector without a sintered bin has the depth $R_p$: the whole sector is porous.
- With `measure = :phi_gt_04`, a marker counts as sintered when `phi_d` $\leq 0.4$.
  The depth is then that of the region with $\phi > 0.4$, the measure of Fig. 8 of Lichtenberg et al. (2016a).
  This measure gives a result only when the cold-pressed porosity at the surface is above 0.4.

On a simulation state, the bins are half the smaller grid spacing.
Here `cfg` is a configuration with `porosity.model = "henke2012"`:

```julia
using Erebus
st = simulation_loop(cfg; output_path=mktempdir())
shell = porous_shell_depth(st, GridCoordinates(cfg.grid), cfg.geometry)
shell.median, shell.q25, shell.q75   # depths in m
```

A model without the dry porosity group (`porosity.model = "none"`) gives an `ArgumentError`.

**Why the shell can be thin.**
The porous layer conducts heat much less well than compact rock.
At $\phi = 0.25$ the conductivity is 0.21 W/(m K), against 3 W/(m K) at $\phi = 0$.
A thin layer can therefore hold the drop from the sintering temperature of 700 K to the cold surface.
Below it, the rock is above the sintering temperature of 700 K and sinters.
Take a body with $R_p = 20$ km and $\phi_{init} = 0.25$ in the 1-D reference model, with steps of 2 kyr.
It keeps a shell of 380 m when it forms 1.5 Myr after the first solids.
It keeps a shell of 4.55 km when it forms at 2.5 Myr, because it has less $^{26}$Al and stays cooler.

---

## Which tests check this?

The validation case [Porous Shell](../validation/porous_shell.md) runs this model in the complete time loop.
It compares the shell depth and the maximum temperature with an independent 1-D model.
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
    A rock reference permeability `kphim0` of 1e-6 $\text{m}^2$ passes.
    A `kphim0` of 0, below 0 or NaN is rejected, and so are a `kphim0` of 1e-305 $\text{m}^2$ and a `phimin` of the smallest positive float, whose Darcy resistance overflows.
    The matrix-free and the iterative solver are rejected with `henke2012`; the direct solver passes.
    They check that the marker step and `dry_porosity_update` reject a field of a directly built `PorosityConfig` outside its domain.
    The marker step also rejects a bad `dt` or `psurface` before it changes a marker, also when no marker is rock.
    It checks `dt` and `psurface` in Float64, so a BigFloat of 1e400 is rejected and a BigFloat of 1e5 gives the Float64 result.
    The yield stress, `dry_porosity_update` and the cold pressing law check their pressure the same way.
    Testsets of `Dry porosity marker group and step` check the marker arrays and the setup pressure against the uniform-density limit.
    They check the residual of the hydrostatic balance and the setup pressure against the Stokes pressure of the first step, to 5 %.
    They check the setup pressure `pfm0` of the melting law and the 2-D Poisson profile of compressible powder.
    A 3-step run with melting checks that `pfm0` follows the total pressure, to 1e-6, and that the dry solve gives $P_f = P_t$.
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
    A water porosity `phim0` of 1e-4 and a permeability `kphim0` of 1e-12 $\text{m}^2$ give the same step, to 1e-10.
    A test helper repeats the step functions of `simulation_loop`, bit for bit.
    With `phimin` = 1e-3 and `kphim0` of 1e-9 or 1e-12 $\text{m}^2$, a 4-step run keeps the time step and the temperature to 1e-12, the shear and adiabatic heating stay below 1e-6 $\text{W/m}^3$, and each step takes 1 outer iteration.
    A 4-step run of a 200 km body changes the solid velocity, the temperature and $\phi_d$ by less than 1e-10 between 2 values of a water parameter.
    These are a rock `kphim0` of 1e-30 and 1e-18 $\text{m}^2$, an air `kphim0` of 1e-17 and 1e-9 $\text{m}^2$, and `betafluid` of 0 and 1e-5 $\text{Pa}^{-1}$.
    The runs use the direct 6-variable solve, the 4-variable solve with eliminated Darcy fluxes and the spherical slice.
    A rock `kphim0` of 1e-6 $\text{m}^2$ runs without a warning.
    With model none `betafluid` changes the solid velocity by more than 1e-6.
    Assembled alone, the dry rows give $P_f = P_t$ and zero Darcy flux on the faces of a dry cell, and the water of a dry cell changes no velocity.
    A mask without a dry cell gives the matrix of model none bit for bit.
    The plastic iteration gives $P_f$ the weight 0 in its Anderson mixing and leaves $P_f$ and the Darcy fluxes out of its cycle stop.
    The dry solve gives $P_f = P_t$, so no cell fractures and no fracture flag flips; a fluid pressure 2e5 Pa above $P_t$ fractures cells only with pore water.
    Without pore water the fluid pressure is carried with the solid velocity, and `time.dxymax_fluid` does not change the step of a melting run.
    Fluid and solid pressures 1 MPa per column off the solution leave the heating, the temperature and the convergence of a dry thermal step unchanged; with pore water they change them.
    Active yielding gives the same yielding nodes and step at `psurface` of 1e3 Pa and 1e8 Pa, in both gravity modes.
    The telemetry porosity counts crust markers and is 0 without rock markers.
    The plastic yield stress takes the pore pressure `psurface`: a 3-step run of a 200 km body with a cohesion of 1e4 Pa keeps the step of the default cohesion.
    A 6-step run at 900 K checks the rise of the rock mass, 37 % (3-D) and 46 % (2-D).
    A restart that changes `porosity.model` is refused, with `force_restart_config` false and true.
    A restart from the checkpoint of step 4 repeats an 8-step run bit for bit.
    A forced restart with `phimin` = 1e-8 and `kphim0` = 1e-9 $\text{m}^2$ gives the velocity, the step and $\phi_d$ of a restart with the checkpoint values, to 1e-10.

    Command: `julia --project -t 2 -e 'using Test; include("test/test_porosity_henke.jl")'`
    With 2 threads the test also runs the threaded marker-to-grid mapping.

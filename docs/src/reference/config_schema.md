# Configuration Schema Reference

This reference documents every parameter in `Erebus.jl` configuration files (`.toml`).

---

## `[grid]`

> [!NOTE]
> In the current release, grid resolution and domain dimensions are compiled into static array stencils. `validate_config` asserts that `xsize`, `ysize`, `Nx`, and `Ny` equal the compiled values in `src/constants.jl`.

| Parameter | Type | Default | Units | Description | Bounds / Invariant |
|:---|:---|:---|:---|:---|:---|
| `xsize` | `Float64` | `140000.0` | m | Total horizontal domain size | Must match compiled constant |
| `ysize` | `Float64` | `140000.0` | m | Total vertical domain size | Must match compiled constant |
| `Nx` | `Int` | `33` | - | Number of basic grid points in x | Must match compiled constant |
| `Ny` | `Int` | `33` | - | Number of basic grid points in y | Must match compiled constant |

---

## `[geometry]`

| Parameter | Type | Default | Units | Description | Bounds / Invariant |
|:---|:---|:---|:---|:---|:---|
| `rplanet` | `Float64` | `50000.0` | m | Outer radius of the planetesimal | Must match compiled constant |
| `rcrust` | `Float64` | `50000.0` | m | Boundary radius between core/mantle and crust | Must match compiled constant |
| `xcenter` | `Float64` | `70000.0` | m | Horizontal position of planetesimal center | $\in [0, \text{xsize}]$ |
| `ycenter` | `Float64` | `70000.0` | m | Vertical position of planetesimal center | $\in [0, \text{ysize}]$ |
| `psurface` | `Float64` | `1.0e+3` | Pa | Surface pressure anchor | $\ge 0$ |
| `spherical_metric` | `Bool` | `false` | - | Enable 3D spherical geometric metric heat source term in 2D Cartesian solver | `true` / `false` |
| `metric_regularization_cells` | `Float64` | `0.5` | grid units | Regularization radius at center $r \to 0$ in units of grid cell spacing $\Delta$ | $> 0$ |

---

## `[time]`

> [!NOTE]
> `TimeConfig` defines `yearlength = 31557600.0` s (Julian year: $365.25 \times 86400\text{ s}$). Simulation time parameters (`dt_initial`, `dt_longest`, `start_time`, `endtime`) are specified in years [yr] in configuration files. Default `start_time = 2.25e6` yr ($2.25\text{ Ma}$) and `endtime = 15.0e6` yr ($15.0\text{ Ma}$) reflect canonical CAI formation timing in the early Solar System.

| Parameter | Type | Default | Units | Description | Bounds |
|:---|:---|:---|:---|:---|:---|
| `dt_initial` | `Float64` | `3168.80878` | yr | Initial computational timestep (~$10^{11}\text{ s}$) | $> 0$ |
| `dt_longest` | `Float64` | `3168.80878` | yr | Maximum allowed computational timestep (~$10^{11}\text{ s}$) | $\ge \text{dt\_initial}$ |
| `dtcoefdn` | `Float64` | `0.5` | - | Factor to reduce timestep on non-convergence | $\in (0, 1)$ |
| `dtcoefup` | `Float64` | `1.2` | - | Factor to increase timestep on convergence | $> 1$ |
| `dtstep` | `Int` | `200` | - | Steps between increasing timestep | $\ge 1$ |
| `dxymax` | `Float64` | `0.05` | grid units | Maximum marker displacement per timestep | $> 0$ |
| `vpratio` | `Float64` | `0.333333333333` | - | Velocity weighting parameter | $\in [0, 1]$ |
| `DTmax` | `Float64` | `20.0` | K | Maximum allowed temperature change per step | $> 0$ |
| `yearlength` | `Float64` | `31557600.0` | s | Length of one Julian year ($365.25 \times 86400\text{ s}$) | $> 0$ |
| `start_time` | `Float64` | `2.25e+6` | yr | Simulation start time after CAIs (2.25 Ma) | $\ge 0$ |
| `endtime` | `Float64` | `15.0e+6` | yr | Total simulation end time (15 Ma) | $> \text{start\_time}$ |
| `start_step` | `Int` | `1` | - | Initial timestep counter index | $\ge 1$ |
| `n_steps` | `Int` | `10` | - | Total number of computational timesteps | $\ge 1$ |

---

## `[solver]`

| Parameter | Type | Default | Units | Description | Bounds |
|:---|:---|:---|:---|:---|:---|
| `titermax` | `Int` | `10000` | - | Maximum global thermochemical iterations | $\ge 1$ |
| `nplast` | `Int` | `100000` | - | Maximum plastic yielding iterations | $\ge \text{titermax}$ |
| `yerrmax` | `Float64` | `100.0` | - | Plastic yielding relative error tolerance | $> 0$ |
| `etawt` | `Float64` | `0.0` | - | Viscosity relaxation weight | $\in [0, 1)$ |
| `dphimax` | `Float64` | `100.01` | - | Maximum porosity change ratio per step | $> 1$ |
| `seed` | `Int` | `42` | - | Random seed for marker initialization | Any `Int` |
| `use_pardiso` | `Bool` | `false` | - | Enable Pardiso solver instead of UMFPACK | `true` / `false` |
| `etaphikoef` | `Float64` | `1.0` | - | Bulk viscosity scaling factor | $> 0$ |
| `etamin` | `Float64` | `1.0e+12` | Pa s | Lower shear viscosity cutoff | $> 0$ |
| `etamax` | `Float64` | `1.0e+23` | Pa s | Upper shear viscosity cutoff | $\ge \text{etamin}$ |

---

## `[poroelasticity]`

| Parameter | Type | Default | Units | Description | Bounds |
|:---|:---|:---|:---|:---|:---|
| `betasolid` | `Float64` | `0.0` | $\text{Pa}^{-1}$ | Solid silicate matrix compressibility (production: `2.5e-11`) | $\ge 0$ |
| `betafluid` | `Float64` | `0.0` | $\text{Pa}^{-1}$ | Pore fluid compressibility (production: `4.0e-10`) | $\ge 0$ |
| `phimin` | `Float64` | `1.0e-4` | - | Minimum porosity floor threshold | $0 < \phi_{\text{min}} < \phi_{\text{max}}$ |
| `phimax` | `Float64` | `0.9999` | - | Maximum porosity ceiling threshold | $\phi_{\text{min}} < \phi_{\text{max}} < 1$ |
| `hydrofracture` | `Bool` | `false` | - | Enable dynamic hydrofracturing permeability enhancement | `true` / `false` |
| `kappa_frac` | `Float64` | `1.0e+3` | - | Dimensionless hydrofracture permeability multiplier | $\ge 0$ |
| `gamma_frac` | `Float64` | `1.0` | - | Power-law exponent for overpressure scaling | $> 0$ |
| `k_frac_max` | `Float64` | `1.0e-9` | $\text{m}^2$ | Maximum fractured permeability ceiling | $> 0$ |

---

## `[thermodynamics]`

| Parameter | Type | Default | Units | Description | Bounds |
|:---|:---|:---|:---|:---|:---|
| `hr_al` | `Bool` | `true` | - | Enable 26Al radiogenic decay heating | `true` / `false` |
| `hr_fe` | `Bool` | `false` | - | Enable 60Fe radiogenic decay heating | `true` / `false` |
| `ratio_al` | `Float64` | `5.0e-5` | - | Initial 26Al/27Al isotope ratio at CAIs | $\in [0, 1]$ |
| `ratio_fe` | `Float64` | `1.0e-6` | - | Initial 60Fe/56Fe isotope ratio at CAIs | $\in [0, 1]$ |
| `E_al` | `Float64` | `5.0470e-13` | J | Decay energy per 26Al atom | $> 0$ |
| `f_al` | `Float64` | `1.9e+23` | atoms/kg | Abundance of 27Al atoms per unit mass | $> 0$ |
| `t_half_al` | `Float64` | `2.2614e+13` | s | 26Al half-life (717,000 yr) | $> 0$ |
| `E_fe` | `Float64` | `4.34e-13` | J | Decay energy per 60Fe atom | $> 0$ |
| `f_fe` | `Float64` | `1.957e+24` | atoms/kg | Abundance of 56Fe atoms per unit mass | $> 0$ |
| `t_half_fe` | `Float64` | `8.2635e+13` | s | 60Fe half-life (2.62 Myr) | $> 0$ |
| `tmsolidphase` | `Float64` | `1416.0` | K | Silicate rock solidus temperature | $> \text{tmfluidphase}$ |
| `tmfluidphase` | `Float64` | `273.0` | K | Water ice melting temperature | $> 0$ |
| `Lᶠ` | `Float64` | `333.55e+3` | J/kg | Latent heat of melting for water ice | $> 0$ |
| `phim0` | `Float64` | `0.2` | - | Standard reference porosity | $\in (0, 1)$ |
| `thermal_buoyancy` | `Bool` | `true` | - | Enable temperature-dependent fluid thermal buoyancy in Darcy flow | `true` / `false` |
| `fluid_viscosity_mode` | `Symbol` | `:arrhenius` | - | Pore fluid viscosity calculation mode | `:arrhenius` / `:constant` |
| `fluid_viscosity_Ea` | `Float64` | `15.0e+3` | J/mol | Activation energy for fluid viscous flow | $\ge 0$ |
| `fluid_viscosity_T0` | `Float64` | `293.15` | K | Reference temperature for fluid viscosity | $> 0$ |
| `fluid_viscosity_eta0` | `Float64` | `1.0e-3` | Pa s | Reference liquid water viscosity at T0 | $> 0$ |
| `surface_radiation` | `Bool` | `false` | - | Enable linearized Stefan-Boltzmann radiative cooling boundary condition | `true` / `false` |
| `emissivity` | `Float64` | `0.9` | - | Surface thermal emissivity $\epsilon$ | $\in [0, 1]$ |
| `sigma_sb` | `Float64` | `5.670374419e-8` | $\text{W}/(\text{m}^2\cdot\text{K}^4)$ | Stefan-Boltzmann radiation constant $\sigma$ | $> 0$ |

---

## `[materials]`

3-element vectors representing `[Index 1: Core, Index 2: Crust, Index 3: Sticky Air]`.

> [!WARNING]
> Eight material arrays are compiled into numerical stencils and cannot be modified without recompiling: `rhosolidm`, `rhofluidm`, `etasolidm`, `etasolidmm`, `etafluidm`, `etafluidmm`, `ksolidm`, and `kfluidm`. `validate_config` throws an `ArgumentError` if custom values differ from `src/constants.jl`. The remaining ten arrays can be modified freely.

| Parameter | Type | Default | Units | Status | Bounds | Description |
|:---|:---|:---|:---|:---|:---|:---|
| `rhosolidm` | `SVector{3}` | `[3300.0, 3300.0, 1.0]` | $\text{kg/m}^3$ | Compiled Constant | All $> 0$ | Solid matrix density |
| `rhofluidm` | `SVector{3}` | `[1000.0, 1000.0, 1.0]` | $\text{kg/m}^3$ | Compiled Constant | All $> 0$ | Pore fluid density |
| `etasolidm` | `SVector{3}` | `[1.0e19, 1.0e19, 1.0e16]` | Pa s | Compiled Constant | All $> 0$ | Solid reference shear viscosity |
| `etasolidmm` | `SVector{3}` | `[1.0e19, 1.0e19, 1.0e16]` | Pa s | Compiled Constant | All $> 0$ | Molten solid shear viscosity |
| `etafluidm` | `SVector{3}` | `[1.0e12, 1.0e12, 1.0e-3]` | Pa s | Compiled Constant | All $> 0$ | Unmelted fluid phase viscosity |
| `etafluidmm` | `SVector{3}` | `[1.0e-3, 1.0e-3, 1.0e-3]` | Pa s | Compiled Constant | All $> 0$ | Liquid water dynamic viscosity |
| `ksolidm` | `SVector{3}` | `[3.0, 3.0, 3000.0]` | $\text{W}/(\text{m}\cdot\text{K})$ | Compiled Constant | All $> 0$ | Solid thermal conductivity |
| `kfluidm` | `SVector{3}` | `[50.0, 50.0, 3000.0]` | $\text{W}/(\text{m}\cdot\text{K})$ | Compiled Constant | All $> 0$ | Fluid thermal conductivity |
| `rhocpsolidm` | `SVector{3}` | `[3.3e6, 3.3e6, 3.0e6]` | $\text{J}/(\text{m}^3\cdot\text{K})$ | Configurable | All $> 0$ | Volumetric solid heat capacity |
| `rhocpfluidm` | `SVector{3}` | `[1.0e6, 1.0e6, 3.0e6]` | $\text{J}/(\text{m}^3\cdot\text{K})$ | Configurable | All $> 0$ | Volumetric fluid heat capacity |
| `alphasolidm` | `SVector{3}` | `[3.0e-5, 3.0e-5, 0.0]` | $1/\text{K}$ | Configurable | All $\ge 0$ | Solid thermal expansion |
| `alphafluidm` | `SVector{3}` | `[5.0e-5, 5.0e-5, 0.0]` | $1/\text{K}$ | Configurable | All $\ge 0$ | Fluid thermal expansion |
| `gggsolidm` | `SVector{3}` | `[1.0e10, 1.0e10, 1.0e10]` | Pa | Configurable | All $> 0$ | Solid shear elastic modulus |
| `frictsolidm` | `SVector{3}` | `[0.6, 0.6, 0.0]` | - | Configurable | All $\ge 0$ | Internal friction coefficient |
| `cohessolidm` | `SVector{3}` | `[1.0e8, 1.0e8, 1.0e8]` | Pa | Configurable | All $> 0$ | Cohesion |
| `tenssolidm` | `SVector{3}` | `[6.0e7, 6.0e7, 6.0e7]` | Pa | Configurable | All $> 0$ | Tensile strength |
| `kphim0` | `SVector{3}` | `[1.0e-13, 1.0e-13, 1.0e-17]` | $\text{m}^2$ | Configurable | All $> 0$ | Reference permeability |
| `tkm0` | `SVector{3}` | `[170.0, 170.0, 170.0]` | K | Configurable | All $> 0$ | Initial temperature |

---

## `[output]`

| Parameter | Type | Default | Units | Description | Bounds |
|:---|:---|:---|:---|:---|:---|
| `output_dir` | `String` | `"output"` | - | Output directory path | Non-empty string |
| `savematstep` | `Int` | `10` | - | Checkpoint saving frequency | $\ge 1$ |
| `visstep` | `Int` | `1` | - | Visualization step cadence | $\ge 1$ |

---

## `[disk]`

Parameters controlling protoplanetary disk ambient temperature evolution and astronomical host star scalings.

| Parameter | Type | Default | Units | Description | Bounds |
|:---|:---|:---|:---|:---|:---|
| `enabled` | `Bool` | `false` | - | Enable disk ambient temperature evolution | `true` / `false` |
| `model` | `Symbol` | `:fixed` | - | Evolution model (`:fixed`, `:monotonic`, `:class1_to_class2`; `:class0_to_class2` alias) | `:fixed` / `:monotonic` / `:class1_to_class2` / `:class0_to_class2` |
| `t_ambient` | `Float64` | `170.0` | K | Constant background temperature for `:fixed` mode and ambient sink temperature when `surface_radiation = true` with `disk.enabled = false` | $> 0$ |
| `orbital_distance_au` | `Float64` | `2.5` | AU | Planetesimal heliocentric orbital distance $r$ | $> 0$ |
| `stellar_mass_msun` | `Float64` | `1.0` | $M_\odot$ | Host star mass $M_\star$ | $> 0$ |
| `t_cloud` | `Float64` | `30.0` | K | Molecular cloud background temperature floor $T_{\text{cloud}}$ | $> 0$ |
| `t_irr_1au` | `Float64` | `150.0` | K | Flared disk irradiation temperature at 1 AU for $1\,M_\odot$ | $> 0$ |
| `t_peak_1au` | `Float64` | `520.0` | K | Peak viscous dissipation temperature at 1 AU for $1\,M_\odot$ | $> 0$ |
| `t_peak_time_1au_myr` | `Float64` | `0.12` | Myr | Time of peak viscous dissipation at 1 AU for $1\,M_\odot$ | $> 0$ |
| `t_visc_0_myr` | `Float64` | `0.25` | Myr | Viscous dissipation power-law reference time $t_0$ | $> 0$ |
| `gamma` | `Float64` | `1.4` | - | Viscous dissipation power-law decay index $\gamma$ ($t > t_0$) | $> 0$ |
| `alpha` | `Float64` | `2.0` | - | Early infall heating rise power-law index $\alpha$ ($t \le t_{\text{peak}}$) | $> 0$ |
| `q_irr` | `Float64` | `0.42857142857142855` | - | Flared disk irradiation radial scaling exponent $q_{\text{irr}}$ ($3/7$) | $> 0$ |
| `q_visc` | `Float64` | `0.75` | - | Viscous dissipation radial scaling exponent $q_{\text{visc}}$ ($3/4$) | $> 0$ |
| `p_r_t` | `Float64` | `0.25` | - | Radial scaling exponent for peak heating time $p_{r,t}$ | $\ge 0$ |
| `p_m_irr` | `Float64` | `0.25` | - | Stellar mass scaling exponent for irradiation temperature $p_{M,\text{irr}}$ | $\ge 0$ |
| `p_m_visc` | `Float64` | `0.30` | - | Stellar mass scaling exponent for peak viscous temperature $p_{M,\text{visc}}$ | $\ge 0$ |
| `p_m_t` | `Float64` | `0.40` | - | Stellar mass scaling exponent for peak heating time $p_{M,t}$ | $\ge 0$ |
| `p_m_visc_decay` | `Float64` | `0.30` | - | Stellar mass scaling exponent for viscous dissipation time $p_{M,\text{visc,decay}}$ | $\ge 0$ |

---

## `[reaction]`

Parameters controlling hydrothermal water-rock hydration and dehydration reactions.

| Parameter | Type | Default | Units | Description | Bounds |
|:---|:---|:---|:---|:---|:---|
| `active` | `Bool` | `true` | - | Enable two-way hydrothermal reaction coupling | `true` / `false` |
| `hydration_active` | `Bool` | `true` | - | Enable serpentine hydration reaction pathway | `true` / `false` |
| `dehydration_active` | `Bool` | `true` | - | Enable serpentine dehydration reaction pathway | `true` / `false` |
| `hydration_mode` | `Int` | `1` | - | Hydration kinetics formulation mode | `1` |
| `dehydration_mode` | `Int` | `2` | - | Dehydration kinetics formulation mode | `2` |
| `dtreaction_hydration` | `Float64` | `1.0e10` | s | Timescale for serpentine hydration kinetics | $> 0$ |
| `dtreaction_dehydration` | `Float64` | `1.0e8` | s | Timescale for serpentine dehydration kinetics | $> 0$ |
| `delta_H` | `Float64` | `40000.0` | J/mol | Enthalpy of reaction | $> 0$ |
| `delta_S` | `Float64` | `60.0` | J/(mol K) | Entropy of reaction | $> 0$ |
| `A_I` | `Float64` | `1.0e-11` | $\text{s}^{-1}$ | Kinetic rate pre-factor | $> 0$ |
| `b_I` | `Float64` | `2.5e-4` | $1/\text{K}$ | Temperature sensitivity coefficient | $> 0$ |
| `c_I` | `Float64` | `543.0` | K | Equilibrium temperature parameter | $> 0$ |
| `Sxo_B` | `Float64` | `2.0e-11` | $\text{s}^{-1}$ | Reaction scale pre-factor | $> 0$ |
| `Tscl_B` | `Float64` | `10.0` | K | Temperature scale factor | $> 0$ |
| `To_B` | `Float64` | `293.0` | K | Reference temperature | $> 0$ |
| `alpha_relaxation` | `Float64` | `0.5` | - | Reaction rate under-relaxation factor | $\in (0, 1]$ |
| `pfcoeff` | `Float64` | `0.5` | - | Fluid pressure relaxation coefficient | $\in (0, 1]$ |
| `pferrmax` | `Float64` | `1.0e5` | Pa | Maximum fluid pressure iteration residual | $> 0$ |
| `p_cavitation` | `Float64` | `1.0e7` | Pa | Cavitation pressure limit | $> 0$ |

---

## `[melting]`

Parameters controlling silicate rock melting, latent heat buffering, and melt-weakened rheology.

| Parameter | Type | Default | Units | Description | Bounds |
|:---|:---|:---|:---|:---|:---|
| `active` | `Bool` | `false` | - | Enable silicate rock melting and magma rheology | `true` / `false` |
| `T_solidus` | `SVector{3}` | `[1400.0, 1400.0, NaN]` | K | Material solidus temperatures | All $> 0$ or `NaN` |
| `T_liquidus` | `SVector{3}` | `[1800.0, 1800.0, NaN]` | K | Material liquidus temperatures | $> T_{\text{solidus}}$ or `NaN` |
| `L_melt` | `Float64` | `4.0e5` | J/kg | Latent heat of silicate melting | $> 0$ |
| `rho_melt` | `Float64` | `2800.0` | $\text{kg/m}^3$ | Molten silicate magma density | $> 0$ |
| `alpha_eta` | `Float64` | `28.0` | - | Melt weakening exponential coefficient | $\ge 0$ |
| `phi_crit` | `Float64` | `0.4` | - | Critical melt fraction for crystal suspension transition | $\in (0, 1)$ |
| `eta_melt` | `Float64` | `10.0` | Pa s | Dynamic viscosity of pure silicate melt | $> 0$ |
| `dpdt_clapeyron` | `Float64` | `0.0` | K/Pa | Clapeyron slope for pressure-dependent solidus and liquidus | $\ge 0$ |
| `latent_heat_mode` | `Symbol` | `:apparent_cp` | - | Latent heat formulation mode | `:apparent_cp` |
| `soft_turbulence` | `Bool` | `false` | - | Enable regularized sub-grid soft turbulence thermal conductivity enhancement | `true` / `false` |
| `turb_exponent` | `Float64` | `0.333333333333` | - | Power-law exponent for viscosity ratio ($1/3$ for Solomatov 2007, $1/2$ for boundary layer scaling) | $> 0$ |
| `eta_fluid_silicate` | `Float64` | `100.0` | Pa s | Dynamic viscosity of turbulent silicate fluid magma | $> 0$ |
| `F_turb_start` | `Float64` | `0.30` | - | Melt fraction threshold for onset of turbulent conductivity enhancement | $0 \le F_{\text{start}} < F_{\text{end}} \le 1$ |
| `F_turb_end` | `Float64` | `0.50` | - | Melt fraction threshold for fully developed turbulent conductivity | $F_{\text{start}} < F_{\text{end}} \le 1$ |
| `dT_turb_min` | `Float64` | `10.0` | K | Minimum temperature contrast scale for thermal regularization | $> 0$ |
| `T_surface_ref` | `Float64` | `300.0` | K | Reference ambient/surface temperature for contrast evaluation | $> 0$ |
| `k_turb_cutoff` | `Float64` | `1.0e+6` | W/(m K) | Upper cutoff for turbulent thermal conductivity | $> k_{\text{turb,floor}}$ |
| `k_turb_floor` | `Float64` | `1.0e-3` | W/(m K) | Lower cutoff floor for regularized thermal conductivity | $> 0$ |



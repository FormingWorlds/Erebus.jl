# Configuration Schema Reference

This reference documents every parameter in `Erebus.jl` configuration files (`.toml`).

---

## `[grid]`

Grid resolution and domain dimensions are configured per simulation run and constructed dynamically via `GridCoordinates(cfg.grid)`.

| Parameter | Type | Default | Units | Description | Bounds / Invariant |
|:---|:---|:---|:---|:---|:---|
| `xsize` | `Float64` | `140000.0` | m | Total horizontal domain size | $> 0$ |
| `ysize` | `Float64` | `140000.0` | m | Total vertical domain size | $> 0$ |
| `Nx` | `Int` | `33` | - | Number of basic grid points in x | $\ge 3$ |
| `Ny` | `Int` | `33` | - | Number of basic grid points in y | $\ge 3$ |

---

## `[geometry]`

| Parameter | Type | Default | Units | Description | Bounds / Invariant |
|:---|:---|:---|:---|:---|:---|
| `rplanet` | `Float64` | `50000.0` | m | Outer radius of the planetesimal | $> 0$ |
| `rcrust` | `Float64` | `50000.0` | m | Boundary radius between core/mantle and crust | $\in (0, \text{rplanet}]$ |
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
> Eight material arrays are compiled into numerical stencils and cannot be modified without recompiling: `rhosolidm`, `rhofluidm`, `etasolidm`, `etasolidmm`, `etafluidm`, `etafluidmm`, `ksolidm`, and `kfluidm`. `validate_config` throws an `ArgumentError` if custom values differ from `src/constants.jl`. The remaining eleven property arrays can be configured freely.

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
| `XWsolidm_init` | `SVector{3}` | `[0.5, 0.5, NaN]` | - | Configurable | All $\ge 0$ or `NaN` | Initial solid water fraction |

---

## `[output]`

| Parameter | Type | Default | Units | Description | Bounds |
|:---|:---|:---|:---|:---|:---|
| `output_dir` | `String` | `"output"` | - | Output directory path | Non-empty string |
| `savematstep` | `Int` | `10` | - | Checkpoint saving frequency | $\ge 1$ |
| `visstep` | `Int` | `1` | - | Visualization step cadence | $\ge 1$ |
| `restart_from` | `String` | `""` | - | Checkpoint JLD2 file path to resume simulation from | File path or empty string |

---

## `[disk]`

Parameters controlling protoplanetary disk ambient temperature evolution, gas dispersal, and astronomical host star scalings.

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
| `t_dispersal_myr` | `Float64` | `3.0` | Myr | Gas disk dispersal midpoint time | $> 0$ |
| `dt_dispersal_myr` | `Float64` | `0.1` | Myr | Gas disk dispersal transition half-width | $> 0$ |
| `p_amb_disk` | `Float64` | `10.0` | Pa | Nebular gas ambient pressure before dispersal | $> 0$ |
| `p_amb_space` | `Float64` | `1.0e-4` | Pa | Interplanetary vacuum ambient pressure floor | $\ge 0$ |
| `albedo` | `Float64` | `0.06` | - | Planetesimal surface Bond albedo | $\in [0, 1)$ |
| `t_eq_custom` | `Float64` | `NaN` | K | Custom solar equilibrium surface temperature (`NaN` = use solar scaling) | `NaN` or $> 0$ |
| `dispersal_active` | `Bool` | `false` | - | Enable dynamic nebular gas dispersal pressure decay | `true` / `false` |

---

## `[reaction]`

Parameters controlling hydrothermal water-rock hydration and dehydration reactions.

| Parameter | Type | Default | Units | Description | Bounds |
|:---|:---|:---|:---|:---|:---|
| `active` | `Bool` | `true` | - | Enable two-way hydrothermal reaction coupling | `true` / `false` |
| `hydration_active` | `Bool` | `true` | - | Enable serpentine hydration reaction pathway | `true` / `false` |
| `dehydration_active` | `Bool` | `true` | - | Enable serpentine dehydration reaction pathway | `true` / `false` |
| `hydration_mode` | `Int` | `1` | - | Hydration kinetics formulation mode | `1, 2, 3, 9` |
| `dehydration_mode` | `Int` | `2` | - | Dehydration kinetics formulation mode | `1, 2, 3, 9` |
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
| `pfcoeff` | `Float64` | `0.5` | - | Fluid pressure relaxation coefficient | $\in [0, 1]$ |
| `pferrmax` | `Float64` | `1.0e5` | Pa | Maximum fluid pressure iteration residual | $> 0$ |
| `p_cavitation` | `Float64` | `1.0e7` | Pa | Cavitation pressure limit | $> 0$ |

---

<a id="melting"></a>
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

---

## `[venting]`

Parameters controlling planetesimal surface volatile venting, ice sealing, and hydrofracture breaching.

| Parameter | Type | Default | Units | Description | Bounds |
|:---|:---|:---|:---|:---|:---|
| `active` | `Bool` | `false` | - | Enable surface volatile venting sink | `true` / `false` |
| `mode` | `Symbol` | `:darcy_sink` | - | Venting activation mode | `:darcy_sink` / `:hydrofracture_gated` |
| `k_vent` | `Float64` | `1.0e-11` | $\text{m}^2$ | Surface venting boundary permeability | $> 0$ |
| `conductance_factor` | `Float64` | `1.0` | - | Dimensionless boundary conductance multiplier | $> 0$ |
| `L_sublimation` | `Float64` | `2.83e6` | J/kg | Latent heat of ice sublimation | $> 0$ |
| `latent_cooling` | `Bool` | `true` | - | Enable volatile sublimation latent heat cooling | `true` / `false` |
| `ice_sealing` | `Bool` | `false` | - | Enable cryogenic pore ice permeability reduction | `true` / `false` |
| `t_freeze` | `Float64` | `273.15` | K | Water freezing temperature threshold | $> 0$ |
| `dt_seal` | `Float64` | `10.0` | K | Exponential ice sealing temperature scale | $> 0$ |
| `k_seal_min_ratio` | `Float64` | `1.0e-6` | - | Minimum cryogenic permeability residual ratio | $\in (0, 1]$ |
| `species` | `Symbol` | `:H2O` | - | Primary vented volatile gas species | `:H2O`, `:H2`, `:N2`, `:NH3`, `:CO`, `:CO2`, `:CH4`, `:H2S`, `:S2`, `:SO2` |

---

## `[volatiles]`

Parameters controlling multi-species volatile solubility in silicate melt and primordial organic devolatilization.

| Parameter | Type | Default | Units | Description | Bounds |
|:---|:---|:---|:---|:---|:---|
| `active` | `Bool` | `false` | - | Enable multi-species volatile solubility and chemistry | `true` / `false` |
| `fO2_delta_IW` | `Float64` | `-1.0` | log10 units | Redox state relative to iron-wüstite buffer ($\Delta\text{IW}$) | $\in [-50, 50]$ |
| `water_solubility_coeff` | `Float64` | `0.40` | $\text{wt}\% / \text{MPa}^{0.5}$ | Burnham low-pressure water solubility coefficient $A_s$ | $> 0$ |
| `water_law` | `Symbol` | `:burnham_dixon` | - | Water solubility law (`:burnham_dixon`, `:sossi_peridotite`, `:basalt_dixon`, `:newcombe_lunar`) | valid symbol |
| `h2_active` | `Bool` | `false` | - | Enable molecular $\text{H}_2$ dissolution in silicate melt | `true` / `false` |
| `h2_law` | `Symbol` | `:hirschmann2012` | - | Molecular $\text{H}_2$ solubility law (`:hirschmann2012`, `:gaillard2003`) | valid symbol |
| `nitrogen_law` | `Symbol` | `:dasgupta2022` | - | Nitrogen solubility law (`:dasgupta2022`, `:libourel2003`) | valid symbol |
| `nitrogen_henry_coeff` | `Float64` | `0.40` | ppm / bar | Henry coefficient $K_h$ for molecular $\text{N}_2$ dissolution | $> 0$ |
| `nitrogen_nitride_capacity` | `Float64` | `1.0e-3` | $\text{wt}\% / \text{bar}^{0.5}$ | Chemical nitride capacity $C_{\text{nitride}}$ | $> 0$ |
| `t_organic_devol` | `Float64` | `550.0` | K | Characteristic midpoint temperature $T_{\text{devol}}$ for organic devolatilization | $> 0$ |
| `dt_organic_devol` | `Float64` | `50.0` | K | Transition temperature scale $\Delta T$ for organic devolatilization | $> 0$ |
| `organic_n_initial_ppm` | `Float64` | `500.0` | ppm | Initial primordial organic nitrogen concentration in rocky core | $\ge 0$ |
| `initial_water_wtpct` | `Float64` | `1.0` | wt% | Initial water concentration in solid silicate matrix | $\ge 0$ |
| `initial_carbon_ppm` | `Float64` | `500.0` | ppm | Initial carbon concentration in solid silicate matrix | $\ge 0$ |
| `initial_nitrogen_ppm` | `Float64` | `50.0` | ppm | Initial dissolved nitrogen concentration in solid silicate matrix | $\ge 0$ |
| `initial_sulfur_ppm` | `Float64` | `1000.0` | ppm | Initial sulfur concentration in solid silicate matrix | $\ge 0$ |
| `carbon_active` | `Bool` | `false` | - | Enable carbon species solubility ($\text{CO}, \text{CH}_4, \text{CO}_2$) | `true` / `false` |
| `co_law` | `Symbol` | `:armstrong2015` | - | Carbon monoxide solubility law (`:armstrong2015`, `:yoshioka2019_morb`) | valid symbol |
| `ch4_law` | `Symbol` | `:ardia2013` | - | Methane solubility law (`:ardia2013`) | valid symbol |
| `co2_law` | `Symbol` | `:dixon1995` | - | Carbon dioxide carbonate solubility law (`:dixon1995`) | valid symbol |
| `graphite_saturation` | `Bool` | `true` | - | Enforce graphite saturation ceiling on carbon fugacities | `true` / `false` |
| `sulfur_active` | `Bool` | `false` | - | Enable sulfur species solubility in silicate melt | `true` / `false` |
| `sulfide_law` | `Symbol` | `:boulliung2023` | - | Sulfide solubility law (`:boulliung2023`, `:gaillard2022`) | valid symbol |
| `sulfide_melt` | `Symbol` | `:basalt` | - | Silicate melt composition for Boulliung (`:basalt`, `:andesite`, `:trachybasalt`) | valid symbol |
| `include_sulfate` | `Bool` | `false` | - | Include sulfate capacity at oxidizing conditions | `true` / `false` |
| `scss_active` | `Bool` | `true` | - | Cap dissolved sulfur at sulfide saturation (SCSS) | `true` / `false` |
| `scss_law` | `Symbol` | `:smythe2017` | - | SCSS formulation (`:smythe2017`, `:oneill2002`) | valid symbol |
| `melt_feo_wtpct` | `Float64` | `10.0` | wt% | Silicate melt $\text{FeO}$ concentration for SCSS calculation | $\ge 0$ |
| `x_sio2` | `Float64` | `0.56` | - | Silicate melt $\text{SiO}_2$ mole fraction for nitrogen solubility | $\in [0, 1]$ |
| `x_al2o3` | `Float64` | `0.11` | - | Silicate melt $\text{Al}_2\text{O}_3$ mole fraction for nitrogen solubility | $\in [0, 1]$ |
| `x_tio2` | `Float64` | `0.01` | - | Silicate melt $\text{TiO}_2$ mole fraction for nitrogen solubility | $\in [0, 1]$ |

---

## `[escape]`

Parameters controlling planetary atmospheric accumulation, kinetic Jeans escape, and surface pressure feedback.

| Parameter | Type | Default | Units | Description | Bounds |
|:---|:---|:---|:---|:---|:---|
| `active` | `Bool` | `false` | - | Enable atmospheric inventory evolution and Jeans escape | `true` / `false` |
| `M_planet` | `Float64` | `1.309e+18` | kg | Planetesimal mass for gravitational potential | $> 0$ |
| `R_planet` | `Float64` | `50000.0` | m | Planetesimal surface radius for atmospheric surface pressure | $> 0$ |
| `T_exobase` | `Float64` | `200.0` | K | Exobase temperature for Maxwellian thermal velocity | $> 0$ |
| `R_exobase` | `Float64` | `50000.0` | m | Exobase radius for escape flux surface integration | $\ge \text{R\_planet}$ |
| `species` | `Symbol` | `:H2O` | - | Primary outgassed volatile species for kinetic escape | `:H2O`, `:H2`, `:N2`, `:NH3`, `:CO`, `:CO2`, `:CH4`, `:H2S`, `:S2`, `:SO2` |
| `multi_species` | `Bool` | `false` | - | Enable coupled 10-species atmospheric inventory and escape tracking | `true` / `false` |
| `species_list` | `Vector{Symbol}` | `[:H2O, :H2, :CO, :CO2, :CH4, :N2, :NH3, :H2S, :S2, :SO2]` | - | Volatile species tracked for multi-species escape | valid symbols |
| `gamma` | `Float64` | `1.4` | - | Atmospheric adiabatic index (ratio of specific heats) | $> 0$ |
| `hydrodynamic` | `Bool` | `true` | - | Check hydrodynamic energy-limited escape rate ceiling | `true` / `false` |

---

## `[retention]`

Parameters controlling thermodynamic volatile retention floors in nominally anhydrous minerals (NAMs) and low-temperature surface venting drainage coupling.

| Parameter | Type | Default | Units | Description | Bounds |
|:---|:---|:---|:---|:---|:---|
| `active` | `Bool` | `false` | - | Enable thermodynamic volatile retention floors and venting drainage (requires `[volatiles]` active = true) | `true` / `false` |
| `h2o_retention_ppm` | `Float64` | `50.0` | ppmw | Subsolidus water retention floor in nominally anhydrous minerals (NAMs) | $\ge 0$ |
| `carbon_retention_ppm` | `Float64` | `50.0` | ppmw | Subsolidus carbon retention floor in refractory graphite/matrix | $\ge 0$ |
| `nitrogen_retention_ppm` | `Float64` | `5.0` | ppmw | Subsolidus nitrogen retention floor in crystalline silicates | $\ge 0$ |
| `sulfur_retention_ppm` | `Float64` | `100.0` | ppmw | Subsolidus sulfur retention floor in refractory sulfides/matrix | $\ge 0$ |
| `T_solidus_ref` | `Float64` | `1400.0` | K | Reference solidus temperature for retention depletion | $> 0$ |
| `dT_retention` | `Float64` | `200.0` | K | Temperature scale for supersolidus retention floor decay | $> 0$ |
| `retention_law` | `Symbol` | `:nams_exponential` | - | Retention floor law (`:constant_floor`, `:linear_melt_blend`, `:nams_exponential`) | valid symbol |
| `venting_drainage_active` | `Bool` | `true` | - | Drain mobile dissolved marker volatiles during surface venting | `true` / `false` |
| `chi_vent` | `Float64` | `1.0` | - | Volatile venting extraction efficiency factor | $\in [0, 1]$ |

---

## `[coreformation]`

Parameters controlling iron core formation, porous metal percolation, Stokes droplet settling, segregation dissipation heating, and mass conservation.

| Parameter | Type | Default | Units | Description | Bounds |
|:---|:---|:---|:---|:---|:---|
| `percolation_active` | `Bool` | `false` | - | Enable Darcy percolation of liquid Fe-FeS through solid silicate matrix | `true` / `false` |
| `settling_active` | `Bool` | `false` | - | Enable Stokes settling of liquid metal droplets in magma ocean | `true` / `false` |
| `rho_metal` | `Float64` | `5450.0` | $\text{kg/m}^3$ | Density of liquid metal (Fe-FeS) phase | $> \rho_{\text{silicate}}$ |
| `rho_metal_solid` | `Float64` | `5700.0` | $\text{kg/m}^3$ | Density of solid metal (Fe-FeS) phase | $> \rho_{\text{metal}}$ |
| `sulfur_fraction` | `Float64` | `0.31` | - | Sulfur mass fraction in Fe-FeS liquid metal phase ($w_S$) | $\in [0.0, 0.40]$ |
| `metal_density_mode` | `Symbol` | `:sanloup2000` | - | Liquid metal density equation of state (`:sanloup2000`, `:morard2014`, `:constant`) | valid symbol |
| `L_metal` | `Float64` | `2.7e5` | J/kg | Latent heat of melting for Fe-FeS eutectic mixture | $\ge 0$ |
| `eta_metal` | `Float64` | `1.0e-2` | Pa s | Dynamic viscosity of liquid metal phase | $> 0$ |
| `k_metal` | `Float64` | `40.0` | W/(m K) | Thermal conductivity of metallic phase | $> 0$ |
| `rhocp_metal` | `Float64` | `4.0e6` | $\text{J/(m}^3\text{ K)}$ | Volumetric heat capacity of metallic phase | $> 0$ |
| `Xfe_bulk` | `Float64` | `0.20` | - | Initial bulk metal volume fraction in planet interior | $\in [0, \phi_{\text{pack}}]$ |
| `phi_pack` | `Float64` | `0.65` | - | Maximum packing fraction of segregated core metal | $\in (0, 1]$ |
| `T_eutectic` | `Float64` | `1213.0` | K | Fe-FeS binary eutectic melting temperature | $> 0$ |
| `dT_metal` | `Float64` | `50.0` | K | Temperature range for complete metal melting | $> 0$ |
| `k_metal_ref` | `Float64` | `1.0e-9` | $\text{m}^2$ | Reference permeability for liquid metal percolation | $> 0$ |
| `perm_exponent` | `Float64` | `3.0` | - | Power-law exponent for metal permeability function | $> 0$ |
| `phi_crit_perc` | `Float64` | `0.05` | - | Critical porosity threshold for percolation connectivity | $\in [0, \phi_{\text{pack}})$ |
| `phi_residual` | `Float64` | `0.02` | - | Residual trapped metal volume fraction | $\in [0, \phi_{\text{crit\_perc}}]$ |
| `phi0` | `Float64` | `0.1` | - | Reference porosity for liquid metal Kozeny-Carman permeability | $\in (0, 1)$ |
| `droplet_size_mode` | `Symbol` | `:capillary_mean` | - | Metal droplet diameter calculation mode (`:capillary_mean`, `:bond_mean`, `:weber_mean`, `:weber_turbulent`, `:fixed`) | valid symbol |
| `droplet_diameter_fixed` | `Float64` | `5.0e-3` | m | Fixed droplet diameter when `droplet_size_mode = :fixed` | $> 0$ |
| `sigma_metal_silicate` | `Float64` | `1.0` | N/m | Metal-silicate interfacial surface tension | $> 0$ |
| `We_crit` | `Float64` | `10.0` | - | Critical Weber number for droplet breakup | $> 0$ |
| `hindered_exponent` | `Float64` | `4.5` | - | Richardson-Zaki hindered settling power-law exponent | $\ge 0$ |
| `hadamard_rybczynski` | `Bool` | `false` | - | Enable fluid droplet internal circulation correction factor | `true` / `false` |
| `F_settle_start` | `Float64` | `0.40` | - | Silicate melt fraction where settling commences | $\in [0, 1]$ |
| `F_perc_end` | `Float64` | `0.50` | - | Silicate melt fraction where percolation terminates | $\in [F_{\text{settle\_start}}, 1]$ |
| `segregation_heating` | `Bool` | `true` | - | Enable gravitational potential energy dissipation heating | `true` / `false` |
| `cfl_settling` | `Float64` | `0.5` | - | Courant-Friedrichs-Lewy stability safety factor for subcycling | $\in (0, 1]$ |
| `max_subcycles` | `Int` | `2000` | - | Maximum allowed subcycles per hydrodynamic step | $\ge 1$ |

---

### `[metal_partition]` - Metal-Silicate Volatile Partitioning & Transport

The `[metal_partition]` section controls thermodynamic volatile exchange between molten iron alloy and silicate melt, donor-cell advective volatile transport during core segregation, and dynamic liquid metal density coupling.

```toml
[metal_partition]
active = false
model_carbon = "grewal2019"
model_nitrogen = "grewal2019"
model_hydrogen = "clesi2018"
model_sulfur = "boujibar2014"
D_H_const = 0.5
D_C_const = 500.0
D_N_const = 20.0
D_S_const = 200.0
equilibration_rate = 1.0
dynamic_sulfur_density = true
D_min = 1.0e-4
D_max = 1.0e5
initial_metal_h_ppm = 0.0
initial_metal_c_ppm = 0.0
initial_metal_n_ppm = 0.0
initial_metal_s_ppm = 0.0
core_radius_fraction = 0.5
phi_core_threshold = 0.40
```

| Parameter | Type | Default | Units | Description | Bounds / Options |
|:----------|:-----|:--------|:------|:------------|:-----------------|
| `active` | `Bool` | `false` | - | Enable metal-silicate volatile partitioning and core segregation transport | `true` / `false` |
| `model_carbon` | `Symbol` | `:grewal2019` | - | Carbon partition model parameterization | `:constant`, `:grewal2019`, `:fischer2020` |
| `model_nitrogen` | `Symbol` | `:grewal2019` | - | Nitrogen partition model parameterization | `:constant`, `:grewal2019` |
| `model_hydrogen` | `Symbol` | `:clesi2018` | - | Hydrogen partition model parameterization | `:constant`, `:clesi2018` |
| `model_sulfur` | `Symbol` | `:boujibar2014` | - | Sulfur partition model parameterization | `:constant`, `:boujibar2014` |
| `D_H_const` | `Float64` | `0.5` | - | Constant partition coefficient for hydrogen | $\ge 0$ |
| `D_C_const` | `Float64` | `500.0` | - | Constant partition coefficient for carbon | $\ge 0$ |
| `D_N_const` | `Float64` | `20.0` | - | Constant partition coefficient for nitrogen | $\ge 0$ |
| `D_S_const` | `Float64` | `200.0` | - | Constant partition coefficient for sulfur | $\ge 0$ |
| `equilibration_rate` | `Float64` | `1.0` | - | Kinetic equilibration fraction per timestep pass ($\alpha_{\text{eq}}$) | $\in [0, 1]$ |
| `dynamic_sulfur_density` | `Bool` | `true` | - | Evaluate liquid metal density using marker sulfur mass fraction $w_S$ | `true` / `false` |
| `D_min` | `Float64` | `1.0e-4` | - | Numerical floor on partition coefficients | $> 0$ |
| `D_max` | `Float64` | `1.0e5` | - | Numerical ceiling on partition coefficients | $\ge D_{\text{min}}$ |
| `initial_metal_h_ppm` | `Float64` | `0.0` | ppmw | Initial hydrogen concentration in primordial metallic phase | $\ge 0$ |
| `initial_metal_c_ppm` | `Float64` | `0.0` | ppmw | Initial carbon concentration in primordial metallic phase | $\ge 0$ |
| `initial_metal_n_ppm` | `Float64` | `0.0` | ppmw | Initial nitrogen concentration in primordial metallic phase | $\ge 0$ |
| `initial_metal_s_ppm` | `Float64` | `0.0` | ppmw | Initial sulfur concentration in primordial metallic phase | $\ge 0$ |
| `core_radius_fraction` | `Float64` | `0.5` | - | Fractional planet radius defining central core region ($r_c / R_{\text{planet}}$) | $\in [0, 1]$ |
| `phi_core_threshold` | `Float64` | `0.40` | - | Metal volume fraction threshold for core membership | $\in [0, 1]$ |

---

## `[phase_tracking]`

Normative accessory mineral tracking and meteorite diagnostic parameters configure sub-eutectic stoichiometric allocation of S, P, C, and N into solid accessory phases (troilite $\text{FeS}$, schreibersite $(\text{Fe,Ni})_3\text{P}$, cohenite $(\text{Fe,Ni})_3\text{C}$, graphite $\text{C}$, nitrides $\text{Fe}_4\text{N}/\text{CrN}/\text{TiN}$) and residual metallic matrix, as well as thermal dissolution across the eutectic transition ($T_{\text{eutectic}} \approx 1213\text{ K}$). Enabling `[phase_tracking] active = true` requires active core formation (`[coreformation] percolation_active = true` or `settling_active = true`, or `[metal_partition] active = true`) to provide the metallic carrier phase, and enforces that `phase_tracking.T_eutectic` matches `coreformation.T_eutectic`.

```toml
[phase_tracking]
active = false
T_eutectic = 1213.0
dT_transition = 50.0
bulk_P_ppm = 1000.0
schreibersite_ni_frac = 0.25
cohenite_carbide_max = 0.0667
nitride_mode = "roaldite"
track_regional_modes = true
r_core_norm = 0.5
r_mantle_norm = 0.85
```

| Parameter | Type | Default | Units | Description | Bounds / Options |
|:----------|:-----|:--------|:------|:------------|:-----------------|
| `active` | `Bool` | `false` | - | Enable normative accessory mineral tracking | `true` / `false` |
| `T_eutectic` | `Float64` | `1213.0` | K | Metallic Fe-FeS eutectic temperature | $> 0$ |
| `dT_transition` | `Float64` | `50.0` | K | Eutectic phase dissolution temperature interval | $> 0$ |
| `bulk_P_ppm` | `Float64` | `1000.0` | ppmw | Bulk phosphorus concentration in metallic alloy | $\ge 0$ |
| `schreibersite_ni_frac` | `Float64` | `0.25` | - | Nickel molar fraction in schreibersite $(\text{Fe,Ni})_3\text{P}$ | $\in [0, 1]$ |
| `cohenite_carbide_max` | `Float64` | `0.0667` | - | Maximum carbon mass fraction before graphite saturation | $\in (0, 1]$ |
| `nitride_mode` | `Symbol` | `:roaldite` | - | Mineral stoichiometry for accessory nitrides | `:roaldite`, `:carlsbergite`, `:osbornite` |
| `track_regional_modes` | `Bool` | `true` | - | Calculate regional modal mineral distributions and classifications | `true` / `false` |
| `r_core_norm` | `Float64` | `0.5` | - | Normalized radial boundary for planetesimal core ($r / R_{\text{planet}}$) | $\in (0, 1)$ |
| `r_mantle_norm` | `Float64` | `0.85` | - | Normalized radial boundary for planetesimal mantle ($r / R_{\text{planet}}$) | $\in (\text{r\_core\_norm}, 1]$ |

---

## Configuration Loading and Synchronization

### File and String Input

The function `load_config` accepts either a filepath to a TOML configuration file or a raw string containing TOML content:

```julia
# Load from file path
cfg = load_config("configs/hydrothermal_benchmark.toml")

# Load from inline TOML string
cfg = load_config("""
[grid]
xsize = 140000.0
ysize = 140000.0
Nx = 33
Ny = 33
""")
```

### Automatic Radius Synchronization

When `escape.R_planet` is not explicitly defined in an input configuration and `geometry.rplanet` differs from the default radius ($50000.0\text{ m}$), `load_config` automatically initializes both `escape.R_planet` and `escape.R_exobase` to match `geometry.rplanet`. If `escape.R_planet` is explicitly defined without `escape.R_exobase`, `escape.R_exobase` retains its default value. When `escape.active = true`, `validate_config` enforces that `escape.R_planet` matches `geometry.rplanet` within 1% relative tolerance (`isapprox(escape.R_planet, geometry.rplanet; rtol=0.01)`).

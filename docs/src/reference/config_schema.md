# Configuration Schema Reference

This reference documents every parameter in `Erebus.jl` configuration files (`.toml`).

---

## `[grid]`

| Parameter | Type | Default | Units | Description | Bounds |
|:---|:---|:---|:---|:---|:---|
| `xsize` | `Float64` | `14000.0` | m | Total horizontal domain size | $> 0$ |
| `ysize` | `Float64` | `14000.0` | m | Total vertical domain size | $> 0$ |
| `Nx` | `Int` | `15` | - | Number of basic grid points in x | $\ge 3$ |
| `Ny` | `Int` | `15` | - | Number of basic grid points in y | $\ge 3$ |

---

## `[geometry]`

| Parameter | Type | Default | Units | Description | Bounds |
|:---|:---|:---|:---|:---|:---|
| `rplanet` | `Float64` | `5000.0` | m | Outer radius of the planetesimal | $> 0$ |
| `rcrust` | `Float64` | `4800.0` | m | Boundary radius between core/mantle and crust | $0 < r_{\text{crust}} \le r_{\text{planet}}$ |
| `xcenter` | `Float64` | `7000.0` | m | Horizontal position of planetesimal center | $\in [0, \text{xsize}]$ |
| `ycenter` | `Float64` | `7000.0` | m | Vertical position of planetesimal center | $\in [0, \text{ysize}]$ |
| `psurface` | `Float64` | `1.0e+3` | Pa | Surface pressure anchor | $\ge 0$ |

---

## `[time]`

| Parameter | Type | Default | Units | Description | Bounds |
|:---|:---|:---|:---|:---|:---|
| `dt_initial` | `Float64` | `1.0e+11` | s | Initial computational timestep | $> 0$ |
| `dt_longest` | `Float64` | `1.0e+11` | s | Maximum allowed computational timestep | $\ge \text{dt\_initial}$ |
| `dtcoefdn` | `Float64` | `0.5` | - | Factor to reduce timestep on non-convergence | $\in (0, 1)$ |
| `dtcoefup` | `Float64` | `1.2` | - | Factor to increase timestep on convergence | $> 1$ |
| `dtstep` | `Int` | `200` | - | Steps between increasing timestep | $\ge 1$ |
| `dxymax` | `Float64` | `0.05` | grid units | Maximum marker displacement per timestep | $> 0$ |
| `vpratio` | `Float64` | `0.3333` | - | Velocity weighting between grid and markers | $\in [0, 1]$ |
| `DTmax` | `Float64` | `20.0` | K | Maximum allowed temperature change per step | $> 0$ |
| `start_time` | `Float64` | `7.0965e+13` | s | Simulation start time after CAIs (2.25 Ma) | $\ge 0$ |
| `endtime` | `Float64` | `4.731e+14` | s | Total simulation end time (15 Ma) | $> \text{start\_time}$ |
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

---

## `[materials]`

3-element vectors representing `[Index 1: Core, Index 2: Crust, Index 3: Sticky Air]`.

| Parameter | Type | Default | Units | Description |
|:---|:---|:---|:---|:---|
| `rhosolidm` | `SVector{3}` | `[3300.0, 3300.0, 1.0]` | $\text{kg/m}^3$ | Solid matrix density |
| `rhofluidm` | `SVector{3}` | `[1000.0, 1000.0, 1.0]` | $\text{kg/m}^3$ | Pore fluid density |
| `etasolidm` | `SVector{3}` | `[1.0e19, 1.0e19, 1.0e16]` | Pa s | Solid reference shear viscosity |
| `etasolidmm` | `SVector{3}` | `[1.0e19, 1.0e19, 1.0e16]` | Pa s | Molten solid shear viscosity |
| `etafluidm` | `SVector{3}` | `[1.0e12, 1.0e12, 1.0e-3]` | Pa s | Unmelted fluid phase viscosity |
| `etafluidmm` | `SVector{3}` | `[1.0e-3, 1.0e-3, 1.0e-3]` | Pa s | Liquid water dynamic viscosity |
| `rhocpsolidm` | `SVector{3}` | `[3.3e6, 3.3e6, 3.0e6]` | $\text{J}/(\text{m}^3\cdot\text{K})$ | Volumetric solid heat capacity |
| `rhocpfluidm` | `SVector{3}` | `[1.0e6, 1.0e6, 3.0e6]` | $\text{J}/(\text{m}^3\cdot\text{K})$ | Volumetric fluid heat capacity |
| `alphasolidm` | `SVector{3}` | `[3.0e-5, 3.0e-5, 0.0]` | $1/\text{K}$ | Solid thermal expansion |
| `alphafluidm` | `SVector{3}` | `[5.0e-5, 5.0e-5, 0.0]` | $1/\text{K}$ | Fluid thermal expansion |
| `ksolidm` | `SVector{3}` | `[3.0, 3.0, 3000.0]` | $\text{W}/(\text{m}\cdot\text{K})$ | Solid thermal conductivity |
| `kfluidm` | `SVector{3}` | `[50.0, 50.0, 3000.0]` | $\text{W}/(\text{m}\cdot\text{K})$ | Fluid thermal conductivity |
| `gggsolidm` | `SVector{3}` | `[1.0e10, 1.0e10, 1.0e10]` | Pa | Solid shear elastic modulus |
| `frictsolidm` | `SVector{3}` | `[0.6, 0.6, 0.0]` | - | Internal friction coefficient |
| `cohessolidm` | `SVector{3}` | `[1.0e8, 1.0e8, 1.0e8]` | Pa | Cohesion |
| `tenssolidm` | `SVector{3}` | `[6.0e7, 6.0e7, 6.0e7]` | Pa | Tensile strength |
| `kphim0` | `SVector{3}` | `[1.0e-13, 1.0e-13, 1.0e-17]` | $\text{m}^2$ | Reference permeability |
| `tkm0` | `SVector{3}` | `[170.0, 170.0, 170.0]` | K | Initial temperature |

---

## `[output]`

| Parameter | Type | Default | Units | Description | Bounds |
|:---|:---|:---|:---|:---|:---|
| `output_dir` | `String` | `"output"` | - | Output directory path | Non-empty string |
| `savematstep` | `Int` | `10` | - | Checkpoint saving frequency | $\ge 1$ |
| `visstep` | `Int` | `1` | - | Visualization step cadence | $\ge 1$ |

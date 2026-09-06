# Silicate Rock Melting and Magma Rheology

This page documents and validates the silicate rock melting formulation, latent heat buffering, and magma rheology in `Erebus.jl`.

## Physical Formulation

### 1. Melt Fraction Model

The melt fraction $F_m(T, P)$ is computed as a linear function between pressure-dependent solidus $T_s(P)$ and liquidus $T_l(P)$ temperatures:

$$T_s(P) = T_{s,0} + \frac{dT}{dP} P, \quad T_l(P) = T_{l,0} + \frac{dT}{dP} P$$

$$F_m(T, P) = \begin{cases}
0, & T \le T_s(P) \\
\frac{T - T_s(P)}{T_l(P) - T_s(P)}, & T_s(P) < T < T_l(P) \\
1, & T \ge T_l(P)
\end{cases}$$

where $dT/dP$ is the Clapeyron slope (`dpdt_clapeyron`, default 0 Pa/K). For non-melting materials (such as sticky air), $F_m = 0$.

### 2. Latent Heat via Apparent Heat Capacity

Latent heat absorption and release during phase change is accounted for through an apparent volumetric heat capacity $c_{p,\text{eff}}$:

$$\rho c_{p,\text{eff}} = \rho c_p + \rho_s L_m \frac{\partial F_m}{\partial T}$$

where $L_m$ is the latent heat of silicate melting ($4.0 \times 10^5\text{ J/kg}$), $\rho_s$ is the solid silicate density, and $\partial F_m / \partial T = 1 / (T_l - T_s)$ in the melting interval. This formulation conserves enthalpy over the phase change:

$$\int_{T_s}^{T_l} (\rho c_{p,\text{eff}} - \rho c_p) \, dT = \rho_s L_m$$

following Gerya (2019, eq. 16.67).

### 3. Melt-Weakened Matrix Rheology

As silicate rocks melt, liquid melt lubricates crystal boundaries and weakens matrix shear viscosity:
- In the partially molten regime below critical melt fraction $\phi_{\text{crit}} = 0.4$, viscosity decreases exponentially with melt fraction:
  $$\eta(F_m) = \eta_{\text{solid}} \exp(-\alpha_\eta F_m)$$
  where $\alpha_\eta = 28$ is the rheological weakening coefficient (Costa et al., 2009; Gerya, 2019).
- Above the critical melt fraction $\phi_{\text{crit}} = 0.4$, the crystal framework breaks down into a crystal suspension. Viscosity decreases smoothly toward the liquid magma viscosity $\eta_{\text{melt}} = 10\text{ Pa}\cdot\text{s}$:
  $$\eta(F_m) = \eta(\phi_{\text{crit}}) \left( \frac{\eta_{\text{melt}}}{\eta(\phi_{\text{crit}})} \right)^{\frac{F_m - \phi_{\text{crit}}}{1 - \phi_{\text{crit}}}}$$
The resulting viscosity remains bounded within $[\eta_{\text{min}}, \eta_{\text{max}}]$, where $\eta_{\text{min}}$ ($10^{12}\text{ Pa}\cdot\text{s}$ by default) prevents ill-conditioning in the Stokes solver when $\eta_{\text{melt}} < \eta_{\text{min}}$.

### 4. Sub-Grid Soft Turbulence Model

Vigorous thermal convection in molten silicate magma ($Ra \sim 10^{16}\text{ to }10^{20}$) operates at spatial and temporal scales far below planetary grid resolution ($\Delta x \sim 1\text{ to }10\text{ km}$). When `soft_turbulence = true`, `Erebus.jl` parameterizes sub-grid convective heat transport by enhancing effective thermal conductivity following the soft turbulence scaling of Solomatov (2007):

$$k_{\text{turb}} = k_{\text{cond}} \left( \frac{\eta_{\text{num}}}{\eta_{\text{fluid}}} \right)^\beta$$

where $\eta_{\text{num}}$ is the numerical matrix viscosity, $\eta_{\text{fluid}}$ is the fluid silicate viscosity (`eta_fluid_silicate`, default $100\text{ Pa s}$), and $\beta$ is the scaling exponent (`turb_exponent`, default $1/3$ from Solomatov 2007; $\beta = 1/2$ corresponds to classical boundary layer scaling).

To eliminate unphysical conductivity jumps at marker state transitions, the formulation applies a regularized geometric blend in logarithmic space:

$$\log_{10} k_{\text{eff}} = (1 - w) \log_{10} k_{\text{cond}} + w \log_{10} k_{\text{turb}}$$

with composite transition weight $w = w_F \cdot w_T$:
- Melt fraction weighting via cubic smoothstep over window $[F_{\text{start}}, F_{\text{end}}]$ ($[0.30, 0.50]$ by default):
  $$\xi = \text{clamp}\left(\frac{F_m - F_{\text{start}}}{F_{\text{end}} - F_{\text{start}}}, 0, 1\right), \quad w_F = 3\xi^2 - 2\xi^3$$
- Thermal contrast weighting guarding against isothermal artifacts:
  $$w_T = \text{clamp}\left(\frac{T - T_{\text{surface}}}{dT_{\text{min}}}, 0, 1\right)$$

In the mush regime ($F_{\text{start}} \le F_m \le F_{\text{end}}$), the fluid reference viscosity blends smoothly from the current melt-weakened matrix viscosity down to the liquid silicate viscosity:

$$\log \eta_{\text{fluid}}(F_m) = (1 - \xi) \log \eta_{\text{solid}} + \xi \log \eta_{\text{fluid,silicate}}$$

The solver enforces strict lower and upper bounds:

$$k_{\text{eff}} = \max\left(k_{\text{cond}}, \text{clamp}(10^{\log_{10} k_{\text{eff}}}, k_{\text{floor}}, k_{\text{cutoff}})\right)$$

ensuring that $k_{\text{eff}} \ge k_{\text{cond}}$ everywhere with strict monotonicity throughout the melting interval.

## Analytical Verification

The implementation is verified against the following benchmarks:
1. **Stefan Moving-Boundary Problem**: Latent heat buffering matches the analytical similarity solution for 1D phase front propagation.
2. **Enthalpy Conservation**: Numerical integration of apparent heat capacity over the solidus-liquidus interval recovers the theoretical latent heat energy to floating-point precision ($< 10^{-12}$ relative error).
3. **Rheological Invariants**: Viscosity remains monotonic with $F_m$, strictly positive, and continuous at the transition $F_m = \phi_{\text{crit}}$.
4. **Soft Turbulence Monotonicity & Asymptotics**: Effective thermal conductivity matches $k_{\text{cond}}$ exactly for $F_m \le F_{\text{start}}$ or isothermal conditions, recovers $k_{\text{turb}}$ for $F_m \ge F_{\text{end}}$, and satisfies $k_{\text{eff}} \ge k_{\text{cond}}$ strictly monotonically over the transition window.
5. **Planetesimal Magma Ocean Cooling Benchmark**: Effective convective conductivity transports core radiogenic heat to the surface, buffering interior temperatures and preventing runaway super-liquidus overheating.

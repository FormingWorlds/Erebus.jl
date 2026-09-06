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

## Analytical Verification

The implementation is verified against the following benchmarks:
1. **Stefan Moving-Boundary Problem**: Latent heat buffering matches the analytical similarity solution for 1D phase front propagation.
2. **Enthalpy Conservation**: Numerical integration of apparent heat capacity over the solidus-liquidus interval recovers the theoretical latent heat energy to floating-point precision ($< 10^{-12}$ relative error).
3. **Rheological Invariants**: Viscosity remains monotonic with $F_m$, strictly positive, and continuous at the transition $F_m = \phi_{\text{crit}}$.

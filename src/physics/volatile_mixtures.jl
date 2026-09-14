"""
Compute the equilibrium pore fluid mixture freezing point under ammonia and solute depression.

Follows Croft et al. (1988) and Hogenboom et al. (1997) for continuous liquidus depression
from pure water freezing to the ammonia-water eutectic point at 176.0 K.

$(SIGNATURES)

# Parameters
- `X_nh3`: Ammonia mass or mole fraction in pore fluid [0, 1].
- `X_solute`: Dissolved salts or secondary solutes fraction [0, 1] (default: 0.0).

# Keywords
- `cfg`: Optional VolatileMixtureConfig struct supplying default eutectic and depression parameters.
- `T_freeze_pure`: Freezing temperature of pure water [K] (default: 273.15).
- `T_eutectic`: Eutectic temperature of ammonia-water system [K] (default: 176.0).
- `T_freeze_floor`: Absolute floor for liquid stability [K] (default: 176.0).
- `X_eutectic`: Eutectic ammonia composition fraction (default: 0.33).
- `lambda_nh3`: Liquidus depression slope for ammonia [K] (default: derived from eutectic depression (T_freeze_pure - T_eutectic) / X_eutectic ≈ 294.39 K).
- `lambda_solute`: Solute depression slope [K] (default: 50.0).

# Returns
- `Float64`: Equilibrium freezing point temperature [K].

# Raises
- `DomainError`: If volatile fractions violate bounds [0, 1] or are NaN.
"""
function compute_mixture_freezing_point(
    X_nh3::Real,
    X_solute::Real=0.0;
    cfg::Union{Nothing,VolatileMixtureConfig}=nothing,
    T_freeze_pure::Real=273.15,
    T_eutectic::Real=isnothing(cfg) ? 176.0 : cfg.T_eutectic_ammonia,
    T_freeze_floor::Real=isnothing(cfg) ? 176.0 : cfg.T_freeze_floor,
    X_eutectic::Real=0.33,
    lambda_nh3::Real=if isnothing(cfg)
        ((Float64(T_freeze_pure) - Float64(T_eutectic)) / Float64(X_eutectic))
    else
        cfg.lambda_nh3_depression
    end,
    lambda_solute::Real=isnothing(cfg) ? 50.0 : cfg.lambda_solute_depression,
)::Float64
    (isnan(X_nh3) || !(0.0 <= X_nh3 <= 1.0)) &&
        throw(DomainError(X_nh3, "X_nh3 must be in [0, 1]"))
    (isnan(X_solute) || !(0.0 <= X_solute <= 1.0)) &&
        throw(DomainError(X_solute, "X_solute must be in [0, 1]"))

    floor_limit = max(Float64(T_eutectic), Float64(T_freeze_floor))
    # Continuous liquidus depression toward the ammonia eutectic floor
    T_depressed =
        Float64(T_freeze_pure) - Float64(lambda_nh3) * Float64(X_nh3) -
        Float64(lambda_solute) * Float64(X_solute)
    return max(T_depressed, floor_limit)
end

"""
Compute the fluid mixture density as a function of temperature, pressure, and ammonia content.

$(SIGNATURES)

# Parameters
- `T`: Temperature [K].
- `P`: Pore fluid pressure [Pa].
- `X_nh3`: Ammonia fraction in fluid [0, 1].

# Keywords
- `rho_pure_ref`: Reference density of pure water [kg/m^3] (default: 1000.0).
- `T_ref`: Reference temperature [K] (default: 293.15).
- `P_ref`: Reference pressure [Pa] (default: 1.0e5).
- `alpha_th`: Volumetric thermal expansion coefficient [1/K] (default: 2.0e-4).
- `beta_comp`: Isothermal compressibility [1/Pa] (default: 4.0e-10).
- `coeff_nh3`: Density reduction coefficient for dissolved ammonia (default: 0.25).

# Returns
- `Float64`: Fluid mixture density [kg/m^3].

# Raises
- `DomainError`: If temperature or pressure is negative or NaN, or if X_nh3 is outside [0, 1].
"""
function compute_mixture_fluid_density(
    T::Real,
    P::Real,
    X_nh3::Real;
    rho_pure_ref::Real=1000.0,
    T_ref::Real=293.15,
    P_ref::Real=1.0e5,
    alpha_th::Real=2.0e-4,
    beta_comp::Real=4.0e-10,
    coeff_nh3::Real=0.25,
)::Float64
    (isnan(T) || T < 0.0) && throw(DomainError(T, "Temperature must be non-negative"))
    (isnan(P) || P < 0.0) && throw(DomainError(P, "Pressure must be non-negative"))
    (isnan(X_nh3) || !(0.0 <= X_nh3 <= 1.0)) &&
        throw(DomainError(X_nh3, "X_nh3 must be in [0, 1]"))

    rho_pure =
        Float64(rho_pure_ref) * (
            1.0 - Float64(alpha_th) * (Float64(T) - Float64(T_ref)) +
            Float64(beta_comp) * (Float64(P) - Float64(P_ref))
        )
    return rho_pure * (1.0 - Float64(coeff_nh3) * Float64(X_nh3))
end

"""
Compute pore fluid or ice viscosity across subfreezing and hydrothermal temperature regimes.

$(SIGNATURES)

# Parameters
- `T`: Temperature [K].
- `X_nh3`: Ammonia fraction in fluid [0, 1] (default: 0.0).

# Keywords
- `T_melt`: Optional melting point [K]. Computed via mixture depression if omitted.
- `eta_liquid_ref`: Reference liquid viscosity at reference temperature [Pa s] (default: 1.0e-3).
- `eta_ice`: Solid ice effective viscosity [Pa s] (default: 1.0e12).
- `E_act`: Activation energy for liquid viscous flow [J/mol] (default: 1.5e4).
- `R_gas`: Universal gas constant [J/mol/K] (default: 8.31446).
- `T_ref`: Reference temperature [K] (default: 293.15).

# Returns
- `Float64`: Dynamic viscosity of fluid or solid ice phase [Pa s].

# Raises
- `DomainError`: If temperature is negative or NaN, or if X_nh3 is outside [0, 1].
"""
function compute_mixture_fluid_viscosity(
    T::Real,
    X_nh3::Real=0.0;
    T_melt::Union{Nothing,Real}=nothing,
    eta_liquid_ref::Real=1.0e-3,
    eta_ice::Real=1.0e12,
    E_act::Real=1.5e4,
    R_gas::Real=8.31446,
    T_ref::Real=293.15,
)::Float64
    (isnan(T) || T < 0.0) && throw(DomainError(T, "Temperature must be non-negative"))
    (isnan(X_nh3) || !(0.0 <= X_nh3 <= 1.0)) &&
        throw(DomainError(X_nh3, "X_nh3 must be in [0, 1]"))

    Tm = T_melt !== nothing ? Float64(T_melt) : compute_mixture_freezing_point(X_nh3)
    if Float64(T) < Tm
        return Float64(eta_ice)
    end
    # Arrhenius temperature dependence for liquid mobile phase
    return Float64(eta_liquid_ref) *
           exp((Float64(E_act) / Float64(R_gas)) * (inv(Float64(T)) - inv(Float64(T_ref))))
end

"""
Evaluate thermal pyrolysis and devolatilization of refractory organic matter and structural hydrogen.

Simulates thermal breakdown of insoluble organic matter (IOM) and structural OH in
phyllosilicates / nominally anhydrous minerals across planetesimal thermal evolution.

$(SIGNATURES)

# Parameters
- `T`: Temperature [K].
- `C_refr`: Initial refractory carbon concentration [ppmw or mass fraction].
- `N_refr`: Initial refractory nitrogen concentration [ppmw or mass fraction].
- `H_refr`: Optional initial refractory/structural hydrogen concentration (default: 0.0).
- `cfg`: Refractory configuration struct.

# Keywords
- `f_graphite`: Fraction of pyrolyzed carbon retained as solid graphite residue (default: 0.60).
- `DeltaT_pyro`: Temperature interval for pyrolysis progression [K] (default: 250.0).
- `DeltaT_dehydrate`: Temperature interval for structural dehydration [K] (default: 100.0).

# Returns
- Named tuple with fields:
  - `C_refr_remaining`: Remaining unreacted refractory carbon.
  - `C_graphite_residue`: Solid graphite residue formed by pyrolysis.
  - `C_devolatilized_gas`: Carbon converted to volatile gas phase.
  - `N_refr_remaining`: Remaining refractory nitrogen.
  - `N_devolatilized_gas`: Nitrogen devolatilized to volatile gas phase.
  - `H_refr_remaining`: Remaining structural refractory hydrogen.
  - `H_dehydrated_gas`: Structural hydrogen released as volatile water/gas.

# Raises
- `DomainError`: If temperature, initial concentrations, or parameters are negative, NaN, or out of bounds.
"""
function evaluate_refractory_pyrolysis(
    T::Real,
    C_refr::Real,
    N_refr::Real,
    H_refr::Real=0.0,
    cfg::RefractoryConfig=RefractoryConfig();
    f_graphite::Real=0.60,
    DeltaT_pyro::Real=250.0,
    DeltaT_dehydrate::Real=100.0,
)
    (isnan(T) || T < 0.0) && throw(DomainError(T, "Temperature must be non-negative"))
    (isnan(C_refr) || C_refr < 0.0) &&
        throw(DomainError(C_refr, "C_refr must be non-negative"))
    (isnan(N_refr) || N_refr < 0.0) &&
        throw(DomainError(N_refr, "N_refr must be non-negative"))
    (isnan(H_refr) || H_refr < 0.0) &&
        throw(DomainError(H_refr, "H_refr must be non-negative"))
    (isnan(f_graphite) || !(0.0 <= f_graphite <= 1.0)) &&
        throw(DomainError(f_graphite, "f_graphite must be in [0, 1]"))
    (isnan(DeltaT_pyro) || DeltaT_pyro <= 0.0) &&
        throw(DomainError(DeltaT_pyro, "DeltaT_pyro must be positive"))
    (isnan(DeltaT_dehydrate) || DeltaT_dehydrate <= 0.0) &&
        throw(DomainError(DeltaT_dehydrate, "DeltaT_dehydrate must be positive"))

    if Float64(T) <= cfg.T_pyrolysis_C
        C_refr_remaining = Float64(C_refr)
        C_graphite_residue = 0.0
        C_devolatilized_gas = 0.0
        N_refr_remaining = Float64(N_refr)
        N_devolatilized_gas = 0.0
    else
        xi_pyro = min(0.80, (Float64(T) - cfg.T_pyrolysis_C) / Float64(DeltaT_pyro))
        Delta_C = xi_pyro * Float64(C_refr)
        C_graphite_residue = Float64(f_graphite) * Delta_C
        C_devolatilized_gas = (1.0 - Float64(f_graphite)) * Delta_C
        C_refr_remaining = Float64(C_refr) - Delta_C

        Delta_N = xi_pyro * Float64(N_refr)
        N_devolatilized_gas = Delta_N
        N_refr_remaining = Float64(N_refr) - Delta_N
    end

    if Float64(T) <= cfg.T_dehydrate_H || Float64(H_refr) == 0.0
        H_refr_remaining = Float64(H_refr)
        H_dehydrated_gas = 0.0
    else
        xi_deh = min(1.0, (Float64(T) - cfg.T_dehydrate_H) / Float64(DeltaT_dehydrate))
        H_dehydrated_gas = xi_deh * Float64(H_refr)
        H_refr_remaining = Float64(H_refr) - H_dehydrated_gas
    end

    return (;
        C_refr_remaining,
        C_graphite_residue,
        C_devolatilized_gas,
        N_refr_remaining,
        N_devolatilized_gas,
        H_refr_remaining,
        H_dehydrated_gas,
    )
end

function evaluate_refractory_pyrolysis(
    T::Real,
    C_refr::Real,
    N_refr::Real,
    cfg::RefractoryConfig;
    f_graphite::Real=0.60,
    DeltaT_pyro::Real=250.0,
    DeltaT_dehydrate::Real=100.0,
)
    return evaluate_refractory_pyrolysis(
        T,
        C_refr,
        N_refr,
        0.0,
        cfg;
        f_graphite=f_graphite,
        DeltaT_pyro=DeltaT_pyro,
        DeltaT_dehydrate=DeltaT_dehydrate,
    )
end

"""
Evaluate kinetic Arrhenius thermal pyrolysis and devolatilization of refractory C, N, and H.

Integrates first-order Arrhenius decomposition kinetics:
`dX_k/dt = -A_k * exp(-E_{a,k} / (R * T)) * X_k`
over timestep `dt` [s] at temperature `T` [K].

$(SIGNATURES)

# Arguments
- `T`: Temperature [K].
- `dt`: Timestep duration [s].
- `C_refr`: Current refractory carbon concentration [ppmw or mass fraction].
- `N_refr`: Current refractory nitrogen concentration [ppmw or mass fraction].
- `H_refr`: Current refractory/structural hydrogen concentration [ppmw or mass fraction].
- `cfg`: Refractory configuration struct containing Arrhenius rate parameters.

# Keyword Arguments
- `f_graphite`: Fraction of pyrolyzed carbon retained as solid graphite residue (default: `cfg.f_refr_C`).
- `R_gas`: Universal gas constant [J/(mol K)] (default: 8.314462618).

# Returns
- Named tuple with fields:
  - `C_refr_remaining`: Remaining unreacted refractory carbon.
  - `C_graphite_residue`: Solid graphite residue formed by pyrolysis.
  - `C_devolatilized_gas`: Carbon converted to volatile gas phase.
  - `N_refr_remaining`: Remaining refractory nitrogen.
  - `N_devolatilized_gas`: Nitrogen devolatilized to volatile gas phase.
  - `H_refr_remaining`: Remaining structural refractory hydrogen.
  - `H_dehydrated_gas`: Structural hydrogen released as volatile gas.
  - `dH_pyro_J_per_kg`: Specific enthalpy sink from endothermic pyrolysis [J/kg] (<= 0).

# Raises
- `DomainError`: If temperature, timestep, initial concentrations, or parameters are negative, NaN, or non-finite.
"""
function step_refractory_pyrolysis_kinetic(
    T::Real,
    dt::Real,
    C_refr::Real,
    N_refr::Real,
    H_refr::Real,
    cfg::RefractoryConfig=RefractoryConfig();
    f_graphite::Real=cfg.f_refr_C,
    R_gas::Real=8.314462618,
)
    (isnan(T) || T <= 0.0) &&
        throw(DomainError(T, "Temperature must be positive and finite"))
    (isnan(dt) || dt < 0.0) && throw(DomainError(dt, "dt must be non-negative and finite"))
    (isnan(C_refr) || C_refr < 0.0) &&
        throw(DomainError(C_refr, "C_refr must be non-negative"))
    (isnan(N_refr) || N_refr < 0.0) &&
        throw(DomainError(N_refr, "N_refr must be non-negative"))
    (isnan(H_refr) || H_refr < 0.0) &&
        throw(DomainError(H_refr, "H_refr must be non-negative"))
    (isnan(f_graphite) || !(0.0 <= f_graphite <= 1.0)) &&
        throw(DomainError(f_graphite, "f_graphite must be in [0, 1]"))

    T_val = Float64(T)
    dt_val = Float64(dt)
    fg = Float64(f_graphite)
    Rg = Float64(R_gas)

    # Carbon decomposition
    arg_C = cfg.Ea_C / (Rg * T_val)
    k_C = arg_C > 700.0 ? 0.0 : cfg.A_C * exp(-arg_C)
    factor_C = exp(-k_C * dt_val)
    C_refr_remaining = Float64(C_refr) * factor_C
    Delta_C = Float64(C_refr) - C_refr_remaining
    C_graphite_residue = fg * Delta_C
    C_devolatilized_gas = (1.0 - fg) * Delta_C

    # Nitrogen decomposition
    arg_N = cfg.Ea_N / (Rg * T_val)
    k_N = arg_N > 700.0 ? 0.0 : cfg.A_N * exp(-arg_N)
    factor_N = exp(-k_N * dt_val)
    N_refr_remaining = Float64(N_refr) * factor_N
    Delta_N = Float64(N_refr) - N_refr_remaining
    N_devolatilized_gas = Delta_N

    # Hydrogen decomposition
    arg_H = cfg.Ea_H / (Rg * T_val)
    k_H = arg_H > 700.0 ? 0.0 : cfg.A_H * exp(-arg_H)
    factor_H = exp(-k_H * dt_val)
    H_refr_remaining = Float64(H_refr) * factor_H
    Delta_H = Float64(H_refr) - H_refr_remaining
    H_dehydrated_gas = Delta_H

    # Endothermic latent heat sink (dH <= 0)
    dH_pyro_J_per_kg = -(
        Delta_C * cfg.dh_pyro_C + Delta_N * cfg.dh_pyro_N + Delta_H * cfg.dh_pyro_H
    )

    return (;
        C_refr_remaining,
        C_graphite_residue,
        C_devolatilized_gas,
        N_refr_remaining,
        N_devolatilized_gas,
        H_refr_remaining,
        H_dehydrated_gas,
        dH_pyro_J_per_kg,
    )
end

"""
Update refractory concentrations, porosity, and gas release from IOM pyrolysis across rock markers.

$(SIGNATURES)

# Arguments
- `tkm`: Marker temperature array [K].
- `dt`: Timestep duration [s].
- `phim`: Marker porosity array [-].
- `X_refr_C_m`: Marker refractory carbon array.
- `X_refr_N_m`: Marker refractory nitrogen array.
- `X_refr_H_m`: Marker refractory hydrogen array.
- `cfg`: RefractoryConfig.

# Keyword Arguments
- `xm`: Optional marker x-coordinates [m].
- `ym`: Optional marker y-coordinates [m].
- `coords`: Optional staggered grid coordinates.
- `DHP`: Optional P-node latent heating array [W/m^3] to apply localized endothermic heat sink.
- `rhosolid`: Solid density [kg/m^3] (default: 3000.0).
- `rhofluid`: Pore fluid density [kg/m^3] (default: 1000.0).
- `phimax`: Maximum porosity cap (default: 0.9999).
- `ppm_scale`: True if concentrations are in ppmw, false if mass fraction (default: false).

# Returns
- Named tuple `(; total_dC_gas, total_dN_gas, total_dH_gas, total_dC_graphite, total_dH_pyro)`
"""
function update_marker_pyrolysis!(
    tkm::AbstractVector{Float64},
    dt::Real,
    phim::AbstractVector{Float64},
    X_refr_C_m::AbstractVector{Float64},
    X_refr_N_m::AbstractVector{Float64},
    X_refr_H_m::AbstractVector{Float64},
    cfg::RefractoryConfig;
    xm::Union{Nothing,AbstractVector{Float64}}=nothing,
    ym::Union{Nothing,AbstractVector{Float64}}=nothing,
    coords=nothing,
    DHP::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    rhosolid::Real=3000.0,
    rhofluid::Real=1000.0,
    phimax::Real=0.9999,
    ppm_scale::Bool=false,
)
    marknum = length(tkm)
    rho_s = Float64(rhosolid)
    rho_f = max(Float64(rhofluid), 100.0)
    phi_max_val = Float64(phimax)
    dt_val = Float64(dt)
    scale = ppm_scale ? 1.0e-6 : 1.0

    tot_dC_gas = 0.0
    tot_dN_gas = 0.0
    tot_dH_gas = 0.0
    tot_dC_graphite = 0.0
    tot_dH_pyro = 0.0

    apply_dhp =
        DHP !== nothing &&
        xm !== nothing &&
        ym !== nothing &&
        coords !== nothing &&
        dt_val > 0.0
    DHP_pyro_sum = apply_dhp ? zeros(Float64, coords.Ny1, coords.Nx1) : nothing
    WT_pyro_sum = apply_dhp ? zeros(Float64, coords.Ny1, coords.Nx1) : nothing

    for m in 1:marknum
        T = tkm[m]
        T <= cfg.T_pyro_min && continue # Skip cold markers

        c_c = X_refr_C_m[m]
        c_n = X_refr_N_m[m]
        c_h = X_refr_H_m[m]
        (c_c == 0.0 && c_n == 0.0 && c_h == 0.0) && continue

        res = step_refractory_pyrolysis_kinetic(T, dt_val, c_c, c_n, c_h, cfg)

        X_refr_C_m[m] = res.C_refr_remaining
        X_refr_N_m[m] = res.N_refr_remaining
        X_refr_H_m[m] = res.H_refr_remaining

        d_c_gas = res.C_devolatilized_gas
        d_n_gas = res.N_devolatilized_gas
        d_h_gas = res.H_dehydrated_gas
        d_c_gr = res.C_graphite_residue

        tot_dC_gas += d_c_gas
        tot_dN_gas += d_n_gas
        tot_dH_gas += d_h_gas
        tot_dC_graphite += d_c_gr
        tot_dH_pyro += res.dH_pyro_J_per_kg

        dw_gas = (d_c_gas + d_n_gas + d_h_gas) * scale
        if dw_gas > 0.0
            dphi = dw_gas * (rho_s / rho_f)
            phim[m] = min(phi_max_val, phim[m] + dphi)
        end

        if apply_dhp && res.dH_pyro_J_per_kg < 0.0
            q_pyro = rho_s * (res.dH_pyro_J_per_kg / dt_val)
            i, j, weights = fix_weights(
                xm[m],
                ym[m],
                coords.xp,
                coords.yp,
                coords.dx,
                coords.dy,
                coords.jmin_p,
                coords.jmax_p,
                coords.imin_p,
                coords.imax_p,
            )
            interpolate_add_to_grid!(i, j, weights, q_pyro, DHP_pyro_sum)
            interpolate_add_to_grid!(i, j, weights, one(1.0), WT_pyro_sum)
        end
    end

    if apply_dhp
        for j in 1:(coords.Nx1), i in 1:(coords.Ny1)
            if WT_pyro_sum[i, j] > 0.0
                DHP[i, j] += DHP_pyro_sum[i, j] / WT_pyro_sum[i, j]
            end
        end
    end

    return (;
        total_dC_gas=tot_dC_gas,
        total_dN_gas=tot_dN_gas,
        total_dH_gas=tot_dH_gas,
        total_dC_graphite=tot_dC_graphite,
        total_dH_pyro=tot_dH_pyro,
    )
end

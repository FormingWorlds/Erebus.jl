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

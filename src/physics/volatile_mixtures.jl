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
- `R_gas`: Universal gas constant [J/mol/K] (default: R_GAS).
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
    R_gas::Real=R_GAS,
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
- `R_gas`: Universal gas constant [J/(mol K)] (default: R_GAS).

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
    R_gas::Real=R_GAS,
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
    speciate_pyrolysis_carbon_redox!(
        redox_props, redox_cfg, m::Int, T::Float64, scale::Float64,
        d_c_gas::Float64, d_c_gr::Float64, d_h_gas::Float64, d_n_gas::Float64, Xfem
    )

Speciate devolatilized carbon into CO, CO2, and CH4 based on local redox conditions and rock buffer capacity.
"""
function speciate_pyrolysis_carbon_redox!(
    redox_props,
    redox_cfg,
    m::Int,
    T::Float64,
    scale::Float64,
    d_c_gas::Float64,
    d_c_gr::Float64,
    d_h_gas::Float64,
    d_n_gas::Float64,
    Xfem;
    Xfe_bulk::Union{Nothing,AbstractVector{Float64}}=nothing,
    rho_s::Real=3000.0,
    rho_metal_val::Real=7000.0,
    phi_pack::Real=0.65,
)
    nC_gr_m = hasproperty(redox_props, :nC_graphite_m) ? redox_props.nC_graphite_m : nothing
    nCO_m = hasproperty(redox_props, :nCO_m) ? redox_props.nCO_m : nothing
    nCO2_m = hasproperty(redox_props, :nCO2_m) ? redox_props.nCO2_m : nothing
    nCH4_m = hasproperty(redox_props, :nCH4_m) ? redox_props.nCH4_m : nothing

    dn_c_gr = (d_c_gr * scale) / M_C
    dn_c_gas = (d_c_gas * scale) / M_C
    dn_c_iom = dn_c_gr + dn_c_gas

    def_co2 = d_c_gas * (M_CO2 / M_C)
    def_h2 = d_h_gas
    if nC_gr_m === nothing || dn_c_iom <= 0.0
        return (
            d_c_gas=d_c_gas,
            d_c_gr=d_c_gr,
            m_co=0.0,
            m_co2=def_co2,
            m_ch4=0.0,
            m_h2=def_h2,
            m_h2o=0.0,
            m_n2=d_n_gas,
        )
    end

    nFe0_m = hasproperty(redox_props, :nFe0_m) ? redox_props.nFe0_m : nothing
    nFe2_m = hasproperty(redox_props, :nFe2_m) ? redox_props.nFe2_m : nothing
    nFe3_m = hasproperty(redox_props, :nFe3_m) ? redox_props.nFe3_m : nothing

    if nFe0_m === nothing || nFe2_m === nothing || nFe3_m === nothing
        nC_gr_m[m] += dn_c_gr
        if nCO2_m !== nothing
            nCO2_m[m] += dn_c_gas
        end
        return (
            d_c_gas=d_c_gas,
            d_c_gr=d_c_gr,
            m_co=0.0,
            m_co2=def_co2,
            m_ch4=0.0,
            m_h2=def_h2,
            m_h2o=0.0,
            m_n2=d_n_gas,
        )
    end

    c_cur = marker_redox_components(
        nFe0_m[m],
        nFe2_m[m],
        nFe3_m[m];
        n_C_graphite=nC_gr_m[m],
        n_CO=nCO_m !== nothing ? nCO_m[m] : 0.0,
        n_CO2=nCO2_m !== nothing ? nCO2_m[m] : 0.0,
        n_CH4=nCH4_m !== nothing ? nCH4_m[m] : 0.0,
    )

    d_IW = if (hasproperty(redox_props, :deltaIW_m) && redox_props.deltaIW_m !== nothing)
        redox_props.deltaIW_m[m]
    else
        0.0
    end
    T_safe = max(100.0, T)
    log10_fo2 = compute_iron_wustite_fO2(T_safe; delta_IW=d_IW)

    logK_CO2 = 14800.0 / T_safe - 4.58
    r_CO2 = 10.0^clamp(logK_CO2 + 0.5 * log10_fo2, -50.0, 50.0)
    r_CH4 = if (T_safe < 900.0 && d_IW < 0.0)
        10.0^clamp((900.0 - T_safe) / 200.0 - 0.5 * (d_IW + 1.0), -50.0, 50.0)
    else
        0.0
    end

    denom = 1.0 + r_CO2 + r_CH4
    f_CO = 1.0 / denom
    f_CO2 = r_CO2 / denom
    f_CH4 = r_CH4 / denom

    dn_co = dn_c_gas * f_CO
    dn_co2 = dn_c_gas * f_CO2
    dn_ch4 = dn_c_gas * f_CH4

    dn_h_avail = (d_h_gas * scale) / M_H + 2.0 * c_cur.n_H2
    max_ch4_from_h = dn_h_avail / 4.0
    if dn_ch4 > max_ch4_from_h
        dn_ch4 = max(0.0, max_ch4_from_h)
        dn_rem = dn_c_gas - dn_ch4
        denom_co = 1.0 + r_CO2
        dn_co = dn_rem / denom_co
        dn_co2 = dn_rem * r_CO2 / denom_co
    end

    delta_rb_req = 2.0 * dn_co + 4.0 * dn_co2 - 4.0 * dn_ch4
    if delta_rb_req > 0.0
        avail_ox = 3.0 * c_cur.n_Fe3 + 2.0 * c_cur.n_Fe2 + 2.0 * c_cur.n_H2O
        if delta_rb_req > avail_ox
            scale_ox = max(0.0, avail_ox / delta_rb_req)
            unoxidized_gas = dn_c_gas * (1.0 - scale_ox)
            dn_c_gr += unoxidized_gas
            dn_co *= scale_ox
            dn_co2 *= scale_ox
            dn_ch4 *= scale_ox

            unox_mass = (unoxidized_gas * M_C) / scale
            d_c_gas = max(0.0, d_c_gas - unox_mass)
            d_c_gr += unox_mass
        end
    elseif delta_rb_req < 0.0
        avail_red = 3.0 * c_cur.n_Fe0 + 1.0 * c_cur.n_Fe2 + 2.0 * c_cur.n_H2
        req_red = -delta_rb_req
        if req_red > avail_red
            scale_red = max(0.0, avail_red / req_red)
            unreduced_gas = dn_c_gas * (1.0 - scale_red)
            dn_c_gr += unreduced_gas
            dn_co *= scale_red
            dn_co2 *= scale_red
            dn_ch4 *= scale_red

            unred_mass = (unreduced_gas * M_C) / scale
            d_c_gas = max(0.0, d_c_gas - unred_mass)
            d_c_gr += unred_mass
        end
    end

    m_co = (dn_co * M_CO) / scale
    m_co2 = (dn_co2 * M_CO2) / scale
    m_ch4 = (dn_ch4 * M_CH4) / scale
    h_in_ch4 = 4.0 * (dn_ch4 * M_H) / scale
    m_h2 = max(0.0, d_h_gas - h_in_ch4)
    m_h2o = 0.0
    m_n2 = d_n_gas

    ref_sym = redox_cfg !== nothing ? redox_cfg.reference : :mantle
    c_up = pyrolyze_redox_budget(
        c_cur,
        dn_c_iom,
        dn_c_gr,
        dn_co,
        dn_co2,
        dn_ch4;
        auto_balance=true,
        reference=ref_sym,
    )

    dn_fe0_smelted = max(0.0, c_up.n_Fe0 - c_cur.n_Fe0)
    if dn_fe0_smelted > 0.0
        dw_fe0 = dn_fe0_smelted * M_Fe
        phi_pack_val = Float64(phi_pack)
        if Xfe_bulk !== nothing
            w_cur = metal_volume_to_mass_fraction(Xfe_bulk[m], rho_metal_val, rho_s)
            w_new = min(1.0, w_cur + dw_fe0)
            phi_new = metal_mass_to_volume_fraction(w_new, rho_metal_val, rho_s)
            phi_target = max(Xfe_bulk[m], min(phi_pack_val, phi_new))
            dphi_fe0 = max(0.0, phi_target - Xfe_bulk[m])
            Xfe_bulk[m] = phi_target
            F_fe_local = compute_metal_melt_fraction(T)
            if F_fe_local > 0.0 && Xfem !== nothing
                Xfem[m] = max(Xfem[m], min(phi_pack_val, Xfem[m] + dphi_fe0 * F_fe_local))
            end
        elseif Xfem !== nothing
            w_cur = metal_volume_to_mass_fraction(Xfem[m], rho_metal_val, rho_s)
            w_new = min(1.0, w_cur + dw_fe0)
            phi_new = metal_mass_to_volume_fraction(w_new, rho_metal_val, rho_s)
            Xfem[m] = max(Xfem[m], min(phi_pack_val, phi_new))
        end
    end

    nFe0_m[m] = c_up.n_Fe0
    nFe2_m[m] = c_up.n_Fe2
    nFe3_m[m] = c_up.n_Fe3
    nC_gr_m[m] = c_up.n_C_graphite
    if nCO_m !== nothing
        nCO_m[m] = c_up.n_CO
    end
    if nCO2_m !== nothing
        nCO2_m[m] = c_up.n_CO2
    end
    if nCH4_m !== nothing
        nCH4_m[m] = c_up.n_CH4
    end

    return (
        d_c_gas=d_c_gas,
        d_c_gr=d_c_gr,
        m_co=m_co,
        m_co2=m_co2,
        m_ch4=m_ch4,
        m_h2=m_h2,
        m_h2o=m_h2o,
        m_n2=m_n2,
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
- `redox_props`: Optional NamedTuple with marker redox arrays to deposit graphite and gas products.
- `redox_cfg`: Optional RedoxConfig to govern pyrolysis redox coupling.
- `Xfem`: Optional marker molten metallic iron volume fraction array.
- `Xfe_bulk`: Optional marker bulk metallic iron volume fraction array.
- `phi_pack`: Maximum metal volume fraction packing limit (default: 0.65).

# Returns
- Named tuple `(; total_dC_gas, total_dN_gas, total_dH_gas, total_dC_graphite, total_dH_pyro)`

# Notes
- Local heat sink `DHP` captures macromolecular organic decomposition endothermicity;
  secondary high-temperature gas-phase redox recombination enthalpies are neglected.
"""
function update_marker_pyrolysis!(
    tkm::AbstractVector{Float64},
    dt::Real,
    phim::AbstractVector{Float64},
    X_refr_C_m::AbstractVector{Float64},
    X_refr_N_m::AbstractVector{Float64},
    X_refr_H_m::AbstractVector{Float64},
    cfg::RefractoryConfig;
    tm::Union{Nothing,AbstractVector{<:Integer}}=nothing,
    rhosolidm::Union{Nothing,AbstractVector{<:Real}}=nothing,
    xm::Union{Nothing,AbstractVector{Float64}}=nothing,
    ym::Union{Nothing,AbstractVector{Float64}}=nothing,
    coords::GridCoordinates=default_grid_coordinates(),
    DHP::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    rhosolid::Real=3000.0,
    rhofluid::Real=1000.0,
    phimax::Real=0.9999,
    ppm_scale::Bool=false,
    redox_props=nothing,
    redox_cfg=nothing,
    Xfem::Union{Nothing,AbstractVector{Float64}}=nothing,
    Xfe_bulk::Union{Nothing,AbstractVector{Float64}}=nothing,
    rho_metal::Real=7000.0,
    phi_pack::Real=0.65,
)
    marknum = length(tkm)
    Xfe_bulk !== nothing &&
        length(Xfe_bulk) < marknum &&
        throw(DimensionMismatch("length(Xfe_bulk) must be >= marknum"))
    Xfem !== nothing &&
        length(Xfem) < marknum &&
        throw(DimensionMismatch("length(Xfem) must be >= marknum"))

    rho_f = max(Float64(rhofluid), 100.0)
    phi_max_val = Float64(phimax)
    dt_val = Float64(dt)
    scale = ppm_scale ? 1.0e-6 : 1.0

    tot_dC_gas = 0.0
    tot_dN_gas = 0.0
    tot_dH_gas = 0.0
    tot_dCO_gas = 0.0
    tot_dCO2_gas = 0.0
    tot_dCH4_gas = 0.0
    tot_dH2_gas = 0.0
    tot_dH2O_gas = 0.0
    tot_dN2_gas = 0.0
    tot_dC_graphite = 0.0
    tot_dH_pyro = 0.0

    apply_dhp = DHP !== nothing && xm !== nothing && ym !== nothing && dt_val > 0.0
    if apply_dhp && size(DHP) != (coords.Ny1, coords.Nx1)
        throw(
            DimensionMismatch(
                "DHP size $(size(DHP)) must match grid dimensions ($(coords.Ny1), $(coords.Nx1))",
            ),
        )
    end
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

        rho_s = if tm !== nothing && rhosolidm !== nothing && m <= length(tm)
            m_tm = tm[m]
            if (m_tm >= 1 && m_tm <= length(rhosolidm))
                Float64(rhosolidm[m_tm])
            else
                (length(rhosolidm) >= 1 ? Float64(rhosolidm[1]) : Float64(rhosolid))
            end
        elseif rhosolidm !== nothing && length(rhosolidm) >= 1
            Float64(rhosolidm[1])
        else
            Float64(rhosolid)
        end

        if redox_props !== nothing &&
            (redox_cfg === nothing || (redox_cfg.active && redox_cfg.pyrolysis_redox))
            spec = speciate_pyrolysis_carbon_redox!(
                redox_props,
                redox_cfg,
                m,
                T,
                scale,
                d_c_gas,
                d_c_gr,
                d_h_gas,
                d_n_gas,
                Xfem;
                Xfe_bulk=Xfe_bulk,
                rho_s=rho_s,
                rho_metal_val=rho_metal,
                phi_pack=phi_pack,
            )
            d_c_gas = spec.d_c_gas
            d_c_gr = spec.d_c_gr
            m_co = spec.m_co
            m_co2 = spec.m_co2
            m_ch4 = spec.m_ch4
            m_h2 = spec.m_h2
            m_h2o = spec.m_h2o
            m_n2 = spec.m_n2
        else
            m_co = 0.0
            m_co2 = d_c_gas * (M_CO2 / M_C)
            m_ch4 = 0.0
            m_h2 = d_h_gas
            m_h2o = 0.0
            m_n2 = d_n_gas
        end

        tot_dC_gas += d_c_gas
        tot_dN_gas += d_n_gas
        tot_dH_gas += d_h_gas
        tot_dCO_gas += m_co
        tot_dCO2_gas += m_co2
        tot_dCH4_gas += m_ch4
        tot_dH2_gas += m_h2
        tot_dH2O_gas += m_h2o
        tot_dN2_gas += m_n2
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
        total_dCO_gas=tot_dCO_gas,
        total_dCO2_gas=tot_dCO2_gas,
        total_dCH4_gas=tot_dCH4_gas,
        total_dH2_gas=tot_dH2_gas,
        total_dH2O_gas=tot_dH2O_gas,
        total_dN2_gas=tot_dN2_gas,
        total_dC_graphite=tot_dC_graphite,
        total_dH_pyro=tot_dH_pyro,
    )
end

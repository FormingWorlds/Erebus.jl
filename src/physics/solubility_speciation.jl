"""
Compute oxygen fugacity of the iron-wüstite (IW) buffer.

$(SIGNATURES)

Calculates log10(fO2 [bar]) using the empirical 1-bar parameterization (e.g. O'Neill 1988; Campbell et al. 2009):
    log10(fO2) = 6.541 - 28164 / T + ΔIW

# Arguments
- `T_K`: Temperature [K]

# Keyword Arguments
- `delta_IW`: Oxygen fugacity offset relative to IW buffer in log10 units (default: 0.0)

# Notes
The default `delta_IW = 0.0` represents the neutral iron-wüstite buffer. Planetesimal interiors are typically more reduced, for example `delta_IW = -1.0` in `VolatilesConfig`.

# Returns
- `log10_fO2`: log10 of oxygen fugacity in bar
"""
function compute_iron_wustite_fO2(T_K::Real; delta_IW::Real=0.0)::Float64
    T_val = Float64(T_K)
    if T_val <= 0.0 || !isfinite(T_val)
        throw(DomainError(T_val, "Temperature must be > 0 and finite"))
    end
    d_IW = Float64(delta_IW)
    if !isfinite(d_IW)
        throw(DomainError(d_IW, "delta_IW must be finite"))
    end
    return 6.541 - 28164.0 / T_val + d_IW
end

"""
Compute equilibrium dissolved water solubility in silicate melt at low pressure.

$(SIGNATURES)

Follows the low-pressure square-root law (Burnham 1979; Dixon et al. 1995; Sossi et al. 2023)
where water dissolves dominantly as hydroxyl (OH⁻):
- `:burnham_dixon`: Burnham (1979) / Dixon et al. (1995) baseline:
    w_H2O = As * sqrt(max(0, P [MPa]))  [wt%]
- `:sossi_peridotite`: Sossi et al. (2023) peridotitic melt:
    w_H2O = 524.0 * sqrt(max(0, P [bar]))  [ppmw] (converted to wt%)
- `:basalt_dixon`: Dixon et al. (1995) MORB basalt:
    w_H2O = 965.0 * sqrt(max(0, P [bar]))  [ppmw] (converted to wt%)
- `:newcombe_lunar`: Newcombe et al. (2017) lunar glass:
    w_H2O = 683.0 * sqrt(max(0, P [bar]))  [ppmw] (converted to wt%)

# Arguments
- `P_Pa`: Pore fluid pressure [Pa]

# Keyword Arguments
- `As`: Water solubility coefficient [wt% / MPa^0.5] for `:burnham_dixon` (default: 0.40)
- `law`: Solubility formulation (`:burnham_dixon`, `:sossi_peridotite`, `:basalt_dixon`, `:newcombe_lunar`)

# Returns
- `w_H2O`: Equilibrium dissolved water concentration in melt [wt%]
"""
function compute_water_solubility_melt(
    P_Pa::Real; As::Real=0.40, law::Symbol=:burnham_dixon
)::Float64
    P_val = Float64(P_Pa)
    if !isfinite(P_val)
        throw(DomainError(P_val, "Pressure must be finite"))
    end
    if law === :burnham_dixon
        As_val = Float64(As)
        if As_val <= 0.0 || !isfinite(As_val)
            throw(
                DomainError(
                    As_val, "Water solubility coefficient As must be > 0 and finite"
                ),
            )
        end
        if P_val <= 0.0
            return 0.0
        end
        return As_val * sqrt(P_val * 1.0e-6)
    elseif law === :sossi_peridotite
        if P_val <= 0.0
            return 0.0
        end
        p_bar = P_val * 1.0e-5
        return (524.0 * sqrt(p_bar)) * 1.0e-4
    elseif law === :basalt_dixon
        if P_val <= 0.0
            return 0.0
        end
        p_bar = P_val * 1.0e-5
        return (965.0 * sqrt(p_bar)) * 1.0e-4
    elseif law === :newcombe_lunar
        if P_val <= 0.0
            return 0.0
        end
        p_bar = P_val * 1.0e-5
        return (683.0 * sqrt(p_bar)) * 1.0e-4
    else
        throw(ArgumentError("Unknown water solubility law: $law"))
    end
end

"""
Compute equilibrium nitrogen solubility in silicate melt under reducing conditions.

$(SIGNATURES)

Partitions nitrogen into physical molecular dissolution (N2) and chemical nitride dissolution (N³⁻)
following Libourel et al. (2003) and Boulliung et al. (2020):
    w_phys = Kh * f_N2  [ppm]
    w_chem = (C_nitride * 10^4) * sqrt(f_N2) * 10^(-0.75 * ΔIW)  [ppm]
    w_total = w_phys + w_chem  [ppm]

This parameterization is isothermal at reference magmatic temperature (~1673 K).

# Arguments
- `P_Pa`: Pore fluid pressure [Pa]
- `delta_IW`: Oxygen fugacity offset relative to IW buffer [log10 units] (default: 0.0)

# Keyword Arguments
- `Kh`: Henry law coefficient for molecular N2 [ppm / bar] (default: 0.40)
- `C_nitride`: Chemical nitride capacity [wt% / bar^0.5] (default: 1.0e-3)

# Notes
The default `delta_IW = 0.0` corresponds to the neutral iron-wüstite buffer. Planetesimal interiors are typically more reduced, for example `delta_IW = -1.0` in `VolatilesConfig`.

# Returns
- `NamedTuple`: `(; total_ppm, physical_ppm, chemical_ppm)`
"""
function compute_nitrogen_solubility_melt(
    P_Pa::Real, delta_IW::Real; Kh::Real=0.40, C_nitride::Real=1.0e-3
)::@NamedTuple{total_ppm::Float64, physical_ppm::Float64, chemical_ppm::Float64}
    P_val = Float64(P_Pa)
    if !isfinite(P_val)
        throw(DomainError(P_val, "Pressure must be finite"))
    end
    d_IW = Float64(delta_IW)
    if !isfinite(d_IW) || abs(d_IW) > 50.0
        throw(DomainError(d_IW, "delta_IW must be finite and within [-50, 50]"))
    end
    Kh_val = Float64(Kh)
    if Kh_val <= 0.0 || !isfinite(Kh_val)
        throw(DomainError(Kh_val, "Henry coefficient Kh must be > 0 and finite"))
    end
    Cn_val = Float64(C_nitride)
    if Cn_val <= 0.0 || !isfinite(Cn_val)
        throw(DomainError(Cn_val, "Nitride capacity C_nitride must be > 0 and finite"))
    end

    if P_val <= 0.0
        return (total_ppm=0.0, physical_ppm=0.0, chemical_ppm=0.0)
    end

    # Pore fluid pressure converted to bar for gas fugacity
    f_N2 = P_val * 1.0e-5
    physical_ppm = Kh_val * f_N2

    # Chemical nitride scaling relative to iron-wüstite buffer
    fO2_ratio = 10.0^d_IW
    chemical_ppm = (Cn_val * 1.0e4) * sqrt(f_N2) * (fO2_ratio)^(-0.75)
    total_ppm = physical_ppm + chemical_ppm

    return (total_ppm=total_ppm, physical_ppm=physical_ppm, chemical_ppm=chemical_ppm)
end

function compute_nitrogen_solubility_melt(
    P_Pa::Real; delta_IW::Real=0.0, Kh::Real=0.40, C_nitride::Real=1.0e-3
)::@NamedTuple{total_ppm::Float64, physical_ppm::Float64, chemical_ppm::Float64}
    return compute_nitrogen_solubility_melt(P_Pa, delta_IW; Kh=Kh, C_nitride=C_nitride)
end

"""
Compute devolatilization yield of primordial organic nitrogen as a function of temperature.

$(SIGNATURES)

Models thermal decomposition of organic nitrogen matter via a logistic sigmoid:
    yield = inv(1 + exp(-(T - T_devol) / ΔT))

# Arguments
- `T_K`: Temperature [K]

# Keyword Arguments
- `T_devol`: Characteristic devolatilization midpoint temperature [K] (default: 550.0)
- `delta_T`: Transition temperature scale [K] (default: 50.0)

# Returns
- `yield`: Devolatilized nitrogen fraction in [0, 1]
"""
function compute_organic_nitrogen_yield(
    T_K::Real; T_devol::Real=550.0, delta_T::Real=50.0
)::Float64
    T_val = Float64(T_K)
    if T_val <= 0.0 || !isfinite(T_val)
        throw(DomainError(T_val, "Temperature must be > 0 and finite"))
    end
    Td = Float64(T_devol)
    if Td <= 0.0 || !isfinite(Td)
        throw(DomainError(Td, "T_devol must be > 0 and finite"))
    end
    dT = Float64(delta_T)
    if dT <= 0.0 || !isfinite(dT)
        throw(DomainError(dT, "delta_T must be > 0 and finite"))
    end

    arg = (T_val - Td) / dT
    # Clamp argument to prevent numerical underflow/overflow in exp
    if arg > 40.0
        return 1.0
    elseif arg < -40.0
        return 0.0
    end
    return inv(1.0 + exp(-arg))
end

"""
Compute equilibrium dissolved molecular hydrogen (H2) solubility in silicate melt.

$(SIGNATURES)

Calculates dissolved H2 concentration under reducing magmatic conditions:
- `:hirschmann2012`: Hirschmann et al. (2012) synthetic basalt fit:
    log10(X_H2 [ppmw]) = 1.1008 + 0.5241 * log10(p_H2 [bar])
- `:gaillard2003`: Gaillard et al. (2003) power law fit:
    X_H2 [ppmw] = 0.163 * (p_H2 [bar])^1.252

# Arguments
- `p_H2_Pa`: Partial pressure of H2 [Pa]

# Keyword Arguments
- `law`: Formulation (`:hirschmann2012` or `:gaillard2003`)

# Returns
- `ppmw`: Dissolved H2 concentration in melt [ppmw]
"""
function compute_h2_solubility_melt(p_H2_Pa::Real; law::Symbol=:hirschmann2012)::Float64
    p_val = Float64(p_H2_Pa)
    if !isfinite(p_val)
        throw(DomainError(p_val, "Partial pressure of H2 must be finite"))
    end
    if p_val <= 0.0
        return 0.0
    end
    p_bar = p_val * 1.0e-5
    if law === :hirschmann2012
        return 10.0^(1.10083602 + 0.52413928 * log10(p_bar))
    elseif law === :gaillard2003
        return 0.163 * (p_bar^1.252)
    else
        throw(ArgumentError("Unknown H2 solubility law: $law"))
    end
end

"""
Compute equilibrium nitrogen solubility in silicate melt using Dasgupta et al. (2022).

$(SIGNATURES)

Partitions nitrogen into physical molecular dissolution (N2) and chemical nitride dissolution (N3-)
incorporating temperature, total pressure, redox state, and melt composition:
    w_chem [ppmw] = sqrt(p_N2 [GPa]) * exp(5908.0 * sqrt(p_tot [GPa]) / T - 1.6 * ΔIW)
    w_phys [ppmw] = p_N2 [GPa] * exp(4.67 + 7.11 * x_SiO2 - 13.06 * x_Al2O3 - 120.67 * x_TiO2)
    w_total = w_chem + w_phys

# Arguments
- `p_N2_Pa`: Partial pressure of N2 [Pa]
- `p_total_Pa`: Total pressure [Pa]
- `T_K`: Melt temperature [K]
- `delta_IW`: Oxygen fugacity offset relative to IW buffer [log10 units]

# Keyword Arguments
- `x_SiO2`: Silicate melt SiO2 mole fraction (default: 0.56, Earth/chondritic mantle)
- `x_Al2O3`: Silicate melt Al2O3 mole fraction (default: 0.11)
- `x_TiO2`: Silicate melt TiO2 mole fraction (default: 0.01)

# Returns
- `NamedTuple`: `(; total_ppm, physical_ppm, chemical_ppm)`
"""
function compute_nitrogen_solubility_dasgupta(
    p_N2_Pa::Real,
    p_total_Pa::Real,
    T_K::Real,
    delta_IW::Real;
    x_SiO2::Real=0.56,
    x_Al2O3::Real=0.11,
    x_TiO2::Real=0.01,
)::@NamedTuple{total_ppm::Float64, physical_ppm::Float64, chemical_ppm::Float64}
    p_val = Float64(p_N2_Pa)
    if !isfinite(p_val)
        throw(DomainError(p_val, "p_N2 must be finite"))
    end
    p_tot = Float64(p_total_Pa)
    if !isfinite(p_tot)
        throw(DomainError(p_tot, "p_total must be finite"))
    end
    T = Float64(T_K)
    if T <= 0.0 || !isfinite(T)
        throw(DomainError(T, "Temperature must be > 0 and finite"))
    end
    d_IW = Float64(delta_IW)
    if !isfinite(d_IW) || abs(d_IW) > 50.0
        throw(DomainError(d_IW, "delta_IW must be finite and within [-50, 50]"))
    end
    for (nm, v) in (("x_SiO2", x_SiO2), ("x_Al2O3", x_Al2O3), ("x_TiO2", x_TiO2))
        fv = Float64(v)
        if fv < 0.0 || fv > 1.0 || !isfinite(fv)
            throw(DomainError(fv, "$nm must be in [0, 1] and finite"))
        end
    end

    if p_val <= 0.0
        return (total_ppm=0.0, physical_ppm=0.0, chemical_ppm=0.0)
    end

    pN2_GPa = p_val * 1.0e-9
    ptot_GPa = max(p_tot, 0.0) * 1.0e-9

    chem_exp = (5908.0 * sqrt(max(ptot_GPa, 1.0e-15))) / T - 1.6 * d_IW
    chem_exp_clamped = clamp(chem_exp, -100.0, 100.0)
    chemical_ppm = sqrt(pN2_GPa) * exp(chem_exp_clamped)

    phys_prefactor = exp(
        4.67 + 7.11 * Float64(x_SiO2) - 13.06 * Float64(x_Al2O3) - 120.67 * Float64(x_TiO2)
    )
    physical_ppm = pN2_GPa * phys_prefactor

    total_ppm = physical_ppm + chemical_ppm
    return (total_ppm=total_ppm, physical_ppm=physical_ppm, chemical_ppm=chemical_ppm)
end

"""
Compute equilibrium dissolved carbon monoxide (CO) in silicate melt.

$(SIGNATURES)

- `:armstrong2015`: Armstrong et al. (2015) mafic melt:
    log10(X_CO [ppmw]) = -0.738 + 0.876 * log10(p_CO [bar]) - 5.44e-5 * p_tot [bar]
- `:yoshioka2019_morb`: Yoshioka et al. (2019) MORB basalt at graphite saturation:
    X_C [wt%] = 10^(-5.20 + 0.80 * log10(p_CO [bar])) -> X_CO [ppmw] = X_C * 1e4 * (28.0101 / 12.011)

# Arguments
- `p_CO_Pa`: Partial pressure of CO [Pa]
- `p_total_Pa`: Total pressure [Pa]

# Keyword Arguments
- `law`: Formulation (`:armstrong2015` or `:yoshioka2019_morb`)

# Returns
- `ppmw`: Dissolved CO concentration in melt [ppmw]
"""
function compute_co_solubility_melt(
    p_CO_Pa::Real, p_total_Pa::Real; law::Symbol=:armstrong2015
)::Float64
    p_co = Float64(p_CO_Pa)
    if !isfinite(p_co)
        throw(DomainError(p_co, "p_CO must be finite"))
    end
    p_tot = Float64(p_total_Pa)
    if !isfinite(p_tot)
        throw(DomainError(p_tot, "p_total must be finite"))
    end
    if p_co <= 0.0
        return 0.0
    end
    p_co_bar = p_co * 1.0e-5
    p_tot_bar = max(p_tot, 0.0) * 1.0e-5

    if law === :armstrong2015
        log_co = -0.738 + 0.876 * log10(p_co_bar) - 5.44e-5 * p_tot_bar
        return 10.0^log_co
    elseif law === :yoshioka2019_morb
        co_wtp = 10.0^(-5.20 + 0.80 * log10(p_co_bar))
        return co_wtp * 1.0e4 * (28.0101 / 12.011)
    else
        throw(ArgumentError("Unknown CO solubility law: $law"))
    end
end

"""
Compute equilibrium dissolved methane (CH4) in silicate melt.

$(SIGNATURES)

Ardia et al. (2013) haplobasalt fit under strongly reducing conditions:
    X_CH4 [ppmw] = p_CH4 [GPa] * exp(4.93 - 1.93 * p_tot [GPa])

# Arguments
- `p_CH4_Pa`: Partial pressure of CH4 [Pa]
- `p_total_Pa`: Total pressure [Pa]

# Keyword Arguments
- `law`: Formulation (`:ardia2013`)

# Returns
- `ppmw`: Dissolved CH4 concentration in melt [ppmw]
"""
function compute_ch4_solubility_melt(
    p_CH4_Pa::Real, p_total_Pa::Real; law::Symbol=:ardia2013
)::Float64
    p_ch4 = Float64(p_CH4_Pa)
    if !isfinite(p_ch4)
        throw(DomainError(p_ch4, "p_CH4 must be finite"))
    end
    p_tot = Float64(p_total_Pa)
    if !isfinite(p_tot)
        throw(DomainError(p_tot, "p_total must be finite"))
    end
    if p_ch4 <= 0.0
        return 0.0
    end
    p_ch4_gpa = p_ch4 * 1.0e-9
    p_tot_gpa = max(p_tot, 0.0) * 1.0e-9

    if law === :ardia2013
        return p_ch4_gpa * exp(4.93 - 1.93 * p_tot_gpa)
    else
        throw(ArgumentError("Unknown CH4 solubility law: $law"))
    end
end

"""
Compute equilibrium dissolved carbon dioxide (CO2) in silicate melt as carbonate.

$(SIGNATURES)

Dixon et al. (1995) MORB basalt fit:
    x = 3.8e-7 * p_CO2 [bar] * exp(-23.0 * (p_CO2 [bar] - 1.0) / (83.15 * T [K]))
    X_CO2 [ppmw] = 1e4 * (4400.0 * x) / (36.6 - 44.0 * x)

# Arguments
- `p_CO2_Pa`: Partial pressure of CO2 [Pa]
- `T_K`: Melt temperature [K]

# Keyword Arguments
- `law`: Formulation (`:dixon1995`)

# Returns
- `ppmw`: Dissolved CO2 concentration in melt [ppmw]
"""
function compute_co2_solubility_melt(
    p_CO2_Pa::Real, T_K::Real; law::Symbol=:dixon1995
)::Float64
    p_co2 = Float64(p_CO2_Pa)
    if !isfinite(p_co2)
        throw(DomainError(p_co2, "p_CO2 must be finite"))
    end
    T = Float64(T_K)
    if T <= 0.0 || !isfinite(T)
        throw(DomainError(T, "Temperature must be > 0 and finite"))
    end
    if p_co2 <= 0.0
        return 0.0
    end
    p_co2_bar = p_co2 * 1.0e-5

    if law === :dixon1995
        x = 3.8e-7 * p_co2_bar * exp(-23.0 * (p_co2_bar - 1.0) / (83.15 * T))
        denom = 36.6 - 44.0 * x
        if denom <= 0.0
            throw(
                DomainError(
                    denom, "CO2 mole fraction exceeds Dixon (1995) denominator pole"
                ),
            )
        end
        return 1.0e4 * (4400.0 * x) / denom
    else
        throw(ArgumentError("Unknown CO2 solubility law: $law"))
    end
end

"""
Compute equilibrium dissolved carbon in silicate melt under reducing conditions.

$(SIGNATURES)

Partitions dissolved carbon into CO (Armstrong et al. 2015), CH4 (Ardia et al. 2013),
and CO2 (Dixon et al. 1995). If `graphite_saturation` is true, caps CO and CO2 partial
pressures at graphite saturation fugacities calculated via `compute_graphite_saturation_fugacity`.

# Arguments
- `p_CO_Pa`: Partial pressure of CO [Pa]
- `p_CH4_Pa`: Partial pressure of CH4 [Pa]
- `p_CO2_Pa`: Partial pressure of CO2 [Pa]
- `p_total_Pa`: Total pressure [Pa]
- `T_K`: Melt temperature [K]

# Keyword Arguments
- `co_law`: Law for CO (default: `:armstrong2015`)
- `ch4_law`: Law for CH4 (default: `:ardia2013`)
- `co2_law`: Law for CO2 (default: `:dixon1995`)
- `graphite_saturation`: Whether to enforce graphite saturation ceiling (default: false)
- `delta_IW`: Oxygen fugacity offset relative to IW buffer [log10 units] (default: 0.0)

# Returns
- `NamedTuple`: `(; total_ppm, co_ppm, ch4_ppm, co2_ppm)`
"""
function compute_carbon_solubility_melt(
    p_CO_Pa::Real,
    p_CH4_Pa::Real,
    p_CO2_Pa::Real,
    p_total_Pa::Real,
    T_K::Real;
    co_law::Symbol=:armstrong2015,
    ch4_law::Symbol=:ardia2013,
    co2_law::Symbol=:dixon1995,
    graphite_saturation::Bool=false,
    delta_IW::Real=0.0,
)::@NamedTuple{total_ppm::Float64, co_ppm::Float64, ch4_ppm::Float64, co2_ppm::Float64}
    p_co = Float64(p_CO_Pa)
    p_co2 = Float64(p_CO2_Pa)
    if graphite_saturation
        log10_fO2 = compute_iron_wustite_fO2(T_K; delta_IW=delta_IW)
        gr = compute_graphite_saturation_fugacity(T_K, log10_fO2)
        p_co = min(p_co, gr.f_CO_max_bar * 1.0e5)
        p_co2 = min(p_co2, gr.f_CO2_max_bar * 1.0e5)
    end
    co = compute_co_solubility_melt(p_co, p_total_Pa; law=co_law)
    ch4 = compute_ch4_solubility_melt(p_CH4_Pa, p_total_Pa; law=ch4_law)
    co2 = compute_co2_solubility_melt(p_co2, T_K; law=co2_law)
    return (total_ppm=co + ch4 + co2, co_ppm=co, ch4_ppm=ch4, co2_ppm=co2)
end

"""
Compute composite dissolved carbon concentration in silicate melt from total pressure and speciation.

$(SIGNATURES)

# Arguments
- `p_total_Pa`: Total pressure [Pa]
- `T_K`: Temperature [K]
- `delta_IW`: Oxygen fugacity offset relative to IW buffer [log10 units] (default: 0.0)

# Keyword Arguments
- `co_law`: Law for CO solubility (default: `:armstrong2015`)
- `ch4_law`: Law for CH4 solubility (default: `:ardia2013`)
- `co2_law`: Law for CO2 solubility (default: `:dixon1995`)
- `graphite_saturation`: Whether to enforce graphite saturation ceiling (default: false)

# Returns
- `NamedTuple`: `(; total_ppm, co_ppm, ch4_ppm, co2_ppm)`
"""
function compute_carbon_solubility_melt(
    p_total_Pa::Real,
    T_K::Real,
    delta_IW::Real=0.0;
    co_law::Symbol=:armstrong2015,
    ch4_law::Symbol=:ardia2013,
    co2_law::Symbol=:dixon1995,
    graphite_saturation::Bool=false,
)::@NamedTuple{total_ppm::Float64, co_ppm::Float64, ch4_ppm::Float64, co2_ppm::Float64}
    p_tot = Float64(p_total_Pa)
    if p_tot <= 0.0
        return (total_ppm=0.0, co_ppm=0.0, ch4_ppm=0.0, co2_ppm=0.0)
    end
    spec = solve_chnos_speciation(
        p_tot, T_K, delta_IW; graphite_saturation=graphite_saturation
    )
    return compute_carbon_solubility_melt(
        spec.p_CO_Pa,
        spec.p_CH4_Pa,
        spec.p_CO2_Pa,
        p_tot,
        T_K;
        co_law=co_law,
        ch4_law=ch4_law,
        co2_law=co2_law,
        graphite_saturation=graphite_saturation,
        delta_IW=delta_IW,
    )
end

"""
Compute maximum carbon fugacities at graphite saturation (a_C = 1).

$(SIGNATURES)

French (1966) and Holloway et al. (1992) graphite-gas buffer equilibria:
    C(gr) + 1/2 O2 <=> CO  => log10(f_CO_max)  = 5785.0 / T + 4.545 + 0.5 * log10_fO2
    C(gr) + O2     <=> CO2 => log10(f_CO2_max) = 20590.0 / T - 0.043 + log10_fO2

# Arguments
- `T_K`: Melt temperature [K]
- `log10_fO2`: log10 of oxygen fugacity [bar]

# Returns
- `NamedTuple`: `(; f_CO_max_bar, f_CO2_max_bar)`
"""
function compute_graphite_saturation_fugacity(
    T_K::Real, log10_fO2::Real
)::@NamedTuple{f_CO_max_bar::Float64, f_CO2_max_bar::Float64}
    T = Float64(T_K)
    if T <= 0.0 || !isfinite(T)
        throw(DomainError(T, "Temperature must be > 0 and finite"))
    end
    lfO2 = Float64(log10_fO2)
    if !isfinite(lfO2)
        throw(DomainError(lfO2, "log10_fO2 must be finite"))
    end

    log_co = 5785.0 / T + 4.545 + 0.5 * lfO2
    log_co2 = 20590.0 / T - 0.043 + lfO2

    return (f_CO_max_bar=10.0^log_co, f_CO2_max_bar=10.0^log_co2)
end

"""
Compute equilibrium dissolved sulfur in silicate melt.

$(SIGNATURES)

Calculates dissolved sulfur concentration under reducing-to-oxidizing conditions:
- `:boulliung2023`: Boulliung & Wood (2022, 2023) sulfide capacity (and optional sulfate capacity):
    log10(C_S2-) = 0.225 - slope / T
    S_sulfide [wt%] = C_S2- * sqrt(p_S2 [bar] / f_O2 [bar])
- `:gaillard2022`: Gaillard et al. (2022) basaltic melt sulfide capacity:
    ln(S [ppmw]) = 13.8426 - 26476.0 / T + 0.124 * x_FeO + 0.5 * ln(p_S2 [bar] / f_O2 [bar])

# Arguments
- `p_S2_Pa`: Partial pressure of S2 [Pa]
- `T_K`: Melt temperature [K]
- `delta_IW`: Oxygen fugacity offset relative to IW buffer [log10 units]

# Keyword Arguments
- `law`: Formulation (`:boulliung2023` or `:gaillard2022`)
- `sulfide_melt`: Melt composition for Boulliung (`:basalt`, `:andesite`, `:trachybasalt`)
- `include_sulfate`: Whether to add sulfate capacity (relevant above IW+2)
- `x_FeO`: Melt FeO content [wt%] (default: 10.0)
- `scss_active`: Whether to enforce SCSS saturation limit (default: false)
- `p_total_Pa`: Total pressure for SCSS [Pa] (default: 0.0)
- `scss_law`: SCSS formulation (`:smythe2017` or `:oneill2002`)

# Returns
- `S_ppm`: Dissolved sulfur concentration in melt [ppmw]
"""
function compute_sulfur_solubility_melt(
    p_S2_Pa::Real,
    T_K::Real,
    delta_IW::Real;
    law::Symbol=:boulliung2023,
    sulfide_melt::Symbol=:basalt,
    include_sulfate::Bool=false,
    x_FeO::Real=10.0,
    scss_active::Bool=false,
    p_total_Pa::Real=0.0,
    scss_law::Symbol=:smythe2017,
)::Float64
    T = Float64(T_K)
    if T <= 0.0 || !isfinite(T)
        throw(DomainError(T, "Temperature must be > 0 and finite"))
    end
    d_IW = Float64(delta_IW)
    if !isfinite(d_IW) || abs(d_IW) > 50.0
        throw(DomainError(d_IW, "delta_IW must be finite and within [-50, 50]"))
    end
    p_s2 = Float64(p_S2_Pa)
    if !isfinite(p_s2)
        throw(DomainError(p_s2, "p_S2 must be finite"))
    end
    p_s2_bar = max(p_s2, 0.0) * 1.0e-5
    if p_s2_bar < 1.0e-20
        return 0.0
    end
    x_fe = Float64(x_FeO)
    if !isfinite(x_fe)
        throw(DomainError(x_fe, "x_FeO must be finite"))
    end

    log10_fO2 = compute_iron_wustite_fO2(T; delta_IW=d_IW)
    fO2_bar = 10.0^log10_fO2

    s_ppm = if law === :boulliung2023
        slope_s2 = if sulfide_melt === :basalt
            8045.7465
        elseif sulfide_melt === :andesite
            8921.0927
        elseif sulfide_melt === :trachybasalt
            7842.5
        else
            throw(ArgumentError("Unknown sulfide_melt: $sulfide_melt"))
        end
        logC_s2 = 0.225 - slope_s2 / T
        s_wtp = 10.0^(logC_s2 - 0.5 * (log10_fO2 - log10(p_s2_bar)))
        s_base = s_wtp * 1.0e4

        if include_sulfate
            slope_s6 = if sulfide_melt === :basalt
                32333.5635
            elseif sulfide_melt === :andesite
                31586.2393
            elseif sulfide_melt === :trachybasalt
                32446.366
            end
            logC_s6 = -12.948 + slope_s6 / T
            so4_wtp = 10.0^(logC_s6 + 0.5 * log10(p_s2_bar) + 1.5 * log10_fO2)
            s_base += (so4_wtp * (32.065 / 96.06)) * 1.0e4
        end
        s_base
    elseif law === :gaillard2022
        ln_s = 13.8426 - 26476.0 / T + 0.124 * x_fe + 0.5 * log(p_s2_bar / fO2_bar)
        exp(ln_s)
    else
        throw(ArgumentError("Unknown sulfur solubility law: $law"))
    end

    if scss_active
        cap = compute_scss(T, p_total_Pa; x_FeO=x_fe, law=scss_law)
        return min(s_ppm, cap)
    end
    return s_ppm
end

"""
Compute Sulfur Content at Sulfide Saturation (SCSS) in silicate melt.

$(SIGNATURES)

Calculates the maximum dissolved sulfur content before an immiscible Fe-S sulfide liquid
exsolves (O'Neill & Mavrogenes 2002; Fortin et al. 2015; Smythe et al. 2017):
- `:smythe2017`: Smythe et al. (2017) pressure-dependent formulation:
    ln(SCSS [ppmw]) = 7.50 - 4500.0 / T + 0.90 * ln(max(0.1, x_FeO)) - 2.5e-4 * (P_tot [bar] / T)
- `:oneill2002`: O'Neill & Mavrogenes (2002) 1-bar baseline:
    ln(SCSS [ppmw]) = 7.50 - 4500.0 / T + 0.90 * ln(max(0.1, x_FeO))

# Arguments
- `T_K`: Melt temperature [K]
- `p_total_Pa`: Total pressure [Pa]

# Keyword Arguments
- `x_FeO`: Silicate melt FeO content [wt%] (default: 10.0)
- `law`: SCSS formulation (`:smythe2017` or `:oneill2002`)

# Returns
- `scss_ppm`: Maximum dissolved sulfur in melt [ppmw]
"""
function compute_scss(
    T_K::Real, p_total_Pa::Real; x_FeO::Real=10.0, law::Symbol=:smythe2017
)::Float64
    T = Float64(T_K)
    if T <= 0.0 || !isfinite(T)
        throw(DomainError(T, "Temperature must be > 0 and finite"))
    end
    p_tot = Float64(p_total_Pa)
    if !isfinite(p_tot)
        throw(DomainError(p_tot, "p_total must be finite"))
    end
    x_fe = Float64(x_FeO)
    if x_fe < 0.0 || !isfinite(x_fe)
        throw(DomainError(x_fe, "x_FeO must be >= 0 and finite"))
    end
    p_bar = max(p_tot, 0.0) * 1.0e-5
    fe_term = log(max(x_fe, 0.1))

    if law === :smythe2017
        ln_scss = 7.50 - 4500.0 / T + 0.90 * fe_term - 2.5e-4 * (p_bar / T)
        return exp(ln_scss)
    elseif law === :oneill2002
        ln_scss = 7.50 - 4500.0 / T + 0.90 * fe_term
        return exp(ln_scss)
    else
        throw(ArgumentError("Unknown SCSS law: $law"))
    end
end

"""
Compute the thermodynamic volatile retention floor [ppmw] in nominally anhydrous minerals (NAMs)
and refractory solid phases for species `species` (`:H2O`, `:C`, `:N`, `:S`) at temperature `T_val` [K],
melt fraction `F_melt` [-], and pressure `P_val` [Pa].

$(SIGNATURES)

# Details
- For water (`:H2O`): models hydroxyl defect retention in nominally anhydrous minerals (olivine,
  pyroxene) following Hirschmann et al. (2006) and Peslier et al. (2017).
- For carbon (`:C`): models refractory graphite and interstitial carbon retention in the solid
  silicate lattice following Shcheka et al. (2006) and Hirschmann (2018).
- For nitrogen (`:N`): models lattice-bound nitrogen and refractory nitride retention (Li et al. 2013).
- For sulfur (`:S`): models monosulfide solid solution (MSS) and refractory sulfide retention.

# Retention Laws (`cfg.retention_law`):
- `:constant_floor`: returns constant retention floor `cfg.<species>_retention_ppm`.
- `:nams_exponential`: near/below `T_solidus_ref`, returns baseline floor; for `T > T_solidus_ref`,
  decays exponentially as `C_ret0 * exp(-(T - T_solidus_ref) / dT_retention)`.
- `:linear_melt_blend`: scales as `C_ret0 * max(0.0, 1.0 - F_melt)`.

# Returns
- `C_ret_ppm::Float64`: Retained volatile concentration in solid matrix [ppmw]
"""
function compute_volatile_retention_floor(
    T_val::Real, species::Symbol, cfg::RetentionConfig; F_melt::Real=0.0, P_val::Real=0.0
)::Float64
    if !cfg.active
        return 0.0
    end
    if !isfinite(Float64(T_val)) || Float64(T_val) < 0.0
        throw(
            DomainError(
                T_val, "Temperature T_val must be non-negative and finite, got $T_val K"
            ),
        )
    end
    if !isfinite(Float64(F_melt)) || Float64(F_melt) < 0.0
        throw(
            DomainError(
                F_melt, "Melt fraction F_melt must be non-negative and finite, got $F_melt"
            ),
        )
    end
    if !isfinite(Float64(P_val)) || Float64(P_val) < 0.0
        throw(
            DomainError(
                P_val, "Pressure P_val must be non-negative and finite, got $P_val Pa"
            ),
        )
    end
    T_k = Float64(T_val)
    F_m = min(1.0, Float64(F_melt))

    C_base = if species === :H2O
        cfg.h2o_retention_ppm
    elseif species === :C
        cfg.carbon_retention_ppm
    elseif species === :N
        cfg.nitrogen_retention_ppm
    elseif species === :S
        cfg.sulfur_retention_ppm
    else
        throw(
            ArgumentError(
                "Unknown volatile species for retention floor: $species. Expected :H2O, :C, :N, or :S.",
            ),
        )
    end

    if C_base <= 0.0
        return 0.0
    end

    if cfg.retention_law === :constant_floor
        return C_base
    elseif cfg.retention_law === :linear_melt_blend
        return C_base * max(0.0, 1.0 - F_m)
    elseif cfg.retention_law === :nams_exponential
        T_sol = cfg.T_solidus_ref
        dT = cfg.dT_retention
        if T_k <= T_sol
            return C_base
        else
            arg = -(T_k - T_sol) / dT
            return C_base * exp(clamp(arg, -40.0, 0.0))
        end
    else
        throw(ArgumentError("Unknown retention_law: $(cfg.retention_law)"))
    end
end

"""
Compute water retention floor in nominally anhydrous minerals [ppmw].
"""
compute_h2o_retention_floor(
    T_val::Real, cfg::RetentionConfig; F_melt::Real=0.0, P_val::Real=0.0
)::Float64 = compute_volatile_retention_floor(T_val, :H2O, cfg; F_melt=F_melt, P_val=P_val)

"""
Compute carbon retention floor in refractory solid phases [ppmw].
"""
compute_carbon_retention_floor(
    T_val::Real, cfg::RetentionConfig; F_melt::Real=0.0, P_val::Real=0.0
)::Float64 = compute_volatile_retention_floor(T_val, :C, cfg; F_melt=F_melt, P_val=P_val)

"""
Compute nitrogen retention floor in mineral lattice and nitrides [ppmw].
"""
compute_nitrogen_retention_floor(
    T_val::Real, cfg::RetentionConfig; F_melt::Real=0.0, P_val::Real=0.0
)::Float64 = compute_volatile_retention_floor(T_val, :N, cfg; F_melt=F_melt, P_val=P_val)

"""
Compute sulfur retention floor in solid sulfides and MSS [ppmw].
"""
compute_sulfur_retention_floor(
    T_val::Real, cfg::RetentionConfig; F_melt::Real=0.0, P_val::Real=0.0
)::Float64 = compute_volatile_retention_floor(T_val, :S, cfg; F_melt=F_melt, P_val=P_val)

"""
Compute equilibrium volatile exsolution from silicate melt for H-C-N-S volatile species.

$(SIGNATURES)

When silicate melting occurs (`F_melt > 0`), dissolved volatiles partition into the melt phase.
If the volatile concentration in the melt exceeds the saturation solubility at local pore pressure `P_Pa`
and temperature `T_K`, the excess volatile mass exsolves into the pore fluid phase.
When retention floor modeling is active, exsolution is bounded by the mobile excess above the solid
retention floor, preventing unphysical total dehydration or decarbonation.

# Arguments
- `F_melt`: Silicate melt volume/mass fraction in [0, 1]
- `P_Pa`: Local pore fluid / ambient pressure [Pa]
- `T_K`: Local temperature [K]
- `w_H2O_bulk`: Bulk rock water mass fraction [-] (e.g. 0.01 for 1 wt%)
- `C_C_bulk_ppm`: Bulk rock carbon concentration [ppm]
- `C_N_bulk_ppm`: Bulk rock nitrogen concentration [ppm]
- `C_S_bulk_ppm`: Bulk rock sulfur concentration [ppm]
- `delta_IW`: Oxygen fugacity offset relative to IW buffer

# Keyword Arguments
- `water_law`: Water solubility law (default: `:burnham_dixon`)
- `water_As`: Water solubility coefficient (default: 0.40)
- `carbon_active`: Whether carbon solubility is modeled (default: false)
- `co_law`: Carbon monoxide solubility law (default: `:armstrong2015`)
- `ch4_law`: Methane solubility law (default: `:ardia2013`)
- `co2_law`: Carbon dioxide solubility law (default: `:dixon1995`)
- `nitrogen_law`: Nitrogen solubility law (default: `:dasgupta2022`)
- `nitrogen_henry`: Nitrogen Henry coefficient (default: 0.40)
- `nitrogen_nitride`: Nitrogen nitride capacity (default: 1.0e-3)
- `sulfur_active`: Whether sulfur solubility is modeled (default: false)
- `sulfide_law`: Sulfide solubility law (default: `:boulliung2023`)
- `graphite_saturation`: Whether carbon is capped at graphite saturation (default: true)
- `retention_active`: Whether thermodynamic volatile retention floors are active (default: false)
- `retention_cfg`: RetentionConfig struct (default: nothing)

# Returns
- `NamedTuple`:
  - `w_H2O_ex`: Exsolved water mass fraction [-]
  - `w_C_ex`: Exsolved carbon mass fraction [-]
  - `w_N_ex`: Exsolved nitrogen mass fraction [-]
  - `w_S_ex`: Exsolved sulfur mass fraction [-]
  - `w_total_ex`: Total exsolved volatile mass fraction [-]
  - `w_H2O_diss`: Retained dissolved water mass fraction [-]
  - `C_C_diss_ppm`: Retained dissolved carbon concentration [ppm]
  - `C_N_diss_ppm`: Retained dissolved nitrogen concentration [ppm]
  - `C_S_diss_ppm`: Retained dissolved sulfur concentration [ppm]
  - `z_H`, `z_C`, `z_N`, `z_S`: Elemental atom fractions of the exsolved gas
"""
function compute_volatile_exsolution(
    F_melt::Real,
    P_Pa::Real,
    T_K::Real,
    w_H2O_bulk::Real,
    C_C_bulk_ppm::Real,
    C_N_bulk_ppm::Real,
    C_S_bulk_ppm::Real,
    delta_IW::Real;
    water_law::Symbol=:burnham_dixon,
    water_As::Real=0.40,
    carbon_active::Bool=false,
    co_law::Symbol=:armstrong2015,
    ch4_law::Symbol=:ardia2013,
    co2_law::Symbol=:dixon1995,
    nitrogen_law::Symbol=:dasgupta2022,
    nitrogen_henry::Real=0.40,
    nitrogen_nitride::Real=1.0e-3,
    sulfur_active::Bool=false,
    sulfide_law::Symbol=:boulliung2023,
    graphite_saturation::Bool=true,
    retention_active::Bool=false,
    retention_cfg::Union{Nothing,RetentionConfig}=nothing,
)::@NamedTuple{
    w_H2O_ex::Float64,
    w_C_ex::Float64,
    w_N_ex::Float64,
    w_S_ex::Float64,
    w_total_ex::Float64,
    w_H2O_diss::Float64,
    C_C_diss_ppm::Float64,
    C_N_diss_ppm::Float64,
    C_S_diss_ppm::Float64,
    z_H::Float64,
    z_C::Float64,
    z_N::Float64,
    z_S::Float64,
}
    F_m = clamp(Float64(F_melt), 0.0, 1.0)
    P_val = max(Float64(P_Pa), 0.0)
    T_val = Float64(T_K)
    d_IW = Float64(delta_IW)
    w_H2O = max(Float64(w_H2O_bulk), 0.0)
    C_C = max(Float64(C_C_bulk_ppm), 0.0)
    C_N = max(Float64(C_N_bulk_ppm), 0.0)
    C_S = max(Float64(C_S_bulk_ppm), 0.0)

    if F_m <= 0.0 || T_val <= 0.0
        return (
            w_H2O_ex=0.0,
            w_C_ex=0.0,
            w_N_ex=0.0,
            w_S_ex=0.0,
            w_total_ex=0.0,
            w_H2O_diss=w_H2O,
            C_C_diss_ppm=C_C,
            C_N_diss_ppm=C_N,
            C_S_diss_ppm=C_S,
            z_H=0.80,
            z_C=0.15,
            z_N=0.03,
            z_S=0.02,
        )
    end

    # Evaluate retention floors if retention active
    w_ret_act_H2O = 0.0
    w_H2O_mob = w_H2O
    C_ret_act_N = 0.0
    C_N_mob = C_N
    C_ret_act_C = 0.0
    C_C_mob = C_C
    C_ret_act_S = 0.0
    C_S_mob = C_S

    if retention_active && retention_cfg !== nothing && retention_cfg.active
        C_ret_H2O_ppm = compute_h2o_retention_floor(
            T_val, retention_cfg; F_melt=F_m, P_val=P_val
        )
        w_ret_H2O = C_ret_H2O_ppm * 1.0e-6
        w_ret_act_H2O = min(w_H2O, w_ret_H2O)
        w_H2O_mob = max(0.0, w_H2O - w_ret_act_H2O)

        C_ret_N = compute_nitrogen_retention_floor(
            T_val, retention_cfg; F_melt=F_m, P_val=P_val
        )
        C_ret_act_N = min(C_N, C_ret_N)
        C_N_mob = max(0.0, C_N - C_ret_act_N)

        C_ret_C = compute_carbon_retention_floor(
            T_val, retention_cfg; F_melt=F_m, P_val=P_val
        )
        C_ret_act_C = min(C_C, C_ret_C)
        C_C_mob = max(0.0, C_C - C_ret_act_C)

        C_ret_S = compute_sulfur_retention_floor(
            T_val, retention_cfg; F_melt=F_m, P_val=P_val
        )
        C_ret_act_S = min(C_S, C_ret_S)
        C_S_mob = max(0.0, C_S - C_ret_act_S)
    end

    # 1. Water solubility
    S_H2O_wtpct = compute_water_solubility_melt(P_val; As=water_As, law=water_law)
    S_H2O_frac = S_H2O_wtpct * 0.01
    cap_H2O = F_m * S_H2O_frac
    w_H2O_ex = max(0.0, w_H2O_mob - cap_H2O)
    w_H2O_diss = w_ret_act_H2O + min(w_H2O_mob, cap_H2O)

    # 2. Nitrogen solubility
    S_N_res = compute_nitrogen_solubility_melt(
        P_val, d_IW; Kh=nitrogen_henry, C_nitride=nitrogen_nitride
    )
    cap_N = F_m * S_N_res.total_ppm
    C_N_ex = max(0.0, C_N_mob - cap_N)
    C_N_diss = C_ret_act_N + min(C_N_mob, cap_N)
    w_N_ex = C_N_ex * 1.0e-6

    # 3. Carbon solubility
    w_C_ex = 0.0
    C_C_diss = C_C
    if carbon_active
        S_C_res = compute_carbon_solubility_melt(
            P_val,
            T_val,
            d_IW;
            co_law=co_law,
            ch4_law=ch4_law,
            co2_law=co2_law,
            graphite_saturation=graphite_saturation,
        )
        cap_C = F_m * S_C_res.total_ppm
        C_C_ex = max(0.0, C_C_mob - cap_C)
        C_C_diss = C_ret_act_C + min(C_C_mob, cap_C)
        w_C_ex = C_C_ex * 1.0e-6
    end

    # 4. Sulfur solubility
    w_S_ex = 0.0
    C_S_diss = C_S
    if sulfur_active
        S_S_ppm = compute_sulfur_solubility_melt(P_val, T_val, d_IW; law=sulfide_law)
        cap_S = F_m * S_S_ppm
        C_S_ex = max(0.0, C_S_mob - cap_S)
        C_S_diss = C_ret_act_S + min(C_S_mob, cap_S)
        w_S_ex = C_S_ex * 1.0e-6
    end

    w_total_ex = w_H2O_ex + w_C_ex + w_N_ex + w_S_ex

    # Elemental atom counts of exsolved gas
    mol_H = 2.0 * (w_H2O_ex / 0.01801528)
    mol_C = w_C_ex / 0.012011
    mol_N = w_N_ex / 0.014007
    mol_S = w_S_ex / 0.032065
    mol_tot = mol_H + mol_C + mol_N + mol_S

    z_H, z_C, z_N, z_S = if mol_tot > 0.0
        (mol_H / mol_tot, mol_C / mol_tot, mol_N / mol_tot, mol_S / mol_tot)
    else
        (0.80, 0.15, 0.03, 0.02)
    end

    return (
        w_H2O_ex=w_H2O_ex,
        w_C_ex=w_C_ex,
        w_N_ex=w_N_ex,
        w_S_ex=w_S_ex,
        w_total_ex=w_total_ex,
        w_H2O_diss=w_H2O_diss,
        C_C_diss_ppm=C_C_diss,
        C_N_diss_ppm=C_N_diss,
        C_S_diss_ppm=C_S_diss,
        z_H=z_H,
        z_C=z_C,
        z_N=z_N,
        z_S=z_S,
    )
end

"""
Solve homogeneous gas-phase chemical equilibrium for the C-H-O-N-S volatile system.

$(SIGNATURES)

Given elemental gas fractions `z_H, z_C, z_N, z_S`, total pressure `p_total_Pa`, melt temperature `T_K`,
and oxygen fugacity offset `delta_IW`, solves for partial pressures of major outgassed species:
`H2, H2O, CO, CO2, CH4, N2, NH3, H2S, S2, SO2` while enforcing Dalton's law of partial
pressures (∑ p_i = p_total) and simultaneous atomic mass conservation for H, C, N, and S.

# Arguments
- `p_total_Pa`: Total gas pressure [Pa]
- `T_K`: Gas temperature [K]
- `delta_IW`: Oxygen fugacity offset relative to IW buffer [log10 units]

# Keyword Arguments
- `z_H`: Elemental hydrogen fraction (default: 0.80)
- `z_C`: Elemental carbon fraction (default: 0.15)
- `z_N`: Elemental nitrogen fraction (default: 0.03)
- `z_S`: Elemental sulfur fraction (default: 0.02)
- `graphite_saturation`: Whether to cap C fugacities at graphite saturation (default: true)

# Returns
- `NamedTuple`: `(; p_H2_Pa, p_H2O_Pa, p_CO_Pa, p_CO2_Pa, p_CH4_Pa, p_N2_Pa, p_NH3_Pa, p_H2S_Pa, p_S2_Pa, p_SO2_Pa)`
"""
function solve_chnos_speciation(
    p_total_Pa::Real,
    T_K::Real,
    delta_IW::Real;
    z_H::Real=0.80,
    z_C::Real=0.15,
    z_N::Real=0.03,
    z_S::Real=0.02,
    graphite_saturation::Bool=true,
)::@NamedTuple{
    p_H2_Pa::Float64,
    p_H2O_Pa::Float64,
    p_CO_Pa::Float64,
    p_CO2_Pa::Float64,
    p_CH4_Pa::Float64,
    p_N2_Pa::Float64,
    p_NH3_Pa::Float64,
    p_H2S_Pa::Float64,
    p_S2_Pa::Float64,
    p_SO2_Pa::Float64,
}
    p_tot = Float64(p_total_Pa)
    if !isfinite(p_tot)
        throw(DomainError(p_tot, "p_total must be finite"))
    end
    T = Float64(T_K)
    if T <= 0.0 || !isfinite(T)
        throw(DomainError(T, "Temperature must be > 0 and finite"))
    end
    d_IW = Float64(delta_IW)
    if !isfinite(d_IW) || abs(d_IW) > 50.0
        throw(DomainError(d_IW, "delta_IW must be finite and within [-50, 50]"))
    end
    if p_tot <= 0.0
        return (
            p_H2_Pa=0.0,
            p_H2O_Pa=0.0,
            p_CO_Pa=0.0,
            p_CO2_Pa=0.0,
            p_CH4_Pa=0.0,
            p_N2_Pa=0.0,
            p_NH3_Pa=0.0,
            p_H2S_Pa=0.0,
            p_S2_Pa=0.0,
            p_SO2_Pa=0.0,
        )
    end

    sum_z = Float64(z_H) + Float64(z_C) + Float64(z_N) + Float64(z_S)
    if sum_z <= 0.0 || !isfinite(sum_z)
        throw(
            DomainError(
                sum_z, "Sum of volatile elemental abundances must be > 0 and finite"
            ),
        )
    end
    nH = Float64(z_H) / sum_z
    nC = Float64(z_C) / sum_z
    nN = Float64(z_N) / sum_z
    nS = Float64(z_S) / sum_z

    log10_fO2 = compute_iron_wustite_fO2(T; delta_IW=d_IW)
    p_tot_bar = p_tot * 1.0e-5

    logK_H2O = 12700.0 / T - 2.80
    r_H = 10.0^clamp(logK_H2O + 0.5 * log10_fO2, -100.0, 100.0)

    logK_CO2 = 14800.0 / T - 4.58
    r_CO2 = 10.0^clamp(logK_CO2 + 0.5 * log10_fO2, -100.0, 100.0)

    logK_SO2 = 18800.0 / T - 3.80
    r_SO2 = 10.0^clamp(logK_SO2 + log10_fO2, -100.0, 100.0)

    # Initial guess for total atomic pressure
    A_atoms = 2.0 * p_tot_bar
    pH2 = nH > 0.0 ? (nH * A_atoms) / (2.0 * (1.0 + r_H)) : 0.0

    p_CO = 0.0
    p_CO2 = 0.0
    p_CH4 = 0.0
    p_N2 = 0.0
    p_NH3 = 0.0
    p_S2 = 0.0
    p_H2S = 0.0
    p_SO2 = 0.0
    p_H2O = 0.0

    # Simultaneous element conservation and Dalton law iteration
    for outer_iter in 1:100
        for inner_iter in 1:40
            log_pH2 = pH2 > 0.0 ? log10(max(pH2, 1.0e-30)) : -100.0
            r_CH4 = if pH2 > 0.0
                10.0^clamp(
                    11500.0 / T - 12.0 + 2.0 * log_pH2 - log10(max(r_H, 1.0e-30)),
                    -100.0,
                    100.0,
                )
            else
                0.0
            end
            r_NH3 = if pH2 > 0.0
                10.0^clamp(2800.0 / T - 5.80 + 1.5 * log_pH2, -100.0, 100.0)
            else
                0.0
            end
            r_H2S = pH2 > 0.0 ? 10.0^clamp(4800.0 / T - 2.50 + log_pH2, -100.0, 100.0) : 0.0

            if nC > 0.0
                p_CO = (nC * A_atoms) / (1.0 + r_CO2 + r_CH4)
                p_CO2 = r_CO2 * p_CO
                p_CH4 = r_CH4 * p_CO
            else
                p_CO = 0.0
                p_CO2 = 0.0
                p_CH4 = 0.0
            end

            if nN > 0.0
                A_N = nN * A_atoms
                denom_N = r_NH3 + sqrt(r_NH3^2 + 8.0 * A_N)
                u_N = (2.0 * A_N) / max(denom_N, 1.0e-30)
                p_N2 = u_N^2
                p_NH3 = r_NH3 * u_N
            else
                p_N2 = 0.0
                p_NH3 = 0.0
            end

            if nS > 0.0
                A_S = nS * A_atoms
                B_S = r_H2S + r_SO2
                denom_S = B_S + sqrt(B_S^2 + 8.0 * A_S)
                v_S = (2.0 * A_S) / max(denom_S, 1.0e-30)
                p_S2 = v_S^2
                p_H2S = r_H2S * v_S
                p_SO2 = r_SO2 * v_S
            else
                p_S2 = 0.0
                p_H2S = 0.0
                p_SO2 = 0.0
            end

            if nH > 0.0
                H_sequestered = 4.0 * p_CH4 + 3.0 * p_NH3 + 2.0 * p_H2S
                H_avail = max(0.0, nH * A_atoms - H_sequestered)
                pH2_new = H_avail / (2.0 * (1.0 + r_H))
                diff = abs(pH2_new - pH2)
                pH2 = 0.5 * (pH2 + pH2_new)
                if diff < 1.0e-13 * p_tot_bar
                    break
                end
            else
                pH2 = 0.0
                break
            end
        end

        p_H2O = r_H * pH2
        p_calc = pH2 + p_H2O + p_CO + p_CO2 + p_CH4 + p_N2 + p_NH3 + p_H2S + p_S2 + p_SO2
        err = abs(p_calc - p_tot_bar) / p_tot_bar
        if err < 1.0e-12
            break
        end
        A_atoms *= (p_tot_bar / p_calc)
    end

    # If graphite saturation is enabled, verify carbon activity a_C <= 1.
    # When uncapped p_CO exceeds the graphite saturation ceiling, elemental carbon precipitates
    # as solid graphite. The gas phase carbon partial pressures are fixed by equilibrium with graphite,
    # and the remaining pressure is partitioned among volatile elements (H, N, S) preserving their
    # relative abundances and satisfying Dalton's law exactly.
    if graphite_saturation && nC > 0.0
        gr = compute_graphite_saturation_fugacity(T, log10_fO2)
        f_co_max_bar = gr.f_CO_max_bar
        if p_CO > f_co_max_bar
            p_CO_sat_bar = f_co_max_bar
            p_CO2_sat_bar = r_CO2 * p_CO_sat_bar
            z_HNS = Float64(z_H) + Float64(z_N) + Float64(z_S)
            if z_HNS <= 0.0
                p_C_tot_bar = p_CO_sat_bar + p_CO2_sat_bar
                scale_sat = p_tot_bar / max(p_C_tot_bar, 1.0e-30)
                return (
                    p_H2_Pa=0.0,
                    p_H2O_Pa=0.0,
                    p_CO_Pa=p_CO_sat_bar * scale_sat * 1.0e5,
                    p_CO2_Pa=p_CO2_sat_bar * scale_sat * 1.0e5,
                    p_CH4_Pa=0.0,
                    p_N2_Pa=0.0,
                    p_NH3_Pa=0.0,
                    p_H2S_Pa=0.0,
                    p_S2_Pa=0.0,
                    p_SO2_Pa=0.0,
                )
            end

            # Solve for pH2 with monotonic 1D bisection
            function _sat_residual(test_pH2_bar)
                l_pH2 = test_pH2_bar > 0.0 ? log10(max(test_pH2_bar, 1.0e-30)) : -100.0
                r_ch4_test = if test_pH2_bar > 0.0
                    10.0^clamp(
                        11500.0 / T - 12.0 + 2.0 * l_pH2 - log10(max(r_H, 1.0e-30)),
                        -100.0,
                        100.0,
                    )
                else
                    0.0
                end
                p_C_test = p_CO_sat_bar + p_CO2_sat_bar + r_ch4_test * p_CO_sat_bar
                p_rem_test = p_tot_bar - p_C_test
                if p_rem_test <= 0.0
                    return p_C_test - p_tot_bar
                end
                hns_res = solve_chnos_speciation(
                    p_rem_test * 1.0e5,
                    T,
                    d_IW;
                    z_H=z_H,
                    z_C=0.0,
                    z_N=z_N,
                    z_S=z_S,
                    graphite_saturation=false,
                )
                return test_pH2_bar - (hns_res.p_H2_Pa * 1.0e-5)
            end

            lo_sat = 0.0
            hi_sat = p_tot_bar
            for _ in 1:60
                mid_sat = 0.5 * (lo_sat + hi_sat)
                if _sat_residual(mid_sat) > 0.0
                    hi_sat = mid_sat
                else
                    lo_sat = mid_sat
                end
            end
            best_pH2_bar = 0.5 * (lo_sat + hi_sat)
            l_pH2_final = best_pH2_bar > 0.0 ? log10(max(best_pH2_bar, 1.0e-30)) : -100.0
            r_ch4_final = if best_pH2_bar > 0.0
                10.0^clamp(
                    11500.0 / T - 12.0 + 2.0 * l_pH2_final - log10(max(r_H, 1.0e-30)),
                    -100.0,
                    100.0,
                )
            else
                0.0
            end
            p_CH4_sat_bar = r_ch4_final * p_CO_sat_bar
            p_C_final = p_CO_sat_bar + p_CO2_sat_bar + p_CH4_sat_bar

            if p_C_final >= p_tot_bar
                scale_sat = p_tot_bar / p_C_final
                return (
                    p_H2_Pa=0.0,
                    p_H2O_Pa=0.0,
                    p_CO_Pa=p_CO_sat_bar * scale_sat * 1.0e5,
                    p_CO2_Pa=p_CO2_sat_bar * scale_sat * 1.0e5,
                    p_CH4_Pa=p_CH4_sat_bar * scale_sat * 1.0e5,
                    p_N2_Pa=0.0,
                    p_NH3_Pa=0.0,
                    p_H2S_Pa=0.0,
                    p_S2_Pa=0.0,
                    p_SO2_Pa=0.0,
                )
            end

            p_rem_final_bar = p_tot_bar - p_C_final
            hns_final = solve_chnos_speciation(
                p_rem_final_bar * 1.0e5,
                T,
                d_IW;
                z_H=z_H,
                z_C=0.0,
                z_N=z_N,
                z_S=z_S,
                graphite_saturation=false,
            )
            return (
                p_H2_Pa=hns_final.p_H2_Pa,
                p_H2O_Pa=hns_final.p_H2O_Pa,
                p_CO_Pa=p_CO_sat_bar * 1.0e5,
                p_CO2_Pa=p_CO2_sat_bar * 1.0e5,
                p_CH4_Pa=p_CH4_sat_bar * 1.0e5,
                p_N2_Pa=hns_final.p_N2_Pa,
                p_NH3_Pa=hns_final.p_NH3_Pa,
                p_H2S_Pa=hns_final.p_H2S_Pa,
                p_S2_Pa=hns_final.p_S2_Pa,
                p_SO2_Pa=hns_final.p_SO2_Pa,
            )
        end
    end

    p_calc = pH2 + p_H2O + p_CO + p_CO2 + p_CH4 + p_N2 + p_NH3 + p_H2S + p_S2 + p_SO2
    if p_calc > 0.0
        norm = p_tot_bar / p_calc
        pH2 *= norm
        p_H2O *= norm
        p_CO *= norm
        p_CO2 *= norm
        p_CH4 *= norm
        p_N2 *= norm
        p_NH3 *= norm
        p_H2S *= norm
        p_S2 *= norm
        p_SO2 *= norm
    end

    return (
        p_H2_Pa=pH2 * 1.0e5,
        p_H2O_Pa=p_H2O * 1.0e5,
        p_CO_Pa=p_CO * 1.0e5,
        p_CO2_Pa=p_CO2 * 1.0e5,
        p_CH4_Pa=p_CH4 * 1.0e5,
        p_N2_Pa=p_N2 * 1.0e5,
        p_NH3_Pa=p_NH3 * 1.0e5,
        p_H2S_Pa=p_H2S * 1.0e5,
        p_S2_Pa=p_S2 * 1.0e5,
        p_SO2_Pa=p_SO2 * 1.0e5,
    )
end

"""
    speciate_vented_volatiles(
        m_H2O::Real,
        m_C::Real,
        m_N::Real,
        m_S::Real,
        P_amb_Pa::Real,
        T_surf_K::Real,
        delta_IW::Real;
        graphite_saturation::Bool=true,
    )::Dict{Symbol,Float64}

Calculate equilibrium molecular speciation of vented volatile mass fluxes.
Partitions elemental C, N, S and water mass releases into gaseous species (H2, H2O, CO,
CO2, CH4, N2, NH3, H2S, S2, SO2) at exsolution temperature, pressure, and oxygen fugacity.
When graphite saturation precipitates solid carbon under reducing conditions, gas-phase carbon
is governed by graphite equilibrium, and gas moles are scaled to the non-condensing carrier element.
"""
function speciate_vented_volatiles(
    m_H2O::Real,
    m_C::Real,
    m_N::Real,
    m_S::Real,
    P_amb_Pa::Real,
    T_surf_K::Real,
    delta_IW::Real=0.0;
    graphite_saturation::Bool=true,
)::Dict{Symbol,Float64}
    m_h2o = Float64(m_H2O)
    m_c = Float64(m_C)
    m_n = Float64(m_N)
    m_s = Float64(m_S)
    p_amb = Float64(P_amb_Pa)
    t_surf = Float64(T_surf_K)
    d_iw = Float64(delta_IW)

    if m_h2o < 0.0 || !isfinite(m_h2o)
        throw(DomainError(m_H2O, "Vented H2O mass must be non-negative and finite"))
    end
    if m_c < 0.0 || !isfinite(m_c)
        throw(DomainError(m_C, "Vented C mass must be non-negative and finite"))
    end
    if m_n < 0.0 || !isfinite(m_n)
        throw(DomainError(m_N, "Vented N mass must be non-negative and finite"))
    end
    if m_s < 0.0 || !isfinite(m_s)
        throw(DomainError(m_S, "Vented S mass must be non-negative and finite"))
    end

    species_dict = Dict{Symbol,Float64}(sp => 0.0 for sp in SPECIATION_SPECIES)
    if (m_h2o + m_c + m_n + m_s) <= 0.0
        return species_dict
    end

    if p_amb <= 0.0 || !isfinite(p_amb)
        throw(
            DomainError(P_amb_Pa, "Ambient pressure must be strictly positive and finite")
        )
    end
    if t_surf <= 0.0 || !isfinite(t_surf)
        throw(
            DomainError(
                T_surf_K, "Surface temperature must be strictly positive and finite"
            ),
        )
    end
    if !isfinite(d_iw) || abs(d_iw) > 50.0
        throw(DomainError(delta_IW, "delta_IW must be finite and within [-50, 50]"))
    end

    # Molar elemental amounts released
    nH = 2.0 * m_h2o / 18.01528e-3
    nC = m_c / 12.011e-3
    nN = m_n / 14.007e-3
    nS = m_s / 32.06e-3
    n_tot = nH + nC + nN + nS

    if n_tot <= 0.0
        return species_dict
    end

    z_H = nH / n_tot
    z_C = nC / n_tot
    z_N = nN / n_tot
    z_S = nS / n_tot

    p_amb_eval = max(p_amb, 1.0)
    t_surf_eval = max(t_surf, 273.15)

    spec = solve_chnos_speciation(
        p_amb_eval,
        t_surf_eval,
        d_iw;
        z_H=z_H,
        z_C=z_C,
        z_N=z_N,
        z_S=z_S,
        graphite_saturation=graphite_saturation,
    )

    p_sum = (
        spec.p_H2_Pa +
        spec.p_H2O_Pa +
        spec.p_CO_Pa +
        spec.p_CO2_Pa +
        spec.p_CH4_Pa +
        spec.p_N2_Pa +
        spec.p_NH3_Pa +
        spec.p_H2S_Pa +
        spec.p_S2_Pa +
        spec.p_SO2_Pa
    )

    if p_sum <= 0.0
        # Fallback to simple stoichiometric partition if partial pressures degenerate
        species_dict[:H2O] = m_h2o
        species_dict[:CO2] = m_c * (44.0095 / 12.011)
        species_dict[:N2] = m_n
        species_dict[:H2S] = m_s * (34.08 / 32.06)
        return species_dict
    end

    # Moles of elements per mole of gas
    y_H2 = spec.p_H2_Pa / p_sum
    y_H2O = spec.p_H2O_Pa / p_sum
    y_CO = spec.p_CO_Pa / p_sum
    y_CO2 = spec.p_CO2_Pa / p_sum
    y_CH4 = spec.p_CH4_Pa / p_sum
    y_N2 = spec.p_N2_Pa / p_sum
    y_NH3 = spec.p_NH3_Pa / p_sum
    y_H2S = spec.p_H2S_Pa / p_sum
    y_S2 = spec.p_S2_Pa / p_sum
    y_SO2 = spec.p_SO2_Pa / p_sum

    c_H = 2.0 * y_H2 + 2.0 * y_H2O + 4.0 * y_CH4 + 3.0 * y_NH3 + 2.0 * y_H2S
    c_C = y_CO + y_CO2 + y_CH4
    c_N = 2.0 * y_N2 + y_NH3
    c_S = y_H2S + 2.0 * y_S2 + y_SO2
    c_elem = c_H + c_C + c_N + c_S

    # Scale gas moles from non-condensing carrier elements when graphite saturation occurs.
    N_gas = if graphite_saturation
        if nH > 0.0 && c_H > 0.0
            nH / c_H
        elseif nN > 0.0 && c_N > 0.0
            nN / c_N
        elseif nS > 0.0 && c_S > 0.0
            nS / c_S
        elseif nC > 0.0 && c_C > 0.0
            nC / c_C
        else
            0.0
        end
    else
        c_elem > 0.0 ? n_tot / c_elem : 0.0
    end

    species_dict[:H2] = N_gas * y_H2 * 2.01588e-3
    species_dict[:H2O] = N_gas * y_H2O * 18.01528e-3
    species_dict[:CO] = N_gas * y_CO * 28.0101e-3
    species_dict[:CO2] = N_gas * y_CO2 * 44.0095e-3
    species_dict[:CH4] = N_gas * y_CH4 * 16.0425e-3
    species_dict[:N2] = N_gas * y_N2 * 28.0134e-3
    species_dict[:NH3] = N_gas * y_NH3 * 17.0305e-3
    species_dict[:H2S] = N_gas * y_H2S * 34.0809e-3
    species_dict[:S2] = N_gas * y_S2 * 64.12e-3
    species_dict[:SO2] = N_gas * y_SO2 * 64.066e-3

    return species_dict
end

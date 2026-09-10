# Erebus redox equilibrium buffers and electron budget accounting (Evans 2012)

"""
Retrieve thermodynamic coefficients for standard solid-oxide redox buffers.

$(SIGNATURES)

Calculates coefficients (A, B, C) for the linear oxygen fugacity buffer equation:
`log10(fO2) = A/T + B + C * (P_bar - 1) / T`
with temperature T in Kelvin and pressure P_bar in bar.

# Arguments
- `buffer`: Buffer symbol (`:IW`, `:QFM`, `:NNO`, `:MH`, `:WM`, `:QIF`)

# Returns
- `Tuple{Float64, Float64, Float64}`: Coefficients `(A, B, C)`

# References
- Frost, B. R. (1991), "Introduction to oxygen fugacity and its petrologic importance",
  Reviews in Mineralogy, 25, 1-9. DOI: 10.1515/9781501508684-004
- Campbell et al. (2009) / O'Neill (1988) for Iron-Wüstite baseline.
"""
@inline function _buffer_coeffs(buffer::Symbol)::Tuple{Float64,Float64,Float64}
    if buffer === :IW
        return (-28164.0, 6.541, 0.0)
    elseif buffer === :QFM
        return (-25096.3, 8.735, 0.110)
    elseif buffer === :NNO
        return (-24930.0, 9.360, 0.046)
    elseif buffer === :MH
        return (-25497.5, 14.330, 0.019)
    elseif buffer === :WM
        return (-32807.0, 13.012, 0.083)
    elseif buffer === :QIF
        return (-29435.7, 7.391, 0.044)
    else
        throw(ArgumentError("Unknown redox buffer: $buffer"))
    end
end

"""
Compute absolute log10(fO2 [bar]) of a named redox buffer at specified (T, P).

$(SIGNATURES)

Supports solid-oxide buffers (:IW, :QFM, :NNO, :MH, :WM, :QIF) and the graphite-CO-CO2
(:CCO) buffer via exact fugacity inversion.

# Arguments
- `buffer`: Redox buffer symbol (:IW, :QFM, :NNO, :MH, :WM, :QIF, :CCO)
- `T_K`: Temperature [K]
- `P_Pa`: Total pressure [Pa]

# Returns
- `Float64`: log10 of oxygen fugacity [bar]

# Raises
- `DomainError`: If T_K <= 0, non-finite, or P_Pa < 0
- `ArgumentError`: If buffer symbol is not recognized
"""
function log10_fo2_of_buffer(buffer::Symbol, T_K::Real, P_Pa::Real)::Float64
    T = Float64(T_K)
    if T <= 0.0 || !isfinite(T)
        throw(DomainError(T, "Temperature must be > 0 and finite"))
    end
    P = Float64(P_Pa)
    if P < 0.0 || !isfinite(P)
        throw(DomainError(P, "Pressure must be >= 0 and finite"))
    end

    if buffer === :CCO
        P_bar = P * 1.0e-5
        if P_bar <= 0.0
            throw(DomainError(P, "Pressure must be > 0 for CCO graphite equilibrium"))
        end
        # Invert C + 0.5 O2 <=> CO and C + O2 <=> CO2 equilibrium where P_bar = f_CO + f_CO2
        log_k_co = 5785.0 / T + 4.545
        log_k_co2 = 20590.0 / T - 0.043
        K_co = 10.0^log_k_co
        K_co2 = 10.0^log_k_co2

        # Quadratic equation in x = sqrt(fO2): K_co2 * x^2 + K_co * x - P_bar = 0
        discriminant = K_co^2 + 4.0 * K_co2 * P_bar
        x = (-K_co + sqrt(discriminant)) / (2.0 * K_co2)
        return 2.0 * log10(x)
    else
        A, B, C = _buffer_coeffs(buffer)
        P_bar = P * 1.0e-5
        return A / T + B + C * (P_bar - 1.0) / T
    end
end

"""
Convert a relative buffer offset to absolute log10(fO2 [bar]).

$(SIGNATURES)

# Arguments
- `delta_val`: Offset relative to the specified buffer in log10 units
- `buffer`: Buffer symbol
- `T_K`: Temperature [K]
- `P_Pa`: Total pressure [Pa]

# Returns
- `Float64`: Absolute log10(fO2 [bar])
"""
@inline function delta_buffer_to_log10_fo2(
    delta_val::Real, buffer::Symbol, T_K::Real, P_Pa::Real
)::Float64
    d = Float64(delta_val)
    if !isfinite(d)
        throw(DomainError(d, "delta_val must be finite"))
    end
    return log10_fo2_of_buffer(buffer, T_K, P_Pa) + d
end

"""
Convert absolute log10(fO2 [bar]) to a relative buffer offset.

$(SIGNATURES)

# Arguments
- `log10_fo2`: Absolute log10(fO2 [bar])
- `buffer`: Reference buffer symbol
- `T_K`: Temperature [K]
- `P_Pa`: Total pressure [Pa]

# Returns
- `Float64`: Offset relative to buffer in log10 units
"""
@inline function log10_fo2_to_delta_buffer(
    log10_fo2::Real, buffer::Symbol, T_K::Real, P_Pa::Real
)::Float64
    lfo2 = Float64(log10_fo2)
    if !isfinite(lfo2)
        throw(DomainError(lfo2, "log10_fo2 must be finite"))
    end
    return lfo2 - log10_fo2_of_buffer(buffer, T_K, P_Pa)
end

"""
Convert a relative redox offset from one buffer frame to another.

$(SIGNATURES)

Preserves exact round-trip identity to floating-point precision.

# Arguments
- `val`: Offset in log10 units relative to `from_buf`
- `from_buf`: Source buffer symbol
- `to_buf`: Destination buffer symbol
- `T_K`: Temperature [K]
- `P_Pa`: Total pressure [Pa]

# Returns
- `Float64`: Offset in log10 units relative to `to_buf`
"""
@inline function convert_redox_buffer(
    val::Real, from_buf::Symbol, to_buf::Symbol, T_K::Real, P_Pa::Real
)::Float64
    lfo2 = delta_buffer_to_log10_fo2(val, from_buf, T_K, P_Pa)
    return log10_fo2_to_delta_buffer(lfo2, to_buf, T_K, P_Pa)
end

"""
Identify the local controlling redox buffer and regime for a simulation cell.

$(SIGNATURES)

# Arguments
- `w_metal_fe`: Mass fraction of metallic iron (Fe0)
- `w_graphite`: Mass fraction of elemental carbon (graphite/carbide)
- `x_ferric`: Molar fraction of ferric iron (Fe3+ / Sigma Fe)

# Keyword Arguments
- `tol`: Detection threshold for phase presence (default: 1.0e-6)

# Returns
- `Tuple{Symbol, Symbol}`: (controlling_buffer, regime_description)
"""
function local_controlling_buffer(
    w_metal_fe::Real, w_graphite::Real, x_ferric::Real; tol::Real=1.0e-6
)::Tuple{Symbol,Symbol}
    if w_metal_fe > tol
        return (:IW, :metal_saturated)
    elseif w_graphite > tol
        return (:CCO, :graphite_saturated)
    elseif x_ferric > tol
        return (:QFM, :silicate_melt)
    else
        return (:QFM, :gas_ratio_fallback)
    end
end

"""
Redox active component inventory in moles for Evans (2012) electron accounting.

$(SIGNATURES)

# Fields
- `n_Fe0`: Moles of metallic iron (Fe0)
- `n_Fe2`: Moles of ferrous iron (Fe2+)
- `n_Fe3`: Moles of ferric iron (Fe3+)
- `n_H2`: Moles of molecular hydrogen (H2)
- `n_H2O`: Moles of water (H2O)
- `n_C_graphite`: Moles of elemental carbon (C0)
- `n_CO`: Moles of carbon monoxide (CO)
- `n_CO2`: Moles of carbon dioxide (CO2)
- `n_CH4`: Moles of methane (CH4)
- `n_Fe3C`: Moles of cohenite (Fe3C)
- `n_S_sulfide`: Moles of sulfide (S2-)
- `n_S2`: Moles of diatomic sulfur gas (S2)
- `n_SO2`: Moles of sulfur dioxide (SO2)
- `n_SO4`: Moles of sulfate (SO4 2-)
- `n_P_phosphide`: Moles of phosphide (P3-)
- `n_P_phosphate`: Moles of phosphate (P5+)
"""
Base.@kwdef struct RedoxComponents
    n_Fe0::Float64 = 0.0
    n_Fe2::Float64 = 0.0
    n_Fe3::Float64 = 0.0
    n_H2::Float64 = 0.0
    n_H2O::Float64 = 0.0
    n_C_graphite::Float64 = 0.0
    n_CO::Float64 = 0.0
    n_CO2::Float64 = 0.0
    n_CH4::Float64 = 0.0
    n_Fe3C::Float64 = 0.0
    n_S_sulfide::Float64 = 0.0
    n_S2::Float64 = 0.0
    n_SO2::Float64 = 0.0
    n_SO4::Float64 = 0.0
    n_P_phosphide::Float64 = 0.0
    n_P_phosphate::Float64 = 0.0
end

"""
Compute extensive redox budget RB [mol e-] following Evans (2012).

$(SIGNATURES)

Calculates total moles of electrons required to bring the assemblage to a
specified reference state:
`RB = sum_i n_i * nu_i`

# Arguments
- `c`: `RedoxComponents` inventory

# Keyword Arguments
- `reference`: Reference state (`:mantle` for Fe2+, C0, S2-, H+, O2-, P5+;
  `:crust` for Fe3+, C4+, S6+, H+, O2-, P5+)

# Returns
- `Float64`: Extensive redox budget [mol e-]

# References
- Evans, K. A. (2012), "The redox budget of subduction zones",
  Earth-Science Reviews, 113, 11-32. DOI: 10.1016/j.earscirev.2012.03.003
"""
function compute_redox_budget(c::RedoxComponents; reference::Symbol=:mantle)::Float64
    if reference === :mantle
        # Mantle reference: Fe2+ (+2), C0 (0), S2- (-2), H+ (+1), O2- (-2), P5+ (+5)
        # nu_i = z_i - z_ref,i is the electrons added to reach reference state
        return (
            -2.0 * c.n_Fe0 +
            0.0 * c.n_Fe2 +
            1.0 * c.n_Fe3 +
            -2.0 * c.n_H2 +
            0.0 * c.n_H2O +
            0.0 * c.n_C_graphite +
            2.0 * c.n_CO +
            4.0 * c.n_CO2 +
            -4.0 * c.n_CH4 +
            -6.0 * c.n_Fe3C +
            0.0 * c.n_S_sulfide +
            4.0 * c.n_S2 +
            6.0 * c.n_SO2 +
            8.0 * c.n_SO4 +
            -8.0 * c.n_P_phosphide +
            0.0 * c.n_P_phosphate
        )
    elseif reference === :crust
        # Crust reference: Fe3+ (+3), C4+ (+4), S6+ (+6), H+ (+1), O2- (-2), P5+ (+5)
        return (
            -3.0 * c.n_Fe0 +
            -1.0 * c.n_Fe2 +
            0.0 * c.n_Fe3 +
            -2.0 * c.n_H2 +
            0.0 * c.n_H2O +
            -4.0 * c.n_C_graphite +
            -2.0 * c.n_CO +
            0.0 * c.n_CO2 +
            -8.0 * c.n_CH4 +
            -13.0 * c.n_Fe3C +
            -8.0 * c.n_S_sulfide +
            -12.0 * c.n_S2 +
            -2.0 * c.n_SO2 +
            0.0 * c.n_SO4 +
            -8.0 * c.n_P_phosphide +
            0.0 * c.n_P_phosphate
        )
    else
        throw(
            ArgumentError("Unknown reference state: $reference. Must be :mantle or :crust")
        )
    end
end

"""
Compute specific redox budget RB_M [mol e- / kg] following Evans (2012).

$(SIGNATURES)

# Arguments
- `c`: `RedoxComponents` inventory
- `mass_kg`: Total mass of the system [kg]

# Keyword Arguments
- `reference`: Reference state (`:mantle` or `:crust`)

# Returns
- `Float64`: Specific redox budget [mol e- / kg]
"""
function compute_specific_redox_budget(
    c::RedoxComponents, mass_kg::Real; reference::Symbol=:mantle
)::Float64
    m = Float64(mass_kg)
    if m <= 0.0 || !isfinite(m)
        throw(DomainError(m, "Mass must be > 0 and finite"))
    end
    return compute_redox_budget(c; reference=reference) / m
end

"""
Update redox component inventory during serpentinization reaction.

$(SIGNATURES)

Models the reaction: `3 FeO (Fe2+) + H2O -> Fe3O4 (1 Fe2+ + 2 Fe3+) + H2`.
This reaction exactly conserves the extensive redox budget RB.

# Arguments
- `c`: Initial `RedoxComponents`
- `delta_n_h2o`: Moles of H2O reacted

# Returns
- `RedoxComponents`: Updated inventory
"""
function serpentinize_redox_budget(c::RedoxComponents, delta_n_h2o::Real)::RedoxComponents
    dn = Float64(delta_n_h2o)
    return RedoxComponents(;
        n_Fe0=c.n_Fe0,
        n_Fe2=c.n_Fe2 - 2.0 * dn,
        n_Fe3=c.n_Fe3 + 2.0 * dn,
        n_H2=c.n_H2 + dn,
        n_H2O=c.n_H2O - dn,
        n_C_graphite=c.n_C_graphite,
        n_CO=c.n_CO,
        n_CO2=c.n_CO2,
        n_CH4=c.n_CH4,
        n_Fe3C=c.n_Fe3C,
        n_S_sulfide=c.n_S_sulfide,
        n_S2=c.n_S2,
        n_SO2=c.n_SO2,
        n_SO4=c.n_SO4,
        n_P_phosphide=c.n_P_phosphide,
        n_P_phosphate=c.n_P_phosphate,
    )
end

"""
Partition redox components during core segregation.

$(SIGNATURES)

Transfers metallic Fe0, cohenite, and schreibersite into the core reservoir.
Conserves total whole-body electrons: `RB_bulk = RB_mantle + RB_core`.

# Arguments
- `c`: Bulk planetesimal `RedoxComponents`
- `fraction_to_core`: Fraction of metallic phases segregated into the core (in [0, 1])

# Returns
- `Tuple{RedoxComponents, RedoxComponents}`: `(c_mantle, c_core)`
"""
function segregate_core_redox_budget(
    c::RedoxComponents, fraction_to_core::Real
)::Tuple{RedoxComponents,RedoxComponents}
    f = clamp(Float64(fraction_to_core), 0.0, 1.0)

    c_core = RedoxComponents(;
        n_Fe0=c.n_Fe0 * f, n_Fe3C=c.n_Fe3C * f, n_P_phosphide=c.n_P_phosphide * f
    )

    c_mantle = RedoxComponents(;
        n_Fe0=c.n_Fe0 * (1.0 - f),
        n_Fe2=c.n_Fe2,
        n_Fe3=c.n_Fe3,
        n_H2=c.n_H2,
        n_H2O=c.n_H2O,
        n_C_graphite=c.n_C_graphite,
        n_CO=c.n_CO,
        n_CO2=c.n_CO2,
        n_CH4=c.n_CH4,
        n_Fe3C=c.n_Fe3C * (1.0 - f),
        n_S_sulfide=c.n_S_sulfide,
        n_S2=c.n_S2,
        n_SO2=c.n_SO2,
        n_SO4=c.n_SO4,
        n_P_phosphide=c.n_P_phosphide * (1.0 - f),
        n_P_phosphate=c.n_P_phosphate,
    )

    return (c_mantle, c_core)
end

"""
Subtract vented gas volatiles from residual rock redox inventory.

$(SIGNATURES)

Conserves total electrons across rock and gas reservoirs:
`RB_rock_initial = RB_rock_final + RB_gas_vented`.

# Arguments
- `c_rock`: Initial rock `RedoxComponents`
- `c_vent`: Vented gas `RedoxComponents`

# Returns
- `RedoxComponents`: Updated rock inventory
"""
function vent_gas_redox_budget(
    c_rock::RedoxComponents, c_vent::RedoxComponents
)::RedoxComponents
    return RedoxComponents(;
        n_Fe0=c_rock.n_Fe0,
        n_Fe2=c_rock.n_Fe2,
        n_Fe3=c_rock.n_Fe3,
        n_H2=max(0.0, c_rock.n_H2 - c_vent.n_H2),
        n_H2O=max(0.0, c_rock.n_H2O - c_vent.n_H2O),
        n_C_graphite=c_rock.n_C_graphite,
        n_CO=max(0.0, c_rock.n_CO - c_vent.n_CO),
        n_CO2=max(0.0, c_rock.n_CO2 - c_vent.n_CO2),
        n_CH4=max(0.0, c_rock.n_CH4 - c_vent.n_CH4),
        n_Fe3C=c_rock.n_Fe3C,
        n_S_sulfide=c_rock.n_S_sulfide,
        n_S2=max(0.0, c_rock.n_S2 - c_vent.n_S2),
        n_SO2=max(0.0, c_rock.n_SO2 - c_vent.n_SO2),
        n_SO4=c_rock.n_SO4,
        n_P_phosphide=c_rock.n_P_phosphide,
        n_P_phosphate=c_rock.n_P_phosphate,
    )
end

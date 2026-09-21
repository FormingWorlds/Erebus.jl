"""
Compute metal-silicate partition coefficient D_i = C_metal / C_silicate for volatile species i in {:H, :C, :N, :S}.

$(SIGNATURES)

# Arguments
- `species::Symbol`: Volatile element (`:H`, `:C`, `:N`, `:S`)
- `T::Real`: Temperature [K]
- `P::Real`: Pressure [Pa]
- `ΔIW::Real`: Oxygen fugacity relative to Iron-Wüstite buffer [log10 units]
- `w_S::Real`: Sulfur mass fraction in metallic alloy [0.0, 1.0]

# Keyword Arguments
- `model::Symbol`: Parameterization model. Supported:
  - `:constant`: Uses fixed partition coefficient `D_const`.
  - For `:C`: `:grewal2019`, `:fischer2020`.
  - For `:N`: `:grewal2019`.
  - For `:H`: `:clesi2018`.
  - For `:S`: `:boujibar2014`.
- `D_const::Real`: Fixed partition coefficient value (default: 1.0)
- `D_min::Real`: Numerical lower floor (default: 1.0e-4)
- `D_max::Real`: Numerical upper ceiling (default: 1.0e5)

# Returns
- `D_val::Float64`: Metal-silicate partition coefficient [-].

# Raises
- `DomainError`: If `T <= 0.0`, `P < 0.0`, `w_S < 0.0 || w_S > 1.0`, or inputs are non-finite.
- `ArgumentError`: If `species` or `model` is unsupported.
"""
function compute_metal_silicate_partition_coefficient(
    species::Symbol,
    T::Real,
    P::Real,
    ΔIW::Real,
    w_S::Real;
    model::Symbol=:default,
    D_const::Real=1.0,
    D_min::Real=1.0e-4,
    D_max::Real=1.0e5,
)::Float64
    T_val = Float64(T)
    P_val = Float64(P)
    ΔIW_val = Float64(ΔIW)
    w_val = Float64(w_S)
    D_c = Float64(D_const)
    d_min = Float64(D_min)
    d_max = Float64(D_max)

    if T_val <= 0.0 || !isfinite(T_val)
        throw(DomainError(T_val, "Temperature must be positive and finite"))
    end
    if P_val < 0.0 || !isfinite(P_val)
        throw(DomainError(P_val, "Pressure must be non-negative and finite"))
    end
    if !isfinite(ΔIW_val)
        throw(DomainError(ΔIW_val, "Oxygen fugacity ΔIW must be finite"))
    end
    if !(0.0 <= w_val <= 1.0) || !isfinite(w_val)
        throw(DomainError(w_val, "Sulfur mass fraction must be in [0, 1] and finite"))
    end
    if d_min <= 0.0 || !isfinite(d_min)
        throw(DomainError(d_min, "D_min must be positive and finite"))
    end
    if d_max < d_min || !isfinite(d_max)
        throw(DomainError(d_max, "D_max must be >= D_min and finite"))
    end

    if species !== :H && species !== :C && species !== :N && species !== :S
        throw(
            ArgumentError("Unknown volatile species: $species. Supported: :H, :C, :N, :S")
        )
    end

    mod = if model === :default
        if species === :C
            :grewal2019
        elseif species === :N
            :grewal2019
        elseif species === :H
            :clesi2018
        elseif species === :S
            :boujibar2014
        else
            :constant
        end
    else
        model
    end

    if mod === :constant
        return clamp(D_c, d_min, d_max)
    end

    # Fe-S molar conversion for sulfur-alloy interaction terms
    # M_S = 32.065 g/mol, M_Fe = 55.845 g/mol
    n_S = w_val / 32.065
    n_Fe = (1.0 - w_val) / 55.845
    X_S = (n_S + n_Fe) > 0.0 ? n_S / (n_S + n_Fe) : 0.0
    # Guard against singular log(1 - X_S) when alloy approaches pure sulfur
    ln_1_minus_XS = log(max(1.0 - min(X_S, 0.999), 1.0e-6))

    log10_D = if species === :C
        if mod === :grewal2019
            # Grewal et al. (2019, Science Advances 5:eaau3669):
            # Strong siderophile behavior suppressed by dissolved sulfur in metallic melt
            1.80 + 2200.0 / T_val - 1.5e-8 * (P_val / T_val) - 0.25 * ΔIW_val +
            4.2 * ln_1_minus_XS
        elseif mod === :fischer2020
            # Fischer et al. (2020, PNAS 117:8743-8749):
            1.50 + 2500.0 / T_val - 1.2e-8 * (P_val / T_val) - 0.20 * ΔIW_val
        else
            throw(
                ArgumentError(
                    "Unknown carbon partition model: $mod. Supported: :constant, :grewal2019, :fischer2020",
                ),
            )
        end
    elseif species === :N
        if mod === :grewal2019
            # Grewal et al. (2019, GCA 251:87-115; 2019, Sci. Adv. 5:eaau3669):
            # Nitrogen siderophile partitioning is weakly dependent on sulfur compared to carbon
            0.85 + 1200.0 / T_val - 0.25 * ΔIW_val + 0.60 * ln_1_minus_XS
        else
            throw(
                ArgumentError(
                    "Unknown nitrogen partition model: $mod. Supported: :constant, :grewal2019",
                ),
            )
        end
    elseif species === :H
        if mod === :clesi2018
            # Clesi et al. (2018, Science Advances 4:e1701876):
            # Low-pressure planetesimal regime: moderately siderophile to lithophile
            -0.80 + 300.0 / T_val + 5.0e-8 * (P_val / T_val) + 0.05 * ΔIW_val
        else
            throw(
                ArgumentError(
                    "Unknown hydrogen partition model: $mod. Supported: :constant, :clesi2018",
                ),
            )
        end
    elseif species === :S
        if mod === :boujibar2014
            # Boujibar et al. (2014, EPSL 391:42-54):
            # Strong chalcophile/siderophile partitioning of sulfur into liquid metal
            2.80 - 800.0 / T_val + 1.0e-10 * P_val - 0.20 * ΔIW_val
        else
            throw(
                ArgumentError(
                    "Unknown sulfur partition model: $mod. Supported: :constant, :boujibar2014",
                ),
            )
        end
    else
        throw(ArgumentError("Unknown volatile species: $species. Supported: :H, :C, :N, :S"))
    end

    return clamp(10.0^log10_D, d_min, d_max)
end

"""
Compute metal-silicate partition coefficients for H, C, N, and S in a single call.

$(SIGNATURES)

# Arguments
- `T::Real`: Temperature [K]
- `P::Real`: Pressure [Pa]
- `ΔIW::Real`: Oxygen fugacity relative to Iron-Wüstite buffer [log10 units]
- `w_S::Real`: Sulfur mass fraction in metallic alloy [0.0, 1.0]
- `cfg::MetalPartitionConfig`: Metal partition configuration

# Returns
- NamedTuple `(; D_H, D_C, D_N, D_S)`: Partition coefficients [-].
"""
function compute_metal_silicate_partition_coefficients(
    T::Real, P::Real, ΔIW::Real, w_S::Real, cfg::MetalPartitionConfig
)
    D_H = compute_metal_silicate_partition_coefficient(
        :H,
        T,
        P,
        ΔIW,
        w_S;
        model=cfg.model_hydrogen,
        D_const=cfg.D_H_const,
        D_min=cfg.D_min,
        D_max=cfg.D_max,
    )
    D_C = compute_metal_silicate_partition_coefficient(
        :C,
        T,
        P,
        ΔIW,
        w_S;
        model=cfg.model_carbon,
        D_const=cfg.D_C_const,
        D_min=cfg.D_min,
        D_max=cfg.D_max,
    )
    D_N = compute_metal_silicate_partition_coefficient(
        :N,
        T,
        P,
        ΔIW,
        w_S;
        model=cfg.model_nitrogen,
        D_const=cfg.D_N_const,
        D_min=cfg.D_min,
        D_max=cfg.D_max,
    )
    D_S = compute_metal_silicate_partition_coefficient(
        :S,
        T,
        P,
        ΔIW,
        w_S;
        model=cfg.model_sulfur,
        D_const=cfg.D_S_const,
        D_min=cfg.D_min,
        D_max=cfg.D_max,
    )
    return (; D_H, D_C, D_N, D_S)
end

"""
Equilibrate volatile concentrations between molten metallic iron and silicate melt on marker m.

Conserves total elemental mass of H, C, N, and S across the two interacting reservoirs:
    M_i = m_sil_melt * C_i_sil + m_met * C_i_met = const

$(SIGNATURES)

# Arguments
- `m::Integer`: Marker index
- `F_fe::Real`: Metal melt fraction [0, 1]
- `F_melt::Real`: Silicate melt fraction [0, 1]
- `T_val::Real`: Temperature [K]
- `P_val::Real`: Pressure [Pa]
- `ΔIW::Real`: Oxygen fugacity relative to Iron-Wüstite buffer [log10 units]
- `Xfe_bulk::AbstractVector{Float64}`: Marker bulk metal volume fraction [0, 1]
- `Xfem::AbstractVector{Float64}`: Marker molten metal volume fraction [0, 1]
- `XH2Om::AbstractVector{Float64}`: Marker silicate water concentration array [wt%]
- `XCm::Union{Nothing,AbstractVector{Float64}}`: Marker silicate carbon concentration array [ppmw]
- `XNm::Union{Nothing,AbstractVector{Float64}}`: Marker silicate nitrogen concentration array [ppmw]
- `XSm::Union{Nothing,AbstractVector{Float64}}`: Marker silicate sulfur concentration array [ppmw]
- `Xfe_H_m::AbstractVector{Float64}`: Marker metal hydrogen concentration array [ppmw]
- `Xfe_C_m::AbstractVector{Float64}`: Marker metal carbon concentration array [ppmw]
- `Xfe_N_m::AbstractVector{Float64}`: Marker metal nitrogen concentration array [ppmw]
- `Xfe_S_m::AbstractVector{Float64}`: Marker metal sulfur concentration array [ppmw]
- `cfg::MetalPartitionConfig`: Partition configuration

# Keyword Arguments
- `rho_silicate::Real`: Silicate reference density [kg/m^3] (default: 3000.0)
- `rho_metal::Real`: Liquid metal reference density [kg/m^3] (default: 7000.0)
- `equilibration_fraction::Real`: Kinetic equilibration factor in [0, 1] (default: cfg.equilibration_rate)
"""
function equilibrate_metal_silicate_volatiles!(
    m::Integer,
    F_fe::Real,
    F_melt::Real,
    T_val::Real,
    P_val::Real,
    ΔIW::Real,
    Xfe_bulk::AbstractVector{Float64},
    Xfem::AbstractVector{Float64},
    XH2Om::AbstractVector{Float64},
    XCm::Union{Nothing,AbstractVector{Float64}},
    XNm::Union{Nothing,AbstractVector{Float64}},
    XSm::Union{Nothing,AbstractVector{Float64}},
    Xfe_H_m::AbstractVector{Float64},
    Xfe_C_m::AbstractVector{Float64},
    Xfe_N_m::AbstractVector{Float64},
    Xfe_S_m::AbstractVector{Float64},
    cfg::MetalPartitionConfig;
    rho_silicate::Real=3000.0,
    rho_metal::Real=7000.0,
    equilibration_fraction::Real=cfg.equilibration_rate,
)
    if !isfinite(T_val) || T_val <= 0.0
        throw(DomainError(T_val, "Temperature must be positive and finite"))
    end
    if !isfinite(P_val) || P_val < 0.0
        throw(DomainError(P_val, "Pressure must be non-negative and finite"))
    end
    if !isfinite(ΔIW)
        throw(DomainError(ΔIW, "Oxygen fugacity ΔIW must be finite"))
    end

    phi_fe = Xfe_bulk[m]
    phi_sil = max(1.0 - phi_fe, 0.0)
    F_fe_val = Float64(F_fe)
    F_melt_val = Float64(F_melt)
    if phi_fe <= 1.0e-7 || phi_sil <= 1.0e-7 || F_fe_val <= 0.0 || F_melt_val <= 0.0
        return nothing
    end

    alpha_eq = clamp(Float64(equilibration_fraction), 0.0, 1.0)
    if alpha_eq <= 0.0
        return nothing
    end

    # Interacting phase masses per unit marker volume (only molten silicate participates)
    m_met = Xfem[m] * max(Float64(rho_metal), 100.0)
    m_sil_melt = phi_sil * F_melt_val * max(Float64(rho_silicate), 100.0)
    if m_met <= 0.0 || m_sil_melt <= 0.0
        return nothing
    end

    T_m = Float64(T_val)
    P_m = Float64(P_val)
    ΔIW_m = Float64(ΔIW)

    # Current metal sulfur mass fraction
    w_S = clamp(Xfe_S_m[m] * 1.0e-6, 0.0, 0.365)

    # 1. Carbon equilibration (graphite saturation ceiling in liquid Fe: ~7 wt% = 70,000 ppmw)
    if XCm !== nothing
        D_C = compute_metal_silicate_partition_coefficient(
            :C,
            T_m,
            P_m,
            ΔIW_m,
            w_S;
            model=cfg.model_carbon,
            D_const=cfg.D_C_const,
            D_min=cfg.D_min,
            D_max=cfg.D_max,
        )
        C_sil_bulk = XCm[m]
        C_met = Xfe_C_m[m]
        C_sil_melt = C_sil_bulk / F_melt_val
        M_tot = m_sil_melt * C_sil_melt + m_met * C_met
        denom = m_sil_melt + D_C * m_met
        if denom > 0.0
            C_sil_melt_eq = M_tot / denom
            C_met_eq = D_C * C_sil_melt_eq
            C_met_C_max = min(Float64(cfg.D_max), 7.0e4)
            if C_met_eq > C_met_C_max
                C_met_eq = C_met_C_max
                C_sil_melt_eq = max(0.0, (M_tot - m_met * C_met_eq) / m_sil_melt)
            end
            dC_sil_melt = alpha_eq * (C_sil_melt_eq - C_sil_melt)
            dC_sil_melt = max(dC_sil_melt, -C_sil_melt)
            dC_met = -dC_sil_melt * (m_sil_melt / m_met)
            if C_met + dC_met < 0.0
                dC_met = -C_met
                dC_sil_melt = -dC_met * (m_met / m_sil_melt)
            end
            XCm[m] = max(0.0, C_sil_bulk + dC_sil_melt * F_melt_val)
            Xfe_C_m[m] = clamp(C_met + dC_met, 0.0, C_met_C_max)
        end
    end

    # 2. Nitrogen equilibration (nitrogen saturation ceiling in liquid Fe: ~4 wt% = 40,000 ppmw)
    if XNm !== nothing
        D_N = compute_metal_silicate_partition_coefficient(
            :N,
            T_m,
            P_m,
            ΔIW_m,
            w_S;
            model=cfg.model_nitrogen,
            D_const=cfg.D_N_const,
            D_min=cfg.D_min,
            D_max=cfg.D_max,
        )
        C_sil_bulk = XNm[m]
        C_met = Xfe_N_m[m]
        C_sil_melt = C_sil_bulk / F_melt_val
        M_tot = m_sil_melt * C_sil_melt + m_met * C_met
        denom = m_sil_melt + D_N * m_met
        if denom > 0.0
            C_sil_melt_eq = M_tot / denom
            C_met_eq = D_N * C_sil_melt_eq
            C_met_N_max = min(Float64(cfg.D_max), 4.0e4)
            if C_met_eq > C_met_N_max
                C_met_eq = C_met_N_max
                C_sil_melt_eq = max(0.0, (M_tot - m_met * C_met_eq) / m_sil_melt)
            end
            dC_sil_melt = alpha_eq * (C_sil_melt_eq - C_sil_melt)
            dC_sil_melt = max(dC_sil_melt, -C_sil_melt)
            dC_met = -dC_sil_melt * (m_sil_melt / m_met)
            if C_met + dC_met < 0.0
                dC_met = -C_met
                dC_sil_melt = -dC_met * (m_met / m_sil_melt)
            end
            XNm[m] = max(0.0, C_sil_bulk + dC_sil_melt * F_melt_val)
            Xfe_N_m[m] = clamp(C_met + dC_met, 0.0, C_met_N_max)
        end
    end

    # 3. Sulfur equilibration (troilite/FeS saturation ceiling: ~36.5 wt% = 365,000 ppmw)
    if XSm !== nothing
        D_S = compute_metal_silicate_partition_coefficient(
            :S,
            T_m,
            P_m,
            ΔIW_m,
            w_S;
            model=cfg.model_sulfur,
            D_const=cfg.D_S_const,
            D_min=cfg.D_min,
            D_max=cfg.D_max,
        )
        C_sil_bulk = XSm[m]
        C_met = Xfe_S_m[m]
        C_sil_melt = C_sil_bulk / F_melt_val
        M_tot = m_sil_melt * C_sil_melt + m_met * C_met
        denom = m_sil_melt + D_S * m_met
        if denom > 0.0
            C_sil_melt_eq = M_tot / denom
            C_met_eq = D_S * C_sil_melt_eq
            C_met_S_max = min(Float64(cfg.D_max), 3.65e5)
            if C_met_eq > C_met_S_max
                C_met_eq = C_met_S_max
                C_sil_melt_eq = max(0.0, (M_tot - m_met * C_met_eq) / m_sil_melt)
            end
            dC_sil_melt = alpha_eq * (C_sil_melt_eq - C_sil_melt)
            dC_sil_melt = max(dC_sil_melt, -C_sil_melt)
            dC_met = -dC_sil_melt * (m_sil_melt / m_met)
            if C_met + dC_met < 0.0
                dC_met = -C_met
                dC_sil_melt = -dC_met * (m_met / m_sil_melt)
            end
            XSm[m] = max(0.0, C_sil_bulk + dC_sil_melt * F_melt_val)
            Xfe_S_m[m] = clamp(C_met + dC_met, 0.0, C_met_S_max)
        end
    end

    # 4. Hydrogen equilibration (stoichiometric conversion: H2O [wt%] <-> H [ppmw])
    # (2 * 1.00794 / 18.01528) * 1.0e4 = 1118.9834407236524
    f_H = (2.0 * 1.00794 / 18.01528) * 1.0e4
    D_H = compute_metal_silicate_partition_coefficient(
        :H,
        T_m,
        P_m,
        ΔIW_m,
        w_S;
        model=cfg.model_hydrogen,
        D_const=cfg.D_H_const,
        D_min=cfg.D_min,
        D_max=cfg.D_max,
    )
    C_sil_bulk_H = XH2Om[m] * f_H
    C_met_H = Xfe_H_m[m]
    C_sil_melt_H = C_sil_bulk_H / F_melt_val
    M_tot_H = m_sil_melt * C_sil_melt_H + m_met * C_met_H
    denom_H = m_sil_melt + D_H * m_met
    if denom_H > 0.0
        C_sil_melt_H_eq = M_tot_H / denom_H
        C_met_H_eq = D_H * C_sil_melt_H_eq
        C_met_H_max = min(Float64(cfg.D_max), 1.0e4)
        if C_met_H_eq > C_met_H_max
            C_met_H_eq = C_met_H_max
            C_sil_melt_H_eq = max(0.0, (M_tot_H - m_met * C_met_H_eq) / m_sil_melt)
        end
        dC_sil_melt = alpha_eq * (C_sil_melt_H_eq - C_sil_melt_H)
        dC_sil_melt = max(dC_sil_melt, -C_sil_melt_H)
        dC_met = -dC_sil_melt * (m_sil_melt / m_met)
        if C_met_H + dC_met < 0.0
            dC_met = -C_met_H
            dC_sil_melt = -dC_met * (m_met / m_sil_melt)
        end
        new_C_sil_bulk_H = max(0.0, C_sil_bulk_H + dC_sil_melt * F_melt_val)
        XH2Om[m] = clamp(new_C_sil_bulk_H / f_H, 0.0, 100.0)
        Xfe_H_m[m] = clamp(C_met_H + dC_met, 0.0, C_met_H_max)
    end

    return nothing
end

"""
Compute integrated core mass and volatile budgets for comparison with magmatic iron meteorites.

$(SIGNATURES)

# Arguments
- `xm::AbstractVector{Float64}`: Marker x-coordinates [m]
- `ym::AbstractVector{Float64}`: Marker y-coordinates [m]
- `tm::AbstractVector{<:Integer}`: Marker type array
- `Xfe_bulk::AbstractVector{Float64}`: Marker bulk metal volume fraction [0, 1]
- `Xfe_H_m::Union{Nothing,AbstractVector{Float64}}`: Metal hydrogen array [ppmw]
- `Xfe_C_m::Union{Nothing,AbstractVector{Float64}}`: Metal carbon array [ppmw]
- `Xfe_N_m::Union{Nothing,AbstractVector{Float64}}`: Metal nitrogen array [ppmw]
- `Xfe_S_m::Union{Nothing,AbstractVector{Float64}}`: Metal sulfur array [ppmw]
- `marknum::Integer`: Marker count

# Keyword Arguments
- `xcenter::Real`: Planet center x [m] (default: 70000.0)
- `ycenter::Real`: Planet center y [m] (default: 70000.0)
- `rplanet::Real`: Planet radius [m] (default: 50000.0)
- `rho_metal::Real`: Metal density [kg/m^3] (default: 7000.0)
- `core_radius_fraction::Real`: Fractional radius defining central core region (default: 0.5)
- `phi_core_threshold::Real`: Metal volume fraction threshold for core membership (default: 0.40)
- `V_marker::Union{Nothing,Real}`: Explicit marker volume [m³] (default: derived from planetary volume)
- `use_3d_volume::Bool`: If true (default), use 3D spherical equivalent volume (4/3 π R³); if false, use 2D area (π R²)

# Returns
- NamedTuple containing:
  - `M_core_metal`: Total segregated core metal mass [kg]
  - `M_core_H`: Integrated core hydrogen mass [kg]
  - `M_core_C`: Integrated core carbon mass [kg]
  - `M_core_N`: Integrated core nitrogen mass [kg]
  - `M_core_S`: Integrated core sulfur mass [kg]
  - `w_core_H_ppm`: Core hydrogen concentration [ppmw]
  - `w_core_C_ppm`: Core carbon concentration [ppmw]
  - `w_core_N_ppm`: Core nitrogen concentration [ppmw]
  - `w_core_S_wtpct`: Core sulfur concentration [wt%]
  - `w_core_S_ppm`: Core sulfur concentration [ppmw]
  - `M_total_metal`: Total metal mass across entire planet [kg]
  - `M_total_H_met`: Total metal-hosted H mass across planet [kg]
  - `M_total_C_met`: Total metal-hosted C mass across planet [kg]
  - `M_total_N_met`: Total metal-hosted N mass across planet [kg]
  - `M_total_S_met`: Total metal-hosted S mass across planet [kg]
"""
function compute_core_volatile_budgets(
    xm::AbstractVector{Float64},
    ym::AbstractVector{Float64},
    tm::AbstractVector{<:Integer},
    Xfe_bulk::AbstractVector{Float64},
    Xfe_H_m::Union{Nothing,AbstractVector{Float64}},
    Xfe_C_m::Union{Nothing,AbstractVector{Float64}},
    Xfe_N_m::Union{Nothing,AbstractVector{Float64}},
    Xfe_S_m::Union{Nothing,AbstractVector{Float64}},
    marknum::Integer;
    xcenter::Real=70000.0,
    ycenter::Real=70000.0,
    rplanet::Real=50000.0,
    rho_metal::Real=7000.0,
    core_radius_fraction::Real=0.5,
    phi_core_threshold::Real=0.40,
    V_marker::Union{Nothing,Real}=nothing,
    use_3d_volume::Bool=true,
)
    M_core_metal = 0.0
    M_core_H = 0.0
    M_core_C = 0.0
    M_core_N = 0.0
    M_core_S = 0.0

    M_total_metal = 0.0
    M_total_H_met = 0.0
    M_total_C_met = 0.0
    M_total_N_met = 0.0
    M_total_S_met = 0.0

    rc_cut = Float64(rplanet) * clamp(Float64(core_radius_fraction), 0.0, 1.0)
    phi_cut = clamp(Float64(phi_core_threshold), 0.0, 1.0)
    rho_m = max(Float64(rho_metal), 100.0)

    N_planet = 0
    @inbounds for m in 1:marknum
        if tm[m] < 3
            dx = xm[m] - xcenter
            dy = ym[m] - ycenter
            if sqrt(dx^2 + dy^2) <= rplanet
                N_planet += 1
            end
        end
    end

    r_p = Float64(rplanet)
    V_tot = use_3d_volume ? (4.0 / 3.0) * pi * r_p^3 : pi * r_p^2
    V_m = if V_marker !== nothing
        Float64(V_marker)
    elseif N_planet > 0
        V_tot / N_planet
    else
        1.0
    end

    @inbounds for m in 1:marknum
        if tm[m] < 3
            dx = xm[m] - xcenter
            dy = ym[m] - ycenter
            rmark = sqrt(dx^2 + dy^2)
            if rmark <= rplanet
                fe_frac = Xfe_bulk[m]
                if fe_frac > 0.0
                    dM_fe = fe_frac * rho_m * V_m
                    dH = Xfe_H_m !== nothing ? dM_fe * (Xfe_H_m[m] * 1.0e-6) : 0.0
                    dC = Xfe_C_m !== nothing ? dM_fe * (Xfe_C_m[m] * 1.0e-6) : 0.0
                    dN = Xfe_N_m !== nothing ? dM_fe * (Xfe_N_m[m] * 1.0e-6) : 0.0
                    dS = Xfe_S_m !== nothing ? dM_fe * (Xfe_S_m[m] * 1.0e-6) : 0.0

                    M_total_metal += dM_fe
                    M_total_H_met += dH
                    M_total_C_met += dC
                    M_total_N_met += dN
                    M_total_S_met += dS

                    if rmark <= rc_cut || fe_frac >= phi_cut
                        M_core_metal += dM_fe
                        M_core_H += dH
                        M_core_C += dC
                        M_core_N += dN
                        M_core_S += dS
                    end
                end
            end
        end
    end

    w_core_H_ppm = M_core_metal > 0.0 ? (M_core_H / M_core_metal) * 1.0e6 : 0.0
    w_core_C_ppm = M_core_metal > 0.0 ? (M_core_C / M_core_metal) * 1.0e6 : 0.0
    w_core_N_ppm = M_core_metal > 0.0 ? (M_core_N / M_core_metal) * 1.0e6 : 0.0
    w_core_S_ppm = M_core_metal > 0.0 ? (M_core_S / M_core_metal) * 1.0e6 : 0.0
    w_core_S_wtpct = w_core_S_ppm * 1.0e-4

    return (;
        M_core_metal,
        M_core_H,
        M_core_C,
        M_core_N,
        M_core_S,
        w_core_H_ppm,
        w_core_C_ppm,
        w_core_N_ppm,
        w_core_S_wtpct,
        w_core_S_ppm,
        M_total_metal,
        M_total_H_met,
        M_total_C_met,
        M_total_N_met,
        M_total_S_met,
    )
end

"""
    compute_troilite_stoichiometry(w_S::Real)

Compute stoichiometric conversion of sulfur into troilite (FeS).

Parameters
----------
- `w_S::Real`: Mass fraction of sulfur in the metallic alloy [-].

Returns
-------
- `(w_troilite, w_fe_consumed)::Tuple{Float64, Float64}`: Mass fraction of troilite
  formed and iron consumed per unit mass of metallic alloy [-].

Raises
------
- `DomainError`: If `w_S` is negative or non-finite.

Notes
-----
Molar masses: S = 32.065 g/mol, Fe = 55.845 g/mol, FeS = 87.910 g/mol.
"""
function compute_troilite_stoichiometry(w_S::Real)
    (0.0 <= w_S <= 1.0 && isfinite(w_S)) ||
        throw(DomainError(w_S, "w_S must be in [0, 1] and finite"))
    w_S_f = Float64(w_S)
    f_troilite = 87.910 / 32.065
    f_fe = 55.845 / 32.065
    S_fe_limit = (1.0 - w_S_f) / f_fe
    S_troilite = min(w_S_f, S_fe_limit)
    w_troilite = S_troilite * f_troilite
    w_fe_consumed = S_troilite * f_fe
    return (w_troilite, w_fe_consumed)
end

"""
    compute_schreibersite_stoichiometry(w_P::Real; ni_frac::Real=0.25)

Compute stoichiometric conversion of phosphorus into schreibersite ((Fe,Ni)3P).

Parameters
----------
- `w_P::Real`: Mass fraction of phosphorus in the metallic alloy [-].
- `ni_frac::Real`: Molar nickel fraction in the metal matrix Ni/(Fe+Ni) [-] (default: 0.25).

Returns
-------
- `(w_schreibersite, w_metal_consumed)::Tuple{Float64, Float64}`: Mass fraction of schreibersite
  formed and metal (Fe+Ni) consumed per unit mass of metallic alloy [-].

Raises
------
- `DomainError`: If `w_P` is not in [0, 1] or non-finite, or `ni_frac` is not in [0, 1].

Notes
-----
Molar masses: P = 30.97376 g/mol, Fe = 55.845 g/mol, Ni = 58.6934 g/mol.
"""
function compute_schreibersite_stoichiometry(w_P::Real; ni_frac::Real=0.25)
    (0.0 <= w_P <= 1.0 && isfinite(w_P)) ||
        throw(DomainError(w_P, "w_P must be in [0, 1] and finite"))
    (0.0 <= ni_frac <= 1.0 && isfinite(ni_frac)) ||
        throw(DomainError(ni_frac, "ni_frac must be in [0, 1]"))
    w_P_f = Float64(w_P)
    x_ni = Float64(ni_frac)
    M_metal_avg = (1.0 - x_ni) * 55.845 + x_ni * 58.6934
    M_P = 30.97376
    M_schreib = 3.0 * M_metal_avg + M_P
    f_schreib = M_schreib / M_P
    f_metal = (3.0 * M_metal_avg) / M_P
    P_metal_limit = (1.0 - w_P_f) / f_metal
    P_schreib = min(w_P_f, P_metal_limit)
    w_schreib = P_schreib * f_schreib
    w_metal_consumed = P_schreib * f_metal
    return (w_schreib, w_metal_consumed)
end

"""
    compute_cohenite_graphite_stoichiometry(w_C::Real; carbide_max::Real=0.0667)

Compute stoichiometric allocation of carbon into cohenite (Fe3C) and crystalline graphite (C).

Parameters
----------
- `w_C::Real`: Mass fraction of carbon in the metallic alloy [-].
- `carbide_max::Real`: Maximum carbon mass fraction accommodated in carbide [-] (default: 0.0667).

Returns
-------
- `(w_cohenite, w_graphite, w_fe_consumed)::Tuple{Float64, Float64, Float64}`: Mass fraction of
  cohenite, graphite, and iron consumed per unit mass of metallic alloy [-].

Raises
------
- `DomainError`: If `w_C` is negative or non-finite, or `carbide_max` is not in (0, 1].

Notes
-----
Molar masses: C = 12.011 g/mol, Fe = 55.845 g/mol. Stoichiometric cohenite factor ~ 14.948.
When carbon exceeds carbide_max, cohenite saturates and excess carbon precipitates as graphite.
"""
function compute_cohenite_graphite_stoichiometry(w_C::Real; carbide_max::Real=0.0667)
    (0.0 <= w_C <= 1.0 && isfinite(w_C)) ||
        throw(DomainError(w_C, "w_C must be in [0, 1] and finite"))
    (0.0 < carbide_max <= 1.0 && isfinite(carbide_max)) ||
        throw(DomainError(carbide_max, "carbide_max must be in (0, 1]"))
    w_C_f = Float64(w_C)
    c_max = Float64(carbide_max)
    f_fe = (3.0 * 55.845) / 12.011
    f_cohenite = f_fe + 1.0
    C_fe_limit = (1.0 - w_C_f) / f_fe
    C_carbide = min(w_C_f, c_max, C_fe_limit)
    w_cohenite = C_carbide * f_cohenite
    w_graphite = w_C_f - C_carbide
    w_fe_consumed = C_carbide * f_fe
    return (w_cohenite, w_graphite, w_fe_consumed)
end

"""
    compute_nitride_stoichiometry(w_N::Real; mode::Symbol=:roaldite)

Compute stoichiometric conversion of nitrogen into nitride minerals.

Parameters
----------
- `w_N::Real`: Mass fraction of nitrogen in the metallic alloy [-].
- `mode::Symbol`: Nitride mineral model (`:roaldite` for Fe4N, `:carlsbergite` for CrN, `:osbornite` for TiN).

Returns
-------
- `(w_nitride, w_metal_consumed)::Tuple{Float64, Float64}`: Mass fraction of nitride
  formed and metal consumed per unit mass of metallic alloy [-].

Raises
------
- `DomainError`: If `w_N` is not in [0, 1] or non-finite.
- `ArgumentError`: If `mode` is not one of `:roaldite`, `:carlsbergite`, or `:osbornite`.

Notes
-----
Molar masses: N = 14.007 g/mol, Fe = 55.845 g/mol, Cr = 51.996 g/mol, Ti = 47.867 g/mol.
"""
function compute_nitride_stoichiometry(w_N::Real; mode::Symbol=:roaldite)
    (0.0 <= w_N <= 1.0 && isfinite(w_N)) ||
        throw(DomainError(w_N, "w_N must be in [0, 1] and finite"))
    w_N_f = Float64(w_N)
    M_N = 14.007
    f_nitride, f_metal = if mode === :roaldite
        M_Fe = 55.845
        (4.0 * M_Fe + M_N) / M_N, (4.0 * M_Fe) / M_N
    elseif mode === :carlsbergite
        M_Cr = 51.996
        (M_Cr + M_N) / M_N, M_Cr / M_N
    elseif mode === :osbornite
        M_Ti = 47.867
        (M_Ti + M_N) / M_N, M_Ti / M_N
    else
        throw(
            ArgumentError(
                "Unknown nitride_mode: :$mode. Expected :roaldite, :carlsbergite, or :osbornite",
            ),
        )
    end
    N_metal_limit = (1.0 - w_N_f) / f_metal
    N_nitride = min(w_N_f, N_metal_limit)
    w_nitride = N_nitride * f_nitride
    w_metal_consumed = N_nitride * f_metal
    return (w_nitride, w_metal_consumed)
end

"""
    compute_normative_mineral_assemblage(
        T::Real, w_S::Real, w_C::Real, w_N::Real, w_P::Real, cfg::PhaseTrackingConfig
    )

Compute temperature-dependent normative accessory mineral assemblage and eutectic melt fraction.

Parameters
----------
- `T::Real`: Local temperature [K].
- `w_S::Real`: Sulfur mass fraction in metallic alloy [-].
- `w_C::Real`: Carbon mass fraction in metallic alloy [-].
- `w_N::Real`: Nitrogen mass fraction in metallic alloy [-].
- `w_P::Real`: Phosphorus mass fraction in metallic alloy [-].
- `cfg::PhaseTrackingConfig`: Phase tracking configuration struct.

Returns
-------
- Named tuple with fields:
  - `F_solid`: Solid metal fraction in [0, 1] [-].
  - `F_liquid`: Liquid metal fraction in [0, 1] [-].
  - `w_troilite`: Troilite mass fraction in metallic system [-].
  - `w_schreibersite`: Schreibersite mass fraction in metallic system [-].
  - `w_cohenite`: Cohenite mass fraction in metallic system [-].
  - `w_graphite`: Graphite mass fraction in metallic system [-].
  - `w_nitride`: Nitride mass fraction in metallic system [-].
  - `w_metal_matrix`: Solid Fe-Ni metal matrix mass fraction in metallic system [-].
  - `w_liquid_alloy`: Molten Fe-FeS liquid alloy mass fraction in metallic system [-].

Raises
------
- `DomainError`: If `T`, `w_S`, `w_C`, `w_N`, or `w_P` are negative or non-finite,
  if `w_S + w_C + w_N + w_P > 1.0`, or if `cfg.dT_transition` is not strictly positive.

Notes
-----
Sub-eutectic mineral phases dissolve continuously across the eutectic transition interval
[T_eutectic, T_eutectic + dT_transition]. Total phase mass fractions strictly sum to 1.0.
"""
function compute_normative_mineral_assemblage(
    T::Real, w_S::Real, w_C::Real, w_N::Real, w_P::Real, cfg::PhaseTrackingConfig
)
    (T >= 0.0 && isfinite(T)) ||
        throw(DomainError(T, "Temperature must be non-negative and finite"))
    (0.0 <= w_S <= 1.0 && isfinite(w_S)) ||
        throw(DomainError(w_S, "w_S must be in [0, 1] and finite"))
    (0.0 <= w_C <= 1.0 && isfinite(w_C)) ||
        throw(DomainError(w_C, "w_C must be in [0, 1] and finite"))
    (0.0 <= w_N <= 1.0 && isfinite(w_N)) ||
        throw(DomainError(w_N, "w_N must be in [0, 1] and finite"))
    (0.0 <= w_P <= 1.0 && isfinite(w_P)) ||
        throw(DomainError(w_P, "w_P must be in [0, 1] and finite"))
    w_S_f = Float64(w_S)
    w_C_f = Float64(w_C)
    w_N_f = Float64(w_N)
    w_P_f = Float64(w_P)
    w_volatiles = w_S_f + w_C_f + w_N_f + w_P_f
    (w_volatiles <= 1.0) || throw(
        DomainError(w_volatiles, "Sum of volatile mass fractions must not exceed 1.0")
    )
    (cfg.dT_transition > 0.0 && isfinite(cfg.dT_transition)) || throw(
        DomainError(
            cfg.dT_transition, "dT_transition must be strictly positive and finite"
        ),
    )

    T_f = Float64(T)
    T_eut = cfg.T_eutectic
    dT = cfg.dT_transition

    F_solid = clamp(1.0 - (T_f - T_eut) / dT, 0.0, 1.0)
    F_liquid = 1.0 - F_solid

    # Available metallic iron-nickel pool for mineral formation
    w_metal_avail = 1.0 - w_volatiles

    # 1. Troilite (FeS): sulfide has highest affinity for metallic iron
    f_troilite = 87.910 / 32.065
    f_fe_S = 55.845 / 32.065
    S_troilite = min(w_S_f, w_metal_avail / f_fe_S)
    w_troilite_0 = S_troilite * f_troilite
    w_metal_avail = max(0.0, w_metal_avail - S_troilite * f_fe_S)

    # 2. Schreibersite ((Fe,Ni)3P)
    x_ni = Float64(cfg.schreibersite_ni_frac)
    M_metal_avg = (1.0 - x_ni) * 55.845 + x_ni * 58.6934
    M_P = 30.97376
    f_schreib = (3.0 * M_metal_avg + M_P) / M_P
    f_metal_P = (3.0 * M_metal_avg) / M_P
    P_schreib = min(w_P_f, w_metal_avail / f_metal_P)
    w_schreib_0 = P_schreib * f_schreib
    w_metal_avail = max(0.0, w_metal_avail - P_schreib * f_metal_P)

    # 3. Nitride
    M_N = 14.007
    f_nitride, f_metal_N = if cfg.nitride_mode === :roaldite
        M_Fe = 55.845
        (4.0 * M_Fe + M_N) / M_N, (4.0 * M_Fe) / M_N
    elseif cfg.nitride_mode === :carlsbergite
        M_Cr = 51.996
        (M_Cr + M_N) / M_N, M_Cr / M_N
    elseif cfg.nitride_mode === :osbornite
        M_Ti = 47.867
        (M_Ti + M_N) / M_N, M_Ti / M_N
    else
        throw(
            ArgumentError(
                "Unknown nitride_mode: :$(cfg.nitride_mode). Expected :roaldite, :carlsbergite, or :osbornite",
            ),
        )
    end
    N_nit = min(w_N_f, w_metal_avail / f_metal_N)
    w_nit_0 = N_nit * f_nitride
    w_metal_avail = max(0.0, w_metal_avail - N_nit * f_metal_N)

    # 4. Cohenite (Fe3C) and crystalline Graphite (C):
    # Cohenite forms up to carbide saturation and available iron; excess carbon precipitates as graphite
    f_fe_C = (3.0 * 55.845) / 12.011
    f_cohenite = f_fe_C + 1.0
    c_max = Float64(cfg.cohenite_carbide_max)
    C_fe_limit = w_metal_avail / f_fe_C
    C_carbide = min(w_C_f, c_max, C_fe_limit)
    w_coh_0 = C_carbide * f_cohenite
    w_gra_0 = w_C_f - C_carbide
    w_metal_avail = max(0.0, w_metal_avail - C_carbide * f_fe_C)

    # Residual metallic iron-nickel matrix
    w_metal_matrix_0 = w_metal_avail

    w_troilite = F_solid * w_troilite_0
    w_schreibersite = F_solid * w_schreib_0
    w_cohenite = F_solid * w_coh_0
    w_graphite = F_solid * w_gra_0
    w_nitride = F_solid * w_nit_0
    w_metal_matrix = F_solid * w_metal_matrix_0
    w_liquid_alloy = F_liquid

    return (;
        F_solid,
        F_liquid,
        w_troilite,
        w_schreibersite,
        w_cohenite,
        w_graphite,
        w_nitride,
        w_metal_matrix,
        w_liquid_alloy,
    )
end

"""
    compute_regional_mineral_modes(
        xm, ym, tm, tkm, Xfe_bulk, Xfe_S_m, Xfe_C_m, Xfe_N_m, marknum;
        cfg::PhaseTrackingConfig=PhaseTrackingConfig(),
        rplanet::Real=50000.0,
        xcenter::Real=70000.0,
        ycenter::Real=70000.0,
        rho_metal::Real=7800.0,
        V_marker=nothing,
        use_3d_volume::Bool=true,
    )

Aggregate modal accessory mineral abundances across planetesimal core, mantle, and crust regions.

Parameters
----------
- `xm, ym`: Marker coordinate arrays [m].
- `tm`: Marker type array (1=silicate/metal, 2=crust/ice, 3=sticky air).
- `tkm`: Marker temperature array [K].
- `Xfe_bulk`: Marker bulk metal volume fraction array [-].
- `Xfe_S_m, Xfe_C_m, Xfe_N_m`: Marker volatile concentration arrays in metal [ppmw].
- `marknum`: Number of markers.

Returns
-------
- Named tuple containing integrated regional masses [kg] and diagnostic meteorite classification.
"""
function compute_regional_mineral_modes(
    xm,
    ym,
    tm,
    tkm,
    Xfe_bulk,
    Xfe_S_m,
    Xfe_C_m,
    Xfe_N_m,
    marknum;
    cfg::PhaseTrackingConfig=PhaseTrackingConfig(),
    rplanet::Real=50000.0,
    xcenter::Real=70000.0,
    ycenter::Real=70000.0,
    rho_metal::Real=7800.0,
    V_marker=nothing,
    use_3d_volume::Bool=true,
)
    M_total_metal = 0.0
    M_total_troilite = 0.0
    M_total_schreibersite = 0.0
    M_total_cohenite = 0.0
    M_total_graphite = 0.0
    M_total_nitride = 0.0
    M_total_metal_matrix = 0.0
    M_total_liquid_alloy = 0.0

    M_core_metal = 0.0
    M_core_troilite = 0.0
    M_core_schreibersite = 0.0
    M_core_cohenite = 0.0
    M_core_graphite = 0.0
    M_core_nitride = 0.0
    M_core_metal_matrix = 0.0
    M_core_liquid_alloy = 0.0

    M_mantle_metal = 0.0
    M_mantle_troilite = 0.0
    M_mantle_schreibersite = 0.0
    M_mantle_cohenite = 0.0
    M_mantle_graphite = 0.0
    M_mantle_nitride = 0.0
    M_mantle_metal_matrix = 0.0
    M_mantle_liquid_alloy = 0.0

    M_crust_metal = 0.0
    M_crust_troilite = 0.0
    M_crust_schreibersite = 0.0
    M_crust_cohenite = 0.0
    M_crust_graphite = 0.0
    M_crust_nitride = 0.0
    M_crust_metal_matrix = 0.0
    M_crust_liquid_alloy = 0.0

    marknum >= 0 || throw(ArgumentError("marknum must be non-negative, got $marknum"))
    length(xm) >= marknum || throw(DimensionMismatch("length(xm) must be >= marknum"))
    length(ym) >= marknum || throw(DimensionMismatch("length(ym) must be >= marknum"))
    length(tm) >= marknum || throw(DimensionMismatch("length(tm) must be >= marknum"))
    length(tkm) >= marknum || throw(DimensionMismatch("length(tkm) must be >= marknum"))
    (cfg.r_core_norm < cfg.r_mantle_norm) ||
        throw(ArgumentError("r_core_norm must be strictly less than r_mantle_norm"))

    r_p = Float64(rplanet)
    rc_cut = r_p * clamp(cfg.r_core_norm, 0.0, 1.0)
    rm_cut = r_p * clamp(cfg.r_mantle_norm, 0.0, 1.0)
    (rho_metal > 0.0 && isfinite(rho_metal)) ||
        throw(DomainError(rho_metal, "rho_metal must be strictly positive and finite"))
    rho_m = Float64(rho_metal)

    N_planet = 0
    @inbounds for m in 1:marknum
        if tm[m] < 3
            dx = xm[m] - xcenter
            dy = ym[m] - ycenter
            if sqrt(dx^2 + dy^2) <= rplanet
                N_planet += 1
            end
        end
    end

    V_tot = use_3d_volume ? (4.0 / 3.0) * pi * r_p^3 : pi * r_p^2
    V_m = if V_marker !== nothing
        Float64(V_marker)
    elseif N_planet > 0
        V_tot / N_planet
    else
        1.0
    end

    w_P = cfg.bulk_P_ppm * 1.0e-6

    @inbounds for m in 1:marknum
        if tm[m] < 3
            dx = xm[m] - xcenter
            dy = ym[m] - ycenter
            rmark = sqrt(dx^2 + dy^2)
            if rmark <= rplanet
                fe_frac = Xfe_bulk !== nothing ? Xfe_bulk[m] : 0.0
                if fe_frac > 0.0
                    dM_fe = fe_frac * rho_m * V_m
                    w_S = Xfe_S_m !== nothing ? Xfe_S_m[m] * 1.0e-6 : 0.0
                    w_C = Xfe_C_m !== nothing ? Xfe_C_m[m] * 1.0e-6 : 0.0
                    w_N = Xfe_N_m !== nothing ? Xfe_N_m[m] * 1.0e-6 : 0.0

                    res = compute_normative_mineral_assemblage(
                        tkm[m], w_S, w_C, w_N, w_P, cfg
                    )

                    dM_tro = dM_fe * res.w_troilite
                    dM_sch = dM_fe * res.w_schreibersite
                    dM_coh = dM_fe * res.w_cohenite
                    dM_gra = dM_fe * res.w_graphite
                    dM_nit = dM_fe * res.w_nitride
                    dM_mat = dM_fe * res.w_metal_matrix
                    dM_liq = dM_fe * res.w_liquid_alloy

                    M_total_metal += dM_fe
                    M_total_troilite += dM_tro
                    M_total_schreibersite += dM_sch
                    M_total_cohenite += dM_coh
                    M_total_graphite += dM_gra
                    M_total_nitride += dM_nit
                    M_total_metal_matrix += dM_mat
                    M_total_liquid_alloy += dM_liq

                    if rmark <= rc_cut
                        M_core_metal += dM_fe
                        M_core_troilite += dM_tro
                        M_core_schreibersite += dM_sch
                        M_core_cohenite += dM_coh
                        M_core_graphite += dM_gra
                        M_core_nitride += dM_nit
                        M_core_metal_matrix += dM_mat
                        M_core_liquid_alloy += dM_liq
                    elseif rmark <= rm_cut
                        M_mantle_metal += dM_fe
                        M_mantle_troilite += dM_tro
                        M_mantle_schreibersite += dM_sch
                        M_mantle_cohenite += dM_coh
                        M_mantle_graphite += dM_gra
                        M_mantle_nitride += dM_nit
                        M_mantle_metal_matrix += dM_mat
                        M_mantle_liquid_alloy += dM_liq
                    else
                        M_crust_metal += dM_fe
                        M_crust_troilite += dM_tro
                        M_crust_schreibersite += dM_sch
                        M_crust_cohenite += dM_coh
                        M_crust_graphite += dM_gra
                        M_crust_nitride += dM_nit
                        M_crust_metal_matrix += dM_mat
                        M_crust_liquid_alloy += dM_liq
                    end
                end
            end
        end
    end

    f_molten_core = M_core_metal > 0.0 ? M_core_liquid_alloy / M_core_metal : 0.0
    f_crust_solid_acc = if M_crust_metal > 0.0
        (M_crust_troilite + M_crust_schreibersite + M_crust_cohenite) / M_crust_metal
    else
        0.0
    end

    classification =
        if f_molten_core >= 0.8 &&
            (M_core_metal / max(M_total_metal, 1.0e-12)) >= 0.4 &&
            f_crust_solid_acc <= 0.005
            :magmatic_differentiated
        elseif f_crust_solid_acc >= 0.01 && f_molten_core <= 0.6
            :IAB_winonaite_primitive
        else
            :transitional
        end

    return (;
        M_total_metal,
        M_total_troilite,
        M_total_schreibersite,
        M_total_cohenite,
        M_total_graphite,
        M_total_nitride,
        M_total_metal_matrix,
        M_total_liquid_alloy,
        M_core_metal,
        M_core_troilite,
        M_core_schreibersite,
        M_core_cohenite,
        M_core_graphite,
        M_core_nitride,
        M_core_metal_matrix,
        M_core_liquid_alloy,
        M_mantle_metal,
        M_mantle_troilite,
        M_mantle_schreibersite,
        M_mantle_cohenite,
        M_mantle_graphite,
        M_mantle_nitride,
        M_mantle_metal_matrix,
        M_mantle_liquid_alloy,
        M_crust_metal,
        M_crust_troilite,
        M_crust_schreibersite,
        M_crust_cohenite,
        M_crust_graphite,
        M_crust_nitride,
        M_crust_metal_matrix,
        M_crust_liquid_alloy,
        classification,
    )
end

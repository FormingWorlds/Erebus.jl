# =============================================================================
# Magma Ocean Multi-Component Volatile Partitioning & Degassing Engine
# =============================================================================

using DocStringExtensions

"""
    partial_pressures_to_masses(
        p_dict::Dict{Symbol,Float64},
        P_total::Real,
        col_coeff::Real,
        amu_dict::Dict{Symbol,Float64}=SPECIES_AMU,
    )

Convert species partial pressures to atmospheric column masses [kg].
...
"""
function partial_pressures_to_masses(
    p_dict::Dict{Symbol,Float64},
    P_total::Real,
    col_coeff::Real,
    amu_dict::Dict{Symbol,Float64}=SPECIES_AMU,
)
    P_tot = Float64(P_total)
    col_c = Float64(col_coeff)
    if P_tot < 0.0
        throw(DomainError(P_tot, "Total pressure P_total must be non-negative"))
    end
    if col_c < 0.0
        throw(DomainError(col_c, "Column coefficient must be non-negative"))
    end
    for (k, v) in p_dict
        if v < 0.0
            throw(DomainError(v, "Species partial pressure for $k must be non-negative"))
        end
    end
    if P_tot == 0.0
        return Dict{Symbol,Float64}(k => 0.0 for k in keys(p_dict))
    end
    numerator = sum(p_dict[k] * amu_dict[k] for k in keys(p_dict))
    mu_bar = numerator / P_tot
    if mu_bar <= 0.0
        return Dict{Symbol,Float64}(k => 0.0 for k in keys(p_dict))
    end
    return Dict{Symbol,Float64}(
        k => (p_dict[k] * amu_dict[k] / mu_bar) * col_c for k in keys(p_dict)
    )
end

const PICARD_WARNING_COUNTER = Ref{Int}(0)

"""
    get_picard_warning_count()

Return the number of times the magma ocean volatile partitioning solver fell back to Picard iteration.
"""
function get_picard_warning_count()
    return PICARD_WARNING_COUNTER[]
end

"""
    reset_picard_warning_count!()

Reset the Picard fallback warning counter to zero.
"""
function reset_picard_warning_count!()
    PICARD_WARNING_COUNTER[] = 0
    return 0
end

"""
    solve_magma_ocean_volatile_partitioning(
        M_melt::Real,
        M_tot_H::Real,
        M_tot_C::Real,
        M_tot_N::Real,
        M_tot_S::Real,
        R_planet::Real,
        g::Real,
        T_mo::Real,
        delta_IW::Real;
        kwargs...
    )

Mass-conserved volatile partitioning between liquid silicate melt and atmosphere.

Solves the coupled non-linear system of Dalton partial pressures, gas-phase chemical equilibrium,
and silicate melt solubility across ten gas species:
    H2, H2O, CO, CO2, CH4, N2, NH3, H2S, S2, SO2
enforcing strict elemental mass conservation for H, C, N, and S between magma ocean melt (M_melt)
and the overlying atmosphere:
    M_tot,E = M_melt * w_E + (4π R_p^2 / g) * ∑_i p_i * (m_E * ν_{i, E} / m_i)

# Arguments
- `M_melt`: Total silicate melt mass [kg].
- `M_tot_H`: Total elemental hydrogen inventory in magma ocean system [kg].
- `M_tot_C`: Total elemental carbon inventory in magma ocean system [kg].
- `M_tot_N`: Total elemental nitrogen inventory in magma ocean system [kg].
- `M_tot_S`: Total elemental sulfur inventory in magma ocean system [kg].
- `R_planet`: Planetary surface radius [m].
- `g`: Surface gravitational acceleration [m/s^2].
- `T_mo`: Magma ocean reference melt temperature [K].
- `delta_IW`: Oxygen fugacity offset relative to iron-wüstite buffer [log10 units].

# Keyword Arguments
- `water_As`: Burnham/Dixon water solubility coefficient [wt% / MPa^0.5] (default: 0.40).
- `water_law`: Water solubility formulation (default: `:burnham_dixon`).
- `carbon_active`: Whether carbon dissolves in melt (default: true).
- `sulfur_active`: Whether sulfur dissolves in melt (default: true).
- `co_law`: CO solubility formulation (default: `:armstrong2015`).
- `ch4_law`: CH4 solubility formulation (default: `:ardia2013`).
- `co2_law`: CO2 solubility formulation (default: `:dixon1995`).
- `sulfide_law`: Sulfide solubility formulation (default: `:boulliung2023`).
- `nitrogen_henry`: Physical N2 Henry coefficient [ppm / bar] (default: 0.40).
- `nitrogen_nitride`: Chemical nitride capacity [wt% / bar^0.5] (default: 1.0e-3).
- `graphite_saturation`: Whether to cap carbon fugacities at graphite saturation (default: true).
- `max_newton_iter`: Maximum Newton iterations before Picard fallback (default: 100).
- `max_picard_iter`: Maximum Picard iterations before throwing ConvergenceError (default: 500).
- `M_tot_O`: Initial/reference elemental oxygen inventory [kg] (default: 0.0).

# Returns
- `NamedTuple`:
  - `P_surf`: Total equilibrium surface atmospheric pressure [Pa].
  - `p_i`: Dictionary of species partial pressures [Pa].
  - `M_atm_i`: Dictionary of species atmospheric masses [kg].
  - `M_atm_tot`: Total atmospheric mass [kg].
  - `M_atm_H`: Elemental hydrogen mass in atmosphere [kg].
  - `M_atm_C`: Elemental carbon mass in atmosphere [kg].
  - `M_atm_N`: Elemental nitrogen mass in atmosphere [kg].
  - `M_atm_S`: Elemental sulfur mass in atmosphere [kg].
  - `M_melt_H`: Elemental hydrogen mass retained in silicate melt [kg].
  - `M_melt_C`: Elemental carbon mass retained in silicate melt [kg].
  - `M_melt_N`: Elemental nitrogen mass retained in silicate melt [kg].
  - `M_melt_S`: Elemental sulfur mass retained in silicate melt [kg].
  - `M_melt_tot`: Total volatile mass dissolved in melt [kg].
  - `w_diss_H2O_wtpct`: Dissolved water concentration in melt [wt%].
  - `C_diss_C_ppm`: Dissolved carbon concentration in melt [ppmw].
  - `C_diss_N_ppm`: Dissolved nitrogen concentration in melt [ppmw].
  - `C_diss_S_ppm`: Dissolved sulfur concentration in melt [ppmw].
  - `M_graphite`: Precipitated solid graphite mass [kg].
  - `dO_buffer`: Net oxygen exchanged with rock/melt buffer [kg].
  - `M_atm_O`: Elemental oxygen mass in atmosphere [kg].
  - `warning_counter`: Picard fallback count in telemetry.
"""
function solve_magma_ocean_volatile_partitioning(
    M_melt::Real,
    M_tot_H::Real,
    M_tot_C::Real,
    M_tot_N::Real,
    M_tot_S::Real,
    R_planet::Real,
    g::Real,
    T_mo::Real,
    delta_IW::Real;
    water_As::Real=0.40,
    water_law::Symbol=:burnham_dixon,
    carbon_active::Bool=true,
    sulfur_active::Bool=true,
    co_law::Symbol=:armstrong2015,
    ch4_law::Symbol=:ardia2013,
    co2_law::Symbol=:dixon1995,
    sulfide_law::Symbol=:boulliung2023,
    nitrogen_henry::Real=0.40,
    nitrogen_nitride::Real=1.0e-3,
    graphite_saturation::Bool=true,
    max_newton_iter::Int=100,
    max_picard_iter::Int=500,
    M_tot_O::Real=0.0,
)::@NamedTuple{
    P_surf::Float64,
    p_i::Dict{Symbol,Float64},
    M_atm_i::Dict{Symbol,Float64},
    M_atm_tot::Float64,
    M_atm_H::Float64,
    M_atm_C::Float64,
    M_atm_N::Float64,
    M_atm_S::Float64,
    M_melt_H::Float64,
    M_melt_C::Float64,
    M_melt_N::Float64,
    M_melt_S::Float64,
    M_melt_tot::Float64,
    w_diss_H2O_wtpct::Float64,
    C_diss_C_ppm::Float64,
    C_diss_N_ppm::Float64,
    C_diss_S_ppm::Float64,
    M_graphite::Float64,
    is_graphite_sat::Bool,
    dO_buffer::Float64,
    M_atm_O::Float64,
    warning_counter::Int,
}
    M_m = Float64(M_melt)
    mH = Float64(M_tot_H)
    mC = Float64(M_tot_C)
    mN = Float64(M_tot_N)
    mS = Float64(M_tot_S)
    Rp = Float64(R_planet)
    grav = Float64(g)
    T = Float64(T_mo)
    dIW = Float64(delta_IW)

    if M_m < 0.0 || !isfinite(M_m)
        throw(DomainError(M_m, "Melt mass must be non-negative and finite"))
    end
    if mH < 0.0 || !isfinite(mH)
        throw(DomainError(mH, "Hydrogen mass must be non-negative and finite"))
    end
    if mC < 0.0 || !isfinite(mC)
        throw(DomainError(mC, "Carbon mass must be non-negative and finite"))
    end
    if mN < 0.0 || !isfinite(mN)
        throw(DomainError(mN, "Nitrogen mass must be non-negative and finite"))
    end
    if mS < 0.0 || !isfinite(mS)
        throw(DomainError(mS, "Sulfur mass must be non-negative and finite"))
    end
    if Rp <= 0.0 || !isfinite(Rp)
        throw(DomainError(Rp, "Planetary radius must be > 0 and finite"))
    end
    if grav <= 0.0 || !isfinite(grav)
        throw(DomainError(grav, "Surface gravity must be > 0 and finite"))
    end
    if T <= 0.0 || !isfinite(T)
        throw(DomainError(T, "Temperature must be > 0 and finite"))
    end
    if !isfinite(dIW) || abs(dIW) > 50.0
        throw(DomainError(dIW, "delta_IW must be finite and within [-50, 50]"))
    end

    empty_p_i = Dict{Symbol,Float64}(
        :H2 => 0.0,
        :H2O => 0.0,
        :CO => 0.0,
        :CO2 => 0.0,
        :CH4 => 0.0,
        :N2 => 0.0,
        :NH3 => 0.0,
        :H2S => 0.0,
        :S2 => 0.0,
        :SO2 => 0.0,
    )

    total_volatile_mass = mH + mC + mN + mS
    if total_volatile_mass == 0.0
        return (
            P_surf=0.0,
            p_i=empty_p_i,
            M_atm_i=copy(empty_p_i),
            M_atm_tot=0.0,
            M_atm_H=0.0,
            M_atm_C=0.0,
            M_atm_N=0.0,
            M_atm_S=0.0,
            M_melt_H=0.0,
            M_melt_C=0.0,
            M_melt_N=0.0,
            M_melt_S=0.0,
            M_melt_tot=0.0,
            w_diss_H2O_wtpct=0.0,
            C_diss_C_ppm=0.0,
            C_diss_N_ppm=0.0,
            C_diss_S_ppm=0.0,
            M_graphite=0.0,
            is_graphite_sat=false,
            dO_buffer=0.0,
            M_atm_O=0.0,
            warning_counter=PICARD_WARNING_COUNTER[],
        )
    end

    col_coeff = (4.0 * π * (Rp^2)) / grav
    log10_fO2 = compute_iron_wustite_fO2(T; delta_IW=dIW)

    logK_H2O = 12700.0 / T - 2.80
    r_H = 10.0^clamp(logK_H2O + 0.5 * log10_fO2, -100.0, 100.0)

    logK_CO2 = 14800.0 / T - 4.58
    r_CO2 = 10.0^clamp(logK_CO2 + 0.5 * log10_fO2, -100.0, 100.0)

    logK_SO2 = 18800.0 / T - 3.80
    r_SO2 = 10.0^clamp(logK_SO2 + log10_fO2, -100.0, 100.0)

    gr = compute_graphite_saturation_fugacity(T, log10_fO2)
    f_CO_max_Pa = gr.f_CO_max_bar * 1.0e5
    f_CO2_max_Pa = gr.f_CO2_max_bar * 1.0e5

    tol_H = 1.0e-10 * mH + 1.0e-12 * total_volatile_mass
    tol_C = 1.0e-10 * mC + 1.0e-12 * total_volatile_mass
    tol_N = 1.0e-10 * mN + 1.0e-12 * total_volatile_mass
    tol_S = 1.0e-10 * mS + 1.0e-12 * total_volatile_mass
    tols = [tol_H, tol_C, tol_N, tol_S]

    function eval_coupled_state(u_vec::AbstractVector{Float64})
        p_H2O = mH > 0.0 ? exp(u_vec[1]) : 0.0
        p_CO2_raw = mC > 0.0 ? exp(u_vec[2]) : 0.0
        p_N2 = mN > 0.0 ? exp(u_vec[3]) : 0.0
        p_SO2 = mS > 0.0 ? exp(u_vec[4]) : 0.0

        p_H2 = p_H2O / r_H

        p_CO2 = p_CO2_raw
        p_CO = p_CO2 / r_CO2
        is_graphite_sat = false
        if graphite_saturation && mC > 0.0
            if p_CO >= f_CO_max_Pa || p_CO2 >= f_CO2_max_Pa
                is_graphite_sat = true
                p_CO = f_CO_max_Pa
                p_CO2 = f_CO2_max_Pa
            end
        end

        p_H2_bar = p_H2 * 1.0e-5
        p_CH4 = 0.0
        if p_H2 > 0.0 && p_CO > 0.0
            l_pH2 = log10(max(p_H2_bar, 1.0e-30))
            r_CH4 =
                10.0^clamp(
                    11500.0 / T - 12.0 + 2.0 * l_pH2 - log10(max(r_H, 1.0e-30)),
                    -100.0,
                    100.0,
                )
            p_CH4 = r_CH4 * p_CO
        end

        p_N2_bar = p_N2 * 1.0e-5
        p_NH3 = 0.0
        if p_H2 > 0.0 && p_N2 > 0.0
            l_pH2 = log10(max(p_H2_bar, 1.0e-30))
            r_NH3 = 10.0^clamp(2800.0 / T - 5.80 + 1.5 * l_pH2, -100.0, 100.0)
            p_NH3 = r_NH3 * sqrt(max(p_N2_bar, 0.0)) * 1.0e5
        end

        p_SO2_bar = p_SO2 * 1.0e-5
        p_S2 = 0.0
        p_H2S = 0.0
        if p_SO2 > 0.0
            sqrt_pS2_bar = p_SO2_bar / max(r_SO2, 1.0e-30)
            p_S2 = (sqrt_pS2_bar^2) * 1.0e5
            if p_H2 > 0.0
                l_pH2 = log10(max(p_H2_bar, 1.0e-30))
                r_H2S = 10.0^clamp(4800.0 / T - 2.50 + l_pH2, -100.0, 100.0)
                p_H2S = r_H2S * sqrt_pS2_bar * 1.0e5
            end
        end

        p_dict = Dict{Symbol,Float64}(
            :H2 => p_H2,
            :H2O => p_H2O,
            :CO => p_CO,
            :CO2 => p_CO2,
            :CH4 => p_CH4,
            :N2 => p_N2,
            :NH3 => p_NH3,
            :H2S => p_H2S,
            :S2 => p_S2,
            :SO2 => p_SO2,
        )

        P_surf = sum(values(p_dict))
        if P_surf <= 0.0
            return (;
                P_surf=0.0,
                p_dict=p_dict,
                m_atm_dict=copy(empty_p_i),
                M_atm_tot=0.0,
                M_atm_H=0.0,
                M_atm_C=0.0,
                M_atm_N=0.0,
                M_atm_S=0.0,
                M_atm_O=0.0,
                M_melt_H=0.0,
                M_melt_C=0.0,
                M_melt_N=0.0,
                M_melt_S=0.0,
                w_diss_H2O_wtpct=0.0,
                C_diss_C_ppm=0.0,
                C_diss_N_ppm=0.0,
                C_diss_S_ppm=0.0,
                M_graphite=0.0,
                is_graphite_sat=false,
                R_H=(-mH),
                R_C=(-mC),
                R_N=(-mN),
                R_S=(-mS),
            )
        end

        mu_bar = sum(p_dict[k] * SPECIES_AMU[k] for k in keys(p_dict)) / P_surf
        m_atm_dict = Dict{Symbol,Float64}(
            k => col_coeff * p_dict[k] * (SPECIES_AMU[k] / mu_bar) for k in keys(p_dict)
        )
        M_atm_tot = col_coeff * P_surf

        M_atm_H = (
            m_atm_dict[:H2] +
            m_atm_dict[:H2O] * (2.0 * SPECIES_AMU[:H] / SPECIES_AMU[:H2O]) +
            m_atm_dict[:CH4] * (4.0 * SPECIES_AMU[:H] / SPECIES_AMU[:CH4]) +
            m_atm_dict[:NH3] * (3.0 * SPECIES_AMU[:H] / SPECIES_AMU[:NH3]) +
            m_atm_dict[:H2S] * (2.0 * SPECIES_AMU[:H] / SPECIES_AMU[:H2S])
        )
        M_atm_C = (
            m_atm_dict[:CO] * (SPECIES_AMU[:C] / SPECIES_AMU[:CO]) +
            m_atm_dict[:CO2] * (SPECIES_AMU[:C] / SPECIES_AMU[:CO2]) +
            m_atm_dict[:CH4] * (SPECIES_AMU[:C] / SPECIES_AMU[:CH4])
        )
        M_atm_N = (
            m_atm_dict[:N2] + m_atm_dict[:NH3] * (SPECIES_AMU[:N] / SPECIES_AMU[:NH3])
        )
        M_atm_S = (
            m_atm_dict[:S2] +
            m_atm_dict[:SO2] * (SPECIES_AMU[:S] / SPECIES_AMU[:SO2]) +
            m_atm_dict[:H2S] * (SPECIES_AMU[:S] / SPECIES_AMU[:H2S])
        )
        M_atm_O = (
            m_atm_dict[:H2O] * (SPECIES_AMU[:O] / SPECIES_AMU[:H2O]) +
            m_atm_dict[:CO] * (SPECIES_AMU[:O] / SPECIES_AMU[:CO]) +
            m_atm_dict[:CO2] * (2.0 * SPECIES_AMU[:O] / SPECIES_AMU[:CO2]) +
            m_atm_dict[:SO2] * (2.0 * SPECIES_AMU[:O] / SPECIES_AMU[:SO2])
        )

        w_diss_H2O_wtpct = compute_water_solubility_melt(
            p_dict[:H2O]; As=water_As, law=water_law
        )
        w_H2O = w_diss_H2O_wtpct * 0.01
        w_H2 = compute_h2_solubility_melt(p_dict[:H2]) * 1.0e-6
        w_diss_H = w_H2O * (2.0 * SPECIES_AMU[:H] / SPECIES_AMU[:H2O]) + w_H2
        M_melt_H = M_m * w_diss_H

        C_diss_C_ppm = 0.0
        if carbon_active && mC > 0.0
            co_ppm = compute_co_solubility_melt(p_dict[:CO], P_surf; law=co_law)
            ch4_ppm = compute_ch4_solubility_melt(p_dict[:CH4], P_surf; law=ch4_law)
            co2_ppm = compute_co2_solubility_melt(p_dict[:CO2], T; law=co2_law)
            C_diss_C_ppm = co_ppm + ch4_ppm + co2_ppm
        end
        w_diss_C = C_diss_C_ppm * 1.0e-6
        M_melt_C = M_m * w_diss_C

        C_diss_N_ppm = 0.0
        if mN > 0.0
            S_N_res = compute_nitrogen_solubility_melt(
                p_dict[:N2], dIW; Kh=nitrogen_henry, C_nitride=nitrogen_nitride
            )
            C_diss_N_ppm = S_N_res.total_ppm
        end
        w_diss_N = C_diss_N_ppm * 1.0e-6
        M_melt_N = M_m * w_diss_N

        C_diss_S_ppm = 0.0
        if sulfur_active && mS > 0.0
            C_diss_S_ppm = compute_sulfur_solubility_melt(
                p_dict[:S2], T, dIW; law=sulfide_law
            )
        end
        w_diss_S = C_diss_S_ppm * 1.0e-6
        M_melt_S = M_m * w_diss_S

        M_graphite = 0.0
        if graphite_saturation && mC > 0.0
            if is_graphite_sat || (
                (M_melt_C + M_atm_C) < mC &&
                (p_CO >= f_CO_max_Pa * 0.999999 || p_CO2 >= f_CO2_max_Pa * 0.999999)
            )
                is_graphite_sat = true
                M_graphite = max(0.0, mC - (M_melt_C + M_atm_C))
            end
        end

        M_calc_H = M_melt_H + M_atm_H
        M_calc_C = M_melt_C + M_atm_C + M_graphite
        M_calc_N = M_melt_N + M_atm_N
        M_calc_S = M_melt_S + M_atm_S

        R_H = mH > 0.0 ? M_calc_H - mH : 0.0
        R_C = mC > 0.0 ? M_calc_C - mC : 0.0
        R_N = mN > 0.0 ? M_calc_N - mN : 0.0
        R_S = mS > 0.0 ? M_calc_S - mS : 0.0

        return (;
            P_surf=P_surf,
            p_dict=p_dict,
            m_atm_dict=m_atm_dict,
            M_atm_tot=M_atm_tot,
            M_atm_H=M_atm_H,
            M_atm_C=M_atm_C,
            M_atm_N=M_atm_N,
            M_atm_S=M_atm_S,
            M_atm_O=M_atm_O,
            M_melt_H=M_melt_H,
            M_melt_C=M_melt_C,
            M_melt_N=M_melt_N,
            M_melt_S=M_melt_S,
            w_diss_H2O_wtpct=w_diss_H2O_wtpct,
            C_diss_C_ppm=C_diss_C_ppm,
            C_diss_N_ppm=C_diss_N_ppm,
            C_diss_S_ppm=C_diss_S_ppm,
            M_graphite=M_graphite,
            is_graphite_sat=is_graphite_sat,
            R_H=R_H,
            R_C=R_C,
            R_N=R_N,
            R_S=R_S,
        )
    end

    P_max = total_volatile_mass / col_coeff
    p_H2O_0 = if mH > 0.0
        max(1.0e-5, (mH / total_volatile_mass) * P_max * (r_H / (1.0 + r_H)))
    else
        1.0e-20
    end
    p_CO2_0 = if mC > 0.0
        max(1.0e-5, (mC / total_volatile_mass) * P_max * (r_CO2 / (1.0 + r_CO2)))
    else
        1.0e-20
    end
    if graphite_saturation && p_CO2_0 > f_CO2_max_Pa
        p_CO2_0 = f_CO2_max_Pa
    end
    p_N2_0 = mN > 0.0 ? max(1.0e-5, (mN / total_volatile_mass) * P_max) : 1.0e-20
    p_SO2_0 = if mS > 0.0
        max(1.0e-5, (mS / total_volatile_mass) * P_max * (r_SO2 / (1.0 + r_SO2)))
    else
        1.0e-20
    end

    u = [log(p_H2O_0), log(p_CO2_0), log(p_N2_0), log(p_SO2_0)]
    active_indices = Int[]
    mH > 0.0 && push!(active_indices, 1)
    mC > 0.0 && push!(active_indices, 2)
    mN > 0.0 && push!(active_indices, 3)
    mS > 0.0 && push!(active_indices, 4)

    state = eval_coupled_state(u)
    converged = false

    res_vec = [state.R_H, state.R_C, state.R_N, state.R_S]
    norm_0 = maximum(abs(res_vec[i]) / tols[i] for i in 1:4)
    if norm_0 <= 1.0
        converged = true
    end

    if !converged
        for iter in 1:max_newton_iter
            curr_active = Int[]
            for idx in active_indices
                if idx == 2 && state.is_graphite_sat
                    continue
                end
                push!(curr_active, idx)
            end

            n_act = length(curr_active)
            if n_act == 0
                converged = true
                break
            end

            F_curr = [res_vec[idx] for idx in curr_active]
            J = zeros(Float64, n_act, n_act)
            h = 1.0e-6

            for col in 1:n_act
                idx = curr_active[col]
                u_pert = copy(u)
                u_pert[idx] += h
                state_pert = eval_coupled_state(u_pert)
                res_pert = [state_pert.R_H, state_pert.R_C, state_pert.R_N, state_pert.R_S]
                for row in 1:n_act
                    row_idx = curr_active[row]
                    J[row, col] = (res_pert[row_idx] - F_curr[row]) / h
                end
            end

            delta_u_act = try
                J \ (-F_curr)
            catch
                break
            end

            if any(!isfinite, delta_u_act)
                break
            end

            max_step = maximum(abs.(delta_u_act))
            if max_step > 4.0
                delta_u_act .*= (4.0 / max_step)
            end

            alpha = 1.0
            step_halvings = 0
            u_trial = copy(u)
            for (col, idx) in enumerate(curr_active)
                u_trial[idx] += alpha * delta_u_act[col]
            end
            state_trial = eval_coupled_state(u_trial)
            res_trial = [state_trial.R_H, state_trial.R_C, state_trial.R_N, state_trial.R_S]
            norm_trial = maximum(abs(res_trial[i]) / tols[i] for i in 1:4)

            while norm_trial >= norm_0 && step_halvings < 30
                alpha *= 0.5
                step_halvings += 1
                u_trial = copy(u)
                for (col, idx) in enumerate(curr_active)
                    u_trial[idx] += alpha * delta_u_act[col]
                end
                state_trial = eval_coupled_state(u_trial)
                res_trial = [
                    state_trial.R_H, state_trial.R_C, state_trial.R_N, state_trial.R_S
                ]
                norm_trial = maximum(abs(res_trial[i]) / tols[i] for i in 1:4)
            end

            if norm_trial < norm_0
                u = u_trial
                state = state_trial
                res_vec = res_trial
                norm_0 = norm_trial
                if norm_0 <= 1.0
                    converged = true
                    break
                end
            else
                break
            end
        end
    end

    if !converged
        PICARD_WARNING_COUNTER[] += 1
        for picard_iter in 1:max_picard_iter
            res_vec = [state.R_H, state.R_C, state.R_N, state.R_S]
            norm_0 = maximum(abs(res_vec[i]) / tols[i] for i in 1:4)
            if norm_0 <= 1.0
                converged = true
                break
            end

            curr_active = Int[]
            for idx in active_indices
                if idx == 2 && state.is_graphite_sat
                    continue
                end
                push!(curr_active, idx)
            end

            m_tot_targets = [mH, mC, mN, mS]
            m_calc_current = [
                state.M_melt_H + state.M_atm_H,
                state.M_melt_C + state.M_atm_C + state.M_graphite,
                state.M_melt_N + state.M_atm_N,
                state.M_melt_S + state.M_atm_S,
            ]

            for idx in curr_active
                target = m_tot_targets[idx]
                curr = m_calc_current[idx]
                if curr > 0.0 && target > 0.0
                    ratio = target / curr
                    damping = 0.5
                    u[idx] += damping * log(clamp(ratio, 0.05, 20.0))
                end
            end
            state = eval_coupled_state(u)
        end
        res_vec = [state.R_H, state.R_C, state.R_N, state.R_S]
        norm_0 = maximum(abs(res_vec[i]) / tols[i] for i in 1:4)
        if norm_0 <= 1.0
            converged = true
        end
    end

    if !converged
        throw(
            ConvergenceError(
                "Magma ocean volatile partitioning failed to converge after Newton and Picard iterations";
                residuals=Dict(
                    :H => state.R_H, :C => state.R_C, :N => state.R_N, :S => state.R_S
                ),
            ),
        )
    end

    M_atm_O = state.M_atm_O
    dO_buffer = M_atm_O - Float64(M_tot_O)

    return (
        P_surf=state.P_surf,
        p_i=state.p_dict,
        M_atm_i=state.m_atm_dict,
        M_atm_tot=state.M_atm_tot,
        M_atm_H=state.M_atm_H,
        M_atm_C=state.M_atm_C,
        M_atm_N=state.M_atm_N,
        M_atm_S=state.M_atm_S,
        M_melt_H=state.M_melt_H,
        M_melt_C=state.M_melt_C,
        M_melt_N=state.M_melt_N,
        M_melt_S=state.M_melt_S,
        M_melt_tot=state.M_melt_H + state.M_melt_C + state.M_melt_N + state.M_melt_S,
        w_diss_H2O_wtpct=state.w_diss_H2O_wtpct,
        C_diss_C_ppm=state.C_diss_C_ppm,
        C_diss_N_ppm=state.C_diss_N_ppm,
        C_diss_S_ppm=state.C_diss_S_ppm,
        M_graphite=state.M_graphite,
        is_graphite_sat=state.is_graphite_sat,
        dO_buffer=dO_buffer,
        M_atm_O=M_atm_O,
        warning_counter=PICARD_WARNING_COUNTER[],
    )
end

"""
    degas_magma_ocean_markers!(
        xm, ym, tm, tkm, Fm, Fm_old, XH2Om, XCm, XNm, XSm,
        marknum, dt, P_surf, R_planet, cfg;
        kwargs...
    )

Evaluate volatile exsolution and dynamic degassing rates from Lagrangian melt markers.

Identifies molten and ascending markers in the magma ocean or near-surface zone,
extracts supersaturated dissolved volatiles above equilibrium solubility, and generates
species-resolved degassing rates [kg/s].

# Arguments
- `xm`: Marker x-coordinates [m].
- `ym`: Marker y-coordinates [m].
- `tm`: Marker material types (tm < 3 for rock).
- `tkm`: Marker temperatures [K].
- `Fm`: Current marker silicate melt fractions.
- `Fm_old`: Previous step marker silicate melt fractions (retained for signature compatibility).
- `XH2Om`: Marker dissolved water mass fractions.
- `XCm`: Marker dissolved carbon mass fractions.
- `XNm`: Marker dissolved nitrogen mass fractions.
- `XSm`: Marker dissolved sulfur mass fractions.
- `marknum`: Number of active markers.
- `dt`: Simulation timestep [s].
- `P_surf`: Current atmospheric surface pressure [Pa].
- `R_planet`: Planetary surface radius [m].
- `cfg`: Magma ocean degassing configuration (`MagmaOceanDegassingConfig`).
- `T_melt_ref`: Reference melt temperature for saturation evaluation [K] (default: 1500.0).

# Keyword Arguments
- `xcenter`: Planetary center horizontal coordinate [m] (default: 0.0).
- `ycenter`: Planetary center vertical coordinate [m] (default: 0.0).
- `w3d_m`: Optional precomputed out-of-plane spherical integration lengths [m].
- `rho_solid`: Reference solid rock density [kg/m^3] (default: 3000.0).
- `marker_volume`: Marker 2D volume / cross-sectional area [m^2].
- `delta_IW`: Oxygen fugacity offset relative to IW (default: 0.0).
- `retention_cfg`: Volatile retention configuration (optional).
- `step`: Simulation step index for TransferRecord provenance (default: 0).

# Returns
- `NamedTuple`: `(; rates, dM_2D, dM_3D, records)` containing species rates [kg/s], 2D planar masses [kg/m], 3D spherical masses [kg], and transfer audit records.
"""
function degas_magma_ocean_markers!(
    xm::AbstractVector{Float64},
    ym::AbstractVector{Float64},
    tm::AbstractVector{Int},
    tkm::AbstractVector{Float64},
    Fm::AbstractVector{Float64},
    Fm_old::AbstractVector{Float64},
    XH2Om::AbstractVector{Float64},
    XCm::AbstractVector{Float64},
    XNm::AbstractVector{Float64},
    XSm::AbstractVector{Float64},
    marknum::Int,
    dt::Real,
    P_surf::Real,
    R_planet::Real,
    cfg::MagmaOceanDegassingConfig,
    T_melt_ref::Real=1500.0;
    xcenter::Real=0.0,
    ycenter::Real=0.0,
    w3d_m::Union{Nothing,AbstractVector{Float64}}=nothing,
    rho_solid::Real=3000.0,
    marker_volume::Real=1.0,
    delta_IW::Real=0.0,
    retention_cfg=nothing,
    step::Int=0,
)
    empty_rates = Dict{Symbol,Float64}(
        :H2 => 0.0,
        :H2O => 0.0,
        :CO => 0.0,
        :CO2 => 0.0,
        :CH4 => 0.0,
        :N2 => 0.0,
        :NH3 => 0.0,
        :H2S => 0.0,
        :S2 => 0.0,
        :SO2 => 0.0,
    )
    empty_dM = Dict{Symbol,Float64}(:H => 0.0, :C => 0.0, :N => 0.0, :S => 0.0, :H2O => 0.0)
    empty_res = (;
        rates=empty_rates, dM_2D=empty_dM, dM_3D=copy(empty_dM), records=TransferRecord[]
    )

    if !cfg.active || dt <= 0.0 || marknum == 0
        return empty_res
    end

    (w3d_m === nothing || length(w3d_m) >= marknum) || throw(
        DimensionMismatch(
            "length(w3d_m) must be >= marknum (got $(length(w3d_m)), expected $marknum)"
        ),
    )

    Rp = Float64(R_planet)
    r_degas_sq = (cfg.degas_depth_fraction * Rp)^2
    Rp_sq = Rp^2
    F_thresh = cfg.F_melt_threshold
    eff = cfg.efficiency
    psurf_val = max(1.0, Float64(P_surf))
    rho_s = Float64(rho_solid)
    v_m = Float64(marker_volume)
    m_marker = rho_s * v_m
    delta_IW_eff = cfg.redox_coupled ? Float64(delta_IW) : 0.0
    T_ref = Float64(T_melt_ref)
    spec_surf = solve_chnos_speciation(psurf_val, T_ref, delta_IW_eff)

    # Evaluate equilibrium solubilities once at surface ambient partial pressures and reference melt temperature
    S_H2O_wtpct = compute_water_solubility_melt(spec_surf.p_H2O_Pa; As=cfg.water_As)
    w_H2O_sat = S_H2O_wtpct * 0.01

    S_N_res = compute_nitrogen_solubility_melt(spec_surf.p_N2_Pa, delta_IW_eff)
    w_N_sat = S_N_res.total_ppm * 1.0e-6

    S_C_res = compute_carbon_solubility_melt(psurf_val, T_ref, delta_IW_eff)
    w_C_sat = S_C_res.total_ppm * 1.0e-6

    S_S_ppm = compute_sulfur_solubility_melt(spec_surf.p_S2_Pa, T_ref, delta_IW_eff)
    w_S_sat = S_S_ppm * 1.0e-6

    tot_ex_H2O_2D = 0.0
    tot_ex_H2O_3D = 0.0
    tot_ex_H_2D = 0.0
    tot_ex_H_3D = 0.0
    tot_ex_C_2D = 0.0
    tot_ex_C_3D = 0.0
    tot_ex_N_2D = 0.0
    tot_ex_N_3D = 0.0
    tot_ex_S_2D = 0.0
    tot_ex_S_3D = 0.0
    records = TransferRecord[]

    xc_val = Float64(xcenter)
    yc_val = Float64(ycenter)
    h_conv = 2.01588 / 18.01528

    @inbounds for m in 1:marknum
        if tm[m] >= 3
            continue
        end

        dx = xm[m] - xc_val
        dy = ym[m] - yc_val
        r_sq = dx * dx + dy * dy
        if r_sq > Rp_sq
            continue
        end

        F_curr = Fm[m]
        if F_curr <= 0.0
            continue
        end

        # Check degassing activation
        is_degassing_zone = (r_sq >= r_degas_sq) && (F_curr >= F_thresh || F_curr > 0.01)
        if !is_degassing_zone
            continue
        end

        # Retention floor bounds (if enabled)
        ret_H2O = 0.0
        ret_C = 0.0
        ret_N = 0.0
        ret_S = 0.0
        if retention_cfg !== nothing && retention_cfg.active
            ret_H2O =
                compute_h2o_retention_floor(
                    tkm[m], retention_cfg; F_melt=F_curr, P_val=psurf_val
                ) * 1.0e-6
            ret_N =
                compute_nitrogen_retention_floor(
                    tkm[m], retention_cfg; F_melt=F_curr, P_val=psurf_val
                ) * 1.0e-6
            ret_C =
                compute_carbon_retention_floor(
                    tkm[m], retention_cfg; F_melt=F_curr, P_val=psurf_val
                ) * 1.0e-6
            ret_S =
                compute_sulfur_retention_floor(
                    tkm[m], retention_cfg; F_melt=F_curr, P_val=psurf_val
                ) * 1.0e-6
        end

        # Marker volatile concentrations converted to dimensionless mass fractions
        w_H2O_m = XH2Om[m] * 0.01
        w_C_m = XCm[m] * 1.0e-6
        w_N_m = XNm[m] * 1.0e-6
        w_S_m = XSm[m] * 1.0e-6

        # Supersaturated volatile extraction (evaluating supersaturation in the melt volume, bounded by retention floors)
        ex_H2O = min(
            max(0.0, w_H2O_m - ret_H2O),
            max(0.0, w_H2O_m / F_curr - w_H2O_sat) * F_curr * eff,
        )
        ex_C = min(
            max(0.0, w_C_m - ret_C), max(0.0, w_C_m / F_curr - w_C_sat) * F_curr * eff
        )
        ex_N = min(
            max(0.0, w_N_m - ret_N), max(0.0, w_N_m / F_curr - w_N_sat) * F_curr * eff
        )
        ex_S = min(
            max(0.0, w_S_m - ret_S), max(0.0, w_S_m / F_curr - w_S_sat) * F_curr * eff
        )

        w3d = w3d_m !== nothing ? w3d_m[m] : (2.0 * sqrt(r_sq))

        if ex_H2O > 0.0
            XH2Om[m] = max(0.0, (w_H2O_m - ex_H2O) * 100.0)
            dM2_h2o = ex_H2O * m_marker
            dM3_h2o = dM2_h2o * w3d
            tot_ex_H2O_2D += dM2_h2o
            tot_ex_H2O_3D += dM3_h2o
            dM2_h = dM2_h2o * h_conv
            dM3_h = dM3_h2o * h_conv
            tot_ex_H_2D += dM2_h
            tot_ex_H_3D += dM3_h
            push!(
                records, TransferRecord(step, :degassing, :H, m, xm[m], ym[m], dM2_h, dM3_h)
            )
        end
        if ex_C > 0.0
            XCm[m] = max(0.0, (w_C_m - ex_C) * 1.0e6)
            dM2_c = ex_C * m_marker
            dM3_c = dM2_c * w3d
            tot_ex_C_2D += dM2_c
            tot_ex_C_3D += dM3_c
            push!(
                records, TransferRecord(step, :degassing, :C, m, xm[m], ym[m], dM2_c, dM3_c)
            )
        end
        if ex_N > 0.0
            XNm[m] = max(0.0, (w_N_m - ex_N) * 1.0e6)
            dM2_n = ex_N * m_marker
            dM3_n = dM2_n * w3d
            tot_ex_N_2D += dM2_n
            tot_ex_N_3D += dM3_n
            push!(
                records, TransferRecord(step, :degassing, :N, m, xm[m], ym[m], dM2_n, dM3_n)
            )
        end
        if ex_S > 0.0
            XSm[m] = max(0.0, (w_S_m - ex_S) * 1.0e6)
            dM2_s = ex_S * m_marker
            dM3_s = dM2_s * w3d
            tot_ex_S_2D += dM2_s
            tot_ex_S_3D += dM3_s
            push!(
                records, TransferRecord(step, :degassing, :S, m, xm[m], ym[m], dM2_s, dM3_s)
            )
        end
    end

    dM_2D = Dict{Symbol,Float64}(
        :H => tot_ex_H_2D,
        :C => tot_ex_C_2D,
        :N => tot_ex_N_2D,
        :S => tot_ex_S_2D,
        :H2O => tot_ex_H2O_2D,
    )
    dM_3D = Dict{Symbol,Float64}(
        :H => tot_ex_H_3D,
        :C => tot_ex_C_3D,
        :N => tot_ex_N_3D,
        :S => tot_ex_S_3D,
        :H2O => tot_ex_H2O_3D,
    )

    if (tot_ex_H2O_3D + tot_ex_C_3D + tot_ex_N_3D + tot_ex_S_3D) == 0.0
        return (; rates=empty_rates, dM_2D=dM_2D, dM_3D=dM_3D, records=records)
    end

    # Thermodynamic gas speciation of degassed volatile mixture
    spec_dict = speciate_vented_volatiles(
        tot_ex_H2O_3D,
        tot_ex_C_3D,
        tot_ex_N_3D,
        tot_ex_S_3D,
        psurf_val,
        T_ref,
        delta_IW_eff;
        graphite_saturation=true,
    )

    rates = Dict{Symbol,Float64}()
    dt_val = Float64(dt)
    for (sp, mass) in spec_dict
        rates[sp] = mass / dt_val
    end

    return (; rates=rates, dM_2D=dM_2D, dM_3D=dM_3D, records=records)
end

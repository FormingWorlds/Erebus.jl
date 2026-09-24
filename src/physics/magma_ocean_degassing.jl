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
- `T_mo`: Magma ocean / surface interface temperature [K].
- `delta_IW`: Oxygen fugacity offset relative to iron-wüstite buffer [log10 units].

# Keyword Arguments
- `water_As`: Burnham/Dixon water solubility coefficient [wt% / MPa^0.5] (default: 0.40).
- `water_law`: Water solubility formulation (`:burnham_dixon`, `:sossi_peridotite`, etc.) (default: `:burnham_dixon`).
- `carbon_active`: Whether carbon dissolves in melt (default: true).
- `sulfur_active`: Whether sulfur dissolves in melt (default: true).
- `co_law`: CO solubility formulation (default: `:armstrong2015`).
- `ch4_law`: CH4 solubility formulation (default: `:ardia2013`).
- `co2_law`: CO2 solubility formulation (default: `:dixon1995`).
- `sulfide_law`: Sulfide solubility formulation (default: `:boulliung2023`).
- `nitrogen_henry`: Physical N2 Henry coefficient [ppm / bar] (default: 0.40).
- `nitrogen_nitride`: Chemical nitride capacity [wt% / bar^0.5] (default: 1.0e-3).
- `graphite_saturation`: Whether to cap carbon fugacities at graphite saturation (default: true).

# Returns
- `NamedTuple`:
  - `P_surf`: Total equilibrium surface atmospheric pressure [Pa].
  - `p_i`: Dictionary of species partial pressures [Pa].
  - `M_atm_i`: Dictionary of species atmospheric masses [kg].
  - `M_atm_tot`: Total atmospheric mass [kg].
  - `M_melt_H`: Elemental hydrogen mass retained in silicate melt [kg].
  - `M_melt_C`: Elemental carbon mass retained in silicate melt [kg].
  - `M_melt_N`: Elemental nitrogen mass retained in silicate melt [kg].
  - `M_melt_S`: Elemental sulfur mass retained in silicate melt [kg].
  - `M_melt_tot`: Total volatile mass dissolved in melt [kg].
  - `w_diss_H2O_wtpct`: Dissolved water concentration in melt [wt%].
  - `C_diss_C_ppm`: Dissolved carbon concentration in melt [ppmw].
  - `C_diss_N_ppm`: Dissolved nitrogen concentration in melt [ppmw].
  - `C_diss_S_ppm`: Dissolved sulfur concentration in melt [ppmw].
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
)::@NamedTuple{
    P_surf::Float64,
    p_i::Dict{Symbol,Float64},
    M_atm_i::Dict{Symbol,Float64},
    M_atm_tot::Float64,
    M_melt_H::Float64,
    M_melt_C::Float64,
    M_melt_N::Float64,
    M_melt_S::Float64,
    M_melt_tot::Float64,
    w_diss_H2O_wtpct::Float64,
    C_diss_C_ppm::Float64,
    C_diss_N_ppm::Float64,
    C_diss_S_ppm::Float64,
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

    total_volatile_mass = mH + mC + mN + mS
    area = 4.0 * π * (Rp^2)
    col_coeff = area / grav # kg / Pa

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

    if total_volatile_mass == 0.0
        return (
            P_surf=0.0,
            p_i=empty_p_i,
            M_atm_i=copy(empty_p_i),
            M_atm_tot=0.0,
            M_melt_H=0.0,
            M_melt_C=0.0,
            M_melt_N=0.0,
            M_melt_S=0.0,
            M_melt_tot=0.0,
            w_diss_H2O_wtpct=0.0,
            C_diss_C_ppm=0.0,
            C_diss_N_ppm=0.0,
            C_diss_S_ppm=0.0,
        )
    end

    # Elemental molar inventories
    nH_tot = mH / 1.008e-3
    nC_tot = mC / 12.011e-3
    nN_tot = mN / 14.007e-3
    nS_tot = mS / 32.060e-3
    n_sum = nH_tot + nC_tot + nN_tot + nS_tot

    z_H = n_sum > 0.0 ? nH_tot / n_sum : 0.80
    z_C = n_sum > 0.0 ? nC_tot / n_sum : 0.15
    z_N = n_sum > 0.0 ? nN_tot / n_sum : 0.03
    z_S = n_sum > 0.0 ? nS_tot / n_sum : 0.02

    # Asymptotic check: If no melt, atmosphere holds 100% of volatiles
    if M_m == 0.0
        p_surf_pure = total_volatile_mass / col_coeff
        spec = solve_chnos_speciation(
            max(p_surf_pure, 1.0),
            T,
            dIW;
            z_H=z_H,
            z_C=z_C,
            z_N=z_N,
            z_S=z_S,
            graphite_saturation=graphite_saturation,
        )
        p_dict = Dict{Symbol,Float64}(
            :H2 => spec.p_H2_Pa,
            :H2O => spec.p_H2O_Pa,
            :CO => spec.p_CO_Pa,
            :CO2 => spec.p_CO2_Pa,
            :CH4 => spec.p_CH4_Pa,
            :N2 => spec.p_N2_Pa,
            :NH3 => spec.p_NH3_Pa,
            :H2S => spec.p_H2S_Pa,
            :S2 => spec.p_S2_Pa,
            :SO2 => spec.p_SO2_Pa,
        )
        # Normalize partial pressures to sum to p_surf_pure
        p_tot_spec = sum(values(p_dict))
        if p_tot_spec > 0.0
            for k in keys(p_dict)
                p_dict[k] = p_dict[k] * (p_surf_pure / p_tot_spec)
            end
        end
        m_atm_dict = Dict{Symbol,Float64}(k => v * col_coeff for (k, v) in p_dict)
        return (
            P_surf=p_surf_pure,
            p_i=p_dict,
            M_atm_i=m_atm_dict,
            M_atm_tot=total_volatile_mass,
            M_melt_H=0.0,
            M_melt_C=0.0,
            M_melt_N=0.0,
            M_melt_S=0.0,
            M_melt_tot=0.0,
            w_diss_H2O_wtpct=0.0,
            C_diss_C_ppm=0.0,
            C_diss_N_ppm=0.0,
            C_diss_S_ppm=0.0,
        )
    end

    # Helper function: evaluate mass of volatiles at a given trial surface pressure
    function eval_partitioning_at_pressure(
        P_trial::Float64, zH::Float64, zC::Float64, zN::Float64, zS::Float64
    )
        P_eval = max(P_trial, 1.0e-3)
        spec = solve_chnos_speciation(
            P_eval,
            T,
            dIW;
            z_H=zH,
            z_C=zC,
            z_N=zN,
            z_S=zS,
            graphite_saturation=graphite_saturation,
        )
        p_dict = Dict{Symbol,Float64}(
            :H2 => spec.p_H2_Pa,
            :H2O => spec.p_H2O_Pa,
            :CO => spec.p_CO_Pa,
            :CO2 => spec.p_CO2_Pa,
            :CH4 => spec.p_CH4_Pa,
            :N2 => spec.p_N2_Pa,
            :NH3 => spec.p_NH3_Pa,
            :H2S => spec.p_H2S_Pa,
            :S2 => spec.p_S2_Pa,
            :SO2 => spec.p_SO2_Pa,
        )
        # Rescale partial pressures to sum exactly to P_trial
        p_tot_spec = sum(values(p_dict))
        if p_tot_spec > 0.0
            for k in keys(p_dict)
                p_dict[k] = p_dict[k] * (P_trial / p_tot_spec)
            end
        end

        m_atm_dict = Dict{Symbol,Float64}(k => v * col_coeff for (k, v) in p_dict)

        # Atmospheric elemental mass contributions
        m_atm_H = (
            m_atm_dict[:H2] * 1.0 +
            m_atm_dict[:H2O] * (2.01588 / 18.01528) +
            m_atm_dict[:CH4] * (4.03176 / 16.04246) +
            m_atm_dict[:NH3] * (3.02382 / 17.03052) +
            m_atm_dict[:H2S] * (2.01588 / 34.08088)
        )
        m_atm_C = (
            m_atm_dict[:CO] * (12.011 / 28.0101) +
            m_atm_dict[:CO2] * (12.011 / 44.0095) +
            m_atm_dict[:CH4] * (12.011 / 16.04246)
        )
        m_atm_N = (m_atm_dict[:N2] * 1.0 + m_atm_dict[:NH3] * (14.007 / 17.03052))
        m_atm_S = (
            m_atm_dict[:H2S] * (32.060 / 34.08088) +
            m_atm_dict[:SO2] * (32.060 / 64.066) +
            m_atm_dict[:S2] * 1.0
        )

        # Melt dissolved volatile mass fractions
        # 1. Water & H2
        w_H2O_wtpct = compute_water_solubility_melt(
            p_dict[:H2O]; As=water_As, law=water_law
        )
        w_H2O_frac = w_H2O_wtpct * 0.01
        w_H2_frac = compute_h2_solubility_melt(p_dict[:H2]) * 1.0e-6
        w_diss_H = w_H2O_frac * (2.01588 / 18.01528) + w_H2_frac
        m_melt_H = M_m * w_diss_H

        # 2. Nitrogen
        S_N_res = compute_nitrogen_solubility_melt(
            p_dict[:N2], dIW; Kh=nitrogen_henry, C_nitride=nitrogen_nitride
        )
        C_N_ppm = S_N_res.total_ppm
        w_diss_N = C_N_ppm * 1.0e-6
        m_melt_N = M_m * w_diss_N

        # 3. Carbon
        C_C_ppm = 0.0
        if carbon_active
            S_C_res = compute_carbon_solubility_melt(
                P_eval,
                T,
                dIW;
                co_law=co_law,
                ch4_law=ch4_law,
                co2_law=co2_law,
                graphite_saturation=graphite_saturation,
            )
            # Scale dissolved C to the actual carbon gas mole fraction in the atmosphere
            f_C_gas = clamp((p_dict[:CO] + p_dict[:CO2] + p_dict[:CH4]) / P_eval, 0.0, 1.0)
            C_C_ppm = S_C_res.total_ppm * f_C_gas
        end
        w_diss_C = C_C_ppm * 1.0e-6
        m_melt_C = M_m * w_diss_C

        # 4. Sulfur
        C_S_ppm = 0.0
        if sulfur_active
            C_S_ppm = compute_sulfur_solubility_melt(p_dict[:S2], T, dIW; law=sulfide_law)
        end
        w_diss_S = C_S_ppm * 1.0e-6
        m_melt_S = M_m * w_diss_S

        m_tot_calc =
            (m_atm_H + m_melt_H) +
            (m_atm_C + m_melt_C) +
            (m_atm_N + m_melt_N) +
            (m_atm_S + m_melt_S)
        return (
            m_tot_calc=m_tot_calc,
            p_dict=p_dict,
            m_atm_dict=m_atm_dict,
            m_atm_H=m_atm_H,
            m_atm_C=m_atm_C,
            m_atm_N=m_atm_N,
            m_atm_S=m_atm_S,
            m_melt_H=m_melt_H,
            m_melt_C=m_melt_C,
            m_melt_N=m_melt_N,
            m_melt_S=m_melt_S,
            w_H2O_wtpct=w_H2O_wtpct,
            C_C_ppm=C_C_ppm,
            C_N_ppm=C_N_ppm,
            C_S_ppm=C_S_ppm,
        )
    end

    # Outer iteration to converge elemental atmospheric gas fractions
    P_low = 1.0e-4
    P_high = max(1.0e8, 10.0 * (total_volatile_mass / col_coeff))
    best_res = nothing

    for outer_iter in 1:8
        # Monotonic 1D bisection/Brent search for P_surf
        P_a = P_low
        P_b = P_high

        # Ensure bracket
        res_b = eval_partitioning_at_pressure(P_b, z_H, z_C, z_N, z_S)
        while res_b.m_tot_calc < total_volatile_mass && P_b < 1.0e11
            P_b *= 10.0
            res_b = eval_partitioning_at_pressure(P_b, z_H, z_C, z_N, z_S)
        end

        for _ in 1:60
            P_mid = 0.5 * (P_a + P_b)
            res_mid = eval_partitioning_at_pressure(P_mid, z_H, z_C, z_N, z_S)
            diff = res_mid.m_tot_calc - total_volatile_mass
            if abs(diff) / total_volatile_mass < 1.0e-7 || (P_b - P_a) / P_mid < 1.0e-7
                best_res = res_mid
                break
            end
            if diff > 0.0
                P_b = P_mid
            else
                P_a = P_mid
            end
            best_res = res_mid
        end

        # Update elemental fractions in gas for next outer iteration based on target exsolved inventories
        nH_atm = max(0.0, mH - best_res.m_melt_H) / 1.008e-3
        nC_atm = max(0.0, mC - best_res.m_melt_C) / 12.011e-3
        nN_atm = max(0.0, mN - best_res.m_melt_N) / 14.007e-3
        nS_atm = max(0.0, mS - best_res.m_melt_S) / 32.060e-3
        n_atm_tot = nH_atm + nC_atm + nN_atm + nS_atm
        if n_atm_tot > 0.0
            z_H_new = nH_atm / n_atm_tot
            z_C_new = nC_atm / n_atm_tot
            z_N_new = nN_atm / n_atm_tot
            z_S_new = nS_atm / n_atm_tot
            if max(
                abs(z_H_new - z_H),
                abs(z_C_new - z_C),
                abs(z_N_new - z_N),
                abs(z_S_new - z_S),
            ) < 1.0e-5
                break
            end
            z_H = 0.5 * (z_H + z_H_new)
            z_C = 0.5 * (z_C + z_C_new)
            z_N = 0.5 * (z_N + z_N_new)
            z_S = 0.5 * (z_S + z_S_new)
        else
            break
        end
    end

    # Enforce strict conservation of elemental mass across melt and atmosphere
    final_atm_i = copy(best_res.m_atm_dict)
    final_p_i = Dict{Symbol,Float64}(k => v / col_coeff for (k, v) in final_atm_i)
    final_P_surf = sum(values(final_p_i))
    final_atm_tot = sum(values(final_atm_i))

    m_atm_H = (
        get(final_atm_i, :H2, 0.0) * 1.0 +
        get(final_atm_i, :H2O, 0.0) * (2.01588 / 18.01528) +
        get(final_atm_i, :CH4, 0.0) * (4.03176 / 16.04246) +
        get(final_atm_i, :NH3, 0.0) * (3.02382 / 17.03052) +
        get(final_atm_i, :H2S, 0.0) * (2.01588 / 34.08088)
    )
    m_atm_C = (
        get(final_atm_i, :CO, 0.0) * (12.011 / 28.0101) +
        get(final_atm_i, :CO2, 0.0) * (12.011 / 44.0095) +
        get(final_atm_i, :CH4, 0.0) * (12.011 / 16.04246)
    )
    m_atm_N = (
        get(final_atm_i, :N2, 0.0) * 1.0 + get(final_atm_i, :NH3, 0.0) * (14.007 / 17.03052)
    )
    m_atm_S = (
        get(final_atm_i, :H2S, 0.0) * (32.060 / 34.08088) +
        get(final_atm_i, :SO2, 0.0) * (32.060 / 64.066) +
        get(final_atm_i, :S2, 0.0) * 1.0
    )

    scale_H = m_atm_H > mH ? mH / max(m_atm_H, 1.0e-30) : 1.0
    scale_C = m_atm_C > mC ? mC / max(m_atm_C, 1.0e-30) : 1.0
    scale_N = m_atm_N > mN ? mN / max(m_atm_N, 1.0e-30) : 1.0
    scale_S = m_atm_S > mS ? mS / max(m_atm_S, 1.0e-30) : 1.0

    if min(scale_H, scale_C, scale_N, scale_S) < 1.0
        for (sp, mass) in final_atm_i
            if sp in (:H2, :H2O)
                final_atm_i[sp] = mass * scale_H
            elseif sp in (:CO, :CO2)
                final_atm_i[sp] = mass * scale_C
            elseif sp === :CH4
                final_atm_i[sp] = mass * min(scale_H, scale_C)
            elseif sp === :N2
                final_atm_i[sp] = mass * scale_N
            elseif sp === :NH3
                final_atm_i[sp] = mass * min(scale_H, scale_N)
            elseif sp in (:SO2, :S2)
                final_atm_i[sp] = mass * scale_S
            elseif sp === :H2S
                final_atm_i[sp] = mass * min(scale_H, scale_S)
            end
        end
        m_atm_H = (
            get(final_atm_i, :H2, 0.0) * 1.0 +
            get(final_atm_i, :H2O, 0.0) * (2.01588 / 18.01528) +
            get(final_atm_i, :CH4, 0.0) * (4.03176 / 16.04246) +
            get(final_atm_i, :NH3, 0.0) * (3.02382 / 17.03052) +
            get(final_atm_i, :H2S, 0.0) * (2.01588 / 34.08088)
        )
        m_atm_C = (
            get(final_atm_i, :CO, 0.0) * (12.011 / 28.0101) +
            get(final_atm_i, :CO2, 0.0) * (12.011 / 44.0095) +
            get(final_atm_i, :CH4, 0.0) * (12.011 / 16.04246)
        )
        m_atm_N = (
            get(final_atm_i, :N2, 0.0) * 1.0 +
            get(final_atm_i, :NH3, 0.0) * (14.007 / 17.03052)
        )
        m_atm_S = (
            get(final_atm_i, :H2S, 0.0) * (32.060 / 34.08088) +
            get(final_atm_i, :SO2, 0.0) * (32.060 / 64.066) +
            get(final_atm_i, :S2, 0.0) * 1.0
        )
        final_p_i = Dict{Symbol,Float64}(k => v / col_coeff for (k, v) in final_atm_i)
        final_P_surf = sum(values(final_p_i))
        final_atm_tot = sum(values(final_atm_i))
    end

    final_melt_H = max(0.0, mH - m_atm_H)
    final_melt_C = max(0.0, mC - m_atm_C)
    final_melt_N = max(0.0, mN - m_atm_N)
    final_melt_S = max(0.0, mS - m_atm_S)
    final_melt_tot = final_melt_H + final_melt_C + final_melt_N + final_melt_S

    w_diss_H2O = M_m > 0.0 ? (final_melt_H * (18.01528 / 2.01588) / M_m) * 100.0 : 0.0
    C_diss_C = M_m > 0.0 ? (final_melt_C / M_m) * 1.0e6 : 0.0
    C_diss_N = M_m > 0.0 ? (final_melt_N / M_m) * 1.0e6 : 0.0
    C_diss_S = M_m > 0.0 ? (final_melt_S / M_m) * 1.0e6 : 0.0

    return (
        P_surf=final_P_surf,
        p_i=final_p_i,
        M_atm_i=final_atm_i,
        M_atm_tot=final_atm_tot,
        M_melt_H=final_melt_H,
        M_melt_C=final_melt_C,
        M_melt_N=final_melt_N,
        M_melt_S=final_melt_S,
        M_melt_tot=final_melt_tot,
        w_diss_H2O_wtpct=w_diss_H2O,
        C_diss_C_ppm=C_diss_C,
        C_diss_N_ppm=C_diss_N,
        C_diss_S_ppm=C_diss_S,
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
        ex_H2O = min(max(0.0, w_H2O_m - ret_H2O), max(0.0, w_H2O_m / F_curr - w_H2O_sat) * F_curr * eff)
        ex_C = min(max(0.0, w_C_m - ret_C), max(0.0, w_C_m / F_curr - w_C_sat) * F_curr * eff)
        ex_N = min(max(0.0, w_N_m - ret_N), max(0.0, w_N_m / F_curr - w_N_sat) * F_curr * eff)
        ex_S = min(max(0.0, w_S_m - ret_S), max(0.0, w_S_m / F_curr - w_S_sat) * F_curr * eff)

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

# Benchmark generator for Erebus.jl iron core formation and planetesimal differentiation
using JSON
using LinearAlgebra
using StaticArrays

push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..")))
using Erebus
using Erebus.Physics
using Erebus.Config

println("=== Erebus.jl Core Formation Benchmark Generation ===")

# -----------------------------------------------------------------------------
# Planetesimal Physics and Simulation Setup
# -----------------------------------------------------------------------------
const R_PLANET = 50_000.0        # 50 km radius planetesimal
const T_SURF = 150.0             # Surface and initial primordial temperature [K]
const PHI_ICE_INIT = 0.30        # Initial pore ice volume fraction (30 vol%)
const PHI_FE_INIT = 0.12         # Initial chondritic metal volume fraction (12 vol%)
const PHI_ROCK_INIT = 1.0 - PHI_ICE_INIT - PHI_FE_INIT # 0.58 silicate rock matrix

const T_ICE_MELT = 273.15        # Pore ice melting point [K]
const T_EUTECTIC = 1213.0        # Fe-FeS eutectic melting point [K]
const DT_METAL = 50.0            # Melting interval for metal [K]
const T_SOLIDUS = 1416.0         # Silicate solidus [K]
const T_LIQUIDUS = 1800.0        # Silicate liquidus [K]
const L_SILICATE = 4.0e5         # Silicate latent heat of melting [J/kg]
const L_ICE = 3.33e5             # Water ice latent heat of melting [J/kg]
const L_METAL = 2.7e5            # Fe-FeS eutectic latent heat of melting [J/kg]

const RHO_SIL_SOLID = 3300.0     # Silicate solid density [kg/m^3]
const RHO_SIL_MELT = 2700.0      # Silicate melt density [kg/m^3]
const RHO_METAL_SOLID = 5700.0   # Solid Fe-FeS eutectic density [kg/m^3]
const RHO_METAL_MELT = 5450.0    # Molten Fe-FeS eutectic density [kg/m^3]
const RHO_ICE = 920.0            # Water ice density [kg/m^3]

const CP_SIL = 1000.0            # Silicate heat capacity [J/(kg K)]
const K_COND_SIL = 3.0           # Silicate thermal conductivity [W/(m K)]
const ETA_SOLID_MATRIX = 1.0e19  # Cold rock effective viscosity [Pa s]
const ETA_METAL = 1.0e-2         # Liquid Fe-FeS dynamic viscosity [Pa s]
const SIGMA_FE_SIL = 0.5         # Interfacial surface tension [N/m]
const PHI_PACK = 0.65            # Maximum droplet packing fraction
const PHI_RESIDUAL = 0.02        # Residual metal trapped in pores
const PHI_CRIT_PERC = 0.05       # Critical percolation threshold
const K_REF_METAL = 1.0e-9       # Reference permeability [m^2]
const PHI0_PERC = 0.10           # Reference percolation porosity scale [-]
const PERM_EXPONENT = 3.0        # Kozeny-Carman permeability exponent [-]
const HINDERED_EXPONENT = 4.5    # Richardson-Zaki hindered settling exponent [-]
const F_SETTLE_START = 0.40      # Melt fraction at percolation-settling transition start [-]
const F_PERC_END = 0.50          # Melt fraction where percolation ends [-]
const CFL_1D = 0.5               # Advective CFL safety factor [-]
const MAX_SUBCYCLES_1D = 2000    # Maximum segregation substeps per outer step [-]

# Radiogenic heating (26Al)
const TAU_26AL = 1.034e6 * 365.25 * 86400.0 # 1.034 Ma mean lifetime [s]
const Q0_26AL = 1.10e-7                      # Radiogenic heating rate per kg of rock [W/kg] (accretion ~1.4 Ma after CAIs)

# Gravitational acceleration profile g(r) = (4/3) pi G rho_mean r
const G_CONST = 6.6743e-11
const RHO_MEAN = 3200.0
const G_SURF = (4.0 / 3.0) * pi * G_CONST * RHO_MEAN * R_PLANET # ~ 0.0447 m/s^2

"""
Evaluate `metal_segregation_velocity` using shared 1-D benchmark configuration parameters.
"""
function _v_seg_1d(
    phi_m_val::Float64,
    Fm_val::Float64,
    drho::Float64,
    g_val::Float64,
    eta_val::Float64,
    r_drop::Float64,
    cfg::CoreFormationConfig;
    perc::Bool=cfg.percolation_active,
    settle::Bool=cfg.settling_active,
)::Float64
    return metal_segregation_velocity(
        phi_m_val,
        Fm_val,
        drho,
        g_val,
        eta_val;
        percolation_active=perc,
        settling_active=settle,
        k_metal_ref=cfg.k_metal_ref,
        eta_metal=cfg.eta_metal,
        phi_crit_perc=cfg.phi_crit_perc,
        phi_residual=cfg.phi_residual,
        phi0=cfg.phi0,
        perm_exponent=cfg.perm_exponent,
        r_drop=r_drop,
        hindered_exponent=cfg.hindered_exponent,
        phi_pack=cfg.phi_pack,
        hadamard_rybczynski=cfg.hadamard_rybczynski,
        F_settle_start=cfg.F_settle_start,
        F_perc_end=cfg.F_perc_end,
    )
end

"""
Return smooth cubic Hermite regime blend weight (0 in percolation, 1 in settling).
"""
function _settling_blend_1d(F_m::Float64, cfg::CoreFormationConfig)::Float64
    if cfg.percolation_active && !cfg.settling_active
        return 0.0
    elseif !cfg.percolation_active && cfg.settling_active
        return 1.0
    elseif !cfg.percolation_active && !cfg.settling_active
        return 0.0
    elseif F_m <= cfg.F_settle_start
        return 0.0
    elseif F_m >= cfg.F_perc_end
        return 1.0
    else
        s = (F_m - cfg.F_settle_start) / (cfg.F_perc_end - cfg.F_settle_start)
        return s * s * (3.0 - 2.0 * s)
    end
end

"""
Return Darcy percolation characteristic wave speed |dq_perc/dphi_m| [m/s].
"""
function _perc_wave_speed_1d(
    phi_m_val::Float64, v_perc_val::Float64, cfg::CoreFormationConfig
)::Float64
    v_perc_val <= 0.0 && return 0.0
    phi_perm = min(clamp(phi_m_val, 0.0, 1.0), 1.0 - 1.0e-7)
    phi_mob_p = phi_perm - cfg.phi_crit_perc
    phi_mob_p <= 0.0 && return 0.0
    phi_mob_eff = 0.0 < cfg.perm_exponent < 1.0 ? max(phi_mob_p, 1.0e-7) : phi_mob_p
    return v_perc_val *
           phi_perm *
           (cfg.perm_exponent / phi_mob_eff + 2.0 / (1.0 - phi_mob_p))
end

"""
Compute 1-D Godunov hindered-settling flux and maximum characteristic wave speed [m/s].
"""
function _godunov_settling_1d(
    phi_D::Float64, phi_R::Float64, v0_D::Float64, m_exp::Float64, phi_pack::Float64
)::Tuple{Float64,Float64}
    if v0_D <= 0.0 || phi_pack <= 0.0
        return (0.0, 0.0)
    end
    pd = clamp(phi_D, 0.0, phi_pack)
    pr = clamp(phi_R, 0.0, phi_pack)
    f_val = (p) -> (p <= 0.0 || p >= phi_pack) ? 0.0 : v0_D * p * (1.0 - p / phi_pack)^m_exp
    df_abs =
        (p) -> begin
            if p < 0.0 || p > phi_pack
                return 0.0
            elseif iszero(m_exp)
                return p < phi_pack ? v0_D : 0.0
            end
            u = 0.0 < m_exp < 1.0 ? min(p / phi_pack, 1.0 - 1.0e-7) : (p / phi_pack)
            if p >= phi_pack && m_exp >= 1.0
                return isapprox(m_exp, 1.0; atol=1.0e-12) ? v0_D : 0.0
            end
            return v0_D * (max(1.0 - u, 0.0)^(m_exp - 1.0)) * abs(1.0 - (m_exp + 1.0) * u)
        end
    q_d = f_val(pd)
    q_r = f_val(pr)
    q_settle = if pd <= pr
        min(q_d, q_r)
    else
        qm = max(q_d, q_r)
        p_star = phi_pack / (m_exp + 1.0)
        (pr < p_star < pd && 0.0 < p_star < phi_pack) ? max(qm, f_val(p_star)) : qm
    end
    pa = min(pd, pr)
    pb = max(pd, pr)
    wmax = max(df_abs(pa), df_abs(pb))
    if m_exp > 1.0
        p_inflect = 2.0 * phi_pack / (m_exp + 1.0)
        if pa < p_inflect < pb
            wmax = max(wmax, df_abs(p_inflect))
        end
    end
    return (q_settle, wmax)
end

"""
Run a 1D spherical core formation benchmark simulation.
Resolves thermal conduction, 26Al decay, ice/rock melting latent heats,
metal melting, porous percolation, Stokes droplet settling,
viscosity weakening, core ponding, and dissipation heating.
"""
function run_core_formation_sim(
    Nr::Int=128;
    percolation_active::Bool=true,
    settling_active::Bool=true,
    droplet_mode::Symbol=:capillary_mean, # :fixed, :capillary_mean, :bond_mean, :weber_mean, :weber_turbulent
    segregation_heating::Bool=true,
    droplet_diameter_fixed::Float64=1.0e-2,
    t_end_yr::Float64=3.5e6,          # 3.5 Myr
    dt_yr::Float64=500.0,             # 500 yr initial timestep
    save_snapshots::Bool=true,
    n_save_snaps::Int=70,
)
    dr = R_PLANET / Nr
    r_centers = [(i - 0.5) * dr for i in 1:Nr]
    r_faces = [i * dr for i in 0:Nr]
    vol_shells = [(4.0 / 3.0) * pi * (r_faces[i + 1]^3 - r_faces[i]^3) for i in 1:Nr]
    areas = [4.0 * pi * r_faces[i]^2 for i in 1:(Nr + 1)]

    # Initial state: homogeneous cold primordial mixture of ice, rock, and metal
    T = fill(T_SURF, Nr)
    phi_fe = fill(PHI_FE_INIT, Nr)
    phi_ice = fill(PHI_ICE_INIT, Nr)
    phi_rock = fill(PHI_ROCK_INIT, Nr)

    sec_per_yr = 365.25 * 86400.0
    t_sec = 0.0
    t_yr = 0.0
    dt_sec = dt_yr * sec_per_yr

    # Diagnostics history
    times_yr = Float64[]
    T_core_hist = Float64[]
    R_core_hist = Float64[]
    R_magma_hist = Float64[]
    R_ice_hist = Float64[]
    v_seg_peak_hist = Float64[]
    E_diss_cum_hist = Float64[]
    phi_fe_core_hist = Float64[]

    snapshots_T = Dict{String,Vector{Float64}}()
    snapshots_phi_fe = Dict{String,Vector{Float64}}()
    snapshots_phi_ice = Dict{String,Vector{Float64}}()
    snapshots_Fm = Dict{String,Vector{Float64}}()
    snapshots_vseg = Dict{String,Vector{Float64}}()
    snapshots_Qseg = Dict{String,Vector{Float64}}()
    snapshots_rho = Dict{String,Vector{Float64}}()

    t_save_interval = t_end_yr / n_save_snaps
    next_save_t = 0.0
    total_diss_energy = 0.0
    M_fe_init = sum(phi_fe[i] * vol_shells[i] for i in 1:Nr)
    cfg_1d = CoreFormationConfig(;
        percolation_active=percolation_active,
        settling_active=settling_active,
        k_metal_ref=K_REF_METAL,
        eta_metal=ETA_METAL,
        phi_crit_perc=PHI_CRIT_PERC,
        phi_residual=PHI_RESIDUAL,
        phi0=PHI0_PERC,
        perm_exponent=PERM_EXPONENT,
        hindered_exponent=HINDERED_EXPONENT,
        phi_pack=PHI_PACK,
        hadamard_rybczynski=true,
        F_settle_start=F_SETTLE_START,
        F_perc_end=F_PERC_END,
    )

    step = 0
    while t_yr <= t_end_yr
        step += 1

        # 1. Physical evaluations at cell centers
        g_r = [G_SURF * (r_centers[i] / R_PLANET) for i in 1:Nr]

        # Silicate melt fraction Fm
        Fm = zeros(Float64, Nr)
        for i in 1:Nr
            if T[i] < T_SOLIDUS
                Fm[i] = 0.0
            elseif T[i] >= T_LIQUIDUS
                Fm[i] = 1.0
            else
                Fm[i] = (T[i] - T_SOLIDUS) / (T_LIQUIDUS - T_SOLIDUS)
            end
        end

        # Metal melt fraction chi_fe and mobile metal phi_m
        chi_fe = [
            compute_metal_melt_fraction(T[i]; T_eutectic=T_EUTECTIC, dT_metal=DT_METAL) for
            i in 1:Nr
        ]
        phi_m = chi_fe .* phi_fe

        # Dynamic pore ice melting & water desiccation
        for i in 1:Nr
            if T[i] >= T_ICE_MELT
                phi_ice[i] = max(0.0, phi_ice[i] - 0.01)
            end
        end

        # Effective viscosity
        eta_matrix = zeros(Float64, Nr)
        r_drop_arr = zeros(Float64, Nr)
        v_seg = zeros(Float64, Nr)
        v_stokes0_r = zeros(Float64, Nr)
        Q_seg = zeros(Float64, Nr)
        drho_shell = zeros(Float64, Nr)
        rho_bulk = zeros(Float64, Nr)

        for i in 1:Nr
            # Smooth melt-weakened suspension viscosity (Roscoe / Solomatov)
            eta_matrix[i] = compute_silicate_suspension_viscosity(
                ETA_SOLID_MATRIX,
                Fm[i];
                alpha_eta=28.0,
                phi_crit=0.40,
                eta_melt=10.0,
                etamin=0.01,
                etamax=ETA_SOLID_MATRIX,
            )

            rho_sil_eff = Fm[i] * RHO_SIL_MELT + (1.0 - Fm[i]) * RHO_SIL_SOLID
            drho = max(RHO_METAL_MELT - rho_sil_eff, 1.0)

            # Droplet radius
            r_drop = if droplet_mode === :fixed
                droplet_diameter_fixed / 2.0
            elseif droplet_mode === :capillary_mean ||
                droplet_mode === :bond_mean ||
                droplet_mode === :weber_mean
                d_w = sqrt(10.0 * SIGMA_FE_SIL / max(drho * g_r[i], 1.0e-8))
                clamp(d_w / 2.0, 1.0e-4, 0.05)
            else # :weber_turbulent
                v_est = stokes_settling_velocity(
                    droplet_diameter_fixed / 2.0,
                    drho,
                    max(g_r[i], 1.0e-5),
                    eta_matrix[i],
                )
                v_rel = max(v_est, 1.0e-6)
                d_w = weber_equilibrium_diameter(
                    RHO_SIL_SOLID, v_rel, SIGMA_FE_SIL; We_crit=10.0
                )
                clamp(d_w / 2.0, 1.0e-4, 0.05)
            end
            r_drop_arr[i] = r_drop

            v_seg[i] = _v_seg_1d(
                phi_m[i], Fm[i], drho, g_r[i], eta_matrix[i], r_drop, cfg_1d
            )
            v_stokes0_r[i] =
                if settling_active && Fm[i] >= cfg_1d.F_settle_start && g_r[i] > 0.0
                    stokes_settling_velocity(
                        r_drop,
                        drho,
                        g_r[i],
                        eta_matrix[i];
                        hadamard_rybczynski=cfg_1d.hadamard_rybczynski,
                        eta_metal=cfg_1d.eta_metal,
                    )
                else
                    0.0
                end

            drho_shell[i] = drho

            # Bulk density blending
            rho_sil_local = Fm[i] * RHO_SIL_MELT + (1.0 - Fm[i]) * RHO_SIL_SOLID
            rho_metal_local =
                chi_fe[i] * RHO_METAL_MELT + (1.0 - chi_fe[i]) * RHO_METAL_SOLID
            rho_bulk[i] =
                phi_fe[i] * rho_metal_local +
                (1.0 - phi_fe[i] - phi_ice[i]) * rho_sil_local +
                phi_ice[i] * RHO_ICE
        end

        # 2. Segregation transport (subcycled conservative finite-volume inward flux)
        max_v = maximum(v_seg)
        if percolation_active || settling_active
            F_flux = zeros(Float64, Nr + 1)
            req_rate = zeros(Float64, Nr + 1)
            w_s_arr = [_settling_blend_1d(Fm[i], cfg_1d) for i in 1:Nr]
            phi_trap_arr = [
                effective_metal_trap_fraction(
                    Fm[i],
                    PHI_RESIDUAL;
                    percolation_active=percolation_active,
                    settling_active=settling_active,
                    F_settle_start=cfg_1d.F_settle_start,
                    F_perc_end=cfg_1d.F_perc_end,
                ) for i in 1:Nr
            ]
            t_rem = dt_sec
            sub = 0
            while t_rem > 1.0e-12 * dt_sec && sub < MAX_SUBCYCLES_1D
                fill!(req_rate, 0.0)
                max_v_char = 0.0
                for i in 2:Nr
                    donor = i
                    receiver = i - 1
                    phi_m_D = phi_m[donor]
                    phi_mob = max(phi_m_D - phi_trap_arr[donor], 0.0)
                    cap_R = max(PHI_PACK - phi_fe[receiver], 0.0) * vol_shells[receiver]
                    if phi_mob <= 1.0e-15 || cap_R <= 1.0e-12 * vol_shells[receiver]
                        continue
                    end
                    w_s_D = w_s_arr[donor]
                    v_perc_D = if percolation_active && w_s_D < 1.0
                        _v_seg_1d(
                            phi_m_D,
                            Fm[donor],
                            drho_shell[donor],
                            g_r[donor],
                            eta_matrix[donor],
                            r_drop_arr[donor],
                            cfg_1d;
                            perc=true,
                            settle=false,
                        )
                    else
                        0.0
                    end
                    q_perc = (1.0 - w_s_D) * v_perc_D * phi_m_D
                    q_settle, c_settle = if w_s_D > 0.0
                        _godunov_settling_1d(
                            phi_m_D,
                            phi_fe[receiver],
                            v_stokes0_r[donor],
                            cfg_1d.hindered_exponent,
                            cfg_1d.phi_pack,
                        )
                    else
                        (0.0, 0.0)
                    end
                    q_tot = q_perc + w_s_D * q_settle
                    q_tot <= 0.0 && continue
                    c_perc = (1.0 - w_s_D) * _perc_wave_speed_1d(phi_m_D, v_perc_D, cfg_1d)
                    v_seg_D = _v_seg_1d(
                        phi_m_D,
                        Fm[donor],
                        drho_shell[donor],
                        g_r[donor],
                        eta_matrix[donor],
                        r_drop_arr[donor],
                        cfg_1d,
                    )
                    max_v_char = max(max_v_char, max(v_seg_D, c_perc + w_s_D * c_settle))
                    req_rate[i] = areas[i] * q_tot
                end
                (max_v_char <= 0.0 || !any(>(0.0), req_rate)) && break
                sub += 1
                dt_sub = min(t_rem, (CFL_1D * dr) / max_v_char)
                fill!(F_flux, 0.0)
                for i in 2:Nr
                    req_rate[i] <= 0.0 && continue
                    donor = i
                    receiver = i - 1
                    avail = max(phi_m[donor] - phi_trap_arr[donor], 0.0) * vol_shells[donor]
                    cap = max(PHI_PACK - phi_fe[receiver], 0.0) * vol_shells[receiver]
                    F_flux[i] = min(req_rate[i] * dt_sub, cap, avail)
                end

                if segregation_heating
                    for i in 1:Nr
                        moved = 0.5 * (F_flux[i] + F_flux[i + 1])
                        q_d = drho_shell[i] * g_r[i] * moved * dr / vol_shells[i] / dt_sec
                        Q_seg[i] += q_d
                        total_diss_energy += q_d * vol_shells[i] * dt_sec
                    end
                end

                for i in 1:Nr
                    phi_fe[i] += (F_flux[i + 1] - F_flux[i]) / vol_shells[i]
                    @assert -1.0e-12 <= phi_fe[i] <= PHI_PACK + 1.0e-12 "1-D phi_fe[$i] out of bounds: $(phi_fe[i])"
                    phi_fe[i] = clamp(phi_fe[i], 0.0, PHI_PACK)
                    phi_m[i] = chi_fe[i] * phi_fe[i]
                end
                t_rem -= dt_sub
            end
            if sub == MAX_SUBCYCLES_1D && t_rem > 1.0e-9 * dt_sec
                @warn "1-D benchmark reached MAX_SUBCYCLES_1D=$(MAX_SUBCYCLES_1D) with $(round(t_rem / dt_sec; digits=3)) of dt unstepped" maxlog=10
            end
            M_fe_now = sum(phi_fe[i] * vol_shells[i] for i in 1:Nr)
            @assert isapprox(M_fe_now, M_fe_init; rtol=1.0e-11) "1-D metal mass drift: $(M_fe_now) vs $(M_fe_init)"
        end

        # 3. Thermal solve (Implicit spherical diffusion with 26Al decay and Q_seg)
        H_rad_vol = (phi_rock[1] * RHO_SIL_SOLID) * (Q0_26AL * exp(-t_sec / TAU_26AL))

        # Construct tridiagonal system for T^{n+1}
        diag_A = zeros(Float64, Nr)
        sub_A = zeros(Float64, Nr - 1)
        super_A = zeros(Float64, Nr - 1)
        rhs = zeros(Float64, Nr)

        for i in 1:Nr
            rc = r_centers[i]
            rho_eff = rho_bulk[i]
            # Apparent heat capacity buffering
            cp_eff = CP_SIL
            if T[i] >= T_SOLIDUS && T[i] <= T_LIQUIDUS
                cp_eff += L_SILICATE / (T_LIQUIDUS - T_SOLIDUS)
            end
            if T[i] >= T_EUTECTIC && T[i] <= T_EUTECTIC + DT_METAL && phi_fe[i] > 0.0
                cp_eff += (phi_fe[i] * RHO_METAL_SOLID / rho_eff) * (L_METAL / DT_METAL)
            end
            if T[i] >= 270.0 && T[i] <= 276.0 && phi_ice[i] > 0.0
                cp_eff += L_ICE / 6.0
            end

            heat_cap = rho_eff * cp_eff

            vol_shell = vol_shells[i]
            area_w = areas[i]
            area_e = areas[i + 1]

            k_w = K_COND_SIL
            k_e = K_COND_SIL

            # Conduction coefficients
            c_w = (i > 1) ? (area_w * k_w / (dr * vol_shell)) * dt_sec / heat_cap : 0.0
            c_e = (i < Nr) ? (area_e * k_e / (dr * vol_shell)) * dt_sec / heat_cap : 0.0

            diag_A[i] = 1.0 + c_w + c_e
            if i > 1
                sub_A[i - 1] = -c_w
            end
            if i < Nr
                super_A[i] = -c_e
            end

            source_term = (H_rad_vol + Q_seg[i]) * dt_sec / heat_cap
            rhs[i] = T[i] + source_term

            # Surface boundary condition T(R) = T_SURF (distance from cell Nr center to face is dr/2)
            if i == Nr
                c_surf =
                    (areas[Nr + 1] * k_e / ((dr / 2.0) * vol_shell)) * dt_sec / heat_cap
                diag_A[i] += c_surf
                rhs[i] += c_surf * T_SURF
            end
        end

        # Tridiagonal Thomas algorithm
        T_new = copy(rhs)
        c_prime = copy(super_A)
        d_prime = copy(rhs)

        # Forward sweep
        c_prime[1] /= diag_A[1]
        d_prime[1] /= diag_A[1]
        for i in 2:(Nr - 1)
            denom = diag_A[i] - sub_A[i - 1] * c_prime[i - 1]
            c_prime[i] /= denom
            d_prime[i] = (d_prime[i] - sub_A[i - 1] * d_prime[i - 1]) / denom
        end
        denom_end = diag_A[Nr] - sub_A[Nr - 1] * c_prime[Nr - 1]
        d_prime[Nr] = (d_prime[Nr] - sub_A[Nr - 1] * d_prime[Nr - 1]) / denom_end

        # Back substitution
        T_new[Nr] = d_prime[Nr]
        for i in (Nr - 1):-1:1
            T_new[i] = d_prime[i] - c_prime[i] * T_new[i + 1]
        end
        T .= clamp.(T_new, T_SURF, 2500.0)

        # 4. Save diagnostics
        t_sec += dt_sec
        t_yr += dt_yr

        # Core radius: outer boundary where phi_fe >= 0.50
        core_idx = findlast(x -> x >= 0.50, phi_fe)
        r_core_now = core_idx !== nothing ? r_faces[core_idx + 1] / 1000.0 : 0.0

        # Magma ocean radius: where Fm >= 0.40
        magma_idx = findlast(x -> x >= 0.40, Fm)
        r_magma_now = magma_idx !== nothing ? r_faces[magma_idx + 1] / 1000.0 : 0.0

        # Icy shell base radius: innermost radius where phi_ice >= 0.10
        ice_cells = findall(x -> x >= 0.10, phi_ice)
        r_ice_now = !isempty(ice_cells) ? minimum(ice_cells) * dr / 1000.0 : 50.0

        push!(times_yr, t_yr)
        push!(T_core_hist, T[1])
        push!(R_core_hist, r_core_now)
        push!(R_magma_hist, r_magma_now)
        push!(R_ice_hist, r_ice_now)
        push!(v_seg_peak_hist, max_v)
        push!(E_diss_cum_hist, total_diss_energy)
        push!(phi_fe_core_hist, phi_fe[1])

        if save_snapshots && (t_yr >= next_save_t || t_yr >= t_end_yr)
            k_str = string(round(t_yr; digits=1))
            snapshots_T[k_str] = copy(T)
            snapshots_phi_fe[k_str] = copy(phi_fe)
            snapshots_phi_ice[k_str] = copy(phi_ice)
            snapshots_Fm[k_str] = copy(Fm)
            snapshots_vseg[k_str] = copy(v_seg)
            snapshots_Qseg[k_str] = copy(Q_seg)
            snapshots_rho[k_str] = copy(rho_bulk)
            next_save_t += t_save_interval
        end
    end

    return Dict{String,Any}(
        "r" => r_centers,
        "times_yr" => times_yr,
        "T_core_hist" => T_core_hist,
        "R_core_hist" => R_core_hist,
        "R_magma_hist" => R_magma_hist,
        "R_ice_hist" => R_ice_hist,
        "v_seg_peak_hist" => v_seg_peak_hist,
        "E_diss_cum_hist" => E_diss_cum_hist,
        "phi_fe_core_hist" => phi_fe_core_hist,
        "snapshot_T" => snapshots_T,
        "snapshot_phi_fe" => snapshots_phi_fe,
        "snapshot_phi_ice" => snapshots_phi_ice,
        "snapshot_Fm" => snapshots_Fm,
        "snapshot_vseg" => snapshots_vseg,
        "snapshot_Qseg" => snapshots_Qseg,
        "snapshot_rho" => snapshots_rho,
    )
end

# -----------------------------------------------------------------------------
# Execute Benchmark Suite
# -----------------------------------------------------------------------------
println(
    "1. Running Reference Evolutionary Planetesimal Run (Percolation + Settling + Heating)...",
)
res_ref = run_core_formation_sim(
    128;
    percolation_active=true,
    settling_active=true,
    droplet_mode=:capillary_mean,
    segregation_heating=true,
    t_end_yr=3.5e6,
    dt_yr=500.0,
    save_snapshots=true,
    n_save_snaps=80,
)

println("2. Running Regime Sweep: Percolation Only vs Settling Only vs Coupled Hermite...")
res_perc_only = run_core_formation_sim(
    64;
    percolation_active=true,
    settling_active=false,
    droplet_mode=:capillary_mean,
    segregation_heating=true,
    t_end_yr=3.5e6,
    dt_yr=1000.0,
    save_snapshots=false,
)

res_settle_only = run_core_formation_sim(
    64;
    percolation_active=false,
    settling_active=true,
    droplet_mode=:capillary_mean,
    segregation_heating=true,
    t_end_yr=3.5e6,
    dt_yr=1000.0,
    save_snapshots=false,
)

println("3. Running Droplet Physics Sweep: Fixed 1cm vs Weber Mean vs Turbulent Breakup...")
res_drop_fixed = run_core_formation_sim(
    64;
    percolation_active=true,
    settling_active=true,
    droplet_mode=:fixed,
    droplet_diameter_fixed=1.0e-2,
    segregation_heating=true,
    t_end_yr=3.5e6,
    dt_yr=1000.0,
    save_snapshots=false,
)

res_drop_turb = run_core_formation_sim(
    64;
    percolation_active=true,
    settling_active=true,
    droplet_mode=:weber_turbulent,
    segregation_heating=true,
    t_end_yr=3.5e6,
    dt_yr=1000.0,
    save_snapshots=false,
)

println("4. Running Energetics Sweep: Segregation Dissipation Heating ON vs OFF...")
res_no_heating = run_core_formation_sim(
    128;
    percolation_active=true,
    settling_active=true,
    droplet_mode=:capillary_mean,
    segregation_heating=false,
    t_end_yr=3.5e6,
    dt_yr=500.0,
    save_snapshots=false,
)

# -----------------------------------------------------------------------------
# Save Structured Benchmark Dataset
# -----------------------------------------------------------------------------
output_dir = normpath(joinpath(@__DIR__, "..", "output_files"))
mkpath(output_dir)
output_json = joinpath(output_dir, "core_formation_benchmark_data.json")

benchmark_data = Dict{String,Any}(
    "reference" => res_ref,
    "perc_only" => res_perc_only,
    "settle_only" => res_settle_only,
    "drop_fixed" => res_drop_fixed,
    "drop_turb" => res_drop_turb,
    "no_heating" => res_no_heating,
)

println("Writing benchmark dataset to $output_json...")
open(output_json, "w") do io
    return JSON.print(io, benchmark_data)
end

println("=== Core Formation Benchmark Generation Complete ===")

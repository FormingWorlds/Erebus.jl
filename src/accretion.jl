"""
Planetesimal accretion engine: orbital mechanics, pebble accretion, gravitational focusing,
impact heating thermodynamics, and marker boundary expansion.
"""

using DocStringExtensions
using LinearAlgebra

# Gravitational constant [m^3 / (kg s^2)]
const G_GRAV = 6.67430e-11

# Boltzmann constant [J / K]
const K_BOLTZMANN = 1.380649e-23

# Proton mass [kg]
const M_PROTON = 1.67262192e-27

# Astronomical unit [m]
const AU_METERS = 1.495978707e11

# Solar mass [kg]
const M_SUN_KG = 1.98847e30

# Seconds per year [s]
const SEC_PER_YEAR = 3.15576e7

"""
Compute Keplerian orbital angular frequency at semi-major axis `a` around star of mass `M_star`.

$(SIGNATURES)

# Arguments
- `a_m`: Semi-major axis [m]
- `M_star_kg`: Stellar mass [kg] (default: 1 Solar mass)

# Returns
- `Omega_K`: Keplerian angular frequency [rad/s]
"""
function compute_keplerian_frequency(a_m::Real, M_star_kg::Real=M_SUN_KG)::Float64
    a_val = Float64(a_m)
    M_val = Float64(M_star_kg)
    if !isfinite(a_val) || a_val <= 0.0
        throw(DomainError(a_val, "Semi-major axis must be positive and finite"))
    end
    if !isfinite(M_val) || M_val <= 0.0
        throw(DomainError(M_val, "Stellar mass must be positive and finite"))
    end
    return sqrt(G_GRAV * M_val / (a_val^3))
end

"""
Compute Keplerian orbital velocity at semi-major axis `a` around star of mass `M_star`.

$(SIGNATURES)

# Arguments
- `a_m`: Semi-major axis [m]
- `M_star_kg`: Stellar mass [kg] (default: 1 Solar mass)

# Returns
- `v_K`: Keplerian orbital velocity [m/s]
"""
function compute_keplerian_velocity(a_m::Real, M_star_kg::Real=M_SUN_KG)::Float64
    a_val = Float64(a_m)
    M_val = Float64(M_star_kg)
    if !isfinite(a_val) || a_val <= 0.0
        throw(DomainError(a_val, "Semi-major axis must be positive and finite"))
    end
    if !isfinite(M_val) || M_val <= 0.0
        throw(DomainError(M_val, "Stellar mass must be positive and finite"))
    end
    return sqrt(G_GRAV * M_val / a_val)
end

"""
Compute sound speed of protoplanetary disk gas.

$(SIGNATURES)

# Arguments
- `T_gas`: Gas temperature [K]

# Keyword Arguments
- `gamma`: Adiabatic index (default: 1.0 for isothermal disk gas)
- `mu_gas`: Mean molecular weight in proton masses (default: 2.34 for H2/He mixture)

# Returns
- `c_s`: Gas sound speed [m/s]
"""
function compute_sound_speed(T_gas::Real; gamma::Real=1.0, mu_gas::Real=2.34)::Float64
    T_val = Float64(T_gas)
    gam_val = Float64(gamma)
    mu_val = Float64(mu_gas)
    if !isfinite(T_val) || T_val <= 0.0
        throw(DomainError(T_val, "Gas temperature must be positive and finite"))
    end
    if !isfinite(gam_val) || gam_val <= 0.0
        throw(DomainError(gam_val, "Adiabatic index must be positive and finite"))
    end
    if !isfinite(mu_val) || mu_val <= 0.0
        throw(DomainError(mu_val, "Mean molecular weight must be positive and finite"))
    end
    return sqrt(gam_val * K_BOLTZMANN * T_val / (mu_val * M_PROTON))
end

"""
Compute disk gas vertical scale height.

$(SIGNATURES)

# Arguments
- `c_s`: Gas sound speed [m/s]
- `Omega_K`: Keplerian frequency [rad/s]

# Returns
- `H_gas`: Gas scale height [m]
"""
function compute_gas_scale_height(c_s::Real, Omega_K::Real)::Float64
    cs_val = Float64(c_s)
    om_val = Float64(Omega_K)
    if !isfinite(cs_val) || cs_val <= 0.0
        throw(DomainError(cs_val, "Sound speed must be positive and finite"))
    end
    if !isfinite(om_val) || om_val <= 0.0
        throw(DomainError(om_val, "Keplerian frequency must be positive and finite"))
    end
    return cs_val / om_val
end

"""
Compute settled pebble midplane layer vertical scale height.

$(SIGNATURES)

# Arguments
- `H_gas`: Gas scale height [m]
- `St`: Aerodynamic Stokes number [-]
- `alpha_turb`: Dimensionless Shakura-Sunyaev turbulence parameter [-]

# Returns
- `H_peb`: Pebble layer scale height [m]
"""
function compute_pebble_scale_height(H_gas::Real, St::Real, alpha_turb::Real)::Float64
    H_val = Float64(H_gas)
    St_val = Float64(St)
    al_val = Float64(alpha_turb)
    if !isfinite(H_val) || H_val <= 0.0
        throw(DomainError(H_val, "Gas scale height must be positive and finite"))
    end
    if !isfinite(St_val) || St_val <= 0.0
        throw(DomainError(St_val, "Stokes number must be positive and finite"))
    end
    if !isfinite(al_val) || al_val <= 0.0
        throw(DomainError(al_val, "Turbulent alpha must be positive and finite"))
    end
    return H_val * sqrt(al_val / (al_val + St_val))
end

"""
Compute gravitational Bondi radius for gas or pebble capture.

$(SIGNATURES)

# Arguments
- `M`: Planetesimal mass [kg]
- `c_s`: Gas sound speed [m/s]

# Returns
- `R_B`: Bondi radius [m]
"""
function compute_bondi_radius(M::Real, c_s::Real)::Float64
    M_val = Float64(M)
    cs_val = Float64(c_s)
    if !isfinite(M_val) || M_val < 0.0
        throw(DomainError(M_val, "Mass must be non-negative and finite"))
    end
    if !isfinite(cs_val) || cs_val <= 0.0
        throw(DomainError(cs_val, "Sound speed must be positive and finite"))
    end
    return G_GRAV * M_val / (cs_val^2)
end

"""
Compute gravitational Hill radius of planetesimal at semi-major axis `a`.

$(SIGNATURES)

# Arguments
- `M`: Planetesimal mass [kg]
- `a`: Semi-major axis [m]
- `M_star`: Stellar mass [kg] (default: M_SUN_KG)

# Returns
- `R_H`: Hill radius [m]
"""
function compute_hill_radius(M::Real, a::Real, M_star::Real=M_SUN_KG)::Float64
    M_val = Float64(M)
    a_val = Float64(a)
    Ms_val = Float64(M_star)

    if !isfinite(M_val) || M_val < 0.0
        throw(DomainError(M_val, "Mass must be non-negative and finite"))
    end
    if !isfinite(a_val) || a_val <= 0.0
        throw(DomainError(a_val, "Semi-major axis must be positive and finite"))
    end
    if !isfinite(Ms_val) || Ms_val <= 0.0
        throw(DomainError(Ms_val, "Stellar mass must be positive and finite"))
    end
    return a_val * cbrt(M_val / (3.0 * Ms_val))
end

"""
Compute radial power-law pebble surface density across the protoplanetary disk.

$(SIGNATURES)

# Arguments
- `a_au`: Orbital semi-major axis in Astronomical Units [AU]

# Keyword Arguments
- `Sigma_peb_0`: Reference pebble surface density at 1 AU [kg/m^2] (default: 50.0)
- `p_peb`: Radial power-law index (default: 1.0)

# Returns
- `Sigma_peb`: Local pebble surface density [kg/m^2]
"""
function compute_pebble_surface_density(
    a_au::Real; Sigma_peb_0::Real=50.0, p_peb::Real=1.0
)::Float64
    a_val = Float64(a_au)
    sig0_val = Float64(Sigma_peb_0)
    p_val = Float64(p_peb)
    if !isfinite(a_val) || a_val <= 0.0
        throw(DomainError(a_val, "Orbital distance must be positive and finite"))
    end
    if !isfinite(sig0_val) || sig0_val < 0.0
        throw(
            DomainError(
                sig0_val, "Pebble reference density must be non-negative and finite"
            ),
        )
    end
    if !isfinite(p_val)
        throw(DomainError(p_val, "Power-law index must be finite"))
    end
    return sig0_val * (a_val^(-p_val))
end

"""
Compute pebble accretion rate [kg/s] onto planetesimal in Bondi, Hill, or automated transition regime.

$(SIGNATURES)

# Arguments
- `M`: Planetesimal mass [kg]
- `M_star`: Stellar mass [kg]
- `a`: Semi-major axis [m]
- `Sigma_peb`: Local pebble surface density [kg/m^2]
- `St`: Stokes number [-]
- `c_s`: Gas sound speed [m/s]
- `alpha_turb`: Turbulent viscosity parameter [-]

# Keyword Arguments
- `c_bondi`: Calibration prefactor for Bondi regime (default: 1.0)
- `c_hill`: Calibration prefactor for Hill regime (default: 1.0)
- `regime`: Accretion mode (`:auto`, `:pebble_bondi`, or `:pebble_hill`)

# Returns
- `dM_dt`: Pebble accretion mass rate [kg/s]
"""
function compute_pebble_accretion_rate(
    M::Real,
    M_star::Real,
    a::Real,
    Sigma_peb::Real,
    St::Real,
    c_s::Real,
    alpha_turb::Real;
    c_bondi::Real=1.0,
    c_hill::Real=1.0,
    regime::Symbol=:auto,
)::Float64
    M_val = Float64(M)
    if M_val <= 0.0 || Sigma_peb <= 0.0
        return 0.0
    end

    Omega_K = compute_keplerian_frequency(a, M_star)
    v_K = compute_keplerian_velocity(a, M_star)
    H_gas = compute_gas_scale_height(c_s, Omega_K)
    H_peb = compute_pebble_scale_height(H_gas, St, alpha_turb)
    rho_peb = Float64(Sigma_peb) / (sqrt(2.0 * pi) * H_peb)

    # Sub-Keplerian gas headwind velocity: eta * v_K with eta ~ (c_s / v_K)^2
    eta_disk = 1.5 * (Float64(c_s) / v_K)^2
    v_rel = eta_disk * v_K

    # Hill radius and orbital shearing velocity
    R_H = compute_hill_radius(M_val, a, M_star)
    v_H = Omega_K * R_H

    # Bondi regime: capture radius bounded by Hill radius (Lambrechts & Johansen 2012 eq. 6)
    R_acc_B = min(
        R_H, 2.0 * sqrt((Float64(St) / Omega_K) * G_GRAV * M_val / max(v_rel, 1.0e-6))
    )
    dM_B_2D = 2.0 * Float64(c_bondi) * R_acc_B * Float64(Sigma_peb) * v_rel
    dM_B_3D = pi * Float64(c_bondi) * (R_acc_B^2) * rho_peb * v_rel
    dM_Bondi = min(dM_B_2D, dM_B_3D)

    # Hill regime: R_acc,H = R_H * St^(1/3)
    R_acc_H = R_H * cbrt(Float64(St))
    dM_H_2D = 2.0 * Float64(c_hill) * R_acc_H * Float64(Sigma_peb) * v_H
    dM_H_3D = pi * Float64(c_hill) * (R_acc_H^2) * rho_peb * v_H
    dM_Hill = min(dM_H_2D, dM_H_3D)

    if regime === :pebble_bondi
        return max(0.0, dM_Bondi)
    elseif regime === :pebble_hill
        return max(0.0, dM_Hill)
    elseif regime === :auto || regime === :pebble_auto
        # Bondi regime for M < M_trans, Hill regime for M >= M_trans (Lambrechts & Johansen 2012)
        M_trans = sqrt(1.0 / 3.0) * (v_rel^3) / (G_GRAV * Omega_K)
        return M_val < M_trans ? max(0.0, dM_Bondi) : max(0.0, dM_Hill)
    else
        throw(
            ArgumentError(
                "Unrecognized pebble accretion regime: :$regime. Must be :auto, :pebble_auto, :pebble_bondi, or :pebble_hill",
            ),
        )
    end
end

"""
Compute Safronov planetesimal collision accretion rate with gravitational focusing.

$(SIGNATURES)

# Arguments
- `M`: Planetesimal mass [kg]
- `R`: Planetesimal radius [m]
- `Sigma_pl`: Planetesimal swarm surface density [kg/m^2]
- `sigma_v`: Planetesimal velocity dispersion [m/s]
- `Omega_K`: Keplerian orbital frequency [rad/s]

# Returns
- `dM_dt`: Accretion mass rate [kg/s]
"""
function compute_safronov_accretion_rate(
    M::Real, R::Real, Sigma_pl::Real, sigma_v::Real, Omega_K::Real
)::Float64
    M_val = Float64(M)
    R_val = Float64(R)
    sig_val = Float64(Sigma_pl)
    sv_val = Float64(sigma_v)
    om_val = Float64(Omega_K)

    if M_val <= 0.0 || R_val <= 0.0 || sig_val <= 0.0 || om_val <= 0.0
        return 0.0
    end
    if !isfinite(sv_val) || sv_val <= 0.0
        throw(DomainError(sv_val, "Velocity dispersion must be positive and finite"))
    end

    # Gravitational focusing factor: F_g = 1 + 2 * Theta with Theta = G M / (R sigma_v^2)
    v_esc_sq = 2.0 * G_GRAV * M_val / R_val
    Theta = v_esc_sq / (2.0 * (sv_val^2))
    F_g = 1.0 + 2.0 * Theta

    # Geometric cross-section enhancement
    return pi * (R_val^2) * sig_val * om_val * F_g
end

"""
Compute specific impact energy and resulting temperature rise in newly accreted shell.

$(SIGNATURES)

# Arguments
- `M`: Planetesimal mass [kg]
- `R`: Planetesimal radius [m]

# Keyword Arguments
- `h_impact`: Retention fraction of impact kinetic energy as heat (default: 0.5)
- `c_p`: Rock specific heat capacity [J/(kg K)] (default: 1000.0)
- `v_inf`: Impactor approach velocity at infinity [m/s] (default: 0.0)

# Returns
- `(u_acc, Delta_T)`: Specific kinetic energy [J/kg] and temperature rise [K]
"""
function compute_impact_heating(
    M::Real, R::Real; h_impact::Real=0.5, c_p::Real=1000.0, v_inf::Real=0.0
)::Tuple{Float64,Float64}
    M_val = Float64(M)
    R_val = Float64(R)
    cp_val = Float64(c_p)
    h_val = Float64(h_impact)
    v_val = Float64(v_inf)

    if !isfinite(M_val) || M_val < 0.0
        throw(DomainError(M_val, "Mass must be non-negative and finite"))
    end
    if !isfinite(R_val) || R_val <= 0.0
        throw(DomainError(R_val, "Radius must be positive and finite"))
    end
    if !isfinite(cp_val) || cp_val <= 0.0
        throw(DomainError(cp_val, "Heat capacity must be positive and finite"))
    end
    if !isfinite(h_val) || !(0.0 <= h_val <= 1.0)
        throw(DomainError(h_val, "Impact retention fraction must be in [0, 1]"))
    end
    if !isfinite(v_val) || v_val < 0.0
        throw(DomainError(v_val, "Velocity at infinity must be non-negative and finite"))
    end

    # Specific kinetic energy of accretion: u_acc = G M / R + 0.5 v_inf^2 [J/kg]
    u_acc = (G_GRAV * M_val / R_val) + 0.5 * (v_val^2)
    delta_T = Float64(h_val * u_acc / cp_val)
    return (u_acc, delta_T)
end

"""
Compute exact 3D spherical radius increment ΔR from accreted mass ΔM and bulk density ρ.

$(SIGNATURES)

# Arguments
- `R`: Current planetary radius [m]
- `dM`: Accreted mass increment [kg]
- `rho_bulk`: Bulk density of accreted shell [kg/m^3]

# Returns
- `Delta_R`: Radius increment [m]
"""
function compute_radius_increment(R::Real, dM::Real, rho_bulk::Real)::Float64
    R_val = Float64(R)
    dM_val = Float64(dM)
    rho_val = Float64(rho_bulk)

    if !isfinite(R_val) || R_val <= 0.0
        throw(DomainError(R_val, "Radius must be positive and finite"))
    end
    if !isfinite(dM_val) || dM_val < 0.0
        throw(DomainError(dM_val, "Mass increment must be non-negative and finite"))
    end
    if !isfinite(rho_val) || rho_val <= 0.0
        throw(DomainError(rho_val, "Bulk density must be positive and finite"))
    end

    if dM_val == 0.0
        return 0.0
    end

    # Numerically stable spherical shell mapping without catastrophic cancellation:
    # R_new = R * (1 + x)^(1/3), where x = 3 dM / (4 pi rho R^3)
    x = (3.0 * dM_val) / (4.0 * pi * rho_val * (R_val^3))
    return R_val * expm1(log1p(x) / 3.0)
end

"""
Evaluate volatile water fraction from protoplanetary disk water snowline.

$(SIGNATURES)

# Arguments
- `T_disk`: Local disk temperature [K]

# Keyword Arguments
- `T_snowline_cond`: Water condensation temperature threshold [K] (default: 160.0)
- `XW_wet`: Hydrated silicate fraction outside snowline (default: 0.40)
- `XW_dry`: Anhydrous silicate fraction inside snowline (default: 0.0)
- `H2O_wet_wtpct`: Bulk water wt% outside snowline (default: 10.0)
- `H2O_dry_wtpct`: Bulk water wt% inside snowline (default: 0.1)

# Returns
- `(XWsolid, XH2O_wtpct)`: Hydrated rock fraction and bulk water content wt%
"""
function evaluate_snowline_water_content(
    T_disk::Real;
    T_snowline_cond::Real=160.0,
    XW_wet::Real=0.40,
    XW_dry::Real=0.0,
    H2O_wet_wtpct::Real=10.0,
    H2O_dry_wtpct::Real=0.1,
)::Tuple{Float64,Float64}
    T_val = Float64(T_disk)
    T_snow = Float64(T_snowline_cond)
    if T_val <= T_snow
        # Cold region outside snowline: volatile-rich ice/silicate mixture
        return (Float64(XW_wet), Float64(H2O_wet_wtpct))
    else
        # Warm region inside snowline: devolatilized anhydrous silicate rock
        return (Float64(XW_dry), Float64(H2O_dry_wtpct))
    end
end

"""
Advance planetesimal accretion boundary and convert sticky-air markers inside the new radius.

$(SIGNATURES)

# Arguments
- `R_current`: Current planetary radius [m]
- `delta_R`: Radius expansion increment [m]
- `xm`: Marker x-coordinates [m]
- `ym`: Marker y-coordinates [m]
- `tm`: Marker material type indices (1: core, 2: crust/rock, 3: sticky air)
- `tkm`: Marker temperature array [K]
- `phim`: Marker porosity array [-]
- `XWsolidm0`: Marker hydrated solid fraction array [-]
- `Xfe_bulk`: Marker bulk metal fraction array [-] (or nothing)
- `Xfem`: Marker molten metal fraction array [-] (or nothing)

# Keyword Arguments
- `xcenter`: Horizontal center coordinate of planetesimal [m]
- `ycenter`: Vertical center coordinate of planetesimal [m]
- `T_accreted`: Initial temperature of accreted shell [K]
- `phi_accreted`: Initial porosity of accreted shell [-]
- `XWsolid_accreted`: Hydrated solid fraction of accreted shell [-]
- `Xfe_accreted`: Bulk metal fraction of accreted shell [-]
- `t_accreted`: Optional array to record accretion epoch timestamps [s]
- `current_time`: Current simulation time [s]
- `XH2Om`: Optional water tracer array
- `XCm`: Optional carbon tracer array
- `XNm`: Optional nitrogen tracer array
- `XSm`: Optional sulfur tracer array
- `XH2O_accreted`: Water wt% for accreted shell
- `XC_accreted`: Carbon ppm for accreted shell
- `XN_accreted`: Nitrogen ppm for accreted shell
- `XS_accreted`: Sulfur ppm for accreted shell

# Returns
- `n_converted`: Number of sticky-air markers converted to rock
"""
function advance_accretion_boundary!(
    R_current::Real,
    delta_R::Real,
    xm::AbstractVector{<:Real},
    ym::AbstractVector{<:Real},
    tm::AbstractVector{<:Integer},
    tkm::AbstractVector{<:Real},
    phim::AbstractVector{<:Real},
    XWsolidm0::AbstractVector{<:Real},
    Xfe_bulk::Union{Nothing,AbstractVector{<:Real}}=nothing,
    Xfem::Union{Nothing,AbstractVector{<:Real}}=nothing;
    xcenter::Real=70000.0,
    ycenter::Real=70000.0,
    T_accreted::Real=200.0,
    phi_accreted::Real=0.35,
    XWsolid_accreted::Real=0.40,
    Xfe_accreted::Real=0.10,
    t_accreted::Union{Nothing,AbstractVector{<:Real}}=nothing,
    current_time::Real=0.0,
    XWsolidm::Union{Nothing,AbstractVector{<:Real}}=nothing,
    phinewm::Union{Nothing,AbstractVector{<:Real}}=nothing,
    XH2Om::Union{Nothing,AbstractVector{<:Real}}=nothing,
    XCm::Union{Nothing,AbstractVector{<:Real}}=nothing,
    XNm::Union{Nothing,AbstractVector{<:Real}}=nothing,
    XSm::Union{Nothing,AbstractVector{<:Real}}=nothing,
    XH2O_accreted::Real=10.0,
    XC_accreted::Real=1000.0,
    XN_accreted::Real=100.0,
    XS_accreted::Real=10000.0,
    hcnspo_props=nothing,
    disk_state=nothing,
)::Int
    R_new = Float64(R_current) + Float64(delta_R)
    R_new_sq = R_new^2
    xc = Float64(xcenter)
    yc = Float64(ycenter)

    n_converted = Threads.Atomic{Int}(0)
    marknum = length(tm)

    Threads.@threads :static for m in 1:marknum
        # Only convert sticky air / space markers (tm == 3)
        if @inbounds tm[m] == 3
            dx_m = Float64(@inbounds(xm[m])) - xc
            dy_m = Float64(@inbounds(ym[m])) - yc
            r_sq = dx_m^2 + dy_m^2
            if r_sq <= R_new_sq
                @inbounds tm[m] = 2  # Convert to solid crust/rock phase
                @inbounds tkm[m] = Float64(T_accreted)
                @inbounds phim[m] = Float64(phi_accreted)
                @inbounds XWsolidm0[m] = Float64(XWsolid_accreted)
                if phinewm !== nothing
                    @inbounds phinewm[m] = Float64(phi_accreted)
                end
                if XWsolidm !== nothing
                    @inbounds XWsolidm[m] = Float64(XWsolid_accreted)
                end

                if Xfe_bulk !== nothing
                    @inbounds Xfe_bulk[m] = Float64(Xfe_accreted)
                end
                if Xfem !== nothing
                    @inbounds Xfem[m] = 0.0  # Cold accreted metal is solid
                end

                if XH2Om !== nothing
                    @inbounds XH2Om[m] = Float64(XH2O_accreted)
                end
                if XCm !== nothing
                    @inbounds XCm[m] = Float64(XC_accreted)
                end
                if XNm !== nothing
                    @inbounds XNm[m] = Float64(XN_accreted)
                end
                if XSm !== nothing
                    @inbounds XSm[m] = Float64(XS_accreted)
                end

                if hcnspo_props !== nothing && disk_state !== nothing
                    @inbounds hcnspo_props.X_ice_H2O_m[m] = Float64(disk_state.X_ice_H2O)
                    @inbounds hcnspo_props.X_ice_NH3_m[m] = Float64(disk_state.X_ice_NH3)
                    @inbounds hcnspo_props.X_ice_CO2_m[m] = Float64(disk_state.X_ice_CO2)
                    @inbounds hcnspo_props.X_ice_CO_m[m] = Float64(disk_state.X_ice_CO)
                    @inbounds hcnspo_props.X_ice_CH4_m[m] = Float64(disk_state.X_ice_CH4)
                    @inbounds hcnspo_props.X_ice_N2_m[m] = Float64(disk_state.X_ice_N2)
                    @inbounds hcnspo_props.X_ice_H2S_m[m] = Float64(disk_state.X_ice_H2S)
                    @inbounds hcnspo_props.X_ice_PH3_m[m] = Float64(disk_state.X_ice_PH3)

                    @inbounds hcnspo_props.X_refr_C_m[m] = Float64(disk_state.f_refr_C)
                    @inbounds hcnspo_props.X_refr_S_m[m] = Float64(disk_state.f_refr_S)
                    @inbounds hcnspo_props.X_refr_N_m[m] = Float64(disk_state.f_refr_N)
                    @inbounds hcnspo_props.X_refr_P_m[m] = Float64(disk_state.f_refr_P)
                    @inbounds hcnspo_props.X_refr_H_m[m] = Float64(disk_state.f_refr_H)
                end

                if t_accreted !== nothing
                    @inbounds t_accreted[m] = Float64(current_time)
                end

                Threads.atomic_add!(n_converted, 1)
            end
        end
    end

    return n_converted[]
end

"""
Compute total accretion mass rate [kg/s] for the planetesimal according to configured mode.

$(SIGNATURES)

# Arguments
- `time_seconds`: Current simulation time [s]
- `M`: Current planetesimal mass [kg]
- `R`: Current planetesimal radius [m]
- `acc_cfg`: AccretionConfig struct
- `disk_cfg`: DiskConfig struct

# Returns
- `dM_dt`: Accretion mass rate [kg/s]
"""
function compute_accretion_rate(
    time_seconds::Real,
    M::Real,
    R::Real,
    acc_cfg::AccretionConfig,
    disk_cfg::DiskConfig=DiskConfig(),
)::Float64
    if !acc_cfg.active
        return 0.0
    end

    t_sec = Float64(time_seconds)
    t_start_sec = acc_cfg.t_start_myr * 1.0e6 * SEC_PER_YEAR
    t_end_sec = t_start_sec + acc_cfg.t_duration_myr * 1.0e6 * SEC_PER_YEAR

    if t_sec < t_start_sec || t_sec >= t_end_sec
        return 0.0
    end

    M_val = Float64(M)
    R_val = Float64(R)

    if M_val >= acc_cfg.M_target || R_val >= acc_cfg.R_target
        return 0.0
    end

    if acc_cfg.mode === :constant_rate
        return max(0.0, acc_cfg.dM_dt_constant)
    elseif acc_cfg.mode === :linear_radius
        # dM/dt = 4 pi R^2 rho_bulk dR/dt
        return max(0.0, 4.0 * pi * (R_val^2) * acc_cfg.rho_bulk * acc_cfg.dR_dt_constant)
    elseif acc_cfg.mode === :exponential
        tau_sec = acc_cfg.tau_growth_myr * 1.0e6 * SEC_PER_YEAR
        return max(0.0, M_val / tau_sec)
    elseif acc_cfg.mode === :safronov
        a_m = disk_cfg.orbital_distance_au * AU_METERS
        M_star = disk_cfg.stellar_mass_msun * M_SUN_KG
        Omega_K = compute_keplerian_frequency(a_m, M_star)
        sigma_v = acc_cfg.v_disp_kms * 1000.0
        return compute_safronov_accretion_rate(
            M_val, R_val, acc_cfg.Sigma_pl_0, sigma_v, Omega_K
        )
    elseif acc_cfg.mode in Set([:pebble_bondi, :pebble_hill, :pebble_auto])
        a_m = disk_cfg.orbital_distance_au * AU_METERS
        M_star = disk_cfg.stellar_mass_msun * M_SUN_KG
        Sigma_peb = compute_pebble_surface_density(
            disk_cfg.orbital_distance_au;
            Sigma_peb_0=acc_cfg.Sigma_peb_0,
            p_peb=acc_cfg.p_peb,
        )
        T_disk = if disk_cfg.enabled
            compute_disk_temperature(t_sec, disk_cfg)
        else
            disk_cfg.t_ambient
        end
        c_s = compute_sound_speed(T_disk)
        return compute_pebble_accretion_rate(
            M_val,
            M_star,
            a_m,
            Sigma_peb,
            acc_cfg.stokes_number,
            c_s,
            acc_cfg.alpha_turbulence;
            c_bondi=acc_cfg.c_bondi,
            c_hill=acc_cfg.c_hill,
            regime=acc_cfg.mode,
        )
    else
        return 0.0
    end
end

"""
Evaluate multi-snowline volatile condensation and refractory delivery state in the protoplanetary disk.

$(SIGNATURES)

# Parameters
- `T_disk`: Midplane disk temperature [K].
- `P_disk`: Midplane gas pressure [Pa].
- `mix_cfg`: Volatile mixture configuration struct.
- `refr_cfg`: Refractory phase configuration struct.

# Keywords
- `P_ref`: Reference midplane pressure for snowline condensation [Pa] (default: 1.0).
- `alpha_P`: Pressure sensitivity coefficient for Clausius-Clapeyron snowline shifting (default: 0.0).

# Returns
- Named tuple with condensation boolean flags, volatile ice mass fractions, and refractory delivery fractions.

# Raises
- `DomainError`: If disk temperature, pressure, P_ref, or alpha_P is negative or NaN.
"""
function evaluate_disk_volatile_condensation(
    T_disk::Real,
    P_disk::Real=1.0,
    mix_cfg::VolatileMixtureConfig=VolatileMixtureConfig(),
    refr_cfg::RefractoryConfig=RefractoryConfig();
    P_ref::Real=mix_cfg.P_ref,
    alpha_P::Real=mix_cfg.alpha_P,
)
    (isnan(T_disk) || T_disk < 0.0) &&
        throw(DomainError(T_disk, "Disk temperature must be non-negative"))
    (isnan(P_disk) || P_disk < 0.0) &&
        throw(DomainError(P_disk, "Disk pressure must be non-negative"))
    (isnan(P_ref) || P_ref <= 0.0) && throw(DomainError(P_ref, "P_ref must be positive"))
    (isnan(alpha_P) || alpha_P < 0.0) &&
        throw(DomainError(alpha_P, "alpha_P must be non-negative"))

    T = Float64(T_disk)
    p_factor = max(
        0.01, 1.0 + Float64(alpha_P) * log(max(Float64(P_disk), 1.0e-8) / Float64(P_ref))
    )

    condensed_H2O = mix_cfg.active && (T <= mix_cfg.T_cond_H2O * p_factor)
    condensed_NH3 = mix_cfg.active && (T <= mix_cfg.T_cond_NH3 * p_factor)
    condensed_CO2 = mix_cfg.active && (T <= mix_cfg.T_cond_CO2 * p_factor)
    condensed_H2S = mix_cfg.active && (T <= mix_cfg.T_cond_H2S * p_factor)
    condensed_CH4 = mix_cfg.active && (T <= mix_cfg.T_cond_CH4 * p_factor)
    condensed_CO = mix_cfg.active && (T <= mix_cfg.T_cond_CO * p_factor)
    condensed_N2 = mix_cfg.active && (T <= mix_cfg.T_cond_N2 * p_factor)
    condensed_PH3 = mix_cfg.active && (T <= mix_cfg.T_cond_PH3 * p_factor)

    X_ice_H2O = condensed_H2O ? Float64(mix_cfg.X_ice_H2O) : 0.0
    X_ice_NH3 = condensed_NH3 ? Float64(mix_cfg.X_ice_NH3) : 0.0
    X_ice_CO2 = condensed_CO2 ? Float64(mix_cfg.X_ice_CO2) : 0.0
    X_ice_H2S = condensed_H2S ? Float64(mix_cfg.X_ice_H2S) : 0.0
    X_ice_CH4 = condensed_CH4 ? Float64(mix_cfg.X_ice_CH4) : 0.0
    X_ice_CO = condensed_CO ? Float64(mix_cfg.X_ice_CO) : 0.0
    X_ice_N2 = condensed_N2 ? Float64(mix_cfg.X_ice_N2) : 0.0
    X_ice_PH3 = condensed_PH3 ? Float64(mix_cfg.X_ice_PH3) : 0.0

    f_refr_C = refr_cfg.active ? Float64(refr_cfg.f_refr_C) : 0.0
    f_refr_S = refr_cfg.active ? Float64(refr_cfg.f_refr_S) : 0.0
    f_refr_N = refr_cfg.active ? Float64(refr_cfg.f_refr_N) : 0.0
    f_refr_P = refr_cfg.active ? Float64(refr_cfg.f_refr_P) : 0.0
    f_refr_H = refr_cfg.active ? Float64(refr_cfg.f_refr_H) : 0.0

    return (;
        condensed_H2O,
        condensed_NH3,
        condensed_CO2,
        condensed_H2S,
        condensed_CH4,
        condensed_CO,
        condensed_N2,
        condensed_PH3,
        X_ice_H2O,
        X_ice_NH3,
        X_ice_CO2,
        X_ice_H2S,
        X_ice_CH4,
        X_ice_CO,
        X_ice_N2,
        X_ice_PH3,
        f_refr_C,
        f_refr_S,
        f_refr_N,
        f_refr_P,
        f_refr_H,
    )
end

"""
Advance planetesimal accretion boundary with multi-component volatile ices and refractory phases.

$(SIGNATURES)

# Parameters
- `R_current`: Current planetesimal radius [m].
- `delta_R`: Incremental radial growth [m].
- `xm`: Marker x-coordinates [m].
- `ym`: Marker y-coordinates [m].
- `tm`: Marker material type array.
- `tkm`: Marker temperature array [K].
- `phim`: Marker porosity array.
- `XWsolidm0`: Reference solid water content array.
- `hcnspo_props`: Marker HCNSPO properties struct or named tuple.
- `disk_state`: Disk condensation state from `evaluate_disk_volatile_condensation`.

# Keywords
- `xcenter`: Domain x-center coordinate [m] (default: 70000.0).
- `ycenter`: Domain y-center coordinate [m] (default: 70000.0).
- `T_accreted`: Default accreted temperature [K] (default: 150.0).
- `phi_accreted`: Default accreted porosity (default: 0.35).
- `current_time`: Current simulation time [s] (default: 0.0).
- `t_accreted`: Optional marker accretion epoch array.

# Returns
- `Int`: Count of converted markers from sticky air (tm == 3) to rock (tm == 2).
"""
function advance_accretion_boundary_hcnspo!(
    R_current::Real,
    delta_R::Real,
    xm::AbstractVector{<:Real},
    ym::AbstractVector{<:Real},
    tm::AbstractVector{<:Integer},
    tkm::AbstractVector{<:Real},
    phim::AbstractVector{<:Real},
    XWsolidm0::AbstractVector{<:Real},
    hcnspo_props,
    disk_state;
    xcenter::Real=70000.0,
    ycenter::Real=70000.0,
    T_accreted::Real=150.0,
    phi_accreted::Real=0.35,
    current_time::Real=0.0,
    t_accreted::Union{Nothing,AbstractVector{<:Real}}=nothing,
)::Int
    return advance_accretion_boundary!(
        R_current,
        delta_R,
        xm,
        ym,
        tm,
        tkm,
        phim,
        XWsolidm0,
        nothing,
        nothing;
        xcenter=xcenter,
        ycenter=ycenter,
        T_accreted=T_accreted,
        phi_accreted=phi_accreted,
        XWsolid_accreted=disk_state.X_ice_H2O,
        current_time=current_time,
        t_accreted=t_accreted,
        hcnspo_props=hcnspo_props,
        disk_state=disk_state,
    )
end

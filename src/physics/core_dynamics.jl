# -----------------------------------------------------------------------------
# Iron Core Formation Physics (Percolation & Gravitational Settling)
# -----------------------------------------------------------------------------

"""
Compute liquid metal permeability using a modified Kozeny-Carman formulation.

$(SIGNATURES)

References:
- McKenzie (1984), J. Petrol., 25(3), 713-765.
- Gerya (2019), Introduction to Numerical Geodynamic Modelling, Cambridge Univ. Press.

# Arguments
- `phi_m`: Metal volume fraction [-]

# Keyword Arguments
- `k_metal_ref`: Reference permeability at reference porosity [m^2] (default 1.0e-9)
- `phi_crit_perc`: Percolation connectivity threshold [-] (default 0.05)
- `perm_exponent`: Porosity power exponent [-] (default 3.0)
- `phi0`: Reference porosity scale [-] (default 0.1)

# Returns
- Permeability of the interconnected liquid metal network [m^2].

# Raises
- `DomainError`: If `phi_m` is not in [0, 1), or parameters are non-positive/unphysical.
"""
function metal_permeability(
    phi_m::Real;
    k_metal_ref::Real=1.0e-9,
    phi_crit_perc::Real=0.05,
    perm_exponent::Real=3.0,
    phi0::Real=0.1,
)
    if !(0.0 <= phi_m < 1.0) || !isfinite(phi_m)
        throw(DomainError(phi_m, "Metal volume fraction must be in [0, 1) and finite"))
    end
    if k_metal_ref < 0.0 || !isfinite(k_metal_ref)
        throw(
            DomainError(
                k_metal_ref, "Reference permeability must be non-negative and finite"
            ),
        )
    end
    if !(0.0 <= phi_crit_perc < 1.0) || !isfinite(phi_crit_perc)
        throw(
            DomainError(phi_crit_perc, "Percolation threshold must be in [0, 1) and finite")
        )
    end
    if perm_exponent < 0.0 || !isfinite(perm_exponent)
        throw(
            DomainError(
                perm_exponent, "Permeability exponent must be non-negative and finite"
            ),
        )
    end
    if !(0.0 < phi0 < 1.0) || !isfinite(phi0)
        throw(DomainError(phi0, "Reference porosity must be in (0, 1) and finite"))
    end

    if phi_m <= phi_crit_perc
        return 0.0
    end

    phi_mob = phi_m - phi_crit_perc
    return k_metal_ref *
           (phi_mob / phi0)^perm_exponent *
           ((1.0 - phi_mob) / (1.0 - phi0))^(-2.0)
end

"""
Compute terminal Stokes settling velocity for liquid metal droplets in magma.

$(SIGNATURES)

References:
- Rubie et al. (2003), Earth Planet. Sci. Lett., 205(3-4), 239-255.
- Hadamard (1911), C. R. Acad. Sci., 152, 1735-1738.
- Rybczynski (1911), Bull. Acad. Sci. Cracovie, A, 40-46.

# Arguments
- `r_drop`: Droplet radius [m]
- `drho`: Density contrast between metal and magma [kg/m^3]
- `g_acc`: Gravitational acceleration magnitude [m/s^2]
- `eta_susp`: Dynamic viscosity of the silicate suspension [Pa s]

# Keyword Arguments
- `hadamard_rybczynski`: Apply fluid-droplet velocity correction (1.5x) (default false)

# Returns
- Terminal settling velocity [m/s].

# Raises
- `DomainError`: If droplet radius is negative, viscosity is non-positive, or inputs non-finite.
"""
function stokes_settling_velocity(
    r_drop::Real,
    drho::Real,
    g_acc::Real,
    eta_susp::Real;
    hadamard_rybczynski::Bool=false,
    eta_metal::Union{Nothing,Real}=nothing,
)
    if r_drop < 0.0 || !isfinite(r_drop)
        throw(DomainError(r_drop, "Droplet radius must be non-negative and finite"))
    end
    if eta_susp <= 0.0 || !isfinite(eta_susp)
        throw(DomainError(eta_susp, "Suspension viscosity must be positive and finite"))
    end
    if !isfinite(drho) || !isfinite(g_acc)
        throw(DomainError((drho, g_acc), "Density contrast and gravity must be finite"))
    end
    if eta_metal !== nothing && (!isfinite(eta_metal) || eta_metal <= 0.0)
        throw(DomainError(eta_metal, "Metal viscosity must be positive and finite"))
    end

    v_st = (2.0 / 9.0) * drho * g_acc * (r_drop^2) / eta_susp
    if hadamard_rybczynski
        f_hr = if eta_metal !== nothing
            (3.0 * eta_susp + 3.0 * eta_metal) / (2.0 * eta_susp + 3.0 * eta_metal)
        else
            1.5
        end
        v_st *= f_hr
    end
    return v_st
end

"""
Compute hydrodynamic equilibrium droplet diameter from critical Weber number.

$(SIGNATURES)

References:
- Rubie et al. (2003), Earth Planet. Sci. Lett., 205(3-4), 239-255.
- Samuel (2012), Earth Planet. Sci. Lett., 313-314, 105-114.
- Wacheul & Le Bars (2018), J. Fluid Mech., 846, 5-41.

# Arguments
- `rho_silicate`: Magma ocean liquid silicate density [kg/m^3]
- `v_rel`: Relative settling or convective velocity [m/s]
- `sigma`: Interfacial tension between liquid metal and silicate [N/m]

# Keyword Arguments
- `We_crit`: Critical Weber number for hydrodynamic breakup [-] (default 10.0)

# Returns
- Equilibrium stable droplet diameter [m].

# Raises
- `DomainError`: If densities, velocity, surface tension, or We_crit are non-positive.
"""
function weber_equilibrium_diameter(
    rho_silicate::Real, v_rel::Real, sigma::Real; We_crit::Real=10.0
)
    if rho_silicate <= 0.0 || !isfinite(rho_silicate)
        throw(DomainError(rho_silicate, "Silicate density must be positive and finite"))
    end
    if v_rel <= 0.0 || !isfinite(v_rel)
        throw(DomainError(v_rel, "Relative velocity must be positive and finite"))
    end
    if sigma <= 0.0 || !isfinite(sigma)
        throw(DomainError(sigma, "Interfacial surface tension must be positive and finite"))
    end
    if We_crit <= 0.0 || !isfinite(We_crit)
        throw(DomainError(We_crit, "Critical Weber number must be positive and finite"))
    end

    return We_crit * sigma / (rho_silicate * (v_rel^2))
end

"""
Compute Richardson-Zaki hindered settling reduction factor.

$(SIGNATURES)

References:
- Richardson & Zaki (1954), Trans. Inst. Chem. Eng., 32, 35-53.

# Arguments
- `phi_m`: Dispersed metal volume fraction [-]

# Keyword Arguments
- `hindered_exponent`: Richardson-Zaki empirical power exponent [-] (default 4.5)
- `phi_pack`: Maximum packing volume fraction [-] (default 0.65)

# Returns
- Velocity hindrance factor in [0, 1].

# Raises
- `DomainError`: If `phi_m` is not in [0, 1], exponent negative, or `phi_pack` out of range.
"""
function richardson_zaki_hindrance(
    phi_m::Real; hindered_exponent::Real=4.5, phi_pack::Real=0.65
)
    if !(0.0 <= phi_m <= 1.0) || !isfinite(phi_m)
        throw(DomainError(phi_m, "Metal volume fraction must be in [0, 1] and finite"))
    end
    if hindered_exponent < 0.0 || !isfinite(hindered_exponent)
        throw(
            DomainError(
                hindered_exponent, "Hindered exponent must be non-negative and finite"
            ),
        )
    end
    if !(0.0 < phi_pack <= 1.0) || !isfinite(phi_pack)
        throw(DomainError(phi_pack, "Packing fraction must be in (0, 1] and finite"))
    end

    if phi_m >= phi_pack
        return 0.0
    end

    xi = clamp(1.0 - phi_m / phi_pack, 0.0, 1.0)
    return xi^hindered_exponent
end

"""
Compute eutectic liquid metal melt fraction as a function of temperature.

$(SIGNATURES)

References:
- Kubaschewski (1982), Iron-Binary Phase Diagrams, Springer.
- Neumann et al. (2012), Astron. Astrophys., 543, A141.

# Arguments
- `T`: Local temperature [K]

# Keyword Arguments
- `T_eutectic`: Eutectic melting temperature of Fe-FeS [K] (default 1213.0)
- `dT_metal`: Temperature interval between eutectic and complete melting [K] (default 50.0)

# Returns
- Liquid metal fraction in [0, 1].

# Raises
- `DomainError`: If temperatures are non-positive or non-finite.
"""
function compute_metal_melt_fraction(T::Real; T_eutectic::Real=1213.0, dT_metal::Real=50.0)
    if !isfinite(T) || T < 0.0
        throw(DomainError(T, "Temperature must be non-negative and finite"))
    end
    if !isfinite(T_eutectic) || T_eutectic <= 0.0
        throw(DomainError(T_eutectic, "Eutectic temperature must be positive and finite"))
    end
    if !isfinite(dT_metal) || dT_metal <= 0.0
        throw(DomainError(dT_metal, "Melting interval must be positive and finite"))
    end

    if T <= T_eutectic
        return 0.0
    elseif T >= T_eutectic + dT_metal
        return 1.0
    else
        return (T - T_eutectic) / dT_metal
    end
end

@inline function _cubic_hermite_weight(x::Real, x0::Real, x1::Real)::Float64
    xf, x0f, x1f = Float64(x), Float64(x0), Float64(x1)
    xf <= x0f && return 0.0
    xf >= x1f && return 1.0
    xi = (xf - x0f) / (x1f - x0f)
    return xi * xi * (3.0 - 2.0 * xi)
end

@inline function _percolation_interstitial_velocity(
    phi_m::Float64,
    F_m::Float64,
    drho::Float64,
    g_acc::Float64,
    phi_trap::Float64;
    percolation_active::Bool,
    F_perc_end::Float64,
    phi_crit_perc::Float64,
    k_metal_ref::Float64,
    eta_metal::Float64,
    perm_exponent::Float64,
    phi0::Float64,
)::Float64
    if !(
        percolation_active &&
        F_m < F_perc_end &&
        phi_m > phi_crit_perc &&
        phi_m > phi_trap &&
        drho > 0.0 &&
        g_acc > 0.0
    )
        return 0.0
    end
    phi_perm = min(phi_m, 1.0 - 1.0e-7)
    km = metal_permeability(
        phi_perm;
        k_metal_ref=k_metal_ref,
        phi_crit_perc=phi_crit_perc,
        perm_exponent=perm_exponent,
        phi0=phi0,
    )
    return (km / (phi_perm * eta_metal)) * drho * g_acc
end

"""
Compute unified segregation velocity magnitude across percolation and settling regimes.

$(SIGNATURES)

References:
- Solomatov (2000, 2007), Treatise on Geophysics.
- Ricard et al. (2009), J. Geophys. Res., 114, B07404.

# Arguments
- `phi_m`: Metal volume fraction [-]
- `F_m`: Silicate melt fraction [-]
- `drho`: Metal-silicate density contrast [kg/m^3]
- `g_acc`: Gravity acceleration magnitude [m/s^2]
- `eta_susp`: Viscosity of the silicate suspension [Pa s]

# Keyword Arguments
- `percolation_active`: Enable porous flow percolation (default false)
- `settling_active`: Enable Stokes gravitational settling (default false)
- `k_metal_ref`: Reference metal permeability [m^2] (default 1.0e-9)
- `eta_metal`: Dynamic viscosity of molten metal [Pa s] (default 1.0e-2)
- `phi_crit_perc`: Percolation connectivity threshold [-] (default 0.05)
- `phi_residual`: Capillary-trapped residual metal fraction [-] (default 0.02)
- `phi0`: Reference porosity scale [-] (default 0.1)
- `perm_exponent`: Porosity exponent for permeability [-] (default 3.0)
- `r_drop`: Settling droplet radius [m] (default 5.0e-3)
- `hindered_exponent`: Richardson-Zaki exponent [-] (default 4.5)
- `phi_pack`: Packing volume fraction [-] (default 0.65)
- `hadamard_rybczynski`: Fluid droplet correction (default false)
- `F_settle_start`: Lower silicate melt fraction bound for settling (default 0.40)
- `F_perc_end`: Upper silicate melt fraction bound for percolation (default 0.50)

# Returns
- Segregation velocity magnitude [m/s].
"""
function metal_segregation_velocity(
    phi_m::Real,
    F_m::Real,
    drho::Real,
    g_acc::Real,
    eta_susp::Real;
    percolation_active::Bool=false,
    settling_active::Bool=false,
    k_metal_ref::Real=1.0e-9,
    eta_metal::Real=1.0e-2,
    phi_crit_perc::Real=0.05,
    phi_residual::Real=0.02,
    phi0::Real=0.1,
    perm_exponent::Real=3.0,
    r_drop::Real=5.0e-3,
    hindered_exponent::Real=4.5,
    phi_pack::Real=0.65,
    hadamard_rybczynski::Bool=false,
    F_settle_start::Real=0.40,
    F_perc_end::Real=0.50,
)
    if !(0.0 <= phi_m <= 1.0) || !isfinite(phi_m)
        throw(DomainError(phi_m, "Metal volume fraction must be in [0, 1] and finite"))
    end
    if !(0.0 <= F_m <= 1.0) || !isfinite(F_m)
        throw(DomainError(F_m, "Silicate melt fraction must be in [0, 1] and finite"))
    end
    if !isfinite(drho) || drho < 0.0
        throw(DomainError(drho, "Density contrast must be non-negative and finite"))
    end
    if !isfinite(g_acc) || g_acc < 0.0
        throw(DomainError(g_acc, "Gravity magnitude must be non-negative and finite"))
    end
    if !isfinite(eta_susp) || eta_susp <= 0.0
        throw(DomainError(eta_susp, "Suspension viscosity must be positive and finite"))
    end
    if percolation_active
        if !(0.0 <= phi_residual <= phi_crit_perc) || !isfinite(phi_residual)
            throw(
                DomainError(
                    phi_residual,
                    "Residual metal fraction must be non-negative and <= phi_crit_perc",
                ),
            )
        end
    end

    if (!percolation_active && !settling_active) || iszero(phi_m)
        return 0.0
    end

    phi_f = Float64(phi_m)
    F_f = Float64(F_m)
    drho_f = Float64(drho)
    g_f = Float64(g_acc)

    # Percolation component
    v_perc = _percolation_interstitial_velocity(
        phi_f,
        F_f,
        drho_f,
        g_f,
        Float64(phi_residual);
        percolation_active=percolation_active,
        F_perc_end=Float64(F_perc_end),
        phi_crit_perc=Float64(phi_crit_perc),
        k_metal_ref=Float64(k_metal_ref),
        eta_metal=Float64(eta_metal),
        perm_exponent=Float64(perm_exponent),
        phi0=Float64(phi0),
    )

    # Settling component
    v_settle = if settling_active && F_f >= F_settle_start && phi_f > 0.0
        v0 = stokes_settling_velocity(
            r_drop,
            drho_f,
            g_f,
            eta_susp;
            hadamard_rybczynski=hadamard_rybczynski,
            eta_metal=eta_metal,
        )
        h = richardson_zaki_hindrance(
            phi_f; hindered_exponent=hindered_exponent, phi_pack=phi_pack
        )
        v0 * h
    else
        0.0
    end

    if percolation_active && !settling_active
        return v_perc
    elseif !percolation_active && settling_active
        return v_settle
    else
        w = _cubic_hermite_weight(F_f, F_settle_start, F_perc_end)
        return (1.0 - w) * v_perc + w * v_settle
    end
end

"""
Compute volumetric heating rate from gravitational potential energy dissipation.

$(SIGNATURES)

References:
- Stevenson (1990), Origin of the Earth, Oxford Univ. Press.
- Ricard et al. (2009), J. Geophys. Res., 114, B07404.

# Arguments
- `phi_m`: Metal volume fraction [-]
- `drho`: Density contrast [kg/m^3]
- `g_acc`: Gravity acceleration [m/s^2]
- `v_seg`: Segregation velocity magnitude [m/s]

# Returns
- Volumetric heating rate [W/m^3].

# Raises
- `DomainError`: If fractions or inputs are negative or unphysical.
"""
function segregation_dissipation_heating(phi_m::Real, drho::Real, g_acc::Real, v_seg::Real)
    if !(0.0 <= phi_m <= 1.0) || !isfinite(phi_m)
        throw(DomainError(phi_m, "Metal volume fraction must be in [0, 1] and finite"))
    end
    if drho < 0.0 || !isfinite(drho)
        throw(DomainError(drho, "Density contrast must be non-negative and finite"))
    end
    if g_acc < 0.0 || !isfinite(g_acc)
        throw(DomainError(g_acc, "Gravity magnitude must be non-negative and finite"))
    end
    if v_seg < 0.0 || !isfinite(v_seg)
        throw(DomainError(v_seg, "Segregation velocity must be non-negative and finite"))
    end

    return phi_m * drho * g_acc * v_seg
end

"""
Compute volume-weighted blended density including liquid metal.

$(SIGNATURES)

# Arguments
- `rho_silicate`: Silicate matrix density [kg/m^3]
- `rho_metal`: Molten metal density [kg/m^3]
- `phi_m`: Metal volume fraction [-]

# Returns
- Effective composite density [kg/m^3].
"""
function metal_blended_density(rho_silicate::Real, rho_metal::Real, phi_m::Real)
    if rho_silicate <= 0.0 || !isfinite(rho_silicate)
        throw(DomainError(rho_silicate, "Silicate density must be positive and finite"))
    end
    if rho_metal <= 0.0 || !isfinite(rho_metal)
        throw(DomainError(rho_metal, "Metal density must be positive and finite"))
    end
    if !(0.0 <= phi_m <= 1.0) || !isfinite(phi_m)
        throw(DomainError(phi_m, "Metal volume fraction must be in [0, 1] and finite"))
    end

    return (1.0 - phi_m) * rho_silicate + phi_m * rho_metal
end

"""
Compute blended thermal conductivity including metal phase.

$(SIGNATURES)

# Arguments
- `k_silicate`: Silicate thermal conductivity [W/(m K)]
- `k_metal`: Metal thermal conductivity [W/(m K)]
- `phi_m`: Metal volume fraction [-]

# Keyword Arguments
- `mode`: Blending method (`:arithmetic` or `:geometric`) (default `:arithmetic`)

# Returns
- Effective thermal conductivity [W/(m K)].
"""
function metal_blended_conductivity(
    k_silicate::Real, k_metal::Real, phi_m::Real; mode::Symbol=:arithmetic
)
    if k_silicate <= 0.0 || !isfinite(k_silicate)
        throw(DomainError(k_silicate, "Silicate conductivity must be positive and finite"))
    end
    if k_metal <= 0.0 || !isfinite(k_metal)
        throw(DomainError(k_metal, "Metal conductivity must be positive and finite"))
    end
    if !(0.0 <= phi_m <= 1.0) || !isfinite(phi_m)
        throw(DomainError(phi_m, "Metal volume fraction must be in [0, 1] and finite"))
    end

    if mode === :arithmetic
        return (1.0 - phi_m) * k_silicate + phi_m * k_metal
    elseif mode === :geometric
        return (k_silicate^(1.0 - phi_m)) * (k_metal^phi_m)
    else
        throw(ArgumentError("mode must be :arithmetic or :geometric, got $mode"))
    end
end

"""
Compute volume-weighted volumetric heat capacity including metal phase.

$(SIGNATURES)

# Arguments
- `rhocp_silicate`: Silicate volumetric heat capacity [J/(m^3 K)]
- `rhocp_metal`: Metal volumetric heat capacity [J/(m^3 K)]
- `phi_m`: Metal volume fraction [-]

# Returns
- Effective volumetric heat capacity [J/(m^3 K)].
"""
function metal_blended_heat_capacity(rhocp_silicate::Real, rhocp_metal::Real, phi_m::Real)
    if rhocp_silicate <= 0.0 || !isfinite(rhocp_silicate)
        throw(
            DomainError(
                rhocp_silicate, "Silicate heat capacity must be positive and finite"
            ),
        )
    end
    if rhocp_metal <= 0.0 || !isfinite(rhocp_metal)
        throw(DomainError(rhocp_metal, "Metal heat capacity must be positive and finite"))
    end
    if !(0.0 <= phi_m <= 1.0) || !isfinite(phi_m)
        throw(DomainError(phi_m, "Metal volume fraction must be in [0, 1] and finite"))
    end

    return (1.0 - phi_m) * rhocp_silicate + phi_m * rhocp_metal
end

"""
Compute Rouse number comparing Stokes settling to turbulent convective mixing.

$(SIGNATURES)

References:
- Deguen et al. (2011, 2014), Earth Planet. Sci. Lett., 310(3-4), 308-318.

# Arguments
- `v_settle`: Stokes settling velocity [m/s]
- `u_conv`: Convective RMS velocity [m/s]

# Returns
- Rouse number [-] (> 1 indicates rainout, < 1 indicates suspension).
"""
function suspension_rouse_number(v_settle::Real, u_conv::Real)
    if v_settle < 0.0 || !isfinite(v_settle)
        throw(DomainError(v_settle, "Settling velocity must be non-negative and finite"))
    end
    if u_conv <= 0.0 || !isfinite(u_conv)
        throw(DomainError(u_conv, "Convective velocity must be positive and finite"))
    end

    return v_settle / u_conv
end

"""
Compute density of liquid Fe-FeS metallic melt as a function of sulfur mass fraction.

$(SIGNATURES)

References:
- Sanloup et al. (2000), Geophys. Res. Lett., 27(6), 811-814.
- Morard et al. (2014), C. R. Geoscience, 346(5-6), 130-139.

# Arguments
- `w_S::Real`: Sulfur mass fraction in metallic melt [-] (e.g., 0.0 for pure iron, 0.31 for Fe-FeS eutectic, 0.365 for stoichiometric FeS)

# Keyword Arguments
- `T::Real=1500.0`: Melt temperature [K] (default: 1500.0)
- `P::Real=0.0`: Confining pressure [Pa] (default: 0.0)
- `law::Symbol=:sanloup2000`: Density parameterization (`:sanloup2000` or `:morard2014`)
- `alpha_m::Real=1.0e-4`: Volumetric thermal expansion coefficient [1/K] (default: 1.0e-4)
- `K_T::Real=6.5e10`: Isothermal bulk modulus [Pa] (default: 65.0 GPa)
- `T0::Real=1500.0`: Reference temperature for 1-bar density calibration [K] (default: 1500.0)

# Returns
- `rho_metal::Float64`: Liquid metal density [kg/m^3].

# Raises
- `DomainError`: If `w_S` is not in [0.0, 0.40], `T <= 0.0`, or `P < 0.0`.
- `ArgumentError`: If `law` is not `:sanloup2000` or `:morard2014`.
"""
function compute_liquid_metal_density(
    w_S::Real;
    T::Real=1500.0,
    P::Real=0.0,
    law::Symbol=:sanloup2000,
    alpha_m::Real=1.0e-4,
    K_T::Real=6.5e10,
    T0::Real=1500.0,
)::Float64
    w_val = Float64(w_S)
    T_val = Float64(T)
    P_val = Float64(P)
    if !(0.0 <= w_val <= 0.40) || !isfinite(w_val)
        throw(DomainError(w_val, "Sulfur mass fraction must be in [0, 0.40] and finite"))
    end
    if T_val <= 0.0 || !isfinite(T_val)
        throw(DomainError(T_val, "Temperature must be positive and finite"))
    end
    if P_val < 0.0 || !isfinite(P_val)
        throw(DomainError(P_val, "Pressure must be non-negative and finite"))
    end

    rho0 = if law === :sanloup2000
        # Sanloup et al. (2000) linear calibration: pure Fe ~7020 kg/m^3 at 1500 K,
        # dropping to ~5450 kg/m^3 at eutectic w_S = 0.31 (~5050 kg/m^3 per unit w_S).
        7020.0 - 5050.0 * w_val
    elseif law === :morard2014
        # Morard et al. (2014) relative deficit model: rho = rho_Fe * (1 - 0.72 * w_S)
        7020.0 * (1.0 - 0.72 * w_val)
    else
        throw(
            ArgumentError(
                "Unknown liquid metal density law: $law. Supported: :sanloup2000, :morard2014",
            ),
        )
    end

    # Thermal expansion and compressibility correction:
    # Floor of 0.5 prevents unphysical negative/zero density during extreme numerical transients.
    thermal_factor = max(0.5, 1.0 - alpha_m * (T_val - T0) + P_val / K_T)
    return rho0 * thermal_factor
end

"""
Compute silicate suspension viscosity with smooth melt weakening and crystal-mush cutoff.

$(SIGNATURES)

References:
- Roscoe (1952), Br. J. Appl. Phys., 3(8), 267-269.
- Solomatov (2015), Treatise on Geophysics, Ch. 9.04, Eq. 17.

# Arguments
- `eta_solid`: Reference silicate matrix viscosity at `F_m = 0` [Pa s]
- `F_m`: Silicate melt fraction in `[0, 1]` [-]

# Keyword Arguments
- `alpha_eta`: Arrhenius melt-weakening exponent below `phi_crit` [-] (default `28.0`)
- `phi_crit`: Rheological critical melt fraction [-] (default `0.4`)
- `eta_melt`: Pure silicate melt viscosity at `F_m = 1` [Pa s] (default `10.0`)
- `etamin`: Lower viscosity clamp [Pa s] (default `0.01`)
- `etamax`: Upper viscosity clamp [Pa s] (default `1.0e23`)
- `F_blend_width`: Melt-fraction width for Hermite blend into Roscoe branch [-] (default `0.10`)

# Returns
- Effective silicate suspension viscosity [Pa s] clamped to `[etamin, etamax]`.

# Raises
- `DomainError`: If any argument is non-finite or outside its physical domain.
"""
function compute_silicate_suspension_viscosity(
    eta_solid::Real,
    F_m::Real;
    alpha_eta::Real=28.0,
    phi_crit::Real=0.4,
    eta_melt::Real=10.0,
    etamin::Real=0.01,
    etamax::Real=1.0e23,
    F_blend_width::Real=0.10,
)::Float64
    eta_s = Float64(eta_solid)
    Fm = Float64(F_m)
    alpha_e = Float64(alpha_eta)
    phi_c = Float64(phi_crit)
    eta_m = Float64(eta_melt)
    etamin_f = Float64(etamin)
    etamax_f = Float64(etamax)
    F_bw = Float64(F_blend_width)

    if !isfinite(eta_s) || eta_s <= 0.0
        throw(DomainError(eta_s, "Solid viscosity must be positive and finite"))
    end
    if !isfinite(Fm) || !(0.0 <= Fm <= 1.0)
        throw(DomainError(Fm, "Silicate melt fraction must be in [0, 1] and finite"))
    end
    if !isfinite(alpha_e) || alpha_e < 0.0
        throw(
            DomainError(alpha_e, "Melt-weakening exponent must be non-negative and finite")
        )
    end
    if !isfinite(phi_c) || !(0.0 < phi_c < 1.0)
        throw(DomainError(phi_c, "Critical melt fraction must be in (0, 1) and finite"))
    end
    if !isfinite(eta_m) || eta_m <= 0.0
        throw(DomainError(eta_m, "Melt viscosity must be positive and finite"))
    end
    if !isfinite(etamin_f) || etamin_f <= 0.0
        throw(DomainError(etamin_f, "Minimum viscosity must be positive and finite"))
    end
    if !isfinite(etamax_f) || etamax_f < etamin_f
        throw(DomainError(etamax_f, "Maximum viscosity must be >= etamin and finite"))
    end
    if !isfinite(F_bw) || F_bw <= 0.0
        throw(DomainError(F_bw, "Blend width must be positive and finite"))
    end

    eta_raw = if Fm <= phi_c
        exp(log(eta_s) - alpha_e * Fm)
    else
        pack_denom = max(1.0 - (1.0 - Fm) / (1.0 - phi_c), 1.0e-12)
        log_roscoe = log(eta_m) - 2.5 * log(pack_denom)
        log_crit = log(eta_s) - alpha_e * phi_c
        log_target = min(log_roscoe, log_crit)
        F_end = min(phi_c + F_bw, max(0.95, 0.5 * (phi_c + 1.0)))
        if Fm < F_end && F_end > phi_c
            theta = (Fm - phi_c) / (F_end - phi_c)
            w = theta * theta * (3.0 - 2.0 * theta)
            exp((1.0 - w) * log_crit + w * log_target)
        else
            exp(log_target)
        end
    end
    return clamp(eta_raw, etamin_f, etamax_f)
end

"""
Compute the effective trapped metal volume fraction across the percolation-to-settling transition.

$(SIGNATURES)

# Arguments
- `F_m`: Silicate melt fraction in `[0, 1]` [-]
- `phi_residual`: Residual percolation metal fraction in `[0, 1]` [-]

# Keyword Arguments
- `percolation_active`: Enable solid-matrix percolation trapping (default `true`)
- `settling_active`: Enable magma-ocean droplet settling release (default `true`)
- `F_settle_start`: Melt fraction where settling transition starts [-] (default `0.40`)
- `F_perc_end`: Melt fraction where percolation trapping vanishes [-] (default `0.50`)

# Returns
- Effective irreducible trapped metal volume fraction [-].

# Raises
- `DomainError`: If fractions or transition bounds are non-finite or out of order.
"""
function effective_metal_trap_fraction(
    F_m::Real,
    phi_residual::Real;
    percolation_active::Bool=true,
    settling_active::Bool=true,
    F_settle_start::Real=0.40,
    F_perc_end::Real=0.50,
)::Float64
    Fm = Float64(F_m)
    phi_r = Float64(phi_residual)
    F0 = Float64(F_settle_start)
    F1 = Float64(F_perc_end)

    if !isfinite(Fm) || !(0.0 <= Fm <= 1.0)
        throw(DomainError(Fm, "Silicate melt fraction must be in [0, 1] and finite"))
    end
    if !isfinite(phi_r) || !(0.0 <= phi_r <= 1.0)
        throw(DomainError(phi_r, "Residual metal fraction must be in [0, 1] and finite"))
    end
    if !isfinite(F0) || !isfinite(F1) || !(0.0 <= F0 <= F1 <= 1.0)
        throw(
            DomainError(
                (F0, F1),
                "Transition bounds must satisfy 0 <= F_settle_start <= F_perc_end <= 1",
            ),
        )
    end

    !percolation_active && return 0.0
    !settling_active && return phi_r
    w = _cubic_hermite_weight(Fm, F0, F1)
    return (1.0 - w) * phi_r
end

"""
Compute equivalent spherical metallic core radius for a packed core mass.

Converts a compact metallic mass `M_core_metal` at density `rho_metal` and packing
fraction `phi_pack` into an equivalent sphere radius. Run diagnostics
`state.accumulators.rcore` instead compute the volume-equivalent sphere radius from
summed 3D core marker volumes `V_core_total` where `Xfe_bulk >= phi_core_threshold`.

$(SIGNATURES)

# Arguments
- `M_core_metal`: Segregated core metal mass [kg]
- `rho_metal`: Metallic core density [kg/m^3]
- `phi_pack`: Core metal packing volume fraction in `(0, 1]` [-] (default `0.65`)

# Returns
- Equivalent spherical metallic core radius [m].

# Raises
- `DomainError`: If mass is negative, density is non-positive, or packing is outside `(0, 1]`.
"""
function compute_core_radius(
    M_core_metal::Real, rho_metal::Real, phi_pack::Real=0.65
)::Float64
    M_val = Float64(M_core_metal)
    rho_val = Float64(rho_metal)
    phi_val = Float64(phi_pack)

    if !isfinite(M_val) || M_val < 0.0
        throw(DomainError(M_val, "Core metal mass must be non-negative and finite"))
    end
    if !isfinite(rho_val) || rho_val <= 0.0
        throw(DomainError(rho_val, "Core metal density must be positive and finite"))
    end
    if !isfinite(phi_val) || !(0.0 < phi_val <= 1.0)
        throw(DomainError(phi_val, "Core packing fraction must be in (0, 1] and finite"))
    end

    return iszero(M_val) ? 0.0 : cbrt(3.0 * M_val / (4.0 * π * rho_val * phi_val))
end

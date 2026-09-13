"""
Compute porous Rayleigh-Darcy number Ra_m for porous hydrothermal convection.

Calculates the dimensionless Rayleigh-Darcy number (Horton and Rogers, 1945; Lapwood, 1948):
    Ra_m = (rho_f^2 * cp_f * g * alpha_f * K * dT * H) / (mu_f * k_cond)

$(SIGNATURES)

# Arguments
- `rho_f`: Fluid density [kg/m³]
- `cp_f`: Fluid isobaric heat capacity [J/(kg K)]
- `g`: Local gravitational acceleration [m/s²]
- `alpha_f`: Fluid isobaric thermal expansivity [1/K]
- `K`: Medium permeability [m²]
- `dT`: Characteristic temperature contrast ΔT = max(0, T - T_surface) [K]
- `H`: Characteristic layer thickness [m]
- `mu_f`: Dynamic fluid viscosity [Pa s]
- `k_cond`: Bulk conductive thermal conductivity [W/(m K)]

# Returns
- `Ra_m`: Dimensionless Rayleigh-Darcy number [-]
"""
function compute_porous_rayleigh_darcy(
    rho_f::Real,
    cp_f::Real,
    g::Real,
    alpha_f::Real,
    K::Real,
    dT::Real,
    H::Real,
    mu_f::Real,
    k_cond::Real,
)
    if !isfinite(rho_f) ||
        rho_f <= 0.0 ||
        !isfinite(cp_f) ||
        cp_f <= 0.0 ||
        !isfinite(g) ||
        g < 0.0 ||
        !isfinite(alpha_f) ||
        alpha_f < 0.0 ||
        !isfinite(K) ||
        K < 0.0 ||
        !isfinite(dT) ||
        !isfinite(H) ||
        H < 0.0 ||
        !isfinite(mu_f) ||
        mu_f <= 0.0 ||
        !isfinite(k_cond) ||
        k_cond <= 0.0
        throw(
            DomainError(
                (rho_f, cp_f, g, alpha_f, K, dT, H, mu_f, k_cond),
                "Thermodynamic parameters and properties must be positive, finite, and within physical bounds",
            ),
        )
    end
    if dT <= 0.0 || K <= 0.0 || H <= 0.0 || g <= 0.0
        return 0.0
    end
    num = rho_f^2 * cp_f * g * alpha_f * K * dT * H
    den = mu_f * k_cond
    return Float64(num / den)
end

"""
Compute free-fluid thermal Rayleigh number Ra for open fluid convection.

Calculates the dimensionless Rayleigh number (Rayleigh, 1916; Kraichnan, 1962):
    Ra = (rho_f^2 * cp_f * g * alpha_f * dT * H^3) / (mu_f * k_f)

$(SIGNATURES)

# Arguments
- `rho_f`: Fluid density [kg/m³]
- `cp_f`: Fluid isobaric heat capacity [J/(kg K)]
- `g`: Local gravitational acceleration [m/s²]
- `alpha_f`: Fluid isobaric thermal expansivity [1/K]
- `dT`: Characteristic temperature contrast ΔT = max(0, T - T_surface) [K]
- `H`: Characteristic layer thickness [m]
- `mu_f`: Dynamic fluid viscosity [Pa s]
- `k_f`: Fluid thermal conductivity [W/(m K)]

# Returns
- `Ra`: Dimensionless Rayleigh number [-]
"""
function compute_free_fluid_rayleigh(
    rho_f::Real,
    cp_f::Real,
    g::Real,
    alpha_f::Real,
    dT::Real,
    H::Real,
    mu_f::Real,
    k_f::Real,
)
    if !isfinite(rho_f) ||
        rho_f <= 0.0 ||
        !isfinite(cp_f) ||
        cp_f <= 0.0 ||
        !isfinite(g) ||
        g < 0.0 ||
        !isfinite(alpha_f) ||
        alpha_f < 0.0 ||
        !isfinite(dT) ||
        !isfinite(H) ||
        H < 0.0 ||
        !isfinite(mu_f) ||
        mu_f <= 0.0 ||
        !isfinite(k_f) ||
        k_f <= 0.0
        throw(
            DomainError(
                (rho_f, cp_f, g, alpha_f, dT, H, mu_f, k_f),
                "Thermodynamic parameters and properties must be positive, finite, and within physical bounds",
            ),
        )
    end
    if dT <= 0.0 || H <= 0.0 || g <= 0.0
        return 0.0
    end
    num = rho_f^2 * cp_f * g * alpha_f * dT * H^3
    den = mu_f * k_f
    return Float64(num / den)
end

"""
Compute subgrid hydrothermal Nusselt number blending porous and free-fluid regimes.

Interpolates between porous Rayleigh-Darcy scaling (Lapwood, 1948) and boundary-layer
free-fluid Rayleigh scaling (Kraichnan, 1962; Howard, 1966) across a cubic smoothstep
porosity transition [phi_start, phi_end]:
    Nu_porous = 1.0 + c_porous * max(0.0, Ra_m / Ra_m_crit - 1.0)
    Nu_free = max(1.0, c_free * cbrt(Ra))
    xi = clamp((phi - phi_start) / (phi_end - phi_start), 0.0, 1.0)
    w_phi = xi^2 * (3.0 - 2.0 * xi)
    log10(Nu) = (1.0 - w_phi) * log10(Nu_porous) + w_phi * log10(Nu_free)

$(SIGNATURES)

# Arguments
- `Ra_m`: Dimensionless porous Rayleigh-Darcy number [-]
- `Ra`: Dimensionless free-fluid Rayleigh number [-]
- `phi`: Current medium porosity [-]

# Keyword Arguments
- `phi_start`: Lower porosity threshold for transition onset (default: 0.30)
- `phi_end`: Upper porosity threshold for full free-fluid regime (default: 0.70)
- `Ra_m_crit`: Critical Rayleigh-Darcy number for porous convection onset (default: 4π²)
- `Ra_crit`: Critical Rayleigh number for free-fluid convection onset (default: 1100.0)
- `c_porous`: Prefactor for porous convective Nusselt scaling (default: 1.0)
- `c_free`: Prefactor for free-fluid boundary layer Nusselt scaling (default: 0.088)

# Returns
- `Nu`: Dimensionless Nusselt number (Nu >= 1.0) [-]
"""
function compute_hydrothermal_nusselt(
    Ra_m::Real,
    Ra::Real,
    phi::Real;
    phi_start::Real=0.30,
    phi_end::Real=0.70,
    Ra_m_crit::Real=4.0 * pi^2,
    Ra_crit::Real=1100.0,
    c_porous::Real=1.0,
    c_free::Real=0.088,
)
    if !isfinite(Ra_m) ||
        Ra_m < 0.0 ||
        !isfinite(Ra) ||
        Ra < 0.0 ||
        !isfinite(phi) ||
        !(0.0 <= phi <= 1.0)
        throw(
            DomainError(
                (Ra_m, Ra, phi),
                "Rayleigh numbers must be non-negative and porosity in [0, 1]",
            ),
        )
    end
    if !isfinite(phi_start) || !isfinite(phi_end) || phi_end <= phi_start
        throw(
            DomainError(
                (phi_start, phi_end),
                "phi_end must be strictly greater than phi_start and both finite",
            ),
        )
    end
    if !isfinite(Ra_m_crit) ||
        Ra_m_crit <= 0.0 ||
        !isfinite(Ra_crit) ||
        Ra_crit <= 0.0 ||
        !isfinite(c_porous) ||
        c_porous <= 0.0 ||
        !isfinite(c_free) ||
        c_free <= 0.0
        throw(
            DomainError(
                (Ra_m_crit, Ra_crit, c_porous, c_free),
                "Scaling parameters must be strictly positive and finite",
            ),
        )
    end

    # Porous Nusselt scaling (Horton-Rogers-Lapwood criterion)
    Nu_porous = if Ra_m <= Ra_m_crit
        1.0
    else
        1.0 + c_porous * (Ra_m / Ra_m_crit - 1.0)
    end

    # Free-fluid boundary-layer Nusselt scaling (Kraichnan / Howard scaling)
    Nu_free = if Ra <= Ra_crit
        1.0
    else
        max(1.0, c_free * cbrt(Ra))
    end

    # Smoothstep porosity blending
    if phi <= phi_start
        return Float64(Nu_porous)
    elseif phi >= phi_end
        return Float64(Nu_free)
    else
        xi = (phi - phi_start) / (phi_end - phi_start)
        w_phi = xi * xi * (3.0 - 2.0 * xi)
        log_Nu = (1.0 - w_phi) * log10(Nu_porous) + w_phi * log10(Nu_free)
        return Float64(10.0^log_Nu)
    end
end

"""
Compute effective thermal conductivity enhancement from hydrothermal convection.

Applies cell-Péclet resolution weighting to avoid double-counting convective flux
when resolved on the grid, clamps within physical floors and cutoffs, and applies
Picard relaxation damping:
    k_target = max(k_cond, clamp(Nu * k_cond, k_floor, k_cutoff))
    w_res = clamp(Pe_cell / Pe_crit, 0.0, 1.0)
    k_damped = (1.0 - w_res) * k_target + w_res * k_cond
    k_eff = (1.0 - gamma) * k_prev + gamma * k_damped

$(SIGNATURES)

# Arguments
- `k_cond`: Conductive baseline thermal conductivity [W/(m K)]
- `Nu`: Dimensionless hydrothermal Nusselt number [-]

# Keyword Arguments
- `k_floor`: Minimum thermal conductivity floor [W/(m K)] (default: 1.0e-3)
- `k_cutoff`: Maximum enhanced thermal conductivity ceiling [W/(m K)] (default: 1.0e6)
- `Pe_cell`: Cell-Péclet number of resolved Darcy flux [-] (default: 0.0)
- `Pe_crit`: Critical cell-Péclet number for full resolved damping [-] (default: 2.0)
- `resolution_weighting`: Enable grid resolution damping (default: true)
- `picard_damping`: Picard relaxation damping parameter γ in (0, 1] (default: 1.0)
- `k_prev`: Thermal conductivity iterate from previous step/iteration [W/(m K)] (default: k_cond)

# Returns
- `k_eff`: Effective thermal conductivity [W/(m K)]
"""
function compute_effective_hydrothermal_conductivity(
    k_cond::Real,
    Nu::Real;
    k_floor::Real=1.0e-3,
    k_cutoff::Real=1.0e6,
    Pe_cell::Real=0.0,
    Pe_crit::Real=2.0,
    resolution_weighting::Bool=true,
    picard_damping::Real=1.0,
    k_prev::Real=k_cond,
)
    if !isfinite(k_cond) || k_cond <= 0.0 || !isfinite(Nu) || Nu < 1.0
        throw(
            DomainError((k_cond, Nu), "k_cond must be positive and finite, and Nu >= 1.0")
        )
    end
    if !isfinite(k_floor) || !isfinite(k_cutoff) || k_cutoff <= k_floor || k_floor <= 0.0
        throw(
            DomainError(
                (k_floor, k_cutoff),
                "Conductivity bounds must satisfy 0 < k_floor < k_cutoff and be finite",
            ),
        )
    end
    if !isfinite(Pe_cell) || Pe_cell < 0.0 || !isfinite(Pe_crit) || Pe_crit <= 0.0
        throw(
            DomainError(
                (Pe_cell, Pe_crit),
                "Pe_cell must be non-negative and Pe_crit positive and finite",
            ),
        )
    end
    if !isfinite(picard_damping) || !(0.0 < picard_damping <= 1.0)
        throw(DomainError(picard_damping, "picard_damping must be in (0, 1] and finite"))
    end
    if !isfinite(k_prev) || k_prev <= 0.0
        throw(DomainError(k_prev, "k_prev must be positive and finite"))
    end

    # Target convective conductivity
    k_raw = Nu * k_cond
    k_target = max(k_cond, clamp(k_raw, k_floor, k_cutoff))

    # Grid-resolution weighting (damps subgrid enhancement if flow is resolved)
    k_res = if resolution_weighting && Pe_cell > 0.0
        w_res = clamp(Pe_cell / Pe_crit, 0.0, 1.0)
        (1.0 - w_res) * k_target + w_res * k_cond
    else
        k_target
    end

    # Picard relaxation damping
    gamma = Float64(picard_damping)
    k_damped = (1.0 - gamma) * Float64(k_prev) + gamma * Float64(k_res)

    return Float64(max(k_cond, clamp(k_damped, k_floor, k_cutoff)))
end

"""
Apply complete subgrid hydrothermal convection closure to marker conductivity.

Evaluates fluid thermodynamic properties, porous Rayleigh-Darcy and free Rayleigh
numbers, smoothstep Nusselt number, temperature regularization, cell-Péclet damping,
and Picard relaxation damping.

$(SIGNATURES)

# Arguments
- `k_cond`: Bulk conductive thermal conductivity [W/(m K)]
- `T`: Local temperature [K]
- `phi`: Local porosity [-]
- `tm`: Material type index (1: core, 2: crust, 3: sticky air)

# Keyword Arguments
- `cfg`: HydrothermalConfig parameter set
- `K`: Local permeability [m²] (default: -1.0, calculates from Kozeny-Carman)
- `Pe_cell`: Cell-Péclet number for resolved Darcy flow (default: 0.0)
- `k_prev`: Conductivity iterate from previous Picard step [W/(m K)] (default: k_cond)
- `tmfluidphase_val`: Melting temperature of pore fluid [K] (default: 273.15)

# Returns
- `k_eff`: Effective hydrothermal thermal conductivity [W/(m K)]
"""
function apply_hydrothermal_convection_closure(
    k_cond::Real,
    T::Real,
    phi::Real,
    tm::Integer;
    cfg::HydrothermalConfig=HydrothermalConfig(),
    K::Real=-1.0,
    Pe_cell::Real=0.0,
    k_prev::Real=k_cond,
    tmfluidphase_val::Real=273.15,
)
    # Physical domain contracts: validate inputs before checking activity or temperature thresholds
    if !isfinite(T) || T <= 0.0
        throw(DomainError(T, "Temperature must be positive and finite"))
    end
    if !isfinite(phi) || !(0.0 <= phi <= 1.0)
        throw(DomainError(phi, "Porosity must be in [0, 1] and finite"))
    end
    if !isfinite(k_cond) || k_cond <= 0.0
        throw(DomainError(k_cond, "Conductive conductivity must be positive and finite"))
    end

    if !cfg.active || tm >= 3
        return Float64(k_cond)
    end
    if T <= tmfluidphase_val || T <= cfg.T_surface_ref
        return Float64(k_cond)
    end

    dT = max(0.0, T - cfg.T_surface_ref)
    if dT <= 0.0
        return Float64(k_cond)
    end

    # Permeability evaluation: use provided K or compute from Kozeny-Carman
    K_eff = if K >= 0.0
        Float64(K)
    else
        if phi <= 0.0
            0.0
        elseif phi < 1.0
            kphi(cfg.kphi_ref, phi)
        else
            cfg.kphi_ref
        end
    end

    # Thermodynamic fluid properties
    rho_f = compute_rhofluid(T, cfg.rho_fluid_ref, cfg.alpha_fluid, cfg.T_surface_ref)
    mu_f = compute_fluid_viscosity(
        T,
        tm;
        mode=:arrhenius,
        eta0=cfg.mu_fluid_ref,
        T0=cfg.T_surface_ref,
        tmfluidphase=tmfluidphase_val,
    )

    # Rayleigh numbers
    Ra_m = compute_porous_rayleigh_darcy(
        rho_f,
        cfg.cp_fluid,
        cfg.gravity,
        cfg.alpha_fluid,
        K_eff,
        dT,
        cfg.H_layer,
        mu_f,
        k_cond,
    )
    Ra = compute_free_fluid_rayleigh(
        rho_f,
        cfg.cp_fluid,
        cfg.gravity,
        cfg.alpha_fluid,
        dT,
        cfg.H_layer,
        mu_f,
        cfg.k_fluid_ref,
    )

    # Blended Nusselt number
    Nu = compute_hydrothermal_nusselt(
        Ra_m,
        Ra,
        phi;
        phi_start=cfg.phi_start,
        phi_end=cfg.phi_end,
        Ra_m_crit=cfg.Ra_m_crit,
        Ra_crit=cfg.Ra_crit,
        c_porous=cfg.c_porous,
        c_free=cfg.c_free,
    )

    # Quadratic temperature ramp across dT_min for C¹ continuity at surface boundary
    w_T = clamp(dT / cfg.dT_min, 0.0, 1.0)^2
    Nu_ramped = if w_T <= 0.0
        1.0
    elseif w_T >= 1.0
        Nu
    else
        10.0^(w_T * log10(Nu))
    end

    return compute_effective_hydrothermal_conductivity(
        k_cond,
        Nu_ramped;
        k_floor=cfg.k_floor,
        k_cutoff=cfg.k_cutoff,
        Pe_cell=Pe_cell,
        Pe_crit=cfg.Pe_crit,
        resolution_weighting=cfg.resolution_weighting,
        picard_damping=cfg.picard_damping,
        k_prev=k_prev,
    )
end

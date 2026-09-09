
"""
Calculate Euclidean distance between two point coordinates.

$(SIGNATURES)

# Details

    - x1: x-coordinate of point 1 [m]
    - y1: y-coordinate of point 1 [m]
    - x2: x-coordinate of point 2 [m]
    - y2: y-coordinate of point 2 [m]

# Returns

    - Euclidean distance between point 1 and point 2 [m]
"""
function distance(x1, y1, x2, y2)
    return sqrt(abs2(x1-x2) + abs2(y1-y2))
end

"""
Compute convex combination of fluid and solid properties to get total property.

$(SIGNATURES)

# Details

    - fluid: fluid properties
    - solid: solid properties
    - ϕ: porosity (fraction of fluid)

# Returns

    - total: computed total property
"""
function total(solid, fluid, ϕ)
    if !(0.0 <= ϕ <= 1.0)
        throw(DomainError(ϕ, "Porosity must be in [0, 1]"))
    end
    return solid * (1.0 - ϕ) + fluid * ϕ
end

"""
Compute total thermal conductivity of two-phase material.

$(SIGNATURES)

# Details

    - ksolid: solid thermal conductivity [W/m/K]
    - kfluid: fluid thermal conductivity [W/m/K]
    - phi: porosity (fraction of fluid)

# Returns

    - ktotal: total thermal conductivity of mixed phase [W/m/K]
"""
function ktotal(ksolid, kfluid, phi)
    if !(0.0 <= phi <= 1.0)
        throw(DomainError(phi, "Porosity must be in [0, 1]"))
    end
    if ksolid < 0.0 || kfluid < 0.0
        throw(DomainError((ksolid, kfluid), "Thermal conductivities must be non-negative"))
    end
    return (
        sqrt(
            ksolid * kfluid / 2 +
            ((ksolid * (3.0 * phi - 2.0) + kfluid * (1.0 - 3.0 * phi))^2) * inv(16.0),
        ) - 0.25 * (ksolid * (3.0 * phi - 2.0) + kfluid * (1.0 - 3.0 * phi))
    )
end

"""
Compute porosity-dependent permeability (eqn 16.64 in Gerya (2019)).

$(SIGNATURES)

# Details

    - kphim0m: standard (reference) permeability (of marker type) [m^2]
    - phimm: actual (marker) porosity

# Returns

    - kphim: empirical porosity-dependent permeability [m^2]
"""
function kphi(kphim0m, phimm)
    if !(0.0 <= phimm < 1.0)
        throw(DomainError(phimm, "Porosity must be in [0, 1)"))
    end
    if kphim0m < 0.0
        throw(DomainError(kphim0m, "Reference permeability must be non-negative"))
    end
    # phim0 is a global constant defined independent of material type
    return kphim0m * (phimm * inv(phim0))^3.0 * ((1.0 - phimm) * inv(1.0 - phim0))^-2.0
end

"""
Compute inverse of porosity-dependent permeability (eqn 16.64 in Gerya (2019)) 
times current fluid viscosity.

$(SIGNATURES)

# Details

    - kϕᵣ: reference permeability [m^2]
    - ϕ: current porosity
    - ηᶠcur: current fluid viscosity [Pa s]

# Returns

    - etafluidcur_inv_kphi: inverse empirical porosity-dependent permeability 
                            times current fluid viscosity
"""
function ηᶠcur_inv_kᵠ(kϕᵣ, ϕ, ηᶠcur)
    if !(0.0 < ϕ < 1.0)
        throw(DomainError(ϕ, "Porosity must be strictly between 0 and 1"))
    end
    if kϕᵣ <= 0.0 || ηᶠcur < 0.0
        throw(
            DomainError(
                (kϕᵣ, ηᶠcur), "Permeability must be positive and viscosity non-negative"
            ),
        )
    end
    return ηᶠcur * inv(kϕᵣ) * (phim0 * inv(ϕ))^3.0 * ((1.0 - ϕ) * inv(1.0 - phim0))^2.0
end

"""
Compute radiogenic heat production of isotope mixture.

$(SIGNATURES)

# Details

    - f: fraction of radioactive matter [atoms/kg]
    - ratio: initial ratio of radioactive to non-radioactive isotopes
    - E: heat energy [J]
    - tau: exp decay mean lifetime ``\\tau=\\frac{t_{1/2}}{\\log{2}}`` [s]
    - time: time elapsed since start of radioactive decay [s]

# Returns

    - Q: radiogenic heat production [W/kg]
"""
function Q_radiogenic(f, ratio, E, tau, time)
    if time < 0.0 || tau <= 0.0 || f < 0.0 || ratio < 0.0
        throw(
            DomainError(
                (f, ratio, E, tau, time),
                "Decay time, abundance, and ratio must be non-negative, and lifetime must be positive",
            ),
        )
    end
    return f * ratio * E * exp(-time * inv(tau)) * inv(tau)
end

"""
Compute total rocky marker viscosity based on temperature and material type.

$(SIGNATURES)

# Details

    - tkmm: marker temperature [K]
    - tmm: marker type [1, 2]

# Returns
    
    - etatotal: rocky marker temperature-dependent total viscosity 
"""
function etatotal_rocks(tkmm, tmm)
    if tkmm <= 0.0
        throw(DomainError(tkmm, "Absolute temperature must be positive"))
    end
    @inbounds etasolidcur = ifelse(tkmm > tmsolidphase, etasolidmm[tmm], etasolidm[tmm])
    @inbounds etafluidcur = ifelse(tkmm > tmfluidphase, etafluidmm[tmm], etafluidm[tmm])
    return max(etamin, etasolidcur, etafluidcur)
end

"""
Compute radiogenic heat production of 26Al and 60Fe isotopes.

$(SIGNATURES)

# Details

    - al: true if radioactive isotope 26Al is present
    - fe: true if radioactive isotope 60Fe is present
    - timesum: time elapsed since initial conditions at start of simulation

# Keyword Arguments

    - ratio_al: initial isotopic ratio of 26Al/27Al
    - E_al: decay energy of 26Al [J]
    - f_al: mass fraction of 27Al in silicate matrix
    - tau_al: mean lifetime of 26Al [s]
    - ratio_fe: initial isotopic ratio of 60Fe/56Fe
    - E_fe: decay energy of 60Fe [J]
    - f_fe: mass fraction of 56Fe in metallic phase
    - tau_fe: mean lifetime of 60Fe [s]
    - rhosolidm: density vector of solid silicate phases [kg/m^3]
    - rhofluidm: density vector of pore fluid phases [kg/m^3]
    - rho_metal: density of metallic phase [kg/m^3] (default 5450.0)

# Returns

    - hrsolidm: radiogenic heat production of 26Al in silicate rock [W/m^3]
    - hrfluidm: radiogenic heat production in pore fluid [W/m^3] (zero: pore water carries no 26Al/60Fe)
    - hrmetalm: radiogenic heat production of 60Fe in metallic iron phase [W/m^3]
"""
function calculate_radioactive_heating(
    al,
    fe,
    timesum;
    ratio_al=ratio_al,
    E_al=E_al,
    f_al=f_al,
    tau_al=tau_al,
    ratio_fe=ratio_fe,
    E_fe=E_fe,
    f_fe=f_fe,
    tau_fe=tau_fe,
    rhosolidm=rhosolidm,
    rhofluidm=rhofluidm,
    rho_metal=5450.0,
)
    # 26Al: planet ✓, crust ✓, space × (lithophile, in silicate rock phase)
    if al
        # 26Al radiogenic heat production [W/kg]
        Q_al = Q_radiogenic(f_al, ratio_al, E_al, tau_al, timesum)
        # Solid phase 26Al radiogenic heat production [W/m^3]
        @inbounds hrsolidm = @SVector [Q_al * rhosolidm[1], Q_al * rhosolidm[2], 0.0]
    else
        hrsolidm = @SVector zeros(3)
    end
    # Fluid phase: pore fluid (water/ice) carries no radiogenic isotopes
    hrfluidm = @SVector zeros(3)
    # 60Fe: planet ✓, crust ✓, space × (siderophile, resides in metallic iron phase)
    if fe
        # 60Fe radiogenic heat production [W/kg]
        Q_fe = Q_radiogenic(f_fe, ratio_fe, E_fe, tau_fe, timesum)
        # Metallic phase 60Fe radiogenic heat production [W/m^3]
        # Silicate mantle (1) and crust (2) metal carries 60Fe; sticky air/space (3) carries zero.
        @inbounds hrmetalm = @SVector [Q_fe * rho_metal, Q_fe * rho_metal, 0.0]
    else
        hrmetalm = @SVector zeros(3)
    end
    return hrsolidm, hrfluidm, hrmetalm
end

"""
Compute volumetric isobaric heat capacity of H₂O (fluid phase)
based on temperature.

$(SIGNATURES)

# Details
    
        - T: temperature [K]
        - mode: marker property computation mode
            - 1: dynamic, based on (Touloukian, 1970; Hobbs, 1974;
                 Travis and Schubert, 2005)
            - 9: constant parameter rhocpfluidm

# Returns
    
        - ρᶠCₚᶠ: volumetric isobaric heat capacity of fluid
"""
function compute_rhocpfluidm(T, mode)
    if mode == 1
        if T <= 0.0
            throw(DomainError(T, "Absolute temperature must be positive"))
        end
        if T < tmfluidphase - 5.0
            ρᶠCₚᶠ = ρH₂Oᶠⁱ * 7.67T
        elseif T < tmfluidphase
            ρᶠCₚᶠ = ρH₂Oᶠⁱ * (7.67T + 0.1Lᶠ)
        elseif T < tmfluidphase + 5.0
            ρᶠCₚᶠ = ρH₂Oᶠ * (4200.0 + 0.1Lᶠ)
        elseif T < 410.0
            ρᶠCₚᶠ = ρH₂Oᶠ * 4200.0
        else
            ρᶠCₚᶠ = ρH₂Oᶠ * (-4.67e4 + 333T - 0.731T^2 + 5.4e-4T^3)
        end
    elseif mode == 9
        @inbounds ρᶠCₚᶠ = rhocpfluidm[1]
    else
        throw(ArgumentError("unknown mode $mode"))
    end
    return ρᶠCₚᶠ
end # function compute_rhocpfluidm

"""
Compute thermal conductivity of silicate (solid phase) based on temperature.

$(SIGNATURES)

# Details
    
        - T: temperature [K]
        - mode:
            - 1: dynamic, based on (Gerya, 2019)
            - 9: constant parameter ksolidm

# Returns
    
        - kᶠ: thermal conductivity of solid
"""
function compute_ksolidm(T, mode)
    if mode == 1
        if T <= 0.0
            throw(DomainError(T, "Absolute temperature must be positive"))
        end
        kˢ = 0.73 + 1293.0 / (T + 77.0)
    elseif mode == 9
        @inbounds kˢ = ksolidm[1]
    else
        throw(ArgumentError("unknown mode $mode"))
    end
    return kˢ
end # function compute_ksolidm

"""
Compute thermal conductivity of H₂O (fluid phase) based on temperature.

$(SIGNATURES)

# Details
    
        - T: temperature [K]
        - mode:
            - 1: dynamic, based on (Touloukian, 1970; Hobbs, 1974;
                 Grimm & Mcsween, 1989; Bland & Travis, 2017)
            - 9: constant parameter kfluidm

# Returns
    
        - kᶠ: thermal conductivity of fluid
"""
function compute_kfluidm(T, mode)
    if mode == 1
        if T <= 0.0
            throw(DomainError(T, "Absolute temperature must be positive"))
        end
        if T < tmfluidphase
            kᶠ = 0.465 + 488.0 / T
        elseif T < 410.0
            kᶠ = -0.581 + 6.34e-3T - 7.93e-6T^2
        else
            kᶠ = -0.142 + 4.12e-3T - 5.01e-6T^2
        end
    elseif mode == 9
        @inbounds kᶠ = kfluidm[1]
    else
        throw(ArgumentError("unknown mode $mode"))
    end
    return kᶠ
end # function compute_kfluidm

"""
Compute dehydration reaction time Δtreaction based on temperature and porosity 
according to selected method:

$(SIGNATURES)

# Details

    - T: temperature
    - ϕ: porosity
    - mode:
        - 1: Gaussian form reaction rate coefficient, based on (Martin & Fyfe,
             1970; Emmanuel & Berkowitz, 2006; Iyer et al., 2012). 
        - 2: pseudo-Arrhenius form reaction rate coefficient, based on
             (Bland & Travis, 2017).
        - 3: Arrhenius form reaction coefficent, based on (Travis et al., 2018).
        - 9: constant parameter Δtreaction
    
# Returns

    - Δtreaction: dehydration reaction time
"""
function compute_Δtreaction(
    T, ϕ, mode; cfg=nothing, is_hydration::Union{Bool,Nothing}=nothing
)
    if mode in (1, 2, 3)
        if T <= 0.0
            throw(DomainError(T, "Absolute temperature must be positive"))
        end
        if !(0.0 < ϕ <= 1.0)
            throw(DomainError(ϕ, "Porosity must be in (0, 1]"))
        end
    end
    A_I_val = cfg === nothing ? A_I : cfg.A_I
    b_I_val = cfg === nothing ? b_I : cfg.b_I
    c_I_val = cfg === nothing ? c_I : cfg.c_I
    Sxo_B_val = cfg === nothing ? Sxo_B : cfg.Sxo_B
    Tscl_B_val = cfg === nothing ? Tscl_B : cfg.Tscl_B
    To_B_val = cfg === nothing ? To_B : cfg.To_B
    dtreaction_val = if cfg === nothing
        Δtreaction
    else
        hyd = is_hydration !== nothing ? is_hydration : (mode == 1)
        hyd ? cfg.dtreaction_hydration : cfg.dtreaction_dehydration
    end

    if mode == 1
        Δtr = -log_completion_rate / (A_I_val * ϕ) * exp(b_I_val * (T - c_I_val)^2)
    elseif mode == 2
        Δtr = -log_completion_rate / (Sxo_B_val * ϕ) * 2.0^((To_B_val - T) / Tscl_B_val)
    elseif mode == 3
        Δtr =
            -log_completion_rate / (Sxo_B_val * ϕ) * exp(Ea_T / RG * (1.0 / T - 1.0 / To_T))
    elseif mode == 9
        Δtr = dtreaction_val
    else
        throw(ArgumentError("unknown mode $mode"))
    end
    return Δtr
end # function compute_dtreaction

"""
Compute molar Gibbs free energy for single dehydration reaction
Wsilicate = Dsilicate + H₂O (16.144)

$(SIGNATURES)

# Details

    - T: temperature
    - pf: fluid pressure
    - XDˢ: molar fraction of dry solid
    - XWˢ: molar fraction of wet solid
    - Δt: timestep size
    - Δtr: total reaction time Δtreaction

# Returns

    - ΔGWD: molar Gibbs free energy for single dehydration reaction (16.165a/b).
"""
@inline function compute_gibbs_free_energy(T, pf, XDˢ, XWˢ, Δt, Δtr; cfg=nothing)
    if T <= 0.0
        throw(DomainError(T, "Absolute temperature must be positive"))
    end
    p_cav = cfg === nothing ? 0.0 : cfg.p_cavitation
    if pf < -p_cav
        throw(DomainError(pf, "Fluid pressure must be >= -$p_cav"))
    end
    if !(0.0 < XDˢ < 1.0) || !(0.0 < XWˢ < 1.0)
        throw(DomainError((XDˢ, XWˢ), "Molar fractions must be strictly in (0, 1)"))
    end
    if Δt < 0.0 || Δtr <= 0.0
        throw(
            DomainError(
                (Δt, Δtr), "Timestep must be non-negative and reaction time positive"
            ),
        )
    end
    dH = cfg === nothing ? ΔHWD : cfg.delta_H
    dS = cfg === nothing ? ΔSWD : cfg.delta_S
    # compute incomplete reaction for short timestep Δt < Δtreaction
    if Δt < Δtr
        # compute ΔG for dehydration reaction (16.145), (16.165b)
        ΔGWD = (dH - T * dS + pf * ΔVWD + RG * T * log(XDˢ / XWˢ)) * (1.0 - Δt / Δtr)
    else
        # Δt ≥ Δtreaction (16.165a)
        ΔGWD = zero(0.0)
    end
    return ΔGWD
end # function compute_gibbs_free_energy

"""
Compute relative enthalpy of system for single dehydration reaction
Wsilicate = Dsilicate + H₂O (16.144).

$(SIGNATURES)

# Details

    - T: temperature
    - pf: fluid pressure
    - XDsolid: molar fraction of dry solid
    - XWsolid: molar fraction of wet solid
    - Δt: timestep size

# Returns

    - Hᵗ: relative enthalpy of system for single dehydration reaction (16.163)
"""
function compute_relative_enthalpy(Xsolid, XWsolid; cfg=nothing, tol=1e-12)
    if !(-tol <= Xsolid <= 1.0 + tol) || !(-tol <= XWsolid <= 1.0 + tol)
        throw(DomainError((Xsolid, XWsolid), "Solid fractions must be in [0, 1]"))
    end
    Xsolid_cl = clamp(Xsolid, 0.0, 1.0)
    XWsolid_cl = clamp(XWsolid, 0.0, 1.0)
    dH = cfg === nothing ? ΔHWD : cfg.delta_H
    return -Xsolid_cl * XWsolid_cl * dH / (MD + MH₂O)
end # function compute_relative_enthalpy

"""
Compute dehydration reaction constant (16.151).

$(SIGNATURES)

# Details

    - T: temperature
    - pf: fluid pressure
    - ΔGWD: Gibbs free energy for dehydration reaction

# Returns

    - KWD: dehydration reaction constant (16.151)
"""
@inline function compute_reaction_constant(T, pf, ΔGWD; cfg=nothing)
    if T <= 0.0
        throw(DomainError(T, "Absolute temperature must be positive"))
    end
    p_cav = cfg === nothing ? 0.0 : cfg.p_cavitation
    if pf < -p_cav
        throw(DomainError(pf, "Fluid pressure must be >= -$p_cav"))
    end
    dH = cfg === nothing ? ΔHWD : cfg.delta_H
    dS = cfg === nothing ? ΔSWD : cfg.delta_S
    # compute reaction constant (16.151)
    return exp(-(dH - T * dS + ΔVWD * pf - ΔGWD) / (RG * T))
end # function compute_reaction_constant

"""
Compute thermodynamic properties at P nodes based on interpolation arrays.

$(SIGNATURES)

# Details

    - DMPSUM: DMP interpolation array
    - DHPSUM: DHP interpolation array
    - WTPSUM: WTP interpolation array
    - DMP: mass transfer term at P nodes
    - DHP: enthalpy transfer/latent heating term at P nodes

# Returns

    - nothing
"""
function compute_thermodynamic_xfer!(
    DMPSUM, DHPSUM, WTPSUM, DMP, DHP, DQPFSUM=nothing, DQPF=nothing
)
    Ny1, Nx1 = size(DMP)
    @inbounds begin
        for j in 1:1:Nx1, i in 1:1:Ny1
            if WTPSUM[i, j] > 0.0
                inv_wt = inv(WTPSUM[i, j])
                DMP[i, j] = DMPSUM[i, j] * inv_wt
                DHP[i, j] = DHPSUM[i, j] * inv_wt
                if DQPF !== nothing && DQPFSUM !== nothing
                    DQPF[i, j] = DQPFSUM[i, j] * inv_wt
                end
            else
                DMP[i, j] = zero(0.0)
                DHP[i, j] = zero(0.0)
                if DQPF !== nothing
                    DQPF[i, j] = zero(0.0)
                end
            end
        end
    end # @inbounds
    return nothing
end # function compute_thermodynamic_xfer!

"""
Perform hydrothermomechanical iterations to time step thermal field at P nodes.

$(SIGNATURES)

# Details

    - DMP: mass transfer term at P nodes
    - DHP: enthalpy transfer/latent heating term at P nodes
    - DMPSUM: interpolation of DMP (mass transfer term) at P nodes
    - DHPSUM: interpolation of DHP (enthalpy transfer term) at P nodes
    - WTPSUM: interpolation weights at P nodes 
    - pf: fluid pressure at P nodes
    - tk2: next temperature at P nodes 
    - tm: type of markers
    - xm: x-coordinate of markers
    - ym: y-coordinate of markers
    - XWˢm₀: previous marker wet silicate (solid) fraction
    - XWˢm: current marker wet silicate (solid) fraction
    - phim: current marker porosity
    - phinewm: next generation marker porosity
    - pfm₀: previous marker fluid pressure
    - marknum: current total number of markers
    - Δt: current time step length
    - timestep: current time step
    - titer: current thermochemical iteration number

# Returns

    - nothing
"""
function perform_thermochemical_reaction!(
    DMP,
    DHP,
    DMPSUM,
    DHPSUM,
    WTPSUM,
    pf,
    tk2,
    tm,
    xm,
    ym,
    XWˢm₀,
    XWˢm,
    phim,
    phinewm,
    pfm₀,
    marknum,
    Δt,
    timestep,
    titer;
    coords=nothing,
    DQPF=nothing,
    DQPFSUM=nothing,
    cfg=nothing,
    backload_step1::Bool=true,
)
    react_cfg = if cfg === nothing
        ReactionConfig()
    elseif cfg isa ReactionConfig
        cfg
    elseif hasproperty(cfg, :reaction)
        cfg.reaction
    else
        ReactionConfig()
    end

    if !react_cfg.active
        DMP .= zero(0.0)
        DHP .= zero(0.0)
        if DQPF !== nothing
            DQPF .= zero(0.0)
        end
        return nothing
    end

    # reset interpolation arrays
    reset_thermochemical_properties!(DMPSUM, DHPSUM, WTPSUM, DQPFSUM)
    xp_val = coords === nothing ? xp : coords.xp
    yp_val = coords === nothing ? yp : coords.yp
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    jmin_p_val = coords === nothing ? jmin_p : coords.jmin_p
    jmax_p_val = coords === nothing ? jmax_p : coords.jmax_p
    imin_p_val = coords === nothing ? imin_p : coords.imin_p
    imax_p_val = coords === nothing ? imax_p : coords.imax_p

    pfcoeff_val = react_cfg.pfcoeff
    p_cav = react_cfg.p_cavitation
    alpha_rel = react_cfg.alpha_relaxation
    eps_comp = 1.0e-4

    # iterate over markers
    @inbounds begin
        for m in 1:1:marknum
            if tm[m] < 3
                # for rocks only
                i, j, weights = fix_weights(
                    xm[m],
                    ym[m],
                    xp_val,
                    yp_val,
                    dx_val,
                    dy_val,
                    jmin_p_val,
                    jmax_p_val,
                    imin_p_val,
                    imax_p_val,
                )
                # interpolate temperature from P nodes
                tknm = dot4(grid_vector(i, j, tk2), weights)
                # interpolate fluid pressure from P nodes with cavitation floor
                pfnm = max(-p_cav, dot4(grid_vector(i, j, pf), weights))
                # factor in previous iteration marker fluid pressure
                if titer > 2
                    pfnm = pfnm * (1.0 - pfcoeff_val) + pfm₀[m] * pfcoeff_val
                end
                # store current marker fluid pressure for next iteration
                pfm₀[m] = pfnm

                # compute bulk composition of solid and fluid system
                if isnan(XWˢm₀[m]) || !(0.0 <= XWˢm₀[m] <= 1.0)
                    throw(
                        DomainError(
                            XWˢm₀[m],
                            "Marker wet solid fraction must be in [0, 1] and non-NaN",
                        ),
                    )
                end
                # clamp previous wet solid molar fraction to interior (prevent singular logs)
                XWˢm₀_cl = clamp(XWˢm₀[m], eps_comp, 1.0 - eps_comp)
                XDˢm₀ = 1.0 - XWˢm₀_cl

                # get fluid molar volume
                VH₂O = ifelse(tknm > tmfluidphase, VH₂Oᶠ, VH₂Oᶠⁱ)
                # compute previous fluid molar fraction (16.164)
                Xᶠ₀ =
                    phim[m] * (XWˢm₀_cl * VWˢ + XDˢm₀ * VDˢ) /
                    ((1.0 - phim[m]) * VH₂O + phim[m] * (XWˢm₀_cl * VWˢ + XDˢm₀ * VDˢ))
                # compute previous equilibrium solid molar fraction (16.150)
                Xˢ₀ = 1.0 - Xᶠ₀
                # compute previous water molar fraction (16.147)
                XH₂Oᵗ = (XWˢm₀_cl * Xˢ₀ + Xᶠ₀) / (1.0 + XWˢm₀_cl * Xˢ₀)
                # compute dry solid molar fraction (16.149)
                XDᵗ = 1.0 - XH₂Oᵗ
                # compute previous solid density (16.161)
                ρˢ₀ = (MD + MH₂O * XWˢm₀_cl) / (VDˢ * XDˢm₀ + VWˢ * XWˢm₀_cl)
                # compute previous fluid density (16.162)
                ρᶠ₀ = ifelse(tknm > tmfluidphase, ρH₂Oᶠ, ρH₂Oᶠⁱ)

                # Continuous equilibrium state
                KWD_eq = compute_reaction_constant(tknm, pfnm, 0.0; cfg=react_cfg)
                XW_eq = inv(KWD_eq + 1.0)

                # Select kinetic rate timescale based on continuous affinity direction
                is_hydration = XW_eq > XWˢm₀_cl
                is_dehydration = XW_eq < XWˢm₀_cl

                skip_rx = false
                if is_hydration
                    if !react_cfg.hydration_active
                        skip_rx = true
                    else
                        mode = react_cfg.hydration_mode
                        Δtr = compute_Δtreaction(
                            tknm, phim[m], mode; cfg=react_cfg, is_hydration=true
                        )
                    end
                elseif is_dehydration
                    if !react_cfg.dehydration_active
                        skip_rx = true
                    else
                        mode = react_cfg.dehydration_mode
                        Δtr = compute_Δtreaction(
                            tknm, phim[m], mode; cfg=react_cfg, is_hydration=false
                        )
                    end
                else
                    Δtr = 1.0e30
                end

                if skip_rx
                    XWˢm₁ = XWˢm₀_cl
                else
                    # Finite-rate relaxation kinetics
                    ΔGWD_star = compute_gibbs_free_energy(
                        tknm, pfnm, XDˢm₀, XWˢm₀_cl, Δt, Δtr; cfg=react_cfg
                    )
                    KWD_star = compute_reaction_constant(
                        tknm, pfnm, ΔGWD_star; cfg=react_cfg
                    )
                    XWˢm₁ = inv(KWD_star + 1.0)
                end

                # Clamp reacted fraction to interior
                XWˢm₁ = clamp(XWˢm₁, eps_comp, 1.0 - eps_comp)

                # Picard under-relaxation on reacted wet silicate fraction across iterations
                if titer > 1
                    XWˢm₁_rel = alpha_rel * XWˢm₁ + (1.0 - alpha_rel) * XWˢm[m]
                else
                    XWˢm₁_rel = XWˢm₁
                end
                XWˢm₁_star = clamp(XWˢm₁_rel, eps_comp, 1.0 - eps_comp)
                XDˢm₁_star = 1.0 - XWˢm₁_star

                # compute reacted total solid molar fraction (16.154)
                Xˢ₁ = clamp(XDᵗ / (1.0 - XDᵗ * XWˢm₁_star), 0.0, 1.0)
                # compute reacted fluid molar fraction (16.155)
                Xᶠ₁ = 1.0 - Xˢ₁

                # only process fluid-bearing rocks
                if 0.0 < Xᶠ₁ < 1.0
                    # compute reacted equilibrium porosity (16.156)
                    ϕ₁ =
                        Xᶠ₁ * VH₂O /
                        (Xᶠ₁ * VH₂O + Xˢ₁ * (XWˢm₁_star * VWˢ + XDˢm₁_star * VDˢ))
                    # under-relax porosity across iterations
                    if titer > 1
                        ϕ_rel = alpha_rel * ϕ₁ + (1.0 - alpha_rel) * phinewm[m]
                    else
                        ϕ_rel = ϕ₁
                    end
                    ϕ_new = clamp(ϕ_rel, phimin, phimax)

                    # compute equilibrium solid density (16.161)
                    ρˢ₁ = (MD + MH₂O * XWˢm₁_star) / (VDˢ * XDˢm₁_star + VWˢ * XWˢm₁_star)
                    # compute equilibrium fluid density (16.162)
                    ρᶠ₁ = ifelse(tknm > tmfluidphase, ρH₂Oᶠ, ρH₂Oᶠⁱ)
                    # compute previous-to-reacted-equilibrium volume ratio (16.106)
                    RV =
                        (ρˢ₁ * (1.0 - ϕ_new) + ρᶠ₁ * ϕ_new) /
                        (ρˢ₀ * (1.0 - phim[m]) + ρᶠ₀ * phim[m])
                    # compute mass transfer rate (16.103)
                    Γmass = (ρˢ₀ * RV * (1.0 - phim[m]) - ρˢ₁ * (1.0 - ϕ_new)) / Δt
                    # compute total mass continuity term (16.112e)
                    ΔMm = (1.0 - RV) / Δt

                    # compute relative enthalpies (16.163)
                    Hᵗ₀ = compute_relative_enthalpy(Xˢ₀, XWˢm₀_cl; cfg=react_cfg)
                    Hᵗ₁ = compute_relative_enthalpy(Xˢ₁, XWˢm₁_star; cfg=react_cfg)
                    ΔHᵗ = Hᵗ₁ - Hᵗ₀

                    # Latent heat transfer term DHP:
                    # Hydration is exothermic (heat source, ΔHm > 0).
                    # Dehydration is endothermic (heat sink, ΔHm < 0).
                    if XWˢm₁_star > XWˢm₀_cl
                        ΔHm = abs(Γmass * ΔHᵗ)
                    elseif XWˢm₁_star < XWˢm₀_cl
                        ΔHm = -abs(Γmass * ΔHᵗ)
                    else
                        ΔHm = zero(0.0)
                    end

                    # Fluid continuity exchange diagnostic DQPF (water mass exchange into pore space)
                    Γwater = Γmass
                    ΔQmᶠ = Γwater / ρᶠ₁

                    # update wet solid molar fraction
                    XWˢm[m] = XWˢm₁_star
                    # update porosity
                    phinewm[m] = ϕ_new

                    # backload properties during first timestep only when requested
                    if backload_step1 && timestep == 1
                        XWˢm₀[m] = XWˢm[m]
                        phim[m] = phinewm[m]
                    end

                    # interpolate terms to P nodes
                    interpolate_add_to_grid!(i, j, weights, ΔMm, DMPSUM)
                    interpolate_add_to_grid!(i, j, weights, ΔHm, DHPSUM)
                    if DQPFSUM !== nothing
                        interpolate_add_to_grid!(i, j, weights, ΔQmᶠ, DQPFSUM)
                    end
                    interpolate_add_to_grid!(i, j, weights, one(1.0), WTPSUM)
                end
            end # if tm[m] < 3
        end # for m=1:1:marknum

        # compute thermodynamic properties at P nodes
        compute_thermodynamic_xfer!(DMPSUM, DHPSUM, WTPSUM, DMP, DHP, DQPFSUM, DQPF)
    end # @inbounds
    @info "min/max mass transfer term" extrema(DMP)
    @info "min/max enthalpy transfer term" extrema(DHP)
    if DQPF !== nothing
        @info "min/max fluid source term" extrema(DQPF)
    end
    return nothing
end # function perform_thermochemical_reaction!

"""
Compute shear heating based on basic (temperature) and P grids.

$(SIGNATURES)

# Details

    - HS: shear heating
    - ETA: viscoplastic viscosity at basic nodes
    - SXY: σ₀xy XY stress at basic nodes
    - ETAP: viscosity at P nodes
    - SXX: normal stress at P nodes
    - RX: ηfluid/Kϕ at Vx nodes
    - RY: ηfluid/Kϕ at Vy nodes
    - qxD: qx-Darcy flux at Vx nodes
    - qyD: qy-Darcy flux at Vy nodes
    - PHI: porosity at P nodes
    - ETAPHI: bulk viscosity at P nodes
    - pr: total pressure at P nodes
    - pf: fluid pressure at P nodes

# Returns

    - nothing
"""
function compute_shear_heating!(
    HS,
    ETA,
    SXY,
    ETAP,
    SXX,
    RX,
    RY,
    qxD,
    qyD,
    PHI,
    ETAPHI,
    pr,
    pf;
    hydrofracture::Bool=false,
    TEN=nothing,
    KX=nothing,
    KY=nothing,
    kappa_frac::Real=1.0e3,
    gamma_frac::Real=1.0,
    k_frac_max::Real=1.0e-9,
    coords=nothing,
)
    Ny1, Nx1 = size(HS)
    Nx = Nx1 - 1
    Ny = Ny1 - 1
    for j in 2:1:Nx, i in 2:1:Ny
        # average SXY⋅EXY
        SXYEXY = 0.25 * sum(grid_vector(i-1, j-1, SXY) .^ 2 ./ grid_vector(i-1, j-1, ETA))
        rx_jm1 = RX[i, j - 1]
        rx_j = RX[i, j]
        ry_im1 = RY[i - 1, j]
        ry_i = RY[i, j]
        if hydrofracture && TEN !== nothing
            Peff_x1 = 0.5 * (pr[i, j - 1] + pr[i, j] - pf[i, j - 1] - pf[i, j])
            sigma_tx1 = 0.5 * (TEN[i, j - 1] + TEN[i - 1, j - 1])
            kphi_x1 = (KX !== nothing) ? KX[i, j - 1] : 0.0
            if kphi_x1 > 0.0
                keff_x1 = compute_hydrofracture_permeability(
                    kphi_x1,
                    Peff_x1,
                    sigma_tx1;
                    active=true,
                    kappa_frac=kappa_frac,
                    gamma=gamma_frac,
                    kmax=k_frac_max,
                )
                rx_jm1 = RX[i, j - 1] * (kphi_x1 / keff_x1)
            else
                fx1 = compute_hydrofracture_factor(
                    Peff_x1, sigma_tx1; kappa_frac=kappa_frac, gamma=gamma_frac
                )
                rx_jm1 = max(RX[i, j - 1] / fx1, 1.0e-5 / k_frac_max)
            end

            Peff_x2 = 0.5 * (pr[i, j] + pr[i, j + 1] - pf[i, j] - pf[i, j + 1])
            sigma_tx2 = 0.5 * (TEN[i, j] + TEN[i - 1, j])
            kphi_x2 = (KX !== nothing) ? KX[i, j] : 0.0
            if kphi_x2 > 0.0
                keff_x2 = compute_hydrofracture_permeability(
                    kphi_x2,
                    Peff_x2,
                    sigma_tx2;
                    active=true,
                    kappa_frac=kappa_frac,
                    gamma=gamma_frac,
                    kmax=k_frac_max,
                )
                rx_j = RX[i, j] * (kphi_x2 / keff_x2)
            else
                fx2 = compute_hydrofracture_factor(
                    Peff_x2, sigma_tx2; kappa_frac=kappa_frac, gamma=gamma_frac
                )
                rx_j = max(RX[i, j] / fx2, 1.0e-5 / k_frac_max)
            end

            Peff_y1 = 0.5 * (pr[i - 1, j] + pr[i, j] - pf[i - 1, j] - pf[i, j])
            sigma_ty1 = 0.5 * (TEN[i - 1, j] + TEN[i - 1, j - 1])
            kphi_y1 = (KY !== nothing) ? KY[i - 1, j] : 0.0
            if kphi_y1 > 0.0
                keff_y1 = compute_hydrofracture_permeability(
                    kphi_y1,
                    Peff_y1,
                    sigma_ty1;
                    active=true,
                    kappa_frac=kappa_frac,
                    gamma=gamma_frac,
                    kmax=k_frac_max,
                )
                ry_im1 = RY[i - 1, j] * (kphi_y1 / keff_y1)
            else
                fy1 = compute_hydrofracture_factor(
                    Peff_y1, sigma_ty1; kappa_frac=kappa_frac, gamma=gamma_frac
                )
                ry_im1 = max(RY[i - 1, j] / fy1, 1.0e-5 / k_frac_max)
            end

            Peff_y2 = 0.5 * (pr[i, j] + pr[i + 1, j] - pf[i, j] - pf[i + 1, j])
            sigma_ty2 = 0.5 * (TEN[i, j] + TEN[i, j - 1])
            kphi_y2 = (KY !== nothing) ? KY[i, j] : 0.0
            if kphi_y2 > 0.0
                keff_y2 = compute_hydrofracture_permeability(
                    kphi_y2,
                    Peff_y2,
                    sigma_ty2;
                    active=true,
                    kappa_frac=kappa_frac,
                    gamma=gamma_frac,
                    kmax=k_frac_max,
                )
                ry_i = RY[i, j] * (kphi_y2 / keff_y2)
            else
                fy2 = compute_hydrofracture_factor(
                    Peff_y2, sigma_ty2; kappa_frac=kappa_frac, gamma=gamma_frac
                )
                ry_i = max(RY[i, j] / fy2, 1.0e-5 / k_frac_max)
            end
        end
        # compute shear heating HS
        @inbounds HS[i, j] = (
            SXX[i, j]^2 / ETAP[i, j] +
            SXYEXY +
            (pr[i, j]-pf[i, j])^2 / (1-PHI[i, j]) / ETAPHI[i, j] +
            0.5 * (rx_jm1*qxD[i, j - 1]^2 + rx_j*qxD[i, j]^2) +
            0.5 * (ry_im1*qyD[i - 1, j]^2 + ry_i*qyD[i, j]^2)
        )
    end
    return nothing
end # function compute_shear_heating!

"""
Compute adiabatic heating based on basic (temperature) and P grids.

$(SIGNATURES)

# Details

    - HA: adiabatic heating at P nodes
    - tk1: previous temperature at P nodes
    - ALPHA: thermal expansion coefficient at P nodes
    - ALPHAF: fluid thermal expansion coefficient at P nodes
    - PHI: porosity at P nodes
    - vx: solid vx-velocity at Vx nodes
    - vy: solid vy-velocity at Vy nodes
    - vxf: fluid vx-velocity at Vx nodes
    - vyf: fluid vy-velocity at Vy nodes
    - ps: solid pressure at P nodes
    - pf: fluid pressure at P nodes

# Returns

    - nothing
"""
function compute_adiabatic_heating!(
    HA, tk1, ALPHA, ALPHAF, PHI, vx, vy, vxf, vyf, ps, pf; coords=nothing
)
    Ny1, Nx1 = size(HA)
    Nx = Nx1 - 1
    Ny = Ny1 - 1
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    @inbounds begin
        for j in 2:1:Nx, i in 2:1:Ny
            # indirect calculation of DP/Dt ≈ (∂P/∂x)⋅vx + (∂P/∂y)⋅vy (eq. 9.23)
            # average vy, vx, vxf, vyf
            VXP = 0.5 * (vx[i, j]+vx[i, j - 1])
            VYP = 0.5 * (vy[i, j]+vy[i - 1, j])
            VXFP = 0.5 * (vxf[i, j]+vxf[i, j - 1])
            VYFP = 0.5 * (vyf[i, j]+vyf[i - 1, j])
            # evaluate DPsolid/Dt with upwind differences
            if VXP > 0.0
                dpsdx = (ps[i, j] - ps[i, j - 1]) * inv(dx_val)
            else
                dpsdx = (ps[i, j + 1] - ps[i, j]) * inv(dx_val)
            end
            if VYP > 0.0
                dpsdy = (ps[i, j] - ps[i - 1, j]) * inv(dy_val)
            else
                dpsdy = (ps[i + 1, j] - ps[i, j]) * inv(dy_val)
            end
            dpsdt = VXP * dpsdx + VYP * dpsdy
            # evaluate DPfluid/Dt with upwind differences
            if VXFP > 0.0
                dpfdx = (pf[i, j]-pf[i, j - 1]) * inv(dx_val)
            else
                dpfdx = (pf[i, j + 1]-pf[i, j]) * inv(dx_val)
            end
            if VYFP > 0.0
                dpfdy = (pf[i, j]-pf[i - 1, j]) * inv(dy_val)
            else
                dpfdy = (pf[i + 1, j]-pf[i, j]) * inv(dy_val)
            end
            dpfdt = VXFP*dpfdx + VYFP*dpfdy
            # Hₐ = (1-ϕ)Tαˢ⋅DPˢ/Dt + ϕTαᶠ⋅DPᶠ/Dt (eq. 9.23)
            HA[i, j] = (
                (1-PHI[i, j]) * tk1[i, j] * ALPHA[i, j] * dpsdt +
                PHI[i, j] * tk1[i, j] * ALPHAF[i, j] * dpfdt
            )
        end
    end # @inbounds
end # function compute_adiabatic_heating!

"""
Compute drained bulk compressibility of a porous medium.

$(SIGNATURES)

# Details

    β_d = (β_ϕ + β_s) / (1 - ϕ)

where β_ϕ is pore compressibility [1/Pa], β_s is solid matrix compressibility [1/Pa],
and ϕ is porosity [-]. References: Biot (1941), Detournay & Cheng (1993), Gerya (2019).

# Arguments

    - betaphi: pore compressibility β_ϕ [1/Pa]
    - phi: porosity ϕ [-]
    - betasolid: solid matrix compressibility β_s [1/Pa]
    - phimin: minimum porosity limit
    - phimax: maximum porosity limit

# Returns

    - betadrained: drained bulk compressibility β_d [1/Pa]
"""
function compute_drained_compressibility(
    betaphi::Real, phi::Real, betasolid::Real; phimin::Real=phimin, phimax::Real=phimax
)
    bphi = max(betaphi, 0.0)
    bsolid = max(betasolid, 0.0)
    phi_eff = clamp(phi, phimin, phimax)
    return (bphi + bsolid) / (1.0 - phi_eff)
end

"""
Compute Biot-Willis coefficient for poroelastic coupling.

$(SIGNATURES)

# Details

    K_BW = 1 - β_s / β_d

For an incompressible solid matrix (β_s = 0), K_BW = 1.
For intact zero-porosity rock (β_d → β_s), K_BW → 0.
Physical bounds: K_BW ∈ [0, 1]. References: Biot (1941), Wang (2000).

# Arguments

    - betadrained: drained bulk compressibility β_d [1/Pa]
    - betasolid: solid matrix compressibility β_s [1/Pa]

# Returns

    - kbw: Biot-Willis coefficient K_BW [-]
"""
function compute_biot_willis_coefficient(betadrained::Real, betasolid::Real)
    if betasolid <= 0.0
        return 1.0
    end
    if betadrained <= betasolid
        return 0.0
    end
    return clamp(1.0 - betasolid / betadrained, 0.0, 1.0)
end

"""
Compute Skempton coefficient B for pore pressure response to mean stress.

$(SIGNATURES)

# Details

    B = (β_d - β_s) / (β_d - β_s + ϕ * (β_f - β_s))

For incompressible constituents (β_s = 0, β_f = 0), B = 1.
Physical bounds: B ∈ [0, 1]. References: Skempton (1954), Rice & Cleary (1976).

# Arguments

    - betadrained: drained bulk compressibility β_d [1/Pa]
    - phi: porosity ϕ [-]
    - betasolid: solid matrix compressibility β_s [1/Pa]
    - betafluid: pore fluid compressibility β_f [1/Pa]
    - phimin: minimum porosity limit
    - phimax: maximum porosity limit

# Returns

    - ksk: Skempton coefficient B [-]
"""
function compute_skempton_coefficient(
    betadrained::Real,
    phi::Real,
    betasolid::Real,
    betafluid::Real;
    phimin::Real=phimin,
    phimax::Real=phimax,
)
    if betasolid <= 0.0 && betafluid <= 0.0
        return 1.0
    end
    bsolid = max(betasolid, 0.0)
    bfluid = max(betafluid, 0.0)
    phi_eff = clamp(phi, phimin, phimax)
    num = betadrained - bsolid
    denom = num + phi_eff * (bfluid - bsolid)
    if denom <= 0.0 || num <= 0.0
        return 1.0
    end
    return clamp(num / denom, 0.0, 1.0)
end

"""
    compute_rhofluid(T::Real, rho0::Real, alpha::Real, T0::Real; thermal_buoyancy::Bool = true)

Compute temperature-dependent pore fluid density with volumetric thermal expansion:
    ρ_f(T) = ρ_{f0} * max(0.1, 1.0 - α_f * (T - T_0))   for T > T_0

When `thermal_buoyancy = false`, returns reference density `rho0` unmodified.
Density is clamped to a lower bound of `0.1 * rho0` to prevent unphysical negative values
at extreme temperatures.

# Arguments

    - T: temperature [K]
    - rho0: reference fluid density at T0 [kg/m³]
    - alpha: fluid volumetric thermal expansion coefficient α_f [1/K]
    - T0: reference temperature [K]
    - thermal_buoyancy: enable or disable thermal expansion (default: true)

# Returns

    - rhof: temperature-dependent fluid density [kg/m³]
"""
function compute_rhofluid(
    T::Real, rho0::Real, alpha::Real, T0::Real; thermal_buoyancy::Bool=true
)
    if !isfinite(T) || !thermal_buoyancy || alpha <= 0.0 || T <= T0
        return Float64(rho0)
    end
    factor = max(0.1, 1.0 - alpha * (T - T0))
    return isnan(factor) ? Float64(rho0) : Float64(rho0 * factor)
end

"""
    compute_fluid_viscosity(T::Real, tm::Integer;
                            mode::Symbol = :arrhenius,
                            eta0::Real = 1.0e-3,
                            eta_ice::Real = 1.0e12,
                            eta_air::Real = 1.0e-3,
                            Ea::Real = 15.0e3,
                            T0::Real = 293.15,
                            tmfluidphase::Real = 273.0,
                            etamin::Real = 1.0e-5,
                            etamax::Real = 1.0e12)

Compute temperature-dependent dynamic fluid viscosity η_f(T) [Pa s].

For sticky air markers (`tm >= 3`), returns `eta_air`.
For sub-freezing rock markers (`T <= tmfluidphase`), returns `eta_ice`.
For non-finite or corrupt temperatures (`!isfinite(T)`), returns `eta_ice` to prevent runaway mobility.
For liquid fluid markers (`T > tmfluidphase`):
  - `:constant` mode: returns `eta0`.
  - `:arrhenius` mode:
      η_f(T) = η_0 * exp((E_a / R) * (1/T - 1/T_0))
    clamped to [etamin, etamax].

# Arguments

    - T: temperature [K]
    - tm: material type index (1: core, 2: crust, 3: sticky air)
    - mode: `:arrhenius` or `:constant` (default: `:arrhenius`)
    - eta0: reference liquid fluid viscosity at T0 [Pa s]
    - eta_ice: sub-freezing ice viscosity [Pa s]
    - eta_air: sticky air fluid viscosity [Pa s]
    - Ea: activation energy for fluid viscous flow [J/mol]
    - T0: reference temperature [K]
    - tmfluidphase: melting temperature [K]
    - etamin: minimum viscosity floor [Pa s]
    - etamax: maximum viscosity ceiling [Pa s]

# Returns

    - etafluid: dynamic fluid viscosity [Pa s]
"""
function compute_fluid_viscosity(
    T::Real,
    tm::Integer;
    mode::Symbol=:arrhenius,
    eta0::Real=1.0e-3,
    eta_ice::Real=1.0e12,
    eta_air::Real=1.0e-3,
    Ea::Real=15.0e3,
    T0::Real=293.15,
    tmfluidphase::Real=273.0,
    etamin::Real=1.0e-5,
    etamax::Real=1.0e12,
)
    if tm >= 3
        return Float64(eta_air)
    end
    if !isfinite(T)
        return Float64(eta_ice)
    end
    if T <= tmfluidphase
        return Float64(eta_ice)
    end
    if mode === :constant || Ea <= 0.0
        return Float64(eta0)
    elseif mode === :arrhenius
        # Universal gas constant R [J/(mol K)]
        R_gas = 8.31446261815324
        log_ratio = (Ea / R_gas) * (inv(T) - inv(T0))
        val = eta0 * exp(log_ratio)
        return clamp(val, Float64(etamin), Float64(etamax))
    else
        throw(
            ArgumentError(
                "Unknown fluid viscosity mode: $mode (expected :arrhenius or :constant)"
            ),
        )
    end
end

"""
    compute_hydrofracture_factor(Peff::Real, sigma_t::Real;
                                 active::Bool = true,
                                 kappa_frac::Real = 1.0e3,
                                 gamma::Real = 1.0,
                                 max_factor::Real = Inf)

Compute dimensionless permeability enhancement factor from dynamic hydrofracturing.

When pore fluid pressure exceeds total confining pressure plus tensile strength
(Terzaghi effective pressure Peff = Pt - Pf <= -sigma_t), hydraulic tensile
fractures open and increase effective permeability:

    factor = 1.0 + kappa_frac * ((-Peff - sigma_t) / sigma_t)^gamma

clamped to [1.0, max_factor].
"""
function compute_hydrofracture_factor(
    Peff::Real,
    sigma_t::Real;
    active::Bool=true,
    kappa_frac::Real=1.0e3,
    gamma::Real=1.0,
    max_factor::Real=Inf,
)
    if !active || !isfinite(Peff) || !isfinite(sigma_t) || sigma_t <= 0.0
        return 1.0
    end
    overpressure = -Peff - sigma_t
    if overpressure <= 0.0
        return 1.0
    end
    norm_overpressure = overpressure / sigma_t
    factor = 1.0 + kappa_frac * (norm_overpressure ^ gamma)
    return clamp(Float64(factor), 1.0, Float64(max_factor))
end

"""
    compute_hydrofracture_permeability(kphi::Real, Peff::Real, sigma_t::Real;
                                       active::Bool = true,
                                       kappa_frac::Real = 1.0e3,
                                       gamma::Real = 1.0,
                                       kmax::Real = 1.0e-9)

Compute effective permeability k_eff [m²] with dynamic hydrofracturing enhancement.

When pore fluid pressure exceeds total confining pressure plus rock tensile strength:

    Peff = Pt - Pf <= -sigma_t

tensile microcracks open and enhance permeability according to:

    k_eff = min(kphi * compute_hydrofracture_factor(Peff, sigma_t; active, kappa_frac, gamma), kmax)

# Arguments
- `kphi`: baseline matrix permeability [m²]
- `Peff`: Terzaghi effective pressure Pt - Pf [Pa]
- `sigma_t`: rock tensile strength [Pa]
- `active`: enable hydrofracture enhancement (default: true)
- `kappa_frac`: dimensionless enhancement multiplier (default: 1.0e3)
- `gamma`: power-law exponent (default: 1.0)
- `kmax`: maximum permeability ceiling [m²] (default: 1.0e-9)

# Returns
- `k_eff`: effective permeability [m²]
"""
function compute_hydrofracture_permeability(
    kphi::Real,
    Peff::Real,
    sigma_t::Real;
    active::Bool=true,
    kappa_frac::Real=1.0e3,
    gamma::Real=1.0,
    kmax::Real=1.0e-9,
)
    if !active || !isfinite(Peff) || !isfinite(sigma_t) || sigma_t <= 0.0 || kphi <= 0.0
        return Float64(kphi)
    end
    factor = compute_hydrofracture_factor(
        Peff, sigma_t; active=active, kappa_frac=kappa_frac, gamma=gamma
    )
    k_enhanced = kphi * factor
    return clamp(Float64(k_enhanced), Float64(kphi), Float64(kmax))
end

"""
    compute_disk_temperature(t_seconds::Real, cfg::DiskConfig)::Float64

Compute ambient protoplanetary disk midplane temperature T_disk [K] at time t [s].

Supports three operating modes in `cfg.model`:
- `:fixed`: Constant ambient temperature `cfg.t_ambient`.
- `:monotonic`: Power-law viscous clearing (Lynden-Bell & Pringle 1974;
  Johansen et al. 2015).
- `:class1_to_class2`: Cold-to-hot-to-cold evolution from disk buildup
  through peak accretion to flared irradiation floor
  (Drążkowska & Dullemond 2018; Lichtenberg et al. 2021; Williams et al. 2026;
  `:class0_to_class2` supported as alias).

# Arguments
- `t_seconds`: Simulation time [s]
- `cfg`: Protoplanetary disk configuration struct (`DiskConfig`)
- `orbital_distance_au`: Optional orbital distance override [AU]
  (defaults to `cfg.orbital_distance_au`)
- `stellar_mass_msun`: Optional host star mass override [M_sun]
  (defaults to `cfg.stellar_mass_msun`)

# Returns
- `T_disk`: Ambient disk temperature [K]
"""
function compute_disk_temperature(
    t_seconds::Real,
    cfg::DiskConfig;
    orbital_distance_au::Real=cfg.orbital_distance_au,
    stellar_mass_msun::Real=cfg.stellar_mass_msun,
)::Float64
    if !cfg.enabled || cfg.model === :fixed
        return Float64(cfg.t_ambient)
    end

    t_sec_nonneg = max(0.0, Float64(t_seconds))
    t_Myr = t_sec_nonneg / (1.0e6 * (365.25 * 86400.0))
    r_au = Float64(orbital_distance_au)
    m_star = Float64(stellar_mass_msun)

    T_irr = cfg.t_irr_1au * (m_star ^ cfg.p_m_irr) * (r_au ^ (-cfg.q_irr))
    T_peak = cfg.t_peak_1au * (m_star ^ cfg.p_m_visc) * (r_au ^ (-cfg.q_visc))

    T_visc_excess4 = max(0.0, T_peak^4 - T_irr^4)

    if cfg.model === :monotonic
        t_visc = cfg.t_visc_0_myr * (m_star ^ cfg.p_m_visc_decay)
        decay = (1.0 + t_Myr / t_visc) ^ (-cfg.gamma)
        T4 = T_irr^4 + T_visc_excess4 * decay
        return max(cfg.t_cloud, T4 ^ 0.25)
    elseif cfg.model === :class1_to_class2 || cfg.model === :class0_to_class2
        t_peak = cfg.t_peak_time_1au_myr * (m_star ^ cfg.p_m_t) * (r_au ^ cfg.p_r_t)
        tau_star = 0.8 * t_peak

        if t_Myr <= 0.0
            f_acc = 0.0
            g_star = 0.0
        else
            x = t_Myr / t_peak
            f_acc =
                (1.0 + cfg.alpha / cfg.gamma) * (x ^ cfg.alpha) /
                (1.0 + (cfg.alpha / cfg.gamma) * (x ^ (cfg.alpha + cfg.gamma)))
            g_star = 1.0 - exp(-t_Myr / tau_star)
        end

        T_eff_irr4 = cfg.t_cloud^4 + (T_irr^4 - cfg.t_cloud^4) * g_star
        T4 = T_eff_irr4 + T_visc_excess4 * f_acc
        return max(cfg.t_cloud, T4 ^ 0.25)
    else
        throw(ArgumentError("Unknown disk model: $(cfg.model)"))
    end
end

"""
    compute_snowline_radius(
        t_seconds::Real, cfg::DiskConfig;
        T_sub::Real=170.0, r_min::Real=0.05, r_max::Real=100.0,
        tol::Real=1e-4, max_iter::Int=50,
        stellar_mass_msun::Real=cfg.stellar_mass_msun
    )::Float64

Compute heliocentric water snowline radius [AU] at time `t_seconds` where ambient
disk temperature equals volatile sublimation temperature `T_sub` (default: 170.0 K).
Returns `r_min` if the entire disk is below `T_sub`, or `r_max` if the disk remains
above `T_sub`. Employs a coarse-to-fine radial scan from `r_max` inward to guarantee
locating the outermost snowline under non-monotonic profiles.

# Arguments
- `t_seconds`: Simulation time [s]
- `cfg`: Protoplanetary disk configuration struct (`DiskConfig`)
- `T_sub`: Volatile sublimation threshold temperature [K] (default: 170.0 K)
- `r_min`: Minimum search radius [AU] (default: 0.05 AU)
- `r_max`: Maximum search radius [AU] (default: 100.0 AU)
- `tol`: Absolute convergence tolerance in orbital distance [AU] (default: 1e-4 AU)
- `max_iter`: Maximum bisection iterations (default: 50)
- `stellar_mass_msun`: Host star mass [M_sun] (default: `cfg.stellar_mass_msun`)

# Returns
- `r_snow`: Water snowline orbital distance [AU]
"""
function compute_snowline_radius(
    t_seconds::Real,
    cfg::DiskConfig;
    T_sub::Real=170.0,
    r_min::Real=0.05,
    r_max::Real=100.0,
    tol::Real=1e-4,
    max_iter::Int=50,
    stellar_mass_msun::Real=cfg.stellar_mass_msun,
)::Float64
    T_outer = compute_disk_temperature(
        t_seconds, cfg; orbital_distance_au=r_max, stellar_mass_msun=stellar_mass_msun
    )
    if T_outer >= T_sub
        return Float64(r_max)
    end

    # Scan radially inward from r_max to r_min across log-spaced intervals
    # to locate the outermost crossing bracket [r_lo, r_hi]
    n_coarse = 200
    log_rmax = log(Float64(r_max))
    log_rmin = log(Float64(r_min))
    dlog_r = (log_rmax - log_rmin) / n_coarse

    r_hi = Float64(r_max)
    r_lo = Float64(r_min)
    found_bracket = false

    for step in 1:n_coarse
        r_step = exp(log_rmax - step * dlog_r)
        T_step = compute_disk_temperature(
            t_seconds, cfg; orbital_distance_au=r_step, stellar_mass_msun=stellar_mass_msun
        )
        if T_step >= T_sub
            r_lo = r_step
            r_hi = exp(log_rmax - (step - 1) * dlog_r)
            found_bracket = true
            break
        end
    end

    if !found_bracket
        return Float64(r_min)
    end

    for _ in 1:max_iter
        r_mid = 0.5 * (r_lo + r_hi)
        T_mid = compute_disk_temperature(
            t_seconds, cfg; orbital_distance_au=r_mid, stellar_mass_msun=stellar_mass_msun
        )
        if abs(T_mid - T_sub) < tol || (r_hi - r_lo) < tol
            return r_mid
        end
        if T_mid >= T_sub
            r_lo = r_mid
        else
            r_hi = r_mid
        end
    end
    return 0.5 * (r_lo + r_hi)
end

"""
    compute_radiation_htc(
        T_surf::Real, T_amb::Real;
        emissivity::Real=0.9, sigma_sb::Real=5.670374419e-8
    )::Float64

Compute linearized Stefan-Boltzmann radiative heat transfer coefficient h_rad [W/(m² K)]:

    h_rad = ε * σ_SB * (T_surf² + T_amb²) * (T_surf + T_amb)

such that h_rad * (T_surf - T_amb) = ε * σ_SB * (T_surf⁴ - T_amb⁴).

# Arguments
- `T_surf`: Surface temperature [K]
- `T_amb`: Ambient disk temperature [K]
- `emissivity`: Surface thermal emissivity in [0, 1]
- `sigma_sb`: Stefan-Boltzmann constant [W/(m² K⁴)]

# Returns
- `h_rad`: Linearized radiative heat transfer coefficient [W/(m² K)]
"""
function compute_radiation_htc(
    T_surf::Real, T_amb::Real; emissivity::Real=0.9, sigma_sb::Real=5.670374419e-8
)::Float64
    if !(0.0 <= emissivity <= 1.0)
        throw(DomainError(emissivity, "Emissivity must be in [0.0, 1.0]"))
    end
    if !isfinite(T_surf) || !isfinite(T_amb) || T_surf <= 0.0 || T_amb <= 0.0
        return 0.0
    end
    T_s = Float64(T_surf)
    T_a = Float64(T_amb)
    return Float64(emissivity) * Float64(sigma_sb) * (T_s^2 + T_a^2) * (T_s + T_a)
end

"""
    compute_disk_dispersal_weight(time_seconds::Real; t_dispersal_myr::Real=3.0, dt_dispersal_myr::Real=0.1)::Float64

Compute smooth sigmoid transition weight w_disp in [0, 1] representing the fraction
of circumstellar gas disk cleared at time `time_seconds`:

    w_disp = 1 / (1 + exp(-(t_Myr - t_dispersal_myr) / dt_dispersal_myr))

# Arguments
- `time_seconds`: Simulation time [s]
- `t_dispersal_myr`: Epoch of disk gas dispersal [Myr] (default: 3.0 Myr)
- `dt_dispersal_myr`: Characteristic duration of dispersal transition [Myr] (default: 0.1 Myr)

# Returns
- `w_disp`: Sigmoid dispersal weight in [0, 1]
"""
function compute_disk_dispersal_weight(
    time_seconds::Real; t_dispersal_myr::Real=3.0, dt_dispersal_myr::Real=0.1
)::Float64
    if dt_dispersal_myr <= 0.0 || !isfinite(dt_dispersal_myr)
        throw(DomainError(dt_dispersal_myr, "dt_dispersal_myr must be > 0 and finite"))
    end
    t_sec_nonneg = max(0.0, Float64(time_seconds))
    t_Myr = t_sec_nonneg / (1.0e6 * (365.25 * 86400.0))
    t_disp = Float64(t_dispersal_myr)
    dt_disp = Float64(dt_dispersal_myr)
    arg = clamp((t_Myr - t_disp) / dt_disp, -100.0, 100.0)
    w = 1.0 / (1.0 + exp(-arg))
    return clamp(w, 0.0, 1.0)
end

"""
    compute_solar_equilibrium_temperature(orbital_distance_au::Real; albedo::Real=0.06, stellar_luminosity_lsun::Real=1.0)::Float64

Compute vacuum solar radiation equilibrium temperature T_eq [K] at heliocentric distance
`orbital_distance_au` [AU] assuming fast planetary rotation or uniform spherical emission:

    T_eq = ((1 - A) * L_star / (16 * π * σ_SB * d²))^(1/4)

# Arguments
- `orbital_distance_au`: Heliocentric orbital distance [AU]
- `albedo`: Bond albedo in [0, 1) (default: 0.06 for dark carbonaceous planetesimals)
- `stellar_luminosity_lsun`: Host star luminosity in solar units [L_sun] (default: 1.0)

# Returns
- `T_eq`: Solar radiation equilibrium temperature [K]
"""
function compute_solar_equilibrium_temperature(
    orbital_distance_au::Real; albedo::Real=0.06, stellar_luminosity_lsun::Real=1.0
)::Float64
    r_au = Float64(orbital_distance_au)
    if r_au <= 0.0 || !isfinite(r_au)
        throw(DomainError(r_au, "orbital_distance_au must be > 0 and finite"))
    end
    A = Float64(albedo)
    if !(0.0 <= A < 1.0) || !isfinite(A)
        throw(DomainError(A, "albedo must be in [0.0, 1.0)"))
    end
    L_sun = 3.828e26 * Float64(stellar_luminosity_lsun)
    sigma_sb = 5.670374419e-8
    d_m = r_au * 1.495978707e11
    F_sun = L_sun / (4.0 * π * d_m^2)
    T_eq4 = (1.0 - A) * F_sun / (4.0 * sigma_sb)
    return T_eq4^0.25
end

"""
    compute_ambient_conditions(time_seconds::Real, cfg::DiskConfig)::Tuple{Float64,Float64,Float64}

Compute evolving ambient temperature T_amb [K], ambient pressure P_amb [Pa], and disk
dispersal weight w_disp at time `time_seconds`. Transitions smoothly from nebular disk
conditions to solar radiative equilibrium and space vacuum upon disk gas clearing.

# Arguments
- `time_seconds`: Simulation time [s]
- `cfg`: Protoplanetary disk configuration struct (`DiskConfig`)

# Returns
- `(T_amb, P_amb, w_disp)`: Ambient temperature [K], ambient pressure [Pa], and dispersal weight in [0, 1]
"""
function compute_ambient_conditions(
    time_seconds::Real, cfg::DiskConfig
)::Tuple{Float64,Float64,Float64}
    w_disp = if cfg.dispersal_active
        compute_disk_dispersal_weight(
            time_seconds;
            t_dispersal_myr=cfg.t_dispersal_myr,
            dt_dispersal_myr=cfg.dt_dispersal_myr,
        )
    else
        0.0
    end
    T_disk = compute_disk_temperature(time_seconds, cfg)
    T_eq = if isfinite(cfg.t_eq_custom) && cfg.t_eq_custom > 0.0
        Float64(cfg.t_eq_custom)
    else
        compute_solar_equilibrium_temperature(cfg.orbital_distance_au; albedo=cfg.albedo)
    end
    T_amb = (1.0 - w_disp) * T_disk + w_disp * T_eq
    P_amb = (1.0 - w_disp) * cfg.p_amb_disk + w_disp * cfg.p_amb_space
    return (T_amb, P_amb, w_disp)
end

"""
    compute_ice_vapor_pressure(T::Real; P0::Real=611.66, T0::Real=273.16, L_sub::Real=2.83e6, Rv::Real=461.5)::Float64

Compute water ice sublimation equilibrium vapor pressure P_sat,ice [Pa] at temperature `T` [K]
using the integrated Clausius-Clapeyron relation:

    P_sat,ice = P0 * exp(-(L_sub / Rv) * (1/T - 1/T0))

anchored at the water triple point (T0 = 273.16 K, P0 = 611.66 Pa).

# Arguments
- `T`: Temperature [K]
- `P0`: Triple-point water vapor pressure [Pa] (default: 611.66 Pa)
- `T0`: Triple-point temperature [K] (default: 273.16 K)
- `L_sub`: Latent heat of ice sublimation [J/kg] (default: 2.83e6 J/kg)
- `Rv`: Specific gas constant for water vapor [J/(kg K)] (default: 461.5 J/(kg K))

# Returns
- `P_sat`: Equilibrium ice sublimation vapor pressure [Pa]

# Notes
- For temperatures below the triple point (`T < T0`), vapor pressure follows ice sublimation
  via the integrated Clausius-Clapeyron relation with latent heat `L_sub = 2.83e6 J/kg`.
- Above the triple point (`T >= T0`), saturation vapor pressure over liquid water follows
  the Arden Buck (1981) formulation up to the critical point (`T_crit = 647.096 K`), clamped
  at water critical pressure `P_crit = 22.064 MPa`.
"""
function compute_water_vapor_pressure(
    T::Real; P0::Real=611.66, T0::Real=273.16, L_sub::Real=2.83e6, Rv::Real=461.5
)::Float64
    T_val = Float64(T)
    if T_val <= 0.0 || !isfinite(T_val)
        throw(DomainError(T_val, "Temperature must be > 0 and finite"))
    end
    P0_val = Float64(P0)
    if P0_val <= 0.0 || !isfinite(P0_val)
        throw(DomainError(P0_val, "Triple-point pressure P0 must be > 0 and finite"))
    end
    T0_val = Float64(T0)
    if T0_val <= 0.0 || !isfinite(T0_val)
        throw(DomainError(T0_val, "Triple-point temperature T0 must be > 0 and finite"))
    end
    Rv_val = Float64(Rv)
    if Rv_val <= 0.0 || !isfinite(Rv_val)
        throw(DomainError(Rv_val, "Gas constant Rv must be > 0 and finite"))
    end
    L_sub_val = Float64(L_sub)

    if T_val <= T0_val
        return P0_val * exp(-(L_sub_val / Rv_val) * (1.0 / T_val - 1.0 / T0_val))
    end

    # Liquid water saturation vapor pressure via Arden Buck (1981)
    T_crit = 647.096
    P_crit = 22.064e6
    if T_val >= T_crit
        return P_crit
    end
    Tc = T_val - 273.15
    p_buck = 611.21 * exp((18.678 - Tc / 234.5) * (Tc / (Tc + 257.14)))
    return min(P_crit, max(P0_val, p_buck))
end

const compute_ice_vapor_pressure = compute_water_vapor_pressure

"""
Compute equilibrium saturation vapor pressure [Pa] for volatile species at temperature T_K.

Supported species:
- `:H2O`: Water ice / liquid sublimation & vapor pressure (Clausius-Clapeyron)
- `:CO2`: Carbon dioxide sublimation vapor pressure
- `:CH4`: Methane vapor pressure
- `:CO`: Carbon monoxide vapor pressure
- `:N2`: Molecular nitrogen vapor pressure
- `:H2`: Molecular hydrogen (hyper-volatile / supercritical at T >= 33 K)
- `:NH3`: Ammonia vapor pressure
- `:H2S`: Hydrogen sulfide vapor pressure
- `:SO2`: Sulfur dioxide vapor pressure
- `:S2`: Diatomic sulfur vapor pressure
"""
function compute_species_vapor_pressure(species::Symbol, T_K::Real)::Float64
    T_val = Float64(T_K)
    if T_val <= 0.0 || !isfinite(T_val)
        throw(DomainError(T_val, "Temperature must be > 0 and finite"))
    end
    if species === :H2O
        return compute_water_vapor_pressure(T_val)
    elseif species === :CO2
        # Triple point: T0 = 216.58 K, P0 = 5.18e5 Pa, L_sub = 5.71e5 J/kg, Rv = 188.92 J/(kg K)
        arg = clamp(-(5.71e5 / 188.92) * (1.0 / T_val - 1.0 / 216.58), -100.0, 100.0)
        return 5.18e5 * exp(arg)
    elseif species === :CH4
        # Non-condensible / supercritical above Tc = 190.6 K; triple point: T0 = 90.69 K, P0 = 1.17e4 Pa
        if T_val >= 190.6
            return 0.0
        end
        arg = clamp(-(5.10e5 / 518.3) * (1.0 / T_val - 1.0 / 90.69), -100.0, 100.0)
        return 1.17e4 * exp(arg)
    elseif species === :CO
        # Non-condensible / supercritical above Tc = 132.9 K; triple point: T0 = 68.15 K, P0 = 1.54e4 Pa
        if T_val >= 132.9
            return 0.0
        end
        arg = clamp(-(2.97e5 / 296.8) * (1.0 / T_val - 1.0 / 68.15), -100.0, 100.0)
        return 1.54e4 * exp(arg)
    elseif species === :N2
        # Non-condensible / supercritical above Tc = 126.2 K; triple point: T0 = 63.15 K, P0 = 1.25e4 Pa
        if T_val >= 126.2
            return 0.0
        end
        arg = clamp(-(2.54e5 / 296.8) * (1.0 / T_val - 1.0 / 63.15), -100.0, 100.0)
        return 1.25e4 * exp(arg)
    elseif species === :H2
        # Non-condensible / supercritical above Tc = 33.1 K; triple point: T0 = 13.8 K, P0 = 7.04e3 Pa
        if T_val >= 33.0
            return 0.0
        end
        arg = clamp(-(4.54e5 / 4124.0) * (1.0 / T_val - 1.0 / 13.8), -100.0, 100.0)
        return 7.04e3 * exp(arg)
    elseif species === :NH3
        # Triple point: T0 = 195.4 K, P0 = 6060.0 Pa, L_sub = 1.70e6 J/kg, Rv = 488.2 J/(kg K)
        arg = clamp(-(1.70e6 / 488.2) * (1.0 / T_val - 1.0 / 195.4), -100.0, 100.0)
        return 6060.0 * exp(arg)
    elseif species === :H2S
        # Triple point: T0 = 187.6 K, P0 = 2.32e4 Pa, L_sub = 6.98e5 J/kg, Rv = 243.9 J/(kg K)
        arg = clamp(-(6.98e5 / 243.9) * (1.0 / T_val - 1.0 / 187.6), -100.0, 100.0)
        return 2.32e4 * exp(arg)
    elseif species === :SO2
        # Triple point: T0 = 197.7 K, P0 = 1670.0 Pa, L_sub = 5.25e5 J/kg, Rv = 129.8 J/(kg K)
        arg = clamp(-(5.25e5 / 129.8) * (1.0 / T_val - 1.0 / 197.7), -100.0, 100.0)
        return 1670.0 * exp(arg)
    elseif species === :S2
        # Boiling / reference: T0 = 718.0 K, P0 = 1.0e5 Pa, L_sub = 1.45e6 J/kg, Rv = 129.6 J/(kg K)
        arg = clamp(-(1.45e6 / 129.6) * (1.0 / T_val - 1.0 / 718.0), -100.0, 100.0)
        return 1.0e5 * exp(arg)
    else
        throw(ArgumentError("Unknown volatile species: $species"))
    end
end

"""
    compute_venting_pressure(T_surf::Real, P_amb::Real; species::Symbol=:H2O, P0::Real=611.66, T0::Real=273.16, L_sub::Real=2.83e6, Rv::Real=461.5)::Float64

Compute effective boundary venting fluid pressure P_vent [Pa] at a planetesimal surface:

    P_vent = max(P_amb, P_sat(T_surf))

Enforces the physical cold-trap and boiling constraints: if ambient nebular gas pressure exceeds
saturation vapor pressure at surface temperatures, ambient gas confines pore fluid;
if ambient pressure drops below saturation vapor pressure, boiling or flash sublimation
sets the effective boundary venting pressure.

# Arguments
- `T_surf`: Planetesimal surface temperature [K]
- `P_amb`: Ambient surrounding gas pressure [Pa]

# Keyword Arguments
- `species`: Volatile species identifier (default: `:H2O`)

# Returns
- `P_vent`: Effective venting boundary pressure [Pa]
"""
function compute_venting_pressure(
    T_surf::Real,
    P_amb::Real;
    species::Symbol=:H2O,
    P0::Real=611.66,
    T0::Real=273.16,
    L_sub::Real=2.83e6,
    Rv::Real=461.5,
)::Float64
    P_sat = if species === :H2O
        compute_water_vapor_pressure(T_surf; P0=P0, T0=T0, L_sub=L_sub, Rv=Rv)
    else
        compute_species_vapor_pressure(species, T_surf)
    end
    return max(Float64(P_amb), P_sat)
end

"""
    compute_ice_sealed_permeability(
        k0::Real, T::Real;
        T_freeze::Real=273.15, delta_T_seal::Real=10.0, k_min_ratio::Real=1.0e-6
    )::Float64

Compute effective rock permeability k_sealed [m²] reduced by cryogenic pore ice freezing:

    k_sealed = k0 * [ (1 - r_min) * exp(-(T_freeze - T) / ΔT_seal) + r_min ]

For temperatures at or above freezing (`T >= T_freeze`), pore ice melts and permeability
equals `k0`. For sub-freezing temperatures (`T < T_freeze`), pore ice blocks pore throats
and reduces permeability exponentially toward the floor ratio `r_min = k_min_ratio`.

# Arguments
- `k0`: Reference unsealed permeability [m²]
- `T`: Local rock temperature [K]
- `T_freeze`: Water freezing temperature [K] (default: 273.15 K)
- `delta_T_seal`: Temperature scale for ice sealing [K] (default: 10.0 K)
- `k_min_ratio`: Minimum residual permeability floor ratio in (0, 1] (default: 1.0e-6)

# Returns
- `k_sealed`: Effective sealed permeability [m²]
"""
function compute_ice_sealed_permeability(
    k0::Real,
    T::Real;
    T_freeze::Real=273.15,
    delta_T_seal::Real=10.0,
    k_min_ratio::Real=1.0e-6,
)::Float64
    k0_val = Float64(k0)
    if k0_val <= 0.0 || !isfinite(k0_val)
        throw(DomainError(k0_val, "Reference permeability k0 must be > 0 and finite"))
    end
    T_val = Float64(T)
    if T_val <= 0.0 || !isfinite(T_val)
        throw(DomainError(T_val, "Temperature T must be > 0 and finite"))
    end
    T_frz = Float64(T_freeze)
    if T_frz <= 0.0 || !isfinite(T_frz)
        throw(DomainError(T_frz, "T_freeze must be > 0 and finite"))
    end
    dT_seal = Float64(delta_T_seal)
    if dT_seal <= 0.0 || !isfinite(dT_seal)
        throw(DomainError(dT_seal, "delta_T_seal must be > 0 and finite"))
    end
    r_min = Float64(k_min_ratio)
    if !(0.0 < r_min <= 1.0) || !isfinite(r_min)
        throw(DomainError(r_min, "k_min_ratio must be in (0, 1] and finite"))
    end

    if T_val >= T_frz
        return k0_val
    end

    arg = (T_frz - T_val) / dT_seal
    factor = (1.0 - r_min) * exp(-arg) + r_min
    return clamp(k0_val * factor, k0_val * r_min, k0_val)
end

"""
    is_hydrofracture_breached(Peff::Real, sigma_t::Real)::Bool
    is_hydrofracture_breached(Pt::Real, Pf::Real, sigma_t::Real)::Bool

Assess whether hydraulic tensile failure criterion (Peff <= -sigma_t) is satisfied,
breaching the rock matrix or cryogenic ice lid.

# Arguments
- `Peff`: Terzaghi effective pressure Pt - Pf [Pa]
- `Pt`: Total confining pressure [Pa]
- `Pf`: Pore fluid pressure [Pa]
- `sigma_t`: Rock tensile strength [Pa]

# Returns
- `breached::Bool`: `true` if hydraulic tensile fractures open, `false` otherwise.
  Non-finite inputs or non-positive `sigma_t <= 0.0` return `false` as an invalid-input safety guard.
"""
function is_hydrofracture_breached(Peff::Real, sigma_t::Real)::Bool
    if !isfinite(Peff) || !isfinite(sigma_t) || sigma_t <= 0.0
        return false
    end
    return Float64(Peff) <= -Float64(sigma_t)
end

function is_hydrofracture_breached(Pt::Real, Pf::Real, sigma_t::Real)::Bool
    return is_hydrofracture_breached(Pt - Pf, sigma_t)
end

"""
    compute_spherical_metric_heat_source!(Q_metric, tk, KX, KY, coords;
                                          xcenter, ycenter, rplanet, reg_cells=0.5)

Compute geometric metric volumetric heat source Q_metric [W/m³] on P-nodes
to account for 3D spherical divergence on a 2D Cartesian grid:

    Q_metric = (k / r_eff²) * [ (x - xc) * (∂T/∂x) + (y - yc) * (∂T/∂y) ]

inside r <= rplanet, and 0 in sticky air.

# Arguments
- `Q_metric`: Output matrix [W/m³] of size (Ny1, Nx1)
- `tk`: Temperature field [K] of size (Ny1, Nx1)
- `KX`: Thermal conductivity at Vx nodes [W/(m K)]
- `KY`: Thermal conductivity at Vy nodes [W/(m K)]
- `coords`: Grid coordinate descriptors
- `xcenter`: Horizontal center coordinate [m]
- `ycenter`: Vertical center coordinate [m]
- `rplanet`: Planetesimal radius [m]
- `reg_cells`: Radial regularization parameter in grid cell units
"""
function compute_spherical_metric_heat_source!(
    Q_metric::AbstractMatrix{Float64},
    tk::AbstractMatrix{Float64},
    KX::AbstractMatrix{Float64},
    KY::AbstractMatrix{Float64},
    coords::GridCoordinates;
    xcenter::Real,
    ycenter::Real,
    rplanet::Real,
    reg_cells::Real=0.5,
)
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    dx = coords.dx
    dy = coords.dy
    inv_2dx = inv(2.0 * dx)
    inv_2dy = inv(2.0 * dy)
    eps_r = reg_cells * min(dx, dy)
    eps_r2 = eps_r^2
    rplanet2 = rplanet^2

    fill!(Q_metric, 0.0)

    @inbounds for j in 2:(Nx1 - 1)
        xj = coords.xp[j]
        dx_c = xj - xcenter
        for i in 2:(Ny1 - 1)
            yi = coords.yp[i]
            dy_c = yi - ycenter
            r2 = dx_c^2 + dy_c^2
            if r2 <= rplanet2
                k_node = 0.25 * (KX[i, j - 1] + KX[i, j] + KY[i - 1, j] + KY[i, j])
                dT_dx = (tk[i, j + 1] - tk[i, j - 1]) * inv_2dx
                dT_dy = (tk[i + 1, j] - tk[i - 1, j]) * inv_2dy
                r_eff2 = r2 + eps_r2
                Q_metric[i, j] = (k_node / r_eff2) * (dx_c * dT_dx + dy_c * dT_dy)
            end
        end
    end
    return Q_metric
end

"""
Compute silicate partial melt fraction based on temperature, pressure, and rock type.

$(SIGNATURES)

# Details
- `T`: marker temperature [K]
- `P`: marker pressure [Pa] (fluid or lithostatic)
- `tm`: marker phase material type (1: core/mantle, 2: crust, 3: air)
- `T_sol`: solidus temperature [K]
- `T_liq`: liquidus temperature [K]
- `dpdt`: Clapeyron slope dT_sol/dP [K/Pa]

# Returns
- `F_m`: silicate melt mass/volume fraction in [0, 1]
"""
function compute_melt_fraction(
    T::Real, P::Real, tm::Integer; T_sol::Real=1400.0, T_liq::Real=1800.0, dpdt::Real=0.0
)
    if !isfinite(T) || T <= 0.0
        throw(DomainError(T, "Absolute temperature must be positive and finite"))
    end
    if !isfinite(P)
        throw(DomainError(P, "Pressure must be finite"))
    end
    if !isfinite(T_sol) || !isfinite(T_liq) || T_sol >= T_liq
        throw(
            DomainError(
                (T_sol, T_liq),
                "Solidus temperature must be finite and strictly less than liquidus temperature",
            ),
        )
    end
    if tm >= 3
        return 0.0
    end
    T_s = T_sol + dpdt * max(0.0, P)
    T_l = T_liq + dpdt * max(0.0, P)
    if T_s >= T_l
        throw(
            DomainError(
                (T_s, T_l),
                "Pressure-shifted solidus must be strictly less than shifted liquidus",
            ),
        )
    end
    if T <= T_s
        return 0.0
    elseif T >= T_l
        return 1.0
    else
        return (T - T_s) / (T_l - T_s)
    end
end

"""
Compute apparent volumetric heat capacity of silicate rock including latent heat of melting.

$(SIGNATURES)

# Details
- `T`: marker temperature [K]
- `P`: marker pressure [Pa]
- `rhocp_solid`: sensible solid volumetric heat capacity [J/(m^3 K)]
- `rho_solid`: solid density [kg/m^3]
- `tm`: marker phase material type (1: core/mantle, 2: crust, 3: air)
- `T_sol`: solidus temperature [K]
- `T_liq`: liquidus temperature [K]
- `L_melt`: latent heat of silicate melting [J/kg]
- `active`: boolean flag to enable latent heat addition
- `dpdt`: Clapeyron slope dT_sol/dP [K/Pa]

# Returns
- `rhocp_eff`: effective volumetric heat capacity [J/(m^3 K)]
"""
function rhocp_apparent_silicate(
    T::Real,
    P::Real,
    rhocp_solid::Real,
    rho_solid::Real,
    tm::Integer;
    T_sol::Real=1400.0,
    T_liq::Real=1800.0,
    L_melt::Real=4.0e5,
    active::Bool=true,
    dpdt::Real=0.0,
)
    if !isfinite(T) || T <= 0.0
        throw(DomainError(T, "Absolute temperature must be positive and finite"))
    end
    if !isfinite(P)
        throw(DomainError(P, "Pressure must be finite"))
    end
    if !isfinite(T_sol) || !isfinite(T_liq) || T_sol >= T_liq
        throw(
            DomainError(
                (T_sol, T_liq),
                "Solidus temperature must be finite and strictly less than liquidus temperature",
            ),
        )
    end
    if !active || tm >= 3
        return rhocp_solid
    end
    T_s = T_sol + dpdt * max(0.0, P)
    T_l = T_liq + dpdt * max(0.0, P)
    if T_s >= T_l
        throw(
            DomainError(
                (T_s, T_l),
                "Pressure-shifted solidus must be strictly less than shifted liquidus",
            ),
        )
    end
    if T_s < T < T_l
        dFdT = inv(T_l - T_s)
        return rhocp_solid + rho_solid * L_melt * dFdT
    else
        return rhocp_solid
    end
end

"""
Compute melt-weakened matrix viscosity and suspension transition.

$(SIGNATURES)

# Details
- `eta_solid`: solid rock matrix viscosity [Pa s]
- `F_m`: silicate melt fraction in [0, 1]
- `tm`: marker phase material type (1: core/mantle, 2: crust, 3: air)
- `alpha_eta`: melt weakening exponent (Costa et al. 2009; Gerya 2019 Section 16.6.2)
- `phi_crit`: rheologically critical melt fraction for disaggregation
- `eta_melt`: pure liquid melt viscosity limit [Pa s]
- `etamin`: lower viscosity clamp [Pa s]
- `etamax`: upper viscosity clamp [Pa s]

# Returns
- `eta_eff`: effective shear viscosity [Pa s]
"""
function compute_melt_weakened_viscosity(
    eta_solid::Real,
    F_m::Real,
    tm::Integer;
    alpha_eta::Real=28.0,
    phi_crit::Real=0.4,
    eta_melt::Real=10.0,
    etamin::Real=1.0e12,
    etamax::Real=1.0e23,
)
    if !isfinite(F_m)
        throw(DomainError(F_m, "Melt fraction must be finite"))
    end
    if !isfinite(eta_solid) || eta_solid <= 0.0
        throw(DomainError(eta_solid, "Solid viscosity must be positive and finite"))
    end
    if tm >= 3
        return clamp(eta_solid, etamin, etamax)
    end
    F_clamped = clamp(F_m, 0.0, 1.0)
    if iszero(F_clamped)
        return clamp(eta_solid, etamin, etamax)
    elseif F_clamped < phi_crit
        eta_weak = eta_solid * exp(-alpha_eta * F_clamped)
        return clamp(eta_weak, etamin, etamax)
    else
        # Critical disaggregation into magma ocean suspension
        eta_at_crit = eta_solid * exp(-alpha_eta * phi_crit)
        frac = (F_clamped - phi_crit) / (1.0 - phi_crit)
        log_eta = (1.0 - frac) * log(eta_at_crit) + frac * log(eta_melt)
        eta_susp = exp(log_eta)
        return clamp(eta_susp, etamin, etamax)
    end
end

"""
Compute regularized effective thermal conductivity from soft turbulence.

$(SIGNATURES)

Blends conductive thermal conductivity and turbulent convective conductivity
smoothly in logarithmic space across a melt fraction transition window.

# Arguments
- `k_cond`: Conductive thermal conductivity [W/(m K)]
- `eta_num`: Numerical shear viscosity used in momentum solver [Pa s]
- `eta_fluid`: Physical fluid/magma viscosity [Pa s]
- `F_m`: Silicate melt fraction [0, 1]
- `T_marker`: Local marker temperature [K]
- `T_surface`: Reference surface/ambient temperature [K]

# Keyword Arguments
- `turb_exponent`: Power exponent for viscosity ratio (default 1/3 from Solomatov 2007)
- `F_start`: Lower boundary of transition window [0, 1] (default 0.30)
- `F_end`: Upper boundary of transition window [0, 1] (default 0.50)
- `dT_min`: Temperature contrast scale [K] (default 10.0)
- `k_floor`: Minimum thermal conductivity [W/(m K)] (default 1.0e-3)
- `k_cutoff`: Maximum thermal conductivity [W/(m K)] (default 1.0e6)

# Returns
- `k_eff`: Effective thermal conductivity [W/(m K)]
"""
function regularized_soft_turbulence_conductivity(
    k_cond::Real,
    eta_num::Real,
    eta_fluid::Real,
    F_m::Real,
    T_marker::Real,
    T_surface::Real;
    turb_exponent::Real=1.0 / 3.0,
    F_start::Real=0.30,
    F_end::Real=0.50,
    dT_min::Real=10.0,
    k_floor::Real=1.0e-3,
    k_cutoff::Real=1.0e6,
)
    if !isfinite(k_cond) ||
        k_cond <= 0.0 ||
        !isfinite(eta_num) ||
        eta_num <= 0.0 ||
        !isfinite(eta_fluid) ||
        eta_fluid <= 0.0
        throw(
            DomainError(
                (k_cond, eta_num, eta_fluid),
                "k_cond and viscosities must be positive and finite",
            ),
        )
    end
    if !isfinite(F_m) || !isfinite(T_marker) || !isfinite(T_surface)
        throw(
            DomainError(
                (F_m, T_marker, T_surface), "Melt fraction and temperatures must be finite"
            ),
        )
    end
    if F_end <= F_start
        throw(DomainError((F_start, F_end), "F_end must be strictly greater than F_start"))
    end
    if !isfinite(turb_exponent) || turb_exponent <= 0.0
        throw(DomainError(turb_exponent, "Turbulence exponent must be positive and finite"))
    end

    # Melt fraction weight factor using cubic smoothstep
    xi = clamp((F_m - F_start) / (F_end - F_start), 0.0, 1.0)
    w_F = xi * xi * (3.0 - 2.0 * xi)

    # Temperature contrast weight factor
    dT = max(0.0, T_marker - T_surface)
    w_T = clamp(dT / dT_min, 0.0, 1.0)^2

    w_total = w_F * w_T
    if w_total <= 0.0
        return max(k_cond, clamp(k_cond, k_floor, k_cutoff))
    end

    # Target turbulent conductivity from Solomatov (2007) scaling reduction
    k_turb_raw = k_cond * (eta_num / eta_fluid)^turb_exponent
    k_turb = max(k_cond, clamp(k_turb_raw, k_floor, k_cutoff))

    # Geometric blend in logarithmic space
    log_k = (1.0 - w_total) * log10(k_cond) + w_total * log10(k_turb)
    return max(k_cond, clamp(10.0^log_k, k_floor, k_cutoff))
end

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

# =============================================================================
# Atmospheric Jeans Kinetic Escape & Volatile Inventory Dynamics
# =============================================================================

# Fundamental physical constants
"""
Newtonian gravitational constant G in m^3 / (kg s^2), from CODATA 2018.
"""
const GRAVITATIONAL_CONSTANT = 6.67430e-11

"""
Boltzmann constant k_B in J / K (SI exact definition).
"""
const BOLTZMANN_CONSTANT = 1.380649e-23

"""
Avogadro constant N_A in 1/mol (SI exact definition).
"""
const AVOGADRO_CONSTANT = 6.02214076e23

# Volatile molecular masses [kg] (Standard atomic weights divided by Avogadro constant)
"""
Molecular mass of water (H2O) in kilograms [kg], with molar mass 18.01528 g/mol.
"""
const MASS_H2O_KG = 2.991507e-26

"""
Molecular mass of molecular hydrogen (H2) in kilograms [kg], with molar mass 2.01588 g/mol.
"""
const MASS_H2_KG = 3.347447e-27

"""
Molecular mass of molecular nitrogen (N2) in kilograms [kg], with molar mass 28.01340 g/mol.
"""
const MASS_N2_KG = 4.651735e-26

"""
Molecular mass of ammonia (NH3) in kilograms [kg], with molar mass 17.03052 g/mol.
"""
const MASS_NH3_KG = 2.827986e-26

"""
Molecular mass of carbon monoxide (CO) in kilograms [kg], with molar mass 28.01010 g/mol.
"""
const MASS_CO_KG = 4.651187e-26

"""
Molecular mass of carbon dioxide (CO2) in kilograms [kg], with molar mass 44.00950 g/mol.
"""
const MASS_CO2_KG = 7.307950e-26

"""
Molecular mass of methane (CH4) in kilograms [kg], with molar mass 16.04246 g/mol.
"""
const MASS_CH4_KG = 2.663920e-26

"""
Molecular mass of hydrogen sulfide (H2S) in kilograms [kg], with molar mass 34.08088 g/mol.
"""
const MASS_H2S_KG = 5.659267e-26

"""
Molecular mass of disulfur (S2) in kilograms [kg], with molar mass 64.130 g/mol.
"""
const MASS_S2_KG = 1.064904e-25

"""
Molecular mass of sulfur dioxide (SO2) in kilograms [kg], with molar mass 64.066 g/mol.
"""
const MASS_SO2_KG = 1.063841e-25

"""
Upper threshold on the Jeans parameter λ above which kinetic effusion is numerically negligible.
"""
const JEANS_LAMBDA_CUTOFF = 100.0

"""
Lower threshold on the Jeans parameter λ below which escape is pure hydrodynamic sound-speed blow-off.
"""
const HYDRODYNAMIC_ESCAPE_LAMBDA_LOW = 1.0

"""
Upper threshold on the Jeans parameter λ above which escape transitions fully to kinetic effusion.
"""
const HYDRODYNAMIC_ESCAPE_LAMBDA_CUTOFF = 2.0

"""
Retrieve molecular mass in kilograms for standard planetary volatile species.

$(SIGNATURES)

# Arguments
- `species::Symbol`: Volatile species identifier (`:H2O`, `:H2`, `:N2`, `:NH3`, `:CO`, `:CO2`, `:CH4`, `:H2S`, `:S2`, `:SO2`).

# Returns
- `mass::Float64`: Molecular mass [kg].
"""
function get_species_molecular_mass(species::Symbol)::Float64
    s = Symbol(uppercase(String(species)))
    if s === :H2O
        return MASS_H2O_KG
    elseif s === :H2
        return MASS_H2_KG
    elseif s === :N2
        return MASS_N2_KG
    elseif s === :NH3
        return MASS_NH3_KG
    elseif s === :CO
        return MASS_CO_KG
    elseif s === :CO2
        return MASS_CO2_KG
    elseif s === :CH4
        return MASS_CH4_KG
    elseif s === :H2S
        return MASS_H2S_KG
    elseif s === :S2
        return MASS_S2_KG
    elseif s === :SO2
        return MASS_SO2_KG
    else
        throw(
            ArgumentError(
                "Unknown species: $species. Supported species: :H2O, :H2, :N2, :NH3, :CO, :CO2, :CH4, :H2S, :S2, :SO2",
            ),
        )
    end
end

"""
Compute planetary escape velocity at a given radial distance.

$(SIGNATURES)

    v_esc = sqrt(2 * G * M / r)

# Arguments
- `M_planet::Real`: Planetary mass [kg]
- `r::Real`: Radial distance from planetary center [m]

# Returns
- `v_esc`: Escape velocity [m/s]
"""
function compute_escape_velocity(M_planet::Real, r::Real)::Float64
    M_val = Float64(M_planet)
    if M_val <= 0.0 || !isfinite(M_val)
        throw(DomainError(M_val, "Planetary mass must be > 0 and finite"))
    end
    r_val = Float64(r)
    if r_val <= 0.0 || !isfinite(r_val)
        throw(DomainError(r_val, "Radius must be > 0 and finite"))
    end
    return sqrt(2.0 * GRAVITATIONAL_CONSTANT * M_val / r_val)
end

"""
Compute most probable thermal speed of a gas particle.

$(SIGNATURES)

    v_th = sqrt(2 * k_B * T / m)

# Arguments
- `T_K::Real`: Gas temperature [K]
- `m_species_kg::Real`: Particle molecular mass [kg]

# Returns
- `v_th`: Most probable thermal speed [m/s]
"""
function compute_thermal_velocity(T_K::Real, m_species_kg::Real)::Float64
    T_val = Float64(T_K)
    if T_val <= 0.0 || !isfinite(T_val)
        throw(DomainError(T_val, "Temperature must be > 0 and finite"))
    end
    m_val = Float64(m_species_kg)
    if m_val <= 0.0 || !isfinite(m_val)
        throw(DomainError(m_val, "Molecular mass must be > 0 and finite"))
    end
    return sqrt(2.0 * BOLTZMANN_CONSTANT * T_val / m_val)
end

"""
Compute dimensionless Jeans escape parameter λ at the exobase.

$(SIGNATURES)

    λ = (G * M * m) / (k_B * T_exo * R_exo) = (v_esc / v_th)^2

# Arguments
- `M_planet::Real`: Planetary mass [kg]
- `r_exo_m::Real`: Exobase radius [m]
- `T_exo_K::Real`: Exobase temperature [K]
- `m_species_kg::Real`: Particle molecular mass [kg]

# Returns
- `lambda`: Dimensionless Jeans parameter
"""
function compute_jeans_parameter(
    M_planet::Real, r_exo_m::Real, T_exo_K::Real, m_species_kg::Real
)::Float64
    M_val = Float64(M_planet)
    if M_val <= 0.0 || !isfinite(M_val)
        throw(DomainError(M_val, "Planetary mass must be > 0 and finite"))
    end
    r_val = Float64(r_exo_m)
    if r_val <= 0.0 || !isfinite(r_val)
        throw(DomainError(r_val, "Exobase radius must be > 0 and finite"))
    end
    T_val = Float64(T_exo_K)
    if T_val <= 0.0 || !isfinite(T_val)
        throw(DomainError(T_val, "Exobase temperature must be > 0 and finite"))
    end
    m_val = Float64(m_species_kg)
    if m_val <= 0.0 || !isfinite(m_val)
        throw(DomainError(m_val, "Molecular mass must be > 0 and finite"))
    end
    return (GRAVITATIONAL_CONSTANT * M_val * m_val) / (BOLTZMANN_CONSTANT * T_val * r_val)
end

"""
Compute Jeans kinetic escape particle number flux across the exobase.

$(SIGNATURES)

    Φ_Jeans = (n_exo * v_th) / (2 * sqrt(π)) * (1 + λ) * exp(-λ)

# Arguments
- `n_exo::Real`: Particle number density at exobase [m^-3]
- `T_exo_K::Real`: Exobase temperature [K]
- `m_species_kg::Real`: Particle molecular mass [kg]
- `lambda::Real`: Dimensionless Jeans parameter

# Returns
- `Phi_Jeans`: Kinetic escape number flux [m^-2 s^-1]
"""
function compute_jeans_escape_flux(
    n_exo::Real,
    T_exo_K::Real,
    m_species_kg::Real,
    lambda::Real;
    gamma::Real=1.4,
    hydrodynamic::Bool=true,
)::Float64
    n_val = Float64(n_exo)
    if !isfinite(n_val)
        throw(DomainError(n_val, "Number density must be finite"))
    end
    if n_val <= 0.0
        return 0.0
    end
    lam_val = Float64(lambda)
    if !isfinite(lam_val) || lam_val < 0.0
        throw(DomainError(lam_val, "Jeans parameter must be non-negative and finite"))
    end
    # Strong gravitational retention underflow guard
    if lam_val > JEANS_LAMBDA_CUTOFF
        return 0.0
    end
    T_val = Float64(T_exo_K)
    m_val = Float64(m_species_kg)
    if hydrodynamic && lam_val < HYDRODYNAMIC_ESCAPE_LAMBDA_CUTOFF
        c_s = sqrt(Float64(gamma) * BOLTZMANN_CONSTANT * T_val / m_val)
        flux_hydro = n_val * c_s
        v_th = compute_thermal_velocity(T_val, m_val)
        effusion_factor = (1.0 + lam_val) * exp(-lam_val)
        flux_eff = (n_val * v_th / (2.0 * sqrt(π))) * effusion_factor
        s = clamp(
            (lam_val - HYDRODYNAMIC_ESCAPE_LAMBDA_LOW) /
            (HYDRODYNAMIC_ESCAPE_LAMBDA_CUTOFF - HYDRODYNAMIC_ESCAPE_LAMBDA_LOW),
            0.0,
            1.0,
        )
        w = s * s * (3.0 - 2.0 * s)
        return (1.0 - w) * flux_hydro + w * flux_eff
    end
    v_th = compute_thermal_velocity(T_val, m_val)
    effusion_factor = (1.0 + lam_val) * exp(-lam_val)
    return (n_val * v_th / (2.0 * sqrt(π))) * effusion_factor
end

"""
Compute global planetary Jeans mass loss rate for a given volatile species.

$(SIGNATURES)

    dM/dt = 4 * π * R_exo^2 * ρ_exo * (v_th / (2 * sqrt(π))) * (1 + λ) * exp(-λ)

# Arguments
- `M_planet::Real`: Planetary mass [kg]
- `R_exo_m::Real`: Exobase radius [m]
- `T_exo_K::Real`: Exobase temperature [K]
- `m_species_kg::Real`: Particle molecular mass [kg]
- `rho_exo_kg_m3::Real`: Gas mass density at exobase [kg/m^3]

# Keyword Arguments
- `gamma::Real`: Adiabatic index for hydrodynamic blow-off regime (default: 1.4)
- `hydrodynamic::Bool`: Enable hydrodynamic sound speed escape when λ < 2.0 (default: true)

# Returns
- `loss_rate`: Mass escape rate [kg/s]
"""
function compute_jeans_mass_loss_rate(
    M_planet::Real,
    R_exo_m::Real,
    T_exo_K::Real,
    m_species_kg::Real,
    rho_exo_kg_m3::Real;
    gamma::Real=1.4,
    hydrodynamic::Bool=true,
)::Float64
    rho_val = Float64(rho_exo_kg_m3)
    if !isfinite(rho_val)
        throw(DomainError(rho_val, "Exobase density must be finite"))
    end
    if rho_val <= 0.0
        return 0.0
    end
    R_val = Float64(R_exo_m)
    if R_val <= 0.0 || !isfinite(R_val)
        throw(DomainError(R_val, "Exobase radius must be > 0 and finite"))
    end
    lam = compute_jeans_parameter(M_planet, R_val, T_exo_K, m_species_kg)
    if lam > JEANS_LAMBDA_CUTOFF
        return 0.0
    end
    area = 4.0 * π * R_val^2
    T_val = Float64(T_exo_K)
    m_val = Float64(m_species_kg)
    if hydrodynamic && lam < HYDRODYNAMIC_ESCAPE_LAMBDA_CUTOFF
        c_s = sqrt(Float64(gamma) * BOLTZMANN_CONSTANT * T_val / m_val)
        flux_hydro = rho_val * c_s
        v_th = compute_thermal_velocity(T_val, m_val)
        effusion_factor = (1.0 + lam) * exp(-lam)
        flux_eff = rho_val * (v_th / (2.0 * sqrt(π))) * effusion_factor
        s = clamp(
            (lam - HYDRODYNAMIC_ESCAPE_LAMBDA_LOW) /
            (HYDRODYNAMIC_ESCAPE_LAMBDA_CUTOFF - HYDRODYNAMIC_ESCAPE_LAMBDA_LOW),
            0.0,
            1.0,
        )
        w = s * s * (3.0 - 2.0 * s)
        mass_flux = (1.0 - w) * flux_hydro + w * flux_eff
        return area * mass_flux
    end
    v_th = compute_thermal_velocity(T_val, m_val)
    effusion_factor = (1.0 + lam) * exp(-lam)
    mass_flux = rho_val * (v_th / (2.0 * sqrt(π))) * effusion_factor
    return area * mass_flux
end

"""
Compute atmospheric barometric scale height.

$(SIGNATURES)

    H = (k_B * T) / (m * g)

# Arguments
- `M_planet::Real`: Planetary mass [kg]
- `R_planet::Real`: Planetary radius [m]
- `T_exo_K::Real`: Atmospheric temperature [K]
- `m_species_kg::Real`: Particle molecular mass [kg]

# Returns
- `H`: Scale height [m]
"""
function compute_atmospheric_scale_height(
    M_planet::Real, R_planet::Real, T_exo_K::Real, m_species_kg::Real
)::Float64
    M_val = Float64(M_planet)
    if M_val <= 0.0 || !isfinite(M_val)
        throw(DomainError(M_val, "Planetary mass must be > 0 and finite"))
    end
    R_val = Float64(R_planet)
    if R_val <= 0.0 || !isfinite(R_val)
        throw(DomainError(R_val, "Planetary radius must be > 0 and finite"))
    end
    T_val = Float64(T_exo_K)
    if T_val <= 0.0 || !isfinite(T_val)
        throw(DomainError(T_val, "Temperature must be > 0 and finite"))
    end
    m_val = Float64(m_species_kg)
    if m_val <= 0.0 || !isfinite(m_val)
        throw(DomainError(m_val, "Molecular mass must be > 0 and finite"))
    end
    g = GRAVITATIONAL_CONSTANT * M_val / (R_val^2)
    return (BOLTZMANN_CONSTANT * T_val) / (m_val * g)
end

"""
Compute surface atmospheric pressure from total atmospheric mass.

$(SIGNATURES)

    P_surf = (M_atm * g) / (4 * π * R_planet^2)

# Arguments
- `M_atm_total::Real`: Total atmospheric mass [kg]
- `M_planet::Real`: Planetary mass [kg]
- `R_planet::Real`: Planetary radius [m]

# Returns
- `P_surf`: Surface pressure [Pa]
"""
function compute_surface_atmospheric_pressure(
    M_atm_total::Real, M_planet::Real, R_planet::Real
)::Float64
    M_atm_val = Float64(M_atm_total)
    if !isfinite(M_atm_val) || M_atm_val < 0.0
        throw(DomainError(M_atm_val, "Atmospheric mass must be non-negative and finite"))
    end
    M_val = Float64(M_planet)
    if M_val <= 0.0 || !isfinite(M_val)
        throw(DomainError(M_val, "Planetary mass must be > 0 and finite"))
    end
    R_val = Float64(R_planet)
    if R_val <= 0.0 || !isfinite(R_val)
        throw(DomainError(R_val, "Planetary radius must be > 0 and finite"))
    end
    if M_atm_val == 0.0
        return 0.0
    end
    g = GRAVITATIONAL_CONSTANT * M_val / (R_val^2)
    area = 4.0 * π * R_val^2
    return M_atm_val * g / area
end

"""
Evolve atmospheric species mass inventory over a time step dt under venting and Jeans escape.

$(SIGNATURES)

Integrates the first-order ODE:
    dM_atm / dt = M_vent_rate - k_escape * M_atm

where k_escape is determined from Jeans kinetic flux and effective atmospheric column scale.
Mass is conserved to machine precision:
    M_atm(t + dt) + M_escaped_step == M_atm(t) + M_vent_rate * dt

# Arguments
- `M_atm_prev::Real`: Initial atmospheric mass of species [kg]
- `M_vent_rate::Real`: Venting mass influx rate from interior [kg/s]
- `dt_s::Real`: Time step duration [s]
- `M_planet::Real`: Planetary mass [kg]
- `R_planet::Real`: Planetary radius [m]
- `T_exo::Real`: Exobase temperature [K]
- `m_species::Real`: Molecular mass of volatile species [kg]
- `R_exobase::Real=R_planet`: Exobase radius for escape evaluation [m]

# Keyword Arguments
- `gamma::Real`: Adiabatic index for hydrodynamic blow-off regime (default: 1.4)
- `hydrodynamic::Bool`: Enable hydrodynamic blow-off when λ < 2.0 (default: true)

# Returns
- `NamedTuple`: `(; M_atm, M_escaped_step, escape_rate)`
"""
function evolve_atmospheric_species_inventory(
    M_atm_prev::Real,
    M_vent_rate::Real,
    dt_s::Real,
    M_planet::Real,
    R_planet::Real,
    T_exo::Real,
    m_species::Real;
    R_exobase::Real=R_planet,
    gamma::Real=1.4,
    hydrodynamic::Bool=true,
)::@NamedTuple{M_atm::Float64, M_escaped_step::Float64, escape_rate::Float64}
    M_prev = Float64(M_atm_prev)
    if !isfinite(M_prev) || M_prev < 0.0
        throw(
            DomainError(M_prev, "Initial atmospheric mass must be non-negative and finite")
        )
    end
    M_dot_vent = Float64(M_vent_rate)
    if !isfinite(M_dot_vent) || M_dot_vent < 0.0
        throw(DomainError(M_dot_vent, "Venting rate must be non-negative and finite"))
    end
    dt = Float64(dt_s)
    if !isfinite(dt) || dt < 0.0
        throw(DomainError(dt, "Time step dt must be non-negative and finite"))
    end
    if dt == 0.0
        return (M_atm=M_prev, M_escaped_step=0.0, escape_rate=0.0)
    end

    R_val = Float64(R_planet)
    if R_val <= 0.0 || !isfinite(R_val)
        throw(DomainError(R_val, "Planetary radius must be > 0 and finite"))
    end
    R_exo_val = Float64(R_exobase)
    if R_exo_val < R_val || !isfinite(R_exo_val)
        throw(
            DomainError(R_exo_val, "Exobase radius must be >= planetary radius and finite")
        )
    end

    lam = compute_jeans_parameter(M_planet, R_exo_val, T_exo, m_species)
    v_th = compute_thermal_velocity(T_exo, m_species)

    # Compute loss rate coefficient k_escape [s^-1]
    k_escape = if lam > JEANS_LAMBDA_CUTOFF
        0.0
    elseif hydrodynamic && lam < HYDRODYNAMIC_ESCAPE_LAMBDA_CUTOFF
        c_s = sqrt(Float64(gamma) * BOLTZMANN_CONSTANT * Float64(T_exo) / Float64(m_species))
        k_hydro = c_s / R_exo_val
        H = compute_atmospheric_scale_height(M_planet, R_exo_val, T_exo, m_species)
        effusion_factor = (1.0 + lam) * exp(-lam)
        k_eff = (v_th / (2.0 * sqrt(π) * H)) * effusion_factor
        s = clamp(
            (lam - HYDRODYNAMIC_ESCAPE_LAMBDA_LOW) /
            (HYDRODYNAMIC_ESCAPE_LAMBDA_CUTOFF - HYDRODYNAMIC_ESCAPE_LAMBDA_LOW),
            0.0,
            1.0,
        )
        w = s * s * (3.0 - 2.0 * s)
        (1.0 - w) * k_hydro + w * k_eff
    else
        H = compute_atmospheric_scale_height(M_planet, R_exo_val, T_exo, m_species)
        effusion_factor = (1.0 + lam) * exp(-lam)
        (v_th / (2.0 * sqrt(π) * H)) * effusion_factor
    end

    x = k_escape * dt
    # Analytical solution of dM/dt = S - k*M
    M_next = if x < 1.0e-6
        # Taylor expansion to prevent cancellation or divide-by-zero
        int_factor = dt * (1.0 - 0.5 * x + (x^2) / 6.0)
        M_prev * exp(-x) + M_dot_vent * int_factor
    elseif x > 40.0
        # Fully decoupled rapid escape regime
        M_dot_vent / k_escape
    else
        M_prev * exp(-x) + (M_dot_vent / k_escape) * (1.0 - exp(-x))
    end

    # Guarantee non-negative mass
    M_next = max(0.0, M_next)

    # Mass balance: escaped mass is exactly influx plus initial minus final
    total_available = M_prev + M_dot_vent * dt
    M_escaped = max(0.0, total_available - M_next)
    current_escape_rate = k_escape * M_next

    return (M_atm=M_next, M_escaped_step=M_escaped, escape_rate=current_escape_rate)
end

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
    if !isfinite(F_m)
        throw(DomainError(F_m, "Silicate melt fraction must be finite"))
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

    # Percolation component
    v_perc =
        if percolation_active &&
            F_m < F_perc_end &&
            phi_m > phi_crit_perc &&
            phi_m > phi_residual
            phi_m_perm = min(phi_m, 1.0 - 1e-7)
            km = metal_permeability(
                phi_m_perm;
                k_metal_ref=k_metal_ref,
                phi_crit_perc=phi_crit_perc,
                perm_exponent=perm_exponent,
                phi0=phi0,
            )
            phi_mobile = phi_m - phi_residual
            (km / (phi_m * eta_metal)) * drho * g_acc * (phi_mobile / phi_m)
        else
            0.0
        end

    # Settling component
    v_settle = if settling_active && F_m >= F_settle_start && phi_m > 0.0
        v0 = stokes_settling_velocity(
            r_drop,
            drho,
            g_acc,
            eta_susp;
            hadamard_rybczynski=hadamard_rybczynski,
            eta_metal=eta_metal,
        )
        h = richardson_zaki_hindrance(
            phi_m; hindered_exponent=hindered_exponent, phi_pack=phi_pack
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
        if F_m <= F_settle_start
            return v_perc
        elseif F_m >= F_perc_end
            return v_settle
        else
            xi = (F_m - F_settle_start) / (F_perc_end - F_settle_start)
            w = xi * xi * (3.0 - 2.0 * xi)
            return (1.0 - w) * v_perc + w * v_settle
        end
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
    M_i = m_sil * C_i_sil + m_met * C_i_met = const

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

    # Interacting phase masses per unit marker volume
    m_met = Xfem[m] * max(Float64(rho_metal), 100.0)
    m_sil = phi_sil * max(Float64(rho_silicate), 100.0)
    if m_met <= 0.0 || m_sil <= 0.0
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
        C_sil = XCm[m]
        C_met = Xfe_C_m[m]
        M_tot = m_sil * C_sil + m_met * C_met
        denom = m_sil + D_C * m_met
        if denom > 0.0
            C_sil_eq = M_tot / denom
            C_met_eq = D_C * C_sil_eq
            C_met_C_max = min(Float64(cfg.D_max), 7.0e4)
            if C_met_eq > C_met_C_max
                C_met_eq = C_met_C_max
                C_sil_eq = max(0.0, (M_tot - m_met * C_met_eq) / m_sil)
            end
            dC_sil = alpha_eq * (C_sil_eq - C_sil)
            dC_sil = max(dC_sil, -C_sil)
            dC_met = -dC_sil * (m_sil / m_met)
            if C_met + dC_met < 0.0
                dC_met = -C_met
                dC_sil = -dC_met * (m_met / m_sil)
            end
            XCm[m] = max(0.0, C_sil + dC_sil)
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
        C_sil = XNm[m]
        C_met = Xfe_N_m[m]
        M_tot = m_sil * C_sil + m_met * C_met
        denom = m_sil + D_N * m_met
        if denom > 0.0
            C_sil_eq = M_tot / denom
            C_met_eq = D_N * C_sil_eq
            C_met_N_max = min(Float64(cfg.D_max), 4.0e4)
            if C_met_eq > C_met_N_max
                C_met_eq = C_met_N_max
                C_sil_eq = max(0.0, (M_tot - m_met * C_met_eq) / m_sil)
            end
            dC_sil = alpha_eq * (C_sil_eq - C_sil)
            dC_sil = max(dC_sil, -C_sil)
            dC_met = -dC_sil * (m_sil / m_met)
            if C_met + dC_met < 0.0
                dC_met = -C_met
                dC_sil = -dC_met * (m_met / m_sil)
            end
            XNm[m] = max(0.0, C_sil + dC_sil)
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
        C_sil = XSm[m]
        C_met = Xfe_S_m[m]
        M_tot = m_sil * C_sil + m_met * C_met
        denom = m_sil + D_S * m_met
        if denom > 0.0
            C_sil_eq = M_tot / denom
            C_met_eq = D_S * C_sil_eq
            C_met_S_max = min(Float64(cfg.D_max), 3.65e5)
            if C_met_eq > C_met_S_max
                C_met_eq = C_met_S_max
                C_sil_eq = max(0.0, (M_tot - m_met * C_met_eq) / m_sil)
            end
            dC_sil = alpha_eq * (C_sil_eq - C_sil)
            dC_sil = max(dC_sil, -C_sil)
            dC_met = -dC_sil * (m_sil / m_met)
            if C_met + dC_met < 0.0
                dC_met = -C_met
                dC_sil = -dC_met * (m_met / m_sil)
            end
            XSm[m] = max(0.0, C_sil + dC_sil)
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
    C_sil_H = XH2Om[m] * f_H
    C_met_H = Xfe_H_m[m]
    M_tot_H = m_sil * C_sil_H + m_met * C_met_H
    denom_H = m_sil + D_H * m_met
    if denom_H > 0.0
        C_sil_H_eq = M_tot_H / denom_H
        C_met_H_eq = D_H * C_sil_H_eq
        C_met_H_max = min(Float64(cfg.D_max), 1.0e4)
        if C_met_H_eq > C_met_H_max
            C_met_H_eq = C_met_H_max
            C_sil_H_eq = max(0.0, (M_tot_H - m_met * C_met_H_eq) / m_sil)
        end
        dC_sil = alpha_eq * (C_sil_H_eq - C_sil_H)
        dC_sil = max(dC_sil, -C_sil_H)
        dC_met = -dC_sil * (m_sil / m_met)
        if C_met_H + dC_met < 0.0
            dC_met = -C_met_H
            dC_sil = -dC_met * (m_met / m_sil)
        end
        new_C_sil_H = max(0.0, C_sil_H + dC_sil)
        XH2Om[m] = clamp(new_C_sil_H / f_H, 0.0, 100.0)
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

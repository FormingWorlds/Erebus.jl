
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

# Returns

    - hrsolidm: radiogenic heat production of 26Al [W/m^3]
    - hrfluidm: radiogenic heat production of 60Fe [W/m^3]
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
)
    #26Al: planet ✓, crust ✓, space ×
    if al
        # 26Al radiogenic heat production [W/kg]
        Q_al = Q_radiogenic(f_al, ratio_al, E_al, tau_al, timesum)
        # Solid phase 26Al radiogenic heat production [W/m^3]
        @inbounds hrsolidm = @SVector [Q_al*rhosolidm[1], Q_al*rhosolidm[2], 0.0]
    else
        hrsolidm = @SVector zeros(3)
    end
    #60Fe: planet ✓, crust ×, space ×
    if fe
        # 60Fe radiogenic heat production [W/kg]
        Q_fe = Q_radiogenic(f_fe, ratio_fe, E_fe, tau_fe, timesum)
        # Fluid phase 60Fe radiogenic heat production [W/m^3]
        @inbounds hrfluidm = @SVector [Q_fe*rhofluidm[1], 0.0, 0.0]
    else
        hrfluidm = @SVector zeros(3)
    end
    return hrsolidm, hrfluidm
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
            - 2: constant parameter rhocpfluidm

# Returns
    
        - ρᶠCₚᶠ: volumetric isobaric heat capacity of fluid
"""
function compute_rhocpfluidm(T, mode)
    # @timeit to "compute_rhocpfluidm" begin
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
    # end # @timeit to "compute_rhocpfluidm"
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
    # @timeit to "compute_ksolidm" begin
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
    # end # @timeit to "compute_ksolidm"
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
            - 9: constant parameter ksolidm

# Returns
    
        - kᶠ: thermal conductivity of fluid
"""
function compute_kfluidm(T, mode)
    # @timeit to "compute_kfluidm" begin
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
    # end # @timeit to "compute_kfluidm"
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
    # @timeit to "compute_Δtreaction" begin
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
    # end # @timeit to "compute_Δtreaction"
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
    # @timeit to "compute_gibbs_free_energy" begin
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
    # end # @timeit to "compute_gibbs_free_energy"
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
    # @timeit to "compute_thermodynamic_xfer!" begin
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
    # end # @timeit to "compute_thermodynamic_xfer!"
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
    # @timeit to "perform_thermochemical_reaction!" begin
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
    # end # @timeit to "perform_thermochemical_reaction!"
end # function perform_thermochemical_reaction!

"""
Compute shear heating based on basic (temperature) and P grids.

$(SIGNATURES)

# Details

    - HS: shear heating
    - ETA: viscoplastic viscosity at basic nodes
    - SXY: σ₀xy XY stress at basic nodes
    - ETAP: viscosity at P nodes
    - SXX: σ₀xy XY stress at basic nodes
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
    # @timeit to "compute_shear_heating!" begin
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
    # end # @timeit to "compute_shear_heating!" 
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
    # @timeit to "compute_adiabatic_heating!" begin
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
            if VXP < 0.0
                dpsdx = (ps[i, j]-ps[i, j - 1]) * inv(dx_val)
            else
                dpsdx = (ps[i, j + 1]-ps[i, j]) * inv(dx_val)
            end
            if VYP < 0.0
                dpsdy = (ps[i, j]-ps[i - 1, j]) * inv(dy_val)
            else
                dpsdy = (ps[i + 1, j]-ps[i, j]) * inv(dy_val)
            end
            dpsdt = VXP*dpsdx + VYP*dpsdy
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
    # end # @timeit to "compute_adiabatic_heating!"
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
- For temperatures at or above the triple point (`T >= T0`), the vapor pressure saturates
  at the triple-point value `P0 = 611.66 Pa` because the bulk ice phase transitions to liquid water.
"""
function compute_ice_vapor_pressure(
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

    if T_val >= T0_val
        return P0_val
    end
    return P0_val * exp(-(L_sub_val / Rv_val) * (1.0 / T_val - 1.0 / T0_val))
end

"""
    compute_venting_pressure(T_surf::Real, P_amb::Real; P0::Real=611.66, T0::Real=273.16, L_sub::Real=2.83e6, Rv::Real=461.5)::Float64

Compute effective boundary venting fluid pressure P_vent [Pa] at a planetesimal surface:

    P_vent = max(P_amb, P_sat,ice(T_surf))

Enforces the physical cold-trap constraint: if ambient nebular gas pressure exceeds
ice sublimation pressure at cold surface temperatures, the ambient gas confines pore fluid;
if ambient pressure drops below sublimation pressure (space vacuum), flash sublimation
sets the effective boundary vapor pressure.

# Arguments
- `T_surf`: Planetesimal surface temperature [K]
- `P_amb`: Ambient surrounding gas pressure [Pa]

# Returns
- `P_vent`: Effective venting boundary pressure [Pa]
"""
function compute_venting_pressure(
    T_surf::Real,
    P_amb::Real;
    P0::Real=611.66,
    T0::Real=273.16,
    L_sub::Real=2.83e6,
    Rv::Real=461.5,
)::Float64
    P_sat = compute_ice_vapor_pressure(T_surf; P0=P0, T0=T0, L_sub=L_sub, Rv=Rv)
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

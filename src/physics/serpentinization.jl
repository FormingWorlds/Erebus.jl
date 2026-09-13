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


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
    return solid*(1.0-ϕ) + fluid*ϕ
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
    return (
        sqrt(
            ksolid * kfluid/2 + ((ksolid*(3.0*phi-2.0) + kfluid*(1.0-3.0*phi))^2)*inv(16.0)
        ) - 0.25 * (ksolid*(3.0*phi-2.0) + kfluid*(1.0-3.0*phi))
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
    # phim0 is a global constant defined independent of material type
    return kphim0m * (phimm*inv(phim0))^3.0 * ((1.0-phimm)*inv(1.0-phim0))^-2.0
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
    return ηᶠcur * inv(kϕᵣ) * (phim0*inv(ϕ))^3.0 * ((1.0-ϕ)*inv(1.0-phim0))^2.0
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
    return f * ratio * E * exp(-time*inv(tau)) * inv(tau)
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
    @inbounds etasolidcur = ifelse(tkmm>tmsolidphase, etasolidmm[tmm], etasolidm[tmm])
    @inbounds etafluidcur = ifelse(tkmm>tmfluidphase, etafluidmm[tmm], etafluidm[tmm])
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
        if T < tmfluidphase-5.0
            ρᶠCₚᶠ = ρH₂Oᶠⁱ * 7.67T
        elseif T < tmfluidphase
            ρᶠCₚᶠ = ρH₂Oᶠⁱ * (7.67T + 0.1Lᶠ)
        elseif T < tmfluidphase+5.0
            ρᶠCₚᶠ = ρH₂Oᶠ * (4200.0 + 0.1Lᶠ)
        elseif T < 410.0
            ρᶠCₚᶠ = ρH₂Oᶠ * 4200.0
        else
            ρᶠCₚᶠ = ρH₂Oᶠ * (-4.67e4 + 333T - 0.731T^2 + 5.4e-4T^3)
        end
    elseif mode == 9
        @inbounds ρᶠCₚᶠ = rhocpfluidm[1]
    else
        throw("unknown mode $mode")
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
        kˢ = 0.73 + 1293.0/(T+77.0)
    elseif mode == 9
        @inbounds kˢ = ksolidm[1]
    else
        throw("unknown mode $mode")
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
        if T < tmfluidphase
            kᶠ = 0.465 + 488.0/T
        elseif T < 410
            kᶠ = -0.581 + 6.34e-3T - 7.93e-6T^2
        else
            kᶠ = -0.142 + 4.12e-3T - 5.01e-6T^2
        end
    elseif mode == 9
        @inbounds kᶠ = kfluidm[1]
    else
        throw("unknown mode $mode")
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
function compute_Δtreaction(T, ϕ, mode)
    # @timeit to "compute_Δtreaction" begin
    if mode == 1
        Δtr = -log_completion_rate / (A_I*ϕ) * exp(b_I*(T-c_I)^2)
    elseif mode == 2
        Δtr = -log_completion_rate / (Sxo_B*ϕ) * 2.0^((To_B-T)/Tscl_B)
    elseif mode == 3
        Δtr = -log_completion_rate / (Sxo_B*ϕ) * exp(Ea_T / RG * (1.0/T - 1.0/To_T))
    elseif mode == 9
        Δtr = Δtreaction
    else
        throw("unknown mode $mode")
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
function compute_gibbs_free_energy(T, pf, XDˢ, XWˢ, Δt, Δtr)
    # @timeit to "compute_gibbs_free_energy" begin
    # compute incomplete reaction for short timestep Δt < Δtreaction
    if Δt < Δtr
        # compute ΔG for dehydration reaction (16.145), (16.165b)
        ΔGWD = (ΔHWD - T*ΔSWD + pf*ΔVWD + RG*T*log(XDˢ/XWˢ)) * (1.0 - Δt/Δtr)
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
function compute_relative_enthalpy(Xsolid, XWsolid)
    return -Xsolid * XWsolid * ΔHWD / (MD+MH₂O)
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
function compute_reaction_constant(T, pf, ΔGWD)
    # compute reaction constant (16.151)
    return exp(-(ΔHWD - T*ΔSWD + ΔVWD*pf - ΔGWD) / (RG*T))
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
function compute_thermodynamic_xfer!(DMPSUM, DHPSUM, WTPSUM, DMP, DHP)
    # @timeit to "compute_thermodynamic_xfer!" begin
    Ny1, Nx1 = size(DMP)
    @inbounds begin
        for j in 1:1:Nx1, i in 1:1:Ny1
            if WTPSUM[i, j] > 0.0
                DMP[i, j] = DMPSUM[i, j] * inv(WTPSUM[i, j])
                DHP[i, j] = DHPSUM[i, j] * inv(WTPSUM[i, j])
            else
                DMP[i, j] = DHP[i, j] = zero(0.0)
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
)
    # @timeit to "perform_thermochemical_reaction!" begin
    # reset interpolation arrays
    reset_thermochemical_properties!(DMPSUM, DHPSUM, WTPSUM)
    xp_val = coords === nothing ? xp : coords.xp
    yp_val = coords === nothing ? yp : coords.yp
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy
    jmin_p_val = coords === nothing ? jmin_p : coords.jmin_p
    jmax_p_val = coords === nothing ? jmax_p : coords.jmax_p
    imin_p_val = coords === nothing ? imin_p : coords.imin_p
    imax_p_val = coords === nothing ? imax_p : coords.imax_p
    # iterate over markers
    @inbounds begin
        for m in 1:1:marknum
            # for 
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
                # interpolate non-negative fluid pressure from P nodes
                pfnm = max(zero(0.0), dot4(grid_vector(i, j, pf), weights))
                # factor in previous iteration marker fluid pressure
                if titer > 2
                    pfnm = pfnm*(1.0-pfcoeff) + pfm₀[m]*pfcoeff
                end
                # store current marker fluid pressure for next iteration
                pfm₀[m] = pfnm
                # compute bulk composition of solid and fluid system:
                # compute previous dry solid molar fraction (16.146)
                XDˢm₀ = 1.0 - XWˢm₀[m]
                # get fluid molar volume
                VH₂O = ifelse(tknm>tmfluidphase, VH₂Oᶠ, VH₂Oᶠⁱ)
                # compute previous fluid molar fraction (16.164)
                Xᶠ₀ =
                    phim[m]*(XWˢm₀[m]*VWˢ + XDˢm₀*VDˢ) /
                    ((1.0-phim[m])*VH₂O + phim[m] * (XWˢm₀[m]*VWˢ + XDˢm₀*VDˢ))
                # compute previous equilibrium solid molar fraction (16.150)
                Xˢ₀ = 1.0 - Xᶠ₀
                # compute previous water molar fraction (16.147)
                XH₂Oᵗ = (XWˢm₀[m]*Xˢ₀ + Xᶠ₀) / (1.0 + XWˢm₀[m]*Xˢ₀)
                # compute dry solid molar fraction (16.149)
                XDᵗ = 1.0 - XH₂Oᵗ
                # compute previous solid density (16.161)
                ρˢ₀ = (MD + MH₂O*XWˢm₀[m]) / (VDˢ*XDˢm₀ + VWˢ*XWˢm₀[m])
                # compute previous fluid density (16.162)
                ρᶠ₀ = ifelse(tknm>tmfluidphase, ρH₂Oᶠ, ρH₂Oᶠⁱ)
                # compute Δtreaction
                Δtr = compute_Δtreaction(tknm, phim[m], reaction_rate_coeff_mode)
                # compute previous relative enthalpy of the system (16.163)
                Hᵗ₀ = compute_relative_enthalpy(Xˢ₀, XWˢm₀[m])
                # compute previous ΔG for dehydration reaction (16.165a/b)
                ΔGWD₀ = compute_gibbs_free_energy(tknm, pfnm, XDˢm₀, XWˢm₀[m], Δt, Δtr)
                # compute dehydration reaction constant (16.151)
                KWD = compute_reaction_constant(tknm, pfnm, ΔGWD₀)
                # compute reacted wet solid molar fraction (16.152)
                XWˢm₁ = inv(KWD + 1.0)
                # compute reacted dry solid molar fraction (16.153)
                XDˢm₁ = 1.0 - XWˢm₁
                # compute reacted total solid molar fraction (16.154)
                Xˢ₁ = XDᵗ / (1.0 - XDᵗ*XWˢm₁)
                # compute reacted fluid molar fraction (16.155)
                Xᶠ₁ = 1.0 - Xˢ₁
                # only process fluid-bearing rocks
                if 0.0 < Xᶠ₁ < 1.0
                    # compute reacted equilibrium porosity (16.156)
                    ϕ₁ = Xᶠ₁*VH₂O / (Xᶠ₁*VH₂O + Xˢ₁*(XWˢm₁*VWˢ+XDˢm₁*VDˢ))
                    # compute equilibrium solid density (16.161)
                    ρˢ₁ = (MD + MH₂O*XWˢm₁) / (VDˢ*XDˢm₁ + VWˢ*XWˢm₁)
                    # compute equilibrium fluid density (16.162)
                    ρᶠ₁ = ifelse(tknm>tmfluidphase, ρH₂Oᶠ, ρH₂Oᶠⁱ)
                    # compute equilibrium relative enthalpy of the system
                    # (16.163)
                    Hᵗ₁ = compute_relative_enthalpy(Xˢ₁, XWˢm₁)
                    # compute enthalpy change
                    ΔHᵗ = Hᵗ₁ - Hᵗ₀
                    # compute previous-to-reacted-equilibrium volume ratio
                    # (16.106)
                    RV = (ρˢ₁*(1.0-ϕ₁) + ρᶠ₁*ϕ₁) / (ρˢ₀*(1.0-phim[m]) + ρᶠ₀*phim[m])
                    # compute mass transfer rate (16.103)
                    Γmass = (ρˢ₀*RV*(1.0-phim[m]) - ρˢ₁*(1.0-ϕ₁)) / Δt
                    # compute mass transfer term (16.112e)
                    ΔMm = (1.0-RV) / Δt
                    # compute enthalpy transfer/latent heating term (16.113)
                    ΔHm = Γmass * ΔHᵗ
                    # update wet solid (melt) molar fraction
                    XWˢm[m]=XWˢm₁
                    # update porosity
                    phinewm[m] = ϕ₁
                    # backload properties during first timestep
                    if timestep==1
                        XWˢm₀[m] = XWˢm[m]
                        phim[m] = phinewm[m]
                    end
                    # interpolate mass, enthalpy transfer terms to P nodes
                    interpolate_add_to_grid!(i, j, weights, ΔMm, DMPSUM)
                    interpolate_add_to_grid!(i, j, weights, ΔHm, DHPSUM)
                    interpolate_add_to_grid!(i, j, weights, one(1.0), WTPSUM)
                end
            end # if tm[m] < 3
        end # for m=1:1:marknum
        # compute thermodynamic properties at P nodes
        compute_thermodynamic_xfer!(DMPSUM, DHPSUM, WTPSUM, DMP, DHP)
    end # @inbounds
    @info "min/max mass transfer term" extrema(DMP)
    @info "min/max enthalpy transfer term" extrema(DHP)
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
    - thermal_buoyancy: toggle thermal expansion (default: true)

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

"""
Assemble the LHS sparse coefficient matrix and fill RHS coefficient vector
of the energy conservation (heat) equation.

$(SIGNATURES)

# Details

	- tk1: current temperature at P nodes
	- RHOCP: volumetric heat capacity at P nodes  
	- KX: thermal conductivity at Vx nodes
	- KY: thermal conductivity at Vy nodes 
	- HR: radioactive heating at P nodes
	- HA: adiabatic heating at P nodes 
	- HS: shear heating at P nodes
    - DHP: latent heating (HL) at P nodes
    - RT: thermal RHS coefficient vector
    - dt: current time step length
    - coords: grid coordinates
    - LT: optional ExtendableSparseMatrix buffer to reuse

# Returns

    - LT: LHS sparse coefficient matrix
"""
function assemble_thermal_lse!(
    tk1,
    RHOCP,
    KX,
    KY,
    HR,
    HA,
    HS,
    DHP,
    RT,
    dt;
    coords=nothing,
    LT=nothing,
    Q_metric=nothing,
    Q_lat=nothing,
    Q_seg=nothing,
    workspace=nothing,
)
    Ny1, Nx1 = size(tk1)
    @unpack_coords coords dx dy

    # fresh or reusable LHS coefficient matrix
    LT_target =
        if workspace !== nothing &&
            hasproperty(workspace, :Ny1) &&
            hasproperty(workspace, :Nx1) &&
            hasproperty(workspace, :LT) &&
            workspace.Ny1 == Ny1 &&
            workspace.Nx1 == Nx1
            workspace.LT
        else
            LT
        end

    LT = if LT_target === nothing
        ExtendableSparseMatrix(Ny1 * Nx1, Ny1 * Nx1)
    else
        if !isempty(LT_target.cscmatrix.nzval)
            nonzeros(LT_target.cscmatrix) .= zero(0.0)
        end
        LT_target
    end
    # reset RHS coefficient vector
    RT .= zero(0.0)
    # compose global thermal matrix LT and coefficient vector RT
    @inbounds begin
        for j in 1:1:Nx1, i in 1:1:Ny1
            # define global index in algebraic space
            gk = (j-1)*Ny1 + i
            # External points
            if i==1 || i==Ny1 || j==1 || j==Nx1
                # thermal equation external points: boundary conditions
                # all locations: ghost unknowns T₃=0 -> 1.0⋅T[i,j]=0.0
                updateindex!(LT, +, 1.0, gk, gk)
                # R[gk] = 0.0 # already done with initialization
                # left boundary: ∂T/∂x=0
                if j == 1
                    updateindex!(LT, +, -1.0, gk, gk+Ny1)
                end
                # right boundary: ∂T/∂x=0
                if j == Nx1
                    updateindex!(LT, +, -1.0, gk, gk-Ny1)
                end
                # top inner boundary: ∂T/∂y=0
                if i==1 && 1<j<Nx1
                    updateindex!(LT, +, -1.0, gk, gk+1)
                end
                # bottom inner boundary: ∂T/∂y=0
                if i==Ny1 && 1<j<Nx1
                    updateindex!(LT, +, -1.0, gk, gk-1)
                end
            else
                # internal points: 2D thermal equation (conservative formulation)
                # extract thermal conductivities
                Kx₁ = KX[i, j - 1]
                Kx₂ = KX[i, j]
                Ky₁ = KY[i - 1, j]
                Ky₂ = KY[i, j]
                # fill system of equations: LHS
                updateindex!(LT, +, -Kx₁*inv(dx_val^2), gk, gk-Ny1) # T₁
                updateindex!(LT, +, -Ky₁*inv(dy_val^2), gk, gk-1) # T₂
                updateindex!(
                    LT,
                    +,
                    (RHOCP[i, j]/dt + (Kx₁+Kx₂)*inv(dx_val^2) + (Ky₁+Ky₂)*inv(dy_val^2)),
                    gk,
                    gk,
                ) # T₃
                updateindex!(LT, +, -Ky₂*inv(dy_val^2), gk, gk+1) # T₄
                updateindex!(LT, +, -Kx₂*inv(dx_val^2), gk, gk+Ny1) # T₅
                # fill system of equations: RHS
                RT[gk] = (
                    RHOCP[i, j]/dt*tk1[i, j] + HR[i, j] + HA[i, j] + HS[i, j] + DHP[i, j]
                )
                if Q_metric !== nothing
                    RT[gk] += Q_metric[i, j]
                end
                if Q_lat !== nothing
                    RT[gk] += Q_lat[i, j]
                end
                if Q_seg !== nothing
                    RT[gk] += Q_seg[i, j]
                end
            end
        end
    end # @inbounds

    flush!(LT) # finalize CSC matrix
    if workspace !== nothing &&
        hasproperty(workspace, :Ny1) &&
        hasproperty(workspace, :Nx1) &&
        hasproperty(workspace, :is_initialized) &&
        workspace.Ny1 == Ny1 &&
        workspace.Nx1 == Nx1
        workspace.is_initialized = true
    end
    return LT
end # function assemble_thermal_lse!

"""
    apply_radiative_surface_boundary!(
        KX, KY, tk, coords, rplanet, xcenter, ycenter, T_amb;
        emissivity=0.9, sigma_sb=5.670374419e-8,
        k_rock=nothing, marker_property_mode=1
    )

Apply linearized Stefan-Boltzmann radiative cooling at planetesimal-sticky air
interface faces.

Modifies interface conductivities in KX and KY using harmonic series resistance:

    1 / U_eff = (Δ / (2 * k_bulk)) + (1 / h_rad)
    k_face = U_eff * Δ = 2 * k_bulk * (h_rad * Δ) / (2 * k_bulk + h_rad * Δ)

where `h_rad = compute_radiation_htc(T_surf, T_amb; emissivity, sigma_sb)`, and `k_bulk`
is local rock thermal conductivity: evaluated dynamically via
`compute_ksolidm(T_surf, marker_property_mode)` when `k_rock === nothing`, or given by
constant scalar `k_rock` when specified.

# Arguments
- `KX`: Horizontal conductivity array at Vx nodes [W/(m K)]
- `KY`: Vertical conductivity array at Vy nodes [W/(m K)]
- `tk`: Temperature field [K] at P nodes
- `coords`: Grid coordinate descriptors
- `rplanet`: Planetesimal radius [m]
- `xcenter`: Center horizontal coordinate [m]
- `ycenter`: Center vertical coordinate [m]
- `T_amb`: Ambient disk temperature [K]
- `emissivity`: Radiative emissivity
- `sigma_sb`: Stefan-Boltzmann constant
- `k_rock`: Rock thermal conductivity [W/(m K)] override
  (default: `nothing` for temperature-dependent calculation)
- `marker_property_mode`: Marker property regime index (default: 1)
- `tau_LW`: Infrared optical depth of atmosphere for greenhouse blanketing (default: 0.0)
"""
function apply_radiative_surface_boundary!(
    KX::AbstractMatrix{Float64},
    KY::AbstractMatrix{Float64},
    tk::AbstractMatrix{Float64},
    coords::GridCoordinates,
    rplanet::Real,
    xcenter::Real,
    ycenter::Real,
    T_amb::Real;
    emissivity::Real=0.9,
    sigma_sb::Real=5.670374419e-8,
    k_rock::Union{Real,Nothing}=nothing,
    marker_property_mode::Int=1,
    phi::Real=0.0,
    kfluid::Real=50.0,
    tau_LW::Real=0.0,
)
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    dx = coords.dx
    dy = coords.dy
    rplanet2 = rplanet^2

    # Horizontal faces KX[i, j] between P(i, j) and P(i, j+1)
    @inbounds for j in 1:(Nx1 - 1)
        xj1 = coords.xp[j] - xcenter
        xj2 = coords.xp[j + 1] - xcenter
        for i in 1:Ny1
            yi = coords.yp[i] - ycenter
            r1_sq = xj1^2 + yi^2
            r2_sq = xj2^2 + yi^2
            is_rock1 = r1_sq <= rplanet2
            is_rock2 = r2_sq <= rplanet2
            if is_rock1 != is_rock2
                T_surf = is_rock1 ? tk[i, j] : tk[i, j + 1]
                h_rad =
                    if T_surf > 0.0 && T_amb > 0.0 && isfinite(T_surf) && isfinite(T_amb)
                        compute_effective_radiation_htc(
                            T_surf, T_amb, tau_LW; emissivity=emissivity, sigma_sb=sigma_sb
                        )
                    else
                        0.0
                    end
                k_rad = h_rad * dx
                k_bulk = if k_rock !== nothing
                    Float64(k_rock)
                elseif phi > 0.0
                    ktotal(
                        compute_ksolidm(T_surf, marker_property_mode),
                        Float64(kfluid),
                        Float64(phi),
                    )
                else
                    compute_ksolidm(T_surf, marker_property_mode)
                end
                KX[i, j] = (2.0 * k_bulk * k_rad) / (2.0 * k_bulk + k_rad)
            end
        end
    end

    # Vertical faces KY[i, j] between P(i, j) and P(i+1, j)
    @inbounds for j in 1:Nx1
        xj = coords.xp[j] - xcenter
        for i in 1:(Ny1 - 1)
            yi1 = coords.yp[i] - ycenter
            yi2 = coords.yp[i + 1] - ycenter
            r1_sq = xj^2 + yi1^2
            r2_sq = xj^2 + yi2^2
            is_rock1 = r1_sq <= rplanet2
            is_rock2 = r2_sq <= rplanet2
            if is_rock1 != is_rock2
                T_surf = is_rock1 ? tk[i, j] : tk[i + 1, j]
                h_rad =
                    if T_surf > 0.0 && T_amb > 0.0 && isfinite(T_surf) && isfinite(T_amb)
                        compute_effective_radiation_htc(
                            T_surf, T_amb, tau_LW; emissivity=emissivity, sigma_sb=sigma_sb
                        )
                    else
                        0.0
                    end
                k_rad = h_rad * dy
                k_bulk = if k_rock !== nothing
                    Float64(k_rock)
                elseif phi > 0.0
                    ktotal(
                        compute_ksolidm(T_surf, marker_property_mode),
                        Float64(kfluid),
                        Float64(phi),
                    )
                else
                    compute_ksolidm(T_surf, marker_property_mode)
                end
                KY[i, j] = (2.0 * k_bulk * k_rad) / (2.0 * k_bulk + k_rad)
            end
        end
    end
    return nothing
end

"""
    compute_face_venting_permeability(
        k_v, breached, ice_sealing, T_surf, peff, sigma_t;
        t_freeze=273.15, dt_seal=10.0, k_seal_min_ratio=1.0e-6,
        kappa_frac=1.0e3, gamma_frac=1.0, k_frac_max=1.0e-9,
    )

Compute effective rock face permeability at the planetesimal surface, accounting for
tensile hydrofracture breaching or cryogenic pore ice sealing.

# Arguments
- `k_v::Real`: Reference matrix permeability [m^2].
- `breached::Bool`: Whether tensile failure has ruptured the rock lid.
- `ice_sealing::Bool`: Whether cryogenic pore ice sealing is active below freezing.
- `T_surf::Real`: Local surface rock temperature [K].
- `peff::Real`: Terzaghi effective stress `P_t - P_f` [Pa].
- `sigma_t::Real`: Rock tensile strength [Pa].

# Keywords
- `t_freeze::Real`: Water freezing temperature [K] (default: 273.15).
- `dt_seal::Real`: Temperature sealing interval [K] (default: 10.0).
- `k_seal_min_ratio::Real`: Minimum residual cryogenic permeability ratio (default: 1.0e-6).
- `kappa_frac::Real`: Hydrofracture multiplier (default: 1.0e3).
- `gamma_frac::Real`: Hydrofracture power-law exponent (default: 1.0).
- `k_frac_max::Real`: Maximum fractured permeability ceiling [m^2] (default: 1.0e-9).

# Returns
- Effective face permeability [m^2].
"""
function compute_face_venting_permeability(
    k_v::Real,
    breached::Bool,
    ice_sealing::Bool,
    T_surf::Real,
    peff::Real,
    sigma_t::Real;
    species::Symbol=:H2O,
    t_freeze::Real=273.15,
    dt_seal::Real=10.0,
    k_seal_min_ratio::Real=1.0e-6,
    kappa_frac::Real=1.0e3,
    gamma_frac::Real=1.0,
    k_frac_max::Real=1.0e-9,
)
    if breached
        return compute_hydrofracture_permeability(
            k_v,
            peff,
            sigma_t;
            active=true,
            kappa_frac=kappa_frac,
            gamma=gamma_frac,
            kmax=k_frac_max,
        )
    elseif ice_sealing && (species === :H2O || species === :water)
        return compute_ice_sealed_permeability(
            k_v,
            T_surf;
            T_freeze=t_freeze,
            delta_T_seal=dt_seal,
            k_min_ratio=k_seal_min_ratio,
        )
    else
        return Float64(k_v)
    end
end

"""
    apply_venting_surface_boundary!(
        L, R, tk, coords, rplanet, xcenter, ycenter, P_amb;
        k_vent=1.0e-11, conductance_factor=1.0, mode=:darcy_sink,
        hydrofracture=false,
        ice_sealing=false, t_freeze=273.15, dt_seal=10.0, k_seal_min_ratio=1.0e-6,
        kappa_frac=1.0e3, gamma_frac=1.0, k_frac_max=1.0e-9,
        pr=nothing, pf=nothing, TEN=nothing, PHI=nothing, phimin=1.0e-4, dt=1.0e10,
        eta_fluid_surf=1.0e-3, L_sub=2.83e6, Kcont=1.0e20, S_vent_out=nothing
    )

Apply permeable venting sink boundary condition at rock-air interface faces (`r = rplanet`).
Computes local venting pressure `P_vent = max(P_amb, P_sat,ice(T_surf))` and assembles Robin
conductance into the fluid continuity row (scaled by `Kcont`).

If `mode === :hydrofracture_gated`, venting requires `pr`, `pf`, and `TEN` to evaluate tensile
failure `Peff <= -sigma_t`. If any pressure array is missing or the lid is unbreached, the face
remains closed (`is_open = false`).
If `ice_sealing === true`, sub-freezing rock faces (`T_surf < t_freeze`) experience exponential
pore ice permeability sealing during unbreached porous flow (`:darcy_sink` mode).
When overpressure breaches the lid (`Peff <= -sigma_t` and `hydrofracture === true` or
`mode === :hydrofracture_gated`), enhanced hydrofracture permeability opens.
Only outward venting is permitted (`pf > P_vent`), and venting is fluid-limited (`phi > phimin`).
"""
function apply_venting_surface_boundary!(
    L,
    R::Union{AbstractVector{Float64},Nothing},
    tk::AbstractMatrix{Float64},
    coords::GridCoordinates,
    rplanet::Real,
    xcenter::Real,
    ycenter::Real,
    P_amb::Real;
    species::Symbol=:H2O,
    k_vent::Real=1.0e-11,
    conductance_factor::Real=1.0,
    mode::Symbol=:darcy_sink,
    hydrofracture::Bool=false,
    ice_sealing::Bool=false,
    t_freeze::Real=273.15,
    dt_seal::Real=10.0,
    k_seal_min_ratio::Real=1.0e-6,
    kappa_frac::Real=1.0e3,
    gamma_frac::Real=1.0,
    k_frac_max::Real=1.0e-9,
    pr::Union{AbstractMatrix{Float64},Nothing}=nothing,
    pf::Union{AbstractMatrix{Float64},Nothing}=nothing,
    TEN::Union{AbstractMatrix{Float64},Nothing}=nothing,
    PHI::Union{AbstractMatrix{Float64},Nothing}=nothing,
    phimin::Real=1.0e-4,
    dt::Real=1.0e10,
    eta_fluid_surf::Real=1.0e-3,
    L_sub::Real=2.83e6,
    Kcont::Real=1.0e20,
    S_vent_out::Union{AbstractMatrix{Float64},Nothing}=nothing,
    dof_stride::Int=6,
)
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    dx = coords.dx
    dy = coords.dy
    if dx <= 0.0 || dy <= 0.0 || !isfinite(dx) || !isfinite(dy)
        throw(DomainError((dx, dy), "Grid spacing must be > 0 and finite"))
    end
    eta_f = Float64(eta_fluid_surf)
    if eta_f <= 0.0 || !isfinite(eta_f)
        throw(DomainError(eta_f, "eta_fluid_surf must be > 0 and finite"))
    end
    rplanet2 = Float64(rplanet)^2
    p_amb_val = Float64(P_amb)
    k_v = Float64(k_vent)
    c_factor = Float64(conductance_factor)
    kcont_val = Float64(Kcont)
    dt_val = max(Float64(dt), 1.0e-12)
    phimin_val = Float64(phimin)

    if S_vent_out !== nothing
        S_vent_out .= 0.0
    end

    # Horizontal faces between P(i, j) and P(i, j+1)
    @inbounds for j in 1:(Nx1 - 1)
        xj1 = coords.xp[j] - xcenter
        xj2 = coords.xp[j + 1] - xcenter
        for i in 1:Ny1
            yi = coords.yp[i] - ycenter
            r1_sq = xj1^2 + yi^2
            r2_sq = xj2^2 + yi^2
            is_rock1 = r1_sq <= rplanet2
            is_rock2 = r2_sq <= rplanet2
            if is_rock1 != is_rock2
                i_rock = i
                j_rock = is_rock1 ? j : (j + 1)

                # Skip domain boundary ghost and anchor nodes
                if i_rock < 2 || i_rock > Ny1 - 1 || j_rock < 2 || j_rock > Nx1 - 1
                    continue
                end

                breached = false
                peff = 0.0
                sigma_t = 0.0
                if (mode === :hydrofracture_gated || hydrofracture) &&
                    pr !== nothing &&
                    pf !== nothing &&
                    TEN !== nothing
                    peff = pr[i_rock, j_rock] - pf[i_rock, j_rock]
                    sigma_t = TEN[i_rock, j_rock]
                    breached = is_hydrofracture_breached(peff, sigma_t)
                end

                is_open = true
                if mode === :hydrofracture_gated && !breached
                    is_open = false
                end

                T_raw = tk[i_rock, j_rock]
                T_surf = isfinite(T_raw) ? max(T_raw, 1.0e-3) : 1.0e-3
                P_vent = compute_venting_pressure(
                    T_surf, p_amb_val; species=species, L_sub=L_sub
                )

                # One-sided venting condition: pore fluid must exceed venting pressure
                if pf !== nothing
                    pf_cur = pf[i_rock, j_rock]
                    if pf_cur <= P_vent
                        is_open = false
                    end
                end

                # Fluid-availability limit: no venting from dry rock
                phi_avail = 1.0
                if PHI !== nothing
                    phi_rock = PHI[i_rock, j_rock]
                    phi_avail = max(0.0, phi_rock - phimin_val)
                    if phi_avail <= 0.0
                        is_open = false
                    end
                end

                if is_open
                    k_face = compute_face_venting_permeability(
                        k_v,
                        breached,
                        ice_sealing,
                        T_surf,
                        peff,
                        sigma_t;
                        species=species,
                        t_freeze=t_freeze,
                        dt_seal=dt_seal,
                        k_seal_min_ratio=k_seal_min_ratio,
                        kappa_frac=kappa_frac,
                        gamma_frac=gamma_frac,
                        k_frac_max=k_frac_max,
                    )

                    C_face = (k_face / (eta_f * dx^2)) * c_factor
                    kpf = ((j_rock - 1) * Ny1 + i_rock - 1) * dof_stride + dof_stride

                    # If pf and PHI are known, check and cap Darcy rate by available fluid
                    C_face_eff = C_face
                    S_vent_actual = 0.0
                    if pf !== nothing
                        pf_cur = pf[i_rock, j_rock]
                        S_darcy = C_face * (pf_cur - P_vent)
                        if PHI !== nothing
                            S_max = phi_avail / dt_val
                            if S_darcy > S_max && S_darcy > 0.0
                                C_face_eff = C_face * (S_max / S_darcy)
                                S_vent_actual = S_max
                            else
                                S_vent_actual = max(0.0, S_darcy)
                            end
                        else
                            S_vent_actual = max(0.0, S_darcy)
                        end
                    end

                    if L !== nothing && R !== nothing
                        updateindex!(L, +, kcont_val * C_face_eff, kpf, kpf)
                        R[kpf] += C_face_eff * P_vent
                    end

                    if S_vent_out !== nothing
                        S_vent_out[i_rock, j_rock] += S_vent_actual
                    end
                end
            end
        end
    end

    # Vertical faces between P(i, j) and P(i+1, j)
    @inbounds for j in 1:Nx1
        xj = coords.xp[j] - xcenter
        for i in 1:(Ny1 - 1)
            yi1 = coords.yp[i] - ycenter
            yi2 = coords.yp[i + 1] - ycenter
            r1_sq = xj^2 + yi1^2
            r2_sq = xj^2 + yi2^2
            is_rock1 = r1_sq <= rplanet2
            is_rock2 = r2_sq <= rplanet2
            if is_rock1 != is_rock2
                i_rock = is_rock1 ? i : (i + 1)
                j_rock = j

                # Skip domain boundary ghost and anchor nodes
                if i_rock < 2 || i_rock > Ny1 - 1 || j_rock < 2 || j_rock > Nx1 - 1
                    continue
                end

                breached = false
                peff = 0.0
                sigma_t = 0.0
                if (mode === :hydrofracture_gated || hydrofracture) &&
                    pr !== nothing &&
                    pf !== nothing &&
                    TEN !== nothing
                    peff = pr[i_rock, j_rock] - pf[i_rock, j_rock]
                    sigma_t = TEN[i_rock, j_rock]
                    breached = is_hydrofracture_breached(peff, sigma_t)
                end

                is_open = true
                if mode === :hydrofracture_gated && !breached
                    is_open = false
                end

                T_raw = tk[i_rock, j_rock]
                T_surf = isfinite(T_raw) ? max(T_raw, 1.0e-3) : 1.0e-3
                P_vent = compute_venting_pressure(
                    T_surf, p_amb_val; species=species, L_sub=L_sub
                )

                # One-sided venting condition: pore fluid must exceed venting pressure
                if pf !== nothing
                    pf_cur = pf[i_rock, j_rock]
                    if pf_cur <= P_vent
                        is_open = false
                    end
                end

                # Fluid-availability limit: no venting from dry rock
                phi_avail = 1.0
                if PHI !== nothing
                    phi_rock = PHI[i_rock, j_rock]
                    phi_avail = max(0.0, phi_rock - phimin_val)
                    if phi_avail <= 0.0
                        is_open = false
                    end
                end

                if is_open
                    k_face = compute_face_venting_permeability(
                        k_v,
                        breached,
                        ice_sealing,
                        T_surf,
                        peff,
                        sigma_t;
                        species=species,
                        t_freeze=t_freeze,
                        dt_seal=dt_seal,
                        k_seal_min_ratio=k_seal_min_ratio,
                        kappa_frac=kappa_frac,
                        gamma_frac=gamma_frac,
                        k_frac_max=k_frac_max,
                    )

                    C_face = (k_face / (eta_f * dy^2)) * c_factor
                    kpf = ((j_rock - 1) * Ny1 + i_rock - 1) * dof_stride + dof_stride

                    # If pf and PHI are known, check and cap Darcy rate by available fluid
                    C_face_eff = C_face
                    S_vent_actual = 0.0
                    if pf !== nothing
                        pf_cur = pf[i_rock, j_rock]
                        S_darcy = C_face * (pf_cur - P_vent)
                        if PHI !== nothing
                            S_max = phi_avail / dt_val
                            if S_darcy > S_max && S_darcy > 0.0
                                C_face_eff = C_face * (S_max / S_darcy)
                                S_vent_actual = S_max
                            else
                                S_vent_actual = max(0.0, S_darcy)
                            end
                        else
                            S_vent_actual = max(0.0, S_darcy)
                        end
                    end

                    if L !== nothing && R !== nothing
                        updateindex!(L, +, kcont_val * C_face_eff, kpf, kpf)
                        R[kpf] += C_face_eff * P_vent
                    end

                    if S_vent_out !== nothing
                        S_vent_out[i_rock, j_rock] += S_vent_actual
                    end
                end
            end
        end
    end
    return nothing
end

"""
    compute_mean_surface_temperature(
        tk::AbstractMatrix{<:Real},
        coords::GridCoordinates,
        rplanet::Real,
        xcenter::Real,
        ycenter::Real;
        T_default::Real=300.0,
    )::Float64

Compute average temperature of planetesimal rock nodes at the surface boundary (r ≈ rplanet).
Returns `T_default` if no surface rock nodes are identified.
"""
function compute_mean_surface_temperature(
    tk::AbstractMatrix{<:Real},
    coords::GridCoordinates,
    rplanet::Real,
    xcenter::Real,
    ycenter::Real;
    T_default::Real=300.0,
)::Float64
    Ny1, Nx1 = coords.Ny1, coords.Nx1
    rplanet2 = Float64(rplanet)^2
    t_sum = 0.0
    count = 0
    visited = falses(Ny1, Nx1)
    # Horizontal rock-air interfaces
    @inbounds for j in 1:(Nx1 - 1)
        xj1 = coords.xp[j] - xcenter
        xj2 = coords.xp[j + 1] - xcenter
        for i in 1:Ny1
            yi = coords.yp[i] - ycenter
            is_rock1 = (xj1^2 + yi^2) <= rplanet2
            is_rock2 = (xj2^2 + yi^2) <= rplanet2
            if is_rock1 != is_rock2
                i_rock = i
                j_rock = is_rock1 ? j : (j + 1)
                if 2 <= i_rock <= Ny1 - 1 &&
                    2 <= j_rock <= Nx1 - 1 &&
                    !visited[i_rock, j_rock]
                    t_val = Float64(tk[i_rock, j_rock])
                    if isfinite(t_val) && t_val > 0.0
                        visited[i_rock, j_rock] = true
                        t_sum += t_val
                        count += 1
                    end
                end
            end
        end
    end
    # Vertical rock-air interfaces
    @inbounds for j in 1:Nx1
        xj = coords.xp[j] - xcenter
        for i in 1:(Ny1 - 1)
            yi1 = coords.yp[i] - ycenter
            yi2 = coords.yp[i + 1] - ycenter
            is_rock1 = (xj^2 + yi1^2) <= rplanet2
            is_rock2 = (xj^2 + yi2^2) <= rplanet2
            if is_rock1 != is_rock2
                i_rock = is_rock1 ? i : (i + 1)
                j_rock = j
                if 2 <= i_rock <= Ny1 - 1 &&
                    2 <= j_rock <= Nx1 - 1 &&
                    !visited[i_rock, j_rock]
                    t_val = Float64(tk[i_rock, j_rock])
                    if isfinite(t_val) && t_val > 0.0
                        visited[i_rock, j_rock] = true
                        t_sum += t_val
                        count += 1
                    end
                end
            end
        end
    end
    return count > 0 ? t_sum / count : Float64(T_default)
end

"""
Perform thermal iterations to time step thermal field at P nodes.

$(SIGNATURES)

# Details

    - tk0: previous temperature at P nodes 
	- tk1: current temperature at P nodes
	- tk2: next temperature at P nodes 
	- DT: calculated temperature difference at P nodes 
	- DT0: previous calculated temperature difference at P nodes
	- RHOCP: volumetric heat capacity at P nodes  
	- KX: thermal conductivity at Vx nodes
	- KY: thermal conductivity at Vy nodes 
	- HR: radioactive heating at P nodes
	- HA: adiabatic heating at P nodes 
	- HS: shear heating at P nodes
    - DHP: latent heating (HL) at P nodes
    - RT: thermal RHS coefficient vector
    - ST: thermal solution vector
	- dt: computational time step

# Returns

    - nothing
"""
function perform_thermal_iterations!(
    tk0,
    tk1,
    tk2,
    DT,
    DT0,
    RHOCP,
    KX,
    KY,
    HR,
    HA,
    HS,
    DHP,
    RT,
    ST,
    dt;
    coords=nothing,
    Q_metric=nothing,
    Q_lat=nothing,
    DTmax_val::Real=20.0,
)
    # set up thermal iterations
    Ny1, Nx1 = size(tk1)
    tk0 .= tk1
    dtt = dt
    dttsum = 0.0
    titer = 1
    # perform thermal iterations until reaching time limit
    while dttsum < dt
        # fresh LHS coefficient matrix
        LT = assemble_thermal_lse!(
            tk1,
            RHOCP,
            KX,
            KY,
            HR,
            HA,
            HS,
            DHP,
            RT,
            dtt;
            coords=coords,
            Q_metric=Q_metric,
            Q_lat=Q_lat,
        )

        # solve system of equations
        ST .= LT \ RT # implicit: flush!(LT)
        # reshape solution vector to 2D array
        tk2 .= reshape(ST, Ny1, Nx1)
        # compute ΔT
        DT .= tk2 .- tk1
        maxDTcurrent = maximum(abs, DT)
        if maxDTcurrent > DTmax_val
            dtt *= DTmax_val * inv(maxDTcurrent)
        else
            dttsum += dtt
            tk1 .= tk2
            if dttsum >= dt || (dt - dttsum) <= 1e-14 * dt
                break
            end
            dtt = min(dtt, dt - dttsum)
        end
        # increase thermal iteration counter
        titer += 1
    end
    # finalize overall temperature change and advance temperature field
    DT .= tk2 .- tk0
    DT0 .= DT
    return nothing
end # function perform_thermal_iterations!

"""
Decide next pass thermochemical iteration time step.

$(SIGNATURES)

# Details:

    - maxDTcurrent: maximum temperature difference between current and
                    previous time step
    - dt: current time step duration
    - titer: current thermochemical iteration counter

# Returns

    - dt: adjusted next time step
"""
function finalize_thermochemical_iteration_pass(
    maxDTcurrent, dt, titer, DTmax_val::Real=20.0
)
    if maxDTcurrent > DTmax_val
        dt *= (DTmax_val * inv(maxDTcurrent))
        @info "titer $titer: reducing dt due to maxDT: dt=$dt s"
    end
    return dt
end # function finalize_thermochemical_iteration_pass

"""
Assess outcome of thermochemical iteration and return thermochemical iterations
completeness status.

$(SIGNATURES)

# Details:

    - DMP: mass transfer term at P nodes
    - pf: fluid pressure at P nodes
    - pf0: previous time step fluid pressure at P nodes
    - titer: current thermochemical iteration counter

# Returns

    - dt: adjusted next time step
"""
function compute_thermochemical_iteration_outcome(DMP, pf, pf0, titer; pferrmax=1.0e5)
    pferrcur = maximum(abs, pf - pf0)
    DMPmax = maximum(abs, DMP)
    @info "end thermochemical iter $titer" pferrcur DMPmax
    return pferrcur < pferrmax && (titer > 2 || DMPmax <= 0.0)
end # function compute_thermochemical_iteration_outcome

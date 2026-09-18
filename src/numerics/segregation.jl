"""
Apply iron metal segregation via sub-cycled conservative drift-flux transport.

$(SIGNATURES)

Solves the conservative drift-flux transport equation for molten iron metal
percolation through a solid silicate matrix and Stokes settling through a magma ocean.
Subcycles the explicit finite-volume transport step using a local CFL criterion.
Guarantees mass conservation of total metal to machine precision when input marker
metal fractions satisfy 0 <= Xfe_bulk <= phi_pack.

References:
- Stevenson (1990), Fluid dynamics of core formation.
- Deguen et al. (2014), Earth Planet. Sci. Lett., 391, 274-287.
- Lichtenberg et al. (2019, 2021), Science / JGR Planets.

# Arguments
- `xm::AbstractVector{Float64}`: Marker x-coordinates [m]
- `ym::AbstractVector{Float64}`: Marker y-coordinates [m]
- `tm::AbstractVector{<:Integer}`: Marker material phase type
- `tkm::AbstractVector{Float64}`: Marker temperature [K]
- `phim::AbstractVector{Float64}`: Marker silicate melt fraction / porosity [-]
- `Xfe_bulk::AbstractVector{Float64}`: Marker bulk metal volume fraction [-]
- `Xfem::AbstractVector{Float64}`: Marker molten metal fraction [-]
- `marknum::Integer`: Number of markers
- `dt::Real`: Timestep duration [s]
- `cfg_core::CoreFormationConfig`: Core formation configuration parameters

# Keyword Arguments
- `coords=nothing`: `GridCoordinates` domain geometry struct
- `xcenter::Real=coords !== nothing ? coords.xcenter : 70000.0`: Planet center x [m]
- `ycenter::Real=coords !== nothing ? coords.ycenter : 70000.0`: Planet center y [m]
- `rplanet::Real=50000.0`: Planet radius [m]
- `g_surf::Real=0.1`: Reference surface gravity magnitude [m/s^2]
- `gx::Union{Nothing,AbstractMatrix{Float64}}=nothing`: Optional x-gravity on grid [m/s^2]
- `gy::Union{Nothing,AbstractMatrix{Float64}}=nothing`: Optional y-gravity on grid [m/s^2]
- `Q_seg_grid::Union{Nothing,AbstractMatrix{Float64}}=nothing`: Optional grid to accumulate dissipation heating [W/m^3]
- `rho_silicate::Real=3000.0`: Reference silicate rock density [kg/m^3]
- `eta_silicate::Real=1.0e18`: Reference silicate rock dynamic viscosity [Pa s]
- `ETA::Union{Nothing,AbstractMatrix{Float64}}=nothing`: Matrix viscosity on grid [Pa s]
- `Fm::Union{Nothing,AbstractVector{Float64}}=nothing`: Marker silicate melt fraction [-]
- `T_solidus_silicate::Real=1400.0`: Reference silicate solidus temperature [K]
- `T_liquidus_silicate::Real=1800.0`: Reference silicate liquidus temperature [K]

# Returns
- NamedTuple `(; max_v_seg, n_subcycles, dt_sub, total_dissipation_energy)` where `total_dissipation_energy` is in [J/m].
"""
function apply_metal_segregation!(
    xm::AbstractVector{Float64},
    ym::AbstractVector{Float64},
    tm::AbstractVector{<:Integer},
    tkm::AbstractVector{Float64},
    phim::AbstractVector{Float64},
    Xfe_bulk::AbstractVector{Float64},
    Xfem::AbstractVector{Float64},
    marknum::Integer,
    dt::Real,
    cfg_core::CoreFormationConfig;
    coords=nothing,
    xcenter::Real=coords !== nothing ? coords.xcenter : 70000.0,
    ycenter::Real=coords !== nothing ? coords.ycenter : 70000.0,
    rplanet::Real=50000.0,
    g_surf::Real=0.1,
    gx::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    gy::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    Q_seg_grid::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    rho_silicate::Real=3000.0,
    eta_silicate::Real=1.0e18,
    ETA::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    Fm::Union{Nothing,AbstractVector{Float64}}=nothing,
    T_solidus_silicate::Real=1400.0,
    T_liquidus_silicate::Real=1800.0,
    Xfe_H_m::Union{Nothing,AbstractVector{Float64}}=nothing,
    Xfe_C_m::Union{Nothing,AbstractVector{Float64}}=nothing,
    Xfe_N_m::Union{Nothing,AbstractVector{Float64}}=nothing,
    Xfe_S_m::Union{Nothing,AbstractVector{Float64}}=nothing,
    cfg_partition::Union{Nothing,MetalPartitionConfig}=nothing,
    workspace=nothing,
)
    if (!cfg_core.percolation_active && !cfg_core.settling_active) ||
        dt <= 0.0 ||
        marknum <= 0
        return (; max_v_seg=0.0, n_subcycles=0, dt_sub=0.0, total_dissipation_energy=0.0)
    end

    track_volatiles = (
        cfg_partition !== nothing &&
        cfg_partition.active &&
        Xfe_H_m !== nothing &&
        Xfe_C_m !== nothing &&
        Xfe_N_m !== nothing &&
        Xfe_S_m !== nothing
    )

    # Validate input marker bounds
    @inbounds for m in 1:marknum
        if tm[m] < 3
            rmark = distance(xm[m], ym[m], xcenter, ycenter)
            if rmark <= rplanet
                xfe = Xfe_bulk[m]
                if !isfinite(xfe) || xfe < 0.0 || xfe > cfg_core.phi_pack + 1.0e-7
                    throw(
                        DomainError(
                            xfe,
                            "Marker bulk metal fraction must be finite, non-negative, and <= phi_pack",
                        ),
                    )
                end
            end
        end
    end

    Nx_val = coords !== nothing ? coords.Nx : 32
    Ny_val = coords !== nothing ? coords.Ny : 32
    dx_val = if coords !== nothing
        coords.dx
    else
        (coords !== nothing ? coords.xsize / Nx_val : 4375.0)
    end
    dy_val = if coords !== nothing
        coords.dy
    else
        (coords !== nothing ? coords.ysize / Ny_val : 4375.0)
    end

    has_ws = (
        workspace !== nothing &&
        size(workspace.M_fe_cell) == (Ny_val, Nx_val) &&
        (!track_volatiles || size(workspace.M_fe_H_cell) == (Ny_val, Nx_val))
    )

    # Allocate or reset cell accumulations
    M_fe_cell = has_ws ? fill!(workspace.M_fe_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)
    M_rock_markers =
        has_ws ? fill!(workspace.M_rock_markers, 0) : zeros(Int, Ny_val, Nx_val)
    v_seg_cell = has_ws ? fill!(workspace.v_seg_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)
    phi_m_cell = has_ws ? fill!(workspace.phi_m_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)
    F_m_cell = has_ws ? fill!(workspace.F_m_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)
    Xfem_cell = has_ws ? fill!(workspace.Xfem_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)
    g_acc_cell = has_ws ? fill!(workspace.g_acc_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)
    cap_cell = has_ws ? fill!(workspace.cap_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)
    phi_fe_cell =
        has_ws ? fill!(workspace.phi_fe_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)
    T_cell = has_ws ? fill!(workspace.T_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)
    drho_cell = has_ws ? fill!(workspace.drho_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)

    M_fe_H_cell = if track_volatiles
        has_ws ? fill!(workspace.M_fe_H_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)
    else
        zeros(Float64, 0, 0)
    end
    M_fe_C_cell = if track_volatiles
        has_ws ? fill!(workspace.M_fe_C_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)
    else
        zeros(Float64, 0, 0)
    end
    M_fe_N_cell = if track_volatiles
        has_ws ? fill!(workspace.M_fe_N_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)
    else
        zeros(Float64, 0, 0)
    end
    M_fe_S_cell = if track_volatiles
        has_ws ? fill!(workspace.M_fe_S_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)
    else
        zeros(Float64, 0, 0)
    end

    # Bin markers into grid cells
    @inbounds for m in 1:marknum
        if tm[m] < 3
            rmark = distance(xm[m], ym[m], xcenter, ycenter)
            if rmark <= rplanet
                j_c = clamp(Int(floor(xm[m] / dx_val)) + 1, 1, Nx_val)
                i_c = clamp(Int(floor(ym[m] / dy_val)) + 1, 1, Ny_val)
                fe_m = Xfe_bulk[m]
                M_fe_cell[i_c, j_c] += fe_m
                M_rock_markers[i_c, j_c] += 1
                phi_m_cell[i_c, j_c] += Xfem[m]
                T_cell[i_c, j_c] += tkm[m]
                F_m_val = if Fm !== nothing
                    Fm[m]
                elseif tkm[m] >= T_solidus_silicate &&
                    T_liquidus_silicate > T_solidus_silicate
                    clamp(
                        (tkm[m] - T_solidus_silicate) /
                        (T_liquidus_silicate - T_solidus_silicate),
                        0.0,
                        1.0,
                    )
                else
                    0.0
                end
                F_m_cell[i_c, j_c] += F_m_val
                cap_cell[i_c, j_c] += max(cfg_core.phi_pack - fe_m, 0.0)

                if track_volatiles
                    M_fe_H_cell[i_c, j_c] += fe_m * Xfe_H_m[m]
                    M_fe_C_cell[i_c, j_c] += fe_m * Xfe_C_m[m]
                    M_fe_N_cell[i_c, j_c] += fe_m * Xfe_N_m[m]
                    M_fe_S_cell[i_c, j_c] += fe_m * Xfe_S_m[m]
                end
            end
        end
    end

    # Normalize cell averages
    @inbounds for j in 1:Nx_val, i in 1:Ny_val
        n_m = M_rock_markers[i, j]
        if n_m > 0
            phi_fe_cell[i, j] = M_fe_cell[i, j] / n_m
            phi_m_cell[i, j] /= n_m
            F_m_cell[i, j] /= n_m
            T_cell[i, j] /= n_m
            m_bulk = phi_fe_cell[i, j]
            Xfem_cell[i, j] =
                m_bulk > 0.0 ? clamp(phi_m_cell[i, j] / m_bulk, 0.0, 1.0) : 0.0
        end
    end

    # Compute cell segregation velocities
    @inbounds for j in 1:Nx_val, i in 1:Ny_val
        n_m = M_rock_markers[i, j]
        if n_m == 0
            continue
        end
        xc = (j - 0.5) * dx_val
        yc = (i - 0.5) * dy_val
        rc = distance(xc, yc, xcenter, ycenter)
        if rc > rplanet
            continue
        end

        g_acc = if gx !== nothing && gy !== nothing && i <= size(gx, 1) && j <= size(gx, 2)
            g_mag = sqrt(gx[i, j]^2 + gy[i, j]^2)
            g_mag > 0.0 ? g_mag : g_surf * min(rc / rplanet, 1.0)
        else
            g_surf * min(rc / rplanet, 1.0)
        end
        g_acc_cell[i, j] = g_acc

        phi_m = phi_m_cell[i, j]
        F_m = F_m_cell[i, j]

        rho_metal_eff = if cfg_core.metal_density_mode !== :constant
            w_S_val =
                if track_volatiles &&
                    cfg_partition.dynamic_sulfur_density &&
                    M_fe_cell[i, j] > 0.0
                    clamp((M_fe_S_cell[i, j] / M_fe_cell[i, j]) * 1.0e-6, 0.0, 0.40)
                else
                    cfg_core.sulfur_fraction
                end
            compute_liquid_metal_density(
                w_S_val; T=max(T_cell[i, j], 100.0), law=cfg_core.metal_density_mode
            )
        else
            cfg_core.rho_metal
        end
        drho = max(rho_metal_eff - rho_silicate, 1.0)
        drho_cell[i, j] = drho

        if phi_m > 0.0 && g_acc > 0.0 && drho > 0.0
            eta_matrix = if ETA !== nothing && i <= size(ETA, 1) && j <= size(ETA, 2)
                ETA[i, j]
            else
                eta_silicate
            end
            eta_susp = compute_melt_weakened_viscosity(
                eta_matrix, F_m, 1; phi_crit=0.4, eta_melt=10.0, etamin=0.1, etamax=1.0e20
            )
            r_drop = if cfg_core.droplet_size_mode === :fixed
                cfg_core.droplet_diameter_fixed / 2.0
            elseif cfg_core.droplet_size_mode === :capillary_mean ||
                cfg_core.droplet_size_mode === :bond_mean ||
                cfg_core.droplet_size_mode === :weber_mean
                # Gravity-capillary (Bond) balance: d = sqrt(We_crit * sigma / (drho * g))
                d_cap = sqrt(
                    cfg_core.We_crit * cfg_core.sigma_metal_silicate /
                    max(drho * g_acc, 1.0e-8),
                )
                clamp(d_cap / 2.0, 1.0e-4, 5.0e-2)
            else # :weber_turbulent
                v_est = stokes_settling_velocity(
                    cfg_core.droplet_diameter_fixed / 2.0,
                    drho,
                    max(g_acc, 1.0e-5),
                    eta_susp,
                )
                v_rel = max(v_est, 1.0e-6)
                d_weber = weber_equilibrium_diameter(
                    rho_silicate,
                    v_rel,
                    cfg_core.sigma_metal_silicate;
                    We_crit=cfg_core.We_crit,
                )
                clamp(d_weber / 2.0, 1.0e-4, 5.0e-2)
            end

            v_seg_cell[i, j] = metal_segregation_velocity(
                phi_m,
                F_m,
                drho,
                g_acc,
                eta_susp;
                percolation_active=cfg_core.percolation_active,
                settling_active=cfg_core.settling_active,
                k_metal_ref=cfg_core.k_metal_ref,
                eta_metal=cfg_core.eta_metal,
                phi_crit_perc=cfg_core.phi_crit_perc,
                phi_residual=cfg_core.phi_residual,
                phi0=cfg_core.phi0,
                perm_exponent=cfg_core.perm_exponent,
                r_drop=r_drop,
                hindered_exponent=cfg_core.hindered_exponent,
                phi_pack=cfg_core.phi_pack,
                hadamard_rybczynski=cfg_core.hadamard_rybczynski,
                F_settle_start=cfg_core.F_settle_start,
                F_perc_end=cfg_core.F_perc_end,
            )
        end
    end

    max_v = maximum(v_seg_cell)
    if max_v <= 0.0
        return (; max_v_seg=0.0, n_subcycles=0, dt_sub=0.0, total_dissipation_energy=0.0)
    end

    # CFL calculation and subcycling
    dt_cfl = cfg_core.cfl_settling * min(dx_val, dy_val) / max_v
    n_sub_raw = Int(ceil(dt / dt_cfl))
    if n_sub_raw > cfg_core.max_subcycles
        @warn "CFL subcycling requires $n_sub_raw steps, capped at max_subcycles $(cfg_core.max_subcycles); segregation transport throttled" maxlog=10
    end
    n_sub = clamp(n_sub_raw, 1, cfg_core.max_subcycles)
    dt_sub = dt / n_sub

    # Working copy of cell metal mass for subcycling
    m_fe = has_ws ? copyto!(workspace.m_fe, M_fe_cell) : copy(M_fe_cell)
    total_diss_energy = 0.0

    m_fe_H = if track_volatiles
        has_ws ? copyto!(workspace.m_fe_H, M_fe_H_cell) : copy(M_fe_H_cell)
    else
        zeros(Float64, 0, 0)
    end
    m_fe_C = if track_volatiles
        has_ws ? copyto!(workspace.m_fe_C, M_fe_C_cell) : copy(M_fe_C_cell)
    else
        zeros(Float64, 0, 0)
    end
    m_fe_N = if track_volatiles
        has_ws ? copyto!(workspace.m_fe_N, M_fe_N_cell) : copy(M_fe_N_cell)
    else
        zeros(Float64, 0, 0)
    end
    m_fe_S = if track_volatiles
        has_ws ? copyto!(workspace.m_fe_S, M_fe_S_cell) : copy(M_fe_S_cell)
    else
        zeros(Float64, 0, 0)
    end

    flux_H_x = if track_volatiles
        has_ws ? fill!(workspace.flux_H_x, 0.0) : zeros(Float64, Ny_val, Nx_val - 1)
    else
        zeros(Float64, 0, 0)
    end
    flux_H_y = if track_volatiles
        has_ws ? fill!(workspace.flux_H_y, 0.0) : zeros(Float64, Ny_val - 1, Nx_val)
    else
        zeros(Float64, 0, 0)
    end
    flux_C_x = if track_volatiles
        has_ws ? fill!(workspace.flux_C_x, 0.0) : zeros(Float64, Ny_val, Nx_val - 1)
    else
        zeros(Float64, 0, 0)
    end
    flux_C_y = if track_volatiles
        has_ws ? fill!(workspace.flux_C_y, 0.0) : zeros(Float64, Ny_val - 1, Nx_val)
    else
        zeros(Float64, 0, 0)
    end
    flux_N_x = if track_volatiles
        has_ws ? fill!(workspace.flux_N_x, 0.0) : zeros(Float64, Ny_val, Nx_val - 1)
    else
        zeros(Float64, 0, 0)
    end
    flux_N_y = if track_volatiles
        has_ws ? fill!(workspace.flux_N_y, 0.0) : zeros(Float64, Ny_val - 1, Nx_val)
    else
        zeros(Float64, 0, 0)
    end
    flux_S_x = if track_volatiles
        has_ws ? fill!(workspace.flux_S_x, 0.0) : zeros(Float64, Ny_val, Nx_val - 1)
    else
        zeros(Float64, 0, 0)
    end
    flux_S_y = if track_volatiles
        has_ws ? fill!(workspace.flux_S_y, 0.0) : zeros(Float64, Ny_val - 1, Nx_val)
    else
        zeros(Float64, 0, 0)
    end

    # Pre-allocated arrays for subcycling fluxes and limiters
    req_flux_x =
        has_ws ? fill!(workspace.req_flux_x, 0.0) : zeros(Float64, Ny_val, Nx_val - 1)
    req_flux_y =
        has_ws ? fill!(workspace.req_flux_y, 0.0) : zeros(Float64, Ny_val - 1, Nx_val)
    flux_x = has_ws ? fill!(workspace.flux_x, 0.0) : zeros(Float64, Ny_val, Nx_val - 1)
    flux_y = has_ws ? fill!(workspace.flux_y, 0.0) : zeros(Float64, Ny_val - 1, Nx_val)
    outflow_tot =
        has_ws ? fill!(workspace.outflow_tot, 0.0) : zeros(Float64, Ny_val, Nx_val)
    inflow_tot = has_ws ? fill!(workspace.inflow_tot, 0.0) : zeros(Float64, Ny_val, Nx_val)
    alpha_out = has_ws ? fill!(workspace.alpha_out, 1.0) : ones(Float64, Ny_val, Nx_val)
    alpha_in = has_ws ? fill!(workspace.alpha_in, 1.0) : ones(Float64, Ny_val, Nx_val)

    # Subcycling loop
    for _ in 1:n_sub
        fill!(outflow_tot, 0.0)
        fill!(inflow_tot, 0.0)
        fill!(req_flux_x, 0.0)
        fill!(req_flux_y, 0.0)

        # 1. Compute unscaled requested fluxes across East-West faces
        @inbounds for j in 1:(Nx_val - 1)
            xf = j * dx_val
            for i in 1:Ny_val
                yf = (i - 0.5) * dy_val
                dxf = xf - xcenter
                dyf = yf - ycenter
                rf = sqrt(dxf^2 + dyf^2)
                if rf > rplanet || rf < 1.0e-3
                    continue
                end

                # Determine transport direction
                nx =
                    if gx !== nothing &&
                        gy !== nothing &&
                        i <= size(gx, 1) &&
                        j <= size(gx, 2)
                        gx_f = gx[i, j]
                        gy_f =
                            0.5 *
                            (gy[i, j] + (j + 1 <= size(gy, 2) ? gy[i, j + 1] : gy[i, j]))
                        g_f = sqrt(gx_f^2 + gy_f^2)
                        g_f > 1.0e-10 ? gx_f / g_f : -dxf / rf
                    else
                        -dxf / rf
                    end

                vf = 0.5 * (v_seg_cell[i, j] + v_seg_cell[i, j + 1])
                uf = vf * nx

                donor_j = uf > 0.0 ? j : j + 1
                rec_j = uf > 0.0 ? j + 1 : j

                n_donor = M_rock_markers[i, donor_j]
                n_rec = M_rock_markers[i, rec_j]
                if n_donor == 0 || n_rec == 0
                    continue
                end

                X_donor = m_fe[i, donor_j] / n_donor
                X_mob = max(X_donor - cfg_core.phi_residual, 0.0) * Xfem_cell[i, donor_j]
                m_avail = X_mob * n_donor

                fx = abs(uf) * (dt_sub / dx_val) * m_avail
                fx_req = uf > 0.0 ? fx : -fx
                req_flux_x[i, j] = fx_req

                if fx_req > 0.0
                    outflow_tot[i, j] += fx_req
                    inflow_tot[i, j + 1] += fx_req
                else
                    outflow_tot[i, j + 1] += -fx_req
                    inflow_tot[i, j] += -fx_req
                end
            end
        end

        # 2. Compute unscaled requested fluxes across North-South faces
        @inbounds for i in 1:(Ny_val - 1)
            yf = i * dy_val
            for j in 1:Nx_val
                xf = (j - 0.5) * dx_val
                dxf = xf - xcenter
                dyf = yf - ycenter
                rf = sqrt(dxf^2 + dyf^2)
                if rf > rplanet || rf < 1.0e-3
                    continue
                end

                ny =
                    if gx !== nothing &&
                        gy !== nothing &&
                        i <= size(gy, 1) &&
                        j <= size(gy, 2)
                        gy_f = gy[i, j]
                        gx_f =
                            0.5 *
                            (gx[i, j] + (i + 1 <= size(gx, 1) ? gx[i + 1, j] : gx[i, j]))
                        g_f = sqrt(gx_f^2 + gy_f^2)
                        g_f > 1.0e-10 ? gy_f / g_f : -dyf / rf
                    else
                        -dyf / rf
                    end

                vf = 0.5 * (v_seg_cell[i, j] + v_seg_cell[i + 1, j])
                wf = vf * ny

                donor_i = wf > 0.0 ? i : i + 1
                rec_i = wf > 0.0 ? i + 1 : i

                n_donor = M_rock_markers[donor_i, j]
                n_rec = M_rock_markers[rec_i, j]
                if n_donor == 0 || n_rec == 0
                    continue
                end

                X_donor = m_fe[donor_i, j] / n_donor
                X_mob = max(X_donor - cfg_core.phi_residual, 0.0) * Xfem_cell[donor_i, j]
                m_avail = X_mob * n_donor

                fy = abs(wf) * (dt_sub / dy_val) * m_avail
                fy_req = wf > 0.0 ? fy : -fy
                req_flux_y[i, j] = fy_req

                if fy_req > 0.0
                    outflow_tot[i, j] += fy_req
                    inflow_tot[i + 1, j] += fy_req
                else
                    outflow_tot[i + 1, j] += -fy_req
                    inflow_tot[i, j] += -fy_req
                end
            end
        end

        # 3. Multi-dimensional flux limiters per cell
        @inbounds for j in 1:Nx_val, i in 1:Ny_val
            n_m = M_rock_markers[i, j]
            if n_m > 0
                X_c = m_fe[i, j] / n_m
                m_avail = max(X_c - cfg_core.phi_residual, 0.0) * Xfem_cell[i, j] * n_m
                m_cap = max(cfg_core.phi_pack - X_c, 0.0) * n_m
                alpha_out[i, j] = if outflow_tot[i, j] > m_avail && m_avail > 0.0
                    m_avail / outflow_tot[i, j]
                else
                    (outflow_tot[i, j] > m_avail ? 0.0 : 1.0)
                end
                alpha_in[i, j] = if inflow_tot[i, j] > m_cap && m_cap > 0.0
                    m_cap / inflow_tot[i, j]
                else
                    (inflow_tot[i, j] > m_cap ? 0.0 : 1.0)
                end
            else
                alpha_out[i, j] = 0.0
                alpha_in[i, j] = 0.0
            end
        end

        # 4. Scale fluxes by joint donor-receiver limiters
        @inbounds for j in 1:(Nx_val - 1), i in 1:Ny_val
            fx_req = req_flux_x[i, j]
            if iszero(fx_req)
                flux_x[i, j] = 0.0
            else
                donor_j = fx_req > 0.0 ? j : j + 1
                rec_j = fx_req > 0.0 ? j + 1 : j
                lim = min(alpha_out[i, donor_j], alpha_in[i, rec_j])
                flux_x[i, j] = fx_req * lim
            end
        end

        @inbounds for j in 1:Nx_val, i in 1:(Ny_val - 1)
            fy_req = req_flux_y[i, j]
            if iszero(fy_req)
                flux_y[i, j] = 0.0
            else
                donor_i = fy_req > 0.0 ? i : i + 1
                rec_i = fy_req > 0.0 ? i + 1 : i
                lim = min(alpha_out[donor_i, j], alpha_in[rec_i, j])
                flux_y[i, j] = fy_req * lim
            end
        end

        if track_volatiles
            @inbounds for j in 1:(Nx_val - 1), i in 1:Ny_val
                fx = flux_x[i, j]
                if iszero(fx)
                    flux_H_x[i, j] = 0.0
                    flux_C_x[i, j] = 0.0
                    flux_N_x[i, j] = 0.0
                    flux_S_x[i, j] = 0.0
                else
                    donor_j = fx > 0.0 ? j : j + 1
                    m_d = m_fe[i, donor_j]
                    if m_d > 0.0
                        flux_H_x[i, j] = fx * (m_fe_H[i, donor_j] / m_d)
                        flux_C_x[i, j] = fx * (m_fe_C[i, donor_j] / m_d)
                        flux_N_x[i, j] = fx * (m_fe_N[i, donor_j] / m_d)
                        flux_S_x[i, j] = fx * (m_fe_S[i, donor_j] / m_d)
                    else
                        flux_H_x[i, j] = 0.0
                        flux_C_x[i, j] = 0.0
                        flux_N_x[i, j] = 0.0
                        flux_S_x[i, j] = 0.0
                    end
                end
            end

            @inbounds for j in 1:Nx_val, i in 1:(Ny_val - 1)
                fy = flux_y[i, j]
                if iszero(fy)
                    flux_H_y[i, j] = 0.0
                    flux_C_y[i, j] = 0.0
                    flux_N_y[i, j] = 0.0
                    flux_S_y[i, j] = 0.0
                else
                    donor_i = fy > 0.0 ? i : i + 1
                    m_d = m_fe[donor_i, j]
                    if m_d > 0.0
                        flux_H_y[i, j] = fy * (m_fe_H[donor_i, j] / m_d)
                        flux_C_y[i, j] = fy * (m_fe_C[donor_i, j] / m_d)
                        flux_N_y[i, j] = fy * (m_fe_N[donor_i, j] / m_d)
                        flux_S_y[i, j] = fy * (m_fe_S[donor_i, j] / m_d)
                    else
                        flux_H_y[i, j] = 0.0
                        flux_C_y[i, j] = 0.0
                        flux_N_y[i, j] = 0.0
                        flux_S_y[i, j] = 0.0
                    end
                end
            end
        end

        # 5. Conservative update of cell metal masses
        @inbounds for j in 1:Nx_val, i in 1:Ny_val
            F_w = (j > 1) ? flux_x[i, j - 1] : 0.0
            F_e = (j < Nx_val) ? flux_x[i, j] : 0.0
            F_n = (i > 1) ? flux_y[i - 1, j] : 0.0
            F_s = (i < Ny_val) ? flux_y[i, j] : 0.0
            m_fe[i, j] += (F_w - F_e + F_n - F_s)

            if track_volatiles
                F_H_w = (j > 1) ? flux_H_x[i, j - 1] : 0.0
                F_H_e = (j < Nx_val) ? flux_H_x[i, j] : 0.0
                F_H_n = (i > 1) ? flux_H_y[i - 1, j] : 0.0
                F_H_s = (i < Ny_val) ? flux_H_y[i, j] : 0.0
                m_fe_H[i, j] = max(0.0, m_fe_H[i, j] + (F_H_w - F_H_e + F_H_n - F_H_s))

                F_C_w = (j > 1) ? flux_C_x[i, j - 1] : 0.0
                F_C_e = (j < Nx_val) ? flux_C_x[i, j] : 0.0
                F_C_n = (i > 1) ? flux_C_y[i - 1, j] : 0.0
                F_C_s = (i < Ny_val) ? flux_C_y[i, j] : 0.0
                m_fe_C[i, j] = max(0.0, m_fe_C[i, j] + (F_C_w - F_C_e + F_C_n - F_C_s))

                F_N_w = (j > 1) ? flux_N_x[i, j - 1] : 0.0
                F_N_e = (j < Nx_val) ? flux_N_x[i, j] : 0.0
                F_N_n = (i > 1) ? flux_N_y[i - 1, j] : 0.0
                F_N_s = (i < Ny_val) ? flux_N_y[i, j] : 0.0
                m_fe_N[i, j] = max(0.0, m_fe_N[i, j] + (F_N_w - F_N_e + F_N_n - F_N_s))

                F_S_w = (j > 1) ? flux_S_x[i, j - 1] : 0.0
                F_S_e = (j < Nx_val) ? flux_S_x[i, j] : 0.0
                F_S_n = (i > 1) ? flux_S_y[i - 1, j] : 0.0
                F_S_s = (i < Ny_val) ? flux_S_y[i, j] : 0.0
                m_fe_S[i, j] = max(0.0, m_fe_S[i, j] + (F_S_w - F_S_e + F_S_n - F_S_s))
            end
        end

        # 6. Gravitational potential energy dissipation heating
        @inbounds for j in 1:Nx_val, i in 1:Ny_val
            n_m = M_rock_markers[i, j]
            v_s = v_seg_cell[i, j]
            if n_m > 0 && v_s > 0.0
                phi_m_curr = (m_fe[i, j] / n_m) * Xfem_cell[i, j]
                Q_diss = segregation_dissipation_heating(
                    min(phi_m_curr, 1.0), drho_cell[i, j], g_acc_cell[i, j], v_s
                )
                total_diss_energy += Q_diss * (dx_val * dy_val) * dt_sub
                if Q_seg_grid !== nothing
                    dQ = 0.25 * Q_diss * (dt_sub / dt)
                    if i <= size(Q_seg_grid, 1) && j <= size(Q_seg_grid, 2)
                        Q_seg_grid[i, j] += dQ
                    end
                    if i <= size(Q_seg_grid, 1) && (j + 1) <= size(Q_seg_grid, 2)
                        Q_seg_grid[i, j + 1] += dQ
                    end
                    if (i + 1) <= size(Q_seg_grid, 1) && j <= size(Q_seg_grid, 2)
                        Q_seg_grid[i + 1, j] += dQ
                    end
                    if (i + 1) <= size(Q_seg_grid, 1) && (j + 1) <= size(Q_seg_grid, 2)
                        Q_seg_grid[i + 1, j + 1] += dQ
                    end
                end
            end
        end
    end

    # Distribute net cell mass changes to markers in each cell
    initial_sum = 0.0
    @inbounds for m in 1:marknum
        initial_sum += Xfe_bulk[m]
    end

    @inbounds for m in 1:marknum
        if tm[m] < 3
            rmark = distance(xm[m], ym[m], xcenter, ycenter)
            if rmark <= rplanet
                j_c = clamp(Int(floor(xm[m] / dx_val)) + 1, 1, Nx_val)
                i_c = clamp(Int(floor(ym[m] / dy_val)) + 1, 1, Ny_val)
                n_m = M_rock_markers[i_c, j_c]
                if n_m > 0
                    m_target = m_fe[i_c, j_c]
                    m_init = M_fe_cell[i_c, j_c]
                    dm_cell = m_target - m_init
                    dX = 0.0
                    if dm_cell > 0.0
                        c_tot = cap_cell[i_c, j_c]
                        if c_tot > 0.0
                            frac_gain = min(dm_cell / c_tot, 1.0)
                            dX = frac_gain * max(cfg_core.phi_pack - Xfe_bulk[m], 0.0)
                            Xfe_bulk[m] = clamp(Xfe_bulk[m] + dX, 0.0, cfg_core.phi_pack)
                        else
                            Xfe_bulk[m] = clamp(Xfe_bulk[m], 0.0, cfg_core.phi_pack)
                        end
                    elseif dm_cell < 0.0
                        if m_init > 0.0
                            scale_loss = max(m_target / m_init, 0.0)
                            Xfe_bulk[m] = clamp(
                                Xfe_bulk[m] * scale_loss, 0.0, cfg_core.phi_pack
                            )
                        else
                            Xfe_bulk[m] = 0.0
                        end
                    else
                        Xfe_bulk[m] = clamp(Xfe_bulk[m], 0.0, cfg_core.phi_pack)
                    end

                    if track_volatiles
                        X_new = Xfe_bulk[m]
                        if X_new > 0.0 && m_target > 0.0
                            Xfe_H_m[m] = m_fe_H[i_c, j_c] / m_target
                            Xfe_C_m[m] = m_fe_C[i_c, j_c] / m_target
                            Xfe_N_m[m] = m_fe_N[i_c, j_c] / m_target
                            Xfe_S_m[m] = m_fe_S[i_c, j_c] / m_target
                        else
                            Xfe_H_m[m] = 0.0
                            Xfe_C_m[m] = 0.0
                            Xfe_N_m[m] = 0.0
                            Xfe_S_m[m] = 0.0
                        end
                    end
                end
            end
        end
    end

    # Enforce floating point conservation without creating out-of-bounds markers
    final_sum = 0.0
    @inbounds for m in 1:marknum
        final_sum += Xfe_bulk[m]
    end

    diff_sum = initial_sum - final_sum
    if abs(diff_sum) > 1.0e-12 * initial_sum
        eligible_count = 0
        @inbounds for m in 1:marknum
            if tm[m] < 3 && distance(xm[m], ym[m], xcenter, ycenter) <= rplanet
                if diff_sum > 0.0 && Xfe_bulk[m] < cfg_core.phi_pack
                    eligible_count += 1
                elseif diff_sum < 0.0 && Xfe_bulk[m] > 0.0
                    eligible_count += 1
                end
            end
        end
        if eligible_count > 0
            corr = diff_sum / eligible_count
            @inbounds for m in 1:marknum
                if tm[m] < 3 && distance(xm[m], ym[m], xcenter, ycenter) <= rplanet
                    if diff_sum > 0.0 && Xfe_bulk[m] < cfg_core.phi_pack
                        Xfe_bulk[m] = min(Xfe_bulk[m] + corr, cfg_core.phi_pack)
                    elseif diff_sum < 0.0 && Xfe_bulk[m] > 0.0
                        Xfe_bulk[m] = max(Xfe_bulk[m] + corr, 0.0)
                    end
                end
            end
        end
    end

    # Keep molten metal fraction Xfem consistent with newly segregated bulk metal
    @inbounds for m in 1:marknum
        if tm[m] < 3
            F_fe = compute_metal_melt_fraction(
                tkm[m]; T_eutectic=cfg_core.T_eutectic, dT_metal=cfg_core.dT_metal
            )
            Xfem[m] = Xfe_bulk[m] * F_fe
        else
            Xfem[m] = 0.0
        end
    end

    return (;
        max_v_seg=max_v,
        n_subcycles=n_sub,
        dt_sub=dt_sub,
        total_dissipation_energy=total_diss_energy,
    )
end

"""
Apply silicate melt segregation via sub-cycled conservative drift-flux transport.

$(SIGNATURES)

Solves the conservative drift-flux transport equation for buoyant silicate melt
percolation through a solid silicate mantle matrix and Stokes crystal settling through
magma mush and suspension regimes. Subcycles explicit finite-volume transport using
a local CFL criterion and applies multi-dimensional donor-receiver flux limiters.
Guarantees mass conservation of total silicate melt to machine precision when input
marker melt fractions satisfy 0 <= Fm <= phi_pack and subsolidus freezing is inactive.

References:
- McKenzie (1984), J. Petrol., 25(3), 713-765.
- Bercovici & Ricard (2003), J. Geophys. Res., 108(B5), 2258.
- Katz (2008), J. Petrol., 49(12), 2099-2121.

# Arguments
- `xm::AbstractVector{Float64}`: Marker x-coordinates [m]
- `ym::AbstractVector{Float64}`: Marker y-coordinates [m]
- `tm::AbstractVector{<:Integer}`: Marker material phase type
- `tkm::AbstractVector{Float64}`: Marker temperature [K]
- `Fm::AbstractVector{Float64}`: Marker silicate melt volume fraction [-]
- `marknum::Integer`: Number of markers
- `dt::Real`: Timestep duration [s]
- `cfg_magma::MagmaTransportConfig`: Magma transport configuration parameters

# Keyword Arguments
- `coords=nothing`: `GridCoordinates` domain geometry struct
- `xcenter::Real=coords !== nothing ? coords.xcenter : 70000.0`: Planet center x [m]
- `ycenter::Real=coords !== nothing ? coords.ycenter : 70000.0`: Planet center y [m]
- `rplanet::Real=50000.0`: Planet radius [m]
- `g_surf::Real=0.1`: Reference surface gravity magnitude [m/s^2]
- `gx::Union{Nothing,AbstractMatrix{Float64}}=nothing`: Optional x-gravity on grid [m/s^2]
- `gy::Union{Nothing,AbstractMatrix{Float64}}=nothing`: Optional y-gravity on grid [m/s^2]
- `Q_seg_grid::Union{Nothing,AbstractMatrix{Float64}}=nothing`: Grid for dissipation heating [W/m^3]
- `Q_lat_grid::Union{Nothing,AbstractMatrix{Float64}}=nothing`: Grid for crystallization latent heat [W/m^3]
- `rho_silicate::Real=3000.0`: Reference solid rock density [kg/m^3]
- `rho_melt::Real=2800.0`: Reference silicate melt density [kg/m^3]
- `eta_silicate::Real=1.0e18`: Reference solid rock dynamic viscosity [Pa s]
- `ETA::Union{Nothing,AbstractMatrix{Float64}}=nothing`: Matrix viscosity on grid [Pa s]
- `T_solidus_silicate::Real=1400.0`: Reference silicate solidus temperature [K]
- `T_liquidus_silicate::Real=1800.0`: Reference silicate liquidus temperature [K]
- `L_melt::Real=4.0e5`: Latent heat of silicate melting [J/kg]
- `F_extract_m::Union{Nothing,AbstractVector{Float64}}=nothing`: Extracted melt fraction array

# Returns
- NamedTuple `(; max_v_seg, n_subcycles, dt_sub, total_dissipation_energy, total_crystallized_mass)`
"""
function exsolve_magma_volatiles!(
    xm::AbstractVector{Float64},
    ym::AbstractVector{Float64},
    tm::AbstractVector{<:Integer},
    tkm::AbstractVector{Float64},
    Fm::AbstractVector{Float64},
    marknum::Integer,
    XH2Om::AbstractVector{Float64},
    XCm::Union{Nothing,AbstractVector{Float64}},
    XNm::Union{Nothing,AbstractVector{Float64}},
    XSm::Union{Nothing,AbstractVector{Float64}},
    phim::AbstractVector{Float64},
    cfg_volatiles::VolatilesConfig;
    pr::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    xcenter::Real=0.0,
    ycenter::Real=0.0,
    rplanet::Real=1.0e6,
    dx_val::Real=1000.0,
    dy_val::Real=1000.0,
    Nx_val::Integer=32,
    Ny_val::Integer=32,
    rho_silicate::Real=3000.0,
    g_surf::Real=1.0,
)
    total_exsolved_vol = 0.0
    @inbounds for m in 1:marknum
        if tm[m] < 3 && distance(xm[m], ym[m], xcenter, ycenter) <= rplanet && Fm[m] > 0.0
            j_c = clamp(Int(floor(xm[m] / dx_val)) + 1, 1, Nx_val)
            i_c = clamp(Int(floor(ym[m] / dy_val)) + 1, 1, Ny_val)
            P_marker = if pr !== nothing
                i_pr = size(pr, 1) >= Ny_val + 1 ? i_c + 1 : min(i_c, size(pr, 1))
                j_pr = size(pr, 2) >= Nx_val + 1 ? j_c + 1 : min(j_c, size(pr, 2))
                max(pr[i_pr, j_pr], 0.0)
            else
                rmark = distance(xm[m], ym[m], xcenter, ycenter)
                rho_silicate * g_surf * max(rplanet - rmark, 0.0)
            end
            T_marker = tkm[m]
            dw = update_single_marker_volatile_exsolution!(
                m,
                Fm[m],
                P_marker,
                T_marker,
                XH2Om,
                XCm,
                XNm,
                XSm,
                phim,
                cfg_volatiles;
                rhosolid=rho_silicate,
                rhofluid=1000.0,
                phimax=0.9999,
            )
            total_exsolved_vol += dw
        end
    end
    return total_exsolved_vol
end

function apply_silicate_melt_segregation!(
    xm::AbstractVector{Float64},
    ym::AbstractVector{Float64},
    tm::AbstractVector{<:Integer},
    tkm::AbstractVector{Float64},
    Fm::AbstractVector{Float64},
    marknum::Integer,
    dt::Real,
    cfg_magma::MagmaTransportConfig;
    coords=nothing,
    xcenter::Real=coords !== nothing ? coords.xcenter : 70000.0,
    ycenter::Real=coords !== nothing ? coords.ycenter : 70000.0,
    rplanet::Real=50000.0,
    g_surf::Real=0.1,
    gx::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    gy::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    Q_seg_grid::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    Q_lat_grid::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    rho_silicate::Real=3000.0,
    rho_melt::Real=2800.0,
    eta_silicate::Real=1.0e18,
    ETA::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    T_solidus_silicate::Real=1400.0,
    T_liquidus_silicate::Real=1800.0,
    L_melt::Real=4.0e5,
    F_extract_m::Union{Nothing,AbstractVector{Float64}}=nothing,
    vx::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    vy::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    pr::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    div_v::Union{Nothing,AbstractMatrix{Float64}}=nothing,
    XH2Om::Union{Nothing,AbstractVector{Float64}}=nothing,
    XCm::Union{Nothing,AbstractVector{Float64}}=nothing,
    XNm::Union{Nothing,AbstractVector{Float64}}=nothing,
    XSm::Union{Nothing,AbstractVector{Float64}}=nothing,
    phim::Union{Nothing,AbstractVector{Float64}}=nothing,
    cfg_volatiles::Union{Nothing,VolatilesConfig}=nothing,
    workspace=nothing,
)
    if !cfg_magma.active || dt <= 0.0 || marknum <= 0
        return (;
            max_v_seg=0.0,
            n_subcycles=0,
            dt_sub=0.0,
            total_dissipation_energy=0.0,
            total_crystallized_mass=0.0,
            max_compaction_pressure=0.0,
            mean_compaction_length=0.0,
            total_exsolved_volatiles=0.0,
        )
    end

    # Validate input marker bounds
    @inbounds for m in 1:marknum
        if tm[m] < 3
            rmark = distance(xm[m], ym[m], xcenter, ycenter)
            if rmark <= rplanet
                fm = Fm[m]
                if !isfinite(fm) || fm < 0.0 || fm > 1.0 + 1.0e-7
                    throw(
                        DomainError(
                            fm, "Marker melt fraction must be finite and within [0, 1]"
                        ),
                    )
                end
            end
        end
    end

    Nx_val = coords !== nothing ? coords.Nx : 32
    Ny_val = coords !== nothing ? coords.Ny : 32
    dx_val = if coords !== nothing
        coords.dx
    else
        (coords !== nothing ? coords.xsize / Nx_val : 4375.0)
    end
    dy_val = if coords !== nothing
        coords.dy
    else
        (coords !== nothing ? coords.ysize / Ny_val : 4375.0)
    end

    has_ws = (workspace !== nothing && size(workspace.M_melt_cell) == (Ny_val, Nx_val))

    # Allocate or reset cell accumulations
    M_melt_cell =
        has_ws ? fill!(workspace.M_melt_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)
    M_rock_markers =
        has_ws ? fill!(workspace.M_rock_markers, 0) : zeros(Int, Ny_val, Nx_val)
    v_seg_cell = has_ws ? fill!(workspace.v_seg_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)
    F_m_cell = has_ws ? fill!(workspace.Fm_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)
    g_acc_cell = has_ws ? fill!(workspace.g_acc_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)
    cap_cell = has_ws ? fill!(workspace.cap_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)
    T_cell = has_ws ? fill!(workspace.T_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)
    drho_cell = has_ws ? fill!(workspace.drho_cell, 0.0) : zeros(Float64, Ny_val, Nx_val)
    P_comp_cell = has_ws ? fill!(workspace.P_comp, 0.0) : zeros(Float64, Ny_val, Nx_val)
    div_v_cell = has_ws ? fill!(workspace.div_v, 0.0) : zeros(Float64, Ny_val, Nx_val)

    # Bin markers into grid cells
    @inbounds for m in 1:marknum
        if tm[m] < 3
            rmark = distance(xm[m], ym[m], xcenter, ycenter)
            if rmark <= rplanet
                j_c = clamp(Int(floor(xm[m] / dx_val)) + 1, 1, Nx_val)
                i_c = clamp(Int(floor(ym[m] / dy_val)) + 1, 1, Ny_val)
                fm_val = Fm[m]
                M_melt_cell[i_c, j_c] += fm_val
                M_rock_markers[i_c, j_c] += 1
                T_cell[i_c, j_c] += tkm[m]
                cap_cell[i_c, j_c] += max(cfg_magma.phi_pack - fm_val, 0.0)
            end
        end
    end

    # Normalize cell averages
    @inbounds for j in 1:Nx_val, i in 1:Ny_val
        n_m = M_rock_markers[i, j]
        if n_m > 0
            F_m_cell[i, j] = M_melt_cell[i, j] / n_m
            T_cell[i, j] /= n_m
        end
    end

    # Compute cell segregation velocities
    @inbounds for j in 1:Nx_val, i in 1:Ny_val
        xc = (j - 0.5) * dx_val
        yc = (i - 0.5) * dy_val
        r_c = distance(xc, yc, xcenter, ycenter)
        g_acc = if gx !== nothing && gy !== nothing
            hypot(gx[i, j], gy[i, j])
        else
            g_surf * (r_c / rplanet)
        end
        g_acc_cell[i, j] = g_acc

        F_m = F_m_cell[i, j]
        drho = rho_silicate - rho_melt
        drho_cell[i, j] = max(drho, 0.0)

        if F_m > cfg_magma.phi_residual && g_acc > 0.0 && drho > 0.0
            v_seg_cell[i, j] = silicate_melt_segregation_velocity(
                F_m,
                drho,
                g_acc,
                cfg_magma.eta_melt;
                k_melt_ref=cfg_magma.k_melt_ref,
                phi0=cfg_magma.phi0,
                perm_exponent=cfg_magma.perm_exponent,
                phi_residual=cfg_magma.phi_residual,
                phi_crit=cfg_magma.phi_crit,
                r_grain=cfg_magma.r_grain,
                hindered_exponent=cfg_magma.hindered_exponent,
                F_perc_end=cfg_magma.F_perc_end,
                F_settle_start=cfg_magma.F_settle_start,
            )
        end
    end

    # Evaluate solid matrix divergence and dynamic compaction pressure
    if cfg_magma.compaction_active
        if div_v !== nothing
            @inbounds for j in 1:Nx_val, i in 1:Ny_val
                div_v_cell[i, j] = div_v[min(i, size(div_v, 1)), min(j, size(div_v, 2))]
            end
        elseif vx !== nothing && vy !== nothing
            @inbounds for j in 1:Nx_val, i in 1:Ny_val
                # On Gerya staggered grid (size Ny1 x Nx1):
                # Cell (i, j) center is at xp[j+1], yp[i+1].
                # Horizontal velocity Vx is along row i+1, across faces j (left) and j+1 (right).
                # Vertical velocity Vy is along col j+1, across faces i (top) and i+1 (bottom).
                dvx = if size(vx, 2) >= Nx_val + 1 && size(vx, 1) >= Ny_val + 1
                    (vx[i + 1, j + 1] - vx[i + 1, j]) / dx_val
                elseif size(vx, 2) >= Nx_val + 1
                    (vx[min(i, size(vx, 1)), j + 1] - vx[min(i, size(vx, 1)), j]) / dx_val
                else
                    j_l = max(j - 1, 1)
                    j_r = min(j + 1, size(vx, 2))
                    (vx[min(i, size(vx, 1)), j_r] - vx[min(i, size(vx, 1)), j_l]) /
                    (max(j_r - j_l, 1) * dx_val)
                end
                dvy = if size(vy, 1) >= Ny_val + 1 && size(vy, 2) >= Nx_val + 1
                    (vy[i + 1, j + 1] - vy[i, j + 1]) / dy_val
                elseif size(vy, 1) >= Ny_val + 1
                    (vy[i + 1, min(j, size(vy, 2))] - vy[i, min(j, size(vy, 2))]) / dy_val
                else
                    i_b = max(i - 1, 1)
                    i_t = min(i + 1, size(vy, 1))
                    (vy[i_t, min(j, size(vy, 2))] - vy[i_b, min(j, size(vy, 2))]) /
                    (max(i_t - i_b, 1) * dy_val)
                end
                div_v_cell[i, j] = dvx + dvy
            end
        end

        @inbounds for j in 1:Nx_val, i in 1:Ny_val
            F_m = F_m_cell[i, j]
            if F_m > cfg_magma.phi_residual
                eta_s = if ETA !== nothing
                    i_eta = size(ETA, 1) >= Ny_val + 1 ? i + 1 : min(i, size(ETA, 1))
                    j_eta = size(ETA, 2) >= Nx_val + 1 ? j + 1 : min(j, size(ETA, 2))
                    ETA[i_eta, j_eta]
                else
                    eta_silicate
                end
                P_comp_cell[i, j] = compaction_pressure(
                    div_v_cell[i, j],
                    max(eta_s, 1.0e-3),
                    F_m;
                    bulk_ratio=cfg_magma.bulk_viscosity_ratio,
                    phi_min=cfg_magma.min_bulk_porosity,
                )
            else
                P_comp_cell[i, j] = 0.0
            end
        end
    end

    max_P_comp = cfg_magma.compaction_active ? maximum(P_comp_cell) : 0.0
    sum_delta_c = 0.0
    count_delta_c = 0
    if cfg_magma.compaction_active
        @inbounds for j in 1:Nx_val, i in 1:Ny_val
            F_m = F_m_cell[i, j]
            if F_m > cfg_magma.phi_residual
                km = silicate_melt_permeability(
                    F_m;
                    k0=cfg_magma.k_melt_ref,
                    phi0=cfg_magma.phi0,
                    n=cfg_magma.perm_exponent,
                    phi_residual=cfg_magma.phi_residual,
                    phi_crit=cfg_magma.phi_crit,
                )
                eta_s = if ETA !== nothing
                    ETA[min(i, size(ETA, 1)), min(j, size(ETA, 2))]
                else
                    eta_silicate
                end
                delta_c = compaction_length(
                    max(eta_s, 1.0e-3),
                    cfg_magma.eta_melt,
                    km,
                    F_m;
                    bulk_ratio=cfg_magma.bulk_viscosity_ratio,
                    phi_min=cfg_magma.min_bulk_porosity,
                    delta_min=cfg_magma.compaction_length_min,
                    delta_max=cfg_magma.compaction_length_max,
                )
                sum_delta_c += delta_c
                count_delta_c += 1
            end
        end
    end

    max_v = maximum(v_seg_cell)
    if cfg_magma.compaction_active
        max_v_comp = 0.0
        @inbounds for j in 1:Nx_val, i in 1:Ny_val
            F_m = F_m_cell[i, j]
            if F_m > cfg_magma.phi_residual
                km = silicate_melt_permeability(
                    F_m;
                    k0=cfg_magma.k_melt_ref,
                    phi0=cfg_magma.phi0,
                    n=cfg_magma.perm_exponent,
                    phi_residual=cfg_magma.phi_residual,
                    phi_crit=cfg_magma.phi_crit,
                )
                dPx = 0.0
                if j < Nx_val
                    dPx = max(dPx, abs(P_comp_cell[i, j] - P_comp_cell[i, j + 1]) / dx_val)
                end
                if j > 1
                    dPx = max(dPx, abs(P_comp_cell[i, j] - P_comp_cell[i, j - 1]) / dx_val)
                end
                dPy = 0.0
                if i < Ny_val
                    dPy = max(dPy, abs(P_comp_cell[i, j] - P_comp_cell[i + 1, j]) / dy_val)
                end
                if i > 1
                    dPy = max(dPy, abs(P_comp_cell[i, j] - P_comp_cell[i - 1, j]) / dy_val)
                end
                v_c = (km / (cfg_magma.eta_melt * F_m)) * hypot(dPx, dPy)
                if v_c > max_v_comp
                    max_v_comp = v_c
                end
            end
        end
        max_v = max(max_v, max_v_comp)
    end

    total_exsolved_vol = 0.0
    if cfg_magma.exsolution_active &&
        cfg_volatiles !== nothing &&
        cfg_volatiles.active &&
        XH2Om !== nothing &&
        phim !== nothing
        total_exsolved_vol = exsolve_magma_volatiles!(
            xm,
            ym,
            tm,
            tkm,
            Fm,
            marknum,
            XH2Om,
            XCm,
            XNm,
            XSm,
            phim,
            cfg_volatiles;
            pr=pr,
            xcenter=xcenter,
            ycenter=ycenter,
            rplanet=rplanet,
            dx_val=dx_val,
            dy_val=dy_val,
            Nx_val=Nx_val,
            Ny_val=Ny_val,
            rho_silicate=rho_silicate,
            g_surf=g_surf,
        )
    end

    if max_v <= 0.0
        return (;
            max_v_seg=0.0,
            n_subcycles=0,
            dt_sub=0.0,
            total_dissipation_energy=0.0,
            total_crystallized_mass=0.0,
            max_compaction_pressure=max_P_comp,
            mean_compaction_length=count_delta_c > 0 ? sum_delta_c / count_delta_c : 0.0,
            total_exsolved_volatiles=total_exsolved_vol,
        )
    end

    # CFL calculation and subcycling
    dt_cfl = cfg_magma.cfl_melt * min(dx_val, dy_val) / max_v
    n_sub_raw = Int(ceil(dt / dt_cfl))
    if n_sub_raw > cfg_magma.max_subcycles
        @warn "CFL subcycling requires $n_sub_raw steps, capped at max_subcycles $(cfg_magma.max_subcycles); magma segregation throttled" maxlog=10
    end
    n_sub = clamp(n_sub_raw, 1, cfg_magma.max_subcycles)
    dt_sub = dt / n_sub

    # Working copy of cell melt mass for subcycling
    m_melt = has_ws ? copyto!(workspace.m_melt, M_melt_cell) : copy(M_melt_cell)
    total_diss_energy = 0.0
    total_cryst_mass = 0.0

    req_flux_x =
        has_ws ? fill!(workspace.req_flux_x, 0.0) : zeros(Float64, Ny_val, Nx_val - 1)
    req_flux_y =
        has_ws ? fill!(workspace.req_flux_y, 0.0) : zeros(Float64, Ny_val - 1, Nx_val)
    flux_x = has_ws ? fill!(workspace.flux_x, 0.0) : zeros(Float64, Ny_val, Nx_val - 1)
    flux_y = has_ws ? fill!(workspace.flux_y, 0.0) : zeros(Float64, Ny_val - 1, Nx_val)
    outflow_tot =
        has_ws ? fill!(workspace.outflow_tot, 0.0) : zeros(Float64, Ny_val, Nx_val)
    inflow_tot = has_ws ? fill!(workspace.inflow_tot, 0.0) : zeros(Float64, Ny_val, Nx_val)
    alpha_out = has_ws ? fill!(workspace.alpha_out, 1.0) : ones(Float64, Ny_val, Nx_val)
    alpha_in = has_ws ? fill!(workspace.alpha_in, 1.0) : ones(Float64, Ny_val, Nx_val)

    # Subcycling loop
    for _ in 1:n_sub
        fill!(outflow_tot, 0.0)
        fill!(inflow_tot, 0.0)
        fill!(req_flux_x, 0.0)
        fill!(req_flux_y, 0.0)

        # 1. Compute unscaled requested fluxes across East-West faces
        @inbounds for j in 1:(Nx_val - 1)
            xf = j * dx_val
            for i in 1:Ny_val
                yf = (i - 0.5) * dy_val
                dxf = xf - xcenter
                dyf = yf - ycenter
                rf = sqrt(dxf^2 + dyf^2)
                if rf > rplanet || rf < 1.0e-3
                    continue
                end

                # Outward buoyancy unit normal: opposite to gravity
                nx =
                    if gx !== nothing &&
                        gy !== nothing &&
                        i <= size(gx, 1) &&
                        j <= size(gx, 2)
                        gx_f = gx[i, j]
                        gy_f =
                            0.5 *
                            (gy[i, j] + (j + 1 <= size(gy, 2) ? gy[i, j + 1] : gy[i, j]))
                        g_f = sqrt(gx_f^2 + gy_f^2)
                        g_f > 1.0e-10 ? -gx_f / g_f : dxf / rf
                    else
                        dxf / rf
                    end

                vf = 0.5 * (v_seg_cell[i, j] + v_seg_cell[i, j + 1])
                uf = vf * nx
                if cfg_magma.compaction_active
                    dP_x = (P_comp_cell[i, j] - P_comp_cell[i, j + 1]) / dx_val
                    F_face_x = 0.5 * (F_m_cell[i, j] + F_m_cell[i, j + 1])
                    if F_face_x > cfg_magma.phi_residual
                        km_face_x = silicate_melt_permeability(
                            F_face_x;
                            k0=cfg_magma.k_melt_ref,
                            phi0=cfg_magma.phi0,
                            n=cfg_magma.perm_exponent,
                            phi_residual=cfg_magma.phi_residual,
                            phi_crit=cfg_magma.phi_crit,
                        )
                        uf += (km_face_x / (cfg_magma.eta_melt * F_face_x)) * dP_x
                    end
                end
                if cfg_magma.ponding_active
                    if (uf > 0.0 && T_cell[i, j + 1] < T_solidus_silicate) ||
                        (uf < 0.0 && T_cell[i, j] < T_solidus_silicate)
                        P_donor = uf > 0.0 ? P_comp_cell[i, j] : P_comp_cell[i, j + 1]
                        if !cfg_magma.eruption_active ||
                            P_donor <= cfg_magma.tensile_strength
                            uf = 0.0
                        end
                    end
                end

                donor_j = uf > 0.0 ? j : j + 1
                rec_j = uf > 0.0 ? j + 1 : j

                n_donor = M_rock_markers[i, donor_j]
                n_rec = M_rock_markers[i, rec_j]
                if n_donor == 0 || n_rec == 0
                    continue
                end

                F_donor = m_melt[i, donor_j] / n_donor
                F_mob = max(F_donor - cfg_magma.phi_residual, 0.0)
                m_avail = F_mob * n_donor

                fx = abs(uf) * (dt_sub / dx_val) * m_avail
                fx_req = uf > 0.0 ? fx : -fx
                req_flux_x[i, j] = fx_req

                if fx_req > 0.0
                    outflow_tot[i, j] += fx_req
                    inflow_tot[i, j + 1] += fx_req
                else
                    outflow_tot[i, j + 1] += -fx_req
                    inflow_tot[i, j] += -fx_req
                end
            end
        end

        # 2. Compute unscaled requested fluxes across North-South faces
        @inbounds for i in 1:(Ny_val - 1)
            yf = i * dy_val
            for j in 1:Nx_val
                xf = (j - 0.5) * dx_val
                dxf = xf - xcenter
                dyf = yf - ycenter
                rf = sqrt(dxf^2 + dyf^2)
                if rf > rplanet || rf < 1.0e-3
                    continue
                end

                ny =
                    if gx !== nothing &&
                        gy !== nothing &&
                        i <= size(gy, 1) &&
                        j <= size(gy, 2)
                        gy_f = gy[i, j]
                        gx_f =
                            0.5 *
                            (gx[i, j] + (i + 1 <= size(gx, 1) ? gx[i + 1, j] : gx[i, j]))
                        g_f = sqrt(gx_f^2 + gy_f^2)
                        g_f > 1.0e-10 ? -gy_f / g_f : dyf / rf
                    else
                        dyf / rf
                    end

                vf = 0.5 * (v_seg_cell[i, j] + v_seg_cell[i + 1, j])
                wf = vf * ny
                if cfg_magma.compaction_active
                    dP_y = (P_comp_cell[i, j] - P_comp_cell[i + 1, j]) / dy_val
                    F_face_y = 0.5 * (F_m_cell[i, j] + F_m_cell[i + 1, j])
                    if F_face_y > cfg_magma.phi_residual
                        km_face_y = silicate_melt_permeability(
                            F_face_y;
                            k0=cfg_magma.k_melt_ref,
                            phi0=cfg_magma.phi0,
                            n=cfg_magma.perm_exponent,
                            phi_residual=cfg_magma.phi_residual,
                            phi_crit=cfg_magma.phi_crit,
                        )
                        wf += (km_face_y / (cfg_magma.eta_melt * F_face_y)) * dP_y
                    end
                end
                if cfg_magma.ponding_active
                    if (wf > 0.0 && T_cell[i + 1, j] < T_solidus_silicate) ||
                        (wf < 0.0 && T_cell[i, j] < T_solidus_silicate)
                        P_donor = wf > 0.0 ? P_comp_cell[i, j] : P_comp_cell[i + 1, j]
                        if !cfg_magma.eruption_active ||
                            P_donor <= cfg_magma.tensile_strength
                            wf = 0.0
                        end
                    end
                end

                donor_i = wf > 0.0 ? i : i + 1
                rec_i = wf > 0.0 ? i + 1 : i

                n_donor = M_rock_markers[donor_i, j]
                n_rec = M_rock_markers[rec_i, j]
                if n_donor == 0 || n_rec == 0
                    continue
                end

                F_donor = m_melt[donor_i, j] / n_donor
                F_mob = max(F_donor - cfg_magma.phi_residual, 0.0)
                m_avail = F_mob * n_donor

                fy = abs(wf) * (dt_sub / dy_val) * m_avail
                fy_req = wf > 0.0 ? fy : -fy
                req_flux_y[i, j] = fy_req

                if fy_req > 0.0
                    outflow_tot[i, j] += fy_req
                    inflow_tot[i + 1, j] += fy_req
                else
                    outflow_tot[i + 1, j] += -fy_req
                    inflow_tot[i, j] += -fy_req
                end
            end
        end

        # 3. Multi-dimensional flux limiters per cell
        @inbounds for j in 1:Nx_val, i in 1:Ny_val
            n_m = M_rock_markers[i, j]
            if n_m > 0
                F_c = m_melt[i, j] / n_m
                m_avail = max(F_c - cfg_magma.phi_residual, 0.0) * n_m
                m_cap = max(cfg_magma.phi_pack - F_c, 0.0) * n_m
                alpha_out[i, j] = if outflow_tot[i, j] > m_avail && m_avail > 0.0
                    m_avail / outflow_tot[i, j]
                else
                    (outflow_tot[i, j] > m_avail ? 0.0 : 1.0)
                end
                alpha_in[i, j] = if inflow_tot[i, j] > m_cap && m_cap > 0.0
                    m_cap / inflow_tot[i, j]
                else
                    (inflow_tot[i, j] > m_cap ? 0.0 : 1.0)
                end
            else
                alpha_out[i, j] = 0.0
                alpha_in[i, j] = 0.0
            end
        end

        # 4. Scale fluxes by joint donor-receiver limiters
        @inbounds for j in 1:(Nx_val - 1), i in 1:Ny_val
            fx_req = req_flux_x[i, j]
            if iszero(fx_req)
                flux_x[i, j] = 0.0
            else
                donor_j = fx_req > 0.0 ? j : j + 1
                rec_j = fx_req > 0.0 ? j + 1 : j
                lim = min(alpha_out[i, donor_j], alpha_in[i, rec_j])
                flux_x[i, j] = fx_req * lim
            end
        end

        @inbounds for j in 1:Nx_val, i in 1:(Ny_val - 1)
            fy_req = req_flux_y[i, j]
            if iszero(fy_req)
                flux_y[i, j] = 0.0
            else
                donor_i = fy_req > 0.0 ? i : i + 1
                rec_i = fy_req > 0.0 ? i + 1 : i
                lim = min(alpha_out[donor_i, j], alpha_in[rec_i, j])
                flux_y[i, j] = fy_req * lim
            end
        end

        # 5. Conservative update of cell melt masses
        @inbounds for j in 1:Nx_val, i in 1:Ny_val
            F_w = (j > 1) ? flux_x[i, j - 1] : 0.0
            F_e = (j < Nx_val) ? flux_x[i, j] : 0.0
            F_n = (i > 1) ? flux_y[i - 1, j] : 0.0
            F_s = (i < Ny_val) ? flux_y[i, j] : 0.0
            m_melt[i, j] += (F_w - F_e + F_n - F_s)
        end

        # 6. Gravitational potential energy dissipation heating
        @inbounds for j in 1:Nx_val, i in 1:Ny_val
            n_m = M_rock_markers[i, j]
            v_s = v_seg_cell[i, j]
            F_w = (j > 1) ? flux_x[i, j - 1] : 0.0
            F_e = (j < Nx_val) ? flux_x[i, j] : 0.0
            F_n = (i > 1) ? flux_y[i - 1, j] : 0.0
            F_s = (i < Ny_val) ? flux_y[i, j] : 0.0
            flux_thru = 0.5 * (abs(F_w) + abs(F_e) + abs(F_n) + abs(F_s))
            if n_m > 0 && v_s > 0.0 && flux_thru > 0.0
                F_curr = m_melt[i, j] / n_m
                Q_diss = silicate_melt_dissipation_heating(
                    min(F_curr, 1.0), drho_cell[i, j], g_acc_cell[i, j], v_s
                )
                total_diss_energy += Q_diss * (dx_val * dy_val) * dt_sub
                if Q_seg_grid !== nothing
                    dQ = 0.25 * Q_diss * (dt_sub / dt)
                    if i <= size(Q_seg_grid, 1) && j <= size(Q_seg_grid, 2)
                        Q_seg_grid[i, j] += dQ
                    end
                    if i <= size(Q_seg_grid, 1) && (j + 1) <= size(Q_seg_grid, 2)
                        Q_seg_grid[i, j + 1] += dQ
                    end
                    if (i + 1) <= size(Q_seg_grid, 1) && j <= size(Q_seg_grid, 2)
                        Q_seg_grid[i + 1, j] += dQ
                    end
                    if (i + 1) <= size(Q_seg_grid, 1) && (j + 1) <= size(Q_seg_grid, 2)
                        Q_seg_grid[i + 1, j + 1] += dQ
                    end
                end
            end
        end

        # 7. Subsolidus/subliquidus crystallization and latent heat release
        if cfg_magma.latent_crystallization
            @inbounds for j in 1:Nx_val, i in 1:Ny_val
                n_m = M_rock_markers[i, j]
                if n_m > 0
                    t_cell_val = T_cell[i, j]
                    if t_cell_val < T_liquidus_silicate
                        F_eq_cell = if t_cell_val <= T_solidus_silicate
                            0.0
                        else
                            (t_cell_val - T_solidus_silicate) /
                            max(T_liquidus_silicate - T_solidus_silicate, 1.0)
                        end
                        m_eq = n_m * F_eq_cell
                        dm_net = m_melt[i, j] - M_melt_cell[i, j]
                        if dm_net > 0.0 && m_melt[i, j] > m_eq
                            # Newly arrived melt exceeding thermodynamic equilibrium crystallizes
                            m_freeze = min(dm_net, m_melt[i, j] - m_eq)
                            m_melt[i, j] -= m_freeze
                            total_cryst_mass += m_freeze
                            if Q_lat_grid !== nothing
                                F_freeze = m_freeze / n_m
                                Q_cryst = (rho_melt * F_freeze * L_melt) / dt_sub
                                dQ_lat = 0.25 * Q_cryst * (dt_sub / dt)
                                if i <= size(Q_lat_grid, 1) && j <= size(Q_lat_grid, 2)
                                    Q_lat_grid[i, j] += dQ_lat
                                end
                                if i <= size(Q_lat_grid, 1) &&
                                    (j + 1) <= size(Q_lat_grid, 2)
                                    Q_lat_grid[i, j + 1] += dQ_lat
                                end
                                if (i + 1) <= size(Q_lat_grid, 1) &&
                                    j <= size(Q_lat_grid, 2)
                                    Q_lat_grid[i + 1, j] += dQ_lat
                                end
                                if (i + 1) <= size(Q_lat_grid, 1) &&
                                    (j + 1) <= size(Q_lat_grid, 2)
                                    Q_lat_grid[i + 1, j + 1] += dQ_lat
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # Distribute net cell mass changes to markers in each cell
    initial_sum = 0.0
    @inbounds for m in 1:marknum
        initial_sum += Fm[m]
    end

    @inbounds for m in 1:marknum
        if tm[m] < 3
            rmark = distance(xm[m], ym[m], xcenter, ycenter)
            if rmark <= rplanet
                j_c = clamp(Int(floor(xm[m] / dx_val)) + 1, 1, Nx_val)
                i_c = clamp(Int(floor(ym[m] / dy_val)) + 1, 1, Ny_val)
                n_m = M_rock_markers[i_c, j_c]
                if n_m > 0
                    m_target = m_melt[i_c, j_c]
                    m_init = M_melt_cell[i_c, j_c]
                    dm_cell = m_target - m_init
                    if dm_cell > 0.0
                        c_tot = cap_cell[i_c, j_c]
                        if c_tot > 0.0
                            frac_gain = min(dm_cell / c_tot, 1.0)
                            dX = frac_gain * max(cfg_magma.phi_pack - Fm[m], 0.0)
                            Fm[m] = clamp(Fm[m] + dX, 0.0, 1.0)
                        else
                            Fm[m] = clamp(Fm[m], 0.0, 1.0)
                        end
                    elseif dm_cell < 0.0
                        if m_init > 0.0
                            scale_loss = max(m_target / m_init, 0.0)
                            dF_lost = Fm[m] * (1.0 - scale_loss)
                            Fm[m] = clamp(Fm[m] * scale_loss, 0.0, 1.0)
                            if F_extract_m !== nothing
                                F_extract_m[m] = clamp(F_extract_m[m] + dF_lost, 0.0, 1.0)
                            end
                        else
                            if F_extract_m !== nothing
                                F_extract_m[m] = clamp(F_extract_m[m] + Fm[m], 0.0, 1.0)
                            end
                            Fm[m] = 0.0
                        end
                    else
                        Fm[m] = clamp(Fm[m], 0.0, 1.0)
                    end
                end
            end
        end
    end

    # Enforce floating point conservation without creating out-of-bounds markers
    if iszero(total_cryst_mass)
        final_sum = 0.0
        @inbounds for m in 1:marknum
            final_sum += Fm[m]
        end

        diff_sum = initial_sum - final_sum
        if abs(diff_sum) > 1.0e-12 * initial_sum && initial_sum > 0.0
            eligible_count = 0
            @inbounds for m in 1:marknum
                if tm[m] < 3 && distance(xm[m], ym[m], xcenter, ycenter) <= rplanet
                    if diff_sum > 0.0 && Fm[m] < 1.0
                        eligible_count += 1
                    elseif diff_sum < 0.0 && Fm[m] > 0.0
                        eligible_count += 1
                    end
                end
            end
            if eligible_count > 0
                corr = diff_sum / eligible_count
                @inbounds for m in 1:marknum
                    if tm[m] < 3 && distance(xm[m], ym[m], xcenter, ycenter) <= rplanet
                        if diff_sum > 0.0 && Fm[m] < 1.0
                            Fm[m] = min(Fm[m] + corr, 1.0)
                        elseif diff_sum < 0.0 && Fm[m] > 0.0
                            Fm[m] = max(Fm[m] + corr, 0.0)
                        end
                    end
                end
            end
        end
    end

    mean_delta_c = count_delta_c > 0 ? sum_delta_c / count_delta_c : 0.0
    return (;
        max_v_seg=max_v,
        n_subcycles=n_sub,
        dt_sub=dt_sub,
        total_dissipation_energy=total_diss_energy,
        total_crystallized_mass=total_cryst_mass,
        max_compaction_pressure=max_P_comp,
        mean_compaction_length=mean_delta_c,
        total_exsolved_volatiles=total_exsolved_vol,
    )
end

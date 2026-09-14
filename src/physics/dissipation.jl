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
    @unpack_coords coords dx dy
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

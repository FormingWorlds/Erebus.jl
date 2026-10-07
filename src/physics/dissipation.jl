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
    PHIX=nothing,
    PHIY=nothing,
    kphim0=nothing,
    phim0_val=nothing,
    phimin_val::Real=1.0e-4,
    kappa_frac::Real=1.0e3,
    gamma_frac::Real=1.0,
    k_frac_max::Real=1.0e-9,
    coords::GridCoordinates=default_grid_coordinates(),
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
            sigma_tx1 = 0.5 * (TEN[i, j - 1] + TEN[i, j])
            phi_x1 = if PHIX !== nothing
                PHIX[i, j - 1]
            elseif PHI !== nothing
                0.5 * (PHI[i, j - 1] + PHI[i, j])
            else
                phimin_val
            end
            rx_jm1 = evaluate_hydrofracture_resistance(
                RX[i, j - 1],
                Peff_x1,
                sigma_tx1,
                phi_x1;
                kphim0=kphim0,
                phim0_val=phim0_val,
                phimin_val=phimin_val,
                kappa_frac=kappa_frac,
                gamma_frac=gamma_frac,
                k_frac_max=k_frac_max,
            )

            Peff_x2 = 0.5 * (pr[i, j] + pr[i, j + 1] - pf[i, j] - pf[i, j + 1])
            sigma_tx2 = 0.5 * (TEN[i, j] + TEN[i, j + 1])
            phi_x2 = if PHIX !== nothing
                PHIX[i, j]
            elseif PHI !== nothing
                0.5 * (PHI[i, j] + PHI[i, j + 1])
            else
                phimin_val
            end
            rx_j = evaluate_hydrofracture_resistance(
                RX[i, j],
                Peff_x2,
                sigma_tx2,
                phi_x2;
                kphim0=kphim0,
                phim0_val=phim0_val,
                phimin_val=phimin_val,
                kappa_frac=kappa_frac,
                gamma_frac=gamma_frac,
                k_frac_max=k_frac_max,
            )

            Peff_y1 = 0.5 * (pr[i - 1, j] + pr[i, j] - pf[i - 1, j] - pf[i, j])
            sigma_ty1 = 0.5 * (TEN[i - 1, j] + TEN[i, j])
            phi_y1 = if PHIY !== nothing
                PHIY[i - 1, j]
            elseif PHI !== nothing
                0.5 * (PHI[i - 1, j] + PHI[i, j])
            else
                phimin_val
            end
            ry_im1 = evaluate_hydrofracture_resistance(
                RY[i - 1, j],
                Peff_y1,
                sigma_ty1,
                phi_y1;
                kphim0=kphim0,
                phim0_val=phim0_val,
                phimin_val=phimin_val,
                kappa_frac=kappa_frac,
                gamma_frac=gamma_frac,
                k_frac_max=k_frac_max,
            )

            Peff_y2 = 0.5 * (pr[i, j] + pr[i + 1, j] - pf[i, j] - pf[i + 1, j])
            sigma_ty2 = 0.5 * (TEN[i, j] + TEN[i + 1, j])
            phi_y2 = if PHIY !== nothing
                PHIY[i, j]
            elseif PHI !== nothing
                0.5 * (PHI[i, j] + PHI[i + 1, j])
            else
                phimin_val
            end
            ry_i = evaluate_hydrofracture_resistance(
                RY[i, j],
                Peff_y2,
                sigma_ty2,
                phi_y2;
                kphim0=kphim0,
                phim0_val=phim0_val,
                phimin_val=phimin_val,
                kappa_frac=kappa_frac,
                gamma_frac=gamma_frac,
                k_frac_max=k_frac_max,
            )
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
    HA,
    tk1,
    ALPHA,
    ALPHAF,
    PHI,
    vx,
    vy,
    vxf,
    vyf,
    ps,
    pf;
    coords::GridCoordinates=default_grid_coordinates(),
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

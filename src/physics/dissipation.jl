"""
Compute shear heating based on basic (temperature) and P grids.

$(SIGNATURES)

# Details

    - HS: shear heating
    - ETA: viscoplastic viscosity at basic nodes
    - SXY: σ₀xy XY stress at basic nodes
    - ETAP: viscosity at P nodes
    - SXX: normal stress at P nodes
    - RX: face resistance ηfluid/Kϕ at Vx nodes that `qxD` was solved with [Pa s m⁻²]
    - RY: face resistance ηfluid/Kϕ at Vy nodes that `qyD` was solved with [Pa s m⁻²]
    - qxD: qx-Darcy flux at Vx nodes
    - qyD: qy-Darcy flux at Vy nodes
    - PHI: porosity at P nodes
    - ETAPHI: bulk viscosity at P nodes
    - pr: total pressure at P nodes
    - pf: fluid pressure at P nodes

All P-node arrays and `RX`, `RY`, `qxD`, `qyD` must have the size of `HS`, and `ETA`,
`SXY` one row and column less; otherwise an `ArgumentError` is thrown.

# Returns

    - nothing
"""
function compute_shear_heating!(
    HS, ETA, SXY, ETAP, SXX, RX, RY, qxD, qyD, PHI, ETAPHI, pr, pf
)
    Ny1, Nx1 = size(HS)
    for (name, A) in (
        ("ETAP", ETAP),
        ("SXX", SXX),
        ("RX", RX),
        ("RY", RY),
        ("qxD", qxD),
        ("qyD", qyD),
        ("PHI", PHI),
        ("ETAPHI", ETAPHI),
        ("pr", pr),
        ("pf", pf),
    )
        size(A) == (Ny1, Nx1) ||
            throw(ArgumentError("$name has size $(size(A)), HS has size $((Ny1, Nx1))"))
    end
    for (name, A) in (("ETA", ETA), ("SXY", SXY))
        size(A) == (Ny1 - 1, Nx1 - 1) || throw(
            ArgumentError("$name has size $(size(A)), expected $((Ny1 - 1, Nx1 - 1))")
        )
    end
    Nx = Nx1 - 1
    Ny = Ny1 - 1
    for j in 2:1:Nx, i in 2:1:Ny
        # average SXY⋅EXY
        SXYEXY = 0.25 * sum(grid_vector(i-1, j-1, SXY) .^ 2 ./ grid_vector(i-1, j-1, ETA))
        rx_jm1 = RX[i, j - 1]
        rx_j = RX[i, j]
        ry_im1 = RY[i - 1, j]
        ry_i = RY[i, j]
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
        @views HA[1, :] .= 0.0
        @views HA[Ny1, :] .= 0.0
        @views HA[:, 1] .= 0.0
        @views HA[:, Nx1] .= 0.0
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
    return nothing
end # function compute_adiabatic_heating!

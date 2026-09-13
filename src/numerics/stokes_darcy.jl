"""
Assemble hydromechanical system of equations.

$(SIGNATURES)

# Details

    - ETA: viscosity at basic nodes
    - ETAP: viscosity at P nodes
    - GGG: shear modulus at basic nodes
    - GGGP: shear modulus at P nodes
    - SXY0: previous XY stress at basic nodes
    - SXX0: previous XX stress at P nodes
    - RHOX: total density at Vx nodes
    - RHOY: total density at Vy nodes
    - RHOFX: fluid density at Vx nodes
    - RHOFY: fluid density at Vy nodes
    - RX: ηfluid/Kϕ at Vx nodes
    - RY: ηfluid/Kϕ at Vy nodes
    - ETAPHI: bulk viscosity at P nodes
    - BETAPHI: bulk compressibility at P nodes
    - PHI: porosity at P nodes
    - gx: x gravitational acceleration at Vx nodes
    - gy: y gravitational acceleration at Vy nodes
    - pr0: previous total pressure at P nodes
    - pf0: previous fluid pressure at P nodes
    - DMP: mass transfer term at P nodes
    - dt: time step
    - R: vector to store RHS coefficients
    - coords: grid coordinates
    - L: optional ExtendableSparseMatrix buffer to reuse

# Returns

    - L: LHS coefficient matrix (SparseMatrixCSC)
"""
function assemble_hydromechanical_lse!(
    ETA,
    ETAP,
    GGG,
    GGGP,
    SXY0,
    SXX0,
    RHOX,
    RHOY,
    RHOFX,
    RHOFY,
    RX,
    RY,
    ETAPHI,
    BETAPHI,
    PHI,
    gx,
    gy,
    pr0,
    pf0,
    DMP,
    dt,
    R;
    betasolid=betasolid,
    betafluid=betafluid,
    phimin=phimin,
    phimax=phimax,
    hydrofracture::Bool=false,
    pr=nothing,
    pf=nothing,
    TEN=nothing,
    KX=nothing,
    KY=nothing,
    kappa_frac::Real=1.0e3,
    gamma_frac::Real=1.0,
    k_frac_max::Real=1.0e-9,
    coords=nothing,
    L=nothing,
    venting::Bool=false,
    venting_mode::Symbol=:darcy_sink,
    k_vent::Real=1.0e-11,
    conductance_factor::Real=1.0,
    ice_sealing::Bool=false,
    t_freeze::Real=273.15,
    dt_seal::Real=10.0,
    k_seal_min_ratio::Real=1.0e-6,
    rplanet::Real=50000.0,
    xcenter::Real=70000.0,
    ycenter::Real=70000.0,
    P_amb::Real=10.0,
    venting_species::Symbol=:H2O,
    tk=nothing,
    eta_fluid_surf::Real=1.0e-3,
    L_sub::Real=2.83e6,
    S_vent_out=nothing,
    DQPF::Union{AbstractMatrix{<:Real},Nothing}=nothing,
    workspace=nothing,
)
    Ny1, Nx1 = size(ETAP)
    Nx_val = Nx1 - 1
    Ny_val = Ny1 - 1
    dx_val = coords === nothing ? dx : coords.dx
    dy_val = coords === nothing ? dy : coords.dy

    # initialize or reuse LHS sparse coefficient matrix
    L_target =
        if workspace !== nothing &&
            hasproperty(workspace, :Ny1) &&
            hasproperty(workspace, :Nx1) &&
            hasproperty(workspace, :L) &&
            workspace.Ny1 == Ny1 &&
            workspace.Nx1 == Nx1
            workspace.L
        else
            L
        end

    L = if L_target === nothing
        ExtendableSparseMatrix(Nx1 * Ny1 * 6, Nx1 * Ny1 * 6)
    else
        if !isempty(L_target.cscmatrix.nzval)
            nonzeros(L_target.cscmatrix) .= zero(0.0)
        end
        L_target
    end
    # reset RHS coefficient vector
    R .= 0.0
    @inbounds begin
        for j in 1:1:Nx1, i in 1:1:Ny1
            # define global indices in algebraic space
            kvx = ((j-1)*Ny1 + i-1) * 6 + 1 # Vx solid
            kvy = kvx + 1 # Vy solid
            kpm = kvx + 2 # P total
            kqx = kvx + 3 # qx Darcy
            kqy = kvx + 4 # qy Darcy
            kpf = kvx + 5 # P fluid
            # Vx equation
            if i==1 || i==Ny1 || j==1 || j==Nx_val || j==Nx1
                # Vx equation external points: boundary conditions
                # all locations: ghost unknowns Vx₃=0 -> 1.0⋅Vx[i,j]=0.0
                updateindex!(L, +, 1.0, kvx, kvx)
                # R[kvx] = 0.0 # already done with initialization
                # left boundary
                if j == 1
                    R[kvx] = vxleft
                end
                # right boundary
                if j == Nx_val
                    R[kvx] = vxright
                end
                # top boundary
                if i==1 && 1<j<Nx_val
                    updateindex!(L, +, bctop, kvx, kvx+6)
                end
                # bottom boundary
                if i==Ny1 && 1<j<Nx_val
                    updateindex!(L, +, bcbottom, kvx, kvx-6)
                end
            else
                # Vx momentum internal stencil (see docs/src/explanations/discretization_numerics.md)
                # computational viscosity
                ETA₁ =
                    ETA[i - 1, j] * GGG[i - 1, j] * dt / (GGG[i - 1, j]*dt + ETA[i - 1, j])
                ETA₂ = ETA[i, j] * GGG[i, j] * dt / (GGG[i, j]*dt + ETA[i, j])
                ETAP₁ = ETAP[i, j] * GGGP[i, j] * dt / (GGGP[i, j]*dt + ETAP[i, j])
                ETAP₂ =
                    ETAP[i, j + 1] * GGGP[i, j + 1] * dt /
                    (GGGP[i, j + 1]*dt + ETAP[i, j + 1])
                # previous stresses
                SXY₁ = SXY0[i - 1, j] * ETA[i - 1, j] / (GGG[i - 1, j]*dt + ETA[i - 1, j])
                SXY₂ = SXY0[i, j] * ETA[i, j] / (GGG[i, j]*dt + ETA[i, j])
                SXX₁ = SXX0[i, j] * ETAP[i, j] / (GGGP[i, j]*dt + ETAP[i, j])
                SXX₂ =
                    SXX0[i, j + 1] * ETAP[i, j + 1] / (GGGP[i, j + 1]*dt + ETAP[i, j + 1])
                # density gradients
                ∂RHO∂x = 0.5 * (RHOX[i, j + 1] - RHOX[i, j - 1]) * inv(dx_val)
                ∂RHO∂y = 0.5 * (RHOX[i + 1, j] - RHOX[i - 1, j]) * inv(dy_val)
                # LHS coefficient matrix
                updateindex!(L, +, ETAP₁/dx_val^2, kvx, kvx-6*Ny1) # Vx₁
                updateindex!(L, +, ETA₁/dy_val^2, kvx, kvx-6) # Vx₂
                updateindex!(
                    L,
                    +,
                    (
                        -(ETAP₁+ETAP₂) * inv(dx_val^2) - (ETA₁+ETA₂) * inv(dy_val^2) -
                        ∂RHO∂x * gx[i, j] * dt
                    ),
                    kvx,
                    kvx,
                ) # Vx₃
                updateindex!(L, +, ETA₂/dy_val^2, kvx, kvx+6) # Vx₄
                updateindex!(L, +, ETAP₂/dx_val^2, kvx, kvx+6*Ny1) # Vx₅
                updateindex!(
                    L,
                    +,
                    (
                        ETAP₁ * inv(dx_val) * inv(dy_val) -
                        ETA₂ * inv(dx_val) * inv(dy_val) - ∂RHO∂y * gx[i, j] * dt * 0.25
                    ),
                    kvx,
                    kvy,
                ) # Vy₂
                updateindex!(
                    L,
                    +,
                    (
                        -ETAP₂ * inv(dx_val) * inv(dy_val) +
                        ETA₂ * inv(dx_val) * inv(dy_val) - ∂RHO∂y * gx[i, j] * dt * 0.25
                    ),
                    kvx,
                    kvy+6*Ny1,
                ) # Vy₄
                updateindex!(
                    L,
                    +,
                    (
                        -ETAP₁ * inv(dx_val) * inv(dy_val) +
                        ETA₁ * inv(dx_val) * inv(dy_val) - ∂RHO∂y * gx[i, j] * dt * 0.25
                    ),
                    kvx,
                    kvy-6,
                ) # Vy₁
                updateindex!(
                    L,
                    +,
                    (
                        ETAP₂ * inv(dx_val) * inv(dy_val) -
                        ETA₁ * inv(dx_val) * inv(dy_val) - ∂RHO∂y * gx[i, j] * dt * 0.25
                    ),
                    kvx,
                    kvy+6*Ny1-6,
                ) # Vy₃
                updateindex!(L, +, Kcont*inv(dx_val), kvx, kpm) # P₁
                updateindex!(L, +, -Kcont*inv(dx_val), kvx, kpm+6*Ny1) # P₂
                # RHS coefficient vector
                R[kvx] = (
                    -RHOX[i, j] * gx[i, j] - (SXY₂-SXY₁) * inv(dy_val) -
                    (SXX₂-SXX₁) * inv(dx_val)
                )
            end # Vx equation
            # Vy equation
            if i==1 || i==Ny_val || i==Ny1 || j==1 || j==Nx1
                # Vy equation external points: boundary conditions
                # all locations: ghost unknowns Vy₃=0 -> 1.0⋅Vy[i,j]=0.0
                updateindex!(L, +, 1.0, kvy, kvy)
                # R[kvy] = 0.0 # already done with initialization
                # top boundary
                if i == 1
                    R[kvy] = vytop
                end
                # bottom boundary
                if i == Ny_val
                    R[kvy] = vybottom
                end
                # left boundary
                if j==1 && 1<i<Ny_val
                    updateindex!(L, +, bcleft, kvy, kvy+6*Ny1)
                end
                # right boundary
                if j==Nx1 && 1<i<Ny_val
                    updateindex!(L, +, bcright, kvy, kvy-6*Ny1)
                end
            else
                # Vy momentum internal stencil (see docs/src/explanations/discretization_numerics.md)
                # computational viscosity
                ETA₁ =
                    ETA[i, j - 1] * GGG[i, j - 1] * dt / (GGG[i, j - 1]*dt + ETA[i, j - 1])
                ETA₂ = ETA[i, j] * GGG[i, j] * dt / (GGG[i, j]*dt + ETA[i, j])
                ETAP₁ = ETAP[i, j] * GGGP[i, j] * dt / (GGGP[i, j]*dt + ETAP[i, j])
                ETAP₂ =
                    ETAP[i + 1, j] * GGGP[i + 1, j] * dt /
                    (GGGP[i + 1, j]*dt + ETAP[i + 1, j])
                # previous stresses
                SXY₁ = SXY0[i, j - 1] * ETA[i, j - 1] / (GGG[i, j - 1]*dt + ETA[i, j - 1])
                SXY₂ = SXY0[i, j] * ETA[i, j] / (GGG[i, j]*dt + ETA[i, j])
                SYY₁ = -SXX0[i, j] * ETAP[i, j] / (GGGP[i, j]*dt + ETAP[i, j])
                SYY₂ =
                    -SXX0[i + 1, j] * ETAP[i + 1, j] / (GGGP[i + 1, j]*dt + ETAP[i + 1, j])
                # density gradients
                ∂RHO∂x = 0.5 * (RHOY[i, j + 1]-RHOY[i, j - 1]) / dx_val
                ∂RHO∂y = 0.5 * (RHOY[i + 1, j]-RHOY[i - 1, j]) / dy_val
                # LHS coefficient matrix
                updateindex!(L, +, ETA₁/dx_val^2, kvy, kvy-6*Ny1) # Vy₁
                updateindex!(L, +, ETAP₁/dy_val^2, kvy, kvy-6) # Vy₂
                updateindex!(
                    L,
                    +,
                    (
                        -(ETAP₁+ETAP₂) * inv(dy_val^2) - (ETA₁+ETA₂) * inv(dx_val^2) -
                        ∂RHO∂y * gy[i, j] * dt
                    ),
                    kvy,
                    kvy,
                ) # Vy₃
                updateindex!(L, +, ETAP₂ * inv(dy_val^2), kvy, kvy+6) # Vy₄
                updateindex!(L, +, ETA₂ * inv(dx_val^2), kvy, kvy+6*Ny1) # Vy₅
                updateindex!(
                    L,
                    +,
                    (
                        ETAP₁ * inv(dx_val) * inv(dy_val) -
                        ETA₂ * inv(dx_val) * inv(dy_val) - ∂RHO∂x * gy[i, j] * dt * 0.25
                    ),
                    kvy,
                    kvx,
                ) # Vx₃
                updateindex!(
                    L,
                    +,
                    (
                        -ETAP₂ * inv(dx_val) * inv(dy_val) +
                        ETA₂ * inv(dx_val) * inv(dy_val) - ∂RHO∂x * gy[i, j] * dt * 0.25
                    ),
                    kvy,
                    kvx+6,
                ) # Vx₄
                updateindex!(
                    L,
                    +,
                    (
                        -ETAP₁ * inv(dx_val) * inv(dy_val) +
                        ETA₁ * inv(dx_val) * inv(dy_val) - ∂RHO∂x * gy[i, j] * dt * 0.25
                    ),
                    kvy,
                    kvx-6*Ny1,
                ) # Vx₁
                updateindex!(
                    L,
                    +,
                    (
                        ETAP₂ * inv(dx_val) * inv(dy_val) -
                        ETA₁ * inv(dx_val) * inv(dy_val) - ∂RHO∂x * gy[i, j] * dt * 0.25
                    ),
                    kvy,
                    kvx+6-6*Ny1,
                ) # Vx₂
                updateindex!(L, +, Kcont*inv(dy_val), kvy, kpm) # P₁
                updateindex!(L, +, -Kcont*inv(dy_val), kvy, kpm+6) # P₂
                R[kvy] = (
                    -RHOY[i, j] * gy[i, j] - (SXY₂-SXY₁) * inv(dx_val) -
                    (SYY₂-SYY₁) * inv(dy_val)
                ) # RHS
            end # Vy equation
            # P equation
            if i==1 || i==Ny1 || j==1 || j==Nx1
                # P equation external points: boundary conditions
                # all locations: ghost unknowns P=0 -> 1.0⋅P[i,j]=0.0
                updateindex!(L, +, 1.0, kpm, kpm)
                # R[kpm] = 0.0 # already done with initialization
                # elseif i==j==2
            elseif (
                (i==2 && 2<=j<=Nx_val) ||
                (j==2 && 2<i<Ny_val) ||
                (i==Ny_val && 2<=j<=Nx_val) ||
                (j==Nx_val && 2<i<Ny_val)
            )
                # Ptotal/Pfluid real pressure boundary condition 'anchor'
                updateindex!(L, +, Kcont, kpm, kpm)
                R[kpm] = psurface
            else
                # Solid continuity internal stencil: ∂Vx/∂x + ∂Vy/∂y = 0
                updateindex!(L, +, -1.0/dx_val, kpm, kvx-6*Ny1) # Vx₁
                updateindex!(L, +, 1.0/dx_val, kpm, kvx) # Vx₂
                updateindex!(L, +, -1.0/dy_val, kpm, kvy-6) # Vy₁
                updateindex!(L, +, 1.0/dy_val, kpm, kvy) # Vy₂
                # Poroelastic continuity stencils based on simple3anpfl.m (Taras Gerya, pers. comm.)
                betadrained = compute_drained_compressibility(
                    BETAPHI[i, j], PHI[i, j], betasolid; phimin=phimin, phimax=phimax
                )
                kbw = compute_biot_willis_coefficient(betadrained, betasolid)
                ksk = compute_skempton_coefficient(
                    betadrained,
                    PHI[i, j],
                    betasolid,
                    betafluid;
                    phimin=phimin,
                    phimax=phimax,
                )

                # LHS coefficient matrix
                updateindex!(
                    L,
                    +,
                    Kcont * (inv(ETAPHI[i, j]) / (1.0 - PHI[i, j]) + betadrained / dt),
                    kpm,
                    kpm,
                ) # P: Ptotal
                updateindex!(
                    L,
                    +,
                    -Kcont *
                    (inv(ETAPHI[i, j]) / (1.0 - PHI[i, j]) + betadrained * kbw / dt),
                    kpm,
                    kpf,
                ) # P: Pfluid
                # RHS coefficient vector
                R[kpm] = (betadrained * (pr0[i, j] - kbw * pf0[i, j]) / dt + DMP[i, j])
            end # P equation
            # qxDarcy equation
            if i==1 || i==Ny1 || j==1 || j==Nx_val || j==Nx1
                # qxDarcy equation external points: boundary conditions
                # all locations: ghost unknowns qyD = 0 -> 1.0⋅qxD[i, j] = 0.0
                updateindex!(L, +, 1.0, kqx, kqx)
                # R[kqx] = 0.0 # already done with initialization
                # top boundary
                if i==1 && 1<j<Nx_val
                    updateindex!(L, +, bcftop, kqx, kqx+6)
                end
                # bottom boundary
                if i==Ny1 && 1<j<Nx_val
                    updateindex!(L, +, bcfbottom, kqx, kqx-6)
                end
            else
                # x-Darcy flux internal stencil: (η_f/k_ϕx) * qxD + ∂P/∂x = ρ_f * gx
                # See Stencil Topologies in docs/src/explanations/discretization_numerics.md
                # LHS coefficient matrix
                rx_val = RX[i, j]
                if hydrofracture && pr !== nothing && pf !== nothing && TEN !== nothing
                    Peff_x = 0.5 * (pr[i, j] + pr[i, j + 1] - pf[i, j] - pf[i, j + 1])
                    sigma_t_x = 0.5 * (TEN[i, j] + TEN[i - 1, j])
                    kphi_x = (KX !== nothing) ? KX[i, j] : 0.0
                    if kphi_x > 0.0
                        keff_x = compute_hydrofracture_permeability(
                            kphi_x,
                            Peff_x,
                            sigma_t_x;
                            active=true,
                            kappa_frac=kappa_frac,
                            gamma=gamma_frac,
                            kmax=k_frac_max,
                        )
                        rx_val = RX[i, j] * (kphi_x / keff_x)
                    else
                        ffrac_x = compute_hydrofracture_factor(
                            Peff_x,
                            sigma_t_x;
                            active=true,
                            kappa_frac=kappa_frac,
                            gamma=gamma_frac,
                        )
                        rx_floor = 1.0e-5 / k_frac_max
                        rx_val = max(RX[i, j] / ffrac_x, rx_floor)
                    end
                end
                updateindex!(L, +, rx_val, kqx, kqx) # qxD
                updateindex!(L, +, -Kcont*inv(dx_val), kqx, kpf) # P₁
                updateindex!(L, +, Kcont*inv(dx_val), kqx, kpf+6*Ny1) # P₂
                # RHS coefficient vector
                R[kqx] = RHOFX[i, j] * gx[i, j]
            end # qxDarcy equation
            # qyDarcy equation
            if i==1 || i==Ny_val || i==Ny1 || j==1 || j==Nx1
                # qyDarcy equation external points: boundary conditions
                # all locations: ghost unknowns qyD = 0 -> 1.0⋅qyD[i, j] = 0.0
                updateindex!(L, +, 1.0, kqy, kqy)
                # R[kqy] = 0.0 # already done with initialization
                # left boundary
                if j==1 && 1<i<Ny_val
                    updateindex!(L, +, bcfleft, kqy, kqy+6*Ny1)
                end
                # right boundary
                if j==Nx1 && 1<i<Ny_val
                    updateindex!(L, +, bcfright, kqy, kqy-6*Ny1)
                end
            else
                # y-Darcy flux internal stencil: (η_f/k_ϕy) * qyD + ∂P/∂y = ρ_f * gy
                # LHS coefficient matrix
                ry_val = RY[i, j]
                if hydrofracture && pr !== nothing && pf !== nothing && TEN !== nothing
                    Peff_y = 0.5 * (pr[i, j] + pr[i + 1, j] - pf[i, j] - pf[i + 1, j])
                    sigma_t_y = 0.5 * (TEN[i, j] + TEN[i, j - 1])
                    kphi_y = (KY !== nothing) ? KY[i, j] : 0.0
                    if kphi_y > 0.0
                        keff_y = compute_hydrofracture_permeability(
                            kphi_y,
                            Peff_y,
                            sigma_t_y;
                            active=true,
                            kappa_frac=kappa_frac,
                            gamma=gamma_frac,
                            kmax=k_frac_max,
                        )
                        ry_val = RY[i, j] * (kphi_y / keff_y)
                    else
                        ffrac_y = compute_hydrofracture_factor(
                            Peff_y,
                            sigma_t_y;
                            active=true,
                            kappa_frac=kappa_frac,
                            gamma=gamma_frac,
                        )
                        ry_floor = 1.0e-5 / k_frac_max
                        ry_val = max(RY[i, j] / ffrac_y, ry_floor)
                    end
                end
                updateindex!(L, +, ry_val, kqy, kqy) # qyD
                updateindex!(L, +, -Kcont*inv(dy_val), kqy, kpf) # P₁
                updateindex!(L, +, Kcont*inv(dy_val), kqy, kpf+6) # P₂
                # RHS coefficient vector
                R[kqy] = RHOFY[i, j] * gy[i, j]
            end # qyDarcy equation
            # Ptotal/Pfluid equation 
            if i==1 || i==Ny1 || j==1 || j==Nx1
                # Ptotal/Pfluid equation external points: boundary conditions
                # all locations: ghost unknowns P = 0 -> 1.0⋅P[i, j] = 0.0
                updateindex!(L, +, 1.0, kpf, kpf)
                # R[kpf] = 0.0 # already done with initialization
                # elseif i==j==2
            elseif (
                (i==2 && 2<=j<=Nx_val) ||
                (j==2 && 2<i<Ny_val) ||
                (i==Ny_val && 2<=j<=Nx_val) ||
                (j==Nx_val && 2<i<Ny_val)
            )
                # Ptotal/Pfluid real pressure boundary condition 'anchor'
                updateindex!(L, +, Kcont, kpf, kpf)
                R[kpf] = psurface
            else
                # Fluid continuity internal stencil: ∂qxD/∂x + ∂qyD/∂y - (Pt - Pf)/η_ϕ = 0
                # LHS coefficient matrix
                updateindex!(L, +, -inv(dx_val), kpf, kqx-6*Ny1) # qxD₁
                updateindex!(L, +, inv(dx_val), kpf, kqx) # qxD₂
                updateindex!(L, +, -inv(dy_val), kpf, kqy-6) # qyD₁
                updateindex!(L, +, inv(dy_val), kpf, kqy) # qyD₂

                # LHS coefficient matrix
                updateindex!(
                    L,
                    +,
                    -Kcont *
                    (inv(ETAPHI[i, j]) / (1.0 - PHI[i, j]) + betadrained * kbw / dt),
                    kpf,
                    kpm,
                ) # Ptotal
                updateindex!(
                    L,
                    +,
                    Kcont *
                    (inv(ETAPHI[i, j]) / (1.0 - PHI[i, j]) + betadrained * kbw / ksk / dt),
                    kpf,
                    kpf,
                ) # Pfluid
                # RHS coefficient vector
                R[kpf] = -betadrained * kbw * (pr0[i, j] - (1.0 / ksk) * pf0[i, j]) / dt
                if DQPF !== nothing
                    R[kpf] += DQPF[i, j]
                end
            end # Ptotal/Pfluid equation
        end # for j=1:1:Nx1, i=1:1:Ny1
    end # @inbounds 

    if venting && tk !== nothing && coords !== nothing
        apply_venting_surface_boundary!(
            L,
            R,
            tk,
            coords,
            rplanet,
            xcenter,
            ycenter,
            P_amb;
            species=venting_species,
            k_vent=k_vent,
            conductance_factor=conductance_factor,
            mode=venting_mode,
            hydrofracture=hydrofracture,
            ice_sealing=ice_sealing,
            t_freeze=t_freeze,
            dt_seal=dt_seal,
            k_seal_min_ratio=k_seal_min_ratio,
            kappa_frac=kappa_frac,
            gamma_frac=gamma_frac,
            k_frac_max=k_frac_max,
            pr=pr,
            pf=pf,
            TEN=TEN,
            PHI=PHI,
            phimin=phimin,
            dt=dt,
            eta_fluid_surf=eta_fluid_surf,
            L_sub=L_sub,
            Kcont=Kcont,
            S_vent_out=S_vent_out,
        )
    end

    flush!(L) # finalize CSC matrix
    if workspace !== nothing &&
        hasproperty(workspace, :Ny1) &&
        hasproperty(workspace, :Nx1) &&
        hasproperty(workspace, :is_initialized) &&
        workspace.Ny1 == Ny1 &&
        workspace.Nx1 == Nx1
        workspace.is_initialized = true
    end
    # return L
    return L.cscmatrix
end # function assemble_hydromechanical_lse!

"""
Process hydromechanical solution vector to output physical observables.

$(SIGNATURES)

# Details

    - S: hydromechanical solution vector
    - vx: solid velocity at Vx nodes
    - vy: solid velocity at Vy nodes
    - pr: total pressure at P nodes
    - qxD: qx-Darcy flux at Vx nodes
    - qyD: qy-Darcy flux at Vy nodes
    - pf: fluid pressure at P nodes

# Returns

    - nothing
"""
function process_hydromechanical_solution!(S, vx, vy, pr, qxD, qyD, pf; coords=nothing)
    Ny1, Nx1 = size(vx)
    S_mat = reshape(S, (:, Ny1, Nx1))
    @inbounds begin
        @views @. vx = S_mat[1, :, :]
        @views @. vy = S_mat[2, :, :]
        @views @. pr = S_mat[3, :, :] .* Kcont
        @views @. qxD = S_mat[4, :, :]
        @views @. qyD = S_mat[5, :, :]
        @views @. pf = S_mat[6, :, :] .* Kcont
    end # @inbounds
    # Δp = 0.25 * (pf[2, 2]+pf[2, Nx]+pf[Ny, 2]+pf[Ny, Nx]) - psurface
    # pr .-= Δp
    # pf .-= Δp
    return nothing
end # function process_hydromechanical_solution!

"""
Recompute bulk viscosity at P nodes.

# Details

    - ETA: viscoplastic viscosity at basic nodes
    - ETAP: viscosity at P nodes
    - ETAPHI: bulk viscosity at P nodes
    - PHI: porosity at P nodes
    - etaphikoef: coefficient: shear viscosity -> compaction viscosity

# Returns

    - nothing
"""
function recompute_bulk_viscosity!(ETA, ETAP, ETAPHI, PHI, etaphikoef)
    @inbounds begin
        @views @. ETAP[2:(end - 1), 2:(end - 1)] =
            4.0 / (
                inv(ETA[1:(end - 1), 1:(end - 1)]) +
                inv(ETA[2:end, 1:(end - 1)]) +
                inv(ETA[1:(end - 1), 2:end]) +
                inv(ETA[2:end, 2:end])
            )
        @views @. ETAPHI = etaphikoef * ETAP * inv(PHI)
    end # @inbounds
    return nothing
end

"""
Compute porosity coefficient Aϕ = Dln[(1-ϕ)/ϕ]/Dt

$(SIGNATURES)

# Details

## In

    - ETAPHI: bulk viscosity at P Nodes
    - BETAPHI: bulk compressibility at P nodes
    - PHI: porosity at P Nodes
    - pr: total pressure at P nodes
    - pf: fluid pressure at P nodes
    - pr0: previous step total pressure at P nodes
    - pf0: previous step fluid pressure at P nodes
    - dt: time step

## Out

    - APHI: porosity coefficient at P nodes

# Returns

    - aphimax: maximum absolute porosity coefficient
"""
function compute_Aϕ!(
    APHI,
    ETAPHI,
    BETAPHI,
    PHI,
    pr,
    pf,
    pr0,
    pf0,
    dt;
    coords=nothing,
    betasolid=betasolid,
    phimin=phimin,
    phimax=phimax,
    S_vent::Union{AbstractMatrix{Float64},Nothing}=nothing,
)
    # APHI .= 0.0
    Ny1, Nx1 = size(APHI)
    Nx = Nx1 - 1
    Ny = Ny1 - 1
    @inbounds begin
        for j in 2:Nx, i in 2:Ny
            betadrained = compute_drained_compressibility(
                BETAPHI[i, j], PHI[i, j], betasolid; phimin=phimin, phimax=phimax
            )
            kbw = compute_biot_willis_coefficient(betadrained, betasolid)
            compaction = (
                (pr[i, j] - pf[i, j]) / (ETAPHI[i, j] * (1.0 - PHI[i, j])) +
                betadrained * ((pr[i, j] - pr0[i, j]) - kbw * (pf[i, j] - pf0[i, j])) / dt
            )
            if S_vent !== nothing
                s_v = S_vent[i, j]
                if s_v > 0.0
                    # Subtract venting sink to avoid double drainage with explicit marker sink.
                    compaction = max(0.0, compaction - s_v)
                end
            end
            APHI[i, j] = compaction / PHI[i, j]
        end
        return maximum(abs, @view APHI[2:Ny, 2:Nx]) # includes [2, 2] anchor abberation
    end # @inbounds
    # return maximum(abs, APHI[3:Ny-1, 3:Nx-1]) # no abberation
end # function compute_Aϕ!

"""
Compute current fluid velocities.

$(SIGNATURES)

# Details

## In

    - PHIX: porosity at Vx nodes
    - PHIY: porosity at Vy nodes
    - qxD: qx-Darcy flux at Vx nodes
    - qyD: qy-Darcy flux at Vy nodes
    - vx: solid velocity at Vx nodes
    - vy: solid velocity at Vy nodes

## Out 

    - vxf: fluid vx velocity at Vx nodes
    - vyf: fluid vy velocity at Vy nodes

# Returns

    - nothing
"""
function compute_fluid_velocities!(PHIX, PHIY, qxD, qyD, vx, vy, vxf, vyf; coords=nothing)
    Ny1, Nx1 = size(vxf)
    Nx = Nx1 - 1
    Ny = Ny1 - 1
    @inbounds begin
        # vx velocity
        @views @. vxf[2:Ny, 1:Nx] = qxD[2:Ny, 1:Nx] / PHIX[2:Ny, 1:Nx]
        # top boundary
        @views @. vxf[1, :] = -bcftop * vxf[2, :]
        # bottom boundary
        @views @. vxf[Ny1, :] = -bcfbottom * vxf[Ny, :]
        # vy velocity
        @views @. vyf[1:Ny, 2:Nx] = qyD[1:Ny, 2:Nx] / PHIY[1:Ny, 2:Nx]
        # left boundary
        @views @. vyf[:, 1] = -bcfleft * vyf[:, 2]
        # right boundary
        @views @. vyf[:, Nx1] = -bcfright * vyf[:, Nx]
        # adding solid velocity
        @views @. vxf += vx
        @views @. vyf += vy
    end # @inbounds

    return nothing
end # function compute_fluid_velocities!

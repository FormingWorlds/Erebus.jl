"""
Compute viscosities, stresses, and density gradients
for hydromechanical solver.

$(SIGNATURES)

# Details

## In

    - ETA: viscoplastic viscosity at basic nodes
    - ETAP: viscosity at P nodes
    - GGG: shear modulus at basic nodes
    - GGGP: shear modulus at P nodes
    - SXY0: σ₀xy XY stress at basic nodes
    - SXX0:σ₀xy XY stress at basic nodes
    - RHOX: density at Vx nodes
    - RHOY: density at Vy nodes
    - dt: time step

## Out 

    - ETAcomp: computational viscosity at basic nodes
    - ETAPcomp: computational viscosity at P nodes
    - SXYcomp: previous XY stresses at basic nodes
    - SXXcomp: previous XX stresses at P nodes
    - SYYcomp: previous YY stresses at P nodes
    - dRHOXdx: density gradient at Vx nodes in x direction
    - dRHOXdy: density gradient at Vx nodes in y direction
    - dRHOYdx: density gradient at Vy nodes in x direction
    - dRHOYdy: density gradient at Vy nodes in y direction

# Returns

    - nothing
"""
function get_viscosities_stresses_density_gradients!(
    ETA,
    ETAP,
    GGG,
    GGGP,
    SXY0,
    SXX0,
    RHOX,
    RHOY,
    dt,
    ETAcomp,
    ETAPcomp,
    SXYcomp,
    SXXcomp,
    SYYcomp,
    dRHOXdx,
    dRHOXdy,
    dRHOYdx,
    dRHOYdy;
    coords=nothing,
)
    Ny1, Nx1 = size(RHOX)
    Nx_val = Nx1 - 1
    Ny_val = Ny1 - 1
    @unpack_coords coords dx dy
    # computational viscosity
    @views @. ETAcomp = ETA*GGG*dt / (GGG*dt + ETA)
    @views @. ETAPcomp = ETAP*GGGP*dt / (GGGP*dt + ETAP)
    # previous stresses
    @views @. SXYcomp = SXY0*ETA / (GGG*dt+ETA)
    @views @. SXXcomp = SXX0 * ETAP / (GGGP * dt + ETAP)
    @views @. SYYcomp = -SXX0 * ETAP / (GGGP * dt + ETAP)
    # density gradients
    @inbounds begin
        @views @. dRHOXdx[:, 2:Nx_val] =
            0.5 * (RHOX[:, 3:Nx1]-RHOX[:, 1:(Nx1 - 2)]) * inv(dx_val)
        @views @. dRHOXdy[2:Ny_val, :] =
            0.5 * (RHOX[3:Ny1, :]-RHOX[1:(Ny1 - 2), :]) * inv(dy_val)
        @views @. dRHOYdx[:, 2:Nx_val] =
            0.5 * (RHOY[:, 3:Nx1]-RHOY[:, 1:(Nx1 - 2)]) * inv(dx_val)
        @views @. dRHOYdy[2:Ny_val, :] =
            0.5 * (RHOY[3:Ny1, :]-RHOY[1:(Ny1 - 2), :]) * inv(dy_val)
    end # @inbounds
    return nothing
end # function get_viscosities_stresses_density_gradients!

"""
Check if point (i, j) lies on or outside domain boundary for gravitational Poisson solver.

$(SIGNATURES)

# Details

    - i: y-index
    - j: x-index
    - Ny1: total number of P nodes in y
    - Nx1: total number of P nodes in x
    - xp_val: x coordinates of P nodes
    - yp_val: y coordinates of P nodes
    - xc_val: center coordinate in x
    - yc_val: center coordinate in y
    - r_limit: radius limit for circular planetary domain

# Returns

    - true if node is boundary or outside domain, false if interior
"""
@inline function is_gravitational_boundary(
    i::Integer,
    j::Integer,
    Ny1::Integer,
    Nx1::Integer,
    xp_val,
    yp_val,
    xc_val,
    yc_val,
    r_limit,
)
    return (
        i == 1 ||
        i == Ny1 ||
        j == 1 ||
        j == Nx1 ||
        distance(xp_val[j], yp_val[i], xc_val, yc_val) > r_limit
    )
end

"""
Compute stress, stress change, and strain rate components.

$(SIGNATURES)

# Details

## In

    - vx: solid vx velocity at Vx nodes
    - vy: solid vy velocity at Vy nodes
    - ETA: viscosity at basic nodes
    - GGG: shear modulus at basic nodes
    - ETAP: viscosity at P nodes
    - GGGP: shear modulus at P nodes
    - SXX0: previous time step σ₀′xx at P nodes
    - SXY0: previous time step σ₀xy at basic nodes
    - dt: computational time step

## Out

    - EXX: ϵxx at P nodes
    - EXY: ϵxy at basic nodes
    - SXX: σ′xx P nodes
    - SXY: σxy at basic nodes
    - DSXX: stress change Δσ′xx at P nodes
    - DSXY: stress change Δσxy at basic nodes
    - EII: second strain rate invariant ϵᴵᴵ at P nodes
    - SII: second stress invariant σᴵᴵ at P nodes

# Returns
    
        - nothing
"""
function compute_stress_strainrate!(
    vx,
    vy,
    ETA,
    GGG,
    ETAP,
    GGGP,
    SXX0,
    SXY0,
    EXX,
    EXY,
    SXX,
    SXY,
    DSXX,
    DSXY,
    EII,
    SII,
    dt;
    coords=nothing,
)
    Ny, Nx = size(EXY)
    @unpack_coords coords dx dy
    @inbounds begin
        # ϵxy, σxy, Δσxy at basic nodes
        for j in 1:1:Nx, i in 1:1:Ny
            EXY[i, j] =
                0.5 * ((vx[i + 1, j]-vx[i, j]) / dy_val + (vy[i, j + 1]-vy[i, j])/dx_val)
            SXY[i, j] = (
                2*ETA[i, j]*EXY[i, j]*GGG[i, j]*dt / (GGG[i, j]*dt+ETA[i, j]) +
                SXY0[i, j]*ETA[i, j] / (GGG[i, j]*dt+ETA[i, j])
            )
            DSXY[i, j] = SXY[i, j] - SXY0[i, j]
        end
        # ϵxx, σ′xx, Δσ'xx and Eᴵᴵ, Sᴵᴵ at P nodes
        for j in 2:1:Nx, i in 2:1:Ny
            EXX[i, j] =
                0.5 * ((vx[i, j]-vx[i, j - 1]) / dx_val - (vy[i, j]-vy[i - 1, j]) / dy_val)
            SXX[i, j] = (
                2*ETAP[i, j]*EXX[i, j]*GGGP[i, j]*dt / (GGGP[i, j]*dt+ETAP[i, j]) +
                SXX0[i, j]*ETAP[i, j] / (GGGP[i, j]*dt+ETAP[i, j])
            )
            DSXX[i, j] = SXX[i, j] - SXX0[i, j]
            EII[i, j] = sqrt(EXX[i, j]^2 + grid_average(i-1, j-1, EXY)^2)
            SII[i, j] = sqrt(SXX[i, j]^2 + grid_average(i-1, j-1, SXY)^2)
        end
    end # @inbounds
    return nothing
end # function compute_stress_strainrate!

"""
Apply symmetry to P node observables.

$(SIGNATURES)

# Details

    - SXX: σ′xx at P nodes
    - APHI: Aϕ = Dln[(1-ϕ)/ϕ]/Dt at P nodes
    - PHI: porosity at P nodes
    - pr: total pressure at P nodes
    - pf: fluid pressure at P nodes
    - ps: solid pressure at P nodes

# Returns

    - nothing
"""
function symmetrize_p_node_observables!(SXX, APHI, PHI, pr, pf, ps)
    Ny1, Nx1 = size(SXX)
    Nx = Nx1 - 1
    Ny = Ny1 - 1
    # top boundary
    @inbounds @views @. begin
        SXX[1, 2:Nx] = SXX[2, 2:Nx]
        APHI[1, 2:Nx] = APHI[2, 2:Nx]
        PHI[1, 2:Nx] = PHI[2, 2:Nx]
        pr[1, 2:Nx] = pr[2, 2:Nx]
        pf[1, 2:Nx] = pf[2, 2:Nx]
        # bottom boundary
        SXX[Ny1, 2:Nx] = SXX[Ny, 2:Nx]
        APHI[Ny1, 2:Nx] = APHI[Ny, 2:Nx]
        PHI[Ny1, 2:Nx] = PHI[Ny, 2:Nx]
        pr[Ny1, 2:Nx] = pr[Ny, 2:Nx]
        pf[Ny1, 2:Nx] = pf[Ny, 2:Nx]
        # left boundary
        SXX[:, 1] = SXX[:, 2]
        APHI[:, 1] = APHI[:, 2]
        PHI[:, 1] = PHI[:, 2]
        pr[:, 1] = pr[:, 2]
        pf[:, 1] = pf[:, 2]
        # right boundary
        SXX[:, Nx1] = SXX[:, Nx]
        APHI[:, Nx1] = APHI[:, Nx]
        PHI[:, Nx1] = PHI[:, Nx]
        pr[:, Nx1] = pr[:, Nx]
        pf[:, Nx1] = pf[:, Nx]
        # solid pressure
        ps = (pr-pf*PHI) * inv(1-PHI)
    end
    return nothing
end # function symmetrize_p_node_observables!

"""
Compute nodal adjustment and return plastic iterations completeness status.

$(SIGNATURES)

# Details

    - ETA: viscoplastic viscosity at basic nodes
    - ETA0: previous time step viscoplastic viscosity at basic nodes
    - ETA5: plastic iterations viscoplastic viscosity at basic nodes
    - GGG: shear modulus at basic nodes
    - SXX: σ′xx at P nodes
    - SXY: σxy at basic nodes
    - pr: total pressure at P nodes
    - pf: fluid pressure at P nodes
    - COH: compressive strength at basic nodes 
    - TEN: tensile strength at basic nodes 
    - FRI: friction at basic nodes
    - YNY: plastic yielding status at basic nodes 
    - YNY5: plastic iterations plastic yielding status at basic nodes
    - YERRNOD: vector of summed yielding errors of nodes over plastic iterations
    - DSY: (SIIB-syield) at basic nodes
    - dt: time set
    - iplast: plastic iteration step 

# Returns

    - plastic_iterations_complete: true if plastic iterations complete
"""
function compute_nodal_adjustment!(
    ETA,
    ETA0,
    ETA5,
    GGG,
    SXX,
    SXY,
    pr,
    pf,
    COH,
    TEN,
    FRI,
    YNY,
    YNY5,
    YERRNOD,
    DSY,
    dt,
    iplast;
    etawt::Real=0.0,
    etamax::Real=1e23,
    etamin::Real=1e12,
    yerrmax::Real=1e2,
    nplast::Int=100_000,
)
    # reset / setup
    Ny, Nx = size(ETA)
    ETA5 .= ETA0
    YNY5 .= 0
    DSY .= 0.0
    ynpl = 0
    ddd = 0.0
    @inbounds begin
        for j in 1:1:Nx, i in 1:1:Ny
            # second stress invariant at basic nodes
            SIIB = sqrt(SXY[i, j]^2 + grid_average(i, j, SXX)^2)
            # second invariant for purely elastic stress buildup at basic nodes
            siiel = SIIB * (GGG[i, j]*dt+ETA[i, j]) / ETA[i, j]
            # interpolate total and fluid pressure at basic nodes
            prB = grid_average(i, j, pr)
            pfB = grid_average(i, j, pf)
            # Yielding stress: confined and tensile fracture.
            # Note: uses Terzaghi effective stress (prB - pfB) for frictional shear and tensile failure,
            # following standard rock mechanics (Terzaghi 1943, Handin et al. 1963).
            syieldc = COH[i, j] + FRI[i, j] * (prB-pfB)
            syieldt = TEN[i, j] + (prB-pfB)
            # non-negative yielding stress requirement
            syield = max(min(syieldc, syieldt), 0.0)
            # update error for previous yielding nodes
            ynn = false
            if YNY[i, j] > 0
                ynn = true
                DSY[i, j] = SIIB - syield
                ddd += DSY[i, j]^2
                ynpl += 1
            end
            # correcting viscosity for yielding
            if syield < siiel
                # update viscosity for basic node
                etapl = dt * GGG[i, j] * syield/(siiel-syield)
                if etapl < ETA0[i, j]
                    # recompute nodal viscosity, apply min/max viscosity cutoffs
                    ETA5[i, j] = etapl^(1.0-etawt) * ETA[i, j]^etawt
                    if ETA5[i, j] > etamax
                        ETA5[i, j] = etamax
                    elseif ETA5[i, j] < etamin
                        ETA5[i, j] = etamin
                    end
                    # mark yielding nodes
                    YNY5[i, j] = 1
                    # update error for new yielding nodes
                    if ynn == false
                        DSY[i, j] = SIIB - syield
                        ddd += DSY[i, j]^2
                        ynpl += 1
                    end
                end
            end
        end
        if ynpl > 0
            YERRNOD[iplast] = sqrt(ddd/ynpl)
        end
        # return plastic iteration completeness
        @info "end plastic iter $iplast: ynpl=$ynpl, YERRNOD=$(YERRNOD[iplast])"
        return ynpl==0 || YERRNOD[iplast]<yerrmax || iplast==nplast
    end # @inbounds
end # function compute_nodal_adjustment!

"""
Compare two arrays of identical sizes element-wise and fill a third array with the larger value if positive and zero otherwise.

$(SIGNATURES)

# Details

    - A: first array
    - B: second array
    - C: result array

# Returns

    - nothing
"""
function positive_max!(A, B, C)
    @inbounds for i in eachindex(A)
        C[i] = max(0, ifelse(A[i] > B[i], A[i], B[i]))
    end
    return nothing
end # function positive_max

"""
Decide next pass plastic iteration time step, viscoplastic viscosity,
and basic node yielding status.

$(SIGNATURES)

# Details:

    - ETA: viscoplastic viscosity at basic nodes
    - ETA5: plastic iterations viscoplastic viscosity at basic nodes
    - ETA00: previous time step viscoplastic viscosity at basic nodes
    - YNY: plastic yielding status at basic nodes 
    - YNY5: plastic iterations plastic yielding status at basic nodes
    - YNY00: previous time step plastic yielding status at basic nodes
    - YNY_inv_ETA: inverse of plastic viscosity at yielding basic nodes
    - dt: current time step
    - iplast: current plastic iteration counter

# Returns

    - dt: adjusted next time step
"""
function finalize_plastic_iteration_pass!(
    ETA, ETA5, ETA00, YNY, YNY5, YNY00, YNY_inv_ETA, dt, iplast;
    dtstep::Int=200, dtcoefdn::Real=0.5
)
    if iplast % dtstep == 0
        # dtstep plastic iterations performed without reaching targets:
        # decrease time step and reset to previous viscoplastic viscosity
        dt *= dtcoefdn
        @info "reducing dt due to plastic iteration limit: dt=$dt s"
        ETA .= ETA00
        YNY .= YNY00
    else
        # perform next plastic iteration pass with new viscoplastic viscosity
        ETA .= ETA5
        YNY .= YNY5
    end
    @views @. YNY_inv_ETA = YNY / ETA
    return dt
end # function finalize_plastic_iteration_pass

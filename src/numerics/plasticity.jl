"""
Exception thrown when plastic iterations fail to converge within `max_plastic_iterations`
after exceeding the maximum allowable timestep reductions (`max_dt_reductions`).

$(FIELDS)
"""
struct PlasticConvergenceError <: Exception
    step::Int
    residual::Float64
    dt::Float64
    msg::String
end

function PlasticConvergenceError(
    step::Int,
    residual::Float64,
    dt::Float64;
    msg::String="Plastic iterations failed to converge",
)
    return PlasticConvergenceError(step, residual, dt, msg)
end

function Base.showerror(io::IO, e::PlasticConvergenceError)
    return print(
        io,
        "PlasticConvergenceError: ",
        e.msg,
        " at step ",
        e.step,
        " (residual: ",
        e.residual,
        ", dt: ",
        e.dt,
        " s)",
    )
end

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
    coords::GridCoordinates=default_grid_coordinates(),
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
    coords::GridCoordinates=default_grid_coordinates(),
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
    max_plastic_iterations::Int=10_000,
    nplast::Union{Int,Nothing}=nothing,
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
        return ynpl == 0 || YERRNOD[iplast] < yerrmax
    end # @inbounds
end # function compute_nodal_adjustment!

"""
    PlasticCycleTracker(; max_period=10)

Per-pass record of the plastic loop state for the cycle stop of the plastic loop. A ring of
state copies long enough to hold a whole run, the scale of each field-scaled array, the
residual and the time step of each pass, and one run counter per period `p`: the number of
consecutive passes that match the earliest pass of the same phase in the run.
[`reset_plastic_cycle!`](@ref) sizes the ring for a window.

# Arguments

    - max_period: longest period checked, in passes

# Returns

    - an empty tracker

# Throws

    - DomainError: `max_period < 1`
"""
mutable struct PlasticCycleTracker
    ring::Vector{Vector{Float64}}
    scale::Matrix{Float64}
    lens::Vector{Int}
    nnode::Int
    Y::Vector{Float64}
    dt::Vector{Float64}
    run::Vector{Int}
    npass::Int
    window::Int
end
function PlasticCycleTracker(; max_period::Integer=10)
    max_period >= 1 || throw(DomainError(max_period, "max_period must be at least 1"))
    return PlasticCycleTracker(
        Vector{Float64}[],
        zeros(0, 0),
        Int[],
        0,
        Float64[],
        Float64[],
        zeros(Int, max_period),
        0,
        0,
    )
end

"""
    reset_plastic_cycle!(tr, fields, scaled, window) -> tr

Clear the pass count and the run counters of `tr` and size its ring for the arrays and the
window of a plastic loop, before the first pass of the loop. The ring holds
`max(window, 3 max_period) + max_period + 1` state copies, so the earliest pass of a run is
still kept when the run ends; it holds none when `window` is 0.

# Arguments

    - tr: [`PlasticCycleTracker`](@ref)
    - fields: tuple of arrays tested node by node
    - scaled: tuple of arrays tested against their field scale
    - window: minimum run length in passes; 0 disables the test

# Returns

    - tr

# Throws

    - DomainError: `window < 0`
"""
function reset_plastic_cycle!(
    tr::PlasticCycleTracker, fields::Tuple, scaled::Tuple, window::Integer
)
    window >= 0 || throw(DomainError(window, "window must be >= 0"))
    np = length(tr.run)
    nslot = window > 0 ? max(window, 3 * np) + np + 1 : 0
    lens = Int[length(f) for f in (fields..., scaled...)]
    nstate = sum(lens; init=0)
    nold = length(tr.ring)
    resize!(tr.ring, nslot)
    for i in 1:nslot
        i > nold && (tr.ring[i] = Float64[])
        length(tr.ring[i]) == nstate || resize!(tr.ring[i], nstate)
    end
    size(tr.scale) == (length(scaled), nslot) || (tr.scale = zeros(length(scaled), nslot))
    resize!(tr.Y, nslot)
    resize!(tr.dt, nslot)
    tr.lens = lens
    tr.nnode = sum(length, fields; init=0)
    fill!(tr.run, 0)
    tr.npass = 0
    tr.window = window
    return tr
end

"""
    _cycle_match(a, b, rtol, s=0.0) -> Bool

True when `a == b`, or when both are finite and `abs(a - b) <= rtol * max(abs(a), abs(b), s)`.

# Arguments

    - a, b: values of one quantity at two passes
    - rtol: relative tolerance
    - s: scale of the field the values belong to (same unit as `a`); 0 for a per-node test

# Returns

    - true when the values match
"""
function _cycle_match(a, b, rtol, s=0.0)
    return a == b ||
           (isfinite(a) && isfinite(b) && abs(a - b) <= rtol * max(abs(a), abs(b), s))
end

"""
    _finite_absmax(v, r) -> Float64

Largest finite magnitude of `v[r]`, or 0 when `r` is empty or holds no finite value.

# Arguments

    - v: vector
    - r: index range into `v`

# Returns

    - the largest `abs(v[i])` over finite `v[i]`, `i in r`
"""
function _finite_absmax(v, r)
    m = 0.0
    @inbounds for i in r
        isfinite(v[i]) && (m = max(m, abs(v[i])))
    end
    return m
end

"""
    _pack_cycle_state!(cur, off, fields) -> Int

Copy the arrays of `fields` one after the other into `cur`, starting after index `off`.

# Arguments

    - cur: ring slot of the current pass
    - off: number of values already in `cur`
    - fields: tuple of arrays

# Returns

    - the number of values in `cur` after the copy
"""
_pack_cycle_state!(cur, off, ::Tuple{}) = off
function _pack_cycle_state!(cur, off, fields::Tuple)
    f = first(fields)
    copyto!(cur, off + 1, f, 1, length(f))
    return _pack_cycle_state!(cur, off + length(f), Base.tail(fields))
end

"""
    _cycle_lag_match(tr, cur, prev, ic, ip, rtol) -> Bool

True when every per-node value of `cur` matches `prev` node by node and every value of a
field-scaled array matches with the larger scale of that array in slots `ic` and `ip`.

# Arguments

    - tr: [`PlasticCycleTracker`](@ref)
    - cur, prev: ring slots of the current pass and of the pass `p` before
    - ic, ip: slot indices of `cur` and `prev`
    - rtol: relative tolerance

# Returns

    - true when the two states match
"""
function _cycle_lag_match(tr::PlasticCycleTracker, cur, prev, ic, ip, rtol)
    @inbounds for i in 1:tr.nnode
        _cycle_match(cur[i], prev[i], rtol) || return false
    end
    off = tr.nnode
    nfield = length(tr.lens) - size(tr.scale, 1)
    @inbounds for q in axes(tr.scale, 1)
        s = max(tr.scale[q, ic], tr.scale[q, ip])
        for i in (off + 1):(off + tr.lens[nfield + q])
            _cycle_match(cur[i], prev[i], rtol, s) || return false
        end
        off += tr.lens[nfield + q]
    end
    return true
end

"""
    plastic_cycle_record!(tr, fields, scaled, Y, dt, window; rtol=1.0e-9) -> Bool

Record one plastic pass and return true when the loop state is periodic. The state of a
pass is the arrays in `fields` (viscosity, yield flags, Darcy resistances), the arrays in
`scaled` (velocities, pressures, Darcy fluxes), the residual `Y` and the time step `dt`. A
value `a` of `fields`, `Y` or `dt` matches its value `b` of an earlier pass as `a == b` or
`abs(a - b) <= rtol * max(abs(a), abs(b))`. A value of an array in `scaled` matches as
`a == b` or `abs(a - b) <= rtol * max(abs(a), abs(b), s)`, with `s` the largest finite
magnitude of that array in the two passes, so round-off on values near zero does not break
the match; a value below about 1e-5 of the array scale can differ from the anchor pass by
1e-4 of itself without breaking it. NaN never matches; Inf matches only itself.

For each period `p` in `1:max_period`, a pass counts towards the run of `p` when it matches
the earliest pass of the same phase in the run and the pass `p` before. Every pass of a
phase is thus within `rtol` of one fixed pass, so a drift cannot add up over the run. The
result is true when, for one period, the run reaches `max(window, 3p)` passes.

# Arguments

    - tr: [`PlasticCycleTracker`](@ref), reset for this loop
    - fields: tuple of arrays tested node by node
    - scaled: tuple of arrays tested against their field scale
    - Y: plastic residual of the pass (Pa)
    - dt: time step of the pass (s)
    - window: minimum run length in passes; 0 disables the test
    - rtol: relative tolerance per value

# Returns

    - true when the state is periodic, else false

# Throws

    - ArgumentError: the arrays do not hold the lengths set by [`reset_plastic_cycle!`](@ref),
      or `window` exceeds the window the ring was sized for
    - DomainError: `rtol` negative or not finite, or `window < 0`
"""
function plastic_cycle_record!(
    tr::PlasticCycleTracker,
    fields::Tuple,
    scaled::Tuple,
    Y::Real,
    dt::Real,
    window::Integer;
    rtol::Real=1.0e-9,
)
    (isfinite(rtol) && rtol >= 0) ||
        throw(DomainError(rtol, "rtol must be finite and >= 0"))
    window >= 0 || throw(DomainError(window, "window must be >= 0"))
    window > 0 || return false
    window <= tr.window ||
        throw(ArgumentError("window $window exceeds the window $(tr.window) of the reset"))
    _cycle_lengths_match(tr, fields, scaled) || throw(
        ArgumentError("plastic cycle state does not match the arrays set at the reset")
    )
    nslot = length(tr.ring)
    tr.npass += 1
    k = tr.npass
    ic = mod1(k, nslot)
    cur = tr.ring[ic]
    _pack_cycle_state!(cur, _pack_cycle_state!(cur, 0, fields), scaled)
    off = tr.nnode
    for q in axes(tr.scale, 1)
        n = tr.lens[length(fields) + q]
        tr.scale[q, ic] = _finite_absmax(cur, (off + 1):(off + n))
        off += n
    end
    tr.Y[ic] = Y
    tr.dt[ic] = dt
    cycled = false
    for p in eachindex(tr.run)
        k <= p && continue
        r = tr.run[p]
        # Earliest pass of the same phase in the run, kept in the ring.
        L = p * min((r + p) ÷ p, (nslot - 1) ÷ p)
        ok = _cycle_slot_match(tr, cur, ic, mod1(k - L, nslot), Y, dt, rtol)
        if ok && L > p
            ok = _cycle_slot_match(tr, cur, ic, mod1(k - p, nslot), Y, dt, rtol)
        end
        tr.run[p] = ok ? r + 1 : 0
        cycled |= tr.run[p] >= max(window, 3 * p)
    end
    return cycled
end # function plastic_cycle_record!

"""
    _cycle_slot_match(tr, cur, ic, j, Y, dt, rtol) -> Bool

True when the residual, the time step and the arrays of the current pass match those
stored in ring slot `j`.

# Arguments

    - tr: [`PlasticCycleTracker`](@ref)
    - cur: ring slot of the current pass
    - ic: index of the slot of the current pass
    - j: index of the slot of the earlier pass
    - Y: plastic residual of the current pass (Pa)
    - dt: time step of the current pass (s)
    - rtol: relative tolerance

# Returns

    - true when the two passes match
"""
function _cycle_slot_match(tr::PlasticCycleTracker, cur, ic, j, Y, dt, rtol)
    return _cycle_match(Y, tr.Y[j], rtol) &&
           _cycle_match(dt, tr.dt[j], rtol) &&
           _cycle_lag_match(tr, cur, tr.ring[j], ic, j, rtol)
end

"""
    _cycle_lengths_match(tr, fields, scaled) -> Bool

True when `fields` and `scaled` hold the number of arrays and the array lengths set by
[`reset_plastic_cycle!`](@ref).

# Arguments

    - tr: [`PlasticCycleTracker`](@ref)
    - fields, scaled: tuples of arrays of the pass

# Returns

    - true when the lengths match
"""
function _cycle_lengths_match(tr::PlasticCycleTracker, fields::Tuple, scaled::Tuple)
    length(fields) + length(scaled) == length(tr.lens) || return false
    length(scaled) == size(tr.scale, 1) || return false
    return _lengths_from(tr.lens, 1, fields) &&
           _lengths_from(tr.lens, length(fields) + 1, scaled)
end

"""
    _lengths_from(lens, i, t) -> Bool

True when the arrays of tuple `t` have the lengths `lens[i]`, `lens[i + 1]`, and so on.

# Arguments

    - lens: array lengths set at the reset
    - i: index in `lens` of the first array of `t`
    - t: tuple of arrays

# Returns

    - true when every length matches
"""
_lengths_from(lens, i, ::Tuple{}) = true
function _lengths_from(lens, i, t::Tuple)
    return length(first(t)) == lens[i] && _lengths_from(lens, i + 1, Base.tail(t))
end

"""
    AndersonHistory(depth, n)

History of the type II Anderson acceleration of the plastic loop (Walker and Ni, 2011):
ring buffers with the last `depth` differences of the packed state `x` and of the
fixed-point residual `f = G(x) - x`, the work vectors of the plastic loop (iterate `x`, map
value `g`, next iterate `xn`, least-squares weights `w`, set when `wset`), and the counters of
the stagnation stop.

# Arguments

    - depth: number of difference columns kept; 0 gives plain Picard updates
    - n: length of the packed state

# Returns

    - an empty history

# Throws

    - DomainError: `depth < 0`, `n < 0` or `depth > n`
"""
mutable struct AndersonHistory
    depth::Int
    dX::Matrix{Float64}
    dF::Matrix{Float64}
    A::Matrix{Float64}
    xprev::Vector{Float64}
    fprev::Vector{Float64}
    f::Vector{Float64}
    x::Vector{Float64}
    g::Vector{Float64}
    xn::Vector{Float64}
    w::Vector{Float64}
    wset::Bool
    ncol::Int
    head::Int
    hasprev::Bool
    best::Float64
    nstall::Int
    dt::Float64
end
function AndersonHistory(depth::Integer, n::Integer)
    depth >= 0 || throw(DomainError(depth, "depth must be at least 0"))
    n >= 0 || throw(DomainError(n, "state length must be at least 0"))
    depth <= n || throw(DomainError(depth, "depth must not exceed the state length $n"))
    m = depth == 0 ? 0 : n
    return AndersonHistory(
        depth,
        zeros(m, depth),
        zeros(m, depth),
        zeros(m, depth),
        zeros(m),
        zeros(m),
        zeros(m),
        zeros(m),
        zeros(m),
        zeros(m),
        zeros(m),
        false,
        0,
        1,
        false,
        Inf,
        0,
        NaN,
    )
end

"""
    reset_anderson!(h) -> h

Clear the difference columns and the stagnation counters of an [`AndersonHistory`](@ref).

# Arguments

    - h: history

# Returns

    - the history
"""
function reset_anderson!(h::AndersonHistory)
    _clear_anderson_columns!(h)
    h.wset = false
    h.best = Inf
    h.nstall = 0
    h.dt = NaN
    return h
end

"""
    _clear_anderson_columns!(h) -> h

Drop the difference columns and the previous iterate of an [`AndersonHistory`](@ref).

# Arguments

    - h: history

# Returns

    - the history
"""
function _clear_anderson_columns!(h::AndersonHistory)
    h.ncol = 0
    h.head = 1
    h.hasprev = false
    return h
end

"""
    anderson_mix!(xnext, x, G, h, w) -> Bool

Type II Anderson update of a fixed-point iteration (Walker and Ni, 2011, Section 3):
`xnext = G - (dX + dF) gamma` with `gamma` the minimiser of `|| w .* (f - dF gamma) ||`,
`f = G - x`, solved by column-pivoted QR. Columns with a diagonal of `R` below `1e-12`
of the largest are dropped. The update falls back to `xnext = G` without history, when
all columns vanish, and (clearing the history) when `f`, `gamma` or `xnext` is not finite
or `maximum(abs, gamma) >= 1e3`.

# Arguments

    - xnext: next iterate, overwritten
    - x: current iterate
    - G: map value `G(x)`
    - h: [`AndersonHistory`](@ref) of this iteration, updated with `x` and `f`
    - w: non-negative weight of each entry in the least-squares problem

# Returns

    - true when the Anderson update was applied, false for the fallback `xnext = G`

# Throws

    - DimensionMismatch: arrays of different lengths, or a history of another length
    - DomainError: a negative or non-finite weight
"""
function anderson_mix!(
    xnext::AbstractVector{Float64},
    x::AbstractVector{Float64},
    G::AbstractVector{Float64},
    h::AndersonHistory,
    w::AbstractVector{Float64},
)
    n = length(x)
    length(xnext) == n && length(G) == n && length(w) == n ||
        throw(DimensionMismatch("xnext, x, G and w must have the same length"))
    all(wi -> isfinite(wi) && wi >= 0.0, w) ||
        throw(DomainError(w, "weights must be finite and non-negative"))
    copyto!(xnext, G)
    h.depth == 0 && return false
    length(h.f) == n ||
        throw(DimensionMismatch("history length $(length(h.f)) differs from state $n"))
    f = h.f
    @. f = G - x
    if !all(isfinite, f)
        _clear_anderson_columns!(h)
        return false
    end
    if h.hasprev
        j = h.head
        @views @. h.dX[:, j] = x - h.xprev
        @views @. h.dF[:, j] = f - h.fprev
        h.head = j % h.depth + 1
        h.ncol = min(h.ncol + 1, h.depth)
    end
    copyto!(h.xprev, x)
    copyto!(h.fprev, f)
    h.hasprev = true
    return h.ncol > 0 && _anderson_apply!(xnext, G, h, w)
end

"""
    _anderson_apply!(xnext, G, h, w) -> Bool

Solve the weighted least-squares problem of [`anderson_mix!`](@ref) on the columns of `h`
and apply the update to `xnext`, which holds `G` on entry.

# Arguments

    - xnext: next iterate, holds `G` on entry
    - G: map value
    - h: history with at least one column and the residual in `h.f`
    - w: least-squares weights

# Returns

    - true when the update was applied, false for the fallback `xnext = G`
"""
function _anderson_apply!(xnext, G, h::AndersonHistory, w)
    m = h.ncol
    A = view(h.A, :, 1:m)
    @views @. A = w * h.dF[:, 1:m]
    F = qr!(A, ColumnNorm())
    r11 = abs(F.R[1, 1])
    r11 > 0.0 || return false
    k = count(j -> abs(F.R[j, j]) > 1.0e-12 * r11, 1:m)
    b = F.Q' * (w .* h.f)
    gamma = zeros(m)
    gamma[F.p[1:k]] .= UpperTriangular(F.R[1:k, 1:k]) \ b[1:k]
    if !all(isfinite, gamma) || maximum(abs, gamma) >= 1.0e3
        _clear_anderson_columns!(h)
        return false
    end
    for j in 1:m
        @views @. xnext -= gamma[j] * (h.dX[:, j] + h.dF[:, j])
    end
    if !all(isfinite, xnext)
        copyto!(xnext, G)
        _clear_anderson_columns!(h)
        return false
    end
    return true
end

"""
    anderson_stall_record!(h, Y, dt, npasses) -> Bool

Stagnation stop of the accelerated plastic loop: count the passes since the yield residual
`Y` last fell below `(1 - 1e-12)` times its smallest value. A change of the time step `dt`
restarts the count.

# Arguments

    - h: [`AndersonHistory`](@ref) of this loop
    - Y: RMS yield residual of the pass (Pa)
    - dt: time step of the pass (s)
    - npasses: number of passes without progress that stops the loop

# Returns

    - true when the loop must stop

# Throws

    - DomainError: `npasses < 1`
"""
function anderson_stall_record!(h::AndersonHistory, Y::Real, dt::Real, npasses::Integer)
    npasses >= 1 || throw(DomainError(npasses, "npasses must be at least 1"))
    if !isequal(dt, h.dt)
        h.dt = dt
        h.best = Inf
        h.nstall = 0
    end
    if Y < h.best * (1.0 - 1.0e-12)
        h.best = Y
        h.nstall = 0
        return false
    end
    h.nstall += 1
    return h.nstall >= npasses
end

"""
    anderson_pack_viscosity!(s, ETA, ETA0) -> s

Pack nodal viscosities as `s = log(min(1, ETA / ETA0))`: exactly 0 on a node that does
not yield, negative on a yielding node.

# Arguments

    - s: packed values, overwritten
    - ETA: viscosity (Pa s)
    - ETA0: viscosity without yielding (Pa s)

# Returns

    - s

# Throws

    - DimensionMismatch: arrays of different lengths
"""
function anderson_pack_viscosity!(s::AbstractVector{Float64}, ETA, ETA0)
    length(s) == length(ETA) == length(ETA0) ||
        throw(DimensionMismatch("s, ETA and ETA0 must have the same length"))
    @inbounds for k in 1:length(s)
        s[k] = ETA[k] < ETA0[k] ? log(ETA[k] / ETA0[k]) : 0.0
    end
    return s
end

"""
    anderson_unpack_viscosity!(ETA, YNY, YNY_inv_ETA, s, ETA0, etamin, etamax) -> ETA

Inverse of [`anderson_pack_viscosity!`](@ref) with the yield rule of the plastic loop: a node
with `s >= -1e-12` does not yield and gets `ETA0`; any other node yields, with
`ETA = clamp(ETA0 exp(s), etamin, min(etamax, ETA0))`.

# Arguments

    - ETA: viscosity (Pa s), overwritten
    - YNY: yield flag, overwritten
    - YNY_inv_ETA: `YNY / ETA` (1 / (Pa s)), overwritten
    - s: packed values
    - ETA0: viscosity without yielding (Pa s)
    - etamin: lower viscosity limit (Pa s)
    - etamax: upper viscosity limit (Pa s)

# Returns

    - ETA

# Throws

    - DimensionMismatch: arrays of different lengths
"""
function anderson_unpack_viscosity!(ETA, YNY, YNY_inv_ETA, s, ETA0, etamin, etamax)
    n = length(s)
    length(ETA) == n && length(YNY) == n && length(YNY_inv_ETA) == n && length(ETA0) == n ||
        throw(DimensionMismatch("viscosity arrays must have the length of s"))
    @inbounds for k in 1:n
        if s[k] < -1.0e-12
            hi = min(etamax, ETA0[k])
            e = clamp(ETA0[k] * exp(s[k]), min(etamin, hi), hi)
            ETA[k] = e
            YNY[k] = true
            YNY_inv_ETA[k] = 1.0 / e
        else
            ETA[k] = ETA0[k]
            YNY[k] = false
            YNY_inv_ETA[k] = 0.0
        end
    end
    return ETA
end

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
    ETA,
    ETA5,
    ETA00,
    YNY,
    YNY5,
    YNY00,
    YNY_inv_ETA,
    dt,
    iplast;
    dtstep::Int=200,
    dtcoefdn::Real=0.5,
)
    ETA .= ETA5
    YNY .= YNY5
    @views @. YNY_inv_ETA = YNY / ETA
    return dt
end # function finalize_plastic_iteration_pass

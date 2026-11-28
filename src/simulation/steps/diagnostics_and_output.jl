# Diagnostics, streaming telemetry, checkpoint serialization, and progress reporting

"""
    telemetry_porosity(markers::MarkerArrays, PHI::AbstractMatrix) -> (phi_max, phi_mean)

Porosity of the `phi_max` and `phi_mean` columns of `telemetry.csv`. With pore water
(`porosity.model = :none`) these are the maximum and the node mean of the water porosity
`PHI` at the P nodes. Without it (`:henke2012`) they are the maximum and the marker mean
of the dry porosity `phi_d` over the rock markers (`tm < 3`), since `PHI` stays at
`poroelasticity.phimin` there.

# Arguments
- `markers`: marker arrays.
- `PHI`: water porosity at the P nodes [-].

# Returns
- `(phi_max, phi_mean)` [-]; `(0.0, 0.0)` for `:henke2012` without rock markers.
"""
function telemetry_porosity(markers::MarkerArrays, PHI::AbstractMatrix)
    has_pore_water(markers) && return (maximum(PHI), sum(PHI) / length(PHI))
    phi_d = markers.groups.porosity.phi_d
    tm = markers.core.tm
    phi_max = 0.0
    phi_sum = 0.0
    n_rock = 0
    @inbounds for m in eachindex(tm)
        tm[m] < 3 || continue
        phi_max = max(phi_max, phi_d[m])
        phi_sum += phi_d[m]
        n_rock += 1
    end
    return n_rock == 0 ? (0.0, 0.0) : (phi_max, phi_sum / n_rock)
end

# Type-7 quantile (linear between order statistics) of a sorted, non-empty vector.
function _sorted_quantile(v::AbstractVector{Float64}, p::Float64)
    h = 1.0 + p * (length(v) - 1)
    lo = floor(Int, h)
    hi = min(lo + 1, length(v))
    return v[lo] + (h - lo) * (v[hi] - v[lo])
end

"""
    porous_shell_depth(xm, ym, tm, pg; xcenter, ycenter, rplanet, dr, nsector=16,
                       measure=:unsintered) -> (; median, q25, q75, nsector_used)

Depth [m] of the outer porous shell of a planetesimal, measured on the rock markers (`tm`
1 or 2) of the dry porosity group `pg` of `porosity.model = :henke2012`.

The slice is split into `nsector` equal angular sectors about `(xcenter, ycenter)`, and
each sector into radial bins of width `dr`, `(rplanet - k dr, rplanet - (k - 1) dr]`; rock
markers above `rplanet` are not counted (the pressure there is 0, so they never sinter).
From the surface inward, the first non-empty bin
whose median marker is sintered ends the shell, and the depth of the sector is the outer
edge of that bin, `(k - 1) dr`. A sector without a sintered bin has the depth `rplanet`.
A sector with few markers per bin skips empty outer bins and can give a deeper shell, so
keep `nsector` well below the number of rock markers per ring of width `dr`.

A marker is sintered for `measure = :unsintered` when `dphi_sinter >= 1e-3`, the porosity
decrease by hot sintering (cold pressing does not count). For `measure = :phi_gt_04` it is
sintered when `phi_d <= 0.4`, so that the depth is that of the region with `phi_d > 0.4`;
this measure is meaningful only when the cold-pressed porosity at the surface exceeds 0.4.

# Arguments
- `xm`, `ym`: marker coordinates [m].
- `tm`: marker material (1 mantle, 2 crust, 3 air).
- `pg`: dry porosity group (`phi_d`, `dphi_sinter` are read), or `nothing`.
- `xcenter`, `ycenter`: centre of the planetesimal [m].
- `rplanet`: planetesimal radius [m].
- `dr`: width of the radial bins [m].
- `nsector`: number of angular sectors.
- `measure`: `:unsintered` or `:phi_gt_04`.

# Returns
- `(; median, q25, q75, nsector_used)`: median and type-7 quartiles of the sector depths
  [m] over the `nsector_used` sectors that hold rock markers; NaN depths and 0 sectors
  without rock markers.

# Throws
- `ArgumentError`: `pg` is `nothing` (`porosity.model = :none`), `nsector < 1`, an unknown
  `measure`, or marker arrays of different lengths.
- `DomainError`: `rplanet` or `dr` not positive and finite, also after the conversion to
  `Float64`; a centre not finite; or, for a rock marker, a non-finite `xm` or `ym` or a
  non-finite value of the measure in use (`dphi_sinter` or `phi_d`).
"""
function porous_shell_depth(
    xm::AbstractVector{<:Real},
    ym::AbstractVector{<:Real},
    tm::AbstractVector{<:Integer},
    pg::Union{PorosityGroup,Nothing};
    xcenter::Real,
    ycenter::Real,
    rplanet::Real,
    dr::Real,
    nsector::Integer=16,
    measure::Symbol=:unsintered,
)
    pg === nothing &&
        throw(ArgumentError("porous_shell_depth needs the dry porosity group of henke2012"))
    nsector >= 1 || throw(ArgumentError("nsector must be >= 1, got $nsector"))
    measure in (:unsintered, :phi_gt_04) ||
        throw(ArgumentError("measure must be :unsintered or :phi_gt_04, got $measure"))
    n = length(xm)
    lengths = (length(ym), length(tm), length(pg.phi_d), length(pg.dphi_sinter))
    all(==(n), lengths) || throw(ArgumentError("marker arrays differ in length"))
    R = Float64(require_positive_finite(rplanet, "rplanet"))
    h = Float64(require_positive_finite(dr, "dr"))
    # A value beyond the Float64 range is Inf or 0 after the conversion
    require_positive_finite(R, "rplanet as a Float64")
    require_positive_finite(h, "dr as a Float64")
    isfinite(xcenter) || throw(DomainError(xcenter, "xcenter must be finite"))
    isfinite(ycenter) || throw(DomainError(ycenter, "ycenter must be finite"))
    nbin = max(1, ceil(Int, R / h))
    vals = [Float64[] for _ in 1:nbin, _ in 1:nsector]
    for m in 1:n
        (tm[m] == 1 || tm[m] == 2) || continue
        isfinite(xm[m]) || throw(DomainError(xm[m], "xm of rock marker $m must be finite"))
        isfinite(ym[m]) || throw(DomainError(ym[m], "ym of rock marker $m must be finite"))
        v = measure === :unsintered ? pg.dphi_sinter[m] : pg.phi_d[m]
        isfinite(v) || throw(
            DomainError(
                v,
                string(
                    measure === :unsintered ? "dphi_sinter" : "phi_d",
                    " of rock marker $m must be finite",
                ),
            ),
        )
        dx = xm[m] - xcenter
        dy = ym[m] - ycenter
        r = hypot(dx, dy)
        # Above the surface P = 0, so a rock marker there never sinters: not counted
        r > R && continue
        s = min(floor(Int, (atan(dy, dx) + pi) / (2pi) * nsector) + 1, nsector)
        k = min(floor(Int, (R - r) / h) + 1, nbin)
        push!(vals[k, s], v)
    end
    depths = Float64[]
    for s in 1:nsector
        all(isempty, view(vals, :, s)) && continue
        d = R
        for k in 1:nbin
            v = vals[k, s]
            isempty(v) && continue
            med = _sorted_quantile(sort!(v), 0.5)
            if (measure === :unsintered ? med >= 1.0e-3 : med <= 0.4)
                d = (k - 1) * h
                break
            end
        end
        push!(depths, d)
    end
    isempty(depths) && return (; median=NaN, q25=NaN, q75=NaN, nsector_used=0)
    sort!(depths)
    return (;
        median=_sorted_quantile(depths, 0.5),
        q25=_sorted_quantile(depths, 0.25),
        q75=_sorted_quantile(depths, 0.75),
        nsector_used=length(depths),
    )
end

"""
    porous_shell_depth(st::SimulationState, coords::GridCoordinates, geometry::GeometryConfig;
                       nsector=16, measure=:unsintered) -> (; median, q25, q75, nsector_used)

Porous shell depth of a simulation state: the marker method with the planet centre and
radius of `geometry` and radial bins of half the smaller grid spacing.

# Throws
- `ArgumentError`: the state has no dry porosity group (`porosity.model = :none`).
- `DomainError`: from the marker method, for example a rock marker with a non-finite
  position or value of the measure.
"""
function porous_shell_depth(
    st::SimulationState,
    coords::GridCoordinates,
    geometry::GeometryConfig;
    nsector::Integer=16,
    measure::Symbol=:unsintered,
)
    groups = st.markers.groups
    pg = haskey(groups, :porosity) ? groups.porosity : nothing
    core = st.markers.core
    return porous_shell_depth(
        core.xm,
        core.ym,
        core.tm,
        pg;
        xcenter=geometry.xcenter,
        ycenter=geometry.ycenter,
        rplanet=geometry.rplanet,
        dr=0.5 * min(coords.dx, coords.dy),
        nsector=nsector,
        measure=measure,
    )
end

"""
Evaluate diagnostics, write telemetry rows, and serialize simulation checkpoints.

$(SIGNATURES)
"""
function advance_step_diagnostics!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig,
    ws::SimulationWorkspaces;
    output_path::String,
    telemetry_io=nothing,
    dt_aphimax_step_max::Float64=0.0,
    n_flips_last::Int=0,
    n_flips_step_total::Int=0,
    timestep_begin=nothing,
    progress_bar=nothing,
)
    timestep = state.timestep
    n_steps_val = cfg.time.n_steps
    savematstep_val = cfg.output.savematstep

    should_save_snapshot = if cfg.output.mode in (:snapshots, :both)
        savematstep_val > 0 && (timestep % savematstep_val == 0)
    elseif cfg.output.mode == :telemetry
        cfg.output.save_final && (timestep == n_steps_val)
    else
        false
    end

    need_telemetry = (
        telemetry_io !== nothing && (
            (cfg.output.telemetrystep > 0 && timestep % cfg.output.telemetrystep == 0) ||
            timestep == n_steps_val
        )
    )

    core_diag_active =
        cfg.coreformation.percolation_active ||
        cfg.coreformation.settling_active ||
        cfg.metal_partition.active
    if (core_diag_active && state.markers.Xfe_bulk !== nothing) &&
        (need_telemetry || should_save_snapshot || cfg.hydrothermal.active)
        core_budgets = compute_core_volatile_budgets(
            state.markers.xm,
            state.markers.ym,
            state.markers.tm,
            state.markers.Xfe_bulk,
            state.markers.Xfe_H_m,
            state.markers.Xfe_C_m,
            state.markers.Xfe_N_m,
            state.markers.Xfe_S_m,
            length(state.markers);
            coords=coords,
            w3d_m=state.markers.w3d_m,
            xcenter=state.accumulators.xcenter,
            ycenter=state.accumulators.ycenter,
            rplanet=state.accumulators.rplanet,
            rho_metal=cfg.coreformation.rho_metal,
            core_radius_fraction=cfg.metal_partition.core_radius_fraction,
            phi_core_threshold=cfg.metal_partition.phi_core_threshold,
        )
        if cfg.metal_partition.active
            state.accumulators.core_budgets = core_budgets
        end
        state.accumulators.rcore =
            if core_budgets !== nothing && core_budgets.V_core_total > 0.0
                cbrt(3.0 * core_budgets.V_core_total / (4.0 * π))
            else
                0.0
            end
    end

    if need_telemetry
        core_radius_current = state.accumulators.rcore
        tk2 = state.grids.tk2
        phi_max, phi_mean = telemetry_porosity(state.markers, state.grids.PHI)
        Fm = state.markers.Fm
        M_vent_all = (
            state.accumulators.M_vent_H2O_total +
            state.accumulators.M_vent_C_total +
            state.accumulators.M_vent_N_total +
            state.accumulators.M_vent_S_total
        )
        stream_telemetry_row!(
            telemetry_io,
            timestep,
            s_to_Ma(state.timesum; yearlength=cfg.time.yearlength),
            state.dt / cfg.time.yearlength,
            state.accumulators.rplanet,
            core_radius_current,
            maximum(tk2),
            sum(tk2) / length(tk2),
            phi_max,
            phi_mean,
            M_vent_all,
            state.accumulators.M_vent_H2O_total,
            state.accumulators.M_atm_total,
            state.accumulators.M_escaped_total,
            Fm !== nothing ? maximum(Fm) : 0.0,
            Fm !== nothing ? (sum(Fm) / length(Fm)) : 0.0,
            dt_aphimax_step_max,
            n_flips_last,
            n_flips_step_total,
        )
    end

    if should_save_snapshot
        if cfg.phase_tracking.active &&
            cfg.phase_tracking.track_regional_modes &&
            state.markers.Xfe_bulk !== nothing
            state.accumulators.regional_mineral_modes = compute_regional_mineral_modes(
                state.markers.xm,
                state.markers.ym,
                state.markers.tm,
                state.markers.tkm,
                state.markers.Xfe_bulk,
                state.markers.Xfe_S_m,
                state.markers.Xfe_C_m,
                state.markers.Xfe_N_m,
                length(state.markers);
                coords=coords,
                w3d_m=state.markers.w3d_m,
                cfg=cfg.phase_tracking,
                xcenter=state.accumulators.xcenter,
                ycenter=state.accumulators.ycenter,
                rplanet=state.accumulators.rplanet,
                rho_metal=cfg.coreformation.rho_metal,
            )
        end
        save_state(output_path, state, coords, cfg)
    end

    maxT = maximum(state.grids.tk2)
    marknum = length(state.markers)
    if timestep_begin !== nothing
        duration_str = Dates.canonicalize(
            Dates.CompoundPeriod(Dates.now() - timestep_begin)
        )
        @info "timestep $timestep computed in $duration_str"
    else
        @info "timestep $timestep computed"
    end
    @info "total time = $(s_to_Ma(state.timesum; yearlength=cfg.time.yearlength)) Ma"
    @info "markers in use = $marknum"
    @info "max T = $maxT K"

    if progress_bar !== nothing
        endtime_val = cfg.time.endtime * cfg.time.yearlength
        showvals =
            () -> [
                (:timestep, timestep),
                (:markers, marknum),
                (:maxT_K, maxT),
                (:dt_s, state.dt),
                (:timesum_Ma, s_to_Ma(state.timesum; yearlength=cfg.time.yearlength)),
                (
                    :to_go_Ma,
                    s_to_Ma(endtime_val - state.timesum; yearlength=cfg.time.yearlength),
                ),
            ]
        next!(progress_bar; showvalues=showvals)
    end

    return nothing
end

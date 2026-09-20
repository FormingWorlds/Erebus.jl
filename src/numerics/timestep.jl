"""
Compute velocity-/displacement-limited time step.

$(SIGNATURES)

# Details

    - vx: solid vx velocity at Vx nodes
    - vy: solid vy velocity at Vy nodes
    - vxf: fluid vx velocity at Vx nodes
    - vyf: fluid vy velocity at Vy nodes
    - dt: current time step
    - aphimax: maximum observed porosity coefficient
   
# Returns

    - dt: displacement time step
"""
function compute_displacement_timestep(
    vx,
    vy,
    vxf,
    vyf,
    dt,
    aphimax;
    coords=nothing,
    dx_val=coords === nothing ? dx : coords.dx,
    dy_val=coords === nothing ? dy : coords.dy,
    dxymax_val::Real=0.05,
    dphimax_val::Real=100.01,
)
    maxvx = maximum(abs, vx)
    maxvy = maximum(abs, vy)
    maxvxf = maximum(abs, vxf)
    maxvyf = maximum(abs, vyf)
    @info "dt before velocity limitations = $dt s"
    dt = ifelse(dt*maxvx > dxymax_val*dx_val, dxymax_val*dx_val*inv(maxvx), dt)
    @info "dt after vx limitation = $dt s"
    dt = ifelse(dt*maxvy > dxymax_val*dy_val, dxymax_val*dy_val*inv(maxvy), dt)
    @info "dt after vy limitation = $dt s"
    dt = ifelse(dt*maxvxf > dxymax_val*dx_val, dxymax_val*dx_val*inv(maxvxf), dt)
    @info "dt after vxf limitation = $dt s"
    dt = ifelse(dt*maxvyf > dxymax_val*dy_val, dxymax_val*dy_val*inv(maxvyf), dt)
    @info "dt after vyf limitation = $dt s"
    dt = ifelse(dt*aphimax > dphimax_val, dphimax_val*inv(aphimax), dt)
    @info "dt after aphimax limitation = $dt s"
    return dt
end # function compute_displacement_timestep

"""
    compute_adaptive_timestep(
        vx, vy, vxf, vyf, dt, aphimax;
        coords=nothing,
        dx_val=coords === nothing ? dx : coords.dx,
        dy_val=coords === nothing ? dy : coords.dy,
        dxymax_val::Real=0.05,
        dphimax_val::Real=100.01,
        maxDTcurrent=0.0,
        DTmax_val::Real=20.0,
        dt_longest_val::Real=1.0e11,
        dt_min=1.0,
    )

Compute multi-criterion adaptive timestep constrained by velocity CFL, porosity compaction,
thermal variation, and stability bounds.
"""
function compute_adaptive_timestep(
    vx,
    vy,
    vxf,
    vyf,
    dt,
    aphimax;
    coords=nothing,
    dx_val=coords === nothing ? dx : coords.dx,
    dy_val=coords === nothing ? dy : coords.dy,
    dxymax_val::Real=0.05,
    dphimax_val::Real=100.01,
    dt_ref=nothing,
    maxDTcurrent=0.0,
    DTmax_val::Real=20.0,
    dt_longest_val::Real=1.0e11,
    dt_min=1.0,
    max_v_seg::Real=0.0,
    max_subcycles::Integer=2000,
    cfl_settling::Real=0.5,
    DQPF::Union{AbstractMatrix{<:Real},Nothing}=nothing,
    cfl_reaction::Real=0.5,
    dphi_reaction_max::Real=0.01,
)
    dt_cand = compute_displacement_timestep(
        vx,
        vy,
        vxf,
        vyf,
        dt,
        aphimax;
        coords=coords,
        dx_val=dx_val,
        dy_val=dy_val,
        dxymax_val=dxymax_val,
        dphimax_val=dphimax_val,
    )
    ref_dt = dt_ref === nothing ? dt : dt_ref
    if maxDTcurrent > DTmax_val && maxDTcurrent > 0.0
        dt_cand = min(dt_cand, ref_dt * (DTmax_val * inv(maxDTcurrent)))
    end
    if max_v_seg > 0.0
        min_dx = min(dx_val, dy_val)
        dt_cfl_seg = cfl_settling * min_dx / max_v_seg
        dt_seg_bound = max_subcycles * dt_cfl_seg
        dt_cand = min(dt_cand, dt_seg_bound)
    end
    if DQPF !== nothing
        max_dqpf = maximum(abs, DQPF)
        if max_dqpf > 0.0
            dt_rxn = Float64(cfl_reaction) * Float64(dphi_reaction_max) / max_dqpf
            dt_cand = min(dt_cand, dt_rxn)
        end
    end
    dt_cand = clamp(dt_cand, dt_min, dt_longest_val)
    return dt_cand
end

# Hydromechanical Stokes-Darcy assembly, solve, and plastic iteration step

"""
Assemble and solve the coupled Stokes-Darcy linear system of equations.

$(SIGNATURES)
"""
function assemble_and_solve_hydromechanical!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig,
    ws::SimulationWorkspaces;
    titer::Int,
    iplast::Int,
    cur_betasolid::Float64,
    cur_betafluid::Float64,
)
    recompute_bulk_viscosity!(
        state.grids.ETA,
        state.grids.ETAP,
        state.grids.ETAPHI,
        state.grids.PHI,
        cfg.solver.etaphikoef,
    )
    fill!(state.grids.S_vent_grid, 0.0)

    darcy_elim_val = cfg.solver.darcy_elimination
    if darcy_elim_val
        L = assemble_hydromechanical_4var_lse!(
            state.grids.ETA,
            state.grids.ETAP,
            state.grids.GGG,
            state.grids.GGGP,
            state.grids.SXY0,
            state.grids.SXX0,
            state.grids.RHOX,
            state.grids.RHOY,
            state.grids.RHOFX,
            state.grids.RHOFY,
            state.grids.RX,
            state.grids.RY,
            state.grids.ETAPHI,
            state.grids.BETAPHI,
            state.grids.PHI,
            state.grids.gx,
            state.grids.gy,
            state.grids.pr0,
            state.grids.pf0,
            state.grids.DMP,
            state.dt,
            ws.R;
            coords=coords,
            betasolid=cur_betasolid,
            betafluid=cur_betafluid,
            phimin=cfg.poroelasticity.phimin,
            phimax=cfg.poroelasticity.phimax,
            hydrofracture=cfg.poroelasticity.hydrofracture,
            pr=state.grids.pr,
            pf=state.grids.pf,
            TEN=state.grids.TEN,
            KX=state.grids.KX,
            KY=state.grids.KY,
            PHIX=state.grids.PHIX,
            PHIY=state.grids.PHIY,
            kphim0=cfg.materials.kphim0,
            phim0_val=cfg.thermodynamics.phim0,
            kappa_frac=cfg.poroelasticity.kappa_frac,
            gamma_frac=cfg.poroelasticity.gamma_frac,
            k_frac_max=cfg.poroelasticity.k_frac_max,
            ramp_width=cfg.poroelasticity.ramp_width,
            theta_frac=cfg.poroelasticity.theta_frac,
            rx_floor_prefactor=cfg.poroelasticity.rx_floor_prefactor,
            rx_eff_prev=ws.hydromech.rx_eff_prev,
            ry_eff_prev=ws.hydromech.ry_eff_prev,
            rx_eff_out=ws.hydromech.rx_eff,
            ry_eff_out=ws.hydromech.ry_eff,
            L=ws.hydromech.L,
            venting=cfg.venting.active,
            venting_mode=cfg.venting.mode,
            k_vent=cfg.venting.k_vent,
            conductance_factor=cfg.venting.conductance_factor,
            ice_sealing=cfg.venting.ice_sealing,
            t_freeze=cfg.venting.t_freeze,
            dt_seal=cfg.venting.dt_seal,
            k_seal_min_ratio=cfg.venting.k_seal_min_ratio,
            rplanet=state.accumulators.rplanet,
            xcenter=state.accumulators.xcenter,
            ycenter=state.accumulators.ycenter,
            P_amb=state.accumulators.P_amb,
            venting_species=cfg.venting.species,
            tk=state.grids.tk1,
            eta_fluid_surf=cfg.materials.etafluidmm[2],
            L_sub=cfg.venting.L_sublimation,
            S_vent_out=state.grids.S_vent_grid,
            DQPF=state.grids.DQPF,
            fluid_overpressure_coupling=cfg.reaction.fluid_overpressure_coupling,
            psurface=cfg.geometry.psurface,
            workspace=ws.hydromech,
        )
    else
        L = assemble_hydromechanical_lse!(
            state.grids.ETA,
            state.grids.ETAP,
            state.grids.GGG,
            state.grids.GGGP,
            state.grids.SXY0,
            state.grids.SXX0,
            state.grids.RHOX,
            state.grids.RHOY,
            state.grids.RHOFX,
            state.grids.RHOFY,
            state.grids.RX,
            state.grids.RY,
            state.grids.ETAPHI,
            state.grids.BETAPHI,
            state.grids.PHI,
            state.grids.gx,
            state.grids.gy,
            state.grids.pr0,
            state.grids.pf0,
            state.grids.DMP,
            state.dt,
            ws.R;
            coords=coords,
            betasolid=cur_betasolid,
            betafluid=cur_betafluid,
            phimin=cfg.poroelasticity.phimin,
            phimax=cfg.poroelasticity.phimax,
            hydrofracture=cfg.poroelasticity.hydrofracture,
            pr=state.grids.pr,
            pf=state.grids.pf,
            TEN=state.grids.TEN,
            KX=state.grids.KX,
            KY=state.grids.KY,
            PHIX=state.grids.PHIX,
            PHIY=state.grids.PHIY,
            kphim0=cfg.materials.kphim0,
            phim0_val=cfg.thermodynamics.phim0,
            kappa_frac=cfg.poroelasticity.kappa_frac,
            gamma_frac=cfg.poroelasticity.gamma_frac,
            k_frac_max=cfg.poroelasticity.k_frac_max,
            ramp_width=cfg.poroelasticity.ramp_width,
            theta_frac=cfg.poroelasticity.theta_frac,
            rx_floor_prefactor=cfg.poroelasticity.rx_floor_prefactor,
            rx_eff_prev=ws.hydromech.rx_eff_prev,
            ry_eff_prev=ws.hydromech.ry_eff_prev,
            rx_eff_out=ws.hydromech.rx_eff,
            ry_eff_out=ws.hydromech.ry_eff,
            L=ws.hydromech.L,
            venting=cfg.venting.active,
            venting_mode=cfg.venting.mode,
            k_vent=cfg.venting.k_vent,
            conductance_factor=cfg.venting.conductance_factor,
            ice_sealing=cfg.venting.ice_sealing,
            t_freeze=cfg.venting.t_freeze,
            dt_seal=cfg.venting.dt_seal,
            k_seal_min_ratio=cfg.venting.k_seal_min_ratio,
            rplanet=state.accumulators.rplanet,
            xcenter=state.accumulators.xcenter,
            ycenter=state.accumulators.ycenter,
            P_amb=state.accumulators.P_amb,
            venting_species=cfg.venting.species,
            tk=state.grids.tk1,
            eta_fluid_surf=cfg.materials.etafluidmm[2],
            L_sub=cfg.venting.L_sublimation,
            S_vent_out=state.grids.S_vent_grid,
            DQPF=state.grids.DQPF,
            fluid_overpressure_coupling=cfg.reaction.fluid_overpressure_coupling,
            psurface=cfg.geometry.psurface,
            workspace=ws.hydromech,
        )
    end

    if darcy_elim_val
        ws.hydromech.pr_presolve .= state.grids.pr
        ws.hydromech.pf_presolve .= state.grids.pf
    end

    if cfg.solver.hydromech_solver == :matrix_free
        rx_mf = if cfg.poroelasticity.hydrofracture && ws.hydromech.rx_eff !== nothing
            ws.hydromech.rx_eff
        else
            state.grids.RX
        end
        ry_mf = if cfg.poroelasticity.hydrofracture && ws.hydromech.ry_eff !== nothing
            ws.hydromech.ry_eff
        else
            state.grids.RY
        end
        op_mf = MatrixFreeStokesDarcyOperator(
            state.grids.ETA,
            state.grids.ETAP,
            state.grids.GGG,
            state.grids.GGGP,
            state.grids.RHOX,
            state.grids.RHOY,
            state.grids.RHOFX,
            state.grids.RHOFY,
            rx_mf,
            ry_mf,
            state.grids.ETAPHI,
            state.grids.BETAPHI,
            state.grids.PHI,
            state.grids.gx,
            state.grids.gy,
            state.dt;
            coords=coords,
            betasolid=cur_betasolid,
            betafluid=cur_betafluid,
            phimin=cfg.poroelasticity.phimin,
            phimax=cfg.poroelasticity.phimax,
        )
        _, stats = solve_hydromechanical_iterative!(
            op_mf,
            ws.R,
            ws.S;
            coords=coords,
            method=cfg.solver.krylov_method,
            rtol=cfg.solver.krylov_rtol,
            atol=cfg.solver.krylov_atol,
            maxiter=cfg.solver.krylov_maxiter,
            restart=cfg.solver.krylov_restart,
            preconditioner=cfg.solver.preconditioner,
            mg_levels=cfg.solver.mg_levels,
            mg_pre_smooth=cfg.solver.mg_pre_smooth,
            mg_post_smooth=cfg.solver.mg_post_smooth,
            mg_smoother=cfg.solver.mg_smoother,
            mg_omega=cfg.solver.mg_omega,
        )
        if !stats.solved
            error("Matrix-free Krylov solver failed: $(stats.status)")
        end
    elseif cfg.solver.hydromech_solver == :iterative
        _, stats = solve_hydromechanical_iterative!(
            L,
            ws.R,
            ws.S;
            coords=coords,
            method=cfg.solver.krylov_method,
            rtol=cfg.solver.krylov_rtol,
            atol=cfg.solver.krylov_atol,
            maxiter=cfg.solver.krylov_maxiter,
            restart=cfg.solver.krylov_restart,
            preconditioner=cfg.solver.preconditioner,
            mg_levels=cfg.solver.mg_levels,
            mg_pre_smooth=cfg.solver.mg_pre_smooth,
            mg_post_smooth=cfg.solver.mg_post_smooth,
            mg_smoother=cfg.solver.mg_smoother,
            mg_omega=cfg.solver.mg_omega,
        )
        if !stats.solved
            error("Iterative Krylov solver failed: $(stats.status)")
        end
    else
        if ws.hydromech_cache === nothing
            hydromech_prob = LinearProblem(L, ws.R)
            ws.hydromech_cache = init(
                hydromech_prob, UMFPACKFactorization(; reuse_symbolic=true)
            )
        else
            ws.hydromech_cache.A = L
            ws.hydromech_cache.b = ws.R
        end
        hydromech_sol = solve!(ws.hydromech_cache)
        if !LinearSolve.SciMLBase.successful_retcode(hydromech_sol) ||
            !all(isfinite, hydromech_sol.u)
            error("LinearSolve failed with retcode $(hydromech_sol.retcode)")
        end
        ws.S .= hydromech_sol.u
    end

    return nothing
end

"""
Postprocess hydromechanical solution, adapt timestep, and check plastic convergence.

$(SIGNATURES)
"""
function postprocess_hydromechanical_solution!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig,
    ws::SimulationWorkspaces;
    titer::Int,
    iplast::Int,
    dt_step_initial::Float64,
    cur_betasolid::Float64,
    cur_betafluid::Float64,
)
    g = state.grids
    darcy_elim_val = cfg.solver.darcy_elimination
    if darcy_elim_val
        process_hydromechanical_4var_solution!(ws.S, g.vx, g.vy, g.pr, g.pf; coords=coords)
        reconstruct_darcy_fluxes!(
            g.qxD,
            g.qyD,
            g.pf,
            g.RHOFX,
            g.RHOFY,
            g.RX,
            g.RY,
            g.gx,
            g.gy,
            coords;
            hydrofracture=cfg.poroelasticity.hydrofracture,
            pr=ws.hydromech.pr_presolve,
            pf_eff=ws.hydromech.pf_presolve,
            TEN=g.TEN,
            KX=g.KX,
            KY=g.KY,
            PHIX=g.PHIX,
            PHIY=g.PHIY,
            PHI=g.PHI,
            kphim0=cfg.materials.kphim0,
            phim0_val=cfg.thermodynamics.phim0,
            kappa_frac=cfg.poroelasticity.kappa_frac,
            gamma_frac=cfg.poroelasticity.gamma_frac,
            k_frac_max=cfg.poroelasticity.k_frac_max,
            ramp_width=cfg.poroelasticity.ramp_width,
            rx_floor_prefactor=cfg.poroelasticity.rx_floor_prefactor,
            rx_eff=ws.hydromech.rx_eff,
            ry_eff=ws.hydromech.ry_eff,
        )
    else
        process_hydromechanical_solution!(
            ws.S, g.vx, g.vy, g.pr, g.qxD, g.qyD, g.pf; coords=coords
        )
    end

    n_flips_iter = 0
    for j in 1:coords.Nx, i in 1:coords.Ny
        peff_c =
            0.25 * (
                g.pr[i, j] + g.pr[i + 1, j] + g.pr[i, j + 1] + g.pr[i + 1, j + 1] -
                g.pf[i, j] - g.pf[i + 1, j] - g.pf[i, j + 1] - g.pf[i + 1, j + 1]
            )
        is_breached = is_hydrofracture_breached(peff_c, g.TEN[i, j])
        if is_breached != ws.fractured_cells_prev[i, j]
            n_flips_iter += 1
        end
        ws.fractured_cells[i, j] = is_breached
    end
    ws.fractured_cells_prev .= ws.fractured_cells
    ws.hydromech.rx_eff_prev .= ws.hydromech.rx_eff
    ws.hydromech.ry_eff_prev .= ws.hydromech.ry_eff

    aphimax = compute_Aϕ!(
        g.APHI,
        g.ETAPHI,
        g.BETAPHI,
        g.PHI,
        g.pr,
        g.pf,
        g.pr0,
        g.pf0,
        state.dt;
        coords=coords,
        betasolid=cur_betasolid,
        phimin=cfg.poroelasticity.phimin,
        phimax=cfg.poroelasticity.phimax,
        S_vent=cfg.venting.active ? g.S_vent_grid : nothing,
    )

    compute_fluid_velocities!(
        g.PHIX, g.PHIY, g.qxD, g.qyD, g.vx, g.vy, g.vxf, g.vyf; coords=coords
    )

    maxDTcurrent = maximum(abs, g.DT0)
    state.dt = compute_adaptive_timestep(
        g.vx,
        g.vy,
        g.vxf,
        g.vyf,
        state.dt,
        aphimax;
        coords=coords,
        dxymax_val=cfg.time.dxymax,
        dphimax_val=cfg.solver.dphimax,
        dt_ref=dt_step_initial,
        maxDTcurrent=maxDTcurrent,
        DTmax_val=cfg.time.DTmax,
        dt_longest_val=cfg.time.dt_longest * cfg.time.yearlength,
        max_v_seg=state.accumulators.max_v_seg_prev,
        max_subcycles=cfg.coreformation.max_subcycles,
        cfl_settling=cfg.coreformation.cfl_settling,
        DQPF=g.DQPF,
        cfl_reaction=cfg.reaction.cfl_reaction,
        dphi_reaction_max=cfg.reaction.dphi_reaction_max,
    )

    compute_stress_strainrate!(
        g.vx,
        g.vy,
        g.ETA,
        g.GGG,
        g.ETAP,
        g.GGGP,
        g.SXX0,
        g.SXY0,
        g.EXX,
        g.EXY,
        g.SXX,
        g.SXY,
        g.DSXX,
        g.DSXY,
        g.EII,
        g.SII,
        state.dt;
        coords=coords,
    )

    _ = compute_Aϕ!(
        g.APHI,
        g.ETAPHI,
        g.BETAPHI,
        g.PHI,
        g.pr,
        g.pf,
        g.pr0,
        g.pf0,
        state.dt;
        coords=coords,
        betasolid=cur_betasolid,
        phimin=cfg.poroelasticity.phimin,
        phimax=cfg.poroelasticity.phimax,
        S_vent=cfg.venting.active ? g.S_vent_grid : nothing,
    )
    symmetrize_p_node_observables!(g.SXX, g.APHI, g.PHI, g.pr, g.pf, g.ps)

    adjustment_ok = compute_nodal_adjustment!(
        g.ETA,
        g.ETA0,
        g.ETA5,
        g.GGG,
        g.SXX,
        g.SXY,
        g.pr,
        g.pf,
        g.COH,
        g.TEN,
        g.FRI,
        g.YNY,
        g.YNY5,
        ws.YERRNOD,
        g.DSY,
        state.dt,
        iplast;
        etawt=cfg.solver.etawt,
        etamax=cfg.solver.etamax,
        etamin=cfg.solver.etamin,
        yerrmax=cfg.solver.yerrmax,
        max_plastic_iterations=cfg.solver.max_plastic_iterations,
    )

    if !adjustment_ok
        state.dt = finalize_plastic_iteration_pass!(
            g.ETA,
            g.ETA5,
            g.ETA00,
            g.YNY,
            g.YNY5,
            g.YNY00,
            g.YNY_inv_ETA,
            state.dt,
            iplast;
            dtstep=cfg.time.dtstep,
            dtcoefdn=cfg.time.dtcoefdn,
        )
    end

    return (; adjustment_ok, aphimax, n_flips_iter)
end

"""
Solve Stokes-Darcy fluid-solid equations with non-linear plastic iterations.

$(SIGNATURES)
"""
function solve_stokes_darcy!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig,
    ws::SimulationWorkspaces;
    titer::Int,
    dt_step_initial::Float64,
)
    g = state.grids
    g.ETA00 .= g.ETA
    g.YNY00 .= g.YNY
    cur_betasolid = state.timestep == 1 ? 0.0 : cfg.poroelasticity.betasolid
    cur_betafluid = state.timestep == 1 ? 0.0 : cfg.poroelasticity.betafluid
    if state.timestep == 1
        g.BETAPHI .= 0.0
    end

    ws.hydromech.pf_prev_iter .= g.pf
    ws.hydromech.rx_eff .= g.RX
    ws.hydromech.ry_eff .= g.RY
    ws.hydromech.rx_eff_prev .= g.RX
    ws.hydromech.ry_eff_prev .= g.RY

    for j in 1:coords.Nx, i in 1:coords.Ny
        peff_c =
            0.25 * (
                g.pr[i, j] + g.pr[i + 1, j] + g.pr[i, j + 1] + g.pr[i + 1, j + 1] -
                g.pf[i, j] - g.pf[i + 1, j] - g.pf[i, j + 1] - g.pf[i + 1, j + 1]
            )
        ws.fractured_cells_prev[i, j] = is_hydrofracture_breached(peff_c, g.TEN[i, j])
    end

    dt_aphimax_step_max = 0.0
    n_flips_last = 0
    n_flips_step_total = 0
    plastic_converged = false
    max_plastic_iterations_val = cfg.solver.max_plastic_iterations

    for iplast in 1:max_plastic_iterations_val
        assemble_and_solve_hydromechanical!(
            state,
            coords,
            cfg,
            ws;
            titer=titer,
            iplast=iplast,
            cur_betasolid=cur_betasolid,
            cur_betafluid=cur_betafluid,
        )
        res = postprocess_hydromechanical_solution!(
            state,
            coords,
            cfg,
            ws;
            titer=titer,
            iplast=iplast,
            dt_step_initial=dt_step_initial,
            cur_betasolid=cur_betasolid,
            cur_betafluid=cur_betafluid,
        )
        dt_aphimax_step_max = max(dt_aphimax_step_max, state.dt * res.aphimax)
        n_flips_last = res.n_flips_iter
        n_flips_step_total += res.n_flips_iter

        if res.adjustment_ok
            plastic_converged = true
            break
        end
    end

    if !plastic_converged
        last_plastic_residual = ws.YERRNOD[min(
            max_plastic_iterations_val, length(ws.YERRNOD)
        )]
        return (;
            plastic_converged=false,
            last_plastic_residual=last_plastic_residual,
            dt_aphimax_step_max=dt_aphimax_step_max,
            n_flips_last=n_flips_last,
            n_flips_step_total=n_flips_step_total,
        )
    end

    if cfg.venting.active
        apply_venting_surface_boundary!(
            nothing,
            nothing,
            g.tk1,
            coords,
            state.accumulators.rplanet,
            state.accumulators.xcenter,
            state.accumulators.ycenter,
            state.accumulators.P_amb;
            species=cfg.venting.species,
            k_vent=cfg.venting.k_vent,
            conductance_factor=cfg.venting.conductance_factor,
            mode=cfg.venting.mode,
            hydrofracture=cfg.poroelasticity.hydrofracture,
            ice_sealing=cfg.venting.ice_sealing,
            t_freeze=cfg.venting.t_freeze,
            dt_seal=cfg.venting.dt_seal,
            k_seal_min_ratio=cfg.venting.k_seal_min_ratio,
            kappa_frac=cfg.poroelasticity.kappa_frac,
            gamma_frac=cfg.poroelasticity.gamma_frac,
            k_frac_max=cfg.poroelasticity.k_frac_max,
            ramp_width=cfg.poroelasticity.ramp_width,
            pr=g.pr,
            pf=g.pf,
            TEN=g.TEN,
            PHI=g.PHI,
            phimin=cfg.poroelasticity.phimin,
            dt=state.dt,
            eta_fluid_surf=cfg.materials.etafluidmm[2],
            L_sub=cfg.venting.L_sublimation,
            S_vent_out=g.S_vent_grid,
        )
    end

    return (;
        plastic_converged=true,
        last_plastic_residual=0.0,
        dt_aphimax_step_max=dt_aphimax_step_max,
        n_flips_last=n_flips_last,
        n_flips_step_total=n_flips_step_total,
    )
end

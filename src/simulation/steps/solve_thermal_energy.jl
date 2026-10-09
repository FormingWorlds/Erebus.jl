# Thermal energy equation assembly, solve, and convergence assessment step

"""
Solve the thermal energy conservation equation with segregation sources and convection.

$(SIGNATURES)
"""
function solve_thermal_energy!(
    state::SimulationState,
    coords::GridCoordinates,
    cfg::SimulationConfig,
    ws::SimulationWorkspaces;
    titer::Int,
    DHP_pyro=nothing,
)
    g = state.grids
    core = state.markers.core
    marknum = length(state.markers)

    compute_shear_heating!(
        g.HS,
        g.ETA,
        g.SXY,
        g.ETAP,
        g.SXX,
        g.RX,
        g.RY,
        g.qxD,
        g.qyD,
        g.PHI,
        g.ETAPHI,
        g.pr,
        g.pf;
        hydrofracture=cfg.poroelasticity.hydrofracture,
        TEN=g.TEN,
        KX=g.KX,
        KY=g.KY,
        PHIX=g.PHIX,
        PHIY=g.PHIY,
        kphim0=cfg.materials.kphim0,
        phim0_val=cfg.thermodynamics.phim0,
        kappa_frac=cfg.poroelasticity.kappa_frac,
        gamma_frac=cfg.poroelasticity.gamma_frac,
        k_frac_max=cfg.poroelasticity.k_frac_max,
        coords=coords,
    )

    compute_adiabatic_heating!(
        g.HA,
        g.tk1,
        g.ALPHA,
        g.ALPHAF,
        g.PHI,
        g.vx,
        g.vy,
        g.vxf,
        g.vyf,
        g.ps,
        g.pf;
        coords=coords,
    )

    if cfg.geometry.spherical_metric && g.Q_metric !== nothing
        compute_spherical_metric_heat_source!(
            g.Q_metric,
            g.tk1,
            g.KX,
            g.KY,
            coords;
            xcenter=state.accumulators.xcenter,
            ycenter=state.accumulators.ycenter,
            rplanet=state.accumulators.rplanet,
            reg_cells=cfg.geometry.metric_regularization_cells,
        )
    end

    fill!(g.Q_seg_grid, 0.0)

    coreformation_active_val =
        cfg.coreformation.percolation_active ||
        cfg.coreformation.settling_active ||
        cfg.metal_partition.active

    if coreformation_active_val && haskey(state.markers.groups, :metal)
        metal = state.markers.groups.metal
        if ws.Xfe_bulk_step_start !== nothing
            if length(metal.Xfe_bulk) != length(ws.Xfe_bulk_step_start)
                resize!(metal.Xfe_bulk, length(ws.Xfe_bulk_step_start))
            end
            copyto!(metal.Xfe_bulk, ws.Xfe_bulk_step_start)
        end
        if ws.Xfem_step_start !== nothing
            if length(metal.Xfem) != length(ws.Xfem_step_start)
                resize!(metal.Xfem, length(ws.Xfem_step_start))
            end
            copyto!(metal.Xfem, ws.Xfem_step_start)
        end

        if cfg.metal_partition.active
            for (dst, src) in (
                (metal.Xfe_H_m, ws.Xfe_H_m_step_start),
                (metal.Xfe_C_m, ws.Xfe_C_m_step_start),
                (metal.Xfe_N_m, ws.Xfe_N_m_step_start),
                (metal.Xfe_S_m, ws.Xfe_S_m_step_start),
            )
                if src !== nothing && dst !== nothing
                    if length(dst) != length(src)
                        resize!(dst, length(src))
                    end
                    copyto!(dst, src)
                end
            end
        end

        seg_res = apply_metal_segregation!(
            core.xm,
            core.ym,
            core.tm,
            core.tkm,
            core.phim,
            metal.Xfe_bulk,
            metal.Xfem,
            marknum,
            state.dt,
            cfg.coreformation;
            coords=coords,
            xcenter=state.accumulators.xcenter,
            ycenter=state.accumulators.ycenter,
            rplanet=state.accumulators.rplanet,
            gx=g.gx,
            gy=g.gy,
            Q_seg_grid=cfg.coreformation.segregation_heating ? g.Q_seg_grid : nothing,
            rho_silicate=cfg.materials.rhosolidm[1],
            eta_silicate=cfg.materials.etasolidm[1],
            ETA=g.ETA,
            Fm=core.Fm,
            T_solidus_silicate=cfg.melting.T_solidus[1],
            T_liquidus_silicate=cfg.melting.T_liquidus[1],
            cfg_partition=cfg.metal_partition,
            Xfe_H_m=metal.Xfe_H_m,
            Xfe_C_m=metal.Xfe_C_m,
            Xfe_N_m=metal.Xfe_N_m,
            Xfe_S_m=metal.Xfe_S_m,
            workspace=ws.metal_segregation,
        )
        state.accumulators.max_v_seg_prev = seg_res.max_v_seg
    end

    if cfg.venting.active && cfg.venting.latent_cooling
        @. g.Q_lat_grid =
            -cfg.venting.L_sublimation * cfg.materials.rhofluidm[2] * g.S_vent_grid
    else
        fill!(g.Q_lat_grid, 0.0)
    end

    magma_active_val = cfg.magma_transport.active
    if magma_active_val
        if ws.Fm_step_start !== nothing
            if length(core.Fm) != length(ws.Fm_step_start)
                resize!(core.Fm, length(ws.Fm_step_start))
            end
            copyto!(core.Fm, ws.Fm_step_start)
        end

        F_extract_m = nothing
        if haskey(state.markers.groups, :volatiles)
            vols = state.markers.groups.volatiles
            if ws.F_extract_m_step_start !== nothing && vols.F_extract_m !== nothing
                if length(vols.F_extract_m) != length(ws.F_extract_m_step_start)
                    resize!(vols.F_extract_m, length(ws.F_extract_m_step_start))
                end
                copyto!(vols.F_extract_m, ws.F_extract_m_step_start)
            end
            F_extract_m = vols.F_extract_m
        end

        XH2Om = if haskey(state.markers.groups, :volatiles)
            state.markers.groups.volatiles.XH2Om
        else
            nothing
        end
        XCm = if haskey(state.markers.groups, :volatiles)
            state.markers.groups.volatiles.XCm
        else
            nothing
        end
        XNm = if haskey(state.markers.groups, :volatiles)
            state.markers.groups.volatiles.XNm
        else
            nothing
        end
        XSm = if haskey(state.markers.groups, :volatiles)
            state.markers.groups.volatiles.XSm
        else
            nothing
        end

        apply_silicate_melt_segregation!(
            core.xm,
            core.ym,
            core.tm,
            core.tkm,
            core.Fm,
            marknum,
            state.dt,
            cfg.magma_transport;
            coords=coords,
            xcenter=state.accumulators.xcenter,
            ycenter=state.accumulators.ycenter,
            rplanet=state.accumulators.rplanet,
            gx=g.gx,
            gy=g.gy,
            Q_seg_grid=if (
                cfg.magma_transport.segregation_heating ||
                cfg.magma_transport.sensible_heat_transport
            )
                g.Q_seg_grid
            else
                nothing
            end,
            Q_lat_grid=cfg.magma_transport.latent_crystallization ? g.Q_lat_grid : nothing,
            rho_silicate=cfg.materials.rhosolidm[1],
            rho_melt=cfg.melting.rho_melt,
            eta_silicate=cfg.materials.etasolidm[1],
            ETA=g.ETA,
            T_solidus_silicate=cfg.melting.T_solidus[1],
            T_liquidus_silicate=cfg.melting.T_liquidus[1],
            L_melt=cfg.melting.L_melt,
            F_extract_m=F_extract_m,
            vx=g.vx,
            vy=g.vy,
            pr=g.pr,
            XH2Om=cfg.volatiles.active ? XH2Om : nothing,
            XCm=cfg.volatiles.active ? XCm : nothing,
            XNm=cfg.volatiles.active ? XNm : nothing,
            XSm=cfg.volatiles.active ? XSm : nothing,
            phim=core.phim,
            cfg_volatiles=cfg.volatiles,
            workspace=ws.magma_segregation,
        )
    end

    Q_seg_val =
        if (coreformation_active_val && cfg.coreformation.segregation_heating) || (
            magma_active_val && (
                cfg.magma_transport.segregation_heating ||
                cfg.magma_transport.sensible_heat_transport
            )
        )
            g.Q_seg_grid
        else
            nothing
        end

    LT = assemble_thermal_lse!(
        g.tk1,
        g.RHOCP,
        g.KX,
        g.KY,
        g.HR,
        g.HA,
        g.HS,
        g.DHP,
        ws.RT,
        state.dt;
        coords=coords,
        LT=ws.thermal.LT,
        Q_metric=g.Q_metric,
        Q_lat=g.Q_lat_grid,
        Q_seg=Q_seg_val,
        workspace=ws.thermal,
    )

    if ws.thermal_cache === nothing
        thermal_prob = LinearProblem(LT.cscmatrix, ws.RT)
        ws.thermal_cache = init(thermal_prob, UMFPACKFactorization(; reuse_symbolic=true))
    else
        ws.thermal_cache.A = LT.cscmatrix
        ws.thermal_cache.b = ws.RT
    end
    thermal_sol = solve!(ws.thermal_cache)
    if !LinearSolve.SciMLBase.successful_retcode(thermal_sol) ||
        !all(isfinite, thermal_sol.u)
        error("Thermal solver failed with retcode $(thermal_sol.retcode)")
    end

    ws.ST .= thermal_sol.u
    g.tk2 .= reshape(ws.ST, coords.Ny1, coords.Nx1)
    @. g.DT = g.tk2 - g.tk1
    maxDTcurrent = maximum(abs, g.DT)

    dt_next = finalize_thermochemical_iteration_pass(
        maxDTcurrent, state.dt, titer, cfg.time.DTmax
    )
    dt_reduced_by_maxDT = (dt_next < state.dt)

    thermochemical_converged = compute_thermochemical_iteration_outcome(
        g.DMP,
        g.pf,
        ws.hydromech.pf_prev_iter,
        titer;
        pferrmax=cfg.reaction.pferrmax,
        maxDTcurrent=maxDTcurrent,
        DTmax=cfg.time.DTmax,
    )

    return (; thermochemical_converged, maxDTcurrent, dt_next, dt_reduced_by_maxDT)
end

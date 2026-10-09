# Domain expansion and telescoping grid/marker remapping step

"""
Execute domain telescoping expansion and remap grid arrays, markers, and workspaces.

$(SIGNATURES)
"""
function telescope_domain!(
    state::SimulationState,
    coords_ref::Ref{GridCoordinates},
    cfg::SimulationConfig,
    ws::SimulationWorkspaces,
)::Bool
    if !cfg.telescoping.active ||
        !should_telescope_domain(
        state.accumulators.rplanet,
        coords_ref[],
        cfg.telescoping;
        level=state.accumulators.telescope_level,
    )
        return false
    end

    old_coords = coords_ref[]
    coords_new = compute_telescoped_coordinates(old_coords)
    @info "Triggering domain telescoping expansion" level=state.accumulators.telescope_level +
                                                          1 rplanet=state.accumulators.rplanet old_xsize=old_coords.xsize new_xsize=coords_new.xsize

    g = state.grids
    dim_cell = (coords_new.Ny, coords_new.Nx)
    dim_node = (coords_new.Ny1, coords_new.Nx1)

    ETA = remap_staggered_grid_array(
        g.ETA, dim_cell; background_val=cfg.materials.etasolidm[3]
    )
    ETA0 = remap_staggered_grid_array(
        g.ETA0, dim_cell; background_val=cfg.materials.etasolidm[3]
    )
    GGG = remap_staggered_grid_array(
        g.GGG, dim_cell; background_val=cfg.materials.gggsolidm[3]
    )
    EXY = remap_staggered_grid_array(g.EXY, dim_cell)
    SXY = remap_staggered_grid_array(g.SXY, dim_cell)
    SXY0 = remap_staggered_grid_array(g.SXY0, dim_cell)
    wyx = remap_staggered_grid_array(g.wyx, dim_cell)
    COH = remap_staggered_grid_array(
        g.COH, dim_cell; background_val=cfg.materials.cohessolidm[3]
    )
    TEN = remap_staggered_grid_array(
        g.TEN, dim_cell; background_val=cfg.materials.tenssolidm[3]
    )
    FRI = remap_staggered_grid_array(
        g.FRI, dim_cell; background_val=cfg.materials.frictsolidm[3]
    )
    YNY = remap_staggered_grid_array(g.YNY, dim_cell)
    ETA5 = remap_staggered_grid_array(
        g.ETA5, dim_cell; background_val=cfg.materials.etasolidm[3]
    )
    ETA00 = remap_staggered_grid_array(
        g.ETA00, dim_cell; background_val=cfg.materials.etasolidm[3]
    )
    YNY5 = remap_staggered_grid_array(g.YNY5, dim_cell)
    YNY00 = remap_staggered_grid_array(g.YNY00, dim_cell)
    YNY_inv_ETA = remap_staggered_grid_array(g.YNY_inv_ETA, dim_cell)
    DSXY = remap_staggered_grid_array(g.DSXY, dim_cell)
    DSY = remap_staggered_grid_array(g.DSY, dim_cell)

    RHOX = remap_staggered_grid_array(
        g.RHOX, dim_node; background_val=cfg.materials.rhosolidm[3]
    )
    RHOFX = remap_staggered_grid_array(
        g.RHOFX, dim_node; background_val=cfg.materials.rhofluidm[3]
    )
    KX = remap_staggered_grid_array(g.KX, dim_node; background_val=cfg.materials.ksolidm[3])
    PHIX = remap_staggered_grid_array(
        g.PHIX, dim_node; background_val=cfg.poroelasticity.phimin
    )
    vx = remap_staggered_grid_array(g.vx, dim_node)
    vxf = remap_staggered_grid_array(g.vxf, dim_node)
    RX = remap_staggered_grid_array(g.RX, dim_node)
    qxD = remap_staggered_grid_array(g.qxD, dim_node)
    gx = remap_staggered_grid_array(g.gx, dim_node)

    RHOY = remap_staggered_grid_array(
        g.RHOY, dim_node; background_val=cfg.materials.rhosolidm[3]
    )
    RHOFY = remap_staggered_grid_array(
        g.RHOFY, dim_node; background_val=cfg.materials.rhofluidm[3]
    )
    KY = remap_staggered_grid_array(g.KY, dim_node; background_val=cfg.materials.ksolidm[3])
    PHIY = remap_staggered_grid_array(
        g.PHIY, dim_node; background_val=cfg.poroelasticity.phimin
    )
    vy = remap_staggered_grid_array(g.vy, dim_node)
    vyf = remap_staggered_grid_array(g.vyf, dim_node)
    RY = remap_staggered_grid_array(g.RY, dim_node)
    qyD = remap_staggered_grid_array(g.qyD, dim_node)
    gy = remap_staggered_grid_array(g.gy, dim_node)

    RHO = remap_staggered_grid_array(
        g.RHO, dim_node; background_val=cfg.materials.rhosolidm[3]
    )
    RHOCP = remap_staggered_grid_array(
        g.RHOCP, dim_node; background_val=cfg.materials.rhocpsolidm[3]
    )
    ALPHA = remap_staggered_grid_array(
        g.ALPHA, dim_node; background_val=cfg.materials.alphasolidm[3]
    )
    ALPHAF = remap_staggered_grid_array(
        g.ALPHAF, dim_node; background_val=cfg.materials.alphafluidm[3]
    )
    HR = remap_staggered_grid_array(g.HR, dim_node)
    HA = remap_staggered_grid_array(g.HA, dim_node)
    HS = remap_staggered_grid_array(g.HS, dim_node)
    ETAP = remap_staggered_grid_array(
        g.ETAP, dim_node; background_val=cfg.materials.etasolidm[3]
    )
    GGGP = remap_staggered_grid_array(
        g.GGGP, dim_node; background_val=cfg.materials.gggsolidm[3]
    )
    EXX = remap_staggered_grid_array(g.EXX, dim_node)
    SXX = remap_staggered_grid_array(g.SXX, dim_node)
    SXX0 = remap_staggered_grid_array(g.SXX0, dim_node)
    tk1 = remap_staggered_grid_array(g.tk1, dim_node; background_val=cfg.materials.tkm0[3])
    tk2 = remap_staggered_grid_array(g.tk2, dim_node; background_val=cfg.materials.tkm0[3])
    DT = remap_staggered_grid_array(g.DT, dim_node)
    DT0 = remap_staggered_grid_array(g.DT0, dim_node)
    vxp = remap_staggered_grid_array(g.vxp, dim_node)
    vyp = remap_staggered_grid_array(g.vyp, dim_node)
    vxpf = remap_staggered_grid_array(g.vxpf, dim_node)
    vypf = remap_staggered_grid_array(g.vypf, dim_node)
    pr = remap_staggered_grid_array(g.pr, dim_node)
    pf = remap_staggered_grid_array(g.pf, dim_node)
    ps = remap_staggered_grid_array(g.ps, dim_node)
    pr0 = remap_staggered_grid_array(g.pr0, dim_node)
    pf0 = remap_staggered_grid_array(g.pf0, dim_node)
    ps0 = remap_staggered_grid_array(g.ps0, dim_node)
    ETAPHI = remap_staggered_grid_array(
        g.ETAPHI, dim_node; background_val=cfg.materials.etasolidm[3]
    )
    BETAPHI = remap_staggered_grid_array(g.BETAPHI, dim_node)
    PHI = remap_staggered_grid_array(
        g.PHI, dim_node; background_val=cfg.poroelasticity.phimin
    )
    APHI = remap_staggered_grid_array(g.APHI, dim_node)
    FI = remap_staggered_grid_array(g.FI, dim_node)
    DMP = remap_staggered_grid_array(g.DMP, dim_node)
    DHP = remap_staggered_grid_array(g.DHP, dim_node)
    XWS = remap_staggered_grid_array(g.XWS, dim_node)
    EII = remap_staggered_grid_array(g.EII, dim_node)
    SII = remap_staggered_grid_array(g.SII, dim_node)
    DSXX = remap_staggered_grid_array(g.DSXX, dim_node)
    tk0 = remap_staggered_grid_array(g.tk0, dim_node; background_val=cfg.materials.tkm0[3])
    Q_metric =
        g.Q_metric !== nothing ? remap_staggered_grid_array(g.Q_metric, dim_node) : nothing
    DQPF = remap_staggered_grid_array(g.DQPF, dim_node)
    DQPFSUM = remap_staggered_grid_array(g.DQPFSUM, dim_node)
    S_vent_grid = remap_staggered_grid_array(g.S_vent_grid, dim_node)
    Q_lat_grid = remap_staggered_grid_array(g.Q_lat_grid, dim_node)
    Q_seg_grid = remap_staggered_grid_array(g.Q_seg_grid, dim_node)

    state.grids = GridArrays(
        ETA,
        ETA0,
        GGG,
        EXY,
        SXY,
        SXY0,
        wyx,
        COH,
        TEN,
        FRI,
        YNY,
        RHOX,
        RHOFX,
        KX,
        PHIX,
        vx,
        vxf,
        RX,
        qxD,
        gx,
        RHOY,
        RHOFY,
        KY,
        PHIY,
        vy,
        vyf,
        RY,
        qyD,
        gy,
        RHO,
        RHOCP,
        ALPHA,
        ALPHAF,
        HR,
        HA,
        HS,
        ETAP,
        GGGP,
        EXX,
        SXX,
        SXX0,
        tk1,
        tk2,
        DT,
        DT0,
        vxp,
        vyp,
        vxpf,
        vypf,
        pr,
        pf,
        ps,
        pr0,
        pf0,
        ps0,
        ETAPHI,
        BETAPHI,
        PHI,
        APHI,
        FI,
        DMP,
        DHP,
        XWS,
        ETA5,
        ETA00,
        YNY5,
        YNY00,
        YNY_inv_ETA,
        DSXY,
        DSY,
        EII,
        SII,
        DSXX,
        tk0,
        DQPF,
        DQPFSUM,
        S_vent_grid,
        Q_lat_grid,
        Q_seg_grid,
        Q_metric,
    )

    telescope_marker_arrays!(
        state.markers;
        old_coords=old_coords,
        new_coords=coords_new,
        buffer_markers_per_cell=cfg.telescoping.buffer_markers_per_cell,
        cfg=cfg,
    )

    reset_workspaces_for_grid!(ws, coords_new, cfg, length(state.markers))

    coords_ref[] = coords_new
    state.accumulators.telescope_level += 1
    state.accumulators.xcenter = coords_new.xcenter
    state.accumulators.ycenter = coords_new.ycenter

    @info "Telescoping complete" level=state.accumulators.telescope_level Nx=coords_new.Nx Ny=coords_new.Ny marknum=length(
        state.markers
    )
    return true
end

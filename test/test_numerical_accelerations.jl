using Test
using Random
using SparseArrays
using ExtendableSparse
using Erebus
using Erebus.Numerics: assemble_hydromechanical_lse!, assemble_thermal_lse!
using Erebus.Particles:
    update_marker_temperature!,
    update_marker_porosity!,
    sink_vented_marker_porosity!,
    drain_vented_marker_volatiles!,
    advance_marker_thermo_porosity_venting!

@testset "Numerical Accelerations & Workspaces" begin
    @testset "Workspace Submodule Exports" begin
        @test :HydromechanicalLSEWorkspace in names(Erebus.Numerics)
        @test :ThermalLSEWorkspace in names(Erebus.Numerics)
        @test :HydromechanicalLSEWorkspace in names(Erebus.Simulation)
        @test :ThermalLSEWorkspace in names(Erebus.Simulation)
        @test :HydromechanicalLSEWorkspace in names(Erebus)
        @test :ThermalLSEWorkspace in names(Erebus)
        @test :advance_marker_thermo_porosity_venting! in names(Erebus.Particles)
        @test :advance_marker_thermo_porosity_venting! in names(Erebus)
    end

    @testset "LSE Workspace Constructors & Dimensions" begin
        Ny1, Nx1 = 15, 15
        ws_h = HydromechanicalLSEWorkspace(Ny1, Nx1)
        @test ws_h.Ny1 == Ny1
        @test ws_h.Nx1 == Nx1
        @test size(ws_h.L) == (Ny1 * Nx1 * 6, Ny1 * Nx1 * 6)
        @test !ws_h.is_initialized

        ws_t = ThermalLSEWorkspace(Ny1, Nx1)
        @test ws_t.Ny1 == Ny1
        @test ws_t.Nx1 == Nx1
        @test size(ws_t.LT) == (Ny1 * Nx1, Ny1 * Nx1)
        @test !ws_t.is_initialized

        coords = GridCoordinates(GridConfig(; Nx=12, Ny=12, xsize=1.0e5, ysize=1.0e5))
        ws_coords_h = HydromechanicalLSEWorkspace(coords)
        ws_coords_t = ThermalLSEWorkspace(coords)
        @test ws_coords_h.Ny1 == coords.Ny1
        @test ws_coords_h.Nx1 == coords.Nx1
        @test ws_coords_t.Ny1 == coords.Ny1
        @test ws_coords_t.Nx1 == coords.Nx1
    end

    @testset "Stokes-Darcy LSE Assembly Equivalence & Workspace Reuse" begin
        coords = GridCoordinates(GridConfig(; Nx=8, Ny=8, xsize=5.0e4, ysize=5.0e4))
        Ny1 = coords.Ny1
        Nx1 = coords.Nx1

        ETA = fill(1.0e20, Ny1, Nx1)
        ETAP = fill(1.0e20, Ny1, Nx1)
        GGG = fill(1.0e10, Ny1, Nx1)
        GGGP = fill(1.0e10, Ny1, Nx1)
        SXY0 = zeros(Ny1, Nx1)
        SXX0 = zeros(Ny1, Nx1)
        RHOX = fill(3000.0, Ny1, Nx1)
        RHOY = fill(3000.0, Ny1, Nx1)
        RHOFX = fill(1000.0, Ny1, Nx1)
        RHOFY = fill(1000.0, Ny1, Nx1)
        RX = zeros(Ny1, Nx1)
        RY = zeros(Ny1, Nx1)
        ETAPHI = fill(1.0e18, Ny1, Nx1)
        BETAPHI = fill(1.0e-11, Ny1, Nx1)
        PHI = fill(0.05, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = fill(1.62, Ny1, Nx1)
        pr0 = fill(1.0e7, Ny1, Nx1)
        pf0 = fill(8.0e6, Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)
        dt = 3.15e7
        R1 = zeros(Nx1 * Ny1 * 6)
        R2 = zeros(Nx1 * Ny1 * 6)
        R3 = zeros(Nx1 * Ny1 * 6)

        # Baseline: assembly without workspace
        L1 = assemble_hydromechanical_lse!(
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
            R1;
            coords=coords,
            workspace=nothing,
        )

        # Accelerated pass 1: initialize workspace
        ws = HydromechanicalLSEWorkspace(coords)
        @test !ws.is_initialized
        L2 = assemble_hydromechanical_lse!(
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
            R2;
            coords=coords,
            workspace=ws,
        )
        @test ws.is_initialized
        @test L1.rowval == L2.rowval
        @test L1.colptr == L2.colptr
        @test maximum(abs.(L1.nzval .- L2.nzval)) ≈ 0.0 atol=1e-12
        @test maximum(abs.(R1 .- R2)) ≈ 0.0 atol=1e-12

        # Accelerated pass 2: reuse pre-allocated matrix pattern with updated dt
        dt_new = 5.0e7
        R_fresh = zeros(Nx1 * Ny1 * 6)
        L_fresh = assemble_hydromechanical_lse!(
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
            dt_new,
            R_fresh;
            coords=coords,
            workspace=nothing,
        )

        L3 = assemble_hydromechanical_lse!(
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
            dt_new,
            R3;
            coords=coords,
            workspace=ws,
        )
        @test L3.rowval == L_fresh.rowval
        @test L3.colptr == L_fresh.colptr
        @test maximum(abs.(L3.nzval .- L_fresh.nzval)) ≈ 0.0 atol=1e-12
        @test maximum(abs.(R3 .- R_fresh)) ≈ 0.0 atol=1e-12
    end

    @testset "Stokes-Darcy Venting Boundary Assembly Equivalence & Reuse" begin
        coords = GridCoordinates(GridConfig(; Nx=8, Ny=8, xsize=5.0e4, ysize=5.0e4))
        Ny1, Nx1 = coords.Ny1, coords.Nx1
        ETA = fill(1.0e20, Ny1, Nx1)
        ETAP = fill(1.0e20, Ny1, Nx1)
        GGG = fill(1.0e10, Ny1, Nx1)
        GGGP = fill(1.0e10, Ny1, Nx1)
        SXY0 = zeros(Ny1, Nx1)
        SXX0 = zeros(Ny1, Nx1)
        RHOX = fill(3000.0, Ny1, Nx1)
        RHOY = fill(3000.0, Ny1, Nx1)
        RHOFX = fill(1000.0, Ny1, Nx1)
        RHOFY = fill(1000.0, Ny1, Nx1)
        RX = zeros(Ny1, Nx1)
        RY = zeros(Ny1, Nx1)
        ETAPHI = fill(1.0e18, Ny1, Nx1)
        BETAPHI = fill(1.0e-11, Ny1, Nx1)
        PHI = fill(0.05, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = fill(1.62, Ny1, Nx1)
        pr0 = fill(1.0e7, Ny1, Nx1)
        pf0 = fill(8.0e6, Ny1, Nx1)
        pr = fill(1.0e7, Ny1, Nx1)
        pf = fill(8.0e6, Ny1, Nx1)
        TEN = fill(1.0e6, Ny1, Nx1)
        tk = fill(300.0, Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)
        dt = 3.15e7
        R_ref = zeros(Nx1 * Ny1 * 6)
        R_ws = zeros(Nx1 * Ny1 * 6)

        L_ref = assemble_hydromechanical_lse!(
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
            R_ref;
            coords=coords,
            pr=pr,
            pf=pf,
            TEN=TEN,
            tk=tk,
            venting=true,
            rplanet=2.0e4,
            xcenter=2.5e4,
            ycenter=2.5e4,
            P_amb=1.0e5,
            workspace=nothing,
        )

        ws_v = HydromechanicalLSEWorkspace(coords)
        L_ws = assemble_hydromechanical_lse!(
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
            R_ws;
            coords=coords,
            pr=pr,
            pf=pf,
            TEN=TEN,
            tk=tk,
            venting=true,
            rplanet=2.0e4,
            xcenter=2.5e4,
            ycenter=2.5e4,
            P_amb=1.0e5,
            workspace=ws_v,
        )

        @test L_ws.rowval == L_ref.rowval
        @test L_ws.colptr == L_ref.colptr
        @test maximum(abs.(L_ws.nzval .- L_ref.nzval)) ≈ 0.0 atol=1e-12
        @test maximum(abs.(R_ws .- R_ref)) ≈ 0.0 atol=1e-12
    end

    @testset "Thermal LSE Assembly Equivalence & Workspace Reuse" begin
        coords = GridCoordinates(GridConfig(; Nx=8, Ny=8, xsize=5.0e4, ysize=5.0e4))
        Ny1 = coords.Ny1
        Nx1 = coords.Nx1

        tk1 = fill(800.0, Ny1, Nx1)
        RHOCP = fill(3.0e6, Ny1, Nx1)
        KX = fill(3.0, Ny1, Nx1)
        KY = fill(3.0, Ny1, Nx1)
        HR = fill(1.0e-8, Ny1, Nx1)
        HA = zeros(Ny1, Nx1)
        HS = zeros(Ny1, Nx1)
        DHP = zeros(Ny1, Nx1)
        dt = 3.15e7

        RT1 = zeros(Ny1 * Nx1)
        RT2 = zeros(Ny1 * Nx1)
        RT3 = zeros(Ny1 * Nx1)

        # Baseline: assembly without workspace
        LT1 = assemble_thermal_lse!(
            tk1, RHOCP, KX, KY, HR, HA, HS, DHP, RT1, dt; coords=coords, workspace=nothing
        )

        # Accelerated pass 1: initialize workspace
        ws_t = ThermalLSEWorkspace(coords)
        @test !ws_t.is_initialized
        LT2 = assemble_thermal_lse!(
            tk1, RHOCP, KX, KY, HR, HA, HS, DHP, RT2, dt; coords=coords, workspace=ws_t
        )
        @test ws_t.is_initialized
        @test LT1.cscmatrix.rowval == LT2.cscmatrix.rowval
        @test LT1.cscmatrix.colptr == LT2.cscmatrix.colptr
        @test maximum(abs.(LT1.cscmatrix.nzval .- LT2.cscmatrix.nzval)) ≈ 0.0 atol=1e-12
        @test maximum(abs.(RT1 .- RT2)) ≈ 0.0 atol=1e-12

        # Accelerated pass 2: reuse pre-allocated matrix pattern with updated temperature
        tk1_new = fill(950.0, Ny1, Nx1)
        RT_fresh = zeros(Ny1 * Nx1)
        LT_fresh = assemble_thermal_lse!(
            tk1_new,
            RHOCP,
            KX,
            KY,
            HR,
            HA,
            HS,
            DHP,
            RT_fresh,
            dt;
            coords=coords,
            workspace=nothing,
        )

        LT3 = assemble_thermal_lse!(
            tk1_new, RHOCP, KX, KY, HR, HA, HS, DHP, RT3, dt; coords=coords, workspace=ws_t
        )
        @test LT3.cscmatrix.rowval == LT_fresh.cscmatrix.rowval
        @test LT3.cscmatrix.colptr == LT_fresh.cscmatrix.colptr
        @test maximum(abs.(LT3.cscmatrix.nzval .- LT_fresh.cscmatrix.nzval)) ≈ 0.0 atol=1e-12
        @test maximum(abs.(RT3 .- RT_fresh)) ≈ 0.0 atol=1e-12
    end

    @testset "Fused Marker Thermo-Porosity-Venting Equivalence" begin
        coords = GridCoordinates(GridConfig(; Nx=12, Ny=12, xsize=6.0e4, ysize=6.0e4))
        marknum = 120
        rng = MersenneTwister(123)

        xm_base = 5.0e3 .+ 5.0e4 .* rand(rng, marknum)
        ym_base = 5.0e3 .+ 5.0e4 .* rand(rng, marknum)
        tm_base = rand(rng, [1, 2, 3], marknum)
        tkm_base = 700.0 .+ 200.0 .* rand(rng, marknum)
        phim_base = 0.01 .+ 0.15 .* rand(rng, marknum)
        XH2Om_base = 0.5 .+ 1.5 .* rand(rng, marknum)
        XCm_base = 100.0 .+ 200.0 .* rand(rng, marknum)
        XNm_base = 10.0 .+ 30.0 .* rand(rng, marknum)
        XSm_base = 50.0 .+ 100.0 .* rand(rng, marknum)
        Fm_base = 0.05 .* rand(rng, marknum)

        DT = fill(5.0, coords.Ny1, coords.Nx1)
        tk2 = fill(750.0, coords.Ny1, coords.Nx1)
        APHI = fill(-1.0e-12, coords.Ny1, coords.Nx1)
        S_vent_grid = zeros(coords.Ny1, coords.Nx1)
        # Active venting in surface region
        S_vent_grid[1:3, :] .= 1.0e-10

        ret_cfg = RetentionConfig(; active=true, venting_drainage_active=true, chi_vent=1.0)
        dt = 3.15e7

        # Serial reference pass
        xm1 = copy(xm_base)
        ym1 = copy(ym_base)
        tm1 = copy(tm_base)
        tkm1 = copy(tkm_base)
        phim1 = copy(phim_base)
        XH2Om1 = copy(XH2Om_base)
        XCm1 = copy(XCm_base)
        XNm1 = copy(XNm_base)
        XSm1 = copy(XSm_base)
        Fm1 = copy(Fm_base)

        update_marker_temperature!(xm1, ym1, tkm1, DT, tk2, 2, marknum; coords=coords)
        update_marker_porosity!(
            xm1,
            ym1,
            tm1,
            phim1,
            APHI,
            dt,
            marknum;
            phimin=1.0e-4,
            phimax=1.0,
            coords=coords,
        )
        m_vent_ref = sink_vented_marker_porosity!(
            xm1,
            ym1,
            tm1,
            phim1,
            S_vent_grid,
            dt,
            marknum;
            coords=coords,
            phimin=1.0e-4,
            rhofluidcur=1000.0,
        )
        vols_ref = drain_vented_marker_volatiles!(
            xm1,
            ym1,
            tm1,
            tkm1,
            XH2Om1,
            XCm1,
            XNm1,
            XSm1,
            S_vent_grid,
            dt,
            marknum,
            ret_cfg;
            coords=coords,
            rhosolid=3000.0,
            phim=phim1,
            Fm=Fm1,
        )

        # Fused single-pass
        xm2 = copy(xm_base)
        ym2 = copy(ym_base)
        tm2 = copy(tm_base)
        tkm2 = copy(tkm_base)
        phim2 = copy(phim_base)
        XH2Om2 = copy(XH2Om_base)
        XCm2 = copy(XCm_base)
        XNm2 = copy(XNm_base)
        XSm2 = copy(XSm_base)
        Fm2 = copy(Fm_base)

        fused_res = advance_marker_thermo_porosity_venting!(
            xm2,
            ym2,
            tm2,
            tkm2,
            phim2,
            DT,
            tk2,
            APHI,
            dt,
            2,
            marknum;
            coords=coords,
            phimin=1.0e-4,
            phimax=1.0,
            venting=true,
            S_vent_grid=S_vent_grid,
            rhofluidcur=1000.0,
            ret_cfg=ret_cfg,
            XH2Om=XH2Om2,
            XCm=XCm2,
            XNm=XNm2,
            XSm=XSm2,
            rhosolid=3000.0,
            Fm=Fm2,
        )

        @test tkm1 ≈ tkm2
        @test phim1 ≈ phim2
        @test fused_res.delta_m_vent ≈ m_vent_ref
        @test fused_res.vented_vols.M_vent_H2O ≈ vols_ref.M_vent_H2O
        @test fused_res.vented_vols.M_vent_C ≈ vols_ref.M_vent_C
        @test fused_res.vented_vols.M_vent_N ≈ vols_ref.M_vent_N
        @test fused_res.vented_vols.M_vent_S ≈ vols_ref.M_vent_S
        @test XH2Om1 ≈ XH2Om2
        @test XCm1 ≈ XCm2
        @test XNm1 ≈ XNm2
        @test XSm1 ≈ XSm2

        # Venting without volatiles (venting=true, XH2Om=nothing)
        phim_novol1 = copy(phim_base)
        update_marker_porosity!(
            xm1,
            ym1,
            tm1,
            phim_novol1,
            APHI,
            dt,
            marknum;
            phimin=1.0e-4,
            phimax=1.0,
            coords=coords,
        )
        m_vent_novol_ref = sink_vented_marker_porosity!(
            xm1,
            ym1,
            tm1,
            phim_novol1,
            S_vent_grid,
            dt,
            marknum;
            coords=coords,
            phimin=1.0e-4,
            rhofluidcur=1000.0,
        )
        phim_novol2 = copy(phim_base)
        tkm_novol2 = copy(tkm_base)
        fused_novol_res = advance_marker_thermo_porosity_venting!(
            xm2,
            ym2,
            tm2,
            tkm_novol2,
            phim_novol2,
            DT,
            tk2,
            APHI,
            dt,
            2,
            marknum;
            coords=coords,
            phimin=1.0e-4,
            phimax=1.0,
            venting=true,
            S_vent_grid=S_vent_grid,
            rhofluidcur=1000.0,
            ret_cfg=ret_cfg,
            XH2Om=nothing,
        )
        @test phim_novol1 ≈ phim_novol2
        @test fused_novol_res.delta_m_vent ≈ m_vent_novol_ref
        @test fused_novol_res.vented_vols === nothing

        # Initial timestep 1 check
        tkm_init1 = copy(tkm_base)
        tkm_init2 = copy(tkm_base)
        update_marker_temperature!(xm1, ym1, tkm_init1, DT, tk2, 1, marknum; coords=coords)
        advance_marker_thermo_porosity_venting!(
            xm2,
            ym2,
            tm2,
            tkm_init2,
            phim2,
            DT,
            tk2,
            APHI,
            dt,
            1,
            marknum;
            coords=coords,
            venting=false,
        )
        @test tkm_init1 ≈ tkm_init2

        # Empty markers edge case
        empty_res = advance_marker_thermo_porosity_venting!(
            Float64[],
            Float64[],
            Int64[],
            Float64[],
            Float64[],
            DT,
            tk2,
            APHI,
            dt,
            1,
            0;
            coords=coords,
        )
        @test empty_res.delta_m_vent ≈ 0.0
        @test empty_res.vented_vols === nothing
    end

    @testset "Telescoping Domain Reallocation Safety" begin
        coords1 = GridCoordinates(GridConfig(; Nx=8, Ny=8, xsize=5.0e4, ysize=5.0e4))
        ws_h1 = HydromechanicalLSEWorkspace(coords1)
        ws_t1 = ThermalLSEWorkspace(coords1)

        coords2 = GridCoordinates(GridConfig(; Nx=16, Ny=16, xsize=1.0e5, ysize=1.0e5))
        ws_h2 = HydromechanicalLSEWorkspace(coords2)
        ws_t2 = ThermalLSEWorkspace(coords2)

        @test size(ws_h1.L) ==
            (coords1.Ny1 * coords1.Nx1 * 6, coords1.Ny1 * coords1.Nx1 * 6)
        @test size(ws_h2.L) ==
            (coords2.Ny1 * coords2.Nx1 * 6, coords2.Ny1 * coords2.Nx1 * 6)
        @test size(ws_t1.LT) == (coords1.Ny1 * coords1.Nx1, coords1.Ny1 * coords1.Nx1)
        @test size(ws_t2.LT) == (coords2.Ny1 * coords2.Nx1, coords2.Ny1 * coords2.Nx1)
    end
end

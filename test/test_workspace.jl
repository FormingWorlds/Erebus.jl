using Test
using Random
using StaticArrays
using Erebus

@testset "Segregation Workspaces" begin
    @testset "Workspace Submodule Exports" begin
        @test :MetalSegregationWorkspace in names(Erebus.Numerics)
        @test :MagmaSegregationWorkspace in names(Erebus.Numerics)
        @test :MetalSegregationWorkspace in names(Erebus.Simulation)
        @test :MagmaSegregationWorkspace in names(Erebus.Simulation)
        @test :MetalSegregationWorkspace in names(Erebus)
        @test :MagmaSegregationWorkspace in names(Erebus)
    end

    @testset "Workspace Constructors & Dimensions" begin
        Ny, Nx = 16, 20
        ws_metal = MetalSegregationWorkspace(Ny, Nx; track_volatiles=true)
        @test ws_metal.Ny == Ny
        @test ws_metal.Nx == Nx
        @test size(ws_metal.M_fe_cell) == (Ny, Nx)
        @test size(ws_metal.M_rock_markers) == (Ny, Nx)
        @test size(ws_metal.v_seg_cell) == (Ny, Nx)
        @test size(ws_metal.req_flux_x) == (Ny, Nx - 1)
        @test size(ws_metal.req_flux_y) == (Ny - 1, Nx)
        @test size(ws_metal.flux_H_x) == (Ny, Nx - 1)
        @test size(ws_metal.flux_H_y) == (Ny - 1, Nx)

        ws_metal_novol = MetalSegregationWorkspace(Ny, Nx; track_volatiles=false)
        @test size(ws_metal_novol.M_fe_H_cell) == (0, 0)
        @test size(ws_metal_novol.flux_H_x) == (0, 0)

        ws_magma = MagmaSegregationWorkspace(Ny, Nx)
        @test ws_magma.Ny == Ny
        @test ws_magma.Nx == Nx
        @test size(ws_magma.M_melt_cell) == (Ny, Nx)
        @test size(ws_magma.M_rock_markers) == (Ny, Nx)
        @test size(ws_magma.v_seg_cell) == (Ny, Nx)
        @test size(ws_magma.drho_cell) == (Ny, Nx)
        @test size(ws_magma.req_flux_x) == (Ny, Nx - 1)
        @test size(ws_magma.req_flux_y) == (Ny - 1, Nx)
    end

    @testset "Metal Segregation Workspace Equivalence & Reuse" begin
        cfg_core = CoreFormationConfig(;
            percolation_active=true,
            settling_active=true,
            cfl_settling=0.5,
            phi_residual=0.01,
            phi_pack=0.6,
        )
        coords = GridCoordinates(GridConfig(; Nx=20, Ny=20, xsize=140000.0, ysize=140000.0))

        marknum = 200
        rng1 = MersenneTwister(42)
        xm_orig = 20000.0 .+ 60000.0 .* rand(rng1, marknum)
        ym_orig = 20000.0 .+ 60000.0 .* rand(rng1, marknum)
        tm_orig = ones(Int, marknum)
        tkm_orig = fill(1500.0, marknum)
        phim_orig = fill(0.1, marknum)
        Xfe_bulk_orig = fill(0.2, marknum)
        Xfem_orig = fill(0.15, marknum)

        xm1 = copy(xm_orig)
        ym1 = copy(ym_orig)
        tm1 = copy(tm_orig)
        tkm1 = copy(tkm_orig)
        phim1 = copy(phim_orig)
        Xfe_bulk1 = copy(Xfe_bulk_orig)
        Xfem1 = copy(Xfem_orig)

        xm2 = copy(xm_orig)
        ym2 = copy(ym_orig)
        tm2 = copy(tm_orig)
        tkm2 = copy(tkm_orig)
        phim2 = copy(phim_orig)
        Xfe_bulk2 = copy(Xfe_bulk_orig)
        Xfem2 = copy(Xfem_orig)

        dt = 1.0e8

        res_no_ws = apply_metal_segregation!(
            xm1,
            ym1,
            tm1,
            tkm1,
            phim1,
            Xfe_bulk1,
            Xfem1,
            marknum,
            dt,
            cfg_core;
            coords=coords,
            workspace=nothing,
        )

        ws = MetalSegregationWorkspace(coords.Ny, coords.Nx; track_volatiles=false)
        res_ws = apply_metal_segregation!(
            xm2,
            ym2,
            tm2,
            tkm2,
            phim2,
            Xfe_bulk2,
            Xfem2,
            marknum,
            dt,
            cfg_core;
            coords=coords,
            workspace=ws,
        )

        @test res_no_ws.max_v_seg ≈ res_ws.max_v_seg
        @test res_no_ws.n_subcycles == res_ws.n_subcycles
        @test res_no_ws.dt_sub ≈ res_ws.dt_sub
        @test res_no_ws.total_dissipation_energy ≈ res_ws.total_dissipation_energy
        @test Xfe_bulk1 ≈ Xfe_bulk2
        @test Xfem1 ≈ Xfem2

        # Sequential reuse of same workspace with distinct second dataset
        rng2 = MersenneTwister(101)
        xm_b = 25000.0 .+ 50000.0 .* rand(rng2, marknum)
        ym_b = 25000.0 .+ 50000.0 .* rand(rng2, marknum)
        tm_b = ones(Int, marknum)
        tkm_b = fill(1550.0, marknum)
        phim_b = fill(0.12, marknum)
        Xfe_bulk_b1 = fill(0.22, marknum)
        Xfe_bulk_b2 = copy(Xfe_bulk_b1)
        Xfem_b1 = fill(0.18, marknum)
        Xfem_b2 = copy(Xfem_b1)

        res_no_ws_b = apply_metal_segregation!(
            copy(xm_b),
            copy(ym_b),
            copy(tm_b),
            copy(tkm_b),
            copy(phim_b),
            Xfe_bulk_b1,
            Xfem_b1,
            marknum,
            dt,
            cfg_core;
            coords=coords,
            workspace=nothing,
        )
        res_ws_reuse = apply_metal_segregation!(
            copy(xm_b),
            copy(ym_b),
            copy(tm_b),
            copy(tkm_b),
            copy(phim_b),
            Xfe_bulk_b2,
            Xfem_b2,
            marknum,
            dt,
            cfg_core;
            coords=coords,
            workspace=ws,
        )
        @test res_ws_reuse.max_v_seg ≈ res_no_ws_b.max_v_seg
        @test res_ws_reuse.n_subcycles == res_no_ws_b.n_subcycles
        @test res_ws_reuse.dt_sub ≈ res_no_ws_b.dt_sub
        @test res_ws_reuse.total_dissipation_energy ≈ res_no_ws_b.total_dissipation_energy
        @test Xfe_bulk_b1 ≈ Xfe_bulk_b2

        # Oversized workspace mismatch fallback
        xm_over = copy(xm_orig)
        ym_over = copy(ym_orig)
        tm_over = copy(tm_orig)
        tkm_over = copy(tkm_orig)
        phim_over = copy(phim_orig)
        Xfe_bulk_over = copy(Xfe_bulk_orig)
        Xfem_over = copy(Xfem_orig)
        ws_over = MetalSegregationWorkspace(
            coords.Ny + 4, coords.Nx + 4; track_volatiles=false
        )
        res_over = apply_metal_segregation!(
            xm_over,
            ym_over,
            tm_over,
            tkm_over,
            phim_over,
            Xfe_bulk_over,
            Xfem_over,
            marknum,
            dt,
            cfg_core;
            coords=coords,
            workspace=ws_over,
        )
        @test res_over.max_v_seg ≈ res_no_ws.max_v_seg
        @test res_over.n_subcycles == res_no_ws.n_subcycles
        @test res_over.total_dissipation_energy ≈ res_no_ws.total_dissipation_energy

        # Undersized workspace mismatch fallback (ensures bounds check prevents OOB access)
        xm_under = copy(xm_orig)
        ym_under = copy(ym_orig)
        tm_under = copy(tm_orig)
        tkm_under = copy(tkm_orig)
        phim_under = copy(phim_orig)
        Xfe_bulk_under = copy(Xfe_bulk_orig)
        Xfem_under = copy(Xfem_orig)
        ws_under = MetalSegregationWorkspace(
            coords.Ny - 4, coords.Nx - 4; track_volatiles=false
        )
        res_under = apply_metal_segregation!(
            xm_under,
            ym_under,
            tm_under,
            tkm_under,
            phim_under,
            Xfe_bulk_under,
            Xfem_under,
            marknum,
            dt,
            cfg_core;
            coords=coords,
            workspace=ws_under,
        )
        @test res_under.max_v_seg ≈ res_no_ws.max_v_seg
        @test res_under.n_subcycles == res_no_ws.n_subcycles
        @test res_under.total_dissipation_energy ≈ res_no_ws.total_dissipation_energy
        @test Xfe_bulk_under ≈ Xfe_bulk1
    end

    @testset "Metal Segregation with Volatiles Workspace Equivalence & Reuse" begin
        cfg_core = CoreFormationConfig(;
            percolation_active=true,
            settling_active=true,
            cfl_settling=0.5,
            phi_residual=0.01,
            phi_pack=0.6,
        )
        cfg_part = MetalPartitionConfig(; active=true)
        coords = GridCoordinates(GridConfig(; Nx=20, Ny=20, xsize=140000.0, ysize=140000.0))

        marknum = 200
        rng1 = MersenneTwister(77)
        xm_orig = 20000.0 .+ 60000.0 .* rand(rng1, marknum)
        ym_orig = 20000.0 .+ 60000.0 .* rand(rng1, marknum)
        tm_orig = ones(Int, marknum)
        tkm_orig = fill(1550.0, marknum)
        phim_orig = fill(0.12, marknum)
        Xfe_bulk_orig = fill(0.25, marknum)
        Xfem_orig = fill(0.18, marknum)
        Xfe_H_orig = fill(10.0, marknum)
        Xfe_C_orig = fill(500.0, marknum)
        Xfe_N_orig = fill(20.0, marknum)
        Xfe_S_orig = fill(2000.0, marknum)

        xm1 = copy(xm_orig)
        ym1 = copy(ym_orig)
        tm1 = copy(tm_orig)
        tkm1 = copy(tkm_orig)
        phim1 = copy(phim_orig)
        Xfe_bulk1 = copy(Xfe_bulk_orig)
        Xfem1 = copy(Xfem_orig)
        Xfe_H_1 = copy(Xfe_H_orig)
        Xfe_C_1 = copy(Xfe_C_orig)
        Xfe_N_1 = copy(Xfe_N_orig)
        Xfe_S_1 = copy(Xfe_S_orig)

        xm2 = copy(xm_orig)
        ym2 = copy(ym_orig)
        tm2 = copy(tm_orig)
        tkm2 = copy(tkm_orig)
        phim2 = copy(phim_orig)
        Xfe_bulk2 = copy(Xfe_bulk_orig)
        Xfem2 = copy(Xfem_orig)
        Xfe_H_2 = copy(Xfe_H_orig)
        Xfe_C_2 = copy(Xfe_C_orig)
        Xfe_N_2 = copy(Xfe_N_orig)
        Xfe_S_2 = copy(Xfe_S_orig)

        dt = 1.0e8

        res_no_ws = apply_metal_segregation!(
            xm1,
            ym1,
            tm1,
            tkm1,
            phim1,
            Xfe_bulk1,
            Xfem1,
            marknum,
            dt,
            cfg_core;
            coords=coords,
            workspace=nothing,
            cfg_partition=cfg_part,
            Xfe_H_m=Xfe_H_1,
            Xfe_C_m=Xfe_C_1,
            Xfe_N_m=Xfe_N_1,
            Xfe_S_m=Xfe_S_1,
        )

        ws = MetalSegregationWorkspace(coords.Ny, coords.Nx; track_volatiles=true)
        res_ws = apply_metal_segregation!(
            xm2,
            ym2,
            tm2,
            tkm2,
            phim2,
            Xfe_bulk2,
            Xfem2,
            marknum,
            dt,
            cfg_core;
            coords=coords,
            workspace=ws,
            cfg_partition=cfg_part,
            Xfe_H_m=Xfe_H_2,
            Xfe_C_m=Xfe_C_2,
            Xfe_N_m=Xfe_N_2,
            Xfe_S_m=Xfe_S_2,
        )

        @test res_no_ws.max_v_seg ≈ res_ws.max_v_seg
        @test res_no_ws.n_subcycles == res_ws.n_subcycles
        @test res_no_ws.dt_sub ≈ res_ws.dt_sub
        @test res_no_ws.total_dissipation_energy ≈ res_ws.total_dissipation_energy
        @test Xfe_bulk1 ≈ Xfe_bulk2
        @test Xfe_H_1 ≈ Xfe_H_2
        @test Xfe_C_1 ≈ Xfe_C_2
        @test Xfe_N_1 ≈ Xfe_N_2
        @test Xfe_S_1 ≈ Xfe_S_2

        # Sequential reuse of same volatile workspace with distinct second dataset
        rng3 = MersenneTwister(202)
        xm_v2 = 22000.0 .+ 55000.0 .* rand(rng3, marknum)
        ym_v2 = 22000.0 .+ 55000.0 .* rand(rng3, marknum)
        tm_v2 = ones(Int, marknum)
        tkm_v2 = fill(1580.0, marknum)
        phim_v2 = fill(0.15, marknum)
        Xfe_bulk_v2a = fill(0.28, marknum)
        Xfe_bulk_v2b = copy(Xfe_bulk_v2a)
        Xfem_v2a = fill(0.20, marknum)
        Xfem_v2b = copy(Xfem_v2a)
        Xfe_H_v2a = fill(15.0, marknum)
        Xfe_H_v2b = copy(Xfe_H_v2a)
        Xfe_C_v2a = fill(750.0, marknum)
        Xfe_C_v2b = copy(Xfe_C_v2a)
        Xfe_N_v2a = fill(35.0, marknum)
        Xfe_N_v2b = copy(Xfe_N_v2a)
        Xfe_S_v2a = fill(3000.0, marknum)
        Xfe_S_v2b = copy(Xfe_S_v2a)

        res_no_ws_v2 = apply_metal_segregation!(
            copy(xm_v2),
            copy(ym_v2),
            copy(tm_v2),
            copy(tkm_v2),
            copy(phim_v2),
            Xfe_bulk_v2a,
            Xfem_v2a,
            marknum,
            dt,
            cfg_core;
            coords=coords,
            workspace=nothing,
            cfg_partition=cfg_part,
            Xfe_H_m=Xfe_H_v2a,
            Xfe_C_m=Xfe_C_v2a,
            Xfe_N_m=Xfe_N_v2a,
            Xfe_S_m=Xfe_S_v2a,
        )
        res_ws_reuse_v2 = apply_metal_segregation!(
            copy(xm_v2),
            copy(ym_v2),
            copy(tm_v2),
            copy(tkm_v2),
            copy(phim_v2),
            Xfe_bulk_v2b,
            Xfem_v2b,
            marknum,
            dt,
            cfg_core;
            coords=coords,
            workspace=ws,
            cfg_partition=cfg_part,
            Xfe_H_m=Xfe_H_v2b,
            Xfe_C_m=Xfe_C_v2b,
            Xfe_N_m=Xfe_N_v2b,
            Xfe_S_m=Xfe_S_v2b,
        )
        @test res_ws_reuse_v2.max_v_seg ≈ res_no_ws_v2.max_v_seg
        @test res_ws_reuse_v2.n_subcycles == res_no_ws_v2.n_subcycles
        @test res_ws_reuse_v2.total_dissipation_energy ≈
            res_no_ws_v2.total_dissipation_energy
        @test Xfe_bulk_v2a ≈ Xfe_bulk_v2b
        @test Xfe_H_v2a ≈ Xfe_H_v2b
        @test Xfe_C_v2a ≈ Xfe_C_v2b
        @test Xfe_N_v2a ≈ Xfe_N_v2b
        @test Xfe_S_v2a ≈ Xfe_S_v2b

        # Flag mismatch test 1: track_volatiles active but workspace without volatile buffers
        ws_no_vol = MetalSegregationWorkspace(coords.Ny, coords.Nx; track_volatiles=false)
        xm_mism = copy(xm_orig)
        ym_mism = copy(ym_orig)
        tm_mism = copy(tm_orig)
        tkm_mism = copy(tkm_orig)
        phim_mism = copy(phim_orig)
        Xfe_bulk_mism = copy(Xfe_bulk_orig)
        Xfem_mism = copy(Xfem_orig)
        Xfe_H_mism = copy(Xfe_H_orig)
        Xfe_C_mism = copy(Xfe_C_orig)
        Xfe_N_mism = copy(Xfe_N_orig)
        Xfe_S_mism = copy(Xfe_S_orig)
        res_mism = apply_metal_segregation!(
            xm_mism,
            ym_mism,
            tm_mism,
            tkm_mism,
            phim_mism,
            Xfe_bulk_mism,
            Xfem_mism,
            marknum,
            dt,
            cfg_core;
            coords=coords,
            workspace=ws_no_vol,
            cfg_partition=cfg_part,
            Xfe_H_m=Xfe_H_mism,
            Xfe_C_m=Xfe_C_mism,
            Xfe_N_m=Xfe_N_mism,
            Xfe_S_m=Xfe_S_mism,
        )
        @test res_mism.max_v_seg ≈ res_no_ws.max_v_seg
        @test res_mism.n_subcycles == res_no_ws.n_subcycles
        @test res_mism.total_dissipation_energy ≈ res_no_ws.total_dissipation_energy
        @test Xfe_bulk_mism ≈ Xfe_bulk1
        @test Xfe_H_mism ≈ Xfe_H_1

        # Flag mismatch test 2: track_volatiles false but workspace with volatile buffers
        xm_novol_a = copy(xm_orig)
        ym_novol_a = copy(ym_orig)
        tm_novol_a = copy(tm_orig)
        tkm_novol_a = copy(tkm_orig)
        phim_novol_a = copy(phim_orig)
        Xfe_bulk_novol_a = copy(Xfe_bulk_orig)
        Xfem_novol_a = copy(Xfem_orig)
        xm_novol_b = copy(xm_orig)
        ym_novol_b = copy(ym_orig)
        tm_novol_b = copy(tm_orig)
        tkm_novol_b = copy(tkm_orig)
        phim_novol_b = copy(phim_orig)
        Xfe_bulk_novol_b = copy(Xfe_bulk_orig)
        Xfem_novol_b = copy(Xfem_orig)
        res_novol_ref = apply_metal_segregation!(
            xm_novol_a,
            ym_novol_a,
            tm_novol_a,
            tkm_novol_a,
            phim_novol_a,
            Xfe_bulk_novol_a,
            Xfem_novol_a,
            marknum,
            dt,
            cfg_core;
            coords=coords,
            workspace=nothing,
        )
        res_novol_ws = apply_metal_segregation!(
            xm_novol_b,
            ym_novol_b,
            tm_novol_b,
            tkm_novol_b,
            phim_novol_b,
            Xfe_bulk_novol_b,
            Xfem_novol_b,
            marknum,
            dt,
            cfg_core;
            coords=coords,
            workspace=ws,
        )
        @test res_novol_ws.max_v_seg ≈ res_novol_ref.max_v_seg
        @test res_novol_ws.n_subcycles == res_novol_ref.n_subcycles
        @test res_novol_ws.total_dissipation_energy ≈ res_novol_ref.total_dissipation_energy
        @test Xfe_bulk_novol_b ≈ Xfe_bulk_novol_a
    end

    @testset "Magma Segregation Workspace Equivalence & Reuse" begin
        cfg_magma = MagmaTransportConfig(;
            active=true, cfl_melt=0.5, phi_residual=0.01, phi_pack=0.6
        )
        coords = GridCoordinates(GridConfig(; Nx=20, Ny=20, xsize=140000.0, ysize=140000.0))

        marknum = 200
        rng1 = MersenneTwister(99)
        xm_orig = 20000.0 .+ 60000.0 .* rand(rng1, marknum)
        ym_orig = 20000.0 .+ 60000.0 .* rand(rng1, marknum)
        tm_orig = ones(Int, marknum)
        tkm_orig = fill(1600.0, marknum)
        Fm_orig = fill(0.25, marknum)

        xm1 = copy(xm_orig)
        ym1 = copy(ym_orig)
        tm1 = copy(tm_orig)
        tkm1 = copy(tkm_orig)
        Fm1 = copy(Fm_orig)
        xm2 = copy(xm_orig)
        ym2 = copy(ym_orig)
        tm2 = copy(tm_orig)
        tkm2 = copy(tkm_orig)
        Fm2 = copy(Fm_orig)

        dt = 1.0e8

        res_no_ws = apply_silicate_melt_segregation!(
            xm1,
            ym1,
            tm1,
            tkm1,
            Fm1,
            marknum,
            dt,
            cfg_magma;
            coords=coords,
            workspace=nothing,
        )

        ws = MagmaSegregationWorkspace(coords.Ny, coords.Nx)
        res_ws = apply_silicate_melt_segregation!(
            xm2, ym2, tm2, tkm2, Fm2, marknum, dt, cfg_magma; coords=coords, workspace=ws
        )

        @test res_no_ws.max_v_seg ≈ res_ws.max_v_seg
        @test res_no_ws.n_subcycles == res_ws.n_subcycles
        @test res_no_ws.dt_sub ≈ res_ws.dt_sub
        @test res_no_ws.total_dissipation_energy ≈ res_ws.total_dissipation_energy
        @test res_no_ws.total_crystallized_mass ≈ res_ws.total_crystallized_mass
        @test Fm1 ≈ Fm2

        # Sequential reuse of same magma workspace on distinct second dataset
        rng_m2 = MersenneTwister(303)
        xm_m2 = 25000.0 .+ 50000.0 .* rand(rng_m2, marknum)
        ym_m2 = 25000.0 .+ 50000.0 .* rand(rng_m2, marknum)
        tm_m2 = ones(Int, marknum)
        tkm_m2 = fill(1650.0, marknum)
        Fm_m2a = fill(0.30, marknum)
        Fm_m2b = copy(Fm_m2a)

        res_no_ws_m2 = apply_silicate_melt_segregation!(
            copy(xm_m2),
            copy(ym_m2),
            copy(tm_m2),
            copy(tkm_m2),
            Fm_m2a,
            marknum,
            dt,
            cfg_magma;
            coords=coords,
            workspace=nothing,
        )
        res_ws_reuse = apply_silicate_melt_segregation!(
            copy(xm_m2),
            copy(ym_m2),
            copy(tm_m2),
            copy(tkm_m2),
            Fm_m2b,
            marknum,
            dt,
            cfg_magma;
            coords=coords,
            workspace=ws,
        )
        @test res_ws_reuse.max_v_seg ≈ res_no_ws_m2.max_v_seg
        @test res_ws_reuse.n_subcycles == res_no_ws_m2.n_subcycles
        @test res_ws_reuse.dt_sub ≈ res_no_ws_m2.dt_sub
        @test res_ws_reuse.total_dissipation_energy ≈ res_no_ws_m2.total_dissipation_energy
        @test res_ws_reuse.total_crystallized_mass ≈ res_no_ws_m2.total_crystallized_mass
        @test Fm_m2b ≈ Fm_m2a

        # Oversized workspace mismatch fallback tolerance
        xm_over = copy(xm_orig)
        ym_over = copy(ym_orig)
        tm_over = copy(tm_orig)
        tkm_over = copy(tkm_orig)
        Fm_over = copy(Fm_orig)
        ws_over = MagmaSegregationWorkspace(coords.Ny + 2, coords.Nx + 2)
        res_over = apply_silicate_melt_segregation!(
            xm_over,
            ym_over,
            tm_over,
            tkm_over,
            Fm_over,
            marknum,
            dt,
            cfg_magma;
            coords=coords,
            workspace=ws_over,
        )
        @test res_over.max_v_seg ≈ res_no_ws.max_v_seg
        @test res_over.n_subcycles == res_no_ws.n_subcycles
        @test res_over.total_dissipation_energy ≈ res_no_ws.total_dissipation_energy
        @test res_over.total_crystallized_mass ≈ res_no_ws.total_crystallized_mass
        @test Fm_over ≈ Fm1

        # Undersized workspace mismatch fallback tolerance
        xm_under = copy(xm_orig)
        ym_under = copy(ym_orig)
        tm_under = copy(tm_orig)
        tkm_under = copy(tkm_orig)
        Fm_under = copy(Fm_orig)
        ws_under = MagmaSegregationWorkspace(coords.Ny - 2, coords.Nx - 2)
        res_under = apply_silicate_melt_segregation!(
            xm_under,
            ym_under,
            tm_under,
            tkm_under,
            Fm_under,
            marknum,
            dt,
            cfg_magma;
            coords=coords,
            workspace=ws_under,
        )
        @test res_under.max_v_seg ≈ res_no_ws.max_v_seg
        @test res_under.n_subcycles == res_no_ws.n_subcycles
        @test res_under.total_dissipation_energy ≈ res_no_ws.total_dissipation_energy
        @test res_under.total_crystallized_mass ≈ res_no_ws.total_crystallized_mass
        @test Fm_under ≈ Fm1
    end
end

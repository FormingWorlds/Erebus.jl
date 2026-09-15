using Erebus
using ExtendableSparse
using LinearAlgebra
using SparseArrays
using Test

@testset "Stokes-Darcy Analytical Darcy Elimination" begin
    cfg = load_config("configs/test_quick.toml")
    coords = GridCoordinates(
        cfg.grid.Nx, cfg.grid.Ny; xsize=cfg.grid.xsize, ysize=cfg.grid.ysize
    )
    Nx1 = coords.Nx1
    Ny1 = coords.Ny1
    Nx_val = coords.Nx
    Ny_val = coords.Ny
    dx_val = coords.dx
    dy_val = coords.dy

    # Base material fields
    ETA = fill(1.0e19, Ny1, Nx1)
    ETAP = fill(1.0e19, Ny1, Nx1)
    GGG = fill(1.0e10, Ny1, Nx1)
    GGGP = fill(1.0e10, Ny1, Nx1)
    SXY0 = zeros(Ny1, Nx1)
    SXX0 = zeros(Ny1, Nx1)
    RHOX = fill(3000.0, Ny1, Nx1)
    RHOY = fill(3000.0, Ny1, Nx1)
    RHOFX = fill(1000.0, Ny1, Nx1)
    RHOFY = fill(1000.0, Ny1, Nx1)
    RX = fill(1.0e8, Ny1, Nx1)
    RY = fill(1.0e8, Ny1, Nx1)
    ETAPHI = fill(1.0e19, Ny1, Nx1)
    BETAPHI = fill(1.0e-11, Ny1, Nx1)
    PHI = fill(0.1, Ny1, Nx1)
    gx = zeros(Ny1, Nx1)
    gy = fill(1.0, Ny1, Nx1)
    pr0 = fill(1.0e6, Ny1, Nx1)
    pf0 = fill(1.0e6, Ny1, Nx1)
    DMP = zeros(Ny1, Nx1)
    dt = 1.0e10

    @testset "Mathematical equivalence of 6-var and 4-var solutions" begin
        # 6-variable assembly and solve
        R6 = zeros(Ny1 * Nx1 * 6)
        L6 = assemble_hydromechanical_lse!(
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
            R6;
            coords=coords,
        )
        S6 = L6 \ R6

        vx6 = zeros(Ny1, Nx1)
        vy6 = zeros(Ny1, Nx1)
        pr6 = zeros(Ny1, Nx1)
        qxD6 = zeros(Ny1, Nx1)
        qyD6 = zeros(Ny1, Nx1)
        pf6 = zeros(Ny1, Nx1)
        process_hydromechanical_solution!(S6, vx6, vy6, pr6, qxD6, qyD6, pf6; coords=coords)

        # 4-variable assembly and solve
        R4 = zeros(Ny1 * Nx1 * 4)
        L4 = assemble_hydromechanical_4var_lse!(
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
            R4;
            coords=coords,
        )
        S4 = L4 \ R4

        vx4 = zeros(Ny1, Nx1)
        vy4 = zeros(Ny1, Nx1)
        pr4 = zeros(Ny1, Nx1)
        pf4 = zeros(Ny1, Nx1)
        process_hydromechanical_4var_solution!(S4, vx4, vy4, pr4, pf4; coords=coords)

        qxD4 = zeros(Ny1, Nx1)
        qyD4 = zeros(Ny1, Nx1)
        reconstruct_darcy_fluxes!(qxD4, qyD4, pf4, RHOFX, RHOFY, RX, RY, gx, gy, coords)

        # Assert mathematical equivalence between 6-var and 4-var
        rel_diff_vx = norm(vx4 - vx6) / norm(vx6)
        rel_diff_vy = norm(vy4 - vy6) / norm(vy6)
        rel_diff_pr = norm(pr4 - pr6) / norm(pr6)
        rel_diff_pf = norm(pf4 - pf6) / norm(pf6)
        rel_diff_qxD = norm(qxD4 - qxD6) / (norm(qxD6) + 1.0e-30)
        rel_diff_qyD = norm(qyD4 - qyD6) / norm(qyD6)

        @test rel_diff_vx < 1.0e-12
        @test rel_diff_vy < 1.0e-12
        @test rel_diff_pr < 1.0e-8
        @test rel_diff_pf < 1.0e-8
        @test rel_diff_qxD < 1.0e-10
        @test rel_diff_qyD < 1.0e-10
    end

    @testset "Darcy flux reconstruction properties" begin
        R4 = zeros(Ny1 * Nx1 * 4)
        L4 = assemble_hydromechanical_4var_lse!(
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
            R4;
            coords=coords,
        )
        S4 = L4 \ R4
        vx4 = zeros(Ny1, Nx1)
        vy4 = zeros(Ny1, Nx1)
        pr4 = zeros(Ny1, Nx1)
        pf4 = zeros(Ny1, Nx1)
        process_hydromechanical_4var_solution!(S4, vx4, vy4, pr4, pf4; coords=coords)

        qxD4 = zeros(Ny1, Nx1)
        qyD4 = zeros(Ny1, Nx1)
        reconstruct_darcy_fluxes!(qxD4, qyD4, pf4, RHOFX, RHOFY, RX, RY, gx, gy, coords)

        # Boundary ghost relations
        bcftop = Erebus.bcftop
        bcfbottom = Erebus.bcfbottom
        bcfleft = Erebus.bcfleft
        bcfright = Erebus.bcfright

        top_check = maximum(abs, qxD4[1, 2:Nx_val] .+ bcftop .* qxD4[2, 2:Nx_val])
        bottom_check = maximum(
            abs, qxD4[Ny1, 2:Nx_val] .+ bcfbottom .* qxD4[Ny_val, 2:Nx_val]
        )
        left_check = maximum(abs, qyD4[2:Ny_val, 1] .+ bcfleft .* qyD4[2:Ny_val, 2])
        right_check = maximum(
            abs, qyD4[2:Ny_val, Nx1] .+ bcfright .* qyD4[2:Ny_val, Nx_val]
        )

        @test top_check < 1.0e-14
        @test bottom_check < 1.0e-14
        @test left_check < 1.0e-14
        @test right_check < 1.0e-14
    end

    @testset "Workspace allocation and reuse" begin
        ws4 = HydromechanicalLSEWorkspace(coords; dof_per_node=4)
        @test ws4.Ny1 == Ny1
        @test ws4.Nx1 == Nx1
        @test ws4.dof_per_node == 4
        @test size(ws4.L, 1) == Ny1 * Nx1 * 4
        @test !ws4.is_initialized

        R4 = zeros(Ny1 * Nx1 * 4)
        L4 = assemble_hydromechanical_4var_lse!(
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
            R4;
            coords=coords,
            L=ws4.L,
            workspace=ws4,
        )
        @test ws4.is_initialized
        @test size(L4) == (Ny1 * Nx1 * 4, Ny1 * Nx1 * 4)
    end

    @testset "Hydrofracture permeability preservation in 4-variable system" begin
        pr_test = fill(2.0e6, Ny1, Nx1)
        pf_test = fill(2.5e6, Ny1, Nx1) # pore pressure exceeds total pressure
        TEN_test = fill(1.0e5, Ny1, Nx1)
        KX_test = fill(1.0e-15, Ny1, Nx1)
        KY_test = fill(1.0e-15, Ny1, Nx1)

        R6 = zeros(Ny1 * Nx1 * 6)
        L6 = assemble_hydromechanical_lse!(
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
            R6;
            coords=coords,
            hydrofracture=true,
            pr=pr_test,
            pf=pf_test,
            TEN=TEN_test,
            KX=KX_test,
            KY=KY_test,
        )
        S6 = L6 \ R6
        vx6 = zeros(Ny1, Nx1)
        vy6 = zeros(Ny1, Nx1)
        pr6 = zeros(Ny1, Nx1)
        qxD6 = zeros(Ny1, Nx1)
        qyD6 = zeros(Ny1, Nx1)
        pf6 = zeros(Ny1, Nx1)
        process_hydromechanical_solution!(S6, vx6, vy6, pr6, qxD6, qyD6, pf6; coords=coords)

        R4 = zeros(Ny1 * Nx1 * 4)
        L4 = assemble_hydromechanical_4var_lse!(
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
            R4;
            coords=coords,
            hydrofracture=true,
            pr=pr_test,
            pf=pf_test,
            TEN=TEN_test,
            KX=KX_test,
            KY=KY_test,
        )
        S4 = L4 \ R4
        vx4 = zeros(Ny1, Nx1)
        vy4 = zeros(Ny1, Nx1)
        pr4 = zeros(Ny1, Nx1)
        pf4 = zeros(Ny1, Nx1)
        process_hydromechanical_4var_solution!(S4, vx4, vy4, pr4, pf4; coords=coords)

        qxD4 = zeros(Ny1, Nx1)
        qyD4 = zeros(Ny1, Nx1)
        reconstruct_darcy_fluxes!(
            qxD4,
            qyD4,
            pf4,
            RHOFX,
            RHOFY,
            RX,
            RY,
            gx,
            gy,
            coords;
            hydrofracture=true,
            pr=pr_test,
            pf_eff=pf_test,
            TEN=TEN_test,
            KX=KX_test,
            KY=KY_test,
        )

        @test norm(vx4 - vx6) / norm(vx6) < 1.0e-9
        @test norm(vy4 - vy6) / norm(vy6) < 1.0e-9
        @test norm(pr4 - pr6) / norm(pr6) < 1.0e-8
        @test norm(pf4 - pf6) / norm(pf6) < 1.0e-5
        @test norm(qxD4 - qxD6) / (norm(qxD6) + 1.0e-30) < 1.0e-3
        @test norm(qyD4 - qyD6) / norm(qyD6) < 1.0e-9
    end

    @testset "Heterogeneous Viscosity 6-var vs 4-var Equivalence" begin
        # Heterogeneous viscosity and shear modulus to verify shear coupling
        ETA_het = [
            1.0e19 * (1.0 + 0.5 * sin(2π * i / Ny1) * cos(2π * j / Nx1)) for
            i in 1:Ny1, j in 1:Nx1
        ]
        ETAP_het = [
            1.0e19 * (1.0 + 0.5 * cos(2π * i / Ny1) * sin(2π * j / Nx1)) for
            i in 1:Ny1, j in 1:Nx1
        ]
        GGG_het = [
            1.0e10 * (1.0 + 0.3 * sin(2π * i / Ny1) * sin(2π * j / Nx1)) for
            i in 1:Ny1, j in 1:Nx1
        ]
        GGGP_het = [
            1.0e10 * (1.0 + 0.3 * cos(2π * i / Ny1) * cos(2π * j / Nx1)) for
            i in 1:Ny1, j in 1:Nx1
        ]

        R6 = zeros(Ny1 * Nx1 * 6)
        L6 = assemble_hydromechanical_lse!(
            ETA_het,
            ETAP_het,
            GGG_het,
            GGGP_het,
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
            R6;
            coords=coords,
        )
        S6 = L6 \ R6
        vx6 = zeros(Ny1, Nx1)
        vy6 = zeros(Ny1, Nx1)
        pr6 = zeros(Ny1, Nx1)
        qxD6 = zeros(Ny1, Nx1)
        qyD6 = zeros(Ny1, Nx1)
        pf6 = zeros(Ny1, Nx1)
        process_hydromechanical_solution!(S6, vx6, vy6, pr6, qxD6, qyD6, pf6; coords=coords)

        R4 = zeros(Ny1 * Nx1 * 4)
        L4 = assemble_hydromechanical_4var_lse!(
            ETA_het,
            ETAP_het,
            GGG_het,
            GGGP_het,
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
            R4;
            coords=coords,
        )
        S4 = L4 \ R4
        vx4 = zeros(Ny1, Nx1)
        vy4 = zeros(Ny1, Nx1)
        pr4 = zeros(Ny1, Nx1)
        pf4 = zeros(Ny1, Nx1)
        process_hydromechanical_4var_solution!(S4, vx4, vy4, pr4, pf4; coords=coords)

        qxD4 = zeros(Ny1, Nx1)
        qyD4 = zeros(Ny1, Nx1)
        reconstruct_darcy_fluxes!(qxD4, qyD4, pf4, RHOFX, RHOFY, RX, RY, gx, gy, coords)

        @test norm(vx4 - vx6) / norm(vx6) < 1.0e-9
        @test norm(vy4 - vy6) / norm(vy6) < 1.0e-9
        @test norm(pr4 - pr6) / norm(pr6) < 1.0e-8
        @test norm(pf4 - pf6) / norm(pf6) < 1.0e-8
        @test norm(qxD4 - qxD6) / (norm(qxD6) + 1.0e-30) < 1.0e-10
        @test norm(qyD4 - qyD6) / norm(qyD6) < 1.0e-10
    end
end

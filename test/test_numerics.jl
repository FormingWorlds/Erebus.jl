@testset "Numerics" begin
    @testset "compute_gravity_solution!()" begin
        SP = zeros(Float64, Nx1*Ny1)
        RP = zeros(Float64, Nx1*Ny1)
        FI = zeros(Float64, Ny1, Nx1)
        gx = zeros(Float64, Ny1, Nx1)
        gy = zeros(Float64, Ny1, Nx1)
        LP_ver = zeros(Nx1*Ny1, Nx1*Ny1)
        RP_ver = zeros(Float64, Nx1*Ny1)
        FI_ver = zeros(Float64, Ny1, Nx1)
        gx_ver = zeros(Float64, Ny1, Nx1)
        gy_ver = zeros(Float64, Ny1, Nx1)
        # simulate density field RHO
        RHO = rand(rgen, Ny1, Nx1) * 7e3
        # compute gravity solution
        Erebus.compute_gravity_solution!(SP, RP, RHO, FI, gx, gy)
        # verification, from HTM-planetary.m, line 680ff
        for j in 1:1:Nx1
            for i in 1:1:Ny1
                # Define global index in algebraic space
                gk=(j-1)*Ny1+i
                # Distance from the model centre
                rnode=((xp[j]-xsize/2)^2+(yp[i]-ysize/2)^2)^0.5
                # External points
                if rnode>xsize/2 || i==1 || i==Ny1 || j==1 || j==Nx1
                    # Boundary Condition
                    # PHI=0
                    LP_ver[gk, gk]=1 # Left part
                    RP_ver[gk]=0 # Right part
                else
                    # Internal points: Temperature eq.
                    # d2PHI/dx^2+d2PHI/dy^2=2/3*4*G*pi*RHO
                    #          PHI2
                    #           |
                    #           |
                    #  PHI1----PHI3----PHI5
                    #           |
                    #           |
                    #          PHI4
                    #
                    # Density gradients
                    dRHOdx=(RHO[i, j + 1]-RHO[i, j - 1])/2/dx
                    dRHOdy=(RHO[i + 1, j]-RHO[i - 1, j])/2/dy
                    # Left part
                    LP_ver[gk, gk - Ny1]=1/dx^2 # PHI1
                    LP_ver[gk, gk - 1]=1/dy^2 # PHI2
                    LP_ver[gk, gk]=-2/dx^2-2/dy^2 # PHI3
                    LP_ver[gk, gk + 1]=1/dy^2 # PHI4
                    LP_ver[gk, gk + Ny1]=1/dx^2 # PHI5
                    # Right part
                    RP_ver[gk]=2/3*4*G*pi*RHO[i, j]
                end
            end
        end
        # Solving matrixes
        SP_ver=LP_ver\RP_ver # Obtaining algebraic vector of solutions SP[]
        # Reload solutions SP[] to geometrical array PHI[]
        # Going through all grid points
        for j in 1:1:Nx1
            for i in 1:1:Ny1
                # Compute global index
                gk=(j-1)*Ny1+i
                # Reload solution
                FI_ver[i, j]=SP_ver[gk]
            end
        end
        # Compute gravity acceleration
        # gx
        for j in 1:1:Nx
            for i in 1:1:Ny1
                # gx=-dPHI/dx
                gx_ver[i, j]=-(FI_ver[i, j + 1]-FI_ver[i, j])/dx
            end
        end
        # gy
        for j in 1:1:Nx1
            for i in 1:1:Ny
                # gy=-dPHI/dy
                gy_ver[i, j]=-(FI_ver[i + 1, j]-FI_ver[i, j])/dy
            end
        end
        # test
        for j in 1:1:Nx, i in 1:1:Ny
            @test FI[i, j] ≈ FI_ver[i, j] rtol=1e-9
            @test gx[i, j] ≈ gx_ver[i, j] rtol=1e-9
            @test gy[i, j] ≈ gy_ver[i, j] rtol=1e-9
        end
    end # testset "compute_gravity_solution!()"

    @testset "assemble_gravitational_lse()!" begin
        RP = rand(rgen, Nx1*Ny1)
        RP_ver = deepcopy(RP)
        LP_ver = zeros(Nx1*Ny1, Nx1*Ny1)
        # simulate density field RHO
        RHO = rand(rgen, Ny1, Nx1) * 7e3
        LP = Erebus.assemble_gravitational_lse!(RHO, RP)
        # verification, from HTM-planetary.m, line 680ff
        for j in 1:1:Nx1
            for i in 1:1:Ny1
                # Define global index in algebraic space
                gk=(j-1)*Ny1+i
                # Distance from the model centre
                rnode=((xp[j]-xsize/2)^2+(yp[i]-ysize/2)^2)^0.5
                # External points
                if rnode>xsize/2 || i==1 || i==Ny1 || j==1 || j==Nx1
                    # Boundary Condition
                    # PHI=0
                    LP_ver[gk, gk]=1 # Left part
                    RP_ver[gk]=0 # Right part
                else
                    # Internal points: Temperature eq.
                    # d2PHI/dx^2+d2PHI/dy^2=2/3*4*G*pi*RHO
                    #          PHI2
                    #           |
                    #           |
                    #  PHI1----PHI3----PHI5
                    #           |
                    #           |
                    #          PHI4
                    #
                    # Density gradients
                    dRHOdx=(RHO[i, j + 1]-RHO[i, j - 1])/2/dx
                    dRHOdy=(RHO[i + 1, j]-RHO[i - 1, j])/2/dy
                    # Left part
                    LP_ver[gk, gk - Ny1]=1/dx^2 # PHI1
                    LP_ver[gk, gk - 1]=1/dy^2 # PHI2
                    LP_ver[gk, gk]=-2/dx^2-2/dy^2 # PHI3
                    LP_ver[gk, gk + 1]=1/dy^2 # PHI4
                    LP_ver[gk, gk + Ny1]=1/dx^2 # PHI5
                    # Right part
                    RP_ver[gk]=2/3*4*G*pi*RHO[i, j]
                end
            end
        end
        # test
        # for j=1:1:Nx1*Ny1, i=1:1:Ny1*Nx1
        #     @test LP[i, j] ≈ LP_ver[i, j] rtol=1e-12
        #     @test RP[i] ≈ RP_ver[i] rtol=1e-12
        # end
        @test LP ≈ LP_ver rtol=1e-12
        @test RP ≈ RP_ver rtol=1e-12
    end # testset "assemble_gravitational_lse!()"

    @testset "process_gravitational_solution" begin
        SP = rand(rgen, Nx1*Ny1)
        FI = zeros(Float64, Ny1, Nx1)
        gx = zeros(Float64, Ny1, Nx1)
        gy = zeros(Float64, Ny1, Nx1)
        FI_ver = zeros(Float64, Ny1, Nx1)
        gx_ver = zeros(Float64, Ny1, Nx1)
        gy_ver = zeros(Float64, Ny1, Nx1)
        Erebus.process_gravitational_solution!(SP, FI, gx, gy)
        # verification, from HTM-planetary.m, line 680ff
        for j in 1:1:Nx1
            for i in 1:1:Ny1
                # Compute global index
                gk=(j-1)*Ny1+i
                # Reload solution
                FI_ver[i, j]=SP[gk]
            end
        end
        # Compute gravity acceleration
        # gx
        for j in 1:1:Nx
            for i in 1:1:Ny1
                # gx=-dPHI/dx
                gx_ver[i, j]=-(FI_ver[i, j + 1]-FI_ver[i, j])/dx
            end
        end
        # gy
        for j in 1:1:Nx1
            for i in 1:1:Ny
                # gy=-dPHI/dy
                gy_ver[i, j]=-(FI_ver[i + 1, j]-FI_ver[i, j])/dy
            end
        end
        # test
        for j in 1:1:Nx, i in 1:1:Ny
            @test FI[i, j] ≈ FI_ver[i, j] rtol=1e-9
            @test gx[i, j] ≈ gx_ver[i, j] rtol=1e-9
            @test gy[i, j] ≈ gy_ver[i, j] rtol=1e-9
        end
    end # testset "process_gravitational_solution"

    @testset "recompute_bulk_viscosity!()" begin
        ETAP = zeros(Ny1, Nx1)
        ETAPHI = zeros(Ny1, Nx1)
        ETAP_ver = zeros(Ny1, Nx1)
        ETAPHI_ver = zeros(Ny1, Nx1)
        # simulate data
        ETA = rand(rgen, Ny, Nx)
        PHI = rand(rgen, Ny1, Nx1)
        # compute bulk viscosity
        Erebus.recompute_bulk_viscosity!(ETA, ETAP, ETAPHI, PHI, etaphikoef)
        # verification, from HTM-planetary.m, line 771ff
        for i in 2:1:Ny
            for j in 2:1:Nx
                ETAP_ver[i, j]=1/(
                    (1/ETA[i - 1, j - 1]+1/ETA[i, j - 1]+1/ETA[i - 1, j]+1/ETA[i, j])/4
                )
                ETAPHI_ver[i, j]=etaphikoef*ETAP_ver[i, j]/PHI[i, j]
            end
        end
        # test
        for j in 1:1:Nx, i in 1:1:Ny
            @test ETAP[i, j] ≈ ETAP_ver[i, j] rtol=1e-9
            @test ETAPHI[i, j] ≈ ETAPHI_ver[i, j] rtol=1e-9
        end
    end # testset "recompute_bulk_viscosity!()"

    @testset "get_viscosities_stresses_density_gradients()" begin
        dt = dt_longest
        # simulate data
        ETA = rand(rgen, Ny, Nx)
        ETAP = rand(rgen, Ny1, Nx1)
        GGG = rand(rgen, Ny, Nx)
        GGGP = rand(rgen, Ny1, Nx1)
        SXY0 = rand(rgen, Ny, Nx)
        SXX0 = rand(rgen, Ny1, Nx1)
        RHOX = rand(rgen, Ny1, Nx1)
        RHOY = rand(rgen, Ny1, Nx1)
        ETAcomp = zeros(Ny, Nx)
        ETAPcomp = zeros(Ny1, Nx1)
        SXYcomp = zeros(Ny, Nx)
        SXXcomp = zeros(Ny1, Nx1)
        SYYcomp = zeros(Ny1, Nx1)
        dRHOXdx = zeros(Ny1, Nx1)
        dRHOXdy = zeros(Ny1, Nx1)
        dRHOYdx = zeros(Ny1, Nx1)
        dRHOYdy = zeros(Ny1, Nx1)
        # compute viscosities, stresses, density gradients
        Erebus.get_viscosities_stresses_density_gradients!(
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
            dRHOYdy,
        )
        # verification, from HTM-planetary.m, line 832ff, 905ff
        for j in 1:1:Nx, i in 1:1:Ny
            # x-Stokes
            if i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1
                # pass: external points
            else
                # x-Stokes internal points
                # Computational viscosity
                ETA1=ETA[i - 1, j]*GGG[i - 1, j]*dt/(GGG[i - 1, j]*dt+ETA[i - 1, j])
                ETA2=ETA[i, j]*GGG[i, j]*dt/(GGG[i, j]*dt+ETA[i, j])
                ETAP1=ETAP[i, j]*GGGP[i, j]*dt/(GGGP[i, j]*dt+ETAP[i, j])
                ETAP2=ETAP[i, j + 1]*GGGP[i, j + 1]*dt/(GGGP[i, j + 1]*dt+ETAP[i, j + 1])
                # Old stresses
                SXY1=SXY0[i - 1, j]*ETA[i - 1, j]/(GGG[i - 1, j]*dt+ETA[i - 1, j])
                SXY2=SXY0[i, j]*ETA[i, j]/(GGG[i, j]*dt+ETA[i, j])
                SXX1=SXX0[i, j]*ETAP[i, j]/(GGGP[i, j]*dt+ETAP[i, j])
                SXX2=SXX0[i, j + 1]*ETAP[i, j + 1]/(GGGP[i, j + 1]*dt+ETAP[i, j + 1])
                # Density gradients
                dRHOdx=(RHOX[i, j + 1]-RHOX[i, j - 1])/2/dx
                dRHOdy=(RHOX[i + 1, j]-RHOX[i - 1, j])/2/dy
                # test
                @test ETAcomp[i - 1, j] ≈ ETA1 rtol=1e-9
                @test ETAcomp[i, j] ≈ ETA2 rtol=1e-9
                @test ETAPcomp[i, j] ≈ ETAP1 rtol=1e-9
                @test ETAPcomp[i, j + 1] ≈ ETAP2 rtol=1e-9
                @test SXYcomp[i - 1, j] ≈ SXY1 rtol=1e-9
                @test SXYcomp[i, j] ≈ SXY2 rtol=1e-9
                @test SXXcomp[i, j] ≈ SXX1 rtol=1e-9
                @test SXXcomp[i, j + 1] ≈ SXX2 rtol=1e-9
                @test dRHOXdx[i, j] ≈ dRHOdx rtol=1e-9
                @test dRHOXdy[i, j] ≈ dRHOdy rtol=1e-9
            end
            # y-Stokes
            if j==1 || j==Nx1 || i==1 || i==Ny || i==Ny1
                # pass: external points
            else
                # Computational viscosity
                ETA1=ETA[i, j - 1]*GGG[i, j - 1]*dt/(GGG[i, j - 1]*dt+ETA[i, j - 1])
                ETA2=ETA[i, j]*GGG[i, j]*dt/(GGG[i, j]*dt+ETA[i, j])
                ETAP1=ETAP[i, j]*GGGP[i, j]*dt/(GGGP[i, j]*dt+ETAP[i, j])
                ETAP2=ETAP[i + 1, j]*GGGP[i + 1, j]*dt/(GGGP[i + 1, j]*dt+ETAP[i + 1, j])
                # Old stresses
                SXY1=SXY0[i, j - 1]*ETA[i, j - 1]/(GGG[i, j - 1]*dt+ETA[i, j - 1])
                SXY2=SXY0[i, j]*ETA[i, j]/(GGG[i, j]*dt+ETA[i, j])
                SYY1=-SXX0[i, j]*ETAP[i, j]/(GGGP[i, j]*dt+ETAP[i, j])
                SYY2=-SXX0[i + 1, j]*ETAP[i + 1, j]/(GGGP[i + 1, j]*dt+ETAP[i + 1, j])
                # Density gradients
                dRHOdx=(RHOY[i, j + 1]-RHOY[i, j - 1])/2/dx
                dRHOdy=(RHOY[i + 1, j]-RHOY[i - 1, j])/2/dy
                # test
                @test ETAcomp[i, j - 1] ≈ ETA1 rtol=1e-9
                @test ETAcomp[i, j] ≈ ETA2 rtol=1e-9
                @test ETAPcomp[i, j] ≈ ETAP1 rtol=1e-9
                @test ETAPcomp[i + 1, j] ≈ ETAP2 rtol=1e-9
                @test SXYcomp[i, j - 1] ≈ SXY1 rtol=1e-9
                @test SXYcomp[i, j] ≈ SXY2 rtol=1e-9
                @test SYYcomp[i, j] ≈ SYY1 rtol=1e-9
                @test SYYcomp[i + 1, j] ≈ SYY2 rtol=1e-9
                @test dRHOYdx[i, j] ≈ dRHOdx rtol=1e-9
                @test dRHOYdy[i, j] ≈ dRHOdy rtol=1e-9
            end
        end
    end # testset "get_viscosities_stresses_density_gradients()"

    @testset "setup_hydromechanical_lse()" begin
        # setup hydromechanical LSE
        R, S = Erebus.setup_hydromechanical_lse()
        # test
        # @test typeof(L) == ExtendableSparseMatrix{Float64, Int64}
        # @test size(L) == (Nx1*Ny1*6, Nx1*Ny1*6)
        @test typeof(R) == Vector{Float64}
        @test size(R) == (Nx1*Ny1*6,)
        @test typeof(S) == Vector{Float64}
        @test size(S) == (Nx1*Ny1*6,)
    end

    @testset "setup_thermal_lse()" begin
        # setup thermal LSE
        RP, SP = Erebus.setup_thermal_lse()
        # test
        # @test typeof(LP) == ExtendableSparseMatrix{Float64, Int64}
        # @test size(LP) == (Nx1*Ny1, Nx1*Ny1)
        @test typeof(RP) == Vector{Float64}
        @test size(RP) == (Nx1*Ny1,)
        @test typeof(SP) == Vector{Float64}
        @test size(SP) == (Nx1*Ny1,)
    end

    @testset "setup_gravitational_lse()" begin
        # setup gravitational LSE
        RT, ST = Erebus.setup_gravitational_lse()
        # test
        # @test typeof(LT) == ExtendableSparseMatrix{Float64, Int64}
        # @test size(LT) == (Nx1*Ny1, Nx1*Ny1)
        @test typeof(RT) == Vector{Float64}
        @test size(RT) == (Nx1*Ny1,)
        @test typeof(ST) == Vector{Float64}
        @test size(ST) == (Nx1*Ny1,)
    end

    @testset "assemble_hydromechanical_lse()" begin
        dt = dt_longest
        # simulate data
        ETA = rand(rgen, Ny, Nx) * 1e16
        ETAP = rand(rgen, Ny1, Nx1) * 1e16
        GGG = rand(rgen, Ny, Nx) * 1e10
        GGGP = rand(rgen, Ny1, Nx1) * 1e10
        SXY0 = rand(rgen, Ny, Nx) * 1e4
        SXX0 = rand(rgen, Ny, Nx) * 1e4
        RHOX = rand(rgen, Ny1, Nx1) * 1e4
        RHOY = rand(rgen, Ny1, Nx1) * 1e4
        RHOFX = rand(rgen, Ny1, Nx1) * 1e4
        RHOFY = rand(rgen, Ny1, Nx1) * 1e4
        RX = rand(rgen, Ny1, Nx1) * 1e39
        RY = rand(rgen, Ny1, Nx1) * 1e39
        ETAPHI = rand(rgen, Ny1, Nx1) * 1e15
        BETTAPHI = rand(rgen, Ny1, Nx1) * 1e-10
        PHI = rand(rgen, Ny1, Nx1) * 1.0
        gx = rand(rgen, Ny1, Nx1) * 1.0
        gy = rand(rgen, Ny1, Nx1) * 1.0
        pr0 = rand(rgen, Ny1, Nx1) * 1e5
        pf0 = rand(rgen, Ny1, Nx1) * 1e5
        DMP = rand(rgen, Ny1, Nx1) * 1e5
        # LSE
        R = zeros(Nx1*Ny1*6)
        L_ver = zeros(Nx1*Ny1*6, Nx1*Ny1*6)
        R_ver = zeros(Nx1*Ny1*6)
        # assemble hydromechanical LSE
        L = Erebus.assemble_hydromechanical_lse!(
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
            BETTAPHI,
            PHI,
            gx,
            gy,
            pr0,
            pf0,
            DMP,
            dt,
            R,
        )
        L = collect(L)
        # verification, from HTM-planetary.m, line 779ff;
        # HTM-hydration.m, line 707ff
        # Hydro-Mechanical Solution
        # Composing global matrixes L_ver[], R_ver[] for Stokes & continuity equations
        for j in 1:1:Nx1
            for i in 1:1:Ny1
                # Define global indexes in algebraic space
                kvx=((j-1)*Ny1+i-1)*6+1 # Vx solid
                kvy=kvx+1 # Vy solid
                kpm=kvx+2 # Ptotal
                kqx=kvx+3 # qx Darcy
                kqy=kvx+4 # qy Darcy
                kpf=kvx+5 # P fluid

                # Vx equation External points
                if i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1
                    # Boundary Condition 
                    # Ghost unknowns 1*Vx=0
                    if j==Nx1
                        L_ver[kvx, kvx]=1 # Left part
                        R_ver[kvx]=0 # Right part
                    end
                    # Left Boundary
                    if j==1
                        L_ver[kvx, kvx]=1 # Left part
                        R_ver[kvx]=vxleft # Right part
                    end
                    # Right Boundary
                    if j==Nx
                        L_ver[kvx, kvx]=1 # Left part
                        R_ver[kvx]=vxright # Right part
                    end
                    # Top boundary
                    if i==1 && j>1 && j<Nx
                        L_ver[kvx, kvx]=1 # Left part
                        L_ver[kvx, kvx + 6]=bctop # Left part
                        R_ver[kvx]=0 # Right part
                    end
                    # Top boundary
                    if i==Ny1 && j>1 && j<Nx
                        L_ver[kvx, kvx]=1 # Left part
                        L_ver[kvx, kvx - 6]=bcbottom # Left part
                        R_ver[kvx]=0 # Right part
                    end
                else
                    # Internal points: x-Stokes eq.
                    #            Vx2
                    #             |
                    #        Vy1  |  Vy3
                    #             |
                    #     Vx1-P1-Vx3-P2-Vx5
                    #             |
                    #        Vy2  |  Vy4
                    #             |
                    #            Vx4
                    #
                    # Computational viscosity
                    ETA1=ETA[i - 1, j]*GGG[i - 1, j]*dt/(GGG[i - 1, j]*dt+ETA[i - 1, j])
                    ETA2=ETA[i, j]*GGG[i, j]*dt/(GGG[i, j]*dt+ETA[i, j])
                    ETAP1=ETAP[i, j]*GGGP[i, j]*dt/(GGGP[i, j]*dt+ETAP[i, j])
                    ETAP2=ETAP[i, j + 1]*GGGP[i, j + 1]*dt/(
                        GGGP[i, j + 1]*dt+ETAP[i, j + 1]
                    )
                    # Old stresses
                    SXY1=SXY0[i - 1, j]*ETA[i - 1, j]/(GGG[i - 1, j]*dt+ETA[i - 1, j])
                    SXY2=SXY0[i, j]*ETA[i, j]/(GGG[i, j]*dt+ETA[i, j])
                    SXX1=SXX0[i, j]*ETAP[i, j]/(GGGP[i, j]*dt+ETAP[i, j])
                    SXX2=SXX0[i, j + 1]*ETAP[i, j + 1]/(GGGP[i, j + 1]*dt+ETAP[i, j + 1])
                    # Density gradients
                    dRHOdx=(RHOX[i, j + 1]-RHOX[i, j - 1])/2/dx
                    dRHOdy=(RHOX[i + 1, j]-RHOX[i - 1, j])/2/dy
                    # Left part
                    L_ver[kvx, kvx - Ny1 * 6]=ETAP1/dx^2 # Vx1
                    L_ver[kvx, kvx - 6]=ETA1/dy^2 # Vx2
                    L_ver[kvx, kvx]=-(ETAP1+ETAP2)/dx^2 - (ETA1+ETA2)/dy^2 -
                                    dRHOdx*gx[i, j]*dt # Vx3
                    L_ver[kvx, kvx + 6]=ETA2/dy^2 # Vx4
                    L_ver[kvx, kvx + Ny1 * 6]=ETAP2/dx^2 # Vx5
                    L_ver[kvx, kvy]=ETAP1/dx/dy-ETA2/dx/dy-dRHOdy*gx[i, j]*dt/4  # Vy2
                    L_ver[kvx, kvy + Ny1 * 6]=-ETAP2/dx/dy+ETA2/dx/dy-dRHOdy*gx[i, j]*dt/4  # Vy4
                    L_ver[kvx, kvy - 6]=-ETAP1/dx/dy+ETA1/dx/dy-dRHOdy*gx[i, j]*dt/4  # Vy1
                    L_ver[kvx, kvy + Ny1 * 6 - 6]=ETAP2/dx/dy-ETA1/dx/dy-dRHOdy*gx[i, j]*dt/4  # Vy3
                    L_ver[kvx, kpm]=Kcont/dx # P1
                    L_ver[kvx, kpm + Ny1 * 6]=-Kcont/dx # P2
                    # Right part
                    R_ver[kvx]=-RHOX[i, j]*gx[i, j]-(SXY2-SXY1)/dy-(SXX2-SXX1)/dx
                end

                # Vy equation External points
                if j==1 || j==Nx1 || i==1 || i==Ny || i==Ny1
                    # Boundary Condition
                    # Ghost unknowns 1*Vx=0
                    if i==Ny1
                        L_ver[kvy, kvy]=1 # Left part
                        R_ver[kvy]=0 # Right part
                    end
                    # Top boundary
                    if i==1
                        L_ver[kvy, kvy]=1 # Left part
                        R_ver[kvy]=vytop # Right part
                    end
                    # Bottom boundary
                    if i==Ny
                        L_ver[kvy, kvy]=1 # Left part
                        R_ver[kvy]=vybottom # Right part
                    end
                    # Left boundary
                    if j==1 && i>1 && i<Ny
                        L_ver[kvy, kvy]=1 # Left part
                        L_ver[kvy, kvy + 6 * Ny1]=bcleft # Left part
                        R_ver[kvy]=0 # Right part
                    end
                    # Right boundary
                    if j==Nx1 && i>1 && i<Ny
                        L_ver[kvy, kvy]=1 # Left part
                        L_ver[kvy, kvy - 6 * Ny1]=bcright # Left part
                        R_ver[kvy]=0 # Right part
                    end
                else
                    # Internal points: y-Stokes eq.
                    #            Vy2
                    #             |
                    #         Vx1 P1 Vx3
                    #             |
                    #     Vy1----Vy3----Vy5
                    #             |
                    #         Vx2 P2 Vx4
                    #             |
                    #            Vy4
                    #
                    # Computational viscosity
                    ETA1=ETA[i, j - 1]*GGG[i, j - 1]*dt/(GGG[i, j - 1]*dt+ETA[i, j - 1])
                    ETA2=ETA[i, j]*GGG[i, j]*dt/(GGG[i, j]*dt+ETA[i, j])
                    ETAP1=ETAP[i, j]*GGGP[i, j]*dt/(GGGP[i, j]*dt+ETAP[i, j])
                    ETAP2=ETAP[i + 1, j]*GGGP[i + 1, j]*dt/(
                        GGGP[i + 1, j]*dt+ETAP[i + 1, j]
                    )
                    # Old stresses
                    SXY1=SXY0[i, j - 1]*ETA[i, j - 1]/(GGG[i, j - 1]*dt+ETA[i, j - 1])
                    SXY2=SXY0[i, j]*ETA[i, j]/(GGG[i, j]*dt+ETA[i, j])
                    SYY1=-SXX0[i, j]*ETAP[i, j]/(GGGP[i, j]*dt+ETAP[i, j])
                    SYY2=-SXX0[i + 1, j]*ETAP[i + 1, j]/(GGGP[i + 1, j]*dt+ETAP[i + 1, j])
                    # Density gradients
                    dRHOdx=(RHOY[i, j + 1]-RHOY[i, j - 1])/2/dx
                    dRHOdy=(RHOY[i + 1, j]-RHOY[i - 1, j])/2/dy
                    # Left part
                    L_ver[kvy, kvy - Ny1 * 6]=ETA1/dx^2 # Vy1
                    L_ver[kvy, kvy - 6]=ETAP1/dy^2 # Vy2
                    L_ver[kvy, kvy]=-(ETAP1+ETAP2)/dy^2 - (ETA1+ETA2)/dx^2 -
                                    dRHOdy*gy[i, j]*dt # Vy3
                    L_ver[kvy, kvy + 6]=ETAP2/dy^2 # Vy4
                    L_ver[kvy, kvy + Ny1 * 6]=ETA2/dx^2 # Vy5
                    L_ver[kvy, kvx]=ETAP1/dx/dy-ETA2/dx/dy-dRHOdx*gy[i, j]*dt/4 #Vx3
                    L_ver[kvy, kvx + 6]=-ETAP2/dx/dy+ETA2/dx/dy-dRHOdx*gy[i, j]*dt/4 #Vx4
                    L_ver[kvy, kvx - Ny1 * 6]=-ETAP1/dx/dy+ETA1/dx/dy-dRHOdx*gy[i, j]*dt/4 #Vx1
                    L_ver[kvy, kvx + 6 - Ny1 * 6]=ETAP2/dx/dy-ETA1/dx/dy-dRHOdx*gy[i, j]*dt/4 #Vx2
                    L_ver[kvy, kpm]=Kcont/dy # P1
                    L_ver[kvy, kpm + 6]=-Kcont/dy # P2

                    # Right part
                    R_ver[kvy]=-RHOY[i, j]*gy[i, j]-(SXY2-SXY1)/dx-(SYY2-SYY1)/dy
                end

                # P equation External points
                if (
                    i==1 ||
                    j==1 ||
                    i==Ny1 ||
                    j==Nx1 ||
                    (i==2 && j>=2 && j<=Nx) ||
                    (i==Ny && j>=2 && j<=Nx) ||
                    (j==2 && i>=2 && i<=Ny) ||
                    (j==Nx && i>=2 && i<=Ny)
                )
                    # Boundary Condition
                    # 1*P=0
                    L_ver[kpm, kpm]=1 # Left part
                    R_ver[kpm]=0 # Right part
                    # Real BC
                    if (
                        (i==2 && j>=2 && j<=Nx) ||
                        (i==Ny && j>=2 && j<=Nx) ||
                        (j==2 && i>=2 && i<=Ny) ||
                        (j==Nx && i>=2 && i<=Ny)
                    )
                        L_ver[kpm, kpm]=1*Kcont # Left part
                        R_ver[kpm]=psurface # Right part
                    end
                else
                    # Internal points: continuity eq.
                    # dVx/dx+dVy/dy=0
                    #            Vy1
                    #             |
                    #        Vx1--P--Vx2
                    #             |
                    #            Vy2
                    #
                    # Left part
                    L_ver[kpm, kvx - Ny1 * 6]=-1/dx # Vx1
                    L_ver[kpm, kvx]=1/dx # Vx2
                    L_ver[kpm, kvy - 6]=-1/dy # Vy1
                    L_ver[kpm, kvy]=1/dy # Vy2
                    L_ver[kpm, kpm] = Kcont/(1-PHI[i, j])*(1/ETAPHI[i, j]+BETTAPHI[i, j]/dt) # Ptotal
                    L_ver[kpm, kpf]=-Kcont/(1-PHI[i, j])*(1/ETAPHI[i, j]+BETTAPHI[i, j]/dt) # Pfluid
                    # Right part
                    R_ver[kpm]=(pr0[i, j]-pf0[i, j])/(1-PHI[i, j])*BETTAPHI[i, j]/dt+DMP[
                        i, j
                    ]
                end

                # qxDarcy equation External points
                if i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1
                    # Boundary Condition
                    # 1*qx=0
                    L_ver[kqx, kqx]=1 # Left part
                    R_ver[kqx]=0 # Right part
                    # Top boundary
                    if i==1 && j>1 && j<Nx
                        L_ver[kqx, kqx + 6]=bcftop # Left part
                    end
                    # Bottom boundary
                    if i==Ny1 && j>1 && j<Nx
                        L_ver[kqx, kqx - 6]=bcfbottom # Left part
                    end
                else
                    # Internal points: x-Darcy eq.
                    # Rx*qxDarcy+dP/dx=RHOfluid*gx
                    #     P1-qxD-P2
                    # Left part
                    L_ver[kqx, kqx]=RX[i, j] # qxD
                    L_ver[kqx, kpf]=-Kcont/dx # P1
                    L_ver[kqx, kpf + Ny1 * 6]=Kcont/dx # P2
                    # Right part
                    R_ver[kqx]=RHOFX[i, j]*gx[i, j]
                end

                # qyDarcy equation External points
                if j==1 || j==Nx1 || i==1 || i==Ny || i==Ny1
                    # Boundary Condition
                    # 1*Vy=0
                    L_ver[kqy, kqy]=1 # Left part
                    R_ver[kqy]=0 # Right part
                    # Left boundary
                    if j==1 && i>1 && i<Ny
                        L_ver[kqy, kqy + 6 * Ny1]=bcfleft # Left part
                    end
                    # Right boundary
                    if j==Nx1 && i>1 && i<Ny
                        L_ver[kqy, kqy - 6 * Ny1]=bcfright # Left part
                    end
                else
                    # Internal points: y-Stokes eq.
                    # Internal points: x-Darcy eq.
                    # Rx*qxDarcy+dP/dx=RHOfluid*gx
                    #      P1
                    #      |
                    #     qxD
                    #      |
                    #      P2
                    # Left part
                    L_ver[kqy, kqy]=RY[i, j] # qxD
                    L_ver[kqy, kpf]=-Kcont/dy # P1
                    L_ver[kqy, kpf + 6]=Kcont/dy # P
                    # Right part
                    R_ver[kqy]=RHOFY[i, j]*gy[i, j]
                end

                # Pfluid equation External points
                if (
                    i==1 ||
                    j==1 ||
                    i==Ny1 ||
                    j==Nx1 ||
                    (i==2 && j>=2 && j<=Nx) ||
                    (i==Ny && j>=2 && j<=Nx) ||
                    (j==2 && i>=2 && i<=Ny) ||
                    (j==Nx && i>=2 && i<=Ny)
                )
                    # Boundary Condition
                    # 1*Pfluid=0
                    L_ver[kpf, kpf]=1 # Left part
                    R_ver[kpf]=0 # Right part
                    # Real BC
                    if (
                        (i==2 && j>=2 && j<=Nx) ||
                        (i==Ny && j>=2 && j<=Nx) ||
                        (j==2 && i>2 && i<Ny) ||
                        (j==Nx && i>2 && i<Ny)
                    )
                        L_ver[kpf, kpf]=1*Kcont #Left part
                        R_ver[kpf]=psurface # Right part
                    end
                else
                    # Internal points: continuity eq.
                    # dqxD/dx+dqyD/dy-(Ptotal-Pfluid)/ETHAphi=0
                    #            qyD1
                    #              |
                    #        qxD1--P--qxD2
                    #              |
                    #            qyD2
                    #
                    # Left part
                    L_ver[kpf, kqx - Ny1 * 6]=-1/dx # qxD1
                    L_ver[kpf, kqx]=1/dx # qxD2
                    L_ver[kpf, kqy - 6]=-1/dy # qyD1
                    L_ver[kpf, kqy]=1/dy # qyD2
                    L_ver[kpf, kpm]=-Kcont/(1-PHI[i, j])*(1/ETAPHI[i, j]+BETTAPHI[i, j]/dt) # Ptotal
                    L_ver[kpf, kpf] = Kcont/(1-PHI[i, j])*(1/ETAPHI[i, j]+BETTAPHI[i, j]/dt) # Pfluid
                    # Right part
                    R_ver[kpf]=-(pr0[i, j]-pf0[i, j])/(1-PHI[i, j])*BETTAPHI[i, j]/dt
                end
            end
        end
        # test
        # for j=1:1:Nx1*Ny1*6, i=1:1:Nx1*Ny1*6
        #     @test L[i, j] ≈ L_ver[i, j] rtol=1e-9
        # @test R[i] ≈ R_ver[i] rtol=1e-9
        # end
        @test L ≈ L_ver rtol=1e-9
        @test R ≈ R_ver rtol=1e-9
    end # testset "assemble_hydromechanical_lse()"

    @testset "process_hydromechanical_solution!()" begin
        # simulate data
        S = rand(rgen, Nx1*Ny1*6) * 2e-10 .- 1e-10
        vx = zeros(Ny1, Nx1)
        vy = zeros(Ny1, Nx1)
        pr = zeros(Ny1, Nx1)
        qxD = zeros(Ny1, Nx1)
        qyD = zeros(Ny1, Nx1)
        pf = zeros(Ny1, Nx1)
        vx_ver = zeros(Ny1, Nx1)
        vy_ver = zeros(Ny1, Nx1)
        pr_ver = zeros(Ny1, Nx1)
        qxD_ver = zeros(Ny1, Nx1)
        qyD_ver = zeros(Ny1, Nx1)
        pf_ver = zeros(Ny1, Nx1)
        # process solution
        Erebus.process_hydromechanical_solution!(S, vx, vy, pr, qxD, qyD, pf)
        # verification, from HTM-planetary.m, line 1058ff
        for j in 1:1:Nx1
            for i in 1:1:Ny1
                # Define global indexes in algebraic space
                kvx=((j-1)*Ny1+i-1)*6+1 # Vx solid
                kvy=kvx+1 # Vy solid
                kpm=kvx+2 # Ptotal
                kqx=kvx+3 # qx Darcy
                kqy=kvx+4 # qy Darcy
                kpf=kvx+5 # P fluid
                # Reload solution
                vx_ver[i, j]=S[kvx]
                vy_ver[i, j]=S[kvy]
                pr_ver[i, j]=S[kpm]*Kcont
                qxD_ver[i, j]=S[kqx]
                qyD_ver[i, j]=S[kqy]
                pf_ver[i, j]=S[kpf]*Kcont
            end
        end
        # pavr=(pf_ver[2,2]+pf_ver[2,Nx]+pf_ver[Ny,2]+pf_ver[Ny,Nx])/4;
        # @. pr_ver=pr_ver-pavr+psurface;
        # @. pf_ver=pf_ver-pavr+psurface;
        # test
        for j in 1:1:Nx1, i in 1:1:Ny1
            @test vx[i, j] ≈ vx_ver[i, j] rtol=1e-9
            @test vy[i, j] ≈ vy_ver[i, j] rtol=1e-9
            @test pr[i, j] ≈ pr_ver[i, j] rtol=1e-9
            @test qxD[i, j] ≈ qxD_ver[i, j] rtol=1e-9
            @test qyD[i, j] ≈ qyD_ver[i, j] rtol=1e-9
            @test pf[i, j] ≈ pf_ver[i, j] rtol=1e-9
        end
    end # testset "process_hydromechanical_solution!()"

    @testset "compute_Aϕ!()" begin
        dt = dt_longest
        # simulate data
        APHI = rand(rgen, Ny1, Nx1)
        APHI_ver = rand(rgen, Ny1, Nx1)
        ETAPHI = rand(rgen, Ny1, Nx1) * 1e15
        BETTAPHI = rand(rgen, Ny1, Nx1) * 1e-10
        PHI = rand(rgen, Ny1, Nx1)
        pr = rand(rgen, Ny1, Nx1) * 1e4
        pf = rand(rgen, Ny1, Nx1) * 1e4
        pr0 = rand(rgen, Ny1, Nx1) * 1e4
        pf0 = rand(rgen, Ny1, Nx1) * 1e4
        # compute Aϕ
        aphimax = Erebus.compute_Aϕ!(APHI, ETAPHI, BETTAPHI, PHI, pr, pf, pr0, pf0, dt)
        # verification, from HTM-planetary.m, line 1078ff
        APHI_ver = zeros(Ny1, Nx1)
        aphimax_ver=0
        for j in 2:1:Nx
            for i in 2:1:Ny
                APHI_ver[i, j]=(
                    (pr[i, j]-pf[i, j])/ETAPHI[i, j] +
                    ((pr[i, j]-pr0[i, j])-(pf[i, j]-pf0[i, j]))/dt*BETTAPHI[i, j]
                )/(1-PHI[i, j])/PHI[i, j]
                aphimax_ver=max(aphimax_ver, abs(APHI_ver[i, j]))
            end
        end
        # test incompressible baseline
        for j in 2:1:Nx, i in 2:1:Ny
            @test APHI[i, j] ≈ APHI_ver[i, j] rtol=1e-9
        end
        @test aphimax ≈ aphimax_ver rtol=1e-9

        # Test poroelastic compaction with nonzero betasolid
        betasolid = 2.5e-11
        aphimax_poro = Erebus.compute_Aϕ!(
            APHI, ETAPHI, BETTAPHI, PHI, pr, pf, pr0, pf0, dt; betasolid=betasolid
        )
        for j in 2:1:Nx, i in 2:1:Ny
            bd = (BETTAPHI[i, j] + betasolid) / (1.0 - PHI[i, j])
            kbw = 1.0 - betasolid / bd
            comp =
                (pr[i, j] - pf[i, j]) / (ETAPHI[i, j] * (1.0 - PHI[i, j])) +
                bd * ((pr[i, j] - pr0[i, j]) - kbw * (pf[i, j] - pf0[i, j])) / dt
            @test APHI[i, j] ≈ (comp / PHI[i, j]) rtol=1e-9
        end
        @test aphimax_poro != aphimax
    end # testset "compute_Aϕ!()"

    @testset "compute_fluid_velocities!()" begin
        # simulate data
        PHIX = rand(rgen, Ny1, Nx1)
        PHIY = rand(rgen, Ny1, Nx1)
        qxD = rand(rgen, Ny1, Nx1)
        qyD = rand(rgen, Ny1, Nx1)
        vx = rand(rgen, Ny1, Nx1)
        vy = rand(rgen, Ny1, Nx1)
        vxf = zeros(Ny1, Nx1)
        vyf = zeros(Ny1, Nx1)
        vxf_ver = zeros(Ny1, Nx1)
        vyf_ver = zeros(Ny1, Nx1)
        # compute fluid velocities
        Erebus.compute_fluid_velocities!(PHIX, PHIY, qxD, qyD, vx, vy, vxf, vyf)
        # verification, from HTM-planetary.m line 1090ff
        for j in 1:1:Nx
            for i in 2:1:Ny
                vxf_ver[i, j]=qxD[i, j]/PHIX[i, j]
            end
        end
        # Apply BC
        # Top
        vxf_ver[1, :] = -bcftop*vxf_ver[2, :]
        # Bottom
        vxf_ver[Ny1, :] = -bcfbottom*vxf_ver[Ny, :]
        # Vy fluid
        for j in 2:1:Nx
            for i in 1:1:Ny
                vyf_ver[i, j]=qyD[i, j]/PHIY[i, j]
            end
        end
        # Apply BC
        # Left
        vyf_ver[:, 1] = -bcfleft*vyf_ver[:, 2]
        # Right
        vyf_ver[:, Nx1] = -bcfright*vyf_ver[:, Nx]
        # Add solid velocity
        # vxf0=vxf; vxf=vxf+vx
        vxf_ver.=vxf_ver .+ vx
        # vyf0=vyf; vyf=vyf+vy
        vyf_ver.=vyf_ver .+ vy
        # test
        for j in 1:1:Nx1, i in 1:1:Ny1
            @test vxf[i, j] ≈ vxf_ver[i, j] rtol=1e-9
            @test vyf[i, j] ≈ vyf_ver[i, j] rtol=1e-9
        end
    end # testset "compute_fluid_velocity!()"

    @testset "compute_displacement_timestep(): CFL and porosity-rate constraints" begin
        dt = dt_longest

        # 1. Zero velocity and zero porosity rate: timestep preserved exactly
        dt_zero = Erebus.compute_displacement_timestep(
            zeros(Ny1, Nx1), zeros(Ny1, Nx1), zeros(Ny1, Nx1), zeros(Ny1, Nx1), dt, 0.0
        )
        @test dt_zero ≈ dt rtol=1e-12

        # 2. Solid velocity CFL constraint: dt * max(|vx|) <= dxymax * dx
        U_fast = 100.0 * (dxymax * dx / dt)
        vx_fast = fill(U_fast, Ny1, Nx1)
        dt_vx = Erebus.compute_displacement_timestep(
            vx_fast, zeros(Ny1, Nx1), zeros(Ny1, Nx1), zeros(Ny1, Nx1), dt, 0.0
        )
        @test dt_vx * U_fast ≈ dxymax * dx rtol=1e-12
        @test dt_vx < dt

        # 3. Fluid velocity CFL constraint: dt * max(|vxf|) <= dxymax * dx
        V_fast = 200.0 * (dxymax * dx / dt)
        vxf_fast = fill(V_fast, Ny1, Nx1)
        dt_vxf = Erebus.compute_displacement_timestep(
            zeros(Ny1, Nx1), zeros(Ny1, Nx1), vxf_fast, zeros(Ny1, Nx1), dt, 0.0
        )
        @test dt_vxf * V_fast ≈ dxymax * dx rtol=1e-12
        @test dt_vxf < dt

        # 4. Porosity rate constraint: dt * aphimax <= dphimax
        aphi_fast = 50.0 * (dphimax / dt)
        dt_phi = Erebus.compute_displacement_timestep(
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            zeros(Ny1, Nx1),
            dt,
            aphi_fast,
        )
        @test dt_phi * aphi_fast ≈ dphimax rtol=1e-12
        @test dt_phi < dt

        # 5. Monotonicity and positivity under combined random loads
        aphimax = rand(rgen)
        vx = rand(rgen, Ny1, Nx1)
        vy = rand(rgen, Ny1, Nx1)
        vxf = rand(rgen, Ny1, Nx1)
        vyf = rand(rgen, Ny1, Nx1)
        dtm = Erebus.compute_displacement_timestep(vx, vy, vxf, vyf, dt, aphimax)
        @test 0.0 < dtm <= dt
        @test dtm * maximum(abs, vx) <= dxymax * dx + 1e-12
        @test dtm * maximum(abs, vy) <= dxymax * dy + 1e-12
        @test dtm * maximum(abs, vxf) <= dxymax * dx + 1e-12
        @test dtm * maximum(abs, vyf) <= dxymax * dy + 1e-12
        @test dtm * aphimax <= dphimax + 1e-12
    end # testset "compute_displacement_timestep()"

    @testset "compute_stress_strainrate!()" begin
        dtm = dt_longest
        # simulate data
        vx = rand(rgen, Ny1, Nx1)
        vy = rand(rgen, Ny1, Nx1)
        ETA = rand(rgen, Ny, Nx)
        GGG = rand(rgen, Ny, Nx)
        ETAP = rand(rgen, Ny1, Nx1)
        GGGP = rand(rgen, Ny1, Nx1)
        SXX0 = rand(rgen, Ny1, Nx1)
        SXY0 = rand(rgen, Ny, Nx)
        EXY = rand(rgen, Ny, Nx)
        SXY = rand(rgen, Ny, Nx)
        DSXY = rand(rgen, Ny, Nx)
        EXX = rand(rgen, Ny1, Nx1)
        SXX = rand(rgen, Ny1, Nx1)
        DSXX = rand(rgen, Ny1, Nx1)
        EII = zeros(Ny1, Nx1)
        SII = zeros(Ny1, Nx1)
        # compute stress, strainrate
        Erebus.compute_stress_strainrate!(
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
            dtm,
        )
        # verification, from HTM-planetary.m, line 1144ff
        EXY_ver = zeros(Ny, Nx) # Strain rate EPSILONxy, 1/s
        SXY_ver = zeros(Ny, Nx) # Stress SIGMAxy, Pa
        DSXY_ver = zeros(Ny, Nx) # Stress change SIGMAxy, Pa
        for j in 1:1:Nx
            for i in 1:1:Ny
                # EXY;SXY; DSXY
                EXY_ver[i, j]=0.5*((vx[i + 1, j]-vx[i, j])/dy+(vy[i, j + 1]-vy[i, j])/dx)
                SXY_ver[i, j]=2*ETA[i, j]*EXY_ver[i, j]*GGG[i, j]*dtm/(
                    GGG[i, j]*dtm+ETA[i, j]
                )+SXY0[i, j]*ETA[i, j]/(GGG[i, j]*dtm+ETA[i, j])
                DSXY_ver[i, j]=SXY_ver[i, j]-SXY0[i, j]
            end
        end
        # Compute EPSILONxx; SIGMA'xx in pressure nodes
        EXX_ver = zeros(Ny1, Nx1) # Strain rate EPSILONxx, 1/s
        EII_ver = zeros(Ny1, Nx1) # Second strain rate invariant, 1/s
        SXX_ver = zeros(Ny1, Nx1) # Stress SIGMA'xx, Pa
        SII_ver = zeros(Ny1, Nx1) # Second stress invariant, Pa
        DSXX_ver = zeros(Ny1, Nx1) # Stress change SIGMA'xx, Pa
        DIVV_ver = zeros(Ny1, Nx1) # div[v]
        for j in 2:1:Nx
            for i in 2:1:Ny
                # DIVV
                DIVV_ver[i, j]=(vx[i, j]-vx[i, j - 1])/dx+(vy[i, j]-vy[i - 1, j])/dy
                # EXX
                EXX_ver[i, j]=((vx[i, j]-vx[i, j - 1])/dx-(vy[i, j]-vy[i - 1, j])/dy)/2
                # SXX
                SXX_ver[i, j]=2*ETAP[i, j]*EXX_ver[i, j]*GGGP[i, j]*dtm/(
                    GGGP[i, j]*dtm+ETAP[i, j]
                )+SXX0[i, j]*ETAP[i, j]/(GGGP[i, j]*dtm+ETAP[i, j])
                DSXX_ver[i, j]=SXX_ver[i, j]-SXX0[i, j]
                # EII
                EII_ver[i, j]=(
                    EXX_ver[i, j]^2+(
                        (
                            EXY_ver[i, j]+EXY_ver[i - 1, j]+EXY_ver[i, j - 1]+EXY_ver[
                                i - 1, j - 1
                            ]
                        )/4
                    )^2
                )^0.5
                # SII
                SII_ver[i, j]=(
                    SXX_ver[i, j]^2+(
                        (
                            SXY_ver[i, j]+SXY_ver[i - 1, j]+SXY_ver[i, j - 1]+SXY_ver[
                                i - 1, j - 1
                            ]
                        )/4
                    )^2
                )^0.5
            end
        end
        # test
        for j in 1:1:Nx, i in 1:1:Ny
            @test EXY[i, j] ≈ EXY_ver[i, j] rtol=1e-9
            @test SXY[i, j] ≈ SXY_ver[i, j] rtol=1e-9
            @test DSXY[i, j] ≈ DSXY_ver[i, j] rtol=1e-9
        end
        for j in 2:1:Nx, i in 2:1:Ny
            @test EXX[i, j] ≈ EXX_ver[i, j] rtol=1e-9
            @test SXX[i, j] ≈ SXX_ver[i, j] rtol=1e-9
            @test EII[i, j] ≈ EII_ver[i, j] rtol=1e-9
            @test SII[i, j] ≈ SII_ver[i, j] rtol=1e-9
        end
    end # testset "compute_stress_strainrate!()"

    @testset "symmetrize_p_node_observables!()" begin
        # simulate data
        SXX = rand(rgen, Ny1, Nx1)
        APHI = rand(rgen, Ny1, Nx1)
        PHI = rand(rgen, Ny1, Nx1)
        pr = rand(rgen, Ny1, Nx1)
        pf = rand(rgen, Ny1, Nx1)
        ps = zeros(Ny1, Nx1)
        SXX_ver = deepcopy(SXX)
        APHI_ver = deepcopy(APHI)
        PHI_ver = deepcopy(PHI)
        pr_ver = deepcopy(pr)
        pf_ver = deepcopy(pf)
        ps_ver = zeros(Ny1, Nx1)
        # symmetrize p node variables
        Erebus.symmetrize_p_node_observables!(SXX, APHI, PHI, pr, pf, ps)
        # verification, from HTM-planetary.m, line 1196ff
        # Apply Symmetry to Pressure nodes
        # External P-nodes: symmetry
        # Top
        SXX_ver[1, 2:Nx]=SXX_ver[2, 2:Nx]
        APHI_ver[1, 2:Nx]=APHI_ver[2, 2:Nx]
        PHI_ver[1, 2:Nx]=PHI_ver[2, 2:Nx]
        pr_ver[1, 2:Nx]=pr_ver[2, 2:Nx]
        pf_ver[1, 2:Nx]=pf_ver[2, 2:Nx]
        # Bottom
        SXX_ver[Ny1, 2:Nx]=SXX_ver[Ny, 2:Nx]
        APHI_ver[Ny1, 2:Nx]=APHI_ver[Ny, 2:Nx]
        PHI_ver[Ny1, 2:Nx]=PHI_ver[Ny, 2:Nx]
        pr_ver[Ny1, 2:Nx]=pr_ver[Ny, 2:Nx]
        pf_ver[Ny1, 2:Nx]=pf_ver[Ny, 2:Nx]
        # Left
        SXX_ver[:, 1]=SXX_ver[:, 2]
        APHI_ver[:, 1]=APHI_ver[:, 2]
        PHI_ver[:, 1]=PHI_ver[:, 2]
        pr_ver[:, 1]=pr_ver[:, 2]
        pf_ver[:, 1]=pf_ver[:, 2]
        # Right
        SXX_ver[:, Nx1]=SXX_ver[:, Nx]
        APHI_ver[:, Nx1]=APHI_ver[:, Nx]
        PHI_ver[:, Nx1]=PHI_ver[:, Nx]
        pr_ver[:, Nx1]=pr_ver[:, Nx]
        pf_ver[:, Nx1]=pf_ver[:, Nx]
        # Compute solid pressure
        ps_ver=(pr_ver .- pf_ver .* PHI_ver) ./ (1 .- PHI_ver)
        # test
        @test SXX[1, 2:Nx] == SXX_ver[1, 2:Nx]
        @test APHI[1, 2:Nx] == APHI_ver[1, 2:Nx]
        @test PHI[1, 2:Nx] == PHI_ver[1, 2:Nx]
        @test pr[1, 2:Nx] == pr_ver[1, 2:Nx]
        @test pf[1, 2:Nx] == pf_ver[1, 2:Nx]
        @test SXX[Ny1, 2:Nx] == SXX_ver[Ny1, 2:Nx]
        @test APHI[Ny1, 2:Nx] == APHI_ver[Ny1, 2:Nx]
        @test PHI[Ny1, 2:Nx] == PHI_ver[Ny1, 2:Nx]
        @test pr[Ny1, 2:Nx] == pr_ver[Ny1, 2:Nx]
        @test pf[Ny1, 2:Nx] == pf_ver[Ny1, 2:Nx]
        @test SXX[:, 1] == SXX_ver[:, 1]
        @test APHI[:, 1] == APHI_ver[:, 1]
        @test PHI[:, 1] == PHI_ver[:, 1]
        @test pr[:, 1] == pr_ver[:, 1]
        @test pf[:, 1] == pf_ver[:, 1]
        @test SXX[:, Nx1] == SXX_ver[:, Nx1]
        @test APHI[:, Nx1] == APHI_ver[:, Nx1]
        @test PHI[:, Nx1] == PHI_ver[:, Nx1]
        @test pr[:, Nx1] == pr_ver[:, Nx1]
        @test pf[:, Nx1] == pf_ver[:, Nx1]
        @test ps ≈ ps_ver rtol=1e-9
    end # testset "symmetrize_p_node_observables!()"

    @testset "positive_max(): non-negativity and upper bound axioms" begin
        A = rand(rgen, -100:0.1:100, 50, 50)
        B = rand(rgen, -100:0.1:100, 50, 50)
        R = zeros(50, 50)
        Erebus.positive_max!(A, B, R)

        # Axiom 1: Non-negativity everywhere
        @test all(r -> r >= 0.0, R)

        # Axiom 2: Upper bound (R >= A and R >= B everywhere)
        @test all(R .>= A)
        @test all(R .>= B)

        # Axiom 3: Negative saturation: strictly negative inputs produce exact zero
        A_neg = fill(-10.0, 10, 10)
        B_neg = fill(-25.0, 10, 10)
        R_neg = fill(99.0, 10, 10)
        Erebus.positive_max!(A_neg, B_neg, R_neg)
        @test all(iszero, R_neg)

        # Axiom 4: Positive maximum selection
        A_pos = [2.0 5.0; 8.0 1.0]
        B_pos = [3.0 1.0; 6.0 4.0]
        R_pos = zeros(2, 2)
        Erebus.positive_max!(A_pos, B_pos, R_pos)
        @test R_pos ≈ [3.0 5.0; 8.0 4.0] rtol=1e-12

        # Axiom 5: Commutativity / symmetry: positive_max(A, B) == positive_max(B, A)
        R_sym = zeros(50, 50)
        Erebus.positive_max!(B, A, R_sym)
        @test R ≈ R_sym rtol=1e-12
    end # testset "positive_max()"

    @testset "compute_nodal_adjustment!()" begin
        dt = dt_longest
        iplast = 1
        # simulate data
        ETA = rand(rgen, Ny, Nx)
        ETA0 = rand(rgen, Ny, Nx)
        ETA5 = zeros(Ny, Nx)
        GGG = rand(rgen, Ny, Nx)
        SXX = rand(rgen, Ny1, Nx1)
        SXY = rand(rgen, Ny, Nx)
        pr = rand(rgen, Ny1, Nx1)
        pf = rand(rgen, Ny1, Nx1)
        COH = rand(rgen, Ny, Nx)
        TEN = rand(rgen, Ny, Nx)
        FRI = rand(rgen, Ny, Nx)
        YNY = zeros(Bool, Ny, Nx)
        YNY5 = zeros(Bool, Ny, Nx)
        DSY = rand(rgen, Ny, Nx)
        YERRNOD = zeros(nplast)
        YERRNOD_ver = zeros(nplast)
        # compute nodal adjustment
        complete = Erebus.compute_nodal_adjustment!(
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
            iplast,
        )
        # verification, from HTM-planetary.m, line 1232ff
        ETA5_ver=copy(ETA0)
        YNY5_ver = zeros(Bool, Ny, Nx)
        DSY_ver = zeros(Ny, Nx)
        ynpl=0
        ddd=0
        for i in 1:1:Ny
            for j in 1:1:Nx
                # Compute second stress invariant
                SIIB_ver=(
                    SXY[i, j]^2+(
                        (SXX[i, j]+SXX[i + 1, j]+SXX[i, j + 1]+SXX[i + 1, j + 1])/4
                    )^2
                )^0.5
                # Compute second invariant for a purely elastic stress build-up
                siiel_ver=SIIB_ver*(GGG[i, j]*dt+ETA[i, j])/ETA[i, j]
                # Compute total & fluid pressure
                prB_ver=(pr[i, j]+pr[i + 1, j]+pr[i, j + 1]+pr[i + 1, j + 1])/4
                pfB_ver=(pf[i, j]+pf[i + 1, j]+pf[i, j + 1]+pf[i + 1, j + 1])/4
                # Compute yielding stress
                syieldc_ver=COH[i, j]+FRI[i, j]*(prB_ver-pfB_ver) # Confined fracture
                syieldt_ver=TEN[i, j]+(prB_ver-pfB_ver) # Tensile fracture
                syield_ver=max(min(syieldt_ver, syieldc_ver), 0) # Non-negative strength requirement
                # Update error for old yielding nodes
                ynn=0
                if (YNY[i, j]>0)
                    ynn=1
                    DSY_ver[i, j]=SIIB_ver-syield_ver
                    ddd=ddd+DSY_ver[i, j]^2
                    ynpl=ynpl+1
                end
                # Correcting viscosity for yielding
                if syield_ver<siiel_ver
                    # New viscosity for the basic node
                    etapl_ver=dt*GGG[i, j]*syield_ver/(siiel_ver-syield_ver)
                    if etapl_ver<ETA0[i, j]
                        # Recompute nodal visocity
                        ETA5_ver[i, j]=etapl_ver^(1-etawt)*ETA[i, j]^etawt
                        # Mark yielding nodes
                        YNY5_ver[i, j]=1
                        # Apply viscosity cutoff values
                        if ETA5_ver[i, j]<etamin
                            ETA5_ver[i, j]=etamin
                        elseif ETA5_ver[i, j]>etamax
                            ETA5_ver[i, j]=etamax
                        end
                        # Update Error for new yielding nodes
                        if ynn==0
                            DSY_ver[i, j]=SIIB_ver-syield_ver
                            ddd=ddd+DSY_ver[i, j]^2
                            ynpl=ynpl+1
                        end
                    else
                        ETA5_ver[i, j]=ETA0[i, j]
                        YNY5_ver[i, j]=0
                    end
                else
                    ETA5_ver[i, j]=ETA0[i, j]
                    YNY5_ver[i, j]=0
                end
            end
        end
        # Compute yielding error for markers
        if (ynpl>0)
            YERRNOD_ver[iplast]=(ddd/ynpl)^0.5
        end
        # test
        @test ETA5 ≈ ETA5_ver rtol=1e-9
        @test YNY5 == YNY5_ver
        @test YERRNOD[iplast] ≈ YERRNOD_ver[iplast] rtol=1e-9
        @test complete == (ynpl==0 || iplast==nplast || YERRNOD[iplast]<yerrmax)
    end # testset "compute_nodal_adjustment!()"

    @testset "finalize_plastic_iteration_pass!()" begin
        dt = dt_longest
        iplast = 1
        # simulate data
        ETA = zeros(Ny, Nx)
        ETA5 = rand(rgen, Ny, Nx)
        ETA00 = rand(rgen, Ny, Nx)
        YNY = zeros(Bool, Ny, Nx)
        YNY5 = rand(rgen, Bool, Ny, Nx)
        YNY00 = rand(rgen, Bool, Ny, Nx)
        YNY_inv_ETA = zeros(Ny, Nx)
        ETA_ver = zeros(Ny, Nx)
        YNY_ver = zeros(Bool, Ny, Nx)
        YNY_inv_ETA_ver = zeros(Ny, Nx)
        dt_ver = dt_longest
        # finalize_plastic_iteration_pass
        dt = Erebus.finalize_plastic_iteration_pass!(
            ETA, ETA5, ETA00, YNY, YNY5, YNY00, YNY_inv_ETA, dt, iplast
        )
        # verification, from HTM-planetary.m, line 1301ff
        if trunc(Int, iplast/dtstep)*dtstep==iplast
            # Decrease timestep
            dt_ver=dt_ver/dtkoef
            # Reset old viscoplastic viscosity
            ETA_ver=ETA00
            YNY_ver=YNY00
        else
            # Use new viscoplastic viscosity
            ETA_ver=ETA5
            YNY_ver=YNY5
        end
        YNY_inv_ETA_ver.=YNY_ver ./ ETA_ver
        # test
        @test dt == dt_ver
        @test ETA == ETA_ver
        @test YNY == YNY_ver
        @test YNY_inv_ETA == YNY_inv_ETA_ver
    end # testset "finalize_plastic_iteration_pass!()"

    @testset "finalize_thermochemical_iteration_pass(): thermal relaxation step control" begin
        # 1. First iteration (titer == 1) with sub-threshold temperature change: dt unchanged
        dt_sub = Erebus.finalize_thermochemical_iteration_pass(DTmax * 0.5, dt_longest, 1)
        @test dt_sub ≈ dt_longest rtol=1e-12

        # 2. First iteration (titer == 1) with excessive temperature change: dt scaled proportionally
        maxDT_large = DTmax * 4.0
        dt_cut = Erebus.finalize_thermochemical_iteration_pass(maxDT_large, dt_longest, 1)
        @test dt_cut ≈ dt_longest * (DTmax / maxDT_large) rtol=1e-12
        @test dt_cut < dt_longest

        # 3. Subsequent iterations (titer > 1): dt is preserved regardless of DT magnitude
        dt_iter2 = Erebus.finalize_thermochemical_iteration_pass(maxDT_large, dt_longest, 2)
        dt_iter3 = Erebus.finalize_thermochemical_iteration_pass(maxDT_large, dt_longest, 3)
        @test dt_iter2 ≈ dt_longest rtol=1e-12
        @test dt_iter3 ≈ dt_longest rtol=1e-12

        # 4. Positivity and monotonicity: larger excess DT produces strictly smaller dt in titer 1
        dt_cut_larger = Erebus.finalize_thermochemical_iteration_pass(
            DTmax * 8.0, dt_longest, 1
        )
        @test dt_cut_larger < dt_cut
        @test dt_cut_larger > 0.0
    end # testset "finalize_thermochemical_iteration_pass()"

    @testset "compute_thermochemical_iteration_outcome: convergence decision boundaries" begin
        # 1. Converged case: small pressure error and past titer threshold (titer > 2)
        pf_converged = rand(rgen, Ny1, Nx1)
        pf0 = copy(pf_converged) # zero pressure error
        DMP_active = fill(1.0e-5, Ny1, Nx1)
        @test Erebus.compute_thermochemical_iteration_outcome(
            DMP_active, pf_converged, pf0, 3
        ) == true

        # 2. Converged case: small pressure error and zero mass reaction (DMP <= 0) even at titer 1
        DMP_zero = zeros(Ny1, Nx1)
        @test Erebus.compute_thermochemical_iteration_outcome(
            DMP_zero, pf_converged, pf0, 1
        ) == true

        # 3. Non-converged case: excessive pressure error (pferrcur >= pferrmax)
        pf_diverged = copy(pf0)
        pf_diverged[2, 2] += pferrmax * 2.0
        @test Erebus.compute_thermochemical_iteration_outcome(
            DMP_active, pf_diverged, pf0, 3
        ) == false
        @test Erebus.compute_thermochemical_iteration_outcome(
            DMP_zero, pf_diverged, pf0, 1
        ) == false

        # 4. Non-converged case: active reaction (DMP > 0) at early iterations (titer <= 2)
        @test Erebus.compute_thermochemical_iteration_outcome(
            DMP_active, pf_converged, pf0, 1
        ) == false
        @test Erebus.compute_thermochemical_iteration_outcome(
            DMP_active, pf_converged, pf0, 2
        ) == false
    end # testset "compute_thermochemical_iteration_outcome"

    @testset "assemble_thermal_lse!" begin
        dt = rand(rgen)
        tk1 = rand(rgen, Ny1, Nx1)
        RHOCP = rand(rgen, Ny1, Nx1)
        KX = rand(rgen, Ny1, Nx1)
        KY = rand(rgen, Ny1, Nx1)
        HR = rand(rgen, Ny1, Nx1)
        HA = rand(rgen, Ny1, Nx1)
        HS = rand(rgen, Ny1, Nx1)
        DHP = rand(rgen, Ny1, Nx1)
        RT = rand(rgen, Ny1*Nx1)
        LT = Erebus.assemble_thermal_lse!(tk1, RHOCP, KX, KY, HR, HA, HS, DHP, RT, dt)
        # verification, from HTM-planetary.m, line 1618ff
        LT_ver = zeros(Ny1*Nx1, Ny1*Nx1)
        RT_ver = zeros(Ny1*Nx1)
        # Composing global matrixes LT[], RT[]
        # Going through all points of the 2D grid &
        # composing respective equations
        for j in 1:1:Nx1
            for i in 1:1:Ny1
                # Define global index in algebraic space
                gk=(j-1)*Ny1+i
                # External points
                if i==1 || i==Ny1 || j==1 || j==Nx1
                    # Boundary Condition
                    # Top BC: T=273
                    if i==1 && j>1 && j<Nx1
                        LT_ver[gk, gk]=1 # Left part
                        LT_ver[gk, gk + 1]=-1 # Left part
                        RT_ver[gk]=0 # Right part
                    end
                    # Bottom BC: T=1500
                    if i==Ny1 && j>1 && j<Nx1
                        LT_ver[gk, gk]=1 # Left part
                        LT_ver[gk, gk - 1]=-1 # Left part
                        RT_ver[gk]=0 # Right part
                    end
                    # Left BC: dT/dx=0
                    if j==1
                        LT_ver[gk, gk]=1 # Left part
                        LT_ver[gk, gk + Ny1]=-1 # Left part
                        RT_ver[gk]=0 # Right part
                    end
                    # Right BC: dT/dx=0
                    if j==Nx1
                        LT_ver[gk, gk]=1 # Left part
                        LT_ver[gk, gk - Ny1]=-1 # Left part
                        RT_ver[gk]=0 # Right part
                    end
                else
                    # Internal points: Temperature eq.
                    # RHO*CP*dT/dt=-dqx/dx-dqy/dy+Hr+Hs+Ha
                    #          Tdt2
                    #           |
                    #          Ky1
                    #           |
                    #Tdt1-Kx1-T03;Tdt3-Kx2-Tdt5
                    #           |
                    #          Ky2
                    #           |
                    #          Tdt4
                    #
                    # Left part
                    Kx1=KX[i, j - 1]
                    Kx2=KX[i, j]
                    Ky1=KY[i - 1, j]
                    Ky2=KY[i, j]
                    LT_ver[gk, gk - Ny1]=-Kx1/dx^2 # T1
                    LT_ver[gk, gk - 1]=-Ky1/dy^2 # FI2
                    LT_ver[gk, gk]=RHOCP[i, j]/dt+(Kx1+Kx2)/dx^2+(Ky1+Ky2)/dy^2 # FI3
                    LT_ver[gk, gk + 1]=-Ky2/dy^2 # FI4
                    LT_ver[gk, gk + Ny1]=-Kx2/dx^2 # FI5
                    # Right part
                    RT_ver[gk]=RHOCP[i, j]/dt*tk1[i, j]+HR[i, j]+HA[i, j]+HS[i, j]+DHP[i, j]
                end
            end
        end
        # testing 
        @test LT ≈ LT_ver rtol=1e-9
        @test RT ≈ RT_ver rtol=1e-9
    end # testset "assemble_thermal_lse!"

    @testset "poroelastic hydromechanical coupling" begin
        Ny, Nx = Erebus.Ny, Erebus.Nx
        Ny1, Nx1 = Erebus.Ny1, Erebus.Nx1
        dt = 10.0

        ETA = fill(1e22, Ny, Nx)
        ETAP = fill(1e22, Ny1, Nx1)
        GGG = fill(1e10, Ny, Nx)
        GGGP = fill(1e10, Ny1, Nx1)
        SXY0 = zeros(Ny, Nx)
        SXX0 = zeros(Ny, Nx)
        RHOX = fill(3000.0, Ny1, Nx1)
        RHOY = fill(3000.0, Ny1, Nx1)
        RHOFX = fill(1000.0, Ny1, Nx1)
        RHOFY = fill(1000.0, Ny1, Nx1)
        RX = fill(1e-3 / 1e-13, Ny1, Nx1)
        RY = fill(1e-3 / 1e-13, Ny1, Nx1)
        PHI = fill(0.1, Ny1, Nx1)
        ETAPHI = fill(1e22 / 0.1, Ny1, Nx1)
        BETAPHI = fill(0.1 / 1e10, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = zeros(Ny1, Nx1)
        pr0 = zeros(Ny1, Nx1)
        pf0 = fill(1e6, Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)
        R = zeros(Nx1*Ny1*6)

        betasolid = 2.5e-11
        betafluid = 4.0e-10
        L = Erebus.assemble_hydromechanical_lse!(
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
            R;
            betasolid=betasolid,
            betafluid=betafluid,
        )

        # Independent theoretical coefficients (calculated without calling solver functions)
        bd_th = (0.1 / 1e10 + betasolid) / (1.0 - 0.1)
        kbw_th = 1.0 - betasolid / bd_th
        ksk_th = (bd_th - betasolid) / ((bd_th - betasolid) + 0.1 * (betafluid - betasolid))
        C_expected = -Kcont * (inv(1e22 / 0.1) / (1.0 - 0.1) + bd_th * kbw_th / dt)
        D_pm_expected = Kcont * (inv(1e22 / 0.1) / (1.0 - 0.1) + bd_th / dt)
        D_pf_expected =
            Kcont * (inv(1e22 / 0.1) / (1.0 - 0.1) + bd_th * kbw_th / ksk_th / dt)

        # Verify each matrix block independently against theoretical values
        for j in 4:(Nx - 2), i in 4:(Ny - 2)
            kvx = ((j-1)*Ny1 + i-1) * 6 + 1
            kpm = kvx + 2
            kpf = kvx + 5
            @test L[kpm, kpf] ≈ C_expected rtol=1e-10
            @test L[kpf, kpm] ≈ C_expected rtol=1e-10
            @test L[kpm, kpm] ≈ D_pm_expected rtol=1e-10
            @test L[kpf, kpf] ≈ D_pf_expected rtol=1e-10
        end

        # Verify solution solvability and finiteness
        prob = LinearProblem(L, R)
        sol = solve(prob, UMFPACKFactorization())
        @test !any(isnan, sol.u)
        @test !any(isinf, sol.u)
    end # testset "poroelastic hydromechanical coupling"

    @testset "Terzaghi 1D consolidation numerical simulation verification" begin
        Ny, Nx = Erebus.Ny, Erebus.Nx
        Ny1, Nx1 = Erebus.Ny1, Erebus.Nx1
        dy = Erebus.dy

        # Height between draining boundary anchors i=2 and i=Ny
        H = (Ny - 2) * dy
        k_perm = 1e-13
        eta_f = 1e-3
        betasolid = 2.5e-11
        betafluid = 4.0e-10
        G_p = 1e10
        phi_0 = 0.1
        beta_phi = phi_0 / G_p

        bd = Erebus.compute_drained_compressibility(beta_phi, phi_0, betasolid)
        kbw = Erebus.compute_biot_willis_coefficient(bd, betasolid)
        ksk = Erebus.compute_skempton_coefficient(bd, phi_0, betasolid, betafluid)
        S = bd * kbw / ksk
        c_v = k_perm / (eta_f * S)

        u0 = 1.0e6
        dt = 2.0e7 # Timestep [s] scaled with column height (H^2/cv)

        ETA = fill(1e25, Ny, Nx)
        ETAP = fill(1e25, Ny1, Nx1)
        GGG = fill(1e10, Ny, Nx)
        GGGP = fill(1e10, Ny1, Nx1)
        SXY0 = zeros(Ny, Nx)
        SXX0 = zeros(Ny, Nx)
        RHOX = fill(3000.0, Ny1, Nx1)
        RHOY = fill(3000.0, Ny1, Nx1)
        RHOFX = fill(1000.0, Ny1, Nx1)
        RHOFY = fill(1000.0, Ny1, Nx1)
        RX = fill(eta_f / k_perm, Ny1, Nx1)
        RY = fill(eta_f / k_perm, Ny1, Nx1)
        PHI = fill(phi_0, Ny1, Nx1)
        ETAPHI = fill(1e25, Ny1, Nx1)
        BETAPHI = fill(beta_phi, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = zeros(Ny1, Nx1)
        pr0 = zeros(Ny1, Nx1)
        pf0 = fill(u0, Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)
        R = zeros(Nx1*Ny1*6)

        vx = zeros(Ny1, Nx1)
        vy = zeros(Ny1, Nx1)
        pr = zeros(Ny1, Nx1)
        qxD = zeros(Ny1, Nx1)
        qyD = zeros(Ny1, Nx1)
        pf = zeros(Ny1, Nx1)

        nsteps = 3
        t_total = nsteps * dt

        # Step Erebus numerical solver forward in time
        for step in 1:nsteps
            R .= 0.0
            L = Erebus.assemble_hydromechanical_lse!(
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
                R;
                betasolid=betasolid,
                betafluid=betafluid,
            )
            prob = LinearProblem(L, R)
            sol = solve(prob, UMFPACKFactorization())
            Erebus.process_hydromechanical_solution!(sol.u, vx, vy, pr, qxD, qyD, pf)
            pf0 .= pf
        end

        # Analytical 1D consolidation Fourier series solution
        function analytical_2drain(y, t, H_col, c_coeff, p0; nterms=100)
            val = 0.0
            for m in 0:nterms
                M = (2m + 1) * pi
                val +=
                    (4.0 * p0 / M) * sin(M * y / H_col) * exp(-M^2 * c_coeff * t / H_col^2)
            end
            return val
        end

        # Verify numerical solution against analytical solution across column depth
        # Account for psurface boundary condition at draining anchors
        mid_col = div(Nx1, 2)
        for i in 3:(Ny - 1)
            y = (i - 2) * dy
            u_ana = analytical_2drain(y, t_total, H, c_v, u0 - psurface) + psurface
            u_num = pf[i, mid_col]
            rel_err = abs(u_num - u_ana) / u0
            @test rel_err < 0.035
        end

        # Physical Invariant 1: Peak excess pore pressure decays monotonically over consolidation
        @test maximum(pf[:, mid_col]) < u0

        # Physical Invariant 2: Positivity and bounded pressure profile inside column
        @test all(pf[3:(Ny - 1), mid_col] .>= psurface - 1e-6)
        @test all(isfinite, pf)

        # Physical Invariant 3: Spatial symmetry of pressure profile across draining boundaries
        top_idx = 3
        bot_idx = Ny - 1
        @test isapprox(pf[top_idx, mid_col], pf[bot_idx, mid_col]; rtol=0.05)
    end # testset "Terzaghi 1D consolidation numerical simulation verification"

    @testset "assemble_hydromechanical_lse!() dynamic hydrofracture integration" begin
        eta_f = 1.0e-3
        k_perm = 1.0e-12
        k_frac_max = 1.0e-9
        kappa_frac = 100.0
        gamma_frac = 1.0
        sigma_t_val = 5.0e6

        ETA = fill(1e22, Ny, Nx)
        ETAP = fill(1e22, Ny1, Nx1)
        GGG = fill(1e10, Ny, Nx)
        GGGP = fill(1e10, Ny1, Nx1)
        SXY0 = zeros(Ny, Nx)
        SXX0 = zeros(Ny, Nx)
        RHOX = fill(3000.0, Ny1, Nx1)
        RHOY = fill(3000.0, Ny1, Nx1)
        RHOFX = fill(1000.0, Ny1, Nx1)
        RHOFY = fill(1000.0, Ny1, Nx1)
        RX = fill(eta_f / k_perm, Ny1, Nx1)
        RY = fill(eta_f / k_perm, Ny1, Nx1)
        KX = fill(k_perm, Ny1, Nx1)
        KY = fill(k_perm, Ny1, Nx1)
        PHI = fill(0.1, Ny1, Nx1)
        ETAPHI = fill(1e25, Ny1, Nx1)
        BETAPHI = fill(1e-10, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = zeros(Ny1, Nx1)
        pr0 = zeros(Ny1, Nx1)
        pf0 = zeros(Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)
        dt = 1e8
        R = zeros(Nx1*Ny1*6)

        pr = fill(1.0e7, Ny1, Nx1)
        TEN = fill(sigma_t_val, Ny, Nx)

        L_base = Erebus.assemble_hydromechanical_lse!(
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
            R;
            hydrofracture=false,
        )

        pf_intact = fill(5.0e6, Ny1, Nx1)
        L_intact = Erebus.assemble_hydromechanical_lse!(
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
            R;
            hydrofracture=true,
            pr=pr,
            pf=pf_intact,
            TEN=TEN,
            KX=KX,
            KY=KY,
            kappa_frac=kappa_frac,
            gamma_frac=gamma_frac,
            k_frac_max=k_frac_max,
        )
        i_test, j_test = 5, 5
        kqx = ((j_test - 1)*Ny1 + i_test - 1) * 6 + 4
        kqy = kqx + 1
        @test L_intact[kqx, kqx] ≈ L_base[kqx, kqx]
        @test L_intact[kqy, kqy] ≈ L_base[kqy, kqy]

        pf_frac = fill(2.0e7, Ny1, Nx1)
        L_frac = Erebus.assemble_hydromechanical_lse!(
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
            R;
            hydrofracture=true,
            pr=pr,
            pf=pf_frac,
            TEN=TEN,
            KX=KX,
            KY=KY,
            kappa_frac=kappa_frac,
            gamma_frac=gamma_frac,
            k_frac_max=k_frac_max,
        )
        expected_keff = k_perm * 101.0
        expected_rx = RX[i_test, j_test] * (k_perm / expected_keff)
        @test L_frac[kqx, kqx] < L_base[kqx, kqx]
        @test L_frac[kqx, kqx] ≈ expected_rx rtol=1e-12
        @test L_frac[kqy, kqy] ≈ expected_rx rtol=1e-12

        pf_extreme = fill(2.0e8, Ny1, Nx1)
        L_extreme = Erebus.assemble_hydromechanical_lse!(
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
            R;
            hydrofracture=true,
            pr=pr,
            pf=pf_extreme,
            TEN=TEN,
            KX=KX,
            KY=KY,
            kappa_frac=kappa_frac,
            gamma_frac=gamma_frac,
            k_frac_max=k_frac_max,
        )
        expected_ceiling_rx = eta_f / k_frac_max
        @test L_extreme[kqx, kqx] ≈ expected_ceiling_rx rtol=1e-12
        @test L_extreme[kqy, kqy] ≈ expected_ceiling_rx rtol=1e-12

        HS_intact = zeros(Ny, Nx)
        HS_frac = zeros(Ny, Nx)
        SXY_test = fill(1e5, Ny, Nx)
        SXX_test = fill(1e5, Ny1, Nx1)
        qxD_test = fill(1e-6, Ny1, Nx1)
        qyD_test = fill(1e-6, Ny1, Nx1)
        Erebus.compute_shear_heating!(
            HS_intact,
            ETA,
            SXY_test,
            ETAP,
            SXX_test,
            RX,
            RY,
            qxD_test,
            qyD_test,
            PHI,
            ETAPHI,
            pr,
            pf_intact;
            hydrofracture=true,
            TEN=TEN,
            KX=KX,
            KY=KY,
            kappa_frac=kappa_frac,
            gamma_frac=gamma_frac,
            k_frac_max=k_frac_max,
        )
        Erebus.compute_shear_heating!(
            HS_frac,
            ETA,
            SXY_test,
            ETAP,
            SXX_test,
            RX,
            RY,
            qxD_test,
            qyD_test,
            PHI,
            ETAPHI,
            pr,
            pf_frac;
            hydrofracture=true,
            TEN=TEN,
            KX=KX,
            KY=KY,
            kappa_frac=kappa_frac,
            gamma_frac=gamma_frac,
            k_frac_max=k_frac_max,
        )
        @test HS_frac[5, 5] < HS_intact[5, 5]
    end # testset "assemble_hydromechanical_lse!() dynamic hydrofracture integration"

    @testset "assemble_hydromechanical_lse!() Darcy thermal buoyancy integration" begin
        g_accel = 10.0
        rho0_f = 1000.0
        alpha_f = 2.0e-4
        delta_T = 50.0
        eta_f = 1.0e-3
        k_perm = 1.0e-13

        ETA = fill(1e22, Ny, Nx)
        ETAP = fill(1e22, Ny1, Nx1)
        GGG = fill(1e10, Ny, Nx)
        GGGP = fill(1e10, Ny1, Nx1)
        SXY0 = zeros(Ny, Nx)
        SXX0 = zeros(Ny, Nx)
        RHOX = fill(3000.0, Ny1, Nx1)
        RHOY = fill(3000.0, Ny1, Nx1)
        PHI = fill(0.1, Ny1, Nx1)
        ETAPHI = fill(1e25, Ny1, Nx1)
        BETAPHI = fill(1e-10, Ny1, Nx1)
        gx = zeros(Ny1, Nx1)
        gy = fill(g_accel, Ny1, Nx1)
        pr0 = zeros(Ny1, Nx1)
        pf0 = zeros(Ny1, Nx1)
        DMP = zeros(Ny1, Nx1)
        dt = 1e8
        R_const = zeros(Nx1*Ny1*6)
        R_buoy = zeros(Nx1*Ny1*6)

        RX = fill(eta_f / k_perm, Ny1, Nx1)
        RY = fill(eta_f / k_perm, Ny1, Nx1)

        RHOFX_const = fill(rho0_f, Ny1, Nx1)
        RHOFY_const = fill(rho0_f, Ny1, Nx1)
        Erebus.assemble_hydromechanical_lse!(
            ETA,
            ETAP,
            GGG,
            GGGP,
            SXY0,
            SXX0,
            RHOX,
            RHOY,
            RHOFX_const,
            RHOFY_const,
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
            R_const,
        )

        rho_f_hot = rho0_f * (1.0 - alpha_f * delta_T)
        RHOFX_buoy = copy(RHOFX_const)
        RHOFY_buoy = copy(RHOFY_const)
        RHOFY_buoy[Ny - 2, 5] = rho_f_hot

        Erebus.assemble_hydromechanical_lse!(
            ETA,
            ETAP,
            GGG,
            GGGP,
            SXY0,
            SXX0,
            RHOX,
            RHOY,
            RHOFX_buoy,
            RHOFY_buoy,
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
            R_buoy,
        )

        kqy_test = ((5 - 1)*Ny1 + (Ny-2) - 1) * 6 + 5
        delta_R = R_const[kqy_test] - R_buoy[kqy_test]
        expected_delta_R = (rho0_f - rho_f_hot) * g_accel
        @test delta_R ≈ expected_delta_R rtol=1e-12
        @test delta_R ≈ (rho0_f * alpha_f * delta_T * g_accel) rtol=1e-12
    end # testset "assemble_hydromechanical_lse!() Darcy thermal buoyancy integration"
end

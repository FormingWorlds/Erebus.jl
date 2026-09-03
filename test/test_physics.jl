@testset "Physics" begin
    @testset "distance()" begin
        @test Erebus.distance(0, 0, 0, 0) == 0
        @test Erebus.distance(1, 0, 0, 0) == 1
        @test Erebus.distance(0, 1, 0, 0) == 1
        @test Erebus.distance(0, 0, 1, 0) == 1
        @test Erebus.distance(0, 0, 0, 1) == 1
        @test Erebus.distance(0, 0, 1, 1) ≈ sqrt(2) rtol=1e-9
        @test Erebus.distance(1, 1, 0, 0) ≈ sqrt(2) rtol=1e-9
        @test Erebus.distance(-1, -1, 1, 1) ≈ sqrt(8) rtol=1e-9
        # random tests
        num_samples = 1000
        p1 = rand(rgen, 0:0.001:1000, 2, num_samples)
        p2 = rand(rgen, 0:0.001:1000, 2, num_samples)
        for i in 1:num_samples
            @test Erebus.distance(
                p1[:, i]..., p2[:, i]...) ≈ (
                    (p1[1, i]-p2[1, i])^2+(p1[2, i]-p2[2, i])^2)^0.5 rtol=1e-9
        end
    end # testset "distance()"

    @testset "total()" begin
        @test Erebus.total(0, 0, 0) == 0
        @test Erebus.total(1, 0, 0) == 1
        @test Erebus.total(0, 1, 0) == 0
        @test Erebus.total(0, 0, 1) == 0
        @test Erebus.total(1, 2, 0.5) == 1.5
        for i=1:1:1000
            s, f, ϕ = rand(rgen, 3)
            @test Erebus.total(s, f, ϕ) ≈ s*(1-ϕ)+f*ϕ rtol=1e-9
        end
    end # testset "total()"

    @testset "dot4()" begin
        for i=1:1:1000
            v1 = rand(rgen, 4)
            v2 = rand(rgen, 4)
            @test Erebus.dot4(v1, v2) ≈ sum(v1.*v2) rtol=1e-9
        end
    end

    @testset "ktotal()" begin
        # verification, from HTM-planetary.m, line 1761
        for i=1:1000
            ksolid, kfluid, phi = rand(rgen, 3)
            @test Erebus.ktotal(ksolid, kfluid, phi) ≈ (
                ksolid*kfluid/2+((ksolid*(3*phi-2)+kfluid*(1-3*phi))^2)/16
                )^0.5-(ksolid*(3*phi-2)+kfluid*(1-3*phi))/4 rtol=1e-9
        end
    end # testset "ktotal()"

    @testset "kphi()" begin
        # verification, from HTM-planetary.m, line 333
        for i=1:1000
            kphim0m, phimm = rand(rgen, 2)
            @test Erebus.kphi(kphim0m, phimm) ≈ (
                kphim0m*(phimm/phim0)^3/((1-phimm)/(1-phim0))^2
            ) rtol=1e-9
        end
    end # testset "kphi()"

    @testset "ηᶠcur_inv_kᵠ(kϕᵣ, ϕ, ηᶠcur)" begin
        for i=1:1000
            kϕᵣ, ϕ, ηᶠcur = rand(rgen, 3)
            @test Erebus.ηᶠcur_inv_kᵠ(kϕᵣ, ϕ, ηᶠcur) ≈ (
                ηᶠcur*inv(kϕᵣ*(ϕ/phim0)^3/((1-ϕ)/(1-phim0))^2)
            ) rtol=1e-9
        end
    end # testset "ηᶠcur_inv_kᵠ(kϕᵣ, ϕ, ηᶠcur)"

    @testset "etatotal_rocks()" begin
        tmin = min(tmfluidphase, tmsolidphase) - 10
        tmid = min(tmfluidphase, tmsolidphase) + 0.5*abs(tmfluidphase-tmsolidphase)
        tmax = max(tmfluidphase, tmsolidphase) + 10
        for type in 1:2
            @test Erebus.etatotal_rocks(tmin, type) == max(
                etamin, etasolidm[type], etafluidm[type])
            if tmfluidphase <= tmsolidphase
                @test Erebus.etatotal_rocks(tmid, type) ==
                max(etamin, etasolidm[type], etafluidmm[type])
            else
                @test Erebus.etatotal_rocks(tmid, type) ==
                max(etamin, etasolidmm[type], etafluidm[type])
            end
            @test Erebus.etatotal_rocks(tmax, type) == max(
                etamin, etasolidmm[type], etafluidmm[type])
        end
    end # testset "etatotal_rocks()"

    @testset "Q_radiogenic()" begin
        # verification, from HTM-planetary.m, line 276
        Q(f, ratio, E, tau, timesum)=f*ratio*E*exp(-timesum/tau)/tau
        @test Erebus.Q_radiogenic(1., 2., 3., 4., 5.) == Q(
            1., 2., 3., 4., 5.)
        @test Erebus.Q_radiogenic(1., 2., 3., 4., 0.) == Q(
            1., 2., 3., 4., 0.)
    end # testset "Q_radiogenic()"

    @testset "calculate_radioactive_heating()" begin
        al = false
        fe = false
        v = @SVector [0., 0., 0.]
        @test Erebus.calculate_radioactive_heating(
            al, fe, 1000.) == (v, v)

        al = true
        fe = false        
        Q_al = Erebus.Q_radiogenic(
            f_al, ratio_al, E_al, tau_al, 1000.)
        u = @SVector [Q_al * rhosolidm[1], Q_al * rhosolidm[2], 0.]
        @test Erebus.calculate_radioactive_heating(
            al, fe, 1000.) == (u, v)

        al = false
        fe = true
        Q_fe = Erebus.Q_radiogenic(
            f_fe, ratio_fe, E_fe, tau_fe, 1000.)
        w = @SVector [Q_fe * rhofluidm[1], 0., 0.]
        @test Erebus.calculate_radioactive_heating(
            al, fe, 1000.) == (v, w)
    end # testset "calculate_radioactive_heating()"

    @testset "compute_gibbs_free_energy()" begin
        dHWD = ΔHWD
        dSWD = ΔSWD
        dVWD = ΔVWD
        m = 1
        phi = rand(rgen)
        XWsolidm0 = rand(rgen, 1)
        tknm, pfnm, XDsolidm0 = rand(rgen, 3)
        Δtr = dtreaction = Erebus.compute_Δtreaction(
            tknm, phi, reaction_rate_coeff_mode)
        Δt₁ = dt1 = 0.5 * Δtr
        Δt₂ = dt2 = 1.5 * Δtr    
        ΔGWD₁ = Erebus.compute_gibbs_free_energy(
            tknm, pfnm, XDsolidm0, XWsolidm0[m], Δt₁, Δtr)
        ΔGWD₂ = Erebus.compute_gibbs_free_energy(
            tknm, pfnm, XDsolidm0, XWsolidm0[m], Δt₂, Δtr)
        # verification, from i2visHTM_hydration.m, line 614ff
        # Compute old dG for dehydration reaction: Wsilicate=Dsilicate+H2O
        dGWD0=dHWD-tknm*dSWD+dVWD*pfnm+8.314*tknm*log(XDsolidm0/XWsolidm0[m]);
        # Compute incomplete reaction for too short timestep
        dGWD1=0;
        if dt1<dtreaction
            dGWD1=dGWD0*(1-dt1/dtreaction);
        end
        # Compute old dG for dehydration reaction: Wsilicate=Dsilicate+H2O
        dGWD0=dHWD-tknm*dSWD+dVWD*pfnm+8.314*tknm*log(XDsolidm0/XWsolidm0[m]);
        # Compute incomplete reaction for too short timestep
        dGWD2=0;
        if dt2<dtreaction
            dGWD2=dGWD0*(1-dt2/dtreaction);
        end
        # test
        @test ΔGWD₁ ≈ dGWD1 rtol=1e-9
        @test ΔGWD₂ ≈ dGWD2 rtol=1e-9
    end # testset "compute_gibbs_free_energy()"

    @testset "compute_relative_enthalpy()" begin
        dHWD = ΔHWD
        MH2O = MH₂O
        Xsolid0 = rand(rgen)
        XWsolidm0 = rand(rgen, 1)
        m = 1
        Hᵗ₀ = Erebus.compute_relative_enthalpy(
            Xsolid0, XWsolidm0[m])
        # verification, from i2visHTM_hydration.m, line 612ff
        # Compute old relative ehthalpy of the system
        Htotal0=-Xsolid0*XWsolidm0[m]*dHWD/(MD+MH2O);
        # test
        @test Hᵗ₀ ≈ Htotal0 rtol=1e-9
    end # testset "compute_relative_enthalpy()"

    @testset "compute_reaction_constant()" begin
        dHWD = ΔHWD
        dSWD = ΔSWD
        dVWD = ΔVWD
        ΔGWD = dGWD = rand(rgen)
        tknm, pfnm = rand(rgen, 2)
        KWD = Erebus.compute_reaction_constant(
            tknm, pfnm, dGWD)
        # verification, from i2visHTM_hydration.m, line 623ff
        KWD_ver=exp(-(dHWD-tknm*dSWD+dVWD*pfnm-dGWD)/8.314/tknm);
        # test
        @test KWD ≈ KWD_ver rtol=1e-9
    end # testset "compute_reaction_constant()"

    @testset "compute_rhocpfluidm()" begin
        Ts = [170.0, tmfluidphase-10.0, tmfluidphase-3.0, tmfluidphase, tmfluidphase+3.0, tmfluidphase+10.0, 370.0, 470.0]
        # verification
        function rhocpfluidm_Hobbs(T)
            if T >= 410
                return ρH₂Oᶠ * (-4.67e4 + 333T - 0.731T^2 + 5.4e-4T^3) 
            elseif 410 > T >= tmfluidphase+5
                return ρH₂Oᶠ * 4200
            elseif tmfluidphase+5 > T >= tmfluidphase
                return ρH₂Oᶠ * (4200.0 + 0.1Lᶠ)
            elseif tmfluidphase > T >= tmfluidphase-5
                return ρH₂Oᶠⁱ * (7.67T + 0.1Lᶠ)
            elseif tmfluidphase-5 > T
                return ρH₂Oᶠⁱ * 7.67T
            end
        end
        # test
        for T in Ts
            @test Erebus.compute_rhocpfluidm(
                T, 1) ≈ rhocpfluidm_Hobbs(T) rtol=1e-9
            @test Erebus.compute_rhocpfluidm(
                T, 9) ≈ rhocpfluidm[1] rtol=1e-9
        end
    end # testset "compute_rhocpfluidm()"

    @testset "compute_ksolidm()" begin
        Ts = [170.0, tmfluidphase-10.0, tmfluidphase-3.0, tmfluidphase, tmfluidphase+3.0, tmfluidphase+10.0, 370.0, 470.0]
        # test
        for T in Ts
            @test Erebus.compute_ksolidm(
                T, 1) ≈ 0.73 + 1293.0/(T+77.0) rtol=1e-9
            @test Erebus.compute_ksolidm(
                T, 9) ≈ ksolidm[1] rtol=1e-9
        end
    end # testset "compute_ksolidm()"

    @testset "compute_kfluidm()" begin
        Ts = [170.0, tmfluidphase-10.0, tmfluidphase-3.0, tmfluidphase, tmfluidphase+3.0, tmfluidphase+10.0, 370.0, 470.0]
        # verification
        function kfluidm_Hobbs(T)
            if T >= 410
                return -0.142 + 4.12e-3T - 5.01e-6T^2
            elseif 410 > T >= tmfluidphase
                return -0.581 + 6.34e-3T - 7.93e-6T^2
            elseif tmfluidphase > T
                return  0.465 + 488.0/T
            end
        end
        # test
        for T in Ts
            @test Erebus.compute_kfluidm(
                T, 1) ≈ kfluidm_Hobbs(T) rtol=1e-9
            @test Erebus.compute_kfluidm(
                T, 9) ≈ kfluidm[1] rtol=1e-9
        end
    end # testset "compute_kfluidm()"

    @testset "compute_Δtreaction()" begin
        Ts = [170.0, tmfluidphase-10.0, tmfluidphase-3.0, tmfluidphase, tmfluidphase+3.0, tmfluidphase+10.0, 370.0, 470.0]
        ϕ = rand(rgen)
        # test
        for T in Ts
            @test Erebus.compute_Δtreaction(
                T, ϕ, 1) ≈ -log_completion_rate / (A_I*ϕ) * exp(b_I*(T-c_I)^2)
            @test Erebus.compute_Δtreaction(
                T, ϕ, 2) ≈ -log_completion_rate / (Sxo_B*ϕ) * 2.0^(
                    -(T-To_B)/Tscl_B)
            @test Erebus.compute_Δtreaction(
                T, ϕ, 3) ≈ -log_completion_rate / (
                    Sxo_T*ϕ*exp(Ea_T/RG*(1/To_T-1/T))
                )
            @test Erebus.compute_Δtreaction(
                T, ϕ, 9) ≈ Δtreaction
        end
    end # testset "compute_Δtreaction()"

    @testset "perform_thermochemical_reaction!()" begin
        marknum = start_marknum
        DMP = rand(rgen, Ny1, Nx1)
        DHP = rand(rgen, Ny1, Nx1)
        DMPSUM = rand(rgen, Ny1, Nx1)
        DHPSUM = rand(rgen, Ny1, Nx1)
        WTPSUM = rand(rgen, Ny1, Nx1)
        pf = rand(rgen, 2700:0.1:3300, Ny1, Nx1)
        tk2 = rand(rgen, 263:0.1:310, Ny1, Nx1)
        tm = rand(rgen, 1:3, marknum)
        xm = rand(rgen, -dx:0.1:x[end]+dx, marknum)
        ym = rand(rgen, -dy:0.1:y[end]+dy, marknum)
        XWˢm₀ = rand(rgen, marknum)
        XWˢm = rand(rgen, marknum)
        phim = rand(rgen, marknum)
        phinewm = rand(rgen, marknum)
        pfm₀ = rand(rgen, 2700:0.1:3300, marknum)
        dt = Δt = 0.9 * Δtreaction
        timestep = 1
        titer = 3
        DMP_ver = copy(DMP)
        DHP_ver = copy(DHP)
        DMPSUM_ver = copy(DMPSUM)
        DHPSUM_ver = copy(DHPSUM)
        WTPSUM_ver = copy(WTPSUM)
        phim_ver = copy(phim)
        phinewm_ver = copy(phinewm)
        pfm₀_ver = copy(pfm₀)
        XWˢm₀_ver = copy(XWˢm₀)
        XWˢm_ver = copy(XWˢm)
        VDsolid = VDˢ
        VWsolid = VWˢ
        dHWD = ΔHWD
        dSWD = ΔSWD
        dVWD = ΔVWD
        rhoH2Ofluid = ρH₂Oᶠ
        rhoH2Ofluidice = ρH₂Oᶠⁱ
        VH2Ofluid = VH₂Oᶠ
        VH2Ofluidice = VH₂Oᶠⁱ
        MH2O = MH₂O
        mode = reaction_rate_coeff_mode
        Erebus.perform_thermochemical_reaction!(
            DMP,
            DHP,
            DMPSUM,
            DHPSUM,
            WTPSUM,
            pf,
            tk2,
            tm,
            xm,
            ym,
            XWˢm₀,
            XWˢm,
            phim,
            phinewm,
            pfm₀,
            marknum,
            Δt,
            timestep,
            titer
        )
        # verification, from HTM-hydration.m, line 558ff
        # Compute mass transfer rate
        DMPSUM_ver=zeros(Ny1,Nx1);
        DHPSUM_ver=zeros(Ny1,Nx1);
        WTPSUM_ver=zeros(Ny1,Nx1);
        for m=1:1:marknum
        if tm[m]< (4-1)
            # Rocks
            # Interpolate fluid pressure and temperature to the marker
            # Define [i,j] indexes for the upper left node
            j=trunc(Int, (xm[m]-xp[1])/dx)+1;
            i=trunc(Int, (ym[m]-yp[1])/dy)+1;
            if j<1
                j=1;
            elseif j>Nx
                j=Nx;
            end
            if i<1
                i=1;
            elseif i>Ny
                i=Ny;
            end
            # Compute distances
            dxmj=xm[m]-xp[j];
            dymi=ym[m]-yp[i];
            # Compute weights
            wtmij=(1-dxmj/dx)*(1-dymi/dy);
            wtmi1j=(1-dxmj/dx)*(dymi/dy);    
            wtmij1=(dxmj/dx)*(1-dymi/dy);
            wtmi1j1=(dxmj/dx)*(dymi/dy);
            # Interpolate nodal fluid pressure and temperature to the marker
            pfnm=pf[i,j]*wtmij+pf[i+1,j]*wtmi1j+pf[i,j+1]*wtmij1+pf[i+1,j+1]*wtmi1j1;
            tknm=tk2[i,j]*wtmij+tk2[i+1,j]*wtmi1j+tk2[i,j+1]*wtmij1+tk2[i+1,j+1]*wtmi1j1;
            pfnm=max(pfnm,0);
            # Use pressure from the previous iteration
            if titer>2
                pfnm=pfnm*(1-pfcoeff)+pfm₀_ver[m]*pfcoeff;
            end
            pfm₀_ver[m]=pfnm;
            if tknm>tmfluidphase
                VH2O = VH2Ofluid;
            else
                VH2O = VH2Ofluidice;
            end
            #
            # Thermodynamic computations
            # Compute bulk composition of the solid+fluid system
            XDsolidm0=1-XWˢm₀_ver[m];
            Xfluid0=phim_ver[m]*(XWˢm₀_ver[m]*VWsolid+XDsolidm0*VDsolid)/(
                (1-phim_ver[m])*VH2O+ phim_ver[m]*(
                    XWˢm₀_ver[m]*VWsolid+XDsolidm0*VDsolid));
            Xsolid0=1-Xfluid0;
            XH2Ototal=(XWˢm₀_ver[m]*Xsolid0+Xfluid0)/(1+XWˢm₀_ver[m]*Xsolid0);
            XDtotal=1-XH2Ototal;
            # Compute old density of the solid and fluid
            rhosolid0=(MD+MH2O*XWˢm₀_ver[m])/(
                VWsolid*XWˢm₀_ver[m]+VDsolid*XDsolidm0);
            if tknm>tmfluidphase
                rhofluid=rhoH2Ofluid;
            else
                rhofluid=rhoH2Ofluidice;
            end
            # Compute old relative ehthalpy of the system
            Htotal0=-Xsolid0*XWˢm₀_ver[m]*dHWD/(MD+MH2O);
            # Compute old dG for dehydration reaction: Wsilicate=Dsilicate+H2O
            dGWD0=dHWD-tknm*dSWD+dVWD*pfnm+8.314*tknm*log(
                XDsolidm0/XWˢm₀_ver[m]);
            # Compute incomplete reaction for too short timestep
            dGWD=0;
            dtreaction = Erebus.compute_Δtreaction(
                tknm, phim_ver[m], mode)
            if dt<dtreaction
                dGWD=dGWD0*(1-dt/dtreaction);
            end

            # Compute equilibrium compositions and fluid fraction
            # Dehydration reaction: Wsilicate=Dsilicate+H2O
            KWD=exp(-(dHWD-tknm*dSWD+dVWD*pfnm-dGWD)/8.314/tknm);
            # Solid composition
            XWsolidm1=1/(KWD+1);
            XDsolidm1=1-XWsolidm1;
            # Fluid, Solid molar fraction
            Xsolid1=XDtotal/(1-XDtotal*XWsolidm1);
            Xfluid1=1-Xsolid1;
            # Process fluid-bearing rocks only 
            if Xfluid1>0 && Xfluid1<1 
                # Compute equilibrium Porosity
                phinew1=Xfluid1*VH2O/(
                    Xfluid1*VH2O+Xsolid1*(
                        XWsolidm1*VWsolid+XDsolidm1*VDsolid));
                # Compute equilibrium density of the solid and fluid
                rhosolid=(MD+MH2O*XWsolidm1)/(
                    VWsolid*XWsolidm1+VDsolid*XDsolidm1);
                if tknm>tmfluidphase
                    rhofluid0=rhoH2Ofluid;
                else
                    rhofluid0=rhoH2Ofluidice;
                end
                # Compute equilibrium relative ehthalpy of the system
                Htotal1=-Xsolid1*XWsolidm1*dHWD/(MD+MH2O);
                # Compute ehthalpy change
                dHtotal=Htotal1-Htotal0;
                # Compute old/equilibrium volume ratio
                RV=(rhosolid*(1-phinew1)+rhofluid*phinew1)/(
                    rhosolid0*(1-phim_ver[m])+rhofluid0*phim_ver[m]);
                # Compute mass transfer rate
                Gmass=(rhosolid0*RV*(1-phim_ver[m])-rhosolid*(1-phinew1))/dt;
                DMm=(1-RV)/dt;
                # Compute enthalpy transfer term
                DHm=Gmass*dHtotal;
                # Save new composition and melt fraction
                XWˢm_ver[m]=XWsolidm1; # solid composition
                phinewm_ver[m]=phinew1; # porosity
                # Reset properties at the first timestep
                if timestep==1
                    XWˢm₀_ver[m]=XWˢm_ver[m];
                    phim_ver[m]=phinewm_ver[m];
                end
                # Interpolation to pressure nodes 
                # Update subgrid diffusion on nodes
                # [i,j] Node
                DMPSUM_ver[i,j]=DMPSUM_ver[i,j]+DMm*wtmij;
                DHPSUM_ver[i,j]=DHPSUM_ver[i,j]+DHm*wtmij;
                WTPSUM_ver[i,j]=WTPSUM_ver[i,j]+wtmij;
                # i+1,j Node
                DMPSUM_ver[i+1,j]=DMPSUM_ver[i+1,j]+DMm*wtmi1j;
                DHPSUM_ver[i+1,j]=DHPSUM_ver[i+1,j]+DHm*wtmi1j;
                WTPSUM_ver[i+1,j]=WTPSUM_ver[i+1,j]+wtmi1j;
                # [i,j]+1 Node
                DMPSUM_ver[i,j+1]=DMPSUM_ver[i,j+1]+DMm*wtmij1;
                DHPSUM_ver[i,j+1]=DHPSUM_ver[i,j+1]+DHm*wtmij1;
                WTPSUM_ver[i,j+1]=WTPSUM_ver[i,j+1]+wtmij1;
                # i+1,j+1 Node
                DMPSUM_ver[i+1,j+1]=DMPSUM_ver[i+1,j+1]+DMm*wtmi1j1;
                DHPSUM_ver[i+1,j+1]=DHPSUM_ver[i+1,j+1]+DHm*wtmi1j1;
                WTPSUM_ver[i+1,j+1]=WTPSUM_ver[i+1,j+1]+wtmi1j1;
            end
        end
        end
        # P-nodes
        DMP_ver=zeros(Ny1,Nx1)
        DHP_ver=zeros(Ny1,Nx1);
        for j=1:1:Nx1
            for i=1:1:Ny1
                if WTPSUM_ver[i,j]>0
                    DMP_ver[i,j]=DMPSUM_ver[i,j]/WTPSUM_ver[i,j];
                    DHP_ver[i,j]=DHPSUM_ver[i,j]/WTPSUM_ver[i,j];
                end
            end
        end
        # testing
        @test DMP ≈ DMP_ver rtol=1e-6
        @test DHP ≈ DHP_ver rtol=1e-6
        @test pfm₀ ≈ pfm₀_ver rtol=1e-6
        @test XWˢm₀ ≈ XWˢm₀_ver rtol=1e-6
        @test phim ≈ phim_ver rtol=1e-6
        @test phinewm ≈ phinewm_ver rtol=1e-6
    end # testset "perform_thermochemical_reaction!()"

    @testset "compute_shear_heating!()" begin
        HS = zeros(Ny1, Nx1)
        ETA = rand(rgen, Ny, Nx)
        SXY = rand(rgen, Ny, Nx)
        ETAP = rand(rgen, Ny1, Nx1)
        SXX = rand(rgen, Ny1, Nx1)
        RX = rand(rgen, Ny1, Nx1)
        RY = rand(rgen, Ny1, Nx1)
        qxD = rand(rgen, Ny1, Nx1)
        qyD = rand(rgen, Ny1, Nx1)
        PHI = rand(rgen, Ny1, Nx1)
        ETAPHI = rand(rgen, Ny1, Nx1)
        pr = rand(rgen, Ny1, Nx1)
        pf = rand(rgen, Ny1, Nx1)
        # compute shear compute shear heating
        Erebus.compute_shear_heating!(
            HS,
            ETA,
            SXY,
            ETAP,
            SXX,
            RX,
            RY,
            qxD,
            qyD,
            PHI,
            ETAPHI,
            pr,
            pf
        )
        # verification, from HTM-planetary.m, line 1551ff
        HS_ver = zeros(Ny1, Nx1); # Adiabatic heating, W/m^3
        for j=2:1:Nx
            for i=2:1:Ny
                # Average SXY*EXY
                SXYEXY_ver=(SXY[i,j]^2/ETA[i,j]+SXY[i-1,j]^2/ETA[i-1,j]+ SXY[i,j-1]^2/ETA[i,j-1]+SXY[i-1,j-1]^2/ETA[i-1,j-1])/4
                # HS
                HS_ver[i,j]=SXX[i,j]^2/ETAP[i,j]+SXYEXY_ver+ (pr[i,j]-pf[i,j])^2/(1-PHI[i,j])/ETAPHI[i,j]+ (RX[i,j-1]*qxD[i,j-1]^2+RX[i,j]*qxD[i,j]^2)/2+ (RY[i-1,j]*qyD[i-1,j]^2+RY[i,j]*qyD[i,j]^2)/2
            end
        end
        # test
        @test HS ≈ HS_ver rtol=1e-9
    end # testset "compute_shear_heating!()"

    @testset "compute_adiabatic_heating!()" begin
        HA = zeros(Ny1, Nx1)
        tk1 = rand(rgen, Ny1, Nx1)
        ALPHA = rand(rgen, Ny1, Nx1)
        ALPHAF = rand(rgen, Ny1, Nx1)
        PHI = rand(rgen, Ny1, Nx1)
        vx = rand(rgen, Ny1, Nx1)
        vy = rand(rgen, Ny1, Nx1)
        vxf = rand(rgen, Ny1, Nx1)
        vyf = rand(rgen, Ny1, Nx1)
        ps = rand(rgen, Ny1, Nx1)
        pf = rand(rgen, Ny1, Nx1)
        HA_ver = zero(HA)
        # compute adiabatic heating
        Erebus.compute_adiabatic_heating!(
            HA, tk1, ALPHA, ALPHAF, PHI, vx, vy, vxf, vyf, ps, pf)
        # verification, from HTM-planetary.m, line 1573ff
        for j=2:1:Nx
            for i=2:1:Ny
                # HA
                # Indirect calculation of dpdt
                # Average vy; vx; vxf; vyf
                VXP=(vx[i,j]+vx[i,j-1])/2
                VYP=(vy[i,j]+vy[i-1,j])/2
                VXFP=(vxf[i,j]+vxf[i,j-1])/2
                VYFP=(vyf[i,j]+vyf[i-1,j])/2
                # Evaluate DPsolid/Dt with upwind differences
                if VXP<0
                    dpsdx=(ps[i,j]-ps[i,j-1])/dx
                else
                    dpsdx=(ps[i,j+1]-ps[i,j])/dx
                end
                if VYP<0
                    dpsdy=(ps[i,j]-ps[i-1,j])/dy
                else
                    dpsdy=(ps[i+1,j]-ps[i,j])/dy
                end
                dpsdt=VXP*dpsdx+VYP*dpsdy
                # Evaluate DPfluid/Dt with upwind differences
                if VXFP>0
                    dpfdx=(pf[i,j]-pf[i,j-1])/dx
                else
                    dpfdx=(pf[i,j+1]-pf[i,j])/dx
                end
                if VYFP>0
                    dpfdy=(pf[i,j]-pf[i-1,j])/dy
                else
                    dpfdy=(pf[i+1,j]-pf[i,j])/dy
                end
                dpfdt=VXFP*dpsdx+VYFP*dpsdy
        #         # Direct calculation of dpdt
        #         dpsdt=(ps[i,j]-ps0[i,j])/dt
        #         dpfdt=(pf[i,j]-pf0[i,j])/dt
                # HA
                HA_ver[i,j]=(1-PHI[i,j])*tk1[i,j]*ALPHA[i,j]*dpsdt+ PHI[i,j]*tk1[i,j]*ALPHAF[i,j]*dpfdt
            end
        end
        # test
        @test HA ≈ HA_ver rtol=1e-9
    end # testset "compute_adiabatic_heating!()"

end

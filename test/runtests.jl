using HydrologyPlanetesimals
using Parameters
using StaticArrays
using Test

@testset verbose = true "HydrologyPlanetesimals.jl" begin

    @testset "staggered grid" begin
        xsize=140_000.0
        ysize=140_000.0
        rplanet=50_000.0
        rcrust=48_000.0
        Nx=141
        Ny=141
        Nxmc=4
        Nymc=4
        sp = HydrologyPlanetesimals.StaticParameters(
            xsize=xsize,
            ysize=ysize,
            rplanet=rplanet,
            rcrust=rcrust,
            Nx=Nx,
            Ny=Ny,
            Nxmc=Nxmc,
            Nymc=Nymc
        )
        function setup(sp)
            @unpack xsize, ysize,
                Nx, Ny,
                Nx1, Ny1,
                dx, dy,
                jmin_basic, jmax_basic,
                imin_basic, imax_basic,
                jmin_vx, jmax_vx,
                imin_vx, imax_vx,
                jmin_vy, jmax_vy,
                imin_vy, imax_vy,
                jmin_p, jmax_p,
                imin_p, imax_p,
                rhosolidm,
                rhofluidm,
                etasolidm,
                etasolidmm,
                etafluidm,
                etafluidmm,
                rhocpsolidm,
                rhocpfluidm,
                alphasolidm,
                alphafluidm,
                ksolidm,
                kfluidm,
                start_hrsolidm,
                start_hrfluidm,
                gggsolidm,
                frictsolidm,
                cohessolidm,
                tenssolidm,
                kphim0,
                etaphikoef,
                phim0,
                tmsilicate,
                tmiron,
                etamin,
                nplast,
                dtelastic,
                start_step,
                nsteps,
                start_time, 
                endtime,
                start_marknum = sp
            x = SVector{Nx, Float64}([j for j = 0:dx:xsize])
            y = SVector{Ny, Float64}([j for j = 0:dy:ysize])
            xvx = SVector{Ny1, Float64}([j for j = 0:dx:xsize+dy])
            yvx = SVector{Nx1, Float64}([i for i = -dy/2:dy:ysize+dy/2])
            xvy = SVector{Nx1, Float64}([j for j = -dx/2:dx:xsize+dx/2])
            yvy = SVector{Ny1, Float64}([i for i = 0:dy:ysize+dy])
            xp = SVector{Nx1, Float64}([j for j = -dx/2:dx:xsize+dx/2])
            yp = SVector{Ny1, Float64}([i for i = -dy/2:dy:ysize+dy/2])
            return x, y, xvx, yvx, xvy, yvy, xp, yp
        end
        x, y, xvx, yvx, xvy, yvy, xp, yp = setup(sp)
        # from madcph.m lines 24ff
        xsize=140000
        ysize=140000
        Nx=141
        Ny=141
        Nx1=Nx+1
        Ny1=Ny+1
        dx=xsize/(Nx-1)
        dy=ysize/(Ny-1)
        x_ver=0:dx:xsize
        y_ver=0:dy:ysize
        xvx_ver=0:dx:xsize+dy
        yvx_ver=-dy/2:dy:ysize+dy/2
        xvy_ver=-dx/2:dx:xsize+dx/2
        yvy_ver=0:dy:ysize+dy
        xp_ver=-dx/2:dx:xsize+dx/2
        yp_ver=-dy/2:dy:ysize+dy/2
        # tests
        @test x == collect(x_ver)
        @test y == collect(y_ver)
        @test xvx == collect(xvx_ver)
        @test yvx == collect(yvx_ver)
        @test xvy == collect(xvy_ver)
        @test yvy == collect(yvy_ver)
        @test xp == collect(xp_ver)
        @test yp == collect(yp_ver)
    end # testset "staggered grid"

    @testset "distance()" begin
        @test HydrologyPlanetesimals.distance(0, 0, 0, 0) == 0
        @test HydrologyPlanetesimals.distance(1, 0, 0, 0) == 1
        @test HydrologyPlanetesimals.distance(0, 1, 0, 0) == 1
        @test HydrologyPlanetesimals.distance(0, 0, 1, 0) == 1
        @test HydrologyPlanetesimals.distance(0, 0, 0, 1) == 1
        @test HydrologyPlanetesimals.distance(0, 0, 1, 1) ≈ sqrt(2)
        @test HydrologyPlanetesimals.distance(1, 1, 0, 0) ≈ sqrt(2)
        @test HydrologyPlanetesimals.distance(-1, -1, 1, 1) ≈ sqrt(8)
        # random tests
        num_samples = 100
        x = rand(0:0.001:1000, 2, num_samples)
        y = rand(0:0.001:1000, 2, num_samples)
        for i in 1:num_samples
            @test HydrologyPlanetesimals.distance(
                x[:, i]..., y[:, i]...) ≈ sqrt(sum((x[:, i] - y[:, i]).^2))
        end
    end

    @testset "total()" begin
        @test HydrologyPlanetesimals.total(0, 0, 0) == 0
        @test HydrologyPlanetesimals.total(1, 0, 0) == 1
        @test HydrologyPlanetesimals.total(0, 1, 0) == 0
        @test HydrologyPlanetesimals.total(0, 0, 1) == 0
        @test HydrologyPlanetesimals.total(1, 2, 0.5) == 1.5
    end

    @testset "ktotal()" begin
        # from madcph.m, line 1761
        ktotalm(ksolidm, kfluid, phim)=(ksolidm*kfluid/2+((ksolidm*(3*phim-2)+kfluid*(1-3*phim))^2)/16)^0.5-(ksolidm*(3*phim-2)+ kfluid*(1-3*phim))/4
        @test HydrologyPlanetesimals.ktotal(1., 2., 3.) == ktotalm(1., 2., 3.)
    end

    @testset "kphi()" begin
        # from madcph.m, line 333
        kphim(kphim0, phim, phim0)=kphim0*(phim/phim0)^3/((1-phim)/(1-phim0))^2
        @test HydrologyPlanetesimals.kphi(1., 2., 3.) == kphim(1., 2., 3.)
    end

    @testset "Q_radiogenic()" begin
        # from madcph.m, line 276
        Q(f, ratio, E, tau, timesum)=f*ratio*E*exp(-timesum/tau)/tau
        @test HydrologyPlanetesimals.Q_radiogenic(1., 2., 3., 4., 5.) == Q(
            1., 2., 3., 4., 5.)
        @test HydrologyPlanetesimals.Q_radiogenic(1., 2., 3., 4., 0.) == Q(
            1., 2., 3., 4., 0.)
    end

    @testset "calculate_radioactive_heating()" begin
        hr_al = false
        hr_fe = false
        v = @SVector [0., 0., 0.]
        sp = HydrologyPlanetesimals.StaticParameters(hr_al=hr_al, hr_fe=hr_fe)
        @test HydrologyPlanetesimals.calculate_radioactive_heating(
            1000., sp) == (v, v)

        hr_al = true
        hr_fe = false        
        sp2 = HydrologyPlanetesimals.StaticParameters(hr_al=hr_al, hr_fe=hr_fe)
        Q_al = HydrologyPlanetesimals.Q_radiogenic(
            sp2.f_al, sp2.ratio_al, sp2.E_al, sp2.tau_al, 1000.)
        u = @SVector [Q_al * sp2.rhosolidm[1], Q_al * sp2.rhosolidm[2], 0.]
        @test HydrologyPlanetesimals.calculate_radioactive_heating(
            1000., sp2) == (u, v)

        hr_al = false
        hr_fe = true
        sp3 = HydrologyPlanetesimals.StaticParameters(hr_al=hr_al, hr_fe=hr_fe)
        Q_fe = HydrologyPlanetesimals.Q_radiogenic(
            sp3.f_fe, sp3.ratio_fe, sp3.E_fe, sp3.tau_fe, 1000.)
        w = @SVector [Q_fe * sp2.rhofluidm[1], 0., 0.]
        @test HydrologyPlanetesimals.calculate_radioactive_heating(
            1000., sp3) == (v, w)
    end

    @testset "fix_weights() elementary" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        dx, dy = sp.dx, sp.dy
        xsize, ysize = sp.xsize, sp.ysize
        # from madcph.m, line 38ff
        # Basic nodes
        x=0:dx:xsize
        y=0:dy:ysize
        # Vx-Nodes
        xvx=0:dx:xsize+dy
        yvx=-dy/2:dy:ysize+dy/2
        # Vy-nodes
        xvy=-dx/2:dx:xsize+dx/2
        yvy=0:dy:ysize+dy
        # P-Nodes
        xp=-dx/2:dx:xsize+dx/2
        yp=-dy/2:dy:ysize+dy/2

        @testset "basic nodes" begin
        # from madcph.m, line 373ff
        jmin, jmax = sp.jmin_basic, sp.jmax_basic
        imin, imax = sp.imin_basic, sp.imax_basic
        function fix_basic(xm, ym, x_axis, y_axis, dx, dy)
            j=trunc(Int, (xm-x_axis[1])/dx)+1;
            i=trunc(Int, (ym-y_axis[1])/dy)+1;
            if j<1
                j=1
            elseif j>Nx-1 
                j=Nx-1
            end
            if i<1 
                i=1
            elseif i>Ny-1
                i=Ny-1
            end
            dxmj=xm-x[j];
            dymi=ym-y[i];
            wtmij=(1-dxmj/dx)*(1-dymi/dy);
            wtmi1j=(1-dxmj/dx)*(dymi/dy);    
            wtmij1=(dxmj/dx)*(1-dymi/dy);
            wtmi1j1=(dxmj/dx)*(dymi/dy);
            return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
        end
        # top left
        xm = -x[1]
        ym = -y[1]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, x, y, dx, dy, jmin, jmax, imin, imax) == fix_basic(
            xm, ym, x, y, dx, dy)
        # bottom left
        xm = x[1]
        ym = y[end]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, x, y, dx, dy, jmin, jmax, imin, imax) == fix_basic(
            xm, ym, x, y, dx, dy)
        # top right
        xm = x[end]
        ym = y[1]
        j=trunc(Int, (xm-x[1])/dx)+1;
        i=trunc(Int, (ym-y[1])/dy)+1;
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, x, y, dx, dy, jmin, jmax, imin, imax) == fix_basic(
            xm, ym, x, y, dx, dy)
        # bottom right
        xm = x[end]
        ym = y[end]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, x, y, dx, dy, jmin, jmax, imin, imax) == fix_basic(
            xm, ym, x, y, dx, dy)
        end # testset "basic nodes"

        @testset "Vx nodes" begin
        # from madcph.m, line 434ff
        jmin, jmax = sp.jmin_vx, sp.jmax_vx
        imin, imax = sp.imin_vx, sp.imax_vx
        function fix_vx(xm, ym, x_axis, y_axis, dx, dy)
            j=trunc(Int, (xm-x_axis[1])/dx)+1;
            i=trunc(Int, (ym-y_axis[1])/dy)+1;
            if j<1
                j=1
            elseif j>Nx-1 
                j=Nx-1
            end
            if i<1 
                i=1
            elseif i>Ny
                i=Ny
            end
            dxmj=xm-xvx[j];
            dymi=ym-yvx[i];
            wtmij=(1-dxmj/dx)*(1-dymi/dy);
            wtmi1j=(1-dxmj/dx)*(dymi/dy);    
            wtmij1=(dxmj/dx)*(1-dymi/dy);
            wtmi1j1=(dxmj/dx)*(dymi/dy);
            return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
        end
        # top left
        xm = -xvx[1]
        ym = -yvx[1]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xvx, yvx, dx, dy, jmin, jmax, imin, imax) == fix_vx(
            xm, ym, xvx, yvx, dx, dy)
        # bottom left
        xm = xvx[1]
        ym = yvx[end]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xvx, yvx, dx, dy, jmin, jmax, imin, imax) == fix_vx(
            xm, ym, xvx, yvx, dx, dy)
        # top right
        xm = xvx[end]
        ym = yvx[1]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xvx, yvx, dx, dy, jmin, jmax, imin, imax) == fix_vx(
            xm, ym, xvx, yvx, dx, dy)
        # bottom right
        xm = xvx[end]
        ym = yvx[end]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xvx, yvx, dx, dy, jmin, jmax, imin, imax) == fix_vx(
            xm, ym, xvx, yvx, dx, dy)
        end # testset "Vx nodes"

        @testset "Vy nodes" begin
        # from madcph.m, line 484ff
        jmin, jmax = sp.jmin_vy, sp.jmax_vy
        imin, imax = sp.imin_vy, sp.imax_vy
        function fix_vy(xm, ym, x_axis, y_axis, dx, dy)
            j=trunc(Int, (xm-x_axis[1])/dx)+1;
            i=trunc(Int, (ym-y_axis[1])/dy)+1;
            if j<1
                j=1
            elseif j>Nx 
                j=Nx
            end
            if i<1 
                i=1
            elseif i>Ny-1
                i=Ny-1
            end
            dxmj=xm-xvy[j];
            dymi=ym-yvy[i];
            wtmij=(1-dxmj/dx)*(1-dymi/dy);
            wtmi1j=(1-dxmj/dx)*(dymi/dy);    
            wtmij1=(dxmj/dx)*(1-dymi/dy);
            wtmi1j1=(dxmj/dx)*(dymi/dy);
            return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
        end
        # top left
        xm = -xvy[1]
        ym = -yvy[1]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xvy, yvy, dx, dy, jmin, jmax, imin, imax) == fix_vy(
            xm, ym, xvy, yvy, dx, dy)
        # bottom left
        xm = xvy[1]
        ym = yvy[end]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xvy, yvy, dx, dy, jmin, jmax, imin, imax) == fix_vy(
            xm, ym, xvy, yvy, dx, dy)
        # top right
        xm = xvy[end]
        ym = yvy[1]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xvy, yvy, dx, dy, jmin, jmax, imin, imax) == fix_vy(
            xm, ym, xvy, yvy, dx, dy)
        # bottom right
        xm = xvy[end]
        ym = yvy[end]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xvy, yvy, dx, dy, jmin, jmax, imin, imax) == fix_vy(
            xm, ym, xvy, yvy, dx, dy)
        end # testset "Vy nodes"
    
        @testset "P nodes" begin
        # from madcph.m, line 538ff
        jmin, jmax = sp.jmin_p, sp.jmax_p
        imin, imax = sp.imin_p, sp.imax_p
        function fix_p(xm, ym, x_axis, y_axis, dx, dy)
            j=trunc(Int, (xm-x_axis[1])/dx)+1;
            i=trunc(Int, (ym-y_axis[1])/dy)+1;
            if j<1
                j=1
            elseif j>Nx 
                j=Nx
            end
            if i<1 
                i=1
            elseif i>Ny
                i=Ny
            end
            dxmj=xm-xp[j];
            dymi=ym-yp[i];
            wtmij=(1-dxmj/dx)*(1-dymi/dy);
            wtmi1j=(1-dxmj/dx)*(dymi/dy);    
            wtmij1=(dxmj/dx)*(1-dymi/dy);
            wtmi1j1=(dxmj/dx)*(dymi/dy);
            return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
        end
        # top left
        xm = -xp[1]
        ym = -yp[1]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xp, yp, dx, dy, jmin, jmax, imin, imax) == fix_p(
            xm, ym, xp, yp, dx, dy)
        # bottom left
        xm = xp[1]
        ym = yp[end]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xp, yp, dx, dy, jmin, jmax, imin, imax) == fix_p(
            xm, ym, xp, yp, dx, dy)
        # top right
        xm = xp[end]
        ym = yp[1]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xp, yp, dx, dy, jmin, jmax, imin, imax) == fix_p(
            xm, ym, xp, yp, dx, dy)
        # bottom right
        xm = xp[end]
        ym = yp[end]
        @test HydrologyPlanetesimals.fix_weights(
            xm, ym, xp, yp, dx, dy, jmin, jmax, imin, imax) == fix_p(
            xm, ym, xp, yp, dx, dy)
        end # testset "P nodes"    
    end # testset "fix_weights() elementary"

    @testset "fix_weights() advanced" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        dx, dy = sp.dx, sp.dy
        xsize, ysize = sp.xsize, sp.ysize
        # from madcph.m, line 38ff
        # basic nodes
        x=0:dx:xsize
        y=0:dy:ysize
        # Vx nodes
        xvx=0:dx:xsize+dy
        yvx=-dy/2:dy:ysize+dy/2
        # Vy nodes
        xvy=-dx/2:dx:xsize+dx/2
        yvy=0:dy:ysize+dy
        # P Nodes
        xp=-dx/2:dx:xsize+dx/2
        yp=-dy/2:dy:ysize+dy/2
        # simulating markers
        num_markers = 10_000
        @testset "basic nodes" begin
            # from madcph.m, line 373ff
            jmin, jmax = sp.jmin_basic, sp.jmax_basic
            imin, imax = sp.imin_basic, sp.imax_basic
            function fix_basic(xm, ym, x_axis, y_axis, dx, dy)
                j=trunc(Int, (xm-x_axis[1])/dx)+1;
                i=trunc(Int, (ym-y_axis[1])/dy)+1;
                if j<1
                    j=1
                elseif j>Nx-1 
                    j=Nx-1
                end
                if i<1 
                    i=1
                elseif i>Ny-1
                    i=Ny-1
                end
                dxmj=xm-x[j];
                dymi=ym-y[i];
                wtmij=(1-dxmj/dx)*(1-dymi/dy);
                wtmi1j=(1-dxmj/dx)*(dymi/dy);    
                wtmij1=(dxmj/dx)*(1-dymi/dy);
                wtmi1j1=(dxmj/dx)*(dymi/dy);
                return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
            end
            # simulating markers
            xm = rand(-x[1]:0.1:x[end]+dx, num_markers)
            ym = rand(-y[1]:0.1:y[end]+dy, num_markers)
            for m in 1:num_markers
                i, j, weights = HydrologyPlanetesimals.fix_weights(
                    xm[m],
                    ym[m],
                    x,
                    y,
                    dx,
                    dy,
                    jmin,
                    jmax,
                    imin,
                    imax
                )
                i_ver, j_ver, weights_ver = fix_basic(
                    xm[m], ym[m], x, y, dx, dy)
                @debug "fix_weights basic" i i_ver j j_ver weights weights_ver
                @test i == i_ver
                @test j == j_ver
                @test weights == weights_ver
            end
        end # testset "basic nodes"
        
        @testset "Vx nodes" begin
            # from madcph.m, line 434ff
            jmin, jmax = sp.jmin_vx, sp.jmax_vx
            imin, imax = sp.imin_vx, sp.imax_vx
            function fix_vx(xm, ym, x_axis, y_axis, dx, dy)
                j=trunc(Int, (xm-x_axis[1])/dx)+1;
                i=trunc(Int, (ym-y_axis[1])/dy)+1;
                if j<1
                    j=1
                elseif j>Nx-1 
                    j=Nx-1
                end
                if i<1 
                    i=1
                elseif i>Ny
                    i=Ny
                end
                dxmj=xm-xvx[j];
                dymi=ym-yvx[i];
                wtmij=(1-dxmj/dx)*(1-dymi/dy);
                wtmi1j=(1-dxmj/dx)*(dymi/dy);    
                wtmij1=(dxmj/dx)*(1-dymi/dy);
                wtmi1j1=(dxmj/dx)*(dymi/dy);
                return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
            end
            # simulating markers
            xm = rand(-xvx[1]:0.1:xvx[end]+dx, num_markers)
            ym = rand(-yvx[1]:0.1:yvx[end]+dy, num_markers)
            for m in 1:num_markers
                i, j, weights = HydrologyPlanetesimals.fix_weights(
                    xm[m],
                    ym[m],
                    xvx,
                    yvx,
                    dx,
                    dy,
                    jmin,
                    jmax,
                    imin,
                    imax
                )
                i_ver, j_ver, weights_ver = fix_vx(
                    xm[m], ym[m], xvx, yvx, dx, dy)
                @debug "fix_weights Vx" i i_ver j j_ver weights weights_ver
                @test i == i_ver
                @test j == j_ver
                @test weights == weights_ver
            end
        end # testset "Vx nodes"

        @testset "Vy nodes" begin
            # from madcph.m, line 484ff
            jmin, jmax = sp.jmin_vy, sp.jmax_vy
            imin, imax = sp.imin_vy, sp.imax_vy
            function fix_vy(xm, ym, x_axis, y_axis, dx, dy)
                j=trunc(Int, (xm-x_axis[1])/dx)+1;
                i=trunc(Int, (ym-y_axis[1])/dy)+1;
                if j<1
                    j=1
                elseif j>Nx 
                    j=Nx
                end
                if i<1 
                    i=1
                elseif i>Ny-1
                    i=Ny-1
                end
                dxmj=xm-xvy[j];
                dymi=ym-yvy[i];
                wtmij=(1-dxmj/dx)*(1-dymi/dy);
                wtmi1j=(1-dxmj/dx)*(dymi/dy);    
                wtmij1=(dxmj/dx)*(1-dymi/dy);
                wtmi1j1=(dxmj/dx)*(dymi/dy);
                return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
            end
            # simulating markers
            xm = rand(-xvy[1]:0.1:xvy[end]+dx, num_markers)
            ym = rand(-yvy[1]:0.1:yvy[end]+dy, num_markers)
            for m in 1:num_markers
                i, j, weights = HydrologyPlanetesimals.fix_weights(
                    xm[m],
                    ym[m],
                    xvy,
                    yvy,
                    dx,
                    dy,
                    jmin,
                    jmax,
                    imin,
                    imax
                )
                i_ver, j_ver, weights_ver = fix_vy(
                    xm[m], ym[m], xvy, yvy, dx, dy)
                @debug "fix_weights Vy" i i_ver j j_ver weights weights_ver
                @test i == i_ver
                @test j == j_ver
                @test weights == weights_ver
            end
        end # testset "Vy nodes"
    
        @testset "P nodes" begin
            # from madcph.m, line 538ff
            jmin, jmax = sp.jmin_p, sp.jmax_p
            imin, imax = sp.imin_p, sp.imax_p
            function fix_p(xm, ym, x_axis, y_axis, dx, dy)
                j=trunc(Int, (xm-x_axis[1])/dx)+1;
                i=trunc(Int, (ym-y_axis[1])/dy)+1;
                if j<1
                    j=1
                elseif j>Nx 
                    j=Nx
                end
                if i<1 
                    i=1
                elseif i>Ny
                    i=Ny
                end
                dxmj=xm-xp[j];
                dymi=ym-yp[i];
                wtmij=(1-dxmj/dx)*(1-dymi/dy);
                wtmi1j=(1-dxmj/dx)*(dymi/dy);    
                wtmij1=(dxmj/dx)*(1-dymi/dy);
                wtmi1j1=(dxmj/dx)*(dymi/dy);
                return i, j, @SVector [wtmij, wtmi1j, wtmij1, wtmi1j1]
            end
            # simulating markers
            xm = rand(-xp[1]:0.1:xp[end]+dx, num_markers)
            ym = rand(-yp[1]:0.1:yp[end]+dy, num_markers)
            for m in 1:num_markers
                i, j, weights = HydrologyPlanetesimals.fix_weights(
                    xm[m],
                    ym[m],
                    xp,
                    yp,
                    dx,
                    dy,
                    jmin,
                    jmax,
                    imin,
                    imax
                )
                i_ver, j_ver, weights_ver = fix_p(
                    xm[m], ym[m], xp, yp, dx, dy)
                @debug "fix_weights P" xm[m] ym[m] i i_ver j j_ver
                @test i == i_ver
                @test j == j_ver
                @test weights == weights_ver
            end
        end # testset "P nodes"    
    end # testset "fix_weights() advanced"
   
    @testset "interpolate!(i,j,weights,property,grid)" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        dx, dy = sp.dx, sp.dy
        xsize, ysize = sp.xsize, sp.ysize
        jmin, jmax = sp.jmin_basic, sp.jmax_basic
        imin, imax = sp.imin_basic, sp.imax_basic
        x=0:dx:xsize
        y=0:dy:ysize
        # simulate markers
        num_markers = 10_000
        xm = rand(-x[1]:0.1:x[end]+dx, num_markers)
        ym = rand(-y[1]:0.1:y[end]+dy, num_markers)
        property = rand(num_markers)
        # sample interpolation array
        for m=1:1:num_markers
            grid = zeros(Ny, Nx, Base.Threads.nthreads())
            i, j, weights = HydrologyPlanetesimals.fix_weights(
                xm[m], ym[m], x, y, dx, dy, jmin, jmax, imin, imax)
            HydrologyPlanetesimals.interpolate!(
                i, j, weights, property[m], grid)
            @test grid[i, j] == property[m] * weights[1] 
            @test grid[i+1, j] == property[m] * weights[2]
            @test grid[i, j+1] == property[m] * weights[3]
            @test grid[i+1, j+1] == property[m] * weights[4]
        end
    end # testset "interpolate!"

    @testset "compute node properties" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        Nx1, Ny1 = sp.Nx1, sp.Ny1
        dx, dy = sp.dx, sp.dy
        xsize, ysize = sp.xsize, sp.ysize
        # from madcph.m, line 38ff
        # basic nodes
        x=0:dx:xsize
        y=0:dy:ysize
        # Vx nodes
        xvx=0:dx:xsize+dy
        yvx=-dy/2:dy:ysize+dy/2
        # Vy nodes
        xvy=-dx/2:dx:xsize+dx/2
        yvy=0:dy:ysize+dy
        # P Nodes
        xp=-dx/2:dx:xsize+dx/2
        yp=-dy/2:dy:ysize+dy/2
        # simulating markers
        num_markers = 10_000
        @testset "compute_basic_node_properties!()" begin    
            jmin, jmax = sp.jmin_basic, sp.jmax_basic
            imin, imax = sp.imin_basic, sp.imax_basic
            ETA0SUM = zeros(Ny, Nx, Base.Threads.nthreads())
            ETASUM = zeros(Ny, Nx, Base.Threads.nthreads())
            GGGSUM = zeros(Ny, Nx, Base.Threads.nthreads())
            SXYSUM = zeros(Ny, Nx, Base.Threads.nthreads())
            COHSUM = zeros(Ny, Nx, Base.Threads.nthreads())
            TENSUM = zeros(Ny, Nx, Base.Threads.nthreads())
            FRISUM = zeros(Ny, Nx, Base.Threads.nthreads())
            WTSUM = zeros(Ny, Nx, Base.Threads.nthreads())
            ETA0 = zeros(Float64, Ny, Nx)
            ETA = zeros(Float64, Ny, Nx)
            GGG = zeros(Float64, Ny, Nx)
            SXY0 = zeros(Float64, Ny, Nx)
            COH = zeros(Float64, Ny, Nx)
            TEN = zeros(Float64, Ny, Nx)
            FRI = zeros(Float64, Ny, Nx)
            YNY = zeros(Bool, Ny, Nx)
            ETA0SUM_ver = zeros(Ny, Nx, Base.Threads.nthreads())
            ETASUM_ver = zeros(Ny, Nx, Base.Threads.nthreads())
            GGGSUM_ver = zeros(Ny, Nx, Base.Threads.nthreads())
            SXYSUM_ver = zeros(Ny, Nx, Base.Threads.nthreads())
            COHSUM_ver = zeros(Ny, Nx, Base.Threads.nthreads())
            TENSUM_ver = zeros(Ny, Nx, Base.Threads.nthreads())
            FRISUM_ver = zeros(Ny, Nx, Base.Threads.nthreads())
            WTSUM_ver = zeros(Ny, Nx, Base.Threads.nthreads())
            ETA0_ver = zeros(Float64, Ny, Nx)
            ETA_ver = zeros(Float64, Ny, Nx)
            GGG_ver = zeros(Float64, Ny, Nx)
            SXY0_ver = zeros(Float64, Ny, Nx)
            COH_ver = zeros(Float64, Ny, Nx)
            TEN_ver = zeros(Float64, Ny, Nx)
            FRI_ver = zeros(Float64, Ny, Nx)
            YNY_ver = zeros(Bool, Ny, Nx)
            # simulate markers
            xm = rand(-x[1]:0.1:x[end]+dx, num_markers)
            ym = rand(-y[1]:0.1:y[end]+dy, num_markers)
            property = rand(7, num_markers)
            # calculate grid properties
            for m=1:1:num_markers
                i, j, weights = HydrologyPlanetesimals.fix_weights(
                    xm[m], ym[m], x, y, dx, dy, jmin, jmax, imin, imax)
                HydrologyPlanetesimals.interpolate!(
                    i, j, weights, property[1, m], ETA0SUM)
                HydrologyPlanetesimals.interpolate!(
                    i, j, weights, property[2, m], ETASUM)
                HydrologyPlanetesimals.interpolate!(
                    i, j, weights, inv(property[3, m]), GGGSUM)
                HydrologyPlanetesimals.interpolate!(
                    i, j, weights, property[4, m], SXYSUM)
                HydrologyPlanetesimals.interpolate!(
                    i, j, weights, property[5, m], COHSUM)
                HydrologyPlanetesimals.interpolate!(
                    i, j, weights, property[6, m], TENSUM)
                HydrologyPlanetesimals.interpolate!(
                    i, j, weights, property[7, m], FRISUM)
                HydrologyPlanetesimals.interpolate!(
                    i, j, weights, 1.0, WTSUM)
            end
            HydrologyPlanetesimals.compute_basic_node_properties!(
                ETA0SUM,
                ETASUM,
                GGGSUM,
                SXYSUM,
                COHSUM,
                TENSUM,
                FRISUM,
                WTSUM,
                ETA0,
                ETA,
                GGG,
                SXY0,
                COH,
                TEN,
                FRI,
                YNY
            )
            # verification properties, from madcph.m, lines 373ff, 606ff
            for m=1:1:num_markers
                j=trunc(Int, (xm[m]-x[1])/dx)+1
                i=trunc(Int, (ym[m]-y[1])/dy)+1
                if j<1
                    j=1
                elseif j>Nx-1
                    j=Nx-1
                end
                if i<1
                    i=1
                elseif i>Ny-1 
                    i=Ny-1
                end
                # Compute distances
                dxmj=xm[m]-x[j]
                dymi=ym[m]-y[i]
                # Compute weights
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy);    
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Update properties
                # i;j Node
                ETA0SUM_ver[i,j]=ETA0SUM_ver[i,j]+property[1,m]*wtmij
                ETASUM_ver[i,j]=ETASUM_ver[i,j]+property[2,m]*wtmij
                GGGSUM_ver[i,j]=GGGSUM_ver[i,j]+1/property[3,m]*wtmij
                SXYSUM_ver[i,j]=SXYSUM_ver[i,j]+property[4,m]*wtmij
                COHSUM_ver[i,j]=COHSUM_ver[i,j]+property[5,m]*wtmij
                TENSUM_ver[i,j]=TENSUM_ver[i,j]+property[6,m]*wtmij
                FRISUM_ver[i,j]=FRISUM_ver[i,j]+property[7,m]*wtmij
                WTSUM_ver[i,j]=WTSUM_ver[i,j]+wtmij
                # i+1;j Node
                ETA0SUM_ver[i+1,j]=ETA0SUM_ver[i+1,j]+property[1,m]*wtmi1j
                ETASUM_ver[i+1,j]=ETASUM_ver[i+1,j]+property[2,m]*wtmi1j
                GGGSUM_ver[i+1,j]=GGGSUM_ver[i+1,j]+1/property[3,m]*wtmi1j
                SXYSUM_ver[i+1,j]=SXYSUM_ver[i+1,j]+property[4,m]*wtmi1j
                COHSUM_ver[i+1,j]=COHSUM_ver[i+1,j]+property[5,m]*wtmi1j
                TENSUM_ver[i+1,j]=TENSUM_ver[i+1,j]+property[6,m]*wtmi1j
                FRISUM_ver[i+1,j]=FRISUM_ver[i+1,j]+property[7,m]*wtmi1j
                WTSUM_ver[i+1,j]=WTSUM_ver[i+1,j]+wtmi1j
                # i;j+1 Node
                ETA0SUM_ver[i,j+1]=ETA0SUM_ver[i,j+1]+property[1,m]*wtmij1
                ETASUM_ver[i,j+1]=ETASUM_ver[i,j+1]+property[2,m]*wtmij1
                GGGSUM_ver[i,j+1]=GGGSUM_ver[i,j+1]+1/property[3,m]*wtmij1
                SXYSUM_ver[i,j+1]=SXYSUM_ver[i,j+1]+property[4,m]*wtmij1
                COHSUM_ver[i,j+1]=COHSUM_ver[i,j+1]+property[5,m]*wtmij1
                TENSUM_ver[i,j+1]=TENSUM_ver[i,j+1]+property[6,m]*wtmij1
                FRISUM_ver[i,j+1]=FRISUM_ver[i,j+1]+property[7,m]*wtmij1
                WTSUM_ver[i,j+1]=WTSUM_ver[i,j+1]+wtmij1
                # i+1;j+1 Node
                ETA0SUM_ver[i+1,j+1]=ETA0SUM_ver[i+1,j+1]+property[1,m]*wtmi1j1
                ETASUM_ver[i+1,j+1]=ETASUM_ver[i+1,j+1]+property[2,m]*wtmi1j1
                GGGSUM_ver[i+1,j+1]=GGGSUM_ver[i+1,j+1]+1/property[3,m]*wtmi1j1
                SXYSUM_ver[i+1,j+1]=SXYSUM_ver[i+1,j+1]+property[4,m]*wtmi1j1
                COHSUM_ver[i+1,j+1]=COHSUM_ver[i+1,j+1]+property[5,m]*wtmi1j1
                TENSUM_ver[i+1,j+1]=TENSUM_ver[i+1,j+1]+property[6,m]*wtmi1j1
                FRISUM_ver[i+1,j+1]=FRISUM_ver[i+1,j+1]+property[7,m]*wtmi1j1
                WTSUM_ver[i+1,j+1]=WTSUM_ver[i+1,j+1]+wtmi1j1
            end
            for j=1:1:Nx
                for i=1:1:Ny
                    if WTSUM_ver[i,j]>0 
                        ETA0_ver[i,j]=ETA0SUM_ver[i,j]/WTSUM_ver[i,j]
                        ETA_ver[i,j]=ETASUM_ver[i,j]/WTSUM_ver[i,j]
                        if(ETA_ver[i,j]<ETA0_ver[i,j])
                            YNY_ver[i,j]=1
                        end
                        GGG_ver[i,j]=1/(GGGSUM_ver[i,j]/WTSUM_ver[i,j])
                        SXY0_ver[i,j]=SXYSUM_ver[i,j]/WTSUM_ver[i,j]
                        COH_ver[i,j]=COHSUM_ver[i,j]/WTSUM_ver[i,j]
                        TEN_ver[i,j]=TENSUM_ver[i,j]/WTSUM_ver[i,j]
                        FRI_ver[i,j]=FRISUM_ver[i,j]/WTSUM_ver[i,j]
                    end
                end
            end 
            # test
            for j=1:1:Nx, i=1:1:Ny
                @test ETA0[i, j] == ETA0_ver[i, j]
                @test ETA[i, j] == ETA_ver[i, j]
                @test GGG[i, j] ≈ GGG_ver[i, j]
                @test SXY0[i, j] == SXY0_ver[i, j]
                @test COH[i, j] == COH_ver[i, j]
                @test TEN[i, j] == TEN_ver[i, j]
                @test FRI[i, j] == FRI_ver[i, j]
                @test YNY[i, j] == YNY_ver[i, j]
            end
        end # testset "compute_basic_node_properties!"

        @testset "compute_vx_node_properties!()" begin
            jmin, jmax = sp.jmin_vx, sp.jmax_vx
            imin, imax = sp.imin_vx, sp.imax_vx 
            RHOXSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            RHOFXSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            KXSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            PHIXSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            RXSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            WTXSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
            RHOX = zeros(Float64, Ny1, Nx1)
            RHOFX = zeros(Float64, Ny1, Nx1)
            KX = zeros(Float64, Ny1, Nx1)
            PHIX = zeros(Float64, Ny1, Nx1)
            RX = zeros(Float64, Ny1, Nx1)
            RHOXSUM_ver = zeros(Ny1, Nx1)
            RHOFXSUM_ver = zeros(Ny1, Nx1)
            KXSUM_ver = zeros(Ny1, Nx1)
            PHIXSUM_ver = zeros(Ny1, Nx1)
            RXSUM_ver = zeros(Ny1, Nx1)
            WTXSUM_ver = zeros(Ny1, Nx1)
            RHOX_ver = zeros(Float64, Ny1, Nx1)
            RHOFX_ver = zeros(Float64, Ny1, Nx1)
            KX_ver = zeros(Float64, Ny1, Nx1)
            PHIX_ver = zeros(Float64, Ny1, Nx1)
            RX_ver = zeros(Float64, Ny1, Nx1)
            # simulate markers
            xm = rand(-xvx[1]:0.1:xvx[end]+dx, num_markers)
            ym = rand(-yvx[1]:0.1:yvx[end]+dy, num_markers)
            property = rand(5, num_markers)
            # calculate grid properties
            for m=1:1:num_markers
                i, j, weights = HydrologyPlanetesimals.fix_weights(
                    xm[m], ym[m], xvx, yvx, dx, dy, jmin, jmax, imin, imax)
                HydrologyPlanetesimals.interpolate!(
                    i, j, weights, property[1, m], RHOXSUM)
                HydrologyPlanetesimals.interpolate!(
                    i, j, weights, property[2, m], RHOFXSUM)
                HydrologyPlanetesimals.interpolate!(
                    i, j, weights, property[3, m], KXSUM)
                HydrologyPlanetesimals.interpolate!(
                    i, j, weights, property[4, m], PHIXSUM)
                HydrologyPlanetesimals.interpolate!(
                    i, j, weights, property[5, m], RXSUM)
                HydrologyPlanetesimals.interpolate!(
                    i, j, weights, 1.0, WTXSUM)
            end
            HydrologyPlanetesimals.compute_vx_node_properties!(
                RHOXSUM,
                RHOFXSUM,
                KXSUM,
                PHIXSUM,
                RXSUM,
                WTXSUM,
                RHOX,
                RHOFX,
                KX,
                PHIX,
                RX
            )
            # verification properties, from madcph.m, lines 434ff, 624ff
            for m=1:1:num_markers
                j=trunc(Int, (xm[m]-xvx[1])/dx)+1
                i=trunc(Int, (ym[m]-yvx[1])/dy)+1
                if j<1
                    j=1
                elseif j>Nx-1
                    j=Nx-1
                end
                if(i<1)
                    i=1
                elseif i>Ny
                    i=Ny
                end
                # Compute distances
                dxmj=xm[m]-xvx[j]
                dymi=ym[m]-yvx[i]
                # Compute weights
                wtmij=(1-dxmj/dx)*(1-dymi/dy)
                wtmi1j=(1-dxmj/dx)*(dymi/dy);    
                wtmij1=(dxmj/dx)*(1-dymi/dy)
                wtmi1j1=(dxmj/dx)*(dymi/dy)
                # Update properties
                # i;j Node
                RHOXSUM_ver[i,j]=RHOXSUM_ver[i,j]+property[1, m]*wtmij
                RHOFXSUM_ver[i,j]=RHOFXSUM_ver[i,j]+property[2, m]*wtmij
                KXSUM_ver[i,j]=KXSUM_ver[i,j]+property[3, m]*wtmij
                PHIXSUM_ver[i,j]=PHIXSUM_ver[i,j]+property[4, m]*wtmij
                RXSUM_ver[i,j]=RXSUM_ver[i,j]+property[5, m]*wtmij
                WTXSUM_ver[i,j]=WTXSUM_ver[i,j]+wtmij
                # i+1;j Node
                RHOXSUM_ver[i+1,j]=RHOXSUM_ver[i+1,j]+property[1, m]*wtmi1j
                RHOFXSUM_ver[i+1,j]=RHOFXSUM_ver[i+1,j]+property[2, m]*wtmi1j
                KXSUM_ver[i+1,j]=KXSUM_ver[i+1,j]+property[3, m]*wtmi1j
                PHIXSUM_ver[i+1,j]=PHIXSUM_ver[i+1,j]+property[4, m]*wtmi1j
                RXSUM_ver[i+1,j]=RXSUM_ver[i+1,j]+property[5, m]*wtmi1j
                WTXSUM_ver[i+1,j]=WTXSUM_ver[i+1,j]+wtmi1j
                # i;j+1 Node
                RHOXSUM_ver[i,j+1]=RHOXSUM_ver[i,j+1]+property[1, m]*wtmij1
                RHOFXSUM_ver[i,j+1]=RHOFXSUM_ver[i,j+1]+property[2, m]*wtmij1
                KXSUM_ver[i,j+1]=KXSUM_ver[i,j+1]+property[3, m]*wtmij1
                PHIXSUM_ver[i,j+1]=PHIXSUM_ver[i,j+1]+property[4, m]*wtmij1
                RXSUM_ver[i,j+1]=RXSUM_ver[i,j+1]+property[5, m]*wtmij1
                WTXSUM_ver[i,j+1]=WTXSUM_ver[i,j+1]+wtmij1
                # i+1;j+1 Node
                RHOXSUM_ver[i+1,j+1]=RHOXSUM_ver[i+1,j+1]+property[1, m]*wtmi1j1
                RHOFXSUM_ver[i+1,j+1]=RHOFXSUM_ver[i+1,j+1]+property[2, m]*wtmi1j1
                KXSUM_ver[i+1,j+1]=KXSUM_ver[i+1,j+1]+property[3, m]*wtmi1j1
                PHIXSUM_ver[i+1,j+1]=PHIXSUM_ver[i+1,j+1]+property[4, m]*wtmi1j1
                RXSUM_ver[i+1,j+1]=RXSUM_ver[i+1,j+1]+property[5, m]*wtmi1j1
                WTXSUM_ver[i+1,j+1]=WTXSUM_ver[i+1,j+1]+wtmi1j1
            end
            for j=1:1:Nx1
                for i=1:1:Ny1
                    if(WTXSUM_ver[i,j]>0)
                        RHOX_ver[i,j]=RHOXSUM_ver[i,j]/WTXSUM_ver[i,j]
                        RHOFX_ver[i,j]=RHOFXSUM_ver[i,j]/WTXSUM_ver[i,j]
                        KX_ver[i,j]=KXSUM_ver[i,j]/WTXSUM_ver[i,j]
                        PHIX_ver[i,j]=PHIXSUM_ver[i,j]/WTXSUM_ver[i,j]
                        RX_ver[i,j]=RXSUM_ver[i,j]/WTXSUM_ver[i,j]
                    end
                end
            end
            # test
            for j=1:1:Nx1, i=1:1:Ny1
                @test RHOX[i, j] == RHOX_ver[i, j]
                @test RHOFX[i, j] == RHOFX_ver[i, j]
                @test KX[i, j] == KX_ver[i, j]
                @test PHIX[i, j] == PHIX_ver[i, j]
                @test RX[i, j] == RX_ver[i, j]
            end
        end # testset "compute_vx_node_properties!"

        @testset "compute_vy_node_properties!()" begin
            
        end # testset "compute_vy_node_properties!"

        @testset "compute_p_node_properties!()" begin
            
        end # testset "compute_p_node_properties!
    end # testset "compute node properties" 
end

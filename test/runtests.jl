using HydrologyPlanetesimals
using Test
using StaticArrays

@testset verbose = true "HydrologyPlanetesimals.jl" begin

    @testset "distance(x1,y1,x2,y2)" begin
        @test HydrologyPlanetesimals.distance(0, 0, 0, 0) == 0
        @test HydrologyPlanetesimals.distance(1, 0, 0, 0) == 1
        @test HydrologyPlanetesimals.distance(0, 1, 0, 0) == 1
        @test HydrologyPlanetesimals.distance(0, 0, 1, 0) == 1
        @test HydrologyPlanetesimals.distance(0, 0, 0, 1) == 1
        @test HydrologyPlanetesimals.distance(0, 0, 1, 1) ≈ sqrt(2)
        @test HydrologyPlanetesimals.distance(1, 1, 0, 0) ≈ sqrt(2)
        @test HydrologyPlanetesimals.distance(-1, -1, 1, 1) ≈ sqrt(8)
    end

    @testset "total(solid,fluid,phi)" begin
        @test HydrologyPlanetesimals.total(0, 0, 0) == 0
        @test HydrologyPlanetesimals.total(1, 0, 0) == 1
        @test HydrologyPlanetesimals.total(0, 1, 0) == 0
        @test HydrologyPlanetesimals.total(0, 0, 1) == 0
        @test HydrologyPlanetesimals.total(1, 2, 0.5) == 1.5
    end

    @testset "ktotal(ksolid,kfluid,phi)" begin
        # from madcph.m, line 1761
        ktotalm(ksolidm, kfluid, phim)=(ksolidm*kfluid/2+((ksolidm*(3*phim-2)+kfluid*(1-3*phim))^2)/16)^0.5-(ksolidm*(3*phim-2)+ kfluid*(1-3*phim))/4
        @test HydrologyPlanetesimals.ktotal(1., 2., 3.) == ktotalm(1., 2., 3.)
    end

    @testset "kphi(kphim0,phim,phim0)" begin
        # from madcph.m, line 333
        kphim(kphim0, phim, phim0)=kphim0*(phim/phim0)^3/((1-phim)/(1-phim0))^2
        @test HydrologyPlanetesimals.kphi(1., 2., 3.) == kphim(1., 2., 3.)
    end

    @testset "Q_radiogenic(f,ratio,E,tau,time)" begin
        # from madcph.m, line 276
        Q(f, ratio, E, tau, timesum)=f*ratio*E*exp(-timesum/tau)/tau
        @test HydrologyPlanetesimals.Q_radiogenic(1., 2., 3., 4., 5.) == Q(
            1., 2., 3., 4., 5.)
        @test HydrologyPlanetesimals.Q_radiogenic(1., 2., 3., 4., 0.) == Q(
            1., 2., 3., 4., 0.)
    end

    @testset "calculate_radioactive_heating(timesum,sp)" begin
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

    @testset "fix_weights(x,y,x_axis,y_axis,dx,dy,jmin,jmax,imin,imax)" begin
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
    end # testset "fix_weights"
   
    @testset "interpolate!(i,j,weights,property,grid)" begin
        sp = HydrologyPlanetesimals.StaticParameters()
        Nx, Ny = sp.Nx, sp.Ny
        dx, dy = sp.dx, sp.dy
        xsize, ysize = sp.xsize, sp.ysize
        jmin, jmax = sp.jmin_basic, sp.jmax_basic
        imin, imax = sp.imin_basic, sp.imax_basic
        x=0:dx:xsize
        y=0:dy:ysize    
        xm = -x[1]
        ym = -y[1]
        # sample interpolation array
        grid = zeros(Ny, Nx, Base.Threads.nthreads())
        property_1 = rand()
        property_2 = rand()
        i, j, weights = HydrologyPlanetesimals.fix_weights(
            xm, ym, x, y, dx, dy, jmin, jmax, imin, imax)
        # interpolate 1
        HydrologyPlanetesimals.interpolate!(i, j, weights, property_1, grid)
        # check
        @test grid[i, j] == property_1 * weights[1]
        @test grid[i+1, j] == property_1 * weights[2]
        @test grid[i, j+1] == property_1 * weights[3]
        @test grid[i+1, j+1] == property_1 * weights[4]
        # interpolate 2
        HydrologyPlanetesimals.interpolate!(i, j, weights, property_2, grid)
        # check
        @test grid[i, j] == (property_1 + property_2)  * weights[1]
        @test grid[i+1, j] == (property_1 + property_2) * weights[2]
        @test grid[i, j+1] == (property_1 + property_2) * weights[3]
        @test grid[i+1, j+1] == (property_1 + property_2) * weights[4]
    end # testset "interpolate!"
end

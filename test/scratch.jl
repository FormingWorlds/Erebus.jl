using Base.Threads
# using SparseArrays
# using MAT
# using DocStringExtensions
using StaticArrays
using BenchmarkTools
# using TimerOutputs

using HydrologyPlanetesimals
run_simulation()


# xsize=140000.0
# ysize=140000.0
# rplanet=50000.0
# rcrust=48000.0
# Nx=141
# Ny=141
# Nxmc=4
# Nymc=4

# sp = HydrologyPlanetesimals.StaticParameters(
#     xsize=xsize,
#     ysize=ysize,
#     rplanet=rplanet,
#     rcrust=rcrust,
#     Nx=Nx,
#     Ny=Ny,
#     Nxmc=Nxmc,
#     Nymc=Nymc
#     )

function f(vx, vy, Nx, Ny, Nx1, Ny1)
    VXP = zeros(Ny1, Nx1)
    @views @. VXP[2:Ny, 2:Nx]= 0.5 * (vx[2:Ny, 2:Nx] + vx[2:Ny,1:Nx-1])
end

function f2(VXP, vx, vy, Nx, Ny, Nx1, Ny1)
    @views @. VXP[2:Ny, 2:Nx]= 0.5 * (vx[2:Ny, 2:Nx] + vx[2:Ny,1:Nx-1])
end


Nx, Ny = 140, 140
Nx1, Ny1 = 141, 141
vx, vy = rand(Ny1, Nx1), rand(Ny1, Nx1)
VXP = zeros(Ny1, Nx1)

@btime f($vx,$vy,$Nx,$Ny,$Nx1,$Ny1)
@btime f2($VXP,$vx,$vy,$Nx,$Ny,$Nx1,$Ny1)

function f_HA(HA, tk1, ALPHA, ALPHAF, PHI, vx, vy, vxf, vyf, ps, pf, Nx, Ny, dx, dy)
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
            HA[i,j]=(1-PHI[i,j])*tk1[i,j]*ALPHA[i,j]*dpsdt+ PHI[i,j]*tk1[i,j]*ALPHAF[i,j]*dpfdt
        end
    end
end


function f2_HA(HA, tk1, ALPHA, ALPHAF, PHI, vx, vy, vxf, vyf, ps, pf, Nx, Ny, dx, dy)
    for j=2:1:Nx, i=2:1:Ny
        # indirect calculation of dpdt
        # average vy, vx, vxf, vyf
        VXP = 0.5 * (vx[i, j]+vx[i, j-1])
        VYP = 0.5 * (vy[i, j]+vy[i-1, j])
        VXFP = 0.5 * (vxf[i, j]+vxf[i, j-1])
        VYFP = 0.5 * (vyf[i, j]+vyf[i-1, j])
        # evaluate DPsolid/Dt with upwind differences
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
        # dpsdx = ifelse(
        #     VXP<0, (ps[i, j]-ps[i, j-1])/dx, (ps[i, j+1]-ps[i, j])/dx)
        # dpsdy = ifelse(
        #     VYP<0, (ps[i, j]-ps[i-1, j])/dy, (ps[i+1, j]-ps[i, j])/dy)
        #     dpsdy=(ps[i,j]-ps[i-1,j])/dy
        # dpsdt = VXP*dpsdx + VYP*dpsdy
        # # evaluate DPfluid/Dt with upwind differences
        # dpfdx = ifelse(
        #     VXFP>0, (pf[i, j]-pf[i, j-1])/dx, (pf[i, j+1]-pf[i, j])/dx)
        # dpfdy = ifelse(
        #     VYFP>0, (pf[i, j]-pf[i-1, j])/dy, (pf[i+1, j]-pf[i, j])/dy)
        # dpfdt = VXFP*dpfdx + VYFP*dpfdy
        # compute adiabatic heating HA
        HA[i, j] = (
            (1-PHI[i, j]) * tk1[i, j] * ALPHA[i, j] * dpsdt
            + PHI[i, j] * tk1[i, j] * ALPHAF[i, j] * dpfdt
        )
    end
end












HA_ = rand(142,142)
tk1 = rand(142,142)
ALPHA = rand(142,142)
ALPHAF = rand(142,142)
PHI = rand(142,142)
vx = rand(142,142)
vy = rand(142,142)
vxf = rand(142,142)
vyf = rand(142,142)
ps = rand(142,142)
pf = rand(142,142)
VXP = zeros(142,142)
VYP = zeros(142,142)
VXFP = zeros(142,142)
VYFP = zeros(142,142)

∂ps∂xₗ = zeros(142,142)
∂ps∂xᵣ = zeros(142,142)
∂ps∂yₗ = zeros(142,142)
∂ps∂yᵣ = zeros(142,142)
∂pf∂xₗ = zeros(142,142)
∂pf∂xᵣ = zeros(142,142)
∂pf∂yₗ = zeros(142,142)
∂pf∂yᵣ = zeros(142,142)
dpsdx = zeros(142,142)
dpsdy = zeros(142,142)
dpfdx = zeros(142,142)
dpfdy = zeros(142,142)
dpsdt = zeros(142,142)
dpfdt = zeros(142,142)

sp = HydrologyPlanetesimals.StaticParameters()
Nx, Ny, Nx1, Ny1 = sp.Nx, sp.Ny, sp.Nx1, sp.Ny1
tk0 = rand(Ny1, Nx1)
tk1 = rand(Ny1, Nx1)
tk2 = rand(Ny1, Nx1)
DT = rand(Ny1, Nx1)
DT0 = rand(Ny1, Nx1)
RHOCP = rand(Ny1, Nx1)
KX = rand(Ny1, Nx1)
KY = rand(Ny1, Nx1)
HR = rand(Ny1, Nx1)
HA = rand(Ny1, Nx1)
HS = rand(Ny1, Nx1)
dtm = 1.0


function red(A)
    for k in 2:size(A, 3), j in 1:size(A, 2), i in 1:size(A, 1)
        @inbounds A[i, j, 1] += A[i, j, k]
    end
end



sp = HydrologyPlanetesimals.StaticParameters(
    Nxmc=4, Nymc=4, dsubgridt=0.5)
Nx, Ny = sp.Nx, sp.Ny
Nx1, Ny1 = sp.Nx1, sp.Ny1
x, y = sp.x, sp.y
xp, yp = sp.xp, sp.yp
dx, dy = sp.dx, sp.dy
marknum = sp.start_marknum
dtm = sp.dtelastic
dsubgridt = sp.dsubgridt
# simulate markers
xm = rand(-dx:0.1:x[end]+dx, marknum)
ym = rand(-dy:0.1:y[end]+dy, marknum)
tm = rand(1:3, marknum)
tkm = rand(marknum)
phim = rand(marknum)
tk1 = rand(Ny1, Nx1)
DT = rand(Ny1, Nx1)
TKSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
RHOCPSUM = zeros(Ny1, Nx1, Base.Threads.nthreads())
dtm = 0.9


sp = HydrologyPlanetesimals.StaticParameters(
    Nxmc=1, Nymc=1, dtelastic=0.9)
Nx, Ny = sp.Nx, sp.Ny
Nx1, Ny1 = sp.Nx1, sp.Ny1
marknum = sp.start_marknum
x, y = sp.x, sp.y
xp, yp = sp.xp, sp.yp
xvx, yvx = sp.xvx, sp.yvx
xvy, yvy = sp.xvy, sp.yvy
dx, dy = sp.dx, sp.dy
# simulate markers
x = 0:1000:140000
y = 0:1000:140000
dx,dy = 1000,1000
Nx, Ny = 140, 140
Nx1, Ny1 = Nx+1, Ny+1
Nxmc, Nymc =4,4
marknum = Nx*Ny*Nxmc*Nymc
xm = rand(-dx:0.1:x[end]+dx, marknum)
ym = rand(-dy:0.1:y[end]+dy, marknum)
sxym = rand(marknum)
sxxm = rand(marknum)
tkm = rand(marknum)
phim = rand(marknum)
tm = rand(1:3, marknum)
DT = rand(Ny1, Nx1)
tk2 = rand(Ny1, Nx1)
wyx = rand(Ny, Nx)
dtm = dtelastic
APHI = rand(Ny1, Nx1)
pr = rand(Ny1, Nx1)
pf = rand(Ny1, Nx1)
ps = rand(Ny1, Nx1)
pr0 = rand(Ny1, Nx1)
pf0 = rand(Ny1, Nx1)
ps0 = rand(Ny1, Nx1)

vx = rand(Ny1, Nx1)
vy = rand(Ny1, Nx1)
vxf = rand(Ny1, Nx1)
vyf = rand(Ny1, Nx1)
vxp = zeros(Ny1, Nx1)
vyp = zeros(Ny1, Nx1)
vxpf = zeros(Ny1, Nx1)
vypf = zeros(Ny1, Nx1)
    
    function fix1(x, y, x_axis, y_axis, dx, dy, jmin, jmax, imin, imax)
        # @timeit to "fix_weights" begin
            @inbounds begin
            j = min(max(trunc(Int, (x-x_axis[1])/dx)+1, jmin), jmax)
            i = min(max(trunc(Int, (y-y_axis[1])/dy)+1, imin), imax)
            dxmj = x - x_axis[j]
            dymi = y - y_axis[i]
            end # @inbounds
            weights = SVector(
                (1.0-dymi/dy) * (1.0-dxmj/dx),
                (dymi/dy) * (1.0-dxmj/dx),
                (1.0-dymi/dy) * (dxmj/dx),
                (dymi/dy) * (dxmj/dx)
            )
        # end # @timeit to "fix_weights"
            return i, j, weights
        end # function fix_weights

        function fix2(x, y, x_axis, y_axis, dx, dy, jmin, jmax, imin, imax)
            # @timeit to "fix_weights" begin
                @inbounds begin
                j = trunc(Int, (x-x_axis[1])/dx)+1
                if j < jmin
                    j = jmin
                elseif j > jmax
                    j = jmax
                end
                i = trunc(Int, (y-y_axis[1])/dy)+1
                if i < imin
                    i = imin
                elseif i > imax
                    i = imax
                end
                dxmj = x - x_axis[j]
                dymi = y - y_axis[i]
                end # @inbounds
                weights = SVector(
                    (1.0-dymi/dy) * (1.0-dxmj/dx),
                    (dymi/dy) * (1.0-dxmj/dx),
                    (1.0-dymi/dy) * (dxmj/dx),
                    (dymi/dy) * (dxmj/dx)
                )
            # end # @timeit to "fix_weights"
                return i, j, weights
            end # function fix_weights

function rk4(
    xm_ver,
    ym_ver,
    tm,
    tkm_ver,
    phim,
    sxym_ver,
    sxxm_ver,
    vx,
    vy,
    vxf,
    vyf,
    wyx,
    tk2,
    marknum,
    dtm,
    sp
)
@unpack Nx,
Ny,
dx,
dy,
x,
y,
xvx,
yvx,
xvy,
yvy,
xp,
yp,
jmin_basic,
jmax_basic,
imin_basic,
imax_basic,
jmin_vx,
jmax_vx,
imin_vx,
imax_vx,
jmin_vy,
jmax_vy,
imin_vy,
imax_vy,
jmin_p,
jmax_p,
imin_p,
imax_p,
brk4,
crk4,
rhocpsolidm,
rhocpfluidm = sp
vxm = zeros(4,1)
vym = zeros(4,1)
for m=1:1:marknum
    # Interpolate solid temperature for the initial marker location
    # Define i;j indexes for the upper left node
    j=trunc(Int, (xm_ver[m]-xp[1])/dx)+1
    i=trunc(Int, (ym_ver[m]-yp[1])/dy)+1
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
    # Compute distances
    dxmj=xm_ver[m]-xp[j]
    dymi=ym_ver[m]-yp[i]
    # Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy)
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy)
    wtmi1j1=(dxmj/dx)*(dymi/dy)
    # Compute Tsolid
    tksm0=tk2[i,j]*wtmij+tk2[i+1,j]*wtmi1j+ tk2[i,j+1]*wtmij1+tk2[i+1,j+1]*wtmi1j1;     
    # Interpolate local rotation rate
    # Define i;j indexes for the upper left node
    j=trunc(Int, (xm_ver[m]-x[1])/dx)+1
    i=trunc(Int, (ym_ver[m]-y[1])/dy)+1
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
    dxmj=xm_ver[m]-x[j]
    dymi=ym_ver[m]-y[i]
    # Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy)
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy)
    wtmi1j1=(dxmj/dx)*(dymi/dy)
    # Compute vx velocity
    omegam=wyx[i,j]*wtmij+wyx[i+1,j]*wtmi1j+ wyx[i,j+1]*wtmij1+wyx[i+1,j+1]*wtmi1j1
    # Analytical stress rotation using SIGMA"xx=-SIGMA"yy
    THETA=dtm*omegam; # Incremental rotation angle()
    sxxmnew=sxxm_ver[m]*cos(THETA)^2-sxxm_ver[m]*sin(THETA)^2-sxym_ver[m]*sin(2*THETA)
    sxymnew=sxxm_ver[m]*sin(2*THETA)+sxym_ver[m]*cos(2*THETA)
    sxxm_ver[m]=sxxmnew; sxym_ver[m]=sxymnew;    
    
    # Save initial marker coordinates
    xA=xm_ver[m]
    yA=ym_ver[m]
    for rk=1:1:4
        # Interpolate vx
        # Define i;j indexes for the upper left node
        j=trunc(Int, (xm_ver[m]-xvx[1])/dx)+1
        i=trunc(Int, (ym_ver[m]-yvx[1])/dy)+1
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
        # Compute distances
        dxmj=xm_ver[m]-xvx[j]
        dymi=ym_ver[m]-yvx[i]
        # Compute weights
        # Compute vx velocity for the top & bottom of the cell()
        vxm13=vx[i,j]*(1-dxmj/dx)+vx[i,j+1]*dxmj/dx
        vxm24=vx[i+1,j]*(1-dxmj/dx)+vx[i+1,j+1]*dxmj/dx
        # Compute correction
        if dxmj/dx>=0.5
            if j<Nx-1
                vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(vx[i,j]-2*vx[i,j+1]+vx[i,j+2])
                vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(vx[i+1,j]-2*vx[i+1,j+1]+vx[i+1,j+2])
            end
        else
            if j>1
                vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(vx[i,j-1]-2*vx[i,j]+vx[i,j+1])
                vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(vx[i+1,j-1]-2*vx[i+1,j]+vx[i+1,j+1])
            end
        end
        # Compute vx
        vxm[rk]=(1-dymi/dy)*vxm13+(dymi/dy)*vxm24
        
        # Interpolate vy
        # Define i;j indexes for the upper left node
        j=trunc(Int, (xm_ver[m]-xvy[1])/dx)+1
        i=trunc(Int, (ym_ver[m]-yvy[1])/dy)+1
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
        # Compute distances
        dxmj=xm_ver[m]-xvy[j]
        dymi=ym_ver[m]-yvy[i]
        # Compute weights
        # Compute vy velocity for the left & right of the cell()
        vym12=vy[i,j]*(1-dymi/dy)+vy[i+1,j]*dymi/dy
        vym34=vy[i,j+1]*(1-dymi/dy)+vy[i+1,j+1]*dymi/dy
        # Compute correction
        if dymi/dy>=0.5
            if i<Ny-1
                vym12=vym12+1/2*((dymi/dy-0.5)^2)*(vy[i,j]-2*vy[i+1,j]+vy[i+2,j])
                vym34=vym34+1/2*((dymi/dy-0.5)^2)*(vy[i,j+1]-2*vy[i+1,j+1]+vy[i+2,j+1])
            end      
        else
            if i>1
                vym12=vym12+1/2*((dymi/dy-0.5)^2)*(vy[i-1,j]-2*vy[i,j]+vy[i+1,j])
                vym34=vym34+1/2*((dymi/dy-0.5)^2)*(vy[i-1,j+1]-2*vy[i,j+1]+vy[i+1,j+1])
            end
        end
        # Compute vy
        vym[rk]=(1-dxmj/dx)*vym12+(dxmj/dx)*vym34
        
        # Change coordinates to obtain B;C;D points
        if rk==1 || rk==2
           xm_ver[m]=xA+dtm/2*vxm[rk]
           ym_ver[m]=yA+dtm/2*vym[rk]
        elseif rk==3
           xm_ver[m]=xA+dtm*vxm[rk]
           ym_ver[m]=yA+dtm*vym[rk]
        end
    end
    # Restore initial coordinates
    xm_ver[m]=xA
    ym_ver[m]=yA
    # Compute effective velocity
    vxmeff=1/6*(vxm[1]+2*vxm[2]+2*vxm[3]+vxm[4])
    vymeff=1/6*(vym[1]+2*vym[2]+2*vym[3]+vym[4])
    # Move markers
    xm_ver[m]=xm_ver[m]+dtm*vxmeff
    ym_ver[m]=ym_ver[m]+dtm*vymeff
    
    # Backtracing markers with fluid velocity
    xcur=xm_ver[m]
    ycur=ym_ver[m]
    xA=xcur
    yA=ycur
    for rk=1:1:4
        # Interpolate vx
        # Define i;j indexes for the upper left node
        j=trunc(Int, (xcur-xvx[1])/dx)+1
        i=trunc(Int, (ycur-yvx[1])/dy)+1
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
        # Compute distances
        dxmj=xcur-xvx[j]
        dymi=ycur-yvx[i]
        # Compute weights
        # Compute vx velocity for the top & bottom of the cell()
        vxm13=vxf[i,j]*(1-dxmj/dx)+vxf[i,j+1]*dxmj/dx
        vxm24=vxf[i+1,j]*(1-dxmj/dx)+vxf[i+1,j+1]*dxmj/dx
        # Compute correction
        if dxmj/dx>=0.5
            if j<Nx-1
                vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(vxf[i,j]-2*vxf[i,j+1]+vxf[i,j+2])
                vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(vxf[i+1,j]-2*vxf[i+1,j+1]+vxf[i+1,j+2])
            end
        else
            if j>1
                vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(vxf[i,j-1]-2*vxf[i,j]+vxf[i,j+1])
                vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(vxf[i+1,j-1]-2*vxf[i+1,j]+vxf[i+1,j+1])
            end
        end
        # Compute vx
        vxm[rk]=(1-dymi/dy)*vxm13+(dymi/dy)*vxm24
        
        # Interpolate vy
        # Define i;j indexes for the upper left node
        j=trunc(Int, (xcur-xvy[1])/dx)+1
        i=trunc(Int, (ycur-yvy[1])/dy)+1
        if j<1
            j=1
        elseif j>Nx
            j=Nx
        end
        if i<1
            i=1
        elseif  i>Ny-1
            i=Ny-1
        end
        # Compute distances
        dxmj=xcur-xvy[j]
        dymi=ycur-yvy[i]
        # Compute weights
        # Compute vy velocity for the left & right of the cell()
        vym12=vyf[i,j]*(1-dymi/dy)+vyf[i+1,j]*dymi/dy
        vym34=vyf[i,j+1]*(1-dymi/dy)+vyf[i+1,j+1]*dymi/dy
        # Compute correction
        if dymi/dy>=0.5
            if i<Ny-1
                vym12=vym12+1/2*((dymi/dy-0.5)^2)*(vyf[i,j]-2*vyf[i+1,j]+vyf[i+2,j])
                vym34=vym34+1/2*((dymi/dy-0.5)^2)*(vyf[i,j+1]-2*vyf[i+1,j+1]+vyf[i+2,j+1])
            end      
        else
            if i>1
                vym12=vym12+1/2*((dymi/dy-0.5)^2)*(vyf[i-1,j]-2*vyf[i,j]+vyf[i+1,j])
                vym34=vym34+1/2*((dymi/dy-0.5)^2)*(vyf[i-1,j+1]-2*vyf[i,j+1]+vyf[i+1,j+1])
            end
        end
        # Compute vy
        vym[rk]=(1-dxmj/dx)*vym12+(dxmj/dx)*vym34
        
        # Change coordinates to obtain B;C;D points
        if rk==1 || rk==2
            xcur=xA-dtm/2*vxm[rk]
            ycur=yA-dtm/2*vym[rk]
        elseif  rk==3
            xcur=xA-dtm*vxm[rk]
            ycur=yA-dtm*vym[rk]
        end
    end
    # Compute effective velocity
    vxmeff=1/6*(vxm[1]+2*vxm[2]+2*vxm[3]+vxm[4])
    vymeff=1/6*(vym[1]+2*vym[2]+2*vym[3]+vym[4])
    # Trace the node backward
    xcur=xA-dtm*vxmeff
    ycur=yA-dtm*vymeff
    # Interpolate fluid temperature
    # Define i;j indexes for the upper left node
    j=trunc(Int, (xcur-xp[1])/dx)+1
    i=trunc(Int, (ycur-yp[1])/dy)+1
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
    # Compute distances
    dxmj=xcur-xp[j]
    dymi=ycur-yp[i]
    # Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy)
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy)
    wtmi1j1=(dxmj/dx)*(dymi/dy)
    # Compute nodal Tfluid
    tkfm0=tk2[i,j]*wtmij+tk2[i+1,j]*wtmi1j+ tk2[i,j+1]*wtmij1+tk2[i+1,j+1]*wtmi1j1
    # Compute Tfluid-Tsolid for the marker
    dtkfsm=tkfm0-tksm0
    # Correct marker temperature
    tkm_ver[m]=((1-phim[m])*tkm_ver[m]*rhocpsolidm[tm[m]]+ phim[m]*(tkm_ver[m]+dtkfsm)*rhocpfluidm[tm[m]])/ ((1-phim[m])*rhocpsolidm[tm[m]]+phim[m]*rhocpfluidm[tm[m]])
end  # end of marker loop
end





sp = HydrologyPlanetesimals.StaticParameters()
Nx, Ny = sp.Nx, sp.Ny
Nx1, Ny1 = sp.Nx1, sp.Ny1
Nxm, Nym = sp.Nxm, sp.Nym
marknum = sp.start_marknum
x, y = sp.x, sp.y
xp, yp = sp.xp, sp.yp
xvx, yvx = sp.xvx, sp.yvx
xvy, yvy = sp.xvy, sp.yvy
xxm, yym = sp.xxm, sp.yym
dx, dy = sp.dx, sp.dy
jmin_m, jmax_m = sp.jmin_m, sp.jmax_m
imin_m, imax_m = sp.imin_m, sp.imax_m
# simulate markers
Nx, Ny = 140, 140
Nx1, Ny1 = Nx+1, Ny+1
dx, dy = 1000, 1000
x = 0:dx:140_000
y = 0:dy:140_000
Nxmc, Nymc = 4,4
dtm = 0.9
dsubgridt = 0.9
marknum = Nx*Ny*Nxmc*Nymc
xm = rand(-dx:0.1:x[end]+dx, marknum)
ym = rand(-dy:0.1:y[end]+dy, marknum)
xm2 = rand(-dx:0.1:x[end]/2, marknum)
ym2 = rand(-dy:0.1:y[end]/2, marknum)

tm = rand(1:3, marknum)
sxym = rand(marknum)
sxxm = rand(marknum)
tkm = rand(marknum)
phim = rand(marknum)
sxxm = rand(marknum)
sxym = rand(marknum)
etavpm = rand(marknum)
phim = rand(marknum)
tk1 = rand(marknum)
ktotalm = rand(marknum)
rhocptotalm = rand(marknum)
DT = zeros(Ny1, Nx1)
TKSUM = zeros(Ny1, Nx1, nthreads())
RHOCPSUM = zeros(Ny1, Nx1, nthreads())
mdis, mnum = HydrologyPlanetesimals.setup_marker_geometry_helpers()


function δy(i, j, grid)
    return grid[i, j] - 2.0*grid[i+1, j] + grid[i+2, j]
end

function loop(result)
    # for xpjj in SVector{140}(xp[2:1:Nx]), ypii in SVector{140}(yp[2:1:Ny])
    # for xpjj in xp[2:1:Nx], ypii in yp[2:1:Ny]
        # ypii = xpjj
    
        # for ypii in SVector{140}(yp[2:1:Ny])
    for jj=2:1:Nx, ii=2:1:Ny
    # setup RK4 scheme
    result = 0.0
    # xrk4 = [xp[jj], 0.0, 0.0, 0.0]
    # yrk4 = @MVector [yp[ii], 0.0, 0.0, 0.0]
    # RK4 velocities va, vb, vc, vd
    # vxrk4 = @MVector zeros(4)
    # vyrk4 = @MVector zeros(4)
    # backtrace P node using RK4 scheme on solid velocity
    for rk=1:1:4
#         # interpolate vx

        i, j, dxmj, dymi = HydrologyPlanetesimals.fix_distances(
            # xpjj,
            # ypii,
            xp[jj],
            yp[ii],
            # xrk4[rk],
            # yrk4[rk],
            xvx,
            yvx,
            dx,
            dy,
            jmin_vx,
            jmax_vx,
            imin_vx,
            imax_vx
        )
        rr[ii,jj] = dxmj+dymi
        # result += dxmj+dymi
    end
# end
end
return result
end


function loop2(t)
    result = 0.0
    for m=1:1:marknum
    # for jj=2:1:Nx
    #     for ii=2:1:Ny
        # # setup RK4 scheme
        # xrk4 = [xp[jj], 0.0, 0.0, 0.0]
        # yrk4 = @MVector [yp[ii], 0.0, 0.0, 0.0]
        # RK4 velocities va, vb, vc, vd
        # vxrk4 = @MVector zeros(4)
        # vyrk4 = @MVector zeros(4)
        # backtrace P node using RK4 scheme on solid velocity
        for rk=1:1:4
    #         # interpolate vx
    
    
            result += rhocpsolidm[t[m]] 
        end
    end
    return result
    end

    function loop3(t)
        result = 0.0
        for m=1:1:marknum::Int64
        # for m=1:1:marknum
        # for jj=2:1:Nx
        #     for ii=2:1:Ny
            # # setup RK4 scheme
            # xrk4 = [xp[jj], 0.0, 0.0, 0.0]
            # yrk4 = @MVector [yp[ii], 0.0, 0.0, 0.0]
            # RK4 velocities va, vb, vc, vd
            # vxrk4 = @MVector zeros(4)
            # vyrk4 = @MVector zeros(4)
            # backtrace P node using RK4 scheme on solid velocity
            for rk=1:1:4
        #         # interpolate vx
        
        
                result += t[m]
            end
        end
        return result
        end
    

        function add_vrk4(vrk4, v, rk)
            if rk == 1
                return vrk4 + @SVector [v, 0.0, 0.0, 0.0]
            elseif rk == 2
                return vrk4 + @SVector [0.0, v, 0.0, 0.0]
            elseif rk == 3
                return vrk4 + @SVector [0.0, 0.0, v, 0.0]
            elseif rk == 4
                return vrk4 + @SVector [0.0, 0.0, 0.0, v]
            else
                return vrk4
            end
        end

        function test_add_vrk4(vrk4, v)
            for rk=1:1:4
                vrk4 = add_vrk4(vrk4, v, rk)
            end
            return vrk4
        end


    function ∂vx∂x(dxmj, i, j, grid)
        ∂vx∂x₁₃ = grid[i, j-1] - 2.0*grid[i, j] + grid[i, j+1]
        ∂vx∂x₂₄ = grid[i+1, j-1] - 2.0*grid[i+1, j] + grid[i+1, j+1]
    end



using DifferentialEquations, Plots

function lorenz(du,u,p,t)
    du[1] = p[1]*(u[2]-u[1])
    du[2] = u[1]*(p[2]-u[3]) - u[2]
    du[3] = u[1]*u[2] - p[3]*u[3]
end

u0 = [1., 5., 10.]
tspan = (0., 100.)
p = (10.0,28.0,8/3)
prob = ODEProblem(lorenz, u0, tspan,p)
sol = solve(prob);

t = range(sol.prob.tspan...; length=10^4+1)
X, Y, Z = ((t -> sol(t)[i]).(t) for i in 1:3)
xlim, ylim, zlim = extrema.((X, Y, Z));

gr(fmt = :png)
anim = @animate for i in 1:100:length(X)
    @views x, y, z = X[1:i], Y[1:i], Z[1:i]
    A = plot(x, y, z; label="", lw=0.5, xlim, ylim, zlim)
    B = plot(x, y; label="", lw=0.2, title="x, y", titlefontsize=8, xlim=xlim, ylim=ylim)
    C = plot(x, z; label="", lw=0.2, title="x, y", titlefontsize=8, xlim=xlim, ylim=zlim)
    # D = plot(y, z; label="", lw=0.2, title="x, y", titlefontsize=8, xlim=ylim, ylim=zlim)
    D = heatmap(rand(100,100))
    layout = @layout [
        a{0.7h}
        [b c d]
    ]
    plot(A, B, C, D; layout, size=(640, 640))
end
gif(anim, "lorenz2.gif")


using LinearAlgebra
using Plots
gr(fmt = :png)

function arrow_!(x, y, u, v; lar=0.1, lc=:black, linewidth=1.0, la=1.0)
    v1, v2 = normalize([u; v]), normalize([-v; u])
    v3 = normalize(3*v1 + v2) * norm((u, v)) * lar
    v4 = normalize(3*v1 - v2) * norm((u, v)) * lar
    plot!([x, x+u], [y, y+v], lc=lc, la=la)
    plot!([x+u, x+u-v4[1]], [y+v, y+v-v4[2]], lc=lc, linewidth=linewidth, la=la)
    plot!([x+u, x+u-v3[1]], [y+v, y+v-v3[2]], lc=lc, linewidth=linewidth, la=la)
end

N = 30;
xx, yy = 1 .+ 2* rand(N), 1 .+ 2*rand(N);  # points
r, θ = rand(N), LinRange(0,2π,N)
u, v = r .* cos.(θ), r .* sin.(θ)    # arrows 

# plot points and arrows with 10% head sizes
scatter(xx, yy, mc=:red, ms=2.5, ratio=1, ma=0.5, legend=false)
for (x,y,u,v) in zip(xx,yy,u,v)
    display(arrow_!(x, y, u, v; as=0.1, lc=:blue, la=1))
end


meshgrid(x, y) = (repeat(x, outer=length(y)), repeat(y, inner=length(x)))

function custom_quiver(
    x, y; quiver=(u, v), mc=:black, ms=2.5, ma=0.5, lar=0.1, lc=:black, linewidth=1.0, la=1.0)
    scatter(x, y, mc=mc, ms=ms, ma=ma, legend=false)
    for (x, y, u, v) in zip(x, y, vec(u'), vec(v'))
        display(arrow_!(x, y, u, v; lar=lar, lc=lc, linewidth=linewidth, la=la))
    end
end

function f(x, y)
    return x, y
end


using AbstractPlotting, Makie
using ImageFiltering, LinearAlgebra

x = range(-2, stop = 2, length = 21)
y = x
z = x .* exp.(-x .^ 2 .- (y') .^ 2)
scene = Makie.contour(x, y, z, levels = 10, linewidth = 3)
u, v = ImageFiltering.imgradients(z, KernelFactors.ando3)
n = vec(norm.(Vec2f0.(u,v)))
Makie.arrows!(x, y, u, v, arrowsize = n, arrowcolor = n)


using GLMakie

function mandelbrot(x, y)
    z = c = x + y*im
    for i in 1:30.0; abs(z) > 2 && return i; z = z^2 + c; end; 0
end

x = LinRange(-2, 1, 200)
y = LinRange(-1.1, 1.1, 200)
matrix = mandelbrot.(x, y')
fig, ax, hm = heatmap(x, y, matrix)

N = 50
xmin = LinRange(-2.0, -0.72, N)
xmax = LinRange(1, -0.6, N)
ymin = LinRange(-1.1, -0.51, N)
ymax = LinRange(1, -0.42, N)

# we use `record` to show the resulting video in the docs.
# If one doesn't need to record a video, a normal loop works as well.
# Just don't forget to call `display(fig)` before the loop
# and without record, one needs to insert a yield to yield to the render task
record(fig, "heatmap_mandelbrot.mp4", 1:7:N) do i
    _x = LinRange(xmin[i], xmax[i], 200)
    _y = LinRange(ymin[i], ymax[i], 200)
    hm[1] = _x # update x coordinates
    hm[2] = _y # update y coordinates
    hm[3] = mandelbrot.(_x, _y') # update data
    autolimits!(ax) # update limits
    # yield() -> not required with record
end
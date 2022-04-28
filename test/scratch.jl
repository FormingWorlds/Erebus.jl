# using Base.Threads
# using SparseArrays
# using MAT
# using DocStringExtensions
# using Parameters
# using StaticArrays
# using BenchmarkTools
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


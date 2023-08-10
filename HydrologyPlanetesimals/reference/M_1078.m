function [xm_jl,ym_jl,tk0_jl,tk1_jl,tk2_jl,DT_jl,APHI_jl,vxf_jl,vyf_jl,pscale_jl,ETA_jl,ETAP_jl,...
    GGG_jl,GGGP_jl,SXY0_jl,SXX0_jl,RHOX_jl,RHOY_jl,RHOFX_jl,RHOFY_jl,...
    RX_jl,RY_jl,ETAPHI_jl,BETTAPHI_jl,PHI_jl,gx_jl,gy_jl,pr0_jl,pf0_jl,...
    dt_jl,R_jl,L_jl,S_jl,vx_jl,vy_jl,qxD_jl,qyD_jl,pr_jl,pf_jl] = M_1078(step)
% fpath = 'C:\Users\ich\outTest\M_1078_';
fpath = '/Users/z7717/test/M_1078_';
xm_jl = h5read([fpath num2str(step) '.jld2'], '/xm')';
ym_jl = h5read([fpath num2str(step) '.jld2'], '/ym')';
tk0_jl = h5read([fpath num2str(step) '.jld2'], '/tk0');
tk1_jl = h5read([fpath num2str(step) '.jld2'], '/tk1');
tk2_jl = h5read([fpath num2str(step) '.jld2'], '/tk2');
DT_jl = h5read([fpath num2str(step) '.jld2'], '/DT');
APHI_jl = h5read([fpath num2str(step) '.jld2'], '/APHI');
vxf_jl = h5read([fpath num2str(step) '.jld2'], '/vxf');
vyf_jl = h5read([fpath num2str(step) '.jld2'], '/vyf');
pscale_jl =  h5read([fpath num2str(step) '.jld2'], '/Kcont');
ETA_jl = h5read([fpath num2str(step) '.jld2'], '/ETA');
ETAP_jl = h5read([fpath num2str(step) '.jld2'], '/ETAP');
GGG_jl = h5read([fpath num2str(step) '.jld2'], '/GGG');
GGGP_jl = h5read([fpath num2str(step) '.jld2'], '/GGGP');
SXY0_jl = h5read([fpath num2str(step) '.jld2'], '/SXY0');
SXX0_jl = h5read([fpath num2str(step) '.jld2'], '/SXX0');
RHOX_jl = h5read([fpath num2str(step) '.jld2'], '/RHOX');
RHOY_jl = h5read([fpath num2str(step) '.jld2'], '/RHOY');
RHOFX_jl = h5read([fpath num2str(step) '.jld2'], '/RHOFX');
RHOFY_jl = h5read([fpath num2str(step) '.jld2'], '/RHOFY');
RX_jl = h5read([fpath num2str(step) '.jld2'], '/RX');
RY_jl = h5read([fpath num2str(step) '.jld2'], '/RY');
ETAPHI_jl = h5read([fpath num2str(step) '.jld2'], '/ETAPHI');
BETTAPHI_jl = h5read([fpath num2str(step) '.jld2'], '/BETTAPHI');
PHI_jl = h5read([fpath num2str(step) '.jld2'], '/PHI');
gx_jl = h5read([fpath num2str(step) '.jld2'], '/gx');
gy_jl = h5read([fpath num2str(step) '.jld2'], '/gy');
pr0_jl = h5read([fpath num2str(step) '.jld2'], '/pr0');
pf0_jl = h5read([fpath num2str(step) '.jld2'], '/pf0');
dt_jl = h5read([fpath num2str(step) '.jld2'], '/dt');
R_jl = h5read([fpath num2str(step) '.jld2'], '/R');
L_jl = [];
S_jl = h5read([fpath num2str(step) '.jld2'], '/S');
vx_jl = h5read([fpath num2str(step) '.jld2'], '/vx');
vy_jl = h5read([fpath num2str(step) '.jld2'], '/vy');
qxD_jl = h5read([fpath num2str(step) '.jld2'], '/qxD');
qyD_jl = h5read([fpath num2str(step) '.jld2'], '/qyD');
pr_jl = h5read([fpath num2str(step) '.jld2'], '/pr');
pf_jl = h5read([fpath num2str(step) '.jld2'], '/pf');
end
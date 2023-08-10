% Visco-elasto-plastic hydro-thermomechanical (HTM) planetary code
% Solving Poisson, momentum, mass and energy conservation eqs.
% for self-gravitating coupled fluid-solid system
% in primitive variable formulation
% with variable viscosity and thermal conductivity
% using FD with staggered grid
% Clearing memory and figures

% % Load mat file
% fdata=fopen('file.txt','rt');
% timestep=fscanf(fdata,'%d',1);
% fclose(fdata);
% if(timestep>0)
% namemat    =  ['planetary-NGD_',num2str(timestep)];
% load(namemat);
% else % if uncommented, uncoment end in line 275
clear all;

%Switch for radioactive heating
hr_al=1;  %if 1 radioactive heating from 26Al active
hr_fe=0;  %if 1 radioactive heating from 60Fe active

% Define Numerical model
xsize=140000; % Horizontal model size, m
ysize=140000; % Vertical model size, m
Nx=141; % Horizontal grid resolution
Ny=141; % Vertical grid resolution
Nx1=Nx+1;
Ny1=Ny+1;
dx=xsize/(Nx-1); % Horizontal grid step, m
dy=ysize/(Ny-1); % Vertical grid step, m

% Define Gravity
G=6.672e-11; % Gravity constant, N*m^2/kg^2

% Coordinates of different nodal points
% Basic nodes
x=0:dx:xsize; % Horizontal coordinates of basic grid points, m
y=0:dy:ysize; % Vertical coordinates of basic grid points, m
% Vx-Nodes
xvx=0:dx:xsize+dy; % Horizontal coordinates of vx grid points, m
yvx=-dy/2:dy:ysize+dy/2; % Vertical coordinates of vx grid points, m
% Vy-nodes
xvy=-dx/2:dx:xsize+dx/2; % Horizontal coordinates of vy grid points, m
yvy=0:dy:ysize+dy; % Vertical coordinates of vy grid points, m
% P-Nodes
xp=-dx/2:dx:xsize+dx/2; % Horizontal coordinates of P grid points, m
yp=-dy/2:dy:ysize+dy/2; % Vertical coordinates of P grid points, m

% Nodal arrays
% Basic nodes
ETA=zeros(Ny,Nx); % Viscoplastic Viscosity, Pa*s
ETA0=zeros(Ny,Nx); % Viscous Viscosity, Pa*s
GGG=zeros(Ny,Nx); % Shear modulus, Pa
EXY=zeros(Ny,Nx); % EPSILONxy, 1/s
SXY=zeros(Ny,Nx); % SIGMAxy, 1/s
SXY0=zeros(Ny,Nx); % SIGMA0xy, 1/s
wyx=zeros(Ny,Nx); % Rotation rate, 1/s
COH=zeros(Ny,Nx); % Compressive strength, Pa
TEN=zeros(Ny,Nx); % Tensile strength, Pa
FRI=zeros(Ny,Nx); % Friction
YNY=zeros(Ny,Nx); % Plastic yielding mark, 1=yes,0=no
% Vx-Nodes
RHOX=zeros(Ny1,Nx1); % Density, kg/m^3
RHOFX=zeros(Ny1,Nx1); % Fluid Density, kg/m^3
KX=zeros(Ny1,Nx1); % Thermal conductivity, W/m/K
PHIX=zeros(Ny1,Nx1); % Porosity
vx=zeros(Ny1,Nx1); % Solid vx-velocity m/s
vxf=zeros(Ny1,Nx1); % Fluid vx-velocity m/s
RX=zeros(Ny1,Nx1); % ETAfluid/Kphi ratio , m^2
qxD=zeros(Ny1,Nx1); % qx-Darcy flux m/s
gx=zeros(Ny1,Nx1); % gx-gravity, m/s^2
% Vy-Nodes
RHOY=zeros(Ny1,Nx1); % Density, kg/m^3
RHOFY=zeros(Ny1,Nx1); % Fluid Density, kg/m^3
KY=zeros(Ny1,Nx1); % Thermal conductivity, W/m/K
PHIY=zeros(Ny1,Nx1); % Porosity
vy=zeros(Ny1,Nx1); % Solid vy-velocity m/s
vyf=zeros(Ny1,Nx1); % Fluid vy-velocity m/s
RY=zeros(Ny1,Nx1); % ETAfluid/Kphi ratio , m^2
qyD=zeros(Ny1,Nx1); % qy-Darcy flux m/s
gy=zeros(Ny1,Nx1); % gy-gravity, m/s^2
% P-nodes
RHO=zeros(Ny1,Nx1); % Density, kg/m^3
RHOCP=zeros(Ny1,Nx1); % Volumetric heat capacity, J/m^3/K
ALPHA=zeros(Ny1,Nx1); % Thermal expansion, J/m^3/K
ALPHAF=zeros(Ny1,Nx1); % Fluid Thermal expansion, J/m^3/K
HR=zeros(Ny1,Nx1); % Radioactive heating, W/m^3
HA=zeros(Ny1,Nx1); % Adiabatic heating, W/m^3
HS=zeros(Ny1,Nx1); % Shear heating, W/m^3
ETAP=zeros(Ny1,Nx1); % Viscosity, Pa*s
GGGP=zeros(Ny1,Nx1); % Shear modulus, Pa
EXX=zeros(Ny,Nx); % EPSILONxx, 1/s
SXX=zeros(Ny,Nx); % SIGMA'xx, 1/s
SXX0=zeros(Ny,Nx); % SIGMA0'xx, 1/s
tk1=zeros(Ny1,Nx1); % Old temperature, K
tk2=zeros(Ny1,Nx1); % New temperature, K
vxp=zeros(Ny1,Nx1); % Solid Vx in pressure nodes, m/s
vyp=zeros(Ny1,Nx1); % Solid Vy in pressure nodes, m/s
vxpf=zeros(Ny1,Nx1); % Fluid Vx in pressure nodes, m/s
vypf=zeros(Ny1,Nx1); % Fluid Vy in pressure nodes, m/s
pr=zeros(Ny1,Nx1); % Total Pressure, Pa
pf=zeros(Ny1,Nx1); % Fluid Pressure, Pa
ps=zeros(Ny1,Nx1); % Solid Pressure, Pa
pr0=zeros(Ny1,Nx1); % Old Total Pressure, Pa
pf0=zeros(Ny1,Nx1); % Old Fluid Pressure, Pa
ps0=zeros(Ny1,Nx1); % Old Solid Pressure, Pa
ETAPHI=zeros(Ny1,Nx1); % Bulk Viscosity, Pa*s
BETTAPHI=zeros(Ny1,Nx1); % Bulk compresibility, Pa*s
PHI=zeros(Ny1,Nx1); % porosity
APHI=zeros(Ny1,Nx1); % Dln((1-PHI)/PHI)/Dt
FI=zeros(Ny1,Nx1); % Gravity potential, J/kg

% Define markers
Nxmc=4; % Number of markers per cell in horizontal direction
Nymc=4; % Number of markers per cell in vertical direction
Nxm=(Nx-1)*Nxmc; % Marker grid resolution in horizontal direction
Nym=(Ny-1)*Nymc; % Marker grid resolution in vertical direction
dxm=xsize/Nxm; % Marker grid step in horizontal direction,m
dym=ysize/Nym; % Marker grid step in vertical direction,m
marknum=Nxm*Nym; % Number of markers
xm=zeros(1,marknum); % Horizontal coordinates, m
ym=zeros(1,marknum); % Vertical coordinates, m
tm=zeros(1,marknum); % Material type
tkm=zeros(1,marknum); % Marker temperature, K
sxxm=zeros(1,marknum); % SIGMA'xx, Pa
sxym=zeros(1,marknum); % SIGMAxy, Pa
etavpm=zeros(1,marknum); % Visco-plastic viscosity, Pa
phim=zeros(1,marknum); % Marker porosity

% Define properties of materials: 
%            Planet  Crust Space
rhosolidm   = [3300   3300   1     ]; % Solid Density, kg/m^3
rhofluidm   = [1000   1000   1     ]; % Fluid Density, kg/m^3
etasolidm   = [1e+19  1e+19  1e+16 ]; % Solid Viscosity, Pa s
etasolidmm  = [1e+19  1e+19  1e+16 ]; % Molten Solid Viscosity, Pa s
% CHANGE!!!
etafluidm   = [1e+12  1e+12  1e-3 ]; % Fluid Viscosity, Pa s
etafluidmm  = [1e-3   1e-3   1e-3  ]; % Molten Fluid Viscosity, Pa s
rhocpsolidm = [3.3e+6 3.3e+6 3.0e+6]; % Solid Volumetric heat capacity, kg/m^3
rhocpfluidm = [1.0e+6 1.0e+6 3.0e+6]; % Fluid Volumetric heat capacity, kg/m^3
alphasolidm = [3e-5   3e-5   0     ]; % Solid Thermal expansion, 1/K
alphafluidm = [5e-5   5e-5   0     ]; % Fluid Thermal expansion, 1/K
ksolidm     = [3      3      3000  ]; % Solid Thermal conductivity, W/m/K
kfluidm     = [50     50     3000  ]; % Fluid Thermal conductivity, W/m/K
hrsolidm    = [0      0      0     ]; % Solid Radiogenic heat production, W/m^3
hrfluidm    = [0      0      0     ]; % Fluid Radiogenic heat production, W/m^3
gggsolidm   = [1e+10  1e+10  1e+10 ]; % Solid Shear Modulus, Pa
frictsolidm = [0.6    0.6      0   ]; % Solid Friction coefficient
cohessolidm = [1e+8   1e+8   1e+8  ]; % Solid compressive strength, Pa
tenssolidm  = [6e+7   6e+7   6e+7  ]; % Solid tensile strength, Pa
kphim0      = [1e-13  1e-13  1e-17 ]; % Standard permeability, m^2
etaphikoef=1; % Koefficient to compute compaction viscosity from shear viscosity

% Constants for 26Al decay
t_half_al=717000*31540000; %s
tau_al=t_half_al/log(2);
ratio_al=5.0e-5; % ratio of 26Al and 27Al Isotopes 
E_al=5.0470e-13; % [J]
f_al=1.9e23; %atoms/kg

%Constants for 60Fe decay
t_half_fe=2620000*31540000; %s
tau_fe=t_half_fe/log(2);
ratio_fe=1e-6; %Initial 60Fe/56Fe
E_fe=4.34e-13; %[J]
f_fe=1.957e24; %atoms/kg

% Melting temperatures
tmsilicate=1e+6;%1416; 
% CHANGE!!!
tmiron=273; 

phim0=0.2; % standard iron fraction (porosity)
phimin=1e-4; % Min porosity
phimax=1-phimin; % Max porosity

% Define marker coordinates, temperature and material type
rplanet=50000; % Planetary radius
rcrust=48000; % Crust radius
psurface=1e+3; % Surface pressure
m=1; % Marker counter
for jm=1:1:Nxm
    for im=1:1:Nym
        % Define marker coordinates
        xm(m)=dxm/2+(jm-1)*dxm;%+(rand-0.5)*dxm;
        ym(m)=dym/2+(im-1)*dym;%+(rand-0.5)*dym;
        % Marker properties
        rmark=((xm(m)-xsize/2)^2+(ym(m)-ysize/2)^2)^0.5;
        if(rmark<rplanet)
            % Planet
            tm(m)=1; % mantle
            if(rmark>rcrust) 
                tm(m)=2; % crust
            end
% CHANGE!!!
            tkm(m)=170; % Temperature
            phim(m)=phim0;%*(1+1.0*(rand-0.5)); % Porosity
            etavpm(m)=etasolidm(tm(m));%*exp(-28*phim(m)); % Matrix viscosity
        else
            % Sticky space (to have internal free surface)
            tm(m)=3; % Material type
% CHANGE!!!
            tkm(m)=170; % Temperature
            phim(m)=phimin; % Porosity
            etavpm(m)=etasolidm(tm(m)); % Matrix viscosity
        end
        % Update marker counter
        m=m+1;
    end
end

% Introducing scaled pressure
pscale=1e+23/dx;

% Define global matrixes 
% Hydro-Mechanical solution: L(), R()
N=Nx1*Ny1*6; % Global number of unknowns
L=sparse(N,N); % Matrix of coefficients (left part)
R=zeros(N,1); % Vector of right parts
% Thermal solution: LT(), RT()
N=Nx1*Ny1; % Global number of unknowns
LT=sparse(N,N); % Matrix of coefficients (left part)
RT=zeros(N,1); % Vector of right parts
% Gravity solution: LP(), RP()
N=Nx1*Ny1; % Global number of unknowns
LP=sparse(N,N); % Matrix of coefficients (left part)
RP=zeros(N,1); % Vector of right parts

% Mechanical boundary conditions: free slip=-1; No Slip=1
bcleft=-1;
bcright=-1;
bctop=-1;
bcbottom=-1;
% Hydraulic boundary conditions: free slip=-1; No Slip=1
bcfleft=-1;
bcfright=-1;
bcftop=-1;
bcfbottom=-1;
% Extension/shortening velocities
strainrate=0e-13; % Shortening strain rate
vxleft=strainrate*xsize/2;
vxright=-strainrate*xsize/2;
vytop=-strainrate*ysize/2;
vybottom=strainrate*ysize/2;

% Thermal boundary conditions: insulation at all boundaries

% Timestepping
nname='planetary-NGD_'; %mat filename
nnamefile='planetary-NGD_'; % jpeg filename Fig1

savematstep=10; %.mat storage periodicity
dtelastic=1e+11; % Maximal computational timestep, s
dt=dtelastic; % Current computational timestep, s
dtkoef=2; % Koefficient to decrese computational timestep
dtkoefup=1.1; % Koefficient to increase computational timestep
dtstep=200; % Number of iterations before changing computational timestep
dxymax=0.05; % Max marker movement per time step, grid steps
vpratio=1/3; % Weight of averaged velocity for moving markers
DTmax=20; % Max temperature change per time step, K
dsubgridt=0; % Subgrid temperature diffusion parameter
dsubgrids=0; % Subgrid stress diffusion parameter
timesum=2.25e6*365.25*24*3600; % Time sum, s
etamin=1e+12; % Lower viscosity cut-off, Pa s
etamax=1e+23; % Upper viscosity cut-off, Pa s
nplast=100000; % Number of plastic iterations
visstep=1; % Periodicity of visualization
yerrmax=1e+2; % Tolerance level for yielding error
YERRNOD=zeros(1,nplast); % Yielding error of nodes
etawt=0; % Weight for old viscosity
dphimax=0.01; % max porosity ratio change per time step
nsteps=20; % number of timesteps
timestep=1;
% end % comment out if restarts disabled at beginning of code
savematstep=1; %.mat storage periodicity
% convergence evaluation
maxT_m = zeros(nsteps,1);
maxT_j = zeros(nsteps,1);
meanT_10_m = zeros(nsteps,1);
meanT_10_j = zeros(nsteps,1);
meanT_40_m = zeros(nsteps,1);
meanT_40_j = zeros(nsteps,1);
dRHO_m = zeros(nsteps,1);
dRHO_j = zeros(nsteps,1);
meanRHO_10_m = zeros(nsteps,1);
meanRHO_10_j = zeros(nsteps,1);
meanRHO_40_m = zeros(nsteps,1);
meanRHO_40_j = zeros(nsteps,1);

for timestep=timestep:1:nsteps

% Updating radioactive heating
%26Al
if hr_al==1
    Q_al=f_al*ratio_al*E_al*exp(-timesum/tau_al)/tau_al; % [W/kg]
    hrsolidm=Q_al*rhosolidm;    %radiogenic heat production [W/m^3]
    hrsolidm(1,3)=0; %no radioactive heating from space
end

% CHANGE!!!
% %60Fe 
% if hr_fe==1
%     Q_fe=f_fe*ratio_fe*E_fe*exp(-timesum/tau_fe)/tau_fe; %[W/kg]
%     hrfluidm(1,1)=Q_fe*rhofluidm(1,1); %[w/m^3]radioactive heatproduction only in planet
% end
% No iron in water
hrfluidm(1,1)=0;

% Save old stresses
sxxm00=sxxm; 
sxym00=sxym;    
    
% Interpolate properties from markers to nodes
% Basic nodes
ETA0SUM=zeros(Ny,Nx);
ETASUM=zeros(Ny,Nx);
GGGSUM=zeros(Ny,Nx);
SXYSUM=zeros(Ny,Nx);
COHSUM=zeros(Ny,Nx);
TENSUM=zeros(Ny,Nx);
FRISUM=zeros(Ny,Nx);
WTSUM=zeros(Ny,Nx);
% Vx-nodes
RHOXSUM=zeros(Ny1,Nx1);
RHOFXSUM=zeros(Ny1,Nx1);
KXSUM=zeros(Ny1,Nx1);
PHIXSUM=zeros(Ny1,Nx1);
RXSUM=zeros(Ny1,Nx1);
WTXSUM=zeros(Ny1,Nx1);
% Vy-nodes
RHOYSUM=zeros(Ny1,Nx1);
RHOFYSUM=zeros(Ny1,Nx1);
KYSUM=zeros(Ny1,Nx1);
PHIYSUM=zeros(Ny1,Nx1);
RYSUM=zeros(Ny1,Nx1);
WTYSUM=zeros(Ny1,Nx1);
% P-Nodes
GGGPSUM=zeros(Ny1,Nx1);
SXXSUM=zeros(Ny1,Nx1);
RHOSUM=zeros(Ny1,Nx1);
RHOCPSUM=zeros(Ny1,Nx1);
ALPHASUM=zeros(Ny1,Nx1);
ALPHAFSUM=zeros(Ny1,Nx1);
HRSUM=zeros(Ny1,Nx1);
TKSUM=zeros(Ny1,Nx1);
PHISUM=zeros(Ny1,Nx1);
WTPSUM=zeros(Ny1,Nx1);

for m=1:1:marknum


        % Compute marker parameters
    if(tm(m)<3)
        % Rocks
        kphim=kphim0(tm(m))*(phim(m)/phim0)^3/((1-phim(m))/(1-phim0))^2; %Permeability
        rhototalm=rhosolidm(tm(m))*(1-phim(m))+rhofluidm(tm(m))*phim(m);
        rhocptotalm=rhocpsolidm(tm(m))*(1-phim(m))+rhocpfluidm(tm(m))*phim(m);
        etasolidcur=etasolidm(tm(m));
        if(tkm(m)>tmsilicate)
            etasolidcur=etasolidmm(tm(m));
        end
        hrtotalm=hrsolidm(tm(m))*(1-phim(m))+hrfluidm(tm(m))*phim(m);
        ktotalm=(ksolidm(tm(m))*kfluidm(tm(m))/2+((ksolidm(tm(m))*(3*phim(m)-2)+...
            kfluidm(tm(m))*(1-3*phim(m)))^2)/16)^0.5-(ksolidm(tm(m))*(3*phim(m)-2)+...
            kfluidm(tm(m))*(1-3*phim(m)))/4;
        gggtotalm=gggsolidm(tm(m));
        fricttotalm=frictsolidm(tm(m));
        cohestotalm=cohessolidm(tm(m));
        tenstotalm=tenssolidm(tm(m));
        etafluidcur=etafluidm(tm(m));
        rhofluidcur=rhofluidm(tm(m));
        if(tkm(m)>tmiron)
            etafluidcur=etafluidmm(tm(m));
        end
        etatotalm=max(etamin,max(etafluidcur,etasolidcur));%*exp(-28*phim(m))));
    else
        % Sticky air
        kphim=kphim0(tm(m))*(phim(m)/phim0)^3/((1-phim(m))/(1-phim0))^2; %Permeability
        rhototalm=rhosolidm(tm(m));
        rhocptotalm=rhocpsolidm(tm(m));
        etatotalm=etasolidm(tm(m));
        hrtotalm=hrsolidm(tm(m));
        ktotalm=ksolidm(tm(m));
        gggtotalm=gggsolidm(tm(m));
        fricttotalm=frictsolidm(tm(m));
        cohestotalm=cohessolidm(tm(m));
        tenstotalm=tenssolidm(tm(m));
        rhofluidcur=rhofluidm(tm(m));
        etafluidcur=etafluidm(tm(m));
    end

  
    % Interpolation to basic nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-x(1))/dx)+1;
    i=fix((ym(m)-y(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-x(j);
    dymi=ym(m)-y(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    % i,j Node
    ETA0SUM(i,j)=ETA0SUM(i,j)+etatotalm*wtmij;
    ETASUM(i,j)=ETASUM(i,j)+etavpm(m)*wtmij;
    GGGSUM(i,j)=GGGSUM(i,j)+1/gggtotalm*wtmij;
    SXYSUM(i,j)=SXYSUM(i,j)+sxym(m)*wtmij;
    COHSUM(i,j)=COHSUM(i,j)+cohestotalm*wtmij;
    TENSUM(i,j)=TENSUM(i,j)+tenstotalm*wtmij;
    FRISUM(i,j)=FRISUM(i,j)+fricttotalm*wtmij;
    WTSUM(i,j)=WTSUM(i,j)+wtmij;
    % i+1,j Node
    ETA0SUM(i+1,j)=ETA0SUM(i+1,j)+etatotalm*wtmi1j;
    ETASUM(i+1,j)=ETASUM(i+1,j)+etavpm(m)*wtmi1j;
    GGGSUM(i+1,j)=GGGSUM(i+1,j)+1/gggtotalm*wtmi1j;
    SXYSUM(i+1,j)=SXYSUM(i+1,j)+sxym(m)*wtmi1j;
    COHSUM(i+1,j)=COHSUM(i+1,j)+cohestotalm*wtmi1j;
    TENSUM(i+1,j)=TENSUM(i+1,j)+tenstotalm*wtmi1j;
    FRISUM(i+1,j)=FRISUM(i+1,j)+fricttotalm*wtmi1j;
    WTSUM(i+1,j)=WTSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    ETA0SUM(i,j+1)=ETA0SUM(i,j+1)+etatotalm*wtmij1;
    ETASUM(i,j+1)=ETASUM(i,j+1)+etavpm(m)*wtmij1;
    GGGSUM(i,j+1)=GGGSUM(i,j+1)+1/gggtotalm*wtmij1;
    SXYSUM(i,j+1)=SXYSUM(i,j+1)+sxym(m)*wtmij1;
    COHSUM(i,j+1)=COHSUM(i,j+1)+cohestotalm*wtmij1;
    TENSUM(i,j+1)=TENSUM(i,j+1)+tenstotalm*wtmij1;
    FRISUM(i,j+1)=FRISUM(i,j+1)+fricttotalm*wtmij1;
    WTSUM(i,j+1)=WTSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    ETA0SUM(i+1,j+1)=ETA0SUM(i+1,j+1)+etatotalm*wtmi1j1;
    ETASUM(i+1,j+1)=ETASUM(i+1,j+1)+etavpm(m)*wtmi1j1;
    GGGSUM(i+1,j+1)=GGGSUM(i+1,j+1)+1/gggtotalm*wtmi1j1;
    SXYSUM(i+1,j+1)=SXYSUM(i+1,j+1)+sxym(m)*wtmi1j1;
    COHSUM(i+1,j+1)=COHSUM(i+1,j+1)+cohestotalm*wtmi1j1;
    TENSUM(i+1,j+1)=TENSUM(i+1,j+1)+tenstotalm*wtmi1j1;
    FRISUM(i+1,j+1)=FRISUM(i+1,j+1)+fricttotalm*wtmi1j1;
    WTSUM(i+1,j+1)=WTSUM(i+1,j+1)+wtmi1j1;    
    
    
    % Interpolation to vx-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xvx(1))/dx)+1;
    i=fix((ym(m)-yvx(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xvx(j);
    dymi=ym(m)-yvx(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    % i,j Node
    RHOXSUM(i,j)=RHOXSUM(i,j)+rhototalm*wtmij;
    RHOFXSUM(i,j)=RHOFXSUM(i,j)+rhofluidcur*wtmij;
    KXSUM(i,j)=KXSUM(i,j)+ktotalm*wtmij;
    PHIXSUM(i,j)=PHIXSUM(i,j)+phim(m)*wtmij;
    RXSUM(i,j)=RXSUM(i,j)+etafluidcur/kphim*wtmij;
    WTXSUM(i,j)=WTXSUM(i,j)+wtmij;
    % i+1,j Node
    RHOXSUM(i+1,j)=RHOXSUM(i+1,j)+rhototalm*wtmi1j;
    RHOFXSUM(i+1,j)=RHOFXSUM(i+1,j)+rhofluidcur*wtmi1j;
    KXSUM(i+1,j)=KXSUM(i+1,j)+ktotalm*wtmi1j;
    PHIXSUM(i+1,j)=PHIXSUM(i+1,j)+phim(m)*wtmi1j;
    RXSUM(i+1,j)=RXSUM(i+1,j)+etafluidcur/kphim*wtmi1j;
    WTXSUM(i+1,j)=WTXSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    RHOXSUM(i,j+1)=RHOXSUM(i,j+1)+rhototalm*wtmij1;
    RHOFXSUM(i,j+1)=RHOFXSUM(i,j+1)+rhofluidcur*wtmij1;
    KXSUM(i,j+1)=KXSUM(i,j+1)+ktotalm*wtmij1;
    PHIXSUM(i,j+1)=PHIXSUM(i,j+1)+phim(m)*wtmij1;
    RXSUM(i,j+1)=RXSUM(i,j+1)+etafluidcur/kphim*wtmij1;
    WTXSUM(i,j+1)=WTXSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    RHOXSUM(i+1,j+1)=RHOXSUM(i+1,j+1)+rhototalm*wtmi1j1;
    RHOFXSUM(i+1,j+1)=RHOFXSUM(i+1,j+1)+rhofluidcur*wtmi1j1;
    KXSUM(i+1,j+1)=KXSUM(i+1,j+1)+ktotalm*wtmi1j1;
    PHIXSUM(i+1,j+1)=PHIXSUM(i+1,j+1)+phim(m)*wtmi1j1;
    RXSUM(i+1,j+1)=RXSUM(i+1,j+1)+etafluidcur/kphim*wtmi1j1;
    WTXSUM(i+1,j+1)=WTXSUM(i+1,j+1)+wtmi1j1;

    % Interpolation to vy-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xvy(1))/dx)+1;
    i=fix((ym(m)-yvy(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-xvy(j);
    dymi=ym(m)-yvy(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    % i,j Node
    RHOYSUM(i,j)=RHOYSUM(i,j)+rhototalm*wtmij;
    RHOFYSUM(i,j)=RHOFYSUM(i,j)+rhofluidcur*wtmij;
    KYSUM(i,j)=KYSUM(i,j)+ktotalm*wtmij;
    PHIYSUM(i,j)=PHIYSUM(i,j)+phim(m)*wtmij;
    RYSUM(i,j)=RYSUM(i,j)+etafluidcur/kphim*wtmij;
    WTYSUM(i,j)=WTYSUM(i,j)+wtmij;
    % i+1,j Node
    RHOYSUM(i+1,j)=RHOYSUM(i+1,j)+rhototalm*wtmi1j;
    RHOFYSUM(i+1,j)=RHOFYSUM(i+1,j)+rhofluidcur*wtmi1j;
    KYSUM(i+1,j)=KYSUM(i+1,j)+ktotalm*wtmi1j;
    PHIYSUM(i+1,j)=PHIYSUM(i+1,j)+phim(m)*wtmi1j;
    RYSUM(i+1,j)=RYSUM(i+1,j)+etafluidcur/kphim*wtmi1j;
    WTYSUM(i+1,j)=WTYSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    RHOYSUM(i,j+1)=RHOYSUM(i,j+1)+rhototalm*wtmij1;
    RHOFYSUM(i,j+1)=RHOFYSUM(i,j+1)+rhofluidcur*wtmij1;
    KYSUM(i,j+1)=KYSUM(i,j+1)+ktotalm*wtmij1;
    PHIYSUM(i,j+1)=PHIYSUM(i,j+1)+phim(m)*wtmij1;
    RYSUM(i,j+1)=RYSUM(i,j+1)+etafluidcur/kphim*wtmij1;
    WTYSUM(i,j+1)=WTYSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    RHOYSUM(i+1,j+1)=RHOYSUM(i+1,j+1)+rhototalm*wtmi1j1;
    RHOFYSUM(i+1,j+1)=RHOFYSUM(i+1,j+1)+rhofluidcur*wtmi1j1;
    KYSUM(i+1,j+1)=KYSUM(i+1,j+1)+ktotalm*wtmi1j1;
    PHIYSUM(i+1,j+1)=PHIYSUM(i+1,j+1)+phim(m)*wtmi1j1;
    RYSUM(i+1,j+1)=RYSUM(i+1,j+1)+etafluidcur/kphim*wtmi1j1;
    WTYSUM(i+1,j+1)=WTYSUM(i+1,j+1)+wtmi1j1;

    % Interpolation to P-nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    % i,j Node
    GGGPSUM(i,j)=GGGPSUM(i,j)+1/gggtotalm*wtmij;
    SXXSUM(i,j)=SXXSUM(i,j)+sxxm(m)*wtmij;
    RHOSUM(i,j)=RHOSUM(i,j)+rhototalm*wtmij;
    RHOCPSUM(i,j)=RHOCPSUM(i,j)+rhocptotalm*wtmij;
    ALPHASUM(i,j)=ALPHASUM(i,j)+alphasolidm(tm(m))*wtmij;
    ALPHAFSUM(i,j)=ALPHAFSUM(i,j)+alphafluidm(tm(m))*wtmij;
    HRSUM(i,j)=HRSUM(i,j)+hrtotalm*wtmij;
    TKSUM(i,j)=TKSUM(i,j)+tkm(m)*rhocptotalm*wtmij;
    PHISUM(i,j)=PHISUM(i,j)+phim(m)*wtmij;
    WTPSUM(i,j)=WTPSUM(i,j)+wtmij;
    % i+1,j Node
    GGGPSUM(i+1,j)=GGGPSUM(i+1,j)+1/gggtotalm*wtmi1j;
    SXXSUM(i+1,j)=SXXSUM(i+1,j)+sxxm(m)*wtmi1j;
    RHOSUM(i+1,j)=RHOSUM(i+1,j)+rhototalm*wtmi1j;
    RHOCPSUM(i+1,j)=RHOCPSUM(i+1,j)+rhocptotalm*wtmi1j;
    ALPHASUM(i+1,j)=ALPHASUM(i+1,j)+alphasolidm(tm(m))*wtmi1j;
    ALPHAFSUM(i+1,j)=ALPHAFSUM(i+1,j)+alphafluidm(tm(m))*wtmi1j;
    HRSUM(i+1,j)=HRSUM(i+1,j)+hrtotalm*wtmi1j;
    TKSUM(i+1,j)=TKSUM(i+1,j)+tkm(m)*rhocptotalm*wtmi1j;
    PHISUM(i+1,j)=PHISUM(i+1,j)+phim(m)*wtmi1j;
    WTPSUM(i+1,j)=WTPSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    GGGPSUM(i,j+1)=GGGPSUM(i,j+1)+1/gggtotalm*wtmij1;
    SXXSUM(i,j+1)=SXXSUM(i,j+1)+sxxm(m)*wtmij1;
    RHOSUM(i,j+1)=RHOSUM(i,j+1)+rhototalm*wtmij1;
    RHOCPSUM(i,j+1)=RHOCPSUM(i,j+1)+rhocptotalm*wtmij1;
    ALPHASUM(i,j+1)=ALPHASUM(i,j+1)+alphasolidm(tm(m))*wtmij1;
    ALPHAFSUM(i,j+1)=ALPHAFSUM(i,j+1)+alphafluidm(tm(m))*wtmij1;
    HRSUM(i,j+1)=HRSUM(i,j+1)+hrtotalm*wtmij1;
    TKSUM(i,j+1)=TKSUM(i,j+1)+tkm(m)*rhocptotalm*wtmij1;
    PHISUM(i,j+1)=PHISUM(i,j+1)+phim(m)*wtmij1;
    WTPSUM(i,j+1)=WTPSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    GGGPSUM(i+1,j+1)=GGGPSUM(i+1,j+1)+1/gggtotalm*wtmi1j1;
    SXXSUM(i+1,j+1)=SXXSUM(i+1,j+1)+sxxm(m)*wtmi1j1;
    RHOSUM(i+1,j+1)=RHOSUM(i+1,j+1)+rhototalm*wtmi1j1;
    RHOCPSUM(i+1,j+1)=RHOCPSUM(i+1,j+1)+rhocptotalm*wtmi1j1;
    ALPHASUM(i+1,j+1)=ALPHASUM(i+1,j+1)+alphasolidm(tm(m))*wtmi1j1;
    ALPHAFSUM(i+1,j+1)=ALPHAFSUM(i+1,j+1)+alphafluidm(tm(m))*wtmi1j1;
    HRSUM(i+1,j+1)=HRSUM(i+1,j+1)+hrtotalm*wtmi1j1;
    TKSUM(i+1,j+1)=TKSUM(i+1,j+1)+tkm(m)*rhocptotalm*wtmi1j1;
    PHISUM(i+1,j+1)=PHISUM(i+1,j+1)+phim(m)*wtmi1j1;
    WTPSUM(i+1,j+1)=WTPSUM(i+1,j+1)+wtmi1j1;
end
% Compute physical properties
% Basic nodes
YNY=zeros(Ny,Nx);
for j=1:1:Nx
    for i=1:1:Ny
        if(WTSUM(i,j)>0)
            ETA0(i,j)=ETA0SUM(i,j)/WTSUM(i,j);
            ETA(i,j)=ETASUM(i,j)/WTSUM(i,j);
            if(ETA(i,j)<ETA0(i,j))
                YNY(i,j)=1;
            end
            GGG(i,j)=1/(GGGSUM(i,j)/WTSUM(i,j));
            SXY0(i,j)=SXYSUM(i,j)/WTSUM(i,j);
            COH(i,j)=COHSUM(i,j)/WTSUM(i,j);
            TEN(i,j)=TENSUM(i,j)/WTSUM(i,j);
            FRI(i,j)=FRISUM(i,j)/WTSUM(i,j);
        end
    end
end
% Vx-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTXSUM(i,j)>0)
            RHOX(i,j)=RHOXSUM(i,j)/WTXSUM(i,j);
            RHOFX(i,j)=RHOFXSUM(i,j)/WTXSUM(i,j);
            KX(i,j)=KXSUM(i,j)/WTXSUM(i,j);
            PHIX(i,j)=PHIXSUM(i,j)/WTXSUM(i,j);
            RX(i,j)=RXSUM(i,j)/WTXSUM(i,j);
        end
    end
end
% Vy-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTYSUM(i,j)>0)
            RHOY(i,j)=RHOYSUM(i,j)/WTYSUM(i,j);
            RHOFY(i,j)=RHOFYSUM(i,j)/WTYSUM(i,j);
            KY(i,j)=KYSUM(i,j)/WTYSUM(i,j);
            PHIY(i,j)=PHIYSUM(i,j)/WTYSUM(i,j);
            RY(i,j)=RYSUM(i,j)/WTYSUM(i,j);
        end
    end
end
% P-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTPSUM(i,j)>0)
            GGGP(i,j)=1/(GGGPSUM(i,j)/WTPSUM(i,j));
            SXX0(i,j)=SXXSUM(i,j)/WTPSUM(i,j);
            RHO(i,j)=RHOSUM(i,j)/WTPSUM(i,j);
            RHOCP(i,j)=RHOCPSUM(i,j)/WTPSUM(i,j);
            ALPHA(i,j)=ALPHASUM(i,j)/WTPSUM(i,j);
            ALPHAF(i,j)=ALPHAFSUM(i,j)/WTPSUM(i,j);
            HR(i,j)=HRSUM(i,j)/WTPSUM(i,j);
            PHI(i,j)=PHISUM(i,j)/WTPSUM(i,j);
            BETTAPHI(i,j)=1/GGGP(i,j)*PHI(i,j);
            tk1(i,j)=TKSUM(i,j)/RHOCPSUM(i,j);
        end
    end
end
% Applying thermal boundary conditions for interpolated temperature
% Upper boundary 
tk1(1,2:Nx)=tk1(2,2:Nx); % Insulating boundary
% Lower boundary 
tk1(Ny1,2:Nx)=tk1(Ny,2:Nx); % Insulating boundary
% Left boundary
tk1(:,1)=tk1(:,2); % Insulating boundary
% Right boundary
tk1(:,Nx1)=tk1(:,Nx); % Insulating boundary



% Gravity solution
% Composing global matrixes LT(), RT()
% Going through all points of the 2D grid and
% composing respective equations
for j=1:1:Nx1
    for i=1:1:Ny1
        % Define global index in algebraic space
        gk=(j-1)*Ny1+i;
        % Distance from the model centre
        rnode=((xp(j)-xsize/2)^2+(yp(i)-ysize/2)^2)^0.5;
        % External points
        if(rnode>xsize/2 || i==1 || i==Ny1 || j==1 || j==Nx1)
            % Boundary Condition
            % PHI=0
            LP(gk,gk)=1; % Left part
            RP(gk)=0; % Right part
        else
        % Internal points: Temperature eq.
        % d2PHI/dx^2+d2PHI/dy^2=2/3*4*G*pi*RHO
        %          PHI2
        %           |
        %           |
        %  PHI1----PHI3----PHI5
        %           |
        %           |
        %          PHI4
        %
        % Density gradients
        dRHOdx=(RHO(i,j+1)-RHO(i,j-1))/2/dx;
        dRHOdy=(RHO(i+1,j)-RHO(i-1,j))/2/dy;
        % Left part
        LP(gk,gk-Ny1)=1/dx^2; % PHI1
        LP(gk,gk-1)=1/dy^2; % PHI2
        LP(gk,gk)=-2/dx^2-2/dy^2; % PHI3
        LP(gk,gk+1)=1/dy^2; % PHI4
        LP(gk,gk+Ny1)=1/dx^2; % PHI5
        % Right part
        RP(gk)=2/3*4*G*pi*RHO(i,j);
        end
    end
end

% Solving matrixes
SP=LP\RP; % Obtaining algebraic vector of solutions SP()

% Reload solutions SP() to geometrical array PHI()
% Going through all grid points
for j=1:1:Nx1
    for i=1:1:Ny1
        % Compute global index
        gk=(j-1)*Ny1+i;
        % Reload solution
        FI(i,j)=SP(gk);
    end
end
% Compute gravity acceleration
% gx
for j=1:1:Nx
    for i=1:1:Ny1
        % gx=-dPHI/dx
        gx(i,j)=-(FI(i,j+1)-FI(i,j))/dx;
    end
end
% gy
for j=1:1:Nx1
    for i=1:1:Ny
        % gy=-dPHI/dy
        gy(i,j)=-(FI(i+1,j)-FI(i,j))/dy;
    end
end





% Try to increase computational Timestep
dt=min(dt*dtkoefup,dtelastic);

% % Set initial viscoplastic viscosity
% if(timestep==1)
%     ETA=ETA0;
% end

% Save initial viscoplastic viscosity
ETA00=ETA;
% Save initial yielding nodes
YNY00=YNY;

% Start Plastic iterations on Nodes
if (timestep==1)
    BETTAPHI=zeros(Ny1,Nx1); % No elastic compaction for the first timestep
end
for iplast=1:1:nplast

% Recompute viscosity at pressure nodes
for i=2:1:Ny
    for j=2:1:Nx
        ETAP(i,j)=1/((1/ETA(i-1,j-1)+1/ETA(i,j-1)+1/ETA(i-1,j)+1/ETA(i,j))/4);
        ETAPHI(i,j)=etaphikoef*ETAP(i,j)/PHI(i,j);
    end
end
% Hydro-Mechanical Solution
% Composing global matrixes L(), R() for Stokes and continuity equations
for j=1:1:Nx1
    for i=1:1:Ny1
        % Define global indexes in algebraic space
        kvx=((j-1)*Ny1+i-1)*6+1; % Vx solid
        kvy=kvx+1; % Vy solid
        kpm=kvx+2; % Ptotal
        kqx=kvx+3; % qx Darcy
        kqy=kvx+4; % qy Darcy
        kpf=kvx+5; % P fluid
        
        % Vx equation External points
        if(i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1)
            % Boundary Condition 
            % Ghost unknowns 1*Vx=0
            if(j==Nx1)
                L(kvx,kvx)=1; % Left part
                R(kvx)=0; % Right part
            end
            % Left Boundary
            if(j==1)
                L(kvx,kvx)=1; % Left part
                R(kvx)=vxleft; % Right part
            end
            % Right Boundary
            if(j==Nx)
                L(kvx,kvx)=1; % Left part
                R(kvx)=vxright; % Right part
            end
            % Top boundary
            if(i==1 && j>1 && j<Nx)
                L(kvx,kvx)=1; % Left part
                L(kvx,kvx+6)=bctop; % Left part
                R(kvx)=0; % Right part
            end
            % Top boundary
            if(i==Ny1 && j>1 && j<Nx)
                L(kvx,kvx)=1; % Left part
                L(kvx,kvx-6)=bcbottom; % Left part
                R(kvx)=0; % Right part
            end
        else
        % Internal points: x-Stokes eq.
        %            Vx2
        %             |
        %        Vy1  |  Vy3
        %             |
        %     Vx1-P1-Vx3-P2-Vx5
        %             |
        %        Vy2  |  Vy4
        %             |
        %            Vx4
        %
        % Computational viscosity
        ETA1=ETA(i-1,j)*GGG(i-1,j)*dt/(GGG(i-1,j)*dt+ETA(i-1,j));
        ETA2=ETA(i,j)*GGG(i,j)*dt/(GGG(i,j)*dt+ETA(i,j));
        ETAP1=ETAP(i,j)*GGGP(i,j)*dt/(GGGP(i,j)*dt+ETAP(i,j));
        ETAP2=ETAP(i,j+1)*GGGP(i,j+1)*dt/(GGGP(i,j+1)*dt+ETAP(i,j+1));
        % Old stresses
        SXY1=SXY0(i-1,j)*ETA(i-1,j)/(GGG(i-1,j)*dt+ETA(i-1,j));
        SXY2=SXY0(i,j)*ETA(i,j)/(GGG(i,j)*dt+ETA(i,j));
        SXX1=SXX0(i,j)*ETAP(i,j)/(GGGP(i,j)*dt+ETAP(i,j));
        SXX2=SXX0(i,j+1)*ETAP(i,j+1)/(GGGP(i,j+1)*dt+ETAP(i,j+1));
        % Density gradients
        dRHOdx=(RHOX(i,j+1)-RHOX(i,j-1))/2/dx;
        dRHOdy=(RHOX(i+1,j)-RHOX(i-1,j))/2/dy;
        % Left part
        L(kvx,kvx-Ny1*6)=ETAP1/dx^2; % Vx1
        L(kvx,kvx-6)=ETA1/dy^2; % Vx2
        L(kvx,kvx)=-(ETAP1+ETAP2)/dx^2-...
                      (ETA1+ETA2)/dy^2-...
                      dRHOdx*gx(i,j)*dt; % Vx3
        L(kvx,kvx+6)=ETA2/dy^2; % Vx4
        L(kvx,kvx+Ny1*6)=ETAP2/dx^2; % Vx5
        L(kvx,kvy)=ETAP1/dx/dy-ETA2/dx/dy-dRHOdy*gx(i,j)*dt/4;  % Vy2
        L(kvx,kvy+Ny1*6)=-ETAP2/dx/dy+ETA2/dx/dy-dRHOdy*gx(i,j)*dt/4;  % Vy4
        L(kvx,kvy-6)=-ETAP1/dx/dy+ETA1/dx/dy-dRHOdy*gx(i,j)*dt/4;  % Vy1
        L(kvx,kvy+Ny1*6-6)=ETAP2/dx/dy-ETA1/dx/dy-dRHOdy*gx(i,j)*dt/4;  % Vy3
        L(kvx,kpm)=pscale/dx; % P1
        L(kvx,kpm+Ny1*6)=-pscale/dx; % P2
        % Right part
        R(kvx)=-RHOX(i,j)*gx(i,j)-(SXY2-SXY1)/dy-(SXX2-SXX1)/dx;
        end
        
        % Vy equation External points
        if(j==1 || j==Nx1 || i==1 || i==Ny || i==Ny1)
            % Boundary Condition
            % Ghost unknowns 1*Vx=0
            if(i==Ny1)
                L(kvy,kvy)=1; % Left part
                R(kvy)=0; % Right part
            end
            % Top boundary
            if(i==1)
                L(kvy,kvy)=1; % Left part
                R(kvy)=vytop; % Right part
            end
            % Bottom boundary
            if(i==Ny)
                L(kvy,kvy)=1; % Left part
                R(kvy)=vybottom; % Right part
            end
            % Left boundary
            if(j==1 && i>1 && i<Ny)
                L(kvy,kvy)=1; % Left part
                L(kvy,kvy+6*Ny1)=bcleft; % Left part
                R(kvy)=0; % Right part
            end
            % Right boundary
            if(j==Nx1 && i>1 && i<Ny)
                L(kvy,kvy)=1; % Left part
                L(kvy,kvy-6*Ny1)=bcright; % Left part
                R(kvy)=0; % Right part
            end
        else
        % Internal points: y-Stokes eq.
        %            Vy2
        %             |
        %         Vx1 P1 Vx3
        %             |
        %     Vy1----Vy3----Vy5
        %             |
        %         Vx2 P2 Vx4
        %             |
        %            Vy4
        %
        % Computational viscosity
        ETA1=ETA(i,j-1)*GGG(i,j-1)*dt/(GGG(i,j-1)*dt+ETA(i,j-1));
        ETA2=ETA(i,j)*GGG(i,j)*dt/(GGG(i,j)*dt+ETA(i,j));
        ETAP1=ETAP(i,j)*GGGP(i,j)*dt/(GGGP(i,j)*dt+ETAP(i,j));
        ETAP2=ETAP(i+1,j)*GGGP(i+1,j)*dt/(GGGP(i+1,j)*dt+ETAP(i+1,j));
        % Old stresses
        SXY1=SXY0(i,j-1)*ETA(i,j-1)/(GGG(i,j-1)*dt+ETA(i,j-1));
        SXY2=SXY0(i,j)*ETA(i,j)/(GGG(i,j)*dt+ETA(i,j));
        SYY1=-SXX0(i,j)*ETAP(i,j)/(GGGP(i,j)*dt+ETAP(i,j));
        SYY2=-SXX0(i+1,j)*ETAP(i+1,j)/(GGGP(i+1,j)*dt+ETAP(i+1,j));
        % Density gradients
        dRHOdx=(RHOY(i,j+1)-RHOY(i,j-1))/2/dx;
        dRHOdy=(RHOY(i+1,j)-RHOY(i-1,j))/2/dy;
        % Left part
        L(kvy,kvy-Ny1*6)=ETA1/dx^2; % Vy1
        L(kvy,kvy-6)=ETAP1/dy^2; % Vy2
        L(kvy,kvy)=-(ETAP1+ETAP2)/dy^2-...
                      (ETA1+ETA2)/dx^2-...
                      dRHOdy*gy(i,j)*dt; % Vy3
        L(kvy,kvy+6)=ETAP2/dy^2; % Vy4
        L(kvy,kvy+Ny1*6)=ETA2/dx^2; % Vy5
        L(kvy,kvx)=ETAP1/dx/dy-ETA2/dx/dy-dRHOdx*gy(i,j)*dt/4; %Vx3
         L(kvy,kvx+6)=-ETAP2/dx/dy+ETA2/dx/dy-dRHOdx*gy(i,j)*dt/4; %Vx4
        L(kvy,kvx-Ny1*6)=-ETAP1/dx/dy+ETA1/dx/dy-dRHOdx*gy(i,j)*dt/4; %Vx1
        L(kvy,kvx+6-Ny1*6)=ETAP2/dx/dy-ETA1/dx/dy-dRHOdx*gy(i,j)*dt/4; %Vx2
        L(kvy,kpm)=pscale/dy; % P1
        L(kvy,kpm+6)=-pscale/dy; % P2
        
        % Right part
        R(kvy)=-RHOY(i,j)*gy(i,j)-(SXY2-SXY1)/dx-(SYY2-SYY1)/dy;
        end
        
        % P equation External points
        if(i==1 || j==1 || i==Ny1 || j==Nx1)
            % Boundary Condition
            % 1*P=0
            L(kpm,kpm)=1; % Left part
            R(kpm)=0; % Right part
        else
        % Internal points: continuity eq.
        % dVx/dx+dVy/dy=0
        %            Vy1
        %             |
        %        Vx1--P--Vx2
        %             |
        %            Vy2
        %
        % Left part
        L(kpm,kvx-Ny1*6)=-1/dx; % Vx1
        L(kpm,kvx)=1/dx; % Vx2
        L(kpm,kvy-6)=-1/dy; % Vy1
        L(kpm,kvy)=1/dy; % Vy2
        L(kpm,kpm)= pscale/(1-PHI(i,j))*(1/ETAPHI(i,j)+BETTAPHI(i,j)/dt); % Ptotal
        L(kpm,kpf)=-pscale/(1-PHI(i,j))*(1/ETAPHI(i,j)+BETTAPHI(i,j)/dt); % Pfluid
        % Right part
        R(kpm)=(pr0(i,j)-pf0(i,j))/(1-PHI(i,j))*BETTAPHI(i,j)/dt;
        end

        % qxDarcy equation External points
        if(i==1 || i==Ny1 || j==1 || j==Nx || j==Nx1)
            % Boundary Condition
            % 1*qx=0
            L(kqx,kqx)=1; % Left part
            R(kqx)=0; % Right part
            % Top boundary
            if(i==1 && j>1 && j<Nx)
                L(kqx,kqx+6)=bcftop; % Left part
            end
            % Bottom boundary
            if(i==Ny1 && j>1 && j<Nx)
                L(kqx,kqx-6)=bcfbottom; % Left part
            end
        else
        % Internal points: x-Darcy eq.
        % Rx*qxDarcy+dP/dx=RHOfluid*gx
        %     P1-qxD-P2
        % Left part
        L(kqx,kqx)=RX(i,j); % qxD
        L(kqx,kpf)=-pscale/dx; % P1
        L(kqx,kpf+Ny1*6)=pscale/dx; % P2
        % Right part
        R(kqx)=RHOFX(i,j)*gx(i,j);
        end
        
        % qyDarcy equation External points
        if(j==1 || j==Nx1 || i==1 || i==Ny || i==Ny1)
            % Boundary Condition
            % 1*Vy=0
            L(kqy,kqy)=1; % Left part
            R(kqy)=0; % Right part
            % Left boundary
            if(j==1 && i>1 && i<Ny)
                L(kqy,kqy+6*Ny1)=bcfleft; % Left part
            end
            % Right boundary
            if(j==Nx1 && i>1 && i<Ny)
                L(kqy,kqy-6*Ny1)=bcfright; % Left part
            end
        else
        % Internal points: y-Stokes eq.
        % Internal points: x-Darcy eq.
        % Rx*qxDarcy+dP/dx=RHOfluid*gx
        %      P1
        %      |
        %     qxD
        %      |
        %      P2
        % Left part
        L(kqy,kqy)=RY(i,j); % qxD
        L(kqy,kpf)=-pscale/dy; % P1
        L(kqy,kpf+6)=pscale/dy; % P
        % Right part
        R(kqy)=RHOFY(i,j)*gy(i,j);
        end
        
        % Pfluid equation External points
        if(i==1 || j==1 || i==Ny1 || j==Nx1 ||...
          (i==2 && j==2))
            % Boundary Condition
            % 1*Pfluid=0
            L(kpf,kpf)=1; % Left part
            R(kpf)=0; % Right part
            % Real BC
            if(i==2 && j==2)
                L(kpf,kpf)=1*pscale; %Left part
                R(kpf)=psurface; % Right part
            end
        else
        % Internal points: continuity eq.
        % dqxD/dx+dqyD/dy-(Ptotal-Pfluid)/ETHAphi=0
        %            qyD1
        %              |
        %        qxD1--P--qxD2
        %              |
        %            qyD2
        %
        % Left part
        L(kpf,kqx-Ny1*6)=-1/dx; % qxD1
        L(kpf,kqx)=1/dx; % qxD2
        L(kpf,kqy-6)=-1/dy; % qyD1
        L(kpf,kqy)=1/dy; % qyD2
        L(kpf,kpm)=-pscale/(1-PHI(i,j))*(1/ETAPHI(i,j)+BETTAPHI(i,j)/dt); % Ptotal
        L(kpf,kpf)= pscale/(1-PHI(i,j))*(1/ETAPHI(i,j)+BETTAPHI(i,j)/dt); % Pfluid
        % Right part
        R(kpf)=-(pr0(i,j)-pf0(i,j))/(1-PHI(i,j))*BETTAPHI(i,j)/dt;
        end

    end
end

% 4) Solving matrixes, reloading solution
S=L\R; % Obtaining algebraic vector of solutions S()
% Reload solutions S() to vx(), vy(), p()
% Going through all grid points
for j=1:1:Nx1
    for i=1:1:Ny1
        % Define global indexes in algebraic space
        kvx=((j-1)*Ny1+i-1)*6+1; % Vx solid
        kvy=kvx+1; % Vy solid
        kpm=kvx+2; % Ptotal
        kqx=kvx+3; % qx Darcy
        kqy=kvx+4; % qy Darcy
        kpf=kvx+5; % P fluid
        % Reload solution
        vx(i,j)=S(kvx);
        vy(i,j)=S(kvy);
        pr(i,j)=S(kpm)*pscale;
        qxD(i,j)=S(kqx);
        qyD(i,j)=S(kqy);
        pf(i,j)=S(kpf)*pscale;
    end
end
% CHANGE!!!: 
% pavr=(pf(2,2)+pf(2,Nx)+pf(Ny,2)+pf(Ny,Nx))/4;
% pr=pr-pavr+0;
% pf=pf-pavr+0;

% Compute Dln((1-PHI)/PHI)/Dt
APHI=zeros(Ny1,Nx1);
aphimax=0;
for j=2:1:Nx
    for i=2:1:Ny
        APHI(i,j)=((pr(i,j)-pf(i,j))/ETAPHI(i,j)+...
         ((pr(i,j)-pr0(i,j))-(pf(i,j)-pf0(i,j)))/dt*BETTAPHI(i,j))/(1-PHI(i,j))/PHI(i,j);
        aphimax=max(aphimax,abs(APHI(i,j)));
    end
end

% Compute fluid velocity
% Vx fluid
for j=1:1:Nx
    for i=2:1:Ny
        vxf(i,j)=qxD(i,j)/PHIX(i,j);
    end
end
% Apply BC
% Top
vxf(1,:)=-bcftop*vxf(2,:);    
% Bottom
vxf(Ny1,:)=-bcfbottom*vxf(Ny,:);    
% Vy fluid
for j=2:1:Nx
    for i=1:1:Ny
        vyf(i,j)=qyD(i,j)/PHIY(i,j);
    end
end
% Apply BC
% Left
vyf(:,1)=-bcfleft*vyf(:,2);    
% Right
vyf(:,Nx1)=-bcfright*vyf(:,Nx);     
% Add solid velocity
vxf0=vxf; vxf=vxf+vx;
vyf0=vyf; vyf=vyf+vy;

% checking convergence
if timestep==2 || timestep==4 || timestep==10 || timestep==20
    check_hydromech(timestep,xm,ym,tk0,tk1,tk2,DT,APHI,vxf,vyf,pscale,ETA,ETAP,GGG,...
    GGGP,SXY0,SXX0,RHOX,RHOY,RHOFX,RHOFY,RX,RY,ETAPHI,BETTAPHI,PHI,gx,gy,...
    pr0,pf0,dt,R,L,S,vx,vy,qxD,qyD,pr,pf);
end
tk2_jl = load_jld2('tk2', timestep);
RHOX_jl = load_jld2('RHOX', timestep);
maxT_m(timestep) = max(tk2, [], 'all');
maxT_j(timestep) = max(tk2_jl, [], 'all');
meanT_10_m(timestep) = ravg(tk2, 10);
meanT_10_j(timestep) = ravg(tk2_jl, 10);
meanT_40_m(timestep) = ravg(tk2, 40);
meanT_40_j(timestep) = ravg(tk2_jl, 40);
dRHO_m(timestep) = abs(ravg(RHOX, 47)-ravg(RHOX, 1));
dRHO_j(timestep) = abs(ravg(RHOX_jl, 47)-ravg(RHOX_jl, 1));
meanRHO_10_m(timestep) = ravg(RHOX, 10);
meanRHO_10_j(timestep) = ravg(RHOX_jl, 10);
meanRHO_40_m(timestep) = ravg(RHOX, 40);
meanRHO_40_j(timestep) = ravg(RHOX_jl, 40);
% Define displacement timestep dtm
dtm=dt;
maxvx=max(max(abs(vx)));
maxvy=max(max(abs(vy)));
if(dtm*maxvx>dxymax*dx)
    dtm=dxymax*dx/maxvx;
end
if(dtm*maxvy>dxymax*dy)
    dtm=dxymax*dy/maxvy;
end
% Fluid velocity
maxvxf=max(max(abs(vxf)));
maxvyf=max(max(abs(vyf)));
if(dtm*maxvxf>dxymax*dx)
    dtm=dxymax*dx/maxvxf;
end
if(dtm*maxvyf>dxymax*dy)
    dtm=dxymax*dy/maxvyf;
end
% Porosity change
if(aphimax*dtm>dphimax)
    dtm=dphimax/aphimax;
end



% Compute Stress, stress change and strain rate components
% Compute EPSILONxy, SIGMAxy in basic nodes
EXY=zeros(Ny,Nx); % Strain rate EPSILONxy, 1/s
SXY=zeros(Ny,Nx); % Stress SIGMAxy, Pa
DSXY=zeros(Ny,Nx); % Stress change SIGMAxy, Pa
for j=1:1:Nx
    for i=1:1:Ny
        % EXY,SXY, DSXY
        EXY(i,j)=0.5*((vx(i+1,j)-vx(i,j))/dy+...
            (vy(i,j+1)-vy(i,j))/dx);
        SXY(i,j)=2*ETA(i,j)*EXY(i,j)*GGG(i,j)*dtm/(GGG(i,j)*dtm+ETA(i,j))+...
            SXY0(i,j)*ETA(i,j)/(GGG(i,j)*dtm+ETA(i,j));
        DSXY(i,j)=SXY(i,j)-SXY0(i,j);
    end
end
% Compute EPSILONxx, SIGMA'xx in pressure nodes
EXX=zeros(Ny1,Nx1); % Strain rate EPSILONxx, 1/s
EII=zeros(Ny1,Nx1); % Second strain rate invariant, 1/s
SXX=zeros(Ny1,Nx1); % Stress SIGMA'xx, Pa
SII=zeros(Ny1,Nx1); % Second stress invariant, Pa
DSXX=zeros(Ny1,Nx1); % Stress change SIGMA'xx, Pa
DIVV=zeros(Ny1,Nx1); % div(v)
for j=2:1:Nx
    for i=2:1:Ny
        % DIVV
        DIVV(i,j)=(vx(i,j)-vx(i,j-1))/dx+(vy(i,j)-vy(i-1,j))/dy;
        % EXX
        EXX(i,j)=((vx(i,j)-vx(i,j-1))/dx-(vy(i,j)-vy(i-1,j))/dy)/2;
        % SXX
        SXX(i,j)=2*ETAP(i,j)*EXX(i,j)*GGGP(i,j)*dtm/(GGGP(i,j)*dtm+ETAP(i,j))+...
            SXX0(i,j)*ETAP(i,j)/(GGGP(i,j)*dtm+ETAP(i,j));
        DSXX(i,j)=SXX(i,j)-SXX0(i,j);
        % EII
        EII(i,j)=(EXX(i,j)^2+((EXY(i,j)+EXY(i-1,j)+...
            EXY(i,j-1)+EXY(i-1,j-1))/4)^2)^0.5;
        % SII
        SII(i,j)=(SXX(i,j)^2+((SXY(i,j)+SXY(i-1,j)+...
            SXY(i,j-1)+SXY(i-1,j-1))/4)^2)^0.5;
    end
end

% Recompute Dln((1-PHI)/PHI)/Dt
APHI=zeros(Ny1,Nx1);
for j=2:1:Nx
    for i=2:1:Ny
        APHI(i,j)=((pr(i,j)-pf(i,j))/ETAPHI(i,j)+...
         ((pr(i,j)-pr0(i,j))-(pf(i,j)-pf0(i,j)))/dt*BETTAPHI(i,j))/(1-PHI(i,j))/PHI(i,j);
    end
end


% Apply Symmetry to Pressure nodes
% External P-nodes: symmetry
% Top
SXX(1,2:Nx)=SXX(2,2:Nx);
APHI(1,2:Nx)=APHI(2,2:Nx);    
PHI(1,2:Nx)=PHI(2,2:Nx);    
pr(1,2:Nx)=pr(2,2:Nx);    
pf(1,2:Nx)=pf(2,2:Nx);    
% Bottom
SXX(Ny1,2:Nx)=SXX(Ny,2:Nx);
APHI(Ny1,2:Nx)=APHI(Ny,2:Nx);    
PHI(Ny1,2:Nx)=PHI(Ny,2:Nx);    
pr(Ny1,2:Nx)=pr(Ny,2:Nx);    
pf(Ny1,2:Nx)=pf(Ny,2:Nx);    
% Left
SXX(:,1)=SXX(:,2);
APHI(:,1)=APHI(:,2);    
PHI(:,1)=PHI(:,2);    
pr(:,1)=pr(:,2);    
pf(:,1)=pf(:,2);    
% Right
SXX(:,Nx1)=SXX(:,Nx);
APHI(:,Nx1)=APHI(:,Nx);    
PHI(:,Nx1)=PHI(:,Nx);    
pr(:,Nx1)=pr(:,Nx);    
pf(:,Nx1)=pf(:,Nx); 

% Compute solid pressure
ps=(pr-pf.*PHI)./(1-PHI);


% Save nodal stress changes
DSXX0=DSXX;
DSXY0=DSXY;


% Nodal adjusment
% Update viscosity for yielding
% Basic nodes
ETA5=ETA0;
YNY5=zeros(Ny,Nx);
DSY=zeros(Ny,Nx);
ynpl=0;
ddd=0;
for i=1:1:Ny
  for j=1:1:Nx
    % Compute second stress invariant
    % SXX, pt are averaged from four surrounding pressure nodes
%     SIIB=(SXY(i,j)^2+(SXX(i,j)^2+SXX(i+1,j)^2+SXX(i,j+1)^2+SXX(i+1,j+1)^2)/4)^0.5;
    SIIB=(SXY(i,j)^2+((SXX(i,j)+SXX(i+1,j)+SXX(i,j+1)+SXX(i+1,j+1))/4)^2)^0.5;
    % Compute second invariant for a purely elastic stress build-up
    siiel=SIIB*(GGG(i,j)*dt+ETA(i,j))/ETA(i,j);
    % Compute total and fluid pressure
    prB=(pr(i,j)+pr(i+1,j)+pr(i,j+1)+pr(i+1,j+1))/4;
    pfB=(pf(i,j)+pf(i+1,j)+pf(i,j+1)+pf(i+1,j+1))/4;
    % Compute yielding stress
    syieldc=COH(i,j)+FRI(i,j)*(prB-pfB); % Confined fracture
    syieldt=TEN(i,j)+(prB-pfB); % Tensile fracture
    syield=max(min(syieldt,syieldc),0); % Non-negative strength requirement
    % Update error for old yielding nodes
    ynn=0;
    if(YNY(i,j)>0)
        ynn=1;
        DSY(i,j)=SIIB-syield;
        ddd=ddd+DSY(i,j)^2;
        ynpl=ynpl+1;
    end
    % Correcting viscosity for yielding
    if(syield<siiel)
        % New viscosity for the basic node
        etapl=dt*GGG(i,j)*syield/(siiel-syield);
        if(etapl<ETA0(i,j))
            % Recompute nodal visocity
            ETA5(i,j)=etapl^(1-etawt)*ETA(i,j)^etawt;
            % Mark yielding nodes
            YNY5(i,j)=1;
            % Apply viscosity cutoff values
            if(ETA5(i,j)<etamin)
                ETA5(i,j)=etamin;
            elseif(ETA5(i,j)>etamax)
                ETA5(i,j)=etamax;
            end
            % Update Error for new yielding nodes
            if(ynn==0)
                DSY(i,j)=SIIB-syield;
                ddd=ddd+DSY(i,j)^2;
                ynpl=ynpl+1;
            end
        else
            ETA5(i,j)=ETA0(i,j);
            YNY5(i,j)=0;
        end
    else
        ETA5(i,j)=ETA0(i,j);
        YNY5(i,j)=0;
    end
  end
end
% Compute yielding error for markers
if(ynpl>0)
    YERRNOD(iplast)=(ddd/ynpl)^0.5;
end


% Stop iteration
if(ynpl==0 || iplast==nplast || YERRNOD(iplast)<yerrmax)
    break
% Repeat iteration
else
    % Decrease computational timestep if too many iterations
    if(fix(iplast/dtstep)*dtstep==iplast)
        % Decrease timestep
        dt=dt/dtkoef;
        % Reset old viscoplastic viscosity
        ETA=ETA00;
        YNY=YNY00;
    else
        % Use new viscoplastic viscosity
        ETA=ETA5;
        YNY=YNY5;
    end
end
% End Plastic iterations on Nodes
end



% Interpolate new viscoplastic viscosity to markers
for m=1:1:marknum
    % Interpolation viscosity from basic nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-x(1))/dx)+1;
    i=fix((ym(m)-y(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-x(j);
    dymi=ym(m)-y(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Matrix viscosity
    if(tm(m)<3)
        % Rocks
        etasolidcur=etasolidm(tm(m));
        if(tkm(m)>tmsilicate)
            etasolidcur=etasolidmm(tm(m));
        end
        etatotalm=etasolidcur;%*exp(-28*phim(m));
    else
        % Sticky air
        etatotalm=etasolidm(tm(m));
    end
    if(YNY(i,j)>0 || YNY(i+1,j)>0 || YNY(i,j+1)>0 || YNY(i+1,j+1)>0)
%         etavpm(m)=ETA(i,j)*wtmij+ETA(i+1,j)*wtmi1j+...
%                 ETA(i,j+1)*wtmij1+ETA(i+1,j+1)*wtmi1j1;
%         etavpm(m)=1/(1/ETA(i,j)*wtmij+1/ETA(i+1,j)*wtmi1j+...
%                 1/ETA(i,j+1)*wtmij1+1/ETA(i+1,j+1)*wtmi1j1);
        etavpm(m)=1/(YNY(i,j)/ETA(i,j)*wtmij+YNY(i+1,j)/ETA(i+1,j)*wtmi1j+...
                YNY(i,j+1)/ETA(i,j+1)*wtmij1+YNY(i+1,j+1)/ETA(i+1,j+1)*wtmi1j1);
        if(etavpm(m)>=etatotalm)
            etavpm(m)=etatotalm;
        end
    else
        etavpm(m)=etatotalm;
    end
end



% Apply subgrid stress diffusion to markers
if(dsubgrids>0)
SXYSUM=zeros(Ny,Nx);
WTSUM=zeros(Ny,Nx);
SXXSUM=zeros(Ny1,Nx1);
WTPSUM=zeros(Ny1,Nx1);
for m=1:1:marknum
    % SIGMA'xx
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute marker-node SIGMA'xx difference
    dsxxm0=sxxm(m)-(SXX0(i,j)*wtmij+SXX0(i+1,j)*wtmi1j+...
            SXX0(i,j+1)*wtmij1+SXX0(i+1,j+1)*wtmi1j1);
    % Relax stress difference
    dsxxm1=dsxxm0*exp(-dsubgrids*dtm/(etam(tm(m))/gggm(tm(m))));
    % Correct marker stress
    ddsxxm=dsxxm1-dsxxm0;
    sxxm(m)=sxxm(m)+ddsxxm;
    % Update subgrid diffusion on nodes
    % i,j Node
    SXXSUM(i,j)=SXXSUM(i,j)+ddsxxm*wtmij;
    WTPSUM(i,j)=WTPSUM(i,j)+wtmij;
    % i+1,j Node
    SXXSUM(i+1,j)=SXXSUM(i+1,j)+ddsxxm*wtmi1j;
    WTPSUM(i+1,j)=WTPSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    SXXSUM(i,j+1)=SXXSUM(i,j+1)+ddsxxm*wtmij1;
    WTPSUM(i,j+1)=WTPSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    SXXSUM(i+1,j+1)=SXXSUM(i+1,j+1)+ddsxxm*wtmi1j1;
    WTPSUM(i+1,j+1)=WTPSUM(i+1,j+1)+wtmi1j1;

    % SIGMAxy
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-x(1))/dx)+1;
    i=fix((ym(m)-y(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-x(j);
    dymi=ym(m)-y(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute marker-node SIGMAxy difference
    dsxym0=sxym(m)-(SXY0(i,j)*wtmij+SXY0(i+1,j)*wtmi1j+...
            SXY0(i,j+1)*wtmij1+SXY0(i+1,j+1)*wtmi1j1);
    % Relax stress difference
    dsxym1=dsxym0*exp(-dsubgrids*dtm/(etam(tm(m))/gggm(tm(m))));
    % Correct marker stress
    ddsxym=dsxym1-dsxym0;
    sxym(m)=sxym(m)+ddsxym;
    % Update subgrid diffusion on nodes
    % i,j Node
    SXYSUM(i,j)=SXYSUM(i,j)+ddsxym*wtmij;
    WTSUM(i,j)=WTSUM(i,j)+wtmij;
    % i+1,j Node
    SXYSUM(i+1,j)=SXYSUM(i+1,j)+ddsxym*wtmi1j;
    WTSUM(i+1,j)=WTSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    SXYSUM(i,j+1)=SXYSUM(i,j+1)+ddsxym*wtmij1;
    WTSUM(i,j+1)=WTSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    SXYSUM(i+1,j+1)=SXYSUM(i+1,j+1)+ddsxym*wtmi1j1;
    WTSUM(i+1,j+1)=WTSUM(i+1,j+1)+wtmi1j1;
end
% Compute DSXXsubgrid
DSXXsubgrid=zeros(Ny1,Nx1);
% P-nodes
for j=2:1:Nx
    for i=2:1:Ny
        if(WTPSUM(i,j)>0)
            DSXXsubgrid(i,j)=SXXSUM(i,j)/WTPSUM(i,j);
        end
    end
end
% Correct DSXX
DSXX=DSXX-DSXXsubgrid;
% Compute DSXYsubgrid
DSXYsubgrid=zeros(Ny,Nx);
% Basic nodes
for j=1:1:Nx
    for i=1:1:Ny
        if(WTSUM(i,j)>0)
            DSXYsubgrid(i,j)=SXYSUM(i,j)/WTSUM(i,j);
        end
    end
end
% Correct DSXY
DSXY=DSXY-DSXYsubgrid;
end

% Interpolate DSXX, DSXY to markers
for m=1:1:marknum
    % SIGMA'xx
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<2)
        j=2;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<2)
        i=2;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update marker by SIGMA'xx change 
    sxxm(m)=sxxm(m)+(DSXX(i,j)*wtmij+DSXX(i+1,j)*wtmi1j+...
            DSXX(i,j+1)*wtmij1+DSXX(i+1,j+1)*wtmi1j1);

    % SIGMAxy
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-x(1))/dx)+1;
    i=fix((ym(m)-y(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-x(j);
    dymi=ym(m)-y(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update marker by SIGMA'xx change 
    sxym(m)=sxym(m)+(DSXY(i,j)*wtmij+DSXY(i+1,j)*wtmi1j+...
            DSXY(i,j+1)*wtmij1+DSXY(i+1,j+1)*wtmi1j1);
end


% Compute shear heating HS in Temperature/Pressure nodes
HS=zeros(Ny1,Nx1); % Adiabatic heating, W/m^3
for j=2:1:Nx
    for i=2:1:Ny
        % Average SXY*EXY
        SXYEXY=(SXY(i,j)^2/ETA(i,j)+SXY(i-1,j)^2/ETA(i-1,j)+...
            SXY(i,j-1)^2/ETA(i,j-1)+SXY(i-1,j-1)^2/ETA(i-1,j-1))/4;
        % HS
        HS(i,j)=SXX(i,j)^2/ETAP(i,j)+SXYEXY+...
            (pr(i,j)-pf(i,j))^2/(1-PHI(i,j))/ETAPHI(i,j)+...
            (RX(i,j-1)*qxD(i,j-1)^2+RX(i,j)*qxD(i,j)^2)/2+...
            (RY(i-1,j)*qyD(i-1,j)^2+RY(i,j)*qyD(i,j)^2)/2;
    end
end

% Compute adiabatic heating HA in Temperature/Pressure nodes
% Compute solid pressure
if (timestep==1)
    pr0=pr; % No total pressure change for the first timestep
    pf0=pf; % No fluid pressure change for the first timestep
    ps0=pf; % No solid pressure change for the first timestep
end
% Old solid pressure
HA=zeros(Ny1,Nx1); % Shear heating, W/m^3
for j=2:1:Nx
    for i=2:1:Ny
        % HA
        % Indirect calculation of dpdt
        % Average vy, vx, vxf, vyf
        VXP=(vx(i,j)+vx(i,j-1))/2;
        VYP=(vy(i,j)+vy(i-1,j))/2;
        VXFP=(vxf(i,j)+vxf(i,j-1))/2;
        VYFP=(vyf(i,j)+vyf(i-1,j))/2;
        % Evaluate DPsolid/Dt with upwind differences
        if(VXP<0)
            dpsdx=(ps(i,j)-ps(i,j-1))/dx;
        else
            dpsdx=(ps(i,j+1)-ps(i,j))/dx;
        end
        if(VYP<0)
            dpsdy=(ps(i,j)-ps(i-1,j))/dy;
        else
            dpsdy=(ps(i+1,j)-ps(i,j))/dy;
        end
        dpsdt=VXP*dpsdx+VYP*dpsdy;
        % Evaluate DPfluid/Dt with upwind differences
        if(VXFP>0)
            dpfdx=(pf(i,j)-pf(i,j-1))/dx;
        else
            dpfdx=(pf(i,j+1)-pf(i,j))/dx;
        end
        if(VYFP>0)
            dpfdy=(pf(i,j)-pf(i-1,j))/dy;
        else
            dpfdy=(pf(i+1,j)-pf(i,j))/dy;
        end
        dpfdt=VXFP*dpsdx+VYFP*dpsdy;
%         % Direct calculation of dpdt
%         dpsdt=(ps(i,j)-ps0(i,j))/dt;
%         dpfdt=(pf(i,j)-pf0(i,j))/dt;
        % HA
        HA(i,j)=(1-PHI(i,j))*tk1(i,j)*ALPHA(i,j)*dpsdt+...
            PHI(i,j)*tk1(i,j)*ALPHAF(i,j)*dpfdt;
    end
end


% Thermal iterations
tk0=tk1;
dtt=dtm;
dttsum=0;
titer=1;
while(dttsum<dtm)
% Composing global matrixes LT(), RT()
% Going through all points of the 2D grid and
% composing respective equations
for j=1:1:Nx1
    for i=1:1:Ny1
        % Define global index in algebraic space
        gk=(j-1)*Ny1+i;
        % External points
        if(i==1 || i==Ny1 || j==1 || j==Nx1)
            % Boundary Condition
            % Top BC: T=273
            if(i==1 && j>1 && j<Nx1)
                LT(gk,gk)=1; % Left part
                LT(gk,gk+1)=-1; % Left part
                RT(gk)=0; % Right part
            end
            % Bottom BC: T=1500
            if(i==Ny1 && j>1 && j<Nx1)
                LT(gk,gk)=1; % Left part
                LT(gk,gk-1)=-1; % Left part
                RT(gk)=0; % Right part
            end
            % Left BC: dT/dx=0
            if(j==1)
                LT(gk,gk)=1; % Left part
                LT(gk,gk+Ny1)=-1; % Left part
                RT(gk)=0; % Right part
            end
            % Right BC: dT/dx=0
            if(j==Nx1)
                LT(gk,gk)=1; % Left part
                LT(gk,gk-Ny1)=-1; % Left part
                RT(gk)=0; % Right part
            end
        else
        % Internal points: Temperature eq.
        % RHO*CP*dT/dt=-dqx/dx-dqy/dy+Hr+Hs+Ha
        %          Tdt2
        %           |
        %          Ky1
        %           |
        %Tdt1-Kx1-T03,Tdt3-Kx2-Tdt5
        %           |
        %          Ky2
        %           |
        %          Tdt4
        %
        % Left part
        Kx1=KX(i,j-1); 
        Kx2=KX(i,j); 
        Ky1=KY(i-1,j); 
        Ky2=KY(i,j); 
        LT(gk,gk-Ny1)=-Kx1/dx^2; % T1
        LT(gk,gk-1)=-Ky1/dy^2; % FI2
        LT(gk,gk)=RHOCP(i,j)/dtt+(Kx1+Kx2)/dx^2+(Ky1+Ky2)/dy^2; % FI3
        LT(gk,gk+1)=-Ky2/dy^2; % FI4
        LT(gk,gk+Ny1)=-Kx2/dx^2; % FI5
        % Right part
        RT(gk)=RHOCP(i,j)/dtt*tk1(i,j)+HR(i,j)+HA(i,j)+HS(i,j);
        end
    end
end

% Solving matrixes
ST=LT\RT; % Obtaining algebraic vector of solutions ST()

% Reload solutions ST() to geometrical array Tdt()
% Going through all grid points
for j=1:1:Nx1
    for i=1:1:Ny1
        % Compute global index
        gk=(j-1)*Ny1+i;
        % Reload solution
        tk2(i,j)=ST(gk);
    end
end
% Compute DT
DT=tk2-tk1;
titer
dtt
if(titer==1)
    % Apply thermal timestepping condition
    maxDTcurrent=max(max(abs(DT)));
    if(maxDTcurrent>DTmax)
        dtt=dtt/maxDTcurrent*DTmax;
    else
        dttsum=dttsum+dtt; % Update dttsum
    end
else
    dttsum=dttsum+dtt; % Update dttsum
    % Adjust timestep
    if(dtt>dtm-dttsum)
        dtt=dtm-dttsum;
    end
end

dttsum

titer=titer+1; % Update iteration counter
end
% Compute/save overall temperature changes
DT=tk2-tk0;
DT0=DT;

% Apply subgrid temperature diffusion on markers
if(dsubgridt>0)
TKSUM=zeros(Ny1,Nx1);
RHOCPSUM=zeros(Ny1,Nx1);
for m=1:1:marknum
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute marker-node T difference
    dtkm0=tkm(m)-(tk1(i,j)*wtmij+tk1(i+1,j)*wtmi1j+...
            tk1(i,j+1)*wtmij1+tk1(i+1,j+1)*wtmi1j1);
    % Compute marker parameters
    if(tm(m)<3)
        % Rocks
        rhocptotalm=rhocpsolidm(tm(m))*(1-phim(m))+rhocpfluid*phim(m);
        ktotalm=(ksolidm(tm(m))*kfluid/2+((ksolidm(tm(m))*(3*phim(m)-2)+...
            kfluid*(1-3*phim(m)))^2)/16)^0.5-(ksolidm(tm(m))*(3*phim(m)-2)+...
            kfluid*(1-3*phim(m)))/4;
    else
        % Sticky air
        rhocptotalm=rhocpsolidm(tm(m));
        ktotalm=ksolidm(tm(m));
    end    % Relax temperature difference
    dtkm1=dtkm0*exp(-dsubgridt*ktotalm*dtm/rhocptotalm*(2/dx^2+2/dy^2));
    % Correct marker temperature
    ddtkm=dtkm1-dtkm0;
    tkm(m)=tkm(m)+ddtkm;
    % Update subgrid diffusion on nodes
    % i,j Node
    TKSUM(i,j)=TKSUM(i,j)+ddtkm*rhocpm(tm(m))*wtmij;
    RHOCPSUM(i,j)=RHOCPSUM(i,j)+rhocpm(tm(m))*wtmij;
    % i+1,j Node
    TKSUM(i+1,j)=TKSUM(i+1,j)+ddtkm*rhocpm(tm(m))*wtmi1j;
    RHOCPSUM(i+1,j)=RHOCPSUM(i+1,j)+rhocpm(tm(m))*wtmi1j;
    % i,j+1 Node
    TKSUM(i,j+1)=TKSUM(i,j+1)+ddtkm*rhocpm(tm(m))*wtmij1;
    RHOCPSUM(i,j+1)=RHOCPSUM(i,j+1)+rhocpm(tm(m))*wtmij1;
    % i+1,j+1 Node
    TKSUM(i+1,j+1)=TKSUM(i+1,j+1)+ddtkm*rhocpm(tm(m))*wtmi1j1;
    RHOCPSUM(i+1,j+1)=RHOCPSUM(i+1,j+1)+rhocpm(tm(m))*wtmi1j1;
end
% Compute DTsubgrid
DTsubgrid=zeros(Ny1,Nx1);
% P-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(RHOCPSUM(i,j)>0)
            DTsubgrid(i,j)=TKSUM(i,j)/RHOCPSUM(i,j);
        end
    end
end
% Correct DT
DT=DT-DTsubgrid;
end


% Interpolate DT to markers
for m=1:1:marknum
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Update properties
    tkm(m)=tkm(m)+DT(i,j)*wtmij+DT(i+1,j)*wtmi1j+...
            DT(i,j+1)*wtmij1+DT(i+1,j+1)*wtmi1j1;
    % Interpolate tk2 at 1st timestep
    if(timestep==1)
        tkm(m)=tk2(i,j)*wtmij+tk2(i+1,j)*wtmi1j+...
            tk2(i,j+1)*wtmij1+tk2(i+1,j+1)*wtmi1j1;
    end
end


    
% Update porosity on markers (not for sticky air)
for m=1:1:marknum
    if(tm(m)<3)
        % Interpolate APHI
        % Define i,j indexes for the upper left node
        j=fix((xm(m)-xp(1))/dx)+1;
        i=fix((ym(m)-yp(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx)
            j=Nx;
        end
        if(i<1)
            i=1;
        elseif(i>Ny)
            i=Ny;
        end
        % Compute distances
        dxmj=xm(m)-xp(j);
        dymi=ym(m)-yp(i);
        % Compute weights
        wtmij=(1-dxmj/dx)*(1-dymi/dy);
        wtmi1j=(1-dxmj/dx)*(dymi/dy);    
        wtmij1=(dxmj/dx)*(1-dymi/dy);
        wtmi1j1=(dxmj/dx)*(dymi/dy);
        % Compute Dln((1-PHI)/PHI)/Dt
        aphim=APHI(i,j)*wtmij+APHI(i+1,j)*wtmi1j+...
            APHI(i,j+1)*wtmij1+APHI(i+1,j+1)*wtmi1j1;
        % Change Porosity
        phim(m)=phim(m)/((1-phim(m))*exp(aphim*dtm)+phim(m));
        if(phim(m)<phimin)
            phim(m)=phimin;
        elseif(phim(m)>phimax)
            phim(m)=phimax;
        end
    end
end

% Compute fluid velocity in pressure nodes
% vxpf
for j=2:1:Nx
    for i=2:1:Ny
        vxpf(i,j)=(vxf(i,j)+vxf(i,j-1))/2;
    end
end
% Apply BC
% Top
vxpf(1,2:Nx-1)=-bcftop*vxpf(2,2:Nx-1);    
% Bottom
vxpf(Ny1,2:Nx-1)=-bcfbottom*vxpf(Ny,2:Nx-1);    
% Left
vxpf(:,1)=2*vxleft-vxpf(:,2);
% Right
vxpf(:,Nx1)=2*vxright-vxpf(:,Nx);
% vypf
for j=2:1:Nx
    for i=2:1:Ny
        vypf(i,j)=(vyf(i,j)+vyf(i-1,j))/2;
    end
end    
% Apply BC
% Left
vypf(2:Ny-1,1)=-bcfleft*vypf(2:Ny-1,2);    
% Right
vypf(2:Ny-1,Nx1)=-bcfright*vypf(2:Ny-1,Nx); % Free slip    
% Top
vypf(1,:)=2*vytop-vypf(2,:);
% Bottom
vypf(Ny1,:)=2*vybottom-vypf(Ny,:);

% Compute velocity in pressure nodes
% vx
for j=2:1:Nx
    for i=2:1:Ny
        vxp(i,j)=(vx(i,j)+vx(i,j-1))/2;
    end
end
% Apply BC
% Top
vxp(1,2:Nx-1)=-bctop*vxp(2,2:Nx-1);    
% Bottom
vxp(Ny1,2:Nx-1)=-bcbottom*vxp(Ny,2:Nx-1);    
% Left
vxp(:,1)=2*vxleft-vxp(:,2);
% Right
vxp(:,Nx1)=2*vxright-vxp(:,Nx);
% vy
for j=2:1:Nx
    for i=2:1:Ny
        vyp(i,j)=(vy(i,j)+vy(i-1,j))/2;
    end
end    
% Apply BC
% Left
vyp(2:Ny-1,1)=-bcleft*vyp(2:Ny-1,2);    
% Right
vyp(2:Ny-1,Nx1)=-bcright*vyp(2:Ny-1,Nx); % Free slip    
% Top
vyp(1,:)=2*vytop-vyp(2,:);
% Bottom
vyp(Ny1,:)=2*vybottom-vyp(Ny,:);

% Compute rotation rate wyx=1/2(dVy/dx-dVx/dy) for basic nodes
for i=1:1:Ny
    for j=1:1:Nx
        wyx(i,j)=0.5*((vy(i,j+1)-vy(i,j))/dx-(vx(i+1,j)-vx(i,j))/dy);
    end
end

% Move markers with 4th order Runge-Kutta
vxm=zeros(4,1);
vym=zeros(4,1);
for m=1:1:marknum
    
    % Interpolate solid temperature for the initial marker location
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xp(1))/dx)+1;
    i=fix((ym(m)-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xm(m)-xp(j);
    dymi=ym(m)-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute Tsolid
    tksm0=tk2(i,j)*wtmij+tk2(i+1,j)*wtmi1j+...
            tk2(i,j+1)*wtmij1+tk2(i+1,j+1)*wtmi1j1;    
        
    % Interpolate local rotation rate
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-x(1))/dx)+1;
    i=fix((ym(m)-y(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx-1)
        j=Nx-1;
    end
    if(i<1)
        i=1;
    elseif(i>Ny-1)
        i=Ny-1;
    end
    % Compute distances
    dxmj=xm(m)-x(j);
    dymi=ym(m)-y(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute vx velocity
    omegam=wyx(i,j)*wtmij+wyx(i+1,j)*wtmi1j+...
    wyx(i,j+1)*wtmij1+wyx(i+1,j+1)*wtmi1j1;
    % Analytical stress rotation using SIGMA'xx=-SIGMA'yy
    THETA=dtm*omegam; % Incremental rotation angle
    sxxmnew=sxxm(m)*cos(THETA)^2-sxxm(m)*sin(THETA)^2-sxym(m)*sin(2*THETA);
    sxymnew=sxxm(m)*sin(2*THETA)+sxym(m)*cos(2*THETA);
    sxxm(m)=sxxmnew; sxym(m)=sxymnew;    
    
    % Save initial marker coordinates
    xA=xm(m);
    yA=ym(m);
    for rk=1:1:4
        % Interpolate vx
        % Define i,j indexes for the upper left node
        j=fix((xm(m)-xvx(1))/dx)+1;
        i=fix((ym(m)-yvx(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx-1)
            j=Nx-1;
        end
        if(i<1)
            i=1;
        elseif(i>Ny)
            i=Ny;
        end
        % Compute distances
        dxmj=xm(m)-xvx(j);
        dymi=ym(m)-yvx(i);
        % Compute weights
        % Compute vx velocity for the top and bottom of the cell
        vxm13=vx(i,j)*(1-dxmj/dx)+vx(i,j+1)*dxmj/dx;
        vxm24=vx(i+1,j)*(1-dxmj/dx)+vx(i+1,j+1)*dxmj/dx;
        % Compute correction
        if(dxmj/dx>=0.5)
            if(j<Nx-1)
                vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(vx(i,j)-2*vx(i,j+1)+vx(i,j+2));
                vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(vx(i+1,j)-2*vx(i+1,j+1)+vx(i+1,j+2));
            end
        else
            if(j>1)
                vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(vx(i,j-1)-2*vx(i,j)+vx(i,j+1));
                vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(vx(i+1,j-1)-2*vx(i+1,j)+vx(i+1,j+1));
            end
        end
        % Compute vx
        vxm(rk)=(1-dymi/dy)*vxm13+(dymi/dy)*vxm24;
        
        % Interpolate vy
        % Define i,j indexes for the upper left node
        j=fix((xm(m)-xvy(1))/dx)+1;
        i=fix((ym(m)-yvy(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx)
            j=Nx;
        end
        if(i<1)
            i=1;
        elseif(i>Ny-1)
            i=Ny-1;
        end
        % Compute distances
        dxmj=xm(m)-xvy(j);
        dymi=ym(m)-yvy(i);
        % Compute weights
        % Compute vy velocity for the left and right of the cell
        vym12=vy(i,j)*(1-dymi/dy)+vy(i+1,j)*dymi/dy;
        vym34=vy(i,j+1)*(1-dymi/dy)+vy(i+1,j+1)*dymi/dy;
        % Compute correction
        if(dymi/dy>=0.5)
            if(i<Ny-1)
                vym12=vym12+1/2*((dymi/dy-0.5)^2)*(vy(i,j)-2*vy(i+1,j)+vy(i+2,j));
                vym34=vym34+1/2*((dymi/dy-0.5)^2)*(vy(i,j+1)-2*vy(i+1,j+1)+vy(i+2,j+1));
            end      
        else
            if(i>1)
                vym12=vym12+1/2*((dymi/dy-0.5)^2)*(vy(i-1,j)-2*vy(i,j)+vy(i+1,j));
                vym34=vym34+1/2*((dymi/dy-0.5)^2)*(vy(i-1,j+1)-2*vy(i,j+1)+vy(i+1,j+1));
            end
        end
        % Compute vy
        vym(rk)=(1-dxmj/dx)*vym12+(dxmj/dx)*vym34;
        
        % Change coordinates to obtain B,C,D points
        if(rk==1 || rk==2)
            xm(m)=xA+dtm/2*vxm(rk);
            ym(m)=yA+dtm/2*vym(rk);
        elseif(rk==3)
            xm(m)=xA+dtm*vxm(rk);
            ym(m)=yA+dtm*vym(rk);
        end
    end
    % Restore initial coordinates
    xm(m)=xA;
    ym(m)=yA;
    % Compute effective velocity
    vxmeff=1/6*(vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4));
    vymeff=1/6*(vym(1)+2*vym(2)+2*vym(3)+vym(4));
    % Move markers
    xm(m)=xm(m)+dtm*vxmeff;
    ym(m)=ym(m)+dtm*vymeff;
    
    % Backtracing markers with fluid velocity
    xcur=xm(m);
    ycur=ym(m);
    xA=xcur;
    yA=ycur;
    for rk=1:1:4
        % Interpolate vx
        % Define i,j indexes for the upper left node
        j=fix((xcur-xvx(1))/dx)+1;
        i=fix((ycur-yvx(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx-1)
            j=Nx-1;
        end
        if(i<1)
            i=1;
        elseif(i>Ny)
            i=Ny;
        end
        % Compute distances
        dxmj=xcur-xvx(j);
        dymi=ycur-yvx(i);
        % Compute weights
        % Compute vx velocity for the top and bottom of the cell
        vxm13=vxf(i,j)*(1-dxmj/dx)+vxf(i,j+1)*dxmj/dx;
        vxm24=vxf(i+1,j)*(1-dxmj/dx)+vxf(i+1,j+1)*dxmj/dx;
        % Compute correction
        if(dxmj/dx>=0.5)
            if(j<Nx-1)
                vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(vxf(i,j)-2*vxf(i,j+1)+vxf(i,j+2));
                vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(vxf(i+1,j)-2*vxf(i+1,j+1)+vxf(i+1,j+2));
            end
        else
            if(j>1)
                vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(vxf(i,j-1)-2*vxf(i,j)+vxf(i,j+1));
                vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(vxf(i+1,j-1)-2*vxf(i+1,j)+vxf(i+1,j+1));
            end
        end
        % Compute vx
        vxm(rk)=(1-dymi/dy)*vxm13+(dymi/dy)*vxm24;
        
        % Interpolate vy
        % Define i,j indexes for the upper left node
        j=fix((xcur-xvy(1))/dx)+1;
        i=fix((ycur-yvy(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx)
            j=Nx;
        end
        if(i<1)
            i=1;
        elseif(i>Ny-1)
            i=Ny-1;
        end
        % Compute distances
        dxmj=xcur-xvy(j);
        dymi=ycur-yvy(i);
        % Compute weights
        % Compute vy velocity for the left and right of the cell
        vym12=vyf(i,j)*(1-dymi/dy)+vyf(i+1,j)*dymi/dy;
        vym34=vyf(i,j+1)*(1-dymi/dy)+vyf(i+1,j+1)*dymi/dy;
        % Compute correction
        if(dymi/dy>=0.5)
            if(i<Ny-1)
                vym12=vym12+1/2*((dymi/dy-0.5)^2)*(vyf(i,j)-2*vyf(i+1,j)+vyf(i+2,j));
                vym34=vym34+1/2*((dymi/dy-0.5)^2)*(vyf(i,j+1)-2*vyf(i+1,j+1)+vyf(i+2,j+1));
            end      
        else
            if(i>1)
                vym12=vym12+1/2*((dymi/dy-0.5)^2)*(vyf(i-1,j)-2*vyf(i,j)+vyf(i+1,j));
                vym34=vym34+1/2*((dymi/dy-0.5)^2)*(vyf(i-1,j+1)-2*vyf(i,j+1)+vyf(i+1,j+1));
            end
        end
        % Compute vy
        vym(rk)=(1-dxmj/dx)*vym12+(dxmj/dx)*vym34;
        
        % Change coordinates to obtain B,C,D points
        if(rk==1 || rk==2)
            xcur=xA-dtm/2*vxm(rk);
            ycur=yA-dtm/2*vym(rk);
        elseif(rk==3)
            xcur=xA-dtm*vxm(rk);
            ycur=yA-dtm*vym(rk);
        end
    end
    % Compute effective velocity
    vxmeff=1/6*(vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4));
    vymeff=1/6*(vym(1)+2*vym(2)+2*vym(3)+vym(4));
    % Trace the node backward
    xcur=xA-dtm*vxmeff;
    ycur=yA-dtm*vymeff;
    % Interpolate fluid temperature
    % Define i,j indexes for the upper left node
    j=fix((xcur-xp(1))/dx)+1;
    i=fix((ycur-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xcur-xp(j);
    dymi=ycur-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute nodal Tfluid
    tkfm0=tk2(i,j)*wtmij+tk2(i+1,j)*wtmi1j+...
            tk2(i,j+1)*wtmij1+tk2(i+1,j+1)*wtmi1j1;
    % Compute Tfluid-Tsolid for the marker
    dtkfsm=tkfm0-tksm0;
    % Correct marker temperature
    tkm(m)=((1-phim(m))*tkm(m)*rhocpsolidm(tm(m))+...
        phim(m)*(tkm(m)+dtkfsm)*rhocpfluidm(tm(m)))/...
        ((1-phim(m))*rhocpsolidm(tm(m))+phim(m)*rhocpfluidm(tm(m)));
end  

% Backtracing Pressure nodes: Ptotal
% Backtracing is based on 4th order Runge-Kutta
vxm=zeros(4,1);
vym=zeros(4,1);
pr0=pr;
ps0=ps;
for jj=2:1:Nx
for ii=2:1:Ny
    % Save initial nodal coordinates
    xcur=xp(jj);
    ycur=yp(ii);
    xA=xcur;
    yA=ycur;
    for rk=1:1:4
        % Interpolate vx
        % Define i,j indexes for the upper left node
        j=fix((xcur-xvx(1))/dx)+1;
        i=fix((ycur-yvx(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx-1)
            j=Nx-1;
        end
        if(i<1)
            i=1;
        elseif(i>Ny)
            i=Ny;
        end
        % Compute distances
        dxmj=xcur-xvx(j);
        dymi=ycur-yvx(i);
        % Compute weights
        % Compute vx velocity for the top and bottom of the cell
        vxm13=vx(i,j)*(1-dxmj/dx)+vx(i,j+1)*dxmj/dx;
        vxm24=vx(i+1,j)*(1-dxmj/dx)+vx(i+1,j+1)*dxmj/dx;
        % Compute correction
        if(dxmj/dx>=0.5)
            if(j<Nx-1)
                vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(vx(i,j)-2*vx(i,j+1)+vx(i,j+2));
                vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(vx(i+1,j)-2*vx(i+1,j+1)+vx(i+1,j+2));
            end
        else
            if(j>1)
                vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(vx(i,j-1)-2*vx(i,j)+vx(i,j+1));
                vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(vx(i+1,j-1)-2*vx(i+1,j)+vx(i+1,j+1));
            end
        end
        % Compute vx
        vxm(rk)=(1-dymi/dy)*vxm13+(dymi/dy)*vxm24;
        
        % Interpolate vy
        % Define i,j indexes for the upper left node
        j=fix((xcur-xvy(1))/dx)+1;
        i=fix((ycur-yvy(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx)
            j=Nx;
        end
        if(i<1)
            i=1;
        elseif(i>Ny-1)
            i=Ny-1;
        end
        % Compute distances
        dxmj=xcur-xvy(j);
        dymi=ycur-yvy(i);
        % Compute weights
        % Compute vy velocity for the left and right of the cell
        vym12=vy(i,j)*(1-dymi/dy)+vy(i+1,j)*dymi/dy;
        vym34=vy(i,j+1)*(1-dymi/dy)+vy(i+1,j+1)*dymi/dy;
        % Compute correction
        if(dymi/dy>=0.5)
            if(i<Ny-1)
                vym12=vym12+1/2*((dymi/dy-0.5)^2)*(vy(i,j)-2*vy(i+1,j)+vy(i+2,j));
                vym34=vym34+1/2*((dymi/dy-0.5)^2)*(vy(i,j+1)-2*vy(i+1,j+1)+vy(i+2,j+1));
            end      
        else
            if(i>1)
                vym12=vym12+1/2*((dymi/dy-0.5)^2)*(vy(i-1,j)-2*vy(i,j)+vy(i+1,j));
                vym34=vym34+1/2*((dymi/dy-0.5)^2)*(vy(i-1,j+1)-2*vy(i,j+1)+vy(i+1,j+1));
            end
        end
        % Compute vy
        vym(rk)=(1-dxmj/dx)*vym12+(dxmj/dx)*vym34;
        
        % Change coordinates to obtain B,C,D points
        if(rk==1 || rk==2)
            xcur=xA-dtm/2*vxm(rk);
            ycur=yA-dtm/2*vym(rk);
        elseif(rk==3)
            xcur=xA-dtm*vxm(rk);
            ycur=yA-dtm*vym(rk);
        end
    end
    % Compute effective velocity
    vxmeff=1/6*(vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4));
    vymeff=1/6*(vym(1)+2*vym(2)+2*vym(3)+vym(4));
    % Trace the node backward
    xcur=xA-dtm*vxmeff;
    ycur=yA-dtm*vymeff;
    % Interpolate nodal property
    % SIGMA'xx, P
    % Define i,j indexes for the upper left node
    j=fix((xcur-xp(1))/dx)+1;
    i=fix((ycur-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xcur-xp(j);
    dymi=ycur-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute nodal total pressure
    pr0(ii,jj)=pr(i,j)*wtmij+pr(i+1,j)*wtmi1j+...
            pr(i,j+1)*wtmij1+pr(i+1,j+1)*wtmi1j1;
    ps0(ii,jj)=ps(i,j)*wtmij+ps(i+1,j)*wtmi1j+...
            ps(i,j+1)*wtmij1+ps(i+1,j+1)*wtmi1j1;
end
end

% Backtracing Pressure nodes: Pfluid
pf0=pf;
for jj=2:1:Nx
for ii=2:1:Ny
    % Save initial nodal coordinates
    xcur=xp(jj);
    ycur=yp(ii);
    xA=xcur;
    yA=ycur;
    for rk=1:1:4
        % Interpolate vx
        % Define i,j indexes for the upper left node
        j=fix((xcur-xvx(1))/dx)+1;
        i=fix((ycur-yvx(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx-1)
            j=Nx-1;
        end
        if(i<1)
            i=1;
        elseif(i>Ny)
            i=Ny;
        end
        % Compute distances
        dxmj=xcur-xvx(j);
        dymi=ycur-yvx(i);
        % Compute weights
        % Compute vx velocity for the top and bottom of the cell
        vxm13=vxf(i,j)*(1-dxmj/dx)+vxf(i,j+1)*dxmj/dx;
        vxm24=vxf(i+1,j)*(1-dxmj/dx)+vxf(i+1,j+1)*dxmj/dx;
        % Compute correction
        if(dxmj/dx>=0.5)
            if(j<Nx-1)
                vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(vxf(i,j)-2*vxf(i,j+1)+vxf(i,j+2));
                vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(vxf(i+1,j)-2*vxf(i+1,j+1)+vxf(i+1,j+2));
            end
        else
            if(j>1)
                vxm13=vxm13+1/2*((dxmj/dx-0.5)^2)*(vxf(i,j-1)-2*vxf(i,j)+vxf(i,j+1));
                vxm24=vxm24+1/2*((dxmj/dx-0.5)^2)*(vxf(i+1,j-1)-2*vxf(i+1,j)+vxf(i+1,j+1));
            end
        end
        % Compute vx
        vxm(rk)=(1-dymi/dy)*vxm13+(dymi/dy)*vxm24;
        
        % Interpolate vy
        % Define i,j indexes for the upper left node
        j=fix((xcur-xvy(1))/dx)+1;
        i=fix((ycur-yvy(1))/dy)+1;
        if(j<1)
            j=1;
        elseif(j>Nx)
            j=Nx;
        end
        if(i<1)
            i=1;
        elseif(i>Ny-1)
            i=Ny-1;
        end
        % Compute distances
        dxmj=xcur-xvy(j);
        dymi=ycur-yvy(i);
        % Compute weights
        % Compute vy velocity for the left and right of the cell
        vym12=vyf(i,j)*(1-dymi/dy)+vyf(i+1,j)*dymi/dy;
        vym34=vyf(i,j+1)*(1-dymi/dy)+vyf(i+1,j+1)*dymi/dy;
        % Compute correction
        if(dymi/dy>=0.5)
            if(i<Ny-1)
                vym12=vym12+1/2*((dymi/dy-0.5)^2)*(vyf(i,j)-2*vyf(i+1,j)+vyf(i+2,j));
                vym34=vym34+1/2*((dymi/dy-0.5)^2)*(vyf(i,j+1)-2*vyf(i+1,j+1)+vyf(i+2,j+1));
            end      
        else
            if(i>1)
                vym12=vym12+1/2*((dymi/dy-0.5)^2)*(vyf(i-1,j)-2*vyf(i,j)+vyf(i+1,j));
                vym34=vym34+1/2*((dymi/dy-0.5)^2)*(vyf(i-1,j+1)-2*vyf(i,j+1)+vyf(i+1,j+1));
            end
        end
        % Compute vy
        vym(rk)=(1-dxmj/dx)*vym12+(dxmj/dx)*vym34;
        
        % Change coordinates to obtain B,C,D points
        if(rk==1 || rk==2)
            xcur=xA-dtm/2*vxm(rk);
            ycur=yA-dtm/2*vym(rk);
        elseif(rk==3)
            xcur=xA-dtm*vxm(rk);
            ycur=yA-dtm*vym(rk);
        end
    end

    % Compute effective velocity
    vxmeff=1/6*(vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4));
    vymeff=1/6*(vym(1)+2*vym(2)+2*vym(3)+vym(4));
    % Trace the node backward
    xcur=xA-dtm*vxmeff;
    ycur=yA-dtm*vymeff;
    % Interpolate nodal property
    % SIGMA'xx, P
    % Define i,j indexes for the upper left node
    j=fix((xcur-xp(1))/dx)+1;
    i=fix((ycur-yp(1))/dy)+1;
    if(j<1)
        j=1;
    elseif(j>Nx)
        j=Nx;
    end
    if(i<1)
        i=1;
    elseif(i>Ny)
        i=Ny;
    end
    % Compute distances
    dxmj=xcur-xp(j);
    dymi=ycur-yp(i);
    % Compute weights
    wtmij=(1-dxmj/dx)*(1-dymi/dy);
    wtmi1j=(1-dxmj/dx)*(dymi/dy);    
    wtmij1=(dxmj/dx)*(1-dymi/dy);
    wtmi1j1=(dxmj/dx)*(dymi/dy);
    % Compute nodal total pressure
    pf0(ii,jj)=pf(i,j)*wtmij+pf(i+1,j)*wtmi1j+...
            pf(i,j+1)*wtmij1+pf(i+1,j+1)*wtmi1j1;
end
end


% Add markers to empty areas
marknumold=marknum
mdis=1e30*ones(Nym,Nxm);
mnum=zeros(Nym,Nxm);
mtyp=zeros(Nym,Nxm);
mpor=zeros(Nym,Nxm);
xxm=dxm/2:dxm:xsize-dxm/2;
yym=dym/2:dym:ysize-dym/2;
for m=1:1:marknum
    
    % Check markers with the nearest nodes
    % Define i,j indexes for the upper left node
    j=fix((xm(m)-xxm(1))/dxm)+1;
    i=fix((ym(m)-yym(1))/dym)+1;
    if(j<1)
        j=1;
    elseif(j>Nxm-1)
        j=Nxm-1;
    end
    if(i<1)
        i=1;
    elseif(i>Nym-1)
        i=Nym-1;
    end
    
    % Check nodes
    % i,j Node
    % Compute distance
    dxmj=xm(m)-xxm(j);
    dymi=ym(m)-yym(i);
    dismij=(dxmj^2+dymi^2)^0.5;
    if(dismij<mdis(i,j))
        mdis(i,j)=dismij;
        mnum(i,j)=m;
        mtyp(i,j)=tm(m);
        mpor(i,j)=phim(m);
    end
    % i+1,j Node
    % Compute distance
    dxmj=xm(m)-xxm(j);
    dymi=ym(m)-yym(i+1);
    dismi1j=(dxmj^2+dymi^2)^0.5;
    if(dismi1j<mdis(i+1,j))
        mdis(i+1,j)=dismi1j;
        mnum(i+1,j)=m;
        mtyp(i+1,j)=tm(m);
        mpor(i+1,j)=phim(m);
    end
    % i,j+1 Node
    % Compute distance
    dxmj=xm(m)-xxm(j+1);
    dymi=ym(m)-yym(i);
    dismij1=(dxmj^2+dymi^2)^0.5;
    if(dismij1<mdis(i,j+1))
        mdis(i,j+1)=dismij1;
        mnum(i,j+1)=m;
        mtyp(i,j+1)=tm(m);
        mpor(i,j+1)=phim(m);
    end
    % i+1,j+1 Node
    % Compute distance
    dxmj=xm(m)-xxm(j+1);
    dymi=ym(m)-yym(i+1);
    dismi1j1=(dxmj^2+dymi^2)^0.5;
    if(dismi1j1<mdis(i+1,j+1))
        mdis(i+1,j+1)=dismi1j1;
        mnum(i+1,j+1)=m;
        mtyp(i+1,j+1)=tm(m);
        mpor(i+1,j+1)=phim(m);
    end
end

dii=5*Nxmc;
djj=5*Nymc;

for j=1:1:Nxm
    for i=1:1:Nym
        if(mnum(i,j)==0)
            % Serch surrounding nodes
            for jj=j-djj:1:j+djj
                for ii=i-dii:1:i+dii
                    if(ii>=1 && ii<=Nym && jj>=1 && jj<=Nxm && mnum(ii,jj)>0)
                        % Compute distance
                        m=mnum(ii,jj);
                        dxmj=xm(m)-xxm(j);
                        dymi=ym(m)-yym(i);
                        dismij=(dxmj^2+dymi^2)^0.5;
                        if(dismij<mdis(i,j))
                            mdis(i,j)=dismij;
                            mnum(i,j)=-m;
                            mtyp(i,j)=-tm(m);
                            mpor(i,j)=phim(m);
                        end
                    end
                end
            end
            % Add New marker
            if(mnum(i,j)<0)
                % Add marker number
                marknum=marknum+1;
                % Assign marker coordinates
                xm(marknum)=xxm(j);%+(rand-0.5)*dxm;
                ym(marknum)=yym(i);%+(rand-0.5)*dym;
                % Copy marker properties
                m=-mnum(i,j);
                tm(marknum)=tm(m); % Material type
                tkm(marknum)=tkm(m); % Temperature
                phim(marknum)=phim(m); % Porosity
                sxxm(marknum)=sxxm(m); % SIGMA'xx, Pa
                sxym(marknum)=sxym(m); % SIGMAxy, Pa
                etavpm(marknum)=etavpm(m); % Visco-plastic viscosity, Pa
            end
        end
    end
end           
marknumnew=marknum



% Update timesum
timesum=timesum+dtm;

% Translate vx,vy and qxD,qyD into polar coordinates for visualization
vrp=sqrt((vxp+vxp((Ny+1)*(Nx+1)/2)).^2+(vyp+vyp((Ny+1)*(Nx+1)/2)).^2);
qrD=sqrt((qxD+qxD((Ny+1)*(Nx+1)/2)).^2+(qyD+qyD((Ny+1)*(Nx+1)/2)).^2);

if(fix(timestep/savematstep)*savematstep==timestep)

nname1 =[nnamefile,'a_'];
nname2 =[nnamefile,'b_'];
nname3 =[nnamefile,'c_'];
DIVV=zeros(Ny1,Nx1);
gmag=zeros(Ny1,Nx1);
for i=2:1:Ny
  for j=2:1:Nx
    % EXX, SXX, DSXX
    DIVV(i,j)=(vx(i,j)-vx(i,j-1))/dx+(vy(i,j)-vy(i-1,j))/dy;
    gmag(i,j)=(((gx(i,j)+gx(i,j-1))/2)^2+((gy(i,j)+gy(i-1,j))/2)^2)^0.5;
  end
end

% 8) Visialize RHO() P(), vx(), vy() 
figure(1);clf;colormap('Jet')
subplot(3,3,1)
pcolor(xp,yp,RHO);
colorbar; 
shading interp;
axis ij image;
title(['RHO    step = ',num2str(timestep),'   time = ',num2str((timesum+dt)/(365.25*24*3600*1e+6)),' Myr'])

subplot(3,3,2)
pcolor(x,y,log10(ETA));
colorbar; % 
shading interp; 
axis ij image; 
title(['log10 viscosity, Pa*s  step = ',num2str(timestep),'   dtelastic = ',num2str(dt/(365.25*24*3600)),' yr','   dt = ',num2str(dtm/(365.25*24*3600)),' yr'])

subplot(3,3,3)
pcolor(xvx,yvx,vx);
colorbar; 
shading interp;
axis ij image;
title('vxs, m/s')

subplot(3,3,4)
pcolor(xvy,yvy,vy); 
colorbar;
shading interp;
axis ij image;
title('vys, m/s')

subplot(3,3,5)
pcolor(xp/1000,yp/1000,log10(PHI./(1-PHI)))
shading interp;
axis ij image;
colorbar
title('logPHI/(1-PHI)')


subplot(3,3,6)
pcolor(xp(2:Nx),yp(2:Ny),log10(EII(2:Ny,2:Nx)));
colorbar;
shading interp;
axis ij image;
title('log10 EII, 1/s')

subplot(3,3,7)
pcolor(xp(2:Nx),yp(2:Ny),SII(2:Ny,2:Nx));
colorbar;
shading interp;
axis ij image;
title('SII, Pa')

subplot(3,3,8)
pcolor(xp(2:Nx),yp(2:Ny),gmag(2:Ny,2:Nx))
shading interp;
axis ij image;
colorbar
title('gmag, m/s^2')

subplot(3,3,9)
pcolor(xp/1000,yp/1000,tk2)
shading interp;
axis ij image;
colorbar
title('T, K')


    namejpg    =  [nname1,num2str(timestep)];
    print ('-djpeg', '-r300',namejpg);

    
% 8) Visialize RHO() P(), vx(), vy() 
figure(2);clf;colormap('Jet')
subplot(4,4,1)
pcolor(x/1000,y/1000,log10(ETA));% caxis([17 21])
shading flat;
axis ij image;
colorbar
title(['log10 visc,Pa*s step=',num2str(timestep),'   dtelastic = ',num2str(dt/(365.25*24*3600)),' yr','   dt = ',num2str(dtm/(365.25*24*3600)),' yr'])

hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'w')

subplot(4,4,2)
pcolor(xp/1000,yp/1000,pr)
shading interp;
axis ij image;
colorbar
title(['Ptotal,Pa time=',num2str((timesum+dt)/(365.25*24*3600*1e+6)),' Myr'])
hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(4,4,3)
pcolor(xp/1000,yp/1000,vxp)
shading interp;
axis ij image;
colorbar
title('vx, m/s')
hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(4,4,4)
pcolor(xp/1000,yp/1000,vyp)
shading interp;
axis ij image;
colorbar
title('vy, m/s')
hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(4,4,5)
pcolor(xp/1000,yp/1000,HS)
shading interp;
axis ij image;
colorbar
title('HS, W/m^3')

subplot(4,4,6)
pcolor(xp/1000,yp/1000,HA)
shading interp;
axis ij image;
colorbar
title('HA, W/m^3')

subplot(4,4,7)
pcolor(xp/1000,yp/1000,pr-pf)
shading interp;
axis ij image;
colorbar
title('Ptotal-Pfluid, kg/m^3')

subplot(4,4,8)
pcolor(xvx/1000,yvx/1000,log10(KX))
shading interp;
axis ij image;
colorbar
title('logK, W/m/K')

subplot(4,4,9)
pcolor(xp/1000,yp/1000,tk2)
shading interp;
axis ij image;
colorbar
title('T, K')

subplot(4,4,10)
pcolor(xp/1000,yp/1000,log10(PHI./(1-PHI)))
shading interp;
axis ij image;
colorbar
title('logPHI/(1-PHI)')

subplot(4,4,11)
pcolor(xp/1000,yp/1000,pf)
shading interp;
axis ij image;
colorbar
title('Pfluid, Pa')

subplot(4,4,12)
pcolor(xvx/1000,yvx/1000,qxD)
shading interp;
axis ij image;
colorbar
title('qxD, m/s')

subplot(4,4,13)
pcolor(xvy/1000,yvy/1000,qyD)
shading interp;
axis ij image;
colorbar
title('qyD, m/s')

subplot(4,4,14)
pcolor(xvx/1000,yvx/1000,log10(RX))
shading interp;
axis ij image;
colorbar
title('logRX')

subplot(4,4,15)
pcolor(xp/1000,yp/1000,log10(ETAPHI))
shading interp;
axis ij image;
colorbar
title('logETAphi, Pa*s')

subplot(4,4,16)
pcolor(xp,yp,RHO);
colorbar; 
shading interp;
axis ij image;
title(['RHO, kg/m^3'])


    namejpg    =  [nname2,num2str(timestep)];
    print ('-djpeg', '-r300',namejpg);



    namemat    =  [nname,num2str(timestep)];
    save(namemat);
    fdata=fopen('file.txt','wt');
    fprintf(fdata,'%d',timestep);
    fclose(fdata);





end

if timesum > 15*3600*24*365.25*1000000
    break
end
end

range = [1:nsteps];
s1 = figure;
s1.Position =[10 300 1800 800];
p1 = subplot(1,2,1);
axis('square');
plot(range, maxT_m, '-*');
hold on
plot(range, maxT_j, '-o');
hold off
xlabel('timestep');
xlim([1 nsteps]);
ylabel('max(T) [K]');
ylim([170 200]);
title('maximum temperature');
legend('MATLAB planetary', 'Julia HydrologyPlanetesimals');
p2 = subplot(1,2,2);
axis('square');
plot(range, dRHO_m, '-*');
hold on
plot(range, dRHO_j, '-o');
hold off
xlabel('timestep');
xlim([1 nsteps]);
ylabel('\Delta RHO [kg m^{-3}]');
title('density contrast center to crust');
legend('MATLAB planetary', 'Julia HydrologyPlanetesimals');
orient(s1,'landscape');
sgtitle('Comparing maximum temperature, density constrast');
% saveas(s1, 'comparison_1.pdf');
set(gcf,'color','w'); export_fig('comparison_1', '-pdf','-m2');

s2 = figure;
s2.Position =[10 300 1800 800];
p1 = subplot(1,2,1);
axis('square');
plot(range, meanT_10_m, '-*');
hold on
plot(range, meanT_10_j, '-o');
hold off
xlabel('timestep');
xlim([1 nsteps]);
ylabel('mean(T) at R=10km [K]');
ylim([170 200]);
title('mean temperature at radius R=10km');
legend('MATLAB planetary', 'Julia HydrologyPlanetesimals');
p2 = subplot(1,2,2);
axis('square');
plot(range, meanT_40_m, '-*');
hold on
plot(range, meanT_40_j, '-o');
hold off
xlabel('timestep');
xlim([1 nsteps]);
ylim([170 200]);    
ylabel('mean(T) at R=40km [K]');
title('mean temperature at radius R=40km');
legend('MATLAB planetary', 'Julia HydrologyPlanetesimals');
sgtitle('Comparing mean temperature at radii R=10km, R=40km');
orient(s2,'landscape');
% saveas(s2, 'comparison_2.pdf');
set(gcf,'color','w'); export_fig('comparison_2', '-pdf','-m2');

s3 = figure;
s3.Position =[10 300 1800 800];
p1 = subplot(1,2,1);
axis('square');
plot(range, meanRHO_10_m, '-*');
hold on
plot(range, meanRHO_10_j, '-o');
hold off
xlabel('timestep');
xlim([1 nsteps]);
ylabel('mean(RHO) at R=10km [kg m^{-3}]');
ylim([2750 3000]);
title('mean density at radius R=10km');
legend('MATLAB planetary', 'Julia HydrologyPlanetesimals');
p2 = subplot(1,2,2);
axis('square');
plot(range, meanRHO_40_m, '-*');
hold on
plot(range, meanRHO_40_j, '-o');
hold off
xlabel('timestep');
xlim([1 nsteps]);
ylabel('mean(RHO) at R=40km [kg m^{-3}]');
ylim([2750 3000]);
title('mean density at radius R=40km');
legend('MATLAB planetary', 'Julia HydrologyPlanetesimals');
sgtitle('Comparing mean density at radii R=10km, R=40km');
orient(s3,'landscape');
% saveas(s3, 'comparison_3.pdf');
set(gcf,'color','w'); export_fig('comparison_3', '-pdf','-m2');

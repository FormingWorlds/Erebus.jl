% Viscous hydro-thermomechanical (HTM) code
% Solving momentum, mass and energy conservation eqs.
% for coupled fluid-solid system with melting
% in primitive variable formulation
% with variable viscosity and thermal conductivity
% using FD with staggered grid

% Clearing memory and figures
clear all; clf

% Define Numerical model
xsize=700000; % Horizontal model size, m
ysize=400000; % Vertical model size, m
Nx=71; % Horizontal grid resolution
Ny=41; % Vertical grid resolution
Nx1=Nx+1;
Ny1=Ny+1;
dx=xsize/(Nx-1); % Horizontal grid step, m
dy=ysize/(Ny-1); % Vertical grid step, m

% Define Gravity
gx=0; % Horizontal gravity acceleration, m/s^2
gy=10; % Vertical gravity acceleration, m/s^2

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
ETA=zeros(Ny,Nx); % Viscosity, Pa*s
EXY=zeros(Ny,Nx); % EPSILONxy, 1/s
SXY=zeros(Ny,Nx); % SIGMAxy, 1/s
% Vx-Nodes
RHOX=zeros(Ny1,Nx1); % Total density, kg/m^3
RHOFX=zeros(Ny1,Nx1); % Fluid density, kg/m^3
KX=zeros(Ny1,Nx1); % Thermal conductivity, W/m/K
PHIX=zeros(Ny1,Nx1); % Porosity
vx=zeros(Ny1,Nx1); % Solid vx-velocity m/s
vxf=zeros(Ny1,Nx1); % Fluid vx-velocity m/s
RX=zeros(Ny1,Nx1); % ETAfluid/Kphi ratio , m^2
qxD=zeros(Ny1,Nx1); % qx-Darcy flux m/s
% Vy-Nodes
RHOY=zeros(Ny1,Nx1); % Total density, kg/m^3
RHOFY=zeros(Ny1,Nx1); % Fluid density, kg/m^3
KY=zeros(Ny1,Nx1); % Thermal conductivity, W/m/K
PHIY=zeros(Ny1,Nx1); % Porosity
vy=zeros(Ny1,Nx1); % Solid vy-velocity m/s
vyf=zeros(Ny1,Nx1); % Fluid vy-velocity m/s
RY=zeros(Ny1,Nx1); % ETAfluid/Kphi ratio , m^2
qyD=zeros(Ny1,Nx1); % qy-Darcy flux m/s
% P-nodes
RHO=zeros(Ny1,Nx1); % Density, kg/m^3
RHOCP=zeros(Ny1,Nx1); % Volumetric heat capacity, J/m^3/K
ALPHA=zeros(Ny1,Nx1); % Thermal expansion, J/m^3/K
HR=zeros(Ny1,Nx1); % Radioactive heating, W/m^3
HA=zeros(Ny1,Nx1); % Adiabatic heating, W/m^3
HS=zeros(Ny1,Nx1); % Shear heating, W/m^3
ETAP=zeros(Ny1,Nx1); % Viscosity, Pa*s
EXX=zeros(Ny,Nx); % EPSILONxx, 1/s
SXX=zeros(Ny,Nx); % SIGMAxx, 1/s
tk1=zeros(Ny1,Nx1); % Old temperature, K
tk2=zeros(Ny1,Nx1); % New temperature, K
vxp=zeros(Ny1,Nx1); % Solid Vx in pressure nodes, m/s
vyp=zeros(Ny1,Nx1); % Solid Vy in pressure nodes, m/s
vxpf=zeros(Ny1,Nx1); % Fluid Vx in pressure nodes, m/s
vypf=zeros(Ny1,Nx1); % Fluid Vy in pressure nodes, m/s
pr=zeros(Ny1,Nx1); % Total Pressure, Pa
pf=zeros(Ny1,Nx1); % Fluid Pressure, Pa
ETAPHI=zeros(Ny1,Nx1); % Bulk Viscosity, Pa*s
PHI=zeros(Ny1,Nx1); % porosity
APHI=zeros(Ny1,Nx1); % Dln((1-PHI)/PHI)/Dt

% Define markers
Nxmc=4;
Nymc=4;
Nxm=(Nx-1)*Nxmc; % Marker grid resolution in horizontal direction
Nym=(Ny-1)*Nymc; % Marker grid resolution in vertical direction
dxm=xsize/Nxm; % Marker grid step in horizontal direction,m
dym=ysize/Nym; % Marker grid step in vertical direction,m
marknum=Nxm*Nym; % Number of markers
xm=zeros(1,marknum); % Horizontal coordinates, m
ym=zeros(1,marknum); % Vertical coordinates, m
tm=zeros(1,marknum); % Material type
tkm=zeros(1,marknum); % Marker temperature, K
phim=zeros(1,marknum); % Marker porosity
phinewm=zeros(1,marknum); % New marker porosity
XWsolidm0=zeros(1,marknum); % Old wet silicate content in solid
XWsolidm=zeros(1,marknum); % New VWsolid content in solid
pfm0=zeros(1,marknum); % Old fluid pressure, Pa

% Define properties of materials: 
%             mantle  crust  lith.   air
rhosolidm   = [3300   2600   3300    1     ]; % Solid Density, kg/m^3
etasolidm   = [1e+19  1e+18  1e+23   1e+17 ]; % Viscosity, Pa s
rhocpsolidm = [3.3e+6 3.2e+6 3.3e+6  3.3e+6]; % Solid volumetric heat capacity, kg/m^3
alphasolidm = [3e-5   3e-5   3e-5    0     ]; % Solid Thermal expansion, 1/K
ksolidm     = [3      2      3       3000  ]; % Thermal conductivity, W/m/K
hrsolidm    = [2e-8   3e-8   2e-8    0     ]; % Radiogenic heat production, W/m^3
kphim0      = [1e-14  1e-14  1e-14   1e-14 ]; % Standard permeability, m^2

% Fluid properties
etafluid=0.001; % Fluid viscosity
rhocpfluid=3e+6; % Fluid/Melt volumetric isobaric heat capacity, J/K/m^3
alphafluid=5e-5; % Thermal expansion of the fluid
hrfluid=2e-7; % Fluid/Melt heat production, W/m^3
kfluid=5; % Fluid/Melt conductivity

phim0=0.01; % standard porosity
phimin=1e-4; % Min porosity
phimax=1-phimin; % Max porosity

% Thermodynamic model
% Dehydration reaction Wsilicate=Dsilicate+H2O
MH2O=0.018; % Molar mass of water, kg/mol
MD=0.120; % Molar mass of dry silicate, kg/mol 
rhoDsolid=3300; % Density of dry silicate, kg/m^3
rhoWsolid=2600; % Density of wet silicate, kg/m^3
rhoH2Ofluid=1000; % Density of liquid water, kg/m^3
VDsolid=MD/rhoDsolid; % Molar volume of dry solid, m^3/mol
VWsolid=(MD+MH2O)/rhoWsolid; % Molar volume of wet solid, m^3/mol
VH2Ofluid=MH2O/rhoH2Ofluid; % Molar volume of liquid water, m^3/mol
dHWD=40000; % Enthalpy change for dehydration of the wet silicate, J/mol, 
dSWD=60; % Entropy change for dehydration of the wet silicate, J/K/mol, 
dVWD=VDsolid+VH2Ofluid-VWsolid; % Volume change for dehydration of the wet silicate, J/Pa/mol, 
dtreaction=1e+10; % Timescale to complete dehydration reaction, s

% Define marker coordinates, temperature and material type
rbending=750000; % initial slab bending radius, m
tplate=50000; % plate thickness, m
tcrust=20000; % crustal thickness, m
xislab=100000; % slab initial horizontal position, m
lbending=200000; % initial bending length, m
m=1; % Marker counter
for jm=1:1:Nxm
    for im=1:1:Nym
        % Define marker coordinates
        xm(m)=dxm/2+(jm-1)*dxm+(rand-0.5)*dxm;
        ym(m)=dym/2+(im-1)*dym+(rand-0.5)*dym;
        % Mantle
        tm(m)=1; % Material type
        tkm(m)=1500; % Temperature
        phim(m)=phimin;
        XWsolidm0(m)=0.01;
        % Lithosphere window temperature
        if(ym(m)<110000 && ym(m)>50000 && xm(m)<=xbending)
            tkm(m)=273+(1500-273)*(ym(m)-50000)/60000; % Temperature
        end
        % Subducting plate temperature
        if(ym(m)>50000 && ym(m)<50000+tcrust+tplate && xm(m)>xbending)
            tkm(m)=273+(1500-273)*(ym(m)-50000)/(tcrust+tplate); % Temperature
        end
        % Bent slab
        xbending=xislab+lbending;
        ybending=50000+rbending;
        rm=((xm(m)-xbending)^2+(ym(m)-ybending)^2)^0.5;
        abending=pi*lbending/(pi*rbending);
        % Bent slab
        if((xm(m)<=xbending && rm>rbending-tcrust-tplate && rm<rbending-tcrust && xm(m)>xbending-(ybending-ym(m))*tan(abending)) ||...  
            (xm(m)>xbending && ym(m)>50000+tcrust && ym(m)<50000+tcrust+tplate))
            tm(m)=3; % Material type
            if(xm(m)<=xbending)
            tkm(m)=1500+(273-1500)*(rm-(rbending-tcrust-tplate))/(tcrust+tplate); % Temperature
            end
        end
       % Crust
        if(ym(m)>50000 && ym(m)<=50000+tcrust)
            tm(m)=2; % Material type
            XWsolidm0(m)=0.5;
        end
        % Sticky air (to have internal free surface)
        if(ym(m)<=50000)
            tm(m)=4; % Material type
            tkm(m)=273; % Temperature
        end
        % Update marker counter
        m=m+1;
    end
end
XWsolidm=XWsolidm0;
phinewm=phim;


% Introducing scaled pressure
pscale=1e+21/dx;

% Define global matrixes 
% Hydro-Mechanical solution: L(), R()
N=Nx1*Ny1*6; % Global number of unknowns
L=sparse(N,N); % Matrix of coefficients (left part)
R=zeros(N,1); % Vector of right parts
% Thermal solution: LT(), RT()
N=Nx1*Ny1; % Global number of unknowns
LT=sparse(N,N); % Matrix of coefficients (left part)
RT=zeros(N,1); % Vector of right parts

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
% Thermal boundary conditions
tktop=273;
tkbottom=1500;

% Timestepping
dt=1e+10; % initial timestep, s
timesum=0; % Time counter, s
dtkoef=1.2; % timestep increment
dxymax=0.1; % max marker movement per time step, grid steps
vpratio=1/3; % Weight of averaged velocity for moving markers
DTmax=20; % max temperature change per time step, K
dsubgridt=1; % subgrid diffusion parameter
dphimax=100.1; % max porosity ratio change per time step
pferrmax=1e+5; % Error limit to exit eterations, Pa
pfkoef=0.5; % Koeeficient for using pressure from the previous iteration
titermax=10000; % Maximal number of global iterations
timestep=1;
for timestep=timestep:1:1000
    
% Interpolate properties from markers to nodes
% Basic nodes
ETASUM=zeros(Ny,Nx);
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
ETAPSUM=zeros(Ny1,Nx1);
RHOSUM=zeros(Ny1,Nx1);
RHOCPSUM=zeros(Ny1,Nx1);
ALPHASUM=zeros(Ny1,Nx1);
HRSUM=zeros(Ny1,Nx1);
TKSUM=zeros(Ny1,Nx1);
PHISUM=zeros(Ny1,Nx1);
WTPSUM=zeros(Ny1,Nx1);

for m=1:1:marknum
       
    % Compute marker parameters
    if(tm(m)<4)
        % Rocks
        % Compute density of solid and fluid
        XDsolidm0=1-XWsolidm0(m);
        rhosolidm0=(MD+MH2O*XWsolidm0(m))/...
                  (VWsolid*XWsolidm0(m)+VDsolid*XDsolidm0)*...
                  exp(-alphasolidm(tm(m))*(tkm(m)-273));
        rhofluidm0=rhoH2Ofluid*...
                  exp(-alphafluid*(tkm(m)-273));
        % Compute properties
        kphim=kphim0(tm(m))*(phim(m)/phim0)^3/((1-phim(m))/(1-phim0))^2; %Permeability
        rhototalm=rhosolidm0*(1-phim(m))+rhofluidm0*phim(m);
        rhocptotalm=rhocpsolidm(tm(m))*(1-phim(m))+rhocpfluid*phim(m);
        etatotalm=etasolidm(tm(m))*exp(-28*phim(m));
        hrtotalm=hrsolidm(tm(m))*(1-phim(m))+hrfluid*phim(m);
        ktotalm=(ksolidm(tm(m))*kfluid/2+((ksolidm(tm(m))*(3*phim(m)-2)+...
            kfluid*(1-3*phim(m)))^2)/16)^0.5-(ksolidm(tm(m))*(3*phim(m)-2)+...
            kfluid*(1-3*phim(m)))/4;
    else
        % Sticky air
        kphim=kphim0(tm(m))*(phim(m)/phim0)^3/((1-phim(m))/(1-phim0))^2; %Permeability
        rhofluidm0=rhosolidm(tm(m));
        rhototalm=rhosolidm(tm(m));
        rhocptotalm=rhocpsolidm(tm(m));
        etatotalm=etasolidm(tm(m));
        hrtotalm=hrsolidm(tm(m));
        ktotalm=ksolidm(tm(m));
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
    ETASUM(i,j)=ETASUM(i,j)+etatotalm*wtmij;
    WTSUM(i,j)=WTSUM(i,j)+wtmij;
    % i+1,j Node
    ETASUM(i+1,j)=ETASUM(i+1,j)+etatotalm*wtmi1j;
    WTSUM(i+1,j)=WTSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    ETASUM(i,j+1)=ETASUM(i,j+1)+etatotalm*wtmij1;
    WTSUM(i,j+1)=WTSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    ETASUM(i+1,j+1)=ETASUM(i+1,j+1)+etatotalm*wtmi1j1;
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
    RHOFXSUM(i,j)=RHOFXSUM(i,j)+rhofluidm0*wtmij;
    KXSUM(i,j)=KXSUM(i,j)+ktotalm*wtmij;
    PHIXSUM(i,j)=PHIXSUM(i,j)+phim(m)*wtmij;
    RXSUM(i,j)=RXSUM(i,j)+etafluid/kphim*wtmij;
    WTXSUM(i,j)=WTXSUM(i,j)+wtmij;
    % i+1,j Node
    RHOXSUM(i+1,j)=RHOXSUM(i+1,j)+rhototalm*wtmi1j;
    RHOFXSUM(i+1,j)=RHOFXSUM(i+1,j)+rhofluidm0*wtmi1j;
    KXSUM(i+1,j)=KXSUM(i+1,j)+ktotalm*wtmi1j;
    PHIXSUM(i+1,j)=PHIXSUM(i+1,j)+phim(m)*wtmi1j;
    RXSUM(i+1,j)=RXSUM(i+1,j)+etafluid/kphim*wtmi1j;
    WTXSUM(i+1,j)=WTXSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    RHOXSUM(i,j+1)=RHOXSUM(i,j+1)+rhototalm*wtmij1;
    RHOFXSUM(i,j+1)=RHOFXSUM(i,j+1)+rhofluidm0*wtmij1;
    KXSUM(i,j+1)=KXSUM(i,j+1)+ktotalm*wtmij1;
    PHIXSUM(i,j+1)=PHIXSUM(i,j+1)+phim(m)*wtmij1;
    RXSUM(i,j+1)=RXSUM(i,j+1)+etafluid/kphim*wtmij1;
    WTXSUM(i,j+1)=WTXSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    RHOXSUM(i+1,j+1)=RHOXSUM(i+1,j+1)+rhototalm*wtmi1j1;
    RHOFXSUM(i+1,j+1)=RHOFXSUM(i+1,j+1)+rhofluidm0*wtmi1j1;
    KXSUM(i+1,j+1)=KXSUM(i+1,j+1)+ktotalm*wtmi1j1;
    RXSUM(i+1,j+1)=RXSUM(i+1,j+1)+etafluid/kphim*wtmi1j1;
    PHIXSUM(i+1,j+1)=PHIXSUM(i+1,j+1)+phim(m)*wtmi1j1;
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
    RHOFYSUM(i,j)=RHOFYSUM(i,j)+rhofluidm0*wtmij;
    KYSUM(i,j)=KYSUM(i,j)+ktotalm*wtmij;
    PHIYSUM(i,j)=PHIYSUM(i,j)+phim(m)*wtmij;
    RYSUM(i,j)=RYSUM(i,j)+etafluid/kphim*wtmij;
    WTYSUM(i,j)=WTYSUM(i,j)+wtmij;
    % i+1,j Node
    RHOYSUM(i+1,j)=RHOYSUM(i+1,j)+rhototalm*wtmi1j;
    RHOFYSUM(i+1,j)=RHOFYSUM(i+1,j)+rhofluidm0*wtmi1j;
    KYSUM(i+1,j)=KYSUM(i+1,j)+ktotalm*wtmi1j;
    PHIYSUM(i+1,j)=PHIYSUM(i+1,j)+phim(m)*wtmi1j;
    RYSUM(i+1,j)=RYSUM(i+1,j)+etafluid/kphim*wtmi1j;
    WTYSUM(i+1,j)=WTYSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    RHOYSUM(i,j+1)=RHOYSUM(i,j+1)+rhototalm*wtmij1;
    RHOFYSUM(i,j+1)=RHOFYSUM(i,j+1)+rhofluidm0*wtmij1;
    KYSUM(i,j+1)=KYSUM(i,j+1)+ktotalm*wtmij1;
    PHIYSUM(i,j+1)=PHIYSUM(i,j+1)+phim(m)*wtmij1;
    RYSUM(i,j+1)=RYSUM(i,j+1)+etafluid/kphim*wtmij1;
    WTYSUM(i,j+1)=WTYSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    RHOYSUM(i+1,j+1)=RHOYSUM(i+1,j+1)+rhototalm*wtmi1j1;
    RHOFYSUM(i+1,j+1)=RHOFYSUM(i+1,j+1)+rhofluidm0*wtmi1j1;
    KYSUM(i+1,j+1)=KYSUM(i+1,j+1)+ktotalm*wtmi1j1;
    PHIYSUM(i+1,j+1)=PHIYSUM(i+1,j+1)+phim(m)*wtmi1j1;
    RYSUM(i+1,j+1)=RYSUM(i+1,j+1)+etafluid/kphim*wtmi1j1;
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
    ETAPSUM(i,j)=ETAPSUM(i,j)+etatotalm*wtmij;
    RHOSUM(i,j)=RHOSUM(i,j)+rhototalm*wtmij;
    RHOCPSUM(i,j)=RHOCPSUM(i,j)+rhocptotalm*wtmij;
    ALPHASUM(i,j)=ALPHASUM(i,j)+alphasolidm(tm(m))*wtmij;
    HRSUM(i,j)=HRSUM(i,j)+hrtotalm*wtmij;
    TKSUM(i,j)=TKSUM(i,j)+tkm(m)*rhocptotalm*wtmij;
    PHISUM(i,j)=PHISUM(i,j)+phim(m)*wtmij;
    WTPSUM(i,j)=WTPSUM(i,j)+wtmij;
    % i+1,j Node
    ETAPSUM(i+1,j)=ETAPSUM(i+1,j)+etatotalm*wtmi1j;
    RHOSUM(i+1,j)=RHOSUM(i+1,j)+rhototalm*wtmi1j;
    RHOCPSUM(i+1,j)=RHOCPSUM(i+1,j)+rhocptotalm*wtmi1j;
    ALPHASUM(i+1,j)=ALPHASUM(i+1,j)+alphasolidm(tm(m))*wtmi1j;
    HRSUM(i+1,j)=HRSUM(i+1,j)+hrtotalm*wtmi1j;
    TKSUM(i+1,j)=TKSUM(i+1,j)+tkm(m)*rhocptotalm*wtmi1j;
    PHISUM(i+1,j)=PHISUM(i+1,j)+phim(m)*wtmi1j;
    WTPSUM(i+1,j)=WTPSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    ETAPSUM(i,j+1)=ETAPSUM(i,j+1)+etatotalm*wtmij1;
    RHOSUM(i,j+1)=RHOSUM(i,j+1)+rhototalm*wtmij1;
    RHOCPSUM(i,j+1)=RHOCPSUM(i,j+1)+rhocptotalm*wtmij1;
    ALPHASUM(i,j+1)=ALPHASUM(i,j+1)+alphasolidm(tm(m))*wtmij1;
    HRSUM(i,j+1)=HRSUM(i,j+1)+hrtotalm*wtmij1;
    TKSUM(i,j+1)=TKSUM(i,j+1)+tkm(m)*rhocptotalm*wtmij1;
    PHISUM(i,j+1)=PHISUM(i,j+1)+phim(m)*wtmij1;
    WTPSUM(i,j+1)=WTPSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    ETAPSUM(i+1,j+1)=ETAPSUM(i+1,j+1)+etatotalm*wtmi1j1;
    RHOSUM(i+1,j+1)=RHOSUM(i+1,j+1)+rhototalm*wtmi1j1;
    RHOCPSUM(i+1,j+1)=RHOCPSUM(i+1,j+1)+rhocptotalm*wtmi1j1;
    ALPHASUM(i+1,j+1)=ALPHASUM(i+1,j+1)+alphasolidm(tm(m))*wtmi1j1;
    HRSUM(i+1,j+1)=HRSUM(i+1,j+1)+hrtotalm*wtmi1j1;
    TKSUM(i+1,j+1)=TKSUM(i+1,j+1)+tkm(m)*rhocptotalm*wtmi1j1;
    PHISUM(i+1,j+1)=PHISUM(i+1,j+1)+phim(m)*wtmi1j1;
    WTPSUM(i+1,j+1)=WTPSUM(i+1,j+1)+wtmi1j1;
end
% Compute physical properties
% Basic nodes
for j=1:1:Nx
    for i=1:1:Ny
        if(WTSUM(i,j)>0)
            ETA(i,j)=ETASUM(i,j)/WTSUM(i,j);
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
            PHI(i,j)=PHISUM(i,j)/WTPSUM(i,j);
            ETAP(i,j)=ETAPSUM(i,j)/WTPSUM(i,j);
            ETAPHI(i,j)=ETAP(i,j)/PHI(i,j);
            RHO(i,j)=RHOSUM(i,j)/WTPSUM(i,j);
            RHOCP(i,j)=RHOCPSUM(i,j)/WTPSUM(i,j);
            ALPHA(i,j)=ALPHASUM(i,j)/WTPSUM(i,j);
            HR(i,j)=HRSUM(i,j)/WTPSUM(i,j);
            tk1(i,j)=TKSUM(i,j)/RHOCPSUM(i,j);
        end
    end
end
% Applying thermal boundary conditions for interpolated temperature
% Upper boundary 
tk1(1,2:Nx)=2*tktop-tk1(2,2:Nx); % Constant temperature
% Lower boundary 
tk1(Ny1,2:Nx)=2*tkbottom-tk1(Ny,2:Nx); % Constant temperature
% Left boundary
tk1(:,1)=tk1(:,2); % Insulating boundary
% Right boundary
tk1(:,Nx1)=tk1(:,Nx); % Insulating boundary


% HTM iteration
for titer=1:1:titermax
% Compute mass transfer rate
DMPSUM=zeros(Ny1,Nx1);
DHPSUM=zeros(Ny1,Nx1);
WTPSUM=zeros(Ny1,Nx1);
for m=1:1:marknum
if(tm(m)<4)
    % Rocks
    
    % Interpolate fluid pressure and temperature to the marker
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
    % Interpolate nodal fluid pressure and temperature to the marker
    pfnm=pf(i,j)*wtmij+pf(i+1,j)*wtmi1j+...
            pf(i,j+1)*wtmij1+pf(i+1,j+1)*wtmi1j1;
    tknm=tk2(i,j)*wtmij+tk2(i+1,j)*wtmi1j+...
            tk2(i,j+1)*wtmij1+tk2(i+1,j+1)*wtmi1j1;
    pfnm=max(pfnm,0);
    % Use pressure from the previous iteration
    if(titer>2)
        pfnm=pfnm*(1-pfkoef)+pfm0(m)*pfkoef;
    end
    pfm0(m)=pfnm;

    % Thermodynamic computations
    % Compute bulk composition of the solid+fluid system
    XDsolidm0=1-XWsolidm0(m);
    Xfluid0=phim(m)*(XWsolidm0(m)*VWsolid+XDsolidm0*VDsolid)/...
       ((1-phim(m))*VH2Ofluid+...
            phim(m)*(XWsolidm0(m)*VWsolid+XDsolidm0*VDsolid));
    Xsolid0=1-Xfluid0;
    XH2Ototal=(XWsolidm0(m)*Xsolid0+Xfluid0)/(1+XWsolidm0(m)*Xsolid0);
    XDtotal=1-XH2Ototal;
    % Compute old density of the solid and fluid
    rhosolid0=(MD+MH2O*XWsolidm0(m))/...
              (VWsolid*XWsolidm0(m)+VDsolid*XDsolidm0);
    rhofluid0=rhoH2Ofluid;
    % Compute old relative ehthalpy of the system
    Htotal0=-Xsolid0*XWsolidm0(m)*dHWD/(MD+MH2O);
    % Compute old dG for dehydration reaction: Wsilicate=Dsilicate+H2O
    dGWD0=dHWD-tknm*dSWD+dVWD*pfnm+8.314*tknm*log(XDsolidm0/XWsolidm0(m));
    % Compute incomplete reaction for too short timestep
    dGWD=0;
    if(dt<dtreaction)
        dGWD=dGWD0*(1-dt/dtreaction);
    end

    % Compute equilibrium compositions and fluid fraction
    % Dehydration reaction: Wsilicate=Dsilicate+H2O
    KWD=exp(-(dHWD-tknm*dSWD+dVWD*pfnm-dGWD)/8.314/tknm);
    % Solid composition
    XWsolidm1=1/(KWD+1);
    XDsolidm1=1-XWsolidm1;
    % Fluid, Solid molar fraction
    Xsolid1=XDtotal/(1-XDtotal*XWsolidm1);
    Xfluid1=1-Xsolid1;
    % Process fluid-bearing rocks only 
    if(Xfluid1>0 && Xfluid1<1)
        % Compute equilibrium Porosity
        phinew1=Xfluid1*VH2Ofluid/(Xfluid1*VH2Ofluid+...
                Xsolid1*(XWsolidm1*VWsolid+XDsolidm1*VDsolid));
        % Compute equilibrium density of the solid and fluid
        rhosolid=(MD+MH2O*XWsolidm1)/...
                  (VWsolid*XWsolidm1+VDsolid*XDsolidm1);
        rhofluid=rhoH2Ofluid;
        % Compute equilibrium relative ehthalpy of the system
        Htotal1=-Xsolid1*XWsolidm1*dHWD/(MD+MH2O);
        % Compute ehthalpy change
        dHtotal=Htotal1-Htotal0;
        % Compute old/equilibrium volume ratio
        RV=(rhosolid*(1-phinew1)+rhofluid*phinew1)/...
           (rhosolid0*(1-phim(m))+rhofluid0*phim(m));
        % Compute mass transfer rate
        Gmass=(rhosolid0*RV*(1-phim(m))-rhosolid*(1-phinew1))/dt;
%         % Compute density of the transferred component C in the solid
%         rhoCsolid=(rhosolid0*RV*(1-phim(m))-rhosolid*(1-phinew1))/...
%                   (RV*(1-phim(m))-(1-phinew1));
%         % Compute density of the transferred component C in the fluid
%         rhoCfluid=(rhosolid0*RV*(1-phim(m))-rhosolid*(1-phinew1))/...
%                   (phinew1-RV*phim(m));
        % Compute mass transfer term
%         DMm=Gmass*(rhoCsolid-rhoCfluid)/rhoCsolid/rhoCfluid;
        DMm=(1-RV)/dt;
        % Compute enthalpy transfer term
        DHm=Gmass*dHtotal;
        % Save new composition and melt fraction
        XWsolidm(m)=XWsolidm1; % solid composition
        phinewm(m)=phinew1; % porosity
        % Reset properties at the first timestep
        if(timestep==1)
            XWsolidm0(m)=XWsolidm(m);
            phim(m)=phinewm(m);
        end

        % Interpolation to pressure nodes 
        % Update subgrid diffusion on nodes
        % i,j Node
        DMPSUM(i,j)=DMPSUM(i,j)+DMm*wtmij;
        DHPSUM(i,j)=DHPSUM(i,j)+DHm*wtmij;
        WTPSUM(i,j)=WTPSUM(i,j)+wtmij;
        % i+1,j Node
        DMPSUM(i+1,j)=DMPSUM(i+1,j)+DMm*wtmi1j;
        DHPSUM(i+1,j)=DHPSUM(i+1,j)+DHm*wtmi1j;
        WTPSUM(i+1,j)=WTPSUM(i+1,j)+wtmi1j;
        % i,j+1 Node
        DMPSUM(i,j+1)=DMPSUM(i,j+1)+DMm*wtmij1;
        DHPSUM(i,j+1)=DHPSUM(i,j+1)+DHm*wtmij1;
        WTPSUM(i,j+1)=WTPSUM(i,j+1)+wtmij1;
        % i+1,j+1 Node
        DMPSUM(i+1,j+1)=DMPSUM(i+1,j+1)+DMm*wtmi1j1;
        DHPSUM(i+1,j+1)=DHPSUM(i+1,j+1)+DHm*wtmi1j1;
        WTPSUM(i+1,j+1)=WTPSUM(i+1,j+1)+wtmi1j1;
    end
end
end
% P-nodes
DMP=zeros(Ny1,Nx1);
DHP=zeros(Ny1,Nx1);
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTPSUM(i,j)>0)
            DMP(i,j)=DMPSUM(i,j)/WTPSUM(i,j);
            DHP(i,j)=DHPSUM(i,j)/WTPSUM(i,j);
        end
    end
end
    
    
    
% Hydro-Mechanical Solution
pf0=pf;
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
            % 1*Vx=0
            L(kvx,kvx)=1; % Left part
            R(kvx)=0; % Right part
            % Top boundary
            if(i==1 && j>1 && j<Nx)
                L(kvx,kvx+6)=bctop; % Left part
            end
            % Bottom boundary
            if(i==Ny1 && j>1 && j<Nx)
                L(kvx,kvx-6)=bcbottom; % Left part
            end
        else
        % Internal points: x-Stokes eq.
        % ETA*(d2Vx/dx^2+d2Vx/dy^2)-dP/dx=0
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
        % Viscosity points
        ETA1=ETA(i-1,j);
        ETA2=ETA(i,j);
        ETAP1=ETAP(i,j);
        ETAP2=ETAP(i,j+1);
        % Density gradients
        dRHOdx=(RHOX(i,j+1)-RHOX(i,j-1))/2/dx;
        dRHOdy=(RHOX(i+1,j)-RHOX(i-1,j))/2/dy;
        % Left part
        L(kvx,kvx-Ny1*6)=ETAP1/dx^2; % Vx1
        L(kvx,kvx-6)=ETA1/dy^2; % Vx2
        L(kvx,kvx)=-(ETAP1+ETAP2)/dx^2-...
                      (ETA1+ETA2)/dy^2-...
                      dRHOdx*gx*dt; % Vx3
        L(kvx,kvx+6)=ETA2/dy^2; % Vx4
        L(kvx,kvx+Ny1*6)=ETAP2/dx^2; % Vx5
        L(kvx,kvy)=ETAP1/dx/dy-ETA2/dx/dy-dRHOdy*gx*dt/4;  % Vy2
        L(kvx,kvy+Ny1*6)=-ETAP2/dx/dy+ETA2/dx/dy-dRHOdy*gx*dt/4;  % Vy4
        L(kvx,kvy-6)=-ETAP1/dx/dy+ETA1/dx/dy-dRHOdy*gx*dt/4;  % Vy1
        L(kvx,kvy+Ny1*6-6)=ETAP2/dx/dy-ETA1/dx/dy-dRHOdy*gx*dt/4;  % Vy3
        L(kvx,kpm)=pscale/dx; % P1
        L(kvx,kpm+Ny1*6)=-pscale/dx; % P2
        % Right part
        R(kvx)=-RHOX(i,j)*gx;
        end
        
        % Vy equation External points
        if(j==1 || j==Nx1 || i==1 || i==Ny || i==Ny1)
            % Boundary Condition
            % 1*Vy=0
            L(kvy,kvy)=1; % Left part
            R(kvy)=0; % Right part
            % Left boundary
            if(j==1 && i>1 && i<Ny)
                L(kvy,kvy+6*Ny1)=bcleft; % Left part
            end
            % Right boundary
            if(j==Nx1 && i>1 && i<Ny)
                L(kvy,kvy-6*Ny1)=bcright; % Left part
            end
        else
        % Internal points: y-Stokes eq.
        % ETA*(d2Vy/dx^2+d2Vy/dy^2)-dP/dy=-RHO*gy
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
        % Viscosity points
        % Viscosity points
        ETA1=ETA(i,j-1);
        ETA2=ETA(i,j);
        ETAP1=ETAP(i,j);
        ETAP2=ETAP(i+1,j);
        % Density gradients
        dRHOdx=(RHOY(i,j+1)-RHOY(i,j-1))/2/dx;
        dRHOdy=(RHOY(i+1,j)-RHOY(i-1,j))/2/dy;
        % Left part
        L(kvy,kvy-Ny1*6)=ETA1/dx^2; % Vy1
        L(kvy,kvy-6)=ETAP1/dy^2; % Vy2
        L(kvy,kvy)=-(ETAP1+ETAP2)/dy^2-...
                      (ETA1+ETA2)/dx^2-...
                      dRHOdy*gy*dt; % Vy3
        L(kvy,kvy+6)=ETAP2/dy^2; % Vy4
        L(kvy,kvy+Ny1*6)=ETA2/dx^2; % Vy5
        L(kvy,kvx)=ETAP1/dx/dy-ETA2/dx/dy-dRHOdx*gy*dt/4; %Vx3
        L(kvy,kvx+6)=-ETAP2/dx/dy+ETA2/dx/dy-dRHOdx*gy*dt/4; %Vx4
        L(kvy,kvx-Ny1*6)=-ETAP1/dx/dy+ETA1/dx/dy-dRHOdx*gy*dt/4; %Vx1
        L(kvy,kvx+6-Ny1*6)=ETAP2/dx/dy-ETA1/dx/dy-dRHOdx*gy*dt/4; %Vx2
        L(kvy,kpm)=pscale/dy; % P1
        L(kvy,kpm+6)=-pscale/dy; % P2
        
        % Right part
        R(kvy)=-RHOY(i,j)*gy;
        end
        
        % Ptotal equation External points
        if(i==1 || j==1 || i==Ny1 || j==Nx1 ||...
          (i==2 && j>=2 && j<=Nx))
            % Boundary Condition
            % 1*P=0
            L(kpm,kpm)=1; % Left part
            R(kpm)=0; % Right part
            % Real BC
            if(i==2 && j>=2 && j<=Nx)
                L(kpm,kpm)=1*pscale; %Left part
                R(kpm)=1e+5; % Right part
            end
        else
        % Internal points: continuity eq.
        % dVx/dx+dVy/dy+(Ptotal-Pfluid)/ETHAphi=0
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
        L(kpm,kpm)=pscale/(1-PHI(i,j))/ETAPHI(i,j); % Ptotal
        L(kpm,kpf)=-pscale/(1-PHI(i,j))/ETAPHI(i,j); % Pfluid
        % Right part
        R(kpm)=DMP(i,j);
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
        R(kqx)=min(RHOFX(i,j),RHOX(i,j))*gx;
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
        L(kqy,kpf+6)=pscale/dy; % P2
        % Right part
        R(kqy)=min(RHOFY(i,j),RHOY(i,j))*gy;
        end
        
        % Pfluid equation External points
        if(i==1 || j==1 || i==Ny1 || j==Nx1 ||...
          (i==2 && j>=2 && j<=Nx))
            % Boundary Condition
            % 1*Pfluid=0
            L(kpf,kpf)=1; % Left part
            R(kpf)=0; % Right part
            % Real BC
            if(i==2 && j>=2 && j<=Nx)
                L(kpf,kpf)=1*pscale; %Left part
                R(kpf)=1e+5; % Right part
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
        L(kpf,kpm)=-pscale/(1-PHI(i,j))/ETAPHI(i,j); % Ptotal
        L(kpf,kpf)=pscale/(1-PHI(i,j))/ETAPHI(i,j); % Pfluid
        % Right part
        R(kpf)=0;
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

% Compute Stress and strain rate components
% Compute EPSILONxy, SIGMAxy in basic nodes
EXY=zeros(Ny,Nx); % Strain rate EPSILONxy, 1/s
SXY=zeros(Ny,Nx); % Strain rate SIGMAxy, Pa
for j=1:1:Nx
    for i=1:1:Ny
        % EXY,SXY
        EXY(i,j)=0.5*((vx(i+1,j)-vx(i,j))/dy+...
            (vy(i,j+1)-vy(i,j))/dx);
        SXY(i,j)=2*ETA(i,j)*EXY(i,j);
    end
end
% Compute EPSILONxx, SIGMA'xx in pressure nodes
EXX=zeros(Ny1,Nx1); % Strain rate EPSILONxx, 1/s
SXX=zeros(Ny1,Nx1); % Strain rate SIGMAxx, Pa
DIVV=zeros(Ny1,Nx1); % div(v)
for j=2:1:Nx
    for i=2:1:Ny
        % DIVV
        DIVV(i,j)=(vx(i,j)-vx(i,j-1))/dx+(vy(i,j)-vy(i-1,j))/dy;
        % EXX
        EXX(i,j)=((vx(i,j)-vx(i,j-1))/dx-(vy(i,j)-vy(i-1,j))/dy)/2;
        % SXX
        SXX(i,j)=2*ETAP(i,j)*EXX(i,j);
    end
end

% Compute Dln((1-PHI)/PHI)/Dt
APHI=zeros(Ny1,Nx1);
aphimax=0;
for j=2:1:Nx
    for i=2:1:Ny
        APHI(i,j)=((pr(i,j)-pf(i,j))/(1-PHI(i,j))/ETAPHI(i,j))/PHI(i,j);
        aphimax=max(aphimax,abs(APHI(i,j)));
    end
end

% Apply Symmetry
% Top
APHI(1,2:Nx)=APHI(2,2:Nx);    
PHI(1,2:Nx)=PHI(2,2:Nx);    
pr(1,2:Nx)=pr(2,2:Nx);    
pf(1,2:Nx)=pf(2,2:Nx);    
% Bottom
APHI(Ny1,2:Nx)=APHI(Ny,2:Nx);    
PHI(Ny1,2:Nx)=PHI(Ny,2:Nx);    
pr(Ny1,2:Nx)=pr(Ny,2:Nx);    
pf(Ny1,2:Nx)=pf(Ny,2:Nx);    
% Left
APHI(:,1)=APHI(:,2);    
PHI(:,1)=PHI(:,2);    
pr(:,1)=pr(:,2);    
pf(:,1)=pf(:,2);    
% Right
APHI(:,Nx1)=APHI(:,Nx);    
PHI(:,Nx1)=PHI(:,Nx);    
pr(:,Nx1)=pr(:,Nx);    
pf(:,Nx1)=pf(:,Nx); 

% Compute solid pressure
ps=(pr-pf.*PHI)./(1-PHI);

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

% Compute shear heating HS in pressure nodes
HS=zeros(Ny1,Nx1); % Adiabatic heating, W/m^3
for j=2:1:Nx
    for i=2:1:Ny
        % Average SXY*EXY
        SXYEXY=(SXY(i,j)*EXY(i,j)+SXY(i-1,j)*EXY(i-1,j)+...
            SXY(i,j-1)*EXY(i,j-1)+SXY(i-1,j-1)*EXY(i-1,j-1))/4;
        % HS
        HS(i,j)=2*SXX(i,j)*EXX(i,j)+2*SXYEXY+...
            (pr(i,j)-pf(i,j))^2/(1-PHI(i,j))/ETAPHI(i,j)+...
            (RX(i,j-1)*qxD(i,j-1)^2+RX(i,j)*qxD(i,j)^2)/2+...
            (RY(i-1,j)*qyD(i-1,j)^2+RY(i,j)*qyD(i,j)^2)/2;
    end
end

% Compute total adiabatic heating HA*dt in pressure nodes
HA=zeros(Ny1,Nx1); % Shear heating, W/m^3
for j=2:1:Nx
    for i=2:1:Ny
        % HA
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
        % HA
        HA(i,j)=(1-PHI(i,j))*tk1(i,j)*ALPHA(i,j)*dpsdt+...
            PHI(i,j)*tk1(i,j)*alphafluid*dpfdt;
    end
end




% Define timestep
if(titer==1)
    dt=dt*dtkoef;
    % Solid velocity
    maxvx=max(max(abs(vx)));
    maxvy=max(max(abs(vy)));
    if(dt*maxvx>dxymax*dx)
        dt=dxymax*dx/maxvx;
    end
    if(dt*maxvy>dxymax*dy)
        dt=dxymax*dy/maxvy;
    end
    % Fluid velocity
    maxvxf=max(max(abs(vxf)));
    maxvyf=max(max(abs(vyf)));
    if(dt*maxvxf>dxymax*dx)
        dt=dxymax*dx/maxvxf;
    end
    if(dt*maxvyf>dxymax*dy)
        dt=dxymax*dy/maxvyf;
    end
    % Porosity change
    if(aphimax*dt>dphimax)
        dt=dphimax/aphimax;
    end
    dt
end


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
                LT(gk,gk+1)=1; % Left part
                RT(gk)=273*2; % Right part
            end
            % Bottom BC: T=1500
            if(i==Ny1 && j>1 && j<Nx1)
                LT(gk,gk)=1; % Left part
                LT(gk,gk-1)=1; % Left part
                RT(gk)=1500*2; % Right part
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
        LT(gk,gk)=RHOCP(i,j)/dt+(Kx1+Kx2)/dx^2+(Ky1+Ky2)/dy^2; % FI3
        LT(gk,gk+1)=-Ky2/dy^2; % FI4
        LT(gk,gk+Ny1)=-Kx2/dx^2; % FI5
        % Right part
        RT(gk)=RHOCP(i,j)/dt*tk1(i,j)+HR(i,j)+HA(i,j)+HS(i,j)+DHP(i,j);
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
dt
% Apply thermal timestepping condition
if(titer==1)
    maxDTcurrent=max(max(abs(DT)));
    if(maxDTcurrent>DTmax)
        dt=dt/maxDTcurrent*DTmax;
    end
    dt
end


figure(3);colormap('Jet');clf
subplot(4,4,1)
pcolor(x/1000,y/1000,log10(ETA));% caxis([17 21])
shading flat;
axis ij image;
colorbar
title(['logETA,Pa*s iter=',num2str(titer)])
hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'w')

subplot(4,4,2)
pcolor(xp/1000,yp/1000,pr)
shading interp;
axis ij image;
colorbar
title(['Ptotal,Pa time(My)=',num2str(timesum/1e+6/(365.25*24*33600))])
hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(4,4,3)
pcolor(xvx/1000,yvx/1000,vx)
shading interp;
axis ij image;
colorbar
title('vx, m/s')
title(['vx,m/s dt(y)=',num2str(dt/(365.25*24*33600))])
hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'k')

subplot(4,4,4)
pcolor(xvy/1000,yvy/1000,vy)
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
pcolor(xp/1000,yp/1000,pf-pf0)
shading interp;
axis ij image;
colorbar
title('pf-pf0,Pa')

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
pcolor(xp/1000,yp/1000,log10(PHI))
shading interp;
axis ij image;
colorbar
title('logPHI')

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
title('qxDarcy, m/s')

subplot(4,4,13)
pcolor(xvy/1000,yvy/1000,qyD)
shading interp;
axis ij image;
colorbar
title('qyDarcy, m/s')


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

% subplot(4,4,16)
% pcolor(xp/1000,yp/1000,APHI)
% shading interp;
% axis ij image;
% colorbar
% title('Dln((1-PHI)/PHI)/dt')
subplot(4,4,16)
pcolor(xp/1000,yp/1000,DMP)
shading interp;
axis ij image;
colorbar
title('DMP, 1/s')
% subplot(4,4,16)
% pcolor(xp/1000,yp/1000,DHP)
% shading interp;
% axis ij image;
% colorbar
% title('DHP, J/m^3/s')


% figure(2);colormap('Jet');clf
% plot(xm,ym,'. k');
% axis ij image

pause(0.01)




% Stop iterations
pferrcur=max(max(abs(pf-pf0)));
DMPmax=max(max(abs(DMP)));
if(pferrcur<pferrmax && (titer>2 || DMPmax<=0))
    break;
end

end
DT0=DT;


% Apply subgrid diffusion on markers
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
    if(tm(m)<4)
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
    dtkm1=dtkm0*exp(-dsubgridt*ktotalm*dt/rhocptotalm*(2/dx^2+2/dy^2));
    % Correct marker temperature
    ddtkm=dtkm1-dtkm0;
    tkm(m)=tkm(m)+ddtkm;
    % Update subgrid diffusion on nodes
    % i,j Node
    TKSUM(i,j)=TKSUM(i,j)+ddtkm*rhocptotalm*wtmij;
    RHOCPSUM(i,j)=RHOCPSUM(i,j)+rhocptotalm*wtmij;
    % i+1,j Node
    TKSUM(i+1,j)=TKSUM(i+1,j)+ddtkm*rhocptotalm*wtmi1j;
    RHOCPSUM(i+1,j)=RHOCPSUM(i+1,j)+rhocptotalm*wtmi1j;
    % i,j+1 Node
    TKSUM(i,j+1)=TKSUM(i,j+1)+ddtkm*rhocptotalm*wtmij1;
    RHOCPSUM(i,j+1)=RHOCPSUM(i,j+1)+rhocptotalm*wtmij1;
    % i+1,j+1 Node
    TKSUM(i+1,j+1)=TKSUM(i+1,j+1)+ddtkm*rhocptotalm*wtmi1j1;
    RHOCPSUM(i+1,j+1)=RHOCPSUM(i+1,j+1)+rhocptotalm*wtmi1j1;
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
% Update composition and porosity for melting
XWsolidm0=XWsolidm;
phim=phinewm;
% Update porosity for compaction
for m=1:1:marknum
    if(tm(m)<4)
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
        phim(m)=phim(m)/((1-phim(m))*exp(aphim*dt)+phim(m));
        if(phim(m)<phimin)
            phim(m)=phimin;
        elseif(phim(m)>phimax)
            phim(m)=phimax;
        end
    end
end
phinewm=phim;



% Interpolate melt composition to pressure nodes
XWSSUM=zeros(Ny1,Nx1);
WTPSUM=zeros(Ny1,Nx1);
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
    % Update fluid composition on nodes
    % i,j Node
    XWSSUM(i,j)=XWSSUM(i,j)+XWsolidm0(m)*wtmij;
    WTPSUM(i,j)=WTPSUM(i,j)+wtmij;
    % i+1,j Node
    XWSSUM(i+1,j)=XWSSUM(i+1,j)+XWsolidm0(m)*wtmi1j;
    WTPSUM(i+1,j)=WTPSUM(i+1,j)+wtmi1j;
    % i,j+1 Node
    XWSSUM(i,j+1)=XWSSUM(i,j+1)+XWsolidm0(m)*wtmij1;
    WTPSUM(i,j+1)=WTPSUM(i,j+1)+wtmij1;
    % i+1,j+1 Node
    XWSSUM(i+1,j+1)=XWSSUM(i+1,j+1)+XWsolidm0(m)*wtmi1j1;
    WTPSUM(i+1,j+1)=WTPSUM(i+1,j+1)+wtmi1j1;
end
% Compute DTsubgrid
XWS=zeros(Ny1,Nx1);
% P-nodes
for j=1:1:Nx1
    for i=1:1:Ny1
        if(WTPSUM(i,j)>0)
            XWS(i,j)=XWSSUM(i,j)/WTPSUM(i,j);
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
vxpf(:,1)=-vxpf(:,2);
% Right
vxpf(:,Nx1)=-vxpf(:,Nx);
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
vypf(1,:)=-vypf(2,:);
% Bottom
vypf(Ny1,:)=-vypf(Ny,:);



% Compute solid velocity in pressure nodes
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
vxp(:,1)=-vxp(:,2);
% Right
vxp(:,Nx1)=-vxp(:,Nx);
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
vyp(2:Ny-1,Nx1)=-bcright*vyp(2:Ny-1,Nx);     
% Top
vyp(1,:)=-vyp(2,:);
% Bottom
vyp(Ny1,:)=-vyp(Ny,:);

% Compute fluid velocity in pressure nodes
% vx
for j=2:1:Nx
    for i=2:1:Ny
        vxpf(i,j)=(vxf0(i,j)+vxf0(i,j-1))/2;
    end
end
% Apply BC
% Top
vxpf(1,2:Nx-1)=-bcftop*vxpf(2,2:Nx-1);    
% Bottom
vxpf(Ny1,2:Nx-1)=-bcfbottom*vxpf(Ny,2:Nx-1);    
% Left
vxpf(:,1)=-vxpf(:,2);
% Right
vxpf(:,Nx1)=-vxpf(:,Nx);
% vy
for j=2:1:Nx
    for i=2:1:Ny
        vypf(i,j)=(vyf0(i,j)+vyf0(i-1,j))/2;
    end
end    
% Apply BC
% Left
vypf(2:Ny-1,1)=-bcfleft*vypf(2:Ny-1,2);    
% Right
vypf(2:Ny-1,Nx1)=-bcfright*vypf(2:Ny-1,Nx);    
% Top
vypf(1,:)=-vypf(2,:);
% Bottom
vypf(Ny1,:)=-vypf(Ny,:);
% Add solid velocity
vxpf=vxpf+vxp;
vypf=vypf+vyp;

% Move markers with 4th order Runge-Kutta
vxm=zeros(4,1);
vym=zeros(4,1);
for m=1:1:marknum
    
    % Interpolate solid temperature, XBfluid for the initial marker location
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
    % Compute Tsolid, XBfluid
    tksm0=tk2(i,j)*wtmij+tk2(i+1,j)*wtmi1j+...
            tk2(i,j+1)*wtmij1+tk2(i+1,j+1)*wtmi1j1;
    
    
    % Save initial marker coordinates
    xA=xm(m);
    yA=ym(m);
    for rk=1:1:4
        % Interpolate vxp,vyp
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
        % Compute vx, vy velocity
        vxm(rk)=vxp(i,j)*wtmij+vxp(i+1,j)*wtmi1j+...
            vxp(i,j+1)*wtmij1+vxp(i+1,j+1)*wtmi1j1;
        vym(rk)=vyp(i,j)*wtmij+vyp(i+1,j)*wtmi1j+...
            vyp(i,j+1)*wtmij1+vyp(i+1,j+1)*wtmi1j1;
        
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
        wtmij=(1-dxmj/dx)*(1-dymi/dy);
        wtmi1j=(1-dxmj/dx)*(dymi/dy);    
        wtmij1=(dxmj/dx)*(1-dymi/dy);
        wtmi1j1=(dxmj/dx)*(dymi/dy);
        % Compute vx velocity
        vxm(rk)=vpratio*vxm(rk)+(1-vpratio)*(vx(i,j)*wtmij+vx(i+1,j)*wtmi1j+...
            vx(i,j+1)*wtmij1+vx(i+1,j+1)*wtmi1j1);
        
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
        wtmij=(1-dxmj/dx)*(1-dymi/dy);
        wtmi1j=(1-dxmj/dx)*(dymi/dy);    
        wtmij1=(dxmj/dx)*(1-dymi/dy);
        wtmi1j1=(dxmj/dx)*(dymi/dy);
        % Compute vx velocity
        vym(rk)=vpratio*vym(rk)+(1-vpratio)*(vy(i,j)*wtmij+vy(i+1,j)*wtmi1j+...
            vy(i,j+1)*wtmij1+vy(i+1,j+1)*wtmi1j1);        
        
        % Change coordinates to obtain B,C,D points
        if(rk==1 || rk==2)
            xm(m)=xA+dt/2*vxm(rk);
            ym(m)=yA+dt/2*vym(rk);
        elseif(rk==3)
            xm(m)=xA+dt*vxm(rk);
            ym(m)=yA+dt*vym(rk);
        end
    end
    % Restore initial coordinates
    xm(m)=xA;
    ym(m)=yA;
    % Compute effective velocity
    vxmeff=1/6*(vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4));
    vymeff=1/6*(vym(1)+2*vym(2)+2*vym(3)+vym(4));
    % Move markers
    xm(m)=xm(m)+dt*vxmeff;
    ym(m)=ym(m)+dt*vymeff;
    
    % Backtracing markers with fluid velocity
    xcur=xm(m);
    ycur=ym(m);
    xA=xcur;
    yA=ycur;
    for rk=1:1:4
        % Interpolate vxpf,vypf
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
        % Compute vx, vy velocity
        vxm(rk)=vxpf(i,j)*wtmij+vxpf(i+1,j)*wtmi1j+...
            vxpf(i,j+1)*wtmij1+vxpf(i+1,j+1)*wtmi1j1;
        vym(rk)=vypf(i,j)*wtmij+vypf(i+1,j)*wtmi1j+...
            vypf(i,j+1)*wtmij1+vypf(i+1,j+1)*wtmi1j1;
        
        % Interpolate vxf
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
        wtmij=(1-dxmj/dx)*(1-dymi/dy);
        wtmi1j=(1-dxmj/dx)*(dymi/dy);    
        wtmij1=(dxmj/dx)*(1-dymi/dy);
        wtmi1j1=(dxmj/dx)*(dymi/dy);
        % Compute vx velocity
        vxm(rk)=vpratio*vxm(rk)+(1-vpratio)*(vxf(i,j)*wtmij+vxf(i+1,j)*wtmi1j+...
            vxf(i,j+1)*wtmij1+vxf(i+1,j+1)*wtmi1j1);
        
        % Interpolate vyf
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
        wtmij=(1-dxmj/dx)*(1-dymi/dy);
        wtmi1j=(1-dxmj/dx)*(dymi/dy);    
        wtmij1=(dxmj/dx)*(1-dymi/dy);
        wtmi1j1=(dxmj/dx)*(dymi/dy);
        % Compute vx velocity
        vym(rk)=vpratio*vym(rk)+(1-vpratio)*(vyf(i,j)*wtmij+vyf(i+1,j)*wtmi1j+...
            vyf(i,j+1)*wtmij1+vyf(i+1,j+1)*wtmi1j1);        
        
        % Change coordinates to obtain B,C,D points
        if(rk==1 || rk==2)
            xcur=xA-dt/2*vxm(rk);
            ycur=yA-dt/2*vym(rk);
        elseif(rk==3)
            xcur=xA-dt*vxm(rk);
            ycur=yA-dt*vym(rk);
        end
    end
    % Compute effective velocity
    vxmeff=1/6*(vxm(1)+2*vxm(2)+2*vxm(3)+vxm(4));
    vymeff=1/6*(vym(1)+2*vym(2)+2*vym(3)+vym(4));
    % Trace the node backward
    xcur=xA-dt*vxmeff;
    ycur=yA-dt*vymeff;
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
        phim(m)*(tkm(m)+dtkfsm)*rhocpfluid)/...
        ((1-phim(m))*rhocpsolidm(tm(m))+phim(m)*rhocpfluid);
end  

% Update timesum
timesum=timesum+dt;



figure(1);colormap('Jet');clf
subplot(4,4,1)
pcolor(x/1000,y/1000,log10(ETA));% caxis([17 21])
shading flat;
axis ij image;
colorbar
title('log10ETA, Pa*s')
hold on
quiver(xp(3:5:Nx1)/1000,yp(3:5:Ny1)/1000,vxp(3:5:Ny,3:5:Nx1),vyp(3:5:Ny1,3:5:Nx1),'w')

subplot(4,4,2)
pcolor(xp/1000,yp/1000,pr)
shading interp;
axis ij image;
colorbar
title('Ptotal, Pa')
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
pcolor(xp/1000,yp/1000,RHO)
shading interp;
axis ij image;
colorbar
title('RHOtotal, kg/m^3')

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
pcolor(xp/1000,yp/1000,log10(PHI))
shading interp;
axis ij image;
colorbar
title('logPHI')

subplot(4,4,11)
pcolor(xp/1000,yp/1000,pf)
shading interp;
axis ij image;
colorbar
title('Pfluid, Pa')

% subplot(4,4,12)
% pcolor(xvx/1000,yvx/1000,qxD)
% shading interp;
% axis ij image;
% colorbar
% title('qxDarcy, m/s')
subplot(4,4,12)
pcolor(xvx/1000,yvx/1000,vxf)
shading interp;
axis ij image;
colorbar
title('Vxfluid, m/s')
 
% subplot(4,4,13)
% pcolor(xvy/1000,yvy/1000,qyD)
% shading interp;
% axis ij image;
% colorbar
% title('qyDarcy, m/s')
subplot(4,4,13)
pcolor(xvy/1000,yvy/1000,vyf)
shading interp;
axis ij image;
colorbar
title('VyFluid, m/s')

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

% subplot(4,4,16)
% pcolor(xp/1000,yp/1000,APHI)
% shading interp;
% axis ij image;
% colorbar
% title('Dln((1-PHI)/PHI)/dt')
% subplot(4,4,16)
% pcolor(xp/1000,yp/1000,DMP)
% shading interp;
% axis ij image;
% colorbar
% title('DMP, 1/s')
% subplot(4,4,16)
% pcolor(xp/1000,yp/1000,DHP)
% shading interp;
% axis ij image;
% colorbar
% title('DHP, J/m^3/s')
subplot(4,4,16)
pcolor(xp/1000,yp/1000,log10(XWS))
shading interp;
axis ij image;
colorbar
title('logXWS, 1/s')


% figure(2);colormap('Jet');clf
% plot(xm,ym,'. k');
% axis ij image

pause(0.01)


end

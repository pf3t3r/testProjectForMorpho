clear; close all; clc;

% Solving the shallow water equations in shallow basins and estuaries.
% Based on Friedrichs (2011), Chapter 3 in Contemporary Issues in Estuarine
% Physics.

%**************************************************************************
%**************************************************************************
%*              Parameter settings
%**************************************************************************
%**************************************************************************

deltaT=45;               % time step in seconds. Choose appropriate time step yourself based on Courant number. 
deltaX=1000;             % spatial step in meters
Lbasin=1.1e5;             % Length of the basin or estuary in meters
Lb=4e4;                  % e-folding length scale for width.
B0=1e3;                  % Width of the basin in meters at seaward side.
H0=10;                   % Depth of basin.
M2amp=1;                 % Amplitude of M2 tide at seaward side.
discharge=0;             % Constant river discharge at landward boundary. 
Cd=1e-3*[0.2, 0.5, 1, 2, 5];               % Drag coefficients

%**************************************************************************
%**************************************************************************
%*                 End of parameter setting
%**************************************************************************
%**************************************************************************

Td1=24*3600+50*60;
Tm2=12*3600+25*60;       % M2 tidal period in seconds
time=0:deltaT:15*Tm2;    % time in seconds
Nt=length(time);

% Define frequencies to be analysed.
global wn
wn(1)=2*pi/Td1;
wn(2)=2*pi/Tm2;
wn(3)=2*wn(2);
wn(4)=3*wn(2);

% No. of different possible drag coefficients 
N_Cd=length(Cd);

for i=1:N_Cd
x=0:deltaX:Lbasin;
Nx=length(x);
%B(1:Nx)=B0*exp(-x/Lb);
B(1:Nx)=B0;                % when basin width has to be constant.
H(1:Nx)=H0;
Z=zeros(Nx-1,Nt);           % Z points shifted half a grid point to the right with respect to Q points. we start with Q point, therefore 1 Z point less than Q points.      
Q=zeros(Nx,Nt);
A=(B.*H)'*ones(1,Nt);       % A at Q points
P=B'*ones(1,Nt);            % Wetted perimeter at Q points.
Inertia=zeros(Nx,Nt);       % Initalize Inertia, Pressure Gradient and Friction for further analysis later on.
PG=zeros(Nx,Nt);
Fric=zeros(Nx,Nt);

% Boundary conditions
Z(1,:)=M2amp*sin(2*pi*time/Tm2);         % prescribed water levels
Q(Nx,:)=-discharge;                      % river discharge; most often river discharge is taken zero.

courant=sqrt(9.8*max(H0))*deltaT/deltaX;

% Courant criteria
if courant<1
end 

if courant>=1
    disp('Courant criteria not met');
    break; 
end 

% For numerical part, follow thesis of Speer (1984). Staggered grid. Z points shifted half deltaX to
% the right of U points: 
% Q1 Z1 Q2 Z2 ........ ZN QN+1

% solve Bdz/dt + dQ/dx=0
% solve dQ/dt + d/dx Q^2/A = -gA dZ/dx - Cd Q|Q|P/A*A

% Z is water level, B=estuary width, Q is discharge, A is cross-sectional
% area, P is wetted perimeter (~channel width when H/B is small)

% First start simple: rectangular basin, no advection, no river flow.
% Numerical scheme from Speer (1984).

for pt=1:Nt-1
    for px=2:Nx-1
        Z(px,pt+1)=Z(px,pt)-(deltaT/(0.5*(B(px)+B(px+1))))*(Q(px+1,pt)-Q(px,pt))/deltaX;
    end
    for px=2:Nx-1
        A(px,pt+1)=B(px)*(H(px)+0.5*Z(px,pt+1)+0.5*Z(px-1,pt+1));           % A at Q points
        P(px,pt+1)=B(px)+2*H(px)+Z(px,pt+1) +Z(px-1,pt+1);                  % P at Q points
    end
    for px=2:Nx-1
        Q(px,pt+1)=Q(px,pt) ...                                            % Inertia.
            -9.81*A(px,pt+1)*(deltaT/deltaX)*(Z(px,pt+1)-Z(px-1,pt+1)) ...  % Pressure gradient
            -Cd(i)*deltaT*abs(Q(px,pt))*Q(px,pt)*P(px,pt)/(A(px,pt)*A(px,pt)); % Friction
        Inertia(px,pt+1)=(Q(px,pt+1)-Q(px,pt))/deltaT;
        PG(px,pt+1)=-9.81*A(px,pt+1)*(1/deltaX)*(Z(px,pt+1)-Z(px-1,pt+1));
        Fric(px,pt+1)=-Cd(i)*abs(Q(px,pt))*Q(px,pt)*P(px,pt)/(A(px,pt)*A(px,pt));
    end
    Q(1,pt+1)=Q(2,pt+1)+B(1)*deltaX*(Z(1,pt+1)-Z(1,pt))/deltaT;
end

U=Q./A;         % Flow velocity in m/s

% A3. Tidally averaged flow velocities
U_mean=0;
for pm=1:Nt-1
U_mean=U_mean+U(pm);
end
U_Mean=U_mean/Nt;

% Analyse last tidal period only. For example determine amplitude and phase of M2, M4, M6 and mean
% of water level and flow velocity. Design you own code here. I used my code (harmfit). You
% can determine HW level, LW level, moments of LW and HW, propagation speed of LW wave
% and HW wave... You can determine propagation speed of wave by determining
% the phase as a function of space. The phase (in radians) also equals k*x, so a linear
% fit to the phase determines wave number k. 

Nsteps=floor(Td1/deltaT);
for px=1:Nx-1
coefin=[0.1, 0.3, 1, 0.2, 0.1, 0.2, 1, 0.2, 0.1];
coefout=nlinfit(time(end-Nsteps:end),Z(px,end-Nsteps:end),@harmfit,coefin);
Z0(i,px)=coefout(1);
ZM2(i,px)=sqrt(coefout(3).^2+coefout(7).^2);
ZM4(i,px)=sqrt(coefout(4).^2+coefout(8).^2);
ZM6(i,px)=sqrt(coefout(5).^2+coefout(9).^2);
phaseZM2(i,px)=atan(coefout(3)/coefout(7));
phaseZM4(i,px)=atan(coefout(4)/coefout(8));
phaseZM6(i,px)=atan(coefout(5)/coefout(9));
coefin=[0.1, 0.3, 1, 0.2, 0.1, 0.2, 1, 0.2, 0.1];
coefout=nlinfit(time(end-Nsteps:end),U(px,end-Nsteps:end),@harmfit,coefin);
U0(px)=coefout(1);
UM2(i,px)=sqrt(coefout(3).^2+coefout(7).^2);
UM4(i,px)=sqrt(coefout(4).^2+coefout(8).^2);
UM6(i,px)=sqrt(coefout(5).^2+coefout(9).^2);
phaseUM2(i,px)=atan(coefout(3)/coefout(7));
phaseUM4(i,px)=atan(coefout(4)/coefout(8));
phaseUM6(i,px)=atan(coefout(5)/coefout(9));
end
%We analyse the size of kL - lenght of estuary*length scale over which
%tidal phase varies. 
LkM2(i)=(phaseUM2(1)+phaseUM2(Nx-1));
%Looking at the size of kL, we predict that the pumping model is
%representative for M2 in this example, except for when the water height
%<=2m. 
end

% Make a plot of how SSE varies across the length of the basin when the
% using different values for the drag coefficient Cd.

Matlab2_B2=figure;
yyaxis left;
plot(x(2:end),ZM2);
ylabel('ZM2 [m]');
hold on
yyaxis right;
plot(x(2:end),UM2);
hold off
title('M2: Dependency of SSE on c_d (L_{Basin} = 110 km)');
xlabel('L_{Basin} [m]');
ylabel('UM2 [m/s]');
legend('C_d = 0.0002','C_d = 0.0005','C_d = 0.0010','C_d = 0.0020','C_d = 0.0050');
grid on;
saveas(gcf,'Matlab2_B2.png');

% B2. What happens to the wavelength of the tide due to friction?
% -> Higher friction leads to shorter wavelengths (noticible from plot). 

% B2. In case of moderate friction, is the M2 tide resonant for larger basin
% lengths or smaller basin lengths compared to the case of negligible friction?
% -> Let's assume that moderate friction and negligible friction are
% represented by the cases where c_d = 0.0010 and c_d = 0.0002
% respectively. 
% For moderate friction, it appears as though the M2 tide is resonant for
% shorter basin lengths. Specifically, in the case of negligible friction,
% resonance is reached after approx. 104 km. Then, in the case of moderate
% friction, resonance is reached after approx. 90 km. This is due to the
% fact that friction shortens the wavelenght. 

% Give two reasons why the amplitude of the M2 water levels at the end of
% the basin are smaller for increasing values of Cd.
% -> 1. Higher friction damps the signal (energy loss). Both decreased
% amplitude and speed
% -> 2. In this specific case, friction impacts resonance too, as now
% resonance is only present in a shorter basin. 
% -> 3. Much of the energy of M2 has been transferred to M4. [true, but not
% complete].


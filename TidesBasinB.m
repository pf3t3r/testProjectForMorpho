clear all; close all;
clc;

% Solving the shallow water equations in shallow basins and estuaries.
% Based on Friedrichs (2011), Chapter 3 in Contemporary Issues in Estuarine
% Physics.

%**************************************************************************
%**************************************************************************
%*              Paremeter settings
%**************************************************************************
%**************************************************************************

deltaT=45;               % time step in seconds. Choose appropriate time step yourself based on Courant number. 
deltaX=1000;             % spatial step in meters
Lbasin=linspace(1e4,15e4,8); % Length of the basin or estuary in meters
Lb=4e4;                  % e-folding length scale for width.
B0=1e3;                  % Width of the basin in meters at seaward side.
H0=10;                   % Depth of basin.
M2amp=1;                 % Amplitude of M2 tide at seaward side.
discharge=0;             % Constant river discharge at landward boundary. 
Cd=0.2e-3;               % Drag coefficient

%**************************************************************************
%**************************************************************************
%*                 End of parameter setting
%**************************************************************************
%**************************************************************************

Td1=24*3600+50*60;
Tm2=12*3600+25*60;       % M2 tidal period in seconds
time=0:deltaT:15*Tm2;    % time in seconds
Nt=length(time);

%Define frequencies to be analysed. To determine the amplitudes and phase
%use the code you designed in the first Matlab exercise. 

global wn

wn(1)=2*pi/Td1;
wn(2)=2*pi/Tm2;
wn(3)=2*wn(2);
wn(4)=3*wn(2);

%Length of basin lenght variation 
NH=length(Lbasin);

f1=figure; 
for i=1:NH
x=0:deltaX:Lbasin(i);
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

%Courant criteria
if courant<1
end 

if courant>=1
    display('Courant criteria not met');
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
            -Cd*deltaT*abs(Q(px,pt))*Q(px,pt)*P(px,pt)/(A(px,pt)*A(px,pt)); % Friction
        Inertia(px,pt+1)=(Q(px,pt+1)-Q(px,pt))/deltaT;
        PG(px,pt+1)=-9.81*A(px,pt+1)*(1/deltaX)*(Z(px,pt+1)-Z(px-1,pt+1));
        Fric(px,pt+1)=-Cd*abs(Q(px,pt))*Q(px,pt)*P(px,pt)/(A(px,pt)*A(px,pt));
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

subplot(3,2,1)
plot(x(2:end),ZM2);
title('M2: Along-basin Change in SSE');
xlabel('L_{Basin} [m]');
ylabel('SSE [m]');
legend('L = 10km','L = 30km','L = 50km','L = 70km','L = 90km','L = 110km','L = 130km','L = 150km');
grid on;
end

% %sgtitle('Part A');
% subplot(3,2,1)
% plot(x(2:end),ZM2);
% title('M2: Along-basin Change in SSE');
% xlabel('L_{Basin} [m]');
% ylabel('SSE [m]');
% legend('H = 2m','H = 3m','H = 4m','H = 5m','H = 6m','H = 7m','H = 8m','H= 9m','H = 10m');
% grid on;

% subplot(3,2,3)
% plot(x(2:end),UM2);
% title('M2: Along-basin Change in U');
% xlabel('L_{Basin} [m]');
% ylabel('U [m/s]');
% legend('H = 2m','H = 3m','H = 4m','H = 5m','H = 6m','H = 7m','H = 8m','H= 9m','H = 10m');
% 
% grid on;
% 
% subplot(3,2,2)
% plot(x(2:end),phaseUM2-phaseZM2);
% title('M2: Phase difference UM2 - ZM2');
% xlabel('L_{Basin} [m]');
% ylabel('\Phi [rad]');
% legend('H = 2m','H = 3m','H = 4m','H = 5m','H = 6m','H = 7m','H = 8m','H= 9m','H = 10m');
% grid on;
% 
% subplot(3,2,4)
% plot(x(2:end),phaseZM2);
% title('M2: Along-basin Change in Phase');
% xlabel('L_{Basin} [m]');
% ylabel('\Phi [rad]');
% legend('H = 2m','H = 3m','H = 4m','H = 5m','H = 6m','H = 7m','H = 8m','H= 9m','H = 10m');
% grid on;
% 
% subplot(3,2,6)
% scatter(H0,LkM2);
% title('M2: Along-basin Change in Phase');
% xlabel('Lk');
% ylabel('Average water depth');
% grid on;

% subplot(3,2,2)
% yyaxis left;
% plot(x(2:end),ZM4(1,:));
% hold on
% plot(x(2:end),ZM4(end,:));
% ylabel('SSE [m]');
% yyaxis right;
% plot(x(2:end),ZM2(1,:));
% plot(x(2:end),ZM2(end,:));
% hold off
% title('M4 and M2: Deepest and Shallowest Case');
% xlabel('L_{Basin} [m]');
% ylabel('SSE [m]');
% legend('M4 (H = 2m)','M4 (H = 10m)','M2 (H = 2m)','M2 (H = 10m)');
% grid on;
% 
% subplot(3,2,4)
% plot(x(2:end),UM2(1,:));
% hold on
% plot(x(2:end),UM4(1,:));
% hold off
% title('Flow Velocities for M2 and M4');
% xlabel('L_{Basin} [m]');
% ylabel('U [m/s]');
% legend('U_{M2}','U_{M4}')
% grid on;


%A1. Explain the dependence of tidal amplitude on depth. 
%M2 tidal amplitude increases landward for heights above 3m. This is due to
%the negative correlation between friction and tidal amplitude. The
%shallower the water, the higher the levels of friction and so the more
%dampened is M2. This is opposite for M4 and M6 which are max in shallow
%waters. As can be seen in subplot(3,2,5), M4 diminishes in size as H
%growth, this is due to the decreased amount of interactions amongst M2
%tides due to the decreased friction. 

%A1. For which depths are you close to the pumping model solution? 
%Also remember that for pumping model the phase of M2 water 
%level should be uniform in the basin.

%As computed above, this linearisation is only possible for water heights 
%over 2m. This can also be observed in the graph above by comparing the
%amplitude and tidal speed phases. From equation (3.14) it follows that in 
%short estuaries, the linearized component of tidal velocity is 90° out of 
%phase with tidal elevation. This means that max U is when Z minimum.

% A2. Determine the deformation of the tide for the shallowest case (2 m)
% and deepest case (10 m) by determining the amplitude of the M4 water
% level as function of position in the basin. Explain the differences
% between a deep and shallow basin.

% The tide deforms significantly in the shallowest case (2m). Between 
% 0 and 800m, the M2 SSE declines in height by about 20%, while the 
% M4 SSE increases from zero to a height of 0.1 m. It continues thereafter,
% which seems to indicate that there is a phase lag between the loss in
% amplitude of M2 and the corresponding amplitude gain by M4.
% The M4 water level increases closer to the landward end of the basin in 
% both cases, though the amplitudes are much higher in the case of a
% shallow basin. For instance, in the case of a shallow basin (H = 2 m),
% the M4 water level increases from 0 m to 2 m, while in the case of a
% deep basin, the M4 water level increases only barely, from 0 m to 3 mm.

% A3. For the shallowest case: Determine the tidally averaged flow 
% velocities, and the amplitude of the M2 and the M4 tidal currents in 
% the basin and determine the relative phase difference between the M2
% and M4 tidal currents. 

% Average flow velocities are maximum at the seaward end of the basin.
% Velocity declines continuously for M2, but M4 begins to increase for 
% the first part of the basin.
% Is this correct? I don't think the M4 should have a velocity at zero (?)
clear; close all; clc;

%**************************************************************************
%**************************************************************************
%*              Parameter settings
%**************************************************************************
%**************************************************************************

deltaT=25;               % Time step in seconds. Chosen such that the
                         % the Courant number <= 0.9. 
                         % Warning: If deltaT is too high (but still 
                         % satisfies Courant), the harmfit function may
                         % break.
xr = 110e3;               % 50 = Point in river at which estuary stops narrowing
                         % kind of arbitrarily set to the location of Ile 
                         % de Patiras.
                         % 73: pt where rivers join
                         % 110: est splits into two rivers. In the two
                         % cases the width becomes constant after 107 and
                         % 111 km respectively (see fig), so we pick 110 as
                         % an average. 
Lbasin=1.5*xr;           % Length of the basin or estuary in meters
B0=5.8e3;                % Width of the basin in meters at seaward side.
                         % This is based on Google Earth observations of
                         % the mouth of the Gironde.
Briver=0.25e3;           % Width of the river where estuary stops narrowing
                         % 0.2 km width for Dordogne, 0.3 km for Garonne,
                         % so avg taken
deltaX=Lbasin/80;       % spatial step in meters.
Lb = -xr/log(Briver/B0); % e-folding length scale    
M2amp=1.55;              % Amplitude of M2 tide at seaward side <=> This is 
                         % the water level at the mouth.
discharge=0;             % Constant river discharge at landward boundary. 
                         % Prescribed as zero for simplicity in the project
                         % description.
Cd=2.5e-3;               % Drag coefficient
x=0:deltaX:Lbasin;       % With deltax = c.1.3e3, we have 80 grid points.

H0 = 8.5;                % Water Height 
NH0 = length(H0);
%**************************************************************************
%**************************************************************************
%*                 End of parameter setting
%**************************************************************************
%**************************************************************************

% Note that D1 tide is not used in this analysis.

% Td1=24*3600+50*60;       % M1 tidal period [s]
Tm2=12*3600+25*60;       % M2 tidal period [s]
time=0:deltaT:30*Tm2;    % Time [s]
Nt=length(time);         % Length of time array

% Define frequencies to be analysed.
global wn
% wn(1)=2*pi/Td1;         % M1
wn(1)=2*pi/Tm2;         % M2
wn(2)=2*wn(1);          % M4
wn(3)=3*wn(1);          % M6

% This loop must be set up such that it varies over bed level.
for i=1:NH0
Nx=length(x);
B(1:54)=B0*exp(-x(1:54)/Lb);      % Basin width for a converging channel until xr  
B(55:Nx)=B(54);                % basin width constant after xr
H(1:Nx)=H0(i);
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
            -Cd*deltaT*abs(Q(px,pt))*Q(px,pt)*P(px,pt)/(A(px,pt)*A(px,pt)); % Friction
        Inertia(px,pt+1)=(Q(px,pt+1)-Q(px,pt))/deltaT;
        PG(px,pt+1)=-9.81*A(px,pt+1)*(1/deltaX)*(Z(px,pt+1)-Z(px-1,pt+1));
        Fric(px,pt+1)=-Cd*abs(Q(px,pt))*Q(px,pt)*P(px,pt)/(A(px,pt)*A(px,pt));
    end
    Q(1,pt+1)=Q(2,pt+1)+B(1)*deltaX*(Z(1,pt+1)-Z(1,pt))/deltaT;
end

U=Q./A;         % Flow velocity in m/s

% D3. Tidally averaged flow velocities
U_mean=0;
for pm=1:Nt-1
U_mean=U_mean+U(:,pm);
end
U_Mean=(abs(U_mean/Nt))';

% Analyse last tidal period only. For example determine amplitude and phase of M2, M4, M6 and mean
% of water level and flow velocity. Design you own code here. I used my code (harmfit). You
% can determine HW level, LW level, moments of LW and HW, propagation speed of LW wave
% and HW wave... You can determine propagation speed of wave by determining
% the phase as a function of space. The phase (in radians) also equals k*x, so a linear
% fit to the phase determines wave number k. 

% Nsteps=floor(Td1/deltaT);
Nsteps=floor(Tm2/deltaT);
for px=1:Nx-1
    
% coefin=[0.1, 0.3, 1, 0.2, 0.1, 0.2, 1, 0.2, 0.1];
coefin=[0.1, 1, 0.2, 0.1, 1, 0.2, 0.1];
coefout=nlinfit(time(end-Nsteps:end),Z(px,end-Nsteps:end),@harmfit,coefin);
Z0(i,px)=coefout(1);
ZM2(i,px)=sqrt(coefout(2).^2+coefout(5).^2);
ZM4(i,px)=sqrt(coefout(3).^2+coefout(6).^2);
ZM6(i,px)=sqrt(coefout(4).^2+coefout(7).^2);
phaseZM2(i,px)=atan(coefout(2)/coefout(5));   
phaseZM4(i,px)=atan(coefout(3)/coefout(6));
phaseZM6(i,px)=atan(coefout(4)/coefout(7));
% coefin=[0.1, 0.3, 1, 0.2, 0.1, 0.2, 1, 0.2, 0.1];
coefin=[0.1, 1, 0.2, 0.1, 1, 0.2, 0.1];
coefout=nlinfit(time(end-Nsteps:end),U(px,end-Nsteps:end),@harmfit,coefin);
U0(px)=coefout(1);
UM2(i,px)=sqrt(coefout(2).^2+coefout(5).^2);
UM4(i,px)=sqrt(coefout(3).^2+coefout(6).^2);
UM6(i,px)=sqrt(coefout(4).^2+coefout(7).^2);
phaseUM2(i,px)=atan(coefout(2)/coefout(5));  
phaseUM4(i,px)=atan(coefout(3)/coefout(6));
phaseUM6(i,px)=atan(coefout(4)/coefout(7));
end
%We analyse the size of kL - lenght of estuary*length scale over which
%tidal phase varies. 
LkM2(i)=(phaseUM2(1)+phaseUM2(Nx-1));
%Looking at the size of kL, we predict that the pumping model is
%representative for M2 in this example, except for when the water height
%<=2m. 

%Linear fit for k
coefs(i,:) = polyfit(x(2:end-1), phaseZM2(i,2:end), 1);
end

%WE SHOULD USE ONLY 75% BUT GET WAY BETTER FIT FOR K,C IF WHOLE LENGTH, WHY?


%Fix for behaviour of atan for values >pi/2
phaseZM2(1,66:end)=phaseZM2(1,66:end)+pi; 
phaseUM2(1,35:end)=phaseUM2(1,35:end)+pi;
%coefs(1,:) = polyfit(x(2:end-1), phaseZM2(1,2:end), 1);

% D2
% Modeled values

k = (coefs(:,1))';
c = wn(1)./k;
%k = 1.3e-5;
%c = wn(1)./k;
UM2 = UM2(:,2:41);
phi = phaseZM2(:,2:41)-phaseUM2(:,2:41);

% Theoretical values
g = 9.81;
U_t=ones(40);
phi_t=ones(40);
r=Cd;
c_t = (g*H0)^(1/2);
k_t = wn(1)/c_t;
UM2_t = U_t*M2amp*c_t/((H0*(1+(wn(1)/(Cd/H0))^(-2))^(1/2))); %check a
phi_t = phi_t.*-atan((Cd/H0)/wn(1)); 

display(k);
display(k_t);
display(c);
display(c_t);

figure
plot(x(2:41),UM2_t);
hold on
plot(x(2:41),UM2);
hold off
title('Gironde Estuary: M2 Speed');
xlabel('Length [m]');
ylabel('UM2 [m/s]');
grid on;
saveas(gcf,'Matlab2_D2_i.png');

figure
plot(x(2:41),phi_t);
hold on
plot(x(2:41),phi);
hold off
title('Gironde Estuary: phase difference');
xlabel('Length [m]');
ylabel('phase [rad]');
grid on;
saveas(gcf,'Matlab2_D2_ii.png');

% D3. 
figure
plot(x(2:54),U_Mean(2:54));
title('Gironde Estuary: average flow velocities');
xlabel('Length [m]');
ylabel('Flow velocity [m/s]');
grid on;
saveas(gcf,'Matlab2_D2_iii.png');

% Average speed starts from zero at sea boundary and increases up to around 90% of 
% our estuary, and then slows back down to zero as we enter the river of 
% constant width. This means that across time, tidal flow is fastest at
% around 95km inland. This is not surprising??? Also understand why
% phaseUM2 behaves like that jumping. This value is the ratio of discharge
% to cross sectional area. This figure suggests that the ratio increases up
% to almost the end of the e-folding length scale. This means that there is
% higher water flow per volume at the landward portion of the estuary which
% makes sense as we have higher amplitudes? Or as we are at high tide??

% D4. 

clear; close all; clc;

% Solution to Part D of Project 2.
% In order to answer this question we must determine two parameters based
% on observations of the Gironde Estuary:
% i. B_0. This is the width at the entrance of the estuary (defined as a
% relatively narrow constriction point before reaching the open sea).
% ii. Briver. This is the width of the river where the estuary stops
% narrowing. It is supposed to represent the location where the width is
% determined by the river and not by the tides. Unfortunately in our case
% it seems like the estuary splits in two before this happens. 


% A2 avg = 1.5481
% A1 avg = 0.1079
% The yearly averaged semi-diurnal tides are approximately 15x the value of
% the diurnal tides. Therefore, we can use the shortcut mentioned in the
% project description, i.e. we use only semi-diurnal tides in this
% analysis.

% Solving the shallow water equations in shallow basins and estuaries.
% Based on Friedrichs (2011), Chapter 3 in Contemporary Issues in Estuarine
% Physics.

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
xr = 110e3;              % Point in river at which estuary stops narrowing
                         % the estuary splits into two rivers. In the two 
                         % cases the width becomes constant after 107 and 
                         % 111 km respectively (see fig), so we pick 110 
                         % as an average.
Lbasin=1.5*xr;           % Length of the basin or estuary in meters
B0=5.8e3;                % Width of the basin in meters at seaward side.
                         % This is based on Google Earth observations of
                         % the mouth of the Gironde.
Briver=0.25e3;           % Width of the river where estuary stops narrowing
                         % 0.2 km width for Dordogne, 0.3 km for Garonne,
                         % so average taken. 
deltaX=Lbasin/80;        % spatial step in meters.
Lb = -xr/log(Briver/B0); % e-folding length scale    
M2amp=1.55;              % Amplitude of M2 tide at seaward side <=> This is 
                         % the water level at the mouth.
discharge=0;             % Constant river discharge at landward boundary. 
                         % Prescribed as zero for simplicity in the project
                         % description.
Cd=2.5e-3;               % Drag coefficient
x=0:deltaX:Lbasin;       % With deltax = c.1.3e3, we have 80 grid points.

H0 = [0.66*8.5 8.5 1.5*8.5]; % Water Height sensitivity analysis
% H0 = 8.0:0.1:8.5;        % Define the bed level and so the depth of the basin.
                         % we observe that 8.5 is a good measure for water
                         % depth, since we see a variation of 10 cm across
                         % 85 km, corresponding to 6% change.
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

Matlab2_D1 = figure
plot(x(2:end),ZM2);
hold on
%plot(x(2:end),ZM4);
%plot(x(2:end),ZM6);
hold off
title('Gironde Estuary: M2 Amplitude');
xlabel('Length [m]');
ylabel('M2 Amplitude [m]');
grid on;
legend('Shallow water','Equilibrium depth','Deep water');
saveas(gcf,'Matlab2_D1_i.png');

Matlab2_D2 = figure
plot(x(2:end),phaseZM2);
grid on;
saveas(gcf,'Matlab2_D1_ii.png');

% figure
% yyaxis left;
% plot(x(2:end)/1000,ZM2);
% ylabel('SSE [m]');
% hold on
% yyaxis right;
% % plot(x(2:end)/1000,ZM4);
% hold off
% title('Gironde Estuary: M2 and M4 Tides');
% xlabel('x [km]');
% ylabel('SSE [m]');
% % legend('M2 a','b','c','d','e','M4 a','b','c','d','e');
% grid on;

% D1i. Shallow waters: equilibrium is not reached as friction is dominant
% over convergence and dampens the amplitudes across the basin. Equilibrium
% figure: the two factors are in balance such that the amplitude is
% constant up to 75% of the length of xr. Deep waters: equilibrium is not
% reached either because the tidal convergence dominates over friction.

% D1ii.  After xr, the width is constant, and river dynamics play a
% dominant role. Hence, we expect higher water levels due to water
% displacement. [!!!!] This is due to resonance (look at M4, M6) - this is 
% explained in previous parts as you look at length x-xr.

% D1iii. 

%Theoretical calculation of speed of propagation 
g = 9.81;
c_theory = sqrt(g*H0);

k = (coefs(:,1))';
c = wn(1)./abs(k);
%k = [1.934e-5 1.418e-5 9.849e-6];
%c = wn(1)./k;
d = c-c_theory; %diversion from theoretical value

display(c);
display(c_theory);
display(d);

% C = 7.2680    9.9128   14.2718
% C_th =  7.4185    9.1315   11.1838
% The propagation speed according to -kx is similar to that which was 
% calculated theoretically from the sqrt(gh). It seems like it is slower
% than the theoretical value, possibly due to friction effects. On the
% other hand, it is faster than the the theoretical value for the
% equilbrium value and deep value. In the former case, it is just slightly
% faster (approx 10%), while for the latter, it is significantly faster
% (approx 30%). This seems to suggest that the channel convergence has a
% quickening effect on the velocity.
% Compare to paper's function for c and explain the crossing by looking at
% Lb and La. 


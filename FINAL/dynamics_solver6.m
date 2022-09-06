% re-entry trajectory solver
% nadja radovic
clc
clear all

%%
% next steps:
% extend orbit apogee
% code in AoA - think about which angle it's trimmed/stable at
% think about how Cl Cd vary with AoA
% when and where (at which altitude) to angle down the re-entry path- i.e.
% -- make gamma non-zero

% input historical values for AoA, mass, Cd, Cl, gamma, V0, h0 etc
% see what happens 
% plot horizontal lines to show which part of atmosphere you're in

%% solver
tic
% define global parameters
global RE C_D C_L A m omega_E

RE      = 6371e3;
%C_D     = 1.5;
%C_L     = 0.3;
alpha = deg2rad(-12.6);
C_L = 0.2477;
C_D = 1.436;
%C_L = 0.275;
%C_D = 1.414;
A       = pi*(4.45/2)^2; 
m       = 7000;
omega_E = 7.2921159e-5;
%omega_E = 0;
%mu = 3.986004418e14;

%{
r       = y(1);
lambda  = y(2);
delta   = y(3);
v       = y(4);
chi     = y(5);
gamma   = y(6);
range   = y(7);
h       = y(1) - RE;
%}

% initial conditions

% Rotating earth frame
h0       = 120e3; 
r0       = h0 + RE; % m, initial altitude
lambda0  = deg2rad(0); % Latitude
delta0   = deg2rad(0); % Longitude
range0   = 0;

% Body orientation frame
%
V0      = 1000; % m/s, initial velocity
chi0    = deg2rad(15); % Roll angle?
gamma0  = deg2rad(-9); % Flight path angle
%}

%{
V0      = 1.098961773126173e+04; % m/s, initial velocity
chi0    = deg2rad(90-5.145); % Roll angle?
%gamma0  = deg2rad(-3.1781); % Flight path angle -other code input
gamma0  = deg2rad(-6); % Flight path angle -very funky
%gamma0  = deg2rad(-8); % Flight path angle
%}

% assign event for terminating integration (when reached altitude of zero)
opts=odeset('Events',@events);

tstep = [0 5e3];

% solve the differential equations
% note, ode45 gets funky sometimes with the timesteps, if the code isn't
% running all the way until re-entry, try increasing the '3000' to a
% greater t_final
%



% Constants
global G_const M_earth mu R_Earth
G_const = 6.673*10^-11; % Gravitational Constant in SI
M_earth = 5.972*10^24; % Mass of Earth in kg
mu = G_const*M_earth; % Gravitational parameter of Earth
R_Earth = 6371e3; % Radius of Earth in m

% Inputs
alt_max = 382260e3; % Lunar orbit
alt_VPA = 30.8e3; % Obtain nominal value from Design Envelope code
arg_p   = 0;
inc     = 5.145; % Lunar orbit inclination
asc_n   = 0;

% Atmospheric Entry Conditions obtained from orbital mechanics
[r_ent0, long_ent0, lat_ent0, V_ent0, azi_ent0, fpa_ent0, h_ent0,...
 delV, r_range, nu] = atmosEntry(alt_max, alt_VPA, arg_p, inc, asc_n);

[t, y] = ode45( @dynamics6, tstep, [r_ent0, long_ent0, lat_ent0, V_ent0, azi_ent0, fpa_ent0, range0, h_ent0], opts);

%[t, y] = ode45( @dynamics6, tstep, [6491000,-3.58596584221759,0.0897971900151083,10989.5989812087,1.48099913681576,-0.123291695753393,0,120000], opts);



%y(:,2) = rem(y(:,2),2*pi);
%y(:,3) = rem(y(:,3),2*pi);


%y(:,9) = y(:,4)./sqrt(1.4*287*)

%{
%check latitude and longitude angles
for n = 1:length(y)
    if y(n,2) > deg2rad(90)
        y(n,2) = deg2rad(180) - y(n,2);
    elseif y(n,2) < deg2rad(-90)
        holding = abs(deg2rad(180)-abs(y(n,2)));
        y(n,2) = -holding;
    end

    if y(n,3) > deg2rad(180)
        y(n,3) = y(n,3) - deg2rad(360);
    elseif y(n,3) < deg2rad(-180)
        y(n,3) = y(n,3) + deg2rad(360);
    end
end
%}


%% Convective Heat Transfer Rate
[T,P,rho]=atmosphere(y(:,8));

%hr(1:length(t)) = 1.1*ones(1,length(t));
for i=1:length(t)

if y(i,4)/(heat_ratio(T(i,1))*287*T(i,1))^0.5 <= 0.4
    %break;
end
    
% Engineering Method

% Determine convective heat flux using Sutton-Graves
R_eff=5.6162;   % Effective radius of heat shield (approx)
q_conv(i,1)=1.74e-4*y(i,4)^3*(rho(i)/R_eff)^0.5;  % Convective heat flux


% Determine peak wall temperature T_w_peak ASSUMING RADIATIVE EQUILIBRIUM,
% NO CONDUCTION!
T_w(i,1)=(q_conv(i,1)/(0.9*5.67e-8))^0.25;


% Boundary Layer Method

hr(i)=heat_ratio(T(i,1));  % Ratio of specific heats
hr_Cp=hr(i)*287/(hr(i)-1);
M(i,1)=y(i,4)/(hr(i)*287*T(i,1))^0.5;  % Mach number


% Calculate post shock flow properties (use Euler equations from datasheet)
% ASSUMES PERFECT GAS

    
% Note: At beginning of trajectory T2=29000K, how to model
% dissociation+ionisation?


p2      = P(i,1)*(1+2*hr(i)*(M(i,1)^2-1)/(hr(i)+1));
rho2    = rho(i,1)*(hr(i)+1)*M(i,1)^2/((hr(i)-1)*M(i,1)^2+2);
T2      = p2/(rho2*287);
M2      = (1+0.5*((i)-1)*M(i,1)^2)/(hr(i)*M(i,1)^2-0.5*(hr(i)-1));
q       = 0.5*rho(i,1)*y(i,4)^2;

% Stagnation values
T02     = T2*(1+0.5*(hr(i)-1)*M2^2);

% Calculate pressure coefficient for each panel
p02     = P(i,1)*((hr(i)+1)^2*M(i,1)^2/(4*hr(i)*M(i,1)^2-2*(hr(i)-1)))^(hr(i)/(hr(i)-1))*((1-hr(i)+2*hr(i)*M(i,1)^2)/(hr(i)+1));
Cpmax(i,1) = (p02-P(i,1))/q;

[theta,n]  = cone();

for j=1:n    
    Cp(i,j) = Cpmax(i,1)*sin(theta(j))^2;
    Cp_upper(i,j)=Cpmax(i,1)*(sin(theta(j)+abs(alpha)))^2;
end

% Calculate flow properties at edge of BL (assuming thin BL)
for j = 1:n
    q_nadja_upper(i,j) = real(q*Cp_upper(i,j));%+P(i,1);
    Pe(i,j) = q*Cp(i,j)+P(i,1);
    Te(i,j) = T02*((Pe(i,j))/(q*Cpmax(i,1)+P(i,1)))^((hr(i)-1)/hr(i));
    rhoe(i,j) = Pe(i,j)/(287*Te(i,j));
    Ve(i,j) = (2*hr_Cp*(T02-Te(i,j)))^0.5;
    Me(i,j) = Ve(i,j)/(hr(i)*287*Te(i,j))^0.5;
end

% Use reference temperature method to model BL properties

% Test streamline function

end
toc
tic
%% Plotting
%
% earth parameters
theta_earth = deg2rad(-180:1: 180);
r_earth = (RE)*ones(1,length(theta_earth));

%{
figure(1) % polar plot of earth + trajectory
polarplot(theta_earth,r_earth)
hold on
polarplot(y(:,3),y(:,1))
legend('Earth','Trajectory')
hold off
%}

figure(1)
%title('3D orbit in km scale')

[X1,Y1,Z1] = sphere;

X2 = X1 * RE;
Y2 = Y1 * RE;
Z2 = Z1 * RE;

surf(X2,Y2,Z2)
hold off


a = y(:,1).*sin(pi/2 - rem(y(:,3),2*pi)).*cos(rem(y(:,2),2*pi));
b = y(:,1).*sin(pi/2 - rem(y(:,3),2*pi)).*sin(rem(y(:,2),2*pi));
c = y(:,1).*cos(pi/2 - rem(y(:,3),2*pi));
xyz_prime = [a,b,c]';
for i = 2:length(xyz_prime)
    hold on
    %fnplt(cscvn(xyz_prime(:,[i-1:i i])),'m',2)
    plot3(a,b,c,'m-')
    %if i == length(xyz_prime)
     %   fnplt(cscvn(xyz_prime(:,[[end,1] i])),'m',2)
    %end
end 
axis equal


% recalculating parameters to allow for plotting
q_load = 0.5*rho.*(y(:,4)).^2;
D = q_load*A*C_D;
g = (9.81*(RE./(RE+y(:,8))).^2); % m/s^2, acceleration due to gravity
a_n = (-D./m - g.*sin(y(:,6)))/9.81; % acceleration in gs in normal direction 

figure(2) % acceleration vs time
plot(t, -a_n)
xlabel('Time (s)')
ylabel('Deceleration in gs (m/s^2)')
grid on

figure(3) % altitude vs time
plot(t, y(:,8)./1000)
xlabel('Time (s)')
ylabel('Altitude (km)')
ylim([0 max(y(:,8)./1000)+10])
grid on

figure(4) % velocity vs time
plot(t, y(:,4))
xlabel('Time (s)')
ylabel('Velocity (m/s)')
grid on

figure(5) % altitude vs velocity
plot(y(:,4), y(:,8)./1000)
xlabel('Velocity (m/s)')
ylabel('Altitude (km)')
ylim([0 max(y(:,8))/1000])
grid on

figure(6)
plot3(rad2deg(y(:,2)),rad2deg(y(:,3)),y(:,8)./1000)
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
zlabel('Altitude (km)')
zlim([0 max(y(:,8)/1000)])
grid on

figure(7)
plot3(t,y(:,7)./1000,y(:,8)./1000)
ylabel('Range (km)')
xlabel('time (s)')
zlabel('Altitude (km)')
zlim([0 max(y(:,8))/1000])
grid on

figure(8)
plot(rad2deg(y(:,2)),rad2deg(y(:,3)))
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
grid on

figure(9)
plot(y(:,7)./1000,y(:,8)./1000)
xlabel('Range (km)')
ylabel('Altitude (km)')
ylim([0 max(y(:,8)./1000)+10])
grid on

figure(10)

[~,index] = min(Cpmax);

plot(M(1:index),Cpmax(1:index))
xlabel('Mach Number [-]')
ylabel('Cp_{max} [-]')
ylim([floor(min(Cpmax(1:index))) ceil(max(Cpmax(1:index)))])
grid on

figure(11)
plot(T,y(:,8)./1000,'x')
xlabel('Atmosphere Temperature (K)')
ylabel('Altitude (km)')
ylim([ceil(min(y(:,8))/1000) max(y(:,8))/1000])
grid on

figure(12)
plot(T_w,y(1:length(T_w),8)./1000)
xlabel('Wall Temperature (K)')
ylabel('Altitude (km)')
%ylim([y(1,8)/1000 max(y(:,8))/1000])
grid on

%
figure(13)
plot(q_conv/(100)^2,y(1:length(q_conv),8)./1000)
xlabel('Heat Flux (W/cm^2)')
ylabel('Altitude (km)')
%ylim([y(length(q_conv),8)/1000 y(1,8)/1000])
grid on
%}

figure(14)
plot(rad2deg(y(:,6)),y(:,8)./1000)
xlabel('Flight Path Angle (degrees)')
ylabel('Altitude (km)')
grid on

figure(15)
%plot(q_nadja_upper/1000,y(:,8)./1000)
plot(q_nadja_upper(1:(0.9*length(y)),1)/1000,y(1:(0.9*length(y)),8)./1000)
xlabel('Dynamic Load (kPa)')
ylabel('Altitude (km)')
ylim([20 120])
grid on
%}
toc

%% Functions

function [value,terminate,direction] = events(t,y)
% check altitude  = 0 and stop
value = y(8);
terminate = 0;
if value <= 0
    terminate = 1; % stop the integration
end
direction = 0; % all events
end

% Cone discretisation
function [theta,n]=cone()
theta = [1.488324054880543,1.371729588672399,1.311009043609386,1.262456550274142,1.220249371580128,1.048213346583875,0.545455742565992,0.172956295308091];
%L = 0.609243337650500;
n=8;
end
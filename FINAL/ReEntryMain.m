%% Simplified 6 d.o.f. code accounting for spatial coordinate, and body trajectory
% Tanmay Ubgade 210608

%% housekeeping

clear
clc
%clf(fig1);

%% Overview

% All calculations are in SI units

% Run Orbital Mechanics calculations to obtain delta-V requirements,
% trajectory and velocity profile until Re-Entry 
% --> Moon orbit + max DeltaV point

% Obtain flight path angle, azimuth angle, local latitude & longitude,
% velocity and altitude 

% Compute entry trajectory until skip out altitude (find by trial and
% error?) 
% --> Design Envelope Code

% Compute Descent stage dynamics with atmospheric entry values 

% Compute Landing stage accounting for parachute deployment
% --> Parachute sizing code and CD change from Nadja

% Concatenate data sets of the above stages together to compile for overall
% plots

%% Entry - Orbital Mechanics calculations

% Design Delta-V requiremenet from maximum Delta-V point

% Constants
global G_const M_earth mu R_Earth
G_const = 6.673*10^-11; % Gravitational Constant in SI
M_earth = 5.972*10^24; % Mass of Earth in kg
mu = G_const*M_earth; % Gravitational parameter of Earth
R_Earth = 6371e3; % Radius of Earth in m

% Inputs
alt_max = 382260e3; % Lunar orbit
alt_VPA = 60e3; % Obtain nominal value from Design Envelope code
arg_p   = 0;
inc     = 5.145; % Lunar orbit inclination
asc_n   = 0;

% Atmospheric Entry Conditions obtained from orbital mechanics
[r_ent0, long_ent0, lat_ent0, V_ent0, azi_ent0, fpa_ent0, h_ent0,...
 delV, r_range, nu] = atmosEntry(alt_max, alt_VPA, arg_p, inc, asc_n);

% Orbit plotting (Comment out if not plotting return orbit)
%[xyz_rotated] = orbitRot(arg_p, inc, asc_n, r_range, nu); 

%% EDL Overall test case
%{
% define global parameters
global R_Earth C_D C_L A m omega_E

RE      = 6371e3;
%C_D     = 1.5;
%C_L     = 0.3;
AoA = -12.6;
C_D = 1.436;
C_L = 0.2477;
A       = pi*(4.45/2)^2; 
m       = 12383;
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
%{
V0      = 10600; % m/s, initial velocity
chi0    = deg2rad(15); % Roll angle?
gamma0  = deg2rad(-14.7); % Flight path angle
%}


%
V0      = 7900; % m/s, initial velocity
chi0    = deg2rad(75); % Roll angle?
gamma0  = deg2rad(-4.2); % Flight path angle
%gamma0  = deg2rad(-15); % Flight path angle
%

% assign event for terminating integration (when reached altitude of zero)
opts=odeset('Events',@events);

tstep = [0:1:10e5];

% solve the differential equations
% note, ode45 gets funky sometimes with the timesteps, if the code isn't
% running all the way until re-entry, try increasing the '3000' to a
% greater t_final
[t, y] = ode45( @ReEntryDynamics, tstep, [r0 lambda0 delta0 V0 chi0 gamma0 range0 h0], opts);

%y(:,2) = rem(y(:,2),2*pi);
%y(:,3) = rem(y(:,3),2*pi);

%}

%% Descent constants

% Obtain CL CD A and m from vehicle geometry function (to be made)

global C_D C_L A m omega_E
AoA = -12.6;
C_L = 0.2477;
C_D = 1.436;
A   = pi*(4.45/2)^2; 
m   = 7600;
omega_E = 7.2921159e-5;
range0 = 0;

tstep = [0 1e6]; % Time step

%% Descent Part 1 (Entry till drogue deployment)

% Inputs for entry - atmospheric entry taken from Orbit Mechanics

% ODE stop event - Drouge deployment speed of 137.2  m/s
opts_drog = odeset('Events',@drogue);

% Solve dynamics
[t_drog, y_drog] = ode45(@ReEntryDynamics, tstep, [r_ent0, long_ent0,...
                 lat_ent0, V_ent0, azi_ent0, fpa_ent0, range0, h_ent0], opts_drog);

%% Descent Part 2 (Drogue deployment till main deployment)
% ODE stop event - Main deployment speed of 58 m/s
opts_main = odeset('Events',@main);

% New inputs for drogue deployment
[r_drog, long_drog, lat_drog, V_drog, azi_drog, fpa_drog, range_drog, h_drog] = y_drog(end,:);

C_L = 0;
C_D = 0.6;
A   = 100; 

[t_main, y_main] = ode45(@ReEntryDynamics, tstep, [r_drog, long_drog,...
                   lat_drog, V_drog, azi_drog, fpa_drog, range_drog, h_drog], opts_main);


%% Landing (Main deployment to land)

% ODE stop event - Landed
opts_land = odeset('Events',@land);

% New inputs for main deployment
[r_main, long_main, lat_main, V_main, azi_main, fpa_main, range_main, h_main] = y_main(end,:);

C_L = 0;
C_D = 0.96;
A   = 2160;

[t_land, y_land] = ode45(@ReEntryDynamics, tstep, [r_main, long_main,...
                   lat_main, V_main, azi_main, fpa_main, range_main, h_main], opts_land);

%% Compiling data

% Overall time
t_main = t_main + t_drog(end);
t_land = t_land + t_main(end);
t_overall = [t_drog ; t_main ; t_land];

% Overall trajectory
y_overall = [y_drog ; y_main ; y_land];

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


%% Plotting


%% Local functions

function [value,terminate,direction] = drogue(t,y)
% Check 
value = y(4);
terminate = 0;
if value <= 140
    terminate = 1; % stop the integration
end
direction = 0; % all events
end

function [value,terminate,direction] = main(t,y)
% Check 
value = y(4);
terminate = 0;
if value <= 60
    terminate = 1; % stop the integration
end
direction = 0; % all events
end

function [value,terminate,direction] = land(t,y)
% Check 
value = y(8);
terminate = 0;
if value <= 0
    terminate = 1; % stop the integration
end
direction = 0; % all events
end
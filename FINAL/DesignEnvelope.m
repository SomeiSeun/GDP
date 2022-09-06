%% Code to acquire limits for Lifeboat entry corridor with a Design Envelope
% Tanmay Ubgade 210608

%% housekeeping

clear
clc
close all

%% Design envelope loop

% Initialise matrices
%success = zeros(1,101); skip_fail = success; acc_fail  = success;
j = 1;
startVPA = 35e3;

for i = 0:(230/5)
    tic
    %% Entry stage - Orbital Mechanics calculations

    % Design Delta-V requiremenet from maximum Delta-V point

    % Constants
    global G_const M_earth mu R_Earth
    G_const = 6.673*10^-11; % Gravitational Constant in SI
    M_earth = 5.972*10^24; % Mass of Earth in kg
    mu = G_const*M_earth; % Gravitational parameter of Earth
    R_Earth = 6371e3; % Radius of Earth in m

    % Inputs
    alt_max = 382260e3; % Lunar orbit
    alt_VPA = startVPA - 5*i*10^2; % Obtain nominal value from Design Envelope code
    arg_p   = 0;
    inc     = 5.145; % Lunar orbit inclination
    asc_n   = 0;

    % Atmospheric Entry Conditions obtained from orbital mechanics
    [r_ent0, long_ent0, lat_ent0, V_ent0, azi_ent0, fpa_ent0, h_ent0,...
     delV, ~, ~] = atmosEntry(alt_max, alt_VPA, arg_p, inc, asc_n);

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

    %tstep = 0:1:1e4; % Time step
    tstep = [0 1e4];

    %% Descent-Land stage

    % Inputs for entry - atmospheric entry taken from Orbit Mechanics

    % ODE stop event - Drouge deployment speed of 137.2  m/s
    opts = odeset('Events',@events);

    % Solve dynamics
    [t, y] = ode45(@ReEntryDynamics, tstep, [r_ent0, long_ent0,...
                     lat_ent0, V_ent0, azi_ent0, fpa_ent0, range0, h_ent0], opts);              

    %% Dynamic load and deceleration
    % recalculating parameters to allow for plotting
    [~,~,rho]=atmosphere(y(:,8));
    q_load = 0.5*rho.*(y(:,4)).^2;
    D = q_load*A*C_D;
    g = gravity(y(:,8)); % m/s^2, acceleration due to gravity
    acc = (-D./m - g.*sin(y(:,6)))/9.81; % acceleration in gs in normal direction 
    
    %% What kind of failure?

    % Checking for skip out
    check_skip  = y(:,8) > 200e3;
    if max(check_skip) == logical(1)
        skip_fail(i+1) = 1;
    else
        skip_fail(i+1) = 0;
    end

    % Checking for too high deceleration
    check_acc   = abs(acc(:)) > 10;
    if max(check_acc) == logical(1)
        acc_fail(i+1) = 1;
    else
        acc_fail(i+1) = 0;
    end

    if skip_fail(i+1) == 1 || acc_fail(i+1) == 1
        success(i+1) = 0;
    else
        success(i+1) = 1;
    end

    %% Storing important parameters
    % VPA
    VPA(i+1) = (startVPA - 5*i*10^2)./1000;
    if success(i+1) == 1
        good_start(j,:) = [r_ent0, long_ent0, lat_ent0, V_ent0, azi_ent0, fpa_ent0, range0, h_ent0, VPA(i+1), t(end)];
        j = j+1;
    end
    toc
    
    clear r_ent0 long_ent0 lat_ent0 V_ent0 azi_ent0 fpa_ent0 h_ent0 delV alt_VPA t y acc D q_load rho
end

%% Which ones work

figure(1)
hold on
title('Successful VPAs')
plot(VPA, success,'bx')
xlabel('VPA (km)')
ylabel('Success')
xlim([VPA(end) VPA(1)])
ylim([-0.2 1.2])
grid on
hold off

figure(2)
hold on
title('Max G-loading exceeded')
plot(VPA, acc_fail,'rx')
xlabel('VPA (km)')
ylabel('G failure')
xlim([VPA(end) VPA(1)])
ylim([-0.2 1.2])
grid on
hold off

figure(3)
hold on
title('Skip outs occur')
plot(VPA, skip_fail,'rx')
xlabel('VPA (km)')
ylabel('Skip out')
xlim([VPA(end) VPA(1)])
ylim([-0.2 1.2])
grid on
hold off

figure(4)
hold on
title('Entry Corridor')
plot(rad2deg(good_start(:,6)),good_start(:,9),'k-')
plot([rad2deg(good_start(1,6)) rad2deg(good_start(1,6))],[0 100],'b--')
plot([rad2deg(good_start(end,6)) rad2deg(good_start(end,6))],[0 100],'r--')
legend('Entry Parameters','Skip-out limit','G-loading limit','Location','Southeast')
xlabel('Entry Flight Path Angle (degrees)')
ylabel('Vaccuum Periapsis Altitude (km)')
ylim([good_start(end,9) good_start(1,9)])

%{
figure
plot3(skip_fail,acc_fail, (130e3 - 0.5.*(0:300).*10^3)./1000)
xlabel('Skip Fail')
ylabel('Deceleration Fail')
zlabel('Altitude VPA (km)')
grid on
%}
%% Local Functions

function [value,terminate,direction] = events(t,y)
% check altitude  = 0 and stop
value = y(8);
terminate = 0;
if value <= 0
    terminate = 1; % stop the integration
end
direction = 0; % all events
end

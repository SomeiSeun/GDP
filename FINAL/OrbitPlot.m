%% Code to plot 3D Orbits
% Tanmay Ubgade 210611

%% housekeeping
clear all
close all
clc

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
arg_p   = 240; %L4 point
inc     = 5.145; % Lunar orbit inclination
asc_n   = 0;

% Atmospheric Entry Conditions obtained from orbital mechanics
[r_ent0, long_ent0, lat_ent0, V_ent0, azi_ent0, fpa_ent0, h_ent0,...
 delV, r_range, nu] = atmosEntry(alt_max, alt_VPA, arg_p, inc, asc_n);

% Circular
[~, ~, ~, ~, ~, ~, ~,~, r_circ, ~] = atmosEntry(alt_max, alt_max, arg_p, inc, asc_n);

% Orbit plotting (Comment out if not plotting return orbit)
[xyz_rotated] = orbitRot(arg_p, inc, asc_n, r_range, nu); 
[xyz_rotated_circ] = orbitRot(arg_p, inc, asc_n, r_circ, nu); 
[moon]        = orbitRot(0, inc, asc_n, 382260e3, 0); 

xyz_rotated = xyz_rotated./1000; 
xyz_rotated_circ = xyz_rotated_circ./1000; 
%% Plotting

figure(1)
title('3D orbit in km scale','FontSize',24)
hold on
planet3D('Earth Cloudy')
planet3D('Moon',moon./1000)

plot3(xyz_rotated_circ(1,:),xyz_rotated_circ(2,:),xyz_rotated_circ(3,:),'b-')
plot3(xyz_rotated(1,:),xyz_rotated(2,:),xyz_rotated(3,:),'m-')
axis equal
xlim([-10000-R_Earth/1000 390000])
hold off

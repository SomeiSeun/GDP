function [xyz] = earth2xyz(r,long,lat)
% Function to convert radius, longitude and latitude to xyz for plotting
% trajectory.
% Radius in metres, Longitude and Latitude in rads. All as same sized
% arrays linked to same timestep.

x = r.*sin(pi/2 - rem(lat,2*pi)).*cos(rem(long,2*pi));
y = r.*sin(pi/2 - rem(lat,2*pi)).*sin(rem(long,2*pi));
z = r.*cos(pi/2 - rem(lat,2*pi));

xyz = [x,y,z]';
end
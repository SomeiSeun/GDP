% test
clear all
clc
close all


%% Plotting Atmosphere Model
%{
h = 0:100:120e3;
[T,~,rho] = atmosphere(h);

figure(1)
hold on
title('Density Model','FontSize',16)
plot(rho,h./1000);
xlabel('Density (g/cm^3)','FontSize',14)
ylabel('Altitude (km)','FontSize',14)
box on
grid on
hold off

figure(2)
hold on
title('Temperature Model','FontSize',16)
plot(T,h./1000);
xlabel('Temperature (K)','FontSize',14)
ylabel('Altitude (km)','FontSize',14)
box on
grid on
hold off

h = 0:100:500e3;
[T_ext,~,~] = atmosphere(h);

figure(3)
hold on
title('Extended Temperature Model','FontSize',16)
plot(T_ext,h./1000);
xlabel('Temperature (K)','FontSize',14)
ylabel('Altitude (km)','FontSize',14)
box on
grid on
hold off
%}
%% Earth

figure(1)
title('3D orbit in km scale','FontSize',24)
hold on
planet3D('Moon')
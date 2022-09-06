%% Finding the Optimal deltaV value from re-entry perspective
% Assuming no atmospheric effects as it gives upper constraint of delV
% requirement
% Tanmay Ubgade 210517

%% housekeeping
clear
clc
close all

%% Constants
G_const = 6.673*10^-11; % Gravitational Constant in SI
M_earth = 5.972*10^24; % Mass of Earth in kg
mu = G_const*M_earth/((10^3)^3); % Gravitational parameter of Earth
r_earth = 6371; % Radius of Earth in km

%% Inputs

% Apoapsis
alt_max = 300:100:500000; 
% 384400
r_a = r_earth + alt_max;
r_p = r_a;

%% Calculations

a= 0; % initialise
for i = 1:length(r_a)
    a(i)    = (r_a(i) + r_p(i))/2; % semi major axis
    V1(i)   = VisViva(r_a(i),a(i),mu); % Before burn
    V2(i)   = VisViva(r_a(i),(r_a(i)+r_earth+100)/2,mu); % After burn
    delV(i) = V1(i)-V2(i);
    Vp(i)   = VisViva(r_earth+100,(r_a(i)+r_earth+100)/2,mu); % Velocity at atmospheric entry
end
delV = delV*1000;
[delV_min, index] = min(delV);
Vp_min = Vp(index)
V2_min = V2(index)
V_entry = 0;
[delV_max, index_max] = max(delV);

%% Plotting
figure(1)
hold on
plot(alt_max,delV,'b-')
plot([0 500000],[2*delV(3841) 2*delV(3841)],'m-')
plot(300+index_max*100,delV(index_max),'rx','MarkerSize',30)
plot(384400,delV(3841),'kx','MarkerSize',30)
grid minor
xlabel('Altitude (km)','FontSize',24)
ylabel('Delta V (m/s)','FontSize',24)
title('Delta V Variation with Altitude for Re-entry','FontSize',24)
xlim([300 500000])
ylim([0 2*delV(3841)+100])
title('Delta V variation with altitude')
legend('Delta V variation','Saftey Factor Fuel','Max Delta V','Moon Lagrange Point','FontSize',22,'Location','Southeast')
box on
hold off

%% Functions 

function [v_sqrt] = VisViva(r,a,mu)

v_sqrt = sqrt(mu.*((2./r) - (1./a)));

end

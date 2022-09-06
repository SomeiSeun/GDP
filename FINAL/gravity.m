function [g] = gravity(h)
% This function receives in input the altitude H and provides as output the
% corresponding gravitational acceleration g for earth
%   g: gravitational acceleration at altitude h (m/s^2)
%   h: altitude (m)

g0 = 9.80665; % Standard gravitational acceleration of earth at sea leverl (m/s^2)
R_E = 6.371e6; % Mean earth radius (m)

g = g0*(R_E./(R_E + h)).^2;

end
function [T] = Kepler3(a,mu)
% Outputs orbit time period using Kepler's 3rd law
T = 2*pi*sqrt((a^3)/(mu));

end
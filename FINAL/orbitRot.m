%% Function that takes orbit characteristics and outputs cartesian points for plot
% Tanmay Ubgade 210609

function [xyz_rotated] = orbitRot(arg_p, inc, asc_n, r_range, nu)
% Orbital rotations (deg), radius (m) and true anamoly matrix (rad)
x_transit = r_range.*cos(nu);
y_transit = r_range.*sin(nu);
z_transit = zeros(1,length(x_transit));

xyz = [x_transit;y_transit;z_transit]; % Base coordinates

syms arg incl an
R_arg       =  [cos(arg) -sin(arg) 0;...
                sin(arg) cos(arg)  0;...
                0           0      1]; % Rotation related to argument of periapsis
R_incl      =  [cos(incl)   0      -sin(incl);...
                0           1      0;...
                sin(incl)   0      cos(incl)]; % Rotation related to inclination of orbit
R_ascnode   =  [cos(an)  -sin(an)  0;...
                sin(an)   cos(an)  0;...
                0           0      1]; % Rotation related to ascending node
R_tot       = R_ascnode*R_incl*R_arg; % Full rotation matrix

% Angles for rotations
arg_p_rad = deg2rad(arg_p);
inc_rad   = deg2rad(inc);
asc_n_rad = deg2rad(asc_n);

% Rotated coordinates
for i = 1:length(x_transit)
    xyz_rotated(:,i) = double(subs(R_tot*xyz(:,i),[arg, incl, an],[arg_p_rad, inc_rad, asc_n_rad]));
end

end
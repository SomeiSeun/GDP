%% Function to provide atmospheric entry parameters using orbital mechanics
% Tanmay Ubgade 210609

function [r_entry, long_entry, lat_entry, V_entry, azi_entry, fpa_entry, h_entry, delV, r_range, nu] = atmosEntry(alt_max, alt_VPA, arg_p, inc, asc_n)

global mu R_Earth

% Initial orbit assumption - circular
r_a     = R_Earth + alt_max; % Also semi-major axis before burn
V1      = VisViva(r_a,r_a,mu);

% Post-burn
r_VPA   = R_Earth + alt_VPA;
V2      = VisViva(r_a,(r_a+r_VPA)/2,mu); % After burn

a       = (r_a + r_VPA)/2; % semi major axis
b       = sqrt(r_a * r_VPA); % semi minor axis
e       = sqrt(1 - ((b^2)/(a^2))); % eccentricity 

% Transit time to VPA
T_VPA   = Kepler3((r_a+r_VPA)/2,mu)/(2*60*60);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                Outputs                              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Delta V - De-orbit burn (m/s)
delV    = V1-V2;

% True anamoly array (rads)
nu      = linspace(0,360,2000).*pi/180; 
%nu     = [linspace(180,360,5000), linspace(0,180,5000)].*pi/180;

% Radius Range array (m)
r_range = a.*(1-e.^2)./(1+e.*cos(nu));

% Entry Altitude (m)
h_entry = 120e3;

% Entry Radius (m)
r_entry = R_Earth + h_entry;

% Entry Velocity (m/s)
V_entry = VisViva(r_entry,a,mu);

% Entry FPA (rads)
[~,index]   = min(abs(r_range-r_entry));
radius_1    = r_range(index);
radius_2    = r_range(index+1);
delta_nu    = nu(index)-nu(index+1);
if radius_1 > radius_2
    holder = radius_1;
    radius_1 = radius_2;
    radius_2 = holder;
    clear holder
    delta_nu    = nu(index+1)-nu(index);
end
delta_entry = sqrt(radius_1^2 + radius_2^2 -...
              2*radius_1*radius_2*cos(delta_nu));
fpa_entry   = pi/2 - acos((delta_entry^2 + radius_1^2 - radius_2^2)...
              /(2*delta_entry*radius_1));

          
          
% Entry Azimuthal Angle (rads)
xyz         = orbitRot(arg_p, inc, asc_n, r_range(1:2), nu(1:2));
azi_entry   = pi/2 - atan(max([(xyz(3,1)-xyz(3,2))/((xyz(1,1)-xyz(1,2))),...
                  (xyz(3,1)-xyz(3,2))/((xyz(2,1)-xyz(2,2))),...
                  (xyz(3,1)-xyz(3,2))/(sqrt(((xyz(1,1)-xyz(1,2)))^2 + ((xyz(2,1)-xyz(2,2)))^2))]));

% Entry Longitude (rads)
long_init   = rem(deg2rad(arg_p + asc_n),pi);
Earth_rot   = 23.93; % hours
long_change = (T_VPA/Earth_rot)*2*pi; 
long_entry  = rem(long_init-long_change,2*pi)-pi;

% Entry Latitude (rads)
lat_entry = deg2rad(inc); % bad approximation

end
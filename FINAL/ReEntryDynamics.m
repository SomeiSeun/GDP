%% Governing equations for re-entry model with rotating Earth
% Tanmay Ubgade 210603

function dy = ReEntryDynamics(t, y)
% this function integrates the dynamic system dy = y' = f(y, t)

% global parameters
global R_Earth C_D C_L A m omega_E

% Specify the inputs to the function
r       = y(1); % radius
lambda  = y(2); % longitude
delta   = y(3); % latitude
v       = y(4); % velocity
chi     = y(5); % azimuth/roll angle
gamma   = y(6); % flight path angle
range   = y(7); % range of travel
h       = y(1) - R_Earth; % altitude

%R_Earth = 6371e3; % m, earth radius

% calculating gravity
g = 9.81*(R_Earth/(R_Earth+h))^2; % m/s^2, acceleration due to gravity

% calculating aerodynamic forces
%rho = exp_density(h);
%[TEMP, RHO] = atmosnrlmsise00((r-RE),lambda,delta,2007,4,0);
%rho = RHO(6);
[~,~,rho] = atmosphere(h);
q = 0.5*rho*(v)^2;
D = q*A*C_D;
L = q*A*C_L;


% governing equations
dy(1) = v*sin(gamma);

dy(2) = (v*cos(gamma)*sin(chi))/(r*cos(delta));

dy(3) = (v*cos(gamma)*cos(chi))/r;

dy(4) = -g*sin(gamma) - D/m + ...
        (omega_E^2)*r*cos(delta)*(cos(delta)*sin(gamma) - ...
        sin(delta)*cos(gamma)*cos(chi));

dy(5) = v*cos(gamma)*sin(chi)*sin(delta)/(r*cos(delta))+ ...
        2*omega_E*(sin(delta)-tan(gamma)*cos(delta)*cos(chi)) + ...
        (omega_E^2)*r*cos(delta)*sin(chi)*sin(delta)/(v*cos(gamma));

dy(6) = -g*cos(gamma)/v + L/(m*v) + v*cos(gamma)/r + ...
        2*omega_E*sin(chi)*cos(delta) + ...
        (omega_E^2)*r*cos(delta)*(cos(delta)*cos(gamma) + ...
        sin(gamma)*sin(delta)*cos(chi))/v;

dy(7) = R_Earth*v*cos(gamma)/(R_Earth+h);
    
dy(8) = dy(1);

% return solution
dy = dy';
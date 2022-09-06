function gamma = heat_ratio(T)
% gamma - Specific Heat Ratio of air varying by temperature.
% T - Temperature of air
% Data taken from https://www.engineeringtoolbox.com/air-specific-heat-capacity-d_705.html
% Assumes pressure of 1 bar and extrapolates from Temp varying dataset 

temp_range  = [60, 78.79, 81.61, 100:20:260, 273.2, 280, 288.7,...
               300:20:400, 500:100:900, 1100:400:1900];
gamma_range = [1.621, 1.839, 1.452, 1.428, 1.415, 1.411, 1.410, 1.407,...
               1.406, 1.404, 1.404, 1.403, 1.403, 1.402, 1.402, 1.402,...
               1.400, 1.400, 1.398, 1.397, 1.396, 1.387, 1.375, 1.365,...
               1.354, 1.344, 1.329, 1.311, 1.301];

gamma       = interp1(temp_range, gamma_range, T,'linear','extrap');

end
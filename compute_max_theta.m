function [max_theta,max_alpha,max_beta] = compute_max_theta(LEO_altitude,UAV_altitude)
    R_E = 6371;  
    %% ----- Step 1: Define constants -----
          % Earth's radius in km
    z_ant = UAV_altitude/1000;          % Antenna height in km 
    r_A = R_E + z_ant;     % Distance from Earth center to Antenna
    
    h_sat = LEO_altitude;           % Satellite altitude in km
    r_sat = R_E + h_sat;   % Distance from Earth center to Satellite
    
    %% ----- Step 2: Geometric Calculation -----
    % alpha: Angle from antenna to tangent point
    max_alpha = acos(R_E / r_A); 
    
    % beta: Angle from tangent point to satellite
    max_beta = acos(R_E / r_sat);
    
    % Total Maximum Visibility Angle
    max_theta = max_alpha + max_beta;
end
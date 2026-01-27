%% Maximum Satellite Visibility Angle - Vertical Orientation
% Calculates the maximum central angle (theta_max) and visualizes 
% the geometry with the antenna at the top (Y-axis).

clc; clear; close all;

%% ----- Step 1: Define constants -----
          % Earth's radius in km
z_ant = 0.15;          % Antenna height in km (150 meters)
r_A = R_E + z_ant;     % Distance from Earth center to Antenna

h_sat = 500;           % Satellite altitude in km
r_sat = R_E + h_sat;   % Distance from Earth center to Satellite

%% ----- Step 2: Geometric Calculation -----
% alpha: Angle from antenna to tangent point
alpha = acos(R_E / r_A); 

% beta: Angle from tangent point to satellite
beta = acos(R_E / r_sat);

% Total Maximum Visibility Angle
theta_max_rad = alpha + beta;
theta_max_deg = rad2deg(theta_max_rad);

%% ----- Step 3: Visualization (Rotated) -----
figure('Color', 'w', 'Position', [100, 100, 800, 800]); 
hold on; axis equal; grid on;

% Offset for the antenna to be at the top (90 degrees)
offset = pi/2;

% Draw Earth (2D Slice)
t = linspace(0, 2*pi, 1000);
plot(R_E*cos(t), R_E*sin(t), 'k--', 'LineWidth', 1.2); 

% 1. Antenna Position (Fixed at the top)
A = [r_A * cos(offset); r_A * sin(offset)];
plot(A(1), A(2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
text(A(1), A(2) + 100, 'Antenna', 'Horiz','center', 'FontWeight', 'bold');

% 2. Satellite Position (at the visibility limit)
% We subtract theta_max to show it on the right side of the antenna
S_angle = offset - theta_max_rad;
S = [r_sat * cos(S_angle); r_sat * sin(S_angle)];
plot(S(1), S(2), 'bs', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
text(S(1) + 100, S(2), 'Satellite (Limit)', 'FontWeight', 'bold');

% 3. Tangent Point (where LOS touches Earth surface)
T_angle = offset - alpha;
T = [R_E * cos(T_angle); R_E * sin(T_angle)];
plot(T(1), T(2), 'gx', 'MarkerSize', 12, 'LineWidth', 2);
text(T(1) + 50, T(2) + 50, 'Tangent Point', 'Color', [0 0.5 0]);

% Draw Line of Sight (LOS)
plot([A(1) T(1) S(1)], [A(2) T(2) S(2)], 'm', 'LineWidth', 2.5);

% Draw Radial Lines for Clarity
plot([0 A(1)], [0 A(2)], 'k:', 'LineWidth', 1); % To Antenna
plot([0 S(1)], [0 S(2)], 'k:', 'LineWidth', 1); % To Satellite
plot([0 T(1)], [0 T(2)], 'g--', 'LineWidth', 1); % To Tangent Point (Right angle)

% Adjust view
title('Satellite Visibility Geometry (Antenna at Zenith)');
xlabel('km'); ylabel('km');
legend('Earth Surface', 'Antenna', 'Satellite', 'Tangent Point', 'Line of Sight');
xlim([-r_sat-200, r_sat+200]);
ylim([-R_E-500, r_sat+500]);
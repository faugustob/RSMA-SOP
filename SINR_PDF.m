%% Monte Carlo Verification: PDF Transformation of a Mixture Model
% Transformation Function: Y = (a*X) / (b*X + c)
% Input X: A mixture of Log-Normal and Gamma distributions
% Purpose: Validate the Change of Variables formula f_Y(y) = f_X(g^-1(y)) * |dx/dy|

clear; clc; close all;

%% 1. Parameters for the Transformation and Mixture
% Transformation coefficients
a = 2.0; 
b = 1.5; 
c = 3.0;

% Mixing proportion (alpha = weight of the first component)
alpha = 0.4; 

% Component 1: Log-Normal distribution parameters (mu and sigma in log-scale)
mu_ln = 0.5; 
sigma_ln = 0.5; 

% Component 2: Gamma distribution parameters (shape k and scale theta)
k_g = 2.0; 
theta_g = 1.0; 

% Number of Monte Carlo samples
N = 1e6; 

%% 2. Generate Random Samples from the Mixture Model
% We first determine which component each sample belongs to using a Bernoulli trial
component_idx = rand(N, 1) < alpha; % Logic array: true for Log-Normal, false for Gamma

x_samples = zeros(N, 1);

% Sample from Log-Normal for indices where component_idx is true
x_samples(component_idx) = lognrnd(mu_ln, sigma_ln, [sum(component_idx), 1]);

% Sample from Gamma for indices where component_idx is false
x_samples(~component_idx) = gamrnd(k_g, theta_g, [sum(~component_idx), 1]);

%% 3. Apply the Transformation to the Samples
% This creates the empirical data set for Y
y_samples = (a .* x_samples) ./ (b .* x_samples + c);

%% 4. Define the Theoretical PDF Expression
% The inverse function: x = g^-1(y) = (c*y) / (a - b*y)
inv_g = @(y) (c .* y) ./ (a - b .* y);

% The Jacobian: |dx/dy| = |(ac) / (a - b*y)^2|
jacobian = @(y) abs(a * c) ./ (a - b .* y).^2;

% Individual PDF of the Log-Normal component
ln_pdf = @(x) (1./(x .* sigma_ln .* sqrt(2*pi))) .* exp(-(log(x) - mu_ln).^2 ./ (2*sigma_ln^2));

% Individual PDF of the Gamma component
g_pdf = @(x) (x.^(k_g-1) .* exp(-x ./ theta_g)) ./ (gamma(k_g) * theta_g^k_g);

% Combined Mixture PDF for the input variable X
mix_pdf_X = @(x) alpha * ln_pdf(x) + (1 - alpha) * g_pdf(x);

% The Final Theoretical PDF for Y using the Transformation Theorem
% f_Y(y) = f_X(inv_g(y)) * jacobian(y)
f_Y_theoretical = @(y) mix_pdf_X(inv_g(y)) .* jacobian(y);

%% 5. Visualization and Comparison
figure('Name', 'PDF Transformation Verification', 'Color', 'w', 'Units', 'normalized', 'Position', [0.2 0.2 0.6 0.5]);
hold on;

% Plot Normalized Histogram (Empirical Density)
% We use 'pdf' normalization so the total area under the bars equals 1
histogram(y_samples, 100, 'Normalization', 'pdf', ...
    'FaceColor', [0.4, 0.6, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.7);

% Define the range for the theoretical curve
% The support of Y is bounded by the horizontal asymptote at y = a/b
y_limit = a/b;
y_grid = linspace(0.001, y_limit - 0.01, 1000); % Avoid division by zero at the asymptote

% Plot the Theoretical PDF Curve
plot(y_grid, f_Y_theoretical(y_grid), 'r-', 'LineWidth', 2.5);

% Formatting the Plot
grid on;
set(gca, 'FontSize', 12);
xlabel('y = (aX)/(bX+c)', 'Interpreter', 'latex');
ylabel('Probability Density', 'Interpreter', 'latex');
title('Verification of PDF Transformation (Mixture Model)', 'FontSize', 14);
legend('Monte Carlo Samples (Empirical)', 'Theoretical Formula (Analytical)', 'Location', 'best');

% Print summary to command window
fprintf('Simulation complete.\n');
fprintf('Transformation: Y = (%.1f*X) / (%.1f*X + %.1f)\n', a, b, c);
fprintf('Theoretical Support of Y: [0, %.2f)\n', y_limit);

%% Analytical CDF calculation
y_grid = linspace(0.001, a/b - 0.01, 1000);
x_inv = (c .* y_grid) ./ (a - b .* y_grid);

% Component CDFs
F_ln = logncdf(x_inv, mu_ln, sigma_ln);
F_g  = gamcdf(x_inv, k_g, theta_g);
F_y_theoretical = alpha * F_ln + (1-alpha) * F_g;

% Plotting
figure;
[f_emp, y_emp] = ecdf(y_samples); % Empirical CDF from samples
plot(y_emp, f_emp, 'b', 'LineWidth', 2); hold on;
plot(y_grid, F_y_theoretical, 'r--', 'LineWidth', 2);
title('CDF Verification');
legend('Empirical CDF', 'Theoretical CDF');
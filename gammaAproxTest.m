%% === Gamma approximation test for the marginal PDF of Z elements ===

% Choose which element of Z to inspect (e.g., first component)
idx = 1;

% Resimulate one big batch (or reuse from loop above if you stored samples)
nPlot = 5e5;
U_a = copularnd('Gaussian', (K_a./var_a), nPlot);
U_b = copularnd('Gaussian', (K_b./var_b), nPlot);

theta1 = Omega1 / m1; theta2 = Omega2 / m2;
A2 = gaminv(U_a, m1, theta1);
B2 = gaminv(U_b, m2, theta2);
A = sqrt(A2);
B = sqrt(B2);
Z = A .* B;

% Take the column of interest
Zi = Z(:, idx);

% Empirical mean and variance of that element
mu_Zi  = mean(Zi);
var_Zi = var(Zi);

% Gamma parameters (moment matching)
alpha = mu_Zi^2 / var_Zi;       % shape
beta  = var_Zi / mu_Zi;         % scale

% Plot histogram and Gamma PDF overlay
figure;
hold on;
histogram(Zi, 'Normalization', 'pdf');
xgrid = linspace(min(Zi), max(Zi), 400);
plot(xgrid, gampdf(xgrid, alpha, beta), 'r', 'LineWidth', 1.8);
hold off;
xlabel('z value');
ylabel('PDF');
title(sprintf('Gamma approximation for Z(%d)\nshape=%.2f, scale=%.2f', idx, alpha, beta));
legend('Empirical histogram', 'Gamma fit');
grid on;

%% Parameters
k1 = 3;      % shape of first gamma
lambda1 = 0.5; % rate of first gamma  (theta1 = 1/lambda1 = 2)
k2 = 1;      % shape of second gamma
lambda2 = 100;   % rate of second gamma  (theta2 = 1/lambda2 = 1)

N = 1e6; % Monte Carlo sample size

%% Generate Gamma random variables
% Note: MATLAB's gamrnd uses scale = 1/lambda
X1 = gamrnd(k1, 1/lambda1, [N, 1]);
X2 = gamrnd(k2, 1/lambda2, [N, 1]);
S = X1 + X2;

%% Moment matching for the sum (Gamma approximation)
% Theoretical mean and variance of the sum
meanS = k1/lambda1 + k2/lambda2;
varS  = k1/(lambda1^2) + k2/(lambda2^2);

% Moment-matched parameters for the approximating Gamma
lambda_s = meanS / varS;      % rate parameter
k_s      = meanS * lambda_s;  % shape parameter

fprintf('Moment-matched parameters:\n');
fprintf('   k_s = %.4f\n', k_s);
fprintf('   lambda_s = %.4f\n', lambda_s);

%% Compare with Monte Carlo histogram
figure;
histogram(S, 'Normalization', 'pdf', 'EdgeColor', 'none');
hold on;

x = linspace(min(S), max(S), 1000);
approx_pdf = gampdf(x, k_s, 1/lambda_s); % MATLAB uses scale = 1/lambda

plot(x, approx_pdf, 'r', 'LineWidth', 2);
title('Sum of Two Gamma RVs — Moment Matching (Rate Form)');
legend('Monte Carlo Simulation', 'Moment-Matched Gamma Approximation');
xlabel('S'); ylabel('PDF');
grid on;

%% Optional: Compare CDFs
[f_emp, x_emp] = ecdf(S);
F_approx = gamcdf(x_emp, k_s, 1/lambda_s);

figure;
plot(x_emp, f_emp, 'b', 'LineWidth', 1.5); hold on;
plot(x_emp, F_approx, 'r--', 'LineWidth', 2);
legend('Empirical CDF', 'Approximated Gamma CDF');
xlabel('S'); ylabel('CDF');
title('CDF Comparison');
grid on;

% K–S statistic
ks_stat = max(abs(f_emp - F_approx));
fprintf('\nKolmogorov–Smirnov statistic: %.4f\n', ks_stat);

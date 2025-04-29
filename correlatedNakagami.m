% Optimized simulation of correlated Nakagami-m random variables
% for IRS-aided wireless system (based on Ajayan et al., IEEE Access, 2023)

%% 1. Initialization and Parameters
rng(42);                        % Fixed seed for reproducibility
N         = 25;                 % Number of IRS elements
m1        = 1;                  % Nakagami-m shape for h_i
m2        = 2;                  % Nakagami-m shape for g_i
Omega1    = 5;                  % Spread for h_i
Omega2    = 5;                  % Spread for g_i
kappa   = 0.4;                % Exponential correlation coefficient
N_samples = 1e6;                % Monte Carlo sample count

% Derived gamma distribution scales
theta1 = Omega1 / m1;           % scale for |h_i|^2
theta2 = Omega2 / m2;           % scale for |g_i|^2

%% 2. Exponential Correlation Matrix Construction
% Use vectorized indexing to build |i-j| matrix
i = (1:N)'; j = 1:N;
R = kappa .^ abs(i - j);      % Exponential correlation: R(i,j) = epsilon^|i-j|

% Ensure symmetry and unit diagonal
R = (R + R') / 2;               % force symmetry
R(1:N+1:end) = 1;               % set diagonal elements to 1

% Enforce positive definiteness if needed
minEig = min(eig(R));
if minEig <= 0
    warning('Adjusting correlation matrix to nearest SPD (min eig = %.3e)...', minEig);
    R = nearestSPD(R);
    % normalize diagonals
    D = sqrt(diag(R));
    R = R ./ (D * D');
end

%% 3. Generate Correlated Samples via Gaussian Copula
% Generate correlated uniforms without loops
U_h = copularnd('Gaussian', R, N_samples);
U_g = copularnd('Gaussian', R, N_samples);

% Transform uniforms to gamma variates (|h|^2 and |g|^2)
H2 = gaminv(U_h, m1, theta1);   % size: [N_samples x N]
G2 = gaminv(U_g, m2, theta2);

%% 4. Form Nakagami-m Variables and Composite Channel
H = sqrt(H2);                   % |h_i|
G = sqrt(G2);                   % |g_i|
X = sum(H .* G, 2) / sqrt(Omega1 * Omega2);

%% 5. Visualization: PDF of Composite Coefficient
figure;
histogram(X, 'Normalization', 'pdf', 'EdgeColor', 'none');
xlabel('Composite Channel Coefficient X');
ylabel('Probability Density');
title('Simulated PDF of Composite Channel Coefficient');
grid on;

%% 6. Covariance Check for H
H_cov_est = cov(H);             % estimated covariance matrix
disp('Estimated covariance (top-left 5x5 block):');
disp(H_cov_est(1:5, 1:5));

%% Subfunction: Nearest Symmetric Positive Definite Matrix
function A_spd = nearestSPD(A)
    % NearestSPD: find nearest symmetric positive definite matrix
    B   = (A + A') / 2;
    [U,S,V] = svd(B);
    Hm  = V * max(S, 1e-12) * V';
    A2  = (B + Hm) / 2;
    A_spd = (A2 + A2') / 2;
    % final check and perturb if necessary
    [~,p] = chol(A_spd);
    if p > 0
        A_spd = A_spd + eye(size(A_spd)) * 1e-8;
    end
end

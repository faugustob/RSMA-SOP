% correlated_nakagami_sim.m
% Simulate correlated Nakagami-m vectors using a Gaussian copula with
% exponential correlation, compute covariances and compare empirical
% correlation to target R.
%
% RNG seed fixed for reproducibility.

%% 1. Initialization and Parameters
rng(42);
N         = 6;                % number of IRS elements (vector dimension)
m1        = 1;                 % Nakagami-m for h_i
m2        = 2;                 % Nakagami-m for g_i
Omega1    = 5;                 % spread for h_i
Omega2    = 5;                 % spread for g_i
kappa     = 0.4;               % exponential correlation coefficient (0<kappa<1)
N_samples = 10e6;               % total Monte Carlo samples
batchSize = 1e5;               % process samples in batches to save memory

% Derived gamma scale parameters (theta)
theta1 = Omega1 / m1;          % scale for |h|^2
theta2 = Omega2 / m2;          % scale for |g|^2

%% 2. Build exponential correlation matrix R (NxN)
i = (1:N)'; j = 1:N;
R = kappa .^ abs(i - j);
R = (R + R')/2;
R(1:N+1:end) = 1;

% Ensure positive definite (tiny regularization if necessary)
minEig = min(eig(R));
if minEig <= 0
    warning('R not PD (min eig = %.3e). Regularizing...', minEig);
    R = nearestSPD(R);
    D = sqrt(diag(R));
    R = R ./ (D * D');   % re-normalize diag to 1
end

%% 3. Allocate accumulators (for means / covariances)
% We'll compute sample mean and covariance in online fashion (Welford-like).
sumH = zeros(1,N);
sumG = zeros(1,N);
sumHHt = zeros(N,N);   % accumulate H'*H across samples
sumGGt = zeros(N,N);

nDone = 0;
while nDone < N_samples
    nBatch = min(batchSize, N_samples - nDone);
    
    % 3a. Generate correlated uniforms using Gaussian copula
    % copularnd returns (nBatch x N) matrix of uniforms
    U_h = copularnd('Gaussian', R, nBatch);    % for H2 = |h|^2
    U_g = copularnd('Gaussian', R, nBatch);    % for G2 = |g|^2
    
    % 3b. Transform uniforms -> Gamma variates (|h|^2 and |g|^2)
    % gaminv(u, a, b) where a = shape (m), b = scale (theta)
    H2_batch = gaminv(U_h, m1, theta1);   % size nBatch x N
    G2_batch = gaminv(U_g, m2, theta2);
    
    % 3c. Form Nakagami magnitudes
    H_batch = sqrt(H2_batch);             % |h_i|
    G_batch = sqrt(G2_batch);             % |g_i|
    
    % 3d. Update accumulators
    sumH = sumH + sum(H_batch, 1);
    sumG = sumG + sum(G_batch, 1);
    sumHHt = sumHHt + (H_batch' * H_batch);   % N x N
    sumGGt = sumGGt + (G_batch' * G_batch);
    
    nDone = nDone + nBatch;
    fprintf('Processed %d / %d samples (%.1f%%)\n', nDone, N_samples, 100*nDone/N_samples);
end

%% 4. Finalize sample means and covariances
muH = sumH / N_samples;                 % 1 x N
muG = sumG / N_samples;

% Sample covariance: cov(H) = E[H H^T] - muH^T * muH
EHHt = sumHHt / N_samples;             % empirical E[H H^T]
EGGt = sumGGt / N_samples;

covH_emp = EHHt - (muH' * muH);        % N x N
covG_emp = EGGt - (muG' * muG);

% Empirical correlation matrices
stdH = sqrt(diag(covH_emp));
stdG = sqrt(diag(covG_emp));
corrH_emp = covH_emp ./ (stdH * stdH');   % elementwise division by outer product
corrG_emp = covG_emp ./ (stdG * stdG');

%% 5. Diagnostics and comparisons
fprintf('\n--- Diagnostics ---\n');
fprintf('Target Gaussian copula correlation (top-left 5x5):\n');
disp(R(1:5,1:5));
fprintf('Empirical correlation of |H| (top-left 5x5):\n');
disp(corrH_emp(1:5,1:5));
fprintf('Empirical correlation of |G| (top-left 5x5):\n');
disp(corrG_emp(1:5,1:5));

% Error metrics between empirical corr and target R
errH_frob = norm(corrH_emp - R, 'fro') / norm(R, 'fro');
errG_frob = norm(corrG_emp - R, 'fro') / norm(R, 'fro');
fprintf('Relative Frobenius error (corr H vs R): %.4f\n', errH_frob);
fprintf('Relative Frobenius error (corr G vs R): %.4f\n', errG_frob);

% Display covariance (top-left 5x5) for H
fprintf('\nEstimated covariance of H (top-left 5x5):\n');
disp(covH_emp(1:5,1:5));

%% 6. (Optional) Composite channel X = sum(H .* G) / sqrt(Omega1*Omega2)
% NOTE: recomputing X in batches to avoid storing all samples
nDone = 0;
sumX = 0; sumX2 = 0;
while nDone < N_samples
    nBatch = min(batchSize, N_samples - nDone);
    U_h = copularnd('Gaussian', R, nBatch);
    U_g = copularnd('Gaussian', R, nBatch);
    H2_batch = gaminv(U_h, m1, theta1);
    G2_batch = gaminv(U_g, m2, theta2);
    H_batch = sqrt(H2_batch);
    G_batch = sqrt(G2_batch);
    X_batch = sum(H_batch .* G_batch, 2) / sqrt(Omega1 * Omega2);  % nBatch x 1
    sumX = sumX + sum(X_batch);
    sumX2 = sumX2 + sum(X_batch.^2);
    nDone = nDone + nBatch;
end
muX = sumX / N_samples;
varX = (sumX2 / N_samples) - muX^2;
fprintf('\nComposite X: mean = %.6f, var = %.6f\n', muX, varX);

%% 7. Plot small diagnostic: histogram of X (subsampled)
subPlotSamples = min(50000, N_samples);  % limit plotted points
U_h = copularnd('Gaussian', R, subPlotSamples);
U_g = copularnd('Gaussian', R, subPlotSamples);
H2_s = gaminv(U_h, m1, theta1); H_s = sqrt(H2_s);
G2_s = gaminv(U_g, m2, theta2); G_s = sqrt(G2_s);
X_s = sum(H_s .* G_s, 2) / sqrt(Omega1 * Omega2);
figure;
histogram(X_s, 'Normalization', 'pdf', 'EdgeColor', 'none');
xlabel('Composite Channel Coefficient X');
ylabel('PDF');
title('Simulated PDF of Composite Channel Coefficient (subsample)');
grid on;

%% Subfunction: nearestSPD (same as your provided)
function A_spd = nearestSPD(A)
    B   = (A + A') / 2;
    [U,S,V] = svd(B);
    Hm  = V * max(S, 1e-12) * V';
    A2  = (B + Hm) / 2;
    A_spd = (A2 + A2') / 2;
    [~,p] = chol(A_spd);
    if p > 0
        A_spd = A_spd + eye(size(A_spd)) * 1e-8;
    end
end

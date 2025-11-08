% correlated_nakagami_sim_full.m
% Simulate correlated Nakagami-m vectors using a Gaussian copula with
% exponential correlation, compute covariances, and compare empirical
% correlation to theoretical one.

%% 1. Initialization and Parameters
rng(42);
N         = 5;                % number of IRS elements (vector dimension)
m1        = 1;                % Nakagami-m for h_i
m2        = 2;                % Nakagami-m for g_i
Omega1    = 5;                % spread for h_i
Omega2    = 5;                % spread for g_i
kappa     = 0.4;               % exponential correlation coefficient (0<kappa<1)
N_samples = 1e5;               % total Monte Carlo samples
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

%% 3. Allocate accumulators
sumH = zeros(1,N);
sumG = zeros(1,N);
sumZ = zeros(1,N);
sumHHt = zeros(N,N);
sumGGt = zeros(N,N);
sumZZt = zeros(N,N);

a1 = randn(1,N);
a2 = randn(1,N);

A1_batch = ones(batchSize,N).*a1;
A2_batch = ones(batchSize,N).*a2;


nDone = 0;
while nDone < N_samples
    nBatch = min(batchSize, N_samples - nDone);
    
    % 3a. Generate correlated uniforms using Gaussian copula
    U_h = copularnd('Gaussian', R, nBatch);    % for H^2
    U_g = copularnd('Gaussian', R, nBatch);    % for G^2
    
    % 3b. Transform uniforms -> Gamma variates (|h|^2 and |g|^2)
    H2_batch = gaminv(U_h, m1, theta1);
    G2_batch = gaminv(U_g, m2, theta2);
    
    % 3c. Form Nakagami magnitudes
    H_batch = sqrt(H2_batch);
    G_batch = sqrt(G2_batch);

    Z_batch = A1_batch.*H_batch .* G_batch.*A2_batch;

    % 3d. Update accumulators
    sumH = sumH + sum(H_batch, 1);
    sumG = sumG + sum(G_batch, 1);
    sumZ = sumZ + sum(Z_batch, 1);

    sumHHt = sumHHt + (H_batch' * H_batch);
    sumGGt = sumGGt + (G_batch' * G_batch);
    sumZZt = sumZZt + (Z_batch' * Z_batch);
    
    nDone = nDone + nBatch;
    fprintf('Processed %d / %d samples (%.1f%%)\n', nDone, N_samples, 100*nDone/N_samples);

end

%% 4. Finalize sample means and covariances
muH = sumH / N_samples;
muG = sumG / N_samples;
muZ = sumZ / N_samples;

EHHt = sumHHt / N_samples;
EGGt = sumGGt / N_samples;
EZZt = sumZZt / N_samples;

covH_emp = EHHt - (muH' * muH);
covG_emp = EGGt - (muG' * muG);
covZ_emp = EZZt - (muZ' * muZ);

stdH = sqrt(diag(covH_emp));
stdG = sqrt(diag(covG_emp));
stdZ = sqrt(diag(covZ_emp));

corrH_emp = covH_emp ./ (stdH * stdH');
corrG_emp = covG_emp ./ (stdG * stdG');
corrZ_emp = covZ_emp ./ (stdZ * stdZ');

%% 5. Diagnostics for H and G
fprintf('\n--- Diagnostics ---\n');
fprintf('Target Gaussian copula correlation (top-left 5x5):\n');
disp(R(1:5,1:5));
fprintf('Empirical correlation of |H| (top-left 5x5):\n');
disp(corrH_emp(1:5,1:5));
fprintf('Empirical correlation of |G| (top-left 5x5):\n');
disp(corrG_emp(1:5,1:5));

errH_frob = norm(corrH_emp - R, 'fro') / norm(R, 'fro');
errG_frob = norm(corrG_emp - R, 'fro') / norm(R, 'fro');
fprintf('Relative Frobenius error (corr H vs R): %.4f\n', errH_frob);
fprintf('Relative Frobenius error (corr G vs R): %.4f\n', errG_frob);

%% 6. Theoretical covariance of Z = H .* G
% Nakagami moments
muH_theo = gamma(m1 + 0.5)/gamma(m1) * sqrt(Omega1/m1);
muG_theo = gamma(m2 + 0.5)/gamma(m2) * sqrt(Omega2/m2);

varH_theo = Omega1 - (Omega1/m1) * (gamma(m1 + 0.5)/gamma(m1))^2;
varG_theo = Omega2 - (Omega2/m2) * (gamma(m2 + 0.5)/gamma(m2))^2;

% Ensure column vectors for means
muH_vec = muH_theo * ones(N,1);   % N x 1 column vector
muG_vec = muG_theo * ones(N,1);   % N x 1 column vector


covH_theo = varH_theo * R;        % N x N
covG_theo = varG_theo * R;        % N x N

covZ_theo = diag(a1)*diag(a2)*(covH_theo .* covG_theo ...
            + covH_theo .* (muG_vec * muG_vec') ...
            + covG_theo .* (muH_vec * muH_vec'))*diag(a1)*diag(a2);


% Theoretical correlation
stdZ_theo = sqrt(diag(covZ_theo));
corrZ_theo = covZ_theo ./ (stdZ_theo * stdZ_theo');

%% 7. Compare empirical vs theoretical Z
fprintf('\n--- Comparison of Z ---\n');
fprintf('Empirical covariance of Z (top-left 5x5):\n');
disp(covZ_emp(1:5,1:5));
fprintf('Theoretical covariance of Z (top-left 5x5):\n');
disp(covZ_theo(1:5,1:5));

fprintf('Empirical correlation of Z (top-left 5x5):\n');
disp(corrZ_emp(1:5,1:5));
fprintf('Theoretical correlation of Z (top-left 5x5):\n');
disp(corrZ_theo(1:5,1:5));

errZ_cov = norm(covZ_emp - covZ_theo,'fro') / norm(covZ_theo,'fro');
errZ_corr = norm(corrZ_emp - corrZ_theo,'fro') / norm(corrZ_theo,'fro');
fprintf('Relative Frobenius error (cov Z vs theoretical): %.4f\n', errZ_cov);
fprintf('Relative Frobenius error (corr Z vs theoretical): %.4f\n', errZ_corr);


%% Subfunction: nearestSPD
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

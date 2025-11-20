%% Covariance test for S = sum_{p,q} a_{p,q} .*( h_p .* h_q )
clear; clc; rng(2);

%% Parameters
N = 3;                 % vector dimension
Pp = 3;                % P_r (number of h_p)
Pq = 2;                % P_{r,j} (number of h_q)
m_val = 2.0;           % Nakagami m
Omega = 1.0;           % average power parameter for gamma
N_samples = 2e5;       % Monte Carlo samples
batchSize = 5e4;       % for copula generation

%% Correlation matrix R (exponential)
kappa = 0.6;
i = (1:N)'; j = 1:N;
R = kappa .^ abs(i - j);
R = (R + R')/2;
R(1:N+1:end) = 1;
% regularize if necessary
minEig = min(eig(R));
if minEig <= 0
    warning('R not PD (%.3e), regularizing...', minEig);
    R = nearestSPD(R);
    D = sqrt(diag(R));
    R = R ./ (D*D');
end

%% helper: correlated Nakagami sample generator (N_samples x N)
% Use Gaussian copula -> uniform -> gamma inverse -> r = sqrt(gamma)
nakagami_corr = @(m,Omega,R,nSamples) ...
    sqrt( gaminv( copularnd('Gaussian', R, nSamples), m, Omega/m ) );

%% 1) Construct deterministic a_{p,q} vectors
% Example: let g_{p,q} be scalar and hat_h_{P}^{(p,q)} be a deterministic vector:
% choose random deterministic vectors for illustration:
a = cell(Pp,Pq);
g_scalar = 0.8 + 0.4*rand(Pp,Pq);            % example scalars g_{p,q}
hatH = cell(Pp,Pq);
for p=1:Pp
    for q=1:Pq
        hatH{p,q} = 0.5 + 0.5*rand(1,N);     % deterministic vector (1xN) example
        a{p,q} = (g_scalar(p,q)) .* hatH{p,q};% a_{p,q} is 1xN
    end
end

%% 2) Generate families h_p and h_q (independent across p and q)
% For reproducibility we generate them separately; they share the same R
hP = cell(1,Pp);
hQ = cell(1,Pq);
for p = 1:Pp
    hP{p} = nakagami_corr(m_val,Omega,R,N_samples);   % N_samples x N
    % add a nonzero mean offset (elementwise) to each h_p if desired:
    hP{p} = hP{p} + 0.2*p;  % shifts all components by scalar 0.2*p
end
for q = 1:Pq
    hQ{q} = nakagami_corr(m_val,Omega,R,N_samples);
    hQ{q} = hQ{q} + 0.15*q;
end

%% 3) Empirical S and empirical covariance
S = zeros(N_samples, N);
for p = 1:Pp
    for q = 1:Pq
        % a{p,q} is 1xN, hP{p} and hQ{q} are N_samples x N
        S = S + (hP{p} .* hQ{q}) .* (ones(N_samples,1) * a{p,q});
    end
end
empCovS = cov(S);

%% 4) Empirical per-vector moments (means and covariances)
muP = cell(1,Pp); KP = cell(1,Pp);
muQ = cell(1,Pq); KQ = cell(1,Pq);
for p = 1:Pp
    muP{p} = mean(hP{p},1)';         % N x 1
    KP{p}  = cov(hP{p});             % N x N
end
for q = 1:Pq
    muQ{q} = mean(hQ{q},1)';         % N x 1
    KQ{q}  = cov(hQ{q});             % N x N
end

%% 5) Theoretical Cov(S) via the formula from the message
theoCovS = zeros(N,N);

% (A) diagonal terms p=p', q=q'
for p = 1:Pp
    for q = 1:Pq
        apq = a{p,q}(:); % make Nx1
        Aouter = apq * apq.'; % NxN
        T = KP{p} .* KQ{q} + KP{p} .* (muQ{q}*muQ{q}.' ) + (muP{p}*muP{p}.') .* KQ{q};
        theoCovS = theoCovS + (Aouter .* T);
    end
end

% (B) share p only: p=p', q != q'
for p = 1:Pp
    for q = 1:Pq
        for q2 = 1:Pq
            if q2==q, continue; end
            apq = a{p,q}(:); apq2 = a{p,q2}(:);
            A12 = apq * apq2.'; % NxN
            T = KP{p} .* ( muQ{q} * muQ{q2}.' );
            theoCovS = theoCovS + (A12 .* T);
        end
    end
end

% (C) share q only: q=q', p != p'
for q = 1:Pq
    for p = 1:Pp
        for p2 = 1:Pp
            if p2==p, continue; end
            apq = a{p,q}(:); ap2q = a{p2,q}(:);
            A12 = apq * ap2q.'; % NxN
            T = ( muP{p} * muP{p2}.' ) .* KQ{q};
            theoCovS = theoCovS + (A12 .* T);
        end
    end
end

%% 6) Display and compare
disp('Empirical Cov(S):'); disp(empCovS);
disp('Theoretical Cov(S):'); disp(theoCovS);
relErr = norm(empCovS - theoCovS,'fro') / norm(theoCovS,'fro');
fprintf('Relative Frobenius error = %.3e\n', relErr);

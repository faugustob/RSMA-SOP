%% CVX + MOSEK : 10 Canonical Convex Optimization Problems
% This script contains 10 different optimization problems.
% Each problem is independent and can be run by executing its section.
% Requirements:
%   - CVX installed
%   - MOSEK installed and licensed
%   - CVX configured to use MOSEK

clear all; % Clear workspace to start fresh
cvx_clear; % Reset CVX state to avoid any lingering model construction issues

cvx_solver mosek
cvx_precision high

% 1. Create a local 'vec' function on the fly


rng(1); % For reproducibility

%% 1. Mean-Variance Portfolio Optimization (Quadratic Program)
n = 10;
Sigma = randn(n); Sigma = Sigma'*Sigma;
mu = randn(n,1);
r_target = 0.1;

cvx_begin
    variable w(n)
    minimize( quad_form(w, Sigma) )
    subject to
        mu'*w >= r_target
        sum(w) == 1
        w >= 0
cvx_end

% Display results if successful
if strcmp(cvx_status, 'Solved')
    disp('Optimal weights:');
    disp(w);
    disp('Achieved return:');
    disp(mu'*w);
    disp('Portfolio variance:');
    disp(cvx_optval);
else
    disp('Optimization failed with status:');
    disp(cvx_status);
end

%% 2. LASSO Regression (SOCP)
m = 50; n = 20;
A = randn(m,n);
b = randn(m,1);
lambda = 0.5;

cvx_begin
    variable x(n)
    minimize( sum_square(A*x - b) + lambda*norm(x,1) )
cvx_end


%% 3. Matrix Completion via Nuclear Norm Minimization (SDP)
n = 8;
M = randn(n);
Omega = rand(n) > 0.5;

cvx_begin sdp
    variable X(n,n)
    minimize( norm_nuc(X) )
    subject to
        X(Omega) == M(Omega)
cvx_end


%% 4. Robust Least Squares (SOCP)
m = 40; n = 10;
A = randn(m,n);
b = randn(m,1);
eps = 0.1;

cvx_begin
    variables x(n) t
    minimize(t)
    subject to
        norm(A*x - b,2) + eps*norm(x,2) <= t
cvx_end


%% 5. Max-Cut Semidefinite Relaxation (SDP)
n = 6;
W = rand(n); W = (W + W')/2;

cvx_begin sdp
    variable X(n,n) symmetric
    maximize( 0.25*sum(sum(W.*(1 - X))) )
    subject to
        diag(X) == 1
        X == semidefinite(n)
cvx_end


%% 6. Chebyshev Center of a Polytope (SOCP)
m = 20; n = 5;
A = randn(m,n);
b = ones(m,1);

cvx_begin
    variables x(n) r
    maximize(r)
    subject to
        for i = 1:m
            A(i,:)*x + r*norm(A(i,:),2) <= b(i);
        end
        r >= 0
cvx_end


%% 7. Optimal Transport (Earth Mover's Distance) (LP)
n = 10;
p = rand(n,1); p = p/sum(p);
q = rand(n,1); q = q/sum(q);
C = rand(n,n);

cvx_begin
    variable T(n,n)
    minimize( sum(sum(C .* T)) )
    subject to
        T >= 0
        sum(T,2) == p
        sum(T,1)' == q
cvx_end


%% 8. Geometric Programming: Power Control (GP)
n = 4;
G = rand(n,1);
H = rand(n,n);
sigma = 0.1*ones(n,1);
gamma = ones(n,1);

cvx_begin gp
    variables pwr(n)
    minimize( sum(pwr) )
    subject to
        for i = 1:n
            G(i)*pwr(i) >= gamma(i)*( sum(H(i,[1:i-1,i+1:end]).*pwr([1:i-1,i+1:end])) + sigma(i) );
        end
cvx_end


%% 9. H-infinity Control via LMIs (SDP)
n = 3; m = 2;
A = randn(n);
B = randn(n,m);
C = randn(m,n);
D = zeros(m,m);

cvx_begin sdp
    variables P(n,n) gamma
    minimize(gamma)
    subject to
        P == semidefinite(n)
        [A'*P+P*A, P*B, C';
         B'*P, -gamma*eye(m), D';
         C, D, -gamma*eye(m)] <= 0
        gamma >= 0
cvx_end


%% 10. Log-Determinant Maximization (Maximum Entropy Covariance) (SDP)
n = 5;
S = randn(n); S = S'*S + eye(n);

cvx_begin sdp
    variable X(n,n) symmetric
    maximize( log_det(X) )
    subject to
        X <= S
        X == semidefinite(n)
cvx_end

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
%% 2D SCA Animation with CVX
% Problem: Minimize 0.5*||w||^2 - 15*log(1 + ||w||^2)
% This creates a 2D landscape with a "volcano" peak at the origin.

clear; clc; close all;

% --- 1. Parameters ---
wk = [8; 8];        % Starting point in 2D space
rho = 10;           % High proximal weight to force 100+ iterations
tol = 1e-4;
max_iters = 300;
history = wk;       % Store coordinates [x; y]

% Function definitions
f = @(x, y) 0.5*(x.^2 + y.^2) - 15*log(1 + x.^2 + y.^2);
% Gradient of the non-convex part v(w) = 15*log(1 + x^2 + y^2)
grad_v = @(w) (15 * 2 * w) / (1 + norm(w)^2);

fprintf('Iter |    x    |    y    |  Step Size\n');
fprintf('------------------------------------\n');

% --- 2. The SCA Loop ---
for k = 1:max_iters
    % Current constants for this iteration
    vk_val = 15 * log(1 + norm(wk)^2);
    gk_vec = grad_v(wk);
    
    cvx_begin quiet
        variable w_next(2)
        
        % The 2D Surrogate:
        % Convex Part: 0.5 * quad_over_lin(w_next, 1)
        % Linearized Part: (vk_val + gk_vec' * (w_next - wk))
        % Proximal Term: (rho/2) * quad_over_lin(w_next - wk, 1)
        
        minimize( 0.5*sum_square(w_next) - (vk_val + gk_vec'*(w_next - wk)) + (rho/2)*sum_square(w_next - wk) )
        
        subject to
            w_next >= -10;
            w_next <= 10;
    cvx_end
    
    step_dist = norm(w_next - wk);
    history = [history, w_next];
    
    if mod(k, 10) == 0
        fprintf('%4d | %7.3f | %7.3f | %9.6f\n', k, w_next(1), w_next(2), step_dist);
    end
    
    if step_dist < tol
        fprintf('Converged at iteration %d!\n', k);
        break;
    end
    
    wk = w_next;
end
% =====================================================
% 2D Contour + 3D Surface Animation (Two Figures)
% =====================================================

% --- Figures ---
figContour = figure('Color','w','Position',[100 100 520 520]);
figSurface = figure('Color','w','Position',[650 100 650 520]);

% --- Grid and surface ---
[X, Y] = meshgrid(linspace(-10,10,80));
Z = f(X,Y);

for i = 1:2:size(history,2)

    % =================================================
    % FIGURE 1: 2D CONTOUR (Top View)
    % =================================================
    figure(figContour)
    clf

    contourf(X, Y, Z, 30, 'LineColor','none')
    colormap(parula)
    colorbar
    hold on

    % Optimization path
    plot(history(1,1:i), history(2,1:i), ...
        'r-', 'LineWidth',2)

    % Current point
    plot(history(1,i), history(2,i), ...
        'ko', 'MarkerFaceColor','y', 'MarkerSize',8)

    title(['SCA Optimization Path (Iteration ', num2str(i), ')'], ...
        'FontWeight','bold')
    xlabel('x'); ylabel('y')
    axis square tight
    set(gca,'FontSize',11)
    hold off


    % =================================================
    % FIGURE 2: 3D SURFACE VIEW
    % =================================================
    figure(figSurface)
    clf

    surf(X, Y, Z, ...
        'EdgeColor','none', ...
        'FaceAlpha',0.85)
    hold on

    shading interp
    colormap(turbo)
    camlight headlight
    lighting gouraud
    set(gca,'SortMethod','childorder')

    % Path floating above the surface
    plot3(history(1,1:i), history(2,1:i), ...
          f(history(1,1:i), history(2,1:i)) + 0.6, ...
          'r-', 'LineWidth',2)

    % Current / optimum point (always visible)
    z_opt = f(history(1,i), history(2,i));
    plot3(history(1,i), history(2,i), z_opt + 2, ...
          'kp', 'MarkerFaceColor','y', ...
          'MarkerSize',14, 'LineWidth',1.5)

    % Vertical beacon line
    plot3([history(1,i) history(1,i)], ...
          [history(2,i) history(2,i)], ...
          [min(Z(:)) z_opt + 2], ...
          'k--', 'LineWidth',1.2)

    view(-35,40)
    xlabel('x'); ylabel('y'); zlabel('f(x,y)')
    title('Objective Surface','FontWeight','bold')
    set(gca,'FontSize',11)
    hold off

    drawnow
    pause(0.2)   % Animation speed control
end

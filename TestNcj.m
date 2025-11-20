%% Monte Carlo Simulation for N_{c,j} with correlated Nakagami
clc; clear; close all;

%% Simulation parameters
numIter = 1e6; % Monte Carlo realizations

% Nakagami parameters
m_r = 2; Omega_r = 1; % RIS contamination
m_j = 2; Omega_j = 1; % RIS reflected
m_g = 2; Omega_g = 1; % LOS

% Number of elements
N_R = 30;  % RIS elements
P_e = 2;  % LOS interference

P_r = 2; % Number of paths from LEO to RIS
P_rj = 4; %Number of paths from RIS to user j.

M=1;
N=1;

% Deterministic coefficients (example)
alpha_c = 1; P_LJ = 1; P_LOS = 1; P_t = 1;
beta_r = rand([1 N_R]);

% Deterministic matrices (identity for illustration)
H_RIS = eye(N_R);
H_LOS = eye(N_R);
g_pq = randn([P_r,P_rj]);

% Initialize the cell array
A = cell(P_r, P_rj);
B = cell(P_e);

% Fill each cell with a MN x MN random matrix
for i = 1:P_r
    for j = 1:P_rj
        A{i,j} = randn(M*N, M*N);
    end
end

for j = 1:P_e
    B{j} = randn(M*N, M*N);
end

%% Correlation matrix R (exponential)
kappa = 0.6;
i = (1:N_R)'; j = 1:N_R;
R = kappa .^ abs(i - j);
R = (R + R')/2;
R(1:N_R+1:end) = 1;

% Regularize if not PD
minEig = min(eig(R));
if minEig <= 0
    warning('R not PD (%.3e), regularizing...', minEig);
    R = nearestSPD(R);
    D = sqrt(diag(R));
    R = R ./ (D*D');
end

%% Correlated Nakagami sample generator (Gaussian copula)
nakagami_corr = @(m,Omega,R,nSamples) ...
    sqrt( gaminv( copularnd('Gaussian', R, nSamples), m, Omega/m ) );

%% Monte Carlo
Ncj_MC = zeros(numIter,1);

for n = 1:numIter
    % Generate correlated Nakagami vectors
    h_r_M = []; %N_R x P_r
    h_j_M = []; %N_R x P_rj
    g_eu = sqrt(gamrnd(m_g, Omega_g/m_g, P_e, 1)); % LOS independent

    for path_p=1:P_r
        h_r_v = nakagami_corr(m_r, Omega_r, R, 1)'; % 1xN_R
        h_r_M = [h_r_M,h_r_v];
    end

    for path_q=1:P_rj
        h_j_v = nakagami_corr(m_j, Omega_j, R, 1)'; % 1xN_R
        h_j_M = [h_j_M,h_j_v];
    end

    %% ----------------Numerical simulation For-loop method ----------------

     result_NR = 0;
     
     for r = 1:N_R
        result_loop = 0;
        for p = 1:P_r
            for q = 1:P_rj
                result_loop = result_loop + h_r_M(r,p) .* h_j_M(r,q)*g_pq(p,q)* A{p,q};
            end
        end
               result_NR = result_NR + beta_r(r)*result_loop;
     end

     result_NR = sqrt(P_LJ)*result_NR;

     resut_E=0;

     for e=1:P_e
         resut_E = resut_E+g_eu(e)*B{e};
     end
     resut_E = sqrt(P_LOS)*resut_E;

     Ncj_MC(n)=norm(result_NR+resut_E,"fro")^2;

    
end

var_NCj = var(Ncj_MC);
mean_NCj = mean(Ncj_MC);

% Gamma parameters using lambda (rate)
k = mean_NCj^2 / var_NCj;          % shape
lambda = mean_NCj / var_NCj;       % rate

% Define gamma PDF using lambda
f_N = @(x) (lambda^k / gamma(k)) * x.^(k-1) .* exp(-lambda * x);

% Example: evaluate PDF at some points
x = linspace(0, 2*mean_NCj, 100);
y = f_N(x);

%% Plot results
figure;
histogram(Ncj_MC,'Normalization','pdf','DisplayName','Monte Carlo');
hold on;
plot(x, y, 'r','LineWidth',2,'DisplayName','Approx. PDF');
xlabel('N_{c,j}');
ylabel('PDF');
title('Empirical vs Theoretical PDF of N_{c,j} (Correlated Nakagami)');
legend;
grid on;



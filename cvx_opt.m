%% SCA for RSMA Power Allocation (Subproblem 1) - Corrected Version
clear; clc;
cvx_clear;

% --- Simulation Parameters ---
K = 3;              
ln2 = log(2);
sigma2 = 1e-2;      
P_total = 1;        
R_min = 0.5;        
R_c0 = 0.2;         

% --- Dummy Channel Gains (Constants for Subproblem 1) ---
P_legit = abs(randn(K,1)).^2 + 5; 
P_eve = abs(randn(K,1)).^2 + 1;   
A_legit = 0.5;                     
A_eve = 0.5;                       

% --- Initialization ---
alpha_k = ones(K+1, 1) / (K+1); 
prev_val = 0; % FIXED: Initialized prev_val
max_iters = 20;
tolerance = 1e-5;

fprintf('Iter | Avg Secrecy Rate | alpha_c \n');
fprintf('------------------------------------\n');

for k = 1:max_iters
    % Pre-calculate interference and total power at anchor alpha_k
    % These constants are used for the Taylor expansions
    I_legit_k = zeros(K,1);
    I_eve_k = zeros(K,1);
    total_eve_k = zeros(K,1);
    
    for j = 1:K
        I_legit_k(j) = (sum(alpha_k(2:end)) - alpha_k(j+1))*P_legit(j) + A_legit + sigma2;
        I_eve_k(j)   = (sum(alpha_k(2:end)) - alpha_k(j+1))*P_eve(j) + A_eve + sigma2;
        total_eve_k(j) = (sum(alpha_k(2:end)))*P_eve(j) + A_eve + sigma2;
    end
    
    cvx_begin quiet
        cvx_solver mosek 
        variable vecAlpha(K+1) nonnegative 
        variable s_j(K)     
        
        % Auxiliary variables for denominators to keep code clean
        variable I_legit(K)
        variable I_eve(K)

        maximize( sum(s_j) / K )

        subject to
            sum(vecAlpha) <= P_total;

            for j = 1:K
                % Define denominators as linear constraints
                I_legit(j) == (sum(vecAlpha(2:end)) - vecAlpha(j+1))*P_legit(j) + A_legit + sigma2;
                I_eve(j)   == (sum(vecAlpha(2:end)) - vecAlpha(j+1))*P_eve(j) + A_eve + sigma2;
                
                % --- Legitimate Rate (Lower Bound) ---
                % Concave part: log(Total Power)
                log_total_legit = log( (sum(vecAlpha(2:end)))*P_legit(j) + A_legit + sigma2 )/ln2;
                % Linear part: Taylor expansion of log(Interference)
                log_I_legit_linear = log(I_legit_k(j))/ln2 + (1/(ln2*I_legit_k(j)))*(I_legit(j) - I_legit_k(j));
                
                R_legit = log_total_legit - log_I_legit_linear;

                % --- Eavesdropper Rate (Upper Bound) ---
                % Linear part: Taylor expansion of log(Total Power)
                log_total_eve_linear = log(total_eve_k(j))/ln2 + (1/(ln2*total_eve_k(j)))*((sum(vecAlpha(2:end)))*P_eve(j) - (sum(alpha_k(2:end)))*P_eve(j));
                % Concave part: log(Interference)
                log_I_eve = log(I_eve(j))/ln2;
                
                R_eve_approx = log_total_eve_linear - log_I_eve;

                % Secrecy Slack
                s_j(j) <= R_legit - R_eve_approx;
                s_j(j) >= 0;
            end
            
            % --- QoS Constraints ---
            for j = 1:K
                vecAlpha(1)*P_legit(j) >= (2^R_c0 - 1) * I_legit(j);
                (vecAlpha(1) + vecAlpha(j+1))*P_legit(j) >= (2^R_min - 1) * I_legit(j);
            end
    cvx_end

    % 3. Check for convergence
    if k > 1 && abs(cvx_optval - prev_val) < tolerance
        fprintf('Converged at iteration %d\n', k);
        break;
    end
    
    prev_val = cvx_optval;
    alpha_k = vecAlpha;
    %fprintf('%4d | %16.4f | %7.4f\n', k, cvx_optval, vecAlpha(1));
end
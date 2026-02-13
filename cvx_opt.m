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

        maximize( sum(s_j) / (K) )

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

%%
display('Convex approximatin with AO');

Num_agents  = Num_agents;
Max_iteration = Max_iteration;
Rmin=0.1;

% Check if more than one STAR-RIS side is being used.
any_reflect = any(reflect > 0) && any(reflect < 0);

% Problem bounds and dimensionality
dim = K+1+Nr;
ub=[ones(1,K+1),2*pi*ones(1,Nr)];
alpha_min = 1e-4;
lb = [alpha_min * ones(1,K+1),zeros(1,Nr)];
zeta_k_St = ones(1,Nr); % RIS amplitude coefficients, we may use it to boost for active RIS

Active_Gain_dB = 0; 
zeta_k_St = (10^(Active_Gain_dB/10)) * ones(1, Nr);


% zeta_k_Sr = rand(Num_agents,Nr); % reflection coefficients
phi_Sr = 2*pi*rand(Num_agents,Nr);
phi_St = 2*pi*rand(Num_agents,Nr);% transmission phases


alpha = rand(Num_agents, K+1); 
% random values
alpha = alpha ./ sum(alpha, 2);      % divide each row by its row sum

alpha = alpha - (sum(alpha,2)-1)/(K+1);
alpha = alpha - (sum(alpha,2)-1)/(K+1);
alpha = alpha - (sum(alpha,2)-1)/(K+1);


X = [alpha,phi_St];
% phi_init = -angle(h_jq(:, 1, 1) .* h_rp(:, 1, 1, 1)); 
% 
% % phi_init is now (Nr x 1), we transpose it to fit the agent row (1 x Nr)
% X(1, K+2:K+1+Nr) = phi_init';

if any_reflect
    dim = K+1+3*Nr;
    ub=[ones(1,K+1),2*pi*ones(1,2*Nr),ones(1,Nr)];
    alpha_min = 1e-4;
    lb = [alpha_min * ones(1,K+1),zeros(1,3*Nr)];
    zeta_k_St = (10^(Active_Gain_dB/10)) *rand(Num_agents,Nr);
    X = [alpha,phi_Sr,phi_St,zeta_k_St];
end


% --- Problem Dimensions and Bounds ---
dim_pso = dim;
alpha_min_pso = alpha_min;
lb_pso =lb;
ub_pso = ub;



Destination_position=zeros(1,dim);
Destination_fitness=inf;
best_fake_secrecy_rate=0;
best_real_secrecy_rate = 0;

Convergence_curve_AO=zeros(1,Max_iteration);
sum_rate_curve=zeros(1,Max_iteration);

min_sum_secrecy = zeros(1,Max_iteration);

Objective_values = zeros(1,size(X,1));
%All_objective_values=zeros(Max_iteration,size(X,1));


% Calculate the fitness of the first set and find the best one
for i=1:size(X,1)
    C_k = zeros(K,1);

    [sc_c_lk,sc_p_lk,sc_p_kk,rate_c,rate_k,R_k,~] = compute_sinr_sc_an(Pe,P,Q_j,nF+L,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,Rmin,h_rp,h_jq,h_e,zeta_k_St,Active_Gain_dB,X(i,:));

    %sum_secrecy = sc_c_lk+sc_p_lk; %Private + Common secrecy capacities.

    mean_fake_p_secrecy = mean(mean(sc_p_lk(1:nF,:)));
    mean_p_secrecy = mean(mean(sc_p_lk(nF+1:end,:)));

    
    penalty = 0;
    violation = max(Rmin - R_k, 0);
    penalty = penalty + sum(violation.^2);

    
    Objective_values(1,i) = -mean_fake_p_secrecy + 1e3 * penalty;

    if i==1
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
        best_fake_secrecy_rate = mean_fake_p_secrecy;
        best_real_secrecy_rate = mean_p_secrecy;
    elseif Objective_values(1,i)<Destination_fitness
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
        best_fake_secrecy_rate = mean_fake_p_secrecy;
        best_real_secrecy_rate = mean_p_secrecy;
    end

    All_objective_values(1,i)=Objective_values(1,i);
end

alpha_idx = 1:(K+1);
ris_idx = (K+2):size(X,2);

% Main loop - Alternating Optimization (alpha first)
t = 2;
while t <= Max_iteration
    

    % === Substep 1: Optimize alpha first (phases/amplitudes fixed) ===
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
     X(:,alpha_idx) = alpha_k;
    alpha_fixed = X(:,alpha_idx);

    % === Substep 2: Optimize non-alpha variables (phases/amplitudes/zeta) with new alpha fixed ===
    for i = 1:size(X, 1)
        X(i,alpha_idx) = alpha_fixed(i,:);
        for j = ris_idx  % phase and zeta dimensions
            r2 = (2*pi) * rand();
            r3 = 2 * rand();
            r4 = rand();
            if r4 < 0.5
                X(i,j) = X(i,j) + (r1 * sin(r2) * abs(r3 * Destination_position(j) - X(i,j)));
            else
                X(i,j) = X(i,j) + (r1 * cos(r2) * abs(r3 * Destination_position(j) - X(i,j)));
            end
        end

        % Bound check
        Flag4ub = X(i,ris_idx) > ub(ris_idx);
        Flag4lb = X(i,ris_idx) < lb(ris_idx);
        
        X(i,ris_idx) = ...
            X(i,ris_idx).*(~(Flag4ub+Flag4lb)) + ...
            ub(ris_idx).*Flag4ub + ...
            lb(ris_idx).*Flag4lb;

       
        % Evaluate objective again
        [sc_c_lk, sc_p_lk, sc_p_kk, rate_c, rate_k, R_k, ~] = compute_sinr_sc_an(Pe, P, Q_j, nF+L, K, delta_f, Plos, PLj, Nr, HB, HA, g_pq, Nsymb, reflect, Rmin, h_rp, h_jq, h_e, zeta_k_St, Active_Gain_dB, X(i,:));
        mean_fake_p_secrecy = mean(mean(sc_p_lk(1:nF,:)));
        mean_p_secrecy = mean(mean(sc_p_lk(nF+1:end,:)));

        penalty = 0;
        violation = max(Rmin - R_k, 0);
        penalty = penalty + sum(violation.^2);

        Objective_values(1,i) = -mean_fake_p_secrecy + 1e3 * penalty;

        % Update global best
        if Objective_values(1,i) < Destination_fitness
            Destination_position = X(i,:);
            Destination_fitness = Objective_values(1,i);
            best_fake_secrecy_rate = mean_fake_p_secrecy;
            best_real_secrecy_rate = mean_p_secrecy;
        end
    end

    % Record curves
    Convergence_curve_AO(t) = -Destination_fitness;
    Fake_secrecy_rate_curve_AO(t) = best_fake_secrecy_rate;
    Real_secrecy_rate_curve_AO(t) = best_real_secrecy_rate;

    if mod(t,1) == 0
        display(['AO At iteration ', num2str(t), ' the optimum fake sc is ', num2str(best_fake_secrecy_rate), ' the optimum real sc is ', num2str(best_real_secrecy_rate)]);
    end

    t = t + 1;
end

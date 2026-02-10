%% Successive Convex Approximation (SCA) with CVX
% Problem: Minimize f(x) = -log2(1 + x) + log2(1 + x/2)
% Constraints: x >= 1;
clear all; % Clear workspace to start fresh
cvx_clear; % Reset CVX state to avoid any lingering model construction issues
cvx_solver mosek
cvx_precision high

% --- Initial Setup ---
ln2 = log(2); % Natural log of 2 for log2 emulation
xk = 0.5; % Our starting point (The "Anchor")
tolerance = 1e-4; % Stopping criteria
max_iters = 50; % Safety exit
history = xk; % Store for visualization

% Define the actual objective for printing (using log/ln2)
f = @(x) -log(1 + x)/ln2 + log(1 + x/2)/ln2;

fprintf('Iteration | x_k | f(x_k) \n');
fprintf('-----------------------------\n');

for k = 1:max_iters
    % 1. Calculate values at the current anchor point xk
    % These are CONSTANTS for the CVX block
    vk = @(x) log(1 + x/2)/ln2; % Value of the non-convex part at xk
    grad_vk = @(x) 1/(2*ln2*(1 + x/2)); % Gradient of the non-convex part at x

    % 2. Solve the CONVEX sub-problem using CVX
    cvx_begin quiet % 'quiet' hides the solver output each loop
        variable x_next
        % THE SURROGATE OBJECTIVE:
        % We keep the "nice" convex part -log2(1 + x) as is.
        % We replace the "nasty" part +log2(1 + x/2) with its linear Taylor expansion.
        % Linearization: f(x) \approx log2(1 + xk/2) + f'(xk)*(x - xk)
        minimize( -log(1 + x_next)/ln2 + (vk(xk) + grad_vk(xk) * (x_next - xk)) )
        subject to
            x_next >= 0;
    cvx_end

    % 3. Check for convergence
    % If the solver found a point very close to where we already are, stop.
    if abs(x_next - xk) < tolerance
        break;
    end

    % 4. Update the anchor for the next iteration
    xk = x_next;
    history = [history, xk];
    fprintf('%8d | %7.4f | %8.4f\n', k, xk, f(xk));
end

%--- Animation Logic ---
% This visualizes the "Bowl" (Surrogate) moving as the algorithm iterates.

f = @(x) -log2( 1 + x) + log(1+x/2); % Original function
x_range = linspace(-3, 3, 200);

figure('Color', 'w');
for i = 1:length(history)
    curr_xk = history(i);
    
    % Build the surrogate for plotting
    % g(x) = x^4 - (v(xk) + grad_v(xk)*(x-xk))
    v_val = log2(1 + curr_xk/2);
    g_val = 1/(2*log(2)*(1 + curr_xk/2));
    surrogate_plot = log2(1 + x_range) - (v_val + g_val * (x_range - curr_xk));
    
    clf; hold on; grid on;
    plot(x_range, f(x_range), 'k', 'LineWidth', 2); % The real landscape
    plot(x_range, surrogate_plot, 'r--', 'LineWidth', 1.5); % The "Fake" landscape
    plot(curr_xk, f(curr_xk), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    
    title(['SCA Iteration: ', num2str(i)]);
    xlabel('x'); ylabel('f(x)');
    legend('True Objective', 'CVX Surrogate', 'Current x_k');
    ylim([-10, 30]);
    pause(0.5);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
display('SCA is optimizing your problem AO');

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
    % Eq. (3.4)
    a = 3;
    r1 = a - t * (a / Max_iteration);  % r1 decreases linearly

    % === Substep 1: Optimize alpha first (phases/amplitudes fixed) ===
    for i = 1:size(X, 1)
        for j = alpha_idx  % only alpha dimensions
            r2 = (2*pi) * rand();
            r3 = 2 * rand();
            r4 = rand();
            if r4 < 0.5
                X(i,j) = X(i,j) + (r1 * sin(r2) * abs(r3 * Destination_position(j) - X(i,j)));
            else
                X(i,j) = X(i,j) + (r1 * cos(r2) * abs(r3 * Destination_position(j) - X(i,j)));
            end
        end

        % Project alpha onto simplex (sum = 1, >= alpha_min)
             % Bound check
        Flag4ub = X(i,alpha_idx) > ub(alpha_idx);
        Flag4lb = X(i,alpha_idx) < lb(alpha_idx);
        
        X(i,alpha_idx) = ...
            X(i,alpha_idx).*(~(Flag4ub+Flag4lb)) + ...
            ub(alpha_idx).*Flag4ub + ...
            lb(alpha_idx).*Flag4lb;

        X(i,alpha_idx) = X(i,alpha_idx) ./ sum(X(i,alpha_idx));


        % Optional repeated correction for floating-point precision (keep your original style)
        X(i,alpha_idx) = X(i,alpha_idx) - (sum(X(i,alpha_idx))-1)/(K+1);
        X(i,alpha_idx) = X(i,alpha_idx) - (sum(X(i,alpha_idx))-1)/(K+1);
        X(i,alpha_idx) = X(i,alpha_idx) - (sum(X(i,alpha_idx))-1)/(K+1);

        % Evaluate objective with updated alpha
        [sc_c_lk, sc_p_lk, sc_p_kk, rate_c, rate_k, R_k, ~] = compute_sinr_sc_an(Pe, P, Q_j, nF+L, K, delta_f, Plos, PLj, Nr, HB, HA, g_pq, Nsymb, reflect, Rmin, h_rp, h_jq, h_e, zeta_k_St, Active_Gain_dB, X(i,:));
        mean_fake_p_secrecy = mean(mean(sc_p_lk(1:nF,:)));
        mean_p_secrecy = mean(mean(sc_p_lk(nF+1:end,:)));

        penalty = 0;
        violation = max(Rmin - R_k, 0);
        penalty = penalty + sum(violation.^2);

        Objective_values(1,i) = -mean_fake_p_secrecy + 1e3 * penalty;

        % Update global best (alpha-optimized solution)
        if Objective_values(1,i) < Destination_fitness
            Destination_position = X(i,:);
            Destination_fitness = Objective_values(1,i);
            best_fake_secrecy_rate = mean_fake_p_secrecy;
            best_real_secrecy_rate = mean_p_secrecy;
        end
    end
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
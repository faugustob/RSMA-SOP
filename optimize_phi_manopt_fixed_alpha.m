function [phi_St,cost_opt] = optimize_phi_manopt_fixed_alpha(Rmin,L_node,E_node,problem,b0,alpha,K, nF, sigma2, Pw, AN_P_ratio, Ck)

     
        %% ========================= MANOPT =========================
        




        s_param = 50; % Smoothing parameter for max-min
        lambda_penalty = 1e3;

        % 
        %   X = [alpha,beta.'];
        % 
        % [~,sc_p_lk,~] = compute_sinr_sc_an_manopt(Pe,P,Q_j,nF,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,Rmin,h_rp,h_jq,h_e,zeta_k_St,Active_Gain_dB,X);
        % 
        % f_min =  -(1/s_param) * log(sum(exp(-s_param * (R_sec)),'all'));
        %min_Rsec = min(min(R_sec));

        % grad_val = grad_func(beta, L_node, E_node, alpha, s_param, K, nF, sigma2, Pw, AN_P_ratio);
         
        Nr = length(b0);

    
        problem.cost = @(b) cost_func(b, L_node, E_node, alpha, s_param, K, nF, sigma2, Pw, AN_P_ratio, Rmin, lambda_penalty,Ck);
        %problem = manoptAD(problem);

        problem.egrad = @(b) grad_func(b, L_node, E_node, alpha, s_param, K, nF, sigma2, Pw, AN_P_ratio, Rmin, lambda_penalty,Ck);

        problem.ehess = @(b,v) hess_func(b,v, L_node, E_node, alpha, s_param, K, nF, sigma2, Pw, AN_P_ratio, Rmin, lambda_penalty,Ck);
        % 
        % % Initial guess
        % b0 = manifold.rand();
        % checkgradient(problem);
       % checkhessian(problem);

        % Solve
       % fprintf('Starting Manifold Optimization...\n');
        % Define options structure
        options.tolgradnorm = 1e-4;    % Stop when the gradient norm is very small
        options.maxiter = 30;         % Maximum outer iterations
        options.verbosity = 0;          % 2 shows summary, 3 shows detailed inner steps
        options.linesearch = @linesearch; % Trust-regions usually manages step size via the radius
        
        % Inner iteration control (Krylov steps)
        % Since your Hessian is perfect, we can allow more inner iterations 
        % to solve the sub-problem more accurately.
        options.maxinner = Nr/2; 
        
        % Execute with options
        [beta_opt, cost_opt, info] = trustregions(problem, b0, options);
        
        phi_St = angle(beta_opt).';
        
end



function f_val = cost_func(beta, L_node, E_node, alpha, s_param, K, nF, sigma2, Pw, AN_P_ratio, Rmin, lambda_penalty,Ck)

        %  X = [alpha,beta.'];
        % 
        % [~,sc_p_lk,~] = compute_sinr_sc_an_manopt(Pe,P,Q_j,nF,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,Rmin,h_rp,h_jq,h_e,zeta_k_St,Active_Gain_dB,X);

       
    [R_sec,rate_p,~] = get_Secrecy_matrix(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);

    % Log-Sum-Exp smoothing for max(-R)
    %f_val =  (1/s_param) * log(sum(exp(-s_param * (R_sec)),'all'));
    Z = -s_param * R_sec;
    m = max(Z, [], 'all');                 % stabilization shift
    f_lse = (1/s_param) * ( m + log(sum(exp(Z - m), 'all')) );
    
    % --- NEW: Penalty Term ---
    % Quadratic penalty for violating Rmin
    penalty = lambda_penalty * sum(max(0, Rmin - (rate_p+Ck)).^2, 'all');  
    f_val = f_lse + penalty;
    
end

function [grad_Rsec_lk,g_user] = grad_R_sec(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio)
    Nr = length(beta);
    ln2 = log(2);
    noise_total = (sigma2/Pw);
   
    alpha_pi = alpha(2:end);    
   
    % 1. Compute all current rates and power terms
  
   
   [~,~,~,T_k,T_l,I_k,I_l,Ak,Pk,Al,Pl] = get_Secrecy_matrix(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);

    % 2. Log-Sum-Exp Weights    
   
    grad_Rsec_lk = zeros(Nr,nF,K);
    g_user        = zeros(Nr, K);          % <-- NEW
   
    % 3. Gradient Accumulation
    g = zeros(Nr, 1);   % already correct
    
    for k = 1:K
        % ==================== FIXED GRADIENT LINES ====================
        grad_Pk = 2 * (L_node(k).V1 * beta) + 2 * L_node(k).V2;           
        grad_Ak = 2 * (L_node(k).V1_AN * beta) + 2 * L_node(k).V2_AN;    
        
                    
        g_user_k = (alpha_pi(k)/(ln2*T_k(k)*I_k(k))) * ((Ak(k)+noise_total)*grad_Pk-Pk(k)*grad_Ak);
        g_user(:,k) = g_user_k;               % <-- NEW
        
        for l = 1:nF
            grad_Pl = 2 * (E_node(l).V1_l * beta) + 2 * E_node(l).V2_l;          
            grad_Al = 2 * (E_node(l).V1_AN_l * beta) + 2 * E_node(l).V2_AN_l;                
                        
            g_eav_l = (alpha_pi(k)/(ln2*T_l(l,k)*I_l(l,k))) * ((Al(l)+noise_total)*grad_Pl-Pl(l)*grad_Al);

            grad_Rsec_lk(:,l,k) = g_user_k - g_eav_l;
            
        end
    end
end

function [g, grad_sec_lk] = grad_func(beta, L_node, E_node, alpha, s_param, K, nF, sigma2, Pw, AN_P_ratio, Rmin, lambda_penalty,Ck)
    Nr = length(beta);   
    [R_sec,rate_p,~] = get_Secrecy_matrix(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);

    % 2. Log-Sum-Exp Weights      
    Z = -s_param * R_sec;

    m = max(Z, [], "all");        % stabilization shift
    weights = exp(Z - m);
    
    weights = weights / sum(weights, "all");
    [grad_sec_lk,g_log_user] = grad_R_sec(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);
   
    g_lse = zeros(Nr, 1);
    for k = 1:K   
        for l = 1:nF           
            g_lse = g_lse + weights(l,k) * (-(grad_sec_lk(:,l,k)));
        end
    end

    % === PENALTY GRADIENT ON rate_p(k) + Ck (per user only) ===
    g_penalty = zeros(Nr, 1);
    violation = Rmin - (rate_p + Ck);      % works if Rmin scalar or K×1
    mask = violation > 0;
    for k = 1:K
        if mask(k)
            g_rate_k = log(2) * rate_p(k) * g_log_user(:,k);
            g_penalty = g_penalty + (-2 * lambda_penalty * violation(k)) * g_rate_k;
        end
    end

    g = g_lse + g_penalty;
end

function h = hess_func(beta, v, L_node, E_node, alpha, s_param, K, nF, sigma2, Pw, AN_P_ratio, Rmin, lambda_penalty, Ck)
    Nr = length(beta);
    ln2 = log(2);
    noise_total = (sigma2 / Pw);
    alpha_pi = alpha(2:end);
    sum_alpha_pi = sum(alpha_pi);

    % 1. Get current rates and power terms
    [R_sec, rate_p, rate_c, T_k, T_l, I_k, I_l, Ak, Pk, Al, Pl] = get_Secrecy_matrix(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);

    % 2. Log-Sum-Exp weights (identical to grad_func)
    Z = -s_param * R_sec;
    m = max(Z, [], "all");
    weights = exp(Z - m);
    weights = weights / sum(weights, "all");

    % 3. Initialize accumulators
    g_final     = zeros(Nr, 1);
    h_sum_parts = zeros(Nr, 1);

    % Precompute violation for the penalty (on rate_p(k) + Ck)
    violation = Rmin - (rate_p + Ck);

    for k = 1:K
        % ===================== USER k PRECOMPUTATIONS (once per k) =====================
        grad_Pk = 2 * (L_node(k).V1 * beta) + 2 * L_node(k).V2;
        grad_Ak = 2 * (L_node(k).V1_AN * beta) + 2 * L_node(k).V2_AN;
        H_Pk_v  = 2 * (L_node(k).V1 * v);
        H_Ak_v  = 2 * (L_node(k).V1_AN * v);

        dPk = real(grad_Pk' * v);
        dAk = real(grad_Ak' * v);

        u_k      = (Ak(k) + noise_total) * grad_Pk - Pk(k) * grad_Ak;
        du_k     = dAk * grad_Pk + (Ak(k) + noise_total) * H_Pk_v - dPk * grad_Ak - Pk(k) * H_Ak_v;
        denom_k  = ln2 * T_k(k) * I_k(k);
        dIk      = (sum_alpha_pi - alpha_pi(k)) * dPk + AN_P_ratio * dAk;
        dTk      = alpha_pi(k) * dPk + dIk;
        ddenom_k = ln2 * (dTk * I_k(k) + T_k(k) * dIk);

        h_user_k = alpha_pi(k) * (du_k * denom_k - u_k * ddenom_k) / (denom_k^2);
        g_user_k = (alpha_pi(k) / denom_k) * u_k;

        % ===================== PENALTY HESSIAN-VECTOR (on rate_p(k) + Ck) =====================
        % rate_p(k) = linear SINR ratio → we need ∇rate and Hrate*v
        if violation(k) > 0
            g_rate_k  = ln2 * rate_p(k) * g_user_k;                    % ∇(linear rate)
            df_rate_v = real(g_rate_k' * v);

            % Correct second derivative of linear rate (derived from log2 rate)
            h_rate_k = ln2 * rate_p(k) * h_user_k ...
                     + (ln2)^2 * rate_p(k) * (real(g_user_k' * v))^2;

            % Quadratic penalty Hessian-vector product
            h_penalty_k = 2 * lambda_penalty * (df_rate_v * g_rate_k - violation(k) * h_rate_k);
            h_sum_parts = h_sum_parts + h_penalty_k;
        end

        % ===================== EAVESDROPPER LOOP (per l, for LSE part only) =====================
        for l = 1:nF
            % --- Precompute Eavesdropper l components ---
            grad_Pl = 2 * (E_node(l).V1_l * beta) + 2 * E_node(l).V2_l;
            grad_Al = 2 * (E_node(l).V1_AN_l * beta) + 2 * E_node(l).V2_AN_l;
            H_Pl_v  = 2 * (E_node(l).V1_l * v);
            H_Al_v  = 2 * (E_node(l).V1_AN_l * v);

            dPl = real(grad_Pl' * v);
            dAl = real(grad_Al' * v);

            u_l      = (Al(l) + noise_total) * grad_Pl - Pl(l) * grad_Al;
            du_l     = dAl * grad_Pl + (Al(l) + noise_total) * H_Pl_v - dPl * grad_Al - Pl(l) * H_Al_v;
            denom_l  = ln2 * T_l(l,k) * I_l(l,k);
            dIl      = (sum_alpha_pi - alpha_pi(k)) * dPl + AN_P_ratio * dAl;
            dTl      = alpha_pi(k) * dPl + dIl;
            ddenom_l = ln2 * (dTl * I_l(l,k) + T_l(l,k) * dIl);

            h_eav_l = alpha_pi(k) * (du_l * denom_l - u_l * ddenom_l) / (denom_l^2);
            g_eav_l = (alpha_pi(k) / denom_l) * u_l;

            % --- Component g_lk and h_lk for f_lk = -R_sec(l,k) ---
            g_lk = -(g_user_k - g_eav_l);
            h_lk = -(h_user_k - h_eav_l);

            df_lk_v = real(g_lk' * v);

            % Accumulate for LSE Hessian
            g_final     = g_final     + weights(l,k) * g_lk;
            h_sum_parts = h_sum_parts + weights(l,k) * (h_lk + s_param * df_lk_v * g_lk);
        end
    end

    % 4. Final LSE Hessian assembly
    % ehess = sum(w_i * (H_i v + s (g_i' v) g_i)) - s (grad' v) grad
    df_final_v = real(g_final' * v);
    h = h_sum_parts - s_param * df_final_v * g_final;
end
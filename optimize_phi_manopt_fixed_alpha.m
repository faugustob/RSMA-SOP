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
        options.maxiter = 10;         % Maximum outer iterations
        options.verbosity = 0;          % 2 shows summary, 3 shows detailed inner steps
        options.linesearch = @linesearch; % Trust-regions usually manages step size via the radius
        
        % Inner iteration control (Krylov steps)     
        options.maxinner = 60; 
        
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

function [g, grad_sec_lk] = grad_func(beta, L_node, E_node, alpha, s_param, K, nF, sigma2, Pw, AN_P_ratio, Rmin, lambda_penalty, Ck)
    Nr = length(beta);   
    [R_sec, rate_p, ~] = get_Secrecy_matrix(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);
    
    % 1. Log-Sum-Exp Weights for smoothing max(-R_sec)
    Z = -s_param * R_sec;
    m = max(Z, [], "all");
    weights = exp(Z - m);
    weights = weights / sum(weights, "all");
    
    % 2. Get Raw Gradients
    [grad_sec_lk, g_user_all] = grad_R_sec(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);
   
    % 3. LSE Gradient Component
    g_lse = zeros(Nr, 1);
    for k = 1:K   
        for l = 1:nF           
            g_lse = g_lse + weights(l,k) * (-(grad_sec_lk(:,l,k)));
        end
    end

    % 4. Penalty Gradient Component (Chain Rule: 2 * lambda * violation * -grad_R)
    g_penalty = zeros(Nr, 1);
    violation = Rmin - (rate_p + Ck); 
    mask = violation > 0;
    
    for k = 1:K
        if mask(k)
            % Corrected: The gradient of the penalty (Rmin - R)^2 is -2*penalty*grad_R
            % Note: g_user_all is already dR/dbeta
            g_penalty = g_penalty + (-2 * lambda_penalty * violation(k)) * g_user_all(:,k);
        end
    end
    
    g = g_lse + g_penalty;
end

function [grad_Rsec_lk, g_user_all] = grad_R_sec(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio)
    Nr = length(beta);
    ln2 = log(2);
    noise_total = (sigma2/Pw);
    alpha_pi = alpha(2:end);    
   
    [~,~,T_k,T_l,I_k,I_l,Ak,Pk,Al,Pl] = get_Secrecy_matrix(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);
    
    grad_Rsec_lk = zeros(Nr, nF, K);
    g_user_all   = zeros(Nr, K); 
   
    for k = 1:K
        % Legitimate User Gradients
        grad_Pk = 2 * (L_node(k).V1 * beta) + 2 * L_node(k).V2;           
        grad_Ak = 2 * (L_node(k).V1_AN * beta) + 2 * L_node(k).V2_AN;    
        
        % CORRECTED: Include AN_P_ratio in the Quotient Rule numerator
        u_k = (AN_P_ratio * Ak(k) + noise_total) * grad_Pk - Pk(k) * (AN_P_ratio * grad_Ak);
        denom_k = ln2 * T_k(k) * I_k(k);
        g_user_k = (alpha_pi(k) / denom_k) * u_k;
        
        g_user_all(:,k) = g_user_k; 
        
        for l = 1:nF
            % Eavesdropper Gradients
            grad_Pl = 2 * (E_node(l).V1_l * beta) + 2 * E_node(l).V2_l;          
            grad_Al = 2 * (E_node(l).V1_AN_l * beta) + 2 * E_node(l).V2_AN_l;                
            
            % CORRECTED: Include AN_P_ratio here as well
            u_l = (AN_P_ratio * Al(l) + noise_total) * grad_Pl - Pl(l) * (AN_P_ratio * grad_Al);
            denom_l = ln2 * T_l(l,k) * I_l(l,k);
            g_eav_l = (alpha_pi(k) / denom_l) * u_l;
            
            grad_Rsec_lk(:,l,k) = g_user_k - g_eav_l;
        end
    end
end

function h = hess_func(beta, v, L_node, E_node, alpha, s_param, K, nF, sigma2, Pw, AN_P_ratio, Rmin, lambda_penalty, Ck)
    Nr = length(beta);
    ln2 = log(2);
    noise_total = (sigma2 / Pw);
    alpha_pi = alpha(2:end);
    
    [R_sec,rate_p,T_k,T_l,I_k,I_l,Ak,Pk,Al,Pl] = get_Secrecy_matrix(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);
    
    Z = -s_param * R_sec;
    m = max(Z, [], "all");
    weights = exp(Z - m);
    weights = weights / sum(weights, "all");

    g_final = zeros(Nr, 1);
    h_sum_parts = zeros(Nr, 1);
    violation = Rmin - (rate_p + Ck);

    for k = 1:K
        % USER k PRECOMPUTATIONS
        grad_Pk = 2 * (L_node(k).V1 * beta) + 2 * L_node(k).V2;
        grad_Ak = 2 * (L_node(k).V1_AN * beta) + 2 * L_node(k).V2_AN;
        H_Pk_v  = 2 * (L_node(k).V1 * v);
        H_Ak_v  = 2 * (L_node(k).V1_AN * v);
        
        dPk = real(grad_Pk' * v);
        dAk = real(grad_Ak' * v);
        
        % CORRECTED: AN_P_ratio in Hessian-vector product components
        u_k  = (AN_P_ratio * Ak(k) + noise_total) * grad_Pk - Pk(k) * (AN_P_ratio * grad_Ak);
        du_k = (AN_P_ratio * dAk) * grad_Pk + (AN_P_ratio * Ak(k) + noise_total) * H_Pk_v ...
               - dPk * (AN_P_ratio * grad_Ak) - Pk(k) * (AN_P_ratio * H_Ak_v);
        
        dIk      = (sum(alpha_pi) - alpha_pi(k)) * dPk + AN_P_ratio * dAk;
        dTk      = alpha_pi(k) * dPk + dIk;
        denom_k  = ln2 * T_k(k) * I_k(k);
        ddenom_k = ln2 * (dTk * I_k(k) + T_k(k) * dIk);
        
        h_user_k = alpha_pi(k) * (du_k * denom_k - u_k * ddenom_k) / (denom_k^2);
        g_user_k = (alpha_pi(k) / denom_k) * u_k;

        % PENALTY HESSIAN
        if violation(k) > 0
            df_rate_v = real(g_user_k' * v);
            % H_penalty = 2 * lambda * (grad_R * grad_R' - violation * Hess_R)
            h_penalty_k = 2 * lambda_penalty * (df_rate_v * g_user_k - violation(k) * h_user_k);
            h_sum_parts = h_sum_parts + h_penalty_k;
        end

        for l = 1:nF
            grad_Pl = 2 * (E_node(l).V1_l * beta) + 2 * E_node(l).V2_l;
            grad_Al = 2 * (E_node(l).V1_AN_l * beta) + 2 * E_node(l).V2_AN_l;
            H_Pl_v  = 2 * (E_node(l).V1_l * v);
            H_Al_v  = 2 * (E_node(l).V1_AN_l * v);
            
            dPl = real(grad_Pl' * v);
            dAl = real(grad_Al' * v);
            
            % CORRECTED: AN_P_ratio for Eavesdropper
            u_l  = (AN_P_ratio * Al(l) + noise_total) * grad_Pl - Pl(l) * (AN_P_ratio * grad_Al);
            du_l = (AN_P_ratio * dAl) * grad_Pl + (AN_P_ratio * Al(l) + noise_total) * H_Pl_v ...
                   - dPl * (AN_P_ratio * grad_Al) - Pl(l) * (AN_P_ratio * H_Al_v);
            
            dIl      = (sum(alpha_pi) - alpha_pi(k)) * dPl + AN_P_ratio * dAl;
            dTl      = alpha_pi(k) * dPl + dIl;
            denom_l  = ln2 * T_l(l,k) * I_l(l,k);
            ddenom_l = ln2 * (dTl * I_l(l,k) + T_l(l,k) * dIl);
            
            h_eav_l = alpha_pi(k) * (du_l * denom_l - u_l * ddenom_l) / (denom_l^2);
            g_eav_l = (alpha_pi(k) / denom_l) * u_l;

            g_lk = -(g_user_k - g_eav_l);
            h_lk = -(h_user_k - h_eav_l);
            df_lk_v = real(g_lk' * v);

            g_final     = g_final     + weights(l,k) * g_lk;
            h_sum_parts = h_sum_parts + weights(l,k) * (h_lk + s_param * df_lk_v * g_lk);
        end
    end

    df_final_v = real(g_final' * v);
    h = h_sum_parts - s_param * df_final_v * g_final;
end
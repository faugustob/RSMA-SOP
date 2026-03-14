function [phi_St,cost_opt] = optimize_phi_manopt_fixed_alpha(Rmin,L_node,E_node,problem,b0,alpha,K, nF, sigma2, Pw, AN_P_ratio)

     
        %% ========================= MANOPT =========================
        




        s_param = 100; % Smoothing parameter for max-min
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

    
        problem.cost = @(b) cost_func(b, L_node, E_node, alpha, s_param, K, nF, sigma2, Pw, AN_P_ratio, Rmin, lambda_penalty);
        %problem = manoptAD(problem);

        problem.egrad = @(b) grad_func(b, L_node, E_node, alpha, s_param, K, nF, sigma2, Pw, AN_P_ratio, Rmin, lambda_penalty);

        problem.ehess = @(b,v) hess_func(b,v, L_node, E_node, alpha, s_param, K, nF, sigma2, Pw, AN_P_ratio, Rmin, lambda_penalty);
        % 
        % % Initial guess
        % b0 = manifold.rand();
        % checkgradient(problem);
        % checkhessian(problem);

        % Solve
       % fprintf('Starting Manifold Optimization...\n');
        % Define options structure
        options.tolgradnorm = 1e-10;    % Stop when the gradient norm is very small
        options.maxiter = 1e4;         % Maximum outer iterations
        options.verbosity = 0;          % 2 shows summary, 3 shows detailed inner steps
        options.linesearch = @linesearch; % Trust-regions usually manages step size via the radius
        
        % Inner iteration control (Krylov steps)
        % Since your Hessian is perfect, we can allow more inner iterations 
        % to solve the sub-problem more accurately.
        options.maxinner = Nr * 2; 
        
        % Execute with options
        [beta_opt, cost_opt, info] = trustregions(problem, b0, options);
        
        phi_St = angle(beta_opt).';
        
end



function f_val = cost_func(beta, L_node, E_node, alpha, s_param, K, nF, sigma2, Pw, AN_P_ratio, Rmin, lambda_penalty)

        %  X = [alpha,beta.'];
        % 
        % [~,sc_p_lk,~] = compute_sinr_sc_an_manopt(Pe,P,Q_j,nF,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,Rmin,h_rp,h_jq,h_e,zeta_k_St,Active_Gain_dB,X);

       
    [R_sec,~] = get_Secrecy_matrix(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);

    % Log-Sum-Exp smoothing for max(-R)
    %f_val =  (1/s_param) * log(sum(exp(-s_param * (R_sec)),'all'));
    Z = -s_param * R_sec;
    m = max(Z, [], 'all');                 % stabilization shift
    f_lse = (1/s_param) * ( m + log(sum(exp(Z - m), 'all')) );
    
    % --- NEW: Penalty Term ---
    % Quadratic penalty for violating Rmin
    penalty = lambda_penalty * sum(max(0, Rmin - R_sec).^2, 'all');
    penalty = 0;
    f_val = f_lse + penalty;
    
end

function [grad_Rsec_lk] = grad_R_sec(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio)
    Nr = length(beta);
    ln2 = log(2);
    noise_total = (sigma2/Pw);
   
    alpha_pi = alpha(2:end);    
   
    % 1. Compute all current rates and power terms
  
   
    [~,T_k,T_l,I_k,I_l,Ak,Pk,Al,Pl] = get_Secrecy_matrix(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);

    % 2. Log-Sum-Exp Weights    
   
    grad_Rsec_lk = zeros(Nr,nF,K);
   
    % 3. Gradient Accumulation
    g = zeros(Nr, 1);   % already correct
    
    for k = 1:K
        % ==================== FIXED GRADIENT LINES ====================
        grad_Pk = 2 * (L_node(k).V1 * beta) + 2 * L_node(k).V2;           
        grad_Ak = 2 * (L_node(k).V1_AN * beta) + 2 * L_node(k).V2_AN;    
        
                    
        g_user_k = (alpha_pi(k)/(ln2*T_k(k)*I_k(k))) * ((Ak(k)+noise_total)*grad_Pk-Pk(k)*grad_Ak);
        
        for l = 1:nF
            grad_Pl = 2 * (E_node(l).V1_l * beta) + 2 * E_node(l).V2_l;          
            grad_Al = 2 * (E_node(l).V1_AN_l * beta) + 2 * E_node(l).V2_AN_l;                
                        
            g_eav_l = (alpha_pi(k)/(ln2*T_l(l,k)*I_l(l,k))) * ((Al(l)+noise_total)*grad_Pl-Pl(l)*grad_Al);

            grad_Rsec_lk(:,l,k) = g_user_k - g_eav_l;
            
        end
    end
end

function [g, grad_sec_lk] = grad_func(beta, L_node, E_node, alpha, s_param, K, nF, sigma2, Pw, AN_P_ratio, Rmin, lambda_penalty)
    Nr = length(beta);   
    [R_sec,~] = get_Secrecy_matrix(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);

    % 2. Log-Sum-Exp Weights      
    Z = -s_param * R_sec;

    m = max(Z, [], "all");        % stabilization shift
    weights = exp(Z - m);
    
    weights = weights / sum(weights, "all");
    grad_sec_lk = grad_R_sec(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);
   
    g_lse = zeros(Nr, 1);
    for k = 1:K   
        for l = 1:nF           
            g_lse = g_lse + weights(l,k) * (-(grad_sec_lk(:,l,k)));
        end
    end

    % 2. Gradient of Penalty Part
    % d/dR [ lambda * (Rmin - R)^2 ] = -2 * lambda * (Rmin - R)
    g_penalty = zeros(Nr, 1);
    violation = Rmin - R_sec;
    mask = violation > 0; % Only penalize if R_sec < Rmin
    
    for k = 1:K
        for l = 1:nF
            if mask(l,k)
                g_penalty = g_penalty + (-2 * lambda_penalty * violation(l,k)) * grad_sec_lk(:,l,k);
            end
        end
    end

    g_penalty = zeros(Nr, 1);

    g = g_lse + g_penalty;
end

function h = hess_func(beta, v, L_node, E_node, alpha, s_param, K, nF, sigma2, Pw, AN_P_ratio, Rmin, lambda_penalty)
    Nr = length(beta);
    ln2 = log(2);
    noise_total = (sigma2/Pw);
    alpha_pi = alpha(2:end);
    sum_alpha_pi = sum(alpha_pi);
   
    % 1. Get current rates and power terms
    [R_sec, T_k, T_l, I_k, I_l, Ak, Pk, Al, Pl] = get_Secrecy_matrix(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);
    
    % 2. Compute Log-Sum-Exp Weights (identical to grad_func)
    % We are minimizing the LSE of (-R_sec)
        % 2. Log-Sum-Exp Weights      
    Z = -s_param * R_sec;

    m = max(Z, [], "all");        % stabilization shift
    weights = exp(Z - m);
    
    weights = weights / sum(weights, "all");
    
    % 3. Initialize accumulators
    g_final = zeros(Nr, 1);  % To store the total gradient
    h_sum_parts = zeros(Nr, 1); % To store sum(w_i * (H_i*v + s*(g_i'*v)*g_i))
    
    for k = 1:K
        % --- Precompute User K components ---
        grad_Pk = 2 * (L_node(k).V1 * beta) + 2 * L_node(k).V2;           
        grad_Ak = 2 * (L_node(k).V1_AN * beta) + 2 * L_node(k).V2_AN;
        H_Pk_v  = 2 * (L_node(k).V1 * v);
        H_Ak_v  = 2 * (L_node(k).V1_AN * v);
        
        % Directional derivatives of scalar power terms: real(grad' * v)
        dPk = real(grad_Pk' * v);
        dAk = real(grad_Ak' * v);
        
        % Intermediate terms for User Rate Gradient
        u_k = (Ak(k) + noise_total) * grad_Pk - Pk(k) * grad_Ak;
        du_k = dAk * grad_Pk + (Ak(k) + noise_total) * H_Pk_v - dPk * grad_Ak - Pk(k) * H_Ak_v;
        
        denom_k = ln2 * T_k(k) * I_k(k);
        dIk = (sum_alpha_pi - alpha_pi(k)) * dPk + AN_P_ratio * dAk;
        dTk = alpha_pi(k) * dPk + dIk;
        ddenom_k = ln2 * (dTk * I_k(k) + T_k(k) * dIk);
        
        % Individual Hessian-vector product for User K rate
        % Using quotient rule: (du*v - u*dv)/v^2
        h_user_k = alpha_pi(k) * (du_k * denom_k - u_k * ddenom_k) / (denom_k^2);
        g_user_k = (alpha_pi(k) / denom_k) * u_k;

        for l = 1:nF
            % --- Precompute Eavesdropper L components ---
            grad_Pl = 2 * (E_node(l).V1_l * beta) + 2 * E_node(l).V2_l;          
            grad_Al = 2 * (E_node(l).V1_AN_l * beta) + 2 * E_node(l).V2_AN_l;
            H_Pl_v  = 2 * (E_node(l).V1_l * v);
            H_Al_v  = 2 * (E_node(l).V1_AN_l * v);
            
            dPl = real(grad_Pl' * v);
            dAl = real(grad_Al' * v);
            
            u_l = (Al(l) + noise_total) * grad_Pl - Pl(l) * grad_Al;
            du_l = dAl * grad_Pl + (Al(l) + noise_total) * H_Pl_v - dPl * grad_Al - Pl(l) * H_Al_v;
            
            denom_l = ln2 * T_l(l,k) * I_l(l,k);
            dIl = (sum_alpha_pi - alpha_pi(k)) * dPl + AN_P_ratio * dAl;
            dTl = alpha_pi(k) * dPl + dIl;
            ddenom_l = ln2 * (dTl * I_l(l,k) + T_l(l,k) * dIl);
            
            h_eav_l = alpha_pi(k) * (du_l * denom_l - u_l * ddenom_l) / (denom_l^2);
            g_eav_l = (alpha_pi(k) / denom_l) * u_l;
            
            % --- Component Gradient and Hessian ---
            % Since objective component is f_lk = -R_sec:
            g_lk = -(g_user_k - g_eav_l);
            h_lk = -(h_user_k - h_eav_l);
            
            % Directional derivative of the component for the LSE second term
            df_lk_v = real(g_lk' * v);
            
            % Accumulate parts
            g_final = g_final + weights(l,k) * g_lk;
            h_sum_parts = h_sum_parts + weights(l,k) * (h_lk + s_param * df_lk_v * g_lk);
            % --- NEW: Penalty Hessian accumulation ---
            if (Rmin - R_sec(l,k)) > 0
                violation = Rmin - R_sec(l,k);
                % d/dbeta [ -2 * lambda * (Rmin - R) * grad_R ]
                % = 2 * lambda * (grad_R * grad_R')*v - 2 * lambda * (Rmin - R) * Hessian_R*v
                
                % Note: h_lk here is already Hessian of (-R), which is exactly what we need
                h_penalty_lk = 2 * lambda_penalty * (real((-g_lk)' * v) * (-g_lk) + violation * h_lk);
                h_penalty_lk = zeros(size(h_penalty_lk));
                h_sum_parts = h_sum_parts + h_penalty_lk;
            end
        end
    end
    
    % 4. Final LSE Hessian Assembly
    % ehess = sum(w_i * (H_i*v + s*(g_i'*v)*g_i)) - s*(grad'*v)*grad
    df_final_v = real(g_final' * v);
    h = h_sum_parts - s_param * df_final_v * g_final;
end
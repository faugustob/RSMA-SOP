function [alpha_final] = new_optimize_alpha_cvx_fixed_phi(sigma2,Pw,L_node, E_node, beta, K, nF, max_SCA)
    %% ========================= CONSTANTS =========================
    AN_P_ratio = 1;         
    noise_total = (sigma2/Pw);
    
    %% ========================= PRECOMPUTE CHANNELS =========================
    Pk = zeros(K, 1); Ak = zeros(K, 1);
    Pl = zeros(nF, 1); Al = zeros(nF, 1);
    for k = 1:K
        Pk(k) = real(beta' * L_node(k).V1 * beta + 2*real(beta' * L_node(k).V2) + L_node(k).term3);
        Ak(k) = real(beta' * L_node(k).V1_AN * beta + 2*real(beta' * L_node(k).V2_AN) + L_node(k).term3_AN);      
    end
    for l = 1:nF
        Pl(l) = real(beta' * E_node(l).V1_l * beta + 2*real(beta' * E_node(l).V2_l) + E_node(l).term3_l);
        Al(l) = real(beta' * E_node(l).V1_AN_l * beta + 2*real(beta' * E_node(l).V2_AN_l) + E_node(l).term3_AN_l);      
    end

    %% ========================= INITIALIZATION =========================
    % Start with a feasible uniform distribution
    alpha_prev = ones(K+1, 1) / (K+1); 

    %% ========================= SCA LOOP =========================
    for sca_iter = 1:max_SCA
        cvx_begin quiet
            cvx_solver mosek
            
            variable vecAlpha(K+1) nonnegative
            variable s_fake(nF, K) 
            
            % Expressions for readability
            alpha_c = vecAlpha(1);
            alpha_pi = vecAlpha(2:end);
            sum_alpha_pi = sum(alpha_pi);
            
            % Previous values for linearization
            alpha_prev_c = alpha_prev(1);
            alpha_prev_pi = alpha_prev(2:end);
            sum_alpha_prev_pi = sum(alpha_prev_pi);

            maximize( min(min(s_fake)) )
            
            subject to
                sum(vecAlpha) <= 1;
                vecAlpha >= 1e-3; % Small floor to prevent log(0)
                
                for k = 1:K
                    % --- User Rate (Concave) ---
                    % R_k = log2(1 + (alpha_k * Pk) / (Interference))
                    % We use the log_det or log-sum-inv forms, but for SCA, 
                    % we can simplify if we assume sum_alpha is the main interference.
                    
                    I_k = (sum_alpha_pi - alpha_pi(k)) * Pk(k) + AN_P_ratio * Ak(k) + noise_total;
                    T_k = alpha_pi(k) * Pk(k) + I_k;
                    
                    % Linearizing the Eavesdropper part for Secrecy
                    for l = 1:nF
                        I_lk = (sum_alpha_pi - alpha_pi(k)) * Pl(l) + AN_P_ratio * Al(l) + noise_total;
                        T_lk = alpha_pi(k) * Pl(l) + I_lk;
                        
                        % Previous values for the Eavesdropper log term
                        I_lk_prev = (sum_alpha_prev_pi - alpha_prev_pi(k)) * Pl(l) + AN_P_ratio * Al(l) + noise_total;
                        T_lk_prev = alpha_prev_pi(k) * Pl(l) + I_lk_prev;
                        
                        % --- Secrecy Rate Approximation ---
                        % Rate_k - Rate_eve_lk
                        % Use First-order Taylor: -log(T_lk) <= -log(T_lk_prev) - (T_lk - T_lk_prev)/T_lk_prev
                        R_user = log(T_k)/log(2); 
                        R_eve_approx = (log(T_lk_prev) + (T_lk - T_lk_prev)/T_lk_prev)/log(2) - log(I_lk)/log(2);
                        
                        s_fake(l,k) <= R_user - R_eve_approx;
                    end
                end
        cvx_end
        
        if contains(cvx_status, 'Solved')
            % Check for convergence
            if norm(vecAlpha - alpha_prev) < 1e-4
                alpha_prev = double(vecAlpha);
                break;
            end
            alpha_prev = double(vecAlpha);
        else
            fprintf('SCA failed at iter %d: %s\n', sca_iter, cvx_status);
            break;
        end
    end
    alpha_final = alpha_prev;
end
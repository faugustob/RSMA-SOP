function [phi_St] = optimize_phi_sca_fixed_alpha_backup(alpha, phi_St, phi_Sr, zeta_k_St, ...
              K, Nr, nF, Pe, P, Q_j, Plos, PLj, HB, HA, g_pq, Nsymb, ...
              reflect, h_rp, h_jq, h_e, delta_f, Active_Gain_dB, max_sca)
        cvx_clear;
        Rmin = 1e-8;
        Nr = length(phi_St);
        zeta_k_Sr = (10^(Active_Gain_dB/10)) - zeta_k_St;
        phase_St = exp(1j .* phi_St);
        phase_Sr = exp(1j .* phi_Sr);
        beta_St = sqrt(zeta_k_St) .* phase_St;
        beta_Sr = sqrt(zeta_k_Sr) .* phase_Sr;
        BW = delta_f;
        N0_dBm = -174;
        sigma2 = 10^((N0_dBm + 10*log10(BW) - 30)/10);
        % Define the scaling factor
        scaling_factor = 1 / sigma2;
        Pw_dBm = 46;
        Pw = 10^((Pw_dBm - 30)/10);
        AN_P_ratio = 1;          % Increase this (e.g. 5-10) if eavesdroppers are too strong
        noise = (sigma2/Pw)*scaling_factor;
        
        %% ========================= PRECOMPUTE CHANNELS =========================
        Nc_k_all = zeros(K,1);
        Nc_k_AN_all = zeros(K,1);
        Nc_l_all = zeros(nF,K);
        Nc_l_AN_all = zeros(nF,K);
        
        for k = 1:K
             reflect_coeff = reflect(k);
                beta_r = (reflect_coeff == 1) * beta_Sr + (reflect_coeff == -1) * beta_St;
        
            % Nc_k_all(k) = compute_SDR( ...
            % 0, Pe, P, Q_j, Plos(k,1), PLj(k,1), Nr, ...
            % HB(:,:,:,k), HA(:,:,:,:,k), g_pq(:,:,k), ...
            % beta_r, Nsymb, h_rp(:,:,k,1), h_jq(:,:,k), ...
            % h_e(:,k,1), 'loop')* scaling_factor;
            % 
            % Nc_k_all(k) = compute_SDR( ...
            % 0, Pe, P, Q_j, Plos(k,1), PLj(k,1), Nr, ...
            % HB(:,:,:,k), HA(:,:,:,:,k), g_pq(:,:,k), ...
            % beta_r, Nsymb, h_rp(:,:,k,1), h_jq(:,:,k), ...
            % h_e(:,k,1), 'vectorized')* scaling_factor;

             Nc_k_all(k) = compute_OTFS_static_channel( ...
            0, Pe, P, Q_j, Plos(k,1), PLj(k,1), Nr, ...
            HB(:,:,:,k), HA(:,:,:,:,k), g_pq(:,:,k), ...
            beta_r, Nsymb, h_rp(:,:,k,1), h_jq(:,:,k), ...
            h_e(:,k,1), 'vectorized')* scaling_factor;
    
            Nc_k_AN_all(k) = compute_OTFS_static_channel( ...
                0, Pe, P, Q_j, Plos(k,2), PLj(k,2), Nr, ...
                HB(:,:,:,k), HA(:,:,:,:,k), g_pq(:,:,k), ...
                beta_r, Nsymb, h_rp(:,:,k,2), h_jq(:,:,k), ...
                h_e(:,k,2), 'vectorized')* scaling_factor;;
        end
        for l = 1:nF
            for k = 1:K
                reflect_coeff = reflect(k);
                beta_r = (reflect_coeff == 1) * beta_Sr + (reflect_coeff == -1) * beta_St;
                Nc_l_all(l,k) = compute_OTFS_static_channel( ...
                    1, Pe, P, Q_j, Plos(K+l,1), PLj(K+l,1), Nr, ...
                    HB(:,:,:,K+l), HA(:,:,:,:,K+l), g_pq(:,:,K+l), ...
                    beta_r, Nsymb, h_rp(:,:,K+l,1), h_jq(:,:,K+l), ...
                    h_e(:,K+l,1), 'vectorized')* scaling_factor;
                Nc_l_AN_all(l,k) = compute_OTFS_static_channel( ...
                    1, Pe, P, Q_j, Plos(K+l,2), PLj(K+l,2), Nr, ...
                    HB(:,:,:,K+l), HA(:,:,:,:,K+l), g_pq(:,:,K+l), ...
                    beta_r, Nsymb, h_rp(:,:,K+l,2), h_jq(:,:,K+l), ...
                    h_e(:,K+l,2), 'vectorized')* scaling_factor;
            end
        end
        
        %% ========================= INITIALIZATION =========================
        alpha_pi = alpha(2:end);
        sum_alpha_pi = sum(alpha_pi);
      
        I_j_prev = zeros(K,1); 
        gamma_j_prev = zeros(K,1); 
        I_l_prev = zeros(nF,K); gamma_l_prev = zeros(nF,K);
        
        for k = 1:K
            I_j_prev(k) = (sum_alpha_pi - alpha_pi(k)) * Nc_k_all(k) + AN_P_ratio * Nc_k_AN_all(k);
            gamma_j_prev(k) = (alpha_pi(k) * Nc_k_all(k)) / (I_j_prev(k) + noise);
        
            
            for l = 1:nF
                I_l_prev(l,k) = (sum_alpha_pi - alpha_pi(k)) * Nc_l_all(l,k) + AN_P_ratio * Nc_l_AN_all(l,k);
                gamma_l_prev(l,k) = (alpha_pi(k) * Nc_l_all(l,k)) / (I_l_prev(l,k) + noise);
            end
        end
        gamma_j_prev = max(gamma_j_prev, 1e-6);
        gamma_l_prev = max(gamma_l_prev, 1e-6);
        
        %% ========================= SCA LOOP =========================
        for sca_iter = 1:max_sca
            cvx_begin quiet
                cvx_solver mosek
                
                % 1. Define the Augmented Matrix Variable
               % variable W(Nr+1, Nr+1) complex hermitian
                variable beta_r(1,Nr) complex  
                variable gamma_j(K) nonnegative
                variable gamma_l(nF,K) nonnegative
                variable s_fake(nF,K)
                
                maximize( min(min(s_fake)) )
        
                subject to

                   % 2. SDR Constraints
                    % W == semidefinite(Nr+1); % Must be Positive Semidefinite
                    % W(Nr+1, Nr+1) == 1;      % The "1" in the augmented corner
                    % diag(W(1:Nr, 1:Nr)) == 1; % Equivalent to abs(beta_r) <= 1
                    
                    abs(beta_r) <=1;
                    
                    % ---------- USER-LEVEL CONSTRAINTS ----------
                    for k = 1:K    
                      
                       [V1, V2_gpu, term3_gpu] = compute_V(    0, Pe, P, Q_j, Plos(k,1), PLj(k,1), Nr, HB(:,:,:,k), HA(:,:,:,:,k), g_pq(:,:,k), ...
                                                       Nsymb, h_rp(:,:,k,1), h_jq(:,:,k), h_e(:,k,1));

                       [V1_AN, V2_AN_gpu, term3_AN_gpu] = compute_V(    0, Pe, P, Q_j, Plos(k,2), PLj(k,2), Nr, HB(:,:,:,k), HA(:,:,:,:,k), g_pq(:,:,k), ...
                                                       Nsymb, h_rp(:,:,k,2), h_jq(:,:,k), h_e(:,k,2));
                        


                       V2 = gather(V2_gpu)*scaling_factor;

                       V2_AN = gather(V2_AN_gpu)*scaling_factor;
                       
                       term3 = gather(term3_gpu)*scaling_factor;

                       term3_AN = gather(term3_AN_gpu)*scaling_factor;

                       % V_aug = [V1, V2'; V2, term3]
                        % V_k_aug = [V1, V2.'; conj(V2), term3]*scaling_factor;
                        % 
                        % V_AN_k_aug = [V1_AN, V2_AN.'; conj(V2_AN), term3_AN]*scaling_factor;
    
   
                       % 4. Quadratic expression in beta_r              
                        Nc_k =    quad_form(beta_r', V1)  + 2*real(V2 * beta_r.') + term3;
                        AN_k =    quad_form(beta_r', V1_AN) + 2*real(V2_AN * beta_r.') + term3_AN;

                        % Nc_k = real(trace(V_k_aug * W));
                        % AN_k = real(trace(V_AN_k_aug * W));
        
                        S_j = alpha_pi(k) * Nc_k;
                        I_j = (sum_alpha_pi - alpha_pi(k)) * Nc_k+ AN_P_ratio * AN_k;
        
                       % Standard SCA linearization logic continues...
                        % (Note: gamma_j(k)*(I_j + noise) <= S_j is still non-convex,
                        %  usually handled by linearizing the product)
                        gamma_j_prev(k)*(I_j + noise) + I_j_prev(k)*gamma_j(k) - gamma_j_prev(k)*I_j_prev(k) <= S_j;
                    end
        
                    % ---------- EAVESDROPPER / SECRECY CONSTRAINTS (SDR VERSION) ----------
            for l = 1:nF
                for k = 1:K
                    % 1. Get the raw channel components
                    [V1_lk, V2_lk_gpu, term3_lk_gpu] = compute_V(1, Pe, P, Q_j, Plos(K+l,1), PLj(K+l,1), Nr, ...
                                                HB(:,:,:,K+l), HA(:,:,:,:,K+l), g_pq(:,:,K+l), ...
                                                Nsymb, h_rp(:,:,K+l,1), h_jq(:,:,K+l), h_e(:,K+l,1));
                                                
                    [V1_AN_lk, V2_AN_lk_gpu, term3_AN_lk_gpu] = compute_V(1, Pe, P, Q_j, Plos(K+l,2), PLj(K+l,2), Nr, ...
                                                        HB(:,:,:,K+l), HA(:,:,:,:,K+l), g_pq(:,:,K+l), ...
                                                        Nsymb, h_rp(:,:,K+l,2), h_jq(:,:,K+l), h_e(:,K+l,2));

                    V2_lk = gather(V2_lk_gpu)*scaling_factor;
                    V2_AN_lk= gather(V2_AN_lk_gpu)*scaling_factor;
                    term3_lk = gather(term3_lk_gpu)*scaling_factor;
                    term3_AN_lk  = gather(term3_AN_lk_gpu)*scaling_factor;

                    % 2. Construct the Augmented Matrices for SDR
                    % These matrices contain all channel information (quadratic, linear, and constant)
                    % V_lk_aug = [V1_lk, V2_lk.'; conj(V2_lk), term3_lk]*scaling_factor;
                    % V_AN_lk_aug = [V1_AN_lk, V2_AN_lk.'; conj(V2_AN_lk), term3_AN_lk]*scaling_factor;
            
                    % 3. Quadratic expression in beta_r 
                    
                    Nc_lk =   quad_form(beta_r', V1_lk) + 2*real(V2_lk * beta_r.') + term3_lk;    
                    AN_lk =   quad_form(beta_r', V1_AN_lk) + 2*real(V2_AN_lk * beta_r.') + term3_AN_lk;
            
                    S_l = alpha(k+1) * Nc_lk;
                    I_l = (sum(alpha(2:end)) - alpha(k+1)) * Nc_lk + AN_P_ratio * AN_lk;
            
                    % 4. SCA Linearization: Upper bound on the product gamma_l * (I_l + noise)
                    % We use the identity: xy <= x_prev*y + y_prev*x - x_prev*y_prev
                    gamma_l_prev(l,k)*(I_l + noise) + I_l_prev(l,k)*gamma_l(l,k) - ...
                    gamma_l_prev(l,k)*I_l_prev(l,k) <= S_l;
            
                    % 5. Secrecy Rate Approximation
                    % Upper bound on log(1 + gamma_l) to ensure a lower bound on (Rate_j - Rate_l)
                    log_l_approx = log(1 + gamma_l_prev(l,k))/log(2) + ...
                                   (1/((1 + gamma_l_prev(l,k))*log(2))) * (gamma_l(l,k) - gamma_l_prev(l,k));
            
                    % Legitimate rate is concave (log is handled by CVX), eavesdropper rate is linearized
                    s_fake(l,k) <= log(1 + gamma_j(k))/log(2) - log_l_approx;
                end
            end
        
            cvx_end

            
        
            if ~strcmp(cvx_status,'Solved') && ~strcmp(cvx_status,'Inaccurate/Solved')
                fprintf('SCA failed at iter %d: %s\n', sca_iter, cvx_status);
                fprintf('Try increasing AN_P_ratio (e.g. to 5-10) or check channel strengths.\n');
                break;
            end

            % % 5. Extract beta_r from W (Rank-1 Approximation)
            % [V_eig, D_eig] = eig(double(W));
            % beta_full = V_eig(:, end) * sqrt(D_eig(end, end)); % Principal eigenvector
            % beta_r_res = beta_full(1:Nr) / beta_full(Nr+1);    % Normalize by the augmented '1'

            beta_r_res = double(beta_r);
        
            % ---------- UPDATE FOR NEXT ITERATION ----------

           for k = 1:K               
            
                Nc_k_all(k) = compute_OTFS_static_channel( ...
                0, Pe, P, Q_j, Plos(k,1), PLj(k,1), Nr, ...
                HB(:,:,:,k), HA(:,:,:,:,k), g_pq(:,:,k), ...
                beta_r_res, Nsymb, h_rp(:,:,k,1), h_jq(:,:,k), ...
                h_e(:,k,1), 'vectorized')* scaling_factor;
        
                Nc_k_AN_all(k) = compute_OTFS_static_channel( ...
                    0, Pe, P, Q_j, Plos(k,2), PLj(k,2), Nr, ...
                    HB(:,:,:,k), HA(:,:,:,:,k), g_pq(:,:,k), ...
                    beta_r_res, Nsymb, h_rp(:,:,k,2), h_jq(:,:,k), ...
                    h_e(:,k,2), 'vectorized')* scaling_factor;
            end
            for l = 1:nF
                for k = 1:K
               
                    Nc_l_all(l,k) = compute_OTFS_static_channel( ...
                        1, Pe, P, Q_j, Plos(K+l,1), PLj(K+l,1), Nr, ...
                        HB(:,:,:,K+l), HA(:,:,:,:,K+l), g_pq(:,:,K+l), ...
                        beta_r_res, Nsymb, h_rp(:,:,K+l,1), h_jq(:,:,K+l), ...
                        h_e(:,K+l,1), 'vectorized')* scaling_factor;
                    Nc_l_AN_all(l,k) = compute_OTFS_static_channel( ...
                        1, Pe, P, Q_j, Plos(K+l,2), PLj(K+l,2), Nr, ...
                        HB(:,:,:,K+l), HA(:,:,:,:,K+l), g_pq(:,:,K+l), ...
                        beta_r_res, Nsymb, h_rp(:,:,K+l,2), h_jq(:,:,K+l), ...
                        h_e(:,K+l,2), 'vectorized')* scaling_factor;
                end
            end
            
            for k = 1:K
                gamma_j_prev(k) = double(gamma_j(k));
                I_j_prev(k) = (sum_alpha_pi - alpha_pi(k)) * Nc_k_all(k) + AN_P_ratio * Nc_k_AN_all(k);
        
                for l = 1:nF
                    gamma_l_prev(l,k) = double(gamma_l(l,k));
                    I_l_prev(l,k) = (sum_alpha_pi - alpha_pi(k)) * Nc_l_all(l,k) + AN_P_ratio * Nc_l_AN_all(l,k);
                end
            end
        
        end
        phi_St = angle(double(beta_r_res));

        cvx_clear;
  end
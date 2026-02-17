function [phi_St, phi_Sr, zeta_k_St] = optimize_phi_sca_fixed_alpha_test(alpha, phi_St, phi_Sr, zeta_k_St, ...
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
        Pw_dBm = 46;
        Pw = 10^((Pw_dBm - 30)/10);
        AN_P_ratio = 1;          % Increase this (e.g. 5-10) if eavesdroppers are too strong
        noise = max(sigma2/Pw, 1e-10);
        
        %% ========================= PRECOMPUTE CHANNELS =========================
        Nc_k_all = zeros(K,1);
        Nc_k_AN_all = zeros(K,1);
        Nc_l_all = zeros(nF,K);
        Nc_l_AN_all = zeros(nF,K);
        
        for k = 1:K
             reflect_coeff = reflect(k);
                beta_r = (reflect_coeff == 1) * beta_Sr + (reflect_coeff == -1) * beta_St;
        
            Nc_k_all(k) = compute_OTFS_static_channel( ...
            0, Pe, P, Q_j, Plos(k,1), PLj(k,1), Nr, ...
            HB(:,:,:,k), HA(:,:,:,:,k), g_pq(:,:,k), ...
            beta_r, Nsymb, h_rp(:,:,k,1), h_jq(:,:,k), ...
            h_e(:,k,1), 'vectorized');
    
            Nc_k_AN_all(k) = compute_OTFS_static_channel( ...
                0, Pe, P, Q_j, Plos(k,2), PLj(k,2), Nr, ...
                HB(:,:,:,k), HA(:,:,:,:,k), g_pq(:,:,k), ...
                beta_r, Nsymb, h_rp(:,:,k,2), h_jq(:,:,k), ...
                h_e(:,k,2), 'vectorized');
        end
        for l = 1:nF
            for k = 1:K
                reflect_coeff = reflect(k);
                beta_r = (reflect_coeff == 1) * beta_Sr + (reflect_coeff == -1) * beta_St;
                Nc_l_all(l,k) = compute_OTFS_static_channel( ...
                    1, Pe, P, Q_j, Plos(K+l,1), PLj(K+l,1), Nr, ...
                    HB(:,:,:,K+l), HA(:,:,:,:,K+l), g_pq(:,:,K+l), ...
                    beta_r, Nsymb, h_rp(:,:,K+l,1), h_jq(:,:,K+l), ...
                    h_e(:,K+l,1), 'vectorized');
                Nc_l_AN_all(l,k) = compute_OTFS_static_channel( ...
                    1, Pe, P, Q_j, Plos(K+l,2), PLj(K+l,2), Nr, ...
                    HB(:,:,:,K+l), HA(:,:,:,:,K+l), g_pq(:,:,K+l), ...
                    beta_r, Nsymb, h_rp(:,:,K+l,2), h_jq(:,:,K+l), ...
                    h_e(:,K+l,2), 'vectorized');
            end
        end
        
        %% ========================= INITIALIZATION =========================
        alpha_pi_init = alpha(2:end);
        sum_pi_init = sum(alpha_pi_init);
        
        I_j_prev = zeros(K,1); I_c_prev = zeros(K,1);
        gamma_j_prev = zeros(K,1); gamma_c_prev = zeros(K,1);
        I_l_prev = zeros(nF,K); gamma_l_prev = zeros(nF,K);
        
        for k = 1:K
            I_j_prev(k) = (sum_pi_init - alpha_pi_init(k)) * Nc_k_all(k) + AN_P_ratio * Nc_k_AN_all(k);
            gamma_j_prev(k) = (alpha_pi_init(k) * Nc_k_all(k)) / (I_j_prev(k) + noise);
        
            
            for l = 1:nF
                I_l_prev(l,k) = (sum_pi_init - alpha_pi_init(k)) * Nc_l_all(l,k) + AN_P_ratio * Nc_l_AN_all(l,k);
                gamma_l_prev(l,k) = (alpha_pi_init(k) * Nc_l_all(l,k)) / (I_l_prev(l,k) + noise);
            end
        end
        gamma_j_prev = max(gamma_j_prev, 1e-6);
        gamma_l_prev = max(gamma_l_prev, 1e-6);
        
        %% ========================= SCA LOOP =========================
        for sca_iter = 1:max_sca
            cvx_begin quiet
                cvx_solver mosek
                
                variable gamma_j(K) nonnegative
                variable gamma_l(nF,K) nonnegative
                variable beta_r(Nr,1) complex          
                variable s_fake(nF,K)            % can be negative (we take max(0,.) later if needed)
        
                %maximize( (1/(nF*K)) * sum(sum(s_fake)) )
        
                maximize( min(min(s_fake)) )
        
                subject to
                   abs(beta_r) <= 1;
               
                    % ---------- USER-LEVEL CONSTRAINTS ----------
                    for k = 1:K    
                      
                       [V_1, V_2, term3] = compute_V(    0, Pe, P, Q_j, Plos(k,1), PLj(k,1), Nr, HB(:,:,:,k), HA(:,:,:,:,k), g_pq(:,:,k), ...
                                                       Nsymb, h_rp(:,:,k,1), h_jq(:,:,k), h_e(:,k,1));
                       [V_1_AN, V_2_AN, term3_AN] = compute_V(    0, Pe, P, Q_j, Plos(k,2), PLj(k,2), Nr, HB(:,:,:,k), HA(:,:,:,:,k), g_pq(:,:,k), ...
                                                       Nsymb, h_rp(:,:,k,2), h_jq(:,:,k), h_e(:,k,2));
    
   
                        Nc_k =    beta_r*V_1*beta_r' + 2*real(trace(V_2 * beta_r.')) + term3;
                        AN_k =    beta_r*V_1_AN*beta_r' + 2*real(trace(V_2_AN * beta_r.')) + term3_AN;
        
                        S_j = alpha_pi(k) * Nc_k;
                        I_j = (sum_alpha_pi - alpha_pi(k)) * Nc_k+ AN_P_ratio * AN_k;
        
                        % Correct SCA linearization (upper bound on bilinear)
                        gamma_j(k)<=S_j/(I_j+noise);
                    end
        
                    % ---------- EAVESDROPPER / SECRECY CONSTRAINTS ----------
                    for l = 1:nF
                        for k = 1:K
                           [V_1_lk, V_2_lk, term3_lk] = compute_V(    1, Pe, P, Q_j, Plos(K+l,1), PLj(K+l,1), Nr, HB(:,:,:,K+l), HA(:,:,:,:,K+l), g_pq(:,:,K+l), ...
                                                       Nsymb, h_rp(:,:,K+l,1), h_jq(:,:,K+l), h_e(:,K+l,1));
                           [V_1_AN_lk, V_2_AN_lk, term3_AN_lk] = compute_V(    1, Pe, P, Q_j, Plos(K+l,2), PLj(K+l,2), Nr, HB(:,:,:,K+l), HA(:,:,:,:,K+l), g_pq(:,:,K+l), ...
                                                       Nsymb, h_rp(:,:,K+l,2), h_jq(:,:,K+l), h_e(:,K+l,2));
    
   
                            Nc_lk =    beta_r*V_1_lk*beta_r' + 2*real(trace(V_2_lk * beta_r.')) + term3_lk;
    
                            AN_lk =    beta_r*V_1_AN_lk*beta_r' + 2*real(trace(V_2_AN_lk * beta_r.')) + term3_AN_lk;
                            S_l = alpha_pi(k) * Nc_lk;
                            I_l = (sum_alpha_pi - alpha_pi(k)) * Nc_lk + AN_P_ratio * AN_lk;
                            gamma_l(l,k)>=S_l/(I_l+noise);
        
                            % Corrected: upper bound on gamma_l (same direction as legitimate users)
                            gamma_l_prev(l,k)*I_l + I_l_prev(l,k)*gamma_l(l,k) - gamma_l_prev(l,k)*I_l_prev(l,k) + gamma_l(l,k)*noise <= S_l;
        
                            % First-order upper bound on log(1+gamma_l) → lower bound on secrecy
                            log_l_approx = log(1 + gamma_l_prev(l,k))/log(2) + ...
                                           (1/((1 + gamma_l_prev(l,k))*log(2))) * (gamma_l(l,k) - gamma_l_prev(l,k));
        
                            s_fake(l,k) <= log(1 + gamma_j(k))/log(2) - log_l_approx;
                        end
                    end
        
            cvx_end
        
            if ~strcmp(cvx_status,'Solved') && ~strcmp(cvx_status,'Inaccurate/Solved')
                fprintf('SCA failed at iter %d: %s\n', sca_iter, cvx_status);
                fprintf('Try increasing AN_P_ratio (e.g. to 5-10) or check channel strengths.\n');
                break;
            end
        
            % ---------- UPDATE FOR NEXT ITERATION ----------
            
            for k = 1:K
                gamma_j_prev(k) = double(gamma_j(k));
                I_j_prev(k) = (sum_pi_val - alpha_pi_val(k)) * Nc_k_all(k) + AN_P_ratio * Nc_k_AN_all(k);
        
                for l = 1:nF
                    gamma_l_prev(l,k) = double(gamma_l(l,k));
                    I_l_prev(l,k) = (sum_pi_val - alpha_pi_val(k)) * Nc_l_all(l,k) + AN_P_ratio * Nc_l_AN_all(l,k);
                end
            end
        
        end
        phi_st = angle(doble(beta_r));
        cvx_clear;
  end
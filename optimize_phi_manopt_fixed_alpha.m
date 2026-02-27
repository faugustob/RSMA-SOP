function [phi_St] = optimize_phi_manopt_fixed_alpha(alpha, phi_St, phi_Sr, zeta_k_St, ...
              K, Nr, nF, Pe, P, Q_j, Plos, PLj, HB, HA, g_pq, Nsymb, ...
              reflect, h_rp, h_jq, h_e, delta_f, Active_Gain_dB)

        Rmin = 1e-8;
        Nr = length(phi_St);
        BW = delta_f;
        N0_dBm = -174;
        sigma2 = 10^((N0_dBm + 10*log10(BW) - 30)/10);
        % Define the scaling factor
        Pw_dBm = 46;
        Pw = 10^((Pw_dBm - 30)/10);
        AN_P_ratio = 1;         % Increase this (e.g. 5-10) if eavesdroppers are too strong      
        
               
        %% ========================= INITIALIZATION =========================
        

            % ---------- CHANNELS ----------
                for k = 1:K    
                  
                   [V1, V2_gpu, term3_gpu] = compute_V(    0, Pe, P, Q_j, Plos(k,1), PLj(k,1), Nr, HB(:,:,:,k), HA(:,:,:,:,k), g_pq(:,:,k), ...
                                                   Nsymb, h_rp(:,:,k,1), h_jq(:,:,k), h_e(:,k,1));

                   [V1_AN, V2_AN_gpu, term3_AN_gpu] = compute_V(    0, Pe, P, Q_j, Plos(k,2), PLj(k,2), Nr, HB(:,:,:,k), HA(:,:,:,:,k), g_pq(:,:,k), ...
                                                   Nsymb, h_rp(:,:,k,2), h_jq(:,:,k), h_e(:,k,2));
                    

                   L_node(k).V1 = V1;
                   L_node(k).V1_AN = V1_AN;

                   L_node(k).V2 = gather(V2_gpu.');
                   L_node(k).V2_AN = gather(V2_AN_gpu.');
                   
                   L_node(k).term3 = gather(term3_gpu);
                   L_node(k).term3_AN = gather(term3_AN_gpu);

                 % Nc_lk =   quad_form(beta_r', V1) + 2*real(V2 * beta_r.') + term3;    
                 % AN_lk =   quad_form(beta_r', V1_AN) + 2*real(V2_AN * beta_r.') + term3_AN;
    
                    
                end
        
                    % ---------- EAVESDROPPER / SECRECY CONSTRAINTS (SDR VERSION) ----------
                for l = 1:nF
                        % 1. Get the raw channel components
                        [V1_l, V2_l_gpu, term3_l_gpu] = compute_V(1, Pe, P, Q_j, Plos(K+l,1), PLj(K+l,1), Nr, ...
                                                    HB(:,:,:,K+l), HA(:,:,:,:,K+l), g_pq(:,:,K+l), ...
                                                    Nsymb, h_rp(:,:,K+l,1), h_jq(:,:,K+l), h_e(:,K+l,1));
                                                    
                        [V1_AN_l, V2_AN_l_gpu, term3_AN_l_gpu] = compute_V(1, Pe, P, Q_j, Plos(K+l,2), PLj(K+l,2), Nr, ...
                                                            HB(:,:,:,K+l), HA(:,:,:,:,K+l), g_pq(:,:,K+l), ...
                                                            Nsymb, h_rp(:,:,K+l,2), h_jq(:,:,K+l), h_e(:,K+l,2));


                        E_node(l).V1_l = V1_l;
                        E_node(l).V1_AN_l = V1_AN_l;

                        E_node(l).V2_l = gather(V2_l_gpu.');
                        E_node(l).V2_AN_l= gather(V2_AN_l_gpu.');

                        E_node(l).term3_l = gather(term3_l_gpu);
                        E_node(l).term3_AN_l  = gather(term3_AN_l_gpu);
    
                                      
                        % 3. Quadratic expression in beta_r 
                        
                        % Nc_lk =   quad_form(beta_r', V1_lk) + 2*real(V2_lk * beta_r.') + term3_lk;    
                        % AN_lk =   quad_form(beta_r', V1_AN_lk) + 2*real(V2_AN_lk * beta_r.') + term3_AN_lk;   
                      
                           
                 end
                
       

      
        
        %% ========================= MANOPT =========================
        manifold = complexcirclefactory(Nr,1);
        problem.M = manifold;

        s_param = 20; % Smoothing parameter for max-min

       


        problem.cost = @(b) cost_func(b, L_node, E_node, alpha, s_param, K, nF, sigma2, Pw, AN_P_ratio);
        problem = manoptAD(problem);

        % problem.egrad = @(b) grad_func(b, L_node, E_node, alpha, s_param, K, nF, sigma2, Pw, AN_P_ratio);
        % 
        % % Initial guess
        b0 = manifold.rand();
        % checkgradient(problem);

        % Solve
        fprintf('Starting Manifold Optimization...\n');
        [beta_opt, cost_opt] = trustregions(problem, b0);
        
        phi_St = angle(beta_opt).';
        
end

function R_sec = get_Secrecy_matrix(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_ratio)
 noise_total = (sigma2/Pw);
        alpha_pi = alpha(2:end);
        sum_alpha_pi = sum(alpha_pi);
        % 1. Compute all current rates and power terms
        R_sec = zeros(nF, K);
        T_k = zeros(K, 1); I_k = zeros(K, 1);
        T_l = zeros(nF, K); I_l = zeros(nF, K);

        for k = 1:K
            Pk = real(beta' * L_node(k).V1 * beta + 2*real(beta' * L_node(k).V2) + L_node(k).term3);   % you already fixed this earlier
            Ak = real(beta' * L_node(k).V1_AN * beta + 2*real(beta' * L_node(k).V2_AN) + L_node(k).term3_AN);
           
            I_k(k) = (sum_alpha_pi - alpha_pi(k)) * Pk + AN_ratio * Ak + noise_total;
            T_k(k) = alpha_pi(k) * Pk + I_k(k);
       
            for l = 1:nF
                Pl = real(beta' * E_node(l).V1_l * beta + 2*real(beta' * E_node(l).V2_l) + E_node(l).term3_l);
                Al = real(beta' * E_node(l).V1_AN_l * beta + 2*real(beta' * E_node(l).V2_AN_l) + E_node(l).term3_AN_l);
               
                I_cl =  (sum_alpha_pi - alpha_pi(k)) * Pl + AN_ratio * Al + noise_total;
                T_cl = alpha_pi(k) * Pl + I_cl;
               
                I_l(l,k) = (sum_alpha_pi - alpha_pi(k)) * Pl + AN_ratio * Al + noise_total;
                T_l(l,k) = alpha_pi(k) * Pl + I_l(l,k);
               
                R_sec(l,k) = log2(T_k(k)/I_k(k)) - log2(T_l(l,k)/I_l(l,k));
            end
        end
end

function f_val = cost_func(beta, L_node, E_node, alpha, s_param, K, nF, sigma2, Pw, AN_ratio)

        %  X = [alpha,beta.'];
        % 
        % [~,sc_p_lk,~] = compute_sinr_sc_an_manopt(Pe,P,Q_j,nF,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,Rmin,h_rp,h_jq,h_e,zeta_k_St,Active_Gain_dB,X);

       
    R_sec = get_Secrecy_matrix(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_ratio);
    % Log-Sum-Exp smoothing for max(-R)
    f_val =  (1/s_param) * log(sum(exp(-s_param * (R_sec)),'all'));
end

function g = grad_func(beta, L_node, E_node, alpha, s_param, K, nF, sigma2, Pw, AN_ratio)
    Nr = length(beta);
    ln2 = log(2);
    noise_total = (sigma2/Pw);
   
    alpha_pi = alpha(2:end);
    sum_alpha_pi = sum(alpha_pi);
   
    % 1. Compute all current rates and power terms
    % R_sec = zeros(nF, K);
    T_k = zeros(K, 1); I_k = zeros(K, 1);
    T_l = zeros(nF, K); I_l = zeros(nF, K);
   
    R_sec = get_Secrecy_matrix(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_ratio);

    % 2. Log-Sum-Exp Weights      
    weights = exp(-s_param * (R_sec));
    weights = weights / sum(weights,"all");
    
   
    % 3. Gradient Accumulation
    g = zeros(Nr, 1);   % already correct
    
    for k = 1:K
        % ==================== FIXED GRADIENT LINES ====================
        grad_Pk = 2 * (L_node(k).V1 * beta) + 2 * L_node(k).V2;          % ← removed the '
        grad_Ak = 2 * (L_node(k).V1_AN * beta) + 2 * L_node(k).V2_AN;    % ← removed the '
        
        grad_Ik = (sum_alpha_pi - alpha_pi(k)) * grad_Pk + AN_ratio * grad_Ak;
        grad_Tk = alpha_pi(k) * grad_Pk + grad_Ik;
        
        g_user_k = (1/ln2) * (grad_Tk/T_k(k) - grad_Ik/I_k(k));
        
        for l = 1:nF
            grad_Pl = 2 * (E_node(l).V1_l * beta) + 2 * E_node(l).V2_l;          % ← removed the '
            grad_Al = 2 * (E_node(l).V1_AN_l * beta) + 2 * E_node(l).V2_AN_l;    % ← removed the '
            
            grad_Il = (sum_alpha_pi - alpha_pi(k)) * grad_Pl + AN_ratio * grad_Al;
            grad_Tl = alpha_pi(k) * grad_Pl + grad_Il;
            
            g_eav_l = (1/ln2) * (grad_Tl/T_l(l,k) - grad_Il/I_l(l,k));
            
            g = g + weights(l,k) * (-(g_user_k - g_eav_l));
        end
    end
end
function [R_sec,T_k,T_l,I_k,I_l,Ak,Pk,Al,Pl] = get_Secrecy_matrix(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio)
        noise_total = (sigma2/Pw);
        alpha_pi = alpha(2:end);
        sum_alpha_pi = sum(alpha_pi);
        % 1. Compute all current rates and power terms
        R_sec = zeros(nF, K);
        T_k = zeros(K, 1); I_k = zeros(K, 1); Ak = zeros(K, 1); Pk = zeros(K, 1);
        T_l = zeros(nF, K); I_l = zeros(nF, K); Pl = zeros(nF,1); Al = zeros(nF,1);

         for k = 1:K
            Pk(k) = real(beta' * L_node(k).V1 * beta + 2*real(beta' * L_node(k).V2) + L_node(k).term3);   % you already fixed this earlier
            Ak(k) = real(beta' * L_node(k).V1_AN * beta + 2*real(beta' * L_node(k).V2_AN) + L_node(k).term3_AN);      
         end

         for l = 1:nF
             Pl(l) = real(beta' * E_node(l).V1_l * beta + 2*real(beta' * E_node(l).V2_l) + E_node(l).term3_l);
             Al(l) = real(beta' * E_node(l).V1_AN_l * beta + 2*real(beta' * E_node(l).V2_AN_l) + E_node(l).term3_AN_l);      
         end

        for k = 1:K           
           
            I_k(k) = (sum_alpha_pi - alpha_pi(k)) * Pk(k) + AN_P_ratio * Ak(k) + noise_total;
            T_k(k) = alpha_pi(k) * Pk(k) + I_k(k);
       
            for l = 1:nF                             
                          
                I_l(l,k) = (sum_alpha_pi - alpha_pi(k)) * Pl(l) + AN_P_ratio * Al(l) + noise_total;
                T_l(l,k) = alpha_pi(k) * Pl(l) + I_l(l,k);
               
                R_sec(l,k) = log2(T_k(k)/I_k(k)) - log2(T_l(l,k)/I_l(l,k));
            end
        end
end
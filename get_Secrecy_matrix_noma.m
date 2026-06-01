function [R_sec,rate_p,interf_power_coeff,T_k,T_l,I_k,I_l,Ak,Pk,Al,Pl] = get_Secrecy_matrix_noma(beta, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio)
        
        T_k = zeros(K, 1); I_k = zeros(K, 1); Ak = zeros(K, 1); Pk = zeros(K, 1);
        T_l = zeros(nF, K); I_l = zeros(nF, K); Pl = zeros(nF,1); Al = zeros(nF,1);
        interf_power_coeff=zeros(K,1);

        %Compute channel gains and order them.
         for k = 1:K
            Pk(k) = real(beta' * L_node(k).V1 * beta + 2*real(beta' * L_node(k).V2) + L_node(k).term3);   % you already fixed this earlier
            Ak(k) = real(beta' * L_node(k).V1_AN * beta + 2*real(beta' * L_node(k).V2_AN) + L_node(k).term3_AN);      
         end

         %Sort channels from weakest to strongest
        [~, order] = sort(Pk);
         alpha_sorted = alpha(order);
        
        noise_total = (sigma2/Pw);       
        
        % 1. Compute all current rates and power terms
        R_sec = zeros(nF, K);
      

         for l = 1:nF
             Pl(l) = real(beta' * E_node(l).V1_l * beta + 2*real(beta' * E_node(l).V2_l) + E_node(l).term3_l);
             Al(l) = real(beta' * E_node(l).V1_AN_l * beta + 2*real(beta' * E_node(l).V2_AN_l) + E_node(l).term3_AN_l);      
         end

        for j = 1:K     
           
            k = order(j);
            
            % Intra-cluster interference comes from all users stronger than the j-th user
            if j < K
                interf_power_coeff(k) = sum(alpha_sorted(j+1:K));
            else
                interf_power_coeff(k) = 0; % The strongest user decodes perfectly (SIC)
            end

            I_k(k) = interf_power_coeff(k) * Pk(k) + AN_P_ratio * Ak(k) + noise_total;
            T_k(k) = alpha(k) * Pk(k) + I_k(k);

            rate_p(k) = log2(T_k(k)./I_k(k));
       
            for l = 1:nF                             
                          
                I_l(l,k) = interf_power_coeff(k) * Pl(l) + AN_P_ratio * Al(l) + noise_total;
                T_l(l,k) = alpha(k) * Pl(l) + I_l(l,k);
               
                R_sec(l,k) = log2(T_k(k)/I_k(k)) - log2(T_l(l,k)/I_l(l,k));
            end
        end
end
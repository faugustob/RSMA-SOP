function [L_node,E_node] = compute_channels( K, Nr, nF, Pe, P, Q_j, Plos, PLj, HB, HA, g_pq, Nsymb, ...
              reflect, h_rp, h_jq, h_e,  Active_Gain_dB)
    for k = 1:K    
                  
       [V1, V2_gpu, term3_gpu] = compute_V(    0, Pe, P, Q_j, Plos(k,1), PLj(k,1), Nr, HB(:,:,:,k,1), HA(:,:,:,:,k,1), g_pq(:,:,k,1), ...
                                       Nsymb, h_rp(:,:,1), h_jq(:,:,k), h_e(:,k,1),'vectorized');   

       [V1_AN, V2_AN_gpu, term3_AN_gpu] = compute_V(    0, Pe, P, Q_j, Plos(k,2), PLj(k,2), Nr, HB(:,:,:,k,2), HA(:,:,:,:,k,2), g_pq(:,:,k,2), ...
                                       Nsymb, h_rp(:,:,2), h_jq(:,:,k), h_e(:,k,2),'vectorized');
        

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
                                        Nsymb, h_rp(:,:,1), h_jq(:,:,K+l), h_e(:,K+l,1),'vectorized');
                                        
            [V1_AN_l, V2_AN_l_gpu, term3_AN_l_gpu] = compute_V(1, Pe, P, Q_j, Plos(K+l,2), PLj(K+l,2), Nr, ...
                                                HB(:,:,:,K+l), HA(:,:,:,:,K+l), g_pq(:,:,K+l), ...
                                                Nsymb, h_rp(:,:,2), h_jq(:,:,K+l), h_e(:,K+l,2),'vectorized');


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
end
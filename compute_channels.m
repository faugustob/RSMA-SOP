function [L_node,E_node] = compute_channels( K, Nr, nF, Pe, P, Q_j, Plos, PLj, HB, HA, g_pq, Nsymb, ...
              reflect, h_rp, h_jq, h_e,  Active_Gain_dB)
    for k = 1:K    
                  
       [V1, V2_gpu, term3_gpu] = compute_V(    0, Pe, P, Q_j, Plos(k,1), PLj(k,1), Nr, HB(:,:,:,k,1), HA(:,:,:,:,k,1), g_pq(:,:,k,1), ...
                                       Nsymb, h_rp(:,:,1), h_jq(:,:,k,1), h_e(:,k,1),'vectorized');   

       [V1_AN, V2_AN_gpu, term3_AN_gpu] = compute_V(    0, Pe, P, Q_j, Plos(k,2), PLj(k,2), Nr, HB(:,:,:,k,2), HA(:,:,:,:,k,2), g_pq(:,:,k,2), ...
                                       Nsymb, h_rp(:,:,2), h_jq(:,:,k,2), h_e(:,k,2),'vectorized');
        

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
                                        Nsymb, h_rp(:,:,1), h_jq(:,:,K+l,1), h_e(:,K+l,1),'vectorized');
                                        
            [V1_AN_l, V2_AN_l_gpu, term3_AN_l_gpu] = compute_V(1, Pe, P, Q_j, Plos(K+l,2), PLj(K+l,2), Nr, ...
                                                HB(:,:,:,K+l), HA(:,:,:,:,K+l), g_pq(:,:,K+l), ...
                                                Nsymb, h_rp(:,:,2), h_jq(:,:,K+l,2), h_e(:,K+l,2),'vectorized');


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

function [V_1, V_2, term3] = compute_V(I,Pe,P,Q_j,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,h_rp,h_jq,h_e, mode)
    % mode: 'loop' or 'vectorized'
    if nargin < 15, mode = 'vectorized'; end % Default to fast mode
    
    K = P * Q_j;
    index = @(p,q) (p-1)*Q_j + q;

    switch lower(mode)
        case 'loop'
            %% ===================== Original Loop-Based Computation =====================
            Gram_1 = zeros(K,K);
            for p = 1:P
                for q = 1:Q_j
                    k = index(p,q);
                    for pp = 1:P
                        for qp = 1:Q_j
                            l = index(pp,qp);
                            Gram_1(k,l) = trace(HA(:,:,p,q)' * HA(:,:,pp,qp));
                        end
                    end
                end
            end
            
            M_1 = zeros(K,Nr);
            for r = 1:Nr
                for p = 1:P
                    for q = 1:Q_j
                        k = index(p,q);
                        M_1(k,r) = h_rp(r,p)*h_jq(r,q)*g_pq(p,q);
                    end
                end
            end
            
            V_1 = zeros(Nr,Nr);
            for r = 1:Nr
                for s = 1:Nr
                    V_1(r,s) = (M_1(:,r)' * Gram_1 * M_1(:,s));
                end
            end
            V_1 = PLj*V_1;
            
            B1 = zeros(Nsymb,Nsymb);
            for u = 1:Pe
                B1 = B1 + I *h_e(u) * HB(:,:,u);
            end
            term3 =  Plos * norm(B1, 'fro').^2;    
            
            G_BA = zeros(Pe, K);
            for u = 1:Pe
                for p = 1:P
                    for q = 1:Q_j
                        k = index(p,q);
                        G_BA(u,k) = trace(HB(:,:,u)' * HA(:,:,p,q));
                    end
                end
            end
            
            V_2 = sqrt(PLj)*sqrt(Plos)*(I * conj(h_e(:)))' * G_BA * M_1 ;

        case 'vectorized'
            %% ===================== Optimized Vectorized Computation =====================
            HA_flat = reshape(HA, [], K); 
            % Gram_1 is A' * A, so we force it to be real if it's supposed to be
            Gram_1 = HA_flat' * HA_flat;
        
            [Q_idx, P_idx] = meshgrid(1:Q_j, 1:P);
            P_idx = P_idx(:); Q_idx = Q_idx(:);
            g_indices = sub2ind([P, Q_j], P_idx, Q_idx);
            
            M_1 = h_rp(:, P_idx).' .* h_jq(:, Q_idx).' .* g_pq(g_indices);
        
            % Force V_1 to be real to clean up 10^-37 noise
            V_1 = PLj * M_1' * Gram_1 * M_1; 
        
            HB_flat = reshape(HB, [], Pe);
            b_coeffs = I * h_e(:);
            B1_vec = HB_flat * b_coeffs;
            
            % term3 is a norm squared, it MUST be real
            term3 = Plos * real(B1_vec' * B1_vec); 
            
            G_BA = HB_flat' * HA_flat;
            V_2 = sqrt(PLj * Plos) * (conj(b_coeffs)' * G_BA * M_1);

        otherwise
            error('Invalid mode. Choose ''loop'' or ''vectorized''.');
    end
end
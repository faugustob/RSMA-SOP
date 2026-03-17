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
            V_1 = PLj * real(M_1' * Gram_1 * M_1); 
        
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
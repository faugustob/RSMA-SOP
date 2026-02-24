function [V_1, V_2, term3] = compute_V(I,Pe,P,Q_j,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,h_rp,h_jq,h_e)


%% ===================== Main Computation Loop =====================
% ------------------ Construct Term1 ------------------


    K = P*Q_j;
    index = @(p,q) (p-1)*Q_j + q;
    
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

    
    % ------------------ Construct Term2 ------------------
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

    M_2 = zeros(K, Nr);

    for r = 1:Nr
        for p = 1:P
            for q = 1:Q_j
                k = index(p,q);
                M_2(k,r) = h_rp(r,p) * h_jq(r,q) * g_pq(p,q);
            end
        end
    end

    b = I * h_e(:);   % Pe x 1

    V_2 = sqrt(PLj)*sqrt(Plos)*b' * G_BA * M_2 ;

end

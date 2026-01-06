function [Nc] = compute_OTFS_static_channel(I,Pe,P,Q_j,Plos,PLj,Nr,HB,HA,g_pq,beta_r,Nsymb,h_rp,h_jq,h_e)


%% ===================== Main Computation Loop =====================
% ------------------ Construct Term1 ------------------
    term1 = zeros(Nsymb,Nsymb);
    for r = 1:Nr
        beta = beta_r(r);
        for p = 1:P
            for q = 1:Q_j
                term1 = term1 + beta * h_rp(r,p) * h_jq(r,q) * g_pq(p,q) * HA(:,:,p,q);
            end
        end
    end
    term1 = sqrt(PLj)*term1;
    
    % ------------------ Construct Term2 ------------------
    B1 = zeros(Nsymb,Nsymb);
    for u = 1:Pe
        B1 = B1 + h_e(u) * HB(:,:,u);
    end
    term2 = I * sqrt(Plos) * B1;    

    Nc = norm(term1 + term2,'fro').^2;
end
    
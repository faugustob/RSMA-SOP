function Nc = compute_OTFS_static_channel( ...
    I,Pe,P,Q_j,Plos,PLj,Nr,HB,HA,g_pq, ...
    beta_r,Nsymb,h_rp,h_jq,h_e,method)

        if nargin < 16
            method = 'loop';   % default behavior
        end
        
        switch lower(method)
            case 'loop'
                Nc = compute_loop( ...
                    I,Pe,P,Q_j,Plos,PLj,Nr,HB,HA,g_pq, ...
                    beta_r,Nsymb,h_rp,h_jq,h_e);
        
            case 'vectorized'
                Nc = compute_vectorized( ...
                    I,Pe,P,Q_j,Plos,PLj,Nr,HB,HA,g_pq, ...
                    beta_r,Nsymb,h_rp,h_jq,h_e);
        
            otherwise
                error('Unknown method. Use ''loop'' or ''vectorized''.');
        end

end


function [Nc] = compute_loop(I,Pe,P,Q_j,Plos,PLj,Nr,HB,HA,g_pq,beta_r,Nsymb,h_rp,h_jq,h_e)

W = eye(Nsymb);   % matched filter / identity kernel

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
    % Heff = term1 + term2;
    % Nc = abs(trace(W' * Heff)).^2;

end
% function Nc = compute_vectorized( ...
%     I,Pe,P,Q_j,Plos,PLj,Nr,HB,HA,g_pq, ...
%     beta_r,Nsymb,h_rp,h_jq,h_e)
% 
% 
% 
% 
%     %% -------- Term 1 --------
%     term1 = zeros(Nsymb,Nsymb);
% 
%     for r = 1:Nr
%         beta = beta_r(r);
%         Hpq = beta * (h_rp(r,:).' * h_jq(r,:)) .* g_pq;
%         term1 = term1 + sum(HA .* reshape(Hpq,1,1,P,Q_j), [3 4]);
%     end
% 
%     term1 = sqrt(PLj) * term1;
% 
%     %% -------- Term 2 --------
%     B1 = sum(HB .* reshape(h_e,1,1,[]), 3);
%     term2 = I * sqrt(Plos) * B1;
% 
%     %% -------- Output --------
%     Heff = term1 + term2;
%     Nc = real(Heff(:)' * Heff(:));
% 
% 
% end
% function Nc = compute_vectorized(I, Pe, P, Q_j, Plos, PLj, Nr, HB, HA, g_pq, beta_r, Nsymb, h_rp, h_jq, h_e)
% 
%     % -------- Optimized Term 1 --------
%     % 1. Pre-calculate the combined weights for all r, p, and q
%     % The loop was calculating: sum_r(beta_r * (h_rp' * h_jq) .* g_pq)
%     % This can be rewritten as: (h_rp' * (beta_r .* h_jq)) .* g_pq
%     combined_weights = (h_rp.' * (beta_r(:) .* h_jq)) .* g_pq;
% 
%     % 2. Perform weighted sum of HA using a single matrix-vector multiply
%     % Reshape HA to (Nsymb^2 x P*Q_j) and multiply by flattened weights
%     HA_reshaped = reshape(HA, Nsymb^2, []);
%     term1 = reshape(HA_reshaped * combined_weights(:), Nsymb, Nsymb);
%     term1 = sqrt(PLj) * term1;
% 
%     % -------- Optimized Term 2 --------
%     % Use matrix-vector multiplication for the weighted sum of HB
%     HB_reshaped = reshape(HB, Nsymb^2, []);
%     B1 = reshape(HB_reshaped * h_e(:), Nsymb, Nsymb);
%     term2 = (I * sqrt(Plos)) * B1;
% 
%     % -------- Output --------
%     Heff = term1 + term2;
% 
%     % Use norm squared for speed (equivalent to real(Heff(:)' * Heff(:)))
%     Nc = norm(Heff, 'fro')^2;
% end
function Nc = compute_vectorized(I, Pe, P, Q_j, Plos, PLj, Nr, HB, HA, g_pq, beta_r, Nsymb, h_rp, h_jq, h_e)
    % Move inputs to GPU (do this once outside the function if possible)
    HA_gpu = gpuArray(reshape(HA, Nsymb^2, []));
    HB_gpu = gpuArray(reshape(HB, Nsymb^2, []));
    
    % Term 1
    combined_weights = (h_rp.' * (beta_r(:) .* h_jq)) .* g_pq;
    w1_gpu = gpuArray(combined_weights(:));
    term1 = reshape(HA_gpu * w1_gpu, Nsymb, Nsymb);
    
    % Term 2
    w2_gpu = gpuArray(h_e(:));
    term2 = (I * sqrt(Plos)) * reshape(HB_gpu * w2_gpu, Nsymb, Nsymb);
    
    % Final Calculation
    Heff = (sqrt(PLj) * term1) + term2;
    Nc = gather(norm(Heff, 'fro')^2); % Move result back to CPU
end
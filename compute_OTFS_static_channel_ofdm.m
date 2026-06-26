function [Nc, ICI, Heff] = compute_OTFS_static_channel_ofdm( ...
    I,Pe,P,Q_j,Plos,PLj,Nr,HB,HA,g_pq, ...
    beta_r,Nsymb,h_rp,h_jq,h_e,method)
        if nargin < 16
            method = 'loop';   % default behavior
        end
        
        switch lower(method)
            case 'loop'
                [Nc, ICI, Heff] = compute_loop( ...
                    I,Pe,P,Q_j,Plos,PLj,Nr,HB,HA,g_pq, ...
                    beta_r,Nsymb,h_rp,h_jq,h_e);
        
            case 'vectorized'
                [Nc, ICI, Heff] = compute_vectorized( ...
                    I,Pe,P,Q_j,Plos,PLj,Nr,HB,HA,g_pq, ...
                    beta_r,Nsymb,h_rp,h_jq,h_e);
        
            otherwise
                error('Unknown method. Use ''loop'' or ''vectorized''.');
        end
end

function [Nc, ICI, Heff] = compute_loop(I,Pe,P,Q_j,Plos,PLj,Nr,HB,HA,g_pq,beta_r,Nsymb,h_rp,h_jq,h_e)
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
    
    Heff = term1 + term2;
    
    % Use sum(abs(...).^2) to avoid vector Frobenius norm error and boost speed
    total_norm_sq = sum(abs(Heff(:)).^2);
    diag_norm_sq  = sum(abs(diag(Heff)).^2);
    
    ICI = total_norm_sq - diag_norm_sq;
    Nc  = diag_norm_sq; 
end

function [Nc, ICI, Heff] = compute_vectorized(I, Pe, P, Q_j, Plos, PLj, Nr, HB, HA, g_pq, beta_r, Nsymb, h_rp, h_jq, h_e)
    % 1. Pre-process shapes on CPU first
    HA_proc = reshape(HA, Nsymb^2, []);
    HB_proc = reshape(HB, Nsymb^2, []);
    
    combined_weights = (h_rp.' * (beta_r(:) .* h_jq)) .* g_pq;
    w1 = combined_weights(:);
    w2 = h_e(:);
    
    % 2. Move to GPU ONLY if available
    if gpuDeviceCount > 0
        HA_proc = gpuArray(HA_proc);
        HB_proc = gpuArray(HB_proc);
        w1 = gpuArray(w1);
        w2 = gpuArray(w2);
        I = gpuArray(I);
    end
    
    % 3. Core Math 
    term1 = reshape(HA_proc * w1, Nsymb, Nsymb);
    term2 = (I * sqrt(Plos)) * reshape(HB_proc * w2, Nsymb, Nsymb);
    
    Heff = (sqrt(PLj) * term1) + term2;
    
    % 4. Final Calculation (Optimized for GPU/CPU consistency)
    total_norm_sq = sum(abs(Heff(:)).^2);
    diag_norm_sq  = sum(abs(diag(Heff)).^2);
    
    ICI_dev = total_norm_sq - diag_norm_sq;
    Nc_dev  = diag_norm_sq;
    
    % 5. Move back to CPU
    Nc  = gather(Nc_dev); 
    ICI = gather(ICI_dev);
end
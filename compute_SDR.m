function Nc = compute_SDR( ...
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


%% ===================== Main Computation Loop =====================
% ------------------ Construct Term1 ------------------
   [V_1, V_2, term3] = compute_V(I,Pe,P,Q_j,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,h_rp,h_jq,h_e);
    
   
   Nc =    beta_r*V_1*beta_r' + 2*real(trace(V_2 * beta_r.')) + term3;


end

function Nc = compute_vectorized(I, Pe, P, Q_j, Plos, PLj, Nr, HB, HA, g_pq, beta_r, Nsymb, h_rp, h_jq, h_e)
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
        % Also move constants if they are used in large matrix operations
        I = gpuArray(I);
    end
    
    % 3. Core Math (This code runs on whichever device the data is currently on)
    term1 = reshape(HA_proc * w1, Nsymb, Nsymb);
    term2 = (I * sqrt(Plos)) * reshape(HB_proc * w2, Nsymb, Nsymb);
    
    Heff = (sqrt(PLj) * term1) + term2;
    
    % 4. Final Calculation & Move back to CPU
    % norm(..., 'fro')^2 is equivalent to sum(abs(x(:)).^2) which is often faster
    result = norm(Heff, 'fro')^2;
    Nc = gather(result); 
end
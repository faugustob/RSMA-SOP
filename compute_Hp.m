function Hp = compute_Hp(tau, nu, M, N, T, Deltaf, method,OTFStype)
    
    if nargin < 8
        OTFStype = 'MC';
    end
    switch lower(OTFStype)
        case 'mc'
            Hp = compute_Hp_mc(tau, nu, M, N, T, Deltaf,method);
        case 'zakr'
            Hp = compute_Hp_zakr(tau, nu, M, N, T, Deltaf,method);
        otherwise
            error('Unknown method. Use ''loop'' or ''blocked''.');
    end 


end

 function Hp = compute_Hp_mc(tau, nu, M, N, T, Deltaf, method)
    % compute_Hp  Wrapper to compute Hp using either 'loop' or 'blocked'
    %
    % Usage:
    %   Hp = compute_Hp(tau, nu, M, N, T, Deltaf, 'loop');
    %   Hp = compute_Hp(tau, nu, M, N, T, Deltaf, 'blocked');
    %
    % The file also contains:
    %   - Hp_loop(tau,nu,M,N,T,Deltaf)
    %   - Hp_blocked(tau,nu,M,N,T,Deltaf)
    
    if nargin < 7
        method = 'loop';
    end
    
    switch lower(method)
        case 'loop'
            Hp = Hp_loop(tau, nu, M, N, T, Deltaf);
        case 'blocked'
            Hp = Hp_blocked(tau, nu, M, N, T, Deltaf);
        otherwise
            error('Unknown method. Use ''loop'' or ''blocked''.');
    end
    
end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 1) CLEAR NESTED-LOOP VERSION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function Hp = Hp_loop(tau, nu, M, N, T, Deltaf)
    % Hp_loop  Build Hp (MN x MN) using straightforward nested loops
    
    MN = M * N;
    Hp = complex(zeros(MN, MN));
    
    nvec = 0:(N-1);
    mvec = 0:(M-1);
    mprimevec = 0:(M-1);
    j2pi = 1j * 2*pi;
    factor_common_exp = exp(-j2pi * nu * tau);
    
    % Prebuild m,m' matrices (MxM)
    [Mmat, Mprime] = meshgrid(mvec, mprimevec);
    delta = Mprime - Mmat;
    
    for kp = 0:(N-1)
        for k = 0:(N-1)
    
            % n-sum
            alpha_n = (kp-k)/N - nu/Deltaf;
            S_n = sum(exp(-j2pi * (nvec.' * alpha_n))); % scalar
    
            for lp = 0:(M-1)
                for l = 0:(M-1)
    
                    r = lp + kp*M + 1;
                    c = l  + k *M + 1;
    
                    % ---- h_{P,1} ----
                    pref1 = factor_common_exp * (1/M) * (1 - tau/T) * (1/N) * S_n;
                    term_exp1 = exp(1j*pi * ((tau/T + 1) * (delta + nu/Deltaf)));
                    phase1 = exp(1j*2*pi * ( (Mmat*lp)/M - (Mprime*l)/M - Mprime*(tau/T) ));
                    sinc1 = sinc((delta + nu/Deltaf)*(1 - tau/T));
                    hP1 = pref1 * sum( (term_exp1 .* phase1 .* sinc1), 'all' );
    
                    % ---- h_{P,2} ----
                    pref2 = factor_common_exp * exp(-j2pi*(k/N)) * (1/M) * (tau/T) * (1/N) * S_n;
                    term_exp2 = exp(1j*2*pi * ( (Mmat*lp)/M - (Mprime*l)/M - Mprime*(tau/T) ));
                    phase2 = exp(1j*pi*(tau/T)*(delta + nu/Deltaf));
                    sinc2 = sinc((delta + nu/Deltaf)*(tau/T));
                    hP2 = pref2 * sum( (term_exp2 .* phase2 .* sinc2), 'all' );
    
                    Hp(r,c) = hP1 + hP2;
                end
            end
    
        end
    end
    
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 2) BLOCK/VECTORIZED VERSION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function Hp = Hp_blocked(tau, nu, M, N, T, Deltaf)
    % Hp_blocked  Build Hp (MN x MN) using block computation for each (k',k)
    
    MN = M*N;
    Hp = complex(zeros(MN, MN));
    
    nvec = 0:(N-1);
    mvec = 0:(M-1);
    mprimevec = 0:(M-1);
    
    j2pi = 1j*2*pi;
    factor_common_exp = exp(-j2pi * nu * tau);
    
    % Preconstruct m,m' matrices
    [Mmat, Mprime] = meshgrid(mvec, mprimevec);
    Mmat = Mmat.';       % Mmat(i,j)=m_i
    Mprime = Mprime.';   % Mprime(i,j)=m'_j
    delta = Mprime - Mmat;
    
    % phase parts that do not depend on l,l'
    phase_base1 = exp(1j*pi * ((tau/T + 1) * (delta + nu/Deltaf)));
    phase_base2 = exp(1j*pi * ((tau/T)     * (delta + nu/Deltaf)));
    
    sinc1 = sinc((delta + nu/Deltaf) * (1 - tau/T));
    sinc2 = sinc((delta + nu/Deltaf) * (tau/T));
    
    A1 = phase_base1 .* sinc1;   % MxM
    A2 = phase_base2 .* sinc2;   % MxM
    
    % Precompute vector forms
    exp_mprime_tau = exp(-j2pi * (mprimevec.' * (tau/T)));    % Mx1
    
    for kp = 0:(N-1)
        for k = 0:(N-1)
    
            % --- n-sum ---
            alpha_n = (kp-k)/N - nu/Deltaf;
            S_n = sum(exp(-j2pi * (nvec.' * alpha_n)));
    
            pref1 = factor_common_exp * (1/M) * (1 - tau/T) * (1/N) * S_n;
            pref2 = factor_common_exp * exp(-j2pi*(k/N)) * (1/M) * (tau/T) * (1/N) * S_n;
    
            % --- Build MxM block for all (l',l) ---
    
            % exp(j2π m l'/M) matrix: rows m, cols l'
            lpvec = 0:(M-1);
            lvec  = (0:(M-1)).'; % column
    
            exp_m_lp = exp(1j*2*pi*((mvec.') * lpvec / M));   % MxM (m x l')
            exp_mprime_l = exp(-1j*2*pi*((mprimevec.') * lvec.' / M)); % MxM (m' x l)
    
            % Step 1: C1 = A1.' * exp_m_lp
            C1 = A1.' * exp_m_lp;   % (m' x l')
            C2 = A2.' * exp_m_lp;
    
            % Multiply by exp(-j2π m' τ/T)
            C1 = C1 .* (exp_mprime_tau * ones(1,M));
            C2 = C2 .* (exp_mprime_tau * ones(1,M));
    
            % Combine with exp(-j2π m' l/M) and sum over m'
            Block1 = (C1.') * exp_mprime_l;  % (l' x l)
            Block2 = (C2.') * exp_mprime_l;
    
            Block_total = pref1 * Block1 + pref2 * Block2;
    
            % Insert block
            row_start = kp*M + 1;
            col_start = k *M + 1;
    
            Hp(row_start:row_start+M-1, col_start:col_start+M-1) = Block_total;
    
        end
    end
    
    end
function Hp = compute_Hp_zakr(tau, nu, M, N, T, Deltaf, method)
% compute_Hp  Wrapper to compute Hp using either 'loop' or 'blocked'
tau_frac = tau - floor(tau / T) * T;
if nargin < 7
    method = 'loop';
end

switch lower(method)
    case 'loop'
        Hp = Hp_loop_zakr(tau, nu, M, N, T, Deltaf);
    case 'blocked'
        Hp = Hp_blocked_zakr(tau, nu, M, N, T, Deltaf);
    otherwise
        error('Unknown method. Use ''loop'' or ''blocked''.');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) NESTED-LOOP VERSION (Direct implementation of formula)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Hp = Hp_loop_zakr(tau, nu, M, N, T, Deltaf)
    MN = M * N;
    Hp = complex(zeros(MN, MN));
    n_vec = 0:N-1;
    m_vec = 0:M-1;
    j2pi = 1i * 2 * pi;

    for kp = 0:N-1
        for k = 0:N-1
            % Term C: Summation over n (Depends only on k' and k)
            % 1/N * sum_{n=0}^{N-1} e^{-j 2pi n T ((k'-k)Delta_f/N - nu)}
            phi_n = T * (((kp - k) * Deltaf / N) - nu);
            val_C = (1/N) * sum(exp(-j2pi * n_vec * phi_n));
            
            for lp = 0:M-1
                % Leading terms (Depends on l' and k')
                term1 = exp(j2pi * (nu / Deltaf) * (lp / M - tau / T));
                term2 = exp(j2pi * (kp / N - nu / Deltaf) * floor(lp / M - tau / T));
                
                for l = 0:M-1
                    % Term D: Summation over m (Depends on l' and l)
                    % 1/M * sum_{m=0}^{M-1} e^{j 2pi m Delta_f ((l'-l)T/M - tau)}
                    phi_m = Deltaf * (((lp - l) * T / M) - tau);
                    val_D = (1/M) * sum(exp(j2pi * m_vec * phi_m));
                    
                    % Map to matrix indices (1-based)
                    row = lp + kp * M + 1;
                    col = l + k * M + 1;
                    Hp(row, col) = term1 * term2 * val_C * val_D;
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) BLOCK/VECTORIZED VERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Hp = Hp_blocked_zakr(tau, nu, M, N, T, Deltaf)
    MN = M * N;
    Hp = complex(zeros(MN, MN));
    n_vec = (0:N-1).'; % Column vector
    m_vec = 0:M-1;      % Row vector
    lp_vec = (0:M-1).'; % Column vector for rows
    l_vec = 0:M-1;      % Row vector for columns
    j2pi = 1i * 2 * pi;

    % Precompute Leading Terms for all lp (M x N matrix)
    % Each column represents a different kp, each row a different lp
    [KP, LP] = meshgrid(0:N-1, 0:M-1);
    Term1 = exp(j2pi * (nu / Deltaf) * (LP / M - tau / T));
    Term2 = exp(j2pi * (KP / N - nu / Deltaf) .* floor(LP / M - tau / T));
    LeadingTerms = Term1 .* Term2; % (M x N)

    % Precompute Term D matrix (M x M)
    % This depends on (lp - l) and is independent of k, kp
    % Uses implicit expansion (broadcasting) for (lp_vec - l_vec)
    phi_m_mat = Deltaf * (((lp_vec - l_vec) * T / M) - tau); % (M x M)
    % We need to sum over m for each element of phi_m_mat
    % Resulting Val_D is M x M
    Val_D = (1/M) * reshape(sum(exp(j2pi * (phi_m_mat(:) * m_vec)), 2), M, M);

    for kp = 0:N-1
        % Extract leading terms for this specific kp (M x 1)
        L_kp = LeadingTerms(:, kp + 1);
        
        for k = 0:N-1
            % Term C: Scalar for this (kp, k) pair
            phi_n = T * (((kp - k) * Deltaf / N) - nu);
            val_C = (1/N) * sum(exp(-j2pi * n_vec * phi_n));
            
            % Compute the M x M block
            % The leading terms apply to the rows (lp)
            % Val_D covers all (lp, l)
            % val_C is a constant scalar for the block
            Block = (L_kp .* Val_D) * val_C;
            
            % Insert block into the large matrix
            row_idx = kp * M + (1:M);
            col_idx = k * M + (1:M);
            Hp(row_idx, col_idx) = Block;
        end
    end
end
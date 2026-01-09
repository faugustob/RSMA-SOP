function [M, N] = computeOTFSgrid(max_tau, max_nu, method_used, B, delta_f, T, Tf, L_tau, L_nu)
% Computes OTFS grid parameters M, N given max delay and Doppler
%
% Inputs:
%   max_tau: maximum delay spread [s] (direct + RIS)
%   max_nu:  maximum Doppler spread [Hz] (direct + RIS)
%   method_used: 'numerology' (default), 'channel_driven', or 'both'
%   B:       OTFS bandwidth [Hz] (needed for 'numerology')
%   delta_f: subcarrier spacing [Hz]
%   T:       OFDM symbol duration [s] (usually 1/delta_f)
%   Tf:      OTFS frame duration [s] (needed for 'numerology')
%   L_tau:   # delay bins over [0,max_tau] (needed for 'channel_driven')
%   L_nu:    # Doppler bins over [-max_nu,max_nu] (needed for 'channel_driven')
%
% Outputs:
%   M, N: integer grid parameters

% Defaults
if nargin < 3 || isempty(method_used), method_used = 'numerology'; end
if nargin < 8 || isempty(L_tau), L_tau = 4; end
if nargin < 9 || isempty(L_nu), L_nu = 4; end

if strcmpi(method_used, 'numerology')
    % METHOD 1: Standard numerology-driven (most common)
    M_real = B / delta_f;
    N_real = Tf / T;
    M = round(M_real);
    N = round(N_real);
    
    % Verify channel support fits
    tau_span = 1 / delta_f;
    nu_span  = 1 / T;
    
    if max_tau > tau_span
        warning('Delay aliasing: max_tau (%.3e s) > 1/delta_f (%.3e s)', max_tau, tau_span);
    end
    if abs(max_nu) > nu_span
        warning('Doppler aliasing: max_nu (%.3e Hz) > 1/T (%.3e Hz)', max_nu, nu_span);
    end
    
elseif strcmpi(method_used, 'channel_driven')
    % METHOD 2: Drive M,N from max_tau, max_nu via # bins
    Delta_tau = max_tau / L_tau;
    Delta_nu  = max_nu  / L_nu;
    
    M_real = 1 / (delta_f * Delta_tau);   % = L_tau / (delta_f * max_tau)
    N_real = 1 / (T       * Delta_nu);    % = L_nu  / (T       * max_nu)
    
    M = round(M_real);
    N = round(N_real);
    
elseif strcmpi(method_used, 'both')
    % METHOD 3: Numerology first, override if channel doesn't fit
    tau_span = 1 / delta_f;
    nu_span  = 1 / T;
    
    if max_tau <= tau_span && abs(max_nu) <= nu_span
        % Numerology works
        M = round(B / delta_f);
        N = round(Tf / T);
        fprintf('Using numerology: M=%d, N=%d (channel fits)\n', M, N);
    else
        % Fallback to channel-driven
        Delta_tau = max_tau / L_tau;
        Delta_nu  = max_nu  / L_nu;
        M = round(1 / (delta_f * Delta_tau));
        N = round(1 / (T       * Delta_nu));
        fprintf('Using channel-driven: M=%d, N=%d (numerology failed)\n', M, N);
    end
    
else
    error('method_used must be ''numerology'', ''channel_driven'', or ''both''');
end

fprintf('[OTFS Grid] Method: %s, M = %d, N = %d\n', method_used, M, N);
end

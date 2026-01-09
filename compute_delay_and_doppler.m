function [taus, nus, u_paths] = compute_delay_and_doppler( ...
    c, r_tx, v_tx, r_rx, v_rx, f_c, P, sigma_ang)

% ============================================================
% General delay–Doppler model for any TX → RX link
%
% Inputs:
%   c         : speed of light [m/s]
%   r_tx      : transmitter position [3x1]
%   v_tx      : transmitter velocity [3x1]
%   r_rx      : receiver position [3x1]
%   v_rx      : receiver velocity [3x1]
%   f_c       : carrier frequency [Hz]
%   P         : number of propagation paths
%   sigma_ang : angular spread (rad) for NLOS paths
%
% Outputs:
%   taus      : path delays [P x 1]
%   nus       : path Doppler shifts [P x 1]
%   u_paths   : path direction unit vectors [3 x P]
%
% Notes:
%   • Path 1 is LOS
%   • Paths 2..P are clustered NLOS paths
%   • Same delay & Doppler per path (shared by RIS elements if applicable)
% ============================================================

    taus    = zeros(P,1);
    nus     = zeros(P,1);
    u_paths = zeros(3,P);

    % ============================================================
    % LOS path (p = 1)
    % ============================================================
    d_LOS = norm(r_rx - r_tx);
    u_LOS = (r_rx - r_tx) / d_LOS;

    taus(1)      = d_LOS / c;
    nus(1)       = -(f_c / c) * dot(v_tx - v_rx, u_LOS);
    u_paths(:,1) = u_LOS;

    % ============================================================
    % NLOS paths (p = 2,...,P)
    % ============================================================
    for p = 2:P

        % --- Excess delay ---
        Delta_tau = 1e-9 + 1e-5 * rand();
        taus(p) = taus(1) + Delta_tau;

        % --- Direction (clustered around LOS) ---
        delta = sigma_ang * randn(3,1);
        u_p = u_LOS + delta;
        u_p = u_p / norm(u_p);

        u_paths(:,p) = u_p;

        % --- Doppler ---
        nus(p) = -(f_c / c) * dot(v_tx - v_rx, u_p);
    end
end

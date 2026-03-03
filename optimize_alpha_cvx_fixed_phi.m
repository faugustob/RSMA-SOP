function [alpha] = optimize_alpha_cvx_fixed_phi(phi_St, phi_Sr, zeta_k_St, ...
    K, nF, L, Rmin, Pe, P, Q_j, Plos, PLj, HB, HA, g_pq, Nsymb, ...
    reflect, h_rp, h_jq, h_e, delta_f, Active_Gain_dB, max_SCA)
%% ========================= CONSTANTS =========================
Rmin = 1e-8;
Nr = length(phi_St);
zeta_k_Sr = (10^(Active_Gain_dB/10)) - zeta_k_St;
phase_St = exp(1j .* phi_St);
phase_Sr = exp(1j .* phi_Sr);
beta_St = sqrt(zeta_k_St) .* phase_St;
beta_Sr = sqrt(zeta_k_Sr) .* phase_Sr;
BW = delta_f;
N0_dBm = -174;
sigma2 = 10^((N0_dBm + 10*log10(BW) - 30)/10);
Pw_dBm = 46;
Pw = 10^((Pw_dBm - 30)/10);
AN_P_ratio = 1;          % Increase this (e.g. 5-10) if eavesdroppers are too strong
noise = max(sigma2/Pw, 1e-10);

%% ========================= PRECOMPUTE CHANNELS =========================
Nc_k_all = zeros(K,1);
Nc_k_AN_all = zeros(K,1);
Nc_l_all = zeros(nF,K);
Nc_l_AN_all = zeros(nF,K);

for k = 1:K
     reflect_coeff = reflect(k);
        beta_r = (reflect_coeff == 1) * beta_Sr + (reflect_coeff == -1) * beta_St;

    Nc_k_all(k) = compute_OTFS_static_channel( ...
        0, Pe, P, Q_j, Plos(k,1), PLj(k,1), Nr, ...
        HB(:,:,:,k), HA(:,:,:,:,k), g_pq(:,:,k), ...
        beta_r, Nsymb, h_rp(:,:,k,1), h_jq(:,:,k), ...
        h_e(:,k,1), 'vectorized');
    Nc_k_AN_all(k) = compute_OTFS_static_channel( ...
        0, Pe, P, Q_j, Plos(k,2), PLj(k,2), Nr, ...
        HB(:,:,:,k), HA(:,:,:,:,k), g_pq(:,:,k), ...
        beta_r, Nsymb, h_rp(:,:,k,2), h_jq(:,:,k), ...
        h_e(:,k,2), 'vectorized');
end
for l = 1:nF
    for k = 1:K
        reflect_coeff = reflect(k);
        beta_r = (reflect_coeff == 1) * beta_Sr + (reflect_coeff == -1) * beta_St;
        Nc_l_all(l,k) = compute_OTFS_static_channel( ...
            1, Pe, P, Q_j, Plos(K+l,1), PLj(K+l,1), Nr, ...
            HB(:,:,:,K+l), HA(:,:,:,:,K+l), g_pq(:,:,K+l), ...
            beta_r, Nsymb, h_rp(:,:,K+l,1), h_jq(:,:,K+l), ...
            h_e(:,K+l,1), 'vectorized');
        Nc_l_AN_all(l,k) = compute_OTFS_static_channel( ...
            1, Pe, P, Q_j, Plos(K+l,2), PLj(K+l,2), Nr, ...
            HB(:,:,:,K+l), HA(:,:,:,:,K+l), g_pq(:,:,K+l), ...
            beta_r, Nsymb, h_rp(:,:,K+l,2), h_jq(:,:,K+l), ...
            h_e(:,K+l,2), 'vectorized');
    end
end

%% ========================= INITIALIZATION =========================
alpha_init = ones(K+1,1) / (K+1);
alpha_c_init = alpha_init(1);
alpha_pi = alpha_init(2:end);
sum_alpha_pi = sum(alpha_pi);

I_j_prev = zeros(K,1); I_c_prev = zeros(K,1);
gamma_j_prev = zeros(K,1); gamma_c_prev = zeros(K,1);
I_l_prev = zeros(nF,K); gamma_l_prev = zeros(nF,K);

for k = 1:K
    I_j_prev(k) = (sum_alpha_pi - alpha_pi(k)) * Nc_k_all(k) + AN_P_ratio * Nc_k_AN_all(k);
    gamma_j_prev(k) = (alpha_pi(k) * Nc_k_all(k)) / (I_j_prev(k) + noise);

    I_c_prev(k) = sum_alpha_pi * Nc_k_all(k) + AN_P_ratio * Nc_k_AN_all(k);
    gamma_c_prev(k) = (alpha_c_init * Nc_k_all(k)) / (I_c_prev(k) + noise);

    for l = 1:nF
        I_l_prev(l,k) = (sum_alpha_pi - alpha_pi(k)) * Nc_l_all(l,k) + AN_P_ratio * Nc_l_AN_all(l,k);
        gamma_l_prev(l,k) = (alpha_pi(k) * Nc_l_all(l,k)) / (I_l_prev(l,k) + noise);
    end
end
gamma_j_prev = max(gamma_j_prev, 1e-6);
gamma_c_prev = max(gamma_c_prev, 1e-6);
gamma_l_prev = max(gamma_l_prev, 1e-6);

%% ========================= SCA LOOP =========================
for sca_iter = 1:max_SCA
    cvx_begin quiet
        cvx_solver mosek
        
        variable vecAlpha(K+1) nonnegative
        variable gamma_j(K) nonnegative
        variable gamma_l(nF,K) nonnegative
        variable gamma_c(K) nonnegative
        variable Rc nonnegative
        variable s_fake(nF,K)            % can be negative (we take max(0,.) later if needed)

        %maximize( (1/(nF*K)) * sum(sum(s_fake)) )

        maximize( min(min(s_fake)) )

        subject to
            sum(vecAlpha) <= 1;
            vecAlpha >= 1e-2;

            alpha_c = vecAlpha(1);
            alpha_pi = vecAlpha(2:end);
            sum_alpha_pi = sum(alpha_pi);

            % ---------- USER-LEVEL CONSTRAINTS ----------
            for k = 1:K
 
                Rc <= log(1 + gamma_c(k))/log(2);

                S_j = alpha_pi(k) * Nc_k_all(k);
                I_j = (sum_alpha_pi - alpha_pi(k)) * Nc_k_all(k) + AN_P_ratio * Nc_k_AN_all(k);
                S_c = alpha_c * Nc_k_all(k);
                I_c = sum_alpha_pi * Nc_k_all(k) + AN_P_ratio * Nc_k_AN_all(k);

                % Correct SCA linearization (upper bound on bilinear)
                gamma_j_prev(k)*I_j + I_j_prev(k)*gamma_j(k) - gamma_j_prev(k)*I_j_prev(k) + gamma_j(k)*noise <= S_j;
                gamma_c_prev(k)*I_c + I_c_prev(k)*gamma_c(k) - gamma_c_prev(k)*I_c_prev(k) + gamma_c(k)*noise <= S_c;
            end

            % ---------- EAVESDROPPER / SECRECY CONSTRAINTS ----------
            for l = 1:nF
                for k = 1:K
                    S_l = alpha_pi(k) * Nc_l_all(l,k);
                    I_l = (sum_alpha_pi - alpha_pi(k)) * Nc_l_all(l,k) + AN_P_ratio * Nc_l_AN_all(l,k);

                    % Corrected: upper bound on gamma_l (same direction as legitimate users)
                    gamma_l_prev(l,k)*I_l + I_l_prev(l,k)*gamma_l(l,k) - gamma_l_prev(l,k)*I_l_prev(l,k) + gamma_l(l,k)*noise <= S_l;

                    % First-order upper bound on log(1+gamma_l) → lower bound on secrecy
                    log_l_approx = log(1 + gamma_l_prev(l,k))/log(2) + ...
                                   (1/((1 + gamma_l_prev(l,k))*log(2))) * (gamma_l(l,k) - gamma_l_prev(l,k));

                    s_fake(l,k) <= log(1 + gamma_j(k))/log(2) - log_l_approx;
                end
            end

    cvx_end

    if ~strcmp(cvx_status,'Solved') && ~strcmp(cvx_status,'Inaccurate/Solved')
        fprintf('SCA failed at iter %d: %s\n', sca_iter, cvx_status);
        fprintf('Try increasing AN_P_ratio (e.g. to 5-10) or check channel strengths.\n');
        break;
    end

    % ---------- UPDATE FOR NEXT ITERATION ----------
    alpha_val = double(vecAlpha);
    alpha_pi_val = alpha_val(2:end);
    sum_pi_val = sum(alpha_pi_val);

    for k = 1:K
        gamma_j_prev(k) = double(gamma_j(k));
        gamma_c_prev(k) = double(gamma_c(k));
        I_j_prev(k) = (sum_pi_val - alpha_pi_val(k)) * Nc_k_all(k) + AN_P_ratio * Nc_k_AN_all(k);
        I_c_prev(k) = sum_pi_val * Nc_k_all(k) + AN_P_ratio * Nc_k_AN_all(k);

        for l = 1:nF
            gamma_l_prev(l,k) = double(gamma_l(l,k));
            I_l_prev(l,k) = (sum_pi_val - alpha_pi_val(k)) * Nc_l_all(l,k) + AN_P_ratio * Nc_l_AN_all(l,k);
        end
    end

    alpha = alpha_val;
end
cvx_clear;
end
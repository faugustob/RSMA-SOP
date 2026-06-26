function [sc_c_lk, sc_p_lk, rate_c_min, rate_p_vec, R_k, ...
          sinr_c_k, sinr_p_k, sinr_c_l, sinr_p_l] = ...
    compute_sinr_sc_an_ofdm(Pe, P, Q_j, L, K, delta_f, Plos, PLj, Nr, ...
                            HB, HA, g_pq, Nsymb, reflect, Rmin, ...
                            h_rp, h_jq, h_e, zeta_k_St, Active_Gain_dB, ...
                            AN_P_ratio, x, eq_mode)

if nargin < 23 || isempty(eq_mode)
    eq_mode = 'none';
end

alpha = x(1:K+1);
phi_St = wrapToPi(x(K+2:K+1+Nr));
phi_Sr = 0;

if sum(alpha) - 1 > 1e-10
    error('Illegal RSMA power allocation');
end

any_reflect = any(reflect > 0) && any(reflect < 0);
if any_reflect
    phi_Sr = wrapToPi(x(K+2+Nr:K+1+2*Nr));
    zeta_k_St = x(K+2+2*Nr:K+1+3*Nr);
end

zeta_k_Sr = (10^(Active_Gain_dB/10)) - zeta_k_St;
phase_St = exp(1j .* phi_St);
phase_Sr = exp(1j .* phi_Sr);
beta_St = sqrt(zeta_k_St) .* phase_St;
beta_Sr = sqrt(zeta_k_Sr) .* phase_Sr;

BW = delta_f;
N0_dBm = -174;
sigma2 = 10^((N0_dBm + 10*log10(BW) - 30)/10);
Pw_dBm = 46;
Pw = 10^((Pw_dBm-30)/10);

ris_noise_term = 0;

alpha_c = alpha(1);
alpha_pi_v = alpha(2:end);
sum_alpha_pi = sum(alpha_pi_v);

% Preallocate
sinr_c_k = zeros(K,1);
sinr_p_k = zeros(K,1);
sinr_c_l = zeros(L,1);
sinr_p_l = zeros(L,K);

% ====================== Legitimate Users ======================
for k = 1:K
    reflect_coeff = reflect(k);
    beta_r = (reflect_coeff == 1) * beta_Sr + (reflect_coeff == -1) * beta_St;

    [~, ~, Heff_k]     = compute_OTFS_static_channel_ofdm(0, Pe, P, Q_j, Plos(k,1), PLj(k,1), Nr, ...
        HB(:,:,:,k), HA(:,:,:,:,k), g_pq(:,:,k), beta_r, Nsymb, ...
        h_rp(:,:,1), h_jq(:,:,k,1), h_e(:,k,1), 'vectorized');

    [~, ~, Heff_k_AN]  = compute_OTFS_static_channel_ofdm(0, Pe, P, Q_j, Plos(k,2), PLj(k,2), Nr, ...
        HB(:,:,:,k), HA(:,:,:,:,k), g_pq(:,:,k), beta_r, Nsymb, ...
        h_rp(:,:,2), h_jq(:,:,k,2), h_e(:,k,2), 'vectorized');

    sinr_c_k(k) = compute_post_eq_sinr(Heff_k, Heff_k_AN, alpha_c, ...
        sum_alpha_pi, AN_P_ratio, sigma2/Pw, ris_noise_term, Nsymb, eq_mode);

    sinr_p_k(k) = compute_post_eq_sinr(Heff_k, Heff_k_AN, alpha_pi_v(k), ...
        sum_alpha_pi - alpha_pi_v(k), AN_P_ratio, sigma2/Pw, ...
        ris_noise_term, Nsymb, eq_mode);
end

% ====================== Eavesdroppers ======================
for l = 1:L
    reflect_coeff = reflect(K+l);
    beta_r = (reflect_coeff == 1) * beta_Sr + (reflect_coeff == -1) * beta_St;

    [~, ~, Heff_l]    = compute_OTFS_static_channel_ofdm(1, Pe, P, Q_j, Plos(K+l,1), PLj(K+l,1), Nr, ...
        HB(:,:,:,K+l), HA(:,:,:,:,K+l), g_pq(:,:,K+l), beta_r, Nsymb, ...
        h_rp(:,:,1), h_jq(:,:,K+l,1), h_e(:,K+l,1), 'vectorized');

    [~, ~, Heff_l_AN] = compute_OTFS_static_channel_ofdm(1, Pe, P, Q_j, Plos(K+l,2), PLj(K+l,2), Nr, ...
        HB(:,:,:,K+l), HA(:,:,:,:,K+l), g_pq(:,:,K+l), beta_r, Nsymb, ...
        h_rp(:,:,2), h_jq(:,:,K+l,2), h_e(:,K+l,2), 'vectorized');

    sinr_c_l(l) = compute_post_eq_sinr(Heff_l, Heff_l_AN, alpha_c, ...
        sum_alpha_pi, AN_P_ratio, sigma2/Pw, ris_noise_term, Nsymb, eq_mode);

    for kk = 1:K
        sinr_p_l(l,kk) = compute_post_eq_sinr(Heff_l, Heff_l_AN, ...
            alpha_pi_v(kk), sum_alpha_pi - alpha_pi_v(kk), AN_P_ratio, ...
            sigma2/Pw, ris_noise_term, Nsymb, eq_mode);
    end
end

% ====================== Rate Calculations ======================
rate_c_min = min(log2(1 + sinr_c_k));
rate_p_vec = log2(1 + sinr_p_k);

sc_c_lk = zeros(L,1);
for l = 1:L
    sc_c_lk(l) = rate_c_min - log2(1 + sinr_c_l(l));
end

sc_p_lk = zeros(L,K);
for l = 1:L
    for k = 1:K
        sc_p_lk(l,k) = rate_p_vec(k) - log2(1 + sinr_p_l(l,k));
    end
end

% RSMA rate allocation
C_k = zeros(K,1);
rate_c_available = rate_c_min;
for k = 1:K
    deficit = max(Rmin - rate_p_vec(k), 0);
    C_k(k) = min(deficit, rate_c_available);
    rate_c_available = max(rate_c_available - C_k(k), 0);
end
if rate_c_available > 0
    C_k = C_k + (rate_c_available / K);
end
R_k = rate_p_vec(:) + C_k;
end


function sinr = compute_post_eq_sinr(Heff, Heff_AN, alpha_sig, alpha_interf, ...
                                     AN_ratio, noise_var, ris_noise, Nsymb, eq_mode)
    N = Nsymb;
    total_noise = noise_var + ris_noise;

    if strcmpi(eq_mode, 'none')
        % === No Equalization: Evaluate average power per resource element ===
        % Normalizing by N protects against non-physical dimension scaling
        sig_pwr   = alpha_sig * mean(abs(Heff(:)).^2);
        int_pwr   = alpha_interf * mean(abs(Heff(:)).^2) + AN_ratio * mean(abs(Heff_AN(:)).^2);
        sinr      = sig_pwr / (int_pwr + total_noise + 1e-12);
        return;
    end

    % ====================== ZF / MMSE Equalization ======================
    if strcmpi(eq_mode, 'zf')
        W = pinv(Heff);
    elseif strcmpi(eq_mode, 'mmse')
        W = (Heff' * Heff + total_noise * eye(N)) \ Heff';
    else
        error('Unknown eq_mode: %s', eq_mode);
    end

    % Effective equalized channel matrices
    Heff_eq    = W * Heff;
    Heff_AN_eq = W * Heff_AN;

    % Extract desired signal (Main Diagonal components only)
    diag_main    = diag(Heff_eq);
    signal_power = alpha_sig * mean(abs(diag_main).^2);

    % Residual ISI/ICI (Off-diagonal energy of the desired stream)
    total_sig_stream_power = alpha_sig * sum(abs(Heff_eq(:)).^2) / N;
    residual_isi_power     = total_sig_stream_power - signal_power;

    % Co-channel RSMA Interference Power (Incoherent)
    interf_data_power = alpha_interf * sum(abs(Heff_eq(:)).^2) / N;

    % Artificial Noise Power (Incoherent)
    interf_an_power   = AN_ratio * sum(abs(Heff_AN_eq(:)).^2) / N;

    % Noise Enhancement Calculation (Average per symbol)
    noise_enhancement = sum(abs(W(:)).^2) / N; 
    enhanced_noise    = total_noise * noise_enhancement;

    % Total SINR calculation
    total_interference = residual_isi_power + interf_data_power + interf_an_power;
    sinr = signal_power / (total_interference + enhanced_noise + 1e-12);
end
function [sc_c_lk,sc_p_lk,sc_p_kk,rate_c_min,rate_p_vec,R_k,sinr_c_k, sinr_p_k, sinr_c_l, sinr_p_l] = compute_sinr_sc_an(Pe,P,Q_j,L,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,Rmin,h_rp,h_jq,h_e,zeta_k_St,Active_Gain_dB,x)
    %
    % phi=[];
    % alpha=[];
    alpha = x(1:K+1);
    phi_St = x(K+2:K+1+Nr);
    phi_Sr=0;


    if sum(alpha) - 1> 1e-15
        error('Illegal RSMA power allocation');
    end
    any_reflect = any(reflect > 0) && any(reflect < 0);

    if any_reflect
        phi_Sr =  x(K+2+Nr:K+1+2*Nr);
        zeta_k_St = x(K+2+2*Nr:K+1+3*Nr);
    end

    if gpuDeviceCount > 0
         % --- GPU CHANGES: Convert x-derived vars to GPU if not already ---
        alpha = gpuArray(alpha);  % Low-cost; ensures RSMA ops on GPU
        phi_St = gpuArray(phi_St);
        phi_Sr = gpuArray(phi_Sr);

        zeta_k_St = gpuArray(zeta_k_St);
    end
    zeta_k_Sr = (10^(Active_Gain_dB/10))  - zeta_k_St;
   
    phase_Sr = exp(1j .* phi_Sr);  % Complex exp on GPU
    phase_St = exp(1j .* phi_St);  % Complex exp on GPU

    beta_St = sqrt(zeta_k_St) .* phase_St;  % Assuming zeta_k_Sr pre-GPU'd outside; else gpuArray it
    beta_Sr = sqrt(zeta_k_Sr) .* phase_Sr;

    % beta_St = sqrt(zeta_k_St).*phase;  % Uncomment if needed

    % Use the realistic SNR values we discussed
    % snr_e_db = 5;
    % snr_db = 20;
    %
    % b_a = 1 / 10^(snr_db/10);
    % b_e = 1 / 10^(snr_e_db/10);

    BW = delta_f; % Subcarrier bandwidth (Hz)
    N0_dBm = -174; % Thermal noise density (dBm/Hz)
    sigma2 = 10^((N0_dBm + 10*log10(BW) - 30)/10);

    Pw_dBm = 46; % LEO RF transmit power (10 W)
    Pw = 10^((Pw_dBm-30)/10); % Watts
    AN_P_ratio = 1;

    % 1. Define RIS Noise parameters
    NF_ris = 0; % Noise Figure of the RIS amplifiers (in dB)
    sigma2_ris = 10^((N0_dBm + 10*log10(BW) + NF_ris - 30)/10); 
    
    % 2. Calculate the Total Amplified RIS Noise (normalized by Pw)
    % Note: zeta_k_Sr is the per-element gain (e.g., 1000 for 30dB)
    ris_noise_term = (Nr * (10^(Active_Gain_dB/10)) * sigma2_ris) / Pw;
    ris_noise_term = 0;
   
    %Pw=80;

    alpha_c = alpha(1);
    alpha_pi_v = alpha(2:end);

    % --- GPU CHANGES: Precompute sums once on GPU ---
    sum_alpha_pi = sum(alpha_pi_v);  % Scalar, but on GPU

    sinr_c_k = zeros(K,1);  % Allocate output
    sinr_p_k = zeros(K,1);

    % --- Legitimate Users ---
    for k = 1:K
        reflect_coeff = reflect(k);
        beta_r = (reflect_coeff == 1) * beta_Sr + (reflect_coeff == -1) * beta_St;

        % Channel call: Inputs already GPU'd outside, so stays on GPU
        [Nc_k] = compute_OTFS_static_channel(0,Pe,P,Q_j,Plos(k,1),PLj(k,1),Nr,HB(:,:,:,k),HA(:,:,:,:,k),g_pq(:,:,k),beta_r,Nsymb,h_rp(:,:,k,1),h_jq(:,:,k),h_e(:,k,1),'vectorized'); %Data Channel
        [Nc_k_AN] = compute_OTFS_static_channel(0,Pe,P,Q_j,Plos(k,2),PLj(k,2),Nr,HB(:,:,:,k),HA(:,:,:,:,k),g_pq(:,:,k),beta_r,Nsymb,h_rp(:,:,k,2),h_jq(:,:,k),h_e(:,k,2),'vectorized'); %Data Channel

        % RSMA SINR Logic: Common message sees all private power as interference
        signal_c = alpha_c * Nc_k;
        interf_c = sum_alpha_pi * Nc_k+AN_P_ratio*Nc_k_AN;  % Use precomputed sum
        sinr_c_k(k) = signal_c / (interf_c +  sigma2/Pw + ris_noise_term);

        % Private message only sees OTHER private messages as interference (after SIC)
        signal_p = alpha_pi_v(k) * Nc_k;
        interf_p = (sum_alpha_pi - alpha_pi_v(k)) * Nc_k+AN_P_ratio*Nc_k_AN;  % Use precomputed
        sinr_p_k(k) = signal_p / (interf_p + sigma2/Pw + ris_noise_term);
    end

    % --- Eavesdroppers ---
    sinr_c_l = zeros(L,1);
    sinr_p_l = zeros(L,K);

    for l=1:L
        reflect_coeff = reflect(K+l);
        beta_r = (reflect_coeff == 1) * beta_Sr + (reflect_coeff == -1) * beta_St;

        [Nc_l] = compute_OTFS_static_channel(1,Pe,P,Q_j,Plos(K+l,1),PLj(K+l,1),Nr,HB(:,:,:,K+l),HA(:,:,:,:,K+l),g_pq(:,:,K+l),beta_r,Nsymb,h_rp(:,:,K+l,1),h_jq(:,:,K+l),h_e(:,K+l,1),'vectorized');
        [Nc_l_AN] = compute_OTFS_static_channel(1,Pe,P,Q_j,Plos(K+l,2),PLj(K+l,2),Nr,HB(:,:,:,K+l),HA(:,:,:,:,K+l),g_pq(:,:,K+l),beta_r,Nsymb,h_rp(:,:,K+l,2),h_jq(:,:,K+l),h_e(:,K+l,2),'vectorized');

        signal_c_l = alpha_c * Nc_l;
        interf_c_l = sum_alpha_pi * Nc_l + AN_P_ratio*Nc_l_AN;
        sinr_c_l(l) = signal_c_l / (interf_c_l + sigma2/Pw + ris_noise_term);

        %sinr_c_l(l) = (alpha_c * Nc_l) / ((1 - alpha_c) * Nc_l + b_e);
        for k=1:K
            signal_p_kl = alpha_pi_v(k) * Nc_l;
            interf_p_kl = (sum_alpha_pi - alpha_pi_v(k))  * Nc_l + AN_P_ratio*Nc_l_AN;  % Use precomputed
            sinr_p_l(l,k) = (signal_p_kl) / (interf_p_kl + sigma2/Pw + ris_noise_term);
        end
    end

    % --- RSMA Rate Calculation ---
    rate_c_min = min(log2(1 + sinr_c_k)); % The common rate is limited by the worst user
    rate_p_vec = log2(1 + sinr_p_k);      % Private rates for each user

    % --- Secrecy Capacities ---
    sc_c_lk = zeros(L,1); 
    for l = 1:L
        % Common Secrecy: Shared rate minus what the eavesdropper can see
        sc_c_lk(l) = rate_c_min - log2(1 + sinr_c_l(l));
    end
    
    sc_p_lk = zeros(L,K);
    for l = 1:L
        for k = 1:K
            %sc_p_lk(l,k) = max(rate_p_vec(k) - log2(1 + sinr_p_l(l,k)), 0);
            sc_p_lk(l,k) = rate_p_vec(k) - log2(1 + sinr_p_l(l,k));
        end
    end

     C_k = zeros(K,1);
        rate_c_available = rate_c_min;
        for k = 1:K
            deficit = max(Rmin - rate_p_vec(k), 0);
            C_k(k) = min(deficit, rate_c_available);
            rate_c_available = max(rate_c_available - C_k(k), 0);
        end
        
        R_k = rate_p_vec(:) + C_k;

        for k = 1:K
            for kp = 1:K
                %sc_p_lk(l,k) = max(rate_p_vec(k) - log2(1 + sinr_p_l(l,k)), 0);
                sc_p_kk(l,k) = rate_p_vec(k) - rate_p_vec(kp);
            end
        end
end
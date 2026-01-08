%function [sc_c_lk,sc_p_lk,rate_c_min,rate_p_vec,sinr_c_k, sinr_p_k, sinr_c_l, sinr_p_l] = compute_sinr_sc(Pe,P,Q_j,L,K,m_e,m_q,m_p,omega_e,omega_p,omega_q,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,h_rp,h_jq,h_e,phi,zeta_k_Sr,x)
function [sc_c_lk,sc_p_lk,rate_c_min,rate_p_vec,sinr_c_k, sinr_p_k, sinr_c_l, sinr_p_l] = compute_sinr_sc(Pe,P,Q_j,L,K,m_e,m_q,m_p,omega_e,omega_p,omega_q,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,h_rp,h_jq,h_e,alpha,phi,zeta_k_Sr,x)
    % 
    % phi=[];
    % alpha=[];
    alpha = x(1:K+1);
    phi = x(K+2:K+1+Nr);
     % phi = x;
    
    
    phase = exp(1j.*phi);
    beta_r = sqrt(zeta_k_Sr).*phase;
    % beta_St = sqrt(zeta_k_St).*phase;
    
    % Use the realistic SNR values we discussed
    % snr_e_db = 5;
    % snr_db = 20; 
    % 
    % b_a = 1 / 10^(snr_db/10);
    % b_e = 1 / 10^(snr_e_db/10);

    BW = delta_f;                 % Subcarrier bandwidth (Hz)
    N0_dBm = -174;             % Thermal noise density (dBm/Hz)
    sigma2 = 10^((N0_dBm + 10*log10(BW) - 30)/10);
    
    Pw_dBm = 46;               % LEO RF transmit power (10 W)
    Pw = 10^((Pw_dBm-30)/10);  % Watts
    %Pw=80;

    
    alpha_c = alpha(1);
    alpha_pi_v = alpha(2:end); 

    sinr_c_k = zeros(K,1);
    sinr_p_k = zeros(K,1);

    % --- Legitimate Users ---
    for k = 1:K  
        % reflect_coeff = reflect(k);
        % beta_r = (reflect_coeff == 1) * beta_Sr + (reflect_coeff == 0) * beta_St;
        
        [Nc_k] = compute_OTFS_static_channel(0,Pe,P,Q_j,Plos(k),PLj(k),Nr,HB(:,:,:,k),HA(:,:,:,:,k),g_pq(:,:,k),beta_r,Nsymb,h_rp(:,:,k),h_jq(:,:,k),h_e(:,k),'vectorized');
        
        % RSMA SINR Logic: Common message sees all private power as interference
        signal_c = alpha_c * Nc_k;
        interf_c = sum(alpha_pi_v) * Nc_k;
        
        sinr_c_k(k) = signal_c / (interf_c + sigma2/Pw); 
        
        % Private message only sees OTHER private messages as interference (after SIC)
        signal_p = alpha_pi_v(k) * Nc_k;
        interf_p = (sum(alpha_pi_v)-alpha_pi_v(k)) * Nc_k;
        
        sinr_p_k(k) = signal_p / (interf_p + sigma2/Pw);
    end
    
    % --- Eavesdroppers ---
    sinr_c_l = zeros(L,1);
    sinr_p_l = zeros(L,K);
    for l=1:L
        % reflect_coeff = reflect(K+l);
        % beta_r = (reflect_coeff == 1) * beta_Sr + (reflect_coeff == 0) * beta_St;
        
        [Nc_l] = compute_OTFS_static_channel(1,Pe,P,Q_j,Plos(K+l),PLj(K+l),Nr,HB(:,:,:,K+l),HA(:,:,:,:,K+l),g_pq(:,:,K+l),beta_r,Nsymb,h_rp(:,:,K+l),h_jq(:,:,K+l),h_e(:,K+l),'vectorized');
        

        signal_c_l = alpha_c * Nc_l;
        interf_c_l = sum(alpha_pi_v) * Nc_l;
        
        sinr_c_l(l) = signal_c_l / (interf_c_l + sigma2/Pw); 

        %sinr_c_l(l) = (alpha_c * Nc_l) / ((1 - alpha_c) * Nc_l + b_e);
        for k=1:K
            signal_p_kl = alpha_pi_v(k) * Nc_l;
            interf_p_kl = (sum(alpha_pi_v)-alpha_pi_v(k))  * Nc_l;
            
          
            sinr_p_l(l,k) = (signal_p_kl) / (interf_p_kl + sigma2/Pw);
        end
    end

    % --- RSMA Rate Calculation ---
    rate_c_min = min(log2(1 + sinr_c_k)); % The common rate is limited by the worst user
    rate_p_vec = log2(1 + sinr_p_k);      % Private rates for each user

    % --- Secrecy Capacities ---
    sc_c_lk = zeros(L,1); 
    for l = 1:L
        % Common Secrecy: Shared rate minus what the eavesdropper can see
        sc_c_lk(l) = max(rate_c_min - log2(1 + sinr_c_l(l)), 0);
    end
    
    sc_p_lk = zeros(L,K);
    for l = 1:L
        for k = 1:K
            sc_p_lk(l,k) = max(rate_p_vec(k) - log2(1 + sinr_p_l(l,k)), 0);
        end
    end
end
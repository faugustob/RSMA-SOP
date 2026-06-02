function [sc_lk_noma,R_k_noma, sinr_k_noma, sinr_l_noma] = compute_sinr_sc_an_noma(Pe,P,Q_j,L,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,h_rp,h_jq,h_e,zeta_k_St,Active_Gain_dB,AN_P_ratio,x)
 
    alpha = x(1:K);
    phi_St = wrapToPi(x(K+1:K+Nr));
    phi_Sr=0;


    if sum(alpha) - 1> 1e-10
        error('Illegal RSMA power allocation');
    end
    any_reflect = any(reflect > 0) && any(reflect < 0);

    if any_reflect
        phi_Sr =  wrapToPi(x(K+1+Nr:K+2*Nr));
        zeta_k_St = x(K+1+2*Nr:K+3*Nr);
    end


    zeta_k_Sr = (10^(Active_Gain_dB/10))  - zeta_k_St;
   
    phase_Sr = exp(1j .* phi_Sr);                  % Complex exp on GPU
    phase_St = exp(1j .* phi_St);  
    beta_St = sqrt(zeta_k_St) .* phase_St;         % Assuming zeta_k_Sr pre-GPU'd outside; else gpuArray it
    beta_Sr = sqrt(zeta_k_Sr) .* phase_Sr;



    BW = delta_f;                                   % Subcarrier bandwidth (Hz)
    N0_dBm = -174;                                  % Thermal noise density (dBm/Hz)
    sigma2 = 10^((N0_dBm + 10*log10(BW) - 30)/10);

    Pw_dBm = 46;                                    % LEO RF transmit power (10 W)
    Pw = 10^((Pw_dBm-30)/10);                       % Watts
 

    % 1. Define RIS Noise parameters
    NF_ris = 0;                                     % Noise Figure of the RIS amplifiers (in dB)
    sigma2_ris = 10^((N0_dBm + 10*log10(BW) + NF_ris - 30)/10); 
    
    % 2. Calculate the Total Amplified RIS Noise (normalized by Pw)
    ris_noise_term = (Nr * (10^(Active_Gain_dB/10)) * sigma2_ris) / Pw;
    ris_noise_term = 0; 

     

    Nc_k = zeros(K,1);
    Nc_k_AN = zeros(K,1);
    interf_coeff = zeros(K,1);


    for k = 1:K
        reflect_coeff = reflect(k);
        beta_r = (reflect_coeff == 1) * beta_Sr + (reflect_coeff == -1) * beta_St;

        % Channel call: Inputs already GPU'd outside, so stays on GPU
        Nc_k(k) = compute_OTFS_static_channel(0,Pe,P,Q_j,Plos(k,1),PLj(k,1),Nr,HB(:,:,:,k),HA(:,:,:,:,k),g_pq(:,:,k),beta_r,Nsymb,h_rp(:,:,1),h_jq(:,:,k,1),h_e(:,k,1),'vectorized'); %Data Channel
        Nc_k_AN(k) = compute_OTFS_static_channel(0,Pe,P,Q_j,Plos(k,2),PLj(k,2),Nr,HB(:,:,:,k),HA(:,:,:,:,k),g_pq(:,:,k),beta_r,Nsymb,h_rp(:,:,2),h_jq(:,:,k,2),h_e(:,k,2),'vectorized'); %Data Channel

    end


    % Sort users: weakest to strongest channel (important for SIC)
    [~, order] = sort(Nc_k);   
    alpha_sorted = alpha(order);     
    sinr_k_noma = zeros(K,1);

    % --- Legitimate Users ---
    for j = 1:K  

        k = order(j);
       
       if j < K
            interf_coeff(k) = sum(alpha_sorted(j+1:end));   
        else
            interf_coeff(k) = 0; % Strongest user cancels all intra-cluster interference
        end

       % Private message only sees OTHER private messages as interference (after SIC)
        signal= alpha(k) * Nc_k(k);
        interf = interf_coeff(k) * Nc_k(k)+AN_P_ratio*Nc_k_AN(k);  % Use precomputed
        sinr_k_noma(k) = signal / (interf + sigma2/Pw + ris_noise_term);
    end

    % --- Eavesdroppers ---
    sinr_l_noma = zeros(L,K);

    for l=1:L
        reflect_coeff = reflect(K+l);
        beta_r = (reflect_coeff == 1) * beta_Sr + (reflect_coeff == -1) * beta_St;

        [Nc_l] = compute_OTFS_static_channel(1,Pe,P,Q_j,Plos(K+l,1),PLj(K+l,1),Nr,HB(:,:,:,K+l),HA(:,:,:,:,K+l),g_pq(:,:,K+l),beta_r,Nsymb,h_rp(:,:,1),h_jq(:,:,K+l,1),h_e(:,K+l,1),'vectorized');
        [Nc_l_AN] = compute_OTFS_static_channel(1,Pe,P,Q_j,Plos(K+l,2),PLj(K+l,2),Nr,HB(:,:,:,K+l),HA(:,:,:,:,K+l),g_pq(:,:,K+l),beta_r,Nsymb,h_rp(:,:,2),h_jq(:,:,K+l,2),h_e(:,K+l,2),'vectorized');

       
        for k=1:K
            signal_p_kl = alpha(k) * Nc_l;
            interf_p_kl = interf_coeff(k)  * Nc_l + AN_P_ratio*Nc_l_AN;  % Use precomputed
            sinr_l_noma(l,k) = (signal_p_kl) / (interf_p_kl + sigma2/Pw + ris_noise_term);
        end
    end

   
    rate_vec = log2(1 + sinr_k_noma);      % Private rates for each user

    % --- Secrecy Capacities ---
  
    
    sc_lk_noma = zeros(L,K);
    for l = 1:L
        for k = 1:K
            %sc_p_lk(l,k) = max(rate_p_vec(k) - log2(1 + sinr_p_l(l,k)), 0);
            sc_lk_noma(l,k) = rate_vec(k) - log2(1 + sinr_l_noma(l,k));
        end
    end

   
        
        R_k_noma = rate_vec(:);
     
end
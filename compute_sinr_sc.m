function [sc_c_lk,sc_p_lk,rate_c,rate_k,sinr_c_k, sinr_p_k, sinr_c_l, sinr_p_l] = compute_sinr_sc(Pe,P,Q_j,L,K,m_e,m_q,m_p,omega_e,omega_p,omega_q,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,h_rp,h_jq,h_e,x)
    


   

    alpha = x(1:K+1);
    phi = x(K+2:K+1+Nr);
    zeta_k_Sr = x(K+2+Nr:K+1+2*Nr);

    zeta_k_St = 1 - zeta_k_Sr; % transmit coefficients

    

    phase = exp(1j.*phi.*2.*pi);


    beta_Sr = sqrt(zeta_k_Sr).*phase;
    beta_St = sqrt(zeta_k_St).*phase;

    

    snr_e_db = 50;
    
    snr_e = 10^(snr_e_db/10);
    
    
    snr_db =  90;
    
    snr_linear = 10^(snr_db/10);
    
    b_a = 1/snr_linear;
    b_e = 1/snr_e;

    alpha_c = alpha(1); % power allocation factor for common message
    alpha_pi_v = alpha(2:end); % power allocation factor for private message of each user
    % power allocation coefficients must satisfy constraint:
    % alpha_c+alpha_pi_v*K <= 1
    sinr_c_k = zeros(K,1);
    sinr_p_k = zeros(K,1);
    for k =1:K  

       
        % h_rp = zeros(Nr, P);
        % for p = 1:P
        %     % Scale parameter theta = omega/m
        %     h_rp(:, p) = sqrt(gamrnd(m_p(p), omega_p(p)/m_p(p), [Nr, 1])) .* exp(1i*2*pi*rand(Nr, 1));
        % end 
        % 
        % 
        % 
        % h_jq = zeros(Nr, Q_j);
        % for q = 1:Q_j
        %     h_jq(:, q) = sqrt(gamrnd(m_q(k,q), omega_q(q)/m_q(k,q), [Nr, 1])) .* exp(1i*2*pi*rand(Nr,1));
        % end
        % 
        % 
        % h_e = zeros(Pe,1);
        % for u = 1:Pe
        %     h_e(u) = sqrt(gamrnd(m_e(k,u), omega_e(u)/m_e(k,u))) .* exp(1i*2*pi*rand());
        % end
      

        reflect_coeff = reflect(k);

        if reflect_coeff == 1
            beta_r = beta_Sr;
        else
            beta_r = beta_St;
        end
        [Nc_k] = compute_OTFS_static_channel(1,Pe,P,Q_j,Plos(k),PLj(k),Nr,HB(:,:,:,k),HA(:,:,:,:,k),g_pq(:,:,k),beta_r,Nsymb,h_rp(:,:,k),h_jq(:,:,k),h_e(:,k));
        % sinr_c_k(k) = (alpha_c*Nc_k)/(sum(sqrt(alpha_pi_v))^2*Nc_k+b_a);
        % sinr_p_k(k) = (alpha_pi_v(k)*Nc_k)/((sum(sqrt(alpha_pi_v))-sqrt(alpha_pi_v(k)))^2*Nc_k+b_a);
        sinr_c_k(k) = (alpha_c * Nc_k) / (sum(alpha_pi_v) * Nc_k + b_a);  % Linear sum alpha_pi
        interference_p_k = sum(alpha_pi_v) - alpha_pi_v(k);  % sum_{v≠k} alpha_pi(v)
        sinr_p_k(k) = (alpha_pi_v(k) * Nc_k) / (interference_p_k * Nc_k + b_a);
    
    end
    
    sinr_c_l = zeros(L,1);
    sinr_p_l = zeros(L,K);
    for l=1:L
        
        % h_rp = zeros(Nr, P);
        % for p = 1:P
        %     % Scale parameter theta = omega/m
        %     h_rp(:, p) = sqrt(gamrnd(m_p(p), omega_p(p)/m_p(p), [Nr, 1])) .* exp(1i*2*pi*rand(Nr, 1));
        % end
        % 
        % h_jq = zeros(Nr, Q_j);
        % for q = 1:Q_j
        %     h_jq(:, q) = sqrt(gamrnd(m_q(K+l,q), omega_q(q)/m_q(K+l,q), [Nr, 1])) .* exp(1i*2*pi*rand(Nr,1));
        % end
        % 
        % h_e = zeros(Pe,1);
        % for u = 1:Pe
        %     h_e(u) = sqrt(gamrnd(m_e(K+l,u), omega_e(u)/m_e(K+l,u))) .* exp(1i*2*pi*rand());
        % end

        reflect_coeff = reflect(K+l);

        if reflect_coeff == 1
            beta_r = beta_Sr;
        else
            beta_r = beta_St;
        end

        [Nc_l] = compute_OTFS_static_channel(1,Pe,P,Q_j,Plos(K+l),PLj(K+l),Nr,HB(:,:,:,K+l),HA(:,:,:,:,K+l),g_pq(:,:,K+l),beta_r,Nsymb,h_rp(:,:,K+l),h_jq(:,:,K+l),h_e(:,K+l));
       sinr_c_l(l) = (alpha_c * Nc_l) / (sum(alpha_pi_v) * Nc_l + b_e);
        for k=1:K
            interference_p_l = sum(alpha_pi_v) - alpha_pi_v(k);
            sinr_p_l(l,k) = (alpha_pi_v(k) * Nc_l) / (interference_p_l * Nc_l + b_e);
        end
    end
    sc_c_lk= zeros(L,K);
    sc_p_lk= zeros(L,K);

    for l=1:L
        for k=1:K
                sc_c_lk(l,k) = log2(1+sinr_c_k(k))-log2(1+sinr_c_l(l));
                sc_p_lk(l,k) = log2(1+sinr_p_k(k))-log2(1+sinr_p_l(l,k));
        end
    end

    rate_c = log2(1+sinr_c_k);
    rate_k = log2(1+sinr_p_k(k));
end
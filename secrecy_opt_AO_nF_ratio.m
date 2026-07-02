clear; clc;
cvx_clear;

% --- Choose how many workers (cores) you want ---
numWorkers = 2;          
pool = gcp('nocreate');
if ~isempty(pool)
    delete(pool);   
end
parpool('local', numWorkers);  

Ns = 20000; % number of samples for Monte Carlo simulation
transmissionType = 'mc';
c = physconst('LightSpeed'); 
f_c = 10e9; 
lambda = c / f_c; 

% ---- Earth-centered conversion ----
R_earth = 6371e3;
nSat = 2;
M = 16;
N = 16;
Nsymb = M*N; 

K_h = 1;  % high speed legit users
K_s = 1;  % slow speed legit users
K = K_h+K_s; 

% --- OTFS System Parameters ---
delta_f = 15e3;      
T = 1/delta_f;       
R = 10;
m_rician = (R+1)^2/(2*R+1);
N_V = 30; 
N_H = 30; 
Nr = N_V * N_H; 
d_x = floor(lambda/2 * 1000) / 1000; 
d_y = floor(lambda/2 * 1000) / 1000; 

%% Radiation patterns
F      = @(theta,phi) cos(theta).^(1/64);
F_tx   = @(theta,phi) cos(theta).^1100;
F_rx   = @(theta,phi) cos(theta).^15; 

%% Gains
G   = 4*pi / integral2(@(t,p) F(t,p).*sin(t),    0, pi/2, 0, 2*pi);
G_t = 4*pi / integral2(@(t,p) F_tx(t,p).*sin(t), 0, pi/2, 0, 2*pi);
G_r = 4*pi / integral2(@(t,p) F_rx(t,p).*sin(t), 0, pi/2, 0, 2*pi);

P = 2; 
Q_j = 3; 
Pe = 2; 

% Coordinates
HAP_altitude = 10e3;
LEO_altitude = 250e3;
RIS_normal = [0;0;1];
[max_theta,max_alpha,max_beta] = compute_max_theta(LEO_altitude,HAP_altitude);

S_elevation = pi/2-max_theta*(0.009);
S_azimuth  = pi/4;
S_radius  = R_earth+LEO_altitude;
S_sph = [ S_azimuth, S_elevation, S_radius]; 

AN_elevation = pi/2-max_theta*(0.1);
AN_azimuth  = pi/4 -pi;
AN_radius  = S_radius;
AN_sph = [ AN_azimuth, AN_elevation, AN_radius]; 

[S_x, S_y, S_z] = sph2cart(S_sph(1), S_sph(2), S_sph(3));
S_xyz = [S_x;S_y;S_z]; 
[AN_x, AN_y, AN_z] = sph2cart(AN_sph(1), AN_sph(2), AN_sph(3));
AN_xyz = [AN_x;AN_y;AN_z];
S_v = 7800;
R_xyz = [0; 0; R_earth+HAP_altitude]; 

RIS_size_x = N_H * d_x;
RIS_size_y = N_V * d_y;
minAngle=(lambda/max(RIS_size_x,RIS_size_y));

m_p = [m_rician; 1*ones(P-1,1)]; 
omega_p = (1/P)*ones(1,P); 
nF_ratio_vec = 0.5:0.5:8;

%% ===================== PRE-ALLOCATION =====================
num_NF = length(nF_ratio_vec);   
feasible_record              = false(Ns, num_NF);
Convex_min_Rk                = zeros(Ns, num_NF);
Convex_Convergence_curve_AO  = zeros(Ns, num_NF);
Convex_Fake_Convergence_curve_AO = zeros(Ns, num_NF);
Convex_Real_Convergence_curve_AO = zeros(Ns, num_NF);
L = 2; 

orbit_normal = [1; 0; 0];
Rs = norm(S_xyz);
omega_orb = (S_v / Rs) * orbit_normal;
vS = cross(omega_orb, S_xyz);
vAN = cross(omega_orb, AN_xyz);
vR = [0;0;0];
sigma_ang = deg2rad(30);   

%% ===================== PARFOR ENGINE =====================
for mc_iter = 1:Ns
    % 1. Determine maximum dimensions for slicing
    max_nF = max(nF_ratio_vec) * L;
    total_nodes_max = K + max_nF + L;
    
    % 2. Pre-allocate maximum sized matrices for this worker iteration
    g_pq_max   = zeros(P, Q_j, total_nodes_max, nSat);
    Plos_max   = zeros(total_nodes_max, nSat);
    PLj_max    = zeros(total_nodes_max, nSat);
    h_jq_max   = zeros(Nr, Q_j, total_nodes_max, nSat);
    h_e_max    = zeros(Pe, total_nodes_max, nSat);
    HA_max     = zeros(Nsymb, Nsymb, P, Q_j, total_nodes_max, nSat);
    HB_max     = zeros(Nsymb, Nsymb, Pe, total_nodes_max, nSat);
    
    taus_kq    = zeros(K, P, Q_j, nSat);
    nus_kq     = zeros(K, P, Q_j, nSat);
    taus_ku    = zeros(Pe, K, nSat);
    nus_ku     = zeros(Pe, K, nSat);
    h_rp       = zeros(Nr, P, nSat);
    
    % Generate Locations for Maximum Scenario Upfront
    el = pi/2 - (0.05)*max_alpha * rand(1, K);
    az = (2*pi) * rand(1, K);
    r = R_earth * ones(1, K);
    [x, y, z] = sph2cart(az, el, r);
    ground_users_cart = [x; y; z];   
    
    x_f = 1000*rand(1, max_nF)+1000;
    y_f = 1000*rand(1, max_nF)+1000;
    z_f = R_earth + 50 + 950*rand(1, max_nF);
    fake_eavesdroppers_xyz = [x_f; y_f; z_f];
    
    x_e = 1000*rand(1, L)+1000;
    y_e = 1000*rand(1, L)+1000;
    z_e = R_earth + 50 + 950*rand(1, L);
    eavesdroppers_xyz = [x_e; y_e; z_e];
    
    % Master coordinate topology matrix
    rho_j_xyz_max = [ground_users_cart, fake_eavesdroppers_xyz, eavesdroppers_xyz];
    reflect_max = sign(RIS_normal.' * (rho_j_xyz_max - R_xyz));
    
    delay_res = 1/(M*delta_f);
    tau_rms = 0.25*delay_res;
    
    [taus_R, nus_R, ~] = compute_delay_and_doppler(c, S_xyz, vS, R_xyz, vR, f_c, P, tau_rms, sigma_ang);
    [taus_R_AN, nus_R_AN, ~] = compute_delay_and_doppler(c, AN_xyz, vAN, R_xyz, vR, f_c, P, tau_rms, sigma_ang);
    
    for p = 1:P
        h_rp(:, p, 1) = sqrt(gamrnd(m_p(p), omega_p(p)/m_p(p), [Nr, 1])) .* exp(1i*2*pi*rand(Nr, 1)); 
        h_rp(:, p, 2) = sqrt(gamrnd(m_p(p), omega_p(p)/m_p(p), [Nr, 1])) .* exp(1i*2*pi*rand(Nr, 1)); 
    end
    
    % Construct Master Statistical Parameters
    m_q_max = ones(total_nodes_max, Q_j);
    omega_q_max = (1/Q_j) * ones(total_nodes_max, Q_j);
    m_e_max = ones(total_nodes_max, Pe);
    omega_e_max = (1/Pe) * ones(total_nodes_max, Pe);
    
    %% --- Compute Channels for Legitimate Users ---
    for k = 1:K
        if k <= K_h, vk_ms = 50 / 3.6; else, vk_ms = 1.2; end
        User_k_loc = rho_j_xyz_max(:, k);
        
        d_los = sqrt(sum((User_k_loc-S_xyz).^2));
        d_los_AN = sqrt(sum((User_k_loc-AN_xyz).^2));
        
        Plos_max(k, 1) = ((lambda/(4*pi))^2*G_t*G_r)/(d_los^2);
        Plos_max(k, 2) = ((lambda/(4*pi))^2*G_t*G_r)/(d_los_AN^2);
        
        PLj_max(k, 1) = compute_ris_PL(lambda,N_V,N_H,S_xyz,User_k_loc,R_xyz,RIS_normal,F,F_tx,F_rx,G,G_t,G_r);
        PLj_max(k, 2) = compute_ris_PL(lambda,N_V,N_H,AN_xyz,User_k_loc,R_xyz,RIS_normal,F,F_tx,F_rx,G,G_t,G_r);
       
        if norm(User_k_loc(1:2)-R_xyz(1:2)) < 1e-6
            d_ru = [R_xyz(1)+1; R_xyz(1)+1];
        else
            d_ru = User_k_loc(1:2) - R_xyz(1:2);
        end
        d_ru = d_ru / norm(d_ru);
        v_l  = vk_ms * [d_ru; 0];
        
        [taus_k, nus_k, ~] = compute_delay_and_doppler(c, R_xyz, vR, User_k_loc, v_l, f_c, Q_j, tau_rms, sigma_ang);
      
        for p=1:P
            for q=1:Q_j
               taus_kq(k,p,q,1) = taus_R(p)+taus_k(q);
               nus_kq(k,p,q,1) = nus_R(p)+nus_k(q);
               taus_kq(k,p,q,2) = taus_R_AN(p)+taus_k(q);
               nus_kq(k,p,q,2) = nus_R_AN(p)+nus_k(q);
            end 
        end
        
        [taus_u, nus_u, ~] = compute_delay_and_doppler(c, S_xyz, vS, User_k_loc, v_l, f_c, Pe, tau_rms, sigma_ang);
        [taus_u_AN, nus_u_AN, ~] = compute_delay_and_doppler(c, AN_xyz, vAN, User_k_loc, v_l, f_c, Pe, tau_rms, sigma_ang);
        
        g_pq_max(:,:,k,1) = exp(1i*2*pi*(taus_R*nus_k'));    
        g_pq_max(:,:,k,2) = exp(1i*2*pi*(taus_R_AN*nus_k'));    
        taus_ku(:,k,1) = taus_u; 
        nus_ku(:,k,1) = nus_u;   
        taus_ku(:,k,2) = taus_u_AN; 
        nus_ku(:,k,2) = nus_u_AN;  
       
        for q = 1:Q_j
            h_jq_max(:, q, k, 1) = sqrt(gamrnd(m_q_max(k,q), omega_q_max(k,q)/m_q_max(k,q), [Nr, 1])) .* exp(1i*2*pi*rand(Nr,1));
            h_jq_max(:, q, k, 2) = sqrt(gamrnd(m_q_max(k,q), omega_q_max(k,q)/m_q_max(k,q), [Nr, 1])) .* exp(1i*2*pi*rand(Nr,1));
        end
        
        for u = 1:Pe
            h_e_max(u, k, 1) = sqrt(gamrnd(m_e_max(k,u), omega_e_max(k,u)/m_e_max(k,u))) .* exp(1i*2*pi*rand());
            h_e_max(u, k, 2) = sqrt(gamrnd(m_e_max(k,u), omega_e_max(k,u)/m_e_max(k,u))) .* exp(1i*2*pi*rand());
        end
        
        for p=1:P
            for q=1:Q_j
                HA_max(:,:,p,q,k,1) = compute_Hp(taus_kq(k,p,q,1), nus_kq(k,p,q,1), M, N, T, delta_f, 'blocked', transmissionType);
                HA_max(:,:,p,q,k,2) = compute_Hp(taus_kq(k,p,q,2), nus_kq(k,p,q,2), M, N, T, delta_f, 'blocked', transmissionType);
            end
        end
       
        for u = 1:Pe
            HB_max(:,:,u,k,1) = compute_Hp(taus_u(u), nus_u(u), M, N, T, delta_f, 'blocked', transmissionType);
            HB_max(:,:,u,k,2) = compute_Hp(taus_u_AN(u), nus_u_AN(u), M, N, T, delta_f, 'blocked', transmissionType);
        end
    end

    %% --- Compute Channels for All Eavesdroppers (Fake + Real) ---
    for l = 1:(max_nF + L)
        node_idx = K + l;
        vl_ms = 0;        
        User_l_loc = rho_j_xyz_max(:, node_idx);
        
        d_los = sqrt(sum((User_l_loc-S_xyz).^2));
        d_los_AN = sqrt(sum((User_l_loc-AN_xyz).^2));
        
        Plos_max(node_idx, 1) = ((lambda/(4*pi))^2*G_t*G_r)/(d_los^2);
        Plos_max(node_idx, 2) = ((lambda/(4*pi))^2*G_t*G_r)/(d_los_AN^2);
        
        PLj_max(node_idx, 1) = compute_ris_PL(lambda,N_V,N_H,S_xyz,User_l_loc,R_xyz,RIS_normal,F,F_tx,F_rx,G,G_t,G_r);
        PLj_max(node_idx, 2) = compute_ris_PL(lambda,N_V,N_H,AN_xyz,User_l_loc,R_xyz,RIS_normal,F,F_tx,F_rx,G,G_t,G_r);
    
        u_re = User_l_loc - R_xyz;
        u_re = u_re / norm(u_re);
        
        a_hat = [0;0;1];
        if abs(dot(a_hat,u_re)) > 0.99, a_hat = [1;0;0]; end
        
        v_dir = cross(a_hat, u_re);
        v_dir = v_dir / norm(v_dir);
        v_l = vl_ms * v_dir;
    
        [taus_l, nus_l, ~] = compute_delay_and_doppler(c, R_xyz, vR, User_l_loc, v_l, f_c, Q_j, tau_rms, sigma_ang);
        [taus_u_l, nus_u_l, ~] = compute_delay_and_doppler(c, S_xyz, vS, User_l_loc, v_l, f_c, Pe, tau_rms, sigma_ang);    
        [taus_u_l_AN, nus_u_l_AN, ~] = compute_delay_and_doppler(c, AN_xyz, vAN, User_l_loc, v_l, f_c, Pe, tau_rms, sigma_ang);
    
        for q = 1:Q_j
            h_jq_max(:, q, node_idx, 1) = sqrt(gamrnd(m_q_max(node_idx,q), omega_q_max(node_idx,q)/m_q_max(node_idx,q), [Nr, 1])) .* exp(1i*2*pi*rand(Nr,1));
            h_jq_max(:, q, node_idx, 2) = sqrt(gamrnd(m_q_max(node_idx,q), omega_q_max(node_idx,q)/m_q_max(node_idx,q), [Nr, 1])) .* exp(1i*2*pi*rand(Nr,1));
        end
        
        for u = 1:Pe
            h_e_max(u, node_idx, 1) = sqrt(gamrnd(m_e_max(node_idx,u), omega_e_max(node_idx,u)/m_e_max(node_idx,u))) .* exp(1i*2*pi*rand());
            h_e_max(u, node_idx, 2) = sqrt(gamrnd(m_e_max(node_idx,u), omega_e_max(node_idx,u)/m_e_max(node_idx,u))) .* exp(1i*2*pi*rand());
        end
    
        for p=1:P
            for q=1:Q_j
                HA_max(:,:,p,q,node_idx,1) = compute_Hp(taus_R(p)+taus_l(q), nus_R(p)+nus_l(q), M, N, T, delta_f, 'blocked', transmissionType);
                HA_max(:,:,p,q,node_idx,2) = compute_Hp(taus_R_AN(p)+taus_l(q), nus_R_AN(p)+nus_l(q), M, N, T, delta_f, 'blocked', transmissionType);
            end
        end
        for u = 1:Pe
            HB_max(:,:,u,node_idx,1) = compute_Hp(taus_u_l(u), nus_u_l(u), M, N, T, delta_f, 'blocked', transmissionType);
            HB_max(:,:,u,node_idx,2) = compute_Hp(taus_u_l_AN(u), nus_u_l_AN(u), M, N, T, delta_f, 'blocked', transmissionType);
        end
        g_pq_max(:,:,node_idx,1) = exp(1i*2*pi*(taus_R*nus_l'));  
        g_pq_max(:,:,node_idx,2) = exp(1i*2*pi*(taus_R_AN*nus_l'));            
    end

    %% ===================== INNER RATIO LOOP (ZERO ALLOCATIONS) =====================
    parfor nFr_idx = 1:num_NF
        nF = nF_ratio_vec(nFr_idx) * L;    
        
        % Formulate slice matrix selectors
        % Keeps: [Legit Users, Slice of Active Fake Eves, True Eves at the very end]
        slice_idx = [1:K, (K+1):(K+nF), (total_nodes_max-L+1):total_nodes_max];
        
        % Slice arrays down without generating heap allocation overhead
        HA   = HA_max(:, :, :, :, slice_idx, :);
        HB   = HB_max(:, :, :, slice_idx, :);
        g_pq = g_pq_max(:, :, slice_idx, :);
        Plos = Plos_max(slice_idx, :);
        PLj  = PLj_max(slice_idx, :);
        h_jq = h_jq_max(:, :, slice_idx, :);
        h_e  = h_e_max(:, slice_idx, :);
        reflect = reflect_max(slice_idx);
        
        %% Optimization Initializations
        Num_agents  = 100;
        Max_iteration = 5;
        Rmin = 1e-5;
        any_reflect = any(reflect > 0) && any(reflect < 0);
        
        zeta_k_St = ones(1, Nr); 
        Active_Gain_dB = 0; 
        zeta_k_St = (10^(Active_Gain_dB/10)) * ones(1, Nr);
        
        phi_Sr = 2*pi*rand(Num_agents, Nr);
        phi_St = 2*pi*rand(Num_agents, Nr);
        alpha = rand(Num_agents, K+1); 
        alpha = alpha ./ sum(alpha, 2);      
        alpha = alpha - (sum(alpha,2)-1)/(K+1);
        alpha = alpha - (sum(alpha,2)-1)/(K+1);
        alpha = alpha - (sum(alpha,2)-1)/(K+1);
        AN_P_ratio = 1;  
        
        %% --- Alternating Optimization Execution ---
        max_AO_iter = Max_iteration;           
        max_SCA = 3;         
        tol = 1e-3;
        
        if any_reflect
            phi_Sr = 2*pi*rand(1, Nr);
            zeta_k_St = (10^(Active_Gain_dB/10)) * rand(Num_agents, Nr);
        else
            phi_Sr = zeros(1, Nr);
            zeta_k_St = (10^(Active_Gain_dB/10)) * ones(1, Nr);
        end
        
        manifold = complexcirclefactory(Nr, 1);
        problem = struct('M', manifold);
        num_agents  = 1;
        
        BW = delta_f;
        N0_dBm = -174;
        sigma2 = 10^((N0_dBm + 10*log10(BW) - 30)/10);
        Pw_dBm = 46;
        Pw = 10^((Pw_dBm - 30)/10);
        
        [L_node, E_node] = compute_channels(K, Nr, nF, Pe, P, Q_j, Plos, PLj, HB, HA, g_pq, Nsymb, reflect, h_rp, h_jq, h_e, Active_Gain_dB);    
        
        beta = zeros(Nr, num_agents);
        min_Rsec = zeros(num_agents, 1);
        for i = 1:num_agents
            beta(:, i) = manifold.rand();
            [R_sec, ~] = get_Secrecy_matrix(beta(:, i), L_node, E_node, alpha(1, :), K, nF, sigma2, Pw, AN_P_ratio);
            min_Rsec(i, 1) = min(min(R_sec));
        end
        
        b0 = beta(:, min_Rsec == max(max(min_Rsec)));
        for i = 1:num_agents
            beta(:, i) = manifold.rand();
            [R_sec, ~] = get_Secrecy_matrix(b0, L_node, E_node, alpha(i, :), K, nF, sigma2, Pw, AN_P_ratio);
            min_Rsec(i, 1) = min(min(R_sec));
        end
        
        alpha_prev = alpha(min_Rsec == max(max(min_Rsec)), :);
        alpha = alpha_prev;
        [R_sec_prev, rate_p, ~, ~] = get_Secrecy_matrix(b0, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);
        phi_St = wrapToPi(angle(b0)).';
        
        if any_reflect
            X = [alpha, phi_Sr, phi_St, zeta_k_St];
        else
            X = [alpha, phi_St(:).'];
        end
        
        [~, sc_p_lk, ~, ~, R_k, ~, ~, ~] = compute_sinr_sc_an(Pe, P, Q_j, nF+L, K, delta_f, Plos, PLj, Nr, HB, HA, g_pq, Nsymb, reflect, Rmin, h_rp, h_jq, h_e, zeta_k_St, Active_Gain_dB, AN_P_ratio, X);
        
        prev_cost = min(min(R_sec_prev));
        best_fake_secrecy = prev_cost;
        best_real_secrecy = min(min(sc_p_lk(nF+1:end, :)));
        Ck = max(0, Rmin - rate_p);
        feasible_flag = false;
        
        for ao = 1:max_AO_iter
            prev_fake = prev_cost;
          
            [phi_St, ~] = optimize_phi_manopt_fixed_alpha(Rmin, L_node, E_node, problem, b0, alpha, K, nF, sigma2, Pw, AN_P_ratio, Ck);
            b0 = exp(1i*phi_St(:));
            
            [cost, alpha_prev, Ck, feasible_flag, xi_val] = new_optimize_alpha_cvx_fixed_phi(Rmin, alpha_prev, L_node, E_node, phi_St, phi_Sr, zeta_k_St, K, nF, reflect, delta_f, Active_Gain_dB, AN_P_ratio, max_SCA);
            alpha = alpha_prev;
           
            if any_reflect
                X = [alpha, phi_Sr, phi_St, zeta_k_St];
            else
                X = [alpha, phi_St(:).'];
            end
            
            [~, sc_p_lk, ~, ~, R_k, ~, ~, ~] = compute_sinr_sc_an(Pe, P, Q_j, nF+L, K, delta_f, Plos, PLj, Nr, HB, HA, g_pq, Nsymb, reflect, Rmin, h_rp, h_jq, h_e, zeta_k_St, Active_Gain_dB, AN_P_ratio, X);
            
            if feasible_flag
                current_fake = min(min(sc_p_lk(1:nF, :)));
                current_real = min(min(sc_p_lk(nF+1:end, :)));
            else
                R_k = zeros(K, 1);
                current_fake = 0;
                current_real = 0;
                prev_min_Rk = 0;        
            end
            
            if feasible_flag && cost > prev_cost
                best_fake_secrecy = current_fake;
                best_real_secrecy = current_real;
                prev_cost = cost;
                prev_min_Rk = min(R_k);
            end
         
            fprintf(['AO Iter %2d | Feasible = %d | Fake Sec = %.6f | Δ = %.6f | ' ...
                     'max(xi)=%.2e |L = %2d| Nf_ratio = %2d | Ns=%2d\n'], ...
                    ao, feasible_flag, best_fake_secrecy, ...
                    best_fake_secrecy - prev_fake, max(xi_val), L, nF_ratio_vec(nFr_idx), mc_iter);
                
            if ~feasible_flag || abs(best_fake_secrecy)<1e-8 
                break;
            end
            
            if feasible_flag && abs(best_fake_secrecy - prev_fake) < tol && ao >= 5
                break;
            end
        end
        
        feasible_record(mc_iter, nFr_idx) = feasible_flag;
        Convex_min_Rk(mc_iter, nFr_idx) = prev_min_Rk;
        Convex_Convergence_curve_AO(mc_iter, nFr_idx) = prev_cost;
        Convex_Fake_Convergence_curve_AO(mc_iter, nFr_idx) = best_fake_secrecy;
        Convex_Real_Convergence_curve_AO(mc_iter, nFr_idx) = best_real_secrecy;
    end
end

%% ===================== PLOTTING & STRUCT SAVING =====================
colors = [0, 0.4470, 0.7410; 0.8500, 0.3250, 0.0980; 0.4660, 0.6740, 0.1880; 0.4940, 0.1840, 0.5560];
markerInterval = 1;
eps_val = 1e-12;

rel_diff = ((Convex_Real_Convergence_curve_AO - Convex_Fake_Convergence_curve_AO).^2) ./ max(abs(Convex_Real_Convergence_curve_AO), eps_val);
ratio = Convex_Fake_Convergence_curve_AO ./ max(Convex_Real_Convergence_curve_AO, eps_val);

valid_records_Qtd = sum(feasible_record, 1);
valid_records_Qtd_safe = max(valid_records_Qtd, 1);

Convex_min_Rk_mean = sum(Convex_min_Rk.*feasible_record,1)./valid_records_Qtd_safe;
Convex_Convergence_curve_AO_mean = sum(Convex_Convergence_curve_AO.*feasible_record,1)./valid_records_Qtd_safe; 
Convex_Fake_Convergence_curve_AO_mean = sum(Convex_Fake_Convergence_curve_AO.*feasible_record,1)./valid_records_Qtd_safe;  
Convex_Real_Convergence_curve_AO_mean = sum(Convex_Real_Convergence_curve_AO.*feasible_record,1)./valid_records_Qtd_safe; 
Convex_Real_Fake_diff_Convergence_curve_AO_mean = sum(rel_diff.*feasible_record,1)./valid_records_Qtd_safe; 
Convex_Real_Fake_ratio_Convergence_curve_AO_mean = sum(ratio .* feasible_record,1) ./ valid_records_Qtd_safe;
mean_diff = abs(Convex_Real_Convergence_curve_AO_mean - Convex_Fake_Convergence_curve_AO_mean);

figure('Color','w'); hold on;
plot(nF_ratio_vec, Convex_min_Rk_mean, '--s','Color',colors(1,:), 'LineWidth',1.5);
plot(nF_ratio_vec, Convex_Fake_Convergence_curve_AO_mean, '--s','Color',colors(2,:), 'LineWidth',1.5);
plot(nF_ratio_vec, Convex_Real_Convergence_curve_AO_mean, '-s','Color',colors(3,:), 'LineWidth',1.5);
xlabel('$V_E/E$','Interpreter','latex');
ylabel('Minimum SC (b/s/Hz)');
legend('Min Rate','Virtual SC','Real SC','Location','best');
grid on; box on;

NF_ratio_Results = struct();
NF_ratio_Results.Convex_min_Rk = Convex_min_Rk;
NF_ratio_Results.Convex_Convergence_curve_AO = Convex_Convergence_curve_AO;
NF_ratio_Results.Convex_Fake_Convergence_curve_AO = Convex_Fake_Convergence_curve_AO;
NF_ratio_Results.Convex_Real_Convergence_curve_AO = Convex_Real_Convergence_curve_AO;
NF_ratio_Results.feasible_record = feasible_record;
NF_ratio_Results.valid_records_Qtd = valid_records_Qtd;
save('N_F_ratio_Results.mat', 'NF_ratio_Results', '-v7.3');
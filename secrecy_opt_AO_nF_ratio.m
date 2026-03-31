clear; clc;
cvx_clear;

Ns = 10000; % number of samples for Monte Carlo simulation
%rng(3);

transmissionType = 'mc';

c = physconst('LightSpeed'); % speed of light
f_c = 10e9; % frequency
lambda = c / f_c; % wavelength
 % ---- Earth-centered conversion ----
R_earth = 6371e3;

% number of satellites
nSat = 2;

%Number of nodes
K_h = 1;  % number of high speed legit users % 50 km/h
K_s = 1;  % number of slow speed legit users % 1.2 m/s
K = K_h+K_s; % number of legit users





% --- OTFS System Parameters ---
delta_f = 100e3;      % Subcarrier spacing (Hz)
T = 1/delta_f;       % Symbol duration
B = 10e6;        % [Hz] ← Use this
Tf      = 14*T;      % 14-symbol frame (~1 ms)


L_tau = 8;   % 8 delay taps over max_tau (covers multipath + RIS)
L_nu  = 8;   % 8 Doppler taps over [-max_nu, max_nu]

R = 10;

m_rician = (R+1)^2/(2*R+1);

N_V = 20; % number of rows of regularly arranged unit cells of RIS
N_H = 20; % number of columns of regularly arranged unit cells of RIS
Nr = N_V * N_H; % total number of unit cells of RIS

d_x = floor(lambda/2 * 1000) / 1000; % horizontal size of RIS element
d_y = floor(lambda/2 * 1000) / 1000; % vertical size of RIS element


  %% Radiation patterns (Tang paper)
F      = @(theta,phi) cos(theta).^(1/64);
F_tx   = @(theta,phi) cos(theta).^1100;%33 2 degree
F_rx   = @(theta,phi) cos(theta).^15; %15 17.28 degree



%% Gains (Eq. 5–6)
G   = 4*pi / integral2(@(t,p) F(t,p).*sin(t),    0, pi/2, 0, 2*pi);
G_t = 4*pi / integral2(@(t,p) F_tx(t,p).*sin(t), 0, pi/2, 0, 2*pi);
G_r = 4*pi / integral2(@(t,p) F_rx(t,p).*sin(t), 0, pi/2, 0, 2*pi);


P = 2; % number of propagation paths arriving at the r-th RIS element
Q_j = 3; % number of propagation paths departing from the r-th RIS element
Pe = 2; % number of propagation paths departing from the LEO satellite

% location coordinates
HAP_altitude = 10e3;
LEO_altitude = 250e3;

RIS_normal = [0;0;1];


[max_theta,max_alpha,max_beta] = compute_max_theta(LEO_altitude,HAP_altitude);

S_elevation = pi/2-max_theta*(0.009);
S_azimuth  = pi/4;
S_radius  = R_earth+LEO_altitude;


S_sph = [ S_azimuth, S_elevation, S_radius]; % location of LEO satellite in spherical coordinates

AN_elevation = pi/2-max_theta*(0.1);
AN_azimuth  = pi/4 -pi;
AN_radius  = S_radius;


AN_sph = [ AN_azimuth, AN_elevation, AN_radius]; % location of AN LEO satellite in spherical coordinates

[S_x, S_y, S_z] = sph2cart(S_sph(1), S_sph(2), S_sph(3));
S_xyz = [S_x;S_y;S_z]; %[-1.4e4,-1.4e4,6.6e6]


[AN_x, AN_y, AN_z] = sph2cart(AN_sph(1), AN_sph(2), AN_sph(3));
AN_xyz = [AN_x;AN_y;AN_z];


S_v = 7800;
R_xyz = [0; 0; R_earth+HAP_altitude]; % location of STAR-RIS; code assumes this to be origin;
% note: this code assumes surface is on x-y plane (surface normal points in
% z-axis direction)



RIS_size_x = N_H * d_x;
RIS_size_y = N_V * d_y;

minAngle=(lambda/max(RIS_size_x,RIS_size_y));

receiving_ang = acos( ...
    dot( ...
        (AN_xyz - R_xyz) / norm(AN_xyz - R_xyz,'fro'), ...
        (S_xyz  - R_xyz) / norm(S_xyz  - R_xyz,'fro') ...
    ) ...
);

% Nakagami Parameters:

% from LEO satellite to STAR-RIS (one value for each path and assume same
% for each RIS element):
m_p = [m_rician;1*ones(P-1,1)]; % shape parameter
omega_p = (1/P)*ones(1,P); % spread parameter

nF_ratio_vec = 0.5:0.5:6;

% Convex_min_Rk= zeros(Ns,10,20);
% Convex_Convergence_curve_AO = zeros(Ns,10,20);
% Convex_Fake_Convergence_curve_AO = zeros(Ns,10,20);
% Convex_Real_Convergence_curve_AO = zeros(Ns,10,20);f


for mc_iter = 1:Ns
for nFr_idx = 1:length(nF_ratio_vec)

    L = 2*randi([1, 3]); % number of eavesdroppers
    nF = nF_ratio_vec(nFr_idx)*L;
    

% from STAR-RIS to users and eavesdroppers (one value for each path and 
% each receiver, and assume same for each RIS element):
m_q = 1*ones(K+nF+L,Q_j); % shape parameter
omega_q = (1/Q_j)*ones(K+nF+L,Q_j); % spread parameter
% note: m_j_1(1:K,:) and omega_j_q(1:K,:) are for legit users
% and m_j_1(K+1:end,:) and omega_j_q(K+1:end,:) are for eavesdroppers
% direct channel from LEO satellite to Eve (one value for each path and
% each eavesdropper)

m_e = 1*ones(K+nF+L,Pe); % shape parameter
omega_e = (1/Pe)*ones(K+nF+L,Pe); % spread parameter

% ============================================================
% Earth-centered positions and velocities (LEO → RIS)
% ============================================================

% Orbit geometry
orbit_normal = [1; 0; 0];
Rs = norm(S_xyz);

omega_orb = (S_v / Rs) * orbit_normal;

% Satellite velocity (ECI)
vS = cross(omega_orb, S_xyz);

% RIS velocity due to Earth rotation (ECI)
vR = [0;0;0];

sigma_ang = deg2rad(30);   % angular spread

g_pq = zeros(P,Q_j,K+nF+L);
Plos = zeros(K+nF+L,nSat);
PLj = zeros(K+nF+L,nSat);

h_rp = zeros(Nr, P,K+nF+L,nSat);
h_jq = zeros(Nr, Q_j,K+nF+L);
h_e = zeros(Pe,K+nF+L,nSat);
taus_ku = zeros(Pe,K);
nus_ku = zeros(Pe,K);


% ELEVATION (UNCHANGED)
el = pi/2 - (0.05)*max_alpha * rand(1, K);


az = (2*pi) * rand(1, K);


% RADIUS (FIXED ON EARTH SURFACE)
r = R_earth * ones(1, K);

% SPHERICAL → CARTESIAN
[x, y, z] = sph2cart(az, el, r);

ground_users_cart = [x; y; z];   % 3 x K

%% Eavesdroppers position
x = 1000*rand(1, L)+1000;
y = 1000*rand(1, L)+1000;
z = R_earth + 50 + 950*rand(1, L);


%% Fake Eavesdroppers position
x_f = 1000*rand(1, nF)+1000;
y_f = 1000*rand(1, nF)+1000;
z_f = R_earth + 50 + 950*rand(1, nF);

%% SHIFT TO GLOBAL COORDINATES

eavesdroppers_xyz =  [x; y; z];
fake_eavesdroppers_xyz =  [x_f; y_f; z_f];

rho_j_xyz = [ground_users_cart,eavesdroppers_xyz,fake_eavesdroppers_xyz];

% find out whether each receiver is on the reflect side or transmit side
reflect = sign(RIS_normal.' * (rho_j_xyz - R_xyz));




% Satellite to RIS delays and doppler coefficients.
[taus_R, nus_R, u_paths_R] = compute_delay_and_doppler( ...
    c, S_xyz, vS, R_xyz, vR, f_c, P, sigma_ang);

%Channels
for p = 1:P
    % Scale parameter theta = omega/m
    h_rp(:, p,1) = sqrt(gamrnd(m_p(p), omega_p(p)/m_p(p), [Nr, 1])) .* exp(1i*2*pi*rand(Nr, 1)); % Data carrying channel
    h_rp(:, p,2) = sqrt(gamrnd(m_p(p), omega_p(p)/m_p(p), [Nr, 1])) .* exp(1i*2*pi*rand(Nr, 1)); % Noise carrying channel
end




%OTFS Per User channel conditions
for k =1:K

  
    if k <= K_h
        vk_ms = 50 / 3.6;
    else
        vk_ms = 1.2;
    end


    User_k_loc = rho_j_xyz(:,k);

    d_los = sqrt(sum((User_k_loc-S_xyz).^2));
    d_los_AN = sqrt(sum((User_k_loc-AN_xyz).^2));
    % d_ris = sqrt(sum((User_k_loc-R_xyz).^2));


    Plos(k,1) = ((lambda/(4*pi))^2*G_t*G_r)/(d_los^2);
    Plos(k,2) = ((lambda/(4*pi))^2*G_t*G_r)/(d_los_AN^2);
    % Plos_ris(k) = ((lambda/(4*pi))^2*G_t*G_r)/(d_ris^2);
    % 
    % PL = Plos(k)*Plos_ris(k) ;
    
    PLj(k,1) = compute_ris_PL(lambda,N_V,N_H,S_xyz,User_k_loc,R_xyz,RIS_normal,F,F_tx,F_rx,G,G_t,G_r);
    PLj(k,2) = compute_ris_PL(lambda,N_V,N_H,AN_xyz,User_k_loc,R_xyz,RIS_normal,F,F_tx,F_rx,G,G_t,G_r);
   
  

     % direction in xy (away from ris)
    if User_k_loc(1)==R_xyz(1) && User_k_loc(2)==R_xyz(2)
        d_ru = [R_xyz(1)+1;R_xyz(1)+1];
    else
        d_ru = User_k_loc(1:2) - R_xyz(1:2);
    end
    d_ru = d_ru / norm(d_ru);


    % receiver velocity (guaranteed norm)
    v_l  = vk_ms * [d_ru; 0];

    % RIS to legitimate users delays and doppler coefficients.
    [taus_k, nus_k, u_paths_k] = compute_delay_and_doppler( ...
    c, R_xyz, vR, User_k_loc, v_l, f_c, Q_j, sigma_ang);
  

    for p=1:P
        for q=1:Q_j
           taus_kq(k,p,q) = taus_R(p)+taus_k(q);
           nus_kq(k,p,q) = nus_R(p)+nus_k(q);
        end 
    end
    

    % sat to legitimate users delays and doppler coefficients.
    [taus_u, nus_u, u_paths_u] = compute_delay_and_doppler( ...
    c, S_xyz, vS, User_k_loc, v_l, f_c, Pe, sigma_ang);

    g_pq(:,:,k) = exp(1i*2*pi*(taus_R*nus_k'));    

    taus_ku(:,k) = taus_u; 
    nus_ku(:,k) = nus_u; 
   
     
    for q = 1:Q_j
        h_jq(:, q,k) = sqrt(gamrnd(m_q(k,q), omega_q(q)/m_q(k,q), [Nr, 1])) .* exp(1i*2*pi*rand(Nr,1));
    end

    
    for u = 1:Pe
        h_e(u,k,1) = sqrt(gamrnd(m_e(k,u), omega_e(u)/m_e(k,u))) .* exp(1i*2*pi*rand());% Direct Data carrying channel
        h_e(u,k,2) = sqrt(gamrnd(m_e(k,u), omega_e(u)/m_e(k,u))) .* exp(1i*2*pi*rand());% Direct Noise carrying channel
    end
            
end

% % % % Max tau and nu
% max_tau = max([taus_kq(:);taus_ku(:)])-min([taus_kq(:);taus_ku(:)]); 
% max_nu  = max([nus_kq(:);nus_ku(:)])-min([nus_kq(:);nus_ku(:)]);    


% % Compute M and N based on the parameters
% [M, N] = computeOTFSgrid(max_tau, max_nu, 'numerology', B, delta_f, T, Tf);
% M = max(M, 64); N = max(N, 20);  % Minimum practical size

M = 18;
N = 18;

Nsymb = M*N; 

HA = zeros(Nsymb,Nsymb,P,Q_j,K+nF+L); % Relay link
HB = zeros(Nsymb,Nsymb,Pe,K+nF+L);


for k=1:K

    for p=1:P
        for q=1:Q_j
            HA(:,:,p,q,k) =  compute_Hp(taus_kq(k,p,q), nus_kq(k,p,q), M, N, T, delta_f, 'blocked',transmissionType);
        end
    end
   

    for u = 1:Pe
        HB(:,:,u,k) =  compute_Hp(taus_u(u), nus_u(u), M, N, T, delta_f, 'blocked',transmissionType);
    end

end



for l=1:nF+L

        vl_ms = 0;        

        User_l_loc = rho_j_xyz(:,K+l);

        d_los = sqrt(sum((User_l_loc-S_xyz).^2));
        d_los_AN = sqrt(sum((User_l_loc-AN_xyz).^2));

        Plos(K+l,1) = ((lambda/(4*pi))^2*G_t*G_r)/(d_los^2);
        Plos(K+l,2) = ((lambda/(4*pi))^2*G_t*G_r)/(d_los_AN^2);


        PLj(K+l,1) = compute_ris_PL(lambda,N_V,N_H,S_xyz,User_l_loc,R_xyz,RIS_normal,F,F_tx,F_rx,G,G_t,G_r);
        PLj(K+l,2) = compute_ris_PL(lambda,N_V,N_H,AN_xyz,User_l_loc,R_xyz,RIS_normal,F,F_tx,F_rx,G,G_t,G_r);
    

        % RIS → UAV radial direction
        u_re = User_l_loc - R_xyz;
        u_re = u_re / norm(u_re);
        
        % rotation axis (choose vertical)
        a_hat = [0;0;1];
        if abs(dot(a_hat,u_re)) > 0.99
            a_hat = [1;0;0];
        end
        
        % tangential velocity (motion on sphere)
        v_dir = cross(a_hat, u_re);
        v_dir = v_dir / norm(v_dir);
        
        % UAV velocity
        v_l = vl_ms * v_dir;
    
        % RIS to eavesdropper users delays and doppler coefficients.
        [taus_l, nus_l, u_paths_l] = compute_delay_and_doppler( ...
        c, R_xyz, vR, User_l_loc, v_l, f_c, Q_j, sigma_ang);
    
        % sat to eavesdropper users users delays and doppler coefficients.
        [taus_u_l, nus_u_l, u_paths_u] = compute_delay_and_doppler( ...
        c, S_xyz, vS, User_l_loc, v_l, f_c, Pe, sigma_ang);         
    
        
        for q = 1:Q_j
            h_jq(:, q,K+l) = sqrt(gamrnd(m_q(K+l,q), omega_q(q)/m_q(K+l,q), [Nr, 1])) .* exp(1i*2*pi*rand(Nr,1));
        end
    
        
        for u = 1:Pe
            h_e(u,K+l,1) = sqrt(gamrnd(m_e(K+l,u), omega_e(u)/m_e(K+l,u))) .* exp(1i*2*pi*rand());
            h_e(u,K+l,2) = sqrt(gamrnd(m_e(K+l,u), omega_e(u)/m_e(K+l,u))) .* exp(1i*2*pi*rand());
        end
    

       for p=1:P
            for q=1:Q_j
                HA(:,:,p,q,K+l) =  compute_Hp(taus_R(p)+taus_l(q), nus_R(p)+nus_l(q), M, N, T, delta_f, 'blocked',transmissionType);
            end
       end

        for u = 1:Pe
            HB(:,:,u,K+l) =  compute_Hp(taus_u_l(u) ,nus_u_l(u), M, N, T, delta_f, 'blocked',transmissionType);
        end

        g_pq(:,:,K+l) = exp(1i*2*pi*(taus_R*nus_l'));           
end



% if gpuDeviceCount > 0
%     % GPU is available
%     HA = gpuArray(HA); HB = gpuArray(HB);  % After precompute  
%     h_rp = gpuArray(h_rp); h_jq = gpuArray(h_jq); g_pq = gpuArray(g_pq); h_e = gpuArray(h_e);
% end
 
%% Sine-Cosine optimization

display('SCA is optimizing your problem');

Num_agents  = 100;
Max_iteration = 5;
Rmin=1e-5;

% Check if more than one STAR-RIS side is being used.
any_reflect = any(reflect > 0) && any(reflect < 0);

% Problem bounds and dimensionality
dim = K+1+Nr;
ub=[ones(1,K+1),2*pi*ones(1,Nr)];
alpha_min = 1e-4;
lb = [alpha_min * ones(1,K+1),zeros(1,Nr)];
zeta_k_St = ones(1,Nr); % RIS amplitude coefficients, we may use it to boost for active RIS

Active_Gain_dB = 0; 
zeta_k_St = (10^(Active_Gain_dB/10)) * ones(1, Nr);


% zeta_k_Sr = rand(Num_agents,Nr); % reflection coefficients
phi_Sr = 2*pi*rand(Num_agents,Nr);
phi_St = 2*pi*rand(Num_agents,Nr);% transmission phases


alpha = rand(Num_agents, K+1); 
% random values
alpha = alpha ./ sum(alpha, 2);      % divide each row by its row sum

alpha = alpha - (sum(alpha,2)-1)/(K+1);
alpha = alpha - (sum(alpha,2)-1)/(K+1);
alpha = alpha - (sum(alpha,2)-1)/(K+1);


% X = [alpha,phi_St];
% 
% if any_reflect
%     dim = K+1+3*Nr;
%     ub=[ones(1,K+1),2*pi*ones(1,2*Nr),ones(1,Nr)];
%     alpha_min = 1e-4;
%     lb = [alpha_min * ones(1,K+1),zeros(1,3*Nr)];
%     zeta_k_St = (10^(Active_Gain_dB/10)) *rand(Num_agents,Nr);
%     X = [alpha,phi_Sr,phi_St,zeta_k_St];
% end
% 
% 
% % --- Problem Dimensions and Bounds ---
% dim_pso = dim;
% alpha_min_pso = alpha_min;
% lb_pso =lb;
% ub_pso = ub;

AN_P_ratio = 1;  



%% ===================== CONVEX ALTERNATING OPTIMIZATION (AO) =====================
display('Convex Approximation with AO');

max_AO_iter = Max_iteration;           % Outer AO iterations
max_SCA = 3;         % Inner SCA iterations for alpha subproblem
tol = 1e-3;

Active_Gain_dB = 0; 


if any_reflect
    phi_Sr = 2*pi*rand(1,Nr);
    zeta_k_St = (10^(Active_Gain_dB/10)) *rand(Num_agents,Nr);
else
    phi_Sr = zeros(1,Nr);
    zeta_k_St = (10^(Active_Gain_dB/10)) * ones(1, Nr);
end


Convergence_curve_AO = zeros(1, max_AO_iter);

fprintf('\n=== Starting Convex AO ===\n');

manifold = complexcirclefactory(Nr,1);
problem.M = manifold;
num_agents  = 1;



%Parameters
BW = delta_f;
N0_dBm = -174;
sigma2 = 10^((N0_dBm + 10*log10(BW) - 30)/10);
% Define the scaling factor
Pw_dBm = 46;
Pw = 10^((Pw_dBm - 30)/10);



 % ---------- CHANNELS ----------
[L_node,E_node] = compute_channels( K, Nr, nF, Pe, P, Q_j, Plos, PLj, HB, HA, g_pq, Nsymb, ...
reflect, h_rp, h_jq, h_e,  Active_Gain_dB);    


beta = zeros(Nr,num_agents);

min_Rsec = zeros(num_agents,1);
for i = 1:num_agents
    beta(:,i) = manifold.rand();
   [R_sec,~] = get_Secrecy_matrix(beta(:,i), L_node, E_node, alpha(1,:), K, nF, sigma2, Pw, AN_P_ratio);
    min_Rsec(i,1) = min(min(R_sec));
end

b0 =  beta(:,min_Rsec == max(max(min_Rsec)));


for i = 1:num_agents
    beta(:,i) = manifold.rand();
   [R_sec,~] = get_Secrecy_matrix(b0, L_node, E_node, alpha(i,:), K, nF, sigma2, Pw, AN_P_ratio);
    min_Rsec(i,1) = min(min(R_sec));
end

alpha_prev = alpha(min_Rsec == max(max(min_Rsec)),:);

alpha = alpha_prev;

[R_sec_prev,rate_p,rate_c,~] = get_Secrecy_matrix(b0, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);
phi_St = wrapToPi(angle(b0)).';

 % ================================================================
    % Build X
    % ================================================================
    if any_reflect
        X = [alpha, phi_Sr, phi_St, zeta_k_St];
    else
        X = [alpha, phi_St(:).'];
    end

[~, sc_p_lk, ~, ~, ~, R_k,sinr_c_k, sinr_p_k, ~] = compute_sinr_sc_an(...
        Pe, P, Q_j, nF+L, K, delta_f, Plos, PLj, Nr, HB, HA, g_pq, ...
        Nsymb, reflect, Rmin, h_rp, h_jq, h_e, ...
        zeta_k_St, Active_Gain_dB,AN_P_ratio, X);

prev_cost = min(min(R_sec_prev));
best_fake_secrecy = prev_cost;
best_real_secrecy = min(min(sc_p_lk(nF+1:end,:)));

Ck = max(0, Rmin - rate_p);

feasible_ao = true;      % track feasibility of this realization

for ao = 1:max_AO_iter

    prev_fake = prev_cost;
  

   

    % ================================================================
    % 2. SUBPROBLEM 2: Optimize RIS Phases Φ
    % ================================================================
    [phi_St,cost_opt] = optimize_phi_manopt_fixed_alpha(...
        Rmin,L_node,E_node,problem,b0,alpha,K, nF, sigma2, Pw, AN_P_ratio,Ck);

    b0 = exp(1i*phi_St(:));

    % [R_sec_next,~] = get_Secrecy_matrix(b0, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);
    % 
    % min_next = min(min(R_sec_next));


    % ================================================================
    % 1. SUBPROBLEM 1: Optimize Power Allocation α
    % ================================================================
    [cost,alpha_prev,Ck,feasible_flag,xi_val] = new_optimize_alpha_cvx_fixed_phi(...
        Rmin,alpha_prev,L_node,E_node,phi_St, phi_Sr, zeta_k_St, ...
        K, nF, reflect, delta_f, Active_Gain_dB,AN_P_ratio, max_SCA);

    alpha = alpha_prev;
   

   [R_sec_next2,~] = get_Secrecy_matrix(b0, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);
    
    cost_opt = min(min(R_sec_next2));

    % ================================================================
    % Build X
    % ================================================================
    if any_reflect
        X = [alpha, phi_Sr, phi_St, zeta_k_St];
    else
        X = [alpha, phi_St(:).'];
    end

    % ================================================================
    % Final evaluation
    % ================================================================
    [~, sc_p_lk, ~, ~, ~, R_k,sinr_c_k, sinr_p_k, ~] = compute_sinr_sc_an(...
        Pe, P, Q_j, nF+L, K, delta_f, Plos, PLj, Nr, HB, HA, g_pq, ...
        Nsymb, reflect, Rmin, h_rp, h_jq, h_e, ...
        zeta_k_St, Active_Gain_dB,AN_P_ratio, X);

    rate_p_vec = log2(1 + sinr_p_k);

    % ================================================================
    % Handle feasibility properly
    % ================================================================
    if feasible_flag
        Rk = rate_p_vec(:) + Ck;

        current_fake = min(min(sc_p_lk(1:nF,:)));
        current_real = min(min(sc_p_lk(nF+1:end,:)));
    else
        % Treat as outage
        Rk = zeros(K,1);
        current_fake = 0;
        current_real = 0;
        prev_min_Rk=0;        
    end

    % ================================================================
    % Update best solution
    % ================================================================
    if feasible_flag && cost > prev_cost
        best_fake_secrecy = current_fake;
        best_real_secrecy = current_real;
        Destination_position = X;
        prev_cost = cost;
        prev_min_Rk = min(Rk);
    end

 
 
    % ================================================================
    % Logging
    % ================================================================
    fprintf(['AO Iter %2d | Feasible = %d | Fake Sec = %.6f | Δ = %.6f | ' ...
             'max(xi)=%.2e |L = %2d| Nf_ratio = %2d | Ns=%2d\n'], ...
            ao, feasible_flag, best_fake_secrecy, ...
            best_fake_secrecy - prev_fake, max(xi_val),L, nF_ratio_vec(nFr_idx), mc_iter);

    if ~feasible_flag || abs(best_fake_secrecy)<1e-8 
        break;
    end

    
    % ================================================================
    % Convergence
    % ================================================================
    if feasible_flag && abs(best_fake_secrecy - prev_fake) < tol && ao >= 5
        fprintf('→ AO Converged at iteration %d\n', ao);
        break;
    end

end

 % Track feasibility
 feasible_record(mc_iter,nFr_idx) = feasible_flag;

 Convex_min_Rk(mc_iter,nFr_idx) = prev_min_Rk;
 Convex_Convergence_curve_AO(mc_iter,nFr_idx) = prev_cost;
 Convex_Fake_Convergence_curve_AO(mc_iter,nFr_idx) = best_fake_secrecy;
 Convex_Real_Convergence_curve_AO(mc_iter,nFr_idx) = best_real_secrecy;

fprintf('\nConvex AO Finished! Best Fake Secrecy Rate = %.8f\n', best_fake_secrecy);

end
end



%% Plot
figure('Color','w');

colors = [0, 0.4470, 0.7410;      
          0.8500, 0.3250, 0.0980; 
          0.4660, 0.6740, 0.1880; 
          0.4940, 0.1840, 0.5560];

markerInterval = 1;
eps_val = 1e-12;

% ================= NORMALIZATION =================
rel_diff = ((Convex_Real_Convergence_curve_AO - ...
            Convex_Fake_Convergence_curve_AO).^2) ./ ...
            max(abs(Convex_Real_Convergence_curve_AO), eps_val);

ratio = Convex_Fake_Convergence_curve_AO ./ ...
        max(Convex_Real_Convergence_curve_AO, eps_val);

% ================= FEASIBILITY =================
valid_records_Qtd = sum(feasible_record, 1);
valid_records_Qtd_safe = max(valid_records_Qtd, 1);

% ================= MEANS =================
Convex_min_Rk_mean = sum(Convex_min_Rk.*feasible_record,1)./valid_records_Qtd_safe;
Convex_Convergence_curve_AO_mean = sum(Convex_Convergence_curve_AO.*feasible_record,1)./valid_records_Qtd_safe; 
Convex_Fake_Convergence_curve_AO_mean = sum(Convex_Fake_Convergence_curve_AO.*feasible_record,1)./valid_records_Qtd_safe;  
Convex_Real_Convergence_curve_AO_mean = sum(Convex_Real_Convergence_curve_AO.*feasible_record,1)./valid_records_Qtd_safe; 
Convex_Real_Fake_diff_Convergence_curve_AO_mean = sum(rel_diff.*feasible_record,1)./valid_records_Qtd_safe; 
Convex_Real_Fake_ratio_Convergence_curve_AO_mean = sum(ratio .* feasible_record,1) ./ valid_records_Qtd_safe;

mean_diff = (Convex_Real_Convergence_curve_AO_mean - Convex_Fake_Convergence_curve_AO_mean).^2;

x = 1:length(Convex_Convergence_curve_AO_mean);

%% ================= FIGURE 1: Convergence =================
figure('Color','w'); hold on;

plot(x, Convex_Convergence_curve_AO_mean, '-o', ...
    'Color', colors(3,:), 'LineWidth',2, ...
    'MarkerIndices',1:markerInterval:length(x), ...
    'MarkerFaceColor',colors(3,:));

title('Convergence Curve','FontWeight','bold');
xlabel('Nf');
ylabel('Best Fake Secrecy Rate');
grid on; box on;

% saveas(gcf,'ConvergenceCurve.png');
% savefig('ConvergenceCurve.fig');

%% ================= FIGURE 2: Real vs Fake =================
figure('Color','w'); hold on;

plot(nF_ratio_vec, Convex_min_Rk_mean, '--s','Color',colors(1,:), 'LineWidth',1.5);
plot(nF_ratio_vec, Convex_Fake_Convergence_curve_AO_mean, '--s','Color',colors(2,:), 'LineWidth',1.5);
plot(nF_ratio_vec, Convex_Real_Convergence_curve_AO_mean, '-^','Color',colors(3,:), 'LineWidth',1.5);

title('Real vs Fake Secrecy Rate');
xlabel('Number of Fake Eves');
ylabel('Secrecy Rate (b/s/Hz)');
legend('Min Rate','Fake','Real','Location','best');

grid on; box on;

% saveas(gcf,'Real_vs_Fake.png');
% savefig('Real_vs_Fake.fig');

%% ================= FIGURE 3: Relative Difference =================
figure('Color','w'); hold on;

plot(nF_ratio_vec, Convex_Real_Fake_diff_Convergence_curve_AO_mean, '-d', ...
    'Color', colors(4,:), 'LineWidth',1.8);

yline(0,'--k');

title('Relative Difference (Real - Fake) / Real');
xlabel('Number of Fake Eves ratio');
ylabel('Relative Gap');
grid on; box on;

% saveas(gcf,'Relative_Difference.png');
% savefig('Relative_Difference.fig');

%% ================= FIGURE 4: Ratio =================
figure('Color','w'); hold on;

plot(nF_ratio_vec, Convex_Real_Fake_ratio_Convergence_curve_AO_mean, '-o', ...
    'Color', colors(1,:), 'LineWidth',1.8);

yline(1,'--k');

title('Performance Ratio (Fake / Real)');
xlabel('Fake/Ratio Eves ratio');
ylabel('Ratio');
grid on; box on;

% saveas(gcf,'Ratio.png');
% savefig('Ratio.fig');

%% ================= FIGURE 5: Mean diff =================
figure('Color','w'); hold on;

plot(nF_ratio_vec, mean_diff, '-o', ...
    'Color', colors(1,:), 'LineWidth',1.8);

yline(1,'--k');

title('(Mean real - Mean fake) squared');
xlabel('Fake/Ratio Eves ratio');
ylabel('Ratio');
grid on; box on;

%% ================= SAVE STRUCT =================
NF_ratio_Results = struct();

NF_ratio_Results.Convex_min_Rk = Convex_min_Rk;
NF_ratio_Results.Convex_Convergence_curve_AO = Convex_Convergence_curve_AO;
NF_ratio_Results.Convex_Fake_Convergence_curve_AO = Convex_Fake_Convergence_curve_AO;
NF_ratio_Results.Convex_Real_Convergence_curve_AO = Convex_Real_Convergence_curve_AO;
NF_ratio_Results.feasible_record = feasible_record;

NF_ratio_Results.valid_records_Qtd = valid_records_Qtd;
NF_ratio_Results.Convex_min_Rk_mean = Convex_min_Rk_mean;
NF_ratio_Results.Convex_Convergence_curve_AO_mean = Convex_Convergence_curve_AO_mean;
NF_ratio_Results.Convex_Fake_Convergence_curve_AO_mean = Convex_Fake_Convergence_curve_AO_mean;
NF_ratio_Results.Convex_Real_Convergence_curve_AO_mean = Convex_Real_Convergence_curve_AO_mean;
NF_ratio_Results.rel_diff_mean = Convex_Real_Fake_diff_Convergence_curve_AO_mean;
NF_ratio_Results.ratio_mean = Convex_Real_Fake_ratio_Convergence_curve_AO_mean;
NF_ratio_Results.mean_diff = mean_diff;

save('N_F_ratio_Results.mat', 'NF_ratio_Results', '-v7.3');
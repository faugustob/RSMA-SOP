clear; clc;
cvx_clear;

Ns = 200; % number of samples for Monte Carlo simulation
%rng(1);

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

nF = 4; % Number of fake eavesdroppers

L = 1; % number of eavesdroppers

% --- OTFS System Parameters ---
delta_f = 100e3;      % Subcarrier spacing (Hz)
T = 1/delta_f;       % Symbol duration
B = 10e6;        % [Hz] ← Use this
Tf      = 14*T;      % 14-symbol frame (~1 ms)


L_tau = 8;   % 8 delay taps over max_tau (covers multipath + RIS)
L_nu  = 8;   % 8 Doppler taps over [-max_nu, max_nu]

R = 10;

m_rician = (R+1)^2/(2*R+1);



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




% Nakagami Parameters:

% from LEO satellite to STAR-RIS (one value for each path and assume same
% for each RIS element):
m_p = [m_rician;1*ones(P-1,1)]; % shape parameter
omega_p = (1/P)*ones(1,P); % spread parameter

nV_vec = 20:1:30;
nH_vec = 20:2:40;
N_H=20;

% Convex_min_Rk= zeros(Ns,10,20);
% Convex_Convergence_curve_AO = zeros(Ns,10,20);
% Convex_Fake_Convergence_curve_AO = zeros(Ns,10,20);
% Convex_Real_Convergence_curve_AO = zeros(Ns,10,20);

for mc_iter = 1:Ns
for nV_idx = 1:length(nV_vec)
%for nH_idx = 1:length(nH_vec)
    N_V = nV_vec(nV_idx);
    %N_H = nH_vec(nH_idx);

%N_V = 20; % number of rows of regularly arranged unit cells of RIS
%N_H = 20; % number of columns of regularly arranged unit cells of RIS
Nr = N_V * N_H; % total number of unit cells of RIS



RIS_size_x = N_H * d_x;
RIS_size_y = N_V * d_y;

minAngle=(lambda/max(RIS_size_x,RIS_size_y));

receiving_ang = acos( ...
    dot( ...
        (AN_xyz - R_xyz) / norm(AN_xyz - R_xyz,'fro'), ...
        (S_xyz  - R_xyz) / norm(S_xyz  - R_xyz,'fro') ...
    ) ...
);

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
   
    %Channels
    for p = 1:P
        % Scale parameter theta = omega/m
        h_rp(:, p,k,1) = sqrt(gamrnd(m_p(p), omega_p(p)/m_p(p), [Nr, 1])) .* exp(1i*2*pi*rand(Nr, 1)); % Data carrying channel
        h_rp(:, p,k,2) = sqrt(gamrnd(m_p(p), omega_p(p)/m_p(p), [Nr, 1])) .* exp(1i*2*pi*rand(Nr, 1)); % Noise carrying channel
    end

    
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

M = 100;
N = 16;

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


        %Channels
        for p = 1:P
            % Scale parameter theta = omega/m
            h_rp(:, p,K+l,1) = sqrt(gamrnd(m_p(p), omega_p(p)/m_p(p), [Nr, 1])) .* exp(1i*2*pi*rand(Nr, 1));
            h_rp(:, p,K+l,2) = sqrt(gamrnd(m_p(p), omega_p(p)/m_p(p), [Nr, 1])) .* exp(1i*2*pi*rand(Nr, 1));
        end
    
    
        
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



if gpuDeviceCount > 0
    % GPU is available
    HA = gpuArray(HA); HB = gpuArray(HB);  % After precompute  
    h_rp = gpuArray(h_rp); h_jq = gpuArray(h_jq); g_pq = gpuArray(g_pq); h_e = gpuArray(h_e);
end
 
%% Sine-Cosine optimization

display('SCA is optimizing your problem');

Num_agents  = 100;
Max_iteration = 20;
Rmin=1e-2;

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
tol = 1e-6;

Active_Gain_dB = 0; 


if any_reflect
    phi_Sr = 2*pi*rand(1,Nr);
    zeta_k_St = (10^(Active_Gain_dB/10)) *rand(Num_agents,Nr);
else
    phi_Sr = zeros(1,Nr);
    zeta_k_St = (10^(Active_Gain_dB/10)) * ones(1, Nr);
end

best_fake_secrecy = -5;
best_real_secrecy = -5;
Convergence_curve_AO = zeros(1, max_AO_iter);

fprintf('\n=== Starting Convex AO ===\n');

manifold = complexcirclefactory(Nr,1);
problem.M = manifold;
num_agents  = Num_agents;
num_agents  = 1;

prev_cost  = 10;

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


for i = 1:num_agents
    beta(:,i) = manifold.rand();
   [R_sec,rate_p,~] = get_Secrecy_matrix(beta(:,i), L_node, E_node, alpha(1,:), K, nF, sigma2, Pw, AN_P_ratio);
    min_Rsec(i,1) = min(min(R_sec));
end

b0 =  beta(:,min_Rsec == max(max(min_Rsec)));


for i = 1:num_agents   
   [R_sec,rate_p,~] = get_Secrecy_matrix(b0, L_node, E_node, alpha(i,:), K, nF, sigma2, Pw, AN_P_ratio);
    min_Rsec(i,1) = min(min(R_sec));
end

alpha_prev = alpha(min_Rsec == max(max(min_Rsec)),:);

alpha = alpha_prev;

[R_sec_prev,rate_p,rate_c,~] = get_Secrecy_matrix(b0, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);
phi_St = wrapToPi(angle(b0)).';

min_prev = min(min(R_sec_prev));

Ck = max(0, Rmin - rate_p);

feasible_ao = true;      % track feasibility of this realization

for ao = 1:max_AO_iter

    prev_fake = best_fake_secrecy;
  

   

    % ================================================================
    % 2. SUBPROBLEM 2: Optimize RIS Phases Φ
    % ================================================================
    [phi_St,cost_opt] = optimize_phi_manopt_fixed_alpha(...
        Rmin,L_node,E_node,problem,b0,alpha,K, nF, sigma2, Pw, AN_P_ratio,Ck);

    b0 = exp(1i*phi_St(:));

    [R_sec_next,~] = get_Secrecy_matrix(b0, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);
    
    min_next = min(min(R_sec_next));


    % ================================================================
    % 1. SUBPROBLEM 1: Optimize Power Allocation α
    % ================================================================
    [alpha_prev,Ck,feasible_flag,xi_val] = new_optimize_alpha_cvx_fixed_phi(...
        Rmin,alpha_prev,L_node,E_node,phi_St, phi_Sr, zeta_k_St, ...
        K, nF, reflect, delta_f, Active_Gain_dB,AN_P_ratio, max_SCA);

    alpha = alpha_prev;

    % Track feasibility
    xi_record(mc_iter,ao) = feasible_flag;

   [R_sec_next2,~] = get_Secrecy_matrix(b0, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);
    
    min_next2 = min(min(R_sec_next2));

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
    end

    % ================================================================
    % Update best solution
    % ================================================================
    if feasible_flag && cost_opt < prev_cost
        best_fake_secrecy = current_fake;
        best_real_secrecy = current_real;
        Destination_position = X;
        prev_cost = cost_opt;
        prev_min_Rk = min(Rk);
    end



    % ================================================================
    % Logging
    % ================================================================
    fprintf(['AO Iter %2d | Feasible = %d | Fake Sec = %.6f | Δ = %.6f | ' ...
             'max(xi)=%.2e | Ns=%2d\n'], ...
            ao, feasible_flag, best_fake_secrecy, ...
            best_fake_secrecy - prev_fake, max(xi_val), mc_iter);

    % ================================================================
    % Convergence
    % ================================================================
    if feasible_flag && abs(best_fake_secrecy - prev_fake) < tol && ao >= 5
        fprintf('→ AO Converged at iteration %d\n', ao);
        break;
    end

end

Convex_min_Rk(mc_iter,nV_idx) = prev_min_Rk;
Convex_Convergence_curve_AO(mc_iter,nV_idx,ao) = -prev_cost;
Convex_Fake_Convergence_curve_AO(mc_iter,nV_idx) = best_fake_secrecy;
Convex_Real_Convergence_curve_AO(mc_iter,nV_idx) = best_real_secrecy;

fprintf('\nConvex AO Finished! Best Fake Secrecy Rate = %.8f\n', best_fake_secrecy);

end
end
%end



%  %% Plot
% % Convergence Curve with Markers
% figure('Color','w'); % White background
% 
% % Define color palette 
% colors = [0, 0.4470, 0.7410;      % Blue
%           0.8500, 0.3250, 0.0980; % Orange/Red
%           0.4660, 0.6740, 0.1880; % Green
%           0.4940, 0.1840, 0.5560]; % Purple
% 
% 
% % Marker interval
% markerInterval = 50;
% 
% % Data Processing for 3D Plot
% % Averaging over mc_iter (dim 1) and selecting last ao (dim 4)
% % Resulting matrices will be [nV_idx x nH_idx]
% Z_R = squeeze(mean(Convex_min_Rk(:,:,:,end), 1));
% Z_conv = squeeze(mean(Convex_Convergence_curve_AO(:,:,:,end), 1));
% Z_fake = squeeze(mean(Convex_Fake_Convergence_curve_AO(:,:,:,end), 1));
% Z_real = squeeze(mean(Convex_Real_Convergence_curve_AO(:,:,:,end), 1));
% 
% % Define your X and Y axis scales based on your simulation parameters
% % Replace nV_vectors and nH_vectors with your actual range (e.g., 1:10)
% nV_axis = 1:size(Z_conv, 1); 
% nH_axis = 1:size(Z_conv, 2);
% [X, Y] = meshgrid(nH_axis, nV_axis);
% 
% 
% %% 3D Visualization
% figure('Color','w','Name','3D Secrecy Analysis');
% 
% % Define color map
% colormap(parula); 
% 
% % Plotting the Real Secrecy Rate
% surf(X, Y, Z_real, 'FaceAlpha', 0.8, 'EdgeColor', 'interp');
% hold on;
% 
% % Optional: Plot the Fake Secrecy Rate as a mesh or transparent surface
% % mesh(X, Y, Z_fake, 'EdgeColor', 'r', 'FaceAlpha', 0.1); 
% 
% % Aesthetics
% title('Average Secrecy Rate vs. nV and nH','FontWeight','bold','FontSize',12);
% xlabel('nH Index','FontWeight','bold','FontSize',11);
% ylabel('nV Index','FontWeight','bold','FontSize',11);
% zlabel('Secrecy Rate (b/s/Hz)','FontWeight','bold','FontSize',11);
% 
% % Enhancing the view
% grid on;
% view(45, 30); % Adjust camera angle: view(azimuth, elevation)
% colorbar;      % Shows the scale of the Z-axis
% ax = gca;
% ax.LineWidth = 1.1;
% box on;
% 
% % Add a legend if you plot multiple surfaces
% % legend('Real Secrecy Rate', 'Fake Secrecy Rate', 'Location', 'best');
% 
% figure;
% % Plotting the Real Secrecy Rate
% surf(X, Y, Z_R, 'FaceAlpha', 0.8, 'EdgeColor', 'interp');
% hold on;
% 
% % Optional: Plot the Fake Secrecy Rate as a mesh or transparent surface
% % mesh(X, Y, Z_fake, 'EdgeColor', 'r', 'FaceAlpha', 0.1); 
% 
% % Aesthetics
% title('Average Min-Rate vs. nV and nH','FontWeight','bold','FontSize',12);
% xlabel('nH Index','FontWeight','bold','FontSize',11);
% ylabel('nV Index','FontWeight','bold','FontSize',11);
% zlabel('Min Rate (b/s/Hz)','FontWeight','bold','FontSize',11);
% 
% % Enhancing the view
% grid on;
% view(45, 30); % Adjust camera angle: view(azimuth, elevation)
% colorbar;      % Shows the scale of the Z-axis
% ax = gca;
% ax.LineWidth = 1.1;
% box on;

 %% Plot
% Convergence Curve with Markers
figure('Color','w'); % White background

% Define color palette 
colors = [0, 0.4470, 0.7410;      % Blue
          0.8500, 0.3250, 0.0980; % Orange/Red
          0.4660, 0.6740, 0.1880; % Green
          0.4940, 0.1840, 0.5560]; % Purple


% Marker interval
markerInterval = 50;

Convex_min_Rk_mean = mean(Convex_min_Rk,1);
Convex_Convergence_curve_AO_mean = mean(Convex_Convergence_curve_AO,1);
Convex_Fake_Convergence_curve_AO_mean = mean(Convex_Fake_Convergence_curve_AO,1);
Convex_Real_Convergence_curve_AO_mean = mean(Convex_Real_Convergence_curve_AO,1);


hold on;
plot(Convex_Convergence_curve_AO_mean(2:end), 'Color', colors(3,:), 'LineStyle','-.', 'LineWidth',2, 'Marker','o', 'MarkerIndices',1:markerInterval:length(Convex_Convergence_curve_AO_mean), 'MarkerFaceColor',colors(3,:))


title('Convergence Curve','FontWeight','bold','FontSize',12);
xlabel('N_V','FontWeight','bold','FontSize',11);
ylabel('Best Fake Secrecy Rate','FontWeight','bold','FontSize',11);
legend('Convex-Manifold','Location','best','FontSize',10);

grid on;
ax = gca;
ax.GridAlpha = 0.3; % Lighter grid
ax.LineWidth = 1.1; % Thicker axes
box on;

% Fake & Real Secrecy Rate Curve with Markers
figure('Color','w');

% Marker definitions
markers = {'o','s','d','^','v','>','<','p','h','x'};
markerCount = 1;

% Helper function for plotting with marker cycling
plotWithMarker = @(y, color, style) plot(y, 'Color', color, 'LineStyle', style, 'LineWidth',1.5, 'Marker', markers{markerCount}, 'MarkerIndices',1:markerInterval:length(y), 'MarkerFaceColor',color);

hold on;


% Convex + Manopt
plot(Convex_min_Rk_mean(2:end), 'Color', colors(1,:), 'LineStyle','--', 'LineWidth',1.5, 'Marker','s', 'MarkerIndices',1:markerInterval:length(Convex_Fake_Convergence_curve_AO(2:end)), 'MarkerFaceColor',colors(1,:));
plot(Convex_Fake_Convergence_curve_AO_mean(2:end), 'Color', colors(3,:), 'LineStyle','--', 'LineWidth',1.5, 'Marker','s', 'MarkerIndices',1:markerInterval:length(Convex_Fake_Convergence_curve_AO(2:end)), 'MarkerFaceColor',colors(3,:));
plot(Convex_Real_Convergence_curve_AO_mean(2:end), 'Color', colors(3,:), 'LineStyle','-', 'LineWidth',1.5, 'Marker','^', 'MarkerIndices',1:markerInterval:length(Convex_Real_Convergence_curve_AO(2:end)), 'MarkerFaceColor',colors(3,:));



title('Best Fake & Real Private Secrecy Rate','FontWeight','bold','FontSize',12);
xlabel('Number of fake Eve','FontWeight','bold','FontSize',11);
ylabel('Minimum secrecy rate (b/s/Hz)','FontWeight','bold','FontSize',11);

legend('Min_rate','Convex-fake','Convex-real', ...
    'Location','best','FontSize',10);

grid on;
ax = gca;
ax.GridAlpha = 0.3;
ax.LineWidth = 1.1;
box on;

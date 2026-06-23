clear; clc;
cvx_clear;

% % %--- Choose how many workers (cores) you want ---
numWorkers =8;          % ←←← CHANGE THIS TO YOUR PREFERRED NUMBER
                          % Recommended: feature('numcores') or feature('numcores')-1

pool = gcp('nocreate');
if ~isempty(pool)
    delete(pool);   % Stop existing pool
end

parpool('local', numWorkers);  % Start new one with desired workers

Ns = 200; % number of samples for Monte Carlo simulation
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

nF = 2; % Number of fake eavesdroppers

L = 2; % number of eavesdroppers

% --- OTFS System Parameters ---
delta_f = 15e3;      % Subcarrier spacing (Hz)
T = 1/delta_f;       % Symbol duration


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


N_H = 40; % number of rows of regularly arranged unit cells of RIS
N_V = 40; % number of columns of regularly arranged unit cells of RIS

SV_vec = 0:800:8000;
% ADD THIS RIGHT BEFORE: for mc_iter = 1:Ns
N_SV = length(SV_vec);
feasible_record = zeros(Ns, N_SV);
Convex_min_Rk = zeros(Ns, N_SV);

Convex_Convergence_curve_AO = zeros(Ns, N_SV);
Convex_Fake_Convergence_curve_AO = zeros(Ns, N_SV);
Convex_Real_Convergence_curve_AO = zeros(Ns, N_SV);

feasible_record_ofdm = zeros(Ns, N_SV);
Convex_min_Rk_ofdm = zeros(Ns, N_SV);
Convex_Convergence_curve_AO_ofdm = zeros(Ns, N_SV);
Convex_Fake_Convergence_curve_AO_ofdm = zeros(Ns, N_SV);
Convex_Real_Convergence_curve_AO_ofdm = zeros(Ns, N_SV);





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

% Nakagami Parameters:

% from LEO satellite to STAR-RIS (one value for each path and assume same
% for each RIS element):
m_p = [m_rician;1*ones(P-1,1)]; % shape parameter
omega_p = (1/P)*ones(1,P); % spread parameter


    


% from STAR-RIS to users and eavesdroppers (one value for each path and 
% each receiver, and assume same for each RIS element):
m_q = ones(K+nF+L,Q_j); % shape parameter
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


parfor mc_iter = 1:Ns

    % ================================================================
% INITIALIZE TEMPORARIES (fixes uninitialized warnings)
% ================================================================
taus_u     = zeros(Pe, 1);
nus_u      = zeros(Pe, 1);
taus_u_AN  = zeros(Pe, 1);
nus_u_AN   = zeros(Pe, 1);
feasible_flag = false;
feasible_flag_ofdm = false;


for sv_idx = 1:N_SV

S_v = SV_vec(sv_idx);

omega_orb = (S_v / Rs) * orbit_normal;

% Satellite velocity (ECI)
vS = cross(omega_orb, S_xyz);
vAN = cross(omega_orb, AN_xyz);

% RIS velocity due to Earth rotation (ECI)
vR = [0;0;0];

sigma_ang = deg2rad(30);   % angular spread


g_pq = zeros(P,Q_j,K+nF+L,nSat);
Plos = zeros(K+nF+L,nSat);
PLj = zeros(K+nF+L,nSat);

h_rp = zeros(Nr, P,nSat);
h_jq = zeros(Nr, Q_j,K+nF+L,nSat);
h_e = zeros(Pe,K+nF+L,nSat);
taus_ku = zeros(Pe,K,nSat);
nus_ku = zeros(Pe,K,nSat);
taus_kq = zeros(K, P, Q_j, nSat);
nus_kq = zeros(K, P, Q_j, nSat);


% % ELEVATION (UNCHANGED)
% el = pi/2 - (0.05)*max_alpha * rand(1, K);
% 
% 
% az = (2*pi) * rand(1, K);
% 
% 
% % RADIUS (FIXED ON EARTH SURFACE)
% r = R_earth * ones(1, K);
% 
% % SPHERICAL → CARTESIAN
% [x, y, z] = sph2cart(az, el, r);
% 
% ground_users_cart = [x; y; z];   % 3 x K

%% ============================================================
% Legitimate User Locations
%
% Users are uniformly distributed inside the HAP coverage area.
%
% Geometry:
%   - Earth-centered coordinate system.
%   - HAP/STAR-RIS located above the North Pole:
%
%       R_xyz = [0;0;R_earth+HAP_altitude]
%
%   - HAP nadir point on Earth:
%
%       [0;0;R_earth]
%
% Users are generated uniformly over a circular footprint
% on the Earth's surface.
%
% This is a realistic model for HAP-assisted communications.
%% ============================================================

% Coverage radius of HAP cell [m]
coverage_radius = 20e3;      % 20 km

% ------------------------------------------------------------
% Uniform distribution over a disk
%
% sqrt(rand) is important:
% without sqrt(), users cluster near the center.
% ------------------------------------------------------------
r_user = coverage_radius * sqrt(rand(1,K));

% Random angle
theta_user = 2*pi*rand(1,K);

% ------------------------------------------------------------
% Local East-North coordinates relative to HAP nadir
%
% Since the HAP is located above the North Pole in the current
% coordinate system, the tangent plane at the nadir is the
% global x-y plane.
% ------------------------------------------------------------
x_local = r_user .* cos(theta_user);
y_local = r_user .* sin(theta_user);

% ------------------------------------------------------------
% Project users onto Earth's surface
%
% x^2 + y^2 + z^2 = R_earth^2
% ------------------------------------------------------------
z_local = sqrt(R_earth^2 - x_local.^2 - y_local.^2);

% Earth-centered coordinates
ground_users_cart = [x_local;
                     y_local;
                     z_local];

%% Optional sanity check
%
% Verify users are within desired footprint
%
% footprint_distance = sqrt(x_local.^2 + y_local.^2);
% fprintf('Max user distance from nadir = %.2f km\n', ...
%          max(footprint_distance)/1000);
% 
% d_HAP_user = vecnorm(ground_users_cart - R_xyz);
% 
% fprintf('Min HAP-user distance = %.2f km\n', min(d_HAP_user)/1000);
% fprintf('Max HAP-user distance = %.2f km\n', max(d_HAP_user)/1000);
% 
% disp(R_xyz')
% disp(ground_users_cart')

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

rho_j_xyz = [ground_users_cart,fake_eavesdroppers_xyz,eavesdroppers_xyz];

% find out whether each receiver is on the reflect side or transmit side
reflect = sign(RIS_normal.' * (rho_j_xyz - R_xyz));

M = 16;
N = 16;

delay_res = 1/(M*delta_f);
tau_rms = 0.25*delay_res;


% Satellite to RIS delays and doppler coefficients.
[taus_R, nus_R, u_paths_R] = compute_delay_and_doppler( ...
    c, S_xyz, vS, R_xyz, vR, f_c, P, tau_rms, sigma_ang);


[taus_R_AN, nus_R_AN, u_paths_R_AN] = compute_delay_and_doppler( ...
    c, AN_xyz, vAN, R_xyz, vR, f_c, P, tau_rms, sigma_ang);

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
    if norm(User_k_loc(1:2)-R_xyz(1:2)) < 1e-6
        d_ru = [R_xyz(1)+1;R_xyz(1)+1];
    else
        d_ru = User_k_loc(1:2) - R_xyz(1:2);
    end
    d_ru = d_ru / norm(d_ru);


    % receiver velocity (guaranteed norm)
    v_l  = vk_ms * [d_ru; 0];

    % RIS to legitimate users delays and doppler coefficients.
    [taus_k, nus_k, u_paths_k] = compute_delay_and_doppler( ...
    c, R_xyz, vR, User_k_loc, v_l, f_c, Q_j, tau_rms, sigma_ang);



     for p=1:P
        for q=1:Q_j
           taus_kq(k,p,q,1) = taus_R(p)+taus_k(q);
           nus_kq(k,p,q,1) = nus_R(p)+nus_k(q);

           taus_kq(k,p,q,2) = taus_R_AN(p)+taus_k(q);
           nus_kq(k,p,q,2) = nus_R_AN(p)+nus_k(q);
        end 
      end
    

    % sat to legitimate users delays and doppler coefficients.
    [taus_u, nus_u, u_paths_u] = compute_delay_and_doppler( ...
    c, S_xyz, vS, User_k_loc, v_l, f_c, Pe, tau_rms, sigma_ang);

    [taus_u_AN, nus_u_AN, u_paths_AN] = compute_delay_and_doppler( ...
    c, AN_xyz, vAN, User_k_loc, v_l, f_c, Pe, tau_rms, sigma_ang);

    g_pq(:,:,k,1) = exp(1i*2*pi*(taus_R*nus_k'));    
    g_pq(:,:,k,2) = exp(1i*2*pi*(taus_R_AN*nus_k'));    

    taus_ku(:,k,1) = taus_u; 
    nus_ku(:,k,1) = nus_u;   

    taus_ku(:,k,2) = taus_u_AN; 
    nus_ku(:,k,2) = nus_u_AN;  
   
     
    for q = 1:Q_j
        h_jq(:, q,k,1) = sqrt(gamrnd(m_q(k,q), omega_q(k,q)/m_q(k,q), [Nr, 1])) .* exp(1i*2*pi*rand(Nr,1));
        h_jq(:, q,k,2) = sqrt(gamrnd(m_q(k,q), omega_q(k,q)/m_q(k,q), [Nr, 1])) .* exp(1i*2*pi*rand(Nr,1));
    end

    
    for u = 1:Pe
        h_e(u,k,1) = sqrt(gamrnd(m_e(k,u), omega_e(k,u)/m_e(k,u))) .* exp(1i*2*pi*rand());% Direct Data carrying channel
        h_e(u,k,2) = sqrt(gamrnd(m_e(k,u), omega_e(k,u)/m_e(k,u))) .* exp(1i*2*pi*rand());% Direct Noise carrying channel
    end
            
end

% % % % Max tau and nu
% max_tau = max([taus_kq(:);taus_ku(:)])-min([taus_kq(:);taus_ku(:)]); 
% max_nu  = max([nus_kq(:);nus_ku(:)])-min([nus_kq(:);nus_ku(:)]);    



% taus_kq_rel = taus_kq - 0;
% taus_ku_rel = taus_ku - 0;
% 
% fprintf('\nSatellite S\n');
% fprintf('min taus_kq(:,:,:,1) = %.3f us\n', ...
%     min(taus_kq(:,:,:,1),[],'all')/T);
% fprintf('max taus_kq(:,:,:,1) = %.3f us\n', ...
%     max(taus_kq(:,:,:,1),[],'all')/T);
% 
% fprintf('\nSatellite AN\n');
% fprintf('min taus_kq(:,:,:,2) = %.3f us\n', ...
%     min(taus_kq(:,:,:,2),[],'all')/T);
% fprintf('max taus_kq(:,:,:,2) = %.3f us\n', ...
%     max(taus_kq(:,:,:,2),[],'all')/T);
% 
% fprintf('\nDirect S\n');
% fprintf('min taus_ku(:,:,1) = %.3f us\n', ...
%     min(taus_ku(:,:,1),[],'all')/T);
% fprintf('max taus_ku(:,:,1) = %.3f us\n', ...
%     max(taus_ku(:,:,1),[],'all')/T);
% 
% fprintf('\nDirect AN\n');
% fprintf('min taus_ku(:,:,2) = %.3f us\n', ...
%     min(taus_ku(:,:,2),[],'all')/T);
% fprintf('max taus_ku(:,:,2) = %.3f us\n', ...
%     max(taus_ku(:,:,2),[],'all')/T);


Nsymb = M*N; 

HA = zeros(Nsymb,Nsymb,P,Q_j,K+nF+L,nSat); % Relay link
HB = zeros(Nsymb,Nsymb,Pe,K+nF+L,nSat);


HA_ofdm = zeros(Nsymb,Nsymb,P,Q_j,K+nF+L,nSat); % Relay link
HB_ofdm = zeros(Nsymb,Nsymb,Pe,K+nF+L,nSat);


for k=1:K

    for p=1:P
        for q=1:Q_j
            HA(:,:,p,q,k,1) =  compute_Hp(taus_kq(k,p,q,1), nus_kq(k,p,q,1), M, N, T, delta_f, 'blocked',transmissionType);
            HA(:,:,p,q,k,2) =  compute_Hp(taus_kq(k,p,q,2), nus_kq(k,p,q,2), M, N, T, delta_f, 'blocked',transmissionType);

             HA_ofdm(:,:,p,q,k,1) =  compute_Hp(taus_kq(k,p,q,1), nus_kq(k,p,q,1), M, N, T, delta_f, 'blocked','ofdm');
             HA_ofdm(:,:,p,q,k,2) =  compute_Hp(taus_kq(k,p,q,2), nus_kq(k,p,q,2), M, N, T, delta_f, 'blocked','ofdm');
        end
    end
   

    for u = 1:Pe
        HB(:,:,u,k,1) =  compute_Hp(taus_u(u), nus_u(u), M, N, T, delta_f, 'blocked',transmissionType);
        HB(:,:,u,k,2) =  compute_Hp(taus_u_AN(u), nus_u_AN(u), M, N, T, delta_f, 'blocked',transmissionType);

        HB_ofdm(:,:,u,k,1) =  compute_Hp(taus_u(u), nus_u(u), M, N, T, delta_f, 'blocked','ofdm');
        HB_ofdm(:,:,u,k,2) =  compute_Hp(taus_u_AN(u), nus_u_AN(u), M, N, T, delta_f, 'blocked','ofdm');
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
        c, R_xyz, vR, User_l_loc, v_l, f_c, Q_j, tau_rms, sigma_ang);
    
        % sat to eavesdropper users users delays and doppler coefficients.
        [taus_u_l, nus_u_l, u_paths_u] = compute_delay_and_doppler( ...
        c, S_xyz, vS, User_l_loc, v_l, f_c, Pe, tau_rms, sigma_ang);    

           % sat to eavesdropper users users delays and doppler coefficients.
        [taus_u_l_AN, nus_u_l_AN, u_paths_u_AN] = compute_delay_and_doppler( ...
        c, AN_xyz, vAN, User_l_loc, v_l, f_c, Pe, tau_rms, sigma_ang);
    
        
        for q = 1:Q_j
            h_jq(:, q,K+l,1) = sqrt(gamrnd(m_q(K+l,q), omega_q(K+l,q)/m_q(K+l,q), [Nr, 1])) .* exp(1i*2*pi*rand(Nr,1));
            h_jq(:, q,K+l,2) = sqrt(gamrnd(m_q(K+l,q), omega_q(K+l,q)/m_q(K+l,q), [Nr, 1])) .* exp(1i*2*pi*rand(Nr,1));
        end
    
        
        for u = 1:Pe
            h_e(u,K+l,1) = sqrt(gamrnd(m_e(K+l,u), omega_e(K+l,u)/m_e(K+l,u))) .* exp(1i*2*pi*rand());
            h_e(u,K+l,2) = sqrt(gamrnd(m_e(K+l,u), omega_e(K+l,u)/m_e(K+l,u))) .* exp(1i*2*pi*rand());
        end
    

       for p=1:P
            for q=1:Q_j
                HA(:,:,p,q,K+l,1) =  compute_Hp(taus_R(p)+taus_l(q), nus_R(p)+nus_l(q), M, N, T, delta_f, 'blocked',transmissionType);
                HA(:,:,p,q,K+l,2) =  compute_Hp(taus_R_AN(p)+taus_l(q), nus_R_AN(p)+nus_l(q), M, N, T, delta_f, 'blocked',transmissionType);

                HA_ofdm(:,:,p,q,K+l,1) =  compute_Hp(taus_R(p)+taus_l(q), nus_R(p)+nus_l(q), M, N, T, delta_f, 'blocked','ofdm');
                HA_ofdm(:,:,p,q,K+l,2) =  compute_Hp(taus_R_AN(p)+taus_l(q), nus_R_AN(p)+nus_l(q), M, N, T, delta_f, 'blocked','ofdm');
            end
       end

        for u = 1:Pe
            HB(:,:,u,K+l,1) =  compute_Hp(taus_u_l(u) ,nus_u_l(u), M, N, T, delta_f, 'blocked',transmissionType);
            HB(:,:,u,K+l,2) =  compute_Hp(taus_u_l_AN(u) ,nus_u_l_AN(u), M, N, T, delta_f, 'blocked',transmissionType);

            HB_ofdm(:,:,u,K+l,1) =  compute_Hp(taus_u_l(u) ,nus_u_l(u), M, N, T, delta_f, 'blocked','ofdm');
            HB_ofdm(:,:,u,K+l,2) =  compute_Hp(taus_u_l_AN(u) ,nus_u_l_AN(u), M, N, T, delta_f, 'blocked','ofdm');
        end

        g_pq(:,:,K+l,1) = exp(1i*2*pi*(taus_R*nus_l'));  
        g_pq(:,:,K+l,2) = exp(1i*2*pi*(taus_R_AN*nus_l'));            
end



%% Sine-Cosine optimization

display('SCA is optimizing your problem');

Num_agents  = 100;
Max_iteration = 5;
Rmin=0;

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





alpha_ofdm = rand(Num_agents, K+1); 
% random values
alpha_ofdm = alpha_ofdm ./ sum(alpha_ofdm, 2);      % divide each row by its row sum

alpha_ofdm = alpha_ofdm - (sum(alpha_ofdm,2)-1)/(K+1);
alpha_ofdm = alpha_ofdm - (sum(alpha_ofdm,2)-1)/(K+1);
alpha_ofdm = alpha_ofdm - (sum(alpha_ofdm,2)-1)/(K+1);







AN_P_ratio = 0.5;  



%% ===================== CONVEX ALTERNATING OPTIMIZATION (AO) =====================
display('Convex Approximation with AO');

max_AO_iter = Max_iteration;           % Outer AO iterations
max_SCA = 3;         % Inner SCA iterations for alpha subproblem
tol = 1e-3;

Active_Gain_dB = 0; 


if any_reflect
    phi_Sr = 2*pi*rand(1,Nr);
    zeta_k_St = (10^(Active_Gain_dB/10)) *rand(1,Nr);
else
    phi_Sr = zeros(1,Nr);
    zeta_k_St = (10^(Active_Gain_dB/10)) * ones(1, Nr);
end


Convergence_curve_AO = zeros(1, max_AO_iter);

fprintf('\n=== Starting Convex AO ===\n');

manifold = complexcirclefactory(Nr,1);
problem  = struct('M', manifold);
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

[L_node_ofdm,E_node_ofdm] = compute_channels( K, Nr, nF, Pe, P, Q_j, Plos, PLj, HB_ofdm, HA_ofdm, g_pq, Nsymb, ...
reflect, h_rp, h_jq, h_e,  Active_Gain_dB);  


beta = zeros(Nr,num_agents);

min_Rsec = zeros(num_agents,1);
min_Rsec_ofdm = zeros(num_agents,1);

for i = 1:num_agents
    beta(:,i) = manifold.rand();
   [R_sec,~] = get_Secrecy_matrix(beta(:,i), L_node, E_node, alpha(1,:), K, nF, sigma2, Pw, AN_P_ratio);
   [R_sec_ofdm,~] = get_Secrecy_matrix(beta(:,i), L_node_ofdm, E_node_ofdm, alpha_ofdm(1,:), K, nF, sigma2, Pw, AN_P_ratio);

    min_Rsec(i,1) = min(min(R_sec));
    min_Rsec_ofdm(i,1) = min(min(R_sec_ofdm));
end

b0 =  beta(:,min_Rsec == max(max(min_Rsec)));
b0_ofdm =  beta(:,min_Rsec_ofdm == max(max(min_Rsec_ofdm)));


for i = 1:num_agents
    beta(:,i) = manifold.rand();
   [R_sec,~] = get_Secrecy_matrix(b0, L_node, E_node, alpha(i,:), K, nF, sigma2, Pw, AN_P_ratio);
   [R_sec_ofdm,~] = get_Secrecy_matrix(b0_ofdm , L_node_ofdm, E_node_ofdm, alpha_ofdm(i,:), K, nF, sigma2, Pw, AN_P_ratio);


    min_Rsec(i,1) = min(min(R_sec));
    min_Rsec_ofdm(i,1) = min(min(R_sec_ofdm));
end

alpha_prev = alpha(min_Rsec == max(max(min_Rsec)),:);
alpha_prev_ofdm = alpha_ofdm(min_Rsec_ofdm == max(max(min_Rsec_ofdm)),:);

alpha = alpha_prev;
alpha_ofdm = alpha_prev_ofdm;

[R_sec_prev,rate_p,~] = get_Secrecy_matrix(b0, L_node, E_node, alpha, K, nF, sigma2, Pw, AN_P_ratio);
[R_sec_prev_ofdm ,rate_ofdm ,~] = get_Secrecy_matrix(b0_ofdm , L_node_ofdm , E_node_ofdm , alpha_ofdm, K, nF, sigma2, Pw, AN_P_ratio);

phi_St = wrapToPi(angle(b0)).';
phi_St_ofdm = wrapToPi(angle(b0_ofdm)).';

 % ================================================================
    % Build X
    % ================================================================
    if any_reflect
        X = [alpha, phi_Sr, phi_St, zeta_k_St];
        X_ofdm = [alpha_ofdm, phi_Sr, phi_St_ofdm, zeta_k_St];
    else
        X = [alpha, phi_St(:).'];
        X_ofdm = [alpha_ofdm, phi_St_ofdm(:).'];
    end



[~, sc_p_lk, ~, ~, R_k,sinr_c_k, sinr_p_k, ~] = compute_sinr_sc_an(...
        Pe, P, Q_j, nF+L, K, delta_f, Plos, PLj, Nr, HB, HA, g_pq, ...
        Nsymb, reflect, Rmin, h_rp, h_jq, h_e, ...
        zeta_k_St, Active_Gain_dB,AN_P_ratio, X);


[~, sc_p_lk_ofdm, ~, ~, Rk_ofdm,sinr_c_k_ofdm, sinr_p_k_ofdm, ~] = compute_sinr_sc_an(...
        Pe, P, Q_j, nF+L, K, delta_f, Plos, PLj, Nr, HB_ofdm, HA_ofdm, g_pq, ...
        Nsymb, reflect, Rmin, h_rp, h_jq, h_e, ...
        zeta_k_St, Active_Gain_dB,AN_P_ratio, X_ofdm);

prev_cost = min(min(R_sec_prev));
prev_cost_ofdm = min(min(R_sec_prev_ofdm));

prev_min_Rk = 0;
prev_min_Rk_ofdm = 0;

best_fake_secrecy = prev_cost;
best_real_secrecy = min(min(sc_p_lk(nF+1:end,:)));


best_fake_secrecy_ofdm = prev_cost_ofdm;
best_real_secrecy_ofdm = min(min(sc_p_lk_ofdm(nF+1:end,:)));

Ck = max(0, Rmin - rate_p);
Ck_ofdm = max(0, Rmin - rate_ofdm);

feasible_ao = true;      % track feasibility of this realization

for ao = 1:max_AO_iter

    prev_fake = prev_cost;
    prev_fake_ofdm = prev_cost_ofdm;
   

    % ================================================================
    % 2. SUBPROBLEM 2: Optimize RIS Phases Φ
    % ================================================================
    [phi_St,cost_opt] = optimize_phi_manopt_fixed_alpha(...
        Rmin,L_node,E_node,problem,b0,alpha,K, nF, sigma2, Pw, AN_P_ratio,Ck);

    [phi_St_ofdm,cost_opt_ofdm] = optimize_phi_manopt_fixed_alpha(...
        Rmin,L_node_ofdm,E_node_ofdm,problem,b0_ofdm,alpha_ofdm,K, nF, sigma2, Pw, AN_P_ratio,Ck);

    b0 = exp(1i*phi_St(:));
    b0_ofdm = exp(1i*phi_St_ofdm(:));

  
    % ================================================================
    % 1. SUBPROBLEM 1: Optimize Power Allocation α
    % ================================================================
    [cost,alpha_prev,Ck,feasible_flag,xi_val] = new_optimize_alpha_cvx_fixed_phi(...
        Rmin,alpha_prev,L_node,E_node,phi_St, phi_Sr, zeta_k_St, ...
        K, nF, reflect, delta_f, Active_Gain_dB,AN_P_ratio, max_SCA);

    [cost_ofdm,alpha_prev_ofdm,Ck_ofdm,feasible_flag_ofdm,xi_val_ofdm] = new_optimize_alpha_cvx_fixed_phi(...
        Rmin,alpha_prev_ofdm,L_node_ofdm,E_node_ofdm,phi_St_ofdm, phi_Sr, zeta_k_St, ...
        K, nF, reflect, delta_f, Active_Gain_dB,AN_P_ratio, max_SCA);

    alpha = alpha_prev;
    alpha_ofdm = alpha_prev_ofdm;  
    

    % ================================================================
    % Build X
    % ================================================================
    if any_reflect
        X = [alpha, phi_Sr, phi_St, zeta_k_St];
        X_ofdm = [alpha_ofdm, phi_Sr, phi_St_ofdm, zeta_k_St];
    else
        X = [alpha, phi_St(:).'];
        X_ofdm = [alpha_ofdm, phi_St_ofdm(:).'];
    end

    % ================================================================
    % Final evaluation
    % ================================================================
    [~, sc_p_lk, ~, ~, R_k,sinr_c_k, sinr_p_k, ~] = compute_sinr_sc_an(...
        Pe, P, Q_j, nF+L, K, delta_f, Plos, PLj, Nr, HB, HA, g_pq, ...
        Nsymb, reflect, Rmin, h_rp, h_jq, h_e, ...
        zeta_k_St, Active_Gain_dB,AN_P_ratio, X);



    [~, sc_p_lk_ofdm, ~, ~, Rk_ofdm,sinr_c_k_ofdm, sinr_p_k_ofdm, ~] = compute_sinr_sc_an(...
    Pe, P, Q_j, nF+L, K, delta_f, Plos, PLj, Nr, HB_ofdm, HA_ofdm, g_pq, ...
    Nsymb, reflect, Rmin, h_rp, h_jq, h_e, ...
    zeta_k_St, Active_Gain_dB,AN_P_ratio, X_ofdm);

    rate_p_vec_ofdm = log2(1 + sinr_p_k_ofdm);

 

    % ================================================================
    % Handle feasibility properly
    % ================================================================
    if feasible_flag
        current_fake = min(min(sc_p_lk(1:nF,:)));
        current_real = min(min(sc_p_lk(nF+1:end,:)));
    else
        % Treat as outage
        R_k = zeros(K,1);
        current_fake = 0;
        current_real = 0;
        prev_min_Rk=0;        
    end

   
    if feasible_flag_ofdm        
        current_fake_ofdm = min(min(sc_p_lk_ofdm(1:nF,:)));
        current_real_ofdm = min(min(sc_p_lk_ofdm(nF+1:end,:)));
    else
        % Treat as outage
        Rk_ofdm = zeros(K,1);
        current_fake_ofdm = 0;
        current_real_ofdm = 0;
        prev_min_Rk_ofdm=0;        
     end

     
    % ================================================================
    % Update best solution
    % ================================================================
    if feasible_flag && cost > prev_cost
        best_fake_secrecy = current_fake;
        best_real_secrecy = current_real;
        Destination_position = X;
        prev_cost = cost;
        prev_min_Rk = min(R_k);
    end

 
    if feasible_flag_ofdm && cost_ofdm  > prev_cost_ofdm  
        best_fake_secrecy_ofdm  = current_fake_ofdm ;
        best_real_secrecy_ofdm = current_real_ofdm ;
        Destination_position_ofdm  = X_ofdm;
        prev_cost_ofdm = cost_ofdm ;
        prev_min_Rk_ofdm  = min(Rk_ofdm);
    end
    
  
    
    
     
 
    % ================================================================
    % Logging
    % ================================================================
      fprintf(['AO Iter %2d | Feasible = %d | Fake Sec = %.6f | Δ = %.6f | ' ...
             'max(xi)=%.2e | V_s = %2.1f | Ns=%2d\n'], ...
            ao, feasible_flag, best_fake_secrecy, ...
            best_fake_secrecy - prev_fake, max(xi_val), S_v, mc_iter);


        fprintf(['OFDM: AO Iter %2d | Feasible = %d | Fake Sec = %.6f | Δ = %.6f | ' ...
             'max(xi)=%.2e | V_s = %2.1f | Ns=%2d\n'], ...
            ao, feasible_flag_ofdm, best_fake_secrecy_ofdm, ...
            best_fake_secrecy_ofdm - prev_fake_ofdm, max(xi_val_ofdm), S_v , mc_iter);

    % if ~feasible_flag || abs(best_fake_secrecy)<1e-8 
    %     break;
    % end


    % ================================================================
    % Convergence
    % ================================================================
    % if feasible_flag && abs(best_fake_secrecy - prev_fake) < tol && ao >= 5
    %     fprintf('→ AO Converged at iteration %d\n', ao);
    %     break;
    % end

end



% fprintf('\nConvex AO Finished! Best Fake Secrecy Rate = %.8f\n', best_fake_secrecy);
 feasible_record(mc_iter,sv_idx) = feasible_flag;
 feasible_record_ofdm(mc_iter,sv_idx) = feasible_flag_ofdm;

 Convex_min_Rk(mc_iter,sv_idx) = prev_min_Rk;
 Convex_min_Rk_ofdm(mc_iter,sv_idx) = prev_min_Rk_ofdm;


 Convex_Convergence_curve_AO(mc_iter,sv_idx) = prev_cost;
 Convex_Fake_Convergence_curve_AO(mc_iter,sv_idx) = best_fake_secrecy;
 Convex_Real_Convergence_curve_AO(mc_iter,sv_idx) = best_real_secrecy;

 Convex_Convergence_curve_AO_ofdm(mc_iter,sv_idx) = prev_cost_ofdm;
 Convex_Fake_Convergence_curve_AO_ofdm(mc_iter,sv_idx) = best_fake_secrecy_ofdm;
 Convex_Real_Convergence_curve_AO_ofdm(mc_iter,sv_idx) = best_real_secrecy_ofdm;


% fprintf('alpha=[%.4f %.4f %.4f] AN_idx=%d iter=%d\n', alpha, AN_idx,mc_iter);

 
end
end



%% Plot
% Convergence Curve with Markers
figure('Color','w'); % White background

% Define color palette 
colors = [0, 0.4470, 0.7410;      % 1: Blue
          0.8500, 0.3250, 0.0980; % 2: Orange/Red
          0.4660, 0.6740, 0.1880; % 3: Green
          0.4940, 0.1840, 0.5560; % 4: Purple
          0.9290, 0.6940, 0.1250]; % 5: Yellow/Gold (Added for distinct OFDM representation)

% Marker interval
markerInterval = 2;

valid_records_Qtd = sum(feasible_record, 1); % returns a logical column vector
valid_records_Qtd_ofdm = sum(feasible_record_ofdm, 1); % returns a logical column vector

% Step 2: Compute means only for the selected rows
% RSMA (Oversight label 'Convex' in your initial code)
Convex_min_Rk_mean = sum(Convex_min_Rk.*feasible_record,1)./valid_records_Qtd;
Convex_Convergence_curve_AO_mean = sum(Convex_Convergence_curve_AO.*feasible_record,1)./valid_records_Qtd; 
Convex_Fake_Convergence_curve_AO_mean = sum(Convex_Fake_Convergence_curve_AO.*feasible_record,1)./valid_records_Qtd;  
Convex_Real_Convergence_curve_AO_mean = sum(Convex_Real_Convergence_curve_AO.*feasible_record,1)./valid_records_Qtd; 


% OFDM
Convex_min_Rk_ofdm_mean = sum(Convex_min_Rk_ofdm.*feasible_record_ofdm,1)./valid_records_Qtd_ofdm;
Convex_Convergence_curve_AO_ofdm_mean = sum(Convex_Convergence_curve_AO_ofdm.*feasible_record_ofdm,1)./valid_records_Qtd_ofdm; 
Convex_Fake_Convergence_curve_AO_ofdm_mean = sum(Convex_Fake_Convergence_curve_AO_ofdm.*feasible_record_ofdm,1)./valid_records_Qtd_ofdm;  
Convex_Real_Convergence_curve_AO_ofdm_mean = sum(Convex_Real_Convergence_curve_AO_ofdm.*feasible_record_ofdm,1)./valid_records_Qtd_ofdm;

hold on;
plot(SV_vec, Convex_Convergence_curve_AO_mean(1:end), 'Color', colors(1,:), 'LineStyle','-.', 'LineWidth',2, 'Marker','o', 'MarkerIndices',1:markerInterval:length(SV_vec), 'MarkerFaceColor',colors(1,:))
plot(SV_vec, Convex_Convergence_curve_AO_ofdm_mean(1:end), 'Color', colors(5,:), 'LineStyle','-.', 'LineWidth',2, 'Marker','s', 'MarkerIndices',1:markerInterval:length(SV_vec), 'MarkerFaceColor',colors(5,:))

xlabel('$v_s$','FontWeight','bold','FontSize',12,'Interpreter','latex');
ylabel('min. SC (b/s/Hz)','FontWeight','bold','FontSize',12,'Interpreter','latex');
legend('OTFS', 'OFDM', 'Location','best','FontSize',10);
grid on;
ax = gca;
ax.GridAlpha = 0.3; % Lighter grid
ax.LineWidth = 1.1; % Thicker axes
box on;

% Fake & Real Secrecy Rate Curve with Markers
figure('Color','w');

hold on;
% RSMA Curves
plot(SV_vec, Convex_min_Rk_mean(1:end), 'Color', colors(1,:), 'LineStyle','--', 'LineWidth',1.5, 'Marker','s', 'MarkerIndices',1:markerInterval:length(SV_vec), 'MarkerFaceColor',colors(1,:));
plot(SV_vec, Convex_Fake_Convergence_curve_AO_mean(1:end), 'Color', colors(3,:), 'LineStyle','--', 'LineWidth',1.5, 'Marker','s', 'MarkerIndices',1:markerInterval:length(SV_vec), 'MarkerFaceColor',colors(3,:));
plot(SV_vec, Convex_Real_Convergence_curve_AO_mean(1:end), 'Color', colors(3,:), 'LineStyle','-', 'LineWidth',1.5, 'Marker','^', 'MarkerIndices',1:markerInterval:length(SV_vec), 'MarkerFaceColor',colors(3,:));


% OFDM Curves
plot(SV_vec, Convex_min_Rk_ofdm_mean, 'Color', colors(5,:), 'LineStyle', '--', 'LineWidth', 1.5, 'Marker', 'o');
plot(SV_vec, Convex_Fake_Convergence_curve_AO_ofdm_mean, 'Color', colors(5,:), 'LineStyle', ':', 'LineWidth', 1.5, 'Marker', 'o');
plot(SV_vec, Convex_Real_Convergence_curve_AO_ofdm_mean, 'Color', colors(5,:), 'LineStyle', '-', 'LineWidth', 1.5, 'Marker', 'p', 'MarkerFaceColor',colors(5,:));

legend('Min Rate (OTFS)', 'Virtual SC (OTFS)', 'Real SC (OTFS)', ...
       'Min Rate (OFDM)', 'Virtual SC (OFDM)', 'Real SC (OFDM)', ...
       'Location', 'best', 'FontSize', 11);
xlabel('$v_s (m/s)$','FontWeight','bold','FontSize',12,'Interpreter','latex');
ylabel('min. SC (b/s/Hz)','FontWeight','bold','FontSize',12,'Interpreter','latex');
grid on;
ax = gca;
ax.GridAlpha = 0.3;
ax.LineWidth = 1.1;
box on;

% Store in struct
OFDMResults = struct();
% Original data
OFDMResults.Convex_min_Rk = Convex_min_Rk;
OFDMResults.Convex_min_Rk_ofdm = Convex_min_Rk_ofdm; % <--- Add
OFDMResults.Convex_Convergence_curve_AO = Convex_Convergence_curve_AO;
OFDMResults.Convex_Convergence_curve_AO_ofdm = Convex_Convergence_curve_AO_ofdm; % <--- Add
OFDMResults.Convex_Fake_Convergence_curve_AO = Convex_Fake_Convergence_curve_AO;
OFDMResults.Convex_Fake_Convergence_curve_AO_ofdm = Convex_Fake_Convergence_curve_AO_ofdm; % <--- Add
OFDMResults.Convex_Real_Convergence_curve_AO = Convex_Real_Convergence_curve_AO;
OFDMResults.Convex_Real_Convergence_curve_AO_ofdm = Convex_Real_Convergence_curve_AO_ofdm; % <--- Add
OFDMResults.feasible_record = feasible_record;
OFDMResults.feasible_record_ofdm = feasible_record_ofdm; % <--- Add

% Processed data
OFDMResults.valid_records_Qtd = valid_records_Qtd;
OFDMResults.valid_records_Qtd_ofdm = valid_records_Qtd_ofdm; % <--- Add
OFDMResults.Convex_min_Rk_mean = Convex_min_Rk_mean;
OFDMResults.Convex_min_Rk_ofdm_mean = Convex_min_Rk_ofdm_mean; % <--- Add
OFDMResults.Convex_Convergence_curve_AO_mean = Convex_Convergence_curve_AO_mean;
OFDMResults.Convex_Convergence_curve_AO_ofdm_mean = Convex_Convergence_curve_AO_ofdm_mean; % <--- Add
OFDMResults.Convex_Fake_Convergence_curve_AO_mean = Convex_Fake_Convergence_curve_AO_mean;
OFDMResults.Convex_Fake_Convergence_curve_AO_ofdm_mean = Convex_Fake_Convergence_curve_AO_ofdm_mean; % <--- Add
OFDMResults.Convex_Real_Convergence_curve_AO_mean = Convex_Real_Convergence_curve_AO_mean;
OFDMResults.Convex_Real_Convergence_curve_AO_ofdm_mean = Convex_Real_Convergence_curve_AO_ofdm_mean; % <--- Add

% Save
save('OFDMResults.mat', 'OFDMResults', '-v7.3');
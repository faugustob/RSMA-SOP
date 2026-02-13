clear; clc;
cvx_clear;

Ns = 1e6; % number of samples for Monte Carlo simulation
rng(3);

transmissionType = 'zakr';

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

nF = 5; % Number of fake eavesdroppers

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

N_V = 40; % number of rows of regularly arranged unit cells of RIS
N_H = 40; % number of columns of regularly arranged unit cells of RIS
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

AN_elevation = pi/2-max_theta*(0.009);
AN_azimuth  = pi/4 -pi;
AN_radius  = S_radius;


AN_sph = [ AN_azimuth, AN_elevation, AN_radius]; % location of AN LEO satellite in spherical coordinates

[S_x, S_y, S_z] = sph2cart(S_sph(1), S_sph(2), S_sph(3));
S_xyz = [S_x;S_y;S_z];


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

% ELEVATION (UNCHANGED)
el = pi/2 - (0.1)*max_alpha * rand(1, K);


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



% Nakagami Parameters:

% from LEO satellite to STAR-RIS (one value for each path and assume same
% for each RIS element):
m_p = [m_rician;1*ones(P-1,1)]; % shape parameter
omega_p = (1/P)*ones(1,P); % spread parameter

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
% note: if we want different paths to have different parameters, or legit
% users and eavesdroppers to have different parameters, modify the above
% matrices correspondingly, but ensure the shape doesn't change


% STAR-RIS phase and magnitude of each RIS element


% note: beta_r will have shape (Ns,Nr,K+nF+L)



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

% Satellite to RIS delays and doppler coefficients.
[taus_R, nus_R, u_paths_R] = compute_delay_and_doppler( ...
    c, S_xyz, vS, R_xyz, vR, f_c, P, sigma_ang);


g_pq = zeros(P,Q_j,K+nF+L);
Plos = zeros(K+nF+L,nSat);
PLj = zeros(K+nF+L,nSat);

h_rp = zeros(Nr, P,K+nF+L,nSat);
h_jq = zeros(Nr, Q_j,K+nF+L);
h_e = zeros(Pe,K+nF+L,nSat);
taus_ku = zeros(Pe,K);
nus_ku = zeros(Pe,K);

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
    Ns=1;
  

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

% % % Max tau and nu
max_tau = max([taus_kq(:);taus_ku(:)])-min([taus_kq(:);taus_ku(:)]); 
max_nu  = max([nus_kq(:);nus_ku(:)])-min([nus_kq(:);nus_ku(:)]);    


% Compute M and N based on the parameters
[M, N] = computeOTFSgrid(max_tau, max_nu, 'numerology', B, delta_f, T, Tf);
M = max(M, 64); N = max(N, 12);  % Minimum practical size



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

         
        Ns=1;

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

Num_agents  = 60;
Max_iteration = 10;
Rmin=0.1;

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


X = [alpha,phi_St];

if any_reflect
    dim = K+1+3*Nr;
    ub=[ones(1,K+1),2*pi*ones(1,2*Nr),ones(1,Nr)];
    alpha_min = 1e-4;
    lb = [alpha_min * ones(1,K+1),zeros(1,3*Nr)];
    zeta_k_St = (10^(Active_Gain_dB/10)) *rand(Num_agents,Nr);
    X = [alpha,phi_Sr,phi_St,zeta_k_St];
end


% --- Problem Dimensions and Bounds ---
dim_pso = dim;
alpha_min_pso = alpha_min;
lb_pso =lb;
ub_pso = ub;



Destination_position=zeros(1,dim);
Destination_fitness=inf;
best_fake_secrecy_rate=0;
best_real_secrecy_rate = 0;

Convergence_curve=zeros(1,Max_iteration);
sum_rate_curve=zeros(1,Max_iteration);

min_sum_secrecy = zeros(1,Max_iteration);

Objective_values = zeros(1,size(X,1));
%All_objective_values=zeros(Max_iteration,size(X,1));


% Calculate the fitness of the first set and find the best one
for i=1:size(X,1)
    C_k = zeros(K,1);

    [sc_c_lk,sc_p_lk,sc_p_kk,rate_c,rate_k,R_k,~] = compute_sinr_sc_an(Pe,P,Q_j,nF+L,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,Rmin,h_rp,h_jq,h_e,zeta_k_St,Active_Gain_dB,X(i,:));

    %sum_secrecy = sc_c_lk+sc_p_lk; %Private + Common secrecy capacities.

    mean_fake_p_secrecy = mean(mean(sc_p_lk(1:nF,:)));
    mean_p_secrecy = mean(mean(sc_p_lk(nF+1:end,:)));

    
    penalty = 0;
    violation = max(Rmin - R_k, 0);
    penalty = penalty + sum(violation.^2);

    
    Objective_values(1,i) = -mean_fake_p_secrecy + 1e3 * penalty;

    if i==1
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
        best_fake_secrecy_rate = mean_fake_p_secrecy;
        best_real_secrecy_rate = mean_p_secrecy;
    elseif Objective_values(1,i)<Destination_fitness
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
        best_fake_secrecy_rate = mean_fake_p_secrecy;
        best_real_secrecy_rate = mean_p_secrecy;
    end

    All_objective_values(1,i)=Objective_values(1,i);
end



%Main loop
t=2; % start from the second iteration since the first iteration was dedicated to calculating the fitness
while t<=Max_iteration

    % Eq. (3.4)
    a = 3;
    Max_iteration = Max_iteration;
    r1=a-t*((a)/Max_iteration); % r1 decreases linearly from a to 0

    % Update the position of solutions with respect to destination
    for i=1:size(X,1) % in i-th solution
       

        for j=1:size(X,2) % in j-th dimension

            % Update r2, r3, and r4 for Eq. (3.3)
            r2=(2*pi)*rand();
            r3=2*rand;
            r4=rand();

            % Eq. (3.3)
            if r4<0.5
                % Eq. (3.1)
                X(i,j)= X(i,j)+(r1*sin(r2)*abs(r3*Destination_position(j)-X(i,j)));
            else
                % Eq. (3.2)
                X(i,j)= X(i,j)+(r1*cos(r2)*abs(r3*Destination_position(j)-X(i,j)));
            end

        end
    end
    for i=1:size(X,1)
        



        % Check if solutions go outside the search spaceand bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

     
        X(i,1:K+1) = X(i,1:K+1) ./ (sum(X(i,1:K+1), 2)); % Normalize to ensure sum alpha = 1;

        X(i,1:K+1) = X(i,1:K+1) - (sum(X(i,1:K+1))-1)/(K+1);
        X(i,1:K+1) = X(i,1:K+1) - (sum(X(i,1:K+1))-1)/(K+1);

        % Calculate the objective values

        [sc_c_lk,sc_p_lk,sc_p_kk,rate_c,rate_k,R_k,~] = compute_sinr_sc_an(Pe,P,Q_j,nF+L,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,Rmin,h_rp,h_jq,h_e, zeta_k_St,Active_Gain_dB, X(i,:));
        sum_secrecy = sc_c_lk+sc_p_lk; %Private + Common secrecy capacities.
        mean_fake_p_secrecy = mean(mean(sc_p_lk(1:nF,:)));
        mean_p_secrecy = mean(mean(sc_p_lk(nF+1:end,:)));

      
       
          
        penalty = 0;
        violation = max(Rmin - R_k, 0);
        penalty = penalty + sum(violation.^2);
    
        
        Objective_values(1,i) = -mean_fake_p_secrecy + 1e3 * penalty;
     


        % Update the destination if there is a better solution
        if Objective_values(1,i)<Destination_fitness
            Destination_position=X(i,:);
            Destination_fitness=Objective_values(1,i);
            best_fake_secrecy_rate = mean_fake_p_secrecy;
            best_real_secrecy_rate = mean_p_secrecy;
        end
    end

    Convergence_curve(t)=-Destination_fitness;
    Fake_secrecy_rate_curve(t)=best_fake_secrecy_rate;
    Real_secrecy_rate_curve(t)=best_real_secrecy_rate;


    % Display the iteration and best optimum obtained so far
    if mod(t,1)==0
        display(['At iteration ', num2str(t), ' the optimum fake sc is ', num2str(best_fake_secrecy_rate), ' the optimum real sc is ', num2str(best_real_secrecy_rate)]);
    end

    % Increase the iteration counter
    t=t+1;

    

end

%% ===================== CONVEX ALTERNATING OPTIMIZATION (AO) =====================
display('Convex Approximation with AO');

max_AO_iter = 15;           % Outer AO iterations
max_SCA_inner = 10;         % Inner SCA iterations for alpha subproblem
tol = 1e-3;

% === Start from the best solution found so far (CRITICAL) ===
X = Destination_position;   


phi_St = 2*pi*rand(1,Nr);% transmission phases

Active_Gain_dB = 0; 


if any_reflect
    phi_Sr = 2*pi*rand(1,Nr);
    zeta_k_St = (10^(Active_Gain_dB/10)) *rand(Num_agents,Nr);
else
    phi_Sr = zeros(1,Nr);
    zeta_k_St = (10^(Active_Gain_dB/10)) * ones(1, Nr);
end

best_fake_secrecy = 0;
best_real_secrecy = 0;
Convergence_curve_AO = zeros(1, max_AO_iter);

fprintf('\n=== Starting Convex AO ===\n');

for ao = 1:max_AO_iter
    
    prev_fake = best_fake_secrecy;
    
    % ================================================================
    % 1. SUBPROBLEM 1: Optimize Power Allocation α  (CVX + SCA)
    % ================================================================
    alpha = optimize_alpha_cvx_fixed_phi(phi_St, phi_Sr, zeta_k_St, ...
              K, nF, L, Rmin, Pe, P, Q_j, Plos, PLj, HB, HA, g_pq, Nsymb, ...
              reflect, h_rp, h_jq, h_e, delta_f, Active_Gain_dB, max_SCA_inner);
    alpha = alpha.';

          % Rebuild X
        if any_reflect
            X = [alpha, phi_Sr, phi_St, zeta_k_St];
        else
            X = [alpha, phi_St];
        end

     [sc_c_lk,sc_p_lk,sc_p_kk,rate_c,rate_k,R_k,~] = compute_sinr_sc_an(Pe,P,Q_j,nF+L,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,Rmin,h_rp,h_jq,h_e,zeta_k_St,Active_Gain_dB,X);


    % ================================================================
    % 2. SUBPROBLEM 2: Optimize RIS Phases Φ  (Strong SCA heuristic)
    % ================================================================
    [phi_St, phi_Sr, zeta_k_St] = optimize_phi_sca_fixed_alpha(alpha, phi_St, phi_Sr, zeta_k_St, ...
              K, Nr, nF, L, Rmin, Pe, P, Q_j, Plos, PLj, HB, HA, g_pq, Nsymb, ...
              reflect, h_rp, h_jq, h_e, delta_f, Active_Gain_dB, 8);

    % Rebuild X
    if any_reflect
        X = [alpha, phi_Sr, phi_St, zeta_k_St];
    else
        X = [alpha, phi_St];
    end

    % Final evaluation
    [~, sc_p_lk, ~, ~, ~, R_k, ~] = compute_sinr_sc_an(Pe, P, Q_j, nF+L, K, delta_f, ...
        Plos, PLj, Nr, HB, HA, g_pq, Nsymb, reflect, Rmin, h_rp, h_jq, h_e, ...
        zeta_k_St, Active_Gain_dB, X);

    current_fake = mean(mean(sc_p_lk(1:nF,:)));
    current_real = mean(mean(sc_p_lk(nF+1:end,:)));

    if current_fake > best_fake_secrecy
        best_fake_secrecy = current_fake;
        best_real_secrecy = current_real;
        Destination_position = X;
    end

    Convergence_curve_AO(ao) = best_fake_secrecy;

    fprintf('AO Iter %2d | Fake Secrecy = %.4f | Real = %.4f | Δ = %.4f\n', ...
            ao, best_fake_secrecy, best_real_secrecy, best_fake_secrecy - prev_fake);

    if abs(best_fake_secrecy - prev_fake) < tol && ao >= 5
        fprintf('→ AO Converged at iteration %d\n', ao);
        break;
    end
end

fprintf('\nConvex AO Finished! Best Fake Secrecy Rate = %.4f\n', best_fake_secrecy);


%% Functions

function alpha = optimize_alpha_cvx_fixed_phi(phi_St, phi_Sr, zeta_k_St, ...
    K, nF, L, Rmin, Pe, P, Q_j, Plos, PLj, HB, HA, g_pq, Nsymb, ...
    reflect, h_rp, h_jq, h_e, delta_f, Active_Gain_dB, max_SCA)

%% ========================= CONSTANTS =========================

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

AN_P_ratio = 1;
noise = max(sigma2/Pw, 1e-10);

%% ========================= PRECOMPUTE CHANNELS =========================

Nc_k_all     = zeros(K,1);
Nc_k_AN_all  = zeros(K,1);
Nc_l_all     = zeros(nF,K);
Nc_l_AN_all  = zeros(nF,K);

for k = 1:K
    
    beta_r = beta_Sr; % adjust if reflect logic needed
    
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
        
        beta_r = beta_Sr;
        
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



gamma_j_prev = ones(K,1)*0.1;
gamma_l_prev = ones(nF,K)*0.1;
I_j_prev     = ones(K,1)*0.1;
I_l_prev     = ones(nF,K)*0.1;

gamma_c_prev = ones(K,1)*0.1;
I_c_prev     = ones(K,1)*0.1;


%% ========================= SCA LOOP =========================

for sca_iter = 1:max_SCA
    
    cvx_begin quiet
        cvx_solver mosek 
        variable vecAlpha(K+1) nonnegative
        variable gamma_j(K) nonnegative
        variable gamma_l(nF,K) nonnegative
        variable gamma_c(K) nonnegative

        variable Rc nonnegative
        variable Rp(K) nonnegative
        variable C_k(K) nonnegative



        variable s_fake(nF,K) nonnegative
        
        maximize( (1/(nF*K)) * sum(sum(s_fake)) )
        
        subject to
        
        sum(vecAlpha) <= 1;
        vecAlpha <= 1;
        sum(C_k) <= Rc;

        
        alpha_pi = vecAlpha(2:end);
        sum_alpha_pi = sum(alpha_pi);
        alpha_c = vecAlpha(1);

        %% ===== QoS =====
        for k = 1:K
            Rp(k) + C_k(k) >= Rmin;
        end

        
        for l = 1:nF
            for k = 1:K
                
                %% ===== SIGNAL & INTERFERENCE =====
                
                S_j = alpha_pi(k) * Nc_k_all(k);
                I_j = (sum_alpha_pi - alpha_pi(k)) * Nc_k_all(k) ...
                      + AN_P_ratio * Nc_k_AN_all(k);
                
                S_l = alpha_pi(k) * Nc_l_all(l,k);
                I_l = (sum_alpha_pi - alpha_pi(k)) * Nc_l_all(l,k) ...
                      + AN_P_ratio * Nc_l_AN_all(l,k);

                

                S_c = alpha_c * Nc_k_all(k);
                
                I_c = sum_alpha_pi * Nc_k_all(k) ...
                      + AN_P_ratio * Nc_k_AN_all(k);

                
                %% ===== SINR LINEARIZATION =====
                
                bilinear_j = gamma_j_prev(k)*I_j ...
                           + I_j_prev(k)*gamma_j(k) ...
                           - gamma_j_prev(k)*I_j_prev(k);
                
                bilinear_j + gamma_j(k)*noise <= S_j;
                
                
                bilinear_l = gamma_l_prev(l,k)*I_l ...
                           + I_l_prev(l,k)*gamma_l(l,k) ...
                           - gamma_l_prev(l,k)*I_l_prev(l,k);
                
                bilinear_l + gamma_l(l,k)*noise >= S_l;


                bilinear_c = gamma_c_prev(k)*I_c ...
           + I_c_prev(k)*gamma_c(k) ...
           - gamma_c_prev(k)*I_c_prev(k);

                bilinear_c + gamma_c(k)*noise <= S_c;

                
                
                %% ===== LOG LINEARIZATION =====
                
                log_l_approx = log(1 + gamma_l_prev(l,k))/log(2) ...
                    + (1/(1 + gamma_l_prev(l,k))) ...
                    * (gamma_l(l,k) - gamma_l_prev(l,k))/log(2);
                
                s_fake(l,k) <= log(1 + gamma_j(k))/log(2) - log_l_approx;
                s_fake(l,k) >= 0;               

                 
            end
            
        end

        %% Common part decodability constraint.
        for k = 1:K
            Rc <= log(1 + gamma_c(k))/log(2);
        end

        
    cvx_end
    
    if ~strcmp(cvx_status,'Solved') && ~strcmp(cvx_status,'Inaccurate/Solved')
        fprintf('SCA failed at iteration %d (%s)\n', sca_iter, cvx_status);
        break;
    end
    
    %% ===== UPDATE FOR NEXT ITERATION =====
    
    alpha_pi_val = vecAlpha(2:end);
    sum_alpha_val = sum(alpha_pi_val);
    
    for l = 1:nF
        for k = 1:K
            
            I_j_prev(k) = (sum_alpha_val - alpha_pi_val(k)) ...
                * Nc_k_all(k) + AN_P_ratio*Nc_k_AN_all(k);
            
            I_l_prev(l,k) = (sum_alpha_val - alpha_pi_val(k)) ...
                * Nc_l_all(l,k) + AN_P_ratio*Nc_l_AN_all(l,k);
        end
    end
    
    gamma_j_prev = gamma_j;
    gamma_l_prev = gamma_l;
    alpha = vecAlpha;
    
end

end




function [phi_St, phi_Sr, zeta_k_St] = optimize_phi_sca_fixed_alpha(alpha, phi_St_init, phi_Sr_init, zeta_k_St_init, ...
    K, Nr, nF, L, Rmin, Pe, P, Q_j, Plos, PLj, HB, HA, g_pq, Nsymb, ...
    reflect, h_rp, h_jq, h_e, delta_f, Active_Gain_dB, max_inner)

    % Use Particle Swarm for phase subproblem (SDP too heavy for large Nr)
    any_reflect = any(reflect > 0) && any(reflect < 0);
    
    % Build variable vector for phi: [phi_St, phi_Sr, zeta_k_St] if any_reflect
    if any_reflect
        dim = 2*Nr + Nr;  % phi_St, phi_Sr, zeta
        lb = zeros(1, dim);
        ub = [2*pi*ones(1, 2*Nr), 10^(Active_Gain_dB/10) * ones(1, Nr)];
        phi_var_init = [phi_St_init, phi_Sr_init, zeta_k_St_init];
    else
        dim = Nr;
        lb = zeros(1, dim);
        ub = 2*pi*ones(1, dim);
        phi_var_init = phi_St_init;
    end
    
    % Objective: negative mean fake secrecy + penalty for QoS
    f_phi = @(phi_var) ao_objective_phi(phi_var, alpha, K, Nr, Rmin, Pe, P, Q_j, ...
        nF, L, delta_f, Plos, PLj, HB, HA, g_pq, Nsymb, ...
        reflect, h_rp, h_jq, h_e, Active_Gain_dB, any_reflect);
    
    opts = optimoptions('particleswarm', ...
        'SwarmSize', 40, ...
        'MaxIterations', 80, ...
        'Display', 'off');
    
    phi_var_opt = particleswarm(f_phi, dim, lb, ub, opts);
    
    if any_reflect
        phi_St = phi_var_opt(1:Nr);
        phi_Sr = phi_var_opt(Nr+1:2*Nr);
        zeta_k_St = phi_var_opt(2*Nr+1:end);
    else
        phi_St = phi_var_opt;
        phi_Sr = [];
        zeta_k_St = [];
    end
end

function fitness = ao_objective_phi(phi_var, alpha, K, Nr, Rmin, Pe, P, Q_j, ...
    nF, L, delta_f, Plos, PLj, HB, HA, g_pq, Nsymb, ...
    reflect, h_rp, h_jq, h_e, Active_Gain_dB, any_reflect)
    
    if any_reflect
        phi_St = phi_var(1:Nr);
        phi_Sr = phi_var(Nr+1:2*Nr);
        zeta_k_St = phi_var(2*Nr+1:end);
        x = [alpha, phi_St, phi_Sr, zeta_k_St];
    else
        phi_St = phi_var;
        phi_Sr = zeros(1, Nr);  % Dummy
        zeta_k_St = (10^(Active_Gain_dB/10)) * ones(1, Nr);  % Fixed if not mixed
        x = [alpha, phi_St];
    end
    
    [~, sc_p_lk, ~, ~, ~, R_k, ~] = compute_sinr_sc_an_cvx( ...
        Pe, P, Q_j, nF+L, K, delta_f, Plos, PLj, Nr, HB, HA, ...
        g_pq, Nsymb, reflect, Rmin, h_rp, h_jq, h_e, ...
        zeta_k_St, Active_Gain_dB, x);
    
    mean_fake = mean(sc_p_lk(1:nF, :), 'all');
    penalty = 1e3 * sum(max(Rmin - R_k, 0).^2);
    
    fitness = -mean_fake + penalty;
end




%%
display('SCA is optimizing your problem AO');

Num_agents  = Num_agents;
Max_iteration = Max_iteration;
Rmin=0.1;

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


X = [alpha,phi_St];
% phi_init = -angle(h_jq(:, 1, 1) .* h_rp(:, 1, 1, 1)); 
% 
% % phi_init is now (Nr x 1), we transpose it to fit the agent row (1 x Nr)
% X(1, K+2:K+1+Nr) = phi_init';

if any_reflect
    dim = K+1+3*Nr;
    ub=[ones(1,K+1),2*pi*ones(1,2*Nr),ones(1,Nr)];
    alpha_min = 1e-4;
    lb = [alpha_min * ones(1,K+1),zeros(1,3*Nr)];
    zeta_k_St = (10^(Active_Gain_dB/10)) *rand(Num_agents,Nr);
    X = [alpha,phi_Sr,phi_St,zeta_k_St];
end


% --- Problem Dimensions and Bounds ---
dim_pso = dim;
alpha_min_pso = alpha_min;
lb_pso =lb;
ub_pso = ub;



Destination_position=zeros(1,dim);
Destination_fitness=inf;
best_fake_secrecy_rate=0;
best_real_secrecy_rate = 0;

Convergence_curve_AO=zeros(1,Max_iteration);
sum_rate_curve=zeros(1,Max_iteration);

min_sum_secrecy = zeros(1,Max_iteration);

Objective_values = zeros(1,size(X,1));
%All_objective_values=zeros(Max_iteration,size(X,1));


% Calculate the fitness of the first set and find the best one
for i=1:size(X,1)
    C_k = zeros(K,1);

    [sc_c_lk,sc_p_lk,sc_p_kk,rate_c,rate_k,R_k,~] = compute_sinr_sc_an(Pe,P,Q_j,nF+L,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,Rmin,h_rp,h_jq,h_e,zeta_k_St,Active_Gain_dB,X(i,:));

    %sum_secrecy = sc_c_lk+sc_p_lk; %Private + Common secrecy capacities.

    mean_fake_p_secrecy = mean(mean(sc_p_lk(1:nF,:)));
    mean_p_secrecy = mean(mean(sc_p_lk(nF+1:end,:)));

    
    penalty = 0;
    violation = max(Rmin - R_k, 0);
    penalty = penalty + sum(violation.^2);

    
    Objective_values(1,i) = -mean_fake_p_secrecy + 1e3 * penalty;

    if i==1
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
        best_fake_secrecy_rate = mean_fake_p_secrecy;
        best_real_secrecy_rate = mean_p_secrecy;
    elseif Objective_values(1,i)<Destination_fitness
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
        best_fake_secrecy_rate = mean_fake_p_secrecy;
        best_real_secrecy_rate = mean_p_secrecy;
    end

    All_objective_values(1,i)=Objective_values(1,i);
end

alpha_idx = 1:(K+1);
ris_idx = (K+2):size(X,2);

% Main loop - Alternating Optimization (alpha first)
t = 2;
while t <= Max_iteration
    % Eq. (3.4)
    a = 3;
    r1 = a - t * (a / Max_iteration);  % r1 decreases linearly

    % === Substep 1: Optimize alpha first (phases/amplitudes fixed) ===
    for i = 1:size(X, 1)
        for j = alpha_idx  % only alpha dimensions
            r2 = (2*pi) * rand();
            r3 = 2 * rand();
            r4 = rand();
            if r4 < 0.5
                X(i,j) = X(i,j) + (r1 * sin(r2) * abs(r3 * Destination_position(j) - X(i,j)));
            else
                X(i,j) = X(i,j) + (r1 * cos(r2) * abs(r3 * Destination_position(j) - X(i,j)));
            end
        end

        % Project alpha onto simplex (sum = 1, >= alpha_min)
             % Bound check
        Flag4ub = X(i,alpha_idx) > ub(alpha_idx);
        Flag4lb = X(i,alpha_idx) < lb(alpha_idx);
        
        X(i,alpha_idx) = ...
            X(i,alpha_idx).*(~(Flag4ub+Flag4lb)) + ...
            ub(alpha_idx).*Flag4ub + ...
            lb(alpha_idx).*Flag4lb;

        X(i,alpha_idx) = X(i,alpha_idx) ./ sum(X(i,alpha_idx));


        % Optional repeated correction for floating-point precision (keep your original style)
        X(i,alpha_idx) = X(i,alpha_idx) - (sum(X(i,alpha_idx))-1)/(K+1);
        X(i,alpha_idx) = X(i,alpha_idx) - (sum(X(i,alpha_idx))-1)/(K+1);
        X(i,alpha_idx) = X(i,alpha_idx) - (sum(X(i,alpha_idx))-1)/(K+1);

        % Evaluate objective with updated alpha
        [sc_c_lk, sc_p_lk, sc_p_kk, rate_c, rate_k, R_k, ~] = compute_sinr_sc_an(Pe, P, Q_j, nF+L, K, delta_f, Plos, PLj, Nr, HB, HA, g_pq, Nsymb, reflect, Rmin, h_rp, h_jq, h_e, zeta_k_St, Active_Gain_dB, X(i,:));
        mean_fake_p_secrecy = mean(mean(sc_p_lk(1:nF,:)));
        mean_p_secrecy = mean(mean(sc_p_lk(nF+1:end,:)));

        penalty = 0;
        violation = max(Rmin - R_k, 0);
        penalty = penalty + sum(violation.^2);

        Objective_values(1,i) = -mean_fake_p_secrecy + 1e3 * penalty;

        % Update global best (alpha-optimized solution)
        if Objective_values(1,i) < Destination_fitness
            Destination_position = X(i,:);
            Destination_fitness = Objective_values(1,i);
            best_fake_secrecy_rate = mean_fake_p_secrecy;
            best_real_secrecy_rate = mean_p_secrecy;
        end
    end
    alpha_fixed = X(:,alpha_idx);

    % === Substep 2: Optimize non-alpha variables (phases/amplitudes/zeta) with new alpha fixed ===
    for i = 1:size(X, 1)
        X(i,alpha_idx) = alpha_fixed(i,:);
        for j = ris_idx  % phase and zeta dimensions
            r2 = (2*pi) * rand();
            r3 = 2 * rand();
            r4 = rand();
            if r4 < 0.5
                X(i,j) = X(i,j) + (r1 * sin(r2) * abs(r3 * Destination_position(j) - X(i,j)));
            else
                X(i,j) = X(i,j) + (r1 * cos(r2) * abs(r3 * Destination_position(j) - X(i,j)));
            end
        end

        % Bound check
        Flag4ub = X(i,ris_idx) > ub(ris_idx);
        Flag4lb = X(i,ris_idx) < lb(ris_idx);
        
        X(i,ris_idx) = ...
            X(i,ris_idx).*(~(Flag4ub+Flag4lb)) + ...
            ub(ris_idx).*Flag4ub + ...
            lb(ris_idx).*Flag4lb;

       
        % Evaluate objective again
        [sc_c_lk, sc_p_lk, sc_p_kk, rate_c, rate_k, R_k, ~] = compute_sinr_sc_an(Pe, P, Q_j, nF+L, K, delta_f, Plos, PLj, Nr, HB, HA, g_pq, Nsymb, reflect, Rmin, h_rp, h_jq, h_e, zeta_k_St, Active_Gain_dB, X(i,:));
        mean_fake_p_secrecy = mean(mean(sc_p_lk(1:nF,:)));
        mean_p_secrecy = mean(mean(sc_p_lk(nF+1:end,:)));

        penalty = 0;
        violation = max(Rmin - R_k, 0);
        penalty = penalty + sum(violation.^2);

        Objective_values(1,i) = -mean_fake_p_secrecy + 1e3 * penalty;

        % Update global best
        if Objective_values(1,i) < Destination_fitness
            Destination_position = X(i,:);
            Destination_fitness = Objective_values(1,i);
            best_fake_secrecy_rate = mean_fake_p_secrecy;
            best_real_secrecy_rate = mean_p_secrecy;
        end
    end

    % Record curves
    Convergence_curve_AO(t) = -Destination_fitness;
    Fake_secrecy_rate_curve_AO(t) = best_fake_secrecy_rate;
    Real_secrecy_rate_curve_AO(t) = best_real_secrecy_rate;

    if mod(t,1) == 0
        display(['AO At iteration ', num2str(t), ' the optimum fake sc is ', num2str(best_fake_secrecy_rate), ' the optimum real sc is ', num2str(best_real_secrecy_rate)]);
    end

    t = t + 1;
end
    



%%
% phi1 = rand(1,Nr)*2*pi;
% phi2 = phi1; 
% 
% Nc1 = compute_OTFS_static_channel(0,Pe,P,Q_j,Plos(1),PLj(1),Nr,HB(:,:,:,1),HA(:,:,:,:,1),g_pq(:,:,1),phi1,Nsymb,h_rp(:,:,1),h_jq(:,:,1),h_e(:,1),'loop')
% Nc2 = compute_OTFS_static_channel(0,Pe,P,Q_j,Plos(1),PLj(1),Nr,HB(:,:,:,1),HA(:,:,:,:,1),g_pq(:,:,1),phi2,Nsymb,h_rp(:,:,1),h_jq(:,:,1),h_e(:,1),'vectorized')


%%

%% Hybrid SCA optimization (Brajevic et al)
display('Hybrid SCA is optimizing your problem');

SP  = Num_agents;
MNI = Max_iteration;

any_reflect = any(reflect > 0) && any(reflect < 0);

Active_Gain_dB = 0;
zeta_k_St = (10^(Active_Gain_dB/10)) * ones(1, Nr);

phi_Sr = 2*pi*rand(SP,Nr);
phi_St = 2*pi*rand(SP,Nr);

alpha_hsca = rand(SP, K+1);
alpha_hsca = alpha_hsca ./ sum(alpha_hsca,2);
for k = 1:3
    alpha_hsca = alpha_hsca - (sum(alpha_hsca,2)-1)/(K+1);
end

X = [alpha_hsca , phi_St];
dim_hsca = size(X,2);

if any_reflect
    zeta_k_St = (10^(Active_Gain_dB/10)) * rand(SP,Nr);
    X = [alpha_hsca , phi_Sr , phi_St , zeta_k_St];
    dim_hsca = size(X,2);
end

lb_hsca = lb;
ub_hsca = ub;
alpha_min_hsca = alpha_min;

y_optimal = zeros(1,dim_hsca);

objective = -inf(SP,1);
population_fake_secrecy_rate = zeros(SP,1);
population_real_secrecy_rate = zeros(SP,1);

best_objective = -inf;
best_mean_fake_secrecy_rate = 0;
best_mean_real_secrecy_rate = 0;

t = 1;
a = 0.75;
r1 = a;
MR = 0.1;
MR_max = 0.9;

while t <= MNI

    %% ================= Evaluate Population =================
    for i = 1:SP
        [~,sc_p_lk,~,~,~,R_k,~] = compute_sinr_sc_an( ...
            Pe,P,Q_j,nF+L,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq, ...
            Nsymb,reflect,Rmin,h_rp,h_jq,h_e,zeta_k_St, ...
            Active_Gain_dB,X(i,:));

        mean_fake = mean(mean(sc_p_lk(1:nF,:)));
        mean_real = mean(mean(sc_p_lk(nF+1:end,:)));

        violation = max(Rmin - R_k,0);
        penalty = -sum(violation.^2);

        objective(i) = mean_fake + 1e3*penalty;
        population_fake_secrecy_rate(i) = mean_fake;
        population_real_secrecy_rate(i) = mean_real;
    end

    %% ================= Global Best (ONLY PLACE IT UPDATES) =================
    [iter_best_obj, idx] = max(objective);

    if iter_best_obj > best_objective
        best_objective = iter_best_obj;
        y_optimal = X(idx,:);
        best_mean_fake_secrecy_rate = population_fake_secrecy_rate(idx);
        best_mean_real_secrecy_rate = population_real_secrecy_rate(idx);
    end

    display(['HSCA Iter ',num2str(t), ...
             ' | Fake = ',num2str(best_mean_fake_secrecy_rate), ...
             ' | Real = ',num2str(best_mean_real_secrecy_rate)]);

    %% ================= Population Update =================
    if mod(t,2) == 0
        % -------- Modified SCA --------
        for i = 1:SP
            neigh = randsample([1:i-1,i+1:SP],2);
            r2 = 2*pi*rand(1,dim_hsca);
            R_ij = rand(1,dim_hsca);
            rand_i = rand;

            under = R_ij < 0.5;
            over  = ~under;

            v_i = X(neigh(1),:) ...
                + rand_i .* abs(y_optimal - X(neigh(2),:)) ...
                + r1 .* abs(y_optimal - X(i,:)) .* ...
                  (under.*sin(r2) + over.*cos(r2));

            v_i = (v_i<lb_hsca).*(2*lb_hsca-v_i)+(v_i>=lb_hsca).*v_i;
            v_i = (v_i>ub_hsca).*(2*ub_hsca-v_i)+(v_i<=ub_hsca).*v_i;

            alpha = max(v_i(1:K+1),alpha_min_hsca);
            alpha = alpha/sum(alpha);
            for k=1:3
                alpha = alpha - (sum(alpha)-1)/(K+1);
            end
            v_i(1:K+1) = alpha;

            if sum(alpha) > 1, continue; end

            [~,sc_p_lk,~,~,~,R_k,~] = compute_sinr_sc_an( ...
                Pe,P,Q_j,nF+L,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq, ...
                Nsymb,reflect,Rmin,h_rp,h_jq,h_e,zeta_k_St, ...
                Active_Gain_dB,v_i);

            mean_fake = mean(mean(sc_p_lk(1:nF,:)));
            violation = max(Rmin - R_k,0);
            penalty = -sum(violation.^2);
            new_obj = mean_fake + 1e3*penalty;

            if new_obj > objective(i)
                objective(i) = new_obj;
                X(i,:) = v_i;
                population_fake_secrecy_rate(i) = mean_fake;
                population_real_secrecy_rate(i) = mean(mean(sc_p_lk(nF+1:end,:)));
            end
        end
    else
        % -------- ABC Strategies --------
        for i = 1:SP
            neigh = randsample([1:i-1,i+1:SP],2);
            phi_i = 2*rand-1;

            if rand < 0.5
                v_i = X(i,:) + (rand(1,dim_hsca)<MR).*phi_i.*(X(i,:)-X(neigh(1),:));
            else
                v_i = X(i,:) + phi_i*(X(neigh(1),:)-X(neigh(2),:));
            end

            v_i = (v_i<lb_hsca).*(2*lb_hsca-v_i)+(v_i>=lb_hsca).*v_i;
            v_i = (v_i>ub_hsca).*(2*ub_hsca-v_i)+(v_i<=ub_hsca).*v_i;

            alpha = max(v_i(1:K+1),alpha_min_hsca);
            alpha = alpha/sum(alpha);
            for k=1:3
                alpha = alpha - (sum(alpha)-1)/(K+1);
            end
            v_i(1:K+1) = alpha;

            if sum(alpha) > 1, continue; end

            [~,sc_p_lk,~,~,~,R_k,~] = compute_sinr_sc_an( ...
                Pe,P,Q_j,nF+L,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq, ...
                Nsymb,reflect,Rmin,h_rp,h_jq,h_e,zeta_k_St, ...
                Active_Gain_dB,v_i);

            mean_fake = mean(mean(sc_p_lk(1:nF,:)));
            violation = max(Rmin - R_k,0);
            penalty = -sum(violation.^2);
            new_obj = mean_fake + 1e3*penalty;

            if new_obj > objective(i)
                objective(i) = new_obj;
                X(i,:) = v_i;
                population_fake_secrecy_rate(i) = mean_fake;
                population_real_secrecy_rate(i) = mean(mean(sc_p_lk(nF+1:end,:)));
            end
        end
    end

    %% ================= Logging =================
    HSCA_Convergence_curve(t)        = best_objective;
    HSCA_Fake_secrecy_rate_curve(t)  = best_mean_fake_secrecy_rate;
    HSCA_Real_secrecy_rate_curve(t)  = best_mean_real_secrecy_rate;

    r1 = a*(1 - t/MNI);
    MR = min(MR_max, MR + (MR_max-0.1)/(max(rand,1e-3)*MNI));
    t = t + 1;
end

best_HSCA = best_mean_fake_secrecy_rate;



%% Hybrid SCA optimization (Brajevic et al) AO
display('Hybrid AO SCA is optimizing your problem AO');

SP = Num_agents;
MNI = Max_iteration;

% Check if more than one STAR-RIS side is being used.
any_reflect = any(reflect > 0) && any(reflect < 0);

% Problem bounds and dimensionality

Active_Gain_dB = 0; 
zeta_k_St = (10^(Active_Gain_dB/10)) * ones(1, Nr);


phi_Sr = 2*pi*rand(SP,Nr);
phi_St = 2*pi*rand(SP,Nr);% transmission phases


alpha_hsca = rand(SP, K+1); 
% random values
alpha_hsca  = alpha_hsca  ./ sum(alpha_hsca , 2);      % divide each row by its row sum

alpha_hsca  = alpha_hsca  - (sum(alpha_hsca ,2)-1)/(K+1);
alpha_hsca  = alpha_hsca  - (sum(alpha_hsca ,2)-1)/(K+1);
alpha_hsca  = alpha_hsca  - (sum(alpha_hsca ,2)-1)/(K+1);

dim_hsca = dim;
alpha_min_hsca = alpha_min;
lb_hsca =lb;
ub_hsca = ub;

X = [alpha_hsca ,phi_St];

if any_reflect  
   
    zeta_k_St = (10^(Active_Gain_dB/10)) *rand(SP,Nr);
    X = [alpha_hsca,phi_Sr,phi_St,zeta_k_St];
    dim_hsca = size(X,2);
end

y_optimal = zeros(1, dim_hsca); % best solution reached so far

best_mean_fake_secrecy_rate=-Inf;
best_mean_real_secrecy_rate=-Inf;
best_objective = -10;

objective = zeros(SP,1); % sum rate of each agent in population
population_real_secrecy_rate = zeros(SP,1); % sum rate of each agent in population
population_fake_secrecy_rate = zeros(SP,1);

t = 1; % iteration counter
a = 0.75; % constant
r1 = a; % parameter that regulates the location of the novel solution
MR = 0.1; % modification rate control parameter
MR_max = 0.9;
P_control = 0.3;

alpha_idx = 1:(K+1);
ris_idx = (K+2):size(X,2);

t = 1;
while t <= MNI
        display([' AO HSCA At iteration ', num2str(t), ' the optimum fake sc is ', num2str(best_mean_fake_secrecy_rate), ' the optimum real sc is ', num2str(best_mean_real_secrecy_rate)]);

    if mod(t, 2) == 0
        % Even: Modified SCA branch with AO
        for i = 1:SP
            neighbours = randsample([1:i-1 i+1:SP], 2);  % Compute once per agent
            % Substep 1: alpha only
            r2 = 2*pi*rand(1, length(alpha_idx));
            rand_i = rand();
            R_ij = rand(1, length(alpha_idx));
            under = R_ij < 0.5; over = 1-under;
            v_i_alpha = X(i,alpha_idx) + rand_i.*abs(y_optimal(alpha_idx)-X(neighbours(2),alpha_idx)) + r1.*abs(y_optimal(alpha_idx)-X(i,alpha_idx)).*(under.*sin(r2)+over.*cos(r2));
            v_i_alpha = (v_i_alpha<lb_hsca(alpha_idx)).*(2*lb_hsca(alpha_idx)-v_i_alpha) + (v_i_alpha>=lb_hsca(alpha_idx)).*v_i_alpha;
            v_i_alpha = (v_i_alpha>ub_hsca(alpha_idx)).*(2*ub_hsca(alpha_idx)-v_i_alpha) + (v_i_alpha<=ub_hsca(alpha_idx)).*v_i_alpha;
            alpha = v_i_alpha;
           
            alpha = alpha / sum(alpha);
            alpha = alpha - (sum(alpha)-1)/(K+1);
            alpha = alpha - (sum(alpha)-1)/(K+1);
            alpha = alpha - (sum(alpha)-1)/(K+1);
            alpha = max(alpha, alpha_min_hsca);
             alpha = max(alpha, alpha_min_hsca);
            alpha = alpha / sum(alpha);
            v_i = X(i,:);
            v_i(alpha_idx) = alpha;

            % Evaluate after alpha substep
            [sc_c_lk,sc_p_lk,sc_p_kk,rate_c,rate_k,R_k,~] = compute_sinr_sc_an(Pe,P,Q_j,nF+L,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,Rmin,h_rp,h_jq,h_e, zeta_k_St, Active_Gain_dB, v_i);
            mean_fake = mean(mean(sc_p_lk(1:nF,:)));
            penalty = sum(max(Rmin - R_k, 0).^2);
            new_obj = mean_fake - 1e3 * penalty;
            if new_obj > objective(i)
                objective(i) = new_obj;
                population_fake_secrecy_rate(i) = mean_fake;
                population_real_secrecy_rate(i) = mean(mean(sc_p_lk(nF+1:end,:)));
                X(i,alpha_idx) = alpha;
                % Immediate global best check
                if new_obj > best_objective
                    best_objective = new_obj;
                    y_optimal = X(i,:);
                    best_mean_fake_secrecy_rate = population_fake_secrecy_rate(i);
                    best_mean_real_secrecy_rate = population_real_secrecy_rate(i);
                end
            end

            % Substep 2: RIS params only (alpha fixed)
            r2 = 2*pi*rand(1, length(ris_idx));
            rand_i = rand();
            R_ij = rand(1, length(ris_idx));
            under = R_ij < 0.5; over = 1-under;
            v_i_ris = X(i,ris_idx) + rand_i.*abs(y_optimal(ris_idx)-X(neighbours(2),ris_idx)) + r1.*abs(y_optimal(ris_idx)-X(i,ris_idx)).*(under.*sin(r2)+over.*cos(r2));
            v_i_ris = (v_i_ris<lb_hsca(ris_idx)).*(2*lb_hsca(ris_idx)-v_i_ris) + (v_i_ris>=lb_hsca(ris_idx)).*v_i_ris;
            v_i_ris = (v_i_ris>ub_hsca(ris_idx)).*(2*ub_hsca(ris_idx)-v_i_ris) + (v_i_ris<=ub_hsca(ris_idx)).*v_i_ris;
            v_i = X(i,:);
            v_i(ris_idx) = v_i_ris;

            % Evaluate after RIS substep
            [sc_c_lk,sc_p_lk,sc_p_kk,rate_c,rate_k,R_k,~] = compute_sinr_sc_an(Pe,P,Q_j,nF+L,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,Rmin,h_rp,h_jq,h_e, zeta_k_St, Active_Gain_dB, v_i);
            mean_fake = mean(mean(sc_p_lk(1:nF,:)));
            penalty = sum(max(Rmin - R_k, 0).^2);
            new_obj = mean_fake - 1e3 * penalty;
            if new_obj > objective(i)
                objective(i) = new_obj;
                population_fake_secrecy_rate(i) = mean_fake;
                population_real_secrecy_rate(i) = mean(mean(sc_p_lk(nF+1:end,:)));
                X(i,ris_idx) = v_i_ris;
                if new_obj > best_objective
                    best_objective = new_obj;
                    y_optimal = X(i,:);
                    best_mean_fake_secrecy_rate = population_fake_secrecy_rate(i);
                    best_mean_real_secrecy_rate = population_real_secrecy_rate(i);
                end
            end
        end

else
        % Odd: Multi-strategy ABC branch with AO
        for i = 1:SP
            % Strategy Selection: S_i = 1 (Exploration), S_i = 2 (Exploitation)
            S_i = (rand > P_control) + 1; 
            neighbours = randsample([1:i-1 i+1:SP], 2);
            
            % --- SUBSTEP 1: Alpha (Power Allocation) Only ---
            v_i_alpha = X(i, alpha_idx);
            phi_rand = (2*rand(1, length(alpha_idx)) - 1); % ABC step size
            
            if S_i == 1
                % Strategy 1: Search based on random neighbor
                v_i_alpha = v_i_alpha + phi_rand .* (v_i_alpha - X(neighbours(1), alpha_idx));
            else
                % Strategy 2: Search based on global best (Exploitation)
                v_i_alpha = v_i_alpha + phi_rand .* (y_optimal(alpha_idx) - v_i_alpha);
            end
            
            % Projection to Simplex & Minimum Constraint
            alpha = max(v_i_alpha, alpha_min_hsca);                  
            alpha = alpha / sum(alpha);
            alpha = alpha - (sum(alpha)-1)/(K+1);
            alpha = alpha - (sum(alpha)-1)/(K+1);
            alpha = alpha - (sum(alpha)-1)/(K+1);
            alpha = max(alpha, alpha_min_hsca);
            alpha = alpha / sum(alpha);
            % Prepare agent for evaluation
            v_i = X(i,:);
            v_i(alpha_idx) = alpha;
            
            
            [sc_c_lk, sc_p_lk, ~, ~, ~, R_k, ~] = compute_sinr_sc_an(Pe,P,Q_j,nF+L,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,Rmin,h_rp,h_jq,h_e, zeta_k_St, Active_Gain_dB, v_i);
            
            mean_fake = mean(mean(sc_p_lk(1:nF,:)));
            penalty = sum(max(Rmin - R_k, 0).^2);
            new_obj = mean_fake - 1e3 * penalty;

            % Greedy Selection for Substep 1
            if new_obj > objective(i)
                objective(i) = new_obj;
                X(i, alpha_idx) = alpha;
                population_fake_secrecy_rate(i) = mean_fake;
                population_real_secrecy_rate(i) = mean(mean(sc_p_lk(nF+1:end,:)));
                if new_obj > best_objective
                    best_objective = new_obj;
                    y_optimal = X(i,:);
                     best_mean_fake_secrecy_rate = population_fake_secrecy_rate(i);
                    best_mean_real_secrecy_rate = population_real_secrecy_rate(i);
                end
            end

            % --- SUBSTEP 2: RIS Parameters Only (Alpha Fixed) ---
            v_i_ris = X(i, ris_idx);
            phi_rand = (2*rand(1, length(ris_idx)) - 1);
            
            if S_i == 1
                v_i_ris = v_i_ris + phi_rand .* (v_i_ris - X(neighbours(2), ris_idx));
            else
                v_i_ris = v_i_ris + phi_rand .* (y_optimal(ris_idx) - v_i_ris);
            end
            
            % Boundary handling (Reflection)
            v_i_ris = (v_i_ris < lb_hsca(ris_idx)).*(2*lb_hsca(ris_idx) - v_i_ris) + (v_i_ris >= lb_hsca(ris_idx)).*v_i_ris;
            v_i_ris = (v_i_ris > ub_hsca(ris_idx)).*(2*ub_hsca(ris_idx) - v_i_ris) + (v_i_ris <= ub_hsca(ris_idx)).*v_i_ris;

            v_i = X(i,:);
            v_i(ris_idx) = v_i_ris;
            
     
            [sc_c_lk, sc_p_lk, ~, ~, ~, R_k, ~] = compute_sinr_sc_an(Pe,P,Q_j,nF+L,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,Rmin,h_rp,h_jq,h_e, zeta_k_St, Active_Gain_dB, v_i);
            
            mean_fake = mean(mean(sc_p_lk(1:nF,:)));
            penalty = sum(max(Rmin - R_k, 0).^2);
            new_obj = mean_fake - 1e3 * penalty;

            % Greedy Selection for Substep 2
            if new_obj > objective(i)
                objective(i) = new_obj;
                X(i, ris_idx) = v_i_ris;
                population_fake_secrecy_rate(i) = mean_fake;
                population_real_secrecy_rate(i) = mean(mean(sc_p_lk(nF+1:end,:)));
                if new_obj > best_objective
                    best_objective = new_obj;
                    y_optimal = X(i,:);
                    best_mean_fake_secrecy_rate = population_fake_secrecy_rate(i);
                    best_mean_real_secrecy_rate = population_real_secrecy_rate(i);
                end
            end
        end

    end

    HSCA_Convergence_curve_AO(t) = best_objective;
    HSCA_Fake_secrecy_rate_curve_AO(t) = best_mean_fake_secrecy_rate;
    HSCA_Real_secrecy_rate_curve_AO(t) = best_mean_real_secrecy_rate;

    r1 = a - t*a/MNI;
    P = rand();
    if MR < MR_max
        MR = MR + (MR_max-0.1)/(P*MNI);
    else
        MR = MR_max;
    end
    t = t + 1;

end
best_HSCA = best_mean_fake_secrecy_rate;



%% Particle Swarm Optimization

display('Using MATLAB built-in particleswarm for optimization...');


PSO_Convergence_curve = nan(1, Max_iteration);
PSO_Fake_secrecy_rate_curve = nan(1, Max_iteration);
PSO_Real_secrecy_rate_curve = nan(1, Max_iteration);

params.K = K;
params.Nr = Nr;
params.Rmin = Rmin;
params.Pe = Pe;
params.P = P;
params.Q_j = Q_j;
params.nF = nF;
params.L = L;
params.delta_f = delta_f;
params.Plos = Plos;
params.PLj = PLj;
params.HB = HB;
params.HA = HA;
params.g_pq = g_pq;
params.Nsymb = Nsymb;
params.reflect = reflect;
params.h_rp = h_rp;
params.h_jq = h_jq;
params.h_e = h_e;
params.Active_Gain_dB = Active_Gain_dB;


% --- PSO Options ---
% You can adjust SwarmSize and MaxIterations to match your original SCA settings
options = optimoptions('particleswarm', ...
    'SwarmSize', 50, ...
    'MaxIterations', Max_iteration, ...
    'Display', 'iter', ...
     'OutputFcn', @psoOutputFcn, ...
    'PlotFcn', @(optimValues,state) myCustomPlot(optimValues,state));

% --- Define the Objective Function Wrapper ---
% We pass all your environment variables into the function handle
fitness_func = @(x) objective_wrapper(x, params);

% --- Run Built-in PSO ---
[best_x, best_fval] = particleswarm(fitness_func, dim_pso, lb_pso, ub_pso, options);

% --- Post-Processing ---
% Extract final best sum rate from the best position found
[~, final_sec_rate,mean_p_secrecy] = objective_wrapper(best_x, params);

display(['Optimization complete. Best min. fake private sec. Rate: ', num2str(final_sec_rate),'Best min. real private sec. Rate: ', num2str(mean_p_secrecy)]);

% --- The Objective Function Wrapper (Local Function) ---
function [fitness, min_fake, min_real] = objective_wrapper(x, params)
    K = params.K;
    Nr = params.Nr;
    Rmin = params.Rmin;
    Pe = params.Pe;
    P = params.P;
    Q_j = params.Q_j;
    nF = params.nF;
    L = params.L;
    delta_f = params.delta_f;
    Plos = params.Plos;
    PLj = params.PLj;
    HB = params.HB;
    HA = params.HA;
    g_pq = params.g_pq;
    Nsymb = params.Nsymb;
    reflect = params.reflect;
    h_rp = params.h_rp;
    h_jq = params.h_jq;
    h_e = params.h_e;
    Active_Gain_dB = params.Active_Gain_dB;

    % Extract alpha, phi, zeta, etc. from x
    alpha = x(1:K+1);
    alpha = alpha ./ sum(alpha);
    phi_St = x(K+2:K+1+Nr);
    zeta_k_St = 10^(Active_Gain_dB/10) * ones(1,Nr);

    any_reflect = any(reflect>0) && any(reflect<0);
    X = [alpha, phi_St];
    if any_reflect
        phi_Sr = x(K+2+Nr:K+1+2*Nr);
        zeta_k_St = 10^(Active_Gain_dB/10) * x(K+2+2*Nr:K+1+3*Nr);
        X = [alpha, phi_St, phi_Sr, zeta_k_St];
    end

    [sc_c_lk, sc_p_lk, sc_p_kk, rate_c, rate_k, R_k, ~] = ...
        compute_sinr_sc_an(Pe,P,Q_j,nF+L,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,Rmin,h_rp,h_jq,h_e,zeta_k_St,Active_Gain_dB,X);

    min_fake = mean(mean(sc_p_lk(1:nF,:)));
    min_real = mean(mean(sc_p_lk(nF+1:end,:)));

    penalty = 1e3*sum(max(Rmin-R_k,0).^2);
    fitness = -min_fake + penalty;
end


function stop = myCustomPlot(optimValues, state)
    stop = false;
    % We plot the negative of the best fval to see the actual Sum Rate
    plot(optimValues.iteration, -optimValues.bestfval, 'bo', 'MarkerFaceColor', 'b');
    xlabel('Iteration');
    ylabel('Actual min priv. secrecy');
    grid on;
    title('Maximization Progress');
    hold on;
end


function stop = psoOutputFcn(optimValues, state)
    stop = false;

    persistent iter PSO_Convergence_curve PSO_Fake_secrecy_rate_curve PSO_Real_secrecy_rate_curve params_local

    if strcmp(state, 'init')
        iter = 1;
        % Initialize the curves
        PSO_Convergence_curve = nan(1, 100); % or Max_iteration
        PSO_Fake_secrecy_rate_curve = nan(1, 100);
        PSO_Real_secrecy_rate_curve = nan(1, 100);
        % Copy params from base workspace
        if evalin('base', 'exist(''params'',''var'')')
            params_local = evalin('base','params');
        else
            error('params struct not found in base workspace');
        end

    elseif strcmp(state, 'iter')
        PSO_Convergence_curve(iter) = -optimValues.bestfval;
        best_x = optimValues.bestx;

        % Call objective_wrapper with the params struct
        [~, fake_sec, real_sec] = objective_wrapper(best_x, params_local);

        PSO_Fake_secrecy_rate_curve(iter) = fake_sec;
        PSO_Real_secrecy_rate_curve(iter) = real_sec;

        iter = iter + 1;

    elseif strcmp(state, 'done')
        % You can save the curves to base workspace
        assignin('base','PSO_Convergence_curve', PSO_Convergence_curve);
        assignin('base','PSO_Fake_secrecy_rate_curve', PSO_Fake_secrecy_rate_curve);
        assignin('base','PSO_Real_secrecy_rate_curve', PSO_Real_secrecy_rate_curve);
    end
end




%% Plot
%% Convergence Curve with Markers
figure('Color','w'); % White background

% Define color palette (colorblind-friendly)
colors = lines(5);

% Marker interval
markerInterval = 50;

% Plot each curve with markers every 5 points
plot(Convergence_curve(2:end), 'Color', colors(1,:), 'LineStyle','-', 'LineWidth',1.8, 'Marker','o', 'MarkerIndices',1:markerInterval:length(Convergence_curve(2:end)), 'MarkerFaceColor',colors(1,:))
hold on;
plot(Convergence_curve_AO(2:end), 'Color', colors(1,:), 'LineStyle','--', 'LineWidth',1.8, 'Marker','s', 'MarkerIndices',1:markerInterval:length(Convergence_curve_AO(2:end)), 'MarkerFaceColor',colors(1,:))
plot(HSCA_Convergence_curve(2:end), 'Color', colors(2,:), 'LineStyle','-', 'LineWidth',1.8, 'Marker','d', 'MarkerIndices',1:markerInterval:length(HSCA_Convergence_curve(2:end)), 'MarkerFaceColor',colors(2,:))
plot(HSCA_Convergence_curve_AO(2:end), 'Color', colors(2,:), 'LineStyle','--', 'LineWidth',1.8, 'Marker','^', 'MarkerIndices',1:markerInterval:length(HSCA_Convergence_curve_AO(2:end)), 'MarkerFaceColor',colors(2,:))
plot(PSO_Convergence_curve(1:end), 'Color', colors(5,:), 'LineStyle','-.', 'LineWidth',2, 'Marker','v', 'MarkerIndices',1:markerInterval:length(PSO_Convergence_curve), 'MarkerFaceColor',colors(5,:))

title('Convergence Curve','FontWeight','bold','FontSize',12);
xlabel('Iteration','FontWeight','bold','FontSize',11);
ylabel('Best Fake Secrecy Rate','FontWeight','bold','FontSize',11);
legend('SCA','SCA-AO','HSCA','HSCA-AO','PSO','Location','best','FontSize',10);
grid on;
ax = gca;
ax.GridAlpha = 0.3; % Lighter grid
ax.LineWidth = 1.1; % Thicker axes
box on;

%% Fake & Real Secrecy Rate Curve with Markers
figure('Color','w');

% Marker definitions
markers = {'o','s','d','^','v','>','<','p','h','x'};
markerCount = 1;

% Helper function for plotting with marker cycling
plotWithMarker = @(y, color, style) plot(y, 'Color', color, 'LineStyle', style, 'LineWidth',1.5, 'Marker', markers{markerCount}, 'MarkerIndices',1:markerInterval:length(y), 'MarkerFaceColor',color);

hold on;

% SCA & SCA-AO
plot(Fake_secrecy_rate_curve(2:end), 'Color', colors(1,:), 'LineStyle','-', 'LineWidth',1.5, 'Marker','o', 'MarkerIndices',1:markerInterval:length(Fake_secrecy_rate_curve(2:end)), 'MarkerFaceColor',colors(1,:));
plot(Fake_secrecy_rate_curve_AO(2:end), 'Color', colors(1,:), 'LineStyle','--', 'LineWidth',1.5, 'Marker','s', 'MarkerIndices',1:markerInterval:length(Fake_secrecy_rate_curve_AO(2:end)), 'MarkerFaceColor',colors(1,:));
plot(Real_secrecy_rate_curve(2:end), 'Color', colors(1,:), 'LineStyle','-.', 'LineWidth',1.5, 'Marker','d', 'MarkerIndices',1:markerInterval:length(Real_secrecy_rate_curve(2:end)), 'MarkerFaceColor',colors(1,:));
plot(Real_secrecy_rate_curve_AO(2:end), 'Color', colors(1,:), 'LineStyle',':', 'LineWidth',1.5, 'Marker','^', 'MarkerIndices',1:markerInterval:length(Real_secrecy_rate_curve_AO(2:end)), 'MarkerFaceColor',colors(1,:));

% HSCA & HSCA-AO
plot(HSCA_Fake_secrecy_rate_curve(2:end), 'Color', colors(2,:), 'LineStyle','-', 'LineWidth',1.5, 'Marker','o', 'MarkerIndices',1:markerInterval:length(HSCA_Fake_secrecy_rate_curve(2:end)), 'MarkerFaceColor',colors(2,:));
plot(HSCA_Fake_secrecy_rate_curve_AO(2:end), 'Color', colors(2,:), 'LineStyle','--', 'LineWidth',1.5, 'Marker','s', 'MarkerIndices',1:markerInterval:length(HSCA_Fake_secrecy_rate_curve_AO(2:end)), 'MarkerFaceColor',colors(2,:));
plot(HSCA_Real_secrecy_rate_curve(2:end), 'Color', colors(2,:), 'LineStyle','-.', 'LineWidth',1.5, 'Marker','d', 'MarkerIndices',1:markerInterval:length(HSCA_Real_secrecy_rate_curve(2:end)), 'MarkerFaceColor',colors(2,:));
plot(HSCA_Real_secrecy_rate_curve_AO(2:end), 'Color', colors(2,:), 'LineStyle',':', 'LineWidth',1.5, 'Marker','^', 'MarkerIndices',1:markerInterval:length(HSCA_Real_secrecy_rate_curve_AO(2:end)), 'MarkerFaceColor',colors(2,:));

% PSO
plot(PSO_Fake_secrecy_rate_curve, 'k--', 'LineWidth',1.8, 'Marker','v', 'MarkerIndices',1:markerInterval:length(PSO_Fake_secrecy_rate_curve), 'MarkerFaceColor','k');
plot(PSO_Real_secrecy_rate_curve, 'k-.', 'LineWidth',1.8, 'Marker','>', 'MarkerIndices',1:markerInterval:length(PSO_Real_secrecy_rate_curve), 'MarkerFaceColor','k');

title('Best Fake & Real Private Secrecy Rate','FontWeight','bold','FontSize',12);
xlabel('Iteration','FontWeight','bold','FontSize',11);
ylabel('Secrecy Rate','FontWeight','bold','FontSize',11);
legend( ...
    'SCA-fake','SCA-fake-AO','SCA-real','SCA-real-AO', ...
    'HSCA-fake','HSCA-fake-AO','HSCA-real','HSCA-real-AO', ...
    'PSO-fake','PSO-real', ...
    'Location','best','FontSize',10);

grid on;
ax = gca;
ax.GridAlpha = 0.3;
ax.LineWidth = 1.1;
box on;



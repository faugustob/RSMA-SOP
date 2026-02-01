clc; clear;

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

Num_agents  = 50;
Max_iteration = 1000;
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

%%
% phi1 = rand(1,Nr)*2*pi;
% phi2 = phi1; 
% 
% Nc1 = compute_OTFS_static_channel(0,Pe,P,Q_j,Plos(1),PLj(1),Nr,HB(:,:,:,1),HA(:,:,:,:,1),g_pq(:,:,1),phi1,Nsymb,h_rp(:,:,1),h_jq(:,:,1),h_e(:,1),'loop')
% Nc2 = compute_OTFS_static_channel(0,Pe,P,Q_j,Plos(1),PLj(1),Nr,HB(:,:,:,1),HA(:,:,:,:,1),g_pq(:,:,1),phi2,Nsymb,h_rp(:,:,1),h_jq(:,:,1),h_e(:,1),'vectorized')


%%
figure;
plot(Convergence_curve(2:end),'Color','b','LineWidth',1.5)
title('SCA Convergence curve')
xlabel('Iteration');
ylabel('Best flame (score) obtained so far');

figure;
plot(Fake_secrecy_rate_curve(2:end),'Color','b','LineWidth',1.5)
title('Best SCA fake private secrecy rate curve')
xlabel('Iteration');
ylabel('Best fake mean private secrecy rate curve');

figure;
plot(Fake_secrecy_rate_curve(2:end),'Color','b','LineWidth',1.5)
title('Best SCA  real private secrecy rate curve')
xlabel('Iteration');
ylabel('Best real fake mean private secrecy rate curve');



%% Hybrid SCA optimization (Brajevic et al)
display('Hybrid SCA is optimizing your problem');

SP = 20;
MNI = 400;

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

best_mean_fake_secrecy_rate=-inf;
best_mean_real_secrecy_rate=-inf;
best_objective = -inf;

objective = zeros(SP,1); % sum rate of each agent in population
population_real_secrecy_rate = zeros(SP,1); % sum rate of each agent in population
population_fake_secrecy_rate = zeros(SP,1);

t = 1; % iteration counter
a = 0.75; % constant
r1 = a; % parameter that regulates the location of the novel solution
MR = 0.1; % modification rate control parameter
MR_max = 0.9;
P_control = 0.3;

while t<=MNI
    % find and record the best solution y reached so far
    for i=1:SP
        [sc_c_lk,sc_p_lk,sc_p_kk,rate_c,rate_k,R_k,~] = compute_sinr_sc_an(Pe,P,Q_j,nF+L,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,Rmin,h_rp,h_jq,h_e, zeta_k_St,Active_Gain_dB, X(i,:));
        mean_fake_p_secrecy = mean(mean(sc_p_lk(1:nF,:)));
        mean_p_secrecy = mean(mean(sc_p_lk(nF+1:end,:)));

        penalty = 0;
        violation = max(Rmin - R_k, 0);
        penalty = penalty - sum(violation.^2);

        objective(i) = mean_fake_p_secrecy + 1e3 * penalty;
        population_fake_secrecy_rate(i) = mean_fake_p_secrecy;
        population_real_secrecy_rate(i) = mean_p_secrecy;

        

        if objective(i)>best_objective
            y_optimal=X(i,:);
            best_objective = objective(i);
            best_mean_fake_secrecy_rate = population_fake_secrecy_rate(i);
            best_mean_real_secrecy_rate = population_real_secrecy_rate(i);
        end
    end
     display(['HSCA At iteration ', num2str(t), ' the optimum fake sc is ', num2str(best_mean_fake_secrecy_rate), ' the optimum real sc is ', num2str(best_mean_real_secrecy_rate)]);
    if mod(t,2) == 0
        % Produce new solution with modified SCA algorithm
        for i=1:SP
            % sample the random parameters
            neighbours = randsample([1:i-1,i+1:SP], 2);
            r2 = 2*pi*rand(1,size(X,2));
            rand_i = rand();
            R_ij = rand(1,size(X,2));
            under = R_ij < 0.5;
            over = 1 - under;     
            % produce new solution
            % equation (5) from Hybrid SCA paper
            v_i = X(neighbours(1),:)+rand_i.*abs(y_optimal-X(neighbours(2),:))+r1.*abs(y_optimal-X(i,:)).*(under.*sin(r2)+over.*cos(r2));
            
            % to ensure solutions are within the lower and upper bounds:
            % (if this is changed, also change section below)
            % 1. clamp new values to the lower and upper bounds
            % v_i = min(v_i, ub);
            % v_i = max(v_i, lb);
            % or 2. handle with technique from [31] "Hybrid Firefly and
            % Multi-Strategy ABC Algorithm"
            v_i = (v_i<lb_hsca).*(2*lb_hsca-v_i)+(v_i>=lb_hsca).*v_i;
            v_i = (v_i>ub_hsca).*(2*ub_hsca-v_i)+(v_i<=ub_hsca).*v_i;
            % note: there may be issues if new solution is too far from the
            % boundaries, causing the above technique to, for example, 
            % bring the new solution from violating the lower bound to
            % now violating the upper bound

            alpha_hsca  = v_i(1:K+1);
            alpha_hsca  = max(alpha_hsca , alpha_min_hsca);
            alpha_hsca  = alpha_hsca  / sum(alpha_hsca );
            alpha_hsca  = alpha_hsca  - (sum(alpha_hsca ,2)-1)/(K+1);
            alpha_hsca  = alpha_hsca  - (sum(alpha_hsca ,2)-1)/(K+1);
            alpha_hsca  = alpha_hsca  - (sum(alpha_hsca ,2)-1)/(K+1);
             v_i(1:K+1) = alpha_hsca;

            % calculate constraint violations
            % in our case, our only constraint is sum(alpha)<=1
            CV = sum(v_i(1,1:K+1)) > 1;
            if CV > 0
                % if new solution violates constraint, we keep old solution
                continue;
            end
            % for the sum(alpha)<=1 constraint, we can try either:
            % 1. normalize, or
            % 2. just discard solution if it is infeasible

            % evaluate the sum rate
            [sc_c_lk,sc_p_lk,sc_p_kk,rate_c,rate_k,R_k,~] = compute_sinr_sc_an(Pe,P,Q_j,nF+L,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,Rmin,h_rp,h_jq,h_e, zeta_k_St,Active_Gain_dB, v_i);
            mean_fake_p_secrecy = mean(mean(sc_p_lk(1:nF,:)));
            mean_p_secrecy = mean(mean(sc_p_lk(nF+1:end,:)));

            penalty = 0;
            violation = max(Rmin - R_k, 0);
            penalty = penalty - sum(violation.^2);
    
            new_objective = mean_fake_p_secrecy + 1e3 * penalty;
        


            if new_objective > objective(i)
                % if new solution doesn't violate constraint and is better
                % than old solution, we replace old solution with new
                % solution
                objective(i) = new_objective;
                population_fake_secrecy_rate(i) = mean_fake_p_secrecy;
                population_real_secrecy_rate(i) = mean_p_secrecy;
                X(i,:) = v_i;
            end
        end
    else
        % Produce new solution with multi-strategy ABC algorithm
        R_j = rand(1,size(X,2));
        for i=1:SP
            phi_i = rand()*2-1;
            neighbours = randsample([1:i-1,i+1:SP], 2);
            % randomly assign strategy 1 or 2 to each agent in the population
            S_i = rand() < 0.5; % S_i = True is strategy 1 and S_i = False is strategy 2
            if S_i
                % equation (7) from Hybrid SCA paper
                v_i = X(i,:) + (R_j<MR).*phi_i.*(X(i,:)-X(neighbours(1),:));
            else
                % equation (9) from Hybrid SCA paper
                v_i = X(i,:) + phi_i*(X(neighbours(1),:)-X(neighbours(2),:));
            end

            % to ensure solutions are within the lower and upper bounds:
            % (if this is changed, also change section above)
            % 1. clamp new values to the lower and upper bounds
            % v_i = min(v_i, ub);
            % v_i = max(v_i, lb);
            % or 2. handle with technique from [31] "Hybrid Firefly and
            % Multi-Strategy ABC Algorithm"
            v_i = (v_i<lb_hsca).*(2*lb_hsca-v_i)+(v_i>=lb_hsca).*v_i;
            v_i = (v_i>ub_hsca).*(2*ub_hsca-v_i)+(v_i<=ub_hsca).*v_i;
            % note: there may be issues if new solution is too far from the
            % boundaries, causing the above technique to, for example, 
            % bring the new solution from violating the lower bound to
            % now violating the upper bound

            alpha_hsca  = v_i(1:K+1);
            alpha_hsca  = max(alpha_hsca , alpha_min_hsca);
            alpha_hsca  = alpha_hsca  / sum(alpha_hsca );
            alpha_hsca  = alpha_hsca  - (sum(alpha_hsca ,2)-1)/(K+1);
            alpha_hsca  = alpha_hsca  - (sum(alpha_hsca ,2)-1)/(K+1);
            alpha_hsca  = alpha_hsca  - (sum(alpha_hsca ,2)-1)/(K+1);

            v_i(1:K+1) = alpha_hsca;

            % calculate constraint violations
            % in our case, our only constraint is sum(alpha)<=1
            CV = sum(v_i(1,1:K+1)) > 1;
            if CV > 0
                % if new solution violates constraint, we keep old solution
                continue;
            end
            % for the sum(alpha)<=1 constraint, we can try either:
            % 1. normalize, or
            % 2. just discard solution if it is infeasible

            % evaluate the sum rate
            [sc_c_lk,sc_p_lk,sc_p_kk,rate_c,rate_k,R_k,~] = compute_sinr_sc_an(Pe,P,Q_j,nF+L,K,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,Rmin,h_rp,h_jq,h_e, zeta_k_St,Active_Gain_dB, v_i);
            mean_fake_p_secrecy = mean(mean(sc_p_lk(1:nF,:)));
            mean_p_secrecy = mean(mean(sc_p_lk(nF+1:end,:)));

            penalty = 0;
            violation = max(Rmin - R_k, 0);
            penalty = penalty - sum(violation.^2);
    
            new_objective = mean_fake_p_secrecy + 1e3 * penalty;
           

            if new_objective > objective(i)
                % if new solution doesn't violate constraint and is better
                % than old solution, we replace old solution with new
                % solution
                objective(i) = new_objective;
                population_fake_secrecy_rate(i) = mean_fake_p_secrecy;
                population_real_secrecy_rate(i) = mean_p_secrecy;
                X(i,:) = v_i;
            end
        end
    end

    HSCA_Convergence_curve(t)=best_objective;
    HSCA_Fake_secrecy_rate_curve(t)=best_mean_fake_secrecy_rate;
    HSCA_Real_secrecy_rate_curve(t)=best_mean_real_secrecy_rate;


    r1 = a-t*a/MNI;
    P = rand();
    % equation (8) from Hybrid SCA paper
    if MR < MR_max
        MR = MR + (MR_max-0.1)/(P*MNI);
        % note: MR may exceed MR_max for one iteration
    else
        MR = MR_max;
    end
    t = t + 1;

    
end

best_HSCA = best_mean_fake_secrecy_rate;

figure;
plot(HSCA_Convergence_curve(2:end),'Color','b','LineWidth',1.5)
title('HSCA Convergence curve')
xlabel('Iteration');
ylabel('Best flame (score) obtained so far');

figure;
plot(HSCA_Fake_secrecy_rate_curve(2:end),'Color','b','LineWidth',1.5)
title('Best HSCA fake private secrecy rate curve')
xlabel('Iteration');
ylabel('Best fake mean private secrecy rate curve');

figure;
plot(HSCA_Real_secrecy_rate_curve(2:end),'Color','b','LineWidth',1.5)
title('Best HSCA  real private secrecy rate curve')
xlabel('Iteration');
ylabel('Best real fake mean private secrecy rate curve');

%% Particle Swarm Optimization

display('Using MATLAB built-in particleswarm for optimization...');



% --- PSO Options ---
% You can adjust SwarmSize and MaxIterations to match your original SCA settings
options = optimoptions('particleswarm', ...
    'SwarmSize', 50, ...
    'MaxIterations', 1000, ...
    'Display', 'iter', ...
    'PlotFcn', @(optimValues,state) myCustomPlot(optimValues,state));

% --- Define the Objective Function Wrapper ---
% We pass all your environment variables into the function handle
fitness_func = @(x) objective_wrapper(x, K, Nr, Rmin, Pe, P, Q_j, nF,L, delta_f, Plos, PLj, HB, HA, g_pq, Nsymb, reflect, h_rp, h_jq, h_e,Active_Gain_dB);

% --- Run Built-in PSO ---
[best_x, best_fval] = particleswarm(fitness_func, dim_pso, lb_pso, ub_pso, options);

% --- Post-Processing ---
% Extract final best sum rate from the best position found
[~, final_sec_rate,mean_p_secrecy] = objective_wrapper(best_x, K, Nr, Rmin, Pe, P, Q_j,nF, L, delta_f, Plos, PLj, HB, HA, g_pq, Nsymb, reflect, h_rp, h_jq, h_e,Active_Gain_dB);

display(['Optimization complete. Best min. fake private sec. Rate: ', num2str(final_sec_rate),'Best min. real private sec. Rate: ', num2str(mean_p_secrecy)]);

% --- The Objective Function Wrapper (Local Function) ---
function [fitness, min_fake_p_secrecy,min_p_secrecy] = objective_wrapper(x, K, Nr, Rmin, Pe, P, Q_j,nF, L,  delta_f, Plos, PLj, HB, HA, g_pq, Nsymb, reflect, h_rp, h_jq, h_e, Active_Gain_dB)

     % 1. Extract and Normalize Alpha (ensure sum = 1)
    alpha = x(1:K+1);
    alpha = alpha ./ (sum(alpha));

    alpha = alpha - (sum(alpha)-1)/(K+1);
    
    % 2. Extract Phi
    
    phi_St = x(K+2:K+1+Nr);
    zeta_k_St = (10^(Active_Gain_dB/10)) * ones(1, Nr);

   X =  [alpha, phi_St];
 
   any_reflect = any(reflect > 0) && any(reflect < 0);


   if any_reflect
        phi_Sr =  x(K+2+Nr:K+1+2*Nr);
        zeta_k_St = (10^(Active_Gain_dB/10)) * x(K+2+2*Nr:K+1+3*Nr);
        X = [alpha,phi_St,phi_Sr,zeta_k_St];
    end
    % 3. Call your SINR function
    % We pass normalized alpha and phi back into the compute function
     [sc_c_lk,sc_p_lk,sc_p_kk,rate_c,rate_k,R_k,~]  = compute_sinr_sc_an(Pe, P, Q_j,nF + L, K,  delta_f, Plos, PLj, Nr, HB, HA, g_pq, Nsymb, ...
        reflect,Rmin,h_rp, h_jq, h_e, zeta_k_St,Active_Gain_dB, X);
    min_fake_p_secrecy = mean(mean(sc_p_lk(1:nF,:)));
    min_p_secrecy = mean(mean(sc_p_lk(nF+1:end,:)));

    % 5. Penalty for Rmin constraint violation
    penalty = 1e3 * sum(max(Rmin - R_k, 0).^2);
    
    % Objective: Minimize -SumRate + Penalty
    fitness = -min_fake_p_secrecy + penalty;
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
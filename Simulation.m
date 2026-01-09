clc; clear;

Ns = 1e6; % number of samples for Monte Carlo simulation
rng(22);

c = physconst('LightSpeed'); % speed of light
f_c = 10e9; % frequency
lambda = c / f_c; % wavelength
 % ---- Earth-centered conversion ----
R_earth = 6371e3;


%Number of nodes
K_h = 1;  % number of high speed legit users % 50 km/h
K_s = 2;  % number of high speed legit users % 1.2 m/s
K = K_h+K_s; % number of legit users

L = 4; % number of eavesdroppers

% --- OTFS System Parameters ---
delta_f = 100e3;      % Subcarrier spacing (Hz)
T = 1/delta_f;       % Symbol duration
B = 10e6;        % [Hz] ← Use this
Tf      = 14*T;      % 14-symbol frame (~1 ms)


L_tau = 8;   % 8 delay taps over max_tau (covers multipath + RIS)
L_nu  = 8;   % 8 Doppler taps over [-max_nu, max_nu]

R = 10;

m_rician = (R+1)^2/(2*R+1);

N_V = 16; % number of rows of regularly arranged unit cells of RIS
N_H = 12; % number of columns of regularly arranged unit cells of RIS
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
HAP_altitude = 20e3;
LEO_altitude = 250e3;

RIS_normal = [0;0;1];


[max_theta,max_alpha,max_beta] = compute_max_theta(LEO_altitude,HAP_altitude);

S_elevation = pi/2-max_theta*(0.009);
S_azimuth  = pi/4;
S_radius  = R_earth+LEO_altitude;

S_sph = [ S_azimuth, S_elevation, S_radius]; % location of LEO satellite in spherical coordinates

[S_x, S_y, S_z] = sph2cart(S_sph(1), S_sph(2), S_sph(3));
S_xyz = [S_x;S_y;S_z];


S_v = 7800;
R_xyz = [0; 0; R_earth+HAP_altitude]; % location of STAR-RIS; code assumes this to be origin;
% note: this code assumes surface is on x-y plane (surface normal points in
% z-axis direction)
R_v_xyz = [0;0;0];


% ELEVATION (UNCHANGED)
el = pi/2 - (0.1)*max_alpha * rand(1, K);


az = (2*pi) * rand(1, K);


% RADIUS (FIXED ON EARTH SURFACE)
r = R_earth * ones(1, K);

% SPHERICAL → CARTESIAN
[x, y, z] = sph2cart(az, el, r);

ground_users_cart = [x; y; z];   % 3 x K

%% Eavesdroppers position
x = 1000*rand(1, L);
y = 1000*rand(1, L);
z = R_earth + 50 + 950*rand(1, L);


%% SHIFT TO GLOBAL COORDINATES
% Add the UAV center position to each point
eavesdroppers_xyz =  [x; y; z];

rho_j_xyz = [ground_users_cart,eavesdroppers_xyz];

% find out whether each receiver is on the reflect side or transmit side
reflect = sign(RIS_normal.' * (rho_j_xyz - R_xyz));



% Nakagami Parameters:

% from LEO satellite to STAR-RIS (one value for each path and assume same
% for each RIS element):
m_p = [m_rician;1*ones(P-1,1)]; % shape parameter
omega_p = (1/P)*ones(1,P); % spread parameter

% from STAR-RIS to users and eavesdroppers (one value for each path and 
% each receiver, and assume same for each RIS element):
m_q = 1*ones(K+L,Q_j); % shape parameter
omega_q = (1/Q_j)*ones(K+L,Q_j); % spread parameter
% note: m_j_1(1:K,:) and omega_j_q(1:K,:) are for legit users
% and m_j_1(K+1:end,:) and omega_j_q(K+1:end,:) are for eavesdroppers
% direct channel from LEO satellite to Eve (one value for each path and
% each eavesdropper)

m_e = 1*ones(K+L,Pe); % shape parameter
omega_e = (1/Pe)*ones(K+L,Pe); % spread parameter
% note: if we want different paths to have different parameters, or legit
% users and eavesdroppers to have different parameters, modify the above
% matrices correspondingly, but ensure the shape doesn't change


% STAR-RIS phase and magnitude of each RIS element


% note: beta_r will have shape (Ns,Nr,K+L)



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


g_pq = zeros(P,Q_j,K+L);
Plos = zeros(K+L,1);
PLj = zeros(K+L,1);

h_rp = zeros(Nr, P,K+L);
h_jq = zeros(Nr, Q_j,K+L);
h_e = zeros(Pe,K+L);
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
    d_ris = sqrt(sum((User_k_loc-R_xyz).^2));


    Plos(k) = ((lambda/(4*pi))^2*G_t*G_r)/(d_los^2);
    Plos_ris(k) = ((lambda/(4*pi))^2*G_t*G_r)/(d_ris^2);

    PL = Plos(k)*Plos_ris(k) ;
    
    PLj(k) = compute_ris_PL(lambda,N_V,N_H,S_xyz,User_k_loc,R_xyz,RIS_normal,F,F_tx,F_rx,G,G_t,G_r);
    Ns=1;
  

    % direction in xy (away from ris)
    d_ru = User_k_loc(1:2) - R_xyz(1:2);
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
        h_rp(:, p,k) = sqrt(gamrnd(m_p(p), omega_p(p)/m_p(p), [Nr, 1])) .* exp(1i*2*pi*rand(Nr, 1));
    end

    
    for q = 1:Q_j
        h_jq(:, q,k) = sqrt(gamrnd(m_q(k,q), omega_q(q)/m_q(k,q), [Nr, 1])) .* exp(1i*2*pi*rand(Nr,1));
    end

    
    for u = 1:Pe
        h_e(u,k) = sqrt(gamrnd(m_e(k,u), omega_e(u)/m_e(k,u))) .* exp(1i*2*pi*rand());
    end
            
end

% % % Max tau and nu
max_tau = max([taus_kq(:);taus_ku(:)])-min([taus_kq(:);taus_ku(:)]); 
max_nu  = max([nus_kq(:);nus_ku(:)])-min([nus_kq(:);nus_ku(:)]);    


% Compute M and N based on the parameters
[M, N] = computeOTFSgrid(max_tau, max_nu, 'numerology', B, delta_f, T, Tf);
M = max(M, 64); N = max(N, 12);  % Minimum practical size



Nsymb = M*N; 

HA = zeros(Nsymb,Nsymb,P,Q_j,K+L); % Relay link
HB = zeros(Nsymb,Nsymb,Pe,K+L);


for k=1:K

    for p=1:P
        for q=1:Q_j
            HA(:,:,p,q,k) =  compute_Hp(taus_kq(k,p,q), nus_kq(k,p,q), M, N, T, delta_f, 'blocked');
        end
    end
   

    for u = 1:Pe
        HB(:,:,u,k) =  compute_Hp(taus_u(u), nus_u(u), M, N, T, delta_f, 'blocked');
    end

end



for l=1:L

        vl_ms = 14;        

        User_l_loc = rho_j_xyz(:,K+l);

        d_los = sqrt(sum((User_l_loc-S_xyz).^2));

        Plos(K+l) = ((lambda/(4*pi))^2*G_t*G_r)/(d_los^2);


        PLj(K+l) = compute_ris_PL(lambda,N_V,N_H,S_xyz,User_l_loc,R_xyz,RIS_normal,F,F_tx,F_rx,G,G_t,G_r);

         
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
            h_rp(:, p,K+l) = sqrt(gamrnd(m_p(p), omega_p(p)/m_p(p), [Nr, 1])) .* exp(1i*2*pi*rand(Nr, 1));
        end
    
    
        
        for q = 1:Q_j
            h_jq(:, q,K+l) = sqrt(gamrnd(m_q(K+l,q), omega_q(q)/m_q(K+l,q), [Nr, 1])) .* exp(1i*2*pi*rand(Nr,1));
        end
    
        
        for u = 1:Pe
            h_e(u,K+l) = sqrt(gamrnd(m_e(K+l,u), omega_e(u)/m_e(K+l,u))) .* exp(1i*2*pi*rand());
        end
    

       for p=1:P
            for q=1:Q_j
                HA(:,:,p,q,K+l) =  compute_Hp(taus_R(p)+taus_l(q), nus_R(p)+nus_l(q), M, N, T, delta_f, 'blocked');
            end
       end

        for u = 1:Pe
            HB(:,:,u,K+l) =  compute_Hp(taus_u_l(u) ,nus_u_l(u), M, N, T, delta_f, 'blocked');
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

Num_agents  = 30;
Max_iteration = 100;
Rmin=0.01;

% Check if more than one STAR-RIS side is being used.
any_reflect = any(reflect > 0) && any(reflect < 0);

% Problem bounds and dimensionality
dim = K+1+Nr;
ub=[ones(1,K+1),2*pi*ones(1,Nr)];
alpha_min = 1e-4;
lb = [alpha_min * ones(1,K+1),zeros(1,Nr)];
zeta_k_St = ones(1,Nr); % RIS amplitude coefficients, we may use it to boost for active RIS


% zeta_k_Sr = rand(Num_agents,Nr); % reflection coefficients
phi = 2*pi*rand(Num_agents,Nr);


alpha = (1/(K+1))*rand(Num_agents, K+1); 
% random values
alpha = alpha ./ sum(alpha, 2);      % divide each row by its row sum
X = [alpha,phi];

if any_reflect
    dim = K+1+2*Nr;
    ub=[ones(1,K+1),2*pi*ones(1,Nr),ones(1,Nr)];
    alpha_min = 1e-4;
    lb = [alpha_min * ones(1,K+1),zeros(1,2*Nr)];
    zeta_k_St = rand(Num_agents,Nr);
    X = [alpha,phi,zeta_k_St];
end

% --- Problem Dimensions and Bounds ---
dim_pso = dim;
alpha_min_pso = alpha_min;
lb_pso =lb;
ub_pso = ub;



Destination_position=zeros(1,dim);
Destination_fitness=inf;
best_sum_rate=0;

Convergence_curve=zeros(1,Max_iteration);
sum_rate_curve=zeros(1,Max_iteration);

min_sum_secrecy = zeros(1,Max_iteration);

Objective_values = zeros(1,size(X,1));
%All_objective_values=zeros(Max_iteration,size(X,1));


% Calculate the fitness of the first set and find the best one
for i=1:size(X,1)
    C_k = zeros(K,1);

    [sc_c_lk,sc_p_lk,rate_c,rate_k,~] = compute_sinr_sc(Pe,P,Q_j,L,K,m_e,m_q,m_p,omega_e,omega_p,omega_q,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,h_rp,h_jq,h_e,alpha,phi,zeta_k_St,X(i,:));

    sum_secrecy = sc_c_lk+sc_p_lk; %Private + Common secrecy capacities.
    sum_rate_k = rate_c + sum(rate_k(:));
    % Objective_values(1,i)=-mean(sum_secrecy(:));

    rate_c_available = rate_c;

    for k = 1:K
        deficit = max(Rmin - rate_k(k), 0);
        C_k(k) = min(deficit, rate_c_available);
        rate_c_available = max(rate_c_available - C_k(k), 0);
    end
    
    R_k = rate_k(:) + C_k;
    sum_rate_k = sum(R_k);


    penalty = 0;
    violation = max(Rmin - R_k, 0);
    penalty = penalty + sum(violation.^2);

    
    Objective_values(1,i) = -sum_rate_k + 1e3 * penalty;

    if i==1
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
        best_sum_rate = sum_rate_k;
    elseif Objective_values(1,i)<Destination_fitness
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
        best_sum_rate = sum_rate_k;
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
        C_k = zeros(K,1);

        % Check if solutions go outside the search spaceand bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

       X(i,1:K+1) = X(i,1:K+1) ./ (sum(X(i,1:K+1), 2)+ 1e-15); % Normalize to ensure sum alpha = 1;



        % Calculate the objective values

        [sc_c_lk,sc_p_lk,rate_c,rate_k,~] = compute_sinr_sc(Pe,P,Q_j,L,K,m_e,m_q,m_p,omega_e,omega_p,omega_q,delta_f,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,h_rp,h_jq,h_e,alpha,phi, zeta_k_St, X(i,:));
        sum_secrecy = sc_c_lk+sc_p_lk; %Private + Common secrecy capacities.
        sum_rate_k = rate_c + sum(rate_k(:));

        %Objective_values(1,i)=-mean(sum_secrecy(:));

         rate_c_available = rate_c;

        for k = 1:K
            deficit = max(Rmin - rate_k(k), 0);
            C_k(k) = min(deficit, rate_c_available);
            rate_c_available = max(rate_c_available - C_k(k), 0);
        end
        
        R_k = rate_k(:) + C_k;
        sum_rate_k = sum(R_k);

    
        penalty = 0;
        violation = max(Rmin - R_k, 0);
        penalty = penalty + sum(violation.^2);
    
        
        Objective_values(1,i) = -sum_rate_k + 1e3 * penalty;
        %Objective_values(1,i)=-sum_rate_k;


        % Update the destination if there is a better solution
        if Objective_values(1,i)<Destination_fitness
            Destination_position=X(i,:);
            Destination_fitness=Objective_values(1,i);
            best_sum_rate = sum_rate_k;
        end
    end

    Convergence_curve(t)=-Destination_fitness;
    Sum_rate_curve(t)=best_sum_rate;


    % Display the iteration and best optimum obtained so far
    if mod(t,1)==0
        display(['At iteration ', num2str(t), ' the optimum is ', num2str(best_sum_rate)]);
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
title('Convergence curve')
xlabel('Iteration');
ylabel('Best flame (score) obtained so far');

figure;
plot(Sum_rate_curve(2:end),'Color','b','LineWidth',1.5)
title('Best sum rate curve')
xlabel('Iteration');
ylabel('Best flame (score) obtained so far');

%% Particle Swarm Optimization

display('Using MATLAB built-in particleswarm for optimization...');



% --- PSO Options ---
% You can adjust SwarmSize and MaxIterations to match your original SCA settings
options = optimoptions('particleswarm', ...
    'SwarmSize', 50, ...
    'MaxIterations', 50, ...
    'Display', 'iter', ...
    'PlotFcn', @(optimValues,state) myCustomPlot(optimValues,state));

% --- Define the Objective Function Wrapper ---
% We pass all your environment variables into the function handle
fitness_func = @(x) objective_wrapper(x, K, Nr, Rmin, Pe, P, Q_j, L, m_e, m_q, m_p, ...
    omega_e, omega_p, omega_q, delta_f, Plos, PLj, HB, HA, g_pq, Nsymb, reflect, h_rp, h_jq, h_e);

% --- Run Built-in PSO ---
[best_x, best_fval] = particleswarm(fitness_func, dim_pso, lb_pso, ub_pso, options);

% --- Post-Processing ---
% Extract final best sum rate from the best position found
[~, final_sum_rate] = objective_wrapper(best_x, K, Nr, Rmin, Pe, P, Q_j, L, m_e, m_q, m_p, ...
    omega_e, omega_p, omega_q, delta_f, Plos, PLj, HB, HA, g_pq, Nsymb, reflect, h_rp, h_jq, h_e);

display(['Optimization complete. Best Sum Rate: ', num2str(final_sum_rate)]);

% --- The Objective Function Wrapper (Local Function) ---
function [fitness, sum_rate_k] = objective_wrapper(x, K, Nr, Rmin, Pe, P, Q_j, L, m_e, m_q, m_p, ...
    omega_e, omega_p, omega_q, delta_f, Plos, PLj, HB, HA, g_pq, Nsymb, reflect, h_rp, h_jq, h_e)

    % 1. Extract and Normalize Alpha (ensure sum = 1)
    alpha = x(1:K+1);
    alpha = alpha ./ (sum(alpha) + 1e-15);
    
    % 2. Extract Phi
    phi = x(K+2:end);
    zeta_k_Sr = ones(1, Nr); % Constant RIS amplitude
    
    % 3. Call your SINR function
    % We pass normalized alpha and phi back into the compute function
    [~, ~, rate_c, rate_k, ~] = compute_sinr_sc(Pe, P, Q_j, L, K, m_e, m_q, m_p, ...
        omega_e, omega_p, omega_q, delta_f, Plos, PLj, Nr, HB, HA, g_pq, Nsymb, ...
        reflect,h_rp, h_jq, h_e, alpha, phi, zeta_k_Sr, [alpha, phi]);

    % 4. Calculate Rate with Common Rate Allocation
    rate_c_available = rate_c;
    C_k = zeros(K, 1);
    for k = 1:K
        deficit = max(Rmin - rate_k(k), 0);
        C_k(k) = min(deficit, rate_c_available);
        rate_c_available = max(rate_c_available - C_k(k), 0);
    end
    
    R_k = rate_k(:) + C_k;
    sum_rate_k = sum(R_k);
    
    % 5. Penalty for Rmin constraint violation
    penalty = 1e3 * sum(max(Rmin - R_k, 0).^2);
    
    % Objective: Minimize -SumRate + Penalty
    fitness = -sum_rate_k + penalty;
end

function stop = myCustomPlot(optimValues, state)
    stop = false;
    % We plot the negative of the best fval to see the actual Sum Rate
    plot(optimValues.iteration, -optimValues.bestfval, 'bo', 'MarkerFaceColor', 'b');
    xlabel('Iteration');
    ylabel('Actual Sum Rate');
    grid on;
    title('Maximization Progress');
    hold on;
end
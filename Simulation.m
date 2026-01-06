clc; clear;

Ns = 1e6; % number of samples for Monte Carlo simulation
rng(2);

c = physconst('LightSpeed'); % speed of light
f_c = 10.5e9; % frequency
lambda = c / f_c; % wavelength
 % ---- Earth-centered conversion ----
R_earth = 6371e3;


%Number of nodes
K_h = 1;  % number of high speed legit users % 50 km/h
K_s = 1;  % number of high speed legit users % 1.2 m/s
K = K_h+K_s; % number of legit users

L = 1; % number of eavesdroppers

% --- System Parameters ---
delta_f = 15e3;      % Subcarrier spacing (Hz)
T = 1/delta_f;       % Symbol duration

M=8;
N=8;

I=1;

R = 10;

m_rician = (R+1)^2/(2*R+1);

N_V = 4; % number of rows of regularly arranged unit cells of RIS
N_H = 4; % number of columns of regularly arranged unit cells of RIS
Nr = N_V * N_H; % total number of unit cells of RIS

d_x = floor(lambda/2 * 1000) / 1000; % horizontal size of RIS element
d_y = floor(lambda/2 * 1000) / 1000; % vertical size of RIS element


G_tx = 10^(10/10); %10 dB
G_rx = 10^(10/10); %10 dB




P = 2; % number of propagation paths arriving at the r-th RIS element
Q_j = 3; % number of propagation paths departing from the r-th RIS element
Pe = 2; % number of propagation paths departing from the LEO satellite

% location coordinates
UAV_altitude = 150;
LEO_altitude = 480e3;


[max_theta,max_alpha,max_beta] = compute_max_theta(LEO_altitude,UAV_altitude);

S_elevation = max_theta/2;
S_azimuth  = pi/4;
S_radius  = R_earth+LEO_altitude;

S_sph = [ S_azimuth, S_elevation, S_radius]; % location of LEO satellite in spherical coordinates

[S_x, S_y, S_z] = sph2cart(S_sph(1), S_sph(2), S_sph(3));
S_xyz = [S_x;S_y;S_z];


S_v = 7800;
R_xyz = [0; 0; R_earth+UAV_altitude]; % location of STAR-RIS; code assumes this to be origin;
% note: this code assumes surface is on x-y plane (surface normal points in
% z-axis direction)
R_v_xyz = [0;0;0];


% ELEVATION (UNCHANGED)
el = max_alpha * rand(1, K);

% AZIMUTH (FORCE HALF x>0, HALF x<0)
Khalf = floor(K/2);

if Khalf == 0
    Khalf=1;
end

% x > 0  → cos(az) > 0
az_pos = (pi/2) * rand(1, Khalf) - pi/2;

% x < 0  → cos(az) < 0
az_neg = (pi/2) * rand(1, K - Khalf) + pi/2;

% Combine and shuffle
az = [az_pos, az_neg];
az = az(randperm(K));

% RADIUS (FIXED ON EARTH SURFACE)
r = R_earth * ones(1, K);

% SPHERICAL → CARTESIAN
[x, y, z] = sph2cart(az, el, r);

ground_users_cart = [x; y; z];   % 3 x K

% PARAMETERS
Rmin = 20;                 % Inner radius of shell [m]
Rmax = 100;                 % Outer radius of shell [m]


% Center of the sphere (UAV position in ECEF-like XYZ)
center = [0; 0; R_earth + UAV_altitude];


% RANDOM VARIABLES
u = rand(1, L);   % radius
w = rand(1, L);   % polar angle

% RADIUS SAMPLING (UNIFORM VOLUME)
r = (Rmin^3 + (Rmax^3 - Rmin^3).*u).^(1/3);

% POLAR ANGLE (ISOTROPIC)
costheta = 2*w - 1;
theta = acos(costheta);

% AZIMUTH ANGLE (FORCE HALF x>0, HALF x<0)

Lhalf = floor(L/2);

if Lhalf==0
    Lhalf=1;
end

% x > 0  → cos(phi) > 0
phi_pos = (pi/2) * rand(1, Lhalf) - pi/2;

% x < 0  → cos(phi) < 0
phi_neg = (pi/2) * rand(1, L - Lhalf) + pi/2;

% Combine and shuffle
phi = [phi_pos, phi_neg];
phi = phi(randperm(L));

%% SPHERICAL → CARTESIAN
x = r .* sin(theta) .* cos(phi);
y = r .* sin(theta) .* sin(phi);
z = r .* cos(theta);


%% SHIFT TO GLOBAL COORDINATES
% Add the UAV center position to each point
eavesdroppers_xyz = R_xyz + [x; y; z];

rho_j_xyz = [ground_users_cart,eavesdroppers_xyz];

% find out whether each receiver is on the reflect side or transmit side
reflect = sign(rho_j_xyz(1,:)) == sign(S_xyz(1)); % shape (1,K+L,Ns)
%transmit = 1 - reflect;  % shape (1,K+L,Ns)


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

Nsymb = M*N; 

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

sigma_ang = deg2rad(20);   % angular spread

% Satellite to RIS delays and doppler coefficients.
[taus_R, nus_R, u_paths_R] = compute_delay_and_doppler( ...
    c, S_xyz, vS, R_xyz, vR, f_c, P, sigma_ang);


HA = zeros(Nsymb,Nsymb,P,Q_j,K+L); % Relay link
HB = zeros(Nsymb,Nsymb,Pe,K+L);
g_pq = zeros(P,Q_j,K+L);
Plos = zeros(K+L,1);
PLj = zeros(K+L,1);

h_rp = zeros(Nr, P,K+L);
h_jq = zeros(Nr, Q_j,K+L);
h_e = zeros(Pe,K+L);

for k =1:K

  
    if k <= K_h
        vk_ms = 50 / 3.6;
    else
        vk_ms = 1.2;
    end


    User_k_loc = rho_j_xyz(:,k);

    d_los = sqrt(sum((User_k_loc-S_xyz).^2));


    Plos(k) = ((lambda/(4*pi))^2*G_tx*G_rx)/(d_los^2);
    
    PLj(k) = compute_ris_PL(lambda,N_V,N_H,S_xyz,User_k_loc,R_xyz);
    Ns=1;
  

    % direction in xy (away from ris)
    d_ru = User_k_loc(1:2) - R_xyz(1:2);
    d_ru = d_ru / norm(d_ru);

    % receiver velocity (guaranteed norm)
    v_l  = vk_ms * [d_ru; 0];

    % RIS to legitimate users delays and doppler coefficients.
    [taus_k, nus_k, u_paths_k] = compute_delay_and_doppler( ...
    c, R_xyz, vR, User_k_loc, v_l, f_c, Q_j, sigma_ang);

    % sat to legitimate users delays and doppler coefficients.
    [taus_u, nus_u, u_paths_u] = compute_delay_and_doppler( ...
    c, S_xyz, vS, User_k_loc, v_l, f_c, Pe, sigma_ang);

    g_pq(:,:,k) = exp(1i*2*pi*(taus_R*nus_k'));    



   
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
       

    for p=1:P
        for q=1:Q_j
            HA(:,:,p,q,k) =  compute_Hp(taus_R(p)+taus_k(q), nus_R(p)+nus_k(q), M, N, T, delta_f, 'blocked');
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

        Plos(K+l) = ((lambda/(4*pi))^2*G_tx*G_rx)/(d_los^2);


        PLj(K+l) = compute_ris_PL(lambda,N_V,N_H,S_xyz,User_l_loc,R_xyz);

         
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
        [taus_u, nus_u, u_paths_u] = compute_delay_and_doppler( ...
        c, S_xyz, vS, User_k_loc, v_l, f_c, Pe, sigma_ang);


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
            HB(:,:,u,K+l) =  compute_Hp(taus_u(u), nus_u(u), M, N, T, delta_f, 'blocked');
        end

        g_pq(:,:,K+l) = exp(1i*2*pi*(taus_R*nus_l'));           
end
 
%  % SCA optimization SINR, Secrecy capacity
% 
% display('SCA is optimizing your problem');
% 
% Num_agents  = 30;
% Max_iteration = 500;
% 
% %Initialize the set of random solutions
% %X=initialization(N,dim,ub,lb);
% 
% dim = K+1;
% 
% ub=ones(1,K+1);
% lb=zeros(1,K+1);
% 
% zeta_k_Sr = rand(Nr,1); % reflection coefficients
% phi = rand(Nr,1);
% 
% alpha = rand(Num_agents, K+1);       % random values
% alpha = alpha ./ sum(alpha, 2);      % divide each row by its row sum
% 
% X = [alpha];
% 
% Destination_position=zeros(1,dim);
% Destination_fitness=inf;
% 
% Convergence_curve=zeros(1,Max_iteration);
% Objective_values = zeros(1,size(X,1));
% %All_objective_values=zeros(Max_iteration,size(X,1));
% 
% % Calculate the fitness of the first set and find the best one
% for i=1:size(X,1)
%     [sc_c_lk,sc_p_lk,rate_c,rate_k,~] = compute_sinr_sc(Pe,P,Q_j,L,K,m_e,m_j_q,m_p,omega_e,omega_p,omega_j_q,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,zeta_k_Sr,phi,X(i,:));
%     sum_secrecy = sc_c_lk+sc_p_lk; %Private + Common secrecy capacities.
%     Objective_values(1,i)=-min(sum_secrecy(:));
%     if i==1
%         Destination_position=X(i,:);
%         Destination_fitness=Objective_values(1,i);
%     elseif Objective_values(1,i)<Destination_fitness
%         Destination_position=X(i,:);
%         Destination_fitness=Objective_values(1,i);
%     end
% 
%     All_objective_values(1,i)=Objective_values(1,i);
% end
% 
% %Main loop
% t=2; % start from the second iteration since the first iteration was dedicated to calculating the fitness
% while t<=Max_iteration
% 
%     % Eq. (3.4)
%     a = 2;
%     Max_iteration = Max_iteration;
%     r1=a-t*((a)/Max_iteration); % r1 decreases linearly from a to 0
% 
%     % Update the position of solutions with respect to destination
%     for i=1:size(X,1) % in i-th solution
%         for j=1:size(X,2) % in j-th dimension
% 
%             % Update r2, r3, and r4 for Eq. (3.3)
%             r2=(2*pi)*rand();
%             r3=2*rand;
%             r4=rand();
% 
%             % Eq. (3.3)
%             if r4<0.5
%                 % Eq. (3.1)
%                 X(i,j)= X(i,j)+(r1*sin(r2)*abs(r3*Destination_position(j)-X(i,j)));
%             else
%                 % Eq. (3.2)
%                 X(i,j)= X(i,j)+(r1*cos(r2)*abs(r3*Destination_position(j)-X(i,j)));
%             end
% 
%         end
%     end
% 
%     for i=1:size(X,1)
% 
%         % Check if solutions go outside the search spaceand bring them back
%         Flag4ub=X(i,:)>ub;
%         Flag4lb=X(i,:)<lb;
%         X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
% 
%        X(i,:) = X(i,:) ./ sum(X(i,:), 2); % Normalize to ensure sum alpha = 1;
% 
% 
% 
%         % Calculate the objective values
% 
%         [sc_c_lk,sc_p_lk,rate_c,rate_k,~] = compute_sinr_sc(Pe,P,Q_j,L,K,m_e,m_j_q,m_p,omega_e,omega_p,omega_j_q,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,zeta_k_Sr,phi,X(i,:));
%         sum_secrecy = sc_c_lk+sc_p_lk; %Private + Common secrecy capacities.
%         Objective_values(1,i)=-min(sum_secrecy(:));
% 
%         % Update the destination if there is a better solution
%         if Objective_values(1,i)<Destination_fitness
%             Destination_position=X(i,:);
%             Destination_fitness=Objective_values(1,i);
%         end
%     end
% 
%     Convergence_curve(t)=-Destination_fitness;
% 
%     % Display the iteration and best optimum obtained so far
%     if mod(t,1)==0
%         display(['At iteration ', num2str(t), ' the optimum is ', num2str(-Destination_fitness)]);
%     end
% 
%     % Increase the iteration counter
%     t=t+1;
% end
%%
% dim = K+1+2*Nr;
% fun = @(x) -min(compute_sinr_sc(Pe,P,Q_j,L,K,m_e,m_q,m_p,omega_e,omega_p,omega_q,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,h_rp,h_jq,h_e,x));  % Wrap with min
% lb = zeros(1,dim); ub = [ones(1,K+1), 2*pi*ones(1,Nr), ones(1,Nr)];
% options = optimoptions('particleswarm','SwarmSize',100,'MaxIterations',1000,'Display','iter');
% [x_opt, fval] = particleswarm(fun, dim, lb, ub, options);

%%
display('SCA is optimizing your problem');

Num_agents  = 5;
Max_iteration = 50;

%Initialize the set of random solutions

dim = K+1+2*Nr;

ub=[ones(1,K+1),2*pi*ones(1,Nr),ones(1,Nr)];
lb=zeros(1,K+1+2*Nr);

zeta_k_Sr = rand(Num_agents,Nr); % reflection coefficients
phi = 2*pi*rand(Num_agents,Nr);

alpha_vec = [1/2, (1/(2*K))*ones(1,K)];
%alpha = rand(Num_agents, K+1);  
alpha = ones(Num_agents, K+1).*alpha_vec; 
% random values
alpha = alpha ./ sum(alpha, 2);      % divide each row by its row sum

X = [alpha,phi,zeta_k_Sr];

Destination_position=zeros(1,dim);
Destination_fitness=inf;

Convergence_curve=zeros(1,Max_iteration);
min_sum_secrecy = zeros(1,Max_iteration);

Objective_values = zeros(1,size(X,1));
%All_objective_values=zeros(Max_iteration,size(X,1));

% Calculate the fitness of the first set and find the best one
for i=1:size(X,1)
    [sc_c_lk,sc_p_lk,rate_c,rate_k,~] = compute_sinr_sc(Pe,P,Q_j,L,K,m_e,m_q,m_p,omega_e,omega_p,omega_q,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,h_rp,h_jq,h_e,X(i,:));
    sum_secrecy = sc_c_lk+sc_p_lk; %Private + Common secrecy capacities.
    sum_rate_k = sum(rate_k(:));
    Objective_values(1,i)=-mean(sum_secrecy(:));
    if i==1
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    elseif Objective_values(1,i)<Destination_fitness
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    end

    All_objective_values(1,i)=Objective_values(1,i);
end




%Main loop
t=2; % start from the second iteration since the first iteration was dedicated to calculating the fitness
while t<=Max_iteration

    % Eq. (3.4)
    a = 2;
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

       X(i,1:K+1) = X(i,1:K+1) ./ (sum(X(i,1:K+1), 2)+ 1e-12); % Normalize to ensure sum alpha = 1;



        % Calculate the objective values

        [sc_c_lk,sc_p_lk,rate_c,rate_k,~] = compute_sinr_sc(Pe,P,Q_j,L,K,m_e,m_q,m_p,omega_e,omega_p,omega_q,Plos,PLj,Nr,HB,HA,g_pq,Nsymb,reflect,h_rp,h_jq,h_e,X(i,:));
        sum_secrecy = sc_c_lk+sc_p_lk; %Private + Common secrecy capacities.
        sum_rate_k = sum(rate_k(:));
        Objective_values(1,i)=-mean(sum_secrecy(:));


        % Update the destination if there is a better solution
        if Objective_values(1,i)<Destination_fitness
            Destination_position=X(i,:);
            Destination_fitness=Objective_values(1,i);
        end
    end

    Convergence_curve(t)=-Destination_fitness;

    % Display the iteration and best optimum obtained so far
    if mod(t,1)==0
        display(['At iteration ', num2str(t), ' the optimum is ', num2str(-Destination_fitness)]);
    end

    % Increase the iteration counter
    t=t+1;

    

end



%%
figure;
plot(Convergence_curve(2:end),'Color','b')
title('Convergence curve')
xlabel('Iteration');
ylabel('Best flame (score) obtained so far');




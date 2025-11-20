num_samples = 100; % number of samples for Monte Carlo simulation

c = physconst('LightSpeed'); % speed of light
f = 10.5e9; % frequency
lambda = c / f; % wavelength

K = 2; % number of legit users
L = 1; % number of eavesdroppers

N_V = 50; % number of rows of regularly arranged unit cells of RIS
N_H = 34; % number of columns of regularly arranged unit cells of RIS
N_R = N_V * N_H; % total number of unit cells of RIS

d_x = floor(lambda/2 * 1000) / 1000; % horizontal size of RIS element
d_y = floor(lambda/2 * 1000) / 1000; % vertical size of RIS element

alpha_c = 0.5; % power allocation factor for common message
alpha_pi_v = (1-alpha_c)/K; % power allocation factor for private message of each user
% power allocation coefficients must satisfy constraint:
% alpha_c+alpha_pi_v*K <= 1
% note: the above assumes all legit users have same allocation factor

P_r = 2; % number of propagation paths arriving at the r-th RIS element
P_r_j = 3; % number of propagation paths departing from the r-th RIS element
P_e = 4; % number of propagation paths departing from the LEO satellite

% location coordinates
UAV_altitude = 300;
LEO_altitude = 150e3;
S_xyz = [50; LEO_altitude-UAV_altitude; 200]; % location of LEO satellite
R_xyz = [0; 0; 0]; % location of STAR-RIS; code assumes this to be origin;
% note: this code assumes surface is on x-y plane (surface normal points in
% z-axis direction)

user_lower_bound = [-100;-UAV_altitude;-100]; % lower boundary for users in the x, y and z axes
user_upper_bound = [100;-UAV_altitude;100]; % upper boundary for users in the x, y and z axes
% note: user_lower_bound(1)==user_upper_bound(1)==(-UAV_altitude) means
% users are always on the ground

eve_lower_bound = [-20;-5;-20]; % boundary for eve in the x direction
eve_upper_bound = [20;5;20]; % boundary for eve in the y direction
% note:assumes UAV_altitude > -eve_lower_bound(1) to ensure eve altitude
% isn't negative (below ground)

% randomly generate user locations with uniform distribution with the above
% boundary ranges
rho_j_xyz = rand(3,K+L,num_samples); % location of each receiver (user and eavesdropper)
rho_j_xyz(:,1:K,:) = rho_j_xyz(:,1:K,:) .* (user_upper_bound-user_lower_bound) + user_lower_bound;
rho_j_xyz(:,K+1:end,:) = rho_j_xyz(:,K+1:end,:) .* (eve_upper_bound-eve_lower_bound) + eve_lower_bound;
% note: the first K elements of rho_j_xyz are for legit users and the next
% L elements are for eavesdroppers

% find out whether each receiver is on the reflect side or transmit side
reflect = sign(rho_j_xyz(3,:,:)) == sign(S_xyz(3)); % shape (1,K+L,num_samples)
transmit = 1 - reflect;  % shape (1,K+L,num_samples)
% how it works: because RIS surface normal is (0,0,1) which points in the z
% direction, receiver is on the reflecting side if the polarity of its z
% coordinate equals the polarity of the LEO's z coordinate. If polarity is
% opposite, then the receiver is on the transmit side.

% Nakagami Parameters:

% from LEO satellite to STAR-RIS (one value for each path and assume same
% for each RIS element):
m_q = 1*ones(1,P_r); % shape parameter
omega_q = 10*ones(1,P_r); % spread parameter

% from STAR-RIS to users and eavesdroppers (one value for each path and 
% each receiver, and assume same for each RIS element):
m_j_q = 0.5*ones(K+L,P_r_j); % shape parameter
omega_j_q = 1*ones(K+L,P_r_j); % spread parameter
% note: m_j_1(1:K,:) and omega_j_q(1:K,:) are for legit users
% and m_j_1(K+1:end,:) and omega_j_q(K+1:end,:) are for eavesdroppers
% direct channel from LEO satellite to Eve (one value for each path and
% each eavesdropper)

m_g = 1*ones(L,P_e); % shape parameter
omega_g = 10*ones(L,P_e); % spread parameter
% note: if we want different paths to have different parameters, or legit
% users and eavesdroppers to have different parameters, modify the above
% matrices correspondingly, but ensure the shape doesn't change

% Parameters for Correlation Matrix:
kappa = 0.4; % exponential correlation coefficient

% derived gamma distribution scales
theta_q = omega_q ./ m_q; % from LEO satellite to STAR-RIS
theta_j_q = omega_j_q ./ m_j_q; % from STAR-RIS to legit users

% Use vectorized indexing to build |i-j| matrix
i = (1:N_R)'; j = 1:N_R;
R = kappa .^ abs(i - j);      % Exponential correlation: R(i,j) = epsilon^|i-j|

% Ensure symmetry and unit diagonal
R = (R + R') / 2;               % force symmetry
R(1:N_R+1:end) = 1;               % set diagonal elements to 1

% Enforce positive definiteness if needed
minEig = min(eig(R));
if minEig <= 0
    warning('Adjusting correlation matrix to nearest SPD (min eig = %.3e)...', minEig);
    R = nearestSPD(R);
    % normalize diagonals
    D = sqrt(diag(R));
    R = R ./ (D * D');
end

% Generate Correlated Samples via Gaussian Copula
H_r_p_2 = zeros(num_samples, N_R, P_r);
H_r_j_q_2 = zeros(num_samples, N_R, P_r_j, K+L);
% Generate correlated uniforms for every path
for path = 1:P_r
    U = copularnd('Gaussian', R, num_samples); % LEO to STAR-RIS
    % Transform uniforms to gamma variates (|h|^2); (nakagami-m variables squared)
    H_r_p_2(:,:,path) = gaminv(U, m_q(path), theta_q(path));
end
for path = 1:P_r_j
    for user = 1:K
        U = copularnd('Gaussian', R, num_samples); % STAR-RIS to the legit user
        % Transform uniforms to gamma variates (|h|^2); (nakagami-m variables squared)
        H_r_j_q_2(:,:,path,user) = gaminv(U, m_j_q(user,path), theta_j_q(user,path));
    end
    for eve = K+1:K+L
        U = copularnd('Gaussian', R, num_samples); % STAR-RIS to the eavesdropper
        % Transform uniforms to gamma variates (|h|^2); (nakagami-m variables squared)
        H_r_j_q_2(:,:,path,eve) = gaminv(U, m_j_q(eve,path), theta_j_q(eve,path));
    end
end
% Form Nakagami-m variable magnitudes
h_r_p = sqrt(H_r_p_2); % magnitudes for channels from LEO to STAR-RIS
h_r_j_q = sqrt(H_r_j_q_2); % magnitudes for channels from STAR-RIS to users and eavesdroppers

% Add uniform random phase
h_r_p = h_r_p.*exp(1j.*rand(size(h_r_p)).*2.*pi);
h_r_j_q = h_r_j_q.*exp(1j.*rand(size(h_r_j_q)).*2.*pi);

% Form Nakagami-m variable magnitudes and uniform random phase for direct
% path from LEO to Eve
g_e_u = zeros(num_samples, P_e);
for path = 1:P_e
    nakagami_dist = makedist('Nakagami','mu',m_g(path),'omega',omega_g(path));
    % direct channel from LEO satellite to Eve (no spatial correlation)
    g_e_u(:,path) = random(nakagami_dist,num_samples,1).*exp(1j.*rand(num_samples,1).*2.*pi);
end

% STAR-RIS phase and magnitude of each RIS element
zeta_k_Sr = rand(num_samples, N_R); % reflection coefficients
zeta_k_St = 1 - zeta_k_Sr; % transmit coefficients
% each receiver either takes the reflection coefficients or the transmit
% coefficients depending on which side it is on
% reminder: boolean matrix 'reflect' have shape (1,K+L,num_samples) and 
% will have value of 1 if user is on reflect side and 0 if not. Vice versa
% for 'transmit' matrix. By multiplying each matrix by its respective
% coefficient zeta, and adding element-wise, we get a matrix where
% reflecting users get the reflecting coefficient and transmitting user
% gets the transmitting coefficient
beta_r = (sqrt(zeta_k_Sr).*permute(reflect,[3,1,2])+sqrt(zeta_k_St).*permute(transmit,[3,1,2])).*exp(1j.*rand(num_samples,N_R,K+L).*2.*pi);
% note: beta_r will have shape (num_samples,N_R,K+L)
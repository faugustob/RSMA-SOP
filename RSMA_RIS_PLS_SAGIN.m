clear;

% Simulation parameters
M = 1e2; % number of samples
N_path = 2; % number of paths

% Parameters

Light_speed = physconst('LightSpeed');
freq=10.5e9;
lambda = Light_speed/freq; % Compute the wavelength based on the frequency
half_wave_length = lambda/2; % Compute the half_wave_length

% location coordinates

UAV_altitude = 300;
LEO_altitude = 150e3;

w_R = [0, 0, 0]; % location of STAR-RIS; code assumes this to be origin;
% note: surface is on x-y plane; surface normal points in z-axis direction
w_U_r = [-50, -UAV_altitude,50]; % location of U_r; z value should be same polarity as z value of satellite (reflecting user)
w_U_t = [50, -UAV_altitude, -50]; % location of U_t; z value should be opposite polarity as z value of satellite (transmitting user)
w_S = [50, LEO_altitude-UAV_altitude, 200]; % location of LEO satellite
w_E = [10, -10, 20]; % location of Eve (assuming UAV_altitude > 10); on the reflecting side if z value is positive

%RIS_geometry

Ny = 50;
Nx = 34;
N = Nx * Ny;

RIS_element_size = floor(half_wave_length * 1000) / 1000;
padding_x = RIS_element_size + RIS_element_size * 0.2;
padding_y = padding_x;

RIS_block_size_x = 2 * padding_x + RIS_element_size;
RIS_block_size_y = RIS_block_size_x;

RIS_size_x = RIS_block_size_x * Nx;
RIS_size_y = RIS_block_size_y * Ny;

rn = cell(Ny, Nx);

r_S = cell(Ny, Nx);
r_U_r = cell(Ny, Nx);
r_U_t = cell(Ny, Nx);
r_E = cell(Ny, Nx);

r_t_nm_S = zeros(Ny, Nx);
r_r_nm_U_r = zeros(Ny, Nx);
r_r_nm_U_t = zeros(Ny, Nx);
r_r_nm_E = zeros(Ny, Nx);

theta_t_nm_S = zeros(Ny, Nx);
theta_r_nm_U_r = zeros(Ny, Nx);
theta_r_nm_U_t = zeros(Ny, Nx);
theta_r_nm_E = zeros(Ny, Nx);

phi_t_nm_S = zeros(Ny, Nx);
phi_r_nm_U_r = zeros(Ny, Nx);
phi_r_nm_U_t = zeros(Ny, Nx);
phi_r_nm_E = zeros(Ny, Nx);

theta_tx_nm_S = zeros(Ny, Nx);
theta_rx_nm_U_r = zeros(Ny, Nx);
theta_rx_nm_U_t = zeros(Ny, Nx);
theta_rx_nm_E = zeros(Ny, Nx);

phi_tx_nm_S = zeros(Ny, Nx);
phi_rx_nm_U_r = zeros(Ny, Nx);
phi_rx_nm_U_t = zeros(Ny, Nx);
phi_rx_nm_E = zeros(Ny, Nx);

distances = zeros(Ny, Nx);

% finding the center positions, as well as distances to
% TX/RX, elevation and azimuth angles for each RIS element
% note: Ny and Nx are assumed to be even numbers
for j = 1:Ny
    for i = 1:Nx
        x = -RIS_size_x / 2 + (i - 0.5) * RIS_block_size_x;
        y = -RIS_size_y / 2 + (j - 0.5) * RIS_block_size_y;
        z = 0;
        % Note: the normal of the RIS surface is in the z direction
        
        % Store position vector
        rn{j, i} = [x, y, z];
        % Compute distance relative to RIS center
        distances(j,i) = sqrt(x^2 + y^2); % or norm(rn{j,i})

        % Vector from satellite to RIS element
        r_S{j, i} = -w_S + rn{j, i};
        % Vector from RIS element to reflecting user
        r_U_r{j, i} = -rn{j, i} + w_U_r;
        % Vector from RIS element to transmitting user
        r_U_t{j, i} = -rn{j, i} + w_U_t;
        % Vector from RIS element to eavesdropper
        r_E{j, i} = -rn{j, i} + w_E;

        % Distance from satellite to RIS element
        r_t_nm_S(j, i) = norm(r_S{j, i});
        % Distance from RIS element to reflecting user
        r_r_nm_U_r(j, i) = norm(r_U_r{j, i});
        % Distance from RIS element to transmitting user
        r_r_nm_U_t(j, i) = norm(r_U_t{j, i});
        % Distance from RIS element to eavesdropper
        r_r_nm_E(j, i) = norm(r_E{j, i});

        % Elevation angle of satellite with RIS element (angle between vector from satellite to RIS element and RIS surface normal)
        theta_t_nm_S(j, i) = acos(norm(dot(r_S{j, i},[0, 0, 1]))/(norm(r_S{j, i})*norm([0, 0, 1])));
        % Elevation angle of reflecting user with RIS element (angle between vector from RIS element to reflecting user and RIS surface normal)
        theta_r_nm_U_r(j, i) = acos(norm(dot(r_U_r{j, i},[0, 0, 1]))/(norm(r_U_r{j, i})*norm([0, 0, 1])));
        % Elevation angle of transmitting user with RIS element (angle between vector from RIS element to transmitting user and RIS surface normal)
        theta_r_nm_U_t(j, i) = acos(norm(dot(r_U_t{j, i},[0, 0, 1]))/(norm(r_U_t{j, i})*norm([0, 0, 1])));        
        % Elevation angle of eavesdropper with RIS element (angle between vector from RIS element to eavesdropper and RIS surface normal)
        theta_r_nm_E(j, i) = acos(norm(dot(r_E{j, i},[0, 0, 1]))/(norm(r_E{j, i})*norm([0, 0, 1])));

        % Azimuth angle of satellite with RIS element (angle between vector from satellite to RIS element projected on x-y plane and x-axis)
        phi_t_nm_S(j, i) = acos(norm(dot([r_S{j, i}(1:2), 0],[1, 0, 0]))/(norm([r_S{j, i}(1:2), 0])*norm([1, 0, 0])));
        % Azimuth angle of reflecting user with RIS element (angle between vector from RIS element to reflecting user projected on x-y plane and x-axis)
        phi_r_nm_U_r(j, i) = acos(norm(dot([r_U_r{j, i}(1:2), 0],[1, 0, 0]))/(norm([r_U_r{j, i}(1:2), 0])*norm([1, 0, 0])));   
        % Azimuth angle of transmitting user with RIS element (angle between vector from RIS element to transmitting user projected on x-y plane and x-axis)
        phi_r_nm_U_t(j, i) = acos(norm(dot([r_U_t{j, i}(1:2), 0],[1, 0, 0]))/(norm([r_U_t{j, i}(1:2), 0])*norm([1, 0, 0])));
        % Azimuth angle of eavesdropper with RIS element (angle between vector from RIS element to eavesdropper projected on x-y plane and x-axis)
        phi_r_nm_E(j, i) = acos(norm(dot([r_E{j, i}(1:2), 0],[1, 0, 0]))/(norm([r_E{j, i}(1:2), 0])*norm([1, 0, 0])));    

        % angle between vector from satellite to RIS element and satellite to RIS center
        theta_tx_nm_S(j, i) = acos(norm(dot(r_S{j, i},[0, 0, 1]))/(norm(r_S{j, i})*norm([0, 0, 1])));
        % angle between vector from RIS element to reflecting user and RIS center to reflecting user
        theta_rx_nm_U_r(j, i) = acos(norm(dot(r_U_r{j, i},[0, 0, 1]))/(norm(r_U_r{j, i})*norm([0, 0, 1])));
        % angle between vector from RIS element to transmitting user and RIS center to transmitting user
        theta_rx_nm_U_t(j, i) = acos(norm(dot(r_U_t{j, i},[0, 0, 1]))/(norm(r_U_t{j, i})*norm([0, 0, 1])));        
        % angle between vector from RIS element to eavesdropper and RIS center to eavesdropper
        theta_rx_nm_E(j, i) = acos(norm(dot(r_E{j, i},[0, 0, 1]))/(norm(r_E{j, i})*norm([0, 0, 1])));
    end
end

% Antenna angles for each TX/RX
% Elevation angle of satellite with RIS center (angle between vector from satellite to RIS center and RIS surface normal)
theta_t_S = acos(norm(dot(w_S,[0, 0, 1]))/(norm(w_S)*norm([0, 0, 1])));
% Elevation angle of reflecting user with RIS center (angle between vector from RIS center to reflecting user and RIS surface normal)
theta_r_U_r = acos(norm(dot(w_U_r,[0, 0, 1]))/(norm(w_U_r)*norm([0, 0, 1])));
% Elevation angle of transmitting user with RIS center (angle between vector from RIS center to transmitting user and RIS surface normal)
theta_r_U_t = acos(norm(dot(w_U_t,[0, 0, 1]))/(norm(w_U_t)*norm([0, 0, 1])));        
% Elevation angle of eavesdropper with RIS center (angle between vector from RIS center to eavesdropper and RIS surface normal)
theta_r_E = acos(norm(dot(w_E,[0, 0, 1]))/(norm(w_E)*norm([0, 0, 1])));
% Azimuth angle of satellite with RIS center (angle between vector from satellite to RIS center projected on x-y plane and x-axis)
phi_t_S = acos(norm(dot([w_S(1:2), 0],[1, 0, 0]))/(norm([w_S(1:2), 0])*norm([1, 0, 0])));
% Azimuth angle of reflecting user with RIS center (angle between vector from RIS center to reflecting user projected on x-y plane and x-axis)
phi_r_U_r = acos(norm(dot([w_U_r(1:2), 0],[1, 0, 0]))/(norm([w_U_r(1:2), 0])*norm([1, 0, 0])));   
% Azimuth angle of transmitting user with RIS center (angle between vector from RIS center to transmitting user projected on x-y plane and x-axis)
phi_r_U_t = acos(norm(dot([w_U_t(1:2), 0],[1, 0, 0]))/(norm([w_U_t(1:2), 0])*norm([1, 0, 0])));
% Azimuth angle of eavesdropper with RIS center (angle between vector from RIS center to eavesdropper projected on x-y plane and x-axis)
phi_r_E = acos(norm(dot([w_E(1:2), 0],[1, 0, 0]))/(norm([w_E(1:2), 0])*norm([1, 0, 0])));  

% The paper has one proposition for specular reflection (incident angle =
% reflection angle) and one proposition for intelligent reflection (the
% reflected signal can be in any desired direction.

% According to the paper, the direction of the reflected signal is
% determined by the reflection coefficient (phase shift). I don't 
% understand how this will work with our paper. The paper calls this
% intelligent reflection.

% note: near-field case considers each RIS element individually while the
% far-field case doesn't

% Additional path loss parameters from W. Tang paper
near_field_threshold = 2*N*RIS_element_size^2/lambda;
F = @(theta,phi) cos(theta).^3; % normalized power radiation pattern function for RIS unit cell (fixed once RIS is designed and fabricated)
F_tx = @(theta,phi) cos(theta).^62; % normalized power radiation pattern function for transmitter antenna
F_rx_reflect = @(theta,phi) cos(theta).^62; % normalized power radiation pattern function for reflecting receiver antenna
F_rx_transmit = @(theta,phi) cos(theta).^62; % normalized power radiation pattern function for transmitting receiver antenna
F_rx_eve = @(theta,phi) cos(theta).^62; % normalized power radiation pattern function for eavesdropper antenna
% The antenna gains and unit cell gains depend on the normalized power
% radiation pattern function F().
% note: theta should be integrated over [0,pi] but it is [0,pi/2] in the
% below expressions because our F() function equals 0 for theta ~ (pi/2,pi]
G_unit_cell = 4*pi / integral2(@(theta,phi) F(theta,phi) .* sin(theta), 0, pi/2, 0, 2*pi); % Gain of unit cell on RIS surface
G_tx = 4*pi / integral2(@(theta,phi) F_tx(theta,phi) .* sin(theta), 0, pi/2, 0, 2*pi); % Gain of the transmitter antenna
G_rx_reflect = 4*pi / integral2(@(theta,phi) F_rx_reflect(theta,phi) .* sin(theta), 0, pi/2, 0, 2*pi); % Gain of the reflecting receiver antenna
G_rx_transmit = 4*pi / integral2(@(theta,phi) F_rx_transmit(theta,phi) .* sin(theta), 0, pi/2, 0, 2*pi); % Gain of the transmitting receiver antenna
G_rx_eve = 4*pi / integral2(@(theta,phi) F_rx_transmit(theta,phi) .* sin(theta), 0, pi/2, 0, 2*pi); % Gain of the eavesdropper antenna

% distances
d_t_c = norm(w_S-w_R); % distance from LEO to STAR-RIS
d_r_c_r = norm(w_U_r-w_R); % distance from STAR-RIS to reflect UE
d_r_c_t = norm(w_U_t-w_R); % distance from STAR-RIS to transmit UE
d_r_c_e = norm(w_E-w_R); % distance from STAR-RIS to Eve
d_SE = norm(w_S-w_E); % distance from LEO to Eve;
% note: paper assumes distance between TX/RX and RIS elements to be larger
% than 5*lambda, which is regarded as the lowerbound of near-field

% RIS elements' amplitude and phase coefficients
A_r = 0.9; % reflection coefficients
A_t = 0.9; % transmission (refraction) coefficients
% note: amplitude A_r and A_t assumed to be constant between RIS elements
uniform_dist = makedist('Uniform','lower',0,'upper',2*pi);
ris_angle_r = random(uniform_dist,Ny,Nx,M); % angle of RIS elements for reflection
ris_angle_t = random(uniform_dist,Ny,Nx,M); % angle of RIS elements for transmission

% note: path loss model assumes constant amplitude and different phases accross RIS elements

% path loss experienced by the signal reaching reflect UE
if d_t_c < near_field_threshold || d_r_c_r < near_field_threshold
    P_L_r = near_field_PL(G_tx,G_rx_reflect,G_unit_cell,RIS_element_size,RIS_element_size,lambda,A_r,F_tx,F,F_rx_reflect,theta_tx_nm_S,phi_tx_nm_S,theta_t_nm_S,phi_t_nm_S,theta_r_nm_U_r,phi_r_nm_U_r,theta_rx_nm_U_r,phi_rx_nm_U_r,r_t_nm_S,r_r_nm_U_r,ris_angle_r);
else
    P_L_r = far_field_PL(G_tx,G_rx_reflect,G_unit_cell,RIS_element_size,RIS_element_size,lambda,F,theta_t_S,phi_t_S,theta_r_U_r,phi_r_U_r,A_r,d_t_c,d_r_c_r,Ny,Nx,ris_angle_r);
end
% path loss experienced by the signal reaching transmit UE
if d_t_c < near_field_threshold || d_r_c_t < near_field_threshold
    P_L_t = near_field_PL(G_tx,G_rx_transmit,G_unit_cell,RIS_element_size,RIS_element_size,lambda,A_t,F_tx,F,F_rx_transmit,theta_tx_nm_S,phi_tx_nm_S,theta_t_nm_S,phi_t_nm_S,theta_r_nm_U_t,phi_r_nm_U_t,theta_rx_nm_U_t,phi_rx_nm_U_t,r_t_nm_S,r_r_nm_U_t,ris_angle_t);
else
    P_L_t = far_field_PL(G_tx,G_rx_transmit,G_unit_cell,RIS_element_size,RIS_element_size,lambda,F,theta_t_S,phi_t_S,theta_r_U_t,phi_r_U_t,A_t,d_t_c,d_r_c_t,Ny,Nx,ris_angle_t);
end
% path loss experienced by the signal reaching Eve, direct path from satellite
P_L_E_1 = (lambda/(4*pi))^4*G_tx*G_rx_eve/(d_SE^2); % still using the old path loss model
% path loss experienced by the signal reaching Eve, path via STAR-RIS
if d_t_c < near_field_threshold || d_r_c_e < near_field_threshold
    P_L_E_2 = near_field_PL(G_tx,G_rx_eve,G_unit_cell,RIS_element_size,RIS_element_size,lambda,A_r,F_tx,F,F_rx_eve,theta_tx_nm_S,phi_tx_nm_S,theta_t_nm_S,phi_t_nm_S,theta_r_nm_E,phi_r_nm_E,theta_rx_nm_E,phi_rx_nm_E,r_t_nm_S,r_r_nm_E,ris_angle_r);
else
    P_L_E_2 = far_field_PL(G_tx,G_rx_eve,G_unit_cell,RIS_element_size,RIS_element_size,lambda,F,theta_t_S,phi_t_S,theta_r_E,phi_r_E,A_r,d_t_c,d_r_c_e,Ny,Nx,ris_angle_r);
end

% note: in calculating the path loss, individual RIS element parameters are
% in (Ny,Nx,M) tensors, but in the code below, they are in (N,M) matrices
% (where Ny*Nx=N)

sigma_r = 1; % standard deviation for gaussian noise at reflect UE
sigma_t = 1; % standard deviation for gaussian noise at transmit UE
sigma_E = 1; % standard deviation for gaussian noise at Eve

% RSMA parameters
alpha_c = 2/3; % power allocation factor for common message
alpha_p_r = 1/6; % power allocation factor for private message of reflecting user
alpha_p_t = 1/6; % power allocation factor for private message of transmitting user
eta = 0; % error factor associated with imperfect SIC; not used yet

% nakagami parameters
m_SR = 1; % shape parameter; from LEO satellite to STAR-RIS
omega_SR = 10; % spread parameter; from LEO satellite to STAR-RIS
m_r = 0.5; % shape parameter; from STAR-RIS to U_r
omega_r = 1; % spread parameter; from STAR-RIS to U_r
m_t = 0.5; % shape parameter; from STAR-RIS to U_t
omega_t = 1; % spread parameter; from STAR-RIS to U_t
m_E = 0.5; % shape parameter; from STAR-RIS to Eve
omega_E = 1; % spread parameter; from STAR-RIS to Eve
m_SE = 1; % shape parameter; direct channel from LEO satellite to Eve
omega_SE = 10; % spread parameter; direct channel from LEO satellite to Eve

% parameters for correlation matrix
kappa = 0.4; % exponential correlation coefficient
% derived gamma distribution scales
theta_SR = omega_SR / m_SR;
theta_r = omega_r / m_r;
theta_t = omega_t / m_t;
theta_E = omega_E / m_E;

% Use vectorized indexing to build |i-j| matrix
i = (1:N)'; j = 1:N;
R = kappa .^ abs(i - j);      % Exponential correlation: R(i,j) = epsilon^|i-j|
% if this is spacial correlation, we need to consider the RIS as a Nx by Ny
% matrix rather than a flattened vector???

% Ensure symmetry and unit diagonal
R = (R + R') / 2;               % force symmetry
R(1:N+1:end) = 1;               % set diagonal elements to 1

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
% Generate correlated uniforms without loops
U_h = copularnd('Gaussian', R, M); % LEO to i-th STAR-RIS element
U_h_r = copularnd('Gaussian', R, M); % i-th RIS element to reflect UE
U_h_t = copularnd('Gaussian', R, M); % i-th RIS element to transmit UE
U_g = copularnd('Gaussian', R, M); % i-th STAR-RIS element to Eve

% Transform uniforms to gamma variates (|h|^2 and |g|^2); (nakagami-m
% variables squared)
H2 = gaminv(U_h, m_SR, theta_SR);   % size: [M x N]
H_r2 = gaminv(U_h_r, m_r, theta_r);
H_t2 = gaminv(U_h_t, m_t, theta_t);
G2 = gaminv(U_g, m_E, theta_E);

% Element-wise RIS magnitude and phase shift for each element and sample
zeta_r_k = rand(N,M); % reflection coefficients
zeta_t_k = 1-zeta_r_k; % transmission (refraction) coefficients
phi_r = sqrt(zeta_r_k).*exp(1j.*reshape(ris_angle_r,N,M)); % [N x M]
phi_t = sqrt(zeta_t_k).*exp(1j.*reshape(ris_angle_t,N,M)); % [N x M]
% note: (N,M) is (number of RIS elements,number of samples). The path loss
% model above considers (Ny,Nx,M) tensors instead of (N,M) matrices.

% Form Nakagami-m Variables
H = sqrt(H2); % magnitude for channel from LEO to i-th STAR-RIS element
H_r = sqrt(H_r2); % magnitude for channel from i-th RIS element to reflect UE
H_t = sqrt(H_t2); % magnitude for channel from i-th RIS element to transmit UE
G = sqrt(G2); % magnitude for channel from i-th STAR-RIS element to Eve

nakagami_dist_SE = makedist('Nakagami','mu',m_SE,'omega',omega_SE);

% direct channel from LEO satellite to Eve (no correlation matrix because
% no RIS elements; size: [1 x M])
g_1 = random(nakagami_dist_SE,1,M).*exp(1j.*random(uniform_dist,1,M));

% add random uniform phase (0 to 2*pi)
h = H'.*exp(1j.*random(uniform_dist,N,M)); % size: [N x M]
h_r = H_r'.*exp(1j.*random(uniform_dist,N,M));
h_t = H_t'.*exp(1j.*random(uniform_dist,N,M));
g_2 = G'.*exp(1j.*random(uniform_dist,N,M));

% vary the transmit power P_S and plot SINR vs transmit power
P_S_dB_vector = 100:5:400;%dB; range from 100 to 400 dB
P_S_vector = 10.^(P_S_dB_vector./10); % transmit power, in watts (P_S)

gamma_c_r_vector = zeros(size(P_S_vector));
gamma_c_t_vector = zeros(size(P_S_vector));
gamma_p_r_vector = zeros(size(P_S_vector));
gamma_p_t_vector = zeros(size(P_S_vector));
gamma_c_E_vector = zeros(size(P_S_vector));
gamma_p_r_E_vector = zeros(size(P_S_vector));
gamma_p_t_E_vector = zeros(size(P_S_vector));

iteration = 1;
for P_S = P_S_vector
    gamma_r_bar = P_S/sigma_r^2;
    gamma_t_bar = P_S/sigma_t^2;
    % channel coefficients and RIS coefficients and phases for users
    channel_r = abs(sum(conj(h_r).*phi_r.*h,1)).^2;
    
    channel_t = abs(sum(conj(h_t).*phi_t.*h,1)).^2;
    % SINR of the common signal part at the reflecting user
    gamma_c_r = alpha_c.*P_L_r.*channel_r.*gamma_r_bar./(gamma_r_bar.*channel_r.*P_L_r.*(alpha_p_r+alpha_p_t)+1);
    % SINR of the common signal part at the transmitting user
    gamma_c_t = alpha_c.*P_L_t.*channel_t.*gamma_t_bar./(gamma_t_bar.*channel_t.*P_L_t.*(alpha_p_r+alpha_p_t)+1);
    
    % SINR of the private signal part at the reflecting user
    gamma_p_r = alpha_p_r.*P_L_r.*channel_r.*gamma_r_bar./(gamma_r_bar.*channel_r.*P_L_r.*(alpha_p_t+eta.*alpha_c)+1);
    % SINR of the private signal part at the transmitting user
    gamma_p_t = alpha_p_t.*P_L_t.*channel_t.*gamma_t_bar./(gamma_t_bar.*channel_t.*P_L_t.*(alpha_p_r+eta.*alpha_c)+1);
    
    gamma_E_bar = P_S/sigma_E^2;
    % channel coefficients and RIS coefficients and phases, including path loss, for Eve
    channel_E = abs(sum(sqrt(P_L_E_2).*conj(g_2).*phi_r.*h,1)+sqrt(P_L_E_1).*g_1).^2;
    
    % SINR of the common signal part at Eve
    gamma_c_E = alpha_c.*channel_E.*gamma_E_bar./(gamma_E_bar.*channel_E.*(alpha_p_r+alpha_p_t)+1);
    
    % SINR of the reflecting user's private signal part at Eve
    gamma_p_r_E = alpha_p_r.*channel_E.*gamma_E_bar./(gamma_E_bar.*channel_E.*(alpha_p_r+alpha_p_t)+1);
    % SINR of the transmitting user's private signal part at Eve
    gamma_p_t_E = alpha_p_t.*channel_E.*gamma_E_bar./(gamma_E_bar.*channel_E.*(alpha_p_r+alpha_p_t)+1);
    
    C_c_r = max(log2(1+gamma_c_r)-log2(1+gamma_c_E), 0); % secrecy capacity of common message vs reflecting user?
    C_p_r = max(log2(1+gamma_p_r)-log2(1+gamma_p_r_E), 0); % secrecy capacity for private message of reflecting user?
    
    C_c_t = max(log2(1+gamma_c_t)-log2(1+gamma_c_E), 0); % secrecy capacity of common message vs transmitting user?
    C_p_t = max(log2(1+gamma_p_t)-log2(1+gamma_p_t_E), 0); % secrecy capacity for private message of transmitting user?
    
    % secrecy capacity thresholds
    R_c_r = 1e-3; % for common signal vs reflecting user
    R_p_r = 1e-3; % for private signal vs reflecting user
    R_c_t = 1e-3; % for common signal vs transimitting user
    R_p_t = 1e-3; % for private signal vs transmitting user
    
    % secrecy outage only happens if BOTH common and private signal fall below
    % threshold (or does a secrecy outage happen if only one of the two fall
    % below threshold?)
    P_SOP_r = mean((C_c_r<R_c_r)&(C_p_r<R_p_r));
    P_SOP_t = mean((C_c_t<R_c_t)&(C_p_t<R_p_t));

    gamma_c_r_vector(iteration) = mean(gamma_c_r);
    gamma_c_t_vector(iteration) = mean(gamma_c_t);
    gamma_p_r_vector(iteration) = mean(gamma_p_r);
    gamma_p_t_vector(iteration) = mean(gamma_p_t);
    gamma_c_E_vector(iteration) = mean(gamma_c_E);
    gamma_p_r_E_vector(iteration) = mean(gamma_p_r_E);
    gamma_p_t_E_vector(iteration) = mean(gamma_p_t_E);

    iteration = iteration + 1;
end

figure
plot(P_S_dB_vector, gamma_c_r_vector)
title('gamma_c_r'); 
figure
plot(P_S_dB_vector, gamma_c_t_vector)
title('gamma_c_t'); 
figure
plot(P_S_dB_vector, gamma_p_r_vector)
title('gamma_p_r'); 
figure
plot(P_S_dB_vector, gamma_p_t_vector)
title('gamma_p_t');
figure
plot(P_S_dB_vector, gamma_c_E_vector)
title('gamma_c_E'); 
figure
plot(P_S_dB_vector, gamma_p_r_E_vector)
title('gamma_p_r_E'); 
figure
plot(P_S_dB_vector, gamma_p_t_E_vector)
title('gamma_p_t_E');
% At legit user, SINR converges at 2 for common signal and 1 for private signal
% At Eve, SINR converges at 2 for common signal and 0.5 for private signal
% SINR converges at around 200dB for Eve but around 300dB for legit user
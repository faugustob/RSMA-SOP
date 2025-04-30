clear

% location coordinates
w_U_r = [0, 0, 0]; % location of U_r
w_U_t = [100, 0, 0]; % location of U_t
w_R = [50, 0, 100]; % location of STAR-RIS
w_S = [50, 100, 400000]; % location of LEO satellite
w_E = [100, 50, 100]; % location of Eve

% parameters
N = 100; % number of passive reflecting elements
P_S = 50; % transmit power, in watts
lambda = 0.1; % wavelength
G_S_dB = 30; % transmit gain, in dB
G_S = 10^(G_S_dB/10); % transmit gain, linear
G_U_r_dB = 160; % receiving gains for U_r, in dB
G_U_t_dB = 160; % receiving gains for U_t, in dB
G_U_r = 10^(G_U_r_dB/10); % receiving gains for U_r, linear
G_U_t = 10^(G_U_t_dB/10); % receiving gains for U_t, linear
G_E_dB = 160; % receiving gain at Eve, in dB
G_E = 10^(G_E_dB/10); % receiving gain at Eve, linear
epsilon_p = 1; % RIS efficiency
d_SR = norm(w_S-w_R); % distance from LEO to STAR-RIS
d_SE = norm(w_S-w_E); % distance from LEO to Eve
d_U_r = norm(w_R-w_U_r); % distance from STAR-RIS to reflect UE
d_U_t = norm(w_R-w_U_t); % distance from STAR-RIS to transmit UE
d_E = norm(w_R-w_E);
% path loss experienced by the signal reaching reflect/transmit UE
P_L_r = (lambda/(4*pi))^4*G_S*G_U_r/(d_SR^2*d_U_r^2)*epsilon_p;
P_L_t = (lambda/(4*pi))^4*G_S*G_U_t/(d_SR^2*d_U_t^2)*epsilon_p;
% path loss experienced by the signal reaching Eve, direct path from satellite
P_L_E_1 = (lambda/(4*pi))^4*G_S*G_E/(d_SE^2)*epsilon_p;
% path loss experienced by the signal reaching Eve, path via STAR-RIS
P_L_E_2 = (lambda/(4*pi))^4*G_S*G_E/(d_SR^2*d_E^2)*epsilon_p;
sigma_r = 1; % standard deviation for gaussian noise at reflect UE
sigma_t = 1; % standard deviation for gaussian noise at transmit UE
sigma_e = 1; % standard deviation for gaussian noise at Eve

M = 10000; % number of samples

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

% RSMA parameters
alpha_c = 1/3; % power allocation factor for common message
alpha_p_r = 1/3; % power allocation factor for private message of reflecting user
alpha_p_t = 1/3; % power allocation factor for private message of transmitting user
eta = 0; % error factor associated with imperfect SIC; not used yet

% spatial correlation
% not implemented yet

% generate random samples
uniform_dist = makedist('Uniform','lower',0,'upper',2*pi);
nakagami_dist_SR = makedist('Nakagami','mu',m_SR,'omega',omega_SR);
nakagami_dist_E = makedist('Nakagami','mu',m_E,'omega',omega_E);
nakagami_dist_SE = makedist('Nakagami','mu',m_SE,'omega',omega_SE);
% channel from LEO to i-th STAR-RIS element
h = random(nakagami_dist_SR,N,M).*exp(1j.*random(uniform_dist,N,M));
% direct channel from LEO satellite to Eve
g_1 = random(nakagami_dist_SE,1,M).*exp(1j.*random(uniform_dist,1,M));
% channel from i-th STAR-RIS element to Eve
g_2 = random(nakagami_dist_E,N,M).*exp(1j.*random(uniform_dist,N,M));
% channel from i-th RIS element to reflect UE
nakagami_dist_r = makedist('Nakagami','mu',m_r,'omega',omega_r);
h_r = random(nakagami_dist_r,N,M).*exp(1j.*random(uniform_dist,N,M));
% channel from i-th RIS element to transmit UE
nakagami_dist_t = makedist('Nakagami','mu',m_t,'omega',omega_t);
h_t = random(nakagami_dist_t,N,M).*exp(1j.*random(uniform_dist,N,M));

zeta_r_k = rand(N,M); % reflection coefficients
zeta_t_k = 1-zeta_r_k; % transmission (refraction) coefficients
ris_angle_r = random(uniform_dist,N,M); % angle of RIS elements for reflection
ris_angle_t = random(uniform_dist,N,M); % angle of RIS elements for transmission
ris_angle_E = random(uniform_dist,N,M); % angle of RIS elements for Eve
% Element-wise RIS magnitude and phase shift for each element and sample
phi_r = sqrt(zeta_r_k).*exp(1j.*ris_angle_r); % [N x M]
phi_t = sqrt(zeta_t_k).*exp(1j.*ris_angle_t); % [N x M]
phi_E = sqrt(zeta_r_k).*exp(1j.*ris_angle_E); % [N x M]; on the reflecting side so takes reflection coefficient

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
channel_E = abs(sum(sqrt(P_L_E_2).*conj(g_2).*phi_E.*h,1)+sqrt(P_L_E_1).*g_1).^2;

% SINR of the common signal part at Eve
gamma_c_E = alpha_c.*channel_E.*gamma_E_bar./(gamma_E_bar.*channel_E.*(alpha_p_r+alpha_p_t)+1);

% SINR of the reflecting user's private signal part at Eve
gamma_p_E_r = alpha_p_r.*channel_E.*gamma_E_bar./(gamma_E_bar.*channel_E.*(alpha_p_r+alpha_p_t)+1);
% SINR of the transmitting user's private signal part at Eve
gamma_p_E_t = alpha_p_t.*channel_E.*gamma_E_bar./(gamma_E_bar.*channel_E.*(alpha_p_r+alpha_p_t)+1);
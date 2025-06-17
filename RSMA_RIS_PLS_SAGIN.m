clear

% location coordinates
w_U_r = [0, 0, 0]; % location of U_r
w_U_t = [100, 0, 0]; % location of U_t
w_R = [50, -50, 100]; % location of STAR-RIS
w_S = [50, 100, 400000]; % location of LEO satellite
w_E = [100, 50, 100]; % location of Eve
% need to add a parameter for orientation of RIS surface. currently
% assuming fixed

% parameters
N = 100; % number of passive reflecting elements
P_S = 1000; % transmit power, in watts
lambda = 0.1; % wavelength

% additional path loss parameters from J. Jeong et al. paper
d_t_c = norm(w_S-w_R); % distance from LEO to STAR-RIS
d_r_c_r = norm(w_R-w_U_r); % distance from STAR-RIS to reflect UE
d_r_c_t = norm(w_R-w_U_t); % distance from STAR-RIS to transmit UE
d_r_c_e = norm(w_R-w_E); % distance from STAR-RIS to Eve
d_SE = norm(w_S-w_E); % distance from LEO to Eve; not part of J. Jeong et al. path loss model
L_P_E = 1; % phase-error loss due to quantization
L_S_R = 1; % specular reflection loss
r_t = 0.2; % transmitter antenna's effective radius
r_d = 0.2; % receiver antenna's effective radius
G_t_a = (r_t*2*pi/lambda)^2; % average Tx antenna gain
G_r_a_r = (r_d*2*pi/lambda)^2; % average Rx antenna gain - reflecting user
G_r_a_t = (r_d*2*pi/lambda)^2; % average Rx antenna gain - transmitting user
G_r_a_e = (r_d*2*pi/lambda)^2; % average Rx antenna gain - eavesdropper
A_p = 5; % % physical aperture (surface area of RIS I assume)
RIS_normal = [0, 1, 0]; % normal of the RIS surface

incident_wave = w_R - w_S;
reflected_wave_r = w_U_r - w_R;
reflected_wave_t = w_U_t - w_R;
reflected_wave_e = w_E - w_R;

theta_i_c = acos(norm(dot(incident_wave,RIS_normal))/(norm(RIS_normal)*norm(incident_wave))); % angle between incident wave and the normal of surface (wave points to center of surface)
theta_r_c_r = acos(norm(dot(reflected_wave_r,RIS_normal))/(norm(RIS_normal)*norm(reflected_wave_r))); % angle between reflected wave to reflecting user and the normal of surface (wave points to center of surface)
theta_r_c_t = acos(norm(dot(reflected_wave_t,RIS_normal))/(norm(RIS_normal)*norm(reflected_wave_t))); % angle between reflected wave to reflecting user and the normal of surface (wave points to center of surface)
theta_r_c_e = acos(norm(dot(reflected_wave_e,RIS_normal))/(norm(RIS_normal)*norm(reflected_wave_e))); % angle between reflected wave to reflecting user and the normal of surface (wave points to center of surface)

% path loss experienced by the signal reaching reflect/transmit UE
P_L_r = N^2*L_P_E*L_S_R*cos(theta_i_c)*cos(theta_r_c_r)*G_t_a*G_r_a_r*A_p^2/(d_t_c^2*d_r_c_r^2*16*pi^2);
P_L_t = N^2*L_P_E*L_S_R*cos(theta_i_c)*cos(theta_r_c_t)*G_t_a*G_r_a_t*A_p^2/(d_t_c^2*d_r_c_t^2*16*pi^2);
% path loss experienced by the signal reaching Eve, direct path from satellite
P_L_E_1 = (lambda/(4*pi))^4*G_t_a*G_r_a_e/(d_SE^2); % still using the old path loss model
% path loss experienced by the signal reaching Eve, path via STAR-RIS
P_L_E_2 = N^2*L_P_E*L_S_R*cos(theta_i_c)*cos(theta_r_c_e)*G_t_a*G_r_a_e*A_p^2/(d_t_c^2*d_r_c_e^2*16*pi^2);

sigma_r = 1; % standard deviation for gaussian noise at reflect UE
sigma_t = 1; % standard deviation for gaussian noise at transmit UE
sigma_E = 1; % standard deviation for gaussian noise at Eve

M = 1e4; % number of samples

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

% Form Nakagami-m Variables
H = sqrt(H2); % magnitude for channel from LEO to i-th STAR-RIS element
H_r = sqrt(H_r2); % magnitude for channel from i-th RIS element to reflect UE
H_t = sqrt(H_t2); % magnitude for channel from i-th RIS element to transmit UE
G = sqrt(G2); % magnitude for channel from i-th STAR-RIS element to Eve

uniform_dist = makedist('Uniform','lower',0,'upper',2*pi);
nakagami_dist_SE = makedist('Nakagami','mu',m_SE,'omega',omega_SE);

% direct channel from LEO satellite to Eve (no correlation matrix because
% no RIS elements; size: [1 x M])
g_1 = random(nakagami_dist_SE,1,M).*exp(1j.*random(uniform_dist,1,M));

% add random uniform phase (0 to 2*pi)
h = H'.*exp(1j.*random(uniform_dist,N,M)); % size: [N x M]
h_r = H_r'.*exp(1j.*random(uniform_dist,N,M));
h_t = H_t'.*exp(1j.*random(uniform_dist,N,M));
g_2 = G'.*exp(1j.*random(uniform_dist,N,M));

zeta_r_k = rand(N,M); % reflection coefficients
zeta_t_k = 1-zeta_r_k; % transmission (refraction) coefficients
ris_angle_r = random(uniform_dist,N,M); % angle of RIS elements for reflection
ris_angle_t = random(uniform_dist,N,M); % angle of RIS elements for transmission
ris_angle_E = -angle(g_1)-angle(conj(g_2).*h); % angle of RIS elements for Eve (optimal for PLS?)
% notes: g_1 has size [1xM] and conj(g_2).*h has size [NxM];
% in channel calculation, conj(g_2).*phi_E.*h gets summed along the all the
% first dimension (N RIS elements), yet angle(sum(conj(g_2).*phi_E.*h,1))
% still equal -angle(g_1), which is good.

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

gamma_E_bar = 1; % P_S/sigma_E^2;
% channel coefficients and RIS coefficients and phases, including path loss, for Eve
channel_E = abs(sum(sqrt(P_L_E_2).*conj(g_2).*phi_E.*h,1)+sqrt(P_L_E_1).*g_1).^2;

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
R_c_r = 1; % for common signal vs reflecting user
R_p_r = 1; % for private signal vs reflecting user
R_c_t = 1; % for common signal vs transimitting user
R_p_t = 1; % for private signal vs transmitting user

% secrecy outage only happens if BOTH common and private signal fall below
% threshold (or does a secrecy outage happen if only one of the two fall
% below threshold?)
P_SOP_r = mean((C_c_r<R_c_r)&(C_p_r<R_p_r));
P_SOP_t = mean((C_c_t<R_c_t)&(C_p_t<R_p_t));
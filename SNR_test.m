%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monte-Carlo vs Closed-Form PDF for R = aX / (Z + b)
% X ~ Gamma(alpha, theta1), Z ~ Gamma(beta, theta2)
% Requires Symbolic Math Toolbox (for meijerG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%% ------------------- PARAMETERS -----------------------
alpha  = 2.5;       % shape of X
beta   = 3.2;       % shape of Z
theta1 = 1.0;       % scale of X
theta2 = 0.8;       % scale of Z
a      = 1.7;       % numerator multiplier
b      = 2.0;       % denominator shift (must be >0 for this code)

N = 2e6;            % number of Monte-Carlo samples

%% ------------------- MONTE-CARLO SIMULATION -----------

% MATLAB gamrnd uses shape/scale
X = gamrnd(alpha, theta1, [N, 1]);
Z = gamrnd(beta,  theta2, [N, 1]);

R_sim = a * X ./ (Z + b);

%% ------------------- THEORETICAL PDF (Meijer-G) --------

syms r positive

c = r/(a*theta1) + 1/theta2;

% Closed form PDF
f_r = ( r^(alpha-1) / (a^alpha * gamma(alpha) * gamma(beta) * theta1^alpha * theta2^beta) ) ...
      * exp( -b * r/(a*theta1) ) ...
      * ( b^(alpha+beta) / gamma(-alpha) ) ...
      * meijerG([1-beta], [], [0, -alpha-beta], [], b*c);

% Convert to numeric function handle
f_r_fun = matlabFunction(f_r, 'Vars', r);

%% ------------------- COMPARE HISTOGRAM & PDF ----------

figure; hold on;

% Histogram (normalized to pdf)
histogram(R_sim, 'Normalization', 'pdf', 'EdgeColor','none', 'FaceAlpha',0.35);

% Domain for plotting
r_vals = linspace(0, prctile(R_sim, 99.9), 600);

% Plot theoretical PDF
plot(r_vals, f_r_fun(r_vals), 'r', 'LineWidth', 2);

xlabel('r');
ylabel('Density');
title('Monte-Carlo vs Closed-Form PDF of R = aX/(Z+b)');
legend('Monte-Carlo Histogram','Closed-Form PDF (Meijer-G)');
grid on;


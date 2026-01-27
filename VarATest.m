clear; clc; rng(2);


%% ===================== Parameters =====================

Ns = 1e6; %Number of monte-carlo repetitions.

I = 1; % Indicator function
Pe = 2; Pr = 2; Pq = 1; 
m_e = 0.6; m_q = 1; m_p = 0.9; Omega_e = 1/Pe; Omega_p = 1/Pr; Omega_q = 1/Pq;

M = 2; N = 2; 
K=2; % number of legitimate users
L=1; % number of eavesdroppers

 alpha_array = rand(1, K);   % generate K random numbers in (0,1)
 alpha_array = alpha_array / sum(alpha_array);  % normalize so the sum is 1

Nsymb = M*N; 
Plos = 0.9; PLj = 0.9; Nr = 5;



%% ===================== Generate channel matrices and channels=====================
HB = randn(Nsymb,Nsymb,Pe) + 1i*randn(Nsymb,Nsymb,Pe); % Direct link
HA = randn(Nsymb,Nsymb,Pr,Pq) + 1i*randn(Nsymb,Nsymb,Pr,Pq); % Relay link

%% ===================== Allocate arrays =====================
Nc_opt_term1 = zeros(Ns,1);
Nc_opt_term2 = zeros(Ns,1);
Nc_opt_term3 = zeros(Ns,1);
Nc_opt= zeros(Ns,1);

Nc_trace = zeros(Ns,1);
Term1 = zeros(Ns,1);
Term2 = zeros(Ns,1);
Term3 = zeros(Ns,1);

g_pq = rand(Pr,Pq) + 1i*rand(Pr,Pq);
beta_opt = ones(Nr,1);
beta_non_opt = rand(Nr,1) + 1i*rand(Nr,1);

% % h_rp = sqrt(gamrnd(m_p, Omega_p/m_p, [Nr,Pr,Ns]));
% % h_jq = sqrt(gamrnd(m_q, Omega_q/m_q, [Nr,Pq,Ns]));
% % 
% % h_e = sqrt(gamrnd(m_e, Omega_e/m_e, [Pe,Ns]));
% 
%  for n = 1:Ns
% 
%     % ------------------ Construct Term1 ------------------
% 
%      term1 = zeros(Nsymb,Nsymb);
%     for r = 1:Nr
%        beta = beta_opt(r);
%         for p = 1:Pr
%             for q = 1:Pq
% 
%                 term1 = term1 + beta * h_rp(r,p,n) * h_jq(r,q,n) * g_pq(p,q) * HA(:,:,p,q);
%             end
%         end
%     end
%     term1 = sqrt(PLj)*term1;
% 
% 
%     % ------------------ Construct Term2 ------------------
%     B1 = zeros(Nsymb,Nsymb);
%     for i = 1:Pe
%         B1 = B1 + h_e(i,n) * HB(:,:,i);
%     end
%     term2 = I * sqrt(Plos) * B1;
% 
%     norm_term2 = norm(term2,'fro').^2;
% 
% 
%     % norm_term2_expanded = 0;
%     % for u=1:Pe
%     %     for up=1:Pe
%     %         norm_term2_expanded = norm_term2_expanded +  conj(h_e(u,n))*h_e(up,n)*trace(HB(:,:,u)'*HB(:,:,up));
%     %     end
%     % end
%     % 
%     % norm_term2_expanded = I * Plos * norm_term2_expanded;
% 
% 
%     % norm_Term3_exp = 0;
%     % 
%     % for r = 1: Nr
%     %     beta = beta_r(r);
%     %     for u=1:Pe
%     %         for p=1:Pr
%     %             for q=1:Pq
%     %                 norm_Term3_exp = norm_Term3_exp  + beta*h_e(u,n)* h_rp(r,p,n) * h_jq(r,q,n) * g_pq(p,q) * trace(HB(:,:,u)'*HA(:,:,p,q));
%     %             end
%     %         end
%     %     end
%     % end
%     % 
%     % norm_Term3_exp = 2*I*sqrt(Plos)*sqrt(PLj)*norm_Term3_exp;
%     norm_Term3 = 2*real(trace(term2'*term1));
%     % norm_Term3_diff = norm(term1+term2,'fro').^2-norm(term1,'fro').^2-norm_term2
% 
% 
% 
% 
%     Nc_opt_term1(n) = norm(term1,'fro').^2;
% 
%     Nc_opt_term2(n) = norm_term2;
% 
%     Nc_opt_term3(n) =  norm_Term3;
% 
%     Nc_opt(n) = norm(term1 + term2,'fro').^2;
% 
% 
%     % ------------------ Alternative normA_clean ------------------
%     % normA_clean = 0;
%     % for r=1:Nr
%     %     beta_s = conj(beta_r(r));
%     %     for rp=1:Nr
%     %          beta = beta_r(rp);
%     %         for p = 1:Pr
%     %             for q = 1:Pq
%     %                 for pp=1:Pr
%     %                     for qp =1:Pq
%     % 
%     %                         normA_clean = normA_clean+beta_s*beta*conj(h_rp(r,p,n))*h_rp(rp,pp,n)*conj(h_jq(r,q,n))*h_jq(rp,qp,n)...
%     %                         *conj(g_pq(p,q))*g_pq(pp,qp)*trace(HA(:,:,p,q)'*HA(:,:,pp,qp));
%     % 
%     %                     end
%     %                 end
%     %             end
%     %         end
%     %     end
%     % end
%     % normA_clean = real(normA_clean);
% 
%     %Nc_exp(n) = normA_clean;
%  end
% 
% 
% 
% mean_NC_term1 = mean(Nc_opt_term1)
% 
% 
% var_NC_term1 = var(Nc_opt_term1)
% 
% mean_NC_term2 = mean(Nc_opt_term2)
% 
% 
% 
% var_NC_term2 = var(Nc_opt_term2)
% 
% mean_NC_term3 = mean(Nc_opt_term3)
% 
% 
% var_NC_term3 = var(Nc_opt_term3)
% 
% mean_NC = mean(Nc_opt)
% 
% 
% var_NC = var(Nc_opt);

% [Dc,meanDMC,varDMC,meansD_MC,varsD_MC] = computeMCNonOptDMoments(I,Pe,Pr,Pq,m_e,m_q,m_p,Omega_e,Omega_p,Omega_q,Plos,PLj,Nr,HB,HA,g_pq,beta_non_opt,alpha_array,Ns,Nsymb);
% [meanD,varD,meansD,varsD] = computeNonOptDMoments(I,Pe,Pr,Pq,m_e,m_q,m_p,Omega_e,Omega_p,Omega_q,Plos,PLj,Nr,HB,HA,g_pq,beta_non_opt,alpha_array);


[Nc_opt,meanNMC_opt,varNMC_opt,meansMC_opt,varsMC_opt] = computeMCOptNMoments(I,Pe,Pr,Pq,m_e,m_q,m_p,Omega_e,Omega_p,Omega_q,Plos,PLj,Nr,HB,HA,g_pq,beta_opt,Ns,Nsymb);
[meanN_opt,varN_opt,mean_N_opt,vars_N_Opt,Covs] = computeOptNMoments(I,Pe,Pr,Pq,m_e,m_q,m_p,Omega_e,Omega_p,Omega_q,Plos,PLj,Nr,HB,HA,g_pq,beta_opt);
[Nc,meanNMC,varNMC,meansMC,varsMC] = computeMCNonOptNMoments(I,Pe,Pr,Pq,m_e,m_q,m_p,Omega_e,Omega_p,Omega_q,Plos,PLj,Nr,HB,HA,g_pq,beta_non_opt,Ns,Nsymb);

[meanN,varN,means,vars] = computeNonOptNMoments(I,Pe,Pr,Pq,m_e,m_q,m_p,Omega_e,Omega_p,Omega_q,Plos,PLj,Nr,HB,HA,g_pq,beta_non_opt);

%% =========================
x = linspace(1e-6, 600, 2000);
Nc = Nc(:);   % Monte-Carlo only for histogram

figure; hold on;

% ============================================================
% Non Opt
% ============================================================
histogram(Nc, 'Normalization', 'pdf', ...
          'FaceColor',[0.7 0.7 0.7], ...
          'EdgeColor','none', ...
          'DisplayName','Monte Carlo');

% ============================================================
% Closed-form moments
% ============================================================
mu = meanN;
v  = varN;

% ============================================================
% Gamma (MM)
% ============================================================
k_g     = mu^2 / v;
theta_g = v / mu;

g_pdf = @(x) gampdf(x, k_g, theta_g);

plot(x, g_pdf(x), 'r-', 'LineWidth',2, ...
     'DisplayName','Gamma (MM)');

% ============================================================
% Lognormal (MM)
% ============================================================
sigma2 = log(1 + v/mu^2);
sigma  = sqrt(sigma2);
mu_l   = log(mu) - 0.5*sigma2;

ln_pdf = @(x) lognpdf(x, mu_l, sigma);

plot(x, ln_pdf(x), 'k--', 'LineWidth',2, ...
     'DisplayName','Lognormal (MM)');

% % ============================================================
% % Generalized Gamma (moment surrogate)
% % ============================================================
% c_gg = 0.7;  % tail control parameter
% 
% a_gg = (mu^2 / v) / c_gg;
% b_gg = mu * gamma(a_gg) / gamma(a_gg + 1/c_gg);
% 
% genGammaPDF = @(x,a,b,c) (abs(c) ./ (b.^(a.*c) .* gamma(a))) .* ...
%                          x.^(a.*c-1) .* exp(-(x./b).^c);
% 
% gg_pdf = @(x) genGammaPDF(x, a_gg, b_gg, c_gg);
% 
% plot(x, gg_pdf(x), 'm:', 'LineWidth',2, ...
%      'DisplayName','Generalized Gamma (MM)');

% ============================================================
% ✅ CORRECT: Moment-based Lognormal–Gamma MIXTURE
% ============================================================

% Mixture weight (closed-form, no data)
alpha = (vars(1) / varN)*(1 - 1/k_g);   % physically reasonable; can be tied to skewness

mix_pdf = @(x) alpha * ln_pdf(x) + (1-alpha) * g_pdf(x);

plot(x, mix_pdf(x), 'g-', 'LineWidth',2.5, ...
     'DisplayName','MM Lognormal–Gamma Mixture');

% ============================================================
% Formatting
% ============================================================
grid on;
xlabel('x');
ylabel('PDF');
title('N Non-Opt');
xlim([0 1000])
legend('show','Location','best');
%%
figure; hold on;




% ============================================================
% Optimized
% ============================================================
histogram(Nc_opt, 'Normalization', 'pdf', ...
          'FaceColor',[0.7 0.7 0.7], ...
          'EdgeColor','none', ...
          'DisplayName','Monte Carlo');



% ============================================================
% Closed-form moments
% ============================================================
mu = meanN_opt;
v  = varN_opt;

% ============================================================
% Gamma (MM)
% ============================================================
k_g     = mu^2 / v;
theta_g = v / mu;

g_pdf = @(x) gampdf(x, k_g, theta_g);

plot(x, g_pdf(x), 'r-', 'LineWidth',2, ...
     'DisplayName','Gamma (MM)');

% ============================================================
% Lognormal (MM)
% ============================================================
sigma2 = log(1 + v/mu^2);
sigma  = sqrt(sigma2);
mu_l   = log(mu) - 0.5*sigma2;

ln_pdf = @(x) lognpdf(x, mu_l, sigma);

plot(x, ln_pdf(x), 'k--', 'LineWidth',2, ...
     'DisplayName','Lognormal (MM)');



% ============================================================
% ✅ CORRECT: Moment-based Lognormal–Gamma MIXTURE
% ============================================================

% Mixture weight (closed-form, no data)
alpha = (vars_N_Opt(1) / varN_opt)*(1 - 1/k_g);   % physically reasonable; can be tied to skewness

mix_pdf = @(x) alpha * ln_pdf(x) + (1-alpha) * g_pdf(x);

plot(x, mix_pdf(x), 'g-', 'LineWidth',2.5, ...
     'DisplayName','MM Lognormal–Gamma Mixture');

% ============================================================
% Formatting
% ============================================================
grid on;
xlabel('x');
ylabel('PDF');
title('N Opt');
xlim([0 1000])
legend('show','Location','best');

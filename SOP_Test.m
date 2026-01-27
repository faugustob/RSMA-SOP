%% ============================================================
%  SOP & ASC VALIDATION: MONTE-CARLO vs GAUSS-LAGUERRE vs INTEGRAL
%% ============================================================
clear; clc; rng(1);

%% ---------------- PARAMETERS ----------------
Rs = 0.5;       
Nmc = 1e6;      
Ngl = 40;       

% Legitimate Link (L)
alpha = 1.5; beta = 0.2; bL = 1;
a = 0.6; muL = 0; sigmaL = 0.6; kL = 2; thetaL = 1;

% Eavesdropper Link (E)
alphaE = 1.5; betaE = 0.4; bE = 1;
aE = 0.5; muE = 0; sigmaE = 0.7; kE = 2.5; thetaE = 1;

%% ============================================================
%  1. MONTE-CARLO (GROUND TRUTH)
%% ============================================================
uL = rand(Nmc,1); XL = (uL<=a).*lognrnd(muL,sigmaL,Nmc,1) + (uL>a).*gamrnd(kL,thetaL,Nmc,1);
uE = rand(Nmc,1); XE = (uE<=aE).*lognrnd(muE,sigmaE,Nmc,1) + (uE>aE).*gamrnd(kE,thetaE,Nmc,1);

gL = (alpha*XL)./(beta*XL + bL);
gE = (alphaE*XE)./(betaE*XE + bE);
Cs_samples = max(log2(1+gL) - log2(1+gE), 0);

SOP_MC = mean(Cs_samples < Rs);
ASC_MC = mean(Cs_samples);

%% ============================================================
%  2. EXPANDED FORM USING 'INTEGRAL' FUNCTION
%% ============================================================
% Component CDFs and PDFs
FXL = @(x) a*normcdf((log(x)-muL)/sigmaL) + (1-a)*gammainc(x/thetaL, kL);
FXE = @(x) aE*normcdf((log(x)-muE)/sigmaE) + (1-aE)*gammainc(x/thetaE, kE);

fXE = @(x) aE*(1./(x*sigmaE*sqrt(2*pi)).*exp(-(log(x)-muE).^2/(2*sigmaE^2))) + ...
           (1-aE)*(x.^(kE-1).*exp(-x/thetaE)./(thetaE^kE*gamma(kE)));

% --- Expanded SOP Integral ---
% SOP = \int F_XL( T(gamma_E(x)) ) * f_XE(x) dx
SOP_kernel = @(x) arrayfun(@(xe) FXL(max(eps, ...
    (bL*(2^Rs*(1+(alphaE*xe/(betaE*xe+bE)))-1)) / ...
    (alpha - beta*(2^Rs*(1+(alphaE*xe/(betaE*xe+bE)))-1)))), x) .* fXE(x);

% We handle the SINR ceiling (\alpha/\beta) by capping the mapping
SOP_INT = integral(@(x) SOP_kernel(x), 0, Inf);

% --- Expanded ASC Integral ---
% ASC = 1/ln2 * \int (FE(t)*(1-FL(t)))/(1+t) dt
t_max = alpha/beta; % Legitimate ceiling
ASC_kernel = @(t) (arrayfun(@(ti) FXE(max(eps, (bE*ti)/(alphaE-betaE*ti))), t) .* ...
                  (1 - arrayfun(@(ti) FXL(max(eps, (bL*ti)/(alpha-beta*ti))), t))) ./ (1+t);

ASC_INT = (1/log(2)) * integral(@(t) ASC_kernel(t), 0, t_max);

%% ============================================================
%  3. GAUSS-LAGUERRE (YOUR CURRENT IMPLEMENTATION)
%% ============================================================
[xi, wi] = GaussLaguerreReal(Ngl);
SOP_GL = 0; ASC_GL = 0;
for i = 1:Ngl
    u = xi(i); w = wi(i);
    
    % SOP Part
    xE_sop = thetaE * u;
    gE_sop = (alphaE * xE_sop) / (betaE * xE_sop + bE);
    Lambda = 2^Rs * (1 + gE_sop) - 1;
    yL_sop = (alpha > beta*Lambda) * (bL*Lambda/(alpha-beta*Lambda)) + (alpha <= beta*Lambda)*1e15;
    valFXL = FXL(max(yL_sop, eps));
    
    term_sop = w * ((1-aE)*(u^(kE-1)/gamma(kE)) + aE*exp(u)*(1/(xE_sop*sigmaE*sqrt(2*pi))*exp(-(log(xE_sop)-muE)^2/(2*sigmaE^2)))) * valFXL;
    SOP_GL = SOP_GL + term_sop;
    
    % ASC Part
    t = exp(u) - 1;
    FE = (t < alphaE/betaE) * FXE(max(eps, bE*t/(alphaE-betaE*t))) + (t >= alphaE/betaE);
    FL = (t < alpha/beta) * FXL(max(eps, bL*t/(alpha-beta*t))) + (t >= alpha/beta);
    ASC_GL = ASC_GL + w * exp(u) * FE * (1 - FL);
end
ASC_GL = ASC_GL / log(2);

%% ============================================================
%  RESULTS
%% ============================================================
fprintf('METRIC |    MC      |  INTEGRAL  |     GL     \n');
fprintf('-------|------------|------------|------------\n');
fprintf('SOP    |  %.6f  |  %.6f  |  %.6f  \n', SOP_MC, SOP_INT, SOP_GL);
fprintf('ASC    |  %.6f  |  %.6f  |  %.6f  \n', ASC_MC, ASC_INT, ASC_GL);

function [x, w] = GaussLaguerreReal(n)
    syms z; L = laguerreL(n, z); r = vpasolve(L == 0, z, [0, inf]);
    x = double(sort(r)); L_next = laguerreL(n + 1, z);
    w = zeros(size(x));
    for i = 1:length(x)
        w(i) = x(i) / ((n + 1)^2 * double(subs(L_next, z, x(i)))^2);
    end
end
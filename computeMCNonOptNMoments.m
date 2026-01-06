function [Nc,meanNMC,varNMC,meansMC,varsMC] = computeMCNonOptNMoments(I,Pe,Pr,Pq,m_e,m_q,m_p,Omega_e,Omega_p,Omega_q,Plos,PLj,Nr,HB,HA,g_pq,beta_r,Ns,Nsymb)

%% ===================== Parameter Size Enforcement =====================
% Force m_p and Omega_p to be [Pr x 1]
if isscalar(m_p),     m_p = repmat(m_p, Pr, 1);         end
if isscalar(Omega_p), Omega_p = repmat(Omega_p, Pr, 1); end

% Force m_q and Omega_q to be [Pq x 1]
if isscalar(m_q),     m_q = repmat(m_q, Pq, 1);         end
if isscalar(Omega_q), Omega_q = repmat(Omega_q, Pq, 1); end

% Force m_e and Omega_e to be [Pe x 1]
if isscalar(m_e),     m_e = repmat(m_e, Pe, 1);         end
if isscalar(Omega_e), Omega_e = repmat(Omega_e, Pe, 1); end

rng(2);

%% ===================== Generate channel matrices =====================
h_rp = zeros(Nr, Pr, Ns);
for p = 1:Pr
    % Scale parameter theta = Omega/m
    h_rp(:, p, :) = sqrt(gamrnd(m_p(p), Omega_p(p)/m_p(p), [Nr, 1, Ns])) .* exp(1i*2*pi*rand(Nr, 1, Ns));
end

h_jq = zeros(Nr, Pq, Ns);
for q = 1:Pq
    h_jq(:, q, :) = sqrt(gamrnd(m_q(q), Omega_q(q)/m_q(q), [Nr, 1, Ns])) .* exp(1i*2*pi*rand(Nr, 1, Ns));
end

h_e = zeros(Pe, Ns);
for u = 1:Pe
    h_e(u, :) = sqrt(gamrnd(m_e(u), Omega_e(u)/m_e(u), [1, Ns])) .* exp(1i*2*pi*rand(1, Ns));
end

%% ===================== Allocate arrays =====================
Nc_term1 = zeros(Ns,1);
Nc_term2 = zeros(Ns,1);
Nc_term3 = zeros(Ns,1);
Nc = zeros(Ns,1);

%% ===================== Main Computation Loop =====================
for n = 1:Ns
    % ------------------ Construct Term1 ------------------
    term1 = zeros(Nsymb,Nsymb);
    for r = 1:Nr
        beta = beta_r(r);
        for p = 1:Pr
            for q = 1:Pq
                term1 = term1 + beta * h_rp(r,p,n) * h_jq(r,q,n) * g_pq(p,q) * HA(:,:,p,q);
            end
        end
    end
    term1 = sqrt(PLj)*term1;
    
    % ------------------ Construct Term2 ------------------
    B1 = zeros(Nsymb,Nsymb);
    for u = 1:Pe
        B1 = B1 + h_e(u,n) * HB(:,:,u);
    end
    term2 = I * sqrt(Plos) * B1;
    
    % ------------------ Power Calculations ------------------
    Nc_term1(n) = norm(term1,'fro').^2;
    Nc_term2(n) = norm(term2,'fro').^2;
    Nc_term3(n) = 2*real(trace(term2'*term1));
    Nc(n) = norm(term1 + term2,'fro').^2;
end

meansMC = [mean(Nc_term1); mean(Nc_term2); mean(Nc_term3)];
varsMC = [var(Nc_term1); var(Nc_term2); var(Nc_term3)];
meanNMC = mean(Nc);
varNMC = var(Nc);

end
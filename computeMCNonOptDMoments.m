function [Dc,meanDMC,varDMC,meansD_MC,varsD_MC] = computeMCNonOptDMoments(I,Pe,Pr,Pq,m_e,m_q,m_p,Omega_e,Omega_p,Omega_q,Plos,PLj,Nr,HB,HA,g_pq,beta_r,alpha_array,Ns,Nsymb)
rng(2);
%% ===================== Allocate arrays =====================
Nc_term1 = zeros(Ns,1);
Nc_term2 = zeros(Ns,1);
Nc_term3 = zeros(Ns,1);
Dc= zeros(Ns,1);
K = length(alpha_array);

%% ===================== Generate channel matrices =====================


h_rp = sqrt(gamrnd(m_p, Omega_p/m_p, [Nr,Pr,Ns])) .* exp(1i*2*pi*rand(Nr,Pr,Ns));
h_jq = sqrt(gamrnd(m_q, Omega_q/m_q, [Nr,Pq,Ns])) .* exp(1i*2*pi*rand(Nr,Pq,Ns));

h_e = sqrt(gamrnd(m_e, Omega_e/m_e, [Pe,Ns])) .* exp(1i*2*pi*rand(Pe,Ns));

 for n = 1:Ns

    % ------------------ Construct Term1 ------------------

     term1 = zeros(Nsymb,Nsymb);

    for k = 1:K
        for r = 1:Nr
           beta = beta_r(r);
            for p = 1:Pr
                for q = 1:Pq
       
                    term1 = term1 + sqrt(alpha_array(k))*beta * h_rp(r,p,n) * h_jq(r,q,n) * g_pq(p,q) * HA(:,:,p,q);
                end
            end
        end
    end
    
    term1 = sqrt(PLj)*term1;

    norm_term1 = norm(term1,'fro').^2;

   

    % ------------------ Construct Term2 ------------------
    B1 = zeros(Nsymb,Nsymb);

    for k=1:K
        for i = 1:Pe
          B1 = B1 + sqrt(alpha_array(k))*h_e(i,n) * HB(:,:,i);
        end
    end
    
    term2 = I * sqrt(Plos) * B1;

    norm_term2 = norm(term2,'fro').^2;


  
    norm_Term3 = 2*real(trace(term2'*term1));    
    

    Nc_term1(n) = norm_term1;

    Nc_term2(n) = norm_term2;

    Nc_term3(n) =  norm_Term3;

    Dc(n) = norm(term1 + term2,'fro').^2;

 end

 meansD_MC=[mean(Nc_term1);mean(Nc_term2);mean(Nc_term3)];
 varsD_MC = [var(Nc_term1);var(Nc_term2);var(Nc_term3)];

 meanDMC = mean(Dc);
 varDMC = var(Dc);
 


end

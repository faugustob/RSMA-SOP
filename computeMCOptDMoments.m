function [Dc_opt,meanDMC_opt,varDMC_opt,meansD_MC_opt,varsD_MC_opt] = computeMCOptDMoments(I,Pe,Pr,Pq,m_e,m_q,m_p,Omega_e,Omega_p,Omega_q,Plos,PLj,Nr,HB,HA,g_pq,beta_opt,alpha_array,Ns,Nsymb)
rng(2);
%% ===================== Allocate arrays =====================

Dc_opt= zeros(Ns,1);
meansD_MC_opt = zeros(3,1);
varsD_MC_opt = zeros(3,1);
K = length(alpha_array);

Nc_opt_term1 = zeros(Ns,1);
Nc_opt_term2 = zeros(Ns,1);
Nc_opt_term3 = zeros(Ns,1);




h_rp = sqrt(gamrnd(m_p, Omega_p/m_p, [Nr,Pr,Ns]));
h_jq = sqrt(gamrnd(m_q, Omega_q/m_q, [Nr,Pq,Ns]));

h_e = sqrt(gamrnd(m_e, Omega_e/m_e, [Pe,Ns]));

 for n = 1:Ns

    % ------------------ Construct Term1 ------------------

     term1 = zeros(Nsymb,Nsymb);
    for r = 1:Nr
       beta = beta_opt(r);
        for p = 1:Pr
            for q = 1:Pq
   
                term1 = term1 + beta * h_rp(r,p,n) * h_jq(r,q,n) * g_pq(p,q) * HA(:,:,p,q);
            end
        end
    end
    term1 = sqrt(PLj)*term1;


    % ------------------ Construct Term2 ------------------
    B1 = zeros(Nsymb,Nsymb);
    for i = 1:Pe
        B1 = B1 + h_e(i,n) * HB(:,:,i);
    end
    term2 = I * sqrt(Plos) * B1;

    norm_term2 = norm(term2,'fro').^2;


    % norm_term2_expanded = 0;
    % for u=1:Pe
    %     for up=1:Pe
    %         norm_term2_expanded = norm_term2_expanded +  conj(h_e(u,n))*h_e(up,n)*trace(HB(:,:,u)'*HB(:,:,up));
    %     end
    % end
    % 
    % norm_term2_expanded = I * Plos * norm_term2_expanded;


    % norm_Term3_exp = 0;
    % 
    % for r = 1: Nr
    %     beta = beta_r(r);
    %     for u=1:Pe
    %         for p=1:Pr
    %             for q=1:Pq
    %                 norm_Term3_exp = norm_Term3_exp  + beta*h_e(u,n)* h_rp(r,p,n) * h_jq(r,q,n) * g_pq(p,q) * trace(HB(:,:,u)'*HA(:,:,p,q));
    %             end
    %         end
    %     end
    % end
    % 
    % norm_Term3_exp = 2*I*sqrt(Plos)*sqrt(PLj)*norm_Term3_exp;
    norm_Term3 = 2*real(trace(term2'*term1));
    % norm_Term3_diff = norm(term1+term2,'fro').^2-norm(term1,'fro').^2-norm_term2

    
    

    Nc_opt_term1(n) = norm(term1,'fro').^2;

    Nc_opt_term2(n) = norm_term2;

    Nc_opt_term3(n) =  norm_Term3;

    Dc_opt(n) = norm(term1 + term2,'fro').^2;


    % ------------------ Alternative normA_clean ------------------
    % normA_clean = 0;
    % for r=1:Nr
    %     beta_s = conj(beta_r(r));
    %     for rp=1:Nr
    %          beta = beta_r(rp);
    %         for p = 1:Pr
    %             for q = 1:Pq
    %                 for pp=1:Pr
    %                     for qp =1:Pq
    % 
    %                         normA_clean = normA_clean+beta_s*beta*conj(h_rp(r,p,n))*h_rp(rp,pp,n)*conj(h_jq(r,q,n))*h_jq(rp,qp,n)...
    %                         *conj(g_pq(p,q))*g_pq(pp,qp)*trace(HA(:,:,p,q)'*HA(:,:,pp,qp));
    % 
    %                     end
    %                 end
    %             end
    %         end
    %     end
    % end
    % normA_clean = real(normA_clean);

    %Nc_exp(n) = normA_clean;
 end



meansD_MC_opt(1) = mean(Nc_opt_term1);
meansD_MC_opt(2) = mean(Nc_opt_term2);
meansD_MC_opt(3) = mean(Nc_opt_term3);


varsD_MC_opt(1) = var(Nc_opt_term1);
varsD_MC_opt(2) = var(Nc_opt_term2);
varsD_MC_opt(3) = var(Nc_opt_term3);


meanDMC_opt = mean(Dc_opt);
varDMC_opt = var(Dc_opt);


end

function [meanN_opt,varN_opt,mean_N_opt,vars_N_Opt,Covs] = computeOptNMoments(I,Pe,Pr,Pq,m_e,m_q,m_p,Omega_e,Omega_p,Omega_q,Plos,PLj,Nr,HB,HA,g_pq,beta_opt)


        vars_N_Opt = zeros(3,1);
        mean_N_opt = zeros(3,1);
        Covs = zeros(3,1);

        Covs(1) = 0;
    
        Theo_term1_mean = 0;
        for r=1:Nr
            beta_s = conj(beta_opt(r));
            for rp=1:Nr
                 beta = beta_opt(rp);
                for p = 1:Pr
                    for q = 1:Pq
                        for pp=1:Pr
                            for qp =1:Pq                          
    
                                if r==rp && p == pp
                                    mean_p = nakagami_moment(2,m_p,Omega_p);
                                else
                                    mean_p = nakagami_moment(1,m_p,Omega_p)^2;
                                end
    
                                if r==rp && q == qp
                                    mean_q = nakagami_moment(2,m_q,Omega_q);
                                else
                                    mean_q = nakagami_moment(1,m_q,Omega_q)^2;
                                end
                                                         
    
                                Theo_term1_mean = Theo_term1_mean+beta_s*beta*mean_p*mean_q...
                                *conj(g_pq(p,q))*g_pq(pp,qp)*trace(HA(:,:,p,q)'*HA(:,:,pp,qp));
                                
                            end
                        end
                    end
                end
            end
        end
        mean_N_opt(1) = PLj*real(Theo_term1_mean);
    
    
        Theo_term2_mean = 0;
        for u=1:Pe
            for up=1:Pe
                 % mean_e = mean(conj(h_e(u,:)).*h_e(up,:));
                if u==up 
                    mean_e = nakagami_moment(2,m_e,Omega_e);
                else
                    mean_e = nakagami_moment(1,m_e,Omega_e)^2;
                end
                Theo_term2_mean = Theo_term2_mean +  mean_e*trace(HB(:,:,u)'*HB(:,:,up));
            end
        end
        mean_N_opt(2) = I*Plos*real(Theo_term2_mean);
    
        
    
    
        Theo_term3_mean = 0;
    
      
        for r = 1: Nr
            beta = beta_opt(r);
            for u=1:Pe
                for p=1:Pr
                    for q=1:Pq
               
    
                        mean_e = nakagami_moment(1,m_e,Omega_e);
                        mean_p =  nakagami_moment(1,m_p,Omega_p);
                        mean_q =  nakagami_moment(1,m_q,Omega_q);
    
    
    
                        Theo_term3_mean = Theo_term3_mean  + beta* mean_e.* mean_p .* mean_q * g_pq(p,q) * trace(HB(:,:,u)'*HA(:,:,p,q));
                    end
                end
            end
        end
    
        mean_N_opt(3) = 2*I*sqrt(Plos)*sqrt(PLj)*real(Theo_term3_mean);
    
      
    
    
    
    Term1_Theo_2nd = 0;
    
    for r = 1:Nr
        beta = conj(beta_opt(r));
        for rp = 1:Nr
            beta_rp = beta_opt(rp);
            for rpp = 1:Nr
                beta_rpp = conj(beta_opt(rpp));
                for rppp = 1:Nr
                    beta_rppp = beta_opt(rppp);
    
                    for p = 1:Pr
                        for q = 1:Pq
                            for pp = 1:Pr
                                for qp = 1:Pq
                                    for ppp = 1:Pr
                                        for qpp = 1:Pq
                                            for pppp = 1:Pr
                                                for qppp = 1:Pq
    
                                                 mean_p = 0;
                                                 mean_q = 0;
    
                                                  crp = classify4( ...
                                                    sub2ind([Nr,Pr],r,p), ...
                                                    sub2ind([Nr,Pr],rp,pp), ...
                                                    sub2ind([Nr,Pr],rpp,ppp), ...
                                                    sub2ind([Nr,Pr],rppp,pppp));
                                                
                                                switch crp
                                                    case 4      % h^4
                                                        mean_p = nakagami_moment(4,m_p,Omega_p);
                                                
                                                    case 31     % h^3 h
                                                        mean_p = nakagami_moment(3,m_p,Omega_p) * ...
                                                                 nakagami_moment(1,m_p,Omega_p);
                                                
                                                    case 22     % h^2 h'^2
                                                        mean_p = nakagami_moment(2,m_p,Omega_p)^2;
                                                
                                                    case 211    % h^2 h h
                                                        mean_p = nakagami_moment(2,m_p,Omega_p) * ...
                                                                 nakagami_moment(1,m_p,Omega_p)^2;
                                                
                                                    otherwise   % 1111
                                                        mean_p = nakagami_moment(1,m_p,Omega_p)^4;
                                                end
    
                                                cjq = classify4( ...
                                                    sub2ind([Nr,Pq],r,q), ...
                                                    sub2ind([Nr,Pq],rp,qp), ...
                                                    sub2ind([Nr,Pq],rpp,qpp), ...
                                                    sub2ind([Nr,Pq],rppp,qppp));
                                                
                                                switch cjq
                                                    case 4
                                                        mean_q = nakagami_moment(4,m_q,Omega_q);
                                                    case 31
                                                        mean_q = nakagami_moment(3,m_q,Omega_q) * ...
                                                                 nakagami_moment(1,m_q,Omega_q);
                                                    case 22
                                                        mean_q = nakagami_moment(2,m_q,Omega_q)^2;
                                                    case 211
                                                        mean_q = nakagami_moment(2,m_q,Omega_q) * ...
                                                                 nakagami_moment(1,m_q,Omega_q)^2;
                                                    otherwise
                                                        mean_q = nakagami_moment(1,m_q,Omega_q)^4;
                                                end
    
    
                                                    Term1_Theo_2nd = Term1_Theo_2nd + ...
                                                        beta * beta_rp * beta_rpp * beta_rppp *mean_p*mean_q*...
                                                        conj(g_pq(p,q))*g_pq(pp,qp) *conj(g_pq(ppp,qpp))*g_pq(pppp,qppp)*...
                                                        trace(HA(:,:,p,q)'*HA(:,:,pp,qp))*trace(HA(:,:,ppp,qpp)'*HA(:,:,pppp,qppp));
    
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
    
                end
            end
        end
    end
    
    Term1_Theo_2nd = PLj^2*real(Term1_Theo_2nd);
    
    
    
    
    
    
    
    
    Term2_Theo_2nd = 0;
    
    for u = 1:Pe
        for up=1:Pe
            for upp=1:Pe
                for uppp=1:Pe
    
                    % mean_e = mean(conj(h_e(u,:)).*h_e(up,:).*conj(h_e(upp,:)).*h_e(uppp,:));
    
    
                    ce = classify4(u, up, upp, uppp);
    
                    switch ce
                        case 4      % h^4
                            mean_e = nakagami_moment(4,m_e,Omega_e);
                    
                        case 31     % h^3 h
                            mean_e = nakagami_moment(3,m_e,Omega_e) * ...
                                     nakagami_moment(1,m_e,Omega_e);
                    
                        case 22     % h^2 h'^2
                            mean_e = nakagami_moment(2,m_e,Omega_e)^2;
                    
                        case 211    % h^2 h h
                            mean_e = nakagami_moment(2,m_e,Omega_e) * ...
                                     nakagami_moment(1,m_e,Omega_e)^2;
                    
                        otherwise   % 1111
                            mean_e = nakagami_moment(1,m_e,Omega_e)^4;
                    end
                    Term2_Theo_2nd = Term2_Theo_2nd + mean_e*...
                        trace(HB(:,:,u)'*HB(:,:,up))*trace(HB(:,:,upp)'*HB(:,:,uppp));
                end
            end
        end
    end
    Term2_Theo_2nd = I*Plos^2*real(Term2_Theo_2nd);
    
    
    
    
            Term3_Theo_2nd = 0;
    
            for r = 1:Nr
                for rp = 1:Nr
                    for u = 1:Pe
                        for up = 1:Pe
                            for p = 1:Pr
                                for pp = 1:Pr
                                    for q = 1:Pq
                                        for qp = 1:Pq
                                            
                                            % 1. Expectations of Real Nakagami Variables
                                            % mean_e = E[h_eu * h_eup]
                                            if u == up
                                                mean_e = nakagami_moment(2,m_e,Omega_e);
                                            else
                                                mean_e = nakagami_moment(1,m_e,Omega_e)^2;
                                            end
                                            
                                            % mean_p = E[h_rp * h_rpp]
                                            if r == rp && p == pp
                                                mean_p = nakagami_moment(2,m_p,Omega_p);
                                            else
                                                mean_p = nakagami_moment(1,m_p,Omega_p)^2;
                                            end
                                            
                                            % mean_q = E[h_jq * h_rqp] (Note: index is j_q and j_qp)
                                            if r == rp && q == qp
                                                mean_q = nakagami_moment(2,m_q,Omega_q);
                                            else
                                                mean_q = nakagami_moment(1,m_q,Omega_q)^2;
                                            end
                                            
                                            % 2. Complex Constants
                                            T1 = trace(HB(:,:,u)' * HA(:,:,p,q));
                                            T2 = trace(HB(:,:,up)' * HA(:,:,pp,qp));
                                            
                                            C1 = beta_opt(r) * g_pq(p,q) * T1;
                                            C2 = beta_opt(rp) * g_pq(pp,qp) * T2;
                                            
                                            % 3. The Core Expansion: E[ (X+X*)^2 ] = 2*|X|^2 + 2*Re(X^2)
                                            % We sum (conj(C1)*C2 + C1*C2)
                                            complex_sum = (conj(C1) * C2) + (C1 * C2);
                                            
                                            Term3_Theo_2nd = Term3_Theo_2nd + ...
                                                mean_e * mean_p * mean_q * complex_sum;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            % Final scaling: 2 * (I * sqrt(Plos) * sqrt(PLj))^2 
            % The '2' comes from the identity (X+X*)^2 = 2*Re(...)
            Term3_Theo_2nd = 2 * I^2 * Plos * PLj * real(Term3_Theo_2nd);
    
    
            
    
    % Cov_T1_T3 = 0;
    % 
    % for r = 1:Nr
    %     beta = beta_opt(r);
    %     for rp = 1:Nr
    %         beta_rp = beta_opt(rp);
    %         for rpp = 1:Nr
    %             beta_rpp = beta_opt(rpp);
    %             for u = 1:Pe
    %                 for p = 1:Pr
    %                     for pp = 1:Pr
    %                         for ppp = 1:Pr
    %                             for q = 1:Pq
    %                                 for qp = 1:Pq
    %                                     for qpp = 1:Pq
    % 
    % 
    %                                         E = mean(reshape(h_e(u,:),1,1,[]).*h_rp(r,p,:).*h_rp(rp,pp,:).*h_rp(rpp,ppp,:)...
    %                                             .*h_jq(r,q,:).*h_jq(rp,qp,:).*h_jq(rpp,qpp,:));
    % 
    %                                         Y = conj(beta)*beta_rp*conj(g_pq(p,q))*g_pq(pp,qp)*trace(HA(:,:,p,q)'*HA(:,:,pp,qp));
    % 
    %                                         X=beta_rpp*g_pq(ppp,qpp)*trace(HB(:,:,u)' * HA(:,:,ppp,qpp));
    % 
    %                                         X = X+conj(X);
    % 
    %                                         C = E*Y*X;
    % 
    %                                         % ===== BOTH contractions =====
    %                                         Cov_T1_T3 = Cov_T1_T3 + C;
    %                                     end
    %                                 end
    %                             end
    %                         end
    %                     end
    %                 end
    %             end
    %         end
    %     end
    % end
    % 
    % Cov_T1_T3 = ...
    %     I * sqrt(Plos) * PLj^(3/2) * real(Cov_T1_T3) ...
    %     - mean_N_opt(1) * mean_N_opt(3)
    
    
    Cov_T1_T3 = 0;
    
    for r = 1:Nr
        beta = beta_opt(r);
        for rp = 1:Nr
            beta_rp = beta_opt(rp);
            for rpp = 1:Nr
                beta_rpp = beta_opt(rpp);
                for u = 1:Pe
                    for p = 1:Pr
                        for pp = 1:Pr
                            for ppp = 1:Pr
                                for q = 1:Pq
                                    for qp = 1:Pq
                                        for qpp = 1:Pq
    
                                            cond_p1 = (r == rp  && p == pp);
                                            cond_p2 = (rp == rpp && pp == ppp);
                                            cond_p3 = (r == rpp && p == ppp);
    
                                            mean_e = nakagami_moment(1,m_e,Omega_e);
    
                                             
                                            % mean_p = mean(h_rp(r,p,:).*h_rp(rp,pp,:).*h_rp(rpp,ppp,:));
                                            if r==rp && rp==rpp && p==pp && pp==ppp
                                              mean_p = nakagami_moment(3,m_p,Omega_p);
                                            elseif (cond_p1 + cond_p2 + cond_p3) == 1
                                              mean_p = nakagami_moment(2,m_p,Omega_p)*nakagami_moment(1,m_p,Omega_p);
                                            else
                                              mean_p = nakagami_moment(1,m_p,Omega_p)^3;
                                            end
                                           
                                            % mean_q = mean(h_jq(r,q,:).*h_jq(rp,qp,:).*h_jq(rpp,qpp,:));
    
                                             if r==rp && rp==rpp && q==qp && qp==qpp
                                                % Triple match
                                                mean_q = nakagami_moment(3,m_q,Omega_q);
                                            elseif (r==rp && q==qp) || (rp==rpp && qp==qpp) || (r==rpp && q==qpp)
                                                % Pair match
                                                mean_q = nakagami_moment(2,m_q,Omega_q) * nakagami_moment(1,m_q,Omega_q);
                                            else
                                                % All independent
                                                mean_q = nakagami_moment(1,m_q,Omega_q)^3;
                                            end
                                          
                                            E = mean_e*mean_p*mean_q;
    
                                            Y = conj(beta)*beta_rp*conj(g_pq(p,q))*g_pq(pp,qp)*trace(HA(:,:,p,q)'*HA(:,:,pp,qp));
    
                                            X=beta_rpp*g_pq(ppp,qpp)*trace(HB(:,:,u)' * HA(:,:,ppp,qpp));
    
                                            X = X+conj(X);
    
                                            C = E*Y*X;
    
                                            % ===== BOTH contractions =====
                                            Cov_T1_T3 = Cov_T1_T3 + C;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    Cov_T1_T3 = ...
        I * sqrt(Plos) * PLj^(3/2) * real(Cov_T1_T3) ...
        - mean_N_opt(1) * mean_N_opt(3);
    
 Covs(2) = Cov_T1_T3;
    
    
    Cov_T2_T3 = 0;
    
    for u = 1:Pe
        for up = 1:Pe
            for upp = 1:Pe
                for r = 1:Nr
                    beta=beta_opt(r);
                    for p = 1:Pr
                        for q = 1:Pq
                         
                            % mean_e = mean(h_e(u,:).*h_e(up,:).*h_e(upp,:));
    
                             if u==up && up==upp
                                % Triple match
                                mean_e = nakagami_moment(3,m_e,Omega_e);
                            elseif ( u==up) || ( up==upp) || (u==upp)
                                % Pair match
                                mean_e = nakagami_moment(2,m_e,Omega_e) * nakagami_moment(1,m_e,Omega_e);
                            else
                                % All independent
                                mean_e = nakagami_moment(1,m_e,Omega_e)^3;
                            end
    
                            mean_p =  nakagami_moment(1,m_p,Omega_p);
                            mean_q = nakagami_moment(1,m_q,Omega_q);
    
                            E = mean_e*mean_p*mean_q;
    
                            Y = trace(HB(:,:,up)' * HB(:,:,upp));
    
                           X=beta*g_pq(p,q)*trace(HB(:,:,u)' * HA(:,:,p,q));
                           X = X+conj(X);
    
                            
                            % ---------- Accumulate ----------
                            Cov_T2_T3 = Cov_T2_T3 + ...
                                E*Y*X;
                        end
                    end
                end
            end
        end
    end
    
    Cov_T2_T3 = ...
        I * Plos^(3/2) * sqrt(PLj) * real(Cov_T2_T3) ...
        - mean_N_opt(2) * mean_N_opt(3)
     Covs(3) = Cov_T2_T3;
    
       
   
    vars_N_Opt(1) = (Term1_Theo_2nd - mean_N_opt(1)^2);    
    vars_N_Opt(2)  = (Term2_Theo_2nd-mean_N_opt(2)^2);
    vars_N_Opt(3) = (Term3_Theo_2nd-mean_N_opt(3)^2);  
    varN_opt = vars_N_Opt(1)+ vars_N_Opt(2) + vars_N_Opt(3) + 2*(Cov_T1_T3+Cov_T2_T3);
    meanN_opt = mean_N_opt(1)+mean_N_opt(2)+mean_N_opt(3);

 end



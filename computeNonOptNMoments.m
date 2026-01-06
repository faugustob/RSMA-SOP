function [meanN,varN,means,vars] = computeNonOptNMoments(I,Pe,Pr,Pq,m_e,m_q,m_p,Omega_e,Omega_p,Omega_q,Plos,PLj,Nr,HB,HA,g_pq,beta_r)


    vars = zeros(3,1);
    means = zeros(3,1);

    %% ===================== Parameter Size Enforcement =====================
    if isscalar(m_p),     m_p = repmat(m_p, Pr, 1);         end
    if isscalar(Omega_p), Omega_p = repmat(Omega_p, Pr, 1); end
    if isscalar(m_q),     m_q = repmat(m_q, Pq, 1);         end
    if isscalar(Omega_q), Omega_q = repmat(Omega_q, Pq, 1); end
    if isscalar(m_e),     m_e = repmat(m_e, Pe, 1);         end
    if isscalar(Omega_e), Omega_e = repmat(Omega_e, Pe, 1); end



  % Theo_term1_mean = 0;
  %   for r=1:Nr
  %       beta_s = conj(beta_r(r));
  %       for rp=1:Nr
  %            beta = beta_r(rp);
  %           for p = 1:Pr
  %               for q = 1:Pq
  %                   for pp=1:Pr
  %                       for qp =1:Pq
  % 
  % 
  % 
  %                           Theo_term1_mean = Theo_term1_mean+beta_s*beta*mean(conj(h_rp(r,p,:)).*h_rp(rp,pp,:).*conj(h_jq(r,q,:)).*h_jq(rp,qp,:))...
  %                           *conj(g_pq(p,q))*g_pq(pp,qp)*trace(HA(:,:,p,q)'*HA(:,:,pp,qp));
  % 
  %                       end
  %                   end
  %               end
  %           end
  %       end
  %   end
  %   mean_N_opt(1) = PLj*real(Theo_term1_mean);

   

     Theo_term1_mean = 0;
    for r=1:Nr
        beta_s = conj(beta_r(r));
        for rp=1:Nr
             beta = beta_r(rp);
            for p = 1:Pr
                for q = 1:Pq
                    for pp=1:Pr
                        for qp =1:Pq

                            mean_p = 0;
                            mean_q = 0;

                            if(r==rp && p==pp)
                                mean_p = Omega_p(p);
                            end

                            if(r==rp && q==qp)
                                mean_q = Omega_q(q);
                            end

                            

                            Theo_term1_mean = Theo_term1_mean+beta_s*beta*mean_p*mean_q...
                            *conj(g_pq(p,q))*g_pq(pp,qp)*trace(HA(:,:,p,q)'*HA(:,:,pp,qp));
                            
                        end
                    end
                end
            end
        end
    end
    means(1) = PLj*real(Theo_term1_mean);


    mean_NC_term2_theo = 0;
    % for u=1:Pe
    %     for up=1:Pe
    %         mean_NC_term2_theo = mean_NC_term2_theo +  mean(conj(h_e(u,:)).*h_e(up,:))*trace(HB(:,:,u)'*HB(:,:,up));
    %     end
    % end

    for u=1:Pe
        for up=1:Pe
            mean_e=0;
            if u==up
                mean_e = Omega_e(u);
            end
            mean_NC_term2_theo = mean_NC_term2_theo +  mean_e*trace(HB(:,:,u)'*HB(:,:,up));
        end
    end

    means(2) = I * Plos * real(mean_NC_term2_theo);


    mean_NC_term3_theo = 0;

    % for r = 1: Nr
    %     beta = beta_r(r);
    %     for u=1:Pe
    %         for p=1:Pr
    %             for q=1:Pq
    %                 mean_NC_term3_theo = mean_NC_term3_theo  + beta*mean(reshape(conj(h_e(u,:)),1,1,[]).* h_rp(r,p,:) .* h_jq(r,q,:)) * g_pq(p,q) * trace(HB(:,:,u)'*HA(:,:,p,q));
    %             end
    %         end
    %     end
    % end

    % for r = 1: Nr
    %     beta = beta_r(r);
    %     for u=1:Pe
    %         for p=1:Pr
    %             for q=1:Pq
    %                 mean_NC_term3_theo = mean_NC_term3_theo  + beta*mean(reshape(conj(h_e(u,:)),1,1,[]).* h_rp(r,p,:) .* h_jq(r,q,:)) * g_pq(p,q) * trace(HB(:,:,u)'*HA(:,:,p,q));
    %             end
    %         end
    %     end
    % end
    % 
    % mean_NC_term3_theo = 2*I*sqrt(Plos)*sqrt(PLj)*real(mean_NC_term3_theo);

    means(3) = 0;

% Term1_Theo_2nd = 0;
% 
% for r = 1:Nr
%     beta = conj(beta_r(r));
%     for rp = 1:Nr
%         beta_rp = beta_r(rp);
%         for rpp = 1:Nr
%             beta_rpp = conj(beta_r(rpp));
%             for rppp = 1:Nr
%                 beta_rppp = beta_r(rppp);
% 
%                 for p = 1:Pr
%                     for q = 1:Pq
%                         for pp = 1:Pr
%                             for qp = 1:Pq
%                                 for ppp = 1:Pr
%                                     for qpp = 1:Pq
%                                         for pppp = 1:Pr
%                                             for qppp = 1:Pq
% 
% 
% 
%                                                  % mean_p = 0;
%                                                  %  mean_q = 0;
%                                                  % 
%                                                  %   if(r == rp && rp == rpp && rpp == rppp && p == pp && pp == ppp && ppp == pppp)
%                                                  %       mean_p =  mean(abs(h_rp(r,p,:)).^4);      
%                                                  %   elseif(r==rpp && p == ppp && rp == rppp && pp==pppp)
%                                                  %       mean_p=mean(.^2).*mean(h_rp(rp,pp,:).^2);
%                                                  %   elseif(r==rp&&rpp==rppp&&p==pp&&ppp==pppp)
%                                                  %       mean_p=mean(abs(h_rp(r,p,:).^2))*mean(abs(h_rp(rpp,ppp,:).^2));
%                                                  %   elseif(r==rppp && rp==rpp && p==pppp && pp==ppp)
%                                                  %       mean_p=mean(abs(h_rp(r,p,:).^2))*mean(abs(h_rp(rp,pp,:).^2));
%                                                  %   end
%                                                  % 
%                                                  % 
%                                                  % 
%                                                  %   if(r == rp && rp == rpp && rpp == rppp && q == qp && qp == qpp && qpp == qppp)                                                       
%                                                  %       mean_q =  mean(abs(h_jq(r,q,:)).^4);
%                                                  %   elseif(r==rpp && q == qpp && rp == rppp && qp==qppp)
%                                                  %       mean_q= mean(conj(h_jq(r,q,:)).^2).*mean(h_jq(rp,qp,:).^2);
%                                                  %   elseif(r==rp&&rpp==rppp&&q==qp&&qpp==qppp)
%                                                  %       mean_q=mean(abs(h_jq(r,q,:).^2))*mean(abs(h_jq(rpp,qpp,:).^2));
%                                                  %    elseif(r==rppp && rp==rpp && q==qppp && qp==qpp)
%                                                  %       mean_q=mean(abs(h_jq(r,q,:).^2))*mean(abs(h_jq(rp,qp,:).^2));
%                                                  %   end
% 
% 
% 
%                                                 Term1_Theo_2nd = Term1_Theo_2nd + ...
%                                                     beta * beta_rp * beta_rpp * beta_rppp *mean(conj(h_rp(r,p,:)).*h_rp(rp,pp,:).*conj(h_rp(rpp,ppp,:)).*h_rp(rppp,pppp,:))...
%                                                     *mean(conj(h_jq(r,q,:)).*h_jq(rp,qp,:).*conj(h_jq(rpp,qpp,:)).*h_jq(rppp,qppp,:))*...
%                                                     conj(g_pq(p,q))*g_pq(pp,qp) *conj(g_pq(ppp,qpp))*g_pq(pppp,qppp)*...
%                                                     trace(HA(:,:,p,q)'*HA(:,:,pp,qp))*trace(HA(:,:,ppp,qpp)'*HA(:,:,pppp,qppp));
% 
%                                             end
%                                         end
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 end
% 
%             end
%         end
%     end
% end
% 
% Term1_Theo_2nd = PLj^2*real(Term1_Theo_2nd);




Term1_Theo_2nd = 0;

for r = 1:Nr
    beta = conj(beta_r(r));
    for rp = 1:Nr
        beta_rp = beta_r(rp);
        for rpp = 1:Nr
            beta_rpp = conj(beta_r(rpp));
            for rppp = 1:Nr
                beta_rppp = beta_r(rppp);

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

                                                   if(r == rp && rp == rpp && rpp == rppp && p == pp && pp == ppp && ppp == pppp)                                                           
                                                       mean_p = Omega_p(p)^2+Omega_p(p)^2/m_p(p);
                                                   elseif(r==rpp && p == ppp && rp == rppp && pp==pppp)
                                                       mean_p=0;
                                                   elseif(r==rp&&rpp==rppp&&p==pp&&ppp==pppp)
                                                       mean_p=Omega_p(p)*Omega_p(ppp);
                                                   elseif(r==rppp && rp==rpp && p==pppp && pp==ppp)
                                                       mean_p=Omega_p(p)*Omega_p(pp);
                                                   end



                                                   if(r == rp && rp == rpp && rpp == rppp && q == qp && qp == qpp && qpp == qppp)                                                       
                                                       mean_q = Omega_q(q)^2+Omega_q(q)^2/m_q(q);
                                                   elseif(r==rpp && q == qpp && rp == rppp && qp==qppp)
                                                       mean_q= 0;
                                                   elseif(r==rp&&rpp==rppp&&q==qp&&qpp==qppp)
                                                       mean_q=Omega_q(q)*Omega_q(qpp);
                                                    elseif(r==rppp && rp==rpp && q==qppp && qp==qpp)
                                                       mean_q=Omega_q(q)*Omega_q(qp);
                                                   end

                     

                                                Term1_Theo_2nd = Term1_Theo_2nd + ...
                                                    beta * beta_rp * beta_rpp * beta_rppp *mean_p.*mean_q*...
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

% for u = 1:Pe
%     for up=1:Pe
%         for upp=1:Pe
%             for uppp=1:Pe
%                 Term2_Theo_2nd = Term2_Theo_2nd + mean(conj(h_e(u,:)).*h_e(up,:).*conj(h_e(upp,:)).*h_e(uppp,:))*...
%                     trace(HB(:,:,u)'*HB(:,:,up))*trace(HB(:,:,upp)'*HB(:,:,uppp));
%             end
%         end
%     end
% end

for u = 1:Pe
    for up=1:Pe
        for upp=1:Pe
            for uppp=1:Pe

                mean_e = 0;
           
        
                if(u == up && up == upp && upp == uppp)                                                           
                    mean_e = Omega_e(u)^2+Omega_e(u)^2/m_e(u);
                elseif(u == upp && up == uppp && up==uppp)
                   mean_e=0;
                elseif(u==up&&upp==uppp)
                   mean_e=Omega_e(u)*Omega_e(upp);
                elseif(u==uppp && up==upp)
                   mean_e=Omega_e(u)*Omega_e(up);
                end

                Term2_Theo_2nd = Term2_Theo_2nd + mean_e*...
                    trace(HB(:,:,u)'*HB(:,:,up))*trace(HB(:,:,upp)'*HB(:,:,uppp));
            end
        end
    end
end
Term2_Theo_2nd = I^2 * Plos^2 * real(Term2_Theo_2nd);


%  Term3_Theo_2nd = 0;
% 
% for r = 1: Nr
%     beta = beta_r(r);
%     for rp = 1: Nr
%         beta_p = beta_r(rp);
%         for u = 1:Pe
%             for up = 1:Pe
%                 for p = 1:Pr
%                     for pp = 1:Pr
%                         for q = 1:Pq
%                             for qp = 1:Pq
% 
%                                 mean_u=0;
%                                 mean_p=0;
%                                 mean_q = 0;
% 
%                                 if r==rp && u == up
%                                     mean_u = mean(abs(h_e(u,:)).^2);
%                                 else
%                                     mean_u = mean(reshape(h_e(u,:),1,1,[]) .* reshape(conj(h_e(up,:)),1,1,[]) );
% 
%                                 end
% 
%                                 if r==rp && p == pp
%                                     mean_p = mean(abs(h_rp(r,p,:)).^2);
%                                 else
%                                     mean_p = mean(conj(h_rp(r,p,:)) .* h_rp(rp,pp,:) );
%                                 end
% 
%                                  if r==rp && q == qp
%                                     mean_q = mean(abs(h_jq(r,q,:)).^2);
%                                 else
%                                     mean_q = mean(conj(h_jq(r,q,:)) .* h_jq(rp,qp,:) );
%                                 end
% 
%                                 Term3_Theo_2nd = Term3_Theo_2nd + ...
%                                     conj(beta)*beta_p * ...
%                                     mean_u*mean_p*mean_q.* ...
%                                     conj(g_pq(p,q)) .* g_pq(pp,qp) .* ...
%                                     conj(trace(HB(:,:,u)' * HA(:,:,p,q))) .* ...
%                                     trace(HB(:,:,up)' * HA(:,:,pp,qp));
% 
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% % ---- Correct variance scaling (KEEP mean symbolic)
% Term3_Theo_2nd = 2 * I^2 * Plos * PLj * real(Term3_Theo_2nd);


    Term3_Theo_2nd = 0;

    for r = 1:Nr
        for u = 1:Pe
            for p = 1:Pr
                for q = 1:Pq
                    Term3_Theo_2nd = Term3_Theo_2nd + ...
                        abs(beta_r(r))^2 * ...
                        Omega_e(u) * Omega_p(p) * Omega_q(q) * ...
                        abs(g_pq(p,q))^2 * ...
                        abs(trace(HB(:,:,u)'*HA(:,:,p,q)))^2;
                end
            end
        end
    end
    
    Term3_Theo_2nd = 2 * I^2 * Plos * PLj * Term3_Theo_2nd;


    

    vars(1) = (Term1_Theo_2nd - means(1)^2);
    
    
    
    vars(2)  = (Term2_Theo_2nd-means(2)^2);
    
    
    
    
    vars(3) = (Term3_Theo_2nd-means(3)^2);
    
    
     meanN = means(1)+ means(2) + means(3)
    
     varN = vars(1)+ vars(2) + vars(3)
end



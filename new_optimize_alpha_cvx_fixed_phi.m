function [alpha] = new_optimize_alpha_cvx_fixed_phi(Rmin,alpha_prev,L_node,E_node,phi_St, phi_Sr, zeta_k_St, ...
    K, nF, reflect,  delta_f, Active_Gain_dB, max_SCA)
%% ========================= CONSTANTS =========================
zeta_k_Sr = (10^(Active_Gain_dB/10)) - zeta_k_St;
phase_St = exp(1j .* phi_St);
phase_Sr = exp(1j .* phi_Sr);
beta_St = sqrt(zeta_k_St) .* phase_St;
beta_Sr = sqrt(zeta_k_Sr) .* phase_Sr;
BW = delta_f;
N0_dBm = -174;
sigma2 = 10^((N0_dBm + 10*log10(BW) - 30)/10);
Pw_dBm = 46;
Pw = 10^((Pw_dBm - 30)/10);
AN_P_ratio = 1;          % Increase this (e.g. 5-10) if eavesdroppers are too strong
noise = sigma2/Pw;


%% ========================= PRECOMPUTE CHANNELS =========================
Pk = zeros(K,1);
Ak = zeros(K,1);
Pl = zeros(nF,1);
Al = zeros(nF,1);

for k = 1:K
     reflect_coeff = reflect(k);
        beta_r = (reflect_coeff == 1) * beta_Sr + (reflect_coeff == -1) * beta_St;
 
     beta = beta_r.';
     Pk(k) = real(beta' * L_node(k).V1 * beta + 2*real(beta' * L_node(k).V2) + L_node(k).term3);   % you already fixed this earlier
     Ak(k) = real(beta' * L_node(k).V1_AN * beta + 2*real(beta' * L_node(k).V2_AN) + L_node(k).term3_AN);    
end
beta =  beta_St.';
for l = 1:nF       
 
  Pl(l) = real(beta' * E_node(l).V1_l * beta + 2*real(beta' * E_node(l).V2_l) + E_node(l).term3_l);
  Al(l) = real(beta' * E_node(l).V1_AN_l * beta + 2*real(beta' * E_node(l).V2_AN_l) + E_node(l).term3_AN_l);        
end


%% ========================= INITIALIZATION =========================
tol = 1e-6;           % Convergence tolerance
obj_prev = -inf;      % Track previous objective value
alpha_prev = alpha_prev.';
lambda_penalty = 1e3; % Adjust based on how strictly you want to enforce Rmin

A_pos=diag(ones(size(alpha_prev)));
A_neg = 1 - A_pos;



A_neg_pi = A_neg;
A_neg_pi(1,:)=0;

scale = 1e10;

noise = noise*scale;

Pk = Pk * scale;
Pl = Pl * scale;
Ak = Ak * scale;
Al = Al * scale;  

%% ========================= SCA LOOP =========================
for sca_iter = 1:max_SCA

    % Scaling to avoid ILL Posed issue
     

   cvx_clear;
    cvx_begin quiet
        cvx_solver mosek
        
        variable vecAlpha(K+1) nonnegative   
        variable Rc nonnegative
        variable Ck(K) nonnegative
        variable Rk(K) nonnegative       
        %variable s_fake(nF,K)            % can be negative (we take max(0,.) later if needed)
        variable t

        % --- NEW: Define Penalty Expression ---
        % We use 'pos' because max(0, Rmin - s_fake) is convex
        % We square it to create a quadratic penalty (as in your example)
        % penalty_term = sum(sum(square_pos(Rmin - s_fake)));
        penalty_term = 0;

        % --- UPDATED: Modified Objective ---
        % We subtract the penalty because we are maximizing
        maximize( t - lambda_penalty * penalty_term )

        subject to
            sum(vecAlpha) <= 1;
            vecAlpha >= 1e-2;

            
            % ---------- USER-LEVEL CONSTRAINTS ----------
            for k = 1:K             

                               
                I_c_prev = Pk(k)*A_neg(:,1).'*alpha_prev + AN_P_ratio * Ak(k)+noise;

                S_c = Pk(k)*A_pos(:,1).'*vecAlpha;
                I_c = Pk(k)*A_neg(:,1).'*vecAlpha + AN_P_ratio * Ak(k)+noise;

                %Original : Rc <= log(S_c+I_c)/log(2)-log(I_c)/log(2);
 
                Rc <= log(S_c+I_c)/log(2) ...
                      - log(I_c_prev)/log(2) ...
                      - (1/log(2)) * ((Pk(k)*A_neg(:,1).')/(I_c_prev)) ...
                        * (vecAlpha-alpha_prev);  

                

                 S_k = Pk(k)*A_pos(:,k+1).'*vecAlpha;  
                 I_k = Pk(k)*A_neg_pi(:,k+1).'*vecAlpha + AN_P_ratio * Ak(k)+noise;

                 I_k_prev = Pk(k)*A_neg_pi(:,k+1).'*alpha_prev + AN_P_ratio * Ak(k)+noise; 
                 %S_k_prev(k) = Pk(k)*A_pos(:,k+1).'*alpha_prev; 

                 Rk(k) <= log(S_k+I_k)/log(2) ...
                      - log(I_k_prev)/log(2) ...
                      - (1/log(2)) * ((Pk(k)*A_neg_pi(:,k+1).')/(I_k_prev)) ...
                        * (vecAlpha-alpha_prev);

                 

             
            end
            sum(Ck)<=Rc;

            for k=1:K
                Ck(k) + Rk(k)>=Rmin; %QoS.
            end
        
            % ---------- EAVESDROPPER / SECRECY CONSTRAINTS ----------
            for l = 1:nF

                   

                for k = 1:K
                                     

                     S_l_prev = Pl(l)*A_pos(:,k+1).'*alpha_prev;
                    I_l_prev = Pl(l)*A_neg_pi(:,k+1).'*alpha_prev + AN_P_ratio * Al(l)+noise;
                  
                    

                    S_l = Pl(l)*A_pos(:,k+1).'*vecAlpha;
                    I_l = Pl(l)*A_neg_pi(:,k+1).'*vecAlpha + AN_P_ratio * Al(l)+noise;  

                                     
                    % original
                    % s_fake(l,k) <= log(I_k+S_k)/log(2) - log(I_k)/log(2) - log(I_l+S_l)/log(2) + log(I_l)/log(2);
                    % s_fake(l,k) <= Rk(k) - log(I_l+S_l)/log(2) + log(I_l)/log(2);

                   grad_IlSl = Pl(l)*(A_pos(:,k+1).' + A_neg_pi(:,k+1).');

                  % check below
                   %s_fake(l,k) <= Rk(k)+ log(I_l)/log(2)- log(I_l_prev+S_l_prev)/log(2)-(1/log(2))*(grad_IlSl/(I_l_prev+S_l_prev)) * (vecAlpha-alpha_prev);
                   t <= Rk(k)+ log(I_l)/log(2)- log(I_l_prev+S_l_prev)/log(2)-(1/log(2))*(grad_IlSl/(I_l_prev+S_l_prev)) * (vecAlpha-alpha_prev);
                  % t<= s_fake(l,k);
                end
            end

    cvx_end

    if ~strcmp(cvx_status,'Solved') && ~strcmp(cvx_status,'Inaccurate/Solved')
        fprintf('SCA failed at iter %d: %s\n', sca_iter, cvx_status);
        break;
    end

    % ---------- CONVERGENCE CHECK ----------
    current_obj = double(t);
    % Check relative change in objective and variable
    obj_change = abs(current_obj - obj_prev) / (abs(obj_prev) + 1);
    alpha_change = norm(double(vecAlpha) - alpha_prev) / (norm(alpha_prev) + 1);

    %fprintf('SCA Iter %d: Obj = %.6f, Alpha Delta = %.6e\n', sca_iter, current_obj, alpha_change);

    if obj_change < tol && alpha_change < tol
        fprintf('SCA converged at iteration %d.\n', sca_iter);
        alpha_prev = double(vecAlpha);
        break;
    end

    % ---------- UPDATE FOR NEXT ITERATION ----------
    alpha_prev = double(vecAlpha);
    obj_prev = current_obj;
 
end
alpha = double(vecAlpha);
alpha = alpha.';

end
function [obj_prev,alpha,Ck_out,feasible_flag,xi_val] = new_optimize_alpha_cvx_fixed_phi(Rmin,alpha_prev,L_node,E_node,phi_St, phi_Sr, zeta_k_St, ...
    K, nF, reflect, delta_f, Active_Gain_dB,AN_P_ratio, max_SCA)
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
noise = sigma2/Pw;

%% ========================= PRECOMPUTE CHANNELS =========================
Pk = zeros(K,1); Ak = zeros(K,1); Pl = zeros(nF,1); Al = zeros(nF,1);
for k = 1:K
    reflect_coeff = reflect(k);
    beta_r = (reflect_coeff == 1) * beta_Sr + (reflect_coeff == -1) * beta_St;
    beta = beta_r.';
    Pk(k) = real(beta' * L_node(k).V1 * beta + 2*real(beta' * L_node(k).V2) + L_node(k).term3);
    Ak(k) = real(beta' * L_node(k).V1_AN * beta + 2*real(beta' * L_node(k).V2_AN) + L_node(k).term3_AN);
end
beta = beta_St.';
for l = 1:nF
    Pl(l) = real(beta' * E_node(l).V1_l * beta + 2*real(beta' * E_node(l).V2_l) + E_node(l).term3_l);
    Al(l) = real(beta' * E_node(l).V1_AN_l * beta + 2*real(beta' * E_node(l).V2_AN_l) + E_node(l).term3_AN_l);
end

%% ========================= INITIALIZATION =========================
tol = 1e-6;
obj_prev = -inf;
alpha_prev = alpha_prev.';                  % column vector
lambda_penalty = 10;
A_pos = diag(ones(size(alpha_prev)));
A_neg = 1 - A_pos;
A_neg_pi = A_neg;
A_neg_pi(1,:) = 0;

%% ========================= ADAPTIVE RANGE-CENTERING =========================
all_terms = [Pk(:); Ak(:); Pl(:); Al(:); noise];
valid_terms = all_terms(all_terms > 0);
if ~isempty(valid_terms)
    scale = 10^(-(log10(min(valid_terms)) + log10(max(valid_terms)))/2);
else
    scale = 1;
end
noise = noise * scale; Pk = Pk*scale; Pl = Pl*scale; Ak = Ak*scale; Al = Al*scale;


%% ========================= SCA LOOP (safe retry version) =========================
sca_iter = max_SCA;


while sca_iter > 0
    cvx_clear;
    cvx_begin quiet
        cvx_solver mosek
        variable vecAlpha(K+1) nonnegative
        variable Rc nonnegative
        variable Ck(K) nonnegative
        variable Rk(K) nonnegative
        variable xi(K) nonnegative
        variable t

        penalty_term = sum(xi);
        maximize( t - lambda_penalty * penalty_term )

        subject to
            sum(vecAlpha) <= 1;
            %vecAlpha >= 0;

            % USER-LEVEL CONSTRAINTS
            for k = 1:K
                I_c_prev = Pk(k)*A_neg(:,1).'*alpha_prev + AN_P_ratio * Ak(k)+noise;
                S_c = Pk(k)*A_pos(:,1).'*vecAlpha;
                I_c = Pk(k)*A_neg(:,1).'*vecAlpha + AN_P_ratio * Ak(k)+noise;

                Rc <= log(S_c+I_c)/log(2) ...
                      - log(I_c_prev)/log(2) ...
                      - (1/log(2)) * ((Pk(k)*A_neg(:,1).')/(I_c_prev)) ...
                        * (vecAlpha-alpha_prev);

                S_k = Pk(k)*A_pos(:,k+1).'*vecAlpha;
                I_k = Pk(k)*A_neg_pi(:,k+1).'*vecAlpha + AN_P_ratio * Ak(k)+noise;
                I_k_prev = Pk(k)*A_neg_pi(:,k+1).'*alpha_prev + AN_P_ratio * Ak(k)+noise;

                Rk(k) <= log(S_k+I_k)/log(2) ...
                      - log(I_k_prev)/log(2) ...
                      - (1/log(2)) * ((Pk(k)*A_neg_pi(:,k+1).')/(I_k_prev)) ...
                        * (vecAlpha-alpha_prev);
            end

            sum(Ck) <= Rc;
            for k=1:K
                Ck(k) + Rk(k) + xi(k) >= Rmin;
            end

            % SECRECY CONSTRAINTS
            for l = 1:nF
                for k = 1:K
                    S_l_prev = Pl(l)*A_pos(:,k+1).'*alpha_prev;
                    I_l_prev = Pl(l)*A_neg_pi(:,k+1).'*alpha_prev + AN_P_ratio * Al(l)+noise;
                    S_l = Pl(l)*A_pos(:,k+1).'*vecAlpha;
                    I_l = Pl(l)*A_neg_pi(:,k+1).'*vecAlpha + AN_P_ratio * Al(l)+noise;

                    grad_IlSl = Pl(l)*(A_pos(:,k+1).' + A_neg_pi(:,k+1).');

                    t <= Rk(k) + log(I_l)/log(2) ...
                         - log(I_l_prev + S_l_prev)/log(2) ...
                         - (1/log(2))*(grad_IlSl/(I_l_prev + S_l_prev)) * (vecAlpha-alpha_prev);
                end
            end
    cvx_end

    status_ok = strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved');

    if ~status_ok
        % fallback: keep previous alpha and continue
       % warning('Subproblem failed — keeping previous alpha');
        sca_iter = sca_iter - 1;
        continue;
    end
    
 
    

    % ---------- successful solve → convergence check & update ----------
    current_obj = double(t);
    obj_change   = abs(current_obj - obj_prev) / (abs(obj_prev) + 1);
    alpha_change = norm(double(vecAlpha) - alpha_prev) / (norm(alpha_prev) + 1);

    if obj_change < tol && alpha_change < tol        
        alpha_prev = double(vecAlpha);
        break;
    end

    alpha_prev = double(vecAlpha);
    obj_prev   = current_obj;
    sca_iter   = sca_iter - 1;
end

    
       
    % ========== FINAL SAFETY & REPORTING ==========
    status_ok = strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved');
    
     if status_ok
        alpha    = double(vecAlpha).';
        Ck_out   = double(Ck);       
    
        xi_val = double(xi);
        feasible_flag = max(xi_val) < 1e-4;
    else
        fprintf('WARNING: Final solve failed.\n');
        alpha    = alpha_prev.';
        Ck_out   = zeros(K,1);
        xi_val = 10;
        feasible_flag = false;
    end 
   
end
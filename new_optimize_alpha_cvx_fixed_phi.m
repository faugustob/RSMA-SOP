function [alpha,Ck_out] = new_optimize_alpha_cvx_fixed_phi(Rmin,alpha_prev,L_node,E_node,phi_St, phi_Sr, zeta_k_St, ...
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
lambda_penalty = 1e3;
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

Rmin_orig = Rmin;
Rmin_current = Rmin;
reduction_factor = 0.95;

%% ========================= SCA LOOP (safe retry version) =========================
sca_iter = max_SCA;
adapt_count_rmin = 0;               % renamed for clarity
max_adapt_rmin    = 500;

alpha_lb = 1e-8;                    % initial lower bound
adapt_count_lb   = 0;
max_adapt_lb     = 50;               % usually 3–6 is more than enough
lb_reduction     = 0.1;             % aggressive: 1e-8 → 1e-9 → 1e-10 → ...

adapted_lb_once = false;

while sca_iter > 0
    cvx_clear;
    cvx_begin quiet
        cvx_solver mosek
        variable vecAlpha(K+1) nonnegative
        variable Rc nonnegative
        variable Ck(K) nonnegative
        variable Rk(K) nonnegative
        variable t

        penalty_term = 0;
        maximize( t - lambda_penalty * penalty_term )

        subject to
            sum(vecAlpha) <= 1;
            vecAlpha >= 1e-8;

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
                Ck(k) + Rk(k) >= Rmin_current;
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
        % fprintf('SCA iter %d failed (status: %s)\n', sca_iter, cvx_status);

        if adapt_count_rmin < max_adapt_rmin
            % Stage 1: reduce Rmin (your existing logic)
            Rmin_current = max(0, Rmin_current * reduction_factor);
            adapt_count_rmin = adapt_count_rmin + 1;
           % fprintf('  → Rmin lowered to %.6f (attempt %d/%d)\n', ...
                    %Rmin_current, adapt_count_rmin, max_adapt_rmin);
            continue;
        else
            % Stage 2: Rmin already very low → try relaxing vecAlpha lower bound
            if adapt_count_lb < max_adapt_lb
                old_lb = alpha_lb;
                alpha_lb = alpha_lb * lb_reduction;
                adapt_count_lb = adapt_count_lb + 1;
                adapted_lb_once = true;

                %fprintf('  → vecAlpha lower bound relaxed: %.2e → %.2e (attempt %d/%d)\n', ...
                       % old_lb, alpha_lb, adapt_count_lb, max_adapt_lb);
                continue;                   % retry same sca_iter slot
            else
                fprintf('  Max adaptations (Rmin + alpha_lb) reached. Giving up.\n');
                break;
            end
        end
    end

    % ---------- successful solve → convergence check & update ----------
    current_obj = double(t);
    obj_change   = abs(current_obj - obj_prev) / (abs(obj_prev) + 1);
    alpha_change = norm(double(vecAlpha) - alpha_prev) / (norm(alpha_prev) + 1);

    if obj_change < tol && alpha_change < tol
        fprintf('Converged at iter %d (Rmin=%.6f, alpha_lb=%.2e)\n', ...
                sca_iter, Rmin_current, alpha_lb);
        alpha_prev = double(vecAlpha);
        break;
    end

    alpha_prev = double(vecAlpha);
    obj_prev   = current_obj;
    sca_iter   = sca_iter - 1;
end
    
    % ========== FINAL SAFETY & REPORTING ==========
    status_ok = strcmp(cvx_status,'Solved') || strcmp(cvx_status,'Inaccurate/Solved');
    
    if ~status_ok
        fprintf('WARNING: Final solve still failed after all adaptations.\n');
        alpha    = alpha_prev.';           % last known good point
        Ck_out   = zeros(K,1);
    else
        alpha    = double(vecAlpha).';
        Ck_out   = double(Ck);
    end
    
    % Summary print
    if Rmin_current < Rmin_orig || adapted_lb_once
        fprintf('Adapted solution: Rmin = %.6f (orig %.6f), alpha_lb = %.2e\n', ...
                Rmin_current, Rmin_orig, alpha_lb);
        if ~status_ok
            fprintf('  (but final solve was infeasible → using previous feasible alpha)\n');
        end
    else
        %fprintf('Solution with original Rmin and default alpha lower bound.\n');
    end
end
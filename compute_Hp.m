function Hp = compute_Hp(tau, nu, M, N, T, Deltaf, method)
% compute_Hp
% Builds the MN-by-MN OTFS delay–Doppler channel matrix Hp
%
% USAGE:
%   Hp = compute_Hp(tau, nu, M, N, T, Deltaf, 'loop');
%   Hp = compute_Hp(tau, nu, M, N, T, Deltaf, 'blocked');
%
% INPUTS:
%   tau, nu   : delay and Doppler (tau_{p,q}^{r,j}, nu_{p,q}^{r,j})
%   M, N      : grid sizes
%   T         : symbol duration
%   Deltaf    : subcarrier spacing
%   method    : 'loop' or 'blocked'
%
% OUTPUT:
%   Hp        : complex MN-by-MN matrix

if nargin < 7
    method = 'loop';
end

switch lower(method)
    case 'loop'
        Hp = Hp_loop(tau, nu, M, N, T, Deltaf);
    case 'blocked'
        Hp = Hp_blocked(tau, nu, M, N, T, Deltaf);
    otherwise
        error('Unknown method. Use ''loop'' or ''blocked''.');
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) CLEAR NESTED-LOOP IMPLEMENTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Hp = Hp_loop(tau, nu, M, N, T, Deltaf)

MN = M*N;
Hp = complex(zeros(MN, MN));

n = 0:N-1;
m = 0:M-1;
mp = 0:M-1;

j2pi = 1j*2*pi;
common_phase = exp(-j2pi*nu*tau);

[Mmat, MpMat] = meshgrid(m, mp);
delta = MpMat - Mmat;

for kp = 0:N-1
    for k = 0:N-1

        Sn = sum(exp(-j2pi*n*((kp-k)/N - nu/Deltaf)));

        for lp = 0:M-1
            for l = 0:M-1

                r = lp + kp*M + 1;
                c = l  + k *M + 1;

                % ---- h_{P,1} ----
                pref1 = common_phase * (1/M) * (1 - tau/T) * (1/N) * Sn;

                exp1 = exp(1j*pi*((tau/T+1)*(delta + nu/Deltaf)));
                phase1 = exp(1j*2*pi*((Mmat*lp)/M - (MpMat*l)/M - MpMat*(tau/T)));
                sinc1 = sinc((delta + nu/Deltaf)*(1 - tau/T));

                hP1 = pref1 * sum(exp1 .* phase1 .* sinc1, 'all');

                % ---- h_{P,2} ----
                pref2 = common_phase * exp(-j2pi*k/N) * ...
                        (1/M) * (tau/T) * (1/N) * Sn;

                phase2 = exp(1j*2*pi*((Mmat*lp)/M - (MpMat*l)/M - MpMat*(tau/T)));
                exp2 = exp(1j*pi*(tau/T)*(delta + nu/Deltaf));
                sinc2 = sinc((delta + nu/Deltaf)*(tau/T));

                hP2 = pref2 * sum(phase2 .* exp2 .* sinc2, 'all');

                Hp(r,c) = hP1 + hP2;
            end
        end
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) BLOCK / VECTORIZED IMPLEMENTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Hp = Hp_blocked(tau, nu, M, N, T, Deltaf)

MN = M*N;
Hp = complex(zeros(MN, MN));

n = 0:N-1;
m = 0:M-1;
mp = 0:M-1;

j2pi = 1j*2*pi;
common_phase = exp(-j2pi*nu*tau);

[Mmat, MpMat] = meshgrid(m, mp);
Mmat  = Mmat.';
MpMat = MpMat.';
delta = MpMat - Mmat;

A1 = exp(1j*pi*((tau/T+1)*(delta + nu/Deltaf))) .* ...
     sinc((delta + nu/Deltaf)*(1 - tau/T));

A2 = exp(1j*pi*(tau/T)*(delta + nu/Deltaf)) .* ...
     sinc((delta + nu/Deltaf)*(tau/T));

exp_mptau = exp(-j2pi*(mp.'*(tau/T)));

for kp = 0:N-1
    for k = 0:N-1

        Sn = sum(exp(-j2pi*n*((kp-k)/N - nu/Deltaf)));

        pref1 = common_phase * (1/M)*(1 - tau/T)*(1/N)*Sn;
        pref2 = common_phase * exp(-j2pi*k/N) * ...
                (1/M)*(tau/T)*(1/N)*Sn;

        lp = 0:M-1;
        l  = (0:M-1).';

        Em_lp  = exp(1j*2*pi*(m.'*lp/M));
        Emp_l  = exp(-1j*2*pi*(mp.'*l.'/M));

        C1 = A1.' * Em_lp;
        C2 = A2.' * Em_lp;

        C1 = C1 .* (exp_mptau * ones(1,M));
        C2 = C2 .* (exp_mptau * ones(1,M));

        Block1 = C1.' * Emp_l;
        Block2 = C2.' * Emp_l;

        rs = kp*M + 1;
        cs = k *M + 1;

        Hp(rs:rs+M-1, cs:cs+M-1) = pref1*Block1 + pref2*Block2;
    end
end

end

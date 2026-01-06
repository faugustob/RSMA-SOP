function Mk = nakagami_moment(k,m,Omega)
% k-th raw moment of a Nakagami-m random variable
% k     : moment order (k > -2m)
% m     : shape parameter (m >= 0.5)
% Omega : spread parameter

Mk = (Omega./m).^(k./2) .* gamma(m + k./2) ./ gamma(m);
end

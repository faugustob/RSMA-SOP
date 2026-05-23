function output = far_field_PL(G_t, G_r, G, d_x, d_y, lambda, F, theta_t, phi_t, theta_r, phi_r, A, d_1, d_2, Ny, Nx, phi_nm)
    n_y = (1-Ny/2 : Ny/2)'; % [Ny x 1]
    n_x = 1-Nx/2 : Nx/2; % [1 x Nx]
    % note: Ny and Nx assumed to be even numbers
    Beta = exp((1j*2*pi.*((sin(theta_t).*cos(phi_t)+sin(theta_r).*cos(phi_r)).*(n_x-1/2).*d_x+(sin(theta_t).*sin(phi_t)+sin(theta_r).*sin(phi_r)).*(n_y-1/2).*d_y+(lambda.*phi_nm./2./pi)))./lambda);
    Beta = squeeze(sum(sum(Beta, 1), 2))'; % [1 x M] vector
    output = G_t.*G_r.*G.*d_x.*d_y.*lambda.^2.*F(theta_t,phi_t).*F(theta_r,phi_r).*A.^2./(64*pi^3.*d_1^2.*d_2^2).*abs(Beta).^2; % [1 x M] vector
end
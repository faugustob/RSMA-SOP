function output = near_field_PL(G_t, G_r, G, d_x, d_y, lambda, A, F_tx, F, F_rx, theta_tx_nm, phi_tx_nm, theta_t_nm, phi_t_nm, theta_r_nm, phi_r_nm, theta_rx_nm, phi_rx_nm, r_t_nm, r_r_nm, phi_nm)
    F_combine_nm = @(theta_tx_nm, phi_tx_nm, theta_t_nm, phi_t_nm, theta_r_nm, phi_r_nm, theta_rx_nm, phi_rx_nm) F_tx(theta_tx_nm,phi_tx_nm).*F(theta_t_nm,phi_t_nm).*F(theta_r_nm,phi_r_nm).*F_rx(theta_rx_nm,phi_rx_nm);
    Beta = sqrt(F_combine_nm(theta_tx_nm, phi_tx_nm, theta_t_nm, phi_t_nm, theta_r_nm, phi_r_nm, theta_rx_nm, phi_rx_nm)).*exp(-1j.*(2*pi.*(r_t_nm+r_r_nm)-lambda.*phi_nm)./lambda)./(r_t_nm.*r_r_nm);
    Beta = squeeze(sum(sum(Beta, 1), 2))'; % [1 x M] vector
    output = G_t.*G_r.*G.*d_x.*d_y.*lambda.^2.*A.^2./(64*pi^3).*abs(Beta).^2; % [1 x M] vector
end
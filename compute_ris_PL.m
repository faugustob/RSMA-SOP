function path_loss = compute_ris_PL(lambda, N_V, N_H, S_xyz, User_loc, R_xyz,normal,F,F_tx,F_rx,G,G_t,G_r)

%% Basic parameters
d_x = lambda/2;
d_y = lambda/2;

RIS_size_x = N_H * d_x;
RIS_size_y = N_V * d_y;




%% RIS element positions
rn = zeros(3, N_V, N_H);
for j = 1:N_V
    for i = 1:N_H
        x = -RIS_size_x/2 + (i-0.5)*d_x;
        y = -RIS_size_y/2 + (j-0.5)*d_y;
        rn(:,j,i) = [x; y; 0];
    end
end

%% Per-element distances and angles
r_t_nm = zeros(N_V,N_H);
r_r_nm = zeros(N_V,N_H);

theta_t_nm  = zeros(N_V,N_H);
phi_t_nm    = zeros(N_V,N_H);
theta_r_nm  = zeros(N_V,N_H);
phi_r_nm    = zeros(N_V,N_H);
theta_tx_nm = zeros(N_V,N_H);
theta_rx_nm = zeros(N_V,N_H);

for j = 1:N_V
    for i = 1:N_H
        elem_pos = rn(:,j,i) + R_xyz;

        v_t = S_xyz - elem_pos;      % TX → RIS element
        v_r = User_loc - elem_pos;   % RIS element → RX

        r_t_nm(j,i) = norm(v_t);
        r_r_nm(j,i) = norm(v_r);

        theta_t_nm(j,i) = acos(dot(v_t,normal)/r_t_nm(j,i)); % angle between vector and x
        phi_t_nm(j,i)   = atan2(v_t(2), v_t(1)); % azimuth

        theta_r_nm(j,i) = acos(dot(v_r,normal)/r_r_nm(j,i));
        phi_r_nm(j,i)   = atan2(v_r(2), v_r(1)); % azimuth

        theta_tx_nm(j,i) = acos( dot(v_t, (S_xyz - R_xyz)) / ( norm(v_t) * norm(S_xyz - R_xyz) ));       
        theta_rx_nm(j,i) = acos(dot(v_r,(User_loc-R_xyz))/(norm(v_r)*norm((User_loc-R_xyz))));
    end
end


%% Distances TX/RX to RIS center
% d_1 = norm(R_xyz - S_xyz);
% d_2 = norm(User_loc - R_xyz);

%% Near / far-field threshold
% D = max(RIS_size_x, RIS_size_y);
% near_field_threshold = 2*D^2/lambda;


%% Reflection amplitude


   % Far-field using general expression (Eq. 3)
    path_loss = PL_eq3( ...
    F,F_tx,F_rx,G_t, G_r, G, d_x, d_y, lambda, ...
    theta_tx_nm, theta_t_nm, ...
    theta_r_nm, theta_rx_nm, r_t_nm, r_r_nm );
    x=1;

end
function output = PL_eq3( ...
    F,F_tx,F_rx,G_t, G_r, G, d_x, d_y, lambda, ...
    theta_tx_nm, theta_t_nm, ...
    theta_r_nm, theta_rx_nm, r_t_nm, r_r_nm )

    % %% Radiation patterns (Tang paper)
    % F      = @(theta) cos(theta).^3;
    % F_tx   = @(theta) cos(theta).^62;
    % F_rx   = @(theta) cos(theta).^62;

    % Combined radiation pattern (Eq. 3)
    F_nm = F_tx(theta_tx_nm,0) ...
         .* F(theta_t_nm,0) ...
         .* F(theta_r_nm,0) ...
         .* F_rx(theta_rx_nm,0);

    % Field contribution per RIS element
    E_nm = sqrt(F_nm) ...
        .* exp(-1j * (2*pi.* (r_t_nm + r_r_nm)/lambda) ) ...
        ./ (r_t_nm .* r_r_nm);

    % Coherent sum
    E_total = sum(E_nm(:));

    % Path loss
    output = (G_t .* G_r .* G ...
           .* d_x .* d_y .* lambda.^2 ...
           ./ (64*pi^3)) ...
           .* abs(E_total).^2;
end

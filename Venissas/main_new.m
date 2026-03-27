 
clear; clc; close all

%%%%%%%%%%%%%%%%%%%%%% (NLoS) (SR) (SCA Algorithm) %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SearchAgents_no = 2;
Function_name = 'F_RIS_NLoS_SR'; 

Max_iteration = 1000; % Maximum number of iterations
Total_runs = 50; % Number of times to repeat the experiment


Elements = 50;        %%%%%%%%%%%%%%%%%%%%% Number of elements in the mirror array

%%%%%%%Load details of the selected benchmark function
lb=[repelem(-90,Elements),repelem(-90,Elements)];
ub=[repelem( 90,Elements),repelem( 90,Elements)];
dim=2*Elements; % [ First 2k is for the oriantation angles + power allocation coefficient of NOMA]
 
Best_scores = zeros(1, Total_runs);
Best_positions = zeros(Total_runs, dim);
Secrecy_Rates = zeros(1, Total_runs);
trace_PT_sums = zeros(1, Total_runs);
trace_hp_sums = zeros(1, Total_runs);



for run = 1:Total_runs
    [Best_score,Best_pos,cg_curve_RIS_NLoS_SR,t,early_stop] = SCA(SearchAgents_no,Max_iteration,lb,ub,dim);
    Best_scores(run) = Best_score;
    Best_positions(run, :) = Best_pos;
    Secrecy_Rates(run) = max(-1*Best_score,0);
end

Average_Secrecy_Rate = mean(Secrecy_Rates);
Max_Secrecy_Rate = max(Secrecy_Rates);




function [Destination_fitness,Destination_position,Convergence_curve,Objective_values,t,X,early_stop]=SCA(N,Max_iteration,lb,ub,dim)

% disp('SCA is optimizing your problem');

%Initialize the set of random solutions
X=initialization(N,dim,ub,lb);

Destination_position=zeros(1,dim);
Destination_fitness=inf;

Convergence_curve=zeros(1,Max_iteration);
Objective_values = zeros(1,size(X,1));

% Parameters for early stopping
tolerance = 1e-6;       % Minimum improvement threshold
stall_limit = 100;      % Number of iterations allowed without improvement
stall_counter = 0;      
best_val_old = inf;
early_stop = false;     % Initialize early stop flag

% Calculate the fitness of the first set and find the best one
for i=1:size(X,1) %i=1:30
    Objective_values(1,i)=F_RIS_NLoS_SR(X(i,:));
    if i==1
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    elseif Objective_values(1,i)<Destination_fitness
        Destination_position=X(i,:);
        Destination_fitness=Objective_values(1,i);
    end
    
    All_objective_values(1,i)=Objective_values(1,i);
end

%Main loop
t=2; % start from the second iteration since the first iteration was dedicated to calculating the fitness
while t<=Max_iteration
    
    % Eq. (3.4)
    a = 2;
    Max_iteration = Max_iteration;
    r1=a-t*((a)/Max_iteration); % r1 decreases linearly from a to 0
    
    % Update the position of solutions with respect to destination
    for i=1:size(X,1) % in i-th solution
        for j=1:size(X,2) % in j-th dimension
            
            % Update r2, r3, and r4 for Eq. (3.3)
            r2=(2*pi)*rand();
            r3=2*rand;
            r4=rand();
            
            % Eq. (3.3)
            if r4<0.5
                % Eq. (3.1)
                X(i,j)= X(i,j)+(r1*sin(r2)*abs(r3*Destination_position(j)-X(i,j)));
            else
                % Eq. (3.2)
                X(i,j)= X(i,j)+(r1*cos(r2)*abs(r3*Destination_position(j)-X(i,j)));
            end
            
        end
    end
    
    for i=1:size(X,1)
         
        % Check if solutions go outside the search spaceand bring them back
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        % Calculate the objective values
        Objective_values(1,i)=F_RIS_NLoS_SR(X(i,:));
        
        % Update the destination if there is a better solution
        if Objective_values(1,i)<Destination_fitness
            Destination_position=X(i,:);
            Destination_fitness=Objective_values(1,i);
        end
    end
    
    Convergence_curve(t)=Destination_fitness;
    
    % Display the iteration and best optimum obtained so far
    if mod(t,1)==0
        display(['At iteration ', num2str(t), ' the optimum is ', num2str(-1*Destination_fitness)]);
    end
    
    %Increase the iteration counter

     % ---------------- Early stopping check ----------------
    % ---------------- Early stopping check ----------------
    if abs(Destination_fitness - best_val_old) < tolerance
        stall_counter = stall_counter + 1;
    else
        stall_counter = 0;
    end
    best_val_old = Destination_fitness;
    
    if stall_counter >= stall_limit
        early_stop = true;
        Convergence_curve = Convergence_curve(1:t); % trim unused values
        break;
    end
    % ------------------------------------------------------
    
    t = t + 1;
end
end

function o = F_RIS_NLoS_SR(x)
% VLC parameters

FOV=85;                            %Half-angle FOV
FOVr=(FOV*pi)/180;                 %FOV of a receiver
index=1.5;                         %refractive index of a lens at a PD; ignore if no lens is used
G_Con=(index^2)/(sin(FOVr))^2;     %gain of an optical concentrator
rho_irs=0.95;                      %reflection coeeficient of IRS
lx=4;lz=1;                         %wall/IRS dimension of interest in x-z plane (meters)
w_m=0.1;                           %size of each smaller mirror (omar:Width) 
h_m=0.1;                           %size of each smaller mirror (Omar:Height)
Nx=lx/w_m; Nz=lz/h_m;              %number of grid in each surface 
dA=lx*lz/(Nx*Nz);                  %calculation of grid area
k = 50;                              %%%%%%%%%%%%%%%%%%%%% Number of elements in the mirror array
w_x=1.1:w_m:1.5;                     %%%%%%%% 50 elements w_x=1:w_m:1.5; // 300 elements w_x=1:w_m:4; // 400 elements w_x=1:w_m:5;
w_z=1.1:w_m:2;
theta_sem_ang=70;                  %semi-angle at half power
m_ue=-log10(2)/log10(cosd(theta_sem_ang)); %Lambertian order of emission
Adet=1e-4;                         %detector physical area of a PD
P_opt=4;                         %%%%%%%%%%%%%%%%%%%%% Changes in the simulations: 2:2:14 [10] // 0.5:0.5:2.5 [20] // 1:0.5:4 [17]
rho=0.53;                          %photodetector responsivity (Amp/Watt)
q=3;
mu=41;                         %mean in degrees
sigma=9;                       %standard deviation in degrees
m=1;n=1;                       %The dimension of y
y=laprnd(m,n,mu,sigma);        %polar angle pola Laplace distribution, mean and standard deviation of 41 and 9,
pola=y;                        %uniform distribution [-pi,pi]
azi=unifrnd(0,180);            %azimuth angle azi

TP1 = [0 0 3];  % Transmitter position
RP = [0 0 2];  % Receiver position
EVE = [0 0 1+rand()];  % Eve position

    %% Pathloss channel gain, h_PL
wavelength = 450e-9; % Wavelength in meters
sigma_lambda = 0.056; % Total attenuation coefficient (Np/m) for the given wavelength
l = TP1(3)-RP(3);              % Link distance (m)
h_PL = exp(-sigma_lambda * l);

%% Oceanic turbulence channel gain, h_PT

omega = 0.721;
lambda = 0.1479;
a = 0.0121;
b = 7.4189;	
c = 65.6983;
r=1;
AverageSNR = 1;
eta = 0.9;

%% Parameter initialization
sigma_rd = 0.5;              % Pointing error jitter standard deviation
w_L = 4;                     % Beam waist at Rx
Dr = 0.075;                      % Receiver aperture diameter
nu = (Dr * sqrt(pi)) / (2 * sqrt(2) * w_L);
A0 = (erf(nu))^2;            % Maximum power fraction
w_Leq = w_L * sqrt((sqrt(pi) * erf(nu)) / (2 * nu * exp(-nu^2)));  % Eq (24)


% a_l = 0.02;               % Detector radius (m)
% wL = 4;              % Beam waist at receiver (m)
% sigma_s = 0.8;         % Pointing error std dev (m)

%%
    inc_ang_MA=zeros(length(w_x)-1,length(w_z)-1); %(Omar: MA= mirror array)

% %%% Calculate the channel gain for H_1
H_1_NLoS_dummy=zeros(length(w_x)-1,length(w_z)-1);
H_1_NLoS=0;

H_Eve_1_NLoS_dummy=zeros(length(w_x)-1,length(w_z)-1);
H_Eve_1_NLoS=0;

trace_PT = zeros(length(w_x)-1,length(w_z)-1);
trace_hp = zeros(length(w_x)-1,length(w_z)-1);

% Preallocations
num_x = length(w_x) - 1;
num_z = length(w_z) - 1;

% ====================== MAIN LOOP ======================

    for ii=1:num_x
                for jj=1:num_z
                    WP1 = [w_x(ii)+w_m/2 2 w_z(jj)+h_m/2]; % point of incidence in MA (Omar: WP1=IRS/wall position) 
                    D1=sqrt(dot(TP1-WP1,TP1-WP1)); % distance from transmitter to (WP1==MA)
                    cos_phi= abs(WP1(3)-TP1(3))/D1; % irradiance angle from Tx
                    cos_alpha_irs = abs(TP1(2)-WP1(2))/D1; % incidence angle on MA  (compare both)
                    D2=sqrt(dot(WP1-RP,WP1-RP)); % distance from MA (ii,jj) to receiver
                    cos_psi_Or_irs_U1=abs(((WP1(1)-RP(1))/D2)*sind(pola)*sind(azi) + ((WP1(2)-RP(2))/D2)*sind(pola)*sind(azi) + ((WP1(3)-RP(3))/D2)*cosd(pola));% with orientation CHECK
                    % phi_U1=acosd(cos_psi_Or_irs_U1); %angle of incidence from MA to receiver phi == inc_ang_MA(ii,jj)  [we need one only] 
                    
                    % % Ensemble mean of mEGG turbulence
                    % I_bar = omega*lambda + (1-omega)*b*gamma(a + 1/c)/gamma(a);

                     % === Turbulence and pointing error for this realization ===
                    is_exp = rand < omega; % Boolean selector
                    f_lambda = exprnd(lambda, 1);  
                    gamma_samples = gamrnd(a, 1);  % Match the first code
                    g_alpha = (gamma_samples.^(1/c)) * b; % Apply scaling correctly           
                    I_samples = is_exp .* f_lambda + (~is_exp) .* g_alpha; % Correct mixture sampling
                    
                    I = I_samples;
                    % Compute mu_r based on detection technique
                    if r == 1
                        mu_r = AverageSNR;  % For heterodyne detection
                    elseif r == 2
                        mu_r = AverageSNR / (2 * omega * lambda^2 + b^2 * (1 - omega) * gamma(a + 2/c) / gamma(a)); 
                    else
                        fprintf('r should be either 1 for heterodyne detection or 2 for IM/DD\n');
                        return; 
                    end
        
                    % Transform I to SNR_UWOC

                    h_PT = (I) * mu_r;
                    trace_PT(ii,jj) = h_PT;


                    %Pointing error             
                    rd_samples = sigma_rd * sqrt(-2*log(rand(1, 1)));  % Rayleigh samples
                    hp = A0 * exp(-2 * rd_samples.^2 / w_Leq^2); % Eq (23)
                    trace_hp(ii,jj) = hp;

         
                    % === Eavesdropper ===
                    inc_ang_MA(ii,jj)=acosd(cos_psi_Or_irs_U1);
                    if abs(acosd(cos_psi_Or_irs_U1))<=FOV % Omar: In the next line we calculate H_RIS
                       H_1_NLoS_dummy(ii,jj)=h_PL*hp*h_PT*(m_ue+1)*Adet*rho_irs*dA*...
                       (cos_phi^m_ue)*cos_alpha_irs...
                       *((abs(WP1(1)-RP(1))/D2)*sind(x(ii))*cosd(x(k+jj)) + (abs(WP1(2)-RP(2))/D2)*cosd(x(ii))*cosd(x(k+jj)) + (abs(WP1(3)-RP(3))/D2)*sind(x(k+jj)))*... --> irradiance angle from the k-th reflecting element towards the intended user
                       G_Con*cos_psi_Or_irs_U1/(2*pi^2*D1^2*D2^2);
                        %H_1_NLoS_dummy(ii,jj)=rho_irs*(m_ue+1)*Adet*G_Con*(cos_phi^m_ue)*cos_psi_Or_irs_U1/(2*pi^2*(D1+D2)^2);
                    end

                     RP_EVE = EVE; % receiver position vector (EVE);
                     D2_EVE = sqrt(dot(WP1-RP_EVE,WP1-RP_EVE)); % distance from MA (ii,jj) to receiver
                     cos_psi_Or_irs_Eve = abs(((WP1(1)-RP_EVE(1))/D2_EVE)*sind(pola)*sind(azi) + ((WP1(2)-RP_EVE(2))/D2_EVE)*sind(pola)*sind(azi) + ((WP1(3)-RP_EVE(3))/D2_EVE)*cosd(pola));% with orientation CHECK
                    
                    if abs(acosd(cos_psi_Or_irs_Eve))<=FOV % Omar: In the next line we calculate H_RIS
                       H_Eve_1_NLoS_dummy(ii,jj)=h_PL*hp*h_PT*(m_ue+1)*Adet*rho_irs*dA*...
                       (cos_phi^m_ue)*cos_alpha_irs...
                       *((abs(WP1(1)-RP_EVE(1))/D2_EVE)*sind(x(ii))*cosd(x(k+jj)) + (abs(WP1(2)-RP_EVE(2))/D2_EVE)*cosd(x(ii))*cosd(x(k+jj)) + (abs(WP1(3)-RP_EVE(3))/D2_EVE)*sind(x(k+jj)))*... --> irradiance angle from the k-th reflecting element towards the intended user
                       G_Con*cos_psi_Or_irs_Eve/(2*pi^2*D1^2*D2_EVE^2);
                       %H_Eve_1_NLoS_dummy(ii,jj)=rho_irs*(m_ue+1)*Adet*G_Con*(cos_phi^m_ue)*cos_psi_Or_irs_Eve1/(2*pi^2*(D1+D2)^2);
                    end

                end
    end

save('channel_logs_1.mat', 'h_PT', 'hp')
save('channel_logs.mat', 'trace_PT', 'trace_hp')

H_1_NLoS=    ((eta)^r) * (sum(H_1_NLoS_dummy(:)))^r;
H_Eve_1_NLoS=sum(H_Eve_1_NLoS_dummy(:));


%%  

%%%Calculate the sumR
BW=200e6;                        %system bandwidth
No_vlc=(BW*(10^(-21)));          %noise power 
P_Signal=(P_opt/q)^2;            % Electrical transmit power

R1=0;
R1_Eve=0;

R1 = BW*log2(1+(exp(1)/(2*pi))*(P_Signal*(rho*H_1_NLoS).^2)/(((rho*H_1_NLoS).^2*P_Signal)+No_vlc)); 
R1_Eve = BW*log2(1+(exp(1)/(2*pi))*(P_Signal*(rho*H_Eve_1_NLoS).^2)/(((rho*H_Eve_1_NLoS).^2*P_Signal)+No_vlc)); 


o = R1-R1_Eve;

end

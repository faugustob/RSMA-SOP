% Matlab Simulation of the Zak Transform for ECE 4806 - Team 14 (AFRL SDR University Challenge) Fall'23 - Spring'24

%% This code simulates a 3x5 matrix through the Discrete Zak Transform (DZT) at the Baseband level. 
% The code will take a complex matrix at the input, perform the Inverse
% Fast-Fourier transforom (IFFT) of the rows and then seralize the matrix so that
% it can be sent throught the channel. In the channel the data (transmitted
% symbols) are filtered by the doppler-shift before the channel estimation
% is completed. The doppler-shift was obtained by taking an impulse response 
% and shifted it by the doppler frequency (The calculations are shown in
% the 'Parameters' section below). The filtered data (received signals) 
% is parallelized back to the shape of a matrix. The Fast-Fourier Transform (FFT) is taken on
% the rows and then observed. The sections are arranged below and the plots will be generated 
% after running the code. 
%% Parameters and Delay-Doppler Calculations
clc; clear;
% Carrier Frequency
fc = 2.4e9; % In LABview fc = 915MHz for the USRP 2901
% Sampling frequency
Fs = 1e6;
% Impulse response 
h = [0 0 0 0.3 0 0 0 1.5 0 0];
% Time vector
t = 0:length(h)-1;
% Doppler shift frequency calculation
v = 10; % velocity in m/s
lambda = physconst('LightSpeed')/fc;
fd = speed2dop(v,lambda); 

% Create doppler-shifted signal
doppler_signal = exp(-1j*2*pi*fd*t);

% Apply the Doppler shift
h_doppler_shift = h.*doppler_signal;
% Convert to freq domain
H_doppler = fft(h_doppler_shift); 

% Apply the Delay the Doppler shift
h_delayed = [zeros(1, 2), h_doppler_shift(1:end-2)];

%% Matrix Definition & DZT
% Create MxN (3x5) matrix
matrix_A = zeros(3,5);
% Find the center of the matrix
center_row = ceil(3/2); 
center_col = ceil(5/2);
% Make the center symbol complex
matrix_A(center_row,center_col) = 5+2j; 
matrix_A; % display A
% Take the inverse fast-fourier transform of each row
ifft_A = ifft(matrix_A,[], 2);

% P/S Conversion
serial_data = reshape(ifft_A.', 1,[]);
% Simualte Channel transmission (Tx)
tx_signal = serial_data;

% Received signal (Rx)
h_doppler_shift = [0 1 zeros(1,8)];
rx_signal = filter(h_doppler_shift,1,tx_signal);

% Find the peaks of the rx_signal
[~, peaks] = findpeaks(abs(rx_signal));

% Initialize estimated impulse response
h_estimated = zeros(size(h)); 
% Use peaks to estimate the channel
h_estimated(peaks) = rx_signal(peaks); 

% S/P Conversion - demodulation
num_rows = size(ifft_A,1);
parallel_data = reshape(rx_signal, num_rows, []);
%parallel_data = parallel_data.';

% Take the fast-fourier tranform of each row
fft_A = fft(parallel_data,[],2);

%% Plots
% - Time Domain plots -
% Time domain plot of the impulse response
figure;
stem(t, real(h), 'b', 'LineWidth', 2)
title('Impulse Reponse in the Time Domain');
xlabel('Time Index (t)');
ylabel('Amplitude');
grid on;

% Plot of Tx and Rx in the Time Domain
figure;
stem(real(tx_signal), 'b', 'LineWidth', 1.5);
hold on;
stem(imag(tx_signal), 'r', 'LineWidth', 1.5);
hold off;
title('Transmitted Signal (tx) in the Time Domain');
xlabel('Sample Index');
ylabel('Amplitude');
legend('Real Part', 'Imaginary Part');
grid on; 

figure;
stem(real(rx_signal), 'b', 'LineWidth', 1.5);
hold on;
stem(imag(rx_signal), 'r', 'LineWidth', 1.5);
hold off;
title('Received Signal (rx) with a Delay of 1 in the Time Domain');
xlabel('Sample Index');
ylabel('Amplitude');
legend('Real Part', 'Imaginary Part');
grid on;

% Plot the Doppler Shift of the impulse response
figure;
stem(t, real(h_doppler_shift), 'r', 'LineWidth', 2);
hold on;
stem(t, imag(h_doppler_shift), 'bx', 'LineWidth',2);
hold off;
title('Doppler-Shifted Impulse Response in the Time Domain'); %doppler shift
xlabel('Time');
ylabel('Amplitude');
grid on;
legend('Real','Imag');

% - Delay-Doppler Domain plots -
% Plot of the original impulse response using bar3
zeros_grid = zeros(30,30);
figure;
bar3(zeros_grid); 
hold on;
bar3([real(h) ; imag(h)]);
hold off;
title('Impulse Response in the Delay-Doppler Domain');
xlabel('Delay');
ylabel('Doppler');
zlabel('Amplitude');
colorbar;
grid on;

% Plot of the matrix sent through the DZT using bar3
figure; 
bar3(abs(matrix_A));
xlabel('Delay');
ylabel('Doppler');
zlabel('Amplitude')
title('The plot of the matrix with a complex value in the center');
colorbar;

% Plot of the matrix after iift
figure; 
bar3(abs(ifft_A));
xlabel('Delay');
ylabel('Doppler');
zlabel('Amplitude')
title('IFFT of each Row in the Delay-Doppler Domain');
colorbar;

% Plot of tx_signal and rx_signal in the Delay-Doppler Domain
[TX,RX] = meshgrid(rx_signal,tx_signal);
figure;
% Transpose the Tx and Rx signal
tx_signal_transposed = tx_signal.';
rx_signal_transposed = rx_signal.';
% Transmitted signal
bar3(zeros_grid);
hold on;
bar3([real(tx_signal_transposed), imag(tx_signal_transposed)]);
hold off;
title('Transmitted Symbol in the Delay-Doppler Domain')
xlabel('Delay');
ylabel('Doppler');
zlabel('Amplitude');
colorbar;
% Received signal
figure;
bar3(zeros_grid);
hold on;
bar3([real(rx_signal_transposed), imag(rx_signal_transposed)]);
hold off;
title('Received Symbol in the Delay-Doppler Domain')
xlabel('Delay');
ylabel('Doppler');
zlabel('Amplitude');
colorbar;

% Plot the estimated impulse response in the Delay-Doppler Domain
figure;
% Delay max in seconds
delay_max = 100/Fs; 
% Doppler max
doppler_max = 1:300;
[X_delay, Y_doppler] = meshgrid(delay_max*(1:size(h_estimated,2)), doppler_max);

% Create matrix with the same number of rows as X_delay and Y_doppler
h_matrix = [real(h_estimated(peaks)) ; imag(h_estimated(peaks))]; % Transposed to correct the orientation
% Reshape h_matrix into a column vector
h_vector = reshape(h_matrix, [], 1);

% Reshape Y_doppler into a vector
Y_doppler_vector = repmat(Y_doppler(:), size(h_matrix, 1), 1);
Y_doppler_vector = Y_doppler_vector(1:length(h_vector));

bar3(zeros_grid); 
hold on;
bar3([[real(Y_doppler_vector); imag(Y_doppler_vector)], [real(h_vector); imag(h_vector)]]);
hold off;
title('Estimated Impulse Response in the Delay-Doppler Domain');
xlabel('Delay');
ylabel('Doppler');
zlabel('Amplitude');
colorbar;

% Plot of the matrix after fft
figure; 
bar3(abs(fft_A));
xlabel('Delay');
ylabel('Doppler');
zlabel('Amplitude')
title('FFT of each Row in the Delay-Doppler Domain');
colorbar;
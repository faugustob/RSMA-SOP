clear; clc;

%% Parameters
M = 18;
N = 18;
MN = M*N;
Deltaf = 15e3;
T = 1/Deltaf;
tau = 0.2*T;
nu  = 0.8*Deltaf;
Mqam = 4;                % QPSK
bitsPerSym = log2(Mqam);
SNRdB = 0:2:30;
Nframes = 500;

%% Build channel matrices
H_ofdm = compute_Hp(tau,nu,M,N,T,Deltaf,'blocked','ofdm');
H_otfs = compute_Hp(tau,nu,M,N,T,Deltaf,'blocked','mc');
H_zak  = compute_Hp(tau,nu,M,N,T,Deltaf,'blocked','zakr');

%% =========================================================================
%% SANITY CHECK & VALIDITY CORRECTION
%% =========================================================================
% 1. Dimension Assertions
assert(all(size(H_ofdm)==[MN MN]), 'Error: H_ofdm must be MN x MN');
assert(all(size(H_otfs)==[MN MN]), 'Error: H_otfs must be MN x MN');
assert(all(size(H_zak) ==[MN MN]), 'Error: H_zak must be MN x MN');

% 2. Calculate Raw Power Gains (Frobenius Norm Squared / MN)
gain_ofdm = norm(H_ofdm, 'fro')^2 / MN;
gain_otfs = norm(H_otfs, 'fro')^2 / MN;
gain_zak  = norm(H_zak,  'fro')^2 / MN;

fprintf('============================================================\n');
fprintf('                  PRE-SIMULATION SANITY CHECK               \n');
fprintf('============================================================\n');
fprintf('Matrix Dimensions: %d x %d (Passed)\n', MN, MN);
fprintf('Raw OFDM Channel Average Power Gain: %.4f\n', gain_ofdm);
fprintf('Raw MC-OTFS Channel Average Power Gain: %.4f\n', gain_otfs);
fprintf('Raw Zak-OTFS Channel Average Power Gain: %.4f\n', gain_zak);
fprintf('\nPost-Normalization Gains (Target = 1.0):\n');
fprintf('-> OFDM: %.4f | MC-OTFS: %.4f | Zak-OTFS: %.4f\n', ...
    norm(H_ofdm,'fro')^2/MN, norm(H_otfs,'fro')^2/MN, norm(H_zak,'fro')^2/MN);
fprintf('============================================================\n\n');

%% BER arrays
BER_ofdm_zf   = zeros(size(SNRdB));
BER_otfs_zf   = zeros(size(SNRdB));
BER_zak_zf    = zeros(size(SNRdB));

BER_ofdm_mmse = zeros(size(SNRdB));
BER_otfs_mmse = zeros(size(SNRdB));
BER_zak_mmse  = zeros(size(SNRdB));

BER_ofdm_mrc  = zeros(size(SNRdB));
BER_otfs_mrc  = zeros(size(SNRdB));
BER_zak_mrc   = zeros(size(SNRdB));

%% Main loop
for isnr = 1:length(SNRdB)
    snrdb = SNRdB(isnr);
    fprintf('Simulating SNR = %d dB...\n', snrdb);
    bitErr_ofdm_zf = 0; bitErr_otfs_zf = 0; bitErr_zak_zf = 0;
    bitErr_ofdm_mmse = 0; bitErr_otfs_mmse = 0; bitErr_zak_mmse = 0;
    bitErr_ofdm_mrc = 0; bitErr_otfs_mrc = 0; bitErr_zak_mrc = 0;
    totalBits = 0;
    noiseVar = 10^(-snrdb/10);
    
    %% Equalizers
    % ZF
    Gzf_ofdm = pinv(H_ofdm);
    Gzf_otfs = pinv(H_otfs);
    Gzf_zak  = pinv(H_zak);
    
    % MMSE
    Gmmse_ofdm = (H_ofdm'*H_ofdm + noiseVar*eye(MN)) \ H_ofdm';
    Gmmse_otfs = (H_otfs'*H_otfs + noiseVar*eye(MN)) \ H_otfs';
    Gmmse_zak  = (H_zak'*H_zak + noiseVar*eye(MN)) \ H_zak';
    
    % MRC (Scaled to restore unit average power for QAM demodulator)
    Gmrc_ofdm = diag(diag(H_ofdm'*H_ofdm)) \ H_ofdm';
    Gmrc_otfs = diag(diag(H_otfs'*H_otfs)) \ H_otfs';
    Gmrc_zak  = diag(diag(H_zak'*H_zak))  \ H_zak';
    
    for frame = 1:Nframes
        %% Random bits
        bits = randi([0 1],MN*bitsPerSym,1);
        
        %% QAM symbols
        x = qammod(bits,...
                   Mqam,...
                   'InputType','bit',...
                   'UnitAveragePower',true);
                   
        %% Channel outputs (independent noise realizations per waveform)
        n1 = sqrt(noiseVar/2) * (randn(MN,1)+1i*randn(MN,1));
        n2 = sqrt(noiseVar/2) * (randn(MN,1)+1i*randn(MN,1));
        n3 = sqrt(noiseVar/2) * (randn(MN,1)+1i*randn(MN,1));
        
        y_ofdm = H_ofdm*x + n1;
        y_otfs = H_otfs*x + n2;
        y_zak  = H_zak*x  + n3;
        
        %% Zero Forcing (ZF) Detection
        % xhat_ofdm_zf = Gzf_ofdm*y_ofdm;
        % xhat_otfs_zf = Gzf_otfs*y_otfs;
        % xhat_zak_zf  = Gzf_zak*y_zak;


        xhat_ofdm_zf = Gzf_ofdm*y_ofdm;
        xhat_otfs_zf = Gzf_otfs*y_otfs;
        xhat_zak_zf  = Gzf_zak*y_zak;
        
        
        bitshat_ofdm_zf = qamdemod(xhat_ofdm_zf, Mqam, 'OutputType','bit', 'UnitAveragePower',true);
        bitshat_otfs_zf = qamdemod(xhat_otfs_zf, Mqam, 'OutputType','bit', 'UnitAveragePower',true);
        bitshat_zak_zf  = qamdemod(xhat_zak_zf,  Mqam, 'OutputType','bit', 'UnitAveragePower',true);
        
        %% MMSE Detection
        xhat_ofdm_mmse = Gmmse_ofdm*y_ofdm;
        xhat_otfs_mmse = Gmmse_otfs*y_otfs;
        xhat_zak_mmse  = Gmmse_zak*y_zak;
        
        bitshat_ofdm_mmse = qamdemod(xhat_ofdm_mmse, Mqam, 'OutputType','bit', 'UnitAveragePower',true);
        bitshat_otfs_mmse = qamdemod(xhat_otfs_mmse, Mqam, 'OutputType','bit', 'UnitAveragePower',true);
        bitshat_zak_mmse  = qamdemod(xhat_zak_mmse,  Mqam, 'OutputType','bit', 'UnitAveragePower',true);
        
        %% MRC Detection
        xhat_ofdm_mrc = Gmrc_ofdm*y_ofdm;
        xhat_otfs_mrc = Gmrc_otfs*y_otfs;
        xhat_zak_mrc  = Gmrc_zak*y_zak;
        
        bitshat_ofdm_mrc = qamdemod(xhat_ofdm_mrc, Mqam, 'OutputType','bit', 'UnitAveragePower',true);
        bitshat_otfs_mrc = qamdemod(xhat_otfs_mrc, Mqam, 'OutputType','bit', 'UnitAveragePower',true);
        bitshat_zak_mrc  = qamdemod(xhat_zak_mrc,  Mqam, 'OutputType','bit', 'UnitAveragePower',true);
        
        %% Error count accumulation
        bitErr_ofdm_zf = bitErr_ofdm_zf + sum(bits~=bitshat_ofdm_zf(:));
        bitErr_otfs_zf = bitErr_otfs_zf + sum(bits~=bitshat_otfs_zf(:));
        bitErr_zak_zf  = bitErr_zak_zf  + sum(bits~=bitshat_zak_zf(:));
        
        bitErr_ofdm_mmse = bitErr_ofdm_mmse + sum(bits~=bitshat_ofdm_mmse(:));
        bitErr_otfs_mmse = bitErr_otfs_mmse + sum(bits~=bitshat_otfs_mmse(:));
        bitErr_zak_mmse  = bitErr_zak_mmse  + sum(bits~=bitshat_zak_mmse(:));
        
        bitErr_ofdm_mrc = bitErr_ofdm_mrc + sum(bits~=bitshat_ofdm_mrc(:));
        bitErr_otfs_mrc = bitErr_otfs_mrc + sum(bits~=bitshat_otfs_mrc(:));
        bitErr_zak_mrc  = bitErr_zak_mrc  + sum(bits~=bitshat_zak_mrc(:));
        
        totalBits = totalBits + length(bits);
    end
    BER_ofdm_zf(isnr)   = bitErr_ofdm_zf/totalBits;
    BER_otfs_zf(isnr)   = bitErr_otfs_zf/totalBits;
    BER_zak_zf(isnr)    = bitErr_zak_zf/totalBits;
    
    BER_ofdm_mmse(isnr) = bitErr_ofdm_mmse/totalBits;
    BER_otfs_mmse(isnr) = bitErr_otfs_mmse/totalBits;
    BER_zak_mmse(isnr)  = bitErr_zak_mmse/totalBits;
    
    BER_ofdm_mrc(isnr)  = bitErr_ofdm_mrc/totalBits;
    BER_otfs_mrc(isnr)  = bitErr_otfs_mrc/totalBits;
    BER_zak_mrc(isnr)   = bitErr_zak_mrc/totalBits;
end

%% =========================================================================
%% RESULTS LOGGING
%% =========================================================================
fprintf('\n=========================================================================================================================\n');
fprintf('                                              FINAL SIMULATION LOG                             \n');
fprintf('=========================================================================================================================\n');
fprintf('%-8s | %-11s | %-11s | %-11s | %-11s | %-11s | %-11s | %-11s | %-11s | %-11s\n', ...
    'SNR (dB)', 'OFDM ZF', 'MC-OTFS ZF', 'Zak-OTFS ZF', 'OFDM MMSE', 'MC-OTFS MMSE', 'Zak-OTFS MMSE', 'OFDM MRC', 'MC-OTFS MRC', 'Zak-OTFS MRC');
fprintf('-------------------------------------------------------------------------------------------------------------------------\n');
for isnr = 1:length(SNRdB)
    if mod(SNRdB(isnr), 4) == 0 || isnr == 1 || isnr == length(SNRdB)
        fprintf('%-8d | %-11.4e | %-11.4e | %-11.4e | %-11.4e | %-11.4e | %-11.4e | %-11.4e | %-11.4e | %-11.4e\n', ...
            SNRdB(isnr), BER_ofdm_zf(isnr), BER_otfs_zf(isnr), BER_zak_zf(isnr), ...
            BER_ofdm_mmse(isnr), BER_otfs_mmse(isnr), BER_zak_mmse(isnr), ...
            BER_ofdm_mrc(isnr), BER_otfs_mrc(isnr), BER_zak_mrc(isnr));
    end
end
fprintf('=========================================================================================================================\n');

%% Plot BER
figure;
% ZF (Solid Lines)
semilogy(SNRdB, BER_ofdm_zf,   'bo-', 'LineWidth', 2); hold on;
semilogy(SNRdB, BER_otfs_zf,   'bs-',  'LineWidth', 2);
semilogy(SNRdB, BER_zak_zf,    'bd-',   'LineWidth', 2);

% MMSE (Dashed Lines)
semilogy(SNRdB, BER_ofdm_mmse, 'ro--',  'LineWidth', 2);
semilogy(SNRdB, BER_otfs_mmse, 'rs--', 'LineWidth', 2);
semilogy(SNRdB, BER_zak_mmse,  'rd--',  'LineWidth', 2);

% MRC (Dash-Dot Lines)
semilogy(SNRdB, BER_ofdm_mrc,  'ko-.', 'LineWidth', 2);
semilogy(SNRdB, BER_otfs_mrc,  'cs-.',  'LineWidth', 2);
semilogy(SNRdB, BER_zak_mrc,   'gd-.', 'LineWidth', 2);

grid on;
xlabel('SNR (dB)');
ylabel('BER');
legend(...
    'OFDM ZF', 'MC-OTFS ZF', 'Zak-OTFS ZF', ...
    'OFDM MMSE', 'MC-OTFS MMSE', 'Zak-OTFS MMSE', ...
    'OFDM MRC', 'MC-OTFS MRC', 'Zak-OTFS MRC', ...
    'Location', 'southwest');
title('Waveform Performance Comparison (ZF, MMSE, & MRC)');
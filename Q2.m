% Q2: SIMO (1 x nr) capacity: one fixed Rayleigh realization vs AWGN SIMO
clear; clc; close all;
rng(1);                                    % reproducible "fixed" channel

% Parameters
nr      = 4;                               % number of receive antennas
SNR_dB  = 0:5:40;                          % SNR in dB
SNR_lin = 10.^(SNR_dB/10);                 % SNR in linear scale

% One fixed channel realization h ~ CN(0, I), size nr x 1
h = (randn(nr,1) + 1i*randn(nr,1))/sqrt(2);
h_norm2 = norm(h)^2; % ||h||^2 of THIS realization
fprintf('||h||^2 of the fixed realization = %.3f (E[||h||^2] = nr = %d)\n', h_norm2, nr);

% (a) Fixed (single-realization) SIMO capacity, receiver combines optimally (MRC)
%     C = log2(1 + ||h||^2 * SNR). It is the INSTANTANEOUS capacity of one draw.
C_FIXED_SIMO = log2(1 + h_norm2 * SNR_lin);

% (b) AWGN SIMO: deterministic channel with unit gain on every branch (h = ones),
%     so ||h||^2 = nr exactly:  C = log2(1 + nr*SNR)  (array gain of 10*log10(nr) dB)
C_AWGN_SIMO = log2(1 + nr * SNR_lin);

% Display results
disp(table(SNR_dB.', SNR_lin.', C_FIXED_SIMO.', C_AWGN_SIMO.', 'VariableNames', ...
    {'SNR_dB','SNR_linear','Capacity_FIXED_SIMO','Capacity_AWGN_SIMO'}));

% Plot
figure(2);
plot(SNR_dB, C_FIXED_SIMO, '-o', 'LineWidth', 1.5, 'MarkerSize', 5);
hold on;
plot(SNR_dB, C_AWGN_SIMO, '-s', 'LineWidth', 1.5, 'MarkerSize', 5, 'Color', [1 0 0]);
hold off;
grid on;
xlabel('SNR (dB)');
ylabel('Capacity (bits/s/Hz)');
title('SIMO Channel Capacity: Fixed Rayleigh Draw vs AWGN');
legend('Fixed SIMO channel (one Rayleigh realization)', 'AWGN SIMO (unit-gain branches)', 'Location', 'best');

% Q1: Capacity of the AWGN SISO channel
clear; clc; close all;

SNR_dB  = 0:5:40;              % SNR in dB
SNR_lin = 10.^(SNR_dB/10);     % SNR in linear scale

% Capacity (bits/s/Hz): C = log2(1 + SNR)
C = log2(1 + SNR_lin);

% Display results
disp(table(SNR_dB.', SNR_lin.', C.', 'VariableNames', ...
    {'SNR_dB','SNR_linear','Capacity_bits_per_Hz'}));

% Plot
figure(1);
plot(SNR_lin, C, '-o', 'LineWidth', 1.5, 'MarkerSize', 5, 'Color', [1 0 0]);
grid on;
xlabel('SNR (dB)');
ylabel('Capacity (bits/s/Hz)');
title('AWGN SISO Channel Capacity');

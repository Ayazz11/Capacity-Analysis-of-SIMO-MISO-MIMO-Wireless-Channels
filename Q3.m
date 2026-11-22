% Q3: Ergodic capacity of the fast-fading (Rayleigh) SIMO channel
clear; clc; close all;

% Parameters
SNR_dB  = 0:5:40;                % SNR in dB
SNR_lin = 10.^(SNR_dB/10);       % SNR in linear scale
nr      = 4;                     % number of receive antennas
n       = 1000;                  % number of channel realizations

% Generate n Rayleigh channel vectors h ~ CN(0, I_nr); column t is realization t
Hs = (randn(nr,n) + 1i*randn(nr,n))/sqrt(2);
G  = sum(abs(Hs).^2, 1);         % ||h||^2 for each realization (1 x n)
Eh=mean(G);
disp(Eh);
C_ergodic_SIMO = zeros(size(SNR_lin));
for k = 1:length(SNR_lin)
    % instantaneous capacity of every realization at this SNR, then average
    C_inst = log2(1 + G * SNR_lin(k));
    C_ergodic_SIMO(k) = mean(C_inst);          % ergodic capacity
end

% AWGN SIMO upper reference (Jensen: E[log2(1+||h||^2 SNR)] <= log2(1+E||h||^2 SNR))
C_AWGN_SIMO = log2(1 + nr * SNR_lin);

% Display results
disp(table(SNR_dB.', SNR_lin.', C_ergodic_SIMO.', C_AWGN_SIMO.', 'VariableNames', ...
    {'SNR_dB','SNR_linear','Ergodic_SIMO_Rayleigh','AWGN_SIMO'}));

% Plot
figure(3);
plot(SNR_dB, C_ergodic_SIMO, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
plot(SNR_dB, C_AWGN_SIMO, '--s', 'LineWidth', 1.5, 'MarkerSize', 5, 'Color', [1 0 0]);
hold off;
grid on;
xlabel('SNR (dB)');
ylabel('Ergodic Capacity (bits/s/Hz)');
title('Fast-Fading (Rayleigh) SIMO Channel Capacity');
legend('Ergodic SIMO (Rayleigh)', 'AWGN SIMO bound', 'Location', 'best');

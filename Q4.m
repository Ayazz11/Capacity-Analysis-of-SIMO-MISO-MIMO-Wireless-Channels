% Q4: Ergodic capacity of a 4x4 Rayleigh MIMO channel vs 1x4 SIMO
clear; clc; close all;
rng(1);

% Parameters
SNR_dB  = 0:5:40;                % SNR in dB
SNR_lin = 10.^(SNR_dB/10);       % SNR in linear scale
nr      = 4;                     % number of receive antennas
nt      = 4;                     % number of TRANSMIT antennas
n       = 1000;                  % number of channel realizations

% Generate channels once and store eigenvalues / SIMO gains
L = zeros(n, min(nr,nt));        % eigenvalues of H*H' for each realization
G = zeros(n, 1);                 % ||h||^2 of the 1x4 SIMO link (first Tx antenna)
for t = 1:n
    H = (randn(nr,nt) + 1i*randn(nr,nt))/sqrt(2);   % H ~ CN(0,1) entries
    L(t,:) = sort(real(eig(H*H')), 'descend').';    % real(): removes numerical imag. residue
    G(t)   = norm(H(:,1))^2;                        % SIMO channel = first column of H
end

C_ergodic_MIMO = zeros(size(SNR_lin));
C_ergodic_SIMO = zeros(size(SNR_lin));
for k = 1:length(SNR_lin)
    snr = SNR_lin(k);
    % MIMO: C = sum_i log2(1 + (SNR/nt) * lambda_i), total Tx power split equally over nt
    C_ergodic_MIMO(k) = mean(sum(log2(1 + (snr/nt) * L), 2));
    % SIMO: C = log2(1 + ||h||^2 * SNR)
    C_ergodic_SIMO(k) = mean(log2(1 + G * snr));
end

% Display results
disp(table(SNR_dB.', SNR_lin.', C_ergodic_MIMO.', C_ergodic_SIMO.', 'VariableNames', ...
    {'SNR_dB','SNR_linear','Ergodic_MIMO_4x4','Ergodic_SIMO_1x4'}));

% Sanity check: average ordered eigenvalues of H*H'; their sum should be ~ nt*nr = 16
disp('Mean ordered eigenvalues of H*H'' (largest to smallest):');
disp(mean(L, 1));
fprintf('Sum of mean eigenvalues = %.2f (theory: nt*nr = %d)\n', sum(mean(L,1)), nt*nr);

% Plot
figure(4);
plot(SNR_dB, C_ergodic_MIMO, '-o', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [1 0 0]);
hold on;
plot(SNR_dB, C_ergodic_SIMO, '-s', 'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0 0 1]);
hold off;
grid on;
xlabel('SNR (dB)');
ylabel('Ergodic Capacity (bits/s/Hz)');
title('Fast-Fading 4x4 MIMO vs 1x4 SIMO Capacity');
legend('Ergodic MIMO 4x4', 'Ergodic SIMO 1x4', 'Location', 'best');

clear; clc; close all;
rng(1);

% Parameters
SNR_dB  = 0:5:40;
SNR_lin = 10.^(SNR_dB/10);
nr = 4;  nt = 4;
n  = 5000;                               % Rayleigh realizations

G1 = zeros(n,1);                         % |h_SISO|^2   (SISO: H(1,1))
G4 = zeros(n,1);                         % ||h_SIMO||^2 (SIMO 1x4: H(:,1)
L  = zeros(n, min(nr,nt));               % eigenvalues of H*H' 
for t = 1:n
    H = (randn(nr,nt) + 1i*randn(nr,nt))/sqrt(2);
    G1(t)  = abs(H(1,1))^2;
    G4(t)  = norm(H(:,1))^2;
    L(t,:) = real(eig(H*H')).';
end

C_SISO_AWGN = log2(1 + SNR_lin);
C_SIMO_AWGN = log2(1 + nr*SNR_lin);                    % h = ones(nr,1)
C_MIMO_AWGN = min(nr,nt) * log2(1 + SNR_lin/nt);       % H = I: 4 parallel unit-gain channels

C_SISO_RAY = zeros(size(SNR_lin));
C_SIMO_RAY = zeros(size(SNR_lin));
C_MIMO_RAY = zeros(size(SNR_lin));          % equal power, no CSIT
C_MIMO_WF  = zeros(size(SNR_lin));          % water-filling, CSIT available
for k = 1:length(SNR_lin)
    snr = SNR_lin(k);
    C_SISO_RAY(k) = mean(log2(1 + G1*snr));
    C_SIMO_RAY(k) = mean(log2(1 + G4*snr));
    C_MIMO_RAY(k) = mean(sum(log2(1 + (snr/nt)*L), 2));
    cwf = zeros(n,1);
    for t = 1:n
        cwf(t) = waterfill_capacity(L(t,:), snr);
    end
    C_MIMO_WF(k) = mean(cwf);
end

disp(table(SNR_dB.', C_SISO_AWGN.', C_SISO_RAY.', C_SIMO_AWGN.', C_SIMO_RAY.', ...
    C_MIMO_AWGN.', C_MIMO_RAY.', C_MIMO_WF.', 'VariableNames', ...
    {'SNR_dB','SISO_AWGN','SISO_Rayleigh','SIMO_AWGN','SIMO_Rayleigh', ...
     'MIMO_AWGN','MIMO_Rayleigh','MIMO_Rayleigh_WF'}));

figure(5);
plot(SNR_dB, C_SISO_AWGN, 'k-o',  'LineWidth', 1.5, 'MarkerSize', 5); hold on;
plot(SNR_dB, C_SISO_RAY,  'k--o', 'LineWidth', 1.5, 'MarkerSize', 5);
plot(SNR_dB, C_SIMO_AWGN, 'b-s',  'LineWidth', 1.5, 'MarkerSize', 5);
plot(SNR_dB, C_SIMO_RAY,  'b--s', 'LineWidth', 1.5, 'MarkerSize', 5);
plot(SNR_dB, C_MIMO_AWGN, 'r-^',  'LineWidth', 1.5, 'MarkerSize', 5);
plot(SNR_dB, C_MIMO_RAY,  'r--^', 'LineWidth', 1.5, 'MarkerSize', 5);
plot(SNR_dB, C_MIMO_WF,   'm:d',  'LineWidth', 1.5, 'MarkerSize', 5);
hold off; grid on;
xlabel('SNR (dB)'); ylabel('Capacity (bits/s/Hz)');
title('Capacity of SISO / SIMO (1x4) / MIMO (4x4): AWGN vs Rayleigh');
legend('SISO AWGN','SISO Rayleigh (ergodic)','SIMO 1x4 AWGN','SIMO 1x4 Rayleigh (ergodic)', ...
       'MIMO 4x4 AWGN (H = I)','MIMO 4x4 Rayleigh (equal power)','MIMO 4x4 Rayleigh (water-filling)', ...
       'Location', 'northwest');

i30 = find(SNR_dB == 30);  i40 = find(SNR_dB == 40);   % 10 dB apart
fprintf('\nHigh-SNR slope (bits/s/Hz per 10 dB, from 30 to 40 dB):\n');
fprintf('  SISO Rayleigh : %.2f  (theory 3.32)\n',  C_SISO_RAY(i40)-C_SISO_RAY(i30));
fprintf('  SIMO Rayleigh : %.2f  (theory 3.32, same slope, shifted)\n', C_SIMO_RAY(i40)-C_SIMO_RAY(i30));
fprintf('  MIMO Rayleigh : %.2f  (theory %.2f = min(nt,nr)*3.32)\n', ...
        C_MIMO_RAY(i40)-C_MIMO_RAY(i30), min(nt,nr)*log2(10));
fprintf('Array gain of SIMO over SISO (AWGN): %.2f dB\n', 10*log10(nr));
fprintf('Rayleigh loss vs AWGN at 20 dB (SISO): %.2f bits/s/Hz (Jensen gap)\n', ...
        C_SISO_AWGN(SNR_dB==20) - C_SISO_RAY(SNR_dB==20));

% ---------- Figure 6: outage CDF at a fixed SNR ----------
snr_o  = 10^(10/10);                          % 10 dB
Cs = sort(log2(1 + G1*snr_o));
Ca = sort(log2(1 + G4*snr_o));
Cm = sort(sum(log2(1 + (snr_o/nt)*L), 2));
cdf = (1:n)'/n;
figure(6);
plot(Cs, cdf, 'k', 'LineWidth', 1.5); hold on;
plot(Ca, cdf, 'b', 'LineWidth', 1.5);
plot(Cm, cdf, 'r', 'LineWidth', 1.5);
hold off; grid on;
xlabel('Rate R (bits/s/Hz)'); ylabel('Outage probability  P(C < R)');
title('Outage CDF of instantaneous capacity at 10 dB SNR');
legend('SISO','SIMO 1x4','MIMO 4x4', 'Location', 'southeast');
idx = ceil(0.1*n);
fprintf('\n10%% outage capacity at 10 dB: SISO %.2f | SIMO %.2f | MIMO %.2f bits/s/Hz\n', ...
        Cs(idx), Ca(idx), Cm(idx));

% ---------- local function ----------
function C = waterfill_capacity(lambda, snr)
% Capacity with water-filling over eigenmodes lambda (noise power 1, total power snr)
    lambda = sort(lambda(:), 'descend');
    lambda = lambda(lambda > 1e-12);
    for r = numel(lambda):-1:1
        mu = (snr + sum(1 ./ lambda(1:r))) / r;      % water level for r active modes
        if mu > 1/lambda(r), break; end
    end
    p = max(mu - 1 ./ lambda(1:r), 0);
    C = sum(log2(1 + lambda(1:r) .* p));
end
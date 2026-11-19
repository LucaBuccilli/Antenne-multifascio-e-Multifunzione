%% p06 - Multifascio: fasci ortogonali (matrice di Butler = FFT spaziale)
% Array di N = 8 elementi a lambda/2. Fasci in sin(theta_m) = (m - 1/2)/(N d).
% Si verificano posizione, ortogonalita', incrocio a -3,92 dB ed equivalenza con la FFT.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

N = 8; d = 0.5;
u = linspace(-1, 1, 40001);
ms = -N/2+1 : N/2;
u_beams = (ms - 0.5) / (N*d);
n = 0:N-1;
fprintf('p06 - Multifascio di Butler (N = 8, d = lambda/2)\n');
fig = figure('Position', [100 100 900 450]); hold on;
for k = 1:N
    w = exp(-1j*2*pi*d*u_beams(k)*n) / N;
    p = abs(af(w, u, d));
    [~, im] = max(p);
    riga(sprintf('fascio %+d: direzione', ms(k)), asind(u(im)), asind(u_beams(k)), 'deg');
    plot(asind(u), 20*log10(p + 1e-9), 'LineWidth', 1.2, 'DisplayName', sprintf('fascio %d', ms(k)));
end
G = zeros(N);
for k = 1:N
    w = exp(-1j*2*pi*d*u_beams(k)*n) / N;
    G(k, :) = abs(af(w, u_beams, d)).';
end
riga('massimo fuori diagonale (ortogonalita'')', 20*log10(max(max(G - diag(diag(G)))) + 1e-15), -300, 'dB');
uc = (u_beams(4) + u_beams(5)) / 2;
w = exp(-1j*2*pi*d*u_beams(4)*n) / N;
riga('incrocio tra fasci adiacenti (teoria N grande)', 20*log10(abs(af(w, uc, d))), -3.92, 'dB');
x = exp(1j*2*pi*d*u_beams(6)*n);
y = abs(fftshift(fft(x .* exp(1j*pi*n/N)))) / N;
fprintf('  FFT: uscite dei canali per un''onda dal fascio %+d: %s\n', ms(6), mat2str(round(y, 3)));
yline(-3.92, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
xlabel('\theta [gradi]'); ylabel('guadagno [dB]'); ylim([-30 1]); xlim([-90 90]); grid on;
legend('NumColumns', 4, 'Location', 'south');
title('8 fasci simultanei ortogonali: ognuno ha il massimo sugli zeri degli altri');
fprintf(' figura: %s\n', salva_fig(fig, 'p06_multifascio'));

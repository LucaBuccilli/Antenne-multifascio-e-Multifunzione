%% p07 - Monopulse con array: fasci somma e differenza, stima dell'angolo
% N = 32 elementi a lambda/2. Sigma = S + D, Delta = S - D (meta' sinistra e destra).
% Im{-Delta/Sigma} = tan(N pi (u - u0)/4), pendenza km = N pi / 4.
% Monte Carlo: sigma_u = 1/(km sqrt(2 SNR_Sigma)).
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

N = 32; d = 0.5; n = 0:N-1;
km = N*pi/4;
sx = double(n < N/2); dx = 1 - sx;
rng(7);
canale_S = @(X) sum(X, 2);
canale_D = @(X) X*sx.' - X*dx.';

uu = linspace(-0.04, 0.04, 801).';
X = exp(1j*2*pi*d*uu*n);
r = imag(-canale_D(X) ./ canale_S(X));
fprintf('p07 - Monopulse (N = 32, d = lambda/2)\n');
sel = abs(uu) < 0.003;
cf = polyfit(uu(sel), r(sel), 1);
riga('pendenza monopulse km', cf(1), km);
riga('rapporto in u = 0,03 vs tan(N pi u / 4)', interp1(uu, r, 0.03), tan(N*pi*0.03/4));

u_t = 0.003; prove = 4000;
fprintf(' Monte Carlo: bersaglio in u = %.3f, %d prove per ogni SNR\n', u_t, prove);
snr_list = 10:5:30; sig_mis = zeros(size(snr_list)); sig_teo = sig_mis;
s = exp(1j*2*pi*d*u_t*n);
for k = 1:numel(snr_list)
    snr_sigma = 10^(snr_list(k)/10);
    sigma_n = sqrt(N/snr_sigma);
    rum = sigma_n/sqrt(2) * (randn(prove, N) + 1j*randn(prove, N));
    Xs = s + rum;
    u_hat = imag(-canale_D(Xs) ./ canale_S(Xs)) / km;
    sig_mis(k) = std(u_hat);
    sig_teo(k) = 1/(km*sqrt(2*snr_sigma));
    riga(sprintf('sigma_u con SNR_Sigma = %d dB (x 1e-3)', snr_list(k)), 1e3*sig_mis(k), 1e3*sig_teo(k));
end
bw_u = 0.886/(N*d);
fprintf('  con SNR = 20 dB la precisione e'' 1/%.0f della larghezza del fascio\n', bw_u/sig_teo(3));

fig = figure('Position', [100 100 1100 420]);
subplot(1, 2, 1); plot(uu, r, 'LineWidth', 1.6); hold on; plot(uu, km*uu, '--', 'Color', [0.5 0.5 0.5]);
ylim([-3 3]); grid on; xlabel('u - u_0'); ylabel('rapporto monopulse'); legend('Im\{-\Delta/\Sigma\}', 'k_m (u - u_0)');
title('Rapporto monopulse: lineare vicino al centro');
subplot(1, 2, 2); semilogy(snr_list, sig_mis, 'o', snr_list, sig_teo, '-'); grid on;
xlabel('SNR canale somma [dB]'); ylabel('\sigma della stima di u'); legend('Monte Carlo', 'teoria'); title('Precisione angolare');
fprintf(' figura: %s\n', salva_fig(fig, 'p07_monopulse'));

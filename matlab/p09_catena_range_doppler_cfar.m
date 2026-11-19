%% p09 - Catena di elaborazione di un fascio: compressione, MTI, Doppler, CA-CFAR
% Banda S (3 GHz), M = 32 impulsi chirp (B = 5 MHz, T = 10 us), due bersagli e clutter forte.
% Catena: filtro adattato -> MTI a 2 impulsi -> FFT Doppler (Chebyshev 60 dB) -> |.|^2 -> CA-CFAR.
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

c = 3e8; fc = 3e9; lam = c/fc;
B = 5e6; T = 10e-6; fs = 10e6;
PRF = 4000; M = 32;
Nr = round(fs/PRF);
rng(11);
t = (0:round(T*fs)-1)/fs;
chirp_s = exp(1j*pi*(B/T)*t.^2);
L = numel(chirp_s);
bersagli = [12000 40 0.3; 20010 -25 0.2];          % distanza [m], velocita' [m/s], ampiezza

fprintf('p09 - Catena range-Doppler e CA-CFAR (3 GHz, B = 5 MHz, M = 32)\n');

% 1) compressione con chirp sovracampionato (x20)
fs_os = 20*B;
t_os = (0:round(T*fs_os)-1)/fs_os;
ch_os = exp(1j*pi*(B/T)*t_os.^2);
prof = abs(conv(ch_os, conj(fliplr(ch_os)))).^2;
prof_db = 10*log10(prof/max(prof));
[~, i0] = max(prof_db);
x_m = ((1:numel(prof_db)) - i0)/fs_os*c/2;
riga('risoluzione in distanza (larghezza -3 dB)', larghezza_3db(x_m, prof_db), 0.886*c/(2*B), 'm');
k1 = round(1.5*fs_os/B);
lob = prof_db; lob(i0-k1:i0+k1) = -200;
riga('lobo laterale di compressione', max(lob), -13.26, 'dB');
g = max(abs(conv(chirp_s, conj(fliplr(chirp_s)))))^2 / L;
riga('guadagno di compressione (campioni T*fs)', 10*log10(g), 10*log10(T*fs), 'dB');

% 2) segnale ricevuto ed elaborazione
rx = genera_eco(true, true, bersagli, M, Nr, PRF, lam, c, fs, chirp_s);
[mf, P] = elabora(rx, chirp_s, M, Nr);
[riv, alpha] = ca_cfar(P, 32, 4, 1e-6);
dv = lam*PRF/(2*(M-1));
vel = ((0:M-2) - floor((M-1)/2)) * PRF/(M-1) * lam/2;
fprintf('  soglia CA-CFAR alpha = %.2f (%.1f dB), risoluzione in velocita'' %.2f m/s\n', alpha, 10*log10(alpha), dv);
for b = 1:size(bersagli, 1)
    R = bersagli(b, 1); v = bersagli(b, 2);
    k = round(2*R/c*fs) + 1;                            % indice MATLAB (base 1)
    sub = P(:, k-5:k+5) .* riv(:, k-5:k+5);
    if max(sub(:)) == 0
        fprintf('  bersaglio a %.2f km NON rivelato\n', R/1e3); continue;
    end
    [~, idx] = max(sub(:)); [j, i] = ind2sub(size(sub), idx);
    riga(sprintf('bersaglio %.2f km: distanza', R/1e3), (k-6+i-1)/fs*c/2, R, 'm');
    riga(sprintf('bersaglio %.2f km: velocita'' (entro %.1f m/s)', R/1e3, dv), vel(j), v, 'm/s');
end
fprintf('  celle oltre soglia in totale: %d\n', nnz(riv));

% 3) falsi allarmi su solo rumore, Pfa nominale 1e-3
fa = 0; celle = 0; fa2 = 0; celle2 = 0;
for prova = 1:5
    [~, Pn] = elabora(genera_eco(false, false, bersagli, M, Nr, PRF, lam, c, fs, chirp_s), chirp_s, M, Nr);
    dn = ca_cfar(Pn, 32, 4, 1e-3);   fa = fa + nnz(dn);   celle = celle + numel(dn(:, 21:end-20));
    dn2 = ca_cfar(Pn(:, 1:2:end), 32, 4, 1e-3); fa2 = fa2 + nnz(dn2); celle2 = celle2 + numel(dn2(:, 21:end-20));
end
riga('Pfa con tutte le celle (x 1e-3)', 1e3*fa/celle, 1);
riga('Pfa con celle indipendenti (x 1e-3)', 1e3*fa2/celle2, 1);

fig = figure('Position', [100 100 1200 440]);
rng_km = (0:Nr-1)/fs*c/2/1e3;
subplot(1, 2, 1); plot(rng_km, 20*log10(abs(mf(1, :)) + 1e-9), 'Color', [0.4 0.4 0.4]); grid on;
xlabel('distanza [km]'); ylabel('ampiezza [dB]'); title('Dopo la compressione: il clutter copre i bersagli');
subplot(1, 2, 2);
imagesc(rng_km, vel, 10*log10(P/median(P(:))), [0 40]); axis xy; colorbar; hold on;
[dj, di] = find(riv);
plot(rng_km(di), vel(dj), 'r.', 'MarkerSize', 8);
xlabel('distanza [km]'); ylabel('velocità radiale [m/s]'); title('Mappa range-Doppler dopo MTI e FFT');
fprintf(' figura: %s\n', salva_fig(fig, 'p09_range_doppler'));

%% funzioni locali
function rx = genera_eco(con_bersagli, con_clutter, bersagli, M, Nr, PRF, lam, c, fs, chirp_s)
refl = zeros(M, Nr);
m = (0:M-1).';
if con_bersagli
    for b = 1:size(bersagli, 1)
        k = round(2*bersagli(b, 1)/c*fs) + 1;
        fd = 2*bersagli(b, 2)/lam;
        refl(:, k) = refl(:, k) + bersagli(b, 3) * exp(1j*2*pi*fd*m/PRF) * exp(1j*2*pi*rand);
    end
end
if con_clutter
    celle = 151:(Nr-200);
    amp = sqrt(10^(30/10)/2) * (randn(1, numel(celle)) + 1j*randn(1, numel(celle)));
    fd_c = 4.0 * randn(1, numel(celle));
    refl(:, celle) = refl(:, celle) + amp .* exp(1j*2*pi*m*fd_c/PRF);
end
rx = zeros(M, Nr);
for k = 1:M
    y = conv(refl(k, :), chirp_s);
    rx(k, :) = y(1:Nr);
end
rx = rx + (randn(M, Nr) + 1j*randn(M, Nr))/sqrt(2);
end

function [mf, P] = elabora(rx, chirp_s, M, Nr)
L = numel(chirp_s);
mf = zeros(M, Nr);
for k = 1:M
    y = conv(rx(k, :), conj(fliplr(chirp_s)));
    mf(k, :) = y(L:L-1+Nr) / sqrt(L);
end
mti = mf(2:end, :) - mf(1:end-1, :);
win = cheb_win(M-1, 60);
rd = fftshift(fft(mti .* win, [], 1), 1);
P = abs(rd).^2;
end

function [riv, alpha] = ca_cfar(P, n_ref, n_guard, pfa)
alpha = n_ref * (pfa^(-1/n_ref) - 1);
mezzo = [ones(1, n_ref/2), zeros(1, n_guard)];
kern = [mezzo, 0, fliplr(mezzo)];
rif = zeros(size(P));
for k = 1:size(P, 1)
    rif(k, :) = conv(P(k, :), kern, 'same') / n_ref;
end
riv = P > alpha*rif;
bordo = n_ref/2 + n_guard;
riv(:, 1:bordo) = false;
riv(:, end-bordo+1:end) = false;
end

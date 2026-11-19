%% p04 - Sfasatori digitali a B bit: perdita di guadagno e lobi di quantizzazione
% Array di N = 64 elementi a lambda/2, theta0 = 10 deg, ampiezze Taylor -40 dB.
%   - perdita di guadagno media: teoria (pi^2/3) 2^(-2B)
%   - lobo di quantizzazione con errore periodico: circa -6B dB
%   - con una fase casuale nota aggiunta e poi sottratta il lobo spurio sparisce
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

N = 64; d = 0.5; theta0 = 10;
u0 = sind(theta0);
u = linspace(-1, 1, 80001);
n = 0:N-1;
fase_ideale = -2*pi*d*u0*n;
amp = taylor_win(N, 6, 40).';
rng(1);
dither = 2*pi*rand(1, N);                        % fase fissa nota di ogni elemento
quantizza = @(phi, B) round(phi / (2*pi/2^B)) * (2*pi/2^B);

ideale = abs(af(amp .* exp(1j*fase_ideale), u, d));
G0 = max(ideale);
escl = 8 / (N*d);
fuori = abs(u - u0) > escl;

fprintf('p04 - Quantizzazione degli sfasatori (N = 64, theta0 = 10 deg)\n');
fig = figure('Position', [100 100 900 450]); hold on;
plot(asind(u), 20*log10(ideale/G0 + 1e-12), 'Color', [0.7 0.7 0.7], 'DisplayName', 'fasi ideali (Taylor -40 dB)');
for B = 3:5
    w_q = amp .* exp(1j*quantizza(fase_ideale, B));
    w_d = amp .* exp(1j*(quantizza(fase_ideale + dither, B) - dither));
    pq = abs(af(w_q, u, d)); pd = abs(af(w_d, u, d));
    fprintf(' B = %d bit (passo %.2f deg)\n', B, 360/2^B);
    riga('perdita di guadagno con dithering', -20*log10(max(pd)/G0), -10*log10(1 - pi^2/3*2^(-2*B)), 'dB');
    dbq = 20*log10(pq/max(pq) + 1e-12); dbd = 20*log10(pd/max(pd) + 1e-12);
    riga('lobo di quantizzazione (errore periodico)', max(dbq(fuori)), -6.02*B, 'dB');
    fprintf('  %-42s misurato = %10.3f dB\n', 'lobo spurio piu alto con dithering', max(dbd(fuori)));
    if B == 3
        plot(asind(u), dbq, 'Color', [0.85 0.2 0.2], 'DisplayName', '3 bit, errore periodico');
        plot(asind(u), dbd, 'Color', [0 0.45 0.74], 'DisplayName', '3 bit, con fase casuale');
    end
end
fprintf('  nota: il livello esatto del lobo di quantizzazione dipende da theta0; -6B dB e'' un valore indicativo.\n');
xlabel('\theta [gradi]'); ylabel('|AF| [dB]'); ylim([-60 2]); xlim([-90 90]); grid on; legend;
title('Sfasatori a 3 bit: l''errore periodico crea lobi spuri, la fase casuale li disperde');
fprintf(' figura: %s\n', salva_fig(fig, 'p04_quantizzazione'));

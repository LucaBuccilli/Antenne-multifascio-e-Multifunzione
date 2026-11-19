%% p02 - Rastremazione delle ampiezze (tapering)
% Array di N = 32 elementi a lambda/2: uniforme, Hamming, Taylor (-30 dB, nbar = 4)
% e Dolph-Chebyshev (-30 dB). Si misurano lobo laterale massimo, allargamento del
% fascio rispetto all'uniforme ed efficienza |sum w|^2 / (N sum w^2).
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

N = 32; d = 0.5;
u = linspace(-1, 1, 40001);
nomi = {'Uniforme', 'Hamming', 'Taylor -30 dB, nbar=4', 'Chebyshev -30 dB'};
pesi = {ones(N,1), hamming_win(N), taylor_win(N, 4, 30), cheb_win(N, 30)};
sll_teo = [-13.26, -42.7, -30, -30];

fprintf('p02 - Rastremazione (N = 32, d = lambda/2)\n');
fig = figure('Position', [100 100 900 450]); hold on;
bw_uni = [];
for k = 1:4
    w = pesi{k};
    p = db_norm(af(w, u, d));
    bw = larghezza_3db(u, p);
    if isempty(bw_uni), bw_uni = bw; end
    eff = abs(sum(w))^2 / (N * sum(abs(w).^2));
    fprintf(' %s\n', nomi{k});
    riga('lobo laterale massimo', lobo_laterale_max(p), sll_teo(k), 'dB');
    fprintf('  %-42s misurato = %10.3f\n', 'allargamento del fascio (rispetto uniforme)', bw/bw_uni);
    fprintf('  %-42s misurato = %10.3f       perdita = %.2f dB\n', 'efficienza di rastremazione', eff, -10*log10(eff));
    plot(asind(u), p, 'LineWidth', 1.3, 'DisplayName', nomi{k});
end
xlabel('angolo \theta [gradi]'); ylabel('|AF| [dB]'); ylim([-60 2]); xlim([-90 90]);
title('Lobi laterali più bassi = fascio più largo e meno guadagno'); grid on; legend;
fprintf(' figura: %s\n', salva_fig(fig, 'p02_rastremazione'));

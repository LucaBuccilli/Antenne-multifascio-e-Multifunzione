%% p03 - Grating lobes e spaziatura degli elementi
% Array di 16 elementi puntato a theta0 = 50 deg. Teoria: massimi in
% u_m = u0 + m/d; il grating lobe entra nello spazio visibile se d > 1/(1 + sin theta0).
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

N = 16; theta0 = 50;
u0 = sind(theta0);
u = linspace(-1, 1, 40001);
fprintf('p03 - Grating lobes (N = 16, theta0 = 50 deg)\n');
fprintf('  spaziatura massima senza grating lobes: d < %.3f lambda\n', 1/(1+u0));
fig = figure('Position', [100 100 1000 600]);
spaziature = [0.5 0.566 0.7 1.0];
for k = 1:4
    d = spaziature(k);
    w = exp(-1j*2*pi*d*u0*(0:N-1));
    p = db_norm(af(w, u, d));
    ug = u0 - 1/d;                                   % grating lobe m = -1
    visibile = abs(ug) <= 1;
    if visibile, s = 'VISIBILE'; else, s = 'non visibile'; end
    fprintf(' d = %.3f lambda: grating lobe teorico in u = %+.3f -> %s\n', d, ug, s);
    if visibile
        mask = abs(u - ug) < 0.1;
        uu = u(mask); pp = p(mask);
        [pmax, im] = max(pp);
        riga('posizione grating lobe (gradi)', asind(uu(im)), asind(ug), 'deg');
        riga('livello grating lobe', pmax, 0, 'dB');
    end
    subplot(2, 2, k);
    plot(asind(u), p, 'LineWidth', 1.2); ylim([-40 2]); grid on;
    title(sprintf('d = %.3g \\lambda', d)); xlabel('\theta [gradi]');
end
sgtitle('Puntamento a 50°: oltre 0,566 \lambda compare un secondo lobo alto quanto il principale');
fprintf(' figura: %s\n', salva_fig(fig, 'p03_grating_lobes'));

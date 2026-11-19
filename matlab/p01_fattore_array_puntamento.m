%% p01 - Fattore di array e puntamento elettronico
% Array lineare uniforme di N dipoli a d = lambda/2. Lo sfasamento progressivo
% alpha = -k d sin(theta0) punta il fascio in theta0. Si verifica:
%   - direzione del massimo = theta0
%   - larghezza a -3 dB = 0.886 lambda / (N d cos theta0)
%   - primo lobo laterale = -13.26 dB (N grande)
%   - perdita di scansione con diagramma dell'elemento cos(theta)
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

N = 16; d = 0.5;
theta = linspace(-90, 90, 36001);
u = sind(theta);

fprintf('p01 - Fattore di array e puntamento (N = 16, d = lambda/2)\n');
fig = figure('Position', [100 100 900 450]); hold on;
colori = {[0.6 0.6 0.6], [0 0.45 0.74], [0.85 0.2 0.2]};
angoli = [0 30 60];
for k = 1:3
    theta0 = angoli(k);
    alpha = -2*pi*d*sind(theta0);                 % sfasamento progressivo
    w = exp(1j*alpha*(0:N-1));                     % pesi = solo fase
    p = db_norm(af(w, u, d));
    fprintf(' theta0 = %d deg  (alpha = %.1f deg tra elementi)\n', theta0, rad2deg(alpha));
    [~, im] = max(p);
    riga('direzione del massimo', theta(im), theta0, 'deg');
    riga('larghezza a -3 dB', larghezza_3db(theta, p), rad2deg(0.886/(N*d*cosd(theta0))), 'deg');
    riga('lobo laterale massimo (teoria per N grande)', lobo_laterale_max(p), -13.26, 'dB');
    plot(theta, p, 'Color', colori{k}, 'LineWidth', 1.6, 'DisplayName', sprintf('puntato a %d°', theta0));
end

fprintf(' Perdita di scansione (diagramma dell''elemento cos(theta))\n');
p_ref = [];
for theta0 = [0 30 45 60]
    w = exp(-1j*2*pi*d*sind(theta0)*(0:N-1));
    tot = cosd(theta(:)) .* abs(af(w, u, d));
    if isempty(p_ref), p_ref = max(tot); end
    riga(sprintf('perdita a %d deg (sola andata)', theta0), 20*log10(max(tot)/p_ref), 20*log10(cosd(theta0)), 'dB');
end

yline(-13.26, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility', 'off');
xlabel('angolo \theta dalla normale [gradi]'); ylabel('|AF| normalizzato [dB]'); ylim([-50 2]);
title('Puntamento elettronico: il fascio si sposta e si allarga'); grid on; legend('Location', 'southwest');
fprintf(' figura: %s\n', salva_fig(fig, 'p01_puntamento'));

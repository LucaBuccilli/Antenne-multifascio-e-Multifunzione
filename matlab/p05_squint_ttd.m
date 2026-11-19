%% p05 - Squint di frequenza: sfasatori contro ritardi veri (TTD)
% Array di N = 64 elementi a d = lambda0/2 puntato a 45 deg a f0; f tra 0,9 e 1,1 f0.
%   sfasatori: sin(theta) = (f0/f) sin(theta0)      TTD: theta = theta0 sempre
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

N = 64; theta0 = 45; d0 = 0.5;
n = 0:N-1;
u = linspace(-1, 1, 200001);
u0 = sind(theta0);
rapporti = linspace(0.9, 1.1, 21);
patt = @(fr, ttd) abs(exp(1j*2*pi*d0*fr*u(:)*n) * exp(1j*(-2*pi*(d0*(ttd*fr + (1-ttd)))*u0*n)).');

fprintf('p05 - Squint (N = 64, theta0 = 45 deg)\n');
dir_ps = zeros(size(rapporti)); dir_ttd = dir_ps;
for k = 1:numel(rapporti)
    [~, i1] = max(patt(rapporti(k), 0)); dir_ps(k) = asind(u(i1));
    [~, i2] = max(patt(rapporti(k), 1)); dir_ttd(k) = asind(u(i2));
end
for fr = [0.95 1.05]
    k = find(abs(rapporti - fr) < 1e-9);
    riga(sprintf('direzione con sfasatori a f = %.2f f0', fr), dir_ps(k), asind(u0/fr), 'deg');
    riga(sprintf('direzione con TTD a f = %.2f f0', fr), dir_ttd(k), theta0, 'deg');
end
k = find(abs(rapporti - 1.02) < 1e-9);
riga('squint per df/f0 = +2% (linearizzato)', dir_ps(k) - theta0, -rad2deg(0.02*tand(theta0)), 'deg');
[~, i0] = min(abs(u - u0));
p1 = patt(1.02, 0); p2 = patt(1.02, 1);
fprintf('  perdita in theta0 a f = 1,02 f0 con sfasatori: %.2f dB\n', 20*log10(p1(i0)/max(p1)));
fprintf('  perdita in theta0 a f = 1,02 f0 con TTD:       %.2f dB\n', 20*log10(p2(i0)/max(p2)));
bw = rad2deg(0.886/(N*d0*cosd(theta0)));
fprintf('  larghezza del fascio a f0: %.2f deg -> banda utile con sfasatori ~ %.1f %%\n', bw, deg2rad(bw)/tand(theta0)*100);

fig = figure('Position', [100 100 800 420]); hold on;
plot(rapporti, dir_ps, 'o-', 'Color', [0.85 0.2 0.2], 'DisplayName', 'sfasatori');
plot(rapporti, dir_ttd, 's-', 'Color', [0 0.45 0.74], 'DisplayName', 'ritardi veri (TTD)');
plot(rapporti, asind(u0./rapporti), 'k--', 'DisplayName', 'teoria sin\theta = (f_0/f) sin\theta_0');
xlabel('f / f_0'); ylabel('direzione del fascio [gradi]'); grid on; legend;
title('Con gli sfasatori il fascio si sposta con la frequenza');
fprintf(' figura: %s\n', salva_fig(fig, 'p05_squint'));

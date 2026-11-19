%% simula_beamformer - Diagramma Sigma e Delta misurato con il modello Simulink
% Punta il beamformer_ttd a theta0 = 20 deg e fa variare la direzione dell'onda
% in arrivo da -90 a +90 deg. Per ogni direzione misura l'ampiezza delle uscite
% Sigma e Delta (ultimo secondo di simulazione) e la confronta con la teoria:
%   |Sigma| = |sin(N psi/2) / (N sin(psi/2))|
%   |Delta| = 2 sin^2(N psi/4) / (N |sin(psi/2)|),   psi = pi (sin theta - sin theta0)
cartella = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(cartella), 'matlab', 'utils'));
mdl = 'beamformer_ttd';
if ~exist(fullfile(cartella, [mdl '.slx']), 'file'), run(fullfile(cartella, 'build_models.m')); end
load_system(fullfile(cartella, [mdl '.slx']));

N = 8; theta0 = 20;
thetas = -90:1:90;
S = zeros(size(thetas)); D = S;
for k = 1:numel(thetas)
    in = Simulink.SimulationInput(mdl);
    in = in.setVariable('theta0', theta0);
    in = in.setVariable('theta_src', thetas(k));
    out = sim(in);
    s = out.sigma_out; dlt = out.delta_out;
    S(k) = max(abs(s(end-100:end)));             % ultimo periodo (1 s)
    D(k) = max(abs(dlt(end-100:end)));
end
psi = pi*(sind(thetas) - sind(theta0));
S_teo = abs(sin(N*psi/2) ./ (N*sin(psi/2))); S_teo(abs(psi) < 1e-9) = 1;
D_teo = abs(2*sin(N*psi/4).^2 ./ (N*sin(psi/2))); D_teo(abs(psi) < 1e-9) = 0;

fprintf('simula_beamformer - Simulink, N = 8, theta0 = %d deg\n', theta0);
[~, im] = max(S);
riga('direzione del massimo di Sigma', thetas(im), theta0, 'deg');
[~, iz] = min(abs(thetas - theta0));
riga('Delta nella direzione di puntamento', D(iz), 0);
riga('errore massimo |Sigma| simulato - teorico', max(abs(S - S_teo)), 0);
riga('errore massimo |Delta| simulato - teorico', max(abs(D - D_teo)), 0);

fig = figure('Position', [100 100 900 420]); hold on;
plot(thetas, 20*log10(S + 1e-6), 'o', 'MarkerSize', 4, 'DisplayName', '\Sigma Simulink');
plot(thetas, 20*log10(S_teo + 1e-6), '-', 'DisplayName', '\Sigma teoria');
plot(thetas, 20*log10(D + 1e-6), 's', 'MarkerSize', 4, 'DisplayName', '\Delta Simulink');
plot(thetas, 20*log10(D_teo + 1e-6), '-', 'DisplayName', '\Delta teoria');
xline(theta0, ':', 'HandleVisibility', 'off');
ylim([-40 2]); grid on; legend('Location', 'southwest');
xlabel('direzione dell''onda in arrivo \theta [gradi]'); ylabel('ampiezza [dB]');
title('Beamformer a ritardo vero in Simulink: fasci somma e differenza');
fprintf(' figura: %s\n', salva_fig(fig, 'simulink_beamformer'));

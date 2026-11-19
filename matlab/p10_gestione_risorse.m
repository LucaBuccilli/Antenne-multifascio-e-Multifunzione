%% p10 - Radar multifunzione: gestione del tempo tra sorveglianza e inseguimento
% 420 posizioni di sorveglianza da 4 ms; Nt tracce (dwell 2 ms; 3 su 4 ogni 1 s,
% 1 su 4 in manovra ogni 0,25 s). Scheduler a priorita' con degrado controllato.
%   carico L = sum(T_dwell/T_update),  rivisita T_s = N_beam T_dwell_s / (1 - L)
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

par.N_BEAM = 420; par.T_SURV = 4e-3; par.T_TRK = 2e-3; par.DURATA = 30;
fprintf('p10 - Gestione delle risorse (420 posizioni di sorveglianza da 4 ms)\n');
nts = 0:50:300;
Ts_mis = nan(size(nts)); Ts_teo = nan(size(nts)); deg_nt = []; deg_Ts = [];
for k = 1:numel(nts)
    nt = nts(k);
    [Lc, Ts, rit] = simula(nt, false, par);
    Ts_mis(k) = Ts;
    fprintf(' %3d tracce: carico di inseguimento %5.1f %%\n', nt, 100*Lc);
    if Lc < 1
        Ts_teo(k) = par.N_BEAM*par.T_SURV/(1 - Lc);
        riga('  tempo di rivisita della sorveglianza', Ts, Ts_teo(k), 's');
        fprintf('    ritardo degli aggiornamenti: medio %.2f ms, massimo %.2f ms\n', 1e3*mean(rit), 1e3*max(rit));
        fprintf('    portata di sorveglianza a rivisita fissa: %.1f %% di quella senza tracce\n', 100*(1 - Lc)^0.25);
    else
        fprintf('    SATURAZIONE: la sorveglianza non completa la scansione, ritardo massimo tracce %.0f ms\n', 1e3*max(rit));
    end
    if Lc > 0.8
        [Ld, Tsd, ritd] = simula(nt, true, par);
        deg_nt(end+1) = nt; deg_Ts(end+1) = Tsd; %#ok<SAGROW>
        upd = periodi(nt, true, par);
        fprintf('    con degrado controllato: carico %.1f %%, rivisita %.2f s, ritardo massimo %.1f ms (tracce lente ogni %.2f s)\n', ...
            100*Ld, Tsd, 1e3*max(ritd), max(upd));
    end
end

fig = figure('Position', [100 100 800 420]); hold on;
ok = isfinite(Ts_mis);
plot(nts(ok), Ts_mis(ok), 'o', 'DisplayName', 'simulazione');
plot(nts, Ts_teo, '-', 'DisplayName', 'teoria N T_d/(1-L)');
plot(deg_nt, deg_Ts, 's', 'DisplayName', 'con degrado controllato');
ylim([0 15]); grid on; legend('Location', 'northwest');
xlabel('numero di tracce'); ylabel('rivisita della sorveglianza [s]'); title('Più tracce = sorveglianza più lenta');
fprintf(' figura: %s\n', salva_fig(fig, 'p10_gestione_risorse'));

%% funzioni locali
function upd = periodi(nt, degrado, par)
upd = ones(1, nt);
upd(mod(0:nt-1, 4) == 3) = 0.25;
if nt == 0, return; end
Lc = sum(par.T_TRK ./ upd);
if degrado && Lc > 0.8
    man = upd < 1;
    L_man = sum(par.T_TRK ./ upd(man));
    upd(~man) = par.T_TRK * nnz(~man) / max(0.7 - L_man, 1e-3);
end
end

function [Lc, Ts, rit] = simula(nt, degrado, par)
upd = periodi(nt, degrado, par);
rng(nt);
due = rand(1, nt) .* upd;
t = 0; beam = 0; scansioni = []; inizio = 0; rit = [];
while t < par.DURATA
    if nt > 0 && min(due) <= t
        [~, i] = min(due - t);                           % traccia piu' in ritardo
        rit(end+1) = t - due(i); %#ok<AGROW>
        t = t + par.T_TRK;
        due(i) = due(i) + upd(i);
    else
        t = t + par.T_SURV;
        beam = beam + 1;
        if beam == par.N_BEAM
            scansioni(end+1) = t - inizio; %#ok<AGROW>
            beam = 0; inizio = t;
        end
    end
end
if nt > 0, Lc = sum(par.T_TRK ./ upd); else, Lc = 0; end
if isempty(scansioni), Ts = inf; else, Ts = mean(scansioni); end
if isempty(rit), rit = 0; end
end

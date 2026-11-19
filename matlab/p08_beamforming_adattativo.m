%% p08 - Beamforming adattativo MVDR e regola di Reed-Mallett-Brennan
% N = 16 elementi a lambda/2, bersaglio a 0 deg, disturbatori a +35 e -50 deg (JNR 40 dB).
%   w = R^-1 a / (a^H R^-1 a)          perdita SMI media = (K + 2 - N)/(K + 1)
clear; close all;
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

N = 16; d = 0.5; n = 0:N-1;
rng(3);
a_t = exp(1j*2*pi*d*sind(0)*n).';
jammers = [35 -50];
jnr = 10^(40/10);
A_j = exp(1j*2*pi*d*sind(jammers(:))*n);            % 2 x N
R = eye(N) + jnr * (A_j.' * conj(A_j));
w_opt = (R \ a_t) / (a_t' * (R \ a_t));
w_conv = a_t / N;
sinr = @(w) abs(w'*a_t)^2 / real(w'*R*w);
sinr_opt = real(a_t' * (R \ a_t));

fprintf('p08 - Beamforming adattativo (N = 16, 2 disturbatori a 40 dB)\n');
riga('SINR ottimo (MVDR) rispetto a N', 10*log10(sinr(w_opt)/N), 10*log10(sinr_opt/N), 'dB');
fprintf('  %-42s misurato = %10.3f dB\n', 'SINR beamformer convenzionale (rispetto a N)', 10*log10(sinr(w_conv)/N));
u = linspace(-1, 1, 20001);
p_opt = 20*log10(abs(af(conj(w_opt), u, d)) / abs(w_opt'*a_t) + 1e-12);
p_conv = 20*log10(abs(af(conj(w_conv), u, d)) / abs(w_conv'*a_t) + 1e-12);
for th = jammers
    [~, i] = min(abs(u - sind(th)));
    fprintf('  guadagno verso il disturbatore a %+d deg: convenzionale %.1f dB, MVDR %.1f dB\n', th, p_conv(i), p_opt(i));
end

fprintf(' SMI - perdita media di SINR (200 prove per K)\n');
L = chol(R, 'lower');
Ks = [N, round(1.5*N), 2*N, 4*N, 8*N];
perd_mis = zeros(size(Ks)); perd_teo = perd_mis;
for k = 1:numel(Ks)
    K = Ks(k); rho = zeros(1, 200);
    for prova = 1:200
        Z = L * ((randn(N, K) + 1j*randn(N, K)) / sqrt(2));
        Rh = Z*Z' / K;
        w = Rh \ a_t;
        rho(prova) = sinr(w) / sinr_opt;
    end
    perd_mis(k) = 10*log10(mean(rho));
    perd_teo(k) = 10*log10((K + 2 - N)/(K + 1));
    riga(sprintf('K = %3d campioni', K), perd_mis(k), perd_teo(k), 'dB');
end

fig = figure('Position', [100 100 1100 420]);
subplot(1, 2, 1); plot(asind(u), p_conv, 'Color', [0.6 0.6 0.6]); hold on; plot(asind(u), p_opt, 'LineWidth', 1.4);
xline(jammers(1), 'r--'); xline(jammers(2), 'r--'); ylim([-80 5]); grid on;
xlabel('\theta [gradi]'); ylabel('guadagno [dB]'); legend('convenzionale', 'MVDR'); title('Il beamformer MVDR mette nulli sui disturbatori');
subplot(1, 2, 2); plot(Ks, perd_mis, 'o', Ks, perd_teo, '-'); grid on; xline(2*N, ':');
xlabel('campioni K'); ylabel('perdita di SINR [dB]'); legend('Monte Carlo', 'teoria (K+2-N)/(K+1)'); title('Regola RMB: K = 2N -> circa -3 dB');
fprintf(' figura: %s\n', salva_fig(fig, 'p08_mvdr'));

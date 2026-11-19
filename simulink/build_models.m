%% build_models - Costruisce il modello Simulink "beamformer_ttd" con i soli blocchi base
% Array di 8 elementi a lambda/2 che riceve un'onda piana sinusoidale (periodo 1 s,
% quindi lambda/c = 1 s) dalla direzione theta_src. Ogni elemento e' una Sine Wave
% con la fase dovuta al ritardo di arrivo; un Transport Delay per elemento compensa
% il ritardo verso la direzione di puntamento theta0 (ritardo vero, non sfasatore).
% Due blocchi Sum formano:
%   Sigma = somma di tutti gli elementi           (fascio somma)
%   Delta = meta' sinistra - meta' destra          (fascio differenza del monopulse)
% Variabili richieste nel workspace base: theta0, theta_src (gradi).
cartella = fileparts(mfilename('fullpath'));
mdl = 'beamformer_ttd';
N = 8;
if bdIsLoaded(mdl), close_system(mdl, 0); end
new_system(mdl);
open_system(mdl);
assignin('base', 'theta0', 0);
assignin('base', 'theta_src', 0);

for n = 0:N-1
    y = 30 + n*60;
    add_block('simulink/Sources/Sine Wave', sprintf('%s/Elemento %d', mdl, n), ...
        'Amplitude', '1', 'Frequency', '2*pi', ...
        'Phase', sprintf('-2*pi*0.5*%d*sind(theta_src)', n), ...
        'Position', [30 y 70 y+30]);
    add_block('simulink/Continuous/Transport Delay', sprintf('%s/Ritardo %d', mdl, n), ...
        'DelayTime', sprintf('0.5*(max(sind(theta0),0)*(%d-1-%d) + max(-sind(theta0),0)*%d)', N, n, n), ...
        'InitialOutput', '0', 'Position', [140 y 200 y+30]);
    add_line(mdl, sprintf('Elemento %d/1', n), sprintf('Ritardo %d/1', n));
end

add_block('simulink/Math Operations/Sum', [mdl '/Somma'], 'Inputs', repmat('+', 1, N), ...
    'IconShape', 'rectangular', 'Position', [300 40 330 480]);
add_block('simulink/Math Operations/Sum', [mdl '/Differenza'], ...
    'Inputs', [repmat('+', 1, N/2) repmat('-', 1, N/2)], ...
    'IconShape', 'rectangular', 'Position', [300 520 330 960]);
for n = 0:N-1
    add_line(mdl, sprintf('Ritardo %d/1', n), sprintf('Somma/%d', n+1), 'autorouting', 'on');
    add_line(mdl, sprintf('Ritardo %d/1', n), sprintf('Differenza/%d', n+1), 'autorouting', 'on');
end

add_block('simulink/Math Operations/Gain', [mdl '/Norm Sigma'], 'Gain', sprintf('1/%d', N), 'Position', [390 245 430 275]);
add_block('simulink/Math Operations/Gain', [mdl '/Norm Delta'], 'Gain', sprintf('1/%d', N), 'Position', [390 725 430 755]);
add_line(mdl, 'Somma/1', 'Norm Sigma/1');
add_line(mdl, 'Differenza/1', 'Norm Delta/1');

add_block('simulink/Sinks/To Workspace', [mdl '/sigma_out'], 'VariableName', 'sigma_out', ...
    'SaveFormat', 'Array', 'Position', [500 240 570 280]);
add_block('simulink/Sinks/To Workspace', [mdl '/delta_out'], 'VariableName', 'delta_out', ...
    'SaveFormat', 'Array', 'Position', [500 720 570 760]);
add_block('simulink/Sinks/Scope', [mdl '/Scope'], 'NumInputPorts', '2', 'Position', [500 450 540 530]);
add_line(mdl, 'Norm Sigma/1', 'sigma_out/1');
add_line(mdl, 'Norm Delta/1', 'delta_out/1');
add_line(mdl, 'Norm Sigma/1', 'Scope/1', 'autorouting', 'on');
add_line(mdl, 'Norm Delta/1', 'Scope/2', 'autorouting', 'on');

set_param(mdl, 'Solver', 'ode3', 'FixedStep', '0.01', 'StopTime', '10');
save_system(mdl, fullfile(cartella, [mdl '.slx']));
fprintf('Modello %s creato in %s\n', mdl, cartella);

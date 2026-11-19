function run_all()
%RUN_ALL Esegue tutte le simulazioni MATLAB in ordine.
%   Ogni script gira nel proprio spazio di lavoro, cosi' il "clear" iniziale
%   degli script non cancella le variabili di questo ciclo.
cartella = fileparts(mfilename('fullpath'));
script = dir(fullfile(cartella, 'p*.m'));
t0 = tic;
for k = 1:numel(script)
    fprintf('%s\n', repmat('=', 1, 90));
    esegui(fullfile(cartella, script(k).name));
end
fprintf('%s\nCompletato in %.1f s\n', repmat('=', 1, 90), toc(t0));
end

function esegui(file)
run(file);
end

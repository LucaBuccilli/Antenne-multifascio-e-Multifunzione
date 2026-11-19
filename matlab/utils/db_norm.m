function p = db_norm(x, floor_db)
%DB_NORM Modulo in dB normalizzato al massimo, con pavimento (default -80 dB).
if nargin < 2, floor_db = -80; end
a = abs(x); a = a / max(a);
p = max(20*log10(max(a, 1e-12)), floor_db);
end

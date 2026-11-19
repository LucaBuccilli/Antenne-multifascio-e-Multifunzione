function sll = lobo_laterale_max(p)
%LOBO_LATERALE_MAX Livello del lobo laterale piu' alto rispetto al massimo (dB).
p = p(:); [pk, i0] = max(p);
l = i0; while l > 1 && p(l-1) < p(l), l = l - 1; end
r = i0; while r < numel(p) && p(r+1) < p(r), r = r + 1; end
resto = [p(1:l-1); p(r+1:end)];
sll = max(resto) - pk;
end

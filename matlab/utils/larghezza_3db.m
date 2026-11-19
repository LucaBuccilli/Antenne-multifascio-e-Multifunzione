function bw = larghezza_3db(x, p)
%LARGHEZZA_3DB Larghezza del lobo principale a -3 dB (interpolata), nelle unita' di x.
x = x(:); p = p(:);
[pk, i0] = max(p); lvl = pk - 3;
i = i0; while i > 1 && p(i) > lvl, i = i - 1; end
xl = interp1([p(i) p(i+1)], [x(i) x(i+1)], lvl);
j = i0; while j < numel(p) && p(j) > lvl, j = j + 1; end
xr = interp1([p(j) p(j-1)], [x(j) x(j-1)], lvl);
bw = xr - xl;
end

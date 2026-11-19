function y = af(w, u, d)
%AF Fattore di array AF(u) = sum_n w_n exp(j 2 pi d n u), d in lunghezze d'onda.
n = 0:numel(w)-1;
y = exp(1j*2*pi*d*u(:)*n) * w(:);
end

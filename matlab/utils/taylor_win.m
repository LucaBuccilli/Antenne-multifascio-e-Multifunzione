function w = taylor_win(N, nbar, sll)
%TAYLOR_WIN Pesi di Taylor (stesso algoritmo di scipy.signal.windows.taylor, norm=False).
%   N elementi, nbar lobi quasi costanti, sll livello dei lobi laterali in dB (positivo).
B = 10^(sll/20);
A = acosh(B)/pi;
s2 = nbar^2 / (A^2 + (nbar - 0.5)^2);
ma = 1:nbar-1;
Fm = zeros(1, nbar-1);
for mi = 1:nbar-1
    m = ma(mi);
    segno = (-1)^(mi+1);
    numer = segno * prod(1 - m^2/s2 ./ (A^2 + (ma - 0.5).^2));
    altri = ma([1:mi-1, mi+1:end]);
    denom = 2 * prod(1 - m^2 ./ altri.^2);
    Fm(mi) = numer / denom;
end
n = 0:N-1;
w = 1 + 2 * (Fm * cos(2*pi*ma(:) * (n - N/2 + 0.5) / N));
w = w(:);
end

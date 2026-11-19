function w = cheb_win(N, at)
%CHEB_WIN Finestra di Dolph-Chebyshev con lobi laterali a -at dB
%   (stesso algoritmo di scipy.signal.windows.chebwin), senza toolbox.
order = N - 1;
beta = cosh(acosh(10^(at/20)) / order);
k = 0:N-1;
x = beta * cos(pi*k/N);
p = zeros(1, N);
i1 = x > 1;  p(i1) = cosh(order * acosh(x(i1)));
i2 = x < -1; p(i2) = (2*mod(N,2) - 1) * cosh(order * acosh(-x(i2)));
i3 = abs(x) <= 1; p(i3) = cos(order * acos(x(i3)));
if mod(N, 2)
    w = real(fft(p));
    n = (N + 1) / 2;
    w = w(1:n);
    w = [w(n:-1:2), w];
else
    p = p .* exp(1j*pi/N*(0:N-1));
    w = real(fft(p));
    n = N/2 + 1;
    w = [w(n:-1:2), w(2:n)];
end
w = w(:) / max(w);
end

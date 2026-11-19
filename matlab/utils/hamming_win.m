function w = hamming_win(N)
%HAMMING_WIN Finestra di Hamming simmetrica (senza toolbox).
n = (0:N-1)';
w = 0.54 - 0.46*cos(2*pi*n/(N-1));
end

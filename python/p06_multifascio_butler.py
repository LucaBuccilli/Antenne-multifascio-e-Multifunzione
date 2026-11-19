"""p06 - Multifascio: fasci ortogonali (matrice di Butler = FFT spaziale).

Array di N = 8 elementi a lambda/2. La matrice di Butler genera N fasci con
sin(theta_m) = (m - 1/2) / (N d), m = -N/2+1 ... N/2. Si verifica:
  - posizione dei fasci
  - ortogonalita': il massimo di ogni fascio cade sugli zeri degli altri
  - livello di incrocio tra fasci adiacenti: -3,92 dB
  - lo stesso insieme di fasci si ottiene con una FFT dei segnali degli elementi
"""
import numpy as np
import matplotlib.pyplot as plt
from utils import fattore_array, riga, salva

N, d = 8, 0.5
u = np.linspace(-1, 1, 40001)
ms = np.arange(-N // 2 + 1, N // 2 + 1)
u_beams = (ms - 0.5) / (N * d)
n = np.arange(N)

print("p06 - Multifascio di Butler (N = 8, d = lambda/2)")
fig, ax = plt.subplots(figsize=(9, 4.5))
patterns = []
for m, ub in zip(ms, u_beams):
    w = np.exp(-1j * 2 * np.pi * d * ub * n) / N
    p = np.abs(fattore_array(w, u, d))
    patterns.append(p)
    ax.plot(np.degrees(np.arcsin(u)), 20 * np.log10(p + 1e-9), lw=1.2, label=f"fascio {m}")
patterns = np.array(patterns)
for m, ub, p in zip(ms, u_beams, patterns):
    riga(f"fascio {m:+d}: direzione", np.degrees(np.arcsin(u[np.argmax(p)])), np.degrees(np.arcsin(ub)), "deg")

# ortogonalita': guadagno del fascio k nella direzione del fascio j
G = np.array([[abs(fattore_array(np.exp(-1j * 2 * np.pi * d * uk * n) / N, np.array([uj]), d)[0]) for uj in u_beams] for uk in u_beams])
riga("massimo fuori diagonale (ortogonalita')", 20 * np.log10(np.max(G - np.diag(np.diag(G))) + 1e-15), -300, "dB")
# incrocio tra fasci adiacenti (a meta' strada)
uc = (u_beams[3] + u_beams[4]) / 2
w = np.exp(-1j * 2 * np.pi * d * u_beams[3] * n) / N
riga("incrocio tra fasci adiacenti (teoria N grande)", 20 * np.log10(abs(fattore_array(w, np.array([uc]), d)[0])), -3.92, "dB")

# equivalenza con la FFT: un'onda piana da u = u_beams[5] eccita un solo canale della FFT
x = np.exp(1j * 2 * np.pi * d * u_beams[5] * n)
y = np.abs(np.fft.fftshift(np.fft.fft(x * np.exp(1j * np.pi * n / N)))) / N   # mezzo bin di offset come la Butler
print(f"  FFT: uscite dei canali per un'onda dal fascio {ms[5]:+d}: {np.round(y, 3)}")

ax.set(xlabel="θ [gradi]", ylabel="guadagno [dB]", ylim=(-30, 1), xlim=(-90, 90),
       title="8 fasci simultanei ortogonali: ognuno ha il massimo sugli zeri degli altri")
ax.axhline(-3.92, ls="--", color="0.5", lw=0.8)
ax.grid(alpha=0.3); ax.legend(fontsize=7, ncol=4)
print(" figura:", salva(fig, "p06_multifascio"))

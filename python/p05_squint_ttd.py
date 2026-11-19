"""p05 - Squint di frequenza: sfasatori contro ritardi veri (TTD).

Array di N = 64 elementi a d = lambda0/2, puntato a theta0 = 45 deg alla
frequenza f0. Si varia la frequenza tra 0,9 f0 e 1,1 f0.
  - con sfasatori (fase costante): sin(theta) = (f0/f) sin(theta0)
  - con ritardi veri: la direzione resta theta0 a tutte le frequenze
Si misura anche la perdita di guadagno in theta0 al bordo di una banda del 4%.
"""
import numpy as np
import matplotlib.pyplot as plt
from utils import riga, salva

N, theta0 = 64, 45.0
d0 = 0.5                        # spaziatura in lunghezze d'onda alla frequenza f0
n = np.arange(N)
u = np.linspace(-1, 1, 200001)
u0 = np.sin(np.radians(theta0))
rapporti = np.linspace(0.9, 1.1, 21)          # f / f0


def pattern(fr, ttd):
    d = d0 * fr                                 # la spaziatura elettrica cresce con f
    if ttd:
        fase = -2 * np.pi * d * u0 * n          # ritardo vero: fase proporzionale a f
    else:
        fase = -2 * np.pi * d0 * u0 * n         # sfasatore: fase fissa calcolata a f0
    a = np.exp(1j * 2 * np.pi * d * np.outer(u, n)) @ np.exp(1j * fase)
    return np.abs(a)


print("p05 - Squint (N = 64, theta0 = 45 deg)")
dir_ps, dir_ttd = [], []
for fr in rapporti:
    dir_ps.append(np.degrees(np.arcsin(u[np.argmax(pattern(fr, False))])))
    dir_ttd.append(np.degrees(np.arcsin(u[np.argmax(pattern(fr, True))])))
for fr in [0.95, 1.05]:
    teo = np.degrees(np.arcsin(u0 / fr))
    riga(f"direzione con sfasatori a f = {fr:.2f} f0", dir_ps[list(np.round(rapporti, 3)).index(round(fr, 3))], teo, "deg")
    riga(f"direzione con TTD a f = {fr:.2f} f0", dir_ttd[list(np.round(rapporti, 3)).index(round(fr, 3))], theta0, "deg")
# formula linearizzata
riga("squint per df/f0 = +2% (linearizzato)", dir_ps[12] - theta0, -np.degrees(0.02 * np.tan(np.radians(theta0))), "deg")

# perdita in theta0 al bordo di una banda del 4% (f = 1,02 f0)
i0 = np.argmin(np.abs(u - u0))
p_ref = pattern(1.0, False)[i0]
print(f"  perdita in theta0 a f = 1,02 f0 con sfasatori: {20*np.log10(pattern(1.02, False)[i0] / pattern(1.02, False).max()):.2f} dB")
print(f"  perdita in theta0 a f = 1,02 f0 con TTD:       {20*np.log10(pattern(1.02, True)[i0] / pattern(1.02, True).max()):.2f} dB")
bw = np.degrees(0.886 / (N * d0 * np.cos(np.radians(theta0))))
print(f"  larghezza del fascio a f0: {bw:.2f} deg -> banda utile con sfasatori ~ theta3dB/tan(theta0) = {np.radians(bw)/np.tan(np.radians(theta0))*100:.1f} %")

fig, ax = plt.subplots(figsize=(8, 4.2))
ax.plot(rapporti, dir_ps, "o-", color="C3", label="sfasatori")
ax.plot(rapporti, dir_ttd, "s-", color="C0", label="ritardi veri (TTD)")
ax.plot(rapporti, np.degrees(np.arcsin(u0 / rapporti)), "k--", lw=0.8, label="teoria sin θ = (f0/f) sin θ0")
ax.set(xlabel="f / f0", ylabel="direzione del fascio [gradi]", title="Con gli sfasatori il fascio si sposta con la frequenza")
ax.grid(alpha=0.3); ax.legend()
print(" figura:", salva(fig, "p05_squint"))

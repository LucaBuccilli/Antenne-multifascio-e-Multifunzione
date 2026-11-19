"""p03 - Grating lobes e spaziatura degli elementi.

Array di 16 elementi puntato a theta0 = 50 deg. Si varia la spaziatura d.
Teoria: massimi in u_m = u0 + m / d (d in lambda). Il grating lobe entra
nello spazio visibile se d > 1 / (1 + sin theta0) = 0.566 lambda.
"""
import numpy as np
import matplotlib.pyplot as plt
from utils import fattore_array, db, riga, salva

N, theta0 = 16, 50.0
u0 = np.sin(np.radians(theta0))
u = np.linspace(-1, 1, 40001)
d_max = 1 / (1 + u0)

print("p03 - Grating lobes (N = 16, theta0 = 50 deg)")
print(f"  spaziatura massima senza grating lobes: d < {d_max:.3f} lambda")
fig, axs = plt.subplots(2, 2, figsize=(10, 6), sharey=True)
for ax, d in zip(axs.ravel(), [0.5, 0.566, 0.7, 1.0]):
    w = np.exp(-1j * 2 * np.pi * d * u0 * np.arange(N))
    p = db(fattore_array(w, u, d))
    # grating lobe teorico m = -1
    ug = u0 - 1 / d
    visibile = abs(ug) <= 1
    print(f" d = {d:.3f} lambda: grating lobe teorico in u = {ug:+.3f} -> {'VISIBILE' if visibile else 'non visibile'}")
    if visibile:
        # cerca il massimo nella zona del grating lobe
        mask = np.abs(u - ug) < 0.1
        ug_mis = u[mask][np.argmax(p[mask])]
        riga("posizione grating lobe (gradi)", np.degrees(np.arcsin(ug_mis)), np.degrees(np.arcsin(ug)), "deg")
        riga("livello grating lobe", p[mask].max(), 0.0, "dB")
    ax.plot(np.degrees(np.arcsin(u)), p, lw=1.2)
    ax.set_title(f"d = {d} λ" + ("  (grating lobe!)" if visibile and abs(ug) < 0.999 else ""))
    ax.set_ylim(-40, 2); ax.grid(alpha=0.3)
for ax in axs[1]:
    ax.set_xlabel("θ [gradi]")
fig.suptitle("Puntamento a 50°: oltre 0,566 λ compare un secondo lobo alto quanto il principale")
fig.tight_layout()
print(" figura:", salva(fig, "p03_grating_lobes"))

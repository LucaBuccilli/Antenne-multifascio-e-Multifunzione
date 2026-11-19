"""p01 - Fattore di array e puntamento elettronico.

Array lineare uniforme di N dipoli a d = lambda/2. Lo sfasamento progressivo
alpha = -k d sin(theta0) punta il fascio in theta0. Si verifica:
  - direzione del massimo = theta0
  - larghezza a -3 dB = 0.886 lambda / (N d cos theta0)
  - primo lobo laterale = -13.26 dB
  - guadagno con diagramma dell'elemento ~ cos(theta0) (perdita di scansione)
"""
import numpy as np
import matplotlib.pyplot as plt
from utils import fattore_array, db, larghezza_3db, lobo_laterale_max, riga, salva

N, d = 16, 0.5
theta = np.linspace(-90, 90, 36001)
u = np.sin(np.radians(theta))

print("p01 - Fattore di array e puntamento (N = 16, d = lambda/2)")
fig, ax = plt.subplots(figsize=(9, 4.5))
for theta0, col in [(0, "0.6"), (30, "C0"), (60, "C3")]:
    alpha = -2 * np.pi * d * np.sin(np.radians(theta0))      # sfasamento progressivo
    w = np.exp(1j * alpha * np.arange(N))                      # pesi = solo fase
    p = db(fattore_array(w, u, d))
    print(f" theta0 = {theta0} deg  (alpha = {np.degrees(alpha):.1f} deg tra elementi)")
    riga("direzione del massimo", theta[np.argmax(p)], theta0, "deg")
    bw_teo = np.degrees(0.886 / (N * d * np.cos(np.radians(theta0))))
    riga("larghezza a -3 dB", larghezza_3db(theta, p), bw_teo, "deg")
    riga("lobo laterale massimo (teoria per N grande)", lobo_laterale_max(p), -13.26, "dB")
    ax.plot(theta, p, color=col, lw=1.6, label=f"puntato a {theta0}°")

# perdita di scansione: con un diagramma dell'elemento ~ cos(theta) (dipolo davanti a
# un piano di massa, approssimato) il massimo del diagramma totale E_el * AF cala come cos(theta0)
print(" Perdita di scansione (diagramma dell'elemento cos(theta))")
p_ref = None
for theta0 in [0, 30, 45, 60]:
    w = np.exp(-1j * 2 * np.pi * d * np.sin(np.radians(theta0)) * np.arange(N))
    tot = np.cos(np.radians(theta)) * np.abs(fattore_array(w, u, d))
    if p_ref is None:
        p_ref = tot.max()
    riga(f"perdita a {theta0} deg (sola andata)", 20 * np.log10(tot.max() / p_ref),
         20 * np.log10(np.cos(np.radians(theta0))), "dB")

ax.set(xlabel="angolo θ dalla normale [gradi]", ylabel="|AF| normalizzato [dB]", ylim=(-50, 2),
       title="Puntamento elettronico: il fascio si sposta e si allarga")
ax.axhline(-13.26, ls="--", color="0.5", lw=0.8)
ax.grid(alpha=0.3); ax.legend()
print(" figura:", salva(fig, "p01_puntamento"))

"""p02 - Rastremazione delle ampiezze (tapering).

Array di N = 32 elementi a lambda/2. Confronto tra alimentazione uniforme,
Hamming, Taylor (-30 dB, nbar = 4) e Dolph-Chebyshev (-30 dB).
Per ogni finestra si misura: lobo laterale massimo, allargamento del fascio
rispetto all'uniforme, efficienza di rastremazione |sum w|^2 / (N sum w^2).
"""
import numpy as np
import matplotlib.pyplot as plt
import warnings
from scipy.signal.windows import taylor, chebwin, hamming
warnings.filterwarnings("ignore", message="This window is not suitable")
from utils import fattore_array, db, larghezza_3db, lobo_laterale_max, riga, salva

N, d = 32, 0.5
u = np.linspace(-1, 1, 40001)
finestre = {
    "Uniforme": (np.ones(N), -13.26),
    "Hamming": (hamming(N), -42.7),
    "Taylor -30 dB, nbar=4": (taylor(N, nbar=4, sll=30, norm=False), -30.0),
    "Chebyshev -30 dB": (chebwin(N, at=30), -30.0),
}

print("p02 - Rastremazione (N = 32, d = lambda/2)")
bw_uni = None
fig, ax = plt.subplots(figsize=(9, 4.5))
for nome, (w, sll_teo) in finestre.items():
    p = db(fattore_array(w, u, d))
    bw = larghezza_3db(u, p)
    if bw_uni is None:
        bw_uni = bw
    eff = abs(w.sum()) ** 2 / (N * (abs(w) ** 2).sum())
    print(f" {nome}")
    riga("lobo laterale massimo", lobo_laterale_max(p), sll_teo, "dB")
    print(f"  {'allargamento del fascio (rispetto uniforme)':<42s} misurato = {bw / bw_uni:>10.3f}")
    print(f"  {'efficienza di rastremazione':<42s} misurato = {eff:>10.3f}       perdita = {-10*np.log10(eff):.2f} dB")
    ax.plot(np.degrees(np.arcsin(u)), p, lw=1.3, label=nome)

ax.set(xlabel="angolo θ [gradi]", ylabel="|AF| [dB]", ylim=(-60, 2), xlim=(-90, 90),
       title="Lobi laterali piu' bassi = fascio piu' largo e meno guadagno")
ax.grid(alpha=0.3); ax.legend(fontsize=8)
print(" figura:", salva(fig, "p02_rastremazione"))

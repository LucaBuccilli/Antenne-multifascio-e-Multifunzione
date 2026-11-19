"""p04 - Sfasatori digitali a B bit: perdita di guadagno e lobi di quantizzazione.

Array di N = 64 elementi a lambda/2 puntato a theta0 = 10 deg, con ampiezze
Taylor -40 dB (cosi' i lobi laterali normali non coprono quelli spuri). Le fasi ideali
vengono arrotondate al passo 360/2^B. Si confronta:
  - perdita di guadagno media: teoria ~ (pi^2/3) * 2^(-2B)
  - lobo di quantizzazione piu' alto con errore periodico: teoria ~ 1/2^(2B) (-6 dB per bit)
  - lo stesso con una fase casuale nota aggiunta e poi sottratta (dithering):
    l'errore diventa casuale e il lobo spurio sparisce nel fondo dei lobi laterali
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal.windows import taylor
from utils import fattore_array, riga, salva

N, d, theta0 = 64, 0.5, 10.0
u0 = np.sin(np.radians(theta0))
u = np.linspace(-1, 1, 80001)
n = np.arange(N)
fase_ideale = -2 * np.pi * d * u0 * n
amp = taylor(N, nbar=6, sll=40, norm=False)   # rastremazione in ampiezza
rng = np.random.default_rng(1)
dither = rng.uniform(0, 2 * np.pi, N)          # fase fissa nota di ogni elemento (es. lunghezza dei cavi)


def quantizza(phi, B):
    passo = 2 * np.pi / 2 ** B
    return np.round(phi / passo) * passo


def picco_spurio(p_db, u, u0, esclusione):
    """Livello del lobo piu' alto fuori dal lobo principale (|u - u0| > esclusione)."""
    return p_db[np.abs(u - u0) > esclusione].max()


ideale = np.abs(fattore_array(amp * np.exp(1j * fase_ideale), u, d))
G0 = ideale.max()
escl = 8.0 / (N * d)                            # fuori dai lobi laterali vicini (Taylor)

print("p04 - Quantizzazione degli sfasatori (N = 64, theta0 = 10 deg)")
fig, ax = plt.subplots(figsize=(9, 4.5))
ax.plot(np.degrees(np.arcsin(u)), 20 * np.log10(ideale / G0 + 1e-12), color="0.7", lw=1, label="fasi ideali (Taylor -40 dB)")
for B in [3, 4, 5]:
    w_q = amp * np.exp(1j * quantizza(fase_ideale, B))
    # con dithering: il comando e' quantizzato su (fase + dither); la fase fisica del cavo (-dither) si somma
    w_d = amp * np.exp(1j * (quantizza(fase_ideale + dither, B) - dither))
    pq = np.abs(fattore_array(w_q, u, d))
    pd = np.abs(fattore_array(w_d, u, d))
    print(f" B = {B} bit (passo {360 / 2 ** B:.2f} deg)")
    perdita_teo = -10 * np.log10(1 - np.pi ** 2 / 3 * 2.0 ** (-2 * B))
    riga("perdita di guadagno con dithering", -20 * np.log10(pd.max() / G0), perdita_teo, "dB")
    dbq = 20 * np.log10(pq / pq.max() + 1e-12)
    dbd = 20 * np.log10(pd / pd.max() + 1e-12)
    riga("lobo di quantizzazione (errore periodico)", picco_spurio(dbq, u, u0, escl), -6.02 * B, "dB")
    print(f"  {'lobo spurio piu alto con dithering':<42s} misurato = {picco_spurio(dbd, u, u0, escl):>10.3f} dB")
    if B == 3:
        ax.plot(np.degrees(np.arcsin(u)), dbq, color="C3", lw=1.2, label="3 bit, errore periodico")
        ax.plot(np.degrees(np.arcsin(u)), dbd, color="C0", lw=1.0, label="3 bit, con fase casuale")
print("  nota: il livello esatto del lobo di quantizzazione dipende da theta0; -6B dB e' un valore indicativo (circa 1-2 dB).")
print("        Con il dithering il lobo spurio scompare e resta solo un fondo di lobi laterali piu' alto.")
ax.set(xlabel="θ [gradi]", ylabel="|AF| [dB]", ylim=(-60, 2), xlim=(-90, 90),
       title="Sfasatori a 3 bit: l'errore periodico crea lobi spuri, la fase casuale li disperde")
ax.grid(alpha=0.3); ax.legend(fontsize=8)
print(" figura:", salva(fig, "p04_quantizzazione"))

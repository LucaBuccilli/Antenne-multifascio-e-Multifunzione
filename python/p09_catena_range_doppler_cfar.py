"""p09 - Catena di elaborazione di un fascio: compressione, MTI, Doppler, CA-CFAR.

Radar in banda S (3 GHz). Un treno di M = 32 impulsi chirp (B = 5 MHz,
T = 10 us) illumina due bersagli in movimento e un clutter fisso molto forte.
Catena: filtro adattato -> cancellatore MTI a 2 impulsi -> FFT Doppler con
finestra di Chebyshev -> rivelatore quadratico -> CA-CFAR in distanza.
Verifiche:
  - risoluzione in distanza c/(2B) e lobi laterali di compressione (-13,3 dB)
  - distanza e velocita' misurate dei bersagli
  - probabilita' di falso allarme del CA-CFAR su solo rumore
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal.windows import chebwin
from utils import riga, salva, larghezza_3db

c, fc = 3e8, 3e9
lam = c / fc
B, T, fs = 5e6, 10e-6, 10e6
PRF, M = 4000.0, 32
Nr = int(fs / PRF)                      # campioni per PRI (2500 -> 37,5 km)
rng = np.random.default_rng(11)

t = np.arange(int(T * fs)) / fs
chirp = np.exp(1j * np.pi * (B / T) * t ** 2)
bersagli = [  # distanza [m], velocita' radiale [m/s], ampiezza per campione
    (12000.0, 40.0, 0.3),
    (20010.0, -25.0, 0.2),
]


def eco(con_bersagli=True, con_clutter=True):
    """Matrice M x Nr dei campioni ricevuti (prima della compressione)."""
    refl = np.zeros((M, Nr), dtype=complex)
    m = np.arange(M)[:, None]
    if con_bersagli:
        for R, v, a in bersagli:
            k = int(round(2 * R / c * fs))
            fd = 2 * v / lam
            refl[:, k] += a * np.exp(1j * 2 * np.pi * fd * m / PRF).ravel() * np.exp(1j * rng.uniform(0, 2 * np.pi))
    if con_clutter:
        celle = np.arange(150, Nr - 200)
        amp = np.sqrt(10 ** (30 / 10) / 2) * (rng.standard_normal(len(celle)) + 1j * rng.standard_normal(len(celle)))
        fd_c = rng.normal(0, 4.0, len(celle))                       # moto interno del clutter (Hz)
        refl[:, celle] += amp * np.exp(1j * 2 * np.pi * fd_c * m / PRF)
    rx = np.array([np.convolve(r, chirp)[:Nr] for r in refl])
    rx += (rng.standard_normal((M, Nr)) + 1j * rng.standard_normal((M, Nr))) / np.sqrt(2)
    return rx


def elabora(rx):
    mf = np.array([np.correlate(r, chirp, mode="full")[len(chirp) - 1:len(chirp) - 1 + Nr] for r in rx]) / np.sqrt(len(chirp))
    mti = mf[1:] - mf[:-1]                                           # cancellatore a 2 impulsi
    win = chebwin(M - 1, at=60)[:, None]
    rd = np.fft.fftshift(np.fft.fft(mti * win, axis=0), axes=0)
    return mf, np.abs(rd) ** 2


def ca_cfar(P, n_ref=32, n_guard=4, pfa=1e-6):
    alpha = n_ref * (pfa ** (-1 / n_ref) - 1)
    k = np.ones(n_ref // 2 + n_guard) ; k[-n_guard:] = 0
    kern = np.concatenate([k, [0], k[::-1]])
    rif = np.array([np.convolve(row, kern, mode="same") for row in P]) / n_ref
    det = P > alpha * rif
    det[:, : n_ref // 2 + n_guard] = False
    det[:, -(n_ref // 2 + n_guard):] = False
    return det, alpha


print("p09 - Catena range-Doppler e CA-CFAR (3 GHz, B = 5 MHz, M = 32)")
rx = eco()
mf, P = elabora(rx)

# 1) compressione: impulso compresso di un bersaglio isolato, chirp sovracampionato (x20) per misurare bene
fs_os = 20 * B
t_os = np.arange(int(T * fs_os)) / fs_os
ch_os = np.exp(1j * np.pi * (B / T) * t_os ** 2)
prof = np.abs(np.correlate(ch_os, ch_os, mode="full")) ** 2
prof_db = 10 * np.log10(prof / prof.max())
i0 = np.argmax(prof_db)
x_m = (np.arange(len(prof_db)) - i0) / fs_os * c / 2
riga("risoluzione in distanza (larghezza -3 dB)", larghezza_3db(x_m, prof_db), 0.886 * c / (2 * B), "m")
lob = prof_db.copy()
k1 = int(1.5 * fs_os / B)                                   # esclude il lobo principale
lob[i0 - k1:i0 + k1] = -200
riga("lobo laterale di compressione", lob.max(), -13.26, "dB")
g = np.abs(np.correlate(chirp, chirp, mode="full")).max() ** 2 / len(chirp)
riga("guadagno di compressione (campioni T*fs)", 10 * np.log10(g), 10 * np.log10(T * fs), "dB")
print(f"  rapporto di compressione B*T = {B*T:.0f}: l'impulso da {T*1e6:.0f} us diventa lungo circa {1/B*1e9:.0f} ns")

# 2) rivelazione
det, alpha = ca_cfar(P)
dv = lam * PRF / (2 * (M - 1))
vel = (np.arange(M - 1) - (M - 1) // 2) * PRF / (M - 1) * lam / 2
print(f"  soglia CA-CFAR alpha = {alpha:.2f} ({10*np.log10(alpha):.1f} dB), risoluzione in velocita' {dv:.2f} m/s")
for R, v, _ in bersagli:
    k = int(round(2 * R / c * fs))
    zona = (slice(None), slice(k - 5, k + 6))
    sub = np.where(det[zona], P[zona], 0)
    if sub.max() == 0:
        print(f"  bersaglio a {R/1e3:.2f} km NON rivelato")
        continue
    j, i = np.unravel_index(np.argmax(sub), sub.shape)
    riga(f"bersaglio {R/1e3:.2f} km: distanza", (k - 5 + i) / fs * c / 2, R, "m")
    riga(f"bersaglio {R/1e3:.2f} km: velocita' (entro {dv:.1f} m/s)", vel[j], v, "m/s")
n_det = det.sum()
print(f"  celle oltre soglia in totale: {n_det}: i due bersagli piu' qualche residuo di clutter vicino a velocita' zero")

# 3) falsi allarmi su solo rumore, Pfa nominale = 1e-3
fa, celle, fa2, celle2 = 0, 0, 0, 0
for _ in range(5):
    _, Pn = elabora(eco(False, False))
    dn, _ = ca_cfar(Pn, pfa=1e-3)
    fa += dn.sum(); celle += dn[:, 20:-20].size
    dn2, _ = ca_cfar(Pn[:, ::2], pfa=1e-3)                  # una cella si' e una no: celle indipendenti
    fa2 += dn2.sum(); celle2 += dn2[:, 20:-20].size
riga("Pfa con tutte le celle (x 1e-3)", 1e3 * fa / celle, 1.0)
riga("Pfa con celle indipendenti (x 1e-3)", 1e3 * fa2 / celle2, 1.0)
print("  le celle a passo 1/fs = 1/(2B) sono correlate: la media del CA-CFAR stima peggio il rumore e i falsi allarmi crescono")

fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
rng_km = np.arange(Nr) / fs * c / 2 / 1e3
axs[0].plot(rng_km, 20 * np.log10(np.abs(mf[0]) + 1e-9), lw=0.6, color="0.4")
axs[0].set(xlabel="distanza [km]", ylabel="ampiezza [dB]", title="Dopo la compressione: il clutter copre i bersagli")
axs[0].grid(alpha=0.3)
im = axs[1].imshow(10 * np.log10(P / np.median(P)), aspect="auto", origin="lower", cmap="viridis", vmin=0, vmax=40,
                   extent=[rng_km[0], rng_km[-1], vel[0] - dv / 2, vel[-1] + dv / 2])
dj, di = np.where(det)
axs[1].plot(rng_km[di], vel[dj], "r.", ms=4, label="rivelazioni CFAR")
axs[1].set(xlabel="distanza [km]", ylabel="velocità radiale [m/s]", title="Mappa range-Doppler dopo MTI e FFT")
axs[1].legend(loc="upper right", fontsize=8)
fig.colorbar(im, ax=axs[1], label="dB sul rumore")
print(" figura:", salva(fig, "p09_range_doppler"))

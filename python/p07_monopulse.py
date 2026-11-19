"""p07 - Monopulse con array: fasci somma e differenza, stima dell'angolo.

Array di N = 32 elementi a lambda/2 puntato alla normale.
  Sigma = somma delle due meta', Delta = meta' sinistra - meta' destra
  Rapporto monopulse: Im{-Delta/Sigma} = tan(N pi (u - u0) / 4)  ->  pendenza km = N pi / 4
Monte Carlo: bersaglio a u = 0,003 con rumore su ogni elemento; si confronta
la deviazione standard della stima con sigma_u = 1 / (km sqrt(2 SNR_Sigma)).
"""
import numpy as np
import matplotlib.pyplot as plt
from utils import riga, salva

N, d = 32, 0.5
n = np.arange(N)
km = N * np.pi / 4
sx = (n < N // 2).astype(float)          # meta' sinistra
dx = 1 - sx                               # meta' destra
rng = np.random.default_rng(7)


def canali(x):
    """Uscite somma e differenza per il vettore x dei segnali degli elementi."""
    return x.sum(axis=-1), (x * sx).sum(axis=-1) - (x * dx).sum(axis=-1)


# curva del rapporto monopulse senza rumore
uu = np.linspace(-0.04, 0.04, 801)
X = np.exp(1j * 2 * np.pi * d * np.outer(uu, n))
S, D = canali(X)
r = np.imag(-D / S)

print("p07 - Monopulse (N = 32, d = lambda/2)")
pend = np.polyfit(uu[np.abs(uu) < 0.003], r[np.abs(uu) < 0.003], 1)[0]
riga("pendenza monopulse km", pend, km)
riga("rapporto in u = 0,03 vs tan(N pi u / 4)", np.interp(0.03, uu, r), np.tan(N * np.pi * 0.03 / 4))

# Monte Carlo
u_t, prove = 0.003, 4000
print(f" Monte Carlo: bersaglio in u = {u_t}, {prove} prove per ogni SNR")
snr_list, sig_mis, sig_teo = [], [], []
for snr_db in [10, 15, 20, 25, 30]:
    snr_sigma = 10 ** (snr_db / 10)                  # SNR all'uscita del canale somma
    sigma_n = np.sqrt(N / snr_sigma)                 # rumore per elemento (ampiezza segnale = 1)
    s = np.exp(1j * 2 * np.pi * d * u_t * n)
    rum = sigma_n / np.sqrt(2) * (rng.standard_normal((prove, N)) + 1j * rng.standard_normal((prove, N)))
    S, D = canali(s + rum)
    u_hat = np.imag(-D / S) / km
    sig_mis.append(np.std(u_hat)); sig_teo.append(1 / (km * np.sqrt(2 * snr_sigma))); snr_list.append(snr_db)
    riga(f"sigma_u con SNR_Sigma = {snr_db} dB (x 1e-3)", 1e3 * sig_mis[-1], 1e3 * sig_teo[-1])
bw_u = 0.886 / (N * d)
print(f"  larghezza del fascio in u: {bw_u:.4f}; con SNR = 20 dB la precisione e' 1/{bw_u / sig_teo[2]:.0f} del fascio")

fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
axs[0].plot(uu, r, lw=1.6, label="Im{−Δ/Σ}")
axs[0].plot(uu, km * uu, "--", color="0.5", label="km · (u − u0)")
axs[0].set(xlabel="u − u0", ylabel="rapporto monopulse", ylim=(-3, 3), title="Rapporto monopulse: lineare vicino al centro")
axs[0].grid(alpha=0.3); axs[0].legend()
axs[1].semilogy(snr_list, sig_mis, "o", label="Monte Carlo")
axs[1].semilogy(snr_list, sig_teo, "-", label="teoria 1/(km √(2 SNR))")
axs[1].set(xlabel="SNR canale somma [dB]", ylabel="σ della stima di u", title="Precisione angolare")
axs[1].grid(alpha=0.3, which="both"); axs[1].legend()
print(" figura:", salva(fig, "p07_monopulse"))

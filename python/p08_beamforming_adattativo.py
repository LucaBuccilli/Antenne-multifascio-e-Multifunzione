"""p08 - Beamforming adattativo MVDR e regola di Reed-Mallett-Brennan.

Array di N = 16 elementi a lambda/2, bersaglio a 0 deg, due disturbatori a
+35 deg e -50 deg con rapporto disturbo/rumore 40 dB per elemento.
  - pesi ottimi: w = R^-1 a / (a^H R^-1 a)  -> nulli sui disturbatori
  - SMI: R stimata con K campioni; perdita media di SINR = (K + 2 - N) / (K + 1)
    (regola RMB: con K = 2N la perdita e' circa 3 dB)
"""
import numpy as np
import matplotlib.pyplot as plt
from utils import steering, fattore_array, riga, salva

N, d = 16, 0.5
rng = np.random.default_rng(3)
a_t = steering(N, np.sin(np.radians(0.0)), d)[0]
jammers = [35.0, -50.0]
jnr = 10 ** (40 / 10)
A_j = steering(N, np.sin(np.radians(jammers)), d)            # 2 x N
R = np.eye(N, dtype=complex) + jnr * (A_j.T @ A_j.conj())     # covarianza disturbo + rumore (sigma^2 = 1)

Ri = np.linalg.inv(R)
w_opt = Ri @ a_t / (a_t.conj() @ Ri @ a_t)
w_conv = a_t / N


def sinr(w):
    return abs(w.conj() @ a_t) ** 2 / np.real(w.conj() @ R @ w)


sinr_opt = np.real(a_t.conj() @ Ri @ a_t)
print("p08 - Beamforming adattativo (N = 16, 2 disturbatori a 40 dB)")
riga("SINR ottimo (pesi MVDR) rispetto a SNR di N elementi", 10 * np.log10(sinr(w_opt) / N), 10 * np.log10(sinr_opt / N), "dB")
print(f"  {'SINR beamformer convenzionale (rispetto a N)':<42s} misurato = {10*np.log10(sinr(w_conv) / N):>10.3f} dB")
u = np.linspace(-1, 1, 20001)
p_opt = 20 * np.log10(np.abs(fattore_array(w_opt.conj(), u, d)) / abs(w_opt.conj() @ a_t) + 1e-12)
p_conv = 20 * np.log10(np.abs(fattore_array(w_conv.conj(), u, d)) / abs(w_conv.conj() @ a_t) + 1e-12)
for th in jammers:
    i = np.argmin(np.abs(u - np.sin(np.radians(th))))
    print(f"  guadagno verso il disturbatore a {th:+.0f} deg: convenzionale {p_conv[i]:.1f} dB, MVDR {p_opt[i]:.1f} dB")

# SMI: perdita di SINR in funzione del numero di campioni K
print(" SMI - perdita media di SINR (200 prove per K)")
L = np.linalg.cholesky(R)
Ks, perd_mis, perd_teo = [N, int(1.5 * N), 2 * N, 4 * N, 8 * N], [], []
for K in Ks:
    rho = []
    for _ in range(200):
        Z = L @ ((rng.standard_normal((N, K)) + 1j * rng.standard_normal((N, K))) / np.sqrt(2))
        Rh = Z @ Z.conj().T / K
        w = np.linalg.solve(Rh, a_t)
        rho.append(sinr(w) / sinr_opt)
    perd_mis.append(10 * np.log10(np.mean(rho)))
    perd_teo.append(10 * np.log10((K + 2 - N) / (K + 1)))
    riga(f"K = {K:3d} campioni", perd_mis[-1], perd_teo[-1], "dB")

fig, axs = plt.subplots(1, 2, figsize=(11, 4.2))
th = np.degrees(np.arcsin(u))
axs[0].plot(th, p_conv, color="0.6", lw=1, label="convenzionale")
axs[0].plot(th, p_opt, color="C0", lw=1.4, label="MVDR")
for j in jammers:
    axs[0].axvline(j, color="C3", ls="--", lw=0.8)
axs[0].set(xlabel="θ [gradi]", ylabel="guadagno [dB]", ylim=(-80, 5), title="Il beamformer MVDR mette nulli sui disturbatori")
axs[0].grid(alpha=0.3); axs[0].legend()
axs[1].plot(Ks, perd_mis, "o", label="Monte Carlo")
axs[1].plot(Ks, perd_teo, "-", label="teoria (K+2−N)/(K+1)")
axs[1].axvline(2 * N, color="0.5", ls=":", lw=0.8)
axs[1].set(xlabel="campioni K", ylabel="perdita di SINR [dB]", title="Regola RMB: K = 2N → circa −3 dB")
axs[1].grid(alpha=0.3); axs[1].legend()
print(" figura:", salva(fig, "p08_mvdr"))

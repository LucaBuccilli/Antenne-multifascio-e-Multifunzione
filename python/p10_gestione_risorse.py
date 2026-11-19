"""p10 - Radar multifunzione: gestione del tempo tra sorveglianza e inseguimento.

Il radar ha un settore di sorveglianza di 420 posizioni di fascio (dwell 4 ms)
e Nt tracce da aggiornare (dwell 2 ms; 3 su 4 ogni 1 s, 1 su 4 in manovra ogni 0,25 s).
Se il carico supera l'80%, il radar allunga l'intervallo delle tracce non in
manovra (degrado controllato) per riportare il carico al 70%.
Scheduler a priorita': una traccia scaduta viene servita per prima (la piu'
in ritardo), altrimenti il radar esegue la prossima posizione di sorveglianza.
Teoria:
  carico di inseguimento  L = sum(T_dwell / T_update)
  tempo di rivisita della sorveglianza  T_s = N_beam * T_dwell_s / (1 - L)
  portata di sorveglianza a T_s fisso: R / R0 = (1 - L)^(1/4)
"""
import numpy as np
import matplotlib.pyplot as plt
from utils import riga, salva

N_BEAM, T_SURV = 420, 4e-3
T_TRK = 2e-3
DURATA = 30.0


def periodi(nt, degrado):
    upd = np.where(np.arange(nt) % 4 == 3, 0.25, 1.0)             # 1 traccia su 4 in manovra
    L = np.sum(T_TRK / upd) if nt else 0.0
    if degrado and L > 0.8:
        man = upd < 1
        L_man = np.sum(T_TRK / upd[man])
        upd[~man] = T_TRK * np.sum(~man) / max(0.7 - L_man, 1e-3)   # nuovo periodo per le tracce lente
    return upd


def simula(nt, registra=False, degrado=False):
    upd = periodi(nt, degrado)
    due = np.random.default_rng(nt).uniform(0, upd)              # fasi iniziali casuali
    t, beam, scansioni, inizio_scan = 0.0, 0, [], 0.0
    ritardi, log = [], []
    while t < DURATA:
        if nt and due.min() <= t:
            i = int(np.argmin(due - t))                           # traccia piu' in ritardo
            ritardi.append(t - due[i])
            if registra:
                log.append((t, T_TRK, "traccia in manovra" if upd[i] < 1 else "traccia"))
            t += T_TRK
            due[i] += upd[i]
        else:
            if registra:
                log.append((t, T_SURV, "sorveglianza"))
            t += T_SURV
            beam += 1
            if beam == N_BEAM:
                scansioni.append(t - inizio_scan)
                beam, inizio_scan = 0, t
    carico = np.sum(T_TRK / upd) if nt else 0.0
    return carico, (np.mean(scansioni) if scansioni else np.inf), (np.array(ritardi) if ritardi else np.zeros(1)), log


print("p10 - Gestione delle risorse (420 posizioni di sorveglianza da 4 ms)")
nts = [0, 50, 100, 150, 200, 250, 300]
cariche, Ts_mis, Ts_teo, Ts_deg = [], [], [], []
for nt in nts:
    L, Ts, rit, _ = simula(nt)
    cariche.append(L); Ts_mis.append(Ts)
    print(f" {nt:3d} tracce: carico di inseguimento {100*L:5.1f} %")
    if L < 1:
        teo = N_BEAM * T_SURV / (1 - L)
        Ts_teo.append(teo)
        riga("  tempo di rivisita della sorveglianza", Ts, teo, "s")
        print(f"    ritardo degli aggiornamenti di traccia: medio {1e3*rit.mean():.2f} ms, massimo {1e3*rit.max():.2f} ms")
        print(f"    portata di sorveglianza a rivisita fissa: {100*(1-L)**0.25:.1f} % di quella senza tracce")
    else:
        Ts_teo.append(np.nan)
        print(f"    SATURAZIONE: la sorveglianza non completa mai la scansione, ritardo massimo delle tracce {1e3*rit.max():.0f} ms")
    if L > 0.8:
        Ld, Tsd, ritd, _ = simula(nt, degrado=True)
        Ts_deg.append((nt, Tsd))
        print(f"    con degrado controllato: carico {100*Ld:.1f} %, rivisita {Tsd:.2f} s, ritardo massimo tracce {1e3*ritd.max():.1f} ms")
        print(f"      (le tracce non in manovra passano da 1 s a {periodi(nt, True).max():.2f} s)")

fig, axs = plt.subplots(1, 2, figsize=(12, 4.3))
ok = np.isfinite(Ts_mis)
axs[0].plot(np.array(nts)[ok], np.array(Ts_mis)[ok], "o", label="simulazione")
axs[0].plot(nts, Ts_teo, "-", label="teoria N·Td/(1−L)")
if Ts_deg:
    axs[0].plot([a for a, _ in Ts_deg], [b for _, b in Ts_deg], "s", color="C2", label="con degrado controllato")
axs[0].set_ylim(0, 15)
axs[0].set(xlabel="numero di tracce", ylabel="rivisita della sorveglianza [s]", title="Piu' tracce = sorveglianza piu' lenta")
axs[0].grid(alpha=0.3); axs[0].legend()
_, _, _, log = simula(100, registra=True)
col = {"sorveglianza": "0.75", "traccia": "C0", "traccia in manovra": "C3"}
for t0, dt, tipo in log:
    if t0 > 0.3:
        break
    axs[1].barh(0, dt * 1e3, left=t0 * 1e3, color=col[tipo], edgecolor="none")
for tipo, c in col.items():
    axs[1].barh(0, 0, color=c, label=tipo)
axs[1].set(xlabel="tempo [ms]", yticks=[], xlim=(0, 300), title="Primi 300 ms con 100 tracce: dwell interfogliati")
axs[1].legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=3, fontsize=8)
print(" figura:", salva(fig, "p10_gestione_risorse"))

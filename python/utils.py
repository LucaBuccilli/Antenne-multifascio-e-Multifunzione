"""Funzioni condivise dalle simulazioni sugli array a fase."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIG_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(FIG_DIR, exist_ok=True)


def salva(fig, nome):
    """Salva una figura in python/figures/<nome>.png."""
    path = os.path.join(FIG_DIR, nome + ".png")
    fig.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(fig)
    return os.path.join("python", "figures", nome + ".png")


def steering(N, u, d=0.5):
    """Vettore/i di puntamento a(u) per un array lineare (d in lunghezze d'onda).
    u puo' essere uno scalare o un array: ritorna matrice len(u) x N."""
    n = np.arange(N)
    return np.exp(1j * 2 * np.pi * d * np.outer(np.atleast_1d(u), n))


def fattore_array(w, u, d=0.5):
    """Fattore di array AF(u) = sum_n w_n * exp(j 2 pi d n u)."""
    return steering(len(w), u, d) @ w


def db(x, floor=-80.0):
    """Modulo in dB normalizzato al massimo."""
    a = np.abs(x)
    a = a / a.max()
    return np.maximum(20 * np.log10(np.maximum(a, 1e-12)), floor)


def larghezza_3db(x_axis, pattern_db):
    """Larghezza del lobo principale a -3 dB (interpolata), nelle unita' di x_axis."""
    i0 = int(np.argmax(pattern_db))
    lvl = pattern_db[i0] - 3.0
    i = i0
    while i > 0 and pattern_db[i] > lvl:
        i -= 1
    xl = np.interp(lvl, [pattern_db[i], pattern_db[i + 1]], [x_axis[i], x_axis[i + 1]])
    j = i0
    while j < len(pattern_db) - 1 and pattern_db[j] > lvl:
        j += 1
    xr = np.interp(lvl, [pattern_db[j], pattern_db[j - 1]], [x_axis[j], x_axis[j - 1]])
    return xr - xl


def lobo_laterale_max(pattern_db):
    """Livello del lobo laterale piu' alto (dB rispetto al massimo)."""
    i0 = int(np.argmax(pattern_db))
    # scende dal massimo fino al primo minimo a destra e a sinistra
    l = i0
    while l > 0 and pattern_db[l - 1] < pattern_db[l]:
        l -= 1
    r = i0
    while r < len(pattern_db) - 1 and pattern_db[r + 1] < pattern_db[r]:
        r += 1
    resto = np.concatenate([pattern_db[:l], pattern_db[r + 1:]])
    return resto.max() - pattern_db[i0]


def riga(nome, misurato, teorico, unita=""):
    """Stampa una riga 'misurato vs teorico'."""
    print(f"  {nome:<42s} misurato = {misurato:>10.3f} {unita:<5s} teorico = {teorico:>10.3f} {unita}")

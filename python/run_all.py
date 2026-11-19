"""Esegue tutte le simulazioni Python in ordine e salva le figure in python/figures/."""
import glob
import os
import runpy
import time

qui = os.path.dirname(os.path.abspath(__file__))
os.chdir(qui)
t0 = time.time()
for f in sorted(glob.glob("p[0-9][0-9]_*.py")):
    print("=" * 90)
    runpy.run_path(f, run_name="__main__")
print("=" * 90)
print(f"Completato in {time.time() - t0:.1f} s")

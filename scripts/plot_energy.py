#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

prefix = 'gc_txt'
outdir = Path('plot_folder')
energy_file = outdir / f'{prefix}_energy_table.dat'
hist_file = outdir / f'{prefix}_energy_histogram.dat'

if not energy_file.exists():
    raise FileNotFoundError(energy_file)

#data = np.loadtxt(energy_file, comments='#')
data = np.loadtxt(energy_file, comments="#", usecols=(0, 1, 2))
# Format: id E_initial E_final active
E_initial = data[:, 1]
E_final = data[:, 2]
mask = (E_initial >= 0.0) & (E_final >= 0.0)
E_initial = E_initial[mask]
E_final = E_final[mask]

print(f'Particles used: {len(E_initial)}')
print(f'E_initial mean = {E_initial.mean():.8e}')
print(f'E_final   mean = {E_final.mean():.8e}')
print(f'mean dE/E      = {np.mean((E_final-E_initial)/E_initial):.8e}')

nbins = 80
emin = min(E_initial.min(), E_final.min())
emax = max(E_initial.max(), E_final.max())
if emin == emax:
    delta = 0.5*abs(emin) if emin != 0.0 else 0.5
    emin -= delta
    emax += delta
bins = np.linspace(emin, emax, nbins+1)

plt.figure(figsize=(7,5))
plt.hist(E_initial, bins=bins, histtype='step', linewidth=2, label='Inicial')
plt.hist(E_final, bins=bins, histtype='step', linewidth=2, label='Final')
plt.xlabel(r'$E_{\rm kin}$')
plt.ylabel(r'$N(E)$')
plt.legend()
plt.tight_layout()
fig = outdir / f'{prefix}_energy_distribution.png'
plt.savefig(fig, dpi=200)
print(f'Wrote {fig}')

if hist_file.exists():
    h = np.loadtxt(hist_file, comments='#')
    plt.figure(figsize=(7,5))
    plt.step(h[:,0], h[:,1], where='mid', linewidth=2, label='Inicial')
    plt.step(h[:,0], h[:,2], where='mid', linewidth=2, label='Final')
    plt.xlabel(r'$E_{\rm kin}$')
    plt.ylabel(r'$N(E)$')
    plt.legend()
    plt.tight_layout()
    fig2 = outdir / f'{prefix}_energy_histogram_from_fortran.png'
    plt.savefig(fig2, dpi=200)
    print(f'Wrote {fig2}')

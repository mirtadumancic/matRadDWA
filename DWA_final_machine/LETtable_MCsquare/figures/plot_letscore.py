import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib import cm, colors
F, OUT = sys.argv[1:3]
plt.style.use('style.mplstyle')
S = pd.read_csv(F); E = sorted(S.energy_MeV.unique())
norm = colors.Normalize(min(E), max(E)); cmap = cm.viridis
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.2))
rows = []
for e in E:
    g = S[S.energy_MeV == e]; z = g.depth_mm.values; I = g.IDD_r60.values; L = g.LETd_r60.values.copy()
    ip = I.argmax(); m = I.max(); k = ip + np.argmax(I[ip:] < 0.8 * m)
    R80 = np.interp(0.8 * m, [I[k], I[k-1]], [z[k], z[k-1]])
    keep = (I >= 0.01 * m) | (np.arange(len(z)) <= ip)          # drop the near-zero-dose tail
    c = cmap(norm(e))
    ax[0].plot(z[keep], L[keep], color=c, lw=0.9)
    ax[1].plot(z[keep] / R80, L[keep], color=c, lw=0.9)
    rows.append(dict(energy_MeV=e, R80_mm=round(R80, 1), LETd_entrance=round(L[(z > 1) & (z < 3)].mean(), 2),
                     LETd_at_R80=round(np.interp(R80, z, L), 2)))
ax[0].set_xlabel('Depth in water (mm)'); ax[0].set_xlim(0, 350)
ax[1].set_xlabel('Relative depth z / R$_{80}$'); ax[1].set_xlim(0, 1.1)
for a in ax: a.set_ylabel('LET$_d$ (keV/µm)'); a.set_ylim(0, 16)
for a, s in zip(ax, 'AB'): a.text(-0.17, 1.03, s, transform=a.transAxes, fontsize=10, fontweight='bold')
sm = cm.ScalarMappable(norm=norm, cmap=cmap)
cb = fig.colorbar(sm, ax=ax, fraction=0.03, pad=0.02); cb.set_label('Energy (MeV)')
for f in ('png', 'pdf'): fig.savefig(f'{OUT}/LETd_vsDepth_MCsquare_DWA_20energies.{f}', dpi=300, bbox_inches='tight')
pd.DataFrame(rows).to_csv(f'{OUT}/LETd_vsDepth_MCsquare_DWA_20energies_summary.csv', index=False)
print(pd.DataFrame(rows).to_string(index=False))

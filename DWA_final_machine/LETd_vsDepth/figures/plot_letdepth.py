import sys, numpy as np, pandas as pd, scipy.io as sio, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
PB, MC, MACH, OUT = sys.argv[1:5]
plt.style.use('style.mplstyle')
P = pd.read_csv(PB); M = pd.read_csv(MC)
mach = sio.loadmat(MACH, squeeze_me=True, struct_as_record=False)['machine'].data
RED, BLUE, GREY = '#C62828', '#0072B2', '#9E9E9E'
E = sorted(P.energy_MeV.unique())
fig, ax = plt.subplots(2, 3, figsize=(7.2, 4.6), sharex='col', gridspec_kw=dict(height_ratios=[2.2, 1]))
rows = []
for j, e in enumerate(E):
    p = P[np.isclose(P.energy_MeV, e)].copy(); m = M[np.isclose(M.energy_MeV, e)].copy()
    for d in (p, m):                       # no LET where there is no dose (avoid drop lines)
        for c in ('LETd_r60', 'LET_axis'): d.loc[d[c] <= 0, c] = np.nan
    k = int(np.argmin([abs(x.energy - e) for x in mach])); x = mach[k]
    z = x.depths; Z = x.Z / x.Z.max()
    ipk = np.argmax(Z); r80 = np.interp(0.8, Z[ipk:][::-1], z[ipk:][::-1])
    a = ax[0, j]; a2 = a.twinx()
    a2.fill_between(z, Z * 100, color=GREY, alpha=0.18, lw=0); a2.set_ylim(0, 110)
    a2.set_yticks([]) if j < 2 else a2.set_ylabel('IDD (% of max)', color=GREY)
    a2.spines['right'].set_visible(j == 2)
    a.set_zorder(a2.get_zorder() + 1); a.patch.set_visible(False)
    a.plot(p.depth_mm, p.LETd_r60, color=RED, lw=1.3, label='matRad PB (LET table)')
    a.plot(m.depth_mm, m.LETd_r60, color=BLUE, lw=1.3, label='MCsquare, r ≤ 60 mm')
    a.plot(m.depth_mm, m.LET_axis, color=BLUE, lw=0.9, ls='--', label='MCsquare, beam axis')
    a.set_xlim(0, r80 + 12)
    a.plot([], []); a.set_ylim(0, 16 if e < 100 else 8)
    a.set_title(f'{e:.0f} MeV (R$_{{80}}$ = {r80:.0f} mm)', loc='left', fontsize=8)
    if j == 0: a.set_ylabel('LET$_d$ (keV/µm)')
    a.axvline(r80, color=GREY, lw=0.6, ls=':')
    zz = m.depth_mm.values; ok = (zz > 2) & (zz < r80 + 3)
    rat = np.interp(zz, p.depth_mm, p.LETd_r60)
    b = ax[1, j]; b.plot(zz[ok], m.LETd_r60.values[ok] / rat[ok], color=BLUE, lw=1.2)
    b.axhline(1, color=GREY, lw=0.6); b.axvline(r80, color=GREY, lw=0.6, ls=':')
    b.set_xlabel('Depth in water (mm)'); b.set_ylim(0, 3.5)
    if j == 0: b.set_ylabel('MCsquare / PB')
    for f in [0.2, 0.5, 0.8, 0.95, 1.0]:
        zi = f * r80
        rows.append(dict(energy=e, R80=r80, z_over_R80=f, depth=zi,
            PB=np.interp(zi, p.depth_mm, p.LETd_r60), MC=np.interp(zi, m.depth_mm, m.LETd_r60),
            MC_axis=np.interp(zi, m.depth_mm, m.LET_axis)))
for a, s in zip(ax[0], 'ABC'): a.text(-0.2, 1.06, s, transform=a.transAxes, fontsize=10, fontweight='bold')
h, l = ax[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=3, fontsize=7, bbox_to_anchor=(0.5, -0.01))
fig.tight_layout(rect=(0, 0.06, 1, 1))
for f in ('png', 'pdf'): fig.savefig(f'{OUT}/LETd_vsDepth_PB_vs_MCsquare_DWA_proton_v3.{f}', dpi=300)
T = pd.DataFrame(rows); T['ratio'] = T.MC / T.PB
T.round(3).to_csv(f'{OUT}/LETd_vsDepth_points.csv', index=False); print(T.round(2).to_string(index=False))

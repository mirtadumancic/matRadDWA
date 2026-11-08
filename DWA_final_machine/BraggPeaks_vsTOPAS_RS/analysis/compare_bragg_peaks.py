"""Bragg-peak (IDD) comparison: TOPAS ground truth vs matRad vs RayStation, pencil beam vs pencil beam
and Monte Carlo vs Monte Carlo, for the DWA machine (single spot, water, r = 60 mm detector).

usage: python3 compare_bragg_peaks.py TOPAS.csv RS_PB.csv RS_MC.csv [matRad_PB.csv|-] [matRad_MC.csv|-] OUTDIR [TAG]

All curves are normalised to their own maximum (TOPAS input and RayStation exports are relative).
Metrics per energy: depth of the maximum, distal R90, R80, R50, distal fall-off R20-R80 (mm),
entrance dose at 2 mm (% of max), and the mean absolute difference from TOPAS between 2 mm and R90
(% of max). Style: scientific-visualization skill (MCP Market Hub) publication preset.
"""
import sys, os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from parse_bp import parse

TOP, RSPB, RSMC, MRPB, MRMC, OUT = sys.argv[1:7]
TAG = ('_' + sys.argv[7]) if len(sys.argv) > 7 else ''   # optional file-name suffix
os.makedirs(OUT, exist_ok=True)
plt.style.use(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'style.mplstyle'))
RED, BLK, GREY = '#C62828', 'black', '#7A7A7A'

def load_matrad(fn):
    if fn in ('-', '') or not os.path.exists(fn): return None
    T = pd.read_csv(fn); return {round(e, 2): (g.depth_mm.values, g.idd_r60.values) for e, g in T.groupby('energy_MeV')}

src = {'TOPAS': parse(TOP)[0], 'RS_PB': parse(RSPB)[0], 'RS_MC': parse(RSMC)[0],
       'mR_PB': load_matrad(MRPB), 'mR_MC': load_matrad(MRMC)}
src = {k: v for k, v in src.items() if v}
E = sorted(src['TOPAS'])
zf = np.arange(-0.5, 420.0, 0.05)

def norm_curve(z, d):
    o = np.argsort(z); z, d = np.asarray(z)[o], np.asarray(d, float)[o]
    f = np.interp(zf, z, d, left=np.nan, right=0.0)
    return 100 * f / np.nanmax(f)

def metrics(f, ref=None):
    ok = np.isfinite(f); ipk = np.nanargmax(f)
    dist = lambda lvl: zf[ipk + np.argmax(f[ipk:] < lvl)]
    r90, r80, r50, r20 = dist(90), dist(80), dist(50), dist(20)
    m = dict(peak=zf[ipk], R90=r90, R80=r80, R50=r50, falloff=r20 - r80, entrance=np.interp(2.0, zf[ok], f[ok]))
    if ref is not None:
        sel = (zf >= 2) & (zf <= metrics(ref)['R90']) & ok & np.isfinite(ref)
        m['plateauMAD'] = np.mean(np.abs(f[sel] - ref[sel]))
    return m

def match(dic, e):
    k = min(dic, key=lambda x: abs(x - e)); return dic[k] if abs(k - e) < 0.06 else None

curves, rows = {}, []
for e in E:
    ref = norm_curve(*src['TOPAS'][e]); curves[('TOPAS', e)] = ref
    rows.append(dict(source='TOPAS', energy=e, **metrics(ref)))
    for s in [k for k in src if k != 'TOPAS']:
        c = match(src[s], e)
        if c is None: continue
        f = norm_curve(*c); curves[(s, e)] = f
        rows.append(dict(source=s, energy=e, **metrics(f, ref)))
M = pd.DataFrame(rows)
T0 = M[M.source == 'TOPAS'].set_index('energy')
for c in ['peak', 'R90', 'R80', 'R50', 'falloff', 'entrance']:
    M['d_' + c] = M[c].values - T0.loc[M.energy, c].values
M.to_csv(f'{OUT}/BraggPeak_metrics_all_energies{TAG}.csv', index=False)

LAB = {'mR_PB': 'matRad pencil beam', 'RS_PB': 'RayStation pencil beam',
       'mR_MC': 'matRad MC (MCsquare)', 'RS_MC': 'RayStation MC'}
if TAG: LAB['mR_PB'] += f' ({TAG[1:]})'
STY = {'mR_PB': dict(color=RED, ls='-'), 'RS_PB': dict(color=BLK, ls='--'),
       'mR_MC': dict(color=RED, ls='-'), 'RS_MC': dict(color=BLK, ls='--')}

# ---------- figure 1: IDD overlays, PB row and MC row, difference panels below ----------
show = [e for e in [50.03, 80.05, 110.11, 149.16, 190.30, 228.84] if e in E]
fig, axs = plt.subplots(4, 1, figsize=(7.0, 7.4), sharex=True, gridspec_kw=dict(height_ratios=[3, 1.3, 3, 1.3]))
for r, (kind, pair) in enumerate([('Pencil beam', ['mR_PB', 'RS_PB']), ('Monte Carlo', ['mR_MC', 'RS_MC'])]):
    a, ad = axs[2 * r], axs[2 * r + 1]
    for i, e in enumerate(show):
        ref = curves[('TOPAS', e)]
        a.plot(zf, ref, color=GREY, lw=2.6, alpha=0.55, label='TOPAS (input)' if i == 0 else None)
        for s in pair:
            if (s, e) in curves:
                a.plot(zf, curves[(s, e)], lw=1.1, label=LAB[s] if i == 0 else None, **STY[s])
                ad.plot(zf, curves[(s, e)] - ref, lw=0.9, **STY[s])
        pk = zf[np.nanargmax(ref)]
        a.text(pk, 104, f'{e:.0f}', ha='center', fontsize=6.5, color=GREY)
    a.set_ylim(0, 112); a.set_ylabel('IDD (% of max)')
    a.set_title(kind, loc='left')
    a.legend(loc='lower right', bbox_to_anchor=(1.0, 1.0), ncol=3, fontsize=6.5)
    ad.axhline(0, color=GREY, lw=0.6); ad.set_ylim(-6, 6); ad.set_ylabel('Δ vs TOPAS\n(% of max)')
axs[-1].set_xlabel('Depth in water (mm)'); axs[-1].set_xlim(0, 340)
for a, s in zip(axs[[0, 2]], 'AB'):
    a.text(-0.09, 1.04, s, transform=a.transAxes, fontsize=10, fontweight='bold', va='bottom')
for x in ('pdf', 'png'): fig.savefig(f'{OUT}/BraggPeaks_IDD_vs_TOPAS{TAG}.{x}', dpi=300)

# ---------- figure 2: metric differences vs energy ----------
MK = {'mR_PB': dict(color=RED, marker='o', mfc=RED, ls='-'), 'RS_PB': dict(color=BLK, marker='s', mfc=BLK, ls='--'),
      'mR_MC': dict(color=RED, marker='o', mfc='white', ls='-'), 'RS_MC': dict(color=BLK, marker='s', mfc='white', ls='--')}
panels = [('d_R80', 'ΔR80 vs TOPAS (mm)', None), ('d_falloff', 'Δ fall-off R20−R80\nvs TOPAS (mm)', None),
          ('d_entrance', 'Δ entrance dose at 2 mm\nvs TOPAS (% of max)', None), ('plateauMAD', 'Mean |Δ| vs TOPAS, 2 mm to R90\n(% of max)', None)]
fig2, ax2 = plt.subplots(2, 2, figsize=(7.0, 5.0), sharex=True)
for a, (c, lab, tol) in zip(ax2.flat, panels):
    for s in [k for k in ['mR_PB', 'RS_PB', 'mR_MC', 'RS_MC'] if k in src]:
        d = M[M.source == s]
        a.plot(d.energy, d[c], ms=3.2, lw=0.8, label=LAB[s], **MK[s])
    if tol: a.axhspan(-tol, tol, color=GREY, alpha=0.12, lw=0)
    if c != 'plateauMAD': a.axhline(0, color=GREY, lw=0.6)
    a.set_ylabel(lab)
for a in ax2[1]: a.set_xlabel('Energy (MeV)')
ax2[0, 0].legend(loc='best', fontsize=6.5)
for a, s in zip(ax2.flat, 'ABCD'):
    a.text(-0.2, 1.04, s, transform=a.transAxes, fontsize=10, fontweight='bold', va='bottom')
for x in ('pdf', 'png'): fig2.savefig(f'{OUT}/BraggPeaks_metrics_vs_energy{TAG}.{x}', dpi=300)

# ---------- summary table ----------
S = (M[M.source != 'TOPAS'].groupby('source')
     .agg(n=('energy', 'size'), dR80_mean=('d_R80', 'mean'), dR80_maxabs=('d_R80', lambda x: np.max(np.abs(x))),
          dPeak_mean=('d_peak', 'mean'), dFalloff_mean=('d_falloff', 'mean'), dEntr_mean=('d_entrance', 'mean'),
          plateauMAD_mean=('plateauMAD', 'mean'), plateauMAD_max=('plateauMAD', 'max')).round(2))
S.to_csv(f'{OUT}/BraggPeak_summary{TAG}.csv'); print(S.to_string())

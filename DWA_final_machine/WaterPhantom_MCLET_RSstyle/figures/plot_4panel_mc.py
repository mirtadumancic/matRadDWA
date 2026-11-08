import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
DVH, LVH, RSD, RSL, OUT = sys.argv[1:6]
plt.style.use('/home/claude/bp/style.mplstyle')
BLUE, BLK = '#0072B2', 'black'
def parse(f):
    R, cur = {}, None
    for l in open(f, encoding='utf-8', errors='ignore'):
        l = l.strip()
        if l.startswith('#RoiName:'): cur = l.split(':', 1)[1]; R[cur] = []
        elif l and not l.startswith('#') and cur: R[cur].append([float(x) for x in l.split()])
    return {k: np.array(v) for k, v in R.items() if len(v)}
rsD, rsL = parse(RSD), parse(RSL)
rs = {('CTV', 'D'): rsD['CTV_Centre_5cm'] * [0.01, 1], ('Ext', 'D'): rsD['Ext-CTVcen5'] * [0.01, 1],
      ('CTV', 'L'): rsL['CTV_Centre_5cm'], ('Ext', 'L'): rsL['Ext-CTVcen5']}
Dv = pd.read_csv(DVH); Lv = pd.read_csv(LVH)
R = {'CTV': 'CTV', 'Ext': 'PhantomMinusCTV'}
d95 = lambda x, v: np.interp(-95, -np.asarray(v), np.asarray(x))
c = Dv[Dv.region == 'CTV']; d95_rs = d95(rs[('CTV', 'D')][:, 0], rs[('CTV', 'D')][:, 1])
sc = d95_rs / d95(c.dose_Gy.values, c.volume_pct.values); print(f'RS D95 {d95_rs:.2f}; MC scale {sc:.4f}')
mc = {}
for r in R:
    t = Dv[Dv.region == R[r]]; mc[(r, 'D')] = (t.dose_Gy.values * sc, t.volume_pct.values)
    t = Lv[(Lv.source == 'MC_StopPow') & (Lv.region == R[r]) & (Lv.doseCut_fracRx == 0)]
    mc[(r, 'L')] = (t.LETd_keVum.values, t.volume_pctROI.values)
fig, ax = plt.subplots(2, 2, figsize=(7.2, 5.6))
panels = [(0, 0, 'CTV', 'D', 'CTV (5 cm, centre): DVH'), (0, 1, 'Ext', 'D', 'Phantom − CTV: DVH'),
          (1, 0, 'CTV', 'L', 'CTV (5 cm, centre): LVH'), (1, 1, 'Ext', 'L', 'Phantom − CTV: LVH')]
for i, j, r, q, t in panels:
    a = ax[i, j]; x, v = mc[(r, q)]
    a.plot(x, v, color=BLUE, lw=1.5, label='DWA 1 mm matRad MCsquare')
    a.plot(rs[(r, q)][:, 0], rs[(r, q)][:, 1], color=BLK, lw=1.4, ls='--', label='DWA 1 mm RayStation MC')
    a.set_title(t, loc='left'); a.set_ylabel('Volume (%)')
    a.set_xlabel('Dose (Gy(RBE))' if q == 'D' else 'LET$_d$ (keV/µm)')
    if r == 'Ext': a.set_yscale('log'); a.set_ylim(1e-2, 110); a.set_xlim(0, 60 if q == 'D' else 25)
    else: a.set_ylim(0, 105); a.set_xlim((45, 60) if q == 'D' else (1.5, 4.5))
for a, s in zip(ax.flat, 'ABCD'): a.text(-0.16, 1.06, s, transform=a.transAxes, fontsize=10, fontweight='bold')
ax[0, 0].legend(loc='lower left', fontsize=7)
fig.tight_layout()
for f in ('png', 'pdf'): fig.savefig(f'{OUT}/DWA_matRadMCsquare_vs_RayStationMC_DVH_LVH.{f}', dpi=300)
# metrics
def lat(x, v, q): return np.interp(-q, -np.asarray(v), np.asarray(x))
def mean(x, v): o = np.argsort(x); return np.trapezoid(np.asarray(v)[o] / 100, np.asarray(x)[o])
rows = []
for name, d in [('matRad MCsquare', mc), ('RayStation MC', {k: (rs[k][:, 0], rs[k][:, 1]) for k in rs})]:
    xd, vd = d[('CTV', 'D')]; xl, vl = d[('CTV', 'L')]
    rows.append(dict(source=name, D98=lat(xd, vd, 98), D95=lat(xd, vd, 95), D50=lat(xd, vd, 50), D2=lat(xd, vd, 2),
                     LETd_mean=mean(xl, vl), L98=lat(xl, vl, 98), L50=lat(xl, vl, 50), L2=lat(xl, vl, 2)))
T = pd.DataFrame(rows).set_index('source').round(3); T.to_csv(f'{OUT}/DWA_matRadMCsquare_vs_RayStationMC_CTV_metrics.csv'); print(T.to_string())

import sys, os, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from load import load
OLD, NEWCSV, RSD, RSL, OUT = sys.argv[1:6]
plt.style.use('/home/claude/bp/style.mplstyle')
RED, BLUE, BLK = '#C62828', '#0072B2', 'black'
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
D, L, m, _ = load(OLD)
reg = {'CTV': m['tMask'], 'Ext': m['outsideEll']}
def cum(v, g): v = np.sort(v); return g, 100 * (1 - np.searchsorted(v, g, side='left') / v.size)
N = pd.read_csv(NEWCSV)
new = lambda r, q: N[(N.region == {'CTV': 'CTV', 'Ext': 'PhantomMinusCTV'}[r]) & (N.quantity == q)]
# scale and match: both matRad doses normalised to the RayStation CTV D95 (LET unaffected)
d95 = lambda x, v: np.interp(-95, -np.asarray(v), np.asarray(x))
rsC = rs[('CTV', 'D')]; d95_rs = d95(rsC[:, 0], rsC[:, 1])
o = np.sort(D[reg['CTV']]); scOld = d95_rs / o[int(np.ceil(0.05 * o.size)) - 1]
nC = N[(N.region == 'CTV') & (N.quantity == 'dose_Gy')]; scNew = d95_rs / d95(nC.x.values, nC.volume_pct.values)
D = D * scOld; N.loc[N.quantity == 'dose_Gy', 'x'] *= scNew
print(f'RS CTV D95 {d95_rs:.2f} Gy; scale old {scOld:.4f}, new {scNew:.4f}')
fig, ax = plt.subplots(2, 2, figsize=(7.2, 5.6))
panels = [(0, 0, 'CTV', 'D', 'CTV (5 cm, centre): DVH'), (0, 1, 'Ext', 'D', 'Phantom − CTV: DVH'),
          (1, 0, 'CTV', 'L', 'CTV (5 cm, centre): LVH'), (1, 1, 'Ext', 'L', 'Phantom − CTV: LVH')]
for i, j, r, q, t in panels:
    a = ax[i, j]
    g = np.linspace(0, 65, 3000) if q == 'D' else np.linspace(0, 35, 3000)
    x, v = cum((D if q == 'D' else L)[reg[r]], g)
    a.plot(x, v, color=RED, lw=1.5, label='DWA 1 mm matRad')
    n = new(r, 'dose_Gy' if q == 'D' else 'LETd_keVum')
    a.plot(n.x, n.volume_pct, color=BLUE, lw=1.5, label='DWA 1 mm matRad NEW')
    a.plot(rs[(r, q)][:, 0], rs[(r, q)][:, 1], color=BLK, lw=1.4, ls='--', label='DWA 1 mm RayStation')
    a.set_title(t, loc='left'); a.set_ylabel('Volume (%)')
    a.set_xlabel('Dose (Gy(RBE))' if q == 'D' else 'LET$_d$ (keV/µm)')
    if r == 'Ext': a.set_yscale('log'); a.set_ylim(1e-2, 110); a.set_xlim(0, 60 if q == 'D' else 25)
    else:
        a.set_ylim(0, 105); a.set_xlim((45, 60) if q == 'D' else (1.5, 4.5))
for a, s in zip(ax.flat, 'ABCD'): a.text(-0.16, 1.06, s, transform=a.transAxes, fontsize=10, fontweight='bold')
ax[0, 0].legend(loc='lower left', fontsize=7)
fig.tight_layout()
for f in ('png', 'pdf'): fig.savefig(f'{OUT}/DWA_matRad_vs_matRadNEW_vs_RayStation_DVH_LVH.{f}', dpi=300)

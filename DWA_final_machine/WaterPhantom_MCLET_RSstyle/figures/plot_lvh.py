import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
LVH, RSL, OUT = sys.argv[1:4]
plt.style.use('style.mplstyle')
H = pd.read_csv(LVH)
def parse(f):
    R, cur = {}, None
    for l in open(f, encoding='utf-8', errors='ignore'):
        l = l.strip()
        if l.startswith('#RoiName:'): cur = l.split(':', 1)[1]; R[cur] = []
        elif l and not l.startswith('#') and cur: R[cur].append([float(x) for x in l.split()])
    return {k: np.array(v) for k, v in R.items() if len(v)}
rs = parse(RSL)
rsC, rsE = rs['CTV_Centre_5cm'], rs['Ext-CTVcen5']
def Lat(x, v, q):  # LET at volume q % (cumulative)
    o = np.argsort(x); x, v = x[o], v[o]
    return np.interp(-q, -v, x)
def meanL(x, v):
    o = np.argsort(x); x, v = x[o], v[o] / 100
    return np.trapezoid(v, x)
S = {'PB_table_cut1pct': ('matRad PB, table, 1 % cut', dict(color='#C62828', ls='-')),
     'PB_table_noCut':   ('matRad PB, table, no cut', dict(color='#C62828', ls=':')),
     'MC_StopPow':       ('matRad MCsquare, StopPow, no cut', dict(color='#0072B2', ls='-')),
     'MC_DepositedEnergy': ('matRad MCsquare, DepositedEnergy, no cut', dict(color='#0072B2', ls='-.'))}
RSsty = dict(color='black', ls='--')
rows = []
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.1))
for k, (lab, st) in S.items():
    c = H[(H.source == k) & (H.region == 'CTV') & (H.doseCut_fracRx == 0)]
    ax[0].plot(c.LETd_keVum, c.volume_pctROI, lw=1.3, label=lab, **st)
    x, v = c.LETd_keVum.values, c.volume_pctROI.values
    rows.append(dict(source=lab, mean=meanL(x, v), L98=Lat(x, v, 98), L50=Lat(x, v, 50), L2=Lat(x, v, 2)))
    e = H[(H.source == k) & (H.region == 'PhantomMinusCTV') & (H.doseCut_fracRx == 0)]
    ax[1].plot(e.LETd_keVum, e.volume_pctROI, lw=1.3, **st)
ax[0].plot(rsC[:, 0], rsC[:, 1], lw=1.3, label='RayStation MC', **RSsty)
ax[1].plot(rsE[:, 0], rsE[:, 1], lw=1.3, **RSsty)
x, v = rsC[:, 0], rsC[:, 1]
rows.append(dict(source='RayStation MC', mean=meanL(x, v), L98=Lat(x, v, 98), L50=Lat(x, v, 50), L2=Lat(x, v, 2)))
ax[0].set_xlim(1.5, 4.5); ax[0].set_ylim(0, 105); ax[0].set_title('CTV (5 cm, centre)', loc='left')
ax[1].set_yscale('log'); ax[1].set_ylim(1e-3, 110); ax[1].set_xlim(0, 30); ax[1].set_title('Phantom − CTV (no dose cut-off)', loc='left')
for a in ax: a.set_xlabel('LET$_d$ (keV/µm)'); a.set_ylabel('Volume (%)')
for a, s in zip(ax, 'AB'): a.text(-0.16, 1.04, s, transform=a.transAxes, fontsize=10, fontweight='bold', va='bottom')
h, l = ax[0].get_legend_handles_labels()
fig.legend(h, l, loc='lower center', ncol=3, fontsize=6.5, bbox_to_anchor=(0.5, -0.02))
fig.tight_layout(rect=(0, 0.12, 1, 1))
for x in ('png', 'pdf'): fig.savefig(f'{OUT}/LVH_CTV_matRadMC_vs_RayStation.{x}', dpi=300)
T = pd.DataFrame(rows).set_index('source').round(3)
rsrow = T.loc['RayStation MC']
for c in ['mean', 'L98', 'L50', 'L2']: T[c + '_vsRS_pct'] = (100 * (T[c] / rsrow[c] - 1)).round(1)
T.to_csv(f'{OUT}/LVH_CTV_metrics.csv'); print(T.to_string())
